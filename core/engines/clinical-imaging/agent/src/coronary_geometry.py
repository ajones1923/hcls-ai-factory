"""Geometric stenosis analysis computed from a segmented coronary surface.

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
`ct_coronary_angiography.infer()` previously returned fixed numbers — CAD-RADS 4A, Agatston 385,
72 % stenosis — regardless of input. This module replaces the stenosis half of that with an actual
measurement.

It is NOT coronary segmentation from a CT volume. That is not possible with what the repo holds:
there is no cardiac CTA volume (the only CT series are an abdomen study and a 64-slice chest
series), and no coronary model weights — the SegResNet the workflow docstring names is absent, and
VISTA-3D is a general organ segmenter, not a coronary tool.

What IS available is `data/cardiac_ct/CoronariesNC6/`: real coronary artery surfaces with two
independent manual rater ground truths per case. Those are segmentation OUTPUT. So the honest
computation is the step that follows segmentation in a real pipeline — quantitative stenosis
analysis of the vessel lumen:

    surface mesh -> per-point lumen radius -> radius profile -> stenosis -> CAD-RADS

Every number this returns is derived from the geometry. Nothing is hardcoded.

METHOD
------
Local lumen radius comes from the **shrinking-ball medial axis** (Ma, Bae & Choi 2012): for a
surface point p with inward normal -n, the maximal inscribed ball tangent at p has its centre on
the medial axis, and its radius is the local vessel radius. Iterating

    c = p - r*n ;  q = nearest other surface point to c ;  r = |p-q|^2 / (2 * (p-q).n)

converges in a handful of steps and needs only a KD-tree.

Diameter stenosis is then measured against a *local proximal reference*, the way a reader grades a
CTA: compare the narrowest point to healthy calibre nearby, not to the whole-vessel mean, so a
tapering distal vessel is not scored as a lesion.

CAD-RADS follows deterministically from maximum diameter stenosis per the published categories
(Cury et al., JACC Cardiovasc Imaging 2016; 2022 update), so that grade is a lookup on a measured
value rather than an assertion.

LIMITS — stated, not buried
---------------------------
* Calcium scoring (Agatston) is NOT computable here. Agatston needs CT Hounsfield units and slice
  geometry; a surface mesh carries neither. Any calcium figure remains representative and is
  labelled as such by the caller.
* Branch identity (LAD / LCx / RCA) is NOT resolved. The meshes carry no anatomical labels, so
  lesions are reported by segment index, not by vessel name.
* Rater variability is real: the two manual ground truths per case differ, and comparing them is a
  legitimate uncertainty estimate rather than noise to hide.
"""
from __future__ import annotations

import math
import pathlib
from dataclasses import dataclass, field, asdict
from typing import Optional

import numpy as np

try:
    from scipy.spatial import cKDTree
    _HAVE_SCIPY = True
except Exception:                                            # pragma: no cover
    _HAVE_SCIPY = False


# CAD-RADS categories by maximum diameter stenosis (Cury et al.). 4B and 5 additionally require
# vessel-level or occlusion context we cannot resolve from an unlabelled mesh, so the ceiling
# reachable here is 4A; the caller is told which categories are unreachable.
CAD_RADS_BANDS = [
    (0.0,   0.0,  "0",  "No plaque"),
    (0.01, 24.0,  "1",  "Minimal non-obstructive"),
    (25.0, 49.0,  "2",  "Mild non-obstructive"),
    (50.0, 69.0,  "3",  "Moderate stenosis"),
    (70.0, 99.0,  "4A", "Severe stenosis"),
    (100.0, 100.0, "5", "Total occlusion"),
]
UNREACHABLE_CATEGORIES = ["4B (left main >50% or 3-vessel obstructive) — needs vessel labels"]


@dataclass
class Lesion:
    segment_index: int
    min_radius_mm: float
    reference_radius_mm: float
    diameter_stenosis_pct: float
    area_stenosis_pct: float
    position_mm: list = field(default_factory=list)


@dataclass
class CoronaryGeometryResult:
    source: str
    n_surface_points: int
    n_points_sampled: int
    total_vessel_length_mm: float
    mean_radius_mm: float
    min_radius_mm: float
    max_radius_mm: float
    max_diameter_stenosis_pct: float
    cad_rads: str
    cad_rads_label: str
    lesions: list
    method: str = "shrinking-ball medial axis (Ma et al. 2012); CAD-RADS per Cury et al."
    computed: bool = True
    unreachable_categories: list = field(default_factory=lambda: list(UNREACHABLE_CATEGORIES))
    caveats: list = field(default_factory=lambda: [
        "Geometric analysis of a segmented surface — not segmentation of a patient CT volume.",
        "Agatston calcium score is not computable from a surface mesh (needs CT Hounsfield units).",
        "Branch identity (LAD/LCx/RCA) is not resolved; lesions are reported by segment index.",
        "Decision support only — not a diagnosis.",
        "Lesions must be interior to the vessel (healthy calibre both sides); tips and end caps "
        "are excluded because the medial-axis radius degenerates there.",
    ])

    def to_dict(self):
        return asdict(self)


def read_stl_vertices(path: pathlib.Path):
    """Unique surface points and per-point normals from a binary STL."""
    b = pathlib.Path(path).read_bytes()
    if b[:5].lower() == b"solid" and b"facet" in b[:512]:
        raise ValueError(f"{path}: ASCII STL not supported")
    n = int.from_bytes(b[80:84], "little")
    rec = np.frombuffer(b[84:84 + n * 50],
                        dtype=np.dtype([("n", "<3f4"), ("v", "<3,3f4"), ("attr", "<u2")]),
                        count=n)
    tri = rec["v"].astype(np.float64)
    fn = rec["n"].astype(np.float64)
    pts = tri.reshape(-1, 3)
    nrm = np.repeat(fn, 3, axis=0)
    # collapse duplicate vertices, averaging their face normals into a point normal
    key = np.round(pts, 4)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    P = pts[idx]
    N = np.zeros_like(P)
    np.add.at(N, inv, nrm)
    ln = np.linalg.norm(N, axis=1, keepdims=True)
    N = np.divide(N, np.where(ln == 0, 1, ln))
    return P, N


def local_radii(P, N, sample=20000, iters=24, seed=0):
    """Per-point lumen radius via the shrinking-ball medial axis. Returns (points, radii)."""
    if not _HAVE_SCIPY:
        raise RuntimeError("scipy is required for medial-axis radius estimation")
    rng = np.random.default_rng(seed)
    if len(P) > sample:
        sel = rng.choice(len(P), sample, replace=False)
    else:
        sel = np.arange(len(P))
    Q, NQ = P[sel], N[sel]
    tree = cKDTree(P)

    span = float(np.linalg.norm(P.max(0) - P.min(0)))
    r = np.full(len(Q), span / 8.0)
    for _ in range(iters):
        c = Q - NQ * r[:, None]
        d, j = tree.query(c, k=2)
        # ignore the tangent point itself
        qd = np.where(np.isclose(d[:, 0], r, atol=1e-6), d[:, 1], d[:, 0])
        qi = np.where(np.isclose(d[:, 0], r, atol=1e-6), j[:, 1], j[:, 0])
        nearest = P[qi]
        v = Q - nearest
        denom = 2.0 * np.einsum("ij,ij->i", v, NQ)
        with np.errstate(divide="ignore", invalid="ignore"):
            r_new = np.einsum("ij,ij->i", v, v) / denom
        ok = np.isfinite(r_new) & (r_new > 1e-6) & (r_new < r)
        if not ok.any():
            break
        r = np.where(ok, r_new, r)
    good = np.isfinite(r) & (r > 1e-3) & (r < span)
    return Q[good], r[good]


def analyse(stl_path, sample=20000, ref_neighbourhood_mm=12.0, seed=0) -> CoronaryGeometryResult:
    """Measure lumen radii and grade the worst stenosis."""
    stl_path = pathlib.Path(stl_path)
    P, N = read_stl_vertices(stl_path)
    C, R = local_radii(P, N, sample=sample, seed=seed)
    if len(R) < 100:
        raise RuntimeError(f"{stl_path.name}: too few medial samples ({len(R)})")

    tree = cKDTree(C)

    # Measure over a SEGMENT, not a point.
    #
    # The raw medial radii are anatomically right in bulk (median 0.97 mm, p99 1.98 mm for a
    # coronary tree) but carry a ~3 % tail below 0.3 mm from mesh noise and branch-junction
    # degeneracy. A point minimum therefore finds noise every time -- the first version reported
    # four "lesions" all pinned to the plausibility floor. Quantitative coronary angiography grades
    # a lesion over a length of vessel, so smooth the radius field over ~2 mm first; a true
    # narrowing survives that, a one-point artifact does not.
    smooth_ids = tree.query_ball_point(C, r=2.0)
    R = np.array([np.median(R[ids]) if len(ids) >= 5 else R[i]
                  for i, ids in enumerate(smooth_ids)])

    # A local reference calibre: the 75th-percentile radius among nearby medial samples. Using a
    # local window rather than a global mean is what stops normal distal tapering being graded as
    # disease -- the same reason a reader compares a lesion to the adjacent healthy segment.
    ref = np.empty(len(R))
    nbrs = tree.query_ball_point(C, r=ref_neighbourhood_mm)
    for i, ids in enumerate(nbrs):
        ref[i] = np.percentile(R[ids], 75) if len(ids) >= 8 else R[i]

    with np.errstate(divide="ignore", invalid="ignore"):
        dia_sten = np.clip((1.0 - R / ref) * 100.0, 0.0, 100.0)
    dia_sten = np.nan_to_num(dia_sten)

    # A narrowing is only a LESION if the vessel continues past it on BOTH sides.
    #
    # Without this, the top "stenoses" are all artifacts: the shrinking ball degenerates at vessel
    # tips and end caps, where the maximal inscribed sphere collapses toward zero. The first run of
    # this module duly reported 98.5 % stenosis at a 0.026 mm radius -- a mesh endpoint, not
    # disease. A true stenosis has healthy calibre proximal AND distal to it, which is also how a
    # reader grades one, so require the local medial samples to extend both ways along the vessel
    # axis and reject radii below any plausible lumen.
    MIN_PLAUSIBLE_RADIUS_MM = 0.20

    def is_interior(i):
        ids = nbrs[i]
        if len(ids) < 12:
            return False
        local = C[ids] - C[i]
        axis_l = np.linalg.svd(local - local.mean(0), full_matrices=False)[2][0]
        t = local @ axis_l
        reach = ref_neighbourhood_mm / 3.0
        return t.min() < -reach and t.max() > reach

    # Report the worst lesion plus any other well-separated significant narrowings.
    order = np.argsort(-dia_sten)
    lesions, taken = [], []
    for i in order:
        if dia_sten[i] < 25.0 or len(lesions) >= 4:
            break
        if R[i] < MIN_PLAUSIBLE_RADIUS_MM or not is_interior(i):
            continue
        if any(np.linalg.norm(C[i] - C[j]) < ref_neighbourhood_mm for j in taken):
            continue
        taken.append(i)
        lesions.append(Lesion(
            segment_index=len(lesions) + 1,
            min_radius_mm=round(float(R[i]), 3),
            reference_radius_mm=round(float(ref[i]), 3),
            diameter_stenosis_pct=round(float(dia_sten[i]), 1),
            area_stenosis_pct=round(float((1.0 - (R[i] / ref[i]) ** 2) * 100.0), 1),
            position_mm=[round(float(x), 2) for x in C[i]],
        ))

    max_sten = float(lesions[0].diameter_stenosis_pct) if lesions else 0.0
    grade, label = "0", "No plaque"
    for lo, hi, g, lab in CAD_RADS_BANDS:
        if lo <= max_sten <= hi:
            grade, label = g, lab
    # crude vessel length: extent of the medial samples along their principal axis
    ctr = C - C.mean(0)
    axis = np.linalg.svd(ctr, full_matrices=False)[2][0]
    proj = ctr @ axis
    length = float(proj.max() - proj.min())

    return CoronaryGeometryResult(
        source=str(stl_path.relative_to(stl_path.parents[3]) if len(stl_path.parents) > 3 else stl_path),
        n_surface_points=int(len(P)),
        n_points_sampled=int(len(R)),
        total_vessel_length_mm=round(length, 1),
        mean_radius_mm=round(float(np.mean(R)), 3),
        min_radius_mm=round(float(np.min(R)), 3),
        max_radius_mm=round(float(np.max(R)), 3),
        max_diameter_stenosis_pct=round(max_sten, 1),
        cad_rads=grade,
        cad_rads_label=label,
        lesions=[asdict(x) for x in lesions],
    )


def rater_agreement(stl_a, stl_b, **kw):
    """Compare two raters' surfaces — a real uncertainty estimate, not a hidden one."""
    a, b = analyse(stl_a, **kw), analyse(stl_b, **kw)
    return {
        "rater_a_max_stenosis_pct": a.max_diameter_stenosis_pct,
        "rater_b_max_stenosis_pct": b.max_diameter_stenosis_pct,
        "abs_difference_pct": round(abs(a.max_diameter_stenosis_pct
                                        - b.max_diameter_stenosis_pct), 1),
        "cad_rads_agree": a.cad_rads == b.cad_rads,
    }


if __name__ == "__main__":                                    # pragma: no cover
    import argparse, json
    ap = argparse.ArgumentParser()
    ap.add_argument("stl")
    ap.add_argument("--sample", type=int, default=20000)
    a = ap.parse_args()
    print(json.dumps(analyse(a.stl, sample=a.sample).to_dict(), indent=2))
