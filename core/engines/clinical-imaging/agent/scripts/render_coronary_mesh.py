#!/usr/bin/env python3
"""Render the CoronariesNC6 coronary-artery meshes to demo imagery.

WHY THIS EXISTS
---------------
The cardiac demo previously illustrated "CTA Coronary — AI Stenosis Analysis" with an axial CT
slice of a VERTEBRA, and its "V2 — Real Clinical CT" fallback with an upper-abdominal slice; a
third asset titled "Real Clinical CT — Heart Level" showed the kidneys. All three were drawn from
`data/cardiac_ct/abdomen_ct/`, which is honestly named. A clinician spots that in about a second,
and it costs exactly the credibility the demo is meant to build.

`data/cardiac_ct/CoronariesNC6/` holds what the demo actually needed the whole time: real coronary
artery surface meshes, six cases, each with a machine segmentation (`ML.stl`) and two independent
manual rater ground truths. This renders them, so the panel shows genuine coronary anatomy that can
be labelled truthfully — a 3-D vessel model, not a CT slice it is not.

No mesh library is required (no trimesh/pyvista/vtk in this environment): binary STL is parsed with
numpy and rasterised with a painter's-algorithm fill, which also keeps full control of the palette.

Usage:
    python3 scripts/render_coronary_mesh.py                 # default case, all outputs
    python3 scripts/render_coronary_mesh.py --case V1P2 --spin
"""
from __future__ import annotations
import argparse, math, pathlib, sys
import numpy as np
from PIL import Image, ImageDraw, ImageFont

BASE = pathlib.Path(__file__).resolve().parents[1]
MESH_DIR = BASE / "data" / "cardiac_ct" / "CoronariesNC6"
OUT_DIR = BASE / "data" / "demo" / "coronary"

# Portal palette (the imaging UI is dark; these match its accents)
BG = (10, 12, 16)
VESSEL = (214, 226, 236)
VESSEL_DEEP = (120, 148, 172)
NVIDIA_GREEN = (118, 185, 0)
AMBER = (255, 176, 32)
RED = (255, 72, 72)
MUTED = (150, 162, 175)

SS = 2                      # supersample factor for antialiasing
# One title size for every panel. A3 and B1 must MATCH, so the size lives in one place; both
# panels pass it through _fit(), which can only shrink it if a string somehow does not fit.
# 72 in 1024-space (the file's `n * SS` idiom), i.e. ~53 px on a 760 px panel -- about 1.8x the
# 30 px the one-line title was width-limited to. 80 * SS is the ceiling before the longer line
# stops fitting, so this leaves headroom for _fit() rather than relying on it.
TITLE_PX = 72 * SS


def read_stl(path: pathlib.Path):
    """Vertices (T,3,3) and face normals (T,3) from a binary STL."""
    b = path.read_bytes()
    if b[:5].lower() == b"solid" and b"facet" in b[:512]:
        raise ValueError(f"{path.name}: ASCII STL not supported")
    n = int.from_bytes(b[80:84], "little")
    rec = np.frombuffer(b[84:84 + n * 50],
                        dtype=np.dtype([("n", "<3f4"), ("v", "<3,3f4"), ("attr", "<u2")]),
                        count=n)
    return rec["v"].astype(np.float64), rec["n"].astype(np.float64)


def rotation(rx, ry, rz):
    cx, sx, cy, sy, cz, sz = (math.cos(rx), math.sin(rx), math.cos(ry),
                              math.sin(ry), math.cos(rz), math.sin(rz))
    Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    return Rz @ Ry @ Rx


def render(tris, normals, size=1024, rx=-1.15, ry=0.55, rz=0.0, margin=0.10,
           focus_world=None, span_mm=None, highlight_world=None, highlight_mm=3.0,
           shell=None, shell_alpha=132):
    """Painter's-algorithm render with Lambertian shading. Returns (image, projector).

    focus_world / span_mm let the camera frame a specific point at a specific field of view --
    used for the stenosis close-up, so the detail view is generated from the SAME geometry and the
    SAME measurement as the overview rather than being drawn by hand.
    """
    R = rotation(rx, ry, rz)
    v = tris.reshape(-1, 3) @ R.T
    nrm = normals @ R.T

    lo, hi = v.min(0), v.max(0)
    if shell is not None and shell[0] is not None:
        sv = shell[0].reshape(-1, 3) @ R.T          # include the heart in the framing
        lo = np.minimum(lo, sv.min(0)); hi = np.maximum(hi, sv.max(0))
    ctr = (lo + hi) / 2.0
    span = float(max(hi[0] - lo[0], hi[1] - lo[1])) or 1.0
    if focus_world is not None:
        ctr = np.asarray(focus_world, float) @ R.T
    if span_mm:
        span = float(span_mm)
    W = size * SS
    scale = (W * (1 - 2 * margin)) / span

    def project(p3):
        q = (np.atleast_2d(p3) - ctr) * scale
        x = q[:, 0] + W / 2.0
        y = W / 2.0 - q[:, 1]                       # screen y grows downward
        return np.stack([x, y], 1)

    p = project(v).reshape(-1, 3, 2)
    depth = v.reshape(-1, 3, 3)[:, :, 2].mean(1)

    # Lambertian from a front-left-ish key light, plus a little ambient
    L = np.array([-0.45, 0.55, 0.75]); L /= np.linalg.norm(L)
    ndl = np.clip(nrm @ L, 0, 1)
    shade = 0.24 + 0.76 * ndl

    # Tint the geometry AT the lesion rather than drawing a marker over it. Because the colour
    # belongs to the triangles, it rotates with the vessel and is correctly hidden when that part
    # of the artery turns away -- an overlay ring would float in front of the mesh and read as a
    # sticker. Falls off smoothly so the patch has no hard edge.
    tint = None
    if highlight_world is not None:
        centres = tris.reshape(-1, 3, 3).mean(1) @ R.T
        dist = np.linalg.norm(centres - (np.asarray(highlight_world, float) @ R.T), axis=1)
        # gamma < 1 widens the visible core instead of feathering it away
        tint = np.clip(1.0 - dist / float(highlight_mm), 0.0, 1.0) ** 0.5
        if tint.max() > 0:
            tint = tint / tint.max()          # the closest triangles reach FULL red

    img = Image.new("RGB", (W, W), BG)
    if shell is not None and shell[0] is not None:
        img = _draw_shell(img, shell[0], shell[1], R, project, front=False,
                          alpha=shell_alpha,
                          ao=shell[2] if len(shell) > 2 else None).convert("RGB")
    d = ImageDraw.Draw(img)
    order = np.argsort(depth)                        # far to near
    base = np.array(VESSEL_DEEP, float)
    top = np.array(VESSEL, float)
    hot = np.array(RED, float)
    for i in order:
        c = base + (top - base) * shade[i]
        if tint is not None and tint[i] > 0.01:
            c = c + (hot * (0.85 + 0.15 * shade[i]) - c) * tint[i]
        d.polygon([tuple(pt) for pt in p[i]], fill=tuple(np.clip(c, 0, 255).astype(int)))
    if shell is not None and shell[0] is not None:
        # front half last, at lower alpha: the vessels then read as sitting INSIDE the shell
        img = _draw_shell(img, shell[0], shell[1], R, project, front=True,
                          alpha=max(28, int(shell_alpha * 0.42)),
                          ao=shell[2] if len(shell) > 2 else None).convert("RGB")
    return img, project, R, ctr, scale, W



def mark_lesion(img, project, world_pos, W, R, label=None):
    """Draw a red ring at the projected lesion so it is unmissable in every frame.

    NOTE the `@ R.T`: project() expects coordinates that have ALREADY been rotated into camera
    space -- it is used internally on the rotated vertex array. Passing raw world coordinates put
    the ring in empty space beside the vessel instead of on it.
    """
    d = ImageDraw.Draw(img)
    x, y = project(np.asarray(world_pos, float) @ R.T)[0]
    r = max(int(W * 0.030), 10)
    d.ellipse([x - r, y - r, x + r, y + r], outline=RED, width=max(3, int(W * 0.006)))
    r2 = int(r * 1.7)
    d.ellipse([x - r2, y - r2, x + r2, y + r2], outline=(255, 120, 120), width=max(2, int(W * 0.003)))
    if label:
        f = _font(max(18, int(W * 0.046)))
        tw, th = d.textbbox((0, 0), label, font=f)[2:]
        lx = min(max(8, x - tw / 2), W - tw - 8)
        ly = y + r2 + 6
        # Scrim, translucent: the label tracks the lesion around the rotation, so on some frames
        # it lands on the reference heart, where red on dark red is barely readable.
        pill = Image.new("L", img.size, 0)
        ImageDraw.Draw(pill).rounded_rectangle(
            [lx - 10, ly - 6, lx + tw + 10, ly + th + 8], radius=10, fill=190)
        img.paste((10, 12, 16), mask=pill)
        d.text((lx, ly), label, font=f, fill=RED)
    return img




# ---------------------------------------------------------------------------------------------
# Reference heart (BodyParts3D FMA7274, "wall of heart")
#
# Real anatomy, not a shape I invented: BodyParts3D is the Database Center for Life Science's
# anatomical model set, and FMA7274 is the myocardial wall -- 331,592 triangles, 115 x 110 x 105 mm.
# Licence: CC BY-SA 2.1 Japan (see data/cardiac_ct/bodyparts3d/LICENSE_content). It stays in data/
# and is not redistributed by this repo.
#
# IMPORTANT, and said on the panel: this is a REFERENCE heart shown for orientation. It is a
# different person's anatomy from the coronary tree, and the two are in unrelated coordinate frames
# (whole-body vs patient CT), so the fit below is a similarity transform for context -- NOT a
# registration. It must never be captioned as this patient's myocardium.
# ---------------------------------------------------------------------------------------------
BP3D_DIR = BASE / "data" / "cardiac_ct" / "bodyparts3d"
HEART_STL = BP3D_DIR / "FMA7274.stl"

# Great vessels, so the heart reads as a heart rather than a lump: the aortic root and arch
# emerging from the base, the venae cavae returning to the right atrium, the coronary sinus in the
# posterior AV groove. All BodyParts3D, same source and licence as the myocardium, and all fitted
# by the SAME transform as the heart so their relative anatomy stays correct.
# Myocardium ONLY. Adding the great vessels was tried and reverted: the venae cavae and the aortic
# arch are full-length torso vessels, so the combined bounding box became dominated by a metre of
# tubing, the fit shrank the heart to accommodate it, and the coronaries no longer read as lying on
# the epicardium. The extra STLs are kept on disk in case a future view frames the mediastinum
# rather than the heart, but the demo shell is the heart wall alone.
HEART_PARTS = [
    ("FMA7274.stl", "wall of heart"),
]
HEART_TINT = (188, 104, 100)   # myocardium: browner and less pink than a pure red


def load_reference_heart(coronary_tris, scale_pad=1.10):
    """Reference heart wall, fitted around the coronary tree by PCA + uniform scale."""
    tri_list, nrm_list, used = [], [], []
    for fn_, label in HEART_PARTS:
        fp = BP3D_DIR / fn_
        if not fp.exists():
            continue
        t, n = read_stl(fp)
        tri_list.append(t); nrm_list.append(n); used.append(label)
    if not tri_list:
        return None, None
    ht = np.concatenate(tri_list, 0)
    hn = np.concatenate(nrm_list, 0)
    load_reference_heart.parts_used = used
    H = ht.reshape(-1, 3)
    C = coronary_tris.reshape(-1, 3)

    def axes(X):
        Xc = X - X.mean(0)
        u = np.linalg.svd(Xc[np.random.default_rng(0).choice(len(Xc),
                          min(20000, len(Xc)), replace=False)], full_matrices=False)[2]
        return u

    Ah, Ac = axes(H), axes(C)
    Rot = Ac.T @ Ah                      # heart principal frame -> coronary principal frame
    if np.linalg.det(Rot) < 0:
        Ah = Ah.copy(); Ah[2] *= -1; Rot = Ac.T @ Ah

    Hc = H - H.mean(0)
    Hr = Hc @ Rot.T
    # scale so the heart comfortably contains the vessel tree
    diag = lambda X: float(np.linalg.norm(X.max(0) - X.min(0)))
    sc = (diag(C) / max(diag(Hr), 1e-6)) * scale_pad
    Hf = Hr * sc + C.mean(0)

    ht_f = Hf.reshape(-1, 3, 3)
    hn_f = (hn @ Rot.T)
    ln = np.linalg.norm(hn_f, axis=1, keepdims=True)
    hn_f = np.divide(hn_f, np.where(ln == 0, 1, ln))
    ao = ambient_occlusion(ht_f, hn_f)
    return ht_f, hn_f, ao


# ---------------------------------------------------------------------------------------------
# Epicardial envelope
#
# There is NO segmented myocardium in this repo -- CoronariesNC6 ships vessel surfaces only, and
# nothing else here contains a heart mesh. Rather than model a heart from imagination and caption
# it as anatomy (the exact failure this whole rewrite exists to correct), the shell below is
# DERIVED FROM THE PATIENT'S OWN VESSELS: the coronaries run on the epicardial surface, so the
# convex envelope of the tree approximates the outside of the heart. For case V1P2 that envelope
# measures ~178 cm3 at 89 x 65 x 79 mm, against ~250-350 cm3 and ~12 x 9 x 6 cm for an adult heart
# -- the right order, and smaller because the arteries do not wrap the full posterior surface.
#
# It is therefore a SCHEMATIC context shell, and every caption says so. It is not a segmentation,
# it is not diagnostic, and it must never be labelled "anatomically correct".
# ---------------------------------------------------------------------------------------------
SHELL = (86, 132, 150)


def epicardial_shell(tris, sample=20000, seed=0):
    """Convex envelope of the coronary tree -> (facet vertices, facet normals)."""
    try:
        from scipy.spatial import ConvexHull
    except Exception:
        return None, None
    P = tris.reshape(-1, 3)
    rng = np.random.default_rng(seed)
    sub = P[rng.choice(len(P), min(sample, len(P)), replace=False)]
    h = ConvexHull(sub)
    facets = sub[h.simplices]                       # (F, 3, 3)
    n = np.cross(facets[:, 1] - facets[:, 0], facets[:, 2] - facets[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, np.where(ln == 0, 1, ln))
    ctr = sub.mean(0)
    out = np.einsum("ij,ij->i", facets.mean(1) - ctr, n) < 0
    n[out] *= -1                                    # orient outward
    return facets, n



def ambient_occlusion(facets, nrm, radius_mm=9.0, samples=30000, seed=0):
    """Per-facet ambient occlusion, computed once from the geometry.

    Creases -- the atrioventricular sulcus, the interventricular groove, the beds the coronaries
    lie in -- are geometrically enclosed, so less ambient light reaches them. Without this every
    surface receives the same fill and the heart reads as a smooth blob; with it the grooves that
    actually define cardiac shape appear.

    Point-based approximation: for each facet, neighbours lying in the hemisphere its normal points
    into are what would block incoming light, weighted by 1/distance so close geometry counts most.
    AO depends only on shape, never on the camera, so it is computed ONCE and reused for all 30
    rotation frames -- otherwise this would be the most expensive thing in the render.

    Cost control: computed on a subsample and propagated to every facet by nearest neighbour. At
    330k facets the difference from a full evaluation is invisible at these output sizes.
    """
    try:
        from scipy.spatial import cKDTree
    except Exception:
        return np.ones(len(facets))
    c = facets.mean(1)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(c), min(samples, len(c)), replace=False)
    cs, ns = c[idx], nrm[idx]
    tree = cKDTree(c)
    nb = tree.query_ball_point(cs, r=radius_mm)
    occ = np.zeros(len(cs))
    for k, ids in enumerate(nb):
        if len(ids) < 2:
            continue
        d = c[ids] - cs[k]
        dist = np.linalg.norm(d, axis=1)
        m = dist > 1e-6
        if not m.any():
            continue
        d, dist = d[m], dist[m]
        cosang = (d @ ns[k]) / dist              # >0 == in front of this facet, i.e. blocking
        w = np.clip(cosang, 0, 1) * (1.0 - dist / radius_mm)
        occ[k] = w.sum()
    if occ.max() > 0:
        occ = occ / np.percentile(occ, 97)
    ao_s = np.clip(1.0 - 0.85 * np.clip(occ, 0, 1.4), 0.15, 1.0)
    # propagate to every facet
    back = cKDTree(cs)
    _, near = back.query(c, k=1)
    return ao_s[near]


def _draw_shell(img, facets, nrm, R, project, front, alpha, ao=None):
    """Paint one side of the myocardium so its FORM reads, not just its silhouette.

    A single flat alpha across every facet renders the heart as a featureless slab: the chamber
    bulges, the atrioventricular groove and the apex all wash out, because a surface angled away
    from the viewer is drawn exactly as opaquely as one facing it. Three things fix that, and all
    three are cheap:

      * Fresnel-weighted opacity -- facets seen edge-on are far more opaque than facets facing the
        camera. This is what makes a translucent surface read as a solid object: the rim and every
        fold darken while the flat front stays see-through, so curvature becomes visible.
      * A specular term, so the lit shoulders catch a highlight and the surface reads as wet tissue
        rather than matte plastic.
      * Depth cue -- facets further from the camera desaturate slightly, which separates the far
        wall from the near one instead of letting them merge into one flat mass.
    """
    v = facets.reshape(-1, 3) @ R.T
    fn = nrm @ R.T
    p = project(v).reshape(-1, 3, 2)
    depth = v.reshape(-1, 3, 3)[:, :, 2].mean(1)
    facing = fn[:, 2] > 0                            # +z is toward the camera
    keep = np.where(facing if front else ~facing)[0]
    if len(keep) == 0:
        return img
    layer = Image.new("RGBA", img.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)

    L = np.array([-0.45, 0.55, 0.75]); L /= np.linalg.norm(L)
    ndl = np.clip(fn @ L, 0, 1)
    sh = 0.32 + 0.68 * ndl
    # specular: half-vector between light and the viewer (+z)
    H = L + np.array([0.0, 0.0, 1.0]); H /= np.linalg.norm(H)
    spec = np.clip(fn @ H, 0, 1) ** 22
    # broad, weak second lobe: wet tissue has a soft sheen as well as a tight highlight
    sheen = np.clip(fn @ H, 0, 1) ** 4

    # Fresnel: 1 when edge-on, 0 when facing us
    fres = 1.0 - np.abs(np.clip(fn[:, 2], -1, 1))
    a_face = float(alpha)
    a_edge = min(255.0, alpha * 3.4)
    av = a_face + (a_edge - a_face) * (fres ** 1.5)

    # Subsurface scattering, approximated. Muscle is translucent: where the wall is thin or turned
    # away from the key light, light that entered elsewhere leaks back out warm and red. Two cheap
    # proxies for that -- back-lighting (surfaces facing AWAY from the key) and grazing angles --
    # do most of the perceptual work, and it is the single biggest step from "plastic" to "tissue".
    back = np.clip(-(fn @ L), 0, 1) ** 1.6
    sss = np.clip(0.55 * back + 0.45 * (fres ** 2.2), 0, 1)
    SSS_TINT = np.array([236.0, 92.0, 84.0])

    # Rim light from behind-left, so the silhouette separates from the black background.
    RL = np.array([0.55, 0.25, -0.8]); RL /= np.linalg.norm(RL)
    rim = (np.clip(fn @ RL, 0, 1) ** 3) * (fres ** 0.8)

    # Tissue mottling: myocardium is not one flat colour. A deterministic low-frequency variation
    # keyed to facet position breaks up the plastic uniformity without inventing structure.
    ctr_f = facets.mean(1)
    mot = (np.sin(ctr_f[:, 0] * 0.9) * np.sin(ctr_f[:, 1] * 0.7)
           * np.sin(ctr_f[:, 2] * 1.1))
    mot = 1.0 + 0.085 * mot
    aov = np.ones(len(facets)) if ao is None else ao

    dmin, dmax = float(depth[keep].min()), float(depth[keep].max())
    drange = max(dmax - dmin, 1e-6)
    base = np.array(HEART_TINT, float)

    for i in keep[np.argsort(depth[keep])]:
        far = 1.0 - (depth[i] - dmin) / drange           # 1 = nearest
        c = base * sh[i] * (0.72 + 0.28 * far) * mot[i] * aov[i]
        c = c + (SSS_TINT - c) * (0.42 * sss[i])         # warm bleed through thin tissue
        c = c + 235.0 * spec[i] + 46.0 * sheen[i] + 150.0 * rim[i] * np.array([1.0, 0.55, 0.5])
        d.polygon([tuple(pt) for pt in p[i]],
                  fill=tuple(np.clip(c, 0, 255).astype(int))
                       + (int(np.clip(av[i] * (1.9 - 0.9 * aov[i]), 0, 255)),))
    return Image.alpha_composite(img.convert("RGBA"), layer)


def _font(px):
    for f in ("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
              "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"):
        if pathlib.Path(f).exists():
            return ImageFont.truetype(f, px)
    return ImageFont.load_default()



def _fit(d, text, max_w, px, min_px=12):
    """Largest font <= px whose rendering of `text` fits max_w.

    Fixed font sizes plus fixed x-positions is why enlarging the type clipped the title, ran the
    subtitle off the canvas and collided the metric labels. Measure, then size.
    """
    size = int(px)
    while size > min_px:
        f = _font(size)
        if d.textbbox((0, 0), text, font=f)[2] <= max_w:
            return f
        size = int(size * 0.94)
    return _font(min_px)


def _metric_row(d, items, y_label, y_value, x0, x1, f_sm, f_md, min_gap=None):
    """Lay metric columns out from MEASURED widths so they cannot overlap.

    Measuring alone is not enough: when the type is enlarged the columns can want more room than
    the strip has, and a fixed minimum gap then leaves them technically non-overlapping but
    visually touching ("CAD-RADS Max stenosis Agatston"). So shrink the type until the columns
    plus a readable gutter genuinely fit, then lay out.
    """
    avail = x1 - x0
    min_gap = min_gap if min_gap is not None else int(26 * SS)
    lo, md = f_sm.size, f_md.size
    for _ in range(24):
        f_s, f_m = _font(lo), _font(md)
        widths = [max(d.textbbox((0, 0), k, font=f_s)[2],
                      d.textbbox((0, 0), v, font=f_m)[2]) for k, v, _ in items]
        if sum(widths) + min_gap * (len(items) - 1) <= avail or lo <= 14:
            break
        lo, md = int(lo * 0.94), int(md * 0.94)
    gap = max(min_gap, (avail - sum(widths)) / max(1, len(items) - 1))
    x = x0
    for (k, v, col), w in zip(items, widths):
        d.text((x, y_label), k, font=f_s, fill=MUTED)
        d.text((x, y_value), v, font=f_m, fill=col)
        x += w + gap


def annotate(img, project, verts, size, case, stenosis_pct=72,
             agatston=385, cad_rads="4A", shell_note=False,
             lesion_world=None, R=None):
    # NOTE: stenosis_pct and cad_rads are MEASURED (passed in from coronary_analysis.json).
    # agatston is NOT -- calcium scoring needs CT Hounsfield units, which a surface mesh does not
    # carry, so it stays representative and is rendered in muted type and marked "(repr.)".
    """Mark a stenosis and label the panel truthfully.

    The marker is deliberately NOT tied to a named vessel. CoronariesNC6 meshes carry no
    anatomical branch labels, so writing "LAD" beside an arbitrary branch would be the same
    class of error as captioning a vertebra "CTA Coronary" — an unverifiable anatomical claim
    presented as fact. Vessel-level findings belong in the case narrative, which is where the
    clinical story already states them.
    """
    d = ImageDraw.Draw(img)
    W = img.size[0]
    # ~1.65x: at a 460 px panel these land near 25 / 18 / 14 px on screen.
    # Sized for a ~460 px panel: roughly 34 / 25 / 19 px on screen.
    f_lg, f_md, f_sm = _font(int(76 * SS)), _font(int(56 * SS)), _font(int(43 * SS))

    # Mark the MEASURED lesion. This used to be a heuristic -- the median of the top 34 % of
    # projected vertices -- which put the ring "somewhere proximal" while render() tinted the
    # geometry at the real lesion and the rotating model ringed the real lesion. Three panels,
    # two different places, one of them decorative. The ring now comes from the same
    # position_mm in coronary_analysis.json as everything else, so all three agree.
    if lesion_world is not None and R is not None:
        mx, my = project(np.asarray(lesion_world, float) @ R.T)[0]
    else:
        pts = project(verts)
        y_thresh = np.percentile(pts[:, 1], 34)
        cand = pts[pts[:, 1] < y_thresh]
        mx, my = float(np.median(cand[:, 0])), float(np.median(cand[:, 1]))

    r = int(30 * SS)
    d.ellipse([mx - r, my - r, mx + r, my + r], outline=RED, width=int(5 * SS))
    # Callout is the measurement alone. The "illustrative marker" qualifier that used to sit under
    # it moved to the provenance subtitle below, which is where the panel's other caveats live --
    # the claim is unchanged, it is just made once instead of twice.
    l1 = f"{stenosis_pct}% stenosis"
    f1 = _fit(d, l1, W * 0.40, 66 * SS)
    tw = d.textbbox((0, 0), l1, font=f1)[2]
    th = int(62 * SS)
    # Keep the callout inside the frame: flip it to the left of the marker when it would
    # overflow, and clamp vertically under the title block.
    lx, ly = mx + r * 3.4, my - r * 2.9
    if lx + tw + 20 * SS > W:
        lx = mx - r * 3.4 - tw
    lx = max(24 * SS, min(lx, W - tw - 24 * SS))
    ly = max(int(268 * SS), ly)
    d.line([mx + (r * 0.9 if lx > mx else -r * 0.9), my - r * 0.9,
            lx + (0 if lx > mx else tw), ly + 18 * SS], fill=RED, width=int(4 * SS))
    d.rounded_rectangle([lx - 12 * SS, ly - 10 * SS, lx + tw + 14 * SS, ly + th],
                        radius=8 * SS, fill=(28, 12, 12), outline=RED, width=int(3 * SS))
    d.text((lx, ly), l1, font=f1, fill=RED)

    # Title + provenance. This is the part that keeps the panel honest.
    #
    # Wrapped onto two lines deliberately. On one line the title is width-bound, and _fit caps it
    # at 80 px on this canvas -- about 30 px on screen, which is why it read small next to
    # "Stenosis Detail" (a short string that could reach 111 px). Wrapping removes the width
    # constraint, so both panels can share TITLE_PX and actually look the same size.
    d.rectangle([0, 0, W, int(TITLE_PX * 2.34) + int(60 * SS)], fill=(9, 11, 15))
    y_t = 22 * SS
    for i, line in enumerate(("Coronary Artery Tree", "3D Segmentation")):
        d.text((36 * SS, y_t + i * int(TITLE_PX * 1.12)), line,
               font=_fit(d, line, W - 72 * SS, TITLE_PX), fill=NVIDIA_GREEN)
    y_sub = y_t + int(TITLE_PX * 2.34)
    marker_note = (" · marker at the measured lesion" if lesion_world is not None
                   else " · lesion marker is illustrative")
    sub = (f"CoronariesNC6 {case} · manual rater ground truth · surface mesh, not a CT slice"
           + marker_note
           + (" · reference heart BodyParts3D FMA7274" if shell_note else ""))
    d.text((36 * SS, y_sub), sub,
           font=_fit(d, sub, W - 72 * SS, 43 * SS), fill=MUTED)

    # Findings strip
    y0 = W - int(212 * SS)
    d.rectangle([0, y0, W, W], fill=(16, 19, 24))
    d.line([0, y0, W, y0], fill=NVIDIA_GREEN, width=int(3 * SS))
    _metric_row(d, [("CAD-RADS", cad_rads, AMBER),
                    ("Max stenosis", f"{stenosis_pct}%", RED),
                    ("Agatston (repr.)", str(agatston), MUTED)],
                y0 + 42 * SS, y0 + 106 * SS, 40 * SS, W - 40 * SS, f_sm, f_md)
    # No per-tile status badge here. The same statement is made once, above the panel grid,
    # where it covers all four panels; repeating it on every tile added visual noise without
    # adding a claim. If these renders are ever used outside that page, restore it here.
    return img.resize((size, size), Image.LANCZOS)



def render_detail(tris, nrm, lesion, size=1024, fov_mm=20.0):
    """Close-up of the measured lesion, with calipers showing lumen vs reference calibre.

    Everything drawn here comes from coronary_analysis.json: the position the camera centres on,
    the minimum lumen radius, the local reference radius and the resulting stenosis. Nothing is
    positioned or captioned by hand, so the close-up cannot drift away from the measurement the
    rest of the demo reports.
    """
    pos = lesion["position_mm"]
    img, project, R, ctr, scale, W = render(tris, nrm, size=size, margin=0.06,
                                            focus_world=pos, span_mm=fov_mm)
    d = ImageDraw.Draw(img)
    # ~1.65x -- see the matching note in annotate().
    f_lg, f_md, f_sm = _font(int(72 * SS)), _font(int(54 * SS)), _font(int(42 * SS))

    cx, cy = W / 2.0, W / 2.0
    px_per_mm = scale

    r_min = float(lesion["min_radius_mm"])
    r_ref = float(lesion["reference_radius_mm"])
    sten = float(lesion["diameter_stenosis_pct"])
    area = float(lesion["area_stenosis_pct"])

    # Ring at the true measured lumen calibre, and a dashed ring at the reference calibre.
    rp = max(4 * SS, r_min * px_per_mm)
    rr = max(rp + 6 * SS, r_ref * px_per_mm)
    d.ellipse([cx - rp, cy - rp, cx + rp, cy + rp], outline=RED, width=int(4 * SS))
    for k in range(0, 360, 14):                       # dashed reference ring
        a0, a1 = math.radians(k), math.radians(k + 7)
        d.arc([cx - rr, cy - rr, cx + rr, cy + rr], math.degrees(a0), math.degrees(a1),
              fill=AMBER, width=int(3 * SS))

    # Calipers
    d.line([cx - rp, cy, cx + rp, cy], fill=RED, width=int(3 * SS))
    d.line([cx - rr, cy + rr + 14 * SS, cx + rr, cy + rr + 14 * SS], fill=AMBER, width=int(3 * SS))
    # Caliper labels live in a legend block in the emptiest corner, with a leader to the rings,
    # rather than floating beside them. Beside the rings each label needed a dark pill to stay
    # readable over the near-white vessel, and once the field of view widened those pills became
    # bars lying across the artery -- obscuring the very geometry the panel exists to show.
    # Each row carries a sample of the line style it names -- solid for the measured lumen, dashed
    # for the reference calibre. Red and amber are a poor pair under the common red-green colour
    # vision deficiencies, so the key must not rely on hue alone to say which ring is which.
    lines = [(f"lumen {r_min * 2:.2f} mm", RED, False),
             (f"reference {r_ref * 2:.2f} mm", AMBER, True)]
    chip = int(52 * SS)
    lw = max(d.textbbox((0, 0), t, font=f_sm)[2] for t, _, _ in lines) + chip
    lh = int(58 * SS)
    # Pick the corner with the least geometry in it, measured rather than assumed, so the legend
    # cannot land on the vessel just because the lesion sits somewhere unusual this time.
    arr = np.asarray(img.convert("L"))
    pad, bw_, bh_ = int(30 * SS), lw + int(34 * SS), lh * len(lines) + int(24 * SS)
    # Bound the placement by what is actually free canvas. The first version used the full image
    # height and put the legend at y = W - pad - bh, i.e. underneath the findings strip that gets
    # drawn afterwards -- the legend vanished and only its leader line survived.
    top_free = int(190 * SS)                       # below the header scrim
    bot_free = W - int(212 * SS) - pad - bh_       # above the findings strip
    corners = {
        "bl": (pad, bot_free), "br": (W - pad - bw_, bot_free),
        "tl": (pad, top_free), "tr": (W - pad - bw_, top_free),
    }
    def _busy(xy):
        x0, y0 = xy
        return float(arr[int(y0):int(y0 + bh_), int(x0):int(x0 + bw_)].mean())
    lx0, ly0 = min(corners.values(), key=_busy)

    scrim = Image.new("L", img.size, 0)
    ImageDraw.Draw(scrim).rounded_rectangle([lx0, ly0, lx0 + bw_, ly0 + bh_],
                                            radius=12 * SS, fill=188)
    img.paste((10, 12, 16), mask=scrim)
    for i, (t, col, dashed) in enumerate(lines):
        ty = ly0 + 12 * SS + i * lh
        my = ty + int(22 * SS)
        x_a, x_b = lx0 + 17 * SS, lx0 + 17 * SS + int(chip * 0.72)
        if dashed:
            step = int(11 * SS)
            for xs in range(int(x_a), int(x_b), step * 2):
                d.line([xs, my, min(xs + step, x_b), my], fill=col, width=int(4 * SS))
        else:
            d.line([x_a, my, x_b, my], fill=col, width=int(5 * SS))
        d.text((lx0 + 17 * SS + chip, ty), t, font=f_sm, fill=col)
    # Leader from the legend to the measured rings.
    d.line([lx0 + bw_ / 2, ly0 + (0 if ly0 > cy else bh_), cx, cy + (rr if ly0 > cy else -rr)],
           fill=(120, 132, 146), width=int(2 * SS))

    # Scrim behind the header: the close-up fills the frame with a near-white vessel, and muted
    # grey type on top of it was unreadable. The tree panel does not need this -- its header sits
    # over background -- but at this field of view the vessel reaches the top edge.
    d.rectangle([0, 0, W, int(180 * SS)], fill=(8, 9, 12))
    t2 = "Stenosis Detail"
    d.text((34 * SS, 20 * SS), t2, font=_fit(d, t2, W - 68 * SS, TITLE_PX), fill=NVIDIA_GREEN)
    s2 = f"{fov_mm:.0f} mm field of view \u00b7 centred on the measured lumen"
    d.text((34 * SS, 118 * SS), s2, font=_fit(d, s2, W - 68 * SS, 42 * SS), fill=MUTED)

    y0 = W - int(212 * SS)
    d.rectangle([0, y0, W, W], fill=(16, 19, 24))
    d.line([0, y0, W, y0], fill=NVIDIA_GREEN, width=int(3 * SS))
    _metric_row(d, [("Diameter stenosis", f"{sten:.1f}%", RED),
                    ("Area stenosis", f"{area:.1f}%", RED),
                    ("Lumen / reference", f"{r_min*2:.2f} / {r_ref*2:.2f} mm", AMBER)],
                y0 + 42 * SS, y0 + 106 * SS, 40 * SS, W - 40 * SS, f_sm, f_md)
    # No per-tile status badge here. The same statement is made once, above the panel grid,
    # where it covers all four panels; repeating it on every tile added visual noise without
    # adding a claim. If these renders are ever used outside that page, restore it here.
    return img.resize((size, size), Image.LANCZOS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="TEV2P2")
    ap.add_argument("--mesh", default="manualGT_rater1.stl")
    ap.add_argument("--size", type=int, default=1024)
    ap.add_argument("--spin", action="store_true", help="also write a rotating GIF")
    # Shell is ON by default: the demo assets are expected to carry the heart, and a
    # regeneration that silently dropped it would leave the panels disagreeing with the
    # walkthrough. --no-shell is there for the vessel-only view.
    ap.add_argument("--no-shell", dest="shell", action="store_false",
                    help="render the vessels without the reference heart")
    ap.set_defaults(shell=True)
    ap.add_argument("--frames", type=int, default=36)
    a = ap.parse_args()

    src = MESH_DIR / a.case / a.mesh
    if not src.exists():
        sys.exit(f"missing mesh: {src}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tris, nrm = read_stl(src)
    print(f"  {src.relative_to(BASE)}: {len(tris)} triangles")
    shell = load_reference_heart(tris) if a.shell else None
    if a.shell:
        if shell and shell[0] is not None:
            print(f"  reference heart: {len(shell[0])} facets (BodyParts3D FMA7274, fitted for context)")
        else:
            print("  !! reference heart mesh not found — rendering without it")

    # Draw the MEASURED numbers, not literals. coronary_analysis.json is written by
    # scripts/precompute_coronary_analysis.py from src/coronary_geometry, so the panel and the
    # workflow can never quietly disagree with each other.
    import json as _json
    meas = {}
    mp = OUT_DIR / "coronary_analysis.json"
    if mp.exists():
        try:
            meas = _json.loads(mp.read_text())
        except Exception:
            meas = {}
    sten = float(meas.get("max_diameter_stenosis_pct", 72))
    # Prefer the full CAD-RADS 2.0 report ("4A/P3/HRP") over the bare stenosis category.
    grade = str(meas.get("cadrads_report", "") or f"CAD-RADS {meas.get('cad_rads', '4A')}")
    grade = grade.replace("CAD-RADS ", "")

    _hl = meas["lesions"][0]["position_mm"] if meas.get("lesions") else None
    img, project, R, ctr, scale, W = render(tris, nrm, size=a.size,
                                            highlight_world=_hl, highlight_mm=6.5,
                                            shell=shell)
    verts = tris.reshape(-1, 3) @ R.T
    out = annotate(img, project, tris.reshape(-1, 3) @ R.T, a.size, a.case,
                   stenosis_pct=round(sten, 1), cad_rads=grade,
                   shell_note=bool(a.shell), lesion_world=_hl, R=R)
    p = OUT_DIR / "coronary_tree_annotated.png"
    out.save(p)
    print(f"  wrote {p.relative_to(BASE)}")

    # clean version, no annotation, for slides that supply their own labels
    clean, *_ = render(tris, nrm, size=a.size)
    if meas.get("lesions"):
        det = render_detail(tris, nrm, meas["lesions"][0], size=a.size)
        pd_ = OUT_DIR / "coronary_stenosis_detail.png"
        det.save(pd_)
        print(f"  wrote {pd_.relative_to(BASE)}")

    pc = OUT_DIR / "coronary_tree_clean.png"
    clean.resize((a.size, a.size), Image.LANCZOS).save(pc)
    print(f"  wrote {pc.relative_to(BASE)}")

    if a.spin:
        frames = []
        hl = meas["lesions"][0]["position_mm"] if meas.get("lesions") else None
        lbl = f"{sten:.1f}% stenosis" if meas.get("lesions") else None
        for k in range(a.frames):
            f, proj, _R, _c, _s, Wf = render(tris, nrm, size=512,
                                             ry=0.55 + 2 * math.pi * k / a.frames,
                                             highlight_world=hl, highlight_mm=6.5,
                                             shell=shell)
            if hl is not None:
                f = mark_lesion(f, proj, hl, Wf, _R, label=lbl)
            frames.append(f.resize((512, 512), Image.LANCZOS))
        pg = OUT_DIR / "coronary_tree_spin.gif"
        # optimize=True lets Pillow drop the small red region into a neighbouring palette entry;
        # an explicit adaptive palette per frame keeps it.
        frames = [f.convert("RGB").quantize(colors=128, method=Image.MEDIANCUT) for f in frames]
        # 60 frames x 280 ms = 16.8 s per revolution. GIF delays are stored in centiseconds,
        # so only multiples of 10 ms exist -- 280 is exact, 275 would silently become 280.
        #
        # GIF delays are stored in CENTISECONDS, so only multiples of 10 ms exist -- 77 and 85 both
        # silently became 80. Slowing by delay alone would step visibly, so the frame COUNT rose to
        # 60 as well: 6 deg per frame reads as smooth motion rather than a stutter.
        frames[0].save(pg, save_all=True, append_images=frames[1:], duration=280, loop=0,
                       optimize=False, disposal=2)
        print(f"  wrote {pg.relative_to(BASE)} ({a.frames} frames)")


if __name__ == "__main__":
    main()
