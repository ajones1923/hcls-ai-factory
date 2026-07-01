"""
projection -> SceneSpec for the lesion-trajectory twin (FR-VZ-2..9, FR-VZ-16).

Reads the engine's trajectory quantities and authors, per lesion, a keyframe at each
observed month (radius = measured size, no uncertainty) and each forecast horizon (radius =
forecast mean, with the 50%/90% envelope radii set to the engine's prediction-interval upper
bounds). The envelope IS the interval — the render is exactly as crisp as the engine is
certain. The SYNTHETIC watermark and per-lesion provenance travel into the spec.
"""
from __future__ import annotations

import math

from config.settings import settings
from src.viz.atlas import HPO_ORGAN, ORGANS, anchor_for
from src.viz.palette import STATE_RGB, class_rgb, state_for_probability, state_for_value
from src.viz.scene_spec import (
    AtlasScene, BodyFigure, Cell, Keyframe, LesionScene, MosaicScene, OrganMark,
    PopulationScene, SceneSpec,
)

_SCALE = 1.0   # cm -> scene unit (stylized 1:1; lesion size maps to sphere radius)


def _lesion_scene(q: dict, provenance: dict) -> LesionScene:
    direction = q.get("direction", "increasing")
    threshold = float(q["threshold_cm"])
    a = anchor_for(q.get("location", ""))
    kfs: list[Keyframe] = []

    # observed history: radius = measured value, no uncertainty shell (env = radius)
    for month, value in zip(q["observed_months"], q["observed_cm"]):
        st = state_for_value(value, threshold, direction)
        r = round(value * _SCALE, 4)
        kfs.append(Keyframe(time_code=float(month), month=float(month), radius=r, state=st,
                            color=STATE_RGB[st], env50_radius=r, env90_radius=r, forecast=False))

    # forecast horizons: radius = mean; envelope radii = the prediction-interval UPPER bounds
    cp = q.get("crossing_probability", {})
    for h in (6, 12, 18):
        f = (q.get("forecast") or {}).get(f"m{h}")
        if not f:
            continue
        mean = f["mean_cm"]
        st = state_for_probability(cp.get(f"m{h}", 0.0))
        # for a decreasing quantity the "outer" interval bound toward threshold is the LOWER bound
        pi50 = f["pi50"][1] if direction == "increasing" else f["pi50"][0]
        pi90 = f["pi90"][1] if direction == "increasing" else f["pi90"][0]
        kfs.append(Keyframe(time_code=float(h), month=float(h), radius=round(mean * _SCALE, 4),
                            state=st, color=STATE_RGB[st],
                            env50_radius=round(pi50 * _SCALE, 4),
                            env90_radius=round(pi90 * _SCALE, 4), forecast=True))

    return LesionScene(
        label=q.get("lesion", q.get("quantity", "lesion")), location=q.get("location", ""),
        unit=q.get("unit", "cm"), threshold=threshold, anchor=a["pos"],
        context=a["context"], context_radius=a["context_radius"],
        threshold_radius=round(threshold * _SCALE, 4),
        crossing_grade=q.get("crossing_grade", ""), keyframes=kfs, provenance=provenance,
    )


def build_lesion_scene(patient_id: str, projection: dict, only_lesions: bool = True) -> SceneSpec:
    """Author the lesion-trajectory SceneSpec for one patient from the projection."""
    traj = projection.get("trajectory") or {}
    quantities = traj.get("lesions") if only_lesions else (traj.get("quantities") or traj.get("lesions"))
    quantities = quantities or []
    prov_records = []
    for e_section in ("trajectory",):
        st = (projection.get("staleness") or {}).get(e_section, {})
        if st:
            prov_records.append({"section": e_section, **st})

    lesions = [_lesion_scene(q, {"section": "trajectory"}) for q in quantities if q.get("observed_cm")]
    tcs = [k.time_code for le in lesions for k in le.keyframes] or [-24.0, 18.0]
    return SceneSpec(
        patient_id=patient_id, scene_kind="lesion_trajectory", watermark=settings.WATERMARK,
        scale=_SCALE, start_tc=min(tcs), end_tc=max(tcs), time_codes_per_second=6.0,
        lesions=lesions, provenance=prov_records,
    )


_GRID = 12   # 12x12 = 144 cells; resolution at which the VAF fraction is rendered


def build_mosaic_scene(patient_id: str, projection: dict, grid: int = _GRID) -> MosaicScene:
    """Author the mosaic 'powers-of-ten' SceneSpec: a cellular field in which exactly the
    recovered VAF fraction of cells carries the variant (FR-VZ-17). Variant-carrying cells
    are spread deterministically so the rendered fraction equals the engine's VAF."""
    p = (projection.get("variant_interp") or {}).get("primary") or {}
    vaf = float(p.get("vaf") or 0.0)
    n = grid * grid
    n_var = round(vaf * n)
    # evenly distribute the variant-carrying cells across the field (deterministic)
    var_idx = {round(i * n / n_var) for i in range(n_var)} if n_var else set()
    cells = []
    for i in range(n):
        row, col = divmod(i, grid)
        cells.append(Cell(pos=(float(col - grid / 2), float(row - grid / 2), 0.0),
                          variant=i in var_idx))
    return MosaicScene(
        patient_id=patient_id, watermark=settings.WATERMARK,
        gene=p.get("gene") or "—",
        variant=f"{p.get('hgvsc') or ''} {p.get('hgvsp') or ''}".strip() or "—",
        vaf=vaf, classification=p.get("acmg_classification") or "—",
        acmg_rule=p.get("acmg_rule") or "", recovered=bool(p.get("recovered")),
        criteria=[{"code": c.get("code"), "bucket": c.get("bucket"), "rationale": c.get("rationale")}
                  for c in p.get("acmg_criteria", [])],
        n_cells=n, n_variant=sum(1 for c in cells if c.variant), cells=cells,
    )


# ── whole-child organ atlas + population (Scene 3) ──────────────────────────
def _figure_for(patient_id: str, featured, projection: dict, pos: tuple) -> BodyFigure:
    vi = (projection.get("variant_interp") or {}).get("primary") or {}
    cls = vi.get("acmg_classification") or "No variant identified"
    prof = projection.get("hpo_profile") or {}
    involved = {HPO_ORGAN.get(t["hpo_id"]) for t in prof.get("hpo_terms", [])}
    involved.discard(None)
    organs = [OrganMark(organ=name, pos=spec["pos"], rgb=spec["rgb"], involved=name in involved)
              for name, spec in ORGANS.items()]
    overdue = [g["type"] for g in prof.get("surveillance_gaps", []) if g.get("status") == "overdue"]
    return BodyFigure(
        patient_id=patient_id, featured=featured, pos=pos, body_rgb=class_rgb(cls),
        classification=cls, recovered=bool(vi.get("recovered")), organs=organs,
        overdue=overdue, burden=len(involved),
    )


def build_atlas_scene(patient_id: str, projection: dict) -> AtlasScene:
    """Per-patient whole-child organ atlas (FR-VZ-18): organs lit by the phenome profile."""
    return AtlasScene(patient_id=patient_id, watermark=settings.WATERMARK,
                      figure=_figure_for(patient_id, None, projection, (0.0, 0.0, 0.0)))


def build_population_scene(patients: list, spacing: float = 14.0) -> PopulationScene:
    """Cohort population array (FR-VZ-14/18): a grid of whole-child figures coloured by ACMG
    class, organs lit by phenome, mosaic recoveries haloed. `patients` = [(pid, featured, proj)]."""
    n = len(patients)
    cols = max(1, round(math.sqrt(n)))
    figures, cls_dist, recovered = [], {}, 0
    for i, (pid, featured, proj) in enumerate(patients):
        r, c = divmod(i, cols)
        fig = _figure_for(pid, featured, proj, (c * spacing, -r * spacing, 0.0))
        figures.append(fig)
        cls_dist[fig.classification] = cls_dist.get(fig.classification, 0) + 1
        recovered += int(fig.recovered)
    return PopulationScene(
        watermark=settings.WATERMARK, n_patients=n, n_recovered=recovered, cols=cols,
        spacing=spacing, figures=figures,
        distributions={"classification": cls_dist},
    )
