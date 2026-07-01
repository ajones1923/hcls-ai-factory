"""
Cohort / population endpoint (PRD §2.4, §11 — "the engine runs the whole panel, not a
demo of three"). Aggregates every patient's projection into a population view: genotype
and classification distributions, the count of mosaic variants recovered below the
standard calling threshold, lesions at/approaching an intervention threshold, and TAND
flags — plus a per-patient row table. This is what makes "50 patients on one Spark"
legible at a glance and lets a reviewer drill into any single case.
"""
from __future__ import annotations

from collections import Counter

from fastapi import APIRouter, Request

router = APIRouter(prefix="/cohort", tags=["cohort"])


def _summarize(eng, manifest, featured: dict) -> dict:
    rev = {v: k for k, v in (featured or {}).items()}
    rows = []
    cls_dist: Counter = Counter()
    gene_dist: Counter = Counter()
    zyg_dist: Counter = Counter()
    recovered = 0
    lesions_threshold = 0
    lesions_window = 0
    tand_flagged = 0
    overdue_total = 0

    for p in manifest["patients"]:
        pid = p["patient_id"]
        proj = eng.store.projection(pid)
        vi = (proj.get("variant_interp") or {}).get("primary") or {}
        tr = proj.get("trajectory") or {}
        td = proj.get("tand_briefing") or {}
        prof = proj.get("hpo_profile") or {}

        cls = vi.get("acmg_classification") or "—"
        gene = vi.get("gene") or p.get("gene") or "NMI"
        zyg = "mosaic" if vi.get("mosaic") else (p.get("zygosity") or "germline")
        cls_dist[cls] += 1
        gene_dist[gene] += 1
        zyg_dist[zyg] += 1
        if vi.get("recovered"):
            recovered += 1
        les_flags = []
        for l in tr.get("lesions", []):
            if l.get("at_or_above_threshold"):
                lesions_threshold += 1
                les_flags.append(f"{l['lesion']}≥{l.get('threshold_cm')}cm")
            elif l.get("crosses_in_12_18mo_window"):
                lesions_window += 1
                les_flags.append(f"{l['lesion']}→thr {l.get('months_to_threshold')}mo")
        flagged = td.get("flagged_clusters") or []
        if flagged:
            tand_flagged += 1
        overdue = [g for g in prof.get("surveillance_gaps", []) if g.get("status") == "overdue"]
        overdue_total += len(overdue)

        rows.append({
            "patient_id": pid,
            "featured": rev.get(pid),
            "gene": gene,
            "zygosity": zyg,
            "vaf": vi.get("vaf", p.get("vaf")),
            "variant": (vi.get("hgvsc") or "") + (f" {vi.get('hgvsp')}" if vi.get("hgvsp") else ""),
            "classification": cls,
            "recovered": bool(vi.get("recovered")),
            "lesion_flags": les_flags,
            "tand": flagged,
            "overdue": [g["type"] for g in overdue],
        })

    # featured first, then by classification severity, then id — stable, demo-friendly order
    sev = {"Pathogenic": 0, "Likely Pathogenic": 1, "Variant of Uncertain Significance": 2}
    rows.sort(key=lambda r: (r["featured"] is None, sev.get(r["classification"], 9), r["patient_id"]))

    return {
        "n_patients": manifest["n_patients"],
        "distributions": {
            "classification": dict(cls_dist),
            "gene": dict(gene_dist),
            "zygosity": dict(zyg_dist),
        },
        "highlights": {
            "mosaic_recovered": recovered,
            "lesions_at_threshold": lesions_threshold,
            "lesions_crossing_window": lesions_window,
            "tand_flagged": tand_flagged,
            "surveillance_overdue": overdue_total,
        },
        "patients": rows,
    }


@router.get("")
def cohort(request: Request):
    app = request.app
    return _summarize(app.state.engine, app.state.manifest, app.state.featured)
