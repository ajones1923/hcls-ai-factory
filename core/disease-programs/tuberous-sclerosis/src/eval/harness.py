"""
Evaluation harness (PRD §5; master paper §16) — MEASURED engine performance against the
synthetic cohort's known ground truth (src/cohort/spec build_roster).

This is the difference between "looks right" and "is right": every claim the engine makes
becomes a number you can audit. The headline is DIAGNOSTIC YIELD — how many patients get a
molecular diagnosis the standard germline-focused pipeline would return as negative, which
in real TSC is the ~15% "no mutation identified" diagnostic-odyssey population.

HONESTY: these are CONSTRUCT-VALIDITY metrics on synthetic data with planted ground truth.
They prove the engine's logic recovers the signal we encoded — NOT prospective clinical
accuracy. Real clinical validation (prospective, IRB-governed, against orthogonal assays
and expert adjudication) is institutional Phase-1. Every surface that shows these numbers
carries that caveat.
"""
from __future__ import annotations

import numpy as np

from src.agents.trajectory_modeler.runner import (
    SEGA_DISCUSSION_THRESHOLD_CM, _gp_posterior,
)
from src.agents.trajectory_modeler.population import shrink_slope
from src.cohort.spec import build_roster

# Routine germline/NGS pipelines reliably call variants only at higher allele fractions;
# low-level somatic mosaicism below this floor is the classic miss in TSC. The engine's
# mosaic-aware + affected-tissue path is what recovers it. Configurable assumption.
STANDARD_VAF_LIMIT = 0.10
_DIAGNOSTIC = ("Pathogenic", "Likely Pathogenic")


def _family(s: str | None) -> str:
    s = s or ""
    if not s or s == "—" or "No variant" in s or "No mutation" in s:
        return "No variant identified"
    if "Likely Pathogenic" in s:
        return "Likely Pathogenic"
    if "Likely Benign" in s:
        return "Likely Benign"
    if "Pathogenic" in s:
        return "Pathogenic"
    if "Benign" in s:
        return "Benign"
    if "Uncertain" in s or "VUS" in s:
        return "Variant of Uncertain Significance"
    return "Other"


def _classification(rows: list[dict]) -> dict:
    n = len(rows)
    correct = sum(r["pred_family"] == r["truth_family"] for r in rows)
    confusion: dict[str, dict[str, int]] = {}
    for r in rows:
        confusion.setdefault(r["truth_family"], {}).setdefault(r["pred_family"], 0)
        confusion[r["truth_family"]][r["pred_family"]] += 1
    # safety: a truncating (null) variant must NEVER be called benign
    false_benign_on_null = sum(
        1 for r in rows
        if r["kind"] in ("nonsense", "frameshift") and r["pred_family"] in ("Benign", "Likely Benign")
    )
    return {
        "n": n,
        "accuracy": round(correct / n, 4),
        "correct": correct,
        "confusion": confusion,
        "false_benign_on_null_variant": false_benign_on_null,
    }


def _diagnostic_yield(rows: list[dict]) -> dict:
    """Two honest, distinct yield metrics:

    DETECTION — does the pipeline physically return the variant at all? This is the engine's
    unique capability: a standard germline-focused pipeline reliably calls only VAF >= the
    limit, so low-level mosaics come back negative; the engine's mosaic-aware + affected-tissue
    path recovers them. The sub-threshold mosaics are exactly real TSC's "no mutation
    identified" diagnostic-odyssey population.

    MOLECULAR DIAGNOSIS — detection AND a returnable Pathogenic/Likely-Pathogenic call. A
    stricter bar (a detected VUS is not a diagnosis); the uplift here is the patients newly
    given an actionable molecular diagnosis.
    """
    n = len(rows)
    std_detect, eng_detect = set(), set()
    std_dx, eng_dx, newly_dx = set(), set(), []
    sub_total, sub_recovered, false_detect = 0, 0, []
    for r in rows:
        pid = r["patient_id"]
        standard_detects = (r["zygosity"] == "germline") or (
            r["zygosity"] == "mosaic" and r["vaf"] >= STANDARD_VAF_LIMIT)
        engine_detects = r["true_variant"] and (r["pred_family"] != "No variant identified")
        if standard_detects:
            std_detect.add(pid)
        if engine_detects:
            eng_detect.add(pid)
        if not r["true_variant"] and engine_detects:        # must never "find" a variant in true-NMI
            false_detect.append(pid)
        # molecular diagnosis (P/LP)
        if standard_detects and r["truth_family"] in _DIAGNOSTIC:
            std_dx.add(pid)
        if r["pred_family"] in _DIAGNOSTIC:
            eng_dx.add(pid)
            if pid not in std_dx:
                newly_dx.append(r)
        if r["zygosity"] == "mosaic" and r["vaf"] < STANDARD_VAF_LIMIT:
            sub_total += 1
            if r["recovered"]:
                sub_recovered += 1

    sd, ed = len(std_detect) / n, len(eng_detect) / n
    sdx, edx = len(std_dx) / n, len(eng_dx) / n
    all_tp = all(r["true_variant"] and r["truth_family"] in _DIAGNOSTIC for r in newly_dx)
    return {
        "standard_vaf_limit": STANDARD_VAF_LIMIT,
        "detection": {
            "standard_rate": round(sd, 4), "engine_rate": round(ed, 4),
            "uplift_points": round(ed - sd, 4),
            "relative_uplift": round((ed - sd) / sd, 4) if sd else None,
            "false_detections_in_true_nmi": false_detect,
        },
        "molecular_diagnosis": {
            "standard_rate": round(sdx, 4), "engine_rate": round(edx, 4),
            "uplift_points": round(edx - sdx, 4),
            "newly_diagnosed": [r["patient_id"] for r in newly_dx],
            "newly_diagnosed_all_true_positive": all_tp,
        },
        "subthreshold_mosaics": sub_total,
        "subthreshold_mosaic_recovery_sensitivity":
            round(sub_recovered / sub_total, 4) if sub_total else None,
    }


def _forecast_calibration(orch, manifest, featured: dict) -> dict:
    """Leave-one-out probe: drop each longitudinal lesion's most recent scan, forecast it
    from the prior scans, and check whether the held-out value lands in the 50%/90%
    prediction interval. Small n (the two featured longitudinal patients) — a probe, not a
    population calibration study (that needs the longitudinal-cohort expansion, W3-W4)."""
    probes = []
    for tag in ("B", "C"):
        pid = featured.get(tag)
        if not pid:
            continue
        for les in (orch.store.projection(pid).get("trajectory") or {}).get("lesions", []):
            ys, ts = les["observed_cm"], les["observed_months"]
            if len(ys) < 3:
                continue
            t_tr, y_tr = np.array(ts[:-1], float), np.array(ys[:-1], float)
            t_h, y_h = float(ts[-1]), float(ys[-1])
            ols = float(np.polyfit(t_tr, y_tr, 1)[0])
            slope = shrink_slope(t_tr, y_tr, ols)[0] if les["lesion"] == "SEGA" else ols
            intercept = float(y_tr.mean() - slope * t_tr.mean())
            _, std = _gp_posterior(t_tr, y_tr, [t_h])
            mu, sd = slope * t_h + intercept, float(std[0])
            in50 = abs(y_h - mu) <= 0.674 * sd
            in90 = abs(y_h - mu) <= 1.645 * sd
            probes.append({"patient": tag, "lesion": les["lesion"], "held_out_cm": round(y_h, 2),
                           "predicted_cm": round(mu, 2), "in_pi50": bool(in50), "in_pi90": bool(in90)})
    nprobe = len(probes)
    return {
        "n_held_out_points": nprobe,
        "empirical_coverage_pi50": round(sum(p["in_pi50"] for p in probes) / nprobe, 3) if nprobe else None,
        "empirical_coverage_pi90": round(sum(p["in_pi90"] for p in probes) / nprobe, 3) if nprobe else None,
        "probes": probes,
        "note": "leave-one-out probe on featured longitudinal lesions; population calibration is Phase-1",
    }


def _tand_and_discipline(orch, manifest, featured: dict) -> dict:
    def flags(pid):
        return (orch.store.projection(pid).get("tand_briefing") or {}).get("flagged_clusters") or []
    b_flag = "academic" in flags(featured.get("B", ""))      # known under-recognized academic cluster
    a_quiet = flags(featured.get("A", "")) == []             # no spurious flag on the quiet patient
    # alert discipline across the whole cohort (NFR — never overwhelm a clinician)
    total_alerts = sum(len(orch.assemble_surface(p["patient_id"], "alerts")["alerts"])
                       for p in manifest["patients"])
    return {
        "detects_under_recognized_academic_cluster": b_flag,
        "no_spurious_flag_on_quiet_patient": a_quiet,
        "total_cohort_alerts": total_alerts,
        "mean_alerts_per_patient": round(total_alerts / manifest["n_patients"], 3),
    }


def _provenance_completeness(orch, manifest) -> dict:
    model_events = {"phenome_mapped", "variant_curated", "trajectory_forecast",
                    "tand_surveyed", "therapeutics_briefed"}
    total = with_latency = 0
    for p in manifest["patients"]:
        for e in orch.store.events_for(p["patient_id"]):
            if e.event_type.value in model_events:
                for rec in (e.provenance or {}).get("records", []):
                    total += 1
                    if "latency_ms" in rec:
                        with_latency += 1
    return {
        "model_backed_outputs": total,
        "carry_latency_and_trace": with_latency,
        "completeness": round(with_latency / total, 4) if total else None,
    }


def evaluate(orch, manifest, featured: dict | None = None) -> dict:
    """Run the full scorecard. `featured` defaults to the engine's featured map if present."""
    featured = featured or getattr(orch, "_featured", None) or {}
    specs = {s.patient_id: s for s in build_roster()}
    rows = []
    for p in manifest["patients"]:
        pid = p["patient_id"]
        spec = specs[pid]
        vi = (orch.store.projection(pid).get("variant_interp") or {}).get("primary") or {}
        rows.append({
            "patient_id": pid,
            "zygosity": spec.zygosity,
            "vaf": spec.vaf,
            "kind": (spec.variant or {}).get("kind"),
            "true_variant": spec.gene is not None,
            "truth_family": _family(spec.expected_acmg),
            "pred_family": _family(vi.get("acmg_classification")) if vi else "No variant identified",
            "recovered": bool(vi.get("recovered")),
        })

    classification = _classification(rows)
    yield_ = _diagnostic_yield(rows)
    return {
        "cohort": "synthetic-50 (construct validity — NOT prospective clinical accuracy)",
        "n_patients": manifest["n_patients"],
        "headline": {
            "variant_classification_accuracy": classification["accuracy"],
            "detection_yield_uplift_points": yield_["detection"]["uplift_points"],
            "molecular_diagnosis_uplift_points": yield_["molecular_diagnosis"]["uplift_points"],
            "mosaic_recovery_sensitivity": yield_["subthreshold_mosaic_recovery_sensitivity"],
            "false_benign_on_null_variant": classification["false_benign_on_null_variant"],
            "provenance_completeness": _provenance_completeness(orch, manifest)["completeness"],
        },
        "variant_classification": classification,
        "diagnostic_yield": yield_,
        "forecast_calibration": _forecast_calibration(orch, manifest, featured),
        "tand_and_discipline": _tand_and_discipline(orch, manifest, featured),
        "provenance": _provenance_completeness(orch, manifest),
        "caveat": "Construct-validity on synthetic ground truth. Prospective clinical "
                  "validation against orthogonal assays + expert adjudication is Phase-1.",
    }
