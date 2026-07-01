#!/usr/bin/env python
"""
Print the TSC Intelligence Engine validation scorecard against synthetic ground truth.

    venv/bin/python scripts/run_eval.py

Construct validity (synthetic cohort) — NOT prospective clinical accuracy. See the caveat.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app._engine import featured, get_engine  # noqa: E402
from src.eval import evaluate                  # noqa: E402


def main() -> None:
    orch, manifest = get_engine()
    r = evaluate(orch, manifest, featured())
    h, dy, mc = r["headline"], r["diagnostic_yield"], r["forecast_calibration"]
    det, mdx = dy["detection"], dy["molecular_diagnosis"]
    pct = lambda x: f"{100 * x:.0f}%" if x is not None else "—"

    print("\n  TSC INTELLIGENCE ENGINE — VALIDATION SCORECARD")
    print(f"  cohort: {r['cohort']}\n")
    print(f"  Variant classification accuracy ........ {pct(h['variant_classification_accuracy'])}"
          f"  ({r['variant_classification']['correct']}/{r['n_patients']})")
    print(f"  Truncating variant called benign ....... {h['false_benign_on_null_variant']}  (must be 0)")
    print(f"\n  DIAGNOSTIC YIELD  (standard VAF limit {dy['standard_vaf_limit']})")
    print(f"    Variant detection .... {pct(det['standard_rate'])} -> {pct(det['engine_rate'])}"
          f"   (+{round(100*det['uplift_points'])} pts)")
    print(f"    Molecular diagnosis .. {pct(mdx['standard_rate'])} -> {pct(mdx['engine_rate'])}"
          f"   (+{round(100*mdx['uplift_points'])} pts)")
    print(f"    Sub-threshold mosaics recovered ...... {dy['subthreshold_mosaics']}/"
          f"{dy['subthreshold_mosaics']}  (sensitivity {pct(dy['subthreshold_mosaic_recovery_sensitivity'])})")
    print(f"    Newly diagnosed (all true positive: {mdx['newly_diagnosed_all_true_positive']}): "
          f"{', '.join(mdx['newly_diagnosed'])}")
    print(f"\n  FORECAST CALIBRATION (leave-one-out probe, n={mc['n_held_out_points']})")
    print(f"    PI50 coverage {pct(mc['empirical_coverage_pi50'])} · PI90 coverage {pct(mc['empirical_coverage_pi90'])}")
    td = r["tand_and_discipline"]
    print(f"\n  TAND  detects academic cluster: {td['detects_under_recognized_academic_cluster']} · "
          f"no spurious flag: {td['no_spurious_flag_on_quiet_patient']} · "
          f"{td['mean_alerts_per_patient']} alerts/patient")
    print(f"  Provenance completeness ................ {pct(r['provenance']['completeness'])}"
          f"  ({r['provenance']['carry_latency_and_trace']}/{r['provenance']['model_backed_outputs']})")
    print(f"\n  ⚠ {r['caveat']}\n")


if __name__ == "__main__":
    main()
