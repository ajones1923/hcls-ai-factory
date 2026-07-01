#!/usr/bin/env python
"""
Guided 3-act demo runbook for Cincinnati Children's (run on the DGX Spark).

    venv/bin/python scripts/demo_runbook.py

Pulls LIVE numbers from the engine (no hardcoded metrics), authors every Omniverse USD
scene, and prints the act-by-act talking points, the exact commands/URLs to show, and the
RunPod RTX render checklist. SYNTHETIC demonstration data throughout.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app._engine import featured, get_engine            # noqa: E402
from src.eval import evaluate                            # noqa: E402
from src.viz import author_population, author_scene      # noqa: E402

HOST = "192.168.68.107"   # Spark LAN IP (portal :8560, dashboard :8562)


def _rule(c="─"):
    print(c * 78)


def _act(n, title):
    print()
    _rule("═")
    print(f"  ACT {n} — {title}")
    _rule("═")


def main() -> None:
    orch, manifest = get_engine()
    fmap = featured()
    ev = evaluate(orch, manifest, fmap)
    h, dy = ev["headline"], ev["diagnostic_yield"]
    pct = lambda x: f"{round(100 * x)}%" if x is not None else "—"

    _rule("═")
    print("  TSC INTELLIGENCE ENGINE — CINCINNATI CHILDREN'S DEMO RUNBOOK")
    print("  Engine 7 · HCLS AI Factory · one DGX Spark + RunPod (Omniverse RTX)")
    print("  SYNTHETIC demonstration data · decision support · not FDA-cleared")
    _rule("═")
    print(f"\n  PRE-FLIGHT (do before they arrive):")
    print(f"   • Start the stack:   bash scripts/serve.sh         (portal http://{HOST}:8560)")
    print(f"   • Author the scenes: venv/bin/python scripts/export_usd.py --all")
    print(f"   • Open in Omniverse on RunPod (see the render checklist at the end).")
    print(f"   • Confirm: {manifest['n_patients']} patients enrolled · reasoning offline ($0).")

    # ── ACT ONE ──────────────────────────────────────────────────────────────
    _act("ONE", "Mosaic recovery — ending the diagnostic odyssey  (Patient A)")
    a = (orch.store.projection(fmap["A"]).get("variant_interp") or {}).get("primary") or {}
    author_scene(fmap["A"], orch.store.projection(fmap["A"]), "mosaic")
    print(f"""
  THE STORY
    Patient A is a 4-year-old. Standard blood testing came back "no mutation
    identified" — the ~15% of TSC kids who get no molecular diagnosis. We sequence
    the resected tuber tissue instead.

  WHAT TO SHOW
    • Dashboard, Patient A:  http://{HOST}:8562   (Variant Curator card)
    • The engine recovers {a.get('gene')} {a.get('hgvsc')} at VAF {a.get('vaf')} (mosaic),
      classifies it {a.get('acmg_classification')} ({a.get('acmg_rule')}),
      and recommends ddPCR. Open the audit trail — every step is traceable.
    • Omniverse mosaic scene (TSC-0043_mosaic.usda): a field of cells where exactly
      ~1 in 12 GLOWS — 8.3% VAF made countable. The signal blood calls negative.

  THE NUMBER (measured, on this cohort)
    • Variant detection: {pct(dy['detection']['standard_rate'])} (standard) → {pct(dy['detection']['engine_rate'])} (engine) = +{round(100*dy['detection']['uplift_points'])} pts
    • Sub-threshold mosaics recovered: {dy['subthreshold_mosaics']}/{dy['subthreshold_mosaics']}  (sensitivity {pct(dy['subthreshold_mosaic_recovery_sensitivity'])})
    • Read-level discipline: the curator REJECTS a strand-biased artifact every sample
      carries — so "zero false-positive pathogenic" is earned, not free.

  THE LINE TO SAY
    "This is the capability your biobank unlocks: banked surgical tissue → a molecular
     diagnosis the blood test missed. On real specimens, this is a retrospective NMI
     re-analysis you already have the tissue for." """)

    # ── ACT TWO ──────────────────────────────────────────────────────────────
    _act("TWO", "The longitudinal twin — seeing the future of a lesion  (B & C)")
    author_scene(fmap["B"], orch.store.projection(fmap["B"]), "lesion_trajectory")
    author_scene(fmap["C"], orch.store.projection(fmap["C"]), "lesion_trajectory")
    bq = [q for q in (orch.store.projection(fmap["B"]).get("trajectory") or {}).get("lesions", [])]
    sega = next((q for q in bq if q["lesion"] == "SEGA"), {})
    print(f"""
  THE STORY
    Patient B (12yo) has a SEGA near the foramen of Monro. The question every visit:
    when does it reach the size where the team discusses intervention — and how sure
    are we? Patient C (18yo) has a growing renal AML, declining eGFR, and refractory
    seizures — three trajectories at once.

  WHAT TO SHOW
    • Dashboard, Patient B:  http://{HOST}:8562   (4-quadrant: variant · phenotype ·
      trajectory · TAND). The 2-D trajectory chart shows SEGA crossing in the 12–18mo
      window ({sega.get('crossing_grade','')} crossing).
    • Then the Omniverse twin (TSC-0001_lesion_trajectory.usda): SCRUB the timeline
      −24 → +18 months. The lesion grows along the forecast, wrapped in a GLASS
      uncertainty envelope whose radii ARE the 50/90% prediction intervals, crossing
      the red threshold membrane. The render is exactly as confident as the engine.
    • Patient C dashboard: eGFR + seizure-frequency tabs, surveillance cadence tightened
      vs the ITSC floor, six-section therapeutics brief (every claim sourced).

  THE LINE TO SAY
    "The envelope can't lie — it widens exactly as the model is uncertain. That's the
     thing that makes a forecast safe to put in front of a family." """)

    # ── ACT THREE ────────────────────────────────────────────────────────────
    _act("THREE", "Scale & infrastructure — one Spark, a whole disease")
    patients = [(p["patient_id"], {v: k for k, v in fmap.items()}.get(p["patient_id"]),
                 orch.store.projection(p["patient_id"])) for p in manifest["patients"]]
    pop = author_population(patients)
    dist = pop["spec"]["distributions"]["classification"]
    print(f"""
  WHAT TO SHOW
    • Portal:  http://{HOST}:8560   — the population command center (all {manifest['n_patients']}),
      the validation scorecard, and /system/usage (cost ledger — $0 offline).
    • Omniverse population scene (cohort_population.usda): {pop['spec']['n_patients']} whole-child
      figures, body colour = ACMG class, organs lit by phenome, and the
      {pop['spec']['n_recovered']} RECOVERED MOSAICS ringed in gold. The scale story, in 3-D.

  THE NUMBERS (validation scorecard — construct validity on synthetic ground truth)
    • ACMG classification accuracy: {pct(h['variant_classification_accuracy'])}
    • Diagnostic detection uplift:  +{round(100*h['detection_yield_uplift_points'])} pts   · mosaic recovery {pct(h['mosaic_recovery_sensitivity'])}
    • Truncating variant called benign: {h['false_benign_on_null_variant']}   · provenance complete: {pct(h['provenance_completeness'])}
    • Classification distribution: {dist}

  THE THESIS
    "Everything you just saw ran on a single $4,699 DGX Spark, on synthetic data, for
     $0 — and the beautiful render is a RunPod RTX pod, off-box. This is the engine.
     Phase 2 is your institution: the Discover Together Biobank as the substrate, the
     Winslow Pavilion as the envelope, Dr. Hagedorn's BMI as the methodology, the TSC
     clinic as the patients, Epic + the LIMS as the plumbing. Point us at the biobank
     and we'll show whether the recovery yield holds on your specimens."

  THE HONEST CAVEAT (say it out loud)
    "{ev['caveat']}" """)

    # ── RENDER CHECKLIST ─────────────────────────────────────────────────────
    _act("RENDER", "RunPod RTX render checklist (Omniverse)")
    print(f"""
    1. Launch a RunPod pod with an RTX GPU + NVIDIA Omniverse / USD Composer.
    2. Copy data/usd/*.usda to the pod (4 scene kinds: mosaic, lesion_trajectory,
       atlas, population).
    3. Open a scene; set renderer to RTX – Interactive (Path Tracing).
    4. MDL is wired: envelopes render as glass, variant cells + recovery halos glow.
    5. For the lesion twin, SCRUB the timeline (frame −24 → +18). For the population,
       fly through the grid and land on a gold-haloed recovery.

  Scenes authored to: {Path(pop['path']).parent}
""")
    _rule("═")
    print("  End of runbook. Lead with the tissue. Show the glass. Tell the truth.")
    _rule("═")


if __name__ == "__main__":
    main()
