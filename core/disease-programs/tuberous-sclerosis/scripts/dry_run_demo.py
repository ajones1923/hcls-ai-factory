"""
Three-act demonstration runbook driver (PRD §8; master paper §19). Rehearses the live
demo against the built synthetic cohort: Act One (Patient A mosaic recovery), Act Two
(Patient B dashboard + Patient C therapeutics), Act Three (infrastructure/cost). Failure
of any act is caught and reported (conservative handling), never a silent break.

Run from the engine root:  python3 scripts/dry_run_demo.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.cohort.build import build_cohort                  # noqa: E402
from src.cohort.loader import featured_map, load_manifest, load_patient  # noqa: E402
from src.orchestrator.graph import Orchestrator            # noqa: E402
from src.orchestrator.state import SqliteEventStore        # noqa: E402


def _rule(title: str) -> None:
    print("\n" + "=" * 72 + f"\n{title}\n" + "=" * 72)


def main() -> None:
    if not (Path("data/cohort/manifest.json").exists()):
        build_cohort()
    manifest = load_manifest()
    fmap = featured_map()
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    for p in manifest["patients"]:
        orch.enroll(p["patient_id"], load_patient(p["patient_id"]))

    print("SYNTHETIC demonstration — not real patients. Decision support; clinician review required.")
    print(f"Cohort: {manifest['n_patients']} patients · hash {manifest['cohort_hash'][:23]}…")

    # ── Act One — mosaic recovery (Patient A) ──────────────────────────────
    try:
        _rule("ACT ONE — Diagnostic Recovery (Patient A)")
        a = orch.store.projection(fmap["A"])["variant_interp"]["primary"]
        print(f"  {a['gene']} {a['hgvsc']}  VAF={a['vaf']}  mosaic={a['mosaic']} recovered={a['recovered']}")
        print(f"  ACMG: {a['acmg_classification']}  via {a['acmg_rule']}")
        print(f"  criteria: {[c['code'] + '/' + c['bucket'] for c in a['acmg_criteria']]}")
        print(f"  ddPCR recommended: {a['ddpcr_recommended']} · review: pending_molecular_geneticist")
        print("  briefing:", json.dumps(orch.assemble_surface(fmap["A"], "briefing")["action_items"]))
    except Exception as exc:  # noqa: BLE001
        print(f"  [act one issue — open audit log, switch to pre-computed] {exc!r}")

    # ── Act Two — longitudinal patient (B) + therapeutics (C) ──────────────
    try:
        _rule("ACT TWO — Longitudinal Dashboard (Patient B) + Therapeutics (Patient C)")
        db = orch.assemble_surface(fmap["B"], "dashboard")
        traj = db["quadrants"]["trajectory"]["lesions"][0]
        print(f"  [B] SEGA {traj['observed_cm']} cm → crosses threshold in "
              f"~{traj['months_to_threshold']} mo (window={traj['crosses_in_12_18mo_window']})")
        tand = db["quadrants"]["tand_and_therapeutics"]["tand"]
        print(f"  [B] TAND flagged: {tand['flagged_clusters']} (surfaced as {tand['surfaced_as']})")
        th = orch.store.projection(fmap["C"])["therapeutics"]
        print(f"  [C] therapeutics: {len(th['sections'])} sections, "
              f"{len(th['trial_matches'])} trial match(es), {len(th['sources'])} cited sources")
        print("  [C] trial matches:", json.dumps([{m['nct']: m['match']} for m in th['trial_matches']]))
    except Exception as exc:  # noqa: BLE001
        print(f"  [act two issue] {exc!r}")

    # ── Act Three — infrastructure & cost ──────────────────────────────────
    _rule("ACT THREE — Infrastructure & Scale")
    n_events = sum(len(orch.store.events_for(p["patient_id"])) for p in manifest["patients"])
    print(f"  {manifest['n_patients']} patients · {n_events} audited events · runs on one DGX Spark")
    print("  LLM tiers via Claude API (Haiku/Sonnet/Opus); local Llama fallback; RunPod burst for cohort build.")
    print("  Apache 2.0 · synthetic-data demo · Epic/LIMS/imaging-AI = institutional Phase-1 (not built).")


if __name__ == "__main__":
    main()
