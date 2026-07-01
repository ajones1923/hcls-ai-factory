"""
Real-reasoning demo on the THREE featured patients only (A/B/C) — the cheap way to
exercise the API key. ~20-25 billed Claude calls total (not the ~300 a full 50-patient
run would cost). Requires TSC_ANTHROPIC_API_KEY in .env and TSC_OFFLINE unset/0.

  venv/bin/python scripts/real_demo.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.cohort.build import build_cohort                  # noqa: E402
from src.cohort.loader import featured_map, load_patient   # noqa: E402
from src.orchestrator.graph import Orchestrator            # noqa: E402
from src.orchestrator.state import SqliteEventStore        # noqa: E402
from src.utils.model_router import get_router              # noqa: E402


def main() -> None:
    if not get_router().online:
        print("Router is OFFLINE. Set TSC_ANTHROPIC_API_KEY in .env and do NOT set TSC_OFFLINE=1.")
        return
    if not (Path("data/cohort/manifest.json").exists()):
        build_cohort()
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    fmap = featured_map()
    print("Running REAL Claude reasoning on featured patients A/B/C (this bills your key)...\n")
    for pid in fmap.values():
        orch.enroll(pid, load_patient(pid))

    A = fmap["A"]
    vi = orch.store.projection(A)["variant_interp"]["primary"]
    nar = orch.store.projection(A)["variant_interp"].get("synthesis_narrative")
    print(f"[A {A}] {vi['gene']} {vi['hgvsc']}  ->  {vi['acmg_classification']}  (validator)")
    print("  Opus ACMG narrative (live):")
    print("   ", json.dumps(nar, indent=2)[:700] if nar else "(none)")

    C = fmap["C"]
    th = orch.store.projection(C)["therapeutics"]
    print(f"\n[C {C}] therapeutics — six sections (live):")
    for k, v in th["sections"].items():
        print(f"   {k}: {str(v)[:140]}")
    print(f"\n  cited sources: {[s['source_uri'] for s in th['sources']]}")


if __name__ == "__main__":
    main()
