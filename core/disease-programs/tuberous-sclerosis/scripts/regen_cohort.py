"""
Regenerate the synthetic cohort and prove orchestrator integration against the built
artifacts (PRD §10 AC-1). Run from the engine root:  python3 scripts/regen_cohort.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.cohort.build import build_cohort           # noqa: E402
from src.cohort.loader import featured_map, load_patient  # noqa: E402
from src.orchestrator.graph import Orchestrator      # noqa: E402
from src.orchestrator.state import SqliteEventStore  # noqa: E402


def main() -> None:
    manifest = build_cohort()
    print("Built cohort:", manifest["n_patients"], "patients")
    print("  composition:", manifest["composition"])
    print("  featured:", manifest["featured"])
    print("  cohort_hash:", manifest["cohort_hash"])

    # integration: run the engine against the built (loaded-from-disk) featured patients
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    for tag, pid in featured_map().items():
        cohort = load_patient(pid)
        orch.enroll(pid, cohort)
        briefing = orch.assemble_surface(pid, "briefing")
        print(f"\n[{tag}={pid}] briefing action_items:")
        print(json.dumps(briefing["action_items"], indent=2, default=str))


if __name__ == "__main__":
    main()
