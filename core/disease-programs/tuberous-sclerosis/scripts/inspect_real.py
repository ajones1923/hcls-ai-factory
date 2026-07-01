"""
Capture live-Claude output across all LLM steps on the featured patients, for quality
enrichment. Bills the key (~20-25 calls). Requires the key in .env, TSC_OFFLINE unset.
  venv/bin/python scripts/inspect_real.py
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


def show(title: str, obj, n: int = 1600) -> None:
    print(f"\n========== {title} ==========")
    print(json.dumps(obj, indent=2, default=str)[:n])


def main() -> None:
    r = get_router()
    print("router.online:", r.online)
    if not r.online:
        print("OFFLINE — set TSC_ANTHROPIC_API_KEY in .env and do NOT set TSC_OFFLINE=1.")
        return
    if not Path("data/cohort/manifest.json").exists():
        build_cohort()
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    fmap = featured_map()
    for pid in fmap.values():
        orch.enroll(pid, load_patient(pid))
    A, B, C = fmap["A"], fmap["B"], fmap["C"]
    pa, pb, pc = (orch.store.projection(x) for x in (A, B, C))

    show("A · Variant Curator — synthesis_narrative (Opus)", pa["variant_interp"].get("synthesis_narrative"))
    show("A · Phenome Mapper — note-extracted HPO terms (Sonnet)",
         [t for t in pa["hpo_profile"]["hpo_terms"] if t.get("source") == "note"])
    show("B · TAND — flagged + discourse highlights (Sonnet/Opus)",
         {"flagged": pb["tand_briefing"]["flagged_clusters"],
          "highlights": pb["tand_briefing"]["discourse_highlights"]})
    show("C · Therapeutics — six sections (Opus)", pc["therapeutics"]["sections"])


if __name__ == "__main__":
    main()
