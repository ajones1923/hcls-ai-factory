"""Production-environment backends (PRD §2.2-2.5): LangGraph runtime, Redis ephemeral
state, MinIO artifact store. Each tested via its runs-here fallback / equivalence."""
from __future__ import annotations

import config.settings as cs
import pytest

from src.cohort.build import build_cohort
from src.cohort.loader import featured_map, load_patient
from src.orchestrator.ephemeral import InMemoryEphemeral
from src.orchestrator.graph import Orchestrator
from src.orchestrator.state import SqliteEventStore
from src.utils.artifacts import FsArtifactStore


def test_langgraph_runtime_equivalent_to_plain_dispatcher(tmp_path):
    pytest.importorskip("langgraph")
    from src.orchestrator.langgraph_graph import LangGraphRunner

    build_cohort(out_dir=tmp_path, seed=42)
    cs.settings.COHORT_DIR = tmp_path
    pid = featured_map(tmp_path)["B"]
    cohort = load_patient(pid, tmp_path)

    plain = Orchestrator(store=SqliteEventStore(":memory:"))
    plain.enroll(pid, cohort)
    lg = LangGraphRunner(Orchestrator(store=SqliteEventStore(":memory:")))
    lg.enroll(pid, cohort)

    pj, lj = plain.store.projection(pid), lg.orch.store.projection(pid)
    assert all(lj[s] for s in ("hpo_profile", "variant_interp", "trajectory",
                               "tand_briefing", "therapeutics"))
    assert (pj["variant_interp"]["primary"]["acmg_classification"]
            == lj["variant_interp"]["primary"]["acmg_classification"])
    assert (pj["trajectory"]["lesions"][0]["crosses_in_12_18mo_window"]
            == lj["trajectory"]["lesions"][0]["crosses_in_12_18mo_window"])
    assert pj["tand_briefing"]["flagged_clusters"] == lj["tand_briefing"]["flagged_clusters"]


def test_ephemeral_locks_inflight_cache_alerts():
    e = InMemoryEphemeral()
    assert e.acquire_lock("P", "variant_curator") is True
    assert e.acquire_lock("P", "variant_curator") is False           # already held
    e.set_inflight("P", "variant_curator", True)
    assert e.inflight("P") == ["variant_curator"]
    e.set_inflight("P", "variant_curator", False)
    assert e.inflight("P") == []
    e.cache_surface("briefing", "P", {"x": 1})
    assert e.get_surface("briefing", "P") == {"x": 1}
    assert e.alert_window_incr("doc") == 1 and e.alert_window_incr("doc") == 2


def test_artifact_store_roundtrip(tmp_path):
    a = FsArtifactStore(root=tmp_path)
    uri = a.put("tsc-reports", "P/variant_draft.md", b"draft report")
    assert a.get("tsc-reports", "P/variant_draft.md") == b"draft report"
    assert uri.startswith("file://")
