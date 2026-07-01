"""Tests for the RAG retriever and the Therapeutics Strategist (PRD §3 FR-TR-*; §2.5.4)."""
from __future__ import annotations

import config.settings as cs
from src.agents.therapeutics_strategist.runner import SECTION_KEYS
from src.agents.therapeutics_strategist.trials import match_trials
from src.cohort.build import build_cohort
from src.cohort.loader import featured_map, load_patient
from src.orchestrator.graph import Orchestrator
from src.orchestrator.state import SqliteEventStore
from src.rag.retriever import get_retriever


def test_rag_returns_sourced_chunks():
    hits = get_retriever().retrieve("everolimus reduces seizures in TSC", k=3)
    assert hits and all(h.get("source_uri") for h in hits)
    assert all("score" in h for h in hits)


def test_trial_matcher_resolves_each_criterion():
    matches = match_trials({"molecular_dx_confirmed": True, "age": 18,
                            "refractory_epilepsy": True, "aml_max_cm": 4.0})
    verdicts = {m["nct"]: m["match"] for m in matches}
    # aml 4cm trial: dx confirmed + aml>=3 -> eligible
    assert verdicts["NCT-SYN-002"] == "eligible"
    # epilepsy trial: stable-regimen criterion unknown -> requires_clarification
    assert verdicts["NCT-SYN-001"] == "requires_clarification"


def test_therapeutics_six_section_brief_via_orchestrator(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    cs.settings.COHORT_DIR = tmp_path
    pid = featured_map(tmp_path)["C"]
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    orch.enroll(pid, load_patient(pid, tmp_path))
    th = orch.store.projection(pid)["therapeutics"]
    assert set(th["sections"]) == set(SECTION_KEYS)
    assert th["trial_matches"]
    assert th["sources"] and all(s.get("source_uri") for s in th["sources"])
    assert "decision-support" in th["framing"]
