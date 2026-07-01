"""
Evaluation harness (PRD §5 eval criteria; AC-2). Per-agent demonstration targets against
synthetic ground truth. These are demonstration-config targets, NOT clinical-validation
thresholds (real-data validation is institutional Phase-1).
"""
from __future__ import annotations

import config.settings as cs
import pytest

from src.cohort.build import build_cohort
from src.cohort.loader import featured_map, load_manifest, load_patient
from src.orchestrator.graph import Orchestrator
from src.orchestrator.state import SqliteEventStore


@pytest.fixture(scope="module")
def engine(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("cohort")
    build_cohort(out_dir=tmp, seed=42)
    cs.settings.COHORT_DIR = tmp
    manifest = load_manifest(tmp)
    orch = Orchestrator(store=SqliteEventStore(":memory:"))
    for p in manifest["patients"]:
        orch.enroll(p["patient_id"], load_patient(p["patient_id"], tmp))
    return orch, manifest, tmp


def _primary(orch, pid):
    return (orch.store.projection(pid)["variant_interp"] or {}).get("primary")


# AC: Variant Curator — recover featured mosaic + classify nulls correctly, no FP Benign
def test_variant_curator_targets(engine):
    orch, manifest, tmp = engine
    fmap = featured_map(tmp)
    a = _primary(orch, fmap["A"])
    assert a["recovered"] and a["acmg_classification"] == "Likely Pathogenic"
    assert _primary(orch, fmap["B"])["acmg_classification"] == "Pathogenic"
    for p in manifest["patients"]:
        vi = _primary(orch, p["patient_id"])
        if vi and vi["consequence"] in ("frameshift", "nonsense"):
            assert vi["acmg_classification"] in ("Pathogenic", "Likely Pathogenic")
            assert "Benign" not in vi["acmg_classification"]


def test_mosaic_recovery_at_or_above_5pct(engine):
    orch, manifest, tmp = engine
    above = [p for p in manifest["patients"] if p["zygosity"] == "mosaic" and p["vaf"] >= 0.05]
    assert above
    for p in above:
        assert _primary(orch, p["patient_id"])["recovered"] is True


# AC: Phenome Mapper — recall >= 90% on structured ground truth
def test_phenome_mapper_recall(engine):
    orch, manifest, tmp = engine
    pid = featured_map(tmp)["B"]
    want = {p["hpo_id"] for p in load_patient(pid, tmp)["phenotypes"]}
    got = {t["hpo_id"] for t in orch.store.projection(pid)["hpo_profile"]["hpo_terms"]}
    assert len(got & want) / len(want) >= 0.9


# AC: Trajectory — Patient B SEGA crosses the 12-18 month window
def test_trajectory_target(engine):
    orch, manifest, tmp = engine
    les = orch.store.projection(featured_map(tmp)["B"])["trajectory"]["lesions"][0]
    assert les["crosses_in_12_18mo_window"] is True


# AC: TAND — flags Patient B academic; no spurious flag on a quiet patient
def test_tand_targets(engine):
    orch, manifest, tmp = engine
    assert "academic" in orch.store.projection(featured_map(tmp)["B"])["tand_briefing"]["flagged_clusters"]
    assert orch.store.projection(featured_map(tmp)["A"])["tand_briefing"]["flagged_clusters"] == []


# AC: Therapeutics — six-section, source-attributed brief
def test_therapeutics_target(engine):
    orch, manifest, tmp = engine
    th = orch.store.projection(featured_map(tmp)["C"])["therapeutics"]
    assert len(th["sections"]) == 6 and th["sources"]


# AC: every agent output carries provenance with latency (NFR-R-1)
def test_all_outputs_have_provenance(engine):
    orch, manifest, tmp = engine
    for e in orch.store.events_for(featured_map(tmp)["B"]):
        if e.event_type.value in ("phenome_mapped", "variant_curated", "trajectory_forecast",
                                  "tand_surveyed", "therapeutics_briefed"):
            recs = e.provenance["records"]
            assert recs and all("latency_ms" in r for r in recs)
