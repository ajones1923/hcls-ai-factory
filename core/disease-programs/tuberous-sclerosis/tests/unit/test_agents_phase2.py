"""Tests for Phenome Mapper, Trajectory Modeler (GP), and TAND Surveillance (PRD §3)."""
from __future__ import annotations

import config.settings as cs
from src.agents.base import AgentContext
from src.agents.phenome_mapper.runner import PhenomeMapper
from src.agents.tand_surveillance.runner import TANDSurveillance
from src.agents.trajectory_modeler.runner import TrajectoryModeler
from src.cohort.build import build_cohort
from src.cohort.loader import featured_map, load_patient


def _cohort_b(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    cs.settings.COHORT_DIR = tmp_path
    return load_patient(featured_map(tmp_path)["B"], tmp_path)


# ── Phenome Mapper ─────────────────────────────────────────────────────────
def test_phenome_mapper_profile_and_surveillance(tmp_path):
    cohort = _cohort_b(tmp_path)
    out = PhenomeMapper().run(AgentContext(patient_id="B", cohort=cohort))
    d = out.data
    assert d["n_terms"] >= 1 and d["hpo_terms"]
    gap_types = {g["type"] for g in d["surveillance_gaps"]}
    assert "brain MRI" in gap_types                     # SEGA/tuber triggers it
    assert "neuropsychiatric screen (TAND-L)" in gap_types  # recommended for all


# ── Trajectory Modeler (Gaussian-process intervals) ────────────────────────
def test_trajectory_gp_intervals_and_crossing(tmp_path):
    cohort = _cohort_b(tmp_path)
    out = TrajectoryModeler().run(AgentContext(patient_id="B", cohort=cohort))
    les = out.data["lesions"][0]
    assert les["crosses_in_12_18mo_window"] is True
    m12 = les["forecast"]["m12"]
    assert len(m12["pi50"]) == 2 and len(m12["pi90"]) == 2
    # 90% interval is wider than 50% interval (proper posterior)
    assert (m12["pi90"][1] - m12["pi90"][0]) > (m12["pi50"][1] - m12["pi50"][0])
    assert "Gaussian-process" in les["model"]


# ── TAND Surveillance (Marshall-Hagedorn taxonomy) ─────────────────────────
def test_tand_flags_academic_with_discourse_markers(tmp_path):
    cohort = _cohort_b(tmp_path)
    out = TANDSurveillance().run(AgentContext(patient_id="B", cohort=cohort))
    d = out.data
    assert "academic" in d["flagged_clusters"]
    assert d["surfaced_as"] == "pre_visit_briefing"     # never an alert
    # the "Mother reports ... will reassess next visit" note yields real markers
    markers = {m for h in d["discourse_highlights"] for m in h["markers"]}
    assert "third_party_attribution" in markers
    assert "deferral" in markers or "hedging" in markers
    assert d["taxonomy_version"].startswith("marshall-hagedorn")


def test_tand_no_spurious_flag_for_quiet_patient(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    cs.settings.COHORT_DIR = tmp_path
    # Patient A has no TAND signals/notes of concern -> no flagged clusters
    cohort = load_patient(featured_map(tmp_path)["A"], tmp_path)
    out = TANDSurveillance().run(AgentContext(patient_id="A", cohort=cohort))
    assert out.data["flagged_clusters"] == []
