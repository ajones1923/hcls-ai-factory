"""Trajectory population prior, BAMSurgeon tool-gate, and optional Postgres backend."""
from __future__ import annotations

import config.settings as cs
import pytest

from src.agents.base import AgentContext
from src.agents.trajectory_modeler.runner import TrajectoryModeler
from src.cohort.bamsurgeon.pipeline import available, ensure_tools
from src.cohort.build import build_cohort
from src.cohort.loader import featured_map, load_patient


def test_trajectory_population_informed_and_preserves_crossing(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    cs.settings.COHORT_DIR = tmp_path
    pid = featured_map(tmp_path)["B"]
    out = TrajectoryModeler().run(AgentContext(patient_id=pid, cohort=load_patient(pid, tmp_path)))
    les = out.data["lesions"][0]
    # This asserted `population_informed is True` and `"mixed-effects" in model` until the
    # 2026-11-04 claim-integrity review. Both were the defect, not the contract: the prior is
    # fitted to 24 simulated patients, so nothing was population-informed and naming a
    # mixed-effects population fit implied the TSC Alliance NHD. The assertions are corrected
    # rather than dropped -- what they were really guarding is the SHAPE of the payload and
    # Patient B's crossing, both of which still hold.
    assert les["population_prior"]["population_informed"] is False
    assert les["population_prior"]["prior_is_synthetic"] is True
    assert "prior_slope" in les["population_prior"]
    # the prior still borrows strength: Patient B's documented 12-18mo crossing is unchanged
    assert les["crosses_in_12_18mo_window"] is True
    assert "synthetic-prior shrinkage" in les["model"]


def test_bamsurgeon_gate_fails_clearly_without_tools():
    avail = available()
    assert set(avail) >= {"samtools", "bcftools", "bwa"}
    if not all(avail.values()):           # this host lacks the genomics tools
        with pytest.raises(RuntimeError, match="Faithful blinded calling requires"):
            ensure_tools()


def test_postgres_backend_optional():
    if not cs.settings.POSTGRES_DSN:
        pytest.skip("no TSC_POSTGRES_DSN configured (Postgres is the prod backend)")
    from src.orchestrator.events import Event, EventType
    from src.orchestrator.state import PostgresEventStore

    store = PostgresEventStore()
    eid = store.append(Event(patient_id="PGTEST", event_type=EventType.patient_enrolled))
    assert eid and store.events_for("PGTEST")
