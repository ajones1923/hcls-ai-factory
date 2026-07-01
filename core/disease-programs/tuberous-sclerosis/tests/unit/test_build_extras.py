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
    assert les["population_prior"]["population_informed"] is True
    assert "prior_slope" in les["population_prior"]
    # population prior borrows strength but doesn't break Patient B's documented 12-18mo crossing
    assert les["crosses_in_12_18mo_window"] is True
    assert "mixed-effects" in les["model"]


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
