"""Unit tests for the event-sourced core and the deterministic orchestrator (PRD §2.5-2.6)."""
from __future__ import annotations

from src.agents.base import Agent, AgentContext, AgentOutput
from src.orchestrator.events import AGENT_DEPENDS_ON, AGENT_EMITS, Event, EventType
from src.orchestrator.graph import Orchestrator
from src.orchestrator.policies import DISPATCH_ORDER, default_agents
from src.orchestrator.state import SqliteEventStore, materialize_projection
from src.utils.provenance import stable_hash


def fresh_orch():
    return Orchestrator(store=SqliteEventStore(":memory:"))


def cohort_b():
    return {
        "genomics": {"gene": "TSC2", "expected_acmg": "Pathogenic"},
        "sega_series_cm": [0.8, 1.1, 1.3], "sega_months": [-24, -12, -6],
        "tand_signals": [{"cluster": "academic", "source": f"n{i}"} for i in range(3)],
    }


def test_thirteen_event_types():
    assert len(list(EventType)) == 13


def test_provenance_hash_is_deterministic():
    a = stable_hash({"x": 1, "y": [2, 3]})
    b = stable_hash({"y": [2, 3], "x": 1})
    assert a == b and a.startswith("sha256:")


def test_sqlite_append_and_readback():
    store = SqliteEventStore(":memory:")
    eid = store.append(Event(patient_id="P1", event_type=EventType.patient_enrolled, payload={"k": 1}))
    assert eid == 1
    evs = store.events_for("P1")
    assert len(evs) == 1 and evs[0].event_type == EventType.patient_enrolled


def test_projection_folds_agent_outputs():
    evs = [
        Event(event_id=1, patient_id="P1", event_type=EventType.patient_enrolled),
        Event(event_id=2, patient_id="P1", event_type=EventType.phenome_mapped, payload={"hpo_terms": []}),
        Event(event_id=3, patient_id="P1", event_type=EventType.variant_curated, payload={"gene": "TSC2"}),
    ]
    proj = materialize_projection(evs)
    assert proj["variant_interp"] == {"gene": "TSC2"}
    assert proj["staleness"]["variant_interp"]["status"] == "ok"
    assert proj["last_event_id"] == 3


def test_enroll_runs_all_agents_in_order_no_staleness():
    orch = fresh_orch()
    orch.enroll("DEMO-B", cohort_b())
    evs = [e.event_type.value for e in orch.store.events_for("DEMO-B")]
    for emitted in AGENT_EMITS.values():
        assert emitted.value in evs, f"missing {emitted}"
    proj = orch.store.projection("DEMO-B")
    assert all(v["status"] == "ok" for v in proj["staleness"].values())


def test_trajectory_crosses_window_for_patient_b():
    orch = fresh_orch()
    orch.enroll("DEMO-B", cohort_b())
    traj = orch.store.projection("DEMO-B")["trajectory"]["lesions"][0]
    # 0.8 -> 1.1 -> 1.3 over 18 months; linear extrapolation crosses 1.5 cm in the 12-18mo window
    assert traj["crosses_in_12_18mo_window"] is True
    assert traj["slope_cm_per_month"] > 0


def test_surfaces_assemble():
    orch = fresh_orch()
    orch.enroll("DEMO-B", cohort_b())
    for kind in ("briefing", "dashboard", "alerts"):
        s = orch.assemble_surface("DEMO-B", kind)
        assert s["kind"] == kind and s["watermark"] == "SYNTHETIC"
    # Patient B's SEGA forecast should produce at least one action item / alert
    assert orch.assemble_surface("DEMO-B", "alerts")["alerts"]


def test_conservative_failure_marks_pending_not_silent():
    class Boom(Agent):
        name = "tand_surveillance"

        def run(self, ctx: AgentContext) -> AgentOutput:
            raise RuntimeError("boom")

    agents = default_agents()
    agents = [Boom() if a.name == "tand_surveillance" else a for a in agents]
    orch = Orchestrator(store=SqliteEventStore(":memory:"), agents=agents)
    orch.enroll("DEMO-B", cohort_b())
    proj = orch.store.projection("DEMO-B")
    assert proj["staleness"]["tand_briefing"]["status"] == "pending"
    # other sections still ok — failure is isolated
    assert proj["staleness"]["variant_interp"]["status"] == "ok"


def test_dispatch_order_is_topological():
    assert DISPATCH_ORDER[0] == "phenome_mapper"
    assert DISPATCH_ORDER[-1] == "therapeutics_strategist"
    for name, deps in AGENT_DEPENDS_ON.items():
        if EventType.phenome_mapped in deps and name != "phenome_mapper":
            assert DISPATCH_ORDER.index(name) > DISPATCH_ORDER.index("phenome_mapper")
