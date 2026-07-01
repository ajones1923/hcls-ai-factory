"""
LangGraph orchestrator runtime (PRD §2.2, §2.6). The PRD names LangGraph as the
production runtime for the deterministic router. This wraps the same agents and the
same event store as the plain Orchestrator in a compiled LangGraph StateGraph, with
explicit nodes (one per agent) and deterministic edges in dependency order — so the
workflow topology is itself inspectable demonstration material. It produces byte-identical
state to the plain dispatcher (verified in tests); LangGraph adds the explicit graph, not
new behavior.

Requires `langgraph` (in the engine venv). The default Orchestrator does not import it.
"""
from __future__ import annotations

from typing import TypedDict

from langgraph.graph import END, START, StateGraph

from src.orchestrator.events import Event, EventType
from src.orchestrator.graph import Orchestrator


class _GState(TypedDict):
    patient_id: str
    cohort: dict


class LangGraphRunner:
    def __init__(self, orchestrator: Orchestrator | None = None) -> None:
        self.orch = orchestrator or Orchestrator()
        self.app = self._build()

    def _build(self):
        order = [a.name for a in self.orch.agents]
        by_name = {a.name: a for a in self.orch.agents}
        sg = StateGraph(_GState)

        def make_node(agent):
            def node(state: _GState) -> dict:
                self.orch.run_agent(state["patient_id"], agent, state["cohort"])
                return {}
            return node

        for name in order:
            sg.add_node(name, make_node(by_name[name]))
        sg.add_edge(START, order[0])
        for a, b in zip(order, order[1:]):
            sg.add_edge(a, b)
        sg.add_edge(order[-1], END)
        return sg.compile()

    def enroll(self, patient_id: str, cohort: dict | None = None) -> None:
        cohort = cohort or {}
        self.orch.store.append(Event(patient_id=patient_id, event_type=EventType.patient_enrolled,
                                     payload={"cohort_keys": sorted(cohort.keys())}))
        self.orch.store.append(Event(patient_id=patient_id, event_type=EventType.genomic_substrate_ready,
                                     payload={"source": "synthetic", "ref": f"tsc-cohort/{patient_id}"}))
        self.app.invoke({"patient_id": patient_id, "cohort": cohort})

    def assemble_surface(self, patient_id: str, kind: str) -> dict:
        return self.orch.assemble_surface(patient_id, kind)

    def mermaid(self) -> str:
        """The graph topology as Mermaid (demonstration material)."""
        try:
            return self.app.get_graph().draw_mermaid()
        except Exception:
            return "graph TD; " + " --> ".join(a.name for a in self.orch.agents)
