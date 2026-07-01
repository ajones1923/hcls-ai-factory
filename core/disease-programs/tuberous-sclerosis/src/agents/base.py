"""
Uniform agent contract (PRD §2.3).

Every agent — classical statistics for Trajectory, Opus-class synthesis for
Therapeutics — exposes the same surface so the orchestrator can treat them as
interchangeable event consumers: a name, the event it emits, its event
prerequisites, and run(ctx) -> AgentOutput. Agents never call each other; they
consume the projection and emit an output the orchestrator appends to the log.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from pydantic import BaseModel, Field

from src.orchestrator.events import AGENT_DEPENDS_ON, AGENT_EMITS, EventType
from src.utils.provenance import Provenance


@dataclass
class AgentContext:
    """What an agent sees: the patient, the current projection, and (for the demo)
    the patient's synthetic source data. Agents read, never mutate, the projection."""

    patient_id: str
    projection: dict = field(default_factory=dict)
    cohort: dict = field(default_factory=dict)   # synthetic source data for this patient


class AgentOutput(BaseModel):
    agent: str
    status: str = "ok"                            # ok | failed
    data: dict = Field(default_factory=dict)      # the agent's schema-validated payload
    provenance: list[dict] = Field(default_factory=list)
    error: str | None = None


class Agent(ABC):
    name: str = "agent"

    @property
    def emits(self) -> EventType:
        return AGENT_EMITS[self.name]

    @property
    def depends_on(self) -> list[EventType]:
        return AGENT_DEPENDS_ON[self.name]

    @abstractmethod
    def run(self, ctx: AgentContext) -> AgentOutput:
        """Produce this agent's output for the patient."""

    # convenience for runners
    @staticmethod
    def ok(name: str, data: dict, provenance: list[Provenance]) -> AgentOutput:
        return AgentOutput(
            agent=name, status="ok", data=data,
            provenance=[p.to_dict() for p in provenance],
        )
