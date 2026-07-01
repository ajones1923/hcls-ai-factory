"""
Orchestration policies (PRD §2.6): the agent registry and the dependency-ordered
dispatch plan. Phenome Mapper runs first; the others gate on its output; the
Therapeutics Strategist fans in from all upstream agents.
"""
from __future__ import annotations

from src.agents.base import Agent
from src.agents.phenome_mapper.runner import PhenomeMapper
from src.agents.variant_curator.runner import VariantCurator
from src.agents.trajectory_modeler.runner import TrajectoryModeler
from src.agents.tand_surveillance.runner import TANDSurveillance
from src.agents.therapeutics_strategist.runner import TherapeuticsStrategist

# Topological dispatch order satisfying AGENT_DEPENDS_ON (src/orchestrator/events.py).
DISPATCH_ORDER = [
    "phenome_mapper",
    "variant_curator",
    "trajectory_modeler",
    "tand_surveillance",
    "therapeutics_strategist",
]


def default_agents() -> list[Agent]:
    agents = {
        "phenome_mapper": PhenomeMapper(),
        "variant_curator": VariantCurator(),
        "trajectory_modeler": TrajectoryModeler(),
        "tand_surveillance": TANDSurveillance(),
        "therapeutics_strategist": TherapeuticsStrategist(),
    }
    return [agents[name] for name in DISPATCH_ORDER]
