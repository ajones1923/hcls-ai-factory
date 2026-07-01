"""
Event taxonomy (PRD §2.5.1). The whole engine is event-sourced: every ingest,
every agent output, and every orchestrator decision is an append-only event.
The router keys all decisions off this enum.
"""
from __future__ import annotations

from datetime import datetime, timezone
from enum import Enum

from pydantic import BaseModel, Field


class EventType(str, Enum):
    cohort_loaded = "cohort_loaded"
    patient_enrolled = "patient_enrolled"
    genomic_substrate_ready = "genomic_substrate_ready"
    variant_curated = "variant_curated"
    phenome_mapped = "phenome_mapped"
    trajectory_forecast = "trajectory_forecast"
    tand_surveyed = "tand_surveyed"
    therapeutics_briefed = "therapeutics_briefed"
    surface_requested = "surface_requested"
    surface_assembled = "surface_assembled"
    alert_emitted = "alert_emitted"
    agent_failed = "agent_failed"
    provenance_logged = "provenance_logged"


# Which event an agent emits on success, and the prerequisite events that gate it
# (PRD §2.6 dependency-ordered enrollment: Phenome Mapper runs first).
AGENT_EMITS: dict[str, EventType] = {
    "phenome_mapper": EventType.phenome_mapped,
    "variant_curator": EventType.variant_curated,
    "trajectory_modeler": EventType.trajectory_forecast,
    "tand_surveillance": EventType.tand_surveyed,
    "therapeutics_strategist": EventType.therapeutics_briefed,
}

AGENT_DEPENDS_ON: dict[str, list[EventType]] = {
    "phenome_mapper": [EventType.patient_enrolled],
    "variant_curator": [EventType.genomic_substrate_ready],
    "trajectory_modeler": [EventType.phenome_mapped],
    "tand_surveillance": [EventType.phenome_mapped],
    "therapeutics_strategist": [
        EventType.phenome_mapped,
        EventType.variant_curated,
        EventType.trajectory_forecast,
        EventType.tand_surveyed,
    ],
}


def utcnow() -> datetime:
    return datetime.now(timezone.utc)


class Event(BaseModel):
    """The append-only event envelope (mirrors the engine_events row)."""

    event_id: int | None = None
    patient_id: str
    event_type: EventType
    payload: dict = Field(default_factory=dict)
    provenance: dict = Field(default_factory=dict)
    created_at: datetime = Field(default_factory=utcnow)
    parent_event: int | None = None
