"""Output schema for the TSC-Phenome Mapper (PRD §3 FR-PM-*)."""
from __future__ import annotations

from pydantic import BaseModel, Field


class HPOAssertion(BaseModel):
    hpo_id: str
    label: str
    onset: str | None = None
    evidence_span: str | None = None       # source text supporting the assertion
    source: str = "structured"             # structured | note | imaging
    confidence: str = "high"
    validation: str = "unverified"         # verified|relabeled|remapped|recovered|unknown|unverified
    polarity: str = "present"              # present | absent (negated mentions are not admitted)
    temporality: str = "current"           # onset | current | historical


class SurveillanceItem(BaseModel):
    type: str
    interval_months: int
    status: str                            # due_per_schedule | overdue
    rationale: str = ""


class PhenomeProfile(BaseModel):
    hpo_terms: list[HPOAssertion] = Field(default_factory=list)
    discordances: list[dict] = Field(default_factory=list)
    surveillance_gaps: list[SurveillanceItem] = Field(default_factory=list)
    n_terms: int = 0
    ontology_release: str = "unavailable"  # HPO release status; set by the runtime validator
    n_dropped_unverified: int = 0          # model-emitted terms rejected as not in the ontology
    build_status: str = "structured HPO + ITSC gap analyzer real; per-note NLP live when keyed"
