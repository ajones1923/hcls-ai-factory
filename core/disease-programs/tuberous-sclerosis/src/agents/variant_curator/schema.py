"""Output schema for the TSC-Variant Curator (PRD §3 FR-VC-8; ClinVar-spec shape)."""
from __future__ import annotations

from pydantic import BaseModel, Field


class CriterionOut(BaseModel):
    code: str
    bucket: str
    rationale: str = ""


class VariantInterpretation(BaseModel):
    gene: str | None = None
    hgvsc: str | None = None
    hgvsp: str | None = None
    consequence: str | None = None
    vaf: float | None = None
    mosaic: bool = False
    recovered: bool = False                 # found in tissue where standard blood testing was negative
    acmg_classification: str = "Variant of Uncertain Significance"
    acmg_rule: str = ""
    acmg_criteria: list[CriterionOut] = Field(default_factory=list)
    confidence: str = "low"
    ddpcr_recommended: bool = False
    # read-level evidence (FR-VC-4)
    depth: int | None = None
    alt_reads: int | None = None
    strand_balance: float | None = None
    artifact_assessment: str = "not assessed"   # pass | strand-bias-suspected | not assessed


class VariantCuratorOutput(BaseModel):
    primary: VariantInterpretation | None = None
    all_variants: list[VariantInterpretation] = Field(default_factory=list)
    review_status: str = "pending_molecular_geneticist"   # human gate — non-negotiable
    synthesis_narrative: dict | None = None               # Opus per-criterion reasoning (when keyed)
    note: str | None = None
    draft_report_ref: str | None = None
    source: str = "vcf"                                    # vcf | inline | none
    artifacts_filtered: int = 0                            # strand-biased candidates rejected (FR-VC-1)
    build_status: str = "real ACMG-AMP combinatorial classifier; live ClinVar/gnomAD/LOVD = W2"
