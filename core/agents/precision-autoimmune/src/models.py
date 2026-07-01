"""
Data models for the Autoimmune Intelligence Agent.

Defines structures for autoantibody panels, disease activity scores,
HLA associations, biologics pharmacogenomics, and flare prediction.

Part of the HCLS AI Factory: Patient DNA -> Drug Candidates pipeline.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from enum import Enum, unique
from typing import Any, Dict, List, Optional

from pydantic import BaseModel, Field, model_validator


@unique
class AutoimmuneDisease(str, Enum):
    """Autoimmune diseases supported by the agent."""
    RHEUMATOID_ARTHRITIS = "rheumatoid_arthritis"
    SYSTEMIC_LUPUS = "systemic_lupus_erythematosus"
    MULTIPLE_SCLEROSIS = "multiple_sclerosis"
    TYPE_1_DIABETES = "type_1_diabetes"
    INFLAMMATORY_BOWEL = "inflammatory_bowel_disease"
    PSORIASIS = "psoriasis"
    ANKYLOSING_SPONDYLITIS = "ankylosing_spondylitis"
    SJOGRENS_SYNDROME = "sjogrens_syndrome"
    SYSTEMIC_SCLEROSIS = "systemic_sclerosis"
    MYASTHENIA_GRAVIS = "myasthenia_gravis"
    CELIAC_DISEASE = "celiac_disease"
    GRAVES_DISEASE = "graves_disease"
    HASHIMOTO_THYROIDITIS = "hashimoto_thyroiditis"


@unique
class DiseaseActivityLevel(str, Enum):
    """Standardized disease activity levels."""
    REMISSION = "remission"
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    VERY_HIGH = "very_high"


@unique
class FlareRisk(str, Enum):
    """Predicted flare risk level."""
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    IMMINENT = "imminent"


class AutoantibodyResult(BaseModel):
    """Single autoantibody test result."""
    antibody: str = Field(..., description="Autoantibody name (e.g., 'ANA', 'anti-dsDNA', 'RF')")
    value: float = Field(..., description="Measured value")
    unit: str = Field(default="", description="Unit of measurement")
    reference_range: str = Field(default="", description="Normal reference range")
    positive: bool = Field(default=False, description="Whether result is positive")
    titer: Optional[str] = Field(None, description="Titer if applicable (e.g., '1:320')")
    pattern: Optional[str] = Field(None, description="Staining pattern for ANA (e.g., 'homogeneous', 'speckled')")


class AutoantibodyPanel(BaseModel):
    """Complete autoantibody panel results."""
    patient_id: str
    collection_date: str = Field(..., description="ISO-8601 date")
    results: List[AutoantibodyResult] = Field(default_factory=list)

    @property
    def positive_antibodies(self) -> List[str]:
        return [r.antibody for r in self.results if r.positive]

    @property
    def positive_count(self) -> int:
        return len(self.positive_antibodies)


class HLAProfile(BaseModel):
    """HLA typing results for disease association analysis."""
    hla_a: List[str] = Field(default_factory=list, description="HLA-A alleles")
    hla_b: List[str] = Field(default_factory=list, description="HLA-B alleles")
    hla_c: List[str] = Field(default_factory=list, description="HLA-C alleles")
    hla_drb1: List[str] = Field(default_factory=list, description="HLA-DRB1 alleles")
    hla_dqb1: List[str] = Field(default_factory=list, description="HLA-DQB1 alleles")

    @property
    def all_alleles(self) -> List[str]:
        return self.hla_a + self.hla_b + self.hla_c + self.hla_drb1 + self.hla_dqb1


class DiseaseActivityScore(BaseModel):
    """Standardized disease activity score."""
    disease: AutoimmuneDisease
    score_name: str = Field(..., description="Score system (e.g., 'DAS28-CRP', 'SLEDAI-2K', 'CDAI')")
    score_value: float
    level: DiseaseActivityLevel
    components: Dict[str, float] = Field(default_factory=dict, description="Individual score components")
    thresholds: Dict[str, float] = Field(
        default_factory=dict,
        description="Thresholds for each level (e.g., {'remission': 2.6, 'low': 3.2, ...})",
    )


class FlarePredictor(BaseModel):
    """Flare prediction result."""
    disease: AutoimmuneDisease
    current_activity: DiseaseActivityLevel
    predicted_risk: FlareRisk
    risk_score: float = Field(..., ge=0, le=1, description="Numeric risk 0-1")
    contributing_factors: List[str] = Field(default_factory=list)
    protective_factors: List[str] = Field(default_factory=list)
    recommended_monitoring: List[str] = Field(default_factory=list)
    time_horizon_days: int = Field(default=90, description="Prediction window in days")


class BiologicTherapy(BaseModel):
    """Biologic therapy recommendation with PGx context."""
    drug_name: str
    drug_class: str = Field(..., description="e.g., 'TNF inhibitor', 'IL-6 inhibitor', 'B-cell depleter'")
    mechanism: str = ""
    indicated_diseases: List[AutoimmuneDisease] = Field(default_factory=list)
    pgx_considerations: List[str] = Field(default_factory=list)
    contraindications: List[str] = Field(default_factory=list)
    monitoring_requirements: List[str] = Field(default_factory=list)
    efficacy_evidence: str = ""


class AutoimmunePatientProfile(BaseModel):
    """Complete patient profile for autoimmune analysis."""
    patient_id: str = Field(..., description="Unique patient identifier")
    age: int = Field(..., ge=0, le=150)
    sex: str = Field(..., pattern="^[MFmf]$")

    # Autoantibody panel
    autoantibody_panel: Optional[AutoantibodyPanel] = None

    # HLA typing
    hla_profile: Optional[HLAProfile] = None

    # Biomarkers (from Biomarker Agent)
    biomarkers: Dict[str, float] = Field(default_factory=dict)

    # Genotypes (from genomics pipeline)
    genotypes: Dict[str, str] = Field(default_factory=dict)

    # Clinical data
    diagnosed_conditions: List[AutoimmuneDisease] = Field(default_factory=list)
    current_medications: List[str] = Field(default_factory=list)
    symptom_duration_months: Optional[int] = None
    family_history: List[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def _check_has_data(self) -> "AutoimmunePatientProfile":
        if not self.autoantibody_panel and not self.hla_profile and not self.biomarkers:
            raise ValueError("At least one of autoantibody_panel, hla_profile, or biomarkers must be provided")
        return self


class AutoimmuneAnalysisResult(BaseModel):
    """Complete analysis result from the Autoimmune Intelligence Agent."""
    patient_id: str
    disease_activity_scores: List[DiseaseActivityScore] = Field(default_factory=list)
    flare_predictions: List[FlarePredictor] = Field(default_factory=list)
    hla_associations: List[Dict[str, Any]] = Field(default_factory=list)
    biologic_recommendations: List[BiologicTherapy] = Field(default_factory=list)
    critical_alerts: List[str] = Field(default_factory=list)
    cross_agent_findings: List[Dict[str, Any]] = Field(default_factory=list)


# =====================================================================
# Collection record models (for embedding text generation)
# =====================================================================

class ClinicalDocumentRecord(BaseModel):
    """Record for autoimmune_clinical_documents collection."""
    id: str
    text_chunk: str
    patient_id: str = ""
    doc_type: str = ""
    specialty: str = ""
    provider: str = ""
    visit_date: str = ""
    source_file: str = ""
    page_number: int = 0
    chunk_index: int = 0

    def to_embedding_text(self) -> str:
        parts = [self.text_chunk]
        if self.doc_type:
            parts.append(f"Document type: {self.doc_type}")
        if self.specialty:
            parts.append(f"Specialty: {self.specialty}")
        return " ".join(parts)


class LabResultRecord(BaseModel):
    """Record for autoimmune_patient_labs collection."""
    id: str
    text_chunk: str
    patient_id: str = ""
    test_name: str = ""
    value: float = 0.0
    unit: str = ""
    reference_range: str = ""
    flag: str = "normal"
    collection_date: str = ""
    panel_name: str = ""

    def to_embedding_text(self) -> str:
        return (
            f"{self.test_name}: {self.value} {self.unit} ({self.flag}). "
            f"Reference range: {self.reference_range}. {self.text_chunk}"
        )


class TimelineEventRecord(BaseModel):
    """Record for autoimmune_patient_timelines collection."""
    id: str
    text_chunk: str
    patient_id: str = ""
    event_type: str = ""
    event_date: str = ""
    description: str = ""
    provider: str = ""
    specialty: str = ""
    days_from_first_symptom: int = 0

    def to_embedding_text(self) -> str:
        return (
            f"[{self.event_date}] {self.event_type}: {self.description}. "
            f"Specialty: {self.specialty}. {self.text_chunk}"
        )
