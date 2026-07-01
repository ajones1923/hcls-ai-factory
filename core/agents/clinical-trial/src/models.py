"""Pydantic data models for the Clinical Trial Intelligence Agent.

Comprehensive enums and models for a clinical trial RAG-based decision
support system covering protocol design, patient matching, site selection,
eligibility optimization, adaptive design, safety signal detection,
regulatory document generation, and competitive intelligence.

Follows the same dataclass/Pydantic pattern as:
  - cardiology_intelligence_agent/src/models.py
  - pharmacogenomics_intelligence_agent/src/models.py

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


# ===================================================================
# ENUMS
# ===================================================================


class TrialWorkflowType(str, Enum):
    """Types of clinical trial query workflows."""
    PROTOCOL_DESIGN = "protocol_design"
    PATIENT_MATCHING = "patient_matching"
    SITE_SELECTION = "site_selection"
    ELIGIBILITY_OPTIMIZATION = "eligibility_optimization"
    ELIGIBILITY_ANALYSIS = "eligibility_analysis"
    ENDPOINT_STRATEGY = "endpoint_strategy"
    ADAPTIVE_DESIGN = "adaptive_design"
    SAFETY_SIGNAL = "safety_signal"
    SAFETY_MONITORING = "safety_monitoring"
    REGULATORY_DOCS = "regulatory_docs"
    REGULATORY_STRATEGY = "regulatory_strategy"
    COMPETITIVE_INTEL = "competitive_intel"
    COMPETITIVE_INTELLIGENCE = "competitive_intelligence"
    BIOMARKER_STRATEGY = "biomarker_strategy"
    RWE_ANALYSIS = "rwe_analysis"
    RECRUITMENT_OPTIMIZATION = "recruitment_optimization"
    DIVERSITY_ASSESSMENT = "diversity_assessment"
    DECENTRALIZED_PLANNING = "decentralized_planning"
    GENERAL = "general"


class TrialPhase(str, Enum):
    """Clinical trial phases."""
    PHASE_I = "phase_i"
    PHASE_I_II = "phase_i_ii"
    PHASE_II = "phase_ii"
    PHASE_II_III = "phase_ii_iii"
    PHASE_III = "phase_iii"
    PHASE_IV = "phase_iv"
    NOT_APPLICABLE = "not_applicable"


class TrialStatus(str, Enum):
    """Clinical trial recruitment/activity status."""
    RECRUITING = "recruiting"
    ACTIVE_NOT_RECRUITING = "active_not_recruiting"
    COMPLETED = "completed"
    TERMINATED = "terminated"
    WITHDRAWN = "withdrawn"
    SUSPENDED = "suspended"
    NOT_YET_RECRUITING = "not_yet_recruiting"


class EvidenceLevel(str, Enum):
    """Level of evidence classification for clinical trial data.

    A1: Systematic review of RCTs
    A2: High-quality RCT
    B: Non-randomized controlled study
    C: Observational study
    D: Case series / case report
    E: Expert opinion / consensus
    """
    A1 = "a1"
    A2 = "a2"
    B = "b"
    C = "c"
    D = "d"
    E = "e"


class CriterionType(str, Enum):
    """Type of eligibility criterion."""
    INCLUSION = "inclusion"
    EXCLUSION = "exclusion"


class EndpointType(str, Enum):
    """Type of clinical trial endpoint."""
    PRIMARY = "primary"
    SECONDARY = "secondary"
    EXPLORATORY = "exploratory"
    SAFETY = "safety"


class RegulatoryAgency(str, Enum):
    """Major regulatory agencies."""
    FDA = "fda"
    EMA = "ema"
    PMDA = "pmda"
    HEALTH_CANADA = "health_canada"
    TGA = "tga"
    MHRA = "mhra"


class DocumentType(str, Enum):
    """Regulatory document types."""
    IND = "ind"
    CSR = "csr"
    BRIEFING = "briefing"
    PSP = "psp"
    RMP = "rmp"
    DSUR = "dsur"


class SeverityLevel(str, Enum):
    """Clinical finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


class TherapeuticArea(str, Enum):
    """Therapeutic area classification."""
    ONCOLOGY = "oncology"
    CARDIOLOGY = "cardiology"
    NEUROLOGY = "neurology"
    IMMUNOLOGY = "immunology"
    RARE_DISEASE = "rare_disease"
    INFECTIOUS_DISEASE = "infectious_disease"
    METABOLIC = "metabolic"
    RESPIRATORY = "respiratory"
    DERMATOLOGY = "dermatology"
    OPHTHALMOLOGY = "ophthalmology"
    HEMATOLOGY = "hematology"
    GASTROENTEROLOGY = "gastroenterology"
    OTHER = "other"


class DCTComponent(str, Enum):
    """Decentralized clinical trial technology components."""
    ECONSENT = "econsent"
    TELEMEDICINE = "telemedicine"
    HOME_HEALTH = "home_health"
    LOCAL_LABS = "local_labs"
    WEARABLES = "wearables"
    EPRO_ECOA = "epro_ecoa"
    DIRECT_TO_PATIENT = "direct_to_patient"


# ===================================================================
# PYDANTIC MODELS - QUERY & SEARCH
# ===================================================================


class TrialQuery(BaseModel):
    """Input query to the Clinical Trial Intelligence Agent."""
    question: str = Field(..., min_length=1, description="Clinical trial question")
    workflow_type: Optional[TrialWorkflowType] = Field(
        default=None,
        description="Specific workflow to route the query to; auto-detected if omitted",
    )
    patient_context: Optional[Dict] = Field(
        default=None,
        description="Patient demographics, biomarkers, genomic data, and clinical history",
    )
    top_k: int = Field(
        default=5, ge=1, le=50,
        description="Number of results to return per collection",
    )
    include_guidelines: bool = Field(
        default=True,
        description="Whether to include guideline references in the response",
    )


class TrialSearchResult(BaseModel):
    """A single search result from any trial knowledge collection."""
    collection: str = Field(..., description="Source Milvus collection name")
    content: str = Field(..., description="Retrieved text content")
    score: float = Field(..., ge=0.0, description="Similarity score")
    metadata: Dict = Field(default_factory=dict, description="Source metadata")


# ===================================================================
# PYDANTIC MODELS - PATIENT MATCHING
# ===================================================================


class MatchScore(BaseModel):
    """Score for a single eligibility criterion match."""
    criterion_text: str = Field(..., description="The eligibility criterion text")
    criterion_type: CriterionType = Field(..., description="Inclusion or exclusion")
    met: bool = Field(..., description="Whether the criterion is met")
    confidence: float = Field(
        ..., ge=0.0, le=1.0,
        description="Confidence in the match assessment",
    )
    evidence: str = Field(
        default="",
        description="Supporting evidence for the match decision",
    )


class OverallMatch(BaseModel):
    """Overall patient-trial match assessment."""
    trial_id: str = Field(..., description="ClinicalTrials.gov NCT identifier")
    trial_title: str = Field(..., description="Trial title")
    phase: TrialPhase = Field(..., description="Trial phase")
    status: TrialStatus = Field(..., description="Trial recruitment status")
    inclusion_met: int = Field(
        ..., ge=0,
        description="Number of inclusion criteria met",
    )
    inclusion_total: int = Field(
        ..., ge=0,
        description="Total number of inclusion criteria",
    )
    exclusion_clear: int = Field(
        ..., ge=0,
        description="Number of exclusion criteria cleared (not triggered)",
    )
    exclusion_total: int = Field(
        ..., ge=0,
        description="Total number of exclusion criteria",
    )
    overall_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Overall match score (0.0 - 1.0)",
    )
    confidence: float = Field(
        ..., ge=0.0, le=1.0,
        description="Confidence in the overall match",
    )
    site_distances: List[Dict] = Field(
        default_factory=list,
        description="Nearby trial sites with distances",
    )


class PatientProfile(BaseModel):
    """Patient profile for trial matching."""
    age: Optional[int] = Field(
        default=None, ge=0, le=120,
        description="Patient age in years",
    )
    sex: Optional[str] = Field(
        default=None,
        description="Patient sex: male or female",
    )
    diagnosis: Optional[str] = Field(
        default=None, max_length=512,
        description="Primary diagnosis",
    )
    biomarkers: List[str] = Field(
        default_factory=list,
        description="List of biomarker results (e.g., HER2+, EGFR T790M)",
    )
    medications: List[str] = Field(
        default_factory=list,
        description="Current medications",
    )
    genomic_variants: List[str] = Field(
        default_factory=list,
        description="Known genomic variants",
    )
    comorbidities: List[str] = Field(
        default_factory=list,
        description="Comorbid conditions",
    )
    geographic_location: Optional[str] = Field(
        default=None, max_length=256,
        description="Patient geographic location for site matching",
    )


# ===================================================================
# PYDANTIC MODELS - ELIGIBILITY & SITE
# ===================================================================


class EligibilityAnalysis(BaseModel):
    """Analysis of a single eligibility criterion."""
    criterion: str = Field(..., description="Eligibility criterion text")
    population_impact: float = Field(
        ..., ge=0.0, le=1.0,
        description="Estimated fraction of target population excluded by this criterion",
    )
    scientific_justification_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Strength of scientific justification for the criterion",
    )
    competitor_comparison: str = Field(
        default="",
        description="How competing trials handle this criterion",
    )
    recommendation: str = Field(
        default="",
        description="Recommendation to broaden, narrow, or retain the criterion",
    )


class SiteScore(BaseModel):
    """Site evaluation score for site selection workflow."""
    site_id: str = Field(..., description="Site identifier")
    facility_name: str = Field(..., description="Facility name")
    city: str = Field(..., description="City")
    country: str = Field(..., description="Country")
    enrollment_rate: float = Field(
        ..., ge=0.0,
        description="Historical enrollment rate (patients/month)",
    )
    screen_failure_rate: float = Field(
        ..., ge=0.0, le=1.0,
        description="Historical screen failure rate (0.0 - 1.0)",
    )
    diversity_index: float = Field(
        ..., ge=0.0, le=1.0,
        description="Site diversity index (0.0 - 1.0)",
    )
    overall_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Overall site suitability score",
    )


# ===================================================================
# PYDANTIC MODELS - SAFETY & COMPETITIVE
# ===================================================================


class SafetySignal(BaseModel):
    """Detected safety signal from trial data."""
    event_type: str = Field(..., description="Adverse event type (MedDRA preferred term)")
    severity: SeverityLevel = Field(..., description="Severity classification")
    frequency: float = Field(
        ..., ge=0.0, le=1.0,
        description="Observed frequency in the study population",
    )
    prr: Optional[float] = Field(
        default=None, ge=0.0,
        description="Proportional reporting ratio",
    )
    ror: Optional[float] = Field(
        default=None, ge=0.0,
        description="Reporting odds ratio",
    )
    causality_assessment: str = Field(
        default="",
        description="Causality assessment (certain, probable, possible, unlikely)",
    )


class CompetitorProfile(BaseModel):
    """Competitive intelligence profile for a rival trial."""
    trial_id: str = Field(..., description="ClinicalTrials.gov NCT identifier")
    sponsor: str = Field(..., description="Sponsoring organization")
    phase: TrialPhase = Field(..., description="Trial phase")
    indication: str = Field(..., description="Target indication")
    mechanism: str = Field(..., description="Mechanism of action")
    enrollment_progress: float = Field(
        ..., ge=0.0, le=1.0,
        description="Enrollment progress (0.0 - 1.0)",
    )
    estimated_completion: str = Field(
        default="",
        description="Estimated primary completion date",
    )
    threat_level: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Competitive threat level",
    )


# ===================================================================
# PYDANTIC MODELS - PROTOCOL COMPLEXITY
# ===================================================================


class ProtocolComplexity(BaseModel):
    """Protocol complexity assessment."""
    procedure_count: int = Field(
        ..., ge=0,
        description="Total number of distinct procedures",
    )
    visit_count: int = Field(
        ..., ge=0,
        description="Total number of scheduled visits",
    )
    endpoint_count: int = Field(
        ..., ge=0,
        description="Total number of endpoints",
    )
    eligibility_criteria_count: int = Field(
        ..., ge=0,
        description="Total number of eligibility criteria",
    )
    complexity_score: float = Field(
        ..., ge=0.0, le=1.0,
        description="Normalized complexity score (0.0 - 1.0)",
    )
    percentile_rank: float = Field(
        ..., ge=0.0, le=100.0,
        description="Percentile rank relative to similar trials",
    )


# ===================================================================
# PYDANTIC MODELS - WORKFLOW & RESPONSE
# ===================================================================


class WorkflowResult(BaseModel):
    """Output from a single clinical trial workflow execution."""
    workflow_type: TrialWorkflowType = Field(
        ...,
        description="Workflow that generated these results",
    )
    findings: List[str] = Field(
        default_factory=list,
        description="Key findings from the workflow",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Actionable recommendations",
    )
    guideline_references: List[str] = Field(
        default_factory=list,
        description="Supporting guideline or regulatory citations",
    )
    severity: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Overall severity of findings",
    )
    cross_agent_triggers: List[str] = Field(
        default_factory=list,
        description="Triggers for cross-agent consultation (e.g., genomics, safety)",
    )
    confidence: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Workflow confidence score (0.0 - 1.0)",
    )


class TrialResponse(BaseModel):
    """Top-level output from the Clinical Trial Intelligence Agent."""
    answer: str = Field(..., description="Synthesized answer to the trial question")
    citations: List[Dict] = Field(
        default_factory=list,
        description="Source citations (collection, id, title, score)",
    )
    workflow_results: List[WorkflowResult] = Field(
        default_factory=list,
        description="Results from each executed workflow",
    )
    matches: List[OverallMatch] = Field(
        default_factory=list,
        description="Patient-trial match results if applicable",
    )
    confidence: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Agent confidence in the response (0.0 - 1.0)",
    )


# ===================================================================
# DATACLASS - SEARCH PLAN
# ===================================================================


@dataclass
class SearchPlan:
    """Pre-retrieval search planning for a clinical trial query.

    Built by the query analyzer to guide collection routing, sub-question
    decomposition, and search-strategy selection before retrieval.
    """
    question: str = ""
    therapeutic_areas: List[str] = field(default_factory=list)
    drugs: List[str] = field(default_factory=list)
    biomarkers: List[str] = field(default_factory=list)
    relevant_workflows: List[TrialWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad | targeted | comparative | regulatory
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)
