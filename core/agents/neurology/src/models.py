"""Pydantic data models for the Neurology Intelligence Agent.

Comprehensive enums and models for a neurology RAG-based clinical
decision support system covering acute stroke triage, dementia evaluation,
epilepsy focus localization, brain tumor grading, MS monitoring,
Parkinson's assessment, headache classification, neuromuscular
evaluation, and general neurology workflows.

Follows the same dataclass/Pydantic pattern as:
  - clinical_trial_intelligence_agent/src/models.py
  - rare_disease_diagnostic_agent/src/models.py

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


class NeuroWorkflowType(str, Enum):
    """Types of neurology query workflows."""
    ACUTE_STROKE = "acute_stroke"
    DEMENTIA_EVALUATION = "dementia_evaluation"
    EPILEPSY_FOCUS = "epilepsy_focus"
    BRAIN_TUMOR = "brain_tumor"
    MS_MONITORING = "ms_monitoring"
    PARKINSONS_ASSESSMENT = "parkinsons_assessment"
    HEADACHE_CLASSIFICATION = "headache_classification"
    NEUROMUSCULAR_EVALUATION = "neuromuscular_evaluation"
    GENERAL = "general"


class NeuroDomain(str, Enum):
    """Primary neurology sub-domains for query routing."""
    CEREBROVASCULAR = "cerebrovascular"
    NEURODEGENERATIVE = "neurodegenerative"
    EPILEPSY = "epilepsy"
    MS = "multiple_sclerosis"
    NEURO_ONCOLOGY = "neuro_oncology"
    MOVEMENT_DISORDERS = "movement_disorders"
    HEADACHE = "headache"
    NEUROMUSCULAR = "neuromuscular"
    NEUROGENETICS = "neurogenetics"
    NEUROIMMUNOLOGY = "neuroimmunology"


class StrokeType(str, Enum):
    """Stroke classification."""
    ISCHEMIC = "ischemic"
    HEMORRHAGIC = "hemorrhagic"
    TIA = "tia"
    SAH = "subarachnoid_hemorrhage"


class DementiaSubtype(str, Enum):
    """Dementia differential subtypes."""
    ALZHEIMERS = "alzheimers"
    FRONTOTEMPORAL = "frontotemporal"
    LEWY_BODY = "lewy_body"
    VASCULAR = "vascular"
    PSP = "progressive_supranuclear_palsy"
    MSA = "multiple_system_atrophy"
    MIXED = "mixed"
    NORMAL_PRESSURE_HYDROCEPHALUS = "nph"
    CJD = "creutzfeldt_jakob"


class ATNStage(str, Enum):
    """ATN biomarker staging for Alzheimer's disease."""
    A_NEG_T_NEG_N_NEG = "A-T-N-"
    A_POS_T_NEG_N_NEG = "A+T-N-"
    A_POS_T_POS_N_NEG = "A+T+N-"
    A_POS_T_POS_N_POS = "A+T+N+"
    A_POS_T_NEG_N_POS = "A+T-N+"
    A_NEG_T_POS_N_NEG = "A-T+N-"
    A_NEG_T_NEG_N_POS = "A-T-N+"
    A_NEG_T_POS_N_POS = "A-T+N+"


class SeizureType(str, Enum):
    """ILAE 2017 seizure classification."""
    FOCAL_AWARE = "focal_aware"
    FOCAL_IMPAIRED = "focal_impaired_awareness"
    FOCAL_TO_BILATERAL = "focal_to_bilateral_tonic_clonic"
    GENERALIZED_TONIC_CLONIC = "generalized_tonic_clonic"
    ABSENCE = "generalized_absence"
    MYOCLONIC = "generalized_myoclonic"
    GENERALIZED_ATONIC = "generalized_atonic"
    GENERALIZED_TONIC = "generalized_tonic"
    UNKNOWN_ONSET = "unknown_onset"


class EpilepsySyndrome(str, Enum):
    """Common epilepsy syndromes."""
    JUVENILE_MYOCLONIC = "juvenile_myoclonic"
    CHILDHOOD_ABSENCE = "childhood_absence"
    TEMPORAL_LOBE = "temporal_lobe"
    FRONTAL_LOBE = "frontal_lobe"
    DRAVET = "dravet"
    LENNOX_GASTAUT = "lennox_gastaut"
    WEST = "west"
    BENIGN_ROLANDIC = "benign_rolandic"
    PROGRESSIVE_MYOCLONIC = "progressive_myoclonic"
    UNCLASSIFIED = "unclassified"


class MSPhenotype(str, Enum):
    """Multiple sclerosis clinical phenotypes (2013 revised criteria)."""
    CIS = "clinically_isolated_syndrome"
    RRMS = "relapsing_remitting"
    SPMS = "secondary_progressive"
    PPMS = "primary_progressive"


class DMTCategory(str, Enum):
    """Disease-modifying therapy categories for MS."""
    PLATFORM = "platform"
    MODERATE_EFFICACY = "moderate_efficacy"
    HIGH_EFFICACY = "high_efficacy"


class TumorGrade(str, Enum):
    """WHO 2021 CNS tumor grading."""
    GRADE_1 = "grade_1"
    GRADE_2 = "grade_2"
    GRADE_3 = "grade_3"
    GRADE_4 = "grade_4"


class TumorMolecularMarker(str, Enum):
    """Key molecular markers for CNS tumors."""
    IDH_MUTANT = "idh_mutant"
    IDH_WILDTYPE = "idh_wildtype"
    MGMT_METHYLATED = "mgmt_methylated"
    MGMT_UNMETHYLATED = "mgmt_unmethylated"
    CODEL_1P19Q = "1p19q_codeleted"
    NO_CODEL_1P19Q = "1p19q_intact"
    H3K27M_MUTANT = "h3k27m_mutant"
    TERT_MUTANT = "tert_mutant"
    ATRX_LOSS = "atrx_loss"
    BRAF_V600E = "braf_v600e"
    EGFR_AMPLIFIED = "egfr_amplified"


class ParkinsonsSubtype(str, Enum):
    """Parkinson's disease motor subtypes."""
    TREMOR_DOMINANT = "tremor_dominant"
    PIGD = "postural_instability_gait_difficulty"
    INDETERMINATE = "indeterminate"


class HeadacheType(str, Enum):
    """ICHD-3 primary headache classification."""
    MIGRAINE_WITHOUT_AURA = "migraine_without_aura"
    MIGRAINE_WITH_AURA = "migraine_with_aura"
    CHRONIC_MIGRAINE = "chronic_migraine"
    TENSION_TYPE_EPISODIC = "tension_type_episodic"
    TENSION_TYPE_CHRONIC = "tension_type_chronic"
    CLUSTER = "cluster"
    TRIGEMINAL_AUTONOMIC = "trigeminal_autonomic"
    MEDICATION_OVERUSE = "medication_overuse"
    NEW_DAILY_PERSISTENT = "new_daily_persistent"
    SECONDARY = "secondary"


class NMJPattern(str, Enum):
    """Neuromuscular junction / EMG pattern classification."""
    AXONAL_MOTOR = "axonal_motor"
    AXONAL_SENSORY = "axonal_sensory"
    AXONAL_SENSORIMOTOR = "axonal_sensorimotor"
    DEMYELINATING = "demyelinating"
    NMJ_PRESYNAPTIC = "nmj_presynaptic"
    NMJ_POSTSYNAPTIC = "nmj_postsynaptic"
    MYOPATHIC = "myopathic"
    MIXED = "mixed"


class SeverityLevel(str, Enum):
    """Clinical finding severity classification."""
    CRITICAL = "critical"
    HIGH = "high"
    MODERATE = "moderate"
    LOW = "low"
    INFORMATIONAL = "informational"


class EvidenceLevel(str, Enum):
    """AAN evidence classification scheme for clinical recommendations."""
    CLASS_I = "class_I"
    CLASS_II = "class_II"
    CLASS_III = "class_III"
    CLASS_IV = "class_IV"


class GuidelineClass(str, Enum):
    """AAN/AHA/ILAE guideline recommendation classification."""
    CLASS_I = "class_i"
    CLASS_IIA = "class_iia"
    CLASS_IIB = "class_iib"
    CLASS_III_NO_BENEFIT = "class_iii_no_benefit"
    CLASS_III_HARM = "class_iii_harm"


class ClinicalScaleType(str, Enum):
    """Validated neurological scale calculators."""
    NIHSS = "nihss"
    GCS = "gcs"
    MOCA = "moca"
    UPDRS = "updrs_part_iii"
    EDSS = "edss"
    MRS = "mrs"
    HIT6 = "hit6"
    ALSFRS = "alsfrs_r"
    ASPECTS = "aspects"
    HOEHN_YAHR = "hoehn_yahr"


# ===================================================================
# PYDANTIC MODELS
# ===================================================================


class NeuroQuery(BaseModel):
    """Incoming neurology query with clinical context."""

    query: str = Field(
        ...,
        description="Free-text clinical query or question",
        min_length=1,
        max_length=4096,
    )
    workflow: NeuroWorkflowType = Field(
        default=NeuroWorkflowType.GENERAL,
        description="Workflow type for collection weight routing",
    )
    domain: Optional[NeuroDomain] = Field(
        default=None,
        description="Primary neurology domain for filtering",
    )
    patient_age: Optional[int] = Field(
        default=None,
        ge=0,
        le=150,
        description="Patient age in years",
    )
    patient_sex: Optional[str] = Field(
        default=None,
        description="Patient biological sex (M/F/Other)",
    )
    top_k: int = Field(
        default=20,
        ge=1,
        le=200,
        description="Maximum number of results to return",
    )
    include_genomic: bool = Field(
        default=False,
        description="Whether to include genomic evidence collection in search",
    )
    conversation_id: Optional[str] = Field(
        default=None,
        description="Conversation session ID for context tracking",
    )


class NeuroSearchResult(BaseModel):
    """A single search result from a neurology collection."""

    collection: str = Field(
        ..., description="Source Milvus collection name"
    )
    score: float = Field(
        ..., ge=0.0, le=1.0, description="Cosine similarity score"
    )
    content: str = Field(
        ..., description="Retrieved text content"
    )
    metadata: Dict = Field(
        default_factory=dict,
        description="Additional metadata fields from the collection record",
    )
    evidence_level: Optional[EvidenceLevel] = Field(
        default=None,
        description="Evidence classification if available",
    )
    guideline_class: Optional[GuidelineClass] = Field(
        default=None,
        description="Guideline recommendation class if applicable",
    )
    citation: Optional[str] = Field(
        default=None,
        description="Formatted citation string",
    )


class StrokeAssessment(BaseModel):
    """Structured acute stroke evaluation data."""

    stroke_type: StrokeType = Field(
        ..., description="Classification of stroke subtype"
    )
    nihss_score: Optional[int] = Field(
        default=None,
        ge=0,
        le=42,
        description="NIH Stroke Scale score (0-42)",
    )
    onset_time_hours: Optional[float] = Field(
        default=None,
        ge=0,
        description="Hours since symptom onset (0 = unknown/wake-up stroke)",
    )
    large_vessel_occlusion: Optional[bool] = Field(
        default=None,
        description="Whether large vessel occlusion is suspected or confirmed",
    )
    tpa_eligible: Optional[bool] = Field(
        default=None,
        description="Whether patient is eligible for IV tPA",
    )
    thrombectomy_eligible: Optional[bool] = Field(
        default=None,
        description="Whether patient is eligible for mechanical thrombectomy",
    )
    aspects_score: Optional[int] = Field(
        default=None,
        ge=0,
        le=10,
        description="Alberta Stroke Program Early CT Score (0-10)",
    )
    toast_classification: Optional[str] = Field(
        default=None,
        description="TOAST classification of ischemic stroke etiology",
    )
    severity: SeverityLevel = Field(
        default=SeverityLevel.HIGH,
        description="Clinical severity level",
    )


class DementiaAssessment(BaseModel):
    """Structured dementia / cognitive decline evaluation."""

    suspected_subtype: Optional[DementiaSubtype] = Field(
        default=None,
        description="Suspected dementia subtype",
    )
    atn_stage: Optional[ATNStage] = Field(
        default=None,
        description="ATN biomarker framework staging",
    )
    mmse_score: Optional[int] = Field(
        default=None,
        ge=0,
        le=30,
        description="Mini-Mental State Examination score (0-30)",
    )
    moca_score: Optional[int] = Field(
        default=None,
        ge=0,
        le=30,
        description="Montreal Cognitive Assessment score (0-30)",
    )
    cdr_global: Optional[float] = Field(
        default=None,
        ge=0,
        le=3,
        description="Clinical Dementia Rating global score (0-3)",
    )
    apoe_genotype: Optional[str] = Field(
        default=None,
        description="APOE genotype (e.g., e3/e4, e4/e4)",
    )
    functional_decline: bool = Field(
        default=False,
        description="Whether functional decline in ADLs is documented",
    )
    duration_months: Optional[int] = Field(
        default=None,
        ge=0,
        description="Duration of cognitive symptoms in months",
    )
    neuroimaging_findings: Optional[str] = Field(
        default=None,
        description="Key neuroimaging findings (e.g., hippocampal atrophy, WMH)",
    )


class SeizureClassification(BaseModel):
    """Structured seizure / epilepsy classification and focus assessment."""

    seizure_type: SeizureType = Field(
        ..., description="ILAE seizure classification"
    )
    epilepsy_syndrome: Optional[EpilepsySyndrome] = Field(
        default=None,
        description="Epilepsy syndrome if identified",
    )
    eeg_findings: Optional[str] = Field(
        default=None,
        description="Key EEG findings (e.g., temporal sharp waves, 3Hz spike-wave)",
    )
    mri_findings: Optional[str] = Field(
        default=None,
        description="Relevant MRI findings (e.g., mesial temporal sclerosis, FCD)",
    )
    seizure_frequency: Optional[str] = Field(
        default=None,
        description="Seizure frequency description (e.g., 2-3/month)",
    )
    drug_resistant: bool = Field(
        default=False,
        description="Whether epilepsy meets ILAE criteria for drug resistance",
    )
    surgical_candidate: Optional[bool] = Field(
        default=None,
        description="Whether patient is a candidate for epilepsy surgery",
    )
    current_aeds: List[str] = Field(
        default_factory=list,
        description="List of current anti-epileptic drugs",
    )


class MSAssessment(BaseModel):
    """Structured multiple sclerosis evaluation and monitoring."""

    phenotype: MSPhenotype = Field(
        ..., description="MS clinical phenotype classification"
    )
    edss_score: Optional[float] = Field(
        default=None,
        ge=0,
        le=10,
        description="Expanded Disability Status Scale score (0-10, 0.5 steps)",
    )
    relapse_count_2yr: Optional[int] = Field(
        default=None,
        ge=0,
        description="Number of relapses in the past 2 years",
    )
    new_t2_lesions: Optional[int] = Field(
        default=None,
        ge=0,
        description="Number of new or enlarging T2 lesions on most recent MRI",
    )
    gadolinium_enhancing: Optional[int] = Field(
        default=None,
        ge=0,
        description="Number of gadolinium-enhancing lesions on most recent MRI",
    )
    current_dmt: Optional[str] = Field(
        default=None,
        description="Current disease-modifying therapy",
    )
    dmt_category: Optional[DMTCategory] = Field(
        default=None,
        description="Category of current DMT",
    )
    jcv_antibody_status: Optional[str] = Field(
        default=None,
        description="JC virus antibody status (positive/negative/unknown)",
    )
    nfl_level: Optional[float] = Field(
        default=None,
        ge=0,
        description="Serum neurofilament light chain level (pg/mL)",
    )
    disease_duration_years: Optional[float] = Field(
        default=None,
        ge=0,
        description="Duration since MS diagnosis in years",
    )


class TumorAssessment(BaseModel):
    """Structured CNS tumor evaluation."""

    tumor_type: Optional[str] = Field(
        default=None,
        description="Histological tumor type (e.g., glioblastoma, meningioma, astrocytoma)",
    )
    who_grade: Optional[TumorGrade] = Field(
        default=None,
        description="WHO CNS tumor grade (2021 classification)",
    )
    location: Optional[str] = Field(
        default=None,
        description="Anatomical tumor location (e.g., left temporal lobe, posterior fossa)",
    )
    molecular_markers: List[TumorMolecularMarker] = Field(
        default_factory=list,
        description="Identified molecular markers",
    )
    size_cm: Optional[float] = Field(
        default=None,
        ge=0,
        description="Maximum tumor dimension in centimeters",
    )
    kps_score: Optional[int] = Field(
        default=None,
        ge=0,
        le=100,
        description="Karnofsky Performance Status (0-100, step 10)",
    )
    prior_treatments: List[str] = Field(
        default_factory=list,
        description="List of prior treatments (surgery, radiation, chemo agents)",
    )
    recurrence: bool = Field(
        default=False,
        description="Whether this is a recurrent tumor",
    )


class HeadacheClassification(BaseModel):
    """Structured headache classification per ICHD-3 criteria."""

    headache_type: Optional[HeadacheType] = Field(
        default=None,
        description="ICHD-3 classification",
    )
    frequency_per_month: Optional[int] = Field(
        default=None,
        ge=0,
        description="Average headache days per month",
    )
    duration_hours: Optional[float] = Field(
        default=None,
        ge=0,
        description="Typical attack duration in hours",
    )
    aura_present: bool = Field(
        default=False,
        description="Whether aura symptoms are present",
    )
    medication_overuse: bool = Field(
        default=False,
        description="Whether medication overuse headache criteria are met",
    )
    red_flags: List[str] = Field(
        default_factory=list,
        description="Red flag symptoms present (e.g., thunderclap, papilledema, fever)",
    )
    preventive_medications: List[str] = Field(
        default_factory=list,
        description="Current or prior preventive medications tried",
    )
    disability_score: Optional[int] = Field(
        default=None,
        ge=0,
        description="MIDAS or HIT-6 disability score",
    )


class NeuromuscularAssessment(BaseModel):
    """Structured neuromuscular disease evaluation."""

    suspected_diagnosis: Optional[str] = Field(
        default=None,
        description="Suspected neuromuscular diagnosis (e.g., ALS, MG, GBS, CIDP)",
    )
    emg_pattern: Optional[NMJPattern] = Field(
        default=None,
        description="EMG/NCS pattern classification",
    )
    weakness_pattern: Optional[str] = Field(
        default=None,
        description="Pattern of weakness (e.g., proximal, distal, bulbar, respiratory)",
    )
    emg_findings: Optional[str] = Field(
        default=None,
        description="Key EMG/NCS findings",
    )
    antibody_panel: Dict[str, str] = Field(
        default_factory=dict,
        description="Relevant antibody results (e.g., AChR: positive, anti-MuSK: negative)",
    )
    ck_level: Optional[float] = Field(
        default=None,
        ge=0,
        description="Creatine kinase level (U/L)",
    )
    respiratory_function: Optional[str] = Field(
        default=None,
        description="Forced vital capacity or other respiratory function metrics",
    )
    progression_rate: Optional[str] = Field(
        default=None,
        description="Rate of disease progression (e.g., rapid, slow, stable)",
    )
    genetic_testing: Optional[str] = Field(
        default=None,
        description="Genetic testing results if applicable",
    )


class ScaleResult(BaseModel):
    """Result of a validated neurological assessment scale."""

    scale_type: ClinicalScaleType = Field(
        ..., description="Calculator used"
    )
    score: float = Field(
        ..., description="Computed score"
    )
    max_score: float = Field(
        ..., description="Maximum possible score"
    )
    interpretation: str = Field(
        ..., description="Clinical interpretation of the score"
    )
    severity_category: str = Field(
        ..., description="Severity or impairment category"
    )
    thresholds: Dict[str, float] = Field(
        default_factory=dict,
        description="Clinical decision thresholds",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Clinical recommendations based on score",
    )


class WorkflowResult(BaseModel):
    """Complete result from a neurology workflow execution."""

    workflow_type: NeuroWorkflowType = Field(
        ..., description="Workflow that generated these results"
    )
    domain: Optional[NeuroDomain] = Field(
        default=None,
        description="Primary neurology domain",
    )
    search_results: List[NeuroSearchResult] = Field(
        default_factory=list,
        description="Retrieved evidence from vector collections",
    )
    findings: List[str] = Field(
        default_factory=list,
        description="Key clinical findings",
    )
    scale_results: List[ScaleResult] = Field(
        default_factory=list,
        description="Computed clinical scale scores if applicable",
    )
    assessment: Optional[str] = Field(
        default=None,
        description="LLM-generated clinical assessment narrative",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Prioritized clinical recommendations",
    )
    guideline_references: List[str] = Field(
        default_factory=list,
        description="Supporting guideline citations",
    )
    evidence_summary: Optional[str] = Field(
        default=None,
        description="Summary of evidence quality and key references",
    )
    severity: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Overall severity of findings",
    )
    cross_modal_triggers: List[str] = Field(
        default_factory=list,
        description="Cross-modal flags for multi-system correlation",
    )
    confidence: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Model confidence in the assessment (0.0-1.0)",
    )


class NeuroResponse(BaseModel):
    """Top-level API response wrapper for neurology queries."""

    query: str = Field(
        ..., description="Original query text"
    )
    workflow: NeuroWorkflowType = Field(
        ..., description="Workflow used for processing"
    )
    result: WorkflowResult = Field(
        ..., description="Complete workflow result"
    )
    citations: List[str] = Field(
        default_factory=list,
        description="Formatted reference citations",
    )
    warnings: List[str] = Field(
        default_factory=list,
        description="Clinical warnings or caveats",
    )
    processing_time_ms: float = Field(
        default=0.0,
        ge=0.0,
        description="Total processing time in milliseconds",
    )


# ===================================================================
# DATACLASSES
# ===================================================================


@dataclass
class SearchPlan:
    """Plan for multi-collection search with per-collection parameters.

    Created by the workflow router to specify which collections to search,
    with what weights, and how many results to retrieve from each.
    """

    collections: List[str] = field(default_factory=list)
    weights: Dict[str, float] = field(default_factory=dict)
    top_k_per_collection: Dict[str, int] = field(default_factory=dict)
    filters: Dict[str, str] = field(default_factory=dict)
    workflow: NeuroWorkflowType = NeuroWorkflowType.GENERAL
    domain: Optional[NeuroDomain] = None
