"""Pydantic data models for the Cardiology Intelligence Agent.

Comprehensive enums and models for a cardiovascular RAG-based clinical
decision support system covering CAD assessment, heart failure management,
valvular disease, arrhythmia evaluation, cardiac imaging, stress testing,
preventive risk scoring, and cardio-oncology surveillance.

Follows the same dataclass/Pydantic pattern as:
  - pharmacogenomics_intelligence_agent/src/models.py
  - cart_intelligence_agent/src/models.py (CARTLiterature, ClinicalTrial)
  - rag-chat-pipeline/src/vcf_parser.py (VariantEvidence)

Author: Adam Jones
Date: March 2026
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional

from pydantic import BaseModel, Field, field_validator


# ===================================================================
# ENUMS
# ===================================================================


class CardioWorkflowType(str, Enum):
    """Types of cardiology query workflows."""
    CAD_ASSESSMENT = "cad_assessment"
    HEART_FAILURE = "heart_failure"
    VALVULAR_DISEASE = "valvular_disease"
    ARRHYTHMIA = "arrhythmia"
    CARDIAC_MRI = "cardiac_mri"
    STRESS_TEST = "stress_test"
    PREVENTIVE_RISK = "preventive_risk"
    CARDIO_ONCOLOGY = "cardio_oncology"
    ACUTE_DECOMPENSATED_HF = "acute_decompensated_hf"
    POST_MI = "post_mi"
    MYOCARDITIS_PERICARDITIS = "myocarditis_pericarditis"
    GENERAL = "general"


class RiskScoreType(str, Enum):
    """Validated cardiovascular risk calculators."""
    ASCVD = "ascvd"
    HEART = "heart"
    CHA2DS2_VASC = "cha2ds2_vasc"
    HAS_BLED = "has_bled"
    MAGGIC = "maggic"
    EUROSCORE_II = "euroscore_ii"


class SeverityLevel(str, Enum):
    """Clinical finding severity classification."""
    CRITICAL = "critical"
    VERY_HIGH = "very_high"
    HIGH = "high"
    INTERMEDIATE = "intermediate"
    MODERATE = "moderate"
    BORDERLINE = "borderline"
    LOW = "low"
    INFORMATIONAL = "informational"


class ImagingModality(str, Enum):
    """Cardiac imaging modalities."""
    ECHOCARDIOGRAPHY = "echocardiography"
    CARDIAC_CT = "cardiac_ct"
    CARDIAC_MRI = "cardiac_mri"
    NUCLEAR_CARDIOLOGY = "nuclear_cardiology"
    CORONARY_ANGIOGRAPHY = "coronary_angiography"


class HeartFailureClass(str, Enum):
    """NYHA functional classification of heart failure symptoms."""
    NYHA_I = "nyha_i"
    NYHA_II = "nyha_ii"
    NYHA_III = "nyha_iii"
    NYHA_IV = "nyha_iv"


class HeartFailureStage(str, Enum):
    """ACC/AHA stages of heart failure progression."""
    STAGE_A = "stage_a"
    STAGE_B = "stage_b"
    STAGE_C = "stage_c"
    STAGE_D = "stage_d"


class EjectionFractionCategory(str, Enum):
    """Heart failure classification by left ventricular ejection fraction.

    HFrEF: LVEF <= 40%
    HFmrEF: LVEF 41-49%
    HFpEF: LVEF >= 50%
    HFimpEF: Previously reduced EF that has improved to > 40%
    """
    HFrEF = "hfref"
    HFmrEF = "hfmref"
    HFpEF = "hfpef"
    HFimpEF = "hfimpef"


class ValveSeverity(str, Enum):
    """Valvular heart disease severity grading."""
    MILD = "mild"
    MODERATE = "moderate"
    SEVERE = "severe"
    CRITICAL = "critical"


class CADRADSScore(str, Enum):
    """CAD-RADS coronary artery disease reporting and data system scores.

    0: No plaque or stenosis (0%)
    1: Minimal stenosis (1-24%)
    2: Mild stenosis (25-49%)
    3: Moderate stenosis (50-69%)
    4A: Severe stenosis (70-99%) in one or two vessels
    4B: Severe stenosis (70-99%) in three vessels or left main >= 50%
    5: Total occlusion (100%)
    """
    CAD_RADS_0 = "cad_rads_0"
    CAD_RADS_1 = "cad_rads_1"
    CAD_RADS_2 = "cad_rads_2"
    CAD_RADS_3 = "cad_rads_3"
    CAD_RADS_4A = "cad_rads_4a"
    CAD_RADS_4B = "cad_rads_4b"
    CAD_RADS_5 = "cad_rads_5"


class AnticoagulationRecommendation(str, Enum):
    """Anticoagulation therapy recommendation status."""
    NO_ANTICOAG = "no_anticoag"
    CONSIDER_ANTICOAG = "consider_anticoag"
    ANTICOAG_RECOMMENDED = "anticoag_recommended"
    ANTICOAG_CONTRAINDICATED = "anticoag_contraindicated"


class GDMTPillar(str, Enum):
    """Four pillars of guideline-directed medical therapy for HFrEF."""
    BETA_BLOCKER = "beta_blocker"
    ARNI_ACEI_ARB = "arni_acei_arb"
    MRA = "mra"
    SGLT2I = "sglt2i"


class GDMTStatus(str, Enum):
    """Current status of a GDMT medication for a patient."""
    NOT_STARTED = "not_started"
    INITIATED = "initiated"
    UPTITRATING = "uptitrating"
    AT_TARGET = "at_target"
    CONTRAINDICATED = "contraindicated"
    INTOLERANT = "intolerant"


class CardiotoxicityRisk(str, Enum):
    """Risk level for cancer-therapy-related cardiac dysfunction."""
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"
    VERY_HIGH = "very_high"


class LGEPattern(str, Enum):
    """Late gadolinium enhancement patterns on cardiac MRI."""
    SUBENDOCARDIAL = "subendocardial"
    MID_WALL = "mid_wall"
    EPICARDIAL = "epicardial"
    TRANSMURAL = "transmural"
    RV_INSERTION = "rv_insertion"
    PATCHY = "patchy"
    NONE = "none"


class GuidelineClass(str, Enum):
    """ACC/AHA/ESC guideline recommendation classification.

    Class I: Benefit >>> Risk (IS recommended)
    Class IIa: Benefit >> Risk (IS reasonable)
    Class IIb: Benefit >= Risk (MAY be considered)
    Class III No Benefit: Benefit = Risk (NOT recommended)
    Class III Harm: Risk > Benefit (potentially harmful)
    """
    CLASS_I = "class_i"
    CLASS_IIA = "class_iia"
    CLASS_IIB = "class_iib"
    CLASS_III_NO_BENEFIT = "class_iii_no_benefit"
    CLASS_III_HARM = "class_iii_harm"


class EvidenceLevel(str, Enum):
    """Level of evidence supporting a guideline recommendation.

    A: High-quality evidence from multiple RCTs or meta-analyses
    B-R: Moderate-quality evidence from one or more RCTs
    B-NR: Moderate-quality evidence from non-randomized studies
    C-LD: Limited data from observational or registry studies
    C-EO: Expert opinion / consensus
    """
    LEVEL_A = "level_a"
    LEVEL_B_R = "level_b_r"
    LEVEL_B_NR = "level_b_nr"
    LEVEL_C_LD = "level_c_ld"
    LEVEL_C_EO = "level_c_eo"


# ===================================================================
# PYDANTIC MODELS - QUERY & SEARCH
# ===================================================================


class CardioQuery(BaseModel):
    """Input query to the Cardiology Intelligence Agent."""
    question: str = Field(..., min_length=1, description="Clinical question")
    workflow_type: Optional[CardioWorkflowType] = Field(
        default=None,
        description="Specific workflow to route the query to; auto-detected if omitted",
    )
    patient_context: Optional[Dict] = Field(
        default=None,
        description="Patient demographics, vitals, labs, and medications",
    )


class CardioSearchResult(BaseModel):
    """A single search result from any cardiology knowledge collection."""
    collection: str = Field(..., description="Source Milvus collection name")
    content: str = Field(..., description="Retrieved text content")
    score: float = Field(..., ge=0.0, description="Similarity score")
    metadata: Dict = Field(default_factory=dict, description="Source metadata")


# ===================================================================
# PYDANTIC MODELS - RISK SCORING
# ===================================================================


class RiskScoreInput(BaseModel):
    """Input parameters for cardiovascular risk score calculators.

    Not all fields are required for every calculator. Each score_type uses
    the subset of fields relevant to its formula.
    """
    score_type: RiskScoreType = Field(..., description="Which risk calculator to use")

    # Demographics
    age: Optional[int] = Field(
        default=None, ge=18, le=120,
        description="Patient age in years",
    )
    sex: Optional[str] = Field(
        default=None,
        description="Patient sex: male or female",
    )
    race: Optional[str] = Field(
        default=None,
        description="Race/ethnicity for ASCVD calculator",
    )

    # Vitals
    systolic_bp: Optional[float] = Field(
        default=None, ge=60, le=300,
        description="Systolic blood pressure (mmHg)",
    )
    heart_rate: Optional[int] = Field(
        default=None, ge=20, le=300,
        description="Heart rate (bpm)",
    )

    # Lipids
    total_cholesterol: Optional[float] = Field(
        default=None, ge=50, le=600,
        description="Total cholesterol (mg/dL)",
    )
    hdl: Optional[float] = Field(
        default=None, ge=10, le=200,
        description="HDL cholesterol (mg/dL)",
    )
    ldl: Optional[float] = Field(
        default=None, ge=20, le=500,
        description="LDL cholesterol (mg/dL)",
    )

    # Labs
    creatinine: Optional[float] = Field(
        default=None, ge=0.1, le=20.0,
        description="Serum creatinine (mg/dL)",
    )
    bnp: Optional[float] = Field(
        default=None, ge=0.0,
        description="BNP (pg/mL)",
    )
    nt_pro_bnp: Optional[float] = Field(
        default=None, ge=0.0,
        description="NT-proBNP (pg/mL)",
    )
    troponin: Optional[float] = Field(
        default=None, ge=0.0,
        description="Troponin (ng/mL)",
    )
    hba1c: Optional[float] = Field(
        default=None, ge=3.0, le=20.0,
        description="HbA1c (%)",
    )
    hemoglobin: Optional[float] = Field(
        default=None, ge=3.0, le=25.0,
        description="Hemoglobin (g/dL)",
    )
    inr: Optional[float] = Field(
        default=None, ge=0.5, le=15.0,
        description="INR",
    )

    # Comorbidities and risk factors
    diabetes: Optional[bool] = Field(default=None, description="Diabetes mellitus")
    smoker: Optional[bool] = Field(default=None, description="Current smoker")
    hypertension_treatment: Optional[bool] = Field(
        default=None,
        description="On antihypertensive therapy",
    )
    history_of_stroke: Optional[bool] = Field(
        default=None,
        description="Prior stroke or TIA",
    )
    history_of_bleeding: Optional[bool] = Field(
        default=None,
        description="Prior major bleeding event",
    )
    congestive_heart_failure: Optional[bool] = Field(
        default=None,
        description="CHF diagnosis",
    )
    vascular_disease: Optional[bool] = Field(
        default=None,
        description="Peripheral or aortic vascular disease",
    )
    renal_disease: Optional[bool] = Field(
        default=None,
        description="Renal dysfunction",
    )
    liver_disease: Optional[bool] = Field(
        default=None,
        description="Liver dysfunction",
    )
    labile_inr: Optional[bool] = Field(
        default=None,
        description="Labile INR (time in therapeutic range < 60%)",
    )
    alcohol_excess: Optional[bool] = Field(
        default=None,
        description="Excess alcohol intake",
    )
    antiplatelet_nsaid: Optional[bool] = Field(
        default=None,
        description="On antiplatelets or NSAIDs",
    )

    # Heart failure specific
    lvef: Optional[float] = Field(
        default=None, ge=5.0, le=90.0,
        description="Left ventricular ejection fraction (%)",
    )
    nyha_class: Optional[HeartFailureClass] = Field(
        default=None,
        description="NYHA functional class",
    )
    beta_blocker_use: Optional[bool] = Field(
        default=None,
        description="Currently on beta-blocker",
    )
    acei_arb_use: Optional[bool] = Field(
        default=None,
        description="Currently on ACEi/ARB/ARNI",
    )
    bmi: Optional[float] = Field(
        default=None, ge=10.0, le=80.0,
        description="Body mass index (kg/m2)",
    )

    # Surgical risk (EuroSCORE II)
    procedure_type: Optional[str] = Field(
        default=None,
        description="Surgical procedure type",
    )
    urgency: Optional[str] = Field(
        default=None,
        description="Procedure urgency: elective, urgent, emergent, salvage",
    )
    redo_surgery: Optional[bool] = Field(
        default=None,
        description="Redo cardiac surgery",
    )

    @field_validator("sex")
    @classmethod
    def validate_sex(cls, v: Optional[str]) -> Optional[str]:
        if v is not None and v.lower() not in ("male", "female"):
            raise ValueError("sex must be 'male' or 'female'")
        return v.lower() if v else v


class RiskScoreResult(BaseModel):
    """Output from a cardiovascular risk score calculation."""
    score_type: RiskScoreType = Field(..., description="Calculator used")
    score_value: float = Field(..., description="Computed risk score or percentage")
    risk_category: str = Field(
        ...,
        description="Risk stratification category (e.g., low, intermediate, high)",
    )
    interpretation: str = Field(
        ...,
        description="Clinical interpretation of the score",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Actionable clinical recommendations",
    )
    guideline_reference: str = Field(
        default="",
        description="Source guideline citation",
    )


# ===================================================================
# PYDANTIC MODELS - HEART FAILURE & GDMT
# ===================================================================


class GDMTMedication(BaseModel):
    """A single GDMT pillar medication and its current status."""
    pillar: GDMTPillar = Field(..., description="GDMT pillar category")
    drug_name: str = Field(
        ..., max_length=200,
        description="Specific drug name (e.g., sacubitril/valsartan)",
    )
    current_dose: Optional[str] = Field(
        default=None, max_length=100,
        description="Current dose if started",
    )
    target_dose: str = Field(
        ..., max_length=100,
        description="Guideline-recommended target dose",
    )
    status: GDMTStatus = Field(
        default=GDMTStatus.NOT_STARTED,
        description="Current titration status",
    )
    contraindications: List[str] = Field(
        default_factory=list,
        description="Active contraindications for this medication",
    )


class GDMTRecommendation(BaseModel):
    """GDMT optimization recommendation for a heart failure patient."""
    ef_category: EjectionFractionCategory = Field(
        ...,
        description="Patient's EF classification",
    )
    current_meds: List[GDMTMedication] = Field(
        default_factory=list,
        description="Current GDMT medication list with titration status",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Specific GDMT optimization steps",
    )
    next_steps: List[str] = Field(
        default_factory=list,
        description="Follow-up actions and monitoring plan",
    )
    guideline_references: List[str] = Field(
        default_factory=list,
        description="Supporting guideline citations",
    )


# ===================================================================
# PYDANTIC MODELS - VALVULAR DISEASE
# ===================================================================


class ValveAssessment(BaseModel):
    """Assessment of a single cardiac valve lesion."""
    valve: str = Field(
        ..., max_length=50,
        description="Valve name (e.g., aortic, mitral, tricuspid, pulmonic)",
    )
    pathology: str = Field(
        ..., max_length=200,
        description="Valve pathology (e.g., aortic stenosis, mitral regurgitation)",
    )
    severity: ValveSeverity = Field(..., description="Severity grade")
    measurements: Dict = Field(
        default_factory=dict,
        description="Key measurements (e.g., AVA, mean gradient, EROA, regurgitant volume)",
    )
    intervention_criteria_met: bool = Field(
        False,
        description="Whether guideline criteria for intervention are met",
    )
    recommendation: str = Field(
        default="",
        max_length=2000,
        description="Clinical recommendation based on severity and criteria",
    )


# ===================================================================
# PYDANTIC MODELS - ECG & DIAGNOSTICS
# ===================================================================


class ECGInterpretation(BaseModel):
    """Structured ECG interpretation result."""
    rhythm: str = Field(
        ..., max_length=200,
        description="Rhythm interpretation (e.g., normal sinus rhythm, atrial fibrillation)",
    )
    rate: int = Field(..., ge=0, le=400, description="Ventricular rate (bpm)")
    intervals: Dict = Field(
        default_factory=dict,
        description="Key intervals: PR (ms), QRS (ms), QTc (ms)",
    )
    axis: str = Field(
        default="", max_length=100,
        description="Electrical axis (e.g., normal, LAD, RAD)",
    )
    findings: List[str] = Field(
        default_factory=list,
        description="List of ECG findings (e.g., ST elevation, LVH, LBBB)",
    )
    urgency: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Clinical urgency of the findings",
    )


# ===================================================================
# PYDANTIC MODELS - IMAGING
# ===================================================================


class ImagingResult(BaseModel):
    """Structured cardiac imaging result across modalities."""
    modality: ImagingModality = Field(..., description="Imaging modality used")
    findings: List[str] = Field(
        default_factory=list,
        description="Key imaging findings",
    )
    measurements: Dict = Field(
        default_factory=dict,
        description="Quantitative measurements (e.g., LVEF, LV dimensions, wall motion)",
    )
    impression: str = Field(
        default="",
        max_length=3000,
        description="Overall impression / summary",
    )
    cross_modal_triggers: List[str] = Field(
        default_factory=list,
        description="Findings that should trigger additional imaging or genomic workup",
    )


# ===================================================================
# PYDANTIC MODELS - CARDIO-ONCOLOGY
# ===================================================================


class CardiotoxicityAssessment(BaseModel):
    """Cardiotoxicity surveillance assessment for oncology patients."""
    agent: str = Field(
        ..., max_length=200,
        description="Cancer therapy agent (e.g., doxorubicin, trastuzumab)",
    )
    risk_level: CardiotoxicityRisk = Field(
        ...,
        description="Baseline cardiotoxicity risk level",
    )
    baseline_lvef: Optional[float] = Field(
        default=None, ge=5.0, le=90.0,
        description="Baseline LVEF (%) before treatment",
    )
    current_lvef: Optional[float] = Field(
        default=None, ge=5.0, le=90.0,
        description="Most recent LVEF (%)",
    )
    gls_decline_percent: Optional[float] = Field(
        default=None, ge=0.0, le=100.0,
        description="Relative decline in global longitudinal strain (%)",
    )
    troponin_trend: Optional[str] = Field(
        default=None, max_length=200,
        description="Troponin trajectory (e.g., stable, rising, elevated)",
    )
    recommendation: str = Field(
        default="",
        max_length=2000,
        description="Cardioprotective and monitoring recommendations",
    )
    monitoring_schedule: str = Field(
        default="",
        max_length=500,
        description="Recommended follow-up imaging and lab schedule",
    )


# ===================================================================
# PYDANTIC MODELS - WORKFLOW & CROSS-MODAL
# ===================================================================


class CrossModalTrigger(BaseModel):
    """A trigger linking cardiac findings to genomic or multi-modal workup."""
    trigger_source: str = Field(
        ..., max_length=200,
        description="Source of the trigger (e.g., cardiac MRI, echocardiogram)",
    )
    finding: str = Field(
        ..., max_length=500,
        description="The clinical finding that initiated the trigger",
    )
    gene_panel: List[str] = Field(
        default_factory=list,
        description="Recommended gene panel (e.g., MYH7, MYBPC3, TNNT2, LMNA, TTN)",
    )
    conditions: List[str] = Field(
        default_factory=list,
        description="Associated conditions (e.g., HCM, DCM, ARVC, LQTS)",
    )
    rationale: str = Field(
        default="",
        max_length=2000,
        description="Clinical rationale for the cross-modal referral",
    )


class WorkflowResult(BaseModel):
    """Output from a single cardiology workflow execution."""
    workflow_type: CardioWorkflowType = Field(
        ...,
        description="Workflow that generated these results",
    )
    findings: List[str] = Field(
        default_factory=list,
        description="Key clinical findings",
    )
    risk_scores: List[RiskScoreResult] = Field(
        default_factory=list,
        description="Computed risk scores if applicable",
    )
    recommendations: List[str] = Field(
        default_factory=list,
        description="Clinical recommendations",
    )
    guideline_references: List[str] = Field(
        default_factory=list,
        description="Supporting guideline citations",
    )
    severity: SeverityLevel = Field(
        default=SeverityLevel.INFORMATIONAL,
        description="Overall severity of findings",
    )
    cross_modal_triggers: List[str] = Field(
        default_factory=list,
        description="Triggers for additional cross-modal workup",
    )


# ===================================================================
# DATACLASS - SEARCH PLAN
# ===================================================================


@dataclass
class SearchPlan:
    """Pre-retrieval search planning for a cardiology query.

    Built by the query analyzer to guide collection routing, sub-question
    decomposition, and search-strategy selection before retrieval.
    """
    question: str = ""
    conditions: List[str] = field(default_factory=list)
    drugs: List[str] = field(default_factory=list)
    imaging_modalities: List[str] = field(default_factory=list)
    relevant_workflows: List[CardioWorkflowType] = field(default_factory=list)
    search_strategy: str = "broad"  # broad | targeted | comparative | clinical
    sub_questions: List[str] = field(default_factory=list)
    identified_topics: List[str] = field(default_factory=list)


# ===================================================================
# PYDANTIC MODELS - AGENT RESPONSE
# ===================================================================


class CardioResponse(BaseModel):
    """Top-level output from the Cardiology Intelligence Agent."""
    answer: str = Field(..., description="Synthesized clinical answer")
    citations: List[Dict] = Field(
        default_factory=list,
        description="Source citations (collection, id, title, score)",
    )
    risk_scores: List[RiskScoreResult] = Field(
        default_factory=list,
        description="Computed risk scores",
    )
    workflow_results: List[WorkflowResult] = Field(
        default_factory=list,
        description="Results from each executed workflow",
    )
    cross_modal_triggers: List[CrossModalTrigger] = Field(
        default_factory=list,
        description="Cross-modal referral triggers",
    )
    confidence: float = Field(
        default=0.0, ge=0.0, le=1.0,
        description="Agent confidence in the response (0.0 - 1.0)",
    )
