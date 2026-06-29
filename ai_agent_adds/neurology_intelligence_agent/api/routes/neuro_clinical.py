"""Neurology clinical API routes.

Provides endpoints for RAG-powered neurology queries, clinical scale
calculators (NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6, ALSFRS-R,
ASPECTS, Hoehn-Yahr), acute stroke triage, dementia evaluation,
epilepsy classification, brain tumor grading, MS assessment,
Parkinson's assessment, headache classification, neuromuscular
evaluation, and reference catalogues.

Author: Adam Jones
Date: March 2026
"""

from typing import List, Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

from src.models import (
    NeuroWorkflowType,
    ClinicalScaleType,
)
from src.knowledge import KNOWLEDGE_VERSION

router = APIRouter(prefix="/v1/neuro", tags=["neurology"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/integrated-assessment")
async def integrated_assessment(request: dict, req: Request):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries imaging, cardiology, biomarker, oncology, and rare disease agents
    for a comprehensive neurological assessment.
    """
    try:
        from src.cross_modal import (
            query_imaging_agent,
            query_cardiology_agent,
            query_biomarker_agent,
            query_oncology_agent,
            query_rare_disease_agent,
            integrate_cross_agent_results,
        )

        patient_profile = request.get("patient_profile", {})
        imaging_data = request.get("imaging_data", {})
        biomarkers = request.get("biomarkers", [])
        tumor_data = request.get("tumor_data", {})

        results = []

        # Query imaging agent for neuroimaging correlation
        if imaging_data:
            results.append(query_imaging_agent(imaging_data))

        # Query cardiology agent for cardiac/stroke overlap
        if patient_profile:
            results.append(query_cardiology_agent(patient_profile))

        # Query biomarker agent for fluid biomarker enrichment
        if biomarkers:
            results.append(query_biomarker_agent(biomarkers))

        # Query oncology agent for neuro-oncology context
        if tumor_data:
            results.append(query_oncology_agent(tumor_data))

        # Query rare disease agent for neurogenetic correlation
        if patient_profile.get("phenotypes") or patient_profile.get("genomic_variants"):
            results.append(query_rare_disease_agent(patient_profile))

        integrated = integrate_cross_agent_results(results)
        return {
            "status": "completed",
            "assessment": integrated,
            "agents_consulted": integrated.get("agents_consulted", []),
        }
    except Exception as exc:
        logger.error(f"Integrated assessment failed: {exc}")
        return {"status": "partial", "assessment": {}, "error": "Cross-agent integration unavailable"}


# =====================================================================
# Request / Response Schemas
# =====================================================================

# -- Query --

class QueryRequest(BaseModel):
    """Free-text RAG query with optional domain and patient context."""
    question: str = Field(..., min_length=3, description="Neurology question")
    domain: Optional[str] = Field(
        None,
        description=(
            "Domain hint: stroke | dementia | epilepsy | tumors | ms | "
            "parkinsons | headache | neuromuscular | general"
        ),
    )
    patient_context: Optional[dict] = Field(None, description="Demographics, imaging, labs")
    top_k: int = Field(5, ge=1, le=50, description="Number of evidence passages")
    include_guidelines: bool = Field(True, description="Include guideline citations")


class QueryResponse(BaseModel):
    answer: str
    evidence: List[dict]
    guidelines_cited: List[str] = []
    confidence: float
    domain_applied: Optional[str] = None


class SearchRequest(BaseModel):
    """Multi-collection semantic search."""
    question: str = Field(..., min_length=3)
    collections: Optional[List[str]] = None
    top_k: int = Field(5, ge=1, le=100)
    threshold: float = Field(0.0, ge=0.0, le=1.0)


class SearchResult(BaseModel):
    collection: str
    text: str
    score: float
    metadata: dict = {}


class SearchResponse(BaseModel):
    results: List[SearchResult]
    total: int
    collections_searched: List[str]


# -- Clinical Scales --

class ScaleCalculateRequest(BaseModel):
    """Clinical scale calculation request."""
    scale_name: str = Field(
        ...,
        description=(
            "Scale identifier: nihss | gcs | moca | updrs | edss | "
            "mrs | hit6 | alsfrs | aspects | hoehn_yahr"
        ),
    )
    items: dict = Field(
        ...,
        description="Item-level scores as key-value pairs specific to each scale",
    )


class ScaleCalculateResponse(BaseModel):
    scale_name: str
    total_score: float
    max_score: float
    interpretation: str
    severity_category: str
    thresholds: dict = {}
    recommendations: List[str] = []
    items_received: int


# -- Stroke Triage --

class StrokeTriageRequest(BaseModel):
    """Acute stroke triage request."""
    nihss_items: Optional[dict] = Field(None, description="NIHSS item scores")
    nihss_total: Optional[int] = Field(None, ge=0, le=42, description="Pre-computed NIHSS total")
    onset_time_hours: Optional[float] = Field(None, ge=0, description="Hours since symptom onset")
    ct_aspects_score: Optional[int] = Field(None, ge=0, le=10, description="CT ASPECTS score")
    lvo_suspected: Optional[bool] = Field(None, description="Large vessel occlusion suspected")
    age: Optional[int] = Field(None, ge=0, le=120)
    anticoagulant_use: Optional[bool] = Field(None, description="Current anticoagulant use")
    blood_pressure_systolic: Optional[int] = Field(None, ge=0, le=300)
    blood_glucose: Optional[float] = Field(None, ge=0, description="Blood glucose mg/dL")
    prior_mrs: Optional[int] = Field(None, ge=0, le=5, description="Pre-stroke mRS score")
    clinical_notes: Optional[str] = Field(None, description="Free-text clinical notes")


class StrokeTriageResponse(BaseModel):
    nihss_score: int
    nihss_severity: str
    tpa_eligible: bool
    tpa_reasoning: List[str]
    thrombectomy_eligible: bool
    thrombectomy_reasoning: List[str]
    aspects_interpretation: Optional[str] = None
    stroke_type_suggestion: str
    urgency_level: str
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Dementia Evaluation --

class DementiaEvaluateRequest(BaseModel):
    """Dementia evaluation request."""
    moca_items: Optional[dict] = Field(None, description="MoCA item scores")
    moca_total: Optional[int] = Field(None, ge=0, le=30, description="Pre-computed MoCA total")
    age: Optional[int] = Field(None, ge=0, le=120)
    education_years: Optional[int] = Field(None, ge=0, le=30)
    symptom_duration_months: Optional[int] = Field(None, ge=0)
    dominant_symptoms: List[str] = Field(default_factory=list, description="e.g., memory_loss, behavioral_change, visual_hallucinations")
    motor_features: Optional[str] = Field(None, description="Parkinsonism, gait, tremor")
    biomarkers: Optional[dict] = Field(None, description="CSF amyloid/tau, PET, MRI findings")
    family_history: Optional[str] = Field(None)
    medications: List[str] = Field(default_factory=list)
    clinical_notes: Optional[str] = Field(None)


class DementiaEvaluateResponse(BaseModel):
    moca_score: Optional[int] = None
    moca_interpretation: Optional[str] = None
    cognitive_domain_analysis: dict
    differential_diagnosis: List[dict]
    atn_staging: Optional[str] = None
    recommended_workup: List[str]
    treatment_options: List[str]
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Epilepsy Classification --

class EpilepsyClassifyRequest(BaseModel):
    """Epilepsy classification request."""
    seizure_description: str = Field(..., min_length=5, description="Seizure semiology description")
    age_at_onset: Optional[int] = Field(None, ge=0, le=120)
    eeg_findings: Optional[str] = Field(None, description="EEG description or conclusion")
    mri_findings: Optional[str] = Field(None, description="Brain MRI findings")
    seizure_frequency: Optional[str] = Field(None, description="Frequency (e.g., 2/month)")
    current_aeds: List[str] = Field(default_factory=list, description="Current anti-epileptic drugs")
    family_history: Optional[str] = Field(None)
    genetic_results: Optional[str] = Field(None, description="Known genetic variants")
    clinical_notes: Optional[str] = Field(None)


class EpilepsyClassifyResponse(BaseModel):
    seizure_type: str
    epilepsy_syndrome: Optional[str] = None
    focal_features: List[str]
    generalized_features: List[str]
    eeg_mri_concordance: Optional[str] = None
    surgical_candidacy: Optional[str] = None
    aed_recommendations: List[str]
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Tumor Grading --

class TumorGradeRequest(BaseModel):
    """Brain tumor grading request."""
    histology: str = Field(..., min_length=3, description="Histological diagnosis")
    molecular_markers: List[str] = Field(default_factory=list, description="e.g., idh_mutant, mgmt_methylated, 1p19q_codeleted")
    location: Optional[str] = Field(None, description="Tumor location")
    age: Optional[int] = Field(None, ge=0, le=120)
    kps: Optional[int] = Field(None, ge=0, le=100, description="Karnofsky Performance Status")
    extent_of_resection: Optional[str] = Field(None, description="GTR | STR | biopsy")
    imaging_features: Optional[str] = Field(None, description="MRI enhancement, edema, etc.")
    clinical_notes: Optional[str] = Field(None)


class TumorGradeResponse(BaseModel):
    who_grade: str
    integrated_diagnosis: str
    molecular_profile: dict
    treatment_recommendations: List[str]
    prognosis_summary: str
    nccn_category: Optional[str] = None
    clinical_trial_relevance: List[str]
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- MS Assessment --

class MSAssessRequest(BaseModel):
    """MS disease assessment request."""
    ms_course: Optional[str] = Field(None, description="rrms | spms | ppms | cis | ris")
    edss_score: Optional[float] = Field(None, ge=0.0, le=10.0, description="EDSS score")
    edss_items: Optional[dict] = Field(None, description="EDSS functional system scores")
    relapse_count_2yr: Optional[int] = Field(None, ge=0, description="Relapses in past 2 years")
    new_t2_lesions: Optional[int] = Field(None, ge=0, description="New T2 lesions on MRI")
    gad_enhancing_lesions: Optional[int] = Field(None, ge=0, description="Gad-enhancing lesions")
    current_dmt: Optional[str] = Field(None, description="Current DMT name")
    disease_duration_years: Optional[float] = Field(None, ge=0)
    age: Optional[int] = Field(None, ge=0, le=120)
    jcv_status: Optional[str] = Field(None, description="JCV antibody status: positive | negative | unknown")
    clinical_notes: Optional[str] = Field(None)


class MSAssessResponse(BaseModel):
    edss_score: Optional[float] = None
    edss_interpretation: Optional[str] = None
    disease_activity: str
    neda_status: dict
    dmt_assessment: dict
    escalation_recommendation: Optional[str] = None
    mri_disease_burden: Optional[str] = None
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Parkinson's Assessment --

class ParkinsonsAssessRequest(BaseModel):
    """Parkinson's disease assessment request."""
    updrs_items: Optional[dict] = Field(None, description="MDS-UPDRS Part III item scores")
    updrs_total: Optional[int] = Field(None, ge=0, le=132, description="Pre-computed UPDRS-III total")
    hoehn_yahr: Optional[float] = Field(None, ge=0, le=5, description="Hoehn-Yahr stage")
    age: Optional[int] = Field(None, ge=0, le=120)
    disease_duration_years: Optional[float] = Field(None, ge=0)
    current_medications: List[str] = Field(default_factory=list)
    motor_fluctuations: Optional[bool] = Field(None, description="Motor fluctuations present")
    dyskinesia: Optional[bool] = Field(None, description="Dyskinesia present")
    non_motor_symptoms: List[str] = Field(default_factory=list, description="e.g., RBD, constipation, orthostatic_hypotension")
    dat_scan: Optional[str] = Field(None, description="DaT-SPECT result")
    clinical_notes: Optional[str] = Field(None)


class ParkinsonsAssessResponse(BaseModel):
    updrs_score: Optional[int] = None
    updrs_interpretation: Optional[str] = None
    hoehn_yahr_stage: Optional[float] = None
    hoehn_yahr_description: Optional[str] = None
    motor_subtype: str
    medication_assessment: dict
    dbs_candidacy: Optional[str] = None
    non_motor_management: List[str]
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Headache Classification --

class HeadacheClassifyRequest(BaseModel):
    """Headache classification request."""
    headache_description: str = Field(..., min_length=5, description="Headache characteristics")
    duration_hours: Optional[float] = Field(None, ge=0, description="Typical duration in hours")
    frequency_per_month: Optional[int] = Field(None, ge=0, description="Episodes per month")
    location: Optional[str] = Field(None, description="Unilateral | bilateral | occipital etc.")
    quality: Optional[str] = Field(None, description="Throbbing | pressure | stabbing etc.")
    associated_symptoms: List[str] = Field(default_factory=list, description="e.g., nausea, photophobia, tearing")
    aura: Optional[bool] = Field(None, description="Aura present")
    triggers: List[str] = Field(default_factory=list, description="Known triggers")
    current_medications: List[str] = Field(default_factory=list)
    analgesic_days_per_month: Optional[int] = Field(None, ge=0, description="Days/month using acute meds")
    red_flags: List[str] = Field(default_factory=list, description="Red flag symptoms present")
    hit6_items: Optional[dict] = Field(None, description="HIT-6 item scores")
    clinical_notes: Optional[str] = Field(None)


class HeadacheClassifyResponse(BaseModel):
    headache_type: str
    ichd3_criteria_met: List[str]
    hit6_score: Optional[int] = None
    hit6_interpretation: Optional[str] = None
    red_flag_assessment: dict
    medication_overuse_risk: Optional[str] = None
    acute_treatment: List[str]
    preventive_treatment: List[str]
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Neuromuscular Evaluation --

class NeuromuscularEvaluateRequest(BaseModel):
    """Neuromuscular evaluation request."""
    presentation: str = Field(..., min_length=5, description="Clinical presentation")
    weakness_pattern: Optional[str] = Field(None, description="Proximal | distal | bulbar | respiratory")
    sensory_involvement: Optional[bool] = Field(None)
    emg_findings: Optional[str] = Field(None, description="EMG/NCS summary")
    ck_level: Optional[float] = Field(None, ge=0, description="Creatine kinase level IU/L")
    antibodies: List[str] = Field(default_factory=list, description="e.g., AChR, MuSK, anti-GM1")
    alsfrs_items: Optional[dict] = Field(None, description="ALSFRS-R item scores")
    age: Optional[int] = Field(None, ge=0, le=120)
    progression_rate: Optional[str] = Field(None, description="Rapid | moderate | slow")
    genetic_results: Optional[str] = Field(None)
    clinical_notes: Optional[str] = Field(None)


class NeuromuscularEvaluateResponse(BaseModel):
    localization: str
    pattern_classification: str
    differential_diagnosis: List[dict]
    alsfrs_score: Optional[int] = None
    alsfrs_interpretation: Optional[str] = None
    recommended_workup: List[str]
    treatment_options: List[str]
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []


# -- Generic Workflow --

class WorkflowRequest(BaseModel):
    """Generic workflow dispatch request."""
    data: dict = Field(default_factory=dict, description="Workflow input data")
    question: Optional[str] = Field(None, description="Optional guiding question")


class WorkflowResponse(BaseModel):
    workflow_type: str
    status: str
    result: str
    evidence_used: bool = False
    note: Optional[str] = None


# =====================================================================
# Clinical Scale Calculators
# =====================================================================

def _calculate_nihss(items: dict) -> dict:
    """Calculate NIHSS total from item scores.

    Items: 1a_loc, 1b_loc_questions, 1c_loc_commands, 2_gaze,
    3_visual, 4_facial_palsy, 5a_motor_left_arm, 5b_motor_right_arm,
    6a_motor_left_leg, 6b_motor_right_leg, 7_limb_ataxia,
    8_sensory, 9_language, 10_dysarthria, 11_extinction
    """
    total = sum(int(v) for v in items.values() if str(v).isdigit())
    if total == 0:
        severity, category = "No stroke symptoms", "none"
    elif total <= 4:
        severity, category = "Minor stroke", "minor"
    elif total <= 15:
        severity, category = "Moderate stroke", "moderate"
    elif total <= 20:
        severity, category = "Moderate-to-severe stroke", "moderate_severe"
    elif total <= 25:
        severity, category = "Severe stroke", "severe"
    else:
        severity, category = "Very severe stroke", "very_severe"

    recs = []
    if total >= 6:
        recs.append("Consider IV tPA if within 4.5h window and no contraindications (AHA/ASA Class I)")
    if total >= 6 and total <= 25:
        recs.append("Evaluate for mechanical thrombectomy if LVO suspected and within 24h (AHA/ASA Class I)")
    if total >= 1:
        recs.append("Admit to stroke unit for monitoring (AHA/ASA Class I)")

    return {
        "total_score": total, "max_score": 42,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"minor": 4, "moderate": 15, "severe": 25, "tpa_threshold": 6},
        "recommendations": recs,
    }


def _calculate_gcs(items: dict) -> dict:
    """Calculate GCS from eye, verbal, motor components."""
    eye = int(items.get("eye", items.get("e", 1)))
    verbal = int(items.get("verbal", items.get("v", 1)))
    motor = int(items.get("motor", items.get("m", 1)))
    total = eye + verbal + motor

    if total <= 8:
        severity, category = "Severe brain injury -- intubation likely required", "severe"
    elif total <= 12:
        severity, category = "Moderate brain injury", "moderate"
    else:
        severity, category = "Mild brain injury", "mild"

    recs = []
    if total <= 8:
        recs.append("Secure airway, consider intubation")
        recs.append("Urgent neurosurgical consultation")
        recs.append("CT head stat")
    elif total <= 12:
        recs.append("Close neurological monitoring")
        recs.append("CT head to evaluate for intracranial pathology")

    return {
        "total_score": total, "max_score": 15,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"severe": 8, "moderate": 12, "mild": 15},
        "recommendations": recs,
    }


def _calculate_moca(items: dict) -> dict:
    """Calculate MoCA total from domain item scores.

    Domains: visuospatial, naming, attention, language, abstraction,
    delayed_recall, orientation. Education correction applied.
    """
    total = sum(int(v) for v in items.values() if str(v).isdigit())
    # Education correction: +1 point if education <= 12 years
    education = items.get("education_years")
    if education is not None and int(education) <= 12:
        total = min(total + 1, 30)

    if total >= 26:
        severity, category = "Normal cognition", "normal"
    elif total >= 18:
        severity, category = "Mild cognitive impairment (MCI)", "mci"
    elif total >= 10:
        severity, category = "Moderate cognitive impairment", "moderate"
    else:
        severity, category = "Severe cognitive impairment", "severe"

    recs = []
    if total < 26:
        recs.append("Comprehensive neuropsychological testing recommended")
        recs.append("Consider neuroimaging (MRI brain with volumetrics)")
        recs.append("Screen for reversible causes: B12, TSH, RPR, HIV")
    if total < 18:
        recs.append("Evaluate for cholinesterase inhibitor therapy")
        recs.append("Assess functional independence and safety")

    return {
        "total_score": total, "max_score": 30,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"normal": 26, "mci": 18, "moderate": 10},
        "recommendations": recs,
    }


def _calculate_updrs(items: dict) -> dict:
    """Calculate MDS-UPDRS Part III motor score from items."""
    total = sum(int(v) for v in items.values() if str(v).isdigit())

    if total <= 10:
        severity, category = "Minimal motor impairment", "minimal"
    elif total <= 32:
        severity, category = "Mild motor impairment", "mild"
    elif total <= 58:
        severity, category = "Moderate motor impairment", "moderate"
    else:
        severity, category = "Severe motor impairment", "severe"

    # Determine motor subtype from item patterns
    sum(int(items.get(k, 0)) for k in ["3.15a", "3.15b", "3.16a", "3.16b", "3.17a", "3.17b", "3.17c", "3.17d", "3.17e", "3.18"] if k in items)
    sum(int(items.get(k, 0)) for k in ["3.10", "3.11", "3.12", "3.13"] if k in items)

    recs = []
    if total > 32:
        recs.append("Consider medication adjustment or addition")
    if total > 58:
        recs.append("Evaluate for deep brain stimulation candidacy")
        recs.append("Refer to movement disorders specialist")

    return {
        "total_score": total, "max_score": 132,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"minimal": 10, "mild": 32, "moderate": 58},
        "recommendations": recs,
    }


def _calculate_edss(items: dict) -> dict:
    """Calculate EDSS from functional system scores.

    Items: pyramidal, cerebellar, brainstem, sensory, bowel_bladder,
    visual, cerebral, ambulation
    """
    fs_scores = {k: int(v) for k, v in items.items() if k != "ambulation" and str(v).isdigit()}
    max_fs = max(fs_scores.values()) if fs_scores else 0
    ambulation = float(items.get("ambulation", 0))

    # Simplified EDSS estimation
    if ambulation > 0:
        edss = ambulation
    elif max_fs == 0:
        edss = 0.0
    elif max_fs == 1:
        edss = 1.0 if sum(1 for v in fs_scores.values() if v >= 1) <= 1 else 1.5
    elif max_fs == 2:
        edss = 2.0 if sum(1 for v in fs_scores.values() if v >= 2) <= 1 else 2.5
    elif max_fs == 3:
        edss = 3.0 if sum(1 for v in fs_scores.values() if v >= 3) <= 1 else 3.5
    elif max_fs == 4:
        edss = 4.0
    else:
        edss = min(float(max_fs), 9.5)

    if edss <= 1.5:
        severity, category = "Minimal disability -- fully ambulatory", "minimal"
    elif edss <= 3.5:
        severity, category = "Moderate disability -- fully ambulatory", "moderate"
    elif edss <= 5.5:
        severity, category = "Significant disability -- ambulatory with limitations", "significant"
    elif edss <= 7.5:
        severity, category = "Severe disability -- restricted to wheelchair", "severe"
    else:
        severity, category = "Very severe disability -- bed-bound", "very_severe"

    recs = []
    if edss >= 3.0:
        recs.append("Consider high-efficacy DMT if not already on one")
    if edss >= 6.0:
        recs.append("Evaluate for rehabilitation and assistive devices")
        recs.append("Assess for symptom management (spasticity, fatigue, pain)")

    return {
        "total_score": edss, "max_score": 10.0,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"minimal": 1.5, "moderate": 3.5, "significant": 5.5, "severe": 7.5},
        "recommendations": recs,
    }


def _calculate_mrs(items: dict) -> dict:
    """Calculate modified Rankin Scale score."""
    score = int(items.get("score", items.get("mrs", 0)))
    descriptions = {
        0: ("No symptoms at all", "no_disability"),
        1: ("No significant disability despite symptoms", "no_significant_disability"),
        2: ("Slight disability -- unable to carry out all prior activities but independent", "slight"),
        3: ("Moderate disability -- requires some help but walks without assistance", "moderate"),
        4: ("Moderately severe disability -- unable to walk without assistance", "moderately_severe"),
        5: ("Severe disability -- bedridden, incontinent, requires constant care", "severe"),
        6: ("Dead", "dead"),
    }
    score = max(0, min(score, 6))
    interpretation, category = descriptions[score]

    recs = []
    if score >= 3:
        recs.append("Arrange rehabilitation services")
        recs.append("Assess for secondary stroke prevention")
    if score >= 4:
        recs.append("Evaluate caregiver support needs")
        recs.append("Consider skilled nursing facility if unable to manage at home")

    return {
        "total_score": score, "max_score": 6,
        "interpretation": interpretation, "severity_category": category,
        "thresholds": {"good_outcome": 2, "moderate": 3, "severe": 5},
        "recommendations": recs,
    }


def _calculate_hit6(items: dict) -> dict:
    """Calculate HIT-6 (Headache Impact Test) score.

    6 items, each scored 6=never, 8=rarely, 10=sometimes, 11=very_often, 13=always.
    """
    total = sum(int(v) for v in items.values() if str(v).isdigit())

    if total <= 49:
        severity, category = "No significant impact", "minimal"
    elif total <= 55:
        severity, category = "Some impact -- occasional headache interference", "some"
    elif total <= 59:
        severity, category = "Substantial impact -- frequent interference with daily activities", "substantial"
    else:
        severity, category = "Severe impact -- headaches significantly limit daily functioning", "severe"

    recs = []
    if total >= 56:
        recs.append("Consider preventive therapy initiation")
    if total >= 60:
        recs.append("Evaluate for prophylactic medication (CGRP mAb, topiramate, or beta-blocker)")
        recs.append("Screen for medication overuse headache")
        recs.append("Consider headache specialist referral")

    return {
        "total_score": total, "max_score": 78,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"minimal": 49, "some": 55, "substantial": 59},
        "recommendations": recs,
    }


def _calculate_alsfrs(items: dict) -> dict:
    """Calculate ALSFRS-R (Revised ALS Functional Rating Scale).

    12 items, each scored 0-4 (4=normal function).
    """
    total = sum(int(v) for v in items.values() if str(v).isdigit())

    if total >= 40:
        severity, category = "Mild functional impairment", "mild"
    elif total >= 30:
        severity, category = "Moderate functional impairment", "moderate"
    elif total >= 20:
        severity, category = "Significant functional impairment", "significant"
    else:
        severity, category = "Severe functional impairment", "severe"

    recs = []
    if total < 40:
        recs.append("Multidisciplinary ALS clinic follow-up recommended")
    if total < 30:
        recs.append("Evaluate for non-invasive ventilation (FVC < 50%)")
        recs.append("Discuss advanced care planning")
    if total < 20:
        recs.append("Assess for PEG tube placement if swallowing affected")
        recs.append("Intensify respiratory monitoring")

    # Bulbar sub-score
    bulbar_items = sum(int(items.get(k, 4)) for k in ["speech", "salivation", "swallowing"] if k in items)
    if bulbar_items < 9:
        recs.append("Prominent bulbar involvement -- speech/swallow therapy referral")

    return {
        "total_score": total, "max_score": 48,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"mild": 40, "moderate": 30, "significant": 20},
        "recommendations": recs,
    }


def _calculate_aspects(items: dict) -> dict:
    """Calculate Alberta Stroke Program Early CT Score (ASPECTS).

    10 regions scored 1=normal, 0=early ischemic change. Total starts at 10.
    """
    affected_regions = sum(1 for v in items.values() if str(v) == "0")
    score = 10 - affected_regions

    if score >= 8:
        severity, category = "Small ischemic core -- favorable for intervention", "small_core"
    elif score >= 6:
        severity, category = "Moderate ischemic core", "moderate_core"
    else:
        severity, category = "Large ischemic core -- higher hemorrhage risk with thrombolysis", "large_core"

    recs = []
    if score >= 6:
        recs.append("Favorable for thrombolysis and/or thrombectomy (AHA/ASA)")
    else:
        recs.append("Large established infarct -- weigh risks/benefits of reperfusion therapy carefully")
        recs.append("Consider advanced imaging (CT perfusion) for salvageable tissue assessment")

    return {
        "total_score": score, "max_score": 10,
        "interpretation": severity, "severity_category": category,
        "thresholds": {"favorable": 8, "moderate": 6, "unfavorable": 5},
        "recommendations": recs,
    }


def _calculate_hoehn_yahr(items: dict) -> dict:
    """Calculate Hoehn-Yahr staging."""
    stage = float(items.get("stage", items.get("score", 0)))
    descriptions = {
        0.0: ("No signs of disease", "no_disease"),
        1.0: ("Unilateral involvement only", "stage_1"),
        1.5: ("Unilateral and axial involvement", "stage_1_5"),
        2.0: ("Bilateral involvement without balance impairment", "stage_2"),
        2.5: ("Mild bilateral disease with recovery on pull test", "stage_2_5"),
        3.0: ("Mild-to-moderate bilateral disease, postural instability, physically independent", "stage_3"),
        4.0: ("Severe disability but still able to walk or stand unassisted", "stage_4"),
        5.0: ("Wheelchair-bound or bedridden unless aided", "stage_5"),
    }
    stage = max(0.0, min(stage, 5.0))
    nearest = min(descriptions.keys(), key=lambda k: abs(k - stage))
    interpretation, category = descriptions[nearest]

    recs = []
    if stage >= 3.0:
        recs.append("Consider DBS evaluation if medication response is waning")
        recs.append("Physical therapy for fall prevention")
    if stage >= 4.0:
        recs.append("Evaluate home safety and caregiver support")
        recs.append("Assess for palliative care integration")

    return {
        "total_score": stage, "max_score": 5.0,
        "interpretation": interpretation, "severity_category": category,
        "thresholds": {"mild": 2.0, "moderate": 3.0, "severe": 4.0},
        "recommendations": recs,
    }


_SCALE_CALCULATORS = {
    "nihss": _calculate_nihss,
    "gcs": _calculate_gcs,
    "moca": _calculate_moca,
    "updrs": _calculate_updrs,
    "updrs_part_iii": _calculate_updrs,
    "edss": _calculate_edss,
    "mrs": _calculate_mrs,
    "hit6": _calculate_hit6,
    "alsfrs": _calculate_alsfrs,
    "alsfrs_r": _calculate_alsfrs,
    "aspects": _calculate_aspects,
    "hoehn_yahr": _calculate_hoehn_yahr,
}


# =====================================================================
# Helper: get engine and manager from app state
# =====================================================================

def _get_engine(request: Request):
    """Get RAG engine from app state, raise 503 if unavailable."""
    engine = getattr(request.app.state, "engine", None)
    if engine is None:
        raise HTTPException(
            status_code=503,
            detail="RAG engine is unavailable. Service is in degraded mode.",
        )
    return engine


def _get_workflow_engine(request: Request):
    """Get workflow engine from app state."""
    return getattr(request.app.state, "workflow_engine", None)


def _get_llm(request: Request):
    """Get LLM client from app state."""
    return getattr(request.app.state, "llm_client", None)


def _increment_metric(request: Request, metric: str):
    """Thread-safe metric increment."""
    metrics = getattr(request.app.state, "metrics", None)
    lock = getattr(request.app.state, "metrics_lock", None)
    if metrics and lock:
        with lock:
            metrics[metric] = metrics.get(metric, 0) + 1


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/query", response_model=QueryResponse)
async def neuro_query(request: QueryRequest, req: Request):
    """RAG-powered neurology Q&A.

    Searches across all neurology knowledge collections and synthesizes
    an evidence-based answer with citations and guideline references.
    """
    _increment_metric(req, "query_requests_total")
    engine = _get_engine(req)

    try:
        results = engine.search(request.question, top_k=request.top_k)
        evidence = [
            {
                "collection": r.get("collection", "unknown"),
                "text": r.get("content", r.get("text", "")),
                "score": r.get("score", 0.0),
                "metadata": r.get("metadata", {}),
            }
            for r in results
        ]
    except Exception as exc:
        logger.warning(f"Search failed: {exc}")
        evidence = []

    # Generate answer via LLM
    llm = _get_llm(req)
    answer = "Search completed. See evidence passages below."
    confidence = 0.5
    guidelines_cited = []

    if llm:
        context = "\n\n".join(e["text"] for e in evidence if e["text"])
        prompt = (
            f"Neurology question: {request.question}\n\n"
            f"Retrieved evidence:\n{context}\n\n"
            f"Provide a comprehensive answer citing specific evidence. "
            f"Include relevant AAN, AHA/ASA, ILAE, or IHS guideline references."
        )
        if request.patient_context:
            prompt += f"\n\nPatient context: {request.patient_context}"

        try:
            answer = llm.generate(prompt)
            confidence = min(0.95, 0.5 + len(evidence) * 0.05)
            # Extract guideline references from answer
            for gline in ["AAN Guideline", "AHA/ASA", "ILAE", "IHS", "ICHD-3",
                          "NCCN", "AAN Practice", "McDonald Criteria", "MDS-UPDRS"]:
                if gline.lower() in answer.lower():
                    guidelines_cited.append(gline)
        except Exception as exc:
            logger.warning(f"LLM generation failed: {exc}")

    return QueryResponse(
        answer=answer,
        evidence=evidence,
        guidelines_cited=guidelines_cited,
        confidence=confidence,
        domain_applied=request.domain,
    )


@router.post("/search", response_model=SearchResponse)
async def neuro_search(request: SearchRequest, req: Request):
    """Multi-collection semantic search across neurology knowledge base."""
    _increment_metric(req, "search_requests_total")
    engine = _get_engine(req)

    try:
        results = engine.search(
            request.question,
            top_k=request.top_k,
            collections=request.collections,
        )
        search_results = [
            SearchResult(
                collection=r.get("collection", "unknown"),
                text=r.get("content", r.get("text", "")),
                score=r.get("score", 0.0),
                metadata=r.get("metadata", {}),
            )
            for r in results
            if r.get("score", 0.0) >= request.threshold
        ]
    except Exception as exc:
        logger.warning(f"Search failed: {exc}")
        search_results = []

    collections_searched = list(set(r.collection for r in search_results))

    return SearchResponse(
        results=search_results,
        total=len(search_results),
        collections_searched=collections_searched,
    )


@router.post("/scale/calculate", response_model=ScaleCalculateResponse)
async def scale_calculate(request: ScaleCalculateRequest, req: Request):
    """Validated neurological clinical scale calculator.

    Supports: NIHSS, GCS, MoCA, UPDRS Part III, EDSS, mRS, HIT-6,
    ALSFRS-R, ASPECTS, Hoehn-Yahr.
    """
    _increment_metric(req, "scale_requests_total")

    scale_key = request.scale_name.lower().replace("-", "_")
    calculator = _SCALE_CALCULATORS.get(scale_key)
    if calculator is None:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Unknown scale '{request.scale_name}'. "
                f"Supported scales: {sorted(set(_SCALE_CALCULATORS.keys()))}"
            ),
        )

    try:
        result = calculator(request.items)
    except Exception as exc:
        logger.error(f"Scale calculation failed: {exc}")
        raise HTTPException(status_code=422, detail="Scale calculation error")

    return ScaleCalculateResponse(
        scale_name=request.scale_name,
        total_score=result["total_score"],
        max_score=result["max_score"],
        interpretation=result["interpretation"],
        severity_category=result["severity_category"],
        thresholds=result.get("thresholds", {}),
        recommendations=result.get("recommendations", []),
        items_received=len(request.items),
    )


@router.post("/stroke/triage", response_model=StrokeTriageResponse)
async def stroke_triage(request: StrokeTriageRequest, req: Request):
    """Acute stroke triage with NIHSS, thrombolysis/thrombectomy eligibility."""
    _increment_metric(req, "workflow_requests_total")

    # Calculate NIHSS
    nihss_score = request.nihss_total or 0
    nihss_severity = "unknown"
    if request.nihss_items:
        result = _calculate_nihss(request.nihss_items)
        nihss_score = int(result["total_score"])
        nihss_severity = result["severity_category"]
    elif request.nihss_total is not None:
        if nihss_score <= 4:
            nihss_severity = "minor"
        elif nihss_score <= 15:
            nihss_severity = "moderate"
        elif nihss_score <= 25:
            nihss_severity = "severe"
        else:
            nihss_severity = "very_severe"

    # tPA eligibility
    tpa_eligible = True
    tpa_reasoning = []
    if request.onset_time_hours is not None and request.onset_time_hours > 4.5:
        tpa_eligible = False
        tpa_reasoning.append(f"Onset {request.onset_time_hours:.1f}h ago -- exceeds 4.5h tPA window")
    elif request.onset_time_hours is not None:
        tpa_reasoning.append(f"Within tPA window ({request.onset_time_hours:.1f}h since onset)")
    if nihss_score < 4:
        tpa_reasoning.append(f"NIHSS {nihss_score} -- mild deficit, weigh risk/benefit")
    if request.anticoagulant_use:
        tpa_eligible = False
        tpa_reasoning.append("Anticoagulant use -- contraindication for tPA")
    if request.blood_pressure_systolic and request.blood_pressure_systolic > 185:
        tpa_reasoning.append(f"SBP {request.blood_pressure_systolic} -- must lower to <185 before tPA")
    if request.blood_glucose and (request.blood_glucose < 50 or request.blood_glucose > 400):
        tpa_reasoning.append(f"Blood glucose {request.blood_glucose} -- correct before tPA")

    # Thrombectomy eligibility
    thrombectomy_eligible = False
    thrombectomy_reasoning = []
    if nihss_score >= 6 and request.lvo_suspected:
        thrombectomy_eligible = True
        thrombectomy_reasoning.append("NIHSS >= 6 with suspected LVO -- thrombectomy candidate")
    elif request.lvo_suspected:
        thrombectomy_reasoning.append("LVO suspected but NIHSS may be low -- consider advanced imaging")
    else:
        thrombectomy_reasoning.append("No suspected LVO or NIHSS < 6")

    if request.onset_time_hours is not None and request.onset_time_hours > 24:
        thrombectomy_eligible = False
        thrombectomy_reasoning.append("Beyond 24h window for thrombectomy")
    if request.ct_aspects_score is not None and request.ct_aspects_score < 6:
        thrombectomy_eligible = False
        thrombectomy_reasoning.append(f"ASPECTS {request.ct_aspects_score} < 6 -- large established infarct")

    # ASPECTS interpretation
    aspects_interp = None
    if request.ct_aspects_score is not None:
        if request.ct_aspects_score >= 8:
            aspects_interp = f"ASPECTS {request.ct_aspects_score}: Small ischemic core, favorable"
        elif request.ct_aspects_score >= 6:
            aspects_interp = f"ASPECTS {request.ct_aspects_score}: Moderate ischemic core"
        else:
            aspects_interp = f"ASPECTS {request.ct_aspects_score}: Large ischemic core, unfavorable"

    # Determine urgency
    if nihss_score >= 15 or (request.lvo_suspected and nihss_score >= 6):
        urgency = "emergent"
    elif nihss_score >= 5:
        urgency = "urgent"
    else:
        urgency = "standard"

    recommendations = [
        "Stat CT head to rule out hemorrhage",
        "CBC, BMP, coagulation studies, troponin",
        "12-lead ECG for atrial fibrillation screening",
        "Continuous telemetry monitoring",
    ]
    if tpa_eligible and nihss_score >= 4:
        recommendations.insert(0, "Administer IV alteplase 0.9 mg/kg (max 90 mg)")
    if thrombectomy_eligible:
        recommendations.insert(0, "Activate neurointerventional team for thrombectomy")

    # LLM evidence enrichment
    evidence = []
    llm = _get_llm(req)
    if llm and request.clinical_notes:
        try:
            prompt = (
                f"Acute stroke triage:\n"
                f"NIHSS: {nihss_score}, Onset: {request.onset_time_hours}h\n"
                f"LVO suspected: {request.lvo_suspected}\n"
                f"ASPECTS: {request.ct_aspects_score}\n"
                f"Notes: {request.clinical_notes}\n\n"
                f"Provide additional clinical considerations."
            )
            llm_result = llm.generate(prompt)
            evidence.append({"source": "llm", "text": llm_result, "score": 0.8})
        except Exception as exc:
            logger.warning(f"Stroke triage LLM failed: {exc}")

    return StrokeTriageResponse(
        nihss_score=nihss_score,
        nihss_severity=nihss_severity,
        tpa_eligible=tpa_eligible,
        tpa_reasoning=tpa_reasoning,
        thrombectomy_eligible=thrombectomy_eligible,
        thrombectomy_reasoning=thrombectomy_reasoning,
        aspects_interpretation=aspects_interp,
        stroke_type_suggestion="ischemic" if request.lvo_suspected else "undetermined",
        urgency_level=urgency,
        recommendations=recommendations,
        guidelines_cited=[
            "AHA/ASA 2019 Acute Ischemic Stroke Guidelines",
            "AHA/ASA 2019 Endovascular Treatment Guidelines",
        ],
        evidence=evidence,
    )


@router.post("/dementia/evaluate", response_model=DementiaEvaluateResponse)
async def dementia_evaluate(request: DementiaEvaluateRequest, req: Request):
    """Comprehensive dementia evaluation with cognitive screening and differential."""
    _increment_metric(req, "workflow_requests_total")

    # MoCA calculation
    moca_score = request.moca_total
    moca_interp = None
    if request.moca_items:
        result = _calculate_moca(request.moca_items)
        moca_score = int(result["total_score"])
        moca_interp = result["interpretation"]
    elif moca_score is not None:
        if moca_score >= 26:
            moca_interp = "Normal cognition"
        elif moca_score >= 18:
            moca_interp = "Mild cognitive impairment (MCI)"
        elif moca_score >= 10:
            moca_interp = "Moderate cognitive impairment"
        else:
            moca_interp = "Severe cognitive impairment"

    # Cognitive domain analysis
    cognitive_domains = {
        "memory": "not_assessed",
        "executive_function": "not_assessed",
        "language": "not_assessed",
        "visuospatial": "not_assessed",
        "attention": "not_assessed",
        "behavior": "not_assessed",
    }
    for symptom in request.dominant_symptoms:
        s = symptom.lower()
        if "memory" in s:
            cognitive_domains["memory"] = "impaired"
        if "executive" in s or "planning" in s:
            cognitive_domains["executive_function"] = "impaired"
        if "language" in s or "aphasia" in s or "word" in s:
            cognitive_domains["language"] = "impaired"
        if "visual" in s or "spatial" in s:
            cognitive_domains["visuospatial"] = "impaired"
        if "attention" in s:
            cognitive_domains["attention"] = "impaired"
        if "behavior" in s or "personality" in s or "disinhibition" in s:
            cognitive_domains["behavior"] = "impaired"

    # Differential diagnosis
    differential = []
    impaired_domains = [d for d, v in cognitive_domains.items() if v == "impaired"]
    if "memory" in impaired_domains:
        differential.append({"diagnosis": "Alzheimer's disease", "likelihood": "high", "rationale": "Prominent memory impairment is hallmark of AD"})
    if "behavior" in impaired_domains or "executive_function" in impaired_domains:
        differential.append({"diagnosis": "Frontotemporal dementia", "likelihood": "moderate", "rationale": "Behavioral/executive changes suggest FTD"})
    if request.motor_features and "parkinson" in request.motor_features.lower():
        differential.append({"diagnosis": "Lewy body dementia", "likelihood": "moderate", "rationale": "Parkinsonism with cognitive decline suggests DLB"})
    if "visual" in " ".join(request.dominant_symptoms).lower():
        differential.append({"diagnosis": "Lewy body dementia", "likelihood": "moderate", "rationale": "Visual hallucinations suggest DLB"})
    differential.append({"diagnosis": "Vascular dementia", "likelihood": "consider", "rationale": "Always consider vascular contributions"})

    if not differential:
        differential.append({"diagnosis": "Neurodegenerative dementia NOS", "likelihood": "pending", "rationale": "Further workup needed"})

    # ATN staging
    atn_staging = None
    if request.biomarkers:
        amyloid = request.biomarkers.get("amyloid", "unknown")
        tau = request.biomarkers.get("tau", "unknown")
        neurodegeneration = request.biomarkers.get("neurodegeneration", "unknown")
        a = "+" if amyloid in ("positive", "+", "abnormal") else ("-" if amyloid in ("negative", "-", "normal") else "?")
        t = "+" if tau in ("positive", "+", "abnormal") else ("-" if tau in ("negative", "-", "normal") else "?")
        n = "+" if neurodegeneration in ("positive", "+", "abnormal") else ("-" if neurodegeneration in ("negative", "-", "normal") else "?")
        atn_staging = f"A{a}T{t}N{n}"

    workup = [
        "MRI brain with volumetric analysis",
        "Comprehensive metabolic panel, CBC, TSH, B12, folate",
        "RPR, HIV screening",
        "Neuropsychological testing",
    ]
    if "memory" in impaired_domains:
        workup.append("Consider amyloid PET or CSF AD biomarkers")
    if request.age and request.age < 65:
        workup.append("Young-onset dementia panel: autoimmune encephalitis antibodies, genetic testing")

    treatment = []
    if moca_score is not None and moca_score < 26:
        treatment.append("Cholinesterase inhibitor (donepezil 5-10mg, rivastigmine, or galantamine)")
    if moca_score is not None and moca_score < 18:
        treatment.append("Consider memantine addition")

    recommendations = []
    if moca_score is not None and moca_score < 26:
        recommendations.append("Assess driving safety and functional independence")
        recommendations.append("Connect with dementia support resources and caregiver education")
    recommendations.append("Follow-up cognitive assessment in 6-12 months")

    # LLM enrichment
    evidence = []
    llm = _get_llm(req)
    if llm and request.clinical_notes:
        try:
            prompt = (
                f"Dementia evaluation:\n"
                f"MoCA: {moca_score}, Age: {request.age}\n"
                f"Symptoms: {', '.join(request.dominant_symptoms)}\n"
                f"Motor: {request.motor_features}\n"
                f"Notes: {request.clinical_notes}\n\n"
                f"Provide differential diagnosis considerations."
            )
            llm_result = llm.generate(prompt)
            evidence.append({"source": "llm", "text": llm_result, "score": 0.8})
        except Exception as exc:
            logger.warning("Dementia evaluation LLM failed: %s", exc)

    return DementiaEvaluateResponse(
        moca_score=moca_score,
        moca_interpretation=moca_interp,
        cognitive_domain_analysis=cognitive_domains,
        differential_diagnosis=differential,
        atn_staging=atn_staging,
        recommended_workup=workup,
        treatment_options=treatment,
        recommendations=recommendations,
        guidelines_cited=[
            "AAN 2018 MCI Practice Guideline",
            "NIA-AA 2018 ATN Research Framework",
            "AAN 2001 Dementia Practice Parameter",
        ],
        evidence=evidence,
    )


@router.post("/epilepsy/classify", response_model=EpilepsyClassifyResponse)
async def epilepsy_classify(request: EpilepsyClassifyRequest, req: Request):
    """Epilepsy seizure classification and syndrome identification."""
    _increment_metric(req, "workflow_requests_total")

    desc = request.seizure_description.lower()

    # Seizure type classification
    focal_features = []
    generalized_features = []

    if any(w in desc for w in ["aura", "focal", "one side", "unilateral", "deja vu", "epigastric"]):
        focal_features.append("Focal onset features identified")
    if any(w in desc for w in ["awareness", "staring", "automatisms", "lip smacking"]):
        focal_features.append("Impaired awareness features")
    if any(w in desc for w in ["bilateral", "generalized", "both sides", "whole body"]):
        generalized_features.append("Bilateral/generalized involvement")
    if any(w in desc for w in ["tonic", "stiffening"]):
        generalized_features.append("Tonic component")
    if any(w in desc for w in ["clonic", "jerking", "shaking"]):
        generalized_features.append("Clonic component")
    if any(w in desc for w in ["myoclonic", "jerk"]):
        generalized_features.append("Myoclonic component")
    if any(w in desc for w in ["absence", "blank stare", "brief staring"]):
        generalized_features.append("Absence features")

    # Determine seizure type
    if focal_features and generalized_features:
        seizure_type = "focal_to_bilateral_tonic_clonic"
    elif focal_features:
        if "impaired awareness" in " ".join(focal_features).lower():
            seizure_type = "focal_impaired_awareness"
        else:
            seizure_type = "focal_aware"
    elif "absence" in " ".join(generalized_features).lower():
        seizure_type = "generalized_absence"
    elif "myoclonic" in " ".join(generalized_features).lower():
        seizure_type = "generalized_myoclonic"
    elif generalized_features:
        seizure_type = "generalized_tonic_clonic"
    else:
        seizure_type = "unknown_onset"

    # Syndrome identification
    syndrome = None
    if request.age_at_onset is not None:
        if request.age_at_onset < 1 and "spasm" in desc:
            syndrome = "west"
        elif 4 <= request.age_at_onset <= 10 and seizure_type == "generalized_absence":
            syndrome = "childhood_absence"
        elif 12 <= request.age_at_onset <= 18 and "myoclonic" in desc:
            syndrome = "juvenile_myoclonic"
        elif request.age_at_onset < 2 and "prolonged" in desc and "febrile" in desc:
            syndrome = "dravet"
    if request.eeg_findings and "temporal" in request.eeg_findings.lower():
        syndrome = syndrome or "temporal_lobe"

    # EEG-MRI concordance
    concordance = None
    if request.eeg_findings and request.mri_findings:
        concordance = "EEG and MRI findings available for concordance analysis"
        if "temporal" in request.eeg_findings.lower() and "temporal" in request.mri_findings.lower():
            concordance = "Concordant: Both EEG and MRI localize to temporal lobe"
        elif "normal" in request.mri_findings.lower():
            concordance = "MRI-negative epilepsy -- EEG findings guide localization"

    # Surgical candidacy
    surgical_candidacy = None
    if len(request.current_aeds) >= 2 and request.seizure_frequency:
        surgical_candidacy = "Drug-resistant epilepsy criteria may be met -- consider epilepsy surgery evaluation"

    # AED recommendations
    aed_recs = []
    if seizure_type.startswith("focal"):
        aed_recs.extend(["Levetiracetam (first-line focal)", "Lamotrigine", "Oxcarbazepine"])
    elif seizure_type.startswith("generalized"):
        aed_recs.extend(["Levetiracetam (first-line generalized)", "Valproate (if not childbearing potential)", "Lamotrigine"])
    if syndrome == "juvenile_myoclonic":
        aed_recs = ["Valproate or levetiracetam (first-line JME)", "Lamotrigine (may worsen myoclonus)"]
    if syndrome == "dravet":
        aed_recs = ["Stiripentol + valproate + clobazam", "Cannabidiol (Epidiolex)", "Fenfluramine"]

    recommendations = [
        "Routine and/or prolonged EEG monitoring",
        "MRI brain with epilepsy protocol",
    ]
    if len(request.current_aeds) >= 2:
        recommendations.append("Refer to comprehensive epilepsy center")
    if surgical_candidacy:
        recommendations.append(surgical_candidacy)

    # LLM enrichment
    evidence = []
    llm = _get_llm(req)
    if llm and request.clinical_notes:
        try:
            prompt = (
                f"Epilepsy classification:\n"
                f"Description: {request.seizure_description}\n"
                f"Age at onset: {request.age_at_onset}\n"
                f"EEG: {request.eeg_findings}\n"
                f"MRI: {request.mri_findings}\n"
                f"Notes: {request.clinical_notes}\n\n"
                f"Provide classification and management considerations."
            )
            llm_result = llm.generate(prompt)
            evidence.append({"source": "llm", "text": llm_result, "score": 0.8})
        except Exception as exc:
            logger.warning("Epilepsy classification LLM failed: %s", exc)

    return EpilepsyClassifyResponse(
        seizure_type=seizure_type,
        epilepsy_syndrome=syndrome,
        focal_features=focal_features,
        generalized_features=generalized_features,
        eeg_mri_concordance=concordance,
        surgical_candidacy=surgical_candidacy,
        aed_recommendations=aed_recs,
        recommendations=recommendations,
        guidelines_cited=[
            "ILAE 2017 Seizure Classification",
            "ILAE 2017 Epilepsy Classification",
            "AAN/AES 2018 Drug-Resistant Epilepsy Guideline",
        ],
        evidence=evidence,
    )


@router.post("/tumor/grade", response_model=TumorGradeResponse)
async def tumor_grade(request: TumorGradeRequest, req: Request):
    """Brain tumor WHO 2021 grading with molecular integration."""
    _increment_metric(req, "workflow_requests_total")

    markers = [m.lower() for m in request.molecular_markers]
    histology_lower = request.histology.lower()

    # Determine WHO grade and integrated diagnosis
    who_grade = "grade_2"
    integrated_dx = request.histology
    molecular_profile = {m: "detected" for m in request.molecular_markers}

    if "glioblastoma" in histology_lower or "gbm" in histology_lower:
        who_grade = "grade_4"
        if "idh_mutant" in markers:
            integrated_dx = "Astrocytoma, IDH-mutant, WHO grade 4"
        else:
            integrated_dx = "Glioblastoma, IDH-wildtype, WHO grade 4"
    elif "astrocytoma" in histology_lower:
        if "idh_mutant" in markers:
            who_grade = "grade_2" if "atrx_loss" not in markers else "grade_3"
            integrated_dx = f"Astrocytoma, IDH-mutant, WHO {who_grade.replace('_', ' ')}"
        else:
            who_grade = "grade_4"
            integrated_dx = "Glioblastoma, IDH-wildtype, WHO grade 4"
    elif "oligodendroglioma" in histology_lower:
        if "1p19q_codeleted" in markers and "idh_mutant" in markers:
            who_grade = "grade_2"
            integrated_dx = "Oligodendroglioma, IDH-mutant, 1p/19q-codeleted, WHO grade 2"
    elif "meningioma" in histology_lower:
        who_grade = "grade_1"
        integrated_dx = "Meningioma, WHO grade 1"
    elif "medulloblastoma" in histology_lower:
        who_grade = "grade_4"
        integrated_dx = "Medulloblastoma, WHO grade 4"

    # Treatment recommendations
    treatment = []
    if who_grade == "grade_4":
        treatment.extend([
            "Maximum safe surgical resection",
            "Concurrent temozolomide + radiation (Stupp protocol)",
        ])
        if "mgmt_methylated" in markers:
            treatment.append("MGMT methylated -- favorable response to temozolomide expected")
        else:
            treatment.append("MGMT unmethylated -- consider clinical trial alternatives")
        treatment.append("Adjuvant temozolomide x6 cycles")
        treatment.append("Consider tumor treating fields (TTFields/Optune)")
    elif who_grade in ("grade_2", "grade_3"):
        treatment.extend([
            "Maximum safe surgical resection",
            "Post-operative MRI surveillance",
        ])
        if who_grade == "grade_3" or (request.age and request.age >= 40):
            treatment.append("Radiation + PCV or temozolomide")

    # Prognosis
    if who_grade == "grade_4" and "idh_wildtype" in markers:
        prognosis = "Median survival 14-16 months with standard therapy"
    elif who_grade == "grade_4" and "idh_mutant" in markers:
        prognosis = "Median survival 3-5 years (better than IDH-wildtype GBM)"
    elif "1p19q_codeleted" in markers:
        prognosis = "Favorable prognosis with median survival >10 years"
    elif who_grade in ("grade_1", "grade_2"):
        prognosis = "Generally favorable prognosis with appropriate management"
    else:
        prognosis = "Prognosis depends on molecular profile, extent of resection, and performance status"

    trial_relevance = []
    if who_grade == "grade_4":
        trial_relevance.append("Eligible for glioblastoma immunotherapy trials")
    if "idh_mutant" in markers:
        trial_relevance.append("IDH-mutant targeted therapy trials (vorasidenib, ivosidenib)")

    recommendations = [
        "Multidisciplinary tumor board review",
        "Serial MRI surveillance per NCCN guidelines",
        "Neuropsychological assessment",
    ]
    if request.kps is not None and request.kps < 70:
        recommendations.append("Supportive/palliative care focus given low KPS")

    return TumorGradeResponse(
        who_grade=who_grade,
        integrated_diagnosis=integrated_dx,
        molecular_profile=molecular_profile,
        treatment_recommendations=treatment,
        prognosis_summary=prognosis,
        nccn_category="NCCN CNS Cancers",
        clinical_trial_relevance=trial_relevance,
        recommendations=recommendations,
        guidelines_cited=[
            "WHO 2021 CNS Tumor Classification",
            "NCCN CNS Cancers v2.2025",
            "Stupp et al. NEJM 2005 (temozolomide + RT)",
        ],
    )


@router.post("/ms/assess", response_model=MSAssessResponse)
async def ms_assess(request: MSAssessRequest, req: Request):
    """MS disease monitoring with EDSS, NEDA, and DMT assessment."""
    _increment_metric(req, "workflow_requests_total")

    # EDSS
    edss_score = request.edss_score
    edss_interp = None
    if request.edss_items:
        result = _calculate_edss(request.edss_items)
        edss_score = result["total_score"]
        edss_interp = result["interpretation"]
    elif edss_score is not None:
        if edss_score <= 1.5:
            edss_interp = "Minimal disability"
        elif edss_score <= 3.5:
            edss_interp = "Moderate disability"
        elif edss_score <= 5.5:
            edss_interp = "Significant disability"
        else:
            edss_interp = "Severe disability"

    # Disease activity
    relapses = request.relapse_count_2yr or 0
    new_lesions = (request.new_t2_lesions or 0) + (request.gad_enhancing_lesions or 0)
    if relapses >= 2 or new_lesions >= 3:
        disease_activity = "highly_active"
    elif relapses >= 1 or new_lesions >= 1:
        disease_activity = "active"
    else:
        disease_activity = "stable"

    # NEDA status
    neda_3 = relapses == 0 and new_lesions == 0
    edss_progression = False
    if edss_score is not None and edss_score > 3.0:
        neda_3 = False
        edss_progression = True
    neda_status = {
        "neda_3_achieved": neda_3,
        "no_relapses": relapses == 0,
        "no_new_lesions": new_lesions == 0,
        "no_edss_progression": not edss_progression,
    }

    # DMT assessment
    dmt_assessment = {
        "current_dmt": request.current_dmt or "none",
        "treatment_response": "adequate" if neda_3 else "suboptimal",
    }

    escalation = None
    if disease_activity in ("active", "highly_active") and request.current_dmt:
        escalation = "Consider DMT escalation to higher-efficacy therapy"
    elif disease_activity == "highly_active" and not request.current_dmt:
        escalation = "Initiate high-efficacy DMT (natalizumab, ocrelizumab, or ofatumumab)"

    # JCV risk for natalizumab consideration
    if request.jcv_status == "positive":
        dmt_assessment["jcv_warning"] = "JCV antibody positive -- avoid natalizumab if index >1.5"

    mri_burden = None
    if new_lesions > 0:
        mri_burden = f"{request.new_t2_lesions or 0} new T2, {request.gad_enhancing_lesions or 0} Gd+ lesions"

    recommendations = [
        "Serial MRI brain and spine per McDonald 2017 protocol",
        "Annual neurological examination with EDSS",
        "Monitor for treatment-specific adverse effects",
    ]
    if disease_activity != "stable":
        recommendations.insert(0, "Discuss DMT optimization with patient")
    if request.ms_course == "spms":
        recommendations.append("Consider siponimod if active SPMS")

    return MSAssessResponse(
        edss_score=edss_score,
        edss_interpretation=edss_interp,
        disease_activity=disease_activity,
        neda_status=neda_status,
        dmt_assessment=dmt_assessment,
        escalation_recommendation=escalation,
        mri_disease_burden=mri_burden,
        recommendations=recommendations,
        guidelines_cited=[
            "McDonald 2017 Diagnostic Criteria",
            "AAN 2018 MS DMT Practice Guideline",
            "ECTRIMS/EAN 2018 MS Treatment Guidelines",
        ],
    )


@router.post("/parkinsons/assess", response_model=ParkinsonsAssessResponse)
async def parkinsons_assess(request: ParkinsonsAssessRequest, req: Request):
    """Parkinson's disease assessment with motor scoring and medication guidance."""
    _increment_metric(req, "workflow_requests_total")

    # UPDRS
    updrs_score = request.updrs_total
    updrs_interp = None
    if request.updrs_items:
        result = _calculate_updrs(request.updrs_items)
        updrs_score = int(result["total_score"])
        updrs_interp = result["interpretation"]
    elif updrs_score is not None:
        if updrs_score <= 10:
            updrs_interp = "Minimal motor impairment"
        elif updrs_score <= 32:
            updrs_interp = "Mild motor impairment"
        elif updrs_score <= 58:
            updrs_interp = "Moderate motor impairment"
        else:
            updrs_interp = "Severe motor impairment"

    # Hoehn-Yahr
    hy_stage = request.hoehn_yahr
    hy_desc = None
    if hy_stage is not None:
        hy_result = _calculate_hoehn_yahr({"stage": hy_stage})
        hy_desc = hy_result["interpretation"]

    # Motor subtype
    motor_subtype = "indeterminate"
    if updrs_score is not None:
        if request.updrs_items:
            tremor_keys = [k for k in request.updrs_items if "tremor" in k.lower() or "3.15" in k or "3.16" in k or "3.17" in k]
            pigd_keys = [k for k in request.updrs_items if "gait" in k.lower() or "postural" in k.lower() or "3.10" in k or "3.11" in k or "3.12" in k]
            tremor_sum = sum(int(request.updrs_items.get(k, 0)) for k in tremor_keys)
            pigd_sum = sum(int(request.updrs_items.get(k, 0)) for k in pigd_keys)
            if tremor_sum > pigd_sum:
                motor_subtype = "tremor_dominant"
            elif pigd_sum > tremor_sum:
                motor_subtype = "postural_instability_gait_difficulty"

    # Medication assessment
    med_assessment = {
        "current_medications": request.current_medications,
        "motor_fluctuations": request.motor_fluctuations,
        "dyskinesia": request.dyskinesia,
    }
    if request.motor_fluctuations:
        med_assessment["adjustment_needed"] = "Consider levodopa dose fractionation or COMT/MAO-B inhibitor addition"
    if request.dyskinesia:
        med_assessment["dyskinesia_management"] = "Consider amantadine or levodopa dose reduction"

    # DBS candidacy
    dbs_candidacy = None
    if request.motor_fluctuations and request.disease_duration_years and request.disease_duration_years >= 4:
        dbs_candidacy = "Potential DBS candidate -- motor fluctuations with adequate levodopa response required"
    if request.age and request.age > 75:
        dbs_candidacy = "Advanced age -- DBS benefit may be limited; careful risk-benefit assessment needed"

    # Non-motor management
    non_motor = []
    for symptom in request.non_motor_symptoms:
        s = symptom.lower()
        if "rbd" in s or "rem" in s:
            non_motor.append("RBD: Melatonin 3-12mg, clonazepam 0.25-0.5mg")
        if "depression" in s:
            non_motor.append("Depression: SNRI or SSRI; pramipexole has antidepressant effect")
        if "constipation" in s:
            non_motor.append("Constipation: Fiber, PEG 3350, exercise")
        if "orthostatic" in s:
            non_motor.append("Orthostatic hypotension: Fludrocortisone, midodrine, droxidopa")
        if "dementia" in s or "cognitive" in s:
            non_motor.append("Cognitive: Rivastigmine (only cholinesterase inhibitor FDA-approved for PDD)")

    recommendations = [
        "Regular movement disorders follow-up every 3-6 months",
        "Physical therapy and exercise program",
        "Speech therapy referral if dysarthria present",
    ]
    if dbs_candidacy and "Potential" in dbs_candidacy:
        recommendations.append("Neurosurgical consultation for DBS evaluation")

    return ParkinsonsAssessResponse(
        updrs_score=updrs_score,
        updrs_interpretation=updrs_interp,
        hoehn_yahr_stage=hy_stage,
        hoehn_yahr_description=hy_desc,
        motor_subtype=motor_subtype,
        medication_assessment=med_assessment,
        dbs_candidacy=dbs_candidacy,
        non_motor_management=non_motor,
        recommendations=recommendations,
        guidelines_cited=[
            "MDS Clinical Diagnostic Criteria for PD (2015)",
            "AAN Quality Measure: PD Annual Assessment",
            "MDS Evidence-Based Review of PD Treatments",
        ],
    )


@router.post("/headache/classify", response_model=HeadacheClassifyResponse)
async def headache_classify(request: HeadacheClassifyRequest, req: Request):
    """Headache classification per ICHD-3 criteria with treatment guidance."""
    _increment_metric(req, "workflow_requests_total")

    request.headache_description.lower()
    symptoms = [s.lower() for s in request.associated_symptoms]

    # Classify headache type
    ichd3_criteria = []
    headache_type = "unclassified"

    # Migraine criteria
    migraine_features = 0
    if request.location and "unilateral" in request.location.lower():
        migraine_features += 1
        ichd3_criteria.append("Unilateral location")
    if request.quality and "throb" in request.quality.lower():
        migraine_features += 1
        ichd3_criteria.append("Pulsating/throbbing quality")
    if request.duration_hours and 4 <= request.duration_hours <= 72:
        migraine_features += 1
        ichd3_criteria.append("Duration 4-72 hours")
    if any(s in symptoms for s in ["nausea", "vomiting"]):
        migraine_features += 1
        ichd3_criteria.append("Nausea/vomiting")
    if any(s in symptoms for s in ["photophobia", "phonophobia"]):
        migraine_features += 1
        ichd3_criteria.append("Photophobia/phonophobia")

    if migraine_features >= 3:
        if request.aura:
            headache_type = "migraine_with_aura"
        elif request.frequency_per_month and request.frequency_per_month >= 15:
            headache_type = "chronic_migraine"
        else:
            headache_type = "migraine_without_aura"
    elif request.location and "bilateral" in request.location.lower() and request.quality and "pressure" in request.quality.lower():
        headache_type = "tension_type_episodic"
        ichd3_criteria.append("Bilateral pressing/tightening quality")
        if request.frequency_per_month and request.frequency_per_month >= 15:
            headache_type = "tension_type_chronic"
    elif any(s in symptoms for s in ["tearing", "rhinorrhea", "ptosis", "miosis"]) and request.location and "unilateral" in request.location.lower():
        headache_type = "cluster"
        ichd3_criteria.append("Unilateral with autonomic features")

    # HIT-6
    hit6_score = None
    hit6_interp = None
    if request.hit6_items:
        result = _calculate_hit6(request.hit6_items)
        hit6_score = int(result["total_score"])
        hit6_interp = result["interpretation"]

    # Red flag assessment
    red_flag_assessment = {
        "flags_present": request.red_flags,
        "urgent_workup_needed": len(request.red_flags) > 0,
    }
    if request.red_flags:
        red_flag_assessment["action"] = "Urgent neuroimaging and further evaluation required"

    # Medication overuse risk
    moh_risk = None
    if request.analgesic_days_per_month is not None:
        if request.analgesic_days_per_month >= 15:
            moh_risk = "High risk -- meets MOH criteria (>=15 days/month analgesic use)"
        elif request.analgesic_days_per_month >= 10:
            moh_risk = "Moderate risk -- approaching MOH threshold"

    # Acute treatment
    acute = []
    if headache_type.startswith("migraine"):
        acute = ["Triptans (sumatriptan 50-100mg PO, or SC 6mg)", "NSAIDs (ibuprofen 400-600mg, naproxen 500mg)", "Antiemetics (metoclopramide 10mg) if nausea"]
    elif headache_type.startswith("tension"):
        acute = ["NSAIDs (ibuprofen 400mg)", "Acetaminophen 1000mg"]
    elif headache_type == "cluster":
        acute = ["Sumatriptan 6mg SC", "High-flow oxygen 12-15 L/min x 15 min"]

    # Preventive treatment
    preventive = []
    if headache_type.startswith("migraine") and request.frequency_per_month and request.frequency_per_month >= 4:
        preventive = [
            "CGRP monoclonal antibody (erenumab, fremanezumab, galcanezumab)",
            "Topiramate 50-100mg BID",
            "Propranolol 80-240mg/day",
            "Amitriptyline 25-75mg QHS",
        ]
    elif headache_type == "cluster":
        preventive = ["Verapamil 240-960mg/day", "Galcanezumab (FDA-approved for episodic cluster)"]
    elif headache_type == "chronic_migraine":
        preventive = [
            "OnabotulinumtoxinA (Botox) 155-195 units",
            "CGRP monoclonal antibody",
            "Topiramate 100-200mg/day",
        ]

    recommendations = []
    if request.red_flags:
        recommendations.append("Urgent evaluation: CT/MRI, LP if indicated")
    if moh_risk and "High" in moh_risk:
        recommendations.append("Initiate medication withdrawal protocol")
    recommendations.append("Maintain headache diary")
    if request.frequency_per_month and request.frequency_per_month >= 4:
        recommendations.append("Preventive therapy indicated")

    return HeadacheClassifyResponse(
        headache_type=headache_type,
        ichd3_criteria_met=ichd3_criteria,
        hit6_score=hit6_score,
        hit6_interpretation=hit6_interp,
        red_flag_assessment=red_flag_assessment,
        medication_overuse_risk=moh_risk,
        acute_treatment=acute,
        preventive_treatment=preventive,
        recommendations=recommendations,
        guidelines_cited=[
            "ICHD-3 (International Classification of Headache Disorders 3rd edition)",
            "AAN/AHS 2021 Migraine Prevention Guideline",
            "AHS 2019 Acute Migraine Treatment Consensus Statement",
        ],
    )


@router.post("/neuromuscular/evaluate", response_model=NeuromuscularEvaluateResponse)
async def neuromuscular_evaluate(request: NeuromuscularEvaluateRequest, req: Request):
    """Neuromuscular evaluation with localization and differential diagnosis."""
    _increment_metric(req, "workflow_requests_total")

    desc = request.presentation.lower()

    # Localization
    localization = "undetermined"
    pattern = "mixed"
    if request.emg_findings:
        emg = request.emg_findings.lower()
        if "demyelinating" in emg:
            localization = "peripheral_nerve"
            pattern = "demyelinating"
        elif "axonal" in emg:
            localization = "peripheral_nerve"
            pattern = "axonal_sensorimotor" if request.sensory_involvement else "axonal_motor"
        elif "myopathic" in emg or "short duration" in emg:
            localization = "muscle"
            pattern = "myopathic"
        elif "decrement" in emg or "jitter" in emg:
            localization = "neuromuscular_junction"
            pattern = "nmj_postsynaptic"
        elif "increment" in emg:
            localization = "neuromuscular_junction"
            pattern = "nmj_presynaptic"
    elif request.weakness_pattern:
        wp = request.weakness_pattern.lower()
        if "proximal" in wp:
            localization = "muscle_or_nmj"
        elif "distal" in wp:
            localization = "peripheral_nerve"
        elif "bulbar" in wp:
            localization = "bulbar"

    # Differential diagnosis
    differential = []
    if "achr" in " ".join(request.antibodies).lower() or "musk" in " ".join(request.antibodies).lower():
        differential.append({"diagnosis": "Myasthenia gravis", "likelihood": "high", "rationale": "Positive NMJ antibodies"})
    if localization == "peripheral_nerve" and pattern == "demyelinating":
        differential.append({"diagnosis": "CIDP", "likelihood": "high", "rationale": "Demyelinating neuropathy pattern"})
        differential.append({"diagnosis": "GBS (if acute)", "likelihood": "consider", "rationale": "Acute demyelinating polyradiculoneuropathy"})
    if localization == "muscle" or (request.ck_level and request.ck_level > 1000):
        differential.append({"diagnosis": "Inflammatory myopathy", "likelihood": "moderate", "rationale": "Myopathic pattern with elevated CK"})
    if any(w in desc for w in ["fasciculation", "wasting", "upper motor", "lower motor"]):
        differential.append({"diagnosis": "ALS", "likelihood": "consider", "rationale": "Combined UMN/LMN features"})
    if not differential:
        differential.append({"diagnosis": "Neuromuscular disorder NOS", "likelihood": "pending", "rationale": "Further workup needed"})

    # ALSFRS-R
    alsfrs_score = None
    alsfrs_interp = None
    if request.alsfrs_items:
        result = _calculate_alsfrs(request.alsfrs_items)
        alsfrs_score = int(result["total_score"])
        alsfrs_interp = result["interpretation"]

    # Workup
    workup = [
        "EMG/NCS if not yet performed",
        "CK, aldolase levels",
        "Comprehensive metabolic panel",
    ]
    if localization in ("neuromuscular_junction", "muscle_or_nmj"):
        workup.extend(["AChR antibodies", "MuSK antibodies", "CT chest (thymoma screening)"])
    if localization == "peripheral_nerve":
        workup.extend(["B12, folate, TSH", "SPEP/UPEP with immunofixation", "HbA1c"])
    if pattern == "demyelinating":
        workup.append("Lumbar puncture with CSF protein (albuminocytologic dissociation)")
    if request.ck_level and request.ck_level > 1000:
        workup.append("Muscle biopsy")
    if any(w in desc for w in ["family", "hereditary", "genetic"]) or request.genetic_results:
        workup.append("Genetic testing panel for hereditary neuropathy/myopathy")

    # Treatment
    treatment = []
    for dx in differential:
        if dx["diagnosis"] == "Myasthenia gravis":
            treatment.extend(["Pyridostigmine 60mg TID", "Prednisone or steroid-sparing immunotherapy"])
        elif dx["diagnosis"] == "CIDP":
            treatment.extend(["IVIg 2g/kg loading then 1g/kg monthly", "Plasma exchange alternative"])
        elif dx["diagnosis"] == "ALS":
            treatment.extend(["Riluzole 50mg BID", "Edaravone or tofersen (if SOD1+)", "Multidisciplinary ALS clinic"])

    recommendations = [
        "Neuromuscular specialist referral",
        "Pulmonary function testing (FVC)",
        "Swallow evaluation if bulbar symptoms",
    ]
    if request.progression_rate and "rapid" in request.progression_rate.lower():
        recommendations.insert(0, "Urgent evaluation given rapid progression")

    return NeuromuscularEvaluateResponse(
        localization=localization,
        pattern_classification=pattern,
        differential_diagnosis=differential,
        alsfrs_score=alsfrs_score,
        alsfrs_interpretation=alsfrs_interp,
        recommended_workup=workup,
        treatment_options=treatment,
        recommendations=recommendations,
        guidelines_cited=[
            "AAN 2009 Evaluation of Distal Symmetric Polyneuropathy",
            "AAN 2021 ALS Practice Guideline",
            "EFNS/PNS 2021 CIDP Guideline",
        ],
    )


@router.post("/workflow/{workflow_type}", response_model=WorkflowResponse)
async def generic_workflow(workflow_type: str, request: WorkflowRequest, req: Request):
    """Generic workflow dispatch for any supported workflow type."""
    _increment_metric(req, "workflow_requests_total")

    valid_types = [wt.value for wt in NeuroWorkflowType]
    if workflow_type not in valid_types:
        raise HTTPException(
            status_code=400,
            detail=f"Unknown workflow_type '{workflow_type}'. Valid types: {valid_types}",
        )

    engine = _get_workflow_engine(req)
    if engine is None:
        return WorkflowResponse(
            workflow_type=workflow_type,
            status="completed",
            result="Workflow engine unavailable. Service in degraded mode.",
            note="Workflow engine not initialized.",
        )

    data = request.data
    if request.question:
        data["question"] = request.question

    try:
        result = engine.execute(workflow_type, data)
        return WorkflowResponse(
            workflow_type=result.get("workflow_type", workflow_type),
            status=result.get("status", "completed"),
            result=result.get("result", ""),
            evidence_used=result.get("evidence_used", False),
            note=result.get("note"),
        )
    except Exception as exc:
        logger.error(f"Workflow execution failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


# =====================================================================
# Reference Endpoints
# =====================================================================

@router.get("/domains")
async def list_domains():
    """Reference catalogue of neurology clinical domains."""
    return {
        "domains": [
            {"id": "stroke", "name": "Cerebrovascular / Stroke", "description": "Acute stroke triage, secondary prevention, neurovascular"},
            {"id": "dementia", "name": "Dementia / Cognitive", "description": "Alzheimer's, FTD, LBD, vascular, cognitive screening"},
            {"id": "epilepsy", "name": "Epilepsy / Seizures", "description": "Seizure classification, syndromes, AED management, surgical evaluation"},
            {"id": "tumors", "name": "Neuro-Oncology", "description": "Brain tumor grading, molecular classification, treatment planning"},
            {"id": "ms", "name": "Multiple Sclerosis", "description": "MS diagnosis, DMT management, NEDA monitoring, relapse management"},
            {"id": "parkinsons", "name": "Movement Disorders / Parkinson's", "description": "PD assessment, medication optimization, DBS evaluation"},
            {"id": "headache", "name": "Headache / Migraine", "description": "ICHD-3 classification, preventive therapy, acute management"},
            {"id": "neuromuscular", "name": "Neuromuscular", "description": "Neuropathy, NMJ disorders, myopathy, ALS, EMG interpretation"},
        ]
    }


@router.get("/scales")
async def list_scales():
    """Reference catalogue of validated neurological clinical scales."""
    return {
        "scales": [
            {"id": "nihss", "name": "NIH Stroke Scale", "max_score": 42, "domain": "stroke", "description": "Stroke severity assessment (15 items)"},
            {"id": "gcs", "name": "Glasgow Coma Scale", "max_score": 15, "domain": "consciousness", "description": "Level of consciousness (eye, verbal, motor)"},
            {"id": "moca", "name": "Montreal Cognitive Assessment", "max_score": 30, "domain": "cognition", "description": "Cognitive screening (8 domains)"},
            {"id": "updrs", "name": "MDS-UPDRS Part III", "max_score": 132, "domain": "parkinsons", "description": "Motor examination (33 items)"},
            {"id": "edss", "name": "Expanded Disability Status Scale", "max_score": 10.0, "domain": "ms", "description": "MS disability (8 functional systems)"},
            {"id": "mrs", "name": "Modified Rankin Scale", "max_score": 6, "domain": "stroke_outcome", "description": "Global disability/dependence (7 levels)"},
            {"id": "hit6", "name": "Headache Impact Test-6", "max_score": 78, "domain": "headache", "description": "Headache disability (6 items)"},
            {"id": "alsfrs", "name": "ALS Functional Rating Scale - Revised", "max_score": 48, "domain": "als", "description": "ALS functional status (12 items)"},
            {"id": "aspects", "name": "ASPECTS", "max_score": 10, "domain": "stroke_imaging", "description": "CT early ischemic changes (10 MCA regions)"},
            {"id": "hoehn_yahr", "name": "Hoehn and Yahr Scale", "max_score": 5.0, "domain": "parkinsons", "description": "PD staging (0-5)"},
        ]
    }


@router.get("/guidelines")
async def list_guidelines():
    """Reference catalogue of neurology practice guidelines."""
    return {
        "guidelines": [
            {"id": "aha_ais_2019", "name": "AHA/ASA Acute Ischemic Stroke 2019", "description": "Early management of acute ischemic stroke"},
            {"id": "aha_evt_2019", "name": "AHA/ASA Endovascular Treatment 2019", "description": "Mechanical thrombectomy guidelines"},
            {"id": "aan_mci_2018", "name": "AAN MCI Practice Guideline 2018", "description": "Mild cognitive impairment evaluation and management"},
            {"id": "nia_aa_atn_2018", "name": "NIA-AA ATN Framework 2018", "description": "Alzheimer's disease biomarker framework"},
            {"id": "ilae_2017", "name": "ILAE 2017 Classification", "description": "Seizure and epilepsy classification"},
            {"id": "who_cns_2021", "name": "WHO CNS Tumors 2021", "description": "5th edition CNS tumor classification"},
            {"id": "nccn_cns", "name": "NCCN CNS Cancers", "description": "CNS cancer treatment guidelines"},
            {"id": "mcdonald_2017", "name": "McDonald Criteria 2017", "description": "MS diagnostic criteria"},
            {"id": "aan_ms_dmt_2018", "name": "AAN MS DMT 2018", "description": "Disease-modifying therapy guidelines"},
            {"id": "mds_pd_2015", "name": "MDS PD Criteria 2015", "description": "Parkinson's disease clinical diagnostic criteria"},
            {"id": "ichd_3", "name": "ICHD-3", "description": "International Classification of Headache Disorders 3rd edition"},
            {"id": "aan_headache_2021", "name": "AAN/AHS Migraine Prevention 2021", "description": "Migraine preventive treatment guideline"},
        ]
    }


@router.get("/knowledge-version")
async def knowledge_version():
    """Version metadata for the neurology knowledge base."""
    return {
        "agent": "neurology-intelligence-agent",
        "neurology_domains": len(NeuroWorkflowType),
        "clinical_scales": len(ClinicalScaleType),
        **KNOWLEDGE_VERSION,
    }
