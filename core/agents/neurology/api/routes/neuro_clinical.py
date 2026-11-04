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

from src.clinical_scales import (
    SUPPORTED_SCALES,
    ScaleDataError,
    UPDRSCalculator,
    calculate_scale,
)
from src.models import (
    NeuroWorkflowType,
    ClinicalScaleType,
    ScaleResult,
)
from src.knowledge import KNOWLEDGE_VERSION

router = APIRouter(prefix="/v1/neuro", tags=["neurology"])




def _as_row(r) -> dict:
    """Normalise a search hit to a plain dict.

    The RAG engines return dataclass results (NeuroSearchResult / RareDiseaseSearchResult /
    ...), but these routes were written against dicts. Calling .get() on a dataclass raises
    AttributeError, which the surrounding try/except swallowed -- so retrieval looked empty
    from a service that was otherwise healthy. One helper rather than duck-typing per site.
    """
    if isinstance(r, dict):
        out = dict(r)
    else:
        out = {
            "collection": getattr(r, "collection", "unknown"),
            "record_id": getattr(r, "record_id", ""),
            "score": getattr(r, "score", 0.0),
            "text": getattr(r, "text", ""),
            "metadata": getattr(r, "metadata", {}) or {},
            "relevance": getattr(r, "relevance", ""),
        }
    out.setdefault("metadata", {})
    out.setdefault("text", out.get("content", ""))
    out.setdefault("score", 0.0)
    return out

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
    # False when retrieval returned nothing: the prose is then the model's own general
    # knowledge, not an answer grounded in this knowledge base. Confidence is 0.0 there.
    grounded: bool = True
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
    # A total computed from a partial instrument is not the patient's score. These say so
    # on the wire, so a caller is not left to infer completeness from `items_received`.
    items_scored: Optional[int] = None
    items_expected: Optional[int] = None
    data_complete: bool = True
    data_warnings: List[str] = []


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
    # `tpa_eligible: bool` cannot distinguish "excluded" from "nobody measured the BP", and
    # the endpoint used to answer `true` for both. These carry the third state explicitly:
    # "eligible" | "not_eligible" | "indeterminate". `tpa_eligible` is true only for
    # "eligible".
    tpa_determination: str = "indeterminate"
    thrombectomy_determination: str = "indeterminate"
    missing_inputs: List[str] = []
    nihss_complete: bool = True
    # The model's own prose. It is NOT retrieved evidence and must not sit under a key the
    # honesty gate skips, so it lands here, in a key the gate scans and annotates.
    discussion: Optional[str] = None


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
    moca_complete: bool = True
    data_warnings: List[str] = []
    # Model prose, moved out of `evidence` -- see the note in stroke_triage (G3).
    discussion: Optional[str] = None


# -- Epilepsy Classification --

class EpilepsyClassifyRequest(BaseModel):
    """Epilepsy classification request."""
    seizure_description: str = Field(..., min_length=5, description="Seizure semiology description")
    age_at_onset: Optional[int] = Field(None, ge=0, le=120)
    eeg_findings: Optional[str] = Field(None, description="EEG description or conclusion")
    mri_findings: Optional[str] = Field(None, description="Brain MRI findings")
    seizure_frequency: Optional[str] = Field(None, description="Frequency (e.g., 2/month)")
    current_aeds: List[str] = Field(default_factory=list, description="Current anti-epileptic drugs")
    failed_aeds: List[str] = Field(
        default_factory=list,
        description=(
            "ASMs that failed an adequate trial. Drug-resistant epilepsy requires FAILURE of "
            "two adequately chosen, tolerated and dosed trials (ILAE), not two concurrent "
            "prescriptions."
        ),
    )
    clinical_diagnosis: Optional[str] = Field(
        None,
        description=(
            "Clinically established epilepsy syndrome, if one is already diagnosed "
            "(e.g. 'Dravet syndrome'). Used for syndrome-specific drug contraindications "
            "that must not wait on genetic testing."
        ),
    )
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
    contraindication_warnings: List[str] = []
    recommendations: List[str]
    guidelines_cited: List[str]
    evidence: List[dict] = []
    # Model prose, moved out of `evidence` -- see the note in stroke_triage (G3).
    discussion: Optional[str] = None


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
    edss_prior: Optional[float] = Field(
        None, ge=0.0, le=10.0,
        description=(
            "EDSS at the previous assessment. EDSS progression is a CHANGE, so it cannot be "
            "determined without it."
        ),
    )
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
    data_warnings: List[str] = []
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

def _scale_result_to_dict(result: ScaleResult) -> dict:
    """Map a ScaleResult onto the dict shape these routes return."""
    return {
        "total_score": result.score,
        "max_score": result.max_score,
        "interpretation": result.interpretation,
        "severity_category": result.severity_category,
        "thresholds": dict(result.thresholds),
        "recommendations": list(result.recommendations),
        "items_scored": result.items_scored,
        "items_expected": result.items_expected,
        "complete": result.complete,
        "warnings": list(result.warnings),
    }


def _calculate(scale: str, items: dict) -> dict:
    """Score *scale* through the single shared implementation.

    STRUCTURAL FIX. Each `_calculate_<scale>` below used to be a second, independent
    implementation written inline in this file, and the Streamlit UI calls this file -- so the
    route-local versions were what a clinician actually saw, and they disagreed with
    `src/clinical_scales.py` on bandings, on defaults and on eligibility. That is why the
    landed tPA fix did not take: it went into the other copy. Each of these is now a thin
    adapter over `src.clinical_scales.calculate_scale`, so there is one answer per scale.

    The defects that lived only in the route copies, now gone with them:

    * MoCA summed *every* numeric value in `items`, and the UI puts `education_years`
      (default 16) in that dict -- a true 12/30 exam reported 28/30 "Normal cognition" with
      no recommendations at all (N1).
    * NIHSS applied no per-item maximum, so one mistyped field (`4_facial: 30`) went straight
      into the total, and unexamined items scored 0 (N20).
    * EDSS invented a step from `max(fs_scores)`, assigning EDSS 6.0 to a blind but ambulatory
      patient, and mislabelled its 6.0-6.5 band (N8).
    * HIT-6 summed whatever numbers it was sent with no check that they were HIT-6 weights,
      so six 1-5 Likert values read as "No significant impact" (N9).
    * ALSFRS-R summed only present items while its bulbar check `sum(... get(k, 4) ...) < 9`
      reported "Prominent bulbar involvement" when no bulbar item was submitted at all (N5).
    * ASPECTS counted only values equal to the string "0", so `{"M1": true}` -- a region
      marked affected -- scored 10/10, and a 12-key payload could go negative (N17).
    """
    return _scale_result_to_dict(calculate_scale(scale, items))


def _calculate_nihss(items: dict) -> dict:
    """NIHSS total from item scores (see NIHSSCalculator for the item names)."""
    return _calculate("nihss", items)


def _calculate_gcs(items: dict) -> dict:
    """GCS from eye, verbal and motor components; all three are required."""
    return _calculate("gcs", items)


def _calculate_moca(items: dict) -> dict:
    """MoCA from the seven domain sub-scores, with the published education adjustment."""
    return _calculate("moca", items)


def _calculate_updrs(items: dict) -> dict:
    """MDS-UPDRS Part III motor score from item scores."""
    return _calculate("updrs", items)


def _calculate_edss(items: dict) -> dict:
    """EDSS: interprets the examiner's step, cross-checked against the FS grades."""
    return _calculate("edss", items)


def _calculate_mrs(items: dict) -> dict:
    """Modified Rankin Scale grade."""
    return _calculate("mrs", items)


def _calculate_hit6(items: dict) -> dict:
    """HIT-6 from exactly six item responses."""
    return _calculate("hit6", items)


def _calculate_alsfrs(items: dict) -> dict:
    """ALSFRS-R from the twelve item scores."""
    return _calculate("alsfrs", items)


def _calculate_aspects(items: dict) -> dict:
    """ASPECTS from a per-region map (1 = normal, 0 = early ischemic change) or a list."""
    return _calculate("aspects", items)


def _calculate_hoehn_yahr(items: dict) -> dict:
    """Hoehn & Yahr stage."""
    return _calculate("hoehn_yahr", items)


_SCALE_CALCULATORS = {name: (lambda items, _n=name: _calculate(_n, items)) for name in SUPPORTED_SCALES}


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
                "collection": _as_row(r).get("collection", "unknown"),
                "text": _as_row(r).get("content", _as_row(r).get("text", "")),
                "score": _as_row(r).get("score", 0.0),
                "metadata": _as_row(r).get("metadata", {}),
            }
            for r in results
        ]
    except Exception as exc:
        logger.warning(f"Search failed: {exc}")
        evidence = []

    # Generate answer via LLM
    llm = _get_llm(req)
    answer = "Search completed. See evidence passages below."
    # Confidence used to start at 0.5 and rise 0.05 per hit, so an answer with NO evidence
    # still reported 0.5 — a coin flip presented as moderate confidence. It is now 0.0
    # until there is something behind the answer.
    confidence = 0.0
    grounded = bool(evidence)
    guidelines_cited = []

    if llm:
        context = "\n\n".join(e["text"] for e in evidence if e["text"])
        if evidence:
            prompt = (
                f"Neurology question: {request.question}\n\n"
                f"Retrieved evidence:\n{context}\n\n"
                f"Provide a comprehensive answer citing specific evidence. "
                f"Include relevant AAN, AHA/ASA, ILAE, or IHS guideline references."
            )
        else:
            # Asking for "specific evidence" with an empty context invited invented citations.
            prompt = (
                f"Neurology question: {request.question}\n\n"
                f"NO evidence was retrieved from the knowledge base for this question. "
                f"Answer from general knowledge only. Do NOT cite specific studies, "
                f"guidelines or figures as though they had been retrieved, and state plainly "
                f"that no source material was available."
            )
        if request.patient_context:
            prompt += f"\n\nPatient context: {request.patient_context}"

        try:
            answer = llm.generate(prompt)
            if evidence:
                confidence = min(0.95, 0.5 + len(evidence) * 0.05)
            else:
                confidence = 0.0
                answer = (
                    "No evidence was retrieved from the knowledge base for this question, so "
                    "the following is the model's general knowledge and is not grounded in "
                    "this platform's sources.\n\n" + answer
                )
            # Extract guideline references from answer
            for gline in ["AAN Guideline", "AHA/ASA", "ILAE", "IHS", "ICHD-3",
                          "NCCN", "AAN Practice", "McDonald Criteria", "MDS-UPDRS"]:
                # A name the model happened to write is not a citation. It counts only if
                # it also appears in the retrieved evidence this answer was built from.
                if gline.lower() in answer.lower() and gline.lower() in context.lower():
                    guidelines_cited.append(gline)
        except Exception as exc:
            logger.warning(f"LLM generation failed: {exc}")

    return QueryResponse(
        answer=answer,
        evidence=evidence,
        guidelines_cited=guidelines_cited,
        confidence=confidence,
        grounded=grounded,
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
                collection=_as_row(r).get("collection", "unknown"),
                text=_as_row(r).get("content", _as_row(r).get("text", "")),
                score=_as_row(r).get("score", 0.0),
                metadata=_as_row(r).get("metadata", {}),
            )
            for r in results
            if _as_row(r).get("score", 0.0) >= request.threshold
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

    try:
        result = _calculate(request.scale_name, request.items)
    except KeyError:
        raise HTTPException(
            status_code=400,
            detail=(
                f"Unknown scale '{request.scale_name}'. "
                f"Supported scales: {SUPPORTED_SCALES}"
            ),
        )
    except ScaleDataError as exc:
        # The items cannot produce a score for this instrument. The message names what is
        # missing or unusable; it is not swallowed into a generic "calculation error",
        # because the caller needs to know which input to supply.
        logger.info("Scale %s refused: %s", request.scale_name, exc)
        raise HTTPException(status_code=422, detail=str(exc))
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
        items_scored=result.get("items_scored"),
        items_expected=result.get("items_expected"),
        data_complete=result.get("complete", True),
        data_warnings=result.get("warnings", []),
    )


@router.post("/stroke/triage", response_model=StrokeTriageResponse)
async def stroke_triage(request: StrokeTriageRequest, req: Request):
    """Acute stroke triage with NIHSS, thrombolysis/thrombectomy eligibility.

    Four defects were fixed here. All four were reachable from the shipped UI.

    **N2 -- exclusions that only appended a sentence.** `SBP > 185` and `glucose < 50 or
    > 400` each appended a reasoning string and left `tpa_eligible` at `True`. A request with
    `{nihss_total: 12, onset_time_hours: 2.0, blood_pressure_systolic: 228, blood_glucose: 38}`
    returned `tpa_eligible: true` whose own `tpa_reasoning` read "SBP 228 -- must lower to
    <185 before tPA". SBP > 185/110 and glucose < 50 mg/dL are exclusions in the AHA/ASA 2019
    acute ischemic stroke guideline this endpoint already names in `guidelines_cited`; the
    code now applies the threshold it was already printing.

    **N2 -- the thrombolytic was ordered ahead of the CT.** `recommendations.insert(0, ...)`
    put "Administer IV alteplase" above "Stat CT head to rule out hemorrhage", which is read
    top-down by whoever is holding the phone. Alteplase into an undiagnosed hemorrhage is
    lethal. The CT is now unconditionally first and the thrombolytic is appended, explicitly
    conditioned on the CT having excluded hemorrhage.

    **N2 -- unknown onset cleared tPA.** `if request.onset_time_hours is not None and ... >
    4.5` fired nothing for `onset_time_hours: null`. The fix that landed in
    `clinical_workflows.py` was never applied to this copy, which is the one the UI calls.

    **N3/N20 -- unmeasured contraindications passed.** An absent BP, glucose or NIHSS is not a
    normal BP, glucose or NIHSS. Each missing input that gates thrombolysis now yields an
    `indeterminate` determination naming the input, rather than an eligible verdict.

    **N6 -- `prior_mrs` was collected, validated `ge=0, le=5`, documented, and never read.**
    `{nihss_total: 10, lvo_suspected: true, prior_mrs: 5}` -- a patient bedridden before the
    stroke -- returned `thrombectomy_eligible: true`. `age` was likewise unused.
    """
    _increment_metric(req, "workflow_requests_total")

    missing_inputs: List[str] = []

    # -- NIHSS -----------------------------------------------------------------
    nihss_complete = True
    nihss_warnings: List[str] = []
    if request.nihss_items:
        result = _calculate_nihss(request.nihss_items)
        nihss_score = int(result["total_score"])
        nihss_severity = result["severity_category"]
        nihss_complete = bool(result.get("complete", True))
        nihss_warnings = list(result.get("warnings", []))
    elif request.nihss_total is not None:
        nihss_score = request.nihss_total
        if nihss_score <= 4:
            nihss_severity = "minor"
        elif nihss_score <= 15:
            nihss_severity = "moderate"
        elif nihss_score <= 20:
            # The 16-20 "moderate-to-severe" band exists in _calculate_nihss and in
            # knowledge.py and was missing only here, so the same NIHSS 18 came back
            # "moderate_severe" when submitted as items and "severe" when submitted as a
            # total (N18).
            nihss_severity = "moderate_severe"
        elif nihss_score <= 25:
            nihss_severity = "severe"
        else:
            nihss_severity = "very_severe"
    else:
        nihss_score = 0
        nihss_severity = "not_assessed"
        nihss_complete = False
        missing_inputs.append("nihss_items or nihss_total")

    # -- tPA eligibility -------------------------------------------------------
    # Starts indeterminate. It becomes "eligible" only when every gating input was supplied
    # and none of them excludes thrombolysis.
    tpa_reasoning: List[str] = []
    tpa_excluded = False
    tpa_blocked_on_missing_data = False

    def _exclude(reason: str) -> None:
        nonlocal tpa_excluded
        tpa_excluded = True
        tpa_reasoning.append(reason)

    def _needs(field: str, reason: str) -> None:
        nonlocal tpa_blocked_on_missing_data
        tpa_blocked_on_missing_data = True
        missing_inputs.append(field)
        tpa_reasoning.append(f"INSUFFICIENT DATA: {reason}")

    if request.onset_time_hours is None:
        _needs(
            "onset_time_hours",
            "time of symptom onset not established. tPA requires a known onset within 4.5h, "
            "or DWI-FLAIR mismatch on MRI for wake-up/unwitnessed stroke (WAKE-UP, NEJM 2018)",
        )
    elif request.onset_time_hours > 4.5:
        _exclude(f"Onset {request.onset_time_hours:.1f}h ago -- exceeds 4.5h tPA window")
    else:
        tpa_reasoning.append(f"Within tPA window ({request.onset_time_hours:.1f}h since onset)")

    if not nihss_complete and request.nihss_items:
        _needs(
            "nihss_items (incomplete)",
            f"NIHSS examination incomplete -- the alteplase threshold (>= 4) cannot be "
            f"applied to a partial exam",
        )
    elif nihss_score == 0 and nihss_severity != "not_assessed":
        _exclude("NIHSS 0 -- no measurable deficit")
    elif nihss_severity == "not_assessed":
        _needs("nihss_items or nihss_total", "no NIHSS recorded")
    elif nihss_score < 4:
        tpa_reasoning.append(f"NIHSS {nihss_score} -- mild deficit, weigh risk/benefit")

    if request.anticoagulant_use:
        # INR / DOAC level is not collected by this endpoint, so anticoagulation cannot be
        # cleared here at all.
        _exclude(
            "On anticoagulation and no INR/DOAC level available -- tPA contraindicated "
            "pending coagulation results"
        )
    elif request.anticoagulant_use is None:
        _needs("anticoagulant_use", "anticoagulant use not recorded")

    if request.blood_pressure_systolic is None:
        _needs("blood_pressure_systolic", "blood pressure not measured")
    elif request.blood_pressure_systolic > 185:
        # AHA/ASA 2019 Early Management of Acute Ischemic Stroke: BP must be lowered to
        # < 185/110 mmHg before IV alteplase. Above that it is an exclusion, not a note.
        _exclude(
            f"SBP {request.blood_pressure_systolic} mmHg exceeds 185 mmHg -- NOT eligible "
            f"until BP is lowered to <185/110; re-triage after treatment"
        )

    if request.blood_glucose is None:
        _needs("blood_glucose", "blood glucose not measured")
    elif request.blood_glucose < 50:
        _exclude(
            f"Blood glucose {request.blood_glucose} mg/dL < 50 -- hypoglycaemia can mimic "
            f"stroke; correct and reassess before tPA"
        )
    elif request.blood_glucose > 400:
        _exclude(
            f"Blood glucose {request.blood_glucose} mg/dL > 400 -- correct before tPA"
        )

    if tpa_excluded:
        tpa_determination = "not_eligible"
    elif tpa_blocked_on_missing_data:
        tpa_determination = "indeterminate"
    else:
        tpa_determination = "eligible"
    tpa_eligible = tpa_determination == "eligible"

    # -- Thrombectomy eligibility ---------------------------------------------
    thrombectomy_reasoning: List[str] = []
    thrombectomy_excluded = False
    thrombectomy_indeterminate = False

    if request.lvo_suspected is None:
        thrombectomy_indeterminate = True
        missing_inputs.append("lvo_suspected")
        thrombectomy_reasoning.append(
            "INSUFFICIENT DATA: no vessel imaging result. Thrombectomy requires a "
            "demonstrated ICA/M1 (or basilar) occlusion -- obtain CTA head/neck"
        )
    elif not request.lvo_suspected:
        thrombectomy_excluded = True
        thrombectomy_reasoning.append("No large vessel occlusion -- not a thrombectomy candidate")
    else:
        thrombectomy_reasoning.append(
            "LVO suspected on the submitted information -- confirm ICA/M1 occlusion on CTA/MRA "
            "before activating the neurointerventional team"
        )

    if not nihss_complete:
        thrombectomy_indeterminate = True
        thrombectomy_reasoning.append(
            "INSUFFICIENT DATA: the thrombectomy threshold (NIHSS >= 6) cannot be applied to "
            "an incomplete or absent NIHSS"
        )
    elif nihss_score < 6:
        thrombectomy_excluded = True
        thrombectomy_reasoning.append(
            f"NIHSS {nihss_score} < 6 -- below the threshold used by this endpoint; "
            f"consider advanced imaging if clinical-imaging mismatch is suspected"
        )

    if request.onset_time_hours is None:
        thrombectomy_indeterminate = True
        thrombectomy_reasoning.append(
            "INSUFFICIENT DATA: onset time unknown -- extended-window selection requires "
            "perfusion imaging (DAWN/DEFUSE-3 criteria)"
        )
    elif request.onset_time_hours > 24:
        thrombectomy_excluded = True
        thrombectomy_reasoning.append("Beyond 24h window for thrombectomy")

    if request.ct_aspects_score is None:
        thrombectomy_indeterminate = True
        missing_inputs.append("ct_aspects_score")
        thrombectomy_reasoning.append(
            "INSUFFICIENT DATA: no ASPECTS recorded -- an unscored CT is not a favourable CT"
        )
    elif request.ct_aspects_score < 6:
        thrombectomy_excluded = True
        thrombectomy_reasoning.append(
            f"ASPECTS {request.ct_aspects_score} < 6 -- large established infarct"
        )

    # Pre-stroke functional status. Declared, validated and never read until now (N6).
    if request.prior_mrs is None:
        thrombectomy_indeterminate = True
        missing_inputs.append("prior_mrs")
        thrombectomy_reasoning.append(
            "INSUFFICIENT DATA: pre-stroke mRS not recorded -- every thrombectomy trial "
            "enrolled on pre-stroke functional independence"
        )
    elif request.prior_mrs >= 3:
        # Outside the enrollment range of every trial behind the recommendation: the 2015
        # trials enrolled pre-stroke mRS 0-1, DAWN and DEFUSE-3 enrolled 0-2.
        thrombectomy_excluded = True
        thrombectomy_reasoning.append(
            f"Pre-stroke mRS {request.prior_mrs} -- outside the enrollment range of the "
            f"trials behind this recommendation (0-1 in the 2015 trials, 0-2 in "
            f"DAWN/DEFUSE-3)"
        )
    elif request.prior_mrs == 2:
        thrombectomy_reasoning.append(
            "Pre-stroke mRS 2 -- within the DAWN/DEFUSE-3 range but outside the 2015 trials' "
            "(0-1); confirm with the treating neurointerventionalist"
        )

    if thrombectomy_excluded:
        thrombectomy_determination = "not_eligible"
    elif thrombectomy_indeterminate:
        thrombectomy_determination = "indeterminate"
    else:
        thrombectomy_determination = "eligible"
    thrombectomy_eligible = thrombectomy_determination == "eligible"

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

    # The non-contrast CT is first in this list and stays first: nothing that could be given
    # to a patient with an intracranial hemorrhage may be ordered above the study that rules
    # one out.
    recommendations = [
        "Stat CT head to rule out hemorrhage -- BEFORE any thrombolytic",
        "CBC, BMP, coagulation studies, troponin",
        "12-lead ECG for atrial fibrillation screening",
        "Continuous telemetry monitoring",
    ]
    if not nihss_complete and request.nihss_items:
        recommendations.append(
            "Complete the NIHSS examination: the recorded items give a lower bound only"
        )
    if tpa_eligible and nihss_score >= 4:
        recommendations.append(
            "AFTER the CT excludes hemorrhage: administer IV alteplase 0.9 mg/kg (max 90 mg), "
            "10% as bolus and the remainder over 60 minutes"
        )
    elif tpa_determination == "indeterminate":
        recommendations.append(
            "tPA eligibility INDETERMINATE -- obtain: " + ", ".join(sorted(set(missing_inputs)))
        )
    if thrombectomy_eligible:
        recommendations.append(
            "AFTER the CT excludes hemorrhage and the occlusion is confirmed on vessel "
            "imaging: activate the neurointerventional team for mechanical thrombectomy"
        )

    # Model-generated prose. It was previously appended to `evidence` as
    # `{"source": "llm", "text": ..., "score": 0.8}` (G3): `evidence` is a key the honesty
    # gate deliberately skips, because retrieved passages must not be rewritten -- so this
    # text was never scanned, never got the decision-support disclaimer, and was handed to
    # the claim verifier as the ground truth against which to check the answer that produced
    # it. The UI rendered it as "Evidence (1 passages) -- [llm] Score: 0.800", which is
    # indistinguishable from a retrieved PubMed passage. It is model knowledge, it was never
    # retrieved and it was never scored, so it carries no relevance number and lives under a
    # key the gate reads.
    evidence: List[dict] = []
    discussion = None
    llm = _get_llm(req)
    if llm and request.clinical_notes:
        try:
            prompt = (
                f"Acute stroke triage:\n"
                f"NIHSS: {nihss_score}, Onset: {request.onset_time_hours}h\n"
                f"LVO suspected: {request.lvo_suspected}\n"
                f"ASPECTS: {request.ct_aspects_score}\n"
                f"Notes: {request.clinical_notes}\n\n"
                f"Provide additional clinical considerations. Nothing was retrieved from the "
                f"knowledge base for this case, so do not cite studies, guidelines or figures "
                f"as though they had been looked up."
            )
            discussion = (
                "Model knowledge, not retrieved from this knowledge base:\n\n"
                + llm.generate(prompt)
            )
        except Exception as exc:
            logger.warning(f"Stroke triage LLM failed: {exc}")

    return StrokeTriageResponse(
        nihss_score=nihss_score,
        nihss_severity=nihss_severity,
        tpa_eligible=tpa_eligible,
        tpa_reasoning=tpa_reasoning + nihss_warnings,
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
        tpa_determination=tpa_determination,
        thrombectomy_determination=thrombectomy_determination,
        missing_inputs=sorted(set(missing_inputs)),
        nihss_complete=nihss_complete,
        discussion=discussion,
    )


@router.post("/dementia/evaluate", response_model=DementiaEvaluateResponse)
async def dementia_evaluate(request: DementiaEvaluateRequest, req: Request):
    """Comprehensive dementia evaluation with cognitive screening and differential."""
    _increment_metric(req, "workflow_requests_total")

    # MoCA calculation. `education_years` is a declared request field that this endpoint
    # never passed to the calculator -- it relied on the client smuggling it inside
    # `moca_items`, which is exactly how N1 inflated the total. It is now passed properly and
    # cannot reach the sum.
    moca_score = request.moca_total
    moca_interp = None
    moca_complete = True
    data_warnings: List[str] = []
    if request.moca_items:
        items = dict(request.moca_items)
        if request.education_years is not None:
            items["education_years"] = request.education_years
        result = _calculate_moca(items)
        moca_score = int(result["total_score"])
        moca_interp = result["interpretation"]
        moca_complete = bool(result.get("complete", True))
        data_warnings.extend(result.get("warnings", []))
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
    # A cholinesterase inhibitor used to be recommended for any MoCA < 26 -- i.e. for MCI --
    # in a response citing the AAN 2018 MCI Practice Guideline, whose MCI recommendation is
    # that clinicians *may choose not to* prescribe cholinesterase inhibitors (Level B) and
    # must first discuss the absence of evidence (Level A). The guidance contradicted its own
    # citation. Scored MCI now gets the guideline's actual position. (Petersen et al. 2018,
    # PMID 29282327; PubMed labels that guideline RETIRED -- confirm its current status before
    # quoting it elsewhere.)
    if moca_score is not None and not moca_complete:
        treatment.append(
            "No pharmacotherapy recommendation: the MoCA was incomplete, so the patient has "
            "no interpretable cognitive stage"
        )
    elif moca_score is not None and moca_score < 18:
        treatment.append("Cholinesterase inhibitor (donepezil 5-10mg, rivastigmine, or galantamine)")
        treatment.append("Consider memantine addition")
    elif moca_score is not None and moca_score < 26:
        treatment.append(
            "MCI range: cholinesterase inhibitors are not approved for MCI. Per the AAN 2018 "
            "MCI guideline, clinicians may choose not to prescribe them and should discuss "
            "the absence of evidence of benefit with the patient first"
        )

    recommendations = []
    if moca_score is not None and moca_score < 26:
        recommendations.append("Assess driving safety and functional independence")
        recommendations.append("Connect with dementia support resources and caregiver education")
    if not moca_complete:
        recommendations.insert(
            0,
            "Complete the MoCA before staging cognition or selecting therapy -- "
            "unrecorded domains are not normal domains",
        )
    recommendations.append("Follow-up cognitive assessment in 6-12 months")

    # Model prose -- see the G3 note in stroke_triage. Not evidence, no relevance score.
    evidence: List[dict] = []
    discussion = None
    llm = _get_llm(req)
    if llm and request.clinical_notes:
        try:
            prompt = (
                f"Dementia evaluation:\n"
                f"MoCA: {moca_score}, Age: {request.age}\n"
                f"Symptoms: {', '.join(request.dominant_symptoms)}\n"
                f"Motor: {request.motor_features}\n"
                f"Notes: {request.clinical_notes}\n\n"
                f"Provide differential diagnosis considerations. Nothing was retrieved from "
                f"the knowledge base for this case, so do not cite studies, guidelines or "
                f"figures as though they had been looked up."
            )
            discussion = (
                "Model knowledge, not retrieved from this knowledge base:\n\n"
                + llm.generate(prompt)
            )
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
        moca_complete=moca_complete,
        data_warnings=data_warnings,
        discussion=discussion,
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
    if request.clinical_diagnosis:
        # N19(a): `syndrome = "dravet"` was reachable only through a pathogenic SCN1A result,
        # so a clinically diagnosed Dravet patient on carbamazepine with no genetic testing in
        # the payload got no contraindication warning at all. Sodium-channel blockers worsen
        # Dravet seizures. A clinical diagnosis is now sufficient to name the syndrome.
        dx = request.clinical_diagnosis.lower()
        for name, key in (
            ("dravet", "dravet"),
            ("juvenile myoclonic", "juvenile_myoclonic"),
            ("childhood absence", "childhood_absence"),
            ("west", "west"),
            ("temporal lobe", "temporal_lobe"),
        ):
            if name in dx:
                syndrome = key
                break
    if request.genetic_results and "scn1a" in request.genetic_results.lower():
        if any(w in request.genetic_results.lower() for w in ("pathogenic", "mutation", "variant")):
            syndrome = syndrome or "dravet"
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

    # Surgical candidacy.
    # N19(c): this fired on `len(current_aeds) >= 2 and seizure_frequency`, so a seizure-free
    # patient on two concurrent ASMs was told "Drug-resistant epilepsy criteria may be met"
    # under a cited drug-resistant-epilepsy guideline. The ILAE definition requires FAILURE of
    # two adequate trials, which is what `failed_aeds` records.
    surgical_candidacy = None
    if len(request.failed_aeds) >= 2 and request.seizure_frequency:
        surgical_candidacy = (
            f"Drug-resistant epilepsy criteria may be met ({len(request.failed_aeds)} ASMs "
            f"failed an adequate trial, ongoing seizures) -- consider epilepsy surgery "
            f"evaluation"
        )
    elif len(request.current_aeds) >= 2 and request.seizure_frequency:
        surgical_candidacy = (
            "Drug-resistance not established: two or more ASMs are currently prescribed, but "
            "the ILAE definition requires failure of two adequately chosen, tolerated and "
            "dosed trials. Record them in `failed_aeds`."
        )

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

    # Syndrome-specific drug contraindications. This endpoint previously only *replaced* the
    # ASM list and never warned that a drug the patient was already taking was contraindicated
    # (N19a).
    _CONTRAINDICATED_ASMS = {
        "dravet": (
            ["carbamazepine", "oxcarbazepine", "phenytoin", "lamotrigine", "vigabatrin"],
            "sodium-channel blockers and vigabatrin worsen seizures in Dravet syndrome",
        ),
        "juvenile_myoclonic": (
            ["carbamazepine", "oxcarbazepine", "phenytoin", "gabapentin", "pregabalin",
             "vigabatrin", "tiagabine"],
            "these agents can aggravate myoclonic and absence seizures in JME",
        ),
        "childhood_absence": (
            ["carbamazepine", "oxcarbazepine", "phenytoin", "gabapentin", "vigabatrin",
             "tiagabine"],
            "these agents can aggravate absence seizures",
        ),
    }
    contraindication_warnings: List[str] = []
    if syndrome in _CONTRAINDICATED_ASMS:
        banned, why = _CONTRAINDICATED_ASMS[syndrome]
        current = " ".join(request.current_aeds).lower()
        hits = [drug for drug in banned if drug in current]
        if hits:
            contraindication_warnings.append(
                f"DANGER: {', '.join(hits)} CONTRAINDICATED in {syndrome.replace('_', ' ')} "
                f"-- {why}. Review urgently."
            )
        contraindication_warnings.append(
            f"Avoid in {syndrome.replace('_', ' ')}: {', '.join(banned)} -- {why}"
        )

    recommendations = [
        "Routine and/or prolonged EEG monitoring",
        "MRI brain with epilepsy protocol",
    ]
    recommendations.extend(contraindication_warnings)
    if len(request.current_aeds) >= 2:
        recommendations.append("Refer to comprehensive epilepsy center")
    if surgical_candidacy:
        recommendations.append(surgical_candidacy)

    # Model prose -- see the G3 note in stroke_triage. Not evidence, no relevance score.
    evidence: List[dict] = []
    discussion = None
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
                f"Provide classification and management considerations. Nothing was retrieved "
                f"from the knowledge base for this case, so do not cite studies, guidelines or "
                f"figures as though they had been looked up."
            )
            discussion = (
                "Model knowledge, not retrieved from this knowledge base:\n\n"
                + llm.generate(prompt)
            )
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
        contraindication_warnings=contraindication_warnings,
        recommendations=recommendations,
        guidelines_cited=[
            "ILAE 2017 Seizure Classification",
            "ILAE 2017 Epilepsy Classification",
            "AAN/AES 2018 Drug-Resistant Epilepsy Guideline",
        ],
        evidence=evidence,
        discussion=discussion,
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
        # N12: `else: "MGMT unmethylated"` reported a tumor whose MGMT was never tested as
        # unmethylated, and routed it down the branch that lowers temozolomide expectations.
        # "Not tested" is now its own answer.
        if "mgmt_methylated" in markers:
            treatment.append("MGMT methylated -- favorable response to temozolomide expected")
        elif "mgmt_unmethylated" in markers:
            treatment.append("MGMT unmethylated -- consider clinical trial alternatives")
        else:
            treatment.append(
                "MGMT promoter methylation NOT TESTED -- order it; temozolomide benefit "
                "cannot be predicted without it, and absence of a result is not an "
                "unmethylated result"
            )
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
    if request.kps is None:
        recommendations.append(
            "Karnofsky Performance Status not recorded -- record it; treatment intensity and "
            "the palliative-care decision both turn on it"
        )
    elif request.kps < 70:
        # NOTE: this endpoint uses KPS < 70 while clinical_workflows.py uses KPS < 60 for the
        # same decision. Neither cut-point is cited in this repo and no primary source was
        # retrieved, so both values are left as they are and the divergence is recorded here
        # rather than silently resolved in one direction.
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
    missing_inputs: List[str] = []
    if request.edss_items:
        items = dict(request.edss_items)
        if request.edss_score is not None:
            items.setdefault("edss_step", request.edss_score)
        result = _calculate_edss(items)
        if result.get("complete", True):
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

    # NEDA-3 status.
    #
    # N7(b): this read `if edss_score is not None and edss_score > 3.0: neda_3 = False;
    # edss_progression = True`. Progression is a CHANGE, and this endpoint had no prior-EDSS
    # field at all, so it was wrong in both directions: a patient stable at EDSS 4.0 for ten
    # years with no relapses and no new lesions was reported as progressing, with
    # `treatment_response: "suboptimal"` and a DMT-escalation suggestion; while a patient who
    # genuinely went 1.0 -> 3.0 got `edss_progression: false` and `neda_3_achieved: true`.
    # The absolute test is gone and the confirmed-progression rule already coded in
    # clinical_workflows.py (>= 1.0 when baseline <= 5.5, >= 0.5 when baseline >= 6.0) is
    # applied to the recorded prior EDSS.
    #
    # N7: NEDA-3 is also no longer asserted from silence. "No evidence of disease activity"
    # requires evidence: with no relapse count, no MRI result and no prior EDSS, the status is
    # not assessable, which is not the same as achieved.
    if request.relapse_count_2yr is None:
        missing_inputs.append("relapse_count_2yr")
    if request.new_t2_lesions is None and request.gad_enhancing_lesions is None:
        missing_inputs.append("new_t2_lesions / gad_enhancing_lesions")
    if request.edss_prior is None:
        missing_inputs.append("edss_prior")

    edss_progression: Optional[bool] = None
    if request.edss_prior is not None and edss_score is not None:
        delta = edss_score - request.edss_prior
        if request.edss_prior <= 5.5:
            edss_progression = delta >= 1.0
        else:
            edss_progression = delta >= 0.5

    neda_components = {
        "no_relapses": None if request.relapse_count_2yr is None else relapses == 0,
        "no_new_lesions": (
            None
            if (request.new_t2_lesions is None and request.gad_enhancing_lesions is None)
            else new_lesions == 0
        ),
        "no_edss_progression": None if edss_progression is None else not edss_progression,
    }
    if any(v is False for v in neda_components.values()):
        neda_3: Optional[bool] = False
        neda_assessable = True
    elif all(v is True for v in neda_components.values()):
        neda_3 = True
        neda_assessable = True
    else:
        neda_3 = None
        neda_assessable = False

    neda_status = {
        "neda_3_achieved": neda_3,
        "neda_3_assessable": neda_assessable,
        **neda_components,
    }
    if not neda_assessable:
        neda_status["note"] = (
            "NEDA-3 NOT ASSESSABLE -- missing: " + ", ".join(missing_inputs)
            + ". NEDA-3 means no evidence of disease activity, which requires evidence; "
            "unreported relapses, unperformed MRI and an unknown prior EDSS are not absence "
            "of activity."
        )

    # DMT assessment
    dmt_assessment = {
        "current_dmt": request.current_dmt or "none",
        "treatment_response": (
            "not_assessable" if neda_3 is None else ("adequate" if neda_3 else "suboptimal")
        ),
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
    if missing_inputs:
        recommendations.insert(
            0,
            "Disease-activity assessment is incomplete -- record: "
            + ", ".join(missing_inputs),
        )

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

    # Motor subtype.
    #
    # N13(b): the PIGD key set was matched by substring -- `"postural" in k.lower()` -- which
    # also matched `3.15a_postural_tremor_r` and `3.15b_postural_tremor_l`. A pure-tremor
    # patient therefore contributed tremor points to the PIGD sum: with
    # `{3.15a_postural_tremor_r: 4, 3.17a_rest_tremor_rue: 4, 3.10_gait: 0,
    # 3.12_postural_stability: 0}` the PIGD keys came out as
    # `['3.15a_postural_tremor_r', '3.10_gait', '3.12_postural_stability']`. The item sets are
    # now named explicitly on UPDRSCalculator and shared with the workflow.
    #
    # OPEN (uncited, unchanged): this compares raw sums of 10 tremor items against 5 PIGD
    # items, while clinical_workflows.py compares means with 1.5 / 0.67 ratio cut-offs. The
    # MDS-UPDRS-based cut-offs (>= 1.15 / <= 0.90) and the canonical PIGD item set come from
    # Stebbins et al., Mov Disord 2013, which was not retrieved, so neither algorithm's
    # numbers are changed here.
    motor_subtype = "indeterminate"
    if updrs_score is not None and request.updrs_items:
        from src.clinical_scales import _resolve_items  # one item-name resolver, not two

        resolved, _ = _resolve_items(
            request.updrs_items, UPDRSCalculator.ITEMS, UPDRSCalculator.ALIASES
        )
        tremor_sum = sum(resolved.get(k, 0) for k in UPDRSCalculator.TREMOR_ITEMS)
        pigd_sum = sum(resolved.get(k, 0) for k in UPDRSCalculator.PIGD_ITEMS)
        tremor_seen = any(k in resolved for k in UPDRSCalculator.TREMOR_ITEMS)
        pigd_seen = any(k in resolved for k in UPDRSCalculator.PIGD_ITEMS)
        if not (tremor_seen and pigd_seen):
            # Subtyping compares two item groups; it cannot be done when one group was never
            # scored, and the group that was not scored is not a group of zeros.
            motor_subtype = "not_assessable (tremor and PIGD items both required)"
        elif tremor_sum > pigd_sum:
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

    # HIT-6. The calculator now refuses a short or unrecognised response set rather than
    # padding it with "never" -- the best possible answer on every item (N9) -- so the refusal
    # is surfaced as a 422 naming the problem rather than as a reassuring score.
    hit6_score = None
    hit6_interp = None
    if request.hit6_items:
        try:
            result = _calculate_hit6(request.hit6_items)
        except ScaleDataError as exc:
            raise HTTPException(status_code=422, detail=f"HIT-6: {exc}")
        hit6_score = int(result["total_score"])
        hit6_interp = result["interpretation"]

    # Red flag assessment
    red_flag_assessment = {
        "flags_present": request.red_flags,
        "urgent_workup_needed": len(request.red_flags) > 0,
    }
    if request.red_flags:
        red_flag_assessment["action"] = "Urgent neuroimaging and further evaluation required"

    # Medication-overuse headache risk.
    #
    # N15: the drug-class-specific threshold was collapsed to 15 days for everything, so a
    # patient on sumatriptan 12 days/month for 6 months -- who MEETS ICHD-3 8.2.2 -- was
    # reported as merely "approaching" the threshold, and the withdrawal-protocol
    # recommendation (gated on "High" in moh_risk) never fired.
    #
    # ICHD-3, the issuing body's own criteria
    # (https://ichd-3.org/8-headache-attributed-to-a-substance-or-its-withdrawal/
    #  8-2-medication-overuse-headache-moh/, sections 8.2.2, 8.2.4, 8.2.6):
    #   - simple analgesics / paracetamol / NSAIDs: >= 15 days/month for > 3 months
    #   - triptans, ergotamine, opioids:            >= 10 days/month for > 3 months
    #   - any combination of these classes:         >= 10 days/month in total
    _MOH_10_DAY_CLASSES = {
        "triptan": ["triptan", "sumatriptan", "rizatriptan", "zolmitriptan", "naratriptan",
                    "almotriptan", "frovatriptan", "eletriptan"],
        "ergotamine": ["ergot", "ergotamine", "dihydroergotamine", "dhe", "cafergot"],
        "opioid": ["opioid", "codeine", "tramadol", "oxycodone", "hydrocodone", "morphine",
                   "butalbital", "fioricet", "percocet", "vicodin"],
        "combination analgesic": ["excedrin", "combination analgesic", "butalbital"],
    }
    moh_risk = None
    if request.analgesic_days_per_month is not None:
        meds = " ".join(request.current_medications).lower()
        classes_present = sorted(
            label for label, names in _MOH_10_DAY_CLASSES.items()
            if any(name in meds for name in names)
        )
        if classes_present:
            threshold, basis = 10, f"{', '.join(classes_present)} (ICHD-3: >= 10 days/month)"
        elif request.current_medications:
            threshold, basis = 15, "simple analgesics/NSAIDs (ICHD-3: >= 15 days/month)"
        else:
            # The acute-medication class determines the threshold and was not supplied.
            # The lower threshold is not assumed; the ambiguity is reported.
            threshold, basis = None, None

        days = request.analgesic_days_per_month
        if threshold is None:
            moh_risk = (
                f"Not determinable -- {days} days/month of acute medication, but no "
                f"medication list was supplied. The ICHD-3 threshold is 10 days/month for "
                f"triptans, ergots, opioids and combination analgesics and 15 days/month for "
                f"simple analgesics, so the class is required."
            )
        elif days >= threshold:
            moh_risk = (
                f"High risk -- meets the ICHD-3 MOH day threshold ({days} >= {threshold} "
                f"days/month, {basis}). ICHD-3 also requires regular overuse for > 3 months; "
                f"confirm duration."
            )
        elif days >= threshold - 5:
            moh_risk = (
                f"Moderate risk -- {days} days/month, approaching the {threshold}-day "
                f"threshold for {basis}"
            )

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
        # N14: the differential was labelled "GBS (if acute)", which matched none of the
        # treatment branches below, so a GBS differential produced no treatment at all. The
        # tempo decides which of the two it is, and when it was not recorded the acute
        # possibility must stay on the list with its airway monitoring: GBS respiratory
        # failure is the non-recoverable error in this domain.
        tempo = (request.progression_rate or "").lower()
        acute_tempo = any(w in tempo for w in ("acute", "rapid", "days", "hyperacute"))
        chronic_tempo = any(w in tempo for w in ("chronic", "slow", "insidious", "months", "years"))
        if acute_tempo or not tempo:
            differential.append({
                "diagnosis": "Guillain-Barre syndrome",
                "likelihood": "high" if acute_tempo else "cannot be excluded",
                "rationale": (
                    "Acute demyelinating polyradiculoneuropathy"
                    if acute_tempo else
                    "Demyelinating pattern with no recorded progression rate -- the tempo "
                    "that distinguishes GBS from CIDP was not supplied, so GBS cannot be "
                    "excluded"
                ),
            })
        if chronic_tempo or not acute_tempo:
            differential.append({
                "diagnosis": "CIDP",
                "likelihood": "high" if chronic_tempo else "consider",
                "rationale": "Demyelinating neuropathy pattern",
            })
    if localization == "muscle" or (request.ck_level and request.ck_level > 1000):
        differential.append({"diagnosis": "Inflammatory myopathy", "likelihood": "moderate", "rationale": "Myopathic pattern with elevated CK"})
    if any(w in desc for w in ["fasciculation", "wasting", "upper motor", "lower motor"]):
        differential.append({"diagnosis": "ALS", "likelihood": "consider", "rationale": "Combined UMN/LMN features"})
    if not differential:
        differential.append({"diagnosis": "Neuromuscular disorder NOS", "likelihood": "pending", "rationale": "Further workup needed"})

    # ALSFRS-R
    alsfrs_score = None
    alsfrs_interp = None
    alsfrs_warnings: List[str] = []
    if request.alsfrs_items:
        result = _calculate_alsfrs(request.alsfrs_items)
        alsfrs_score = int(result["total_score"])
        alsfrs_interp = result["interpretation"]
        alsfrs_warnings = list(result.get("warnings", []))
        if not result.get("complete", True):
            # The route twin reported "Prominent bulbar involvement" when NO bulbar item had
            # been submitted, because `sum(int(items.get(k, 4)) ... if k in items)` summed to
            # 0 over an empty selection (N5). An incomplete instrument is now said to be one.
            alsfrs_warnings.insert(
                0,
                "ALSFRS-R incomplete: the total above is a lower bound, not the patient's "
                "score, and the functional-impairment bands do not apply to it.",
            )

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
        elif dx["diagnosis"] == "Guillain-Barre syndrome":
            treatment.extend([
                "IVIg 0.4 g/kg/day x 5 days, or plasmapheresis (5 exchanges) -- do not give both",
                "Admit for serial respiratory monitoring; corticosteroids are NOT effective in GBS",
            ])
        elif dx["diagnosis"] == "ALS":
            treatment.extend(["Riluzole 50mg BID", "Edaravone or tofersen (if SOD1+)", "Multidisciplinary ALS clinic"])

    recommendations = [
        "Neuromuscular specialist referral",
        "Pulmonary function testing (FVC)",
        "Swallow evaluation if bulbar symptoms",
    ]
    # N14: escalation used to require the literal string "rapid", so a patient described as
    # "acute" got no urgent recommendation at all -- and a patient whose tempo was never
    # recorded got none either.
    gbs_on_differential = any(
        dx["diagnosis"] == "Guillain-Barre syndrome" for dx in differential
    )
    tempo_text = (request.progression_rate or "").lower()
    if any(w in tempo_text for w in ("rapid", "acute", "days", "hyperacute")):
        recommendations.insert(0, "Urgent evaluation given rapid/acute progression")
    elif not tempo_text:
        recommendations.insert(
            0,
            "Progression rate NOT RECORDED -- record it. It is the feature that separates an "
            "acute (GBS) from a chronic (CIDP) demyelinating neuropathy and it changes the "
            "monitoring plan.",
        )
    if gbs_on_differential:
        recommendations.insert(
            0,
            "URGENT: monitor respiratory function (FVC q4h, negative inspiratory force, "
            "single-breath count) -- intubate if FVC < 20 mL/kg or declining rapidly. Also "
            "monitor for autonomic instability.",
        )

    return NeuromuscularEvaluateResponse(
        localization=localization,
        pattern_classification=pattern,
        differential_diagnosis=differential,
        alsfrs_score=alsfrs_score,
        alsfrs_interpretation=alsfrs_interp,
        data_warnings=alsfrs_warnings,
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
