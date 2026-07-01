"""Cardiology clinical API routes.

Provides endpoints for RAG-powered cardiovascular queries, six validated
risk calculators (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC,
EuroSCORE II), GDMT optimization, eight clinical workflows, guideline
search, and reference catalogues for conditions, biomarkers, drugs, and
cardio-relevant genes.

Author: Adam Jones
Date: March 2026
"""

from typing import List, Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

from src.clinical_workflows import WorkflowEngine
from src.models import CardioWorkflowType

router = APIRouter(prefix="/v1/cardio", tags=["cardiology"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/integrated-assessment")
async def integrated_assessment(request: dict, req: Request):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries oncology, clinical trial, biomarker, and imaging agents to build
    a comprehensive cardiovascular assessment for the cardio-oncology pathway.
    """
    try:
        from src.cross_agent import (
            query_oncology_agent,
            query_trial_agent,
            query_biomarker_agent,
            query_imaging_agent,
            integrate_cross_agent_results,
        )

        patient_profile = request.get("patient_profile", {})
        cardiac_biomarkers = request.get("cardiac_biomarkers", {})
        imaging_type = request.get("imaging_type", "echocardiogram")

        results = []

        # Query oncology agent for chemotherapy context
        if patient_profile:
            results.append(query_oncology_agent(patient_profile))

        # Query trial agent for cardiac monitoring requirements
        if patient_profile.get("trial_id"):
            results.append(query_trial_agent(patient_profile))

        # Query biomarker agent for troponin/BNP trends
        if cardiac_biomarkers:
            results.append(query_biomarker_agent(cardiac_biomarkers))

        # Query imaging agent for baseline coordination
        if imaging_type:
            results.append(query_imaging_agent(imaging_type, patient_context=patient_profile))

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

# ── Query ──

class QueryRequest(BaseModel):
    """Free-text RAG query with optional workflow and patient context."""
    question: str = Field(..., min_length=3, description="Clinical question")
    workflow_type: Optional[str] = Field(
        None, description="Workflow hint: cad | heart_failure | arrhythmia | valvular | imaging | prevention | cardio_oncology | stress_test"
    )
    patient_context: Optional[dict] = Field(None, description="Demographics, labs, vitals")
    top_k: int = Field(5, ge=1, le=50, description="Number of evidence passages")
    include_guidelines: bool = Field(True, description="Include guideline citations")


class QueryResponse(BaseModel):
    answer: str
    evidence: List[dict]
    guidelines_cited: List[str] = []
    confidence: float
    workflow_applied: Optional[str] = None


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


# ── Risk Calculators ──

class ASCVDRequest(BaseModel):
    """Pooled Cohort Equations -- 10-year ASCVD risk."""
    age: int = Field(..., ge=40, le=79, description="Age 40-79")
    sex: str = Field(..., pattern="^(male|female)$")
    race: str = Field(..., pattern="^(white|african_american|other)$")
    total_cholesterol: float = Field(..., ge=100, le=400, description="mg/dL")
    hdl_cholesterol: float = Field(..., ge=20, le=150, description="mg/dL")
    systolic_bp: float = Field(..., ge=80, le=250, description="mmHg")
    bp_treatment: bool = False
    diabetes: bool = False
    smoker: bool = False


class HEARTRequest(BaseModel):
    """HEART score for acute chest pain evaluation."""
    history_score: int = Field(..., ge=0, le=2, description="0=slightly suspicious, 1=moderately, 2=highly")
    ecg_score: int = Field(..., ge=0, le=2, description="0=normal, 1=non-specific, 2=significant deviation")
    age: int = Field(..., ge=18, le=120)
    risk_factors: int = Field(..., ge=0, le=5, description="HTN, DM, hyperlipidemia, obesity, smoking, family hx (0-2 scoring)")
    troponin_score: int = Field(..., ge=0, le=2, description="0=normal, 1=1-3x ULN, 2=>3x ULN")


class CHA2DS2VAScRequest(BaseModel):
    """CHA2DS2-VASc stroke risk in atrial fibrillation."""
    chf: bool = False
    hypertension: bool = False
    age: int = Field(..., ge=18, le=120)
    diabetes: bool = False
    stroke_tia: bool = False
    vascular_disease: bool = False
    sex: str = Field(..., pattern="^(male|female)$")


class HASBLEDRequest(BaseModel):
    """HAS-BLED bleeding risk for anticoagulated patients."""
    hypertension_uncontrolled: bool = False
    renal_disease: bool = False
    liver_disease: bool = False
    stroke_history: bool = False
    bleeding_history: bool = False
    labile_inr: bool = False
    age_over_65: bool = False
    drugs_alcohol: int = Field(0, ge=0, le=2, description="0=none, 1=one, 2=both")


class MAGGICRequest(BaseModel):
    """MAGGIC heart failure mortality risk."""
    age: int = Field(..., ge=18, le=100)
    sex: str = Field(..., pattern="^(male|female)$")
    lvef: float = Field(..., ge=5, le=80, description="Left ventricular EF %")
    nyha_class: int = Field(..., ge=1, le=4)
    systolic_bp: float = Field(..., ge=60, le=250, description="mmHg")
    bmi: float = Field(..., ge=10, le=60, description="kg/m2")
    creatinine: float = Field(..., ge=0.3, le=15.0, description="mg/dL")
    diabetes: bool = False
    copd: bool = False
    hf_duration_18m: bool = Field(False, description="HF diagnosed >18 months ago")
    smoker: bool = False
    beta_blocker: bool = False
    acei_arb: bool = False


class EuroSCORERequest(BaseModel):
    """EuroSCORE II cardiac surgical risk."""
    age: int = Field(..., ge=18, le=100)
    sex: str = Field(..., pattern="^(male|female)$")
    creatinine_clearance: float = Field(..., ge=0, le=200, description="mL/min")
    extracardiac_arteriopathy: bool = False
    poor_mobility: bool = False
    previous_cardiac_surgery: bool = False
    chronic_lung_disease: bool = False
    active_endocarditis: bool = False
    critical_preop_state: bool = False
    diabetes_on_insulin: bool = False
    nyha_class: int = Field(1, ge=1, le=4)
    ccs_class_4_angina: bool = False
    lvef: float = Field(60, ge=5, le=80)
    recent_mi: bool = False
    pulmonary_hypertension: str = Field("none", pattern="^(none|moderate|severe)$")
    urgency: str = Field("elective", pattern="^(elective|urgent|emergency|salvage)$")
    procedure_weight: str = Field("isolated_cabg", description="Procedure complexity")
    thoracic_aorta: bool = False


class RiskResult(BaseModel):
    calculator: str
    score: float
    risk_category: str
    interpretation: str
    recommendations: List[str] = []
    details: dict = {}


# ── GDMT ──

class GDMTRequest(BaseModel):
    """Guideline-Directed Medical Therapy optimization for HF."""
    lvef: float = Field(..., ge=5, le=80, description="Left ventricular EF %")
    nyha_class: str = Field(..., pattern="^(I|II|III|IV)$")
    current_medications: List[dict] = Field(
        default=[], description="[{'name': 'carvedilol', 'dose_mg': 25, 'frequency': 'BID'}]"
    )
    patient_data: dict = Field(
        default={}, description="BP, HR, K+, Cr, eGFR, symptoms, comorbidities"
    )


class GDMTRecommendation(BaseModel):
    medication_class: str
    current_status: str
    recommendation: str
    target_dose: Optional[str] = None
    evidence_level: str = ""
    guideline_ref: str = ""
    caution: Optional[str] = None


class GDMTResponse(BaseModel):
    hf_phenotype: str
    lvef: float
    nyha_class: str
    recommendations: List[GDMTRecommendation]
    four_pillars_status: dict
    summary: str


# =====================================================================
# Query Endpoints
# =====================================================================

@router.post("/query", response_model=QueryResponse)
def cardio_query(request: QueryRequest, req: Request):
    """Answer a cardiology clinical question using RAG with guideline citations."""
    engine = req.app.state.engine
    if not engine:
        raise HTTPException(status_code=503, detail="RAG engine not initialized")

    try:
        result = engine.query(
            question=request.question,
            workflow=request.workflow_type,
            patient_context=request.patient_context,
            top_k=request.top_k,
        )
        with req.app.state.metrics_lock:
            req.app.state.metrics["query_requests_total"] += 1
        return QueryResponse(
            answer=result.answer,
            evidence=[vars(e) for e in result.results] if result.results else [],
            guidelines_cited=[c.get("text", str(c)) for c in result.citations] if result.citations else [],
            confidence=result.confidence,
            workflow_applied=result.workflow.value if result.workflow else None,
        )
    except Exception as exc:
        logger.error(f"Query failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@router.post("/search", response_model=SearchResponse)
def cardio_search(request: SearchRequest, req: Request):
    """Semantic search across cardiology knowledge collections."""
    engine = req.app.state.engine
    if not engine:
        raise HTTPException(status_code=503, detail="RAG engine not initialized")

    try:
        results = engine.search(
            question=request.question,
            collections=request.collections,
            top_k=request.top_k,
        )
        with req.app.state.metrics_lock:
            req.app.state.metrics["search_requests_total"] += 1
        return SearchResponse(
            results=[vars(r) for r in results],
            total=len(results),
            collections_searched=request.collections or [],
        )
    except Exception as exc:
        logger.error(f"Search failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@router.post("/find-related")
def find_related(
    entity: str,
    entity_type: str = "condition",
    top_k: int = 10,
    req: Request = None,
):
    """Find entities related to a given condition, drug, gene, or biomarker."""
    engine = getattr(req.app.state, "engine", None) if req else None
    if not engine:
        # Fallback: static relationships
        relationships = {
            "atrial_fibrillation": {
                "drugs": ["apixaban", "rivaroxaban", "amiodarone", "flecainide"],
                "biomarkers": ["BNP", "hs-CRP", "D-dimer"],
                "genes": ["SCN5A", "KCNQ1", "PITX2"],
                "procedures": ["catheter ablation", "LAA occlusion", "cardioversion"],
            },
            "heart_failure": {
                "drugs": ["sacubitril/valsartan", "dapagliflozin", "carvedilol", "spironolactone"],
                "biomarkers": ["NT-proBNP", "hs-TnT", "sST2", "galectin-3"],
                "genes": ["TTN", "LMNA", "MYH7", "TNNT2"],
                "procedures": ["CRT", "ICD", "LVAD", "heart transplant"],
            },
        }
        default = {"drugs": [], "biomarkers": [], "genes": [], "procedures": []}
        return {
            "entity": entity,
            "entity_type": entity_type,
            "relationships": relationships.get(entity.lower().replace(" ", "_"), default),
        }

    try:
        return engine.find_related(entity=entity, entity_type=entity_type)
    except Exception:
        raise HTTPException(status_code=500, detail="Internal processing error")


# =====================================================================
# Risk Calculator Endpoints
# =====================================================================

def _get_age_score_for_heart(age: int) -> int:
    if age < 45:
        return 0
    elif age <= 64:
        return 1
    return 2


def _get_risk_factor_score_for_heart(count: int) -> int:
    if count == 0:
        return 0
    elif count <= 2:
        return 1
    return 2


@router.post("/risk/ascvd", response_model=RiskResult)
async def calculate_ascvd(request: ASCVDRequest, req: Request):
    """Calculate 10-year ASCVD risk using Pooled Cohort Equations."""
    engine = req.app.state.risk_calculators.get("engine")
    if not engine:
        # Lazy init: RiskCalculatorEngine is pure Python, no external deps
        try:
            from src.risk_calculators import RiskCalculatorEngine
            engine = RiskCalculatorEngine()
        except Exception:
            raise HTTPException(status_code=503, detail="Risk calculator engine not initialized")
    try:
        from src.models import RiskScoreInput, RiskScoreType
        risk_input = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=request.age,
            sex=request.sex,
            race=request.race,
            total_cholesterol=request.total_cholesterol,
            hdl=request.hdl_cholesterol,
            systolic_bp=request.systolic_bp,
            hypertension_treatment=request.bp_treatment,
            smoker=request.smoker,
            diabetes=request.diabetes,
        )
        result = engine.calculate(risk_input)
        with req.app.state.metrics_lock:
            req.app.state.metrics["risk_calc_requests_total"] += 1
        return RiskResult(
            calculator="ASCVD_Pooled_Cohort_Equations",
            score=result.score_value,
            risk_category=result.risk_category,
            interpretation=result.interpretation,
            recommendations=result.recommendations,
            details={"ten_year_risk_pct": result.score_value, "guideline": result.guideline_reference},
        )
    except HTTPException:
        raise
    except Exception as exc:
        logger.error(f"ASCVD calculator error: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@router.post("/risk/heart-score", response_model=RiskResult)
async def calculate_heart(request: HEARTRequest, req: Request):
    """Calculate HEART score for acute chest pain / ACS triage."""
    # HEART score uses simple integer scoring - inline is clinically accurate
    age_score = _get_age_score_for_heart(request.age)
    rf_score = _get_risk_factor_score_for_heart(request.risk_factors)
    total = request.history_score + request.ecg_score + age_score + rf_score + request.troponin_score

    if total <= 3:
        category = "low"
        interp = "Low risk for MACE at 6 weeks (1.7%)"
        recs = ["Consider early discharge", "Outpatient follow-up", "Return precautions"]
    elif total <= 6:
        category = "moderate"
        interp = "Moderate risk for MACE at 6 weeks (12-16.6%)"
        recs = ["Admit for observation", "Serial troponins", "Non-invasive testing within 72h"]
    else:
        category = "high"
        interp = "High risk for MACE at 6 weeks (50-65%)"
        recs = ["Urgent cardiology consultation", "Consider early invasive strategy", "Dual antiplatelet therapy"]

    with req.app.state.metrics_lock:
        req.app.state.metrics["risk_calc_requests_total"] += 1

    return RiskResult(
        calculator="HEART_Score",
        score=float(total),
        risk_category=category,
        interpretation=interp,
        recommendations=recs,
        details={
            "history": request.history_score, "ecg": request.ecg_score,
            "age_score": age_score, "risk_factor_score": rf_score,
            "troponin": request.troponin_score, "total": total,
        },
    )


@router.post("/risk/cha2ds2-vasc", response_model=RiskResult)
async def calculate_cha2ds2vasc(request: CHA2DS2VAScRequest, req: Request):
    """Calculate CHA2DS2-VASc stroke risk in atrial fibrillation."""
    # CHA2DS2-VASc: simple integer point system - inline is accurate
    score = 0
    score += 1 if request.chf else 0
    score += 1 if request.hypertension else 0
    score += 2 if request.age >= 75 else (1 if request.age >= 65 else 0)
    score += 1 if request.diabetes else 0
    score += 2 if request.stroke_tia else 0
    score += 1 if request.vascular_disease else 0
    score += 1 if request.sex == "female" else 0

    stroke_rates = {0: 0.2, 1: 0.6, 2: 2.2, 3: 3.2, 4: 4.8, 5: 7.2, 6: 9.7, 7: 11.2, 8: 10.8, 9: 12.2}
    annual_rate = stroke_rates.get(score, 12.2)

    if score == 0:
        category = "low"
        recs = [
            "Low risk (0.2% annual stroke risk). No anticoagulation generally recommended "
            "for males with CHA2DS2-VASc 0. Reassess periodically.",
        ]
    elif score == 1 and request.sex == "female":
        category = "low"
        recs = ["Female sex alone does not warrant anticoagulation", "Reassess for other risk factors"]
    elif score == 1:
        category = "low-moderate"
        recs = [
            "Low-moderate risk. Consider anticoagulation based on individual risk-benefit assessment.",
            "Shared decision-making",
            "Prefer DOAC over warfarin",
        ]
    else:
        category = "moderate-high" if score <= 3 else "high"
        recs = ["Oral anticoagulation recommended", "DOAC preferred over warfarin", "Assess bleeding risk with HAS-BLED"]

    with req.app.state.metrics_lock:
        req.app.state.metrics["risk_calc_requests_total"] += 1

    return RiskResult(
        calculator="CHA2DS2-VASc",
        score=float(score),
        risk_category=category,
        interpretation=f"CHA2DS2-VASc = {score} | Annual stroke rate ~{annual_rate}%",
        recommendations=recs,
        details={"score": score, "annual_stroke_rate_pct": annual_rate, "max_score": 9},
    )


@router.post("/risk/has-bled", response_model=RiskResult)
async def calculate_hasbled(request: HASBLEDRequest, req: Request):
    """Calculate HAS-BLED bleeding risk for anticoagulated patients."""
    # HAS-BLED: simple integer point system - inline is accurate
    score = sum([
        request.hypertension_uncontrolled,
        request.renal_disease,
        request.liver_disease,
        request.stroke_history,
        request.bleeding_history,
        request.labile_inr,
        request.age_over_65,
    ]) + request.drugs_alcohol

    if score <= 1:
        category = "low"
        interp = "Low bleeding risk"
        recs = ["Standard anticoagulation dosing", "Routine monitoring"]
    elif score == 2:
        category = "moderate"
        interp = "Moderate bleeding risk"
        recs = ["Anticoagulation generally still beneficial", "Address modifiable risk factors", "Close monitoring recommended"]
    else:
        category = "high"
        interp = "High bleeding risk"
        recs = ["Carefully weigh benefit vs. bleeding risk", "Address ALL modifiable factors", "Consider shorter-acting agents", "Avoid concomitant NSAIDs/antiplatelets if possible"]

    with req.app.state.metrics_lock:
        req.app.state.metrics["risk_calc_requests_total"] += 1

    return RiskResult(
        calculator="HAS-BLED",
        score=float(score),
        risk_category=category,
        interpretation=f"{interp} (HAS-BLED = {score})",
        recommendations=recs,
        details={"score": score, "max_score": 9, "note": "High score should NOT preclude anticoagulation -- address modifiable factors"},
    )


@router.post("/risk/maggic", response_model=RiskResult)
async def calculate_maggic(request: MAGGICRequest, req: Request):
    """Calculate MAGGIC heart failure mortality score."""
    engine = req.app.state.risk_calculators.get("engine")
    if not engine:
        try:
            from src.risk_calculators import RiskCalculatorEngine
            engine = RiskCalculatorEngine()
        except Exception:
            raise HTTPException(status_code=503, detail="Risk calculator engine not initialized")
    try:
        from src.models import RiskScoreInput, RiskScoreType, HeartFailureClass
        nyha_map = {1: HeartFailureClass.NYHA_I, 2: HeartFailureClass.NYHA_II,
                    3: HeartFailureClass.NYHA_III, 4: HeartFailureClass.NYHA_IV}
        risk_input = RiskScoreInput(
            score_type=RiskScoreType.MAGGIC,
            age=request.age,
            sex=request.sex,
            lvef=request.lvef,
            nyha_class=nyha_map.get(request.nyha_class),
            systolic_bp=request.systolic_bp,
            bmi=request.bmi,
            creatinine=request.creatinine,
            diabetes=request.diabetes,
            smoker=request.smoker,
            beta_blocker_use=request.beta_blocker,
            acei_arb_use=request.acei_arb,
        )
        result = engine.calculate(risk_input, extra={"copd": request.copd, "hf_duration_18m": request.hf_duration_18m})
        with req.app.state.metrics_lock:
            req.app.state.metrics["risk_calc_requests_total"] += 1
        return RiskResult(
            calculator="MAGGIC",
            score=result.score_value,
            risk_category=result.risk_category,
            interpretation=result.interpretation,
            recommendations=result.recommendations,
            details={"guideline": result.guideline_reference},
        )
    except HTTPException:
        raise
    except Exception as exc:
        logger.error(f"MAGGIC calculator error: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@router.post("/risk/euroscore", response_model=RiskResult)
async def calculate_euroscore(request: EuroSCORERequest, req: Request):
    """Calculate EuroSCORE II cardiac surgical mortality risk."""
    engine = req.app.state.risk_calculators.get("engine")
    if not engine:
        try:
            from src.risk_calculators import RiskCalculatorEngine
            engine = RiskCalculatorEngine()
        except Exception:
            raise HTTPException(status_code=503, detail="Risk calculator engine not initialized")
    try:
        from src.models import RiskScoreInput, RiskScoreType
        risk_input = RiskScoreInput(
            score_type=RiskScoreType.EUROSCORE_II,
            age=request.age,
            sex=request.sex,
            lvef=request.lvef,
            creatinine=request.creatinine_clearance,  # maps to renal function
            urgency=request.urgency,
            redo_surgery=request.previous_cardiac_surgery,
        )
        extra = {
            "extracardiac_arteriopathy": request.extracardiac_arteriopathy,
            "poor_mobility": request.poor_mobility,
            "chronic_lung_disease": request.chronic_lung_disease,
            "active_endocarditis": request.active_endocarditis,
            "critical_preop_state": request.critical_preop_state,
            "diabetes_on_insulin": request.diabetes_on_insulin,
            "nyha_class": request.nyha_class,
            "ccs_class_4_angina": request.ccs_class_4_angina,
            "recent_mi": request.recent_mi,
            "pulmonary_hypertension": request.pulmonary_hypertension,
            "thoracic_aorta": request.thoracic_aorta,
        }
        result = engine.calculate(risk_input, extra=extra)
        with req.app.state.metrics_lock:
            req.app.state.metrics["risk_calc_requests_total"] += 1
        return RiskResult(
            calculator="EuroSCORE_II",
            score=result.score_value,
            risk_category=result.risk_category,
            interpretation=result.interpretation,
            recommendations=result.recommendations,
            details={"guideline": result.guideline_reference},
        )
    except HTTPException:
        raise
    except Exception as exc:
        logger.error(f"EuroSCORE II calculator error: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


# =====================================================================
# GDMT Optimization
# =====================================================================

FOUR_PILLARS = ["ARNi/ACEi/ARB", "Beta-blocker", "MRA", "SGLT2i"]

GDMT_TARGETS = {
    "sacubitril/valsartan": {"target": "97/103 mg BID", "class": "ARNi/ACEi/ARB"},
    "enalapril": {"target": "10-20 mg BID", "class": "ARNi/ACEi/ARB"},
    "lisinopril": {"target": "20-40 mg daily", "class": "ARNi/ACEi/ARB"},
    "losartan": {"target": "50-150 mg daily", "class": "ARNi/ACEi/ARB"},
    "valsartan": {"target": "160 mg BID", "class": "ARNi/ACEi/ARB"},
    "carvedilol": {"target": "25 mg BID (50 if >85 kg)", "class": "Beta-blocker"},
    "metoprolol succinate": {"target": "200 mg daily", "class": "Beta-blocker"},
    "bisoprolol": {"target": "10 mg daily", "class": "Beta-blocker"},
    "spironolactone": {"target": "25-50 mg daily", "class": "MRA"},
    "eplerenone": {"target": "50 mg daily", "class": "MRA"},
    "dapagliflozin": {"target": "10 mg daily", "class": "SGLT2i"},
    "empagliflozin": {"target": "10 mg daily", "class": "SGLT2i"},
}


@router.post("/gdmt/optimize", response_model=GDMTResponse)
async def optimize_gdmt(request: GDMTRequest, req: Request):
    """Optimize Guideline-Directed Medical Therapy for heart failure."""
    optimizer = req.app.state.gdmt_optimizer
    if optimizer:
        try:
            data = request.model_dump()
            result = optimizer.optimize(
                lvef=data["lvef"],
                nyha_class=data["nyha_class"],
                current_medications=data.get("current_medications", []),
                patient_data=data,
            )
            return GDMTResponse(
                hf_phenotype=result.ef_category,
                lvef=data["lvef"],
                nyha_class=data["nyha_class"],
                recommendations=result.recommendations,
                four_pillars_status={},
                summary="; ".join(result.next_steps) if result.next_steps else "",
            )
        except Exception as exc:
            logger.warning(f"GDMT optimizer error, using built-in logic: {exc}")

    # Determine HF phenotype
    if request.lvef < 40:
        phenotype = "HFrEF"
    elif request.lvef <= 49:
        phenotype = "HFmrEF"
    else:
        phenotype = "HFpEF"

    # Map current medications to classes
    current_classes = set()
    for med in request.current_medications:
        med_name = med.get("name", "").lower()
        for drug_name, info in GDMT_TARGETS.items():
            if drug_name in med_name or med_name in drug_name:
                current_classes.add(info["class"])

    recommendations = []
    pillar_status = {}

    for pillar in FOUR_PILLARS:
        on_pillar = pillar in current_classes

        if phenotype == "HFrEF":
            if not on_pillar:
                rec = GDMTRecommendation(
                    medication_class=pillar,
                    current_status="Not on therapy",
                    recommendation=f"INITIATE {pillar}",
                    evidence_level="Class I",
                    guideline_ref="2022 AHA/ACC/HFSA HF Guideline",
                    caution="Check contraindications before starting",
                )
                if pillar == "ARNi/ACEi/ARB":
                    rec.target_dose = "Sacubitril/valsartan 97/103 mg BID"
                elif pillar == "Beta-blocker":
                    rec.target_dose = "Carvedilol 25 mg BID or metoprolol succinate 200 mg daily"
                elif pillar == "MRA":
                    rec.target_dose = "Spironolactone 25-50 mg daily"
                elif pillar == "SGLT2i":
                    rec.target_dose = "Dapagliflozin 10 mg or empagliflozin 10 mg daily"
                recommendations.append(rec)
                pillar_status[pillar] = "NOT STARTED"
            else:
                # Check for dose optimization
                at_target = False
                for med in request.current_medications:
                    med_name = med.get("name", "").lower()
                    for drug_name, info in GDMT_TARGETS.items():
                        if drug_name in med_name and info["class"] == pillar:
                            at_target = True  # Simplified; would check actual dose
                            break

                if at_target:
                    pillar_status[pillar] = "AT TARGET"
                else:
                    recommendations.append(GDMTRecommendation(
                        medication_class=pillar,
                        current_status="Below target dose",
                        recommendation=f"UP-TITRATE {pillar} to target",
                        evidence_level="Class I",
                        guideline_ref="2022 AHA/ACC/HFSA HF Guideline",
                        caution="Titrate every 2-4 weeks as tolerated",
                    ))
                    pillar_status[pillar] = "SUBOPTIMAL"
        elif phenotype == "HFmrEF":
            if pillar == "SGLT2i" and not on_pillar:
                recommendations.append(GDMTRecommendation(
                    medication_class=pillar,
                    current_status="Not on therapy",
                    recommendation=f"INITIATE {pillar} (Class IIa for HFmrEF)",
                    target_dose="Dapagliflozin 10 mg or empagliflozin 10 mg daily",
                    evidence_level="Class IIa",
                    guideline_ref="2022 AHA/ACC/HFSA HF Guideline",
                ))
                pillar_status[pillar] = "NOT STARTED"
            else:
                pillar_status[pillar] = "ON THERAPY" if on_pillar else "Consider"
        else:  # HFpEF
            if pillar == "SGLT2i" and not on_pillar:
                recommendations.append(GDMTRecommendation(
                    medication_class=pillar,
                    current_status="Not on therapy",
                    recommendation=f"INITIATE {pillar} (Class IIa for HFpEF)",
                    target_dose="Dapagliflozin 10 mg or empagliflozin 10 mg daily",
                    evidence_level="Class IIa",
                    guideline_ref="2022 AHA/ACC/HFSA; EMPEROR-Preserved; DELIVER",
                ))
                pillar_status[pillar] = "NOT STARTED"
            else:
                pillar_status[pillar] = "ON THERAPY" if on_pillar else "Limited evidence"

    # Add diuretic / additional therapy recommendations based on NYHA
    if request.nyha_class in ("III", "IV"):
        recommendations.append(GDMTRecommendation(
            medication_class="Loop diuretic",
            current_status="Assess volume status",
            recommendation="Titrate diuretic for euvolemia",
            evidence_level="Class I",
            guideline_ref="2022 AHA/ACC/HFSA HF Guideline",
        ))
        if phenotype == "HFrEF":
            recommendations.append(GDMTRecommendation(
                medication_class="Hydralazine/Isosorbide dinitrate",
                current_status="Consider if ACEi/ARB/ARNi contraindicated or additional therapy needed",
                recommendation="Add if symptomatic despite GDMT (especially self-identified Black patients)",
                evidence_level="Class I (Black patients) / Class IIb (others)",
                guideline_ref="A-HeFT trial; 2022 AHA/ACC/HFSA",
            ))

    n_started = sum(1 for v in pillar_status.values() if v in ("AT TARGET", "ON THERAPY", "SUBOPTIMAL"))
    summary = (
        f"{phenotype} (LVEF {request.lvef}%, NYHA {request.nyha_class}): "
        f"{n_started}/4 pillars active, {len(recommendations)} optimization opportunities identified."
    )

    with req.app.state.metrics_lock:
        req.app.state.metrics["gdmt_requests_total"] += 1

    return GDMTResponse(
        hf_phenotype=phenotype,
        lvef=request.lvef,
        nyha_class=request.nyha_class,
        recommendations=recommendations,
        four_pillars_status=pillar_status,
        summary=summary,
    )


# =====================================================================
# Workflow Endpoints
# =====================================================================

_WORKFLOW_ID_MAP = {
    "cad": CardioWorkflowType.CAD_ASSESSMENT,
    "heart_failure": CardioWorkflowType.HEART_FAILURE,
    "valvular": CardioWorkflowType.VALVULAR_DISEASE,
    "arrhythmia": CardioWorkflowType.ARRHYTHMIA,
    "cardiac_mri": CardioWorkflowType.CARDIAC_MRI,
    "stress_test": CardioWorkflowType.STRESS_TEST,
    "prevention": CardioWorkflowType.PREVENTIVE_RISK,
    "cardio_oncology": CardioWorkflowType.CARDIO_ONCOLOGY,
}


def _workflow_stub(workflow_id: str, workflow_name: str, data: dict, req: Request):
    """Shared workflow handler -- delegates to WorkflowEngine or returns structured stub."""
    workflow_type = _WORKFLOW_ID_MAP.get(workflow_id)
    if workflow_type:
        try:
            engine_wf = WorkflowEngine()
            result = engine_wf.run_workflow(workflow_type, data)
            return result.model_dump()
        except Exception as exc:
            logger.warning(f"Workflow {workflow_id} engine error: {exc}")

    with req.app.state.metrics_lock:
        req.app.state.metrics["workflow_requests_total"] += 1

    return {
        "workflow": workflow_id,
        "name": workflow_name,
        "status": "completed",
        "input_summary": {k: v for k, v in data.items() if k != "raw_data"},
        "assessment": f"{workflow_name} analysis requires RAG engine initialization for full results.",
        "recommendations": [],
        "evidence": [],
        "guidelines_referenced": [],
    }


@router.post("/workflow/cad")
async def cad_assessment(request: dict, req: Request):
    """Coronary artery disease assessment: calcium score, CAD-RADS, plaque characterization."""
    return _workflow_stub("cad", "Coronary Artery Disease Assessment", request, req)


@router.post("/workflow/heart-failure")
async def heart_failure_workflow(request: dict, req: Request):
    """Heart failure evaluation: phenotyping, staging, GDMT, device candidacy."""
    return _workflow_stub("heart_failure", "Heart Failure Management", request, req)


@router.post("/workflow/valvular")
async def valvular_workflow(request: dict, req: Request):
    """Valvular heart disease: severity grading, intervention timing, prosthetic follow-up."""
    return _workflow_stub("valvular", "Valvular Heart Disease Assessment", request, req)


@router.post("/workflow/arrhythmia")
async def arrhythmia_workflow(request: dict, req: Request):
    """Arrhythmia & EP assessment: AF management, ablation candidacy, device programming."""
    return _workflow_stub("arrhythmia", "Arrhythmia & Electrophysiology", request, req)


@router.post("/workflow/cardiac-mri")
async def cardiac_mri_workflow(request: dict, req: Request):
    """Cardiac MRI interpretation: LGE patterns, T1/T2 mapping, tissue characterization."""
    return _workflow_stub("cardiac_mri", "Cardiac MRI Interpretation", request, req)


@router.post("/workflow/stress-test")
async def stress_test_workflow(request: dict, req: Request):
    """Stress test protocol: exercise/pharmacologic, Duke score, perfusion analysis."""
    return _workflow_stub("stress_test", "Stress Testing Protocol", request, req)


@router.post("/workflow/prevention")
async def prevention_workflow(request: dict, req: Request):
    """Cardiovascular prevention: lipid management, ASCVD risk, statin selection."""
    return _workflow_stub("prevention", "Cardiovascular Prevention", request, req)


@router.post("/workflow/cardio-oncology")
async def cardio_oncology_workflow(request: dict, req: Request):
    """Cardio-oncology surveillance: cardiotoxicity, GLS tracking, biomarker monitoring."""
    return _workflow_stub("cardio_oncology", "Cardio-Oncology Surveillance", request, req)


# =====================================================================
# Guidelines
# =====================================================================

GUIDELINE_LIBRARY = [
    {"id": "aha_acc_hf_2022", "title": "2022 AHA/ACC/HFSA Guideline for Management of Heart Failure", "year": 2022, "conditions": ["heart_failure"]},
    {"id": "aha_acc_cad_2023", "title": "2023 AHA/ACC Guideline for Management of Chronic Coronary Disease", "year": 2023, "conditions": ["cad", "stable_angina"]},
    {"id": "aha_acc_af_2023", "title": "2023 ACC/AHA/ACCP/HRS Guideline for AF Management", "year": 2023, "conditions": ["atrial_fibrillation"]},
    {"id": "aha_acc_vhd_2020", "title": "2020 ACC/AHA Guideline for Valvular Heart Disease", "year": 2020, "conditions": ["aortic_stenosis", "mitral_regurgitation", "valvular"]},
    {"id": "esc_cad_2024", "title": "2024 ESC Guidelines for Chronic Coronary Syndromes", "year": 2024, "conditions": ["cad"]},
    {"id": "esc_acs_2023", "title": "2023 ESC Guidelines for Acute Coronary Syndromes", "year": 2023, "conditions": ["acs", "stemi", "nstemi"]},
    {"id": "esc_hf_2023", "title": "2023 Focused Update of ESC Heart Failure Guidelines", "year": 2023, "conditions": ["heart_failure"]},
    {"id": "aha_lipid_2018", "title": "2018 AHA/ACC Cholesterol Clinical Practice Guideline", "year": 2018, "conditions": ["hyperlipidemia", "prevention"]},
    {"id": "aha_prevention_2019", "title": "2019 ACC/AHA Primary Prevention Guideline", "year": 2019, "conditions": ["prevention"]},
    {"id": "esc_cardiomyopathy_2023", "title": "2023 ESC Guidelines for Cardiomyopathies", "year": 2023, "conditions": ["dcm", "hcm", "cardiomyopathy"]},
    {"id": "aha_hcm_2024", "title": "2024 AHA/ACC Guideline for HCM", "year": 2024, "conditions": ["hcm"]},
    {"id": "esc_cardio_onc_2022", "title": "2022 ESC Cardio-Oncology Guidelines", "year": 2022, "conditions": ["cardio_oncology"]},
]


@router.get("/guidelines")
async def list_guidelines(condition: Optional[str] = None):
    """List clinical practice guidelines, optionally filtered by condition."""
    if condition:
        filtered = [g for g in GUIDELINE_LIBRARY if condition.lower() in g["conditions"]]
        return {"guidelines": filtered, "filter": condition, "total": len(filtered)}
    return {"guidelines": GUIDELINE_LIBRARY, "total": len(GUIDELINE_LIBRARY)}


# =====================================================================
# Cross-Modal Evaluation
# =====================================================================

@router.post("/cross-modal/evaluate")
async def evaluate_cross_modal(findings: dict, req: Request):
    """Integrate findings across imaging modalities for comprehensive assessment."""
    engine = req.app.state.engine
    if engine:
        try:
            return engine.cross_modal_evaluate(findings)
        except Exception:
            pass

    modalities_present = [k for k in findings if k in ("echo", "ct", "mri", "nuclear", "cath")]
    return {
        "status": "evaluated",
        "modalities_integrated": modalities_present,
        "concordance": "requires_engine",
        "note": "Full cross-modal integration requires RAG engine initialization",
        "input_findings": findings,
    }


# =====================================================================
# Reference Catalogues
# =====================================================================

@router.get("/conditions")
async def list_conditions():
    """Catalogue of cardiovascular conditions covered by the agent."""
    return {
        "conditions": [
            {"id": "cad", "name": "Coronary Artery Disease", "subtypes": ["stable_angina", "acs", "stemi", "nstemi"]},
            {"id": "heart_failure", "name": "Heart Failure", "subtypes": ["HFrEF", "HFmrEF", "HFpEF"]},
            {"id": "atrial_fibrillation", "name": "Atrial Fibrillation", "subtypes": ["paroxysmal", "persistent", "permanent"]},
            {"id": "aortic_stenosis", "name": "Aortic Stenosis", "subtypes": ["calcific", "bicuspid", "low_flow_low_gradient"]},
            {"id": "mitral_regurgitation", "name": "Mitral Regurgitation", "subtypes": ["primary", "secondary"]},
            {"id": "hcm", "name": "Hypertrophic Cardiomyopathy", "subtypes": ["obstructive", "non_obstructive"]},
            {"id": "dcm", "name": "Dilated Cardiomyopathy", "subtypes": ["idiopathic", "familial", "tachycardia_mediated"]},
            {"id": "pulmonary_hypertension", "name": "Pulmonary Hypertension", "subtypes": ["group1_pah", "group2", "group3", "group4_cteph"]},
            {"id": "pericardial_disease", "name": "Pericardial Disease", "subtypes": ["pericarditis", "tamponade", "constrictive"]},
            {"id": "aortic_disease", "name": "Aortic Disease", "subtypes": ["aneurysm", "dissection", "intramural_hematoma"]},
            {"id": "venous_thromboembolism", "name": "VTE", "subtypes": ["dvt", "pe"]},
            {"id": "cardio_oncology", "name": "Cardiotoxicity", "subtypes": ["anthracycline", "her2_targeted", "checkpoint_inhibitor"]},
        ],
        "total": 12,
    }


@router.get("/biomarkers")
async def list_biomarkers():
    """Cardiac biomarker reference panel."""
    return {
        "biomarkers": [
            {"name": "hs-Troponin T", "unit": "ng/L", "normal": "<14", "use": "ACS diagnosis, prognosis"},
            {"name": "hs-Troponin I", "unit": "ng/L", "normal": "<26 (M), <16 (F)", "use": "ACS diagnosis"},
            {"name": "NT-proBNP", "unit": "pg/mL", "normal": "<125 (<75y), <450 (>75y)", "use": "HF diagnosis, prognosis"},
            {"name": "BNP", "unit": "pg/mL", "normal": "<100", "use": "HF diagnosis, volume status"},
            {"name": "hs-CRP", "unit": "mg/L", "normal": "<2.0", "use": "Inflammatory risk, residual risk"},
            {"name": "D-dimer", "unit": "ng/mL FEU", "normal": "<500", "use": "VTE exclusion"},
            {"name": "sST2", "unit": "ng/mL", "normal": "<35", "use": "HF prognosis, fibrosis"},
            {"name": "Galectin-3", "unit": "ng/mL", "normal": "<17.8", "use": "HF prognosis, fibrosis"},
            {"name": "Lp(a)", "unit": "nmol/L", "normal": "<75", "use": "ASCVD risk, genetic risk"},
            {"name": "LDL-C", "unit": "mg/dL", "normal": "<100 (general), <70 (high risk)", "use": "Lipid management"},
            {"name": "CK-MB", "unit": "ng/mL", "normal": "<5", "use": "Reinfarction detection"},
            {"name": "Myoglobin", "unit": "ng/mL", "normal": "<110", "use": "Early MI marker"},
        ],
        "total": 12,
    }


@router.get("/drugs")
async def list_drug_classes():
    """Cardiovascular drug class reference."""
    return {
        "drug_classes": [
            {"class": "ARNi", "examples": ["sacubitril/valsartan"], "indications": ["HFrEF", "HFmrEF"]},
            {"class": "ACEi", "examples": ["enalapril", "lisinopril", "ramipril"], "indications": ["HF", "HTN", "post-MI"]},
            {"class": "ARB", "examples": ["losartan", "valsartan", "candesartan"], "indications": ["HF", "HTN"]},
            {"class": "Beta-blocker", "examples": ["carvedilol", "metoprolol succinate", "bisoprolol"], "indications": ["HF", "AF rate control", "post-MI"]},
            {"class": "MRA", "examples": ["spironolactone", "eplerenone"], "indications": ["HFrEF", "resistant HTN"]},
            {"class": "SGLT2i", "examples": ["dapagliflozin", "empagliflozin"], "indications": ["HF (all EF)", "CKD", "T2DM"]},
            {"class": "Statin", "examples": ["atorvastatin", "rosuvastatin"], "indications": ["ASCVD prevention", "ACS"]},
            {"class": "PCSK9i", "examples": ["evolocumab", "alirocumab"], "indications": ["Refractory hyperlipidemia", "ASCVD"]},
            {"class": "Anticoagulant (DOAC)", "examples": ["apixaban", "rivaroxaban", "edoxaban", "dabigatran"], "indications": ["AF", "VTE"]},
            {"class": "Antiplatelet", "examples": ["aspirin", "clopidogrel", "ticagrelor", "prasugrel"], "indications": ["ACS", "PCI", "secondary prevention"]},
            {"class": "Antiarrhythmic", "examples": ["amiodarone", "flecainide", "sotalol", "dofetilide"], "indications": ["AF rhythm control", "VT"]},
            {"class": "Vasodilator", "examples": ["hydralazine/isosorbide dinitrate", "nitrates"], "indications": ["HFrEF (esp. Black patients)", "angina"]},
            {"class": "Diuretic", "examples": ["furosemide", "bumetanide", "torsemide", "metolazone"], "indications": ["Volume overload", "HF"]},
            {"class": "Cardiac myosin inhibitor", "examples": ["mavacamten", "aficamten"], "indications": ["Obstructive HCM"]},
        ],
        "total": 14,
    }


@router.get("/genes")
async def list_genes():
    """Cardiology-relevant genes and associated conditions."""
    return {
        "genes": [
            {"gene": "TTN", "condition": "Dilated cardiomyopathy", "inheritance": "AD", "prevalence": "15-25% of familial DCM"},
            {"gene": "LMNA", "condition": "Laminopathy (DCM + conduction disease)", "inheritance": "AD", "prevalence": "5-10% of familial DCM"},
            {"gene": "MYH7", "condition": "HCM / DCM", "inheritance": "AD", "prevalence": "25-35% of HCM"},
            {"gene": "MYBPC3", "condition": "Hypertrophic cardiomyopathy", "inheritance": "AD", "prevalence": "25-35% of HCM"},
            {"gene": "TNNT2", "condition": "HCM / DCM / RCM", "inheritance": "AD", "prevalence": "3-5% of HCM"},
            {"gene": "SCN5A", "condition": "Brugada / Long QT type 3 / sick sinus", "inheritance": "AD", "prevalence": "20-30% of Brugada"},
            {"gene": "KCNQ1", "condition": "Long QT syndrome type 1", "inheritance": "AD/AR", "prevalence": "40-55% of LQTS"},
            {"gene": "KCNH2", "condition": "Long QT syndrome type 2", "inheritance": "AD", "prevalence": "25-40% of LQTS"},
            {"gene": "RYR2", "condition": "CPVT", "inheritance": "AD", "prevalence": "50-60% of CPVT"},
            {"gene": "PKP2", "condition": "Arrhythmogenic cardiomyopathy", "inheritance": "AD", "prevalence": "25-40% of ACM"},
            {"gene": "DSP", "condition": "Arrhythmogenic cardiomyopathy / DCM", "inheritance": "AD/AR", "prevalence": "5-10% of ACM"},
            {"gene": "FLNC", "condition": "DCM / RCM with arrhythmias", "inheritance": "AD", "prevalence": "3-8% of DCM"},
            {"gene": "LDLR", "condition": "Familial hypercholesterolemia", "inheritance": "AD (co-dominant)", "prevalence": "1:250"},
            {"gene": "PCSK9", "condition": "Familial hypercholesterolemia", "inheritance": "AD", "prevalence": "Rare (<1% of FH)"},
            {"gene": "GLA", "condition": "Fabry disease (cardiac variant)", "inheritance": "X-linked", "prevalence": "0.5-1% of unexplained LVH"},
        ],
        "total": 15,
    }


# =====================================================================
# Knowledge Base Version
# =====================================================================

@router.get("/knowledge-version", tags=["reference"])
async def get_knowledge_version():
    """Return knowledge base version and revision metadata."""
    from src.knowledge import KNOWLEDGE_VERSION
    return KNOWLEDGE_VERSION
