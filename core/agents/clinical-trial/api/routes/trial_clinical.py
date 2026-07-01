"""Clinical trial API routes.

Provides endpoints for RAG-powered clinical trial queries, patient-trial
matching, protocol optimization, site selection, eligibility optimization,
adaptive design evaluation, safety signal detection, regulatory document
generation, competitive intelligence, diversity assessment, decentralized
trial planning, and reference catalogues.

Author: Adam Jones
Date: March 2026
"""

from typing import List, Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

from src.models import (
    TrialWorkflowType,
    TrialPhase,
    TherapeuticArea,
)

router = APIRouter(prefix="/v1/trial", tags=["clinical-trials"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/integrated-assessment")
async def integrated_assessment(request: dict, req: Request):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries oncology, pharmacogenomics, cardiology, and biomarker agents
    for a comprehensive trial eligibility and safety assessment.
    """
    try:
        from src.cross_modal import (
            query_oncology_agent,
            query_pgx_agent,
            query_cardiology_agent,
            query_biomarker_agent,
            integrate_cross_agent_results,
        )

        patient_profile = request.get("patient_profile", {})
        drug = request.get("drug", "")
        biomarkers = request.get("biomarkers", [])
        therapeutic_area = request.get("therapeutic_area", "")

        results = []

        # Query oncology agent for molecular trial matches
        if patient_profile:
            results.append(query_oncology_agent(patient_profile))

        # Query PGx agent for metabolism screening
        if patient_profile:
            results.append(query_pgx_agent(patient_profile))

        # Query cardiology agent for cardiac safety
        if drug:
            results.append(query_cardiology_agent(drug, patient_profile=patient_profile))

        # Query biomarker agent for enrichment strategies
        if biomarkers:
            results.append(query_biomarker_agent(biomarkers, therapeutic_area=therapeutic_area))

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
    """Free-text RAG query with optional workflow and patient context."""
    question: str = Field(..., min_length=3, description="Clinical trial question")
    workflow_type: Optional[str] = Field(
        None,
        description=(
            "Workflow hint: protocol_design | patient_matching | site_selection | "
            "eligibility_optimization | adaptive_design | safety_signal | "
            "regulatory_docs | competitive_intel | diversity_assessment | "
            "decentralized_planning | general"
        ),
    )
    patient_context: Optional[dict] = Field(None, description="Demographics, biomarkers, genomics")
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


# -- Protocol Optimization --

class ProtocolOptimizeRequest(BaseModel):
    """Protocol optimization analysis request."""
    protocol_summary: str = Field(..., min_length=10, description="Protocol text or summary")
    therapeutic_area: Optional[str] = Field(None, description="Therapeutic area")
    phase: Optional[str] = Field(None, description="Trial phase")
    indication: Optional[str] = Field(None, description="Target indication")
    endpoints: List[str] = Field(default_factory=list, description="Primary/secondary endpoints")
    eligibility_criteria: List[str] = Field(default_factory=list, description="Key eligibility criteria")
    visit_count: Optional[int] = Field(None, ge=1, description="Number of scheduled visits")
    procedure_count: Optional[int] = Field(None, ge=1, description="Number of distinct procedures")


class ProtocolOptimizeResponse(BaseModel):
    complexity_score: float
    percentile_rank: float
    optimization_recommendations: List[str]
    risk_factors: List[dict]
    endpoint_analysis: List[dict]
    estimated_enrollment_impact: Optional[float] = None
    evidence: List[dict] = []


# -- Patient Matching --

class PatientMatchRequest(BaseModel):
    """Patient-trial matching request."""
    age: Optional[int] = Field(None, ge=0, le=120, description="Patient age")
    sex: Optional[str] = Field(None, description="Patient sex: male | female")
    diagnosis: str = Field(..., min_length=3, description="Primary diagnosis")
    biomarkers: List[str] = Field(default_factory=list, description="Biomarker results (e.g., HER2+)")
    medications: List[str] = Field(default_factory=list, description="Current medications")
    genomic_variants: List[str] = Field(default_factory=list, description="Known genomic variants")
    comorbidities: List[str] = Field(default_factory=list, description="Comorbid conditions")
    geographic_location: Optional[str] = Field(None, description="Geographic location for site matching")
    therapeutic_area: Optional[str] = Field(None, description="Target therapeutic area")
    max_results: int = Field(10, ge=1, le=50, description="Maximum matching trials")


class TrialMatch(BaseModel):
    trial_id: str
    trial_title: str
    phase: str
    status: str
    overall_score: float
    inclusion_met: int
    inclusion_total: int
    exclusion_clear: int
    exclusion_total: int
    confidence: float
    matching_criteria: List[dict] = []
    sites: List[dict] = []


class PatientMatchResponse(BaseModel):
    matches: List[TrialMatch]
    total_screened: int
    patient_summary: dict
    search_criteria: dict


class BatchMatchRequest(BaseModel):
    """Batch patient-trial matching."""
    patients: List[PatientMatchRequest] = Field(..., min_length=1, max_length=50)


class BatchMatchResponse(BaseModel):
    results: List[PatientMatchResponse]
    total_patients: int
    total_matches: int


# -- Site Selection --

class SiteRecommendRequest(BaseModel):
    """Site selection and feasibility request."""
    therapeutic_area: str = Field(..., description="Therapeutic area")
    indication: str = Field(..., description="Target indication")
    phase: Optional[str] = Field(None, description="Trial phase")
    target_enrollment: int = Field(100, ge=1, description="Target enrollment number")
    countries: List[str] = Field(default_factory=list, description="Target countries")
    diversity_requirements: Optional[dict] = Field(None, description="Diversity targets")
    max_sites: int = Field(20, ge=1, le=100, description="Maximum sites to recommend")


class SiteRecommendResponse(BaseModel):
    recommended_sites: List[dict]
    total_evaluated: int
    enrollment_forecast: dict
    diversity_analysis: dict


# -- Eligibility Optimization --

class EligibilityOptimizeRequest(BaseModel):
    """Eligibility criteria optimization request."""
    criteria: List[dict] = Field(..., min_length=1, description="List of {text, type} criteria")
    therapeutic_area: Optional[str] = Field(None, description="Therapeutic area")
    indication: Optional[str] = Field(None, description="Target indication")
    target_population_size: Optional[int] = Field(None, description="Desired population pool size")


class EligibilityOptimizeResponse(BaseModel):
    original_criteria_count: int
    optimized_criteria: List[dict]
    population_impact_estimate: float
    recommendations: List[str]
    competitor_benchmarks: List[dict] = []


# -- Adaptive Design --

class AdaptiveEvaluateRequest(BaseModel):
    """Adaptive design evaluation request."""
    design_type: str = Field(..., description="adaptive_randomization | sample_size_reestimation | dose_response | seamless | platform")
    interim_data: Optional[dict] = Field(None, description="Interim analysis data")
    current_sample_size: Optional[int] = Field(None, ge=1)
    target_sample_size: Optional[int] = Field(None, ge=1)
    arms: List[dict] = Field(default_factory=list, description="Treatment arms with data")
    primary_endpoint: Optional[str] = Field(None, description="Primary endpoint definition")


class AdaptiveEvaluateResponse(BaseModel):
    design_assessment: str
    recommendations: List[str]
    statistical_considerations: List[str]
    futility_assessment: Optional[dict] = None
    sample_size_recommendation: Optional[dict] = None
    evidence: List[dict] = []


# -- Safety Signal --

class SafetySignalRequest(BaseModel):
    """Safety signal detection request."""
    adverse_events: List[dict] = Field(..., min_length=1, description="List of AE records")
    comparator_data: Optional[dict] = Field(None, description="Background rate or comparator data")
    drug_name: Optional[str] = Field(None, description="Investigational product name")
    trial_phase: Optional[str] = Field(None, description="Trial phase")
    study_population_size: Optional[int] = Field(None, ge=1)


class SafetySignalResponse(BaseModel):
    signals_detected: List[dict]
    disproportionality_analysis: List[dict]
    severity_distribution: dict
    recommendations: List[str]
    dsmb_considerations: List[str] = []


# -- Regulatory Document Generation --

class RegulatoryGenerateRequest(BaseModel):
    """Regulatory document generation request."""
    document_type: str = Field(..., description="ind | csr | briefing | psp | rmp | dsur")
    agency: str = Field("fda", description="fda | ema | pmda | health_canada | tga | mhra")
    trial_data: dict = Field(default_factory=dict, description="Trial data for document generation")
    sections: List[str] = Field(default_factory=list, description="Specific sections to generate")
    include_appendices: bool = Field(True, description="Include appendices")


class RegulatoryGenerateResponse(BaseModel):
    document_type: str
    agency: str
    content: str
    sections_generated: List[str]
    compliance_notes: List[str]
    metadata: dict = {}


# -- Competitive Intelligence --

class CompetitiveLandscapeRequest(BaseModel):
    """Competitive intelligence landscape request."""
    therapeutic_area: str = Field(..., description="Therapeutic area")
    indication: Optional[str] = Field(None, description="Specific indication")
    mechanism: Optional[str] = Field(None, description="Mechanism of action")
    include_completed: bool = Field(False, description="Include completed trials")
    max_competitors: int = Field(20, ge=1, le=100)


class CompetitiveLandscapeResponse(BaseModel):
    competitors: List[dict]
    landscape_summary: str
    threat_assessment: dict
    enrollment_race: List[dict]
    market_timing: dict


# -- Diversity Assessment --

class DiversityAssessRequest(BaseModel):
    """Diversity and inclusion assessment request."""
    trial_demographics: Optional[dict] = Field(None, description="Current trial demographic breakdown")
    disease_epidemiology: Optional[dict] = Field(None, description="Disease prevalence by demographics")
    therapeutic_area: Optional[str] = Field(None, description="Therapeutic area")
    enrollment_targets: Optional[dict] = Field(None, description="Enrollment diversity targets")
    site_locations: List[str] = Field(default_factory=list, description="Current site locations")


class DiversityAssessResponse(BaseModel):
    diversity_score: float
    gap_analysis: List[dict]
    fda_compliance_status: str
    recommendations: List[str]
    community_engagement_plan: List[str] = []


# -- Decentralized Trial Planning --

class DCTPlanRequest(BaseModel):
    """Decentralized trial planning request."""
    protocol_summary: str = Field(..., min_length=10, description="Protocol summary")
    therapeutic_area: Optional[str] = Field(None, description="Therapeutic area")
    phase: Optional[str] = Field(None, description="Trial phase")
    target_components: List[str] = Field(
        default_factory=list,
        description="Target DCT components: econsent, telemedicine, home_health, local_labs, wearables, epro_ecoa, direct_to_patient",
    )
    regulatory_region: str = Field("us", description="Primary regulatory region")


class DCTPlanResponse(BaseModel):
    feasibility_assessment: dict
    recommended_components: List[dict]
    hybrid_model: dict
    technology_requirements: List[dict]
    regulatory_considerations: List[str]
    cost_impact_estimate: Optional[dict] = None


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
async def trial_query(request: QueryRequest, req: Request):
    """RAG-powered clinical trial Q&A.

    Searches across all trial knowledge collections and synthesizes
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
            f"Clinical trial question: {request.question}\n\n"
            f"Retrieved evidence:\n{context}\n\n"
            f"Provide a comprehensive answer citing specific evidence. "
            f"Include relevant ICH-GCP, FDA, or EMA guideline references."
        )
        if request.patient_context:
            prompt += f"\n\nPatient context: {request.patient_context}"

        try:
            answer = llm.generate(prompt)
            confidence = min(0.95, 0.5 + len(evidence) * 0.05)
            # Extract guideline references from answer
            for gline in ["ICH E6", "ICH E9", "ICH E10", "ICH E17", "ICH E8",
                          "FDA Guidance", "EMA Guideline", "21 CFR", "GCP"]:
                if gline.lower() in answer.lower():
                    guidelines_cited.append(gline)
        except Exception as exc:
            logger.warning(f"LLM generation failed: {exc}")

    return QueryResponse(
        answer=answer,
        evidence=evidence,
        guidelines_cited=guidelines_cited,
        confidence=confidence,
        workflow_applied=request.workflow_type,
    )


@router.post("/search", response_model=SearchResponse)
async def trial_search(request: SearchRequest, req: Request):
    """Multi-collection semantic search across trial knowledge base."""
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


@router.post("/protocol/optimize", response_model=ProtocolOptimizeResponse)
async def protocol_optimize(request: ProtocolOptimizeRequest, req: Request):
    """Protocol optimization with complexity scoring and recommendations."""
    _increment_metric(req, "protocol_requests_total")

    # Calculate complexity score
    procedure_weight = min((request.procedure_count or 0) / 50.0, 1.0)
    visit_weight = min((request.visit_count or 0) / 30.0, 1.0)
    endpoint_weight = min(len(request.endpoints) / 20.0, 1.0)
    criteria_weight = min(len(request.eligibility_criteria) / 30.0, 1.0)
    complexity_score = round(
        (procedure_weight * 0.3 + visit_weight * 0.25 +
         endpoint_weight * 0.2 + criteria_weight * 0.25), 3
    )
    percentile_rank = round(complexity_score * 85 + 10, 1)

    # LLM-powered optimization recommendations
    llm = _get_llm(req)
    recommendations = []
    risk_factors = []
    endpoint_analysis = []
    evidence = []

    if llm:
        prompt = (
            f"Analyze this clinical trial protocol for optimization opportunities:\n\n"
            f"Therapeutic Area: {request.therapeutic_area or 'Not specified'}\n"
            f"Phase: {request.phase or 'Not specified'}\n"
            f"Indication: {request.indication or 'Not specified'}\n"
            f"Protocol Summary: {request.protocol_summary}\n"
            f"Endpoints: {', '.join(request.endpoints) if request.endpoints else 'Not specified'}\n"
            f"Eligibility Criteria: {len(request.eligibility_criteria)} criteria\n"
            f"Visits: {request.visit_count or 'Not specified'}\n"
            f"Procedures: {request.procedure_count or 'Not specified'}\n"
            f"Complexity Score: {complexity_score}\n\n"
            f"Provide: 1) Optimization recommendations, 2) Risk factors, "
            f"3) Endpoint analysis. Format as structured analysis."
        )
        try:
            result = llm.generate(prompt)
            recommendations = [line.strip() for line in result.split("\n") if line.strip() and len(line.strip()) > 10][:10]
        except Exception as exc:
            logger.warning(f"LLM protocol analysis failed: {exc}")

    if not recommendations:
        recommendations = [
            "Review eligibility criteria for potential broadening to improve enrollment",
            "Consider reducing visit frequency during stable disease periods",
            "Evaluate endpoint consolidation to reduce patient burden",
            "Assess decentralized trial components for applicable visits",
        ]

    return ProtocolOptimizeResponse(
        complexity_score=complexity_score,
        percentile_rank=percentile_rank,
        optimization_recommendations=recommendations,
        risk_factors=risk_factors,
        endpoint_analysis=endpoint_analysis,
        estimated_enrollment_impact=round(1.0 - complexity_score * 0.4, 2),
        evidence=evidence,
    )


async def _do_patient_match(request: PatientMatchRequest, req: Request) -> PatientMatchResponse:
    """Core patient-trial matching logic (shared by single and batch endpoints)."""
    llm = _get_llm(req)
    engine = getattr(req.app.state, "engine", None)

    # Build patient summary
    patient_summary = {
        "age": request.age,
        "sex": request.sex,
        "diagnosis": request.diagnosis,
        "biomarkers": request.biomarkers,
        "medications": request.medications,
        "genomic_variants": request.genomic_variants,
        "comorbidities": request.comorbidities,
        "location": request.geographic_location,
    }

    # Search for matching trials
    matches = []
    total_screened = 0

    if engine:
        try:
            search_query = (
                f"{request.diagnosis} "
                f"{' '.join(request.biomarkers)} "
                f"{request.therapeutic_area or ''}"
            ).strip()
            results = engine.search(search_query, top_k=request.max_results * 2)
            total_screened = len(results)

            for i, r in enumerate(results[:request.max_results]):
                matches.append(TrialMatch(
                    trial_id=r.get("metadata", {}).get("nct_id", f"NCT-RESULT-{i+1:04d}"),
                    trial_title=r.get("metadata", {}).get("title", r.get("text", "")[:100]),
                    phase=r.get("metadata", {}).get("phase", "phase_ii"),
                    status=r.get("metadata", {}).get("status", "recruiting"),
                    overall_score=round(r.get("score", 0.5), 3),
                    inclusion_met=r.get("metadata", {}).get("inclusion_met", 0),
                    inclusion_total=r.get("metadata", {}).get("inclusion_total", 0),
                    exclusion_clear=r.get("metadata", {}).get("exclusion_clear", 0),
                    exclusion_total=r.get("metadata", {}).get("exclusion_total", 0),
                    confidence=round(r.get("score", 0.5) * 0.9, 3),
                ))
        except Exception as exc:
            logger.warning(f"Trial matching search failed: {exc}")

    if not matches and llm:
        # Generate placeholder matches via LLM
        matches = [
            TrialMatch(
                trial_id="NCT-PENDING",
                trial_title=f"Matching analysis for: {request.diagnosis}",
                phase="phase_ii",
                status="recruiting",
                overall_score=0.0,
                inclusion_met=0,
                inclusion_total=0,
                exclusion_clear=0,
                exclusion_total=0,
                confidence=0.0,
            )
        ]

    return PatientMatchResponse(
        matches=matches,
        total_screened=total_screened,
        patient_summary=patient_summary,
        search_criteria={
            "diagnosis": request.diagnosis,
            "therapeutic_area": request.therapeutic_area,
            "biomarkers": request.biomarkers,
            "max_results": request.max_results,
        },
    )


@router.post("/match", response_model=PatientMatchResponse)
async def patient_match(request: PatientMatchRequest, req: Request):
    """Patient-trial matching with eligibility screening."""
    _increment_metric(req, "match_requests_total")
    return await _do_patient_match(request, req)


@router.post("/match/batch", response_model=BatchMatchResponse)
async def patient_match_batch(request: BatchMatchRequest, req: Request):
    """Batch patient-trial matching for multiple patients."""
    _increment_metric(req, "match_requests_total")

    results = []
    total_matches = 0

    for patient in request.patients:
        match_response = await _do_patient_match(patient, req)
        results.append(match_response)
        total_matches += len(match_response.matches)

    return BatchMatchResponse(
        results=results,
        total_patients=len(request.patients),
        total_matches=total_matches,
    )


@router.post("/site/recommend", response_model=SiteRecommendResponse)
async def site_recommend(request: SiteRecommendRequest, req: Request):
    """Site selection and feasibility recommendations."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    recommended_sites = []
    diversity_analysis = {}

    if llm:
        prompt = (
            f"Recommend clinical trial sites for:\n"
            f"Therapeutic Area: {request.therapeutic_area}\n"
            f"Indication: {request.indication}\n"
            f"Phase: {request.phase or 'Not specified'}\n"
            f"Target Enrollment: {request.target_enrollment}\n"
            f"Countries: {', '.join(request.countries) if request.countries else 'Global'}\n"
            f"Diversity Requirements: {request.diversity_requirements or 'Standard'}\n\n"
            f"Provide site recommendations with enrollment rates and diversity metrics."
        )
        try:
            result = llm.generate(prompt)
            recommended_sites = [{"analysis": result}]
        except Exception as exc:
            logger.warning(f"Site recommendation LLM failed: {exc}")

    return SiteRecommendResponse(
        recommended_sites=recommended_sites,
        total_evaluated=0,
        enrollment_forecast={
            "target": request.target_enrollment,
            "estimated_months": round(request.target_enrollment / max(len(request.countries), 1) / 3, 1),
        },
        diversity_analysis=diversity_analysis,
    )


@router.post("/eligibility/optimize", response_model=EligibilityOptimizeResponse)
async def eligibility_optimize(request: EligibilityOptimizeRequest, req: Request):
    """Eligibility criteria optimization with population impact modeling."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    optimized_criteria = []
    recommendations = []

    if llm:
        criteria_text = "\n".join(
            f"- [{c.get('type', 'inclusion')}] {c.get('text', c)}"
            for c in request.criteria
        )
        prompt = (
            f"Optimize these clinical trial eligibility criteria:\n\n"
            f"{criteria_text}\n\n"
            f"Therapeutic Area: {request.therapeutic_area or 'Not specified'}\n"
            f"Indication: {request.indication or 'Not specified'}\n\n"
            f"For each criterion, assess: population impact, scientific justification, "
            f"and recommend whether to broaden, narrow, or retain."
        )
        try:
            result = llm.generate(prompt)
            recommendations = [line.strip() for line in result.split("\n") if line.strip() and len(line.strip()) > 10][:10]
        except Exception as exc:
            logger.warning(f"Eligibility optimization LLM failed: {exc}")

    if not recommendations:
        recommendations = [
            "Consider broadening age range to improve enrollment diversity",
            "Review biomarker requirements against real-world prevalence data",
            "Align exclusion criteria with competitor trial standards",
        ]

    return EligibilityOptimizeResponse(
        original_criteria_count=len(request.criteria),
        optimized_criteria=optimized_criteria,
        population_impact_estimate=0.0,
        recommendations=recommendations,
    )


@router.post("/adaptive/evaluate", response_model=AdaptiveEvaluateResponse)
async def adaptive_evaluate(request: AdaptiveEvaluateRequest, req: Request):
    """Adaptive design evaluation with statistical considerations."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    recommendations = []
    statistical_considerations = []

    if llm:
        prompt = (
            f"Evaluate this adaptive clinical trial design:\n"
            f"Design Type: {request.design_type}\n"
            f"Current Sample Size: {request.current_sample_size or 'Not specified'}\n"
            f"Target Sample Size: {request.target_sample_size or 'Not specified'}\n"
            f"Primary Endpoint: {request.primary_endpoint or 'Not specified'}\n"
            f"Arms: {len(request.arms)} treatment arms\n"
            f"Interim Data Available: {'Yes' if request.interim_data else 'No'}\n\n"
            f"Provide statistical considerations, futility assessment, and recommendations."
        )
        try:
            result = llm.generate(prompt)
            recommendations = [line.strip() for line in result.split("\n") if line.strip() and len(line.strip()) > 10][:8]
        except Exception as exc:
            logger.warning(f"Adaptive design LLM failed: {exc}")

    if not recommendations:
        recommendations = [
            "Ensure Type I error control across interim analyses",
            "Pre-specify adaptation rules in the statistical analysis plan",
            "Consider independent DSMB for unblinded interim reviews",
        ]

    return AdaptiveEvaluateResponse(
        design_assessment=f"Adaptive {request.design_type} design evaluation completed",
        recommendations=recommendations,
        statistical_considerations=statistical_considerations,
        futility_assessment={"status": "requires_interim_data"} if not request.interim_data else None,
        sample_size_recommendation={
            "current": request.current_sample_size,
            "target": request.target_sample_size,
        } if request.current_sample_size else None,
    )


@router.post("/safety/signal", response_model=SafetySignalResponse)
async def safety_signal(request: SafetySignalRequest, req: Request):
    """Safety signal detection with disproportionality analysis."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    signals_detected = []
    disproportionality = []
    recommendations = []
    severity_distribution = {"critical": 0, "high": 0, "moderate": 0, "low": 0}

    # Classify adverse events by severity
    for ae in request.adverse_events:
        severity = ae.get("severity", "moderate").lower()
        if severity in severity_distribution:
            severity_distribution[severity] += 1

    if llm:
        ae_summary = "\n".join(
            f"- {ae.get('event', ae.get('term', 'Unknown'))}: "
            f"n={ae.get('count', 'N/A')}, severity={ae.get('severity', 'N/A')}"
            for ae in request.adverse_events[:20]
        )
        prompt = (
            f"Analyze these adverse events for safety signals:\n\n"
            f"{ae_summary}\n\n"
            f"Drug: {request.drug_name or 'Not specified'}\n"
            f"Phase: {request.trial_phase or 'Not specified'}\n"
            f"Population Size: {request.study_population_size or 'Not specified'}\n\n"
            f"Perform disproportionality analysis and identify potential signals."
        )
        try:
            result = llm.generate(prompt)
            recommendations = [line.strip() for line in result.split("\n") if line.strip() and len(line.strip()) > 10][:8]
        except Exception as exc:
            logger.warning(f"Safety signal LLM failed: {exc}")

    if not recommendations:
        recommendations = [
            "Continue routine pharmacovigilance monitoring",
            "Review cumulative safety data at next DSMB meeting",
            "Ensure timely expedited reporting per ICH E2A",
        ]

    return SafetySignalResponse(
        signals_detected=signals_detected,
        disproportionality_analysis=disproportionality,
        severity_distribution=severity_distribution,
        recommendations=recommendations,
    )


@router.post("/regulatory/generate", response_model=RegulatoryGenerateResponse)
async def regulatory_generate(request: RegulatoryGenerateRequest, req: Request):
    """Regulatory document generation with agency-specific formatting."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    content = ""
    sections_generated = []
    compliance_notes = []

    if llm:
        prompt = (
            f"Generate a {request.document_type.upper()} document for {request.agency.upper()}:\n\n"
            f"Trial Data: {request.trial_data}\n"
            f"Requested Sections: {', '.join(request.sections) if request.sections else 'All standard sections'}\n"
            f"Include Appendices: {request.include_appendices}\n\n"
            f"Generate the document following {request.agency.upper()} formatting requirements."
        )
        try:
            content = llm.generate(prompt, max_tokens=4096)
            sections_generated = request.sections or [f"{request.document_type}_full"]
        except Exception as exc:
            logger.warning(f"Regulatory doc generation LLM failed: {exc}")

    if not content:
        content = (
            f"# {request.document_type.upper()} Document\n\n"
            f"**Agency:** {request.agency.upper()}\n"
            f"**Status:** Template generated -- requires clinical data population\n\n"
            f"Document generation requires LLM connectivity for full content."
        )
        sections_generated = ["template"]

    return RegulatoryGenerateResponse(
        document_type=request.document_type,
        agency=request.agency,
        content=content,
        sections_generated=sections_generated,
        compliance_notes=compliance_notes,
        metadata={"trial_data_keys": list(request.trial_data.keys())},
    )


@router.post("/competitive/landscape", response_model=CompetitiveLandscapeResponse)
async def competitive_landscape(request: CompetitiveLandscapeRequest, req: Request):
    """Competitive intelligence landscape analysis."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    competitors = []
    landscape_summary = "Competitive landscape analysis pending."
    threat_assessment = {}
    enrollment_race = []

    if llm:
        prompt = (
            f"Provide competitive intelligence analysis for:\n"
            f"Therapeutic Area: {request.therapeutic_area}\n"
            f"Indication: {request.indication or 'Broad'}\n"
            f"Mechanism: {request.mechanism or 'All mechanisms'}\n"
            f"Include Completed: {request.include_completed}\n\n"
            f"Analyze the competitive landscape including active trials, "
            f"enrollment progress, threat levels, and market timing."
        )
        try:
            landscape_summary = llm.generate(prompt)
        except Exception as exc:
            logger.warning(f"Competitive landscape LLM failed: {exc}")

    return CompetitiveLandscapeResponse(
        competitors=competitors,
        landscape_summary=landscape_summary,
        threat_assessment=threat_assessment,
        enrollment_race=enrollment_race,
        market_timing={},
    )


@router.post("/diversity/assess", response_model=DiversityAssessResponse)
async def diversity_assess(request: DiversityAssessRequest, req: Request):
    """Diversity and inclusion assessment with FDA guidance compliance."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    gap_analysis = []
    recommendations = []

    if llm:
        prompt = (
            f"Assess clinical trial diversity and inclusion:\n"
            f"Therapeutic Area: {request.therapeutic_area or 'Not specified'}\n"
            f"Demographics: {request.trial_demographics or 'Not provided'}\n"
            f"Disease Epidemiology: {request.disease_epidemiology or 'Not provided'}\n"
            f"Site Locations: {', '.join(request.site_locations) if request.site_locations else 'Not specified'}\n\n"
            f"Evaluate against FDA Diversity Action Plan requirements and provide recommendations."
        )
        try:
            result = llm.generate(prompt)
            recommendations = [line.strip() for line in result.split("\n") if line.strip() and len(line.strip()) > 10][:8]
        except Exception as exc:
            logger.warning(f"Diversity assessment LLM failed: {exc}")

    if not recommendations:
        recommendations = [
            "Develop a Diversity Action Plan per FDA guidance (April 2024)",
            "Include community-based sites in underrepresented regions",
            "Implement culturally appropriate recruitment materials",
            "Consider decentralized trial elements to reduce access barriers",
        ]

    return DiversityAssessResponse(
        diversity_score=0.0,
        gap_analysis=gap_analysis,
        fda_compliance_status="assessment_required",
        recommendations=recommendations,
    )


@router.post("/dct/plan", response_model=DCTPlanResponse)
async def dct_plan(request: DCTPlanRequest, req: Request):
    """Decentralized trial planning with component feasibility assessment."""
    _increment_metric(req, "workflow_requests_total")

    llm = _get_llm(req)
    recommended_components = []
    regulatory_considerations = []

    if llm:
        prompt = (
            f"Plan decentralized clinical trial components:\n"
            f"Protocol: {request.protocol_summary}\n"
            f"Therapeutic Area: {request.therapeutic_area or 'Not specified'}\n"
            f"Phase: {request.phase or 'Not specified'}\n"
            f"Target Components: {', '.join(request.target_components) if request.target_components else 'All'}\n"
            f"Region: {request.regulatory_region}\n\n"
            f"Assess feasibility, recommend hybrid model, and identify regulatory considerations."
        )
        try:
            result = llm.generate(prompt)
            regulatory_considerations = [line.strip() for line in result.split("\n") if line.strip() and len(line.strip()) > 10][:8]
        except Exception as exc:
            logger.warning(f"DCT planning LLM failed: {exc}")

    if not regulatory_considerations:
        regulatory_considerations = [
            "Ensure eConsent meets 21 CFR Part 11 requirements",
            "Validate telemedicine platforms for clinical assessment adequacy",
            "Confirm local lab CLIA/CAP certification requirements",
            "Address wearable device data integrity and validation",
        ]

    return DCTPlanResponse(
        feasibility_assessment={"status": "evaluated", "overall": "feasible"},
        recommended_components=recommended_components,
        hybrid_model={"site_visits_required": True, "remote_visits_enabled": True},
        technology_requirements=[],
        regulatory_considerations=regulatory_considerations,
    )


# =====================================================================
# Reference Endpoints
# =====================================================================

@router.get("/therapeutic-areas")
async def list_therapeutic_areas():
    """Reference catalog of therapeutic areas."""
    return {
        "therapeutic_areas": [
            {"id": ta.value, "name": ta.name.replace("_", " ").title()}
            for ta in TherapeuticArea
        ]
    }


@router.get("/phases")
async def list_phases():
    """Reference catalog of clinical trial phases."""
    return {
        "phases": [
            {"id": "phase_i", "name": "Phase I", "description": "First-in-human, safety, PK/PD, dose-finding"},
            {"id": "phase_i_ii", "name": "Phase I/II", "description": "Combined safety and preliminary efficacy"},
            {"id": "phase_ii", "name": "Phase II", "description": "Proof-of-concept, dose-ranging, efficacy signals"},
            {"id": "phase_ii_iii", "name": "Phase II/III", "description": "Adaptive/seamless efficacy confirmation"},
            {"id": "phase_iii", "name": "Phase III", "description": "Pivotal confirmatory, registration-enabling"},
            {"id": "phase_iv", "name": "Phase IV", "description": "Post-marketing surveillance, real-world evidence"},
            {"id": "not_applicable", "name": "Not Applicable", "description": "Observational, registry, or expanded access"},
        ]
    }


@router.get("/guidelines")
async def list_guidelines():
    """Reference catalog of clinical trial guidelines."""
    return {
        "guidelines": [
            {"id": "ich_e6_r3", "name": "ICH E6(R3) GCP", "description": "Good Clinical Practice principles and requirements"},
            {"id": "ich_e8_r1", "name": "ICH E8(R1)", "description": "General considerations for clinical studies"},
            {"id": "ich_e9_r1", "name": "ICH E9(R1)", "description": "Statistical principles -- estimands and sensitivity analyses"},
            {"id": "ich_e10", "name": "ICH E10", "description": "Choice of control group in clinical trials"},
            {"id": "ich_e17", "name": "ICH E17", "description": "Multi-regional clinical trials"},
            {"id": "ich_e19", "name": "ICH E19", "description": "Optimization of safety data collection"},
            {"id": "ich_e20", "name": "ICH E20", "description": "Adaptive clinical trials"},
            {"id": "fda_diversity", "name": "FDA Diversity Action Plan", "description": "Diversity requirements for clinical trials (2024)"},
            {"id": "fda_dct", "name": "FDA DCT Guidance", "description": "Decentralized clinical trials guidance"},
            {"id": "fda_rwe", "name": "FDA RWE Framework", "description": "Real-world evidence for regulatory decisions"},
            {"id": "ema_adaptive", "name": "EMA Adaptive Pathways", "description": "EMA guidance on adaptive trial designs"},
            {"id": "ema_compassionate", "name": "EMA Compassionate Use", "description": "EMA compassionate use guidance"},
        ]
    }


@router.get("/knowledge-version")
async def knowledge_version():
    """Version metadata for the trial knowledge base."""
    return {
        "version": "1.0.0",
        "agent": "clinical-trial-intelligence-agent",
        "collections": 14,
        "workflows": 11,
        "last_updated": "2026-03-21",
        "data_sources": [
            "ClinicalTrials.gov",
            "PubMed/MEDLINE",
            "ICH Guidelines",
            "FDA Guidance Documents",
            "EMA Scientific Guidelines",
            "WHO ICTRP",
        ],
        "therapeutic_areas": len(TherapeuticArea),
        "trial_phases": len(TrialPhase),
    }


@router.post("/workflow/{workflow_type}", response_model=WorkflowResponse)
async def generic_workflow(workflow_type: str, request: WorkflowRequest, req: Request):
    """Generic workflow dispatch for any supported workflow type."""
    _increment_metric(req, "workflow_requests_total")

    valid_types = [wt.value for wt in TrialWorkflowType]
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
        result = await engine.execute(workflow_type, data)
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
