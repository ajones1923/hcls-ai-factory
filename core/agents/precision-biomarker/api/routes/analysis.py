"""Analysis-specific routes for the Biomarker Intelligence Agent.

Provides endpoints for:
  - Full patient analysis (all modules)
  - Biological age calculation
  - Disease risk trajectory analysis
  - Pharmacogenomic mapping
  - RAG Q&A query
  - Health check

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import json
import os
import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import StreamingResponse
from loguru import logger
from pydantic import BaseModel, Field, model_validator

from src.audit import AuditAction, audit_log

router = APIRouter(prefix="/v1", tags=["analysis"])


# =====================================================================
# Cross-Agent Integration Endpoint
# =====================================================================

@router.post("/biomarker/integrated-assessment")
async def integrated_assessment(request: dict, req: Request):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries oncology, CAR-T, pharmacogenomics, and clinical trial agents
    for a comprehensive biomarker-driven assessment.
    """
    try:
        from src.cross_modal import (
            integrate_cross_agent_results,
            query_cart_agent,
            query_oncology_agent,
            query_pgx_agent,
            query_trial_agent,
        )

        biomarker_panel = request.get("biomarker_panel", {})
        target_antigens = request.get("target_antigens", {})
        drug_list = request.get("drug_list", {})

        results = []

        # Query oncology agent for tumor-biomarker correlation
        if biomarker_panel:
            results.append(query_oncology_agent(biomarker_panel))

        # Query CAR-T agent for target validation
        if target_antigens:
            results.append(query_cart_agent(target_antigens))

        # Query PGx agent for drug-gene interactions
        if drug_list:
            results.append(query_pgx_agent(drug_list))

        # Query trial agent for biomarker-driven trials
        if biomarker_panel:
            results.append(query_trial_agent(biomarker_panel))

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

class PatientProfileRequest(BaseModel):
    """Full patient profile for analysis."""
    patient_id: str = Field(..., description="Unique patient identifier")
    age: int = Field(..., ge=0, le=150, description="Patient age in years")
    sex: str = Field(..., description="Patient sex", pattern="^(M|F)$")
    biomarkers: Dict[str, float] = Field(
        default_factory=dict,
        description="Biomarker name -> measured value",
    )
    genotypes: Dict[str, str] = Field(
        default_factory=dict,
        description="rsID -> genotype mapping (e.g., {'rs1801133': 'CT'})",
    )
    star_alleles: Dict[str, str] = Field(
        default_factory=dict,
        description="Gene -> star allele mapping (e.g., {'CYP2D6': '*1/*2'})",
    )

    @model_validator(mode="after")
    def _check_has_data(self) -> "PatientProfileRequest":
        if not self.biomarkers and not self.genotypes and not self.star_alleles:
            raise ValueError("At least one of biomarkers, genotypes, or star_alleles must be provided")
        return self


class BiologicalAgeRequest(BaseModel):
    """Request for biological age calculation only."""
    age: int = Field(..., ge=0, le=150, description="Chronological age")
    biomarkers: Dict[str, float] = Field(
        ..., description="PhenoAge biomarkers (albumin, creatinine, glucose, etc.)",
    )


class DiseaseRiskRequest(BaseModel):
    """Request for disease trajectory analysis."""
    age: int = Field(..., ge=0, le=150)
    sex: str = Field(..., description="Patient sex", pattern="^(M|F)$")
    biomarkers: Dict[str, float] = Field(default_factory=dict)
    genotypes: Dict[str, str] = Field(default_factory=dict)


class PGxRequest(BaseModel):
    """Request for pharmacogenomic mapping."""
    star_alleles: Dict[str, str] = Field(
        default_factory=dict,
        description="Gene -> star allele (e.g., {'CYP2D6': '*4/*4'})",
    )
    genotypes: Dict[str, str] = Field(
        default_factory=dict,
        description="rsID -> genotype for PGx-relevant SNPs",
    )


class QueryRequest(BaseModel):
    """Request for RAG Q&A query."""
    question: str = Field(..., min_length=1, max_length=5000, description="Natural-language question")
    patient_profile: Optional[PatientProfileRequest] = None
    collections: Optional[List[str]] = Field(None, description="Restrict to collections")
    year_min: Optional[int] = Field(None, ge=1990, le=2030)
    year_max: Optional[int] = Field(None, ge=1990, le=2030)

    @model_validator(mode="after")
    def _check_year_range(self) -> "QueryRequest":
        if self.year_min is not None and self.year_max is not None:
            if self.year_min > self.year_max:
                raise ValueError("year_min must be <= year_max")
        return self


class AnalysisResponse(BaseModel):
    """Response from full patient analysis."""
    patient_id: str
    biological_age: float
    age_acceleration: float
    mortality_risk: float
    disease_trajectories: List[Dict[str, Any]]
    pgx_results: List[Dict[str, Any]]
    genotype_adjustments: List[Dict[str, Any]]
    critical_alerts: List[str]
    processing_time_ms: float


class BiologicalAgeResponse(BaseModel):
    """Response from biological age calculation."""
    chronological_age: int
    biological_age: float
    age_acceleration: float
    mortality_risk: str
    top_aging_drivers: List[Dict[str, Any]]
    grimage: Optional[Dict[str, Any]] = None


class DiseaseRiskResponse(BaseModel):
    """Response from disease risk analysis."""
    trajectories: List[Dict[str, Any]]
    processing_time_ms: float


class PGxResponse(BaseModel):
    """Response from pharmacogenomic mapping."""
    results: List[Dict[str, Any]]
    total_genes: int
    critical_findings: List[str]


class QueryResponse(BaseModel):
    """Response from RAG Q&A query."""
    question: str
    answer: str
    evidence_count: int
    collections_searched: int
    search_time_ms: float


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/analyze", response_model=AnalysisResponse)
def full_analysis(request: PatientProfileRequest, req: Request):
    """Run full patient analysis using all modules.

    Executes biological age calculation, disease trajectory analysis,
    pharmacogenomic mapping, and genotype-adjusted reference ranges.
    Returns critical alerts and comprehensive results.
    """
    agent = getattr(req.app.state, "agent", None)
    if not agent:
        raise HTTPException(status_code=503, detail="Agent not initialized")

    t0 = time.perf_counter()

    try:
        audit_log(
            AuditAction.PATIENT_ANALYSIS,
            patient_id=request.patient_id,
            source_ip=req.client.host if req.client else None,
            details={"modules": ["biological_age", "disease_trajectory", "pgx", "genotype_adjustment"]},
        )

        from src.models import PatientProfile
        profile = PatientProfile(
            patient_id=request.patient_id,
            age=request.age,
            sex=request.sex,
            biomarkers=request.biomarkers,
            genotypes=request.genotypes,
            star_alleles=request.star_alleles,
        )
        result = agent.analyze_patient(profile)
        elapsed = (time.perf_counter() - t0) * 1000

        response = AnalysisResponse(
            patient_id=request.patient_id,
            biological_age=result.biological_age.biological_age,
            age_acceleration=result.biological_age.age_acceleration,
            mortality_risk=result.biological_age.mortality_risk,
            disease_trajectories=[t.model_dump() for t in result.disease_trajectories],
            pgx_results=[p.model_dump() for p in result.pgx_results],
            genotype_adjustments=[g.model_dump() for g in result.genotype_adjustments],
            critical_alerts=result.critical_alerts,
            processing_time_ms=round(elapsed, 1),
        )

        # Publish cross-agent events
        try:
            import sys
            _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
            sys.path.insert(0, _lib_path)
            from hcls_common.event_bus import EventPriority, EventType, PipelineStage, publish_event

            publish_event(
                EventType.BIOLOGICAL_AGE_COMPUTED,
                source_stage=PipelineStage.ANNOTATION,
                payload={
                    "biological_age": result.biological_age.biological_age,
                    "age_acceleration": result.biological_age.age_acceleration,
                    "mortality_risk": result.biological_age.mortality_risk,
                },
                patient_id=request.patient_id,
            )

            # Publish alerts for critical findings
            for alert in result.critical_alerts:
                publish_event(
                    EventType.BIOMARKER_ALERT,
                    source_stage=PipelineStage.ANNOTATION,
                    payload={"alert": alert},
                    patient_id=request.patient_id,
                    priority=EventPriority.HIGH,
                )

            # Publish PGx results for drug discovery feedback
            if result.pgx_results:
                publish_event(
                    EventType.PGX_RESULT_READY,
                    source_stage=PipelineStage.ANNOTATION,
                    payload={
                        "pgx_results": [p.model_dump() for p in result.pgx_results],
                        "genes_analyzed": len(result.pgx_results),
                    },
                    patient_id=request.patient_id,
                )
        except Exception:
            logger.debug("Cross-agent event publishing unavailable (non-critical)")

        return response
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Full analysis failed: {e}")
        raise HTTPException(status_code=500, detail="Internal analysis error. Check server logs.")


@router.post("/biological-age", response_model=BiologicalAgeResponse)
def biological_age(request: BiologicalAgeRequest, req: Request):
    """Calculate biological age from blood biomarkers.

    Uses PhenoAge (Levine 2018) algorithm with 9 routine blood biomarkers.
    Also computes GrimAge surrogate if plasma protein markers are available.
    """
    calc = getattr(req.app.state, "bio_age_calc", None)
    if not calc:
        raise HTTPException(status_code=503, detail="Bio age calculator not initialized")

    try:
        audit_log(AuditAction.BIOLOGICAL_AGE, source_ip=req.client.host if req.client else None)

        result = calc.calculate(request.age, request.biomarkers)
        phenoage = result.get("phenoage", {})

        response = BiologicalAgeResponse(
            chronological_age=request.age,
            biological_age=result.get("biological_age", request.age),
            age_acceleration=result.get("age_acceleration", 0.0),
            mortality_risk=phenoage.get("mortality_risk", "UNKNOWN"),
            top_aging_drivers=phenoage.get("top_aging_drivers", []),
            grimage=result.get("grimage"),
        )

        # Publish cross-agent event
        try:
            import sys
            _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
            sys.path.insert(0, _lib_path)
            from hcls_common.event_bus import EventType, PipelineStage, publish_event
            publish_event(
                EventType.BIOLOGICAL_AGE_COMPUTED,
                source_stage=PipelineStage.ANNOTATION,
                payload={
                    "biological_age": result.get("biological_age"),
                    "age_acceleration": result.get("age_acceleration"),
                },
            )
        except Exception:
            logger.debug("Cross-agent event publishing unavailable")

        return response
    except Exception as e:
        logger.exception(f"Biological age calculation failed: {e}")
        raise HTTPException(status_code=500, detail="Internal analysis error. Check server logs.")


@router.post("/disease-risk", response_model=DiseaseRiskResponse)
def disease_risk(request: DiseaseRiskRequest, req: Request):
    """Analyze disease trajectory risk across 6 disease categories.

    Detects pre-symptomatic disease patterns using genotype-stratified
    biomarker thresholds for diabetes, cardiovascular, liver, thyroid,
    iron metabolism, and nutritional conditions.
    """
    analyzer = getattr(req.app.state, "trajectory_analyzer", None)
    if not analyzer:
        raise HTTPException(status_code=503, detail="Trajectory analyzer not initialized")

    t0 = time.perf_counter()

    try:
        audit_log(AuditAction.DISEASE_RISK, source_ip=req.client.host if req.client else None)

        trajectories = analyzer.analyze_all(
            request.biomarkers, request.genotypes, request.age, request.sex,
        )
        elapsed = (time.perf_counter() - t0) * 1000

        response = DiseaseRiskResponse(
            trajectories=[t.model_dump() if hasattr(t, "model_dump") else t for t in trajectories],
            processing_time_ms=round(elapsed, 1),
        )

        # Publish cross-agent event for critical disease trajectories
        try:
            import sys
            _lib_path = os.environ.get("HCLS_LIB_PATH", "/app/lib")
            sys.path.insert(0, _lib_path)
            from hcls_common.event_bus import EventPriority, EventType, PipelineStage, publish_event
            critical_trajectories = []
            for t in trajectories:
                risk = t.risk_level if hasattr(t, "risk_level") else t.get("risk_level") if isinstance(t, dict) else None
                if risk in ("HIGH", "CRITICAL"):
                    critical_trajectories.append(t.model_dump() if hasattr(t, "model_dump") else t)
            if critical_trajectories:
                publish_event(
                    EventType.DISEASE_TRAJECTORY_ALERT,
                    source_stage=PipelineStage.ANNOTATION,
                    payload={"trajectories": critical_trajectories},
                    priority=EventPriority.HIGH,
                )
        except Exception:
            logger.debug("Cross-agent event publishing unavailable")

        return response
    except Exception as e:
        logger.exception(f"Disease risk analysis failed: {e}")
        raise HTTPException(status_code=500, detail="Internal analysis error. Check server logs.")


@router.post("/pgx", response_model=PGxResponse)
def pharmacogenomics(request: PGxRequest, req: Request):
    """Map star alleles and genotypes to CPIC drug recommendations.

    Supports CYP2D6, CYP2C19, CYP2C9, VKORC1, SLCO1B1, DPYD, and CYP3A5.
    Returns metabolizer phenotype and drug-specific dosing guidance.
    """
    mapper = getattr(req.app.state, "pgx_mapper", None)
    if not mapper:
        raise HTTPException(status_code=503, detail="PGx mapper not initialized")

    try:
        audit_log(AuditAction.PGX_MAPPING, source_ip=req.client.host if req.client else None)

        results = mapper.map_all(request.star_alleles, request.genotypes)

        # Extract critical findings
        critical = []
        for pgx_result in results.get("gene_results", []):
            phenotype = pgx_result.get("phenotype", "")
            if phenotype and phenotype.lower().replace(" ", "_").replace("-", "_") in ("poor_metabolizer", "ultra_rapid_metabolizer", "poor", "ultra_rapid"):
                critical.append(
                    f"{pgx_result.get('gene', '')} {pgx_result.get('star_alleles', '')}: {phenotype}"
                )

        return PGxResponse(
            results=results.get("gene_results", []),
            total_genes=results.get("genes_analyzed", 0),
            critical_findings=critical,
        )
    except Exception as e:
        logger.exception(f"PGx mapping failed: {e}")
        raise HTTPException(status_code=500, detail="Internal analysis error. Check server logs.")


@router.post("/query", response_model=QueryResponse)
def rag_query(request: QueryRequest, req: Request):
    """RAG Q&A query with evidence retrieval and LLM synthesis.

    Optionally accepts a patient profile for personalized context injection
    in the LLM prompt.
    """
    engine = getattr(req.app.state, "engine", None)
    if not engine:
        raise HTTPException(status_code=503, detail="RAG engine not initialized")
    if not engine.llm:
        raise HTTPException(status_code=503, detail="LLM client not available")
    if not engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    try:
        audit_log(
            AuditAction.RAG_QUERY,
            patient_id=request.patient_profile.patient_id if request.patient_profile else None,
            source_ip=req.client.host if req.client else None,
        )

        patient_profile = None
        if request.patient_profile:
            from src.models import PatientProfile
            patient_profile = PatientProfile(**request.patient_profile.model_dump())

        # Retrieve evidence once, then generate
        from src.models import AgentQuery
        from src.rag_engine import BIOMARKER_SYSTEM_PROMPT
        agent_query = AgentQuery(question=request.question, patient_profile=patient_profile)
        evidence = engine.retrieve(
            agent_query,
            collections_filter=request.collections,
            year_min=request.year_min,
            year_max=request.year_max,
        )
        prompt = engine._build_prompt(request.question, evidence, patient_profile)
        answer = engine.llm.generate(
            prompt=prompt,
            system_prompt=BIOMARKER_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        return QueryResponse(
            question=request.question,
            answer=answer,
            evidence_count=evidence.hit_count,
            collections_searched=evidence.total_collections_searched,
            search_time_ms=evidence.search_time_ms,
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"RAG query failed: {e}")
        raise HTTPException(status_code=500, detail="Internal analysis error. Check server logs.")


@router.post("/query/stream")
def rag_query_stream(request: QueryRequest, req: Request):
    """Streaming RAG Q&A query. Returns SSE stream with evidence then answer tokens."""
    engine = getattr(req.app.state, "engine", None)
    if not engine:
        raise HTTPException(status_code=503, detail="RAG engine not initialized")
    if not engine.llm:
        raise HTTPException(status_code=503, detail="LLM client not available")
    if not engine.embedder:
        raise HTTPException(status_code=503, detail="Embedding model not loaded")

    patient_profile = None
    if request.patient_profile:
        from src.models import PatientProfile
        patient_profile = PatientProfile(**request.patient_profile.model_dump())

    def event_generator():
        try:
            for chunk in engine.query_stream(
                question=request.question,
                patient_profile=patient_profile,
                collections_filter=request.collections,
                year_min=request.year_min,
                year_max=request.year_max,
            ):
                if chunk["type"] == "token":
                    yield f"data: {json.dumps({'type': 'token', 'content': chunk['content']})}\n\n"
                elif chunk["type"] == "done":
                    yield f"data: {json.dumps({'type': 'done'})}\n\n"
        except Exception as e:
            logger.exception(f"Streaming query failed: {e}")
            yield f"data: {json.dumps({'type': 'error', 'message': str(e)})}\n\n"

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
    )


@router.get("/health")
async def v1_health(req: Request):
    """V1 health check endpoint."""
    agent = getattr(req.app.state, "agent", None)
    engine = getattr(req.app.state, "engine", None)

    return {
        "status": "healthy" if agent else "degraded",
        "agent_ready": agent is not None,
        "rag_engine_ready": engine is not None,
        "llm_available": engine.llm is not None if engine else False,
        "embedder_available": engine.embedder is not None if engine else False,
    }
