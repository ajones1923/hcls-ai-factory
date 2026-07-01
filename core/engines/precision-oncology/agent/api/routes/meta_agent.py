"""
Oncology Intelligence Agent - Meta Agent Router
==============================================
Unified RAG Q&A endpoint that routes through the intelligence agent or
falls back to the RAG engine for evidence retrieval and answer generation.

Author: Adam Jones
Date: February 2026
"""

import time
import logging
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

router = APIRouter(tags=["meta-agent"])


# ---------------------------------------------------------------------------
# Cross-Agent Integration Endpoint
# ---------------------------------------------------------------------------
@router.post("/v1/onco/integrated-assessment")
async def integrated_assessment(request: dict):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries CAR-T, biomarker, clinical trial, cardiology, neurology, PGx,
    imaging, and single-cell agents for a comprehensive oncology assessment.
    """
    try:
        from src.cross_agent import (
            query_cart_agent,
            query_biomarker_agent,
            query_trial_agent,
            query_cardiology_agent,
            query_neurology_agent,
            query_pgx_agent,
            query_imaging_agent,
            query_single_cell_agent,
            integrate_cross_agent_results,
        )

        patient_data = request.get("patient_data", {})
        tumor_profile = request.get("tumor_profile", {})
        therapy_plan = request.get("therapy_plan", {})
        mutations = request.get("mutations", [])
        drug_list = request.get("drug_list", [])

        results = []

        # Query CAR-T agent for immunotherapy eligibility
        target_antigens = request.get("target_antigens", [])
        if target_antigens and patient_data:
            results.append(query_cart_agent(target_antigens, patient_data))

        # Query biomarker agent for risk stratification
        if tumor_profile:
            results.append(query_biomarker_agent(tumor_profile))

        # Query trial agent for precision medicine trials
        if patient_data and mutations:
            results.append(query_trial_agent(patient_data, mutations))

        # Query cardiology agent for cardiotoxicity risk
        if therapy_plan:
            results.append(query_cardiology_agent(therapy_plan))

        # Query neurology agent for neurotoxicity risk
        if therapy_plan:
            results.append(query_neurology_agent(therapy_plan))

        # Query PGx agent for drug-gene interactions
        if drug_list and patient_data.get("patient_id"):
            results.append(query_pgx_agent(drug_list, patient_data["patient_id"]))

        # Query imaging agent for staging protocol
        cancer_type = tumor_profile.get("cancer_type", "")
        stage = tumor_profile.get("stage", "")
        if cancer_type:
            results.append(query_imaging_agent(cancer_type, stage))

        # Query single-cell agent for TME profiling
        tumor_data = request.get("tumor_data", {})
        if tumor_data:
            results.append(query_single_cell_agent(tumor_data))

        integrated = integrate_cross_agent_results(results)
        return {
            "status": "completed",
            "assessment": integrated,
            "agents_consulted": integrated.get("agents_consulted", []),
        }
    except Exception as exc:
        logger.error(f"Integrated assessment failed: {exc}")
        return {"status": "partial", "assessment": {}, "error": "Cross-agent integration unavailable"}


# ---------------------------------------------------------------------------
# Schemas
# ---------------------------------------------------------------------------
class AskRequest(BaseModel):
    question: str = Field(..., min_length=3, description="Clinical question")
    cancer_type: Optional[str] = Field(None, description="Cancer type filter")
    gene: Optional[str] = Field(None, description="Gene symbol filter")
    top_k: int = Field(default=10, ge=1, le=100)
    conversation_history: Optional[List[Dict[str, str]]] = Field(
        default=None,
        description="Prior turns as [{role, content}, ...]",
    )
    include_follow_ups: bool = Field(
        default=True,
        description="Generate follow-up question suggestions",
    )


class SourceItem(BaseModel):
    collection: str
    text: str
    score: float
    metadata: Dict[str, Any] = Field(default_factory=dict)


class AskResponse(BaseModel):
    answer: str
    sources: List[SourceItem] = Field(default_factory=list)
    confidence: float = 0.0
    follow_up_questions: List[str] = Field(default_factory=list)
    processing_time_ms: float = 0.0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_COMPARATIVE_KEYWORDS = [
    "compare",
    "versus",
    " vs ",
    "difference between",
    "better than",
    "prefer",
    "superior",
    "inferior",
    "head-to-head",
    "which is",
    "combination vs",
]


def _is_comparative(question: str) -> bool:
    """Detect whether the question is a comparative query."""
    q_lower = question.lower()
    return any(kw in q_lower for kw in _COMPARATIVE_KEYWORDS)


# ---------------------------------------------------------------------------
# Endpoint
# ---------------------------------------------------------------------------
@router.post("/api/ask", response_model=AskResponse)
async def ask(req: AskRequest):
    """
    Unified clinical Q&A endpoint.

    Routes comparative or complex queries through the OncoIntelligenceAgent
    for multi-step reasoning; simpler evidence queries go directly through
    the OncoRAGEngine.
    """
    from api.main import get_state

    state = get_state()
    agent = state.get("intelligence_agent")
    rag = state.get("rag_engine")

    if agent is None and rag is None:
        raise HTTPException(status_code=503, detail="No query service available")

    t0 = time.time()

    try:
        # Comparative / complex queries -> Intelligence Agent
        if agent is not None and _is_comparative(req.question):
            logger.info("Routing comparative query through intelligence agent")
            result = await agent.run(
                question=req.question,
                cancer_type=req.cancer_type,
                gene=req.gene,
                top_k=req.top_k,
                conversation_history=req.conversation_history,
            )
        elif agent is not None:
            # Attempt full agent pipeline; fall back to RAG on failure
            try:
                result = await agent.run(
                    question=req.question,
                    cancer_type=req.cancer_type,
                    gene=req.gene,
                    top_k=req.top_k,
                    conversation_history=req.conversation_history,
                )
            except Exception:
                logger.warning(
                    "Intelligence agent failed, falling back to RAG engine"
                )
                result = await rag.query(
                    question=req.question,
                    cancer_type=req.cancer_type,
                    gene=req.gene,
                    top_k=req.top_k,
                )
        else:
            result = await rag.query(
                question=req.question,
                cancer_type=req.cancer_type,
                gene=req.gene,
                top_k=req.top_k,
            )
    except Exception as exc:
        logger.error("Query failed: %s", exc, exc_info=True)
        raise HTTPException(status_code=500, detail="Internal processing error")

    elapsed_ms = round((time.time() - t0) * 1000, 1)

    # Normalise result into response schema
    sources = []
    for src in result.get("sources", []):
        sources.append(
            SourceItem(
                collection=src.get("collection", "unknown"),
                text=src.get("text", ""),
                score=float(src.get("score", 0.0)),
                metadata=src.get("metadata", {}),
            )
        )

    follow_ups: List[str] = []
    if req.include_follow_ups:
        follow_ups = result.get("follow_up_questions", [])

    return AskResponse(
        answer=result.get("answer", "No answer generated."),
        sources=sources,
        confidence=float(result.get("confidence", 0.0)),
        follow_up_questions=follow_ups,
        processing_time_ms=elapsed_ms,
    )
