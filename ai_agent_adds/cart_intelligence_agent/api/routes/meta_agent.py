"""Meta-agent route -- unified question-answering endpoint.

Accepts a natural-language question, routes it through the meta-agent
orchestrator, and returns a synthesised answer with provenance, confidence
score, and suggested follow-up questions.

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

import logging
import time
from typing import List, Optional

from fastapi import APIRouter, HTTPException

logger = logging.getLogger(__name__)
from pydantic import BaseModel, Field

from src.metrics import record_query, record_pipeline_stage

router = APIRouter(prefix="/api", tags=["meta-agent"])


# ── Cross-Agent Integration Endpoint ─────────────────────────────────

@router.post("/v1/cart/integrated-assessment")
async def integrated_assessment(request: dict):
    """Multi-agent integrated assessment combining insights from across the HCLS AI Factory.

    Queries biomarker, oncology, single-cell, cardiology, and clinical trial
    agents for a comprehensive CAR-T therapy assessment.
    """
    try:
        from src.cross_modal import (
            query_biomarker_agent,
            query_oncology_agent,
            query_single_cell_agent,
            query_cardiology_agent,
            query_trial_agent,
            integrate_cross_agent_results,
        )

        target_antigens = request.get("target_antigens", {})
        patient_profile = request.get("patient_profile", {})
        tumor_data = request.get("tumor_data", {})
        cart_product = request.get("cart_product", {})
        patient_id = request.get("patient_id", "")

        results = []

        # Query biomarker agent for target antigen expression
        if target_antigens:
            results.append(query_biomarker_agent(target_antigens))

        # Query oncology agent for tumor profile and disease context
        if patient_profile:
            results.append(query_oncology_agent(patient_profile))

        # Query single-cell agent for TME profiling
        if tumor_data:
            results.append(query_single_cell_agent(tumor_data))

        # Query cardiology agent for baseline cardiac assessment
        if patient_id:
            results.append(query_cardiology_agent(patient_id))

        # Query trial agent for CAR-T trial matching
        if cart_product and patient_profile:
            results.append(query_trial_agent(cart_product, patient_profile))

        integrated = integrate_cross_agent_results(results)
        return {
            "status": "completed",
            "assessment": integrated,
            "agents_consulted": integrated.get("agents_consulted", []),
        }
    except Exception as exc:
        logger.error(f"Integrated assessment failed: {exc}")
        return {"status": "partial", "assessment": {}, "error": "Cross-agent integration unavailable"}


# ── Request / Response schemas ───────────────────────────────────────

class AskRequest(BaseModel):
    """Payload for POST /api/ask."""

    question: str = Field(..., min_length=1, description="Natural-language question")
    target_gene: Optional[str] = Field(
        None, description="Optional gene symbol to focus the query (e.g. CD19, BCMA)"
    )
    patient_id: Optional[str] = Field(
        None, description="Optional patient identifier for contextual queries"
    )


class SourceRef(BaseModel):
    """A single evidence source backing the answer."""

    collection: str = Field(..., description="Source collection name")
    doc_id: str = Field(..., description="Document / record identifier")
    title: str = Field("", description="Short title or snippet")
    score: float = Field(0.0, description="Relevance score")


class AskResponse(BaseModel):
    """Structured response from POST /api/ask."""

    answer: str = Field(..., description="Synthesised answer text")
    sources: List[SourceRef] = Field(default_factory=list, description="Provenance trail")
    confidence: float = Field(
        0.0, ge=0.0, le=1.0, description="Confidence score (0-1)"
    )
    follow_up_questions: List[str] = Field(
        default_factory=list, description="Suggested follow-up questions"
    )
    processing_time_ms: float = Field(0.0, description="Server-side latency in ms")


# ── Endpoint ─────────────────────────────────────────────────────────

@router.post("/ask", response_model=AskResponse)
async def ask(request: AskRequest):
    """Accept a question, route through the meta-agent, and return a
    synthesised answer with sources and confidence.

    The meta-agent orchestrates retrieval across all CAR-T collections,
    augments with the knowledge graph, and synthesises via LLM.
    """
    t0 = time.perf_counter()

    # Import engine state from main app module (populated during lifespan)
    try:
        from api.main import _engine
    except ImportError:
        raise HTTPException(status_code=503, detail="Engine module unavailable")

    if not _engine:
        raise HTTPException(status_code=503, detail="Engine not initialized")

    if not _engine.embedder or not _engine.llm:
        raise HTTPException(
            status_code=503,
            detail="Embedding model or LLM client not available",
        )

    try:
        from src.models import AgentQuery
        from src.rag_engine import CART_SYSTEM_PROMPT

        agent_query = AgentQuery(
            question=request.question,
            target_antigen=request.target_gene,
        )

        # Retrieve evidence across collections
        evidence = _engine.retrieve(query=agent_query)

        # Generate LLM synthesis
        prompt_text = _engine._build_prompt(request.question, evidence)
        answer = _engine.llm.generate(
            prompt=prompt_text,
            system_prompt=CART_SYSTEM_PROMPT,
            max_tokens=2048,
            temperature=0.7,
        )

        # Build source references
        sources = [
            SourceRef(
                collection=h.collection,
                doc_id=h.id,
                title=h.text[:120] if h.text else "",
                score=h.score,
            )
            for h in evidence.hits
        ]

        # Derive confidence from mean evidence score
        confidence = 0.0
        if sources:
            confidence = min(1.0, sum(s.score for s in sources) / len(sources))

        elapsed_ms = (time.perf_counter() - t0) * 1000

        # Record metrics
        record_query("meta_agent", elapsed_ms / 1000, len(sources))
        record_pipeline_stage("meta_agent_ask", elapsed_ms / 1000)

        return AskResponse(
            answer=answer,
            sources=sources,
            confidence=round(confidence, 3),
            follow_up_questions=[],
            processing_time_ms=round(elapsed_ms, 1),
        )

    except HTTPException:
        raise
    except Exception as exc:
        logger.error(f"Meta-agent query failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")
