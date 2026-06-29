"""Meta-agent endpoint for Imaging Intelligence Agent.

POST /api/ask -- Accepts a natural language question and returns a
synthesized answer with evidence sources and follow-up questions.

Author: Adam Jones
Date: February 2026
"""

import time
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from loguru import logger
from pydantic import BaseModel, Field

# =====================================================================
# Request / Response Models
# =====================================================================


class AskRequest(BaseModel):
    """Meta-agent question request."""
    question: str = Field(..., min_length=3, max_length=2000, description="Natural language question")
    modality: Optional[str] = Field(None, description="Filter by imaging modality (ct, mri, xray, etc.)")
    body_region: Optional[str] = Field(None, description="Filter by body region (head, chest, etc.)")
    top_k: int = Field(5, ge=1, le=30, description="Results per collection")
    conversation_history: List[Dict[str, str]] = Field(
        default_factory=list,
        description="Prior conversation turns as [{role, content}, ...]",
    )
    include_follow_ups: bool = Field(True, description="Generate follow-up question suggestions")


class SourceReference(BaseModel):
    """A single evidence source used in the answer."""
    collection: str
    id: str
    score: float
    text_snippet: str = Field("", max_length=500)
    metadata: Dict[str, Any] = {}


class AskResponse(BaseModel):
    """Meta-agent answer response."""
    question: str
    answer: str
    sources: List[SourceReference]
    follow_up_questions: List[str] = []
    collections_searched: int
    evidence_count: int
    search_time_ms: float
    total_time_ms: float
    knowledge_context_used: bool = False
    nim_services_available: List[str] = []


# =====================================================================
# Follow-up Question Generation
# =====================================================================

# Template-based follow-up generation keyed to detected query topics
FOLLOW_UP_TEMPLATES = {
    "hemorrhage": [
        "What are the Brain Trauma Foundation guidelines for surgical intervention thresholds?",
        "How does VISTA-3D segment hemorrhage subtypes (epidural vs subdural)?",
        "What FDA-cleared AI devices exist for hemorrhage detection?",
    ],
    "nodule": [
        "What is the ACR Lung-RADS v2022 management algorithm for category 4A nodules?",
        "How does volume doubling time help differentiate benign from malignant nodules?",
        "Which AI models achieve the best sensitivity for sub-6mm nodules?",
    ],
    "ms": [
        "What are the MAGNIMS criteria for MS disease activity classification?",
        "How does lesion volume change correlate with clinical disability (EDSS)?",
        "What role does VISTA-3D play in brain lesion segmentation?",
    ],
    "pneumonia": [
        "What is the sensitivity of DenseNet-121 for consolidation detection on CXR?",
        "How do CT and CXR compare for community-acquired pneumonia diagnosis?",
        "What public datasets are available for training pneumonia detection models?",
    ],
    "segmentation": [
        "How does VISTA-3D compare to nnU-Net for multi-organ segmentation?",
        "What is the Dice score of current MONAI models on the BTCV dataset?",
        "Can MAISI-generated synthetic CT data improve segmentation model training?",
    ],
}

DEFAULT_FOLLOW_UPS = [
    "What FDA-cleared AI devices are available for this imaging modality?",
    "What clinical guidelines apply to this finding?",
    "What are the latest clinical trial results in this area?",
]


def generate_follow_ups(question: str, answer: str) -> List[str]:
    """Generate relevant follow-up questions based on the query topic."""
    question_lower = question.lower()

    for topic, follow_ups in FOLLOW_UP_TEMPLATES.items():
        if topic in question_lower:
            return follow_ups[:3]

    # Generic follow-ups based on answer content
    answer_lower = answer.lower()
    for topic, follow_ups in FOLLOW_UP_TEMPLATES.items():
        if topic in answer_lower:
            return follow_ups[:3]

    return DEFAULT_FOLLOW_UPS


# =====================================================================
# Router
# =====================================================================

router = APIRouter()


def _get_state():
    """Import and access the shared application state from api.main."""
    from api.main import _state
    return _state


@router.post("/ask", response_model=AskResponse)
async def ask(request: AskRequest):
    """Meta-agent endpoint: full RAG question answering with sources and follow-ups.

    Orchestrates:
      1. Multi-collection vector search across 10+ imaging collections
      2. Domain knowledge graph augmentation
      3. LLM synthesis with clinical reasoning
      4. Source extraction and scoring
      5. Follow-up question generation
    """
    start = time.time()
    state = _get_state()

    engine = state.get("engine")
    if engine is None:
        raise HTTPException(status_code=503, detail="Agent engine not initialized")

    nim_manager = state.get("nim_manager")

    # Build search kwargs
    kwargs = {"top_k_per_collection": request.top_k}
    if request.modality:
        kwargs["modality_filter"] = request.modality
    if request.body_region:
        kwargs["body_region_filter"] = request.body_region

    # Build conversation context from history
    conversation_ctx = ""
    if request.conversation_history:
        turns = request.conversation_history[-6:]  # Last 3 exchanges
        conversation_ctx = "\n".join(
            f"{t.get('role', 'user').upper()}: {t.get('content', '')[:500]}"
            for t in turns
        )

    try:
        # Step 1: Retrieve evidence
        evidence = engine.retrieve(request.question, **kwargs)

        # Step 2: Synthesize answer via LLM
        answer = engine.query(
            request.question,
            conversation_context=conversation_ctx,
            **kwargs,
        )

        # Step 3: Build source references
        sources = []
        for hit in evidence.hits[:20]:
            sources.append(SourceReference(
                collection=hit.collection,
                id=hit.id,
                score=hit.score,
                text_snippet=hit.text[:500] if hit.text else "",
                metadata=hit.metadata,
            ))

        # Step 4: Generate follow-up questions
        follow_ups = []
        if request.include_follow_ups:
            follow_ups = generate_follow_ups(request.question, answer)

        # Step 5: Check NIM availability
        nim_available = []
        if nim_manager:
            try:
                status = nim_manager.check_all_services()
                nim_available = [name for name, s in status.items() if s in ("available", "cloud", "anthropic", "mock")]
            except Exception:
                pass

        total_time = (time.time() - start) * 1000

        return AskResponse(
            question=request.question,
            answer=answer,
            sources=sources,
            follow_up_questions=follow_ups,
            collections_searched=evidence.total_collections_searched,
            evidence_count=evidence.hit_count,
            search_time_ms=evidence.search_time_ms,
            total_time_ms=total_time,
            knowledge_context_used=bool(evidence.knowledge_context),
            nim_services_available=nim_available,
        )
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"Meta-agent /ask failed: {e}")
        raise HTTPException(status_code=500, detail=f"Agent query failed: {e}")
