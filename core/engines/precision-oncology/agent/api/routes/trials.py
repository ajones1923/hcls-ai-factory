"""
Oncology Intelligence Agent - Trial Matching Router
==================================================
Match clinical trials to patient profiles and rank therapies based on
molecular evidence, biomarker status, and guideline concordance.

Author: Adam Jones
Date: February 2026
"""

import time
import logging
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

router = APIRouter(tags=["trials"])


# ---------------------------------------------------------------------------
# Schemas
# ---------------------------------------------------------------------------
class TrialMatchRequest(BaseModel):
    cancer_type: str
    biomarkers: Dict[str, Any] = Field(default_factory=dict)
    stage: Optional[str] = None
    age: Optional[int] = Field(None, ge=0, le=120)
    top_k: int = Field(default=10, ge=1, le=50)


class TherapyRankRequest(BaseModel):
    cancer_type: str
    variants: List[Dict[str, Any]] = Field(default_factory=list)
    biomarkers: Dict[str, Any] = Field(default_factory=dict)
    prior_therapies: List[str] = Field(default_factory=list)


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@router.post("/api/trials/match")
async def match_trials(req: TrialMatchRequest):
    """Match clinical trials to a patient profile based on biomarkers and stage."""
    from api.main import get_state

    state = get_state()
    trial_matcher = state.get("trial_matcher")
    if trial_matcher is None:
        raise HTTPException(status_code=503, detail="Trial matcher not initialised")

    t0 = time.time()

    try:
        matches = trial_matcher.match_trials(
            cancer_type=req.cancer_type,
            biomarkers=req.biomarkers,
            stage=req.stage,
            age=req.age,
            top_k=req.top_k,
        )
    except Exception as exc:
        logger.error("Trial matching failed: %s", exc, exc_info=True)
        raise HTTPException(status_code=500, detail="Internal processing error")

    elapsed_ms = round((time.time() - t0) * 1000, 1)

    return {
        "matches": matches,
        "count": len(matches),
        "cancer_type": req.cancer_type,
        "processing_time_ms": elapsed_ms,
    }


@router.post("/api/trials/match-case/{case_id}")
async def match_trials_for_case(case_id: str, top_k: int = 10):
    """Match clinical trials for an existing case by extracting its profile."""
    from api.main import get_state

    state = get_state()
    case_manager = state.get("case_manager")
    trial_matcher = state.get("trial_matcher")

    if case_manager is None or trial_matcher is None:
        raise HTTPException(status_code=503, detail="Required services not available")

    # Retrieve case
    try:
        case = case_manager.get_case(case_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Case {case_id} not found")
    except Exception:
        raise HTTPException(status_code=500, detail="Internal processing error")

    if case is None:
        raise HTTPException(status_code=404, detail=f"Case {case_id} not found")

    t0 = time.time()

    try:
        matches = trial_matcher.match_trials(
            cancer_type=case.cancer_type,
            biomarkers=case.biomarkers or {},
            stage=case.stage,
            top_k=top_k,
        )
    except Exception as exc:
        logger.error(
            "Trial matching for case %s failed: %s", case_id, exc, exc_info=True
        )
        raise HTTPException(status_code=500, detail="Internal processing error")

    elapsed_ms = round((time.time() - t0) * 1000, 1)

    return {
        "case_id": case_id,
        "matches": matches,
        "count": len(matches),
        "processing_time_ms": elapsed_ms,
    }


@router.post("/api/therapies/rank")
async def rank_therapies(req: TherapyRankRequest):
    """Rank therapy options based on molecular profile and prior treatment."""
    from api.main import get_state

    state = get_state()
    therapy_ranker = state.get("therapy_ranker")
    if therapy_ranker is None:
        raise HTTPException(status_code=503, detail="Therapy ranker not initialised")

    t0 = time.time()

    try:
        ranked = therapy_ranker.rank_therapies(
            cancer_type=req.cancer_type,
            variants=req.variants,
            biomarkers=req.biomarkers,
            prior_therapies=req.prior_therapies,
        )
    except Exception as exc:
        logger.error("Therapy ranking failed: %s", exc, exc_info=True)
        raise HTTPException(status_code=500, detail="Internal processing error")

    elapsed_ms = round((time.time() - t0) * 1000, 1)

    return {
        "therapies": ranked,
        "count": len(ranked),
        "cancer_type": req.cancer_type,
        "processing_time_ms": elapsed_ms,
    }
