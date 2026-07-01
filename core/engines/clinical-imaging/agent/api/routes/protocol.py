"""Protocol optimization endpoints.

Provides evidence-based imaging protocol recommendations using ACR
Appropriateness Criteria with patient-specific adjustments.
"""

from fastapi import APIRouter
from pydantic import BaseModel
from typing import Any, Dict, List, Optional

from src.protocol_optimizer import PatientFactors, ProtocolOptimizer, ProtocolRecommendation

router = APIRouter(prefix="/protocol", tags=["Protocol Optimizer"])

# Singleton optimizer instance
_optimizer = ProtocolOptimizer()


# ── Request / Response Models ──────────────────────────────────────────


class ProtocolRequest(BaseModel):
    """Request body for protocol recommendation."""

    indication: str
    patient: Optional[PatientFactors] = None


class ProtocolResponse(BaseModel):
    """Response body for protocol recommendation."""

    indication: str
    recommended_modality: str
    recommended_protocol: str
    acr_appropriateness_rating: int
    parameters: Dict[str, Any] = {}
    contrast: Optional[Dict[str, Any]] = None
    dose_estimate_msv: Optional[float] = None
    warnings: List[str] = []
    alternatives: List[Dict[str, Any]] = []
    rationale: str = ""


# ── Endpoints ──────────────────────────────────────────────────────────


@router.post("/recommend", response_model=ProtocolResponse)
def recommend_protocol(req: ProtocolRequest):
    """Recommend optimal imaging protocol based on clinical indication and patient parameters.

    Accepts a clinical indication string and optional patient factors.
    Returns the top recommendation with ACR rating, dose estimate,
    patient-specific warnings, and ranked alternatives.
    """
    rec = _optimizer.recommend(
        indication=req.indication,
        patient=req.patient,
    )
    return ProtocolResponse(**rec.model_dump())


@router.get("/indications", response_model=List[str])
def list_indications():
    """List all available clinical indications supported by the optimizer."""
    return _optimizer.get_available_indications()
