"""
Oncology Intelligence Agent - Case Management Router
===================================================
Create patient cases, retrieve case details, list variants, and generate
Molecular Tumor Board (MTB) packets for clinical review.

Author: Adam Jones
Date: February 2026
"""

import time
import uuid
import logging
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Union

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

logger = logging.getLogger(__name__)

router = APIRouter(tags=["cases"])


# ---------------------------------------------------------------------------
# Schemas
# ---------------------------------------------------------------------------
class VariantInput(BaseModel):
    gene: str
    variant: str
    variant_type: str = "SNV"
    vaf: Optional[float] = None
    consequence: Optional[str] = None


class CreateCaseRequest(BaseModel):
    patient_id: str = Field(..., min_length=1)
    cancer_type: str
    stage: Optional[str] = None
    variants: Optional[Union[List[VariantInput], List[Dict[str, Any]]]] = None
    vcf_text: Optional[str] = Field(
        None, description="Raw VCF content (alternative to variant list)"
    )
    biomarkers: Optional[Dict[str, Any]] = Field(default_factory=dict)
    prior_therapies: Optional[List[str]] = Field(default_factory=list)


class CreateCaseResponse(BaseModel):
    case_id: str
    patient_id: str
    cancer_type: str
    stage: Optional[str] = None
    variant_count: int = 0
    biomarkers: Dict[str, Any] = Field(default_factory=dict)
    prior_therapies: List[str] = Field(default_factory=list)
    created_at: str


class MTBRequest(BaseModel):
    include_trials: bool = True
    include_therapies: bool = True
    include_resistance: bool = True
    top_k: int = Field(default=10, ge=1, le=50)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _record_event(event_type: str, details: Dict[str, Any]):
    """Append an event to the in-memory event log (from events router)."""
    try:
        from api.routes.events import record_event

        record_event(event_type, details)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@router.post("/api/cases", response_model=CreateCaseResponse)
async def create_case(req: CreateCaseRequest):
    """Create a new oncology case for MTB analysis."""
    from api.main import get_state

    state = get_state()
    case_manager = state.get("case_manager")
    if case_manager is None:
        raise HTTPException(status_code=503, detail="Case manager not initialised")

    t0 = time.time()

    # Normalise variant input
    variants_raw: List[Dict[str, Any]] = []
    if req.variants:
        for v in req.variants:
            if isinstance(v, VariantInput):
                variants_raw.append(v.model_dump())
            elif isinstance(v, dict):
                variants_raw.append(v)

    try:
        # Pass VCF text if provided, otherwise pass parsed variant list
        vcf_or_variants = req.vcf_text if req.vcf_text else variants_raw
        case = case_manager.create_case(
            patient_id=req.patient_id,
            cancer_type=req.cancer_type,
            stage=req.stage or "",
            vcf_content_or_variants=vcf_or_variants,
            biomarkers=req.biomarkers or {},
            prior_therapies=req.prior_therapies or [],
        )
    except Exception as exc:
        logger.error("Failed to create case: %s", exc, exc_info=True)
        raise HTTPException(status_code=500, detail="Internal processing error")

    elapsed_ms = round((time.time() - t0) * 1000, 1)

    # CaseSnapshot is a Pydantic model — access attributes directly
    case_id = getattr(case, "case_id", None) or str(uuid.uuid4())
    variant_count = len(getattr(case, "variants", []))
    created_at = getattr(case, "created_at", None) or datetime.now(timezone.utc).isoformat()

    _record_event(
        "case_created",
        {
            "case_id": case_id,
            "patient_id": req.patient_id,
            "cancer_type": req.cancer_type,
            "processing_time_ms": elapsed_ms,
        },
    )

    return CreateCaseResponse(
        case_id=case_id,
        patient_id=req.patient_id,
        cancer_type=req.cancer_type,
        stage=req.stage,
        variant_count=variant_count,
        biomarkers=req.biomarkers or {},
        prior_therapies=req.prior_therapies or [],
        created_at=created_at,
    )


@router.get("/api/cases/{case_id}")
async def get_case(case_id: str):
    """Retrieve an existing oncology case."""
    from api.main import get_state

    state = get_state()
    case_manager = state.get("case_manager")
    if case_manager is None:
        raise HTTPException(status_code=503, detail="Case manager not initialised")

    try:
        case = await case_manager.get_case(case_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Case {case_id} not found")
    except Exception as exc:
        logger.error("Error retrieving case %s: %s", case_id, exc, exc_info=True)
        raise HTTPException(status_code=500, detail="Internal processing error")

    return case


@router.post("/api/cases/{case_id}/mtb")
async def generate_mtb_packet(case_id: str, req: Optional[MTBRequest] = None):
    """Generate a Molecular Tumor Board packet for a case."""
    from api.main import get_state

    state = get_state()
    case_manager = state.get("case_manager")
    if case_manager is None:
        raise HTTPException(status_code=503, detail="Case manager not initialised")

    if req is None:
        req = MTBRequest()

    t0 = time.time()

    try:
        mtb_packet = await case_manager.generate_mtb_packet(
            case_id=case_id,
            include_trials=req.include_trials,
            include_therapies=req.include_therapies,
            include_resistance=req.include_resistance,
            top_k=req.top_k,
        )
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Case {case_id} not found")
    except Exception as exc:
        logger.error(
            "Failed to generate MTB packet for %s: %s", case_id, exc, exc_info=True
        )
        raise HTTPException(status_code=500, detail="Internal processing error")

    elapsed_ms = round((time.time() - t0) * 1000, 1)

    _record_event(
        "mtb_generated",
        {
            "case_id": case_id,
            "processing_time_ms": elapsed_ms,
            "variant_count": mtb_packet.get("variant_count", 0),
            "therapy_count": len(mtb_packet.get("therapy_ranking", [])),
            "trial_count": len(mtb_packet.get("trial_matches", [])),
        },
    )

    mtb_packet["processing_time_ms"] = elapsed_ms
    return mtb_packet


@router.get("/api/cases/{case_id}/variants")
async def list_variants(case_id: str):
    """List all variants associated with a case."""
    from api.main import get_state

    state = get_state()
    case_manager = state.get("case_manager")
    if case_manager is None:
        raise HTTPException(status_code=503, detail="Case manager not initialised")

    try:
        case = await case_manager.get_case(case_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"Case {case_id} not found")
    except Exception:
        raise HTTPException(status_code=500, detail="Internal processing error")

    variants = case.get("variants", [])
    return {
        "case_id": case_id,
        "variants": variants,
        "count": len(variants),
    }
