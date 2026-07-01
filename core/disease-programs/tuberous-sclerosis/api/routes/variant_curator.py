"""
Variant-curator draft report + sign-off (PRD §3 FR-VC-5). The human gate is non-negotiable:
the engine produces an AI-labeled DRAFT; a board-certified molecular geneticist signs off.
"""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request
from pydantic import BaseModel

from src.agents.variant_curator.report import (
    get_signoff, record_signoff, render_draft, write_draft,
)

router = APIRouter(prefix="/variant-curator", tags=["variant-curator"])


class SignOff(BaseModel):
    signer: str
    decision: str = "approve"      # approve | reject | amend
    comment: str = ""


@router.get("/{patient_id}/draft")
def draft(patient_id: str, request: Request):
    vi = request.app.state.engine.store.projection(patient_id).get("variant_interp") or {}
    path = write_draft(patient_id, vi)
    return {"patient_id": patient_id, "path": str(path),
            "review_status": (get_signoff(patient_id) or {}).get("review_status",
                                                                  "pending_molecular_geneticist"),
            "markdown": render_draft(patient_id, vi)}


@router.post("/{patient_id}/sign-off")
def sign_off(patient_id: str, body: SignOff, request: Request):
    if body.decision not in ("approve", "reject", "amend"):
        raise HTTPException(400, "decision must be approve | reject | amend")
    # ensure a draft exists to sign
    vi = request.app.state.engine.store.projection(patient_id).get("variant_interp") or {}
    write_draft(patient_id, vi)
    return record_signoff(patient_id, body.signer, body.decision, body.comment)


@router.get("/{patient_id}/sign-off")
def sign_off_status(patient_id: str):
    return get_signoff(patient_id) or {"patient_id": patient_id,
                                       "review_status": "pending_molecular_geneticist"}
