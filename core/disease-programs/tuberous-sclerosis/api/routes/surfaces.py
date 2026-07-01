"""Clinician-surface endpoints (PRD §7) — assembled from projections."""
from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request

router = APIRouter(prefix="/surfaces", tags=["surfaces"])
_KINDS = {"briefing", "dashboard", "alerts"}


@router.get("/{kind}/{patient_id}")
def surface(kind: str, patient_id: str, request: Request):
    if kind not in _KINDS:
        raise HTTPException(404, f"unknown surface '{kind}' (expected one of {sorted(_KINDS)})")
    return request.app.state.engine.assemble_surface(patient_id, kind)
