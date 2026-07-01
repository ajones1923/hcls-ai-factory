"""Provenance endpoint (PRD §2.5.1, NFR-R-1) — full audit trail per patient."""
from __future__ import annotations

from fastapi import APIRouter, Request

router = APIRouter(prefix="/provenance", tags=["provenance"])


@router.get("/{patient_id}")
def provenance(patient_id: str, request: Request):
    out = []
    for e in request.app.state.engine.store.events_for(patient_id):
        recs = (e.provenance or {}).get("records")
        if recs:
            out.append({"event": e.event_type.value, "event_id": e.event_id, "records": recs})
    return out
