"""Event endpoints (PRD §2.5.1)."""
from __future__ import annotations

from fastapi import APIRouter, Request

router = APIRouter(prefix="/events", tags=["events"])


@router.get("/{patient_id}")
def events_for(patient_id: str, request: Request):
    evs = request.app.state.engine.store.events_for(patient_id)
    return [{"event_id": e.event_id, "type": e.event_type.value,
             "created_at": e.created_at, "payload_keys": sorted((e.payload or {}).keys())}
            for e in evs]
