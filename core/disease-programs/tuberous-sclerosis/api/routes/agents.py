"""Agent listing + projection access (PRD §2.3; debug/demo control)."""
from __future__ import annotations

from fastapi import APIRouter, Request

router = APIRouter(prefix="/agents", tags=["agents"])


@router.get("")
def list_agents(request: Request):
    eng = request.app.state.engine
    from src.orchestrator.events import AGENT_DEPENDS_ON, AGENT_EMITS
    return [{"name": a.name, "emits": AGENT_EMITS[a.name].value,
             "depends_on": [d.value for d in AGENT_DEPENDS_ON[a.name]]} for a in eng.agents]


@router.get("/projection/{patient_id}")
def projection(patient_id: str, request: Request):
    return request.app.state.engine.store.projection(patient_id)
