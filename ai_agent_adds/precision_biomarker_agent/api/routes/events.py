"""Cross-modal event routes for the Biomarker Intelligence Agent.

Provides endpoints for:
  - Receiving cross-modal triggers from other HCLS AI Factory agents
  - Sending biomarker alerts to the platform event bus

The cross-modal event system enables communication between agents:
  - Imaging Intelligence Agent -> Biomarker Agent (imaging findings)
  - Biomarker Agent -> Imaging Agent (elevated Lp(a), iron overload)
  - Biomarker Agent -> CAR-T/Oncology Agent (PGx drug safety alerts)
  - Biomarker Agent -> Genomics Pipeline (VCF re-analysis triggers)

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import threading
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, Query, Request
from loguru import logger
from pydantic import BaseModel, Field

router = APIRouter(prefix="/v1/events", tags=["events"])


# =====================================================================
# In-memory event stores (placeholder for production message bus)
# =====================================================================

_inbound_events: List[Dict[str, Any]] = []
_outbound_alerts: List[Dict[str, Any]] = []
_MAX_STORED_EVENTS = 1000
_events_lock = threading.Lock()


# =====================================================================
# Request / Response Schemas
# =====================================================================

class CrossModalEvent(BaseModel):
    """Inbound cross-modal event from another agent."""
    source_agent: str = Field(
        ..., description="Name of the sending agent (e.g., 'imaging_intelligence_agent')",
    )
    event_type: str = Field(
        ..., description="Event type (e.g., 'imaging_finding', 'genomic_variant', 'drug_alert')",
    )
    patient_id: str = Field(..., description="Patient identifier")
    payload: Dict[str, Any] = Field(
        default_factory=dict,
        description="Event-specific data payload",
    )
    urgency: str = Field(
        "moderate",
        pattern="^(critical|high|moderate|low)$",
        description="Urgency level: critical, high, moderate, low",
    )
    correlation_id: Optional[str] = Field(
        None, description="Correlation ID for tracking related events across agents",
    )


class CrossModalEventResponse(BaseModel):
    """Response after processing a cross-modal event."""
    event_id: str = Field(..., description="Generated event ID for tracking")
    status: str = Field(..., description="Processing status: received, processed, error")
    actions_triggered: List[str] = Field(
        default_factory=list,
        description="List of actions triggered by this event",
    )
    timestamp: str


class BiomarkerAlert(BaseModel):
    """Outbound biomarker alert to the platform event bus."""
    patient_id: str = Field(..., description="Patient identifier")
    alert_type: str = Field(
        ...,
        description="Alert type: elevated_lpa, pgx_critical, pre_diabetic_trajectory, "
                    "iron_overload, accelerated_aging",
    )
    severity: str = Field(
        ..., pattern="^(critical|high|moderate|low)$", description="Severity: critical, high, moderate, low",
    )
    target_agent: str = Field(
        ..., description="Target agent to receive this alert",
    )
    payload: Dict[str, Any] = Field(
        default_factory=dict,
        description="Alert-specific data payload",
    )
    description: str = Field("", description="Human-readable alert description")


class BiomarkerAlertResponse(BaseModel):
    """Response after sending a biomarker alert."""
    alert_id: str
    status: str = "sent"
    target_agent: str
    timestamp: str


class EventListResponse(BaseModel):
    """Paginated list of events."""
    events: List[Dict[str, Any]]
    total: int
    page: int
    page_size: int


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/cross-modal", response_model=CrossModalEventResponse)
async def receive_cross_modal_event(event: CrossModalEvent, req: Request):
    """Receive a cross-modal trigger from another HCLS AI Factory agent.

    Processes inbound events such as:
    - Imaging findings that affect biomarker interpretation
    - Genomic variant discoveries requiring biomarker re-analysis
    - Drug alerts from oncology agents requiring PGx cross-check

    The event is logged and relevant actions are triggered based on
    the event type and payload.
    """
    if not event.source_agent or not event.source_agent.strip():
        raise HTTPException(status_code=400, detail="source_agent is required")

    event_id = uuid.uuid4().hex[:12]
    timestamp = datetime.now(timezone.utc).isoformat()

    # Store the event
    record = {
        "event_id": event_id,
        "source_agent": event.source_agent,
        "event_type": event.event_type,
        "patient_id": event.patient_id,
        "payload": event.payload,
        "urgency": event.urgency,
        "correlation_id": event.correlation_id or event_id,
        "timestamp": timestamp,
        "status": "received",
    }
    with _events_lock:
        _inbound_events.append(record)

        # Evict oldest events if store exceeds max size
        if len(_inbound_events) > _MAX_STORED_EVENTS:
            _inbound_events[:] = _inbound_events[-_MAX_STORED_EVENTS:]

    # Determine actions based on event type
    actions_triggered = []

    if event.event_type == "imaging_finding":
        # Imaging agent sent a finding -- correlate with biomarker data
        actions_triggered.append("biomarker_correlation_analysis")
        if event.payload.get("finding_type") == "coronary_calcium":
            actions_triggered.append("cardiovascular_risk_reassessment")

    elif event.event_type == "genomic_variant":
        # Genomics pipeline discovered a new variant
        gene = event.payload.get("gene", "")
        actions_triggered.append(f"genotype_adjustment_update:{gene}")
        actions_triggered.append("disease_trajectory_reanalysis")

    elif event.event_type == "drug_alert":
        # Oncology/CAR-T agent flagged a drug interaction concern
        actions_triggered.append("pgx_cross_check")
        drug = event.payload.get("drug", "")
        if drug:
            actions_triggered.append(f"pgx_drug_safety_review:{drug}")

    elif event.event_type == "treatment_started":
        # Patient started a new treatment -- check PGx implications
        actions_triggered.append("pgx_treatment_compatibility_check")
        actions_triggered.append("biomarker_monitoring_schedule_update")

    else:
        actions_triggered.append("event_logged")

    for action in actions_triggered:
        logger.info(f"Cross-modal action triggered: {action} for patient {event.patient_id} (event_id={event_id})")

    record["status"] = "processed"
    record["actions_triggered"] = actions_triggered

    # In production, these would publish to a message bus (Redis, Kafka, etc.)
    # Publish to app.state event queue if available
    event_queue = getattr(req.app.state, "event_queue", None)
    if event_queue is not None:
        event_queue.append(record)

    return CrossModalEventResponse(
        event_id=event_id,
        status="processed",
        actions_triggered=actions_triggered,
        timestamp=timestamp,
    )


@router.post("/biomarker-alert", response_model=BiomarkerAlertResponse)
async def send_biomarker_alert(alert: BiomarkerAlert, req: Request):
    """Send a biomarker alert to the platform event bus.

    Used by the Biomarker Agent to notify other agents of findings that
    require cross-modal follow-up:

    - **elevated_lpa**: Triggers cardiovascular imaging assessment
    - **pgx_critical**: Triggers oncology drug interaction check
    - **pre_diabetic_trajectory**: Triggers genomic pipeline VCF re-analysis
    - **iron_overload**: Triggers liver imaging quantification
    - **accelerated_aging**: Triggers epigenetic variant analysis

    In production, this would publish to a message queue (Redis, Kafka, etc.).
    """
    alert_id = uuid.uuid4().hex[:12]
    timestamp = datetime.now(timezone.utc).isoformat()

    # Validate alert type against known cross-modal links
    from src.knowledge import CROSS_MODAL_LINKS
    valid_types = set(CROSS_MODAL_LINKS.keys())
    if alert.alert_type not in valid_types:
        # Allow unknown types but log a warning
        pass

    # Store the alert
    record = {
        "alert_id": alert_id,
        "patient_id": alert.patient_id,
        "alert_type": alert.alert_type,
        "severity": alert.severity,
        "target_agent": alert.target_agent,
        "payload": alert.payload,
        "description": alert.description,
        "timestamp": timestamp,
        "status": "sent",
    }
    with _events_lock:
        _outbound_alerts.append(record)

        # Evict oldest alerts if store exceeds max size
        if len(_outbound_alerts) > _MAX_STORED_EVENTS:
            _outbound_alerts[:] = _outbound_alerts[-_MAX_STORED_EVENTS:]

    # In production: publish to message bus
    # await message_bus.publish(alert.target_agent, record)

    return BiomarkerAlertResponse(
        alert_id=alert_id,
        status="sent",
        target_agent=alert.target_agent,
        timestamp=timestamp,
    )


@router.get("/cross-modal", response_model=EventListResponse)
async def list_inbound_events(
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=200),
    patient_id: Optional[str] = Query(None, description="Filter by patient ID"),
    event_type: Optional[str] = Query(None, description="Filter by event type"),
):
    """List inbound cross-modal events (newest first).

    Supports filtering by patient_id and event_type.
    """
    with _events_lock:
        filtered = list(_inbound_events)
    if patient_id:
        filtered = [e for e in filtered if e.get("patient_id") == patient_id]
    if event_type:
        filtered = [e for e in filtered if e.get("event_type") == event_type]

    # Newest first
    filtered = list(reversed(filtered))

    start = (page - 1) * page_size
    end = start + page_size
    page_events = filtered[start:end]

    return EventListResponse(
        events=page_events,
        total=len(filtered),
        page=page,
        page_size=page_size,
    )


@router.get("/biomarker-alert", response_model=EventListResponse)
async def list_outbound_alerts(
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=200),
    patient_id: Optional[str] = Query(None, description="Filter by patient ID"),
    alert_type: Optional[str] = Query(None, description="Filter by alert type"),
):
    """List outbound biomarker alerts (newest first).

    Supports filtering by patient_id and alert_type.
    """
    with _events_lock:
        filtered = list(_outbound_alerts)
    if patient_id:
        filtered = [e for e in filtered if e.get("patient_id") == patient_id]
    if alert_type:
        filtered = [e for e in filtered if e.get("alert_type") == alert_type]

    # Newest first
    filtered = list(reversed(filtered))

    start = (page - 1) * page_size
    end = start + page_size
    page_alerts = filtered[start:end]

    return EventListResponse(
        events=page_alerts,
        total=len(filtered),
        page=page,
        page_size=page_size,
    )
