"""Pipeline event / audit-trail routes.

Provides endpoints to query the event log for pipeline activity,
enabling audit trails and operational visibility for the
Pharmacogenomics Intelligence Agent.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import time
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field

router = APIRouter(prefix="/api", tags=["events"])


# -- In-memory event store (placeholder for production persistence) ----
_event_store: List[Dict[str, Any]] = []


# -- Schemas -----------------------------------------------------------

class PipelineEvent(BaseModel):
    """A single pipeline event record."""

    event_id: str = Field(..., description="Unique event identifier")
    event_type: str = Field(..., description="Event category (e.g. query, ingest, report, drug_check, hla_screen)")
    timestamp: str = Field(..., description="ISO-8601 timestamp")
    source: str = Field("", description="Originating service or module")
    summary: str = Field("", description="Human-readable summary")
    metadata: Dict[str, Any] = Field(default_factory=dict)


class EventListResponse(BaseModel):
    """Response for GET /api/events."""

    events: List[PipelineEvent]
    total: int
    page: int
    page_size: int


# -- Public API for other modules to emit events ----------------------

def emit_event(
    event_type: str,
    source: str = "",
    summary: str = "",
    metadata: Optional[Dict[str, Any]] = None,
) -> str:
    """Record a pipeline event and return its event_id.

    Other modules (RAG engine, ingest, report generator, clinical
    endpoints) call this to append entries to the audit trail.

    Args:
        event_type: Category label (e.g. ``"query"``, ``"drug_check"``,
            ``"hla_screen"``, ``"phenoconversion"``, ``"report"``).
        source: Originating module or service name.
        summary: Human-readable description.
        metadata: Arbitrary key-value pairs.

    Returns:
        The generated ``event_id``.
    """
    event_id = uuid.uuid4().hex[:12]
    record = {
        "event_id": event_id,
        "event_type": event_type,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "source": source,
        "summary": summary,
        "metadata": metadata or {},
    }
    _event_store.append(record)
    return event_id


# -- Endpoints ---------------------------------------------------------

@router.get("/events", response_model=EventListResponse)
async def list_events(
    page: int = Query(1, ge=1, description="Page number"),
    page_size: int = Query(50, ge=1, le=200, description="Results per page"),
    event_type: Optional[str] = Query(None, description="Filter by event type"),
):
    """Return paginated pipeline events (audit trail).

    Supports optional filtering by ``event_type``.  Common types for
    the PGx agent include ``query``, ``drug_check``, ``hla_screen``,
    ``phenoconversion``, ``dosing``, ``report``, and ``ingest``.
    """
    filtered = _event_store
    if event_type:
        filtered = [e for e in filtered if e["event_type"] == event_type]

    # Newest first
    filtered = list(reversed(filtered))

    start = (page - 1) * page_size
    end = start + page_size
    page_events = filtered[start:end]

    return EventListResponse(
        events=[PipelineEvent(**e) for e in page_events],
        total=len(filtered),
        page=page,
        page_size=page_size,
    )


@router.get("/events/{event_id}", response_model=PipelineEvent)
async def get_event(event_id: str):
    """Return details for a specific pipeline event."""
    for record in _event_store:
        if record["event_id"] == event_id:
            return PipelineEvent(**record)

    raise HTTPException(status_code=404, detail=f"Event '{event_id}' not found")
