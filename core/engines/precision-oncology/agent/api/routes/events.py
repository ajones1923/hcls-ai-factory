"""
Oncology Intelligence Agent - Audit Event Log Router
===================================================
In-memory audit trail for case creation, MTB generation, report exports,
and other trackable actions within the oncology agent.

Author: Adam Jones
Date: February 2026
"""

import uuid
import logging
from datetime import datetime, timezone
from typing import Any, Dict, List

from fastapi import APIRouter, HTTPException, Query

logger = logging.getLogger(__name__)

router = APIRouter(tags=["events"])

# ---------------------------------------------------------------------------
# In-memory event store
# ---------------------------------------------------------------------------
_MAX_EVENTS = 500
_events: List[Dict[str, Any]] = []


def record_event(
    event_type: str,
    details: Dict[str, Any],
    user: str = "system",
) -> Dict[str, Any]:
    """
    Record an audit event.

    Called by other route handlers when significant actions occur (case
    creation, MTB generation, report export, etc.).
    """
    event = {
        "id": str(uuid.uuid4()),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "event_type": event_type,
        "details": details,
        "user": user,
    }
    _events.append(event)

    # Trim to max size
    if len(_events) > _MAX_EVENTS:
        del _events[: len(_events) - _MAX_EVENTS]

    logger.debug("Event recorded: %s [%s]", event_type, event["id"])
    return event


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
@router.get("/api/events")
async def list_events(
    limit: int = Query(default=50, ge=1, le=500),
    offset: int = Query(default=0, ge=0),
):
    """
    List audit events with pagination.

    Returns events in reverse-chronological order (newest first).
    """
    # Reverse so newest events come first
    ordered = list(reversed(_events))
    page = ordered[offset: offset + limit]

    return {
        "events": page,
        "total": len(_events),
        "limit": limit,
        "offset": offset,
    }


@router.get("/api/events/{event_id}")
async def get_event(event_id: str):
    """Retrieve a specific audit event by ID."""
    for evt in _events:
        if evt["id"] == event_id:
            return evt

    raise HTTPException(status_code=404, detail=f"Event {event_id} not found")
