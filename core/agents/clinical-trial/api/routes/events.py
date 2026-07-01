"""Clinical trial SSE (Server-Sent Events) streaming routes.

Provides real-time event streaming for long-running clinical trial workflows,
patient matching progress, and protocol analysis. Uses FastAPI's
StreamingResponse with text/event-stream content type for SSE protocol
compliance.

Author: Adam Jones
Date: March 2026
"""

import asyncio
import json
import time
from collections import deque
from datetime import datetime, timezone
from typing import AsyncGenerator, List, Optional

from fastapi import APIRouter, Query, Request
from fastapi.responses import StreamingResponse
from loguru import logger

router = APIRouter(prefix="/v1/events", tags=["events"])


# =====================================================================
# Cross-Agent Event Bus
# =====================================================================

_VALID_EVENT_TYPES = {
    "patient_matched",
    "protocol_optimized",
    "safety_signal_detected",
    "workflow_complete",
    "critical_alert",
    "regulatory_update",
    "cross_modal_trigger",
}

_event_queue: deque = deque(maxlen=1000)
_event_subscribers: List[asyncio.Queue] = []
_event_lock = asyncio.Lock()


def publish_event(event_type: str, data: dict) -> dict:
    """Publish a cross-agent event to the event bus.

    This function is importable by other modules:
        from api.routes.events import publish_event
        publish_event("patient_matched", {"trial_id": "NCT001", "score": 0.92})

    Args:
        event_type: One of patient_matched, protocol_optimized,
                    safety_signal_detected, workflow_complete,
                    critical_alert, regulatory_update, cross_modal_trigger.
        data: Event-specific payload dict.

    Returns:
        The full event dict that was published.
    """
    if event_type not in _VALID_EVENT_TYPES:
        logger.warning(f"Unknown event_type '{event_type}'; publishing anyway")

    event = {
        "event_type": event_type,
        "agent": "clinical-trial-intelligence-agent",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "data": data,
    }
    _event_queue.append(event)

    # Notify all active SSE subscribers (non-blocking)
    stale: List[int] = []
    for i, q in enumerate(_event_subscribers):
        try:
            q.put_nowait(event)
        except asyncio.QueueFull:
            stale.append(i)
    # Remove stale subscribers whose queues are full
    for idx in reversed(stale):
        try:
            _event_subscribers.pop(idx)
        except IndexError:
            pass

    logger.debug(f"Published event: {event_type} -> {len(_event_subscribers)} subscribers")
    return event


async def _cross_agent_event_generator(
    req: Request,
    last_n: int = 0,
    max_duration: float = 300.0,
    heartbeat_interval: float = 15.0,
) -> AsyncGenerator[str, None]:
    """SSE generator that yields cross-agent events as they arrive."""
    subscriber_queue: asyncio.Queue = asyncio.Queue(maxsize=256)
    _event_subscribers.append(subscriber_queue)

    try:
        # Optionally replay the last N events from the buffer
        if last_n > 0:
            recent = list(_event_queue)[-last_n:]
            for evt in recent:
                yield _sse_message(evt["event_type"], evt, event_id=f"replay-{evt['timestamp']}")

        start = time.monotonic()
        hb_seq = 0
        while (time.monotonic() - start) < max_duration:
            if await req.is_disconnected():
                break
            try:
                event = await asyncio.wait_for(subscriber_queue.get(), timeout=heartbeat_interval)
                yield _sse_message(
                    event["event_type"],
                    event,
                    event_id=f"evt-{event['timestamp']}",
                )
            except asyncio.TimeoutError:
                hb_seq += 1
                yield _sse_message("heartbeat", {"seq": hb_seq, "agent": "clinical-trial"}, event_id=f"hb-{hb_seq}")

        yield _sse_message("done", {"reason": "max_duration_reached"})
    finally:
        try:
            _event_subscribers.remove(subscriber_queue)
        except ValueError:
            pass


# =====================================================================
# SSE Helpers
# =====================================================================

def _sse_message(event: str, data: dict, event_id: Optional[str] = None) -> str:
    """Format a single SSE message."""
    lines = []
    if event_id:
        lines.append(f"id: {event_id}")
    lines.append(f"event: {event}")
    lines.append(f"data: {json.dumps(data)}")
    lines.append("")  # Trailing blank line required by SSE spec
    return "\n".join(lines) + "\n"


async def _heartbeat_generator(
    interval: float = 15.0,
    max_duration: float = 300.0,
) -> AsyncGenerator[str, None]:
    """Yield SSE heartbeat pings to keep the connection alive."""
    start = time.monotonic()
    seq = 0
    while (time.monotonic() - start) < max_duration:
        seq += 1
        yield _sse_message("heartbeat", {"seq": seq, "agent": "clinical-trial"}, event_id=f"hb-{seq}")
        await asyncio.sleep(interval)
    yield _sse_message("done", {"reason": "max_duration_reached"})


async def _workflow_progress_generator(
    workflow_id: str,
    data: dict,
    req: Request,
) -> AsyncGenerator[str, None]:
    """Simulate workflow progress events for demo / UI integration."""
    steps = {
        "protocol_design": ["Analyzing protocol structure", "Scoring complexity", "Benchmarking against competitors", "Generating optimization recommendations", "Finalizing report"],
        "patient_matching": ["Loading patient profile", "Screening eligibility criteria", "Scoring inclusion criteria", "Checking exclusion criteria", "Ranking matches"],
        "site_selection": ["Evaluating site history", "Calculating enrollment rates", "Assessing diversity metrics", "Scoring feasibility", "Ranking sites"],
        "safety_signal": ["Coding adverse events", "Running disproportionality analysis", "Assessing causality", "Generating DSMB summary", "Finalizing signals"],
        "regulatory_docs": ["Selecting template", "Populating sections", "Checking compliance", "Formatting for agency", "Generating appendices"],
        "competitive_intel": ["Scanning trial registries", "Analyzing enrollment race", "Assessing threat levels", "Mapping landscape", "Summarizing intelligence"],
    }

    workflow_steps = steps.get(workflow_id, ["Processing step 1", "Processing step 2", "Processing step 3", "Finalizing"])
    total = len(workflow_steps)

    yield _sse_message("workflow_start", {
        "workflow_id": workflow_id,
        "total_steps": total,
        "status": "running",
    }, event_id=f"{workflow_id}-start")

    for i, step_name in enumerate(workflow_steps, 1):
        await asyncio.sleep(0.8)  # Simulate processing time
        yield _sse_message("workflow_progress", {
            "workflow_id": workflow_id,
            "step": i,
            "total_steps": total,
            "step_name": step_name,
            "progress_pct": round(i / total * 100, 1),
            "status": "running",
        }, event_id=f"{workflow_id}-step-{i}")

    await asyncio.sleep(0.5)
    yield _sse_message("workflow_complete", {
        "workflow_id": workflow_id,
        "status": "completed",
        "progress_pct": 100.0,
        "result_summary": f"{workflow_id} workflow completed successfully",
    }, event_id=f"{workflow_id}-done")

    yield _sse_message("done", {"workflow_id": workflow_id})


# =====================================================================
# Endpoints
# =====================================================================

@router.get("/stream")
async def event_stream(
    req: Request,
    workflow_id: Optional[str] = Query(None, description="Workflow to stream progress for"),
    heartbeat_only: bool = Query(False, description="Only send heartbeat pings"),
    cross_agent: bool = Query(False, description="Subscribe to cross-agent event bus"),
    last_n: int = Query(0, ge=0, le=1000, description="Replay last N events (cross_agent mode)"),
):
    """SSE endpoint for real-time workflow progress and cross-agent events.

    Connect with EventSource:
        const es = new EventSource('/v1/events/stream?workflow_id=patient_matching');
        es.addEventListener('workflow_progress', (e) => { ... });

    Cross-agent event bus:
        const es = new EventSource('/v1/events/stream?cross_agent=true&last_n=10');
        es.addEventListener('patient_matched', (e) => { ... });
        es.addEventListener('safety_signal_detected', (e) => { ... });
    """
    if cross_agent:
        generator = _cross_agent_event_generator(req, last_n=last_n)
    elif heartbeat_only:
        generator = _heartbeat_generator()
    elif workflow_id:
        generator = _workflow_progress_generator(workflow_id, {}, req)
    else:
        generator = _heartbeat_generator(interval=10.0, max_duration=60.0)

    return StreamingResponse(
        generator,
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no",
        },
    )


@router.get("/health")
async def events_health():
    """SSE subsystem health check."""
    return {
        "status": "healthy",
        "sse_enabled": True,
        "cross_agent_enabled": True,
        "event_queue_size": len(_event_queue),
        "active_subscribers": len(_event_subscribers),
        "supported_events": [
            "heartbeat",
            "workflow_start",
            "workflow_progress",
            "workflow_complete",
            "query_progress",
            "patient_matched",
            "done",
        ],
        "cross_agent_event_types": sorted(_VALID_EVENT_TYPES),
    }
