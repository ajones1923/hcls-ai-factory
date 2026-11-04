"""Rare disease diagnostic SSE (Server-Sent Events) streaming routes.

Provides real-time event streaming for long-running diagnostic workflows,
phenotype analysis progress, and variant interpretation. Uses FastAPI's
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

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import StreamingResponse
from loguru import logger

router = APIRouter(prefix="/v1/events", tags=["events"])

#: Attached to every frame the demo animation emits. The generator sleeps and
#: yields step labels; it calls no engine and computes nothing.
_DEMO_NOTE = (
    "DEMO ANIMATION — no workflow was executed, no collection was queried and "
    "no criterion was evaluated. This stream is not evidence that any step "
    "named in it was performed."
)


# =====================================================================
# Cross-Agent Event Bus
# =====================================================================

_VALID_EVENT_TYPES = {
    "diagnosis_complete",
    "variant_classified",
    "phenotype_matched",
    "therapy_found",
    "workflow_complete",
    "critical_alert",
    "cross_modal_trigger",
}

_event_queue: deque = deque(maxlen=1000)
_event_subscribers: List[asyncio.Queue] = []
_event_lock = asyncio.Lock()


def publish_event(event_type: str, data: dict) -> dict:
    """Publish a cross-agent event to the event bus.

    This function is importable by other modules:
        from api.routes.events import publish_event
        publish_event("diagnosis_complete", {"disease_id": "OMIM:219700", "confidence": 0.92})

    Args:
        event_type: One of diagnosis_complete, variant_classified,
                    phenotype_matched, therapy_found, workflow_complete,
                    critical_alert, cross_modal_trigger.
        data: Event-specific payload dict.

    Returns:
        The full event dict that was published.
    """
    if event_type not in _VALID_EVENT_TYPES:
        logger.warning(f"Unknown event_type '{event_type}'; publishing anyway")

    event = {
        "event_type": event_type,
        "agent": "rare-disease-diagnostic-agent",
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
                yield _sse_message("heartbeat", {"seq": hb_seq, "agent": "rare-disease"}, event_id=f"hb-{hb_seq}")

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
        yield _sse_message("heartbeat", {"seq": seq, "agent": "rare-disease"}, event_id=f"hb-{seq}")
        await asyncio.sleep(interval)
    yield _sse_message("done", {"reason": "max_duration_reached"})


async def _demo_progress_generator(
    workflow_id: str,
    data: dict,
    req: Request,
) -> AsyncGenerator[str, None]:
    """Emit a DEMO ANIMATION. No workflow is executed and nothing is computed.

    This generator sleeps and yields step labels. It does not call any engine,
    query any collection, run any predictor, or evaluate any criterion.

    It used to emit the step names below under the event names
    ``workflow_progress`` / ``workflow_complete``, with ``"status":
    "running"``, then ``"status": "completed"`` and ``"result_summary":
    "<id> workflow completed successfully"``. The word "Simulate" appeared
    only in this docstring and an inline comment -- no part of the wire payload
    said so, and the public endpoint docstring called it "real-time workflow
    progress".

    These streams were the only place in the product that CLAIMED the
    evidence work was done, and the one place guaranteed not to have done it.
    A clinician watching six steps tick to 100% ending "variant_interpretation
    workflow completed successfully" had positive, specific, step-by-step
    assurance that ClinVar was checked, population frequency evaluated,
    in-silico predictors run and ACMG criteria applied -- the four things they
    would otherwise verify by hand. A coordinator watching "Screening
    eligibility criteria -> Scoring inclusion criteria -> Checking exclusion
    criteria -> Ranking matches" concluded the patient's criteria had been
    checked against the trials shown. For ``competitive_intel`` the step
    "Scanning trial registries" asserted a ClinicalTrials.gov query that never
    happened.

    Every payload now carries ``"simulated": true`` and a note, the events are
    named ``demo_progress`` / ``demo_complete``, the terminal status is
    ``"simulated"`` rather than ``"completed"``, and the endpoint refuses to
    serve the animation unless the caller explicitly passes ``demo=true``.
    """
    steps = {
        "differential_diagnosis": ["Parsing phenotype terms", "Querying HPO ontology", "Matching against OMIM", "Matching against Orphanet", "Ranking differential candidates", "Generating recommendations"],
        "variant_interpretation": ["Loading variant data", "Checking ClinVar records", "Evaluating population frequency", "Running in-silico predictors", "Applying ACMG criteria", "Final classification"],
        "phenotype_matching": ["Encoding HPO terms", "Computing IC-weighted similarity", "Screening disease profiles", "Ranking matches", "Generating overlap report"],
        "therapeutic_search": ["Searching orphan drug databases", "Checking gene therapy eligibility", "Querying clinical trials", "Evaluating ERT options", "Compiling recommendations"],
        "gene_panel_analysis": ["Selecting virtual gene panel", "Filtering variants by panel genes", "Prioritizing by pathogenicity", "Checking inheritance patterns", "Generating panel report"],
        "case_similarity": ["Encoding patient phenotype", "Searching case report database", "Computing similarity scores", "Ranking similar cases", "Extracting diagnostic insights"],
    }

    workflow_steps = steps.get(workflow_id, ["Processing step 1", "Processing step 2", "Processing step 3", "Finalizing"])
    total = len(workflow_steps)

    yield _sse_message("demo_start", {
        "workflow_id": workflow_id,
        "total_steps": total,
        "status": "simulated",
        "simulated": True,
        "note": _DEMO_NOTE,
    }, event_id=f"{workflow_id}-start")

    for i, step_name in enumerate(workflow_steps, 1):
        await asyncio.sleep(0.8)  # a sleep, not processing
        yield _sse_message("demo_progress", {
            "workflow_id": workflow_id,
            "step": i,
            "total_steps": total,
            "step_name": f"[SIMULATED] {step_name}",
            "progress_pct": round(i / total * 100, 1),
            "status": "simulated",
            "simulated": True,
            "note": _DEMO_NOTE,
        }, event_id=f"{workflow_id}-step-{i}")

    await asyncio.sleep(0.5)
    yield _sse_message("demo_complete", {
        "workflow_id": workflow_id,
        "status": "simulated",
        "progress_pct": 100.0,
        "simulated": True,
        "note": _DEMO_NOTE,
        "result_summary": (
            f"Demo animation for '{workflow_id}' finished. NOTHING WAS "
            f"COMPUTED: no engine was called, no collection was queried and "
            f"no criterion was evaluated. Do not read this as a completed "
            f"workflow."
        ),
    }, event_id=f"{workflow_id}-done")

    yield _sse_message("done", {
        "workflow_id": workflow_id, "simulated": True,
    })


# =====================================================================
# Endpoints
# =====================================================================

@router.get("/stream")
async def event_stream(
    req: Request,
    workflow_id: Optional[str] = Query(None, description="Workflow to stream progress for"),
    heartbeat_only: bool = Query(False, description="Only send heartbeat pings"),
    cross_agent: bool = Query(False, description="Subscribe to cross-agent event bus"),
    demo: bool = Query(
        False,
        description=(
            "Required to receive the canned demo animation for a workflow_id. "
            "The animation executes no workflow and computes nothing."
        ),
    ),
    last_n: int = Query(0, ge=0, le=1000, description="Replay last N events (cross_agent mode)"),
):
    """SSE endpoint for heartbeats, cross-agent events, and a DEMO animation.

    There is no live progress stream for a workflow. Passing a ``workflow_id``
    without ``demo=true`` returns HTTP 400; with ``demo=true`` it returns a
    canned animation whose every frame is marked ``"simulated": true``.
    Nothing is computed by it.

    Demo animation (nothing is executed):
        const es = new EventSource(
            '/v1/events/stream?workflow_id=differential_diagnosis&demo=true');
        es.addEventListener('demo_progress', (e) => { ... });

    Cross-agent event bus:
        const es = new EventSource('/v1/events/stream?cross_agent=true&last_n=10');
        es.addEventListener('diagnosis_complete', (e) => { ... });
        es.addEventListener('variant_classified', (e) => { ... });
    """
    if cross_agent:
        generator = _cross_agent_event_generator(req, last_n=last_n)
    elif heartbeat_only:
        generator = _heartbeat_generator()
    elif workflow_id:
        # The animation is opt-in. A caller that wants progress for real work
        # has to ask for a simulation by name; it is no longer what you get by
        # default from a workflow_id. Previously any `?workflow_id=` value --
        # unvalidated, so any string at all -- produced a step-by-step success
        # stream for work that never ran.
        if not demo:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"No live progress stream exists for workflow "
                    f"'{workflow_id}'. This endpoint can only replay a "
                    f"canned demo animation, which computes nothing; pass "
                    f"demo=true to receive it. To execute the workflow, POST "
                    f"to the workflow endpoint instead."
                ),
            )
        generator = _demo_progress_generator(workflow_id, {}, req)
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
            "diagnosis_complete",
            "variant_classified",
            "done",
        ],
        "cross_agent_event_types": sorted(_VALID_EVENT_TYPES),
    }
