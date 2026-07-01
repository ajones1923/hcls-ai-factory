"""
Event Bus for the HCLS AI Factory Pipeline.

Event-driven coordination across the three-stage precision medicine pipeline
(Genomics -> RAG/Chat -> Drug Discovery) and the multi-agent constellation
(Biomarker, Oncology, Imaging, CAR-T, Autoimmune).
Implements Section 11.0.3 of the architecture specification.

Design principles:
    - File-based persistence (JSON lines) for audit trail -- no PostgreSQL
      dependency; this platform uses the local filesystem for pipeline state.
    - Priority queue (heapq) for event ordering.
    - Singleton pattern for process-wide event coordination.
    - Prometheus metrics for events emitted and processed.
    - Sync publish for immediate handlers; async publish via background thread.
    - Event replay from audit log for crash recovery.

Hardware target: NVIDIA DGX Spark (GB10 GPU, 128 GB unified memory, $4,699).
"""

from __future__ import annotations

import heapq
import json
import os
import threading
import time
import uuid
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime
from enum import Enum, unique
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

from loguru import logger

# Re-use PipelineStage from enums if available; define locally otherwise.
try:
    from hcls_common.enums import PipelineStage
except ImportError:

    @unique
    class PipelineStage(str, Enum):  # type: ignore[no-redef]
        GENOMICS = "genomics"
        ANNOTATION = "annotation"
        RAG_INGEST = "rag_ingest"
        RAG_CHAT = "rag_chat"
        DRUG_DISCOVERY = "drug_discovery"
        CART_ANALYSIS = "cart_analysis"
        REPORTING = "reporting"


# ---------------------------------------------------------------------------
# Prometheus metrics (optional)
# ---------------------------------------------------------------------------

try:
    from prometheus_client import Counter, Histogram

    EVENTS_EMITTED = Counter(
        "hcls_events_emitted_total",
        "Total pipeline events emitted",
        ["event_type", "source_stage"],
    )
    EVENTS_PROCESSED = Counter(
        "hcls_events_processed_total",
        "Total pipeline events processed by handlers",
        ["event_type", "handler"],
    )
    EVENT_LATENCY = Histogram(
        "hcls_event_processing_seconds",
        "Event handler execution latency",
        ["event_type"],
        buckets=[0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 30.0],
    )
except ImportError:
    EVENTS_EMITTED = None
    EVENTS_PROCESSED = None
    EVENT_LATENCY = None


# ---------------------------------------------------------------------------
# Event types
# ---------------------------------------------------------------------------

@unique
class EventType(str, Enum):
    """Pipeline event types for inter-stage coordination."""

    # Stage 1: Genomics
    FASTQ_ARRIVED = "fastq_arrived"
    VCF_READY = "vcf_ready"

    # Stage 2: RAG / Annotation
    ANNOTATION_COMPLETE = "annotation_complete"
    TARGETS_IDENTIFIED = "targets_identified"

    # Stage 3: Drug Discovery
    DRUG_CANDIDATES_READY = "drug_candidates_ready"

    # CAR-T Intelligence Agent
    CART_ANALYSIS_COMPLETE = "cart_analysis_complete"
    CART_MANUFACTURING_READY = "cart_manufacturing_ready"

    # Precision Biomarker Agent
    BIOMARKER_ALERT = "biomarker_alert"
    BIOLOGICAL_AGE_COMPUTED = "biological_age_computed"
    PGX_RESULT_READY = "pgx_result_ready"
    DISEASE_TRAJECTORY_ALERT = "disease_trajectory_alert"

    # Precision Oncology Agent
    ONCOLOGY_CASE_CREATED = "oncology_case_created"
    THERAPY_RANKED = "therapy_ranked"
    TRIAL_MATCHED = "trial_matched"

    # Imaging Intelligence Agent
    IMAGING_FINDING = "imaging_finding"
    DICOM_INGESTED = "dicom_ingested"

    # Precision Autoimmune Agent
    AUTOIMMUNE_FLARE_PREDICTED = "autoimmune_flare_predicted"
    AUTOANTIBODY_PANEL_READY = "autoantibody_panel_ready"

    # Cross-agent coordination
    PGX_DRUG_FILTER = "pgx_drug_filter"
    CONCORDANCE_CHECK = "concordance_check"

    # Reporting
    REPORT_GENERATED = "report_generated"

    # Cross-stage requests (bidirectional triggers)
    DEEP_ANALYSIS_REQUEST = "deep_analysis_request"
    VARIANT_CONTEXT_REQUEST = "variant_context_request"
    DRUG_DESIGN_REQUEST = "drug_design_request"

    @property
    def display_name(self) -> str:
        return self.value.replace("_", " ").title()


# ---------------------------------------------------------------------------
# Event status
# ---------------------------------------------------------------------------

@unique
class EventStatus(str, Enum):
    """Lifecycle status of a pipeline event."""

    PENDING = "pending"
    DISPATCHED = "dispatched"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"
    EXPIRED = "expired"

    @property
    def is_terminal(self) -> bool:
        return self in (EventStatus.COMPLETED, EventStatus.FAILED, EventStatus.EXPIRED)


# ---------------------------------------------------------------------------
# Priority levels
# ---------------------------------------------------------------------------

@unique
class EventPriority(int, Enum):
    """Event priority levels (lower number = higher priority)."""

    CRITICAL = 0
    HIGH = 1
    NORMAL = 2
    LOW = 3
    BACKGROUND = 4


# ---------------------------------------------------------------------------
# Pipeline Event
# ---------------------------------------------------------------------------

@dataclass(order=False)
class PipelineEvent:
    """
    A single event in the HCLS AI Factory pipeline.

    Attributes
    ----------
    event_type : EventType
        The semantic type of this event.
    source_stage : PipelineStage
        Which pipeline stage emitted this event.
    target_stage : PipelineStage, optional
        Intended recipient stage (None = broadcast).
    patient_id : str, optional
        Patient / sample identifier (e.g., "HG002").
    payload : dict
        Arbitrary event data (VCF path, gene list, candidate SMILES, etc.).
    priority : EventPriority
        Dispatch priority.
    status : EventStatus
        Current lifecycle status.
    event_id : str
        Unique event identifier (UUID).
    created_at : str
        ISO-8601 creation timestamp.
    dispatched_at : str, optional
        ISO-8601 dispatch timestamp.
    completed_at : str, optional
        ISO-8601 completion timestamp.
    error : str, optional
        Error message if status is FAILED.
    """

    event_type: EventType
    source_stage: PipelineStage
    target_stage: PipelineStage | None = None
    patient_id: str | None = None
    payload: dict[str, Any] = field(default_factory=dict)
    priority: EventPriority = EventPriority.NORMAL
    status: EventStatus = EventStatus.PENDING
    event_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    dispatched_at: str | None = None
    completed_at: str | None = None
    error: str | None = None

    # -- comparison for heapq (priority queue) --
    def __lt__(self, other: PipelineEvent) -> bool:
        if self.priority.value != other.priority.value:
            return self.priority.value < other.priority.value
        return self.created_at < other.created_at

    def __le__(self, other: PipelineEvent) -> bool:
        return self == other or self < other

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary (JSON-safe)."""
        d = {
            "event_id": self.event_id,
            "event_type": self.event_type.value,
            "source_stage": self.source_stage.value,
            "target_stage": self.target_stage.value if self.target_stage else None,
            "patient_id": self.patient_id,
            "payload": self.payload,
            "priority": self.priority.value,
            "status": self.status.value,
            "created_at": self.created_at,
            "dispatched_at": self.dispatched_at,
            "completed_at": self.completed_at,
            "error": self.error,
        }
        return d

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> PipelineEvent:
        """Deserialize from dictionary."""
        return cls(
            event_type=EventType(d["event_type"]),
            source_stage=PipelineStage(d["source_stage"]),
            target_stage=PipelineStage(d["target_stage"]) if d.get("target_stage") else None,
            patient_id=d.get("patient_id"),
            payload=d.get("payload", {}),
            priority=EventPriority(d.get("priority", EventPriority.NORMAL.value)),
            status=EventStatus(d.get("status", EventStatus.PENDING.value)),
            event_id=d.get("event_id", str(uuid.uuid4())),
            created_at=d.get("created_at", datetime.now().isoformat()),
            dispatched_at=d.get("dispatched_at"),
            completed_at=d.get("completed_at"),
            error=d.get("error"),
        )

    def mark_dispatched(self) -> None:
        self.status = EventStatus.DISPATCHED
        self.dispatched_at = datetime.now().isoformat()

    def mark_completed(self) -> None:
        self.status = EventStatus.COMPLETED
        self.completed_at = datetime.now().isoformat()

    def mark_failed(self, error: str) -> None:
        self.status = EventStatus.FAILED
        self.completed_at = datetime.now().isoformat()
        self.error = error


# ---------------------------------------------------------------------------
# Patient Context Bus (unified patient case)
# ---------------------------------------------------------------------------

@dataclass
class PatientCase:
    """Unified patient context that propagates across all agents.

    This is the central data carrier for the HCLS AI Factory. As each agent
    processes data, it enriches the PatientCase with its findings. This enables
    the closed-loop architecture where downstream agents benefit from upstream
    intelligence.
    """
    patient_id: str
    case_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())

    # Genomics (Stage 1)
    vcf_path: str | None = None
    variants: list[dict[str, Any]] = field(default_factory=list)
    key_genes: list[str] = field(default_factory=list)

    # Biomarker Agent
    biomarkers: dict[str, float] = field(default_factory=dict)
    biological_age: float | None = None
    age_acceleration: float | None = None
    disease_trajectories: list[dict[str, Any]] = field(default_factory=list)
    pgx_results: list[dict[str, Any]] = field(default_factory=list)
    critical_alerts: list[str] = field(default_factory=list)

    # Oncology Agent
    tumor_profile: dict[str, Any] = field(default_factory=dict)
    therapy_recommendations: list[dict[str, Any]] = field(default_factory=list)
    matched_trials: list[dict[str, Any]] = field(default_factory=list)

    # Imaging Agent
    imaging_findings: list[dict[str, Any]] = field(default_factory=list)
    concordance_scores: dict[str, float] = field(default_factory=dict)

    # CAR-T Agent
    cart_evaluation: dict[str, Any] = field(default_factory=dict)

    # Autoimmune Agent
    autoantibody_panel: dict[str, Any] = field(default_factory=dict)
    disease_activity_scores: dict[str, float] = field(default_factory=dict)
    flare_risk: float | None = None

    # Drug Discovery (Stage 3)
    drug_targets: list[str] = field(default_factory=list)
    drug_candidates: list[dict[str, Any]] = field(default_factory=list)
    pgx_filtered_candidates: list[dict[str, Any]] = field(default_factory=list)

    # Metadata
    ancestry: str | None = None
    age: int | None = None
    sex: str | None = None

    # Provenance
    processing_log: list[dict[str, Any]] = field(default_factory=list)

    def add_processing_step(self, agent: str, action: str, details: dict[str, Any] | None = None) -> None:
        """Record a processing step in the provenance log."""
        self.processing_log.append({
            "agent": agent,
            "action": action,
            "timestamp": datetime.now().isoformat(),
            "details": details or {},
        })

    def to_dict(self) -> dict[str, Any]:
        """Serialize to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "PatientCase":
        """Deserialize from dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ---------------------------------------------------------------------------
# Subscriber type
# ---------------------------------------------------------------------------

# A handler is a callable that takes a PipelineEvent and returns None.
EventHandler = Callable[[PipelineEvent], None]


@dataclass
class Subscription:
    """A registered event subscription."""

    handler: EventHandler
    event_types: set[EventType]
    name: str = ""
    source_filter: PipelineStage | None = None
    target_filter: PipelineStage | None = None
    priority_filter: EventPriority | None = None

    def matches(self, event: PipelineEvent) -> bool:
        """Check whether this subscription matches the given event."""
        if event.event_type not in self.event_types:
            return False
        if self.source_filter and event.source_stage != self.source_filter:
            return False
        if self.target_filter and event.target_stage != self.target_filter:
            return False
        if self.priority_filter and event.priority.value > self.priority_filter.value:
            return False
        return True


# ---------------------------------------------------------------------------
# Event Bus (singleton)
# ---------------------------------------------------------------------------

class EventBus:
    """
    Central event bus for the HCLS AI Factory pipeline.

    Implements publish/subscribe with priority queue dispatch, file-based
    audit persistence, and optional event replay for crash recovery.

    Thread-safe: all mutations are guarded by a reentrant lock.

    Usage::

        bus = EventBus.get_instance()

        # Subscribe
        bus.subscribe(
            handler=my_handler,
            event_types={EventType.VCF_READY},
            name="vcf_processor",
        )

        # Publish
        event = PipelineEvent(
            event_type=EventType.VCF_READY,
            source_stage=PipelineStage.GENOMICS,
            payload={"vcf_path": "/data/output/HG002.vcf.gz"},
        )
        bus.publish(event)

    Parameters
    ----------
    audit_dir : Path, optional
        Directory for JSON-lines audit logs.  Default: ``./data/events/``
    max_queue_size : int
        Maximum events in the priority queue before back-pressure.
    enable_audit : bool
        Whether to write events to the audit log file.
    """

    _instance: EventBus | None = None
    _instance_lock = threading.Lock()

    @classmethod
    def get_instance(
        cls,
        audit_dir: Path | None = None,
        max_queue_size: int = 10000,
        enable_audit: bool = True,
    ) -> EventBus:
        """Return the singleton EventBus instance (create on first call)."""
        with cls._instance_lock:
            if cls._instance is None:
                cls._instance = cls(
                    audit_dir=audit_dir,
                    max_queue_size=max_queue_size,
                    enable_audit=enable_audit,
                )
            return cls._instance

    @classmethod
    def reset_instance(cls) -> None:
        """Reset the singleton (for testing)."""
        with cls._instance_lock:
            if cls._instance is not None:
                cls._instance.shutdown()
            cls._instance = None

    def __init__(
        self,
        audit_dir: Path | None = None,
        max_queue_size: int = 10000,
        enable_audit: bool = True,
    ):
        self._lock = threading.RLock()
        self._subscriptions: list[Subscription] = []
        self._queue: list[PipelineEvent] = []  # heapq
        self._max_queue_size = max_queue_size
        self._enable_audit = enable_audit
        self._running = True

        # Event history (bounded)
        self._history: list[PipelineEvent] = []
        self._max_history = 5000

        # Statistics
        self._stats: dict[str, int] = defaultdict(int)

        # Audit log
        self._audit_dir = audit_dir or Path(
            os.environ.get("HCLS_EVENT_AUDIT_DIR", "data/events")
        )
        self._audit_file: Path | None = None
        self._audit_handle = None
        if self._enable_audit:
            self._init_audit_log()

        # Background dispatch thread
        self._dispatch_event = threading.Event()
        self._dispatch_thread = threading.Thread(
            target=self._dispatch_loop,
            name="hcls-event-dispatch",
            daemon=True,
        )
        self._dispatch_thread.start()

        logger.info(
            f"EventBus initialized: audit={'enabled' if enable_audit else 'disabled'}, "
            f"max_queue={max_queue_size}"
        )

    # ---- lifecycle ---------------------------------------------------------

    def shutdown(self) -> None:
        """Gracefully shut down the event bus."""
        self._running = False
        self._dispatch_event.set()  # wake the dispatch thread
        if self._dispatch_thread.is_alive():
            self._dispatch_thread.join(timeout=5.0)
        if self._audit_handle:
            self._audit_handle.close()
            self._audit_handle = None
        logger.info("EventBus shut down")

    # ---- subscribe / unsubscribe -------------------------------------------

    def subscribe(
        self,
        handler: EventHandler,
        event_types: set[EventType] | EventType | None = None,
        name: str = "",
        source_filter: PipelineStage | None = None,
        target_filter: PipelineStage | None = None,
        priority_filter: EventPriority | None = None,
    ) -> str:
        """
        Register an event handler.

        Parameters
        ----------
        handler : callable
            Function taking a PipelineEvent.
        event_types : set or single EventType
            Event types to subscribe to.  None = all types.
        name : str
            Human-readable name for logging.
        source_filter : PipelineStage, optional
            Only receive events from this source stage.
        target_filter : PipelineStage, optional
            Only receive events targeting this stage.
        priority_filter : EventPriority, optional
            Only receive events at this priority or higher.

        Returns
        -------
        str
            Subscription name (for unsubscribe).
        """
        if event_types is None:
            types_set = set(EventType)
        elif isinstance(event_types, EventType):
            types_set = {event_types}
        else:
            types_set = set(event_types)

        if not name:
            name = f"sub-{handler.__name__}-{len(self._subscriptions)}"

        sub = Subscription(
            handler=handler,
            event_types=types_set,
            name=name,
            source_filter=source_filter,
            target_filter=target_filter,
            priority_filter=priority_filter,
        )

        with self._lock:
            self._subscriptions.append(sub)

        logger.info(
            f"Subscribed '{name}' to events: "
            f"{[t.value for t in types_set]}"
        )
        return name

    def unsubscribe(self, name: str) -> bool:
        """
        Remove a subscription by name.

        Returns True if a subscription was removed.
        """
        with self._lock:
            before = len(self._subscriptions)
            self._subscriptions = [
                s for s in self._subscriptions if s.name != name
            ]
            removed = len(self._subscriptions) < before

        if removed:
            logger.info(f"Unsubscribed '{name}'")
        return removed

    # ---- publish -----------------------------------------------------------

    def publish(self, event: PipelineEvent) -> str:
        """
        Publish an event to the bus (synchronous dispatch).

        The event is dispatched immediately to all matching handlers in the
        calling thread, then written to the audit log.

        Returns the event_id.
        """
        with self._lock:
            self._stats["published"] += 1
            event.mark_dispatched()

            # Write to audit log
            self._audit_write(event)

            # Record Prometheus metric
            if EVENTS_EMITTED is not None:
                EVENTS_EMITTED.labels(
                    event_type=event.event_type.value,
                    source_stage=event.source_stage.value,
                ).inc()

            # Dispatch to handlers
            self._dispatch_event_to_handlers(event)

            # Add to history
            self._add_to_history(event)

        return event.event_id

    def publish_async(self, event: PipelineEvent) -> str:
        """
        Publish an event to the priority queue for background dispatch.

        The event is enqueued and dispatched by the background thread.
        Returns the event_id.
        """
        with self._lock:
            if len(self._queue) >= self._max_queue_size:
                logger.warning(
                    f"Event queue full ({self._max_queue_size}), "
                    f"dropping event {event.event_id}"
                )
                event.mark_failed("Queue full -- back-pressure")
                self._audit_write(event)
                return event.event_id

            heapq.heappush(self._queue, event)
            self._stats["enqueued"] += 1
            self._audit_write(event)

        # Signal the dispatch thread
        self._dispatch_event.set()
        return event.event_id

    # ---- query -------------------------------------------------------------

    def get_event(self, event_id: str) -> PipelineEvent | None:
        """Look up an event by ID from history."""
        with self._lock:
            for ev in reversed(self._history):
                if ev.event_id == event_id:
                    return ev
        return None

    def get_history(
        self,
        event_type: EventType | None = None,
        source_stage: PipelineStage | None = None,
        patient_id: str | None = None,
        limit: int = 100,
    ) -> list[PipelineEvent]:
        """
        Query event history with optional filters.
        """
        with self._lock:
            filtered = self._history[:]

        if event_type:
            filtered = [e for e in filtered if e.event_type == event_type]
        if source_stage:
            filtered = [e for e in filtered if e.source_stage == source_stage]
        if patient_id:
            filtered = [e for e in filtered if e.patient_id == patient_id]

        # Most recent first
        filtered.sort(key=lambda e: e.created_at, reverse=True)
        return filtered[:limit]

    def get_stats(self) -> dict[str, Any]:
        """Return event bus statistics."""
        with self._lock:
            return {
                "published": self._stats.get("published", 0),
                "enqueued": self._stats.get("enqueued", 0),
                "dispatched": self._stats.get("dispatched", 0),
                "handler_errors": self._stats.get("handler_errors", 0),
                "queue_size": len(self._queue),
                "history_size": len(self._history),
                "subscription_count": len(self._subscriptions),
                "audit_enabled": self._enable_audit,
                "running": self._running,
            }

    def get_queue_depth(self) -> int:
        """Current number of events in the priority queue."""
        with self._lock:
            return len(self._queue)

    # ---- replay ------------------------------------------------------------

    def replay_from_audit(
        self,
        audit_file: Path | None = None,
        since: str | None = None,
        event_types: set[EventType] | None = None,
    ) -> int:
        """
        Replay events from the audit log file.

        Parameters
        ----------
        audit_file : Path, optional
            Path to the audit log.  Default: current audit file.
        since : str, optional
            ISO-8601 timestamp -- only replay events after this time.
        event_types : set, optional
            Only replay these event types.

        Returns
        -------
        int
            Number of events replayed.
        """
        target = audit_file or self._audit_file
        if not target or not target.exists():
            logger.warning(f"Audit file not found for replay: {target}")
            return 0

        replayed = 0
        with open(target, "r") as f:
            for line_no, line in enumerate(f, start=1):
                line = line.strip()
                if not line:
                    continue
                try:
                    data = json.loads(line)
                    event = PipelineEvent.from_dict(data)

                    # Apply filters
                    if since and event.created_at < since:
                        continue
                    if event_types and event.event_type not in event_types:
                        continue

                    # Re-dispatch
                    self._dispatch_event_to_handlers(event)
                    replayed += 1
                except json.JSONDecodeError as exc:
                    logger.warning(f"Malformed JSON at line {line_no}: {exc}")
                except Exception as exc:
                    logger.warning(f"Error replaying event at line {line_no}: {exc}")

        logger.info(f"Replayed {replayed} events from {target}")
        return replayed

    # ---- audit log ---------------------------------------------------------

    def _init_audit_log(self) -> None:
        """Initialize the audit log directory and file."""
        try:
            self._audit_dir.mkdir(parents=True, exist_ok=True)
            date_str = datetime.now().strftime("%Y%m%d")
            self._audit_file = self._audit_dir / f"events_{date_str}.jsonl"
            self._audit_handle = open(self._audit_file, "a")
            logger.info(f"Audit log: {self._audit_file}")
        except Exception as exc:
            logger.error(f"Failed to initialize audit log: {exc}")
            self._enable_audit = False

    def _audit_write(self, event: PipelineEvent) -> None:
        """Append an event to the audit log (JSON lines format)."""
        if not self._enable_audit or not self._audit_handle:
            return
        try:
            line = json.dumps(event.to_dict(), default=str)
            self._audit_handle.write(line + "\n")
            self._audit_handle.flush()
        except Exception as exc:
            logger.error(f"Audit write failed: {exc}")

    def get_audit_file(self) -> Path | None:
        """Return the current audit log file path."""
        return self._audit_file

    # ---- dispatch ----------------------------------------------------------

    def _dispatch_event_to_handlers(self, event: PipelineEvent) -> None:
        """Dispatch a single event to all matching subscribers."""
        with self._lock:
            matching_subs = [
                s for s in self._subscriptions if s.matches(event)
            ]

        for sub in matching_subs:
            t0 = time.time()
            try:
                sub.handler(event)
                elapsed = time.time() - t0

                with self._lock:
                    self._stats["dispatched"] += 1

                if EVENTS_PROCESSED is not None:
                    EVENTS_PROCESSED.labels(
                        event_type=event.event_type.value,
                        handler=sub.name,
                    ).inc()
                if EVENT_LATENCY is not None:
                    EVENT_LATENCY.labels(event_type=event.event_type.value).observe(elapsed)

                logger.debug(
                    f"Handler '{sub.name}' processed {event.event_type.value} "
                    f"in {elapsed:.4f}s"
                )
            except Exception as exc:
                elapsed = time.time() - t0
                with self._lock:
                    self._stats["handler_errors"] += 1

                logger.error(
                    f"Handler '{sub.name}' failed for {event.event_type.value}: {exc}"
                )
                # Do not re-raise -- other handlers should still run

    def _dispatch_loop(self) -> None:
        """Background thread that drains the priority queue."""
        logger.debug("Event dispatch thread started")

        while self._running:
            # Wait for work or shutdown signal
            self._dispatch_event.wait(timeout=1.0)
            self._dispatch_event.clear()

            # Drain the queue
            while self._running:
                event: PipelineEvent | None = None
                with self._lock:
                    if self._queue:
                        event = heapq.heappop(self._queue)

                if event is None:
                    break

                event.mark_dispatched()
                self._dispatch_event_to_handlers(event)
                self._audit_write(event)
                self._add_to_history(event)

        logger.debug("Event dispatch thread stopped")

    def _add_to_history(self, event: PipelineEvent) -> None:
        """Add event to bounded history."""
        with self._lock:
            self._history.append(event)
            if len(self._history) > self._max_history:
                self._history = self._history[-self._max_history:]


# ---------------------------------------------------------------------------
# Module-level convenience functions
# ---------------------------------------------------------------------------

def get_event_bus(**kwargs: Any) -> EventBus:
    """
    Get the singleton EventBus instance.

    Keyword arguments are forwarded to ``EventBus.get_instance()`` on first call.
    """
    return EventBus.get_instance(**kwargs)


def publish_event(
    event_type: EventType,
    source_stage: PipelineStage,
    target_stage: PipelineStage | None = None,
    patient_id: str | None = None,
    payload: dict[str, Any] | None = None,
    priority: EventPriority = EventPriority.NORMAL,
) -> str:
    """
    Convenience function to create and publish an event in one call.

    Returns the event_id.
    """
    event = PipelineEvent(
        event_type=event_type,
        source_stage=source_stage,
        target_stage=target_stage,
        patient_id=patient_id,
        payload=payload or {},
        priority=priority,
    )
    bus = get_event_bus()
    return bus.publish(event)


def create_pipeline_event(
    event_type: EventType,
    source_stage: PipelineStage,
    target_stage: PipelineStage | None = None,
    patient_id: str | None = None,
    payload: dict[str, Any] | None = None,
    priority: EventPriority = EventPriority.NORMAL,
) -> PipelineEvent:
    """
    Create a PipelineEvent without publishing it.
    """
    return PipelineEvent(
        event_type=event_type,
        source_stage=source_stage,
        target_stage=target_stage,
        patient_id=patient_id,
        payload=payload or {},
        priority=priority,
    )
