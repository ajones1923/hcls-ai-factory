"""Prometheus metrics for the CAR-T Intelligence Agent.

Exposes counters, histograms, and gauges for query latency, collection hits,
LLM token usage, evidence counts, and ingest freshness.  Scraped by the
Grafana + Prometheus stack on port 9099.

All metrics use the ``cart_`` prefix so they are easily filterable in
Grafana dashboards alongside the existing HCLS AI Factory exporters
(node_exporter:9100, DCGM:9400).

If ``prometheus_client`` is not installed the module silently exports
no-op stubs so the rest of the application can import metrics helpers
without a hard dependency.

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

from typing import Dict

try:
    from prometheus_client import Counter, Gauge, Histogram, generate_latest

    # ── Histograms ────────────────────────────────────────────────────
    QUERY_LATENCY = Histogram(
        "cart_query_latency_seconds",
        "Query processing time",
        ["query_type"],
        buckets=[0.1, 0.5, 1, 2, 5, 10, 30],
    )

    EVIDENCE_COUNT = Histogram(
        "cart_evidence_count",
        "Evidence items per query",
        buckets=[0, 5, 10, 15, 20, 25, 30],
    )

    # ── Counters ──────────────────────────────────────────────────────
    QUERY_COUNT = Counter(
        "cart_queries_total",
        "Total queries processed",
        ["query_type", "status"],
    )

    COLLECTION_HITS = Counter(
        "cart_collection_hits_total",
        "Hits by collection",
        ["collection"],
    )

    LLM_TOKENS = Counter(
        "cart_llm_tokens_total",
        "LLM tokens used",
        ["direction"],
    )

    # ── Cross-collection & pipeline histograms ─────────────────────────
    CROSS_COLLECTION_QUERY_LATENCY = Histogram(
        "cart_cross_collection_query_latency_seconds",
        "Cross-collection query processing time",
        ["query_type"],
        buckets=[0.1, 0.5, 1, 2, 5, 10, 30],
    )

    CROSS_COLLECTION_RESULTS = Histogram(
        "cart_cross_collection_results_count",
        "Number of results returned from cross-collection queries",
        buckets=[0, 1, 5, 10, 20, 50, 100],
    )

    LLM_API_LATENCY = Histogram(
        "cart_llm_api_latency_seconds",
        "LLM API call latency",
        ["provider", "model"],
        buckets=[0.5, 1, 2, 5, 10, 30, 60],
    )

    EMBEDDING_LATENCY = Histogram(
        "cart_embedding_latency_seconds",
        "Embedding generation latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1],
    )

    PIPELINE_STAGE_DURATION = Histogram(
        "cart_pipeline_stage_duration_seconds",
        "Duration of individual pipeline stages",
        ["stage"],
        buckets=[0.1, 0.5, 1, 5, 10, 30, 60, 120],
    )

    MILVUS_SEARCH_LATENCY = Histogram(
        "cart_milvus_search_latency_seconds",
        "Milvus vector search latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1, 2],
    )

    MILVUS_UPSERT_LATENCY = Histogram(
        "cart_milvus_upsert_latency_seconds",
        "Milvus upsert latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1, 5],
    )

    # ── Additional counters ──────────────────────────────────────────
    LLM_COST_ESTIMATE = Counter(
        "cart_llm_cost_estimate_usd",
        "Estimated LLM cost in USD",
        ["model"],
    )

    EMBEDDING_CACHE_HITS = Counter(
        "cart_embedding_cache_hits_total",
        "Embedding cache hits",
    )

    EMBEDDING_CACHE_MISSES = Counter(
        "cart_embedding_cache_misses_total",
        "Embedding cache misses",
    )

    CIRCUIT_BREAKER_TRIPS = Counter(
        "cart_circuit_breaker_trips_total",
        "Circuit breaker trip events",
        ["service"],
    )

    EVENT_BUS_EVENTS_EMITTED = Counter(
        "cart_event_bus_events_emitted_total",
        "Events emitted to the event bus",
        ["event_type"],
    )

    REPORT_GENERATED = Counter(
        "cart_reports_generated_total",
        "Reports generated",
        ["format"],
    )

    # ── Gauges ────────────────────────────────────────────────────────
    ACTIVE_CONNECTIONS = Gauge(
        "cart_active_connections",
        "Active connections",
    )

    COLLECTION_SIZE = Gauge(
        "cart_collection_size",
        "Records per collection",
        ["collection"],
    )

    LAST_INGEST = Gauge(
        "cart_last_ingest_timestamp",
        "Last ingest timestamp",
        ["source"],
    )

    CIRCUIT_BREAKER_STATE = Gauge(
        "cart_circuit_breaker_state",
        "Circuit breaker state (0=closed, 1=open, 2=half-open)",
        ["service"],
    )

    _PROMETHEUS_AVAILABLE = True

except ImportError:
    # ── No-op stubs when prometheus_client is not installed ────────────
    _PROMETHEUS_AVAILABLE = False

    class _NoOpLabeled:
        """Stub that silently ignores .labels().observe/inc/set calls."""

        def labels(self, *args, **kwargs):
            return self

        def observe(self, *args, **kwargs):
            pass

        def inc(self, *args, **kwargs):
            pass

        def dec(self, *args, **kwargs):
            pass

        def set(self, *args, **kwargs):
            pass

    class _NoOpGauge:
        """Stub for label-less Gauge (ACTIVE_CONNECTIONS)."""

        def inc(self, *args, **kwargs):
            pass

        def dec(self, *args, **kwargs):
            pass

        def set(self, *args, **kwargs):
            pass

    QUERY_LATENCY = _NoOpLabeled()        # type: ignore[assignment]
    EVIDENCE_COUNT = _NoOpLabeled()        # type: ignore[assignment]
    QUERY_COUNT = _NoOpLabeled()           # type: ignore[assignment]
    COLLECTION_HITS = _NoOpLabeled()       # type: ignore[assignment]
    LLM_TOKENS = _NoOpLabeled()            # type: ignore[assignment]
    ACTIVE_CONNECTIONS = _NoOpGauge()      # type: ignore[assignment]
    COLLECTION_SIZE = _NoOpLabeled()       # type: ignore[assignment]
    LAST_INGEST = _NoOpLabeled()           # type: ignore[assignment]

    # New metrics — no-op stubs
    CROSS_COLLECTION_QUERY_LATENCY = _NoOpLabeled()   # type: ignore[assignment]
    CROSS_COLLECTION_RESULTS = _NoOpLabeled()          # type: ignore[assignment]
    LLM_API_LATENCY = _NoOpLabeled()                   # type: ignore[assignment]
    EMBEDDING_LATENCY = _NoOpLabeled()                 # type: ignore[assignment]
    PIPELINE_STAGE_DURATION = _NoOpLabeled()           # type: ignore[assignment]
    MILVUS_SEARCH_LATENCY = _NoOpLabeled()             # type: ignore[assignment]
    MILVUS_UPSERT_LATENCY = _NoOpLabeled()             # type: ignore[assignment]
    LLM_COST_ESTIMATE = _NoOpLabeled()                 # type: ignore[assignment]
    EMBEDDING_CACHE_HITS = _NoOpGauge()                # type: ignore[assignment]
    EMBEDDING_CACHE_MISSES = _NoOpGauge()              # type: ignore[assignment]
    CIRCUIT_BREAKER_TRIPS = _NoOpLabeled()             # type: ignore[assignment]
    EVENT_BUS_EVENTS_EMITTED = _NoOpLabeled()          # type: ignore[assignment]
    REPORT_GENERATED = _NoOpLabeled()                  # type: ignore[assignment]
    CIRCUIT_BREAKER_STATE = _NoOpLabeled()             # type: ignore[assignment]

    def generate_latest() -> bytes:  # type: ignore[misc]
        return b""


# ═════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═════════════════════════════════════════════════════════════════════════


def record_query(
    query_type: str,
    latency: float,
    hit_count: int,
    status: str = "success",
) -> None:
    """Record metrics for a completed query.

    Call this once at the end of every RAG query or agent run to update
    latency, hit count, and total-query counters in a single place.

    Args:
        query_type: Query category — e.g. ``"rag"``, ``"agent"``,
            ``"comparative"``, ``"entity_link"``.
        latency: Wall-clock processing time in **seconds**.
        hit_count: Total evidence items returned to the user.
        status: Outcome label — ``"success"`` or ``"error"``.
    """
    QUERY_LATENCY.labels(query_type=query_type).observe(latency)
    QUERY_COUNT.labels(query_type=query_type, status=status).inc()
    EVIDENCE_COUNT.observe(hit_count)


def record_collection_hits(hits_by_collection: Dict[str, int]) -> None:
    """Increment per-collection hit counters.

    Args:
        hits_by_collection: Mapping of collection name (e.g.
            ``"cart_literature"``) to the number of hits returned from
            that collection in a single query.
    """
    for collection, count in hits_by_collection.items():
        COLLECTION_HITS.labels(collection=collection).inc(count)


def update_collection_sizes(stats: Dict[str, int]) -> None:
    """Set the current record count for each collection.

    Typically called after ``CARTCollectionManager.get_collection_stats()``
    returns its ``{collection_name: row_count}`` dict.

    Args:
        stats: Mapping of collection name to entity/row count.
    """
    for collection, size in stats.items():
        COLLECTION_SIZE.labels(collection=collection).set(size)


def record_cross_collection_query(
    query_type: str,
    latency: float,
    result_count: int,
) -> None:
    """Record metrics for a cross-collection query.

    Args:
        query_type: Category (e.g. ``"multi_rag"``, ``"entity_link"``).
        latency: Wall-clock time in **seconds**.
        result_count: Number of results returned.
    """
    CROSS_COLLECTION_QUERY_LATENCY.labels(query_type=query_type).observe(latency)
    CROSS_COLLECTION_RESULTS.observe(result_count)


def record_llm_call(
    provider: str,
    model: str,
    latency: float,
    estimated_cost: float = 0.0,
) -> None:
    """Record metrics for an LLM API call.

    Args:
        provider: LLM provider (e.g. ``"anthropic"``, ``"openai"``).
        model: Model identifier (e.g. ``"claude-sonnet-4-20250514"``).
        latency: Call latency in **seconds**.
        estimated_cost: Estimated cost in USD (default 0).
    """
    LLM_API_LATENCY.labels(provider=provider, model=model).observe(latency)
    if estimated_cost > 0:
        LLM_COST_ESTIMATE.labels(model=model).inc(estimated_cost)


def record_embedding(latency: float, cache_hit: bool = False) -> None:
    """Record metrics for an embedding operation.

    Args:
        latency: Embedding generation time in **seconds**.
        cache_hit: Whether the embedding was served from cache.
    """
    EMBEDDING_LATENCY.observe(latency)
    if cache_hit:
        EMBEDDING_CACHE_HITS.inc()
    else:
        EMBEDDING_CACHE_MISSES.inc()


def record_circuit_breaker(
    service: str,
    state: int,
    tripped: bool = False,
) -> None:
    """Record circuit breaker state change.

    Args:
        service: Service name (e.g. ``"milvus"``, ``"llm"``).
        state: Numeric state (0=closed, 1=open, 2=half-open).
        tripped: Whether the breaker just tripped (increments trip counter).
    """
    CIRCUIT_BREAKER_STATE.labels(service=service).set(state)
    if tripped:
        CIRCUIT_BREAKER_TRIPS.labels(service=service).inc()


def record_pipeline_stage(stage: str, duration: float) -> None:
    """Record duration for a pipeline stage.

    Args:
        stage: Stage name (e.g. ``"embed"``, ``"retrieve"``, ``"synthesize"``).
        duration: Stage duration in **seconds**.
    """
    PIPELINE_STAGE_DURATION.labels(stage=stage).observe(duration)


def record_milvus_search(latency: float) -> None:
    """Record Milvus vector search latency.

    Args:
        latency: Search time in **seconds**.
    """
    MILVUS_SEARCH_LATENCY.observe(latency)


def record_milvus_upsert(latency: float) -> None:
    """Record Milvus upsert latency.

    Args:
        latency: Upsert time in **seconds**.
    """
    MILVUS_UPSERT_LATENCY.observe(latency)


def record_event_emitted(event_type: str) -> None:
    """Record an event emitted to the event bus.

    Args:
        event_type: Event category (e.g. ``"query_complete"``, ``"ingest_done"``).
    """
    EVENT_BUS_EVENTS_EMITTED.labels(event_type=event_type).inc()


def record_report_generated(fmt: str) -> None:
    """Record a report generation.

    Args:
        fmt: Report format (e.g. ``"pdf"``, ``"markdown"``, ``"json"``).
    """
    REPORT_GENERATED.labels(format=fmt).inc()


def get_metrics_text() -> str:
    """Return the current Prometheus metrics exposition in text format.

    This is the string you serve at ``/metrics`` (typically via FastAPI
    or a dedicated ``prometheus_client.start_http_server``).

    Returns:
        UTF-8 decoded Prometheus exposition text, or an empty string if
        ``prometheus_client`` is not installed.
    """
    return generate_latest().decode("utf-8")
