"""
Prometheus metrics for the Oncology Intelligence Agent.

Exposes histograms, counters, and gauges for monitoring query performance,
evidence retrieval, LLM usage, embedding latency, Milvus operations,
circuit-breaker state, and pipeline stages. Falls back to lightweight
no-op stubs when ``prometheus_client`` is not installed.

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)

try:
    from prometheus_client import (
        Counter,
        Gauge,
        Histogram,
        generate_latest,
    )

    _PROMETHEUS_AVAILABLE = True
except ImportError:
    _PROMETHEUS_AVAILABLE = False

# -----------------------------------------------------------------------
# No-op stubs for environments without prometheus_client
# -----------------------------------------------------------------------

if not _PROMETHEUS_AVAILABLE:
    logger.info(
        "prometheus_client not installed – metrics will be no-ops"
    )

    class _NoOpMetric:
        """Drop-in stub that silently ignores all metric operations."""

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

    # Histograms
    QUERY_LATENCY = _NoOpMetric()
    EVIDENCE_COUNT = _NoOpMetric()
    CROSS_COLLECTION_QUERY_LATENCY = _NoOpMetric()
    CROSS_COLLECTION_RESULTS = _NoOpMetric()
    LLM_API_LATENCY = _NoOpMetric()
    EMBEDDING_LATENCY = _NoOpMetric()
    PIPELINE_STAGE_DURATION = _NoOpMetric()
    MILVUS_SEARCH_LATENCY = _NoOpMetric()
    MILVUS_UPSERT_LATENCY = _NoOpMetric()

    # Counters
    QUERY_COUNT = _NoOpMetric()
    COLLECTION_HITS = _NoOpMetric()
    LLM_TOKENS = _NoOpMetric()
    LLM_COST_ESTIMATE = _NoOpMetric()
    EMBEDDING_CACHE_HITS = _NoOpMetric()
    EMBEDDING_CACHE_MISSES = _NoOpMetric()
    CIRCUIT_BREAKER_TRIPS = _NoOpMetric()
    EVENT_BUS_EVENTS_EMITTED = _NoOpMetric()
    REPORT_GENERATED = _NoOpMetric()
    MTB_PACKETS_GENERATED = _NoOpMetric()
    TRIAL_MATCHES_PERFORMED = _NoOpMetric()

    # Gauges
    ACTIVE_CONNECTIONS = _NoOpMetric()
    COLLECTION_SIZE = _NoOpMetric()
    LAST_INGEST = _NoOpMetric()
    CIRCUIT_BREAKER_STATE = _NoOpMetric()

else:
    # -------------------------------------------------------------------
    # Histograms
    # -------------------------------------------------------------------

    QUERY_LATENCY = Histogram(
        "onco_query_latency_seconds",
        "End-to-end latency for a single agent query",
        buckets=(0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0),
    )

    EVIDENCE_COUNT = Histogram(
        "onco_evidence_count",
        "Number of evidence items returned per query",
        buckets=(0, 1, 3, 5, 10, 25, 50, 100),
    )

    CROSS_COLLECTION_QUERY_LATENCY = Histogram(
        "onco_cross_collection_query_latency_seconds",
        "Latency for cross-collection search operations",
        labelnames=["strategy"],
        buckets=(0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0),
    )

    CROSS_COLLECTION_RESULTS = Histogram(
        "onco_cross_collection_results",
        "Number of results returned from cross-collection queries",
        labelnames=["collection"],
        buckets=(0, 1, 5, 10, 25, 50, 100),
    )

    LLM_API_LATENCY = Histogram(
        "onco_llm_api_latency_seconds",
        "Latency for LLM API calls (synthesis, planning)",
        labelnames=["operation"],
        buckets=(0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0),
    )

    EMBEDDING_LATENCY = Histogram(
        "onco_embedding_latency_seconds",
        "Latency for embedding generation",
        buckets=(0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
    )

    PIPELINE_STAGE_DURATION = Histogram(
        "onco_pipeline_stage_duration_seconds",
        "Duration of each pipeline stage (plan, search, evaluate, synthesize)",
        labelnames=["stage"],
        buckets=(0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0),
    )

    MILVUS_SEARCH_LATENCY = Histogram(
        "onco_milvus_search_latency_seconds",
        "Latency for Milvus vector search operations",
        labelnames=["collection"],
        buckets=(0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.5),
    )

    MILVUS_UPSERT_LATENCY = Histogram(
        "onco_milvus_upsert_latency_seconds",
        "Latency for Milvus upsert operations",
        labelnames=["collection"],
        buckets=(0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0),
    )

    # -------------------------------------------------------------------
    # Counters
    # -------------------------------------------------------------------

    QUERY_COUNT = Counter(
        "onco_query_count_total",
        "Total number of queries processed by the agent",
        labelnames=["strategy"],
    )

    COLLECTION_HITS = Counter(
        "onco_collection_hits_total",
        "Total hits per Milvus collection",
        labelnames=["collection"],
    )

    LLM_TOKENS = Counter(
        "onco_llm_tokens_total",
        "Total LLM tokens consumed",
        labelnames=["direction"],  # "input" or "output"
    )

    LLM_COST_ESTIMATE = Counter(
        "onco_llm_cost_estimate_dollars",
        "Estimated cumulative LLM cost in USD",
    )

    EMBEDDING_CACHE_HITS = Counter(
        "onco_embedding_cache_hits_total",
        "Number of embedding cache hits",
    )

    EMBEDDING_CACHE_MISSES = Counter(
        "onco_embedding_cache_misses_total",
        "Number of embedding cache misses",
    )

    CIRCUIT_BREAKER_TRIPS = Counter(
        "onco_circuit_breaker_trips_total",
        "Number of circuit-breaker trip events",
        labelnames=["service"],
    )

    EVENT_BUS_EVENTS_EMITTED = Counter(
        "onco_event_bus_events_emitted_total",
        "Total events emitted on the internal event bus",
        labelnames=["event_type"],
    )

    REPORT_GENERATED = Counter(
        "onco_report_generated_total",
        "Total oncology intelligence reports generated",
    )

    MTB_PACKETS_GENERATED = Counter(
        "onco_mtb_packets_generated_total",
        "Total molecular tumor board packets generated",
    )

    TRIAL_MATCHES_PERFORMED = Counter(
        "onco_trial_matches_performed_total",
        "Total clinical trial matching operations performed",
    )

    # -------------------------------------------------------------------
    # Gauges
    # -------------------------------------------------------------------

    ACTIVE_CONNECTIONS = Gauge(
        "onco_active_connections",
        "Number of currently active client connections",
    )

    COLLECTION_SIZE = Gauge(
        "onco_collection_size",
        "Current number of entities in each Milvus collection",
        labelnames=["collection"],
    )

    LAST_INGEST = Gauge(
        "onco_last_ingest_timestamp",
        "Unix timestamp of the most recent ingest run",
        labelnames=["source"],
    )

    CIRCUIT_BREAKER_STATE = Gauge(
        "onco_circuit_breaker_state",
        "Circuit breaker state (0=closed, 1=half-open, 2=open)",
        labelnames=["service"],
    )


# -----------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------

def record_query(strategy: str, latency: float, evidence_count: int) -> None:
    """Record metrics for a completed agent query."""
    QUERY_COUNT.labels(strategy=strategy).inc()
    QUERY_LATENCY.observe(latency)
    EVIDENCE_COUNT.observe(evidence_count)


def record_collection_hits(collection: str, hits: int) -> None:
    """Record hit counts for a single collection search."""
    COLLECTION_HITS.labels(collection=collection).inc(hits)


def update_collection_sizes(sizes: dict) -> None:
    """Bulk-update the collection size gauge.

    Parameters
    ----------
    sizes:
        Mapping of collection name to entity count.
    """
    for collection, size in sizes.items():
        COLLECTION_SIZE.labels(collection=collection).set(size)


def record_cross_collection_query(
    strategy: str,
    latency: float,
    results_by_collection: Optional[dict] = None,
) -> None:
    """Record cross-collection query latency and per-collection result counts."""
    CROSS_COLLECTION_QUERY_LATENCY.labels(strategy=strategy).observe(latency)
    if results_by_collection:
        for collection, count in results_by_collection.items():
            CROSS_COLLECTION_RESULTS.labels(collection=collection).observe(count)


def record_llm_call(
    operation: str,
    latency: float,
    input_tokens: int = 0,
    output_tokens: int = 0,
    cost: float = 0.0,
) -> None:
    """Record an LLM API call."""
    LLM_API_LATENCY.labels(operation=operation).observe(latency)
    if input_tokens:
        LLM_TOKENS.labels(direction="input").inc(input_tokens)
    if output_tokens:
        LLM_TOKENS.labels(direction="output").inc(output_tokens)
    if cost > 0:
        LLM_COST_ESTIMATE.inc(cost)


def record_embedding(latency: float, cache_hit: bool = False) -> None:
    """Record an embedding operation."""
    EMBEDDING_LATENCY.observe(latency)
    if cache_hit:
        EMBEDDING_CACHE_HITS.inc()
    else:
        EMBEDDING_CACHE_MISSES.inc()


def record_circuit_breaker(service: str, state: int, tripped: bool = False) -> None:
    """Record circuit-breaker state change.

    Parameters
    ----------
    service:
        Name of the downstream service.
    state:
        0=closed, 1=half-open, 2=open.
    tripped:
        Whether this is a new trip event.
    """
    CIRCUIT_BREAKER_STATE.labels(service=service).set(state)
    if tripped:
        CIRCUIT_BREAKER_TRIPS.labels(service=service).inc()


def record_pipeline_stage(stage: str, duration: float) -> None:
    """Record the duration of a pipeline stage."""
    PIPELINE_STAGE_DURATION.labels(stage=stage).observe(duration)


def record_milvus_search(collection: str, latency: float) -> None:
    """Record a Milvus search operation."""
    MILVUS_SEARCH_LATENCY.labels(collection=collection).observe(latency)


def record_milvus_upsert(collection: str, latency: float) -> None:
    """Record a Milvus upsert operation."""
    MILVUS_UPSERT_LATENCY.labels(collection=collection).observe(latency)


def record_event_emitted(event_type: str) -> None:
    """Record an event emitted on the internal event bus."""
    EVENT_BUS_EVENTS_EMITTED.labels(event_type=event_type).inc()


def record_report_generated() -> None:
    """Increment the report-generated counter."""
    REPORT_GENERATED.inc()


def get_metrics_text() -> str:
    """Return the current Prometheus metrics as text.

    Returns an empty string if prometheus_client is not available.
    """
    if not _PROMETHEUS_AVAILABLE:
        return ""
    return generate_latest().decode("utf-8")
