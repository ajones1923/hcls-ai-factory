"""Prometheus metrics for the Clinical Trial Intelligence Agent.

Exposes counters, histograms, gauges, and info metrics for query latency,
collection hits, LLM token usage, workflow executions, patient matching,
safety signal detection, ingest operations, and system health.

Scraped by the Grafana + Prometheus stack alongside existing HCLS AI Factory
exporters (node_exporter:9100, DCGM:9400).

All metrics use the ``trial_`` prefix so they are easily filterable in
Grafana dashboards.

If ``prometheus_client`` is not installed the module silently exports
no-op stubs so the rest of the application can import metrics helpers
without a hard dependency.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import time
from typing import Any, Dict

try:
    from prometheus_client import Counter, Gauge, Histogram, Info, generate_latest

    # ── Query Metrics ──────────────────────────────────────────────────
    QUERY_TOTAL = Counter(
        "trial_queries_total",
        "Total queries processed",
        ["workflow_type"],
    )

    QUERY_LATENCY = Histogram(
        "trial_query_duration_seconds",
        "Query processing time",
        ["workflow_type"],
        buckets=[0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0],
    )

    QUERY_ERRORS = Counter(
        "trial_query_errors_total",
        "Total query errors",
        ["error_type"],
    )

    # ── RAG / Vector Search Metrics ────────────────────────────────────
    SEARCH_TOTAL = Counter(
        "trial_search_total",
        "Total vector searches",
        ["collection"],
    )

    SEARCH_LATENCY = Histogram(
        "trial_search_duration_seconds",
        "Vector search latency",
        ["collection"],
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0],
    )

    SEARCH_RESULTS = Histogram(
        "trial_search_results_count",
        "Number of results per search",
        ["collection"],
        buckets=[0, 1, 5, 10, 20, 50, 100],
    )

    EMBEDDING_LATENCY = Histogram(
        "trial_embedding_duration_seconds",
        "Embedding generation time",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0],
    )

    # ── LLM Metrics ────────────────────────────────────────────────────
    LLM_CALLS = Counter(
        "trial_llm_calls_total",
        "Total LLM calls",
        ["model"],
    )

    LLM_LATENCY = Histogram(
        "trial_llm_duration_seconds",
        "LLM call latency",
        ["model"],
        buckets=[0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 60.0],
    )

    LLM_TOKENS = Counter(
        "trial_llm_tokens_total",
        "Total LLM tokens",
        ["direction"],  # input / output
    )

    # ── Clinical Workflow Metrics ──────────────────────────────────────
    WORKFLOW_TOTAL = Counter(
        "trial_workflow_executions_total",
        "Clinical workflow executions",
        ["workflow_type"],
    )

    WORKFLOW_LATENCY = Histogram(
        "trial_workflow_duration_seconds",
        "Clinical workflow execution time",
        ["workflow_type"],
        buckets=[0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0],
    )

    # ── Patient Matching Metrics ───────────────────────────────────────
    MATCHING_TOTAL = Counter(
        "trial_matching_total",
        "Total patient-trial matching operations",
        ["therapeutic_area"],
    )

    MATCHING_SCORE = Histogram(
        "trial_matching_score",
        "Patient-trial match scores",
        ["phase"],
        buckets=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    )

    # ── Safety Signal Metrics ──────────────────────────────────────────
    SAFETY_SIGNALS = Counter(
        "trial_safety_signals_total",
        "Safety signals detected",
        ["severity"],
    )

    # ── Export Metrics ─────────────────────────────────────────────────
    EXPORT_TOTAL = Counter(
        "trial_exports_total",
        "Report exports",
        ["format"],
    )

    # ── System Metrics ─────────────────────────────────────────────────
    MILVUS_CONNECTED = Gauge(
        "trial_milvus_connected",
        "Milvus connection status (1=connected, 0=disconnected)",
    )

    COLLECTIONS_LOADED = Gauge(
        "trial_collections_loaded",
        "Number of loaded collections",
    )

    COLLECTION_SIZE = Gauge(
        "trial_collection_size",
        "Records per collection",
        ["collection"],
    )

    ACTIVE_CONNECTIONS = Gauge(
        "trial_active_connections",
        "Active client connections",
    )

    AGENT_INFO = Info(
        "trial_agent",
        "Agent version and configuration info",
    )

    # ── Ingest Metrics ─────────────────────────────────────────────────
    INGEST_TOTAL = Counter(
        "trial_ingest_total",
        "Total ingest operations",
        ["source"],
    )

    INGEST_RECORDS = Counter(
        "trial_ingest_records_total",
        "Total records ingested",
        ["collection"],
    )

    INGEST_ERRORS = Counter(
        "trial_ingest_errors_total",
        "Total ingest errors",
        ["source"],
    )

    INGEST_LATENCY = Histogram(
        "trial_ingest_duration_seconds",
        "Ingest operation time",
        ["source"],
        buckets=[1.0, 5.0, 10.0, 30.0, 60.0, 300.0, 600.0],
    )

    LAST_INGEST = Gauge(
        "trial_last_ingest_timestamp",
        "Last ingest timestamp (unix epoch)",
        ["source"],
    )

    # ── Pipeline Stage Metrics ─────────────────────────────────────────
    PIPELINE_STAGE_DURATION = Histogram(
        "trial_pipeline_stage_duration_seconds",
        "Duration of individual pipeline stages",
        ["stage"],
        buckets=[0.1, 0.5, 1.0, 5.0, 10.0, 30.0, 60.0, 120.0],
    )

    MILVUS_SEARCH_LATENCY = Histogram(
        "trial_milvus_search_latency_seconds",
        "Milvus vector search latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0],
    )

    MILVUS_UPSERT_LATENCY = Histogram(
        "trial_milvus_upsert_latency_seconds",
        "Milvus upsert latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 5.0],
    )

    _PROMETHEUS_AVAILABLE = True

except ImportError:
    # ── No-op stubs when prometheus_client is not installed ─────────────
    _PROMETHEUS_AVAILABLE = False

    class _NoOpLabeled:
        """Stub that silently ignores .labels().observe/inc/set calls."""

        def labels(self, *args: Any, **kwargs: Any) -> "_NoOpLabeled":
            return self

        def observe(self, *args: Any, **kwargs: Any) -> None:
            pass

        def inc(self, *args: Any, **kwargs: Any) -> None:
            pass

        def dec(self, *args: Any, **kwargs: Any) -> None:
            pass

        def set(self, *args: Any, **kwargs: Any) -> None:
            pass

    class _NoOpGauge:
        """Stub for label-less Gauge."""

        def inc(self, *args: Any, **kwargs: Any) -> None:
            pass

        def dec(self, *args: Any, **kwargs: Any) -> None:
            pass

        def set(self, *args: Any, **kwargs: Any) -> None:
            pass

    class _NoOpInfo:
        """Stub for Info metric."""

        def info(self, *args: Any, **kwargs: Any) -> None:
            pass

    QUERY_TOTAL = _NoOpLabeled()              # type: ignore[assignment]
    QUERY_LATENCY = _NoOpLabeled()            # type: ignore[assignment]
    QUERY_ERRORS = _NoOpLabeled()             # type: ignore[assignment]
    SEARCH_TOTAL = _NoOpLabeled()             # type: ignore[assignment]
    SEARCH_LATENCY = _NoOpLabeled()           # type: ignore[assignment]
    SEARCH_RESULTS = _NoOpLabeled()           # type: ignore[assignment]
    EMBEDDING_LATENCY = _NoOpLabeled()        # type: ignore[assignment]
    LLM_CALLS = _NoOpLabeled()               # type: ignore[assignment]
    LLM_LATENCY = _NoOpLabeled()             # type: ignore[assignment]
    LLM_TOKENS = _NoOpLabeled()              # type: ignore[assignment]
    WORKFLOW_TOTAL = _NoOpLabeled()           # type: ignore[assignment]
    WORKFLOW_LATENCY = _NoOpLabeled()         # type: ignore[assignment]
    MATCHING_TOTAL = _NoOpLabeled()           # type: ignore[assignment]
    MATCHING_SCORE = _NoOpLabeled()           # type: ignore[assignment]
    SAFETY_SIGNALS = _NoOpLabeled()           # type: ignore[assignment]
    EXPORT_TOTAL = _NoOpLabeled()             # type: ignore[assignment]
    MILVUS_CONNECTED = _NoOpGauge()           # type: ignore[assignment]
    COLLECTIONS_LOADED = _NoOpGauge()         # type: ignore[assignment]
    COLLECTION_SIZE = _NoOpLabeled()          # type: ignore[assignment]
    ACTIVE_CONNECTIONS = _NoOpGauge()         # type: ignore[assignment]
    AGENT_INFO = _NoOpInfo()                  # type: ignore[assignment]
    INGEST_TOTAL = _NoOpLabeled()             # type: ignore[assignment]
    INGEST_RECORDS = _NoOpLabeled()           # type: ignore[assignment]
    INGEST_ERRORS = _NoOpLabeled()            # type: ignore[assignment]
    INGEST_LATENCY = _NoOpLabeled()           # type: ignore[assignment]
    LAST_INGEST = _NoOpLabeled()              # type: ignore[assignment]
    PIPELINE_STAGE_DURATION = _NoOpLabeled()  # type: ignore[assignment]
    MILVUS_SEARCH_LATENCY = _NoOpLabeled()    # type: ignore[assignment]
    MILVUS_UPSERT_LATENCY = _NoOpLabeled()   # type: ignore[assignment]

    def generate_latest() -> bytes:  # type: ignore[misc]
        return b""


# ===================================================================
# METRICS COLLECTOR (CONVENIENCE WRAPPER)
# ===================================================================


class MetricsCollector:
    """Convenience wrapper for recording Clinical Trial Intelligence Agent metrics.

    Provides static methods that bundle related metric updates into single
    calls, reducing boilerplate in the application code.

    Usage::

        from src.metrics import MetricsCollector

        MetricsCollector.record_query("patient_matching", duration=1.23, success=True)
        MetricsCollector.record_search("trial_protocols", duration=0.15, num_results=12)
        MetricsCollector.record_matching("oncology", score=0.85, phase="phase_iii")
    """

    @staticmethod
    def record_query(workflow_type: str, duration: float, success: bool) -> None:
        """Record metrics for a completed query.

        Args:
            workflow_type: Workflow category (e.g. ``"patient_matching"``).
            duration: Processing time in seconds.
            success: Whether the query completed without error.
        """
        QUERY_TOTAL.labels(workflow_type=workflow_type).inc()
        QUERY_LATENCY.labels(workflow_type=workflow_type).observe(duration)
        if not success:
            QUERY_ERRORS.labels(error_type=workflow_type).inc()

    @staticmethod
    def record_search(
        collection: str, duration: float, num_results: int
    ) -> None:
        """Record metrics for a vector search operation.

        Args:
            collection: Milvus collection name.
            duration: Search time in seconds.
            num_results: Number of results returned.
        """
        SEARCH_TOTAL.labels(collection=collection).inc()
        SEARCH_LATENCY.labels(collection=collection).observe(duration)
        SEARCH_RESULTS.labels(collection=collection).observe(num_results)

    @staticmethod
    def record_embedding(duration: float) -> None:
        """Record embedding generation latency.

        Args:
            duration: Embedding generation time in seconds.
        """
        EMBEDDING_LATENCY.observe(duration)

    @staticmethod
    def record_llm_call(
        model: str,
        duration: float,
        input_tokens: int = 0,
        output_tokens: int = 0,
    ) -> None:
        """Record metrics for an LLM API call.

        Args:
            model: Model identifier (e.g. ``"claude-sonnet-4-6"``).
            duration: Call latency in seconds.
            input_tokens: Number of input tokens.
            output_tokens: Number of output tokens.
        """
        LLM_CALLS.labels(model=model).inc()
        LLM_LATENCY.labels(model=model).observe(duration)
        if input_tokens > 0:
            LLM_TOKENS.labels(direction="input").inc(input_tokens)
        if output_tokens > 0:
            LLM_TOKENS.labels(direction="output").inc(output_tokens)

    @staticmethod
    def record_workflow(workflow_type: str, duration: float) -> None:
        """Record a clinical workflow execution.

        Args:
            workflow_type: Workflow type (e.g. ``"patient_matching"``).
            duration: Execution time in seconds.
        """
        WORKFLOW_TOTAL.labels(workflow_type=workflow_type).inc()
        WORKFLOW_LATENCY.labels(workflow_type=workflow_type).observe(duration)

    @staticmethod
    def record_matching(
        therapeutic_area: str, score: float, phase: str = "unknown"
    ) -> None:
        """Record a patient-trial matching operation.

        Args:
            therapeutic_area: Therapeutic area (e.g. ``"oncology"``).
            score: Match score (0.0 - 1.0).
            phase: Trial phase.
        """
        MATCHING_TOTAL.labels(therapeutic_area=therapeutic_area).inc()
        MATCHING_SCORE.labels(phase=phase).observe(score)

    @staticmethod
    def record_safety_signal(severity: str) -> None:
        """Record a detected safety signal.

        Args:
            severity: Signal severity (e.g. ``"critical"``, ``"high"``).
        """
        SAFETY_SIGNALS.labels(severity=severity).inc()

    @staticmethod
    def record_export(format_type: str) -> None:
        """Record a report export.

        Args:
            format_type: Export format (e.g. ``"markdown"``, ``"pdf"``,
                ``"json"``, ``"fhir_r4"``).
        """
        EXPORT_TOTAL.labels(format=format_type).inc()

    @staticmethod
    def record_ingest(
        source: str,
        duration: float,
        record_count: int,
        collection: str,
        success: bool = True,
    ) -> None:
        """Record an ingest operation.

        Args:
            source: Data source (e.g. ``"clinicaltrials"``, ``"pubmed"``).
            duration: Ingest time in seconds.
            record_count: Number of records ingested.
            collection: Target Milvus collection name.
            success: Whether the ingest completed without error.
        """
        INGEST_TOTAL.labels(source=source).inc()
        INGEST_LATENCY.labels(source=source).observe(duration)
        if success:
            INGEST_RECORDS.labels(collection=collection).inc(record_count)
            LAST_INGEST.labels(source=source).set(time.time())
        else:
            INGEST_ERRORS.labels(source=source).inc()

    @staticmethod
    def set_agent_info(
        version: str, collections: int, workflows: int
    ) -> None:
        """Set agent info gauge with version and configuration.

        Args:
            version: Agent version string.
            collections: Number of active Milvus collections.
            workflows: Number of registered clinical workflows.
        """
        AGENT_INFO.info(
            {
                "version": version,
                "collections": str(collections),
                "workflows": str(workflows),
                "agent": "clinical_trial_intelligence_agent",
            }
        )
        COLLECTIONS_LOADED.set(collections)

    @staticmethod
    def set_milvus_status(connected: bool) -> None:
        """Update Milvus connection status gauge.

        Args:
            connected: Whether Milvus is currently reachable.
        """
        MILVUS_CONNECTED.set(1 if connected else 0)

    @staticmethod
    def update_collection_sizes(stats: Dict[str, int]) -> None:
        """Set the current record count for each collection.

        Args:
            stats: Mapping of collection name to row count.
        """
        for collection, size in stats.items():
            COLLECTION_SIZE.labels(collection=collection).set(size)

    @staticmethod
    def record_pipeline_stage(stage: str, duration: float) -> None:
        """Record duration for a pipeline stage.

        Args:
            stage: Stage name (e.g. ``"embed"``, ``"retrieve"``, ``"synthesize"``).
            duration: Stage duration in seconds.
        """
        PIPELINE_STAGE_DURATION.labels(stage=stage).observe(duration)

    @staticmethod
    def record_milvus_search(duration: float) -> None:
        """Record Milvus vector search latency.

        Args:
            duration: Search time in seconds.
        """
        MILVUS_SEARCH_LATENCY.observe(duration)

    @staticmethod
    def record_milvus_upsert(duration: float) -> None:
        """Record Milvus upsert latency.

        Args:
            duration: Upsert time in seconds.
        """
        MILVUS_UPSERT_LATENCY.observe(duration)


# ===================================================================
# CONVENIENCE FUNCTIONS
# ===================================================================


def get_metrics_text() -> str:
    """Return the current Prometheus metrics exposition in text format.

    Serve this at ``/metrics`` via FastAPI or a dedicated HTTP server.

    Returns:
        UTF-8 decoded Prometheus exposition text, or an empty string if
        ``prometheus_client`` is not installed.
    """
    return generate_latest().decode("utf-8")
