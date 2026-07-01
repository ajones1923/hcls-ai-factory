"""Prometheus metrics for the Single-Cell Intelligence Agent.

Exposes counters, histograms, gauges, and info metrics for query latency,
collection hits, LLM token usage, workflow executions, cell type analysis,
spatial transcriptomics, TME profiling, trajectory inference, ingest
operations, and system health.

Scraped by the Grafana + Prometheus stack alongside existing HCLS AI Factory
exporters (node_exporter:9100, DCGM:9400).

All metrics use the ``sc_`` prefix so they are easily filterable in
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

    # -- Query Metrics --
    QUERY_TOTAL = Counter(
        "sc_queries_total",
        "Total queries processed",
        ["workflow_type"],
    )

    QUERY_LATENCY = Histogram(
        "sc_query_duration_seconds",
        "Query processing time",
        ["workflow_type"],
        buckets=[0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0],
    )

    QUERY_ERRORS = Counter(
        "sc_query_errors_total",
        "Total query errors",
        ["error_type"],
    )

    # -- RAG / Vector Search Metrics --
    SEARCH_TOTAL = Counter(
        "sc_search_total",
        "Total vector searches",
        ["collection"],
    )

    SEARCH_LATENCY = Histogram(
        "sc_search_duration_seconds",
        "Vector search latency",
        ["collection"],
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0],
    )

    SEARCH_RESULTS = Histogram(
        "sc_search_results_count",
        "Number of results per search",
        ["collection"],
        buckets=[0, 1, 5, 10, 20, 50, 100],
    )

    EMBEDDING_LATENCY = Histogram(
        "sc_embedding_duration_seconds",
        "Embedding generation time",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0],
    )

    # -- LLM Metrics --
    LLM_CALLS = Counter(
        "sc_llm_calls_total",
        "Total LLM calls",
        ["model"],
    )

    LLM_LATENCY = Histogram(
        "sc_llm_duration_seconds",
        "LLM call latency",
        ["model"],
        buckets=[0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 60.0],
    )

    LLM_TOKENS = Counter(
        "sc_llm_tokens_total",
        "Total LLM tokens",
        ["direction"],  # input / output
    )

    # -- Clinical Workflow Metrics --
    WORKFLOW_TOTAL = Counter(
        "sc_workflow_executions_total",
        "Single-cell workflow executions",
        ["workflow_type"],
    )

    WORKFLOW_LATENCY = Histogram(
        "sc_workflow_duration_seconds",
        "Single-cell workflow execution time",
        ["workflow_type"],
        buckets=[0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0],
    )

    # -- Cell Type Analysis Metrics --
    CELL_TYPE_ANALYSES = Counter(
        "sc_cell_type_analyses_total",
        "Cell type annotation analyses performed",
        ["analysis_type"],
    )

    # -- Spatial Transcriptomics Metrics --
    SPATIAL_ANALYSES = Counter(
        "sc_spatial_analyses_total",
        "Spatial transcriptomics analyses performed",
        ["technology"],
    )

    # -- TME Profiling Metrics --
    TME_PROFILES = Counter(
        "sc_tme_profiles_total",
        "TME profiling analyses performed",
        ["cancer_type"],
    )

    # -- Trajectory Inference Metrics --
    TRAJECTORY_ANALYSES = Counter(
        "sc_trajectory_analyses_total",
        "Trajectory inference analyses performed",
        ["method"],
    )

    # -- Drug Response Metrics --
    DRUG_RESPONSE_ANALYSES = Counter(
        "sc_drug_response_analyses_total",
        "Drug response prediction analyses performed",
        ["prediction_type"],
    )

    # -- Export Metrics --
    EXPORT_TOTAL = Counter(
        "sc_exports_total",
        "Report exports",
        ["format"],
    )

    # -- System Metrics --
    MILVUS_CONNECTED = Gauge(
        "sc_milvus_connected",
        "Milvus connection status (1=connected, 0=disconnected)",
    )

    COLLECTIONS_LOADED = Gauge(
        "sc_collections_loaded",
        "Number of loaded collections",
    )

    COLLECTION_SIZE = Gauge(
        "sc_collection_size",
        "Records per collection",
        ["collection"],
    )

    ACTIVE_CONNECTIONS = Gauge(
        "sc_active_connections",
        "Active client connections",
    )

    AGENT_INFO = Info(
        "sc_agent",
        "Agent version and configuration info",
    )

    # -- Ingest Metrics --
    INGEST_TOTAL = Counter(
        "sc_ingest_total",
        "Total ingest operations",
        ["source"],
    )

    INGEST_RECORDS = Counter(
        "sc_ingest_records_total",
        "Total records ingested",
        ["collection"],
    )

    INGEST_ERRORS = Counter(
        "sc_ingest_errors_total",
        "Total ingest errors",
        ["source"],
    )

    INGEST_LATENCY = Histogram(
        "sc_ingest_duration_seconds",
        "Ingest operation time",
        ["source"],
        buckets=[1.0, 5.0, 10.0, 30.0, 60.0, 300.0, 600.0],
    )

    LAST_INGEST = Gauge(
        "sc_last_ingest_timestamp",
        "Last ingest timestamp (unix epoch)",
        ["source"],
    )

    # -- Pipeline Stage Metrics --
    PIPELINE_STAGE_DURATION = Histogram(
        "sc_pipeline_stage_duration_seconds",
        "Duration of individual pipeline stages",
        ["stage"],
        buckets=[0.1, 0.5, 1.0, 5.0, 10.0, 30.0, 60.0, 120.0],
    )

    MILVUS_SEARCH_LATENCY = Histogram(
        "sc_milvus_search_latency_seconds",
        "Milvus vector search latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0],
    )

    MILVUS_UPSERT_LATENCY = Histogram(
        "sc_milvus_upsert_latency_seconds",
        "Milvus upsert latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 5.0],
    )

    _PROMETHEUS_AVAILABLE = True

except ImportError:
    # -- No-op stubs when prometheus_client is not installed --
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
    CELL_TYPE_ANALYSES = _NoOpLabeled()      # type: ignore[assignment]
    SPATIAL_ANALYSES = _NoOpLabeled()        # type: ignore[assignment]
    TME_PROFILES = _NoOpLabeled()            # type: ignore[assignment]
    TRAJECTORY_ANALYSES = _NoOpLabeled()     # type: ignore[assignment]
    DRUG_RESPONSE_ANALYSES = _NoOpLabeled()  # type: ignore[assignment]
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
    """Convenience wrapper for recording Single-Cell Intelligence Agent metrics.

    Provides static methods that bundle related metric updates into single
    calls, reducing boilerplate in the application code.

    Usage::

        from src.metrics import MetricsCollector

        MetricsCollector.record_query("cell_type_annotation", duration=1.23, success=True)
        MetricsCollector.record_search("sc_cell_types", duration=0.15, num_results=12)
        MetricsCollector.record_cell_type_analysis("clustering")
    """

    @staticmethod
    def record_query(workflow_type: str, duration: float, success: bool) -> None:
        """Record metrics for a completed query."""
        QUERY_TOTAL.labels(workflow_type=workflow_type).inc()
        QUERY_LATENCY.labels(workflow_type=workflow_type).observe(duration)
        if not success:
            QUERY_ERRORS.labels(error_type=workflow_type).inc()

    @staticmethod
    def record_search(
        collection: str, duration: float, num_results: int
    ) -> None:
        """Record metrics for a vector search operation."""
        SEARCH_TOTAL.labels(collection=collection).inc()
        SEARCH_LATENCY.labels(collection=collection).observe(duration)
        SEARCH_RESULTS.labels(collection=collection).observe(num_results)

    @staticmethod
    def record_embedding(duration: float) -> None:
        """Record embedding generation latency."""
        EMBEDDING_LATENCY.observe(duration)

    @staticmethod
    def record_llm_call(
        model: str,
        duration: float,
        input_tokens: int = 0,
        output_tokens: int = 0,
    ) -> None:
        """Record metrics for an LLM API call."""
        LLM_CALLS.labels(model=model).inc()
        LLM_LATENCY.labels(model=model).observe(duration)
        if input_tokens > 0:
            LLM_TOKENS.labels(direction="input").inc(input_tokens)
        if output_tokens > 0:
            LLM_TOKENS.labels(direction="output").inc(output_tokens)

    @staticmethod
    def record_workflow(workflow_type: str, duration: float) -> None:
        """Record a single-cell workflow execution."""
        WORKFLOW_TOTAL.labels(workflow_type=workflow_type).inc()
        WORKFLOW_LATENCY.labels(workflow_type=workflow_type).observe(duration)

    @staticmethod
    def record_cell_type_analysis(analysis_type: str) -> None:
        """Record a cell type annotation analysis."""
        CELL_TYPE_ANALYSES.labels(analysis_type=analysis_type).inc()

    @staticmethod
    def record_spatial_analysis(technology: str) -> None:
        """Record a spatial transcriptomics analysis."""
        SPATIAL_ANALYSES.labels(technology=technology).inc()

    @staticmethod
    def record_tme_profile(cancer_type: str) -> None:
        """Record a TME profiling analysis."""
        TME_PROFILES.labels(cancer_type=cancer_type).inc()

    @staticmethod
    def record_trajectory_analysis(method: str) -> None:
        """Record a trajectory inference analysis."""
        TRAJECTORY_ANALYSES.labels(method=method).inc()

    @staticmethod
    def record_drug_response(prediction_type: str) -> None:
        """Record a drug response prediction analysis."""
        DRUG_RESPONSE_ANALYSES.labels(prediction_type=prediction_type).inc()

    @staticmethod
    def record_export(format_type: str) -> None:
        """Record a report export."""
        EXPORT_TOTAL.labels(format=format_type).inc()

    @staticmethod
    def record_ingest(
        source: str,
        duration: float,
        record_count: int,
        collection: str,
        success: bool = True,
    ) -> None:
        """Record an ingest operation."""
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
        """Set agent info gauge with version and configuration."""
        AGENT_INFO.info(
            {
                "version": version,
                "collections": str(collections),
                "workflows": str(workflows),
                "agent": "single_cell_intelligence_agent",
            }
        )
        COLLECTIONS_LOADED.set(collections)

    @staticmethod
    def set_milvus_status(connected: bool) -> None:
        """Update Milvus connection status gauge."""
        MILVUS_CONNECTED.set(1 if connected else 0)

    @staticmethod
    def update_collection_sizes(stats: Dict[str, int]) -> None:
        """Set the current record count for each collection."""
        for collection, size in stats.items():
            COLLECTION_SIZE.labels(collection=collection).set(size)

    @staticmethod
    def record_pipeline_stage(stage: str, duration: float) -> None:
        """Record duration for a pipeline stage."""
        PIPELINE_STAGE_DURATION.labels(stage=stage).observe(duration)

    @staticmethod
    def record_milvus_search(duration: float) -> None:
        """Record Milvus vector search latency."""
        MILVUS_SEARCH_LATENCY.observe(duration)

    @staticmethod
    def record_milvus_upsert(duration: float) -> None:
        """Record Milvus upsert latency."""
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
