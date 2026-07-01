"""Prometheus metrics for the Rare Disease Diagnostic Agent.

Exposes counters, histograms, gauges, and info metrics for query latency,
collection hits, LLM token usage, workflow executions, diagnostic scoring,
phenotype matching, ingest operations, and system health.

All metrics use the ``rd_`` prefix so they are easily filterable in
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
        "rd_queries_total",
        "Total queries processed",
        ["workflow_type"],
    )

    QUERY_LATENCY = Histogram(
        "rd_query_duration_seconds",
        "Query processing time",
        ["workflow_type"],
        buckets=[0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0],
    )

    QUERY_ERRORS = Counter(
        "rd_query_errors_total",
        "Total query errors",
        ["error_type"],
    )

    # -- RAG / Vector Search Metrics --
    SEARCH_TOTAL = Counter(
        "rd_search_total",
        "Total vector searches",
        ["collection"],
    )

    SEARCH_LATENCY = Histogram(
        "rd_search_duration_seconds",
        "Vector search latency",
        ["collection"],
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0],
    )

    SEARCH_RESULTS = Histogram(
        "rd_search_results_count",
        "Number of results per search",
        ["collection"],
        buckets=[0, 1, 5, 10, 20, 50, 100],
    )

    EMBEDDING_LATENCY = Histogram(
        "rd_embedding_duration_seconds",
        "Embedding generation time",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0],
    )

    # -- LLM Metrics --
    LLM_CALLS = Counter(
        "rd_llm_calls_total",
        "Total LLM calls",
        ["model"],
    )

    LLM_LATENCY = Histogram(
        "rd_llm_duration_seconds",
        "LLM call latency",
        ["model"],
        buckets=[0.5, 1.0, 2.0, 5.0, 10.0, 30.0, 60.0],
    )

    LLM_TOKENS = Counter(
        "rd_llm_tokens_total",
        "Total LLM tokens",
        ["direction"],  # input / output
    )

    # -- Diagnostic Workflow Metrics --
    WORKFLOW_TOTAL = Counter(
        "rd_workflow_executions_total",
        "Diagnostic workflow executions",
        ["workflow_type"],
    )

    WORKFLOW_LATENCY = Histogram(
        "rd_workflow_duration_seconds",
        "Diagnostic workflow execution time",
        ["workflow_type"],
        buckets=[0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0],
    )

    # -- Phenotype Matching Metrics --
    PHENOTYPE_MATCHES = Counter(
        "rd_phenotype_matches_total",
        "Total phenotype matching operations",
        ["disease_category"],
    )

    PHENOTYPE_SCORE = Histogram(
        "rd_phenotype_match_score",
        "Phenotype match similarity scores",
        ["method"],
        buckets=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    )

    # -- ACMG Classification Metrics --
    ACMG_CLASSIFICATIONS = Counter(
        "rd_acmg_classifications_total",
        "ACMG variant classifications",
        ["classification"],
    )

    # -- Gene Therapy Eligibility Metrics --
    THERAPY_ELIGIBILITY = Counter(
        "rd_therapy_eligibility_checks_total",
        "Gene therapy eligibility checks",
        ["therapy"],
    )

    # -- Export Metrics --
    EXPORT_TOTAL = Counter(
        "rd_exports_total",
        "Report exports",
        ["format"],
    )

    # -- System Metrics --
    MILVUS_CONNECTED = Gauge(
        "rd_milvus_connected",
        "Milvus connection status (1=connected, 0=disconnected)",
    )

    COLLECTIONS_LOADED = Gauge(
        "rd_collections_loaded",
        "Number of loaded collections",
    )

    COLLECTION_SIZE = Gauge(
        "rd_collection_size",
        "Records per collection",
        ["collection"],
    )

    ACTIVE_CONNECTIONS = Gauge(
        "rd_active_connections",
        "Active client connections",
    )

    AGENT_INFO = Info(
        "rd_agent",
        "Agent version and configuration info",
    )

    # -- Ingest Metrics --
    INGEST_TOTAL = Counter(
        "rd_ingest_total",
        "Total ingest operations",
        ["source"],
    )

    INGEST_RECORDS = Counter(
        "rd_ingest_records_total",
        "Total records ingested",
        ["collection"],
    )

    INGEST_ERRORS = Counter(
        "rd_ingest_errors_total",
        "Total ingest errors",
        ["source"],
    )

    INGEST_LATENCY = Histogram(
        "rd_ingest_duration_seconds",
        "Ingest operation time",
        ["source"],
        buckets=[1.0, 5.0, 10.0, 30.0, 60.0, 300.0, 600.0],
    )

    LAST_INGEST = Gauge(
        "rd_last_ingest_timestamp",
        "Last ingest timestamp (unix epoch)",
        ["source"],
    )

    # -- Pipeline Stage Metrics --
    PIPELINE_STAGE_DURATION = Histogram(
        "rd_pipeline_stage_duration_seconds",
        "Duration of individual pipeline stages",
        ["stage"],
        buckets=[0.1, 0.5, 1.0, 5.0, 10.0, 30.0, 60.0, 120.0],
    )

    MILVUS_SEARCH_LATENCY = Histogram(
        "rd_milvus_search_latency_seconds",
        "Milvus vector search latency",
        buckets=[0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0],
    )

    MILVUS_UPSERT_LATENCY = Histogram(
        "rd_milvus_upsert_latency_seconds",
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
    PHENOTYPE_MATCHES = _NoOpLabeled()        # type: ignore[assignment]
    PHENOTYPE_SCORE = _NoOpLabeled()          # type: ignore[assignment]
    ACMG_CLASSIFICATIONS = _NoOpLabeled()     # type: ignore[assignment]
    THERAPY_ELIGIBILITY = _NoOpLabeled()      # type: ignore[assignment]
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
    """Convenience wrapper for recording Rare Disease Diagnostic Agent metrics.

    Provides static methods that bundle related metric updates into single
    calls, reducing boilerplate in the application code.

    Usage::

        from src.metrics import MetricsCollector

        MetricsCollector.record_query("phenotype_driven", duration=1.23, success=True)
        MetricsCollector.record_search("rd_phenotypes", duration=0.15, num_results=12)
        MetricsCollector.record_phenotype_match("metabolic", score=0.85, method="ic_score")
    """

    @staticmethod
    def record_query(workflow_type: str, duration: float, success: bool) -> None:
        """Record metrics for a completed query."""
        QUERY_TOTAL.labels(workflow_type=workflow_type).inc()
        QUERY_LATENCY.labels(workflow_type=workflow_type).observe(duration)
        if not success:
            QUERY_ERRORS.labels(error_type=workflow_type).inc()

    @staticmethod
    def record_search(collection: str, duration: float, num_results: int) -> None:
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
        """Record a diagnostic workflow execution."""
        WORKFLOW_TOTAL.labels(workflow_type=workflow_type).inc()
        WORKFLOW_LATENCY.labels(workflow_type=workflow_type).observe(duration)

    @staticmethod
    def record_phenotype_match(
        disease_category: str, score: float, method: str = "semantic"
    ) -> None:
        """Record a phenotype matching operation."""
        PHENOTYPE_MATCHES.labels(disease_category=disease_category).inc()
        PHENOTYPE_SCORE.labels(method=method).observe(score)

    @staticmethod
    def record_acmg_classification(classification: str) -> None:
        """Record an ACMG variant classification."""
        ACMG_CLASSIFICATIONS.labels(classification=classification).inc()

    @staticmethod
    def record_therapy_eligibility(therapy: str) -> None:
        """Record a gene therapy eligibility check."""
        THERAPY_ELIGIBILITY.labels(therapy=therapy).inc()

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
    def set_agent_info(version: str, collections: int, workflows: int) -> None:
        """Set agent info gauge with version and configuration."""
        AGENT_INFO.info(
            {
                "version": version,
                "collections": str(collections),
                "workflows": str(workflows),
                "agent": "rare_disease_diagnostic_agent",
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
    """Return the current Prometheus metrics exposition in text format."""
    return generate_latest().decode("utf-8")
