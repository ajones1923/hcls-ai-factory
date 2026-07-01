"""Tests for the Cardiology Intelligence Agent Prometheus metrics module.

Validates that all Prometheus metric objects are created (real or no-op
stubs), that MetricsCollector convenience methods exist and are callable,
and that the get_metrics_text function returns expected output.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import unittest

from src.metrics import (
    ACTIVE_CONNECTIONS,
    AGENT_INFO,
    COLLECTION_SIZE,
    COLLECTIONS_LOADED,
    CRITICAL_ALERTS,
    CROSS_MODAL_TRIGGERS,
    EMBEDDING_LATENCY,
    EXPORT_TOTAL,
    GDMT_OPTIMIZATIONS,
    INGEST_ERRORS,
    INGEST_LATENCY,
    INGEST_RECORDS,
    INGEST_TOTAL,
    LAST_INGEST,
    LLM_CALLS,
    LLM_LATENCY,
    LLM_TOKENS,
    MILVUS_CONNECTED,
    MILVUS_SEARCH_LATENCY,
    MILVUS_UPSERT_LATENCY,
    PIPELINE_STAGE_DURATION,
    QUERY_ERRORS,
    QUERY_LATENCY,
    QUERY_TOTAL,
    RISK_CALC_TOTAL,
    SEARCH_LATENCY,
    SEARCH_RESULTS,
    SEARCH_TOTAL,
    WORKFLOW_LATENCY,
    WORKFLOW_TOTAL,
    MetricsCollector,
    _PROMETHEUS_AVAILABLE,
    get_metrics_text,
)


# =====================================================================
# Metric objects exist
# =====================================================================


class TestMetricObjectsExist(unittest.TestCase):
    """Verify all expected Prometheus metric objects are defined."""

    def test_query_total_exists(self):
        self.assertIsNotNone(QUERY_TOTAL)

    def test_query_latency_exists(self):
        self.assertIsNotNone(QUERY_LATENCY)

    def test_query_errors_exists(self):
        self.assertIsNotNone(QUERY_ERRORS)

    def test_search_total_exists(self):
        self.assertIsNotNone(SEARCH_TOTAL)

    def test_search_latency_exists(self):
        self.assertIsNotNone(SEARCH_LATENCY)

    def test_search_results_exists(self):
        self.assertIsNotNone(SEARCH_RESULTS)

    def test_embedding_latency_exists(self):
        self.assertIsNotNone(EMBEDDING_LATENCY)

    def test_llm_calls_exists(self):
        self.assertIsNotNone(LLM_CALLS)

    def test_llm_latency_exists(self):
        self.assertIsNotNone(LLM_LATENCY)

    def test_llm_tokens_exists(self):
        self.assertIsNotNone(LLM_TOKENS)

    def test_risk_calc_total_exists(self):
        self.assertIsNotNone(RISK_CALC_TOTAL)

    def test_workflow_total_exists(self):
        self.assertIsNotNone(WORKFLOW_TOTAL)

    def test_workflow_latency_exists(self):
        self.assertIsNotNone(WORKFLOW_LATENCY)

    def test_cross_modal_triggers_exists(self):
        self.assertIsNotNone(CROSS_MODAL_TRIGGERS)

    def test_gdmt_optimizations_exists(self):
        self.assertIsNotNone(GDMT_OPTIMIZATIONS)

    def test_critical_alerts_exists(self):
        self.assertIsNotNone(CRITICAL_ALERTS)

    def test_export_total_exists(self):
        self.assertIsNotNone(EXPORT_TOTAL)

    def test_milvus_connected_exists(self):
        self.assertIsNotNone(MILVUS_CONNECTED)

    def test_collections_loaded_exists(self):
        self.assertIsNotNone(COLLECTIONS_LOADED)

    def test_collection_size_exists(self):
        self.assertIsNotNone(COLLECTION_SIZE)

    def test_active_connections_exists(self):
        self.assertIsNotNone(ACTIVE_CONNECTIONS)

    def test_agent_info_exists(self):
        self.assertIsNotNone(AGENT_INFO)

    def test_ingest_total_exists(self):
        self.assertIsNotNone(INGEST_TOTAL)

    def test_ingest_records_exists(self):
        self.assertIsNotNone(INGEST_RECORDS)

    def test_ingest_errors_exists(self):
        self.assertIsNotNone(INGEST_ERRORS)

    def test_ingest_latency_exists(self):
        self.assertIsNotNone(INGEST_LATENCY)

    def test_last_ingest_exists(self):
        self.assertIsNotNone(LAST_INGEST)

    def test_pipeline_stage_duration_exists(self):
        self.assertIsNotNone(PIPELINE_STAGE_DURATION)

    def test_milvus_search_latency_exists(self):
        self.assertIsNotNone(MILVUS_SEARCH_LATENCY)

    def test_milvus_upsert_latency_exists(self):
        self.assertIsNotNone(MILVUS_UPSERT_LATENCY)


# =====================================================================
# MetricsCollector method existence
# =====================================================================


class TestMetricsCollectorMethodsExist(unittest.TestCase):
    """Verify MetricsCollector has all expected static methods."""

    def test_record_query_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_query"))
        self.assertTrue(callable(MetricsCollector.record_query))

    def test_record_search_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_search"))
        self.assertTrue(callable(MetricsCollector.record_search))

    def test_record_embedding_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_embedding"))
        self.assertTrue(callable(MetricsCollector.record_embedding))

    def test_record_llm_call_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_llm_call"))
        self.assertTrue(callable(MetricsCollector.record_llm_call))

    def test_record_risk_calculation_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_risk_calculation"))
        self.assertTrue(callable(MetricsCollector.record_risk_calculation))

    def test_record_workflow_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_workflow"))
        self.assertTrue(callable(MetricsCollector.record_workflow))

    def test_record_cross_modal_trigger_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_cross_modal_trigger"))
        self.assertTrue(callable(MetricsCollector.record_cross_modal_trigger))

    def test_record_gdmt_optimization_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_gdmt_optimization"))
        self.assertTrue(callable(MetricsCollector.record_gdmt_optimization))

    def test_record_critical_alert_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_critical_alert"))
        self.assertTrue(callable(MetricsCollector.record_critical_alert))

    def test_record_export_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_export"))
        self.assertTrue(callable(MetricsCollector.record_export))

    def test_record_ingest_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_ingest"))
        self.assertTrue(callable(MetricsCollector.record_ingest))

    def test_set_agent_info_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "set_agent_info"))
        self.assertTrue(callable(MetricsCollector.set_agent_info))

    def test_set_milvus_status_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "set_milvus_status"))
        self.assertTrue(callable(MetricsCollector.set_milvus_status))

    def test_update_collection_sizes_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "update_collection_sizes"))
        self.assertTrue(callable(MetricsCollector.update_collection_sizes))

    def test_record_pipeline_stage_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_pipeline_stage"))
        self.assertTrue(callable(MetricsCollector.record_pipeline_stage))

    def test_record_milvus_search_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_milvus_search"))
        self.assertTrue(callable(MetricsCollector.record_milvus_search))

    def test_record_milvus_upsert_exists(self):
        self.assertTrue(hasattr(MetricsCollector, "record_milvus_upsert"))
        self.assertTrue(callable(MetricsCollector.record_milvus_upsert))


# =====================================================================
# MetricsCollector invocations (no-op safe)
# =====================================================================


class TestMetricsCollectorInvocations(unittest.TestCase):
    """Verify MetricsCollector methods can be invoked without error."""

    def test_record_query_no_error(self):
        MetricsCollector.record_query("heart_failure", duration=1.0, success=True)

    def test_record_query_failure_no_error(self):
        MetricsCollector.record_query("test", duration=0.5, success=False)

    def test_record_search_no_error(self):
        MetricsCollector.record_search("cardio_literature", duration=0.1, num_results=5)

    def test_record_embedding_no_error(self):
        MetricsCollector.record_embedding(duration=0.05)

    def test_record_llm_call_no_error(self):
        MetricsCollector.record_llm_call("claude-sonnet-4-20250514", duration=2.0, input_tokens=100, output_tokens=50)

    def test_record_risk_calculation_no_error(self):
        MetricsCollector.record_risk_calculation("ascvd")

    def test_record_workflow_no_error(self):
        MetricsCollector.record_workflow("heart_failure", duration=3.0)

    def test_record_cross_modal_trigger_no_error(self):
        MetricsCollector.record_cross_modal_trigger("imaging")

    def test_record_gdmt_optimization_no_error(self):
        MetricsCollector.record_gdmt_optimization()

    def test_record_critical_alert_no_error(self):
        MetricsCollector.record_critical_alert("stemi")

    def test_record_export_no_error(self):
        MetricsCollector.record_export("markdown")

    def test_record_ingest_success_no_error(self):
        MetricsCollector.record_ingest(
            source="pubmed", duration=10.0, record_count=100,
            collection="cardio_literature", success=True,
        )

    def test_record_ingest_failure_no_error(self):
        MetricsCollector.record_ingest(
            source="pubmed", duration=5.0, record_count=0,
            collection="cardio_literature", success=False,
        )

    def test_set_agent_info_no_error(self):
        MetricsCollector.set_agent_info(version="1.0.0", collections=12, workflows=8)

    def test_set_milvus_status_connected(self):
        MetricsCollector.set_milvus_status(connected=True)

    def test_set_milvus_status_disconnected(self):
        MetricsCollector.set_milvus_status(connected=False)

    def test_update_collection_sizes_no_error(self):
        MetricsCollector.update_collection_sizes({"cardio_literature": 5000, "cardio_trials": 200})

    def test_record_pipeline_stage_no_error(self):
        MetricsCollector.record_pipeline_stage("embed", duration=0.5)

    def test_record_milvus_search_no_error(self):
        MetricsCollector.record_milvus_search(duration=0.1)

    def test_record_milvus_upsert_no_error(self):
        MetricsCollector.record_milvus_upsert(duration=0.2)


# =====================================================================
# get_metrics_text
# =====================================================================


class TestGetMetricsText(unittest.TestCase):
    """Tests for get_metrics_text convenience function."""

    def test_returns_string(self):
        result = get_metrics_text()
        self.assertIsInstance(result, str)

    def test_prometheus_available_flag(self):
        self.assertIsInstance(_PROMETHEUS_AVAILABLE, bool)


if __name__ == "__main__":
    unittest.main()
