"""Tests for src/metrics.py — Prometheus metrics definitions and helpers.

Tests all 22 metric objects exist, helper functions accept correct arguments,
and get_metrics_text returns a string.

Author: Adam Jones
Date: March 2026
"""

from src import metrics


# ═══════════════════════════════════════════════════════════════════════
# METRIC OBJECTS EXIST
# ═══════════════════════════════════════════════════════════════════════

class TestMetricObjectsExist:
    # -- Histograms (10) --
    def test_query_latency(self):
        assert hasattr(metrics, "QUERY_LATENCY")

    def test_evidence_count(self):
        assert hasattr(metrics, "EVIDENCE_COUNT")

    def test_cross_collection_query_latency(self):
        assert hasattr(metrics, "CROSS_COLLECTION_QUERY_LATENCY")

    def test_cross_collection_results(self):
        assert hasattr(metrics, "CROSS_COLLECTION_RESULTS")

    def test_llm_api_latency(self):
        assert hasattr(metrics, "LLM_API_LATENCY")

    def test_embedding_latency(self):
        assert hasattr(metrics, "EMBEDDING_LATENCY")

    def test_pipeline_stage_duration(self):
        assert hasattr(metrics, "PIPELINE_STAGE_DURATION")

    def test_milvus_search_latency(self):
        assert hasattr(metrics, "MILVUS_SEARCH_LATENCY")

    def test_milvus_upsert_latency(self):
        assert hasattr(metrics, "MILVUS_UPSERT_LATENCY")

    def test_vcf_processing(self):
        assert hasattr(metrics, "VCF_PROCESSING")

    # -- Counters (8) --
    def test_query_count(self):
        assert hasattr(metrics, "QUERY_COUNT")

    def test_collection_hits(self):
        assert hasattr(metrics, "COLLECTION_HITS")

    def test_llm_tokens(self):
        assert hasattr(metrics, "LLM_TOKENS")

    def test_alerts_generated(self):
        assert hasattr(metrics, "ALERTS_GENERATED")

    def test_drug_checks(self):
        assert hasattr(metrics, "DRUG_CHECKS")

    def test_phenoconversions_detected(self):
        assert hasattr(metrics, "PHENOCONVERSIONS_DETECTED")

    def test_hla_screens(self):
        assert hasattr(metrics, "HLA_SCREENS")

    def test_reports_generated(self):
        assert hasattr(metrics, "REPORTS_GENERATED")

    # -- Gauges (4) --
    def test_active_connections(self):
        assert hasattr(metrics, "ACTIVE_CONNECTIONS")

    def test_collection_size(self):
        assert hasattr(metrics, "COLLECTION_SIZE")

    def test_last_ingest(self):
        assert hasattr(metrics, "LAST_INGEST")

    def test_patient_profiles_stored(self):
        assert hasattr(metrics, "PATIENT_PROFILES_STORED")


# ═══════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

class TestHelperFunctions:
    def test_record_query(self):
        # Should not raise
        metrics.record_query("rag", 0.5, 10, "success")

    def test_record_collection_hits(self):
        metrics.record_collection_hits({"pgx_drug_guidelines": 5, "pgx_gene_reference": 3})

    def test_update_collection_sizes(self):
        metrics.update_collection_sizes({"pgx_drug_guidelines": 100})

    def test_record_llm_call(self):
        metrics.record_llm_call("anthropic", "claude-sonnet-4-6", 1.5, 500, 200)

    def test_record_embedding(self):
        metrics.record_embedding(0.05)

    def test_record_pipeline_stage(self):
        metrics.record_pipeline_stage("search", 2.3)

    def test_record_alert_generated(self):
        metrics.record_alert_generated("critical")

    def test_record_drug_check(self):
        metrics.record_drug_check("CYP2D6", "codeine")

    def test_record_phenoconversion(self):
        metrics.record_phenoconversion("CYP2D6", "fluoxetine")

    def test_record_hla_screen(self):
        metrics.record_hla_screen("HLA-B*57:01", "abacavir")

    def test_record_report_generated(self):
        metrics.record_report_generated("pdf")

    def test_record_vcf_processing(self):
        metrics.record_vcf_processing("wgs", 30.0)

    def test_get_metrics_text(self):
        result = metrics.get_metrics_text()
        assert isinstance(result, str)


# ═══════════════════════════════════════════════════════════════════════
# PROMETHEUS AVAILABILITY FLAG
# ═══════════════════════════════════════════════════════════════════════

class TestPrometheusFlag:
    def test_flag_is_bool(self):
        assert isinstance(metrics._PROMETHEUS_AVAILABLE, bool)
