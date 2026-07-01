"""Tests for FastAPI endpoints (using TestClient patterns).

Tests the export module's output formats which would be served
by the API, plus TestClient-based endpoint tests for the FastAPI app.

Author: Adam Jones
Date: March 2026
"""

import os
import pytest

# Ensure ANTHROPIC_API_KEY is set so settings module loads cleanly
os.environ.setdefault("ANTHROPIC_API_KEY", "test-key-placeholder")

from fastapi.testclient import TestClient

from api.main import app
from src.export import (
    RareDiseaseReportExporter,
    VERSION,
    REPORT_TEMPLATES,
    _now_iso,
    _now_display,
    _generate_filename,
)
from src.metrics import MetricsCollector, get_metrics_text


# ===================================================================
# FASTAPI TESTCLIENT TESTS
# ===================================================================


@pytest.fixture(scope="module")
def client():
    """Create a TestClient that triggers app lifespan."""
    with TestClient(app) as c:
        yield c


class TestHealth:
    """Test /health endpoint."""

    def test_health_returns_200(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert "status" in data
        assert data["status"] in ("healthy", "degraded")

    def test_health_has_components(self, client):
        resp = client.get("/health")
        data = resp.json()
        assert "components" in data
        assert "milvus" in data["components"]
        assert "rag_engine" in data["components"]
        assert "workflow_engine" in data["components"]

    def test_health_status_value(self, client):
        resp = client.get("/health")
        data = resp.json()
        # Status is either healthy or degraded depending on Milvus availability
        assert data["status"] in ("healthy", "degraded")
        assert data["agent"] == "rare-disease-diagnostic-agent"


class TestWorkflows:
    """Test /workflows endpoint."""

    def test_workflows_returns_list(self, client):
        resp = client.get("/workflows")
        assert resp.status_code == 200
        data = resp.json()
        assert "workflows" in data
        assert isinstance(data["workflows"], list)
        assert len(data["workflows"]) > 0
        assert "id" in data["workflows"][0]


class TestCollections:
    """Test /collections endpoint."""

    def test_collections_returns_response(self, client):
        resp = client.get("/collections")
        # Returns 200 with empty list if Milvus reachable, or 503 if not
        assert resp.status_code in (200, 503)


class TestReferenceEndpoints:
    """Test reference data endpoints."""

    def test_disease_categories(self, client):
        resp = client.get("/v1/diagnostic/disease-categories")
        assert resp.status_code == 200
        data = resp.json()
        assert "categories" in data
        assert len(data["categories"]) > 0

    def test_gene_therapies(self, client):
        resp = client.get("/v1/diagnostic/gene-therapies")
        assert resp.status_code == 200
        data = resp.json()
        assert "approved_therapies" in data

    def test_acmg_criteria(self, client):
        resp = client.get("/v1/diagnostic/acmg-criteria")
        assert resp.status_code == 200
        data = resp.json()
        assert "pathogenic_criteria" in data
        assert "benign_criteria" in data

    def test_hpo_categories(self, client):
        resp = client.get("/v1/diagnostic/hpo-categories")
        assert resp.status_code == 200
        data = resp.json()
        assert "categories" in data

    def test_knowledge_version(self, client):
        resp = client.get("/v1/diagnostic/knowledge-version")
        assert resp.status_code == 200
        data = resp.json()
        assert data["agent"] == "rare-disease-diagnostic-agent"
        assert "knowledge_sources" in data


class TestAuth:
    """Test API key authentication middleware."""

    def test_health_exempt_from_auth(self, client):
        # /health should always return 200 even without API key
        resp = client.get("/health")
        assert resp.status_code == 200

    def test_metrics_exempt_from_auth(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200


class TestMetricsEndpoint:
    """Test /metrics endpoint."""

    def test_metrics_returns_text(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200
        assert "rd_agent_" in resp.text or "rd_" in resp.text


# ===================================================================
# EXPORT FORMAT TESTS
# ===================================================================


class TestExportFormats:
    """Test the four export formats."""

    @pytest.fixture
    def exporter(self):
        return RareDiseaseReportExporter()

    def test_markdown_export(self, exporter):
        response = {"findings": ["F1", "F2"], "confidence": 0.9}
        md = exporter.export_markdown(response)
        assert isinstance(md, str)
        assert "F1" in md
        assert "F2" in md
        assert "Disclaimer" in md

    def test_json_export(self, exporter):
        response = {"findings": ["F1"]}
        result = exporter.export_json(response)
        assert isinstance(result, dict)
        assert result["report_type"] == "rare_disease_workflow"
        assert "generated_at" in result
        assert result["version"] == VERSION

    def test_pdf_export_returns_bytes(self, exporter):
        response = {"findings": ["F1"], "recommendations": ["R1"]}
        result = exporter.export_pdf(response)
        assert isinstance(result, bytes)
        # May be empty if reportlab not installed

    def test_fhir_r4_export_structure(self, exporter):
        response = {"findings": ["Finding 1"]}
        bundle = exporter.export_fhir_r4(response, patient_id="P-001")
        assert bundle["resourceType"] == "Bundle"
        assert bundle["type"] == "document"
        report = bundle["entry"][0]["resource"]
        assert report["resourceType"] == "DiagnosticReport"
        assert report["subject"]["reference"] == "Patient/P-001"
        assert report["status"] == "final"


# ===================================================================
# DIAGNOSTIC REPORT TESTS
# ===================================================================


class TestDiagnosticReport:
    """Test diagnostic report generation."""

    @pytest.fixture
    def exporter(self):
        return RareDiseaseReportExporter()

    def test_diagnostic_report_markdown(self, exporter):
        md = exporter.export_diagnostic_report(
            findings=["Elevated CK levels", "Proximal weakness"],
            differential=[
                {"disease_name": "DMD", "score": 0.92, "gene": "DMD"},
                {"disease_name": "BMD", "score": 0.78, "gene": "DMD"},
            ],
            patient_id="P-001",
            phenotypes=["HP:0003236", "HP:0003701"],
        )
        assert "Diagnostic Report" in md
        assert "DMD" in md
        assert "P-001" in md

    def test_diagnostic_report_json(self, exporter):
        result = exporter.export_diagnostic_report(
            findings=["F1"],
            format_type="json",
        )
        assert isinstance(result, dict)
        assert result["report_type"] == "diagnostic"


# ===================================================================
# VARIANT REPORT TESTS
# ===================================================================


class TestVariantReport:
    """Test variant report generation."""

    @pytest.fixture
    def exporter(self):
        return RareDiseaseReportExporter()

    def test_variant_report_markdown(self, exporter):
        md = exporter.export_variant_report(
            variants=[
                {"gene": "CFTR", "variant": "F508del", "classification": "Pathogenic", "evidence": "PS3, PM2"},
            ],
            patient_id="P-002",
        )
        assert "Variant Interpretation" in md
        assert "CFTR" in md
        assert "F508del" in md
        assert "Pathogenic" in md

    def test_variant_report_json(self, exporter):
        result = exporter.export_variant_report(
            variants=[{"gene": "CFTR", "variant": "F508del", "classification": "Pathogenic"}],
            format_type="json",
        )
        assert isinstance(result, dict)
        assert result["variant_count"] == 1


# ===================================================================
# THERAPY REPORT TESTS
# ===================================================================


class TestTherapyReport:
    """Test therapy report generation."""

    @pytest.fixture
    def exporter(self):
        return RareDiseaseReportExporter()

    def test_therapy_report_markdown(self, exporter):
        md = exporter.export_therapy_report(
            therapies=[
                {"drug_name": "Zolgensma", "approval_status": "FDA approved", "eligible": True},
            ],
            patient_id="P-003",
            diagnosis="SMA Type 1",
            gene="SMN1",
        )
        assert "Therapy Eligibility" in md
        assert "Zolgensma" in md
        assert "SMA Type 1" in md

    def test_therapy_report_json(self, exporter):
        result = exporter.export_therapy_report(
            therapies=[{"drug_name": "Zolgensma"}],
            format_type="json",
        )
        assert isinstance(result, dict)
        assert result["therapy_count"] == 1


# ===================================================================
# HELPER FUNCTION TESTS
# ===================================================================


class TestHelperFunctions:
    """Test export helper functions."""

    def test_now_iso_format(self):
        ts = _now_iso()
        assert "T" in ts
        assert ts.endswith("Z")

    def test_now_display_format(self):
        ts = _now_display()
        assert "UTC" in ts

    def test_generate_filename(self):
        fn = _generate_filename("report", "md")
        assert fn.startswith("report_")
        assert fn.endswith(".md")

    def test_report_templates_exist(self):
        assert "diagnostic" in REPORT_TEMPLATES
        assert "variant" in REPORT_TEMPLATES
        assert "therapy" in REPORT_TEMPLATES
        assert "workflow" in REPORT_TEMPLATES

    def test_version_defined(self):
        assert VERSION == "1.0.0"


# ===================================================================
# METRICS TESTS
# ===================================================================


class TestMetrics:
    """Test metrics collector methods."""

    def test_record_query_no_error(self):
        MetricsCollector.record_query("phenotype_driven", 1.0, True)

    def test_record_search_no_error(self):
        MetricsCollector.record_search("rd_phenotypes", 0.1, 10)

    def test_record_export_no_error(self):
        MetricsCollector.record_export("markdown")

    def test_get_metrics_text(self):
        text = get_metrics_text()
        assert isinstance(text, str)
