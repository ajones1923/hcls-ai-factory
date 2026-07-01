"""Comprehensive API endpoint tests for Imaging Intelligence Agent.

Tests all major endpoint groups using FastAPI TestClient with mocked
backend services (Milvus, NIM, embeddings). No running server or
external services required.

Author: Adam Jones
Date: April 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

# Ensure project root on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =====================================================================
# TestClient — patch heavy init so the app starts without Milvus/NIM
# =====================================================================


@pytest.fixture(scope="module")
def client():
    """Create a TestClient with mocked lifespan dependencies.

    Patches SentenceTransformer, ImagingCollectionManager, NIMServiceManager,
    ImagingRAGEngine, and CrossModalTrigger so the FastAPI app starts cleanly
    without any external services.
    """
    import numpy as np
    from src.models import CrossCollectionResult, SearchHit

    mock_embedder = MagicMock()
    mock_embedder.encode.side_effect = lambda texts, **kw: (
        np.random.randn(384).astype(np.float32)
        if isinstance(texts, str)
        else np.random.randn(len(texts), 384).astype(np.float32)
    )

    mock_manager = MagicMock()
    mock_manager.get_collection_stats.return_value = {
        "imaging_literature": 100,
        "imaging_trials": 80,
        "imaging_findings": 60,
        "imaging_protocols": 50,
        "imaging_devices": 40,
        "imaging_anatomy": 30,
        "imaging_benchmarks": 70,
        "imaging_guidelines": 55,
        "imaging_report_templates": 20,
        "imaging_datasets": 25,
        "genomic_evidence": 90,
    }

    sample_evidence = CrossCollectionResult(
        query="test query",
        hits=[
            SearchHit(
                collection="imaging_literature",
                id="lit-001",
                score=0.85,
                text="AI-based hemorrhage detection achieves 95% sensitivity.",
                metadata={"modality": "ct"},
            ),
            SearchHit(
                collection="imaging_trials",
                id="trial-001",
                score=0.78,
                text="Phase 3 trial of AI-assisted lung cancer screening.",
                metadata={"phase": "Phase 3"},
            ),
        ],
        knowledge_context="## Pathology: Intracranial Hemorrhage",
        total_collections_searched=11,
        search_time_ms=42.0,
    )

    mock_engine = MagicMock()
    mock_engine.retrieve.return_value = sample_evidence
    mock_engine.query.return_value = (
        "Based on the available imaging evidence, the findings are consistent "
        "with normal anatomy. Clinical correlation recommended."
    )
    mock_engine.find_related.return_value = {
        "imaging_literature": sample_evidence.hits[:1],
    }

    mock_nim = MagicMock()
    mock_nim.check_all_services.return_value = {
        "vista3d": "mock",
        "maisi": "mock",
        "vila_m3": "mock",
        "llm": "mock",
    }

    cross_modal_mock = MagicMock()
    cross_modal_mock.evaluate.return_value = None

    patches = [
        patch("sentence_transformers.SentenceTransformer", return_value=mock_embedder),
        patch("src.collections.ImagingCollectionManager", return_value=mock_manager),
        patch("src.nim.service_manager.NIMServiceManager", return_value=mock_nim),
        patch("src.rag_engine.ImagingRAGEngine", return_value=mock_engine),
        patch("src.cross_modal.CrossModalTrigger", return_value=cross_modal_mock),
    ]
    for p in patches:
        p.start()

    from api.main import app, _state

    _state["manager"] = mock_manager
    _state["embedder"] = mock_embedder
    _state["nim_manager"] = mock_nim
    _state["engine"] = mock_engine
    _state["settings"] = MagicMock(
        DICOM_AUTO_INGEST=False,
        ORTHANC_URL="http://localhost:8042",
        DICOM_WATCH_INTERVAL=30,
        NIM_CLOUD_URL="",
        NIM_VISTA3D_URL="",
        NIM_MAISI_URL="",
        NIM_VILAM3_URL="",
        NIM_LLM_URL="",
        API_KEY="",
        MAX_REQUEST_SIZE_MB=10,
        CORS_ORIGINS="*",
    )
    _state["cross_modal_trigger"] = cross_modal_mock

    with TestClient(app) as tc:
        yield tc

    for p in patches:
        p.stop()
    _state.clear()


# =====================================================================
# Health & Status (3 tests)
# =====================================================================


class TestHealthAndStatus:
    """Health, metrics, and collection endpoints."""

    def test_health_returns_200(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] in ("healthy", "degraded", "initializing")

    def test_metrics_returns_200(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200
        assert "imaging_agent" in resp.text or "HELP" in resp.text

    def test_collections_returns_list(self, client):
        resp = client.get("/collections")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        if data:
            assert "name" in data[0]
            assert "count" in data[0]


# =====================================================================
# Workflows (6 tests)
# =====================================================================


class TestWorkflows:
    """Workflow list, info, and run endpoints."""

    def test_list_workflows(self, client):
        resp = client.get("/workflows")
        assert resp.status_code == 200
        data = resp.json()
        assert data["total"] == 9
        assert len(data["workflows"]) == 9

    def test_workflow_info(self, client):
        resp = client.get("/workflow/ct_head_hemorrhage/info")
        assert resp.status_code == 200
        data = resp.json()
        assert data["name"] == "ct_head_hemorrhage"
        assert "modality" in data
        assert "body_region" in data

    def test_run_workflow(self, client):
        resp = client.post(
            "/workflow/ct_head_hemorrhage/run",
            json={"mock_mode": True},
        )
        assert resp.status_code == 200, f"Workflow run failed: {resp.text}"
        data = resp.json()
        assert data["status"] == "completed"
        assert data["is_mock"] is True

    def test_run_all_workflows(self, client):
        from src.workflows import WORKFLOW_REGISTRY

        for name in WORKFLOW_REGISTRY:
            resp = client.post(f"/workflow/{name}/run", json={"mock_mode": True})
            assert resp.status_code == 200, f"Workflow {name} failed: {resp.text}"
            data = resp.json()
            assert data["status"] == "completed", f"Workflow {name} status: {data['status']}"

    def test_workflow_not_found(self, client):
        resp = client.get("/workflow/nonexistent_workflow/info")
        assert resp.status_code == 404

    def test_workflow_result_has_classification(self, client):
        resp = client.post(
            "/workflow/ct_head_hemorrhage/run",
            json={"mock_mode": True},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "classification" in data
        assert len(data["classification"]) > 0


# =====================================================================
# Evidence / RAG (4 tests)
# =====================================================================


class TestEvidenceRAG:
    """Meta-agent /api/ask and /search endpoints."""

    def test_ask_returns_answer(self, client):
        resp = client.post(
            "/api/ask",
            json={"question": "What AI models detect intracranial hemorrhage?"},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "answer" in data
        assert len(data["answer"]) > 0

    def test_ask_has_evidence(self, client):
        resp = client.post(
            "/api/ask",
            json={"question": "What AI models detect intracranial hemorrhage?"},
        )
        data = resp.json()
        assert data["evidence_count"] > 0

    def test_ask_has_follow_ups(self, client):
        resp = client.post(
            "/api/ask",
            json={
                "question": "What AI models detect intracranial hemorrhage?",
                "include_follow_ups": True,
            },
        )
        data = resp.json()
        assert "follow_up_questions" in data
        assert len(data["follow_up_questions"]) > 0

    def test_search_returns_results(self, client):
        resp = client.post(
            "/search",
            json={"question": "Lung nodule CT screening guidelines"},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["total_hits"] > 0
        assert len(data["hits"]) > 0


# =====================================================================
# Demo Cases (3 tests)
# =====================================================================


class TestDemoCases:
    """Demo case list and run endpoints."""

    def test_list_demo_cases(self, client):
        resp = client.get("/demo-cases")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        assert len(data) >= 4

    def test_run_demo_case(self, client):
        resp = client.post("/demo-cases/DEMO-001/run")
        assert resp.status_code == 200
        data = resp.json()
        assert data["case_id"] == "DEMO-001"
        assert "workflow_result" in data

    def test_demo_case_not_found(self, client):
        resp = client.post("/demo-cases/NONEXISTENT/run")
        assert resp.status_code == 404


# =====================================================================
# Protocol Optimizer (3 tests)
# =====================================================================


class TestProtocol:
    """Protocol recommendation endpoints."""

    def test_recommend_protocol(self, client):
        resp = client.post(
            "/protocol/recommend",
            json={"indication": "headache_acute"},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "recommended_modality" in data
        assert "recommended_protocol" in data

    def test_indications_list(self, client):
        resp = client.get("/protocol/indications")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        assert len(data) >= 12

    def test_protocol_has_acr_rating(self, client):
        resp = client.post(
            "/protocol/recommend",
            json={"indication": "headache_acute"},
        )
        data = resp.json()
        assert "acr_appropriateness_rating" in data
        assert isinstance(data["acr_appropriateness_rating"], int)
        assert 1 <= data["acr_appropriateness_rating"] <= 9


# =====================================================================
# Dose Intelligence (4 tests)
# =====================================================================


class TestDoseIntelligence:
    """Dose recording, cumulative tracking, and DRL comparison endpoints."""

    def test_record_dose(self, client):
        resp = client.post(
            "/dose/record",
            json={
                "patient_id": "PAT-TEST-001",
                "study_date": "2026-04-01",
                "modality": "CT",
                "protocol": "ct_head_routine",
                "body_region": "head",
                "effective_dose_msv": 2.0,
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "recorded"

    def test_cumulative_dose(self, client):
        # Record a dose first
        client.post(
            "/dose/record",
            json={
                "patient_id": "PAT-CUMULATIVE-001",
                "study_date": "2026-04-01",
                "modality": "CT",
                "protocol": "ct_head_routine",
                "body_region": "head",
                "effective_dose_msv": 2.0,
            },
        )
        resp = client.get("/dose/cumulative/PAT-CUMULATIVE-001")
        assert resp.status_code == 200
        data = resp.json()
        assert "total_effective_dose_msv" in data
        assert data["patient_id"] == "PAT-CUMULATIVE-001"

    def test_drl_comparison(self, client):
        resp = client.post(
            "/dose/compare-drl",
            json={
                "protocol": "ct_head_routine",
                "body_region": "head",
                "effective_dose_msv": 2.5,
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "ratio" in data
        assert "status" in data

    def test_population_dose(self, client):
        resp = client.get("/dose/population")
        assert resp.status_code == 200
        data = resp.json()
        assert "total_records" in data or "by_protocol" in data or isinstance(data, dict)


# =====================================================================
# Analytics (5 tests)
# =====================================================================


class TestAnalytics:
    """Population analytics endpoints."""

    def test_generate_demo_data(self, client):
        resp = client.post("/api/analytics/generate-demo-data?n_studies=50")
        assert resp.status_code == 200
        data = resp.json()
        assert "generated" in data
        assert data["generated"] >= 50

    def test_population_summary(self, client):
        # Ensure data exists
        client.post("/api/analytics/generate-demo-data?n_studies=50")
        resp = client.get("/api/analytics/population")
        assert resp.status_code == 200
        data = resp.json()
        assert "total_studies" in data

    def test_cohort_query(self, client):
        client.post("/api/analytics/generate-demo-data?n_studies=50")
        resp = client.post(
            "/api/analytics/cohort",
            json={"modality": "CT"},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "matching_studies" in data or "count" in data or isinstance(data, dict)

    def test_trends(self, client):
        client.post("/api/analytics/generate-demo-data?n_studies=100")
        resp = client.get("/api/analytics/trends/study_count")
        assert resp.status_code == 200
        data = resp.json()
        assert "data_points" in data or "metric" in data or isinstance(data, dict)

    def test_severity_by_modality(self, client):
        client.post("/api/analytics/generate-demo-data?n_studies=100")
        resp = client.get("/api/analytics/severity-by-modality")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, dict)


# =====================================================================
# NIM Services (2 tests)
# =====================================================================


class TestNIMServices:
    """NIM service status endpoints."""

    def test_nim_status(self, client):
        resp = client.get("/nim/status")
        assert resp.status_code == 200
        data = resp.json()
        assert "services" in data
        assert isinstance(data["services"], list)

    def test_nim_services_have_status(self, client):
        resp = client.get("/nim/status")
        data = resp.json()
        for svc in data["services"]:
            assert "name" in svc
            assert "status" in svc
            assert svc["status"] in (
                "available", "cloud", "anthropic", "mock", "unavailable",
            )


# =====================================================================
# Reports (2 tests)
# =====================================================================


class TestReports:
    """Report generation endpoints."""

    def test_generate_markdown_report(self, client):
        resp = client.post(
            "/reports/generate",
            json={
                "question": "What is the sensitivity of AI for hemorrhage detection?",
                "format": "markdown",
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "content" in data
        assert len(data["content"]) > 0
        assert "#" in data["content"]  # Markdown heading

    def test_report_has_evidence(self, client):
        resp = client.post(
            "/reports/generate",
            json={
                "question": "What is the sensitivity of AI for hemorrhage detection?",
                "include_evidence": True,
            },
        )
        data = resp.json()
        assert data["evidence_count"] > 0


# =====================================================================
# Knowledge (2 tests)
# =====================================================================


class TestKnowledge:
    """Knowledge graph stats endpoint."""

    def test_knowledge_stats(self, client):
        resp = client.get("/knowledge/stats")
        assert resp.status_code == 200
        data = resp.json()
        assert "pathologies" in data
        assert "modalities" in data

    def test_knowledge_has_pathologies(self, client):
        resp = client.get("/knowledge/stats")
        data = resp.json()
        assert data["pathologies"] >= 20


# =====================================================================
# Events (2 tests)
# =====================================================================


class TestEvents:
    """DICOM event bus endpoints."""

    def test_events_status(self, client):
        resp = client.get("/events/status")
        assert resp.status_code == 200
        data = resp.json()
        assert "active" in data or "supported_events" in data

    def test_events_history(self, client):
        resp = client.get("/events/history")
        assert resp.status_code == 200
        data = resp.json()
        assert "items" in data
        assert isinstance(data["items"], list)
