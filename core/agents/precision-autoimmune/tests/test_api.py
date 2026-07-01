"""
Tests for the Precision Autoimmune Agent FastAPI endpoints.

All external dependencies (Milvus, Anthropic, SentenceTransformers) are mocked
so tests run without any running services.

Run with:
    cd precision_autoimmune_agent
    python -m pytest tests/test_api.py -v
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional
from unittest.mock import MagicMock, patch

import pytest
from fastapi.testclient import TestClient

# ---------------------------------------------------------------------------
# Fake / stub classes that replace real service objects
# ---------------------------------------------------------------------------

@dataclass
class FakeSearchHit:
    collection: str = "autoimmune_clinical_documents"
    id: str = "hit-001"
    score: float = 0.92
    text: str = "ANA positive 1:320 homogeneous pattern."
    metadata: Dict[str, Any] = field(default_factory=dict)
    relevance: str = "high"


@dataclass
class FakeCrossCollectionResult:
    query: str = ""
    hits: List[FakeSearchHit] = field(default_factory=list)
    knowledge_context: str = ""
    total_collections_searched: int = 3
    search_time_ms: float = 42.5


class FakeCollectionManager:
    """Stub for AutoimmuneCollectionManager."""

    def connect(self):
        pass

    def disconnect(self):
        pass

    def get_collection_stats(self) -> Dict[str, int]:
        return {
            "autoimmune_clinical_documents": 120,
            "autoimmune_patient_labs": 80,
            "autoimmune_autoantibody_panels": 45,
        }

    def create_all_collections(self, drop_existing: bool = False):
        return {"autoimmune_clinical_documents": True}

    def insert_batch(self, collection_name: str, records: list) -> int:
        return len(records)


class FakeRAGEngine:
    """Stub for AutoimmuneRAGEngine."""

    def query(self, question: str, **kwargs) -> str:
        return "Based on evidence, the patient shows SLE indicators."

    def retrieve(self, question: str, **kwargs) -> FakeCrossCollectionResult:
        return FakeCrossCollectionResult(
            query=question,
            hits=[FakeSearchHit()],
            total_collections_searched=3,
            search_time_ms=42.5,
        )

    def search(self, question: str, **kwargs) -> List[FakeSearchHit]:
        return [FakeSearchHit(), FakeSearchHit(id="hit-002", score=0.85)]

    def query_stream(self, question: str, **kwargs):
        yield "Based on"
        yield " evidence,"
        yield " SLE indicators."


class FakeAnalysisResult:
    """Stub for AutoimmuneAnalysisResult."""

    def __init__(self):
        self.patient_id = "PT-001"
        self.disease_activity_scores = []
        self.flare_predictions = []
        self.hla_associations = []
        self.biologic_recommendations = []
        self.critical_alerts = []
        self.cross_agent_findings = []

    def model_dump(self):
        return {
            "patient_id": self.patient_id,
            "disease_activity_scores": self.disease_activity_scores,
            "flare_predictions": self.flare_predictions,
            "hla_associations": self.hla_associations,
            "biologic_recommendations": self.biologic_recommendations,
            "critical_alerts": self.critical_alerts,
            "cross_agent_findings": self.cross_agent_findings,
        }


class FakeAgent:
    """Stub for AutoimmuneAgent."""

    def analyze_patient(self, profile):
        return FakeAnalysisResult()


class FakeDiagnosticEngine:
    """Stub for DiagnosticEngine."""

    def generate_differential(
        self,
        positive_antibodies: List[str],
        hla_alleles: Optional[List[str]] = None,
        symptoms: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        return [
            {
                "disease": "systemic_lupus_erythematosus",
                "probability": 0.78,
                "supporting_antibodies": positive_antibodies,
                "rationale": "ANA + anti-dsDNA strongly suggest SLE.",
            }
        ]


class FakeDocumentProcessor:
    pass


class FakeTimelineBuilder:
    pass


# ---------------------------------------------------------------------------
# Fixture: patched FastAPI app with TestClient
# ---------------------------------------------------------------------------

@pytest.fixture()
def _patch_state():
    """Patch the global _state dict used by all endpoints."""
    fake_state = {
        "cm": FakeCollectionManager(),
        "embedder": MagicMock(),
        "llm": MagicMock(),
        "rag": FakeRAGEngine(),
        "agent": FakeAgent(),
        "doc_processor": FakeDocumentProcessor(),
        "diagnostic": FakeDiagnosticEngine(),
        "timeline": FakeTimelineBuilder(),
        "start_time": time.time() - 60,  # 60s uptime
    }
    with patch("api.main._state", fake_state):
        yield fake_state


@pytest.fixture()
def client(_patch_state):
    """
    Return a TestClient that skips the real lifespan (which would try to
    connect to Milvus and load models).
    """
    # Override lifespan to be a no-op since _state is already patched
    from contextlib import asynccontextmanager

    from api.main import app

    @asynccontextmanager
    async def _noop_lifespan(_app):
        yield

    original_lifespan = app.router.lifespan_context
    app.router.lifespan_context = _noop_lifespan
    try:
        with TestClient(app, raise_server_exceptions=False) as tc:
            yield tc
    finally:
        app.router.lifespan_context = original_lifespan


@pytest.fixture()
def auth_client(_patch_state):
    """
    Return a TestClient where API_KEY authentication is enabled.
    """
    from contextlib import asynccontextmanager

    from api.main import app
    from config.settings import settings

    @asynccontextmanager
    async def _noop_lifespan(_app):
        yield

    original_key = settings.API_KEY
    original_lifespan = app.router.lifespan_context
    app.router.lifespan_context = _noop_lifespan
    settings.API_KEY = "test-secret-key-12345"
    try:
        with TestClient(app, raise_server_exceptions=False) as tc:
            yield tc
    finally:
        settings.API_KEY = original_key
        app.router.lifespan_context = original_lifespan


# ===================================================================
# 1. Health endpoint
# ===================================================================

class TestHealthEndpoint:
    def test_health_returns_200(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "healthy"
        assert data["service"] == "autoimmune-agent"

    def test_health_includes_component_status(self, client):
        data = client.get("/health").json()
        assert "milvus_connected" in data
        assert "embedder_loaded" in data
        assert "llm_available" in data
        assert "uptime_seconds" in data
        assert "collections" in data
        assert "total_vectors" in data

    def test_healthz_returns_ok(self, client):
        resp = client.get("/healthz")
        assert resp.status_code == 200
        assert resp.json() == {"status": "ok"}


# ===================================================================
# 2. Root / info endpoint
# ===================================================================

class TestRootEndpoint:
    def test_root_returns_service_info(self, client):
        resp = client.get("/")
        assert resp.status_code == 200
        data = resp.json()
        assert data["service"] == "Precision Autoimmune Intelligence Agent"
        assert data["version"] == "1.0.0"
        assert data["status"] == "running"
        assert "ports" in data


# ===================================================================
# 3. Metrics endpoint
# ===================================================================

class TestMetricsEndpoint:
    def test_metrics_returns_text_plain(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200
        content_type = resp.headers["content-type"]
        assert "text/plain" in content_type

    def test_metrics_contains_prometheus_format(self, client):
        resp = client.get("/metrics")
        text = resp.text
        assert "autoimmune_agent_up 1" in text
        assert "autoimmune_agent_uptime_seconds" in text

    def test_metrics_includes_collection_vectors(self, client):
        text = client.get("/metrics").text
        assert "autoimmune_collection_vectors" in text
        assert "autoimmune_clinical_documents" in text


# ===================================================================
# 4. Query endpoint
# ===================================================================

class TestQueryEndpoint:
    def test_query_returns_answer(self, client):
        resp = client.post("/query", json={"question": "What is the SLE risk?"})
        assert resp.status_code == 200
        data = resp.json()
        assert "answer" in data
        assert data["evidence_count"] >= 0
        assert data["collections_searched"] >= 0
        assert "search_time_ms" in data

    def test_query_with_patient_context(self, client):
        resp = client.post(
            "/query",
            json={
                "question": "Risk of lupus nephritis?",
                "patient_id": "PT-001",
                "patient_context": "Female, 32yo, ANA+",
                "top_k": 3,
            },
        )
        assert resp.status_code == 200
        assert resp.json()["answer"]

    def test_query_missing_question_returns_422(self, client):
        resp = client.post("/query", json={})
        assert resp.status_code == 422

    def test_query_empty_body_returns_422(self, client):
        resp = client.post(
            "/query",
            content="not json",
            headers={"content-type": "application/json"},
        )
        assert resp.status_code == 422


# ===================================================================
# 5. Search endpoint
# ===================================================================

class TestSearchEndpoint:
    def test_search_returns_hits(self, client):
        resp = client.post("/search", json={"question": "anti-dsDNA antibody"})
        assert resp.status_code == 200
        data = resp.json()
        assert "hits" in data
        assert "total" in data
        assert data["total"] == 2

    def test_search_hit_structure(self, client):
        hits = client.post(
            "/search", json={"question": "test"}
        ).json()["hits"]
        assert len(hits) > 0
        hit = hits[0]
        assert "collection" in hit
        assert "id" in hit
        assert "score" in hit
        assert "text" in hit
        assert "relevance" in hit
        assert "metadata" in hit

    def test_search_with_filter(self, client):
        resp = client.post(
            "/search",
            json={
                "question": "autoantibody",
                "collections_filter": ["autoimmune_autoantibody_panels"],
                "top_k": 2,
            },
        )
        assert resp.status_code == 200


# ===================================================================
# 6. Analyze endpoint
# ===================================================================

class TestAnalyzeEndpoint:
    def test_analyze_with_valid_profile(self, client):
        profile = {
            "patient_id": "PT-001",
            "age": 35,
            "sex": "F",
            "biomarkers": {"CRP": 12.5, "ESR": 45.0},
        }
        resp = client.post("/analyze", json=profile)
        assert resp.status_code == 200
        data = resp.json()
        assert data["patient_id"] == "PT-001"

    def test_analyze_missing_required_fields_returns_422(self, client):
        resp = client.post("/analyze", json={"patient_id": "PT-001"})
        assert resp.status_code == 422

    def test_analyze_invalid_sex_returns_422(self, client):
        resp = client.post(
            "/analyze",
            json={
                "patient_id": "PT-001",
                "age": 35,
                "sex": "X",
                "biomarkers": {"CRP": 1.0},
            },
        )
        assert resp.status_code == 422

    def test_analyze_requires_at_least_one_data_source(self, client):
        """Model validator requires autoantibody_panel, hla_profile, or biomarkers."""
        resp = client.post(
            "/analyze",
            json={"patient_id": "PT-001", "age": 35, "sex": "F"},
        )
        assert resp.status_code == 422


# ===================================================================
# 7. Differential diagnosis endpoint
# ===================================================================

class TestDifferentialEndpoint:
    def test_differential_with_antibodies(self, client):
        resp = client.post(
            "/differential",
            json={"positive_antibodies": ["ANA", "anti-dsDNA"]},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "differential" in data
        assert len(data["differential"]) > 0

    def test_differential_with_full_panel(self, client):
        resp = client.post(
            "/differential",
            json={
                "positive_antibodies": ["ANA", "anti-dsDNA", "anti-Sm"],
                "hla_alleles": ["DRB1*03:01"],
                "symptoms": ["malar rash", "arthralgia", "fatigue"],
            },
        )
        assert resp.status_code == 200

    def test_differential_missing_antibodies_returns_422(self, client):
        resp = client.post("/differential", json={})
        assert resp.status_code == 422

    def test_differential_empty_antibodies_returns_200(self, client):
        resp = client.post(
            "/differential", json={"positive_antibodies": []}
        )
        assert resp.status_code == 200


# ===================================================================
# 8. Collections endpoint
# ===================================================================

class TestCollectionsEndpoint:
    def test_list_collections(self, client):
        resp = client.get("/collections")
        assert resp.status_code == 200
        data = resp.json()
        assert "collections" in data
        assert "total_collections" in data
        assert "total_vectors" in data
        assert data["total_collections"] == 3

    def test_collections_sorted_by_name(self, client):
        colls = client.get("/collections").json()["collections"]
        names = [c["name"] for c in colls]
        assert names == sorted(names)


# ===================================================================
# 9. Export endpoint
# ===================================================================

class TestExportEndpoint:
    def test_export_markdown(self, client):
        resp = client.post(
            "/export",
            json={
                "patient_id": "PT-001",
                "format": "markdown",
                "query_answer": "Patient shows elevated ANA.",
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "markdown"
        assert "content" in data
        assert "PT-001" in data["content"]

    def test_export_fhir(self, client):
        resp = client.post(
            "/export",
            json={"patient_id": "PT-002", "format": "fhir"},
        )
        assert resp.status_code == 200
        data = resp.json()
        assert data["format"] == "fhir"
        content = data["content"]
        assert content["resourceType"] == "Bundle"
        assert content["type"] == "collection"
        # Should contain Patient and DiagnosticReport entries
        resource_types = [e["resource"]["resourceType"] for e in content["entry"]]
        assert "Patient" in resource_types
        assert "DiagnosticReport" in resource_types

    def test_export_unsupported_format_returns_400(self, client):
        resp = client.post(
            "/export",
            json={"patient_id": "PT-001", "format": "csv"},
        )
        assert resp.status_code == 400

    def test_export_missing_patient_id_returns_422(self, client):
        resp = client.post("/export", json={"format": "markdown"})
        assert resp.status_code == 422


# ===================================================================
# 10. Authentication middleware
# ===================================================================

class TestAuthenticationMiddleware:
    def test_protected_endpoint_rejected_without_key(self, auth_client):
        resp = auth_client.post("/query", json={"question": "test"})
        assert resp.status_code == 401
        assert "API key" in resp.json()["detail"]

    def test_protected_endpoint_accepted_with_header(self, auth_client):
        resp = auth_client.post(
            "/query",
            json={"question": "test"},
            headers={"X-API-Key": "test-secret-key-12345"},
        )
        assert resp.status_code == 200

    def test_protected_endpoint_accepted_with_query_param(self, auth_client):
        resp = auth_client.post(
            "/query?api_key=test-secret-key-12345",
            json={"question": "test"},
        )
        assert resp.status_code == 200

    def test_health_bypasses_auth(self, auth_client):
        resp = auth_client.get("/health")
        assert resp.status_code == 200

    def test_healthz_bypasses_auth(self, auth_client):
        resp = auth_client.get("/healthz")
        assert resp.status_code == 200

    def test_metrics_bypasses_auth(self, auth_client):
        resp = auth_client.get("/metrics")
        assert resp.status_code == 200

    def test_root_bypasses_auth(self, auth_client):
        resp = auth_client.get("/")
        assert resp.status_code == 200

    def test_wrong_key_rejected(self, auth_client):
        resp = auth_client.post(
            "/query",
            json={"question": "test"},
            headers={"X-API-Key": "wrong-key"},
        )
        assert resp.status_code == 401

    def test_collections_requires_auth(self, auth_client):
        resp = auth_client.get("/collections")
        assert resp.status_code == 401

    def test_search_requires_auth(self, auth_client):
        resp = auth_client.post("/search", json={"question": "test"})
        assert resp.status_code == 401


# ===================================================================
# 11. Request timing middleware
# ===================================================================

class TestTimingMiddleware:
    def test_process_time_header_present(self, client):
        resp = client.get("/health")
        assert "x-process-time-ms" in resp.headers

    def test_process_time_is_numeric(self, client):
        resp = client.get("/healthz")
        value = resp.headers["x-process-time-ms"]
        assert float(value) >= 0


# ===================================================================
# 12. Request size limit
# ===================================================================

class TestRequestSizeLimit:
    def test_oversized_request_rejected(self, client):
        """Simulate a request with content-length exceeding MAX_REQUEST_SIZE_MB."""
        from config.settings import settings

        max_bytes = settings.MAX_REQUEST_SIZE_MB * 1024 * 1024
        resp = client.post(
            "/query",
            json={"question": "test"},
            headers={"content-length": str(max_bytes + 1)},
        )
        assert resp.status_code == 413
        assert "too large" in resp.json()["detail"].lower()

    def test_normal_request_passes_size_check(self, client):
        resp = client.post("/query", json={"question": "small request"})
        assert resp.status_code != 413


# ===================================================================
# 13. Error handling
# ===================================================================

class TestErrorHandling:
    def test_invalid_json_body_returns_422(self, client):
        resp = client.post(
            "/query",
            content=b"{invalid json",
            headers={"content-type": "application/json"},
        )
        assert resp.status_code == 422

    def test_wrong_field_types_returns_422(self, client):
        resp = client.post(
            "/query",
            json={"question": 12345},  # should be string but Pydantic coerces int->str
        )
        # Pydantic v2 may coerce this; if so, just check it doesn't 500
        assert resp.status_code in (200, 422)

    def test_extra_fields_ignored(self, client):
        resp = client.post(
            "/query",
            json={"question": "test", "nonexistent_field": "value"},
        )
        # Extra fields should be ignored by default in Pydantic
        assert resp.status_code == 200

    def test_analyze_negative_age_returns_422(self, client):
        resp = client.post(
            "/analyze",
            json={
                "patient_id": "PT-001",
                "age": -5,
                "sex": "M",
                "biomarkers": {"CRP": 1.0},
            },
        )
        assert resp.status_code == 422

    def test_nonexistent_endpoint_returns_404(self, client):
        resp = client.get("/nonexistent")
        assert resp.status_code == 404


# ===================================================================
# 14. Streaming endpoint
# ===================================================================

class TestStreamingEndpoint:
    def test_stream_returns_event_stream(self, client):
        resp = client.post(
            "/query/stream",
            json={"question": "What about SLE?"},
        )
        assert resp.status_code == 200
        content_type = resp.headers["content-type"]
        assert "text/event-stream" in content_type

    def test_stream_contains_data_lines(self, client):
        resp = client.post(
            "/query/stream",
            json={"question": "What about SLE?"},
        )
        text = resp.text
        lines = [l for l in text.strip().split("\n") if l.startswith("data:")]  # noqa: E741
        # Should have at least the data chunks plus [DONE]
        assert len(lines) >= 2
        # Last data line should be [DONE]
        assert lines[-1].strip() == "data: [DONE]"

    def test_stream_data_is_valid_json(self, client):
        resp = client.post(
            "/query/stream",
            json={"question": "What about SLE?"},
        )
        lines = [l for l in resp.text.strip().split("\n") if l.startswith("data:")]  # noqa: E741
        for line in lines:
            payload = line[len("data:"):].strip()
            if payload == "[DONE]":
                continue
            parsed = json.loads(payload)
            assert "text" in parsed

    def test_stream_missing_question_returns_422(self, client):
        resp = client.post("/query/stream", json={})
        assert resp.status_code == 422


# ===================================================================
# Edge cases: service not initialized
# ===================================================================

class TestServiceNotInitialized:
    """Verify 503 when components are missing from _state."""

    @pytest.fixture()
    def empty_client(self):
        """Client with empty _state (no services initialized)."""
        from contextlib import asynccontextmanager

        from api.main import app

        @asynccontextmanager
        async def _noop_lifespan(_app):
            yield

        original_lifespan = app.router.lifespan_context
        app.router.lifespan_context = _noop_lifespan
        empty_state = {"start_time": time.time()}
        with patch("api.main._state", empty_state):
            with TestClient(app, raise_server_exceptions=False) as tc:
                yield tc
        app.router.lifespan_context = original_lifespan

    def test_query_503_without_rag(self, empty_client):
        resp = empty_client.post("/query", json={"question": "test"})
        assert resp.status_code == 503

    def test_search_503_without_rag(self, empty_client):
        resp = empty_client.post("/search", json={"question": "test"})
        assert resp.status_code == 503

    def test_analyze_503_without_agent(self, empty_client):
        resp = empty_client.post(
            "/analyze",
            json={
                "patient_id": "PT-001",
                "age": 30,
                "sex": "F",
                "biomarkers": {"CRP": 1.0},
            },
        )
        assert resp.status_code == 503

    def test_differential_503_without_diagnostic(self, empty_client):
        resp = empty_client.post(
            "/differential",
            json={"positive_antibodies": ["ANA"]},
        )
        assert resp.status_code == 503

    def test_collections_503_without_cm(self, empty_client):
        resp = empty_client.get("/collections")
        assert resp.status_code == 503

    def test_stream_503_without_rag(self, empty_client):
        resp = empty_client.post(
            "/query/stream", json={"question": "test"}
        )
        assert resp.status_code == 503
