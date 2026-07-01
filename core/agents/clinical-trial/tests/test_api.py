"""FastAPI TestClient tests for the Clinical Trial Intelligence Agent API.

Tests all core endpoints, authentication, rate limiting, and reference
catalogues using the FastAPI TestClient with a no-op lifespan to avoid
Milvus/LLM dependencies.

Author: Adam Jones
Date: March 2026
"""

import sys
import threading
from contextlib import asynccontextmanager
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ---------------------------------------------------------------------------
# No-op lifespan override (skip Milvus / LLM / embedding startup)
# ---------------------------------------------------------------------------
@asynccontextmanager
async def _noop_lifespan(app):
    yield


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def client():
    """Create a TestClient with no-op lifespan and cleared state."""
    from fastapi.testclient import TestClient

    from api.main import _rate_limit_store, app

    _rate_limit_store.clear()
    app.router.lifespan_context = _noop_lifespan
    app.state.engine = None
    app.state.manager = None
    app.state.workflow_engine = None
    app.state.llm_client = None
    app.state.metrics = {
        "requests_total": 0,
        "query_requests_total": 0,
        "search_requests_total": 0,
        "match_requests_total": 0,
        "workflow_requests_total": 0,
        "protocol_requests_total": 0,
        "report_requests_total": 0,
        "errors_total": 0,
    }
    app.state.metrics_lock = threading.Lock()
    with TestClient(app, raise_server_exceptions=False) as c:
        yield c


# ===================================================================
# TestHealth
# ===================================================================
class TestHealth:
    """Tests for GET /health."""

    def test_health_returns_200(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200

    def test_health_degraded_without_milvus(self, client):
        resp = client.get("/health")
        data = resp.json()
        assert data["status"] == "degraded"

    def test_health_has_components(self, client):
        resp = client.get("/health")
        data = resp.json()
        assert "components" in data
        assert "milvus" in data["components"]
        assert "rag_engine" in data["components"]
        assert "workflow_engine" in data["components"]

    def test_health_agent_name(self, client):
        resp = client.get("/health")
        data = resp.json()
        assert data["agent"] == "clinical-trial-intelligence-agent"


# ===================================================================
# TestWorkflows
# ===================================================================
class TestWorkflows:
    """Tests for GET /workflows."""

    def test_workflows_returns_200(self, client):
        resp = client.get("/workflows")
        assert resp.status_code == 200

    def test_workflows_correct_count(self, client):
        resp = client.get("/workflows")
        data = resp.json()
        # 11 workflows: 10 specific + 1 general
        assert len(data["workflows"]) == 11

    def test_workflows_have_required_fields(self, client):
        resp = client.get("/workflows")
        data = resp.json()
        for wf in data["workflows"]:
            assert "id" in wf, f"Missing 'id' in workflow: {wf}"
            assert "name" in wf, f"Missing 'name' in workflow: {wf}"
            assert "description" in wf, f"Missing 'description' in workflow: {wf}"


# ===================================================================
# TestMetrics
# ===================================================================
class TestMetrics:
    """Tests for GET /metrics."""

    def test_metrics_returns_200(self, client):
        resp = client.get("/metrics")
        assert resp.status_code == 200

    def test_metrics_contains_counter_format(self, client):
        resp = client.get("/metrics")
        text = resp.text
        # May use prometheus_client (# TYPE ... counter) or custom format
        assert "# TYPE" in text
        assert "counter" in text


# ===================================================================
# TestCollections
# ===================================================================
class TestCollections:
    """Tests for GET /collections."""

    def test_collections_503_without_milvus(self, client):
        resp = client.get("/collections")
        assert resp.status_code == 503


# ===================================================================
# TestReferenceEndpoints
# ===================================================================
class TestReferenceEndpoints:
    """Tests for reference catalogue endpoints under /v1/trial/."""

    def test_therapeutic_areas_returns_list(self, client):
        resp = client.get("/v1/trial/therapeutic-areas")
        assert resp.status_code == 200
        data = resp.json()
        assert "therapeutic_areas" in data
        assert isinstance(data["therapeutic_areas"], list)
        assert len(data["therapeutic_areas"]) > 0

    def test_phases_returns_list(self, client):
        resp = client.get("/v1/trial/phases")
        assert resp.status_code == 200
        data = resp.json()
        assert "phases" in data
        assert isinstance(data["phases"], list)
        assert len(data["phases"]) >= 6  # at least 6 standard phases

    def test_guidelines_returns_list(self, client):
        resp = client.get("/v1/trial/guidelines")
        assert resp.status_code == 200
        data = resp.json()
        assert "guidelines" in data
        assert isinstance(data["guidelines"], list)
        assert len(data["guidelines"]) > 0

    def test_knowledge_version_returns_dict(self, client):
        resp = client.get("/v1/trial/knowledge-version")
        assert resp.status_code == 200
        data = resp.json()
        assert "version" in data
        assert "agent" in data
        assert data["agent"] == "clinical-trial-intelligence-agent"
        assert "collections" in data
        assert "workflows" in data


# ===================================================================
# TestQueryEndpoints
# ===================================================================
class TestQueryEndpoints:
    """Tests for POST query/search/protocol endpoints."""

    def test_query_503_without_engine(self, client):
        resp = client.post(
            "/v1/trial/query",
            json={"question": "What is the best endpoint for oncology trials?"},
        )
        assert resp.status_code == 503

    def test_search_503_without_engine(self, client):
        resp = client.post(
            "/v1/trial/search",
            json={"question": "PFS endpoint oncology phase III"},
        )
        assert resp.status_code == 503

    def test_protocol_optimize_returns_response(self, client):
        """Protocol optimize uses local logic and should work without engine."""
        resp = client.post(
            "/v1/trial/protocol/optimize",
            json={
                "protocol_summary": "A phase III double-blind trial of Drug X vs placebo in NSCLC",
                "therapeutic_area": "oncology",
                "phase": "phase_iii",
                "indication": "NSCLC",
                "endpoints": ["PFS", "OS"],
                "eligibility_criteria": ["Age >= 18", "ECOG 0-1", "Measurable disease"],
                "visit_count": 12,
                "procedure_count": 15,
            },
        )
        assert resp.status_code == 200
        data = resp.json()
        assert "complexity_score" in data
        assert "optimization_recommendations" in data
        assert isinstance(data["optimization_recommendations"], list)
        assert len(data["optimization_recommendations"]) > 0
        assert data["complexity_score"] >= 0


# ===================================================================
# TestAuth
# ===================================================================
class TestAuth:
    """Tests for API key authentication middleware."""

    def test_unauthenticated_returns_401(self, client):
        """When API_KEY is set, unauthenticated request should return 401."""
        from config.settings import settings

        original = settings.API_KEY
        try:
            settings.API_KEY = "test-secret-key-12345"
            resp = client.get("/workflows")
            assert resp.status_code == 401
        finally:
            settings.API_KEY = original

    def test_health_exempt_from_auth(self, client):
        """Health endpoint should be accessible without API key."""
        from config.settings import settings

        original = settings.API_KEY
        try:
            settings.API_KEY = "test-secret-key-12345"
            resp = client.get("/health")
            assert resp.status_code == 200
        finally:
            settings.API_KEY = original

    def test_valid_key_passes(self, client):
        """Request with correct API key should succeed."""
        from config.settings import settings

        original = settings.API_KEY
        try:
            settings.API_KEY = "test-secret-key-12345"
            resp = client.get(
                "/workflows",
                headers={"X-API-Key": "test-secret-key-12345"},
            )
            assert resp.status_code == 200
        finally:
            settings.API_KEY = original

    def test_metrics_exempt_from_auth(self, client):
        """Metrics endpoint should be accessible without API key."""
        from config.settings import settings

        original = settings.API_KEY
        try:
            settings.API_KEY = "test-secret-key-12345"
            resp = client.get("/metrics")
            assert resp.status_code == 200
        finally:
            settings.API_KEY = original


# ===================================================================
# TestRateLimiting
# ===================================================================
class TestRateLimiting:
    """Tests for rate limiting middleware."""

    def test_rate_limit_triggers_429(self, client):
        """101 rapid requests to a non-exempt path should trigger 429."""
        from api.main import _rate_limit_store

        _rate_limit_store.clear()
        last_status = None
        for i in range(101):
            resp = client.get("/workflows")
            last_status = resp.status_code
            if last_status == 429:
                break
        assert last_status == 429, (
            f"Expected 429 after exceeding rate limit, got {last_status}"
        )
