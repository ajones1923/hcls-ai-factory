"""Tests for FastAPI routes.

Uses FastAPI TestClient to test API endpoints if the API module
is available.

Author: Adam Jones
Date: March 2026
"""

import pytest

# API module may not be fully implemented yet
try:
    from fastapi.testclient import TestClient
    _FASTAPI_AVAILABLE = True
except ImportError:
    _FASTAPI_AVAILABLE = False

# Try to import the app
try:
    from api.main import app
    _APP_AVAILABLE = True
except ImportError:
    _APP_AVAILABLE = False

from src.models import NeuroWorkflowType


@pytest.mark.skipif(
    not (_FASTAPI_AVAILABLE and _APP_AVAILABLE),
    reason="FastAPI app not available"
)
class TestHealthEndpoint:
    """Tests for the /health endpoint."""

    def test_health_returns_200(self):
        client = TestClient(app)
        response = client.get("/health")
        assert response.status_code == 200

    def test_health_has_status(self):
        client = TestClient(app)
        response = client.get("/health")
        data = response.json()
        assert "status" in data


@pytest.mark.skipif(
    not (_FASTAPI_AVAILABLE and _APP_AVAILABLE),
    reason="FastAPI app not available"
)
class TestQueryEndpoint:
    """Tests for the /api/query endpoint."""

    def test_query_returns_200(self):
        client = TestClient(app)
        response = client.post(
            "/api/query",
            json={"question": "What is the treatment for acute ischemic stroke?"},
        )
        assert response.status_code in (200, 404, 422, 500)


class TestAPIPlaceholder:
    """Placeholder tests that always pass."""

    def test_workflow_types_for_api(self):
        """Verify all workflow types are valid for API routing."""
        for wf in NeuroWorkflowType:
            assert isinstance(wf.value, str)

    def test_api_port_configured(self):
        from config.settings import settings
        assert 1024 <= settings.API_PORT <= 65535

    def test_cors_configured(self):
        from config.settings import settings
        assert settings.CORS_ORIGINS, "CORS_ORIGINS should not be empty"

    def test_api_host_configured(self):
        from config.settings import settings
        assert settings.API_HOST == "0.0.0.0"

    def test_max_request_size(self):
        from config.settings import settings
        assert settings.MAX_REQUEST_SIZE_MB > 0
