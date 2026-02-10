"""
Integration smoke tests for HCLS AI Factory.

These tests verify that each portal service starts correctly and responds
to health/readiness probes. They use Flask test clients (no network needed).

Run with: python -m pytest tests/test_integration_smoke.py -v
"""
import json
import pytest
from unittest.mock import patch


@pytest.mark.integration
class TestGenomicsPortalSmoke:
    """Smoke tests for the genomics pipeline portal."""

    def test_health_endpoint(self, genomics_client):
        """Test /health returns 200."""
        response = genomics_client.get('/health')
        assert response.status_code == 200
        data = json.loads(response.data)
        assert 'status' in data

    def test_healthz_endpoint(self, genomics_client):
        """Test /healthz returns 200."""
        response = genomics_client.get('/healthz')
        assert response.status_code == 200

    def test_ready_endpoint_schema(self, genomics_client):
        """Test /api/ready returns correct JSON schema."""
        response = genomics_client.get('/api/ready')
        # May be 200 or 503 depending on environment
        assert response.status_code in (200, 503)
        data = json.loads(response.data)
        assert 'ready' in data
        assert 'checks' in data
        assert isinstance(data['ready'], bool)
        assert isinstance(data['checks'], dict)
        # Verify expected check keys
        assert 'docker' in data['checks']
        assert 'gpu' in data['checks']
        assert 'scripts_dir' in data['checks']

    def test_status_endpoint_returns_json(self, genomics_client, genomics_app):
        """Test /api/status returns JSON with expected structure."""
        with patch.object(genomics_app, 'NVML_AVAILABLE', False):
            response = genomics_client.get('/api/status')
            assert response.status_code == 200
            assert response.content_type == 'application/json'
            data = json.loads(response.data)
            assert 'pipeline' in data
            assert 'system' in data


@pytest.mark.integration
class TestLandingPageSmoke:
    """Smoke tests for the landing page."""

    def test_health_endpoint(self, landing_client):
        """Test /health returns 200."""
        response = landing_client.get('/health')
        assert response.status_code == 200

    def test_ready_endpoint_schema(self, landing_client):
        """Test /api/ready returns correct JSON schema."""
        response = landing_client.get('/api/ready')
        assert response.status_code in (200, 503)
        data = json.loads(response.data)
        assert 'ready' in data
        assert 'checks' in data
        assert isinstance(data['ready'], bool)
        assert isinstance(data['checks'], dict)
