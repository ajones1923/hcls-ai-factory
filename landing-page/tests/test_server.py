"""Tests for AI Factory landing page server."""
import pytest
import json
from pathlib import Path
from unittest.mock import patch, Mock
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from server import app, is_port_open, check_service_health, SERVICES


@pytest.fixture
def client():
    """Create Flask test client."""
    app.config['TESTING'] = True
    with app.test_client() as client:
        yield client


class TestHealthEndpoint:
    """Tests for /health endpoint."""

    def test_returns_200(self, client):
        resp = client.get('/health')
        assert resp.status_code == 200

    def test_returns_healthy_status(self, client):
        resp = client.get('/health')
        data = json.loads(resp.data)
        assert data['status'] == 'healthy'
        assert data['service'] == 'landing-page'


class TestReportStatus:
    """Tests for /api/report-status endpoint."""

    def test_returns_json(self, client):
        resp = client.get('/api/report-status')
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert 'generated_at' in data

    def test_returns_generated_at_key(self, client):
        resp = client.get('/api/report-status')
        data = json.loads(resp.data)
        # Should always return generated_at (may be None if no report exists)
        assert 'generated_at' in data


class TestCheckSingleService:
    """Tests for /api/check-service/<id> endpoint."""

    def test_unknown_service_returns_404(self, client):
        resp = client.get('/api/check-service/nonexistent')
        assert resp.status_code == 404
        data = json.loads(resp.data)
        assert 'error' in data

    @patch('server.get_host_ip', return_value='127.0.0.1')
    @patch('server.check_service_health')
    def test_known_service_returns_result(self, mock_check, mock_ip, client):
        mock_check.return_value = {
            'id': 'genomics', 'name': 'Genomics Pipeline',
            'port': 5000, 'online': False, 'message': 'Connection refused'
        }
        resp = client.get('/api/check-service/genomics')
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert data['id'] == 'genomics'


class TestCheckAllServices:
    """Tests for /api/check-services endpoint."""

    @patch('server.get_host_ip', return_value='127.0.0.1')
    @patch('server.check_service_health')
    def test_returns_summary(self, mock_check, mock_ip, client):
        mock_check.return_value = {
            'id': 'test', 'name': 'Test', 'port': 5000,
            'online': False, 'message': 'Offline'
        }
        resp = client.get('/api/check-services')
        assert resp.status_code == 200
        data = json.loads(resp.data)
        assert 'services' in data
        assert 'summary' in data
        assert 'online' in data['summary']
        assert 'total' in data['summary']
        assert 'all_healthy' in data['summary']


class TestHelperFunctions:
    """Tests for utility functions."""

    def test_is_port_open_closed_port(self):
        # Port 1 should not be open on localhost
        assert is_port_open(1, timeout=0.1) is False

    def test_check_service_health_unreachable(self):
        result = check_service_health(
            'test',
            {'port': 1, 'name': 'Test Service', 'path': '/health'},
            host='127.0.0.1',
            timeout=0.1,
        )
        assert result['online'] is False
        assert result['id'] == 'test'
        assert result['name'] == 'Test Service'

    def test_check_service_health_tcp_only(self):
        result = check_service_health(
            'test-tcp',
            {'port': 1, 'name': 'TCP Service', 'path': None},
            host='127.0.0.1',
            timeout=0.1,
        )
        assert result['online'] is False

    def test_services_dict_has_expected_keys(self):
        assert 'genomics' in SERVICES
        assert 'milvus' in SERVICES
        assert 'grafana' in SERVICES
        for sid, sinfo in SERVICES.items():
            assert 'port' in sinfo
            assert 'name' in sinfo
