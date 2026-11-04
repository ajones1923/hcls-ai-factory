"""
Shared pytest fixtures for RAG Chat Pipeline Portal tests.
"""
import pytest
import tempfile
from pathlib import Path
from unittest.mock import patch
import sys

from flask.testing import FlaskClient

# Add the app directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))


# This portal always enforces an API key: with PORTAL_API_KEY unset it generates an
# ephemeral one rather than opening up. The suite used to send no key at all and so ran
# entirely against 401s. Tests now run with a known key configured; `anon_client` sends
# none, for the tests that assert an unauthenticated caller is refused.
PORTAL_TEST_KEY = 'portal-test-key'


class _KeyedClient(FlaskClient):
    """Test client that sends X-API-Key unless the caller set one explicitly."""

    def open(self, *args, **kwargs):
        headers = kwargs.get('headers') or {}
        if isinstance(headers, dict) and not any(k.lower() == 'x-api-key' for k in headers):
            kwargs['headers'] = {**headers, 'X-API-Key': PORTAL_TEST_KEY}
        return super().open(*args, **kwargs)


@pytest.fixture(autouse=True)
def portal_api_key(monkeypatch):
    """Configure a known portal key for every test."""
    monkeypatch.setenv('PORTAL_API_KEY', PORTAL_TEST_KEY)
    return PORTAL_TEST_KEY


@pytest.fixture
def app():
    """Create Flask test app."""
    from server import app as flask_app
    flask_app.config['TESTING'] = True
    return flask_app


@pytest.fixture
def client(app):
    """Create test client that authenticates."""
    app.test_client_class = _KeyedClient
    try:
        return app.test_client()
    finally:
        app.test_client_class = None


@pytest.fixture
def anon_client(app):
    """Create test client that sends no credential."""
    app.test_client_class = None
    return app.test_client()


@pytest.fixture
def temp_data_dir():
    """Create temporary data directory structure."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        (tmpdir / 'logs').mkdir()
        (tmpdir / 'targets').mkdir()
        yield tmpdir
