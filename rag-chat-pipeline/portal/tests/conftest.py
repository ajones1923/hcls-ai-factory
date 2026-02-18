"""
Shared pytest fixtures for RAG Chat Pipeline Portal tests.
"""
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

# Add the app directory to the path
sys.path.insert(0, str(Path(__file__).parent.parent / 'app'))


@pytest.fixture
def app():
    """Create Flask test app."""
    from server import app as flask_app
    flask_app.config['TESTING'] = True
    return flask_app


@pytest.fixture
def client(app):
    """Create test client."""
    return app.test_client()


@pytest.fixture
def temp_data_dir():
    """Create temporary data directory structure."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        (tmpdir / 'logs').mkdir()
        (tmpdir / 'targets').mkdir()
        yield tmpdir
