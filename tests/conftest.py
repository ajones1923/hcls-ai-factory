"""
Shared fixtures for cross-service integration smoke tests.

Each portal has a server.py, so we use importlib.util to load them
under unique module names to avoid conflicts.
"""
import importlib.util
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

REPO_ROOT = Path(__file__).parent.parent


def _load_server(module_name, server_path):
    """Load a server.py under a unique module name."""
    # Add the server's directory to sys.path for its own imports
    server_dir = str(server_path.parent)
    if server_dir not in sys.path:
        sys.path.insert(0, server_dir)

    spec = importlib.util.spec_from_file_location(module_name, server_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def genomics_app():
    """Load the genomics portal Flask app."""
    server_path = REPO_ROOT / 'genomics-pipeline' / 'web-portal' / 'app' / 'server.py'
    module = _load_server('genomics_server', server_path)
    module.app.config['TESTING'] = True
    return module


@pytest.fixture
def genomics_client(genomics_app):
    """Flask test client for the genomics portal."""
    with genomics_app.app.test_client() as client:
        yield client


@pytest.fixture(scope="session")
def landing_app():
    """Load the landing page Flask app."""
    server_path = REPO_ROOT / 'landing-page' / 'server.py'
    module = _load_server('landing_server', server_path)
    module.app.config['TESTING'] = True
    return module


@pytest.fixture
def landing_client(landing_app):
    """Flask test client for the landing page."""
    with landing_app.app.test_client() as client:
        yield client
