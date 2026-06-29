"""Shared test fixtures for Neurology Intelligence Agent tests."""

import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so that imports like
# ``from src.models import ...`` and ``from config.settings import ...``
# resolve correctly regardless of how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
