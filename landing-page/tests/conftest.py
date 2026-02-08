"""Shared pytest fixtures for Landing Page tests."""
import pytest
from pathlib import Path
from unittest.mock import patch
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))
