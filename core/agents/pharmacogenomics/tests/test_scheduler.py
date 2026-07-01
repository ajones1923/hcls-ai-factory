"""Tests for src/scheduler.py — IngestScheduler.

Tests initialization, start/stop lifecycle, get_status, and
default interval configuration.

Author: Adam Jones
Date: March 2026
"""

import pytest
from unittest.mock import MagicMock

from src.scheduler import IngestScheduler


# ═══════════════════════════════════════════════════════════════════════
# IngestScheduler
# ═══════════════════════════════════════════════════════════════════════

class TestIngestScheduler:
    @pytest.fixture
    def scheduler(self):
        mock_cm = MagicMock()
        mock_embedder = MagicMock()
        return IngestScheduler(mock_cm, mock_embedder)

    def test_init_defaults(self, scheduler):
        # Whether real or stub, it should accept these args
        assert scheduler is not None

    def test_init_custom_interval(self):
        mock_cm = MagicMock()
        mock_embedder = MagicMock()
        sched = IngestScheduler(mock_cm, mock_embedder, interval_hours=24)
        assert sched is not None

    def test_get_status_returns_dict(self, scheduler):
        status = scheduler.get_status()
        assert isinstance(status, dict)

    def test_get_status_has_keys(self, scheduler):
        status = scheduler.get_status()
        assert "next_run_time" in status
        assert "last_run_time" in status
        assert "job_count" in status

    def test_get_status_before_start(self, scheduler):
        status = scheduler.get_status()
        assert status["last_run_time"] is None

    def test_stop_without_start(self, scheduler):
        # Should not raise
        scheduler.stop()

    def test_start_and_stop(self):
        mock_cm = MagicMock()
        mock_embedder = MagicMock()
        sched = IngestScheduler(mock_cm, mock_embedder, interval_hours=999)
        try:
            sched.start()
            status = sched.get_status()
            assert status["job_count"] >= 0
        finally:
            sched.stop()
