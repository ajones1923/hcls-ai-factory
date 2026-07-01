"""Tests for ImagingIngestScheduler and scheduler lifecycle.

Validates:
  - INGEST_ENABLED is False by default
  - INGEST_SCHEDULE_HOURS defaults to 168 (weekly)
  - INGEST_TRIALS_SCHEDULE_HOURS defaults to 720 (monthly)
  - ImagingIngestScheduler imports cleanly
  - Scheduler registers both pubmed and trials jobs
  - Scheduler does not start when INGEST_ENABLED is False
  - Manual ingest cycle runs and returns stats
  - Schedule info is accessible

Author: Adam Jones
Date: April 2026
"""

import pytest
from unittest.mock import MagicMock, patch

from config.settings import ImagingSettings
from src.scheduler import (
    APSCHEDULER_AVAILABLE,
    ImagingIngestScheduler,
    PUBMED_INTERVAL_HOURS,
    TRIALS_INTERVAL_HOURS,
    start_scheduler,
    stop_scheduler,
    get_scheduler,
)


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def mock_pubmed_ingestor():
    """Return a mock PubMed ingestor that returns a record count."""
    ingestor = MagicMock(return_value=42)
    return ingestor


@pytest.fixture
def mock_trials_ingestor():
    """Return a mock trials ingestor that returns a record count."""
    ingestor = MagicMock(return_value=15)
    return ingestor


@pytest.fixture
def scheduler(mock_pubmed_ingestor, mock_trials_ingestor, mock_embedder, mock_collection_manager):
    """Return an ImagingIngestScheduler with mock ingestors."""
    sched = ImagingIngestScheduler(
        pubmed_ingestor=mock_pubmed_ingestor,
        trials_ingestor=mock_trials_ingestor,
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
    )
    yield sched
    if sched.is_running:
        sched.stop()


# ===================================================================
# TESTS
# ===================================================================


def test_scheduler_disabled_by_default():
    """INGEST_ENABLED should be False in default settings."""
    settings = ImagingSettings()
    assert settings.INGEST_ENABLED is False


def test_scheduler_config_hours():
    """INGEST_SCHEDULE_HOURS should default to 168 (weekly)."""
    settings = ImagingSettings()
    assert settings.INGEST_SCHEDULE_HOURS == 168


def test_scheduler_trials_config_hours():
    """INGEST_TRIALS_SCHEDULE_HOURS should default to 720 (monthly)."""
    settings = ImagingSettings()
    assert settings.INGEST_TRIALS_SCHEDULE_HOURS == 720


def test_scheduler_importable():
    """ImagingIngestScheduler should import cleanly."""
    assert ImagingIngestScheduler is not None
    assert callable(ImagingIngestScheduler)


def test_scheduler_default_intervals():
    """Module-level interval constants should match expected values."""
    assert PUBMED_INTERVAL_HOURS == 168
    assert TRIALS_INTERVAL_HOURS == 720


@pytest.mark.skipif(not APSCHEDULER_AVAILABLE, reason="apscheduler not installed")
def test_scheduler_jobs_registered(scheduler):
    """Scheduler should register both pubmed and trials jobs on start."""
    scheduler.start(pubmed_interval_hours=168, trials_interval_hours=720)

    job_ids = scheduler.job_ids
    assert "pubmed_ingest" in job_ids
    assert "trials_ingest" in job_ids
    assert len(job_ids) == 2


@pytest.mark.skipif(not APSCHEDULER_AVAILABLE, reason="apscheduler not installed")
def test_scheduler_does_not_start_when_disabled():
    """start_scheduler() should return None when INGEST_ENABLED is False."""
    with patch("src.scheduler.settings", create=True) as mock_settings:
        mock_settings.INGEST_ENABLED = False
        result = start_scheduler()
        assert result is None


@pytest.mark.skipif(not APSCHEDULER_AVAILABLE, reason="apscheduler not installed")
def test_scheduler_start_stop(scheduler):
    """Scheduler should start and stop cleanly."""
    assert scheduler.is_running is False

    scheduler.start()
    assert scheduler.is_running is True

    scheduler.stop()
    assert scheduler.is_running is False


@pytest.mark.skipif(not APSCHEDULER_AVAILABLE, reason="apscheduler not installed")
def test_scheduler_prevents_double_start(scheduler):
    """Starting an already-running scheduler should log a warning, not crash."""
    scheduler.start()
    assert scheduler.is_running is True

    # Second start should not raise
    scheduler.start()
    assert scheduler.is_running is True


def test_scheduler_run_ingest_cycle(scheduler, mock_pubmed_ingestor, mock_trials_ingestor):
    """Manual ingest cycle should call both ingestors and return stats."""
    stats = scheduler.run_ingest_cycle()

    assert stats["pubmed_records"] == 42
    assert stats["trials_records"] == 15
    assert stats["total_records"] == 57
    assert stats["success"] is True
    assert stats["duration_seconds"] >= 0
    assert len(stats["errors"]) == 0

    mock_pubmed_ingestor.assert_called_once()
    mock_trials_ingestor.assert_called_once()


def test_scheduler_run_ingest_cycle_with_error(mock_embedder, mock_collection_manager):
    """Ingest cycle should capture errors without raising."""
    failing_ingestor = MagicMock(side_effect=RuntimeError("API timeout"))

    sched = ImagingIngestScheduler(
        pubmed_ingestor=failing_ingestor,
        trials_ingestor=MagicMock(return_value=5),
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
    )

    stats = sched.run_ingest_cycle()

    assert stats["pubmed_records"] == 0
    assert stats["trials_records"] == 5
    assert stats["success"] is False
    assert len(stats["errors"]) == 1
    assert "PubMed" in stats["errors"][0]


@pytest.mark.skipif(not APSCHEDULER_AVAILABLE, reason="apscheduler not installed")
def test_scheduler_schedule_info(scheduler):
    """get_schedule_info() should return schedule details."""
    scheduler.start(pubmed_interval_hours=24, trials_interval_hours=48)

    info = scheduler.get_schedule_info()
    assert info["running"] is True
    assert info["pubmed_interval_hours"] == 24
    assert info["trials_interval_hours"] == 48
    assert len(info["jobs"]) == 2


def test_scheduler_last_run_stats_empty(scheduler):
    """last_run_stats should be empty before any run."""
    assert scheduler.last_run_stats == {}


def test_scheduler_last_run_stats_updated(scheduler):
    """last_run_stats should be populated after a run."""
    scheduler.run_ingest_cycle()
    stats = scheduler.last_run_stats
    assert stats["total_records"] == 57
    assert stats["success"] is True
