"""Tests for the Cardiology Intelligence Agent ingest scheduler.

Validates CardioScheduler creation, start/stop lifecycle, job listing,
manual ingest triggering, and CardioSchedulerSettings defaults.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import unittest
from unittest.mock import MagicMock

from src.scheduler import (
    CardioScheduler,
    CardioSchedulerSettings,
    IngestJobStatus,
)


# =====================================================================
# CardioSchedulerSettings
# =====================================================================


class TestCardioSchedulerSettings(unittest.TestCase):
    """Tests for the CardioSchedulerSettings dataclass."""

    def test_create_default_settings(self):
        settings = CardioSchedulerSettings()
        self.assertIsInstance(settings, CardioSchedulerSettings)

    def test_default_ingest_enabled(self):
        settings = CardioSchedulerSettings()
        self.assertTrue(settings.INGEST_ENABLED)

    def test_default_schedule_hours(self):
        settings = CardioSchedulerSettings()
        self.assertEqual(settings.INGEST_SCHEDULE_HOURS, 168)

    def test_default_max_pubmed_results(self):
        settings = CardioSchedulerSettings()
        self.assertEqual(settings.MAX_PUBMED_RESULTS, 500)

    def test_default_max_trials_results(self):
        settings = CardioSchedulerSettings()
        self.assertEqual(settings.MAX_TRIALS_RESULTS, 200)

    def test_default_pubmed_query_not_empty(self):
        settings = CardioSchedulerSettings()
        self.assertTrue(len(settings.PUBMED_QUERY) > 0)

    def test_pubmed_query_contains_cardiovascular(self):
        settings = CardioSchedulerSettings()
        self.assertIn("cardiovascular", settings.PUBMED_QUERY.lower())

    def test_default_trials_conditions(self):
        settings = CardioSchedulerSettings()
        self.assertIsInstance(settings.TRIALS_CONDITIONS, list)
        self.assertGreater(len(settings.TRIALS_CONDITIONS), 0)

    def test_trials_conditions_contains_heart_failure(self):
        settings = CardioSchedulerSettings()
        self.assertIn("heart failure", settings.TRIALS_CONDITIONS)

    def test_trials_conditions_contains_cad(self):
        settings = CardioSchedulerSettings()
        self.assertIn("coronary artery disease", settings.TRIALS_CONDITIONS)

    def test_trials_conditions_contains_af(self):
        settings = CardioSchedulerSettings()
        self.assertIn("atrial fibrillation", settings.TRIALS_CONDITIONS)

    def test_default_guideline_sources(self):
        settings = CardioSchedulerSettings()
        self.assertIsInstance(settings.GUIDELINE_SOURCES, list)
        self.assertIn("ACC/AHA", settings.GUIDELINE_SOURCES)
        self.assertIn("ESC", settings.GUIDELINE_SOURCES)

    def test_custom_schedule_hours(self):
        settings = CardioSchedulerSettings(INGEST_SCHEDULE_HOURS=24)
        self.assertEqual(settings.INGEST_SCHEDULE_HOURS, 24)

    def test_custom_ingest_disabled(self):
        settings = CardioSchedulerSettings(INGEST_ENABLED=False)
        self.assertFalse(settings.INGEST_ENABLED)


# =====================================================================
# IngestJobStatus
# =====================================================================


class TestIngestJobStatus(unittest.TestCase):
    """Tests for the IngestJobStatus dataclass."""

    def test_create_status(self):
        status = IngestJobStatus(job_id="test_1", source="pubmed")
        self.assertIsInstance(status, IngestJobStatus)

    def test_default_status_pending(self):
        status = IngestJobStatus(job_id="test_1", source="pubmed")
        self.assertEqual(status.status, "pending")

    def test_default_records_ingested(self):
        status = IngestJobStatus(job_id="test_1", source="pubmed")
        self.assertEqual(status.records_ingested, 0)

    def test_default_error_message_none(self):
        status = IngestJobStatus(job_id="test_1", source="pubmed")
        self.assertIsNone(status.error_message)

    def test_default_duration(self):
        status = IngestJobStatus(job_id="test_1", source="pubmed")
        self.assertEqual(status.duration_seconds, 0.0)


# =====================================================================
# CardioScheduler creation
# =====================================================================


class TestCardioSchedulerCreation(unittest.TestCase):
    """Tests for CardioScheduler instantiation."""

    def test_create_scheduler_no_args(self):
        scheduler = CardioScheduler()
        self.assertIsInstance(scheduler, CardioScheduler)

    def test_create_scheduler_with_settings(self):
        settings = CardioSchedulerSettings(INGEST_ENABLED=False)
        scheduler = CardioScheduler(settings=settings)
        self.assertIsInstance(scheduler, CardioScheduler)

    def test_create_scheduler_with_mocks(self):
        scheduler = CardioScheduler(
            settings=CardioSchedulerSettings(),
            collection_manager=MagicMock(),
            embedder=MagicMock(),
        )
        self.assertIsInstance(scheduler, CardioScheduler)


# =====================================================================
# CardioScheduler start/stop
# =====================================================================


class TestCardioSchedulerStartStop(unittest.TestCase):
    """Tests for CardioScheduler start and stop methods."""

    def test_start_with_ingest_disabled(self):
        settings = CardioSchedulerSettings(INGEST_ENABLED=False)
        scheduler = CardioScheduler(settings=settings)
        # Should complete without error
        scheduler.start()

    def test_stop_without_start(self):
        scheduler = CardioScheduler()
        # Should complete without error even if not started
        scheduler.stop()

    def test_start_exists(self):
        scheduler = CardioScheduler()
        self.assertTrue(hasattr(scheduler, "start"))
        self.assertTrue(callable(scheduler.start))

    def test_stop_exists(self):
        scheduler = CardioScheduler()
        self.assertTrue(hasattr(scheduler, "stop"))
        self.assertTrue(callable(scheduler.stop))


# =====================================================================
# CardioScheduler get_jobs
# =====================================================================


class TestCardioSchedulerGetJobs(unittest.TestCase):
    """Tests for CardioScheduler.get_jobs."""

    def test_get_jobs_returns_list(self):
        scheduler = CardioScheduler()
        jobs = scheduler.get_jobs()
        self.assertIsInstance(jobs, list)

    def test_get_jobs_empty_before_start(self):
        scheduler = CardioScheduler()
        jobs = scheduler.get_jobs()
        self.assertEqual(len(jobs), 0)

    def test_get_jobs_exists(self):
        scheduler = CardioScheduler()
        self.assertTrue(hasattr(scheduler, "get_jobs"))
        self.assertTrue(callable(scheduler.get_jobs))


# =====================================================================
# CardioScheduler get_status
# =====================================================================


class TestCardioSchedulerGetStatus(unittest.TestCase):
    """Tests for CardioScheduler.get_status."""

    def test_get_status_returns_dict(self):
        scheduler = CardioScheduler()
        status = scheduler.get_status()
        self.assertIsInstance(status, dict)

    def test_get_status_has_running_key(self):
        scheduler = CardioScheduler()
        status = scheduler.get_status()
        self.assertIn("running", status)

    def test_get_status_has_ingest_enabled_key(self):
        scheduler = CardioScheduler()
        status = scheduler.get_status()
        self.assertIn("ingest_enabled", status)

    def test_get_status_has_job_count(self):
        scheduler = CardioScheduler()
        status = scheduler.get_status()
        self.assertIn("job_count", status)

    def test_get_status_has_jobs_list(self):
        scheduler = CardioScheduler()
        status = scheduler.get_status()
        self.assertIn("jobs", status)
        self.assertIsInstance(status["jobs"], list)

    def test_get_status_has_recent_history(self):
        scheduler = CardioScheduler()
        status = scheduler.get_status()
        self.assertIn("recent_history", status)


# =====================================================================
# CardioScheduler trigger_manual_ingest
# =====================================================================


class TestCardioSchedulerManualIngest(unittest.TestCase):
    """Tests for CardioScheduler.trigger_manual_ingest."""

    def test_trigger_manual_ingest_exists(self):
        scheduler = CardioScheduler()
        self.assertTrue(hasattr(scheduler, "trigger_manual_ingest"))
        self.assertTrue(callable(scheduler.trigger_manual_ingest))

    def test_trigger_manual_ingest_returns_dict(self):
        scheduler = CardioScheduler()
        result = scheduler.trigger_manual_ingest("pubmed")
        self.assertIsInstance(result, dict)

    def test_trigger_manual_ingest_has_status(self):
        scheduler = CardioScheduler()
        result = scheduler.trigger_manual_ingest("pubmed")
        self.assertIn("status", result)

    def test_trigger_manual_ingest_has_message(self):
        scheduler = CardioScheduler()
        result = scheduler.trigger_manual_ingest("pubmed")
        self.assertIn("message", result)

    def test_trigger_manual_ingest_invalid_source(self):
        scheduler = CardioScheduler()
        result = scheduler.trigger_manual_ingest("nonexistent_source")
        self.assertEqual(result["status"], "error")
        self.assertIn("Unknown source", result["message"])


if __name__ == "__main__":
    unittest.main()
