"""
Automated ingest scheduler for the Oncology Intelligence Agent.

Wraps APScheduler's BackgroundScheduler to periodically refresh data from
PubMed, ClinicalTrials.gov, and CIViC into the agent's Milvus collections.
Falls back to a no-op stub when ``apscheduler`` is not installed.

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

import logging
import time
from typing import Any, Dict, Optional

logger = logging.getLogger(__name__)

try:
    from apscheduler.schedulers.background import BackgroundScheduler

    _APSCHEDULER_AVAILABLE = True
except ImportError:
    _APSCHEDULER_AVAILABLE = False

# Import metrics (works even as no-ops)
from .metrics import LAST_INGEST


# -----------------------------------------------------------------------
# Main scheduler implementation
# -----------------------------------------------------------------------

class IngestScheduler:
    """Periodic data ingest scheduler.

    Manages recurring background jobs that refresh the agent's knowledge
    base from external sources:

    * **PubMed** -- recent oncology literature
    * **ClinicalTrials.gov** -- active/recruiting clinical trials
    * **CIViC** -- curated clinical interpretation of variants

    Parameters
    ----------
    collection_manager:
        Provides access to Milvus collection CRUD operations.
    embedder:
        Embedding model/service used by the ingest pipelines.
    interval_hours:
        Hours between successive runs for each job. Defaults to 168
        (one week).
    """

    def __init__(
        self,
        collection_manager,
        embedder,
        interval_hours: int = 168,
    ) -> None:
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.interval_hours = interval_hours

        if not _APSCHEDULER_AVAILABLE:
            logger.warning(
                "apscheduler not installed – IngestScheduler will be a no-op. "
                "Install with: pip install apscheduler"
            )
            self._scheduler: Optional[BackgroundScheduler] = None
        else:
            self._scheduler = BackgroundScheduler(
                job_defaults={
                    "coalesce": True,
                    "max_instances": 1,
                    "misfire_grace_time": 3600,
                },
            )

        self._last_run: Dict[str, Optional[float]] = {
            "pubmed": None,
            "clinical_trials": None,
            "civic": None,
        }

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start(self) -> None:
        """Start the scheduler and register all recurring ingest jobs."""
        if self._scheduler is None:
            logger.warning("Scheduler not available – start() is a no-op")
            return

        self._scheduler.add_job(
            self._refresh_pubmed,
            trigger="interval",
            hours=self.interval_hours,
            id="refresh_pubmed",
            name="Refresh PubMed literature",
            next_run_time=None,  # do not fire immediately
        )

        self._scheduler.add_job(
            self._refresh_clinical_trials,
            trigger="interval",
            hours=self.interval_hours,
            id="refresh_clinical_trials",
            name="Refresh ClinicalTrials.gov",
            next_run_time=None,
        )

        self._scheduler.add_job(
            self._refresh_civic,
            trigger="interval",
            hours=self.interval_hours,
            id="refresh_civic",
            name="Refresh CIViC evidence",
            next_run_time=None,
        )

        self._scheduler.start()
        logger.info(
            "IngestScheduler started with %d jobs (interval=%dh)",
            len(self._scheduler.get_jobs()),
            self.interval_hours,
        )

    def stop(self) -> None:
        """Gracefully shut down the scheduler, waiting for running jobs."""
        if self._scheduler is None:
            return

        logger.info("IngestScheduler shutting down...")
        self._scheduler.shutdown(wait=True)
        logger.info("IngestScheduler stopped")

    # ------------------------------------------------------------------
    # Status
    # ------------------------------------------------------------------

    def get_status(self) -> Dict[str, Any]:
        """Return scheduler status information.

        Returns
        -------
        dict
            Keys: ``next_run_time``, ``last_run_time``, ``job_count``.
        """
        if self._scheduler is None:
            return {
                "next_run_time": None,
                "last_run_time": self._last_run,
                "job_count": 0,
            }

        jobs = self._scheduler.get_jobs()

        next_run_times = {}
        for job in jobs:
            nrt = job.next_run_time
            next_run_times[job.id] = nrt.isoformat() if nrt else None

        return {
            "next_run_time": next_run_times,
            "last_run_time": self._last_run,
            "job_count": len(jobs),
        }

    # ------------------------------------------------------------------
    # Ingest jobs
    # ------------------------------------------------------------------

    def _refresh_pubmed(self) -> None:
        """Refresh PubMed oncology literature."""
        logger.info("Starting PubMed ingest refresh")
        try:
            from .ingest.literature_parser import PubMedIngestPipeline

            pipeline = PubMedIngestPipeline(
                collection_manager=self.collection_manager,
                embedder=self.embedder,
            )
            pipeline.run()

            now = time.time()
            self._last_run["pubmed"] = now
            LAST_INGEST.labels(source="pubmed").set(now)
            logger.info("PubMed ingest refresh completed")

        except Exception:
            logger.exception("PubMed ingest refresh failed")

    def _refresh_clinical_trials(self) -> None:
        """Refresh ClinicalTrials.gov data."""
        logger.info("Starting ClinicalTrials.gov ingest refresh")
        try:
            from .ingest.clinical_trials_parser import ClinicalTrialsIngestPipeline

            pipeline = ClinicalTrialsIngestPipeline(
                collection_manager=self.collection_manager,
                embedder=self.embedder,
            )
            pipeline.run()

            now = time.time()
            self._last_run["clinical_trials"] = now
            LAST_INGEST.labels(source="clinical_trials").set(now)
            logger.info("ClinicalTrials.gov ingest refresh completed")

        except Exception:
            logger.exception("ClinicalTrials.gov ingest refresh failed")

    def _refresh_civic(self) -> None:
        """Refresh CIViC (Clinical Interpretation of Variants in Cancer) evidence."""
        logger.info("Starting CIViC ingest refresh")
        try:
            from .ingest.civic_parser import CIViCIngestPipeline

            pipeline = CIViCIngestPipeline(
                collection_manager=self.collection_manager,
                embedder=self.embedder,
            )
            pipeline.run()

            now = time.time()
            self._last_run["civic"] = now
            LAST_INGEST.labels(source="civic").set(now)
            logger.info("CIViC ingest refresh completed")

        except Exception:
            logger.exception("CIViC ingest refresh failed")


# -----------------------------------------------------------------------
# No-op stub when apscheduler is not installed
# -----------------------------------------------------------------------

class _NoOpIngestScheduler:
    """Drop-in replacement that silently ignores all scheduler operations."""

    def __init__(self, *args, **kwargs) -> None:
        logger.info("Using no-op IngestScheduler (apscheduler not installed)")

    def start(self) -> None:
        pass

    def stop(self) -> None:
        pass

    def get_status(self) -> Dict[str, Any]:
        return {
            "next_run_time": None,
            "last_run_time": None,
            "job_count": 0,
        }


# Expose the appropriate class depending on availability
if not _APSCHEDULER_AVAILABLE:
    IngestScheduler = _NoOpIngestScheduler  # type: ignore[misc, assignment]
