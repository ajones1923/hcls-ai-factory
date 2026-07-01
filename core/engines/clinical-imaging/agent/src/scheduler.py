"""APScheduler-based periodic ingest for Imaging Intelligence Agent.

Schedules PubMed literature ingest weekly (168 hours) and
ClinicalTrials.gov trials ingest monthly (720 hours).  Both intervals
are configurable via ImagingSettings.

Author: Adam Jones
Date: February 2026
"""

import time
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from loguru import logger

try:
    from apscheduler.schedulers.background import BackgroundScheduler
    from apscheduler.triggers.interval import IntervalTrigger
    APSCHEDULER_AVAILABLE = True
except ImportError:
    APSCHEDULER_AVAILABLE = False
    logger.warning("apscheduler not installed — scheduler will be unavailable")

# Default intervals (hours)
PUBMED_INTERVAL_HOURS = 168   # 7 days
TRIALS_INTERVAL_HOURS = 720   # 30 days


def _run_pubmed_ingest(collection_manager, embedder) -> int:
    """Run the PubMed imaging literature ingest pipeline.

    Returns:
        Number of records ingested.
    """
    from src.ingest.literature_parser import PubMedImagingIngestPipeline

    pipeline = PubMedImagingIngestPipeline(collection_manager, embedder)
    return pipeline.run()


def _run_trials_ingest(collection_manager, embedder) -> int:
    """Run the ClinicalTrials.gov imaging trials ingest pipeline.

    Returns:
        Number of records ingested.
    """
    from src.ingest.clinical_trials_parser import ImagingTrialsIngestPipeline

    pipeline = ImagingTrialsIngestPipeline(collection_manager, embedder)
    return pipeline.run()


class ImagingIngestScheduler:
    """Manages periodic ingestion of imaging literature and trial data.

    Uses APScheduler BackgroundScheduler to run two independent jobs:
      - PubMed literature ingest (default: every 168 hours / 7 days)
      - ClinicalTrials.gov trials ingest (default: every 720 hours / 30 days)

    Usage:
        scheduler = ImagingIngestScheduler(
            collection_manager=milvus_manager,
            embedder=embedder,
        )
        scheduler.start(
            pubmed_interval_hours=168,
            trials_interval_hours=720,
        )
        # ... application runs ...
        scheduler.stop()
    """

    def __init__(
        self,
        pubmed_ingestor=None,
        trials_ingestor=None,
        collection_manager=None,
        embedder=None,
    ):
        """Initialize the scheduler.

        Args:
            pubmed_ingestor: Callable that ingests PubMed records. Signature:
                fn(collection_manager, embedder) -> int (records ingested).
                If None, the built-in PubMedImagingIngestPipeline is used.
            trials_ingestor: Callable that ingests ClinicalTrials.gov records.
                Signature: fn(collection_manager, embedder) -> int.
                If None, the built-in ImagingTrialsIngestPipeline is used.
            collection_manager: Milvus collection manager instance.
            embedder: Sentence-transformer embedder instance.
        """
        self.pubmed_ingestor = pubmed_ingestor or _run_pubmed_ingest
        self.trials_ingestor = trials_ingestor or _run_trials_ingest
        self.collection_manager = collection_manager
        self.embedder = embedder
        self._scheduler: Optional[Any] = None
        self._running = False
        self._last_run_stats: Dict[str, Any] = {}
        self._pubmed_interval_hours: int = PUBMED_INTERVAL_HOURS
        self._trials_interval_hours: int = TRIALS_INTERVAL_HOURS

    @property
    def is_running(self) -> bool:
        return self._running

    @property
    def last_run_stats(self) -> Dict[str, Any]:
        return self._last_run_stats

    @property
    def job_ids(self) -> List[str]:
        """Return list of registered job IDs."""
        if self._scheduler is None:
            return []
        return [job.id for job in self._scheduler.get_jobs()]

    def get_schedule_info(self) -> Dict[str, Any]:
        """Return human-readable schedule information."""
        info: Dict[str, Any] = {
            "running": self._running,
            "pubmed_interval_hours": self._pubmed_interval_hours,
            "trials_interval_hours": self._trials_interval_hours,
            "jobs": [],
        }
        if self._scheduler:
            for job in self._scheduler.get_jobs():
                info["jobs"].append({
                    "id": job.id,
                    "name": job.name,
                    "next_run": str(job.next_run_time) if job.next_run_time else None,
                })
        return info

    def start(
        self,
        pubmed_interval_hours: int = PUBMED_INTERVAL_HOURS,
        trials_interval_hours: int = TRIALS_INTERVAL_HOURS,
    ) -> None:
        """Start the periodic ingest scheduler with separate PubMed and trials jobs.

        Args:
            pubmed_interval_hours: Hours between PubMed ingest runs.
                Default 168 (weekly / 7 days).
            trials_interval_hours: Hours between ClinicalTrials.gov ingest runs.
                Default 720 (monthly / 30 days).
        """
        if not APSCHEDULER_AVAILABLE:
            logger.error("Cannot start scheduler: apscheduler is not installed")
            return

        if self._running:
            logger.warning("Scheduler is already running")
            return

        self._pubmed_interval_hours = pubmed_interval_hours
        self._trials_interval_hours = trials_interval_hours

        self._scheduler = BackgroundScheduler(daemon=True)

        # Job 1: PubMed literature ingest (weekly)
        self._scheduler.add_job(
            self._run_pubmed_job,
            trigger=IntervalTrigger(hours=pubmed_interval_hours),
            id="pubmed_ingest",
            name="PubMed Imaging Literature Ingest",
            replace_existing=True,
        )

        # Job 2: ClinicalTrials.gov trials ingest (monthly)
        self._scheduler.add_job(
            self._run_trials_job,
            trigger=IntervalTrigger(hours=trials_interval_hours),
            id="trials_ingest",
            name="ClinicalTrials.gov Imaging Trials Ingest",
            replace_existing=True,
        )

        self._scheduler.start()
        self._running = True
        logger.info(
            f"Ingest scheduler started — "
            f"PubMed every {pubmed_interval_hours}h, "
            f"Trials every {trials_interval_hours}h"
        )

    def stop(self) -> None:
        """Stop the periodic ingest scheduler."""
        if self._scheduler and self._running:
            self._scheduler.shutdown(wait=False)
            self._running = False
            logger.info("Ingest scheduler stopped")
        else:
            logger.warning("Scheduler is not running")

    def _run_pubmed_job(self) -> int:
        """Scheduled job: run PubMed ingest and update stats."""
        logger.info("Scheduled PubMed ingest starting...")
        start = time.time()
        count = 0
        try:
            count = self.pubmed_ingestor(
                self.collection_manager,
                self.embedder,
            )
            logger.info(f"Scheduled PubMed ingest complete: {count} records in {time.time() - start:.1f}s")
        except Exception as e:
            logger.error(f"Scheduled PubMed ingest failed: {e}")
            self._last_run_stats["pubmed_error"] = str(e)
        self._last_run_stats["pubmed_records"] = count
        self._last_run_stats["pubmed_last_run"] = datetime.now(timezone.utc).isoformat()
        return count

    def _run_trials_job(self) -> int:
        """Scheduled job: run ClinicalTrials.gov ingest and update stats."""
        logger.info("Scheduled ClinicalTrials.gov ingest starting...")
        start = time.time()
        count = 0
        try:
            count = self.trials_ingestor(
                self.collection_manager,
                self.embedder,
            )
            logger.info(f"Scheduled trials ingest complete: {count} records in {time.time() - start:.1f}s")
        except Exception as e:
            logger.error(f"Scheduled ClinicalTrials.gov ingest failed: {e}")
            self._last_run_stats["trials_error"] = str(e)
        self._last_run_stats["trials_records"] = count
        self._last_run_stats["trials_last_run"] = datetime.now(timezone.utc).isoformat()
        return count

    def run_ingest_cycle(self) -> Dict[str, Any]:
        """Run a single ingest cycle: PubMed + ClinicalTrials.gov.

        This method can be called manually or is invoked by the scheduler.

        Returns:
            Dict with ingest statistics including record counts and timing.
        """
        logger.info("Starting ingest cycle...")
        stats: Dict[str, Any] = {
            "start_time": time.time(),
            "pubmed_records": 0,
            "trials_records": 0,
            "errors": [],
        }

        # PubMed ingest
        if self.pubmed_ingestor:
            try:
                count = self.pubmed_ingestor(
                    self.collection_manager,
                    self.embedder,
                )
                stats["pubmed_records"] = count
                logger.info(f"PubMed ingest complete: {count} records")
            except Exception as e:
                stats["errors"].append(f"PubMed: {e}")
                logger.error(f"PubMed ingest failed: {e}")
        else:
            logger.warning("No PubMed ingestor configured — skipping")

        # ClinicalTrials.gov ingest
        if self.trials_ingestor:
            try:
                count = self.trials_ingestor(
                    self.collection_manager,
                    self.embedder,
                )
                stats["trials_records"] = count
                logger.info(f"ClinicalTrials.gov ingest complete: {count} records")
            except Exception as e:
                stats["errors"].append(f"Trials: {e}")
                logger.error(f"ClinicalTrials.gov ingest failed: {e}")
        else:
            logger.warning("No trials ingestor configured — skipping")

        stats["end_time"] = time.time()
        stats["duration_seconds"] = stats["end_time"] - stats["start_time"]
        stats["total_records"] = stats["pubmed_records"] + stats["trials_records"]
        stats["success"] = len(stats["errors"]) == 0

        self._last_run_stats = stats
        logger.info(
            f"Ingest cycle complete: {stats['total_records']} records "
            f"in {stats['duration_seconds']:.1f}s "
            f"({'OK' if stats['success'] else 'ERRORS: ' + '; '.join(stats['errors'])})"
        )
        return stats


# =====================================================================
# Module-level convenience function for API startup
# =====================================================================

_global_scheduler: Optional[ImagingIngestScheduler] = None


def start_scheduler() -> Optional[ImagingIngestScheduler]:
    """Start the global ingest scheduler using application settings.

    Called from the FastAPI lifespan when INGEST_ENABLED is True.
    Creates the scheduler with the configured collection manager and
    embedder, then starts both PubMed and trials jobs.

    Returns:
        The running ImagingIngestScheduler instance, or None if
        INGEST_ENABLED is False or APScheduler is unavailable.
    """
    global _global_scheduler

    from config.settings import settings

    if not settings.INGEST_ENABLED:
        logger.info("Ingest scheduler disabled (INGEST_ENABLED=False)")
        return None

    if not APSCHEDULER_AVAILABLE:
        logger.error("Cannot start scheduler: apscheduler is not installed")
        return None

    if _global_scheduler and _global_scheduler.is_running:
        logger.warning("Global scheduler already running")
        return _global_scheduler

    # Lazy-import heavy dependencies so this module stays light
    try:
        from sentence_transformers import SentenceTransformer
        from src.collections import ImagingCollectionManager

        manager = ImagingCollectionManager(
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        manager.connect()
        manager.ensure_collections()

        embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
    except Exception as e:
        logger.error(f"Failed to initialize scheduler dependencies: {e}")
        return None

    _global_scheduler = ImagingIngestScheduler(
        collection_manager=manager,
        embedder=embedder,
    )
    _global_scheduler.start(
        pubmed_interval_hours=settings.INGEST_SCHEDULE_HOURS,
        trials_interval_hours=settings.INGEST_TRIALS_SCHEDULE_HOURS,
    )

    schedule_info = _global_scheduler.get_schedule_info()
    logger.info(f"Ingest scheduler schedule: {schedule_info}")
    return _global_scheduler


def stop_scheduler() -> None:
    """Stop the global ingest scheduler if running."""
    global _global_scheduler
    if _global_scheduler and _global_scheduler.is_running:
        _global_scheduler.stop()
        _global_scheduler = None
        logger.info("Global ingest scheduler stopped")


def get_scheduler() -> Optional[ImagingIngestScheduler]:
    """Return the global scheduler instance (may be None)."""
    return _global_scheduler
