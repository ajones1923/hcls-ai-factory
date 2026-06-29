"""APScheduler-based periodic ingest for Imaging Intelligence Agent.

Schedules weekly (or custom interval) ingestion cycles that pull new
records from PubMed and ClinicalTrials.gov into the Milvus collections.
"""

import time
from typing import Any, Dict, Optional

from loguru import logger

try:
    from apscheduler.schedulers.background import BackgroundScheduler
    from apscheduler.triggers.interval import IntervalTrigger
    APSCHEDULER_AVAILABLE = True
except ImportError:
    APSCHEDULER_AVAILABLE = False
    logger.warning("apscheduler not installed — scheduler will be unavailable")


class ImagingIngestScheduler:
    """Manages periodic ingestion of imaging literature and trial data.

    Uses APScheduler BackgroundScheduler to run ingest cycles at a
    configurable interval (default: weekly / 168 hours).

    Usage:
        scheduler = ImagingIngestScheduler(
            pubmed_ingestor=pubmed_ingest_fn,
            trials_ingestor=trials_ingest_fn,
            collection_manager=milvus_manager,
            embedder=embedder,
        )
        scheduler.start(interval_hours=168)
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
            trials_ingestor: Callable that ingests ClinicalTrials.gov records.
                Signature: fn(collection_manager, embedder) -> int.
            collection_manager: Milvus collection manager instance.
            embedder: Sentence-transformer embedder instance.
        """
        self.pubmed_ingestor = pubmed_ingestor
        self.trials_ingestor = trials_ingestor
        self.collection_manager = collection_manager
        self.embedder = embedder
        self._scheduler: Optional[Any] = None
        self._running = False
        self._last_run_stats: Dict[str, Any] = {}

    @property
    def is_running(self) -> bool:
        return self._running

    @property
    def last_run_stats(self) -> Dict[str, Any]:
        return self._last_run_stats

    def start(self, interval_hours: int = 168) -> None:
        """Start the periodic ingest scheduler.

        Args:
            interval_hours: Hours between ingest cycles. Default 168 (weekly).
        """
        if not APSCHEDULER_AVAILABLE:
            logger.error("Cannot start scheduler: apscheduler is not installed")
            return

        if self._running:
            logger.warning("Scheduler is already running")
            return

        self._scheduler = BackgroundScheduler(daemon=True)
        self._scheduler.add_job(
            self.run_ingest_cycle,
            trigger=IntervalTrigger(hours=interval_hours),
            id="imaging_ingest_cycle",
            name="Imaging Intelligence Ingest Cycle",
            replace_existing=True,
        )
        self._scheduler.start()
        self._running = True
        logger.info(f"Ingest scheduler started — interval={interval_hours}h")

    def stop(self) -> None:
        """Stop the periodic ingest scheduler."""
        if self._scheduler and self._running:
            self._scheduler.shutdown(wait=False)
            self._running = False
            logger.info("Ingest scheduler stopped")
        else:
            logger.warning("Scheduler is not running")

    def run_ingest_cycle(self) -> Dict[str, Any]:
        """Run a single ingest cycle: PubMed + ClinicalTrials.gov.

        This method can be called manually or is invoked by the scheduler.

        Returns:
            Dict with ingest statistics including record counts and timing.
        """
        logger.info("Starting ingest cycle...")
        stats = {
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
