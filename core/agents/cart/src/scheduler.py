"""Automated ingest scheduler for the CAR-T Intelligence Agent.

Periodically refreshes the PubMed literature and ClinicalTrials.gov
collections so the knowledge base stays current without manual
intervention.  Uses APScheduler's BackgroundScheduler so the jobs run
in a daemon thread alongside the FastAPI / Streamlit application.

Default cadence: every 168 hours (once per week).

If ``apscheduler`` is not installed the module exports a no-op
``IngestScheduler`` stub so dependent code can import unconditionally.

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

import time
from typing import Any, Dict, Optional

from loguru import logger

# Import metrics (always available — stubs if prometheus_client missing)
from .metrics import LAST_INGEST

try:
    from apscheduler.schedulers.background import BackgroundScheduler

    _APSCHEDULER_AVAILABLE = True
except ImportError:
    _APSCHEDULER_AVAILABLE = False


if _APSCHEDULER_AVAILABLE:

    class IngestScheduler:
        """Background scheduler that keeps PubMed and ClinicalTrials.gov fresh.

        Wraps ``apscheduler.BackgroundScheduler`` with two recurring interval
        jobs — one for each upstream data source.  Each job creates a
        short-lived ingest pipeline, runs it, and updates the
        ``cart_last_ingest_timestamp`` Prometheus gauge.

        Usage::

            from src.scheduler import IngestScheduler

            scheduler = IngestScheduler(collection_manager, embedder)
            scheduler.start()      # non-blocking — runs in background thread
            ...
            scheduler.stop()       # graceful shutdown
        """

        def __init__(
            self,
            collection_manager: Any,
            embedder: Any,
            interval_hours: int = 168,
        ):
            """Initialize the ingest scheduler.

            Args:
                collection_manager: ``CARTCollectionManager`` instance already
                    connected to Milvus.
                embedder: Embedding model / client with an ``encode()`` method
                    (e.g. SentenceTransformer wrapping BGE-small-en-v1.5).
                interval_hours: How often (in hours) each ingest job should
                    run.  Defaults to **168** (once per week).
            """
            self.collection_manager = collection_manager
            self.embedder = embedder
            self.interval_hours = interval_hours
            self._scheduler = BackgroundScheduler(daemon=True)
            self._last_run_time: Optional[float] = None

        # ── Public API ────────────────────────────────────────────────

        def start(self) -> None:
            """Start the background scheduler with PubMed and trials jobs.

            Both jobs execute immediately on first run (``next_run_time="now"``
            is *not* set — the first run happens after one full interval by
            default).  To trigger an immediate refresh call the private
            ``_refresh_*`` methods directly.
            """
            self._scheduler.add_job(
                self._refresh_pubmed,
                trigger="interval",
                hours=self.interval_hours,
                id="refresh_pubmed",
                name="PubMed CAR-T Literature Refresh",
                replace_existing=True,
            )

            self._scheduler.add_job(
                self._refresh_clinical_trials,
                trigger="interval",
                hours=self.interval_hours,
                id="refresh_clinical_trials",
                name="ClinicalTrials.gov CAR-T Refresh",
                replace_existing=True,
            )

            self._scheduler.start()
            logger.info(
                f"IngestScheduler started — refreshing every "
                f"{self.interval_hours}h ({self.interval_hours // 24}d)"
            )

        def stop(self) -> None:
            """Gracefully shut down the background scheduler."""
            if self._scheduler.running:
                self._scheduler.shutdown(wait=False)
                logger.info("IngestScheduler stopped")

        def get_status(self) -> Dict[str, Any]:
            """Return a status summary of all scheduled jobs.

            Returns:
                Dict with ``next_run_time``, ``last_run_time``, and
                ``job_count`` keys.  Times are ISO-8601 strings or
                ``None`` if not yet available.
            """
            jobs = self._scheduler.get_jobs()

            next_run_times = [
                j.next_run_time.isoformat()
                for j in jobs
                if j.next_run_time is not None
            ]

            return {
                "next_run_time": next_run_times[0] if next_run_times else None,
                "last_run_time": self._last_run_time,
                "job_count": len(jobs),
            }

        # ── Private job wrappers ──────────────────────────────────────

        def _refresh_pubmed(self) -> None:
            """Run the PubMed ingest pipeline with default CAR-T query.

            Updates the ``cart_last_ingest_timestamp{source="pubmed"}``
            Prometheus gauge on success.
            """
            from .ingest.literature_parser import PubMedIngestPipeline

            logger.info("Scheduler: starting PubMed refresh")
            start = time.time()

            try:
                pipeline = PubMedIngestPipeline(
                    self.collection_manager,
                    self.embedder,
                )
                count = pipeline.run()
                elapsed = time.time() - start
                self._last_run_time = time.time()

                LAST_INGEST.labels(source="pubmed").set(time.time())
                logger.info(
                    f"Scheduler: PubMed refresh complete — "
                    f"{count} records in {elapsed:.1f}s"
                )

            except Exception as exc:
                logger.error(f"Scheduler: PubMed refresh failed — {exc}")

        def _refresh_clinical_trials(self) -> None:
            """Run the ClinicalTrials.gov ingest pipeline.

            Updates the ``cart_last_ingest_timestamp{source="clinical_trials"}``
            Prometheus gauge on success.
            """
            from .ingest.clinical_trials_parser import (
                ClinicalTrialsIngestPipeline,
            )

            logger.info("Scheduler: starting ClinicalTrials.gov refresh")
            start = time.time()

            try:
                pipeline = ClinicalTrialsIngestPipeline(
                    self.collection_manager,
                    self.embedder,
                )
                count = pipeline.run()
                elapsed = time.time() - start
                self._last_run_time = time.time()

                LAST_INGEST.labels(source="clinical_trials").set(time.time())
                logger.info(
                    f"Scheduler: ClinicalTrials.gov refresh complete — "
                    f"{count} records in {elapsed:.1f}s"
                )

            except Exception as exc:
                logger.error(
                    f"Scheduler: ClinicalTrials.gov refresh failed — {exc}"
                )

else:
    # ── No-op stub when apscheduler is not installed ──────────────────

    class IngestScheduler:  # type: ignore[no-redef]
        """No-op scheduler stub (apscheduler not installed)."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            logger.warning(
                "apscheduler is not installed — IngestScheduler is a no-op. "
                "Install with: pip install apscheduler>=3.10.0"
            )

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
