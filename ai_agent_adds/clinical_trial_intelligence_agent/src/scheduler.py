"""Automated ingest scheduler for the Clinical Trial Intelligence Agent.

Periodically refreshes ClinicalTrials.gov studies, PubMed clinical trial
publications, and FDA regulatory updates so the knowledge base stays
current without manual intervention.

Uses APScheduler's BackgroundScheduler so jobs run in a daemon thread
alongside the FastAPI / Streamlit application.

Default cadence:
  - ClinicalTrials.gov:    every INGEST_SCHEDULE_HOURS (default 24h / daily)
  - PubMed literature:     every INGEST_SCHEDULE_HOURS (default 24h / daily)
  - FDA regulatory:        every INGEST_SCHEDULE_HOURS * 7 (default 168h / weekly)

If ``apscheduler`` is not installed the module exports a no-op
``TrialScheduler`` stub so dependent code can import unconditionally.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import collections
import logging
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

# Import metrics (always available -- stubs if prometheus_client missing)
from .metrics import (
    INGEST_ERRORS,
    MetricsCollector,
)

logger = logging.getLogger(__name__)

try:
    from apscheduler.schedulers.background import BackgroundScheduler

    _APSCHEDULER_AVAILABLE = True
except ImportError:
    _APSCHEDULER_AVAILABLE = False


# ===================================================================
# DEFAULT SETTINGS DATACLASS
# ===================================================================


@dataclass
class TrialSchedulerSettings:
    """Configuration for the clinical trial ingest scheduler.

    Attributes:
        INGEST_ENABLED: Master switch for scheduled ingest jobs.
        INGEST_SCHEDULE_HOURS: Base interval in hours for ClinicalTrials.gov
            and PubMed refreshes. FDA regulatory runs at 7x this interval.
        TRIALS_CONDITIONS: ClinicalTrials.gov condition search terms.
        PUBMED_QUERY: PubMed search query for clinical trial publications.
        MAX_TRIALS_RESULTS: Maximum trials per refresh cycle.
        MAX_PUBMED_RESULTS: Maximum PubMed articles per refresh cycle.
        REGULATORY_DRUGS: Drug names for FDA regulatory monitoring.
    """

    INGEST_ENABLED: bool = True
    INGEST_SCHEDULE_HOURS: int = 24  # daily
    TRIALS_CONDITIONS: List[str] = field(
        default_factory=lambda: [
            "cancer",
            "heart failure",
            "alzheimer",
            "diabetes",
            "rare disease",
            "autoimmune",
            "infectious disease",
        ]
    )
    PUBMED_QUERY: str = (
        '"clinical trial"[Publication Type] AND '
        '("2024"[Date - Publication] : "3000"[Date - Publication])'
    )
    MAX_TRIALS_RESULTS: int = 200
    MAX_PUBMED_RESULTS: int = 500
    REGULATORY_DRUGS: List[str] = field(
        default_factory=lambda: [
            "pembrolizumab",
            "nivolumab",
            "osimertinib",
            "dapagliflozin",
            "empagliflozin",
        ]
    )


# ===================================================================
# INGEST JOB STATUS
# ===================================================================


@dataclass
class IngestJobStatus:
    """Status of a single ingest job execution."""

    job_id: str
    source: str
    status: str = "pending"  # pending | running | success | error
    records_ingested: int = 0
    error_message: Optional[str] = None
    started_at: Optional[str] = None
    completed_at: Optional[str] = None
    duration_seconds: float = 0.0


# ===================================================================
# SCHEDULER IMPLEMENTATION
# ===================================================================


if _APSCHEDULER_AVAILABLE:

    class TrialScheduler:
        """Background scheduler for periodic clinical trial data ingestion.

        Manages three recurring jobs:
          1. ClinicalTrials.gov trial refresh (daily)
          2. PubMed clinical trial literature refresh (daily)
          3. FDA regulatory update check (weekly)

        Usage::

            from src.scheduler import TrialScheduler, TrialSchedulerSettings

            settings = TrialSchedulerSettings(INGEST_ENABLED=True)
            scheduler = TrialScheduler(
                settings=settings,
                collection_manager=cm,
                embedder=embedder,
            )
            scheduler.start()
            ...
            scheduler.stop()
        """

        def __init__(
            self,
            settings: Optional[TrialSchedulerSettings] = None,
            collection_manager: Any = None,
            embedder: Any = None,
        ):
            """Initialize the clinical trial ingest scheduler.

            Args:
                settings: Scheduler configuration. Uses defaults if ``None``.
                collection_manager: Milvus collection manager instance.
                embedder: Embedding model with an ``encode()`` method.
            """
            self.settings = settings or TrialSchedulerSettings()
            self.collection_manager = collection_manager
            self.embedder = embedder
            self.scheduler = BackgroundScheduler(daemon=True)
            self.logger = logging.getLogger(__name__)
            self._job_history: collections.deque = collections.deque(maxlen=100)
            self._last_run_time: Optional[float] = None

        # ── Public API ────────────────────────────────────────────────

        def start(self) -> None:
            """Start the scheduler with configured jobs.

            If ``INGEST_ENABLED`` is ``False`` in settings, logs a message
            and returns without starting any jobs.
            """
            if not self.settings or not self.settings.INGEST_ENABLED:
                self.logger.info("Scheduled ingest disabled.")
                return

            hours = self.settings.INGEST_SCHEDULE_HOURS

            self.scheduler.add_job(
                self._run_trials_ingest,
                "interval",
                hours=hours,
                id="trials_ingest",
                name="ClinicalTrials.gov trial refresh",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_pubmed_ingest,
                "interval",
                hours=hours,
                id="pubmed_ingest",
                name="PubMed clinical trial literature refresh",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_regulatory_check,
                "interval",
                hours=hours * 7,  # Weekly (7x daily)
                id="regulatory_check",
                name="FDA regulatory update check",
                replace_existing=True,
            )

            self.scheduler.start()
            self.logger.info(
                f"TrialScheduler started -- "
                f"Trials/PubMed every {hours}h, "
                f"Regulatory every {hours * 7}h ({hours * 7 // 24}d)"
            )

        def stop(self) -> None:
            """Gracefully shut down the background scheduler."""
            if self.scheduler.running:
                self.scheduler.shutdown(wait=False)
                self.logger.info("TrialScheduler stopped")

        def get_jobs(self) -> list:
            """Return a list of scheduled job summaries.

            Returns:
                List of dicts with ``id``, ``name``, and ``next_run_time``
                for each registered job.
            """
            jobs = self.scheduler.get_jobs()
            return [
                {
                    "id": job.id,
                    "name": job.name,
                    "next_run_time": (
                        job.next_run_time.isoformat()
                        if job.next_run_time
                        else None
                    ),
                }
                for job in jobs
            ]

        def get_status(self) -> Dict[str, Any]:
            """Return a comprehensive status summary.

            Returns:
                Dict with scheduler state, job list, last run time,
                and recent job history.
            """
            jobs = self.get_jobs()
            next_times = [
                j["next_run_time"] for j in jobs if j["next_run_time"]
            ]

            return {
                "running": self.scheduler.running,
                "ingest_enabled": self.settings.INGEST_ENABLED,
                "schedule_hours": self.settings.INGEST_SCHEDULE_HOURS,
                "next_run_time": next_times[0] if next_times else None,
                "last_run_time": self._last_run_time,
                "job_count": len(jobs),
                "jobs": jobs,
                "recent_history": [
                    {
                        "job_id": h.job_id,
                        "source": h.source,
                        "status": h.status,
                        "records": h.records_ingested,
                        "duration_s": round(h.duration_seconds, 1),
                        "completed_at": h.completed_at,
                    }
                    for h in self._job_history[-10:]
                ],
            }

        def trigger_manual_ingest(self, source: str) -> dict:
            """Trigger an immediate manual ingest for the specified source.

            Args:
                source: One of ``"trials"``, ``"pubmed"``, or ``"regulatory"``.

            Returns:
                Dict with ``status`` and ``message`` keys.
            """
            dispatch = {
                "trials": self._run_trials_ingest,
                "pubmed": self._run_pubmed_ingest,
                "regulatory": self._run_regulatory_check,
            }

            runner = dispatch.get(source.lower())
            if runner is None:
                return {
                    "status": "error",
                    "message": (
                        f"Unknown source '{source}'. "
                        f"Valid sources: {', '.join(dispatch.keys())}"
                    ),
                }

            self.logger.info(f"Manual ingest triggered for source: {source}")
            try:
                runner()
                return {
                    "status": "success",
                    "message": f"Manual ingest for '{source}' completed.",
                }
            except Exception as exc:
                return {
                    "status": "error",
                    "message": f"Manual ingest for '{source}' failed: {exc}",
                }

        # ── Private Job Wrappers ──────────────────────────────────────

        def _run_trials_ingest(self) -> None:
            """Run the ClinicalTrials.gov trial ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"trials_{int(time.time())}",
                source="clinicaltrials",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting ClinicalTrials.gov refresh")
            start = time.time()

            try:
                from .ingest.clinicaltrials_parser import ClinicalTrialsParser

                parser = ClinicalTrialsParser(
                    collection_manager=self.collection_manager,
                    embedder=self.embedder,
                )
                records, stats = parser.run(
                    conditions=self.settings.TRIALS_CONDITIONS,
                    max_results=self.settings.MAX_TRIALS_RESULTS,
                )
                elapsed = time.time() - start
                self._last_run_time = time.time()
                count = len(records)

                MetricsCollector.record_ingest(
                    source="clinicaltrials",
                    duration=elapsed,
                    record_count=count,
                    collection="trial_protocols",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    f"Scheduler: ClinicalTrials.gov refresh complete -- "
                    f"{count} records in {elapsed:.1f}s"
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "ClinicalTrialsParser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Trials ingest skipped -- "
                    "clinicaltrials_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="clinicaltrials").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    f"Scheduler: ClinicalTrials.gov refresh failed -- {exc}"
                )

            self._job_history.append(job_status)

        def _run_pubmed_ingest(self) -> None:
            """Run the PubMed clinical trial literature ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"pubmed_{int(time.time())}",
                source="pubmed",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting PubMed trial literature refresh")
            start = time.time()

            try:
                from .ingest.pubmed_parser import PubMedTrialParser

                parser = PubMedTrialParser(
                    collection_manager=self.collection_manager,
                    embedder=self.embedder,
                )
                records, stats = parser.run(
                    query=self.settings.PUBMED_QUERY,
                    max_results=self.settings.MAX_PUBMED_RESULTS,
                )
                elapsed = time.time() - start
                self._last_run_time = time.time()
                count = len(records)

                MetricsCollector.record_ingest(
                    source="pubmed",
                    duration=elapsed,
                    record_count=count,
                    collection="trial_literature",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    f"Scheduler: PubMed refresh complete -- "
                    f"{count} records in {elapsed:.1f}s"
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "PubMedTrialParser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: PubMed ingest skipped -- "
                    "pubmed_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="pubmed").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    f"Scheduler: PubMed refresh failed -- {exc}"
                )

            self._job_history.append(job_status)

        def _run_regulatory_check(self) -> None:
            """Run the FDA regulatory update check."""
            job_status = IngestJobStatus(
                job_id=f"regulatory_{int(time.time())}",
                source="regulatory",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting FDA regulatory update check")
            start = time.time()

            try:
                from .ingest.regulatory_parser import RegulatoryParser

                parser = RegulatoryParser(
                    collection_manager=self.collection_manager,
                    embedder=self.embedder,
                )
                records, stats = parser.run(
                    drug_names=self.settings.REGULATORY_DRUGS,
                    include_milestones=True,
                )
                elapsed = time.time() - start
                self._last_run_time = time.time()
                count = len(records)

                MetricsCollector.record_ingest(
                    source="regulatory",
                    duration=elapsed,
                    record_count=count,
                    collection="trial_regulatory",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    f"Scheduler: Regulatory check complete -- "
                    f"{count} records in {elapsed:.1f}s"
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "RegulatoryParser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Regulatory check skipped -- "
                    "regulatory_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="regulatory").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    f"Scheduler: Regulatory check failed -- {exc}"
                )

            self._job_history.append(job_status)

else:
    # ── No-op stub when apscheduler is not installed ──────────────────

    class TrialScheduler:  # type: ignore[no-redef]
        """No-op scheduler stub (apscheduler not installed).

        All methods are safe to call but perform no work. Install
        apscheduler to enable scheduled ingest::

            pip install apscheduler>=3.10.0
        """

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            logger.warning(
                "apscheduler is not installed -- TrialScheduler is a no-op. "
                "Install with: pip install apscheduler>=3.10.0"
            )

        def start(self) -> None:
            pass

        def stop(self) -> None:
            pass

        def get_jobs(self) -> list:
            return []

        def get_status(self) -> Dict[str, Any]:
            return {
                "running": False,
                "ingest_enabled": False,
                "schedule_hours": 0,
                "next_run_time": None,
                "last_run_time": None,
                "job_count": 0,
                "jobs": [],
                "recent_history": [],
            }

        def trigger_manual_ingest(self, source: str) -> dict:
            return {
                "status": "error",
                "message": (
                    "Scheduler unavailable -- apscheduler is not installed."
                ),
            }
