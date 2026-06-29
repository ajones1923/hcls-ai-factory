"""Automated ingest scheduler for the Cardiology Intelligence Agent.

Periodically refreshes PubMed cardiovascular literature, ClinicalTrials.gov
cardiovascular trials, and ACC/AHA clinical practice guidelines so the
knowledge base stays current without manual intervention.

Uses APScheduler's BackgroundScheduler so jobs run in a daemon thread
alongside the FastAPI / Streamlit application.

Default cadence:
  - PubMed literature:        every INGEST_SCHEDULE_HOURS (default 168h / weekly)
  - ClinicalTrials.gov:       every INGEST_SCHEDULE_HOURS (default 168h / weekly)
  - ACC/AHA guideline check:  every INGEST_SCHEDULE_HOURS * 4 (default 672h / monthly)

If ``apscheduler`` is not installed the module exports a no-op
``CardioScheduler`` stub so dependent code can import unconditionally.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

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


# ═══════════════════════════════════════════════════════════════════════
# DEFAULT SETTINGS DATACLASS
# ═══════════════════════════════════════════════════════════════════════


@dataclass
class CardioSchedulerSettings:
    """Configuration for the cardiology ingest scheduler.

    Attributes:
        INGEST_ENABLED: Master switch for scheduled ingest jobs.
        INGEST_SCHEDULE_HOURS: Base interval in hours for PubMed and
            ClinicalTrials.gov refreshes.  The guideline check runs at
            4x this interval.
        PUBMED_QUERY: PubMed search query for cardiovascular literature.
        TRIALS_CONDITIONS: ClinicalTrials.gov condition search terms.
        MAX_PUBMED_RESULTS: Maximum PubMed articles per refresh cycle.
        MAX_TRIALS_RESULTS: Maximum clinical trials per refresh cycle.
        GUIDELINE_SOURCES: List of guideline sources to monitor.
    """

    INGEST_ENABLED: bool = True
    INGEST_SCHEDULE_HOURS: int = 168  # weekly
    PUBMED_QUERY: str = (
        "(cardiovascular diseases[MeSH] OR heart failure[MeSH] "
        "OR coronary artery disease[MeSH] OR atrial fibrillation[MeSH] "
        "OR valvular heart disease[MeSH]) AND (therapy OR management "
        "OR guidelines OR risk)"
    )
    TRIALS_CONDITIONS: List[str] = field(
        default_factory=lambda: [
            "heart failure",
            "coronary artery disease",
            "atrial fibrillation",
            "valvular heart disease",
            "cardiomyopathy",
            "acute coronary syndrome",
            "cardio-oncology",
        ]
    )
    MAX_PUBMED_RESULTS: int = 500
    MAX_TRIALS_RESULTS: int = 200
    GUIDELINE_SOURCES: List[str] = field(
        default_factory=lambda: [
            "ACC/AHA",
            "ESC",
            "HRS",
        ]
    )


# ═══════════════════════════════════════════════════════════════════════
# INGEST JOB STATUS
# ═══════════════════════════════════════════════════════════════════════


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


# ═══════════════════════════════════════════════════════════════════════
# SCHEDULER IMPLEMENTATION
# ═══════════════════════════════════════════════════════════════════════


if _APSCHEDULER_AVAILABLE:

    class CardioScheduler:
        """Background scheduler for periodic cardiovascular data ingestion.

        Manages three recurring jobs:
          1. PubMed cardiovascular literature refresh
          2. ClinicalTrials.gov cardiovascular trials refresh
          3. ACC/AHA guideline update check

        Usage::

            from src.scheduler import CardioScheduler, CardioSchedulerSettings

            settings = CardioSchedulerSettings(INGEST_ENABLED=True)
            scheduler = CardioScheduler(
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
            settings: Optional[CardioSchedulerSettings] = None,
            collection_manager: Any = None,
            embedder: Any = None,
        ):
            """Initialize the cardiology ingest scheduler.

            Args:
                settings: Scheduler configuration. Uses defaults if ``None``.
                collection_manager: Milvus collection manager instance.
                embedder: Embedding model with an ``encode()`` method.
            """
            self.settings = settings or CardioSchedulerSettings()
            self.collection_manager = collection_manager
            self.embedder = embedder
            self.scheduler = BackgroundScheduler(daemon=True)
            self.logger = logging.getLogger(__name__)
            self._job_history: List[IngestJobStatus] = []
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
                self._run_pubmed_ingest,
                "interval",
                hours=hours,
                id="pubmed_ingest",
                name="PubMed cardiovascular literature ingest",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_trials_ingest,
                "interval",
                hours=hours,
                id="trials_ingest",
                name="ClinicalTrials.gov cardiovascular trials ingest",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_guideline_check,
                "interval",
                hours=hours * 4,  # Monthly (4x weekly)
                id="guideline_check",
                name="ACC/AHA guideline update check",
                replace_existing=True,
            )

            self.scheduler.start()
            self.logger.info(
                f"CardioScheduler started -- "
                f"PubMed/Trials every {hours}h ({hours // 24}d), "
                f"Guidelines every {hours * 4}h ({hours * 4 // 24}d)"
            )

        def stop(self) -> None:
            """Gracefully shut down the background scheduler."""
            if self.scheduler.running:
                self.scheduler.shutdown(wait=False)
                self.logger.info("CardioScheduler stopped")

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
                source: One of ``"pubmed"``, ``"trials"``, or ``"guidelines"``.

            Returns:
                Dict with ``status`` and ``message`` keys.
            """
            dispatch = {
                "pubmed": self._run_pubmed_ingest,
                "trials": self._run_trials_ingest,
                "guidelines": self._run_guideline_check,
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

        def _run_pubmed_ingest(self) -> None:
            """Run the PubMed cardiovascular literature ingest pipeline.

            Searches PubMed for recent cardiovascular publications, parses
            results, generates embeddings, and upserts into the
            ``cardio_literature`` collection.

            Updates Prometheus ``cardio_last_ingest_timestamp{source="pubmed"}``
            on success.
            """
            job_status = IngestJobStatus(
                job_id=f"pubmed_{int(time.time())}",
                source="pubmed",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting PubMed cardiovascular refresh")
            start = time.time()

            try:
                from .ingest.literature_parser import PubMedIngestPipeline

                pipeline = PubMedIngestPipeline(
                    self.collection_manager,
                    self.embedder,
                )
                count = pipeline.run(
                    query=self.settings.PUBMED_QUERY,
                    max_results=self.settings.MAX_PUBMED_RESULTS,
                )
                elapsed = time.time() - start
                self._last_run_time = time.time()

                # Update metrics
                MetricsCollector.record_ingest(
                    source="pubmed",
                    duration=elapsed,
                    record_count=count,
                    collection="cardio_literature",
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
                job_status.error_message = "PubMedIngestPipeline not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: PubMed ingest skipped -- "
                    "literature_parser module not available"
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

        def _run_trials_ingest(self) -> None:
            """Run the ClinicalTrials.gov cardiovascular trials ingest.

            Queries ClinicalTrials.gov for active cardiovascular trials,
            parses results, generates embeddings, and upserts into the
            ``cardio_trials`` collection.

            Updates Prometheus
            ``cardio_last_ingest_timestamp{source="clinical_trials"}``
            on success.
            """
            job_status = IngestJobStatus(
                job_id=f"trials_{int(time.time())}",
                source="clinical_trials",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info(
                "Scheduler: starting ClinicalTrials.gov cardiovascular refresh"
            )
            start = time.time()

            try:
                from .ingest.clinical_trials_parser import (
                    ClinicalTrialsIngestPipeline,
                )

                pipeline = ClinicalTrialsIngestPipeline(
                    self.collection_manager,
                    self.embedder,
                )
                count = pipeline.run(
                    conditions=self.settings.TRIALS_CONDITIONS,
                    max_results=self.settings.MAX_TRIALS_RESULTS,
                )
                elapsed = time.time() - start
                self._last_run_time = time.time()

                MetricsCollector.record_ingest(
                    source="clinical_trials",
                    duration=elapsed,
                    record_count=count,
                    collection="cardio_trials",
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
                job_status.error_message = (
                    "ClinicalTrialsIngestPipeline not available"
                )
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Trials ingest skipped -- "
                    "clinical_trials_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="clinical_trials").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    f"Scheduler: ClinicalTrials.gov refresh failed -- {exc}"
                )

            self._job_history.append(job_status)

        def _run_guideline_check(self) -> None:
            """Check for ACC/AHA/ESC guideline updates.

            Queries guideline sources for recently published or updated
            clinical practice guidelines in cardiovascular medicine.
            New or updated guidelines are ingested into the
            ``cardio_guidelines`` collection.

            Updates Prometheus
            ``cardio_last_ingest_timestamp{source="guidelines"}``
            on success.
            """
            job_status = IngestJobStatus(
                job_id=f"guidelines_{int(time.time())}",
                source="guidelines",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info(
                "Scheduler: starting guideline update check "
                f"(sources: {', '.join(self.settings.GUIDELINE_SOURCES)})"
            )
            start = time.time()

            try:
                from .ingest.guideline_parser import GuidelineIngestPipeline

                pipeline = GuidelineIngestPipeline(
                    self.collection_manager,
                    self.embedder,
                )
                count = pipeline.run(
                    sources=self.settings.GUIDELINE_SOURCES,
                )
                elapsed = time.time() - start
                self._last_run_time = time.time()

                MetricsCollector.record_ingest(
                    source="guidelines",
                    duration=elapsed,
                    record_count=count,
                    collection="cardio_guidelines",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    f"Scheduler: Guideline check complete -- "
                    f"{count} updates in {elapsed:.1f}s"
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = (
                    "GuidelineIngestPipeline not available"
                )
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Guideline check skipped -- "
                    "guideline_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="guidelines").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    f"Scheduler: Guideline check failed -- {exc}"
                )

            self._job_history.append(job_status)

else:
    # ── No-op stub when apscheduler is not installed ──────────────────

    class CardioScheduler:  # type: ignore[no-redef]
        """No-op scheduler stub (apscheduler not installed).

        All methods are safe to call but perform no work.  Install
        apscheduler to enable scheduled ingest::

            pip install apscheduler>=3.10.0
        """

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            logger.warning(
                "apscheduler is not installed -- CardioScheduler is a no-op. "
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
