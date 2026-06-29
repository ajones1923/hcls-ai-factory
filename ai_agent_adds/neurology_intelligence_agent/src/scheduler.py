"""Automated ingest scheduler for the Neurology Intelligence Agent.

Periodically refreshes PubMed neurology literature, ClinicalTrials.gov
neurology studies, and neurology clinical practice guideline updates so
the knowledge base stays current without manual intervention.

Uses APScheduler's BackgroundScheduler so jobs run in a daemon thread
alongside the FastAPI / Streamlit application.

Default cadence:
  - PubMed neurology literature:  weekly (168h)
  - ClinicalTrials.gov neurology: weekly (168h)
  - Guideline updates:            monthly (720h)

If ``apscheduler`` is not installed the module exports a no-op
``NeuroScheduler`` stub so dependent code can import unconditionally.

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
class NeuroSchedulerSettings:
    """Configuration for the neurology ingest scheduler.

    Attributes:
        INGEST_ENABLED: Master switch for scheduled ingest jobs.
        PUBMED_SCHEDULE_HOURS: Interval in hours for PubMed refresh (default: weekly).
        TRIALS_SCHEDULE_HOURS: Interval in hours for ClinicalTrials.gov refresh (default: weekly).
        GUIDELINE_SCHEDULE_HOURS: Interval in hours for guideline refresh (default: monthly).
        PUBMED_QUERY: PubMed search query for neurology publications.
        TRIALS_CONDITIONS: ClinicalTrials.gov condition search terms.
        GUIDELINE_SOURCES: Guideline source organizations.
        MAX_PUBMED_RESULTS: Maximum PubMed articles per refresh cycle.
        MAX_TRIALS_RESULTS: Maximum trials per refresh cycle.
    """

    INGEST_ENABLED: bool = True
    PUBMED_SCHEDULE_HOURS: int = 168  # weekly
    TRIALS_SCHEDULE_HOURS: int = 168  # weekly
    GUIDELINE_SCHEDULE_HOURS: int = 720  # monthly (~30 days)
    PUBMED_QUERY: str = (
        '("neurology"[MeSH Terms] OR "neuroscience"[MeSH Terms]) AND '
        '("2024"[Date - Publication] : "3000"[Date - Publication])'
    )
    TRIALS_CONDITIONS: List[str] = field(
        default_factory=lambda: [
            "stroke",
            "Alzheimer's disease",
            "epilepsy",
            "multiple sclerosis",
            "Parkinson's disease",
            "migraine",
            "glioblastoma",
            "ALS",
            "myasthenia gravis",
            "Huntington disease",
        ]
    )
    GUIDELINE_SOURCES: List[str] = field(
        default_factory=lambda: [
            "AAN",
            "AHA/ASA",
            "ILAE",
            "EFNS/EAN",
            "NICE",
            "IHS",
        ]
    )
    MAX_PUBMED_RESULTS: int = 500
    MAX_TRIALS_RESULTS: int = 200


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

    class NeuroScheduler:
        """Background scheduler for periodic neurology data ingestion.

        Manages three recurring jobs:
          1. PubMed neurology literature refresh (weekly)
          2. ClinicalTrials.gov neurology trial refresh (weekly)
          3. Guideline update check (monthly)

        Usage::

            from src.scheduler import NeuroScheduler, NeuroSchedulerSettings

            settings = NeuroSchedulerSettings(INGEST_ENABLED=True)
            scheduler = NeuroScheduler(
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
            settings: Optional[NeuroSchedulerSettings] = None,
            collection_manager: Any = None,
            embedder: Any = None,
        ):
            self.settings = settings or NeuroSchedulerSettings()
            self.collection_manager = collection_manager
            self.embedder = embedder
            self.scheduler = BackgroundScheduler(daemon=True)
            self.logger = logging.getLogger(__name__)
            self._job_history: collections.deque = collections.deque(maxlen=100)
            self._last_run_time: Optional[float] = None

        # -- Public API --

        def start(self) -> None:
            """Start the scheduler with configured jobs."""
            if not self.settings or not self.settings.INGEST_ENABLED:
                self.logger.info("Scheduled ingest disabled.")
                return

            self.scheduler.add_job(
                self._run_pubmed_ingest,
                "interval",
                hours=self.settings.PUBMED_SCHEDULE_HOURS,
                id="neuro_pubmed_ingest",
                name="PubMed neurology literature refresh",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_trials_ingest,
                "interval",
                hours=self.settings.TRIALS_SCHEDULE_HOURS,
                id="neuro_trials_ingest",
                name="ClinicalTrials.gov neurology refresh",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_guideline_check,
                "interval",
                hours=self.settings.GUIDELINE_SCHEDULE_HOURS,
                id="neuro_guideline_check",
                name="Neurology guideline update check",
                replace_existing=True,
            )

            self.scheduler.start()
            self.logger.info(
                "NeuroScheduler started -- "
                "PubMed every %dh, Trials every %dh, Guidelines every %dh",
                self.settings.PUBMED_SCHEDULE_HOURS,
                self.settings.TRIALS_SCHEDULE_HOURS,
                self.settings.GUIDELINE_SCHEDULE_HOURS,
            )

        def stop(self) -> None:
            """Gracefully shut down the background scheduler."""
            if self.scheduler.running:
                self.scheduler.shutdown(wait=False)
                self.logger.info("NeuroScheduler stopped")

        def get_jobs(self) -> list:
            """Return a list of scheduled job summaries."""
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
            """Return a comprehensive status summary."""
            jobs = self.get_jobs()
            next_times = [
                j["next_run_time"] for j in jobs if j["next_run_time"]
            ]

            return {
                "running": self.scheduler.running,
                "ingest_enabled": self.settings.INGEST_ENABLED,
                "pubmed_schedule_hours": self.settings.PUBMED_SCHEDULE_HOURS,
                "trials_schedule_hours": self.settings.TRIALS_SCHEDULE_HOURS,
                "guideline_schedule_hours": self.settings.GUIDELINE_SCHEDULE_HOURS,
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
            """Trigger an immediate manual ingest for the specified source."""
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

            self.logger.info("Manual ingest triggered for source: %s", source)
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

        # -- Private Job Wrappers --

        def _run_pubmed_ingest(self) -> None:
            """Run the PubMed neurology literature ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"pubmed_{int(time.time())}",
                source="pubmed",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting PubMed neurology refresh")
            start = time.time()

            try:
                from .ingest.pubmed_neuro_parser import PubMedNeuroParser

                parser = PubMedNeuroParser(
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
                    collection="neuro_literature",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    "Scheduler: PubMed refresh complete -- "
                    "%d records in %.1fs",
                    count, elapsed,
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "PubMedNeuroParser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: PubMed ingest skipped -- "
                    "pubmed_neuro_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="pubmed").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    "Scheduler: PubMed refresh failed -- %s", exc
                )

            self._job_history.append(job_status)

        def _run_trials_ingest(self) -> None:
            """Run the ClinicalTrials.gov neurology trial ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"trials_{int(time.time())}",
                source="clinicaltrials",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting ClinicalTrials.gov neurology refresh")
            start = time.time()

            try:
                # Trials parser would be a separate module; gracefully handle ImportError
                raise ImportError("ClinicalTrials neurology parser not yet implemented")

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "ClinicalTrials neurology parser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Trials ingest skipped -- parser not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="clinicaltrials").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    "Scheduler: ClinicalTrials.gov refresh failed -- %s", exc
                )

            self._job_history.append(job_status)

        def _run_guideline_check(self) -> None:
            """Run the neurology guideline update check."""
            job_status = IngestJobStatus(
                job_id=f"guidelines_{int(time.time())}",
                source="guidelines",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting neurology guideline update check")
            start = time.time()

            try:
                # Guideline parser would be a separate module; gracefully handle ImportError
                raise ImportError("Guideline parser not yet implemented")

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "Guideline parser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Guideline check skipped -- parser not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="guidelines").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    "Scheduler: Guideline check failed -- %s", exc
                )

            self._job_history.append(job_status)

else:
    # -- No-op stub when apscheduler is not installed --

    class NeuroScheduler:  # type: ignore[no-redef]
        """No-op scheduler stub (apscheduler not installed).

        All methods are safe to call but perform no work. Install
        apscheduler to enable scheduled ingest::

            pip install apscheduler>=3.10.0
        """

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            logger.warning(
                "apscheduler is not installed -- NeuroScheduler is a no-op. "
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
                "pubmed_schedule_hours": 0,
                "trials_schedule_hours": 0,
                "guideline_schedule_hours": 0,
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
