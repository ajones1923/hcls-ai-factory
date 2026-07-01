"""Automated ingest scheduler for the Rare Disease Diagnostic Agent.

Periodically refreshes disease databases (ClinVar, OMIM, HPO, Orphanet)
so the knowledge base stays current without manual intervention.

Uses APScheduler's BackgroundScheduler so jobs run in a daemon thread
alongside the FastAPI / Streamlit application.

Default cadence:
  - ClinVar variants:    weekly  (INGEST_SCHEDULE_HOURS * 7)
  - OMIM diseases:       weekly  (INGEST_SCHEDULE_HOURS * 7)
  - HPO phenotypes:      monthly (INGEST_SCHEDULE_HOURS * 30)
  - Orphanet diseases:   monthly (INGEST_SCHEDULE_HOURS * 30)

If ``apscheduler`` is not installed the module exports a no-op
``RareDiseaseScheduler`` stub so dependent code can import unconditionally.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import collections
import logging
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Optional

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
class RareDiseaseSchedulerSettings:
    """Configuration for the rare disease ingest scheduler.

    Attributes:
        INGEST_ENABLED: Master switch for scheduled ingest jobs.
        INGEST_SCHEDULE_HOURS: Base interval in hours.
            ClinVar/OMIM = 7x (weekly), HPO/Orphanet = 30x (monthly).
        OMIM_ENABLED: Whether to schedule OMIM refresh.
        HPO_ENABLED: Whether to schedule HPO refresh.
        ORPHANET_ENABLED: Whether to schedule Orphanet refresh.
        CLINVAR_ENABLED: Whether to schedule ClinVar refresh.
    """

    INGEST_ENABLED: bool = True
    INGEST_SCHEDULE_HOURS: int = 24  # base unit
    OMIM_ENABLED: bool = True
    HPO_ENABLED: bool = True
    ORPHANET_ENABLED: bool = True
    CLINVAR_ENABLED: bool = True


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

    class RareDiseaseScheduler:
        """Background scheduler for periodic rare disease data ingestion.

        Manages four recurring jobs:
          1. ClinVar variant refresh (weekly)
          2. OMIM disease refresh (weekly)
          3. HPO phenotype refresh (monthly)
          4. Orphanet disease refresh (monthly)

        Usage::

            from src.scheduler import RareDiseaseScheduler, RareDiseaseSchedulerSettings

            settings = RareDiseaseSchedulerSettings(INGEST_ENABLED=True)
            scheduler = RareDiseaseScheduler(settings=settings)
            scheduler.start()
            ...
            scheduler.stop()
        """

        def __init__(
            self,
            settings: Optional[RareDiseaseSchedulerSettings] = None,
            collection_manager: Any = None,
            embedder: Any = None,
        ):
            self.settings = settings or RareDiseaseSchedulerSettings()
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

            hours = self.settings.INGEST_SCHEDULE_HOURS

            if self.settings.CLINVAR_ENABLED:
                self.scheduler.add_job(
                    self._run_clinvar_ingest,
                    "interval",
                    hours=hours * 7,  # weekly
                    id="clinvar_ingest",
                    name="ClinVar variant refresh (weekly)",
                    replace_existing=True,
                )

            if self.settings.OMIM_ENABLED:
                self.scheduler.add_job(
                    self._run_omim_ingest,
                    "interval",
                    hours=hours * 7,  # weekly
                    id="omim_ingest",
                    name="OMIM disease refresh (weekly)",
                    replace_existing=True,
                )

            if self.settings.HPO_ENABLED:
                self.scheduler.add_job(
                    self._run_hpo_ingest,
                    "interval",
                    hours=hours * 30,  # monthly
                    id="hpo_ingest",
                    name="HPO phenotype refresh (monthly)",
                    replace_existing=True,
                )

            if self.settings.ORPHANET_ENABLED:
                self.scheduler.add_job(
                    self._run_orphanet_ingest,
                    "interval",
                    hours=hours * 30,  # monthly
                    id="orphanet_ingest",
                    name="Orphanet disease refresh (monthly)",
                    replace_existing=True,
                )

            self.scheduler.start()
            self.logger.info(
                "RareDiseaseScheduler started -- "
                "ClinVar/OMIM every %dh (weekly), "
                "HPO/Orphanet every %dh (monthly)",
                hours * 7,
                hours * 30,
            )

        def stop(self) -> None:
            """Gracefully shut down the background scheduler."""
            if self.scheduler.running:
                self.scheduler.shutdown(wait=False)
                self.logger.info("RareDiseaseScheduler stopped")

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
            """Trigger an immediate manual ingest for the specified source."""
            dispatch = {
                "clinvar": self._run_clinvar_ingest,
                "omim": self._run_omim_ingest,
                "hpo": self._run_hpo_ingest,
                "orphanet": self._run_orphanet_ingest,
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

        def _run_ingest_job(self, source: str, parser_factory: Any) -> None:
            """Generic ingest job runner."""
            job_status = IngestJobStatus(
                job_id=f"{source}_{int(time.time())}",
                source=source,
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting %s refresh", source)
            start = time.time()

            try:
                parser = parser_factory(
                    collection_manager=self.collection_manager,
                    embedder=self.embedder,
                )
                records, stats = parser.run()
                elapsed = time.time() - start
                self._last_run_time = time.time()
                count = len(records)

                MetricsCollector.record_ingest(
                    source=source,
                    duration=elapsed,
                    record_count=count,
                    collection=f"rd_{source}",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    "Scheduler: %s refresh complete -- %d records in %.1fs",
                    source, count, elapsed,
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = f"{source} parser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: %s ingest skipped -- parser module not available",
                    source,
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source=source).inc()
                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.error("Scheduler: %s refresh failed -- %s", source, exc)

            self._job_history.append(job_status)

        def _run_clinvar_ingest(self) -> None:
            """Run the ClinVar variant ingest pipeline."""
            self._run_ingest_job("clinvar", _get_clinvar_parser_class)

        def _run_omim_ingest(self) -> None:
            """Run the OMIM disease ingest pipeline."""
            from .ingest.omim_parser import OMIMParser
            self._run_ingest_job("omim", lambda **kw: OMIMParser(**kw))

        def _run_hpo_ingest(self) -> None:
            """Run the HPO phenotype ingest pipeline."""
            from .ingest.hpo_parser import HPOParser
            self._run_ingest_job("hpo", lambda **kw: HPOParser(**kw))

        def _run_orphanet_ingest(self) -> None:
            """Run the Orphanet disease ingest pipeline."""
            from .ingest.orphanet_parser import OrphanetParser
            self._run_ingest_job("orphanet", lambda **kw: OrphanetParser(**kw))

else:
    # -- No-op stub when apscheduler is not installed --

    class RareDiseaseScheduler:  # type: ignore[no-redef]
        """No-op scheduler stub (apscheduler not installed)."""

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            logger.warning(
                "apscheduler is not installed -- RareDiseaseScheduler is a no-op. "
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
                "message": "Scheduler unavailable -- apscheduler is not installed.",
            }


def _get_clinvar_parser_class(**kwargs: Any) -> Any:
    """Stub factory for ClinVar parser (not yet implemented)."""
    raise ImportError("ClinVar parser not yet implemented")
