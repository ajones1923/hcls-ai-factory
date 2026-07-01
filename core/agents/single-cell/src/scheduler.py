"""Automated ingest scheduler for the Single-Cell Intelligence Agent.

Periodically refreshes CellxGene cell type data, CellMarker/PanglaoDB
marker gene databases, and single-cell literature so the knowledge base
stays current without manual intervention.

Uses APScheduler's BackgroundScheduler so jobs run in a daemon thread
alongside the FastAPI / Streamlit application.

Default cadence:
  - CellxGene cell types:    monthly (720h)
  - Marker gene databases:   monthly (720h)
  - Single-cell literature:  weekly (168h)

If ``apscheduler`` is not installed the module exports a no-op
``SingleCellScheduler`` stub so dependent code can import unconditionally.

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
class SingleCellSchedulerSettings:
    """Configuration for the single-cell ingest scheduler.

    Attributes:
        INGEST_ENABLED: Master switch for scheduled ingest jobs.
        CELLXGENE_SCHEDULE_HOURS: Interval for CellxGene refresh (default: monthly).
        MARKERS_SCHEDULE_HOURS: Interval for marker DB refresh (default: monthly).
        LITERATURE_SCHEDULE_HOURS: Interval for literature refresh (default: weekly).
        MAX_CELLXGENE_RESULTS: Maximum cell type records per refresh.
        MAX_MARKER_RESULTS: Maximum marker records per refresh.
        MAX_LITERATURE_RESULTS: Maximum literature records per refresh.
    """

    INGEST_ENABLED: bool = True
    CELLXGENE_SCHEDULE_HOURS: int = 720  # monthly
    MARKERS_SCHEDULE_HOURS: int = 720  # monthly
    LITERATURE_SCHEDULE_HOURS: int = 168  # weekly
    MAX_CELLXGENE_RESULTS: int = 500
    MAX_MARKER_RESULTS: int = 300
    MAX_LITERATURE_RESULTS: int = 200


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

    class SingleCellScheduler:
        """Background scheduler for periodic single-cell data ingestion.

        Manages three recurring jobs:
          1. CellxGene cell type refresh (monthly)
          2. CellMarker/PanglaoDB marker gene refresh (monthly)
          3. Single-cell literature refresh (weekly)

        Usage::

            from src.scheduler import SingleCellScheduler, SingleCellSchedulerSettings

            settings = SingleCellSchedulerSettings(INGEST_ENABLED=True)
            scheduler = SingleCellScheduler(
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
            settings: Optional[SingleCellSchedulerSettings] = None,
            collection_manager: Any = None,
            embedder: Any = None,
        ):
            self.settings = settings or SingleCellSchedulerSettings()
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
                self._run_cellxgene_ingest,
                "interval",
                hours=self.settings.CELLXGENE_SCHEDULE_HOURS,
                id="sc_cellxgene_ingest",
                name="CellxGene cell type refresh",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_marker_ingest,
                "interval",
                hours=self.settings.MARKERS_SCHEDULE_HOURS,
                id="sc_marker_ingest",
                name="CellMarker/PanglaoDB marker gene refresh",
                replace_existing=True,
            )

            self.scheduler.add_job(
                self._run_literature_ingest,
                "interval",
                hours=self.settings.LITERATURE_SCHEDULE_HOURS,
                id="sc_literature_ingest",
                name="Single-cell literature refresh",
                replace_existing=True,
            )

            self.scheduler.start()
            self.logger.info(
                "SingleCellScheduler started -- "
                "CellxGene every %dh, Markers every %dh, Literature every %dh",
                self.settings.CELLXGENE_SCHEDULE_HOURS,
                self.settings.MARKERS_SCHEDULE_HOURS,
                self.settings.LITERATURE_SCHEDULE_HOURS,
            )

        def stop(self) -> None:
            """Gracefully shut down the background scheduler."""
            if self.scheduler.running:
                self.scheduler.shutdown(wait=False)
                self.logger.info("SingleCellScheduler stopped")

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
                "cellxgene_schedule_hours": self.settings.CELLXGENE_SCHEDULE_HOURS,
                "markers_schedule_hours": self.settings.MARKERS_SCHEDULE_HOURS,
                "literature_schedule_hours": self.settings.LITERATURE_SCHEDULE_HOURS,
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
                "cellxgene": self._run_cellxgene_ingest,
                "markers": self._run_marker_ingest,
                "literature": self._run_literature_ingest,
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

        def _run_cellxgene_ingest(self) -> None:
            """Run the CellxGene cell type ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"cellxgene_{int(time.time())}",
                source="cellxgene",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting CellxGene cell type refresh")
            start = time.time()

            try:
                from .ingest.cellxgene_parser import CellxGeneParser

                parser = CellxGeneParser(
                    collection_manager=self.collection_manager,
                    embedder=self.embedder,
                )
                records, stats = parser.run()
                elapsed = time.time() - start
                self._last_run_time = time.time()
                count = len(records)

                MetricsCollector.record_ingest(
                    source="cellxgene",
                    duration=elapsed,
                    record_count=count,
                    collection="sc_cell_types",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    "Scheduler: CellxGene refresh complete -- "
                    "%d records in %.1fs",
                    count, elapsed,
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "CellxGeneParser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: CellxGene ingest skipped -- "
                    "cellxgene_parser module not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="cellxgene").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    "Scheduler: CellxGene refresh failed -- %s", exc
                )

            self._job_history.append(job_status)

        def _run_marker_ingest(self) -> None:
            """Run the CellMarker/PanglaoDB marker gene ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"markers_{int(time.time())}",
                source="cellmarker",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting marker gene database refresh")
            start = time.time()

            try:
                from .ingest.marker_parser import MarkerParser

                parser = MarkerParser(
                    collection_manager=self.collection_manager,
                    embedder=self.embedder,
                )
                records, stats = parser.run()
                elapsed = time.time() - start
                self._last_run_time = time.time()
                count = len(records)

                MetricsCollector.record_ingest(
                    source="cellmarker",
                    duration=elapsed,
                    record_count=count,
                    collection="sc_markers",
                    success=True,
                )

                job_status.status = "success"
                job_status.records_ingested = count
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.info(
                    "Scheduler: Marker gene refresh complete -- "
                    "%d records in %.1fs",
                    count, elapsed,
                )

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "MarkerParser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Marker ingest skipped -- parser not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="cellmarker").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    "Scheduler: Marker gene refresh failed -- %s", exc
                )

            self._job_history.append(job_status)

        def _run_literature_ingest(self) -> None:
            """Run the single-cell literature ingest pipeline."""
            job_status = IngestJobStatus(
                job_id=f"literature_{int(time.time())}",
                source="literature",
                status="running",
                started_at=datetime.now(timezone.utc).isoformat(),
            )

            self.logger.info("Scheduler: starting single-cell literature refresh")
            start = time.time()

            try:
                # Literature parser would be a separate module
                raise ImportError("Single-cell literature parser not yet implemented")

            except ImportError:
                elapsed = time.time() - start
                job_status.status = "error"
                job_status.error_message = "Literature parser not available"
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()
                self.logger.warning(
                    "Scheduler: Literature ingest skipped -- parser not available"
                )

            except Exception as exc:
                elapsed = time.time() - start
                INGEST_ERRORS.labels(source="literature").inc()

                job_status.status = "error"
                job_status.error_message = str(exc)
                job_status.duration_seconds = elapsed
                job_status.completed_at = datetime.now(timezone.utc).isoformat()

                self.logger.error(
                    "Scheduler: Literature refresh failed -- %s", exc
                )

            self._job_history.append(job_status)

else:
    # -- No-op stub when apscheduler is not installed --

    class SingleCellScheduler:  # type: ignore[no-redef]
        """No-op scheduler stub (apscheduler not installed).

        All methods are safe to call but perform no work. Install
        apscheduler to enable scheduled ingest::

            pip install apscheduler>=3.10.0
        """

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            logger.warning(
                "apscheduler is not installed -- SingleCellScheduler is a no-op. "
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
                "cellxgene_schedule_hours": 0,
                "markers_schedule_hours": 0,
                "literature_schedule_hours": 0,
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
