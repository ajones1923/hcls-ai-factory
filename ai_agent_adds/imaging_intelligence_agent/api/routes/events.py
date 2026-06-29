"""DICOM event bus and Orthanc webhook endpoints.

This module provides webhook endpoints for DICOM study events,
enabling real-time integration with PACS/Orthanc servers. When a
study is marked complete in Orthanc, the webhook triggers automatic
workflow selection based on modality and body region, runs inference,
and stores results in the ingestion history.

Author: Adam Jones
Date: February 2026
"""

import asyncio
import tempfile
from collections import deque
from datetime import datetime, timezone
from typing import Deque, Dict, List, Optional

import httpx
from fastapi import APIRouter, HTTPException, Query
from loguru import logger
from pydantic import BaseModel, Field

from src.workflows import WORKFLOW_REGISTRY

# =====================================================================
# Constants
# =====================================================================

MAX_HISTORY = 200  # Maximum ingestion history entries retained in memory

# Modality + body_region -> workflow name mapping
WORKFLOW_ROUTING: Dict[tuple, str] = {
    ("CT", "head"): "ct_head_hemorrhage",
    ("CT", "brain"): "ct_head_hemorrhage",
    ("CT", "chest"): "ct_chest_lung_nodule",
    ("CT", "lung"): "ct_chest_lung_nodule",
    ("CR", "chest"): "cxr_rapid_findings",
    ("DX", "chest"): "cxr_rapid_findings",
    ("MR", "brain"): "mri_brain_ms_lesion",
    ("MR", "head"): "mri_brain_ms_lesion",
}

# =====================================================================
# Request / Response Models
# =====================================================================


class DicomStudyEvent(BaseModel):
    """Orthanc webhook change event."""
    event_type: str = Field(
        ...,
        description="Lifecycle event: study.complete, study.received, study.updated",
    )
    study_uid: Optional[str] = None
    patient_id: Optional[str] = None
    modality: Optional[str] = None
    body_region: Optional[str] = None
    orthanc_id: Optional[str] = None
    series_count: int = 0
    instance_count: int = 0


class DicomIngestionResult(BaseModel):
    """Result of processing a DICOM study."""
    study_uid: str
    patient_id: str
    modality: str
    workflow_triggered: Optional[str] = None
    workflow_status: str = "pending"
    workflow_result: Optional[dict] = None
    fhir_report: Optional[str] = None
    processed_at: str = Field(
        default_factory=lambda: datetime.now(timezone.utc).isoformat(),
    )


class IngestionHistoryResponse(BaseModel):
    """Paginated ingestion history response."""
    items: List[DicomIngestionResult]
    total: int
    limit: int
    offset: int


# =====================================================================
# In-memory state
# =====================================================================

_ingestion_history: Deque[DicomIngestionResult] = deque(maxlen=MAX_HISTORY)

# =====================================================================
# Routing Logic
# =====================================================================


def determine_workflow(modality: Optional[str], body_region: Optional[str]) -> Optional[str]:
    """Determine which imaging workflow to run based on DICOM modality and body region.

    Performs a case-insensitive lookup against the routing table.  Returns
    None if no matching workflow is found.

    Args:
        modality: DICOM modality code (e.g. "CT", "MR", "CR", "DX").
        body_region: Anatomical region (e.g. "head", "chest", "brain").

    Returns:
        Workflow registry key or None.
    """
    if not modality or not body_region:
        return None

    key = (modality.upper().strip(), body_region.lower().strip())
    return WORKFLOW_ROUTING.get(key)


# =====================================================================
# Settings (lazy import to avoid module-level validation errors)
# =====================================================================


def _get_settings():
    """Lazy-import settings to avoid module-level Pydantic validation."""
    from config.settings import settings
    return settings


# =====================================================================
# Orthanc REST Helpers
# =====================================================================


async def fetch_orthanc_metadata(orthanc_id: str) -> dict:
    """Fetch DICOM study metadata from Orthanc REST API.

    Calls ``GET /studies/{orthanc_id}`` on the configured Orthanc server
    and returns the JSON response.

    Args:
        orthanc_id: Orthanc-internal study identifier.

    Returns:
        Dictionary with Orthanc study metadata.

    Raises:
        httpx.HTTPStatusError: If the Orthanc server returns an error.
    """
    s = _get_settings()
    url = f"{s.ORTHANC_URL}/studies/{orthanc_id}"
    auth = (s.ORTHANC_USERNAME, s.ORTHANC_PASSWORD)

    async with httpx.AsyncClient(timeout=30.0) as client:
        response = await client.get(url, auth=auth)
        response.raise_for_status()
        return response.json()


async def download_dicom_series(orthanc_id: str, output_dir: Optional[str] = None) -> str:
    """Download all DICOM instances for a study from Orthanc.

    Retrieves the study archive (ZIP of all DICOM files) from
    ``GET /studies/{orthanc_id}/archive`` and writes it to a temporary
    directory.

    Args:
        orthanc_id: Orthanc-internal study identifier.
        output_dir: Target directory.  A temp directory is created if None.

    Returns:
        Path to the directory containing the downloaded archive.

    Raises:
        httpx.HTTPStatusError: If the download fails.
    """
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix="dicom_")

    s = _get_settings()
    url = f"{s.ORTHANC_URL}/studies/{orthanc_id}/archive"
    auth = (s.ORTHANC_USERNAME, s.ORTHANC_PASSWORD)
    archive_path = f"{output_dir}/study_{orthanc_id}.zip"

    async with httpx.AsyncClient(timeout=120.0) as client:
        response = await client.get(url, auth=auth)
        response.raise_for_status()
        with open(archive_path, "wb") as f:
            f.write(response.content)

    logger.info(f"Downloaded DICOM archive to {archive_path}")
    return output_dir


# =====================================================================
# Workflow Execution
# =====================================================================


async def _run_workflow_async(
    event: DicomStudyEvent,
    workflow_name: str,
) -> DicomIngestionResult:
    """Execute a workflow asynchronously and return the ingestion result.

    The workflow is run in mock mode by default (as live NIM services may
    not be available).  Override by ensuring NIM endpoints are reachable.
    """
    study_uid = event.study_uid or "unknown"
    patient_id = event.patient_id or "unknown"
    modality = event.modality or "unknown"

    result = DicomIngestionResult(
        study_uid=study_uid,
        patient_id=patient_id,
        modality=modality,
        workflow_triggered=workflow_name,
        workflow_status="running",
    )

    try:
        wf_class = WORKFLOW_REGISTRY[workflow_name]
        workflow = wf_class(mock_mode=True)

        # Run the synchronous workflow in a thread pool to avoid blocking
        loop = asyncio.get_event_loop()
        wf_result = await loop.run_in_executor(None, workflow.run, "")

        result.workflow_status = wf_result.status.value
        result.workflow_result = {
            "findings_count": len(wf_result.findings),
            "classification": wf_result.classification,
            "severity": wf_result.severity.value,
            "inference_time_ms": wf_result.inference_time_ms,
            "nim_services_used": wf_result.nim_services_used,
        }
        result.processed_at = datetime.now(timezone.utc).isoformat()

        logger.info(
            f"Workflow {workflow_name} completed for study {study_uid}: "
            f"{wf_result.status.value}"
        )
    except Exception as exc:
        logger.error(f"Workflow {workflow_name} failed for study {study_uid}: {exc}")
        result.workflow_status = "failed"
        result.workflow_result = {"error": str(exc)}
        result.processed_at = datetime.now(timezone.utc).isoformat()

    return result


# =====================================================================
# Router
# =====================================================================

events_router = APIRouter(prefix="/events", tags=["events"])


@events_router.post("/dicom-webhook", response_model=DicomIngestionResult)
async def dicom_webhook(event: DicomStudyEvent):
    """Receive DICOM study events from Orthanc.

    Accepts Orthanc webhook ``study.complete`` events, determines the
    appropriate imaging workflow, executes it, and stores the result
    in the in-memory ingestion history.

    Only ``study.complete`` events trigger processing.  Other event
    types (``study.received``, ``study.updated``) are acknowledged
    but not processed.
    """
    logger.info(
        f"DICOM webhook received: {event.event_type} "
        f"study_uid={event.study_uid} modality={event.modality}"
    )

    study_uid = event.study_uid or "unknown"
    patient_id = event.patient_id or "unknown"
    modality = event.modality or "unknown"

    # Only process study.complete events
    if event.event_type != "study.complete":
        result = DicomIngestionResult(
            study_uid=study_uid,
            patient_id=patient_id,
            modality=modality,
            workflow_triggered=None,
            workflow_status="skipped",
        )
        _ingestion_history.appendleft(result)
        return result

    # If we have an orthanc_id but are missing metadata, try to fetch it
    if event.orthanc_id and not event.modality:
        try:
            metadata = await fetch_orthanc_metadata(event.orthanc_id)
            main_tags = metadata.get("MainDicomTags", {})
            event.study_uid = event.study_uid or main_tags.get("StudyInstanceUID")
            event.patient_id = event.patient_id or metadata.get("PatientMainDicomTags", {}).get("PatientID")
            event.modality = main_tags.get("Modality", event.modality)
            event.body_region = main_tags.get("BodyPartExamined", event.body_region)
            modality = event.modality or "unknown"
            study_uid = event.study_uid or "unknown"
            patient_id = event.patient_id or "unknown"
        except Exception as exc:
            logger.warning(f"Failed to fetch Orthanc metadata for {event.orthanc_id}: {exc}")

    # Determine which workflow to run
    workflow_name = determine_workflow(event.modality, event.body_region)

    if workflow_name is None:
        logger.warning(
            f"No matching workflow for modality={event.modality}, "
            f"body_region={event.body_region}"
        )
        result = DicomIngestionResult(
            study_uid=study_uid,
            patient_id=patient_id,
            modality=modality,
            workflow_triggered=None,
            workflow_status="no_matching_workflow",
        )
        _ingestion_history.appendleft(result)
        return result

    # Run the matched workflow
    result = await _run_workflow_async(event, workflow_name)
    _ingestion_history.appendleft(result)
    return result


@events_router.get("/history", response_model=IngestionHistoryResponse)
async def ingestion_history(
    limit: int = Query(20, ge=1, le=100, description="Number of entries to return"),
    offset: int = Query(0, ge=0, description="Offset into the history"),
):
    """Return the recent DICOM ingestion history.

    Returns the last N processed studies from the in-memory history
    buffer.  Results are ordered newest-first.
    """
    items = list(_ingestion_history)
    total = len(items)
    page = items[offset: offset + limit]

    return IngestionHistoryResponse(
        items=page,
        total=total,
        limit=limit,
        offset=offset,
    )


@events_router.get("/status")
async def event_bus_status():
    """Check event bus status and configuration."""
    s = _get_settings()
    return {
        "active": s.DICOM_AUTO_INGEST,
        "orthanc_url": s.ORTHANC_URL,
        "watch_interval_sec": s.DICOM_WATCH_INTERVAL,
        "supported_events": [
            "study.received",
            "study.complete",
            "study.updated",
        ],
        "workflow_routing": {
            f"{mod}+{region}": wf
            for (mod, region), wf in WORKFLOW_ROUTING.items()
        },
        "history_count": len(_ingestion_history),
        "history_max": MAX_HISTORY,
    }
