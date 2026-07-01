"""Tests for DICOM auto-ingestion pipeline.

Covers:
  - Workflow routing logic (determine_workflow)
  - Webhook event processing (dicom_webhook endpoint)
  - Ingestion history endpoint
  - Event bus status endpoint
  - Orthanc metadata fetch (mocked)
  - Orthanc DICOM download (mocked)
  - DicomWatcher lifecycle (start / stop / polling)
  - Edge cases: unknown modality, missing fields, non-complete events

Author: Adam Jones
Date: February 2026
"""

import threading
import time
from collections import deque
from unittest.mock import AsyncMock, MagicMock, patch

import httpx
import pytest
from fastapi.testclient import TestClient

from api.routes.events import (
    DicomIngestionResult,
    DicomStudyEvent,
    WORKFLOW_ROUTING,
    _ingestion_history,
    determine_workflow,
    download_dicom_series,
    events_router,
    fetch_orthanc_metadata,
)
from src.ingest.dicom_watcher import DicomWatcher


# =====================================================================
# Mock settings for tests (avoids Pydantic env validation)
# =====================================================================


class _MockSettings:
    """Minimal settings stand-in for testing."""
    ORTHANC_URL = "http://localhost:8042"
    ORTHANC_USERNAME = "admin"
    ORTHANC_PASSWORD = "orthanc"
    DICOM_AUTO_INGEST = False
    DICOM_WATCH_INTERVAL = 5


# =====================================================================
# Fixtures
# =====================================================================


@pytest.fixture(autouse=True)
def _clear_history():
    """Ensure ingestion history is empty before and after each test."""
    _ingestion_history.clear()
    yield
    _ingestion_history.clear()


@pytest.fixture
def mock_settings():
    """Patch _get_settings to return a mock settings object."""
    with patch("api.routes.events._get_settings", return_value=_MockSettings()):
        yield _MockSettings()


@pytest.fixture
def client(mock_settings):
    """FastAPI TestClient wired to the events router only."""
    from fastapi import FastAPI

    app = FastAPI()
    app.include_router(events_router)
    return TestClient(app)


@pytest.fixture
def sample_complete_event():
    """A study.complete event for a CT head scan."""
    return {
        "event_type": "study.complete",
        "study_uid": "1.2.840.113619.2.55.3.1234",
        "patient_id": "PAT-001",
        "modality": "CT",
        "body_region": "head",
        "orthanc_id": "abc-123",
        "series_count": 3,
        "instance_count": 450,
    }


@pytest.fixture
def sample_received_event():
    """A study.received event (should be skipped)."""
    return {
        "event_type": "study.received",
        "study_uid": "1.2.840.113619.2.55.3.9999",
        "patient_id": "PAT-002",
        "modality": "MR",
        "body_region": "brain",
    }


# =====================================================================
# Test: determine_workflow routing
# =====================================================================


class TestDetermineWorkflow:
    """Tests for the modality + body_region -> workflow routing logic."""

    def test_ct_head_routes_to_hemorrhage(self):
        assert determine_workflow("CT", "head") == "ct_head_hemorrhage"

    def test_ct_brain_routes_to_hemorrhage(self):
        assert determine_workflow("CT", "brain") == "ct_head_hemorrhage"

    def test_ct_chest_routes_to_lung_nodule(self):
        assert determine_workflow("CT", "chest") == "ct_chest_lung_nodule"

    def test_ct_lung_routes_to_lung_nodule(self):
        assert determine_workflow("CT", "lung") == "ct_chest_lung_nodule"

    def test_cr_chest_routes_to_cxr(self):
        assert determine_workflow("CR", "chest") == "cxr_rapid_findings"

    def test_dx_chest_routes_to_cxr(self):
        assert determine_workflow("DX", "chest") == "cxr_rapid_findings"

    def test_mr_brain_routes_to_ms_lesion(self):
        assert determine_workflow("MR", "brain") == "mri_brain_ms_lesion"

    def test_mr_head_routes_to_ms_lesion(self):
        assert determine_workflow("MR", "head") == "mri_brain_ms_lesion"

    def test_case_insensitive_modality(self):
        assert determine_workflow("ct", "head") == "ct_head_hemorrhage"

    def test_case_insensitive_body_region(self):
        assert determine_workflow("CT", "HEAD") == "ct_head_hemorrhage"

    def test_unknown_modality_returns_none(self):
        assert determine_workflow("US", "abdomen") is None

    def test_unknown_body_region_returns_none(self):
        assert determine_workflow("CT", "pelvis") is None

    def test_none_modality_returns_none(self):
        assert determine_workflow(None, "head") is None

    def test_none_body_region_returns_none(self):
        assert determine_workflow("CT", None) is None

    def test_both_none_returns_none(self):
        assert determine_workflow(None, None) is None

    def test_whitespace_stripped(self):
        assert determine_workflow("  CT  ", "  head  ") == "ct_head_hemorrhage"

    def test_all_routing_entries_map_to_valid_workflows(self):
        """Every entry in WORKFLOW_ROUTING must reference a registered workflow."""
        from src.workflows import WORKFLOW_REGISTRY

        for key, wf_name in WORKFLOW_ROUTING.items():
            assert wf_name in WORKFLOW_REGISTRY, (
                f"Routing entry {key} -> {wf_name} not in WORKFLOW_REGISTRY"
            )


# =====================================================================
# Test: Webhook event processing
# =====================================================================


class TestDicomWebhook:
    """Tests for the POST /events/dicom-webhook endpoint."""

    def test_study_complete_triggers_workflow(self, client, sample_complete_event):
        resp = client.post("/events/dicom-webhook", json=sample_complete_event)
        assert resp.status_code == 200
        data = resp.json()
        assert data["study_uid"] == "1.2.840.113619.2.55.3.1234"
        assert data["workflow_triggered"] == "ct_head_hemorrhage"
        assert data["workflow_status"] == "completed"
        assert data["workflow_result"] is not None

    def test_study_received_is_skipped(self, client, sample_received_event):
        resp = client.post("/events/dicom-webhook", json=sample_received_event)
        assert resp.status_code == 200
        data = resp.json()
        assert data["workflow_status"] == "skipped"
        assert data["workflow_triggered"] is None

    def test_study_updated_is_skipped(self, client):
        event = {
            "event_type": "study.updated",
            "study_uid": "1.2.3.4.5",
            "patient_id": "PAT-X",
            "modality": "CT",
            "body_region": "chest",
        }
        resp = client.post("/events/dicom-webhook", json=event)
        assert resp.status_code == 200
        assert resp.json()["workflow_status"] == "skipped"

    def test_no_matching_workflow(self, client):
        event = {
            "event_type": "study.complete",
            "study_uid": "1.2.3.4.5",
            "patient_id": "PAT-003",
            "modality": "US",
            "body_region": "abdomen",
        }
        resp = client.post("/events/dicom-webhook", json=event)
        assert resp.status_code == 200
        data = resp.json()
        assert data["workflow_status"] == "no_matching_workflow"
        assert data["workflow_triggered"] is None

    def test_missing_modality_no_orthanc_id(self, client):
        """Missing modality with no orthanc_id -> no matching workflow."""
        event = {
            "event_type": "study.complete",
            "study_uid": "1.2.3",
            "patient_id": "PAT-X",
        }
        resp = client.post("/events/dicom-webhook", json=event)
        assert resp.status_code == 200
        assert resp.json()["workflow_status"] == "no_matching_workflow"

    def test_webhook_stores_result_in_history(self, client, sample_complete_event):
        client.post("/events/dicom-webhook", json=sample_complete_event)
        assert len(_ingestion_history) == 1
        assert _ingestion_history[0].study_uid == "1.2.840.113619.2.55.3.1234"

    def test_multiple_events_accumulate_history(self, client, sample_complete_event, sample_received_event):
        client.post("/events/dicom-webhook", json=sample_complete_event)
        client.post("/events/dicom-webhook", json=sample_received_event)
        assert len(_ingestion_history) == 2

    def test_ct_chest_workflow(self, client):
        event = {
            "event_type": "study.complete",
            "study_uid": "1.2.3.4.5.6",
            "patient_id": "PAT-004",
            "modality": "CT",
            "body_region": "chest",
        }
        resp = client.post("/events/dicom-webhook", json=event)
        assert resp.status_code == 200
        data = resp.json()
        assert data["workflow_triggered"] == "ct_chest_lung_nodule"
        assert data["workflow_status"] == "completed"

    def test_cxr_workflow(self, client):
        event = {
            "event_type": "study.complete",
            "study_uid": "1.2.3.4.5.7",
            "patient_id": "PAT-005",
            "modality": "CR",
            "body_region": "chest",
        }
        resp = client.post("/events/dicom-webhook", json=event)
        assert resp.status_code == 200
        assert resp.json()["workflow_triggered"] == "cxr_rapid_findings"

    def test_mri_brain_workflow(self, client):
        event = {
            "event_type": "study.complete",
            "study_uid": "1.2.3.4.5.8",
            "patient_id": "PAT-006",
            "modality": "MR",
            "body_region": "brain",
        }
        resp = client.post("/events/dicom-webhook", json=event)
        assert resp.status_code == 200
        assert resp.json()["workflow_triggered"] == "mri_brain_ms_lesion"


# =====================================================================
# Test: Ingestion history endpoint
# =====================================================================


class TestIngestionHistory:
    """Tests for the GET /events/history endpoint."""

    def test_empty_history(self, client):
        resp = client.get("/events/history")
        assert resp.status_code == 200
        data = resp.json()
        assert data["items"] == []
        assert data["total"] == 0

    def test_history_after_events(self, client, sample_complete_event):
        client.post("/events/dicom-webhook", json=sample_complete_event)
        resp = client.get("/events/history")
        data = resp.json()
        assert data["total"] == 1
        assert len(data["items"]) == 1
        assert data["items"][0]["study_uid"] == "1.2.840.113619.2.55.3.1234"

    def test_history_pagination_limit(self, client, sample_complete_event):
        # Post 3 events
        for _ in range(3):
            client.post("/events/dicom-webhook", json=sample_complete_event)
        resp = client.get("/events/history", params={"limit": 2})
        data = resp.json()
        assert data["total"] == 3
        assert len(data["items"]) == 2
        assert data["limit"] == 2

    def test_history_pagination_offset(self, client, sample_complete_event):
        for _ in range(3):
            client.post("/events/dicom-webhook", json=sample_complete_event)
        resp = client.get("/events/history", params={"limit": 10, "offset": 2})
        data = resp.json()
        assert len(data["items"]) == 1
        assert data["offset"] == 2

    def test_history_newest_first(self, client):
        """Most recent event should be first in history."""
        event1 = {
            "event_type": "study.complete",
            "study_uid": "FIRST",
            "patient_id": "P1",
            "modality": "CT",
            "body_region": "head",
        }
        event2 = {
            "event_type": "study.complete",
            "study_uid": "SECOND",
            "patient_id": "P2",
            "modality": "CT",
            "body_region": "chest",
        }
        client.post("/events/dicom-webhook", json=event1)
        client.post("/events/dicom-webhook", json=event2)
        resp = client.get("/events/history")
        items = resp.json()["items"]
        assert items[0]["study_uid"] == "SECOND"
        assert items[1]["study_uid"] == "FIRST"


# =====================================================================
# Test: Event bus status endpoint
# =====================================================================


class TestEventBusStatus:
    """Tests for the GET /events/status endpoint."""

    def test_status_returns_ok(self, client):
        resp = client.get("/events/status")
        assert resp.status_code == 200

    def test_status_contains_routing(self, client):
        data = client.get("/events/status").json()
        assert "workflow_routing" in data
        assert "CT+head" in data["workflow_routing"]

    def test_status_contains_supported_events(self, client):
        data = client.get("/events/status").json()
        assert "study.complete" in data["supported_events"]

    def test_status_history_count(self, client, sample_complete_event):
        client.post("/events/dicom-webhook", json=sample_complete_event)
        data = client.get("/events/status").json()
        assert data["history_count"] == 1


# =====================================================================
# Test: Orthanc metadata fetch (mocked)
# =====================================================================


class TestOrthancMetadataFetch:
    """Tests for fetch_orthanc_metadata with mocked HTTP."""

    @pytest.mark.asyncio
    async def test_fetch_metadata_success(self):
        mock_json = {
            "ID": "abc-123",
            "MainDicomTags": {
                "StudyInstanceUID": "1.2.3.4.5",
                "Modality": "CT",
                "BodyPartExamined": "HEAD",
            },
            "PatientMainDicomTags": {"PatientID": "PAT-001"},
            "Series": ["series-1", "series-2"],
        }

        mock_response = MagicMock()
        mock_response.json.return_value = mock_json
        mock_response.raise_for_status = MagicMock()

        with patch("api.routes.events._get_settings", return_value=_MockSettings()), \
             patch("api.routes.events.httpx.AsyncClient") as MockClient:
            instance = AsyncMock()
            instance.get = AsyncMock(return_value=mock_response)
            instance.__aenter__ = AsyncMock(return_value=instance)
            instance.__aexit__ = AsyncMock(return_value=False)
            MockClient.return_value = instance

            result = await fetch_orthanc_metadata("abc-123")
            assert result["ID"] == "abc-123"
            assert result["MainDicomTags"]["Modality"] == "CT"

    @pytest.mark.asyncio
    async def test_fetch_metadata_error_raises(self):
        with patch("api.routes.events._get_settings", return_value=_MockSettings()), \
             patch("api.routes.events.httpx.AsyncClient") as MockClient:
            instance = AsyncMock()
            instance.get = AsyncMock(side_effect=httpx.HTTPStatusError(
                "404", request=MagicMock(), response=MagicMock()
            ))
            instance.__aenter__ = AsyncMock(return_value=instance)
            instance.__aexit__ = AsyncMock(return_value=False)
            MockClient.return_value = instance

            with pytest.raises(httpx.HTTPStatusError):
                await fetch_orthanc_metadata("nonexistent")


# =====================================================================
# Test: DICOM download (mocked)
# =====================================================================


class TestDicomDownload:
    """Tests for download_dicom_series with mocked HTTP."""

    @pytest.mark.asyncio
    async def test_download_creates_file(self, tmp_path):
        mock_response = MagicMock()
        mock_response.content = b"PK\x03\x04fake-zip-data"
        mock_response.raise_for_status = MagicMock()

        with patch("api.routes.events._get_settings", return_value=_MockSettings()), \
             patch("api.routes.events.httpx.AsyncClient") as MockClient:
            instance = AsyncMock()
            instance.get = AsyncMock(return_value=mock_response)
            instance.__aenter__ = AsyncMock(return_value=instance)
            instance.__aexit__ = AsyncMock(return_value=False)
            MockClient.return_value = instance

            output = await download_dicom_series("abc-123", str(tmp_path))
            assert output == str(tmp_path)
            archive = tmp_path / "study_abc-123.zip"
            assert archive.exists()
            assert archive.read_bytes() == b"PK\x03\x04fake-zip-data"


# =====================================================================
# Test: DicomWatcher lifecycle
# =====================================================================


class TestDicomWatcher:
    """Tests for the DicomWatcher polling service."""

    def test_init_defaults(self):
        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=lambda c: None,
        )
        assert watcher.poll_interval == 5
        assert watcher.is_running is False
        assert watcher.orthanc_url == "http://localhost:8042"

    def test_init_custom_interval(self):
        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=lambda c: None,
            poll_interval=10,
        )
        assert watcher.poll_interval == 10

    def test_start_stop_lifecycle(self):
        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=lambda c: None,
            poll_interval=1,
        )
        # Patch _poll_changes to avoid real HTTP
        watcher._poll_changes = MagicMock()

        watcher.start()
        assert watcher.is_running is True
        assert watcher._thread is not None
        assert watcher._thread.is_alive()

        watcher.stop()
        assert watcher.is_running is False

    def test_start_twice_raises(self):
        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=lambda c: None,
            poll_interval=1,
        )
        watcher._poll_changes = MagicMock()
        watcher.start()
        with pytest.raises(RuntimeError, match="already running"):
            watcher.start()
        watcher.stop()

    def test_poll_invokes_callback_on_stable_study(self):
        """Verify that a StableStudy change triggers the callback."""
        received = []

        def on_study(change):
            received.append(change)

        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=on_study,
            poll_interval=60,
        )

        mock_response = MagicMock()
        mock_response.json.return_value = {
            "Changes": [
                {
                    "ChangeType": "StableStudy",
                    "ID": "orthanc-study-1",
                    "ResourceType": "Study",
                    "Seq": 42,
                },
                {
                    "ChangeType": "NewInstance",
                    "ID": "orthanc-inst-1",
                    "ResourceType": "Instance",
                    "Seq": 41,
                },
            ],
            "Done": True,
            "Last": 42,
        }
        mock_response.raise_for_status = MagicMock()

        with patch("src.ingest.dicom_watcher.httpx.Client") as MockClient:
            mock_client = MagicMock()
            mock_client.get.return_value = mock_response
            mock_client.__enter__ = MagicMock(return_value=mock_client)
            mock_client.__exit__ = MagicMock(return_value=False)
            MockClient.return_value = mock_client

            watcher._poll_changes()

        assert len(received) == 1
        assert received[0]["ChangeType"] == "StableStudy"
        assert received[0]["ID"] == "orthanc-study-1"

    def test_poll_advances_last_change(self):
        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=lambda c: None,
        )
        assert watcher._last_change == 0

        mock_response = MagicMock()
        mock_response.json.return_value = {
            "Changes": [],
            "Done": True,
            "Last": 99,
        }
        mock_response.raise_for_status = MagicMock()

        with patch("src.ingest.dicom_watcher.httpx.Client") as MockClient:
            mock_client = MagicMock()
            mock_client.get.return_value = mock_response
            mock_client.__enter__ = MagicMock(return_value=mock_client)
            mock_client.__exit__ = MagicMock(return_value=False)
            MockClient.return_value = mock_client

            watcher._poll_changes()

        assert watcher._last_change == 99

    def test_callback_error_does_not_crash_watcher(self):
        """A failing callback should not kill the watcher."""
        def bad_callback(change):
            raise ValueError("Intentional test error")

        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042",
            callback=bad_callback,
        )

        mock_response = MagicMock()
        mock_response.json.return_value = {
            "Changes": [
                {"ChangeType": "StableStudy", "ID": "study-x", "Seq": 10},
            ],
            "Done": True,
            "Last": 10,
        }
        mock_response.raise_for_status = MagicMock()

        with patch("src.ingest.dicom_watcher.httpx.Client") as MockClient:
            mock_client = MagicMock()
            mock_client.get.return_value = mock_response
            mock_client.__enter__ = MagicMock(return_value=mock_client)
            mock_client.__exit__ = MagicMock(return_value=False)
            MockClient.return_value = mock_client

            # Should not raise
            watcher._poll_changes()

        # Sequence should still advance
        assert watcher._last_change == 10

    def test_url_trailing_slash_stripped(self):
        watcher = DicomWatcher(
            orthanc_url="http://localhost:8042/",
            callback=lambda c: None,
        )
        assert watcher.orthanc_url == "http://localhost:8042"


# =====================================================================
# Test: Data models
# =====================================================================


class TestDataModels:
    """Tests for the Pydantic event and result models."""

    def test_dicom_study_event_defaults(self):
        event = DicomStudyEvent(event_type="study.complete")
        assert event.study_uid is None
        assert event.series_count == 0
        assert event.instance_count == 0

    def test_dicom_ingestion_result_has_timestamp(self):
        result = DicomIngestionResult(
            study_uid="1.2.3",
            patient_id="PAT-X",
            modality="CT",
        )
        assert result.processed_at is not None
        assert len(result.processed_at) > 0

    def test_dicom_ingestion_result_default_status(self):
        result = DicomIngestionResult(
            study_uid="1.2.3",
            patient_id="PAT-X",
            modality="CT",
        )
        assert result.workflow_status == "pending"
