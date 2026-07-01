"""Tests for MONAI Label interactive annotation integration.

All tests work without the ``monailabel`` package installed or a
MONAI Label server running.

Author: Adam Jones
Date: April 2026
"""

import time
from unittest.mock import patch

import pytest

from src.annotation import AnnotationConfig, AnnotationSession, MONAILabelManager


# ═══════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def manager():
    """Return a MONAILabelManager with default config."""
    return MONAILabelManager()


@pytest.fixture
def manager_with_session(manager):
    """Return a manager that already has one active session."""
    session = manager.start_session(
        study_id="study-001",
        annotator="dr_smith",
        model="vista3d_interactive",
    )
    return manager, session


# ═══════════════════════════════════════════════════════════════════════
# Graceful Fallback
# ═══════════════════════════════════════════════════════════════════════


class TestGracefulFallback:
    """MONAI Label manager works without the monailabel package."""

    def test_monailabel_graceful_fallback(self, manager):
        """Manager initializes even when monailabel is not installed."""
        assert manager is not None
        assert isinstance(manager.config, AnnotationConfig)
        # _monailabel_available will be False in test env (no monailabel)
        # but the manager should still be fully functional
        assert isinstance(manager._monailabel_available, bool)


# ═══════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════


class TestConfigDefaults:
    """AnnotationConfig has correct default values."""

    def test_config_defaults(self):
        cfg = AnnotationConfig()
        assert cfg.host == "0.0.0.0"
        assert cfg.port == 8527
        assert cfg.app_dir == "data/monai_label_apps"
        assert cfg.studies_dir == "data/monai_label_studies"
        assert cfg.datastore_type == "orthanc"
        assert cfg.orthanc_url == "http://localhost:8042"
        assert "vista3d_interactive" in cfg.models
        assert "segmentation_spleen" in cfg.models
        assert "deepedit" in cfg.models
        assert cfg.auto_segmentation is True
        assert cfg.active_learning is True


# ═══════════════════════════════════════════════════════════════════════
# Session Lifecycle
# ═══════════════════════════════════════════════════════════════════════


class TestSessionManagement:
    """Session creation, tracking, and completion."""

    def test_start_session(self, manager):
        """start_session creates a session with a UUID and active status."""
        session = manager.start_session(
            study_id="study-abc",
            annotator="dr_jones",
            model="deepedit",
        )
        assert session.session_id != ""
        assert len(session.session_id) == 36  # UUID4 format
        assert session.study_id == "study-abc"
        assert session.annotator == "dr_jones"
        assert session.model_used == "deepedit"
        assert session.status == "active"
        assert session.labels_created == 0
        assert session.labels_corrected == 0
        assert session.created_at != ""

    def test_end_session(self, manager_with_session):
        """end_session marks the session completed and records time."""
        mgr, session = manager_with_session
        ended = mgr.end_session(session.session_id, status="completed")
        assert ended.status == "completed"
        assert ended.time_seconds >= 0.0
        # Session should no longer be in active dict
        assert session.session_id not in mgr._sessions
        # But should be in completed list
        assert ended in mgr._completed_sessions

    def test_end_session_abandoned(self, manager_with_session):
        """end_session accepts 'abandoned' status."""
        mgr, session = manager_with_session
        ended = mgr.end_session(session.session_id, status="abandoned")
        assert ended.status == "abandoned"

    def test_end_session_not_found(self, manager):
        """end_session raises KeyError for unknown session."""
        with pytest.raises(KeyError):
            manager.end_session("nonexistent-session-id")

    def test_session_not_found(self, manager):
        """get_session returns None for unknown session ID."""
        result = manager.get_session("does-not-exist")
        assert result is None

    def test_get_session_active(self, manager_with_session):
        """get_session finds active sessions."""
        mgr, session = manager_with_session
        found = mgr.get_session(session.session_id)
        assert found is not None
        assert found.session_id == session.session_id

    def test_get_session_completed(self, manager_with_session):
        """get_session finds completed sessions."""
        mgr, session = manager_with_session
        mgr.end_session(session.session_id, status="completed")
        found = mgr.get_session(session.session_id)
        assert found is not None
        assert found.status == "completed"


# ═══════════════════════════════════════════════════════════════════════
# Label Tracking
# ═══════════════════════════════════════════════════════════════════════


class TestLabelTracking:
    """Label creation and correction recording."""

    def test_record_label_created(self, manager_with_session):
        """record_label increments labels_created."""
        mgr, session = manager_with_session
        assert mgr.record_label(session.session_id, "created") is True
        assert mgr.record_label(session.session_id, "created") is True
        assert mgr.record_label(session.session_id, "created") is True
        found = mgr.get_session(session.session_id)
        assert found.labels_created == 3

    def test_record_label_corrected(self, manager_with_session):
        """record_label increments labels_corrected."""
        mgr, session = manager_with_session
        assert mgr.record_label(session.session_id, "corrected") is True
        assert mgr.record_label(session.session_id, "corrected") is True
        found = mgr.get_session(session.session_id)
        assert found.labels_corrected == 2

    def test_record_label_invalid_session(self, manager):
        """record_label returns False for unknown session."""
        assert manager.record_label("bogus-id", "created") is False

    def test_record_label_invalid_type(self, manager_with_session):
        """record_label returns False for invalid label_type."""
        mgr, session = manager_with_session
        assert mgr.record_label(session.session_id, "deleted") is False


# ═══════════════════════════════════════════════════════════════════════
# Statistics
# ═══════════════════════════════════════════════════════════════════════


class TestStats:
    """Aggregate annotation statistics."""

    def test_get_stats_empty(self, manager):
        """Empty manager returns zero stats."""
        stats = manager.get_stats()
        assert stats.total_sessions == 0
        assert stats.total_labels == 0
        assert stats.total_corrections == 0
        assert stats.correction_rate == 0.0
        assert stats.avg_session_time_seconds == 0.0
        assert stats.models_used == {}
        assert stats.annotators == {}
        assert stats.studies_annotated == 0
        assert stats.labels_for_flare == 0

    def test_get_stats_after_sessions(self, manager):
        """Stats reflect data from multiple completed sessions."""
        # Session 1: dr_smith, vista3d, 5 labels, 1 correction
        s1 = manager.start_session("study-001", "dr_smith", "vista3d_interactive")
        for _ in range(5):
            manager.record_label(s1.session_id, "created")
        manager.record_label(s1.session_id, "corrected")
        manager.end_session(s1.session_id)

        # Session 2: dr_jones, deepedit, 3 labels, 2 corrections
        s2 = manager.start_session("study-002", "dr_jones", "deepedit")
        for _ in range(3):
            manager.record_label(s2.session_id, "created")
        for _ in range(2):
            manager.record_label(s2.session_id, "corrected")
        manager.end_session(s2.session_id)

        stats = manager.get_stats()
        assert stats.total_sessions == 2
        assert stats.total_labels == 8  # 5 + 3
        assert stats.total_corrections == 3  # 1 + 2
        assert stats.studies_annotated == 2
        assert stats.models_used == {"vista3d_interactive": 1, "deepedit": 1}
        assert stats.annotators == {"dr_smith": 1, "dr_jones": 1}
        assert stats.labels_for_flare == 11  # all labels from completed sessions

    def test_correction_rate(self, manager):
        """Correction rate = corrections / (labels + corrections)."""
        s = manager.start_session("study-001", "dr_test", "deepedit")
        # 8 created, 2 corrected -> rate = 2/10 = 0.2
        for _ in range(8):
            manager.record_label(s.session_id, "created")
        for _ in range(2):
            manager.record_label(s.session_id, "corrected")
        manager.end_session(s.session_id)

        stats = manager.get_stats()
        assert stats.correction_rate == pytest.approx(0.2, abs=0.001)


# ═══════════════════════════════════════════════════════════════════════
# Model Catalog
# ═══════════════════════════════════════════════════════════════════════


class TestModels:
    """Available model listing."""

    def test_available_models(self, manager):
        """get_available_models returns 4+ preset models when server is down."""
        models = manager.get_available_models()
        assert len(models) >= 4
        names = [m["name"] for m in models]
        assert "vista3d_interactive" in names
        assert "segmentation_spleen" in names
        assert "deepedit" in names
        assert "segmentation_liver" in names
        # Each model has required keys
        for m in models:
            assert "name" in m
            assert "type" in m
            assert "description" in m
            assert "source" in m


# ═══════════════════════════════════════════════════════════════════════
# FLARE Export
# ═══════════════════════════════════════════════════════════════════════


class TestFLAREExport:
    """Export annotated data for NVIDIA FLARE federated training."""

    def test_export_for_flare(self, manager, tmp_path):
        """export_for_flare returns counts and writes manifest."""
        # Create and complete two sessions with labels
        s1 = manager.start_session("study-001", "dr_a", "vista3d_interactive")
        for _ in range(4):
            manager.record_label(s1.session_id, "created")
        manager.record_label(s1.session_id, "corrected")
        manager.end_session(s1.session_id)

        s2 = manager.start_session("study-002", "dr_b", "deepedit")
        for _ in range(3):
            manager.record_label(s2.session_id, "created")
        manager.end_session(s2.session_id)

        out_dir = str(tmp_path / "flare_out")
        result = manager.export_for_flare(output_dir=out_dir)

        assert result["total_labels"] == 8  # 4+1 + 3
        assert "vista3d_interactive" in result["by_model"]
        assert "deepedit" in result["by_model"]
        assert result["by_model"]["vista3d_interactive"] == 5
        assert result["by_model"]["deepedit"] == 3
        assert "manifest_path" in result

        # Manifest file should exist
        import json
        from pathlib import Path

        manifest = json.loads(Path(result["manifest_path"]).read_text())
        assert manifest["total_labels"] == 8
        assert manifest["total_sessions"] == 2

    def test_export_for_flare_empty(self, manager, tmp_path):
        """export_for_flare with no sessions returns zero counts."""
        result = manager.export_for_flare(output_dir=str(tmp_path / "empty"))
        assert result["total_labels"] == 0
        assert result["by_model"] == {}


# ═══════════════════════════════════════════════════════════════════════
# OHIF Integration
# ═══════════════════════════════════════════════════════════════════════


class TestOHIFIntegration:
    """OHIF Viewer URL generation."""

    def test_ohif_annotation_url(self, manager):
        """get_ohif_annotation_url generates valid URL with study ID."""
        url = manager.get_ohif_annotation_url("1.2.3.4.5.6.789")
        assert "1.2.3.4.5.6.789" in url
        assert "StudyInstanceUIDs=" in url
        assert "8526" in url
        assert "monaiLabel" in url


# ═══════════════════════════════════════════════════════════════════════
# Server Status
# ═══════════════════════════════════════════════════════════════════════


class TestServerStatus:
    """Server availability checks."""

    def test_server_status_when_unavailable(self, manager):
        """get_server_status returns unavailable when server is down."""
        status = manager.get_server_status()
        assert status["available"] is False
        assert status["version"] is None
        assert "message" in status
