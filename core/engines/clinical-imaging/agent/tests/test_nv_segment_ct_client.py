"""Tests for NV-Segment-CT NIM client (132 classes incl. 7 tumors).

All tests use mock mode -- no real NIM service needed.

Author: Adam Jones
Date: April 2026
"""

import pytest

from src.models import SegmentationResult
from src.nim.nv_segment_ct_client import (
    NV_SEGMENT_CT_CLASSES,
    TUMOR_CLASSES,
    NVSegmentCTClient,
)


# ===================================================================
# NVSegmentCTClient
# ===================================================================


class TestNVSegmentCTClient:
    """Tests for NVSegmentCTClient mock mode."""

    def _make_mock_client(self) -> NVSegmentCTClient:
        """Create a client forced into mock mode."""
        client = NVSegmentCTClient("http://localhost:8534", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        return client

    # ── test_segment_mock ──

    def test_segment_mock_returns_segmentation_result(self):
        client = self._make_mock_client()
        result = client.segment("/tmp/test.nii.gz")
        assert isinstance(result, SegmentationResult)
        assert result.is_mock is True

    def test_segment_mock_has_classes_detected(self):
        client = self._make_mock_client()
        result = client.segment("/tmp/test.nii.gz")
        assert len(result.classes_detected) >= 5
        assert len(result.classes_detected) <= 10

    def test_segment_mock_has_volumes(self):
        client = self._make_mock_client()
        result = client.segment("/tmp/test.nii.gz")
        assert len(result.volumes) >= 5
        for cls, vol in result.volumes.items():
            assert vol > 0.0

    def test_segment_mock_has_inference_time(self):
        client = self._make_mock_client()
        result = client.segment("/tmp/test.nii.gz")
        assert result.inference_time_ms > 0

    def test_segment_mock_model_name(self):
        client = self._make_mock_client()
        result = client.segment("/tmp/test.nii.gz")
        assert result.model_name == "nv-segment-ct-mock"

    def test_segment_with_specific_classes(self):
        client = self._make_mock_client()
        result = client.segment("/tmp/test.nii.gz", classes=["liver", "spleen"])
        assert isinstance(result, SegmentationResult)

    # ── test_segment_tumors_mock ──

    def test_segment_tumors_mock_returns_segmentation_result(self):
        client = self._make_mock_client()
        result = client.segment_tumors("/tmp/test.nii.gz")
        assert isinstance(result, SegmentationResult)
        assert result.is_mock is True

    def test_segment_tumors_mock_only_tumor_classes(self):
        client = self._make_mock_client()
        result = client.segment_tumors("/tmp/test.nii.gz")
        for cls in result.classes_detected:
            assert cls in TUMOR_CLASSES, (
                f"Non-tumor class {cls!r} found in tumor-only segmentation"
            )

    def test_segment_tumors_mock_has_small_volumes(self):
        """Tumor volumes should be realistically small (< 25 cm^3)."""
        client = self._make_mock_client()
        result = client.segment_tumors("/tmp/test.nii.gz")
        for cls, vol in result.volumes.items():
            assert vol < 25.0, (
                f"Tumor {cls!r} has unrealistically large volume: {vol} cm^3"
            )

    # ── test_segment_interactive_mock ──

    def test_segment_interactive_mock_returns_result(self):
        client = self._make_mock_client()
        prompts = [
            {"point": [100, 150, 50], "label": 1},
            {"point": [200, 250, 60], "label": 0},
        ]
        result = client.segment_interactive("/tmp/test.nii.gz", prompts)
        assert isinstance(result, SegmentationResult)
        assert result.is_mock is True

    def test_segment_interactive_mock_has_classes(self):
        client = self._make_mock_client()
        prompts = [{"point": [100, 150, 50], "label": 1}]
        result = client.segment_interactive("/tmp/test.nii.gz", prompts)
        assert len(result.classes_detected) >= 1

    # ── test_132_classes ──

    def test_132_classes_count(self):
        """Verify all 132 classes are present in the class list."""
        assert len(NV_SEGMENT_CT_CLASSES) == 132

    def test_132_classes_no_duplicates(self):
        assert len(NV_SEGMENT_CT_CLASSES) == len(set(NV_SEGMENT_CT_CLASSES))

    def test_get_supported_classes_returns_132(self):
        client = NVSegmentCTClient("http://localhost:8534")
        classes = client.get_supported_classes()
        assert isinstance(classes, list)
        assert len(classes) == 132
        assert "liver" in classes
        assert "brain" in classes
        assert "heart" in classes
        assert "liver_tumor" in classes

    # ── test_tumor_classes_subset ──

    def test_tumor_classes_count(self):
        assert len(TUMOR_CLASSES) == 7

    def test_tumor_classes_are_subset_of_all(self):
        """All 7 tumor classes must be present in the full 132-class list."""
        all_classes_set = set(NV_SEGMENT_CT_CLASSES)
        for tumor_cls in TUMOR_CLASSES:
            assert tumor_cls in all_classes_set, (
                f"Tumor class {tumor_cls!r} missing from NV_SEGMENT_CT_CLASSES"
            )

    def test_tumor_classes_expected_names(self):
        expected = {
            "liver_tumor", "lung_tumor", "kidney_tumor",
            "pancreas_tumor", "colon_tumor", "hepatic_vessel_tumor",
            "adrenal_gland_tumor",
        }
        assert set(TUMOR_CLASSES) == expected

    # ── test_service_status ──

    def test_service_status_mock_when_unavailable(self):
        """get_status() returns 'mock' when service is unavailable."""
        client = NVSegmentCTClient("http://localhost:8534", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        assert client.get_status() == "mock"

    def test_service_status_unavailable_when_mock_disabled(self):
        client = NVSegmentCTClient("http://localhost:8534", mock_enabled=False)
        client._available = False
        client._last_check = 1e18
        assert client.get_status() == "unavailable"

    def test_service_name(self):
        client = NVSegmentCTClient("http://localhost:8534")
        assert client.service_name == "nv-segment-ct"

    def test_init_default_mock_enabled(self):
        client = NVSegmentCTClient("http://localhost:8534")
        assert client.mock_enabled is True
