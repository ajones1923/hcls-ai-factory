"""Tests for NV-Generate-CT and NV-Generate-MR NIM clients.

All tests use mock mode -- no real NIM service needed.

Author: Adam Jones
Date: April 2026
"""

import json

import pytest

from src.nim.nv_generate_ct_client import (
    NVGenerateCTClient,
    SyntheticGenerationConfig,
    SyntheticVolumeResult,
)
from src.nim.nv_generate_mr_client import (
    MRContrastType,
    NVGenerateMRClient,
    SyntheticMRConfig,
    SyntheticMRResult,
)


# ===================================================================
# Helpers
# ===================================================================


def _make_ct_client() -> NVGenerateCTClient:
    """Create a CT client forced into mock mode."""
    client = NVGenerateCTClient("http://localhost:8539", mock_enabled=True)
    client._available = False
    client._last_check = 1e18
    return client


def _make_mr_client() -> NVGenerateMRClient:
    """Create an MR client forced into mock mode."""
    client = NVGenerateMRClient("http://localhost:8540", mock_enabled=True)
    client._available = False
    client._last_check = 1e18
    return client


# ===================================================================
# NV-Generate-CT Tests
# ===================================================================


class TestNVGenerateCTClient:
    """Tests for NVGenerateCTClient mock mode."""

    # ── test_ct_generate_mock ──

    def test_ct_generate_mock(self):
        """Default generation returns SyntheticVolumeResult with correct fields."""
        client = _make_ct_client()
        result = client.generate()
        assert isinstance(result, SyntheticVolumeResult)
        assert result.is_mock is True
        assert result.volume_path is not None
        assert result.mask_path is not None
        assert result.body_region == "chest"
        assert result.resolution_mm == 1.0
        assert result.volume_shape == [256, 256, 256]
        assert result.classes_generated == 127
        assert result.generation_time_ms >= 2000.0
        assert result.generation_time_ms <= 15000.0
        assert "nv-generate-ct" in result.model_name

    # ── test_ct_generate_custom_config ──

    def test_ct_generate_custom_config(self):
        """Custom body_region, resolution, and volume_size are respected."""
        client = _make_ct_client()
        config = SyntheticGenerationConfig(
            body_region="abdomen",
            resolution_mm=0.8,
            volume_size=[512, 512, 512],
            num_classes=50,
        )
        result = client.generate(config)
        assert isinstance(result, SyntheticVolumeResult)
        assert result.body_region == "abdomen"
        assert result.resolution_mm == 0.8
        assert result.volume_shape == [512, 512, 512]
        assert result.classes_generated == 50

    # ── test_ct_generate_with_pathology ──

    def test_ct_generate_with_pathology(self):
        """include_pathology=True with pathology_type works."""
        client = _make_ct_client()
        config = SyntheticGenerationConfig(
            body_region="chest",
            include_pathology=True,
            pathology_type="nodule",
        )
        result = client.generate(config)
        assert isinstance(result, SyntheticVolumeResult)
        assert result.is_mock is True
        assert "nodule" in result.volume_path

    # ── test_ct_batch_generation ──

    def test_ct_batch_generation(self):
        """generate_training_batch returns multiple results."""
        client = _make_ct_client()
        config = SyntheticGenerationConfig(body_region="head")
        results = client.generate_training_batch(config, count=3)
        assert isinstance(results, list)
        assert len(results) == 3
        for r in results:
            assert isinstance(r, SyntheticVolumeResult)
            assert r.is_mock is True
            assert r.body_region == "head"

    # ── test_ct_seed_reproducibility ──

    def test_ct_seed_reproducibility(self):
        """Same seed produces the same mock result."""
        client = _make_ct_client()
        config_a = SyntheticGenerationConfig(seed=42)
        config_b = SyntheticGenerationConfig(seed=42)
        result_a = client.generate(config_a)
        result_b = client.generate(config_b)
        assert result_a.volume_path == result_b.volume_path
        assert result_a.mask_path == result_b.mask_path
        assert result_a.generation_time_ms == result_b.generation_time_ms

    # ── test_ct_different_seeds ──

    def test_ct_different_seeds(self):
        """Different seeds produce different mock results."""
        client = _make_ct_client()
        result_a = client.generate(SyntheticGenerationConfig(seed=42))
        result_b = client.generate(SyntheticGenerationConfig(seed=99))
        assert result_a.volume_path != result_b.volume_path

    # ── test_ct_client_status ──

    def test_ct_client_status(self):
        """get_status() returns 'mock' when service is unavailable."""
        client = _make_ct_client()
        assert client.get_status() == "mock"

    def test_ct_client_status_unavailable_when_mock_disabled(self):
        client = NVGenerateCTClient("http://localhost:8539", mock_enabled=False)
        client._available = False
        client._last_check = 1e18
        assert client.get_status() == "unavailable"

    def test_ct_service_name(self):
        client = NVGenerateCTClient("http://localhost:8539")
        assert client.service_name == "nv-generate-ct"

    def test_ct_init_default_mock_enabled(self):
        client = NVGenerateCTClient("http://localhost:8539")
        assert client.mock_enabled is True

    # ── test_ct_volume_size_affects_timing ──

    def test_ct_large_volume_slower(self):
        """Larger volumes should produce longer simulated generation times."""
        client = _make_ct_client()
        small = client.generate(SyntheticGenerationConfig(
            volume_size=[64, 64, 64], seed=100,
        ))
        large = client.generate(SyntheticGenerationConfig(
            volume_size=[512, 512, 768], seed=100,
        ))
        # Large volume should have longer generation time (before clamping)
        # Both are clamped to [2000, 15000] range but large should be >= small
        assert large.generation_time_ms >= small.generation_time_ms

    # ── test_ct_batch_unique_seeds ──

    def test_ct_batch_unique_seeds(self):
        """Batch generation with a base seed produces unique volumes."""
        client = _make_ct_client()
        config = SyntheticGenerationConfig(seed=1000)
        results = client.generate_training_batch(config, count=3)
        paths = [r.volume_path for r in results]
        assert len(set(paths)) == 3, "Each volume in batch should have a unique path"


# ===================================================================
# NV-Generate-MR Tests
# ===================================================================


class TestNVGenerateMRClient:
    """Tests for NVGenerateMRClient mock mode."""

    # ── test_mr_generate_mock ──

    def test_mr_generate_mock(self):
        """Default generation returns SyntheticMRResult."""
        client = _make_mr_client()
        result = client.generate()
        assert isinstance(result, SyntheticMRResult)
        assert result.is_mock is True
        assert result.volume_path is not None
        assert result.mask_path is not None
        assert result.body_region == "brain"
        assert result.contrast == "t1w"
        assert result.resolution_mm == 1.0
        assert result.volume_shape == [256, 256, 256]
        assert result.generation_time_ms >= 2000.0
        assert result.generation_time_ms <= 15000.0
        assert "nv-generate-mr" in result.model_name

    # ── test_mr_t1w_contrast ──

    def test_mr_t1w_contrast(self):
        """T1W contrast type is correctly reflected in result."""
        client = _make_mr_client()
        config = SyntheticMRConfig(contrast=MRContrastType.T1W)
        result = client.generate(config)
        assert result.contrast == "t1w"
        assert "t1w" in result.volume_path

    # ── test_mr_flair_contrast ──

    def test_mr_flair_contrast(self):
        """FLAIR contrast type is correctly reflected in result."""
        client = _make_mr_client()
        config = SyntheticMRConfig(contrast=MRContrastType.FLAIR)
        result = client.generate(config)
        assert result.contrast == "flair"
        assert "flair" in result.volume_path

    # ── test_mr_multi_contrast ──

    def test_mr_multi_contrast(self):
        """generate_multi_contrast produces one volume per contrast."""
        client = _make_mr_client()
        contrasts = [MRContrastType.T1W, MRContrastType.T2W, MRContrastType.FLAIR]
        results = client.generate_multi_contrast(
            body_region="brain",
            contrasts=contrasts,
        )
        assert isinstance(results, list)
        assert len(results) == 3
        contrast_values = [r.contrast for r in results]
        assert "t1w" in contrast_values
        assert "t2w" in contrast_values
        assert "flair" in contrast_values
        for r in results:
            assert isinstance(r, SyntheticMRResult)
            assert r.is_mock is True
            assert r.body_region == "brain"

    # ── test_mr_brain_skull_stripped ──

    def test_mr_brain_skull_stripped(self):
        """skull_stripped option is reflected in the output path."""
        client = _make_mr_client()
        config = SyntheticMRConfig(
            body_region="brain",
            skull_stripped=True,
        )
        result = client.generate(config)
        assert isinstance(result, SyntheticMRResult)
        assert "stripped" in result.volume_path

    # ── test_mr_with_pathology ──

    def test_mr_with_pathology_tumor(self):
        """Tumor pathology injection works."""
        client = _make_mr_client()
        config = SyntheticMRConfig(
            body_region="brain",
            include_pathology=True,
            pathology_type="tumor",
        )
        result = client.generate(config)
        assert result.is_mock is True
        assert "tumor" in result.volume_path

    def test_mr_with_pathology_lesion(self):
        """Lesion pathology injection works."""
        client = _make_mr_client()
        config = SyntheticMRConfig(
            body_region="brain",
            include_pathology=True,
            pathology_type="lesion",
        )
        result = client.generate(config)
        assert "lesion" in result.volume_path

    def test_mr_with_pathology_infarct(self):
        """Infarct pathology injection works."""
        client = _make_mr_client()
        config = SyntheticMRConfig(
            body_region="brain",
            include_pathology=True,
            pathology_type="infarct",
        )
        result = client.generate(config)
        assert "infarct" in result.volume_path

    # ── test_mr_client_status ──

    def test_mr_client_status(self):
        """get_status() returns 'mock' when service is unavailable."""
        client = _make_mr_client()
        assert client.get_status() == "mock"

    def test_mr_client_status_unavailable_when_mock_disabled(self):
        client = NVGenerateMRClient("http://localhost:8540", mock_enabled=False)
        client._available = False
        client._last_check = 1e18
        assert client.get_status() == "unavailable"

    def test_mr_service_name(self):
        client = NVGenerateMRClient("http://localhost:8540")
        assert client.service_name == "nv-generate-mr"

    def test_mr_init_default_mock_enabled(self):
        client = NVGenerateMRClient("http://localhost:8540")
        assert client.mock_enabled is True

    # ── test_mr_dwi_contrast ──

    def test_mr_dwi_contrast(self):
        """DWI contrast type works."""
        client = _make_mr_client()
        config = SyntheticMRConfig(contrast=MRContrastType.DWI)
        result = client.generate(config)
        assert result.contrast == "dwi"

    # ── test_mr_swi_contrast ──

    def test_mr_swi_contrast(self):
        """SWI contrast type works."""
        client = _make_mr_client()
        config = SyntheticMRConfig(contrast=MRContrastType.SWI)
        result = client.generate(config)
        assert result.contrast == "swi"

    # ── test_mr_prostate_region ──

    def test_mr_prostate_region(self):
        """Prostate body region generates correctly."""
        client = _make_mr_client()
        config = SyntheticMRConfig(body_region="prostate", contrast=MRContrastType.T2W)
        result = client.generate(config)
        assert result.body_region == "prostate"
        assert result.contrast == "t2w"

    # ── test_mr_multi_contrast_defaults ──

    def test_mr_multi_contrast_defaults(self):
        """Default multi-contrast generates T1W, T2W, FLAIR."""
        client = _make_mr_client()
        results = client.generate_multi_contrast()
        assert len(results) == 3
        contrast_values = sorted([r.contrast for r in results])
        assert contrast_values == ["flair", "t1w", "t2w"]


# ===================================================================
# Model Serialization Tests
# ===================================================================


class TestModelsSerialization:
    """Verify all Pydantic models serialize to JSON."""

    def test_synthetic_generation_config_serializable(self):
        config = SyntheticGenerationConfig(
            body_region="chest",
            resolution_mm=1.5,
            volume_size=[512, 512, 256],
            include_pathology=True,
            pathology_type="nodule",
            seed=42,
        )
        data = json.loads(config.model_dump_json())
        assert data["body_region"] == "chest"
        assert data["seed"] == 42

    def test_synthetic_volume_result_serializable(self):
        result = SyntheticVolumeResult(
            volume_path="/tmp/test.nii.gz",
            mask_path="/tmp/mask.nii.gz",
            body_region="chest",
            resolution_mm=1.0,
            volume_shape=[256, 256, 256],
            classes_generated=127,
            generation_time_ms=5000.0,
            model_name="nv-generate-ct-mock",
            is_mock=True,
        )
        data = json.loads(result.model_dump_json())
        assert data["body_region"] == "chest"
        assert data["is_mock"] is True

    def test_synthetic_mr_config_serializable(self):
        config = SyntheticMRConfig(
            body_region="brain",
            contrast=MRContrastType.FLAIR,
            skull_stripped=True,
            seed=99,
        )
        data = json.loads(config.model_dump_json())
        assert data["contrast"] == "flair"
        assert data["skull_stripped"] is True

    def test_synthetic_mr_result_serializable(self):
        result = SyntheticMRResult(
            volume_path="/tmp/test_mr.nii.gz",
            mask_path="/tmp/mask_mr.nii.gz",
            body_region="brain",
            contrast="t1w",
            resolution_mm=1.0,
            volume_shape=[256, 256, 256],
            generation_time_ms=6000.0,
            model_name="nv-generate-mr-mock",
            is_mock=True,
        )
        data = json.loads(result.model_dump_json())
        assert data["contrast"] == "t1w"
        assert data["is_mock"] is True

    def test_mr_contrast_type_enum_serializable(self):
        """MRContrastType enum values serialize correctly."""
        for contrast in MRContrastType:
            assert isinstance(contrast.value, str)
            assert len(contrast.value) > 0
