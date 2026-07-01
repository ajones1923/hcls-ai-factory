"""Tests for NVIDIA DALI and cuCIM GPU-accelerated preprocessing.

All tests run WITHOUT DALI or cuCIM installed -- they exercise the CPU
fallback paths using numpy/scipy only.  Synthetic test data is generated
with numpy so no external data files are required.

Author: Adam Jones
Date: April 2026
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.preprocessing.dali_pipeline import DALIPreprocessor
from src.preprocessing.cucim_processor import CuCIMProcessor


# ═══════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════


@pytest.fixture
def dali():
    """Return a DALIPreprocessor instance (CPU fallback expected)."""
    return DALIPreprocessor(device="gpu")


@pytest.fixture
def cucim():
    """Return a CuCIMProcessor instance (CPU fallback expected)."""
    return CuCIMProcessor()


@pytest.fixture
def ct_volume():
    """Synthetic CT volume with HU-like values: shape (32, 64, 64)."""
    rng = np.random.RandomState(42)
    # Simulate HU range roughly -1024 to +3071
    return rng.uniform(-1024, 3071, size=(32, 64, 64)).astype(np.float32)


@pytest.fixture
def mri_volume():
    """Synthetic MRI volume with arbitrary intensity: shape (32, 64, 64)."""
    rng = np.random.RandomState(99)
    return rng.uniform(0, 4096, size=(32, 64, 64)).astype(np.float32)


@pytest.fixture
def cxr_image():
    """Synthetic 2-D chest X-ray grayscale image: shape (512, 512)."""
    rng = np.random.RandomState(7)
    return rng.uniform(0, 255, size=(512, 512)).astype(np.float32)


@pytest.fixture
def binary_mask_3d():
    """Synthetic 3-D binary mask with two blobs and some noise."""
    mask = np.zeros((32, 64, 64), dtype=bool)
    # Blob 1: large sphere at center
    z, y, x = np.ogrid[:32, :64, :64]
    blob1 = ((z - 16) ** 2 + (y - 32) ** 2 + (x - 32) ** 2) < 10 ** 2
    mask[blob1] = True
    # Blob 2: small sphere offset
    blob2 = ((z - 8) ** 2 + (y - 50) ** 2 + (x - 50) ** 2) < 3 ** 2
    mask[blob2] = True
    # Scatter noise: 20 random isolated voxels
    rng = np.random.RandomState(11)
    for _ in range(20):
        mask[rng.randint(0, 32), rng.randint(0, 64), rng.randint(0, 64)] = True
    return mask


@pytest.fixture
def binary_mask_with_hole():
    """3-D binary mask (sphere shell) with an interior hole."""
    mask = np.zeros((32, 64, 64), dtype=bool)
    z, y, x = np.ogrid[:32, :64, :64]
    outer = ((z - 16) ** 2 + (y - 32) ** 2 + (x - 32) ** 2) < 12 ** 2
    inner = ((z - 16) ** 2 + (y - 32) ** 2 + (x - 32) ** 2) < 5 ** 2
    mask[outer & ~inner] = True
    return mask


# ═══════════════════════════════════════════════════════════════════
# DALI PREPROCESSOR TESTS
# ═══════════════════════════════════════════════════════════════════


class TestDALIPreprocessor:
    """Tests for DALIPreprocessor (CPU fallback path)."""

    def test_dali_fallback_when_unavailable(self, dali):
        """DALIPreprocessor initialises and reports CPU mode when DALI is missing."""
        # On CI / dev without DALI, GPU acceleration should be off
        # (This test passes either way -- it just validates the object is usable)
        assert isinstance(dali, DALIPreprocessor)
        assert isinstance(dali.is_gpu_accelerated, bool)

    def test_ct_windowing_brain(self, dali, ct_volume):
        """Brain window (center=40, width=80) clips to 0-80 HU range."""
        result = dali.apply_window(ct_volume, "brain")
        assert result.dtype == np.float32
        assert result.shape == ct_volume.shape
        assert result.min() >= 0.0 - 1e-6
        assert result.max() <= 1.0 + 1e-6

    def test_ct_windowing_lung(self, dali, ct_volume):
        """Lung window (center=-600, width=1500) applied correctly."""
        result = dali.apply_window(ct_volume, "lung")
        assert result.dtype == np.float32
        assert result.shape == ct_volume.shape
        assert result.min() >= 0.0 - 1e-6
        assert result.max() <= 1.0 + 1e-6

    def test_all_window_presets(self, dali):
        """All 8 CT window presets have valid center and width."""
        assert len(dali.CT_WINDOWS) == 8
        for name, preset in dali.CT_WINDOWS.items():
            assert "center" in preset, f"{name} missing center"
            assert "width" in preset, f"{name} missing width"
            assert preset["width"] > 0, f"{name} has non-positive width"

    def test_invalid_window_raises(self, dali, ct_volume):
        """Unknown window name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown CT window"):
            dali.apply_window(ct_volume, "nonexistent_window")

    def test_normalize_zscore(self, dali, ct_volume):
        """Z-score normalization produces mean~0, std~1."""
        result = dali.normalize(ct_volume, method="zscore")
        assert result.dtype == np.float32
        assert abs(result.mean()) < 1e-5
        assert abs(result.std() - 1.0) < 1e-5

    def test_normalize_minmax(self, dali, ct_volume):
        """Min-max normalization produces values in [0, 1]."""
        result = dali.normalize(ct_volume, method="minmax")
        assert result.dtype == np.float32
        assert result.min() >= 0.0 - 1e-6
        assert result.max() <= 1.0 + 1e-6

    def test_normalize_percentile(self, dali, ct_volume):
        """Percentile normalization clips and scales to [0, 1]."""
        result = dali.normalize(ct_volume, method="percentile")
        assert result.dtype == np.float32
        assert result.min() >= 0.0 - 1e-6
        assert result.max() <= 1.0 + 1e-6

    def test_normalize_invalid_method(self, dali, ct_volume):
        """Invalid normalization method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown normalization"):
            dali.normalize(ct_volume, method="invalid")

    def test_resize_3d(self, dali, ct_volume):
        """3-D volume resized to specified target shape."""
        target = (16, 32, 32)
        result = dali.resize(ct_volume, target)
        assert result.shape == target
        assert result.dtype == np.float32

    def test_resize_2d(self, dali, cxr_image):
        """2-D image (CXR) resized to 224x224."""
        target = (224, 224)
        result = dali.resize(cxr_image, target)
        assert result.shape == target
        assert result.dtype == np.float32

    def test_resize_noop(self, dali, cxr_image):
        """Resize with matching shape is a no-op (returns same shape)."""
        result = dali.resize(cxr_image, cxr_image.shape)
        assert result.shape == cxr_image.shape

    def test_resample_isotropic(self, dali):
        """Resampling with known zoom factor produces expected shape."""
        vol = np.random.randn(20, 40, 40).astype(np.float32)
        # spacing (2, 1, 1) -> target 1mm -> zoom (2, 1, 1)
        result = dali.resample_isotropic(vol, current_spacing=(2.0, 1.0, 1.0), target_spacing=1.0)
        # z-axis should approximately double
        assert result.shape[0] == pytest.approx(40, abs=1)
        assert result.shape[1] == pytest.approx(40, abs=1)
        assert result.shape[2] == pytest.approx(40, abs=1)

    def test_preprocess_ct_pipeline(self, dali, ct_volume):
        """Full CT pipeline runs end-to-end without error."""
        result = dali.preprocess_ct(
            ct_volume,
            window="brain",
            target_spacing=1.0,
            current_spacing=(2.0, 1.0, 1.0),
            target_shape=(32, 64, 64),
        )
        assert result.shape == (32, 64, 64)
        assert result.dtype == np.float32

    def test_preprocess_mri_pipeline(self, dali, mri_volume):
        """Full MRI pipeline runs end-to-end without error."""
        result = dali.preprocess_mri(
            mri_volume,
            target_spacing=1.0,
            current_spacing=(1.5, 0.5, 0.5),
            target_shape=(32, 64, 64),
        )
        assert result.shape == (32, 64, 64)
        assert result.dtype == np.float32

    def test_preprocess_cxr_pipeline(self, dali, cxr_image):
        """CXR pipeline runs end-to-end, produces (224, 224) output."""
        result = dali.preprocess_cxr(cxr_image, target_size=(224, 224))
        assert result.shape == (224, 224)
        assert result.dtype == np.float32
        assert result.min() >= 0.0 - 1e-6
        assert result.max() <= 1.0 + 1e-6


# ═══════════════════════════════════════════════════════════════════
# CUCIM PROCESSOR TESTS
# ═══════════════════════════════════════════════════════════════════


class TestCuCIMProcessor:
    """Tests for CuCIMProcessor (CPU fallback path)."""

    def test_cucim_fallback(self, cucim):
        """CuCIMProcessor initialises and reports CPU mode when cuCIM missing."""
        assert isinstance(cucim, CuCIMProcessor)
        assert isinstance(cucim.is_gpu_accelerated, bool)

    def test_remove_small_objects(self, cucim, binary_mask_3d):
        """Small components removed, large blob preserved."""
        # The large sphere has ~4000 voxels; noise dots are 1 voxel each
        result = cucim.remove_small_objects(binary_mask_3d, min_size=50)
        assert result.dtype == bool
        # Large blob should survive
        assert result.sum() > 100
        # Isolated noise voxels should be gone, so fewer True voxels
        assert result.sum() < binary_mask_3d.sum()

    def test_binary_closing(self, cucim):
        """Closing fills small gaps in a mask."""
        # Create a 3-D slab with a 1-voxel gap punched through it
        mask = np.zeros((16, 16, 16), dtype=bool)
        # Solid 4x4x4 block
        mask[6:10, 6:10, 6:10] = True
        # Punch a 1-voxel hole in the middle
        mask[8, 8, 8] = False
        result = cucim.binary_closing(mask, radius=1)
        assert result.dtype == bool
        # After closing, the hole should be filled
        assert result[8, 8, 8]

    def test_binary_opening(self, cucim):
        """Opening removes thin protrusions."""
        mask = np.zeros((16, 16, 16), dtype=bool)
        # Solid 6x6x6 block
        mask[5:11, 5:11, 5:11] = True
        # 1-voxel protrusion
        mask[5, 5, 4] = True
        result = cucim.binary_opening(mask, radius=1)
        assert result.dtype == bool
        # Protrusion should be removed
        assert not result[5, 5, 4]
        # Core block should largely survive
        assert result[8, 8, 8]

    def test_fill_holes(self, cucim, binary_mask_with_hole):
        """Interior holes are filled."""
        result = cucim.fill_holes(binary_mask_with_hole)
        assert result.dtype == bool
        # Filled mask should have more True voxels than the shell
        assert result.sum() > binary_mask_with_hole.sum()
        # Center (which was the hole) should now be True
        assert result[16, 32, 32]

    def test_label_components(self, cucim):
        """Correct number of connected components detected."""
        mask = np.zeros((16, 32, 32), dtype=bool)
        # Component 1
        mask[2:5, 2:5, 2:5] = True
        # Component 2 (well separated)
        mask[10:13, 20:25, 20:25] = True
        # Component 3
        mask[2:4, 25:28, 25:28] = True
        labeled, num = cucim.label_connected_components(mask)
        assert num == 3
        assert labeled.shape == mask.shape
        assert labeled.max() == 3
        # Background is 0
        assert labeled[0, 0, 0] == 0

    def test_region_properties(self, cucim):
        """Area and centroid returned for each region."""
        mask = np.zeros((16, 32, 32), dtype=bool)
        mask[2:6, 2:6, 2:6] = True   # 4x4x4 = 64 voxels
        mask[10:12, 20:24, 20:24] = True  # 2x4x4 = 32 voxels
        labeled, num = cucim.label_connected_components(mask)
        props = cucim.region_properties(labeled)
        assert len(props) == 2
        areas = sorted([p["area"] for p in props])
        assert areas == [32, 64]
        # Each prop should have centroid and bbox
        for p in props:
            assert "centroid" in p
            assert "bbox" in p
            assert len(p["centroid"]) == 3

    def test_region_properties_with_intensity(self, cucim):
        """Mean intensity computed when intensity image provided."""
        mask = np.zeros((16, 32, 32), dtype=bool)
        mask[4:8, 4:8, 4:8] = True
        labeled, _ = cucim.label_connected_components(mask)
        intensity = np.ones((16, 32, 32), dtype=np.float32) * 100.0
        props = cucim.region_properties(labeled, intensity_image=intensity)
        assert len(props) == 1
        assert "mean_intensity" in props[0]
        assert abs(props[0]["mean_intensity"] - 100.0) < 1e-3

    def test_gaussian_filter(self, cucim):
        """Gaussian smoothing produces less noisy output."""
        rng = np.random.RandomState(42)
        noisy = rng.randn(32, 32).astype(np.float32) * 10.0
        smoothed = cucim.gaussian_filter(noisy, sigma=2.0)
        assert smoothed.shape == noisy.shape
        assert smoothed.dtype == np.float32
        # Smoothed should have lower variance than noisy
        assert smoothed.std() < noisy.std()

    def test_median_filter(self, cucim):
        """Median filter reduces impulse noise."""
        rng = np.random.RandomState(42)
        clean = np.ones((32, 32), dtype=np.float32) * 50.0
        noisy = clean.copy()
        # Add salt-and-pepper noise
        for _ in range(30):
            noisy[rng.randint(0, 32), rng.randint(0, 32)] = 255.0
        filtered = cucim.median_filter(noisy, size=3)
        assert filtered.shape == noisy.shape
        # Filtered should be closer to the clean image
        assert np.abs(filtered - clean).mean() < np.abs(noisy - clean).mean()

    def test_volume_computation(self, cucim):
        """mm^3 volume correct given known voxel spacing."""
        # 10x10x10 cube of True voxels = 1000 voxels
        mask = np.zeros((20, 20, 20), dtype=bool)
        mask[5:15, 5:15, 5:15] = True
        spacing = (1.0, 1.0, 1.0)
        vol = cucim.compute_volume_mm3(mask, spacing)
        assert vol == pytest.approx(1000.0)

        # Same mask, 2mm spacing -> 1000 * 8 = 8000 mm^3
        spacing2 = (2.0, 2.0, 2.0)
        vol2 = cucim.compute_volume_mm3(mask, spacing2)
        assert vol2 == pytest.approx(8000.0)

    def test_surface_area(self, cucim):
        """Surface area reasonable for sphere-like object."""
        # Create a sphere of radius 8 in a 32^3 volume
        mask = np.zeros((32, 32, 32), dtype=bool)
        z, y, x = np.ogrid[:32, :32, :32]
        sphere = ((z - 16) ** 2 + (y - 16) ** 2 + (x - 16) ** 2) < 8 ** 2
        mask[sphere] = True

        spacing = (1.0, 1.0, 1.0)
        area = cucim.compute_surface_area_mm2(mask, spacing)

        # Analytic sphere surface area: 4 * pi * r^2 = 4 * pi * 64 ~ 804
        # Discrete approximation will differ but should be in the right ballpark
        expected = 4.0 * np.pi * 8.0 ** 2
        assert area > expected * 0.5, f"Surface area {area} too small (expected ~{expected})"
        assert area < expected * 2.0, f"Surface area {area} too large (expected ~{expected})"
