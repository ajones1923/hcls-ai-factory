"""Tests for src.imaging.preview_generator module.

Uses synthetic NIfTI and DICOM fixtures so no real medical data is required.
MP4 output tests mock ``imageio.mimwrite`` to avoid a hard ffmpeg dependency;
GIF tests use Pillow which is available as a transitive dependency.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import nibabel as nib
import pydicom
import pydicom.uid
import pytest

# Ensure project root is on sys.path so ``from src.…`` imports work.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.imaging.preview_generator import (
    generate_slice_animation,
    generate_thumbnail,
    generate_cine_loop,
    WINDOW_PRESETS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_nifti(tmp_path: Path) -> Path:
    """Create a small synthetic 32x32x20 NIfTI volume."""
    data = np.random.randint(0, 1000, (32, 32, 20), dtype=np.int16)
    img = nib.Nifti1Image(data, affine=np.eye(4))
    path = tmp_path / "test_volume.nii.gz"
    nib.save(img, str(path))
    return path


@pytest.fixture
def sample_dicom_series(tmp_path: Path) -> Path:
    """Create a small synthetic DICOM series (5 slices)."""
    series_dir = tmp_path / "dicom_series"
    series_dir.mkdir()

    study_uid = pydicom.uid.generate_uid()
    series_uid = pydicom.uid.generate_uid()

    for i in range(5):
        file_meta = pydicom.dataset.FileMetaDataset()
        file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian
        file_meta.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()

        ds = pydicom.dataset.FileDataset(
            filename_or_obj="",
            dataset={},
            file_meta=file_meta,
            preamble=b"\x00" * 128,
        )
        ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"
        ds.SOPInstanceUID = file_meta.MediaStorageSOPInstanceUID
        ds.StudyInstanceUID = study_uid
        ds.SeriesInstanceUID = series_uid
        ds.InstanceNumber = i + 1
        ds.Rows = 32
        ds.Columns = 32
        ds.BitsAllocated = 16
        ds.BitsStored = 12
        ds.HighBit = 11
        ds.PixelRepresentation = 0
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.PixelData = np.random.randint(
            0, 4095, (32, 32), dtype=np.uint16
        ).tobytes()

        filepath = str(series_dir / f"slice_{i:03d}.dcm")
        ds.save_as(filepath, write_like_original=False)

    return series_dir


@pytest.fixture
def single_slice_nifti(tmp_path: Path) -> Path:
    """Create a NIfTI volume with only 1 axial slice (32x32x1)."""
    data = np.random.randint(0, 500, (32, 32, 1), dtype=np.int16)
    img = nib.Nifti1Image(data, affine=np.eye(4))
    path = tmp_path / "single_slice.nii.gz"
    nib.save(img, str(path))
    return path


# ---------------------------------------------------------------------------
# Tests — generate_slice_animation
# ---------------------------------------------------------------------------

class TestGenerateSliceAnimation:
    """Tests for the slice animation generator."""

    def test_generate_slice_animation_axial(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Axial sweep produces a GIF file with non-zero size."""
        out = tmp_path / "axial.gif"
        result = generate_slice_animation(sample_nifti, out, axis="axial")
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_generate_slice_animation_sagittal(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Sagittal sweep produces a GIF file with non-zero size."""
        out = tmp_path / "sagittal.gif"
        result = generate_slice_animation(sample_nifti, out, axis="sagittal")
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_generate_slice_animation_coronal(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Coronal sweep produces a GIF file with non-zero size."""
        out = tmp_path / "coronal.gif"
        result = generate_slice_animation(sample_nifti, out, axis="coronal")
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_generate_slice_animation_gif(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Explicit .gif extension uses Pillow-based writer."""
        out = tmp_path / "output.gif"
        result = generate_slice_animation(sample_nifti, out)
        assert result.endswith(".gif")
        assert Path(result).exists()

    def test_generate_slice_animation_mp4(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """MP4 output delegates to imageio.mimwrite (mocked)."""
        out = tmp_path / "output.mp4"
        mock_imageio = MagicMock()
        with patch.dict("sys.modules", {"imageio": mock_imageio}):
            result = generate_slice_animation(sample_nifti, out)
            mock_imageio.mimwrite.assert_called_once()
            assert result == str(out)

    def test_generate_slice_animation_with_brain_window(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Brain window preset applies HU windowing without errors."""
        out = tmp_path / "brain.gif"
        result = generate_slice_animation(sample_nifti, out, window="brain")
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_generate_slice_animation_with_lung_window(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Lung window preset applies HU windowing without errors."""
        out = tmp_path / "lung.gif"
        result = generate_slice_animation(sample_nifti, out, window="lung")
        assert Path(result).exists()

    def test_generate_slice_animation_with_bone_window(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Bone window preset applies HU windowing without errors."""
        out = tmp_path / "bone.gif"
        result = generate_slice_animation(sample_nifti, out, window="bone")
        assert Path(result).exists()

    def test_generate_slice_animation_with_custom_window(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Custom window dict is accepted and applied correctly."""
        out = tmp_path / "custom.gif"
        result = generate_slice_animation(
            sample_nifti, out, window={"center": 100, "width": 200}
        )
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0


# ---------------------------------------------------------------------------
# Tests — generate_thumbnail
# ---------------------------------------------------------------------------

class TestGenerateThumbnail:
    """Tests for the thumbnail generator."""

    def test_generate_thumbnail(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Default thumbnail (middle slice) produces a PNG file."""
        out = tmp_path / "thumb.png"
        result = generate_thumbnail(sample_nifti, out)
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_generate_thumbnail_custom_slice(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """Specifying slice_idx=5 selects that slice without error."""
        out = tmp_path / "thumb_s5.png"
        result = generate_thumbnail(sample_nifti, out, slice_idx=5)
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_generate_thumbnail_auto_middle(
        self, sample_nifti: Path, tmp_path: Path
    ) -> None:
        """When slice_idx is None the middle slice (index 10 for 20 slices) is used."""
        out = tmp_path / "thumb_mid.png"
        with patch(
            "src.imaging.preview_generator._load_volume"
        ) as mock_load:
            # Shape (32, 32, 20) -> middle slice index = 10
            fake_vol = np.zeros((32, 32, 20), dtype=np.float64)
            # Put a distinctive marker at the middle slice
            fake_vol[:, :, 10] = 999.0
            mock_load.return_value = fake_vol

            result = generate_thumbnail(sample_nifti, out, slice_idx=None)
            assert Path(result).exists()


# ---------------------------------------------------------------------------
# Tests — generate_cine_loop
# ---------------------------------------------------------------------------

class TestGenerateCineLoop:
    """Tests for the cine-loop generator."""

    def test_generate_cine_loop(
        self, sample_dicom_series: Path, tmp_path: Path
    ) -> None:
        """Cine loop from a DICOM series directory produces a GIF file."""
        out = tmp_path / "cine.gif"
        result = generate_cine_loop(sample_dicom_series, out, fps=10)
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0


# ---------------------------------------------------------------------------
# Tests — edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge-case and error-handling tests."""

    def test_single_slice_volume(
        self, single_slice_nifti: Path, tmp_path: Path
    ) -> None:
        """A volume with only 1 axial slice should produce a single-frame animation."""
        out = tmp_path / "single.gif"
        result = generate_slice_animation(single_slice_nifti, out, axis="axial")
        assert Path(result).exists()
        assert Path(result).stat().st_size > 0

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        """Passing a non-existent path raises FileNotFoundError."""
        fake = tmp_path / "does_not_exist.nii.gz"
        out = tmp_path / "out.gif"
        with pytest.raises(FileNotFoundError):
            generate_slice_animation(fake, out)

    def test_window_presets_exist(self) -> None:
        """WINDOW_PRESETS contains all five expected presets with correct keys."""
        expected = {"brain", "lung", "bone", "abdomen", "soft_tissue"}
        assert set(WINDOW_PRESETS.keys()) == expected
        for name, preset in WINDOW_PRESETS.items():
            assert "center" in preset, f"Preset '{name}' missing 'center'"
            assert "width" in preset, f"Preset '{name}' missing 'width'"
