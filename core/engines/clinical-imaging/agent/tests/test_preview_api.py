"""Tests for the preview API endpoints.

Verifies slice animation, thumbnail, Orthanc preview, and presets
endpoints using FastAPI TestClient with mocked preview generator
functions and settings.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so imports work regardless of
# how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =====================================================================
# Fixtures
# =====================================================================


@pytest.fixture
def client(tmp_path):
    """TestClient with sample data directory and mocked settings."""
    # Set up sample data directory
    sample_dir = tmp_path / "data" / "sample_images"
    sample_dir.mkdir(parents=True)
    (sample_dir / "sample_ct_chest.nii.gz").write_bytes(b"fake-nifti-data")

    cache_dir = tmp_path / "data" / "cache" / "previews"
    cache_dir.mkdir(parents=True)

    with patch("api.routes.preview.settings") as mock_settings:
        mock_settings.DATA_DIR = tmp_path / "data"
        mock_settings.PREVIEW_CACHE_DIR = str(cache_dir)
        mock_settings.PREVIEW_DEFAULT_FPS = 8
        mock_settings.PREVIEW_DEFAULT_FORMAT = "gif"
        mock_settings.PREVIEW_MAX_FRAMES = 200
        mock_settings.ORTHANC_URL = "http://localhost:8042"
        mock_settings.ORTHANC_USERNAME = "admin"
        mock_settings.ORTHANC_PASSWORD = "orthanc"

        from api.main import app
        from fastapi.testclient import TestClient

        yield TestClient(app, raise_server_exceptions=False)


@pytest.fixture
def cache_dir(tmp_path):
    """Return the cache directory path consistent with the client fixture."""
    return tmp_path / "data" / "cache" / "previews"


# =====================================================================
# Minimal 1-pixel GIF (35 bytes)
# =====================================================================

# GIF89a 1x1 pixel, transparent, minimal valid binary
TINY_GIF = (
    b"GIF89a\x01\x00\x01\x00\x80\x00\x00\xff\xff\xff"
    b"\x00\x00\x00!\xf9\x04\x00\x00\x00\x00\x00,"
    b"\x00\x00\x00\x00\x01\x00\x01\x00\x00\x02\x02D\x01\x00;"
)

# Minimal 1-pixel PNG (67 bytes)
TINY_PNG = (
    b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01"
    b"\x00\x00\x00\x01\x08\x02\x00\x00\x00\x90wS\xde\x00"
    b"\x00\x00\x0cIDATx\x9cc\xf8\x0f\x00\x00\x01\x01\x00"
    b"\x05\x18\xd8N\x00\x00\x00\x00IEND\xaeB`\x82"
)


# =====================================================================
# Tests
# =====================================================================


class TestPreviewSample:
    """Tests for GET /preview/sample/{volume_name}."""

    @patch("src.imaging.preview_generator.generate_slice_animation")
    def test_preview_sample_returns_gif(self, mock_gen, client, tmp_path):
        """Generating a sample preview should return 200 with image/gif."""
        cache_dir = tmp_path / "data" / "cache" / "previews"

        def _fake_generate(volume_path, output_path, **kwargs):
            Path(output_path).write_bytes(TINY_GIF)
            return output_path

        mock_gen.side_effect = _fake_generate

        response = client.get("/preview/sample/sample_ct_chest?axis=axial&fps=8&format=gif")

        assert response.status_code == 200
        assert response.headers["content-type"] == "image/gif"
        assert response.content == TINY_GIF
        mock_gen.assert_called_once()

    def test_preview_sample_not_found(self, client):
        """Requesting a nonexistent volume should return 404."""
        response = client.get("/preview/sample/nonexistent_volume")

        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()

    @patch("src.imaging.preview_generator.generate_slice_animation")
    def test_preview_sample_cached(self, mock_gen, client, tmp_path):
        """A cached preview should be served without calling the generator."""
        cache_dir = tmp_path / "data" / "cache" / "previews"
        cache_file = cache_dir / "sample_ct_chest_axial_auto_8fps.gif"
        cache_file.write_bytes(TINY_GIF)

        response = client.get("/preview/sample/sample_ct_chest?axis=axial&fps=8&format=gif")

        assert response.status_code == 200
        assert response.content == TINY_GIF
        mock_gen.assert_not_called()


class TestPreviewThumbnail:
    """Tests for GET /preview/thumbnail/{volume_name}."""

    @patch("src.imaging.preview_generator.generate_thumbnail")
    def test_preview_thumbnail_returns_png(self, mock_thumb, client, tmp_path):
        """Generating a thumbnail should return 200 with image/png."""

        def _fake_thumbnail(volume_path, output_path, **kwargs):
            Path(output_path).write_bytes(TINY_PNG)
            return output_path

        mock_thumb.side_effect = _fake_thumbnail

        response = client.get("/preview/thumbnail/sample_ct_chest")

        assert response.status_code == 200
        assert response.headers["content-type"] == "image/png"
        assert response.content == TINY_PNG
        mock_thumb.assert_called_once()

    def test_preview_thumbnail_not_found(self, client):
        """Requesting a thumbnail for a nonexistent volume should return 404."""
        response = client.get("/preview/thumbnail/nonexistent_volume")

        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()


class TestPreviewOrthanc:
    """Tests for GET /preview/orthanc/{orthanc_id}."""

    @patch("src.imaging.preview_generator.generate_preview_from_orthanc")
    def test_preview_orthanc_success(self, mock_orthanc, client, tmp_path):
        """A successful Orthanc preview should return 200."""
        cache_dir = tmp_path / "data" / "cache" / "previews"
        output_file = cache_dir / "orthanc_abc123.gif"

        def _fake_orthanc(orthanc_url, orthanc_id, output_dir, **kwargs):
            output_file.write_bytes(TINY_GIF)
            return str(output_file)

        mock_orthanc.side_effect = _fake_orthanc

        response = client.get("/preview/orthanc/abc123?format=gif")

        assert response.status_code == 200
        assert response.headers["content-type"] == "image/gif"
        assert response.content == TINY_GIF
        mock_orthanc.assert_called_once()

    @patch("src.imaging.preview_generator.generate_preview_from_orthanc")
    def test_preview_orthanc_failure(self, mock_orthanc, client):
        """An Orthanc connection failure should return 502."""
        mock_orthanc.side_effect = ConnectionError("Orthanc unreachable")

        response = client.get("/preview/orthanc/bad_id?format=gif")

        assert response.status_code == 502
        assert "orthanc" in response.json()["detail"].lower()


class TestPreviewPresets:
    """Tests for GET /preview/presets."""

    def test_preview_presets_endpoint(self, client):
        """The presets endpoint should return all 5 window presets."""
        response = client.get("/preview/presets")

        assert response.status_code == 200
        data = response.json()
        assert len(data) == 5
        for key in ("brain", "lung", "bone", "abdomen", "soft_tissue"):
            assert key in data
            assert "center" in data[key]
            assert "width" in data[key]
