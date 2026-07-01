"""Preview generation endpoints for medical imaging volumes.

Provides slice animations (GIF/MP4), thumbnails, and Orthanc-sourced
previews for NIfTI volumes and DICOM series.  Generated files are
cached on disk to avoid redundant computation.

Author: Adam Jones
Date: March 2026
"""

from typing import Optional

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import StreamingResponse
from loguru import logger
from pathlib import Path

from config.settings import settings

# =====================================================================
# Router
# =====================================================================

router = APIRouter(prefix="/preview", tags=["Preview"])


# =====================================================================
# Helpers
# =====================================================================


def _media_type(fmt: str) -> str:
    """Return the HTTP media type for a preview format string."""
    return {
        "mp4": "video/mp4",
        "gif": "image/gif",
        "png": "image/png",
    }.get(fmt, "application/octet-stream")


def _stream_file(path: Path, media_type: str) -> StreamingResponse:
    """Return a StreamingResponse that reads *path* in chunks."""

    def _iter():
        with open(path, "rb") as fh:
            while chunk := fh.read(64 * 1024):
                yield chunk

    return StreamingResponse(_iter(), media_type=media_type)


# =====================================================================
# GET /preview/sample/{volume_name}
# =====================================================================


@router.get("/sample/{volume_name}")
async def preview_sample(
    volume_name: str,
    axis: str = Query("axial", pattern="^(axial|sagittal|coronal)$"),
    fps: int = Query(8, ge=1, le=60),
    format: str = Query("gif", pattern="^(mp4|gif)$"),
    window: Optional[str] = Query(None, pattern="^(brain|lung|bone|abdomen|soft_tissue)$"),
):
    """Generate a slice animation preview for a sample NIfTI volume.

    Looks up the volume in ``{DATA_DIR}/sample_images/{volume_name}.nii.gz``,
    generates an animated GIF or MP4, caches the result, and streams it back.
    """
    # Resolve volume path
    volume_path = Path(settings.DATA_DIR) / "sample_images" / f"{volume_name}.nii.gz"
    if not volume_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Sample volume not found: {volume_name}",
        )

    # Ensure cache directory exists
    cache_dir = Path(settings.PREVIEW_CACHE_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Build cache filename
    window_tag = window or "auto"
    cache_name = f"{volume_name}_{axis}_{window_tag}_{fps}fps.{format}"
    cache_path = cache_dir / cache_name

    # Serve from cache if available
    if cache_path.exists():
        logger.debug(f"Serving cached preview: {cache_path}")
        return _stream_file(cache_path, _media_type(format))

    # Generate preview
    from src.imaging.preview_generator import generate_slice_animation

    try:
        output = generate_slice_animation(
            volume_path=str(volume_path),
            output_path=str(cache_path),
            axis=axis,
            fps=fps,
            colormap="gray",
            window=window,
            max_frames=settings.PREVIEW_MAX_FRAMES,
        )
        logger.info(f"Generated preview: {output}")
    except Exception as exc:
        logger.error(f"Preview generation failed for {volume_name}: {exc}")
        raise HTTPException(
            status_code=500,
            detail=f"Preview generation failed: {exc}",
        )

    return _stream_file(cache_path, _media_type(format))


# =====================================================================
# GET /preview/orthanc/{orthanc_id}
# =====================================================================


@router.get("/orthanc/{orthanc_id}")
async def preview_orthanc(
    orthanc_id: str,
    format: str = Query("gif", pattern="^(mp4|gif)$"),
):
    """Generate a preview from an Orthanc DICOM series.

    Fetches the series from the configured Orthanc server, converts it
    to an animated GIF or MP4, caches the result, and streams it back.
    """
    cache_dir = Path(settings.PREVIEW_CACHE_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Serve from cache if available
    cache_path = cache_dir / f"preview_{orthanc_id}.{format}"
    if cache_path.exists():
        logger.debug(f"Serving cached Orthanc preview: {cache_path}")
        return _stream_file(cache_path, _media_type(format))

    from src.imaging.preview_generator import generate_preview_from_orthanc

    try:
        output = generate_preview_from_orthanc(
            orthanc_url=settings.ORTHANC_URL,
            orthanc_id=orthanc_id,
            output_dir=str(cache_dir),
            format=format,
            username=settings.ORTHANC_USERNAME,
            password=settings.ORTHANC_PASSWORD,
        )
        output_path = Path(output)
        logger.info(f"Generated Orthanc preview: {output_path}")
    except Exception as exc:
        logger.error(f"Orthanc preview failed for {orthanc_id}: {exc}")
        raise HTTPException(
            status_code=502,
            detail=f"Failed to generate preview from Orthanc: {exc}",
        )

    return _stream_file(output_path, _media_type(format))


# =====================================================================
# GET /preview/thumbnail/{volume_name}
# =====================================================================


@router.get("/thumbnail/{volume_name}")
async def preview_thumbnail(
    volume_name: str,
):
    """Generate a single-slice PNG thumbnail for a sample NIfTI volume."""
    volume_path = Path(settings.DATA_DIR) / "sample_images" / f"{volume_name}.nii.gz"
    if not volume_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Sample volume not found: {volume_name}",
        )

    cache_dir = Path(settings.PREVIEW_CACHE_DIR)
    cache_dir.mkdir(parents=True, exist_ok=True)

    cache_name = f"{volume_name}_thumb.png"
    cache_path = cache_dir / cache_name

    # Serve from cache if available
    if cache_path.exists():
        logger.debug(f"Serving cached thumbnail: {cache_path}")
        return _stream_file(cache_path, "image/png")

    # Generate thumbnail
    from src.imaging.preview_generator import generate_thumbnail

    try:
        output = generate_thumbnail(
            volume_path=str(volume_path),
            output_path=str(cache_path),
        )
        logger.info(f"Generated thumbnail: {output}")
    except Exception as exc:
        logger.error(f"Thumbnail generation failed for {volume_name}: {exc}")
        raise HTTPException(
            status_code=500,
            detail=f"Thumbnail generation failed: {exc}",
        )

    return _stream_file(cache_path, "image/png")


# =====================================================================
# GET /preview/presets
# =====================================================================


@router.get("/presets")
async def preview_presets():
    """Return the available windowing presets for medical imaging previews.

    Each preset is a ``(window_center, window_width)`` pair commonly used
    in radiology for optimising contrast in different tissue types.
    """
    from src.imaging.preview_generator import WINDOW_PRESETS

    return WINDOW_PRESETS
