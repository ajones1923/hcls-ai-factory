"""Preview generation for medical imaging volumes.

Provides functions to create slice animations (MP4/GIF), cine loops from
DICOM series, single-slice thumbnails, and Orthanc-sourced previews with
configurable windowing presets for CT and MR modalities.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pydicom
import nibabel as nib
from loguru import logger

# ---------------------------------------------------------------------------
# Window presets (HU values for CT)
# ---------------------------------------------------------------------------

WINDOW_PRESETS: Dict[str, Dict[str, int]] = {
    "brain": {"center": 40, "width": 80},
    "lung": {"center": -600, "width": 1500},
    "bone": {"center": 400, "width": 1800},
    "abdomen": {"center": 50, "width": 400},
    "soft_tissue": {"center": 50, "width": 350},
}

# Mapping from human-readable axis name to numpy axis index
_AXIS_MAP: Dict[str, int] = {
    "axial": 2,
    "sagittal": 0,
    "coronal": 1,
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _resolve_window(window: Optional[Union[str, Dict[str, int]]]) -> Optional[Dict[str, int]]:
    """Return a window dict ``{"center": …, "width": …}`` or *None*.

    Parameters
    ----------
    window:
        A preset name (str), an explicit dict, or None for auto-normalisation.

    Raises
    ------
    ValueError
        If the string name is not a known preset.
    """
    if window is None:
        return None
    if isinstance(window, str):
        preset = WINDOW_PRESETS.get(window)
        if preset is None:
            raise ValueError(
                f"Unknown window preset '{window}'. "
                f"Available presets: {list(WINDOW_PRESETS.keys())}"
            )
        return preset
    # Assume dict-like
    return dict(window)


def _apply_window(data: np.ndarray, window: Optional[Dict[str, int]]) -> np.ndarray:
    """Normalize *data* to ``uint8`` using optional windowing.

    Parameters
    ----------
    data:
        Floating-point or integer volume data.
    window:
        ``{"center": …, "width": …}`` or None for min/max auto-scale.

    Returns
    -------
    np.ndarray
        Array scaled to 0-255 as ``uint8``.
    """
    data = data.astype(np.float64)

    if window is not None:
        center = window["center"]
        width = window["width"]
        lower = center - width / 2.0
        upper = center + width / 2.0
        data = np.clip(data, lower, upper)
        if upper != lower:
            data = (data - lower) / (upper - lower) * 255.0
        else:
            data = np.zeros_like(data)
    else:
        dmin = data.min()
        dmax = data.max()
        if dmax != dmin:
            data = (data - dmin) / (dmax - dmin) * 255.0
        else:
            data = np.zeros_like(data)

    return np.clip(data, 0, 255).astype(np.uint8)


def _extract_slices(
    volume: np.ndarray,
    axis: str,
    max_frames: Optional[int] = None,
) -> List[np.ndarray]:
    """Extract 2-D slices along the given anatomical axis.

    Each returned slice is oriented so that the first array dimension
    corresponds to rows (superior-inferior or anterior-posterior) and the
    second to columns.

    Parameters
    ----------
    volume:
        3-D numpy array (X, Y, Z).
    axis:
        One of ``"axial"``, ``"sagittal"``, ``"coronal"``.
    max_frames:
        If set, downsample to at most this many evenly-spaced slices.

    Returns
    -------
    list[np.ndarray]
        List of 2-D slices.
    """
    ax = _AXIS_MAP.get(axis)
    if ax is None:
        raise ValueError(
            f"Unknown axis '{axis}'. Choose from {list(_AXIS_MAP.keys())}."
        )

    n_slices = volume.shape[ax]

    # Determine which slice indices to extract
    if max_frames and n_slices > max_frames:
        indices = np.linspace(0, n_slices - 1, max_frames, dtype=int)
        logger.debug(
            "Downsampling {} slices to {} frames (max_frames={})",
            n_slices,
            max_frames,
            max_frames,
        )
    else:
        indices = range(n_slices)

    slices: List[np.ndarray] = []
    for i in indices:
        slc = np.take(volume, int(i), axis=ax)
        # Rotate so the displayed orientation is conventional
        slc = np.rot90(slc)
        slices.append(slc)
    return slices


def _load_volume(volume_path: Union[str, Path]) -> np.ndarray:
    """Load a 3-D volume from NIfTI file or DICOM series directory.

    Parameters
    ----------
    volume_path:
        Path to a ``.nii`` / ``.nii.gz`` file **or** a directory containing
        DICOM files.

    Returns
    -------
    np.ndarray
        3-D float64 array.

    Raises
    ------
    FileNotFoundError
        If the path does not exist.
    ValueError
        If the path is not a supported format.
    """
    volume_path = Path(volume_path)
    if not volume_path.exists():
        raise FileNotFoundError(f"Volume path does not exist: {volume_path}")

    if volume_path.is_dir():
        return _load_dicom_dir(volume_path)

    suffixes = "".join(volume_path.suffixes).lower()
    if ".nii" in suffixes:
        logger.debug("Loading NIfTI volume: {}", volume_path)
        img = nib.load(str(volume_path))
        data = img.get_fdata()
        # Ensure exactly 3-D
        if data.ndim == 4:
            data = data[..., 0]
        return data

    raise ValueError(
        f"Unsupported volume format: {volume_path.name}. "
        "Expected .nii, .nii.gz, or a directory of DICOM files."
    )


def _load_dicom_dir(dicom_dir: Path) -> np.ndarray:
    """Read a directory of DICOM files and stack into a 3-D array.

    Files are sorted by ``InstanceNumber`` when available.

    Parameters
    ----------
    dicom_dir:
        Directory containing ``.dcm`` (or extensionless) DICOM files.

    Returns
    -------
    np.ndarray
        3-D float64 array (rows, columns, slices).
    """
    dicom_files: List[pydicom.FileDataset] = []
    for p in sorted(dicom_dir.iterdir()):
        if p.is_file():
            try:
                ds = pydicom.dcmread(str(p))
                dicom_files.append(ds)
            except Exception:
                continue

    if not dicom_files:
        raise ValueError(f"No readable DICOM files found in {dicom_dir}")

    # Sort by InstanceNumber if present
    dicom_files.sort(key=lambda d: int(getattr(d, "InstanceNumber", 0)))

    slices = [ds.pixel_array.astype(np.float64) for ds in dicom_files]
    volume = np.stack(slices, axis=-1)
    logger.debug(
        "Loaded DICOM series: {} slices, shape {}",
        len(slices),
        volume.shape,
    )
    return volume


def _write_frames(
    frames: List[np.ndarray],
    output_path: Path,
    fps: int,
) -> None:
    """Write a list of uint8 frames to an MP4 or GIF file.

    Parameters
    ----------
    frames:
        List of 2-D ``uint8`` arrays (H, W).
    output_path:
        Destination file. Extension determines format (``.mp4`` or ``.gif``).
    fps:
        Frames per second.
    """
    if output_path.suffix.lower() == ".gif":
        from PIL import Image  # noqa: F811

        pil_frames = [Image.fromarray(f, mode="L") for f in frames]
        duration_ms = int(1000 / fps)
        pil_frames[0].save(
            str(output_path),
            save_all=True,
            append_images=pil_frames[1:],
            duration=duration_ms,
            loop=0,
        )
        logger.info("Saved GIF ({} frames, {} fps): {}", len(frames), fps, output_path)
    else:
        import imageio

        imageio.mimwrite(str(output_path), frames, fps=fps)
        logger.info("Saved MP4 ({} frames, {} fps): {}", len(frames), fps, output_path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_slice_animation(
    volume_path: Union[str, Path],
    output_path: Union[str, Path],
    axis: str = "axial",
    fps: int = 8,
    colormap: str = "gray",
    window: Optional[Union[str, Dict[str, int]]] = None,
    max_frames: Optional[int] = None,
) -> str:
    """Generate an animated preview (MP4 or GIF) sweeping through volume slices.

    Parameters
    ----------
    volume_path:
        Path to a NIfTI file (``.nii`` / ``.nii.gz``) or a directory of DICOM
        files.
    output_path:
        Destination file path. Use ``.mp4`` or ``.gif`` extension.
    axis:
        Anatomical axis to sweep: ``"axial"`` (z), ``"sagittal"`` (x), or
        ``"coronal"`` (y).
    fps:
        Frames per second for the animation.
    colormap:
        Colormap name (currently only ``"gray"`` is used directly; the
        parameter is reserved for future matplotlib-based colour mapping).
    window:
        Windowing specification — a preset name (``"brain"``, ``"lung"``,
        ``"bone"``, ``"abdomen"``, ``"soft_tissue"``), an explicit dict
        ``{"center": int, "width": int}``, or *None* for auto-normalisation.
    max_frames:
        If set, downsample to at most this many evenly-spaced slices.
        Useful for large volumes where a 200-frame animation is sufficient.

    Returns
    -------
    str
        The ``output_path`` as a string.

    Raises
    ------
    FileNotFoundError
        If *volume_path* does not exist.
    ValueError
        If an unknown axis or window preset is specified.
    """
    volume_path = Path(volume_path)
    output_path = Path(output_path)

    logger.info(
        "Generating slice animation: {} -> {} (axis={}, fps={}, window={})",
        volume_path.name,
        output_path.name,
        axis,
        fps,
        window,
    )

    volume = _load_volume(volume_path)
    win = _resolve_window(window)
    windowed = _apply_window(volume, win)
    slices = _extract_slices(windowed, axis, max_frames=max_frames)

    # Ensure at least one frame
    if not slices:
        raise ValueError(f"No slices extracted along axis '{axis}'.")

    _write_frames(slices, output_path, fps)
    return str(output_path)


def generate_cine_loop(
    dicom_series_dir: Union[str, Path],
    output_path: Union[str, Path],
    fps: int = 15,
    max_frames: Optional[int] = None,
) -> str:
    """Create a cine-loop animation from a temporal DICOM series.

    The function reads all DICOM files in *dicom_series_dir*, sorts them by
    ``InstanceNumber`` (falling back to ``TemporalPositionIndex``), normalises
    pixel data to ``uint8``, and writes an MP4 or GIF.

    Parameters
    ----------
    dicom_series_dir:
        Directory containing the DICOM files of one series.
    output_path:
        Destination file path (``.mp4`` or ``.gif``).
    fps:
        Frames per second.

    Returns
    -------
    str
        The *output_path* as a string.

    Raises
    ------
    FileNotFoundError
        If *dicom_series_dir* does not exist.
    ValueError
        If no DICOM files are found.
    """
    dicom_series_dir = Path(dicom_series_dir)
    output_path = Path(output_path)

    if not dicom_series_dir.exists():
        raise FileNotFoundError(
            f"DICOM series directory does not exist: {dicom_series_dir}"
        )

    logger.info(
        "Generating cine loop: {} -> {} (fps={})",
        dicom_series_dir,
        output_path.name,
        fps,
    )

    datasets: List[pydicom.FileDataset] = []
    for p in sorted(dicom_series_dir.iterdir()):
        if p.is_file():
            try:
                ds = pydicom.dcmread(str(p))
                datasets.append(ds)
            except Exception:
                continue

    if not datasets:
        raise ValueError(
            f"No readable DICOM files found in {dicom_series_dir}"
        )

    # Sort by TemporalPositionIndex first, then InstanceNumber
    def _sort_key(ds: pydicom.FileDataset) -> tuple:
        temporal = int(getattr(ds, "TemporalPositionIndex", 0))
        instance = int(getattr(ds, "InstanceNumber", 0))
        return (temporal, instance)

    datasets.sort(key=_sort_key)

    # Downsample if max_frames is set
    if max_frames and len(datasets) > max_frames:
        indices = np.linspace(0, len(datasets) - 1, max_frames, dtype=int)
        datasets = [datasets[i] for i in indices]
        logger.debug(
            "Downsampled cine loop to {} frames (max_frames={})",
            max_frames,
            max_frames,
        )

    # Extract and normalize frames
    raw_frames = [ds.pixel_array.astype(np.float64) for ds in datasets]
    global_min = min(f.min() for f in raw_frames)
    global_max = max(f.max() for f in raw_frames)

    frames: List[np.ndarray] = []
    for f in raw_frames:
        if global_max != global_min:
            normed = (f - global_min) / (global_max - global_min) * 255.0
        else:
            normed = np.zeros_like(f)
        frames.append(np.clip(normed, 0, 255).astype(np.uint8))

    _write_frames(frames, output_path, fps)
    return str(output_path)


def generate_thumbnail(
    volume_path: Union[str, Path],
    output_path: Union[str, Path],
    slice_idx: Optional[int] = None,
) -> str:
    """Save a single axial slice as a PNG thumbnail.

    Parameters
    ----------
    volume_path:
        Path to a NIfTI file or DICOM series directory.
    output_path:
        Destination PNG path.
    slice_idx:
        Zero-based slice index along the axial (z) axis. If *None*, the
        middle slice is used.

    Returns
    -------
    str
        The *output_path* as a string.

    Raises
    ------
    FileNotFoundError
        If *volume_path* does not exist.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    volume_path = Path(volume_path)
    output_path = Path(output_path)

    volume = _load_volume(volume_path)

    n_slices = volume.shape[2]
    if slice_idx is None:
        slice_idx = n_slices // 2
    slice_idx = max(0, min(slice_idx, n_slices - 1))

    logger.info(
        "Generating thumbnail: {} slice {} -> {}",
        volume_path.name,
        slice_idx,
        output_path.name,
    )

    slc = volume[:, :, slice_idx]
    slc = np.rot90(slc)

    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    ax.imshow(slc, cmap="gray", aspect="equal")
    ax.axis("off")
    fig.savefig(str(output_path), bbox_inches="tight", pad_inches=0, dpi=100)
    plt.close(fig)

    return str(output_path)


def generate_preview_from_orthanc(
    orthanc_url: str,
    orthanc_id: str,
    output_dir: Union[str, Path],
    format: str = "mp4",
    username: str = "admin",
    password: str = "orthanc",
) -> str:
    """Download a series from Orthanc and generate a preview animation.

    The function queries the Orthanc REST API, downloads each instance as raw
    DICOM, saves them to a temporary directory, then delegates to either
    :func:`generate_cine_loop` (for temporal / multi-frame data) or
    :func:`generate_slice_animation` (for spatial stacks).

    Parameters
    ----------
    orthanc_url:
        Base URL of the Orthanc server (e.g. ``http://localhost:8042``).
    orthanc_id:
        Orthanc series identifier.
    output_dir:
        Directory where the output file will be written.
    format:
        Output format extension — ``"mp4"`` or ``"gif"``.
    username:
        Orthanc username for basic auth.
    password:
        Orthanc password for basic auth.

    Returns
    -------
    str
        Path to the generated preview file.

    Raises
    ------
    RuntimeError
        If the Orthanc API request fails.
    """
    import requests

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    orthanc_url = orthanc_url.rstrip("/")
    auth = (username, password)

    logger.info(
        "Fetching series {} from Orthanc at {}",
        orthanc_id,
        orthanc_url,
    )

    # Get list of instances
    resp = requests.get(
        f"{orthanc_url}/series/{orthanc_id}/instances",
        auth=auth,
        timeout=30,
    )
    if resp.status_code != 200:
        raise RuntimeError(
            f"Orthanc API returned status {resp.status_code}: {resp.text}"
        )

    instances = resp.json()
    logger.debug("Found {} instances in series {}", len(instances), orthanc_id)

    # Download each instance to a temp directory (cleaned up after generation)
    tmp_dir_ctx = tempfile.TemporaryDirectory(prefix="orthanc_preview_")
    tmp_dir = Path(tmp_dir_ctx.name)
    temporal_indicators = 0
    _TEMPORAL_CHECK_COUNT = min(5, len(instances))  # Check first N instances

    for idx, inst in enumerate(instances):
        instance_id = inst["ID"]
        dicom_resp = requests.get(
            f"{orthanc_url}/instances/{instance_id}/file",
            auth=auth,
            timeout=60,
        )
        if dicom_resp.status_code != 200:
            logger.warning(
                "Failed to download instance {}: HTTP {}",
                instance_id,
                dicom_resp.status_code,
            )
            continue

        dicom_path = tmp_dir / f"instance_{idx:05d}.dcm"
        dicom_path.write_bytes(dicom_resp.content)

        # Check for temporal data across first few instances
        if idx < _TEMPORAL_CHECK_COUNT:
            try:
                ds = pydicom.dcmread(str(dicom_path))
                if int(getattr(ds, "NumberOfFrames", 1)) > 1:
                    temporal_indicators += 2  # Strong signal
                if int(getattr(ds, "TemporalPositionIndex", 0)) > 0:
                    temporal_indicators += 2
                if int(getattr(ds, "CardiacNumberOfImages", 0)) > 1:
                    temporal_indicators += 2
                if getattr(ds, "TriggerTime", None) is not None:
                    temporal_indicators += 1
                if int(getattr(ds, "NumberOfTemporalPositions", 0)) > 1:
                    temporal_indicators += 2
            except Exception:
                pass

    is_temporal = temporal_indicators >= 2

    output_file = output_dir / f"preview_{orthanc_id}.{format}"

    try:
        if is_temporal:
            logger.info("Detected temporal data — generating cine loop")
            return generate_cine_loop(tmp_dir, output_file, fps=15)
        else:
            logger.info("Detected spatial stack — generating slice animation")
            return generate_slice_animation(tmp_dir, output_file, axis="axial", fps=8)
    finally:
        tmp_dir_ctx.cleanup()
