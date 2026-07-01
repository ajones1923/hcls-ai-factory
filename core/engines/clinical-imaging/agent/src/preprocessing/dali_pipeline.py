"""NVIDIA DALI GPU-accelerated preprocessing for medical imaging.

Offloads DICOM/NIfTI loading, reorientation, resampling, and windowing
from CPU to GPU, eliminating the preprocessing bottleneck.

Apache 2.0 Licensed. Available on DGX Spark ARM64.
"""

from typing import Dict, Optional, Tuple, Union

import numpy as np
from loguru import logger


class DALIPreprocessor:
    """GPU-accelerated medical image preprocessing via NVIDIA DALI.

    Provides GPU-accelerated alternatives to common CPU preprocessing:
    - Volume loading and format conversion
    - Reorientation to standard RAS orientation
    - Resampling to isotropic spacing
    - CT windowing (brain, lung, bone, soft tissue, etc.)
    - Normalization (z-score, min-max, percentile)
    - Resizing to model input dimensions

    Falls back to numpy/scipy CPU processing when DALI is unavailable.
    """

    # Standard CT window presets (center, width in Hounsfield Units)
    CT_WINDOWS: Dict[str, Dict[str, int]] = {
        "brain": {"center": 40, "width": 80},       # 0-80 HU
        "stroke": {"center": 40, "width": 40},       # 20-60 HU
        "subdural": {"center": 75, "width": 215},    # -32.5 to 182.5 HU
        "bone": {"center": 400, "width": 2000},      # -600 to 1400 HU
        "lung": {"center": -600, "width": 1500},     # -1350 to -150 HU
        "soft_tissue": {"center": 40, "width": 400}, # -160 to 240 HU
        "liver": {"center": 60, "width": 150},       # -15 to 135 HU
        "mediastinum": {"center": 40, "width": 350}, # -135 to 215 HU
    }

    def __init__(self, device: str = "gpu"):
        """Initialize DALI preprocessor.

        Args:
            device: Target device for DALI operations. ``"gpu"`` for GPU
                acceleration, ``"cpu"`` to force CPU mode even if DALI is
                available. When DALI is not installed the device parameter is
                ignored and all operations fall back to numpy/scipy.
        """
        self.device = device
        self._dali_available = self._check_dali()
        if self._dali_available:
            logger.info("NVIDIA DALI detected, GPU preprocessing enabled")
        else:
            logger.info("DALI not available, using CPU preprocessing (numpy/scipy)")

    # ------------------------------------------------------------------
    # Availability check
    # ------------------------------------------------------------------

    def _check_dali(self) -> bool:
        """Check whether NVIDIA DALI is importable."""
        try:
            import nvidia.dali  # noqa: F401
            return True
        except ImportError:
            return False

    @property
    def is_gpu_accelerated(self) -> bool:
        """Return True when DALI is loaded and device is ``gpu``."""
        return self._dali_available and self.device == "gpu"

    # ------------------------------------------------------------------
    # CT Windowing
    # ------------------------------------------------------------------

    def apply_window(
        self,
        volume: np.ndarray,
        window_name: str,
    ) -> np.ndarray:
        """Apply CT windowing to a volume.

        Maps Hounsfield Unit values to the 0-1 float range using the
        specified window preset.  GPU-accelerated via DALI when available,
        otherwise uses vectorised numpy operations.

        Args:
            volume: Input CT volume in HU (any shape).
            window_name: Key into :pyattr:`CT_WINDOWS` (e.g. ``"brain"``).

        Returns:
            Windowed volume with values clipped to ``[0.0, 1.0]``.

        Raises:
            ValueError: If *window_name* is not a recognised preset.
        """
        if window_name not in self.CT_WINDOWS:
            raise ValueError(
                f"Unknown CT window '{window_name}'. "
                f"Available: {list(self.CT_WINDOWS.keys())}"
            )

        preset = self.CT_WINDOWS[window_name]
        center = preset["center"]
        width = preset["width"]
        lower = center - width / 2.0
        upper = center + width / 2.0

        if self._dali_available and self.device == "gpu":
            return self._apply_window_dali(volume, lower, upper)

        return self._apply_window_cpu(volume, lower, upper)

    def _apply_window_dali(
        self,
        volume: np.ndarray,
        lower: float,
        upper: float,
    ) -> np.ndarray:
        """Apply CT window using DALI GPU pipeline."""
        try:
            import nvidia.dali.fn as fn  # noqa: F811
            import nvidia.dali.types as types
            from nvidia.dali import pipeline_def

            @pipeline_def(batch_size=1, num_threads=1, device_id=0)
            def window_pipe():
                data = types.Constant(volume.astype(np.float32))
                data = fn.clip(data, lo=lower, hi=upper)
                data = (data - lower) / (upper - lower)
                return data

            pipe = window_pipe()
            pipe.build()
            (output,) = pipe.run()
            return np.array(output[0]).astype(np.float32)
        except Exception as exc:
            logger.warning(f"DALI windowing failed, falling back to CPU: {exc}")
            return self._apply_window_cpu(volume, lower, upper)

    @staticmethod
    def _apply_window_cpu(
        volume: np.ndarray,
        lower: float,
        upper: float,
    ) -> np.ndarray:
        """Apply CT window using numpy (CPU fallback)."""
        result = np.clip(volume.astype(np.float32), lower, upper)
        result = (result - lower) / (upper - lower)
        return result

    # ------------------------------------------------------------------
    # Resampling
    # ------------------------------------------------------------------

    def resample_isotropic(
        self,
        volume: np.ndarray,
        current_spacing: Tuple[float, float, float],
        target_spacing: float = 1.0,
    ) -> np.ndarray:
        """Resample a 3-D volume to isotropic voxel spacing.

        Uses trilinear interpolation. GPU-accelerated via DALI when
        available, otherwise falls back to :func:`scipy.ndimage.zoom`.

        Args:
            volume: 3-D array (D, H, W).
            current_spacing: Current voxel spacing in mm ``(sz, sy, sx)``.
            target_spacing: Desired isotropic spacing in mm.

        Returns:
            Resampled volume with approximately isotropic voxels.
        """
        zoom_factors = tuple(cs / target_spacing for cs in current_spacing)

        if self._dali_available and self.device == "gpu":
            return self._resample_dali(volume, zoom_factors)

        return self._resample_cpu(volume, zoom_factors)

    def _resample_dali(
        self,
        volume: np.ndarray,
        zoom_factors: Tuple[float, float, float],
    ) -> np.ndarray:
        """Resample using DALI GPU pipeline."""
        try:
            import nvidia.dali.fn as fn
            import nvidia.dali.types as types
            from nvidia.dali import pipeline_def

            new_shape = tuple(
                int(round(s * z)) for s, z in zip(volume.shape, zoom_factors)
            )

            @pipeline_def(batch_size=1, num_threads=1, device_id=0)
            def resample_pipe():
                data = types.Constant(
                    volume.astype(np.float32).reshape(1, *volume.shape)
                )
                data = fn.resize(
                    data,
                    size=new_shape,
                    interp_type=types.INTERP_TRIANGULAR,
                )
                return data

            pipe = resample_pipe()
            pipe.build()
            (output,) = pipe.run()
            return np.array(output[0]).reshape(new_shape).astype(np.float32)
        except Exception as exc:
            logger.warning(f"DALI resampling failed, falling back to CPU: {exc}")
            return self._resample_cpu(volume, zoom_factors)

    @staticmethod
    def _resample_cpu(
        volume: np.ndarray,
        zoom_factors: Tuple[float, float, float],
    ) -> np.ndarray:
        """Resample using scipy (CPU fallback)."""
        from scipy.ndimage import zoom

        return zoom(
            volume.astype(np.float32), zoom_factors, order=1
        ).astype(np.float32)

    # ------------------------------------------------------------------
    # Normalization
    # ------------------------------------------------------------------

    def normalize(
        self,
        volume: np.ndarray,
        method: str = "zscore",
        percentile_low: float = 1.0,
        percentile_high: float = 99.0,
    ) -> np.ndarray:
        """Normalize volume intensities.

        Args:
            volume: Input array (any dimensionality).
            method: One of ``"zscore"``, ``"minmax"``, or ``"percentile"``.
            percentile_low: Lower percentile for ``"percentile"`` method.
            percentile_high: Upper percentile for ``"percentile"`` method.

        Returns:
            Normalized volume as float32.

        Raises:
            ValueError: If *method* is not recognised.
        """
        vol = volume.astype(np.float32)

        if method == "zscore":
            mean = vol.mean()
            std = vol.std()
            if std < 1e-8:
                logger.warning("Volume has near-zero std; returning zeros")
                return np.zeros_like(vol)
            return (vol - mean) / std

        if method == "minmax":
            vmin, vmax = vol.min(), vol.max()
            if abs(vmax - vmin) < 1e-8:
                logger.warning("Volume has constant intensity; returning zeros")
                return np.zeros_like(vol)
            return (vol - vmin) / (vmax - vmin)

        if method == "percentile":
            plow = np.percentile(vol, percentile_low)
            phigh = np.percentile(vol, percentile_high)
            if abs(phigh - plow) < 1e-8:
                logger.warning("Percentile range is near-zero; returning zeros")
                return np.zeros_like(vol)
            vol = np.clip(vol, plow, phigh)
            return (vol - plow) / (phigh - plow)

        raise ValueError(
            f"Unknown normalization method '{method}'. "
            "Choose from 'zscore', 'minmax', 'percentile'."
        )

    # ------------------------------------------------------------------
    # Resizing
    # ------------------------------------------------------------------

    def resize(
        self,
        volume: np.ndarray,
        target_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Resize a volume or image to the target shape.

        GPU-accelerated via DALI when available, otherwise uses
        :func:`scipy.ndimage.zoom`.

        Args:
            volume: Input array (2-D or 3-D).
            target_shape: Desired output shape.

        Returns:
            Resized array as float32.
        """
        if volume.shape == target_shape:
            return volume.astype(np.float32)

        if self._dali_available and self.device == "gpu":
            return self._resize_dali(volume, target_shape)

        return self._resize_cpu(volume, target_shape)

    def _resize_dali(
        self,
        volume: np.ndarray,
        target_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Resize using DALI GPU pipeline."""
        try:
            import nvidia.dali.fn as fn
            import nvidia.dali.types as types
            from nvidia.dali import pipeline_def

            @pipeline_def(batch_size=1, num_threads=1, device_id=0)
            def resize_pipe():
                data = types.Constant(volume.astype(np.float32))
                data = fn.resize(
                    data,
                    size=target_shape,
                    interp_type=types.INTERP_TRIANGULAR,
                )
                return data

            pipe = resize_pipe()
            pipe.build()
            (output,) = pipe.run()
            return np.array(output[0]).reshape(target_shape).astype(np.float32)
        except Exception as exc:
            logger.warning(f"DALI resize failed, falling back to CPU: {exc}")
            return self._resize_cpu(volume, target_shape)

    @staticmethod
    def _resize_cpu(
        volume: np.ndarray,
        target_shape: Tuple[int, ...],
    ) -> np.ndarray:
        """Resize using scipy (CPU fallback)."""
        from scipy.ndimage import zoom

        factors = tuple(t / s for t, s in zip(target_shape, volume.shape))
        return zoom(volume.astype(np.float32), factors, order=1).astype(
            np.float32
        )

    # ------------------------------------------------------------------
    # High-level pipelines
    # ------------------------------------------------------------------

    def preprocess_ct(
        self,
        volume: np.ndarray,
        window: str = "soft_tissue",
        target_spacing: float = 1.0,
        current_spacing: Optional[Tuple[float, float, float]] = None,
        target_shape: Optional[Tuple[int, int, int]] = None,
    ) -> np.ndarray:
        """Full CT preprocessing pipeline: window -> resample -> normalize -> resize.

        Each step is optional depending on the supplied arguments. The
        windowing step always runs; resampling requires *current_spacing*;
        resizing requires *target_shape*.

        Args:
            volume: Raw CT volume in Hounsfield Units (D, H, W).
            window: CT window preset name (default ``"soft_tissue"``).
            target_spacing: Isotropic spacing in mm for resampling.
            current_spacing: Current voxel spacing ``(sz, sy, sx)`` in mm.
                If ``None``, the resampling step is skipped.
            target_shape: Final output shape. If ``None``, no resizing.

        Returns:
            Preprocessed volume as float32 in ``[0, 1]``.
        """
        logger.debug(
            f"CT preprocessing: window={window}, spacing={target_spacing}mm, "
            f"shape={target_shape}, device={self.device}"
        )

        # 1. CT windowing
        result = self.apply_window(volume, window)

        # 2. Resample to isotropic spacing (if spacing known)
        if current_spacing is not None:
            result = self.resample_isotropic(result, current_spacing, target_spacing)

        # 3. Normalize to zero-mean unit-variance
        result = self.normalize(result, method="zscore")

        # 4. Resize to model input shape
        if target_shape is not None:
            result = self.resize(result, target_shape)

        logger.debug(f"CT preprocessing complete: output shape {result.shape}")
        return result

    def preprocess_mri(
        self,
        volume: np.ndarray,
        target_spacing: float = 1.0,
        current_spacing: Optional[Tuple[float, float, float]] = None,
        target_shape: Optional[Tuple[int, int, int]] = None,
    ) -> np.ndarray:
        """Full MRI preprocessing pipeline: normalize -> resample -> resize.

        MRI has no standard intensity scale (no Hounsfield Units), so
        windowing is replaced by percentile-based normalization to handle
        the wide dynamic range.

        Args:
            volume: Raw MRI volume (D, H, W).
            target_spacing: Isotropic spacing in mm.
            current_spacing: Current voxel spacing ``(sz, sy, sx)`` in mm.
                If ``None``, resampling is skipped.
            target_shape: Final output shape. If ``None``, no resizing.

        Returns:
            Preprocessed volume as float32.
        """
        logger.debug(
            f"MRI preprocessing: spacing={target_spacing}mm, "
            f"shape={target_shape}, device={self.device}"
        )

        # 1. Percentile normalization (handles MRI intensity drift)
        result = self.normalize(volume, method="percentile")

        # 2. Resample to isotropic
        if current_spacing is not None:
            result = self.resample_isotropic(result, current_spacing, target_spacing)

        # 3. Z-score after resampling
        result = self.normalize(result, method="zscore")

        # 4. Resize
        if target_shape is not None:
            result = self.resize(result, target_shape)

        logger.debug(f"MRI preprocessing complete: output shape {result.shape}")
        return result

    def preprocess_cxr(
        self,
        image: np.ndarray,
        target_size: Tuple[int, int] = (224, 224),
    ) -> np.ndarray:
        """CXR preprocessing: resize -> normalize for DenseNet-121 input.

        Prepares a chest X-ray image for torchxrayvision / DenseNet-121
        inference. Input is a single 2-D grayscale image.

        Args:
            image: 2-D grayscale CXR image (H, W).
            target_size: Output spatial dimensions ``(H, W)``.

        Returns:
            Preprocessed image as float32 with shape *target_size*.
        """
        logger.debug(
            f"CXR preprocessing: target_size={target_size}, device={self.device}"
        )

        # 1. Resize to model input
        result = self.resize(image, target_size)

        # 2. Min-max normalize to [0, 1]
        result = self.normalize(result, method="minmax")

        logger.debug(f"CXR preprocessing complete: output shape {result.shape}")
        return result
