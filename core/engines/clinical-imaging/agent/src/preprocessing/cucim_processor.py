"""NVIDIA cuCIM GPU-accelerated image processing for medical imaging.

Drop-in replacement for scikit-image operations with GPU acceleration.
Used in post-segmentation processing: morphological cleanup, connected
component analysis, region measurement, and filtering.

Apache 2.0 Licensed. Part of RAPIDS.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
from loguru import logger


class CuCIMProcessor:
    """GPU-accelerated image processing via cuCIM.

    Provides GPU-accelerated alternatives to common scikit-image operations
    used in post-segmentation processing:
    - Morphological operations (opening, closing, dilation, erosion)
    - Connected component labeling
    - Region property measurement
    - Image filtering (Gaussian, median)
    - Binary operations (fill holes, remove small objects)

    Falls back to scipy/skimage CPU processing when cuCIM is unavailable.
    """

    def __init__(self):
        """Initialize cuCIM processor with availability check."""
        self._cucim_available = self._check_cucim()
        if self._cucim_available:
            logger.info("NVIDIA cuCIM detected, GPU image processing enabled")
        else:
            logger.info(
                "cuCIM not available, using CPU image processing (scipy/skimage)"
            )

    # ------------------------------------------------------------------
    # Availability
    # ------------------------------------------------------------------

    def _check_cucim(self) -> bool:
        """Check whether cuCIM scikit-image API is importable."""
        try:
            import cucim.skimage  # noqa: F401
            return True
        except ImportError:
            return False

    @property
    def is_gpu_accelerated(self) -> bool:
        """Return True when cuCIM is loaded."""
        return self._cucim_available

    # ------------------------------------------------------------------
    # Morphological operations
    # ------------------------------------------------------------------

    def remove_small_objects(
        self,
        mask: np.ndarray,
        min_size: int = 50,
    ) -> np.ndarray:
        """Remove small connected components from a binary mask.

        Args:
            mask: Binary mask (bool or uint8, any dimensionality).
            min_size: Minimum number of voxels to keep a component.

        Returns:
            Cleaned binary mask with small objects removed.
        """
        binary = mask.astype(bool)

        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.morphology import remove_small_objects as cucim_rso

                gpu_mask = cp.asarray(binary)
                result = cucim_rso(gpu_mask, min_size=min_size)
                return cp.asnumpy(result)
            except Exception as exc:
                logger.warning(
                    f"cuCIM remove_small_objects failed, CPU fallback: {exc}"
                )

        from scipy.ndimage import label

        labeled, num_features = label(binary)
        result = np.zeros_like(binary)
        for i in range(1, num_features + 1):
            component = labeled == i
            if component.sum() >= min_size:
                result[component] = True
        return result

    def binary_closing(
        self,
        mask: np.ndarray,
        radius: int = 3,
    ) -> np.ndarray:
        """Morphological closing to fill small gaps in a binary mask.

        Closing = dilation followed by erosion. Useful for filling narrow
        breaks in segmentation boundaries.

        Args:
            mask: Binary mask (bool or uint8).
            radius: Radius of the spherical/disk structuring element.

        Returns:
            Closed binary mask.
        """
        binary = mask.astype(bool)

        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.morphology import (
                    binary_closing as cucim_closing,
                    ball,
                    disk,
                )

                gpu_mask = cp.asarray(binary)
                if binary.ndim == 3:
                    selem = ball(radius)
                else:
                    selem = disk(radius)
                result = cucim_closing(gpu_mask, selem)
                return cp.asnumpy(result)
            except Exception as exc:
                logger.warning(
                    f"cuCIM binary_closing failed, CPU fallback: {exc}"
                )

        from scipy.ndimage import binary_closing as scipy_closing
        from scipy.ndimage import generate_binary_structure

        ndim = binary.ndim
        struct = generate_binary_structure(ndim, 1)
        return scipy_closing(binary, structure=struct, iterations=radius).astype(
            bool
        )

    def binary_opening(
        self,
        mask: np.ndarray,
        radius: int = 2,
    ) -> np.ndarray:
        """Morphological opening to remove small protrusions.

        Opening = erosion followed by dilation.

        Args:
            mask: Binary mask (bool or uint8).
            radius: Radius of the structuring element.

        Returns:
            Opened binary mask.
        """
        binary = mask.astype(bool)

        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.morphology import (
                    binary_opening as cucim_opening,
                    ball,
                    disk,
                )

                gpu_mask = cp.asarray(binary)
                if binary.ndim == 3:
                    selem = ball(radius)
                else:
                    selem = disk(radius)
                result = cucim_opening(gpu_mask, selem)
                return cp.asnumpy(result)
            except Exception as exc:
                logger.warning(
                    f"cuCIM binary_opening failed, CPU fallback: {exc}"
                )

        from scipy.ndimage import binary_opening as scipy_opening
        from scipy.ndimage import generate_binary_structure

        ndim = binary.ndim
        struct = generate_binary_structure(ndim, 1)
        return scipy_opening(binary, structure=struct, iterations=radius).astype(
            bool
        )

    def fill_holes(self, mask: np.ndarray) -> np.ndarray:
        """Fill holes in a binary mask.

        A hole is a region of background completely surrounded by
        foreground in the binary mask.

        Args:
            mask: Binary mask (bool or uint8).

        Returns:
            Binary mask with interior holes filled.
        """
        binary = mask.astype(bool)

        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.morphology import (
                    remove_small_holes as cucim_rsh,
                )

                gpu_mask = cp.asarray(binary)
                result = cucim_rsh(gpu_mask, area_threshold=binary.size)
                return cp.asnumpy(result)
            except Exception as exc:
                logger.warning(f"cuCIM fill_holes failed, CPU fallback: {exc}")

        from scipy.ndimage import binary_fill_holes

        return binary_fill_holes(binary).astype(bool)

    # ------------------------------------------------------------------
    # Connected component analysis
    # ------------------------------------------------------------------

    def label_connected_components(
        self,
        mask: np.ndarray,
    ) -> Tuple[np.ndarray, int]:
        """Label connected components in a binary mask.

        Args:
            mask: Binary mask (bool or uint8).

        Returns:
            Tuple of ``(labeled_array, num_features)`` where each
            connected component has a unique integer label starting at 1.
        """
        binary = mask.astype(bool)

        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.measure import label as cucim_label

                gpu_mask = cp.asarray(binary)
                labeled = cucim_label(gpu_mask)
                num = int(labeled.max())
                return cp.asnumpy(labeled), num
            except Exception as exc:
                logger.warning(
                    f"cuCIM label failed, CPU fallback: {exc}"
                )

        from scipy.ndimage import label

        return label(binary)

    def region_properties(
        self,
        label_image: np.ndarray,
        intensity_image: Optional[np.ndarray] = None,
    ) -> List[Dict]:
        """Measure region properties for each labeled component.

        Computes area (voxel count), centroid, bounding box, and
        optionally mean intensity for each region.

        Args:
            label_image: Integer-labeled array from :meth:`label_connected_components`.
            intensity_image: Optional grayscale image for intensity measurements.

        Returns:
            List of dicts with keys ``label``, ``area``, ``centroid``,
            ``bbox``, and optionally ``mean_intensity``.
        """
        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.measure import regionprops as cucim_regionprops

                gpu_labels = cp.asarray(label_image)
                gpu_intensity = (
                    cp.asarray(intensity_image) if intensity_image is not None else None
                )
                props = cucim_regionprops(gpu_labels, intensity_image=gpu_intensity)
                return self._props_to_dicts(props, has_intensity=intensity_image is not None)
            except Exception as exc:
                logger.warning(
                    f"cuCIM regionprops failed, CPU fallback: {exc}"
                )

        return self._region_properties_cpu(label_image, intensity_image)

    @staticmethod
    def _props_to_dicts(props, has_intensity: bool = False) -> List[Dict]:
        """Convert regionprops objects to plain dicts."""
        results = []
        for p in props:
            entry: Dict = {
                "label": int(p.label),
                "area": int(p.area),
                "centroid": tuple(float(c) for c in p.centroid),
                "bbox": tuple(int(b) for b in p.bbox),
            }
            if has_intensity:
                try:
                    entry["mean_intensity"] = float(p.mean_intensity)
                except Exception:
                    pass
            results.append(entry)
        return results

    @staticmethod
    def _region_properties_cpu(
        label_image: np.ndarray,
        intensity_image: Optional[np.ndarray] = None,
    ) -> List[Dict]:
        """Compute region properties using scipy (CPU fallback)."""
        from scipy.ndimage import find_objects, label

        results: List[Dict] = []
        num_labels = int(label_image.max())

        slices = find_objects(label_image)
        for i in range(1, num_labels + 1):
            if i - 1 >= len(slices) or slices[i - 1] is None:
                continue
            sl = slices[i - 1]
            component = label_image[sl] == i
            area = int(component.sum())

            # Centroid in global coordinates
            coords = np.argwhere(label_image == i)
            centroid = tuple(float(c) for c in coords.mean(axis=0))

            # Bounding box: (min_row, min_col, ..., max_row, max_col, ...)
            bbox_min = tuple(int(s.start) for s in sl)
            bbox_max = tuple(int(s.stop) for s in sl)
            bbox = bbox_min + bbox_max

            entry: Dict = {
                "label": i,
                "area": area,
                "centroid": centroid,
                "bbox": bbox,
            }
            if intensity_image is not None:
                region_intensities = intensity_image[sl][component]
                entry["mean_intensity"] = float(region_intensities.mean())
            results.append(entry)
        return results

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def gaussian_filter(
        self,
        image: np.ndarray,
        sigma: float = 1.0,
    ) -> np.ndarray:
        """Apply Gaussian smoothing filter.

        Args:
            image: Input image or volume (any dimensionality).
            sigma: Standard deviation of the Gaussian kernel.

        Returns:
            Smoothed image as float32.
        """
        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.filters import gaussian as cucim_gaussian

                gpu_img = cp.asarray(image.astype(np.float32))
                result = cucim_gaussian(gpu_img, sigma=sigma)
                return cp.asnumpy(result).astype(np.float32)
            except Exception as exc:
                logger.warning(
                    f"cuCIM gaussian_filter failed, CPU fallback: {exc}"
                )

        from scipy.ndimage import gaussian_filter as scipy_gaussian

        return scipy_gaussian(image.astype(np.float32), sigma=sigma).astype(
            np.float32
        )

    def median_filter(
        self,
        image: np.ndarray,
        size: int = 3,
    ) -> np.ndarray:
        """Apply median filter for noise reduction.

        Args:
            image: Input image or volume.
            size: Side length of the cubic/square filter kernel.

        Returns:
            Filtered image as float32.
        """
        if self._cucim_available:
            try:
                import cupy as cp
                from cucim.skimage.filters import median as cucim_median

                gpu_img = cp.asarray(image.astype(np.float32))
                footprint = np.ones([size] * image.ndim, dtype=bool)
                result = cucim_median(gpu_img, footprint=cp.asarray(footprint))
                return cp.asnumpy(result).astype(np.float32)
            except Exception as exc:
                logger.warning(
                    f"cuCIM median_filter failed, CPU fallback: {exc}"
                )

        from scipy.ndimage import median_filter as scipy_median

        return scipy_median(image.astype(np.float32), size=size).astype(
            np.float32
        )

    # ------------------------------------------------------------------
    # Volumetric measurements
    # ------------------------------------------------------------------

    def compute_volume_mm3(
        self,
        mask: np.ndarray,
        voxel_spacing: Tuple[float, float, float],
    ) -> float:
        """Compute volume in cubic millimetres from a binary mask.

        Args:
            mask: Binary mask (bool or uint8).
            voxel_spacing: Voxel dimensions in mm ``(sz, sy, sx)``.

        Returns:
            Volume in mm^3.
        """
        voxel_volume = float(np.prod(voxel_spacing))
        num_voxels = int(mask.astype(bool).sum())
        return num_voxels * voxel_volume

    def compute_surface_area_mm2(
        self,
        mask: np.ndarray,
        voxel_spacing: Tuple[float, float, float],
    ) -> float:
        """Estimate surface area in mm^2 using marching cubes.

        Requires a 3-D binary mask. Falls back to a simple voxel-face
        counting approximation if skimage is not available.

        Args:
            mask: 3-D binary mask (D, H, W).
            voxel_spacing: Voxel dimensions in mm ``(sz, sy, sx)``.

        Returns:
            Estimated surface area in mm^2.
        """
        binary = mask.astype(bool).astype(np.float32)

        try:
            from skimage.measure import marching_cubes

            verts, faces, _, _ = marching_cubes(binary, level=0.5, spacing=voxel_spacing)
            # Compute area of each triangle face
            v0 = verts[faces[:, 0]]
            v1 = verts[faces[:, 1]]
            v2 = verts[faces[:, 2]]
            cross = np.cross(v1 - v0, v2 - v0)
            area = 0.5 * np.linalg.norm(cross, axis=1).sum()
            return float(area)
        except ImportError:
            logger.warning(
                "skimage not available for marching_cubes; "
                "using voxel-face approximation"
            )
        except Exception as exc:
            logger.warning(
                f"marching_cubes failed, using voxel-face approximation: {exc}"
            )

        # Fallback: count exposed voxel faces
        return self._surface_area_voxel_faces(binary, voxel_spacing)

    @staticmethod
    def _surface_area_voxel_faces(
        mask: np.ndarray,
        spacing: Tuple[float, float, float],
    ) -> float:
        """Estimate surface area by counting exposed voxel faces."""
        sz, sy, sx = spacing
        face_areas = [sy * sx, sz * sx, sz * sy]  # per axis
        total = 0.0
        padded = np.pad(mask, 1, mode="constant", constant_values=0)
        for axis in range(3):
            diff = np.diff(padded, axis=axis)
            total += float(np.abs(diff).sum()) * face_areas[axis]
        return total
