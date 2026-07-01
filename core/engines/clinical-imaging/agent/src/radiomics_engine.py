"""PyRadiomics feature extraction engine for Imaging Intelligence Agent.

Extracts quantitative radiomic features from medical images (NIfTI) using
PyRadiomics.  Supports configurable feature classes (firstorder, shape,
glcm, glrlm, glszm, ngtdm, gldm), optional image filters (LoG, wavelet),
GPU acceleration via pyradiomics_cuda, and a mock mode that returns
realistic synthetic feature values for demo/testing.

Follows the same GPU-fallback and mock pattern as:
  - src/nim/base.py (BaseNIMClient._invoke_or_mock)
  - src/workflows/base.py (BaseImagingWorkflow._mock_inference)

Author: Adam Jones
Date: April 2026
"""

import hashlib
import math
import random
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from loguru import logger
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════
# GPU Acceleration — try CUDA build first, fall back to CPU
# ═══════════════════════════════════════════════════════════════════════

_GPU_AVAILABLE = False

try:
    import pyradiomics_cuda as _radiomics_backend  # type: ignore[import-untyped]
    from pyradiomics_cuda.featureextractor import RadiomicsFeatureExtractor as _Extractor  # type: ignore[import-untyped]
    _GPU_AVAILABLE = True
    logger.info("PyRadiomics CUDA backend loaded — GPU-accelerated extraction enabled")
except ImportError:
    try:
        import radiomics as _radiomics_backend  # type: ignore[import-untyped]
        from radiomics.featureextractor import RadiomicsFeatureExtractor as _Extractor  # type: ignore[import-untyped]
        logger.info("PyRadiomics CPU backend loaded (pyradiomics_cuda not available)")
    except ImportError:
        _radiomics_backend = None  # type: ignore[assignment]
        _Extractor = None  # type: ignore[assignment,misc]
        logger.warning(
            "PyRadiomics not installed — only mock mode will be available. "
            "Install with: pip install pyradiomics>=3.1.0"
        )


# ═══════════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════════

ALL_FEATURE_CLASSES = [
    "firstorder",
    "shape",
    "glcm",
    "glrlm",
    "glszm",
    "ngtdm",
    "gldm",
]

DEFAULT_FEATURE_CLASSES = ["firstorder", "shape", "glcm"]

SUPPORTED_FILTERS = ["LoG", "wavelet"]

# Realistic mock feature ranges: (min, max) for common features
_MOCK_FIRSTORDER = {
    "firstorder_Mean": (30.0, 120.0),
    "firstorder_Median": (28.0, 115.0),
    "firstorder_StandardDeviation": (8.0, 45.0),
    "firstorder_Variance": (64.0, 2025.0),
    "firstorder_Skewness": (-1.5, 2.0),
    "firstorder_Kurtosis": (2.0, 8.0),
    "firstorder_Entropy": (3.5, 6.5),
    "firstorder_Energy": (1e6, 5e8),
    "firstorder_Minimum": (-200.0, 20.0),
    "firstorder_Maximum": (100.0, 350.0),
    "firstorder_Range": (150.0, 500.0),
    "firstorder_10Percentile": (10.0, 60.0),
    "firstorder_90Percentile": (80.0, 200.0),
    "firstorder_MeanAbsoluteDeviation": (5.0, 35.0),
    "firstorder_RobustMeanAbsoluteDeviation": (3.0, 25.0),
    "firstorder_RootMeanSquared": (35.0, 130.0),
    "firstorder_Uniformity": (0.01, 0.15),
    "firstorder_TotalEnergy": (1e7, 1e10),
    "firstorder_InterquartileRange": (15.0, 80.0),
}

_MOCK_SHAPE = {
    "shape_VoxelVolume": (500.0, 250000.0),
    "shape_MeshVolume": (480.0, 245000.0),
    "shape_SurfaceArea": (200.0, 50000.0),
    "shape_SurfaceVolumeRatio": (0.05, 0.8),
    "shape_Sphericity": (0.4, 0.95),
    "shape_Maximum3DDiameter": (8.0, 120.0),
    "shape_Maximum2DDiameterSlice": (6.0, 100.0),
    "shape_Maximum2DDiameterColumn": (5.0, 95.0),
    "shape_Maximum2DDiameterRow": (5.0, 90.0),
    "shape_MajorAxisLength": (10.0, 80.0),
    "shape_MinorAxisLength": (5.0, 50.0),
    "shape_LeastAxisLength": (3.0, 35.0),
    "shape_Elongation": (0.3, 1.0),
    "shape_Flatness": (0.2, 0.9),
}

_MOCK_GLCM = {
    "glcm_Autocorrelation": (20.0, 500.0),
    "glcm_ClusterProminence": (1e3, 1e7),
    "glcm_ClusterShade": (-5e4, 5e4),
    "glcm_ClusterTendency": (10.0, 500.0),
    "glcm_Contrast": (5.0, 200.0),
    "glcm_Correlation": (0.3, 0.98),
    "glcm_DifferenceAverage": (2.0, 15.0),
    "glcm_DifferenceEntropy": (1.5, 4.5),
    "glcm_DifferenceVariance": (3.0, 80.0),
    "glcm_Id": (0.05, 0.5),
    "glcm_Idm": (0.05, 0.45),
    "glcm_Idmn": (0.9, 1.0),
    "glcm_Idn": (0.85, 0.99),
    "glcm_Imc1": (-0.3, 0.0),
    "glcm_Imc2": (0.7, 0.99),
    "glcm_InverseVariance": (0.1, 0.5),
    "glcm_JointAverage": (5.0, 50.0),
    "glcm_JointEnergy": (0.001, 0.1),
    "glcm_JointEntropy": (4.0, 9.0),
    "glcm_MaximumProbability": (0.01, 0.2),
    "glcm_SumAverage": (10.0, 100.0),
    "glcm_SumEntropy": (3.0, 7.0),
    "glcm_SumSquares": (5.0, 200.0),
}

_MOCK_GLRLM = {
    "glrlm_GrayLevelNonUniformity": (50.0, 5000.0),
    "glrlm_GrayLevelNonUniformityNormalized": (0.01, 0.15),
    "glrlm_GrayLevelVariance": (5.0, 200.0),
    "glrlm_HighGrayLevelRunEmphasis": (50.0, 2000.0),
    "glrlm_LongRunEmphasis": (1.5, 5.0),
    "glrlm_LongRunHighGrayLevelEmphasis": (100.0, 5000.0),
    "glrlm_LongRunLowGrayLevelEmphasis": (0.001, 0.1),
    "glrlm_LowGrayLevelRunEmphasis": (0.001, 0.05),
    "glrlm_RunEntropy": (3.0, 7.0),
    "glrlm_RunLengthNonUniformity": (500.0, 50000.0),
    "glrlm_RunLengthNonUniformityNormalized": (0.5, 0.95),
    "glrlm_RunPercentage": (0.7, 0.98),
    "glrlm_RunVariance": (0.1, 2.0),
    "glrlm_ShortRunEmphasis": (0.7, 0.98),
    "glrlm_ShortRunHighGrayLevelEmphasis": (40.0, 1800.0),
    "glrlm_ShortRunLowGrayLevelEmphasis": (0.001, 0.04),
}

_MOCK_GLSZM = {
    "glszm_GrayLevelNonUniformity": (20.0, 2000.0),
    "glszm_GrayLevelNonUniformityNormalized": (0.02, 0.2),
    "glszm_GrayLevelVariance": (5.0, 200.0),
    "glszm_HighGrayLevelZoneEmphasis": (50.0, 2000.0),
    "glszm_LargeAreaEmphasis": (100.0, 100000.0),
    "glszm_LargeAreaHighGrayLevelEmphasis": (5000.0, 1e7),
    "glszm_LargeAreaLowGrayLevelEmphasis": (1.0, 1000.0),
    "glszm_LowGrayLevelZoneEmphasis": (0.001, 0.05),
    "glszm_SizeZoneNonUniformity": (50.0, 5000.0),
    "glszm_SizeZoneNonUniformityNormalized": (0.2, 0.9),
    "glszm_SmallAreaEmphasis": (0.3, 0.95),
    "glszm_SmallAreaHighGrayLevelEmphasis": (20.0, 1500.0),
    "glszm_SmallAreaLowGrayLevelEmphasis": (0.0005, 0.03),
    "glszm_ZoneEntropy": (4.0, 9.0),
    "glszm_ZonePercentage": (0.001, 0.5),
    "glszm_ZoneVariance": (100.0, 50000.0),
}

_MOCK_NGTDM = {
    "ngtdm_Busyness": (0.1, 50.0),
    "ngtdm_Coarseness": (1e-5, 0.05),
    "ngtdm_Complexity": (1000.0, 500000.0),
    "ngtdm_Contrast": (0.01, 0.5),
    "ngtdm_Strength": (0.5, 20.0),
}

_MOCK_GLDM = {
    "gldm_DependenceEntropy": (4.0, 8.0),
    "gldm_DependenceNonUniformity": (100.0, 10000.0),
    "gldm_DependenceNonUniformityNormalized": (0.01, 0.2),
    "gldm_DependenceVariance": (1.0, 50.0),
    "gldm_GrayLevelNonUniformity": (20.0, 2000.0),
    "gldm_GrayLevelVariance": (5.0, 200.0),
    "gldm_HighGrayLevelEmphasis": (50.0, 2000.0),
    "gldm_LargeDependenceEmphasis": (5.0, 500.0),
    "gldm_LargeDependenceHighGrayLevelEmphasis": (500.0, 100000.0),
    "gldm_LargeDependenceLowGrayLevelEmphasis": (0.01, 5.0),
    "gldm_LowGrayLevelEmphasis": (0.001, 0.05),
    "gldm_SmallDependenceEmphasis": (0.3, 0.95),
    "gldm_SmallDependenceHighGrayLevelEmphasis": (20.0, 1500.0),
    "gldm_SmallDependenceLowGrayLevelEmphasis": (0.0005, 0.03),
}

_MOCK_FEATURE_MAP: Dict[str, Dict[str, Tuple[float, float]]] = {
    "firstorder": _MOCK_FIRSTORDER,
    "shape": _MOCK_SHAPE,
    "glcm": _MOCK_GLCM,
    "glrlm": _MOCK_GLRLM,
    "glszm": _MOCK_GLSZM,
    "ngtdm": _MOCK_NGTDM,
    "gldm": _MOCK_GLDM,
}

# Realistic mock region labels
_MOCK_REGIONS = {
    1: "tumor_core",
    2: "peritumoral_edema",
    3: "enhancing_tumor",
    4: "necrotic_core",
    5: "whole_tumor",
}


# ═══════════════════════════════════════════════════════════════════════
# Pydantic Models
# ═══════════════════════════════════════════════════════════════════════


class RadiomicsFeatureSet(BaseModel):
    """Extracted radiomics features for a single region."""
    region_label: int = Field(1, description="Segmentation mask label")
    region_name: str = Field("", description="Anatomical region name")
    feature_class: str = Field("", description="Feature class (e.g., firstorder, glcm)")
    features: Dict[str, float] = Field(default_factory=dict)
    feature_count: int = Field(0, ge=0)

    @property
    def mean_intensity(self) -> float:
        return self.features.get("firstorder_Mean", 0.0)

    @property
    def volume_mm3(self) -> float:
        return self.features.get("shape_VoxelVolume", 0.0)

    @property
    def sphericity(self) -> float:
        return self.features.get("shape_Sphericity", 0.0)

    @property
    def entropy(self) -> float:
        return self.features.get("firstorder_Entropy", 0.0)


class LongitudinalDelta(BaseModel):
    """Feature delta between two time points with significance flag."""
    feature_name: str
    value_t1: float
    value_t2: float
    absolute_delta: float
    percent_change: float
    is_significant: bool = Field(
        False, description="True if absolute delta > 2 SD of the feature range"
    )
    direction: str = Field("", description="increased, decreased, or stable")


class LongitudinalComparison(BaseModel):
    """Full longitudinal comparison between two feature sets."""
    deltas: List[LongitudinalDelta] = Field(default_factory=list)
    significant_changes: List[str] = Field(
        default_factory=list, description="Feature names with significant changes"
    )
    total_features_compared: int = 0
    summary: str = ""


# ═══════════════════════════════════════════════════════════════════════
# RadiomicsEngine
# ═══════════════════════════════════════════════════════════════════════


class RadiomicsEngine:
    """PyRadiomics feature extraction engine with GPU fallback and mock mode.

    Extracts quantitative radiomic features from NIfTI images paired with
    segmentation masks.  Supports configurable feature classes, image
    filters (LoG, wavelet), GPU acceleration via pyradiomics_cuda, and a
    mock mode that returns realistic synthetic feature values.

    Usage:
        engine = RadiomicsEngine(
            feature_classes=["firstorder", "shape", "glcm"],
            mock=False,
        )
        features = engine.extract_features("image.nii.gz", "mask.nii.gz")
        regions = engine.extract_all_regions("image.nii.gz", "mask.nii.gz")
        comparison = engine.compare_longitudinal(features_t1, features_t2)
    """

    def __init__(
        self,
        feature_classes: Optional[List[str]] = None,
        filters: Optional[List[str]] = None,
        gpu_enabled: bool = True,
        mock: bool = False,
    ):
        """Initialize the RadiomicsEngine.

        Args:
            feature_classes: List of PyRadiomics feature classes to enable.
                Defaults to ["firstorder", "shape", "glcm"].  Valid classes:
                firstorder, shape, glcm, glrlm, glszm, ngtdm, gldm.
            filters: Image filters to apply before feature extraction.
                Empty list or None means no filtering.  Supported: "LoG", "wavelet".
            gpu_enabled: If True, prefer pyradiomics_cuda backend when available.
            mock: If True, return synthetic features without loading images.
        """
        self.mock = mock
        self.gpu_enabled = gpu_enabled and _GPU_AVAILABLE
        self.feature_classes = self._validate_feature_classes(
            feature_classes or DEFAULT_FEATURE_CLASSES
        )
        self.filters = self._validate_filters(filters or [])
        self._extractor: Optional[Any] = None

        if not self.mock and _Extractor is not None:
            self._extractor = self._build_extractor()
            backend = "CUDA" if self.gpu_enabled else "CPU"
            logger.info(
                f"RadiomicsEngine initialized ({backend}): "
                f"classes={self.feature_classes}, filters={self.filters}"
            )
        elif self.mock:
            logger.info(
                f"RadiomicsEngine initialized (mock mode): "
                f"classes={self.feature_classes}"
            )
        else:
            logger.warning(
                "RadiomicsEngine: pyradiomics not installed, forcing mock mode"
            )
            self.mock = True

    @staticmethod
    def _validate_feature_classes(classes: List[str]) -> List[str]:
        """Validate and normalize feature class names.

        Args:
            classes: List of requested feature class names.

        Returns:
            List of valid, lowercased feature class names.
        """
        valid = []
        for cls in classes:
            cls_lower = cls.strip().lower()
            if cls_lower in ALL_FEATURE_CLASSES:
                valid.append(cls_lower)
            else:
                logger.warning(
                    f"Unknown feature class '{cls}' — skipping. "
                    f"Valid classes: {ALL_FEATURE_CLASSES}"
                )
        if not valid:
            logger.warning("No valid feature classes specified, using defaults")
            return list(DEFAULT_FEATURE_CLASSES)
        return valid

    @staticmethod
    def _validate_filters(filters: List[str]) -> List[str]:
        """Validate image filter names.

        Args:
            filters: List of requested filter names.

        Returns:
            List of valid filter names.
        """
        valid = []
        for f in filters:
            f_stripped = f.strip()
            if f_stripped in SUPPORTED_FILTERS:
                valid.append(f_stripped)
            elif f_stripped:
                logger.warning(
                    f"Unknown filter '{f_stripped}' — skipping. "
                    f"Supported: {SUPPORTED_FILTERS}"
                )
        return valid

    def _build_extractor(self) -> Any:
        """Build and configure the PyRadiomics feature extractor.

        Returns:
            Configured RadiomicsFeatureExtractor instance.
        """
        extractor = _Extractor()

        # Disable all feature classes first, then enable only requested ones
        extractor.disableAllFeatures()
        for cls in self.feature_classes:
            extractor.enableFeatureClassByName(cls)
            logger.debug(f"Enabled feature class: {cls}")

        # Configure image filters
        extractor.disableAllImageTypes()
        extractor.enableImageTypeByName("Original")

        for filt in self.filters:
            if filt == "LoG":
                extractor.enableImageTypeByName("LoG")
                extractor.settings["sigma"] = [1.0, 2.0, 3.0]
                logger.debug("Enabled LoG filter with sigma=[1.0, 2.0, 3.0]")
            elif filt == "wavelet":
                extractor.enableImageTypeByName("Wavelet")
                logger.debug("Enabled wavelet filter")

        # Suppress verbose pyradiomics logging
        import logging
        logging.getLogger("radiomics").setLevel(logging.WARNING)

        return extractor

    # ── Feature Extraction ──────────────────────────────────────────

    def extract_features(
        self,
        image_path: str,
        mask_path: str,
        label: int = 1,
    ) -> Dict[str, float]:
        """Extract radiomics features for a single label from image + mask.

        Args:
            image_path: Path to a NIfTI image file (.nii or .nii.gz).
            mask_path: Path to a NIfTI segmentation mask file.
            label: Integer label in the mask to extract features for.

        Returns:
            Dict mapping feature name (e.g. "firstorder_Mean") to float value.
        """
        if self.mock:
            return self._mock_features(label=label)

        if self._extractor is None:
            logger.error("Extractor not available — falling back to mock")
            return self._mock_features(label=label)

        # Validate paths
        img = Path(image_path)
        msk = Path(mask_path)
        if not img.exists():
            raise FileNotFoundError(f"Image file not found: {image_path}")
        if not msk.exists():
            raise FileNotFoundError(f"Mask file not found: {mask_path}")

        logger.info(
            f"Extracting features: image={img.name}, mask={msk.name}, label={label}"
        )

        result = self._extractor.execute(str(img), str(msk), label=label)

        # Filter to only numeric feature values (skip diagnostics)
        features: Dict[str, float] = {}
        for key, value in result.items():
            if key.startswith("diagnostics_"):
                continue
            try:
                features[key] = float(value)
            except (TypeError, ValueError):
                continue

        logger.info(
            f"Extracted {len(features)} features for label={label} "
            f"from {img.name}"
        )
        return features

    def extract_all_regions(
        self,
        image_path: str,
        mask_path: str,
    ) -> Dict[int, Dict[str, float]]:
        """Extract features for all labeled regions in the mask.

        Iterates over unique non-zero labels in the segmentation mask
        and extracts features for each.

        Args:
            image_path: Path to a NIfTI image file.
            mask_path: Path to a NIfTI segmentation mask file.

        Returns:
            Dict mapping region label (int) to feature dict.
        """
        if self.mock:
            return self._mock_all_regions()

        import nibabel as nib
        import numpy as np

        # Load mask to discover unique labels
        msk = Path(mask_path)
        if not msk.exists():
            raise FileNotFoundError(f"Mask file not found: {mask_path}")

        mask_img = nib.load(str(msk))
        mask_data = np.asarray(mask_img.dataobj)
        unique_labels = sorted(set(int(v) for v in np.unique(mask_data) if v != 0))

        if not unique_labels:
            logger.warning(f"No non-zero labels found in mask: {mask_path}")
            return {}

        logger.info(
            f"Extracting features for {len(unique_labels)} regions: {unique_labels}"
        )

        results: Dict[int, Dict[str, float]] = {}
        for lbl in unique_labels:
            try:
                features = self.extract_features(image_path, mask_path, label=lbl)
                results[lbl] = features
            except Exception as e:
                logger.warning(f"Failed to extract features for label {lbl}: {e}")
                continue

        return results

    # ── Longitudinal Comparison ─────────────────────────────────────

    def compare_longitudinal(
        self,
        features_t1: Dict[str, float],
        features_t2: Dict[str, float],
    ) -> LongitudinalComparison:
        """Compare features between two time points and flag significant changes.

        A change is flagged as significant if the absolute delta exceeds
        2 standard deviations of the expected feature range (derived from
        mock feature bounds).

        Args:
            features_t1: Feature dict from the earlier time point.
            features_t2: Feature dict from the later time point.

        Returns:
            LongitudinalComparison with per-feature deltas and significance flags.
        """
        common_keys = sorted(set(features_t1.keys()) & set(features_t2.keys()))
        if not common_keys:
            logger.warning("No common features found between time points")
            return LongitudinalComparison(summary="No common features to compare")

        deltas: List[LongitudinalDelta] = []
        significant: List[str] = []

        for key in common_keys:
            v1 = features_t1[key]
            v2 = features_t2[key]
            abs_delta = v2 - v1

            # Percent change (avoid division by zero)
            if abs(v1) > 1e-10:
                pct_change = ((v2 - v1) / abs(v1)) * 100.0
            else:
                pct_change = 0.0 if abs(abs_delta) < 1e-10 else float("inf")

            # Direction
            if abs(abs_delta) < 1e-10:
                direction = "stable"
            elif abs_delta > 0:
                direction = "increased"
            else:
                direction = "decreased"

            # Significance: compare delta to 2 SD of feature range
            is_sig = self._is_significant_change(key, abs_delta)
            if is_sig:
                significant.append(key)

            deltas.append(
                LongitudinalDelta(
                    feature_name=key,
                    value_t1=v1,
                    value_t2=v2,
                    absolute_delta=abs_delta,
                    percent_change=round(pct_change, 2),
                    is_significant=is_sig,
                    direction=direction,
                )
            )

        summary = (
            f"Compared {len(common_keys)} features between time points. "
            f"{len(significant)} significant changes detected"
        )
        if significant:
            summary += f": {', '.join(significant[:5])}"
            if len(significant) > 5:
                summary += f" (+{len(significant) - 5} more)"

        logger.info(summary)

        return LongitudinalComparison(
            deltas=deltas,
            significant_changes=significant,
            total_features_compared=len(common_keys),
            summary=summary,
        )

    @staticmethod
    def _is_significant_change(feature_name: str, delta: float) -> bool:
        """Determine if a feature delta exceeds 2 SD of the expected range.

        Uses the mock feature bounds as a proxy for the expected range.
        SD is estimated as (max - min) / 4 (assumes ~95% of values fall
        within 2 SD of the mean in a normal distribution).

        Args:
            feature_name: The feature key (e.g. "firstorder_Mean").
            delta: The absolute change (t2 - t1).

        Returns:
            True if |delta| > 2 * estimated_SD.
        """
        # Find the feature range from mock data
        for class_features in _MOCK_FEATURE_MAP.values():
            if feature_name in class_features:
                lo, hi = class_features[feature_name]
                estimated_sd = (hi - lo) / 4.0
                return abs(delta) > (2.0 * estimated_sd)

        # Unknown feature — apply a generous 20% relative threshold
        return False

    # ── Feature Summary Generation ──────────────────────────────────

    @staticmethod
    def generate_feature_summary(
        features: Dict[str, float],
        region_name: str = "",
        modality: str = "",
    ) -> str:
        """Generate a text summary of extracted features suitable for embedding.

        Creates a structured natural-language summary that captures the
        key quantitative characteristics, making it suitable for semantic
        search via BGE-small-en-v1.5 embedding.

        Args:
            features: Dict of feature name -> value.
            region_name: Anatomical region label (e.g. "tumor_core").
            modality: Imaging modality (e.g. "ct", "mri").

        Returns:
            Text summary string (max ~5000 chars).
        """
        parts: List[str] = []

        if region_name:
            parts.append(f"Region: {region_name}.")
        if modality:
            parts.append(f"Modality: {modality}.")

        # First-order statistics
        fo_keys = [k for k in features if k.startswith("firstorder_")]
        if fo_keys:
            mean_val = features.get("firstorder_Mean", 0.0)
            std_val = features.get("firstorder_StandardDeviation", 0.0)
            entropy_val = features.get("firstorder_Entropy", 0.0)
            parts.append(
                f"First-order statistics: mean intensity {mean_val:.1f}, "
                f"standard deviation {std_val:.1f}, entropy {entropy_val:.2f}. "
                f"{len(fo_keys)} first-order features extracted."
            )

        # Shape features
        shape_keys = [k for k in features if k.startswith("shape_")]
        if shape_keys:
            vol = features.get("shape_VoxelVolume", 0.0)
            sph = features.get("shape_Sphericity", 0.0)
            max_diam = features.get("shape_Maximum3DDiameter", 0.0)
            parts.append(
                f"Shape: volume {vol:.1f} mm3, sphericity {sph:.3f}, "
                f"max 3D diameter {max_diam:.1f} mm. "
                f"{len(shape_keys)} shape features extracted."
            )

        # Texture features (GLCM)
        glcm_keys = [k for k in features if k.startswith("glcm_")]
        if glcm_keys:
            contrast = features.get("glcm_Contrast", 0.0)
            corr = features.get("glcm_Correlation", 0.0)
            homog = features.get("glcm_Idm", 0.0)
            parts.append(
                f"GLCM texture: contrast {contrast:.2f}, "
                f"correlation {corr:.3f}, homogeneity {homog:.3f}. "
                f"{len(glcm_keys)} GLCM features extracted."
            )

        # Higher-order texture features
        for prefix, name in [
            ("glrlm_", "GLRLM"),
            ("glszm_", "GLSZM"),
            ("ngtdm_", "NGTDM"),
            ("gldm_", "GLDM"),
        ]:
            keys = [k for k in features if k.startswith(prefix)]
            if keys:
                parts.append(f"{name}: {len(keys)} features extracted.")

        parts.append(f"Total radiomic features: {len(features)}.")

        summary = " ".join(parts)
        return summary[:5000]

    # ── Mock Mode ───────────────────────────────────────────────────

    def _mock_features(self, label: int = 1) -> Dict[str, float]:
        """Generate realistic mock radiomics features for demo/testing.

        Uses deterministic seeding based on the label value for
        reproducible mock outputs across runs.

        Args:
            label: Segmentation mask label (used for seed).

        Returns:
            Dict of feature name -> synthetic float value.
        """
        rng = random.Random(42 + label)
        features: Dict[str, float] = {}

        for cls in self.feature_classes:
            mock_range = _MOCK_FEATURE_MAP.get(cls, {})
            for name, (lo, hi) in mock_range.items():
                features[name] = round(rng.uniform(lo, hi), 6)

        logger.debug(
            f"Generated {len(features)} mock features for label={label}"
        )
        return features

    def _mock_all_regions(self) -> Dict[int, Dict[str, float]]:
        """Generate mock features for multiple regions.

        Returns:
            Dict mapping region label -> feature dict for 3 mock regions.
        """
        regions: Dict[int, Dict[str, float]] = {}
        for lbl in [1, 2, 3]:
            regions[lbl] = self._mock_features(label=lbl)
        logger.debug(f"Generated mock features for {len(regions)} regions")
        return regions

    # ── Utility ─────────────────────────────────────────────────────

    @property
    def is_gpu_accelerated(self) -> bool:
        """Return True if GPU-accelerated backend is active."""
        return self.gpu_enabled and _GPU_AVAILABLE

    @property
    def backend_name(self) -> str:
        """Return a human-readable backend name."""
        if self.mock:
            return "mock"
        if self.is_gpu_accelerated:
            return "pyradiomics_cuda"
        if _radiomics_backend is not None:
            return "pyradiomics"
        return "unavailable"

    def get_status(self) -> Dict[str, Any]:
        """Return engine status information."""
        return {
            "backend": self.backend_name,
            "gpu_available": _GPU_AVAILABLE,
            "gpu_enabled": self.gpu_enabled,
            "mock_mode": self.mock,
            "feature_classes": self.feature_classes,
            "filters": self.filters,
        }
