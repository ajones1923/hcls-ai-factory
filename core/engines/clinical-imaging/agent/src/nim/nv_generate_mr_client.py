"""NV-Generate-MR client for synthetic MRI volume generation.

Generates synthetic MR volumes with multi-contrast support:
T1w, T2w, FLAIR, SWI (brain), and body MRI.

License: Check per release.
"""

import random
import time
from enum import Enum
from typing import Any, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field

from .base import BaseNIMClient


# ═══════════════════════════════════════════════════════════════════════
# CONFIG & RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════

VALID_MR_BODY_REGIONS: List[str] = [
    "brain", "spine", "abdomen", "pelvis", "prostate",
]

VALID_MR_PATHOLOGY_TYPES: List[str] = [
    "tumor", "lesion", "infarct",
]


class MRContrastType(str, Enum):
    """Supported MRI contrast weightings."""

    T1W = "t1w"
    T2W = "t2w"
    FLAIR = "flair"
    SWI = "swi"
    DWI = "dwi"


class SyntheticMRConfig(BaseModel):
    """Configuration for NV-Generate-MR synthetic MRI volume generation."""

    body_region: str = Field("brain", description="Target body region")
    contrast: MRContrastType = Field(MRContrastType.T1W, description="MRI contrast weighting")
    resolution_mm: float = Field(1.0, ge=0.5, le=5.0, description="Voxel spacing in mm")
    volume_size: List[int] = Field(
        default_factory=lambda: [256, 256, 256],
        description="Volume dimensions [W, H, D]",
    )
    skull_stripped: bool = Field(False, description="Remove skull (brain only)")
    include_pathology: bool = Field(False, description="Inject synthetic pathology")
    pathology_type: Optional[str] = Field(None, description="Type: tumor, lesion, infarct")
    seed: Optional[int] = Field(None, description="Random seed for reproducibility")


class SyntheticMRResult(BaseModel):
    """Result from NV-Generate-MR synthetic MRI volume generation."""

    volume_path: Optional[str] = None
    mask_path: Optional[str] = None
    body_region: str = ""
    contrast: str = ""
    resolution_mm: float = 1.0
    volume_shape: List[int] = Field(default_factory=lambda: [256, 256, 256])
    generation_time_ms: float = 0.0
    model_name: str = "nv-generate-mr"
    is_mock: bool = False


# ═══════════════════════════════════════════════════════════════════════
# CLIENT
# ═══════════════════════════════════════════════════════════════════════


class NVGenerateMRClient(BaseNIMClient):
    """Client for NV-Generate-MR synthetic MRI volume generation.

    NV-Generate-MR produces synthetic MR volumes with multi-contrast
    support (T1w, T2w, FLAIR, SWI, DWI) for brain, spine, abdomen,
    pelvis, and prostate imaging. Useful for:
      - Training data augmentation for MRI-based models
      - Privacy-preserving synthetic MRI datasets
      - Multi-contrast paired generation
      - Rare pathology simulation (tumors, lesions, infarcts)
    """

    def __init__(self, base_url: str, mock_enabled: bool = True):
        super().__init__(base_url, service_name="nv-generate-mr", mock_enabled=mock_enabled)

    def generate(
        self,
        config: Optional[SyntheticMRConfig] = None,
    ) -> SyntheticMRResult:
        """Generate a synthetic MR volume.

        Args:
            config: Generation configuration. Uses defaults if None.

        Returns:
            SyntheticMRResult with paths and generation metadata.
        """
        if config is None:
            config = SyntheticMRConfig()

        logger.info(
            f"Generating synthetic MR: region={config.body_region}, "
            f"contrast={config.contrast.value}, size={config.volume_size}"
        )

        payload: Dict[str, Any] = {
            "body_region": config.body_region,
            "contrast": config.contrast.value,
            "resolution_mm": config.resolution_mm,
            "volume_size": config.volume_size,
            "skull_stripped": config.skull_stripped,
            "include_pathology": config.include_pathology,
        }
        if config.pathology_type:
            payload["pathology_type"] = config.pathology_type
        if config.seed is not None:
            payload["seed"] = config.seed

        result = self._invoke_or_mock(
            endpoint="/v1/generate",
            payload=payload,
            timeout=600,
            config=config,
        )

        # If real NIM returned raw dict, parse into SyntheticMRResult
        if isinstance(result, dict) and not isinstance(result, SyntheticMRResult):
            return SyntheticMRResult(
                volume_path=result.get("volume_path"),
                mask_path=result.get("mask_path"),
                body_region=config.body_region,
                contrast=config.contrast.value,
                resolution_mm=config.resolution_mm,
                volume_shape=result.get("volume_shape", config.volume_size),
                generation_time_ms=result.get("generation_time_ms", 0.0),
                model_name="nv-generate-mr",
                is_mock=False,
            )

        return result

    def generate_multi_contrast(
        self,
        body_region: str = "brain",
        contrasts: Optional[List[MRContrastType]] = None,
        resolution_mm: float = 1.0,
        volume_size: Optional[List[int]] = None,
        skull_stripped: bool = False,
        seed: Optional[int] = None,
    ) -> List[SyntheticMRResult]:
        """Generate synthetic MR volumes for multiple contrast types.

        Produces one volume per contrast, all sharing the same body
        region and spatial configuration. Useful for generating
        co-registered multi-contrast training sets.

        Args:
            body_region: Anatomical region for all volumes.
            contrasts: List of contrast types. Defaults to [T1W, T2W, FLAIR].
            resolution_mm: Voxel spacing in mm.
            volume_size: Volume dimensions [W, H, D].
            skull_stripped: Remove skull (brain only).
            seed: Base seed; each contrast gets seed+i.

        Returns:
            List of SyntheticMRResult, one per contrast.
        """
        if contrasts is None:
            contrasts = [MRContrastType.T1W, MRContrastType.T2W, MRContrastType.FLAIR]
        if volume_size is None:
            volume_size = [256, 256, 256]

        logger.info(
            f"Generating multi-contrast MR: region={body_region}, "
            f"contrasts={[c.value for c in contrasts]}"
        )

        results: List[SyntheticMRResult] = []
        for i, contrast in enumerate(contrasts):
            config = SyntheticMRConfig(
                body_region=body_region,
                contrast=contrast,
                resolution_mm=resolution_mm,
                volume_size=volume_size,
                skull_stripped=skull_stripped,
                seed=(seed + i) if seed is not None else None,
            )
            result = self.generate(config)
            results.append(result)

        logger.info(
            f"Multi-contrast generation complete: {len(results)} volumes "
            f"({', '.join(c.value for c in contrasts)})"
        )
        return results

    def _mock_response(self, **kwargs) -> SyntheticMRResult:
        """Return a realistic mock SyntheticMRResult.

        Simulates the metadata that NV-Generate-MR would return
        without actually generating a volume. Generation times scale
        with volume size (2000-15000ms range).
        """
        config: SyntheticMRConfig = kwargs.get("config", SyntheticMRConfig())

        # Use seed for reproducibility if provided
        rng = random.Random(config.seed) if config.seed is not None else random

        # Simulate realistic generation times based on volume size
        total_voxels = 1
        for d in config.volume_size:
            total_voxels *= d
        ref_voxels = 256 * 256 * 256
        base_time_ms = (total_voxels / ref_voxels) * 6000.0  # MR slightly slower
        gen_time_ms = round(base_time_ms * rng.uniform(0.8, 1.2), 1)
        gen_time_ms = max(2000.0, min(gen_time_ms, 15000.0))

        contrast_str = config.contrast.value
        pathology_suffix = ""
        if config.include_pathology and config.pathology_type:
            pathology_suffix = f"_{config.pathology_type}"

        skull_suffix = "_stripped" if config.skull_stripped else ""
        seed_suffix = f"_s{config.seed}" if config.seed is not None else ""

        logger.info(
            f"Mock NV-Generate-MR: {config.body_region} {contrast_str} "
            f"{'x'.join(str(d) for d in config.volume_size)} "
            f"({gen_time_ms:.0f}ms simulated)"
        )

        return SyntheticMRResult(
            volume_path=(
                f"/tmp/nv_generate_mr_mock/"
                f"{config.body_region}_{contrast_str}"
                f"{skull_suffix}{pathology_suffix}{seed_suffix}_synthetic.nii.gz"
            ),
            mask_path=(
                f"/tmp/nv_generate_mr_mock/"
                f"{config.body_region}_{contrast_str}"
                f"{skull_suffix}{pathology_suffix}{seed_suffix}_mask.nii.gz"
            ),
            body_region=config.body_region,
            contrast=contrast_str,
            resolution_mm=config.resolution_mm,
            volume_shape=list(config.volume_size),
            generation_time_ms=gen_time_ms,
            model_name="nv-generate-mr-mock",
            is_mock=True,
        )
