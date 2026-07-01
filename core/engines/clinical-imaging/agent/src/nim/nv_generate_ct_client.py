"""NV-Generate-CT client for synthetic CT volume generation.

Generates high-fidelity synthetic CT volumes with paired segmentation
masks and anatomical annotations up to 512x512x768 voxels.
127 anatomical classes. Rectified flow scheduler (33x faster than DDPM).

License: Check per release (model weights on HuggingFace).
"""

import random
import time
from typing import Any, Dict, List, Optional

from loguru import logger
from pydantic import BaseModel, Field

from .base import BaseNIMClient


# ═══════════════════════════════════════════════════════════════════════
# CONFIG & RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════

VALID_BODY_REGIONS: List[str] = [
    "chest", "abdomen", "head", "pelvis", "whole_body",
]

VALID_PATHOLOGY_TYPES: List[str] = [
    "nodule", "mass", "hemorrhage",
]


class SyntheticGenerationConfig(BaseModel):
    """Configuration for NV-Generate-CT synthetic volume generation."""

    body_region: str = Field("chest", description="Target body region")
    resolution_mm: float = Field(1.0, ge=0.5, le=5.0, description="Voxel spacing in mm")
    volume_size: List[int] = Field(
        default_factory=lambda: [256, 256, 256],
        description="Volume dimensions [W, H, D] up to 512x512x768",
    )
    num_classes: int = Field(127, ge=1, le=127, description="Anatomical classes for paired mask")
    include_pathology: bool = Field(False, description="Inject synthetic pathology")
    pathology_type: Optional[str] = Field(None, description="Type: nodule, mass, hemorrhage")
    seed: Optional[int] = Field(None, description="Random seed for reproducibility")


class SyntheticVolumeResult(BaseModel):
    """Result from NV-Generate-CT synthetic volume generation."""

    volume_path: Optional[str] = None
    mask_path: Optional[str] = None
    body_region: str = ""
    resolution_mm: float = 1.0
    volume_shape: List[int] = Field(default_factory=lambda: [256, 256, 256])
    classes_generated: int = 0
    generation_time_ms: float = 0.0
    model_name: str = "nv-generate-ct"
    is_mock: bool = False


# ═══════════════════════════════════════════════════════════════════════
# CLIENT
# ═══════════════════════════════════════════════════════════════════════


class NVGenerateCTClient(BaseNIMClient):
    """Client for NV-Generate-CT synthetic CT volume generation.

    NV-Generate-CT produces high-fidelity synthetic CT volumes with
    paired segmentation masks across 127 anatomical classes. Uses a
    rectified flow scheduler that is ~33x faster than DDPM. Useful for:
      - Training data augmentation
      - Privacy-preserving synthetic datasets
      - Rare pathology simulation
      - Model validation and benchmarking
    """

    def __init__(self, base_url: str, mock_enabled: bool = True):
        super().__init__(base_url, service_name="nv-generate-ct", mock_enabled=mock_enabled)

    def generate(
        self,
        config: Optional[SyntheticGenerationConfig] = None,
    ) -> SyntheticVolumeResult:
        """Generate a synthetic CT volume with paired segmentation mask.

        Args:
            config: Generation configuration. Uses defaults if None.

        Returns:
            SyntheticVolumeResult with paths and generation metadata.
        """
        if config is None:
            config = SyntheticGenerationConfig()

        logger.info(
            f"Generating synthetic CT: region={config.body_region}, "
            f"size={config.volume_size}, resolution={config.resolution_mm}mm"
        )

        payload: Dict[str, Any] = {
            "body_region": config.body_region,
            "resolution_mm": config.resolution_mm,
            "volume_size": config.volume_size,
            "num_classes": config.num_classes,
            "include_pathology": config.include_pathology,
        }
        if config.pathology_type:
            payload["pathology_type"] = config.pathology_type
        if config.seed is not None:
            payload["seed"] = config.seed

        result = self._invoke_or_mock(
            endpoint="/v1/generate",
            payload=payload,
            timeout=600,  # Generation can take several minutes
            config=config,
        )

        # If real NIM returned raw dict, parse into SyntheticVolumeResult
        if isinstance(result, dict) and not isinstance(result, SyntheticVolumeResult):
            return SyntheticVolumeResult(
                volume_path=result.get("volume_path"),
                mask_path=result.get("mask_path"),
                body_region=config.body_region,
                resolution_mm=config.resolution_mm,
                volume_shape=result.get("volume_shape", config.volume_size),
                classes_generated=result.get("classes_generated", config.num_classes),
                generation_time_ms=result.get("generation_time_ms", 0.0),
                model_name="nv-generate-ct",
                is_mock=False,
            )

        return result

    def generate_training_batch(
        self,
        config: Optional[SyntheticGenerationConfig] = None,
        count: int = 5,
    ) -> List[SyntheticVolumeResult]:
        """Generate multiple synthetic CT volumes for training data augmentation.

        Args:
            config: Generation configuration applied to each volume.
            count: Number of volumes to generate (1-50).

        Returns:
            List of SyntheticVolumeResult, one per generated volume.
        """
        if config is None:
            config = SyntheticGenerationConfig()

        count = max(1, min(count, 50))
        logger.info(
            f"Generating batch of {count} synthetic CT volumes: "
            f"region={config.body_region}, size={config.volume_size}"
        )

        results: List[SyntheticVolumeResult] = []
        for i in range(count):
            # Vary seed per volume if a base seed is provided
            batch_config = config.model_copy()
            if config.seed is not None:
                batch_config.seed = config.seed + i
            else:
                batch_config.seed = random.randint(0, 2**31)

            result = self.generate(batch_config)
            results.append(result)

        logger.info(
            f"Batch generation complete: {len(results)}/{count} volumes generated"
        )
        return results

    def _mock_response(self, **kwargs) -> SyntheticVolumeResult:
        """Return a realistic mock SyntheticVolumeResult.

        Simulates the metadata that NV-Generate-CT would return
        without actually generating a volume. Generation times scale
        with volume size (2000-15000ms range).
        """
        config: SyntheticGenerationConfig = kwargs.get("config", SyntheticGenerationConfig())

        # Use seed for reproducibility if provided
        rng = random.Random(config.seed) if config.seed is not None else random

        # Simulate realistic generation times based on volume size
        total_voxels = 1
        for d in config.volume_size:
            total_voxels *= d
        ref_voxels = 256 * 256 * 256  # reference size
        base_time_ms = (total_voxels / ref_voxels) * 5000.0
        gen_time_ms = round(base_time_ms * rng.uniform(0.8, 1.2), 1)
        # Clamp to realistic range
        gen_time_ms = max(2000.0, min(gen_time_ms, 15000.0))

        pathology_suffix = ""
        if config.include_pathology and config.pathology_type:
            pathology_suffix = f"_{config.pathology_type}"

        seed_suffix = f"_s{config.seed}" if config.seed is not None else ""

        logger.info(
            f"Mock NV-Generate-CT: {config.body_region} "
            f"{'x'.join(str(d) for d in config.volume_size)} "
            f"({gen_time_ms:.0f}ms simulated)"
        )

        return SyntheticVolumeResult(
            volume_path=(
                f"/tmp/nv_generate_ct_mock/"
                f"{config.body_region}{pathology_suffix}{seed_suffix}_synthetic.nii.gz"
            ),
            mask_path=(
                f"/tmp/nv_generate_ct_mock/"
                f"{config.body_region}{pathology_suffix}{seed_suffix}_mask.nii.gz"
            ),
            body_region=config.body_region,
            resolution_mm=config.resolution_mm,
            volume_shape=list(config.volume_size),
            classes_generated=config.num_classes,
            generation_time_ms=gen_time_ms,
            model_name="nv-generate-ct-mock",
            is_mock=True,
        )
