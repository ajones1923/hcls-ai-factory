"""MAISI NIM client for synthetic CT volume generation.

MAISI (Medical AI for Synthetic Imaging) generates high-resolution
synthetic 3D CT volumes with controllable anatomy, useful for training
data augmentation and model validation.
"""

import random
import time
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import SyntheticCTResult

from .base import BaseNIMClient


class MAISIClient(BaseNIMClient):
    """Client for NVIDIA MAISI NIM synthetic CT generation service.

    MAISI generates photorealistic 3D CT scans at up to 512x512x512
    resolution with control over body region, anatomy classes, and
    pathology presence. Useful for:
      - Training data augmentation
      - Privacy-preserving synthetic datasets
      - Rare pathology simulation
      - Model validation and benchmarking
    """

    def __init__(self, base_url: str, mock_enabled: bool = True):
        super().__init__(base_url, service_name="maisi", mock_enabled=mock_enabled)

    def generate(
        self,
        body_region: str = "chest",
        num_classes: int = 127,
        resolution: str = "512x512x512",
    ) -> SyntheticCTResult:
        """Generate a synthetic CT volume.

        Args:
            body_region: Anatomical region to generate. One of:
                "chest", "abdomen", "pelvis", "head", "full_body".
            num_classes: Number of anatomical label classes in the
                generated segmentation mask (up to 127).
            resolution: Volume resolution as "WxHxD" string.

        Returns:
            SyntheticCTResult with paths to generated volume and
            segmentation mask, generation metadata, and timing.
        """
        logger.info(
            f"Generating synthetic CT: region={body_region}, "
            f"resolution={resolution}, classes={num_classes}"
        )

        # Parse resolution
        dims = [int(d) for d in resolution.split("x")]
        if len(dims) != 3:
            raise ValueError(
                f"Resolution must be 'WxHxD' format, got: {resolution}"
            )

        payload: Dict[str, Any] = {
            "body_region": body_region,
            "num_output_channels": num_classes,
            "output_size": dims,
        }

        result = self._invoke_or_mock(
            endpoint="/maisi/generate",
            payload=payload,
            timeout=600,  # Generation can take several minutes
            body_region=body_region,
            num_classes=num_classes,
            resolution=resolution,
        )

        # If real NIM returned raw dict, parse into SyntheticCTResult
        if isinstance(result, dict) and not isinstance(result, SyntheticCTResult):
            return SyntheticCTResult(
                volume_path=result.get("volume_path", ""),
                segmentation_mask_path=result.get("segmentation_mask_path"),
                body_region=body_region,
                resolution=resolution,
                num_classes_annotated=result.get("num_classes", num_classes),
                generation_time_ms=result.get("generation_time_ms", 0.0),
                voxel_spacing_mm=result.get("voxel_spacing_mm", [1.0, 1.0, 1.0]),
                model_name="maisi",
                is_mock=False,
            )

        return result

    def _mock_response(self, **kwargs) -> SyntheticCTResult:
        """Return a realistic mock SyntheticCTResult.

        Simulates the metadata that MAISI would return without
        actually generating a volume.
        """
        body_region = kwargs.get("body_region", "chest")
        resolution = kwargs.get("resolution", "512x512x512")
        num_classes = kwargs.get("num_classes", 127)

        # Simulate realistic generation times based on resolution
        dims = [int(d) for d in resolution.split("x")]
        total_voxels = dims[0] * dims[1] * dims[2]
        base_time_ms = (total_voxels / (512 * 512 * 512)) * 45000.0  # ~45s for 512^3
        gen_time_ms = round(base_time_ms * random.uniform(0.8, 1.2), 1)

        # Standard voxel spacings by body region
        voxel_spacings = {
            "chest": [1.0, 1.0, 1.5],
            "abdomen": [0.8, 0.8, 1.5],
            "pelvis": [1.0, 1.0, 1.5],
            "head": [0.5, 0.5, 0.5],
            "full_body": [1.5, 1.5, 2.0],
        }

        logger.info(
            f"Mock MAISI generation: {body_region} at {resolution} "
            f"({gen_time_ms:.0f}ms simulated)"
        )

        return SyntheticCTResult(
            volume_path=f"/tmp/maisi_mock/{body_region}_synthetic.nii.gz",
            segmentation_mask_path=f"/tmp/maisi_mock/{body_region}_labels.nii.gz",
            body_region=body_region,
            resolution=resolution,
            num_classes_annotated=num_classes,
            generation_time_ms=gen_time_ms,
            voxel_spacing_mm=voxel_spacings.get(body_region, [1.0, 1.0, 1.0]),
            model_name="maisi-mock",
            is_mock=True,
        )
