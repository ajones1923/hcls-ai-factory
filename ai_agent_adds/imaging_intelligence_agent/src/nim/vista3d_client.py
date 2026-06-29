"""VISTA-3D NIM client for 3D medical image segmentation.

VISTA-3D supports 127+ anatomical structures from CT volumes with
both automatic and interactive (point-prompt) segmentation modes.
"""

import base64
import random
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import SegmentationResult

from .base import BaseNIMClient

# 127+ anatomical classes supported by VISTA-3D
VISTA3D_CLASSES: List[str] = [
    "liver", "right_kidney", "spleen", "pancreas", "aorta",
    "inferior_vena_cava", "right_adrenal_gland", "left_adrenal_gland",
    "gallbladder", "esophagus", "stomach", "duodenum",
    "left_kidney", "urinary_bladder", "prostate", "rectum",
    "small_bowel", "colon", "lung_upper_lobe_left", "lung_lower_lobe_left",
    "lung_upper_lobe_right", "lung_middle_lobe_right", "lung_lower_lobe_right",
    "trachea", "heart", "pulmonary_artery", "brain",
    "iliac_artery_left", "iliac_artery_right", "iliac_vena_left",
    "iliac_vena_right", "portal_vein", "hepatic_vein",
    "celiac_trunk", "superior_mesenteric_artery", "inferior_mesenteric_artery",
    "vertebrae_L5", "vertebrae_L4", "vertebrae_L3", "vertebrae_L2", "vertebrae_L1",
    "vertebrae_T12", "vertebrae_T11", "vertebrae_T10", "vertebrae_T9",
    "vertebrae_T8", "vertebrae_T7", "vertebrae_T6", "vertebrae_T5",
    "vertebrae_T4", "vertebrae_T3", "vertebrae_T2", "vertebrae_T1",
    "vertebrae_C7", "vertebrae_C6", "vertebrae_C5", "vertebrae_C4",
    "vertebrae_C3", "vertebrae_C2", "vertebrae_C1",
    "rib_left_1", "rib_left_2", "rib_left_3", "rib_left_4", "rib_left_5",
    "rib_left_6", "rib_left_7", "rib_left_8", "rib_left_9", "rib_left_10",
    "rib_left_11", "rib_left_12",
    "rib_right_1", "rib_right_2", "rib_right_3", "rib_right_4", "rib_right_5",
    "rib_right_6", "rib_right_7", "rib_right_8", "rib_right_9", "rib_right_10",
    "rib_right_11", "rib_right_12",
    "humerus_left", "humerus_right", "scapula_left", "scapula_right",
    "clavicula_left", "clavicula_right", "femur_left", "femur_right",
    "hip_left", "hip_right", "sacrum",
    "gluteus_maximus_left", "gluteus_maximus_right",
    "gluteus_medius_left", "gluteus_medius_right",
    "gluteus_minimus_left", "gluteus_minimus_right",
    "autochthon_left", "autochthon_right",
    "iliopsoas_left", "iliopsoas_right",
    "sternum", "costal_cartilages",
    "thyroid_gland", "face",
    "spinal_cord", "pericardium",
    "adrenal_gland_left", "adrenal_gland_right",
    "left_atrium", "right_atrium", "left_ventricle", "right_ventricle",
    "myocardium", "ascending_aorta", "descending_aorta", "aortic_arch",
    "brachiocephalic_trunk", "subclavian_artery_right", "subclavian_artery_left",
    "common_carotid_artery_right", "common_carotid_artery_left",
    "brachiocephalic_vein_left", "brachiocephalic_vein_right",
    "superior_vena_cava", "skull", "mandible",
]


class VISTA3DClient(BaseNIMClient):
    """Client for NVIDIA VISTA-3D NIM segmentation service.

    VISTA-3D provides automatic and interactive segmentation of 127+
    anatomical structures from CT volumes. Supports both full-volume
    inference and point-click guided segmentation.
    """

    def __init__(self, base_url: str, mock_enabled: bool = True):
        super().__init__(base_url, service_name="vista3d", mock_enabled=mock_enabled)

    def segment(
        self,
        ct_volume_path: str,
        classes: Optional[List[str]] = None,
    ) -> SegmentationResult:
        """Run automatic segmentation on a CT volume.

        Args:
            ct_volume_path: Path to the NIfTI (.nii.gz) CT volume.
            classes: Optional list of anatomical classes to segment.
                     If None, all 127+ classes are segmented.

        Returns:
            SegmentationResult with detected classes, volumes, and timing.
        """
        logger.info(f"Segmenting CT volume: {ct_volume_path}")

        payload: Dict[str, Any] = {
            "image": ct_volume_path,
        }
        if classes:
            payload["classes"] = classes

        result = self._invoke_or_mock(
            endpoint="/vista3d/inference",
            payload=payload,
            timeout=300,
            ct_volume_path=ct_volume_path,
            classes=classes,
        )

        # If real NIM returned raw dict, parse into SegmentationResult
        if isinstance(result, dict) and not isinstance(result, SegmentationResult):
            return SegmentationResult(
                classes_detected=result.get("classes_detected", []),
                volumes=result.get("volumes", {}),
                inference_time_ms=result.get("inference_time_ms", 0.0),
                segmentation_mask_path=result.get("segmentation_mask_path"),
                model_name="vista3d",
                is_mock=False,
            )

        return result

    def segment_interactive(
        self,
        ct_volume_path: str,
        point_prompts: List[Dict],
    ) -> SegmentationResult:
        """Run interactive point-click guided segmentation.

        Args:
            ct_volume_path: Path to the NIfTI (.nii.gz) CT volume.
            point_prompts: List of point prompt dicts, each containing:
                - "point": [x, y, z] coordinates in voxel space
                - "label": 1 for foreground, 0 for background

        Returns:
            SegmentationResult with detected classes and volumes.
        """
        logger.info(
            f"Interactive segmentation with {len(point_prompts)} prompts: "
            f"{ct_volume_path}"
        )

        payload: Dict[str, Any] = {
            "image": ct_volume_path,
            "prompts": point_prompts,
        }

        result = self._invoke_or_mock(
            endpoint="/vista3d/inference",
            payload=payload,
            timeout=300,
            ct_volume_path=ct_volume_path,
            point_prompts=point_prompts,
        )

        if isinstance(result, dict) and not isinstance(result, SegmentationResult):
            return SegmentationResult(
                classes_detected=result.get("classes_detected", []),
                volumes=result.get("volumes", {}),
                inference_time_ms=result.get("inference_time_ms", 0.0),
                segmentation_mask_path=result.get("segmentation_mask_path"),
                model_name="vista3d-interactive",
                is_mock=False,
            )

        return result

    def get_supported_classes(self) -> List[str]:
        """Return the list of 127+ supported anatomical structures."""
        return list(VISTA3D_CLASSES)

    def _mock_response(self, **kwargs) -> SegmentationResult:
        """Return a realistic mock SegmentationResult.

        Generates 5-10 detected anatomical structures with plausible
        volume measurements in cubic centimeters.
        """
        num_classes = random.randint(5, 10)
        selected = random.sample(VISTA3D_CLASSES, min(num_classes, len(VISTA3D_CLASSES)))

        # Plausible organ volumes in cm^3
        volume_ranges: Dict[str, tuple] = {
            "liver": (1200.0, 1800.0),
            "spleen": (100.0, 300.0),
            "right_kidney": (120.0, 200.0),
            "left_kidney": (120.0, 200.0),
            "pancreas": (60.0, 100.0),
            "heart": (250.0, 400.0),
            "aorta": (30.0, 60.0),
            "gallbladder": (30.0, 60.0),
            "stomach": (200.0, 400.0),
            "lung_upper_lobe_left": (400.0, 700.0),
            "lung_lower_lobe_left": (500.0, 800.0),
            "lung_upper_lobe_right": (450.0, 750.0),
            "lung_middle_lobe_right": (200.0, 400.0),
            "lung_lower_lobe_right": (550.0, 850.0),
        }

        volumes: Dict[str, float] = {}
        for cls in selected:
            if cls in volume_ranges:
                lo, hi = volume_ranges[cls]
                volumes[cls] = round(random.uniform(lo, hi), 1)
            else:
                volumes[cls] = round(random.uniform(5.0, 150.0), 1)

        logger.info(
            f"Mock VISTA-3D segmentation: {num_classes} classes detected"
        )

        return SegmentationResult(
            classes_detected=selected,
            volumes=volumes,
            inference_time_ms=round(random.uniform(2000.0, 8000.0), 1),
            segmentation_mask_path=None,
            model_name="vista3d-mock",
            is_mock=True,
        )
