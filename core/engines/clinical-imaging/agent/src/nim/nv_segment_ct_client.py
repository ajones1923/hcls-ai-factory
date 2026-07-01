"""NV-Segment-CT NIM client for 3D CT segmentation with tumor support.

NV-Segment-CT (nvidia/NV-Segment-CT) extends VISTA-3D with 132 anatomical
classes including 7 tumor types for comprehensive CT analysis. Licensed
under NVIDIA Open Model License (commercial use OK). CT volumes only.

HuggingFace weights: nvidia/NV-Segment-CT
License: NVIDIA Open Model License (commercial OK)
Classes: 132 (125 anatomical structures + 7 tumor types)
Modality: CT only
"""

import random
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import SegmentationResult

from .base import BaseNIMClient

# 7 tumor classes added by NV-Segment-CT beyond VISTA-3D
TUMOR_CLASSES: List[str] = [
    "liver_tumor",
    "lung_tumor",
    "kidney_tumor",
    "pancreas_tumor",
    "colon_tumor",
    "hepatic_vessel_tumor",
    "adrenal_gland_tumor",
]

# 132 anatomical classes: 125 anatomical structures + 7 tumor types
NV_SEGMENT_CT_CLASSES: List[str] = [
    # ── Organs ──
    "liver", "right_kidney", "spleen", "pancreas", "aorta",
    "inferior_vena_cava", "right_adrenal_gland", "left_adrenal_gland",
    "gallbladder", "esophagus", "stomach", "duodenum",
    "left_kidney", "urinary_bladder", "prostate", "rectum",
    "small_bowel", "colon", "lung_upper_lobe_left", "lung_lower_lobe_left",
    "lung_upper_lobe_right", "lung_middle_lobe_right", "lung_lower_lobe_right",
    "trachea", "heart", "pulmonary_artery", "brain",
    # ── Vasculature ──
    "iliac_artery_left", "iliac_artery_right", "iliac_vena_left",
    "iliac_vena_right", "portal_vein", "hepatic_vein",
    "celiac_trunk", "superior_mesenteric_artery", "inferior_mesenteric_artery",
    # ── Spine ──
    "vertebrae_L5", "vertebrae_L4", "vertebrae_L3", "vertebrae_L2", "vertebrae_L1",
    "vertebrae_T12", "vertebrae_T11", "vertebrae_T10", "vertebrae_T9",
    "vertebrae_T8", "vertebrae_T7", "vertebrae_T6", "vertebrae_T5",
    "vertebrae_T4", "vertebrae_T3", "vertebrae_T2", "vertebrae_T1",
    "vertebrae_C7", "vertebrae_C6", "vertebrae_C5", "vertebrae_C4",
    "vertebrae_C3", "vertebrae_C2", "vertebrae_C1",
    # ── Ribs ──
    "rib_left_1", "rib_left_2", "rib_left_3", "rib_left_4", "rib_left_5",
    "rib_left_6", "rib_left_7", "rib_left_8", "rib_left_9", "rib_left_10",
    "rib_left_11", "rib_left_12",
    "rib_right_1", "rib_right_2", "rib_right_3", "rib_right_4", "rib_right_5",
    "rib_right_6", "rib_right_7", "rib_right_8", "rib_right_9", "rib_right_10",
    "rib_right_11", "rib_right_12",
    # ── Bones & Musculature ──
    "humerus_left", "humerus_right", "scapula_left", "scapula_right",
    "clavicula_left", "clavicula_right", "femur_left", "femur_right",
    "hip_left", "hip_right", "sacrum",
    "gluteus_maximus_left", "gluteus_maximus_right",
    "gluteus_medius_left", "gluteus_medius_right",
    "gluteus_minimus_left", "gluteus_minimus_right",
    "autochthon_left", "autochthon_right",
    "iliopsoas_left", "iliopsoas_right",
    "sternum",
    # ── Head & Neck ──
    "thyroid_gland",
    # ── CNS ──
    "spinal_cord",
    # ── Cardiac chambers & great vessels ──
    "left_atrium", "right_atrium", "left_ventricle", "right_ventricle",
    "myocardium", "ascending_aorta", "descending_aorta", "aortic_arch",
    "brachiocephalic_trunk", "subclavian_artery_right", "subclavian_artery_left",
    "common_carotid_artery_right", "common_carotid_artery_left",
    "brachiocephalic_vein_left", "brachiocephalic_vein_right",
    "superior_vena_cava", "skull",
    # ── Tumor classes (7) ──
    "liver_tumor",
    "lung_tumor",
    "kidney_tumor",
    "pancreas_tumor",
    "colon_tumor",
    "hepatic_vessel_tumor",
    "adrenal_gland_tumor",
]


class NVSegmentCTClient(BaseNIMClient):
    """Client for NVIDIA NV-Segment-CT NIM segmentation service.

    NV-Segment-CT provides automatic and interactive segmentation of
    132 anatomical structures from CT volumes, including 7 tumor types
    for oncology workflows. Supports both full-volume inference and
    point-click guided segmentation.

    Model: nvidia/NV-Segment-CT (HuggingFace)
    License: NVIDIA Open Model License (commercial use OK)
    Classes: 132 (125 anatomical structures + 7 tumor types)
    Modality: CT only
    """

    def __init__(self, base_url: str, mock_enabled: bool = True):
        super().__init__(base_url, service_name="nv-segment-ct", mock_enabled=mock_enabled)

    def segment(
        self,
        ct_volume_path: str,
        classes: Optional[List[str]] = None,
    ) -> SegmentationResult:
        """Run automatic segmentation on a CT volume.

        Args:
            ct_volume_path: Path to the NIfTI (.nii.gz) CT volume.
            classes: Optional list of anatomical classes to segment.
                     If None, all 132 classes are segmented.

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
            endpoint="/v1/inference",
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
                model_name="nv-segment-ct",
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
            endpoint="/v1/inference",
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
                model_name="nv-segment-ct-interactive",
                is_mock=False,
            )

        return result

    def segment_tumors(
        self,
        ct_volume_path: str,
    ) -> SegmentationResult:
        """Convenience method to segment only the 7 tumor classes.

        Runs segmentation restricted to liver_tumor, lung_tumor,
        kidney_tumor, pancreas_tumor, colon_tumor, hepatic_vessel_tumor,
        and adrenal_gland_tumor.

        Args:
            ct_volume_path: Path to the NIfTI (.nii.gz) CT volume.

        Returns:
            SegmentationResult with only tumor classes detected.
        """
        logger.info(f"Tumor-only segmentation: {ct_volume_path}")
        return self.segment(ct_volume_path, classes=list(TUMOR_CLASSES))

    def get_supported_classes(self) -> List[str]:
        """Return the list of 132 supported anatomical structures."""
        return list(NV_SEGMENT_CT_CLASSES)

    def _mock_response(self, **kwargs) -> SegmentationResult:
        """Return a realistic mock SegmentationResult.

        Generates 5-10 detected anatomical structures with plausible
        volume measurements in cubic centimeters. Tumor volumes are
        realistically small (0.5-25 cm^3) compared to organ volumes.
        """
        classes = kwargs.get("classes")

        if classes:
            # When specific classes requested, return a subset
            pool = classes
            num_classes = min(random.randint(1, max(len(pool), 1)), len(pool))
        else:
            pool = NV_SEGMENT_CT_CLASSES
            num_classes = random.randint(5, 10)

        selected = random.sample(pool, min(num_classes, len(pool)))

        # Plausible organ volumes in cm^3
        volume_ranges: Dict[str, tuple] = {
            # Organs
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
            # Tumors — realistically small compared to organs
            "liver_tumor": (0.5, 15.0),
            "lung_tumor": (0.3, 12.0),
            "kidney_tumor": (0.5, 10.0),
            "pancreas_tumor": (0.8, 8.0),
            "colon_tumor": (1.0, 20.0),
            "hepatic_vessel_tumor": (0.3, 5.0),
            "adrenal_gland_tumor": (0.5, 6.0),
        }

        volumes: Dict[str, float] = {}
        for cls in selected:
            if cls in volume_ranges:
                lo, hi = volume_ranges[cls]
                volumes[cls] = round(random.uniform(lo, hi), 1)
            else:
                volumes[cls] = round(random.uniform(5.0, 150.0), 1)

        logger.info(
            f"Mock NV-Segment-CT segmentation: {len(selected)} classes detected"
        )

        return SegmentationResult(
            classes_detected=selected,
            volumes=volumes,
            inference_time_ms=round(random.uniform(2000.0, 8000.0), 1),
            segmentation_mask_path=None,
            model_name="nv-segment-ct-mock",
            is_mock=True,
        )
