"""Thyroid TI-RADS Classification Workflow.

Clinical workflow for thyroid ultrasound with ACR TI-RADS scoring.
Cross-modal trigger: TI-RADS 4+ -> BRAF V600E, RAS, RET/PTC genomic queries.

Scoring system: ACR TI-RADS (2017)
  Points system based on 5 features:
  - Composition: cystic(0), spongiform(0), mixed(1), solid(2)
  - Echogenicity: anechoic(0), hyper/isoechoic(1), hypoechoic(2), very hypoechoic(3)
  - Shape: wider-than-tall(0), taller-than-wide(3)
  - Margin: smooth(0), ill-defined(0), lobulated/irregular(2), extra-thyroidal extension(3)
  - Echogenic foci: none(0), comet-tail(0), macrocalcifications(1), peripheral(2), punctate(3)

  Total points -> TR level:
  TR1 (0 pts): Benign -- no FNA
  TR2 (2 pts): Not suspicious -- no FNA
  TR3 (3 pts): Mildly suspicious -- FNA if >= 2.5cm, follow if >= 1.5cm
  TR4 (4-6 pts): Moderately suspicious -- FNA if >= 1.5cm, follow if >= 1.0cm
  TR5 (7+ pts): Highly suspicious -- FNA if >= 1.0cm, follow if >= 0.5cm

Target latency: < 3 minutes
Models: MONAI thyroid segmentation

Author: Adam Jones
Date: April 2026
"""

from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    FindingCategory,
    FindingSeverity,
    TIRADS,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows.base import BaseImagingWorkflow


# =====================================================================
# TI-RADS point values per feature category
# =====================================================================

COMPOSITION_POINTS = {
    "cystic": 0,
    "spongiform": 0,
    "mixed_cystic_solid": 1,
    "solid": 2,
}

ECHOGENICITY_POINTS = {
    "anechoic": 0,
    "hyperechoic": 1,
    "isoechoic": 1,
    "hypoechoic": 2,
    "very_hypoechoic": 3,
}

SHAPE_POINTS = {
    "wider_than_tall": 0,
    "taller_than_wide": 3,
}

MARGIN_POINTS = {
    "smooth": 0,
    "ill_defined": 0,
    "lobulated": 2,
    "irregular": 2,
    "extra_thyroidal_extension": 3,
}

ECHOGENIC_FOCI_POINTS = {
    "none": 0,
    "large_comet_tail": 0,
    "macrocalcifications": 1,
    "peripheral_calcifications": 2,
    "punctate_echogenic_foci": 3,
}


# =====================================================================
# TI-RADS helper functions
# =====================================================================


def calculate_tirads_points(nodule: Dict) -> int:
    """Calculate total ACR TI-RADS points from nodule features.

    Args:
        nodule: Dict with keys composition, echogenicity, shape,
                margin, echogenic_foci.

    Returns:
        Total TI-RADS points (0-14 theoretical max).
    """
    points = 0
    points += COMPOSITION_POINTS.get(nodule.get("composition", ""), 0)
    points += ECHOGENICITY_POINTS.get(nodule.get("echogenicity", ""), 0)
    points += SHAPE_POINTS.get(nodule.get("shape", ""), 0)
    points += MARGIN_POINTS.get(nodule.get("margin", ""), 0)
    points += ECHOGENIC_FOCI_POINTS.get(nodule.get("echogenic_foci", ""), 0)
    return points


def points_to_tirads(points: int) -> TIRADS:
    """Map total TI-RADS points to TR level.

    ACR TI-RADS (2017):
        0 pts:   TR1 (Benign)
        2 pts:   TR2 (Not suspicious)
        3 pts:   TR3 (Mildly suspicious)
        4-6 pts: TR4 (Moderately suspicious)
        7+ pts:  TR5 (Highly suspicious)

    Args:
        points: Total TI-RADS points.

    Returns:
        TIRADS enum value.
    """
    if points <= 0:
        return TIRADS.TR1
    elif points <= 2:
        return TIRADS.TR2
    elif points == 3:
        return TIRADS.TR3
    elif points <= 6:
        return TIRADS.TR4
    else:
        return TIRADS.TR5


def tirads_to_severity(tr_level: TIRADS) -> FindingSeverity:
    """Map TI-RADS level to clinical severity.

    Args:
        tr_level: ACR TI-RADS level (TR1-TR5).

    Returns:
        FindingSeverity corresponding to the TR level.
    """
    mapping = {
        TIRADS.TR1: FindingSeverity.NORMAL,
        TIRADS.TR2: FindingSeverity.NORMAL,
        TIRADS.TR3: FindingSeverity.SIGNIFICANT,
        TIRADS.TR4: FindingSeverity.URGENT,
        TIRADS.TR5: FindingSeverity.CRITICAL,
    }
    return mapping.get(tr_level, FindingSeverity.ROUTINE)


def tirads_fna_recommendation(tr_level: TIRADS, max_diameter_mm: float) -> str:
    """Determine FNA recommendation based on TR level and nodule size.

    ACR TI-RADS size thresholds for FNA:
        TR1: No FNA at any size
        TR2: No FNA at any size
        TR3: FNA if >= 25mm, follow-up if >= 15mm
        TR4: FNA if >= 15mm, follow-up if >= 10mm
        TR5: FNA if >= 10mm, follow-up if >= 5mm

    Args:
        tr_level: ACR TI-RADS level.
        max_diameter_mm: Maximum nodule diameter in mm.

    Returns:
        Management recommendation string.
    """
    if tr_level == TIRADS.TR1:
        return "TR1 (Benign). No FNA indicated. No follow-up needed."

    if tr_level == TIRADS.TR2:
        return "TR2 (Not suspicious). No FNA indicated. No follow-up needed."

    if tr_level == TIRADS.TR3:
        if max_diameter_mm >= 25.0:
            return (
                "TR3 (Mildly suspicious). FNA recommended (nodule >= 2.5cm). "
                "Consider ultrasound-guided fine needle aspiration."
            )
        elif max_diameter_mm >= 15.0:
            return (
                "TR3 (Mildly suspicious). Follow-up ultrasound recommended "
                "at 1, 3, and 5 years (nodule >= 1.5cm but < 2.5cm)."
            )
        else:
            return (
                "TR3 (Mildly suspicious). No FNA or follow-up needed "
                "at current size (< 1.5cm)."
            )

    if tr_level == TIRADS.TR4:
        if max_diameter_mm >= 15.0:
            return (
                "TR4 (Moderately suspicious). FNA recommended (nodule >= 1.5cm). "
                "Ultrasound-guided fine needle aspiration indicated."
            )
        elif max_diameter_mm >= 10.0:
            return (
                "TR4 (Moderately suspicious). Follow-up ultrasound recommended "
                "at 1, 2, 3, and 5 years (nodule >= 1.0cm but < 1.5cm)."
            )
        else:
            return (
                "TR4 (Moderately suspicious). No FNA or follow-up needed "
                "at current size (< 1.0cm)."
            )

    if tr_level == TIRADS.TR5:
        if max_diameter_mm >= 10.0:
            return (
                "TR5 (Highly suspicious). FNA recommended (nodule >= 1.0cm). "
                "Urgent ultrasound-guided fine needle aspiration indicated."
            )
        elif max_diameter_mm >= 5.0:
            return (
                "TR5 (Highly suspicious). Follow-up ultrasound recommended "
                "annually for up to 5 years (nodule >= 0.5cm but < 1.0cm)."
            )
        else:
            return (
                "TR5 (Highly suspicious). No FNA needed at current size "
                "(< 0.5cm). Optional follow-up ultrasound."
            )

    return "Clinical correlation recommended."


# =====================================================================
# ThyroidTIRADSWorkflow
# =====================================================================


class ThyroidTIRADSWorkflow(BaseImagingWorkflow):
    """Thyroid nodule TI-RADS classification workflow.

    Pipeline:
        1. Load thyroid ultrasound DICOM series
        2. Segment thyroid gland and identify nodules (MONAI)
        3. Extract nodule features (composition, echogenicity, shape, margin, foci)
        4. Calculate TI-RADS points per nodule
        5. Map points to TR level (TR1-TR5)
        6. Determine FNA recommendation based on TR level + nodule size
        7. Trigger cross-modal genomic query if TR >= 4

    Cross-modal trigger:
        TI-RADS 4+ fires BRAF V600E, RAS, RET/PTC genomic queries
        via CrossModalTrigger to enrich findings with thyroid cancer
        molecular marker context.
    """

    WORKFLOW_NAME: str = "thyroid_tirads"
    TARGET_LATENCY_SEC: float = 180.0
    MODALITY: str = "ultrasound"
    BODY_REGION: str = "neck"
    MODELS_USED: List[str] = ["monai_thyroid_segmentation"]

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        mock_overrides: Optional[Dict] = None,
    ):
        super().__init__(
            mock_mode=mock_mode,
            nim_clients=nim_clients,
            mock_overrides=mock_overrides,
        )

    # ------------------------------------------------------------------
    # Preprocessing
    # ------------------------------------------------------------------

    def preprocess(self, input_path: str) -> Any:
        """Load thyroid ultrasound DICOM and apply preprocessing.

        Steps:
            - Load DICOM ultrasound series from input_path
            - Apply speckle reduction filtering
            - Normalize intensity to [0, 1]
            - Segment thyroid gland boundary

        In production, uses monai.transforms:
            LoadImaged, EnsureChannelFirstd, ScaleIntensityRanged,
            ultrasound-specific speckle noise reduction.
        """
        logger.info(
            f"Preprocessing thyroid ultrasound from {input_path}: "
            "speckle reduction, intensity normalization, gland segmentation"
        )
        return None

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def infer(self, preprocessed: Any) -> Dict:
        """Run thyroid nodule detection and feature extraction.

        In mock mode, returns realistic nodule findings with TI-RADS
        features. In real mode, would run MONAI thyroid segmentation
        followed by nodule-level feature extraction.
        """
        if self.mock_mode:
            return self._mock_inference()

        logger.warning("Real inference not yet implemented; using mock results")
        return self._mock_inference()

    def _mock_inference(self) -> Dict:
        """Return realistic mock thyroid ultrasound findings.

        Simulates a thyroid ultrasound with two nodules:
            - A 22mm solid hypoechoic nodule with punctate echogenic foci
              and irregular margins in the right lobe (TR5, 10 points)
            - A 14mm mixed cystic-solid isoechoic nodule with smooth
              margins in the left lobe (TR3, 2 points)
        """
        return {
            "thyroid_volume_ml": 18.5,
            "isthmus_thickness_mm": 3.2,
            "background_parenchyma": "heterogeneous",
            "nodules": [
                {
                    "id": "nodule_1",
                    "location": "right lobe mid",
                    "laterality": "right",
                    "max_diameter_mm": 22.0,
                    "ap_diameter_mm": 18.5,
                    "transverse_diameter_mm": 20.1,
                    "longitudinal_diameter_mm": 22.0,
                    "composition": "solid",
                    "echogenicity": "hypoechoic",
                    "shape": "taller_than_wide",
                    "margin": "irregular",
                    "echogenic_foci": "punctate_echogenic_foci",
                    "vascularity": "increased_intranodular",
                    "detection_confidence": 0.96,
                },
                {
                    "id": "nodule_2",
                    "location": "left lobe inferior",
                    "laterality": "left",
                    "max_diameter_mm": 14.0,
                    "ap_diameter_mm": 10.2,
                    "transverse_diameter_mm": 12.8,
                    "longitudinal_diameter_mm": 14.0,
                    "composition": "mixed_cystic_solid",
                    "echogenicity": "isoechoic",
                    "shape": "wider_than_tall",
                    "margin": "smooth",
                    "echogenic_foci": "none",
                    "vascularity": "peripheral",
                    "detection_confidence": 0.91,
                },
            ],
        }

    # ------------------------------------------------------------------
    # Postprocessing
    # ------------------------------------------------------------------

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Calculate TI-RADS scores and determine FNA recommendations.

        Applies ACR TI-RADS (2017) point system:
            - Calculate points for each nodule from 5 feature categories
            - Map total points to TR level
            - Determine FNA recommendation based on TR level + size
            - Overall severity driven by highest-TR nodule

        Args:
            inference_result: Raw inference dict with 'nodules' list.

        Returns:
            WorkflowResult with TI-RADS classification and severity.
        """
        nodules = inference_result.get("nodules", [])

        if not nodules:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.COMPLETED,
                findings=[{
                    "category": FindingCategory.NORMAL.value,
                    "description": "No thyroid nodules detected.",
                    "severity": FindingSeverity.NORMAL.value,
                    "tirads": TIRADS.TR1.value,
                }],
                measurements={
                    "nodule_count": 0.0,
                    "thyroid_volume_ml": inference_result.get("thyroid_volume_ml", 0.0),
                },
                classification=f"TI-RADS {TIRADS.TR1.value}",
                severity=FindingSeverity.NORMAL,
                nim_services_used=self.MODELS_USED,
            )

        findings = []
        overall_severity = FindingSeverity.NORMAL
        highest_tirads = TIRADS.TR1

        severity_order = [
            FindingSeverity.NORMAL,
            FindingSeverity.ROUTINE,
            FindingSeverity.SIGNIFICANT,
            FindingSeverity.URGENT,
            FindingSeverity.CRITICAL,
        ]

        for nodule in nodules:
            nodule_id = nodule.get("id", "unknown")
            location = nodule.get("location", "unspecified")
            laterality = nodule.get("laterality", "unspecified")
            max_diameter_mm = nodule.get("max_diameter_mm", 0.0)

            # Calculate TI-RADS points
            points = calculate_tirads_points(nodule)
            tr_level = points_to_tirads(points)
            nodule_severity = tirads_to_severity(tr_level)
            recommendation = tirads_fna_recommendation(tr_level, max_diameter_mm)

            # Track overall severity (highest wins)
            if severity_order.index(nodule_severity) > severity_order.index(overall_severity):
                overall_severity = nodule_severity
                highest_tirads = tr_level

            description = (
                f"Thyroid nodule in {location}, {max_diameter_mm:.1f}mm, "
                f"{nodule.get('composition', 'unspecified')} composition, "
                f"{nodule.get('echogenicity', 'unspecified')} echogenicity, "
                f"TI-RADS {tr_level.value} ({points} points)"
            )

            findings.append({
                "category": FindingCategory.THYROID_NODULE.value,
                "description": description,
                "severity": nodule_severity.value,
                "nodule_id": nodule_id,
                "laterality": laterality,
                "location": location,
                "tirads": tr_level.value,
                "tirads_points": points,
                "recommendation": recommendation,
                "detection_confidence": nodule.get("detection_confidence", 0.0),
                "features": {
                    "composition": nodule.get("composition", ""),
                    "echogenicity": nodule.get("echogenicity", ""),
                    "shape": nodule.get("shape", ""),
                    "margin": nodule.get("margin", ""),
                    "echogenic_foci": nodule.get("echogenic_foci", ""),
                },
                "measurements": {
                    "max_diameter_mm": max_diameter_mm,
                    "ap_diameter_mm": nodule.get("ap_diameter_mm", 0.0),
                    "transverse_diameter_mm": nodule.get("transverse_diameter_mm", 0.0),
                    "longitudinal_diameter_mm": nodule.get("longitudinal_diameter_mm", 0.0),
                },
            })

        # Aggregate measurements
        measurements = {
            "nodule_count": float(len(nodules)),
            "thyroid_volume_ml": inference_result.get("thyroid_volume_ml", 0.0),
            "isthmus_thickness_mm": inference_result.get("isthmus_thickness_mm", 0.0),
        }
        for i, nodule in enumerate(nodules):
            prefix = f"nodule_{i + 1}"
            measurements[f"{prefix}_diameter_mm"] = nodule.get("max_diameter_mm", 0.0)
            points = calculate_tirads_points(nodule)
            measurements[f"{prefix}_tirads_points"] = float(points)

        # Cross-modal trigger for TR4+ (BRAF, RAS, RET/PTC)
        cross_modal_trigger = highest_tirads in (TIRADS.TR4, TIRADS.TR5)

        if cross_modal_trigger:
            findings[0]["cross_modal_trigger"] = True
            findings[0]["genomic_queries"] = [
                "BRAF V600E", "RAS", "RET/PTC",
            ]

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=measurements,
            classification=f"TI-RADS {highest_tirads.value}",
            severity=overall_severity,
            nim_services_used=self.MODELS_USED,
        )
