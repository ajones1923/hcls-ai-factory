"""CT Coronary Angiography (CTA) Workflow.

Reference workflow #5: Coronary artery disease assessment.
Uses SegResNet (MONAI) for coronary segmentation and VISTA-3D for calcium scoring.
Applies CAD-RADS 2.0 classification for standardized reporting.

Pipeline:
    1. Load DICOM CT cardiac series (ECG-gated CTA)
    2. Resample to 0.5mm isotropic, apply cardiac window
    3. Run SegResNet coronary artery segmentation
    4. Run VISTA-3D calcium scoring (Agatston method)
    5. Measure per-vessel stenosis percentage, plaque characterization
    6. Apply CAD-RADS 2.0 classification
"""

from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    CADRADS,
    BodyRegion,
    FindingCategory,
    FindingSeverity,
    ImagingModality,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows.base import BaseImagingWorkflow


# CAD-RADS 2.0 stenosis thresholds (percent)
CADRADS_THRESHOLDS = {
    "cad_0_max": 0,     # 0%: no stenosis
    "cad_1_max": 24,    # 1-24%: minimal
    "cad_2_max": 49,    # 25-49%: mild
    "cad_3_max": 69,    # 50-69%: moderate
    "cad_4a_max": 99,   # 70-99%: severe
    # 100%: total occlusion -> CAD-RADS 5
}


def classify_stenosis_cadrads(
    max_stenosis_pct: float,
    left_main_stenosis_pct: float = 0.0,
    num_vessels_gte50: int = 0,
) -> CADRADS:
    """Classify coronary stenosis using CAD-RADS 2.0 criteria.

    Args:
        max_stenosis_pct: Maximum stenosis percentage across all vessels.
        left_main_stenosis_pct: Left main coronary artery stenosis percentage.
        num_vessels_gte50: Number of vessels with >= 50% stenosis.

    Returns:
        CADRADS category enum value.
    """
    # CAD-RADS 4B: left main > 50% or 3-vessel disease (>= 50% each)
    if left_main_stenosis_pct > 50 or num_vessels_gte50 >= 3:
        return CADRADS.CAD_4B

    if max_stenosis_pct >= 100:
        return CADRADS.CAD_5
    elif max_stenosis_pct >= 70:
        return CADRADS.CAD_4A
    elif max_stenosis_pct >= 50:
        return CADRADS.CAD_3
    elif max_stenosis_pct >= 25:
        return CADRADS.CAD_2
    elif max_stenosis_pct >= 1:
        return CADRADS.CAD_1
    else:
        return CADRADS.CAD_0


def cadrads_to_severity(category: CADRADS) -> FindingSeverity:
    """Map CAD-RADS 2.0 category to clinical severity."""
    mapping = {
        CADRADS.CAD_0: FindingSeverity.ROUTINE,
        CADRADS.CAD_1: FindingSeverity.ROUTINE,
        CADRADS.CAD_2: FindingSeverity.ROUTINE,
        CADRADS.CAD_3: FindingSeverity.SIGNIFICANT,
        CADRADS.CAD_4A: FindingSeverity.URGENT,
        CADRADS.CAD_4B: FindingSeverity.URGENT,
        CADRADS.CAD_5: FindingSeverity.CRITICAL,
        CADRADS.CAD_N: FindingSeverity.ROUTINE,
    }
    return mapping.get(category, FindingSeverity.ROUTINE)


def cadrads_recommendation(category: CADRADS) -> str:
    """Return CAD-RADS 2.0 management recommendation."""
    recommendations = {
        CADRADS.CAD_0: (
            "No stenosis. No further cardiac workup needed. "
            "Consider risk factor modification."
        ),
        CADRADS.CAD_1: (
            "Minimal stenosis (1-24%). No further cardiac workup needed. "
            "Preventive pharmacotherapy and lifestyle modification."
        ),
        CADRADS.CAD_2: (
            "Mild stenosis (25-49%). Consider preventive pharmacotherapy. "
            "Risk factor modification recommended."
        ),
        CADRADS.CAD_3: (
            "Moderate stenosis (50-69%). Consider functional assessment "
            "(stress test, FFR-CT). Guideline-directed medical therapy."
        ),
        CADRADS.CAD_4A: (
            "Severe stenosis (70-99%). Consider invasive coronary angiography "
            "or FFR-CT. Cardiology referral recommended."
        ),
        CADRADS.CAD_4B: (
            "Left main >50% or 3-vessel disease. Urgent cardiology referral. "
            "Consider invasive coronary angiography for revascularization planning."
        ),
        CADRADS.CAD_5: (
            "Total occlusion (100%). Urgent cardiology referral. "
            "Evaluate for chronic total occlusion intervention or bypass surgery."
        ),
        CADRADS.CAD_N: (
            "Non-diagnostic study. Consider repeat CTA with optimized protocol "
            "or alternative functional testing."
        ),
    }
    return recommendations.get(category, "Clinical correlation recommended.")


class CTCoronaryAngiographyWorkflow(BaseImagingWorkflow):
    """Coronary CTA analysis, stenosis grading, and CAD-RADS 2.0 classification.

    Pipeline:
        1. Load DICOM ECG-gated cardiac CT series
        2. Resample to 0.5mm isotropic, apply cardiac window
        3. Run SegResNet coronary artery segmentation (MONAI)
        4. Run VISTA-3D calcium scoring (Agatston method)
        5. Per-vessel stenosis grading, plaque characterization
        6. Apply CAD-RADS 2.0 classification
    """

    WORKFLOW_NAME = "ct_coronary_angiography"
    MODALITY = ImagingModality.CT
    BODY_REGION = BodyRegion.CARDIAC
    TARGET_LATENCY_SEC = 180.0
    MODELS_USED = ["SegResNet (MONAI coronary segmentation)", "VISTA-3D (calcium scoring)"]

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        mock_overrides: Optional[Dict] = None,
    ):
        super().__init__(mock_mode=mock_mode, nim_clients=nim_clients, mock_overrides=mock_overrides)

    # ------------------------------------------------------------------
    # Preprocessing
    # ------------------------------------------------------------------

    def preprocess(self, input_path: str) -> Any:
        """Load DICOM ECG-gated cardiac CT and apply cardiac-optimized preprocessing.

        Steps:
            - Load DICOM series from input_path
            - Reorient to RAS coordinate system
            - Resample to 0.5mm isotropic voxel spacing
            - Apply cardiac window: center=200 HU, width=800 HU
            - Normalize intensity to [0, 1]
        """
        logger.info(
            f"Preprocessing cardiac CTA from {input_path}: "
            "reorient RAS, resample 0.5mm iso, cardiac window"
        )
        return None

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def infer(self, preprocessed: Any) -> Dict:
        """Run coronary segmentation + calcium scoring.

        In mock mode, returns realistic cardiac CTA findings.
        Falls back to mock results if no model is loaded.
        """
        if self.mock_mode:
            return self._mock_inference()

        logger.warning("Real inference not yet implemented — falling back to mock")
        return self._mock_inference()

    def _mock_inference(self) -> Dict:
        """Return realistic mock cardiac CTA result.

        Simulates a patient with significant LAD stenosis, moderate calcium
        burden, and high-risk plaque features.
        """
        return {
            "calcium_score": 385,
            "calcium_score_method": "agatston",
            "calcium_percentile": 82,
            "vessels": [
                {
                    "name": "LAD",
                    "full_name": "Left Anterior Descending",
                    "max_stenosis_pct": 72,
                    "plaque_type": "mixed",
                    "segment": "mid",
                },
                {
                    "name": "LCx",
                    "full_name": "Left Circumflex",
                    "max_stenosis_pct": 30,
                    "plaque_type": "calcified",
                    "segment": "proximal",
                },
                {
                    "name": "RCA",
                    "full_name": "Right Coronary Artery",
                    "max_stenosis_pct": 15,
                    "plaque_type": "none",
                    "segment": "proximal",
                },
                {
                    "name": "Left Main",
                    "full_name": "Left Main Coronary Artery",
                    "max_stenosis_pct": 0,
                    "plaque_type": "none",
                    "segment": "ostial",
                },
            ],
            "high_risk_plaque_features": [
                "low-attenuation plaque",
                "positive remodeling",
                "napkin-ring sign",
            ],
            "ejection_fraction_estimate": 55,
        }

    # ------------------------------------------------------------------
    # Postprocessing
    # ------------------------------------------------------------------

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Apply CAD-RADS 2.0 classification and determine severity.

        CAD-RADS 2.0 stenosis thresholds:
            0%:      CAD-RADS 0 (no stenosis)
            1-24%:   CAD-RADS 1 (minimal)
            25-49%:  CAD-RADS 2 (mild)
            50-69%:  CAD-RADS 3 (moderate)
            70-99%:  CAD-RADS 4A (severe)
            LM>50% or 3-vessel: CAD-RADS 4B
            100%:    CAD-RADS 5 (total occlusion)
        """
        vessels = inference_result.get("vessels", [])
        calcium_score = inference_result.get("calcium_score", 0)
        high_risk_features = inference_result.get("high_risk_plaque_features", [])
        ef_estimate = inference_result.get("ejection_fraction_estimate", 0)

        # Compute per-vessel stenosis stats
        max_stenosis = 0.0
        left_main_stenosis = 0.0
        num_vessels_gte50 = 0

        for vessel in vessels:
            stenosis = vessel.get("max_stenosis_pct", 0)
            if stenosis > max_stenosis:
                max_stenosis = stenosis
            if vessel.get("name") == "Left Main":
                left_main_stenosis = stenosis
            if stenosis >= 50:
                num_vessels_gte50 += 1

        # Apply CAD-RADS 2.0 classification
        cadrads = classify_stenosis_cadrads(
            max_stenosis_pct=max_stenosis,
            left_main_stenosis_pct=left_main_stenosis,
            num_vessels_gte50=num_vessels_gte50,
        )
        severity = cadrads_to_severity(cadrads)
        recommendation = cadrads_recommendation(cadrads)

        # Build findings list
        findings: List[Dict] = []

        # Primary finding: overall CAD-RADS classification
        findings.append({
            "category": FindingCategory.STENOSIS.value,
            "description": (
                f"Coronary CTA: CAD-RADS {cadrads.value}. "
                f"Maximum stenosis {max_stenosis:.0f}% "
                f"(Agatston calcium score: {calcium_score})."
            ),
            "severity": severity.value,
            "cadrads": cadrads.value,
            "recommendation": recommendation,
        })

        # Per-vessel findings
        for vessel in vessels:
            stenosis = vessel.get("max_stenosis_pct", 0)
            if stenosis > 0:
                findings.append({
                    "category": FindingCategory.STENOSIS.value,
                    "description": (
                        f"{vessel['full_name']} ({vessel['name']}): "
                        f"{stenosis}% stenosis ({vessel.get('segment', '')} segment), "
                        f"plaque type: {vessel.get('plaque_type', 'unknown')}"
                    ),
                    "severity": cadrads_to_severity(
                        classify_stenosis_cadrads(stenosis)
                    ).value,
                    "vessel": vessel["name"],
                    "stenosis_pct": stenosis,
                    "plaque_type": vessel.get("plaque_type", "unknown"),
                })

        # High-risk plaque features finding
        if high_risk_features:
            findings.append({
                "category": FindingCategory.CALCIFICATION.value,
                "description": (
                    f"High-risk plaque features identified: "
                    f"{', '.join(high_risk_features)}"
                ),
                "severity": FindingSeverity.SIGNIFICANT.value,
                "features": high_risk_features,
            })

        # Build measurements dict
        measurements: Dict[str, float] = {
            "calcium_score_agatston": float(calcium_score),
            "max_stenosis_pct": max_stenosis,
            "num_vessels_with_stenosis": float(
                sum(1 for v in vessels if v.get("max_stenosis_pct", 0) > 0)
            ),
            "num_vessels_gte50_pct": float(num_vessels_gte50),
            "ejection_fraction_pct": float(ef_estimate),
        }

        # Per-vessel stenosis measurements
        for vessel in vessels:
            key = vessel["name"].lower().replace(" ", "_")
            measurements[f"{key}_stenosis_pct"] = float(
                vessel.get("max_stenosis_pct", 0)
            )

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=measurements,
            classification=f"CAD-RADS {cadrads.value}",
            severity=severity,
            nim_services_used=self.MODELS_USED,
        )
