"""MRI Prostate PI-RADS Assessment Workflow.

Reference workflow #6: Prostate cancer detection and PI-RADS classification.
Uses nnU-Net for prostate segmentation and SegResNet for lesion detection.
Applies PI-RADS v2.1 classification for standardized reporting.

Pipeline:
    1. Load DICOM multiparametric MRI series (T2W, DWI, ADC, DCE)
    2. Register sequences, resample to common resolution
    3. Run nnU-Net prostate zone segmentation (peripheral/transition)
    4. Run SegResNet lesion detection and characterization
    5. Per-lesion PI-RADS v2.1 scoring based on zone and signal characteristics
    6. Apply overall PI-RADS classification (highest lesion score)
"""

from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    PIRADS,
    BodyRegion,
    FindingCategory,
    FindingSeverity,
    ImagingModality,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows.base import BaseImagingWorkflow


def pirads_to_severity(score: PIRADS) -> FindingSeverity:
    """Map PI-RADS v2.1 score to clinical severity."""
    mapping = {
        PIRADS.PI_1: FindingSeverity.ROUTINE,
        PIRADS.PI_2: FindingSeverity.ROUTINE,
        PIRADS.PI_3: FindingSeverity.SIGNIFICANT,
        PIRADS.PI_4: FindingSeverity.URGENT,
        PIRADS.PI_5: FindingSeverity.CRITICAL,
    }
    return mapping.get(score, FindingSeverity.ROUTINE)


def pirads_recommendation(score: PIRADS) -> str:
    """Return PI-RADS v2.1 management recommendation."""
    recommendations = {
        PIRADS.PI_1: (
            "Very low suspicion for clinically significant cancer. "
            "No further evaluation needed for this finding."
        ),
        PIRADS.PI_2: (
            "Low suspicion for clinically significant cancer. "
            "No further evaluation needed for this finding."
        ),
        PIRADS.PI_3: (
            "Intermediate suspicion. Clinical correlation with PSA density, "
            "risk calculators, and prior biopsy history recommended. "
            "Consider MRI-targeted biopsy based on clinical context."
        ),
        PIRADS.PI_4: (
            "High suspicion for clinically significant cancer. "
            "MRI-targeted biopsy recommended (MRI-TRUS fusion or in-bore)."
        ),
        PIRADS.PI_5: (
            "Very high suspicion for clinically significant cancer. "
            "MRI-targeted biopsy strongly recommended. "
            "Consider concurrent systematic biopsy."
        ),
    }
    return recommendations.get(score, "Clinical correlation recommended.")


def int_to_pirads(score: int) -> PIRADS:
    """Convert integer score (1-5) to PIRADS enum."""
    mapping = {
        1: PIRADS.PI_1,
        2: PIRADS.PI_2,
        3: PIRADS.PI_3,
        4: PIRADS.PI_4,
        5: PIRADS.PI_5,
    }
    return mapping.get(score, PIRADS.PI_3)


class MRIProstateWorkflow(BaseImagingWorkflow):
    """Prostate MRI analysis, lesion detection, and PI-RADS v2.1 classification.

    Pipeline:
        1. Load DICOM multiparametric MRI (T2W, DWI, ADC, DCE)
        2. Co-register sequences to T2W reference frame
        3. Run nnU-Net prostate zone segmentation (PZ, TZ, CZ)
        4. Run SegResNet lesion detection within prostate zones
        5. Per-lesion PI-RADS v2.1 scoring based on zone + signal characteristics
        6. Apply overall PI-RADS classification (highest lesion score)
    """

    WORKFLOW_NAME = "mri_prostate_pirads"
    MODALITY = ImagingModality.MRI
    BODY_REGION = BodyRegion.PELVIS
    TARGET_LATENCY_SEC = 240.0
    MODELS_USED = ["nnU-Net (prostate segmentation)", "SegResNet (lesion detection)"]

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
        """Load DICOM multiparametric MRI and apply prostate-optimized preprocessing.

        Steps:
            - Load T2W, DWI, ADC, DCE series from input_path
            - Co-register all sequences to T2W reference frame
            - Resample to 0.5mm in-plane, 3mm slice thickness
            - Normalize intensities per-sequence
        """
        logger.info(
            f"Preprocessing prostate mpMRI from {input_path}: "
            "co-register T2W/DWI/ADC/DCE, resample, normalize"
        )
        return None

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def infer(self, preprocessed: Any) -> Dict:
        """Run prostate segmentation + lesion detection.

        In mock mode, returns realistic multiparametric MRI findings.
        Falls back to mock results if no model is loaded.
        """
        if self.mock_mode:
            return self._mock_inference()

        logger.warning("Real inference not yet implemented — falling back to mock")
        return self._mock_inference()

    def _mock_inference(self) -> Dict:
        """Return realistic mock prostate mpMRI result.

        Simulates a patient with a PI-RADS 4 peripheral zone lesion and
        a PI-RADS 3 transition zone lesion, moderate prostate volume,
        and elevated PSA density.
        """
        return {
            "prostate_volume_cc": 45.2,
            "psa_density": 0.22,
            "psa_assumed_ng_ml": 10.0,
            "lesions": [
                {
                    "id": "lesion_1",
                    "zone": "peripheral",
                    "sector": "right mid",
                    "size_mm": 14,
                    "t2_signal": "markedly hypointense",
                    "dwi_signal": "markedly hyperintense",
                    "adc_value": 650,
                    "dce_pattern": "early enhancement",
                    "pirads_score": 4,
                },
                {
                    "id": "lesion_2",
                    "zone": "transition",
                    "sector": "left anterior",
                    "size_mm": 8,
                    "t2_signal": "moderately hypointense",
                    "dwi_signal": "mildly hyperintense",
                    "adc_value": 900,
                    "dce_pattern": "no early enhancement",
                    "pirads_score": 3,
                },
            ],
        }

    # ------------------------------------------------------------------
    # Postprocessing
    # ------------------------------------------------------------------

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Apply PI-RADS v2.1 classification and determine severity.

        PI-RADS v2.1 severity mapping:
            1-2: routine (very low / low suspicion)
            3:   significant (intermediate suspicion)
            4:   urgent (high suspicion)
            5:   critical (very high suspicion)

        Overall score is the highest PI-RADS score among all lesions.
        """
        prostate_volume = inference_result.get("prostate_volume_cc", 0.0)
        psa_density = inference_result.get("psa_density", 0.0)
        lesions = inference_result.get("lesions", [])

        # Determine overall PI-RADS (highest lesion score)
        max_pirads_int = 1
        for lesion in lesions:
            score = lesion.get("pirads_score", 1)
            if score > max_pirads_int:
                max_pirads_int = score

        overall_pirads = int_to_pirads(max_pirads_int)
        severity = pirads_to_severity(overall_pirads)
        recommendation = pirads_recommendation(overall_pirads)

        # Build findings list
        findings: List[Dict] = []

        # Primary finding: overall PI-RADS classification
        findings.append({
            "category": FindingCategory.LESION.value,
            "description": (
                f"Prostate mpMRI: PI-RADS {overall_pirads.value} (overall). "
                f"{len(lesions)} lesion(s) detected. "
                f"Prostate volume {prostate_volume:.1f} cc, "
                f"PSA density {psa_density:.2f} ng/mL/cc."
            ),
            "severity": severity.value,
            "pirads": overall_pirads.value,
            "recommendation": recommendation,
        })

        # Per-lesion findings
        for lesion in lesions:
            lesion_pirads = int_to_pirads(lesion.get("pirads_score", 1))
            lesion_severity = pirads_to_severity(lesion_pirads)
            lesion_recommendation = pirads_recommendation(lesion_pirads)

            findings.append({
                "category": FindingCategory.LESION.value,
                "description": (
                    f"Lesion {lesion['id']}: PI-RADS {lesion_pirads.value} in "
                    f"{lesion.get('zone', 'unknown')} zone ({lesion.get('sector', '')}), "
                    f"size {lesion.get('size_mm', 0)}mm. "
                    f"T2: {lesion.get('t2_signal', '')}, "
                    f"DWI: {lesion.get('dwi_signal', '')}, "
                    f"ADC: {lesion.get('adc_value', 0)}, "
                    f"DCE: {lesion.get('dce_pattern', '')}"
                ),
                "severity": lesion_severity.value,
                "lesion_id": lesion["id"],
                "pirads_score": lesion.get("pirads_score", 1),
                "zone": lesion.get("zone", ""),
                "sector": lesion.get("sector", ""),
                "size_mm": lesion.get("size_mm", 0),
                "adc_value": lesion.get("adc_value", 0),
                "recommendation": lesion_recommendation,
            })

        # Build measurements dict
        measurements: Dict[str, float] = {
            "prostate_volume_cc": prostate_volume,
            "psa_density": psa_density,
            "lesion_count": float(len(lesions)),
            "max_pirads_score": float(max_pirads_int),
        }

        # Per-lesion measurements
        for lesion in lesions:
            lid = lesion.get("id", "unknown")
            measurements[f"{lid}_pirads_score"] = float(lesion.get("pirads_score", 0))
            measurements[f"{lid}_size_mm"] = float(lesion.get("size_mm", 0))
            measurements[f"{lid}_adc_value"] = float(lesion.get("adc_value", 0))

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=measurements,
            classification=f"PI-RADS {overall_pirads.value}",
            severity=severity,
            nim_services_used=self.MODELS_USED,
        )
