"""Liver LI-RADS Classification Workflow.

Clinical workflow for liver CT/MRI with ACR LI-RADS scoring.
Cross-modal trigger: LI-RADS 4+ -> TP53, CTNNB1, TERT genomic queries.

Scoring system: ACR LI-RADS v2018
  For observations in patients at risk for HCC:
  LR-1: Definitely benign
  LR-2: Probably benign
  LR-3: Intermediate probability of malignancy
  LR-4: Probably HCC
  LR-5: Definitely HCC
  LR-M: Probably or definitely malignant, not HCC specific
  LR-TIV: Tumor in vein

  Major features:
  - Arterial phase hyperenhancement (APHE)
  - Observation size (mm)
  - Enhancing capsule
  - Nonperipheral washout
  - Threshold growth (>= 50% increase in <= 6 months)

Target latency: < 5 minutes
Models: MONAI liver segmentation

Author: Adam Jones
Date: April 2026
"""

from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    FindingCategory,
    FindingSeverity,
    LIRADS,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows.base import BaseImagingWorkflow


# =====================================================================
# LI-RADS helper functions
# =====================================================================


def classify_lirads(observation: Dict) -> LIRADS:
    """Apply ACR LI-RADS v2018 diagnostic table to classify an observation.

    The LI-RADS diagnostic table uses arterial phase hyperenhancement
    (APHE), observation size, and additional major features (washout,
    enhancing capsule) to determine the LR category.

    Special categories:
        LR-TIV: Assigned when tumor in vein is detected, regardless of
                 other features.
        LR-M:   Assigned when features suggest malignancy but are not
                 HCC-specific (e.g., targetoid appearance, peripheral
                 washout, non-enhancement).

    Args:
        observation: Dict with keys: aphe, size_mm, washout, capsule,
                     threshold_growth, tumor_in_vein, targetoid_appearance.

    Returns:
        LIRADS enum value.
    """
    # Check for special categories first
    if observation.get("tumor_in_vein", False):
        return LIRADS.LR_TIV

    if observation.get("targetoid_appearance", False):
        return LIRADS.LR_M

    aphe = observation.get("aphe", False)
    size_mm = observation.get("size_mm", 0.0)
    washout = observation.get("washout", False)
    capsule = observation.get("capsule", False)
    threshold_growth = observation.get("threshold_growth", False)

    # Count additional major features (washout, capsule, threshold_growth)
    additional_features = sum([washout, capsule, threshold_growth])

    # ACR LI-RADS v2018 Diagnostic Table
    # -----------------------------------------------------------------
    # APHE present:
    #   < 10mm:
    #     0 additional features -> LR-3
    #     >= 1 additional feature -> LR-4
    #   10-19mm:
    #     0 additional features -> LR-3
    #     1 additional feature -> LR-4
    #     >= 2 additional features -> LR-5
    #   >= 20mm:
    #     0 additional features -> LR-4
    #     >= 1 additional feature -> LR-5
    #
    # No APHE:
    #   < 20mm:
    #     0 additional features -> LR-3
    #     >= 1 additional feature -> LR-3 (still intermediate)
    #   >= 20mm:
    #     0 additional features -> LR-3
    #     >= 1 additional feature -> LR-4
    # -----------------------------------------------------------------

    if aphe:
        if size_mm < 10.0:
            if additional_features == 0:
                return LIRADS.LR_3
            else:
                return LIRADS.LR_4
        elif size_mm < 20.0:
            if additional_features == 0:
                return LIRADS.LR_3
            elif additional_features == 1:
                return LIRADS.LR_4
            else:
                return LIRADS.LR_5
        else:  # >= 20mm
            if additional_features == 0:
                return LIRADS.LR_4
            else:
                return LIRADS.LR_5
    else:
        # No APHE
        if size_mm < 20.0:
            return LIRADS.LR_3
        else:
            if additional_features >= 1:
                return LIRADS.LR_4
            else:
                return LIRADS.LR_3


def lirads_to_severity(lr_category: LIRADS) -> FindingSeverity:
    """Map LI-RADS category to clinical severity.

    Args:
        lr_category: ACR LI-RADS category.

    Returns:
        FindingSeverity corresponding to the LR category.
    """
    mapping = {
        LIRADS.LR_1: FindingSeverity.NORMAL,
        LIRADS.LR_2: FindingSeverity.ROUTINE,
        LIRADS.LR_3: FindingSeverity.SIGNIFICANT,
        LIRADS.LR_4: FindingSeverity.URGENT,
        LIRADS.LR_5: FindingSeverity.CRITICAL,
        LIRADS.LR_M: FindingSeverity.CRITICAL,
        LIRADS.LR_TIV: FindingSeverity.CRITICAL,
    }
    return mapping.get(lr_category, FindingSeverity.ROUTINE)


def lirads_recommendation(lr_category: LIRADS) -> str:
    """Return ACR LI-RADS v2018 management recommendation.

    Args:
        lr_category: ACR LI-RADS category.

    Returns:
        Management recommendation string.
    """
    recommendations = {
        LIRADS.LR_1: (
            "LR-1 (Definitely benign). No additional workup needed. "
            "Return to routine surveillance."
        ),
        LIRADS.LR_2: (
            "LR-2 (Probably benign). Continue routine surveillance per "
            "clinical guidelines. No additional workup required."
        ),
        LIRADS.LR_3: (
            "LR-3 (Intermediate probability). Repeat diagnostic imaging "
            "in 3-6 months. Consider alternate contrast agent or modality. "
            "Multidisciplinary discussion if features evolve."
        ),
        LIRADS.LR_4: (
            "LR-4 (Probably HCC). Multidisciplinary discussion recommended. "
            "Consider biopsy or additional imaging with alternate contrast "
            "agent for further characterization."
        ),
        LIRADS.LR_5: (
            "LR-5 (Definitely HCC). HCC can be diagnosed without biopsy. "
            "Multidisciplinary discussion for treatment planning. "
            "Stage per BCLC or institution guidelines."
        ),
        LIRADS.LR_M: (
            "LR-M (Probably or definitely malignant, not HCC specific). "
            "Biopsy recommended to determine etiology. "
            "Multidisciplinary discussion for management."
        ),
        LIRADS.LR_TIV: (
            "LR-TIV (Tumor in vein). Definite tumor thrombosis detected. "
            "STAT multidisciplinary discussion. Stage as advanced HCC "
            "(BCLC C). Consider systemic therapy eligibility."
        ),
    }
    return recommendations.get(lr_category, "Clinical correlation recommended.")


# =====================================================================
# LiverLIRADSWorkflow
# =====================================================================


class LiverLIRADSWorkflow(BaseImagingWorkflow):
    """Liver LI-RADS classification workflow.

    Pipeline:
        1. Load contrast-enhanced CT or MRI DICOM series
           (arterial, portal venous, delayed phases)
        2. Segment liver parenchyma and identify observations (MONAI)
        3. Evaluate major features: APHE, washout, capsule, size, growth
        4. Apply ACR LI-RADS v2018 diagnostic table
        5. Check for special categories (LR-TIV, LR-M)
        6. Trigger cross-modal genomic query if LR >= 4

    Cross-modal trigger:
        LI-RADS 4+ fires TP53, CTNNB1, TERT promoter genomic queries
        via CrossModalTrigger to enrich findings with HCC molecular
        marker context.
    """

    WORKFLOW_NAME: str = "liver_lirads"
    TARGET_LATENCY_SEC: float = 300.0
    MODALITY: str = "ct"
    BODY_REGION: str = "abdomen"
    MODELS_USED: List[str] = ["monai_liver_segmentation"]

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
        """Load liver CT/MRI DICOM and apply preprocessing.

        Steps:
            - Load multi-phase DICOM series (arterial, portal venous, delayed)
            - Register phases to common coordinate space
            - Apply liver window: center=80 HU, width=200 HU (CT)
            - Segment liver parenchyma for region-of-interest analysis

        In production, uses monai.transforms:
            LoadImaged, EnsureChannelFirstd, Orientationd, Spacingd,
            ScaleIntensityRanged, liver-specific windowing.
        """
        logger.info(
            f"Preprocessing liver CT/MRI from {input_path}: "
            "multi-phase registration, liver window, parenchyma segmentation"
        )
        return None

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def infer(self, preprocessed: Any) -> Dict:
        """Run liver observation detection and characterization.

        In mock mode, returns realistic liver findings with LI-RADS
        major features. In real mode, would run MONAI liver segmentation
        followed by multi-phase enhancement analysis.
        """
        if self.mock_mode:
            return self._mock_inference()

        logger.warning("Real inference not yet implemented; using mock results")
        return self._mock_inference()

    def _mock_inference(self) -> Dict:
        """Return realistic mock liver CT findings.

        Simulates a contrast-enhanced CT in a cirrhotic patient with:
            - A 25mm arterially enhancing observation in segment 6 with
              washout and capsule (LR-5: definitely HCC)
            - A 12mm observation in segment 8 with APHE but no additional
              features (LR-3: intermediate)
        """
        return {
            "background_liver": "cirrhosis",
            "liver_volume_ml": 1450.0,
            "portal_vein_patent": True,
            "hepatic_veins_patent": True,
            "ascites": False,
            "splenomegaly": True,
            "observations": [
                {
                    "id": "obs_1",
                    "segment": "6",
                    "location": "right hepatic lobe, segment 6",
                    "size_mm": 25.0,
                    "aphe": True,
                    "washout": True,
                    "capsule": True,
                    "threshold_growth": False,
                    "tumor_in_vein": False,
                    "targetoid_appearance": False,
                    "restricted_diffusion": False,
                    "t2_hyperintensity": True,
                    "fat_in_mass": False,
                    "blood_products": False,
                    "detection_confidence": 0.97,
                },
                {
                    "id": "obs_2",
                    "segment": "8",
                    "location": "right hepatic lobe, segment 8",
                    "size_mm": 12.0,
                    "aphe": True,
                    "washout": False,
                    "capsule": False,
                    "threshold_growth": False,
                    "tumor_in_vein": False,
                    "targetoid_appearance": False,
                    "restricted_diffusion": False,
                    "t2_hyperintensity": False,
                    "fat_in_mass": False,
                    "blood_products": False,
                    "detection_confidence": 0.82,
                },
            ],
        }

    # ------------------------------------------------------------------
    # Postprocessing
    # ------------------------------------------------------------------

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Apply LI-RADS diagnostic table and determine severity.

        Applies ACR LI-RADS v2018:
            - Classify each observation via the diagnostic table
            - Check for special categories (LR-TIV, LR-M)
            - Overall severity driven by highest-category observation
            - Cross-modal trigger fires if LR >= 4

        Args:
            inference_result: Raw inference dict with 'observations' list.

        Returns:
            WorkflowResult with LI-RADS classification and severity.
        """
        observations = inference_result.get("observations", [])
        background_liver = inference_result.get("background_liver", "unknown")

        if not observations:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.COMPLETED,
                findings=[{
                    "category": FindingCategory.NORMAL.value,
                    "description": "No suspicious liver observations detected.",
                    "severity": FindingSeverity.NORMAL.value,
                    "lirads": LIRADS.LR_1.value,
                }],
                measurements={
                    "observation_count": 0.0,
                    "liver_volume_ml": inference_result.get("liver_volume_ml", 0.0),
                },
                classification=f"LI-RADS {LIRADS.LR_1.value}",
                severity=FindingSeverity.NORMAL,
                nim_services_used=self.MODELS_USED,
            )

        findings = []
        overall_severity = FindingSeverity.NORMAL
        highest_lirads = LIRADS.LR_1

        severity_order = [
            FindingSeverity.NORMAL,
            FindingSeverity.ROUTINE,
            FindingSeverity.SIGNIFICANT,
            FindingSeverity.URGENT,
            FindingSeverity.CRITICAL,
        ]

        # LI-RADS category ordering for comparison
        lirads_severity_rank = {
            LIRADS.LR_1: 0,
            LIRADS.LR_2: 1,
            LIRADS.LR_3: 2,
            LIRADS.LR_4: 3,
            LIRADS.LR_5: 4,
            LIRADS.LR_M: 5,
            LIRADS.LR_TIV: 6,
        }

        for observation in observations:
            obs_id = observation.get("id", "unknown")
            location = observation.get("location", "unspecified")
            segment = observation.get("segment", "")
            size_mm = observation.get("size_mm", 0.0)

            # Apply LI-RADS diagnostic table
            lr_category = classify_lirads(observation)
            obs_severity = lirads_to_severity(lr_category)
            recommendation = lirads_recommendation(lr_category)

            # Track overall severity (highest wins)
            if severity_order.index(obs_severity) > severity_order.index(overall_severity):
                overall_severity = obs_severity
            if lirads_severity_rank.get(lr_category, 0) > lirads_severity_rank.get(highest_lirads, 0):
                highest_lirads = lr_category

            # Build feature summary
            features_present = []
            if observation.get("aphe"):
                features_present.append("APHE")
            if observation.get("washout"):
                features_present.append("washout")
            if observation.get("capsule"):
                features_present.append("enhancing capsule")
            if observation.get("threshold_growth"):
                features_present.append("threshold growth")
            if observation.get("tumor_in_vein"):
                features_present.append("tumor in vein")

            feature_text = ", ".join(features_present) if features_present else "no major features"

            description = (
                f"Liver observation in {location} (segment {segment}), "
                f"{size_mm:.1f}mm, {feature_text}, "
                f"LI-RADS {lr_category.value}"
            )

            findings.append({
                "category": FindingCategory.LESION.value,
                "description": description,
                "severity": obs_severity.value,
                "observation_id": obs_id,
                "segment": segment,
                "location": location,
                "lirads": lr_category.value,
                "recommendation": recommendation,
                "detection_confidence": observation.get("detection_confidence", 0.0),
                "major_features": {
                    "aphe": observation.get("aphe", False),
                    "washout": observation.get("washout", False),
                    "capsule": observation.get("capsule", False),
                    "threshold_growth": observation.get("threshold_growth", False),
                    "tumor_in_vein": observation.get("tumor_in_vein", False),
                },
                "measurements": {
                    "size_mm": size_mm,
                },
            })

        # Aggregate measurements
        measurements = {
            "observation_count": float(len(observations)),
            "liver_volume_ml": inference_result.get("liver_volume_ml", 0.0),
        }
        for i, obs in enumerate(observations):
            prefix = f"obs_{i + 1}"
            measurements[f"{prefix}_size_mm"] = obs.get("size_mm", 0.0)

        # Background liver context
        if background_liver == "cirrhosis":
            measurements["cirrhosis"] = 1.0
        if inference_result.get("ascites"):
            measurements["ascites"] = 1.0
        if inference_result.get("splenomegaly"):
            measurements["splenomegaly"] = 1.0

        # Cross-modal trigger for LR-4+ (TP53, CTNNB1, TERT)
        cross_modal_trigger = highest_lirads in (
            LIRADS.LR_4, LIRADS.LR_5, LIRADS.LR_M, LIRADS.LR_TIV,
        )

        if cross_modal_trigger:
            findings[0]["cross_modal_trigger"] = True
            findings[0]["genomic_queries"] = [
                "TP53", "CTNNB1", "TERT promoter",
            ]

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=findings,
            measurements=measurements,
            classification=f"LI-RADS {highest_lirads.value}",
            severity=overall_severity,
            nim_services_used=self.MODELS_USED,
        )
