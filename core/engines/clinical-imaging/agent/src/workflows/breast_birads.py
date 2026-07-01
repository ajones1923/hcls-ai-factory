"""Breast Imaging BI-RADS Classification Workflow.

Clinical workflow for mammography and breast MRI with BI-RADS scoring.
Cross-modal trigger: BI-RADS 4+ -> BRCA1, BRCA2, PALB2, CHEK2 genomic queries.

Scoring system: ACR BI-RADS 5th Edition (2013)
  BI-RADS 0: Incomplete -- additional imaging needed
  BI-RADS 1: Negative -- routine screening
  BI-RADS 2: Benign -- routine screening
  BI-RADS 3: Probably benign -- short-interval follow-up (6 months)
  BI-RADS 4: Suspicious -- tissue diagnosis recommended
    4A: Low suspicion (2-10% malignancy)
    4B: Moderate suspicion (10-50%)
    4C: High suspicion (50-95%)
  BI-RADS 5: Highly suggestive of malignancy (>95%)
  BI-RADS 6: Known biopsy-proven malignancy

Target latency: < 5 minutes
Models: MONAI segmentation (breast tissue, lesion detection)

Author: Adam Jones
Date: April 2026
"""

from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import (
    BiRADS,
    FindingCategory,
    FindingSeverity,
    WorkflowResult,
    WorkflowStatus,
)
from src.workflows.base import BaseImagingWorkflow


# =====================================================================
# BI-RADS helper functions
# =====================================================================


def birads_to_severity(category: BiRADS) -> FindingSeverity:
    """Map BI-RADS category to clinical severity.

    Args:
        category: ACR BI-RADS category (0-6).

    Returns:
        FindingSeverity corresponding to the BI-RADS category.

    Mapping:
        BI-RADS 0: ROUTINE (needs additional imaging)
        BI-RADS 1: NORMAL
        BI-RADS 2: NORMAL
        BI-RADS 3: SIGNIFICANT (short-interval follow-up)
        BI-RADS 4: URGENT (tissue diagnosis recommended)
        BI-RADS 5: CRITICAL (highly suggestive of malignancy)
        BI-RADS 6: CRITICAL (known malignancy)
    """
    mapping = {
        BiRADS.CAT_0: FindingSeverity.ROUTINE,
        BiRADS.CAT_1: FindingSeverity.NORMAL,
        BiRADS.CAT_2: FindingSeverity.NORMAL,
        BiRADS.CAT_3: FindingSeverity.SIGNIFICANT,
        BiRADS.CAT_4: FindingSeverity.URGENT,
        BiRADS.CAT_5: FindingSeverity.CRITICAL,
        BiRADS.CAT_6: FindingSeverity.CRITICAL,
    }
    return mapping.get(category, FindingSeverity.ROUTINE)


def birads_recommendation(category: BiRADS) -> str:
    """Return ACR BI-RADS 5th Edition management recommendation.

    Args:
        category: ACR BI-RADS category (0-6).

    Returns:
        Management recommendation string.
    """
    recommendations = {
        BiRADS.CAT_0: (
            "Additional imaging evaluation needed. "
            "Recall for diagnostic mammography, ultrasound, or MRI."
        ),
        BiRADS.CAT_1: (
            "Negative. Continue routine screening mammography."
        ),
        BiRADS.CAT_2: (
            "Benign finding. Continue routine screening mammography."
        ),
        BiRADS.CAT_3: (
            "Probably benign. Short-interval follow-up mammography "
            "at 6 months recommended. <2% likelihood of malignancy."
        ),
        BiRADS.CAT_4: (
            "Suspicious abnormality. Tissue diagnosis recommended. "
            "Biopsy should be considered."
        ),
        BiRADS.CAT_5: (
            "Highly suggestive of malignancy. Tissue diagnosis recommended. "
            "Appropriate action should be taken. >95% likelihood of malignancy."
        ),
        BiRADS.CAT_6: (
            "Known biopsy-proven malignancy. Surgical excision when clinically "
            "appropriate. Ensure tissue diagnosis concordance."
        ),
    }
    return recommendations.get(category, "Clinical correlation recommended.")


def classify_most_suspicious_finding(findings: List[Dict]) -> BiRADS:
    """Determine overall BI-RADS from the most suspicious finding.

    Each finding dict should have a 'birads_score' key with an integer 0-6.
    The overall BI-RADS is the maximum across all findings.

    Args:
        findings: List of finding dicts from mock inference.

    Returns:
        BiRADS enum for the highest-scoring finding.
    """
    birads_order = [
        BiRADS.CAT_0, BiRADS.CAT_1, BiRADS.CAT_2, BiRADS.CAT_3,
        BiRADS.CAT_4, BiRADS.CAT_5, BiRADS.CAT_6,
    ]
    max_score = 1  # default to BI-RADS 1 (negative)
    for finding in findings:
        score = finding.get("birads_score", 1)
        if score > max_score:
            max_score = score
    if max_score < 0 or max_score > 6:
        max_score = 1
    return birads_order[max_score]


# =====================================================================
# BreastBIRADSWorkflow
# =====================================================================


class BreastBIRADSWorkflow(BaseImagingWorkflow):
    """Breast imaging BI-RADS classification workflow.

    Pipeline:
        1. Load mammography DICOM (CC and MLO views) or breast MRI series
        2. Apply breast tissue segmentation (MONAI)
        3. Detect and characterize masses, calcifications, asymmetries
        4. Score each finding per ACR BI-RADS 5th Edition lexicon
        5. Assign overall BI-RADS category based on most suspicious finding
        6. Trigger cross-modal genomic query if BI-RADS >= 4

    Cross-modal trigger:
        BI-RADS 4+ fires BRCA1, BRCA2, PALB2, CHEK2 genomic queries
        via CrossModalTrigger to enrich findings with hereditary breast
        cancer susceptibility context.
    """

    WORKFLOW_NAME: str = "breast_birads"
    TARGET_LATENCY_SEC: float = 300.0
    MODALITY: str = "mammography"
    BODY_REGION: str = "breast"
    MODELS_USED: List[str] = ["monai_breast_segmentation"]

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
        """Load mammography DICOM and apply preprocessing.

        Steps:
            - Load DICOM mammography (CC + MLO views) or breast MRI series
            - Apply breast region segmentation to isolate parenchyma
            - Normalize pixel intensity for model input
            - Apply CLAHE enhancement for calcification visibility

        In production, uses monai.transforms:
            LoadImaged, EnsureChannelFirstd, ScaleIntensityRanged,
            breast-specific windowing and contrast enhancement.
        """
        logger.info(
            f"Preprocessing mammography from {input_path}: "
            "breast segmentation, intensity normalization, CLAHE enhancement"
        )
        return None

    # ------------------------------------------------------------------
    # Inference
    # ------------------------------------------------------------------

    def infer(self, preprocessed: Any) -> Dict:
        """Run breast lesion detection and characterization.

        In mock mode, returns realistic findings with BI-RADS lexicon
        descriptors. In real mode, would run MONAI breast segmentation
        followed by lesion-level feature extraction.
        """
        if self.mock_mode:
            return self._mock_inference()

        logger.warning("Real inference not yet implemented; using mock results")
        return self._mock_inference()

    def _mock_inference(self) -> Dict:
        """Return realistic mock breast imaging findings.

        Simulates a diagnostic mammography with two findings:
            - An irregular mass with spiculated margins in the upper outer
              quadrant of the right breast (BI-RADS 4C)
            - Amorphous calcifications in a grouped distribution in the
              left breast (BI-RADS 4A)
        """
        return {
            "breast_density": "heterogeneously_dense",  # ACR density C
            "density_category": "C",
            "findings": [
                {
                    "finding_type": "mass",
                    "laterality": "right",
                    "location": "upper outer quadrant",
                    "clock_position": "10 o'clock",
                    "depth": "middle third",
                    "distance_from_nipple_mm": 45.0,
                    "size_mm": 18.3,
                    "shape": "irregular",
                    "margin": "spiculated",
                    "density": "high",
                    "associated_calcifications": False,
                    "skin_retraction": False,
                    "architectural_distortion": True,
                    "detection_confidence": 0.95,
                    "birads_score": 4,  # BI-RADS 4C (suspicious)
                    "birads_subcategory": "4C",
                    "malignancy_probability": 0.72,
                },
                {
                    "finding_type": "calcifications",
                    "laterality": "left",
                    "location": "upper inner quadrant",
                    "clock_position": "1 o'clock",
                    "depth": "anterior third",
                    "distance_from_nipple_mm": 32.0,
                    "morphology": "amorphous",
                    "distribution": "grouped",
                    "extent_mm": 12.5,
                    "number_of_calcifications": 8,
                    "associated_mass": False,
                    "detection_confidence": 0.88,
                    "birads_score": 4,  # BI-RADS 4A (low suspicion)
                    "birads_subcategory": "4A",
                    "malignancy_probability": 0.08,
                },
            ],
        }

    # ------------------------------------------------------------------
    # Postprocessing
    # ------------------------------------------------------------------

    def postprocess(self, inference_result: Dict) -> WorkflowResult:
        """Classify BI-RADS 0-6 and determine severity.

        Applies ACR BI-RADS 5th Edition:
            - Each finding is scored individually
            - Overall BI-RADS is determined by the most suspicious finding
            - Cross-modal trigger fires if BI-RADS >= 4

        Args:
            inference_result: Raw inference dict with 'findings' list.

        Returns:
            WorkflowResult with BI-RADS classification and severity.
        """
        raw_findings = inference_result.get("findings", [])
        breast_density = inference_result.get("breast_density", "unknown")
        density_category = inference_result.get("density_category", "")

        if not raw_findings:
            return WorkflowResult(
                workflow_name=self.WORKFLOW_NAME,
                status=WorkflowStatus.COMPLETED,
                findings=[{
                    "category": FindingCategory.NORMAL.value,
                    "description": "No suspicious findings. Negative mammogram.",
                    "severity": FindingSeverity.NORMAL.value,
                    "birads": BiRADS.CAT_1.value,
                }],
                measurements={"breast_density_category": 0.0},
                classification=f"BI-RADS {BiRADS.CAT_1.value}",
                severity=FindingSeverity.NORMAL,
                nim_services_used=self.MODELS_USED,
            )

        # Determine overall BI-RADS from most suspicious finding
        overall_birads = classify_most_suspicious_finding(raw_findings)
        overall_severity = birads_to_severity(overall_birads)
        recommendation = birads_recommendation(overall_birads)

        # Build per-finding output
        processed_findings = []
        for finding in raw_findings:
            finding_type = finding.get("finding_type", "unspecified")
            laterality = finding.get("laterality", "unspecified")
            location = finding.get("location", "unspecified")
            birads_score = finding.get("birads_score", 1)
            birads_sub = finding.get("birads_subcategory", str(birads_score))

            # Map individual finding category
            if finding_type == "mass":
                category = FindingCategory.MASS.value
            elif finding_type == "calcifications":
                category = FindingCategory.CALCIFICATION.value
            elif finding_type == "architectural_distortion":
                category = FindingCategory.DISTORTION.value
            else:
                category = FindingCategory.LESION.value

            # Build description based on finding type
            if finding_type == "mass":
                shape = finding.get("shape", "")
                margin = finding.get("margin", "")
                size_mm = finding.get("size_mm", 0.0)
                description = (
                    f"{shape.title()} mass with {margin} margins in "
                    f"{laterality} breast {location}, {size_mm:.1f}mm, "
                    f"BI-RADS {birads_sub}"
                )
            elif finding_type == "calcifications":
                morphology = finding.get("morphology", "")
                distribution = finding.get("distribution", "")
                extent_mm = finding.get("extent_mm", 0.0)
                description = (
                    f"{morphology.title()} calcifications in {distribution} "
                    f"distribution, {laterality} breast {location}, "
                    f"extent {extent_mm:.1f}mm, BI-RADS {birads_sub}"
                )
            else:
                description = (
                    f"{finding_type.replace('_', ' ').title()} in "
                    f"{laterality} breast {location}, BI-RADS {birads_sub}"
                )

            # Individual finding severity
            birads_enum_order = [
                BiRADS.CAT_0, BiRADS.CAT_1, BiRADS.CAT_2, BiRADS.CAT_3,
                BiRADS.CAT_4, BiRADS.CAT_5, BiRADS.CAT_6,
            ]
            finding_birads = birads_enum_order[min(birads_score, 6)]
            finding_severity = birads_to_severity(finding_birads)

            processed_findings.append({
                "category": category,
                "description": description,
                "severity": finding_severity.value,
                "finding_type": finding_type,
                "laterality": laterality,
                "location": location,
                "birads": birads_sub,
                "recommendation": birads_recommendation(finding_birads),
                "detection_confidence": finding.get("detection_confidence", 0.0),
                "malignancy_probability": finding.get("malignancy_probability", 0.0),
            })

        # Aggregate measurements
        measurements = {
            "finding_count": float(len(raw_findings)),
            "breast_density_category": float(
                ord(density_category) - ord("A") + 1
            ) if density_category in ("A", "B", "C", "D") else 0.0,
        }

        for i, finding in enumerate(raw_findings):
            prefix = f"finding_{i + 1}"
            if "size_mm" in finding:
                measurements[f"{prefix}_size_mm"] = finding["size_mm"]
            if "extent_mm" in finding:
                measurements[f"{prefix}_extent_mm"] = finding["extent_mm"]
            measurements[f"{prefix}_confidence"] = finding.get(
                "detection_confidence", 0.0
            )

        # Determine cross-modal trigger eligibility
        cross_modal_trigger = overall_birads in (
            BiRADS.CAT_4, BiRADS.CAT_5, BiRADS.CAT_6,
        )

        if cross_modal_trigger:
            # Add trigger metadata to first finding
            processed_findings[0]["cross_modal_trigger"] = True
            processed_findings[0]["genomic_queries"] = [
                "BRCA1", "BRCA2", "PALB2", "CHEK2",
            ]

        return WorkflowResult(
            workflow_name=self.WORKFLOW_NAME,
            status=WorkflowStatus.COMPLETED,
            findings=processed_findings,
            measurements=measurements,
            classification=f"BI-RADS {overall_birads.value}",
            severity=overall_severity,
            nim_services_used=self.MODELS_USED,
        )
