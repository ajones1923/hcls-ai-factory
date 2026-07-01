"""Breast BI-RADS Classification — MONAI Application Package.

Clinical workflow:
  DICOM Mammography/MRI Breast -> Mass/calcification detection ->
  Feature characterization -> BI-RADS classification -> DICOM SR output

Scoring:
  - BI-RADS 0: Incomplete — need additional imaging
  - BI-RADS 1: Negative
  - BI-RADS 2: Benign
  - BI-RADS 3: Probably benign (< 2% malignancy)
  - BI-RADS 4: Suspicious (4A: 2-10%, 4B: 10-50%, 4C: 50-95%)
  - BI-RADS 5: Highly suggestive of malignancy (>= 95%)
  - BI-RADS 6: Known biopsy-proven malignancy

Target latency: < 60 seconds
NVIDIA technologies: MONAI Deploy SDK, TensorRT (optional)
License: Apache 2.0
"""

import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

# Resolve project root for imports when running as a MAP
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from maps.common.base_imaging_map import (
    BaseImagingMAP,
    SeriesSelectionRule,
    HIGHDICOM_AVAILABLE,
    PYDICOM_AVAILABLE,
)
from src.models import FindingSeverity, WorkflowResult
from src.workflows.breast_birads import BreastBIRADSWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — Breast BI-RADS
# ═══════════════════════════════════════════════════════════════════════

BREAST_BIRADS_RULES = [
    # Primary: digital mammography MLO view
    SeriesSelectionRule(
        modality="MG",
        body_part="BREAST",
        description_keywords=["mlo", "mediolateral", "oblique", "mammogram"],
        min_instances=1,
    ),
    # Secondary: digital mammography CC view
    SeriesSelectionRule(
        modality="MG",
        body_part="BREAST",
        description_keywords=["cc", "craniocaudal", "mammogram"],
        min_instances=1,
    ),
    # Tertiary: breast MRI (dynamic contrast-enhanced)
    SeriesSelectionRule(
        modality="MR",
        body_part="BREAST",
        description_keywords=["breast", "dce", "dynamic", "contrast", "axial"],
        min_instances=10,
        preferred_slice_thickness_mm=1.5,
    ),
    # Fallback: any mammography
    SeriesSelectionRule(
        modality="MG",
        body_part="BREAST",
        min_instances=1,
    ),
    # Fallback: any breast MRI
    SeriesSelectionRule(
        modality="MR",
        body_part="BREAST",
        min_instances=5,
    ),
]


class BreastBIRADSApp(BaseImagingMAP):
    """MONAI Deploy Application for breast BI-RADS classification.

    Wraps the existing BreastBIRADSWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (mass detection, calcification characterization, BI-RADS
    classification) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Mass description (shape, margin, density)
        - Calcification morphology and distribution
        - BI-RADS classification
        - Laterality and clock-face location
        - Clinical recommendation (biopsy, follow-up, routine)
    """

    APP_NAME: str = "breast_birads"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated breast imaging analysis with BI-RADS classification. "
        "Detects masses and calcifications on mammography and breast MRI, "
        "characterizes features, and assigns BI-RADS assessment category."
    )
    TARGET_LATENCY_SEC: float = 60.0

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(mock_mode=mock_mode, *args, **kwargs)
        self._nim_clients = nim_clients
        self._workflow: Optional[BreastBIRADSWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return breast BI-RADS series selection rules.

        Accepts mammography (MG) with MLO and CC views, and breast MRI
        with dynamic contrast-enhanced sequences. Both modalities are
        supported for BI-RADS assessment.
        """
        return BREAST_BIRADS_RULES

    def get_workflow(self) -> BreastBIRADSWorkflow:
        """Return the breast BI-RADS workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = BreastBIRADSWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from BI-RADS findings.

        Generates a DICOM Structured Report encoding:
            - Mass description (shape, margin, density)
            - Calcification morphology and distribution
            - BI-RADS assessment category
            - Laterality and location
            - Clinical recommendation

        Args:
            workflow_result: WorkflowResult from BreastBIRADSWorkflow.
            source_dicom: Optional source DICOM dataset for patient/study UIDs.

        Returns:
            Serialized DICOM SR bytes, or None if pydicom is unavailable.
        """
        if not PYDICOM_AVAILABLE:
            logger.info("DICOM SR generation skipped — pydicom not available")
            return None

        try:
            return self._build_sr_bytes(workflow_result, source_dicom)
        except Exception as e:
            logger.error(f"DICOM SR generation failed: {e}")
            return None

    def _build_sr_bytes(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any],
    ) -> bytes:
        """Internal SR construction for BI-RADS findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        birads = measurements.get("birads", "0")
        mass_count = int(measurements.get("mass_count", 0))
        calcification_count = int(measurements.get("calcification_count", 0))

        # Per-finding details
        finding_lines = []
        for i, finding in enumerate(workflow_result.findings):
            finding_type = finding.get("type", "mass")
            birads_score = finding.get("birads", "0")
            laterality = finding.get("laterality", "unspecified")
            location = finding.get("location", "unspecified")

            if finding_type == "mass":
                shape = finding.get("shape", "unspecified")
                margin = finding.get("margin", "unspecified")
                density = finding.get("density", "unspecified")
                finding_lines.append(
                    f"  Finding {i + 1} (mass): BI-RADS {birads_score}, "
                    f"{laterality}, {location}, shape={shape}, "
                    f"margin={margin}, density={density}"
                )
            else:
                morphology = finding.get("morphology", "unspecified")
                distribution = finding.get("distribution", "unspecified")
                finding_lines.append(
                    f"  Finding {i + 1} (calcification): BI-RADS {birads_score}, "
                    f"{laterality}, {location}, morphology={morphology}, "
                    f"distribution={distribution}"
                )

        recommendation = workflow_result.findings[0].get("recommendation", "") if workflow_result.findings else ""

        sr_text_lines = [
            f"Breast BI-RADS Assessment Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Masses detected: {mass_count}",
            f"  Calcification clusters: {calcification_count}",
            f"  Overall BI-RADS: {birads}",
            f"",
        ]

        if finding_lines:
            sr_text_lines.extend([
                f"FINDING DETAILS:",
                *finding_lines,
                f"",
            ])

        sr_text_lines.extend([
            f"SEVERITY: {severity.upper()}",
            f"",
            f"RECOMMENDATION:",
            f"  {recommendation}",
            f"",
            f"BI-RADS Classification Applied:",
            f"  0: Incomplete — additional imaging needed",
            f"  1: Negative",
            f"  2: Benign",
            f"  3: Probably benign (< 2% malignancy)",
            f"  4: Suspicious (4A: 2-10%, 4B: 10-50%, 4C: 50-95%)",
            f"  5: Highly suggestive of malignancy (>= 95%)",
            f"  6: Known biopsy-proven malignancy",
        ])

        # Create DICOM SR dataset
        sr_dataset = pydicom.Dataset()
        sr_dataset.SOPClassUID = "1.2.840.10008.5.1.4.1.1.88.33"
        sr_dataset.SOPInstanceUID = generate_uid()
        sr_dataset.SeriesInstanceUID = generate_uid()
        sr_dataset.Modality = "SR"
        sr_dataset.Manufacturer = "HCLS AI Factory"
        sr_dataset.ManufacturerModelName = self.APP_NAME
        sr_dataset.SoftwareVersions = self.APP_VERSION
        sr_dataset.ContentDate = now.strftime("%Y%m%d")
        sr_dataset.ContentTime = now.strftime("%H%M%S")
        sr_dataset.SeriesDescription = (
            f"Breast BI-RADS Assessment - BI-RADS {birads}"
        )

        # Copy patient/study context from source DICOM
        if source_dicom is not None:
            for tag_name in [
                "PatientID",
                "PatientName",
                "PatientBirthDate",
                "PatientSex",
                "StudyInstanceUID",
                "StudyDate",
                "StudyTime",
                "AccessionNumber",
                "ReferringPhysicianName",
            ]:
                if hasattr(source_dicom, tag_name):
                    setattr(sr_dataset, tag_name, getattr(source_dicom, tag_name))
        else:
            sr_dataset.PatientID = "UNKNOWN"
            sr_dataset.PatientName = "UNKNOWN"
            sr_dataset.StudyInstanceUID = generate_uid()

        # Findings TEXT content
        text_content = pydicom.Dataset()
        text_content.RelationshipType = "CONTAINS"
        text_content.ValueType = "TEXT"
        text_content.TextValue = "\n".join(sr_text_lines)

        concept_name = pydicom.Dataset()
        concept_name.CodeValue = "59776-5"
        concept_name.CodingSchemeDesignator = "LN"
        concept_name.CodeMeaning = "Findings"
        text_content.ConceptNameCodeSequence = [concept_name]

        # Severity CODE content
        severity_content = pydicom.Dataset()
        severity_content.RelationshipType = "CONTAINS"
        severity_content.ValueType = "CODE"

        severity_concept = pydicom.Dataset()
        severity_concept.CodeValue = "246112005"
        severity_concept.CodingSchemeDesignator = "SCT"
        severity_concept.CodeMeaning = "Severity"
        severity_content.ConceptNameCodeSequence = [severity_concept]

        severity_code = pydicom.Dataset()
        severity_map = {
            "critical": ("24484000", "Critical"),
            "urgent": ("103391001", "Urgent"),
            "significant": ("371924009", "Significant"),
            "routine": ("17621005", "Routine"),
            "normal": ("17621005", "Normal"),
        }
        code_value, code_meaning = severity_map.get(severity, ("17621005", "Routine"))
        severity_code.CodeValue = code_value
        severity_code.CodingSchemeDesignator = "SCT"
        severity_code.CodeMeaning = code_meaning
        severity_content.ConceptCodeSequence = [severity_code]

        # Build content sequence
        sr_dataset.ContentSequence = [
            text_content,
            severity_content,
        ]

        # Root content item
        root_concept = pydicom.Dataset()
        root_concept.CodeValue = "126000"
        root_concept.CodingSchemeDesignator = "DCM"
        root_concept.CodeMeaning = "Imaging Measurement Report"
        sr_dataset.ConceptNameCodeSequence = [root_concept]
        sr_dataset.ValueType = "CONTAINER"
        sr_dataset.ContinuityOfContent = "SEPARATE"
        sr_dataset.CompletionFlag = "COMPLETE"
        sr_dataset.VerificationFlag = "UNVERIFIED"

        # File meta
        sr_dataset.is_little_endian = True
        sr_dataset.is_implicit_VR = False
        file_meta = pydicom.Dataset()
        file_meta.MediaStorageSOPClassUID = sr_dataset.SOPClassUID
        file_meta.MediaStorageSOPInstanceUID = sr_dataset.SOPInstanceUID
        file_meta.TransferSyntaxUID = "1.2.840.10008.1.2.1"
        file_meta.ImplementationClassUID = "1.2.826.0.1.3680043.8.498.1"
        sr_dataset.file_meta = file_meta
        sr_dataset.preamble = b"\x00" * 128

        # Serialize
        from io import BytesIO

        buffer = BytesIO()
        pydicom.dcmwrite(buffer, sr_dataset)
        sr_bytes = buffer.getvalue()

        logger.info(
            f"Generated DICOM SR: {len(sr_bytes)} bytes, "
            f"severity={severity}, BI-RADS={birads}, "
            f"masses={mass_count}, calcs={calcification_count}"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for breast BI-RADS.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> BreastInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"Breast BI-RADS MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run Breast BI-RADS MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Breast BI-RADS Classification MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM mammography or breast MRI",
    )
    parser.add_argument(
        "--output", "-o",
        default="/output",
        help="Output directory for results and DICOM SR",
    )
    parser.add_argument(
        "--mock",
        action="store_true",
        default=True,
        help="Run in mock mode (no GPU required)",
    )
    parser.add_argument(
        "--real",
        action="store_true",
        help="Run with real model inference (requires GPU)",
    )
    args = parser.parse_args()

    mock_mode = not args.real

    app = BreastBIRADSApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"Breast BI-RADS Classification complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
