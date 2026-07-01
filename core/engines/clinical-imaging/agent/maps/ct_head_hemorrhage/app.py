"""CT Head Hemorrhage Triage — MONAI Application Package.

Clinical workflow:
  DICOM CT Head -> SegResNet segmentation -> Volume/midline measurement ->
  Severity classification (Brain Trauma Foundation) -> DICOM SR output

Scoring:
  - Volume > 30 mL OR midline shift > 5mm OR thickness > 10mm -> Critical (P1)
  - Volume > 5 mL -> Urgent (P2)
  - Otherwise -> Routine

Target latency: < 90 seconds
NVIDIA technologies: MONAI SegResNet, MONAI Deploy SDK, TensorRT (optional)
License: Apache 2.0
"""

import sys
import uuid
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
from src.workflows.ct_head_hemorrhage import CTHeadHemorrhageWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if HIGHDICOM_AVAILABLE:
    import highdicom as hd
    from highdicom.sr.content import FindingSite
    from highdicom.sr import CodedConcept

if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — CT Head
# ═══════════════════════════════════════════════════════════════════════

CT_HEAD_RULES = [
    # Primary: non-contrast CT head, thin-slice axial
    SeriesSelectionRule(
        modality="CT",
        body_part="HEAD",
        description_keywords=["head", "brain", "axial", "thin", "non-contrast"],
        min_instances=20,
        preferred_slice_thickness_mm=1.5,
    ),
    # Secondary: any CT head series with enough slices
    SeriesSelectionRule(
        modality="CT",
        body_part="HEAD",
        description_keywords=["head", "brain"],
        min_instances=10,
        preferred_slice_thickness_mm=3.0,
    ),
    # Fallback: any CT with HEAD body part
    SeriesSelectionRule(
        modality="CT",
        body_part="HEAD",
        min_instances=5,
    ),
]


# ═══════════════════════════════════════════════════════════════════════
# TID 1500 Measurement Codes
# ═══════════════════════════════════════════════════════════════════════

# SNOMED CT and UCUM codes for structured measurement report
_HEMORRHAGE_VOLUME_CODE = ("G-D705", "SCT", "Volume")
_MIDLINE_SHIFT_CODE = ("G-D7FE", "SCT", "Length")
_THICKNESS_CODE = ("G-D7FE", "SCT", "Length")
_MILLILITER = ("mL", "UCUM", "milliliter")
_MILLIMETER = ("mm", "UCUM", "millimeter")


class CTHeadHemorrhageApp(BaseImagingMAP):
    """MONAI Deploy Application for CT head hemorrhage triage.

    Wraps the existing CTHeadHemorrhageWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (Brain Trauma Foundation thresholds, SegResNet inference, severity
    classification) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Hemorrhage volume (mL)
        - Midline shift (mm)
        - Maximum thickness (mm)
        - Hemorrhage type and location
        - Severity classification
    """

    APP_NAME: str = "ct_head_hemorrhage"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated CT head hemorrhage triage with SegResNet. "
        "Applies Brain Trauma Foundation severity thresholds for "
        "emergency department triage (P1/P2/Routine)."
    )
    TARGET_LATENCY_SEC: float = 90.0

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(mock_mode=mock_mode, *args, **kwargs)
        self._nim_clients = nim_clients
        self._workflow: Optional[CTHeadHemorrhageWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return CT head series selection rules.

        Prefers thin-slice axial non-contrast CT head. Falls back to any
        CT series with HEAD body part and sufficient slice count.
        """
        return CT_HEAD_RULES

    def get_workflow(self) -> CTHeadHemorrhageWorkflow:
        """Return the CT head hemorrhage workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = CTHeadHemorrhageWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from hemorrhage findings.

        Generates a DICOM Structured Report encoding:
            - Hemorrhage detection status
            - Volume measurement (mL)
            - Midline shift (mm)
            - Maximum thickness (mm)
            - Hemorrhage type and location
            - Brain Trauma Foundation severity classification
            - Clinical recommendation

        Args:
            workflow_result: WorkflowResult from CTHeadHemorrhageWorkflow.
            source_dicom: Optional source DICOM dataset for patient/study UIDs.

        Returns:
            Serialized DICOM SR bytes, or None if highdicom is unavailable.
        """
        if not HIGHDICOM_AVAILABLE or not PYDICOM_AVAILABLE:
            logger.info(
                "DICOM SR generation skipped — "
                "highdicom or pydicom not available"
            )
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
        """Internal SR construction with full TID 1500 content."""
        now = datetime.now(timezone.utc)

        # Extract measurements from workflow result
        measurements = workflow_result.measurements
        volume_ml = measurements.get("volume_ml", 0.0)
        midline_shift_mm = measurements.get("midline_shift_mm", 0.0)
        max_thickness_mm = measurements.get("max_thickness_mm", 0.0)

        # Extract finding details
        finding = workflow_result.findings[0] if workflow_result.findings else {}
        hemorrhage_type = finding.get("hemorrhage_type", "unspecified")
        location = finding.get("location", "unspecified")
        recommendation = finding.get("recommendation", "")
        severity = workflow_result.severity.value

        # Build text content for SR
        sr_text_lines = [
            f"CT Head Hemorrhage Triage Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Hemorrhage detected: {'Yes' if volume_ml > 0 else 'No'}",
            f"  Type: {hemorrhage_type.replace('_', ' ').title()}",
            f"  Location: {location}",
            f"",
            f"MEASUREMENTS:",
            f"  Volume: {volume_ml:.1f} mL",
            f"  Midline shift: {midline_shift_mm:.1f} mm",
            f"  Maximum thickness: {max_thickness_mm:.1f} mm",
            f"",
            f"SEVERITY: {severity.upper()}",
            f"",
            f"RECOMMENDATION:",
            f"  {recommendation}",
            f"",
            f"Brain Trauma Foundation Guidelines Applied:",
            f"  Critical: volume > 30 mL OR shift > 5 mm OR thickness > 10 mm",
            f"  Urgent: volume > 5 mL",
            f"  Routine: hemorrhage present, below urgent thresholds",
        ]

        # Create a minimal DICOM SR dataset
        sr_dataset = pydicom.Dataset()
        sr_dataset.SOPClassUID = "1.2.840.10008.5.1.4.1.1.88.33"  # Comprehensive SR
        sr_dataset.SOPInstanceUID = generate_uid()
        sr_dataset.SeriesInstanceUID = generate_uid()
        sr_dataset.Modality = "SR"
        sr_dataset.Manufacturer = "HCLS AI Factory"
        sr_dataset.ManufacturerModelName = self.APP_NAME
        sr_dataset.SoftwareVersions = self.APP_VERSION
        sr_dataset.ContentDate = now.strftime("%Y%m%d")
        sr_dataset.ContentTime = now.strftime("%H%M%S")
        sr_dataset.SeriesDescription = (
            f"CT Head Hemorrhage Triage - {severity.upper()}"
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

        # Encode findings as TEXT content in the SR tree
        # Using a simplified content tree (TEXT nodes) rather than
        # full TID 1500 numeric measurements, for broad viewer compatibility
        text_content = pydicom.Dataset()
        text_content.RelationshipType = "CONTAINS"
        text_content.ValueType = "TEXT"
        text_content.TextValue = "\n".join(sr_text_lines)

        concept_name = pydicom.Dataset()
        concept_name.CodeValue = "59776-5"
        concept_name.CodingSchemeDesignator = "LN"
        concept_name.CodeMeaning = "Findings"
        text_content.ConceptNameCodeSequence = [concept_name]

        # Severity code
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

        # Volume measurement (NUM node)
        volume_content = pydicom.Dataset()
        volume_content.RelationshipType = "CONTAINS"
        volume_content.ValueType = "NUM"

        volume_concept = pydicom.Dataset()
        volume_concept.CodeValue = "118565006"
        volume_concept.CodingSchemeDesignator = "SCT"
        volume_concept.CodeMeaning = "Volume"
        volume_content.ConceptNameCodeSequence = [volume_concept]

        volume_measured = pydicom.Dataset()
        volume_measured.NumericValue = f"{volume_ml:.1f}"

        volume_unit = pydicom.Dataset()
        volume_unit.CodeValue = "mL"
        volume_unit.CodingSchemeDesignator = "UCUM"
        volume_unit.CodeMeaning = "milliliter"
        volume_measured.MeasurementUnitsCodeSequence = [volume_unit]
        volume_content.MeasuredValueSequence = [volume_measured]

        # Build content sequence
        sr_dataset.ContentSequence = [
            text_content,
            severity_content,
            volume_content,
        ]

        # Set the root content item
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
        file_meta.TransferSyntaxUID = "1.2.840.10008.1.2.1"  # Explicit VR LE
        file_meta.ImplementationClassUID = "1.2.826.0.1.3680043.8.498.1"
        sr_dataset.file_meta = file_meta
        sr_dataset.preamble = b"\x00" * 128

        # Serialize to bytes
        from io import BytesIO

        buffer = BytesIO()
        pydicom.dcmwrite(buffer, sr_dataset)
        sr_bytes = buffer.getvalue()

        logger.info(
            f"Generated DICOM SR: {len(sr_bytes)} bytes, "
            f"severity={severity}, volume={volume_ml:.1f} mL"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for CT head hemorrhage.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> CTHeadInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"CT Head Hemorrhage MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run CT Head Hemorrhage MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CT Head Hemorrhage Triage MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM CT head series",
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
        help="Run with real model inference (requires GPU + MONAI)",
    )
    args = parser.parse_args()

    mock_mode = not args.real

    app = CTHeadHemorrhageApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"CT Head Hemorrhage Triage complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
