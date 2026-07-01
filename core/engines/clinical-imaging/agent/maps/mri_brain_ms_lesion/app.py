"""MRI Brain MS Lesion Detection — MONAI Application Package.

Clinical workflow:
  DICOM MRI Brain (FLAIR/T2) -> White matter lesion segmentation ->
  Lesion counting and volumetry -> New/enlarging detection ->
  MS activity classification -> DICOM SR output

Scoring:
  - Active: new or enlarging lesions since prior study
  - Stable: no new or enlarging lesions
  - Lesion burden: total volume and count

Target latency: < 120 seconds
NVIDIA technologies: MONAI wholeBrainSeg_Large_UNEST, MONAI Deploy SDK, TensorRT (optional)
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
from src.workflows.mri_brain_ms_lesion import MRIBrainMSLesionWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — MRI Brain MS Lesion
# ═══════════════════════════════════════════════════════════════════════

MRI_BRAIN_MS_RULES = [
    # Primary: FLAIR sequence (preferred for MS lesion detection)
    SeriesSelectionRule(
        modality="MR",
        body_part="BRAIN",
        description_keywords=["flair", "t2_flair", "dark_fluid", "brain"],
        min_instances=20,
        preferred_slice_thickness_mm=1.0,
    ),
    # Secondary: T2-weighted (fallback for MS lesion detection)
    SeriesSelectionRule(
        modality="MR",
        body_part="BRAIN",
        description_keywords=["t2", "t2w", "brain"],
        min_instances=20,
        preferred_slice_thickness_mm=1.5,
    ),
    # Tertiary: HEAD body part with FLAIR
    SeriesSelectionRule(
        modality="MR",
        body_part="HEAD",
        description_keywords=["flair", "brain", "ms"],
        min_instances=15,
        preferred_slice_thickness_mm=2.0,
    ),
    # Fallback: any MR brain/head series
    SeriesSelectionRule(
        modality="MR",
        body_part="BRAIN",
        min_instances=10,
    ),
]


class MRIBrainMSLesionApp(BaseImagingMAP):
    """MONAI Deploy Application for MRI brain MS lesion detection.

    Wraps the existing MRIBrainMSLesionWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (white matter lesion segmentation, volumetry, new/enlarging detection,
    MS activity classification) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Total lesion count
        - Total lesion volume (mL)
        - New lesion count (vs prior)
        - Enlarging lesion count (vs prior)
        - Activity classification (Active/Stable)
        - Per-lesion location and volume
    """

    APP_NAME: str = "mri_brain_ms_lesion"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated MRI brain MS lesion detection and volumetry. "
        "Segments white matter lesions on FLAIR/T2, quantifies lesion "
        "burden, and classifies disease activity for MS monitoring."
    )
    TARGET_LATENCY_SEC: float = 120.0

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(mock_mode=mock_mode, *args, **kwargs)
        self._nim_clients = nim_clients
        self._workflow: Optional[MRIBrainMSLesionWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return MRI brain MS lesion series selection rules.

        Prefers FLAIR sequences for optimal white matter lesion contrast.
        Falls back to T2-weighted sequences. Requires minimum 10 slices
        for volumetric analysis.
        """
        return MRI_BRAIN_MS_RULES

    def get_workflow(self) -> MRIBrainMSLesionWorkflow:
        """Return the MRI brain MS lesion workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = MRIBrainMSLesionWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from MS lesion findings.

        Generates a DICOM Structured Report encoding:
            - Total lesion count and volume
            - New and enlarging lesion counts
            - Activity classification (Active/Stable)
            - Per-lesion location and measurements
            - Clinical recommendation

        Args:
            workflow_result: WorkflowResult from MRIBrainMSLesionWorkflow.
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
        """Internal SR construction for MS lesion findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        lesion_count = int(measurements.get("lesion_count", 0))
        total_volume_ml = measurements.get("total_volume_ml", 0.0)
        new_lesion_count = int(measurements.get("new_lesion_count", 0))
        enlarging_count = int(measurements.get("enlarging_lesion_count", 0))
        activity = measurements.get("activity_classification", "stable")

        # Per-lesion details
        lesion_lines = []
        for i, finding in enumerate(workflow_result.findings):
            loc = finding.get("location", "unspecified")
            vol = finding.get("volume_ml", 0.0)
            status = finding.get("status", "stable")
            lesion_lines.append(
                f"  Lesion {i + 1}: {vol:.2f} mL, location={loc}, status={status}"
            )

        sr_text_lines = [
            f"MRI Brain MS Lesion Detection Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Total lesion count: {lesion_count}",
            f"  Total lesion volume: {total_volume_ml:.2f} mL",
            f"  New lesions: {new_lesion_count}",
            f"  Enlarging lesions: {enlarging_count}",
            f"  Activity classification: {activity.upper()}",
            f"",
        ]

        if lesion_lines:
            sr_text_lines.extend([
                f"LESION DETAILS:",
                *lesion_lines,
                f"",
            ])

        sr_text_lines.extend([
            f"SEVERITY: {severity.upper()}",
            f"",
            f"MS Activity Classification:",
            f"  Active: new or enlarging lesions detected",
            f"  Stable: no new or enlarging lesions since prior",
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
            f"MRI Brain MS Lesion - {activity.upper()}"
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

        # Lesion volume NUM content
        volume_content = pydicom.Dataset()
        volume_content.RelationshipType = "CONTAINS"
        volume_content.ValueType = "NUM"

        volume_concept = pydicom.Dataset()
        volume_concept.CodeValue = "118565006"
        volume_concept.CodingSchemeDesignator = "SCT"
        volume_concept.CodeMeaning = "Volume"
        volume_content.ConceptNameCodeSequence = [volume_concept]

        volume_measured = pydicom.Dataset()
        volume_measured.NumericValue = f"{total_volume_ml:.2f}"

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
            f"severity={severity}, activity={activity}, "
            f"lesions={lesion_count}"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for MRI brain MS lesion detection.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> MSLesionInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"MRI Brain MS Lesion MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run MRI Brain MS Lesion MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="MRI Brain MS Lesion Detection MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM MRI brain FLAIR/T2 series",
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

    app = MRIBrainMSLesionApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"MRI Brain MS Lesion Detection complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
