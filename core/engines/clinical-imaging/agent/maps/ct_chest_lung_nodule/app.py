"""CT Chest Lung Nodule Detection — MONAI Application Package.

Clinical workflow:
  DICOM CT Chest -> Lung nodule detection -> Size/volume measurement ->
  Lung-RADS classification -> DICOM SR output

Scoring:
  - Lung-RADS 1: No nodules / definitely benign
  - Lung-RADS 2: Benign appearance, < 6mm
  - Lung-RADS 3: Probably benign, 6-8mm
  - Lung-RADS 4A: Suspicious, 8-15mm -> cross-modal trigger
  - Lung-RADS 4B: Very suspicious, >= 15mm -> cross-modal trigger

Target latency: < 120 seconds
NVIDIA technologies: MONAI lung_nodule_ct_detection bundle, MONAI Deploy SDK, TensorRT (optional)
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
from src.workflows.ct_chest_lung_nodule import CTChestLungNoduleWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — CT Chest Lung Nodule
# ═══════════════════════════════════════════════════════════════════════

CT_CHEST_LUNG_NODULE_RULES = [
    # Primary: thin-slice axial CT chest with lung kernel
    SeriesSelectionRule(
        modality="CT",
        body_part="CHEST",
        description_keywords=["chest", "lung", "thin", "axial", "lung kernel", "b60"],
        min_instances=50,
        preferred_slice_thickness_mm=1.25,
    ),
    # Secondary: thin-slice CT chest without lung kernel
    SeriesSelectionRule(
        modality="CT",
        body_part="CHEST",
        description_keywords=["chest", "lung", "thin", "axial"],
        min_instances=30,
        preferred_slice_thickness_mm=1.5,
    ),
    # Tertiary: any CT chest with sufficient slices
    SeriesSelectionRule(
        modality="CT",
        body_part="CHEST",
        description_keywords=["chest"],
        min_instances=20,
        preferred_slice_thickness_mm=3.0,
    ),
    # Fallback: any CT with CHEST body part
    SeriesSelectionRule(
        modality="CT",
        body_part="CHEST",
        min_instances=10,
    ),
]


class CTChestLungNoduleApp(BaseImagingMAP):
    """MONAI Deploy Application for CT chest lung nodule detection.

    Wraps the existing CTChestLungNoduleWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (nodule detection, size measurement, Lung-RADS classification,
    volume doubling time estimation) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Nodule count and locations
        - Nodule sizes (long axis, short axis, volume)
        - Lung-RADS classification per nodule
        - Volume doubling time estimate
        - Cross-modal trigger recommendation for Lung-RADS 4A+
    """

    APP_NAME: str = "ct_chest_lung_nodule"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated CT lung nodule detection with Lung-RADS classification. "
        "Measures nodule size and volume, classifies per Lung-RADS v2022, "
        "and triggers cross-modal review for Lung-RADS 4A+ findings."
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
        self._workflow: Optional[CTChestLungNoduleWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return CT chest lung nodule series selection rules.

        Prefers thin-slice axial CT chest with lung kernel reconstruction.
        Falls back to any CT series with CHEST body part and sufficient
        slice count for volumetric nodule analysis.
        """
        return CT_CHEST_LUNG_NODULE_RULES

    def get_workflow(self) -> CTChestLungNoduleWorkflow:
        """Return the CT chest lung nodule workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = CTChestLungNoduleWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from lung nodule findings.

        Generates a DICOM Structured Report encoding:
            - Nodule detection count
            - Per-nodule size measurements (long axis, short axis, volume)
            - Lung-RADS classification per nodule
            - Volume doubling time estimate
            - Cross-modal trigger recommendation
            - Overall severity classification

        Args:
            workflow_result: WorkflowResult from CTChestLungNoduleWorkflow.
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
        """Internal SR construction for lung nodule findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        nodule_count = int(measurements.get("nodule_count", 0))
        lung_rads = measurements.get("lung_rads", "1")
        max_nodule_size_mm = measurements.get("max_nodule_size_mm", 0.0)
        total_volume_mm3 = measurements.get("total_volume_mm3", 0.0)
        volume_doubling_time_days = measurements.get("volume_doubling_time_days", 0.0)

        # Per-nodule findings
        nodule_lines = []
        for i, finding in enumerate(workflow_result.findings):
            size = finding.get("size_mm", 0.0)
            loc = finding.get("location", "unspecified")
            rads = finding.get("lung_rads", "N/A")
            nodule_lines.append(
                f"  Nodule {i + 1}: {size:.1f} mm, location={loc}, Lung-RADS {rads}"
            )

        cross_modal = any(
            finding.get("cross_modal_trigger", False)
            for finding in workflow_result.findings
        )

        sr_text_lines = [
            f"CT Chest Lung Nodule Detection Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Nodule count: {nodule_count}",
            f"  Maximum nodule size: {max_nodule_size_mm:.1f} mm",
            f"  Total nodule volume: {total_volume_mm3:.1f} mm³",
            f"  Overall Lung-RADS: {lung_rads}",
            f"",
        ]

        if nodule_lines:
            sr_text_lines.extend([
                f"NODULE DETAILS:",
                *nodule_lines,
                f"",
            ])

        if volume_doubling_time_days > 0:
            sr_text_lines.append(
                f"  Volume doubling time: {volume_doubling_time_days:.0f} days"
            )
            sr_text_lines.append(f"")

        sr_text_lines.extend([
            f"SEVERITY: {severity.upper()}",
            f"CROSS-MODAL TRIGGER: {'YES — PET/CT or biopsy recommended' if cross_modal else 'No'}",
            f"",
            f"Lung-RADS v2022 Classification Applied:",
            f"  1: No nodules / definitely benign",
            f"  2: Benign appearance, < 6 mm",
            f"  3: Probably benign, 6-8 mm",
            f"  4A: Suspicious, 8-15 mm",
            f"  4B: Very suspicious, >= 15 mm",
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
            f"CT Lung Nodule - Lung-RADS {lung_rads}"
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

        # Nodule size NUM content
        size_content = pydicom.Dataset()
        size_content.RelationshipType = "CONTAINS"
        size_content.ValueType = "NUM"

        size_concept = pydicom.Dataset()
        size_concept.CodeValue = "246120007"
        size_concept.CodingSchemeDesignator = "SCT"
        size_concept.CodeMeaning = "Nodule size"
        size_content.ConceptNameCodeSequence = [size_concept]

        size_measured = pydicom.Dataset()
        size_measured.NumericValue = f"{max_nodule_size_mm:.1f}"

        size_unit = pydicom.Dataset()
        size_unit.CodeValue = "mm"
        size_unit.CodingSchemeDesignator = "UCUM"
        size_unit.CodeMeaning = "millimeter"
        size_measured.MeasurementUnitsCodeSequence = [size_unit]
        size_content.MeasuredValueSequence = [size_measured]

        # Build content sequence
        sr_dataset.ContentSequence = [
            text_content,
            severity_content,
            size_content,
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
            f"severity={severity}, Lung-RADS={lung_rads}, "
            f"nodules={nodule_count}"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for CT lung nodule detection.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> LungNoduleInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"CT Chest Lung Nodule MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run CT Chest Lung Nodule MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CT Chest Lung Nodule Detection MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM CT chest series",
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

    app = CTChestLungNoduleApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"CT Chest Lung Nodule Detection complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
