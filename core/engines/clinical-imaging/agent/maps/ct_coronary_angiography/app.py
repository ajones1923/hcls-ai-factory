"""CT Coronary Angiography — MONAI Application Package.

Clinical workflow:
  DICOM CTA Heart -> Coronary artery segmentation -> Stenosis grading ->
  Calcium scoring -> CAD-RADS classification -> DICOM SR output

Scoring:
  - CAD-RADS 0: No stenosis (0%)
  - CAD-RADS 1: Minimal stenosis (1-24%)
  - CAD-RADS 2: Mild stenosis (25-49%)
  - CAD-RADS 3: Moderate stenosis (50-69%)
  - CAD-RADS 4A: Severe stenosis (70-99%, 1-2 vessels)
  - CAD-RADS 4B: Severe stenosis (70-99%, 3 vessels or left main)
  - CAD-RADS 5: Total occlusion (100%)

Target latency: < 180 seconds
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
from src.workflows.ct_coronary_angiography import CTCoronaryAngiographyWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — CT Coronary Angiography
# ═══════════════════════════════════════════════════════════════════════

CT_CORONARY_RULES = [
    # Primary: ECG-gated CTA cardiac, thin-slice
    SeriesSelectionRule(
        modality="CT",
        body_part="CARDIAC",
        description_keywords=["cta", "coronary", "cardiac", "ecg", "gated", "angio"],
        min_instances=100,
        preferred_slice_thickness_mm=0.75,
    ),
    # Secondary: CT chest CTA protocol
    SeriesSelectionRule(
        modality="CT",
        body_part="CHEST",
        description_keywords=["cta", "coronary", "cardiac", "angio"],
        min_instances=50,
        preferred_slice_thickness_mm=1.0,
    ),
    # Tertiary: any cardiac CT with sufficient slices
    SeriesSelectionRule(
        modality="CT",
        body_part="CARDIAC",
        description_keywords=["cardiac", "heart"],
        min_instances=30,
        preferred_slice_thickness_mm=1.5,
    ),
    # Fallback: any CT with CARDIAC body part
    SeriesSelectionRule(
        modality="CT",
        body_part="CARDIAC",
        min_instances=20,
    ),
]


class CTCoronaryAngiographyApp(BaseImagingMAP):
    """MONAI Deploy Application for CT coronary angiography analysis.

    Wraps the existing CTCoronaryAngiographyWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (coronary artery segmentation, stenosis grading, calcium scoring,
    CAD-RADS classification) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Stenosis percentage per vessel segment
        - Agatston calcium score
        - CAD-RADS classification
        - Per-vessel findings and recommendations
    """

    APP_NAME: str = "ct_coronary_angiography"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated CT coronary angiography analysis with CAD-RADS "
        "classification. Measures stenosis percentage, Agatston calcium "
        "score, and generates per-vessel segment reports."
    )
    TARGET_LATENCY_SEC: float = 180.0

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(mock_mode=mock_mode, *args, **kwargs)
        self._nim_clients = nim_clients
        self._workflow: Optional[CTCoronaryAngiographyWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return CT coronary angiography series selection rules.

        Prefers ECG-gated CTA cardiac acquisitions with thin-slice
        reconstruction. Falls back to any CT series with CARDIAC body
        part and sufficient slice count.
        """
        return CT_CORONARY_RULES

    def get_workflow(self) -> CTCoronaryAngiographyWorkflow:
        """Return the CT coronary angiography workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = CTCoronaryAngiographyWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from coronary findings.

        Generates a DICOM Structured Report encoding:
            - Maximum stenosis percentage
            - Agatston calcium score
            - CAD-RADS classification
            - Per-vessel segment stenosis
            - Clinical recommendation

        Args:
            workflow_result: WorkflowResult from CTCoronaryAngiographyWorkflow.
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
        """Internal SR construction for coronary angiography findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        max_stenosis_pct = measurements.get("max_stenosis_pct", 0.0)
        calcium_score = measurements.get("calcium_score", 0.0)
        cad_rads = measurements.get("cad_rads", "0")

        # Per-vessel findings
        vessel_lines = []
        for finding in workflow_result.findings:
            vessel = finding.get("vessel", "unspecified")
            segment = finding.get("segment", "")
            stenosis = finding.get("stenosis_pct", 0.0)
            plaque_type = finding.get("plaque_type", "none")
            vessel_lines.append(
                f"  {vessel} {segment}: stenosis={stenosis:.0f}%, plaque={plaque_type}"
            )

        recommendation = workflow_result.findings[0].get("recommendation", "") if workflow_result.findings else ""

        sr_text_lines = [
            f"CT Coronary Angiography Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Maximum stenosis: {max_stenosis_pct:.0f}%",
            f"  Agatston calcium score: {calcium_score:.0f}",
            f"  CAD-RADS: {cad_rads}",
            f"",
        ]

        if vessel_lines:
            sr_text_lines.extend([
                f"VESSEL ANALYSIS:",
                *vessel_lines,
                f"",
            ])

        sr_text_lines.extend([
            f"SEVERITY: {severity.upper()}",
            f"",
            f"RECOMMENDATION:",
            f"  {recommendation}",
            f"",
            f"CAD-RADS Classification Applied:",
            f"  0: No stenosis (0%)",
            f"  1: Minimal (1-24%)",
            f"  2: Mild (25-49%)",
            f"  3: Moderate (50-69%)",
            f"  4A: Severe (70-99%, 1-2 vessels)",
            f"  4B: Severe (70-99%, 3 vessels or left main)",
            f"  5: Total occlusion (100%)",
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
            f"CT Coronary Angiography - CAD-RADS {cad_rads}"
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

        # Stenosis NUM content
        stenosis_content = pydicom.Dataset()
        stenosis_content.RelationshipType = "CONTAINS"
        stenosis_content.ValueType = "NUM"

        stenosis_concept = pydicom.Dataset()
        stenosis_concept.CodeValue = "373157009"
        stenosis_concept.CodingSchemeDesignator = "SCT"
        stenosis_concept.CodeMeaning = "Stenosis percentage"
        stenosis_content.ConceptNameCodeSequence = [stenosis_concept]

        stenosis_measured = pydicom.Dataset()
        stenosis_measured.NumericValue = f"{max_stenosis_pct:.0f}"

        stenosis_unit = pydicom.Dataset()
        stenosis_unit.CodeValue = "%"
        stenosis_unit.CodingSchemeDesignator = "UCUM"
        stenosis_unit.CodeMeaning = "percent"
        stenosis_measured.MeasurementUnitsCodeSequence = [stenosis_unit]
        stenosis_content.MeasuredValueSequence = [stenosis_measured]

        # Build content sequence
        sr_dataset.ContentSequence = [
            text_content,
            severity_content,
            stenosis_content,
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
            f"severity={severity}, CAD-RADS={cad_rads}, "
            f"max_stenosis={max_stenosis_pct:.0f}%"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for CT coronary angiography.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> CoronaryInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"CT Coronary Angiography MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run CT Coronary Angiography MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CT Coronary Angiography MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM CTA cardiac series",
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

    app = CTCoronaryAngiographyApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"CT Coronary Angiography complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
