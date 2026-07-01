"""Thyroid TI-RADS Classification — MONAI Application Package.

Clinical workflow:
  DICOM Ultrasound Thyroid -> Nodule detection -> Feature characterization ->
  ACR TI-RADS point scoring -> FNA recommendation -> DICOM SR output

Scoring (ACR TI-RADS):
  - TR1 (0 points): Benign — no FNA
  - TR2 (2 points): Not suspicious — no FNA
  - TR3 (3 points): Mildly suspicious — FNA >= 25mm, follow-up >= 15mm
  - TR4 (4-6 points): Moderately suspicious — FNA >= 15mm, follow-up >= 10mm
  - TR5 (>= 7 points): Highly suspicious — FNA >= 10mm, follow-up >= 5mm

Target latency: < 45 seconds
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
from src.workflows.thyroid_tirads import ThyroidTIRADSWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — Thyroid TI-RADS
# ═══════════════════════════════════════════════════════════════════════

THYROID_TIRADS_RULES = [
    # Primary: ultrasound thyroid
    SeriesSelectionRule(
        modality="US",
        body_part="THYROID",
        description_keywords=["thyroid", "neck", "transverse", "longitudinal"],
        min_instances=1,
    ),
    # Secondary: ultrasound neck
    SeriesSelectionRule(
        modality="US",
        body_part="NECK",
        description_keywords=["thyroid", "nodule"],
        min_instances=1,
    ),
    # Fallback: any US with THYROID body part
    SeriesSelectionRule(
        modality="US",
        body_part="THYROID",
        min_instances=1,
    ),
    # Fallback: any US with NECK body part
    SeriesSelectionRule(
        modality="US",
        body_part="NECK",
        min_instances=1,
    ),
]


class ThyroidTIRADSApp(BaseImagingMAP):
    """MONAI Deploy Application for thyroid TI-RADS classification.

    Wraps the existing ThyroidTIRADSWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (nodule detection, feature characterization, ACR TI-RADS point
    scoring, FNA recommendation) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Nodule size (3 dimensions)
        - Composition (cystic, spongiform, mixed, solid)
        - Echogenicity (anechoic, hyper, iso, hypo, very hypo)
        - Shape (wider-than-tall, taller-than-wide)
        - Margin (smooth, ill-defined, lobulated, irregular, ETE)
        - Echogenic foci (none, large comet-tail, macrocalcifications,
          peripheral rim, punctate echogenic foci)
        - ACR TI-RADS score and category
        - FNA recommendation
    """

    APP_NAME: str = "thyroid_tirads"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated thyroid nodule classification with ACR TI-RADS "
        "point scoring. Characterizes nodule features on ultrasound and "
        "generates FNA biopsy recommendation per ACR guidelines."
    )
    TARGET_LATENCY_SEC: float = 45.0

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(mock_mode=mock_mode, *args, **kwargs)
        self._nim_clients = nim_clients
        self._workflow: Optional[ThyroidTIRADSWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return thyroid TI-RADS series selection rules.

        Accepts ultrasound (US) modality with THYROID or NECK body part.
        Both transverse and longitudinal views are useful for nodule
        characterization.
        """
        return THYROID_TIRADS_RULES

    def get_workflow(self) -> ThyroidTIRADSWorkflow:
        """Return the thyroid TI-RADS workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = ThyroidTIRADSWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from TI-RADS findings.

        Generates a DICOM Structured Report encoding:
            - Nodule size and dimensions
            - Composition, echogenicity, shape, margin, echogenic foci
            - ACR TI-RADS point total and category
            - FNA recommendation with size thresholds
            - Clinical recommendation

        Args:
            workflow_result: WorkflowResult from ThyroidTIRADSWorkflow.
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
        """Internal SR construction for TI-RADS findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        tirads_score = int(measurements.get("tirads_points", 0))
        tirads_category = measurements.get("tirads_category", "TR1")
        nodule_size_mm = measurements.get("nodule_size_mm", 0.0)
        fna_recommended = measurements.get("fna_recommended", False)

        # Per-nodule details
        nodule_lines = []
        for i, finding in enumerate(workflow_result.findings):
            size = finding.get("size_mm", 0.0)
            composition = finding.get("composition", "unspecified")
            echogenicity = finding.get("echogenicity", "unspecified")
            shape = finding.get("shape", "unspecified")
            margin = finding.get("margin", "unspecified")
            foci = finding.get("echogenic_foci", "none")
            points = finding.get("tirads_points", 0)
            nodule_lines.append(
                f"  Nodule {i + 1}: {size:.1f} mm, {tirads_category}"
            )
            nodule_lines.append(
                f"    Composition: {composition}, Echogenicity: {echogenicity}"
            )
            nodule_lines.append(
                f"    Shape: {shape}, Margin: {margin}, Foci: {foci}"
            )
            nodule_lines.append(
                f"    TI-RADS points: {points}"
            )

        recommendation = workflow_result.findings[0].get("recommendation", "") if workflow_result.findings else ""

        sr_text_lines = [
            f"Thyroid TI-RADS Assessment Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Maximum nodule size: {nodule_size_mm:.1f} mm",
            f"  TI-RADS points: {tirads_score}",
            f"  TI-RADS category: {tirads_category}",
            f"  FNA recommended: {'Yes' if fna_recommended else 'No'}",
            f"",
        ]

        if nodule_lines:
            sr_text_lines.extend([
                f"NODULE DETAILS:",
                *nodule_lines,
                f"",
            ])

        sr_text_lines.extend([
            f"SEVERITY: {severity.upper()}",
            f"",
            f"RECOMMENDATION:",
            f"  {recommendation}",
            f"",
            f"ACR TI-RADS Classification Applied:",
            f"  TR1 (0 pts): Benign — no FNA",
            f"  TR2 (2 pts): Not suspicious — no FNA",
            f"  TR3 (3 pts): Mildly suspicious — FNA >= 25mm",
            f"  TR4 (4-6 pts): Moderately suspicious — FNA >= 15mm",
            f"  TR5 (>= 7 pts): Highly suspicious — FNA >= 10mm",
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
            f"Thyroid TI-RADS - {tirads_category}"
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

        # TI-RADS score NUM content
        tirads_content = pydicom.Dataset()
        tirads_content.RelationshipType = "CONTAINS"
        tirads_content.ValueType = "NUM"

        tirads_concept = pydicom.Dataset()
        tirads_concept.CodeValue = "246459004"
        tirads_concept.CodingSchemeDesignator = "SCT"
        tirads_concept.CodeMeaning = "TI-RADS score"
        tirads_content.ConceptNameCodeSequence = [tirads_concept]

        tirads_measured = pydicom.Dataset()
        tirads_measured.NumericValue = f"{tirads_score}"

        tirads_unit = pydicom.Dataset()
        tirads_unit.CodeValue = "1"
        tirads_unit.CodingSchemeDesignator = "UCUM"
        tirads_unit.CodeMeaning = "score"
        tirads_measured.MeasurementUnitsCodeSequence = [tirads_unit]
        tirads_content.MeasuredValueSequence = [tirads_measured]

        # Build content sequence
        sr_dataset.ContentSequence = [
            text_content,
            severity_content,
            tirads_content,
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
            f"severity={severity}, {tirads_category}, "
            f"FNA={'yes' if fna_recommended else 'no'}"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for thyroid TI-RADS.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> ThyroidInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"Thyroid TI-RADS MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run Thyroid TI-RADS MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Thyroid TI-RADS Classification MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM thyroid ultrasound images",
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

    app = ThyroidTIRADSApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"Thyroid TI-RADS Classification complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
