"""CXR Rapid Findings — MONAI Application Package.

Clinical workflow:
  DICOM CXR -> DenseNet-121 multi-label classification ->
  GradCAM heatmap -> Findings with confidence -> DICOM SR output

Detects: consolidation, effusion, pneumothorax, cardiomegaly,
  atelectasis, edema, nodule (18 pathology labels via torchxrayvision)

Target latency: < 30 seconds
NVIDIA technologies: torchxrayvision DenseNet-121, MONAI Deploy SDK
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
from src.workflows.cxr_rapid_findings import (
    CXRRapidFindingsWorkflow,
    CXR_CLASS_NAMES,
    CXR_CLASS_THRESHOLDS,
    CXR_CLASS_DESCRIPTIONS,
)

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — CXR
# ═══════════════════════════════════════════════════════════════════════

CXR_RULES = [
    # Primary: computed radiography / digital radiography chest
    SeriesSelectionRule(
        modality="CR",
        body_part="CHEST",
        description_keywords=["chest", "pa", "ap", "portable"],
        min_instances=1,
    ),
    # Digital X-ray (DX modality)
    SeriesSelectionRule(
        modality="DX",
        body_part="CHEST",
        description_keywords=["chest", "pa", "ap", "portable"],
        min_instances=1,
    ),
    # Fallback: any CR with CHEST body part
    SeriesSelectionRule(
        modality="CR",
        body_part="CHEST",
        min_instances=1,
    ),
    # Fallback: any DX with CHEST body part
    SeriesSelectionRule(
        modality="DX",
        body_part="CHEST",
        min_instances=1,
    ),
]


class CXRRapidFindingsApp(BaseImagingMAP):
    """MONAI Deploy Application for CXR rapid findings triage.

    Wraps the existing CXRRapidFindingsWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (DenseNet-121 multi-label classification, per-class thresholds,
    severity assessment) is fully delegated to the workflow class.

    DICOM SR output follows TID 1500 with:
        - Per-class confidence scores
        - Positive/negative determination per threshold
        - Overall severity classification
        - Clinical recommendations for positive findings
    """

    APP_NAME: str = "cxr_rapid_findings"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "Multi-label chest X-ray classification for emergency triage. "
        "Detects pneumothorax, consolidation, pleural effusion, "
        "cardiomegaly, and fracture with per-class confidence thresholds."
    )
    TARGET_LATENCY_SEC: float = 30.0

    def __init__(
        self,
        mock_mode: bool = True,
        nim_clients: Optional[Dict] = None,
        checkpoint_path: Optional[str] = None,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        super().__init__(mock_mode=mock_mode, *args, **kwargs)
        self._nim_clients = nim_clients
        self._checkpoint_path = checkpoint_path
        self._workflow: Optional[CXRRapidFindingsWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return CXR series selection rules.

        Accepts CR (computed radiography) and DX (digital X-ray) modalities
        with CHEST body part. PA and AP views are both accepted for
        emergency triage, though PA is preferred for cardiothoracic ratio.
        """
        return CXR_RULES

    def get_workflow(self) -> CXRRapidFindingsWorkflow:
        """Return the CXR rapid findings workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = CXRRapidFindingsWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
                checkpoint_path=self._checkpoint_path,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR from CXR multi-label classification findings.

        Generates a DICOM Structured Report encoding:
            - Per-class confidence and threshold status
            - Positive finding descriptions with clinical recommendations
            - Overall severity classification
            - Model metadata (torchxrayvision vs MONAI fallback)

        Args:
            workflow_result: WorkflowResult from CXRRapidFindingsWorkflow.
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
        """Internal SR construction for CXR findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        classification = workflow_result.classification

        # Build per-class results table
        class_lines = []
        for class_name in CXR_CLASS_NAMES:
            confidence = measurements.get(f"{class_name}_confidence", 0.0)
            threshold = CXR_CLASS_THRESHOLDS.get(class_name, 0.50)
            status = "POSITIVE" if confidence >= threshold else "negative"
            class_lines.append(
                f"  {class_name:<20s}  conf={confidence:.3f}  "
                f"thresh={threshold:.2f}  [{status}]"
            )

        # Build positive findings detail
        positive_lines = []
        for finding in workflow_result.findings:
            if finding.get("above_threshold"):
                desc = finding.get("description", "")
                conf = finding.get("confidence", 0.0)
                positive_lines.append(f"  [{conf:.3f}] {desc}")

        sr_text_lines = [
            f"CXR Rapid Findings Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"CLASSIFICATION RESULTS:",
            *class_lines,
            f"",
            f"OVERALL: {classification}",
            f"SEVERITY: {severity.upper()}",
            f"",
        ]

        if positive_lines:
            sr_text_lines.extend([
                f"POSITIVE FINDINGS:",
                *positive_lines,
                f"",
            ])
        else:
            sr_text_lines.extend([
                f"No significant acute findings detected.",
                f"All classification scores below clinical thresholds.",
                f"",
            ])

        using_xrv = measurements.get("using_xrv", 0.0) > 0
        model_backend = (
            "torchxrayvision (densenet121-res224-all)"
            if using_xrv
            else "MONAI DenseNet-121"
        )
        sr_text_lines.extend([
            f"MODEL: {model_backend}",
            f"INPUT: 224x224 normalized",
        ])

        # Create DICOM SR dataset
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
            f"CXR Rapid Findings - {severity.upper()}"
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
        code_value, code_meaning = severity_map.get(
            severity, ("17621005", "Routine")
        )
        severity_code.CodeValue = code_value
        severity_code.CodingSchemeDesignator = "SCT"
        severity_code.CodeMeaning = code_meaning
        severity_content.ConceptCodeSequence = [severity_code]

        # Per-class NUM content items for positive findings
        class_content_items = []
        for class_name in CXR_CLASS_NAMES:
            confidence = measurements.get(f"{class_name}_confidence", 0.0)
            threshold = CXR_CLASS_THRESHOLDS.get(class_name, 0.50)

            if confidence >= threshold:
                num_item = pydicom.Dataset()
                num_item.RelationshipType = "CONTAINS"
                num_item.ValueType = "NUM"

                num_concept = pydicom.Dataset()
                num_concept.CodeValue = "363714003"
                num_concept.CodingSchemeDesignator = "SCT"
                num_concept.CodeMeaning = (
                    f"{class_name.replace('_', ' ').title()} confidence"
                )
                num_item.ConceptNameCodeSequence = [num_concept]

                measured = pydicom.Dataset()
                measured.NumericValue = f"{confidence:.4f}"

                unit = pydicom.Dataset()
                unit.CodeValue = "1"
                unit.CodingSchemeDesignator = "UCUM"
                unit.CodeMeaning = "probability"
                measured.MeasurementUnitsCodeSequence = [unit]
                num_item.MeasuredValueSequence = [measured]

                class_content_items.append(num_item)

        # Build content sequence
        sr_dataset.ContentSequence = [
            text_content,
            severity_content,
            *class_content_items,
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
        file_meta.TransferSyntaxUID = "1.2.840.10008.1.2.1"  # Explicit VR LE
        file_meta.ImplementationClassUID = "1.2.826.0.1.3680043.8.498.1"
        sr_dataset.file_meta = file_meta
        sr_dataset.preamble = b"\x00" * 128

        # Serialize
        from io import BytesIO

        buffer = BytesIO()
        pydicom.dcmwrite(buffer, sr_dataset)
        sr_bytes = buffer.getvalue()

        logger.info(
            f"Generated CXR DICOM SR: {len(sr_bytes)} bytes, "
            f"severity={severity}, classification={classification}"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for CXR rapid findings.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> CXRInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"CXR Rapid Findings MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run CXR Rapid Findings MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="CXR Rapid Findings MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM CXR images",
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
        help="Run with real model inference (requires GPU + torchxrayvision)",
    )
    parser.add_argument(
        "--checkpoint",
        default=None,
        help="Path to fine-tuned DenseNet-121 checkpoint (MONAI fallback path)",
    )
    args = parser.parse_args()

    mock_mode = not args.real

    app = CXRRapidFindingsApp(
        mock_mode=mock_mode,
        checkpoint_path=args.checkpoint,
    )
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"CXR Rapid Findings complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
