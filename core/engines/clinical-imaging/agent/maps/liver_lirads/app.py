"""Liver LI-RADS Classification — MONAI Application Package.

Clinical workflow:
  DICOM CT/MRI Abdomen (multiphase) -> Observation detection ->
  Enhancement pattern analysis (APHE, washout, capsule) ->
  LI-RADS classification -> DICOM SR output

Scoring (LI-RADS v2018):
  - LR-1: Definitely benign
  - LR-2: Probably benign
  - LR-3: Intermediate probability of HCC
  - LR-4: Probably HCC
  - LR-5: Definitely HCC
  - LR-M: Probably or definitely malignant, not HCC specific
  - LR-TIV: Tumor in vein

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
from src.workflows.liver_lirads import LiverLIRADSWorkflow

# ── Optional imports for DICOM SR generation ────────────────────────────
if PYDICOM_AVAILABLE:
    import pydicom
    from pydicom.uid import generate_uid


# ═══════════════════════════════════════════════════════════════════════
# DICOM Series Selection Rules — Liver LI-RADS
# ═══════════════════════════════════════════════════════════════════════

LIVER_LIRADS_RULES = [
    # Primary: CT multiphase abdomen — arterial phase
    SeriesSelectionRule(
        modality="CT",
        body_part="ABDOMEN",
        description_keywords=["arterial", "liver", "multiphase", "dynamic"],
        min_instances=30,
        preferred_slice_thickness_mm=1.5,
    ),
    # Secondary: CT multiphase — portal venous phase
    SeriesSelectionRule(
        modality="CT",
        body_part="ABDOMEN",
        description_keywords=["portal", "venous", "liver"],
        min_instances=30,
        preferred_slice_thickness_mm=1.5,
    ),
    # Tertiary: CT multiphase — delayed phase
    SeriesSelectionRule(
        modality="CT",
        body_part="ABDOMEN",
        description_keywords=["delayed", "equilibrium", "liver"],
        min_instances=20,
        preferred_slice_thickness_mm=2.0,
    ),
    # MRI: liver dynamic contrast-enhanced
    SeriesSelectionRule(
        modality="MR",
        body_part="ABDOMEN",
        description_keywords=["liver", "arterial", "dynamic", "dce", "hepatobiliary"],
        min_instances=15,
        preferred_slice_thickness_mm=3.0,
    ),
    # MRI: liver body part
    SeriesSelectionRule(
        modality="MR",
        body_part="LIVER",
        description_keywords=["arterial", "portal", "delayed", "dynamic"],
        min_instances=10,
        preferred_slice_thickness_mm=3.0,
    ),
    # Fallback: any CT with ABDOMEN body part
    SeriesSelectionRule(
        modality="CT",
        body_part="ABDOMEN",
        description_keywords=["liver"],
        min_instances=10,
    ),
    # Fallback: any CT/MR with LIVER body part
    SeriesSelectionRule(
        modality="CT",
        body_part="LIVER",
        min_instances=10,
    ),
    SeriesSelectionRule(
        modality="MR",
        body_part="LIVER",
        min_instances=5,
    ),
]


class LiverLIRADSApp(BaseImagingMAP):
    """MONAI Deploy Application for liver LI-RADS classification.

    Wraps the existing LiverLIRADSWorkflow in a MONAI Application
    Package with DICOM-native input/output. The clinical logic
    (observation detection, enhancement pattern analysis, LI-RADS
    major and ancillary feature scoring) is fully delegated to the
    workflow class.

    DICOM SR output follows TID 1500 Measurement Report with:
        - Observation size (mm)
        - Arterial phase hyperenhancement (APHE)
        - Washout appearance (portal venous or delayed phase)
        - Enhancing capsule appearance
        - Threshold growth assessment
        - LI-RADS classification (LR-1 through LR-5, LR-M, LR-TIV)
        - Clinical recommendation
    """

    APP_NAME: str = "liver_lirads"
    APP_VERSION: str = "1.0.0"
    APP_DESCRIPTION: str = (
        "GPU-accelerated liver observation classification with LI-RADS v2018. "
        "Analyzes multiphase CT/MRI for arterial phase hyperenhancement, "
        "washout, and capsule appearance to classify hepatocellular carcinoma "
        "probability."
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
        self._workflow: Optional[LiverLIRADSWorkflow] = None

    # ── BaseImagingMAP interface ────────────────────────────────────────

    def get_series_selection_rules(self) -> List[SeriesSelectionRule]:
        """Return liver LI-RADS series selection rules.

        Prefers multiphase CT or MRI with arterial, portal venous, and
        delayed phases. LI-RADS requires multiphasic contrast-enhanced
        imaging for proper characterization of liver observations.
        """
        return LIVER_LIRADS_RULES

    def get_workflow(self) -> LiverLIRADSWorkflow:
        """Return the liver LI-RADS workflow instance.

        Lazily initializes the workflow on first call to avoid loading
        model weights during MAP construction.
        """
        if self._workflow is None:
            self._workflow = LiverLIRADSWorkflow(
                mock_mode=self.mock_mode,
                nim_clients=self._nim_clients,
            )
        return self._workflow

    def build_dicom_sr(
        self,
        workflow_result: WorkflowResult,
        source_dicom: Optional[Any] = None,
    ) -> Optional[bytes]:
        """Build DICOM SR (TID 1500 Measurement Report) from LI-RADS findings.

        Generates a DICOM Structured Report encoding:
            - Observation size
            - APHE (arterial phase hyperenhancement)
            - Washout appearance
            - Enhancing capsule
            - LI-RADS classification
            - Clinical recommendation

        Args:
            workflow_result: WorkflowResult from LiverLIRADSWorkflow.
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
        """Internal SR construction for LI-RADS findings."""
        now = datetime.now(timezone.utc)

        measurements = workflow_result.measurements
        severity = workflow_result.severity.value
        lirads = measurements.get("lirads_category", "LR-1")
        observation_count = int(measurements.get("observation_count", 0))
        max_observation_size_mm = measurements.get("max_observation_size_mm", 0.0)

        # Per-observation details
        obs_lines = []
        for i, finding in enumerate(workflow_result.findings):
            size = finding.get("size_mm", 0.0)
            aphe = finding.get("aphe", "absent")
            washout = finding.get("washout", "absent")
            capsule = finding.get("capsule", "absent")
            lr_cat = finding.get("lirads_category", "LR-1")
            location = finding.get("location", "unspecified")
            obs_lines.append(
                f"  Observation {i + 1}: {size:.1f} mm, {lr_cat}, "
                f"segment {location}"
            )
            obs_lines.append(
                f"    APHE: {aphe}, Washout: {washout}, Capsule: {capsule}"
            )

        recommendation = workflow_result.findings[0].get("recommendation", "") if workflow_result.findings else ""

        sr_text_lines = [
            f"Liver LI-RADS Classification Report",
            f"Generated: {now.strftime('%Y-%m-%d %H:%M:%S UTC')}",
            f"Application: {self.APP_NAME} v{self.APP_VERSION}",
            f"",
            f"FINDINGS:",
            f"  Observation count: {observation_count}",
            f"  Maximum observation size: {max_observation_size_mm:.1f} mm",
            f"  Overall LI-RADS: {lirads}",
            f"",
        ]

        if obs_lines:
            sr_text_lines.extend([
                f"OBSERVATION DETAILS:",
                *obs_lines,
                f"",
            ])

        sr_text_lines.extend([
            f"SEVERITY: {severity.upper()}",
            f"",
            f"RECOMMENDATION:",
            f"  {recommendation}",
            f"",
            f"LI-RADS v2018 Classification Applied:",
            f"  LR-1: Definitely benign",
            f"  LR-2: Probably benign",
            f"  LR-3: Intermediate probability of HCC",
            f"  LR-4: Probably HCC",
            f"  LR-5: Definitely HCC",
            f"  LR-M: Probably/definitely malignant, not HCC specific",
            f"  LR-TIV: Tumor in vein",
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
            f"Liver LI-RADS - {lirads}"
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

        # Observation size NUM content
        size_content = pydicom.Dataset()
        size_content.RelationshipType = "CONTAINS"
        size_content.ValueType = "NUM"

        size_concept = pydicom.Dataset()
        size_concept.CodeValue = "246120007"
        size_concept.CodingSchemeDesignator = "SCT"
        size_concept.CodeMeaning = "Observation size"
        size_content.ConceptNameCodeSequence = [size_concept]

        size_measured = pydicom.Dataset()
        size_measured.NumericValue = f"{max_observation_size_mm:.1f}"

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
            f"severity={severity}, {lirads}, "
            f"observations={observation_count}"
        )
        return sr_bytes

    # ── MONAI Deploy compose ────────────────────────────────────────────

    def compose(self) -> None:
        """Build the MONAI Deploy operator graph for liver LI-RADS.

        Operator chain:
            DICOMDataLoader -> DICOMSeriesSelector -> LiverInference -> SRWriter
        """
        super().compose()
        logger.info(
            f"Liver LI-RADS MAP operator graph composed "
            f"(mock={self.mock_mode})"
        )


# ═══════════════════════════════════════════════════════════════════════
# Standalone entry point
# ═══════════════════════════════════════════════════════════════════════


def main() -> None:
    """Run Liver LI-RADS MAP in standalone mode."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Liver LI-RADS Classification MAP"
    )
    parser.add_argument(
        "--input", "-i",
        default="/input",
        help="Input directory containing DICOM multiphase CT/MRI liver series",
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

    app = LiverLIRADSApp(mock_mode=mock_mode)
    result = app.run_standalone(args.input, args.output)

    workflow_result = result["workflow_result"]
    logger.info(
        f"Liver LI-RADS Classification complete: "
        f"severity={workflow_result.severity.value}, "
        f"classification={workflow_result.classification}, "
        f"elapsed={result['elapsed_sec']:.1f}s"
    )


if __name__ == "__main__":
    main()
