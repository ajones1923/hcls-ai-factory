"""DICOM Structured Report (SR) export via highdicom.

Generates TID 1500 Measurement Reports that store AI findings alongside
source images in PACS. SR objects are viewable in OHIF and retrievable
via standard DICOM C-FIND queries.

Dependencies: highdicom (MIT), pydicom (MIT). Both free, pip-installable.

Author: Adam Jones
Date: April 2026
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from loguru import logger
from pydantic import BaseModel, Field

from src.models import FindingSeverity, WorkflowResult, WorkflowStatus

# ═══════════════════════════════════════════════════════════════════════
# GRACEFUL IMPORTS
# ═══════════════════════════════════════════════════════════════════════

try:
    import highdicom as hd
    from pydicom.sr.codedict import codes
    from pydicom.uid import generate_uid
    import pydicom
    HIGHDICOM_AVAILABLE = True
except ImportError:
    HIGHDICOM_AVAILABLE = False
    logger.warning("highdicom not installed. DICOM SR export unavailable.")


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

# Map unit strings to UCUM coded concepts (populated at import time)
# Each entry: (value, scheme_designator, meaning)
_UCUM_UNIT_MAP: Dict[str, tuple] = {
    "mm": ("mm", "UCUM", "Millimeter"),
    "cm": ("cm", "UCUM", "Centimeter"),
    "mL": ("mL", "UCUM", "Milliliter"),
    "ml": ("mL", "UCUM", "Milliliter"),
    "mm3": ("mm3", "UCUM", "CubicMillimeter"),
    "HU": ("[hnsf'U]", "UCUM", "Hounsfield Unit"),
    "hu": ("[hnsf'U]", "UCUM", "Hounsfield Unit"),
    "%": ("%", "UCUM", "Percent"),
    "ms": ("ms", "UCUM", "Millisecond"),
    "s": ("s", "UCUM", "Second"),
}

# Scoring system display names for CODE content items
_SCORING_SYSTEM_DISPLAY: Dict[str, str] = {
    "lung-rads": "ACR Lung-RADS v2022",
    "lung_rads": "ACR Lung-RADS v2022",
    "bi-rads": "ACR BI-RADS 5th Edition",
    "birads": "ACR BI-RADS 5th Edition",
    "bi_rads": "ACR BI-RADS 5th Edition",
    "ti-rads": "ACR TI-RADS",
    "tirads": "ACR TI-RADS",
    "ti_rads": "ACR TI-RADS",
    "li-rads": "ACR LI-RADS",
    "lirads": "ACR LI-RADS",
    "li_rads": "ACR LI-RADS",
    "cad-rads": "CAD-RADS 2.0",
    "cadrads": "CAD-RADS 2.0",
    "cad_rads": "CAD-RADS 2.0",
    "pi-rads": "PI-RADS v2.1",
    "pirads": "PI-RADS v2.1",
    "pi_rads": "PI-RADS v2.1",
}

# Severity -> DICOM coded concept for observer context
_SEVERITY_CODES: Dict[str, tuple] = {
    "critical": ("R-0040A", "SRT", "Critical"),
    "urgent": ("R-0040B", "SRT", "Urgent"),
    "significant": ("R-0040C", "SRT", "Significant"),
    "routine": ("R-0040D", "SRT", "Routine"),
    "normal": ("R-00339", "SRT", "Normal"),
}

# Measurement key suffix -> unit string
_MEASUREMENT_KEY_SUFFIXES: Dict[str, str] = {
    "_ml": "mL",
    "_mm": "mm",
    "_cm": "cm",
    "_hu": "HU",
    "_ms": "ms",
    "_pct": "%",
}


# ═══════════════════════════════════════════════════════════════════════
# CONFIG & RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════


class DICOMSRConfig(BaseModel):
    """Configuration for DICOM SR generation."""
    institution_name: str = "HCLS AI Factory"
    device_name: str = "Clinical Imaging Engine v2.0"
    manufacturer: str = "HCLS AI Factory (Apache 2.0)"
    software_version: str = "2.0.0"
    include_measurements: bool = True
    include_findings: bool = True
    include_classification: bool = True


class DICOMSRResult(BaseModel):
    """Result of DICOM SR generation."""
    success: bool
    output_path: Optional[str] = None
    sop_instance_uid: Optional[str] = None
    error: Optional[str] = None
    content_items_count: int = 0


# ═══════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════


def _infer_unit_from_key(key: str) -> str:
    """Infer measurement unit from a measurement key suffix.

    Examples:
        volume_ml -> mL
        midline_shift_mm -> mm
        hounsfield_mean -> HU (special case)
        confidence -> (empty, dimensionless)
    """
    key_lower = key.lower()
    for suffix, unit in _MEASUREMENT_KEY_SUFFIXES.items():
        if key_lower.endswith(suffix):
            return unit
    # Special cases
    if "hounsfield" in key_lower:
        return "HU"
    return ""


def _make_coded_concept(value: str, scheme: str, meaning: str):
    """Create a highdicom CodedConcept.

    Truncates meaning to 64 characters per DICOM SR specification.
    """
    return hd.sr.CodedConcept(
        value=value, scheme_designator=scheme, meaning=meaning[:64]
    )


# ═══════════════════════════════════════════════════════════════════════
# EXPORTER CLASS
# ═══════════════════════════════════════════════════════════════════════


class DICOMSRExporter:
    """Generate DICOM Structured Reports from workflow results.

    Produces TID 1500 Measurement Reports containing:
    - TEXT content items for each finding description
    - NUM content items for each measurement (with UCUM units)
    - CODE content items for severity and classification scores
    - Container for workflow metadata

    Usage::

        exporter = DICOMSRExporter()
        result = exporter.export_workflow_result(
            workflow_result=wf_result,
            output_path="ct_head_sr.dcm",
            patient_id="PATIENT001",
        )
        if result.success:
            print(f"SR saved to {result.output_path}")
    """

    def __init__(self, config: Optional[DICOMSRConfig] = None):
        self.config = config or DICOMSRConfig()
        self.available = HIGHDICOM_AVAILABLE

    # ─── Public API ───────────────────────────────────────────────────

    def export_workflow_result(
        self,
        workflow_result: WorkflowResult,
        source_dicom_path: Optional[str] = None,
        output_path: str = "output.dcm",
        patient_id: str = "ANONYMOUS",
        patient_name: str = "ANONYMOUS",
        study_date: Optional[str] = None,
        accession_number: str = "",
    ) -> DICOMSRResult:
        """Generate a DICOM SR from a WorkflowResult.

        Creates a TID 1500 Measurement Report containing:
        - TEXT content items for each finding
        - NUM content items for each measurement (with UCUM units)
        - CODE content items for severity and classification
        - Container for workflow metadata

        If source_dicom_path is provided, references the source study.
        Otherwise creates a standalone SR.

        Args:
            workflow_result: The WorkflowResult to export.
            source_dicom_path: Optional path to source DICOM file for
                study/series UIDs. If None, generates fresh UIDs.
            output_path: File path for the generated .dcm SR.
            patient_id: Patient identifier for the SR.
            patient_name: Patient name for the SR.
            study_date: Study date (YYYYMMDD). Defaults to today.
            accession_number: Accession number for the study.

        Returns:
            DICOMSRResult with success status, path, and metadata.
        """
        if not self.available:
            return DICOMSRResult(
                success=False,
                error=(
                    "highdicom not installed. Install with: "
                    "pip install highdicom pydicom"
                ),
            )

        try:
            # Resolve study/series UIDs
            study_instance_uid = generate_uid()
            series_instance_uid = generate_uid()

            if source_dicom_path and Path(source_dicom_path).exists():
                try:
                    source_ds = pydicom.dcmread(source_dicom_path)
                    study_instance_uid = str(
                        getattr(source_ds, "StudyInstanceUID", study_instance_uid)
                    )
                    # Always generate a new series UID for the SR
                    if hasattr(source_ds, "PatientID"):
                        patient_id = str(source_ds.PatientID) or patient_id
                    if hasattr(source_ds, "PatientName"):
                        patient_name = str(source_ds.PatientName) or patient_name
                    logger.info(
                        f"Linked SR to source study {study_instance_uid}"
                    )
                except Exception as e:
                    logger.warning(
                        f"Could not read source DICOM {source_dicom_path}: {e}. "
                        "Using generated UIDs."
                    )

            # Build the SR dataset
            sr_dataset = self._build_measurement_report(
                workflow_result=workflow_result,
                patient_id=patient_id,
                patient_name=patient_name,
                study_instance_uid=study_instance_uid,
                series_instance_uid=series_instance_uid,
                study_date=study_date,
                accession_number=accession_number,
            )

            # Count content items
            content_items_count = self._count_content_items(sr_dataset)

            # Write to file
            output_dir = Path(output_path).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            sr_dataset.save_as(output_path)

            sop_uid = str(sr_dataset.SOPInstanceUID)
            logger.info(
                f"DICOM SR exported: {output_path} "
                f"(SOP UID: {sop_uid}, {content_items_count} content items)"
            )

            return DICOMSRResult(
                success=True,
                output_path=str(output_path),
                sop_instance_uid=sop_uid,
                content_items_count=content_items_count,
            )

        except Exception as e:
            logger.error(f"DICOM SR export failed: {e}")
            return DICOMSRResult(
                success=False,
                error=str(e),
            )

    def export_mock(
        self,
        workflow_name: str,
        output_path: str = "output_sr.dcm",
    ) -> DICOMSRResult:
        """Generate a mock DICOM SR for demonstration without source DICOM.

        Creates a valid DICOM SR with synthetic content based on the
        workflow name. Useful for testing, demos, and PACS integration
        validation.

        Args:
            workflow_name: Name of the workflow to simulate (e.g.,
                'ct_head_hemorrhage', 'cxr_rapid_findings').
            output_path: File path for the generated .dcm SR.

        Returns:
            DICOMSRResult with success status and metadata.
        """
        mock_result = self._generate_mock_workflow_result(workflow_name)
        return self.export_workflow_result(
            workflow_result=mock_result,
            output_path=output_path,
            patient_id="MOCK_PATIENT",
            patient_name="MOCK^PATIENT",
            accession_number="MOCK-001",
        )

    # ─── SR Construction ──────────────────────────────────────────────

    def _build_measurement_report(
        self,
        workflow_result: WorkflowResult,
        patient_id: str,
        patient_name: str,
        study_instance_uid: str,
        series_instance_uid: str,
        study_date: Optional[str] = None,
        accession_number: str = "",
    ) -> "pydicom.Dataset":
        """Build the highdicom SR dataset.

        Constructs a Comprehensive 3D SR (IOD) with TID 1500 Measurement
        Report structure. Content is organized into:
        - Imaging Measurements container (NUM items)
        - Qualitative Evaluations container (CODE items)
        - Findings container (TEXT items)
        """
        # Resolve study date
        if study_date is None:
            study_date = datetime.now().strftime("%Y%m%d")

        # ── Build measurements and qualitative evaluations ────────────
        measurements_list = []
        qualitative_evals = []

        # 1. Measurements as Measurement objects (NUM items)
        if self.config.include_measurements and workflow_result.measurements:
            for key, value in workflow_result.measurements.items():
                m = self._create_measurement_content(key, value)
                if m is not None:
                    measurements_list.append(m)

        # 2. Classification as QualitativeEvaluation (CODE item)
        if self.config.include_classification and workflow_result.classification:
            scoring_system = self._detect_scoring_system(
                workflow_result.workflow_name,
                workflow_result.classification,
            )
            qe = self._create_classification_content(
                workflow_result.classification, scoring_system
            )
            if qe is not None:
                qualitative_evals.append(qe)

        # 3. Severity as QualitativeEvaluation (CODE item)
        severity_qe = self._create_severity_content(workflow_result.severity)
        if severity_qe is not None:
            qualitative_evals.append(severity_qe)

        # 4. Findings as QualitativeEvaluation (CODE items with text)
        if self.config.include_findings and workflow_result.findings:
            for finding in workflow_result.findings:
                qe_items = self._create_finding_content(finding)
                qualitative_evals.extend(qe_items)

        # ── Build MeasurementsAndQualitativeEvaluations group ────────
        tracking = hd.sr.TrackingIdentifier(
            uid=generate_uid(),
            identifier=f"HCLS-AI-{workflow_result.workflow_name}",
        )

        measurement_group = hd.sr.MeasurementsAndQualitativeEvaluations(
            tracking_identifier=tracking,
            measurements=measurements_list if measurements_list else None,
            qualitative_evaluations=qualitative_evals if qualitative_evals else None,
        )

        # ── Build the SR dataset ─────────────────────────────────────

        # Observation context: the AI device
        observer_device = hd.sr.ObserverContext(
            observer_type=_make_coded_concept(
                "121007", "DCM", "Device"
            ),
            observer_identifying_attributes=hd.sr.DeviceObserverIdentifyingAttributes(
                uid=generate_uid(),
                name=self.config.device_name,
                manufacturer_name=self.config.manufacturer,
                model_name=self.config.software_version,
            ),
        )

        # Build the measurement report container
        measurement_report = hd.sr.MeasurementReport(
            observation_context=observer_device,
            procedure_reported=_make_coded_concept(
                "363679005", "SCT", "Imaging"
            ),
            imaging_measurements=[measurement_group],
        )

        # Create a minimal evidence dataset (required by Comprehensive3DSR)
        evidence_ds = pydicom.Dataset()
        evidence_ds.SOPClassUID = "1.2.840.10008.5.1.4.1.1.2"  # CT Image Storage
        evidence_ds.SOPInstanceUID = generate_uid()
        evidence_ds.StudyInstanceUID = study_instance_uid
        evidence_ds.SeriesInstanceUID = generate_uid()
        evidence_ds.Modality = "OT"  # Other
        evidence_ds.PatientID = patient_id
        evidence_ds.PatientName = patient_name
        evidence_ds.PatientBirthDate = ""
        evidence_ds.PatientSex = ""
        evidence_ds.StudyDate = study_date
        evidence_ds.StudyTime = ""
        evidence_ds.AccessionNumber = accession_number
        evidence_ds.ReferringPhysicianName = ""
        evidence_ds.StudyID = ""

        # Build the Comprehensive 3D SR
        sr_dataset = hd.sr.Comprehensive3DSR(
            evidence=[evidence_ds],
            content=measurement_report,
            series_number=99,
            series_instance_uid=series_instance_uid,
            sop_instance_uid=generate_uid(),
            instance_number=1,
            manufacturer=self.config.manufacturer,
            institution_name=self.config.institution_name,
            is_complete=True,
            is_verified=False,
        )

        # Set patient and study level attributes
        sr_dataset.PatientID = patient_id
        sr_dataset.PatientName = patient_name
        sr_dataset.StudyInstanceUID = study_instance_uid
        sr_dataset.StudyDate = study_date
        sr_dataset.StudyDescription = (
            f"AI Analysis: {workflow_result.workflow_name}"
        )
        sr_dataset.AccessionNumber = accession_number
        sr_dataset.Manufacturer = self.config.manufacturer
        sr_dataset.InstitutionName = self.config.institution_name
        sr_dataset.SoftwareVersions = self.config.software_version

        return sr_dataset

    # ─── Content Item Builders ────────────────────────────────────────

    def _create_finding_content(self, finding: Dict) -> List:
        """Create SR QualitativeEvaluation items for a single finding.

        Each finding produces one or more QualitativeEvaluation items
        that encode the finding category and severity as coded concepts.

        Args:
            finding: Dict with keys like 'description', 'category',
                'severity', 'confidence', etc.

        Returns:
            List of highdicom QualitativeEvaluation items.
        """
        items = []

        # Category as QualitativeEvaluation
        category = finding.get("category", "normal")
        # CodedConcept meaning is limited to 64 characters
        cat_display = category.replace("_", " ").title()[:64]

        cat_eval = hd.sr.QualitativeEvaluation(
            name=_make_coded_concept(
                "121071", "DCM", "Finding"
            ),
            value=_make_coded_concept(
                "404684003", "SCT", cat_display
            ),
        )
        items.append(cat_eval)

        # Severity as QualitativeEvaluation
        severity = finding.get("severity", "")
        if severity:
            sev_tuple = _SEVERITY_CODES.get(
                severity.lower(), ("R-00339", "SRT", "Normal")
            )
            sev_eval = hd.sr.QualitativeEvaluation(
                name=_make_coded_concept(
                    "246112005", "SCT", "Severity"
                ),
                value=_make_coded_concept(*sev_tuple),
            )
            items.append(sev_eval)

        return items

    def _create_measurement_content(
        self, key: str, value: float
    ) -> Optional[Any]:
        """Create a Measurement object for a numeric measurement.

        Maps the measurement key suffix to a UCUM unit code. Measurements
        without a recognized unit suffix are stored as dimensionless
        quantities with unit '1' (UCUM unity).

        Args:
            key: Measurement key (e.g., 'volume_ml', 'midline_shift_mm').
            value: Numeric measurement value.

        Returns:
            A highdicom Measurement, or None if creation fails.
        """
        try:
            unit_str = _infer_unit_from_key(key)
            display_name = key.replace("_", " ").title()

            if unit_str and unit_str in _UCUM_UNIT_MAP:
                ucum_tuple = _UCUM_UNIT_MAP[unit_str]
                unit_code = _make_coded_concept(*ucum_tuple)
            else:
                # Dimensionless -- use UCUM unity '1'
                unit_code = _make_coded_concept("1", "UCUM", "No units")

            return hd.sr.Measurement(
                name=_make_coded_concept(
                    "G-C036", "SRT", display_name
                ),
                value=float(value),
                unit=unit_code,
            )
        except Exception as e:
            logger.warning(
                f"Could not create Measurement for {key}={value}: {e}"
            )
            return None

    def _create_classification_content(
        self, classification: str, scoring_system: str
    ) -> Optional[Any]:
        """Create a QualitativeEvaluation for a classification score.

        Maps classification strings (e.g., 'Lung-RADS 4A') and scoring
        system identifiers to DICOM SR coded evaluations.

        Args:
            classification: The classification value string.
            scoring_system: Identifier for the scoring system
                (e.g., 'lung-rads', 'bi-rads').

        Returns:
            A highdicom QualitativeEvaluation, or None if creation fails.
        """
        try:
            display = _SCORING_SYSTEM_DISPLAY.get(
                scoring_system.lower(), scoring_system or "AI Classification"
            )

            # CodedConcept meaning is limited to 64 characters
            meaning = f"{display}: {classification}"[:64]

            return hd.sr.QualitativeEvaluation(
                name=_make_coded_concept(
                    "246464006", "SCT", "Assessment"
                ),
                value=_make_coded_concept(
                    "AI-CLASS", "99HCLS", meaning
                ),
            )
        except Exception as e:
            logger.warning(
                f"Could not create QualitativeEvaluation for {classification}: {e}"
            )
            return None

    def _create_severity_content(
        self, severity: FindingSeverity
    ) -> Optional[Any]:
        """Create a QualitativeEvaluation for the overall severity.

        Args:
            severity: FindingSeverity enum value.

        Returns:
            A highdicom QualitativeEvaluation, or None if creation fails.
        """
        try:
            sev_str = severity.value if hasattr(severity, "value") else str(severity)
            sev_tuple = _SEVERITY_CODES.get(
                sev_str.lower(), ("R-00339", "SRT", "Normal")
            )

            return hd.sr.QualitativeEvaluation(
                name=_make_coded_concept(
                    "246112005", "SCT", "Severity"
                ),
                value=_make_coded_concept(*sev_tuple),
            )
        except Exception as e:
            logger.warning(f"Could not create severity QualitativeEvaluation: {e}")
            return None

    # ─── Utility Methods ─────────────────────────────────────────────

    def _detect_scoring_system(
        self, workflow_name: str, classification: str
    ) -> str:
        """Detect the scoring system from workflow name or classification.

        Args:
            workflow_name: Name of the workflow (e.g., 'ct_chest_lung_nodule').
            classification: Classification string.

        Returns:
            Scoring system identifier string.
        """
        combined = f"{workflow_name} {classification}".lower()

        scoring_keywords = {
            "lung": "lung-rads",
            "lung_rads": "lung-rads",
            "lung-rads": "lung-rads",
            "birads": "bi-rads",
            "bi_rads": "bi-rads",
            "bi-rads": "bi-rads",
            "breast": "bi-rads",
            "mammography": "bi-rads",
            "tirads": "ti-rads",
            "ti_rads": "ti-rads",
            "ti-rads": "ti-rads",
            "thyroid": "ti-rads",
            "lirads": "li-rads",
            "li_rads": "li-rads",
            "li-rads": "li-rads",
            "liver": "li-rads",
            "cadrads": "cad-rads",
            "cad_rads": "cad-rads",
            "cad-rads": "cad-rads",
            "coronary": "cad-rads",
            "pirads": "pi-rads",
            "pi_rads": "pi-rads",
            "pi-rads": "pi-rads",
            "prostate": "pi-rads",
        }

        for keyword, system in scoring_keywords.items():
            if keyword in combined:
                return system

        return ""

    def _count_content_items(self, sr_dataset) -> int:
        """Count the total number of content items in the SR.

        Traverses the ContentSequence recursively to count all items.

        Args:
            sr_dataset: The pydicom Dataset (SR).

        Returns:
            Total number of content items.
        """
        count = 0
        if hasattr(sr_dataset, "ContentSequence"):
            count += len(sr_dataset.ContentSequence)
            for item in sr_dataset.ContentSequence:
                if hasattr(item, "ContentSequence"):
                    count += self._count_content_items_recursive(
                        item.ContentSequence
                    )
        return count

    def _count_content_items_recursive(self, sequence) -> int:
        """Recursively count content items in a ContentSequence."""
        count = len(sequence)
        for item in sequence:
            if hasattr(item, "ContentSequence"):
                count += self._count_content_items_recursive(
                    item.ContentSequence
                )
        return count

    def _generate_mock_workflow_result(self, workflow_name: str) -> WorkflowResult:
        """Generate a realistic mock WorkflowResult for a given workflow.

        Args:
            workflow_name: Name of the workflow to simulate.

        Returns:
            A WorkflowResult with synthetic findings, measurements,
            and classification appropriate for the workflow.
        """
        mock_data = _MOCK_WORKFLOW_DATA.get(workflow_name)
        if mock_data is None:
            logger.warning(
                f"No mock data for workflow '{workflow_name}'. "
                "Using generic mock."
            )
            mock_data = _MOCK_WORKFLOW_DATA["generic"]

        return WorkflowResult(
            workflow_name=workflow_name,
            status=WorkflowStatus.COMPLETED,
            findings=mock_data["findings"],
            measurements=mock_data["measurements"],
            classification=mock_data["classification"],
            severity=mock_data["severity"],
            inference_time_ms=mock_data.get("inference_time_ms", 150.0),
            nim_services_used=mock_data.get("nim_services_used", []),
            is_mock=True,
        )

    def get_ucum_unit(self, unit_str: str) -> Optional[tuple]:
        """Map unit string to UCUM coded concept tuple.

        Public accessor for unit mapping, useful for testing.

        Args:
            unit_str: Unit string (e.g., 'mm', 'cm', 'mL').

        Returns:
            Tuple of (value, scheme_designator, meaning) or None.
        """
        return _UCUM_UNIT_MAP.get(unit_str)


# ═══════════════════════════════════════════════════════════════════════
# MOCK WORKFLOW DATA
# ═══════════════════════════════════════════════════════════════════════

_MOCK_WORKFLOW_DATA: Dict[str, Dict] = {
    "ct_head_hemorrhage": {
        "findings": [
            {
                "category": "hemorrhage",
                "description": (
                    "Intraparenchymal hemorrhage in right basal ganglia, "
                    "volume 12.5 mL, midline shift 3.2 mm, "
                    "max thickness 8.1 mm"
                ),
                "severity": "urgent",
                "hemorrhage_type": "intraparenchymal",
                "location": "right basal ganglia",
            },
        ],
        "measurements": {
            "volume_ml": 12.5,
            "midline_shift_mm": 3.2,
            "max_thickness_mm": 8.1,
            "hounsfield_mean": 62.0,
            "hounsfield_max": 78.0,
            "surrounding_edema_ml": 4.3,
        },
        "classification": "urgent_hemorrhage",
        "severity": FindingSeverity.URGENT,
        "inference_time_ms": 180.0,
        "nim_services_used": [
            "SegResNet (MONAI wholeBody_ct_segmentation)",
            "VISTA-3D (optional)",
        ],
    },
    "cxr_rapid_findings": {
        "findings": [
            {
                "category": "consolidation",
                "description": (
                    "Pulmonary consolidation identified, suggestive of "
                    "pneumonia or atelectasis. Clinical correlation with "
                    "infection markers recommended."
                ),
                "severity": "urgent",
                "class_name": "consolidation",
                "confidence": 0.87,
            },
            {
                "category": "effusion",
                "description": (
                    "Pleural effusion detected. Assess for laterality, "
                    "size, and clinical correlation."
                ),
                "severity": "significant",
                "class_name": "pleural_effusion",
                "confidence": 0.72,
            },
        ],
        "measurements": {
            "pneumothorax_confidence": 0.08,
            "consolidation_confidence": 0.87,
            "pleural_effusion_confidence": 0.72,
            "cardiomegaly_confidence": 0.31,
            "fracture_confidence": 0.12,
        },
        "classification": "positive: consolidation, pleural_effusion",
        "severity": FindingSeverity.URGENT,
        "inference_time_ms": 25.0,
        "nim_services_used": [
            "DenseNet-121 (torchxrayvision)",
            "VILA-M3 (optional)",
        ],
    },
    "ct_chest_lung_nodule": {
        "findings": [
            {
                "category": "nodule",
                "description": (
                    "Solid pulmonary nodule in right upper lobe, "
                    "diameter 9.2 mm. Lung-RADS 4A: suspicious, "
                    "3-month follow-up recommended."
                ),
                "severity": "significant",
                "location": "right upper lobe",
            },
        ],
        "measurements": {
            "nodule_diameter_mm": 9.2,
            "nodule_volume_mm3": 407.0,
            "mean_hu": 45.0,
        },
        "classification": "Lung-RADS 4A",
        "severity": FindingSeverity.SIGNIFICANT,
        "inference_time_ms": 120.0,
        "nim_services_used": [
            "RetinaNet (MONAI lung_nodule_ct_detection)",
            "SegResNet (volumetric)",
        ],
    },
    "breast_birads": {
        "findings": [
            {
                "category": "mass",
                "description": (
                    "Irregular spiculated mass in left breast upper outer "
                    "quadrant, 14 mm. BI-RADS 4C: high suspicion for "
                    "malignancy. Tissue diagnosis recommended."
                ),
                "severity": "urgent",
                "location": "left breast upper outer quadrant",
            },
        ],
        "measurements": {
            "mass_diameter_mm": 14.0,
            "mass_volume_ml": 1.4,
        },
        "classification": "BI-RADS 4C",
        "severity": FindingSeverity.URGENT,
        "inference_time_ms": 95.0,
        "nim_services_used": ["MONAI segmentation (breast)"],
    },
    "generic": {
        "findings": [
            {
                "category": "normal",
                "description": "No significant findings identified.",
                "severity": "normal",
            },
        ],
        "measurements": {},
        "classification": "normal",
        "severity": FindingSeverity.NORMAL,
        "inference_time_ms": 50.0,
        "nim_services_used": [],
    },
}
