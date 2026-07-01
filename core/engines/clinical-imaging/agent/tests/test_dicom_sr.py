"""Tests for DICOM Structured Report (SR) export module.

Validates DICOMSRExporter generates valid TID 1500 Measurement Reports
from WorkflowResult objects across all clinical workflows.

Author: Adam Jones
Date: April 2026
"""

import json
import os
import tempfile
from pathlib import Path

import pydicom
import pytest

from src.export_dicom_sr import (
    HIGHDICOM_AVAILABLE,
    DICOMSRConfig,
    DICOMSRExporter,
    DICOMSRResult,
    _infer_unit_from_key,
    _UCUM_UNIT_MAP,
)
from src.export import export_dicom_sr
from src.models import FindingSeverity, WorkflowResult, WorkflowStatus


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def exporter():
    """Default DICOMSRExporter instance."""
    return DICOMSRExporter()


@pytest.fixture
def custom_config():
    """Custom DICOMSRConfig."""
    return DICOMSRConfig(
        institution_name="Test Hospital",
        device_name="Test Engine v1.0",
        manufacturer="Test Corp",
        software_version="1.0.0",
    )


@pytest.fixture
def ct_head_result():
    """WorkflowResult for CT head hemorrhage detection."""
    return WorkflowResult(
        workflow_name="ct_head_hemorrhage",
        status=WorkflowStatus.COMPLETED,
        findings=[
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
        measurements={
            "volume_ml": 12.5,
            "midline_shift_mm": 3.2,
            "max_thickness_mm": 8.1,
            "hounsfield_mean": 62.0,
            "hounsfield_max": 78.0,
            "surrounding_edema_ml": 4.3,
        },
        classification="urgent_hemorrhage",
        severity=FindingSeverity.URGENT,
        inference_time_ms=180.0,
        nim_services_used=[
            "SegResNet (MONAI wholeBody_ct_segmentation)",
            "VISTA-3D (optional)",
        ],
        is_mock=True,
    )


@pytest.fixture
def cxr_result():
    """WorkflowResult for CXR rapid findings."""
    return WorkflowResult(
        workflow_name="cxr_rapid_findings",
        status=WorkflowStatus.COMPLETED,
        findings=[
            {
                "category": "consolidation",
                "description": (
                    "Pulmonary consolidation identified, suggestive of "
                    "pneumonia or atelectasis."
                ),
                "severity": "urgent",
                "class_name": "consolidation",
                "confidence": 0.87,
            },
            {
                "category": "effusion",
                "description": "Pleural effusion detected.",
                "severity": "significant",
                "class_name": "pleural_effusion",
                "confidence": 0.72,
            },
        ],
        measurements={
            "consolidation_confidence": 0.87,
            "pleural_effusion_confidence": 0.72,
            "pneumothorax_confidence": 0.08,
        },
        classification="positive: consolidation, pleural_effusion",
        severity=FindingSeverity.URGENT,
        inference_time_ms=25.0,
        nim_services_used=["DenseNet-121 (torchxrayvision)"],
        is_mock=True,
    )


@pytest.fixture
def lung_nodule_result():
    """WorkflowResult for lung nodule with Lung-RADS classification."""
    return WorkflowResult(
        workflow_name="ct_chest_lung_nodule",
        status=WorkflowStatus.COMPLETED,
        findings=[
            {
                "category": "nodule",
                "description": (
                    "Solid pulmonary nodule in right upper lobe, "
                    "diameter 9.2 mm. Lung-RADS 4A."
                ),
                "severity": "significant",
                "location": "right upper lobe",
            },
        ],
        measurements={
            "nodule_diameter_mm": 9.2,
            "nodule_volume_mm3": 407.0,
            "mean_hu": 45.0,
        },
        classification="Lung-RADS 4A",
        severity=FindingSeverity.SIGNIFICANT,
        inference_time_ms=120.0,
        nim_services_used=["RetinaNet (MONAI lung_nodule_ct_detection)"],
        is_mock=True,
    )


@pytest.fixture
def birads_result():
    """WorkflowResult for breast BI-RADS classification."""
    return WorkflowResult(
        workflow_name="breast_birads",
        status=WorkflowStatus.COMPLETED,
        findings=[
            {
                "category": "mass",
                "description": (
                    "Irregular spiculated mass in left breast upper outer "
                    "quadrant, 14 mm. BI-RADS 4C."
                ),
                "severity": "urgent",
                "location": "left breast upper outer quadrant",
            },
        ],
        measurements={
            "mass_diameter_mm": 14.0,
            "mass_volume_ml": 1.4,
        },
        classification="BI-RADS 4C",
        severity=FindingSeverity.URGENT,
        inference_time_ms=95.0,
        nim_services_used=["MONAI segmentation (breast)"],
        is_mock=True,
    )


@pytest.fixture
def output_dir(tmp_path):
    """Temporary output directory for DICOM SR files."""
    return tmp_path


# ===================================================================
# TESTS: MODULE AVAILABILITY
# ===================================================================


def test_exporter_available():
    """HIGHDICOM_AVAILABLE is True since highdicom is in requirements."""
    assert HIGHDICOM_AVAILABLE is True


def test_exporter_available_flag(exporter):
    """Exporter instance reports availability."""
    assert exporter.available is True


# ===================================================================
# TESTS: CONFIG AND RESULT MODELS
# ===================================================================


def test_config_defaults():
    """DICOMSRConfig has correct default values."""
    config = DICOMSRConfig()
    assert config.institution_name == "HCLS AI Factory"
    assert config.device_name == "Clinical Imaging Engine v2.0"
    assert config.manufacturer == "HCLS AI Factory (Apache 2.0)"
    assert config.software_version == "2.0.0"
    assert config.include_measurements is True
    assert config.include_findings is True
    assert config.include_classification is True


def test_config_custom(custom_config):
    """DICOMSRConfig accepts custom values."""
    assert custom_config.institution_name == "Test Hospital"
    assert custom_config.manufacturer == "Test Corp"


def test_result_model():
    """DICOMSRResult serializes correctly to dict and JSON."""
    result = DICOMSRResult(
        success=True,
        output_path="/tmp/test.dcm",
        sop_instance_uid="1.2.3.4.5",
        content_items_count=10,
    )
    d = result.model_dump()
    assert d["success"] is True
    assert d["output_path"] == "/tmp/test.dcm"
    assert d["sop_instance_uid"] == "1.2.3.4.5"
    assert d["content_items_count"] == 10
    assert d["error"] is None

    # JSON round-trip
    j = result.model_dump_json()
    assert "success" in j
    assert "true" in j


def test_result_model_failure():
    """DICOMSRResult captures errors correctly."""
    result = DICOMSRResult(
        success=False,
        error="highdicom not installed",
    )
    assert result.success is False
    assert "highdicom" in result.error
    assert result.output_path is None
    assert result.sop_instance_uid is None


# ===================================================================
# TESTS: MOCK EXPORTS
# ===================================================================


def test_export_mock_ct_head(exporter, output_dir):
    """Mock SR for CT head hemorrhage workflow."""
    path = str(output_dir / "ct_head_sr.dcm")
    result = exporter.export_mock("ct_head_hemorrhage", output_path=path)

    assert result.success is True
    assert result.output_path == path
    assert result.sop_instance_uid is not None
    assert result.content_items_count > 0
    assert Path(path).exists()


def test_export_mock_cxr(exporter, output_dir):
    """Mock SR for CXR rapid findings workflow."""
    path = str(output_dir / "cxr_sr.dcm")
    result = exporter.export_mock("cxr_rapid_findings", output_path=path)

    assert result.success is True
    assert Path(path).exists()
    assert result.content_items_count > 0


def test_export_mock_lung_nodule(exporter, output_dir):
    """Mock SR for lung nodule with Lung-RADS classification."""
    path = str(output_dir / "lung_nodule_sr.dcm")
    result = exporter.export_mock("ct_chest_lung_nodule", output_path=path)

    assert result.success is True
    assert Path(path).exists()


def test_export_mock_birads(exporter, output_dir):
    """Mock SR for breast BI-RADS classification."""
    path = str(output_dir / "breast_sr.dcm")
    result = exporter.export_mock("breast_birads", output_path=path)

    assert result.success is True
    assert Path(path).exists()


def test_export_mock_unknown_workflow(exporter, output_dir):
    """Mock SR for unknown workflow falls back to generic."""
    path = str(output_dir / "unknown_sr.dcm")
    result = exporter.export_mock("unknown_workflow", output_path=path)

    assert result.success is True
    assert Path(path).exists()


# ===================================================================
# TESTS: WORKFLOW RESULT EXPORTS
# ===================================================================


def test_export_with_measurements(exporter, ct_head_result, output_dir):
    """NUM content items are created for measurements."""
    path = str(output_dir / "measurements_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
        patient_id="TEST001",
    )

    assert result.success is True
    assert result.content_items_count > 0

    # Verify the file contains measurement data by reading it
    ds = pydicom.dcmread(path)
    sr_text = str(ds)
    # Should contain measurement-related content
    assert "Volume" in sr_text or "Midline" in sr_text or "Thickness" in sr_text


def test_export_with_classification(exporter, lung_nodule_result, output_dir):
    """CODE content items for scoring system classification."""
    path = str(output_dir / "classification_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=lung_nodule_result,
        output_path=path,
    )

    assert result.success is True

    ds = pydicom.dcmread(path)
    sr_text = str(ds)
    # Should reference the classification
    assert "Lung-RADS" in sr_text or "4A" in sr_text


def test_export_with_findings(exporter, cxr_result, output_dir):
    """TEXT content items for findings."""
    path = str(output_dir / "findings_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=cxr_result,
        output_path=path,
    )

    assert result.success is True

    ds = pydicom.dcmread(path)
    sr_text = str(ds)
    # Should contain finding descriptions
    assert "consolidation" in sr_text.lower() or "effusion" in sr_text.lower()


# ===================================================================
# TESTS: UCUM UNIT MAPPING
# ===================================================================


def test_ucum_unit_mapping():
    """All standard units map correctly."""
    assert _UCUM_UNIT_MAP["mm"] == ("mm", "UCUM", "Millimeter")
    assert _UCUM_UNIT_MAP["cm"] == ("cm", "UCUM", "Centimeter")
    assert _UCUM_UNIT_MAP["mL"] == ("mL", "UCUM", "Milliliter")
    assert _UCUM_UNIT_MAP["ml"] == ("mL", "UCUM", "Milliliter")
    assert _UCUM_UNIT_MAP["mm3"] == ("mm3", "UCUM", "CubicMillimeter")
    assert _UCUM_UNIT_MAP["HU"] == ("[hnsf'U]", "UCUM", "Hounsfield Unit")


def test_ucum_unit_mapping_via_exporter(exporter):
    """Exporter public method maps units correctly."""
    assert exporter.get_ucum_unit("mm") is not None
    assert exporter.get_ucum_unit("cm") is not None
    assert exporter.get_ucum_unit("mL") is not None
    assert exporter.get_ucum_unit("mm3") is not None
    assert exporter.get_ucum_unit("HU") is not None
    assert exporter.get_ucum_unit("unknown_unit") is None


def test_infer_unit_from_key():
    """Unit inference from measurement key suffixes."""
    assert _infer_unit_from_key("volume_ml") == "mL"
    assert _infer_unit_from_key("midline_shift_mm") == "mm"
    assert _infer_unit_from_key("diameter_cm") == "cm"
    assert _infer_unit_from_key("hounsfield_mean") == "HU"
    assert _infer_unit_from_key("hounsfield_max") == "HU"
    assert _infer_unit_from_key("confidence") == ""


# ===================================================================
# TESTS: OUTPUT VALIDATION
# ===================================================================


def test_output_is_valid_dicom(exporter, ct_head_result, output_dir):
    """Output file is readable by pydicom as valid DICOM."""
    path = str(output_dir / "valid_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
    )

    assert result.success is True

    # Must be loadable by pydicom without errors
    ds = pydicom.dcmread(path)
    assert ds is not None

    # Must have SOP Class UID for Comprehensive 3D SR
    assert hasattr(ds, "SOPClassUID")
    # Comprehensive 3D SR: 1.2.840.10008.5.1.4.1.1.88.34
    assert str(ds.SOPClassUID) == "1.2.840.10008.5.1.4.1.1.88.34"

    # Must have a ContentSequence
    assert hasattr(ds, "ContentSequence")
    assert len(ds.ContentSequence) > 0


def test_sr_has_patient_info(exporter, ct_head_result, output_dir):
    """Patient ID, name, and study date are present in SR."""
    path = str(output_dir / "patient_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
        patient_id="PAT12345",
        patient_name="DOE^JOHN",
        study_date="20260413",
    )

    assert result.success is True

    ds = pydicom.dcmread(path)
    assert str(ds.PatientID) == "PAT12345"
    assert str(ds.PatientName) == "DOE^JOHN"
    assert str(ds.StudyDate) == "20260413"


def test_sr_has_device_info(exporter, ct_head_result, output_dir):
    """Manufacturer and software version are present in SR."""
    path = str(output_dir / "device_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
    )

    assert result.success is True

    ds = pydicom.dcmread(path)
    assert "HCLS AI Factory" in str(ds.Manufacturer)
    assert str(ds.InstitutionName) == "HCLS AI Factory"
    assert "2.0.0" in str(ds.SoftwareVersions)


def test_sr_has_custom_device_info(custom_config, ct_head_result, output_dir):
    """Custom config values appear in SR."""
    exp = DICOMSRExporter(config=custom_config)
    path = str(output_dir / "custom_sr.dcm")
    result = exp.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
    )

    assert result.success is True

    ds = pydicom.dcmread(path)
    assert "Test Corp" in str(ds.Manufacturer)
    assert str(ds.InstitutionName) == "Test Hospital"


def test_sr_study_description(exporter, ct_head_result, output_dir):
    """StudyDescription includes workflow name."""
    path = str(output_dir / "study_desc_sr.dcm")
    exporter.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
    )

    ds = pydicom.dcmread(path)
    assert "ct_head_hemorrhage" in str(ds.StudyDescription)


def test_sr_accession_number(exporter, ct_head_result, output_dir):
    """Accession number is stored in SR."""
    path = str(output_dir / "accession_sr.dcm")
    exporter.export_workflow_result(
        workflow_result=ct_head_result,
        output_path=path,
        accession_number="ACC-2026-001",
    )

    ds = pydicom.dcmread(path)
    assert str(ds.AccessionNumber) == "ACC-2026-001"


# ===================================================================
# TESTS: EXPORT FUNCTION IN EXPORT MODULE
# ===================================================================


def test_export_function_in_export_module(ct_head_result, output_dir):
    """export_dicom_sr function in export.py delegates correctly."""
    path = str(output_dir / "export_func_sr.dcm")
    result = export_dicom_sr(
        workflow_result=ct_head_result,
        output_path=path,
        patient_id="FUNC_TEST",
    )

    assert isinstance(result, dict)
    assert result["success"] is True
    assert result["output_path"] == path
    assert result["sop_instance_uid"] is not None
    assert result["content_items_count"] > 0

    # Verify the file is valid DICOM
    ds = pydicom.dcmread(path)
    assert str(ds.PatientID) == "FUNC_TEST"


def test_export_function_returns_dict(ct_head_result, output_dir):
    """export_dicom_sr returns a dict matching DICOMSRResult fields."""
    path = str(output_dir / "dict_sr.dcm")
    result = export_dicom_sr(
        workflow_result=ct_head_result,
        output_path=path,
    )

    expected_keys = {"success", "output_path", "sop_instance_uid", "error", "content_items_count"}
    assert set(result.keys()) == expected_keys


# ===================================================================
# TESTS: EDGE CASES
# ===================================================================


def test_export_empty_workflow(exporter, output_dir):
    """Export a workflow with no findings or measurements."""
    empty_result = WorkflowResult(
        workflow_name="empty_test",
        status=WorkflowStatus.COMPLETED,
        severity=FindingSeverity.NORMAL,
    )
    path = str(output_dir / "empty_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=empty_result,
        output_path=path,
    )

    assert result.success is True
    assert Path(path).exists()

    ds = pydicom.dcmread(path)
    assert hasattr(ds, "ContentSequence")


def test_export_failed_workflow(exporter, output_dir):
    """Export a failed workflow still produces valid SR."""
    failed_result = WorkflowResult(
        workflow_name="failed_test",
        status=WorkflowStatus.FAILED,
        severity=FindingSeverity.ROUTINE,
    )
    path = str(output_dir / "failed_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=failed_result,
        output_path=path,
    )

    assert result.success is True
    assert Path(path).exists()


def test_multiple_findings_single_sr(exporter, cxr_result, output_dir):
    """Multiple findings are all captured in a single SR."""
    path = str(output_dir / "multi_sr.dcm")
    result = exporter.export_workflow_result(
        workflow_result=cxr_result,
        output_path=path,
    )

    assert result.success is True
    # CXR result has 2 findings, each producing multiple content items,
    # plus measurements, classification, severity, and metadata
    assert result.content_items_count >= 5
