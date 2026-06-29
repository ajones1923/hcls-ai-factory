"""Tests for workflow execution patterns.

Verifies WorkflowResult creation and serialization for all workflow types.

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    NeuroWorkflowType,
    SeverityLevel,
    WorkflowResult,
    ScaleResult,
    ClinicalScaleType,
)
from src.export import NeuroReportExporter


class TestWorkflowResultSerialization:
    """Tests for WorkflowResult serialization."""

    def test_model_dump(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.ACUTE_STROKE,
            findings=["LVO detected"],
            recommendations=["Thrombectomy"],
            severity=SeverityLevel.CRITICAL,
        )
        dumped = result.model_dump()
        assert dumped["workflow_type"] == "acute_stroke"
        assert dumped["severity"] == "critical"
        assert len(dumped["findings"]) == 1

    def test_model_dump_with_scales(self):
        scale = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=16,
            max_score=42,
            interpretation="Moderate to severe stroke",
            severity_category="moderate_severe",
        )
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.ACUTE_STROKE,
            findings=["NIHSS 16"],
            scale_results=[scale],
        )
        dumped = result.model_dump()
        assert len(dumped["scale_results"]) == 1
        assert dumped["scale_results"][0]["score"] == 16


class TestWorkflowExportRoundtrip:
    """Tests for workflow result -> export roundtrip."""

    def test_markdown_export(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.DEMENTIA_EVALUATION,
            findings=["MoCA 18/30", "Medial temporal atrophy"],
            recommendations=["Amyloid PET", "CSF biomarkers"],
            severity=SeverityLevel.MODERATE,
        )
        exporter = NeuroReportExporter()
        md = exporter.export_markdown(result)
        assert "Findings" in md
        assert "Recommendations" in md
        assert "Disclaimer" in md

    def test_json_export(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.EPILEPSY_FOCUS,
            findings=["Right temporal spikes"],
            recommendations=["Epilepsy surgery evaluation"],
        )
        exporter = NeuroReportExporter()
        json_data = exporter.export_json(result)
        assert json_data["report_type"] == "neurology_workflow"
        assert "data" in json_data
        assert json_data["version"] == "1.0.0"

    def test_fhir_r4_export(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.MS_MONITORING,
            findings=["New T2 lesions"],
        )
        exporter = NeuroReportExporter()
        bundle = exporter.export_fhir_r4(result, patient_id="NEURO-001")
        assert bundle["resourceType"] == "Bundle"
        assert bundle["type"] == "document"
        assert len(bundle["entry"]) == 1
        report = bundle["entry"][0]["resource"]
        assert report["resourceType"] == "DiagnosticReport"
        assert report["subject"]["reference"] == "Patient/NEURO-001"

    def test_stroke_report_export(self):
        scale = ScaleResult(
            scale_type=ClinicalScaleType.NIHSS,
            score=12,
            max_score=42,
            interpretation="Moderate stroke",
            severity_category="moderate",
        )
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.ACUTE_STROKE,
            findings=["LVO in left MCA"],
            recommendations=["Thrombectomy"],
            severity=SeverityLevel.CRITICAL,
            scale_results=[scale],
        )
        exporter = NeuroReportExporter()
        md = exporter.export_stroke_report(result, patient_id="STK-001")
        assert "Stroke Triage" in md
        assert "NIHSS" in md
        assert "Disclaimer" in md

    def test_scale_report_export(self):
        scale = ScaleResult(
            scale_type=ClinicalScaleType.MOCA,
            score=22,
            max_score=30,
            interpretation="Mild cognitive impairment",
            severity_category="mci",
            recommendations=["Formal neuropsychological testing"],
        )
        exporter = NeuroReportExporter()
        md = exporter.export_scale_report(scale, patient_id="DEM-001")
        assert "MOCA" in md
        assert "22" in md
        assert "Disclaimer" in md
