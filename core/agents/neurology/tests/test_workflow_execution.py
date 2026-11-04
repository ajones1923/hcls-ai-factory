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
        """CORRECTED (2026-11-03): this asserted only that the strings "Stroke Triage",
        "NIHSS" and "Disclaimer" appeared, against a hand-built `WorkflowResult`. It would
        have passed unchanged with N4's "IV tPA indicated" sitting in a hemorrhage report --
        which is exactly what shipped, because `NIHSSCalculator` baked treatment advice into
        the interpretation that this exporter prints verbatim, and the stroke workflow runs
        NIHSS for hemorrhagic stroke too.

        The report is now built from a real workflow run so the exporter is exercised against
        what the agent actually produces, and the hemorrhage case is asserted on in
        `test_scale_data_regressions.py::TestNIHSSCarriesNoTreatmentAdvice`.
        """
        from src.clinical_scales import NIHSSCalculator
        from src.clinical_workflows import AcuteStrokeTriageWorkflow

        scores = {item: 0 for item in NIHSSCalculator.ITEMS}
        remaining = 12
        for item, max_val in NIHSSCalculator.ITEMS.items():
            take = min(max_val, remaining)
            scores[item] = take
            remaining -= take
        result = AcuteStrokeTriageWorkflow().run({
            "onset_hours": 2.0,
            "nihss_scores": scores,
            "large_vessel_occlusion": True,
            "prestroke_mrs": 0,
            "ct_affected_regions": ["M1"],
            "blood_pressure_systolic": 150,
            "blood_glucose": 110.0,
            "on_anticoagulation": False,
            "platelet_count": 250,
            "age": 68,
        })
        exporter = NeuroReportExporter()
        md = exporter.export_stroke_report(result, patient_id="STK-001")
        assert "Stroke Triage" in md
        assert "NIHSS" in md
        assert "Disclaimer" in md
        assert "12.0/42.0" in md
        # The scale's own text must not carry treatment advice; the workflow owns that.
        scale_line = next(line for line in md.splitlines() if "**NIHSS:**" in line)
        assert "tPA" not in scale_line
        assert "thrombectomy" not in scale_line.lower()
        # The workflow's own verdict is what carries it, and it is still there.
        assert "alteplase" in md.lower()

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
