"""Tests for decision support and export functionality.

Tests export formats, FHIR R4 generation, AnnData metadata export,
and cross-modal integration.

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    CellTypeAnnotation,
    CellTypeConfidence,
    DrugResponsePrediction,
    ResistanceRisk,
    SCWorkflowType,
    SeverityLevel,
    TMEClass,
    TMEProfile,
    WorkflowResult,
)
from src.export import (
    SCReportExporter,
    VERSION,
    REPORT_TEMPLATES,
    SEVERITY_COLORS,
    _now_iso,
    _now_display,
)


class TestExportHelpers:
    """Tests for export helper functions."""

    def test_now_iso_format(self):
        ts = _now_iso()
        assert "T" in ts
        assert ts.endswith("Z")

    def test_now_display_format(self):
        ts = _now_display()
        assert "UTC" in ts

    def test_version_string(self):
        assert VERSION == "1.0.0"


class TestReportTemplates:
    """Tests for report template definitions."""

    def test_templates_exist(self):
        assert "cell_type" in REPORT_TEMPLATES
        assert "tme" in REPORT_TEMPLATES
        assert "drug_response" in REPORT_TEMPLATES
        assert "workflow" in REPORT_TEMPLATES

    def test_templates_are_strings(self):
        for key, template in REPORT_TEMPLATES.items():
            assert isinstance(template, str), f"Template {key} is not a string"
            assert len(template) > 20


class TestSeverityColors:
    """Tests for severity color map."""

    def test_all_severities_have_colors(self):
        for sev in SeverityLevel:
            assert sev in SEVERITY_COLORS, f"Missing color for {sev.value}"

    def test_colors_are_hex(self):
        for sev, color in SEVERITY_COLORS.items():
            assert color.startswith("#"), f"Color for {sev.value} is not hex"


class TestCellTypeReport:
    """Tests for cell type annotation report export."""

    def test_cell_type_report_generation(self):
        exporter = SCReportExporter()
        annotations = [
            CellTypeAnnotation(
                cluster_id="0",
                cell_type="CD8+ T cell",
                confidence=CellTypeConfidence.HIGH,
                confidence_score=0.95,
                marker_genes=["CD8A", "GZMB"],
                fraction=0.12,
            ),
            CellTypeAnnotation(
                cluster_id="1",
                cell_type="Macrophage",
                confidence=CellTypeConfidence.MEDIUM,
                confidence_score=0.75,
                marker_genes=["CD68", "CD163"],
                fraction=0.08,
            ),
        ]
        md = exporter.export_cell_type_report(annotations, sample_id="DEMO")
        assert "Cell Type Annotation Report" in md
        assert "CD8+ T cell" in md
        assert "Macrophage" in md
        assert "Disclaimer" in md


class TestTMEReport:
    """Tests for TME profiling report export."""

    def test_tme_report_generation(self):
        exporter = SCReportExporter()
        profile = TMEProfile(
            tme_class=TMEClass.HOT_INFLAMED,
            immune_score=0.82,
            stromal_score=0.25,
            exhaustion_signature=0.4,
            cell_type_fractions={"CD8+ T cell": 0.15, "Macrophage": 0.1},
            predicted_immunotherapy_response="responder",
        )
        md = exporter.export_tme_report(profile, sample_id="S-001")
        assert "Tumor Microenvironment" in md
        assert "hot_inflamed" in md
        assert "Disclaimer" in md


class TestDrugResponseReport:
    """Tests for drug response report export."""

    def test_drug_response_report_generation(self):
        exporter = SCReportExporter()
        predictions = [
            DrugResponsePrediction(
                drug_name="Pembrolizumab",
                drug_class="anti-PD-1",
                predicted_sensitivity=0.85,
                resistance_risk=ResistanceRisk.LOW,
            ),
        ]
        md = exporter.export_drug_response_report(predictions, sample_id="S-002")
        assert "Drug Response" in md
        assert "Pembrolizumab" in md


class TestJSONExport:
    """Tests for JSON export."""

    def test_json_export_workflow_result(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.GENERAL,
            severity=SeverityLevel.INFORMATIONAL,
        )
        data = exporter.export_json(result)
        assert data["report_type"] == "single_cell_analysis"
        assert "generated_at" in data
        assert "data" in data

    def test_json_export_dict(self):
        exporter = SCReportExporter()
        data = exporter.export_json({"key": "value"})
        assert data["data"]["key"] == "value"


class TestFHIRExport:
    """Tests for FHIR R4 export."""

    def test_fhir_bundle_structure(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
        )
        bundle = exporter.export_fhir_r4(result, patient_id="P-001")
        assert bundle["resourceType"] == "Bundle"
        assert bundle["type"] == "document"
        assert len(bundle["entry"]) == 1
        assert bundle["entry"][0]["resource"]["resourceType"] == "DiagnosticReport"
        assert bundle["entry"][0]["resource"]["subject"]["reference"] == "Patient/P-001"


class TestAnnDataExport:
    """Tests for AnnData metadata export."""

    def test_anndata_metadata_basic(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[
                CellTypeAnnotation(
                    cluster_id="0",
                    cell_type="T cell",
                    confidence=CellTypeConfidence.HIGH,
                    confidence_score=0.9,
                    marker_genes=["CD3D"],
                    cell_count=500,
                    fraction=0.15,
                ),
            ],
        )
        meta = exporter.export_anndata_metadata(result)
        assert meta["agent"] == "single_cell_intelligence_agent"
        assert "cell_type_annotations" in meta
        assert "0" in meta["cell_type_annotations"]
        assert meta["cell_type_annotations"]["0"]["cell_type"] == "T cell"

    def test_anndata_metadata_tme(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.TME_PROFILING,
            tme_profile=TMEProfile(
                tme_class=TMEClass.COLD_DESERT,
                immune_score=0.1,
                stromal_score=0.7,
            ),
        )
        meta = exporter.export_anndata_metadata(result)
        assert "tme_profile" in meta
        assert meta["tme_profile"]["tme_class"] == "cold_desert"


class TestMarkdownExport:
    """Tests for generic markdown export."""

    def test_markdown_export_workflow(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.GENERAL,
        )
        md = exporter.export_markdown(result)
        assert "Single-Cell Analysis Report" in md
        assert "Disclaimer" in md
