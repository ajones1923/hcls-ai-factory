"""Tests for workflow execution patterns.

Verifies WorkflowResult creation, serialization, and export roundtrip
for all workflow types.

Author: Adam Jones
Date: March 2026
"""

import pytest

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
from src.export import SCReportExporter


class TestWorkflowResultSerialization:
    """Tests for WorkflowResult serialization."""

    def test_model_dump(self):
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            severity=SeverityLevel.INFORMATIONAL,
        )
        dumped = result.model_dump()
        assert dumped["workflow_type"] == "cell_type_annotation"
        assert dumped["severity"] == "informational"

    def test_model_dump_with_annotations(self):
        ann = CellTypeAnnotation(
            cluster_id="0",
            cell_type="CD8+ T cell",
            confidence=CellTypeConfidence.HIGH,
            confidence_score=0.93,
            marker_genes=["CD8A", "GZMB"],
            cell_count=800,
            fraction=0.12,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[ann],
        )
        dumped = result.model_dump()
        assert len(dumped["cell_annotations"]) == 1
        assert dumped["cell_annotations"][0]["cell_type"] == "CD8+ T cell"
        assert dumped["cell_annotations"][0]["confidence_score"] == 0.93

    def test_model_dump_with_tme(self):
        profile = TMEProfile(
            tme_class=TMEClass.HOT_INFLAMED,
            immune_score=0.8,
            stromal_score=0.3,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.TME_PROFILING,
            tme_profile=profile,
        )
        dumped = result.model_dump()
        assert dumped["tme_profile"]["tme_class"] == "hot_inflamed"

    def test_model_dump_with_drugs(self):
        pred = DrugResponsePrediction(
            drug_name="Nivolumab",
            predicted_sensitivity=0.72,
            resistance_risk=ResistanceRisk.MEDIUM,
        )
        result = WorkflowResult(
            workflow_type=SCWorkflowType.DRUG_RESPONSE,
            drug_predictions=[pred],
        )
        dumped = result.model_dump()
        assert len(dumped["drug_predictions"]) == 1
        assert dumped["drug_predictions"][0]["drug_name"] == "Nivolumab"


class TestWorkflowExportRoundtrip:
    """Tests for workflow result -> export roundtrip."""

    def test_markdown_export(self):
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[
                CellTypeAnnotation(
                    cluster_id="0",
                    cell_type="T cell",
                    confidence=CellTypeConfidence.HIGH,
                    confidence_score=0.9,
                    fraction=0.2,
                ),
            ],
        )
        exporter = SCReportExporter()
        md = exporter.export_markdown(result)
        assert isinstance(md, str)
        assert "T cell" in md

    def test_json_export(self):
        result = WorkflowResult(
            workflow_type=SCWorkflowType.GENERAL,
        )
        exporter = SCReportExporter()
        data = exporter.export_json(result)
        assert "report_type" in data
        assert "data" in data
        assert data["data"]["workflow_type"] == "general"

    def test_fhir_export(self):
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
        )
        exporter = SCReportExporter()
        bundle = exporter.export_fhir_r4(result)
        assert bundle["resourceType"] == "Bundle"

    def test_anndata_metadata_export(self):
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[
                CellTypeAnnotation(
                    cluster_id="0",
                    cell_type="Macrophage",
                    confidence=CellTypeConfidence.HIGH,
                    confidence_score=0.88,
                    cell_count=600,
                    fraction=0.1,
                ),
            ],
        )
        exporter = SCReportExporter()
        meta = exporter.export_anndata_metadata(result)
        assert meta["agent"] == "single_cell_intelligence_agent"
        assert "cell_type_annotations" in meta


class TestAllWorkflowTypes:
    """Verify WorkflowResult can be created for all workflow types."""

    @pytest.mark.parametrize("wf_type", list(SCWorkflowType))
    def test_workflow_result_creation(self, wf_type):
        result = WorkflowResult(workflow_type=wf_type)
        assert result.workflow_type == wf_type
        dumped = result.model_dump()
        assert dumped["workflow_type"] == wf_type.value
