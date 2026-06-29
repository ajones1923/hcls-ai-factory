"""Integration tests for cross-module consistency.

Verifies that data structures, enums, and configurations are consistent
across models, knowledge, collections, ingest, metrics, and export modules.

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    SCWorkflowType,
    WorkflowResult,
    CellTypeAnnotation,
    CellTypeConfidence,
)
from src.collections import (
    WORKFLOW_COLLECTION_WEIGHTS,
    get_all_collection_names,
)
from src.metrics import MetricsCollector, get_metrics_text
from src.export import SCReportExporter, VERSION
from src.ingest.cellxgene_parser import CELL_TYPE_RECORDS, get_cell_type_count
from src.ingest.marker_parser import MARKER_GENE_RECORDS, get_marker_count
from src.ingest.tme_parser import get_tme_profile_count
from config.settings import settings


class TestWorkflowCollectionConsistency:
    """Verify workflow types and collection names are consistent."""

    def test_all_workflows_have_collection_weights(self):
        for wf in SCWorkflowType:
            assert wf in WORKFLOW_COLLECTION_WEIGHTS, (
                f"Workflow {wf.value} missing from WORKFLOW_COLLECTION_WEIGHTS"
            )

    def test_collection_weights_reference_valid_collections(self):
        all_names = set(get_all_collection_names())
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            for coll_name in weights:
                assert coll_name in all_names, (
                    f"Workflow {wf.value} references unknown collection: {coll_name}"
                )

    def test_collection_weight_sums(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 0.05, (
                f"Workflow {wf.value} weights sum to {total}, expected ~1.0"
            )


class TestSettingsCollectionConsistency:
    """Verify settings collection names match defined collections."""

    def test_settings_collection_names_exist(self):
        all_names = set(get_all_collection_names())
        setting_collections = [
            settings.COLLECTION_CELL_TYPES,
            settings.COLLECTION_MARKERS,
            settings.COLLECTION_SPATIAL,
            settings.COLLECTION_TME,
            settings.COLLECTION_DRUG_RESPONSE,
            settings.COLLECTION_LITERATURE,
            settings.COLLECTION_METHODS,
            settings.COLLECTION_DATASETS,
            settings.COLLECTION_TRAJECTORIES,
            settings.COLLECTION_PATHWAYS,
            settings.COLLECTION_CLINICAL,
            settings.COLLECTION_GENOMIC,
        ]
        for name in setting_collections:
            assert name in all_names, (
                f"Settings collection '{name}' not in defined collections"
            )


class TestIngestDataConsistency:
    """Verify ingest seed data is consistent."""

    def test_cell_type_count(self):
        assert get_cell_type_count() >= 30

    def test_marker_count(self):
        assert get_marker_count() >= 50

    def test_tme_profile_count(self):
        assert get_tme_profile_count() >= 10

    def test_cell_types_have_valid_records(self):
        for record in CELL_TYPE_RECORDS:
            assert record.get("cell_type"), "Cell type record missing cell_type"
            assert record.get("markers"), "Cell type record missing markers"

    def test_markers_have_valid_records(self):
        for record in MARKER_GENE_RECORDS:
            assert record.get("gene"), "Marker record missing gene"
            assert record.get("cell_types"), "Marker record missing cell_types"


class TestMetricsConsistency:
    """Verify metrics system works."""

    def test_get_metrics_text(self):
        text = get_metrics_text()
        assert isinstance(text, str)

    def test_metrics_collector_methods(self):
        """All MetricsCollector static methods should be callable."""
        MetricsCollector.record_query("test", 0.1, True)
        MetricsCollector.record_search("sc_cell_types", 0.05, 10)
        MetricsCollector.record_embedding(0.02)
        MetricsCollector.record_workflow("test", 0.5)
        MetricsCollector.record_cell_type_analysis("clustering")
        MetricsCollector.record_spatial_analysis("visium")
        MetricsCollector.record_tme_profile("NSCLC")
        MetricsCollector.record_export("markdown")
        MetricsCollector.set_milvus_status(True)


class TestExportConsistency:
    """Verify export module works with real data."""

    def test_export_version(self):
        assert VERSION == "1.0.0"

    def test_markdown_roundtrip(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.CELL_TYPE_ANNOTATION,
            cell_annotations=[
                CellTypeAnnotation(
                    cluster_id="0",
                    cell_type="T cell",
                    confidence=CellTypeConfidence.HIGH,
                    confidence_score=0.9,
                ),
            ],
        )
        md = exporter.export_markdown(result)
        assert "T cell" in md

    def test_json_roundtrip(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.GENERAL,
        )
        data = exporter.export_json(result)
        assert data["report_type"] == "single_cell_analysis"
        assert data["data"]["workflow_type"] == "general"

    def test_fhir_roundtrip(self):
        exporter = SCReportExporter()
        result = WorkflowResult(
            workflow_type=SCWorkflowType.GENERAL,
        )
        bundle = exporter.export_fhir_r4(result)
        assert bundle["resourceType"] == "Bundle"
