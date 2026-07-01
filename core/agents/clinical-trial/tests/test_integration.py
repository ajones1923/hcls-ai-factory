"""Cross-module consistency tests for the Clinical Trial Intelligence Agent.

Verifies that:
  - All workflow types appear in relevant maps and models
  - Collection names are consistent across settings and schemas
  - Collection weights sum correctly in both settings and schemas
  - Ingest parsers produce records for defined collections
  - Export formats reference valid model types

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    TrialWorkflowType,
    TrialPhase,
    TrialStatus,
    SeverityLevel,
    TherapeuticArea,
    WorkflowResult,
)
from config.settings import TrialSettings
from scripts.setup_collections import COLLECTION_SCHEMAS, get_collection_names
from src.ingest.clinicaltrials_parser import ClinicalTrialsParser, LANDMARK_TRIALS
from src.ingest.regulatory_parser import RegulatoryParser, REGULATORY_MILESTONES
from src.ingest.pubmed_parser import PubMedTrialParser
from src.metrics import MetricsCollector, get_metrics_text
from src.export import TrialReportExporter, SEVERITY_COLORS, REPORT_TEMPLATES


class TestWorkflowConsistency:
    """Verify workflow types are consistent across modules."""

    def test_all_workflow_types_in_model(self):
        """Every workflow type must be usable in WorkflowResult."""
        for wtype in TrialWorkflowType:
            wr = WorkflowResult(
                workflow_type=wtype,
                findings=[f"Test for {wtype.value}"],
            )
            assert wr.workflow_type == wtype

    def test_all_workflow_types_serializable(self):
        """Every workflow type must serialize to its string value."""
        for wtype in TrialWorkflowType:
            wr = WorkflowResult(workflow_type=wtype)
            d = wr.model_dump()
            assert d["workflow_type"] == wtype.value

    def test_workflow_count(self):
        """There should be 19 workflow types (18 specific + general)."""
        assert len(TrialWorkflowType) == 19


class TestCollectionConsistency:
    """Verify collection names are consistent across settings and schemas."""

    def test_settings_collections_match_schemas(self):
        """All collection names in settings should exist in schemas."""
        s = TrialSettings()
        settings_collections = set()
        for attr in dir(s):
            if attr.startswith("COLLECTION_"):
                val = getattr(s, attr)
                if isinstance(val, str):
                    settings_collections.add(val)

        schema_collections = set(get_collection_names())
        assert settings_collections == schema_collections, (
            f"Mismatch: settings has {settings_collections - schema_collections}, "
            f"schemas has {schema_collections - settings_collections}"
        )

    def test_collection_count(self):
        """There should be 14 collections."""
        assert len(COLLECTION_SCHEMAS) == 14
        s = TrialSettings()
        collection_attrs = [
            attr for attr in dir(s)
            if attr.startswith("COLLECTION_")
        ]
        assert len(collection_attrs) == 14


class TestWeightConsistency:
    """Verify collection weights sum correctly in both sources."""

    def test_settings_weights_sum(self):
        s = TrialSettings()
        weight_attrs = [
            attr for attr in dir(s)
            if attr.startswith("WEIGHT_") and isinstance(getattr(s, attr), float)
        ]
        total = sum(getattr(s, attr) for attr in weight_attrs)
        assert abs(total - 1.0) <= 0.05

    def test_schema_weights_sum(self):
        total = sum(s["search_weight"] for s in COLLECTION_SCHEMAS.values())
        assert abs(total - 1.0) <= 0.05

    def test_weight_counts_match(self):
        """Settings and schemas should have the same number of weights."""
        s = TrialSettings()
        weight_attrs = [
            attr for attr in dir(s)
            if attr.startswith("WEIGHT_") and isinstance(getattr(s, attr), float)
        ]
        assert len(weight_attrs) == len(COLLECTION_SCHEMAS)


class TestIngestParserConsistency:
    """Verify ingest parsers produce records for defined collections."""

    def test_clinicaltrials_target_collection(self):
        parser = ClinicalTrialsParser()
        records = parser.parse(LANDMARK_TRIALS[:1])
        assert records[0].collection_name == "trial_protocols"
        assert "trial_protocols" in COLLECTION_SCHEMAS

    def test_regulatory_target_collection(self):
        parser = RegulatoryParser()
        records = parser.parse(REGULATORY_MILESTONES[:1])
        assert records[0].collection_name == "trial_regulatory"
        assert "trial_regulatory" in COLLECTION_SCHEMAS

    def test_pubmed_target_collection(self):
        parser = PubMedTrialParser()
        # Construct a minimal article dict
        article = {
            "pmid": "12345",
            "title": "A Test Clinical Trial Publication",
            "abstract": "This study evaluated the efficacy of drug X.",
            "journal": "NEJM",
            "year": "2024",
            "mesh_terms": ["Clinical Trial"],
            "publication_types": ["Clinical Trial"],
            "authors": ["Smith J"],
        }
        records = parser.parse([article])
        assert records[0].collection_name == "trial_literature"
        assert "trial_literature" in COLLECTION_SCHEMAS


class TestMetricsConsistency:
    """Verify metrics module is consistent and functional."""

    def test_metrics_collector_methods(self):
        """All MetricsCollector static methods should be callable without error."""
        MetricsCollector.record_query("patient_matching", duration=1.0, success=True)
        MetricsCollector.record_search("trial_protocols", duration=0.1, num_results=5)
        MetricsCollector.record_embedding(0.05)
        MetricsCollector.record_llm_call("claude-sonnet-4-6", 2.0, 500, 200)
        MetricsCollector.record_workflow("safety_signal", 3.0)
        MetricsCollector.record_matching("oncology", 0.85, "phase_iii")
        MetricsCollector.record_safety_signal("high")
        MetricsCollector.record_export("markdown")
        MetricsCollector.record_ingest("clinicaltrials", 10.0, 100, "trial_protocols")
        MetricsCollector.set_agent_info("1.0.0", 14, 11)
        MetricsCollector.set_milvus_status(True)
        MetricsCollector.update_collection_sizes({"trial_protocols": 1000})
        MetricsCollector.record_pipeline_stage("embed", 0.5)
        MetricsCollector.record_milvus_search(0.1)
        MetricsCollector.record_milvus_upsert(0.2)

    def test_get_metrics_text(self):
        """get_metrics_text should return a string."""
        text = get_metrics_text()
        assert isinstance(text, str)


class TestExportConsistency:
    """Verify export module references valid model types."""

    def test_severity_colors_cover_all_levels(self):
        for sev in SeverityLevel:
            assert sev in SEVERITY_COLORS, f"Missing color for {sev}"

    def test_report_templates_exist(self):
        assert len(REPORT_TEMPLATES) >= 4
        assert "trial_match" in REPORT_TEMPLATES
        assert "protocol_analysis" in REPORT_TEMPLATES
        assert "competitive_landscape" in REPORT_TEMPLATES
        assert "safety" in REPORT_TEMPLATES

    def test_exporter_trial_match(self):
        exporter = TrialReportExporter()
        md = exporter.export_trial_match_report(
            matches=[],
            patient_id="P-001",
        )
        assert "Patient-Trial Matching Report" in md
        assert "P-001" in md

    def test_exporter_protocol_analysis(self):
        exporter = TrialReportExporter()
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PROTOCOL_DESIGN,
            findings=["Complex protocol"],
            recommendations=["Simplify"],
            confidence=0.8,
        )
        md = exporter.export_protocol_analysis(wr, trial_id="NCT001")
        assert "Protocol Design Analysis" in md
        assert "NCT001" in md

    def test_exporter_json(self):
        exporter = TrialReportExporter()
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.GENERAL,
            findings=["Test finding"],
        )
        result = exporter.export_json(wr)
        assert "generated_at" in result
        assert "data" in result

    def test_exporter_fhir_r4(self):
        exporter = TrialReportExporter()
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            findings=["3 trials matched"],
        )
        bundle = exporter.export_fhir_r4(wr, patient_id="P-001")
        assert bundle["resourceType"] == "Bundle"
        assert bundle["type"] == "document"
        assert len(bundle["entry"]) >= 1
        report = bundle["entry"][0]["resource"]
        assert report["resourceType"] == "DiagnosticReport"
        assert "P-001" in report["subject"]["reference"]


class TestEnumCompleteness:
    """Verify enums have expected member counts."""

    def test_trial_phase_count(self):
        assert len(TrialPhase) == 7

    def test_trial_status_count(self):
        assert len(TrialStatus) == 7

    def test_therapeutic_area_count(self):
        assert len(TherapeuticArea) == 13

    def test_severity_level_count(self):
        assert len(SeverityLevel) == 5
