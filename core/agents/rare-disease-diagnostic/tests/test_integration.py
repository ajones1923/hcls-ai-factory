"""Integration tests for cross-module consistency.

Tests that data flows correctly between parsers, collections, settings,
and export modules.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from config.settings import RareDiseaseSettings
from src.models import DiagnosticWorkflowType
from src.ingest.omim_parser import OMIMParser
from src.ingest.hpo_parser import HPOParser
from src.ingest.orphanet_parser import OrphanetParser
from src.ingest.gene_therapy_parser import GeneTherapyParser
from src.export import RareDiseaseReportExporter
from scripts.setup_collections import COLLECTION_SCHEMAS


class TestCrossModuleConsistency:
    """Test cross-module data consistency."""

    def test_settings_collections_match_schemas(self):
        """Collection names in settings should match setup_collections."""
        s = RareDiseaseSettings()
        settings_collections = {
            getattr(s, attr) for attr in dir(s)
            if attr.startswith("COLLECTION_") and isinstance(getattr(s, attr), str)
        }
        schema_collections = set(COLLECTION_SCHEMAS.keys())
        assert settings_collections == schema_collections

    def test_omim_records_target_valid_collection(self):
        parser = OMIMParser()
        records, _ = parser.run()
        valid = set(COLLECTION_SCHEMAS.keys())
        for r in records:
            assert r.collection_name in valid, f"Invalid collection: {r.collection_name}"

    def test_hpo_records_target_valid_collection(self):
        parser = HPOParser()
        records, _ = parser.run()
        valid = set(COLLECTION_SCHEMAS.keys())
        for r in records:
            assert r.collection_name in valid

    def test_orphanet_records_target_valid_collection(self):
        parser = OrphanetParser()
        records, _ = parser.run()
        valid = set(COLLECTION_SCHEMAS.keys())
        for r in records:
            assert r.collection_name in valid

    def test_gene_therapy_records_target_valid_collection(self):
        parser = GeneTherapyParser()
        records, _ = parser.run()
        valid = set(COLLECTION_SCHEMAS.keys())
        for r in records:
            assert r.collection_name in valid

    def test_ingest_stats_consistent(self):
        """Stats should be consistent across all parsers."""
        parsers = [OMIMParser(), HPOParser(), OrphanetParser(), GeneTherapyParser()]
        for parser in parsers:
            records, stats = parser.run()
            assert stats.total_validated == len(records)
            assert stats.total_errors == 0
            assert stats.duration_seconds >= 0

    def test_export_accepts_dict_response(self):
        exporter = RareDiseaseReportExporter()
        response = {
            "findings": ["Finding 1", "Finding 2"],
            "recommendations": ["Rec 1"],
            "confidence": 0.85,
        }
        md = exporter.export_markdown(response)
        assert "Finding 1" in md
        assert "Rec 1" in md

    def test_export_json_includes_version(self):
        exporter = RareDiseaseReportExporter()
        result = exporter.export_json({"findings": ["test"]})
        assert "version" in result
        assert result["report_type"] == "rare_disease_workflow"

    def test_fhir_r4_export_valid(self):
        exporter = RareDiseaseReportExporter()
        bundle = exporter.export_fhir_r4({"findings": ["test finding"]})
        assert bundle["resourceType"] == "Bundle"
        assert bundle["type"] == "document"
        assert len(bundle["entry"]) == 1
        report = bundle["entry"][0]["resource"]
        assert report["resourceType"] == "DiagnosticReport"

    def test_phenopacket_export(self):
        exporter = RareDiseaseReportExporter()
        pp = exporter.export_phenopacket(
            patient_id="TEST-001",
            phenotypes=["HP:0001250", "HP:0001252"],
            diagnosis="Cystic Fibrosis",
            gene="CFTR",
        )
        assert pp["subject"]["id"] == "TEST-001"
        assert len(pp["phenotypicFeatures"]) == 2
        assert pp["metaData"]["phenopacketSchemaVersion"] == "2.0"

    def test_workflow_types_cover_key_scenarios(self):
        """Verify workflow types cover core rare disease scenarios."""
        values = {wf.value for wf in DiagnosticWorkflowType}
        required = {
            "phenotype_driven", "variant_interpretation",
            "differential_diagnosis", "gene_therapy_eligibility",
            "newborn_screening", "metabolic_workup",
        }
        assert required.issubset(values)
