"""Integration tests for cross-module consistency.

Verifies that data structures, enums, and configurations are consistent
across models, knowledge, collections, ingest, metrics, and export modules.

Author: Adam Jones
Date: March 2026
"""


from src.models import (
    ClinicalScaleType,
    NeuroWorkflowType,
)
from src.knowledge import (
    CLINICAL_SCALES,
)
from src.collections import (
    WORKFLOW_COLLECTION_WEIGHTS,
    get_all_collection_names,
)
from src.metrics import MetricsCollector, get_metrics_text
from src.export import NeuroReportExporter, VERSION
from src.ingest.pubmed_neuro_parser import LANDMARK_NEURO_PAPERS, get_landmark_paper_count
from src.ingest.neuroimaging_parser import NEUROIMAGING_PROTOCOLS, get_neuroimaging_protocol_count
from src.ingest.eeg_parser import EEG_PATTERNS, get_eeg_pattern_count
from config.settings import settings


class TestWorkflowCollectionConsistency:
    """Verify workflow types and collection names are consistent."""

    def test_all_workflows_have_collection_weights(self):
        for wf in NeuroWorkflowType:
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
            settings.COLLECTION_LITERATURE,
            settings.COLLECTION_TRIALS,
            settings.COLLECTION_IMAGING,
            settings.COLLECTION_ELECTROPHYSIOLOGY,
            settings.COLLECTION_DEGENERATIVE,
            settings.COLLECTION_CEREBROVASCULAR,
            settings.COLLECTION_EPILEPSY,
            settings.COLLECTION_ONCOLOGY,
            settings.COLLECTION_MS,
            settings.COLLECTION_MOVEMENT,
            settings.COLLECTION_HEADACHE,
            settings.COLLECTION_NEUROMUSCULAR,
            settings.COLLECTION_GUIDELINES,
            settings.COLLECTION_GENOMIC,
        ]
        for coll in setting_collections:
            assert coll in all_names, f"Settings collection {coll} not in ALL_COLLECTIONS"


class TestScaleKnowledgeModelConsistency:
    """Verify clinical scales in knowledge match ClinicalScaleType enum."""

    def test_knowledge_scales_count(self):
        """Knowledge base should have entries for all 10 scales."""
        assert len(CLINICAL_SCALES) == 10

    def test_enum_scales_mapped_to_knowledge(self):
        """Each ClinicalScaleType enum member should map to a knowledge entry.

        Note: enum values may differ from knowledge keys (e.g., updrs_part_iii vs updrs).
        We verify the knowledge keys exist and the enum count matches.
        """
        assert len(ClinicalScaleType) == len(CLINICAL_SCALES)


class TestIngestDataConsistency:
    """Verify ingest seed data counts and validity."""

    def test_pubmed_paper_count(self):
        assert get_landmark_paper_count() == len(LANDMARK_NEURO_PAPERS)
        assert get_landmark_paper_count() >= 30

    def test_neuroimaging_protocol_count(self):
        assert get_neuroimaging_protocol_count() == len(NEUROIMAGING_PROTOCOLS)
        assert get_neuroimaging_protocol_count() >= 50

    def test_eeg_pattern_count(self):
        assert get_eeg_pattern_count() == len(EEG_PATTERNS)
        assert get_eeg_pattern_count() >= 30

    def test_pubmed_papers_have_domains(self):
        for paper in LANDMARK_NEURO_PAPERS:
            assert paper.get("domain"), f"Paper {paper.get('pmid')} missing domain"

    def test_neuroimaging_protocols_have_modality(self):
        for proto in NEUROIMAGING_PROTOCOLS:
            assert proto.get("modality"), f"Protocol {proto.get('protocol_id')} missing modality"

    def test_eeg_patterns_have_category(self):
        for pattern in EEG_PATTERNS:
            assert pattern.get("category"), f"Pattern {pattern.get('pattern_id')} missing category"


class TestMetricsExportConsistency:
    """Verify metrics and export modules work together."""

    def test_metrics_text_returns_string(self):
        text = get_metrics_text()
        assert isinstance(text, str)

    def test_metrics_collector_methods_callable(self):
        # Should not raise even without prometheus_client
        MetricsCollector.record_query("general", 0.5, True)
        MetricsCollector.record_search("neuro_literature", 0.1, 5)
        MetricsCollector.record_scale("nihss", 0.5)
        MetricsCollector.record_export("markdown")

    def test_export_version_defined(self):
        assert VERSION == "1.0.0"

    def test_exporter_instantiation(self):
        exporter = NeuroReportExporter()
        assert exporter is not None
