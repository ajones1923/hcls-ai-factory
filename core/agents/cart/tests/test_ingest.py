"""Comprehensive tests for CAR-T ingest pipeline parsers.

Covers:
  1. BaseIngestPipeline abstract interface enforcement
  2. Each parser's parse() method with sample input data (mocked external APIs)
  3. Schema validation -- parsed records match expected field names and types
  4. Edge cases: empty responses, malformed data, missing fields
  5. Seed data file loading -- all 13 JSON seed files load without error

Author: Adam Jones
Date: March 2026
"""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Ensure project root is on sys.path
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.ingest.base import BaseIngestPipeline
from src.ingest.literature_parser import PubMedIngestPipeline
from src.ingest.clinical_trials_parser import ClinicalTrialsIngestPipeline
from src.ingest.construct_parser import ConstructIngestPipeline
from src.ingest.assay_parser import AssayIngestPipeline
from src.ingest.manufacturing_parser import ManufacturingIngestPipeline
from src.ingest.safety_parser import SafetyIngestPipeline
from src.ingest.biomarker_parser import BiomarkerIngestPipeline
from src.ingest.regulatory_parser import RegulatoryIngestPipeline
from src.ingest.sequence_parser import SequenceIngestPipeline
from src.ingest.realworld_parser import RealWorldIngestPipeline
from src.ingest.faers_parser import FAERSIngestPipeline
from src.ingest.dailymed_parser import DailyMedIngestPipeline
from src.ingest.uniprot_parser import UniProtIngestPipeline
from src.ingest.cibmtr_parser import CIBMTRIngestPipeline

from src.models import (
    AssayResult,
    AssayType,
    BiomarkerRecord,
    BiomarkerType,
    CARConstruct,
    CARGeneration,
    CARTLiterature,
    CARTStage,
    ClinicalTrial,
    EvidenceLevel,
    FDAStatus,
    ManufacturingRecord,
    ProcessStep,
    RealWorldRecord,
    RegulatoryEvent,
    RegulatoryRecord,
    RWEStudyType,
    SafetyEventType,
    SafetyRecord,
    SequenceRecord,
    SourceType,
    TrialPhase,
    TrialStatus,
)


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════

SEED_DATA_DIR = PROJECT_ROOT / "data" / "reference"

# The project deliberately never publishes data/ — .gitignore excludes it and CLAUDE.md states
# "Data / weights / secrets stay local; only code + docs publish." These tests assert on that
# local-only reference data, so on a fresh clone the files are absent by design. Skip rather than
# fail: a clean clone should report a green run with skips, not 44 failures that look like a broken
# project. Run them on a machine that has the data to exercise them for real.
_DATA_PRESENT = SEED_DATA_DIR.is_dir() and any(SEED_DATA_DIR.glob("*.json"))
pytestmark = pytest.mark.skipif(
    not _DATA_PRESENT,
    reason="local-only reference data under data/ is not published (see .gitignore)",
)


SEED_FILES = [
    "literature_seed_data.json",
    "trials_seed_data.json",
    "constructs_seed_data.json",
    "assay_seed_data.json",
    "manufacturing_seed_data.json",
    "safety_seed_data.json",
    "biomarker_seed_data.json",
    "regulatory_seed_data.json",
    "sequence_seed_data.json",
    "realworld_seed_data.json",
    "patent_seed_data.json",
    "immunogenicity_biomarker_seed.json",
    "immunogenicity_sequence_seed.json",
]


@pytest.fixture
def mock_collection_manager():
    """Return a MagicMock collection manager with sane defaults."""
    manager = MagicMock()
    manager.search.return_value = []
    manager.insert_batch.return_value = None
    manager.connect.return_value = None
    manager.disconnect.return_value = None
    return manager


@pytest.fixture
def mock_embedder():
    """Return a mock embedder that produces 384-dim zero vectors."""
    embedder = MagicMock()
    embedder.encode.return_value = [[0.0] * 384]
    return embedder


@pytest.fixture
def mock_pubmed_client():
    """Return a mock PubMedClient."""
    client = MagicMock()
    client.search.return_value = ["12345678"]
    client.fetch_abstracts.return_value = [
        {
            "pmid": "12345678",
            "title": "CD19 CAR-T cell therapy in B-ALL",
            "abstract": "A phase 1 study of CD19-directed CAR-T cells with 4-1BB costimulatory domain.",
            "authors": ["Smith J", "Doe A"],
            "journal": "Nature Medicine",
            "year": "2023",
            "mesh_terms": ["CAR-T", "CD19", "immunotherapy"],
        }
    ]
    return client


# ═══════════════════════════════════════════════════════════════════════
# 1. BaseIngestPipeline — abstract interface enforcement
# ═══════════════════════════════════════════════════════════════════════


class TestBaseIngestPipeline:
    """Verify that BaseIngestPipeline enforces abstract method contracts."""

    def test_cannot_instantiate_directly(self, mock_collection_manager, mock_embedder):
        """BaseIngestPipeline is abstract and cannot be directly instantiated."""
        with pytest.raises(TypeError):
            BaseIngestPipeline(mock_collection_manager, mock_embedder)

    def test_subclass_missing_fetch_raises(self, mock_collection_manager, mock_embedder):
        """A subclass that omits fetch() cannot be instantiated."""

        class PartialPipeline(BaseIngestPipeline):
            def parse(self, raw_data):
                return []

        with pytest.raises(TypeError):
            PartialPipeline(mock_collection_manager, mock_embedder)

    def test_subclass_missing_parse_raises(self, mock_collection_manager, mock_embedder):
        """A subclass that omits parse() cannot be instantiated."""

        class PartialPipeline(BaseIngestPipeline):
            def fetch(self, **kwargs):
                return []

        with pytest.raises(TypeError):
            PartialPipeline(mock_collection_manager, mock_embedder)

    def test_complete_subclass_instantiates(self, mock_collection_manager, mock_embedder):
        """A subclass implementing both fetch() and parse() can be instantiated."""

        class CompletePipeline(BaseIngestPipeline):
            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = CompletePipeline(mock_collection_manager, mock_embedder)
        assert pipeline.collection_manager is mock_collection_manager
        assert pipeline.embedder is mock_embedder

    def test_embed_and_store_calls_insert_batch(self, mock_collection_manager, mock_embedder):
        """embed_and_store() encodes text and calls insert_batch on the manager."""

        class MinimalPipeline(BaseIngestPipeline):
            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = MinimalPipeline(mock_collection_manager, mock_embedder)

        # Create a mock record with to_embedding_text and model_dump
        mock_record = MagicMock()
        mock_record.to_embedding_text.return_value = "some text"
        mock_record.model_dump.return_value = {"id": "test-1", "text_summary": "test"}

        mock_embedder.encode.return_value = [[0.0] * 384]

        count = pipeline.embed_and_store([mock_record], "test_collection", batch_size=32)
        assert count == 1
        mock_collection_manager.insert_batch.assert_called_once()

    def test_embed_and_store_handles_empty_records(self, mock_collection_manager, mock_embedder):
        """embed_and_store() returns 0 for empty records list."""

        class MinimalPipeline(BaseIngestPipeline):
            def fetch(self, **kwargs):
                return []

            def parse(self, raw_data):
                return []

        pipeline = MinimalPipeline(mock_collection_manager, mock_embedder)
        count = pipeline.embed_and_store([], "test_collection")
        assert count == 0
        mock_collection_manager.insert_batch.assert_not_called()


# ═══════════════════════════════════════════════════════════════════════
# 2. PubMed / Literature Parser
# ═══════════════════════════════════════════════════════════════════════


class TestPubMedParser:
    """Test PubMedIngestPipeline.parse() with sample input."""

    def test_parse_valid_article(self, mock_collection_manager, mock_embedder, mock_pubmed_client):
        pipeline = PubMedIngestPipeline(
            mock_collection_manager, mock_embedder, mock_pubmed_client
        )
        raw = mock_pubmed_client.fetch_abstracts.return_value
        records = pipeline.parse(raw)

        assert len(records) == 1
        record = records[0]
        assert isinstance(record, CARTLiterature)
        assert record.id == "12345678"
        assert record.year == 2023
        assert record.source_type == SourceType.PUBMED
        assert record.target_antigen == "CD19"
        assert record.journal == "Nature Medicine"

    def test_parse_extracts_cart_stage(self, mock_collection_manager, mock_embedder, mock_pubmed_client):
        pipeline = PubMedIngestPipeline(
            mock_collection_manager, mock_embedder, mock_pubmed_client
        )
        raw = [
            {
                "pmid": "99999999",
                "title": "Lentiviral vector production for CAR-T",
                "abstract": "Lentiviral vector transduction efficiency and crispr gene editing.",
                "authors": [],
                "journal": "Mol Ther",
                "year": "2024",
                "mesh_terms": [],
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert records[0].cart_stage == CARTStage.VECTOR_ENG

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder, mock_pubmed_client):
        pipeline = PubMedIngestPipeline(
            mock_collection_manager, mock_embedder, mock_pubmed_client
        )
        records = pipeline.parse([])
        assert records == []

    def test_parse_missing_fields_defaults(self, mock_collection_manager, mock_embedder, mock_pubmed_client):
        """Articles with missing fields should still parse with defaults."""
        pipeline = PubMedIngestPipeline(
            mock_collection_manager, mock_embedder, mock_pubmed_client
        )
        raw = [{"pmid": "00000001"}]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert records[0].year == 2020  # default year
        assert records[0].target_antigen == ""

    def test_parse_invalid_year_defaults(self, mock_collection_manager, mock_embedder, mock_pubmed_client):
        pipeline = PubMedIngestPipeline(
            mock_collection_manager, mock_embedder, mock_pubmed_client
        )
        raw = [
            {
                "pmid": "11111111",
                "title": "Test",
                "abstract": "Abstract",
                "year": "not_a_number",
                "authors": [],
                "journal": "",
                "mesh_terms": [],
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert records[0].year == 2020

    def test_classify_cart_stage_empty_text(self):
        assert PubMedIngestPipeline._classify_cart_stage("") == CARTStage.CLINICAL

    def test_extract_target_antigen_bcma(self):
        text = "This study evaluates BCMA-directed therapy"
        assert PubMedIngestPipeline._extract_target_antigen(text) == "BCMA"

    def test_extract_target_antigen_none(self):
        text = "A general oncology review without specific targets"
        assert PubMedIngestPipeline._extract_target_antigen(text) == ""

    def test_schema_fields(self, mock_collection_manager, mock_embedder, mock_pubmed_client):
        """Verify CARTLiterature record has all expected fields."""
        pipeline = PubMedIngestPipeline(
            mock_collection_manager, mock_embedder, mock_pubmed_client
        )
        raw = mock_pubmed_client.fetch_abstracts.return_value
        records = pipeline.parse(raw)
        records[0]
        expected_fields = {
            "id", "title", "text_chunk", "source_type", "year",
            "cart_stage", "target_antigen", "disease", "keywords", "journal",
        }
        assert expected_fields == set(CARTLiterature.model_fields.keys())


# ═══════════════════════════════════════════════════════════════════════
# 3. ClinicalTrials.gov Parser
# ═══════════════════════════════════════════════════════════════════════


class TestClinicalTrialsParser:
    """Test ClinicalTrialsIngestPipeline.parse() with sample API data."""

    @pytest.fixture
    def sample_study(self):
        """Minimal ClinicalTrials.gov API v2 study JSON."""
        return {
            "protocolSection": {
                "identificationModule": {
                    "nctId": "NCT12345678",
                    "officialTitle": "Phase 1 Study of CD19 CAR-T in DLBCL",
                },
                "descriptionModule": {
                    "briefSummary": "A phase 1 study evaluating CD19-directed CAR-T cells.",
                },
                "designModule": {
                    "phases": ["PHASE1"],
                    "enrollmentInfo": {"count": 50},
                },
                "statusModule": {
                    "overallStatus": "RECRUITING",
                    "startDateStruct": {"date": "2023-01-15"},
                },
                "sponsorCollaboratorsModule": {
                    "leadSponsor": {"name": "Test Pharma"},
                },
                "conditionsModule": {
                    "conditions": ["Diffuse Large B-Cell Lymphoma"],
                },
                "armsInterventionsModule": {
                    "interventions": [
                        {
                            "name": "CD19 CAR-T cells",
                            "description": "Anti-CD19 CAR with 4-1BB domain",
                        }
                    ],
                },
            }
        }

    def test_parse_valid_study(self, mock_collection_manager, mock_embedder, sample_study):
        pipeline = ClinicalTrialsIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse([sample_study])

        assert len(records) == 1
        trial = records[0]
        assert isinstance(trial, ClinicalTrial)
        assert trial.id == "NCT12345678"
        assert trial.phase == TrialPhase.PHASE_1
        assert trial.status == TrialStatus.RECRUITING
        assert trial.sponsor == "Test Pharma"
        assert trial.enrollment == 50
        assert trial.start_year == 2023
        assert trial.target_antigen == "CD19"

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = ClinicalTrialsIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_missing_nct_id_skips(self, mock_collection_manager, mock_embedder):
        """Studies with missing nctId are skipped."""
        pipeline = ClinicalTrialsIngestPipeline(mock_collection_manager, mock_embedder)
        study = {"protocolSection": {"identificationModule": {}}}
        records = pipeline.parse([study])
        assert records == []

    def test_extract_phase_none(self):
        assert ClinicalTrialsIngestPipeline._extract_phase(None) == TrialPhase.NA

    def test_extract_phase_combined(self):
        assert ClinicalTrialsIngestPipeline._extract_phase(["PHASE1", "PHASE2"]) == TrialPhase.PHASE_1_2

    def test_extract_status_unknown(self):
        assert ClinicalTrialsIngestPipeline._extract_status(None) == TrialStatus.UNKNOWN
        assert ClinicalTrialsIngestPipeline._extract_status("GARBAGE") == TrialStatus.UNKNOWN

    def test_parse_malformed_study_skipped(self, mock_collection_manager, mock_embedder):
        """Completely malformed study data is silently skipped."""
        pipeline = ClinicalTrialsIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse([{"garbage": True}])
        assert records == []

    def test_schema_fields(self, mock_collection_manager, mock_embedder, sample_study):
        pipeline = ClinicalTrialsIngestPipeline(mock_collection_manager, mock_embedder)
        pipeline.parse([sample_study])
        expected_fields = {
            "id", "title", "text_summary", "phase", "status", "sponsor",
            "target_antigen", "car_generation", "costimulatory", "disease",
            "enrollment", "start_year", "outcome_summary",
        }
        assert expected_fields == set(ClinicalTrial.model_fields.keys())


# ═══════════════════════════════════════════════════════════════════════
# 4. Construct Parser
# ═══════════════════════════════════════════════════════════════════════


class TestConstructParser:
    """Test ConstructIngestPipeline.parse() and FDA seed data."""

    def test_parse_fda_seed_data(self, mock_collection_manager, mock_embedder):
        ConstructIngestPipeline(mock_collection_manager, mock_embedder)
        fda_constructs = ConstructIngestPipeline._get_fda_approved_constructs()
        assert len(fda_constructs) == 6

        for construct in fda_constructs:
            assert isinstance(construct, CARConstruct)
            assert construct.id != ""
            assert construct.text_summary != ""
            assert construct.fda_status == FDAStatus.APPROVED

    def test_parse_valid_dict(self, mock_collection_manager, mock_embedder):
        pipeline = ConstructIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "test-construct-1",
                "name": "Test CAR",
                "text_summary": "A test CAR construct targeting CD22.",
                "target_antigen": "CD22",
                "generation": CARGeneration.SECOND,
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert records[0].target_antigen == "CD22"

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = ConstructIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_invalid_dict_skipped(self, mock_collection_manager, mock_embedder):
        """Invalid construct data is skipped with a warning."""
        pipeline = ConstructIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [{"garbage": "data"}]
        records = pipeline.parse(raw)
        assert records == []


# ═══════════════════════════════════════════════════════════════════════
# 5. Assay Parser
# ═══════════════════════════════════════════════════════════════════════


class TestAssayParser:
    """Test AssayIngestPipeline.parse() with sample input."""

    def test_parse_valid_assay(self, mock_collection_manager, mock_embedder):
        pipeline = AssayIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "assay-001",
                "text_summary": "Cytotoxicity assay against Nalm-6 cells.",
                "assay_type": "cytotoxicity",
                "construct_id": "fda-tisagenlecleucel",
                "target_antigen": "CD19",
                "cell_line": "Nalm-6",
                "effector_ratio": "10:1",
                "key_metric": "specific_lysis",
                "metric_value": 85.0,
                "outcome": "success",
                "notes": "",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], AssayResult)
        assert records[0].assay_type == AssayType.CYTOTOXICITY
        assert records[0].metric_value == 85.0

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = AssayIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_string_metric_value(self, mock_collection_manager, mock_embedder):
        """metric_value should be converted from string to float."""
        pipeline = AssayIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "assay-002",
                "text_summary": "Test assay.",
                "assay_type": "cytokine",
                "metric_value": "42.5",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert records[0].metric_value == 42.5

    def test_parse_invalid_assay_type_skipped(self, mock_collection_manager, mock_embedder):
        pipeline = AssayIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "assay-bad",
                "text_summary": "Bad assay.",
                "assay_type": "nonexistent_type",
                "metric_value": 0.0,
            }
        ]
        records = pipeline.parse(raw)
        assert records == []


# ═══════════════════════════════════════════════════════════════════════
# 6. Manufacturing Parser
# ═══════════════════════════════════════════════════════════════════════


class TestManufacturingParser:
    """Test ManufacturingIngestPipeline.parse() with sample input."""

    def test_parse_valid_record(self, mock_collection_manager, mock_embedder):
        pipeline = ManufacturingIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "mfg-001",
                "text_summary": "Lentiviral transduction at MOI 5.",
                "process_step": "transduction",
                "vector_type": "lentiviral",
                "parameter": "MOI",
                "parameter_value": "5",
                "target_spec": ">=40%",
                "met_spec": "yes",
                "batch_id": "BATCH-001",
                "notes": "",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], ManufacturingRecord)
        assert records[0].process_step == ProcessStep.TRANSDUCTION

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = ManufacturingIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_invalid_step_skipped(self, mock_collection_manager, mock_embedder):
        pipeline = ManufacturingIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "mfg-bad",
                "text_summary": "Bad record.",
                "process_step": "nonexistent_step",
            }
        ]
        records = pipeline.parse(raw)
        assert records == []


# ═══════════════════════════════════════════════════════════════════════
# 7. Safety Parser
# ═══════════════════════════════════════════════════════════════════════


class TestSafetyParser:
    """Test SafetyIngestPipeline.parse() with sample input."""

    def test_parse_valid_record(self, mock_collection_manager, mock_embedder):
        pipeline = SafetyIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "safety-001",
                "text_summary": "CRS event after Kymriah infusion.",
                "product": "Kymriah",
                "event_type": "CRS",
                "severity_grade": "Grade 3",
                "year": 2023,
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], SafetyRecord)
        assert records[0].event_type == SafetyEventType.CRS
        assert records[0].year == 2023

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = SafetyIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_string_year_converted(self, mock_collection_manager, mock_embedder):
        pipeline = SafetyIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "safety-002",
                "text_summary": "Test safety record.",
                "event_type": "ICANS",
                "year": "2024",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert records[0].year == 2024


# ═══════════════════════════════════════════════════════════════════════
# 8. Biomarker Parser
# ═══════════════════════════════════════════════════════════════════════


class TestBiomarkerParser:
    """Test BiomarkerIngestPipeline.parse() with sample input."""

    def test_parse_valid_record(self, mock_collection_manager, mock_embedder):
        pipeline = BiomarkerIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "bio-001",
                "text_summary": "IL-6 as a predictive biomarker for CRS severity.",
                "biomarker_name": "IL-6",
                "biomarker_type": "predictive",
                "evidence_level": "validated",
                "assay_method": "ELISA",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], BiomarkerRecord)
        assert records[0].biomarker_type == BiomarkerType.PREDICTIVE
        assert records[0].evidence_level == EvidenceLevel.VALIDATED

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = BiomarkerIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []


# ═══════════════════════════════════════════════════════════════════════
# 9. Regulatory Parser
# ═══════════════════════════════════════════════════════════════════════


class TestRegulatoryParser:
    """Test RegulatoryIngestPipeline.parse() with sample input."""

    def test_parse_valid_record(self, mock_collection_manager, mock_embedder):
        pipeline = RegulatoryIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "reg-001",
                "text_summary": "FDA full approval of Kymriah for B-ALL.",
                "product": "Kymriah",
                "regulatory_event": "full_approval",
                "date": "2017-08",
                "agency": "FDA",
                "indication": "r/r B-ALL",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], RegulatoryRecord)
        assert records[0].regulatory_event == RegulatoryEvent.FULL_APPROVAL

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = RegulatoryIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []


# ═══════════════════════════════════════════════════════════════════════
# 10. Sequence Parser
# ═══════════════════════════════════════════════════════════════════════


class TestSequenceParser:
    """Test SequenceIngestPipeline.parse() with sample input."""

    def test_parse_valid_record(self, mock_collection_manager, mock_embedder):
        pipeline = SequenceIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "seq-001",
                "text_summary": "FMC63 scFv targeting CD19.",
                "construct_name": "FMC63",
                "target_antigen": "CD19",
                "scfv_clone": "FMC63",
                "species_origin": "murine",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], SequenceRecord)
        assert records[0].target_antigen == "CD19"

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = SequenceIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []


# ═══════════════════════════════════════════════════════════════════════
# 11. Real-World Evidence Parser
# ═══════════════════════════════════════════════════════════════════════


class TestRealWorldParser:
    """Test RealWorldIngestPipeline.parse() with sample input."""

    def test_parse_valid_record(self, mock_collection_manager, mock_embedder):
        pipeline = RealWorldIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "id": "rwe-001",
                "text_summary": "Registry analysis of Yescarta in LBCL.",
                "study_type": "registry",
                "data_source": "CIBMTR",
                "product": "Yescarta",
                "population_size": "100",
                "median_followup_months": "12.5",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], RealWorldRecord)
        assert records[0].study_type == RWEStudyType.REGISTRY
        assert records[0].population_size == 100
        assert records[0].median_followup_months == 12.5

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = RealWorldIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []


# ═══════════════════════════════════════════════════════════════════════
# 12. FAERS Parser
# ═══════════════════════════════════════════════════════════════════════


class TestFAERSParser:
    """Test FAERSIngestPipeline.parse() with sample openFDA data."""

    @pytest.fixture
    def sample_faers_event(self):
        """Minimal openFDA adverse event JSON."""
        return {
            "safetyreportid": "US-FDA-2023-001",
            "receivedate": "20230615",
            "seriousnesshospitalization": "1",
            "patient": {
                "drug": [
                    {
                        "openfda": {
                            "brand_name": ["KYMRIAH"],
                        },
                        "medicinalproduct": "tisagenlecleucel",
                    }
                ],
                "reaction": [
                    {
                        "reactionmeddrapt": "cytokine release syndrome",
                        "reactionoutcome": "1",
                    },
                    {
                        "reactionmeddrapt": "neutropenia",
                        "reactionoutcome": "2",
                    },
                ],
            },
        }

    def test_parse_valid_event(self, mock_collection_manager, mock_embedder, sample_faers_event):
        pipeline = FAERSIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse([sample_faers_event])

        assert len(records) == 1
        record = records[0]
        assert isinstance(record, SafetyRecord)
        assert record.id.startswith("FAERS-")
        assert record.event_type == SafetyEventType.CRS  # CRS is prioritized
        assert record.year == 2023
        assert "Kymriah" in record.product
        assert record.reporting_source == "FAERS"

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = FAERSIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_event_no_reactions_skipped(self, mock_collection_manager, mock_embedder):
        """Events with no reactions are skipped."""
        pipeline = FAERSIngestPipeline(mock_collection_manager, mock_embedder)
        event = {
            "safetyreportid": "US-FDA-2023-002",
            "patient": {"drug": [], "reaction": []},
        }
        records = pipeline.parse([event])
        assert records == []

    def test_classify_event_type_crs(self):
        assert FAERSIngestPipeline._classify_event_type(
            ["cytokine release syndrome"]
        ) == SafetyEventType.CRS

    def test_classify_event_type_icans(self):
        assert FAERSIngestPipeline._classify_event_type(
            ["neurotoxicity", "neutropenia"]
        ) == SafetyEventType.ICANS

    def test_classify_event_type_fallback(self):
        assert FAERSIngestPipeline._classify_event_type(
            ["unknown_reaction"]
        ) == SafetyEventType.ORGAN_TOXICITY

    def test_extract_severity_fatal(self):
        assert "fatal" in FAERSIngestPipeline._extract_severity(
            {"seriousnessdeath": "1"}
        )

    def test_extract_severity_non_serious(self):
        assert FAERSIngestPipeline._extract_severity({}) == "non-serious"

    def test_extract_cart_product_unspecified(self):
        assert FAERSIngestPipeline._extract_cart_product([]) == "CAR-T product (unspecified)"

    def test_extract_cart_product_by_brand(self):
        drugs = [{"openfda": {"brand_name": ["YESCARTA"]}}]
        result = FAERSIngestPipeline._extract_cart_product(drugs)
        assert "Yescarta" in result

    def test_schema_fields(self, mock_collection_manager, mock_embedder, sample_faers_event):
        pipeline = FAERSIngestPipeline(mock_collection_manager, mock_embedder)
        pipeline.parse([sample_faers_event])
        expected_fields = {
            "id", "text_summary", "product", "event_type", "severity_grade",
            "onset_timing", "incidence_rate", "management_protocol",
            "outcome", "reporting_source", "year",
        }
        assert expected_fields == set(SafetyRecord.model_fields.keys())


# ═══════════════════════════════════════════════════════════════════════
# 13. DailyMed Parser
# ═══════════════════════════════════════════════════════════════════════


class TestDailyMedParser:
    """Test DailyMedIngestPipeline.parse() with sample and fallback data."""

    def test_parse_fallback_data(self, mock_collection_manager, mock_embedder):
        """Verify parsing of the static fallback seed data."""
        from src.ingest.dailymed_parser import _FALLBACK_SEED_DATA

        pipeline = DailyMedIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse(_FALLBACK_SEED_DATA)

        assert len(records) == 6
        for record in records:
            assert isinstance(record, RegulatoryRecord)
            assert record.id != ""
            assert record.text_summary != ""
            assert record.agency == "FDA"

    def test_parse_live_api_format(self, mock_collection_manager, mock_embedder):
        """Verify parsing of live DailyMed API response format."""
        pipeline = DailyMedIngestPipeline(mock_collection_manager, mock_embedder)
        raw = [
            {
                "setid": "abc-123-def",
                "title": "KYMRIAH - tisagenlecleucel suspension",
                "_queried_product": "kymriah",
                "published_date": "20240315",
            }
        ]
        records = pipeline.parse(raw)
        assert len(records) == 1
        assert isinstance(records[0], RegulatoryRecord)
        assert "DAILYMED-" in records[0].id
        assert "Kymriah" in records[0].product

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = DailyMedIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_normalize_date_yyyymmdd(self):
        assert DailyMedIngestPipeline._normalize_date("20240315") == "2024-03"

    def test_normalize_date_empty(self):
        assert DailyMedIngestPipeline._normalize_date("") == ""

    def test_normalize_date_yyyy_mm_dd(self):
        assert DailyMedIngestPipeline._normalize_date("2024-03-15") == "2024-03"

    def test_resolve_product_name_known(self):
        result = DailyMedIngestPipeline._resolve_product_name("kymriah", "")
        assert "Kymriah" in result and "tisagenlecleucel" in result

    def test_resolve_product_name_unknown(self):
        result = DailyMedIngestPipeline._resolve_product_name("unknown", "")
        assert result == "Unknown CAR-T product"


# ═══════════════════════════════════════════════════════════════════════
# 14. UniProt Parser
# ═══════════════════════════════════════════════════════════════════════


class TestUniProtParser:
    """Test UniProtIngestPipeline.parse() with sample UniProt data."""

    @pytest.fixture
    def sample_uniprot_entry(self):
        """Minimal UniProt protein entry JSON."""
        return {
            "primaryAccession": "P15391",
            "proteinDescription": {
                "recommendedName": {
                    "fullName": {"value": "B-lymphocyte antigen CD19"},
                },
            },
            "genes": [
                {"geneName": {"value": "CD19"}},
            ],
            "organism": {
                "scientificName": "Homo sapiens",
            },
            "sequence": {
                "value": "MPPPRLLFFL" * 10,
                "length": 556,
                "molWeight": 61120,
            },
            "features": [
                {
                    "type": "Domain",
                    "description": "Ig-like C2-type 1",
                    "location": {"start": {"value": 20}, "end": {"value": 110}},
                },
            ],
            "comments": [
                {
                    "commentType": "FUNCTION",
                    "texts": [
                        {"value": "Critical signal transduction molecule in B-lymphocytes."}
                    ],
                }
            ],
        }

    def test_parse_valid_entry(self, mock_collection_manager, mock_embedder, sample_uniprot_entry):
        pipeline = UniProtIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse([sample_uniprot_entry])

        assert len(records) == 1
        record = records[0]
        assert isinstance(record, SequenceRecord)
        assert record.id == "UNIPROT-P15391"
        assert record.target_antigen == "CD19"
        assert record.species_origin == "fully_human"
        assert record.immunogenicity_risk == "low"
        assert "CD19" in record.text_summary

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = UniProtIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_missing_accession_skipped(self, mock_collection_manager, mock_embedder):
        pipeline = UniProtIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse([{"no_accession": True}])
        assert records == []

    def test_map_species_murine(self):
        assert UniProtIngestPipeline._map_species_origin("Mus musculus") == "murine"

    def test_map_species_human(self):
        assert UniProtIngestPipeline._map_species_origin("Homo sapiens") == "fully_human"

    def test_extract_next_link_found(self):
        header = '<https://rest.uniprot.org/next?cursor=abc>; rel="next"'
        assert UniProtIngestPipeline._extract_next_link(header) == "https://rest.uniprot.org/next?cursor=abc"

    def test_extract_next_link_not_found(self):
        assert UniProtIngestPipeline._extract_next_link("") is None

    def test_extract_protein_name_recommended(self):
        desc = {"recommendedName": {"fullName": {"value": "Test Protein"}}}
        assert UniProtIngestPipeline._extract_protein_name(desc) == "Test Protein"

    def test_extract_protein_name_fallback(self):
        desc = {}
        assert UniProtIngestPipeline._extract_protein_name(desc) == "Unknown protein"


# ═══════════════════════════════════════════════════════════════════════
# 15. CIBMTR Parser
# ═══════════════════════════════════════════════════════════════════════


class TestCIBMTRParser:
    """Test CIBMTRIngestPipeline.parse() with curated data."""

    def test_parse_curated_data(self, mock_collection_manager, mock_embedder):
        """Verify all curated CIBMTR records parse correctly."""
        from src.ingest.cibmtr_parser import _CURATED_CIBMTR_DATA

        pipeline = CIBMTRIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse(_CURATED_CIBMTR_DATA)

        assert len(records) == 10
        for record in records:
            assert isinstance(record, RealWorldRecord)
            assert record.id != ""
            assert record.text_summary != ""
            assert record.data_source == "CIBMTR"
            assert record.study_type == RWEStudyType.REGISTRY

    def test_parse_empty_input(self, mock_collection_manager, mock_embedder):
        pipeline = CIBMTRIngestPipeline(mock_collection_manager, mock_embedder)
        assert pipeline.parse([]) == []

    def test_parse_missing_text_summary_skipped(self, mock_collection_manager, mock_embedder):
        pipeline = CIBMTRIngestPipeline(mock_collection_manager, mock_embedder)
        records = pipeline.parse([{"id": "test", "text_summary": ""}])
        assert records == []

    @patch("src.ingest.cibmtr_parser.CIBMTRIngestPipeline._request_with_retry")
    def test_fetch_falls_back_on_network_failure(
        self, mock_request, mock_collection_manager, mock_embedder
    ):
        """When CIBMTR website is unreachable, fetch returns curated data."""
        mock_request.return_value = None

        pipeline = CIBMTRIngestPipeline(mock_collection_manager, mock_embedder)
        data = pipeline.fetch()

        assert len(data) == 10  # curated data count


# ═══════════════════════════════════════════════════════════════════════
# 16. Seed Data File Loading
# ═══════════════════════════════════════════════════════════════════════


class TestSeedDataFiles:
    """Verify all 13 JSON seed files load without error and have valid structure."""

    @pytest.mark.parametrize("filename", SEED_FILES)
    def test_seed_file_loads(self, filename):
        """Each seed file should load as valid JSON."""
        filepath = SEED_DATA_DIR / filename
        assert filepath.exists(), f"Seed file missing: {filepath}"

        with open(filepath, "r") as f:
            data = json.load(f)

        assert isinstance(data, list), f"{filename} should contain a JSON array"
        assert len(data) > 0, f"{filename} should not be empty"

    @pytest.mark.parametrize("filename", SEED_FILES)
    def test_seed_records_have_id(self, filename):
        """Every record in every seed file must have a non-empty 'id' field."""
        filepath = SEED_DATA_DIR / filename
        with open(filepath, "r") as f:
            data = json.load(f)

        for i, record in enumerate(data):
            assert "id" in record, f"{filename}[{i}] missing 'id' field"
            assert record["id"], f"{filename}[{i}] has empty 'id' field"

    @pytest.mark.parametrize("filename", [
        f for f in SEED_FILES
        if f not in ("patent_seed_data.json",)  # patents may use text_chunk
    ])
    def test_seed_records_have_text_summary_or_chunk(self, filename):
        """Every record should have a non-empty text field for embedding."""
        filepath = SEED_DATA_DIR / filename
        with open(filepath, "r") as f:
            data = json.load(f)

        for i, record in enumerate(data):
            has_text = bool(
                record.get("text_summary")
                or record.get("text_chunk")
            )
            assert has_text, (
                f"{filename}[{i}] (id={record.get('id', '?')}) "
                "missing both 'text_summary' and 'text_chunk'"
            )
