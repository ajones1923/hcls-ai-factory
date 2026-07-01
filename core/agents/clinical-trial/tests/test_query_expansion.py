"""Tests for ingest parser base functionality (serves as query expansion tests).

Since the clinical trial agent uses ingest parsers rather than a separate
query_expansion module, these tests validate the core parsing and validation
logic that expands raw data into structured records.

Covers:
  - IngestRecord creation and validation
  - IngestStats tracking
  - BaseIngestParser ABC contract
  - Parser-specific field extraction

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.ingest.base import IngestRecord, IngestStats
from src.ingest.clinicaltrials_parser import ClinicalTrialsParser, LANDMARK_TRIALS
from src.ingest.pubmed_parser import PubMedTrialParser, TRIAL_MESH_TERMS
from src.ingest.regulatory_parser import RegulatoryParser, REGULATORY_MILESTONES


class TestIngestRecord:
    """Test IngestRecord dataclass."""

    def test_valid_creation(self):
        r = IngestRecord(
            text="Test clinical trial record",
            metadata={"nct_id": "NCT00000001"},
            collection_name="trial_protocols",
            record_id="NCT00000001",
            source="clinicaltrials",
        )
        assert r.text == "Test clinical trial record"
        assert r.metadata["nct_id"] == "NCT00000001"

    def test_empty_text_raises(self):
        with pytest.raises(ValueError):
            IngestRecord(text="")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValueError):
            IngestRecord(text="   ")

    def test_to_dict(self):
        r = IngestRecord(
            text="Test record",
            metadata={"key": "value"},
            collection_name="test_collection",
            record_id="R001",
            source="test",
        )
        d = r.to_dict()
        assert d["text"] == "Test record"
        assert d["metadata"]["key"] == "value"
        assert d["record_id"] == "R001"

    def test_default_values(self):
        r = IngestRecord(text="Minimal record")
        assert r.metadata == {}
        assert r.collection_name == ""
        assert r.record_id is None
        assert r.source == ""


class TestIngestStats:
    """Test IngestStats dataclass."""

    def test_defaults(self):
        s = IngestStats()
        assert s.source == ""
        assert s.total_fetched == 0
        assert s.total_parsed == 0
        assert s.total_validated == 0
        assert s.total_errors == 0
        assert s.duration_seconds == 0.0
        assert s.error_details == []

    def test_custom_values(self):
        s = IngestStats(
            source="pubmed",
            total_fetched=100,
            total_parsed=95,
            total_validated=90,
            total_errors=5,
            duration_seconds=12.5,
        )
        assert s.source == "pubmed"
        assert s.total_validated == 90


class TestClinicalTrialsParserInit:
    """Test ClinicalTrialsParser initialization."""

    def test_default_init(self):
        parser = ClinicalTrialsParser()
        assert parser.source_name == "clinicaltrials"
        assert parser.api_key is None

    def test_with_api_key(self):
        parser = ClinicalTrialsParser(api_key="test_key")
        assert parser.api_key == "test_key"


class TestClinicalTrialsParserParse:
    """Test ClinicalTrialsParser parse functionality."""

    def test_parse_landmark_trials(self):
        parser = ClinicalTrialsParser()
        records = parser.parse(LANDMARK_TRIALS)
        assert len(records) == len(LANDMARK_TRIALS)

    def test_parse_landmark_record_fields(self):
        parser = ClinicalTrialsParser()
        records = parser.parse(LANDMARK_TRIALS[:1])
        assert len(records) == 1
        r = records[0]
        assert r.collection_name == "trial_protocols"
        assert r.source == "clinicaltrials"
        assert r.record_id.startswith("NCT")
        assert "title" in r.metadata
        assert r.metadata.get("is_landmark") is True

    def test_validate_valid_record(self):
        parser = ClinicalTrialsParser()
        record = IngestRecord(
            text="A sufficiently long clinical trial record text",
            metadata={"title": "Test Trial"},
            record_id="NCT00000001",
        )
        assert parser.validate_record(record) is True

    def test_validate_short_text(self):
        parser = ClinicalTrialsParser()
        record = IngestRecord(
            text="Too short",
            metadata={"title": "Test"},
            record_id="NCT1",
        )
        assert parser.validate_record(record) is False

    def test_validate_no_record_id(self):
        parser = ClinicalTrialsParser()
        record = IngestRecord(
            text="A sufficiently long clinical trial record text",
            metadata={"title": "Test"},
        )
        assert parser.validate_record(record) is False

    def test_validate_no_title(self):
        parser = ClinicalTrialsParser()
        record = IngestRecord(
            text="A sufficiently long clinical trial record text",
            metadata={},
            record_id="NCT1",
        )
        assert parser.validate_record(record) is False


class TestClinicalTrialsParserSeed:
    """Test seed_landmark_trials functionality."""

    def test_seed_returns_list(self):
        parser = ClinicalTrialsParser()
        trials = parser.seed_landmark_trials()
        assert isinstance(trials, list)
        assert len(trials) >= 20

    def test_seed_returns_copy(self):
        parser = ClinicalTrialsParser()
        trials1 = parser.seed_landmark_trials()
        trials2 = parser.seed_landmark_trials()
        assert trials1 is not trials2


class TestPubMedParserInit:
    """Test PubMedTrialParser initialization."""

    def test_default_init(self):
        parser = PubMedTrialParser()
        assert parser.source_name == "pubmed"

    def test_mesh_terms_available(self):
        assert len(TRIAL_MESH_TERMS) >= 5
        assert "Clinical Trial" in TRIAL_MESH_TERMS
        assert "Randomized Controlled Trial" in TRIAL_MESH_TERMS


class TestPubMedParserValidation:
    """Test PubMedTrialParser validation."""

    def test_validate_valid_record(self):
        parser = PubMedTrialParser()
        record = IngestRecord(
            text="Title: A Phase III Randomized Clinical Trial of Drug X",
            metadata={"title": "Test Article"},
            record_id="PMID:12345678",
        )
        assert parser.validate_record(record) is True

    def test_validate_empty_text(self):
        parser = PubMedTrialParser()
        record = IngestRecord(
            text="Short text only",
            metadata={"title": "T"},
            record_id="PMID:1",
        )
        assert parser.validate_record(record) is False


class TestRegulatoryParserInit:
    """Test RegulatoryParser initialization."""

    def test_default_init(self):
        parser = RegulatoryParser()
        assert parser.source_name == "regulatory"


class TestRegulatoryParserParse:
    """Test RegulatoryParser parse functionality."""

    def test_parse_milestones(self):
        parser = RegulatoryParser()
        records = parser.parse(REGULATORY_MILESTONES)
        assert len(records) == len(REGULATORY_MILESTONES)

    def test_parse_record_fields(self):
        parser = RegulatoryParser()
        records = parser.parse(REGULATORY_MILESTONES[:1])
        r = records[0]
        assert r.collection_name == "trial_regulatory"
        assert r.source == "regulatory"
        assert "drug" in r.metadata
        assert "agency" in r.metadata
        assert "decision" in r.metadata

    def test_validate_valid_record(self):
        parser = RegulatoryParser()
        record = IngestRecord(
            text="Regulatory Decision: FDA Approval for Drug X for Indication Y",
            metadata={"drug": "Drug X", "agency": "FDA"},
            record_id="REG:FDA:DrugX:2024-01-01",
        )
        assert parser.validate_record(record) is True

    def test_validate_missing_drug(self):
        parser = RegulatoryParser()
        record = IngestRecord(
            text="Regulatory Decision: FDA Approval for Drug X",
            metadata={"agency": "FDA"},
            record_id="REG:1",
        )
        assert parser.validate_record(record) is False
