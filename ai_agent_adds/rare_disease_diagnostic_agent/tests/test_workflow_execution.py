"""Tests for actual workflow runs using the ingest parsers.

Tests end-to-end parser execution and validates output quality.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.ingest.base import IngestRecord, IngestStats
from src.ingest.omim_parser import OMIMParser
from src.ingest.hpo_parser import HPOParser
from src.ingest.orphanet_parser import OrphanetParser
from src.ingest.gene_therapy_parser import GeneTherapyParser


# ===================================================================
# BASE PARSER TESTS
# ===================================================================


class TestBaseParser:
    """Test BaseIngestParser contract."""

    def test_ingest_record_creation(self):
        record = IngestRecord(text="Test record", source="test")
        assert record.text == "Test record"
        assert record.source == "test"

    def test_ingest_record_empty_text_raises(self):
        with pytest.raises(ValueError, match="text must not be empty"):
            IngestRecord(text="")

    def test_ingest_record_whitespace_only_raises(self):
        with pytest.raises(ValueError, match="text must not be empty"):
            IngestRecord(text="   ")

    def test_ingest_record_to_dict(self):
        record = IngestRecord(
            text="Test",
            metadata={"key": "val"},
            collection_name="test_col",
            record_id="R1",
            source="test",
        )
        d = record.to_dict()
        assert d["text"] == "Test"
        assert d["metadata"]["key"] == "val"
        assert d["collection_name"] == "test_col"
        assert d["record_id"] == "R1"

    def test_ingest_stats_defaults(self):
        stats = IngestStats()
        assert stats.source == ""
        assert stats.total_fetched == 0
        assert stats.total_parsed == 0
        assert stats.total_validated == 0
        assert stats.total_errors == 0
        assert stats.duration_seconds == 0.0
        assert stats.error_details == []


# ===================================================================
# OMIM PARSER EXECUTION
# ===================================================================


class TestOMIMExecution:
    """Test full OMIM parser run."""

    def test_run_returns_records_and_stats(self):
        parser = OMIMParser()
        records, stats = parser.run()
        assert isinstance(records, list)
        assert isinstance(stats, IngestStats)

    def test_run_produces_expected_count(self):
        parser = OMIMParser()
        records, stats = parser.run()
        assert stats.total_validated >= 30
        assert len(records) == stats.total_validated

    def test_run_no_errors(self):
        parser = OMIMParser()
        records, stats = parser.run()
        assert stats.total_errors == 0

    def test_run_duration_positive(self):
        parser = OMIMParser()
        records, stats = parser.run()
        assert stats.duration_seconds > 0

    def test_all_records_are_ingest_records(self):
        parser = OMIMParser()
        records, _ = parser.run()
        for r in records:
            assert isinstance(r, IngestRecord)


# ===================================================================
# HPO PARSER EXECUTION
# ===================================================================


class TestHPOExecution:
    """Test full HPO parser run."""

    def test_run_returns_records(self):
        parser = HPOParser()
        records, stats = parser.run()
        assert len(records) >= 50

    def test_run_stats_consistent(self):
        parser = HPOParser()
        records, stats = parser.run()
        assert stats.total_fetched == stats.total_parsed
        assert stats.total_validated == len(records)
        assert stats.total_errors == 0


# ===================================================================
# ORPHANET PARSER EXECUTION
# ===================================================================


class TestOrphanetExecution:
    """Test full Orphanet parser run."""

    def test_run_returns_records(self):
        parser = OrphanetParser()
        records, stats = parser.run()
        assert len(records) >= 30

    def test_run_stats_consistent(self):
        parser = OrphanetParser()
        records, stats = parser.run()
        assert stats.total_fetched == stats.total_parsed
        assert stats.total_validated == len(records)
        assert stats.total_errors == 0


# ===================================================================
# GENE THERAPY PARSER EXECUTION
# ===================================================================


class TestGeneTherapyExecution:
    """Test full Gene Therapy parser run."""

    def test_run_returns_records(self):
        parser = GeneTherapyParser()
        records, stats = parser.run()
        assert len(records) == 25  # 6 approved + 19 investigational

    def test_run_stats_consistent(self):
        parser = GeneTherapyParser()
        records, stats = parser.run()
        assert stats.total_fetched == 25
        assert stats.total_parsed == 25
        assert stats.total_validated == 25
        assert stats.total_errors == 0
