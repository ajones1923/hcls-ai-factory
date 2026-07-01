"""Tests for the RAG engine integration points.

Tests the ingest parsers' IngestRecord output for search and
conversation readiness.

Author: Adam Jones
Date: March 2026
"""


from src.ingest.omim_parser import OMIMParser
from src.ingest.hpo_parser import HPOParser
from src.ingest.orphanet_parser import OrphanetParser
from src.ingest.gene_therapy_parser import GeneTherapyParser


# ===================================================================
# SEARCH READINESS TESTS
# ===================================================================


class TestSearchReadiness:
    """Test that parser output is ready for vector search."""

    def test_omim_records_have_text(self):
        parser = OMIMParser()
        records, stats = parser.run()
        for r in records:
            assert r.text and len(r.text) > 20

    def test_hpo_records_have_text(self):
        parser = HPOParser()
        records, stats = parser.run()
        for r in records:
            assert r.text and len(r.text) > 20

    def test_orphanet_records_have_text(self):
        parser = OrphanetParser()
        records, stats = parser.run()
        for r in records:
            assert r.text and len(r.text) > 20

    def test_gene_therapy_records_have_text(self):
        parser = GeneTherapyParser()
        records, stats = parser.run()
        for r in records:
            assert r.text and len(r.text) > 20

    def test_all_records_have_collection_name(self):
        """All records must specify a target collection."""
        parsers = [OMIMParser(), HPOParser(), OrphanetParser(), GeneTherapyParser()]
        for parser in parsers:
            records, _ = parser.run()
            for r in records:
                assert r.collection_name, f"Missing collection for {r.record_id}"

    def test_all_records_have_source(self):
        parsers = [OMIMParser(), HPOParser(), OrphanetParser(), GeneTherapyParser()]
        for parser in parsers:
            records, _ = parser.run()
            for r in records:
                assert r.source, f"Missing source for {r.record_id}"

    def test_all_records_have_record_id(self):
        parsers = [OMIMParser(), HPOParser(), OrphanetParser(), GeneTherapyParser()]
        for parser in parsers:
            records, _ = parser.run()
            for r in records:
                assert r.record_id, "Missing record_id"


# ===================================================================
# CONVERSATION READINESS TESTS
# ===================================================================


class TestConversationReadiness:
    """Test that records contain enough context for conversation."""

    def test_omim_text_includes_gene(self):
        parser = OMIMParser()
        records, _ = parser.run()
        for r in records:
            assert "Gene:" in r.text

    def test_omim_text_includes_inheritance(self):
        parser = OMIMParser()
        records, _ = parser.run()
        for r in records:
            assert "Inheritance:" in r.text

    def test_hpo_text_includes_definition(self):
        parser = HPOParser()
        records, _ = parser.run()
        for r in records:
            assert "Definition:" in r.text

    def test_gene_therapy_text_includes_mechanism(self):
        parser = GeneTherapyParser()
        records, _ = parser.run()
        for r in records:
            assert "Mechanism:" in r.text
