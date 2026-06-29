"""Tests for CAR-T Intelligence Agent export module.

Validates Markdown, JSON, and PDF export functions, filename generation,
and handling of comparative results.

Author: Adam Jones
Date: February 2026
"""

import json

import pytest

from src.export import export_json, export_markdown, export_pdf, generate_filename
from src.models import ComparativeResult, CrossCollectionResult, SearchHit


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def basic_evidence():
    """Return a simple CrossCollectionResult for export tests."""
    return CrossCollectionResult(
        query="What is the efficacy of CD19 CAR-T?",
        hits=[
            SearchHit(
                collection="Literature",
                id="12345678",
                score=0.92,
                text="CD19 CAR-T therapy achieves high response rates in B-ALL.",
                metadata={"title": "CD19 study", "year": 2023, "target_antigen": "CD19"},
            ),
            SearchHit(
                collection="Trial",
                id="NCT03958656",
                score=0.87,
                text="Phase 2 tisagenlecleucel in pediatric B-ALL.",
                metadata={"phase": "Phase 2", "status": "Completed"},
            ),
        ],
        knowledge_context="## Target: CD19\nB-cell antigen",
        total_collections_searched=10,
        search_time_ms=42.5,
    )


@pytest.fixture
def comparative_evidence():
    """Return a ComparativeResult for export tests."""
    ev_a = CrossCollectionResult(
        query="CD19 evidence",
        hits=[
            SearchHit(
                collection="Literature",
                id="11111111",
                score=0.90,
                text="CD19 CAR-T data.",
                metadata={"title": "CD19 paper", "year": 2023},
            ),
        ],
        total_collections_searched=10,
        search_time_ms=20.0,
    )
    ev_b = CrossCollectionResult(
        query="BCMA evidence",
        hits=[
            SearchHit(
                collection="Literature",
                id="22222222",
                score=0.88,
                text="BCMA CAR-T data.",
                metadata={"title": "BCMA paper", "year": 2024},
            ),
        ],
        total_collections_searched=10,
        search_time_ms=22.0,
    )
    return ComparativeResult(
        query="Compare CD19 vs BCMA",
        entity_a="CD19",
        entity_b="BCMA",
        evidence_a=ev_a,
        evidence_b=ev_b,
        comparison_context="CD19 targets B-cells; BCMA targets plasma cells.",
        total_search_time_ms=45.0,
    )


# ═══════════════════════════════════════════════════════════════════════
# MARKDOWN EXPORT
# ═══════════════════════════════════════════════════════════════════════


class TestExportMarkdown:
    """Tests for the export_markdown() function."""

    def test_returns_non_empty_string(self, basic_evidence):
        """export_markdown() returns a non-empty string."""
        md = export_markdown(
            query="CD19 CAR-T efficacy",
            response_text="CD19 therapy is effective.",
            evidence=basic_evidence,
        )
        assert isinstance(md, str)
        assert len(md) > 0

    def test_contains_query(self, basic_evidence):
        """The exported Markdown contains the original query."""
        md = export_markdown(
            query="CD19 CAR-T efficacy",
            response_text="Answer here.",
            evidence=basic_evidence,
        )
        assert "CD19 CAR-T efficacy" in md

    def test_contains_response(self, basic_evidence):
        """The exported Markdown contains the LLM response text."""
        md = export_markdown(
            query="test",
            response_text="The response is here.",
            evidence=basic_evidence,
        )
        assert "The response is here." in md

    def test_contains_evidence_section(self, basic_evidence):
        """The exported Markdown includes an evidence section header."""
        md = export_markdown(
            query="test",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        assert "Evidence Sources" in md

    def test_contains_footer(self, basic_evidence):
        """The exported Markdown includes the agent footer."""
        md = export_markdown(
            query="test",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        assert "CAR-T Intelligence Agent" in md

    def test_handles_no_evidence(self):
        """export_markdown() works when evidence is None."""
        md = export_markdown(
            query="test",
            response_text="No evidence available.",
            evidence=None,
        )
        assert isinstance(md, str)
        assert len(md) > 0

    def test_handles_comparative_result(self, comparative_evidence):
        """export_markdown() renders comparative results correctly."""
        md = export_markdown(
            query="Compare CD19 vs BCMA",
            response_text="Comparison analysis.",
            comp_result=comparative_evidence,
        )
        assert "CD19" in md
        assert "BCMA" in md
        assert "Comparative" in md


# ═══════════════════════════════════════════════════════════════════════
# JSON EXPORT
# ═══════════════════════════════════════════════════════════════════════


class TestExportJson:
    """Tests for the export_json() function."""

    def test_returns_valid_json(self, basic_evidence):
        """export_json() returns a valid JSON string."""
        json_str = export_json(
            query="CD19 efficacy",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        data = json.loads(json_str)
        assert isinstance(data, dict)

    def test_json_contains_expected_keys(self, basic_evidence):
        """The JSON output contains required top-level keys."""
        json_str = export_json(
            query="test",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        data = json.loads(json_str)
        assert "report_type" in data
        assert "version" in data
        assert "query" in data
        assert "response" in data
        assert "evidence" in data
        assert "search_metrics" in data

    def test_json_query_value(self, basic_evidence):
        """The JSON 'query' field matches the input."""
        json_str = export_json(
            query="my specific query",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        data = json.loads(json_str)
        assert data["query"] == "my specific query"

    def test_json_is_not_comparative_by_default(self, basic_evidence):
        """Standard queries have is_comparative == False."""
        json_str = export_json(
            query="test",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        data = json.loads(json_str)
        assert data["is_comparative"] is False

    def test_json_comparative(self, comparative_evidence):
        """Comparative results have is_comparative == True and comparative data."""
        json_str = export_json(
            query="Compare CD19 vs BCMA",
            response_text="Comparison.",
            comp_result=comparative_evidence,
        )
        data = json.loads(json_str)
        assert data["is_comparative"] is True
        assert "comparative" in data
        assert data["comparative"]["entity_a"] == "CD19"
        assert data["comparative"]["entity_b"] == "BCMA"


# ═══════════════════════════════════════════════════════════════════════
# PDF EXPORT
# ═══════════════════════════════════════════════════════════════════════


class TestExportPdf:
    """Tests for the export_pdf() function."""

    def test_returns_bytes(self, basic_evidence):
        """export_pdf() returns bytes."""
        pdf_bytes = export_pdf(
            query="CD19 efficacy",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        assert isinstance(pdf_bytes, bytes)

    def test_starts_with_pdf_magic(self, basic_evidence):
        """The PDF output starts with the %PDF magic bytes."""
        pdf_bytes = export_pdf(
            query="test",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        assert pdf_bytes[:5] == b"%PDF-"

    def test_non_empty_output(self, basic_evidence):
        """The PDF output is non-empty."""
        pdf_bytes = export_pdf(
            query="test",
            response_text="Answer.",
            evidence=basic_evidence,
        )
        assert len(pdf_bytes) > 100

    def test_handles_no_evidence(self):
        """export_pdf() works when evidence is None."""
        pdf_bytes = export_pdf(
            query="test",
            response_text="No evidence.",
            evidence=None,
        )
        assert isinstance(pdf_bytes, bytes)
        assert pdf_bytes[:5] == b"%PDF-"

    def test_handles_comparative_result(self, comparative_evidence):
        """export_pdf() renders comparative results without error."""
        pdf_bytes = export_pdf(
            query="Compare CD19 vs BCMA",
            response_text="Comparison analysis.\n\n**Key differences:**\n- CRS rates\n- Persistence",
            comp_result=comparative_evidence,
        )
        assert isinstance(pdf_bytes, bytes)
        assert pdf_bytes[:5] == b"%PDF-"

    def test_handles_markdown_tables_in_response(self, basic_evidence):
        """export_pdf() handles markdown tables in the response text."""
        response_with_table = (
            "## Summary\n\n"
            "| Product | Target | CRS Rate |\n"
            "|---------|--------|----------|\n"
            "| Kymriah | CD19 | 58% |\n"
            "| Yescarta | CD19 | 93% |\n"
        )
        pdf_bytes = export_pdf(
            query="CD19 products CRS comparison",
            response_text=response_with_table,
            evidence=basic_evidence,
        )
        assert isinstance(pdf_bytes, bytes)
        assert len(pdf_bytes) > 100


# ═══════════════════════════════════════════════════════════════════════
# FILENAME GENERATION
# ═══════════════════════════════════════════════════════════════════════


class TestGenerateFilename:
    """Tests for the generate_filename() function."""

    def test_markdown_extension(self):
        """generate_filename('md') returns a .md filename."""
        fname = generate_filename("md")
        assert fname.endswith(".md")
        assert fname.startswith("cart_query_")

    def test_json_extension(self):
        """generate_filename('json') returns a .json filename."""
        fname = generate_filename("json")
        assert fname.endswith(".json")
        assert fname.startswith("cart_query_")

    def test_pdf_extension(self):
        """generate_filename('pdf') returns a .pdf filename."""
        fname = generate_filename("pdf")
        assert fname.endswith(".pdf")

    def test_timestamp_format(self):
        """The filename contains a timestamp in ISO-8601 YYYYMMDDTHHMMSSZ format."""
        fname = generate_filename("md")
        # Expected: cart_query_20260219T143025Z.md
        parts = fname.replace("cart_query_", "").replace(".md", "")
        assert len(parts) == 16  # YYYYMMDDTHHMMSSZ
        assert parts[8] == "T"
        assert parts[-1] == "Z"

    def test_two_filenames_differ(self):
        """Two sequential filenames could be identical (same second) but are valid."""
        fname1 = generate_filename("md")
        fname2 = generate_filename("json")
        assert fname1 != fname2  # Different extensions guarantee difference
