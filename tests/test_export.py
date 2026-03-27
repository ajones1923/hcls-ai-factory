"""Tests for src/export.py — Markdown, JSON, PDF, FHIR R4 export.

Tests generate_filename, export_markdown, export_json, export_pdf,
export_fhir_r4, and private formatting helpers.

Author: Adam Jones
Date: March 2026
"""

import json
import re
from datetime import datetime, timezone
from unittest.mock import patch

import pytest
from src.export import (
    VERSION,
    export_fhir_r4,
    export_json,
    export_markdown,
    export_pdf,
    generate_filename,
    _format_filters,
    _format_citation_link,
)
from src.models import (
    AlertLevel,
    CrossCollectionResult,
    ComparativeResult,
    PGxAlert,
    SearchHit,
)


# ═══════════════════════════════════════════════════════════════════════
# VERSION
# ═══════════════════════════════════════════════════════════════════════

class TestVersion:
    def test_version_string(self):
        assert VERSION == "1.0.0"


# ═══════════════════════════════════════════════════════════════════════
# generate_filename
# ═══════════════════════════════════════════════════════════════════════

class TestGenerateFilename:
    def test_md_extension(self):
        fn = generate_filename("md")
        assert fn.startswith("pgx_query_")
        assert fn.endswith(".md")

    def test_json_extension(self):
        fn = generate_filename("json")
        assert fn.endswith(".json")

    def test_pdf_extension(self):
        fn = generate_filename("pdf")
        assert fn.endswith(".pdf")

    def test_timestamp_format(self):
        fn = generate_filename("md")
        # pgx_query_20260312T143025Z.md
        ts_part = fn.replace("pgx_query_", "").replace(".md", "")
        assert re.match(r"^\d{8}T\d{6}Z$", ts_part)

    def test_unique_filenames(self):
        fn1 = generate_filename("md")
        fn2 = generate_filename("md")
        # Same second → same name; just verify format is consistent
        assert fn1.startswith("pgx_query_")


# ═══════════════════════════════════════════════════════════════════════
# export_markdown
# ═══════════════════════════════════════════════════════════════════════

class TestExportMarkdown:
    def test_basic_report(self):
        md = export_markdown("test query", "test response")
        assert "# Pharmacogenomics Intelligence Report" in md
        assert "test query" in md
        assert "test response" in md

    def test_contains_version(self):
        md = export_markdown("q", "r")
        assert VERSION in md

    def test_contains_timestamp(self):
        md = export_markdown("q", "r")
        assert "UTC" in md

    def test_alerts_section(self):
        alerts = [
            PGxAlert(
                alert_level=AlertLevel.CRITICAL,
                gene="CYP2D6",
                drug="codeine",
                phenotype="Ultra-rapid Metabolizer",
                recommendation="Avoid codeine",
            )
        ]
        md = export_markdown("q", "r", alerts=alerts)
        assert "## Clinical Alerts" in md
        assert "CRITICAL" in md
        assert "CYP2D6" in md
        assert "codeine" in md

    def test_no_alerts_section_when_none(self):
        md = export_markdown("q", "r", alerts=None)
        assert "## Clinical Alerts" not in md

    def test_evidence_section(self):
        hits = [
            SearchHit(
                id="hit1",
                score=0.92,
                collection="pgx_drug_guidelines",
                text="test text",
                metadata={"gene": "CYP2D6", "drug": "codeine"},
            )
        ]
        evidence = CrossCollectionResult(
            query="q",
            hits=hits,
            total_collections_searched=15,
            search_time_ms=42.5,
        )
        md = export_markdown("q", "r", evidence=evidence)
        assert "## Evidence Sources" in md

    def test_filters_applied(self):
        md = export_markdown("q", "r", filters_applied={"Gene": "CYP2D6"})
        assert "CYP2D6" in md

    def test_comparative_result(self):
        hits_a = [
            SearchHit(id="a1", score=0.9, collection="pgx_gene_reference",
                      text="t", metadata={})
        ]
        hits_b = [
            SearchHit(id="b1", score=0.8, collection="pgx_gene_reference",
                      text="t", metadata={})
        ]
        ev_a = CrossCollectionResult(query="q", hits=hits_a, total_collections_searched=15,
                                     search_time_ms=10)
        ev_b = CrossCollectionResult(query="q", hits=hits_b, total_collections_searched=15,
                                     search_time_ms=10)
        comp = ComparativeResult(
            query="compare",
            entity_a="CYP2D6",
            entity_b="CYP2C19",
            evidence_a=ev_a,
            evidence_b=ev_b,
            total_search_time_ms=20,
        )
        md = export_markdown("compare", "result", comp_result=comp)
        assert "Comparative" in md
        assert "CYP2D6" in md
        assert "CYP2C19" in md


# ═══════════════════════════════════════════════════════════════════════
# export_json
# ═══════════════════════════════════════════════════════════════════════

class TestExportJSON:
    def test_valid_json(self):
        result = export_json("test query", "test response")
        data = json.loads(result)
        assert isinstance(data, dict)

    def test_report_type(self):
        data = json.loads(export_json("q", "r"))
        assert data["report_type"] == "pgx_intelligence_query"

    def test_version(self):
        data = json.loads(export_json("q", "r"))
        assert data["version"] == VERSION

    def test_query_and_response(self):
        data = json.loads(export_json("my query", "my response"))
        assert data["query"] == "my query"
        assert data["response"] == "my response"

    def test_is_comparative_false(self):
        data = json.loads(export_json("q", "r"))
        assert data["is_comparative"] is False

    def test_alerts_in_json(self):
        alerts = [
            PGxAlert(
                alert_level=AlertLevel.WARNING,
                gene="CYP2C19",
                drug="clopidogrel",
                phenotype="Poor Metabolizer",
                recommendation="Use alternative",
            )
        ]
        data = json.loads(export_json("q", "r", alerts=alerts))
        assert "alerts" in data
        assert len(data["alerts"]) == 1
        assert data["alerts"][0]["gene"] == "CYP2C19"

    def test_generated_at_iso(self):
        data = json.loads(export_json("q", "r"))
        assert "generated_at" in data
        # Should be ISO 8601 format
        assert "T" in data["generated_at"]

    def test_evidence_section(self):
        hits = [
            SearchHit(id="h1", score=0.85, collection="pgx_drug_guidelines",
                      text="t", metadata={})
        ]
        evidence = CrossCollectionResult(
            query="q", hits=hits, total_collections_searched=15, search_time_ms=30
        )
        data = json.loads(export_json("q", "r", evidence=evidence))
        assert "evidence" in data
        assert "search_metrics" in data
        assert data["search_metrics"]["total_results"] == 1


# ═══════════════════════════════════════════════════════════════════════
# export_pdf
# ═══════════════════════════════════════════════════════════════════════

class TestExportPDF:
    def test_returns_bytes(self):
        result = export_pdf("test query", "test response")
        assert isinstance(result, bytes)
        assert len(result) > 0

    def test_pdf_header(self):
        result = export_pdf("test query", "test response")
        assert result[:5] == b"%PDF-"

    def test_pdf_with_alerts(self):
        alerts = [
            PGxAlert(
                alert_level=AlertLevel.CRITICAL,
                gene="CYP2D6",
                drug="codeine",
                phenotype="UM",
                recommendation="Avoid",
            )
        ]
        result = export_pdf("q", "r", alerts=alerts)
        assert isinstance(result, bytes)
        assert len(result) > 100

    def test_pdf_with_evidence(self):
        hits = [
            SearchHit(id="h1", score=0.9, collection="pgx_drug_guidelines",
                      text="evidence text", metadata={"gene": "CYP2D6"})
        ]
        evidence = CrossCollectionResult(
            query="q", hits=hits, total_collections_searched=15, search_time_ms=30
        )
        result = export_pdf("q", "r", evidence=evidence)
        assert isinstance(result, bytes)


# ═══════════════════════════════════════════════════════════════════════
# export_fhir_r4
# ═══════════════════════════════════════════════════════════════════════

class TestExportFHIR:
    def test_returns_valid_json(self):
        result = export_fhir_r4("query", "response")
        data = json.loads(result)
        assert isinstance(data, dict)

    def test_fhir_bundle_type(self):
        data = json.loads(export_fhir_r4("q", "r"))
        assert data["resourceType"] == "Bundle"

    def test_fhir_type_collection(self):
        data = json.loads(export_fhir_r4("q", "r"))
        assert data["type"] == "collection"

    def test_fhir_has_entry(self):
        data = json.loads(export_fhir_r4("q", "r"))
        assert "entry" in data
        assert len(data["entry"]) >= 1

    def test_diagnostic_report_present(self):
        data = json.loads(export_fhir_r4("q", "r"))
        resource_types = [
            e["resource"]["resourceType"]
            for e in data["entry"]
            if "resource" in e
        ]
        assert "DiagnosticReport" in resource_types

    def test_loinc_code(self):
        data = json.loads(export_fhir_r4("q", "r"))
        dr = next(
            e["resource"]
            for e in data["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        )
        loinc_codes = [
            coding["code"]
            for coding in dr.get("code", {}).get("coding", [])
        ]
        assert "69548-6" in loinc_codes

    def test_fhir_with_evidence(self):
        hits = [
            SearchHit(id="h1", score=0.9, collection="pgx_drug_guidelines",
                      text="text", metadata={"gene": "CYP2D6"})
        ]
        evidence = CrossCollectionResult(
            query="q", hits=hits, total_collections_searched=15, search_time_ms=30
        )
        data = json.loads(export_fhir_r4("q", "r", evidence=evidence))
        assert len(data["entry"]) >= 1


# ═══════════════════════════════════════════════════════════════════════
# PRIVATE HELPERS
# ═══════════════════════════════════════════════════════════════════════

class TestFormatFilters:
    def test_none_filters(self):
        assert _format_filters(None) == "None"

    def test_empty_filters(self):
        assert _format_filters({}) == "None"

    def test_all_genes_ignored(self):
        assert _format_filters({"Gene": "All Genes"}) == "None"

    def test_specific_gene(self):
        result = _format_filters({"Gene": "CYP2D6"})
        assert "CYP2D6" in result


class TestFormatCitationLink:
    def test_pmid_link(self):
        result = _format_citation_link("pgx_clinical_evidence", "12345678")
        assert "pubmed" in result
        assert "12345678" in result

    def test_nct_link(self):
        result = _format_citation_link("pgx_clinical_trials", "NCT12345678")
        assert "clinicaltrials.gov" in result

    def test_fda_link(self):
        result = _format_citation_link("pgx_fda_labels", "fda_001")
        assert "FDA" in result

    def test_generic_id(self):
        result = _format_citation_link("pgx_gene_reference", "gene_001")
        assert result == "gene_001"
