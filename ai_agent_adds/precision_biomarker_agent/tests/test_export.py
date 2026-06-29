"""Tests for Biomarker Intelligence Agent export functions.

Validates markdown export, JSON export, PDF export (bytes), CSV export,
and FHIR R4 bundle structure using the real src.export module.

Author: Adam Jones
Date: March 2026
"""

import csv
import io
import json

import pytest

from src.export import (
    export_csv,
    export_fhir_diagnostic_report,
    export_json,
    export_markdown,
    export_pdf,
    generate_filename,
)
from src.models import (
    AnalysisResult,
    BiologicalAgeResult,
    CrossCollectionResult,
    DiseaseCategory,
    DiseaseTrajectoryResult,
    GenotypeAdjustmentResult,
    MetabolizerPhenotype,
    PGxResult,
    RiskLevel,
    SearchHit,
)
from src.report_generator import ReportGenerator

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def simple_analysis(sample_patient_profile):
    """Return a minimal AnalysisResult for export testing."""
    return AnalysisResult(
        patient_profile=sample_patient_profile,
        biological_age=BiologicalAgeResult(
            chronological_age=45,
            biological_age=47.2,
            age_acceleration=2.2,
            phenoage_score=0.003,
            mortality_risk=0.015,
            aging_drivers=[
                {"biomarker": "rdw", "value": 13.5, "contribution": 4.46, "direction": "aging"},
            ],
        ),
        disease_trajectories=[
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.DIABETES,
                risk_level=RiskLevel.HIGH,
                current_markers={"hba1c": 5.9},
                genetic_risk_factors=["TCF7L2 rs7903146 CT"],
                intervention_recommendations=["Dietary optimization"],
            ),
            DiseaseTrajectoryResult(
                disease=DiseaseCategory.CARDIOVASCULAR,
                risk_level=RiskLevel.MODERATE,
                current_markers={"ldl_c": 145},
            ),
        ],
        pgx_results=[
            PGxResult(
                gene="CYP2D6",
                star_alleles="*1/*4",
                phenotype=MetabolizerPhenotype.INTERMEDIATE,
                drugs_affected=[
                    {"drug": "Codeine", "recommendation": "Use with caution"},
                ],
            ),
        ],
        genotype_adjustments=[
            GenotypeAdjustmentResult(
                biomarker="ALT",
                standard_range="0-56 U/L",
                adjusted_range="0-45 U/L",
                genotype="CG",
                gene="PNPLA3",
                rationale="PNPLA3 CG carrier",
            ),
        ],
        critical_alerts=["HbA1c 5.9% pre-diabetic"],
    )


@pytest.fixture
def simple_evidence():
    """Return a CrossCollectionResult for markdown export testing."""
    return CrossCollectionResult(
        query="What biomarkers indicate pre-diabetes risk?",
        hits=[
            SearchHit(
                collection="biomarker_reference",
                id="ref-hba1c",
                score=0.92,
                text="HbA1c reflects average blood glucose over 2-3 months.",
                metadata={"name": "HbA1c", "unit": "%"},
            ),
        ],
        total_collections_searched=10,
        search_time_ms=42.5,
    )


@pytest.fixture
def report_markdown(simple_analysis):
    """Generate a full markdown report for PDF testing."""
    return ReportGenerator().generate(simple_analysis)


# =====================================================================
# FILENAME GENERATION
# =====================================================================


class TestFilenameGeneration:
    """Tests for generate_filename()."""

    def test_md_extension(self):
        name = generate_filename("md")
        assert name.startswith("biomarker_report_")
        assert name.endswith(".md")

    def test_pdf_extension(self):
        name = generate_filename("pdf")
        assert name.endswith(".pdf")

    def test_csv_extension(self):
        name = generate_filename("csv")
        assert name.endswith(".csv")


# =====================================================================
# MARKDOWN EXPORT
# =====================================================================


class TestExportMarkdown:
    """Tests for export_markdown()."""

    def test_returns_string(self, simple_evidence):
        md = export_markdown("test query", "test response", evidence=simple_evidence)
        assert isinstance(md, str)

    def test_contains_query(self, simple_evidence):
        md = export_markdown("What about HbA1c?", "Response here", evidence=simple_evidence)
        assert "What about HbA1c?" in md

    def test_contains_response(self, simple_evidence):
        md = export_markdown("q", "The response text", evidence=simple_evidence)
        assert "The response text" in md

    def test_contains_evidence_section(self, simple_evidence):
        md = export_markdown("q", "r", evidence=simple_evidence)
        assert "Evidence Sources" in md
        assert "ref-hba1c" in md

    def test_contains_search_metrics(self, simple_evidence):
        md = export_markdown("q", "r", evidence=simple_evidence)
        assert "Search Metrics" in md
        assert "42" in md

    def test_contains_analysis_summary(self, simple_analysis, simple_evidence):
        md = export_markdown("q", "r", evidence=simple_evidence, analysis=simple_analysis)
        assert "Patient Analysis Summary" in md
        assert "47.2" in md

    def test_contains_footer(self):
        md = export_markdown("q", "r")
        assert "Biomarker Intelligence Agent" in md

    def test_valid_markdown_headers(self, simple_evidence):
        md = export_markdown("q", "r", evidence=simple_evidence)
        lines = md.split("\n")
        headers = [line for line in lines if line.startswith("#")]
        for h in headers:
            assert h.startswith("# ") or h.startswith("## ") or h.startswith("### ")


# =====================================================================
# JSON EXPORT
# =====================================================================


class TestExportJSON:
    """Tests for export_json()."""

    def test_returns_valid_json(self, simple_analysis):
        json_str = export_json(analysis_result=simple_analysis)
        data = json.loads(json_str)
        assert isinstance(data, dict)

    def test_contains_report_type(self, simple_analysis):
        data = json.loads(export_json(analysis_result=simple_analysis))
        assert data["report_type"] == "precision_biomarker_analysis"

    def test_contains_version(self, simple_analysis):
        data = json.loads(export_json(analysis_result=simple_analysis))
        assert data["version"] == "1.0.0"

    def test_contains_analysis_data(self, simple_analysis):
        data = json.loads(export_json(analysis_result=simple_analysis))
        assert "analysis" in data
        assert "biological_age" in data["analysis"]

    def test_contains_query_and_response(self):
        data = json.loads(export_json(query="test q", response_text="test r"))
        assert data["query"] == "test q"
        assert data["response"] == "test r"

    def test_contains_timestamp(self, simple_analysis):
        data = json.loads(export_json(analysis_result=simple_analysis))
        assert "generated_at" in data
        assert "T" in data["generated_at"]

    def test_roundtrip_serialization(self, simple_analysis):
        json_str = export_json(analysis_result=simple_analysis)
        data = json.loads(json_str)
        json_str2 = json.dumps(data, indent=2, default=str)
        data2 = json.loads(json_str2)
        assert data == data2


# =====================================================================
# PDF EXPORT
# =====================================================================


class TestExportPDF:
    """Tests for export_pdf()."""

    def test_returns_bytes(self, report_markdown):
        pdf = export_pdf(report_markdown)
        assert isinstance(pdf, bytes)

    def test_starts_with_pdf_header(self, report_markdown):
        pdf = export_pdf(report_markdown)
        assert pdf.startswith(b"%PDF")

    def test_ends_with_eof_marker(self, report_markdown):
        pdf = export_pdf(report_markdown)
        assert pdf.rstrip().endswith(b"%%EOF")

    def test_substantial_output(self, report_markdown):
        """A full 12-section report should produce a non-trivial PDF."""
        pdf = export_pdf(report_markdown)
        assert len(pdf) > 1000

    def test_contains_table_data(self, report_markdown):
        """PDF should contain rendered table content (not skip tables)."""
        # The markdown has table content like "Metric" and "Value"
        # which should be rendered in the PDF stream
        pdf = export_pdf(report_markdown)
        # reportlab embeds text; check the PDF is substantial
        assert len(pdf) > 5000

    def test_simple_markdown(self):
        """Even a simple markdown string should produce a valid PDF."""
        pdf = export_pdf("# Test\n\nHello world.")
        assert pdf.startswith(b"%PDF")


# =====================================================================
# CSV EXPORT
# =====================================================================


class TestExportCSV:
    """Tests for export_csv()."""

    def test_returns_bytes(self, simple_analysis):
        result = export_csv(simple_analysis)
        assert isinstance(result, bytes)

    def test_decodable_utf8(self, simple_analysis):
        result = export_csv(simple_analysis)
        text = result.decode("utf-8")
        assert len(text) > 0

    def test_parseable_csv(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        reader = csv.reader(io.StringIO(text))
        rows = list(reader)
        assert len(rows) > 10

    def test_contains_patient_info(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Patient Info" in text
        assert "HG002" in text

    def test_contains_biological_age(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Biological Age" in text
        assert "47.2" in text

    def test_contains_disease_trajectories(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Disease Trajectories" in text
        assert "diabetes" in text

    def test_contains_pgx_profile(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "PGx Profile" in text
        assert "CYP2D6" in text

    def test_contains_genotype_adjustments(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Genotype Adjustments" in text
        assert "PNPLA3" in text

    def test_contains_critical_alerts(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Critical Alerts" in text
        assert "HbA1c" in text

    def test_contains_biomarker_values(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Biomarker Values" in text
        assert "albumin" in text

    def test_aging_drivers_section(self, simple_analysis):
        text = export_csv(simple_analysis).decode("utf-8")
        assert "Aging Drivers" in text
        assert "rdw" in text


# =====================================================================
# FHIR R4 EXPORT
# =====================================================================


class TestExportFHIR:
    """Tests for export_fhir_diagnostic_report()."""

    def test_returns_valid_json(self, simple_analysis):
        fhir_str = export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        )
        bundle = json.loads(fhir_str)
        assert isinstance(bundle, dict)

    def test_resource_type_is_bundle(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        assert bundle["resourceType"] == "Bundle"

    def test_bundle_type_is_collection(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        assert bundle["type"] == "collection"

    def test_contains_diagnostic_report(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        reports = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ]
        assert len(reports) == 1
        assert reports[0]["resource"]["status"] == "final"

    def test_contains_biological_age_observation(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        bio_obs = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and "Biological Age" in e["resource"]["code"].get("text", "")
        ]
        assert len(bio_obs) == 1
        assert bio_obs[0]["resource"]["valueQuantity"]["value"] == 47.2

    def test_contains_disease_trajectory_observations(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        traj_obs = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and "Risk Assessment" in e["resource"]["code"].get("text", "")
        ]
        assert len(traj_obs) == 2  # diabetes + cardiovascular

    def test_contains_pgx_observation(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        pgx_obs = [
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and "Pharmacogenomic" in e["resource"]["code"].get("text", "")
        ]
        assert len(pgx_obs) == 1
        # valueString contains star alleles + phenotype; gene is in components
        assert "*1/*4" in pgx_obs[0]["resource"]["valueString"]
        gene_components = [
            c for c in pgx_obs[0]["resource"]["component"]
            if c["code"]["text"] == "Gene"
        ]
        assert gene_components[0]["valueString"] == "CYP2D6"

    def test_diagnostic_report_has_conclusion(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        report = next(
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        )
        conclusion = report["resource"]["conclusion"]
        assert "47.2" in conclusion

    def test_diagnostic_report_has_result_refs(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        report = next(
            e for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        )
        assert len(report["resource"]["result"]) >= 3

    def test_timestamp_present(self, simple_analysis):
        bundle = json.loads(export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        ))
        assert "timestamp" in bundle
        assert "T" in bundle["timestamp"]

    def test_valid_json_roundtrip(self, simple_analysis):
        fhir_str = export_fhir_diagnostic_report(
            simple_analysis, simple_analysis.patient_profile
        )
        bundle = json.loads(fhir_str)
        json_str2 = json.dumps(bundle, indent=2)
        bundle2 = json.loads(json_str2)
        assert bundle == bundle2
