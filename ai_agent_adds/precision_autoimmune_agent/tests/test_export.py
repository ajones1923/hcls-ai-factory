"""
Tests for the Precision Autoimmune Agent export module.

Covers:
- Markdown report generation (full and empty results)
- PDF byte generation (with reportlab fallback)
- FHIR R4 bundle creation and JSON serialization

All tests run offline — no external services required.
"""

import json

import pytest

from src.export import AutoimmuneExporter
from src.models import (
    AutoimmuneAnalysisResult,
    AutoimmuneDisease,
    BiologicTherapy,
    DiseaseActivityLevel,
    DiseaseActivityScore,
    FlarePredictor,
    FlareRisk,
)

# ── Fixtures ──────────────────────────────────────────────────────────


@pytest.fixture
def exporter():
    return AutoimmuneExporter(knowledge_version={"version": "2.0.0", "last_updated": "2026-03-10"})


@pytest.fixture
def full_result():
    """Fully populated AutoimmuneAnalysisResult."""
    return AutoimmuneAnalysisResult(
        patient_id="PAT-001",
        disease_activity_scores=[
            DiseaseActivityScore(
                disease=AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
                score_name="DAS28-CRP",
                score_value=4.8,
                level=DiseaseActivityLevel.HIGH,
                components={"tender_joints": 8, "swollen_joints": 5, "CRP": 22.0},
                thresholds={"remission": 2.6, "low": 3.2, "moderate": 5.1},
            ),
        ],
        flare_predictions=[
            FlarePredictor(
                disease=AutoimmuneDisease.RHEUMATOID_ARTHRITIS,
                current_activity=DiseaseActivityLevel.HIGH,
                predicted_risk=FlareRisk.HIGH,
                risk_score=0.82,
                contributing_factors=["Elevated CRP", "Recent medication change"],
                protective_factors=["Adherence to methotrexate"],
                recommended_monitoring=["CRP every 2 weeks", "DAS28 monthly"],
            ),
        ],
        hla_associations=[
            {"allele": "HLA-DRB1*04:01", "disease": "rheumatoid_arthritis", "odds_ratio": 4.2, "pmid": "20301572"},
        ],
        biologic_recommendations=[
            BiologicTherapy(
                drug_name="Adalimumab",
                drug_class="TNF inhibitor",
                mechanism="Binds soluble and membrane-bound TNF-alpha",
                indicated_diseases=[AutoimmuneDisease.RHEUMATOID_ARTHRITIS],
                pgx_considerations=["HLA-DQA1*05 associated with anti-drug antibody formation"],
                contraindications=["Active TB", "Severe heart failure (NYHA III-IV)"],
                monitoring_requirements=["TB screening", "CBC every 3 months", "LFTs quarterly"],
            ),
        ],
        critical_alerts=["Elevated CRP with high disease activity — consider biologic escalation"],
    )


@pytest.fixture
def minimal_result():
    """Minimal AutoimmuneAnalysisResult with empty lists."""
    return AutoimmuneAnalysisResult(patient_id="PAT-EMPTY")


# ── Markdown tests ────────────────────────────────────────────────────


class TestToMarkdown:

    def test_returns_string(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert isinstance(md, str)

    def test_contains_patient_id(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "PAT-001" in md

    def test_contains_disease_activity_section(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "Disease Activity Scores" in md
        assert "DAS28-CRP" in md
        assert "4.8" in md
        assert "HIGH" in md

    def test_contains_flare_predictions_section(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "Flare Risk Predictions" in md
        assert "0.82" in md
        assert "Contributing Factors" in md
        assert "Elevated CRP" in md

    def test_contains_hla_associations_section(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "HLA-Disease Associations" in md
        assert "HLA-DRB1*04:01" in md

    def test_contains_biologic_recommendations_section(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "Biologic Therapy Recommendations" in md
        assert "Adalimumab" in md
        assert "TNF inhibitor" in md
        assert "PGx Considerations" in md
        assert "Contraindications" in md
        assert "Monitoring Requirements" in md

    def test_contains_critical_alerts_section(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "Critical Alerts" in md
        assert "biologic escalation" in md

    def test_contains_footer(self, exporter, full_result):
        md = exporter.to_markdown("PAT-001", analysis_result=full_result)
        assert "clinical validation" in md
        assert "2.0.0" in md

    def test_empty_result_no_sections(self, exporter, minimal_result):
        md = exporter.to_markdown("PAT-EMPTY", analysis_result=minimal_result)
        assert "PAT-EMPTY" in md
        # Should not have body sections since all lists are empty
        assert "Disease Activity Scores" not in md
        assert "Flare Risk Predictions" not in md
        assert "Critical Alerts" not in md
        assert "Biologic Therapy Recommendations" not in md

    def test_no_analysis_result(self, exporter):
        md = exporter.to_markdown("PAT-NONE")
        assert "PAT-NONE" in md
        assert isinstance(md, str)

    def test_with_query_answer(self, exporter):
        md = exporter.to_markdown("PAT-Q", query_answer="TNF inhibitors are recommended.")
        assert "Clinical Query Response" in md
        assert "TNF inhibitors are recommended." in md

    def test_with_evidence_hits(self, exporter):
        hits = [
            {"collection": "clinical_docs", "score": 0.92, "text": "Evidence text here", "relevance": "high"},
        ]
        md = exporter.to_markdown("PAT-E", evidence_hits=hits)
        assert "Evidence Sources" in md
        assert "clinical_docs" in md
        assert "0.920" in md

    def test_knowledge_version_defaults(self):
        exp = AutoimmuneExporter()
        md = exp.to_markdown("PAT-X")
        assert "N/A" in md or "?" in md


# ── PDF tests ─────────────────────────────────────────────────────────


class TestToPdfBytes:

    def test_returns_bytes(self, exporter, full_result):
        result = exporter.to_pdf_bytes("PAT-001", analysis_result=full_result)
        assert isinstance(result, bytes)
        assert len(result) > 0

    def test_empty_result_returns_bytes(self, exporter, minimal_result):
        result = exporter.to_pdf_bytes("PAT-EMPTY", analysis_result=minimal_result)
        assert isinstance(result, bytes)

    def test_no_analysis_returns_bytes(self, exporter):
        result = exporter.to_pdf_bytes("PAT-NONE")
        assert isinstance(result, bytes)

    def test_with_query_answer(self, exporter):
        result = exporter.to_pdf_bytes("PAT-Q", query_answer="Clinical response text.\n\nSecond paragraph.")
        assert isinstance(result, bytes)
        assert len(result) > 0

    def test_fallback_when_no_reportlab(self, exporter, full_result, monkeypatch):
        """When reportlab is not importable, should fall back to markdown bytes."""
        import builtins
        real_import = builtins.__import__

        def mock_import(name, *args, **kwargs):
            if "reportlab" in name:
                raise ImportError("mocked")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", mock_import)
        result = exporter.to_pdf_bytes("PAT-001", analysis_result=full_result)
        assert isinstance(result, bytes)
        # Fallback returns UTF-8 encoded markdown, which should contain the patient ID
        decoded = result.decode("utf-8")
        assert "PAT-001" in decoded


# ── FHIR R4 tests ────────────────────────────────────────────────────


class TestToFhirR4:

    def test_returns_dict(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        assert isinstance(bundle, dict)

    def test_resource_type_is_bundle(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        assert bundle["resourceType"] == "Bundle"

    def test_bundle_type(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        assert bundle["type"] == "collection"

    def test_contains_patient_entry(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        patient_entries = [
            e for e in bundle["entry"] if e["resource"]["resourceType"] == "Patient"
        ]
        assert len(patient_entries) == 1
        assert patient_entries[0]["resource"]["id"] == "PAT-001"

    def test_contains_diagnostic_report(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        reports = [
            e for e in bundle["entry"] if e["resource"]["resourceType"] == "DiagnosticReport"
        ]
        assert len(reports) == 1
        report = reports[0]["resource"]
        assert report["status"] == "final"
        assert report["subject"]["reference"] == "Patient/PAT-001"

    def test_contains_observation_entries(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        observations = [
            e for e in bundle["entry"] if e["resource"]["resourceType"] == "Observation"
        ]
        # 1 disease activity score + 1 flare prediction = 2 observations
        assert len(observations) == 2

    def test_disease_activity_observation(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        obs = [
            e["resource"] for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation" and "activity" in e["resource"]["id"]
        ]
        assert len(obs) == 1
        assert obs[0]["valueQuantity"]["value"] == 4.8

    def test_flare_risk_observation(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        obs = [
            e["resource"] for e in bundle["entry"]
            if e["resource"]["resourceType"] == "Observation" and "flare" in e["resource"]["id"]
        ]
        assert len(obs) == 1
        assert obs[0]["valueQuantity"]["value"] == 0.82

    def test_conclusion_includes_critical_alerts(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        report = [
            e["resource"] for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ][0]
        assert "ALERT" in report["conclusion"]

    def test_minimal_result_has_patient_and_report(self, exporter, minimal_result):
        bundle = exporter.to_fhir_r4("PAT-EMPTY", analysis_result=minimal_result)
        resource_types = [e["resource"]["resourceType"] for e in bundle["entry"]]
        assert "Patient" in resource_types
        assert "DiagnosticReport" in resource_types

    def test_no_analysis_result(self, exporter):
        bundle = exporter.to_fhir_r4("PAT-NONE")
        assert bundle["resourceType"] == "Bundle"
        report = [
            e["resource"] for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ][0]
        assert report["conclusion"] == "Analysis complete. See observations."

    def test_with_query_answer(self, exporter):
        bundle = exporter.to_fhir_r4("PAT-Q", query_answer="The patient should start TNF therapy.")
        report = [
            e["resource"] for e in bundle["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ][0]
        assert "TNF therapy" in report["conclusion"]

    def test_timestamp_present(self, exporter, full_result):
        bundle = exporter.to_fhir_r4("PAT-001", analysis_result=full_result)
        assert "timestamp" in bundle
        assert len(bundle["timestamp"]) > 0


class TestToFhirJson:

    def test_returns_string(self, exporter, full_result):
        result = exporter.to_fhir_json("PAT-001", analysis_result=full_result)
        assert isinstance(result, str)

    def test_valid_json(self, exporter, full_result):
        result = exporter.to_fhir_json("PAT-001", analysis_result=full_result)
        parsed = json.loads(result)
        assert parsed["resourceType"] == "Bundle"

    def test_minimal_result_valid_json(self, exporter, minimal_result):
        result = exporter.to_fhir_json("PAT-EMPTY", analysis_result=minimal_result)
        parsed = json.loads(result)
        assert parsed["resourceType"] == "Bundle"
        assert "entry" in parsed
