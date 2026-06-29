"""
Tests for the export module.
==============================
Validates Markdown, JSON, and FHIR R4 export functions for
oncology reports and MTB packets.
"""

import json
import sys
from pathlib import Path

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.export import (
    EVIDENCE_LEVEL_LABELS,
    FHIR_LOINC_CODES,
    FHIR_SNOMED_CANCER_CODES,
    export_fhir_r4,
    export_json,
    export_markdown,
)


# ═══════════════════════════════════════════════════════════════════════════
# Test Data Fixtures
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def full_mtb_packet():
    """A fully populated MTB packet dict for testing exports."""
    return {
        "patient_id": "PT-001",
        "cancer_type": "nsclc",
        "sample_id": "S-001",
        "title": "Oncology Intelligence Report — PT-001",
        "summary": "Stage IV NSCLC with EGFR L858R. High PD-L1 expression.",
        "variants": [
            {
                "gene": "EGFR",
                "variant_name": "L858R",
                "variant_type": "missense",
                "vaf": 0.35,
                "consequence": "missense_variant",
                "tier": "I",
            },
            {
                "gene": "TP53",
                "variant_name": "R175H",
                "variant_type": "missense",
                "vaf": 0.42,
                "consequence": "missense_variant",
                "tier": "III",
            },
        ],
        "biomarkers": {
            "tmb": 12.5,
            "msi": "MSS",
            "pdl1": 80,
        },
        "evidence": [
            {
                "gene": "EGFR",
                "evidence_level": "level_1",
                "source": "PMID:33096080",
                "summary": "EGFR L858R responds to osimertinib.",
            },
        ],
        "therapies": [
            {
                "name": "osimertinib",
                "targets": ["EGFR"],
                "evidence_level": "A",
                "line_of_therapy": "1L",
                "notes": "FLAURA trial",
            },
            {
                "name": "pembrolizumab",
                "targets": ["PD-1"],
                "evidence_level": "A",
                "line_of_therapy": "1L",
                "notes": "KEYNOTE-024 (PD-L1 >= 50%)",
            },
        ],
        "clinical_trials": [
            {
                "nct_id": "NCT04613596",
                "title": "ADAURA: Adjuvant Osimertinib",
                "phase": "Phase 3",
                "status": "Active",
                "match_rationale": "EGFR L858R match",
            },
        ],
        "resistance_mechanisms": [
            {
                "mechanism": "T790M",
                "drug": "erlotinib",
                "description": "Gatekeeper mutation conferring resistance to 1st-gen EGFR TKIs.",
            },
        ],
        "open_questions": [
            "TP53 R175H is a variant of uncertain significance.",
            "Consider liquid biopsy for ctDNA monitoring.",
        ],
    }


@pytest.fixture
def minimal_mtb_packet():
    """A minimal MTB packet for testing edge cases."""
    return {
        "patient_id": "PT-EMPTY",
        "cancer_type": "breast",
    }


@pytest.fixture
def empty_mtb_packet():
    """An empty dict MTB packet."""
    return {}


# ═══════════════════════════════════════════════════════════════════════════
# Markdown Export Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestExportMarkdown:
    """Test export_markdown returns valid markdown string."""

    def test_returns_string(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert isinstance(result, str)

    def test_contains_title(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "Oncology Intelligence Report" in result

    def test_contains_patient_id(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "PT-001" in result

    def test_contains_cancer_type(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "nsclc" in result.lower() or "NSCLC" in result

    def test_contains_variant_table(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "EGFR" in result
        assert "L858R" in result

    def test_contains_therapy_ranking(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "osimertinib" in result
        assert "Therapy Ranking" in result

    def test_contains_clinical_trials(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "NCT04613596" in result
        assert "Clinical Trial" in result

    def test_contains_biomarkers(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "12.5" in result  # TMB value

    def test_contains_open_questions(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "Open Questions" in result
        assert "TP53" in result

    def test_contains_disclaimer(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "research use only" in result.lower()

    def test_contains_resistance(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet)
        assert "T790M" in result

    def test_custom_title(self, full_mtb_packet):
        result = export_markdown(full_mtb_packet, title="Custom Report Title")
        assert "Custom Report Title" in result

    def test_minimal_packet(self, minimal_mtb_packet):
        """Minimal packet should produce valid markdown without errors."""
        result = export_markdown(minimal_mtb_packet)
        assert isinstance(result, str)
        assert "PT-EMPTY" in result

    def test_empty_packet(self, empty_mtb_packet):
        """Empty packet should not crash."""
        result = export_markdown(empty_mtb_packet)
        assert isinstance(result, str)

    def test_string_input(self):
        """String input should produce minimal report."""
        result = export_markdown("This is a plain text report.")
        assert isinstance(result, str)
        assert "plain text report" in result

    def test_markdown_has_headers(self, full_mtb_packet):
        """Output should use Markdown headers (# or ##)."""
        result = export_markdown(full_mtb_packet)
        assert "#" in result

    def test_variant_table_format(self, full_mtb_packet):
        """Variant table should use Markdown table format."""
        result = export_markdown(full_mtb_packet)
        assert "|" in result  # Markdown table separator


# ═══════════════════════════════════════════════════════════════════════════
# JSON Export Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestExportJSON:
    """Test export_json returns valid JSON-serialisable dict."""

    def test_returns_dict(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert isinstance(result, dict)

    def test_is_json_serializable(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        serialized = json.dumps(result)
        assert isinstance(serialized, str)

    def test_has_meta_section(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert "meta" in result
        assert result["meta"]["format"] == "hcls-ai-factory-oncology-report"

    def test_has_patient_id(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert result["patient_id"] == "PT-001"

    def test_has_cancer_type(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert result["cancer_type"] == "nsclc"

    def test_has_variants(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert "variants" in result
        assert len(result["variants"]) == 2

    def test_has_therapy_ranking(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert "therapy_ranking" in result
        assert len(result["therapy_ranking"]) == 2

    def test_has_clinical_trials(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert "clinical_trials" in result

    def test_has_open_questions(self, full_mtb_packet):
        result = export_json(full_mtb_packet)
        assert "open_questions" in result

    def test_minimal_packet(self, minimal_mtb_packet):
        """Minimal packet should produce valid JSON."""
        result = export_json(minimal_mtb_packet)
        assert isinstance(result, dict)
        assert result["patient_id"] == "PT-EMPTY"

    def test_empty_packet(self, empty_mtb_packet):
        """Empty packet should not crash."""
        result = export_json(empty_mtb_packet)
        assert isinstance(result, dict)
        assert "meta" in result

    def test_none_values_removed(self, empty_mtb_packet):
        """None values at top level should be removed."""
        result = export_json(empty_mtb_packet)
        for key, value in result.items():
            assert value is not None


# ═══════════════════════════════════════════════════════════════════════════
# FHIR R4 Export Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestExportFHIR:
    """Test export_fhir_r4 returns valid FHIR bundle structure."""

    def test_returns_dict(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert isinstance(result, dict)

    def test_bundle_resource_type(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert result["resourceType"] == "Bundle"

    def test_bundle_type_collection(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert result["type"] == "collection"

    def test_has_bundle_id(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert "id" in result
        assert len(result["id"]) > 0

    def test_has_timestamp(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert "timestamp" in result

    def test_has_entries(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert "entry" in result
        assert len(result["entry"]) > 0

    def test_has_patient_resource(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        patient_entries = [
            e for e in result["entry"]
            if e["resource"]["resourceType"] == "Patient"
        ]
        assert len(patient_entries) == 1

    def test_patient_identifier(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        patient = next(
            e["resource"] for e in result["entry"]
            if e["resource"]["resourceType"] == "Patient"
        )
        assert patient["identifier"][0]["value"] == "PT-001"

    def test_has_diagnostic_report(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        report_entries = [
            e for e in result["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        ]
        assert len(report_entries) == 1

    def test_has_variant_observations(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        obs_entries = [
            e for e in result["entry"]
            if e["resource"]["resourceType"] == "Observation"
        ]
        # Should have at least 2 variant observations + biomarker observations
        assert len(obs_entries) >= 2

    def test_observations_have_loinc_coding(self, full_mtb_packet):
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        for entry in result["entry"]:
            if entry["resource"]["resourceType"] == "Observation":
                code = entry["resource"]["code"]
                coding = code["coding"][0]
                assert coding["system"] == "http://loinc.org"

    def test_has_tmb_observation(self, full_mtb_packet):
        """TMB biomarker should generate an Observation."""
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        tmb_obs = [
            e for e in result["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and e["resource"]["code"]["coding"][0]["code"] == FHIR_LOINC_CODES["tumor_mutation_burden"]
        ]
        assert len(tmb_obs) == 1

    def test_has_msi_observation(self, full_mtb_packet):
        """MSI biomarker should generate an Observation."""
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        msi_obs = [
            e for e in result["entry"]
            if e["resource"]["resourceType"] == "Observation"
            and e["resource"]["code"]["coding"][0]["code"] == FHIR_LOINC_CODES["microsatellite_instability"]
        ]
        assert len(msi_obs) == 1

    def test_nsclc_snomed_coding(self, full_mtb_packet):
        """NSCLC should get a SNOMED CT code in the DiagnosticReport."""
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        report = next(
            e["resource"] for e in result["entry"]
            if e["resource"]["resourceType"] == "DiagnosticReport"
        )
        assert "conclusionCode" in report
        coding = report["conclusionCode"][0]["coding"][0]
        assert coding["system"] == "http://snomed.info/sct"
        expected_code = FHIR_SNOMED_CANCER_CODES["nsclc"][0]
        assert coding["code"] == expected_code

    def test_minimal_packet(self, minimal_mtb_packet):
        """Minimal packet should produce valid FHIR bundle."""
        result = export_fhir_r4(minimal_mtb_packet, patient_id="PT-EMPTY")
        assert result["resourceType"] == "Bundle"
        assert len(result["entry"]) >= 1  # At least Patient resource

    def test_empty_packet(self, empty_mtb_packet):
        """Empty packet should produce valid FHIR bundle."""
        result = export_fhir_r4(empty_mtb_packet, patient_id="PT-NONE")
        assert result["resourceType"] == "Bundle"

    def test_meta_profile(self, full_mtb_packet):
        """Bundle should reference genomics-reporting profile."""
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        assert "meta" in result
        profiles = result["meta"]["profile"]
        assert any("genomics-reporting" in p for p in profiles)

    def test_entry_full_urls(self, full_mtb_packet):
        """Each entry should have a fullUrl in urn:uuid format."""
        result = export_fhir_r4(full_mtb_packet, patient_id="PT-001")
        for entry in result["entry"]:
            assert entry["fullUrl"].startswith("urn:uuid:")


# ═══════════════════════════════════════════════════════════════════════════
# Constants Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestExportConstants:
    """Test export constants are correctly defined."""

    def test_evidence_level_labels_has_entries(self):
        assert len(EVIDENCE_LEVEL_LABELS) >= 4

    def test_fhir_loinc_codes_has_genomic_report(self):
        assert "genomic_report" in FHIR_LOINC_CODES

    def test_fhir_loinc_codes_has_variant(self):
        assert "variant" in FHIR_LOINC_CODES

    def test_fhir_snomed_has_nsclc(self):
        assert "nsclc" in FHIR_SNOMED_CANCER_CODES

    def test_fhir_snomed_has_breast(self):
        assert "breast" in FHIR_SNOMED_CANCER_CODES

    def test_fhir_snomed_has_melanoma(self):
        assert "melanoma" in FHIR_SNOMED_CANCER_CODES
