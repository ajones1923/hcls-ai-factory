"""Tests for Clinical Guardrails safety module.

All tests run in lite mode (no NeMo Guardrails dependency required).
Validates PII detection/sanitization, clinical disclaimer injection,
evidence grounding checks, absolute diagnosis flagging, input rejection,
lite-mode fallback, and disabled-mode passthrough.

Author: Adam Jones
Date: April 2026
"""

import sys
from pathlib import Path
from unittest.mock import patch

import pytest

# Ensure project root on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.safety.guardrails_config import (
    CLINICAL_DISCLAIMER,
    DISCLAIMER_SENTINEL,
    ClinicalGuardrails,
    GuardrailResult,
)


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def guardrails():
    """Return a ClinicalGuardrails instance in lite mode (default)."""
    return ClinicalGuardrails(enabled=True, mode="lite")


@pytest.fixture
def disabled_guardrails():
    """Return a ClinicalGuardrails instance with safety checks disabled."""
    return ClinicalGuardrails(enabled=False)


# ===================================================================
# PII DETECTION
# ===================================================================


class TestPIIDetection:
    """Tests for PII pattern detection in input text."""

    def test_ssn_detected(self, guardrails):
        result = guardrails.check_input("Patient SSN is 123-45-6789")
        assert not result.passed
        assert result.blocked
        assert any("ssn" in f for f in result.flags)

    def test_mrn_detected(self, guardrails):
        result = guardrails.check_input("MRN: 12345678 shows hemorrhage")
        assert not result.passed
        assert result.blocked
        assert any("mrn" in f for f in result.flags)

    def test_mrn_no_colon_detected(self, guardrails):
        result = guardrails.check_input("MRN 9876543 imaging study")
        assert not result.passed
        assert any("mrn" in f for f in result.flags)

    def test_email_detected(self, guardrails):
        result = guardrails.check_input(
            "Send results to patient@hospital.org"
        )
        assert not result.passed
        assert any("email" in f for f in result.flags)

    def test_phone_detected(self, guardrails):
        result = guardrails.check_input("Call patient at (555) 123-4567")
        assert not result.passed
        assert any("phone" in f for f in result.flags)

    def test_multiple_pii_types(self, guardrails):
        text = "SSN 123-45-6789, email doctor@clinic.com, MRN: 99887766"
        result = guardrails.check_input(text)
        assert not result.passed
        assert len(result.flags) >= 3

    def test_no_pii_passes(self, guardrails):
        result = guardrails.check_input(
            "What are the Lung-RADS criteria for solid nodules?"
        )
        assert result.passed
        assert not result.blocked


# ===================================================================
# PII SANITIZATION
# ===================================================================


class TestPIISanitization:
    """Tests for PII replacement with [REDACTED]."""

    def test_ssn_redacted(self, guardrails):
        text = "Patient SSN is 123-45-6789"
        sanitized = guardrails.sanitize_pii(text)
        assert "123-45-6789" not in sanitized
        assert "[REDACTED]" in sanitized

    def test_mrn_redacted(self, guardrails):
        sanitized = guardrails.sanitize_pii("MRN: 12345678")
        assert "12345678" not in sanitized
        assert "[REDACTED]" in sanitized

    def test_email_redacted(self, guardrails):
        sanitized = guardrails.sanitize_pii("email: doc@hospital.com")
        assert "doc@hospital.com" not in sanitized
        assert "[REDACTED]" in sanitized

    def test_phone_redacted(self, guardrails):
        sanitized = guardrails.sanitize_pii("Phone: (555) 123-4567")
        assert "(555) 123-4567" not in sanitized
        assert "[REDACTED]" in sanitized

    def test_clean_text_unchanged(self, guardrails):
        text = "CT head shows normal anatomy"
        sanitized = guardrails.sanitize_pii(text)
        assert sanitized == text

    def test_multiple_pii_all_redacted(self, guardrails):
        text = "SSN 123-45-6789 and email test@example.com"
        sanitized = guardrails.sanitize_pii(text)
        assert "123-45-6789" not in sanitized
        assert "test@example.com" not in sanitized
        assert sanitized.count("[REDACTED]") == 2


# ===================================================================
# CLINICAL DISCLAIMER
# ===================================================================


class TestClinicalDisclaimer:
    """Tests for clinical disclaimer injection."""

    def test_disclaimer_added_to_bare_response(self, guardrails):
        response = "CT shows a 12mm solid nodule in the right upper lobe."
        result = guardrails.check_output(response, evidence_count=3)
        assert result.modified_text is not None
        assert DISCLAIMER_SENTINEL in result.modified_text

    def test_disclaimer_not_duplicated(self, guardrails):
        response = (
            "CT shows a 12mm nodule.\n\n"
            "Clinical Decision Support Notice: For informational purposes."
        )
        result = guardrails.check_output(response, evidence_count=3)
        # Disclaimer sentinel already present, so modified_text should
        # only be set if other changes were made (not for disclaimer)
        if result.modified_text:
            assert result.modified_text.count(DISCLAIMER_SENTINEL) == 1

    def test_add_disclaimer_method_idempotent(self, guardrails):
        response = "Normal chest CT."
        with_disclaimer = guardrails.add_disclaimer(response)
        assert DISCLAIMER_SENTINEL in with_disclaimer
        # Calling again should not double-append
        double = guardrails.add_disclaimer(with_disclaimer)
        assert double == with_disclaimer


# ===================================================================
# EVIDENCE GROUNDING
# ===================================================================


class TestEvidenceGrounding:
    """Tests for evidence citation grounding checks."""

    def test_zero_evidence_flagged(self, guardrails):
        result = guardrails.check_output(
            "The CT shows a hemorrhage.", evidence_count=0
        )
        assert "no_evidence_citations" in result.flags

    def test_positive_evidence_not_flagged(self, guardrails):
        result = guardrails.check_output(
            "The CT shows findings consistent with hemorrhage [PMID:123].",
            evidence_count=5,
        )
        assert "no_evidence_citations" not in result.flags


# ===================================================================
# ABSOLUTE DIAGNOSIS FLAGGING
# ===================================================================


class TestAbsoluteDiagnosis:
    """Tests for absolute diagnostic claim detection."""

    def test_this_is_cancer_flagged(self, guardrails):
        result = guardrails.check_output(
            "This is cancer based on imaging features.",
            evidence_count=3,
        )
        assert "absolute_diagnosis_claim" in result.flags
        assert not result.passed

    def test_definitely_flagged(self, guardrails):
        result = guardrails.check_output(
            "This is definitely a malignant lesion requiring surgery.",
            evidence_count=3,
        )
        assert "absolute_diagnosis_claim" in result.flags
        assert not result.passed

    def test_confirm_diagnosis_flagged(self, guardrails):
        result = guardrails.check_output(
            "I can confirm the diagnosis of adenocarcinoma.",
            evidence_count=3,
        )
        assert "absolute_diagnosis_claim" in result.flags
        assert not result.passed

    def test_hedged_language_passes(self, guardrails):
        result = guardrails.check_output(
            "The findings are suggestive of malignancy and warrant "
            "further evaluation with tissue biopsy.",
            evidence_count=3,
        )
        assert "absolute_diagnosis_claim" not in result.flags


# ===================================================================
# SAFE RESPONSE PASSES
# ===================================================================


class TestSafeResponse:
    """Tests that well-formed clinical responses pass all checks."""

    def test_safe_response_passes(self, guardrails):
        response = (
            "Based on the available evidence, the CT findings are suggestive "
            "of a pulmonary nodule measuring 8mm in the right upper lobe. "
            "Per ACR Lung-RADS v2022, this corresponds to category 4A. "
            "[Evidence: PMID:38001234] "
            "Clinical Decision Support Notice: This analysis is for "
            "informational purposes only."
        )
        result = guardrails.check_output(response, evidence_count=5)
        assert result.passed
        assert "absolute_diagnosis_claim" not in result.flags
        assert "contraindicated_recommendation" not in result.flags


# ===================================================================
# INPUT REJECTION
# ===================================================================


class TestInputRejection:
    """Tests for non-medical and unsafe input query rejection."""

    def test_non_medical_query_rejected(self, guardrails):
        result = guardrails.check_input("What is the stock price of NVDA?")
        assert not result.passed
        assert result.blocked
        assert "non_medical_query" in result.flags

    def test_recipe_query_rejected(self, guardrails):
        result = guardrails.check_input("What is the recipe for chocolate cake?")
        assert not result.passed
        assert "non_medical_query" in result.flags

    def test_prescribe_rejected(self, guardrails):
        result = guardrails.check_input(
            "Prescribe amoxicillin for this patient"
        )
        assert not result.passed
        assert "unsafe_request" in result.flags

    def test_stop_medication_rejected(self, guardrails):
        result = guardrails.check_input(
            "Should I stop taking medication based on this scan?"
        )
        assert not result.passed
        assert "unsafe_request" in result.flags

    def test_diagnose_definitively_rejected(self, guardrails):
        result = guardrails.check_input(
            "Diagnose definitively what disease this is"
        )
        assert not result.passed
        assert "unsafe_request" in result.flags

    def test_medical_query_accepted(self, guardrails):
        result = guardrails.check_input(
            "What are the BI-RADS criteria for breast lesion classification?"
        )
        assert result.passed
        assert not result.blocked

    def test_imaging_query_accepted(self, guardrails):
        result = guardrails.check_input(
            "Compare CT vs MRI for detecting brain hemorrhage"
        )
        assert result.passed


# ===================================================================
# LITE MODE FALLBACK
# ===================================================================


class TestLiteModeFallback:
    """Tests that lite mode works when nemoguardrails is not importable."""

    def test_lite_mode_when_import_fails(self):
        """Simulate nemoguardrails not installed — should fall back to lite."""
        with patch.dict("sys.modules", {"nemoguardrails": None}):
            gr = ClinicalGuardrails(enabled=True, mode="full")
            # Should have fallen back to lite mode
            assert gr.mode == "lite"
            assert gr._rails is None

        # Lite mode checks still work
        result = gr.check_input("What is the stock price today?")
        assert not result.passed

    def test_lite_mode_pii_still_works(self):
        """PII detection works in lite mode without NeMo."""
        gr = ClinicalGuardrails(enabled=True, mode="lite")
        result = gr.check_input("Patient SSN 123-45-6789")
        assert not result.passed
        assert any("ssn" in f for f in result.flags)

    def test_lite_mode_output_checks_work(self):
        """Output safety checks work in lite mode."""
        gr = ClinicalGuardrails(enabled=True, mode="lite")
        result = gr.check_output(
            "This is definitely cancer.",
            evidence_count=0,
        )
        assert "absolute_diagnosis_claim" in result.flags
        assert "no_evidence_citations" in result.flags
        assert not result.passed


# ===================================================================
# DISABLED MODE
# ===================================================================


class TestDisabledMode:
    """Tests that all checks pass when guardrails are disabled."""

    def test_disabled_input_always_passes(self, disabled_guardrails):
        result = disabled_guardrails.check_input(
            "Patient SSN 123-45-6789 stock price prescribe"
        )
        assert result.passed
        assert not result.blocked
        assert result.flags == []

    def test_disabled_output_always_passes(self, disabled_guardrails):
        result = disabled_guardrails.check_output(
            "This is definitely cancer. Take 500mg immediately.",
            evidence_count=0,
        )
        assert result.passed
        assert result.flags == []
        assert result.modified_text is None


# ===================================================================
# GUARDRAIL RESULT MODEL
# ===================================================================


class TestGuardrailResult:
    """Tests for the GuardrailResult Pydantic model."""

    def test_default_values(self):
        result = GuardrailResult(passed=True)
        assert result.passed is True
        assert result.flags == []
        assert result.modified_text is None
        assert result.blocked is False

    def test_with_flags(self):
        result = GuardrailResult(
            passed=False,
            flags=["pii_detected:ssn", "non_medical_query"],
            blocked=True,
        )
        assert not result.passed
        assert len(result.flags) == 2
        assert result.blocked

    def test_serialization(self):
        result = GuardrailResult(
            passed=False,
            flags=["test_flag"],
            modified_text="sanitized text",
        )
        data = result.model_dump()
        assert data["passed"] is False
        assert data["flags"] == ["test_flag"]
        assert data["modified_text"] == "sanitized text"


# ===================================================================
# CONTRAINDICATED OUTPUT PATTERNS
# ===================================================================


class TestContraindicatedOutput:
    """Tests for contraindicated recommendation detection in output."""

    def test_prescribing_flagged(self, guardrails):
        result = guardrails.check_output(
            "I am prescribing metformin 500mg.",
            evidence_count=1,
        )
        assert "contraindicated_recommendation" in result.flags
        assert not result.passed

    def test_stop_medication_flagged(self, guardrails):
        result = guardrails.check_output(
            "You should stop your current chemotherapy regimen.",
            evidence_count=2,
        )
        assert "contraindicated_recommendation" in result.flags
        assert not result.passed

    def test_dosage_recommendation_flagged(self, guardrails):
        result = guardrails.check_output(
            "Based on imaging, take 200 mg of ibuprofen daily.",
            evidence_count=3,
        )
        assert "contraindicated_recommendation" in result.flags
        assert not result.passed

    def test_clinical_recommendation_passes(self, guardrails):
        result = guardrails.check_output(
            "Clinical correlation and follow-up imaging in 3 months "
            "is recommended per ACR guidelines.",
            evidence_count=4,
        )
        assert "contraindicated_recommendation" not in result.flags
