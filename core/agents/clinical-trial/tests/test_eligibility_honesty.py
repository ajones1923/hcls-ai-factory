"""Tests that the agent does not claim eligibility work it did not do.

Covers two paths:

* ``TrialRAGEngine.search_eligibility`` -- the LLM path, which derived a
  documented ``eligible (bool)`` from whether the literal string "NOT MET"
  appeared in the first 200 characters of the model's prose (F6).
* ``/v1/trial/match`` -- the only matching path reachable over HTTP, which
  reported a vector-similarity score as "Match Score" and "Confidence" with no
  criterion evaluated, and invented NCT IDs, phase and status (F1).
"""

import pytest

from src.rag_engine import TrialRAGEngine


# ═══════════════════════════════════════════════════════════════════════════
# F6 — search_eligibility must not read a boolean out of prose
# ═══════════════════════════════════════════════════════════════════════════


def _assess(answer):
    """Run the structured-verdict parser over a model response."""
    assessment = {
        "trial_id": "NCT02511106",
        "eligible": None,
        "matched_criteria": [],
        "unmatched_criteria": [],
        "uncertain_criteria": [],
        "assessment": "",
    }
    TrialRAGEngine._apply_structured_eligibility(assessment, answer)
    return assessment


def test_failure_past_character_200_is_not_eligible():
    """The exact prose that measured ``eligible: True`` pre-fix.

    ``assessment["eligible"] = "NOT MET" not in answer.upper()[:200]``. A
    criterion-by-criterion answer puts its first failure well past character
    200, so the model's correct determination was discarded.

    Clinical consequence: a consumer that trusts the boolean and not the prose
    enrols a patient with active brain metastases into a trial that excludes
    them -- a protocol deviation with a plausible CNS-progression safety
    consequence and a reportable GCP violation.
    """
    prose = (
        "Criterion 1: MET - the patient is 62 years old, which is within the "
        "protocol-specified range of 18 to 75 years of age. "
        "Criterion 2: MET - histology is confirmed as non-small cell lung "
        "cancer on the submitted pathology report. "
        "Criterion 3: MET - an EGFR L858R activating mutation is present on "
        "the molecular profile. "
        "Criterion 4: NOT MET - the patient has active brain metastases, "
        "which this protocol lists as an exclusion criterion."
    )
    # The pre-fix boolean was `"NOT MET" not in answer.upper()[:200]`, so the
    # determination had to land inside the first 200 characters to be seen.
    assert prose.upper().find("NOT MET") > 200, "fixture must reproduce the bug"
    assert "NOT MET" not in prose.upper()[:200]
    # Unstructured prose yields no verdict at all -- never a pass.
    assert _assess(prose)["eligible"] is None


PARAPHRASED_FAILURES = [
    "The patient does not meet criterion 3.",
    "This patient is excluded by the brain metastasis criterion.",
    "Patient fails criterion 3.",
    "The patient is ineligible for this study.",
]


@pytest.mark.parametrize("prose", PARAPHRASED_FAILURES)
def test_paraphrased_failure_is_never_reported_eligible(prose):
    """A paraphrase never matched the literal token, so it read as eligible."""
    assert _assess(prose)["eligible"] is not True


def test_all_met_is_eligible():
    """Structured verdicts, all met."""
    answer = """{"criteria": [
        {"criterion": "Age 18-75", "verdict": "met", "rationale": "62"},
        {"criterion": "NSCLC", "verdict": "met", "rationale": "confirmed"}
    ], "summary": "Eligible."}"""
    result = _assess(answer)
    assert result["eligible"] is True
    assert len(result["matched_criteria"]) == 2


def test_one_not_met_is_not_eligible():
    """A single not_met verdict decides the whole assessment."""
    answer = """{"criteria": [
        {"criterion": "Age 18-75", "verdict": "met", "rationale": "62"},
        {"criterion": "No brain mets", "verdict": "not_met",
         "rationale": "active brain metastases"}
    ], "summary": "Not eligible."}"""
    result = _assess(answer)
    assert result["eligible"] is False
    assert len(result["unmatched_criteria"]) == 1


def test_any_uncertain_blocks_an_eligible_verdict():
    """Uncertain is not a pass: eligibility becomes None, not True."""
    answer = """{"criteria": [
        {"criterion": "Age 18-75", "verdict": "met", "rationale": "62"},
        {"criterion": "ECOG 0-1", "verdict": "uncertain",
         "rationale": "performance status not supplied"}
    ], "summary": "Cannot determine."}"""
    result = _assess(answer)
    assert result["eligible"] is None
    assert len(result["uncertain_criteria"]) == 1


def test_json_in_a_fenced_block_is_parsed():
    """Tolerate the model wrapping its JSON in prose or a code fence."""
    answer = (
        "Here is my assessment:\n```json\n"
        '{"criteria": [{"criterion": "Age", "verdict": "met", "rationale": "ok"}],'
        ' "summary": "Eligible."}\n```\nLet me know if you need more.'
    )
    assert _assess(answer)["eligible"] is True


def test_unparseable_response_says_so_in_the_assessment_text():
    """A consumer reading only the prose is told eligibility was not determined."""
    result = _assess("I am not able to assess this patient.")
    assert result["eligible"] is None
    assert "not determined" in result["assessment"].lower()


def test_empty_response_is_not_eligible():
    """No answer at all is not a pass."""
    assert _assess("")["eligible"] is None


def test_malformed_rows_are_not_eligible():
    """A criteria list of the wrong shape yields no verdict."""
    answer = '{"criteria": ["Age 18-75", "NSCLC"], "summary": "Eligible."}'
    assert _assess(answer)["eligible"] is not True


# ═══════════════════════════════════════════════════════════════════════════
# F1 — /v1/trial/match reports retrieval, and says so
# ═══════════════════════════════════════════════════════════════════════════


def test_trial_match_model_has_no_match_score_field():
    """The similarity score is not named "overall_score" or "confidence".

    The route emits one ``engine.search()`` similarity as ``overall_score``
    with ``confidence = round(score * 0.9, 3)``, which the UI rendered as
    ``st.metric("Match Score", ...)`` and ``Confidence: ...``, coloured green
    above 0.7. Measured: a 9-year-old boy with medulloblastoma, on rifampin,
    with interstitial lung disease, scored **0.91 / confidence 0.819** against
    ADAURA -- a closed adult trial of adjuvant osimertinib in resected
    EGFR-mutant NSCLC -- with ``inclusion_met 0/0`` hidden by a
    ``if inc_total > 0 or exc_total > 0`` guard.
    """
    from api.routes.trial_clinical import TrialMatch

    fields = set(TrialMatch.model_fields)
    assert "overall_score" not in fields
    assert "confidence" not in fields
    assert "retrieval_similarity" in fields
    assert "criteria_evaluated" in fields


def test_trial_match_defaults_do_not_invent_registry_values():
    """phase, status and trial_id default to None, not to real clinical values.

    Pre-fix: ``trial_id=md.get("nct_id") or f"NCT-RESULT-{i+1:04d}"``,
    ``phase=md.get("phase", "phase_ii")``, ``status=md.get("status",
    "recruiting")``. A literature review article with no metadata was presented
    as a recruiting Phase II trial under the identifier NCT-RESULT-0002.

    The ClinicalTrials.gov API v2 Status enum includes UNKNOWN as a real
    returned value, so defaulting an absent status to recruiting is never
    justified.
    """
    from api.routes.trial_clinical import TrialMatch

    m = TrialMatch(trial_title="A review article about EGFR inhibitors")
    assert m.trial_id is None
    assert m.phase is None
    assert m.status is None
    assert m.retrieval_similarity is None
    assert m.criteria_evaluated is False


def test_match_response_carries_a_research_use_disclaimer():
    """The machine-readable response states it is not an eligibility finding.

    The agent's only disclaimer lived in ``src/export.py``, which the HTTP
    routes do not use, so the formats designed to be consumed by a clinical
    system without a human reading them carried no qualifier at all.
    """
    from api.routes.trial_clinical import PatientMatchResponse

    r = PatientMatchResponse(
        matches=[], total_screened=0, patient_summary={}, search_criteria={},
    )
    assert r.criteria_evaluated is False
    assert "RESEARCH USE ONLY" in r.disclaimer
    assert "not an eligibility determination" in r.disclaimer
    assert "qualified investigator" in r.disclaimer


def test_match_route_docstring_does_not_claim_eligibility_screening():
    """The docstring said "Patient-trial matching with eligibility screening".

    That is what made the similarity score readable as a match score.
    """
    from api.routes.trial_clinical import patient_match

    doc = patient_match.__doc__ or ""
    assert "does NOT screen eligibility" in doc
    assert "eligibility screening" not in doc.split("does NOT")[0]
