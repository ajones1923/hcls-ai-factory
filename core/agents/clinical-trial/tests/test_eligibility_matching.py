"""Behavioural reference tests for patient-trial eligibility matching.

``PatientMatchingWorkflow._evaluate_criterion`` is the only place in the agent
where a patient is compared to an eligibility criterion, and it had **no direct
test at all** -- ``grep -rn _evaluate_criterion tests/`` returned nothing.
Every defect F3-F10 lived in those 110 lines and all of them were invisible to
a green 781-test suite.

Each test asserts the **verdict, the polarity and the evidence**, in the shape
the audit prescribed: a described patient and a described criterion in, the
arithmetically or clinically correct answer out.

Primary source for the structural claim that an unparsed criterion is UNKNOWN
rather than met: ClinicalTrials.gov API v2 ``GET /studies/metadata`` --
``protocolSection.eligibilityModule.eligibilityCriteria`` is typed ``markup``
(free text), with ``minimumAge``/``maximumAge`` (``NormalizedTime``), ``sex``
(``Sex``) and ``stdAges`` (``StandardAge[]``) as the only structured gates. The
registry publishes no machine-readable per-criterion representation, so a
criterion no rule parses has no verdict.
"""

import pytest

from src.clinical_workflows import PatientMatchingWorkflow
from src.models import CriterionType


@pytest.fixture
def wf():
    return PatientMatchingWorkflow()


def _ev(wf, criterion, patient, is_exclusion=False):
    return wf._evaluate_criterion(criterion, patient, is_exclusion=is_exclusion)


# ═══════════════════════════════════════════════════════════════════════════
# F3 — Unknown must not be met
# ═══════════════════════════════════════════════════════════════════════════


UNPARSEABLE_CRITERIA = [
    "Measurable disease per RECIST v1.1",
    "ECOG performance status 0-1",
    "Platelet count >= 100,000/uL",
    "Left ventricular ejection fraction >= 50% by echocardiogram",
    "Willing and able to comply with study procedures",
    "Adequate hepatic function",
]


@pytest.mark.parametrize("criterion", UNPARSEABLE_CRITERIA)
def test_unparsed_inclusion_criterion_is_unknown_not_met(wf, criterion):
    """A criterion no rule understands has verdict None, not True.

    ``_evaluate_criterion`` started at ``met = True`` with every rule block
    conditional on the patient field being present, so a criterion no rule
    recognised stayed True and was scored as satisfied. A trial whose three
    criteria no rule parses (LVEF, compliance, RECIST) measured **inclusion
    3/3 met, score 0.60**.
    """
    score = _ev(wf, criterion, {})
    assert score.met is None
    assert score.criterion_type is CriterionType.INCLUSION


def test_empty_patient_satisfies_no_inclusion_criterion(wf):
    """An empty patient profile meets nothing.

    Measured pre-fix: ``{}`` against a 5-inclusion/5-exclusion osimertinib
    trial scored **inclusion 5/5 met, score 0.60, severity MODERATE,
    confidence 0.58** -- identical to a correctly-matching 61-year-old
    EGFR-mutant NSCLC patient.

    Clinical consequence: a coordinator pastes a partial referral (name and
    diagnosis only -- no age, no labs, no ECOG) into the matcher for an
    84-year-old with ECOG 3. Every inclusion criterion on an 18-75, ECOG 0-1
    trial is reported as met and the patient is sent for a screening visit they
    cannot pass, consuming a screening slot and exposing them to baseline
    biopsy and contrast CT for a trial they were never eligible for.
    """
    criteria = [
        "Age 18-75 years",
        "Histologically confirmed non-small cell lung cancer",
        "Documented EGFR exon 19 deletion or L858R mutation",
        "ECOG performance status 0-1",
        "Measurable disease per RECIST v1.1",
    ]
    for c in criteria:
        assert _ev(wf, c, {}).met is not True, c


def test_evidence_names_the_absent_fields(wf):
    """An unevaluated criterion says so, rather than reporting silently."""
    score = _ev(wf, "Measurable disease per RECIST v1.1", {})
    assert "not evaluated" in score.evidence.lower()
    assert score.confidence == 0.0


# ═══════════════════════════════════════════════════════════════════════════
# F4 — Exclusion polarity
# ═══════════════════════════════════════════════════════════════════════════


def test_unparsed_exclusion_criterion_is_unknown_not_triggered(wf):
    """For an exclusion, the default must not read as TRIGGERED.

    The code reused ``met`` with the opposite meaning for exclusions -- True
    means triggered -- but kept the ``True`` default, so
    ``exc_clear = sum(1 for s in exc_scores if not s.met)`` was 0 for every
    criterion no rule fired on.
    """
    score = _ev(wf, "Known hypersensitivity to study drug", {}, is_exclusion=True)
    assert score.met is None
    assert score.criterion_type is CriterionType.EXCLUSION


def test_patient_with_no_comorbidities_clears_exclusions(wf):
    """A clean patient CLEARS a comorbidity exclusion, so exc_clear > 0."""
    patient = {"comorbidities": ["well-controlled hypertension"], "medications": []}
    score = _ev(
        wf, "Active brain metastases", patient, is_exclusion=True,
    )
    assert score.met is not True


def test_matching_patient_clears_most_exclusions(wf):
    """The clinically plausible patient is not reported as tripping everything.

    Measured pre-fix: a 61-year-old EGFR-mutant NSCLC patient got
    ``exclusion clear: 0/5``, so the report asserted the patient tripped all
    five exclusion criteria. ``exc_ratio`` was therefore ~always 0 and
    ``overall = 0.6*inc + 0.4*exc`` capped at 0.60, making the
    ``overall_score >= 0.7`` "top match" branch and the ``top_score >= 0.8``
    LOW-severity branch dead code -- every patient got "No high-confidence
    matches found -- consider broadening search criteria or reviewing waiver
    options".

    Clinical consequence: a patient who genuinely qualifies is reported as "0
    of 5 exclusion criteria clear" with a recommendation to review waiver
    options. The coordinator reads five triggered exclusions, concludes a
    protocol waiver is needed for each, and drops the referral.
    """
    patient = {
        "age": 61, "sex": "female",
        "diagnosis": "non-small cell lung cancer",
        "biomarkers": ["EGFR"], "comorbidities": [], "medications": [],
    }
    exclusions = [
        "Active brain metastases",
        "Pregnancy or breastfeeding",
        "Prior treatment with an EGFR tyrosine kinase inhibitor",
        "Concurrent strong CYP3A4 inducers such as rifampin",
        "Interstitial lung disease",
    ]
    cleared = sum(
        1 for c in exclusions
        if _ev(wf, c, patient, is_exclusion=True).met is False
    )
    triggered = sum(
        1 for c in exclusions
        if _ev(wf, c, patient, is_exclusion=True).met is True
    )
    assert triggered == 0
    assert cleared >= 1


# ═══════════════════════════════════════════════════════════════════════════
# F5 — Conjunction, not last-writer-wins
# ═══════════════════════════════════════════════════════════════════════════


def test_age_failure_is_not_erased_by_a_later_biomarker_hit(wf):
    """A 9-year-old does not meet an adults-only criterion.

    The five rule blocks assigned to one shared ``met`` in sequence with no
    conjunction, so the biomarker block overwrote the age block's correct
    ``False``. Measured: age 9 with biomarker EGFR against "Adults aged 18
    years or older with EGFR mutation" returned **met=True, confidence 0.85**,
    evidence "Biomarker 'egfr' found in criterion" -- the audit trail no longer
    mentioned age at all, so the override was invisible.

    Clinical consequence: the child is scored as meeting the adult eligibility
    criterion of an adult-only Phase 3 at *higher* confidence (0.85) than the
    correct age determination would have carried (0.90 for met=False,
    discarded). Referred for screening, the paediatric case risks a
    consent/IRB violation and, if the site trusts the tool, exposure of a minor
    to an adult dosing regimen.
    """
    score = _ev(
        wf, "Adults aged 18 years or older with EGFR mutation",
        {"age": 9, "biomarkers": ["EGFR"]},
    )
    assert score.met is False
    # The evidence must name BOTH dimensions, so the conjunction is auditable.
    assert "9" in score.evidence
    assert "egfr" in score.evidence.lower()


def test_sex_failure_is_not_erased_by_a_later_biomarker_hit(wf):
    """A male BRCA1 carrier does not meet a female-only criterion.

    Measured pre-fix: met=True, confidence 0.85.
    """
    score = _ev(
        wf, "Females only; must carry a BRCA1 mutation",
        {"sex": "male", "biomarkers": ["BRCA1"]},
    )
    assert score.met is False
    assert "male" in score.evidence.lower()


def test_all_satisfied_dimensions_still_yield_met(wf):
    """Conjunction does not break the case where everything is satisfied."""
    score = _ev(
        wf, "Adults aged 18 years or older with EGFR mutation",
        {"age": 61, "biomarkers": ["EGFR"]},
    )
    assert score.met is True


# ═══════════════════════════════════════════════════════════════════════════
# F7 — Diagnosis specificity
# ═══════════════════════════════════════════════════════════════════════════


def test_nsclc_does_not_meet_an_sclc_only_criterion(wf):
    """Non-small cell lung cancer must not satisfy a small cell criterion.

    Matching was ``any(word in crit_lower for word in diagnosis.split() if
    len(word) > 3)`` -- any single shared word over three characters. The word
    "cancer" alone carried it, so this measured **met=True, conf 0.70** with
    evidence asserting the diagnosis "matches criterion keywords".

    Clinical consequence: SCLC and NSCLC have disjoint standard-of-care and
    disjoint trial populations. If the referral goes through and the site's own
    screening is light, the patient receives a platinum-etoposide-based or
    DLL3-targeted investigational regimen developed for a tumour they do not
    have -- no expected benefit, full toxicity exposure, and a contaminated
    trial dataset.
    """
    score = _ev(
        wf, "Histologically confirmed small cell lung cancer",
        {"diagnosis": "non-small cell lung cancer"},
    )
    assert score.met is not True


def test_breast_cancer_does_not_meet_a_colorectal_criterion(wf):
    """Measured pre-fix: met=True, conf 0.70 on the shared word "cancer"."""
    score = _ev(
        wf, "Histologically confirmed metastatic colorectal cancer",
        {"diagnosis": "breast cancer"},
    )
    assert score.met is not True


def test_matching_diagnosis_still_meets_the_criterion(wf):
    """The correct patient still matches."""
    score = _ev(
        wf, "Histologically confirmed non-small cell lung cancer",
        {"diagnosis": "non-small cell lung cancer"},
    )
    assert score.met is True


def test_unrelated_diagnosis_is_not_met(wf):
    """No shared vocabulary at all is a clean False."""
    score = _ev(
        wf, "Histologically confirmed non-small cell lung cancer",
        {"diagnosis": "medulloblastoma"},
    )
    assert score.met is False


# ═══════════════════════════════════════════════════════════════════════════
# F8 — Biomarker polarity
# ═══════════════════════════════════════════════════════════════════════════


NEGATED_BIOMARKER_CRITERIA = [
    "EGFR wild-type tumors only",
    "Patients must NOT harbor an EGFR activating mutation",
    "EGFR-negative disease",
    "Absence of EGFR mutation",
    "Confirmed disease without EGFR alteration",
]


@pytest.mark.parametrize("criterion", NEGATED_BIOMARKER_CRITERIA)
def test_egfr_mutant_patient_does_not_meet_wild_type_criterion(wf, criterion):
    """Mere MENTION of a biomarker is not satisfaction of the criterion.

    ``if bm in crit_lower: met = True`` scored an EGFR-mutant patient as
    meeting "EGFR wild-type tumors only" at 0.85 confidence.

    Clinical consequence: a patient with an EGFR exon 19 deletion is enrolled
    on a trial restricted to EGFR-wild-type disease -- a regimen
    contraindicated by their driver mutation -- while being withheld from the
    osimertinib they should receive. Direct survival harm and a confounded
    primary endpoint.
    """
    score = _ev(wf, criterion, {"biomarkers": ["EGFR"]})
    assert score.met is not True


def test_required_biomarker_present_is_met(wf):
    """A positively-stated biomarker requirement the patient satisfies."""
    score = _ev(
        wf, "Documented EGFR exon 19 deletion or L858R mutation",
        {"biomarkers": ["EGFR"]},
    )
    assert score.met is True


# ═══════════════════════════════════════════════════════════════════════════
# F9 — Age parsing, table-driven
# ═══════════════════════════════════════════════════════════════════════════

AGE_PHRASINGS = [
    "Age 18-75 years",
    "Age 18 to 75",
    "Subjects must be 18 years of age or older and no more than 75 years of age",
    "18 Years to 75 Years",
]

#: age -> expected verdict for an 18-75 criterion.
AGE_EXPECTATIONS = {8: False, 17: False, 18: True, 61: True, 75: True, 76: False, 92: False}


@pytest.mark.parametrize("criterion", AGE_PHRASINGS)
@pytest.mark.parametrize("age,expected", sorted(AGE_EXPECTATIONS.items()))
def test_age_range_is_parsed_from_every_common_phrasing(wf, criterion, age, expected):
    """An 18-75 criterion is enforced whichever way the registry phrases it.

    The regex was ``(\\d+)\\s*(?:years|yrs|y\\.o\\.)`` -- a number was only
    captured if immediately followed by a unit token -- and a single captured
    bound was only used if a comparator keyword was present. So:

      "Age 18-75 years"      -> only "75 years" captured -> single bound, no
                                comparator -> fell through to met = True
      "Age 18 to 75"         -> no unit token -> met = True
      "Patients aged >= 18"  -> no unit token -> met = True
      "Age over 75 years"    -> "over" in neither keyword list -> met = True

    all measured with patient age **92**.

    Clinical consequence: a 92-year-old is referred to a Phase 1
    dose-escalation study with an upper age limit set for a reason (renal
    clearance of the investigational agent, DLT attribution). Screening either
    wastes a slot or, if the site also trusts the tool, admits a patient 17
    years past the protocol ceiling.
    """
    score = _ev(wf, criterion, {"age": age})
    assert score.met is expected, (
        f"age {age} vs '{criterion}' -> {score.met} ({score.evidence})"
    )


LOWER_BOUND_PHRASINGS = ["Patients aged >= 18", "Age >= 18 years", "Age 18 years or older"]


@pytest.mark.parametrize("criterion", LOWER_BOUND_PHRASINGS)
@pytest.mark.parametrize("age,expected", [(8, False), (17, False), (18, True), (92, True)])
def test_lower_bound_only_criteria(wf, criterion, age, expected):
    """A minimum-age criterion is enforced with or without the unit token."""
    score = _ev(wf, criterion, {"age": age})
    assert score.met is expected, f"age {age} vs '{criterion}': {score.evidence}"


@pytest.mark.parametrize("age,expected", [(50, False), (76, True), (92, True)])
def test_over_is_recognised_as_a_lower_bound(wf, age, expected):
    """"Age over 75 years" is a minimum, and "over" was in neither keyword list."""
    score = _ev(wf, "Age over 75 years", {"age": age})
    assert score.met is expected, score.evidence


def test_explicit_range_is_not_overridden_by_a_named_age_band(wf):
    """An explicit numeric range beats the inferred "adult" band.

    ``_AGE_PATTERNS`` substituted its own bands -- adult (18,120),
    elderly/geriatric (65,120), child (0,11), adolescent (12,17) -- and was
    checked BEFORE the numeric path, with a ``break`` out of it. So "Adults
    aged 18 to 75 years" was scored as 18-120 and the upper bound was silently
    discarded.
    """
    score = _ev(wf, "Adults aged 18 to 75 years", {"age": 92})
    assert score.met is False
    assert "75" in score.evidence


def test_age_band_is_still_used_when_no_number_is_given(wf):
    """With no numbers present, a named band is a reasonable inference."""
    score = _ev(wf, "Adult patients", {"age": 8})
    assert score.met is False


def test_bare_number_with_no_comparator_is_unknown(wf):
    """A number with no direction cannot be turned into a bound."""
    score = _ev(wf, "Patients approximately 50 years", {"age": 92})
    assert score.met is None


# ═══════════════════════════════════════════════════════════════════════════
# F10 — Exclusion substring matching
# ═══════════════════════════════════════════════════════════════════════════


def test_controlled_comorbidity_does_not_trip_an_uncontrolled_exclusion(wf):
    """"controlled" is a substring of "uncontrolled", so the negation fired it.

    ``if comor in crit_lower or any(word in crit_lower for word in
    comor.split() if len(word) > 4)`` triggered "controlled hypothyroidism"
    against "Severe UNcontrolled autoimmune thyroid disease" at 0.75
    confidence, with evidence naming the controlled condition.

    Clinical consequence: a patient with well-controlled Hashimoto's
    hypothyroidism on stable levothyroxine -- explicitly permitted by the
    protocol, which excludes only *severe uncontrolled* autoimmune thyroid
    disease -- is reported as triggering that exclusion. They are screened out
    of a trial they qualify for, and the stated evidence is confident enough
    that nobody re-reads the protocol.
    """
    score = _ev(
        wf, "Severe uncontrolled autoimmune thyroid disease",
        {"comorbidities": ["controlled hypothyroidism"]}, is_exclusion=True,
    )
    assert score.met is not True


def test_hay_fever_does_not_trip_a_platinum_allergy_exclusion(wf):
    """"seasonal allergic rhinitis" vs "allergic reaction to platinum agents".

    Measured pre-fix: triggered at 0.75 confidence, so hay fever excluded the
    patient from any platinum trial.
    """
    score = _ev(
        wf, "History of allergic reaction to platinum agents",
        {"comorbidities": ["seasonal allergic rhinitis"]}, is_exclusion=True,
    )
    assert score.met is not True


def test_genuine_comorbidity_still_triggers_its_exclusion(wf):
    """The exclusion fires when the patient record matches it in full."""
    score = _ev(
        wf, "Active brain metastases",
        {"comorbidities": ["active brain metastases"]}, is_exclusion=True,
    )
    assert score.met is True
    assert "brain metastases" in score.evidence


def test_comorbidity_without_the_criterions_qualifier_is_unknown(wf):
    """"brain metastases" against "ACTIVE brain metastases" needs review.

    The two readings of a severity qualifier are not symmetric. Where the
    patient record CONTRADICTS the qualifier -- "controlled hypothyroidism"
    against "severe uncontrolled autoimmune thyroid disease" -- the criterion
    does not apply. Where the record is merely SILENT on it -- "brain
    metastases" against "active brain metastases" -- the patient may well have
    active disease and the record simply not say, so reporting the criterion
    clear would be the unsafe direction. That case is None: it neither
    excludes the patient nor clears them, and it surfaces in the
    ``exclusion_unknown`` count for clinician review.
    """
    score = _ev(
        wf, "Active brain metastases",
        {"comorbidities": ["brain metastases"]}, is_exclusion=True,
    )
    assert score.met is None
    assert "severity" in score.evidence


def test_genuine_medication_still_triggers_its_exclusion(wf):
    """A prohibited concomitant medication trips the exclusion."""
    score = _ev(
        wf, "Concurrent strong CYP3A4 inducers such as rifampin",
        {"medications": ["rifampin"]}, is_exclusion=True,
    )
    assert score.met is True


def test_medication_substring_does_not_trip_an_exclusion(wf):
    """Word-boundary matching on medication names."""
    score = _ev(
        wf, "Treatment with rifampicin within 14 days",
        {"medications": ["rifampin"]}, is_exclusion=True,
    )
    assert score.met is not True


# ═══════════════════════════════════════════════════════════════════════════
# The aggregate test: F3 and F4 together
# ═══════════════════════════════════════════════════════════════════════════

_TRIAL = {
    "trial_id": "NCT02511106",
    "title": "ADAURA: Osimertinib as Adjuvant Therapy in EGFR-Mutant NSCLC",
    "phase": "Phase III",
    "status": "Completed",
    "inclusion_criteria": [
        "Age 18-75 years",
        "Histologically confirmed non-small cell lung cancer",
        "Documented EGFR exon 19 deletion or L858R mutation",
        "ECOG performance status 0-1",
        "Measurable disease per RECIST v1.1",
    ],
    "exclusion_criteria": [
        "Active brain metastases",
        "Pregnancy or breastfeeding",
        "Prior treatment with an EGFR tyrosine kinase inhibitor",
        "Concurrent strong CYP3A4 inducers such as rifampin",
        "Interstitial lung disease",
    ],
}

_MATCHING_PATIENT = {
    "age": 61, "sex": "female",
    "diagnosis": "non-small cell lung cancer",
    "biomarkers": ["EGFR"], "comorbidities": [], "medications": [],
    "genomic_variants": [], "geographic_location": None,
}

_MISMATCHED_PATIENT = {
    "age": 8, "sex": "male", "diagnosis": "medulloblastoma",
    "biomarkers": [], "comorbidities": ["brain metastases"],
    "medications": ["rifampin"], "genomic_variants": [],
    "geographic_location": None,
}


def _run(wf, patient):
    """Run the workflow end to end and return the WorkflowResult."""
    inputs = wf.preprocess({"patient": dict(patient), "trials": [dict(_TRIAL)]})
    return wf.execute(inputs)


def _score(wf, patient, trial=None):
    """Return the structured match dict for a single trial."""
    inputs = wf.preprocess({
        "patient": dict(patient), "trials": [dict(trial or _TRIAL)],
    })
    return wf.score_trials(inputs["patient"], inputs["trials"])[0]


def test_empty_profile_scores_strictly_below_a_matching_profile(wf):
    """The headline aggregate assertion.

    Pre-fix, both profiles scored **0.60** with confidence 0.58 and
    ``exclusion clear: 0/5``, so an empty patient profile was arithmetically
    indistinguishable from a correctly matched one. This is the single test
    that would have caught F3 and F4 together.
    """
    empty_match = _score(wf, {})
    matching_match = _score(wf, _MATCHING_PATIENT)

    assert matching_match["overall_score"] is not None
    assert (
        empty_match["overall_score"] is None
        or empty_match["overall_score"] < matching_match["overall_score"]
    )


def test_empty_profile_meets_no_inclusion_criterion_in_aggregate(wf):
    """Pre-fix the empty profile measured inclusion 5/5 met."""
    match = _score(wf, {})
    assert match["inclusion_met"] == 0


def test_matching_patient_clears_all_exclusions(wf):
    """exclusion_clear == exclusion_total for a patient with nothing to trip."""
    match = _score(wf, _MATCHING_PATIENT)
    assert match["exclusion_triggered"] == 0
    assert match["exclusion_clear"] == match["exclusion_total"]


def test_mismatched_patient_is_excluded_by_a_hard_gate(wf):
    """A triggered exclusion removes the trial rather than deducting 0.4.

    Blending the ratios meant a patient who cleared 9 of 10 exclusions but
    tripped "Active brain metastases" scored 0.6*1.0 + 0.4*0.9 = **0.96** and
    became the single highest-ranked trial.
    """
    match = _score(wf, _MISMATCHED_PATIENT)
    assert match["excluded"] is True
    assert match["exclusion_reason"]


def test_unparseable_trial_gets_no_score(wf):
    """A trial whose criteria no rule parses is not scored at all.

    Measured pre-fix: a trial whose three criteria no rule parses scored a full
    **0.60** inclusion match, and a trial whose criteria failed to ingest
    entirely (``[]``) reported **confidence 0.10 with nothing evaluated**.
    """
    trial = {
        "trial_id": "NCT00000000", "title": "Unparseable",
        "inclusion_criteria": [
            "Left ventricular ejection fraction >= 50% by echocardiogram",
            "Willing and able to comply with study procedures",
            "Measurable disease per RECIST v1.1",
        ],
        "exclusion_criteria": [],
    }
    match = _score(wf, _MATCHING_PATIENT, trial)
    assert match["overall_score"] is None
    assert match["score_withheld_reason"]


def test_findings_report_unevaluated_counts(wf):
    """The output distinguishes "not met" from "not evaluated"."""
    result = _run(wf, _MATCHING_PATIENT)
    joined = " ".join(result.findings)
    assert "not evaluated" in joined


def test_findings_carry_a_research_use_statement(wf):
    """Every matching result states it is not an eligibility determination."""
    result = _run(wf, _MATCHING_PATIENT)
    joined = " ".join(result.findings)
    assert "Research use only" in joined
    assert "qualified investigator" in joined


def test_no_ungrounded_guideline_citations(wf):
    """No guideline is cited, because no guideline method is implemented.

    Every WorkflowResult used to attach a fixed ``guideline_references`` list
    ("FDA Guidance: Enrichment Strategies for Clinical Trials (2019)", "ICH
    E6(R3) GCP -- Eligibility Assessment") regardless of what was computed, and
    nothing in the code derives from either document.
    """
    result = _run(wf, _MATCHING_PATIENT)
    assert result.guideline_references == []


def test_confidence_is_none_when_nothing_was_decided(wf):
    """Confidence is not reported for an assessment that did not happen.

    ``confidence=_clamp(top_score * 0.8 + 0.1)`` reported a floor of 0.10 even
    when the trial had zero criteria.
    """
    trial = {
        "trial_id": "NCT00000001", "title": "No criteria ingested",
        "inclusion_criteria": [], "exclusion_criteria": [],
    }
    inputs = wf.preprocess({"patient": dict(_MATCHING_PATIENT), "trials": [trial]})
    result = wf.execute(inputs)
    assert result.confidence is None
