"""Reference cases for the two acute-stroke implementations, run against both.

`AcuteStrokeTriageWorkflow.execute` is the only function in this agent whose output is an
instruction to administer a thrombolytic, and `/v1/neuro/stroke/triage` is its twin on the
path the shipped Streamlit UI calls. Neither had a single test before this file: the audit's
mechanical check --
`grep -rn "Calculator|\\.run\\(|calculate\\(|_calculate_|stroke/triage" tests/` -- returned
zero matches, and none of the twelve clinical endpoints was covered.

Because the two implementations disagreed, the fixtures are shared and asserted against both
wherever the same clinical question is being asked.
"""
import asyncio
import types

import pytest

from api.routes.neuro_clinical import StrokeTriageRequest, stroke_triage
from src.clinical_scales import NIHSSCalculator
from src.clinical_workflows import AcuteStrokeTriageWorkflow


class _FakeRequest:
    """Minimal stand-in for the FastAPI Request these handlers take."""

    def __init__(self, llm=None):
        self.app = types.SimpleNamespace(
            state=types.SimpleNamespace(llm_client=llm, metrics=None, metrics_lock=None)
        )


def _triage(**kwargs):
    return asyncio.run(stroke_triage(StrokeTriageRequest(**kwargs), _FakeRequest()))


def _full_nihss(total: int) -> dict:
    """A COMPLETE 15-item NIHSS summing to *total*."""
    scores = {item: 0 for item in NIHSSCalculator.ITEMS}
    remaining = total
    for item, max_val in NIHSSCalculator.ITEMS.items():
        take = min(max_val, remaining)
        scores[item] = take
        remaining -= take
        if remaining == 0:
            break
    assert remaining == 0, "fixture cannot reach that total"
    return scores


# Every gating input present and permissive, so each test can vary exactly one thing.
_ELIGIBLE_KWARGS = dict(
    nihss_total=12,
    onset_time_hours=2.0,
    blood_pressure_systolic=150,
    blood_glucose=110.0,
    anticoagulant_use=False,
    lvo_suspected=True,
    prior_mrs=0,
    ct_aspects_score=9,
)


class TestBaselineEligibleCaseStillPasses:
    """The fixes must not simply refuse everything."""

    def test_a_clean_two_hour_case_is_eligible(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        assert response.tpa_determination == "eligible"
        assert response.tpa_eligible is True

    def test_a_clean_lvo_case_is_cleared_for_thrombectomy(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        assert response.thrombectomy_determination == "eligible"
        assert response.thrombectomy_eligible is True

    def test_the_alteplase_order_is_present(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        assert any("alteplase" in r.lower() for r in response.recommendations)


# =====================================================================
# N2 -- SBP 228 and glucose 38 left tpa_eligible True
# =====================================================================

class TestN2ExclusionsActuallyExclude:
    """`{nihss_total: 12, onset_time_hours: 2.0, blood_pressure_systolic: 228,
    blood_glucose: 38}` returned `tpa_eligible: true` whose own `tpa_reasoning` read
    "SBP 228 -- must lower to <185 before tPA". Both branches appended a sentence and left
    `tpa_eligible` untouched. SBP > 185/110 and glucose < 50 are exclusions in the AHA/ASA
    2019 guideline this endpoint names in its own `guidelines_cited`.
    """

    REPRODUCTION = dict(
        nihss_total=12, onset_time_hours=2.0,
        blood_pressure_systolic=228, blood_glucose=38.0,
    )

    def test_the_audit_reproduction_is_no_longer_tpa_eligible(self):
        response = _triage(**self.REPRODUCTION)
        assert response.tpa_eligible is False
        assert response.tpa_determination == "not_eligible"

    def test_sbp_over_185_alone_excludes(self):
        response = _triage(**{**_ELIGIBLE_KWARGS, "blood_pressure_systolic": 228})
        assert response.tpa_eligible is False
        assert any("228" in r and "185" in r for r in response.tpa_reasoning)

    def test_glucose_under_50_alone_excludes(self):
        response = _triage(**{**_ELIGIBLE_KWARGS, "blood_glucose": 38.0})
        assert response.tpa_eligible is False
        assert any("38" in r for r in response.tpa_reasoning)

    def test_glucose_over_400_alone_excludes(self):
        response = _triage(**{**_ELIGIBLE_KWARGS, "blood_glucose": 480.0})
        assert response.tpa_eligible is False

    def test_the_response_no_longer_contradicts_itself(self):
        """It said "eligible" and "must lower BP before tPA" in the same payload."""
        response = _triage(**self.REPRODUCTION)
        says_lower_bp = any("185" in r for r in response.tpa_reasoning)
        assert says_lower_bp and not response.tpa_eligible


class TestN2ThrombolyticIsNotOrderedAheadOfTheCT:
    """`recommendations.insert(0, "Administer IV alteplase ...")` put the thrombolytic above
    "Stat CT head to rule out hemorrhage", which was at index 0 of the base list. An ordered
    recommendation list is read top-down by whoever is holding the phone, and alteplase into
    an undiagnosed intracranial hemorrhage is lethal.
    """

    def test_the_ct_is_the_first_recommendation(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        assert "CT head" in response.recommendations[0]
        assert "hemorrhage" in response.recommendations[0]

    def test_the_alteplase_order_comes_after_the_ct(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        ct_index = next(
            i for i, r in enumerate(response.recommendations) if "CT head" in r
        )
        alteplase_index = next(
            i for i, r in enumerate(response.recommendations) if "alteplase" in r.lower()
        )
        assert ct_index < alteplase_index

    def test_the_alteplase_order_is_conditioned_on_the_ct(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        order = next(r for r in response.recommendations if "alteplase" in r.lower())
        assert "excludes hemorrhage" in order

    def test_the_thrombectomy_activation_also_comes_after_the_ct(self):
        response = _triage(**_ELIGIBLE_KWARGS)
        ct_index = next(
            i for i, r in enumerate(response.recommendations) if "CT head" in r
        )
        activation = next(
            i for i, r in enumerate(response.recommendations)
            if "neurointerventional" in r.lower()
        )
        assert ct_index < activation


class TestN2UnknownOnsetDoesNotClearTPA:
    """`if request.onset_time_hours is not None and ... > 4.5` fired nothing for a null
    onset. The fix that landed in `clinical_workflows.py` was never applied to this copy,
    which is the one the UI calls.
    """

    def test_a_null_onset_is_not_eligible(self):
        kwargs = dict(_ELIGIBLE_KWARGS)
        kwargs.pop("onset_time_hours")
        response = _triage(**kwargs)
        assert response.tpa_eligible is False
        assert response.tpa_determination == "indeterminate"

    def test_the_missing_input_is_named(self):
        kwargs = dict(_ELIGIBLE_KWARGS)
        kwargs.pop("onset_time_hours")
        response = _triage(**kwargs)
        assert "onset_time_hours" in response.missing_inputs
        assert any("DWI-FLAIR" in r for r in response.tpa_reasoning)

    def test_the_workflow_copy_agrees(self):
        result = AcuteStrokeTriageWorkflow().run({
            "nihss_scores": _full_nihss(12),
            "blood_pressure_systolic": 150, "blood_glucose": 110.0,
            "on_anticoagulation": False, "platelet_count": 250, "age": 70,
        })
        assert not any("tPA ELIGIBLE" in f for f in result.findings)
        assert any("tPA NOT eligible" in f for f in result.findings)


# =====================================================================
# N3 -- unmeasured contraindications passed as eligible
# =====================================================================

class TestN3UnmeasuredInputsBlockRatherThanPass:
    """The workflow defaulted `blood_pressure_systolic` to 0, `blood_glucose` to 0.0 and
    `on_anticoagulation` to False, so an unmeasured BP ran no check and emitted no warning,
    while an unmeasured glucose printed the fabricated value "glucose 0.0 mg/dL" as a finding.
    `on_anticoag` true with `inr` None emitted only "On anticoagulation -- verify INR/DOAC
    levels" and stayed eligible -- cleared for a thrombolytic on a verification that had not
    happened.
    """

    BASE = dict(
        onset_hours=2.0, nihss_scores=None, age=70,
        blood_pressure_systolic=150, blood_glucose=110.0,
        on_anticoagulation=False, platelet_count=250,
    )

    def _run(self, **overrides):
        inputs = {**self.BASE, "nihss_scores": _full_nihss(12), **overrides}
        inputs = {k: v for k, v in inputs.items() if v is not None or k in overrides}
        return AcuteStrokeTriageWorkflow().run(inputs)

    def _tpa_line(self, result):
        return next((f for f in result.findings if f.startswith("tPA ")), "")

    def test_the_baseline_case_is_still_eligible(self):
        result = self._run()
        assert "tPA ELIGIBLE" in self._tpa_line(result)

    def test_the_audit_reproduction_sbp_228_is_not_eligible(self):
        result = self._run(blood_pressure_systolic=228)
        assert "tPA NOT eligible" in self._tpa_line(result)

    def test_both_contradictory_instructions_are_no_longer_issued(self):
        result = self._run(blood_pressure_systolic=228)
        says_reduce = any("Reduce BP" in r for r in result.recommendations)
        orders_alteplase = any(
            "Administer IV alteplase" in r for r in result.recommendations
        )
        assert says_reduce and not orders_alteplase

    def test_an_unmeasured_blood_pressure_blocks(self):
        inputs = {k: v for k, v in self.BASE.items() if k != "blood_pressure_systolic"}
        inputs["nihss_scores"] = _full_nihss(12)
        result = AcuteStrokeTriageWorkflow().run(inputs)
        line = self._tpa_line(result)
        assert "tPA NOT eligible" in line
        assert "blood pressure not measured" in line

    def test_an_unmeasured_glucose_blocks_without_printing_a_value(self):
        inputs = {k: v for k, v in self.BASE.items() if k != "blood_glucose"}
        inputs["nihss_scores"] = _full_nihss(12)
        result = AcuteStrokeTriageWorkflow().run(inputs)
        line = self._tpa_line(result)
        assert "blood glucose not measured" in line
        assert "glucose 0.0 mg/dL" not in " ".join(result.findings)

    def test_an_absent_platelet_count_blocks(self):
        inputs = {k: v for k, v in self.BASE.items() if k != "platelet_count"}
        inputs["nihss_scores"] = _full_nihss(12)
        result = AcuteStrokeTriageWorkflow().run(inputs)
        assert "platelet count not available" in self._tpa_line(result)

    def test_anticoagulation_with_unknown_inr_blocks(self):
        result = self._run(on_anticoagulation=True)
        line = self._tpa_line(result)
        assert "tPA NOT eligible" in line
        assert "no INR/DOAC level available" in line

    def test_anticoagulation_with_a_normal_inr_is_still_eligible(self):
        result = self._run(on_anticoagulation=True, inr=1.1)
        assert "tPA ELIGIBLE" in self._tpa_line(result)

    def test_an_unrecognised_contraindication_is_not_discarded(self):
        """`contraindications & _TPA_CONTRAINDICATIONS` dropped `recent_gi_bleed_21d`,
        `recent_major_surgery_14d` and `infective_endocarditis` in silence."""
        result = self._run(contraindications=["recent_gi_bleed_21d"])
        line = self._tpa_line(result)
        assert "tPA NOT eligible" in line
        assert "recent_gi_bleed_21d" in line
        assert "UNRECOGNISED CONTRAINDICATION" in line

    def test_a_recognised_contraindication_still_excludes(self):
        result = self._run(contraindications=["intracranial_hemorrhage"])
        assert "intracranial_hemorrhage" in self._tpa_line(result)


# =====================================================================
# N6 -- thrombectomy cleared with no LVO; prior_mrs never read
# =====================================================================

class TestN6ThrombectomyRequiresAnOcclusionAndFunctionalIndependence:
    """Workflow entry was `if lvo or nihss >= 6` and the standard-window branch re-checked
    only `nihss >= 6`, so `{onset_hours: 3.0, NIHSS 10, age: 60}` with no CTA and no CT
    regions yielded "THROMBECTOMY ELIGIBLE (standard window <6h, NIHSS >= 6)". ASPECTS was
    permissive when absent (`aspects_result is None or score >= 6`).

    In the endpoint, `prior_mrs` was declared, validated `ge=0, le=5`, documented as
    "Pre-stroke mRS score" -- and never referenced in the function body, so
    `{nihss_total: 10, lvo_suspected: true, prior_mrs: 5}`, a patient bedridden before the
    stroke, returned `thrombectomy_eligible: true`.
    """

    def test_the_workflow_reproduction_is_no_longer_cleared(self):
        result = AcuteStrokeTriageWorkflow().run({
            "onset_hours": 3.0, "nihss_scores": _full_nihss(10), "age": 60,
        })
        assert not any("THROMBECTOMY ELIGIBLE" in f for f in result.findings)
        assert not any("Activate neurointerventional" in r for r in result.recommendations)

    def test_the_workflow_names_the_missing_occlusion(self):
        result = AcuteStrokeTriageWorkflow().run({
            "onset_hours": 3.0, "nihss_scores": _full_nihss(10), "age": 60,
        })
        blockers = next(f for f in result.findings if "THROMBECTOMY NOT CLEARED" in f)
        assert "no vessel imaging result" in blockers
        assert "pre-stroke mRS not recorded" in blockers

    def test_an_unscored_ct_is_not_a_favourable_ct(self):
        result = AcuteStrokeTriageWorkflow().run({
            "onset_hours": 3.0, "nihss_scores": _full_nihss(10), "age": 60,
            "large_vessel_occlusion": True, "prestroke_mrs": 0,
        })
        blockers = next(f for f in result.findings if "THROMBECTOMY NOT CLEARED" in f)
        assert "ASPECTS not scored" in blockers

    def test_the_workflow_clears_a_complete_standard_window_case(self):
        result = AcuteStrokeTriageWorkflow().run({
            "onset_hours": 3.0, "nihss_scores": _full_nihss(10), "age": 60,
            "large_vessel_occlusion": True, "prestroke_mrs": 0,
            "ct_affected_regions": ["M1"],
        })
        assert any("THROMBECTOMY ELIGIBLE" in f for f in result.findings)

    def test_the_endpoint_reproduction_bedridden_patient_is_not_cleared(self):
        response = _triage(nihss_total=10, lvo_suspected=True, prior_mrs=5,
                           onset_time_hours=3.0, ct_aspects_score=9)
        assert response.thrombectomy_eligible is False
        assert any("mRS 5" in r for r in response.thrombectomy_reasoning)

    def test_the_endpoint_does_not_clear_thrombectomy_with_no_vessel_imaging(self):
        kwargs = dict(_ELIGIBLE_KWARGS)
        kwargs.pop("lvo_suspected")
        response = _triage(**kwargs)
        assert response.thrombectomy_eligible is False
        assert response.thrombectomy_determination == "indeterminate"

    def test_the_endpoint_does_not_clear_thrombectomy_with_no_aspects(self):
        kwargs = dict(_ELIGIBLE_KWARGS)
        kwargs.pop("ct_aspects_score")
        response = _triage(**kwargs)
        assert response.thrombectomy_eligible is False
        assert "ct_aspects_score" in response.missing_inputs

    def test_prior_mrs_2_is_flagged_as_a_trial_discrepancy_not_silently_cleared(self):
        response = _triage(**{**_ELIGIBLE_KWARGS, "prior_mrs": 2})
        assert any("mRS 2" in r for r in response.thrombectomy_reasoning)


# =====================================================================
# N20 -- a partial NIHSS drove the eligibility gates
# =====================================================================

class TestN20PartialNIHSSDoesNotDriveEligibility:
    """NIHSS >= 4 gates alteplase and >= 6 gates thrombectomy, and unexamined items scored 0,
    so a partial exam systematically under-triaged.
    """

    def test_the_endpoint_reports_the_exam_as_incomplete(self):
        response = _triage(nihss_items={"1a_loc": 2, "9_language": 3},
                           onset_time_hours=2.0)
        assert response.nihss_complete is False

    def test_the_endpoint_will_not_clear_tpa_off_a_partial_exam(self):
        response = _triage(
            nihss_items={"1a_loc": 2, "9_language": 3},
            onset_time_hours=2.0, blood_pressure_systolic=150,
            blood_glucose=110.0, anticoagulant_use=False,
        )
        assert response.tpa_eligible is False
        assert any("NIHSS examination incomplete" in r for r in response.tpa_reasoning)

    def test_the_endpoint_will_not_clear_thrombectomy_off_a_partial_exam(self):
        response = _triage(
            nihss_items={"1a_loc": 2, "9_language": 3},
            onset_time_hours=2.0, lvo_suspected=True, prior_mrs=0, ct_aspects_score=9,
        )
        assert response.thrombectomy_eligible is False

    def test_the_workflow_warns_on_a_partial_exam_not_only_an_absent_one(self):
        """The workflow emitted a warning only when `nihss_scores` was entirely absent."""
        result = AcuteStrokeTriageWorkflow().run({
            "onset_hours": 2.0, "nihss_scores": {"1a_loc": 2, "9_language": 3},
        })
        assert any("NIHSS examination incomplete" in f for f in result.findings)


# =====================================================================
# N10 -- DEFUSE-3 declared out to 24 h with no mismatch-volume or NIHSS gate
# =====================================================================

class TestN10DEFUSE3Criteria:
    """The extended-window branch covered 6-24 h and tested `core_vol < 70 and
    mismatch >= 1.8` only, with no NIHSS re-check after the `lvo or nihss >= 6` entry gate.
    `{onset_hours: 20.0, NIHSS 3, LVO, core 40, mismatch 2.0}` yielded "THROMBECTOMY ELIGIBLE
    (extended window, DEFUSE-3 criteria met)" -- contradicting the guideline string printed
    two lines below it.

    Primary source: DEFUSE-3 protocol / statistical analysis plan, NCT02586415
    (https://cdn.clinicaltrials.gov/large-docs/15/NCT02586415/SAP_000.pdf). Enrollment window
    6-16 h; target mismatch = core < 70 mL AND ratio >= 1.8 AND mismatch volume >= 15 mL;
    baseline NIHSS >= 6; ICA or M1 occlusion; pre-stroke mRS 0-2.
    """

    BASE = dict(
        large_vessel_occlusion=True, prestroke_mrs=0,
        ct_affected_regions=["M1"], age=65,
    )

    def test_the_audit_reproduction_nihss_3_at_20h_is_not_eligible(self):
        result = AcuteStrokeTriageWorkflow().run({
            **self.BASE, "onset_hours": 20.0, "nihss_scores": _full_nihss(3),
            "core_volume_ml": 40.0, "mismatch_ratio": 2.0,
        })
        assert not any("DEFUSE-3 criteria met" in f for f in result.findings)

    def test_defuse3_is_not_applied_beyond_16_hours(self):
        result = AcuteStrokeTriageWorkflow().run({
            **self.BASE, "onset_hours": 20.0, "nihss_scores": _full_nihss(8),
            "core_volume_ml": 40.0, "mismatch_ratio": 2.0, "mismatch_volume_ml": 30.0,
        })
        assert any("DEFUSE-3 not applicable" in f for f in result.findings)
        assert any("6-16h" in f for f in result.findings)

    def test_defuse3_requires_a_mismatch_volume(self):
        result = AcuteStrokeTriageWorkflow().run({
            **self.BASE, "onset_hours": 10.0, "nihss_scores": _full_nihss(8),
            "core_volume_ml": 40.0, "mismatch_ratio": 2.0,
        })
        assert any("mismatch volume (mL) not reported" in f for f in result.findings)
        assert not any("DEFUSE-3 criteria met" in f for f in result.findings)

    def test_defuse3_requires_nihss_at_least_6(self):
        result = AcuteStrokeTriageWorkflow().run({
            **self.BASE, "onset_hours": 10.0, "nihss_scores": _full_nihss(4),
            "core_volume_ml": 40.0, "mismatch_ratio": 2.0, "mismatch_volume_ml": 30.0,
        })
        assert any("below the trial's NIHSS >= 6" in f for f in result.findings)

    def test_a_fully_qualifying_defuse3_case_is_still_eligible(self):
        result = AcuteStrokeTriageWorkflow().run({
            **self.BASE, "onset_hours": 10.0, "nihss_scores": _full_nihss(8),
            "core_volume_ml": 40.0, "mismatch_ratio": 2.0, "mismatch_volume_ml": 30.0,
        })
        assert any("DEFUSE-3" in f and "ELIGIBLE" in f for f in result.findings)

    def test_a_dawn_case_beyond_16h_is_still_eligible(self):
        """DAWN's window is 6-24 h and its three strata were transcribed correctly."""
        result = AcuteStrokeTriageWorkflow().run({
            **self.BASE, "onset_hours": 20.0, "nihss_scores": _full_nihss(12),
            "core_volume_ml": 25.0,
        })
        assert any("DAWN" in f and "ELIGIBLE" in f for f in result.findings)


# =====================================================================
# G3 -- model output filed as retrieved evidence with score 0.8
# =====================================================================

class TestG3ModelOutputIsNotEvidence:
    """`evidence.append({"source": "llm", "text": llm_result, "score": 0.8})`.

    `evidence` is one of the honesty gate's `_EVIDENCE_KEYS`, and the gate deliberately skips
    those subtrees -- retrieved evidence is quoted source material that must not be rewritten.
    So the generated prose was never scanned, never got the decision-support disclaimer, and
    `_evidence_from(payload)` handed it to the claim verifier as the ground truth against
    which to check the answer that produced it. The UI rendered it as
    "Evidence (1 passages) -- [llm] Score: 0.800", indistinguishable from a retrieved PubMed
    passage, and 0.8 was a relevance number attached to something never retrieved.
    """

    class _LLM:
        def generate(self, prompt):
            self.prompt = prompt
            return "...consider idarucizumab reversal..."

    def _with_llm(self):
        llm = self._LLM()
        request = StrokeTriageRequest(
            **_ELIGIBLE_KWARGS, clinical_notes="on apixaban, last dose unknown"
        )
        return asyncio.run(stroke_triage(request, _FakeRequest(llm))), llm

    def test_the_evidence_array_stays_empty(self):
        response, _ = self._with_llm()
        assert response.evidence == []

    def test_no_fabricated_relevance_score_is_attached(self):
        response, _ = self._with_llm()
        assert not any(e.get("score") == 0.8 for e in response.evidence)

    def test_the_prose_lands_in_a_key_the_gate_scans(self):
        from hcls_common.api_gate import _ANSWER_KEYS, _EVIDENCE_KEYS

        response, _ = self._with_llm()
        assert response.discussion is not None
        assert "idarucizumab" in response.discussion
        assert "discussion" in _ANSWER_KEYS
        assert "discussion" not in _EVIDENCE_KEYS

    def test_the_prose_is_labelled_as_model_knowledge(self):
        response, _ = self._with_llm()
        assert "not retrieved from this knowledge base" in response.discussion

    def test_the_prompt_does_not_invite_fabricated_citations(self):
        _, llm = self._with_llm()
        assert "do not cite studies" in llm.prompt.lower()

    @pytest.mark.parametrize("endpoint_name", ["dementia_evaluate", "epilepsy_classify"])
    def test_the_other_two_routes_carry_a_discussion_field_too(self, endpoint_name):
        import api.routes.neuro_clinical as routes

        model = {
            "dementia_evaluate": routes.DementiaEvaluateResponse,
            "epilepsy_classify": routes.EpilepsyClassifyResponse,
        }[endpoint_name]
        assert "discussion" in model.model_fields

    def test_no_route_still_appends_llm_output_to_evidence(self):
        """A source-level guard: the three sites shared one line.

        Comment lines are skipped, since the fix is documented at each site by quoting the
        line it replaced.
        """
        import inspect

        import api.routes.neuro_clinical as routes

        code = [
            line for line in inspect.getsource(routes).splitlines()
            if not line.lstrip().startswith("#")
        ]
        for line in code:
            assert '"source": "llm"' not in line, line
            assert '"score": 0.8' not in line, line
