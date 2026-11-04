"""A risk factor nobody recorded is unknown, not absent.

CHA2DS2-VASc and HAS-BLED read every clinical component as ``inp.<field> or False``, so a
patient record that reached the engine with age and sex populated and the comorbidity flags
null produced a complete-looking score from age and sex alone:

    calculate_all_applicable({"age": 70, "sex": "male"}) ->
        cha2ds2_vasc  1  Low-Moderate  "consider oral anticoagulation (Class IIa)"
        has_bled      1  Low-Moderate  "relatively low bleeding risk. Anticoagulation can
                                        generally be continued."

A 70-year-old man in atrial fibrillation with heart failure, hypertension and diabetes has a
true CHA2DS2-VASc of 5 -- Class I anticoagulation, ~7% annual stroke risk by this module's own
table. He scored 1. At 64 he would have scored 0 and been told "No anticoagulation or
antiplatelet therapy recommended (Class III: No benefit)". The HAS-BLED result was the mirror
image, computed without renal function, liver function, bleeding history, INR lability, or
drug and alcohol history -- which means there is no conservative reading of a defaulted score:
the defaults understate stroke risk AND understate bleeding risk.

``RiskScoreResult`` has no field that can distinguish "this component was absent in the data"
from "this component was absent in the patient", so the fix is the one
``calculate_heart_score`` already applies to a missing troponin: refuse, and name what is
needed. These tests pin that refusal, both at the calculator and at the
``calculate_all_applicable`` applicability gate that would otherwise swallow it.
"""
import pytest

from src.models import RiskScoreInput, RiskScoreType
from src.risk_calculators import (
    RiskCalculatorEngine,
    RiskCalculatorError,
    calculate_cha2ds2_vasc,
    calculate_has_bled,
)


_CHA_COMPONENTS = (
    "congestive_heart_failure", "hypertension_treatment", "diabetes",
    "history_of_stroke", "vascular_disease",
)

_HB_COMPONENTS = (
    "hypertension_treatment", "renal_disease", "liver_disease", "history_of_stroke",
    "history_of_bleeding", "labile_inr", "antiplatelet_nsaid", "alcohol_excess",
)


def _complete_cha(**overrides):
    kw = {c: False for c in _CHA_COMPONENTS}
    kw.update(overrides)
    return RiskScoreInput(score_type=RiskScoreType.CHA2DS2_VASC, **kw)


def _complete_hb(**overrides):
    kw = {c: False for c in _HB_COMPONENTS}
    kw.update(overrides)
    return RiskScoreInput(score_type=RiskScoreType.HAS_BLED, **kw)


class TestCha2ds2VascRefusesPartialData:
    def test_age_and_sex_alone_are_refused(self):
        """The audit's reproduction case: this used to return score 1, "Low-Moderate"."""
        inp = RiskScoreInput(score_type=RiskScoreType.CHA2DS2_VASC, age=70, sex="male")
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_cha2ds2_vasc(inp)
        assert "cannot be scored" in str(exc.value)

    @pytest.mark.parametrize("omitted", _CHA_COMPONENTS)
    def test_each_component_is_named_when_it_is_the_one_missing(self, omitted):
        """The refusal has to say WHICH field to supply, or it is just an error."""
        inp = _complete_cha(age=70, sex="male", **{omitted: None})
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_cha2ds2_vasc(inp)
        expected = "hypertension" if omitted == "hypertension_treatment" else omitted
        assert expected in str(exc.value)

    def test_all_the_missing_components_are_named_at_once(self):
        inp = RiskScoreInput(score_type=RiskScoreType.CHA2DS2_VASC, age=70, sex="male")
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_cha2ds2_vasc(inp)
        message = str(exc.value)
        for field in ("congestive_heart_failure", "diabetes", "history_of_stroke",
                      "vascular_disease", "hypertension"):
            assert field in message

    def test_an_explicit_false_is_accepted(self):
        """False is a clinician saying the patient does not have it; None is not."""
        result = calculate_cha2ds2_vasc(_complete_cha(age=70, sex="male"))
        assert result.score_value == 1  # A(65-74) only

    def test_hypertension_may_come_from_extra_instead(self):
        kw = {c: False for c in _CHA_COMPONENTS if c != "hypertension_treatment"}
        inp = RiskScoreInput(score_type=RiskScoreType.CHA2DS2_VASC, age=70, sex="male", **kw)
        result = calculate_cha2ds2_vasc(inp, extra={"hypertension": True})
        assert result.score_value == 2  # A(65-74) + H

    def test_the_sick_patient_the_default_used_to_hide_scores_five(self):
        """70-year-old man in AF with CHF, hypertension and diabetes: C1+H1+D1+A(65-74)1... = 4.

        Plus a prior stroke (S2=2) he is 6; the audit's stated true score of 5 is this patient
        with a prior TIA and no CHF. Pinned here with the full component set so the score can
        never again come from age and sex.
        """
        result = calculate_cha2ds2_vasc(_complete_cha(
            age=70, sex="male", congestive_heart_failure=True,
            hypertension_treatment=True, diabetes=True, history_of_stroke=True,
        ))
        assert result.score_value == 6
        assert "Class I" in " ".join(result.recommendations)


class TestHasBledRefusesPartialData:
    def test_age_alone_is_refused(self):
        inp = RiskScoreInput(score_type=RiskScoreType.HAS_BLED, age=70)
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_has_bled(inp)
        assert "cannot be scored" in str(exc.value)

    @pytest.mark.parametrize("omitted", [c for c in _HB_COMPONENTS if c != "hypertension_treatment"])
    def test_each_component_is_named_when_it_is_the_one_missing(self, omitted):
        inp = _complete_hb(age=70, **{omitted: None})
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_has_bled(inp)
        assert omitted in str(exc.value)

    def test_the_h_component_accepts_any_of_its_three_inputs(self):
        """Explicit flag, or an SBP to compare against 160, or the treatment flag."""
        kw = {c: False for c in _HB_COMPONENTS if c != "hypertension_treatment"}
        by_extra = calculate_has_bled(
            RiskScoreInput(score_type=RiskScoreType.HAS_BLED, age=70, **kw),
            extra={"hypertension_uncontrolled": True},
        )
        by_sbp = calculate_has_bled(
            RiskScoreInput(score_type=RiskScoreType.HAS_BLED, age=70, systolic_bp=170.0, **kw)
        )
        by_treatment = calculate_has_bled(_complete_hb(age=70, hypertension_treatment=True))
        assert by_extra.score_value == by_sbp.score_value == by_treatment.score_value == 2

    def test_no_h_input_at_all_is_refused_rather_than_called_controlled(self):
        kw = {c: False for c in _HB_COMPONENTS if c != "hypertension_treatment"}
        inp = RiskScoreInput(score_type=RiskScoreType.HAS_BLED, age=70, **kw)
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_has_bled(inp)
        assert "hypertension_uncontrolled" in str(exc.value)

    def test_an_explicit_false_everywhere_is_accepted(self):
        result = calculate_has_bled(_complete_hb(age=70))
        assert result.score_value == 1  # E (age > 65) only


class TestCalculateAllApplicableSkipsRatherThanGuesses:
    """The applicability gate has to agree with the calculators, or it swallows the refusal.

    ``calculate_all_applicable`` catches RiskCalculatorError and continues, so a calculator
    that refuses is simply dropped -- which is the right outcome, but only if the gate does not
    first advertise the calculator as applicable on two fields.
    """

    def test_age_and_sex_alone_produce_no_cha2ds2_vasc_or_has_bled(self):
        results = RiskCalculatorEngine().calculate_all_applicable({"age": 70, "sex": "male"})
        returned = {r.score_type for r in results}
        assert RiskScoreType.CHA2DS2_VASC not in returned
        assert RiskScoreType.HAS_BLED not in returned

    def test_a_complete_record_still_produces_both(self):
        patient = {
            "age": 70, "sex": "male",
            "congestive_heart_failure": True, "hypertension_treatment": True,
            "diabetes": True, "history_of_stroke": False, "vascular_disease": False,
            "systolic_bp": 150.0, "renal_disease": False, "liver_disease": False,
            "history_of_bleeding": False, "labile_inr": False,
            "antiplatelet_nsaid": False, "alcohol_excess": False,
        }
        results = RiskCalculatorEngine().calculate_all_applicable(patient)
        by_type = {r.score_type: r for r in results}
        assert by_type[RiskScoreType.CHA2DS2_VASC].score_value == 4
        assert RiskScoreType.HAS_BLED in by_type

    def test_the_required_field_lists_name_every_component(self):
        """A component missing from the list is a component the gate will let be guessed."""
        required = RiskCalculatorEngine._REQUIRED_FIELDS
        for component in _CHA_COMPONENTS:
            assert component in required[RiskScoreType.CHA2DS2_VASC]
        for component in ("renal_disease", "liver_disease", "history_of_bleeding",
                          "labile_inr", "antiplatelet_nsaid", "alcohol_excess"):
            assert component in required[RiskScoreType.HAS_BLED]
