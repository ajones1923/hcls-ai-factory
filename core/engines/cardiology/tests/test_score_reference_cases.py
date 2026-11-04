"""Reference cases for the point-based scores: patient in, published score out.

Why this file exists. `test_risk_calculators.py` has 117 assertions, of which 13 are range or
bounds checks and 3 name a published source; much of the rest is structural -- that a function
is callable, that the engine dispatches, that bad input raises. That suite is why
CLINICAL_REVIEW §1 can say the EuroSCORE II LVEF coefficients look mis-slotted while every
test passes: nothing compares an output against a number a guideline would give.

`scripts/run_clinical_eval.py` says the same thing in its first line -- the subject suites
"prove the plumbing".

These cases are derived from the published scoring criteria themselves, each cited below, by
adding up the component points for a described patient. Every expected value here is therefore
checkable by hand against the citation, which is the property the existing tests lack. They
deliberately cover the boundaries, because an off-by-one in an age band is the error a range
assertion cannot see.

Scope: the deterministic point scores only. The regression-based calculators (ASCVD Pooled
Cohort Equations, EuroSCORE II, MAGGIC) need published worked examples with their coefficient
tables, which CLINICAL_REVIEW §1-§2 records as outstanding; a case invented from the repo's
own coefficients would only pin current behaviour.
"""
import pytest

from src.models import RiskScoreInput, RiskScoreType
from src.risk_calculators import (
    RiskCalculatorError,
    calculate_cha2ds2_vasc,
    calculate_has_bled,
    calculate_heart_score,
)


# =====================================================================
# CHA2DS2-VASc — Lip GYH et al., Chest 2010;137(2):263-272
#   C 1 · H 1 · A2 (>=75) 2 · D 1 · S2 (stroke/TIA/TE) 2 · V 1 ·
#   A (65-74) 1 · Sc (female) 1   — maximum 9
# =====================================================================
def cha(**kw):
    kw.setdefault("score_type", RiskScoreType.CHA2DS2_VASC)
    extra = kw.pop("extra", None)
    return calculate_cha2ds2_vasc(RiskScoreInput(**kw), extra=extra)


class TestCha2ds2VascReferenceCases:
    def test_healthy_young_man_scores_zero(self):
        assert cha(age=40, sex="male").score_value == 0

    def test_healthy_young_woman_scores_one_for_sex_alone(self):
        """Sc is the only point; this is the published "female with no risk factors" case."""
        assert cha(age=40, sex="female").score_value == 1

    def test_every_component_at_once_is_nine(self):
        """C1+H1+A2(2)+D1+S2(2)+V1+Sc1 = 9. The A(65-74) point cannot co-occur with A2."""
        r = cha(age=80, sex="female", congestive_heart_failure=True,
                hypertension_treatment=True, diabetes=True, history_of_stroke=True,
                vascular_disease=True)
        assert r.score_value == 9

    def test_the_two_age_bands_are_mutually_exclusive(self):
        """A 76-year-old gets A2=2, not A2=2 plus A=1."""
        assert cha(age=76, sex="male").score_value == 2

    @pytest.mark.parametrize("age,expected", [
        (64, 0),    # below both bands
        (65, 1),    # lower edge of A
        (74, 1),    # upper edge of A
        (75, 2),    # lower edge of A2
    ])
    def test_age_band_boundaries(self, age, expected):
        assert cha(age=age, sex="male").score_value == expected

    def test_stroke_history_is_worth_two_not_one(self):
        assert cha(age=50, sex="male", history_of_stroke=True).score_value == 2

    def test_a_worked_case(self):
        """72-year-old man, hypertension, diabetes, prior TIA:
        H1 + D1 + S2(2) + A(65-74)1 = 5."""
        r = cha(age=72, sex="male", hypertension_treatment=True, diabetes=True,
                history_of_stroke=True)
        assert r.score_value == 5

    def test_an_unknown_sex_is_refused_rather_than_scored(self):
        """The guard fires at the model layer -- RiskScoreInput validates `sex` -- before the
        calculator's own check is reached. Either way it is refused, not scored."""
        with pytest.raises(Exception, match="(?i)sex must be"):
            cha(age=70, sex="unspecified")


# =====================================================================
# HAS-BLED — Pisters R et al., Chest 2010;138(5):1093-1100
#   H 1 · A abnormal renal 1 · A abnormal liver 1 · S 1 · B 1 · L 1 ·
#   E (>65) 1 · D drugs 1 · D alcohol 1   — maximum 9
#   (renal and liver score separately; drugs and alcohol score separately)
# =====================================================================
def hb(**kw):
    kw.setdefault("score_type", RiskScoreType.HAS_BLED)
    extra = kw.pop("extra", None)
    return calculate_has_bled(RiskScoreInput(**kw), extra=extra)


class TestHasBledReferenceCases:
    def test_no_risk_factors_scores_zero(self):
        assert hb(age=50, sex="male").score_value == 0

    def test_every_component_at_once_is_nine(self):
        r = hb(age=80, sex="male", renal_disease=True, liver_disease=True,
               history_of_stroke=True, history_of_bleeding=True, labile_inr=True,
               antiplatelet_nsaid=True, alcohol_excess=True,
               extra={"hypertension_uncontrolled": True})
        assert r.score_value == 9

    def test_renal_and_liver_score_separately(self):
        """The 'A' letter covers two conditions worth one point each, not one shared point."""
        both = hb(age=50, sex="male", renal_disease=True, liver_disease=True).score_value
        renal_only = hb(age=50, sex="male", renal_disease=True).score_value
        assert renal_only == 1 and both == 2

    def test_drugs_and_alcohol_score_separately(self):
        both = hb(age=50, sex="male", antiplatelet_nsaid=True,
                  alcohol_excess=True).score_value
        assert both == 2

    @pytest.mark.parametrize("age,expected", [(65, 0), (66, 1)])
    def test_elderly_is_over_65_not_65_and_over(self, age, expected):
        """Published criterion is age > 65."""
        assert hb(age=age, sex="male").score_value == expected

    def test_a_worked_case(self):
        """70-year-old on aspirin with prior bleeding and labile INR:
        E1 + D(drugs)1 + B1 + L1 = 4 — the published high-bleeding-risk threshold."""
        r = hb(age=70, sex="male", history_of_bleeding=True, labile_inr=True,
               antiplatelet_nsaid=True)
        assert r.score_value == 4
        assert "high" in r.risk_category.lower()


# =====================================================================
# HEART — Six AJ et al., Neth Heart J 2008;16(6):191-196;
#         Backus BE et al., Int J Cardiol 2013;168(3):2153-2158
#   History 0-2 · ECG 0-2 · Age (<45 / 45-64 / >=65) 0-2 ·
#   Risk factors (0 / 1-2 / >=3 or known atherosclerosis) 0-2 ·
#   Troponin (<=1x / 1-3x / >3x ULN) 0-2   — maximum 10
#   Bands: 0-3 low, 4-6 moderate, 7-10 high
# =====================================================================
def heart(age, **extra):
    extra.setdefault("history_suspicion", 0)
    extra.setdefault("ecg_finding", 0)
    extra.setdefault("troponin_ratio", 0.5)
    return calculate_heart_score(
        RiskScoreInput(score_type=RiskScoreType.HEART, age=age), extra=extra)


class TestHeartReferenceCases:
    def test_the_lowest_possible_score_is_zero(self):
        assert heart(40).score_value == 0

    def test_the_highest_possible_score_is_ten(self):
        r = heart(70, history_suspicion=2, ecg_finding=2, troponin_ratio=5.0,
                  num_risk_factors=4)
        assert r.score_value == 10

    @pytest.mark.parametrize("age,expected", [(44, 0), (45, 1), (64, 1), (65, 2)])
    def test_age_band_boundaries(self, age, expected):
        assert heart(age).score_value == expected

    @pytest.mark.parametrize("ratio,expected", [
        (1.0, 0),     # at the ULN is still normal
        (1.5, 1),     # 1-3x
        (3.0, 1),     # at 3x is still the 1-3x band
        (3.5, 2),     # above 3x
    ])
    def test_troponin_band_boundaries(self, ratio, expected):
        assert heart(40, troponin_ratio=ratio).score_value == expected

    @pytest.mark.parametrize("n,expected", [(0, 0), (1, 1), (2, 1), (3, 2)])
    def test_risk_factor_band_boundaries(self, n, expected):
        assert heart(40, num_risk_factors=n).score_value == expected

    def test_known_atherosclerosis_scores_two_regardless_of_count(self):
        assert heart(40, num_risk_factors=0, history_atherosclerosis=True).score_value == 2

    @pytest.mark.parametrize("score_inputs,band", [
        ({}, "low"),                                                    # 0
        ({"history_suspicion": 1, "num_risk_factors": 1}, "low"),       # 2
        ({"age_override": None}, "low"),
    ])
    def test_low_band(self, score_inputs, band):
        score_inputs.pop("age_override", None)
        assert heart(40, **score_inputs).risk_category.lower() == band

    def test_the_band_boundaries_are_3_and_6(self):
        three = heart(70, history_suspicion=1, num_risk_factors=0)          # 2+1 = 3
        assert three.score_value == 3 and three.risk_category.lower() == "low"
        four = heart(70, history_suspicion=1, num_risk_factors=1)           # 2+1+1 = 4
        assert four.score_value == 4 and four.risk_category.lower() == "moderate"
        six = heart(70, history_suspicion=2, num_risk_factors=1, troponin_ratio=1.5)
        assert six.score_value == 6 and six.risk_category.lower() == "moderate"
        seven = heart(70, history_suspicion=2, num_risk_factors=3, troponin_ratio=1.5)
        assert seven.score_value == 7 and seven.risk_category.lower() == "high"

    def test_a_worked_case(self):
        """58-year-old, moderately suspicious history, non-specific ECG repolarisation
        change, two risk factors, troponin 2x ULN:
        H1 + E1 + A1 + R1 + T1 = 5, moderate risk -> admit for serial troponins."""
        r = heart(58, history_suspicion=1, ecg_finding=1, num_risk_factors=2,
                  troponin_ratio=2.0)
        assert r.score_value == 5
        assert r.risk_category.lower() == "moderate"


class TestHeartRefusesAnUndrawnTroponin:
    """A troponin that was never drawn is not a normal troponin.

    The T component is worth up to 2 points, which spans the 0-3/4-6 boundary -- discharge
    versus admit for serial troponins. It used to default to a ratio of 0.0 and the
    interpretation then stated "Troponin <= 1x ULN (ratio 0.00)" about a patient with no
    troponin on record.
    """

    def test_a_missing_troponin_raises(self):
        with pytest.raises(RiskCalculatorError, match="requires a troponin"):
            calculate_heart_score(
                RiskScoreInput(score_type=RiskScoreType.HEART, age=58),
                extra={"history_suspicion": 1, "ecg_finding": 1})

    def test_the_error_names_both_ways_to_supply_it(self):
        with pytest.raises(RiskCalculatorError) as exc:
            calculate_heart_score(
                RiskScoreInput(score_type=RiskScoreType.HEART, age=58),
                extra={"history_suspicion": 0, "ecg_finding": 0})
        assert "troponin_ratio" in str(exc.value)
        assert "troponin_uln" in str(exc.value)

    def test_an_explicit_zero_ratio_is_still_accepted(self):
        """A measured, undetectable troponin is a real result and scores 0."""
        r = calculate_heart_score(
            RiskScoreInput(score_type=RiskScoreType.HEART, age=40),
            extra={"history_suspicion": 0, "ecg_finding": 0, "troponin_ratio": 0.0})
        assert r.score_value == 0

    def test_a_measured_troponin_on_the_input_is_used(self):
        """0.12 against the 0.04 default ULN is exactly 3x, the TOP of the 1-3x band."""
        r = calculate_heart_score(
            RiskScoreInput(score_type=RiskScoreType.HEART, age=40, troponin=0.12),
            extra={"history_suspicion": 0, "ecg_finding": 0})
        assert r.score_value == 1

    def test_a_troponin_above_three_times_the_uln_scores_two(self):
        r = calculate_heart_score(
            RiskScoreInput(score_type=RiskScoreType.HEART, age=40, troponin=0.20),
            extra={"history_suspicion": 0, "ecg_finding": 0})
        assert r.score_value == 2
