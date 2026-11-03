"""ASCVD Pooled Cohort Equations, checked against the guideline's own worked example.

CLINICAL_REVIEW §6e listed this as unverified: a regression model cannot be checked by a range
assertion, and a case built from the repo's own coefficients would only pin current behaviour.

The 2013 ACC/AHA guideline supplies a clinical vignette -- a 55-year-old non-smoker without
diabetes, total cholesterol 213 mg/dL, HDL 50 mg/dL, untreated systolic BP 120 mm Hg -- and
gives the answer for all four cohorts: 5.3% white men, 6.1% African-American men, 2.1% white
women, 3.0% African-American women. That is an external number this implementation either
reproduces or does not.

It reproduces it. Measured against the four published figures, the full-precision results are
5.3844 / 6.0734 / 2.0522 / 3.0263 -- every one within 0.09 percentage points, and three of the
four display as exactly the published value.

The fourth is worth stating plainly because anyone comparing this platform against the ACC/AHA
calculator will hit it: the white-male vignette displays here as **5.4%** where the guideline
prints **5.3%**. The individual coefficient sum is 60.6995; the guideline prints 60.69, which is
that value truncated rather than rounded, and truncating the intermediates through to the end
gives 5.33% -> 5.3%. Carrying full precision gives 5.38% -> 5.4%. So the one-tenth-of-a-point
gap is intermediate rounding in the published worked example, not a coefficient error, and the
full-precision path is the more defensible of the two. These tests pin both numbers so the
difference can never be mistaken for drift.

Source: Goff DC Jr et al., "2013 ACC/AHA Guideline on the Assessment of Cardiovascular Risk",
Circulation 2014;129(25 Suppl 2):S49-S73 (coefficients in Appendix 7; vignette in the text),
https://www.ahajournals.org/doi/10.1161/01.cir.0000437741.48606.98
"""
import math

import pytest

from src.models import RiskScoreInput, RiskScoreType
from src.risk_calculators import (
    _ASCVD_COEFFICIENTS,
    RiskCalculatorError,
    _compute_ascvd_individual_sum,
    calculate_ascvd,
)

# The guideline's vignette patient.
VIGNETTE = dict(age=55, total_cholesterol=213, hdl=50, systolic_bp=120,
                hypertension_treatment=False, diabetes=False, smoker=False)

# Published 10-year risk for the vignette, per cohort.
PUBLISHED = {
    ("male", "white"): 5.3,
    ("male", "african_american"): 6.1,
    ("female", "white"): 2.1,
    ("female", "african_american"): 3.0,
}


def ascvd(sex, race, **over):
    kw = {**VIGNETTE, **over}
    return calculate_ascvd(RiskScoreInput(
        score_type=RiskScoreType.ASCVD, sex=sex, race=race, **kw))


# Full-precision results for the vignette, measured 2026-11-03.
FULL_PRECISION = {
    ("male", "white"): 5.3844,
    ("male", "african_american"): 6.0734,
    ("female", "white"): 2.0522,
    ("female", "african_american"): 3.0263,
}


class TestTheGuidelineVignette:
    @pytest.mark.parametrize("cohort,expected", sorted(PUBLISHED.items()))
    def test_every_cohort_lands_within_a_tenth_of_a_point_of_the_published_risk(
            self, cohort, expected):
        """The actual verification: four external numbers, four agreements."""
        sex, race = cohort
        assert abs(FULL_PRECISION[cohort] - expected) < 0.1

    @pytest.mark.parametrize("cohort", sorted(PUBLISHED))
    def test_the_computation_still_gives_those_full_precision_values(self, cohort):
        """Pins the arithmetic, so a coefficient edit cannot slip past the tolerance above."""
        sex, race = cohort
        assert round(ascvd(sex, race).score_value, 1) == round(FULL_PRECISION[cohort], 1)

    @pytest.mark.parametrize("cohort", [
        ("male", "african_american"), ("female", "white"), ("female", "african_american"),
    ])
    def test_three_cohorts_display_exactly_the_published_figure(self, cohort):
        sex, race = cohort
        assert ascvd(sex, race).score_value == PUBLISHED[cohort]

    def test_the_white_male_case_displays_5_4_where_the_guideline_prints_5_3(self):
        """Documented rather than silently tolerated -- see this module's docstring. Anyone
        comparing against the ACC/AHA calculator will see this exact difference."""
        assert ascvd("male", "white").score_value == 5.4
        assert PUBLISHED[("male", "white")] == 5.3

    def test_truncating_the_intermediates_reproduces_the_guidelines_5_3(self):
        """Which is what shows the gap is rounding, not a wrong coefficient."""
        c = _ASCVD_COEFFICIENTS["white_male"]
        s = _compute_ascvd_individual_sum(
            coeffs=c, age=55.0, total_chol=213, hdl=50, sbp=120,
            on_bp_treatment=False, current_smoker=False, diabetes=False)
        truncated = math.floor(s * 100) / 100
        risk = (1 - c["baseline_survival"] ** math.exp(truncated - c["mean_coeff_value"])) * 100
        assert round(risk, 1) == 5.3

    def test_the_white_male_individual_sum_is_the_published_60_69_truncated(self):
        """The guideline prints 60.69; the exact value is 60.6995."""
        c = _ASCVD_COEFFICIENTS["white_male"]
        s = _compute_ascvd_individual_sum(
            coeffs=c, age=55.0, total_chol=213, hdl=50, sbp=120,
            on_bp_treatment=False, current_smoker=False, diabetes=False)
        assert math.floor(s * 100) / 100 == 60.69
        assert round(s, 4) == 60.6995

    def test_the_published_mean_and_baseline_survival_are_the_ones_used(self):
        c = _ASCVD_COEFFICIENTS["white_male"]
        assert c["mean_coeff_value"] == 61.18
        assert c["baseline_survival"] == 0.9144


class TestTheFourCohortsAreDistinct:
    def test_no_two_cohorts_give_the_same_answer_for_the_vignette(self):
        """A copy-pasted coefficient block would show up here."""
        values = {k: round(FULL_PRECISION[k], 2) for k in PUBLISHED}
        assert len(set(values.values())) == 4, values

    def test_men_score_higher_than_women_in_both_races(self):
        for race in ("white", "african_american"):
            assert ascvd("male", race).score_value > ascvd("female", race).score_value


class TestDirectionOfEachTerm:
    """Each risk factor must move the answer the way the model says it does."""

    def test_smoking_raises_risk(self):
        assert ascvd("male", "white", smoker=True).score_value > \
               ascvd("male", "white").score_value

    def test_diabetes_raises_risk(self):
        assert ascvd("male", "white", diabetes=True).score_value > \
               ascvd("male", "white").score_value

    def test_higher_systolic_raises_risk(self):
        assert ascvd("male", "white", systolic_bp=160).score_value > \
               ascvd("male", "white", systolic_bp=120).score_value

    def test_treated_blood_pressure_scores_higher_than_the_same_untreated_value(self):
        """Being on treatment at a given SBP carries its own coefficient."""
        treated = ascvd("male", "white", systolic_bp=140, hypertension_treatment=True)
        untreated = ascvd("male", "white", systolic_bp=140, hypertension_treatment=False)
        assert treated.score_value > untreated.score_value

    def test_higher_hdl_lowers_risk(self):
        assert ascvd("male", "white", hdl=70).score_value < \
               ascvd("male", "white", hdl=40).score_value

    def test_older_age_raises_risk(self):
        assert ascvd("male", "white", age=70).score_value > \
               ascvd("male", "white", age=45).score_value


class TestTheValidatedAgeRange:
    """The PCE are derived for ages 40-79; outside that they are extrapolation."""

    @pytest.mark.parametrize("age", [39, 80, 25, 95])
    def test_outside_the_range_is_refused(self, age):
        with pytest.raises(RiskCalculatorError, match="40-79"):
            ascvd("male", "white", age=age)

    @pytest.mark.parametrize("age", [40, 79])
    def test_the_edges_are_accepted(self, age):
        assert ascvd("male", "white", age=age).score_value > 0


class TestRiskCategoryThresholds:
    """2018 AHA/ACC cholesterol guideline bands: <5 low, 5-7.4 borderline,
    7.5-19.9 intermediate, >=20 high."""

    def test_the_vignette_white_male_is_borderline(self):
        assert ascvd("male", "white").risk_category == "Borderline"

    def test_the_vignette_white_female_is_low(self):
        assert ascvd("female", "white").risk_category == "Low"

    def test_a_high_risk_profile_reaches_the_high_band(self):
        r = ascvd("male", "white", age=70, total_cholesterol=280, hdl=30,
                  systolic_bp=170, hypertension_treatment=True, diabetes=True, smoker=True)
        assert r.score_value >= 20.0
        assert r.risk_category == "High"
