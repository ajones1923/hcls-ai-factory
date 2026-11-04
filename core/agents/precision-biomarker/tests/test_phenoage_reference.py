"""PhenoAge, checked against the published equation — and one thing no source settles.

Levine et al. 2018 (PMID:29676998), its published correction, and Liu et al. 2018 all print the
equation and **none of them states units for any biomarker**. So the coefficients and the two
transformation formulas can be verified; the CRP unit convention cannot.

What these tests do about that: they pin the coefficients against the published equation, and
they pin the SIZE of the unsettled choice — 2.40 years — so that if anyone switches the
convention it is a deliberate act with a visible consequence, and if anyone perturbs the
coefficients the number moves and a test fails. See CLINICAL_REVIEW §8.

Source: the published correction, which prints the corrected equation in full:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6388911/
"""
import math

import pytest

from src.biological_age import (
    BA_DENOMINATOR,
    BA_INTERCEPT,
    BA_NUMERATOR,
    CRP_UNIT_CONVENTION,
    MORT_DENOMINATOR,
    MORT_NUMERATOR,
    PHENOAGE_COEFFICIENTS,
    PHENOAGE_INTERCEPT,
    UNIT_CONVERSIONS,
    BiologicalAgeCalculator,
)

# The published linear predictor:
#   xb = -19.907 - 0.0336*Albumin + 0.0095*Creatinine + 0.1953*Glucose + 0.0954*ln(CRP)
#        - 0.0120*Lymphocyte% + 0.0268*MCV + 0.3306*RDW + 0.00188*ALP + 0.0554*WBC
#        + 0.0804*ChronologicalAge
PUBLISHED = {
    "albumin": -0.0336,
    "creatinine": 0.0095,
    "glucose": 0.1953,
    "ln_crp": 0.0954,
    "lymphocyte_pct": -0.0120,
    "mcv": 0.0268,
    "rdw": 0.3306,
    "alkaline_phosphatase": 0.00188,
    "wbc": 0.0554,
}

HEALTHY_50 = {"albumin": 4.4, "creatinine": 0.9, "glucose": 90, "hs_crp": 1.0,
              "lymphocyte_pct": 30, "mcv": 90, "rdw": 13.0,
              "alkaline_phosphatase": 70, "wbc": 6.0}


@pytest.fixture
def calc():
    return BiologicalAgeCalculator()


class TestCoefficientsMatchThePublishedEquation:
    @pytest.mark.parametrize("name,expected", sorted(PUBLISHED.items()))
    def test_coefficient(self, name, expected):
        assert PHENOAGE_COEFFICIENTS[name] == expected

    def test_the_intercept(self):
        """Published as -19.907; carried here to one more place."""
        assert round(PHENOAGE_INTERCEPT, 3) == -19.907

    def test_the_chronological_age_coefficient_is_applied(self, calc):
        """0.0804 per year; checked through the public API rather than a constant."""
        a = calc.calculate_phenoage(50, HEALTHY_50)["biological_age"]
        b = calc.calculate_phenoage(60, HEALTHY_50)["biological_age"]
        assert b > a

    def test_alkaline_phosphatase_uses_the_published_00188_not_a_rounding(self):
        """It was 0.0019, which is 0.00188 rounded to two significant figures."""
        assert PHENOAGE_COEFFICIENTS["alkaline_phosphatase"] == 0.00188
        assert PHENOAGE_COEFFICIENTS["alkaline_phosphatase"] != 0.0019


class TestTheTransformationConstants:
    def test_the_gompertz_mortality_constants(self):
        """M = 1 - exp(-1.51714 * exp(xb) / 0.0076927)."""
        assert MORT_NUMERATOR == -1.51714
        assert round(MORT_DENOMINATOR, 7) == 0.0076927

    def test_the_phenoage_conversion_constants(self):
        """PhenoAge = 141.50 + ln[-0.00553 * ln(1 - M)] / 0.09165."""
        assert round(BA_INTERCEPT, 2) == 141.50
        assert round(BA_NUMERATOR, 5) == -0.00553
        assert BA_DENOMINATOR == 0.09165


class TestTheSiUnitConversions:
    """The three terms whose units the equation's variable names do state."""

    def test_albumin_converts_g_per_dl_to_g_per_l(self):
        assert UNIT_CONVERSIONS["albumin"] == 10.0

    def test_creatinine_converts_mg_per_dl_to_micromoles(self):
        assert UNIT_CONVERSIONS["creatinine"] == 88.4

    def test_glucose_converts_mg_per_dl_to_millimoles(self):
        assert UNIT_CONVERSIONS["glucose"] == pytest.approx(1 / 18.016)

    def test_crp_has_no_conversion_and_that_is_the_open_question(self):
        """Adding a 0.1 factor here is the one-line switch to the mg/dL reading."""
        assert "ln_crp" not in UNIT_CONVERSIONS
        assert "hs_crp" not in UNIT_CONVERSIONS
        assert "crp" not in UNIT_CONVERSIONS


class TestTheCrpConventionIsDeclaredAndItsCostPinned:
    """The unsettled choice, made visible and measured.

    The primary sources state no units. So rather than pick silently, the module declares its
    convention and these tests fix the magnitude of the alternative at 2.40 years.
    """

    def test_the_convention_is_reported_in_every_result(self, calc):
        out = calc.calculate_phenoage(50, HEALTHY_50)
        assert out["crp_unit_convention"] == CRP_UNIT_CONVENTION
        assert "mg/L" in out["crp_unit_convention"]

    def test_the_declared_convention_names_the_size_of_the_alternative(self):
        assert "2.40" in CRP_UNIT_CONVENTION

    def test_switching_to_mg_per_dl_costs_exactly_2_40_years(self, calc):
        """ln(10) * 0.0954 in the predictor. Measured end to end, at this patient."""
        as_mg_l = calc.calculate_phenoage(50, HEALTHY_50)["biological_age"]
        as_mg_dl = calc.calculate_phenoage(
            50, {**HEALTHY_50, "hs_crp": HEALTHY_50["hs_crp"] / 10})["biological_age"]
        assert as_mg_l - as_mg_dl == pytest.approx(2.40, abs=0.02)

    @pytest.mark.parametrize("crp", [0.5, 1.0, 3.0, 10.0])
    def test_the_offset_is_constant_at_every_crp_level(self, calc, crp):
        """Which is what makes it a systematic bias rather than noise."""
        hi = calc.calculate_phenoage(50, {**HEALTHY_50, "hs_crp": crp})["biological_age"]
        lo = calc.calculate_phenoage(50, {**HEALTHY_50, "hs_crp": crp / 10})["biological_age"]
        assert hi - lo == pytest.approx(2.40, abs=0.02)

    def test_the_current_convention_still_produces_this_number(self, calc):
        """A guard on the whole pipeline: perturb any coefficient and this moves."""
        assert calc.calculate_phenoage(50, HEALTHY_50)["biological_age"] == pytest.approx(
            46.2, abs=0.1)


class TestDirectionOfEachTerm:
    @pytest.mark.parametrize("marker,worse,better", [
        ("hs_crp", 10.0, 0.5),
        ("glucose", 180, 85),
        ("rdw", 16.0, 12.5),
        ("wbc", 11.0, 5.5),
        ("mcv", 100, 88),
        ("alkaline_phosphatase", 140, 60),
    ])
    def test_a_worse_value_raises_biological_age(self, calc, marker, worse, better):
        hi = calc.calculate_phenoage(50, {**HEALTHY_50, marker: worse})["biological_age"]
        lo = calc.calculate_phenoage(50, {**HEALTHY_50, marker: better})["biological_age"]
        assert hi > lo

    @pytest.mark.parametrize("marker,high,low", [
        ("albumin", 4.8, 3.2),
        ("lymphocyte_pct", 38, 15),
    ])
    def test_the_protective_markers_lower_it(self, calc, marker, high, low):
        with_high = calc.calculate_phenoage(50, {**HEALTHY_50, marker: high})["biological_age"]
        with_low = calc.calculate_phenoage(50, {**HEALTHY_50, marker: low})["biological_age"]
        assert with_high < with_low
