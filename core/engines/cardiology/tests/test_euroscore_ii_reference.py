"""EuroSCORE II coefficients, pinned to Table 6 of the paper that defines them.

CLINICAL_REVIEW §1 suspected a copy-paste here: `lvef_21_30` was byte-identical to
`creat_cc_50_85`, and `lvef_31_50` (0.1084443) matched nothing anyone could find. That suspicion
was right, and worse than recorded -- FOUR coefficients were wrong.

  creat_cc_lt_50  0.6521653 -> 0.8592256
  creat_cc_50_85  0.3150652 -> 0.303553     (held the LV-moderate value)
  lvef_31_50      0.1084443 -> 0.3150652
  lvef_21_30      0.3150652 -> 0.8084096    (held the LV-moderate value)

The LV error was the one that mattered. Poor ventricular function (EF 21-30) was scored with
the MODERATE coefficient, 0.3151 instead of 0.8084. For a 72-year-old with poor LV, NYHA III,
urgent, two procedures, predicted operative mortality went from 2.93% to 4.83% once corrected.
EuroSCORE II is used to decide whether to operate, so this understated the risk of surgery for
exactly the patients in whom the number decides.

All 2,046 cardiology tests passed before this fix and after it. Monotonic ordering was
preserved, so every output stayed plausible -- which is precisely why a suite of range
assertions could not see a coefficient wrong by a factor of 2.6. That is the argument for this
file existing.

Source: Nashef SAM, Roques F, Sharples LD, Nilsson J, Smith C, Goldstone AR, Lockowandt U.
"EuroSCORE II", Eur J Cardiothorac Surg 2012;41(4):734-745, Table 6 (final logistic model).
https://academic.oup.com/ejcts/article/41/4/734/646622
Cross-checked against two independent calculator implementations.
"""
import math

import pytest

from src.risk_calculators import (
    _EUROSCORE_II_BETA0,
    _EUROSCORE_II_COEFFICIENTS,
    _euroscore_age_term,
)

# Table 6 of Nashef 2012, final multivariable model.
TABLE_6 = {
    "age": 0.0285181,
    "female": 0.2196434,
    "creat_cc_lt_50": 0.8592256,
    "creat_cc_50_85": 0.303553,
    "creat_on_dialysis": 0.6421508,
    "extracardiac_arteriopathy": 0.5360268,
    "previous_cardiac_surgery": 1.118599,
    "chronic_lung_disease": 0.1886564,
    "active_endocarditis": 0.6194522,
    "diabetes_on_insulin": 0.3542749,
    "nyha_class_ii": 0.1070545,
    "nyha_class_iii": 0.2958358,
    "nyha_class_iv": 0.5597929,
    "ccs_class_4_angina": 0.2226147,
    "lvef_31_50": 0.3150652,
    "lvef_21_30": 0.8084096,
    "lvef_le_20": 0.9346919,
    "recent_mi": 0.1528943,
    "pulm_hyp_moderate": 0.1788899,
    "pulm_hyp_severe": 0.3491475,
    "urgency_urgent": 0.3174673,
    "urgency_emergency": 0.7039121,
    "weight_single_non_cabg": 0.0062118,
    "weight_two_procedures": 0.5521478,
    "weight_three_or_more": 0.9724533,
    "surgery_thoracic_aorta": 0.6527205,
}

# The values whose correction this file exists to lock down.
CORRECTED = {
    "creat_cc_lt_50": (0.6521653, 0.8592256),
    "creat_cc_50_85": (0.3150652, 0.303553),
    "lvef_31_50": (0.1084443, 0.3150652),
    "lvef_21_30": (0.3150652, 0.8084096),
}


class TestEveryCoefficientMatchesTable6:
    @pytest.mark.parametrize("name,expected", sorted(TABLE_6.items()))
    def test_coefficient(self, name, expected):
        assert _EUROSCORE_II_COEFFICIENTS[name] == expected, (
            f"{name} is {_EUROSCORE_II_COEFFICIENTS[name]}, Table 6 gives {expected}")

    def test_the_constant(self):
        assert _EUROSCORE_II_BETA0 == -5.324537


class TestTheCorrectedValuesStayCorrected:
    @pytest.mark.parametrize("name,values", sorted(CORRECTED.items()))
    def test_the_old_wrong_value_is_gone(self, name, values):
        wrong, right = values
        assert _EUROSCORE_II_COEFFICIENTS[name] != wrong
        assert _EUROSCORE_II_COEFFICIENTS[name] == right

    def test_lv_function_is_strictly_ordered_by_severity(self):
        """Moderate < poor < very poor. The bug preserved ordering, so ordering alone is not
        sufficient -- but a future edit that breaks it should still fail here."""
        c = _EUROSCORE_II_COEFFICIENTS
        assert c["lvef_31_50"] < c["lvef_21_30"] < c["lvef_le_20"]

    def test_no_lv_coefficient_shares_a_value_with_a_renal_coefficient(self):
        """The copy-paste signature CLINICAL_REVIEW §1 spotted. It must not come back."""
        c = _EUROSCORE_II_COEFFICIENTS
        lv = {c["lvef_31_50"], c["lvef_21_30"], c["lvef_le_20"]}
        renal = {c["creat_cc_lt_50"], c["creat_cc_50_85"], c["creat_on_dialysis"]}
        assert not (lv & renal)

    def test_severe_renal_impairment_outweighs_moderate(self):
        c = _EUROSCORE_II_COEFFICIENTS
        assert c["creat_cc_50_85"] < c["creat_cc_lt_50"]


class TestTheAgeTerm:
    """Nashef 2012: the age variable is 1 at age <= 60 and rises by one per year after."""

    @pytest.mark.parametrize("age,count", [
        (40, 1), (55, 1), (60, 1),      # the floor is 1, not 0
        (61, 2), (62, 3), (70, 11), (80, 21),
    ])
    def test_the_count_matches_the_published_coding(self, age, count):
        assert _euroscore_age_term(age) == pytest.approx(TABLE_6["age"] * count)

    def test_a_sixty_year_old_contributes_one_unit_not_zero(self):
        """It used to return 0.0 here, making every patient's predictor short by one unit."""
        assert _euroscore_age_term(60) == pytest.approx(TABLE_6["age"])
        assert _euroscore_age_term(60) > 0

    def test_each_additional_year_adds_exactly_one_coefficient(self):
        for age in (61, 65, 70, 75):
            step = _euroscore_age_term(age + 1) - _euroscore_age_term(age)
            assert step == pytest.approx(TABLE_6["age"])


class TestTheClinicalImpactOfTheFix:
    """A regression test on the magnitude, so the understatement cannot silently return."""

    @staticmethod
    def _risk(terms):
        x = _EUROSCORE_II_BETA0 + sum(terms)
        return 100 * math.exp(x) / (1 + math.exp(x))

    def test_poor_lv_patient_predicts_the_corrected_mortality(self):
        """72yo, poor LV (EF 21-30), NYHA III, urgent, two procedures.
        Was 2.93% with the moderate coefficient in the poor slot; 4.83% corrected."""
        c = _EUROSCORE_II_COEFFICIENTS
        risk = self._risk([
            _euroscore_age_term(72), c["lvef_21_30"], c["nyha_class_iii"],
            c["urgency_urgent"], c["weight_two_procedures"],
        ])
        assert risk == pytest.approx(4.83, abs=0.05)

    def test_the_old_coefficients_would_have_understated_it(self):
        c = _EUROSCORE_II_COEFFICIENTS
        with_old_lv = self._risk([
            TABLE_6["age"] * (72 - 60), 0.3150652, c["nyha_class_iii"],
            c["urgency_urgent"], c["weight_two_procedures"],
        ])
        assert with_old_lv == pytest.approx(2.93, abs=0.05)
        assert with_old_lv < 3.0 < 4.0 < self._risk([
            _euroscore_age_term(72), c["lvef_21_30"], c["nyha_class_iii"],
            c["urgency_urgent"], c["weight_two_procedures"]])
