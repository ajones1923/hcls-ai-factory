"""Unit tests for Cardiology Intelligence Agent risk calculators.

Tests RiskCalculatorEngine and all six cardiovascular risk scoring systems:
ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, and EuroSCORE II.

Author: Adam Jones
Date: March 2026
"""


import pytest

from src.models import (
    HeartFailureClass,
    RiskScoreInput,
    RiskScoreResult,
    RiskScoreType,
)
from src.risk_calculators import (
    RiskCalculatorEngine,
    RiskCalculatorError,
    calculate_ascvd,
    calculate_cha2ds2_vasc,
    calculate_has_bled,
    calculate_heart_score,
    calculate_maggic,
    calculate_euroscore_ii,
    _clamp,
    _nyha_to_int,
    _require,
    _require_bool,
    _require_float,
    _require_int,
    _require_str,
    _maggic_age_points,
    _maggic_ef_points,
    _maggic_sbp_points,
    _maggic_bmi_points,
    _maggic_creatinine_points,
    _maggic_nyha_points,
    _heart_age_score,
    _heart_history_score,
    _heart_ecg_score,
    _heart_risk_factors_score,
    _heart_troponin_score,
    _ascvd_cohort_key,
)


# ===================================================================
# HELPER FUNCTION TESTS
# ===================================================================


class TestClamp:
    """Tests for the _clamp helper."""

    def test_within_range(self):
        assert _clamp(5.0, 0.0, 10.0) == 5.0

    def test_below_range(self):
        assert _clamp(-1.0, 0.0, 10.0) == 0.0

    def test_above_range(self):
        assert _clamp(15.0, 0.0, 10.0) == 10.0

    def test_at_lower_bound(self):
        assert _clamp(0.0, 0.0, 10.0) == 0.0

    def test_at_upper_bound(self):
        assert _clamp(10.0, 0.0, 10.0) == 10.0


class TestRequireHelpers:
    """Tests for _require, _require_float, _require_int, _require_bool, _require_str."""

    def test_require_returns_value(self):
        assert _require(42, "test") == 42

    def test_require_none_raises(self):
        with pytest.raises(RiskCalculatorError, match="Missing required field"):
            _require(None, "test_field")

    def test_require_float_returns_float(self):
        assert _require_float(3.14, "val") == 3.14

    def test_require_float_none_raises(self):
        with pytest.raises(RiskCalculatorError):
            _require_float(None, "val")

    def test_require_int_returns_int(self):
        assert _require_int(55, "age") == 55

    def test_require_int_none_raises(self):
        with pytest.raises(RiskCalculatorError):
            _require_int(None, "age")

    def test_require_bool_returns_bool(self):
        assert _require_bool(True, "flag") is True

    def test_require_bool_none_raises(self):
        with pytest.raises(RiskCalculatorError):
            _require_bool(None, "flag")

    def test_require_str_returns_str(self):
        assert _require_str("male", "sex") == "male"

    def test_require_str_none_raises(self):
        with pytest.raises(RiskCalculatorError):
            _require_str(None, "sex")

    def test_require_str_empty_raises(self):
        with pytest.raises(RiskCalculatorError):
            _require_str("", "sex")


class TestNyhaToInt:
    """Tests for _nyha_to_int helper."""

    @pytest.mark.parametrize(
        "nyha_class, expected",
        [
            (HeartFailureClass.NYHA_I, 1),
            (HeartFailureClass.NYHA_II, 2),
            (HeartFailureClass.NYHA_III, 3),
            (HeartFailureClass.NYHA_IV, 4),
        ],
    )
    def test_conversion(self, nyha_class, expected):
        assert _nyha_to_int(nyha_class) == expected

    def test_none_returns_none(self):
        assert _nyha_to_int(None) is None


# ===================================================================
# RiskCalculatorEngine TESTS
# ===================================================================


class TestRiskCalculatorEngine:
    """Tests for the RiskCalculatorEngine class."""

    def test_creation(self):
        engine = RiskCalculatorEngine()
        assert engine is not None

    def test_has_calculate_method(self):
        engine = RiskCalculatorEngine()
        assert callable(getattr(engine, "calculate", None))

    def test_has_named_calculator_methods(self):
        engine = RiskCalculatorEngine()
        assert callable(getattr(engine, "calculate_ascvd", None))
        assert callable(getattr(engine, "calculate_heart", None))
        assert callable(getattr(engine, "calculate_cha2ds2_vasc", None))
        assert callable(getattr(engine, "calculate_has_bled", None))
        assert callable(getattr(engine, "calculate_maggic", None))
        assert callable(getattr(engine, "calculate_euroscore_ii", None))

    def test_calculate_dispatches_ascvd(self):
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=55, sex="male", race="white",
            total_cholesterol=213.0, hdl=50.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            diabetes=False, smoker=False,
        )
        result = engine.calculate(inp)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type == RiskScoreType.ASCVD

    def test_calculate_dispatches_cha2ds2_vasc(self):
        engine = RiskCalculatorEngine()
        inp = RiskScoreInput(
            score_type=RiskScoreType.CHA2DS2_VASC,
            age=72, sex="male",
        )
        result = engine.calculate(inp)
        assert result.score_type == RiskScoreType.CHA2DS2_VASC


# ===================================================================
# ASCVD CALCULATOR TESTS
# ===================================================================


class TestASCVDCalculator:
    """Tests for the ASCVD Pooled Cohort Equations calculator."""

    def _make_input(self, **kwargs):
        defaults = dict(
            score_type=RiskScoreType.ASCVD,
            age=55, sex="male", race="white",
            total_cholesterol=213.0, hdl=50.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            diabetes=False, smoker=False,
        )
        defaults.update(kwargs)
        return RiskScoreInput(**defaults)

    def test_basic_white_male(self):
        inp = self._make_input()
        result = calculate_ascvd(inp)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type == RiskScoreType.ASCVD
        assert 0.0 <= result.score_value <= 100.0

    def test_55yo_white_male_low_risk_range(self):
        inp = self._make_input(
            age=55, sex="male", race="white",
            total_cholesterol=213.0, hdl=50.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            diabetes=False, smoker=False,
        )
        result = calculate_ascvd(inp)
        # Expect somewhere in the low-to-intermediate range
        assert 2.0 <= result.score_value <= 15.0

    def test_risk_category_low(self):
        inp = self._make_input(
            age=40, sex="female", race="white",
            total_cholesterol=170.0, hdl=60.0,
            systolic_bp=110.0,
        )
        result = calculate_ascvd(inp)
        assert result.risk_category == "Low"
        assert result.score_value < 5.0

    def test_risk_category_high(self):
        inp = self._make_input(
            age=70, sex="male", race="white",
            total_cholesterol=280.0, hdl=35.0,
            systolic_bp=160.0,
            hypertension_treatment=True,
            diabetes=True, smoker=True,
        )
        result = calculate_ascvd(inp)
        assert result.risk_category == "High"
        assert result.score_value >= 20.0

    def test_age_below_40_raises(self):
        inp = self._make_input(age=39)
        with pytest.raises(RiskCalculatorError, match="ages 40-79"):
            calculate_ascvd(inp)

    def test_age_above_79_raises(self):
        inp = self._make_input(age=80)
        with pytest.raises(RiskCalculatorError, match="ages 40-79"):
            calculate_ascvd(inp)

    def test_age_boundary_40_valid(self):
        inp = self._make_input(age=40)
        result = calculate_ascvd(inp)
        assert isinstance(result, RiskScoreResult)

    def test_age_boundary_79_valid(self):
        inp = self._make_input(age=79)
        result = calculate_ascvd(inp)
        assert isinstance(result, RiskScoreResult)

    def test_missing_age_raises(self):
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            sex="male", race="white",
            total_cholesterol=213.0, hdl=50.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            diabetes=False, smoker=False,
        )
        with pytest.raises(RiskCalculatorError, match="age"):
            calculate_ascvd(inp)

    def test_missing_sex_raises(self):
        inp = RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=55, race="white",
            total_cholesterol=213.0, hdl=50.0,
            systolic_bp=120.0,
            hypertension_treatment=False,
            diabetes=False, smoker=False,
        )
        with pytest.raises(RiskCalculatorError, match="sex"):
            calculate_ascvd(inp)

    def test_missing_cholesterol_raises(self):
        inp = self._make_input()
        inp.total_cholesterol = None
        with pytest.raises(RiskCalculatorError, match="total_cholesterol"):
            calculate_ascvd(inp)

    def test_african_american_male(self):
        inp = self._make_input(race="african_american")
        result = calculate_ascvd(inp)
        assert isinstance(result, RiskScoreResult)

    def test_african_american_female(self):
        inp = self._make_input(sex="female", race="black")
        result = calculate_ascvd(inp)
        assert isinstance(result, RiskScoreResult)

    def test_other_race_uses_white_coefficients(self):
        inp = self._make_input(race="asian")
        result = calculate_ascvd(inp)
        assert isinstance(result, RiskScoreResult)
        # Should get a note about PCE derivation
        assert any("PCE" in r or "cohort" in r.lower() for r in result.recommendations)

    def test_smoker_increases_risk(self):
        nonsmoker = calculate_ascvd(self._make_input(smoker=False))
        smoker = calculate_ascvd(self._make_input(smoker=True))
        assert smoker.score_value > nonsmoker.score_value

    def test_diabetes_increases_risk(self):
        no_dm = calculate_ascvd(self._make_input(diabetes=False))
        dm = calculate_ascvd(self._make_input(diabetes=True))
        assert dm.score_value > no_dm.score_value

    def test_higher_bp_increases_risk(self):
        low_bp = calculate_ascvd(self._make_input(systolic_bp=110.0))
        high_bp = calculate_ascvd(self._make_input(systolic_bp=160.0))
        assert high_bp.score_value > low_bp.score_value

    def test_result_has_recommendations(self):
        result = calculate_ascvd(self._make_input())
        assert len(result.recommendations) > 0

    def test_result_has_guideline_reference(self):
        result = calculate_ascvd(self._make_input())
        assert "Goff" in result.guideline_reference

    def test_result_has_interpretation(self):
        result = calculate_ascvd(self._make_input())
        assert "ASCVD" in result.interpretation

    def test_cohort_key_white_male(self):
        assert _ascvd_cohort_key("male", "white") == "white_male"

    def test_cohort_key_white_female(self):
        assert _ascvd_cohort_key("female", "white") == "white_female"

    def test_cohort_key_aa_male(self):
        assert _ascvd_cohort_key("male", "african_american") == "african_american_male"

    def test_cohort_key_black_female(self):
        assert _ascvd_cohort_key("female", "black") == "african_american_female"

    def test_cohort_key_invalid_sex_raises(self):
        with pytest.raises(RiskCalculatorError, match="Invalid sex"):
            _ascvd_cohort_key("other", "white")


# ===================================================================
# HEART SCORE TESTS
# ===================================================================


class TestHEARTScore:
    """Tests for the HEART Score calculator."""

    def _make_input(self, age=55, troponin=None):
        return RiskScoreInput(
            score_type=RiskScoreType.HEART,
            age=age, troponin=troponin,
        )

    def _make_extra(self, **kwargs):
        defaults = dict(
            history_suspicion=1, ecg_finding=1,
            num_risk_factors=2, history_atherosclerosis=False,
            troponin_uln=0.04,
        )
        defaults.update(kwargs)
        return defaults

    def test_basic_calculation(self):
        inp = self._make_input(age=55, troponin=0.02)
        extra = self._make_extra()
        result = calculate_heart_score(inp, extra)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type == RiskScoreType.HEART

    def test_low_risk_score_0_to_3(self):
        inp = self._make_input(age=40, troponin=0.01)
        extra = self._make_extra(
            history_suspicion=0, ecg_finding=0,
            num_risk_factors=0,
        )
        result = calculate_heart_score(inp, extra)
        assert result.score_value <= 3
        assert result.risk_category == "Low"

    def test_moderate_risk_score_4_to_6(self):
        inp = self._make_input(age=55, troponin=0.08)
        extra = self._make_extra(
            history_suspicion=1, ecg_finding=1,
            num_risk_factors=2, troponin_uln=0.04,
        )
        result = calculate_heart_score(inp, extra)
        assert 4 <= result.score_value <= 6
        assert result.risk_category == "Moderate"

    def test_high_risk_score_7_to_10(self):
        inp = self._make_input(age=70, troponin=0.5)
        extra = self._make_extra(
            history_suspicion=2, ecg_finding=2,
            num_risk_factors=3, history_atherosclerosis=True,
            troponin_uln=0.04,
        )
        result = calculate_heart_score(inp, extra)
        assert result.score_value >= 7
        assert result.risk_category == "High"

    def test_missing_history_suspicion_raises(self):
        inp = self._make_input()
        extra = {"ecg_finding": 1}
        with pytest.raises(RiskCalculatorError, match="history_suspicion"):
            calculate_heart_score(inp, extra)

    def test_missing_ecg_finding_raises(self):
        inp = self._make_input()
        extra = {"history_suspicion": 1}
        with pytest.raises(RiskCalculatorError, match="ecg_finding"):
            calculate_heart_score(inp, extra)

    def test_missing_age_raises(self):
        inp = RiskScoreInput(score_type=RiskScoreType.HEART)
        extra = self._make_extra()
        with pytest.raises(RiskCalculatorError, match="age"):
            calculate_heart_score(inp, extra)

    def test_result_has_recommendations(self):
        inp = self._make_input(age=55, troponin=0.02)
        result = calculate_heart_score(inp, self._make_extra())
        assert len(result.recommendations) > 0

    def test_result_has_interpretation(self):
        inp = self._make_input(age=55, troponin=0.02)
        result = calculate_heart_score(inp, self._make_extra())
        assert "HEART" in result.interpretation


class TestHEARTScoreComponents:
    """Tests for individual HEART score component helpers."""

    @pytest.mark.parametrize("suspicion, expected_score", [(0, 0), (1, 1), (2, 2)])
    def test_history_score(self, suspicion, expected_score):
        score, _ = _heart_history_score(suspicion)
        assert score == expected_score

    @pytest.mark.parametrize("ecg, expected_score", [(0, 0), (1, 1), (2, 2)])
    def test_ecg_score(self, ecg, expected_score):
        score, _ = _heart_ecg_score(ecg)
        assert score == expected_score

    @pytest.mark.parametrize(
        "age, expected_score",
        [(30, 0), (44, 0), (45, 1), (64, 1), (65, 2), (80, 2)],
    )
    def test_age_score(self, age, expected_score):
        score, _ = _heart_age_score(age)
        assert score == expected_score

    def test_risk_factors_none(self):
        score, _ = _heart_risk_factors_score(0, False)
        assert score == 0

    def test_risk_factors_1_2(self):
        score, _ = _heart_risk_factors_score(2, False)
        assert score == 1

    def test_risk_factors_3_plus(self):
        score, _ = _heart_risk_factors_score(3, False)
        assert score == 2

    def test_risk_factors_atherosclerosis(self):
        score, _ = _heart_risk_factors_score(0, True)
        assert score == 2

    @pytest.mark.parametrize(
        "ratio, expected_score",
        [(0.5, 0), (1.0, 0), (1.5, 1), (3.0, 1), (3.1, 2), (10.0, 2)],
    )
    def test_troponin_score(self, ratio, expected_score):
        score, _ = _heart_troponin_score(ratio)
        assert score == expected_score


# ===================================================================
# CHA2DS2-VASc TESTS
# ===================================================================


class TestCHA2DS2VASc:
    """Tests for the CHA2DS2-VASc score calculator."""

    def _make_input(self, **kwargs):
        defaults = dict(
            score_type=RiskScoreType.CHA2DS2_VASC,
            age=65, sex="male",
        )
        defaults.update(kwargs)
        return RiskScoreInput(**defaults)

    def test_male_score_0(self):
        inp = self._make_input(age=50, sex="male")
        result = calculate_cha2ds2_vasc(inp)
        assert result.score_value == 0
        assert result.risk_category == "Low"

    def test_female_score_1_is_low(self):
        inp = self._make_input(age=50, sex="female")
        result = calculate_cha2ds2_vasc(inp)
        assert result.score_value == 1  # sex category only
        assert result.risk_category == "Low"

    def test_male_score_1_low_moderate(self):
        inp = self._make_input(age=67, sex="male")
        result = calculate_cha2ds2_vasc(inp)
        assert result.score_value == 1  # Age 65-74 = 1
        assert result.risk_category == "Low-Moderate"

    @pytest.mark.parametrize(
        "age, sex, chf, htn, dm, stroke, vasc, expected_score",
        [
            (50, "male", False, False, False, False, False, 0),
            (50, "female", False, False, False, False, False, 1),
            (76, "male", True, True, True, True, True, 8),
            (65, "male", True, False, False, False, False, 2),
            (50, "male", False, False, False, True, False, 2),
        ],
    )
    def test_parametrized_scores(self, age, sex, chf, htn, dm, stroke, vasc, expected_score):
        inp = self._make_input(
            age=age, sex=sex,
            congestive_heart_failure=chf,
            hypertension_treatment=htn,
            diabetes=dm,
            history_of_stroke=stroke,
            vascular_disease=vasc,
        )
        result = calculate_cha2ds2_vasc(inp)
        assert result.score_value == expected_score

    def test_max_score_is_9(self):
        inp = self._make_input(
            age=80, sex="female",
            congestive_heart_failure=True,
            hypertension_treatment=True,
            diabetes=True,
            history_of_stroke=True,
            vascular_disease=True,
        )
        result = calculate_cha2ds2_vasc(inp)
        assert result.score_value == 9

    def test_missing_age_raises(self):
        inp = RiskScoreInput(score_type=RiskScoreType.CHA2DS2_VASC, sex="male")
        with pytest.raises(RiskCalculatorError, match="age"):
            calculate_cha2ds2_vasc(inp)

    def test_missing_sex_raises(self):
        inp = RiskScoreInput(score_type=RiskScoreType.CHA2DS2_VASC, age=65)
        with pytest.raises(RiskCalculatorError, match="sex"):
            calculate_cha2ds2_vasc(inp)

    def test_result_has_guideline_reference(self):
        inp = self._make_input()
        result = calculate_cha2ds2_vasc(inp)
        assert "Lip" in result.guideline_reference

    def test_high_score_mentions_anticoagulation(self):
        inp = self._make_input(
            age=76, sex="male",
            congestive_heart_failure=True,
            hypertension_treatment=True,
            diabetes=True,
        )
        result = calculate_cha2ds2_vasc(inp)
        assert any("anticoagulation" in r.lower() for r in result.recommendations)


# ===================================================================
# HAS-BLED TESTS
# ===================================================================


class TestHASBLED:
    """Tests for the HAS-BLED score calculator."""

    def _make_input(self, **kwargs):
        defaults = dict(
            score_type=RiskScoreType.HAS_BLED,
            age=70,
        )
        defaults.update(kwargs)
        return RiskScoreInput(**defaults)

    def test_basic_calculation(self):
        inp = self._make_input()
        result = calculate_has_bled(inp)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type == RiskScoreType.HAS_BLED

    def test_minimal_score(self):
        inp = self._make_input(age=50)
        result = calculate_has_bled(inp)
        assert result.score_value == 0
        assert result.risk_category == "Low-Moderate"

    def test_elderly_adds_one(self):
        inp = self._make_input(age=66)
        result = calculate_has_bled(inp)
        assert result.score_value >= 1  # E = elderly

    def test_score_with_all_factors(self):
        inp = self._make_input(
            age=70,
            systolic_bp=170.0,
            renal_disease=True,
            liver_disease=True,
            history_of_stroke=True,
            history_of_bleeding=True,
            labile_inr=True,
            antiplatelet_nsaid=True,
            alcohol_excess=True,
        )
        result = calculate_has_bled(inp)
        assert result.score_value == 9

    def test_score_3_is_high(self):
        inp = self._make_input(
            age=70,
            renal_disease=True,
            liver_disease=True,
        )
        result = calculate_has_bled(inp)
        assert result.score_value >= 3
        assert result.risk_category == "High"

    def test_score_below_3_is_low_moderate(self):
        inp = self._make_input(age=60)  # not elderly (>65)
        result = calculate_has_bled(inp)
        assert result.score_value < 3
        assert result.risk_category == "Low-Moderate"

    def test_missing_age_raises(self):
        inp = RiskScoreInput(score_type=RiskScoreType.HAS_BLED)
        with pytest.raises(RiskCalculatorError, match="age"):
            calculate_has_bled(inp)

    def test_uncontrolled_htn_from_sbp(self):
        inp = self._make_input(systolic_bp=165.0)
        result = calculate_has_bled(inp)
        # H component should be 1 because SBP > 160
        # Also E for elderly, so score >= 2
        assert result.score_value >= 2

    def test_result_has_recommendations(self):
        inp = self._make_input()
        result = calculate_has_bled(inp)
        assert len(result.recommendations) > 0

    def test_high_score_does_not_contraindicate_anticoag(self):
        inp = self._make_input(
            age=70, renal_disease=True, liver_disease=True,
            history_of_stroke=True,
        )
        result = calculate_has_bled(inp)
        assert any("NOT contraindicate" in r or "not a reason to withhold" in r.lower()
                    for r in result.recommendations)


# ===================================================================
# MAGGIC SCORE TESTS
# ===================================================================


class TestMAGGIC:
    """Tests for the MAGGIC Heart Failure Risk Score calculator."""

    def _make_input(self, **kwargs):
        defaults = dict(
            score_type=RiskScoreType.MAGGIC,
            age=65, sex="male",
            lvef=30.0,
            nyha_class=HeartFailureClass.NYHA_III,
            systolic_bp=110.0,
            bmi=25.0,
            creatinine=1.2,
            diabetes=True,
            smoker=False,
        )
        defaults.update(kwargs)
        return RiskScoreInput(**defaults)

    def test_basic_calculation(self):
        inp = self._make_input()
        result = calculate_maggic(inp)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type == RiskScoreType.MAGGIC

    def test_result_has_1y_mortality_in_interpretation(self):
        inp = self._make_input()
        result = calculate_maggic(inp)
        assert "1-year" in result.interpretation.lower() or "mortality" in result.interpretation.lower()

    def test_missing_lvef_raises(self):
        inp = self._make_input()
        inp.lvef = None
        with pytest.raises(RiskCalculatorError, match="lvef"):
            calculate_maggic(inp)

    def test_missing_nyha_raises(self):
        inp = self._make_input()
        inp.nyha_class = None
        with pytest.raises(RiskCalculatorError, match="nyha_class"):
            calculate_maggic(inp)

    def test_missing_bmi_raises(self):
        inp = self._make_input()
        inp.bmi = None
        with pytest.raises(RiskCalculatorError, match="bmi"):
            calculate_maggic(inp)

    def test_result_has_recommendations(self):
        inp = self._make_input()
        result = calculate_maggic(inp)
        assert len(result.recommendations) > 0

    def test_result_has_guideline_reference(self):
        inp = self._make_input()
        result = calculate_maggic(inp)
        assert "Pocock" in result.guideline_reference

    def test_higher_nyha_increases_score(self):
        inp_ii = self._make_input(nyha_class=HeartFailureClass.NYHA_II)
        inp_iv = self._make_input(nyha_class=HeartFailureClass.NYHA_IV)
        r_ii = calculate_maggic(inp_ii)
        r_iv = calculate_maggic(inp_iv)
        assert r_iv.score_value > r_ii.score_value

    def test_lower_ef_increases_score(self):
        inp_35 = self._make_input(lvef=35.0)
        inp_15 = self._make_input(lvef=15.0)
        r_35 = calculate_maggic(inp_35)
        r_15 = calculate_maggic(inp_15)
        assert r_15.score_value > r_35.score_value


class TestMAGGICPointFunctions:
    """Tests for MAGGIC point calculation helpers."""

    @pytest.mark.parametrize(
        "age, lvef, expected",
        [
            (50, 25, 0),   # EF<=30, age<55
            (57, 25, 1),   # EF<=30, age 55-59
            (62, 25, 2),   # EF<=30, age 60-64
            (67, 35, 6),   # EF 31-40, age 65-69
            (82, 45, 15),  # EF>40, age>=80
        ],
    )
    def test_age_points(self, age, lvef, expected):
        assert _maggic_age_points(age, lvef) == expected

    @pytest.mark.parametrize(
        "lvef, expected",
        [(15, 7), (22, 6), (27, 5), (32, 3), (37, 2), (42, 1), (50, 0)],
    )
    def test_ef_points(self, lvef, expected):
        assert _maggic_ef_points(lvef) == expected

    @pytest.mark.parametrize(
        "sbp, expected",
        [(100, 5), (115, 3), (125, 2), (135, 1), (145, 0)],
    )
    def test_sbp_points(self, sbp, expected):
        assert _maggic_sbp_points(sbp) == expected

    @pytest.mark.parametrize(
        "bmi, expected",
        [(14, 6), (18, 5), (22, 3), (27, 2), (35, 0)],
    )
    def test_bmi_points(self, bmi, expected):
        assert _maggic_bmi_points(bmi) == expected

    @pytest.mark.parametrize(
        "creat, expected",
        [(0.8, 0), (1.0, 1), (1.2, 2), (1.4, 3), (1.6, 4), (1.8, 5), (2.0, 6), (3.0, 8)],
    )
    def test_creatinine_points(self, creat, expected):
        assert _maggic_creatinine_points(creat) == expected

    @pytest.mark.parametrize(
        "nyha, expected",
        [(1, 0), (2, 2), (3, 6), (4, 8)],
    )
    def test_nyha_points(self, nyha, expected):
        assert _maggic_nyha_points(nyha) == expected


# ===================================================================
# EuroSCORE II TESTS
# ===================================================================


class TestEuroSCOREII:
    """Tests for the EuroSCORE II calculator."""

    def _make_input(self, **kwargs):
        defaults = dict(
            score_type=RiskScoreType.EUROSCORE_II,
            age=65, sex="male",
        )
        defaults.update(kwargs)
        return RiskScoreInput(**defaults)

    def test_basic_calculation(self):
        inp = self._make_input()
        result = calculate_euroscore_ii(inp)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type == RiskScoreType.EUROSCORE_II

    def test_result_has_score_value(self):
        inp = self._make_input()
        result = calculate_euroscore_ii(inp)
        assert result.score_value >= 0.0

    def test_result_has_risk_category(self):
        inp = self._make_input()
        result = calculate_euroscore_ii(inp)
        assert result.risk_category in ("Low", "Intermediate", "High")

    def test_missing_age_raises(self):
        inp = RiskScoreInput(score_type=RiskScoreType.EUROSCORE_II, sex="male")
        with pytest.raises(RiskCalculatorError, match="age"):
            calculate_euroscore_ii(inp)

    def test_missing_sex_raises(self):
        inp = RiskScoreInput(score_type=RiskScoreType.EUROSCORE_II, age=65)
        with pytest.raises(RiskCalculatorError, match="sex"):
            calculate_euroscore_ii(inp)

    def test_result_has_guideline_reference(self):
        inp = self._make_input()
        result = calculate_euroscore_ii(inp)
        assert "Nashef" in result.guideline_reference

    def test_older_patient_higher_risk(self):
        young = calculate_euroscore_ii(self._make_input(age=50))
        old = calculate_euroscore_ii(self._make_input(age=80))
        assert old.score_value > young.score_value

    def test_result_has_recommendations(self):
        inp = self._make_input()
        result = calculate_euroscore_ii(inp)
        assert len(result.recommendations) > 0


# ===================================================================
# ALL CALCULATORS: RiskScoreResult VALIDATION
# ===================================================================


class TestAllCalculatorsReturnCorrectResult:
    """Verify that all calculator methods return a RiskScoreResult with required fields."""

    @pytest.mark.parametrize(
        "calculator_fn, inp_kwargs, extra",
        [
            (
                calculate_ascvd,
                dict(score_type=RiskScoreType.ASCVD, age=55, sex="male", race="white",
                     total_cholesterol=213.0, hdl=50.0, systolic_bp=120.0,
                     hypertension_treatment=False, diabetes=False, smoker=False),
                None,
            ),
            (
                calculate_heart_score,
                dict(score_type=RiskScoreType.HEART, age=55, troponin=0.02),
                {"history_suspicion": 1, "ecg_finding": 1, "num_risk_factors": 2},
            ),
            (
                calculate_cha2ds2_vasc,
                dict(score_type=RiskScoreType.CHA2DS2_VASC, age=70, sex="male"),
                None,
            ),
            (
                calculate_has_bled,
                dict(score_type=RiskScoreType.HAS_BLED, age=70),
                None,
            ),
            (
                calculate_maggic,
                dict(score_type=RiskScoreType.MAGGIC, age=65, sex="male",
                     lvef=30.0, nyha_class=HeartFailureClass.NYHA_III,
                     systolic_bp=110.0, bmi=25.0, creatinine=1.2, diabetes=True),
                None,
            ),
            (
                calculate_euroscore_ii,
                dict(score_type=RiskScoreType.EUROSCORE_II, age=65, sex="male"),
                None,
            ),
        ],
    )
    def test_result_fields(self, calculator_fn, inp_kwargs, extra):
        inp = RiskScoreInput(**inp_kwargs)
        result = calculator_fn(inp, extra)
        assert isinstance(result, RiskScoreResult)
        assert result.score_type is not None
        assert isinstance(result.score_value, float)
        assert isinstance(result.risk_category, str)
        assert len(result.risk_category) > 0
        assert isinstance(result.interpretation, str)
        assert len(result.interpretation) > 0
        assert isinstance(result.recommendations, list)
        assert isinstance(result.guideline_reference, str)
