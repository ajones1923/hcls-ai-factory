"""Validated cardiovascular risk calculators for the Cardiology Intelligence Agent.

Author: Adam Jones
Date: March 2026

Implements six evidence-based cardiovascular risk scoring systems:
  1. ASCVD Pooled Cohort Equations (Goff 2013)
  2. HEART Score for chest pain
  3. CHA2DS2-VASc for AF stroke risk
  4. HAS-BLED for bleeding risk on anticoagulation
  5. MAGGIC Heart Failure Risk Score
  6. EuroSCORE II for cardiac surgical mortality

All calculators accept RiskScoreInput and return RiskScoreResult with
score value, risk category, human-readable interpretation, actionable
recommendations, and guideline references.

The RiskScoreInput model carries the most common cardiovascular fields.
For calculator-specific inputs not present on RiskScoreInput (e.g., HEART
score history_suspicion, MAGGIC COPD flag, EuroSCORE II operative details),
callers pass an ``extra`` dict alongside the input.  The ``RiskCalculatorEngine``
handles this seamlessly through ``calculate_all_applicable``.

Part of the HCLS AI Factory precision medicine platform.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

from src.models import (
    HeartFailureClass,
    RiskScoreInput,
    RiskScoreResult,
    RiskScoreType,
)


# ═══════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════


# ---------------------------------------------------------------------------
# 1. ASCVD Pooled Cohort Equations - Goff 2013 coefficients
#    Source: Goff DC Jr, et al. 2013 ACC/AHA Guideline on the Assessment of
#    Cardiovascular Risk. J Am Coll Cardiol. 2014;63(25 Pt B):2935-2959.
# ---------------------------------------------------------------------------

# Each sub-dictionary stores the beta coefficients for a sex/race cohort.
# Keys:
#   ln_age, ln_age_sq       -- log(age) and its square
#   ln_tc, ln_age_ln_tc     -- log(total cholesterol) terms
#   ln_hdl, ln_age_ln_hdl   -- log(HDL-C) terms
#   ln_treated_sbp, ln_age_ln_treated_sbp   -- if on BP treatment
#   ln_untreated_sbp, ln_age_ln_untreated_sbp -- if not on BP treatment
#   current_smoker, ln_age_smoker            -- smoking terms
#   diabetes                                  -- diabetes term
#   mean_coeff_value                          -- group mean coefficient sum
#   baseline_survival                         -- 10-year baseline survival S0

_ASCVD_COEFFICIENTS: Dict[str, Dict[str, float]] = {
    "white_female": {
        "ln_age": -29.799,
        "ln_age_sq": 4.884,
        "ln_tc": 13.540,
        "ln_age_ln_tc": -3.114,
        "ln_hdl": -13.578,
        "ln_age_ln_hdl": 3.149,
        "ln_treated_sbp": 2.019,
        "ln_age_ln_treated_sbp": 0.0,
        "ln_untreated_sbp": 1.957,
        "ln_age_ln_untreated_sbp": 0.0,
        "current_smoker": 7.574,
        "ln_age_smoker": -1.665,
        "diabetes": 0.661,
        "mean_coeff_value": -29.18,
        "baseline_survival": 0.9665,
    },
    "african_american_female": {
        "ln_age": 17.114,
        "ln_age_sq": 0.0,
        "ln_tc": 0.940,
        "ln_age_ln_tc": 0.0,
        "ln_hdl": -18.920,
        "ln_age_ln_hdl": 4.475,
        "ln_treated_sbp": 29.291,
        "ln_age_ln_treated_sbp": -6.432,
        "ln_untreated_sbp": 27.820,
        "ln_age_ln_untreated_sbp": -6.087,
        "current_smoker": 0.691,
        "ln_age_smoker": 0.0,
        "diabetes": 0.874,
        "mean_coeff_value": 86.61,
        "baseline_survival": 0.9533,
    },
    "white_male": {
        "ln_age": 12.344,
        "ln_age_sq": 0.0,
        "ln_tc": 11.853,
        "ln_age_ln_tc": -2.664,
        "ln_hdl": -7.990,
        "ln_age_ln_hdl": 1.769,
        "ln_treated_sbp": 1.797,
        "ln_age_ln_treated_sbp": 0.0,
        "ln_untreated_sbp": 1.764,
        "ln_age_ln_untreated_sbp": 0.0,
        "current_smoker": 7.837,
        "ln_age_smoker": -1.795,
        "diabetes": 0.658,
        "mean_coeff_value": 61.18,
        "baseline_survival": 0.9144,
    },
    "african_american_male": {
        "ln_age": 2.469,
        "ln_age_sq": 0.0,
        "ln_tc": 0.302,
        "ln_age_ln_tc": 0.0,
        "ln_hdl": -0.307,
        "ln_age_ln_hdl": 0.0,
        "ln_treated_sbp": 1.916,
        "ln_age_ln_treated_sbp": 0.0,
        "ln_untreated_sbp": 1.809,
        "ln_age_ln_untreated_sbp": 0.0,
        "current_smoker": 0.549,
        "ln_age_smoker": 0.0,
        "diabetes": 0.645,
        "mean_coeff_value": 19.54,
        "baseline_survival": 0.8954,
    },
}


# ---------------------------------------------------------------------------
# 2. HEART Score MACE rates by category
# ---------------------------------------------------------------------------

_HEART_RISK_THRESHOLDS = {
    "low": {"min": 0, "max": 3, "mace_rate": 1.7},
    "moderate": {"min": 4, "max": 6, "mace_rate": 16.6},
    "high": {"min": 7, "max": 10, "mace_rate": 50.1},
}


# ---------------------------------------------------------------------------
# 3. CHA2DS2-VASc annual stroke rates by score (%)
#    Source: Lip GYH, et al. Chest 2010; 137(2):263-272.
# ---------------------------------------------------------------------------

_CHA2DS2_VASC_STROKE_RATES: Dict[int, float] = {
    0: 0.2,
    1: 0.6,
    2: 2.2,
    3: 3.2,
    4: 4.8,
    5: 7.2,
    6: 9.7,
    7: 11.2,
    8: 10.8,
    9: 12.2,
}


# ---------------------------------------------------------------------------
# 5. MAGGIC score-to-mortality lookup tables
#    Source: Pocock SJ, et al. Eur Heart J. 2013;34(19):1404-1413.
# ---------------------------------------------------------------------------

_MAGGIC_MORTALITY_1Y: Dict[int, float] = {
    0: 1.5, 1: 1.7, 2: 1.9, 3: 2.1, 4: 2.4, 5: 2.7, 6: 3.0,
    7: 3.4, 8: 3.8, 9: 4.2, 10: 4.8, 11: 5.3, 12: 6.0, 13: 6.7,
    14: 7.5, 15: 8.4, 16: 9.4, 17: 10.5, 18: 11.7, 19: 13.1,
    20: 14.6, 21: 16.2, 22: 18.1, 23: 20.1, 24: 22.3, 25: 24.7,
    26: 27.3, 27: 30.1, 28: 33.1, 29: 36.3, 30: 39.7, 31: 43.3,
    32: 47.0, 33: 50.8, 34: 54.7, 35: 58.6, 36: 62.4, 37: 66.1,
    38: 69.7, 39: 73.1, 40: 76.2, 41: 79.1, 42: 81.7, 43: 84.0,
    44: 86.1, 45: 87.9, 46: 89.5, 47: 90.8, 48: 92.0, 49: 93.0,
    50: 93.8,
}

_MAGGIC_MORTALITY_3Y: Dict[int, float] = {
    0: 4.4, 1: 5.0, 2: 5.6, 3: 6.3, 4: 7.1, 5: 7.9, 6: 8.9,
    7: 9.9, 8: 11.1, 9: 12.4, 10: 13.9, 11: 15.5, 12: 17.2,
    13: 19.2, 14: 21.3, 15: 23.6, 16: 26.1, 17: 28.8, 18: 31.7,
    19: 34.7, 20: 37.9, 21: 41.3, 22: 44.7, 23: 48.3, 24: 51.9,
    25: 55.5, 26: 59.1, 27: 62.7, 28: 66.1, 29: 69.4, 30: 72.6,
    31: 75.5, 32: 78.2, 33: 80.7, 34: 83.0, 35: 85.1, 36: 86.9,
    37: 88.5, 38: 89.9, 39: 91.2, 40: 92.2, 41: 93.1, 42: 93.9,
    43: 94.5, 44: 95.1, 45: 95.5, 46: 95.9, 47: 96.3, 48: 96.6,
    49: 96.8, 50: 97.0,
}


# ---------------------------------------------------------------------------
# 6. EuroSCORE II logistic regression coefficients
#    Source: Nashef SA, et al. Eur J Cardiothorac Surg. 2012;41(4):734-745.
# ---------------------------------------------------------------------------

_EUROSCORE_II_BETA0 = -5.324537

_EUROSCORE_II_COEFFICIENTS: Dict[str, float] = {
    # Patient-related factors
    "age": 0.0285181,              # per year over 60
    "female": 0.2196434,
    "creat_cc_lt_50": 0.6521653,   # CrCl < 50 ml/min (not on dialysis)
    "creat_cc_50_85": 0.3150652,   # CrCl 50-85 ml/min
    "creat_on_dialysis": 0.6421508,
    "extracardiac_arteriopathy": 0.5360268,
    "poor_mobility": 0.2407818,
    "previous_cardiac_surgery": 1.118599,
    "chronic_lung_disease": 0.1886564,
    "active_endocarditis": 0.6194522,
    "critical_preoperative_state": 1.0862550,
    "diabetes_on_insulin": 0.3542749,

    # Cardiac-related factors
    "nyha_class_ii": 0.1070545,
    "nyha_class_iii": 0.2958358,
    "nyha_class_iv": 0.5597929,
    "ccs_class_4_angina": 0.2226147,
    "lvef_21_30": 0.3150652,
    "lvef_31_50": 0.1084443,
    "lvef_le_20": 0.9346919,
    "recent_mi": 0.1528943,
    "pulm_hyp_moderate": 0.1788899,    # systolic PA pressure 31-55 mmHg
    "pulm_hyp_severe": 0.3491475,      # systolic PA pressure > 55 mmHg

    # Operation-related factors
    "urgency_urgent": 0.3174673,
    "urgency_emergency": 0.7039121,
    "urgency_salvage": 1.3626972,
    "weight_single_non_cabg": 0.0062118,
    "weight_two_procedures": 0.5521478,
    "weight_three_or_more": 0.9724533,
    "surgery_thoracic_aorta": 0.6527205,
}


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION HELPERS
# ═══════════════════════════════════════════════════════════════════════════════


class RiskCalculatorError(Exception):
    """Raised when input validation fails for a risk calculator."""
    pass


def _require(value: Optional[Any], field_name: str) -> Any:
    """Raise if a required field is None."""
    if value is None:
        raise RiskCalculatorError(
            f"Missing required field '{field_name}' for this calculator."
        )
    return value


def _require_float(value: Optional[float], field_name: str) -> float:
    """Raise if a required float field is None."""
    return float(_require(value, field_name))


def _require_int(value: Optional[int], field_name: str) -> int:
    """Raise if a required int field is None."""
    return int(_require(value, field_name))


def _require_bool(value: Optional[bool], field_name: str) -> bool:
    """Raise if a required bool field is None."""
    return bool(_require(value, field_name))


def _require_str(value: Optional[str], field_name: str) -> str:
    """Raise if a required string field is None or empty."""
    if not value:
        raise RiskCalculatorError(
            f"Missing required field '{field_name}' for this calculator."
        )
    return value


def _clamp(value: float, lo: float, hi: float) -> float:
    """Clamp a value to [lo, hi]."""
    return max(lo, min(hi, value))


def _nyha_to_int(nyha: Optional[HeartFailureClass]) -> Optional[int]:
    """Convert HeartFailureClass enum to integer 1-4."""
    if nyha is None:
        return None
    mapping = {
        HeartFailureClass.NYHA_I: 1,
        HeartFailureClass.NYHA_II: 2,
        HeartFailureClass.NYHA_III: 3,
        HeartFailureClass.NYHA_IV: 4,
    }
    return mapping.get(nyha, None)


def _extra(extra: Optional[Dict], key: str, default: Any = None) -> Any:
    """Safely get a value from the extras dict, falling back to default."""
    if extra is None:
        return default
    return extra.get(key, default)


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 1: ASCVD POOLED COHORT EQUATIONS
# ═══════════════════════════════════════════════════════════════════════════════


def _ascvd_cohort_key(sex: str, race: str) -> str:
    """Return the coefficient dictionary key for the patient's sex and race.

    The PCE has separate equations for four groups:
      white_female, white_male, african_american_female, african_american_male.
    Patients of other races are calculated using the white coefficients (per
    the original 2013 guideline) with appropriate caveats noted.
    """
    sex_lower = sex.strip().lower()
    race_lower = race.strip().lower() if race else "white"

    if race_lower in ("african_american", "african american", "black", "aa"):
        race_prefix = "african_american"
    else:
        race_prefix = "white"

    if sex_lower in ("female", "f"):
        return f"{race_prefix}_female"
    elif sex_lower in ("male", "m"):
        return f"{race_prefix}_male"
    else:
        raise RiskCalculatorError(
            f"Invalid sex value '{sex}'. Must be 'male' or 'female'."
        )


def _compute_ascvd_individual_sum(
    coeffs: Dict[str, float],
    age: float,
    total_chol: float,
    hdl: float,
    sbp: float,
    on_bp_treatment: bool,
    current_smoker: bool,
    diabetes: bool,
) -> float:
    """Compute the individual sum for the Pooled Cohort Equations.

    Returns the raw sum of beta_i * x_i terms (before applying the
    baseline survival function).
    """
    ln_age = math.log(max(age, 1))
    ln_tc = math.log(max(total_chol, 1))
    ln_hdl = math.log(max(hdl, 1))
    ln_sbp = math.log(max(sbp, 1))

    smoker_val = 1.0 if current_smoker else 0.0
    diabetes_val = 1.0 if diabetes else 0.0

    total = 0.0

    # Age terms
    total += coeffs["ln_age"] * ln_age
    total += coeffs["ln_age_sq"] * (ln_age ** 2)

    # Total cholesterol
    total += coeffs["ln_tc"] * ln_tc
    total += coeffs["ln_age_ln_tc"] * ln_age * ln_tc

    # HDL cholesterol
    total += coeffs["ln_hdl"] * ln_hdl
    total += coeffs["ln_age_ln_hdl"] * ln_age * ln_hdl

    # Systolic blood pressure (treated vs. untreated)
    if on_bp_treatment:
        total += coeffs["ln_treated_sbp"] * ln_sbp
        total += coeffs["ln_age_ln_treated_sbp"] * ln_age * ln_sbp
    else:
        total += coeffs["ln_untreated_sbp"] * ln_sbp
        total += coeffs["ln_age_ln_untreated_sbp"] * ln_age * ln_sbp

    # Smoking
    total += coeffs["current_smoker"] * smoker_val
    total += coeffs["ln_age_smoker"] * ln_age * smoker_val

    # Diabetes
    total += coeffs["diabetes"] * diabetes_val

    return total


def calculate_ascvd(
    inp: RiskScoreInput,
    extra: Optional[Dict] = None,
) -> RiskScoreResult:
    """Calculate 10-year ASCVD risk using the Pooled Cohort Equations.

    Valid for adults aged 40-79 without pre-existing ASCVD.

    Uses ``RiskScoreInput`` fields:
        age, sex, race, total_cholesterol, hdl, systolic_bp,
        hypertension_treatment, diabetes, smoker.

    Reference:
        Goff DC Jr, et al. 2013 ACC/AHA Guideline on the Assessment of
        Cardiovascular Risk. J Am Coll Cardiol. 2014;63(25 Pt B):2935-2959.
    """
    # --- Validate required inputs ---
    age = _require_int(inp.age, "age")
    sex = _require_str(inp.sex, "sex")
    race = inp.race or "white"
    total_chol = _require_float(inp.total_cholesterol, "total_cholesterol")
    hdl = _require_float(inp.hdl, "hdl")
    sbp = _require_float(inp.systolic_bp, "systolic_bp")
    on_bp_treatment = _require_bool(inp.hypertension_treatment, "hypertension_treatment")
    diabetes = _require_bool(inp.diabetes, "diabetes")
    current_smoker = _require_bool(inp.smoker, "smoker")

    if age < 40 or age > 79:
        raise RiskCalculatorError(
            f"ASCVD PCE is validated for ages 40-79. Got age={age}."
        )

    # --- Select coefficient set ---
    cohort_key = _ascvd_cohort_key(sex, race)
    coeffs = _ASCVD_COEFFICIENTS[cohort_key]

    # --- Compute individual sum ---
    individual_sum = _compute_ascvd_individual_sum(
        coeffs=coeffs,
        age=float(age),
        total_chol=total_chol,
        hdl=hdl,
        sbp=sbp,
        on_bp_treatment=on_bp_treatment,
        current_smoker=current_smoker,
        diabetes=diabetes,
    )

    # --- Convert to 10-year risk ---
    mean_coeff = coeffs["mean_coeff_value"]
    baseline_surv = coeffs["baseline_survival"]

    risk_10yr = 1.0 - baseline_surv ** math.exp(individual_sum - mean_coeff)
    risk_10yr_pct = _clamp(risk_10yr * 100.0, 0.0, 100.0)

    # --- Categorize ---
    if risk_10yr_pct < 5.0:
        category = "Low"
    elif risk_10yr_pct < 7.5:
        category = "Borderline"
    elif risk_10yr_pct < 20.0:
        category = "Intermediate"
    else:
        category = "High"

    # --- Build interpretation ---
    interpretation = (
        f"10-year ASCVD risk: {risk_10yr_pct:.1f}% ({category} risk). "
        f"Calculated using the Pooled Cohort Equations for the "
        f"{cohort_key.replace('_', ' ')} cohort. "
        f"Inputs: age {age}, TC {total_chol:.0f}, HDL {hdl:.0f}, "
        f"SBP {sbp:.0f} ({'treated' if on_bp_treatment else 'untreated'}), "
        f"{'diabetic' if diabetes else 'non-diabetic'}, "
        f"{'smoker' if current_smoker else 'non-smoker'}."
    )

    # --- Build recommendations ---
    recommendations: List[str] = []

    if risk_10yr_pct < 5.0:
        recommendations.append(
            "Low 10-year risk. Emphasize lifestyle modification: "
            "heart-healthy diet, regular exercise, smoking cessation if applicable."
        )
        recommendations.append(
            "Reassess risk every 4-6 years or if risk factors change."
        )
    elif risk_10yr_pct < 7.5:
        recommendations.append(
            "Borderline risk. Lifestyle modifications are the primary intervention."
        )
        recommendations.append(
            "Consider coronary artery calcium (CAC) scoring for further "
            "risk stratification if treatment decision is uncertain."
        )
        recommendations.append(
            "If risk-enhancing factors are present (e.g., family history of "
            "premature ASCVD, metabolic syndrome, chronic inflammatory conditions, "
            "South Asian ancestry, elevated Lp(a), elevated hsCRP, ABI < 0.9), "
            "statin therapy may be reasonable."
        )
    elif risk_10yr_pct < 20.0:
        recommendations.append(
            "Intermediate risk. Initiate moderate-intensity statin therapy "
            "(e.g., atorvastatin 10-20 mg or rosuvastatin 5-10 mg)."
        )
        recommendations.append(
            "If risk-enhancing factors are present, a moderate-to-high "
            "intensity statin is favoured."
        )
        recommendations.append(
            "Consider coronary artery calcium (CAC) scoring to refine the "
            "estimate. If CAC = 0, statin may be deferred; if CAC >= 100, "
            "statin is indicated."
        )
        recommendations.append(
            "LDL-C target: >= 30% reduction from baseline."
        )
    else:
        recommendations.append(
            "High risk (>= 20%). Initiate high-intensity statin therapy "
            "(atorvastatin 40-80 mg or rosuvastatin 20-40 mg)."
        )
        recommendations.append(
            "LDL-C target: >= 50% reduction from baseline."
        )
        recommendations.append(
            "If LDL-C remains >= 70 mg/dL on maximally tolerated statin, "
            "consider adding ezetimibe."
        )
        recommendations.append(
            "If LDL-C remains >= 70 mg/dL on statin + ezetimibe, "
            "consider PCSK9 inhibitor therapy."
        )
        recommendations.append(
            "Aggressive lifestyle modification: DASH/Mediterranean diet, "
            ">= 150 min/week moderate-intensity aerobic activity, "
            "smoking cessation, weight management."
        )

    if race and race.strip().lower() not in (
        "white", "african_american", "african american", "black", "aa"
    ):
        recommendations.append(
            "Note: The PCE was derived from white and African American cohorts. "
            "Risk may be over- or under-estimated for other racial/ethnic groups. "
            "Consider CAC scoring for refined risk assessment."
        )

    return RiskScoreResult(
        score_type=RiskScoreType.ASCVD,
        score_value=round(risk_10yr_pct, 1),
        risk_category=category,
        interpretation=interpretation,
        recommendations=recommendations,
        guideline_reference=(
            "Goff DC Jr, et al. 2013 ACC/AHA Guideline on the Assessment of "
            "Cardiovascular Risk. J Am Coll Cardiol. 2014;63(25 Pt B):2935-2959. "
            "Also: 2018 AHA/ACC Cholesterol Guideline (Grundy SM, et al. "
            "J Am Coll Cardiol. 2019;73(24):e285-e350)."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 2: HEART SCORE
# ═══════════════════════════════════════════════════════════════════════════════


def _heart_history_score(suspicion: int) -> Tuple[int, str]:
    """Score the history component (0-2).

    0 = slightly suspicious
    1 = moderately suspicious
    2 = highly suspicious
    """
    score = int(_clamp(suspicion, 0, 2))
    labels = {
        0: "Slightly suspicious history",
        1: "Moderately suspicious history",
        2: "Highly suspicious history",
    }
    return score, labels.get(score, "Unknown")


def _heart_ecg_score(ecg_finding: int) -> Tuple[int, str]:
    """Score the ECG component (0-2).

    0 = normal
    1 = non-specific repolarization disturbance
    2 = significant ST deviation
    """
    score = int(_clamp(ecg_finding, 0, 2))
    labels = {
        0: "Normal ECG",
        1: "Non-specific repolarization disturbance",
        2: "Significant ST deviation",
    }
    return score, labels.get(score, "Unknown")


def _heart_age_score(age: int) -> Tuple[int, str]:
    """Score the age component (0-2).

    0 = age < 45
    1 = age 45-64
    2 = age >= 65
    """
    if age < 45:
        return 0, f"Age {age} (< 45)"
    elif age < 65:
        return 1, f"Age {age} (45-64)"
    else:
        return 2, f"Age {age} (>= 65)"


def _heart_risk_factors_score(
    num_risk_factors: int,
    history_atherosclerosis: bool,
) -> Tuple[int, str]:
    """Score the risk factors component (0-2).

    Risk factors counted: hypertension, diabetes, hypercholesterolaemia,
    obesity (BMI > 30), smoking, family history of CAD.

    0 = no known risk factors
    1 = 1-2 risk factors
    2 = >= 3 risk factors or known atherosclerotic disease
    """
    if history_atherosclerosis:
        return 2, "History of atherosclerotic disease"
    elif num_risk_factors >= 3:
        return 2, f"{num_risk_factors} risk factors (>= 3)"
    elif num_risk_factors >= 1:
        return 1, f"{num_risk_factors} risk factor(s)"
    else:
        return 0, "No known risk factors"


def _heart_troponin_score(troponin_ratio: float) -> Tuple[int, str]:
    """Score the troponin component (0-2).

    troponin_ratio = measured troponin / upper limit of normal (ULN).

    0 = <= 1x ULN (normal)
    1 = 1-3x ULN
    2 = > 3x ULN
    """
    if troponin_ratio <= 1.0:
        return 0, f"Troponin <= 1x ULN (ratio {troponin_ratio:.2f})"
    elif troponin_ratio <= 3.0:
        return 1, f"Troponin 1-3x ULN (ratio {troponin_ratio:.2f})"
    else:
        return 2, f"Troponin > 3x ULN (ratio {troponin_ratio:.2f})"


def calculate_heart_score(
    inp: RiskScoreInput,
    extra: Optional[Dict] = None,
) -> RiskScoreResult:
    """Calculate the HEART Score for emergency department chest pain triage.

    Predicts 6-week risk of Major Adverse Cardiac Events (MACE): death,
    MI, or coronary revascularization.

    Uses ``RiskScoreInput`` fields: age, troponin.

    Additional fields from ``extra`` dict:
        history_suspicion (int 0-2): chest pain history suspicion level
        ecg_finding (int 0-2): ECG classification
        num_risk_factors (int): count of CV risk factors
        history_atherosclerosis (bool): prior known atherosclerotic disease
        troponin_uln (float): upper limit of normal for troponin assay

    If ``history_suspicion`` and ``ecg_finding`` are not provided, the
    calculator will raise a ``RiskCalculatorError``.

    Reference:
        Six AJ, et al. Chest pain in the emergency room: value of the HEART
        score. Neth Heart J. 2008;16(6):191-196.
        Backus BE, et al. Int J Cardiol. 2013;168(3):2153-2158.
    """
    # --- Gather inputs ---
    age = _require_int(inp.age, "age")

    # HEART-specific fields from extra dict
    history_suspicion = _extra(extra, "history_suspicion")
    ecg_finding = _extra(extra, "ecg_finding")
    num_risk_factors = _extra(extra, "num_risk_factors", 0)
    history_athero = _extra(extra, "history_atherosclerosis", False)
    troponin_uln = _extra(extra, "troponin_uln", 0.04)  # typical hs-cTnI ULN

    if history_suspicion is None:
        raise RiskCalculatorError(
            "Missing required extra field 'history_suspicion' (int 0-2) "
            "for HEART Score."
        )
    if ecg_finding is None:
        raise RiskCalculatorError(
            "Missing required extra field 'ecg_finding' (int 0-2) "
            "for HEART Score."
        )

    # Compute troponin ratio
    troponin_val = inp.troponin
    if troponin_val is not None and troponin_uln > 0:
        troponin_ratio = troponin_val / troponin_uln
    else:
        troponin_ratio = _extra(extra, "troponin_ratio", 0.0)

    # --- Compute each component ---
    h_score, h_label = _heart_history_score(int(history_suspicion))
    e_score, e_label = _heart_ecg_score(int(ecg_finding))
    a_score, a_label = _heart_age_score(age)
    r_score, r_label = _heart_risk_factors_score(int(num_risk_factors), bool(history_athero))
    t_score, t_label = _heart_troponin_score(float(troponin_ratio))

    total_score = h_score + e_score + a_score + r_score + t_score

    # --- Categorize ---
    if total_score <= 3:
        category = "Low"
        mace_rate = _HEART_RISK_THRESHOLDS["low"]["mace_rate"]
    elif total_score <= 6:
        category = "Moderate"
        mace_rate = _HEART_RISK_THRESHOLDS["moderate"]["mace_rate"]
    else:
        category = "High"
        mace_rate = _HEART_RISK_THRESHOLDS["high"]["mace_rate"]

    # --- Interpretation ---
    interpretation = (
        f"HEART Score: {total_score}/10 ({category} risk). "
        f"Estimated 6-week MACE rate: {mace_rate}%. "
        f"Components: History={h_score} ({h_label}), ECG={e_score} ({e_label}), "
        f"Age={a_score} ({a_label}), Risk factors={r_score} ({r_label}), "
        f"Troponin={t_score} ({t_label})."
    )

    # --- Recommendations ---
    recommendations: List[str] = []

    if total_score <= 3:
        recommendations.append(
            "Low-risk HEART score (0-3). Consider early discharge with "
            "outpatient follow-up."
        )
        recommendations.append(
            "Non-invasive testing (stress test, CT coronary angiography) "
            "can be arranged as outpatient if clinically indicated."
        )
        recommendations.append(
            "Provide clear discharge instructions including chest pain "
            "return precautions."
        )
    elif total_score <= 6:
        recommendations.append(
            "Moderate-risk HEART score (4-6). Admit for observation and "
            "serial troponin monitoring."
        )
        recommendations.append(
            "Obtain cardiology consultation."
        )
        recommendations.append(
            "Non-invasive stress testing or CT coronary angiography "
            "during index admission."
        )
        recommendations.append(
            "Initiate guideline-directed medical therapy: aspirin, "
            "antianginal agents as appropriate."
        )
    else:
        recommendations.append(
            "High-risk HEART score (7-10). Admit to monitored bed / CCU. "
            "High likelihood of ACS."
        )
        recommendations.append(
            "Urgent cardiology consultation for possible invasive strategy "
            "(coronary angiography +/- PCI)."
        )
        recommendations.append(
            "Initiate dual antiplatelet therapy (aspirin + P2Y12 inhibitor) "
            "and anticoagulation per ACS protocol."
        )
        recommendations.append(
            "Serial ECGs and cardiac biomarkers. Continuous telemetry monitoring."
        )
        recommendations.append(
            "Consider early invasive strategy within 24 hours per ACC/AHA "
            "NSTE-ACS guidelines."
        )

    return RiskScoreResult(
        score_type=RiskScoreType.HEART,
        score_value=float(total_score),
        risk_category=category,
        interpretation=interpretation,
        recommendations=recommendations,
        guideline_reference=(
            "Six AJ, et al. Chest pain in the emergency room: value of the "
            "HEART score. Neth Heart J. 2008;16(6):191-196. "
            "Backus BE, et al. Int J Cardiol. 2013;168(3):2153-2158."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 3: CHA2DS2-VASc
# ═══════════════════════════════════════════════════════════════════════════════


def calculate_cha2ds2_vasc(
    inp: RiskScoreInput,
    extra: Optional[Dict] = None,
) -> RiskScoreResult:
    """Calculate CHA2DS2-VASc score for stroke risk in atrial fibrillation.

    Components (max 9):
      C  Congestive heart failure                +1
      H  Hypertension                            +1
      A2 Age >= 75                               +2
      D  Diabetes mellitus                       +1
      S2 Stroke / TIA / thromboembolism          +2
      V  Vascular disease (prior MI, PAD, aortic plaque) +1
      A  Age 65-74                               +1
      Sc Sex category (female)                   +1

    Uses ``RiskScoreInput`` fields:
        age, sex, congestive_heart_failure, hypertension_treatment (as
        proxy for hypertension), diabetes, history_of_stroke,
        vascular_disease.

    Optional ``extra`` fields:
        hypertension (bool): explicit hypertension flag (overrides
            hypertension_treatment if present).

    Reference:
        Lip GYH, et al. Refining clinical risk stratification for predicting
        stroke and thromboembolism in atrial fibrillation. Chest.
        2010;137(2):263-272.
    """
    # --- Validate required inputs ---
    age = _require_int(inp.age, "age")
    sex = _require_str(inp.sex, "sex")
    chf = inp.congestive_heart_failure or False
    # Use explicit hypertension from extra if available; otherwise infer from
    # hypertension_treatment flag on the model.
    hypertension = _extra(extra, "hypertension", inp.hypertension_treatment or False)
    diabetes = inp.diabetes or False
    stroke_tia = inp.history_of_stroke or False
    vascular = inp.vascular_disease or False

    is_female = sex in ("female", "f")
    is_male = sex in ("male", "m")

    if not (is_female or is_male):
        raise RiskCalculatorError(
            f"Invalid sex '{sex}'. Must be 'male' or 'female'."
        )

    # --- Compute score ---
    score = 0
    detail_parts: List[str] = []

    # C - Congestive heart failure
    c_val = 1 if chf else 0
    score += c_val
    detail_parts.append(f"C(CHF)={c_val}")

    # H - Hypertension
    h_val = 1 if hypertension else 0
    score += h_val
    detail_parts.append(f"H(HTN)={h_val}")

    # A2 - Age >= 75
    a2_val = 2 if age >= 75 else 0
    score += a2_val
    detail_parts.append(f"A2(age>=75)={a2_val}")

    # D - Diabetes
    d_val = 1 if diabetes else 0
    score += d_val
    detail_parts.append(f"D(DM)={d_val}")

    # S2 - Stroke/TIA
    s2_val = 2 if stroke_tia else 0
    score += s2_val
    detail_parts.append(f"S2(stroke)={s2_val}")

    # V - Vascular disease
    v_val = 1 if vascular else 0
    score += v_val
    detail_parts.append(f"V(vasc)={v_val}")

    # A - Age 65-74
    a_val = 1 if 65 <= age < 75 else 0
    score += a_val
    detail_parts.append(f"A(65-74)={a_val}")

    # Sc - Sex category (female)
    sc_val = 1 if is_female else 0
    score += sc_val
    detail_parts.append(f"Sc(female)={sc_val}")

    # --- Annual stroke rate ---
    clamped_score = min(score, 9)
    annual_stroke_rate = _CHA2DS2_VASC_STROKE_RATES.get(clamped_score, 15.2)

    # --- Categorize and recommend ---
    recommendations: List[str] = []

    if is_male and score == 0:
        category = "Low"
        recommendations.append(
            "CHA2DS2-VASc 0 in male: low stroke risk. No anticoagulation "
            "or antiplatelet therapy recommended (Class III: No benefit)."
        )
        recommendations.append(
            "Reassess annually or when clinical status changes."
        )
    elif is_female and score == 1:
        category = "Low"
        recommendations.append(
            "CHA2DS2-VASc 1 in female (sex category point only): low stroke "
            "risk. No anticoagulation recommended."
        )
        recommendations.append(
            "Reassess annually or when clinical status changes."
        )
    elif is_male and score == 1:
        category = "Low-Moderate"
        recommendations.append(
            "CHA2DS2-VASc 1 in male: consider oral anticoagulation "
            "(Class IIa recommendation)."
        )
        recommendations.append(
            "Preferred agents: direct oral anticoagulants (DOACs) -- "
            "apixaban, rivaroxaban, edoxaban, or dabigatran -- "
            "are recommended over warfarin (unless mechanical valve or "
            "moderate-to-severe mitral stenosis)."
        )
        recommendations.append(
            "Discuss bleeding risk (HAS-BLED) and patient preferences."
        )
    elif is_female and score == 2:
        category = "Low-Moderate"
        recommendations.append(
            "CHA2DS2-VASc 2 in female (1 additional risk factor beyond sex): "
            "consider oral anticoagulation (Class IIa)."
        )
        recommendations.append(
            "Preferred agents: DOACs over warfarin (unless valvular AF)."
        )
        recommendations.append(
            "Assess bleeding risk with HAS-BLED before initiating."
        )
    else:
        category = "Moderate-High" if score <= 4 else "High"
        recommendations.append(
            f"CHA2DS2-VASc {score}: oral anticoagulation is recommended "
            f"(Class I). Annual stroke rate ~{annual_stroke_rate}%."
        )
        recommendations.append(
            "Preferred: DOAC therapy (apixaban, rivaroxaban, edoxaban, "
            "or dabigatran) unless contraindicated."
        )
        recommendations.append(
            "Warfarin (target INR 2.0-3.0) is an alternative when DOACs "
            "are contraindicated (mechanical valve, moderate-severe mitral "
            "stenosis, severe CKD)."
        )
        recommendations.append(
            "Assess bleeding risk using HAS-BLED. A high HAS-BLED score "
            "is not a reason to withhold anticoagulation but to identify "
            "and address modifiable bleeding risk factors."
        )
        if score >= 6:
            recommendations.append(
                "Very high stroke risk. Ensure anticoagulation adherence. "
                "Consider left atrial appendage occlusion (LAAO) only if "
                "long-term oral anticoagulation is contraindicated."
            )

    # --- Interpretation ---
    interpretation = (
        f"CHA2DS2-VASc score: {score}/9 ({category} risk). "
        f"Estimated annual stroke rate: {annual_stroke_rate}%. "
        f"Components: {', '.join(detail_parts)}."
    )

    return RiskScoreResult(
        score_type=RiskScoreType.CHA2DS2_VASC,
        score_value=float(score),
        risk_category=category,
        interpretation=interpretation,
        recommendations=recommendations,
        guideline_reference=(
            "Lip GYH, et al. Chest. 2010;137(2):263-272. "
            "2020 ESC Guidelines for AF (Hindricks G, et al. Eur Heart J. "
            "2021;42(5):373-498). "
            "2019 AHA/ACC/HRS Focused Update on AF (January CT, et al. "
            "J Am Coll Cardiol. 2019;74(1):104-132)."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 4: HAS-BLED
# ═══════════════════════════════════════════════════════════════════════════════


def calculate_has_bled(
    inp: RiskScoreInput,
    extra: Optional[Dict] = None,
) -> RiskScoreResult:
    """Calculate HAS-BLED bleeding risk score for patients on anticoagulation.

    Components (max 9):
      H  Hypertension (uncontrolled, SBP > 160 mmHg)      +1
      A  Abnormal renal function                           +1
         Abnormal liver function                           +1
      S  Stroke history                                    +1
      B  Bleeding history or predisposition                +1
      L  Labile INR (TTR < 60%)                           +1
      E  Elderly (age > 65)                               +1
      D  Drugs (antiplatelets / NSAIDs)                    +1
         Alcohol excess                                    +1

    Uses ``RiskScoreInput`` fields:
        age, systolic_bp, renal_disease, liver_disease, history_of_stroke,
        history_of_bleeding, labile_inr, antiplatelet_nsaid, alcohol_excess.

    Optional ``extra`` fields:
        hypertension_uncontrolled (bool): explicit override for the H component.

    Reference:
        Pisters R, et al. A novel user-friendly score (HAS-BLED) to assess
        1-year risk of major bleeding in patients with atrial fibrillation.
        Chest. 2010;138(5):1093-1100.
    """
    # --- Gather inputs ---
    age = _require_int(inp.age, "age")

    # H: uncontrolled hypertension
    explicit_uncontrolled = _extra(extra, "hypertension_uncontrolled")
    if explicit_uncontrolled is not None:
        has_uncontrolled_htn = bool(explicit_uncontrolled)
    elif inp.systolic_bp is not None:
        has_uncontrolled_htn = inp.systolic_bp > 160
    elif inp.hypertension_treatment is not None:
        has_uncontrolled_htn = bool(inp.hypertension_treatment)
    else:
        has_uncontrolled_htn = False

    abnormal_renal = inp.renal_disease or False
    abnormal_liver = inp.liver_disease or False
    stroke = inp.history_of_stroke or False
    bleeding = inp.history_of_bleeding or False
    labile_inr = inp.labile_inr or False
    elderly = age > 65
    drugs = inp.antiplatelet_nsaid or False
    alcohol = inp.alcohol_excess or False

    # --- Compute score ---
    score = 0
    detail_parts: List[str] = []

    h_val = 1 if has_uncontrolled_htn else 0
    score += h_val
    detail_parts.append(f"H(uncontrolled HTN)={h_val}")

    ar_val = 1 if abnormal_renal else 0
    score += ar_val
    detail_parts.append(f"A(renal)={ar_val}")

    al_val = 1 if abnormal_liver else 0
    score += al_val
    detail_parts.append(f"A(liver)={al_val}")

    s_val = 1 if stroke else 0
    score += s_val
    detail_parts.append(f"S(stroke)={s_val}")

    b_val = 1 if bleeding else 0
    score += b_val
    detail_parts.append(f"B(bleeding)={b_val}")

    l_val = 1 if labile_inr else 0
    score += l_val
    detail_parts.append(f"L(labile INR)={l_val}")

    e_val = 1 if elderly else 0
    score += e_val
    detail_parts.append(f"E(elderly)={e_val}")

    d_val = 1 if drugs else 0
    score += d_val
    detail_parts.append(f"D(drugs)={d_val}")

    a_val = 1 if alcohol else 0
    score += a_val
    detail_parts.append(f"D(alcohol)={a_val}")

    # --- Categorize ---
    if score < 3:
        category = "Low-Moderate"
    else:
        category = "High"

    # Approximate annual major bleeding rates by score (%)
    bleeding_rates: Dict[int, float] = {
        0: 1.13, 1: 1.02, 2: 1.88, 3: 3.74, 4: 8.70,
        5: 12.50, 6: 12.50, 7: 12.50, 8: 12.50, 9: 12.50,
    }
    annual_bleed_rate = bleeding_rates.get(min(score, 9), 12.50)

    # --- Recommendations ---
    recommendations: List[str] = []

    if score < 3:
        recommendations.append(
            f"HAS-BLED {score}: relatively low bleeding risk. "
            f"Anticoagulation can generally be continued."
        )
        recommendations.append(
            "Review and address any modifiable bleeding risk factors."
        )
    else:
        recommendations.append(
            f"HAS-BLED >= 3 ({score}): high bleeding risk. "
            f"This does NOT contraindicate anticoagulation if indicated "
            f"for stroke prevention."
        )
        recommendations.append(
            "IMPORTANT: A high HAS-BLED score identifies patients who require "
            "closer follow-up and more frequent monitoring, not necessarily "
            "withdrawal of anticoagulation."
        )

    # Modifiable risk factor advice
    modifiable: List[str] = []
    if has_uncontrolled_htn:
        modifiable.append("Optimize blood pressure control (target SBP < 160 mmHg)")
    if labile_inr:
        modifiable.append(
            "Address labile INR: consider switching from warfarin to a DOAC, "
            "improve INR monitoring frequency, or evaluate drug/diet interactions"
        )
    if drugs:
        modifiable.append(
            "Review and minimize concurrent use of NSAIDs and antiplatelet agents; "
            "discontinue if not clearly indicated"
        )
    if alcohol:
        modifiable.append("Counsel on reducing alcohol intake")
    if bleeding:
        modifiable.append(
            "Investigate and treat underlying bleeding source (e.g., GI evaluation)"
        )

    if modifiable:
        recommendations.append(
            "Address modifiable risk factors: " + "; ".join(modifiable) + "."
        )

    recommendations.append(
        "Schedule regular follow-up (at least every 3 months for high-risk "
        "patients) with CBC, renal function, and liver function monitoring."
    )

    if score >= 3 and not labile_inr:
        recommendations.append(
            "If currently on warfarin, strongly consider switching to a DOAC "
            "for a more favourable bleeding profile."
        )

    # --- Interpretation ---
    interpretation = (
        f"HAS-BLED score: {score}/9 ({category} bleeding risk). "
        f"Estimated annual major bleeding rate: {annual_bleed_rate}%. "
        f"Components: {', '.join(detail_parts)}."
    )

    return RiskScoreResult(
        score_type=RiskScoreType.HAS_BLED,
        score_value=float(score),
        risk_category=category,
        interpretation=interpretation,
        recommendations=recommendations,
        guideline_reference=(
            "Pisters R, et al. A novel user-friendly score (HAS-BLED) to "
            "assess 1-year risk of major bleeding in patients with atrial "
            "fibrillation. Chest. 2010;138(5):1093-1100."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 5: MAGGIC HEART FAILURE RISK SCORE
# ═══════════════════════════════════════════════════════════════════════════════


def _maggic_age_points(age: int, lvef: float) -> int:
    """Calculate MAGGIC age points stratified by EF group.

    The MAGGIC scoring system assigns different age points depending on
    whether LVEF is <= 30%, 31-40%, or > 40%.

    Source: Pocock SJ, et al. Eur Heart J. 2013;34(19):1404-1413.
    """
    if lvef <= 30:
        if age < 55:
            return 0
        elif age < 60:
            return 1
        elif age < 65:
            return 2
        elif age < 70:
            return 4
        elif age < 75:
            return 6
        elif age < 80:
            return 8
        else:
            return 10
    elif lvef <= 40:
        if age < 55:
            return 0
        elif age < 60:
            return 2
        elif age < 65:
            return 4
        elif age < 70:
            return 6
        elif age < 75:
            return 8
        elif age < 80:
            return 10
        else:
            return 13
    else:
        # EF > 40% (HFpEF)
        if age < 55:
            return 0
        elif age < 60:
            return 3
        elif age < 65:
            return 5
        elif age < 70:
            return 7
        elif age < 75:
            return 9
        elif age < 80:
            return 12
        else:
            return 15


def _maggic_ef_points(lvef: float) -> int:
    """Calculate MAGGIC LVEF points. Lower EF = higher risk."""
    if lvef < 20:
        return 7
    elif lvef < 25:
        return 6
    elif lvef < 30:
        return 5
    elif lvef < 35:
        return 3
    elif lvef < 40:
        return 2
    elif lvef < 45:
        return 1
    else:
        return 0


def _maggic_sbp_points(sbp: float) -> int:
    """Calculate MAGGIC systolic BP points. Lower SBP = higher risk in HF."""
    if sbp < 110:
        return 5
    elif sbp < 120:
        return 3
    elif sbp < 130:
        return 2
    elif sbp < 140:
        return 1
    else:
        return 0


def _maggic_bmi_points(bmi: float) -> int:
    """Calculate MAGGIC BMI points. Lower BMI (cachexia) = higher risk."""
    if bmi < 15:
        return 6
    elif bmi < 20:
        return 5
    elif bmi < 25:
        return 3
    elif bmi < 30:
        return 2
    else:
        return 0


def _maggic_creatinine_points(creatinine: float) -> int:
    """Calculate MAGGIC creatinine points (creatinine in mg/dL)."""
    if creatinine < 0.9:
        return 0
    elif creatinine < 1.1:
        return 1
    elif creatinine < 1.3:
        return 2
    elif creatinine < 1.5:
        return 3
    elif creatinine < 1.7:
        return 4
    elif creatinine < 1.9:
        return 5
    elif creatinine < 2.5:
        return 6
    else:
        return 8


def _maggic_nyha_points(nyha_int: int) -> int:
    """Calculate MAGGIC NYHA class points."""
    if nyha_int == 1:
        return 0
    elif nyha_int == 2:
        return 2
    elif nyha_int == 3:
        return 6
    elif nyha_int == 4:
        return 8
    else:
        return 0


def calculate_maggic(
    inp: RiskScoreInput,
    extra: Optional[Dict] = None,
) -> RiskScoreResult:
    """Calculate MAGGIC Heart Failure Risk Score.

    Predicts 1-year and 3-year all-cause mortality in heart failure patients.

    Uses ``RiskScoreInput`` fields:
        age, sex, lvef, nyha_class, systolic_bp, bmi, creatinine,
        diabetes, smoker, beta_blocker_use, acei_arb_use.

    Optional ``extra`` fields:
        copd (bool): chronic obstructive pulmonary disease
        hf_duration_months (int): months since HF first diagnosed
        current_smoker (bool): alternative to inp.smoker

    Reference:
        Pocock SJ, et al. Predicting survival in heart failure: a risk score
        based on 39 372 patients from 30 studies. Eur Heart J.
        2013;34(19):1404-1413.
    """
    # --- Validate required inputs ---
    age = _require_int(inp.age, "age")
    sex = _require_str(inp.sex, "sex")
    lvef = _require_float(inp.lvef, "lvef")
    nyha_enum = _require(inp.nyha_class, "nyha_class")
    nyha_int = _nyha_to_int(nyha_enum)
    if nyha_int is None:
        raise RiskCalculatorError("Unable to convert NYHA class to integer.")
    sbp = _require_float(inp.systolic_bp, "systolic_bp")
    bmi = _require_float(inp.bmi, "bmi")
    creatinine = _require_float(inp.creatinine, "creatinine")
    diabetes = _require_bool(inp.diabetes, "diabetes")
    current_smoker = inp.smoker or _extra(extra, "current_smoker", False)

    is_male = sex in ("male", "m")
    beta_blocker = inp.beta_blocker_use or False
    acei_arb = inp.acei_arb_use or False
    copd = _extra(extra, "copd", False)
    hf_duration = _extra(extra, "hf_duration_months", 0)
    hf_duration_gte_18 = (hf_duration or 0) >= 18

    # --- Compute score ---
    score = 0
    detail_parts: List[str] = []

    age_pts = _maggic_age_points(age, lvef)
    score += age_pts
    detail_parts.append(f"age={age_pts}")

    ef_pts = _maggic_ef_points(lvef)
    score += ef_pts
    detail_parts.append(f"EF={ef_pts}")

    nyha_pts = _maggic_nyha_points(nyha_int)
    score += nyha_pts
    detail_parts.append(f"NYHA={nyha_pts}")

    sbp_pts = _maggic_sbp_points(sbp)
    score += sbp_pts
    detail_parts.append(f"SBP={sbp_pts}")

    bmi_pts = _maggic_bmi_points(bmi)
    score += bmi_pts
    detail_parts.append(f"BMI={bmi_pts}")

    creat_pts = _maggic_creatinine_points(creatinine)
    score += creat_pts
    detail_parts.append(f"creat={creat_pts}")

    male_pts = 1 if is_male else 0
    score += male_pts
    detail_parts.append(f"male={male_pts}")

    diabetes_pts = 3 if diabetes else 0
    score += diabetes_pts
    detail_parts.append(f"DM={diabetes_pts}")

    copd_pts = 2 if copd else 0
    score += copd_pts
    detail_parts.append(f"COPD={copd_pts}")

    # HF first diagnosed within 18 months = +2 (higher risk for new diagnosis)
    hf_new_pts = 2 if not hf_duration_gte_18 else 0
    score += hf_new_pts
    detail_parts.append(f"new_HF={hf_new_pts}")

    smoker_pts = 1 if current_smoker else 0
    score += smoker_pts
    detail_parts.append(f"smoker={smoker_pts}")

    no_bb_pts = 3 if not beta_blocker else 0
    score += no_bb_pts
    detail_parts.append(f"no_BB={no_bb_pts}")

    no_acei_pts = 1 if not acei_arb else 0
    score += no_acei_pts
    detail_parts.append(f"no_ACEi={no_acei_pts}")

    # --- Lookup mortality ---
    clamped_score = int(_clamp(score, 0, 50))
    mortality_1y = _MAGGIC_MORTALITY_1Y.get(clamped_score, 93.8)
    mortality_3y = _MAGGIC_MORTALITY_3Y.get(clamped_score, 97.0)

    # --- Categorize ---
    if mortality_1y < 10:
        category = "Low"
    elif mortality_1y < 20:
        category = "Intermediate"
    elif mortality_1y < 40:
        category = "High"
    else:
        category = "Very High"

    # --- Interpretation ---
    interpretation = (
        f"MAGGIC HF Risk Score: {score} points ({category} risk). "
        f"Predicted 1-year mortality: {mortality_1y}%. "
        f"Predicted 3-year mortality: {mortality_3y}%. "
        f"LVEF: {lvef}%, NYHA class {nyha_int}. "
        f"Point breakdown: {', '.join(detail_parts)}."
    )

    # --- Recommendations ---
    recommendations: List[str] = []

    if lvef <= 40:
        recommendations.append(
            "Heart failure with reduced ejection fraction (HFrEF). "
            "Ensure patient is on all four pillars of GDMT:"
        )
        if not acei_arb:
            recommendations.append(
                "  - START ACEi/ARB/ARNI: sacubitril-valsartan (ARNI) is "
                "preferred over ACEi/ARB (PARADIGM-HF). Titrate to target dose."
            )
        else:
            recommendations.append(
                "  - ACEi/ARB/ARNI: already prescribed. Ensure at target dose. "
                "Consider switch to ARNI if on ACEi/ARB and tolerating."
            )

        if not beta_blocker:
            recommendations.append(
                "  - START beta-blocker: carvedilol, bisoprolol, or "
                "metoprolol succinate. Titrate to target dose."
            )
        else:
            recommendations.append(
                "  - Beta-blocker: already prescribed. Ensure at target dose."
            )

        recommendations.append(
            "  - Mineralocorticoid receptor antagonist (MRA): spironolactone "
            "or eplerenone if K+ < 5.0 and eGFR > 30."
        )
        recommendations.append(
            "  - SGLT2 inhibitor: dapagliflozin or empagliflozin "
            "(benefit independent of diabetes status)."
        )

        if lvef <= 35 and nyha_int >= 2:
            recommendations.append(
                "With LVEF <= 35% and NYHA >= II on optimal GDMT for >= 3 "
                "months: evaluate for ICD (primary prevention of sudden "
                "cardiac death)."
            )
        if lvef <= 35 and nyha_int >= 2 and score >= 20:
            recommendations.append(
                "With persistent symptoms despite optimal GDMT: consider CRT "
                "(if QRS >= 130 ms with LBBB pattern) or advanced HF therapies "
                "(LVAD, heart transplant evaluation)."
            )
    else:
        recommendations.append(
            "Heart failure with preserved ejection fraction (HFpEF). "
            "Evidence-based therapies:"
        )
        recommendations.append(
            "  - SGLT2 inhibitor: empagliflozin (EMPEROR-Preserved) or "
            "dapagliflozin (DELIVER trial)."
        )
        recommendations.append(
            "  - Diuretics for congestion management."
        )
        recommendations.append(
            "  - Treat underlying conditions: hypertension, AF, CAD, obesity."
        )

    if mortality_1y >= 20:
        recommendations.append(
            f"HIGH mortality risk ({mortality_1y}% at 1 year). "
            "Consider referral to advanced heart failure / transplant centre."
        )

    if current_smoker:
        recommendations.append(
            "Smoking cessation: critical for cardiovascular risk reduction."
        )

    if diabetes:
        recommendations.append(
            "Optimise glycaemic control. SGLT2 inhibitors have dual benefit "
            "for HF and diabetes."
        )

    if copd:
        recommendations.append(
            "COPD management: cardioselective beta-blockers (bisoprolol) "
            "are safe and should not be withheld."
        )

    return RiskScoreResult(
        score_type=RiskScoreType.MAGGIC,
        score_value=float(score),
        risk_category=category,
        interpretation=interpretation,
        recommendations=recommendations,
        guideline_reference=(
            "Pocock SJ, et al. Predicting survival in heart failure: a risk "
            "score based on 39 372 patients from 30 studies. Eur Heart J. "
            "2013;34(19):1404-1413. "
            "Also: 2022 AHA/ACC/HFSA Guideline for Management of Heart "
            "Failure (Heidenreich PA, et al. J Am Coll Cardiol. "
            "2022;79(17):e263-e421)."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# CALCULATOR 6: EuroSCORE II
# ═══════════════════════════════════════════════════════════════════════════════


def _euroscore_age_term(age: int) -> float:
    """Compute the EuroSCORE II age contribution.

    For age <= 60, contribution is 0 (reference).
    For age > 60, contribution = coefficient * (age - 60).
    """
    if age <= 60:
        return 0.0
    return _EUROSCORE_II_COEFFICIENTS["age"] * (age - 60)


def _euroscore_renal_term(
    creatinine_clearance: Optional[float],
    on_dialysis: bool = False,
) -> float:
    """Compute the EuroSCORE II renal impairment term.

    Categories:
        CrCl > 85 ml/min         : 0 (reference)
        CrCl 50-85 ml/min        : moderate impairment
        CrCl < 50 ml/min         : severe impairment (not on dialysis)
        On dialysis               : highest risk
    """
    if on_dialysis:
        return _EUROSCORE_II_COEFFICIENTS["creat_on_dialysis"]
    if creatinine_clearance is None:
        return 0.0
    if creatinine_clearance < 50:
        return _EUROSCORE_II_COEFFICIENTS["creat_cc_lt_50"]
    elif creatinine_clearance <= 85:
        return _EUROSCORE_II_COEFFICIENTS["creat_cc_50_85"]
    return 0.0


def _euroscore_lvef_term(lvef: Optional[float]) -> float:
    """Compute the EuroSCORE II LVEF term.

    LVEF > 50% : 0 (good), 31-50%: moderate, 21-30%: poor, <=20%: very poor.
    """
    if lvef is None:
        return 0.0
    if lvef <= 20:
        return _EUROSCORE_II_COEFFICIENTS["lvef_le_20"]
    elif lvef <= 30:
        return _EUROSCORE_II_COEFFICIENTS["lvef_21_30"]
    elif lvef <= 50:
        return _EUROSCORE_II_COEFFICIENTS["lvef_31_50"]
    return 0.0


def _euroscore_nyha_term(nyha_int: Optional[int]) -> float:
    """Compute the EuroSCORE II NYHA class term."""
    if nyha_int is None or nyha_int <= 1:
        return 0.0
    elif nyha_int == 2:
        return _EUROSCORE_II_COEFFICIENTS["nyha_class_ii"]
    elif nyha_int == 3:
        return _EUROSCORE_II_COEFFICIENTS["nyha_class_iii"]
    elif nyha_int >= 4:
        return _EUROSCORE_II_COEFFICIENTS["nyha_class_iv"]
    return 0.0


def _euroscore_pulm_hyp_term(pulm_hyp: Optional[str]) -> float:
    """Compute the EuroSCORE II pulmonary hypertension term.

    'none' / None -> 0, 'moderate' -> 31-55 mmHg, 'severe' -> > 55 mmHg.
    """
    if not pulm_hyp:
        return 0.0
    ph = pulm_hyp.strip().lower()
    if ph == "moderate":
        return _EUROSCORE_II_COEFFICIENTS["pulm_hyp_moderate"]
    elif ph == "severe":
        return _EUROSCORE_II_COEFFICIENTS["pulm_hyp_severe"]
    return 0.0


def _euroscore_urgency_term(urgency: Optional[str]) -> float:
    """Compute the EuroSCORE II surgery urgency term.

    'elective' -> 0, 'urgent', 'emergency', 'salvage'.
    """
    if not urgency:
        return 0.0
    u = urgency.strip().lower()
    if u == "urgent":
        return _EUROSCORE_II_COEFFICIENTS["urgency_urgent"]
    elif u in ("emergency", "emergent"):
        return _EUROSCORE_II_COEFFICIENTS["urgency_emergency"]
    elif u == "salvage":
        return _EUROSCORE_II_COEFFICIENTS["urgency_salvage"]
    return 0.0


def _euroscore_weight_term(weight_of_procedure: Optional[str]) -> float:
    """Compute the EuroSCORE II weight of intervention term.

    'isolated_cabg' -> 0, 'single_non_cabg', 'two_procedures', 'three_or_more'.
    """
    if not weight_of_procedure:
        return 0.0
    w = weight_of_procedure.strip().lower()
    if w == "single_non_cabg":
        return _EUROSCORE_II_COEFFICIENTS["weight_single_non_cabg"]
    elif w == "two_procedures":
        return _EUROSCORE_II_COEFFICIENTS["weight_two_procedures"]
    elif w in ("three_or_more", "three_procedures"):
        return _EUROSCORE_II_COEFFICIENTS["weight_three_or_more"]
    return 0.0


def calculate_euroscore_ii(
    inp: RiskScoreInput,
    extra: Optional[Dict] = None,
) -> RiskScoreResult:
    """Calculate EuroSCORE II predicted operative mortality.

    Uses the published logistic regression model:
        predicted_mortality = exp(Z) / (1 + exp(Z))
    where Z = beta_0 + sum(beta_i * x_i).

    Uses ``RiskScoreInput`` fields:
        age, sex, lvef, nyha_class, urgency, redo_surgery.

    Additional fields from ``extra`` dict (all optional, default False/None):
        creatinine_clearance_ml_min (float): Cockcroft-Gault CrCl
        on_dialysis (bool): currently on dialysis
        extracardiac_arteriopathy (bool): claudication, carotid stenosis, etc.
        poor_mobility (bool): musculoskeletal or neurological impairment
        chronic_lung_disease (bool): long-term bronchodilator/steroid use
        active_endocarditis (bool): on antibiotics at time of surgery
        critical_preoperative_state (bool): VT/VF, preop cardiac massage, etc.
        diabetes_on_insulin (bool)
        ccs_class_4_angina (bool)
        recent_mi (bool): MI within 90 days
        pulmonary_hypertension (str): 'none', 'moderate', 'severe'
        weight_of_procedure (str): 'isolated_cabg', 'single_non_cabg',
            'two_procedures', 'three_or_more'
        surgery_on_thoracic_aorta (bool)

    Reference:
        Nashef SA, et al. EuroSCORE II. Eur J Cardiothorac Surg.
        2012;41(4):734-745.
    """
    # --- Validate required inputs ---
    age = _require_int(inp.age, "age")
    sex = _require_str(inp.sex, "sex")
    is_female = sex in ("female", "f")

    # Optional fields from inp
    nyha_int = _nyha_to_int(inp.nyha_class)
    lvef = inp.lvef
    urgency = inp.urgency
    previous_cardiac_surgery = inp.redo_surgery or False

    # Optional fields from extra
    creatinine_clearance = _extra(extra, "creatinine_clearance_ml_min")
    on_dialysis = _extra(extra, "on_dialysis", False)
    extracardiac_arteriopathy = _extra(extra, "extracardiac_arteriopathy", False)
    poor_mobility = _extra(extra, "poor_mobility", False)
    chronic_lung_disease = _extra(extra, "chronic_lung_disease", False)
    active_endocarditis = _extra(extra, "active_endocarditis", False)
    critical_preop = _extra(extra, "critical_preoperative_state", False)
    diabetes_insulin = _extra(extra, "diabetes_on_insulin", False)
    ccs_4 = _extra(extra, "ccs_class_4_angina", False)
    recent_mi = _extra(extra, "recent_mi", False)
    pulm_hyp = _extra(extra, "pulmonary_hypertension")
    weight_proc = _extra(extra, "weight_of_procedure")
    thoracic_aorta = _extra(extra, "surgery_on_thoracic_aorta", False)

    # --- Compute linear predictor Z ---
    z = _EUROSCORE_II_BETA0
    detail_parts: List[str] = []

    # Age
    age_term = _euroscore_age_term(age)
    z += age_term
    if age_term > 0:
        detail_parts.append(f"age(+{age_term:.3f})")

    # Female sex
    female_term = _EUROSCORE_II_COEFFICIENTS["female"] if is_female else 0.0
    z += female_term
    if female_term > 0:
        detail_parts.append(f"female(+{female_term:.3f})")

    # Renal function
    renal_term = _euroscore_renal_term(creatinine_clearance, on_dialysis)
    z += renal_term
    if renal_term > 0:
        detail_parts.append(f"renal(+{renal_term:.3f})")

    # Extracardiac arteriopathy
    eca_term = (
        _EUROSCORE_II_COEFFICIENTS["extracardiac_arteriopathy"]
        if extracardiac_arteriopathy else 0.0
    )
    z += eca_term
    if eca_term > 0:
        detail_parts.append(f"ECA(+{eca_term:.3f})")

    # Poor mobility
    mob_term = (
        _EUROSCORE_II_COEFFICIENTS["poor_mobility"]
        if poor_mobility else 0.0
    )
    z += mob_term
    if mob_term > 0:
        detail_parts.append(f"mobility(+{mob_term:.3f})")

    # Previous cardiac surgery
    prev_term = (
        _EUROSCORE_II_COEFFICIENTS["previous_cardiac_surgery"]
        if previous_cardiac_surgery else 0.0
    )
    z += prev_term
    if prev_term > 0:
        detail_parts.append(f"redo(+{prev_term:.3f})")

    # Chronic lung disease
    cld_term = (
        _EUROSCORE_II_COEFFICIENTS["chronic_lung_disease"]
        if chronic_lung_disease else 0.0
    )
    z += cld_term
    if cld_term > 0:
        detail_parts.append(f"CLD(+{cld_term:.3f})")

    # Active endocarditis
    endo_term = (
        _EUROSCORE_II_COEFFICIENTS["active_endocarditis"]
        if active_endocarditis else 0.0
    )
    z += endo_term
    if endo_term > 0:
        detail_parts.append(f"endocarditis(+{endo_term:.3f})")

    # Critical preoperative state
    crit_term = (
        _EUROSCORE_II_COEFFICIENTS["critical_preoperative_state"]
        if critical_preop else 0.0
    )
    z += crit_term
    if crit_term > 0:
        detail_parts.append(f"critical(+{crit_term:.3f})")

    # Diabetes on insulin
    diab_term = (
        _EUROSCORE_II_COEFFICIENTS["diabetes_on_insulin"]
        if diabetes_insulin else 0.0
    )
    z += diab_term
    if diab_term > 0:
        detail_parts.append(f"DM-insulin(+{diab_term:.3f})")

    # NYHA class
    nyha_term = _euroscore_nyha_term(nyha_int)
    z += nyha_term
    if nyha_term > 0:
        detail_parts.append(f"NYHA(+{nyha_term:.3f})")

    # CCS class 4 angina
    ccs_term = (
        _EUROSCORE_II_COEFFICIENTS["ccs_class_4_angina"]
        if ccs_4 else 0.0
    )
    z += ccs_term
    if ccs_term > 0:
        detail_parts.append(f"CCS4(+{ccs_term:.3f})")

    # LVEF
    lvef_term = _euroscore_lvef_term(lvef)
    z += lvef_term
    if lvef_term > 0:
        detail_parts.append(f"LVEF(+{lvef_term:.3f})")

    # Recent MI
    mi_term = (
        _EUROSCORE_II_COEFFICIENTS["recent_mi"]
        if recent_mi else 0.0
    )
    z += mi_term
    if mi_term > 0:
        detail_parts.append(f"MI(+{mi_term:.3f})")

    # Pulmonary hypertension
    ph_term = _euroscore_pulm_hyp_term(pulm_hyp)
    z += ph_term
    if ph_term > 0:
        detail_parts.append(f"PH(+{ph_term:.3f})")

    # Urgency
    urg_term = _euroscore_urgency_term(urgency)
    z += urg_term
    if urg_term > 0:
        detail_parts.append(f"urgency(+{urg_term:.3f})")

    # Weight of procedure
    wt_term = _euroscore_weight_term(weight_proc)
    z += wt_term
    if wt_term > 0:
        detail_parts.append(f"weight(+{wt_term:.3f})")

    # Surgery on thoracic aorta
    aorta_term = (
        _EUROSCORE_II_COEFFICIENTS["surgery_thoracic_aorta"]
        if thoracic_aorta else 0.0
    )
    z += aorta_term
    if aorta_term > 0:
        detail_parts.append(f"aorta(+{aorta_term:.3f})")

    # --- Logistic transformation ---
    try:
        exp_z = math.exp(z)
        predicted_mortality = exp_z / (1.0 + exp_z)
    except OverflowError:
        predicted_mortality = 1.0

    mortality_pct = _clamp(predicted_mortality * 100.0, 0.0, 100.0)

    # --- Categorize ---
    if mortality_pct < 1.0:
        category = "Low"
    elif mortality_pct < 3.0:
        category = "Low-Moderate"
    elif mortality_pct < 5.0:
        category = "Moderate"
    elif mortality_pct < 10.0:
        category = "High"
    else:
        category = "Very High"

    # --- Interpretation ---
    contributor_str = (
        f" Key contributors: {', '.join(detail_parts)}." if detail_parts else ""
    )
    interpretation = (
        f"EuroSCORE II: predicted operative mortality {mortality_pct:.2f}% "
        f"({category} risk). "
        f"Age {age}, {'female' if is_female else 'male'}. "
        f"Linear predictor Z = {z:.4f}.{contributor_str}"
    )

    # --- Recommendations ---
    recommendations: List[str] = []

    recommendations.append(
        f"Predicted operative mortality: {mortality_pct:.2f}%. "
        f"Discuss risk with the patient and surgical team as part of "
        f"shared decision-making."
    )

    if mortality_pct < 3.0:
        recommendations.append(
            "Low operative risk. Standard perioperative management appropriate."
        )
    elif mortality_pct < 5.0:
        recommendations.append(
            "Moderate operative risk. Ensure multidisciplinary Heart Team "
            "discussion of surgical vs. interventional options (if applicable)."
        )
    elif mortality_pct < 10.0:
        recommendations.append(
            "High operative risk. Mandatory Heart Team discussion. Consider "
            "transcatheter alternatives (TAVI for aortic stenosis, MitraClip "
            "for mitral regurgitation) where applicable."
        )
        recommendations.append(
            "Optimise modifiable risk factors pre-operatively: nutrition, "
            "pulmonary rehabilitation, glycaemic control, renal function."
        )
    else:
        recommendations.append(
            "Very high operative risk (>= 10% predicted mortality). "
            "Heart Team discussion is essential. Strongly consider "
            "transcatheter or conservative alternatives."
        )
        recommendations.append(
            "If surgery is pursued, ensure ICU bed availability, "
            "consider staged procedures, and plan for possible mechanical "
            "circulatory support."
        )
        recommendations.append(
            "Goals-of-care discussion with patient and family."
        )

    if previous_cardiac_surgery:
        recommendations.append(
            "Redo cardiac surgery increases risk substantially. "
            "CT chest for sternal re-entry planning is recommended."
        )

    if critical_preop:
        recommendations.append(
            "Critical preoperative state identified. Ensure haemodynamic "
            "stabilisation prior to OR. Consider mechanical circulatory "
            "support (IABP, ECMO) if haemodynamically unstable."
        )

    if active_endocarditis:
        recommendations.append(
            "Active endocarditis: continue IV antibiotics. Timing of "
            "surgery per infectious disease and cardiac surgery consensus."
        )

    if urgency and urgency.strip().lower() in ("emergency", "emergent", "salvage"):
        recommendations.append(
            "Emergency/salvage procedure: abbreviated risk discussion. "
            "Ensure blood products are available and perfusion team is "
            "on standby."
        )

    return RiskScoreResult(
        score_type=RiskScoreType.EUROSCORE_II,
        score_value=round(mortality_pct, 2),
        risk_category=category,
        interpretation=interpretation,
        recommendations=recommendations,
        guideline_reference=(
            "Nashef SA, et al. EuroSCORE II. Eur J Cardiothorac Surg. "
            "2012;41(4):734-745. "
            "Also: 2021 ESC/EACTS Guidelines for Valvular Heart Disease "
            "(Vahanian A, et al. Eur Heart J. 2022;43(7):561-632)."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════════════
# RISK CALCULATOR ENGINE
# ═══════════════════════════════════════════════════════════════════════════════


class RiskCalculatorEngine:
    """Unified dispatcher for all cardiovascular risk calculators.

    Usage::

        engine = RiskCalculatorEngine()

        # Single calculator via RiskScoreInput
        result = engine.calculate(RiskScoreInput(
            score_type=RiskScoreType.ASCVD,
            age=55, sex="male", race="white",
            total_cholesterol=213, hdl=50,
            systolic_bp=120, hypertension_treatment=False,
            diabetes=False, smoker=True,
        ))

        # All applicable calculators for a patient
        results = engine.calculate_all_applicable({
            "age": 68, "sex": "male", "race": "white",
            "total_cholesterol": 220, "hdl": 45,
            "systolic_bp": 148, "hypertension_treatment": True,
            "diabetes": True, "smoker": False,
            "congestive_heart_failure": True,
            "lvef": 35, "nyha_class": "nyha_iii",
            "bmi": 29.5, "creatinine": 1.4,
        })
    """

    # Map score type to its calculator function
    _CALCULATORS = {
        RiskScoreType.ASCVD: calculate_ascvd,
        RiskScoreType.HEART: calculate_heart_score,
        RiskScoreType.CHA2DS2_VASC: calculate_cha2ds2_vasc,
        RiskScoreType.HAS_BLED: calculate_has_bled,
        RiskScoreType.MAGGIC: calculate_maggic,
        RiskScoreType.EUROSCORE_II: calculate_euroscore_ii,
    }

    # Minimum required fields per calculator (for applicability check).
    # These are RiskScoreInput field names.
    _REQUIRED_FIELDS: Dict[RiskScoreType, List[str]] = {
        RiskScoreType.ASCVD: [
            "age", "sex", "total_cholesterol", "hdl",
            "systolic_bp", "hypertension_treatment", "diabetes", "smoker",
        ],
        RiskScoreType.HEART: [
            "age",
        ],
        RiskScoreType.CHA2DS2_VASC: [
            "age", "sex",
        ],
        RiskScoreType.HAS_BLED: [
            "age",
        ],
        RiskScoreType.MAGGIC: [
            "age", "sex", "lvef", "nyha_class",
            "systolic_bp", "bmi", "creatinine", "diabetes",
        ],
        RiskScoreType.EUROSCORE_II: [
            "age", "sex",
        ],
    }

    def calculate(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate risk using the specified score type.

        Args:
            inp: Validated RiskScoreInput with score_type set.
            extra: Optional dictionary of calculator-specific fields not
                present on RiskScoreInput (e.g., HEART score
                history_suspicion, EuroSCORE II operative details).

        Returns:
            RiskScoreResult with score, category, interpretation,
            recommendations, and guideline reference.

        Raises:
            RiskCalculatorError: If required fields are missing or
                values are out of range.
            ValueError: If score_type is not supported.
        """
        calculator = self._CALCULATORS.get(inp.score_type)
        if calculator is None:
            raise ValueError(
                f"Unsupported risk score type: {inp.score_type}. "
                f"Supported types: {list(self._CALCULATORS.keys())}"
            )
        return calculator(inp, extra)

    def calculate_ascvd(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate 10-year ASCVD risk (Pooled Cohort Equations)."""
        inp.score_type = RiskScoreType.ASCVD
        return calculate_ascvd(inp, extra)

    def calculate_heart(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate HEART Score for ED chest pain triage."""
        inp.score_type = RiskScoreType.HEART
        return calculate_heart_score(inp, extra)

    def calculate_cha2ds2_vasc(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate CHA2DS2-VASc score for AF stroke risk."""
        inp.score_type = RiskScoreType.CHA2DS2_VASC
        return calculate_cha2ds2_vasc(inp, extra)

    def calculate_has_bled(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate HAS-BLED bleeding risk score."""
        inp.score_type = RiskScoreType.HAS_BLED
        return calculate_has_bled(inp, extra)

    def calculate_maggic(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate MAGGIC heart failure risk score."""
        inp.score_type = RiskScoreType.MAGGIC
        return calculate_maggic(inp, extra)

    def calculate_euroscore_ii(
        self,
        inp: RiskScoreInput,
        extra: Optional[Dict] = None,
    ) -> RiskScoreResult:
        """Calculate EuroSCORE II operative mortality prediction."""
        inp.score_type = RiskScoreType.EUROSCORE_II
        return calculate_euroscore_ii(inp, extra)

    def _is_applicable(
        self,
        score_type: RiskScoreType,
        patient_data: Dict,
    ) -> bool:
        """Check whether a calculator can run given the available data.

        A calculator is considered applicable if all its required fields
        are present (not None) in the patient data dictionary.

        Additional clinical eligibility:
          - ASCVD: age must be 40-79.
        """
        required = self._REQUIRED_FIELDS.get(score_type, [])
        for field in required:
            if field not in patient_data or patient_data[field] is None:
                return False

        if score_type == RiskScoreType.ASCVD:
            age = patient_data.get("age")
            if age is not None and (age < 40 or age > 79):
                return False

        return True

    def calculate_all_applicable(
        self,
        patient_data: Dict,
        extra: Optional[Dict] = None,
    ) -> List[RiskScoreResult]:
        """Run every applicable risk calculator for the given patient data.

        Takes a flat dictionary of patient parameters (matching
        RiskScoreInput field names) and attempts each calculator where
        sufficient data is available.

        Args:
            patient_data: Dictionary of patient parameters. Keys should
                match RiskScoreInput field names. If ``nyha_class`` is
                provided as an integer (1-4), it will be converted to the
                HeartFailureClass enum automatically.
            extra: Optional dictionary of calculator-specific fields not
                present on RiskScoreInput.

        Returns:
            List of RiskScoreResult objects for all applicable and
            successfully computed calculators. Calculators that fail
            due to data issues are silently skipped.
        """
        results: List[RiskScoreResult] = []

        # Pre-process: convert nyha_class int -> HeartFailureClass enum
        data = dict(patient_data)
        nyha_val = data.get("nyha_class")
        if isinstance(nyha_val, int):
            nyha_map = {
                1: HeartFailureClass.NYHA_I,
                2: HeartFailureClass.NYHA_II,
                3: HeartFailureClass.NYHA_III,
                4: HeartFailureClass.NYHA_IV,
            }
            data["nyha_class"] = nyha_map.get(nyha_val)
        elif isinstance(nyha_val, str) and nyha_val.startswith("nyha_"):
            try:
                data["nyha_class"] = HeartFailureClass(nyha_val)
            except ValueError:
                pass

        for score_type in RiskScoreType:
            if not self._is_applicable(score_type, data):
                continue

            try:
                data_with_type = {**data, "score_type": score_type.value}
                inp = RiskScoreInput(**data_with_type)
                result = self.calculate(inp, extra)
                results.append(result)
            except (RiskCalculatorError, ValueError, TypeError):
                continue

        return results

    @staticmethod
    def get_supported_calculators() -> List[Dict[str, str]]:
        """Return metadata for all supported risk calculators."""
        return [
            {
                "id": RiskScoreType.ASCVD.value,
                "name": "ASCVD Pooled Cohort Equations",
                "description": (
                    "10-year atherosclerotic cardiovascular disease risk. "
                    "Validated for adults aged 40-79."
                ),
                "reference": (
                    "Goff DC Jr, et al. J Am Coll Cardiol. "
                    "2014;63(25 Pt B):2935-2959"
                ),
                "max_score": "100 (percentage)",
            },
            {
                "id": RiskScoreType.HEART.value,
                "name": "HEART Score",
                "description": (
                    "Emergency department chest pain risk assessment. "
                    "Predicts 6-week MACE (death, MI, revascularization)."
                ),
                "reference": (
                    "Six AJ, et al. Neth Heart J. 2008;16(6):191-196"
                ),
                "max_score": "10",
            },
            {
                "id": RiskScoreType.CHA2DS2_VASC.value,
                "name": "CHA2DS2-VASc",
                "description": (
                    "Stroke risk in non-valvular atrial fibrillation. "
                    "Guides anticoagulation decisions."
                ),
                "reference": (
                    "Lip GYH, et al. Chest. 2010;137(2):263-272"
                ),
                "max_score": "9",
            },
            {
                "id": RiskScoreType.HAS_BLED.value,
                "name": "HAS-BLED",
                "description": (
                    "Bleeding risk on anticoagulation. "
                    "Identifies modifiable bleeding risk factors."
                ),
                "reference": (
                    "Pisters R, et al. Chest. 2010;138(5):1093-1100"
                ),
                "max_score": "9",
            },
            {
                "id": RiskScoreType.MAGGIC.value,
                "name": "MAGGIC Heart Failure Risk Score",
                "description": (
                    "1-year and 3-year all-cause mortality in heart failure. "
                    "Based on the MAGGIC meta-analysis of 30 studies."
                ),
                "reference": (
                    "Pocock SJ, et al. Eur Heart J. 2013;34(19):1404-1413"
                ),
                "max_score": "50",
            },
            {
                "id": RiskScoreType.EUROSCORE_II.value,
                "name": "EuroSCORE II",
                "description": (
                    "Predicted operative mortality for cardiac surgery. "
                    "Uses logistic regression with 18 risk factors."
                ),
                "reference": (
                    "Nashef SA, et al. Eur J Cardiothorac Surg. "
                    "2012;41(4):734-745"
                ),
                "max_score": "100 (percentage)",
            },
        ]

    @staticmethod
    def get_required_fields(score_type: RiskScoreType) -> List[str]:
        """Return the list of required input fields for a given calculator.

        Useful for building dynamic forms or validating API requests
        before invoking the calculator.
        """
        fields_map: Dict[RiskScoreType, List[str]] = {
            RiskScoreType.ASCVD: [
                "age (int, 40-79)",
                "sex ('male' or 'female')",
                "race ('white', 'african_american', or 'other')",
                "total_cholesterol (float, mg/dL)",
                "hdl (float, mg/dL)",
                "systolic_bp (float, mmHg)",
                "hypertension_treatment (bool)",
                "diabetes (bool)",
                "smoker (bool)",
            ],
            RiskScoreType.HEART: [
                "age (int)",
                "troponin (float, ng/mL; or troponin_ratio in extra)",
                "extra.history_suspicion (int 0-2)",
                "extra.ecg_finding (int 0-2)",
                "extra.num_risk_factors (int)",
                "extra.history_atherosclerosis (bool, optional)",
                "extra.troponin_uln (float, default 0.04)",
            ],
            RiskScoreType.CHA2DS2_VASC: [
                "age (int)",
                "sex ('male' or 'female')",
                "congestive_heart_failure (bool)",
                "hypertension_treatment (bool, proxy for hypertension)",
                "diabetes (bool)",
                "history_of_stroke (bool)",
                "vascular_disease (bool)",
            ],
            RiskScoreType.HAS_BLED: [
                "age (int)",
                "systolic_bp (float, optional -- > 160 = uncontrolled HTN)",
                "renal_disease (bool)",
                "liver_disease (bool)",
                "history_of_stroke (bool)",
                "history_of_bleeding (bool)",
                "labile_inr (bool)",
                "antiplatelet_nsaid (bool)",
                "alcohol_excess (bool)",
            ],
            RiskScoreType.MAGGIC: [
                "age (int)",
                "sex ('male' or 'female')",
                "lvef (float, %)",
                "nyha_class (HeartFailureClass or int 1-4 via engine)",
                "systolic_bp (float, mmHg)",
                "bmi (float, kg/m2)",
                "creatinine (float, mg/dL)",
                "diabetes (bool)",
                "smoker (bool)",
                "beta_blocker_use (bool)",
                "acei_arb_use (bool)",
                "extra.copd (bool, optional)",
                "extra.hf_duration_months (int, optional)",
            ],
            RiskScoreType.EUROSCORE_II: [
                "age (int)",
                "sex ('male' or 'female')",
                "lvef (float, %, optional)",
                "nyha_class (HeartFailureClass, optional)",
                "urgency ('elective', 'urgent', 'emergent', 'salvage')",
                "redo_surgery (bool)",
                "extra.creatinine_clearance_ml_min (float)",
                "extra.extracardiac_arteriopathy (bool)",
                "extra.poor_mobility (bool)",
                "extra.chronic_lung_disease (bool)",
                "extra.active_endocarditis (bool)",
                "extra.critical_preoperative_state (bool)",
                "extra.diabetes_on_insulin (bool)",
                "extra.ccs_class_4_angina (bool)",
                "extra.recent_mi (bool)",
                "extra.pulmonary_hypertension ('none'/'moderate'/'severe')",
                "extra.weight_of_procedure (str)",
                "extra.surgery_on_thoracic_aorta (bool)",
            ],
        }
        return fields_map.get(score_type, [])
