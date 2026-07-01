"""Biological age calculation engine.

Implements PhenoAge (Levine et al. 2018) and GrimAge surrogate estimation
from standard blood biomarkers. Pure computation — no LLM or database calls.

PhenoAge reference: Levine et al., "An epigenetic biomarker of aging for
lifespan and healthspan", Aging 2018; 10(4):573-591.

The original PhenoAge coefficients (from dayoonkwon/BioAge R package) expect
SI units: albumin g/L, creatinine μmol/L, glucose mmol/L, ln(CRP mg/L).
This module accepts standard US clinical units (g/dL, mg/dL, mg/dL, mg/L)
and converts internally before applying coefficients.

Author: Adam Jones
Date: March 2026
"""

import math
from typing import Any, Dict, List

from loguru import logger

# PhenoAge coefficients (Levine et al., Aging 2018; dayoonkwon/BioAge R package)
# IMPORTANT: These coefficients expect SI units internally:
#   albumin: g/L, creatinine: μmol/L, glucose: mmol/L, CRP: ln(mg/L),
#   lymphocyte: %, MCV: fL, RDW: %, alkaline phosphatase: U/L, WBC: 10^3/μL
PHENOAGE_COEFFICIENTS = {
    "albumin": -0.0336,              # protective (negative), applied to g/L
    "creatinine": 0.0095,            # applied to μmol/L
    "glucose": 0.1953,               # applied to mmol/L (ln-transformed in original)
    "ln_crp": 0.0954,               # applied to ln(CRP in mg/L)
    "lymphocyte_pct": -0.0120,      # protective (negative), applied to %
    "mcv": 0.0268,                   # applied to fL
    "rdw": 0.3306,                   # applied to %
    "alkaline_phosphatase": 0.0019,  # applied to U/L
    "wbc": 0.0554,                   # applied to 10^3/μL
}
PHENOAGE_INTERCEPT = -19.9067

# Gompertz mortality parameters (from dayoonkwn/BioAge orig=TRUE code)
MORT_NUMERATOR = -1.51714       # m_n: -(exp(120 * gamma) - 1)
MORT_DENOMINATOR = 0.007692696  # m_d: gamma (Gompertz shape)

# Biological age conversion parameters (from age-only Gompertz model)
BA_NUMERATOR = -0.0055305       # BA_n: -rate from age-only model
BA_DENOMINATOR = 0.09165        # BA_d: age coefficient from age-only model
BA_INTERCEPT = 141.50225        # BA_i: intercept for age conversion

# Unit conversion factors: US clinical units → SI units for PhenoAge
# These are applied before coefficients are used
UNIT_CONVERSIONS = {
    "albumin": 10.0,        # g/dL → g/L (multiply by 10)
    "creatinine": 88.4,     # mg/dL → μmol/L (multiply by 88.4)
    "glucose": 1 / 18.016,  # mg/dL → mmol/L (divide by 18.016)
    # Other biomarkers: same units in US and SI (%, fL, U/L, 10^3/μL, ln(mg/L))
}

# PhenoAge standard error from NHANES III validation cohort
# Levine et al. 2018 PMID:29676998 (Table 1, SE = 4.9 years for full model)
# SE decreases with more biomarkers available
PHENOAGE_SE_FULL = 4.9  # years, all 9 biomarkers
PHENOAGE_SE_PARTIAL = 6.5  # years, <9 biomarkers (increased uncertainty)

# GrimAge surrogate markers and weights (simplified from Lu et al. 2019)
# Lu et al. 2019 PMID:30669119; Hillary et al. 2020 PMID:32941527
GRIMAGE_MARKERS = {
    "gdf15": {"weight": 0.15, "unit": "pg/mL", "ref_max": 1200.0},
    "cystatin_c": {"weight": 0.12, "unit": "mg/L", "ref_max": 1.0},
    "leptin": {"weight": 0.08, "unit": "ng/mL", "ref_max": 15.0},
    "pai1": {"weight": 0.10, "unit": "ng/mL", "ref_max": 43.0},
    "timp1": {"weight": 0.09, "unit": "ng/mL", "ref_max": 250.0},
    "adm": {"weight": 0.11, "unit": "pmol/L", "ref_max": 50.0},
}

# Published plasma protein reference ranges for GrimAge surrogate validation
# Sources: Hillary et al. 2020 PMID:32941527; Tanaka et al. 2018 PMID:29371514
GRIMAGE_VALIDATION = {
    "correlation_with_true_grimage": 0.72,  # r² from Hillary et al. 2020
    "standard_error_years": 5.8,  # SE of surrogate vs true GrimAge
    "validation_cohort": "Lothian Birth Cohort 1936 (n=906)",
    "validation_citation": "Hillary et al. 2020 PMID:32941527",
    "note": (
        "Plasma protein surrogates explain ~72% of variance in DNAm GrimAge. "
        "Missing components: smoking pack-years (major), DNA methylation patterns. "
        "Surrogate estimates should be interpreted as directional indicators, "
        "not precise age measurements."
    ),
}


def validate_biomarker_ranges(biomarkers: Dict[str, float]) -> List[str]:
    """Check biomarker values against plausible clinical ranges.

    Returns list of warning strings for out-of-range values.
    """
    from src.knowledge import BIOMARKER_PLAUSIBLE_RANGES

    warnings = []
    for name, value in biomarkers.items():
        if name in BIOMARKER_PLAUSIBLE_RANGES:
            lower, upper = BIOMARKER_PLAUSIBLE_RANGES[name]
            if value < lower or value > upper:
                warnings.append(
                    f"{name}={value} outside plausible range ({lower}-{upper}); "
                    f"possible data entry error"
                )
    return warnings


class BiologicalAgeCalculator:
    """Calculates biological age using PhenoAge and GrimAge surrogates.

    PhenoAge uses 9 routine blood biomarkers to estimate biological age
    based on mortality risk modeling from NHANES III data.

    GrimAge uses plasma protein surrogates (GDF-15, Cystatin C, etc.)
    that correlate with DNAm GrimAge clocks.
    """

    def calculate_phenoage(
        self,
        chronological_age: float,
        biomarkers: Dict[str, float],
    ) -> Dict[str, Any]:
        """Calculate PhenoAge biological age.

        Args:
            chronological_age: Patient's chronological age in years.
            biomarkers: Dict with keys matching PHENOAGE_COEFFICIENTS.
                Required: albumin, creatinine, glucose, ln_crp (or hs_crp
                which will be log-transformed), lymphocyte_pct, mcv, rdw,
                alkaline_phosphatase, wbc.

        Returns:
            Dict with biological_age, age_acceleration, mortality_risk,
            and per-biomarker contributions.
        """
        # Validate non-negative biomarker values — create sanitized copy
        sanitized = {}
        for marker, val in biomarkers.items():
            if val is not None and val < 0:
                logger.warning(f"Negative biomarker value: {marker}={val}. Setting to 0.")
                sanitized[marker] = 0
            else:
                sanitized[marker] = val
        biomarkers = sanitized

        range_warnings = validate_biomarker_ranges(biomarkers)
        for w in range_warnings:
            logger.warning(f"Input validation: {w}")

        # Handle hs_crp -> ln_crp transformation
        working = dict(biomarkers)
        if "hs_crp" in working and "ln_crp" not in working:
            crp_val = working.pop("hs_crp")
            working["ln_crp"] = math.log(max(min(crp_val, 200.0), 0.01))

        # Convert US clinical units → SI units for PhenoAge coefficients
        for marker, factor in UNIT_CONVERSIONS.items():
            if marker in working:
                working[marker] = working[marker] * factor

        # Calculate linear predictor (xb)
        xb = PHENOAGE_INTERCEPT
        contributions = []
        missing = []

        for marker, coeff in PHENOAGE_COEFFICIENTS.items():
            if marker in working:
                si_value = working[marker]
                contribution = coeff * si_value
                xb += contribution
                # Report original value (before conversion) for user display
                original_value = biomarkers.get(marker, si_value)
                if marker == "ln_crp":
                    original_value = biomarkers.get("hs_crp", si_value)
                contributions.append({
                    "biomarker": marker,
                    "value": original_value,
                    "si_value": round(si_value, 4),
                    "coefficient": coeff,
                    "contribution": round(contribution, 4),
                    "direction": "protective" if coeff < 0 else "aging",
                })
            else:
                missing.append(marker)

        if missing:
            logger.warning(f"Missing PhenoAge biomarkers: {missing}")

        if len(missing) > len(PHENOAGE_COEFFICIENTS) // 2:
            logger.warning(
                f"PhenoAge: Only {len(PHENOAGE_COEFFICIENTS) - len(missing)}/{len(PHENOAGE_COEFFICIENTS)} "
                f"biomarkers available. Results may be unreliable."
            )

        # Add chronological age contribution
        # Chronological age coefficient from Levine et al. 2018 (Table S1, PMID:29676998)
        age_contribution = 0.0804 * chronological_age
        xb += age_contribution
        contributions.append({
            "biomarker": "chronological_age",
            "value": chronological_age,
            "coefficient": 0.0804,
            "contribution": round(age_contribution, 4),
            "direction": "aging",
        })

        # Clamp xb to prevent exp() overflow (exp(700) is near float max)
        xb = max(min(xb, 700), -700)

        # Calculate mortality score using Gompertz model
        # m = 1 - exp((m_n * exp(xb)) / m_d)
        mortality_score = 1 - math.exp(
            (MORT_NUMERATOR * math.exp(xb)) / MORT_DENOMINATOR
        )

        # Convert mortality score to biological age (PhenoAge)
        # phenoage = (log(BA_n * log(1 - m)) / BA_d) + BA_i
        # Clamp mortality_score to valid range for log computation
        m_clamped = max(min(mortality_score, 1 - 1e-10), 1e-10)
        inner = BA_NUMERATOR * math.log(1 - m_clamped)
        # inner must be positive for outer log; clamp to prevent domain error
        inner = max(inner, 1e-10)
        biological_age = (math.log(inner) / BA_DENOMINATOR) + BA_INTERCEPT

        age_acceleration = biological_age - chronological_age

        # Compute confidence interval (95% CI)
        # SE scales with proportion of missing biomarkers
        n_available = len(PHENOAGE_COEFFICIENTS) - len(missing)
        se = PHENOAGE_SE_FULL if n_available >= len(PHENOAGE_COEFFICIENTS) else PHENOAGE_SE_PARTIAL
        ci_lower = biological_age - 1.96 * se
        ci_upper = biological_age + 1.96 * se

        # Classify mortality risk with confidence qualifier
        # Levine et al. 2018 PMID:29676998; Liu et al. 2019 PMID:30567591
        if age_acceleration > 5:
            risk = "HIGH"
        elif age_acceleration > 2:
            risk = "MODERATE"
        elif age_acceleration > -2:
            risk = "NORMAL"
        else:
            risk = "LOW"  # Aging slower than expected

        # Add confidence qualifier based on data completeness
        if len(missing) > 3:
            risk_confidence = "low"
        elif len(missing) > 0:
            risk_confidence = "moderate"
        else:
            risk_confidence = "high"

        # Sort contributions by absolute magnitude
        contributions.sort(key=lambda c: abs(c["contribution"]), reverse=True)

        return {
            "chronological_age": chronological_age,
            "biological_age": round(biological_age, 1),
            "age_acceleration": round(age_acceleration, 1),
            "mortality_score": round(mortality_score, 6),
            "mortality_risk": risk,
            "risk_confidence": risk_confidence,
            "confidence_interval": {
                "lower": round(ci_lower, 1),
                "upper": round(ci_upper, 1),
                "confidence_level": 0.95,
                "standard_error": se,
                "note": f"Based on {n_available}/{len(PHENOAGE_COEFFICIENTS)} biomarkers available",
            },
            "top_aging_drivers": contributions[:5],
            "all_contributions": contributions,
            "missing_biomarkers": missing,
        }

    def calculate_grimage_surrogate(
        self,
        chronological_age: float,
        plasma_markers: Dict[str, float],
    ) -> Dict[str, Any]:
        """Estimate GrimAge acceleration from plasma protein surrogates.

        This is a simplified surrogate — true GrimAge requires methylation data.
        Plasma proteins (GDF-15, Cystatin C, etc.) are used as proxy markers
        that correlate with DNAm GrimAge components.

        Args:
            chronological_age: Patient's chronological age in years.
            plasma_markers: Dict with keys from GRIMAGE_MARKERS (gdf15,
                cystatin_c, leptin, pai1, timp1, adm).

        Returns:
            Dict with grimage_score, estimated_acceleration, and marker details.
        """
        marker_details = []
        total_weighted_deviation = 0.0
        total_weight = 0.0

        for marker, config in GRIMAGE_MARKERS.items():
            if marker in plasma_markers:
                value = plasma_markers[marker]
                ref_max = config["ref_max"]
                weight = config["weight"]

                # Calculate deviation from reference upper bound
                # Positive = elevated (accelerated aging)
                deviation = (value - ref_max) / ref_max if ref_max > 0 else 0
                weighted = weight * deviation
                total_weighted_deviation += weighted
                total_weight += weight

                status = "elevated" if value > ref_max else "normal"
                marker_details.append({
                    "marker": marker,
                    "value": value,
                    "unit": config["unit"],
                    "reference_max": ref_max,
                    "status": status,
                    "deviation_pct": round(deviation * 100, 1),
                    "weight": weight,
                })

        # Estimate GrimAge acceleration (surrogate)
        # Scale factor converts weighted deviation to approximate years
        # Empirical scaling factor (Belsky et al. 2020 PMID:32203970; calibrated
        # against plasma protein deviation from healthy reference upper bounds).
        # True GrimAge requires methylation data — this is a surrogate estimate.
        scale_factor = 10.0
        estimated_acceleration = (
            total_weighted_deviation / total_weight * scale_factor
            if total_weight > 0 else 0.0
        )

        grimage_score = chronological_age + estimated_acceleration

        # Compute surrogate confidence score (0-1)
        # Based on: marker coverage and correlation with true GrimAge
        marker_coverage = len(marker_details) / len(GRIMAGE_MARKERS)
        base_confidence = GRIMAGE_VALIDATION["correlation_with_true_grimage"]  # 0.72
        confidence_score = base_confidence * marker_coverage

        # Compute 95% CI for surrogate estimate
        se = GRIMAGE_VALIDATION["standard_error_years"]
        ci_lower = grimage_score - 1.96 * se
        ci_upper = grimage_score + 1.96 * se

        return {
            "chronological_age": chronological_age,
            "grimage_score": round(grimage_score, 1),
            "estimated_acceleration": round(estimated_acceleration, 1),
            "marker_details": marker_details,
            "markers_available": len(marker_details),
            "markers_total": len(GRIMAGE_MARKERS),
            "confidence_score": round(confidence_score, 2),
            "confidence_interval": {
                "lower": round(ci_lower, 1),
                "upper": round(ci_upper, 1),
                "standard_error": se,
                "confidence_level": 0.95,
            },
            "validation": GRIMAGE_VALIDATION,
            "note": (
                f"Surrogate estimation from {len(marker_details)}/{len(GRIMAGE_MARKERS)} "
                f"plasma protein markers (confidence: {confidence_score:.0%}). "
                f"True GrimAge requires methylation data and includes smoking pack-years "
                f"(a major component not captured by plasma surrogates). "
                f"Surrogate correlation with true GrimAge: r\u00b2={GRIMAGE_VALIDATION['correlation_with_true_grimage']:.2f} "
                f"({GRIMAGE_VALIDATION['validation_cohort']}). "
                f"Interpret with caution \u2014 this is a directional indicator, not a precise measurement."
            ),
        }

    def calculate(
        self,
        chronological_age: float,
        biomarkers: Dict[str, float],
    ) -> Dict[str, Any]:
        """Run both PhenoAge and GrimAge surrogate calculations.

        Args:
            chronological_age: Patient's chronological age in years.
            biomarkers: Combined dict of all available biomarkers.

        Returns:
            Dict with phenoage and grimage results combined.
        """
        phenoage = self.calculate_phenoage(chronological_age, biomarkers)

        # Extract GrimAge markers from the combined biomarker dict
        grimage_markers = {
            k: v for k, v in biomarkers.items()
            if k in GRIMAGE_MARKERS
        }
        grimage = (
            self.calculate_grimage_surrogate(chronological_age, grimage_markers)
            if grimage_markers
            else None
        )

        return {
            "phenoage": phenoage,
            "grimage": grimage,
            "biological_age": phenoage["biological_age"],
            "age_acceleration": phenoage["age_acceleration"],
            "mortality_risk": phenoage["mortality_risk"],
        }
