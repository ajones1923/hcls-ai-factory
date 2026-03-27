"""Genotype-guided dosing calculators for dose-critical drugs.

Implements validated pharmacogenomic dosing algorithms for drugs where
genotype information significantly impacts optimal dose selection.

Algorithms:
    - Warfarin: IWPC algorithm (International Warfarin Pharmacogenetics Consortium)
    - Tacrolimus: CPIC CYP3A5-guided dosing
    - Fluoropyrimidines: DPYD activity score-based dose reduction
    - Thiopurines: TPMT + NUDT15 combined dosing
    - Clopidogrel: CPIC CYP2C19-guided therapy selection
    - Simvastatin: CPIC SLCO1B1-guided dose adjustment
    - SSRIs: CPIC CYP2C19-guided SSRI (citalopram/escitalopram) dosing
    - Phenytoin: CPIC CYP2C9-guided phenytoin dosing
    - TCAs: CPIC CYP2D6/CYP2C19-guided tricyclic antidepressant dosing

References:
    - IWPC: Klein et al., N Engl J Med 2009; 360:753-764
    - CPIC warfarin guideline: Johnson et al., Clin Pharmacol Ther 2017; 102:397
    - CPIC tacrolimus guideline: Birdwell et al., Clin Pharmacol Ther 2015; 98:19
    - CPIC fluoropyrimidine guideline: Amstutz et al., Clin Pharmacol Ther 2018; 103:210
    - CPIC thiopurine guideline: Relling et al., Clin Pharmacol Ther 2019; 105:1095
    - CPIC clopidogrel guideline: Lee et al., Clin Pharmacol Ther 2022; 112:959
    - CPIC simvastatin guideline: Cooper-DeHoff et al., Clin Pharmacol Ther 2022; 112:1125
    - CPIC SSRI guideline: Hicks et al., Clin Pharmacol Ther 2015; 98:127 (2023 Update)
    - CPIC phenytoin guideline: Karnes et al., Clin Pharmacol Ther 2021; 109:302 (2020 Update)
    - CPIC TCA guideline: Hicks et al., Clin Pharmacol Ther 2017; 102:37 (2016 Update)

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


# ═════════════════════════════════════════════════════════════════════════════
# Activity score lookups (shared with PhenotypeTranslator but duplicated
# here for module independence)
# ═════════════════════════════════════════════════════════════════════════════

_DPYD_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,
    "*2A": 0.0,
    "*13": 0.0,
    "c.2846A>T": 0.5,
    # Additional DPYD alleles
    "*7": 0.5,       # decreased function
    "*9A": 0.5,      # decreased function
    "*10": 0.5,      # decreased function
    "HapB3": 0.5,    # rs75017182, decreased function
}

_TPMT_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,
    "*2": 0.0,
    "*3A": 0.0,
    "*3B": 0.0,
    "*3C": 0.0,
    "*4": 0.0,
}

_NUDT15_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,
    "*2": 0.0,
    "*3": 0.0,
    "*4": 0.0,
    "*5": 0.0,
    "*6": 1.0,  # normal function
}

_CYP3A5_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,  # expresser
    "*3": 0.0,  # non-expresser (most common)
    "*6": 0.0,  # non-expresser
    "*7": 0.0,  # non-expresser
}

_CYP2C19_ACTIVITY: Dict[str, float] = {
    "*1": 1.0, "*2": 0.0, "*3": 0.0, "*4": 0.0, "*5": 0.0,
    "*6": 0.0, "*7": 0.0, "*8": 0.0, "*9": 0.0, "*10": 0.0,
    "*17": 1.5,
}

_SLCO1B1_FUNCTION: Dict[str, Dict[str, str]] = {
    "*1/*1": {"phenotype": "Normal Function", "risk": "average"},
    "*1/*5": {"phenotype": "Decreased Function", "risk": "increased"},
    "*5/*5": {"phenotype": "Poor Function", "risk": "high"},
    "*1/*15": {"phenotype": "Decreased Function", "risk": "increased"},
    "*1/*17": {"phenotype": "Decreased Function", "risk": "increased"},
    "*15/*15": {"phenotype": "Poor Function", "risk": "high"},
    "*5/*15": {"phenotype": "Poor Function", "risk": "high"},
    "*5/*17": {"phenotype": "Poor Function", "risk": "high"},
}

_CYP2C9_ACTIVITY: Dict[str, float] = {
    "*1": 1.0, "*2": 0.5, "*3": 0.0, "*5": 0.0, "*6": 0.0,
    "*8": 0.5, "*11": 0.5,
}

_CYP2D6_ACTIVITY: Dict[str, float] = {
    "*1": 1.0, "*2": 1.0, "*3": 0.0, "*4": 0.0, "*5": 0.0,
    "*6": 0.0, "*7": 0.0, "*8": 0.0, "*9": 0.5, "*10": 0.25,
    "*11": 0.0, "*12": 0.0, "*14": 0.0, "*15": 0.0,
    "*17": 0.5, "*29": 0.5, "*35": 1.0, "*36": 0.0,
    "*40": 0.0, "*41": 0.5, "*45": 1.0,
}

_TCA_CLASS: Dict[str, str] = {
    "amitriptyline": "tertiary", "clomipramine": "tertiary",
    "doxepin": "tertiary", "imipramine": "tertiary",
    "trimipramine": "tertiary",
    "nortriptyline": "secondary", "desipramine": "secondary",
}

_TCA_TDM_RANGES: Dict[str, str] = {
    "amitriptyline": "amitriptyline+nortriptyline: 80-200 ng/mL",
    "nortriptyline": "nortriptyline: 50-150 ng/mL",
    "clomipramine": "clomipramine+desmethylclomipramine: 230-450 ng/mL",
    "doxepin": "doxepin+nordoxepin: 150-300 ng/mL",
    "desipramine": "desipramine: 100-300 ng/mL",
    "imipramine": "imipramine+desipramine: 175-300 ng/mL",
    "trimipramine": "trimipramine: 150-300 ng/mL",
}


# ═════════════════════════════════════════════════════════════════════════════
# Dosing Calculator
# ═════════════════════════════════════════════════════════════════════════════


class DosingCalculator:
    """Genotype-guided dosing calculators for dose-critical drugs.

    Each method implements a validated clinical algorithm that incorporates
    genetic and clinical variables to predict optimal drug dosing.
    """

    # ─────────────────────────────────────────────────────────────────────
    # 1. WARFARIN — IWPC pharmacogenomic dosing algorithm
    # ─────────────────────────────────────────────────────────────────────

    def warfarin_dose(
        self,
        cyp2c9_diplotype: str,
        vkorc1_genotype: str,
        cyp4f2_genotype: str = "*1/*1",
        age: int = 60,
        weight: float = 75.0,
        height: float = 170.0,
        race: str = "other",
        amiodarone: bool = False,
        smoker: bool = False,
        enzyme_inducer: bool = False,
    ) -> Dict[str, Any]:
        """Calculate predicted warfarin weekly dose using the IWPC algorithm.

        The IWPC pharmacogenetic algorithm (Klein et al., NEJM 2009) uses
        CYP2C9, VKORC1, and clinical variables to predict the stable
        therapeutic warfarin dose.

        Formula (square-root-transformed dose):
        sqrt(weekly_dose) = 5.6044
            - 0.2614 * age_decades
            + 0.0087 * height_cm
            + 0.0128 * weight_kg
            - 0.8677 * VKORC1_AG
            - 1.6974 * VKORC1_AA
            - 0.5211 * CYP2C9_1_2
            - 0.9357 * CYP2C9_1_3
            - 1.0616 * CYP2C9_2_2
            - 1.9206 * CYP2C9_2_3
            - 2.3312 * CYP2C9_3_3
            - 0.2188 * CYP4F2_1_3
            - 0.2760 * CYP4F2_3_3
            + 1.1816 * african_american
            - 0.1070 * asian
            - 0.2029 * amiodarone
            + 0.2107 * smoker
            + 1.2799 * enzyme_inducer

        Parameters
        ----------
        cyp2c9_diplotype : str
            e.g. "*1/*1", "*1/*2", "*2/*3", "*3/*3"
        vkorc1_genotype : str
            e.g. "*1/*1" (GG), "*1/-1639G>A" (AG), "-1639G>A/-1639G>A" (AA)
        cyp4f2_genotype : str
            e.g. "*1/*1", "*1/*3", "*3/*3"
        age : int
            Patient age in years.
        weight : float
            Weight in kg.
        height : float
            Height in cm.
        race : str
            One of "african_american", "asian", "european", "other".
        amiodarone : bool
            Whether the patient is taking amiodarone.
        smoker : bool
            Whether the patient is a current smoker.
        enzyme_inducer : bool
            Whether the patient is taking a CYP enzyme inducer
            (rifampin, carbamazepine, phenytoin, phenobarbital).

        Returns
        -------
        dict
            predicted_weekly_dose (mg), predicted_daily_dose (mg),
            confidence_interval (tuple), algorithm_name, variables_used (dict),
            dose_category (str), clinical_notes (list).
        """
        # ── Encode genetic variables ──
        cyp2c9_norm = self._normalize_diplotype(cyp2c9_diplotype)
        cyp2c9_1_2 = 1.0 if cyp2c9_norm == "*1/*2" else 0.0
        cyp2c9_1_3 = 1.0 if cyp2c9_norm == "*1/*3" else 0.0
        cyp2c9_2_2 = 1.0 if cyp2c9_norm == "*2/*2" else 0.0
        cyp2c9_2_3 = 1.0 if cyp2c9_norm == "*2/*3" else 0.0
        cyp2c9_3_3 = 1.0 if cyp2c9_norm == "*3/*3" else 0.0

        # VKORC1 encoding
        vkorc1_ag = 1.0 if self._is_vkorc1_ag(vkorc1_genotype) else 0.0
        vkorc1_aa = 1.0 if self._is_vkorc1_aa(vkorc1_genotype) else 0.0

        # CYP4F2 encoding
        cyp4f2_norm = self._normalize_diplotype(cyp4f2_genotype)
        cyp4f2_1_3 = 1.0 if cyp4f2_norm == "*1/*3" else 0.0
        cyp4f2_3_3 = 1.0 if cyp4f2_norm == "*3/*3" else 0.0

        # Race encoding
        race_lower = race.lower().strip()
        african_american = 1.0 if race_lower in ("african_american", "african american", "black") else 0.0
        asian = 1.0 if race_lower in ("asian", "east_asian", "southeast_asian") else 0.0

        # Clinical variables
        amiodarone_val = 1.0 if amiodarone else 0.0
        smoker_val = 1.0 if smoker else 0.0
        enzyme_inducer_val = 1.0 if enzyme_inducer else 0.0
        age_decades = age / 10.0

        # ── IWPC algorithm ──
        sqrt_dose = (
            5.6044
            - 0.2614 * age_decades
            + 0.0087 * height
            + 0.0128 * weight
            - 0.8677 * vkorc1_ag
            - 1.6974 * vkorc1_aa
            - 0.5211 * cyp2c9_1_2
            - 0.9357 * cyp2c9_1_3
            - 1.0616 * cyp2c9_2_2
            - 1.9206 * cyp2c9_2_3
            - 2.3312 * cyp2c9_3_3
            - 0.2188 * cyp4f2_1_3
            - 0.2760 * cyp4f2_3_3
            + 1.1816 * african_american
            - 0.1070 * asian
            - 0.2029 * amiodarone_val
            + 0.2107 * smoker_val
            + 1.2799 * enzyme_inducer_val
        )

        # Square the result to get weekly dose
        weekly_dose = max(sqrt_dose ** 2, 0.5)  # floor at 0.5 mg/week
        daily_dose = weekly_dose / 7.0

        # Approximate confidence interval (IWPC reports RMSE ~8.5 mg/week
        # for the pharmacogenetic algorithm)
        rmse = 8.5
        ci_low = max(weekly_dose - 1.96 * rmse, 0.5)
        ci_high = weekly_dose + 1.96 * rmse

        # Dose category
        dose_category = self._warfarin_dose_category(weekly_dose)

        # Clinical notes
        clinical_notes = self._warfarin_clinical_notes(
            cyp2c9_norm, vkorc1_genotype, weekly_dose, amiodarone, enzyme_inducer,
        )

        result = {
            "predicted_weekly_dose": round(weekly_dose, 1),
            "predicted_daily_dose": round(daily_dose, 1),
            "confidence_interval": (round(ci_low, 1), round(ci_high, 1)),
            "algorithm_name": "IWPC Pharmacogenetic Warfarin Dosing Algorithm",
            "dose_category": dose_category,
            "variables_used": {
                "cyp2c9_diplotype": cyp2c9_diplotype,
                "vkorc1_genotype": vkorc1_genotype,
                "cyp4f2_genotype": cyp4f2_genotype,
                "age": age,
                "weight_kg": weight,
                "height_cm": height,
                "race": race,
                "amiodarone": amiodarone,
                "smoker": smoker,
                "enzyme_inducer": enzyme_inducer,
            },
            "clinical_notes": clinical_notes,
        }

        logger.info(
            "Warfarin dose prediction: %.1f mg/week (CYP2C9 %s, VKORC1 %s)",
            weekly_dose, cyp2c9_diplotype, vkorc1_genotype,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 2. TACROLIMUS — CYP3A5-guided dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def tacrolimus_dose(
        self,
        cyp3a5_diplotype: str,
        weight: float = 70.0,
    ) -> Dict[str, Any]:
        """Calculate CYP3A5 genotype-guided tacrolimus starting dose.

        CPIC guideline (Birdwell et al. 2015):
        - CYP3A5 expressers (*1/*1, *1/*3): increase starting dose 1.5-2x
        - CYP3A5 non-expressers (*3/*3): standard dose (0.15-0.2 mg/kg/day)

        Parameters
        ----------
        cyp3a5_diplotype : str
            e.g. "*1/*1", "*1/*3", "*3/*3"
        weight : float
            Body weight in kg.

        Returns
        -------
        dict
            dose_recommendation, mg_per_kg_per_day, total_daily_mg,
            phenotype, rationale.
        """
        allele1, allele2 = self._parse_diplotype(cyp3a5_diplotype)
        score1 = _CYP3A5_ACTIVITY.get(allele1, 0.0)
        score2 = _CYP3A5_ACTIVITY.get(allele2, 0.0)
        total_score = score1 + score2

        # Determine phenotype and dosing
        if total_score >= 1.5:
            # Expresser (*1/*1)
            phenotype = "CYP3A5 Expresser"
            mg_per_kg = 0.30
            recommendation = (
                "CYP3A5 Expresser: increase starting dose to 0.3 mg/kg/day "
                "(1.5-2x standard). Rapid CYP3A5-mediated metabolism will "
                "result in lower trough levels at standard doses. Titrate "
                "based on trough level monitoring (target 10-15 ng/mL early "
                "post-transplant)."
            )
            rationale = (
                "Functional CYP3A5 protein significantly increases tacrolimus "
                "clearance. Multiple studies show expressers require 30-50% "
                "higher doses to achieve target trough concentrations."
            )
        elif total_score >= 0.5:
            # Intermediate expresser (*1/*3, *1/*6)
            phenotype = "CYP3A5 Intermediate Expresser"
            mg_per_kg = 0.25
            recommendation = (
                "CYP3A5 Intermediate Expresser: increase starting dose to "
                "0.25 mg/kg/day (~1.25-1.5x standard). One functional "
                "CYP3A5 allele confers intermediate clearance. Titrate "
                "based on trough level monitoring."
            )
            rationale = (
                "Heterozygous CYP3A5 function; intermediate metabolism. "
                "CPIC recommends increased starting dose with close "
                "therapeutic drug monitoring."
            )
        else:
            # Non-expresser (*3/*3, *3/*6, *6/*6)
            phenotype = "CYP3A5 Non-expresser"
            mg_per_kg = 0.15
            recommendation = (
                "CYP3A5 Non-expresser: standard dose of 0.15 mg/kg/day. "
                "No functional CYP3A5 — metabolism primarily via CYP3A4. "
                "Standard trough monitoring applies."
            )
            rationale = (
                "No functional CYP3A5 alleles detected. ~80% of European "
                "and ~30% of African patients are CYP3A5 non-expressers. "
                "Standard tacrolimus dosing is appropriate."
            )

        total_daily = round(mg_per_kg * weight, 1)
        # Tacrolimus is typically dosed BID
        per_dose = round(total_daily / 2, 1)

        result = {
            "dose_recommendation": recommendation,
            "mg_per_kg_per_day": mg_per_kg,
            "total_daily_mg": total_daily,
            "per_dose_mg": per_dose,
            "frequency": "BID (every 12 hours)",
            "phenotype": phenotype,
            "cyp3a5_diplotype": cyp3a5_diplotype,
            "activity_score": total_score,
            "rationale": rationale,
            "monitoring": (
                "Therapeutic drug monitoring (trough levels) required. "
                "Target: 10-15 ng/mL early post-transplant, 5-10 ng/mL "
                "maintenance (varies by organ and protocol)."
            ),
        }

        logger.info(
            "Tacrolimus dose: %.1f mg/day (CYP3A5 %s → %s)",
            total_daily, cyp3a5_diplotype, phenotype,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 3. FLUOROPYRIMIDINES — DPYD activity score-based dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def fluoropyrimidine_dose(
        self,
        dpyd_diplotype: str,
    ) -> Dict[str, Any]:
        """Calculate DPYD genotype-guided fluoropyrimidine dose adjustment.

        CPIC guideline (Amstutz et al. 2018):
        - Activity score 2.0: normal, full dose
        - Activity score 1.0-1.5: intermediate, reduce by 50%
        - Activity score 0.5: intermediate, reduce by 50% (or 25-50%)
        - Activity score 0.0: deficient, AVOID (or reduce to 25% with monitoring)

        Parameters
        ----------
        dpyd_diplotype : str
            e.g. "*1/*1", "*1/*2A", "*2A/*2A", "*1/c.2846A>T"

        Returns
        -------
        dict
            activity_score, dose_reduction_percent, recommendation,
            phenotype, affected_drugs (list).
        """
        allele1, allele2 = self._parse_diplotype(dpyd_diplotype)
        score1 = _DPYD_ACTIVITY.get(allele1, 1.0)
        score2 = _DPYD_ACTIVITY.get(allele2, 1.0)
        activity_score = score1 + score2

        affected_drugs = [
            "5-fluorouracil (5-FU)",
            "capecitabine",
            "tegafur",
        ]

        if activity_score == 0.0:
            phenotype = "DPD Deficient"
            dose_reduction = 100
            recommendation = (
                "AVOID all fluoropyrimidines. Complete DPD deficiency — "
                "life-threatening and potentially fatal toxicity (severe "
                "mucositis, myelosuppression, neurotoxicity, hand-foot "
                "syndrome). Select an alternative non-fluoropyrimidine "
                "chemotherapy regimen. If no alternatives exist and "
                "fluoropyrimidine is essential, reduce dose to at most "
                "25% with intensive inpatient monitoring (CPIC)."
            )
        elif activity_score <= 0.5:
            phenotype = "DPD Intermediate Activity (Partial Deficiency)"
            dose_reduction = 50
            recommendation = (
                "Reduce fluoropyrimidine dose by at least 50%. Activity "
                "score 0.5 indicates significant DPD impairment. Titrate "
                "dose based on toxicity and therapeutic drug monitoring "
                "(uracil/dihydrouracil ratio) if available. Monitor "
                "closely for mucositis, cytopenias, and neurotoxicity."
            )
        elif activity_score <= 1.0:
            phenotype = "DPD Intermediate Activity"
            dose_reduction = 50
            recommendation = (
                "Reduce fluoropyrimidine dose by 50%. CPIC recommendation "
                "for activity score 1.0. May titrate up based on toxicity "
                "in subsequent cycles. Monitor CBC, mucositis."
            )
        elif activity_score <= 1.5:
            phenotype = "DPD Intermediate Activity (Mild)"
            dose_reduction = 25
            recommendation = (
                "Reduce fluoropyrimidine dose by 25-50%. Activity score "
                "1.5 indicates one decreased-function allele. Titrate "
                "based on tolerability. Enhanced toxicity monitoring "
                "recommended for first 2 cycles."
            )
        else:
            phenotype = "DPD Normal Activity"
            dose_reduction = 0
            recommendation = (
                "Standard fluoropyrimidine dose. No DPYD-based dose "
                "adjustment needed. Activity score 2.0 (normal). Standard "
                "toxicity monitoring applies."
            )

        result = {
            "activity_score": activity_score,
            "dose_reduction_percent": dose_reduction,
            "recommendation": recommendation,
            "phenotype": phenotype,
            "dpyd_diplotype": dpyd_diplotype,
            "allele_scores": {allele1: score1, allele2: score2},
            "affected_drugs": affected_drugs,
            "monitoring": (
                "CBC with differential prior to each cycle. Monitor for "
                "mucositis, diarrhea, hand-foot syndrome, and neurotoxicity. "
                "Consider uracil test dose or DPD phenotyping for borderline "
                "cases."
            ),
        }

        logger.info(
            "DPYD dosing: activity score %.1f → %d%% dose reduction (%s)",
            activity_score, dose_reduction, dpyd_diplotype,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 4. THIOPURINES — TPMT + NUDT15 combined dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def thiopurine_dose(
        self,
        tpmt_diplotype: str,
        nudt15_diplotype: str = "*1/*1",
    ) -> Dict[str, Any]:
        """Calculate TPMT + NUDT15 genotype-guided thiopurine dose.

        CPIC guideline (Relling et al. 2019) considers both TPMT and NUDT15
        activity scores.  The more restrictive (lower) activity score
        determines the dose recommendation.

        Parameters
        ----------
        tpmt_diplotype : str
            e.g. "*1/*1", "*1/*3A", "*3A/*3C"
        nudt15_diplotype : str
            e.g. "*1/*1", "*1/*3", "*3/*3"

        Returns
        -------
        dict
            dose_recommendation, tpmt_phenotype, nudt15_phenotype,
            combined_recommendation, affected_drugs.
        """
        # TPMT activity
        t_a1, t_a2 = self._parse_diplotype(tpmt_diplotype)
        tpmt_score = _TPMT_ACTIVITY.get(t_a1, 1.0) + _TPMT_ACTIVITY.get(t_a2, 1.0)

        # NUDT15 activity
        n_a1, n_a2 = self._parse_diplotype(nudt15_diplotype)
        nudt15_score = _NUDT15_ACTIVITY.get(n_a1, 1.0) + _NUDT15_ACTIVITY.get(n_a2, 1.0)

        # Determine individual phenotypes
        tpmt_phenotype = self._thiopurine_phenotype(tpmt_score, "TPMT")
        nudt15_phenotype = self._thiopurine_phenotype(nudt15_score, "NUDT15")

        # Combined recommendation: use the more restrictive phenotype
        combined_score = min(tpmt_score, nudt15_score)

        affected_drugs = [
            "azathioprine",
            "6-mercaptopurine (6-MP)",
            "6-thioguanine (6-TG)",
        ]

        if combined_score == 0.0:
            dose_recommendation = (
                "AVOID thiopurines if possible. If essential, reduce dose "
                "to 10% of standard (e.g., azathioprine 0.5 mg/kg/day "
                "thrice weekly) with intensive monitoring. Both TPMT and/or "
                "NUDT15 deficient — extremely high risk of life-threatening "
                "myelosuppression (pancytopenia). Consider alternative "
                "immunosuppressant (mycophenolate mofetil)."
            )
            dose_percent = 10
        elif combined_score <= 1.0:
            dose_recommendation = (
                "Reduce thiopurine dose to 30-70% of standard. Start at "
                "lower end (30-50%) and titrate based on CBC. Heterozygous "
                "for TPMT and/or NUDT15 variant — intermediate risk of "
                "myelosuppression. Monitor CBC weekly for first 8 weeks, "
                "then every 2-4 weeks."
            )
            dose_percent = 50
        else:
            dose_recommendation = (
                "Standard thiopurine dose. Normal TPMT and NUDT15 activity. "
                "Routine CBC monitoring applies (at least every 3 months). "
                "Note: even with normal genotype, idiosyncratic toxicity "
                "is possible — maintain standard monitoring."
            )
            dose_percent = 100

        result = {
            "dose_recommendation": dose_recommendation,
            "dose_percent_of_standard": dose_percent,
            "tpmt_diplotype": tpmt_diplotype,
            "tpmt_activity_score": tpmt_score,
            "tpmt_phenotype": tpmt_phenotype,
            "nudt15_diplotype": nudt15_diplotype,
            "nudt15_activity_score": nudt15_score,
            "nudt15_phenotype": nudt15_phenotype,
            "combined_activity_score": combined_score,
            "affected_drugs": affected_drugs,
            "monitoring": (
                "CBC with differential: weekly x8 weeks, then biweekly x4 "
                "weeks, then at least every 3 months. More frequent "
                "monitoring for intermediate/deficient patients. Monitor "
                "thiopurine metabolites (6-TGN, 6-MMP) to guide dosing."
            ),
        }

        logger.info(
            "Thiopurine dosing: TPMT %s (AS %.1f), NUDT15 %s (AS %.1f) → %d%% dose",
            tpmt_diplotype, tpmt_score, nudt15_diplotype, nudt15_score, dose_percent,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 5. CLOPIDOGREL — CYP2C19-guided therapy selection (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def clopidogrel_dose(
        self,
        cyp2c19_diplotype: str,
    ) -> Dict[str, Any]:
        """Calculate CYP2C19 genotype-guided clopidogrel therapy selection.

        CPIC guideline (Lee et al. 2022):
        - Poor metabolizers: alternative antiplatelet therapy required
        - Intermediate metabolizers: consider alternative therapy
        - Normal/Rapid/Ultrarapid metabolizers: standard clopidogrel

        Parameters
        ----------
        cyp2c19_diplotype : str
            e.g. "*1/*1", "*1/*2", "*2/*2", "*1/*17", "*17/*17"

        Returns
        -------
        dict
            cyp2c19_diplotype, activity_score, phenotype,
            recommended_therapy, clopidogrel_use, rationale,
            warnings, affected_drugs, reference.
        """
        allele1, allele2 = self._parse_diplotype(cyp2c19_diplotype)
        score1 = _CYP2C19_ACTIVITY.get(allele1, 1.0)
        score2 = _CYP2C19_ACTIVITY.get(allele2, 1.0)
        activity_score = score1 + score2

        # Determine phenotype using half-open intervals
        if activity_score == 0.0:
            phenotype = "Poor Metabolizer"
            recommended_therapy = "prasugrel 10mg daily OR ticagrelor 90mg BID"
            clopidogrel_use = "contraindicated"
            rationale = (
                "Alternative antiplatelet therapy REQUIRED. Clopidogrel is "
                "INEFFECTIVE in CYP2C19 poor metabolizers — no conversion "
                "to active metabolite. Use prasugrel 10mg daily or ticagrelor "
                "90mg BID. CPIC Level A."
            )
            warnings = [
                "Clopidogrel provides no antiplatelet benefit in PMs",
                "Elevated risk of MACE, stent thrombosis, and recurrent ACS",
            ]
        elif activity_score < 2.0:
            phenotype = "Intermediate Metabolizer"
            recommended_therapy = (
                "prasugrel 10mg daily OR ticagrelor 90mg BID (preferred); "
                "clopidogrel 75mg with enhanced monitoring "
                "(if alternatives unavailable)"
            )
            clopidogrel_use = "use_with_caution"
            rationale = (
                "Consider alternative antiplatelet therapy. CYP2C19 "
                "intermediate metabolizers have significantly reduced "
                "clopidogrel activation. Prasugrel or ticagrelor preferred "
                "per CPIC guidelines."
            )
            warnings = [
                "30-40% reduction in active metabolite formation",
                "Increased residual platelet reactivity",
            ]
        elif activity_score < 2.5:
            phenotype = "Normal Metabolizer"
            recommended_therapy = "clopidogrel 75mg daily (standard)"
            clopidogrel_use = "standard"
            rationale = (
                "Standard clopidogrel therapy appropriate. Normal CYP2C19 "
                "metabolism expected."
            )
            warnings = []
        elif activity_score < 3.0:
            phenotype = "Rapid Metabolizer"
            recommended_therapy = "clopidogrel 75mg daily (standard)"
            clopidogrel_use = "standard"
            rationale = (
                "Standard clopidogrel therapy appropriate. Enhanced CYP2C19 "
                "metabolism may improve antiplatelet response."
            )
            warnings = [
                "Monitor for increased bleeding risk",
            ]
        else:
            phenotype = "Ultrarapid Metabolizer"
            recommended_therapy = "clopidogrel 75mg daily (standard)"
            clopidogrel_use = "standard"
            rationale = (
                "Standard clopidogrel therapy appropriate. Ultra-rapid "
                "metabolism may increase active metabolite and bleeding risk."
            )
            warnings = [
                "Possible increased bleeding risk — monitor clinically",
                "Consider lower dose if bleeding complications occur",
            ]

        result = {
            "cyp2c19_diplotype": cyp2c19_diplotype,
            "activity_score": activity_score,
            "phenotype": phenotype,
            "recommended_therapy": recommended_therapy,
            "clopidogrel_use": clopidogrel_use,
            "rationale": rationale,
            "warnings": warnings,
            "affected_drugs": ["clopidogrel"],
            "reference": (
                "CPIC Guideline for Clopidogrel and CYP2C19, 2022 Update"
            ),
        }

        logger.info(
            "Clopidogrel therapy: %s (CYP2C19 %s, AS %.1f)",
            phenotype, cyp2c19_diplotype, activity_score,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 6. SIMVASTATIN — SLCO1B1-guided dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def simvastatin_dose(
        self,
        slco1b1_diplotype: str,
    ) -> Dict[str, Any]:
        """Calculate SLCO1B1 genotype-guided simvastatin dose adjustment.

        CPIC guideline (Cooper-DeHoff et al. 2022):
        - Normal function: standard simvastatin dosing
        - Decreased function: limit to 20mg/day
        - Poor function: avoid simvastatin

        Parameters
        ----------
        slco1b1_diplotype : str
            e.g. "*1/*1", "*1/*5", "*5/*5", "*1/*15"

        Returns
        -------
        dict
            slco1b1_diplotype, phenotype, max_simvastatin_dose_mg,
            dose_adjustment, myopathy_risk, rationale, warnings,
            alternatives, affected_drugs, reference.
        """
        normalized = self._normalize_diplotype(slco1b1_diplotype)
        lookup = _SLCO1B1_FUNCTION.get(normalized)

        if lookup is None:
            phenotype = "Unknown"
            risk_level = "unknown"
        else:
            phenotype = lookup["phenotype"]
            risk_level = lookup["risk"]

        if phenotype == "Normal Function":
            max_dose = 80
            dose_adjustment = "none"
            myopathy_risk = "average (population baseline)"
            rationale = (
                "Standard simvastatin dosing. No increased myopathy risk "
                "from SLCO1B1 genotype."
            )
            warnings: List[str] = []
        elif phenotype == "Decreased Function":
            max_dose = 20
            dose_adjustment = "reduce_to_20mg_max"
            myopathy_risk = "increased (4.5-fold with 80mg dose)"
            rationale = (
                "Prescribe simvastatin ≤20mg/day. SLCO1B1 decreased "
                "function increases simvastatin exposure and myopathy "
                "risk. Consider rosuvastatin or pravastatin as safer "
                "alternatives."
            )
            warnings = [
                "Do NOT use simvastatin 80mg",
                "Monitor for muscle pain, weakness, elevated CK",
            ]
        elif phenotype == "Poor Function":
            max_dose = 0
            dose_adjustment = "avoid_simvastatin"
            myopathy_risk = "high (17-fold increased with standard dose)"
            rationale = (
                "AVOID simvastatin. SLCO1B1 poor function confers very "
                "high myopathy/rhabdomyolysis risk. Use rosuvastatin, "
                "pravastatin, or fluvastatin instead."
            )
            warnings = [
                "Simvastatin contraindicated",
                "Also applies to lovastatin (similar transport)",
                "Rosuvastatin, pravastatin, fluvastatin are safer alternatives",
            ]
        else:
            # Unknown genotype fallback
            max_dose = 40
            dose_adjustment = "standard_precautions"
            myopathy_risk = "unknown"
            rationale = (
                "SLCO1B1 genotype not recognized. Use standard prescribing "
                "guidelines with awareness of myopathy risk."
            )
            warnings = []

        result = {
            "slco1b1_diplotype": slco1b1_diplotype,
            "phenotype": phenotype,
            "max_simvastatin_dose_mg": max_dose,
            "dose_adjustment": dose_adjustment,
            "myopathy_risk": myopathy_risk,
            "rationale": rationale,
            "warnings": warnings,
            "alternatives": [
                "rosuvastatin", "pravastatin", "fluvastatin", "pitavastatin",
            ],
            "affected_drugs": ["simvastatin", "lovastatin"],
            "reference": (
                "CPIC Guideline for Simvastatin and SLCO1B1, 2022 Update"
            ),
        }

        logger.info(
            "Simvastatin dosing: %s (SLCO1B1 %s, max %dmg)",
            phenotype, slco1b1_diplotype, max_dose,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 7. SSRIs — CYP2C19-guided SSRI dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def ssri_dose(
        self,
        drug: str,
        cyp2c19_diplotype: str,
    ) -> Dict[str, Any]:
        """Calculate CYP2C19 genotype-guided SSRI dose adjustment.

        CPIC guideline for citalopram and escitalopram (Hicks et al.,
        2023 Update).

        Parameters
        ----------
        drug : str
            One of "citalopram" or "escitalopram".
        cyp2c19_diplotype : str
            e.g. "*1/*1", "*1/*2", "*2/*2", "*1/*17", "*17/*17"

        Returns
        -------
        dict
            drug, cyp2c19_diplotype, activity_score, phenotype,
            max_dose_mg, dose_adjustment, rationale, warnings,
            alternatives, reference.
        """
        drug_lower = drug.lower().strip()

        # Validate drug
        if drug_lower not in ("citalopram", "escitalopram"):
            logger.warning(
                "SSRI dosing: unsupported drug '%s'. Only citalopram and "
                "escitalopram have CPIC CYP2C19 guidance.", drug,
            )
            return {
                "drug": drug,
                "cyp2c19_diplotype": cyp2c19_diplotype,
                "activity_score": None,
                "phenotype": "Unknown",
                "max_dose_mg": None,
                "dose_adjustment": "no_cpic_guidance",
                "rationale": (
                    f"No CPIC CYP2C19 dosing guidance available for {drug}. "
                    "CPIC SSRI guidance covers citalopram and escitalopram only."
                ),
                "warnings": [
                    f"No genotype-guided dosing recommendation for {drug}",
                ],
                "alternatives": ["sertraline", "mirtazapine", "bupropion"],
                "reference": "CPIC Guideline for SSRIs and CYP2C19, 2023 Update",
            }

        # Calculate activity score
        allele1, allele2 = self._parse_diplotype(cyp2c19_diplotype)
        score1 = _CYP2C19_ACTIVITY.get(allele1, 1.0)
        score2 = _CYP2C19_ACTIVITY.get(allele2, 1.0)
        activity_score = score1 + score2

        # Determine phenotype (same thresholds as clopidogrel_dose)
        if activity_score == 0.0:
            phenotype = "Poor Metabolizer"
        elif activity_score < 2.0:
            phenotype = "Intermediate Metabolizer"
        elif activity_score < 2.5:
            phenotype = "Normal Metabolizer"
        elif activity_score < 3.0:
            phenotype = "Rapid Metabolizer"
        else:
            phenotype = "Ultrarapid Metabolizer"

        # Drug-specific dose recommendations
        if drug_lower == "citalopram":
            if phenotype == "Poor Metabolizer":
                max_dose_mg = 20
                dose_adjustment = "50% reduction, max 20mg/day"
                rationale = (
                    "CYP2C19 PM causes ~2-fold increase in citalopram "
                    "exposure. QT prolongation risk increases with higher "
                    "plasma levels. FDA max 20mg/day for PMs."
                )
                warnings: List[str] = [
                    "Maximum 20 mg/day due to QT prolongation risk",
                    "ECG monitoring recommended",
                    "Avoid concomitant QT-prolonging drugs",
                ]
            elif phenotype == "Intermediate Metabolizer":
                max_dose_mg = 40
                dose_adjustment = "consider_dose_reduction"
                rationale = (
                    "CYP2C19 IM may have moderately increased citalopram "
                    "exposure. Standard initial dosing acceptable with "
                    "monitoring."
                )
                warnings = []
            elif phenotype == "Normal Metabolizer":
                max_dose_mg = 40
                dose_adjustment = "none"
                rationale = (
                    "Standard citalopram dosing. Normal CYP2C19 metabolism."
                )
                warnings = []
            elif phenotype == "Rapid Metabolizer":
                max_dose_mg = 40
                dose_adjustment = "none"
                rationale = (
                    "Standard dosing; slightly enhanced metabolism unlikely "
                    "to be clinically significant."
                )
                warnings = []
            else:  # Ultrarapid Metabolizer
                max_dose_mg = 60
                dose_adjustment = "consider_dose_increase"
                rationale = (
                    "CYP2C19 UM may have reduced citalopram efficacy due "
                    "to rapid metabolism. Consider higher dose or alternative "
                    "SSRI (fluoxetine, paroxetine)."
                )
                warnings = [
                    "Monitor for therapeutic response",
                    "Consider alternative not dependent on CYP2C19",
                ]
        else:  # escitalopram
            if phenotype == "Poor Metabolizer":
                max_dose_mg = 10
                dose_adjustment = "50% reduction, max 10mg/day"
                rationale = (
                    "CYP2C19 PM causes ~2-fold increase in escitalopram "
                    "exposure. QT prolongation risk increases with higher "
                    "plasma levels. Max 10mg/day for PMs."
                )
                warnings = [
                    "Maximum 10 mg/day due to QT prolongation risk",
                    "ECG monitoring recommended",
                    "Avoid concomitant QT-prolonging drugs",
                ]
            elif phenotype == "Intermediate Metabolizer":
                max_dose_mg = 20
                dose_adjustment = "consider_dose_reduction"
                rationale = (
                    "CYP2C19 IM may have moderately increased escitalopram "
                    "exposure. Standard initial dosing acceptable with "
                    "monitoring."
                )
                warnings = []
            elif phenotype == "Normal Metabolizer":
                max_dose_mg = 20
                dose_adjustment = "none"
                rationale = (
                    "Standard escitalopram dosing. Normal CYP2C19 metabolism."
                )
                warnings = []
            elif phenotype == "Rapid Metabolizer":
                max_dose_mg = 20
                dose_adjustment = "none"
                rationale = (
                    "Standard dosing; slightly enhanced metabolism unlikely "
                    "to be clinically significant."
                )
                warnings = []
            else:  # Ultrarapid Metabolizer
                max_dose_mg = 30
                dose_adjustment = "consider_dose_increase"
                rationale = (
                    "CYP2C19 UM may have reduced escitalopram efficacy due "
                    "to rapid metabolism. Consider higher dose or alternative "
                    "SSRI (fluoxetine, paroxetine)."
                )
                warnings = [
                    "Monitor for therapeutic response",
                    "Consider alternative not dependent on CYP2C19",
                ]

        result = {
            "drug": drug,
            "cyp2c19_diplotype": cyp2c19_diplotype,
            "activity_score": activity_score,
            "phenotype": phenotype,
            "max_dose_mg": max_dose_mg,
            "dose_adjustment": dose_adjustment,
            "rationale": rationale,
            "warnings": warnings,
            "alternatives": ["sertraline", "mirtazapine", "bupropion"],
            "reference": "CPIC Guideline for SSRIs and CYP2C19, 2023 Update",
        }

        logger.info(
            "SSRI dosing: %s %s (CYP2C19 %s, AS %.1f, max %dmg)",
            drug, phenotype, cyp2c19_diplotype, activity_score, max_dose_mg,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 8. PHENYTOIN — CYP2C9-guided dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def phenytoin_dose(
        self,
        cyp2c9_diplotype: str,
        weight_kg: float = 70.0,
    ) -> Dict[str, Any]:
        """Calculate CYP2C9 genotype-guided phenytoin dose adjustment.

        CPIC guideline (Karnes et al., 2020 Update):
        - PM (AS < 1.0): reduce maintenance by 50%
        - IM (1.0 <= AS < 2.0): reduce maintenance by 25%
        - NM (AS >= 2.0): standard dosing

        Parameters
        ----------
        cyp2c9_diplotype : str
            e.g. "*1/*1", "*1/*2", "*1/*3", "*2/*3", "*3/*3"
        weight_kg : float
            Body weight in kg.

        Returns
        -------
        dict
            cyp2c9_diplotype, activity_score, phenotype,
            dose_reduction_percent, standard_loading_dose_mg,
            suggested_maintenance_mg_per_day, rationale, warnings,
            tdm_target, alternatives, affected_drugs, reference.
        """
        allele1, allele2 = self._parse_diplotype(cyp2c9_diplotype)
        score1 = _CYP2C9_ACTIVITY.get(allele1, 1.0)
        score2 = _CYP2C9_ACTIVITY.get(allele2, 1.0)
        activity_score = score1 + score2

        # Determine phenotype using half-open intervals
        if activity_score < 1.0:
            phenotype = "Poor Metabolizer"
            dose_reduction_percent = 50
            rationale = (
                "CYP2C9 PM causes markedly reduced phenytoin clearance "
                "(50-75% reduction). Administer standard loading dose but "
                "reduce maintenance by 50%. Titrate using therapeutic drug "
                "monitoring (TDM). Target trough: 10-20 mcg/mL."
            )
            warnings = [
                "Reduce maintenance dose by 50%",
                "Mandatory TDM — measure levels at steady state (7-10 days)",
                "Risk of toxicity: nystagmus, ataxia, diplopia, seizure exacerbation",
                "Consider levetiracetam as alternative (no CYP2C9 dependence)",
            ]
        elif activity_score < 2.0:
            phenotype = "Intermediate Metabolizer"
            dose_reduction_percent = 25
            rationale = (
                "CYP2C9 IM has moderately reduced phenytoin clearance. "
                "Reduce maintenance by 25% and use TDM."
            )
            warnings = [
                "Reduce maintenance dose by 25%",
                "TDM recommended",
                "Monitor for dose-related side effects",
            ]
        else:
            phenotype = "Normal Metabolizer"
            dose_reduction_percent = 0
            rationale = (
                "Standard phenytoin dosing. Normal CYP2C9 metabolism."
            )
            warnings = []

        # Calculate doses
        standard_loading_dose_mg = min(round(weight_kg * 18), 1500)
        suggested_maintenance_mg = round(
            weight_kg * 5 * (1 - dose_reduction_percent / 100)
        )

        # Maintenance dose range string
        if phenotype == "Poor Metabolizer":
            max_maintenance_range = "150-200"
        elif phenotype == "Intermediate Metabolizer":
            max_maintenance_range = "225-300"
        else:
            max_maintenance_range = "300-400"

        result = {
            "cyp2c9_diplotype": cyp2c9_diplotype,
            "activity_score": activity_score,
            "phenotype": phenotype,
            "dose_reduction_percent": dose_reduction_percent,
            "standard_loading_dose_mg": standard_loading_dose_mg,
            "suggested_maintenance_mg_per_day": suggested_maintenance_mg,
            "max_maintenance_mg_per_day": max_maintenance_range,
            "rationale": rationale,
            "warnings": warnings,
            "tdm_target": "10-20 mcg/mL",
            "alternatives": ["levetiracetam", "lacosamide", "valproic acid"],
            "affected_drugs": ["phenytoin", "fosphenytoin"],
            "reference": (
                "CPIC Guideline for Phenytoin and CYP2C9, 2020 Update"
            ),
        }

        logger.info(
            "Phenytoin dosing: %s (CYP2C9 %s, AS %.1f, reduction %d%%)",
            phenotype, cyp2c9_diplotype, activity_score,
            dose_reduction_percent,
        )
        return result

    # ─────────────────────────────────────────────────────────────────────
    # 9. TCAs — CYP2D6/CYP2C19-guided dosing (CPIC)
    # ─────────────────────────────────────────────────────────────────────

    def tca_dose(
        self,
        drug: str,
        cyp2d6_diplotype: str,
        cyp2c19_diplotype: str = "*1/*1",
    ) -> Dict[str, Any]:
        """Calculate CYP2D6/CYP2C19 genotype-guided TCA dose adjustment.

        CPIC guideline (Hicks et al., 2016 Update) for tricyclic
        antidepressants.  CYP2D6 is the primary metabolizer; CYP2C19
        provides secondary adjustment for tertiary amine TCAs.

        Parameters
        ----------
        drug : str
            One of: amitriptyline, nortriptyline, clomipramine, doxepin,
            desipramine, imipramine, trimipramine.
        cyp2d6_diplotype : str
            e.g. "*1/*1", "*1/*4", "*4/*4", "*1/*41"
        cyp2c19_diplotype : str
            e.g. "*1/*1", "*1/*2", "*2/*2", "*1/*17"

        Returns
        -------
        dict
            drug, tca_class, cyp2d6_diplotype, cyp2d6_activity_score,
            cyp2d6_phenotype, cyp2c19_diplotype, cyp2c19_phenotype,
            dose_adjustment, rationale, warnings, tdm_required,
            tdm_target_range, alternatives, affected_drugs, reference.
        """
        drug_lower = drug.lower().strip()
        tca_class = _TCA_CLASS.get(drug_lower)

        if tca_class is None:
            logger.warning(
                "TCA dosing: unsupported drug '%s'. Supported: %s",
                drug, ", ".join(sorted(_TCA_CLASS.keys())),
            )
            return {
                "drug": drug,
                "tca_class": None,
                "cyp2d6_diplotype": cyp2d6_diplotype,
                "cyp2d6_activity_score": None,
                "cyp2d6_phenotype": "Unknown",
                "cyp2c19_diplotype": cyp2c19_diplotype,
                "cyp2c19_phenotype": "Unknown",
                "dose_adjustment": "no_cpic_guidance",
                "rationale": (
                    f"No CPIC TCA dosing guidance available for {drug}."
                ),
                "warnings": [f"No genotype-guided dosing for {drug}"],
                "tdm_required": False,
                "tdm_target_range": None,
                "alternatives": ["sertraline", "venlafaxine", "mirtazapine"],
                "affected_drugs": list(_TCA_CLASS.keys()),
                "reference": (
                    "CPIC Guideline for Tricyclic Antidepressants and "
                    "CYP2D6/CYP2C19, 2016 Update"
                ),
            }

        # ── CYP2D6 activity score ──
        d6_a1, d6_a2 = self._parse_diplotype(cyp2d6_diplotype)
        d6_score1 = _CYP2D6_ACTIVITY.get(d6_a1, 1.0)
        d6_score2 = _CYP2D6_ACTIVITY.get(d6_a2, 1.0)
        cyp2d6_as = d6_score1 + d6_score2

        # CYP2D6 phenotype
        if cyp2d6_as < 0.25:
            cyp2d6_phenotype = "Poor Metabolizer"
        elif cyp2d6_as < 1.0:
            cyp2d6_phenotype = "Intermediate Metabolizer"
        elif cyp2d6_as < 2.25:
            cyp2d6_phenotype = "Normal Metabolizer"
        else:
            cyp2d6_phenotype = "Ultrarapid Metabolizer"

        # ── CYP2C19 activity score (for secondary adjustment) ──
        c19_a1, c19_a2 = self._parse_diplotype(cyp2c19_diplotype)
        c19_score1 = _CYP2C19_ACTIVITY.get(c19_a1, 1.0)
        c19_score2 = _CYP2C19_ACTIVITY.get(c19_a2, 1.0)
        cyp2c19_as = c19_score1 + c19_score2

        if cyp2c19_as == 0.0:
            cyp2c19_phenotype = "Poor Metabolizer"
        elif cyp2c19_as < 2.0:
            cyp2c19_phenotype = "Intermediate Metabolizer"
        elif cyp2c19_as < 2.5:
            cyp2c19_phenotype = "Normal Metabolizer"
        elif cyp2c19_as < 3.0:
            cyp2c19_phenotype = "Rapid Metabolizer"
        else:
            cyp2c19_phenotype = "Ultrarapid Metabolizer"

        # ── CYP2D6-based dose recommendation ──
        is_tertiary = tca_class == "tertiary"
        warnings: List[str] = []
        tdm_required = False

        if cyp2d6_phenotype == "Ultrarapid Metabolizer":
            if is_tertiary:
                dose_adjustment = "avoid"
                rationale = (
                    "CYP2D6 UM rapidly metabolizes TCAs. Tertiary amines: "
                    "avoid — subtherapeutic parent and metabolite levels "
                    "expected. Secondary amines (nortriptyline, desipramine): "
                    "may attempt with mandatory TDM and dose titration to "
                    "therapeutic range."
                )
            else:
                dose_adjustment = "use_with_tdm"
                rationale = (
                    "CYP2D6 UM rapidly metabolizes TCAs. Secondary amines "
                    "(nortriptyline, desipramine): may attempt with mandatory "
                    "TDM and dose titration to therapeutic range."
                )
            tdm_required = True
            warnings = [
                "Avoid tertiary amine TCAs — treatment failure expected",
                "If secondary amine used, mandatory TDM",
                "Consider SSRI/SNRI alternative",
            ]
        elif cyp2d6_phenotype == "Poor Metabolizer":
            dose_adjustment = "reduce_50_percent"
            rationale = (
                "CYP2D6 PM causes 2-4 fold increase in TCA exposure. "
                "Reduce starting dose by 50% and titrate with TDM. "
                "Target therapeutic range."
            )
            tdm_required = True
            warnings = [
                "Reduce dose by 50%",
                "Mandatory TDM — target therapeutic range",
                "Risk of cardiac toxicity (QT prolongation), "
                "anticholinergic crisis, seizures",
                "Consider SSRI alternative if TDM unavailable",
            ]
        elif cyp2d6_phenotype == "Intermediate Metabolizer":
            dose_adjustment = "reduce_25_percent"
            rationale = (
                "CYP2D6 IM has moderately reduced TCA clearance. Consider "
                "25% dose reduction with TDM."
            )
            tdm_required = True
            warnings = [
                "Consider 25% dose reduction",
                "TDM recommended",
            ]
        else:  # Normal Metabolizer
            dose_adjustment = "none"
            rationale = (
                "Standard TCA dosing. Normal CYP2D6 metabolism."
            )
            tdm_required = False

        # ── CYP2C19 secondary adjustment (tertiary amines only) ──
        if is_tertiary:
            if cyp2c19_phenotype == "Poor Metabolizer":
                warnings.append(
                    "CYP2C19 PM further reduces demethylation of tertiary "
                    "amine TCA. Additional dose reduction may be needed."
                )
            elif cyp2c19_phenotype == "Ultrarapid Metabolizer":
                warnings.append(
                    "CYP2C19 UM increases demethylation; monitor metabolite "
                    "levels."
                )

        alternatives = ["sertraline", "venlafaxine", "mirtazapine"]
        tdm_target_range = _TCA_TDM_RANGES.get(drug_lower, "See drug-specific reference")

        result = {
            "drug": drug,
            "tca_class": tca_class,
            "cyp2d6_diplotype": cyp2d6_diplotype,
            "cyp2d6_activity_score": cyp2d6_as,
            "cyp2d6_phenotype": cyp2d6_phenotype,
            "cyp2c19_diplotype": cyp2c19_diplotype,
            "cyp2c19_phenotype": cyp2c19_phenotype,
            "dose_adjustment": dose_adjustment,
            "rationale": rationale,
            "warnings": warnings,
            "tdm_required": tdm_required,
            "tdm_target_range": tdm_target_range,
            "alternatives": alternatives,
            "affected_drugs": list(_TCA_CLASS.keys()),
            "reference": (
                "CPIC Guideline for Tricyclic Antidepressants and "
                "CYP2D6/CYP2C19, 2016 Update"
            ),
        }

        logger.info(
            "TCA dosing: %s %s (CYP2D6 %s AS %.2f %s, CYP2C19 %s %s)",
            drug, dose_adjustment, cyp2d6_diplotype, cyp2d6_as,
            cyp2d6_phenotype, cyp2c19_diplotype, cyp2c19_phenotype,
        )
        return result

    # ═════════════════════════════════════════════════════════════════════
    # Private helpers
    # ═════════════════════════════════════════════════════════════════════

    @staticmethod
    def _parse_diplotype(diplotype: str) -> Tuple[str, str]:
        """Split a diplotype string into two alleles."""
        parts = diplotype.split("/")
        if len(parts) == 2:
            return parts[0].strip(), parts[1].strip()
        return "*1", "*1"

    @staticmethod
    def _normalize_diplotype(diplotype: str) -> str:
        """Normalize diplotype to canonical form (lower number first)."""
        parts = diplotype.split("/")
        if len(parts) != 2:
            return diplotype
        a, b = parts[0].strip(), parts[1].strip()

        def _sort_key(s: str) -> float:
            import re as _re
            m = _re.match(r"\*(\d+)", s)
            return float(m.group(1)) if m else 999.0

        if _sort_key(a) > _sort_key(b):
            return f"{b}/{a}"
        return f"{a}/{b}"

    @staticmethod
    def _is_vkorc1_ag(genotype: str) -> bool:
        """Check if VKORC1 genotype is heterozygous (AG)."""
        g = genotype.lower().replace(" ", "")
        return any(x in g for x in [
            "*1/-1639g>a", "-1639g>a/*1",
            "g/a", "a/g", "ag", "ga",
            "*1/-1639", "-1639/*1",
        ]) and not DosingCalculator._is_vkorc1_aa(genotype)

    @staticmethod
    def _is_vkorc1_aa(genotype: str) -> bool:
        """Check if VKORC1 genotype is homozygous variant (AA)."""
        g = genotype.lower().replace(" ", "")
        return any(x in g for x in [
            "-1639g>a/-1639g>a",
            "a/a",
            "aa",
        ]) and "*1/*1" not in g

    @staticmethod
    def _warfarin_dose_category(weekly_dose: float) -> str:
        """Categorize warfarin dose as low, standard, or high."""
        if weekly_dose < 21:
            return "Low Dose (sensitive)"
        elif weekly_dose <= 49:
            return "Standard Dose"
        else:
            return "High Dose (resistant)"

    @staticmethod
    def _warfarin_clinical_notes(
        cyp2c9: str, vkorc1: str, dose: float, amiodarone: bool,
        enzyme_inducer: bool = False,
    ) -> List[str]:
        """Generate clinical notes for warfarin dosing."""
        notes = []

        # CYP2C9 variant notes
        if cyp2c9 in ("*2/*2", "*2/*3", "*3/*3"):
            notes.append(
                "CYP2C9 poor metabolizer — significantly increased bleeding "
                "risk. Longer time to stable dose expected. Consider more "
                "frequent INR monitoring during initiation."
            )
        elif cyp2c9 in ("*1/*2", "*1/*3"):
            notes.append(
                "CYP2C9 intermediate metabolizer — moderate dose reduction "
                "expected. Standard INR monitoring with closer follow-up "
                "during initiation."
            )

        # VKORC1 variant notes
        if DosingCalculator._is_vkorc1_aa(vkorc1):
            notes.append(
                "VKORC1 AA genotype — significantly reduced warfarin "
                "requirement. Low VKORC1 expression results in greater "
                "vitamin K epoxide reductase sensitivity."
            )

        # Dose-specific notes
        if dose < 14:
            notes.append(
                "Predicted dose <2 mg/day — consider alternative "
                "anticoagulant (DOAC) if clinically appropriate, as very "
                "low warfarin doses are difficult to manage."
            )
        elif dose > 63:
            notes.append(
                "Predicted dose >9 mg/day — high warfarin resistance. "
                "Verify CYP4F2 genotype. Ensure adherence and dietary "
                "vitamin K consistency. Consider DOAC alternative."
            )

        if amiodarone:
            notes.append(
                "Amiodarone co-administration — inhibits warfarin "
                "metabolism. The algorithm already accounts for this, but "
                "monitor closely if amiodarone is started/stopped."
            )

        if enzyme_inducer:
            notes.append(
                "Enzyme inducer co-administration (rifampin, carbamazepine, "
                "phenytoin, or phenobarbital) — increases warfarin metabolism. "
                "Higher doses typically required. Monitor INR closely if "
                "inducer is started, stopped, or dose-changed."
            )

        if not notes:
            notes.append(
                "Standard pharmacogenomic warfarin dosing. Initiate at "
                "predicted dose and titrate to target INR (typically 2-3). "
                "Monitor INR at least weekly during initiation."
            )

        return notes

    @staticmethod
    def _thiopurine_phenotype(activity_score: float, gene: str) -> str:
        """Determine thiopurine-related phenotype from activity score."""
        if activity_score == 0.0:
            return f"{gene} Deficient"
        elif activity_score <= 1.0:
            return f"{gene} Intermediate"
        else:
            return f"{gene} Normal"
