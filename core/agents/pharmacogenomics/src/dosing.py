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
from typing import Any, Dict, List, Tuple

logger = logging.getLogger(__name__)


class IndeterminateGenotype(ValueError):
    """A genotype no CPIC table can score, so no dose can be recommended.

    CPIC expresses this as a distinct phenotype value — Indeterminate —
    whose recommendation is "No recommendation".  It is not a normal
    metabolizer, not wild-type, and not a dose.
    """

    def __init__(self, genotype: str, reason: str, alleles: List[str]) -> None:
        super().__init__(reason)
        self.genotype = genotype
        self.reason = reason
        self.alleles = alleles


class CapabilityWithdrawn(RuntimeError):
    """A calculator deliberately withheld for a correctness reason, not a failure.

    Distinct from an exception, and the distinction is the point. A caller that sees HTTP 500
    "Internal processing error" reasonably concludes the server is broken and may retry, report
    a bug, or treat the result as transient. None of those is right for a calculator that has
    been withdrawn because its coefficients are known to be wrong: the correct response is
    permanent-until-fixed and non-retryable, and it should carry the reason.

    Added 2026-11-04. The IWPC refusal below raised a bare RuntimeError, and the API route
    caught `Exception` and replaced it with "Internal processing error" — so a carefully
    written explanation of a race-specific dosing error reached the server log and never the
    caller. The same shape as the rest of this review: the honest words were somewhere nobody
    reads.
    """


# ═════════════════════════════════════════════════════════════════════════════
# Activity score lookups (shared with PhenotypeTranslator but duplicated
# here for module independence)
# ═════════════════════════════════════════════════════════════════════════════

# DPYD activity values, from CPIC's maintained allele-functionality table
# (api.cpicpgx.org/v1/allele?genesymbol=eq.DPYD, retrieved 2026-11-03).
#
# Three of these were wrong, in both directions:
#   *7  (c.295_298delTCAT) was 0.5 "decreased"; CPIC: NO function, 0.0
#   *10 (c.2983G>T)        was 0.5 "decreased"; CPIC: NO function, 0.0
#   *9A (c.85T>C)          was 0.5 "decreased"; CPIC: NORMAL function, 1.0
# *7/*7 -- a DPD-deficient patient CPIC tells you to avoid fluoropyrimidines
# in -- therefore scored AS 1.0 and was handed a half dose of 5-FU; and
# *1/*9A, a variant with no effect on DPD activity and a ~20% allele
# frequency in Europeans, had its curative-intent adjuvant chemotherapy
# reduced.
_DPYD_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,       # Reference / normal function
    "*2A": 0.0,      # c.1905+1G>A, no function
    "*3": 0.0,       # c.1898delC, no function
    "*7": 0.0,       # c.295_298delTCAT, NO function (was 0.5)
    "*8": 0.0,       # c.703C>T, no function
    "*9A": 1.0,      # c.85T>C, NORMAL function (was 0.5)
    "*10": 0.0,      # c.2983G>T, NO function (was 0.5)
    "*12": 0.0,      # c.1156G>T, no function
    "*13": 0.0,      # c.1679T>G, no function
    "c.2846A>T": 0.5,  # decreased function
    "HapB3": 0.5,    # c.1129-5923C>G, c.1236G>A — ONE allele, decreased function
}

# TPMT activity values, from CPIC's maintained allele-functionality table
# (api.cpicpgx.org/v1/allele?genesymbol=eq.TPMT, retrieved 2026-11-03).
# *8 is CPIC's only decreased-function TPMT allele; it was absent, so
# *1/*8 and *2/*8 were unscoreable.
_TPMT_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,
    "*2": 0.0,
    "*3A": 0.0,
    "*3B": 0.0,
    "*3C": 0.0,
    "*4": 0.0,
    "*8": 0.5,   # decreased function
}

# NUDT15 activity values, from CPIC's MAINTAINED allele-functionality table
# (api.cpicpgx.org/v1/allele?genesymbol=eq.NUDT15, retrieved 2026-11-03).
#
# This table previously carried *1/*2/*3 only, sourced to Relling 2019
# Table 1, which in 2018 called *4-*9 "unclear function". **CPIC has since
# assigned them**, and a guideline paper is not the same artefact as the
# table the issuing body maintains: for allele function the maintained
# table governs. CPIC now says *4, *6, *7, *9 are NO function and *5, *8
# are DECREASED function, and its diplotype table follows:
#     *1/*4 -> Intermediate Metabolizer      *2/*4 -> Poor Metabolizer
#     *1/*6 -> Intermediate Metabolizer      *3/*6 -> Poor Metabolizer
#     *1/*5 -> Normal Metabolizer            *6/*6 -> Poor Metabolizer
# Treating them as unscoreable sent a NUDT15 poor metabolizer away with
# "dose on clinical grounds" -- for a thiopurine in a child, a full
# weight-based dose and profound prolonged myelosuppression.
#
# Indeterminate is reserved for the alleles CPIC still calls uncertain,
# below.
_NUDT15_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,                                         # normal function
    "*2": 0.0, "*3": 0.0, "*4": 0.0, "*6": 0.0,        # no function
    "*7": 0.0, "*9": 0.0, "*10": 0.0, "*13": 0.0,      # no function
    "*14": 0.0, "*17": 0.0, "*18": 0.0, "*19": 0.0,    # no function
    "*5": 0.5, "*8": 0.5,                              # decreased function
}

# CPIC's current "Uncertain function" NUDT15 alleles — these, and only
# these, make a NUDT15 diplotype Indeterminate
# (api.cpicpgx.org/v1/diplotype?genesymbol=eq.NUDT15 gives Indeterminate
# for every diplotype carrying one). Previously held *4-*9.
_NUDT15_UNCERTAIN: frozenset = frozenset({
    "*11", "*12", "*15", "*16", "*20", "*21", "*22", "*23",
})

# CYP3A5, from api.cpicpgx.org/v1/allele?genesymbol=eq.CYP3A5
# (retrieved 2026-11-03): *1 normal; *3, *6, *7 no function; *8 and *9
# Unknown function; *10 and *11 carry no assignment at all.
_CYP3A5_ACTIVITY: Dict[str, float] = {
    "*1": 1.0,  # expresser
    "*3": 0.0,  # non-expresser (most common)
    "*6": 0.0,  # non-expresser
    "*7": 0.0,  # non-expresser
}

_CYP3A5_UNCERTAIN: frozenset = frozenset({"*8", "*9", "*10", "*11"})

# CYP2C19, from api.cpicpgx.org/v1/allele?genesymbol=eq.CYP2C19
# (retrieved 2026-11-03). Two corrections and sixteen additions:
#   *9, *10 were scored 0.0 (no function); CPIC calls both DECREASED, 0.5
#   *22, *24, *35, *36, *37 are no-function alleles that were absent, so
#     `.get(allele, 1.0)` scored them as fully functional -- CYP2C19
#     *1/*35 came back "Normal Metabolizer, standard clopidogrel" where
#     CPIC's diplotype table says Intermediate Metabolizer and "avoid
#     standard dose clopidogrel".
# Verified: these values plus _CYP2C19_* thresholds reproduce CPIC's
# phenotype for all 351 CYP2C19 diplotypes built from assigned alleles.
_CYP2C19_ACTIVITY: Dict[str, float] = {
    "*1": 1.0, "*11": 1.0, "*13": 1.0, "*15": 1.0, "*18": 1.0, "*28": 1.0,
    "*38": 1.0,
    "*2": 0.0, "*3": 0.0, "*4": 0.0, "*5": 0.0, "*6": 0.0, "*7": 0.0,
    "*8": 0.0, "*22": 0.0, "*24": 0.0, "*35": 0.0, "*36": 0.0, "*37": 0.0,
    "*9": 0.5, "*10": 0.5, "*16": 0.5, "*19": 0.5, "*25": 0.5, "*26": 0.5,
    "*17": 1.5,
}

_CYP2C19_UNCERTAIN: frozenset = frozenset({
    "*12", "*14", "*23", "*29", "*30", "*31", "*32", "*33", "*34", "*39",
})

# SLCO1B1 is a function-class gene, not an activity-score one. This was a
# hand-listed table of EIGHT diplotypes; CPIC's own diplotype table has 741,
# and two of the eight -- "*1/*17" (Decreased Function) and "*5/*17" (Poor
# Function) -- asserted a phenotype for an allele that appears in neither
# CPIC's SLCO1B1 allele-functionality table nor its diplotype table, so they
# are gone (recorded as an open question in CLINICAL_REVIEW.md).
#
# Replaced by CPIC's allele functions plus CPIC's own function-pair rule,
# which covers every diplotype and adds the two categories the hand-list
# lacked (Possible Decreased Function, Increased Function).
# api.cpicpgx.org/v1/allele?genesymbol=eq.SLCO1B1 and
# api.cpicpgx.org/v1/diplotype?genesymbol=eq.SLCO1B1, retrieved 2026-11-03.
_SLCO1B1_ALLELE_FUNCTION: Dict[str, str] = {
    "*1": "normal", "*37": "normal",
    "*14": "increased", "*20": "increased",
    "*5": "none", "*9": "none", "*15": "none", "*23": "none",
    "*31": "none", "*39": "none", "*41": "none", "*45": "none",
    "*47": "none", "*48": "none", "*49": "none",
}

# CPIC calls these "Uncertain function"; a diplotype carrying one is either
# "Possible Decreased Function" (paired with a no-function allele) or
# Indeterminate.
_SLCO1B1_UNCERTAIN: frozenset = frozenset({
    "*2", "*4", "*6", "*8", "*10", "*11", "*16", "*19", "*24", "*25",
    "*26", "*27", "*28", "*29", "*30", "*32", "*33", "*34", "*36", "*40",
    "*42", "*43", "*44", "*50", "*51", "*52", "*53", "*54", "*55", "*56",
    "*57", "*58",
})

# CPIC's function-pair rule, read off its 741-row diplotype table. Verified:
# reproduces every row built from alleles CPIC assigns a function to.
_SLCO1B1_PAIR_PHENOTYPE: Dict[frozenset, str] = {
    frozenset({"none"}): "Poor Function",
    frozenset({"none", "normal"}): "Decreased Function",
    frozenset({"none", "increased"}): "Decreased Function",
    frozenset({"none", "uncertain"}): "Possible Decreased Function",
    frozenset({"normal"}): "Normal Function",
    frozenset({"normal", "increased"}): "Normal Function",
    frozenset({"increased"}): "Increased Function",
}

_SLCO1B1_MYOPATHY_RISK: Dict[str, str] = {
    "Increased Function": "average or below (population baseline)",
    "Normal Function": "average (population baseline)",
    "Possible Decreased Function": "possibly increased — function uncertain",
    "Decreased Function": "increased (4.5-fold with 80mg dose)",
    "Poor Function": "high (17-fold increased with standard dose)",
}

# CYP2C9, from api.cpicpgx.org/v1/allele?genesymbol=eq.CYP2C9
# (retrieved 2026-11-03).
#   *5 was scored 0.0 (no function); CPIC calls it DECREASED, 0.5, and its
#     diplotype table makes *5/*5 an Intermediate Metabolizer (AS 1.0), not
#     a Poor one -- the repo gave a *5/*5 patient (an allele essentially
#     restricted to African-ancestry populations) a 50% phenytoin
#     maintenance reduction where CPIC says 25%.
#   *4 and *13 and eleven other assigned alleles were absent, so
#     `.get(allele, 1.0)` scored them as fully functional: *1/*13 (CPIC
#     Intermediate, 25% less phenytoin) came back "Normal Metabolizer,
#     standard phenytoin dosing".
# Verified: these values reproduce CPIC's phenotype for all 2,850 CYP2C9
# diplotypes built from assigned alleles.
_CYP2C9_ACTIVITY: Dict[str, float] = {
    "*1": 1.0, "*9": 1.0,
    "*2": 0.5, "*4": 0.5, "*5": 0.5, "*8": 0.5, "*11": 0.5, "*12": 0.5,
    "*14": 0.5, "*16": 0.5, "*23": 0.5, "*26": 0.5, "*28": 0.5, "*29": 0.5,
    "*30": 0.5, "*31": 0.5, "*37": 0.5, "*38": 0.5, "*44": 0.5, "*46": 0.5,
    "*50": 0.5, "*55": 0.5, "*61": 0.5,
    "*3": 0.0, "*6": 0.0, "*13": 0.0, "*15": 0.0, "*24": 0.0, "*25": 0.0,
    "*33": 0.0, "*35": 0.0, "*39": 0.0, "*42": 0.0, "*43": 0.0, "*45": 0.0,
    "*52": 0.0,
}

# CPIC calls these "Uncertain function" or assigns no function at all;
# a diplotype carrying one is Indeterminate, not a metabolizer class.
_CYP2C9_UNCERTAIN: frozenset = frozenset({
    "*7", "*10", "*17", "*18", "*19", "*20", "*21", "*22", "*27", "*32",
    "*34", "*36", "*40", "*41", "*47", "*48", "*49", "*51", "*53", "*54",
    "*56", "*57", "*58", "*59", "*60", "*62", "*63", "*64", "*65", "*66",
    "*67", "*68", "*69", "*70", "*71", "*72", "*73", "*74", "*75",
})

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

        # ─────────────────────────────────────────────────────────────────────────────────
        # THIS DOSE IS NOT RETURNED. The coefficient block above is mis-transcribed, and the
        # error is race-specific, so it must not reach a prescriber.
        #
        # Every number above IS an IWPC coefficient (Klein et al., "Estimation of the Warfarin
        # Dose with Clinical and Pharmacogenetic Data", N Engl J Med 2009;360:753-764), but the
        # column is shifted: the published algorithm assigns
        #     -0.2188  CYP2C9 genotype UNKNOWN      (here: multiplying cyp4f2_1_3)
        #     -0.2760  Black or African American    (here: multiplying cyp4f2_3_3)
        #     +1.1816  enzyme inducer status        (here: multiplying african_american)
        # and IWPC has no CYP4F2 term and no smoking term at all, so `cyp4f2_*` and
        # `smoker_val` are not in the model being cited. The amiodarone coefficient is also not
        # the published one.
        #
        # The consequence is the reason this is a refusal rather than a TODO: a Black patient
        # picks up +1.1816 on the square-root-of-weekly-dose scale where the algorithm applies
        # -0.2760, which inflates the predicted dose by roughly 1.5-1.7x -- a race-specific
        # overdose of a narrow-therapeutic-index anticoagulant, in the one direction that
        # causes a bleed. No other patient group is affected the same way.
        #
        # The fix is a ten-minute re-transcription, and it is deliberately NOT done here,
        # because this repo's rule is that a clinical value changes only against a primary
        # source in hand. Supplementary Appendix 1 of Klein 2009 (and the warfarindosing.org
        # workbook) carries the table; neither could be retrieved tonight. Restoring it also
        # means adding the three "genotype unknown" terms the block never had.
        # Recorded as CLINICAL_REVIEW.md §21 with the full mapping.
        raise CapabilityWithdrawn(
            "IWPC warfarin dosing is disabled: the coefficient block in this function is"
            "mis-transcribed from Klein 2009 and the error is race-specific (a Black patient "
            "receives the enzyme-inducer coefficient, inflating the predicted dose ~1.5-1.7x). "
            "No dose is returned until the published coefficients are restored from "
            "Supplementary Appendix 1. See CLINICAL_REVIEW.md §21."
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

        CPIC's tacrolimus recommendation table gives ONE recommendation to
        CYP3A5 Normal Metabolizer, Possible Intermediate Metabolizer and
        Intermediate Metabolizer alike: "Increase starting dose 1.5 to 2
        times recommended starting dose. Total starting dose should not
        exceed 0.3 mg/kg/day. Use therapeutic drug monitoring to guide dose
        adjustments."  Poor Metabolizer: "Initiate therapy with standard
        recommended dose."  Indeterminate: "No recommendation."
        (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:42316, retrieved
        2026-11-03.)

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
        # `.get(allele, 0.0)` here was the opposite default from the rest of
        # the module, and it called an unrecognised allele a NON-EXPRESSER
        # and gave it the standard 0.15 mg/kg/day. A true expresser dosed as
        # a non-expresser is subtherapeutic tacrolimus in the first week
        # post-transplant — the window in which acute rejection happens.
        # CPIC's answer for Indeterminate is "No recommendation".
        try:
            allele1, allele2, score1, score2 = self._score_diplotype(
                "CYP3A5", cyp3a5_diplotype, _CYP3A5_ACTIVITY, _CYP3A5_UNCERTAIN)
        except IndeterminateGenotype as exc:
            return {
                "dose_recommendation": (
                    f"No genotype-guided tacrolimus dose can be given: "
                    f"{exc.reason} An unrecognised CYP3A5 allele must NOT be "
                    f"dosed as a non-expresser — a true expresser given the "
                    f"standard 0.15 mg/kg/day is subtherapeutic in the first "
                    f"week post-transplant. Dose per local protocol with "
                    f"trough monitoring from day 1."
                ),
                "mg_per_kg_per_day": None,
                "total_daily_mg": None,
                "per_dose_mg": None,
                "frequency": "BID (every 12 hours)",
                "phenotype": "Indeterminate",
                "cyp3a5_diplotype": cyp3a5_diplotype,
                "activity_score": None,
                "indeterminate_alleles": exc.alleles,
                "rationale": exc.reason,
                "monitoring": (
                    "Therapeutic drug monitoring (trough levels) required."
                ),
            }
        total_score = score1 + score2

        # Determine phenotype and dosing
        if total_score >= 0.5:
            # CPIC gives Normal, Possible Intermediate and Intermediate
            # Metabolizers the SAME recommendation. The code used to split
            # *1/*1 (0.30 mg/kg/day) from *1/*3 ("CYP3A5 Intermediate
            # Expresser", 0.25 mg/kg/day, "~1.25-1.5x standard") — an
            # invented clinical distinction presented under a CPIC heading,
            # contradicted by this function's own docstring.
            phenotype = (
                "CYP3A5 Expresser" if total_score >= 1.5
                else "CYP3A5 Intermediate Expresser"
            )
            mg_per_kg = 0.30
            recommendation = (
                f"{phenotype}: increase starting dose 1.5 to 2 times the "
                f"recommended starting dose. Total starting dose should not "
                f"exceed 0.3 mg/kg/day. Use therapeutic drug monitoring to "
                f"guide dose adjustments (target 10-15 ng/mL early "
                f"post-transplant). CPIC gives Normal, Possible Intermediate "
                f"and Intermediate Metabolizers this same recommendation."
            )
            rationale = (
                "Functional CYP3A5 protein significantly increases tacrolimus "
                "clearance; expressers require higher doses to reach target "
                "trough concentrations. CPIC draws no dose distinction "
                "between one and two functional alleles."
            )
        else:
            # Non-expresser (*3/*3, *3/*6, *6/*6)
            phenotype = "CYP3A5 Non-expresser"
            mg_per_kg = 0.15
            recommendation = (
                "CYP3A5 Non-expresser (CPIC Poor Metabolizer): initiate "
                "therapy with the standard recommended dose, 0.15 mg/kg/day. "
                "No functional CYP3A5 — metabolism primarily via CYP3A4. "
                "Use therapeutic drug monitoring to guide dose adjustments."
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

        CPIC's own fluorouracil and capecitabine recommendation rows, keyed
        on DPYD activity score, verbatim
        (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:4492 and
        =eq.RxNorm:194000, retrieved 2026-11-03):

            2.0  "Based on genotype, there is no indication to change dose
                  or therapy."
            1.5  "Reduce starting dose by 50% followed by titration of dose
                  based on toxicity or therapeutic drug monitoring."
            1.0  as 1.5, plus "Patients with the c.2846A>T/c.2846A>T
                  genotype may require >50% reduction in starting dose."
            0.5  "AVOID use of 5-fluorouracil or 5-fluorouracil
                  prodrug-based regimens. In the event ... alternative
                  agents are not considered a suitable therapeutic option,
                  5-fluorouracil should be administered at a strongly
                  reduced dose with early therapeutic drug monitoring."
            0.0  "Avoid use of 5-fluorouracil or 5-fluorouracil
                  prodrug-based regimens."

        CPIC does NOT grade 1.5 and 1.0 differently, and it calls activity
        score 0.5 a **Poor Metabolizer** for whom fluoropyrimidines are to
        be avoided.

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
        try:
            allele1, allele2, score1, score2 = self._score_diplotype(
                "DPYD", dpyd_diplotype, _DPYD_ACTIVITY)
        except IndeterminateGenotype as exc:
            return {
                "activity_score": None,
                "dose_reduction_percent": None,
                "recommendation": (
                    f"No DPYD-guided fluoropyrimidine dose can be given: "
                    f"{exc.reason} Do NOT give a standard 5-FU or capecitabine "
                    f"dose on the strength of an unscoreable DPYD result — "
                    f"check the genotype against CPIC's DPYD "
                    f"allele-functionality table, or consider DPD phenotyping "
                    f"(uracil/dihydrouracil ratio) before the first cycle."
                ),
                "phenotype": "Indeterminate",
                "dpyd_diplotype": dpyd_diplotype,
                "allele_scores": {},
                "indeterminate_alleles": exc.alleles,
                "affected_drugs": [
                    "5-fluorouracil (5-FU)", "capecitabine", "tegafur",
                ],
                "monitoring": (
                    "CBC with differential prior to each cycle. Monitor for "
                    "mucositis, diarrhea, hand-foot syndrome, neurotoxicity."
                ),
            }
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
                "AVOID use of 5-fluorouracil or 5-fluorouracil prodrug-based "
                "regimens. Complete DPD deficiency — life-threatening and "
                "potentially fatal toxicity (severe mucositis, "
                "myelosuppression, neurotoxicity, hand-foot syndrome). "
                "Select an alternative non-fluoropyrimidine regimen."
            )
        elif activity_score <= 0.5:
            # CLINICAL_REVIEW §4, SETTLED 2026-11-03. This branch said "DPD
            # Intermediate Activity (Partial Deficiency)" and a 50%
            # reduction, while pgx_pipeline.py called AS 0.5 "DPD
            # Deficient"; knowledge.py gave a third answer again ("25-50%
            # dose with intensive monitoring"). CPIC settles it: AS 0.5 is a
            # **Poor Metabolizer** and the recommendation is AVOID.
            phenotype = "DPD Deficient (Poor Metabolizer, activity score 0.5)"
            dose_reduction = 100
            recommendation = (
                "AVOID use of 5-fluorouracil or 5-fluorouracil prodrug-based "
                "regimens. Activity score 0.5 is a CPIC Poor Metabolizer. In "
                "the event that, based on clinical advice, alternative agents "
                "are not considered a suitable therapeutic option, "
                "5-fluorouracil should be administered at a strongly reduced "
                "dose with early therapeutic drug monitoring."
            )
        elif activity_score <= 1.0:
            phenotype = "DPD Intermediate Activity"
            dose_reduction = 50
            recommendation = (
                "Reduce starting dose by 50%, followed by titration of dose "
                "based on toxicity or therapeutic drug monitoring (if "
                "available). Patients with the c.2846A>T/c.2846A>T genotype "
                "may require a >50% reduction in starting dose."
            )
        elif activity_score <= 1.5:
            # WAS: 25%, labelled "(Mild)". CPIC gives activity score 1.5 the
            # SAME 50% reduction as 1.0 — and *1/c.2846A>T is the single most
            # frequently reported actionable DPYD genotype in
            # European-ancestry patients. A 25% reduction leaves such a
            # patient on 75% of a standard 5-FU dose where CPIC calls for
            # 50%; DPD heterozygotes near full dose are the population the
            # fatal-toxicity case reports come from.
            phenotype = "DPD Intermediate Activity"
            dose_reduction = 50
            recommendation = (
                "Reduce starting dose by 50%, followed by titration of dose "
                "based on toxicity or therapeutic drug monitoring (if "
                "available). Activity score 1.5 indicates one "
                "decreased-function allele; CPIC gives it the same 50% "
                "starting-dose reduction as activity score 1.0."
            )
        else:
            phenotype = "DPD Normal Activity"
            dose_reduction = 0
            recommendation = (
                "Based on genotype, there is no indication to change dose or "
                "therapy. Use label-recommended dosage and administration. "
                "Activity score 2.0 (normal). Standard toxicity monitoring "
                "applies."
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
        nudt15_diplotype: str = "",
    ) -> Dict[str, Any]:
        """Calculate TPMT + NUDT15 genotype-guided thiopurine dose.

        CPIC's recommendation is a TWO-GENE table over the (TPMT phenotype,
        NUDT15 phenotype) pair, not a function of the lower of the two
        activity scores — see the comment at the combined branch below for
        the rows and the primary source.

        Parameters
        ----------
        tpmt_diplotype : str
            e.g. "*1/*1", "*1/*3A", "*3A/*3C"
        nudt15_diplotype : str
            e.g. "*1/*1", "*1/*3", "*3/*3".  The default used to be
            "*1/*1", so a TPMT-only panel silently reported NUDT15 Normal
            and dosed on it (pgx audit §20h).  It is now empty, which
            returns CPIC's Indeterminate: "Neither TPMT or NUDT15 phenotype
            could be assigned ... if thiopurines are required and either
            TPMT or NUDT15 status is unknown, monitor closely for toxicity."

        Returns
        -------
        dict
            dose_recommendation, tpmt_phenotype, nudt15_phenotype,
            dose_percent_of_standard, dose_percent_range_of_standard,
            recommendations_by_drug, affected_drugs.
        """
        # TPMT activity. Same hazard as NUDT15 below: this table names six alleles, CPIC's
        # allele-functionality table names dozens, and `.get(allele, 1.0)` scored every one it
        # had never seen as fully functional -- a no-function allele reported as a normal
        # metabolizer, with a full thiopurine dose behind it.
        try:
            t_a1, t_a2 = self._parse_diplotype(tpmt_diplotype)
        except IndeterminateGenotype as exc:
            t_a1 = t_a2 = exc.genotype
        t_unscoreable = [a for a in (t_a1, t_a2) if a not in _TPMT_ACTIVITY]
        if t_unscoreable:
            t_unscoreable = list(dict.fromkeys(t_unscoreable))
            return {
                "drug_class": "thiopurines",
                "tpmt_diplotype": tpmt_diplotype,
                "nudt15_diplotype": nudt15_diplotype,
                "tpmt_phenotype": "TPMT Indeterminate",
                "nudt15_phenotype": None,
                "dose_recommendation": (
                    f"No genotype-guided dose can be given: TPMT {', '.join(t_unscoreable)} is "
                    "not in this agent's allele table, so its function is unknown. Check the "
                    "CPIC TPMT allele-functionality table before dosing, and do not read this "
                    "as a normal-metabolizer result."
                ),
                "dose_percent_of_standard": None,
                "dose_percent_range_of_standard": None,
                "recommendations_by_drug": {},
                "indeterminate_alleles": t_unscoreable,
                "affected_drugs": [
                    "azathioprine", "6-mercaptopurine (6-MP)", "6-thioguanine (6-TG)",
                ],
                "guideline": (
                    "CPIC Thiopurine Dosing Based on TPMT and NUDT15 Genotypes, 2018 Update "
                    "(Relling et al., Clin Pharmacol Ther 2019;105:1095-1105), Table 1"
                ),
            }
        tpmt_score = _TPMT_ACTIVITY[t_a1] + _TPMT_ACTIVITY[t_a2]

        # NUDT15 activity. An uncertain-function or unrecognised allele cannot be scored, and
        # defaulting it to 1.0 would recommend a full dose on an allele CPIC declines to call --
        # the direction of error that ends in myelosuppression. Return Indeterminate instead and
        # say which allele caused it.
        try:
            n_a1, n_a2 = self._parse_diplotype(nudt15_diplotype)
        except IndeterminateGenotype as exc:
            n_a1 = n_a2 = exc.genotype
        n_unscoreable = [a for a in (n_a1, n_a2)
                         if a in _NUDT15_UNCERTAIN or a not in _NUDT15_ACTIVITY]
        if n_unscoreable:
            n_unscoreable = list(dict.fromkeys(n_unscoreable))
            return {
                "drug_class": "thiopurines",
                "tpmt_diplotype": tpmt_diplotype,
                "nudt15_diplotype": nudt15_diplotype,
                "tpmt_phenotype": self._thiopurine_phenotype(tpmt_score, "TPMT"),
                "nudt15_phenotype": "NUDT15 Indeterminate",
                "dose_recommendation": (
                    "No genotype-guided dose can be given: NUDT15 "
                    f"{', '.join(n_unscoreable)} has uncertain function in, or is absent "
                    "from, CPIC's current allele-functionality table, which assigns "
                    "Indeterminate rather than a "
                    "metabolizer class. CPIC's own answer here is: if thiopurines are "
                    "required and NUDT15 status is unknown, monitor closely for toxicity; "
                    "consider evaluating TPMT erythrocyte activity. Do NOT read this as a "
                    "normal-metabolizer result and do NOT give a full weight-based dose on "
                    "the strength of it."
                ),
                "dose_percent_of_standard": None,
                "dose_percent_range_of_standard": None,
                "recommendations_by_drug": {},
                "indeterminate_alleles": n_unscoreable,
                "affected_drugs": [
                    "azathioprine", "6-mercaptopurine (6-MP)", "6-thioguanine (6-TG)",
                ],
                "guideline": (
                    "CPIC Thiopurine Dosing Based on TPMT and NUDT15 Genotypes, 2018 Update "
                    "(Relling et al., Clin Pharmacol Ther 2019;105:1095-1105), Table 1"
                ),
            }
        nudt15_score = _NUDT15_ACTIVITY[n_a1] + _NUDT15_ACTIVITY[n_a2]

        # Determine individual phenotypes
        tpmt_phenotype = self._thiopurine_phenotype(tpmt_score, "TPMT")
        nudt15_phenotype = self._thiopurine_phenotype(nudt15_score, "NUDT15")

        # CPIC's thiopurine recommendation is a TWO-GENE table, keyed on the
        # (TPMT phenotype, NUDT15 phenotype) pair. This used to be
        # `combined_score = min(tpmt_score, nudt15_score)`, and min() cannot
        # express a two-gene recommendation: two intermediate metabolizers
        # looked identical to one. CPIC's own azathioprine rows
        # (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:1256,
        # retrieved 2026-11-03):
        #
        #   TPMT IM + NUDT15 NM  -> "reduced starting doses (30-80% ...)"
        #   TPMT NM + NUDT15 IM  -> "reduced starting doses (30-80% ...)"
        #   TPMT IM + NUDT15 IM  -> "reduced starting doses (20%-50% ...)"
        #   either gene PM       -> "Consider alternative nonthiopurine
        #                            immunosuppressant therapy."
        #
        # The old code gave the double heterozygote "30-70% of standard,
        # start at lower end (30-50%)" with dose_percent 50 -- a headline
        # range running 20 points above CPIC's ceiling, and a structured
        # field sitting exactly on that ceiling with titrate-upward advice.
        # Double heterozygotes have additively reduced thiopurine
        # clearance; starting at the high end is a predictable route to
        # grade 3-4 leukopenia at week 2-4.
        #
        # The 10%/thrice-weekly figure the Poor branch used to carry for
        # azathioprine belongs to MERCAPTOPURINE and THIOGUANINE
        # (?drugid=eq.RxNorm:103 -- "Reduce starting dose by 10-fold and
        # reduce frequency to thrice weekly"). CPIC's azathioprine answer
        # for a poor metabolizer at either locus is an alternative drug, not
        # a reduced azathioprine dose. This also settles pharmacogenomics
        # open question #2.
        tpmt_class = tpmt_phenotype.rsplit(" ", 1)[-1]      # Deficient/Intermediate/Normal
        nudt15_class = nudt15_phenotype.rsplit(" ", 1)[-1]
        classes = {tpmt_class, nudt15_class}

        affected_drugs = [
            "azathioprine",
            "6-mercaptopurine (6-MP)",
            "6-thioguanine (6-TG)",
        ]

        _ALTERNATIVE = (
            "Consider alternative nonthiopurine immunosuppressant therapy."
        )
        _MP_PM = (
            "For malignancy: initiate therapy with drastically reduced "
            "starting doses. Reduce starting dose by 10-fold and reduce "
            "frequency to thrice weekly instead of daily (e.g. 10 mg/m2/day "
            "given 3 days/week). Allow 4-6 weeks of stable dosing to reach "
            "steady state after each adjustment."
        )

        if "Deficient" in classes:
            which = " and ".join(
                g for g, c in (("TPMT", tpmt_class), ("NUDT15", nudt15_class))
                if c == "Deficient"
            )
            dose_range = None
            dose_percent = None
            recommendations_by_drug = {
                "azathioprine": _ALTERNATIVE,
                "6-mercaptopurine (6-MP)": _MP_PM,
                "6-thioguanine (6-TG)": _MP_PM,
            }
            dose_recommendation = (
                f"{which} poor metabolizer — extremely high risk of "
                f"life-threatening myelosuppression. Azathioprine: "
                f"{_ALTERNATIVE} Mercaptopurine / thioguanine: {_MP_PM} "
                f"No single percent-of-standard applies across the three "
                f"drugs, so none is reported."
            )
        elif tpmt_class == "Intermediate" and nudt15_class == "Intermediate":
            dose_range = (20, 50)
            dose_percent = 20
            shared = (
                "Initiate therapy with reduced starting doses (20-50% of "
                "standard starting dose) if the standard starting dose is "
                "at or above 2 mg/kg/day (azathioprine) or 75 mg/m2/day or "
                "1.5 mg/kg/day (mercaptopurine). If the starting dose is "
                "already below standard, dose reduction might not be "
                "necessary. Allow 2-4 weeks of stable dosing to reach "
                "steady state after each adjustment."
            )
            recommendations_by_drug = {d: shared for d in affected_drugs}
            dose_recommendation = (
                "TPMT intermediate AND NUDT15 intermediate — reduced "
                "clearance at both loci. " + shared
            )
        elif "Intermediate" in classes:
            which = "TPMT" if tpmt_class == "Intermediate" else "NUDT15"
            dose_range = (30, 80)
            dose_percent = 30
            shared = (
                "Initiate therapy with reduced starting doses (30-80% of "
                "standard starting dose) if the standard starting dose is "
                "at or above 2 mg/kg/day (azathioprine) or 75 mg/m2/day or "
                "1.5 mg/kg/day (mercaptopurine). If the starting dose is "
                "already below standard, dose reduction might not be "
                "necessary. Allow 2-4 weeks of stable dosing to reach "
                "steady state after each adjustment."
            )
            recommendations_by_drug = {d: shared for d in affected_drugs}
            dose_recommendation = (
                f"{which} intermediate metabolizer. " + shared
            )
        else:
            dose_range = (100, 100)
            dose_percent = 100
            shared = (
                "Initiate therapy with the standard starting dose (e.g. 2 "
                "mg/kg/day azathioprine for autoimmune disease; 75 mg/m2/day "
                "mercaptopurine for malignancy). Adjust per disease-specific "
                "guidelines. Even with a normal genotype, idiosyncratic "
                "toxicity is possible — maintain standard monitoring."
            )
            recommendations_by_drug = {d: shared for d in affected_drugs}
            dose_recommendation = (
                "Standard thiopurine dose. Normal TPMT and NUDT15 activity. "
                + shared
            )

        result = {
            "dose_recommendation": dose_recommendation,
            # CPIC gives a RANGE, not a point. dose_percent_of_standard is
            # the conservative (low) end of that range, which is where a
            # starting dose belongs; dose_percent_range_of_standard carries
            # CPIC's full range. None means CPIC recommends a different drug
            # rather than a reduced dose.
            "dose_percent_of_standard": dose_percent,
            "dose_percent_range_of_standard": dose_range,
            "recommendations_by_drug": recommendations_by_drug,
            "tpmt_diplotype": tpmt_diplotype,
            "tpmt_activity_score": tpmt_score,
            "tpmt_phenotype": tpmt_phenotype,
            "nudt15_diplotype": nudt15_diplotype,
            "nudt15_activity_score": nudt15_score,
            "nudt15_phenotype": nudt15_phenotype,
            "affected_drugs": affected_drugs,
            "guideline": (
                "CPIC Thiopurine Dosing Based on TPMT and NUDT15 Genotypes, "
                "2018 Update (Relling et al., Clin Pharmacol Ther "
                "2019;105:1095-1105); dose rows from CPIC's maintained "
                "azathioprine and mercaptopurine recommendation tables, "
                "api.cpicpgx.org/v1/recommendation, retrieved 2026-11-03"
            ),
            "monitoring": (
                "CBC with differential: weekly x8 weeks, then biweekly x4 "
                "weeks, then at least every 3 months. More frequent "
                "monitoring for intermediate/deficient patients. Monitor "
                "thiopurine metabolites (6-TGN, 6-MMP) to guide dosing."
            ),
        }

        logger.info(
            "Thiopurine dosing: TPMT %s (AS %.1f, %s), NUDT15 %s (AS %.1f, %s) → %s",
            tpmt_diplotype, tpmt_score, tpmt_class,
            nudt15_diplotype, nudt15_score, nudt15_class,
            f"{dose_range[0]}-{dose_range[1]}% of standard" if dose_range
            else "alternative drug recommended",
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
        try:
            allele1, allele2, score1, score2 = self._score_diplotype(
                "CYP2C19", cyp2c19_diplotype, _CYP2C19_ACTIVITY,
                _CYP2C19_UNCERTAIN)
        except IndeterminateGenotype as exc:
            return {
                "cyp2c19_diplotype": cyp2c19_diplotype,
                "activity_score": None,
                "phenotype": "Indeterminate",
                "recommended_therapy": None,
                "clopidogrel_use": "indeterminate",
                "rationale": (
                    f"No CYP2C19-guided antiplatelet recommendation can be "
                    f"given: {exc.reason} Do NOT read this as a normal "
                    f"metabolizer — an unrecognised CYP2C19 allele is as "
                    f"likely to be a no-function one, and CPIC's Intermediate "
                    f"and Poor Metabolizer rows both say to avoid standard "
                    f"clopidogrel."
                ),
                "warnings": [
                    "CYP2C19 genotype could not be scored; clopidogrel "
                    "response cannot be predicted from it",
                ],
                "indeterminate_alleles": exc.alleles,
                "affected_drugs": ["clopidogrel"],
                "reference": (
                    "CPIC Guideline for Clopidogrel and CYP2C19, 2022 Update; "
                    "Indeterminate row of api.cpicpgx.org/v1/recommendation"
                ),
            }
        activity_score = score1 + score2

        # Phenotype boundaries follow CPIC's diplotype table: the Poor
        # bucket's top is activity score 1.0 (exclusive), not 0.0 — CPIC
        # calls a no-function + decreased-function pair (AS 0.5, e.g.
        # *2/*10) a Poor Metabolizer.
        if activity_score < 1.0:
            phenotype = "Poor Metabolizer"
            recommended_therapy = "prasugrel 10mg daily OR ticagrelor 90mg BID"
            clopidogrel_use = "contraindicated"
            rationale = (
                "Avoid clopidogrel if possible. Use prasugrel or ticagrelor "
                "at standard dose if no contraindication. Clopidogrel is "
                "INEFFECTIVE in CYP2C19 poor metabolizers — no conversion "
                "to active metabolite. CPIC Level A."
            )
            warnings = [
                "Clopidogrel provides no antiplatelet benefit in PMs",
                "Elevated risk of MACE, stent thrombosis, and recurrent ACS",
            ]
        elif activity_score < 2.0:
            phenotype = "Intermediate Metabolizer"
            recommended_therapy = "prasugrel 10mg daily OR ticagrelor 90mg BID"
            clopidogrel_use = "avoid_standard_dose"
            # CPIC verbatim. This branch used to read "Consider alternative
            # antiplatelet therapy ... clopidogrel 75mg with enhanced
            # monitoring (if alternatives unavailable)" with
            # clopidogrel_use "use_with_caution" — softer than the guideline
            # it cites.
            rationale = (
                "Avoid standard dose (75 mg) clopidogrel if possible. Use "
                "prasugrel or ticagrelor at standard dose if no "
                "contraindication. CYP2C19 intermediate metabolizers have "
                "significantly reduced clopidogrel activation."
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

        CPIC's own simvastatin recommendation rows, verbatim
        (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:36567, retrieved
        2026-11-03):

            Increased / Normal Function
                "Prescribe desired starting dose and adjust doses based on
                 disease-specific guidelines."
            Decreased and Possible Decreased Function
                "Prescribe an alternative statin depending on the desired
                 potency ... If simvastatin therapy is warranted, limit dose
                 to <20mg/day."
            Poor Function
                "Prescribe an alternative statin depending on the desired
                 potency."
            Indeterminate
                "No recommendation."

        Note CPIC's bound is **strictly less than** 20 mg/day; 20 mg is the
        one dose it excludes.

        Parameters
        ----------
        slco1b1_diplotype : str
            e.g. "*1/*1", "*1/*5", "*5/*5", "*1/*15"

        Returns
        -------
        dict
            slco1b1_diplotype, phenotype, max_simvastatin_dose_mg,
            max_simvastatin_dose_bound, dose_adjustment, myopathy_risk,
            rationale, warnings, alternatives, affected_drugs, reference.
        """
        alternatives = [
            "rosuvastatin", "pravastatin", "fluvastatin", "pitavastatin",
        ]
        base = {
            "slco1b1_diplotype": slco1b1_diplotype,
            "alternatives": alternatives,
            "affected_drugs": ["simvastatin", "lovastatin"],
            "reference": (
                "CPIC Guideline for Simvastatin and SLCO1B1, 2022 Update "
                "(Cooper-DeHoff et al., Clin Pharmacol Ther 2022;112:1125); "
                "recommendation rows from api.cpicpgx.org/v1/recommendation, "
                "retrieved 2026-11-03. Labelled dose range from the "
                "FDA-approved ZOCOR label (DailyMed, Organon LLC, setid "
                "8f55d5de-5a4f-4a39-8c84-c53976dd6af9) §2.1-2.2."
            ),
        }

        # An unrecognised diplotype used to be reported as phenotype
        # "Unknown" with max_simvastatin_dose_mg 40 and an EMPTY warnings
        # list, so simvastatin_dose("5/5") -- a lab reporting without
        # asterisks for a genuinely *5/*5 patient, CPIC Poor Function with
        # ~17-fold myopathy risk -- handed back a 40 mg ceiling. CPIC's
        # answer for Indeterminate is "No recommendation", which is not 40
        # mg.
        try:
            allele1, allele2 = self._parse_diplotype(slco1b1_diplotype)
        except IndeterminateGenotype as exc:
            allele1 = allele2 = ""
            unrecognised = [exc.genotype]
        else:
            unrecognised = [
                a for a in (allele1, allele2)
                if a not in _SLCO1B1_ALLELE_FUNCTION
                and a not in _SLCO1B1_UNCERTAIN
            ]
        if unrecognised:
            return {
                **base,
                "phenotype": "Indeterminate",
                "max_simvastatin_dose_mg": None,
                "max_simvastatin_dose_bound": None,
                "dose_adjustment": "no_recommendation",
                "myopathy_risk": "unknown — genotype not recognised",
                "rationale": (
                    f"SLCO1B1 "
                    f"{', '.join(dict.fromkeys(unrecognised))} is not in CPIC's "
                    f"SLCO1B1 allele table, so no phenotype and no dose "
                    f"ceiling can be assigned. CPIC's answer for Indeterminate "
                    f"is 'No recommendation'. Do NOT read this as normal "
                    f"transporter function: an unrecognised label may be a "
                    f"no-function allele, for which CPIC says prescribe an "
                    f"alternative statin."
                ),
                "warnings": [
                    "SLCO1B1 genotype not recognised — myopathy risk cannot "
                    "be assessed from it",
                    "Check the reported genotype against CPIC's SLCO1B1 "
                    "allele-functionality table before prescribing",
                ],
                "indeterminate_alleles": list(dict.fromkeys(unrecognised)),
            }

        def _fn(allele: str) -> str:
            return _SLCO1B1_ALLELE_FUNCTION.get(allele, "uncertain")

        phenotype = _SLCO1B1_PAIR_PHENOTYPE.get(
            frozenset({_fn(allele1), _fn(allele2)}), "Indeterminate")
        myopathy_risk = _SLCO1B1_MYOPATHY_RISK.get(phenotype, "unknown")

        if phenotype in ("Normal Function", "Increased Function"):
            # WAS max_dose = 80. The FDA-approved ZOCOR label: §2.2 "The
            # recommended dosage range of ZOCOR is 20 mg to 40 mg once
            # daily"; §2.1 and §5.1 "An 80 mg daily dosage of ZOCOR is
            # restricted to patients who have been taking simvastatin 80 mg
            # daily chronically (e.g., for 12 months or more) without
            # evidence of muscle toxicity"; §2.1 "The 5 mg and 80 mg
            # strengths of ZOCOR are no longer marketed". CPIC names no
            # number at all. So 80 was double the top of the labelled range,
            # restricted to patients already established on it, and no
            # longer manufactured -- and SLCO1B1 (OATP1B1, named in §5.1) is
            # the very risk factor that restriction exists for.
            max_dose: Any = 40
            bound = "<=40 mg/day (top of the labelled range)"
            dose_adjustment = "none"
            rationale = (
                "Prescribe the desired starting dose and adjust based on "
                "disease-specific guidelines. No increased myopathy risk "
                "from SLCO1B1 genotype. The labelled ZOCOR dosage range is "
                "20-40 mg once daily; 80 mg is restricted to patients "
                "already established on it for 12 months or more without "
                "muscle toxicity, and that strength is no longer marketed."
            )
            warnings: List[str] = []
        elif phenotype in ("Decreased Function", "Possible Decreased Function"):
            # WAS max_dose = 20 with the text "Prescribe simvastatin
            # <=20mg/day". CPIC says **<20 mg/day** -- strictly less -- so
            # 20 mg is the one dose CPIC excludes, and it was the number in
            # the structured field a prescriber reads.
            max_dose = 10
            bound = "<20 mg/day (CPIC); highest marketed strength below it is 10 mg"
            dose_adjustment = "alternative_statin_or_below_20mg"
            rationale = (
                "Prescribe an alternative statin depending on the desired "
                "potency. If simvastatin therapy is warranted, limit the dose "
                "to LESS THAN 20 mg/day — CPIC's bound is strictly below 20 "
                "mg, so 20 mg itself is excluded; the highest marketed "
                "strength below it is 10 mg. SLCO1B1 "
                + ("decreased" if phenotype == "Decreased Function"
                   else "possible decreased")
                + " function increases simvastatin exposure and myopathy risk."
            )
            warnings = [
                "Do NOT use simvastatin 20mg or above",
                "Monitor for muscle pain, weakness, elevated CK",
            ]
            if phenotype == "Possible Decreased Function":
                warnings.append(
                    "One allele is of uncertain function; CPIC gives this "
                    "class the same recommendation as Decreased Function"
                )
        elif phenotype == "Poor Function":
            max_dose = 0
            bound = "no simvastatin dose recommended"
            dose_adjustment = "avoid_simvastatin"
            rationale = (
                "Prescribe an alternative statin depending on the desired "
                "potency. SLCO1B1 poor function confers very high "
                "myopathy/rhabdomyolysis risk. Use rosuvastatin, "
                "pravastatin, or fluvastatin instead."
            )
            warnings = [
                "Simvastatin contraindicated",
                "Also applies to lovastatin (similar transport)",
                "Rosuvastatin, pravastatin, fluvastatin are safer alternatives",
            ]
        else:
            max_dose = None
            bound = None
            dose_adjustment = "no_recommendation"
            rationale = (
                f"SLCO1B1 {slco1b1_diplotype}: at least one allele is of "
                f"uncertain function at CPIC and the pair maps to no "
                f"phenotype. CPIC's answer for Indeterminate is 'No "
                f"recommendation'. Do NOT read this as normal transporter "
                f"function."
            )
            warnings = [
                "SLCO1B1 phenotype Indeterminate — myopathy risk cannot be "
                "assessed from this genotype",
            ]

        result = {
            **base,
            "phenotype": phenotype,
            "max_simvastatin_dose_mg": max_dose,
            "max_simvastatin_dose_bound": bound,
            "dose_adjustment": dose_adjustment,
            "myopathy_risk": myopathy_risk,
            "rationale": rationale,
            "warnings": warnings,
        }

        logger.info(
            "Simvastatin dosing: %s (SLCO1B1 %s, ceiling %s)",
            phenotype, slco1b1_diplotype, bound,
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
        try:
            allele1, allele2, score1, score2 = self._score_diplotype(
                "CYP2C19", cyp2c19_diplotype, _CYP2C19_ACTIVITY,
                _CYP2C19_UNCERTAIN)
        except IndeterminateGenotype as exc:
            return {
                "drug": drug,
                "cyp2c19_diplotype": cyp2c19_diplotype,
                "activity_score": None,
                "phenotype": "Indeterminate",
                "max_dose_mg": None,
                "dose_adjustment": "no_recommendation",
                "rationale": (
                    f"No CYP2C19-guided {drug_lower} dose can be given: "
                    f"{exc.reason} CPIC's answer for Indeterminate is 'No "
                    f"recommendation'. Prescribe within the drug label and do "
                    f"NOT read this as a normal-metabolizer result."
                ),
                "warnings": [
                    "CYP2C19 genotype could not be scored",
                ],
                "indeterminate_alleles": exc.alleles,
                "alternatives": ["sertraline", "mirtazapine", "bupropion"],
                "reference": "CPIC Guideline for SSRIs and CYP2C19, 2023 Update",
            }
        activity_score = score1 + score2

        # Determine phenotype (same thresholds as clopidogrel_dose)
        if activity_score < 1.0:
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
                # WAS max_dose_mg = 60 with "Consider higher dose".
                # FDA-approved CELEXA label (DailyMed, Allergan, setid
                # 4259d9b1-...): §2.1 "Dosages above 40 mg once daily are
                # not recommended due to the risk of QT prolongation";
                # §12.2 the thorough QT study measured 60 mg, described in
                # the label as "1.5 times the maximum recommended dosage",
                # and found a QTcNi increase of 18.5 msec (upper bound
                # 21.0). So the agent was emitting, as a structured field,
                # the precise dose the label reserves for its own QT safety
                # study. CPIC says nothing about 60 mg either; its UM row is
                # "Consider a clinically appropriate alternative
                # antidepressant not predominantly metabolized by CYP2C19.
                # If citalopram or escitalopram are clinically appropriate,
                # and adequate efficacy is not achieved at standard
                # maintenance dosing, consider titrating to a higher
                # maintenance dose" -- within the label, which caps at 40 mg
                # (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:2556,
                # retrieved 2026-11-03).
                max_dose_mg = 40
                dose_adjustment = "consider_alternative_or_titrate_within_label"
                rationale = (
                    "Consider a clinically appropriate alternative "
                    "antidepressant not predominantly metabolized by CYP2C19. "
                    "If citalopram is clinically appropriate and adequate "
                    "efficacy is not achieved at standard maintenance dosing, "
                    "consider titrating to a higher maintenance dose WITHIN "
                    "the label: the CELEXA label does not recommend dosages "
                    "above 40 mg once daily because of QT prolongation."
                )
                warnings = [
                    "Do NOT exceed 40 mg/day — the label's own thorough QT "
                    "study measured +18.5 msec QTcNi at 60 mg",
                    "Monitor for therapeutic response",
                    "Consider an alternative not dependent on CYP2C19",
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
                # WAS max_dose_mg = 30. FDA-approved LEXAPRO label
                # (DailyMed, Allergan, setid
                # 13bb8267-1cab-43e5-acae-55a4d957630a): §2.1 "The
                # recommended dosage of Lexapro in adults is 10 mg once
                # daily", which "may be increased to the maximum
                # recommended dosage of 20 mg once daily". The label names
                # NO recommended dose above 20 mg once daily for adult MDD;
                # it names 30 mg exactly once, in §6.1, as a
                # "supratherapeutic 30 mg" once-daily arm OF A QTcF STUDY.
                # Same shape as the citalopram 60 mg: the number was taken
                # from the pharmacology section, not the dosing section.
                max_dose_mg = 20
                dose_adjustment = "consider_alternative_or_titrate_within_label"
                rationale = (
                    "Consider a clinically appropriate alternative "
                    "antidepressant not predominantly metabolized by CYP2C19. "
                    "If escitalopram is clinically appropriate and adequate "
                    "efficacy is not achieved at standard maintenance dosing, "
                    "consider titrating to a higher maintenance dose WITHIN "
                    "the label: the LEXAPRO label's maximum recommended "
                    "adult dosage is 20 mg once daily."
                )
                warnings = [
                    "Do NOT exceed 20 mg/day — the label names 30 mg only as "
                    "a supratherapeutic arm of a QTcF study",
                    "Monitor for therapeutic response",
                    "Consider an alternative not dependent on CYP2C19",
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
            "reference": (
                "CPIC Guideline for SSRIs and CYP2C19, 2023 Update (Hicks et "
                "al.); recommendation rows from "
                "api.cpicpgx.org/v1/recommendation, retrieved 2026-11-03. "
                "Dose ceilings from the FDA-approved CELEXA label (DailyMed, "
                "Allergan) §2.1/§2.3/§12.2 and LEXAPRO label (DailyMed, "
                "Allergan, setid 13bb8267-1cab-43e5-acae-55a4d957630a) §2.1."
            ),
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

        CPIC's own phenytoin recommendation rows, by CYP2C9 activity score
        (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:8183, retrieved
        2026-11-03):

            2.0  "No adjustments needed from typical dosing strategies."
            1.5  "No adjustments needed from typical dosing strategies.
                  Subsequent doses should be adjusted according to
                  therapeutic drug monitoring, response and side effects."
            1.0  "For first dose, use typical initial or loading dose. For
                  subsequent doses, use approximately 25% less than typical
                  maintenance dose."
            0.5  as 1.0 but "approximately 50% less".
            0.0  as 0.5.
            n/a  (Indeterminate) "No recommendation."

        Every row also carries: "An HLA-B*15:02 negative test does not
        eliminate the risk of phenytoin-induced SJS/TEN", and CPIC's first
        phenytoin recommendation for an HLA-B*15:02-POSITIVE patient is "If
        patient is phenytoin-naive, do not use phenytoin/fosphenytoin".

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
            hla_b_1502_caveat, tdm_target, alternatives, affected_drugs,
            reference.
        """
        # CPIC attaches this to EVERY phenytoin row, and a function that
        # returned "Standard phenytoin dosing" under a heading citing the
        # phenytoin guideline while omitting the guideline's own headline
        # recommendation misrepresented it by omission.
        hla_caveat = (
            "HLA-B*15:02: if the patient is HLA-B*15:02 POSITIVE and "
            "phenytoin-naive, CPIC says do not use phenytoin or "
            "fosphenytoin — risk of Stevens-Johnson syndrome / toxic "
            "epidermal necrolysis. An HLA-B*15:02 NEGATIVE test does not "
            "eliminate that risk. This function does not know the patient's "
            "HLA-B status; screen separately (see hla_screener.py)."
        )

        try:
            allele1, allele2, score1, score2 = self._score_diplotype(
                "CYP2C9", cyp2c9_diplotype, _CYP2C9_ACTIVITY,
                _CYP2C9_UNCERTAIN)
        except IndeterminateGenotype as exc:
            return {
                "cyp2c9_diplotype": cyp2c9_diplotype,
                "activity_score": None,
                "phenotype": "Indeterminate",
                "dose_reduction_percent": None,
                "standard_loading_dose_mg": None,
                "suggested_maintenance_mg_per_day": None,
                "max_maintenance_mg_per_day": None,
                "rationale": (
                    f"No CYP2C9-guided phenytoin dose can be given: "
                    f"{exc.reason} CPIC's answer for Indeterminate is 'No "
                    f"recommendation'. Phenytoin has saturable kinetics and a "
                    f"narrow window; dose by therapeutic drug monitoring, and "
                    f"do NOT read this as a normal-metabolizer result."
                ),
                "warnings": [
                    "CYP2C9 genotype could not be scored",
                    "Therapeutic drug monitoring required",
                ],
                "indeterminate_alleles": exc.alleles,
                "hla_b_1502_caveat": hla_caveat,
                "tdm_target": "10-20 mcg/mL",
                "alternatives": ["levetiracetam", "lacosamide", "valproic acid"],
                "affected_drugs": ["phenytoin", "fosphenytoin"],
                "reference": (
                    "CPIC Guideline for Phenytoin and CYP2C9, 2020 Update"
                ),
            }
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
        elif activity_score < 1.5:
            phenotype = "Intermediate Metabolizer"
            dose_reduction_percent = 25
            rationale = (
                "For the first dose, use the typical initial or loading dose. "
                "For subsequent doses, use approximately 25% less than the "
                "typical maintenance dose, then adjust according to "
                "therapeutic drug monitoring, response and side effects."
            )
            warnings = [
                "Reduce maintenance dose by 25%",
                "TDM recommended",
                "Monitor for dose-related side effects",
            ]
        else:
            # Activity score 1.5 belongs HERE, not in the 25%-reduction
            # branch. The condition used to be `activity_score < 2.0`, which
            # lumped every score in [1.0, 2.0) into one 25% bucket. CPIC
            # splits it and gives AS 1.5 "No adjustments needed from typical
            # dosing strategies" -- and AS 1.5 is CYP2C9 *1/*2, the commonest
            # actionable CYP2C9 genotype in European-ancestry patients
            # (10-15% allele frequency). Starting such a patient 25% below
            # the standard maintenance dose of an antiepileptic with
            # saturable kinetics risks a subtherapeutic level and a
            # breakthrough seizure; under-dosing an anticonvulsant is not
            # the safe direction of error.
            phenotype = (
                "Intermediate Metabolizer" if activity_score < 2.0
                else "Normal Metabolizer"
            )
            dose_reduction_percent = 0
            rationale = (
                f"No adjustments needed from typical dosing strategies "
                f"({phenotype}, activity score {activity_score}). Subsequent "
                f"doses should be adjusted according to therapeutic drug "
                f"monitoring, response and side effects."
            )
            warnings = (
                ["TDM recommended — one decreased-function CYP2C9 allele, but "
                 "CPIC makes no starting-dose adjustment at activity score 1.5"]
                if activity_score < 2.0 else []
            )

        # Calculate doses
        standard_loading_dose_mg = min(round(weight_kg * 18), 1500)
        suggested_maintenance_mg = round(
            weight_kg * 5 * (1 - dose_reduction_percent / 100)
        )

        # Maintenance dose range string
        if phenotype == "Poor Metabolizer":
            max_maintenance_range = "150-200"
        elif dose_reduction_percent == 25:
            max_maintenance_range = "225-300"
        else:
            max_maintenance_range = "300-400"

        warnings = warnings + [hla_caveat]

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
            "hla_b_1502_caveat": hla_caveat,
            "tdm_target": "10-20 mcg/mL",
            "alternatives": ["levetiracetam", "lacosamide", "valproic acid"],
            "affected_drugs": ["phenytoin", "fosphenytoin"],
            "reference": (
                "CPIC Guideline for Phenytoin and CYP2C9, 2020 Update "
                "(Karnes et al., Clin Pharmacol Ther 2021;109:302); "
                "recommendation rows from api.cpicpgx.org/v1/recommendation, "
                "retrieved 2026-11-03"
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

        _indeterminate = {
            "drug": drug,
            "tca_class": tca_class,
            "cyp2d6_diplotype": cyp2d6_diplotype,
            "cyp2c19_diplotype": cyp2c19_diplotype,
            "dose_adjustment": "no_recommendation",
            "tdm_required": True,
            "tdm_target_range": _TCA_TDM_RANGES.get(
                drug_lower, "See drug-specific reference"),
            "alternatives": ["sertraline", "venlafaxine", "mirtazapine"],
            "affected_drugs": list(_TCA_CLASS.keys()),
            "reference": (
                "CPIC Guideline for Tricyclic Antidepressants and "
                "CYP2D6/CYP2C19, 2016 Update"
            ),
        }

        # ── CYP2D6 activity score ──
        try:
            d6_a1, d6_a2, d6_score1, d6_score2 = self._score_diplotype(
                "CYP2D6", cyp2d6_diplotype, _CYP2D6_ACTIVITY)
        except IndeterminateGenotype as exc:
            return {
                **_indeterminate,
                "cyp2d6_activity_score": None,
                "cyp2d6_phenotype": "Indeterminate",
                "cyp2c19_phenotype": "Not evaluated",
                "rationale": (
                    f"No genotype-guided {drug_lower} dose can be given: "
                    f"{exc.reason} TCAs have a narrow therapeutic index; dose "
                    f"by therapeutic drug monitoring and do NOT read this as a "
                    f"normal-metabolizer result."
                ),
                "warnings": ["CYP2D6 genotype could not be scored"],
                "indeterminate_alleles": exc.alleles,
            }
        cyp2d6_as = d6_score1 + d6_score2

        # CYP2D6 phenotype — CPIC/DPWG consensus: PM = 0, IM = 0 < AS < 1.25,
        # NM = 1.25-2.25 inclusive, UM > 2.25. The consensus reclassified AS 1.0 from NM to IM,
        # and 2.25 is the inclusive TOP of NM, not the bottom of UM.
        #
        # This block previously cut IM at < 1.0 and NM at < 2.25, so a patient with AS 1.0
        # (*1/*4) got standard tricyclic dosing where CPIC calls for a 25% reduction, and a
        # patient with AS 2.25 got the ultrarapid recommendation to switch drug class.
        #
        # Source: Caudle KE et al. "Standardizing CYP2D6 Genotype to Phenotype Translation:
        # Consensus Recommendations from the Clinical Pharmacogenetics Implementation
        # Consortium and Dutch Pharmacogenetics Working Group." Clin Transl Sci
        # 2020;13(1):116-124, Table 3. doi:10.1111/cts.12692, PMID 31647186.
        if cyp2d6_as == 0.0:
            cyp2d6_phenotype = "Poor Metabolizer"
        elif cyp2d6_as < 1.25:
            cyp2d6_phenotype = "Intermediate Metabolizer"
        elif cyp2d6_as <= 2.25:
            cyp2d6_phenotype = "Normal Metabolizer"
        else:
            cyp2d6_phenotype = "Ultrarapid Metabolizer"

        # ── CYP2C19 activity score (second axis of CPIC's table) ──
        try:
            c19_a1, c19_a2, c19_score1, c19_score2 = self._score_diplotype(
                "CYP2C19", cyp2c19_diplotype, _CYP2C19_ACTIVITY,
                _CYP2C19_UNCERTAIN)
        except IndeterminateGenotype as exc:
            return {
                **_indeterminate,
                "cyp2d6_activity_score": cyp2d6_as,
                "cyp2d6_phenotype": cyp2d6_phenotype,
                "cyp2c19_phenotype": "Indeterminate",
                "rationale": (
                    f"No genotype-guided {drug_lower} dose can be given: "
                    f"{exc.reason} CPIC's TCA recommendation is a "
                    f"two-dimensional table over CYP2D6 x CYP2C19 and the "
                    f"CYP2C19 axis can make the answer 'avoid', so a CYP2D6 "
                    f"result alone is not enough. Dose by therapeutic drug "
                    f"monitoring."
                ),
                "warnings": ["CYP2C19 genotype could not be scored"],
                "indeterminate_alleles": exc.alleles,
            }
        cyp2c19_as = c19_score1 + c19_score2

        if cyp2c19_as < 1.0:
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
            # CPIC: "Avoid amitriptyline use; if amitriptyline is warranted,
            # utilize therapeutic drug monitoring to guide dose adjustment"
            # for a tertiary amine, and for a secondary amine "avoid ...
            # consider titrating to a higher target dose". The old secondary-
            # amine branch said "may attempt with mandatory TDM", which is
            # softer than the guideline's "avoid".
            dose_adjustment = "avoid"
            rationale = (
                "Avoid " + drug_lower + " use. CYP2D6 ultrarapid metabolism "
                "gives subtherapeutic parent and metabolite levels. If a TCA "
                "is warranted, utilize therapeutic drug monitoring to guide "
                "dose adjustment"
                + ("; a secondary amine may be titrated to a higher target "
                   "dose than in normal metabolizers." if not is_tertiary
                   else ".")
            )
            tdm_required = True
            warnings = [
                "Avoid this TCA — treatment failure expected",
                "If used anyway, mandatory TDM",
                "Consider SSRI/SNRI alternative",
            ]
        elif cyp2d6_phenotype == "Poor Metabolizer":
            # CPIC leads with "Avoid amitriptyline use." before the 50%.
            dose_adjustment = "avoid_or_reduce_50_percent"
            rationale = (
                "Avoid " + drug_lower + " use. If " + drug_lower + " is "
                "warranted, consider a 50% reduction of the recommended "
                "starting dose; utilize therapeutic drug monitoring to guide "
                "dose adjustments. CYP2D6 PM causes a 2-4 fold increase in "
                "TCA exposure."
            )
            tdm_required = True
            warnings = [
                "CPIC says avoid; if used, reduce dose by 50%",
                "Mandatory TDM — target therapeutic range",
                "Risk of cardiac toxicity (QT prolongation), "
                "anticholinergic crisis, seizures",
                "Consider SSRI alternative if TDM unavailable",
            ]
        elif cyp2d6_phenotype == "Intermediate Metabolizer":
            dose_adjustment = "reduce_25_percent"
            rationale = (
                "Consider a 25% reduction of the recommended starting dose. "
                "Utilize therapeutic drug monitoring to guide dose "
                "adjustments. CYP2D6 IM has moderately reduced TCA clearance."
            )
            tdm_required = True
            warnings = [
                "Consider 25% dose reduction",
                "TDM recommended",
            ]
        else:  # Normal Metabolizer
            dose_adjustment = "none"
            rationale = (
                "Initiate therapy with the recommended starting dose. Normal "
                "CYP2D6 metabolism."
            )
            tdm_required = False

        # ── CYP2C19 axis: the SECOND dimension of CPIC's table ──
        # CPIC's TCA guidance is a two-dimensional table over CYP2D6 x
        # CYP2C19. This used to treat CYP2D6 as the decision and CYP2C19 as
        # a footnote: dose_adjustment, tdm_required and rationale were set
        # from CYP2D6 alone and the CYP2C19 contribution was an appended
        # free-text string, so tca_dose('amitriptyline', '*1/*1', '*2/*2')
        # returned dose_adjustment "none", tdm_required False and "Standard
        # TCA dosing. Normal CYP2D6 metabolism." for a CYP2C19 poor
        # metabolizer — 12-23% of East Asian populations. CPIC's row for
        # every CYP2D6-Normal x CYP2C19-Poor pair (CYP2D6 AS 1.25, 1.5,
        # 1.75, 2.0 and 2.25 alike) reads: "Avoid amitriptyline use. If
        # amitriptyline is warranted, consider a 50% reduction of
        # recommended starting dose. Utilize therapeutic drug monitoring to
        # guide dose adjustments."
        # (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:704, retrieved
        # 2026-11-03.)
        #
        # The structured fields now take the MORE RESTRICTIVE of the two
        # genes' recommendations.
        _RESTRICTIVENESS = {
            "none": 0,
            "reduce_25_percent": 1,
            "avoid_or_reduce_50_percent": 2,
            "avoid": 3,
        }
        if is_tertiary:
            c19_adjustment = None
            if cyp2c19_phenotype == "Poor Metabolizer":
                c19_adjustment = "avoid_or_reduce_50_percent"
                c19_rationale = (
                    "CYP2C19 poor metabolism raises tertiary-amine "
                    "parent-drug exposure. CPIC: avoid " + drug_lower + " "
                    "use; if " + drug_lower + " is warranted, consider a 50% "
                    "reduction of the recommended starting dose, with "
                    "therapeutic drug monitoring."
                )
            elif cyp2c19_phenotype in (
                    "Rapid Metabolizer", "Ultrarapid Metabolizer"):
                c19_adjustment = "avoid"
                c19_rationale = (
                    "CPIC: consider an alternative drug not metabolized by "
                    "CYP2C19; if " + drug_lower + " is warranted, utilize "
                    "therapeutic drug monitoring to guide dose adjustment."
                )
            elif cyp2c19_phenotype == "Intermediate Metabolizer":
                c19_adjustment = None
                c19_rationale = ""
                warnings.append(
                    "CYP2C19 intermediate metabolizer; CPIC makes no "
                    "additional starting-dose adjustment on this axis"
                )

            if c19_adjustment and (
                _RESTRICTIVENESS.get(c19_adjustment, 0)
                > _RESTRICTIVENESS.get(dose_adjustment, 0)
            ):
                dose_adjustment = c19_adjustment
                rationale = c19_rationale
                tdm_required = True
                warnings.insert(0, (
                    f"CYP2C19 {cyp2c19_phenotype} drives this recommendation, "
                    f"not CYP2D6"
                ))
            elif c19_adjustment:
                tdm_required = True
                warnings.append(c19_rationale)

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
                "CYP2D6/CYP2C19, 2016 Update (Hicks et al., Clin Pharmacol "
                "Ther 2017;102:37); two-gene recommendation rows from "
                "api.cpicpgx.org/v1/recommendation, retrieved 2026-11-03"
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
        """Split a diplotype string into its two alleles.

        Raises
        ------
        IndeterminateGenotype
            If the string is not exactly two slash-separated non-empty
            fields.  This used to return ``("*1", "*1")`` — homozygous
            wild-type — for an empty string, a single allele label, a
            three-field string, "Unknown", "not reported" or any free text,
            and every one of the nine calculators downstream then reported a
            normal metabolizer and a standard dose with no warning anywhere
            in the result.  ``fluoropyrimidine_dose("*2A")`` — a pharmacy
            reporting a patient as a *2A carrier — printed "Standard
            fluoropyrimidine dose. No DPYD-based dose adjustment needed.
            Activity score 2.0 (normal)", which is the published fatal-5-FU-
            toxicity scenario.  An empty field is the other common case: a
            genotype that was ordered and has not come back yet.
        """
        parts = diplotype.split("/") if diplotype else []
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            return parts[0].strip(), parts[1].strip()
        raise IndeterminateGenotype(
            genotype=diplotype,
            reason=(
                f"{diplotype!r} is not a diplotype — two slash-separated "
                f"alleles were expected, e.g. '*1/*2'. This is NOT wild-type; "
                f"a single allele label or an empty result cannot be scored."
            ),
            alleles=[str(diplotype)],
        )

    @classmethod
    def _score_diplotype(
        cls,
        gene: str,
        diplotype: str,
        table: Dict[str, float],
        uncertain: frozenset = frozenset(),
    ) -> Tuple[str, str, float, float]:
        """Parse a diplotype and look up both CPIC activity values, or refuse.

        Replaces ``_ACTIVITY.get(allele, 1.0)`` — "an allele I have never
        seen is fully functional", the least safe default available.  The
        repo scores a fraction of what CPIC assigns (11 of 39 CYP2C19
        alleles, 7 of ~50 CYP2C9, 21 of >100 CYP2D6, 9 of 84 DPYD), so that
        default was reached often and it reported no-function alleles as
        normal metabolizers: ``clopidogrel_dose("*1/*35")`` returned "Normal
        Metabolizer, standard clopidogrel" where CPIC's diplotype table says
        Intermediate Metabolizer and "avoid standard dose clopidogrel (75
        mg) ... use prasugrel or ticagrelor".  CPIC's own answer for an
        allele it cannot call is **Indeterminate** / "No recommendation".

        Raises
        ------
        IndeterminateGenotype
            If the string is unparseable, or either allele is uncertain-
            function or absent from CPIC's assigned values.
        """
        allele1, allele2 = cls._parse_diplotype(diplotype)
        unscoreable = [
            (a, "uncertain function at CPIC" if a in uncertain
                else "not in this agent's %s allele table" % gene)
            for a in (allele1, allele2)
            if a in uncertain or a not in table
        ]
        if unscoreable:
            raise IndeterminateGenotype(
                genotype=diplotype,
                reason=(
                    f"{gene} "
                    + "; ".join(f"{a} ({why})" for a, why in unscoreable)
                    + ". CPIC assigns Indeterminate and makes no "
                      "recommendation here. Do NOT infer normal function "
                      "from the absence of a call."
                ),
                alleles=[a for a, _ in unscoreable],
            )
        score1, score2 = table[allele1], table[allele2]
        # An increased-function allele does not rescue a reduced-function
        # one — see PhenotypeTranslator.INCREASED_FUNCTION_SCORE for CPIC's
        # own diplotype rows.
        if score1 == 1.5 and score2 < 1.0:
            score1 = 1.0
        elif score2 == 1.5 and score1 < 1.0:
            score2 = 1.0
        return allele1, allele2, score1, score2

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
        """Determine thiopurine-related phenotype from activity score.

        TPMT and NUDT15 do NOT share a rule, and this used to apply one to
        both ("0 -> Deficient, <=1.0 -> Intermediate, else Normal").  From
        CPIC's own diplotype tables (api.cpicpgx.org/v1/diplotype, retrieved
        2026-11-03; 1,225 TPMT and 276 NUDT15 rows, all reproduced by the
        bounds below):

            pair                   TPMT            NUDT15
            no / no                Poor            Poor
            no / decreased         **Intermediate**  **Poor**
            no / normal            Intermediate    Intermediate
            decreased / decreased  Intermediate    Intermediate
            decreased / normal     **Normal**      **Normal**
            normal / normal        Normal          Normal

        The old shared rule put decreased/normal (AS 1.5) in Intermediate
        for both genes and no/decreased (AS 0.5) in Intermediate for both.
        """
        if gene == "NUDT15":
            deficient_below, intermediate_below = 1.0, 1.5
        else:  # TPMT
            deficient_below, intermediate_below = 0.5, 1.5
        if activity_score < deficient_below:
            return f"{gene} Deficient"
        if activity_score < intermediate_below:
            return f"{gene} Intermediate"
        return f"{gene} Normal"
