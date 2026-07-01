"""Disease trajectory detection engine.

Detects pre-symptomatic disease trajectories across 9 disease categories
using genotype-stratified biomarker thresholds. Pure computation — no LLM
or database calls.

Author: Adam Jones
Date: March 2026
"""

import math
from typing import Any, Dict, List, Optional

from loguru import logger

from src.knowledge import GENOTYPE_THRESHOLDS

# ---------------------------------------------------------------------------
# Disease configuration: biomarkers, genetic modifiers, thresholds, stages
# ---------------------------------------------------------------------------

DISEASE_CONFIGS = {
    "type2_diabetes": {
        "display_name": "Type 2 Diabetes",
        "biomarkers": ["hba1c", "fasting_glucose", "fasting_insulin", "homa_ir"],
        "genetic_modifiers": {
            "TCF7L2_rs7903146": {"risk_allele": "T", "effect": "beta_cell_dysfunction"},
            "PPARG_rs1801282": {"risk_allele": "C", "effect": "insulin_sensitivity"},
            "SLC30A8_rs13266634": {"risk_allele": "C", "effect": "zinc_transport"},
            "KCNJ11_rs5219": {"risk_allele": "T", "effect": "potassium_channel"},
            "GCKR_rs780094": {"risk_allele": "T", "effect": "glucokinase_regulation"},
        },
        "stages": ["normal", "early_metabolic_shift", "insulin_resistance", "pre_diabetic", "diabetic"],
    },
    "cardiovascular": {
        "display_name": "Cardiovascular Disease",
        "biomarkers": ["lpa", "ldl_c", "apob", "hs_crp", "total_cholesterol", "hdl_c", "triglycerides"],
        "genetic_modifiers": {
            "APOE": {"risk_allele": "E4", "effect": "lipid_metabolism"},
            "PCSK9_rs11591147": {"risk_allele": "T", "effect": "ldl_clearance"},
            "LPA_genetic": {"risk_allele": "elevated", "effect": "lpa_production"},
            "IL6": {"risk_allele": "C", "effect": "inflammation"},
        },
        "stages": ["optimal", "borderline", "elevated_risk", "high_risk"],
    },
    "liver": {
        "display_name": "Liver Disease (NAFLD/Fibrosis)",
        "biomarkers": ["alt", "ast", "ggt", "ferritin", "platelets", "albumin"],
        "genetic_modifiers": {
            "PNPLA3_rs738409": {"risk_allele": "G", "effect": "hepatic_lipid_accumulation"},
            "TM6SF2_rs58542926": {"risk_allele": "T", "effect": "vldl_secretion"},
            "HSD17B13_rs72613567": {"risk_allele": "TA", "effect": "protective_lipid_droplet"},
        },
        "stages": ["normal", "steatosis_risk", "early_fibrosis", "advanced_fibrosis"],
    },
    "thyroid": {
        "display_name": "Thyroid Dysfunction",
        "biomarkers": ["tsh", "free_t4", "free_t3"],
        "genetic_modifiers": {
            "DIO2_rs225014": {"risk_allele": "A", "effect": "impaired_t4_t3_conversion"},
            "DIO1_rs2235544": {"risk_allele": "A", "effect": "deiodinase_activity"},
        },
        "stages": ["euthyroid", "subclinical", "overt_dysfunction"],
    },
    "iron": {
        "display_name": "Iron Metabolism Disorder",
        "biomarkers": ["ferritin", "transferrin_saturation", "serum_iron", "tibc"],
        "genetic_modifiers": {
            "HFE_rs1800562": {"risk_allele": "A", "effect": "c282y_hemochromatosis"},
            "HFE_rs1799945": {"risk_allele": "G", "effect": "h63d_iron_loading"},
        },
        "stages": ["normal", "early_accumulation", "iron_overload"],
    },
    "nutritional": {
        "display_name": "Nutritional Deficiency",
        "biomarkers": ["omega3_index", "vitamin_d", "vitamin_b12", "folate", "magnesium", "zinc", "selenium"],
        "genetic_modifiers": {
            "FADS1_rs174546": {"risk_allele": "C", "effect": "reduced_omega3_conversion"},
            "FADS2_rs1535": {"risk_allele": "A", "effect": "fatty_acid_desaturation"},
            "VDR": {"risk_allele": "variant", "effect": "vitamin_d_receptor"},
            "MTHFR_rs1801133": {"risk_allele": "T", "effect": "folate_metabolism"},
        },
        "stages": ["optimal", "suboptimal", "deficient"],
    },
    "kidney": {
        "display_name": "Chronic Kidney Disease",
        "biomarkers": ["egfr", "creatinine", "cystatin_c", "bun", "urine_acr", "potassium"],
        "genetic_modifiers": {
            "APOL1": {"risk_allele": "G1/G2", "effect": "focal_segmental_glomerulosclerosis"},
            "UMOD_rs12917707": {"risk_allele": "G", "effect": "uromodulin_ckd_risk"},
        },
        "stages": ["normal", "early_decline", "moderate_ckd", "advanced_ckd"],
    },
    "bone_health": {
        "display_name": "Osteoporosis / Bone Health",
        "biomarkers": ["vitamin_d_25oh", "calcium", "pth", "alkaline_phosphatase", "ctx", "p1np"],
        "genetic_modifiers": {
            "VDR_rs2228570": {"risk_allele": "T", "effect": "vitamin_d_receptor_function"},
            "COL1A1_rs1800012": {"risk_allele": "T", "effect": "collagen_structure"},
            "ESR1_rs2234693": {"risk_allele": "T", "effect": "estrogen_receptor_bone"},
        },
        "stages": ["optimal", "osteopenia_risk", "osteoporosis_risk", "fracture_risk"],
    },
    "cognitive": {
        "display_name": "Cognitive Decline / Alzheimer's Risk",
        "biomarkers": ["homocysteine", "vitamin_b12", "folate", "omega3_index", "hs_crp", "vitamin_d_25oh"],
        "genetic_modifiers": {
            "APOE": {"risk_allele": "E4", "effect": "amyloid_clearance"},
            "MTHFR_rs1801133": {"risk_allele": "T", "effect": "homocysteine_metabolism"},
            "BDNF_rs6265": {"risk_allele": "A", "effect": "neurotrophin_signaling"},
        },
        "stages": ["low_risk", "modifiable_risk_factors", "elevated_risk", "high_risk"],
    },
}


def _count_risk_alleles(genotype: str, risk_allele: str) -> int:
    """Count risk alleles in a genotype string like 'CT' or 'TT'."""
    if not genotype or not risk_allele:
        return 0
    return genotype.upper().count(risk_allele.upper())


def _calculate_homa_ir(fasting_insulin: float, fasting_glucose: float) -> float:
    """Calculate HOMA-IR (Homeostatic Model Assessment of Insulin Resistance).

    Requires: insulin in µIU/mL, glucose in mg/dL.
    Formula: (insulin × glucose) / 405.0

    Note: If glucose is provided in mmol/L (typical values 3-7),
    results will be ~18x too low. This function assumes mg/dL.
    """
    # Note: heuristic-based unit detection — fasting_glucose < 20 assumed mmol/L.
    # Patients on medication may have actual 15-18 mg/dL glucose, causing false conversion.
    # Explicit unit specification recommended for production use.
    if fasting_glucose < 20:
        logger.warning(
            f"HOMA-IR: fasting_glucose={fasting_glucose} appears to be in mmol/L "
            f"(expected mg/dL, typical range 70-130). Auto-converting to mg/dL. "
            f"If this is an actual mg/dL value, the conversion is incorrect."
        )
        fasting_glucose = fasting_glucose * 18.016
    return (fasting_insulin * fasting_glucose) / 405.0


_RISK_ORDER = {"LOW": 0, "MODERATE": 1, "HIGH": 2, "CRITICAL": 3}

def _max_risk(*levels: str) -> str:
    """Return the highest risk level from the given values."""
    return max(levels, key=lambda x: _RISK_ORDER.get(x, -1))


class DiseaseTrajectoryAnalyzer:
    """Detects pre-symptomatic disease trajectories using genotype-stratified
    biomarker thresholds across 6 disease categories.

    Each disease analyzer applies genetic modifiers to shift detection
    thresholds earlier for high-risk genotypes, enabling intervention
    before clinical symptoms manifest.
    """

    def analyze_type2_diabetes(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
    ) -> Dict[str, Any]:
        """Analyze type 2 diabetes trajectory.

        Applies TCF7L2-stratified thresholds for HbA1c and HOMA-IR.
        Calculates HOMA-IR from fasting insulin + glucose if not provided.

        Args:
            biomarkers: Dict with hba1c, fasting_glucose, fasting_insulin,
                and/or homa_ir values.
            genotypes: Dict with genetic marker genotypes (e.g.,
                TCF7L2_rs7903146 -> "CT").

        Returns:
            Dict with disease, risk_level, stage, markers, genetic factors,
            and recommendations.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["type2_diabetes"]["biomarkers"]}

        # Calculate HOMA-IR if insulin and glucose provided but HOMA-IR not
        if "homa_ir" not in available:
            if "fasting_insulin" in available and "fasting_glucose" in available:
                available["homa_ir"] = round(
                    _calculate_homa_ir(available["fasting_insulin"],
                                       available["fasting_glucose"]), 2
                )

        # Determine TCF7L2 genotype and stratified thresholds
        tcf7l2 = genotypes.get("TCF7L2_rs7903146", "")
        risk_allele_count = _count_risk_alleles(tcf7l2, "T")

        # Genotype-stratified thresholds (shared constants from knowledge.py)
        # Florez et al. 2012 PMID:22399527; Grant et al. 2006 PMID:16415884
        thresholds = {
            2: {  # TT: tightest
                "hba1c": GENOTYPE_THRESHOLDS["TCF7L2_hba1c"][2],  # ADA Standards of Care 2024; PMID:38078589
                "homa_ir": 3.0,   # Matthews et al. 1985 PMID:3899825
                "fasting_glucose": GENOTYPE_THRESHOLDS["TCF7L2_fasting_glucose"][2],  # ADA Standards of Care 2024; PMID:38078589
            },
            1: {  # CT: intermediate
                "hba1c": GENOTYPE_THRESHOLDS["TCF7L2_hba1c"][1],  # ADA Standards of Care 2024; PMID:38078589
                "homa_ir": 3.5,   # Matthews et al. 1985 PMID:3899825
                "fasting_glucose": GENOTYPE_THRESHOLDS["TCF7L2_fasting_glucose"][1],  # ADA Standards of Care 2024; PMID:38078589
            },
            0: {  # CC: standard ADA
                "hba1c": GENOTYPE_THRESHOLDS["TCF7L2_hba1c"][0],  # ADA Standards of Care 2024; PMID:38078589
                "homa_ir": 4.0,   # Matthews et al. 1985 PMID:3899825
                "fasting_glucose": GENOTYPE_THRESHOLDS["TCF7L2_fasting_glucose"][0],  # ADA Standards of Care 2024; PMID:38078589
            },
        }
        thresh = thresholds.get(risk_allele_count, thresholds[0])

        # Collect genetic risk factors
        genetic_risk_factors = []
        for gene_key, config in DISEASE_CONFIGS["type2_diabetes"]["genetic_modifiers"].items():
            gt = genotypes.get(gene_key, "")
            if gt:
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key,
                        "genotype": gt,
                        "risk_alleles": count,
                        "effect": config["effect"],
                    })

        # Determine stage
        stage = "normal"
        risk_level = "LOW"
        findings = []
        recommendations = []

        hba1c = available.get("hba1c")
        homa_ir = available.get("homa_ir")
        fasting_glucose = available.get("fasting_glucose")

        if hba1c is not None:
            if hba1c >= 6.5:  # ADA Standards of Care 2024; PMID:38078589
                stage = "diabetic"
                risk_level = "CRITICAL"
                findings.append(f"HbA1c {hba1c}% is in diabetic range (>=6.5%)")
                recommendations.append("Immediate endocrinology referral for diabetes management")
            elif hba1c >= 5.7:  # ADA Standards of Care 2024; PMID:38078589
                stage = "pre_diabetic"
                risk_level = "HIGH"
                findings.append(f"HbA1c {hba1c}% is in pre-diabetic range (5.7-6.4%)")
                recommendations.append("Intensive lifestyle intervention: diet + 150 min/week exercise")
                recommendations.append("Consider metformin if BMI >35 or age <60")
            elif hba1c >= thresh["hba1c"]:
                stage = "early_metabolic_shift"
                risk_level = "MODERATE"
                findings.append(
                    f"HbA1c {hba1c}% exceeds genotype-adjusted threshold "
                    f"({thresh['hba1c']}% for {tcf7l2 or 'unknown'} genotype)"
                )
                recommendations.append("Dietary optimization: reduce refined carbohydrates")
                recommendations.append("Recheck HbA1c and fasting insulin in 3 months")

        if homa_ir is not None:
            if homa_ir >= thresh["homa_ir"]:
                if stage in ("normal", "early_metabolic_shift"):
                    stage = "insulin_resistance" if stage == "normal" else stage
                risk_level = _max_risk(risk_level, "MODERATE")
                findings.append(
                    f"HOMA-IR {homa_ir} exceeds genotype-adjusted threshold "
                    f"({thresh['homa_ir']} for {tcf7l2 or 'unknown'} genotype)"
                )
                recommendations.append("Assess insulin resistance with oral glucose tolerance test")
                recommendations.append("Consider berberine or inositol supplementation")

        if fasting_glucose is not None:
            if fasting_glucose >= 126:
                if stage != "diabetic":
                    stage = "diabetic"
                    risk_level = "CRITICAL"
                findings.append(f"Fasting glucose {fasting_glucose} mg/dL is in diabetic range (>=126)")
            elif fasting_glucose >= thresh["fasting_glucose"]:
                findings.append(
                    f"Fasting glucose {fasting_glucose} mg/dL exceeds genotype-adjusted "
                    f"threshold ({thresh['fasting_glucose']} mg/dL)"
                )

        if genetic_risk_factors and not findings:
            findings.append("Genetic risk factors present; recommend closer monitoring")
            recommendations.append("Annual HbA1c and fasting insulin monitoring")

        return {
            "disease": "type2_diabetes",
            "display_name": "Type 2 Diabetes",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "genotype_adjusted_thresholds": thresh,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_cardiovascular(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
    ) -> Dict[str, Any]:
        """Analyze cardiovascular disease trajectory.

        Applies APOE-stratified LDL thresholds and evaluates Lp(a),
        ApoB, hs-CRP, and lipid panel markers.

        Args:
            biomarkers: Dict with lpa, ldl_c, apob, hs_crp, etc.
            genotypes: Dict with APOE, PCSK9_rs11591147, etc.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["cardiovascular"]["biomarkers"]}

        # Determine APOE genotype and LDL threshold
        apoe = genotypes.get("APOE", "")
        # APOE-stratified LDL targets; Bennet et al. 2007 PMID:17327455
        # ACC/AHA 2019 Cholesterol Guideline PMID:30586774
        ldl_thresholds = {
            "E4/E4": 100, "E3/E4": 115, "E4/E3": 115,
            "E3/E3": 130, "E2/E3": 130, "E3/E2": 130,
            "E2/E2": 130,  # Special case: type III dyslipidemia
        }
        ldl_thresh = ldl_thresholds.get(apoe, 130)
        has_e4 = "E4" in apoe

        genetic_risk_factors = []
        if has_e4:
            genetic_risk_factors.append({
                "gene": "APOE",
                "genotype": apoe,
                "risk_alleles": apoe.count("E4"),
                "effect": "lipid_metabolism",
            })

        for gene_key in ("PCSK9_rs11591147", "LPA_genetic", "IL6"):
            gt = genotypes.get(gene_key, "")
            if gt:
                config = DISEASE_CONFIGS["cardiovascular"]["genetic_modifiers"][gene_key]
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key,
                        "genotype": gt,
                        "risk_alleles": count,
                        "effect": config["effect"],
                    })

        stage = "optimal"
        risk_level = "LOW"
        findings = []
        recommendations = []

        # Lp(a) — genetically determined
        # ESC/EAS 2019 Guidelines PMID:31504418; >50 nmol/L is elevated
        lpa = available.get("lpa")
        if lpa is not None:
            if lpa > 125:  # Nordestgaard et al. 2010 PMID:20031622
                stage = "high_risk"
                risk_level = "HIGH"
                findings.append(f"Lp(a) {lpa} nmol/L is very elevated (>125, strongly genetic)")
                recommendations.append("Lp(a) is largely genetic; aggressive LDL lowering recommended")
                recommendations.append("Consider PCSK9 inhibitor therapy")
            elif lpa > 50:  # ESC/EAS 2019 Guidelines PMID:31504418
                stage = "elevated_risk"
                risk_level = "MODERATE"
                findings.append(f"Lp(a) {lpa} nmol/L is elevated (>50)")
                recommendations.append("Optimize all modifiable cardiovascular risk factors")

        # LDL-C with APOE-adjusted threshold
        ldl = available.get("ldl_c")
        if ldl is not None:
            if ldl >= 190:  # ACC/AHA 2019 Cholesterol Guideline PMID:30586774
                stage = "high_risk"
                risk_level = "HIGH"
                findings.append(f"LDL-C {ldl} mg/dL is very high (>=190)")
                recommendations.append("High-intensity statin therapy recommended")
            elif ldl >= ldl_thresh:
                if stage in ("optimal", "borderline"):
                    stage = "elevated_risk"
                risk_level = _max_risk(risk_level, "MODERATE")
                findings.append(
                    f"LDL-C {ldl} mg/dL exceeds APOE-adjusted threshold "
                    f"({ldl_thresh} mg/dL for {apoe or 'standard'})"
                )
                if has_e4:
                    recommendations.append("APOE E4 carrier: aggressive LDL lowering to <100 mg/dL")
                else:
                    recommendations.append("Lifestyle modifications: diet, exercise, consider statin")

        # ApoB — more predictive than LDL in some cases
        apob = available.get("apob")
        if apob is not None:
            if apob > 130:
                findings.append(f"ApoB {apob} mg/dL is elevated (>130)")
                if stage == "optimal":
                    stage = "borderline"
            elif apob > 90 and has_e4:
                findings.append(f"ApoB {apob} mg/dL may warrant attention in APOE E4 carrier")

        # hs-CRP — inflammatory risk
        crp = available.get("hs_crp")
        if crp is not None:
            if crp > 3.0:
                findings.append(f"hs-CRP {crp} mg/L indicates high inflammatory risk (>3.0)")
                recommendations.append("Evaluate sources of chronic inflammation")
                if stage in ("optimal", "borderline"):
                    stage = "elevated_risk"
            elif crp > 1.0:
                findings.append(f"hs-CRP {crp} mg/L is borderline (1.0-3.0)")

        # APOE E2/E2 special case
        if apoe == "E2/E2":
            trig = available.get("triglycerides")
            if trig is not None and trig > 200:
                findings.append(
                    "APOE E2/E2 with elevated triglycerides: risk of type III dyslipidemia"
                )
                recommendations.append("Screen for type III hyperlipoproteinemia")

        return {
            "disease": "cardiovascular",
            "display_name": "Cardiovascular Disease",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "ldl_threshold": ldl_thresh,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_liver(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
        age: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Analyze liver disease trajectory (NAFLD/fibrosis).

        Calculates FIB-4 index and applies PNPLA3-stratified ALT thresholds.

        Args:
            biomarkers: Dict with alt, ast, ggt, ferritin, platelets, albumin.
            genotypes: Dict with PNPLA3_rs738409, TM6SF2_rs58542926, etc.
            age: Patient age in years (needed for FIB-4 calculation).

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["liver"]["biomarkers"]}

        # PNPLA3 genotype determines ALT upper limit
        pnpla3 = genotypes.get("PNPLA3_rs738409", "")
        pnpla3_risk = _count_risk_alleles(pnpla3, "G")
        # PNPLA3-adjusted ALT limits
        # Romeo et al. 2008 PMID:18820127; Sookoian & Pirola 2011 PMID:21520172
        # ALT thresholds: Prati et al. 2002 PMID:12029600
        alt_upper = GENOTYPE_THRESHOLDS["PNPLA3_alt_upper"].get(
            "GG" if pnpla3_risk == 2 else "CG" if pnpla3_risk == 1 else "CC", 56
        )

        genetic_risk_factors = []
        for gene_key, config in DISEASE_CONFIGS["liver"]["genetic_modifiers"].items():
            gt = genotypes.get(gene_key, "")
            if gt:
                count = _count_risk_alleles(gt, config["risk_allele"])
                is_protective = (gene_key == "HSD17B13_rs72613567" and "TA" in gt)
                if count > 0 or is_protective:
                    genetic_risk_factors.append({
                        "gene": gene_key,
                        "genotype": gt,
                        "risk_alleles": count,
                        "effect": config["effect"],
                        "direction": "protective" if is_protective else "risk",
                    })

        # Calculate FIB-4 if age, AST, platelets, ALT available
        fib4 = None
        alt = available.get("alt")
        ast = available.get("ast")
        platelets = available.get("platelets")

        # FIB-4 = (age × AST) / (platelets × √ALT); requires platelets > 0 and alt > 0
        if age is not None and age > 0 and ast is not None and platelets is not None and platelets > 0 and alt is not None and alt > 0:
            fib4 = round((age * ast) / (platelets * math.sqrt(alt)), 2)
            available["fib4_calculated"] = fib4

        stage = "normal"
        risk_level = "LOW"
        findings = []
        recommendations = []

        # ALT with PNPLA3-adjusted threshold
        if alt is not None:
            if alt > alt_upper:
                stage = "steatosis_risk"
                risk_level = "MODERATE"
                findings.append(
                    f"ALT {alt} U/L exceeds PNPLA3-adjusted upper limit "
                    f"({alt_upper} U/L for {'GG' if pnpla3_risk == 2 else 'CG' if pnpla3_risk == 1 else 'CC'} genotype)"
                )
                recommendations.append("Liver ultrasound to evaluate for hepatic steatosis")

        # AST/ALT ratio
        if ast is not None and alt is not None and alt > 0:
            ratio = ast / alt
            if ratio > 1.0:
                findings.append(f"AST/ALT ratio {ratio:.1f} (>1.0 may indicate fibrosis)")

        # GGT
        ggt = available.get("ggt")
        if ggt is not None and ggt > 60:
            findings.append(f"GGT {ggt} U/L is elevated (>60)")
            recommendations.append("Evaluate alcohol consumption and medication hepatotoxicity")

        # FIB-4 staging
        if fib4 is not None:
            if fib4 > 2.67:  # Sterling et al. 2006 PMID:16702586; AASLD 2023 Practice Guidelines
                stage = "advanced_fibrosis"
                risk_level = "HIGH"
                findings.append(f"FIB-4 score {fib4} suggests advanced fibrosis (>2.67)")
                recommendations.append("Urgent hepatology referral; consider FibroScan or liver biopsy")
            elif fib4 > 1.30:  # Sterling et al. 2006 PMID:16702586; AASLD 2023 Practice Guidelines
                if stage != "advanced_fibrosis":
                    stage = "early_fibrosis"
                risk_level = _max_risk(risk_level, "MODERATE")
                findings.append(f"FIB-4 score {fib4} is indeterminate (1.30-2.67)")
                recommendations.append("FibroScan recommended; recheck in 6 months")

        # Ferritin — elevated in NAFLD
        ferritin = available.get("ferritin")
        if ferritin is not None and ferritin > 300:
            findings.append(f"Ferritin {ferritin} ng/mL is elevated (>300)")
            recommendations.append("Check transferrin saturation to differentiate NAFLD vs hemochromatosis")

        # Albumin — low suggests hepatic synthetic dysfunction
        albumin = available.get("albumin")
        if albumin is not None and albumin < 3.5:
            findings.append(f"Albumin {albumin} g/dL is low (<3.5), suggesting impaired liver synthesis")
            risk_level = "HIGH"

        # PNPLA3 GG specific guidance
        if pnpla3_risk == 2:
            recommendations.append(
                "PNPLA3 GG carrier: 2-3x NAFLD risk. Strict weight management "
                "and avoidance of high-fructose diet strongly recommended"
            )

        # HSD17B13 protective effect
        hsd = genotypes.get("HSD17B13_rs72613567", "")
        if "TA" in hsd:
            findings.append("HSD17B13 TA variant detected: protective against NASH progression")

        return {
            "disease": "liver",
            "display_name": "Liver Disease (NAFLD/Fibrosis)",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "alt_upper_limit": alt_upper,
            "fib4_score": fib4,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_thyroid(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
    ) -> Dict[str, Any]:
        """Analyze thyroid dysfunction trajectory.

        Applies DIO2-stratified TSH/free_T3 thresholds to detect impaired
        T4-to-T3 conversion that standard testing may miss.

        Args:
            biomarkers: Dict with tsh, free_t4, free_t3.
            genotypes: Dict with DIO2_rs225014, DIO1_rs2235544.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["thyroid"]["biomarkers"]}

        # DIO2 genotype adjusts TSH upper limit and free_T3 lower limit
        dio2 = genotypes.get("DIO2_rs225014", "")
        dio2_risk = _count_risk_alleles(dio2, "A")

        # DIO2-adjusted reference ranges
        # Panicker et al. 2009 PMID:19820026; Castagna et al. 2017 PMID:28100792
        dio2_label = "AA" if dio2_risk == 2 else "GA" if dio2_risk == 1 else "GG"
        tsh_upper = GENOTYPE_THRESHOLDS["DIO2_tsh_upper"].get(dio2_label, 4.0)  # ATA 2017 Guidelines PMID:28056690
        ft3_lower = {0: 2.3, 1: 2.5, 2: 2.8}.get(dio2_risk, 2.3)

        genetic_risk_factors = []
        for gene_key, config in DISEASE_CONFIGS["thyroid"]["genetic_modifiers"].items():
            gt = genotypes.get(gene_key, "")
            if gt:
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key,
                        "genotype": gt,
                        "risk_alleles": count,
                        "effect": config["effect"],
                    })

        stage = "euthyroid"
        risk_level = "LOW"
        findings = []
        recommendations = []

        tsh = available.get("tsh")
        ft4 = available.get("free_t4")
        ft3 = available.get("free_t3")

        # TSH evaluation with DIO2 adjustment
        if tsh is not None:
            if tsh > 10.0:
                stage = "overt_dysfunction"
                risk_level = "HIGH"
                findings.append(f"TSH {tsh} mIU/L indicates overt hypothyroidism (>10.0)")
                recommendations.append("Endocrinology referral; levothyroxine initiation likely needed")
            elif tsh < 0.1:
                stage = "overt_dysfunction"
                risk_level = "HIGH"
                findings.append(f"TSH {tsh} mIU/L indicates possible hyperthyroidism (<0.1)")
                recommendations.append("Urgent evaluation: free T4, free T3, TSI/TRAb")
            elif tsh > tsh_upper:
                stage = "subclinical"
                risk_level = "MODERATE"
                findings.append(
                    f"TSH {tsh} mIU/L exceeds DIO2-adjusted upper limit "
                    f"({tsh_upper} mIU/L for {'AA' if dio2_risk == 2 else 'GA' if dio2_risk == 1 else 'GG'} genotype)"
                )
                recommendations.append("Recheck TSH + anti-TPO antibodies in 6-8 weeks")

        # Free T3 — critical for DIO2 variants
        if ft3 is not None:
            if ft3 < ft3_lower:
                if stage == "euthyroid":
                    stage = "subclinical"
                risk_level = _max_risk(risk_level, "MODERATE")
                findings.append(
                    f"Free T3 {ft3} pg/mL is below DIO2-adjusted threshold "
                    f"({ft3_lower} pg/mL for {'AA' if dio2_risk == 2 else 'GA' if dio2_risk == 1 else 'GG'} genotype)"
                )
                if dio2_risk == 2:
                    recommendations.append(
                        "DIO2 AA genotype with low free T3: impaired T4->T3 conversion. "
                        "Consider combination T4+T3 therapy if symptomatic"
                    )
                else:
                    recommendations.append("Monitor free T3 levels; consider T3 supplementation if symptomatic")

        # DIO2 AA with normal TSH but not-checked T3 — flag for monitoring
        if dio2_risk == 2 and tsh is not None and ft3 is None:
            findings.append(
                "DIO2 AA genotype detected but free T3 not measured. "
                "TSH alone may miss impaired T4->T3 conversion"
            )
            recommendations.append("Order free T3 measurement — essential for DIO2 AA genotype")

        # Free T4
        if ft4 is not None:
            if ft4 < 0.8:
                findings.append(f"Free T4 {ft4} ng/dL is low (<0.8)")
            elif ft4 > 1.8:
                findings.append(f"Free T4 {ft4} ng/dL is elevated (>1.8)")

        return {
            "disease": "thyroid",
            "display_name": "Thyroid Dysfunction",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "tsh_upper_limit": tsh_upper,
            "ft3_lower_limit": ft3_lower,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_iron(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
        sex: str = "male",
    ) -> Dict[str, Any]:
        """Analyze iron metabolism disorder trajectory.

        Applies HFE C282Y-stratified ferritin thresholds and evaluates
        transferrin saturation for hemochromatosis risk.

        Args:
            biomarkers: Dict with ferritin, transferrin_saturation,
                serum_iron, tibc.
            genotypes: Dict with HFE_rs1800562, HFE_rs1799945.
            sex: Patient sex ('male' or 'female') for reference ranges.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["iron"]["biomarkers"]}

        # HFE C282Y genotype determines ferritin thresholds
        hfe_c282y = genotypes.get("HFE_rs1800562", "")
        c282y_risk = _count_risk_alleles(hfe_c282y, "A")

        # HFE-adjusted ferritin limits
        # WHO Ferritin Guidelines 2020; EASL HFE Guidelines PMID:20471131
        # Adams et al. 2005 PMID:15729334; European Clinical Practice Guidelines PMID:20471131
        if sex.lower() == "female":
            ferritin_upper = {0: 150, 1: 120, 2: 100}.get(c282y_risk, 150)  # WHO Ferritin Guidelines 2020; EASL HFE Guidelines PMID:20471131
        else:
            ferritin_upper = {0: 400, 1: 300, 2: 200}.get(c282y_risk, 400)  # WHO Ferritin Guidelines 2020; EASL HFE Guidelines PMID:20471131

        genetic_risk_factors = []
        for gene_key, config in DISEASE_CONFIGS["iron"]["genetic_modifiers"].items():
            gt = genotypes.get(gene_key, "")
            if gt:
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key,
                        "genotype": gt,
                        "risk_alleles": count,
                        "effect": config["effect"],
                    })

        stage = "normal"
        risk_level = "LOW"
        findings = []
        recommendations = []

        ferritin = available.get("ferritin")
        tsat = available.get("transferrin_saturation")

        # Ferritin with HFE-adjusted thresholds
        if ferritin is not None:
            if ferritin > ferritin_upper:
                stage = "early_accumulation"
                risk_level = "MODERATE"
                genotype_label = "AA" if c282y_risk == 2 else "GA" if c282y_risk == 1 else "GG"
                findings.append(
                    f"Ferritin {ferritin} ng/mL exceeds HFE-adjusted upper limit "
                    f"({ferritin_upper} ng/mL for {genotype_label} genotype, {sex})"
                )
                recommendations.append("Monitor ferritin every 3-6 months")

        # Transferrin saturation
        if tsat is not None:
            if tsat > 45:
                if stage == "normal":
                    stage = "early_accumulation"
                risk_level = _max_risk(risk_level, "MODERATE")
                findings.append(f"Transferrin saturation {tsat}% is elevated (>45%)")
                recommendations.append("Evaluate for iron overload; consider HFE genotyping if not done")
            if tsat > 60:
                stage = "iron_overload"
                risk_level = "HIGH"
                findings.append(f"Transferrin saturation {tsat}% is very elevated (>60%)")
                recommendations.append("Urgent hematology referral for therapeutic phlebotomy evaluation")

        # C282Y homozygous — highest risk
        if c282y_risk == 2:
            findings.append(
                "HFE C282Y homozygous (AA): hereditary hemochromatosis genotype. "
                "Lifetime risk of iron overload is significant"
            )
            recommendations.append("Annual ferritin + transferrin saturation monitoring")
            recommendations.append("Liver MRI for iron quantification if ferritin >500")
            if ferritin is not None and ferritin > 200:  # Adams et al. 2005 PMID:15729334; European Clinical Practice Guidelines PMID:20471131
                stage = "iron_overload"
                risk_level = "HIGH"
                recommendations.append("Consider therapeutic phlebotomy")

        # H63D compound heterozygote
        hfe_h63d = genotypes.get("HFE_rs1799945", "")
        if c282y_risk >= 1 and _count_risk_alleles(hfe_h63d, "G") >= 1:
            findings.append("C282Y/H63D compound heterozygote: moderate hemochromatosis risk")

        return {
            "disease": "iron",
            "display_name": "Iron Metabolism Disorder",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "ferritin_upper_limit": ferritin_upper,
            "sex": sex,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_nutritional(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
    ) -> Dict[str, Any]:
        """Analyze nutritional deficiency trajectory.

        Applies FADS1/MTHFR-stratified thresholds for omega-3 index
        and folate metabolism.

        Args:
            biomarkers: Dict with omega3_index, vitamin_d, vitamin_b12,
                folate, magnesium, zinc, selenium.
            genotypes: Dict with FADS1_rs174546, MTHFR_rs1801133, etc.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["nutritional"]["biomarkers"]}

        # FADS1 genotype for omega-3 targets
        fads1 = genotypes.get("FADS1_rs174546", "")
        fads1_risk = _count_risk_alleles(fads1, "C")
        # FADS1-adjusted omega-3 targets (Schaeffer 2006 PMID:16825694; Harris 2018 PMID:29215303)
        omega3_target = {0: 8.0, 1: 6.0, 2: 5.0}.get(fads1_risk, 8.0)  # Harris & von Schacky 2004 PMID:14676843

        # MTHFR genotype for folate metabolism
        mthfr = genotypes.get("MTHFR_rs1801133", "")
        mthfr_risk = _count_risk_alleles(mthfr, "T")

        genetic_risk_factors = []
        for gene_key, config in DISEASE_CONFIGS["nutritional"]["genetic_modifiers"].items():
            gt = genotypes.get(gene_key, "")
            if gt:
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key,
                        "genotype": gt,
                        "risk_alleles": count,
                        "effect": config["effect"],
                    })

        stage = "optimal"
        risk_level = "LOW"
        findings = []
        recommendations = []
        deficiencies = []

        # Omega-3 index with FADS1 adjustment
        omega3 = available.get("omega3_index")
        if omega3 is not None:
            if omega3 < omega3_target:
                deficiencies.append("omega3")
                genotype_label = "CC" if fads1_risk == 2 else "CT" if fads1_risk == 1 else "TT"
                findings.append(
                    f"Omega-3 index {omega3}% below FADS1-adjusted target "
                    f"({omega3_target}% for {genotype_label} genotype)"
                )
                if fads1_risk == 2:
                    recommendations.append(
                        "FADS1 CC genotype: severely reduced ALA->EPA/DHA conversion. "
                        "Needs direct EPA/DHA supplementation (fish oil 2-4g/day). "
                        "Plant-based omega-3 sources insufficient"
                    )
                elif fads1_risk == 1:
                    recommendations.append("Increase dietary EPA/DHA; consider fish oil 1-2g/day")
                else:
                    recommendations.append("Increase fatty fish consumption or supplement omega-3")

        # Vitamin D
        vit_d = available.get("vitamin_d")
        if vit_d is not None:
            if vit_d < 20:
                deficiencies.append("vitamin_d")
                findings.append(f"Vitamin D {vit_d} ng/mL is deficient (<20)")
                recommendations.append("Vitamin D3 supplementation: 2000-5000 IU/day with K2")
            elif vit_d < 30:  # Endocrine Society 2024 PMID:38828931
                findings.append(f"Vitamin D {vit_d} ng/mL is insufficient (20-30)")
                recommendations.append("Vitamin D3 supplementation: 1000-2000 IU/day")

        # Vitamin B12
        b12 = available.get("vitamin_b12")
        if b12 is not None:
            if b12 < 200:
                deficiencies.append("vitamin_b12")
                findings.append(f"Vitamin B12 {b12} pg/mL is deficient (<200)")
                recommendations.append("B12 supplementation: methylcobalamin 1000-2000 mcg/day")
            elif b12 < 400:
                findings.append(f"Vitamin B12 {b12} pg/mL is suboptimal (<400)")
                recommendations.append("Consider B12 supplementation: methylcobalamin 500-1000 mcg/day")

        # Folate with MTHFR adjustment; Frosst et al. 1995 PMID:7647779
        folate = available.get("folate")
        if folate is not None:
            if mthfr_risk == 2:
                # MTHFR TT: needs L-methylfolate, not folic acid
                findings.append(
                    f"MTHFR C677T TT genotype detected. Folate level: {folate} ng/mL"
                )
                recommendations.append(
                    "MTHFR TT genotype: use L-methylfolate (1-5mg/day) instead of folic acid. "
                    "Folic acid may accumulate unmetabolized"
                )
                if folate < 10:
                    deficiencies.append("folate")
                    recommendations.append("Prioritize L-methylfolate supplementation")
            elif mthfr_risk == 1:
                if folate < 5:
                    deficiencies.append("folate")
                    findings.append(f"Folate {folate} ng/mL is low (<5)")
                    recommendations.append("MTHFR CT: consider L-methylfolate supplementation")
            else:
                if folate < 3:
                    deficiencies.append("folate")
                    findings.append(f"Folate {folate} ng/mL is deficient (<3)")
                    recommendations.append("Folate supplementation recommended")

        # Magnesium
        mag = available.get("magnesium")
        if mag is not None and mag < 1.8:
            deficiencies.append("magnesium")
            findings.append(f"Magnesium {mag} mg/dL is low (<1.8)")
            recommendations.append("Magnesium glycinate or threonate 200-400mg/day")

        # Zinc
        zinc = available.get("zinc")
        if zinc is not None and zinc < 60:
            deficiencies.append("zinc")
            findings.append(f"Zinc {zinc} mcg/dL is low (<60)")
            recommendations.append("Zinc supplementation 15-30mg/day (take with food, separate from iron)")

        # Selenium
        selenium = available.get("selenium")
        if selenium is not None and selenium < 70:
            deficiencies.append("selenium")
            findings.append(f"Selenium {selenium} mcg/L is low (<70)")
            recommendations.append("Selenium supplementation 100-200mcg/day (or 2-3 Brazil nuts/day)")

        # Determine overall stage
        if len(deficiencies) >= 3:
            stage = "deficient"
            risk_level = "HIGH"
        elif len(deficiencies) >= 1:
            stage = "suboptimal"
            risk_level = "MODERATE"

        return {
            "disease": "nutritional",
            "display_name": "Nutritional Deficiency",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "omega3_target": omega3_target,
            "deficiencies": deficiencies,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_kidney(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
        age: Optional[float] = None,
        sex: str = "male",
    ) -> Dict[str, Any]:
        """Analyze chronic kidney disease trajectory.

        Evaluates eGFR (calculated from creatinine and/or cystatin C),
        urine albumin-to-creatinine ratio, and APOL1 risk variants.

        Args:
            biomarkers: Dict with egfr, creatinine, cystatin_c, bun, urine_acr.
            genotypes: Dict with APOL1, UMOD_rs12917707.
            age: Patient age in years (needed for eGFR calculation).
            sex: Patient sex for eGFR calculation.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["kidney"]["biomarkers"]}

        # Calculate eGFR from creatinine if not provided (CKD-EPI 2021)
        # Inker et al. 2021 PMID:34554658 — race-neutral equation
        if "egfr" not in available and "creatinine" in available and age is not None:
            scr = available["creatinine"]
            if sex.lower() == "female":
                kappa, alpha, female_mult = 0.7, -0.241, 1.012
            else:
                kappa, alpha, female_mult = 0.9, -0.302, 1.0
            egfr = 142 * min(scr / kappa, 1) ** alpha * max(scr / kappa, 1) ** (-1.200) * 0.9938 ** age * female_mult
            available["egfr_calculated"] = round(egfr, 1)

        # APOL1 risk assessment (high-risk = two G1/G2 alleles)
        apol1 = genotypes.get("APOL1", "")
        has_apol1_risk = apol1.lower() in ("g1/g1", "g2/g2", "g1/g2", "high_risk", "two_risk_alleles")

        genetic_risk_factors = []
        if has_apol1_risk:
            genetic_risk_factors.append({
                "gene": "APOL1",
                "genotype": apol1,
                "risk_alleles": 2,
                "effect": "focal_segmental_glomerulosclerosis",
            })
        for gene_key in ("UMOD_rs12917707",):
            gt = genotypes.get(gene_key, "")
            if gt:
                config = DISEASE_CONFIGS["kidney"]["genetic_modifiers"][gene_key]
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key, "genotype": gt,
                        "risk_alleles": count, "effect": config["effect"],
                    })

        stage = "normal"
        risk_level = "LOW"
        findings = []
        recommendations = []

        egfr = available.get("egfr") or available.get("egfr_calculated")
        if egfr is not None:
            if egfr < 30:
                stage = "advanced_ckd"
                risk_level = "CRITICAL"
                findings.append(f"eGFR {egfr} mL/min/1.73m2 — CKD Stage 4-5 (<30)")
                recommendations.append("Urgent nephrology referral; evaluate for renal replacement therapy planning")
            elif egfr < 60:
                stage = "moderate_ckd"
                risk_level = "HIGH"
                findings.append(f"eGFR {egfr} mL/min/1.73m2 — CKD Stage 3 (30-59)")
                recommendations.append("Nephrology consultation; ACEi/ARB for renoprotection")
                recommendations.append("SGLT2 inhibitor (dapagliflozin/empagliflozin) for CKD progression reduction")
            elif egfr < 90:
                stage = "early_decline"
                risk_level = "MODERATE"
                findings.append(f"eGFR {egfr} mL/min/1.73m2 — mildly decreased (60-89)")
                recommendations.append("Monitor creatinine and urine ACR every 6 months")

        # Urine albumin-to-creatinine ratio
        uacr = available.get("urine_acr")
        if uacr is not None:
            if uacr > 300:
                risk_level = _max_risk(risk_level, "HIGH")
                findings.append(f"Urine ACR {uacr} mg/g — severely increased albuminuria (>300, A3)")
                recommendations.append("Optimize blood pressure <130/80; maximize ACEi/ARB dose")
            elif uacr > 30:
                risk_level = _max_risk(risk_level, "MODERATE")
                findings.append(f"Urine ACR {uacr} mg/g — moderately increased albuminuria (30-300, A2)")
                recommendations.append("Start or optimize ACEi/ARB therapy for albuminuria reduction")

        # Cystatin C — more accurate GFR estimation in muscle wasting
        cystatin = available.get("cystatin_c")
        if cystatin is not None and cystatin > 1.0:
            findings.append(f"Cystatin C {cystatin} mg/L elevated (>1.0) — consider cystatin-based eGFR")

        # BUN
        bun = available.get("bun")
        if bun is not None and bun > 20:
            findings.append(f"BUN {bun} mg/dL elevated (>20)")

        # APOL1 high-risk genotype
        if has_apol1_risk:
            risk_level = _max_risk(risk_level, "MODERATE")
            findings.append(
                "APOL1 high-risk genotype (two G1/G2 alleles): 7-10x increased risk of "
                "FSGS and CKD progression, especially in African ancestry populations"
            )
            recommendations.append("Annual eGFR + urine ACR monitoring due to APOL1 risk")
            recommendations.append("Avoid nephrotoxic agents; maintain aggressive blood pressure control")

        # Potassium — important in CKD
        potassium = available.get("potassium")
        if potassium is not None:
            if potassium > 5.5:
                findings.append(f"Potassium {potassium} mEq/L is elevated (>5.5) — hyperkalemia risk in CKD")
                recommendations.append("Review ACEi/ARB dose; consider potassium binder if persistent")

        return {
            "disease": "kidney",
            "display_name": "Chronic Kidney Disease",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_bone_health(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
        sex: str = "male",
        age: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Analyze osteoporosis and bone health trajectory.

        Evaluates vitamin D, calcium, PTH, bone turnover markers, and
        VDR/COL1A1 genetic modifiers.

        Args:
            biomarkers: Dict with vitamin_d_25oh, calcium, pth,
                alkaline_phosphatase, ctx (C-telopeptide), p1np.
            genotypes: Dict with VDR_rs2228570, COL1A1_rs1800012, ESR1_rs2234693.
            sex: Patient sex for risk stratification.
            age: Patient age for risk context.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["bone_health"]["biomarkers"]}

        genetic_risk_factors = []
        for gene_key, config in DISEASE_CONFIGS["bone_health"]["genetic_modifiers"].items():
            gt = genotypes.get(gene_key, "")
            if gt:
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key, "genotype": gt,
                        "risk_alleles": count, "effect": config["effect"],
                    })

        stage = "optimal"
        risk_level = "LOW"
        findings = []
        recommendations = []
        risk_factors_count = 0

        # Vitamin D — critical for calcium absorption
        vit_d = available.get("vitamin_d_25oh")
        if vit_d is not None:
            vdr = genotypes.get("VDR_rs2228570", "")
            vdr_risk = _count_risk_alleles(vdr, "T")
            # VDR variants need higher vitamin D levels
            vit_d_target = 40 if vdr_risk >= 2 else 35 if vdr_risk == 1 else 30

            if vit_d < 20:
                risk_factors_count += 1
                findings.append(f"Vitamin D {vit_d} ng/mL is deficient (<20) — bone mineralization impaired")
                recommendations.append("Vitamin D3 5000 IU/day + K2 for 8-12 weeks, then recheck")
            elif vit_d < vit_d_target:
                risk_factors_count += 1
                findings.append(
                    f"Vitamin D {vit_d} ng/mL below VDR-adjusted target "
                    f"({vit_d_target} ng/mL for {'TT' if vdr_risk >= 2 else 'CT' if vdr_risk == 1 else 'standard'})"
                )
                recommendations.append(f"Supplement vitamin D3 2000-4000 IU/day to target >{vit_d_target} ng/mL")

        # PTH — secondary hyperparathyroidism
        pth = available.get("pth")
        if pth is not None and pth > 65:
            risk_factors_count += 1
            findings.append(f"PTH {pth} pg/mL is elevated (>65) — secondary hyperparathyroidism")
            recommendations.append("Evaluate vitamin D status and calcium intake; rule out primary hyperparathyroidism")

        # Calcium
        calcium = available.get("calcium")
        if calcium is not None:
            if calcium < 8.5:
                risk_factors_count += 1
                findings.append(f"Calcium {calcium} mg/dL is low (<8.5)")
                recommendations.append("Evaluate PTH and vitamin D; consider calcium citrate supplementation")
            elif calcium > 10.5:
                findings.append(f"Calcium {calcium} mg/dL is elevated (>10.5) — evaluate for hyperparathyroidism")

        # Bone turnover markers
        ctx = available.get("ctx")  # C-telopeptide (bone resorption)
        if ctx is not None and ctx > 0.6:
            risk_factors_count += 1
            findings.append(f"CTX {ctx} ng/mL elevated (>0.6) — increased bone resorption")
            recommendations.append("Consider DEXA scan; evaluate for anti-resorptive therapy need")

        p1np = available.get("p1np")  # Procollagen I N-propeptide (bone formation)
        if p1np is not None and p1np < 25:
            risk_factors_count += 1
            findings.append(f"P1NP {p1np} mcg/L is low (<25) — reduced bone formation")

        # ALP elevated with other bone markers
        alp = available.get("alkaline_phosphatase")
        if alp is not None and alp > 120:
            findings.append(f"Alkaline phosphatase {alp} U/L elevated — evaluate bone vs liver source")

        # COL1A1 risk
        col1a1 = genotypes.get("COL1A1_rs1800012", "")
        if _count_risk_alleles(col1a1, "T") > 0:
            risk_factors_count += 1
            findings.append("COL1A1 Sp1 variant: altered collagen structure; increased fracture risk")

        # Postmenopausal female + ESR1 variant
        if sex.lower() == "female" and age is not None and age > 50:
            esr1 = genotypes.get("ESR1_rs2234693", "")
            if _count_risk_alleles(esr1, "T") > 0:
                risk_factors_count += 1
                findings.append("ESR1 variant in postmenopausal female: increased bone loss risk")
                recommendations.append("Consider DEXA scan if not done in last 2 years")

        # Determine stage
        if risk_factors_count >= 3:
            stage = "fracture_risk"
            risk_level = "HIGH"
            recommendations.append("DEXA scan recommended; evaluate for pharmacological intervention")
        elif risk_factors_count >= 2:
            stage = "osteoporosis_risk"
            risk_level = "MODERATE"
        elif risk_factors_count >= 1:
            stage = "osteopenia_risk"
            risk_level = "MODERATE"

        return {
            "disease": "bone_health",
            "display_name": "Osteoporosis / Bone Health",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "genetic_risk_factors": genetic_risk_factors,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_cognitive(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
    ) -> Dict[str, Any]:
        """Analyze cognitive decline and Alzheimer's risk trajectory.

        Evaluates APOE genotype with modifiable risk factors: homocysteine,
        B12, folate, omega-3, inflammation, and vitamin D.

        Args:
            biomarkers: Dict with homocysteine, vitamin_b12, folate,
                omega3_index, hs_crp, vitamin_d_25oh.
            genotypes: Dict with APOE, MTHFR_rs1801133, BDNF_rs6265.

        Returns:
            Dict with disease trajectory analysis.
        """
        available = {k: v for k, v in biomarkers.items()
                     if k in DISEASE_CONFIGS["cognitive"]["biomarkers"]}

        # APOE genotype — strongest genetic risk factor for Alzheimer's
        apoe = genotypes.get("APOE", "")
        has_e4 = "E4" in apoe
        e4_count = apoe.count("E4")

        genetic_risk_factors = []
        if has_e4:
            genetic_risk_factors.append({
                "gene": "APOE", "genotype": apoe,
                "risk_alleles": e4_count, "effect": "amyloid_clearance",
            })
        for gene_key in ("MTHFR_rs1801133", "BDNF_rs6265"):
            gt = genotypes.get(gene_key, "")
            if gt:
                config = DISEASE_CONFIGS["cognitive"]["genetic_modifiers"][gene_key]
                count = _count_risk_alleles(gt, config["risk_allele"])
                if count > 0:
                    genetic_risk_factors.append({
                        "gene": gene_key, "genotype": gt,
                        "risk_alleles": count, "effect": config["effect"],
                    })

        stage = "low_risk"
        risk_level = "LOW"
        findings = []
        recommendations = []
        modifiable_risk_count = 0

        # APOE E4 — non-modifiable risk
        if e4_count >= 2:
            risk_level = "HIGH"
            findings.append(
                "APOE E4/E4 homozygote: ~12x increased Alzheimer's risk. "
                "Aggressive modifiable risk factor optimization critical"
            )
        elif e4_count == 1:
            risk_level = _max_risk(risk_level, "MODERATE")
            findings.append(
                f"APOE {apoe}: ~3x increased Alzheimer's risk. "
                "Focus on modifiable protective factors"
            )

        # Homocysteine — neurotoxic at elevated levels
        hcy = available.get("homocysteine")
        if hcy is not None:
            if hcy > 15:
                modifiable_risk_count += 1
                findings.append(f"Homocysteine {hcy} umol/L elevated (>15) — neurotoxic, accelerates brain atrophy")
                mthfr = genotypes.get("MTHFR_rs1801133", "")
                mthfr_risk = _count_risk_alleles(mthfr, "T")
                if mthfr_risk >= 2:
                    recommendations.append(
                        "MTHFR TT + elevated homocysteine: L-methylfolate 1-5mg + methylB12 1000mcg + P5P 50mg daily"
                    )
                else:
                    recommendations.append("Optimize B-vitamins: folate + B12 + B6 to reduce homocysteine")
            elif hcy > 10 and has_e4:
                modifiable_risk_count += 1
                findings.append(f"Homocysteine {hcy} umol/L is suboptimal for APOE E4 carrier (target <10)")

        # Vitamin B12
        b12 = available.get("vitamin_b12")
        if b12 is not None and b12 < 400:
            modifiable_risk_count += 1
            findings.append(f"Vitamin B12 {b12} pg/mL suboptimal for neuroprotection (<400)")
            recommendations.append("Methylcobalamin 1000-2000 mcg/day for neuroprotection")

        # Omega-3 index — DHA critical for brain structure
        omega3 = available.get("omega3_index")
        if omega3 is not None and omega3 < 8:
            modifiable_risk_count += 1
            findings.append(f"Omega-3 index {omega3}% below neuroprotective threshold (<8%)")
            recommendations.append("EPA/DHA supplementation 2-4g/day; prioritize DHA for brain health")

        # hs-CRP — neuroinflammation marker
        crp = available.get("hs_crp")
        if crp is not None and crp > 2.0:
            modifiable_risk_count += 1
            findings.append(f"hs-CRP {crp} mg/L indicates inflammation — neuroinflammation accelerates neurodegeneration")
            recommendations.append("Anti-inflammatory lifestyle: Mediterranean diet, regular exercise, stress management")

        # Vitamin D
        vit_d = available.get("vitamin_d_25oh")
        if vit_d is not None and vit_d < 30:
            modifiable_risk_count += 1
            findings.append(f"Vitamin D {vit_d} ng/mL insufficient (<30) — associated with cognitive decline")
            recommendations.append("Vitamin D3 2000-5000 IU/day; target >40 ng/mL for neuroprotection")

        # BDNF Val66Met
        bdnf = genotypes.get("BDNF_rs6265", "")
        if _count_risk_alleles(bdnf, "A") > 0:
            findings.append("BDNF Val66Met variant: reduced activity-dependent BDNF secretion")
            recommendations.append("Regular aerobic exercise (150+ min/week) — strongest BDNF inducer")

        # Stage determination
        if has_e4 and modifiable_risk_count >= 2:
            stage = "high_risk"
            risk_level = _max_risk(risk_level, "HIGH")
            recommendations.append("Consider comprehensive cognitive assessment and longitudinal monitoring")
        elif has_e4 or modifiable_risk_count >= 3:
            stage = "elevated_risk"
            risk_level = _max_risk(risk_level, "MODERATE")
        elif modifiable_risk_count >= 1:
            stage = "modifiable_risk_factors"

        return {
            "disease": "cognitive",
            "display_name": "Cognitive Decline / Alzheimer's Risk",
            "risk_level": risk_level,
            "stage": stage,
            "current_markers": available,
            "genetic_risk_factors": genetic_risk_factors,
            "modifiable_risk_factors": modifiable_risk_count,
            "findings": findings,
            "recommendations": recommendations,
        }

    def analyze_all(
        self,
        biomarkers: Dict[str, float],
        genotypes: Dict[str, str],
        age: Optional[float] = None,
        sex: str = "male",
    ) -> List[Dict[str, Any]]:
        """Run all 9 disease trajectory analyses.

        Args:
            biomarkers: Combined dict of all available biomarkers.
            genotypes: Combined dict of all available genotypes.
            age: Patient age in years (needed for FIB-4 in liver analysis).
            sex: Patient sex for sex-stratified reference ranges.

        Returns:
            List of disease trajectory result dicts, sorted by risk level.
        """
        # Sanitize biomarker values — filter out NaN and None
        clean_biomarkers = {
            k: v for k, v in biomarkers.items()
            if v is not None and not (isinstance(v, float) and v != v)  # v != v detects NaN
        }

        results = []

        results.append(self.analyze_type2_diabetes(clean_biomarkers, genotypes))
        results.append(self.analyze_cardiovascular(clean_biomarkers, genotypes))
        results.append(self.analyze_liver(clean_biomarkers, genotypes, age=age))
        results.append(self.analyze_thyroid(clean_biomarkers, genotypes))
        results.append(self.analyze_iron(clean_biomarkers, genotypes, sex=sex))
        results.append(self.analyze_nutritional(clean_biomarkers, genotypes))
        results.append(self.analyze_kidney(clean_biomarkers, genotypes, age=age, sex=sex))
        results.append(self.analyze_bone_health(clean_biomarkers, genotypes, sex=sex, age=age))
        results.append(self.analyze_cognitive(clean_biomarkers, genotypes))

        # Sort by risk level severity
        risk_order = {"CRITICAL": 0, "HIGH": 1, "MODERATE": 2, "LOW": 3}
        results.sort(key=lambda r: risk_order.get(r["risk_level"], 4))

        logger.info(
            f"Disease trajectory analysis complete: "
            f"{sum(1 for r in results if r['risk_level'] in ('HIGH', 'CRITICAL'))} high-risk, "
            f"{sum(1 for r in results if r['risk_level'] == 'MODERATE')} moderate-risk"
        )

        return results
