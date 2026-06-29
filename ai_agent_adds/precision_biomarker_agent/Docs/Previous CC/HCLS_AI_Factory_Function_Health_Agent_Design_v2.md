# HCLS AI Factory: Function Health Biomarker-Genomic Correlation Agent

## Advanced Design Document - Version 2.0

**Author:** Adam Jones
**Date:** January 2026
**Version:** 2.0 - Enhanced Edition

---

## Executive Summary

This document presents a next-generation AI Agent that transforms comprehensive blood biomarker data from Function Health into **life-changing personalized health insights** through deep integration with genomic variant analysis. This agent goes far beyond simple reference range checking to provide:

- **Biological Age Calculation** using blood biomarkers correlated with epigenetic clocks
- **Pre-Symptomatic Disease Detection** - identifying disease trajectories 5-10 years before clinical diagnosis
- **Genotype-Adjusted Optimal Ranges** - personalized biomarker thresholds based on individual genetics
- **Disease Stage Assessment** - determining where someone falls on disease progression continuum
- **Biomarker Trajectory Modeling** - longitudinal analysis predicting future health states
- **Pharmacogenomic Optimization** - drug-gene interactions with dosing guidance
- **Multi-Pathway Risk Integration** - connecting inflammation, metabolism, cardiovascular, and cognitive pathways

**The Vision:** A person uploads their Function Health results and receives insights that would typically require a team of specialists: a cardiologist, endocrinologist, geneticist, pharmacologist, and longevity researcher - all synthesized into actionable intelligence.

---

## Part I: Revolutionary Capabilities

### 1. Biological Age Engine

The agent calculates biological age using biomarkers that correlate with validated epigenetic clocks (GrimAge, PhenoAge, DunedinPACE).

#### Blood Biomarkers That Predict Biological Age

| Biomarker | Epigenetic Clock Correlation | What It Reveals |
|-----------|------------------------------|-----------------|
| **GDF-15** | Strong correlation with GrimAge (r=0.664) | Mitochondrial stress, cellular senescence, all-cause mortality |
| **Cystatin C** | GrimAge component | Kidney function beyond creatinine, CV mortality |
| **hs-CRP** | PhenoAge component | Systemic inflammation, vascular aging |
| **Albumin** | PhenoAge component (inverse) | Liver synthetic function, nutritional status |
| **Lymphocyte %** | PhenoAge component | Immune aging (immunosenescence) |
| **Red Cell Distribution Width (RDW)** | PhenoAge component | Bone marrow aging, chronic disease burden |
| **Mean Cell Volume (MCV)** | PhenoAge component | Nutritional deficiencies, bone marrow health |
| **White Blood Cell Count** | PhenoAge component | Chronic inflammation, immune activation |
| **Alkaline Phosphatase** | PhenoAge component | Liver/bone health, vascular calcification |

```python
class BiologicalAgeCalculator:
    """
    Estimates biological age from blood biomarkers using validated
    correlations with epigenetic clocks (GrimAge, PhenoAge)
    """

    def __init__(self):
        # PhenoAge biomarker weights (Levine et al. 2018)
        self.phenoage_coefficients = {
            "albumin": -0.0336,        # g/L (protective)
            "creatinine": 0.0095,       # µmol/L
            "glucose": 0.1953,          # mmol/L
            "hs_crp": 0.0954,           # ln(mg/L)
            "lymphocyte_pct": -0.0120,  # % (protective)
            "mcv": 0.0268,              # fL
            "rdw": 0.3306,              # %
            "alkaline_phosphatase": 0.0019,  # U/L
            "wbc": 0.0554,              # 10^3 cells/µL
            "chronological_age": 0.0804
        }

        # GrimAge surrogate markers
        self.grimage_markers = {
            "gdf15": {"weight": 0.15, "unit": "pg/mL"},      # Growth differentiation factor 15
            "cystatin_c": {"weight": 0.12, "unit": "mg/L"},  # Kidney function
            "leptin": {"weight": 0.08, "unit": "ng/mL"},     # Adiposity
            "pai1": {"weight": 0.10, "unit": "ng/mL"},       # Plasminogen activator inhibitor
            "timp1": {"weight": 0.09, "unit": "ng/mL"},      # Tissue inhibitor metalloproteinase
            "adm": {"weight": 0.11, "unit": "pmol/L"},       # Adrenomedullin
            "b2m": {"weight": 0.10, "unit": "mg/L"},         # Beta-2-microglobulin
        }

    def calculate_phenoage(self, biomarkers: dict, chronological_age: int) -> dict:
        """
        Calculate PhenoAge and age acceleration

        Returns:
            biological_age: Estimated biological age
            age_acceleration: Difference from chronological age
            risk_interpretation: What this means for health
        """
        # Transform biomarkers
        transformed = {
            "albumin": biomarkers.get("albumin", 45),
            "creatinine": biomarkers.get("creatinine", 80),
            "glucose": biomarkers.get("glucose_fasting", 90) / 18,  # mg/dL to mmol/L
            "hs_crp": np.log(biomarkers.get("hs_crp", 1.0) + 0.001),
            "lymphocyte_pct": biomarkers.get("lymphocyte_pct", 30),
            "mcv": biomarkers.get("mcv", 90),
            "rdw": biomarkers.get("rdw", 13),
            "alkaline_phosphatase": biomarkers.get("alkaline_phosphatase", 70),
            "wbc": biomarkers.get("wbc", 6.0),
        }

        # Calculate mortality score
        mortality_score = sum(
            self.phenoage_coefficients[k] * transformed[k]
            for k in transformed
        )
        mortality_score += self.phenoage_coefficients["chronological_age"] * chronological_age

        # Convert to biological age (inverse of mortality equation)
        biological_age = self._mortality_to_age(mortality_score)
        age_acceleration = biological_age - chronological_age

        return {
            "biological_age": round(biological_age, 1),
            "chronological_age": chronological_age,
            "age_acceleration": round(age_acceleration, 1),
            "interpretation": self._interpret_age_acceleration(age_acceleration),
            "contributing_factors": self._identify_aging_drivers(transformed, biomarkers)
        }

    def _interpret_age_acceleration(self, acceleration: float) -> dict:
        """Interpret what age acceleration means"""
        if acceleration <= -5:
            return {
                "status": "EXCEPTIONAL",
                "description": "Biological age significantly younger than chronological age",
                "mortality_risk": "Substantially reduced all-cause mortality risk",
                "implications": "Current lifestyle/interventions highly effective"
            }
        elif acceleration <= -2:
            return {
                "status": "EXCELLENT",
                "description": "Aging slower than average",
                "mortality_risk": "Reduced mortality risk",
                "implications": "Health optimization working well"
            }
        elif acceleration <= 2:
            return {
                "status": "NORMAL",
                "description": "Aging at expected rate",
                "mortality_risk": "Average population risk",
                "implications": "Room for optimization"
            }
        elif acceleration <= 5:
            return {
                "status": "ACCELERATED",
                "description": "Aging faster than expected",
                "mortality_risk": "Elevated mortality and disease risk",
                "implications": "Intervention recommended"
            }
        else:
            return {
                "status": "SIGNIFICANTLY ACCELERATED",
                "description": "Biological age much older than chronological",
                "mortality_risk": "Substantially elevated mortality risk",
                "implications": "Urgent intervention and medical evaluation needed"
            }

    def _identify_aging_drivers(self, transformed: dict, raw: dict) -> list:
        """Identify which biomarkers are driving age acceleration"""
        drivers = []

        # Check inflammation
        if raw.get("hs_crp", 1) > 2.0:
            drivers.append({
                "factor": "Chronic Inflammation",
                "biomarker": "hs-CRP",
                "value": raw["hs_crp"],
                "impact": "HIGH",
                "intervention": "Anti-inflammatory diet, omega-3s, address root cause"
            })

        # Check metabolic
        if raw.get("glucose_fasting", 90) > 100:
            drivers.append({
                "factor": "Metabolic Dysfunction",
                "biomarker": "Fasting Glucose",
                "value": raw["glucose_fasting"],
                "impact": "HIGH",
                "intervention": "Carbohydrate reduction, exercise, metabolic optimization"
            })

        # Check immune aging
        if raw.get("lymphocyte_pct", 30) < 20:
            drivers.append({
                "factor": "Immunosenescence",
                "biomarker": "Lymphocyte %",
                "value": raw["lymphocyte_pct"],
                "impact": "MODERATE",
                "intervention": "Exercise, sleep optimization, stress reduction"
            })

        # Check RDW (bone marrow/chronic disease)
        if raw.get("rdw", 13) > 14.5:
            drivers.append({
                "factor": "Chronic Disease Burden",
                "biomarker": "RDW",
                "value": raw["rdw"],
                "impact": "MODERATE",
                "intervention": "Investigate underlying conditions, nutritional optimization"
            })

        return drivers
```

---

### 2. Pre-Symptomatic Disease Detection Engine

The agent identifies disease trajectories **years before clinical diagnosis** by correlating biomarker patterns with genetic risk.

#### Type 2 Diabetes: Detectable 5-7 Years Before Diagnosis

```python
class PreDiabeticTrajectoryAnalyzer:
    """
    Identifies pre-diabetic trajectory using biomarkers + genetics
    Studies show diabetes can be predicted 5-7 years before diagnosis
    """

    def __init__(self):
        # Key diabetes susceptibility genes
        self.diabetes_genes = {
            "TCF7L2": {
                "rs7903146": {
                    "risk_allele": "T",
                    "effect": "1.4x risk per allele - most important T2D gene",
                    "mechanism": "Impaired beta cell function and insulin secretion"
                }
            },
            "PPARG": {
                "rs1801282": {
                    "risk_allele": "C",
                    "effect": "Protective Pro12Ala reduces risk 20%",
                    "mechanism": "Improved insulin sensitivity"
                }
            },
            "SLC30A8": {
                "rs13266634": {
                    "risk_allele": "C",
                    "effect": "15% increased risk per allele",
                    "mechanism": "Zinc transporter affecting insulin secretion"
                }
            },
            "KCNJ11": {
                "rs5219": {
                    "risk_allele": "T",
                    "effect": "15% increased risk per allele",
                    "mechanism": "ATP-sensitive potassium channel in beta cells"
                }
            },
            "GCKR": {
                "rs780094": {
                    "risk_allele": "T",
                    "effect": "Increases fasting glucose and triglycerides",
                    "mechanism": "Glucokinase regulatory protein"
                }
            }
        }

    def calculate_diabetes_trajectory(
        self,
        biomarkers: dict,
        genetic_variants: dict,
        age: int
    ) -> dict:
        """
        Calculate personalized diabetes trajectory

        Key predictive biomarkers:
        - Fasting glucose trajectory
        - Fasting insulin / HOMA-IR
        - HbA1c
        - Triglyceride/HDL ratio
        - 1-hour post-glucose (if available)
        """
        result = {
            "current_state": "",
            "genetic_risk": {},
            "trajectory_prediction": {},
            "years_to_diagnosis": None,
            "interventions": [],
            "reversal_potential": ""
        }

        # Calculate HOMA-IR
        glucose = biomarkers.get("glucose_fasting", 90)
        insulin = biomarkers.get("insulin_fasting", 8)
        homa_ir = (glucose * insulin) / 405

        hba1c = biomarkers.get("hba1c", 5.4)
        tg_hdl_ratio = biomarkers.get("triglycerides", 100) / biomarkers.get("hdl_c", 50)

        # Genetic risk score
        genetic_risk_score = self._calculate_genetic_risk(genetic_variants)
        result["genetic_risk"] = genetic_risk_score

        # Current metabolic state
        if hba1c >= 6.5 or glucose >= 126:
            result["current_state"] = "DIABETES"
            result["reversal_potential"] = "May be reversible with aggressive intervention if recent onset"
        elif hba1c >= 5.7 or glucose >= 100 or homa_ir > 2.5:
            result["current_state"] = "PREDIABETES"
            result["reversal_potential"] = "HIGH - 58% risk reduction possible with lifestyle intervention"
            result["years_to_diagnosis"] = self._estimate_years_to_diabetes(
                hba1c, glucose, homa_ir, genetic_risk_score["overall_risk"]
            )
        elif homa_ir > 1.5 or tg_hdl_ratio > 2.5:
            result["current_state"] = "EARLY INSULIN RESISTANCE"
            result["reversal_potential"] = "VERY HIGH - Early intervention highly effective"
            result["trajectory_prediction"] = {
                "without_intervention": "Progression to prediabetes likely within 3-5 years",
                "with_intervention": "Fully reversible with lifestyle modification"
            }
        else:
            result["current_state"] = "METABOLICALLY HEALTHY"

        # Genotype-specific recommendations
        if "TCF7L2" in genetic_variants:
            tcf7l2 = genetic_variants["TCF7L2"].get("rs7903146", {})
            if tcf7l2.get("genotype") in ["TT", "CT"]:
                result["interventions"].append({
                    "priority": "HIGH",
                    "intervention": "Carbohydrate quality focus",
                    "rationale": "TCF7L2 risk carriers show greater benefit from reduced glycemic load",
                    "evidence": "DPP study showed TCF7L2 carriers benefit equally from lifestyle intervention"
                })

        # HOMA-IR specific interventions
        if homa_ir > 2.0:
            result["interventions"].append({
                "priority": "HIGH",
                "intervention": "Time-restricted eating (16:8)",
                "rationale": "Improves insulin sensitivity independent of weight loss",
                "expected_benefit": "20-30% improvement in HOMA-IR within 12 weeks"
            })

        return result

    def _calculate_genetic_risk(self, variants: dict) -> dict:
        """Calculate composite genetic risk score for T2D"""
        risk_allele_count = 0
        total_possible = 0
        variant_contributions = []

        for gene, snps in self.diabetes_genes.items():
            for rsid, info in snps.items():
                total_possible += 2
                if gene in variants and rsid in variants[gene]:
                    genotype = variants[gene][rsid].get("genotype", "")
                    count = genotype.count(info["risk_allele"])
                    risk_allele_count += count
                    if count > 0:
                        variant_contributions.append({
                            "gene": gene,
                            "rsid": rsid,
                            "risk_alleles": count,
                            "effect": info["effect"]
                        })

        genetic_risk_percentile = (risk_allele_count / total_possible) * 100

        return {
            "risk_allele_count": risk_allele_count,
            "percentile": genetic_risk_percentile,
            "overall_risk": "HIGH" if genetic_risk_percentile > 60 else "MODERATE" if genetic_risk_percentile > 40 else "LOW",
            "contributing_variants": variant_contributions
        }

    def _estimate_years_to_diabetes(
        self, hba1c: float, glucose: float, homa_ir: float, genetic_risk: str
    ) -> dict:
        """
        Estimate years until diabetes diagnosis without intervention
        Based on trajectory modeling from DPP and other cohort studies
        """
        # Base trajectory from biomarkers
        if hba1c >= 6.0:
            base_years = 2
        elif hba1c >= 5.7:
            base_years = 4
        else:
            base_years = 7

        # Modify by HOMA-IR
        if homa_ir > 3.0:
            base_years -= 1
        elif homa_ir > 2.5:
            base_years -= 0.5

        # Modify by genetic risk
        if genetic_risk == "HIGH":
            base_years *= 0.75
        elif genetic_risk == "LOW":
            base_years *= 1.25

        return {
            "without_intervention": f"{max(1, int(base_years))}-{max(2, int(base_years + 2))} years",
            "with_intervention": f"Risk reduced by 58% - may never develop diabetes",
            "key_metric_to_track": "HbA1c every 6 months, HOMA-IR annually"
        }
```

---

### 3. Cardiovascular Risk: Beyond Cholesterol

The agent provides sophisticated cardiovascular risk assessment integrating emerging biomarkers with genetics.

```python
class AdvancedCardiovascularAnalyzer:
    """
    Multi-pathway cardiovascular risk assessment integrating:
    - Lipid genetics (PCSK9, LPA, APOE)
    - Inflammation pathway (IL-6, hs-CRP, genetics)
    - Metabolic pathway
    - Vascular aging markers
    """

    def __init__(self):
        # Key cardiovascular genes
        self.cv_genes = {
            "LPA": {
                "description": "Lipoprotein(a) - stable from childhood, 90% genetic",
                "clinical_note": "Lp(a) >50 nmol/L = significantly elevated CV risk"
            },
            "PCSK9": {
                "rs11591147": {
                    "effect": "Loss-of-function: 88% CHD reduction per 37 mg/dL LDL decrease",
                    "clinical": "Predicts exceptional statin response"
                }
            },
            "IL6": {
                "description": "Variants mimicking IL-6 inhibition reduce CV risk",
                "clinical": "IL-6 more predictive than hs-CRP for CV events"
            },
            "CRP": {
                "rs1800947": "Associated with higher CRP and mortality in heart failure",
                "rs11265263": "Associated with higher CRP and CV mortality"
            }
        }

        # Inflammation pathway
        self.inflammation_thresholds = {
            "hs_crp": {"low": 1.0, "moderate": 3.0, "high": 10.0},
            "il6": {"low": 2.0, "moderate": 5.0, "high": 10.0},
            "nlr": {"low": 2.0, "moderate": 3.0, "high": 5.0},  # Neutrophil-lymphocyte ratio
            "gdf15": {"low": 1000, "moderate": 1500, "high": 2000}  # pg/mL
        }

    def comprehensive_cv_assessment(
        self,
        biomarkers: dict,
        genetic_variants: dict,
        age: int,
        sex: str
    ) -> dict:
        """
        Comprehensive cardiovascular risk assessment
        """
        result = {
            "traditional_risk": {},
            "advanced_lipid_risk": {},
            "inflammation_risk": {},
            "genetic_risk": {},
            "vascular_aging": {},
            "integrated_risk": "",
            "personalized_interventions": []
        }

        # === ADVANCED LIPID ASSESSMENT ===
        apoB = biomarkers.get("apoB", 90)
        lp_a = biomarkers.get("lp_a", 30)
        ldl_p = biomarkers.get("ldl_particle_number", 1000)

        result["advanced_lipid_risk"] = {
            "apoB": {
                "value": apoB,
                "interpretation": "ELEVATED" if apoB > 90 else "OPTIMAL" if apoB < 80 else "BORDERLINE",
                "note": "ApoB is the best measure of atherogenic particle burden"
            },
            "lp_a": {
                "value": lp_a,
                "interpretation": self._interpret_lpa(lp_a),
                "genetic_note": "Lp(a) is 90% genetically determined - stable from childhood",
                "action": self._lpa_action(lp_a)
            }
        }

        # === INFLAMMATION PATHWAY ===
        hs_crp = biomarkers.get("hs_crp", 1.0)
        neutrophils = biomarkers.get("neutrophils_absolute", 4.0)
        lymphocytes = biomarkers.get("lymphocytes_absolute", 2.0)
        nlr = neutrophils / lymphocytes if lymphocytes > 0 else 0

        result["inflammation_risk"] = {
            "hs_crp": {
                "value": hs_crp,
                "risk_category": "LOW" if hs_crp < 1 else "MODERATE" if hs_crp < 3 else "HIGH",
                "note": "hs-CRP is downstream marker; IL-6 is upstream driver"
            },
            "nlr": {
                "value": round(nlr, 2),
                "interpretation": self._interpret_nlr(nlr),
                "note": "NLR from routine CBC predicts CV events independently"
            },
            "inflammatory_phenotype": self._determine_inflammatory_phenotype(hs_crp, nlr)
        }

        # === GENETIC CARDIOVASCULAR RISK ===
        result["genetic_risk"] = self._assess_cv_genetic_risk(genetic_variants)

        # === VASCULAR AGING ===
        gdf15 = biomarkers.get("gdf15", None)
        if gdf15:
            result["vascular_aging"] = {
                "gdf15": {
                    "value": gdf15,
                    "interpretation": "Normal" if gdf15 < 1000 else "Elevated" if gdf15 < 1500 else "High",
                    "note": "GDF-15 correlates with arterial stiffness and predicts all-cause mortality"
                }
            }

        # === INTEGRATED RISK AND INTERVENTIONS ===
        result["integrated_risk"] = self._integrate_risk(result)
        result["personalized_interventions"] = self._generate_cv_interventions(result, genetic_variants)

        return result

    def _interpret_lpa(self, lp_a: float) -> dict:
        """Interpret Lp(a) levels"""
        if lp_a < 30:
            return {"risk": "LOW", "percentile": "Below 50th"}
        elif lp_a < 50:
            return {"risk": "BORDERLINE", "percentile": "50th-75th"}
        elif lp_a < 125:
            return {"risk": "ELEVATED", "percentile": "75th-90th"}
        else:
            return {"risk": "VERY HIGH", "percentile": "Above 90th"}

    def _lpa_action(self, lp_a: float) -> dict:
        """Action recommendations for Lp(a)"""
        if lp_a < 50:
            return {"action": "No specific intervention needed", "monitoring": "Recheck in 5 years"}
        elif lp_a < 125:
            return {
                "action": "Aggressive LDL-C lowering recommended",
                "target": "LDL-C <70 mg/dL to offset Lp(a) risk",
                "note": "Lp(a) adds ~6% CV risk per year at this level",
                "emerging": "PCSK9 inhibitors reduce Lp(a) 20-30%"
            }
        else:
            return {
                "action": "Very aggressive risk factor modification",
                "target": "LDL-C <55 mg/dL",
                "note": "Consider cardiology referral",
                "emerging": "Novel Lp(a)-lowering therapies in development (pelacarsen)"
            }

    def _interpret_nlr(self, nlr: float) -> str:
        """Interpret neutrophil-to-lymphocyte ratio"""
        if nlr < 2.0:
            return "LOW INFLAMMATION - Favorable cardiovascular profile"
        elif nlr < 3.0:
            return "NORMAL - Average inflammatory state"
        elif nlr < 5.0:
            return "ELEVATED - Subclinical inflammation, increased CV risk"
        else:
            return "HIGH - Significant inflammatory state, evaluate for underlying conditions"

    def _determine_inflammatory_phenotype(self, hs_crp: float, nlr: float) -> dict:
        """Determine inflammatory phenotype for targeted intervention"""
        if hs_crp < 1.0 and nlr < 2.0:
            return {
                "phenotype": "LOW INFLAMMATION",
                "intervention_priority": "Low",
                "note": "Focus on other risk factors"
            }
        elif hs_crp > 3.0 and nlr > 3.0:
            return {
                "phenotype": "SYSTEMIC INFLAMMATION",
                "intervention_priority": "High",
                "note": "Both hepatic (CRP) and cellular (NLR) inflammation elevated",
                "action": "Investigate sources: visceral adiposity, infection, autoimmune, gut"
            }
        elif hs_crp > 3.0 and nlr < 2.0:
            return {
                "phenotype": "HEPATIC/METABOLIC INFLAMMATION",
                "intervention_priority": "Moderate",
                "note": "Liver-mediated inflammation without immune activation",
                "action": "Focus on metabolic health, weight loss if indicated"
            }
        else:
            return {
                "phenotype": "IMMUNE-MEDIATED INFLAMMATION",
                "intervention_priority": "Moderate",
                "note": "Cellular inflammation without hepatic acute phase response",
                "action": "Evaluate immune function, consider autoimmune screening"
            }

    def _assess_cv_genetic_risk(self, variants: dict) -> dict:
        """Assess genetic cardiovascular risk"""
        risk_factors = []

        # Check APOE
        apoe = variants.get("APOE", {}).get("genotype", "E3/E3")
        if "E4" in apoe:
            risk_factors.append({
                "gene": "APOE",
                "genotype": apoe,
                "effect": "Elevated LDL-C, enhanced statin response",
                "action": "Aggressive LDL management, statins highly effective"
            })

        # Check PCSK9
        if "PCSK9" in variants:
            pcsk9_lof = variants["PCSK9"].get("rs11591147", {})
            if pcsk9_lof.get("genotype", "").count("T") > 0:
                risk_factors.append({
                    "gene": "PCSK9",
                    "variant": "rs11591147 (R46L)",
                    "effect": "Loss-of-function: 15-28% lower LDL, 47-88% lower CHD risk",
                    "note": "Lifelong protection from elevated LDL"
                })

        # Check LPA for genetic context
        # (Lp(a) is 90% genetic, so elevated Lp(a) = genetic risk)

        return {
            "identified_variants": risk_factors,
            "overall_genetic_cv_risk": "ELEVATED" if len(risk_factors) > 0 else "AVERAGE"
        }

    def _generate_cv_interventions(self, risk_result: dict, variants: dict) -> list:
        """Generate personalized cardiovascular interventions"""
        interventions = []

        # High Lp(a)
        if risk_result["advanced_lipid_risk"]["lp_a"]["interpretation"]["risk"] in ["ELEVATED", "VERY HIGH"]:
            interventions.append({
                "target": "Elevated Lp(a)",
                "priority": "HIGH",
                "interventions": [
                    "Aggressive LDL-C lowering (target <70, consider <55 mg/dL)",
                    "PCSK9 inhibitor consideration (reduces Lp(a) 20-30%)",
                    "Niacin may reduce Lp(a) 20-30% (discuss with physician)",
                    "Aspirin for primary prevention may be beneficial (consult physician)"
                ],
                "monitoring": "Annual Lp(a) not necessary (genetically fixed)",
                "emerging": "Pelacarsen (antisense oligonucleotide) in Phase 3 trials"
            })

        # Inflammation
        if risk_result["inflammation_risk"]["hs_crp"]["risk_category"] == "HIGH":
            interventions.append({
                "target": "Elevated Inflammation",
                "priority": "HIGH",
                "interventions": [
                    "Mediterranean diet (shown to reduce hs-CRP 30-40%)",
                    "Omega-3 fatty acids 2-4g EPA+DHA daily",
                    "Regular aerobic exercise (150 min/week)",
                    "Weight loss if BMI >25 (each kg lost reduces CRP 0.13 mg/L)",
                    "Consider colchicine 0.5mg daily (COLCOT trial: 23% CV reduction)"
                ],
                "investigate": "If CRP >10, evaluate for infection or inflammatory condition",
                "monitoring": "Recheck hs-CRP in 3 months after interventions"
            })

        # APOE4 specific
        apoe = variants.get("APOE", {}).get("genotype", "")
        if "E4" in apoe:
            interventions.append({
                "target": "APOE4 Carrier",
                "priority": "MODERATE",
                "interventions": [
                    "Statins highly effective - enhanced LDL response",
                    "Mediterranean diet particularly beneficial",
                    "Limit saturated fat (<7% calories) - E4 carriers more sensitive",
                    "Omega-3s may provide cognitive and CV protection",
                    "Regular aerobic exercise (neuroprotective and cardioprotective)"
                ],
                "note": "APOE4 is both CV and Alzheimer's risk factor - address both"
            })

        return interventions
```

---

### 4. Liver Health Genetics: Silent Disease Detection

```python
class LiverHealthAnalyzer:
    """
    Integrates liver enzymes with genetic risk for fatty liver disease (NAFLD/MASH)

    Key insight: Liver enzymes within "normal" range can still indicate disease
    in genetically susceptible individuals
    """

    def __init__(self):
        # NAFLD susceptibility genes
        self.nafld_genes = {
            "PNPLA3": {
                "rs738409": {
                    "risk_allele": "G",
                    "effect": "Major NAFLD/NASH risk gene - 2-3x risk per allele",
                    "mechanism": "Impaired lipid droplet metabolism",
                    "note": "Associated with fibrosis, not just steatosis"
                }
            },
            "TM6SF2": {
                "rs58542926": {
                    "risk_allele": "T",
                    "effect": "Increased liver fat but REDUCED cardiovascular risk",
                    "mechanism": "Impaired VLDL secretion - keeps fat in liver",
                    "note": "Paradoxical protection from heart disease"
                }
            },
            "HSD17B13": {
                "rs72613567": {
                    "effect": "PROTECTIVE - loss-of-function reduces progression",
                    "note": "Drug development target (HSD17B13 inhibitors)"
                }
            }
        }

    def analyze_liver_health(
        self,
        biomarkers: dict,
        genetic_variants: dict
    ) -> dict:
        """
        Comprehensive liver health analysis integrating enzymes and genetics
        """
        alt = biomarkers.get("alt", 25)
        ast = biomarkers.get("ast", 25)
        ggt = biomarkers.get("ggt", 30)
        ferritin = biomarkers.get("ferritin", 100)
        platelets = biomarkers.get("platelets", 250)
        albumin = biomarkers.get("albumin", 4.5)

        result = {
            "enzyme_analysis": {},
            "genetic_risk": {},
            "fibrosis_assessment": {},
            "integrated_assessment": "",
            "recommendations": []
        }

        # Genetic risk assessment
        pnpla3_risk = self._assess_pnpla3(genetic_variants)
        tm6sf2_risk = self._assess_tm6sf2(genetic_variants)
        result["genetic_risk"] = {
            "PNPLA3": pnpla3_risk,
            "TM6SF2": tm6sf2_risk,
            "overall_genetic_susceptibility": "HIGH" if pnpla3_risk["risk"] == "HIGH" else "MODERATE" if pnpla3_risk["risk"] == "MODERATE" else "LOW"
        }

        # Genotype-adjusted enzyme interpretation
        # Key insight: "Normal" ALT may still be elevated in context of genetics
        result["enzyme_analysis"] = {
            "alt": {
                "value": alt,
                "standard_interpretation": "Normal" if alt < 40 else "Elevated",
                "genotype_adjusted": self._genotype_adjusted_alt(alt, pnpla3_risk["risk"])
            },
            "ast_alt_ratio": {
                "value": round(ast / alt, 2) if alt > 0 else 0,
                "interpretation": "Suggests fibrosis" if ast / alt > 1 else "No fibrosis signal"
            },
            "ggt": {
                "value": ggt,
                "interpretation": "Suggests biliary/alcohol component" if ggt > 50 else "Normal"
            }
        }

        # FIB-4 fibrosis score
        age = biomarkers.get("age", 45)
        fib4 = (age * ast) / (platelets * (alt ** 0.5))
        result["fibrosis_assessment"] = {
            "fib4_score": round(fib4, 2),
            "interpretation": self._interpret_fib4(fib4),
            "genetic_context": self._fib4_genetic_context(fib4, pnpla3_risk["risk"])
        }

        # Integrated assessment
        result["integrated_assessment"] = self._integrated_liver_assessment(result)
        result["recommendations"] = self._liver_recommendations(result, biomarkers)

        return result

    def _assess_pnpla3(self, variants: dict) -> dict:
        """Assess PNPLA3 rs738409 risk"""
        if "PNPLA3" in variants:
            genotype = variants["PNPLA3"].get("rs738409", {}).get("genotype", "CC")
            if genotype == "GG":
                return {"risk": "HIGH", "genotype": "GG", "effect": "3x increased NASH/fibrosis risk"}
            elif genotype == "CG":
                return {"risk": "MODERATE", "genotype": "CG", "effect": "1.7x increased risk"}
        return {"risk": "LOW", "genotype": "CC", "effect": "Average population risk"}

    def _assess_tm6sf2(self, variants: dict) -> dict:
        """Assess TM6SF2 rs58542926"""
        if "TM6SF2" in variants:
            genotype = variants["TM6SF2"].get("rs58542926", {}).get("genotype", "CC")
            if "T" in genotype:
                return {
                    "risk": "ELEVATED for liver disease",
                    "genotype": genotype,
                    "paradox": "PROTECTIVE for cardiovascular disease",
                    "note": "Fat accumulates in liver instead of arteries"
                }
        return {"risk": "LOW", "genotype": "CC"}

    def _genotype_adjusted_alt(self, alt: float, genetic_risk: str) -> dict:
        """
        Genotype-adjusted ALT interpretation

        Key insight: PNPLA3 carriers may have liver damage at lower ALT levels
        """
        if genetic_risk == "HIGH":
            threshold = 25  # Lower threshold for high-risk genotype
        elif genetic_risk == "MODERATE":
            threshold = 30
        else:
            threshold = 40  # Standard threshold

        if alt > threshold:
            return {
                "interpretation": "ELEVATED FOR GENOTYPE",
                "threshold_used": threshold,
                "note": f"Genotype-adjusted threshold is {threshold} U/L"
            }
        else:
            return {
                "interpretation": "Within genotype-adjusted range",
                "threshold_used": threshold
            }

    def _interpret_fib4(self, fib4: float) -> str:
        """Interpret FIB-4 fibrosis score"""
        if fib4 < 1.3:
            return "LOW risk of advanced fibrosis (NPV >90%)"
        elif fib4 < 2.67:
            return "INDETERMINATE - further evaluation recommended"
        else:
            return "HIGH risk of advanced fibrosis - hepatology referral recommended"

    def _fib4_genetic_context(self, fib4: float, genetic_risk: str) -> str:
        """Add genetic context to FIB-4 interpretation"""
        if fib4 < 1.3 and genetic_risk == "HIGH":
            return "Low FIB-4 but HIGH genetic risk - closer monitoring recommended"
        elif fib4 > 1.3 and genetic_risk == "HIGH":
            return "Both FIB-4 and genetics suggest elevated risk - consider elastography"
        return "FIB-4 consistent with genetic risk level"

    def _integrated_liver_assessment(self, result: dict) -> str:
        """Generate integrated liver health assessment"""
        genetic = result["genetic_risk"]["overall_genetic_susceptibility"]
        fibrosis = result["fibrosis_assessment"]["interpretation"]

        if genetic == "HIGH" and "HIGH" in fibrosis:
            return "HIGH RISK - Genetic susceptibility combined with fibrosis signal. Hepatology evaluation recommended."
        elif genetic == "HIGH":
            return "GENETICALLY SUSCEPTIBLE - Aggressive lifestyle intervention recommended even with normal enzymes"
        elif "HIGH" in fibrosis:
            return "FIBROSIS SIGNAL - Despite average genetics, enzymes suggest liver stress. Evaluate further."
        else:
            return "LOW RISK - Normal enzymes and average genetic risk. Maintain healthy lifestyle."

    def _liver_recommendations(self, result: dict, biomarkers: dict) -> list:
        """Generate liver health recommendations"""
        recommendations = []

        if result["genetic_risk"]["overall_genetic_susceptibility"] == "HIGH":
            recommendations.append({
                "priority": "HIGH",
                "intervention": "Lifestyle modification critical",
                "specifics": [
                    "Limit fructose and added sugars (<25g/day)",
                    "Avoid alcohol or strict moderation (<1 drink/day)",
                    "Weight loss if BMI >25 (5-10% body weight reduces liver fat 40%)",
                    "Mediterranean diet pattern",
                    "Coffee 2-3 cups/day (hepatoprotective)"
                ],
                "monitoring": "ALT every 6 months, consider FibroScan annually"
            })

        if biomarkers.get("ferritin", 100) > 300:
            recommendations.append({
                "priority": "MODERATE",
                "intervention": "Evaluate elevated ferritin",
                "rationale": "Elevated ferritin common in NAFLD, but rule out hemochromatosis",
                "action": "Check transferrin saturation, consider HFE genetic testing"
            })

        return recommendations
```

---

### 5. Thyroid Optimization with DIO1/DIO2 Genetics

```python
class ThyroidGeneticOptimizer:
    """
    Optimizes thyroid assessment using DIO1/DIO2 genetic variants

    Key insight: ~20% of patients don't feel well on T4 monotherapy
    due to impaired T4→T3 conversion. Genetics can explain this.
    """

    def __init__(self):
        self.deiodinase_genes = {
            "DIO1": {
                "rs2235544": {
                    "effect": "C allele: Increased DIO1 activity, higher FT3/FT4 ratio",
                    "clinical": "May do better on T4 monotherapy"
                },
                "rs11206244": {
                    "effect": "Affects FT3, rT3, and FT4 levels",
                    "clinical": "Contributes to thyroid hormone variability"
                }
            },
            "DIO2": {
                "rs225014": {
                    "name": "Thr92Ala",
                    "effect": "Ala allele: Reduced DIO2 activity, impaired T4→T3 conversion",
                    "clinical": "May explain persistent symptoms despite 'normal' labs",
                    "prevalence": "~16% homozygous AA, ~50% heterozygous",
                    "associations": "Insulin resistance, obesity, hypertension"
                }
            }
        }

    def analyze_thyroid_with_genetics(
        self,
        biomarkers: dict,
        genetic_variants: dict
    ) -> dict:
        """
        Genotype-informed thyroid analysis
        """
        tsh = biomarkers.get("tsh", 2.0)
        ft4 = biomarkers.get("free_t4", 1.2)
        ft3 = biomarkers.get("free_t3", 3.0)
        rt3 = biomarkers.get("reverse_t3", 15)
        tpo_ab = biomarkers.get("tpo_antibodies", 0)

        result = {
            "standard_assessment": {},
            "genetic_context": {},
            "conversion_analysis": {},
            "personalized_optimal_ranges": {},
            "recommendations": []
        }

        # Standard assessment
        result["standard_assessment"] = {
            "tsh": {"value": tsh, "status": self._interpret_tsh(tsh)},
            "ft4": {"value": ft4, "status": "Normal" if 0.8 <= ft4 <= 1.8 else "Abnormal"},
            "ft3": {"value": ft3, "status": "Normal" if 2.3 <= ft3 <= 4.2 else "Abnormal"},
            "tpo_antibodies": {"value": tpo_ab, "status": "Positive" if tpo_ab > 35 else "Negative"}
        }

        # Genetic context
        dio2_status = self._assess_dio2(genetic_variants)
        dio1_status = self._assess_dio1(genetic_variants)
        result["genetic_context"] = {
            "DIO2": dio2_status,
            "DIO1": dio1_status
        }

        # Conversion analysis
        if ft4 and ft3:
            ft3_ft4_ratio = ft3 / ft4
            result["conversion_analysis"] = {
                "ft3_ft4_ratio": round(ft3_ft4_ratio, 2),
                "interpretation": self._interpret_conversion(ft3_ft4_ratio, dio2_status),
                "rt3_analysis": self._analyze_rt3(rt3, ft3) if rt3 else None
            }

        # Personalized optimal ranges based on genetics
        result["personalized_optimal_ranges"] = self._personalized_thyroid_ranges(
            dio2_status, dio1_status
        )

        # Recommendations
        result["recommendations"] = self._thyroid_recommendations(
            result, biomarkers, genetic_variants
        )

        return result

    def _assess_dio2(self, variants: dict) -> dict:
        """Assess DIO2 Thr92Ala (rs225014)"""
        if "DIO2" in variants:
            genotype = variants["DIO2"].get("rs225014", {}).get("genotype", "CC")
            if genotype == "AA":
                return {
                    "status": "IMPAIRED CONVERSION",
                    "genotype": "AA (Ala/Ala homozygous)",
                    "enzyme_activity": "Reduced ~20-50%",
                    "clinical_significance": [
                        "May have lower tissue T3 despite normal serum levels",
                        "Higher rates of persistent hypothyroid symptoms on T4 therapy",
                        "May benefit from combination T4/T3 therapy",
                        "Associated with insulin resistance and metabolic syndrome"
                    ]
                }
            elif genotype in ["CA", "AC"]:
                return {
                    "status": "MILDLY IMPAIRED",
                    "genotype": "CA (heterozygous)",
                    "enzyme_activity": "Modestly reduced",
                    "clinical_significance": [
                        "May have subtle conversion impairment",
                        "Monitor symptoms on T4 therapy"
                    ]
                }
        return {"status": "NORMAL", "genotype": "CC (Thr/Thr)", "enzyme_activity": "Normal"}

    def _assess_dio1(self, variants: dict) -> dict:
        """Assess DIO1 rs2235544"""
        if "DIO1" in variants:
            genotype = variants["DIO1"].get("rs2235544", {}).get("genotype", "AA")
            if "C" in genotype:
                return {
                    "status": "ENHANCED",
                    "genotype": genotype,
                    "effect": "Higher DIO1 activity, better peripheral T4→T3 conversion"
                }
        return {"status": "NORMAL", "genotype": "AA"}

    def _interpret_conversion(self, ratio: float, dio2_status: dict) -> str:
        """Interpret T4→T3 conversion with genetic context"""
        if ratio < 2.0:
            base = "LOW FT3:FT4 ratio suggests impaired conversion"
        elif ratio > 3.5:
            base = "HIGH FT3:FT4 ratio suggests efficient conversion"
        else:
            base = "Normal conversion ratio"

        if dio2_status["status"] == "IMPAIRED CONVERSION" and ratio < 2.5:
            return f"{base}. GENETICALLY CONFIRMED - DIO2 Ala/Ala explains low conversion."
        elif dio2_status["status"] == "IMPAIRED CONVERSION":
            return f"{base}. Note: DIO2 variant present but ratio is acceptable - tissue levels may still be low."

        return base

    def _analyze_rt3(self, rt3: float, ft3: float) -> dict:
        """Analyze reverse T3 in context"""
        rt3_ft3_ratio = rt3 / ft3 if ft3 > 0 else 0

        return {
            "rt3": rt3,
            "rt3_ft3_ratio": round(rt3_ft3_ratio, 2),
            "interpretation": "Elevated rT3:FT3 suggests stress/illness-related conversion shift" if rt3_ft3_ratio > 10 else "Normal"
        }

    def _personalized_thyroid_ranges(self, dio2: dict, dio1: dict) -> dict:
        """
        Generate personalized optimal ranges based on genetics
        """
        ranges = {
            "tsh": {"optimal": (0.5, 2.5), "note": "Standard optimal range"},
            "ft4": {"optimal": (1.0, 1.5), "note": "Mid-range often optimal"},
            "ft3": {"optimal": (3.0, 4.0), "note": "Upper third of range often optimal"}
        }

        if dio2["status"] == "IMPAIRED CONVERSION":
            ranges["ft3"]["optimal"] = (3.2, 4.2)
            ranges["ft3"]["note"] = "DIO2 variant: aim for upper range FT3"
            ranges["recommendation"] = "Consider checking FT3 in addition to TSH for monitoring"

        return ranges

    def _thyroid_recommendations(self, result: dict, biomarkers: dict, variants: dict) -> list:
        """Generate thyroid recommendations"""
        recommendations = []

        dio2_status = result["genetic_context"]["DIO2"]

        if dio2_status["status"] == "IMPAIRED CONVERSION":
            recommendations.append({
                "priority": "HIGH",
                "category": "Treatment Optimization",
                "recommendations": [
                    "If on levothyroxine (T4) with persistent symptoms despite normal TSH:",
                    "  - Discuss combination T4/T3 therapy with endocrinologist",
                    "  - Studies show DIO2 Ala/Ala carriers have improved well-being on combination therapy",
                    "Monitor Free T3 levels, not just TSH",
                    "Selenium 200mcg daily may support deiodinase function"
                ],
                "evidence": "Panicker V et al. J Clin Endocrinol Metab 2009"
            })

            recommendations.append({
                "priority": "MODERATE",
                "category": "Metabolic Considerations",
                "recommendations": [
                    "DIO2 Thr92Ala associated with insulin resistance",
                    "Monitor fasting glucose and HOMA-IR",
                    "Exercise and weight management particularly important"
                ]
            })

        # Autoimmunity
        if biomarkers.get("tpo_antibodies", 0) > 35:
            recommendations.append({
                "priority": "MODERATE",
                "category": "Autoimmune Management",
                "recommendations": [
                    "Positive TPO antibodies indicate Hashimoto's thyroiditis",
                    "Selenium 200mcg daily may reduce antibody levels",
                    "Gluten-free diet may benefit some patients",
                    "Monitor TSH every 6-12 months for progression"
                ]
            })

        return recommendations
```

---

### 6. Iron Status and HFE Hemochromatosis

```python
class IronMetabolismAnalyzer:
    """
    Comprehensive iron status analysis with HFE genetic testing

    Key insight: Hereditary hemochromatosis affects 1 in 200-300 people
    of Northern European descent but is often undiagnosed until organ damage occurs.
    """

    def __init__(self):
        self.hfe_variants = {
            "C282Y": {
                "rsid": "rs1800562",
                "effect": "Major hemochromatosis mutation",
                "homozygous_risk": "85-90% of clinical hemochromatosis cases",
                "penetrance": "~10% develop clinical disease"
            },
            "H63D": {
                "rsid": "rs1799945",
                "effect": "Minor HFE variant",
                "compound_het_risk": "C282Y/H63D: moderate risk"
            }
        }

    def analyze_iron_status(
        self,
        biomarkers: dict,
        genetic_variants: dict,
        sex: str
    ) -> dict:
        """
        Comprehensive iron analysis with genetic context
        """
        ferritin = biomarkers.get("ferritin", 100)
        iron = biomarkers.get("iron", 100)
        tibc = biomarkers.get("tibc", 300)
        transferrin_sat = (iron / tibc * 100) if tibc > 0 else None

        result = {
            "iron_status": {},
            "hfe_genetics": {},
            "hemochromatosis_risk": {},
            "recommendations": []
        }

        # Calculate transferrin saturation if not provided
        if transferrin_sat is None:
            transferrin_sat = biomarkers.get("transferrin_saturation", 30)

        # Iron status assessment
        result["iron_status"] = self._assess_iron_status(
            ferritin, transferrin_sat, sex
        )

        # HFE genetic analysis
        result["hfe_genetics"] = self._assess_hfe_genetics(genetic_variants)

        # Hemochromatosis risk
        result["hemochromatosis_risk"] = self._assess_hemochromatosis_risk(
            result["iron_status"], result["hfe_genetics"], ferritin
        )

        # Recommendations
        result["recommendations"] = self._iron_recommendations(result, biomarkers, sex)

        return result

    def _assess_iron_status(self, ferritin: float, tsat: float, sex: str) -> dict:
        """Assess iron status"""
        # Sex-specific ferritin thresholds
        ferritin_high = 300 if sex == "male" else 200
        ferritin_concern = 500 if sex == "male" else 300
        ferritin_danger = 1000

        status = {
            "ferritin": {
                "value": ferritin,
                "interpretation": ""
            },
            "transferrin_saturation": {
                "value": tsat,
                "interpretation": ""
            }
        }

        # Ferritin interpretation
        if ferritin < 30:
            status["ferritin"]["interpretation"] = "LOW - Iron deficiency"
        elif ferritin <= ferritin_high:
            status["ferritin"]["interpretation"] = "NORMAL"
        elif ferritin <= ferritin_concern:
            status["ferritin"]["interpretation"] = "ELEVATED - Investigate cause"
        elif ferritin <= ferritin_danger:
            status["ferritin"]["interpretation"] = "HIGH - Possible iron overload"
        else:
            status["ferritin"]["interpretation"] = "VERY HIGH - Cirrhosis risk if >1000"

        # Transferrin saturation interpretation
        if tsat < 20:
            status["transferrin_saturation"]["interpretation"] = "LOW - Suggests iron deficiency"
        elif tsat <= 45:
            status["transferrin_saturation"]["interpretation"] = "NORMAL"
        elif tsat <= 55:
            status["transferrin_saturation"]["interpretation"] = "ELEVATED - Screen for hemochromatosis"
        else:
            status["transferrin_saturation"]["interpretation"] = "HIGH - Hemochromatosis likely if confirmed"

        return status

    def _assess_hfe_genetics(self, variants: dict) -> dict:
        """Assess HFE genotype for hemochromatosis"""
        hfe = variants.get("HFE", {})
        c282y = hfe.get("C282Y", {}).get("genotype", "CC")  # Wild-type is CC
        h63d = hfe.get("H63D", {}).get("genotype", "CC")     # Wild-type is CC

        # Interpret combination
        if c282y == "AA":  # Homozygous C282Y
            return {
                "genotype": "C282Y/C282Y (homozygous)",
                "risk": "HIGH",
                "interpretation": "Classic hemochromatosis genotype",
                "penetrance": "~28% of males, ~1% of females develop clinical disease",
                "action": "REQUIRES MONITORING even if ferritin normal now"
            }
        elif c282y == "CA" and h63d in ["CD", "DD"]:  # Compound heterozygote
            return {
                "genotype": "C282Y/H63D (compound heterozygote)",
                "risk": "MODERATE",
                "interpretation": "Elevated iron overload risk",
                "action": "Monitor ferritin and transferrin saturation"
            }
        elif c282y == "CA":  # C282Y carrier
            return {
                "genotype": "C282Y carrier (heterozygous)",
                "risk": "LOW",
                "interpretation": "Carrier status, low personal risk",
                "action": "Family screening may be indicated"
            }
        elif h63d in ["CD", "DD"]:
            return {
                "genotype": f"H63D {'homozygous' if h63d == 'DD' else 'heterozygous'}",
                "risk": "MINIMAL",
                "interpretation": "H63D alone rarely causes clinical iron overload"
            }

        return {"genotype": "Wild-type", "risk": "LOW", "interpretation": "No HFE mutations detected"}

    def _assess_hemochromatosis_risk(
        self,
        iron_status: dict,
        genetics: dict,
        ferritin: float
    ) -> dict:
        """Integrated hemochromatosis risk assessment"""
        genetic_risk = genetics["risk"]
        tsat_high = "ELEVATED" in iron_status["transferrin_saturation"]["interpretation"] or \
                    "HIGH" in iron_status["transferrin_saturation"]["interpretation"]
        ferritin_high = ferritin > 300

        if genetic_risk == "HIGH" and tsat_high:
            return {
                "risk_level": "CONFIRMED HIGH RISK",
                "interpretation": "HFE C282Y/C282Y with elevated iron markers",
                "cirrhosis_risk": "HIGH if ferritin >1000" if ferritin > 1000 else "Moderate if ferritin rises",
                "organs_at_risk": ["Liver (cirrhosis, HCC)", "Heart (cardiomyopathy)",
                                   "Pancreas (diabetes)", "Joints (arthritis)",
                                   "Pituitary (hypogonadism)"],
                "action": "Therapeutic phlebotomy indicated - hematology referral"
            }
        elif genetic_risk == "HIGH" and ferritin < 300:
            return {
                "risk_level": "HIGH GENETIC RISK - NOT YET EXPRESSED",
                "interpretation": "C282Y homozygote but ferritin still normal",
                "note": "Iron accumulation is gradual - will likely rise without intervention",
                "action": "Monitor ferritin/transferrin saturation every 6-12 months"
            }
        elif genetic_risk == "MODERATE" and (tsat_high or ferritin_high):
            return {
                "risk_level": "MODERATE - EXPRESSING",
                "interpretation": "Compound heterozygote with evidence of iron loading",
                "action": "Consider therapeutic phlebotomy if ferritin continues rising"
            }
        elif ferritin_high and genetic_risk == "LOW":
            return {
                "risk_level": "ELEVATED FERRITIN - NON-GENETIC CAUSE LIKELY",
                "interpretation": "Elevated ferritin without HFE mutations",
                "investigate": [
                    "Inflammation (check CRP)",
                    "Metabolic syndrome/NAFLD (most common cause)",
                    "Alcohol",
                    "Non-HFE hemochromatosis (rare)"
                ]
            }

        return {"risk_level": "LOW", "interpretation": "Normal iron status and genetics"}

    def _iron_recommendations(self, result: dict, biomarkers: dict, sex: str) -> list:
        """Generate iron-related recommendations"""
        recommendations = []

        if result["hemochromatosis_risk"]["risk_level"] in ["CONFIRMED HIGH RISK", "HIGH GENETIC RISK - NOT YET EXPRESSED"]:
            recommendations.append({
                "priority": "URGENT" if "CONFIRMED" in result["hemochromatosis_risk"]["risk_level"] else "HIGH",
                "category": "Hemochromatosis Management",
                "actions": [
                    "Hematology referral for phlebotomy consideration",
                    "Avoid iron supplements and vitamin C supplements",
                    "Limit red meat consumption",
                    "Avoid alcohol (accelerates liver damage)",
                    "Screen first-degree relatives (siblings 25% chance of same genotype)"
                ],
                "monitoring": "Ferritin every 3-6 months until stable"
            })

        # Iron deficiency
        if result["iron_status"]["ferritin"]["interpretation"] == "LOW - Iron deficiency":
            recommendations.append({
                "priority": "MODERATE",
                "category": "Iron Deficiency",
                "actions": [
                    "Investigate cause (GI blood loss, heavy menstruation, celiac)",
                    "Consider iron supplementation with vitamin C",
                    "Recheck ferritin in 3 months"
                ]
            })

        return recommendations
```

---

### 7. Personalized Nutrition: FADS Gene-Diet Interactions

```python
class PersonalizedNutritionAnalyzer:
    """
    Personalized omega-3/omega-6 recommendations based on FADS1/FADS2 genetics

    Key insight: FADS genetic variants affect fatty acid metabolism and
    inflammation, meaning dietary requirements vary by genotype.
    """

    def __init__(self):
        self.fads_genes = {
            "FADS1": {
                "rs174546": {
                    "risk_allele": "C",
                    "effect": "Lower FADS1 activity, less efficient PUFA conversion",
                    "clinical": "May need more preformed EPA/DHA"
                }
            },
            "FADS2": {
                "rs1535": {
                    "risk_allele": "G",
                    "effect": "GG genotype: increased inflammatory conditions, benefit from omega-3s",
                    "clinical": "Greater therapeutic response to fish oil"
                }
            }
        }

    def analyze_fads_and_recommend(
        self,
        biomarkers: dict,
        genetic_variants: dict
    ) -> dict:
        """
        Generate personalized fatty acid recommendations
        """
        omega3_index = biomarkers.get("omega3_index", None)
        arachidonic_acid = biomarkers.get("arachidonic_acid", None)
        epa = biomarkers.get("epa", None)
        dha = biomarkers.get("dha", None)

        result = {
            "fads_genetics": {},
            "omega_status": {},
            "personalized_recommendation": {}
        }

        # Genetic assessment
        fads1_status = self._assess_fads1(genetic_variants)
        fads2_status = self._assess_fads2(genetic_variants)
        result["fads_genetics"] = {
            "FADS1": fads1_status,
            "FADS2": fads2_status,
            "overall": self._overall_fads_phenotype(fads1_status, fads2_status)
        }

        # Omega status if available
        if omega3_index:
            result["omega_status"] = {
                "omega3_index": {
                    "value": omega3_index,
                    "interpretation": self._interpret_omega3_index(omega3_index)
                }
            }

        # Personalized recommendations
        result["personalized_recommendation"] = self._generate_fads_recommendations(
            result["fads_genetics"]["overall"], omega3_index
        )

        return result

    def _assess_fads1(self, variants: dict) -> dict:
        """Assess FADS1 genotype"""
        if "FADS1" in variants:
            genotype = variants["FADS1"].get("rs174546", {}).get("genotype", "TT")
            if genotype == "CC":
                return {
                    "genotype": "CC",
                    "activity": "LOW",
                    "effect": "Less efficient conversion of ALA to EPA/DHA",
                    "recommendation": "Emphasize direct EPA/DHA sources"
                }
            elif genotype == "CT":
                return {"genotype": "CT", "activity": "MODERATE"}
        return {"genotype": "TT", "activity": "NORMAL"}

    def _assess_fads2(self, variants: dict) -> dict:
        """Assess FADS2 genotype"""
        if "FADS2" in variants:
            genotype = variants["FADS2"].get("rs1535", {}).get("genotype", "AA")
            if genotype == "GG":
                return {
                    "genotype": "GG",
                    "effect": "Higher inflammatory risk, enhanced benefit from omega-3s",
                    "therapeutic_response": "ENHANCED response to fish oil supplementation"
                }
        return {"genotype": "AA", "effect": "Standard response to omega-3s"}

    def _overall_fads_phenotype(self, fads1: dict, fads2: dict) -> dict:
        """Determine overall FADS phenotype"""
        if fads1.get("activity") == "LOW" and fads2.get("therapeutic_response") == "ENHANCED":
            return {
                "phenotype": "HIGH RESPONDER - CONVERSION LIMITED",
                "interpretation": [
                    "Low endogenous EPA/DHA synthesis",
                    "But high therapeutic response to direct supplementation",
                    "Direct fish/fish oil particularly important"
                ]
            }
        elif fads1.get("activity") == "LOW":
            return {
                "phenotype": "CONVERSION LIMITED",
                "interpretation": "Dietary omega-3 plant sources (ALA) less effective; need preformed EPA/DHA"
            }
        elif fads2.get("therapeutic_response") == "ENHANCED":
            return {
                "phenotype": "HIGH RESPONDER",
                "interpretation": "Greater-than-average benefit from omega-3 supplementation"
            }
        return {"phenotype": "NORMAL", "interpretation": "Standard fatty acid metabolism"}

    def _interpret_omega3_index(self, index: float) -> dict:
        """Interpret Omega-3 Index"""
        if index < 4:
            return {"risk": "HIGH CV RISK", "recommendation": "Aggressive supplementation needed"}
        elif index < 8:
            return {"risk": "MODERATE", "recommendation": "Increase omega-3 intake"}
        else:
            return {"risk": "LOW (CARDIOPROTECTIVE)", "recommendation": "Maintain current intake"}

    def _generate_fads_recommendations(self, fads_overall: dict, omega3_index: float) -> dict:
        """Generate personalized fatty acid recommendations"""
        phenotype = fads_overall.get("phenotype", "NORMAL")

        base_recommendations = {
            "fatty_fish_servings": 2,
            "fish_oil_dose": "1-2g EPA+DHA daily",
            "plant_omega3": "Moderate (ALA from flax, chia, walnuts)"
        }

        if "CONVERSION LIMITED" in phenotype:
            base_recommendations["fatty_fish_servings"] = 3
            base_recommendations["fish_oil_dose"] = "2-4g EPA+DHA daily"
            base_recommendations["plant_omega3"] = "Less effective due to genetics - prioritize fish/fish oil"
            base_recommendations["note"] = "Your genetics suggest limited ALA→EPA/DHA conversion. Plant omega-3s alone insufficient."

        if "HIGH RESPONDER" in phenotype:
            base_recommendations["fish_oil_dose"] = "2-4g EPA+DHA daily (enhanced benefit expected)"
            base_recommendations["note"] = "Your genetics suggest greater-than-average benefit from omega-3 supplementation"

        if omega3_index and omega3_index < 4:
            base_recommendations["urgency"] = "HIGH - Current omega-3 index indicates elevated CV risk"
            base_recommendations["target"] = "Omega-3 Index >8%"

        return base_recommendations
```

---

## Part II: Complete Agent Implementation

### Main Orchestration Class

```python
class FunctionHealthCorrelationAgent:
    """
    Main orchestration class for the Function Health Biomarker-Genomic Correlation Agent

    This agent integrates all analysis modules to provide comprehensive,
    personalized health insights.
    """

    def __init__(self, vast_db_connection, llm_client):
        self.db = vast_db_connection
        self.llm = llm_client

        # Initialize all analyzers
        self.bio_age = BiologicalAgeCalculator()
        self.diabetes_trajectory = PreDiabeticTrajectoryAnalyzer()
        self.cv_analyzer = AdvancedCardiovascularAnalyzer()
        self.liver_analyzer = LiverHealthAnalyzer()
        self.thyroid_analyzer = ThyroidGeneticOptimizer()
        self.iron_analyzer = IronMetabolismAnalyzer()
        self.nutrition_analyzer = PersonalizedNutritionAnalyzer()
        self.mthfr_analyzer = MTHFRHomocysteineAnalyzer()
        self.apoe_analyzer = ApoEAnalyzer()

    def generate_comprehensive_report(
        self,
        biomarkers: dict,
        genomic_variants: dict,
        patient_info: dict
    ) -> dict:
        """
        Generate comprehensive personalized health report
        """
        age = patient_info.get("age", 45)
        sex = patient_info.get("sex", "male")

        report = {
            "patient_info": patient_info,
            "report_date": datetime.now().isoformat(),
            "biological_age": {},
            "disease_trajectories": {},
            "cardiovascular": {},
            "metabolic": {},
            "liver": {},
            "thyroid": {},
            "iron": {},
            "nutrition": {},
            "pharmacogenomics": {},
            "executive_summary": {},
            "action_items": []
        }

        # 1. Biological Age
        report["biological_age"] = self.bio_age.calculate_phenoage(biomarkers, age)

        # 2. Disease Trajectories
        report["disease_trajectories"]["diabetes"] = self.diabetes_trajectory.calculate_diabetes_trajectory(
            biomarkers, genomic_variants, age
        )

        # 3. Cardiovascular
        report["cardiovascular"] = self.cv_analyzer.comprehensive_cv_assessment(
            biomarkers, genomic_variants, age, sex
        )

        # 4. Liver
        report["liver"] = self.liver_analyzer.analyze_liver_health(
            biomarkers, genomic_variants
        )

        # 5. Thyroid
        report["thyroid"] = self.thyroid_analyzer.analyze_thyroid_with_genetics(
            biomarkers, genomic_variants
        )

        # 6. Iron
        report["iron"] = self.iron_analyzer.analyze_iron_status(
            biomarkers, genomic_variants, sex
        )

        # 7. Nutrition
        report["nutrition"] = self.nutrition_analyzer.analyze_fads_and_recommend(
            biomarkers, genomic_variants
        )

        # 8. MTHFR/Homocysteine
        if "MTHFR" in genomic_variants:
            report["methylation"] = self.mthfr_analyzer.analyze_mthfr_homocysteine(
                genomic_variants["MTHFR"], biomarkers
            )

        # 9. ApoE
        if "APOE" in genomic_variants:
            report["apoe_analysis"] = self.apoe_analyzer.analyze_apoe_comprehensive(
                genomic_variants["APOE"], biomarkers
            )

        # Generate executive summary using LLM
        report["executive_summary"] = self._generate_executive_summary(report)

        # Compile prioritized action items
        report["action_items"] = self._compile_action_items(report)

        return report

    def _generate_executive_summary(self, report: dict) -> dict:
        """
        Use LLM to generate executive summary
        """
        prompt = f"""
        Based on the following comprehensive health analysis, generate an executive summary
        highlighting the TOP 5 most important findings and their implications.

        Focus on:
        1. Most actionable findings
        2. Genetic factors that explain biomarker patterns
        3. Pre-symptomatic disease risks identified
        4. Drug-gene interactions to be aware of
        5. Personalized interventions with highest expected impact

        Analysis Data:
        {json.dumps(report, indent=2, default=str)}

        Format as:
        - Key Finding 1: [Finding] -> [Action]
        - Key Finding 2: [Finding] -> [Action]
        ...
        """

        summary_text = self.llm.generate(prompt, max_tokens=1500)

        return {
            "summary": summary_text,
            "biological_age_delta": report["biological_age"].get("age_acceleration", 0),
            "highest_priority_items": self._extract_priorities(report)
        }

    def _compile_action_items(self, report: dict) -> list:
        """
        Compile all action items across analyses, sorted by priority
        """
        action_items = []

        # Collect from all sections
        sections_with_recommendations = [
            ("cardiovascular", "personalized_interventions"),
            ("liver", "recommendations"),
            ("thyroid", "recommendations"),
            ("iron", "recommendations"),
            ("disease_trajectories", "interventions"),
        ]

        for section, key in sections_with_recommendations:
            if section in report and key in report[section]:
                items = report[section][key]
                if isinstance(items, list):
                    for item in items:
                        item["source"] = section
                        action_items.append(item)
                elif isinstance(items, dict) and "interventions" in items:
                    for item in items["interventions"]:
                        item["source"] = section
                        action_items.append(item)

        # Sort by priority
        priority_order = {"URGENT": 0, "HIGH": 1, "MODERATE": 2, "LOW": 3}
        action_items.sort(key=lambda x: priority_order.get(x.get("priority", "LOW"), 4))

        return action_items

    def _extract_priorities(self, report: dict) -> list:
        """Extract highest priority items"""
        priorities = []

        # Biological age acceleration
        if report["biological_age"].get("age_acceleration", 0) > 5:
            priorities.append({
                "category": "Biological Aging",
                "finding": f"Aging {report['biological_age']['age_acceleration']:.1f} years faster than chronological age",
                "urgency": "HIGH"
            })

        # Pre-diabetes
        if report.get("disease_trajectories", {}).get("diabetes", {}).get("current_state") == "PREDIABETES":
            priorities.append({
                "category": "Metabolic Health",
                "finding": "Pre-diabetic trajectory detected",
                "urgency": "HIGH"
            })

        # High Lp(a)
        if report.get("cardiovascular", {}).get("advanced_lipid_risk", {}).get("lp_a", {}).get("interpretation", {}).get("risk") == "VERY HIGH":
            priorities.append({
                "category": "Cardiovascular",
                "finding": "Very high Lp(a) - genetic cardiovascular risk factor",
                "urgency": "HIGH"
            })

        return priorities
```

---

## Part III: Implementation Roadmap

### Phase 1: Core Infrastructure (Weeks 1-2)
- [ ] VAST DataBase schema for biomarker-gene correlations
- [ ] Function Health data parser (PDF/CSV)
- [ ] Genomic data integration from Pipeline 1

### Phase 2: Analysis Engines (Weeks 3-4)
- [ ] Biological Age Calculator
- [ ] Pre-Diabetic Trajectory Analyzer
- [ ] Cardiovascular Risk Engine
- [ ] Liver Health Analyzer

### Phase 3: Genetic Integrations (Weeks 5-6)
- [ ] MTHFR-Homocysteine correlator
- [ ] ApoE-Lipid-Cognitive correlator
- [ ] DIO2-Thyroid correlator
- [ ] FADS-Nutrition optimizer
- [ ] HFE-Iron correlator

### Phase 4: Agent Orchestration (Weeks 7-8)
- [ ] Main orchestration class
- [ ] LLM report generation
- [ ] Action item prioritization
- [ ] Patient-friendly summary generation

### Phase 5: Testing & Validation (Weeks 9-10)
- [ ] Test with real Function Health data
- [ ] Clinical review of outputs
- [ ] Iterate on recommendations
- [ ] Documentation and GitHub release

---

## References

### Biological Aging
1. [GDF-15 as proxy for epigenetic aging](https://pmc.ncbi.nlm.nih.gov/articles/PMC11625061/)
2. [PhenoAge - An epigenetic biomarker of aging](https://pmc.ncbi.nlm.nih.gov/articles/PMC5940111/)
3. [GrimAge epigenetic clock](https://www.aginganddisease.org/EN/10.14336/AD.2024.1495)

### Cardiovascular
4. [IL-6 and Cardiovascular Risk](https://pmc.ncbi.nlm.nih.gov/articles/PMC11599326/)
5. [IL6 genetic perturbation and cardiometabolic risk](https://www.nature.com/articles/s44161-025-00700-7)
6. [PCSK9 variants and LDL response](https://pmc.ncbi.nlm.nih.gov/articles/PMC4995153/)
7. [Lp(a) and cardiovascular risk](https://www.nejm.org/doi/full/10.1056/NEJMoa054013)

### Metabolic
8. [TCF7L2 and insulin resistance](https://pubmed.ncbi.nlm.nih.gov/19509102/)
9. [Prediabetes reversibility](https://pmc.ncbi.nlm.nih.gov/articles/PMC4116271/)
10. [HbA1c genetic determinants](https://pmc.ncbi.nlm.nih.gov/articles/PMC3207128/)

### Liver
11. [PNPLA3, TM6SF2 in NAFLD](https://pmc.ncbi.nlm.nih.gov/articles/PMC8382644/)
12. [Genetic prediction of NAFLD trajectories](https://www.gastrojournal.org/article/S0016-5085(20)30229-8/fulltext)

### Thyroid
13. [DIO2 Thr92Ala polymorphism](https://www.frontiersin.org/articles/10.3389/fendo.2019.00912/full)
14. [DIO1 genetic variation](https://pmc.ncbi.nlm.nih.gov/articles/PMC2515080/)

### Iron
15. [HFE Hemochromatosis](https://www.ncbi.nlm.nih.gov/books/NBK1440/)
16. [Ferritin and cirrhosis risk](https://pmc.ncbi.nlm.nih.gov/articles/PMC2275006/)

### Inflammation
17. [Neutrophil-Lymphocyte Ratio and CV disease](https://pmc.ncbi.nlm.nih.gov/articles/PMC9687310/)
18. [hs-CRP genetic variants](https://pubmed.ncbi.nlm.nih.gov/29627531/)

### Nutrition
19. [FADS genetic variation and omega-3](https://www.frontiersin.org/journals/nutrition/articles/10.3389/fnut.2025.1538505/full)
20. [FADS polymorphisms and cardiovascular health](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0222061)

### Multi-Omics
21. [Multi-omics profiling for health](https://pmc.ncbi.nlm.nih.gov/articles/PMC10220275/)
22. [AI-driven multi-omics integration](https://www.sciencedirect.com/science/article/pii/S2001037024004513)

---

*HCLS AI Factory - Function Health Biomarker-Genomic Correlation Agent v2.0*
*Precision Medicine Through Integrated Analysis*
*Apache 2.0 Licensed | Open Source*
