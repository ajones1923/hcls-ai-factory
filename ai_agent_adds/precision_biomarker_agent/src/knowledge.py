"""Biomarker Intelligence Agent -- Knowledge Graph.

Extends the Clinker pattern from rag-chat-pipeline/src/knowledge.py and
mirrors cart_intelligence_agent/src/knowledge.py, adapted for precision
biomarker analysis. Contains:

1. BIOMARKER_DOMAINS: 7 disease domains with biomarkers, genetic modifiers, and
   intervention targets
2. PHENOAGE_KNOWLEDGE: PhenoAge clock biomarker descriptions, coefficients,
   and clinical interpretation
3. PGX_KNOWLEDGE: 9 pharmacogenes with key drug interactions and CPIC guidance
4. CROSS_MODAL_LINKS: 8 cross-modal links mapping biomarker findings to triggers
   for other HCLS AI Factory agents (imaging, oncology, genomics)

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict

# Knowledge base version tracking
# Update these dates when clinical guidelines or thresholds change
KNOWLEDGE_VERSION = {
    "version": "1.0.0",
    "last_updated": "2026-03-01",
    "cpic_version": "March 2025",  # CPIC guideline versions used
    "ada_standards": "2024",       # ADA Standards of Medical Care
    "esc_guidelines": "2021",      # ESC/EAS Dyslipidemia Guidelines
    "aasld_guidelines": "2023",    # AASLD NAFLD/NASH Guidance
    "levine_phenoage": "2018",     # PhenoAge algorithm (Levine et al.)
    "lu_grimage": "2019",          # GrimAge algorithm (Lu et al.)
    "sources": [
        "CPIC (cpicpgx.org) - Pharmacogenomic guidelines",
        "ADA Standards of Medical Care in Diabetes 2024",
        "ESC/EAS Guidelines for Dyslipidemia Management 2021",
        "AASLD Practice Guidance on NAFLD 2023",
        "Levine et al. 2018 PMID:29676998 - PhenoAge",
        "Lu et al. 2019 PMID:30669119 - GrimAge",
        "Nordestgaard et al. 2010 PMID:20031622 - Lp(a)",
    ],
}


# =============================================================================
# 0. SHARED CLINICAL THRESHOLDS -- for cross-module consistency
# =============================================================================

# Shared clinical thresholds for cross-module consistency
# Florez et al. 2012 PMID:22399527; Grant et al. 2006 PMID:16415884
GENOTYPE_THRESHOLDS = {
    "TCF7L2_hba1c": {0: 6.0, 1: 5.8, 2: 5.5},  # risk_alleles: threshold (%)
    "TCF7L2_fasting_glucose": {0: 100, 1: 95, 2: 90},  # mg/dL
    "PNPLA3_alt_upper": {"CC": 56, "CG": 45, "GG": 35},  # U/L; Romeo et al. 2008 PMID:18820127; Sookoian & Pirola 2011 PMID:21520172
    "DIO2_tsh_upper": {"GG": 4.0, "GA": 3.5, "AA": 3.0},  # mIU/L; Panicker et al. 2009 PMID:19820026; Castagna et al. 2017 PMID:28100792
    "APOL1_eGFR_threshold": {"low_risk": 60, "high_risk": 75},  # APOL1 two-risk-allele carriers need tighter eGFR monitoring
    "MTHFR_homocysteine_upper": {"CC": 15, "CT": 12, "TT": 10},  # umol/L; Frosst et al. 1995 PMID:7647779
}

# Age- and sex-stratified reference ranges for key biomarkers
# Sources: NHANES III, Harrison's Principles of Internal Medicine 21e,
# Tietz Clinical Guide to Laboratory Tests 6e
AGE_SEX_REFERENCE_RANGES = {
    "creatinine": {
        # mg/dL; Inker et al. 2021 PMID:34554658 (CKD-EPI 2021)
        "M": {"18-49": (0.7, 1.2), "50-69": (0.8, 1.3), "70+": (0.9, 1.5)},
        "F": {"18-49": (0.5, 1.0), "50-69": (0.6, 1.1), "70+": (0.7, 1.3)},
    },
    "alt": {
        # U/L; Prati et al. 2002 PMID:12407587; updated Kwo et al. 2017 PMID:27677441
        "M": {"18-49": (7, 45), "50-69": (7, 50), "70+": (7, 55)},
        "F": {"18-49": (7, 30), "50-69": (7, 35), "70+": (7, 40)},
    },
    "alkaline_phosphatase": {
        # U/L; varies with bone turnover (higher in elderly)
        "M": {"18-49": (44, 147), "50-69": (44, 147), "70+": (44, 165)},
        "F": {"18-49": (44, 147), "50-69": (44, 165), "70+": (44, 187)},
    },
    "ferritin": {
        # ng/mL; Adams et al. 2005 PMID:15818862
        "M": {"18-49": (24, 336), "50-69": (24, 336), "70+": (24, 400)},
        "F": {"18-49": (11, 150), "50-69": (12, 263), "70+": (12, 300)},
    },
    "tsh": {
        # mIU/L; Boucai et al. 2011 PMID:21270357
        "M": {"18-49": (0.4, 4.0), "50-69": (0.4, 4.5), "70+": (0.4, 6.0)},
        "F": {"18-49": (0.4, 4.0), "50-69": (0.4, 5.0), "70+": (0.4, 7.0)},
    },
    "hemoglobin": {
        # g/dL; WHO criteria
        "M": {"18-49": (13.5, 17.5), "50-69": (13.0, 17.0), "70+": (12.0, 16.5)},
        "F": {"18-49": (12.0, 16.0), "50-69": (11.5, 15.5), "70+": (11.0, 15.0)},
    },
    "bun": {
        # mg/dL; Inker et al. 2021 PMID:34554658
        "M": {"18-49": (7, 20), "50-69": (8, 23), "70+": (10, 28)},
        "F": {"18-49": (6, 18), "50-69": (7, 21), "70+": (8, 25)},
    },
    "cystatin_c": {
        # mg/L; Stevens et al. 2008 PMID:18701668
        "M": {"18-49": (0.56, 0.98), "50-69": (0.63, 1.08), "70+": (0.70, 1.21)},
        "F": {"18-49": (0.52, 0.90), "50-69": (0.58, 1.02), "70+": (0.65, 1.15)},
    },
    "homocysteine": {
        # umol/L; Refsum et al. 2004 PMID:15205206
        "M": {"18-49": (5, 15), "50-69": (6, 17), "70+": (7, 20)},
        "F": {"18-49": (4, 13), "50-69": (5, 15), "70+": (6, 18)},
    },
    "vitamin_d_25oh": {
        # ng/mL; Endocrine Society 2024 PMID:38828931
        "M": {"18-49": (30, 100), "50-69": (30, 100), "70+": (30, 100)},
        "F": {"18-49": (30, 100), "50-69": (30, 100), "70+": (30, 100)},
    },
}

# Plausible clinical ranges for input validation
# Values outside these ranges are likely data entry errors
# Sources: Tietz Clinical Guide, clinical experience
BIOMARKER_PLAUSIBLE_RANGES = {
    "albumin": (1.0, 6.0),        # g/dL
    "creatinine": (0.1, 25.0),    # mg/dL
    "glucose": (20, 600),          # mg/dL
    "fasting_glucose": (20, 600),  # mg/dL
    "hs_crp": (0.01, 300.0),      # mg/L
    "wbc": (0.5, 100.0),          # 10^3/μL
    "mcv": (50, 130),              # fL
    "rdw": (8.0, 30.0),           # %
    "lymphocyte_pct": (1.0, 80.0), # %
    "alkaline_phosphatase": (5, 2000), # U/L
    "hba1c": (2.0, 20.0),         # %
    "ldl": (10, 500),              # mg/dL
    "hdl": (5, 150),               # mg/dL
    "triglycerides": (20, 5000),   # mg/dL
    "alt": (1, 5000),              # U/L
    "ast": (1, 5000),              # U/L
    "tsh": (0.01, 100.0),         # mIU/L
    "ferritin": (1, 10000),        # ng/mL
    "hemoglobin": (3.0, 25.0),    # g/dL
    "platelets": (5, 2000),        # K/μL
    "lpa": (0.1, 500),             # nmol/L
    "homocysteine": (1.0, 100.0), # μmol/L
    "vitamin_d_25oh": (1.0, 200.0), # ng/mL
    "vitamin_b12": (50, 5000),     # pg/mL
    "folate_serum": (1.0, 50.0),  # ng/mL
    "egfr": (5, 150),              # mL/min/1.73m2
    "cystatin_c": (0.1, 10.0),    # mg/L
    "bun": (2, 100),               # mg/dL
    "urine_acr": (0, 5000),       # mg/g
    "potassium": (2.0, 8.0),      # mEq/L
    "calcium": (5.0, 15.0),       # mg/dL
    "pth": (5, 500),               # pg/mL
    "ctx": (0.01, 3.0),           # ng/mL (C-telopeptide)
    "p1np": (5, 300),              # mcg/L (Procollagen I N-propeptide)
    "omega3_index": (1.0, 20.0),  # %
    "vitamin_d_25oh_plausible": (1.0, 200.0),  # ng/mL (already exists as vitamin_d_25oh)
    "magnesium": (0.5, 5.0),      # mg/dL
    "zinc": (20, 200),             # mcg/dL
    "selenium": (20, 400),         # mcg/L
    "adiponectin": (1.0, 50.0),   # mcg/mL
    "igf1": (10, 1000),            # ng/mL
}


# =============================================================================
# 0b. ANCESTRY ADJUSTMENTS -- population-specific biomarker reference shifts
# =============================================================================

# Population-specific biomarker adjustments
# Sources: NHANES III, UK Biobank, Multi-Ethnic Study of Atherosclerosis (MESA)
ANCESTRY_ADJUSTMENTS = {
    "african": {
        # African ancestry: higher Lp(a), lower triglycerides, higher creatinine
        # Virani et al. 2020 PMID:31992061; Inker et al. 2021 PMID:34554658
        "lpa": {"threshold_multiplier": 1.0, "note": "Lp(a) levels 2-3x higher; standard thresholds apply but prevalence is higher"},
        "creatinine": {"threshold_multiplier": 1.15, "note": "Higher muscle mass raises baseline creatinine ~15%; eGFR equations now race-neutral (CKD-EPI 2021)"},
        "triglycerides": {"threshold_multiplier": 0.85, "note": "Lower baseline triglycerides; adjusted thresholds recommended"},
        "vitamin_d_25oh": {"threshold_multiplier": 0.80, "note": "Lower 25(OH)D levels but equivalent bone density; adjusted sufficiency threshold"},
        "hba1c": {"threshold_multiplier": 1.0, "note": "HbA1c may overestimate glycemia due to RBC lifespan differences; consider fructosamine"},
    },
    "south_asian": {
        # South Asian: higher cardiovascular risk at lower BMI/lipid levels
        # Gujral et al. 2013 PMID:23404868; Sattar & Gill 2015 PMID:25533203
        "ldl": {"threshold_multiplier": 0.90, "note": "Lower LDL thresholds recommended due to higher CVD risk at equivalent levels"},
        "hba1c": {"threshold_multiplier": 0.95, "note": "Diabetes screening recommended at lower HbA1c (5.5% vs 5.7%)"},
        "triglycerides": {"threshold_multiplier": 0.90, "note": "Lower triglyceride thresholds due to atherogenic dyslipidemia pattern"},
        "lpa": {"threshold_multiplier": 1.0, "note": "Elevated Lp(a) prevalence; standard thresholds apply"},
    },
    "east_asian": {
        # East Asian: different alcohol metabolism, statin sensitivity
        # Dang et al. 2014 PMID:24458560; Lee et al. 2019 PMID:31276099
        "alt": {"threshold_multiplier": 0.85, "note": "Lower ALT upper limits recommended (30 U/L male, 19 U/L female)"},
        "creatinine": {"threshold_multiplier": 0.95, "note": "Lower average muscle mass; slightly lower baseline creatinine"},
        "ldl": {"threshold_multiplier": 0.95, "note": "Statins more potent at equivalent doses; lower LDL targets achievable"},
    },
    "hispanic": {
        # Hispanic/Latino: higher NAFLD prevalence, diabetes risk
        # Lazo et al. 2013 PMID:23532013; Aguayo-Mazzucato et al. 2019 PMID:30737276
        "alt": {"threshold_multiplier": 1.10, "note": "Higher NAFLD prevalence; ALT may be elevated at baseline"},
        "hba1c": {"threshold_multiplier": 0.95, "note": "Earlier diabetes screening recommended; lower HbA1c threshold"},
        "triglycerides": {"threshold_multiplier": 1.05, "note": "Higher baseline triglycerides common"},
    },
}


# =============================================================================
# 1. BIOMARKER_DOMAINS -- 6 disease domain knowledge graphs
# =============================================================================

BIOMARKER_DOMAINS: Dict[str, Dict[str, Any]] = {
    "diabetes": {
        "name": "Diabetes / Metabolic Syndrome",
        "key_biomarkers": {
            "HbA1c": {
                "unit": "%",
                "normal_range": "4.0-5.6",
                "pre_disease": "5.7-6.4",
                "disease": ">=6.5",
                "clinical_note": "Glycated hemoglobin reflects 2-3 month average glucose.",
            },
            "fasting_glucose": {
                "unit": "mg/dL",
                "normal_range": "70-99",
                "pre_disease": "100-125",
                "disease": ">=126",
                "clinical_note": "Fasting plasma glucose; confirm on repeat testing.",
            },
            "fasting_insulin": {
                "unit": "uIU/mL",
                "normal_range": "2.6-24.9",
                "elevated": ">25",
                "clinical_note": "Elevated insulin with normal glucose suggests insulin resistance.",
            },
            "HOMA_IR": {
                "unit": "index",
                "normal_range": "<2.0",
                "pre_disease": "2.0-2.9",
                "disease": ">=3.0",
                "clinical_note": "Homeostatic model assessment; (glucose x insulin) / 405.",
            },
            "triglycerides": {
                "unit": "mg/dL",
                "normal_range": "<150",
                "borderline": "150-199",
                "elevated": ">=200",
                "clinical_note": "Part of metabolic syndrome criteria.",
            },
        },
        "genetic_modifiers": [
            {"gene": "TCF7L2", "rs_id": "rs7903146", "risk_allele": "T",
             "effect": "1.4x increased T2DM risk per allele; strongest common T2DM variant"},
            {"gene": "SLC30A8", "rs_id": "rs13266634", "risk_allele": "C",
             "effect": "Zinc transporter affecting insulin secretion"},
            {"gene": "PPARG", "rs_id": "rs1801282", "risk_allele": "C",
             "effect": "Pro12Ala; affects insulin sensitivity and adipogenesis"},
        ],
        "intervention_targets": [
            "Lifestyle modification (diet, exercise) -- 58% risk reduction (DPP trial)",
            "Metformin for pre-diabetes with HOMA-IR >3.0",
            "Berberine for insulin sensitization (TCF7L2 risk carriers)",
            "Chromium picolinate 200-1000 mcg/day for insulin sensitivity",
        ],
        "clinical_context": (
            "Pre-diabetes is detectable 5-10 years before clinical T2DM diagnosis. "
            "The biomarker trajectory of rising fasting insulin followed by rising "
            "glucose followed by rising HbA1c provides a window for intervention. "
            "Genetic risk from TCF7L2 can accelerate this trajectory."
        ),
    },

    "cardiovascular": {
        "name": "Cardiovascular Risk",
        "key_biomarkers": {
            "Lp(a)": {
                "unit": "nmol/L",
                "normal_range": "<75",
                "elevated": "75-125",
                "high_risk": ">125",
                "clinical_note": "Genetically determined; 90% of variation due to LPA gene.",
            },
            "ApoB": {
                "unit": "mg/dL",
                "normal_range": "<90",
                "borderline": "90-130",
                "elevated": ">130",
                "clinical_note": "Better predictor of cardiovascular events than LDL-C.",
            },
            "hs_CRP": {
                "unit": "mg/L",
                "normal_range": "<1.0",
                "borderline": "1.0-3.0",
                "elevated": ">3.0",
                "clinical_note": "Marker of systemic inflammation and vascular risk.",
            },
            "homocysteine": {
                "unit": "umol/L",
                "normal_range": "5-15",
                "elevated": ">15",
                "clinical_note": "Elevated by MTHFR variants; independent CVD risk factor.",
            },
            "LDL_C": {
                "unit": "mg/dL",
                "normal_range": "<100",
                "borderline": "100-159",
                "elevated": ">=160",
                "clinical_note": "Primary target for statin therapy; PCSK9 variants affect levels.",
            },
        },
        "genetic_modifiers": [
            {"gene": "APOE", "rs_id": "rs429358/rs7412", "risk_allele": "e4",
             "effect": "APOE4 increases LDL-C by 10-15 mg/dL; 2-3x Alzheimer risk"},
            {"gene": "LPA", "rs_id": "rs10455872", "risk_allele": "G",
             "effect": "KIV-2 repeats determine Lp(a) levels; not modifiable by statins"},
            {"gene": "PCSK9", "rs_id": "rs11591147", "risk_allele": "T",
             "effect": "Loss-of-function variants reduce LDL-C and CVD risk by ~50%"},
            {"gene": "MTHFR", "rs_id": "rs1801133", "risk_allele": "T",
             "effect": "C677T reduces enzyme activity 30-70%; raises homocysteine"},
        ],
        "intervention_targets": [
            "PCSK9 inhibitors for high Lp(a) with established CVD",
            "Niacin 1-3g/day may reduce Lp(a) by 20-30%",
            "Methylfolate + B12 for MTHFR-related hyperhomocysteinemia",
            "ApoB-guided statin therapy (target <80 mg/dL for high risk)",
            "Omega-3 fatty acids for triglyceride reduction",
        ],
        "clinical_context": (
            "Lp(a) is the single most heritable cardiovascular risk factor. "
            "ApoB better captures atherogenic particle burden than LDL-C alone. "
            "APOE genotype affects both lipid metabolism and Alzheimer risk, "
            "requiring integrative counseling. MTHFR C677T homozygotes need "
            "methylfolate rather than folic acid."
        ),
    },

    "liver": {
        "name": "Liver Health / NAFLD / MASLD",
        "key_biomarkers": {
            "ALT": {
                "unit": "U/L",
                "normal_range": "7-56",
                "mildly_elevated": "56-100",
                "elevated": ">100",
                "clinical_note": "Primary marker of hepatocyte injury.",
            },
            "AST": {
                "unit": "U/L",
                "normal_range": "10-40",
                "elevated": ">40",
                "clinical_note": "AST:ALT ratio >2 suggests alcoholic liver disease.",
            },
            "GGT": {
                "unit": "U/L",
                "normal_range": "9-48",
                "elevated": ">48",
                "clinical_note": "Sensitive marker for biliary disease and alcohol use.",
            },
            "FIB4_index": {
                "unit": "index",
                "low_risk": "<1.3",
                "intermediate": "1.3-2.67",
                "high_risk": ">2.67",
                "clinical_note": "Non-invasive fibrosis score: (age x AST) / (platelets x sqrt(ALT)).",
            },
            "ferritin": {
                "unit": "ng/mL",
                "normal_range_male": "24-336",
                "normal_range_female": "11-307",
                "elevated": ">500",
                "clinical_note": "Elevated in iron overload, inflammation, and NAFLD.",
            },
        },
        "genetic_modifiers": [
            {"gene": "PNPLA3", "rs_id": "rs738409", "risk_allele": "G",
             "effect": "I148M variant; 2-3x NAFLD/NASH risk; strongest genetic predictor"},
            {"gene": "TM6SF2", "rs_id": "rs58542926", "risk_allele": "T",
             "effect": "E167K variant; increases hepatic steatosis but may protect against CVD"},
            {"gene": "HSD17B13", "rs_id": "rs72613567", "risk_allele": "TA",
             "effect": "Loss-of-function protects against NASH progression (therapeutic target)"},
        ],
        "intervention_targets": [
            "Weight loss 7-10% for NAFLD regression",
            "Mediterranean diet reduces hepatic steatosis",
            "Vitamin E 800 IU/day for non-diabetic NASH (PNPLA3 GG carriers benefit most)",
            "Avoid excess fructose (exacerbates PNPLA3-mediated steatosis)",
            "GLP-1 receptor agonists for NASH with diabetes",
        ],
        "clinical_context": (
            "PNPLA3 I148M is the strongest genetic determinant of NAFLD/NASH risk. "
            "Carriers of the G allele have 2-3x higher risk and progress faster to "
            "fibrosis. FIB-4 index is the recommended first-line non-invasive fibrosis "
            "assessment. The combination of PNPLA3 genotype + FIB-4 provides superior "
            "risk stratification compared to either alone."
        ),
    },

    "thyroid": {
        "name": "Thyroid Function",
        "key_biomarkers": {
            "TSH": {
                "unit": "mIU/L",
                "normal_range": "0.4-4.0",
                "subclinical_hypo": "4.0-10.0",
                "overt_hypo": ">10.0",
                "subclinical_hyper": "0.1-0.4",
                "overt_hyper": "<0.1",
                "clinical_note": "Most sensitive marker of thyroid dysfunction.",
            },
            "free_T4": {
                "unit": "ng/dL",
                "normal_range": "0.8-1.8",
                "clinical_note": "Low in hypothyroidism; high in hyperthyroidism.",
            },
            "free_T3": {
                "unit": "pg/mL",
                "normal_range": "2.3-4.2",
                "clinical_note": "Active thyroid hormone; may be normal in early hypothyroidism.",
            },
            "anti_TPO": {
                "unit": "IU/mL",
                "normal_range": "<35",
                "elevated": ">35",
                "clinical_note": "Thyroid peroxidase antibodies; marker of Hashimoto's thyroiditis.",
            },
            "reverse_T3": {
                "unit": "ng/dL",
                "normal_range": "10-24",
                "elevated": ">24",
                "clinical_note": "Elevated in euthyroid sick syndrome, chronic stress, and T4-to-T3 conversion issues.",
            },
        },
        "genetic_modifiers": [
            {"gene": "DIO2", "rs_id": "rs225014", "risk_allele": "A",
             "effect": "Thr92Ala; impaired T4-to-T3 conversion; may need combination T4/T3 therapy"},
            {"gene": "DIO1", "rs_id": "rs2235544", "risk_allele": "A",
             "effect": "Affects peripheral T4-to-T3 conversion efficiency"},
        ],
        "intervention_targets": [
            "Levothyroxine dose optimization guided by DIO2 genotype",
            "Selenium 200 mcg/day for elevated anti-TPO (reduces antibody titers)",
            "Combination T4/T3 therapy for DIO2 Thr92Ala homozygotes with persistent symptoms",
            "Iron and selenium optimization (cofactors for deiodinase enzymes)",
        ],
        "clinical_context": (
            "DIO2 Thr92Ala polymorphism affects T4-to-T3 conversion and may explain "
            "persistent hypothyroid symptoms despite normal TSH on levothyroxine. "
            "Thyroid autoimmunity (anti-TPO positive) can precede clinical "
            "hypothyroidism by years, providing an intervention window."
        ),
    },

    "iron": {
        "name": "Iron Metabolism",
        "key_biomarkers": {
            "ferritin": {
                "unit": "ng/mL",
                "normal_range_male": "24-336",
                "normal_range_female": "11-307",
                "iron_deficiency": "<15",
                "iron_overload": ">500",
                "clinical_note": "Acute phase reactant; interpret with CRP.",
            },
            "iron_saturation": {
                "unit": "%",
                "normal_range": "20-50",
                "low": "<20",
                "elevated": ">45",
                "high_risk": ">60",
                "clinical_note": "Transferrin saturation >45% prompts HFE genotyping.",
            },
            "serum_iron": {
                "unit": "ug/dL",
                "normal_range": "60-170",
                "clinical_note": "Diurnal variation; best measured fasting AM.",
            },
            "TIBC": {
                "unit": "ug/dL",
                "normal_range": "250-370",
                "clinical_note": "Total iron-binding capacity; elevated in iron deficiency.",
            },
            "hepcidin": {
                "unit": "ng/mL",
                "normal_range": "1-20",
                "elevated": ">20",
                "clinical_note": "Master regulator of iron homeostasis; elevated suppresses absorption.",
            },
        },
        "genetic_modifiers": [
            {"gene": "HFE", "rs_id": "rs1800562", "risk_allele": "A",
             "effect": "C282Y homozygosity causes hereditary hemochromatosis (iron overload)"},
            {"gene": "HFE", "rs_id": "rs1799945", "risk_allele": "G",
             "effect": "H63D; mild iron overload risk, especially compound heterozygotes C282Y/H63D"},
            {"gene": "TMPRSS6", "rs_id": "rs855791", "risk_allele": "T",
             "effect": "Affects hepcidin regulation; associated with iron-refractory iron deficiency anemia"},
        ],
        "intervention_targets": [
            "Therapeutic phlebotomy for HFE C282Y homozygotes (target ferritin 50-100 ng/mL)",
            "Iron chelation therapy (deferasirox) for transfusion-dependent iron overload",
            "IV iron (ferric carboxymaltose) for TMPRSS6-related oral iron resistance",
            "Vitamin C 200 mg with iron supplements to enhance absorption",
            "Avoid iron supplements and vitamin C supplementation in HFE carriers with elevated ferritin",
        ],
        "clinical_context": (
            "HFE C282Y homozygosity is the most common genetic cause of iron overload "
            "(hereditary hemochromatosis) in Northern European populations. Early detection "
            "via transferrin saturation screening followed by HFE genotyping prevents "
            "irreversible organ damage (cirrhosis, cardiomyopathy, diabetes). "
            "TMPRSS6 variants can cause iron deficiency that is refractory to oral iron."
        ),
    },

    "nutritional": {
        "name": "Nutritional Genomics",
        "key_biomarkers": {
            "vitamin_D_25OH": {
                "unit": "ng/mL",
                "deficient": "<20",
                "insufficient": "20-30",
                "optimal": "30-50",
                "clinical_note": "VDR polymorphisms affect vitamin D receptor sensitivity.",
            },
            "folate_serum": {
                "unit": "ng/mL",
                "deficient": "<3.0",
                "normal_range": "3.0-20.0",
                "clinical_note": "MTHFR C677T carriers may need higher levels for homocysteine control.",
            },
            "vitamin_B12": {
                "unit": "pg/mL",
                "deficient": "<200",
                "borderline": "200-300",
                "normal_range": "300-900",
                "clinical_note": "FUT2 non-secretors have higher B12 requirements.",
            },
            "omega3_index": {
                "unit": "%",
                "low_risk": ">8",
                "moderate_risk": "4-8",
                "high_risk": "<4",
                "clinical_note": "FADS1/2 variants affect omega-3 fatty acid metabolism.",
            },
            "methylmalonic_acid": {
                "unit": "nmol/L",
                "normal_range": "<270",
                "elevated": ">270",
                "clinical_note": "Functional B12 marker; more specific than serum B12.",
            },
        },
        "genetic_modifiers": [
            {"gene": "VDR", "rs_id": "rs2228570", "risk_allele": "T",
             "effect": "FokI polymorphism; affects vitamin D receptor function"},
            {"gene": "MTHFR", "rs_id": "rs1801133", "risk_allele": "T",
             "effect": "C677T; 30-70% reduced enzyme activity; need methylfolate not folic acid"},
            {"gene": "FADS1", "rs_id": "rs174537", "risk_allele": "T",
             "effect": "Reduced desaturase activity; impaired conversion of ALA to EPA/DHA"},
            {"gene": "FUT2", "rs_id": "rs602662", "risk_allele": "A",
             "effect": "Non-secretor status; reduced B12 absorption and higher requirements"},
            {"gene": "BCMO1", "rs_id": "rs7501331", "risk_allele": "T",
             "effect": "Reduced beta-carotene to vitamin A conversion by 30-70%"},
        ],
        "intervention_targets": [
            "Methylfolate (not folic acid) for MTHFR C677T carriers",
            "Preformed EPA/DHA for FADS1 poor converters (2-4 g/day)",
            "Vitamin D3 2000-5000 IU/day for VDR FokI T carriers",
            "Methylcobalamin for FUT2 non-secretors (sublingual preferred)",
            "Preformed vitamin A (retinol) for BCMO1 poor converters",
        ],
        "clinical_context": (
            "Nutritional genomics reveals why standard RDAs are insufficient for "
            "many individuals. MTHFR C677T is the most clinically relevant nutritional "
            "variant, affecting folate metabolism and methylation capacity. FADS1/2 "
            "variants are especially important for vegetarians who rely on ALA-to-EPA "
            "conversion. VDR variants affect the threshold for vitamin D sufficiency."
        ),
    },
    "kidney": {
        "name": "Kidney Function / CKD",
        "key_biomarkers": {
            "eGFR": {
                "unit": "mL/min/1.73m2",
                "normal": ">=90",
                "stage_2": "60-89",
                "stage_3": "30-59",
                "stage_4": "15-29",
                "stage_5": "<15",
                "clinical_note": "Estimated from creatinine (CKD-EPI 2021, race-neutral). Confirm with cystatin C if muscle mass is atypical.",
            },
            "cystatin_C": {
                "unit": "mg/L",
                "normal_range": "0.56-1.0",
                "elevated": ">1.0",
                "clinical_note": "Less affected by muscle mass than creatinine. Preferred for sarcopenic and muscular patients.",
            },
            "urine_ACR": {
                "unit": "mg/g",
                "normal": "<30 (A1)",
                "moderately_increased": "30-300 (A2)",
                "severely_increased": ">300 (A3)",
                "clinical_note": "Urine albumin-to-creatinine ratio; early marker of glomerular damage.",
            },
            "BUN": {
                "unit": "mg/dL",
                "normal_range": "7-20",
                "elevated": ">20",
                "clinical_note": "Blood urea nitrogen; affected by diet, hydration, and hepatic function.",
            },
        },
        "genetic_modifiers": [
            {"gene": "APOL1", "rs_id": "G1/G2 variants", "risk_allele": "G1 or G2",
             "effect": "7-10x FSGS/CKD risk with two risk alleles; critical in African ancestry"},
            {"gene": "UMOD", "rs_id": "rs12917707", "risk_allele": "G",
             "effect": "Uromodulin variant; CKD risk through tubular dysfunction"},
            {"gene": "PKD1/PKD2", "rs_id": "various", "risk_allele": "pathogenic",
             "effect": "Autosomal dominant polycystic kidney disease; screen family members"},
        ],
        "intervention_targets": [
            "SGLT2 inhibitors (dapagliflozin, empagliflozin) — slow CKD progression regardless of diabetes status",
            "ACEi/ARB therapy for albuminuria reduction and renoprotection",
            "Blood pressure target <130/80 mmHg in CKD with albuminuria",
            "Finerenone (non-steroidal MRA) for diabetic CKD with albuminuria",
            "Dietary protein restriction (0.6-0.8 g/kg/day) in CKD Stage 3b-5",
            "Avoid nephrotoxic agents (NSAIDs, aminoglycosides) — especially with APOL1 risk",
        ],
        "clinical_context": (
            "CKD affects 15% of US adults. APOL1 high-risk genotype (G1/G2 two-allele) "
            "explains the 3-4x higher CKD prevalence in African Americans compared to "
            "European Americans. SGLT2 inhibitors represent a paradigm shift, showing "
            "renoprotection independent of diabetes status (DAPA-CKD, EMPA-KIDNEY trials). "
            "Early detection via eGFR + urine ACR screening enables intervention before "
            "irreversible nephron loss."
        ),
    },
}


# =============================================================================
# 2. PHENOAGE_KNOWLEDGE -- PhenoAge clock interpretation
# =============================================================================

PHENOAGE_KNOWLEDGE: Dict[str, Any] = {
    "description": (
        "PhenoAge is a biological age estimator developed by Levine et al. (2018) "
        "using 9 routine blood biomarkers and chronological age. It was trained on "
        "NHANES III mortality data and validated against lifespan, healthspan, and "
        "disease outcomes. PhenoAge acceleration (PhenoAge - chronological age) is "
        "associated with all-cause mortality, cardiovascular disease, cancer, and "
        "cognitive decline."
    ),
    "reference": "Levine ME et al. Aging (Albany NY). 2018;10(4):573-591.",
    "biomarkers": {
        "albumin": {
            "coefficient": -0.0336,
            "unit": "g/dL",
            "direction": "protective",
            "interpretation": (
                "Higher albumin is associated with slower aging. Albumin reflects "
                "hepatic synthetic function and nutritional status. Levels <3.5 g/dL "
                "are a strong predictor of mortality."
            ),
        },
        "creatinine": {
            "coefficient": 0.0095,
            "unit": "mg/dL",
            "direction": "aging",
            "interpretation": (
                "Higher creatinine contributes to accelerated aging, reflecting "
                "declining kidney function. Must be interpreted with muscle mass."
            ),
        },
        "glucose": {
            "coefficient": 0.1953,
            "unit": "mg/dL",
            "direction": "aging",
            "interpretation": (
                "Fasting glucose has the largest coefficient. Elevated glucose drives "
                "advanced glycation end-products (AGEs), oxidative stress, and "
                "accelerated aging across multiple organ systems."
            ),
        },
        "ln_crp": {
            "coefficient": 0.0954,
            "unit": "ln(mg/L)",
            "direction": "aging",
            "interpretation": (
                "Log-transformed hs-CRP captures chronic low-grade inflammation "
                "(inflammaging). Elevated CRP accelerates biological aging and "
                "is modifiable through diet, exercise, and anti-inflammatory interventions."
            ),
        },
        "lymphocyte_pct": {
            "coefficient": -0.0120,
            "unit": "%",
            "direction": "protective",
            "interpretation": (
                "Higher lymphocyte percentage is protective. Low lymphocytes "
                "(lymphopenia) reflects immunosenescence and is associated with "
                "increased infection risk and mortality."
            ),
        },
        "mcv": {
            "coefficient": 0.0268,
            "unit": "fL",
            "direction": "aging",
            "interpretation": (
                "Higher MCV (macrocytosis) contributes to accelerated aging. "
                "Can reflect B12/folate deficiency, liver disease, or myelodysplasia."
            ),
        },
        "rdw": {
            "coefficient": 0.3306,
            "unit": "%",
            "direction": "aging",
            "interpretation": (
                "Red cell distribution width has the second-largest coefficient. "
                "Elevated RDW (anisocytosis) is one of the strongest predictors of "
                "mortality and reflects oxidative stress, chronic inflammation, and "
                "impaired erythropoiesis."
            ),
        },
        "alkaline_phosphatase": {
            "coefficient": 0.0019,
            "unit": "U/L",
            "direction": "aging",
            "interpretation": (
                "Higher ALP contributes modestly to accelerated aging. Reflects "
                "hepatobiliary or bone disease. Elevated ALP is associated with "
                "cardiovascular and all-cause mortality."
            ),
        },
        "wbc": {
            "coefficient": 0.0554,
            "unit": "10^3/uL",
            "direction": "aging",
            "interpretation": (
                "Higher WBC count reflects chronic inflammation and immune activation. "
                "Even high-normal WBC (>7.0) is associated with accelerated aging."
            ),
        },
    },
    "clinical_interpretation": {
        "acceleration_positive": (
            "Positive age acceleration (PhenoAge > chronological age) indicates "
            "biological aging faster than expected. Each year of acceleration is "
            "associated with approximately 9% increased mortality risk."
        ),
        "acceleration_negative": (
            "Negative age acceleration (PhenoAge < chronological age) indicates "
            "biological aging slower than expected, suggesting resilience and "
            "favorable metabolic health."
        ),
        "actionable_drivers": (
            "The most modifiable PhenoAge drivers are glucose (diet/exercise), "
            "hs-CRP (anti-inflammatory interventions), RDW (nutritional optimization), "
            "and albumin (protein intake). Targeting these can reduce PhenoAge "
            "by 2-5 years within 6-12 months."
        ),
    },
}


# =============================================================================
# 3. PGX_KNOWLEDGE -- 7 key pharmacogenes
# =============================================================================

PGX_KNOWLEDGE: Dict[str, Dict[str, Any]] = {
    "CYP2D6": {
        "full_name": "Cytochrome P450 2D6",
        "chromosome": "22q13.1",
        "phenotypes": {
            "ultra_rapid": "Gene duplication/amplification (*1xN, *2xN); increased metabolism",
            "normal": "*1/*1, *1/*2; standard enzyme activity",
            "intermediate": "*1/*4, *1/*5, *41/*41; reduced activity (50-75%)",
            "poor": "*4/*4, *5/*5, *4/*5; absent or minimal activity (<10%)",
        },
        "key_drugs": [
            {"drug": "codeine", "effect": "Ultra-rapid: risk of toxicity/respiratory depression; Poor: no analgesic effect",
             "recommendation": "Avoid codeine in ultra-rapid and poor metabolizers; use morphine alternative",
             "cpic_level": "1A"},
            {"drug": "tamoxifen", "effect": "Poor: reduced conversion to endoxifen; compromised efficacy",
             "recommendation": "Consider aromatase inhibitor for poor metabolizers",
             "cpic_level": "1A"},
            {"drug": "tramadol", "effect": "Ultra-rapid: toxicity risk; Poor: reduced efficacy",
             "recommendation": "Avoid in ultra-rapid and poor metabolizers",
             "cpic_level": "1A"},
            {"drug": "ondansetron", "effect": "Ultra-rapid: reduced efficacy",
             "recommendation": "Select alternative antiemetic for ultra-rapid metabolizers",
             "cpic_level": "1B"},
        ],
        "clinical_context": (
            "CYP2D6 is the most polymorphic CYP enzyme with >100 allelic variants. "
            "It metabolizes ~25% of clinically used drugs. Gene deletions (*5) and "
            "duplications (*1xN) are common. Ethnicity-based frequency: PM prevalence "
            "5-10% European, 1-2% East Asian."
        ),
    },
    "CYP2C19": {
        "full_name": "Cytochrome P450 2C19",
        "chromosome": "10q23.33",
        "phenotypes": {
            "ultra_rapid": "*17/*17; increased enzyme activity",
            "normal": "*1/*1; standard activity",
            "intermediate": "*1/*2, *1/*3; reduced activity",
            "poor": "*2/*2, *2/*3, *3/*3; absent activity",
        },
        "key_drugs": [
            {"drug": "clopidogrel", "effect": "Poor/IM: reduced active metabolite; increased CVD event risk",
             "recommendation": "Use prasugrel or ticagrelor for poor/intermediate metabolizers",
             "cpic_level": "1A"},
            {"drug": "omeprazole", "effect": "Ultra-rapid: reduced efficacy; Poor: increased drug levels",
             "recommendation": "Increase dose for ultra-rapid; reduce for poor metabolizers",
             "cpic_level": "1A"},
            {"drug": "escitalopram", "effect": "Poor: increased drug levels and side effects; Ultra-rapid: subtherapeutic",
             "recommendation": "Reduce dose 50% for poor metabolizers; consider alternative for ultra-rapid",
             "cpic_level": "1A"},
            {"drug": "voriconazole", "effect": "Poor: toxicity risk from increased exposure",
             "recommendation": "Choose alternative antifungal or reduce dose with TDM",
             "cpic_level": "1A"},
        ],
        "clinical_context": (
            "CYP2C19 is critical for clopidogrel activation. The *2 allele (loss-of-function) "
            "is found in 15% of Europeans and 30% of East Asians. *17 (gain-of-function) is "
            "present in 20-30% of Europeans. Genotyping before PCI is recommended by AHA/ACC."
        ),
    },
    "CYP2C9": {
        "full_name": "Cytochrome P450 2C9",
        "chromosome": "10q23.33",
        "phenotypes": {
            "normal": "*1/*1; standard activity",
            "intermediate": "*1/*2, *1/*3; reduced activity (50-80%)",
            "poor": "*2/*2, *3/*3, *2/*3; markedly reduced activity (<20%)",
        },
        "key_drugs": [
            {"drug": "warfarin", "effect": "Poor/IM: reduced clearance; increased bleeding risk",
             "recommendation": "Reduce initial dose 25-50%; use pharmacogenomic dosing algorithm",
             "cpic_level": "1A"},
            {"drug": "phenytoin", "effect": "Poor: increased drug levels and toxicity",
             "recommendation": "Reduce dose by 25%; monitor levels closely",
             "cpic_level": "1A"},
            {"drug": "celecoxib", "effect": "Poor: increased exposure; GI bleeding risk",
             "recommendation": "Reduce dose by 50% or use alternative NSAID",
             "cpic_level": "1B"},
        ],
        "clinical_context": (
            "CYP2C9 is critical for warfarin metabolism. Combined CYP2C9 + VKORC1 "
            "genotyping explains ~40% of warfarin dose variability. CYP2C9*2 and *3 "
            "are the most common reduced-function alleles in Europeans (10-15% each)."
        ),
    },
    "VKORC1": {
        "full_name": "Vitamin K Epoxide Reductase Complex Subunit 1",
        "chromosome": "16p11.2",
        "phenotypes": {
            "normal_sensitivity": "GG genotype; standard warfarin dose",
            "intermediate_sensitivity": "GA genotype; reduced dose needed",
            "high_sensitivity": "AA genotype; significantly reduced dose needed",
        },
        "key_drugs": [
            {"drug": "warfarin", "effect": "AA: 2-3 mg/day typical (vs 5-7 mg/day for GG)",
             "recommendation": "Use pharmacogenomic dosing algorithm combining CYP2C9 + VKORC1",
             "cpic_level": "1A"},
        ],
        "clinical_context": (
            "VKORC1 rs9923231 (-1639G>A) is the strongest genetic predictor of warfarin "
            "dose requirement, explaining ~25% of dose variability. The A allele frequency "
            "is ~40% in Europeans and ~90% in East Asians."
        ),
    },
    "SLCO1B1": {
        "full_name": "Solute Carrier Organic Anion Transporter Family Member 1B1",
        "chromosome": "12p12.1",
        "phenotypes": {
            "normal_function": "TT (rs4149056); standard statin transport",
            "intermediate_function": "TC; reduced hepatic uptake of statins",
            "poor_function": "CC; markedly reduced transport; myopathy risk",
        },
        "key_drugs": [
            {"drug": "simvastatin", "effect": "CC: 17x increased myopathy risk; avoid high doses",
             "recommendation": "Use max 20 mg simvastatin or switch to rosuvastatin/pravastatin",
             "cpic_level": "1A"},
            {"drug": "atorvastatin", "effect": "CC: increased myopathy risk (lower than simvastatin)",
             "recommendation": "Monitor for myopathy; consider lower dose",
             "cpic_level": "2A"},
        ],
        "clinical_context": (
            "SLCO1B1*5 (rs4149056 T>C) reduces hepatic statin uptake, increasing systemic "
            "exposure and myopathy risk. The C allele is present in ~15% of Europeans. "
            "The SEARCH trial showed 17x increased myopathy risk with 80 mg simvastatin "
            "in CC homozygotes."
        ),
    },
    "DPYD": {
        "full_name": "Dihydropyrimidine Dehydrogenase",
        "chromosome": "1p21.3",
        "phenotypes": {
            "normal": "*1/*1; standard 5-FU metabolism",
            "intermediate": "*1/*2A, *1/*13; reduced activity (30-70%)",
            "poor": "*2A/*2A; absent activity (<10%); life-threatening toxicity risk",
        },
        "key_drugs": [
            {"drug": "5-fluorouracil", "effect": "Poor: fatal toxicity risk (neutropenia, mucositis, death)",
             "recommendation": "Contraindicated in poor metabolizers; reduce 50% for intermediate",
             "cpic_level": "1A"},
            {"drug": "capecitabine", "effect": "Same as 5-FU (capecitabine is 5-FU prodrug)",
             "recommendation": "Pre-treatment DPYD testing is mandatory in EU (EMA 2020)",
             "cpic_level": "1A"},
        ],
        "clinical_context": (
            "DPYD deficiency causes life-threatening fluoropyrimidine toxicity. "
            "Pre-treatment testing is mandatory in the EU and recommended by CPIC. "
            "DPYD*2A (IVS14+1G>A) completely eliminates enzyme activity. "
            "Approximately 3-5% of the population carries a reduced-function allele."
        ),
    },
    "CYP3A5": {
        "full_name": "Cytochrome P450 3A5",
        "chromosome": "7q22.1",
        "phenotypes": {
            "expresser": "*1/*1 or *1/*3; functional enzyme (more common in African descent)",
            "non_expresser": "*3/*3; absent expression (most Europeans)",
        },
        "key_drugs": [
            {"drug": "tacrolimus", "effect": "Expressers: rapid metabolism; need higher starting dose",
             "recommendation": "Start 1.5-2x standard dose for *1/*1; titrate with TDM",
             "cpic_level": "1A"},
        ],
        "clinical_context": (
            "CYP3A5 expression is highly ethnicity-dependent: *1 (functional) allele "
            "frequency is 60-90% in African populations vs 5-15% in Europeans. This "
            "explains racial disparities in tacrolimus dosing requirements after "
            "organ transplantation."
        ),
    },
    "UGT1A1": {
        "full_name": "UDP-Glucuronosyltransferase 1A1",
        "chromosome": "2q37.1",
        "phenotypes": {
            "normal": "*1/*1; standard glucuronidation activity",
            "intermediate": "*1/*28 or *1/*6; reduced activity (Gilbert's-associated)",
            "poor": "*28/*28 or *6/*6; markedly reduced activity",
        },
        "key_drugs": [
            {"drug": "irinotecan", "effect": "Poor: 3-4x increased severe neutropenia and diarrhea risk",
             "recommendation": "Reduce dose by 30-50% for *28/*28; FDA label includes UGT1A1*28 dosing guidance",
             "cpic_level": "1A"},
            {"drug": "atazanavir", "effect": "Poor: hyperbilirubinemia (cosmetic jaundice); may affect adherence",
             "recommendation": "Consider alternative PI if jaundice is unacceptable to patient",
             "cpic_level": "1B"},
        ],
        "clinical_context": (
            "UGT1A1*28 (TA7/TA7) is the genetic basis of Gilbert's syndrome, present in "
            "~10% of European and ~30-40% of African populations. While benign in isolation, "
            "it becomes critically important when prescribing irinotecan, where impaired "
            "glucuronidation of the active metabolite SN-38 can cause life-threatening toxicity. "
            "UGT1A1*6 is more common in East Asian populations."
        ),
    },
    "NUDT15": {
        "full_name": "Nudix Hydrolase 15",
        "chromosome": "13q14.2",
        "phenotypes": {
            "normal": "*1/*1; standard thiopurine metabolism",
            "intermediate": "*1/*2, *1/*3; reduced nucleotide pool regulation",
            "poor": "*2/*2, *3/*3; absent activity — severe myelosuppression risk",
        },
        "key_drugs": [
            {"drug": "azathioprine", "effect": "Poor: severe myelosuppression risk; critical for East Asian patients",
             "recommendation": "Reduce dose by 75-90% or avoid; *3 frequency is 7-10% in East Asians",
             "cpic_level": "1A"},
            {"drug": "6-mercaptopurine", "effect": "Same as azathioprine — shared thiopurine pathway",
             "recommendation": "Combined TPMT + NUDT15 testing recommended before thiopurines",
             "cpic_level": "1A"},
        ],
        "clinical_context": (
            "NUDT15 was identified as a major determinant of thiopurine toxicity in East Asian "
            "populations, where TPMT variants are rare. NUDT15*3 (p.Arg139Cys) has a frequency "
            "of 7-10% in East Asians but <1% in Europeans. Combined TPMT + NUDT15 testing "
            "is now recommended by CPIC for all patients before starting thiopurines."
        ),
    },
}


# =============================================================================
# 4. CROSS_MODAL_LINKS -- Triggers to other HCLS AI Factory agents
# =============================================================================

CROSS_MODAL_LINKS: Dict[str, Dict[str, Any]] = {
    "elevated_lpa": {
        "trigger_condition": "Lp(a) > 125 nmol/L",
        "target_agent": "imaging_intelligence_agent",
        "target_url": "not_deployed",
        "action": "cardiovascular_imaging_assessment",
        "description": (
            "Elevated Lp(a) is the strongest heritable cardiovascular risk factor. "
            "When detected, trigger the Imaging Intelligence Agent to assess coronary "
            "artery calcium score (CAC), carotid intima-media thickness (CIMT), and "
            "cardiac CT angiography for subclinical atherosclerosis."
        ),
        "urgency": "high",
        "payload": {
            "biomarker": "Lp(a)",
            "threshold_exceeded": ">125 nmol/L",
            "imaging_modalities": ["CAC_score", "CIMT", "CCTA"],
            "clinical_rationale": "Screen for subclinical atherosclerosis in high Lp(a)",
        },
    },
    "pgx_critical": {
        "trigger_condition": "Critical PGx finding (DPYD poor metabolizer or CYP2D6 ultra-rapid on opioids)",
        "target_agent": "cart_intelligence_agent",
        "target_url": "http://localhost:8525/api/v1/drug-interaction-check",
        "action": "drug_interaction_check",
        "description": (
            "When a critical pharmacogenomic finding is detected (e.g., DPYD poor "
            "metabolizer on 5-FU, CYP2D6 ultra-rapid on codeine), trigger the "
            "CAR-T/Oncology Intelligence Agent to check for drug interaction conflicts "
            "in any active treatment protocols."
        ),
        "urgency": "critical",
        "payload": {
            "pgx_gene": "variable",
            "metabolizer_status": "variable",
            "flagged_drugs": [],
            "clinical_rationale": "Cross-check oncology treatment protocols for PGx conflicts",
        },
    },
    "pre_diabetic_trajectory": {
        "trigger_condition": "HOMA-IR > 2.5 AND HbA1c > 5.5% AND TCF7L2 risk carrier",
        "target_agent": "genomics_pipeline",
        "target_url": "http://localhost:5001/api/v1/reanalyze",
        "action": "vcf_reanalysis_diabetes_genes",
        "description": (
            "When biomarker trajectory analysis detects pre-diabetic pattern combined "
            "with TCF7L2 risk genotype, trigger the Genomics Pipeline to re-analyze "
            "the patient VCF for additional diabetes-associated rare variants in "
            "HNF1A, HNF4A, GCK, and other MODY genes."
        ),
        "urgency": "moderate",
        "payload": {
            "gene_panel": ["TCF7L2", "HNF1A", "HNF4A", "GCK", "SLC30A8", "KCNJ11", "ABCC8"],
            "analysis_type": "rare_variant_screening",
            "clinical_rationale": "Screen for monogenic diabetes (MODY) in pre-diabetic trajectory",
        },
    },
    "iron_overload_detected": {
        "trigger_condition": "Ferritin > 500 AND transferrin saturation > 45%",
        "target_agent": "imaging_intelligence_agent",
        "target_url": "not_deployed",
        "action": "liver_iron_quantification",
        "description": (
            "Iron overload detected by biomarkers. Trigger imaging assessment for "
            "hepatic iron concentration (HIC) via MRI R2* and liver stiffness "
            "measurement via FibroScan to assess fibrosis stage."
        ),
        "urgency": "high",
        "payload": {
            "biomarkers": {"ferritin": ">500", "transferrin_saturation": ">45%"},
            "imaging_modalities": ["MRI_R2star", "FibroScan"],
            "hfe_genotype": "check_required",
            "clinical_rationale": "Quantify hepatic iron and assess fibrosis in iron overload",
        },
    },
    "accelerated_aging_detected": {
        "trigger_condition": "PhenoAge acceleration > 5 years",
        "target_agent": "genomics_pipeline",
        "target_url": "http://localhost:5001/api/v1/reanalyze",
        "action": "epigenetic_variant_analysis",
        "description": (
            "Significant biological age acceleration detected. Trigger genomic pipeline "
            "to analyze telomere-related variants (TERT, TERC) and DNA repair genes "
            "(ATM, BRCA1) for germline variants that may explain accelerated aging."
        ),
        "urgency": "moderate",
        "payload": {
            "gene_panel": ["TERT", "TERC", "ATM", "BRCA1", "LMNA", "WRN"],
            "analysis_type": "aging_variant_screening",
            "clinical_rationale": "Screen for germline aging-associated variants",
        },
    },
    "ckd_detected": {
        "trigger_condition": "eGFR < 60 OR urine ACR > 300",
        "target_agent": "imaging_intelligence_agent",
        "target_url": "not_deployed",
        "action": "renal_imaging_assessment",
        "description": (
            "CKD Stage 3+ or severely increased albuminuria detected. Trigger imaging "
            "assessment for renal ultrasound (kidney size, cortical thickness, obstruction), "
            "and Doppler for renal artery stenosis if clinical suspicion."
        ),
        "urgency": "high",
        "payload": {
            "biomarkers": {"eGFR": "variable", "urine_ACR": "variable"},
            "imaging_modalities": ["renal_ultrasound", "renal_doppler"],
            "clinical_rationale": "Structural assessment of kidneys in CKD detection",
        },
    },
    "cognitive_risk_detected": {
        "trigger_condition": "APOE E4 carrier AND (homocysteine > 15 OR omega3_index < 6)",
        "target_agent": "imaging_intelligence_agent",
        "target_url": "not_deployed",
        "action": "brain_imaging_assessment",
        "description": (
            "Elevated Alzheimer's risk with APOE E4 genotype and modifiable risk factors "
            "detected. Trigger imaging assessment for brain MRI volumetry (hippocampal "
            "volume), white matter hyperintensities, and optional amyloid PET if indicated."
        ),
        "urgency": "moderate",
        "payload": {
            "genetic_risk": "APOE_E4",
            "imaging_modalities": ["brain_MRI_volumetry", "WMH_assessment"],
            "clinical_rationale": "Screen for neurodegeneration biomarkers in high-risk genotype",
        },
    },
    "bone_health_alert": {
        "trigger_condition": "Vitamin D < 20 AND PTH > 65 AND (VDR risk OR COL1A1 risk)",
        "target_agent": "imaging_intelligence_agent",
        "target_url": "not_deployed",
        "action": "bone_density_assessment",
        "description": (
            "Osteoporosis risk detected with vitamin D deficiency, secondary "
            "hyperparathyroidism, and genetic bone risk factors. Trigger imaging "
            "for DEXA scan and vertebral fracture assessment."
        ),
        "urgency": "moderate",
        "payload": {
            "biomarkers": {"vitamin_D": "deficient", "PTH": "elevated"},
            "imaging_modalities": ["DEXA_scan", "vertebral_fracture_assessment"],
            "clinical_rationale": "Bone density screening with genetic + biochemical risk factors",
        },
    },
}


# =============================================================================
# 5. HELPER FUNCTIONS
# =============================================================================


def get_domain_context(domain: str) -> str:
    """Return formatted knowledge context for a disease domain.

    Args:
        domain: One of 'diabetes', 'cardiovascular', 'liver', 'thyroid', 'iron', 'nutritional', 'kidney'.

    Returns:
        Formatted string with domain knowledge or empty string if not found.
    """
    info = BIOMARKER_DOMAINS.get(domain.lower())
    if not info:
        return ""

    lines = [f"Domain: {info['name']}"]

    # Key biomarkers
    lines.append("Key Biomarkers:")
    for name, details in info["key_biomarkers"].items():
        unit = details.get("unit", "")
        normal = details.get("normal_range", details.get("normal_range_male", ""))
        note = details.get("clinical_note", "")
        lines.append(f"  - {name} ({unit}): normal {normal}. {note}")

    # Genetic modifiers
    if info.get("genetic_modifiers"):
        lines.append("Genetic Modifiers:")
        for mod in info["genetic_modifiers"]:
            lines.append(
                f"  - {mod['gene']} {mod['rs_id']} ({mod['risk_allele']}): {mod['effect']}"
            )

    # Intervention targets
    if info.get("intervention_targets"):
        lines.append("Intervention Targets:")
        for target in info["intervention_targets"]:
            lines.append(f"  - {target}")

    # Clinical context
    if info.get("clinical_context"):
        lines.append(f"Clinical Context: {info['clinical_context']}")

    return "\n".join(lines)


def get_pgx_context(gene: str) -> str:
    """Return formatted PGx knowledge context for a pharmacogene.

    Args:
        gene: Gene symbol (e.g., 'CYP2D6', 'DPYD').

    Returns:
        Formatted string with PGx knowledge or empty string if not found.
    """
    gene_upper = gene.upper()
    info = PGX_KNOWLEDGE.get(gene_upper)
    if not info:
        # Try case-insensitive lookup
        for key, val in PGX_KNOWLEDGE.items():
            if key.upper() == gene_upper:
                info = val
                break
    if not info:
        return ""

    lines = [f"Pharmacogene: {info['full_name']} ({gene_upper})"]

    # Phenotypes
    lines.append("Phenotypes:")
    for phenotype, desc in info.get("phenotypes", {}).items():
        lines.append(f"  - {phenotype}: {desc}")

    # Key drugs
    lines.append("Key Drug Interactions:")
    for drug in info.get("key_drugs", []):
        lines.append(
            f"  - {drug['drug']}: {drug['effect']} "
            f"(CPIC Level {drug['cpic_level']})"
        )
        lines.append(f"    Recommendation: {drug['recommendation']}")

    # Clinical context
    if info.get("clinical_context"):
        lines.append(f"Clinical Context: {info['clinical_context']}")

    return "\n".join(lines)


def get_biomarker_context(biomarker_name: str) -> str:
    """Return formatted knowledge context for a specific biomarker.

    Searches across all disease domains for the biomarker.

    Args:
        biomarker_name: Biomarker name (e.g., 'HbA1c', 'Lp(a)', 'ferritin').

    Returns:
        Formatted string with biomarker knowledge or empty string if not found.
    """
    name_lower = biomarker_name.lower().replace(" ", "_")

    # Search all domains
    for domain_key, domain_info in BIOMARKER_DOMAINS.items():
        for marker_key, marker_details in domain_info["key_biomarkers"].items():
            if marker_key.lower() == name_lower or biomarker_name.lower() in marker_key.lower():
                lines = [
                    f"Biomarker: {marker_key}",
                    f"Domain: {domain_info['name']}",
                    f"Unit: {marker_details.get('unit', 'N/A')}",
                ]
                # Add all range info
                for rk, rv in marker_details.items():
                    if rk not in ("unit", "clinical_note"):
                        lines.append(f"  {rk}: {rv}")
                if marker_details.get("clinical_note"):
                    lines.append(f"Clinical Note: {marker_details['clinical_note']}")

                # Add relevant genetic modifiers from same domain
                modifiers = [
                    m for m in domain_info.get("genetic_modifiers", [])
                ]
                if modifiers:
                    lines.append("Related Genetic Modifiers:")
                    for mod in modifiers:
                        lines.append(f"  - {mod['gene']} {mod['rs_id']}: {mod['effect']}")

                return "\n".join(lines)

    # Check PhenoAge biomarkers
    for marker_key, marker_details in PHENOAGE_KNOWLEDGE.get("biomarkers", {}).items():
        if marker_key.lower() == name_lower or biomarker_name.lower() in marker_key.lower():
            return (
                f"Biomarker: {marker_key} (PhenoAge Clock)\n"
                f"Coefficient: {marker_details['coefficient']}\n"
                f"Unit: {marker_details['unit']}\n"
                f"Direction: {marker_details['direction']}\n"
                f"Interpretation: {marker_details['interpretation']}"
            )

    return ""


def get_cross_modal_context(trigger_key: str) -> str:
    """Return formatted cross-modal link context.

    Args:
        trigger_key: Key from CROSS_MODAL_LINKS (e.g., 'elevated_lpa').

    Returns:
        Formatted string with cross-modal link details.
    """
    link = CROSS_MODAL_LINKS.get(trigger_key)
    if not link:
        return ""

    return (
        f"Cross-Modal Trigger: {trigger_key}\n"
        f"Condition: {link['trigger_condition']}\n"
        f"Target Agent: {link['target_agent']}\n"
        f"Action: {link['action']}\n"
        f"Urgency: {link['urgency']}\n"
        f"Description: {link['description']}"
    )


def get_knowledge_stats() -> Dict[str, int]:
    """Return statistics about the biomarker knowledge graph.

    Returns:
        Dict with counts of domains, biomarkers, genetic modifiers,
        pharmacogenes, PhenoAge markers, and cross-modal links.
    """
    total_biomarkers = sum(
        len(d["key_biomarkers"]) for d in BIOMARKER_DOMAINS.values()
    )
    total_genetic_modifiers = sum(
        len(d.get("genetic_modifiers", [])) for d in BIOMARKER_DOMAINS.values()
    )
    total_pgx_drugs = sum(
        len(g.get("key_drugs", [])) for g in PGX_KNOWLEDGE.values()
    )

    return {
        "disease_domains": len(BIOMARKER_DOMAINS),
        "total_biomarkers": total_biomarkers,
        "total_genetic_modifiers": total_genetic_modifiers,
        "pharmacogenes": len(PGX_KNOWLEDGE),
        "pgx_drug_interactions": total_pgx_drugs,
        "phenoage_markers": len(PHENOAGE_KNOWLEDGE.get("biomarkers", {})),
        "cross_modal_links": len(CROSS_MODAL_LINKS),
    }
