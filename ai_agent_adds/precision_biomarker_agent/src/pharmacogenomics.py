"""Pharmacogenomic mapping engine.

Maps star alleles and genotypes to drug recommendations following CPIC
Level 1A guidelines. Pure computation — no LLM or database calls.

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List, Optional

from loguru import logger

# CPIC guideline versions used for each gene
# Update when new CPIC publications change recommendations
CPIC_GUIDELINE_VERSIONS = {
    "CYP2D6": {"version": "2019", "pmid": "33387367", "update": "2020-12", "level": "1A"},
    "CYP2C19": {"version": "2022", "pmid": "34697867", "update": "2022-12", "level": "1A"},
    "CYP2C9": {"version": "2017", "pmid": "27441996", "update": "2020-01", "level": "1A"},
    "VKORC1": {"version": "2017", "pmid": "27441996", "update": "2020-01", "level": "1A"},
    "SLCO1B1": {"version": "2022", "pmid": "35152405", "update": "2022-06", "level": "1A"},
    "TPMT": {"version": "2018", "pmid": "30447069", "update": "2023-03", "level": "1A"},
    "DPYD": {"version": "2017", "pmid": "29152729", "update": "2023-12", "level": "1A"},
    "MTHFR": {"version": "N/A", "pmid": "N/A", "update": "N/A", "level": "Informational"},
    "HLA-B*57:01": {"version": "2014", "pmid": "24561393", "update": "2014-01", "level": "1A"},
    "G6PD": {"version": "N/A", "pmid": "N/A", "update": "N/A", "level": "Informational"},
    "HLA-B*58:01": {"version": "2015", "pmid": "23232549", "update": "2015-01", "level": "1A"},
    "CYP3A5": {"version": "2015", "pmid": "25801146", "update": "2022-11", "level": "1A"},
    "UGT1A1": {"version": "2020", "pmid": "32233013", "update": "2020-04", "level": "1A"},
    "NUDT15": {"version": "2018", "pmid": "30447069", "update": "2023-03", "level": "1A"},
}


# ---------------------------------------------------------------------------
# PGx gene configurations: star alleles -> phenotype -> drug recommendations
# All mappings follow CPIC Level 1A guidelines
#
# Phenotype naming convention:
#   - CYP enzymes (CYP2D6, CYP2C19, CYP2C9, TPMT, DPYD): use CPIC standard
#     metabolizer terms: "Normal Metabolizer", "Intermediate Metabolizer",
#     "Poor Metabolizer", "Ultra-rapid Metabolizer", "Rapid Metabolizer"
#   - SLCO1B1: uses transporter function terms ("Normal Function",
#     "Intermediate Function", "Poor Function") -- transport activity, not metabolism
#   - MTHFR: uses enzyme activity terms ("Normal Activity",
#     "Intermediate Activity", "Reduced Activity") -- folate enzyme activity
#   - VKORC1: uses drug sensitivity terms ("Normal Sensitivity",
#     "Intermediate Sensitivity", "High Sensitivity") -- warfarin target sensitivity
#   - HLA genes (HLA-B*57:01, HLA-B*58:01): "Negative" / "Positive"
#   - G6PD: "Normal" / "Intermediate" / "Deficient"
#
# Drug recommendation format:
#   Each drug entry maps phenotype -> {recommendation, action, alert_level}
#   - recommendation: Free-text clinical recommendation string
#   - action: Clinical decision support category, one of:
#       STANDARD_DOSING   -- no change needed
#       DOSE_REDUCTION    -- reduce dose per recommendation
#       DOSE_ADJUSTMENT   -- adjust dose (up or down) based on context
#       CONSIDER_ALTERNATIVE -- current drug may work but alternative preferred
#       AVOID             -- do not use this drug
#       CONTRAINDICATED   -- absolute contraindication (FDA/EMA mandated)
#   - alert_level: INFO (routine), WARNING (clinical review), CRITICAL (immediate action)
# ---------------------------------------------------------------------------

PGX_GENE_CONFIGS = {
    "CYP2D6": {
        "display_name": "CYP2D6",
        "description": "Cytochrome P450 2D6 — metabolizes ~25% of drugs",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Normal Metabolizer",
            "*2/*2": "Normal Metabolizer",
            "*1/*4": "Intermediate Metabolizer",
            "*1/*5": "Intermediate Metabolizer",
            "*1/*10": "Intermediate Metabolizer",
            "*2/*4": "Intermediate Metabolizer",
            "*4/*10": "Intermediate Metabolizer",
            "*10/*10": "Intermediate Metabolizer",
            "*4/*4": "Poor Metabolizer",
            "*4/*5": "Poor Metabolizer",
            "*5/*5": "Poor Metabolizer",
            "*4/*6": "Poor Metabolizer",
            "*1/*1xN": "Ultra-rapid Metabolizer",
            "*1/*2xN": "Ultra-rapid Metabolizer",
            "*2/*2xN": "Ultra-rapid Metabolizer",
        },
        "drug_recommendations": {
            "codeine": {
                "Normal Metabolizer": {
                    "recommendation": "Use codeine per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use codeine with caution; consider alternative analgesics (morphine, non-opioid). Reduced conversion to morphine expected",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "AVOID codeine — no conversion to morphine, will be ineffective. Use morphine, non-opioid analgesics, or non-tramadol opioids",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "AVOID codeine — excess conversion to morphine, risk of fatal respiratory depression. Use non-opioid analgesics or non-tramadol alternatives",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "tramadol": {
                "Normal Metabolizer": {
                    "recommendation": "Use tramadol per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use tramadol with caution; may have reduced efficacy",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "AVOID tramadol — reduced efficacy due to impaired O-demethylation. Use alternative analgesics",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "AVOID tramadol — risk of respiratory depression from excess O-desmethyltramadol. Use non-opioid analgesics",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "tamoxifen": {
                "Normal Metabolizer": {
                    "recommendation": "Use tamoxifen per standard dosing (20 mg/day)",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Consider higher tamoxifen dose (40 mg/day) or switch to aromatase inhibitor if postmenopausal",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "AVOID tamoxifen — greatly reduced conversion to active endoxifen. Use aromatase inhibitor if postmenopausal; consider alternative endocrine therapy if premenopausal",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "Use tamoxifen per standard dosing; may have increased active metabolite levels",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
            },
            "ondansetron": {
                "Normal Metabolizer": {
                    "recommendation": "Use ondansetron per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use ondansetron per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Ondansetron may have increased exposure. Consider dose reduction or alternative antiemetic (granisetron)",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "Ondansetron may have reduced efficacy. Consider alternative antiemetic or higher dose",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
            },
        },
    },
    "CYP2C19": {
        "display_name": "CYP2C19",
        "description": "Cytochrome P450 2C19 — critical for clopidogrel activation",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*2/*17": "Intermediate Metabolizer",
            "*2/*2": "Poor Metabolizer",
            "*2/*3": "Poor Metabolizer",
            "*3/*3": "Poor Metabolizer",
            "*1/*17": "Rapid Metabolizer",
            "*17/*17": "Ultra-rapid Metabolizer",
        },
        "drug_recommendations": {
            "clopidogrel": {
                "Normal Metabolizer": {
                    "recommendation": "Use clopidogrel per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Clopidogrel may have reduced efficacy. Consider alternative antiplatelet (prasugrel or ticagrelor) especially for ACS/PCI",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "AVOID clopidogrel — significantly reduced platelet inhibition. Use prasugrel or ticagrelor (FDA boxed warning)",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
                "Rapid Metabolizer": {
                    "recommendation": "Use clopidogrel per standard dosing; expected normal response",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "Use clopidogrel per standard dosing; may have enhanced response",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
            },
            "citalopram": {
                "Normal Metabolizer": {
                    "recommendation": "Use citalopram per standard dosing (20-40 mg/day)",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use citalopram per standard dosing; monitor for side effects",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce citalopram dose by 50% (max 20 mg/day). QTc prolongation risk at standard doses. Consider alternative SSRI (sertraline)",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Rapid Metabolizer": {
                    "recommendation": "Use citalopram per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "Citalopram may have reduced efficacy. Consider alternative SSRI or dose increase with monitoring",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
            },
            "escitalopram": {
                "Normal Metabolizer": {
                    "recommendation": "Use escitalopram per standard dosing (10-20 mg/day)",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use escitalopram per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce escitalopram dose by 50% (max 10 mg/day). Consider alternative SSRI (sertraline)",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Rapid Metabolizer": {
                    "recommendation": "Use escitalopram per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "Escitalopram may have reduced efficacy. Consider alternative SSRI or dose increase",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
            },
            "omeprazole": {
                "Normal Metabolizer": {
                    "recommendation": "Use omeprazole per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use omeprazole per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Omeprazole exposure increased 4-12x. Consider 50% dose reduction for long-term use. Monitor for adverse effects",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Rapid Metabolizer": {
                    "recommendation": "May need increased omeprazole dose for H. pylori eradication",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "INFO",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "Increase omeprazole dose (2-3x standard) or switch to rabeprazole (less CYP2C19 dependent)",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
            },
            "voriconazole": {
                "Normal Metabolizer": {
                    "recommendation": "Use voriconazole per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use voriconazole per standard dosing; therapeutic drug monitoring recommended",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Voriconazole exposure significantly increased. Reduce dose and perform therapeutic drug monitoring. Consider alternative antifungal",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Rapid Metabolizer": {
                    "recommendation": "Voriconazole may be subtherapeutic. Consider alternative antifungal or therapeutic drug monitoring with dose adjustment",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Ultra-rapid Metabolizer": {
                    "recommendation": "AVOID voriconazole — likely subtherapeutic levels. Use alternative antifungal (isavuconazole, posaconazole)",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    "CYP2C9": {
        "display_name": "CYP2C9",
        "description": "Cytochrome P450 2C9 — metabolizes warfarin, phenytoin, NSAIDs",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*2/*2": "Poor Metabolizer",
            "*2/*3": "Poor Metabolizer",
            "*3/*3": "Poor Metabolizer",
        },
        "drug_recommendations": {
            "warfarin": {
                "Normal Metabolizer": {
                    "recommendation": "Use standard warfarin dosing algorithm with VKORC1 genotype",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce warfarin dose by 20-30%. CYP2C9 intermediate metabolizers clear S-warfarin more slowly. Use pharmacogenomic dosing algorithm with VKORC1",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce warfarin dose by 30-50%. CYP2C9 poor metabolizers have significantly delayed S-warfarin clearance. High bleeding risk. Use pharmacogenomic dosing algorithm",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "CRITICAL",
                },
            },
            "phenytoin": {
                "Normal Metabolizer": {
                    "recommendation": "Use phenytoin per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce phenytoin dose by 25%. Monitor levels closely — CYP2C9 IMs have reduced phenytoin clearance",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce phenytoin dose by 50%. CYP2C9 PMs have markedly reduced clearance with toxicity risk (ataxia, nystagmus). Consider alternative anticonvulsant",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "CRITICAL",
                },
            },
            "celecoxib": {
                "Normal Metabolizer": {
                    "recommendation": "Use celecoxib per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use celecoxib with caution at lowest effective dose",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce celecoxib starting dose by 50%. CYP2C9 PMs have increased exposure. Consider alternative NSAID with less CYP2C9 dependence",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
            },
            "ibuprofen": {
                "Normal Metabolizer": {
                    "recommendation": "Use ibuprofen per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use lowest effective ibuprofen dose. CYP2C9 IMs have reduced NSAID clearance with increased GI bleeding risk",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce ibuprofen dose by 25-50% or use alternative analgesic (acetaminophen). CYP2C9 PMs have increased NSAID exposure and GI/renal toxicity risk",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
            },
            "flurbiprofen": {
                "Normal Metabolizer": {
                    "recommendation": "Use flurbiprofen per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use lowest effective flurbiprofen dose. Monitor for GI side effects",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce flurbiprofen dose by 50% or use alternative. CYP2C9 PMs have significantly increased flurbiprofen exposure",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
            },
            "losartan": {
                "Normal Metabolizer": {
                    "recommendation": "Use losartan per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Losartan may have reduced efficacy — CYP2C9 converts losartan to active E3174 metabolite. Monitor blood pressure closely. Consider alternative ARB (valsartan, irbesartan)",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Losartan likely ineffective — CYP2C9 PMs produce minimal active metabolite E3174. Switch to non-CYP2C9-dependent ARB (valsartan, irbesartan, candesartan)",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    "SLCO1B1": {
        "display_name": "SLCO1B1",
        "description": "Solute carrier organic anion transporter — statin hepatic uptake",
        "genotype_field": "SLCO1B1_rs4149056",
        "genotype_to_phenotype": {
            "TT": "Normal Function",
            "TC": "Intermediate Function",
            "CC": "Poor Function",
        },
        "drug_recommendations": {
            "simvastatin": {
                "Normal Function": {
                    "recommendation": "Use simvastatin per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Function": {
                    "recommendation": "Simvastatin max 20 mg/day. AVOID 80mg dose. Consider switching to rosuvastatin or pitavastatin (less SLCO1B1 dependent)",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Function": {
                    "recommendation": "AVOID simvastatin — high risk of myopathy/rhabdomyolysis. Use rosuvastatin or pitavastatin instead",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "atorvastatin": {
                "Normal Function": {
                    "recommendation": "Use atorvastatin per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Function": {
                    "recommendation": "Use atorvastatin with caution; lower starting dose. Monitor for myopathy symptoms",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
                "Poor Function": {
                    "recommendation": "Atorvastatin at reduced dose with close monitoring, or switch to rosuvastatin/pitavastatin. Monitor CK levels",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
            },
            "rosuvastatin": {
                "Normal Function": {
                    "recommendation": "Use rosuvastatin per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Function": {
                    "recommendation": "Use rosuvastatin per standard dosing; less affected by SLCO1B1 variants",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Function": {
                    "recommendation": "Use rosuvastatin at standard dose; preferred statin for SLCO1B1 poor function. Alternatively consider pitavastatin",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
            },
        },
    },
    "VKORC1": {
        "display_name": "VKORC1",
        "description": "Vitamin K epoxide reductase — warfarin target",
        "genotype_field": "VKORC1_rs9923231",
        "genotype_to_phenotype": {
            "GG": "Normal Sensitivity",
            "AG": "Intermediate Sensitivity",
            "GA": "Intermediate Sensitivity",
            "AA": "High Sensitivity",
        },
        "drug_recommendations": {
            "warfarin": {
                "Normal Sensitivity": {
                    "recommendation": "Standard warfarin dosing range (5-7 mg/day). Titrate based on INR",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Sensitivity": {
                    "recommendation": "Reduce warfarin dose (3-4 mg/day expected maintenance). VKORC1 AG: ~25% dose reduction. Titrate based on INR",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "High Sensitivity": {
                    "recommendation": "Low warfarin dose required (1-2 mg/day expected maintenance). VKORC1 AA: ~50% dose reduction. High bleeding risk — close INR monitoring essential",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
            },
        },
    },
    "MTHFR": {
        "display_name": "MTHFR",
        "description": "Methylenetetrahydrofolate reductase — folate metabolism",
        "genotype_field": "MTHFR_rs1801133",
        "genotype_to_phenotype": {
            "CC": "Normal Activity",
            "CT": "Intermediate Activity",
            "TC": "Intermediate Activity",
            "TT": "Reduced Activity",
        },
        "drug_recommendations": {
            "methotrexate": {
                "Normal Activity": {
                    "recommendation": "Use methotrexate per standard dosing with standard folate supplementation",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Activity": {
                    "recommendation": "Monitor for methotrexate toxicity. Consider dose adjustment if adverse effects occur. Ensure adequate folate supplementation",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
                "Reduced Activity": {
                    "recommendation": "5-7x increased methotrexate toxicity risk. Dose reduction recommended. Use leucovorin rescue. Supplement with L-methylfolate (NOT folic acid). Close CBC monitoring",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "CRITICAL",
                },
            },
            "folate_supplementation": {
                "Normal Activity": {
                    "recommendation": "Standard folic acid supplementation is effective",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Activity": {
                    "recommendation": "Consider L-methylfolate instead of folic acid for optimal folate status",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "INFO",
                },
                "Reduced Activity": {
                    "recommendation": "Use L-methylfolate (1-5 mg/day) instead of folic acid. Folic acid cannot be efficiently converted to active form (5-MTHF). Unmetabolized folic acid may accumulate",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
            },
        },
    },
    "TPMT": {
        "display_name": "TPMT",
        "description": "Thiopurine S-methyltransferase — thiopurine metabolism",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*3A": "Intermediate Metabolizer",
            "*1/*3C": "Intermediate Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*3A/*3A": "Poor Metabolizer",
            "*3A/*3C": "Poor Metabolizer",
            "*3C/*3C": "Poor Metabolizer",
            "*2/*3A": "Poor Metabolizer",
        },
        "drug_recommendations": {
            "azathioprine": {
                "Normal Metabolizer": {
                    "recommendation": "Use azathioprine per standard dosing (2-3 mg/kg/day)",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce azathioprine dose by 50% (1-1.5 mg/kg/day). Monitor CBC weekly for first month, then monthly",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce azathioprine to 10-15% of standard dose or AVOID. Severe myelosuppression risk (life-threatening). If used, start at 0.5 mg/kg/day with intensive CBC monitoring",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "6-mercaptopurine": {
                "Normal Metabolizer": {
                    "recommendation": "Use 6-mercaptopurine per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce 6-mercaptopurine dose by 50%. Monitor CBC closely",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce 6-mercaptopurine to 10-15% of standard dose. Severe myelosuppression risk. Consider alternative therapy",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "thioguanine": {
                "Normal Metabolizer": {
                    "recommendation": "Use thioguanine per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce thioguanine dose by 50%. Monitor for myelosuppression",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce thioguanine to 10-15% of standard dose. Severe myelosuppression risk",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    "HLA-B*57:01": {
        "display_name": "HLA-B*57:01",
        "description": "Human leukocyte antigen B*57:01 — abacavir hypersensitivity",
        "genotype_field": "HLA_B5701",
        "genotype_to_phenotype": {
            "negative": "Negative",
            "Negative": "Negative",
            "positive": "Positive",
            "Positive": "Positive",
            "carrier": "Positive",
            "Carrier": "Positive",
        },
        "drug_recommendations": {
            "abacavir": {
                "Negative": {
                    "recommendation": "Abacavir use is permitted. HLA-B*57:01 negative — low risk of hypersensitivity reaction",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Positive": {
                    "recommendation": "CONTRAINDICATED: Do NOT prescribe abacavir. HLA-B*57:01 positive — 5-8% risk of severe, potentially fatal hypersensitivity reaction. FDA-mandated screening. Use alternative NRTI (tenofovir)",
                    "action": "CONTRAINDICATED",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    # -----------------------------------------------------------------------
    # DPYD — CPIC Level 1A for fluoropyrimidines (5-FU, capecitabine)
    # Key variant: *2A (c.1905+1G>A, rs3918290) — most common pathogenic allele
    # CPIC Guideline: Amstutz et al., Clin Pharmacol Ther, 2018
    # NOTE for agent.py: Add critical alert for DPYD Poor Metabolizer —
    #   fluoropyrimidines are CONTRAINDICATED (life-threatening toxicity).
    #   Alert text: "DPYD *2A/*2A: 5-FU and capecitabine are CONTRAINDICATED
    #   — risk of lethal DPD deficiency toxicity (mucositis, myelosuppression,
    #   neurotoxicity, death)."
    # -----------------------------------------------------------------------
    "DPYD": {
        "display_name": "DPYD",
        "description": "Dihydropyrimidine dehydrogenase — fluoropyrimidine metabolism",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2A": "Intermediate Metabolizer",
            "*2A/*2A": "Poor Metabolizer",
        },
        "drug_recommendations": {
            "5-fluorouracil": {
                "Normal Metabolizer": {
                    "recommendation": "Use 5-fluorouracil per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce 5-fluorouracil dose by 50%. DPYD *1/*2A carriers have partial DPD deficiency with increased risk of severe toxicity (mucositis, diarrhea, myelosuppression). Monitor closely",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "CONTRAINDICATED: Do NOT administer 5-fluorouracil. DPYD *2A/*2A — complete DPD deficiency. Life-threatening toxicity (severe mucositis, myelosuppression, neurotoxicity, death). Use alternative non-fluoropyrimidine regimen",
                    "action": "CONTRAINDICATED",
                    "alert_level": "CRITICAL",
                },
            },
            "capecitabine": {
                "Normal Metabolizer": {
                    "recommendation": "Use capecitabine per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce capecitabine dose by 50%. DPYD *1/*2A carriers have partial DPD deficiency with increased risk of severe toxicity (hand-foot syndrome, diarrhea, myelosuppression). Monitor closely",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "CONTRAINDICATED: Do NOT administer capecitabine. DPYD *2A/*2A — complete DPD deficiency. Life-threatening toxicity (capecitabine is converted to 5-FU in vivo). Use alternative non-fluoropyrimidine regimen",
                    "action": "CONTRAINDICATED",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    # -----------------------------------------------------------------------
    # G6PD — Important for rasburicase, dapsone, primaquine
    # X-linked gene: males are hemizygous (normal or deficient),
    # females can be heterozygous (intermediate).
    # Uses genotype_field approach for simplicity.
    # CPIC Guideline: Relling et al., Clin Pharmacol Ther, 2014
    # NOTE for agent.py: Add critical alert for G6PD Deficient —
    #   rasburicase is CONTRAINDICATED (risk of severe hemolytic anemia).
    #   Alert text: "G6PD Deficient: Rasburicase is CONTRAINDICATED —
    #   risk of severe hemolytic anemia and methemoglobinemia.
    #   Dapsone and primaquine should also be avoided."
    # -----------------------------------------------------------------------
    "G6PD": {
        "display_name": "G6PD",
        "description": "Glucose-6-phosphate dehydrogenase — X-linked enzyme protecting against oxidative damage",
        "genotype_field": "G6PD",
        "genotype_to_phenotype": {
            "normal": "Normal",
            "Normal": "Normal",
            "deficient": "Deficient",
            "Deficient": "Deficient",
            "intermediate": "Intermediate",
            "Intermediate": "Intermediate",
        },
        "drug_recommendations": {
            "rasburicase": {
                "Normal": {
                    "recommendation": "Use rasburicase per standard dosing for tumor lysis syndrome prophylaxis/treatment",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate": {
                    "recommendation": "Use rasburicase with caution in G6PD intermediate patients (heterozygous females). Monitor for hemolysis (LDH, haptoglobin, reticulocyte count, peripheral smear). Consider alternative urate-lowering strategy",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Deficient": {
                    "recommendation": "CONTRAINDICATED: Do NOT administer rasburicase. G6PD deficiency — risk of severe hemolytic anemia and methemoglobinemia. Use alternative urate-lowering therapy (allopurinol, febuxostat with hydration)",
                    "action": "CONTRAINDICATED",
                    "alert_level": "CRITICAL",
                },
            },
            "dapsone": {
                "Normal": {
                    "recommendation": "Use dapsone per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate": {
                    "recommendation": "Use dapsone with caution. Monitor for hemolysis (CBC, reticulocyte count, methemoglobin). Consider alternative prophylaxis (atovaquone, pentamidine)",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
                "Deficient": {
                    "recommendation": "AVOID dapsone — G6PD deficiency increases risk of dose-dependent hemolytic anemia. Use alternative agents (atovaquone for PCP prophylaxis, alternative for dermatitis herpetiformis)",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "primaquine": {
                "Normal": {
                    "recommendation": "Use primaquine per standard dosing for P. vivax/ovale radical cure",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate": {
                    "recommendation": "Use primaquine with caution. Consider extended 8-week regimen (0.75 mg/kg weekly) with close hemolysis monitoring. Alternative: tafenoquine is also contraindicated in G6PD deficiency",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
                "Deficient": {
                    "recommendation": "AVOID primaquine — G6PD deficiency causes severe, potentially fatal hemolytic anemia with primaquine. Tafenoquine is also contraindicated. Consider chloroquine prophylaxis without radical cure, or consult infectious disease specialist",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    # -----------------------------------------------------------------------
    # HLA-B*58:01 — CPIC Level 1A for allopurinol
    # Associated with allopurinol hypersensitivity syndrome (AHS):
    # Stevens-Johnson syndrome (SJS) / toxic epidermal necrolysis (TEN).
    # Prevalence: ~6-8% in Southeast Asian, ~3.8% in African American,
    # ~1-2% in European populations.
    # CPIC Guideline: Hershfield et al., Clin Pharmacol Ther, 2013
    # NOTE for agent.py: Add critical alert for HLA-B*58:01 Positive —
    #   allopurinol is CONTRAINDICATED (SJS/TEN risk).
    #   Alert text: "HLA-B*58:01 Positive: Allopurinol is CONTRAINDICATED —
    #   risk of Stevens-Johnson syndrome / toxic epidermal necrolysis.
    #   Higher prevalence in Southeast Asian and African American populations.
    #   Use febuxostat or probenecid as alternatives."
    # -----------------------------------------------------------------------
    "HLA-B*58:01": {
        "display_name": "HLA-B*58:01",
        "description": "Human leukocyte antigen B*58:01 — allopurinol hypersensitivity (SJS/TEN)",
        "genotype_field": "HLA_B_5801",
        "genotype_to_phenotype": {
            "negative": "Negative",
            "Negative": "Negative",
            "positive": "Positive",
            "Positive": "Positive",
            "carrier": "Positive",
            "Carrier": "Positive",
        },
        "drug_recommendations": {
            "allopurinol": {
                "Negative": {
                    "recommendation": "Allopurinol use is permitted. HLA-B*58:01 negative — low risk of allopurinol hypersensitivity syndrome",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Positive": {
                    "recommendation": "CONTRAINDICATED: Do NOT prescribe allopurinol. HLA-B*58:01 positive — high risk of Stevens-Johnson syndrome (SJS) / toxic epidermal necrolysis (TEN). Especially high prevalence in Southeast Asian and African American populations. Use febuxostat or probenecid as urate-lowering alternatives",
                    "action": "CONTRAINDICATED",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
    "CYP3A5": {
        "display_name": "CYP3A5",
        "description": "Cytochrome P450 3A5 — tacrolimus and cyclosporine metabolism",
        "allele_to_phenotype": {
            "*1/*1": "Extensive Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*3/*3": "Poor Metabolizer",
            "*1/*6": "Intermediate Metabolizer",
            "*3/*6": "Poor Metabolizer",
            "*6/*6": "Poor Metabolizer",
        },
        "drug_recommendations": {
            "tacrolimus": {
                "Extensive Metabolizer": {
                    "recommendation": "Increase tacrolimus starting dose to 1.5-2x standard (0.3 mg/kg/day). CYP3A5 expressers clear tacrolimus rapidly. Therapeutic drug monitoring essential",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "WARNING",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use standard tacrolimus starting dose (0.2 mg/kg/day). Therapeutic drug monitoring recommended. May need dose adjustment based on trough levels",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Use standard or reduced tacrolimus starting dose (0.15-0.2 mg/kg/day). CYP3A5 non-expressers have slower clearance. Therapeutic drug monitoring essential to avoid supratherapeutic levels and nephrotoxicity",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
            },
            "cyclosporine": {
                "Extensive Metabolizer": {
                    "recommendation": "May require higher cyclosporine doses. CYP3A5 expressers metabolize cyclosporine faster. Monitor C2 levels closely",
                    "action": "DOSE_ADJUSTMENT",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Use cyclosporine per standard dosing with therapeutic drug monitoring",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "Standard or reduced cyclosporine dosing. Monitor C2 levels to avoid nephrotoxicity",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
            },
        },
    },
    "UGT1A1": {
        "display_name": "UGT1A1",
        "description": "UDP-glucuronosyltransferase 1A1 — irinotecan and atazanavir metabolism",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*28": "Intermediate Metabolizer",
            "*1/*6": "Intermediate Metabolizer",
            "*28/*28": "Poor Metabolizer",
            "*6/*6": "Poor Metabolizer",
            "*6/*28": "Poor Metabolizer",
        },
        "drug_recommendations": {
            "irinotecan": {
                "Normal Metabolizer": {
                    "recommendation": "Use irinotecan per standard dosing",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce irinotecan dose by 25-30% if using doses >=250 mg/m2. UGT1A1*28 carriers have reduced SN-38 glucuronidation. Monitor for neutropenia and diarrhea",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce irinotecan dose by 30-50%. UGT1A1*28/*28 (Gilbert's-associated) — significantly impaired SN-38 metabolism with high risk of severe neutropenia and diarrhea. FDA label includes UGT1A1*28 dosing guidance",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "CRITICAL",
                },
            },
            "atazanavir": {
                "Normal Metabolizer": {
                    "recommendation": "Use atazanavir per standard dosing (300 mg + ritonavir 100 mg daily)",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Monitor for jaundice/hyperbilirubinemia. UGT1A1*28 carriers have higher bilirubin levels on atazanavir. Usually cosmetic but may affect adherence",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Poor Metabolizer": {
                    "recommendation": "High risk of clinically significant hyperbilirubinemia on atazanavir. Consider alternative protease inhibitor (darunavir) if jaundice is concerning to patient. UGT1A1*28/*28 genotype",
                    "action": "CONSIDER_ALTERNATIVE",
                    "alert_level": "WARNING",
                },
            },
        },
    },
    "NUDT15": {
        "display_name": "NUDT15",
        "description": "Nudix hydrolase 15 — thiopurine metabolite dephosphorylation",
        "allele_to_phenotype": {
            "*1/*1": "Normal Metabolizer",
            "*1/*2": "Intermediate Metabolizer",
            "*1/*3": "Intermediate Metabolizer",
            "*2/*2": "Poor Metabolizer",
            "*2/*3": "Poor Metabolizer",
            "*3/*3": "Poor Metabolizer",
        },
        "drug_recommendations": {
            "azathioprine": {
                "Normal Metabolizer": {
                    "recommendation": "Use azathioprine per standard dosing. NUDT15 normal — standard thiopurine toxicity risk",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce azathioprine dose by 25-50%. NUDT15 intermediate metabolizers have increased thioguanine nucleotide levels and myelosuppression risk. Monitor CBC weekly for first 8 weeks",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce azathioprine to 10% of standard dose or AVOID. NUDT15 poor metabolizers have severe myelosuppression risk. This is especially important in East Asian populations where NUDT15*3 frequency is 7-10%",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
            "6-mercaptopurine": {
                "Normal Metabolizer": {
                    "recommendation": "Use 6-mercaptopurine per standard dosing with NUDT15 normal status",
                    "action": "STANDARD_DOSING",
                    "alert_level": "INFO",
                },
                "Intermediate Metabolizer": {
                    "recommendation": "Reduce 6-mercaptopurine dose by 25-50%. Monitor CBC closely for myelosuppression",
                    "action": "DOSE_REDUCTION",
                    "alert_level": "WARNING",
                },
                "Poor Metabolizer": {
                    "recommendation": "Reduce 6-mercaptopurine to 10% of standard dose. Severe myelosuppression risk — especially in East Asian populations",
                    "action": "AVOID",
                    "alert_level": "CRITICAL",
                },
            },
        },
    },
}


# ---------------------------------------------------------------------------
# Known drug-drug interactions relevant to PGx recommendations
# Sources: FDA drug labels, CPIC guidelines, Lexicomp
# ---------------------------------------------------------------------------
DRUG_INTERACTIONS = [
    {
        "drug_a": "codeine",
        "drug_b": "tramadol",
        "severity": "high",
        "mechanism": "Both metabolized by CYP2D6; combined use increases CNS depression risk",
        "recommendation": "Avoid concurrent use",
    },
    {
        "drug_a": "clopidogrel",
        "drug_b": "omeprazole",
        "severity": "high",
        "mechanism": "Omeprazole inhibits CYP2C19, reducing clopidogrel activation",
        "recommendation": "Use pantoprazole instead of omeprazole with clopidogrel",
    },
    {
        "drug_a": "simvastatin",
        "drug_b": "amiodarone",
        "severity": "high",
        "mechanism": "Amiodarone inhibits SLCO1B1/CYP3A4, increasing simvastatin myopathy risk",
        "recommendation": "Limit simvastatin to 20mg/day or switch to rosuvastatin",
    },
    {
        "drug_a": "warfarin",
        "drug_b": "amiodarone",
        "severity": "high",
        "mechanism": "Amiodarone inhibits CYP2C9, increasing warfarin bleeding risk",
        "recommendation": "Reduce warfarin dose by 30-50% and monitor INR closely",
    },
    {
        "drug_a": "methotrexate",
        "drug_b": "trimethoprim",
        "severity": "high",
        "mechanism": "Both are folate antagonists; combined use increases bone marrow suppression",
        "recommendation": "Avoid concurrent use or monitor CBC closely",
    },
    {
        "drug_a": "carbamazepine",
        "drug_b": "simvastatin",
        "severity": "moderate",
        "mechanism": "Carbamazepine induces CYP3A4, reducing simvastatin efficacy",
        "recommendation": "Consider higher statin dose or switch to rosuvastatin",
    },
    {
        "drug_a": "fluoxetine",
        "drug_b": "codeine",
        "severity": "high",
        "mechanism": "Fluoxetine strongly inhibits CYP2D6, blocking codeine-to-morphine conversion",
        "recommendation": "Codeine will be ineffective; use non-CYP2D6 analgesic",
    },
    {
        "drug_a": "fluoxetine",
        "drug_b": "tramadol",
        "severity": "high",
        "mechanism": "Fluoxetine inhibits CYP2D6 and combined serotonergic effect increases seizure/serotonin syndrome risk",
        "recommendation": "Avoid concurrent use",
    },
    {
        "drug_a": "tacrolimus",
        "drug_b": "itraconazole",
        "severity": "high",
        "mechanism": "Itraconazole potently inhibits CYP3A4/5, dramatically increasing tacrolimus levels",
        "recommendation": "If unavoidable, reduce tacrolimus dose by 60-80% with intensive TDM",
    },
    {
        "drug_a": "irinotecan",
        "drug_b": "atazanavir",
        "severity": "high",
        "mechanism": "Both metabolized by UGT1A1; atazanavir inhibits UGT1A1 increasing irinotecan toxicity",
        "recommendation": "Avoid concurrent use; if required, reduce irinotecan dose significantly",
    },
    {
        "drug_a": "azathioprine",
        "drug_b": "allopurinol",
        "severity": "high",
        "mechanism": "Allopurinol inhibits xanthine oxidase, shifting thiopurine metabolism toward cytotoxic 6-TGN",
        "recommendation": "Reduce azathioprine dose by 67-75% if allopurinol required, with intensive CBC monitoring",
    },
    {
        "drug_a": "warfarin",
        "drug_b": "ibuprofen",
        "severity": "high",
        "mechanism": "Both affected by CYP2C9; ibuprofen inhibits COX-1 platelets and displaces warfarin from albumin",
        "recommendation": "Avoid concurrent use if possible; use acetaminophen for analgesia",
    },
    {
        "drug_a": "tacrolimus",
        "drug_b": "voriconazole",
        "severity": "high",
        "mechanism": "Voriconazole inhibits CYP3A4/5, increasing tacrolimus levels 2-3x",
        "recommendation": "Reduce tacrolimus dose by 60% and monitor trough levels daily",
    },
]


class PharmacogenomicMapper:
    """Maps star alleles and genotypes to drug recommendations.

    Follows CPIC Level 1A guidelines for 14 pharmacogenes covering
    major drug classes including analgesics, antiplatelets, statins,
    anticoagulants, antimetabolites, thiopurines, antiretrovirals,
    fluoropyrimidines, oxidative hemolysis triggers, and xanthine
    oxidase inhibitors.
    """

    def __init__(self) -> None:
        self._gene_configs = PGX_GENE_CONFIGS

    def _resolve_phenotype(
        self,
        gene: str,
        star_alleles: Optional[str] = None,
        genotype: Optional[str] = None,
    ) -> Optional[str]:
        """Resolve phenotype from star alleles or genotype.

        Args:
            gene: Gene name (e.g., 'CYP2D6').
            star_alleles: Star allele notation (e.g., '*1/*4').
            genotype: Genotype notation (e.g., 'TC').

        Returns:
            Phenotype string or None if not resolved.
        """
        config = self._gene_configs.get(gene)
        if not config:
            return None

        # Try star allele lookup first
        if star_alleles and "allele_to_phenotype" in config:
            phenotype = config["allele_to_phenotype"].get(star_alleles)
            if phenotype:
                # Check for conflict if genotype is also provided
                if genotype and "genotype_to_phenotype" in config:
                    gt_phenotype = config["genotype_to_phenotype"].get(genotype)
                    if gt_phenotype and gt_phenotype != phenotype:
                        logger.warning(
                            f"PGx conflict for {gene}: star alleles ({star_alleles}) -> "
                            f"{phenotype}, but genotype ({genotype}) -> {gt_phenotype}. "
                            f"Using star allele result."
                        )
                return phenotype

        # Try genotype lookup
        if genotype and "genotype_to_phenotype" in config:
            phenotype = config["genotype_to_phenotype"].get(genotype)
            if phenotype:
                return phenotype

        return None

    def map_gene(
        self,
        gene: str,
        star_alleles: Optional[str] = None,
        genotype: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Map a single gene to drug recommendations.

        Args:
            gene: Gene name (e.g., 'CYP2D6', 'SLCO1B1').
            star_alleles: Star allele notation (e.g., '*1/*4').
            genotype: Genotype notation (e.g., 'TC' for SLCO1B1).

        Returns:
            Dict with gene, phenotype, affected drugs, and recommendations.
        """
        config = self._gene_configs.get(gene)
        if not config:
            logger.warning(f"Unknown PGx gene: {gene}")
            return {
                "gene": gene,
                "phenotype": None,
                "error": f"Gene {gene} not in PGx database",
                "affected_drugs": [],
                "critical_alerts": [],
            }

        phenotype = self._resolve_phenotype(gene, star_alleles, genotype)
        if not phenotype:
            allele_info = star_alleles or genotype or "not provided"
            logger.warning(f"Could not resolve phenotype for {gene} ({allele_info})")
            return {
                "gene": gene,
                "display_name": config["display_name"],
                "description": config["description"],
                "star_alleles": star_alleles,
                "genotype": genotype,
                "phenotype": None,
                "error": f"Unrecognized allele/genotype: {allele_info}",
                "affected_drugs": [],
                "critical_alerts": [],
            }

        affected_drugs = []
        critical_alerts = []

        for drug, phenotype_recs in config.get("drug_recommendations", {}).items():
            rec = phenotype_recs.get(phenotype)
            if rec:
                drug_result = {
                    "drug": drug,
                    "recommendation": rec["recommendation"],
                    "action": rec["action"],
                    "alert_level": rec["alert_level"],
                }
                affected_drugs.append(drug_result)

                if rec["alert_level"] == "CRITICAL":
                    critical_alerts.append({
                        "drug": drug,
                        "gene": gene,
                        "phenotype": phenotype,
                        "action": rec["action"],
                        "message": rec["recommendation"],
                    })

        return {
            "gene": gene,
            "display_name": config["display_name"],
            "description": config["description"],
            "star_alleles": star_alleles,
            "genotype": genotype,
            "phenotype": phenotype,
            "affected_drugs": affected_drugs,
            "critical_alerts": critical_alerts,
        }

    def get_guideline_versions(self) -> Dict[str, Dict[str, str]]:
        """Return CPIC guideline version info for all supported genes.

        Useful for audit trails and ensuring recommendations are current.
        """
        return dict(CPIC_GUIDELINE_VERSIONS)

    def map_all(
        self,
        star_alleles: Optional[Dict[str, str]] = None,
        genotypes: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """Map all provided PGx genes to drug recommendations.

        Args:
            star_alleles: Dict of gene -> star alleles (e.g.,
                {'CYP2D6': '*1/*4', 'TPMT': '*1/*3A'}).
            genotypes: Dict of rsID or gene -> genotype (e.g.,
                {'SLCO1B1_rs4149056': 'TC', 'VKORC1_rs9923231': 'AG',
                 'MTHFR_rs1801133': 'CT', 'HLA_B5701': 'positive'}).

        Returns:
            Dict with gene_results list, all critical_alerts, and summary.
        """
        star_alleles = star_alleles or {}
        genotypes = genotypes or {}

        gene_results = []
        all_critical_alerts = []
        drugs_to_avoid = []
        drugs_to_adjust = []

        # Process star allele genes
        for gene, alleles in star_alleles.items():
            result = self.map_gene(gene, star_alleles=alleles)
            gene_results.append(result)
            all_critical_alerts.extend(result.get("critical_alerts", []))

        # Process genotype-based genes
        # Map genotype field names to gene names
        genotype_field_to_gene = {}
        for gene, config in self._gene_configs.items():
            field = config.get("genotype_field")
            if field:
                genotype_field_to_gene[field] = gene

        for field, gt in genotypes.items():
            gene = genotype_field_to_gene.get(field)
            if gene:
                # Skip if already processed via star alleles
                if gene not in star_alleles:
                    result = self.map_gene(gene, genotype=gt)
                    gene_results.append(result)
                    all_critical_alerts.extend(result.get("critical_alerts", []))

        # Compile drug action lists
        for result in gene_results:
            for drug_rec in result.get("affected_drugs", []):
                if drug_rec["action"] in ("AVOID", "CONTRAINDICATED"):
                    drugs_to_avoid.append({
                        "drug": drug_rec["drug"],
                        "gene": result["gene"],
                        "phenotype": result["phenotype"],
                        "reason": drug_rec["recommendation"],
                    })
                elif drug_rec["action"] in ("DOSE_REDUCTION", "DOSE_ADJUSTMENT", "CONSIDER_ALTERNATIVE"):
                    drugs_to_adjust.append({
                        "drug": drug_rec["drug"],
                        "gene": result["gene"],
                        "phenotype": result["phenotype"],
                        "action": drug_rec["action"],
                        "reason": drug_rec["recommendation"],
                    })

        # Detect drug-drug interaction warnings (same CYP pathway)
        drug_interaction_warnings = []
        cyp_drugs = {}  # gene -> list of drugs
        for result in gene_results:
            gene = result.get("gene", "")
            for drug_rec in result.get("affected_drugs", []):
                if drug_rec["action"] != "STANDARD_DOSING":
                    cyp_drugs.setdefault(gene, []).append(drug_rec["drug"])
        for gene, drugs in cyp_drugs.items():
            if len(drugs) > 1:
                drug_interaction_warnings.append({
                    "gene": gene,
                    "drugs": drugs,
                    "warning": (
                        f"Multiple drugs affected by {gene} metabolizer status: "
                        f"{', '.join(drugs)}. Concomitant use may compound pharmacokinetic effects."
                    ),
                })

        # Build preliminary result for drug interaction checking
        preliminary_result = {
            "gene_results": gene_results,
            "critical_alerts": all_critical_alerts,
            "drugs_to_avoid": drugs_to_avoid,
            "drugs_to_adjust": drugs_to_adjust,
            "drug_interaction_warnings": drug_interaction_warnings,
            "genes_analyzed": len(gene_results),
            "has_critical_findings": len(all_critical_alerts) > 0,
        }

        # Cross-check for drug-drug interactions across all recommendations
        drug_interactions = self.check_drug_interactions(preliminary_result)
        preliminary_result["drug_interactions"] = drug_interactions

        # Add guideline version info for audit trail
        genes_analyzed = [r["gene"] for r in gene_results if r.get("gene")]
        preliminary_result["guideline_versions"] = {
            gene: CPIC_GUIDELINE_VERSIONS.get(gene, {}).get("update", "unknown")
            for gene in genes_analyzed
        }

        logger.info(
            f"PGx mapping complete: {len(gene_results)} genes analyzed, "
            f"{len(all_critical_alerts)} critical alerts, "
            f"{len(drugs_to_avoid)} drugs to avoid, "
            f"{len(drug_interactions)} drug-drug interactions detected"
        )

        return preliminary_result

    def check_drug_interactions(
        self,
        map_all_result: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        """Check for known drug-drug interactions across PGx recommendations.

        Scans all drugs mentioned in gene_results and cross-checks them
        against the DRUG_INTERACTIONS table. Returns warnings for any
        dangerous combinations found in the patient's recommended drugs.

        Args:
            map_all_result: The dict returned by map_all(), containing
                a "gene_results" list with per-gene drug recommendations.

        Returns:
            List of interaction warning dicts, each containing:
                - drug_a, drug_b: the interacting drug names
                - severity: "high" or "moderate"
                - mechanism: pharmacological explanation
                - recommendation: clinical guidance
                - genes_involved: list of genes that recommended these drugs
        """
        # Collect all recommended drugs and which gene(s) recommended them
        drug_to_genes: Dict[str, List[str]] = {}
        for gene_result in map_all_result.get("gene_results", []):
            gene = gene_result.get("gene", "")
            for drug_rec in gene_result.get("affected_drugs", []):
                drug_name = drug_rec["drug"].lower()
                drug_to_genes.setdefault(drug_name, [])
                if gene not in drug_to_genes[drug_name]:
                    drug_to_genes[drug_name].append(gene)

        # Check each known interaction pair
        warnings: List[Dict[str, Any]] = []
        recommended_drugs = set(drug_to_genes.keys())

        for interaction in DRUG_INTERACTIONS:
            a = interaction["drug_a"].lower()
            b = interaction["drug_b"].lower()
            if a in recommended_drugs and b in recommended_drugs:
                genes_involved = sorted(
                    set(drug_to_genes[a] + drug_to_genes[b])
                )
                warnings.append({
                    "drug_a": interaction["drug_a"],
                    "drug_b": interaction["drug_b"],
                    "severity": interaction["severity"],
                    "mechanism": interaction["mechanism"],
                    "recommendation": interaction["recommendation"],
                    "genes_involved": genes_involved,
                })

        if warnings:
            logger.warning(
                f"Drug-drug interactions detected: {len(warnings)} "
                f"interaction(s) across PGx recommendations"
            )

        return warnings

    def check_drug(
        self,
        drug_name: str,
        star_alleles: Optional[Dict[str, str]] = None,
        genotypes: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        """Check a specific drug against all available PGx data.

        Args:
            drug_name: Drug name to check (case-insensitive).
            star_alleles: Dict of gene -> star alleles.
            genotypes: Dict of rsID/gene -> genotype.

        Returns:
            Dict with drug name, safety assessment, and gene-specific findings.
        """
        drug_lower = drug_name.lower()
        star_alleles = star_alleles or {}
        genotypes = genotypes or {}

        findings = []
        is_safe = True
        overall_action = "STANDARD_DOSING"

        # Check each gene config for this drug
        for gene, config in self._gene_configs.items():
            if drug_lower not in config.get("drug_recommendations", {}):
                continue

            # Resolve phenotype
            alleles = star_alleles.get(gene)
            gt = None
            field = config.get("genotype_field")
            if field:
                gt = genotypes.get(field)

            if not alleles and not gt:
                findings.append({
                    "gene": gene,
                    "status": "NOT_TESTED",
                    "message": f"{gene} genotype not available — cannot assess {drug_name} safety",
                })
                continue

            phenotype = self._resolve_phenotype(gene, alleles, gt)
            if not phenotype:
                continue

            rec = config["drug_recommendations"][drug_lower].get(phenotype)
            if rec:
                findings.append({
                    "gene": gene,
                    "phenotype": phenotype,
                    "alleles": alleles or gt,
                    "action": rec["action"],
                    "alert_level": rec["alert_level"],
                    "recommendation": rec["recommendation"],
                })

                if rec["action"] in ("AVOID", "CONTRAINDICATED"):
                    is_safe = False
                    overall_action = rec["action"]
                elif rec["action"] in ("DOSE_REDUCTION", "DOSE_ADJUSTMENT") and overall_action == "STANDARD_DOSING":
                    overall_action = rec["action"]

        # Flag if the drug isn't in any PGx gene's recommendations
        drug_found = any(
            drug_lower in config.get("drug_recommendations", {})
            for config in self._gene_configs.values()
        )

        if not drug_found:
            findings.append({
                "gene": "N/A",
                "status": "NOT_IN_DATABASE",
                "message": f"{drug_name} is not in the pharmacogenomic database. "
                           "No PGx guidance available.",
            })

        return {
            "drug": drug_name,
            "drug_in_database": drug_found,
            "is_safe": is_safe,
            "overall_action": overall_action,
            "findings": findings,
            "genes_checked": len(findings),
        }
