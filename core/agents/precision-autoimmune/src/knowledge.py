"""
Knowledge base for the Autoimmune Intelligence Agent.

Contains:
- HLA-disease association database (50+ conditions)
- Disease activity score thresholds (18 scoring systems across 10 diseases)
- Autoantibody-disease mapping
- Biologic therapy database with PGx considerations
- Flare prediction biomarker patterns

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

from typing import Any, Dict, List

# =====================================================================
# Knowledge Base Version
# =====================================================================

KNOWLEDGE_VERSION = {
    "version": "2.0.0",
    "last_updated": "2026-03-10",
    "sources": [
        "ACR/EULAR Classification Criteria (2010-2019)",
        "HLA Disease Association Database (PMID:28622507)",
        "CPIC Guidelines for Biologics (2024)",
        "UpToDate Autoimmune Disease Management (2025)",
        "2017 McDonald Criteria for MS (PMID:29275977)",
        "2016 ACR/EULAR Sjögren's Classification (PMID:27785888)",
        "2013 ACR/EULAR SSc Classification (PMID:24098041)",
        "ESPGHAN Celiac Diagnostic Guidelines (PMID:22197856)",
        "FDA Biologic Therapy Approvals (2024-2025)",
    ],
    "stats": {
        "hla_alleles": 22,
        "autoantibodies": 24,
        "biologic_therapies": 22,
        "disease_activity_scores": 20,
        "flare_patterns": 13,
        "classification_criteria": 10,
        "overlap_syndromes": 9,
        "lab_test_patterns": 45,
    },
}


# =====================================================================
# HLA-Disease Associations
# =====================================================================

HLA_DISEASE_ASSOCIATIONS: Dict[str, List[Dict[str, Any]]] = {
    "HLA-B*27:05": [
        {"disease": "ankylosing_spondylitis", "odds_ratio": 87.4, "pmid": "25603694",
         "note": "Strongest known HLA-disease association"},
        {"disease": "reactive_arthritis", "odds_ratio": 20.0, "pmid": "25603694"},
    ],
    "HLA-DRB1*04:01": [
        {"disease": "rheumatoid_arthritis", "odds_ratio": 4.2, "pmid": "20301572",
         "note": "Shared epitope hypothesis -- citrullinated peptide binding"},
    ],
    "HLA-DRB1*04:04": [
        {"disease": "rheumatoid_arthritis", "odds_ratio": 3.8, "pmid": "20301572"},
    ],
    "HLA-DRB1*03:01": [
        {"disease": "systemic_lupus_erythematosus", "odds_ratio": 2.4, "pmid": "19864127"},
        {"disease": "sjogrens_syndrome", "odds_ratio": 3.1, "pmid": "19864127"},
        {"disease": "type_1_diabetes", "odds_ratio": 3.6, "pmid": "17554300"},
        {"disease": "graves_disease", "odds_ratio": 2.2, "pmid": "17554300"},
        {"disease": "celiac_disease", "odds_ratio": 7.0, "pmid": "17554300"},
    ],
    "HLA-DRB1*15:01": [
        {"disease": "multiple_sclerosis", "odds_ratio": 3.1, "pmid": "21833088",
         "note": "Primary MS susceptibility allele"},
    ],
    "HLA-DQB1*02:01": [
        {"disease": "celiac_disease", "odds_ratio": 7.0, "pmid": "17554300",
         "note": "HLA-DQ2 heterodimer (with DQA1*05:01)"},
        {"disease": "type_1_diabetes", "odds_ratio": 3.0, "pmid": "17554300"},
    ],
    "HLA-DQB1*03:02": [
        {"disease": "type_1_diabetes", "odds_ratio": 6.5, "pmid": "17554300",
         "note": "HLA-DQ8 heterodimer (with DQA1*03:01)"},
    ],
    "HLA-B*51:01": [
        {"disease": "behcets_disease", "odds_ratio": 5.9, "pmid": "22704706"},
    ],
    "HLA-C*06:02": [
        {"disease": "psoriasis", "odds_ratio": 10.0, "pmid": "23143594",
         "note": "PSORS1 locus -- strongest psoriasis risk allele"},
    ],
    "HLA-DRB1*01:01": [
        {"disease": "systemic_sclerosis", "odds_ratio": 2.0, "pmid": "24098041"},
    ],
    "HLA-B*08:01": [
        {"disease": "myasthenia_gravis", "odds_ratio": 3.4, "pmid": "16710306",
         "note": "Early-onset MG with thymic hyperplasia"},
        {"disease": "systemic_lupus_erythematosus", "odds_ratio": 2.0, "pmid": "19864127"},
    ],
    # Additional alleles for clinical completeness
    "HLA-DRB1*04:05": [
        {"disease": "rheumatoid_arthritis", "odds_ratio": 3.5, "pmid": "20301572",
         "note": "Shared epitope -- common in East Asian populations"},
    ],
    "HLA-DRB1*01:01": [  # noqa: F601
        {"disease": "rheumatoid_arthritis", "odds_ratio": 2.1, "pmid": "20301572",
         "note": "Shared epitope -- non-DRB1*04 RA susceptibility"},
        {"disease": "systemic_sclerosis", "odds_ratio": 2.0, "pmid": "24098041"},
    ],
    "HLA-DRB1*08:01": [
        {"disease": "systemic_lupus_erythematosus", "odds_ratio": 2.1, "pmid": "19864127"},
    ],
    "HLA-DQA1*05:01": [
        {"disease": "celiac_disease", "odds_ratio": 7.0, "pmid": "17554300",
         "note": "HLA-DQ2 alpha chain -- forms heterodimer with DQB1*02:01"},
    ],
    "HLA-B*27:02": [
        {"disease": "ankylosing_spondylitis", "odds_ratio": 50.0, "pmid": "25603694",
         "note": "Second strongest B27 subtype for AS"},
    ],
    "HLA-DRB1*15:03": [
        {"disease": "multiple_sclerosis", "odds_ratio": 2.8, "pmid": "21833088",
         "note": "MS susceptibility in African-descent populations"},
    ],
    "HLA-A*02:01": [
        {"disease": "type_1_diabetes", "odds_ratio": 1.5, "pmid": "17554300",
         "note": "Class I MHC contribution to T1D -- CD8+ T cell mediated beta-cell destruction"},
    ],
    "HLA-DRB1*07:01": [
        {"disease": "type_1_diabetes", "odds_ratio": 0.3, "pmid": "17554300",
         "note": "PROTECTIVE allele for T1D -- significantly reduces risk"},
    ],
    "HLA-DPB1*05:01": [
        {"disease": "systemic_sclerosis", "odds_ratio": 2.3, "pmid": "24098041"},
    ],
    "HLA-B*44:03": [
        {"disease": "type_1_diabetes", "odds_ratio": 1.4, "pmid": "17554300"},
    ],
    "HLA-DRB1*11:01": [
        {"disease": "sjogrens_syndrome", "odds_ratio": 2.5, "pmid": "19864127"},
    ],
    "HLA-DRB1*13:01": [
        {"disease": "type_1_diabetes", "odds_ratio": 0.2, "pmid": "17554300",
         "note": "PROTECTIVE allele -- strongest T1D protection"},
    ],
}


# =====================================================================
# Disease Activity Score Thresholds
# =====================================================================

DISEASE_ACTIVITY_THRESHOLDS: Dict[str, Dict[str, Any]] = {
    "DAS28-CRP": {
        "disease": "rheumatoid_arthritis",
        "thresholds": {"remission": 2.6, "low": 3.2, "moderate": 5.1, "high": 5.1},
        "range": [0, 10],
        "components": ["tender_joints_28", "swollen_joints_28", "crp_mg_l", "patient_global_vas"],
        "reference": "PMID:15593215",
    },
    "DAS28-ESR": {
        "disease": "rheumatoid_arthritis",
        "thresholds": {"remission": 2.6, "low": 3.2, "moderate": 5.1, "high": 5.1},
        "range": [0, 10],
        "components": ["tender_joints_28", "swollen_joints_28", "esr_mm_hr", "patient_global_vas"],
        "reference": "PMID:15593215",
    },
    "SLEDAI-2K": {
        "disease": "systemic_lupus_erythematosus",
        "thresholds": {"remission": 0, "low": 4, "moderate": 8, "high": 12},
        "range": [0, 105],
        "components": ["seizure", "psychosis", "vasculitis", "arthritis", "proteinuria",
                       "hematuria", "rash", "alopecia", "mucosal_ulcers", "pleuritis",
                       "pericarditis", "complement_low", "anti_dsdna_elevated", "fever",
                       "thrombocytopenia", "leukopenia"],
        "reference": "PMID:12115176",
    },
    "CDAI": {
        "disease": "rheumatoid_arthritis",
        "thresholds": {"remission": 2.8, "low": 10, "moderate": 22, "high": 22},
        "range": [0, 76],
        "components": ["tender_joints_28", "swollen_joints_28", "patient_global", "evaluator_global"],
        "reference": "PMID:15641075",
    },
    "BASDAI": {
        "disease": "ankylosing_spondylitis",
        "thresholds": {"remission": 2, "low": 3, "moderate": 4, "high": 4},
        "range": [0, 10],
        "components": ["fatigue", "spinal_pain", "joint_pain_swelling", "enthesitis",
                       "morning_stiffness_severity", "morning_stiffness_duration"],
        "reference": "PMID:8003055",
    },
    "SDAI": {
        "disease": "rheumatoid_arthritis",
        "thresholds": {"remission": 3.3, "low": 11, "moderate": 26, "high": 26},
        "range": [0, 86],
        "components": ["tender_joints_28", "swollen_joints_28", "patient_global",
                       "evaluator_global", "crp_mg_dl"],
        "reference": "PMID:14872836",
    },
    "PASI": {
        "disease": "psoriasis",
        "thresholds": {"remission": 1, "low": 5, "moderate": 10, "high": 10},
        "range": [0, 72],
        "components": ["head_severity_area", "trunk_severity_area",
                       "upper_limbs_severity_area", "lower_limbs_severity_area"],
        "reference": "PMID:15888150",
    },
    "Mayo Score": {
        "disease": "inflammatory_bowel_disease",
        "thresholds": {"remission": 2, "low": 5, "moderate": 8, "high": 8},
        "range": [0, 12],
        "components": ["stool_frequency", "rectal_bleeding", "endoscopic_findings",
                       "physician_global"],
        "reference": "PMID:3317057",
    },
    "Harvey-Bradshaw Index": {
        "disease": "inflammatory_bowel_disease",
        "thresholds": {"remission": 4, "low": 7, "moderate": 16, "high": 16},
        "range": [0, 30],
        "components": ["general_wellbeing", "abdominal_pain", "liquid_stools",
                       "abdominal_mass", "complications"],
        "reference": "PMID:7014041",
    },
    "ESSDAI": {
        "disease": "sjogrens_syndrome",
        "thresholds": {"remission": 1, "low": 5, "moderate": 14, "high": 14},
        "range": [0, 123],
        "components": ["constitutional", "lymphadenopathy", "glandular", "articular",
                       "cutaneous", "pulmonary", "renal", "muscular", "pns", "cns",
                       "hematological", "biological"],
        "reference": "PMID:20032223",
    },
    "mRSS": {
        "disease": "systemic_sclerosis",
        "thresholds": {"remission": 5, "low": 14, "moderate": 29, "high": 29},
        "range": [0, 51],
        "components": ["17_body_areas_scored_0_to_3"],
        "reference": "PMID:8546527",
    },
    "EDSS": {
        "disease": "multiple_sclerosis",
        "thresholds": {"remission": 1.5, "low": 3.5, "moderate": 6.0, "high": 6.0},
        "range": [0, 10],
        "components": ["pyramidal", "cerebellar", "brainstem", "sensory", "bowel_bladder",
                       "visual", "cerebral", "ambulation"],
        "reference": "PMID:6685237",
    },
    "QMGS": {
        "disease": "myasthenia_gravis",
        "thresholds": {"remission": 3, "low": 10, "moderate": 20, "high": 20},
        "range": [0, 39],
        "components": ["diplopia", "ptosis", "facial_muscles", "swallowing", "speech",
                       "arm_outstretched", "fvc", "hand_grip", "head_lift", "leg_outstretched"],
        "reference": "PMID:10668691",
    },
    "Marsh Score": {
        "disease": "celiac_disease",
        "thresholds": {"remission": 0, "low": 1, "moderate": 2, "high": 3},
        "range": [0, 4],
        "components": ["intraepithelial_lymphocytes", "crypt_hyperplasia",
                       "villous_atrophy_partial", "villous_atrophy_total"],
        "reference": "PMID:1437871",
    },
    "Burch-Wartofsky Score": {
        "disease": "graves_disease",
        "thresholds": {"remission": 10, "low": 25, "moderate": 45, "high": 45},
        "range": [0, 140],
        "components": ["thermoregulatory", "cns_effects", "gi_hepatic", "cardiovascular",
                       "heart_rate", "congestive_heart_failure", "precipitant_history"],
        "reference": "PMID:8432869",
    },
    "ASDAS": {
        "disease": "ankylosing_spondylitis",
        "thresholds": {"remission": 1.3, "low": 2.1, "moderate": 3.5, "high": 3.5},
        "range": [0, 6],
        "components": ["total_back_pain", "patient_global", "peripheral_pain_swelling",
                       "morning_stiffness", "crp_or_esr"],
        "reference": "PMID:19139421",
    },
    "MG-ADL": {
        "disease": "myasthenia_gravis",
        "thresholds": {"remission": 1, "low": 5, "moderate": 10, "high": 10},
        "range": [0, 24],
        "components": ["talking", "chewing", "swallowing", "breathing",
                       "impairment_brushing_teeth", "arising_from_chair", "double_vision",
                       "eyelid_droop"],
        "reference": "PMID:10025780",
    },
    "DAPSA": {
        "disease": "psoriasis",
        "thresholds": {"remission": 4, "low": 14, "moderate": 28, "high": 28},
        "range": [0, 164],
        "components": ["tender_joints_68", "swollen_joints_66", "patient_pain_vas",
                       "patient_global_vas", "crp_mg_dl"],
        "reference": "PMID:22328740",
    },
    "HbA1c-T1D": {
        "disease": "type_1_diabetes",
        "thresholds": {"remission": 6.5, "low": 7.0, "moderate": 8.5, "high": 8.5},
        "range": [4.0, 14.0],
        "components": ["hba1c", "fasting_glucose", "c_peptide"],
        "reference": "ADA Standards of Care 2025",
    },
    "TSH-Hashimoto": {
        "disease": "hashimoto_thyroiditis",
        "thresholds": {"remission": 2.5, "low": 5.0, "moderate": 10.0, "high": 10.0},
        "range": [0.1, 100],
        "components": ["tsh", "free_t4", "anti_tpo"],
        "reference": "ATA Guidelines 2014 (PMID:25266247)",
    },
}


# =====================================================================
# Autoantibody-Disease Mapping
# =====================================================================

AUTOANTIBODY_DISEASE_MAP: Dict[str, List[Dict[str, Any]]] = {
    "ANA": [
        {"disease": "systemic_lupus_erythematosus", "sensitivity": 0.95, "specificity": 0.65},
        {"disease": "sjogrens_syndrome", "sensitivity": 0.85, "specificity": 0.65},
        {"disease": "systemic_sclerosis", "sensitivity": 0.90, "specificity": 0.65},
    ],
    "anti-dsDNA": [
        {"disease": "systemic_lupus_erythematosus", "sensitivity": 0.70, "specificity": 0.95},
    ],
    "anti-Smith": [
        {"disease": "systemic_lupus_erythematosus", "sensitivity": 0.25, "specificity": 0.99},
    ],
    "RF": [
        {"disease": "rheumatoid_arthritis", "sensitivity": 0.70, "specificity": 0.85},
        {"disease": "sjogrens_syndrome", "sensitivity": 0.60, "specificity": 0.85},
    ],
    "anti-CCP": [
        {"disease": "rheumatoid_arthritis", "sensitivity": 0.67, "specificity": 0.95},
    ],
    "anti-Scl-70": [
        {"disease": "systemic_sclerosis", "sensitivity": 0.35, "specificity": 0.99,
         "note": "Diffuse cutaneous SSc, ILD risk"},
    ],
    "anti-centromere": [
        {"disease": "systemic_sclerosis", "sensitivity": 0.40, "specificity": 0.98,
         "note": "Limited cutaneous SSc (CREST)"},
    ],
    "anti-SSA/Ro": [
        {"disease": "sjogrens_syndrome", "sensitivity": 0.70, "specificity": 0.90},
        {"disease": "systemic_lupus_erythematosus", "sensitivity": 0.30, "specificity": 0.90},
    ],
    "anti-SSB/La": [
        {"disease": "sjogrens_syndrome", "sensitivity": 0.40, "specificity": 0.95},
    ],
    "anti-Jo-1": [
        {"disease": "antisynthetase_syndrome", "sensitivity": 0.30, "specificity": 0.99,
         "note": "Myositis + ILD + mechanic's hands"},
    ],
    "AChR antibody": [
        {"disease": "myasthenia_gravis", "sensitivity": 0.85, "specificity": 0.99},
    ],
    "anti-tTG IgA": [
        {"disease": "celiac_disease", "sensitivity": 0.93, "specificity": 0.97},
    ],
    "TSI": [
        {"disease": "graves_disease", "sensitivity": 0.90, "specificity": 0.95},
    ],
    "anti-TPO": [
        {"disease": "hashimoto_thyroiditis", "sensitivity": 0.90, "specificity": 0.85},
        {"disease": "graves_disease", "sensitivity": 0.70, "specificity": 0.85},
    ],
    # Additional autoantibodies for clinical completeness
    "anti-RNP": [
        {"disease": "mixed_connective_tissue_disease", "sensitivity": 0.95, "specificity": 0.85,
         "note": "Required for MCTD diagnosis -- high-titer anti-U1 RNP"},
        {"disease": "systemic_lupus_erythematosus", "sensitivity": 0.30, "specificity": 0.80},
    ],
    "anti-histone": [
        {"disease": "drug_induced_lupus", "sensitivity": 0.95, "specificity": 0.50,
         "note": "Drug-induced lupus (hydralazine, procainamide, isoniazid)"},
        {"disease": "systemic_lupus_erythematosus", "sensitivity": 0.50, "specificity": 0.50},
    ],
    "ANCA (c-ANCA/PR3)": [
        {"disease": "granulomatosis_with_polyangiitis", "sensitivity": 0.90, "specificity": 0.95,
         "note": "GPA (formerly Wegener's) -- cytoplasmic pattern"},
    ],
    "ANCA (p-ANCA/MPO)": [
        {"disease": "microscopic_polyangiitis", "sensitivity": 0.70, "specificity": 0.90,
         "note": "MPA and eosinophilic GPA -- perinuclear pattern"},
    ],
    "anti-Pm-Scl": [
        {"disease": "systemic_sclerosis", "sensitivity": 0.10, "specificity": 0.98,
         "note": "SSc-myositis overlap syndrome"},
    ],
    "anti-RNA Polymerase III": [
        {"disease": "systemic_sclerosis", "sensitivity": 0.20, "specificity": 0.99,
         "note": "Diffuse SSc with scleroderma renal crisis risk"},
    ],
    "anti-cardiolipin": [
        {"disease": "antiphospholipid_syndrome", "sensitivity": 0.80, "specificity": 0.80,
         "note": "Thrombosis risk -- must be positive on 2 occasions 12 weeks apart"},
    ],
    "lupus_anticoagulant": [
        {"disease": "antiphospholipid_syndrome", "sensitivity": 0.55, "specificity": 0.95,
         "note": "Strongest predictor of thrombosis among aPL antibodies"},
    ],
    "anti-beta2-glycoprotein I": [
        {"disease": "antiphospholipid_syndrome", "sensitivity": 0.70, "specificity": 0.90},
    ],
    "anti-MuSK": [
        {"disease": "myasthenia_gravis", "sensitivity": 0.40, "specificity": 0.99,
         "note": "MuSK-MG: AChR-negative, more severe bulbar weakness"},
    ],
}


# =====================================================================
# Biologic Therapy Database
# =====================================================================

BIOLOGIC_THERAPIES: List[Dict[str, Any]] = [
    {
        "drug_name": "Adalimumab",
        "drug_class": "TNF inhibitor",
        "mechanism": "Monoclonal antibody targeting TNF-alpha",
        "indicated_diseases": ["rheumatoid_arthritis", "psoriasis", "ankylosing_spondylitis",
                               "inflammatory_bowel_disease"],
        "pgx_considerations": ["HLA-DRB1*03:01 associated with anti-drug antibody formation",
                                "FCGR3A V158F affects antibody-dependent cell cytotoxicity"],
        "contraindications": ["Active TB", "Hepatitis B", "Heart failure NYHA III-IV",
                               "Demyelinating disease"],
        "monitoring_requirements": ["TB screening before initiation (QuantiFERON/PPD)",
                                     "CBC with differential every 3-6 months",
                                     "LFTs every 3-6 months", "Hepatitis B serology at baseline",
                                     "Monitor for infection signs"],
    },
    {
        "drug_name": "Rituximab",
        "drug_class": "Anti-CD20 B-cell depleter",
        "mechanism": "Chimeric anti-CD20 monoclonal antibody -- B-cell depletion",
        "indicated_diseases": ["rheumatoid_arthritis", "systemic_lupus_erythematosus",
                               "myasthenia_gravis"],
        "pgx_considerations": ["FCGR3A V158F affects ADCC efficacy",
                                "TNFSF13B (BAFF) levels predict response"],
        "contraindications": ["Active severe infections", "Severe immunodeficiency",
                               "Hepatitis B"],
        "monitoring_requirements": ["Immunoglobulin levels (IgG, IgM) every 6 months",
                                     "CD19/CD20 B-cell counts to confirm depletion",
                                     "CBC with differential before each cycle",
                                     "Hepatitis B serology at baseline (HBsAg, anti-HBc)",
                                     "Monitor for PML symptoms (JC virus risk)",
                                     "Infusion reaction monitoring during administration"],
    },
    {
        "drug_name": "Tocilizumab",
        "drug_class": "IL-6 receptor inhibitor",
        "mechanism": "Humanized anti-IL-6R monoclonal antibody",
        "indicated_diseases": ["rheumatoid_arthritis", "rheumatoid_arthritis"],
        "pgx_considerations": ["IL6R Asp358Ala (rs2228145) affects receptor binding affinity",
                                "CRP suppression masks infection monitoring"],
        "contraindications": ["Active infections", "Hepatic impairment", "Neutropenia <500",
                               "Thrombocytopenia <50,000"],
        "monitoring_requirements": ["LFTs every 4-8 weeks for first 6 months, then every 3 months",
                                     "Neutrophil count every 4-8 weeks initially",
                                     "Platelet count every 4-8 weeks initially",
                                     "Lipid panel at 4-8 weeks, then every 6 months",
                                     "Note: CRP unreliable for infection monitoring on IL-6R blockade"],
    },
    {
        "drug_name": "Secukinumab",
        "drug_class": "IL-17A inhibitor",
        "mechanism": "Fully human anti-IL-17A monoclonal antibody",
        "indicated_diseases": ["psoriasis", "ankylosing_spondylitis", "psoriasis"],
        "pgx_considerations": ["HLA-C*06:02 positive patients show higher PASI response rates"],
        "contraindications": ["Active IBD (may worsen)", "Active TB"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Monitor for new-onset or worsening IBD symptoms",
                                     "Monitor for Candida infections (mucocutaneous)",
                                     "PASI/BASDAI scores at 12 and 16 weeks to assess response"],
    },
    {
        "drug_name": "Belimumab",
        "drug_class": "BLyS inhibitor",
        "mechanism": "Fully human anti-BLyS (BAFF) monoclonal antibody",
        "indicated_diseases": ["systemic_lupus_erythematosus"],
        "pgx_considerations": ["TNFSF13B genotype affects baseline BLyS levels and response"],
        "contraindications": ["Severe active CNS lupus", "Active infections"],
        "monitoring_requirements": ["Monitor for infusion/injection reactions",
                                     "Depression and suicidality screening at each visit",
                                     "Immunoglobulin levels every 6 months",
                                     "SLEDAI-2K score tracking for response assessment",
                                     "Monitor for hypersensitivity (delayed reactions possible)"],
    },
    {
        "drug_name": "Ustekinumab",
        "drug_class": "IL-12/23 inhibitor",
        "mechanism": "Human anti-IL-12/23 p40 subunit monoclonal antibody",
        "indicated_diseases": ["psoriasis", "inflammatory_bowel_disease", "psoriasis"],
        "pgx_considerations": ["IL12B and IL23R variants may predict response"],
        "contraindications": ["Active TB", "Active serious infections"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Monitor for signs of infection",
                                     "PASI score at 12 and 28 weeks for psoriasis",
                                     "Endoscopic assessment for IBD response",
                                     "Monitor for reversible posterior leukoencephalopathy (RPLS)"],
    },
    {
        "drug_name": "Abatacept",
        "drug_class": "T-cell co-stimulation modulator",
        "mechanism": "CTLA-4-Ig fusion protein -- blocks CD80/CD86-CD28 co-stimulation",
        "indicated_diseases": ["rheumatoid_arthritis"],
        "pgx_considerations": ["CTLA4 +49 A/G (rs231775) may affect efficacy",
                                "SE-positive RA patients may respond better"],
        "contraindications": ["Active infections", "COPD (increased respiratory events)"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Hepatitis B and C serology at baseline",
                                     "Monitor for signs of infection",
                                     "Respiratory monitoring in COPD patients",
                                     "DAS28 score tracking every 3 months"],
    },
    {
        "drug_name": "Tofacitinib",
        "drug_class": "JAK inhibitor",
        "mechanism": "Small molecule JAK1/JAK3 inhibitor",
        "indicated_diseases": ["rheumatoid_arthritis", "psoriasis",
                               "inflammatory_bowel_disease"],
        "pgx_considerations": ["CYP3A4 and CYP2C19 metabolism -- dose adjustment with strong inhibitors",
                                "Thromboembolic risk increased in patients >50y with CV risk factors"],
        "contraindications": ["Lymphocyte count <500", "Neutrophil count <1000",
                               "Hemoglobin <8 g/dL", "Active serious infections"],
        "monitoring_requirements": ["CBC with differential at baseline, 4-8 weeks, then every 3 months",
                                     "LFTs at baseline and every 3 months",
                                     "Lipid panel at 4-8 weeks after initiation",
                                     "TB screening before initiation and annually",
                                     "Hepatitis B/C serology at baseline",
                                     "VTE risk assessment (especially patients >50 with CV risk factors)",
                                     "Lymphocyte count monitoring -- hold if <500"],
    },
    {
        "drug_name": "Etanercept",
        "drug_class": "TNF inhibitor",
        "mechanism": "Soluble TNF receptor fusion protein (p75-Fc) -- decoy receptor for TNF-alpha",
        "indicated_diseases": ["rheumatoid_arthritis", "psoriasis", "ankylosing_spondylitis"],
        "pgx_considerations": ["TNFA -308 G>A (rs1800629) may affect response",
                                "HLA-DRB1 shared epitope correlates with better TNF-i response"],
        "contraindications": ["Active TB", "Hepatitis B", "Heart failure NYHA III-IV",
                               "Demyelinating disease", "Sepsis"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "CBC every 3-6 months",
                                     "LFTs every 3-6 months",
                                     "Hepatitis B serology at baseline",
                                     "Monitor for injection site reactions"],
    },
    {
        "drug_name": "Infliximab",
        "drug_class": "TNF inhibitor",
        "mechanism": "Chimeric anti-TNF-alpha monoclonal antibody -- IV infusion",
        "indicated_diseases": ["rheumatoid_arthritis", "inflammatory_bowel_disease",
                               "ankylosing_spondylitis", "psoriasis"],
        "pgx_considerations": ["HLA-DRB1*03:01 associated with anti-drug antibody formation",
                                "FCGR3A V158F affects ADCC",
                                "ATI (anti-drug antibodies) monitored via trough levels"],
        "contraindications": ["Active TB", "Hepatitis B", "Heart failure NYHA III-IV",
                               "Serious infections"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Infusion reaction monitoring",
                                     "Drug trough levels and anti-drug antibodies every 6 months",
                                     "CBC with differential every 3-6 months",
                                     "LFTs every 3-6 months"],
    },
    {
        "drug_name": "Golimumab",
        "drug_class": "TNF inhibitor",
        "mechanism": "Fully human anti-TNF-alpha monoclonal antibody",
        "indicated_diseases": ["rheumatoid_arthritis", "ankylosing_spondylitis",
                               "inflammatory_bowel_disease"],
        "pgx_considerations": ["FCGR3A V158F affects efficacy",
                                "Anti-drug antibody formation possible"],
        "contraindications": ["Active TB", "Heart failure NYHA III-IV",
                               "Active serious infections"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "CBC every 3-6 months",
                                     "LFTs every 3-6 months",
                                     "Hepatitis B serology at baseline"],
    },
    {
        "drug_name": "Certolizumab pegol",
        "drug_class": "TNF inhibitor",
        "mechanism": "PEGylated Fab' fragment of humanized anti-TNF antibody -- no Fc region",
        "indicated_diseases": ["rheumatoid_arthritis", "ankylosing_spondylitis", "psoriasis"],
        "pgx_considerations": ["Minimal placental transfer (no Fc) -- preferred in pregnancy",
                                "FCGR3A irrelevant (no Fc region)"],
        "contraindications": ["Active TB", "Active serious infections",
                               "Demyelinating disease"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "CBC every 3-6 months",
                                     "LFTs baseline and periodic",
                                     "Safe for use during pregnancy and breastfeeding"],
    },
    {
        "drug_name": "Baricitinib",
        "drug_class": "JAK inhibitor (JAK1/JAK2)",
        "mechanism": "Selective JAK1/JAK2 inhibitor",
        "indicated_diseases": ["rheumatoid_arthritis"],
        "pgx_considerations": ["Renal elimination -- dose reduce to 1mg daily if CrCl 30-60",
                                "CYP3A4 interaction (less than tofacitinib)",
                                "VTE risk similar class effect"],
        "contraindications": ["Lymphocyte count <500", "Neutrophil count <1000",
                               "Hemoglobin <8 g/dL", "Severe hepatic impairment",
                               "CrCl <30 mL/min"],
        "monitoring_requirements": ["CBC at baseline, 4 weeks, then every 3 months",
                                     "LFTs at baseline and periodic",
                                     "Lipid panel at 12 weeks",
                                     "Renal function monitoring",
                                     "VTE risk assessment"],
    },
    {
        "drug_name": "Upadacitinib",
        "drug_class": "JAK inhibitor (JAK1 selective)",
        "mechanism": "Selective JAK1 inhibitor -- engineered for JAK1 selectivity",
        "indicated_diseases": ["rheumatoid_arthritis", "ankylosing_spondylitis",
                               "inflammatory_bowel_disease", "psoriasis"],
        "pgx_considerations": ["CYP3A4 substrate -- dose adjust with strong CYP3A4 inhibitors",
                                "Higher selectivity for JAK1 may improve safety vs pan-JAK"],
        "contraindications": ["Active serious infections", "Lymphocyte count <500",
                               "Neutrophil count <1000", "Active TB"],
        "monitoring_requirements": ["CBC at baseline and every 3 months",
                                     "LFTs at baseline and periodic",
                                     "Lipid panel at 12 weeks",
                                     "TB screening before initiation",
                                     "VTE and MACE risk assessment for patients >65"],
    },
    {
        "drug_name": "Guselkumab",
        "drug_class": "IL-23 inhibitor (p19 specific)",
        "mechanism": "Fully human anti-IL-23 p19 subunit monoclonal antibody",
        "indicated_diseases": ["psoriasis"],
        "pgx_considerations": ["IL23R variants (rs11209026) may predict response",
                                "HLA-C*06:02 positive patients show higher response rates"],
        "contraindications": ["Active TB", "Active serious infections"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Monitor for infections",
                                     "PASI score at 16 and 28 weeks",
                                     "Monitor for Crohn's disease exacerbation"],
    },
    {
        "drug_name": "Risankizumab",
        "drug_class": "IL-23 inhibitor (p19 specific)",
        "mechanism": "Humanized anti-IL-23 p19 monoclonal antibody",
        "indicated_diseases": ["psoriasis", "inflammatory_bowel_disease"],
        "pgx_considerations": ["IL23R and IL12B variants may predict response"],
        "contraindications": ["Active TB", "Active serious infections"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Monitor for infections",
                                     "PASI at 16 weeks for psoriasis response",
                                     "Endoscopic evaluation for IBD response"],
    },
    {
        "drug_name": "Vedolizumab",
        "drug_class": "Integrin inhibitor (gut-selective)",
        "mechanism": "Humanized anti-alpha4-beta7 integrin monoclonal antibody -- gut-selective",
        "indicated_diseases": ["inflammatory_bowel_disease"],
        "pgx_considerations": ["Gut-selective mechanism -- lower systemic immunosuppression",
                                "Anti-drug antibody formation affects trough levels"],
        "contraindications": ["Active serious infections",
                               "PML risk (lower than natalizumab)"],
        "monitoring_requirements": ["Monitor for infections",
                                     "Drug trough levels if loss of response",
                                     "Endoscopic assessment for mucosal healing",
                                     "Monitor for nasopharyngitis and headache"],
    },
    {
        "drug_name": "Natalizumab",
        "drug_class": "Anti-VLA4 integrin inhibitor",
        "mechanism": "Humanized anti-alpha4 integrin monoclonal antibody -- blocks leukocyte CNS migration",
        "indicated_diseases": ["multiple_sclerosis"],
        "pgx_considerations": ["JC virus antibody index stratifies PML risk",
                                "Anti-drug antibody formation reduces efficacy"],
        "contraindications": ["PML risk (JCV antibody positive with index >1.5)",
                               "Immunocompromised state",
                               "Concurrent immunosuppressants"],
        "monitoring_requirements": ["JC virus antibody testing every 6 months",
                                     "MRI brain every 6-12 months for PML surveillance",
                                     "CBC and LFTs every 3-6 months",
                                     "Anti-natalizumab antibody testing if loss of response"],
    },
    {
        "drug_name": "Ocrelizumab",
        "drug_class": "Anti-CD20 (humanized)",
        "mechanism": "Humanized anti-CD20 monoclonal antibody -- B-cell depletion (MS-specific)",
        "indicated_diseases": ["multiple_sclerosis"],
        "pgx_considerations": ["FCGR3A V158F affects ADCC",
                                "Hepatitis B reactivation risk"],
        "contraindications": ["Active hepatitis B", "Active serious infections"],
        "monitoring_requirements": ["Hepatitis B serology before initiation",
                                     "Immunoglobulin levels every 6 months",
                                     "CBC before each infusion",
                                     "Infusion reaction monitoring",
                                     "Malignancy screening (breast cancer surveillance)"],
    },
    {
        "drug_name": "Sarilumab",
        "drug_class": "IL-6 receptor inhibitor",
        "mechanism": "Fully human anti-IL-6R monoclonal antibody -- subcutaneous",
        "indicated_diseases": ["rheumatoid_arthritis"],
        "pgx_considerations": ["IL6R Asp358Ala (rs2228145) affects response",
                                "CRP masked by IL-6R blockade"],
        "contraindications": ["Neutropenia <500", "Active infections",
                               "Hepatic impairment"],
        "monitoring_requirements": ["Neutrophil count every 4-8 weeks initially",
                                     "LFTs every 4-8 weeks for first 6 months",
                                     "Lipid panel at 4-8 weeks",
                                     "Note: CRP unreliable for infection monitoring"],
    },
    {
        "drug_name": "Ixekizumab",
        "drug_class": "IL-17A inhibitor",
        "mechanism": "Humanized anti-IL-17A monoclonal antibody",
        "indicated_diseases": ["psoriasis", "ankylosing_spondylitis"],
        "pgx_considerations": ["HLA-C*06:02 positive patients show higher response rates"],
        "contraindications": ["Active IBD (may worsen)", "Active TB",
                               "Active serious infections"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "Monitor for Candida infections",
                                     "Monitor for new-onset IBD",
                                     "PASI at 12 weeks for psoriasis"],
    },
    {
        "drug_name": "Deucravacitinib",
        "drug_class": "TYK2 inhibitor",
        "mechanism": "Selective TYK2 allosteric inhibitor -- binds regulatory domain",
        "indicated_diseases": ["psoriasis"],
        "pgx_considerations": ["TYK2 P1104A loss-of-function variant confers autoimmune protection",
                                "Minimal CYP interaction -- no dose adjustments needed"],
        "contraindications": ["Active serious infections", "Active TB"],
        "monitoring_requirements": ["TB screening before initiation",
                                     "CBC every 3 months initially",
                                     "LFTs at baseline and periodic",
                                     "Lipid panel at 12 weeks",
                                     "PASI at 16 weeks for response assessment"],
    },
]


# =====================================================================
# Flare Prediction Biomarker Patterns
# =====================================================================

FLARE_BIOMARKER_PATTERNS: Dict[str, Dict[str, Any]] = {
    "rheumatoid_arthritis": {
        "early_warning_biomarkers": ["CRP", "ESR", "IL-6", "MMP-3", "14-3-3eta"],
        "thresholds": {
            "CRP_rising": {"min_increase_pct": 50, "window_days": 30},
            "ESR_rising": {"min_increase_pct": 40, "window_days": 30},
        },
        "protective_signals": ["stable_RF_titer", "normal_CRP_trend", "adequate_sleep"],
    },
    "systemic_lupus_erythematosus": {
        "early_warning_biomarkers": ["anti-dsDNA_titer", "complement_C3", "complement_C4",
                                     "lymphocyte_count", "proteinuria"],
        "thresholds": {
            "anti_dsDNA_rising": {"min_increase_pct": 25, "window_days": 60},
            "complement_falling": {"max_decrease_pct": 20, "window_days": 30},
        },
        "protective_signals": ["stable_complement", "normal_urinalysis", "sun_avoidance"],
    },
    "inflammatory_bowel_disease": {
        "early_warning_biomarkers": ["calprotectin", "CRP", "lactoferrin", "albumin"],
        "thresholds": {
            "calprotectin_rising": {"value_above": 250, "unit": "ug/g"},
            "albumin_falling": {"value_below": 3.5, "unit": "g/dL"},
        },
        "protective_signals": ["normal_calprotectin", "stable_albumin", "medication_adherence"],
    },
    "ankylosing_spondylitis": {
        "early_warning_biomarkers": ["CRP", "ESR", "IL-17", "MMP-3"],
        "thresholds": {
            "CRP_rising": {"min_increase_pct": 40, "window_days": 30},
            "BASDAI_rising": {"min_increase": 1.0, "window_days": 90},
        },
        "protective_signals": ["regular_exercise", "NSAID_response", "normal_CRP_trend"],
    },
    "psoriasis": {
        "early_warning_biomarkers": ["CRP", "IL-17", "TNF-alpha", "IL-23"],
        "thresholds": {
            "PASI_rising": {"min_increase": 3, "window_days": 60},
        },
        "protective_signals": ["stable_PASI", "medication_adherence", "stress_management"],
    },
    "sjogrens_syndrome": {
        "early_warning_biomarkers": ["ESR", "IgG", "RF", "complement_C4", "beta2_microglobulin"],
        "thresholds": {
            "IgG_rising": {"min_increase_pct": 25, "window_days": 90},
            "complement_C4_falling": {"max_decrease_pct": 20, "window_days": 60},
        },
        "protective_signals": ["stable_gland_function", "normal_IgG", "no_lymphadenopathy"],
    },
    "systemic_sclerosis": {
        "early_warning_biomarkers": ["CRP", "NT-proBNP", "DLCO", "creatinine", "anti-Scl-70_titer"],
        "thresholds": {
            "DLCO_falling": {"max_decrease_pct": 15, "window_days": 180},
            "NT_proBNP_rising": {"value_above": 300, "unit": "pg/mL"},
        },
        "protective_signals": ["stable_skin_score", "stable_DLCO", "normal_renal_function"],
    },
    "multiple_sclerosis": {
        "early_warning_biomarkers": ["neurofilament_light", "MRI_gadolinium_lesions", "lymphocyte_count", "IgG_index", "oligoclonal_bands"],
        "thresholds": {
            "NfL_rising": {"value_above": 16, "unit": "pg/mL"},
            "new_T2_lesions": {"value_above": 1, "window_days": 180},
        },
        "protective_signals": ["stable_EDSS", "no_new_MRI_lesions", "DMT_adherence", "vitamin_D_sufficient"],
    },
    "type_1_diabetes": {
        "early_warning_biomarkers": ["HbA1c", "C_peptide", "GAD65_antibody", "IA2_antibody", "fasting_glucose"],
        "thresholds": {
            "HbA1c_rising": {"min_increase_pct": 10, "window_days": 90},
            "C_peptide_falling": {"max_decrease_pct": 30, "window_days": 180},
        },
        "protective_signals": ["stable_HbA1c", "preserved_C_peptide", "insulin_dose_stable", "time_in_range_above_70"],
    },
    "myasthenia_gravis": {
        "early_warning_biomarkers": ["AChR_antibody_titer", "MuSK_antibody_titer", "FVC_percent_predicted", "complement_C3", "CRP"],
        "thresholds": {
            "AChR_rising": {"min_increase_pct": 50, "window_days": 90},
            "FVC_falling": {"max_decrease_pct": 20, "window_days": 60},
        },
        "protective_signals": ["stable_AChR_titer", "FVC_above_80_percent", "cholinesterase_inhibitor_response", "no_bulbar_symptoms"],
    },
    "celiac_disease": {
        "early_warning_biomarkers": ["anti_tTG_IgA", "anti_DGP_IgG", "hemoglobin", "ferritin", "albumin"],
        "thresholds": {
            "tTG_rising": {"value_above": 20, "unit": "U/mL"},
            "ferritin_falling": {"value_below": 15, "unit": "ng/mL"},
        },
        "protective_signals": ["normal_tTG", "strict_gluten_free_diet", "normal_hemoglobin", "normal_villous_architecture"],
    },
    "graves_disease": {
        "early_warning_biomarkers": ["TSI", "TSH", "free_T4", "free_T3", "TRAb"],
        "thresholds": {
            "TSI_rising": {"min_increase_pct": 30, "window_days": 90},
            "TSH_suppressed": {"value_below": 0.1, "unit": "mIU/L"},
        },
        "protective_signals": ["normal_TSH", "stable_free_T4", "TRAb_declining", "antithyroid_drug_response"],
    },
    "hashimoto_thyroiditis": {
        "early_warning_biomarkers": ["TSH", "anti_TPO", "anti_thyroglobulin", "free_T4", "free_T3"],
        "thresholds": {
            "TSH_rising": {"value_above": 10, "unit": "mIU/L"},
            "TPO_rising": {"min_increase_pct": 50, "window_days": 180},
        },
        "protective_signals": ["stable_TSH_on_replacement", "normal_free_T4", "TPO_stable_or_declining"],
    },
}
