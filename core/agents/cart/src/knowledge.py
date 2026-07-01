"""CAR-T Intelligence Agent — Knowledge Graph.

Extends the Clinker pattern from core/engines/precision-intelligence/src/knowledge.py,
adapted for the CAR-T cell therapy domain. Contains:

1. CART_TARGETS: 34 target antigens with clinical data
2. CART_TOXICITIES: 17 toxicity profiles with grading/management
3. CART_MANUFACTURING: 20 manufacturing process parameters
4. CART_BIOMARKERS: 23 biomarker profiles
5. CART_REGULATORY: 6 approved product regulatory records
6. CART_IMMUNOGENICITY: 6 immunogenicity topic profiles

Author: Adam Jones
Date: February 2026
"""

from typing import Any, Dict, Optional


# =============================================================================
# 1. CART_TARGETS — Target antigen knowledge graph (~33 entries)
# =============================================================================

CART_TARGETS: Dict[str, Dict[str, Any]] = {
    "CD19": {
        "protein": "B-Lymphocyte Antigen CD19",
        "uniprot_id": "P15391",
        "expression": "B-cell lineage (pro-B to mature B, not plasma cells)",
        "diseases": ["B-ALL", "DLBCL", "FL", "MCL", "CLL"],
        "approved_products": [
            "Kymriah (tisagenlecleucel)",
            "Yescarta (axicabtagene ciloleucel)",
            "Tecartus (brexucabtagene autoleucel)",
            "Breyanzi (lisocabtagene maraleucel)",
        ],
        "key_trials": ["ELIANA", "ZUMA-1", "ZUMA-2", "TRANSFORM", "TRANSCEND"],
        "known_resistance": [
            "CD19 loss/mutation", "lineage switch", "trogocytosis",
            "alternative splicing (exon 2 deletion)",
        ],
        "toxicity_profile": {"CRS": "30-90%", "ICANS": "20-65%", "B_cell_aplasia": "expected"},
        "normal_tissue": "B-cells (acceptable on-target/off-tumor)",
    },
    "BCMA": {
        "protein": "B-Cell Maturation Antigen (TNFRSF17)",
        "uniprot_id": "Q02223",
        "expression": "Mature B-cells and plasma cells",
        "diseases": ["Multiple Myeloma", "Waldenstrom's Macroglobulinemia"],
        "approved_products": [
            "Abecma (idecabtagene vicleucel)",
            "Carvykti (ciltacabtagene autoleucel)",
        ],
        "key_trials": ["KarMMa", "CARTITUDE-1", "CARTITUDE-4"],
        "known_resistance": [
            "BCMA downregulation", "biallelic BCMA loss",
            "soluble BCMA shedding (gamma-secretase)",
        ],
        "toxicity_profile": {"CRS": "80-95%", "ICANS": "10-20%", "cytopenia": "common"},
        "normal_tissue": "Mature B-cells, plasma cells",
    },
    "CD22": {
        "protein": "B-cell receptor CD22 (SIGLEC-2)",
        "uniprot_id": "P20273",
        "expression": "B-cell lineage, broader than CD19",
        "diseases": ["B-ALL", "DLBCL", "Hairy cell leukemia"],
        "approved_products": [],
        "key_trials": ["NCT02315612", "NCT04088890"],
        "known_resistance": ["CD22 downregulation", "diminished site density"],
        "toxicity_profile": {"CRS": "moderate", "ICANS": "low"},
        "normal_tissue": "B-cells",
    },
    "CD20": {
        "protein": "B-lymphocyte antigen CD20 (MS4A1)",
        "uniprot_id": "P11836",
        "expression": "Pre-B to mature B-cells",
        "diseases": ["DLBCL", "FL", "MCL", "CLL"],
        "approved_products": [],
        "key_trials": ["NCT03277729"],
        "known_resistance": ["CD20 loss (rare)", "rituximab resistance"],
        "toxicity_profile": {"CRS": "moderate", "B_cell_aplasia": "expected"},
        "normal_tissue": "B-cells",
    },
    "CD30": {
        "protein": "TNFRSF8 (Ki-1 antigen)",
        "uniprot_id": "P28908",
        "expression": "Activated T/B cells, Reed-Sternberg cells",
        "diseases": ["Hodgkin Lymphoma", "ALCL", "PTCL"],
        "approved_products": [],
        "key_trials": ["RELY-30", "NCT02690545"],
        "known_resistance": ["CD30 shedding"],
        "toxicity_profile": {"CRS": "moderate"},
        "normal_tissue": "Activated lymphocytes (limited)",
    },
    "CD33": {
        "protein": "Siglec-3",
        "uniprot_id": "P20138",
        "expression": "Myeloid progenitors, AML blasts",
        "diseases": ["AML", "MDS"],
        "approved_products": [],
        "key_trials": ["NCT03126864", "NCT03971799"],
        "known_resistance": ["CD33 splice variants", "heterogeneous expression"],
        "toxicity_profile": {"CRS": "moderate", "myelosuppression": "significant"},
        "normal_tissue": "Myeloid progenitors (significant on-target/off-tumor)",
    },
    "CD38": {
        "protein": "ADP-ribosyl cyclase/cyclic ADP-ribose hydrolase 1",
        "uniprot_id": "P28907",
        "expression": "Plasma cells, activated lymphocytes",
        "diseases": ["Multiple Myeloma", "T-ALL"],
        "approved_products": [],
        "key_trials": ["NCT03464916"],
        "known_resistance": ["CD38 downregulation"],
        "toxicity_profile": {"CRS": "moderate"},
        "normal_tissue": "T-cells, NK cells, red blood cells (broad)",
    },
    "CD123": {
        "protein": "CD123 (IL-3 receptor alpha chain, IL3RA)",
        "uniprot_id": "P26951",
        "expression": "Leukemic blasts, plasmacytoid dendritic cells, basophils; low on normal HSCs",
        "diseases": ["AML", "BPDCN", "B-ALL (CD123+)", "Systemic mastocytosis", "Hairy cell leukemia"],
        "approved_products": [],
        "key_trials": [
            "NCT02159495 (CD123 CAR-T, City of Hope)",
            "NCT04106076 (UCART123, Cellectis)",
            "NCT03190278 (MB-102, Mustang Bio)",
            "NCT04678336 (UniCAR-T-CD123)",
        ],
        "known_resistance": [
            "CD123 downregulation under selective pressure",
            "Myeloid lineage plasticity",
            "HSC toxicity limiting dosing",
            "TME immunosuppression in AML",
        ],
        "toxicity_profile": {
            "CRS": "40-70%",
            "myelosuppression": "60-90% (on-target HSC toxicity)",
            "ICANS": "10-25%",
            "prolonged_cytopenias": "50-80%",
        },
        "normal_tissue": "Expressed on normal HSCs at low levels — significant myelosuppression risk; may require bridge to transplant",
    },
    "GD2": {
        "protein": "GD2 (disialoganglioside)",
        "uniprot_id": None,
        "expression": "Neuroblastoma, melanoma, retinoblastoma; CNS neurons, peripheral nerves",
        "diseases": ["Neuroblastoma", "Melanoma", "Retinoblastoma", "DIPG/DMG", "Osteosarcoma", "Small cell lung cancer"],
        "approved_products": [],
        "key_trials": [
            "NCT03635632 (GD2 CAR-T, Baylor)",
            "NCT04099797 (GD2-APRIL CAR-T)",
            "NCT04196413 (C7R-GD2 CAR-T armored)",
            "NCT05298995 (GD2 CAR-T DIPG, Stanford)",
        ],
        "known_resistance": [
            "GD2 shedding by tumor cells",
            "TME immunosuppression (TGF-beta, IL-10)",
            "Poor T-cell infiltration in solid tumors",
            "Neuropathic pain from on-target off-tumor",
        ],
        "toxicity_profile": {
            "pain_syndrome": "30-60% (neuropathic)",
            "CRS": "20-50%",
            "ICANS": "5-15%",
            "GvHD": "if allogeneic",
        },
        "normal_tissue": "Expressed on CNS neurons and peripheral nerves — neuropathic pain is dose-limiting; careful dose escalation required",
    },
    "HER2": {
        "protein": "Human Epidermal Growth Factor Receptor 2 (ERBB2)",
        "uniprot_id": "P04626",
        "expression": "Epithelial cells, overexpressed in tumors",
        "diseases": ["Breast Cancer", "Gastric Cancer", "Sarcoma", "GBM"],
        "approved_products": [],
        "key_trials": ["NCT00902044", "NCT03500991"],
        "known_resistance": ["HER2 heterogeneity", "solid tumor barriers"],
        "toxicity_profile": {"CRS": "moderate", "on_target_off_tumor": "cardiac/pulmonary risk"},
        "normal_tissue": "Heart, lung, GI epithelium (safety concern)",
    },
    "GPC3": {
        "protein": "Glypican-3",
        "uniprot_id": "P51654",
        "expression": "Hepatocellular carcinoma, some pediatric tumors",
        "diseases": ["HCC", "Hepatoblastoma", "Wilms tumor"],
        "approved_products": [],
        "key_trials": ["NCT02905188", "NCT03884751"],
        "known_resistance": ["Solid tumor microenvironment"],
        "toxicity_profile": {"CRS": "moderate"},
        "normal_tissue": "Fetal liver (low in adults)",
    },
    "EGFR": {
        "protein": "Epidermal Growth Factor Receptor",
        "uniprot_id": "P00533",
        "expression": "Epithelial cells, overexpressed in many solid tumors",
        "diseases": ["NSCLC", "GBM", "Head and Neck Cancer"],
        "approved_products": [],
        "key_trials": ["NCT03182816"],
        "known_resistance": ["Antigen heterogeneity", "solid tumor barriers"],
        "toxicity_profile": {"on_target_off_tumor": "skin/GI toxicity risk"},
        "normal_tissue": "Skin, GI tract, lung epithelium",
    },
    "EGFRvIII": {
        "protein": "EGFR variant III (tumor-specific mutant)",
        "uniprot_id": None,
        "expression": "GBM (20-30%), some NSCLC",
        "diseases": ["Glioblastoma"],
        "approved_products": [],
        "key_trials": ["NCT02209376", "NCT01454596"],
        "known_resistance": ["Antigen loss (heterogeneous expression)"],
        "toxicity_profile": {"CRS": "low-moderate"},
        "normal_tissue": "None (tumor-specific neoantigen)",
    },
    "Mesothelin": {
        "protein": "Mesothelin (MSLN)",
        "uniprot_id": "Q13421",
        "expression": "Mesothelial cells, overexpressed in mesothelioma/pancreatic/ovarian",
        "diseases": ["Mesothelioma", "Pancreatic Cancer", "Ovarian Cancer"],
        "approved_products": [],
        "key_trials": ["NCT02159716", "NCT03054298"],
        "known_resistance": ["Immunosuppressive TME", "antigen shedding"],
        "toxicity_profile": {"CRS": "moderate", "on_target_off_tumor": "pleuritis risk"},
        "normal_tissue": "Pleura, peritoneum, pericardium",
    },
    "Claudin18.2": {
        "protein": "Claudin 18.2 (CLDN18.2, tight junction protein)",
        "uniprot_id": "P56856",
        "expression": "Gastric mucosa (normal), overexpressed in gastric, pancreatic, esophageal cancers",
        "diseases": ["Gastric/GEJ adenocarcinoma", "Pancreatic ductal adenocarcinoma", "Esophageal cancer", "Mucinous ovarian cancer"],
        "approved_products": [],
        "key_trials": [
            "NCT03874897 (CT041, CARsgen)",
            "NCT04404595 (Claudin18.2 CAR-T, Nanjing Legend)",
            "NCT05393986 (CLDN18.2 CAR-T, Innovent)",
        ],
        "known_resistance": [
            "Claudin18.2 expression heterogeneity",
            "Desmoplastic stroma limiting T-cell access",
            "TGF-beta-rich TME",
            "Tight junction shielding of epitope",
        ],
        "toxicity_profile": {"CRS": "40-70%", "GI_toxicity": "15-30%", "mucositis": "10-25%", "ICANS": "5-10%"},
        "normal_tissue": "Normal gastric mucosa expression — GI toxicity expected but generally manageable; expression becomes accessible only in malignant transformation (tight junction disruption)",
    },
    "MUC1": {
        "protein": "Mucin 1 (episialin)",
        "uniprot_id": "P15941",
        "expression": "Epithelial cells, aberrantly glycosylated in tumors",
        "diseases": ["Breast Cancer", "Pancreatic Cancer", "NSCLC"],
        "approved_products": [],
        "key_trials": ["NCT02587689"],
        "known_resistance": ["Glycosylation heterogeneity"],
        "toxicity_profile": {"on_target_off_tumor": "epithelial toxicity risk"},
        "normal_tissue": "All epithelial surfaces",
    },
    "PSMA": {
        "protein": "Prostate-Specific Membrane Antigen (FOLH1)",
        "uniprot_id": "Q04609",
        "expression": "Prostate epithelium, overexpressed in prostate cancer",
        "diseases": ["Prostate Cancer"],
        "approved_products": [],
        "key_trials": ["NCT04227275", "NCT03089203"],
        "known_resistance": ["Heterogeneous expression", "solid tumor barriers"],
        "toxicity_profile": {"CRS": "moderate"},
        "normal_tissue": "Prostate, kidney, small intestine, brain",
    },
    "ROR1": {
        "protein": "Receptor Tyrosine Kinase-Like Orphan Receptor 1",
        "uniprot_id": "Q01973",
        "expression": "Embryonic, re-expressed in CLL/MCL/TNBC",
        "diseases": ["CLL", "MCL", "TNBC", "NSCLC"],
        "approved_products": [],
        "key_trials": ["NCT02706392"],
        "known_resistance": ["Variable expression levels"],
        "toxicity_profile": {"CRS": "moderate"},
        "normal_tissue": "Minimal in adults (oncofetal antigen)",
    },
    "GPRC5D": {
        "protein": "GPRC5D (G protein-coupled receptor class C group 5 member D)",
        "uniprot_id": "Q9NZD1",
        "expression": "Myeloma cells, hair follicles; limited normal tissue",
        "diseases": ["Multiple myeloma (including BCMA-refractory)", "Light chain amyloidosis"],
        "approved_products": ["Talquetamab (bispecific mAb, not CAR-T)"],
        "key_trials": [
            "NCT04555551 (GPRC5D CAR-T, MSKCC)",
            "NCT05016778 (OriCAR-017)",
            "NCT05739188 (GPRC5D CAR-T, BMS)",
        ],
        "known_resistance": [
            "GPRC5D downregulation (less common than BCMA)",
            "Dual antigen loss (BCMA+GPRC5D)",
            "Immune evasion via PD-L1",
        ],
        "toxicity_profile": {
            "CRS": "60-80%",
            "ICANS": "10-20%",
            "nail_toxicity": "40-60%",
            "skin_toxicity": "30-50%",
            "dysgeusia": "40-70%",
        },
        "normal_tissue": "Hair follicle expression causes reversible nail/skin toxicity and taste changes; manageable and dose-dependent",
    },
    "IL13Ra2": {
        "protein": "Interleukin-13 Receptor Alpha 2",
        "uniprot_id": "Q14627",
        "expression": "GBM (overexpressed in >50% of tumors)",
        "diseases": ["Glioblastoma", "Medulloblastoma"],
        "approved_products": [],
        "key_trials": ["NCT02208362", "NCT04510051"],
        "known_resistance": ["Antigen heterogeneity within tumor"],
        "toxicity_profile": {"CRS": "low-moderate"},
        "normal_tissue": "Testis, minimal elsewhere",
    },
    "DLL3": {
        "protein": "DLL3 (Delta-like ligand 3, Notch pathway)",
        "uniprot_id": "Q9NYJ7",
        "expression": "SCLC (>80% of tumors), large cell neuroendocrine carcinoma; minimal normal tissue",
        "diseases": ["SCLC", "Large cell neuroendocrine carcinoma", "High-grade neuroendocrine tumors", "Neuroendocrine prostate cancer"],
        "approved_products": ["Tarlatamab (bispecific, not CAR-T)"],
        "key_trials": [
            "NCT05680922 (DLL3 CAR-T, BMS)",
            "NCT05471609 (DLL3 x CD3 bispecific)",
            "NCT05680922 (DLL3 armored CAR-T)",
        ],
        "known_resistance": [
            "DLL3 internalization reducing surface density",
            "Neuroendocrine-to-mesenchymal transition",
            "Rapid disease progression outpacing manufacturing",
        ],
        "toxicity_profile": {"CRS": "30-55%", "ICANS": "5-15%", "fatigue": "40-60%"},
        "normal_tissue": "Minimal normal tissue expression (intracellular in Notch pathway) — favorable safety profile expected",
    },
    "B7-H3": {
        "protein": "B7-H3 (CD276, immune checkpoint molecule)",
        "uniprot_id": "Q5ZPR3",
        "expression": "Broadly overexpressed in solid tumors; limited normal tissue",
        "diseases": ["Neuroblastoma", "Rhabdomyosarcoma", "Ewing sarcoma", "Osteosarcoma", "DIPG/DMG", "Prostate cancer", "NSCLC", "Ovarian cancer"],
        "approved_products": [],
        "key_trials": [
            "NCT04483778 (B7-H3 CAR-T, Seattle Children's)",
            "NCT04897321 (B7-H3 CAR-T + checkpoint)",
            "NCT04670068 (B7-H3 intracranial, Stanford)",
        ],
        "known_resistance": [
            "B7-H3 expression heterogeneity in CNS tumors",
            "Immunosuppressive TME",
            "Blood-brain barrier for CNS tumors",
            "T-cell exhaustion in solid tumors",
        ],
        "toxicity_profile": {"CRS": "20-45%", "ICANS": "5-15%", "musculoskeletal_pain": "10-20%"},
        "normal_tissue": "Limited normal tissue expression — favorable safety profile; one of the most promising pan-solid tumor CAR-T targets",
    },
    "NKG2D_ligands": {
        "protein": "NKG2D ligands (MICA, MICB, ULBP1-6, stress-induced surface molecules)",
        "uniprot_id": None,
        "expression": "Upregulated on stressed/transformed cells; minimal on healthy tissue",
        "diseases": ["AML", "Multiple myeloma", "Ovarian cancer", "Colorectal cancer", "Hepatocellular carcinoma"],
        "approved_products": [],
        "key_trials": [
            "NCT03018405 (NKG2D CAR-T, Celyad)",
            "NCT04167696 (CYAD-101, allogeneic NKG2D)",
            "NCT03370198 (NKG2D CAR-T, AML)",
        ],
        "known_resistance": [
            "NKG2D ligand shedding (soluble MICA/MICB)",
            "TGF-beta downregulation of NKG2D",
            "Metalloprotease cleavage of ligands",
            "Immune checkpoint upregulation",
        ],
        "toxicity_profile": {"CRS": "15-40%", "autoimmune_risk": "theoretical (stress ligand on normal cells)", "ICANS": "5-10%"},
        "normal_tissue": "Stress ligands theoretically on any stressed normal cell — autoimmune risk is theoretical but not observed clinically at therapeutic doses",
    },
    "CD7": {
        "protein": "T-cell antigen CD7",
        "uniprot_id": "P09564",
        "expression": "T-cells, NK cells, T-ALL blasts",
        "diseases": ["T-ALL", "T-lymphoma"],
        "approved_products": [],
        "key_trials": ["NCT04572308"],
        "known_resistance": ["Fratricide (T-cell target on T-cell product)"],
        "toxicity_profile": {"CRS": "moderate-high", "T_cell_aplasia": "risk"},
        "normal_tissue": "Normal T-cells and NK cells",
    },
    "CD5": {
        "protein": "Lymphocyte Antigen CD5",
        "uniprot_id": "P06127",
        "expression": "T-cells, subset of B-cells",
        "diseases": ["T-ALL", "T-lymphoma", "CLL"],
        "approved_products": [],
        "key_trials": ["NCT03081910"],
        "known_resistance": ["Fratricide"],
        "toxicity_profile": {"CRS": "moderate"},
        "normal_tissue": "Normal T-cells",
    },
    "FcRH5": {
        "protein": "FcRH5 (Fc receptor-homolog 5, FCRLA)",
        "uniprot_id": "Q96RD9",
        "expression": "B cells, plasma cells, myeloma cells; broader B-cell lineage than BCMA",
        "diseases": ["Multiple myeloma", "B-cell lymphomas"],
        "approved_products": [],
        "key_trials": [
            "NCT05338775 (FcRH5 CAR-T)",
            "MCARH109 (Cevostamab, bispecific mAb)",
        ],
        "known_resistance": [
            "Broader expression may reduce selectivity",
            "B-cell aplasia expected",
            "Potential for antigen masking",
        ],
        "toxicity_profile": {"CRS": "50-70%", "B_cell_aplasia": "near 100%", "ICANS": "10-25%"},
        "normal_tissue": "Expressed on normal B cells — expect B-cell aplasia similar to CD19-targeted therapies",
    },
    "SLAMF7": {
        "protein": "SLAM Family Member 7 (CS1/CD319)",
        "uniprot_id": "Q9NQ25",
        "expression": "Plasma cells, NK cells, activated T-cells",
        "diseases": ["Multiple Myeloma"],
        "approved_products": [],
        "key_trials": ["UCART-CS1 Phase 1", "CARAMBA (Phase 1/2)"],
        "known_resistance": ["SLAMF7 shedding", "NK-cell fratricide risk"],
        "toxicity_profile": {"crs_any": "60-80%", "crs_grade3": "5-10%", "nk_depletion": "Expected"},
        "normal_tissue": "NK cells, activated T-cells — risk of NK/T depletion",
    },
    "CD70": {
        "protein": "CD70 (TNFSF7/CD27L)",
        "uniprot_id": "P32970",
        "expression": "Activated lymphocytes, tumor cells (AML, RCC)",
        "diseases": ["AML", "T-ALL", "Renal Cell Carcinoma", "Glioblastoma"],
        "approved_products": [],
        "key_trials": ["PRGN-3006 Phase 1/2 (AML)", "CD70-targeting Phase 1 (RCC)"],
        "known_resistance": ["CD70 internalization", "Immunosuppressive TME in solid tumors"],
        "toxicity_profile": {"crs_any": "50-70%", "crs_grade3": "5-15%", "t_cell_fratricide": "Risk in allogeneic"},
        "normal_tissue": "Activated lymphocytes — controlled by limited expression kinetics",
    },
    "TROP2": {
        "protein": "TROP2 (Trophoblast cell surface antigen 2, TACSTD2)",
        "uniprot_id": "P09758",
        "expression": "Epithelial cancers broadly (breast, NSCLC, urothelial, HNSCC); low-level normal epithelium",
        "diseases": ["Triple-negative breast cancer", "NSCLC (non-squamous)", "Urothelial carcinoma", "HNSCC", "Endometrial cancer"],
        "approved_products": ["Sacituzumab govitecan (ADC, not CAR-T)"],
        "key_trials": [
            "NCT05451849 (TROP2 CAR-T)",
            "NCT05671029 (TROP2 armored CAR-T)",
        ],
        "known_resistance": [
            "TROP2 expression variability across tumor regions",
            "Normal epithelial toxicity",
            "Solid tumor TME barriers",
            "Poor T-cell persistence in solid tumors",
        ],
        "toxicity_profile": {"CRS": "20-50%", "skin_toxicity": "15-30%", "mucosal_toxicity": "10-25%", "diarrhea": "10-20%"},
        "normal_tissue": "Low-level normal epithelial expression — mucosal/skin toxicity expected; dose optimization critical",
    },
    "FLT3": {
        "protein": "FMS-Like Tyrosine Kinase 3 (CD135)",
        "uniprot_id": "P36888",
        "expression": "Hematopoietic progenitors, AML blasts",
        "diseases": ["AML", "B-ALL", "MDS"],
        "approved_products": [],
        "key_trials": ["FLT3-CAR-T Phase 1 (UPenn)", "AMG 553 (bispecific, validates target)"],
        "known_resistance": ["FLT3 splice variants", "Myeloablative on-target toxicity"],
        "toxicity_profile": {"myeloablation": "Expected (stem cell rescue needed)", "crs_any": "60-80%"},
        "normal_tissue": "Hematopoietic stem cells — myeloablation expected, requires HSCT rescue",
    },
    "CLL1": {
        "protein": "CLL1 (C-type lectin-like molecule 1, CLEC12A)",
        "uniprot_id": "Q5QGZ9",
        "expression": "Leukemic stem cells, monocytes; absent from normal HSCs",
        "diseases": ["AML", "MDS"],
        "approved_products": [],
        "key_trials": [
            "NCT04219163 (CLL1 CAR-T, China)",
            "NCT04789408 (CLL1/CD33 bispecific)",
            "NCT05252572 (CLL1 CAR-T, Wuhan)",
        ],
        "known_resistance": [
            "CLL1 expression heterogeneity",
            "Monocyte/macrophage depletion",
            "AML clonal evolution",
        ],
        "toxicity_profile": {"CRS": "30-60%", "myelosuppression": "20-40%", "ICANS": "5-15%"},
        "normal_tissue": "Absent from normal HSCs — potentially safer AML target than CD33/CD123",
    },
    "CD44v6": {
        "protein": "CD44 variant 6 (CD44v6, adhesion molecule splice variant)",
        "uniprot_id": "P16070",
        "expression": "AML blasts, head/neck SCC, pancreatic cancer, colorectal cancer; normal keratinocytes and monocytes",
        "diseases": ["AML", "Head and neck SCC", "Pancreatic cancer", "Colorectal cancer"],
        "approved_products": [],
        "key_trials": [
            "NCT04097301 (CD44v6 CAR-T, Italy)",
            "NCT04427449 (CD44v6 CAR-T, AML/solid)",
        ],
        "known_resistance": [
            "CD44v6 splice variant switching",
            "Keratinocyte toxicity",
            "Monocyte depletion",
        ],
        "toxicity_profile": {
            "CRS": "30-60%",
            "skin_toxicity": "30-50% (keratinocyte expression)",
            "monocytopenia": "expected",
            "ICANS": "10-20%",
        },
        "normal_tissue": "Keratinocyte and monocyte expression — skin toxicity and monocyte depletion expected; dose-finding critical",
    },
    "EpCAM": {
        "protein": "Epithelial Cell Adhesion Molecule (CD326/TACSTD1)",
        "uniprot_id": "P16422",
        "expression": "Epithelial cells, overexpressed in carcinomas",
        "diseases": ["Colorectal Cancer", "Gastric Cancer", "Ovarian Cancer", "Hepatocellular Carcinoma", "Pancreatic Cancer"],
        "approved_products": [],
        "key_trials": ["EpCAM-CAR-T Phase 1 (HCC)", "EpCAM-CAR-T Phase 1/2 (peritoneal carcinomatosis)"],
        "known_resistance": ["Epithelial-mesenchymal transition", "Antigen heterogeneity", "Solid tumor barriers"],
        "toxicity_profile": {"on_target_off_tumor": "GI epithelial toxicity risk", "crs_any": "30-50%"},
        "normal_tissue": "Normal GI, hepatic, pancreatic epithelia — significant on-target/off-tumor risk",
    },
    "LNGFR": {
        "protein": "LNGFR (Low-affinity nerve growth factor receptor, CD271, p75NTR)",
        "expression": "Cancer stem cells in melanoma, glioma, squamous cancers; normal neural crest derivatives",
        "diseases": ["Melanoma (cancer stem cells)", "Glioblastoma", "Squamous cell carcinoma", "HNSCC"],
        "approved_products": [],
        "key_trials": [
            "Preclinical stage — no active CAR-T trials",
            "Multiple academic preclinical programs",
        ],
        "known_resistance": [
            "Cancer stem cell plasticity",
            "Low surface density",
            "Neural crest normal tissue toxicity risk",
        ],
        "toxicity_profile": {"neurotoxicity": "theoretical (neural crest expression)", "CRS": "estimated 20-40%"},
        "normal_tissue": "Expressed on neural crest-derived cells — neurotoxicity risk; preclinical only",
    },
}


# =============================================================================
# 2. CART_TOXICITIES — Toxicity knowledge graph (~12 profiles)
# =============================================================================

CART_TOXICITIES: Dict[str, Dict[str, Any]] = {
    "CRS": {
        "full_name": "Cytokine Release Syndrome",
        "mechanism": (
            "Massive cytokine release (IL-6, IFN-gamma, IL-2R, TNF-alpha, GM-CSF) "
            "from activated CAR-T cells and bystander monocytes/macrophages"
        ),
        "grading": {
            "1": "Fever only (>=38C)",
            "2": "Hypotension responsive to fluids; O2 by low-flow nasal cannula",
            "3": "Hypotension requiring vasopressors; O2 by high-flow or non-rebreather",
            "4": "Life-threatening; vasopressors + positive pressure ventilation",
        },
        "grading_system": "Lee 2014 / ASTCT 2019 consensus",
        "incidence": "50-95% any grade; 10-25% grade 3+",
        "timing": "Typically onset day 1-14 post-infusion, peak day 3-7",
        "management": [
            "Tocilizumab (IL-6R blockade) — first-line for grade 2+",
            "Corticosteroids (dexamethasone 10mg IV) — second-line or grade 3+",
            "Siltuximab (anti-IL-6) — if tocilizumab-refractory",
            "Anakinra (IL-1R antagonist) — emerging evidence",
            "Ruxolitinib (JAK inhibitor) — refractory cases",
            "Supportive care: fluids, vasopressors, O2",
        ],
        "biomarkers": ["Ferritin", "CRP", "IL-6", "IFN-gamma", "sIL-2R"],
        "risk_factors": [
            "High tumor burden", "High CAR-T cell dose",
            "CD28 costimulation (vs 4-1BB)", "ALL > lymphoma",
        ],
    },
    "ICANS": {
        "full_name": "Immune Effector Cell-Associated Neurotoxicity Syndrome",
        "mechanism": (
            "BBB disruption from endothelial activation by cytokines (IL-6, IFN-gamma); "
            "direct CAR-T cell CNS trafficking; cerebral edema"
        ),
        "grading": {
            "1": "ICE score 7-9 (mild impairment)",
            "2": "ICE score 3-6 (moderate impairment)",
            "3": "ICE score 0-2; any seizure; raised ICP not requiring intervention",
            "4": "Prolonged seizure; cerebral edema; coma",
        },
        "grading_system": "ASTCT 2019 consensus (ICE score)",
        "incidence": "20-65% any grade; 10-30% grade 3+",
        "timing": "Onset typically day 3-10, often after CRS peak",
        "management": [
            "Corticosteroids (dexamethasone 10mg IV q6h) — first-line",
            "Levetiracetam — seizure prophylaxis",
            "Tocilizumab — if concurrent CRS (may not help isolated ICANS)",
            "Siltuximab — if steroid-refractory",
            "ICU monitoring for grade 3+",
        ],
        "biomarkers": ["ICE score", "CSF protein", "CSF IL-6", "MRI (edema)"],
        "risk_factors": [
            "High tumor burden", "Pre-existing neurologic conditions",
            "Severe CRS", "Thrombocytopenia", "CD28 costimulation",
        ],
    },
    "B_CELL_APLASIA": {
        "full_name": "B-Cell Aplasia / Hypogammaglobulinemia",
        "mechanism": (
            "On-target, off-tumor effect of CD19/CD22-directed CAR-T cells; "
            "depletion of normal B-cell compartment leading to agammaglobulinemia"
        ),
        "grading": {
            "expected": "Anticipated pharmacodynamic effect of CD19 CAR-T",
        },
        "incidence": "Near 100% with CD19 CAR-T; duration correlates with persistence",
        "timing": "Onset with CAR-T engraftment; may persist months to years",
        "management": [
            "IVIG replacement (target IgG > 400 mg/dL)",
            "Infection prophylaxis (PJP, fungal, viral)",
            "Monitoring: quantitative immunoglobulins monthly",
            "Vaccination (live vaccines contraindicated)",
        ],
        "biomarkers": ["IgG levels", "B-cell count (CD19+ cells)", "CAR-T persistence (qPCR)"],
        "risk_factors": ["CD19 CAR-T (expected)", "Prolonged CAR-T persistence"],
    },
    "HLH_MAS": {
        "full_name": "Hemophagocytic Lymphohistiocytosis / Macrophage Activation Syndrome",
        "mechanism": (
            "Severe immune hyperactivation with uncontrolled macrophage/histiocyte "
            "proliferation and hemophagocytosis; overlaps with severe CRS"
        ),
        "grading": {
            "diagnostic_criteria": "Ferritin >10,000, cytopenias, hepatosplenomegaly, "
            "coagulopathy, elevated sIL-2R, elevated triglycerides",
        },
        "incidence": "1-5% of CAR-T recipients",
        "timing": "Typically day 5-20 post-infusion, during or after CRS",
        "management": [
            "Etoposide — severe/refractory cases",
            "Anakinra (IL-1R antagonist)",
            "Ruxolitinib (JAK inhibitor)",
            "High-dose corticosteroids",
            "Tocilizumab (if concurrent CRS)",
        ],
        "biomarkers": ["Ferritin (>10,000)", "Fibrinogen (decreased)", "Triglycerides", "sIL-2R"],
        "risk_factors": ["Severe CRS", "High tumor burden", "Prior HLH history"],
    },
    "CYTOPENIAS": {
        "full_name": "Prolonged Cytopenias",
        "mechanism": (
            "Bone marrow suppression from lymphodepletion, CRS-related inflammation, "
            "and direct CAR-T effects on hematopoietic progenitors"
        ),
        "grading": {
            "early": "Within 30 days (expected from lymphodepletion)",
            "prolonged": "Beyond 30 days (concerning)",
            "severe": "Grade 3-4 neutropenia, thrombocytopenia, or anemia >90 days",
        },
        "incidence": "30-50% prolonged cytopenias (>30 days)",
        "timing": "Onset day 0-7; prolonged in 30-50%",
        "management": [
            "G-CSF for severe neutropenia (>14 days post-infusion)",
            "Platelet/RBC transfusion support",
            "TPO-RA for prolonged thrombocytopenia",
            "Infection monitoring and prophylaxis",
        ],
        "biomarkers": ["CBC with differential", "Reticulocyte count", "Bone marrow biopsy if prolonged"],
        "risk_factors": [
            "Severe CRS", "High Flu/Cy dose", "Pre-existing cytopenias",
            "Prior transplant", "Heavy prior therapy",
        ],
    },
    "TLS": {
        "full_name": "Tumor Lysis Syndrome",
        "mechanism": "Rapid tumor cell destruction releasing intracellular contents",
        "grading": {
            "laboratory": "2+ of: uric acid, K+, phosphate elevated; Ca2+ decreased",
            "clinical": "Laboratory TLS + renal/cardiac/neurologic sequelae",
        },
        "incidence": "Rare with CAR-T (<5%); more common in high tumor burden",
        "timing": "Day 1-7 post-infusion",
        "management": [
            "Rasburicase for hyperuricemia",
            "Allopurinol prophylaxis",
            "Aggressive hydration",
            "Electrolyte monitoring and correction",
        ],
        "biomarkers": ["Uric acid", "Potassium", "Phosphate", "Calcium", "LDH", "Creatinine"],
        "risk_factors": ["High tumor burden", "ALL", "Bulky disease"],
    },
    "GVHD": {
        "full_name": "Graft-versus-Host Disease",
        "mechanism": (
            "Allogeneic CAR-T cells with intact TCR attacking host tissues; "
            "relevant primarily for donor-derived or allogeneic CAR-T products"
        ),
        "grading": {
            "acute": "Skin, GI, liver involvement (Glucksberg criteria)",
            "chronic": "Skin, mouth, liver, lung, joints",
        },
        "incidence": "Rare with autologous CAR-T; risk with allogeneic products",
        "timing": "Variable; typically >2 weeks post-infusion",
        "management": [
            "TCR knockout in allogeneic products (prevention)",
            "Corticosteroids",
            "Ruxolitinib for steroid-refractory",
            "Ibrutinib for chronic GVHD",
        ],
        "biomarkers": ["Skin biopsy", "Liver function tests", "GI biopsy"],
        "risk_factors": [
            "Allogeneic (donor-derived) CAR-T",
            "Post-transplant donor lymphocyte-derived CAR-T",
            "Intact TCR on product",
        ],
    },
    "ON_TARGET_OFF_TUMOR": {
        "full_name": "On-Target, Off-Tumor Toxicity",
        "mechanism": (
            "CAR-T cells attacking normal tissues that express the target antigen; "
            "severity depends on target expression pattern in normal tissues"
        ),
        "grading": {
            "low_risk": "B-cell aplasia with CD19 (manageable with IVIG)",
            "moderate_risk": "Skin/nail toxicity with GPRC5D; mucosal with Claudin18.2",
            "high_risk": "Cardiac/pulmonary with HER2; myelosuppression with CD33",
        },
        "examples": [
            "CD19 → B-cell aplasia (manageable)",
            "BCMA → plasma cell depletion (manageable)",
            "HER2 → cardiac toxicity (fatal case reported 2010)",
            "CD33 → myelosuppression",
            "GPRC5D → skin/nail/taste changes",
            "Claudin18.2 → GI mucosal toxicity",
        ],
        "management": [
            "Target selection: tumor-restricted antigens preferred",
            "Affinity tuning: lower-affinity CARs spare low-expression normal tissue",
            "Logic gates: AND-gate CARs requiring two antigens",
            "Safety switches: iCasp9, EGFRt for rapid elimination",
        ],
        "biomarkers": ["Target-antigen dependent"],
        "risk_factors": ["Broad normal tissue expression", "High-affinity scFv"],
    },
    "COAGULOPATHY": {
        "full_name": "Disseminated Intravascular Coagulation / Coagulopathy",
        "mechanism": (
            "CRS-driven endothelial activation releases tissue factor and von Willebrand "
            "factor, triggering coagulation cascade consumption. IL-6 and TNF-alpha promote "
            "a procoagulant state."
        ),
        "grading": {
            "grade_1": "D-dimer elevated <4x ULN, fibrinogen normal",
            "grade_2": "Fibrinogen 100-200 mg/dL, D-dimer 4-10x ULN",
            "grade_3": "Fibrinogen <100 mg/dL, PT/INR elevated, active bleeding",
            "grade_4": "Life-threatening hemorrhage, multi-organ DIC, hemodynamic instability",
        },
        "grading_system": "ISTH DIC scoring system + CTCAE v5.0",
        "incidence": "5-15% with severe CRS; up to 30% subclinical coagulopathy",
        "timing": "Onset day 3-7 post-infusion, concurrent with peak CRS",
        "management": [
            "Monitor fibrinogen q6h during CRS",
            "Cryoprecipitate if fibrinogen <150 mg/dL",
            "Platelet transfusion if <50,000/\u03bcL with bleeding",
            "Heparin generally contraindicated",
            "Treat underlying CRS aggressively",
        ],
        "biomarkers": ["fibrinogen", "d_dimer", "pt_inr", "platelets", "thrombin_time"],
        "risk_factors": [
            "High tumor burden", "Grade 3+ CRS",
            "Baseline low fibrinogen", "Prior anticoagulation",
        ],
    },
    "CARDIAC_TOXICITY": {
        "full_name": "Cardiac Toxicity",
        "mechanism": (
            "Direct myocardial injury from cytokine storm (IL-6, TNF-alpha), catecholamine "
            "surge, and capillary leak. Can progress to cardiomyopathy, arrhythmia, or "
            "cardiac arrest."
        ),
        "grading": {
            "grade_1": "Troponin elevated above 99th percentile, asymptomatic",
            "grade_2": "Troponin elevation with ECG changes (ST-T wave), mild LV dysfunction (LVEF 40-50%)",
            "grade_3": "Symptomatic heart failure, LVEF <40%, hemodynamically significant arrhythmia",
            "grade_4": "Cardiogenic shock, life-threatening arrhythmia, cardiac arrest",
        },
        "grading_system": "CTCAE v5.0 + ACC/AHA cardio-oncology guidelines",
        "incidence": "5-20% troponin elevation; 2-5% clinically significant cardiac events",
        "timing": "Onset day 2-10, often concurrent with CRS peak",
        "management": [
            "Baseline echocardiogram and troponin",
            "Serial troponin and BNP monitoring during CRS",
            "Cardiology consultation if troponin >3x ULN or BNP >500",
            "Standard heart failure management if LVEF depressed",
            "ICU monitoring for hemodynamic instability",
        ],
        "biomarkers": ["troponin_hstni", "bnp_nt_probnp", "ecg", "echocardiogram_lvef"],
        "risk_factors": [
            "Pre-existing cardiac disease", "Grade 3+ CRS",
            "Age >65", "Prior anthracycline exposure", "Hypertension",
        ],
    },
    "RENAL_TOXICITY": {
        "full_name": "Renal Toxicity / Acute Kidney Injury",
        "mechanism": (
            "TLS-associated uric acid nephropathy, cytokine-mediated capillary leak with "
            "prerenal AKI, or direct tubular injury from inflammatory mediators."
        ),
        "grading": {
            "grade_1": "Creatinine 1.5-2x baseline, eGFR mildly reduced",
            "grade_2": "Creatinine 2-3x baseline, oliguria <0.5 mL/kg/h for 6-12h",
            "grade_3": "Creatinine >3x baseline or >4.0 mg/dL, dialysis indicated",
            "grade_4": "Life-threatening; dialysis-dependent, multi-organ failure",
        },
        "grading_system": "KDIGO AKI staging + CTCAE v5.0",
        "incidence": "10-25% any grade AKI; 3-8% dialysis-requiring",
        "timing": "Onset day 2-14, often with CRS or TLS",
        "management": [
            "Aggressive IV hydration during CRS",
            "Monitor electrolytes q6-12h (K+, phosphate, uric acid)",
            "Rasburicase prophylaxis for high tumor burden",
            "Avoid nephrotoxins (NSAIDs, contrast, aminoglycosides)",
            "Early nephrology consultation if creatinine >2x baseline",
        ],
        "biomarkers": ["creatinine", "bun", "potassium", "phosphate", "uric_acid", "urine_output"],
        "risk_factors": [
            "High tumor burden (TLS risk)", "Pre-existing CKD",
            "Severe CRS", "Concomitant nephrotoxic medications",
        ],
    },
    "SECONDARY_MALIGNANCY": {
        "full_name": "Secondary Malignancy (T-Cell Lymphoma)",
        "mechanism": (
            "Insertional mutagenesis from viral vector integration near proto-oncogenes, "
            "or clonal expansion of transduced T-cells. FDA class-wide boxed warning added "
            "November 2023."
        ),
        "grading": {
            "grade_1": "Clonal expansion detected on integration site analysis, no clinical lymphoma",
            "grade_2": "Atypical T-cell population detected, monitoring required",
            "grade_3": "Confirmed T-cell lymphoma diagnosis, treatment-requiring",
            "grade_4": "Aggressive/refractory T-cell lymphoma, life-threatening",
        },
        "grading_system": "FDA Risk Evaluation + WHO lymphoma classification",
        "incidence": "~33 reported cases by Nov 2023 across all approved products; estimated 0.1-0.4% incidence",
        "timing": "Onset 2 months to 5+ years post-infusion (median ~1-2 years)",
        "management": [
            "FDA requires 15-year long-term follow-up for all CAR-T recipients",
            "Annual CBC with differential and flow cytometry",
            "Integration site analysis if unexplained lymphocytosis",
            "Standard lymphoma workup if T-cell malignancy suspected",
            "Report to REMS program and FDA MedWatch",
        ],
        "biomarkers": ["cbc_differential", "flow_cytometry_tcrVb", "integration_site_analysis", "car_transgene_pcr"],
        "risk_factors": [
            "Retroviral vectors (higher insertional risk than lentiviral)",
            "Multiple prior lines of therapy",
            "Prior DNA-damaging agents", "Baseline clonal hematopoiesis",
        ],
    },
    "IMMUNE_RECONSTITUTION": {
        "full_name": "Immune Reconstitution — Long-term Immune Recovery",
        "mechanism": (
            "Prolonged B-cell aplasia and hypogammaglobulinemia following CAR-T therapy, "
            "with delayed recovery of humoral and cellular immunity"
        ),
        "grades": {
            "normal": "Full immune reconstitution by 12 months",
            "delayed": "Persistent hypogammaglobulinemia >12 months",
            "severe": "Recurrent opportunistic infections requiring prophylaxis >24 months",
        },
        "incidence": "Delayed reconstitution in 30-50% at 12 months; 10-20% requiring IVIG >24 months",
        "management": [
            "Monthly IVIG (trough >400 mg/dL)",
            "Antimicrobial prophylaxis (acyclovir, fluconazole, TMP-SMX)",
            "Vaccination schedule (inactivated only, start at 6-12 months post CAR-T)",
            "Serial immunoglobulin level monitoring",
        ],
        "biomarkers": ["IgG levels", "CD4 count", "CD19+ B-cell recovery", "Vaccine titers"],
        "risk_factors": [
            "Prior lines of therapy (>3)",
            "CD19 CAR-T (vs BCMA)",
            "Flu/Cy conditioning intensity",
            "Age >60",
        ],
    },
    "PROLONGED_CYTOPENIAS": {
        "full_name": "Prolonged Cytopenias — Beyond 30 Days",
        "mechanism": (
            "Persistent bone marrow suppression beyond expected lymphodepletion recovery, "
            "driven by inflammatory marrow damage, stromal injury, and impaired hematopoietic recovery"
        ),
        "grades": {
            "grade1": "ANC 1000-1500 or platelets 50-75K at day 30",
            "grade2": "ANC 500-1000 or platelets 25-50K at day 30",
            "grade3": "ANC <500 or platelets <25K persisting >30 days",
            "grade4": "Transfusion-dependent or ANC <100 >60 days",
        },
        "incidence": "Any grade: 30-60%; Grade 3+: 15-30% at day 30; Persistent >90 days: 5-15%",
        "management": [
            "G-CSF (if ANC <500 after day 14)",
            "Thrombopoietin receptor agonists (eltrombopag, romiplostim)",
            "Platelet and RBC transfusion support",
            "Bone marrow biopsy if persistent >90 days",
            "Iron studies and B12/folate supplementation",
        ],
        "biomarkers": ["Baseline ANC/platelets", "Ferritin (HLH screen)", "Reticulocyte count", "EPO level", "Bone marrow cellularity"],
        "risk_factors": [
            "High tumor burden",
            "Prior HSCT",
            "CRS grade >=2",
            ">4 prior lines",
            "Baseline cytopenias",
        ],
    },
    "INFECTIONS": {
        "full_name": "Infections — Post CAR-T",
        "mechanism": (
            "Immunocompromised state from lymphodepletion, B-cell aplasia, hypogammaglobulinemia, "
            "and prolonged cytopenias predisposing to bacterial, viral, and fungal infections"
        ),
        "grades": {
            "mild": "Localized infection, oral antibiotics",
            "moderate": "IV antibiotics, no ICU",
            "severe": "ICU admission, sepsis",
            "fatal": "Treatment-related mortality from infection",
        },
        "incidence": "Any infection: 20-40% in first 90 days; Severe: 10-20%; Fatal: 2-5%",
        "management": [
            "Prophylaxis: acyclovir/valacyclovir (HSV/VZV), fluconazole (fungal), TMP-SMX (PJP)",
            "PJP prophylaxis for 6-12 months post CAR-T",
            "IVIG replacement if IgG <400",
            "COVID-19 pre-exposure prophylaxis (tixagevimab/cilgavimab or equivalent)",
            "Avoid live vaccines for 12+ months",
        ],
        "biomarkers": ["ANC", "IgG trough", "CD4 count", "CMV/EBV PCR monitoring", "Beta-D-glucan", "Galactomannan"],
        "risk_factors": [
            "Prolonged cytopenias",
            "CRS grade >=2",
            "High-dose corticosteroids",
            "Prior HSCT",
            "Hypogammaglobulinemia",
        ],
    },
    "HEPATOTOXICITY": {
        "full_name": "Hepatotoxicity — Liver Injury",
        "mechanism": (
            "Immune-mediated liver injury from CAR-T cell infiltration, cytokine-driven hepatocyte "
            "damage, or on-target off-tumor toxicity with liver-antigen-directed CARs (GPC3, Claudin18.2)"
        ),
        "grades": {
            "grade1": "AST/ALT 1-3x ULN",
            "grade2": "AST/ALT 3-5x ULN",
            "grade3": "AST/ALT 5-20x ULN",
            "grade4": "AST/ALT >20x ULN or liver failure",
        },
        "incidence": "Any grade: 10-25%; Grade 3+: 3-8%; More common with GPC3/Claudin18.2-targeted therapies",
        "management": [
            "Hold therapy if grade 3+",
            "N-acetylcysteine if severe",
            "Corticosteroids if immune-mediated",
            "Hepatology consult if grade >=3",
            "Serial LFT monitoring (daily during CRS)",
        ],
        "biomarkers": ["AST/ALT", "Bilirubin", "INR", "Albumin", "Ferritin (HLH screen)"],
        "risk_factors": [
            "GPC3-targeted CAR-T (liver antigen)",
            "Hepatic tumor burden",
            "Prior hepatotoxic therapy",
            "Viral hepatitis (HBV/HCV)",
        ],
    },
    "MACROPHAGE_ACTIVATION": {
        "full_name": "Macrophage Activation Syndrome — Distinct from CRS",
        "mechanism": (
            "Pathologic macrophage hyperactivation with hemophagocytosis, distinct from but overlapping "
            "with severe CRS. Driven by IFN-gamma, IL-18, and CXCL9 from CAR-T cells activating "
            "tissue-resident macrophages"
        ),
        "grades": {
            "mild": "Ferritin 5000-10000, mild organ dysfunction",
            "moderate": "Ferritin >10000, multi-organ involvement",
            "severe": "Ferritin >50000, DIC, liver failure, renal failure",
        },
        "incidence": "Subclinical: 10-20%; Clinical MAS: 1-5%; Overlaps with severe CRS",
        "management": [
            "Anakinra (IL-1 blockade, preferred first-line for MAS)",
            "Ruxolitinib (JAK1/2 inhibition)",
            "Emapalumab (anti-IFN-gamma, for refractory)",
            "Etoposide (refractory HLH/MAS only)",
            "Avoid tocilizumab monotherapy (may worsen hepatic MAS)",
        ],
        "biomarkers": [
            "Ferritin (>10,000 highly suggestive)",
            "sIL-2R/sCD25",
            "Fibrinogen (falling)",
            "Triglycerides (rising)",
            "NK cell function",
            "Bone marrow hemophagocytosis",
        ],
        "risk_factors": [
            "High tumor burden",
            "CRS grade >=3",
            "Pre-existing immune dysregulation",
            "EBV reactivation",
        ],
    },
}


# =============================================================================
# 3. CART_MANUFACTURING — Manufacturing process knowledge (~15 entries)
# =============================================================================

CART_MANUFACTURING: Dict[str, Dict[str, Any]] = {
    "lentiviral_transduction": {
        "description": "Lentiviral vector transduction of T-cells to introduce CAR transgene",
        "typical_efficiency": "20-60%",
        "target_vcn": "0.5-5.0 copies/cell",
        "critical_parameters": [
            "Multiplicity of infection (MOI)",
            "Viral titer (TU/mL)",
            "Polybrene/RetroNectin concentration",
            "Transduction duration (16-24 hours)",
            "Spinoculation protocol",
        ],
        "failure_modes": [
            "Low viral titer", "Insertional mutagenesis risk",
            "RCL (replication-competent lentivirus) contamination",
            "Low transduction efficiency", "High VCN (safety concern)",
        ],
        "release_criteria": {
            "CAR_expression": ">10-20% CAR+ by flow cytometry",
            "VCN": "<5 copies/cell (FDA guideline)",
            "RCL": "Negative",
        },
        "products_using": ["Kymriah", "Breyanzi", "Abecma", "Carvykti"],
    },
    "retroviral_transduction": {
        "description": "Gamma-retroviral vector transduction (requires actively dividing cells)",
        "typical_efficiency": "40-80%",
        "target_vcn": "1-5 copies/cell",
        "critical_parameters": [
            "Cell activation status (must be actively dividing)",
            "MOI", "Viral titer", "RetroNectin coating",
        ],
        "failure_modes": [
            "RCR (replication-competent retrovirus) risk",
            "Requires cell division for integration",
            "Insertional mutagenesis (historical concern with MLV)",
        ],
        "release_criteria": {
            "CAR_expression": ">40% CAR+",
            "RCR": "Negative",
        },
        "products_using": ["Yescarta", "Tecartus"],
    },
    "t_cell_activation": {
        "description": "T-cell activation using anti-CD3/CD28 stimulation prior to transduction",
        "methods": [
            "Anti-CD3/CD28 magnetic beads (Dynabeads)",
            "TransAct (soluble anti-CD3/CD28)",
            "Plate-bound anti-CD3 + soluble anti-CD28",
            "Artificial APCs",
        ],
        "critical_parameters": [
            "Bead-to-cell ratio (3:1 typical for Dynabeads)",
            "Activation duration (24-72 hours)",
            "Cytokine cocktail (IL-2, IL-7, IL-15)",
            "Activation marker verification (CD25, CD69)",
        ],
        "failure_modes": [
            "Over-activation leading to exhaustion",
            "Poor activation in lymphopenic patients",
            "Contaminating tumor cells activated",
        ],
    },
    "ex_vivo_expansion": {
        "description": "Expansion of CAR-T cells to therapeutic dose (typically 1e8-1e9 cells)",
        "target_dose": "1e6 to 1e8 CAR+ cells/kg",
        "expansion_fold": "100-1000x typical",
        "duration": "7-14 days",
        "platforms": [
            "G-Rex (gas-permeable static culture)",
            "WAVE/Xuri bioreactor (rocking motion)",
            "CliniMACS Prodigy (automated closed system)",
            "Lonza Cocoon (automated)",
            "Standard culture flasks/bags",
        ],
        "critical_parameters": [
            "Cell density", "Media/feed schedule",
            "Cytokine concentration (IL-2: 50-300 IU/mL; IL-7/IL-15: 5-10 ng/mL)",
            "pH and dissolved O2", "Temperature (37C)",
        ],
        "failure_modes": [
            "Insufficient expansion (<50x)", "T-cell exhaustion from over-expansion",
            "Contamination", "Shift to effector phenotype (loss of Tscm/Tcm)",
        ],
    },
    "leukapheresis": {
        "description": "Collection of patient PBMCs via apheresis as starting material",
        "volume_processed": "2-3 blood volumes typical",
        "target_collection": ">2e9 total nucleated cells",
        "critical_parameters": [
            "Timing relative to prior therapy (washout period)",
            "Anti-coagulant (ACD-A)",
            "Flow rate and processing time",
            "T-cell purity (CD3+ percentage)",
        ],
        "failure_modes": [
            "Low T-cell count (lymphopenic patient)",
            "High circulating tumor cell contamination",
            "Poor T-cell quality (exhausted/senescent)",
            "Apheresis access issues",
        ],
        "enrichment_options": [
            "CD4/CD8 positive selection (CliniMACS)",
            "T-cell enrichment by elutriation",
            "Tumor cell depletion",
        ],
    },
    "cryopreservation": {
        "description": "Cryopreservation of final CAR-T product for storage and shipping",
        "media": "CryoStor CS10 (10% DMSO) or equivalent",
        "protocol": [
            "Controlled-rate freezing (-1C/min to -40C, -10C/min to -80C)",
            "Transfer to liquid nitrogen vapor phase (-150 to -196C)",
            "Chain of custody and temperature monitoring",
        ],
        "critical_parameters": [
            "Post-thaw viability (>70% required)",
            "Post-thaw recovery (% of pre-freeze)",
            "Post-thaw potency (functional assay)",
            "DMSO concentration",
            "Cooling rate",
        ],
        "failure_modes": [
            "Freezing too fast/slow", "Temperature excursion during shipping",
            "Thaw protocol deviation", "Low post-thaw viability",
        ],
    },
    "release_testing": {
        "description": "Quality control testing of final product before release",
        "required_tests": {
            "identity": "CAR expression by flow cytometry (Protein L or anti-idiotype)",
            "purity": "CD3+ percentage; tumor cell absence",
            "viability": ">70% viable cells (trypan blue or flow)",
            "potency": "Cytotoxicity assay or IFN-gamma secretion",
            "sterility": "USP <71> or BacT/ALERT (14-day or rapid)",
            "endotoxin": "<3.5 EU/mL (LAL or recombinant Factor C)",
            "mycoplasma": "PCR or culture-based",
            "VCN": "<5 copies/cell (qPCR for integrated transgene)",
            "RCL_RCR": "Negative for replication-competent virus",
        },
        "turnaround": "Sterility is rate-limiting (14-day incubation)",
        "failure_modes": [
            "Sterility failure", "Low viability (<70%)",
            "Low CAR expression", "Potency failure",
            "Out-of-spec VCN",
        ],
    },
    "point_of_care_manufacturing": {
        "description": "Decentralized manufacturing at or near the treatment site",
        "platforms": [
            "CliniMACS Prodigy (Miltenyi)",
            "Lonza Cocoon",
            "Ori Biotech",
        ],
        "advantages": [
            "Reduced vein-to-vein time (3-5 days vs 3-6 weeks)",
            "Eliminated shipping/cold-chain logistics",
            "Fresh (non-cryopreserved) product possible",
            "Lower manufacturing cost potential",
        ],
        "challenges": [
            "GMP compliance at multiple sites",
            "Operator training and standardization",
            "Regulatory framework for decentralized manufacturing",
            "Scale-out (not scale-up) model",
        ],
    },
    "lymphodepletion": {
        "description": "Pre-infusion conditioning chemotherapy to create immunologic space",
        "standard_regimen": "Fludarabine 30mg/m2 x3 days + Cyclophosphamide 500mg/m2 x3 days",
        "alternative_regimens": [
            "Flu 25/Cy 250 (reduced intensity)",
            "Bendamustine (CD19 products)",
            "Flu/Cy/Alemtuzumab (allogeneic products)",
        ],
        "mechanism": [
            "Depletion of regulatory T-cells",
            "Reduction of competing lymphocytes",
            "Induction of homeostatic cytokines (IL-7, IL-15)",
            "Creation of immunologic space for CAR-T expansion",
        ],
        "critical_parameters": [
            "Timing: typically day -5 to day -3 before infusion",
            "Dose adjustments for renal function",
            "ALC (absolute lymphocyte count) at infusion",
        ],
        "failure_modes": [
            "Inadequate lymphodepletion (poor CAR-T expansion)",
            "Excessive myelosuppression",
            "Infection before CAR-T infusion",
        ],
    },
    "vein_to_vein_time": {
        "description": "Total time from leukapheresis to patient infusion",
        "typical_centralized": "3-6 weeks (commercial products)",
        "typical_academic": "2-4 weeks",
        "typical_poc": "3-7 days (point-of-care)",
        "bottlenecks": [
            "Leukapheresis scheduling",
            "Shipping to manufacturing site",
            "Manufacturing slot availability",
            "Release testing (14-day sterility)",
            "Return shipping",
            "Lymphodepletion scheduling",
        ],
        "clinical_impact": [
            "Disease progression during manufacturing",
            "Need for bridging therapy",
            "Patient dropout/death before infusion",
            "30-40% of patients never receive manufactured product",
        ],
    },
    "non_viral_transposon": {
        "description": "Non-viral gene transfer using Sleeping Beauty or PiggyBac transposon systems. DNA plasmid + transposase electroporated into T-cells.",
        "typical_efficiency": "5-15% stable integration",
        "critical_parameters": [
            "Transposon:transposase ratio (10:1 optimal)",
            "Electroporation voltage/pulse (1700V, 20ms)",
            "Plasmid quality (endotoxin-free)",
            "Cell density at electroporation (1e7/mL)",
        ],
        "failure_modes": [
            "Low integration efficiency vs viral",
            "Multi-copy insertions",
            "Transposase re-mobilization risk",
            "Cell death from electroporation (30-50% loss)",
        ],
        "regulatory_note": "Lower manufacturing cost ($20-50K vs $100-300K viral). PACT alliance, Ziopharm programs.",
    },
    "mrna_electroporation": {
        "description": "Transient CAR expression via in vitro transcribed mRNA electroporation. No genomic integration — CAR expression lasts 24-72 hours, enabling repeated dosing.",
        "typical_efficiency": "70-95% transfection, expression decays over 48-72h",
        "critical_parameters": [
            "mRNA cap structure (CleanCap/ARCA)",
            "Poly-A tail length (>100 nt)",
            "Modified nucleosides (pseudouridine, m1\u03a8)",
            "Multiple dosing schedule (weekly x 3-6)",
        ],
        "failure_modes": [
            "Rapid CAR loss requiring repeated infusions",
            "Anti-RNA immune responses",
            "Inconsistent expression kinetics",
            "Manufacturing complexity for multiple doses",
        ],
        "regulatory_note": "Enables iterative dose-finding with safety shutoff. UPenn mesoCAR-T program.",
    },
    "crispr_knock_in": {
        "description": "CRISPR-Cas9 site-specific integration of CAR construct at the TRAC locus, disrupting endogenous TCR while achieving uniform CAR expression under the T-cell receptor promoter.",
        "typical_efficiency": "20-50% knock-in at TRAC locus",
        "critical_parameters": [
            "Guide RNA specificity (off-target <0.1%)",
            "Cas9 delivery (RNP preferred over plasmid)",
            "HDR template design (AAV6 or dsDNA)",
            "Selection/enrichment post-editing",
        ],
        "failure_modes": [
            "Off-target editing (chromosomal translocations)",
            "Low HDR efficiency in primary T-cells",
            "p53-mediated cell death from DSBs",
            "Large-scale Cas9 RNP manufacturing",
        ],
        "regulatory_note": "TRAC-integrated CAR shows superior persistence (Eyquem et al. Nature 2017). Caribou Biosciences CB-010 program.",
    },
    "ipsc_derived": {
        "description": "CAR-T cells derived from induced pluripotent stem cells (iPSCs) — master cell bank enables unlimited, standardized manufacturing of off-the-shelf allogeneic product.",
        "typical_efficiency": "Differentiation to mature T-cells: 40-60% yield",
        "critical_parameters": [
            "iPSC quality (karyotype stability)",
            "Directed differentiation protocol (21-28 days)",
            "TCR/HLA knockout (TRAC/B2M deletion)",
            "NK-resistance engineering (HLA-E/CD47 knock-in)",
        ],
        "failure_modes": [
            "Incomplete T-cell maturation (immature phenotype)",
            "Genomic instability during culture",
            "Residual undifferentiated cells (teratoma risk)",
            "Functional inferiority vs autologous CAR-T",
        ],
        "regulatory_note": "Fate Therapeutics FT819 (Phase 1). Standardized manufacturing eliminates patient-specific variability.",
    },
    "automated_manufacturing": {
        "description": "Automated closed-system manufacturing using platforms like CliniMACS Prodigy (Miltenyi), Cocoon (Lonza), or Sepax (Cytiva). Reduces cleanroom requirements, operator variability, and contamination risk.",
        "typical_efficiency": "Comparable to manual (85-95% viability, 50-200x expansion)",
        "critical_parameters": [
            "Protocol programming and validation",
            "Reagent lot qualification",
            "Sensor calibration (temperature, pH, dissolved O2)",
            "Closed-system integrity testing",
        ],
        "failure_modes": [
            "Software/hardware failures during automated run",
            "Limited flexibility for non-standard protocols",
            "High capital equipment cost ($200-500K per unit)",
            "Batch size limitations",
        ],
        "regulatory_note": "Enables decentralized manufacturing at clinical sites. Reduced vein-to-vein time to 7-12 days. Novartis T-Charge technology.",
    },
    "automated_closed_system": {
        "description": "Fully automated closed manufacturing platforms that perform all steps from T-cell enrichment through harvest in a single device",
        "key_parameters": {
            "platforms": "CliniMACS Prodigy (Miltenyi), Cocoon (Lonza), Sepax (Cytiva)",
            "processing_time": "7-12 days",
            "operator_interventions": "2-3 per run",
            "contamination_risk": "reduced 10-fold vs open process",
            "scalability": "1 device = 1 patient",
        },
        "clinical_relevance": "Enables decentralized/point-of-care manufacturing; reduces vein-to-vein time from 4-6 weeks to 7-14 days; critical for aggressive diseases where patients deteriorate during wait",
    },
    "t_cell_selection": {
        "description": "Defined composition T-cell products using CD4+ and CD8+ selection to achieve optimal 1:1 ratio",
        "key_parameters": {
            "target_ratio": "1:1 CD4:CD8",
            "selection_method": "CliniMACS CD4/CD8 selection or FACS sorting",
            "yield": "60-80% recovery",
            "purity": ">95% CD3+",
            "products_using": "Breyanzi (lisocabtagene maraleucel)",
        },
        "clinical_relevance": "Defined composition reduces inter-patient variability; Breyanzi's defined ratio showed consistent expansion kinetics; enables better dose-response prediction",
    },
    "viral_vector_production": {
        "description": "Large-scale production of lentiviral or retroviral vectors for CAR-T transduction",
        "key_parameters": {
            "cell_line": "HEK293T (transient) or stable producer lines",
            "scale": "10L-200L bioreactors",
            "titer": "1e7-1e9 TU/mL",
            "lot_size": "sufficient for 50-500 patients",
            "cost": "$50,000-$500,000 per lot",
            "timeline": "8-12 weeks for GMP lot",
            "testing": "RCL, sterility, endotoxin, identity, potency, adventitious agents",
        },
        "clinical_relevance": "Vector supply is the primary bottleneck for CAR-T manufacturing scale-up; stable producer lines reduce cost 5-10x but require 12-18 months to develop",
    },
    "quality_control_analytics": {
        "description": "Comprehensive analytical testing required for CAR-T product release under GMP",
        "key_parameters": {
            "identity": "CAR expression by flow cytometry (% CAR+)",
            "purity": "CD3+ %, viability by trypan blue or Nucleocounter",
            "potency": "Cytotoxicity assay (E:T = 1:1 to 10:1) or cytokine secretion (IFN-gamma ELISPOT)",
            "sterility": "14-day direct inoculation (USP <71>)",
            "mycoplasma": "qPCR (result in <4 hours)",
            "endotoxin": "LAL (<5 EU/kg)",
            "RCL_RCR": "Replication competent lentivirus/retrovirus (qPCR + culture)",
            "VCN": "Vector copy number by ddPCR (target 0.5-5 copies/cell)",
            "in_process": "Cell count, viability, transduction efficiency at days 0, 3, 7, harvest",
        },
        "clinical_relevance": "Release testing takes 14+ days (sterility rate-limiting); rapid sterility methods (BacT/ALERT, Bactec) enable earlier release; VCN monitoring critical for insertional mutagenesis risk (FDA secondary malignancy concern)",
    },
    "supply_chain_logistics": {
        "description": "End-to-end logistics from leukapheresis collection through product delivery and infusion",
        "key_parameters": {
            "leukapheresis_shipping": "Fresh (2-8C) within 24-48 hours or cryopreserved (-150C)",
            "manufacturing_site": "Centralized (Novartis Morris Plains, Kite El Segundo, BMS Bothell) or decentralized",
            "product_shipping": "Cryopreserved in vapor-phase LN2 shipper (-150C)",
            "chain_of_identity": "Unique patient identifier, barcode tracking, chain of custody forms",
            "thaw_to_infuse": "Within 30 minutes of thaw",
            "scheduling_coordination": "Lymphodepletion timing, bed availability, ICU capacity, tocilizumab stock",
        },
        "clinical_relevance": "30-40% of patients referred for CAR-T never receive product due to disease progression, manufacturing failure, or logistic delays; supply chain optimization is as critical as CAR design",
    },
}


# =============================================================================
# 4. PUBLIC API — Context retrieval functions
# =============================================================================

def get_target_context(antigen: str) -> str:
    """Get formatted knowledge context for a target antigen.

    Args:
        antigen: Target antigen name (case-insensitive, e.g. 'CD19', 'BCMA').

    Returns:
        Formatted string with target knowledge, or empty string if not found.
    """
    key = antigen.upper()
    if key not in CART_TARGETS:
        # Try case-insensitive match
        for k in CART_TARGETS:
            if k.upper() == key:
                key = k
                break
        else:
            return ""

    data = CART_TARGETS[key]
    lines = [f"## Target Antigen: {key}"]
    lines.append(f"- **Protein:** {data['protein']}")
    lines.append(f"- **Expression:** {data['expression']}")
    lines.append(f"- **Diseases:** {', '.join(data['diseases'])}")

    if data.get("approved_products"):
        lines.append(f"- **Approved Products:** {'; '.join(data['approved_products'])}")
    if data.get("key_trials"):
        lines.append(f"- **Key Trials:** {', '.join(data['key_trials'])}")
    if data.get("known_resistance"):
        lines.append(f"- **Resistance Mechanisms:** {'; '.join(data['known_resistance'])}")
    if data.get("toxicity_profile"):
        tox = ", ".join(f"{k}: {v}" for k, v in data["toxicity_profile"].items())
        lines.append(f"- **Toxicity Profile:** {tox}")
    if data.get("normal_tissue"):
        lines.append(f"- **Normal Tissue Expression:** {data['normal_tissue']}")

    return "\n".join(lines)


def get_toxicity_context(toxicity: str) -> str:
    """Get formatted knowledge context for a toxicity profile.

    Args:
        toxicity: Toxicity identifier (e.g. 'CRS', 'ICANS', 'HLH_MAS').

    Returns:
        Formatted string with toxicity knowledge, or empty string if not found.
    """
    key = toxicity.upper()
    if key not in CART_TOXICITIES:
        for k in CART_TOXICITIES:
            if k.upper() == key:
                key = k
                break
        else:
            return ""

    data = CART_TOXICITIES[key]
    lines = [f"## Toxicity: {data['full_name']} ({key})"]
    lines.append(f"- **Mechanism:** {data['mechanism']}")
    lines.append(f"- **Incidence:** {data.get('incidence', 'Variable')}")
    lines.append(f"- **Timing:** {data.get('timing', 'Variable')}")

    if data.get("management"):
        lines.append("- **Management:**")
        for m in data["management"]:
            lines.append(f"  - {m}")
    if data.get("biomarkers"):
        lines.append(f"- **Biomarkers:** {', '.join(data['biomarkers'])}")
    if data.get("risk_factors"):
        lines.append(f"- **Risk Factors:** {'; '.join(data['risk_factors'])}")

    return "\n".join(lines)


def get_manufacturing_context(process: str) -> str:
    """Get formatted knowledge context for a manufacturing process.

    Args:
        process: Process identifier (e.g. 'lentiviral_transduction', 'expansion').

    Returns:
        Formatted string with manufacturing knowledge, or empty string.
    """
    key = process.lower().replace(" ", "_")
    if key not in CART_MANUFACTURING:
        for k in CART_MANUFACTURING:
            if key in k.lower():
                key = k
                break
        else:
            return ""

    data = CART_MANUFACTURING[key]
    lines = [f"## Manufacturing: {key.replace('_', ' ').title()}"]
    lines.append(f"- **Description:** {data['description']}")

    for field in ["typical_efficiency", "target_vcn", "target_dose",
                  "expansion_fold", "duration", "standard_regimen"]:
        if field in data:
            label = field.replace("_", " ").title()
            lines.append(f"- **{label}:** {data[field]}")

    if data.get("critical_parameters"):
        lines.append("- **Critical Parameters:**")
        for p in data["critical_parameters"]:
            lines.append(f"  - {p}")
    if data.get("failure_modes"):
        lines.append("- **Failure Modes:**")
        for f in data["failure_modes"]:
            lines.append(f"  - {f}")

    return "\n".join(lines)


def get_all_context_for_query(query: str) -> str:
    """Extract all relevant knowledge context from a query string.

    Scans the query for mentions of target antigens, toxicities, and
    manufacturing terms, returning combined context.

    Args:
        query: User question about CAR-T therapy.

    Returns:
        Combined knowledge context string.
    """
    query_upper = query.upper()
    query_lower = query.lower()
    sections = []

    # Check targets
    for antigen in CART_TARGETS:
        if antigen.upper() in query_upper:
            ctx = get_target_context(antigen)
            if ctx:
                sections.append(ctx)

    # Check toxicities
    tox_aliases = {
        "CRS": ["CRS", "CYTOKINE RELEASE"],
        "ICANS": ["ICANS", "NEUROTOXICITY", "NEUROLOGIC"],
        "B_CELL_APLASIA": ["B-CELL APLASIA", "HYPOGAMMAGLOBULINEMIA", "AGAMMAGLOBULINEMIA"],
        "HLH_MAS": ["HLH", "MAS", "HEMOPHAGOCYTIC"],
        "CYTOPENIAS": ["CYTOPENIA", "NEUTROPENIA", "THROMBOCYTOPENIA"],
        "TLS": ["TUMOR LYSIS"],
        "GVHD": ["GVHD", "GRAFT-VERSUS-HOST"],
        "ON_TARGET_OFF_TUMOR": ["ON-TARGET OFF-TUMOR", "ON TARGET OFF TUMOR"],
        "COAGULOPATHY": ["COAGULOPATHY", "DIC", "FIBRINOGEN"],
        "CARDIAC_TOXICITY": ["CARDIAC", "TROPONIN", "CARDIOMYOPATHY"],
        "RENAL_TOXICITY": ["RENAL", "AKI", "KIDNEY"],
        "SECONDARY_MALIGNANCY": ["SECONDARY MALIGNANCY", "T-CELL LYMPHOMA", "INSERTIONAL"],
    }
    for tox_id, aliases in tox_aliases.items():
        if any(a in query_upper for a in aliases):
            ctx = get_toxicity_context(tox_id)
            if ctx:
                sections.append(ctx)

    # Check manufacturing
    mfg_keywords = {
        "lentiviral_transduction": ["TRANSDUCTION", "LENTIVIRAL", "LENTIVIRUS"],
        "retroviral_transduction": ["RETROVIRAL", "RETROVIRUS", "GAMMA-RETROVIRAL"],
        "t_cell_activation": ["ACTIVATION", "ANTI-CD3", "DYNABEADS"],
        "ex_vivo_expansion": ["EXPANSION", "EX VIVO", "BIOREACTOR", "G-REX"],
        "leukapheresis": ["LEUKAPHERESIS", "APHERESIS"],
        "cryopreservation": ["CRYOPRESERVATION", "CRYO", "THAW"],
        "release_testing": ["RELEASE TESTING", "QC", "POTENCY ASSAY"],
        "point_of_care_manufacturing": ["POINT OF CARE", "POC MANUFACTURING", "DECENTRALIZED"],
        "lymphodepletion": ["LYMPHODEPLETION", "FLUDARABINE", "CYCLOPHOSPHAMIDE"],
        "vein_to_vein_time": ["VEIN-TO-VEIN", "TURNAROUND", "MANUFACTURING TIME"],
        "non_viral_transposon": ["TRANSPOSON", "SLEEPING BEAUTY", "PIGGYBACK"],
        "mrna_electroporation": ["MRNA", "ELECTROPORATION"],
        "crispr_knock_in": ["CRISPR", "TRAC", "KNOCK-IN"],
        "ipsc_derived": ["IPSC", "IPS CELL"],
        "automated_manufacturing": ["AUTOMATED", "CLINIMACS", "PRODIGY"],
    }
    for proc_id, keywords in mfg_keywords.items():
        if any(kw in query_upper for kw in keywords):
            ctx = get_manufacturing_context(proc_id)
            if ctx:
                sections.append(ctx)

    # Check biomarkers
    for key in CART_BIOMARKERS:
        aliases = [key, CART_BIOMARKERS[key]["full_name"].lower()]
        if any(a.lower() in query_lower for a in aliases):
            ctx = get_biomarker_context(key)
            if ctx:
                sections.append(ctx)

    # Check regulatory (product names)
    for key in CART_REGULATORY:
        if key.lower() in query_lower or CART_REGULATORY[key]["generic_name"].lower() in query_lower:
            ctx = get_regulatory_context(key)
            if ctx:
                sections.append(ctx)

    # Check immunogenicity
    immuno_keywords = {
        "murine_scfv_immunogenicity": ["IMMUNOGENICITY", "ADA", "ANTI-DRUG ANTIBOD", "HAMA", "ANTI-MURINE"],
        "humanization_strategies": ["HUMANIZATION", "HUMANIZED", "CDR GRAFTING", "DEIMMUNIZ", "FRAMEWORK SHUFFL"],
        "ada_clinical_impact": ["ADA TITER", "NEUTRALIZING ANTIBOD", "ANTI-CAR ANTIBOD"],
        "hla_restricted_epitopes": ["HLA", "MHC", "T-CELL EPITOPE", "NETCHMII", "NETMHC", "EPIMATRIX"],
        "immunogenicity_testing": ["IMMUNOGENICITY TEST", "ELISPOT", "DC-T CELL", "MAPPS"],
        "allogeneic_hla_considerations": ["ALLOGENEIC HLA", "B2M KNOCKOUT", "HLA-E", "TRAC KNOCKOUT", "OFF-THE-SHELF"],
    }
    for immuno_id, keywords in immuno_keywords.items():
        if any(kw in query_upper for kw in keywords):
            ctx = get_immunogenicity_context(immuno_id)
            if ctx:
                sections.append(ctx)

    if not sections:
        return ""

    return "\n\n".join(sections)


def get_knowledge_stats() -> Dict[str, int]:
    """Return statistics about the CAR-T knowledge graph."""
    return {
        "target_antigens": len(CART_TARGETS),
        "targets_with_approved_products": sum(
            1 for t in CART_TARGETS.values() if t.get("approved_products")
        ),
        "toxicity_profiles": len(CART_TOXICITIES),
        "manufacturing_processes": len(CART_MANUFACTURING),
        "biomarkers": len(CART_BIOMARKERS),
        "regulatory_products": len(CART_REGULATORY),
        "immunogenicity_topics": len(CART_IMMUNOGENICITY),
    }


# =============================================================================
# 4. CART_BIOMARKERS — Biomarker knowledge graph (~23 entries)
# =============================================================================

CART_BIOMARKERS: Dict[str, Dict[str, Any]] = {
    "ferritin": {
        "full_name": "Serum Ferritin",
        "type": "predictive",
        "assay_method": "Immunoassay (ELISA/CLIA)",
        "clinical_cutoff": ">500 mg/L pre-infusion",
        "predictive_value": "Elevated baseline ferritin (>500 mg/L) predicts grade 3+ CRS with PPV 72%",
        "associated_outcome": "CRS severity",
        "evidence_level": "validated",
        "key_references": ["PMID:29385376", "PMID:30409105"],
    },
    "crp": {
        "full_name": "C-Reactive Protein",
        "type": "predictive",
        "assay_method": "Immunoturbidimetry",
        "clinical_cutoff": ">200 mg/L within 72h post-infusion",
        "predictive_value": "Peak CRP >200 mg/L within 72h has 85% sensitivity for grade 3+ CRS",
        "associated_outcome": "CRS severity and timing",
        "evidence_level": "validated",
        "key_references": ["PMID:29385376", "PMID:31533922"],
    },
    "il6": {
        "full_name": "Interleukin-6",
        "type": "pharmacodynamic",
        "assay_method": "ELISA / Luminex multiplex",
        "clinical_cutoff": ">1000 pg/mL",
        "predictive_value": "Peak IL-6 correlates with CRS grade (r=0.82). Target of tocilizumab therapy",
        "associated_outcome": "CRS grade, tocilizumab response",
        "evidence_level": "validated",
        "key_references": ["PMID:29385376", "PMID:30409105"],
    },
    "sil2r": {
        "full_name": "Soluble IL-2 Receptor (sCD25)",
        "type": "pharmacodynamic",
        "assay_method": "ELISA",
        "clinical_cutoff": ">5000 pg/mL",
        "predictive_value": "Elevated sIL-2R correlates with T-cell activation and CRS/HLH risk",
        "associated_outcome": "CRS, HLH/MAS",
        "evidence_level": "validated",
        "key_references": ["PMID:30409105"],
    },
    "car_t_expansion": {
        "full_name": "CAR-T Cell Peak Expansion (Cmax)",
        "type": "pharmacodynamic",
        "assay_method": "qPCR (transgene copies/μg DNA) or flow cytometry",
        "clinical_cutoff": "Peak >50,000 copies/μg DNA (varies by product)",
        "predictive_value": "Higher peak expansion correlates with response (OR 3.2) but also CRS risk",
        "associated_outcome": "Response rate, CRS, durability",
        "evidence_level": "validated",
        "key_references": ["PMID:29385376", "PMID:28687837"],
    },
    "tcm_percentage": {
        "full_name": "Central Memory T-cell Percentage (Tcm%)",
        "type": "predictive",
        "assay_method": "Flow cytometry (CD45RA-/CCR7+/CD62L+)",
        "clinical_cutoff": ">40% Tcm in apheresis product",
        "predictive_value": "Higher Tcm% predicts better expansion, persistence, and response",
        "associated_outcome": "Manufacturing success, clinical response",
        "evidence_level": "emerging",
        "key_references": ["PMID:31040380"],
    },
    "cd4_cd8_ratio": {
        "full_name": "CD4:CD8 T-cell Ratio",
        "type": "predictive",
        "assay_method": "Flow cytometry",
        "clinical_cutoff": "Optimal 1:1 (defined CD4/CD8 composition in product)",
        "predictive_value": "Balanced CD4:CD8 ratio in product associated with improved persistence (lisocabtagene)",
        "associated_outcome": "CAR-T persistence, reduced toxicity",
        "evidence_level": "emerging",
        "key_references": ["PMID:33888460"],
    },
    "ldh": {
        "full_name": "Lactate Dehydrogenase",
        "type": "prognostic",
        "assay_method": "Enzymatic assay (serum)",
        "clinical_cutoff": ">ULN (upper limit of normal) at lymphodepletion",
        "predictive_value": "Elevated LDH (>ULN) pre-lymphodepletion predicts inferior PFS (HR 1.8)",
        "associated_outcome": "PFS, OS, tumor burden",
        "evidence_level": "validated",
        "key_references": ["PMID:28687837", "PMID:31040380"],
    },
    "pd1": {
        "full_name": "Programmed Death-1 (PD-1/CD279)",
        "type": "resistance",
        "assay_method": "Flow cytometry / IHC",
        "clinical_cutoff": ">30% PD-1+ on CAR-T cells",
        "predictive_value": "High PD-1 expression on infused product predicts shorter persistence and inferior response",
        "associated_outcome": "T-cell exhaustion, relapse",
        "evidence_level": "emerging",
        "key_references": ["PMID:31040380"],
    },
    "lag3": {
        "full_name": "Lymphocyte Activation Gene-3 (LAG-3/CD223)",
        "type": "resistance",
        "assay_method": "Flow cytometry",
        "clinical_cutoff": ">20% LAG-3+ on CAR-T cells",
        "predictive_value": "Co-expression with PD-1 and TIM-3 marks terminally exhausted T-cells",
        "associated_outcome": "T-cell exhaustion, reduced efficacy",
        "evidence_level": "emerging",
        "key_references": ["PMID:31040380"],
    },
    "tim3": {
        "full_name": "T-cell Immunoglobulin and Mucin-3 (TIM-3)",
        "type": "resistance",
        "assay_method": "Flow cytometry",
        "clinical_cutoff": ">25% TIM-3+ on CAR-T cells",
        "predictive_value": "TIM-3+ CAR-T cells show reduced cytokine production and cytotoxicity",
        "associated_outcome": "T-cell exhaustion, functional impairment",
        "evidence_level": "emerging",
        "key_references": ["PMID:31040380"],
    },
    "mrd_flow": {
        "full_name": "Minimal Residual Disease by Flow Cytometry",
        "type": "monitoring",
        "assay_method": "Multiparameter flow cytometry (≥8-color)",
        "clinical_cutoff": "<10^-4 (0.01%)",
        "predictive_value": "MRD negativity at day 28 post-infusion predicts PFS (HR 0.3 vs MRD+)",
        "associated_outcome": "PFS, relapse-free survival",
        "evidence_level": "validated",
        "key_references": ["PMID:29385376", "PMID:30409105"],
    },
    "ctdna": {
        "full_name": "Circulating Tumor DNA",
        "type": "monitoring",
        "assay_method": "NGS panel / ddPCR",
        "clinical_cutoff": "2-log reduction from baseline by day 28",
        "predictive_value": "Early ctDNA clearance correlates with durable CR in DLBCL post-CAR-T",
        "associated_outcome": "Response durability, early relapse detection",
        "evidence_level": "emerging",
        "key_references": ["PMID:34516408"],
    },
    "sbcma": {
        "full_name": "Soluble BCMA (sBCMA)",
        "type": "resistance",
        "assay_method": "ELISA",
        "clinical_cutoff": ">40 ng/mL baseline",
        "predictive_value": "High sBCMA acts as decoy, sequestering CAR-T cells and reducing on-target activity",
        "associated_outcome": "BCMA CAR-T response, resistance mechanism",
        "evidence_level": "emerging",
        "key_references": ["PMID:33657416"],
    },
    "ifn_gamma": {
        "full_name": "Interferon-gamma (IFN-γ)",
        "type": "pharmacodynamic",
        "assay_method": "ELISA / Luminex / ELISpot",
        "clinical_cutoff": ">500 pg/mL at day 7",
        "predictive_value": "Peak IFN-γ correlates with CAR-T effector function and antitumor activity",
        "associated_outcome": "Antitumor response, immune activation",
        "evidence_level": "validated",
        "key_references": ["PMID:28687837"],
    },
    "tox": {
        "full_name": "TOX (Thymocyte Selection Associated HMG Box)",
        "type": "resistance",
        "assay_method": "Flow cytometry (intracellular), scRNA-seq",
        "clinical_cutoff": ">50% TOX+ in CAR-T product",
        "predictive_value": "PPV 72% for poor persistence at 6 months",
        "associated_outcome": "TOX high expression correlates with terminal T-cell exhaustion and reduced antitumor activity",
        "evidence_level": "emerging",
        "key_references": ["Khan et al. Nature 2019 PMID:31207603", "Scott et al. Nature 2019 PMID:31207604"],
    },
    "nr4a": {
        "full_name": "NR4A Family (NR4A1/NR4A2/NR4A3)",
        "type": "resistance",
        "assay_method": "Flow cytometry, transcriptomic profiling",
        "clinical_cutoff": "NR4A triple-positive >30% of CAR-T cells",
        "predictive_value": "NR4A knockout CAR-T show 3-5x improved tumor rejection in preclinical models",
        "associated_outcome": "NR4A transcription factors enforce T-cell exhaustion program downstream of chronic antigen stimulation",
        "evidence_level": "emerging",
        "key_references": ["Chen et al. Nature 2019 PMID:30814732", "Liu et al. Nature 2019 PMID:30814735"],
    },
    "il10": {
        "full_name": "Interleukin-10 (IL-10)",
        "type": "pharmacodynamic",
        "assay_method": "ELISA, Luminex multiplex",
        "clinical_cutoff": ">100 pg/mL within 72h post-infusion",
        "predictive_value": "PPV 68% for concurrent or subsequent ICANS when IL-10 >100 pg/mL",
        "associated_outcome": "IL-10 elevation predicts neurotoxicity independent of CRS grade; marker of monocyte/macrophage activation",
        "evidence_level": "validated",
        "key_references": ["Santomasso et al. Cancer Discov 2018 PMID:30425106", "Gust et al. Cancer Discov 2017 PMID:29025771"],
    },
    "ang2": {
        "full_name": "Angiopoietin-2 (Ang-2)",
        "type": "pharmacodynamic",
        "assay_method": "ELISA",
        "clinical_cutoff": ">5 ng/mL pre-infusion or >10 ng/mL peak",
        "predictive_value": "PPV 75% for severe ICANS when Ang-2 elevated; Ang-2:Ang-1 ratio >2 correlates with BBB disruption",
        "associated_outcome": "Endothelial activation marker predicting blood-brain barrier breakdown and neurotoxicity severity",
        "evidence_level": "validated",
        "key_references": ["Gust et al. Cancer Discov 2017 PMID:29025771", "Hay KA et al. Blood 2017 PMID:28924019"],
    },
    "antigen_density": {
        "full_name": "Target Antigen Surface Density",
        "type": "predictive",
        "assay_method": "Quantitative flow cytometry (QIFIKIT/BD Quantibrite)",
        "clinical_cutoff": "Antigen-specific: CD19 >1,000 molecules/cell; BCMA >1,500 molecules/cell",
        "predictive_value": "Below threshold associated with 40-60% lower response rates; dose-response relationship for CAR-T killing",
        "associated_outcome": "Low antigen density reduces CAR-T avidity and serial killing capacity, enabling tumor escape",
        "evidence_level": "validated",
        "key_references": ["Majzner et al. Cancer Discov 2020 PMID:31974172", "Walker et al. Mol Ther 2017 PMID:28110864"],
    },
    "tumor_burden": {
        "full_name": "Baseline Tumor Burden (Metabolic Tumor Volume)",
        "type": "prognostic",
        "assay_method": "PET/CT (MTV calculation), LDH, bone marrow blast %",
        "clinical_cutoff": "MTV >80 cm\u00b3 or marrow blasts >50% = high burden; LDH >500 U/L",
        "predictive_value": "High tumor burden: 2-3x higher grade 3+ CRS risk, but paradoxically higher initial response rates",
        "associated_outcome": "Predicts both toxicity severity (CRS/TLS) and response depth; used for risk-stratified management",
        "evidence_level": "validated",
        "key_references": ["Locke et al. Blood 2017 PMID:28924018", "Park et al. NEJM 2018 PMID:29385376"],
    },
    "d_dimer": {
        "full_name": "D-Dimer",
        "type": "monitoring",
        "assay_method": "Immunoturbidimetric assay, latex agglutination",
        "clinical_cutoff": ">4x ULN (typically >2,000 ng/mL FEU)",
        "predictive_value": "D-dimer >4x ULN with grade 2+ CRS has 80% sensitivity for developing coagulopathy",
        "associated_outcome": "Early indicator of CRS-associated DIC; serial monitoring guides cryoprecipitate and CRS management escalation",
        "evidence_level": "emerging",
        "key_references": ["Fried et al. Blood Adv 2019 PMID:31189560", "Lee DW et al. Biol Blood Marrow Transplant 2019 PMID:30359826"],
    },
    "troponin": {
        "full_name": "High-Sensitivity Cardiac Troponin (hs-cTnI/hs-cTnT)",
        "type": "monitoring",
        "assay_method": "High-sensitivity immunoassay (Roche Elecsys, Abbott Architect)",
        "clinical_cutoff": ">99th percentile URL (14 ng/L for hs-cTnI); serial rise >20% = acute injury",
        "predictive_value": "Troponin elevation during CRS has 85% sensitivity for cardiac events; serial monitoring enables early intervention",
        "associated_outcome": "Detects subclinical myocardial injury during CRS; guides cardiology referral and ICU escalation",
        "evidence_level": "emerging",
        "key_references": ["Alvi et al. JACC CardioOncol 2021 PMID:34396344", "Lefebvre et al. Circ Heart Fail 2020 PMID:32634040"],
    },
}


# =============================================================================
# 5. CART_REGULATORY — Regulatory knowledge graph (~6 FDA-approved products)
# =============================================================================

CART_REGULATORY: Dict[str, Dict[str, Any]] = {
    "Kymriah": {
        "generic_name": "tisagenlecleucel",
        "manufacturer": "Novartis",
        "initial_approval": "2017-08-30",
        "initial_indication": "Pediatric/young adult r/r B-ALL",
        "pivotal_trial": "ELIANA (NCT02435849)",
        "designations": ["Breakthrough Therapy (2014)", "RMAT (not applicable, pre-RMAT)"],
        "subsequent_approvals": [
            {"date": "2018-05-01", "indication": "Adult r/r DLBCL (2+ prior lines)", "trial": "JULIET"},
            {"date": "2022-05-27", "indication": "Adult r/r FL (3+ prior lines)", "trial": "ELARA"},
        ],
        "rems": "CAR-T REMS: CRS and neurological toxicity management certification required",
        "post_marketing": ["15-year follow-up study for secondary malignancies", "RCR testing protocol"],
        "ema_approval": "2018-08-23",
    },
    "Yescarta": {
        "generic_name": "axicabtagene ciloleucel",
        "manufacturer": "Kite/Gilead",
        "initial_approval": "2017-10-18",
        "initial_indication": "Adult r/r LBCL (2+ prior lines)",
        "pivotal_trial": "ZUMA-1 (NCT02348216)",
        "designations": ["Breakthrough Therapy (2015)", "Priority Review"],
        "subsequent_approvals": [
            {"date": "2021-03-05", "indication": "Adult r/r LBCL (2L after 1st-line chemoimmunotherapy)", "trial": "ZUMA-7"},
            {"date": "2024-03-15", "indication": "Adult r/r FL (3+ prior lines)", "trial": "ZUMA-5"},
        ],
        "rems": "CAR-T REMS: Shared REMS program with other CAR-T products",
        "post_marketing": ["15-year follow-up study", "Post-marketing observational study"],
        "ema_approval": "2018-08-27",
    },
    "Tecartus": {
        "generic_name": "brexucabtagene autoleucel",
        "manufacturer": "Kite/Gilead",
        "initial_approval": "2020-07-24",
        "initial_indication": "Adult r/r MCL",
        "pivotal_trial": "ZUMA-2 (NCT02601313)",
        "designations": ["Breakthrough Therapy (2018)", "Priority Review", "Orphan Drug"],
        "subsequent_approvals": [
            {"date": "2021-10-01", "indication": "Adult r/r B-ALL", "trial": "ZUMA-3"},
        ],
        "rems": "CAR-T REMS",
        "post_marketing": ["15-year follow-up study"],
        "ema_approval": "2020-12-14",
    },
    "Breyanzi": {
        "generic_name": "lisocabtagene maraleucel",
        "manufacturer": "Bristol Myers Squibb/Juno",
        "initial_approval": "2021-02-05",
        "initial_indication": "Adult r/r LBCL (2+ prior lines)",
        "pivotal_trial": "TRANSCEND (NCT02631044)",
        "designations": ["Breakthrough Therapy (2016)", "RMAT (2018)", "Priority Review"],
        "subsequent_approvals": [
            {"date": "2022-06-24", "indication": "Adult r/r LBCL (2L)", "trial": "TRANSFORM"},
        ],
        "rems": "CAR-T REMS",
        "post_marketing": ["15-year follow-up study", "TRANSCEND WORLD long-term follow-up"],
        "ema_approval": "2022-04-04",
    },
    "Abecma": {
        "generic_name": "idecabtagene vicleucel",
        "manufacturer": "Bristol Myers Squibb/bluebird bio",
        "initial_approval": "2021-03-26",
        "initial_indication": "Adult r/r multiple myeloma (4+ prior lines incl. PI, IMiD, anti-CD38)",
        "pivotal_trial": "KarMMa (NCT03361748)",
        "designations": ["Breakthrough Therapy (2017)", "Priority Review", "Orphan Drug"],
        "subsequent_approvals": [],
        "rems": "CAR-T REMS",
        "post_marketing": ["KarMMa-3 confirmatory trial", "15-year follow-up study"],
        "ema_approval": "2021-08-18",
    },
    "Carvykti": {
        "generic_name": "ciltacabtagene autoleucel",
        "manufacturer": "Janssen/Legend Biotech",
        "initial_approval": "2022-02-28",
        "initial_indication": "Adult r/r multiple myeloma (4+ prior lines incl. PI, IMiD, anti-CD38)",
        "pivotal_trial": "CARTITUDE-1 (NCT03548207)",
        "designations": ["Breakthrough Therapy (2019)", "Priority Review", "Orphan Drug"],
        "subsequent_approvals": [
            {"date": "2024-04-04", "indication": "Adult r/r MM (1-3 prior lines incl. PI and IMiD)", "trial": "CARTITUDE-4"},
        ],
        "rems": "CAR-T REMS",
        "post_marketing": ["CARTITUDE-2 expansion cohorts", "15-year follow-up study"],
        "ema_approval": "2022-05-26",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# COMPARATIVE ANALYSIS — Entity Resolution
# ═══════════════════════════════════════════════════════════════════════

ENTITY_ALIASES: Dict[str, Dict[str, str]] = {
    # FDA-approved products → canonical name + target antigen
    "KYMRIAH": {"type": "product", "canonical": "Kymriah (tisagenlecleucel)", "target": "CD19"},
    "TISAGENLECLEUCEL": {"type": "product", "canonical": "Kymriah (tisagenlecleucel)", "target": "CD19"},
    "YESCARTA": {"type": "product", "canonical": "Yescarta (axicabtagene ciloleucel)", "target": "CD19"},
    "AXICABTAGENE": {"type": "product", "canonical": "Yescarta (axicabtagene ciloleucel)", "target": "CD19"},
    "TECARTUS": {"type": "product", "canonical": "Tecartus (brexucabtagene autoleucel)", "target": "CD19"},
    "BREXUCABTAGENE": {"type": "product", "canonical": "Tecartus (brexucabtagene autoleucel)", "target": "CD19"},
    "BREYANZI": {"type": "product", "canonical": "Breyanzi (lisocabtagene maraleucel)", "target": "CD19"},
    "LISOCABTAGENE": {"type": "product", "canonical": "Breyanzi (lisocabtagene maraleucel)", "target": "CD19"},
    "ABECMA": {"type": "product", "canonical": "Abecma (idecabtagene vicleucel)", "target": "BCMA"},
    "IDECABTAGENE": {"type": "product", "canonical": "Abecma (idecabtagene vicleucel)", "target": "BCMA"},
    "CARVYKTI": {"type": "product", "canonical": "Carvykti (ciltacabtagene autoleucel)", "target": "BCMA"},
    "CILTACABTAGENE": {"type": "product", "canonical": "Carvykti (ciltacabtagene autoleucel)", "target": "BCMA"},
    # Costimulatory domains
    "4-1BB": {"type": "costimulatory", "canonical": "4-1BB (CD137)"},
    "CD137": {"type": "costimulatory", "canonical": "4-1BB (CD137)"},
    "CD28": {"type": "costimulatory", "canonical": "CD28"},
    "OX40": {"type": "costimulatory", "canonical": "OX40 (CD134)"},
    "ICOS": {"type": "costimulatory", "canonical": "ICOS"},
    # Vector types
    "LENTIVIRAL": {"type": "manufacturing", "canonical": "lentiviral_transduction"},
    "RETROVIRAL": {"type": "manufacturing", "canonical": "retroviral_transduction"},
    # Biomarker aliases
    "FERRITIN": {"type": "biomarker", "canonical": "ferritin"},
    "CRP": {"type": "biomarker", "canonical": "crp"},
    "IL-6": {"type": "biomarker", "canonical": "il6"},
    "IL6": {"type": "biomarker", "canonical": "il6"},
    "PD-1": {"type": "biomarker", "canonical": "pd1"},
    "PD1": {"type": "biomarker", "canonical": "pd1"},
    "LAG-3": {"type": "biomarker", "canonical": "lag3"},
    "LAG3": {"type": "biomarker", "canonical": "lag3"},
    "TIM-3": {"type": "biomarker", "canonical": "tim3"},
    "TIM3": {"type": "biomarker", "canonical": "tim3"},
    "MRD": {"type": "biomarker", "canonical": "mrd_flow"},
    "SBCMA": {"type": "biomarker", "canonical": "sbcma"},
    "CTDNA": {"type": "biomarker", "canonical": "ctdna"},
    # Immunogenicity aliases
    "ADA": {"type": "immunogenicity", "canonical": "murine_scfv_immunogenicity"},
    "ANTI-DRUG ANTIBODY": {"type": "immunogenicity", "canonical": "murine_scfv_immunogenicity"},
    "HAMA": {"type": "immunogenicity", "canonical": "murine_scfv_immunogenicity"},
    "HUMANIZATION": {"type": "immunogenicity", "canonical": "humanization_strategies"},
    "DEIMMUNIZATION": {"type": "immunogenicity", "canonical": "humanization_strategies"},
    "HLA": {"type": "immunogenicity", "canonical": "hla_restricted_epitopes"},
    "MHC": {"type": "immunogenicity", "canonical": "hla_restricted_epitopes"},
    # New target aliases
    "FCRH5": {"type": "target", "canonical": "FcRH5"},
    "FCRL5": {"type": "target", "canonical": "FcRH5"},
    "CS1": {"type": "target", "canonical": "SLAMF7"},
    "CD319": {"type": "target", "canonical": "SLAMF7"},
    "CLEC12A": {"type": "target", "canonical": "CLL1"},
    "CD327": {"type": "target", "canonical": "CLL1"},
    "CD135": {"type": "target", "canonical": "FLT3"},
    "TACSTD2": {"type": "target", "canonical": "TROP2"},
    "CD326": {"type": "target", "canonical": "EpCAM"},
    # New biomarker aliases
    "TOX": {"type": "biomarker", "canonical": "tox"},
    "NR4A": {"type": "biomarker", "canonical": "nr4a"},
    "NR4A1": {"type": "biomarker", "canonical": "nr4a"},
    "D-DIMER": {"type": "biomarker", "canonical": "d_dimer"},
    "TROPONIN": {"type": "biomarker", "canonical": "troponin"},
    "ANGIOPOIETIN": {"type": "biomarker", "canonical": "ang2"},
    # Pediatric CAR-T aliases
    "PEDIATRIC ALL": {"type": "pediatric_cart", "canonical": "cd19_pediatric_all"},
    "PEDIATRIC B-ALL": {"type": "pediatric_cart", "canonical": "cd19_pediatric_all"},
    "ELIANA": {"type": "pediatric_cart", "canonical": "cd19_pediatric_all"},
    "PEDIATRIC CD22": {"type": "pediatric_cart", "canonical": "cd22_pediatric_relapse"},
    "PEDIATRIC NEUROBLASTOMA": {"type": "pediatric_cart", "canonical": "gd2_neuroblastoma"},
    "PEDIATRIC HODGKIN": {"type": "pediatric_cart", "canonical": "cd30_pediatric_hodgkin"},
    "PEDIATRIC CRS": {"type": "pediatric_cart", "canonical": "pediatric_crs_management"},
    "PEDIATRIC ICANS": {"type": "pediatric_cart", "canonical": "pediatric_icans_management"},
}


# =============================================================================
# 7. PEDIATRIC_CART — Pediatric CAR-T Knowledge Graph
# =============================================================================

PEDIATRIC_CART: Dict[str, Dict[str, Any]] = {
    "cd19_pediatric_all": {
        "target": "CD19",
        "indication": "Pediatric/Young Adult Relapsed/Refractory B-ALL",
        "age_range": "Up to 25 years",
        "approved_product": "Tisagenlecleucel (Kymriah)",
        "fda_approval": "August 2017 (first-ever CAR-T approval)",
        "pivotal_trial": {
            "name": "ELIANA (NCT02435849)",
            "reference": "Maude et al. NEJM 2018",
            "population": "Pediatric/young adult r/r B-ALL, ages 3-25",
            "n_enrolled": 75,
            "complete_remission_rate": "82% (61/75 evaluable patients)",
            "mrd_negativity": "100% of responders achieved MRD-negative CR",
            "12_month_efs": "50%",
            "12_month_os": "76%",
        },
        "pediatric_crs": {
            "any_grade": "77%",
            "grade_3_4": "47%",
            "median_onset": "Day 3 (range 1-22)",
            "median_duration": "8 days",
            "icu_admission": "~47% of patients",
            "tocilizumab_use": "69% of CRS patients",
            "vasopressor_use": "35%",
        },
        "pediatric_icans": {
            "any_grade": "40%",
            "grade_3_4": "13%",
            "median_onset": "Day 6",
            "seizures": "5%",
        },
        "manufacturing": {
            "turnaround": "~22 days from apheresis to infusion",
            "manufacturing_success": "~92%",
            "bridging_chemotherapy": "Often required during manufacturing wait",
        },
        "special_considerations": [
            "Lymphodepletion: fludarabine 30 mg/m²/day x3 + cyclophosphamide 500 mg/m²/day x2",
            "B-cell aplasia requiring IVIG replacement (target IgG >400 mg/dL)",
            "Prolonged cytopenias: 35% grade 3-4 at day 28",
            "Infection risk: PJP prophylaxis, antifungal prophylaxis, monitor for viral reactivation",
            "Long-term monitoring: 15-year follow-up per REMS requirement",
            "Growth and development monitoring in pediatric patients",
            "Fertility counseling prior to lymphodepletion (cyclophosphamide gonadotoxicity)",
        ],
        "description": "Tisagenlecleucel (Kymriah) was the first FDA-approved CAR-T product, "
                       "indicated for patients up to 25 years with r/r B-ALL. The ELIANA trial "
                       "demonstrated 82% CR rate with durable responses. Pediatric CRS is "
                       "common (77%) and requires specialized management with weight-based "
                       "tocilizumab dosing.",
    },
    "cd22_pediatric_relapse": {
        "target": "CD22",
        "indication": "Pediatric r/r B-ALL with CD19-Negative Relapse",
        "age_range": "Pediatric and young adult",
        "approved_product": None,
        "investigational_status": "Phase I/II trials",
        "rationale": "CD19 antigen loss occurs in 10-20% of patients post-CD19 CAR-T, "
                     "necessitating alternative targets. CD22 is expressed on most B-ALL "
                     "blasts and is retained in many CD19-negative relapses.",
        "key_trial": {
            "name": "NCT02315612 (NIH/NCI)",
            "reference": "Fry et al. Nature Medicine 2018",
            "population": "Pediatric/young adult r/r B-ALL (including post-CD19 CAR-T relapse)",
            "complete_remission_rate": "73% (11/15 patients)",
            "mrd_negativity": "82% of responders",
            "cd22_site_density": "Response correlated with CD22 site density (>2,000 molecules/cell)",
            "relapse_mechanism": "CD22 downregulation (diminished site density) — not complete loss",
        },
        "dual_targeting_approaches": [
            "Sequential CD19 → CD22 CAR-T: treat CD19-negative relapse",
            "Bivalent CD19/CD22 CAR-T: single construct targeting both antigens",
            "Co-infusion: separate CD19 and CD22 CAR-T products",
            "Loop CAR: tandem CD19-CD22 single-chain construct",
        ],
        "special_considerations": [
            "CD22 surface density matters: low density → reduced efficacy",
            "CD22 is not fully lost (unlike CD19) — downregulation rather than deletion",
            "Dose optimization ongoing: higher doses may overcome low CD22 density",
        ],
        "description": "CD22 CAR-T is an emerging approach for pediatric B-ALL patients who "
                       "relapse with CD19 antigen loss after CD19 CAR-T therapy. Fry et al. "
                       "demonstrated 73% CR rate. Dual CD19/CD22 targeting strategies are "
                       "being developed to prevent antigen escape.",
    },
    "dual_cd19_cd22": {
        "target": "CD19/CD22 (Dual)",
        "indication": "Pediatric r/r B-ALL — Prevention of Antigen Escape",
        "age_range": "Pediatric and young adult",
        "approved_product": None,
        "investigational_status": "Multiple Phase I/II trials",
        "approaches": {
            "bivalent": "Single CAR construct with both CD19 and CD22 binding domains",
            "co_administration": "Two separate CAR-T products infused together",
            "sequential": "CD19 CAR-T followed by CD22 CAR-T at relapse",
        },
        "key_trials": [
            "NCT03241940 (bivalent CD19/CD22, Stanford)",
            "NCT03448393 (dual CD19/CD22, CHOP)",
            "NCT03289455 (AUTO3, Autolus — sequential)",
        ],
        "rationale": "Dual antigen targeting reduces risk of antigen-negative relapse. "
                     "Bivalent constructs may provide more durable responses than single-target "
                     "approaches in pediatric ALL.",
        "description": "Dual CD19/CD22 CAR-T targeting is designed to prevent the 10-20% antigen "
                       "escape seen with single-target CD19 CAR-T. Multiple approaches are under "
                       "investigation including bivalent constructs and sequential administration.",
    },
    "cd30_pediatric_hodgkin": {
        "target": "CD30",
        "indication": "Pediatric Relapsed/Refractory Hodgkin Lymphoma",
        "age_range": "Pediatric and adolescent",
        "approved_product": None,
        "investigational_status": "Phase I/II trials",
        "key_trials": [
            "RELY-30 (NCT02690545)",
            "NCT04268706 (CD30 CAR-T, Baylor)",
        ],
        "toxicity_profile": {
            "CRS": "Moderate (typically grade 1-2)",
            "ICANS": "Low",
        },
        "special_considerations": [
            "CD30 is expressed on Reed-Sternberg cells but also on activated T-cells",
            "Brentuximab vedotin (ADC) failures may still respond to CD30 CAR-T",
            "Lower toxicity profile compared to CD19 CAR-T in B-ALL",
        ],
        "description": "CD30 CAR-T is under investigation for pediatric Hodgkin lymphoma "
                       "patients who fail brentuximab vedotin and checkpoint inhibitors. "
                       "Early data show manageable toxicity and promising responses.",
    },
    "gd2_neuroblastoma": {
        "target": "GD2",
        "indication": "High-Risk Neuroblastoma (Investigational CAR-T)",
        "age_range": "Pediatric (median age 2-5 years)",
        "approved_product": None,
        "investigational_status": "Phase I/II trials (CAR-T); dinutuximab (anti-GD2 mAb) is FDA-approved",
        "key_trials": [
            "NCT03635632 (GD2 CAR-T with C7R armoring, Baylor)",
            "NCT04196413 (C7R-GD2 armored CAR-T)",
            "NCT05298995 (GD2 CAR-T DIPG/DMG, Stanford)",
        ],
        "challenges": [
            "Solid tumor microenvironment: immunosuppressive (TGF-beta, IL-10)",
            "GD2 expression on CNS neurons → neuropathic pain (dose-limiting)",
            "Poor T-cell trafficking to solid tumor sites",
            "T-cell exhaustion in hostile TME",
        ],
        "armoring_strategies": [
            "C7R (constitutive IL-7 receptor): enhances persistence in solid tumors",
            "IL-15 co-expression: improves T-cell fitness",
            "Dominant-negative TGF-beta receptor: overcomes TME suppression",
        ],
        "special_considerations": [
            "GD2 mAb (dinutuximab) is standard-of-care for high-risk neuroblastoma maintenance",
            "CAR-T may overcome limitations of antibody therapy (penetration, persistence)",
            "DIPG/DMG: intracavitary delivery being explored to bypass BBB",
            "Neuropathic pain management critical in pediatric patients",
        ],
        "description": "GD2 CAR-T for neuroblastoma is investigational but represents a major "
                       "pediatric solid tumor CAR-T target. Key challenges include TME suppression "
                       "and on-target/off-tumor neuropathic pain. Armored CAR-T constructs (C7R) "
                       "are designed to enhance persistence in the solid tumor setting.",
    },
}


# =============================================================================
# 8. PEDIATRIC CRS/ICANS MANAGEMENT — Toxicity Management for Children
# =============================================================================

PEDIATRIC_CRS_ICANS_MANAGEMENT: Dict[str, Dict[str, Any]] = {
    "pediatric_crs_management": {
        "topic": "Pediatric CRS Management",
        "description": "CRS management in pediatric CAR-T patients requires weight-based "
                       "dosing and earlier intervention thresholds compared to adults.",
        "tocilizumab_dosing": {
            "weight_lt_30kg": "12 mg/kg IV (higher weight-based dose for smaller children)",
            "weight_gte_30kg": "8 mg/kg IV (standard adult dosing)",
            "max_dose": "800 mg per dose",
            "repeat_dosing": "May repeat every 8 hours, maximum 3 doses in 24 hours",
            "max_total_doses": "4 doses total",
        },
        "intervention_thresholds": {
            "grade_1": "Supportive care, antipyretics, monitoring",
            "grade_2": "Tocilizumab (lower threshold than adults — treat at grade 2)",
            "grade_3": "Tocilizumab + corticosteroids (dexamethasone 10 mg/m² or methylprednisolone 2 mg/kg)",
            "grade_4": "Tocilizumab + high-dose corticosteroids + vasopressors + ICU",
        },
        "pediatric_specific": [
            "Earlier intervention threshold: tocilizumab at grade 2 CRS (vs grade 2-3 in adults)",
            "Children may decompensate faster than adults",
            "Fluid resuscitation: 10-20 mL/kg isotonic fluid boluses (weight-based)",
            "Vasopressor dosing: weight-based (mcg/kg/min)",
            "Temperature monitoring: continuous; children may have higher baseline fevers",
            "Parent/caregiver education on CRS recognition critical",
        ],
        "monitoring": [
            "Continuous telemetry and pulse oximetry",
            "Vital signs every 4 hours minimum (every 1-2 hours for grade ≥2)",
            "Daily labs: ferritin, CRP, LDH, comprehensive metabolic panel, CBC",
            "Cytokine panel (IL-6, IFN-gamma, IL-10) if available",
            "Echocardiography for grade ≥3 CRS",
        ],
    },
    "pediatric_icans_management": {
        "topic": "Pediatric ICANS (Immune Effector Cell-Associated Neurotoxicity Syndrome) Management",
        "description": "Neurotoxicity management in pediatric CAR-T patients includes "
                       "seizure prophylaxis and adapted assessment tools for younger children.",
        "assessment": {
            "ice_score": "ICE (Immune Effector Cell-Associated Encephalopathy) score — adapted for age",
            "pediatric_adaptation": "Cornell Assessment of Pediatric Delirium (CAPD) for age <12 years",
            "standard_ice": "Standard ICE score for age ≥12 years",
        },
        "seizure_prophylaxis": {
            "indication": "ICANS grade ≥2",
            "drug": "Levetiracetam (Keppra)",
            "pediatric_dose": "10-20 mg/kg/dose BID (max 1500 mg/dose)",
            "duration": "Continue for ≥72 hours after ICANS resolution",
            "breakthrough": "Benzodiazepines (lorazepam 0.1 mg/kg IV) for breakthrough seizures",
        },
        "management_by_grade": {
            "grade_1": "Observation, seizure prophylaxis, non-sedating environment",
            "grade_2": "Dexamethasone 0.5 mg/kg IV q6h (max 10 mg/dose) + seizure prophylaxis",
            "grade_3": "Dexamethasone 0.5 mg/kg IV q6h + consider methylprednisolone 2 mg/kg/day",
            "grade_4": "High-dose methylprednisolone 30 mg/kg/day (max 1g) + ICU + neurology consult",
        },
        "pediatric_specific": [
            "Seizure prophylaxis with levetiracetam for all ICANS grade ≥2",
            "Non-convulsive status epilepticus: continuous EEG monitoring if altered mental status",
            "Papilledema check: concern for cerebral edema (rare but fatal)",
            "MRI brain if focal deficits or grade ≥3 ICANS",
            "Neurodevelopmental assessment at follow-up (long-term cognitive monitoring)",
        ],
    },
    "pediatric_icu_monitoring": {
        "topic": "Pediatric ICU Monitoring Post-CAR-T Infusion",
        "description": "All pediatric patients should be monitored in a CAR-T-capable center "
                       "with pediatric ICU access for a minimum of 14 days post-infusion.",
        "monitoring_requirements": [
            "All pediatric patients: minimum 14 days post-infusion monitoring",
            "REMS-certified center with pediatric ICU capability",
            "Tocilizumab must be available on-site (minimum 2 doses per patient weight)",
            "Pediatric intensivist available 24/7",
            "Pediatric neurology consultation available within 1 hour",
        ],
        "discharge_criteria": [
            "No CRS ≥ grade 2 for ≥48 hours off treatment",
            "No ICANS ≥ grade 1 for ≥24 hours",
            "Tolerating oral medications and nutrition",
            "Reliable caregiver education completed (CRS/ICANS recognition)",
            "Within 1 hour of REMS-certified center for ≥4 weeks post-infusion",
            "Driving restriction: patient/caregiver informed of 8-week restriction",
        ],
        "infection_prophylaxis": [
            "PJP prophylaxis: TMP-SMX (dose: 5 mg/kg/day TMP divided BID, 3 days/week)",
            "Antiviral: acyclovir 10 mg/kg/dose TID (HSV/VZV prophylaxis)",
            "Antifungal: fluconazole 6 mg/kg/day (if prolonged neutropenia)",
            "IVIG replacement: target IgG >400 mg/dL (B-cell aplasia management)",
            "Vaccination: live vaccines contraindicated; re-immunization per post-SCT schedule",
        ],
    },
}


def get_pediatric_cart_context(topic: str) -> str:
    """Return formatted knowledge for a pediatric CAR-T topic.

    Args:
        topic: Topic key or keyword (e.g. 'cd19_pediatric_all', 'cd22',
               'neuroblastoma', 'pediatric crs').

    Returns:
        Formatted string with pediatric CAR-T knowledge, or empty string.
    """
    key = topic.lower().replace("-", "_").replace(" ", "_")

    # Check PEDIATRIC_CART first
    data = PEDIATRIC_CART.get(key)
    if data:
        lines = [
            f"## Pediatric CAR-T: {data.get('target', '')} — {data.get('indication', '')}",
            f"- Age Range: {data.get('age_range', 'N/A')}",
            f"- Approved Product: {data.get('approved_product') or 'Investigational'}",
        ]
        if data.get("pivotal_trial"):
            trial = data["pivotal_trial"]
            lines.append(f"- Pivotal Trial: {trial.get('name', '')} ({trial.get('reference', '')})")
            lines.append(f"  CR Rate: {trial.get('complete_remission_rate', 'N/A')}")
        if data.get("key_trial"):
            trial = data["key_trial"]
            lines.append(f"- Key Trial: {trial.get('name', '')} ({trial.get('reference', '')})")
            lines.append(f"  CR Rate: {trial.get('complete_remission_rate', 'N/A')}")
        if data.get("pediatric_crs"):
            crs = data["pediatric_crs"]
            lines.append(f"- Pediatric CRS: any grade {crs.get('any_grade', '')}, "
                         f"grade 3-4 {crs.get('grade_3_4', '')}")
        if data.get("pediatric_icans"):
            icans = data["pediatric_icans"]
            lines.append(f"- Pediatric ICANS: any grade {icans.get('any_grade', '')}, "
                         f"grade 3-4 {icans.get('grade_3_4', '')}")
        lines.append(f"- {data.get('description', '')}")
        return "\n".join(lines)

    # Check PEDIATRIC_CRS_ICANS_MANAGEMENT
    data = PEDIATRIC_CRS_ICANS_MANAGEMENT.get(key)
    if not data:
        for k, v in PEDIATRIC_CRS_ICANS_MANAGEMENT.items():
            if key in k or key in v.get("topic", "").lower().replace(" ", "_"):
                data = v
                break
    if not data:
        # Fuzzy search in PEDIATRIC_CART
        for k, v in PEDIATRIC_CART.items():
            if key in k or key in v.get("target", "").lower() or key in v.get("indication", "").lower():
                data = v
                return get_pediatric_cart_context(k)

    if data and "topic" in data:
        lines = [f"## {data['topic']}", f"- {data['description']}"]
        if "tocilizumab_dosing" in data:
            toci = data["tocilizumab_dosing"]
            lines.append("- Tocilizumab Dosing:")
            lines.append(f"  <30 kg: {toci['weight_lt_30kg']}")
            lines.append(f"  ≥30 kg: {toci['weight_gte_30kg']}")
        if "seizure_prophylaxis" in data:
            sp = data["seizure_prophylaxis"]
            lines.append(f"- Seizure Prophylaxis: {sp['drug']} {sp['pediatric_dose']}")
        if "pediatric_specific" in data:
            lines.append("- Pediatric-Specific Considerations:")
            for item in data["pediatric_specific"]:
                lines.append(f"  * {item}")
        return "\n".join(lines)

    return ""


def get_biomarker_context(biomarker: str) -> str:
    """Return formatted knowledge for a biomarker."""
    key = biomarker.lower().replace("-", "").replace(" ", "_")
    data = CART_BIOMARKERS.get(key)
    if not data:
        return ""
    lines = [
        f"Biomarker: {data['full_name']}",
        f"  Type: {data['type']}",
        f"  Assay Method: {data['assay_method']}",
        f"  Clinical Cutoff: {data['clinical_cutoff']}",
        f"  Predictive Value: {data['predictive_value']}",
        f"  Associated Outcome: {data['associated_outcome']}",
        f"  Evidence Level: {data['evidence_level']}",
    ]
    return "\n".join(lines)


def get_regulatory_context(product: str) -> str:
    """Return formatted regulatory knowledge for a product."""
    data = CART_REGULATORY.get(product)
    if not data:
        # Try matching by generic name
        for k, v in CART_REGULATORY.items():
            if v.get("generic_name", "").lower() in product.lower() or product.lower() in k.lower():
                data = v
                break
    if not data:
        return ""
    lines = [
        f"Regulatory Profile: {product}",
        f"  Generic Name: {data['generic_name']}",
        f"  Manufacturer: {data['manufacturer']}",
        f"  Initial FDA Approval: {data['initial_approval']}",
        f"  Initial Indication: {data['initial_indication']}",
        f"  Pivotal Trial: {data['pivotal_trial']}",
        f"  Designations: {', '.join(data['designations'])}",
        f"  REMS: {data['rems']}",
    ]
    if data.get("subsequent_approvals"):
        lines.append("  Subsequent Approvals:")
        for sa in data["subsequent_approvals"]:
            lines.append(f"    - {sa['date']}: {sa['indication']} ({sa['trial']})")
    if data.get("ema_approval"):
        lines.append(f"  EMA Approval: {data['ema_approval']}")
    return "\n".join(lines)


# =============================================================================
# 6. CART_IMMUNOGENICITY — HLA & Immunogenicity knowledge (~6 entries)
# =============================================================================

CART_IMMUNOGENICITY: Dict[str, Dict[str, Any]] = {
    "murine_scfv_immunogenicity": {
        "topic": "Murine scFv Immunogenicity in CAR-T",
        "description": (
            "Murine-derived scFvs (e.g., FMC63 in Kymriah/Yescarta/Breyanzi, "
            "11D5-3 in Abecma) elicit anti-drug antibodies (ADA) in ~3-5% of patients. "
            "ADA can reduce CAR-T persistence, impair re-dosing efficacy, and rarely "
            "cause infusion reactions. ELIANA trial: ~5% ADA against FMC63."
        ),
        "key_constructs": ["FMC63 (murine anti-CD19)", "SJ25C1 (murine anti-CD19)",
                           "11D5-3 (murine anti-BCMA)", "14G2a (murine anti-GD2)"],
        "ada_incidence": "3-8% for murine scFvs, <1% for humanized, <0.5% for fully human",
        "clinical_impact": "Reduced persistence at 12 months; impaired re-dosing response",
        "management": "Monitor ADA titers; consider humanized alternatives for re-dosing",
    },
    "humanization_strategies": {
        "topic": "scFv Humanization Strategies for CAR-T",
        "description": (
            "Humanization reduces immunogenicity by grafting murine CDRs onto human "
            "framework regions. Key approaches: (1) CDR grafting onto closest human "
            "germline (VH3-23, VK1-39 are common acceptors), (2) framework back-mutations "
            "for affinity retention, (3) deimmunization via T-cell epitope removal, "
            "(4) fully human binders from phage/yeast display libraries."
        ),
        "methods": ["CDR grafting", "Framework shuffling", "Back-mutation optimization",
                     "Deimmunization (epitope deletion)", "Fully human library selection"],
        "tools": ["EpiMatrix/EpiVax", "NetMHCIIpan", "IEDB", "AbDesigner", "Rosetta"],
        "tradeoffs": "Humanization may reduce affinity (1.5-5x); back-mutations partially restore",
    },
    "ada_clinical_impact": {
        "topic": "Clinical Impact of Anti-Drug Antibodies on CAR-T",
        "description": (
            "ADA against CAR constructs can: (1) neutralize CAR-T cells (blocking "
            "scFv-antigen interaction), (2) accelerate clearance via Fc-mediated "
            "opsonization, (3) cause infusion reactions on re-dosing, (4) activate "
            "complement-dependent cytotoxicity of CAR-T cells. High-titer ADA "
            "(>1:1000) correlates with 40% reduced persistence at 6 months."
        ),
        "risk_factors": ["Murine scFv origin", "Repeated dosing", "Non-human hinge/spacer",
                         "High transgene immunogenicity", "Intact immune system (less lymphodepletion)"],
        "monitoring": "ELISA (screening) → confirmatory (drug tolerance) → titer → neutralizing assay",
    },
    "hla_restricted_epitopes": {
        "topic": "HLA-Restricted T-Cell Epitopes in CAR Constructs",
        "description": (
            "MHC class II presentation of CAR-derived peptides drives CD4 T-cell help "
            "for ADA production. HLA-DRB1*04:01 and HLA-DRB1*15:01 alleles show "
            "highest risk for FMC63-derived epitope presentation. Computational tools "
            "(NetMHCIIpan, EpiMatrix) predict 15-20 potential epitopes in murine scFvs "
            "vs 2-5 in humanized versions."
        ),
        "high_risk_alleles": ["HLA-DRB1*04:01", "HLA-DRB1*15:01", "HLA-DRB1*07:01"],
        "prediction_tools": ["NetMHCIIpan-4.0", "EpiMatrix", "IEDB MHC-II binding"],
        "epitope_counts": {"murine_FMC63": 18, "humanized_FMC63": 4, "fully_human_CD19": 2},
    },
    "immunogenicity_testing": {
        "topic": "Immunogenicity Testing Paradigm for CAR-T",
        "description": (
            "Three-tier testing: (1) In silico: EpiMatrix/NetMHCIIpan T-cell epitope "
            "prediction, aggregatrix score for epitope clustering. (2) In vitro: DC-T cell "
            "co-culture, ELISpot for IL-4/IFN-γ, MAPPs (MHC-associated peptide proteomics). "
            "(3) Clinical: ADA sampling at baseline, day 28, months 3/6/12; tiered approach "
            "per FDA 2024 guidance."
        ),
        "in_silico": ["EpiMatrix", "NetMHCIIpan", "IEDB", "iTope/TCED"],
        "in_vitro": ["DC-T cell co-culture", "ELISpot (IL-4, IFN-γ)", "MAPPs assay"],
        "clinical": ["Screening ELISA", "Confirmatory assay", "Titer", "Neutralizing antibody"],
        "fda_guidance": "FDA 2024: Immunogenicity Testing of Therapeutic Protein Products",
    },
    "allogeneic_hla_considerations": {
        "topic": "HLA Considerations in Allogeneic CAR-T",
        "description": (
            "Allogeneic (off-the-shelf) CAR-T requires HLA management to prevent "
            "graft-versus-host disease (GvHD) and host-versus-graft rejection. "
            "Strategies: (1) TRAC/TRBC knockout (TALEN or CRISPR) eliminates GvHD, "
            "(2) B2M knockout removes MHC-I to evade host CD8 T-cells, (3) HLA-E "
            "overexpression (B2M-HLA-E fusion) prevents NK cell lysis of MHC-I-negative "
            "cells, (4) CD52 knockout enables alemtuzumab-based lymphodepletion selectivity."
        ),
        "gene_edits": ["TRAC KO (prevents GvHD)", "B2M KO (evades CD8 rejection)",
                        "HLA-E knock-in (NK evasion)", "CD52 KO (alemtuzumab resistance)"],
        "platforms": ["TALEN (Cellectis UCART)", "CRISPR/Cas9 (CRISPR Therapeutics CTX110)",
                      "Base editing (Beam BEAM-201)", "iPSC-derived (Fate FT819)"],
    },
}


def get_immunogenicity_context(topic: str) -> str:
    """Return formatted knowledge for an immunogenicity topic.

    Args:
        topic: Topic key or keyword (e.g. 'humanization', 'ada', 'hla').

    Returns:
        Formatted string with immunogenicity knowledge, or empty string.
    """
    key = topic.lower().replace("-", "_").replace(" ", "_")
    data = CART_IMMUNOGENICITY.get(key)
    if not data:
        # Try partial match
        for k in CART_IMMUNOGENICITY:
            if key in k.lower():
                data = CART_IMMUNOGENICITY[k]
                break
    if not data:
        return ""

    lines = [f"## Immunogenicity: {data['topic']}", f"- {data['description']}"]
    for field in ["ada_incidence", "clinical_impact", "management", "tradeoffs",
                  "fda_guidance"]:
        if field in data:
            label = field.replace("_", " ").title()
            lines.append(f"- **{label}:** {data[field]}")
    return "\n".join(lines)


def resolve_comparison_entity(text: str) -> Optional[Dict[str, str]]:
    """Resolve a raw text string to a known CAR-T entity for comparison.

    Checks target antigens, product aliases, toxicities, and manufacturing
    processes in priority order.

    Returns:
        Dict with 'type', 'canonical', and optionally 'target' keys,
        or None if the entity is not recognized.
    """
    cleaned = text.strip().upper()

    # 1. Target antigens (exact match)
    for key in CART_TARGETS:
        if key.upper() == cleaned:
            return {"type": "target", "canonical": key, "target": key}

    # 2. Product / alias table
    if cleaned in ENTITY_ALIASES:
        return dict(ENTITY_ALIASES[cleaned])

    # 3. Toxicities
    for key in CART_TOXICITIES:
        if key.upper() == cleaned or key.upper().replace("_", " ") == cleaned:
            return {"type": "toxicity", "canonical": key}

    # 4. Manufacturing processes (fuzzy substring)
    for key in CART_MANUFACTURING:
        if cleaned.replace(" ", "_").lower() in key.lower() or key.lower().startswith(cleaned.lower()):
            return {"type": "manufacturing", "canonical": key}

    # 5. Biomarkers
    for key in CART_BIOMARKERS:
        if cleaned.lower() == key or CART_BIOMARKERS[key]["full_name"].upper() == cleaned:
            return {"type": "biomarker", "canonical": key}

    return None


def get_comparison_context(entity_a: Dict[str, str], entity_b: Dict[str, str]) -> str:
    """Build side-by-side knowledge graph context for two entities.

    Reuses existing get_target_context / get_toxicity_context /
    get_manufacturing_context depending on entity type.

    Returns:
        Formatted comparison context string with both entities' data.
    """
    def _get_entity_context(entity: Dict[str, str]) -> str:
        etype = entity["type"]
        canonical = entity["canonical"]
        if etype == "target":
            return get_target_context(canonical)
        elif etype == "product":
            target = entity.get("target", "")
            return get_target_context(target) if target else ""
        elif etype == "toxicity":
            return get_toxicity_context(canonical)
        elif etype == "manufacturing":
            return get_manufacturing_context(canonical)
        elif etype == "biomarker":
            return get_biomarker_context(canonical)
        elif etype == "costimulatory":
            return f"Costimulatory Domain: {canonical}"
        return ""

    sections = []
    ctx_a = _get_entity_context(entity_a)
    ctx_b = _get_entity_context(entity_b)

    if ctx_a:
        sections.append(f"### {entity_a['canonical']}\n{ctx_a}")
    if ctx_b:
        sections.append(f"### {entity_b['canonical']}\n{ctx_b}")

    return "\n\n---\n\n".join(sections)
