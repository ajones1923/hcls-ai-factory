"""CAR-T query expansion maps for multi-collection RAG search.

When a user asks about a CAR-T topic, expand the search to include
related terms, antigens, products, and concepts.  This improves recall
across the 5 CAR-T collections (literature, trials, constructs, assays,
manufacturing).

Pattern from: rag-chat-pipeline/src/rag_engine.py
  - 10 therapeutic-area expansion dictionaries
  - _get_expanded_genes() scans all maps for keyword hits

Author: Adam Jones
Date: February 2026
"""

from typing import Dict, List, Set

from loguru import logger


# ═══════════════════════════════════════════════════════════════════════
# 1. TARGET_ANTIGEN_EXPANSION
#    Maps antigen keywords to related terms (diseases, products, aliases)
# ═══════════════════════════════════════════════════════════════════════

TARGET_ANTIGEN_EXPANSION: Dict[str, List[str]] = {
    # --- B-cell antigens ---
    "cd19": [
        "CD19", "B-ALL", "DLBCL", "B-cell lymphoma", "B-cell leukemia",
        "tisagenlecleucel", "axicabtagene ciloleucel", "lisocabtagene maraleucel",
        "brexucabtagene autoleucel", "Kymriah", "Yescarta", "Breyanzi", "Tecartus",
        "FMC63", "SJ25C1",
    ],
    "cd22": [
        "CD22", "B-ALL", "DLBCL", "hairy cell leukemia",
        "inotuzumab ozogamicin", "moxetumomab pasudotox",
        "CD19/CD22 bispecific", "dual-targeting", "m971",
    ],
    "cd20": [
        "CD20", "DLBCL", "follicular lymphoma", "B-NHL",
        "rituximab", "obinutuzumab", "Leu-16",
    ],

    # --- Myeloma antigens ---
    "bcma": [
        "BCMA", "B-cell maturation antigen", "TNFRSF17",
        "multiple myeloma", "relapsed refractory myeloma",
        "idecabtagene vicleucel", "ciltacabtagene autoleucel",
        "Abecma", "Carvykti", "teclistamab", "APRIL", "BAFF",
        "bispecific T-cell engager",
    ],
    "gprc5d": [
        "GPRC5D", "multiple myeloma", "talquetamab",
        "BCMA/GPRC5D bispecific", "orphan GPCR",
    ],
    "cd38": [
        "CD38", "multiple myeloma", "daratumumab", "isatuximab",
        "Darzalex", "Sarclisa", "NAD+ hydrolase",
    ],

    # --- Lymphoma / Hodgkin ---
    "cd30": [
        "CD30", "TNFRSF8", "Hodgkin lymphoma", "HL",
        "anaplastic large cell lymphoma", "ALCL",
        "brentuximab vedotin", "Adcetris",
    ],

    # --- Myeloid / AML ---
    "cd33": [
        "CD33", "Siglec-3", "AML", "acute myeloid leukemia",
        "myelodysplastic syndrome", "MDS",
        "gemtuzumab ozogamicin", "Mylotarg",
    ],
    "cd123": [
        "CD123", "IL-3R-alpha", "IL3RA", "IL-3 receptor alpha", "interleukin-3 receptor",
        "AML", "AML target", "BPDCN", "blastic plasmacytoid dendritic cell neoplasm",
        "blastic plasmacytoid dendritic cell",
        "tagraxofusp", "Elzonris",
    ],
    "flt3": [
        "FLT3", "CD135", "FMS-like tyrosine kinase 3",
        "AML", "FLT3-ITD", "midostaurin", "gilteritinib",
        "FLT3-CAR-T", "AML FLT3", "AMG 553",
    ],
    "cll-1": [
        "CLL-1", "CLEC12A", "AML", "myeloid leukemia",
        "CLL1/CD33 bispecific",
    ],
    "cll1": [
        "CLL-1", "CLL1", "CLEC12A", "CD371", "C-type lectin", "myeloid CLL-1",
        "AML CLL-1 CAR-T", "leukemic stem cell marker", "AML LSC target",
    ],

    # --- T-cell malignancies ---
    "cd7": [
        "CD7", "T-ALL", "T-cell lymphoma",
        "fratricide", "CD7 knockout", "CRISPR-edited",
        "universal CAR-T", "allogeneic",
    ],
    "cd5": [
        "CD5", "T-ALL", "T-cell lymphoma",
        "fratricide resistance", "CD5 knockout",
    ],

    # --- Solid tumor antigens ---
    "her2": [
        "HER2", "ERBB2", "HER2-positive", "breast cancer",
        "gastric cancer", "glioblastoma", "osteosarcoma",
        "trastuzumab", "Herceptin", "4D5",
    ],
    "gpc3": [
        "GPC3", "glypican-3", "hepatocellular carcinoma", "HCC",
        "liver cancer", "pediatric liver tumors",
    ],
    "egfr": [
        "EGFR", "EGFRvIII", "glioblastoma", "GBM",
        "non-small cell lung cancer", "NSCLC",
        "cetuximab", "variant III",
    ],
    "mesothelin": [
        "mesothelin", "MSLN", "mesothelioma",
        "ovarian cancer", "pancreatic cancer",
        "SS1P", "anetumab ravtansine",
    ],
    "claudin18.2": [
        "Claudin 18.2", "CLDN18.2", "gastric cancer",
        "pancreatic cancer", "gastric adenocarcinoma",
        "zolbetuximab", "tight junction",
    ],
    "gd2": [
        "GD2", "disialoganglioside", "ganglioside", "neuroblastoma", "neuroblastoma target",
        "retinoblastoma", "osteosarcoma", "melanoma",
        "dinutuximab", "Unituxin", "14g2a",
        "DIPG target", "DMG",
    ],
    "psma": [
        "PSMA", "prostate-specific membrane antigen", "FOLH1",
        "prostate cancer", "mCRPC",
        "metastatic castration-resistant",
    ],
    "ror1": [
        "ROR1", "receptor tyrosine kinase-like orphan receptor 1",
        "CLL", "mantle cell lymphoma", "triple-negative breast cancer",
        "TNBC", "Wnt signaling",
    ],
    "dll3": [
        "DLL3", "delta-like ligand 3", "small cell lung cancer", "SCLC",
        "neuroendocrine tumors", "Notch ligand",
        "tarlatamab", "rovalpituzumab tesirine",
    ],
    "b7-h3": [
        "B7-H3", "CD276", "neuroblastoma", "NSCLC",
        "head and neck cancer", "checkpoint",
    ],
    "muc1": [
        "MUC1", "mucin-1", "breast cancer",
        "pancreatic cancer", "ovarian cancer",
        "Tn-MUC1", "tumor-associated glycoform",
    ],
    "il13ra2": [
        "IL13Ralpha2", "IL-13 receptor alpha-2",
        "glioblastoma", "GBM", "intracranial delivery",
    ],
    "epcam": [
        "EpCAM", "CD326", "TACSTD1", "epithelial cell adhesion molecule",
        "EpCAM CAR-T", "colorectal CAR-T", "peritoneal carcinomatosis",
        "colorectal cancer", "gastric cancer",
        "epithelial tumors",
    ],

    # --- Additional antigens ---
    "fcrh5": [
        "FcRH5", "FCRLA", "FCRL5", "Fc receptor-like 5", "MCARH109",
        "multiple myeloma FcRH5", "tecvayli target",
    ],
    "slamf7": [
        "SLAMF7", "CS1", "CD319", "SLAM family member 7", "UCART-CS1",
        "CARAMBA", "elotuzumab target",
    ],
    "cd70": [
        "CD70", "CD27 ligand", "TNFSF7", "CD27L", "PRGN-3006",
        "AML CD70 CAR-T", "renal cell carcinoma CD70",
    ],
    "trop2": [
        "TROP2", "TACSTD2", "trophoblast antigen", "sacituzumab target",
        "triple-negative breast CAR-T", "solid tumor TROP2",
    ],
    "cd44v6": [
        "CD44v6", "CD44 variant 6", "CD44v6 CAR-T", "MLM-CAR44.1",
        "AML CD44v6", "pancreatic CD44v6",
        "adhesion molecule", "hyaluronan receptor", "splice variant",
    ],

    # --- Additional targets (expanded) ---
    "nkg2d": [
        "NKG2D", "NKG2D ligand", "MICA", "MICB", "ULBP", "stress ligand",
        "natural killer group 2D", "innate immune receptor",
    ],
    "lngfr": [
        "LNGFR", "CD271", "p75NTR", "nerve growth factor receptor",
        "cancer stem cell marker", "melanoma stem cell",
    ],
    "prame": [
        "PRAME", "preferentially expressed antigen in melanoma",
        "cancer testis antigen", "melanoma antigen",
    ],
    "wt1": [
        "WT1", "Wilms tumor 1", "transcription factor",
        "leukemia antigen", "solid tumor antigen",
    ],
    "gcc": [
        "GCC", "guanylyl cyclase C", "GUCY2C",
        "colorectal target", "intestinal receptor",
    ],
    "lilrb4": [
        "LILRB4", "ILT3", "leukocyte immunoglobulin-like receptor",
        "monocytic AML target",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 2. DISEASE_EXPANSION
#    Maps disease keywords to related terms, antigens, and therapies
# ═══════════════════════════════════════════════════════════════════════

DISEASE_EXPANSION: Dict[str, List[str]] = {
    # --- Hematologic malignancies ---
    "b-all": [
        "B-ALL", "B-cell acute lymphoblastic leukemia",
        "precursor B-cell ALL", "Philadelphia-positive ALL", "Ph+ ALL",
        "CD19", "CD22", "CD19/CD22 dual",
        "tisagenlecleucel", "Kymriah",
        "blinatumomab", "Blincyto", "MRD",
        "minimal residual disease", "pediatric ALL",
    ],
    "t-all": [
        "T-ALL", "T-cell acute lymphoblastic leukemia",
        "T-lymphoblastic leukemia",
        "CD7", "CD5", "CD1a",
        "fratricide", "nelarabine",
    ],
    "dlbcl": [
        "DLBCL", "diffuse large B-cell lymphoma",
        "relapsed refractory DLBCL", "r/r DLBCL",
        "CD19", "axicabtagene ciloleucel", "lisocabtagene maraleucel",
        "Yescarta", "Breyanzi",
        "GCB subtype", "ABC subtype", "double-hit lymphoma",
        "high-grade B-cell lymphoma", "transformed follicular",
    ],
    "follicular lymphoma": [
        "follicular lymphoma", "FL", "indolent lymphoma",
        "B-NHL", "CD19", "CD20",
        "grade 3B", "transformed FL",
        "axicabtagene ciloleucel",
    ],
    "mantle cell lymphoma": [
        "mantle cell lymphoma", "MCL",
        "CD19", "brexucabtagene autoleucel", "Tecartus",
        "BTK inhibitor", "ibrutinib",
    ],
    "multiple myeloma": [
        "multiple myeloma", "MM", "myeloma",
        "relapsed refractory myeloma", "r/r MM",
        "BCMA", "GPRC5D", "CD38", "FcRH5",
        "idecabtagene vicleucel", "ciltacabtagene autoleucel",
        "Abecma", "Carvykti",
        "APRIL", "plasma cell",
    ],
    "hodgkin": [
        "Hodgkin lymphoma", "HL", "classical Hodgkin",
        "CD30", "Reed-Sternberg",
        "brentuximab vedotin", "nivolumab",
    ],
    "aml": [
        "AML", "acute myeloid leukemia",
        "CD33", "CD123", "CLL-1", "FLT3", "CLEC12A",
        "myeloid antigen", "myeloid leukemia",
        "lineage switch", "mixed phenotype",
    ],
    "cll": [
        "CLL", "chronic lymphocytic leukemia",
        "CD19", "ROR1",
        "ibrutinib", "venetoclax", "Richter transformation",
    ],

    # --- Solid tumors ---
    "neuroblastoma": [
        "neuroblastoma", "GD2", "B7-H3",
        "pediatric solid tumor", "MYCN amplification",
        "dinutuximab", "high-risk neuroblastoma",
    ],
    "glioblastoma": [
        "glioblastoma", "GBM", "glioma",
        "EGFRvIII", "IL13Ralpha2", "HER2", "GD2",
        "intracranial", "blood-brain barrier",
        "tumor microenvironment", "immunosuppressive",
    ],
    "mesothelioma": [
        "mesothelioma", "malignant pleural mesothelioma",
        "mesothelin", "MSLN",
        "asbestos", "pleural",
    ],
    "pancreatic cancer": [
        "pancreatic cancer", "PDAC",
        "mesothelin", "Claudin 18.2", "MUC1", "HER2",
        "pancreatic ductal adenocarcinoma",
        "desmoplastic stroma", "immunosuppressive TME",
    ],
    "liver cancer": [
        "liver cancer", "hepatocellular carcinoma", "HCC",
        "GPC3", "glypican-3", "AFP",
    ],
    "prostate cancer": [
        "prostate cancer", "PSMA", "mCRPC",
        "castration-resistant", "prostate-specific membrane antigen",
    ],
    "ovarian cancer": [
        "ovarian cancer", "mesothelin", "MUC16",
        "HER2", "folate receptor alpha", "FRa",
    ],

    # --- Additional disease expansions ---
    "bpdcn": [
        "BPDCN", "blastic plasmacytoid dendritic cell neoplasm",
        "CD123 positive neoplasm",
    ],
    "mds": [
        "MDS", "myelodysplastic syndrome", "myelodysplasia",
        "pre-leukemia", "refractory anemia",
    ],
    "dipg": [
        "DIPG", "diffuse intrinsic pontine glioma", "DMG",
        "diffuse midline glioma", "H3K27M mutant",
        "pediatric brainstem glioma",
    ],
    "autoimmune_cart": [
        "autoimmune CAR-T", "CAR-T for lupus", "CAR-T for myasthenia gravis",
        "CAR-T for scleroderma", "B-cell depletion autoimmune",
        "anti-CD19 autoimmune",
    ],
    "t_cell_malignancy": [
        "T-ALL", "T-cell ALL", "T-cell lymphoma", "PTCL", "AITL", "ATLL",
        "Sézary syndrome", "mycosis fungoides", "CD7 positive lymphoma",
    ],
    "cns_lymphoma": [
        "CNS lymphoma", "primary CNS lymphoma", "PCNSL",
        "secondary CNS lymphoma", "leptomeningeal lymphoma",
        "intrathecal CAR-T",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 3. TOXICITY_EXPANSION
#    Maps toxicity / safety keywords to related terms
# ═══════════════════════════════════════════════════════════════════════

TOXICITY_EXPANSION: Dict[str, List[str]] = {
    # --- Cytokine release syndrome ---
    "crs": [
        "CRS", "cytokine release syndrome", "cytokine storm",
        "tocilizumab", "Actemra", "siltuximab",
        "IL-6", "IL-6 receptor", "sIL-6R",
        "ferritin", "CRP", "C-reactive protein",
        "fever", "hypotension", "hypoxia",
        "Lee grading", "ASTCT grading",
        "grade 3 CRS", "grade 4 CRS",
        "vasopressors", "dexamethasone",
    ],
    "cytokine release": [
        "CRS", "cytokine release syndrome", "cytokine storm",
        "tocilizumab", "IL-6", "ferritin",
    ],

    # --- Neurotoxicity / ICANS ---
    "icans": [
        "ICANS", "immune effector cell-associated neurotoxicity syndrome",
        "neurotoxicity", "CRES",
        "ICE score", "encephalopathy",
        "cerebral edema", "aphasia", "tremor", "seizure",
        "dexamethasone", "corticosteroids",
        "endothelial activation", "BBB disruption",
        "blood-brain barrier",
    ],
    "neurotoxicity": [
        "ICANS", "neurotoxicity", "CRES",
        "immune effector cell-associated neurotoxicity syndrome",
        "ICE score", "encephalopathy",
        "cerebral edema", "dexamethasone",
    ],

    # --- B-cell aplasia ---
    "b-cell aplasia": [
        "B-cell aplasia", "hypogammaglobulinemia",
        "immunoglobulin replacement", "IVIG",
        "agammaglobulinemia", "B-cell recovery",
        "on-target off-tumor", "CD19 depletion",
        "infection risk", "humoral immunity",
    ],
    "hypogammaglobulinemia": [
        "hypogammaglobulinemia", "B-cell aplasia",
        "IVIG", "immunoglobulin replacement",
        "infection risk",
    ],

    # --- HLH / MAS ---
    "hlh": [
        "HLH", "hemophagocytic lymphohistiocytosis",
        "macrophage activation syndrome", "MAS",
        "ferritin", "triglycerides", "fibrinogen",
        "soluble IL-2 receptor", "sCD25",
        "NK cell dysfunction", "etoposide", "ruxolitinib",
        "anakinra",
    ],
    "mas": [
        "MAS", "macrophage activation syndrome",
        "HLH", "hemophagocytic lymphohistiocytosis",
        "ferritin", "hyperinflammation",
    ],

    # --- Tumor lysis syndrome ---
    "tumor lysis": [
        "tumor lysis syndrome", "TLS",
        "hyperuricemia", "hyperkalemia", "hyperphosphatemia",
        "hypocalcemia", "rasburicase", "allopurinol",
        "high tumor burden", "renal injury",
    ],
    "tls": [
        "tumor lysis syndrome", "TLS",
        "hyperuricemia", "rasburicase", "high tumor burden",
    ],

    # --- Graft-versus-host disease ---
    "gvhd": [
        "GvHD", "graft-versus-host disease",
        "allogeneic", "donor-derived",
        "TRAC knockout", "TCR disruption",
        "CRISPR", "universal CAR-T",
        "alemtuzumab", "cyclophosphamide",
        "skin", "liver", "GI tract",
    ],
    "graft versus host": [
        "GvHD", "graft-versus-host disease",
        "allogeneic", "TRAC knockout", "TCR disruption",
    ],

    # --- Cytopenias ---
    "cytopenia": [
        "cytopenia", "pancytopenia", "neutropenia",
        "thrombocytopenia", "prolonged cytopenia",
        "bone marrow suppression", "GCSF",
        "hematopoietic recovery", "delayed recovery",
    ],

    # --- On-target off-tumor ---
    "on-target off-tumor": [
        "on-target off-tumor", "B-cell aplasia",
        "normal tissue toxicity", "antigen expression",
        "safety switch", "suicide gene",
        "iCasp9", "inducible caspase 9",
        "truncated EGFR", "CD20 safety switch",
    ],

    # --- Coagulopathy ---
    "coagulopathy": [
        "DIC", "disseminated intravascular coagulation", "fibrinogen consumption",
        "coagulopathy", "D-dimer elevation", "coagulation cascade", "ISTH DIC score",
    ],

    # --- Cardiac toxicity ---
    "cardiac toxicity": [
        "cardiac toxicity", "troponin elevation", "cardiomyopathy",
        "LV dysfunction", "cardiac arrest", "myocardial injury",
        "BNP elevation", "arrhythmia", "cardio-oncology",
    ],

    # --- Renal toxicity ---
    "renal toxicity": [
        "renal toxicity", "acute kidney injury", "AKI",
        "creatinine elevation", "dialysis", "TLS nephropathy",
        "uric acid nephropathy", "electrolyte derangement",
    ],

    # --- Secondary malignancy ---
    "secondary malignancy": [
        "secondary malignancy", "T-cell lymphoma", "insertional mutagenesis",
        "insertional oncogenesis", "clonal expansion", "FDA boxed warning",
        "long-term follow-up", "integration site analysis",
    ],

    # --- Additional toxicity expansions ---
    "immune_reconstitution": [
        "immune reconstitution", "immunodeficiency", "hypogammaglobulinemia",
        "IVIG replacement", "vaccination post CAR-T", "immune recovery",
        "IgG levels",
    ],
    "infection_post_cart": [
        "infection", "opportunistic infection", "PJP prophylaxis",
        "CMV reactivation", "EBV reactivation", "fungal infection",
        "bacterial sepsis post CAR-T", "COVID CAR-T",
    ],
    "hepatotoxicity": [
        "hepatotoxicity", "liver injury", "AST ALT elevation",
        "hepatic toxicity", "liver failure CAR-T", "transaminitis",
    ],
    "macrophage_activation": [
        "macrophage activation syndrome", "MAS",
        "hemophagocytic lymphohistiocytosis", "HLH MAS overlap",
        "ferritin spike", "secondary HLH", "hyperferritinemia",
    ],
    "dermatologic": [
        "skin toxicity", "rash", "nail toxicity", "dermatologic adverse event",
        "GPRC5D skin", "alopecia", "Stevens-Johnson", "skin GVHD",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 4. MANUFACTURING_EXPANSION
#    Maps manufacturing / CMC terms to related concepts
# ═══════════════════════════════════════════════════════════════════════

MANUFACTURING_EXPANSION: Dict[str, List[str]] = {
    # --- Transduction ---
    "transduction": [
        "transduction efficiency", "lentiviral vector", "LVV",
        "gamma-retroviral", "retroviral",
        "MOI", "multiplicity of infection",
        "VCN", "vector copy number",
        "viral titer", "functional titer",
        "p24", "FACS transduction",
        "RetroNectin", "polybrene",
        "spinoculation", "static transduction",
    ],
    "lentiviral": [
        "lentiviral vector", "LVV", "lentivirus",
        "HIV-1 derived", "self-inactivating", "SIN vector",
        "VSV-G pseudotyping", "third-generation packaging",
        "transduction efficiency", "VCN",
        "insertional mutagenesis",
    ],
    "viral vector": [
        "lentiviral vector", "gamma-retroviral vector",
        "AAV", "adeno-associated virus",
        "viral titer", "transduction",
        "VSV-G", "pseudotyping",
    ],
    "non-viral": [
        "non-viral delivery", "electroporation",
        "transposon", "Sleeping Beauty", "piggyBac",
        "mRNA electroporation", "lipid nanoparticle", "LNP",
        "CRISPR knock-in", "HDR",
        "transient expression",
    ],

    # --- T-cell expansion ---
    "expansion": [
        "T-cell expansion", "ex vivo expansion",
        "anti-CD3/CD28 beads", "Dynabeads",
        "TransAct", "T Cell TransAct",
        "Expander beads",
        "IL-2", "IL-7", "IL-15", "IL-21",
        "fold expansion", "doubling time",
        "population doublings",
        "G-Rex", "WAVE bioreactor", "Prodigy",
        "CliniMACS Prodigy",
        "fed-batch", "perfusion",
    ],
    "bioreactor": [
        "bioreactor", "G-Rex", "WAVE bioreactor",
        "CliniMACS Prodigy", "Cocoon",
        "Lonza Cocoon", "Miltenyi Prodigy",
        "automated manufacturing", "closed system",
        "fed-batch", "rocking motion",
    ],

    # --- Leukapheresis ---
    "leukapheresis": [
        "leukapheresis", "apheresis",
        "starting material", "PBMC",
        "peripheral blood mononuclear cells",
        "T-cell enrichment", "CD4/CD8 selection",
        "Ficoll", "elutriation",
        "Spectra Optia", "lymphodepletion timing",
    ],
    "apheresis": [
        "apheresis", "leukapheresis",
        "starting material", "PBMC", "T-cell enrichment",
    ],

    # --- Cryopreservation ---
    "cryopreservation": [
        "cryopreservation", "cryopreserved product",
        "CryoStor", "CS10", "DMSO",
        "controlled-rate freezer", "CRF",
        "liquid nitrogen", "LN2", "vapor phase",
        "thaw", "post-thaw viability",
        "post-thaw recovery", "cold chain",
        "bedside thaw",
    ],
    "formulation": [
        "formulation", "final product",
        "infusion bag", "cryopreservation medium",
        "excipients", "human serum albumin", "HSA",
        "PlasmaLyte", "dextran",
    ],

    # --- Release testing ---
    "release testing": [
        "release testing", "certificate of analysis", "CoA",
        "lot release", "quality control", "QC",
        "sterility", "endotoxin", "LAL",
        "mycoplasma", "RCL", "replication-competent lentivirus",
        "identity", "potency", "purity",
        "CAR expression", "percent CAR+",
        "viability", "cell count", "total viable cells",
        "VCN release", "in-process testing",
    ],
    "potency": [
        "potency assay", "functional potency",
        "cytotoxicity assay", "cytokine secretion",
        "IFN-gamma", "specific lysis",
        "potency release criterion",
    ],
    "identity": [
        "identity testing", "flow cytometry",
        "CAR detection", "Protein L",
        "anti-idiotype", "CAR+ percentage",
        "CD3+", "CD4/CD8 ratio",
    ],

    # --- Vein-to-vein time ---
    "vein-to-vein": [
        "vein-to-vein time", "turnaround time",
        "manufacturing slot", "scheduling",
        "bridging therapy", "patient waitlist",
        "out-of-spec", "manufacturing failure",
        "point-of-care manufacturing",
    ],
    "point-of-care": [
        "point-of-care", "POC manufacturing",
        "decentralized manufacturing",
        "bedside manufacturing", "automated",
        "CliniMACS Prodigy", "Cocoon",
        "rapid manufacturing",
    ],

    # --- Lymphodepletion ---
    "lymphodepletion": [
        "lymphodepletion", "conditioning regimen",
        "fludarabine", "cyclophosphamide",
        "Flu/Cy", "bendamustine",
        "T-cell depletion", "lymphopenia",
        "homeostatic expansion",
    ],

    # --- Transposon / non-viral gene transfer ---
    "transposon": [
        "Sleeping Beauty", "PiggyBac", "transposon", "non-viral gene transfer",
        "DNA transposition", "transposase",
    ],

    # --- mRNA electroporation ---
    "mrna electroporation": [
        "mRNA CAR-T", "mRNA electroporation", "transient CAR expression",
        "mRNA transfection", "repeated dosing", "non-integrating",
    ],

    # --- CRISPR knock-in ---
    "crispr knock-in": [
        "CRISPR knock-in", "TRAC integration", "site-specific integration",
        "Cas9 RNP", "HDR template", "gene editing CAR-T", "Caribou CB-010",
    ],

    # --- iPSC manufacturing ---
    "ipsc manufacturing": [
        "iPSC-derived", "iPSC CAR-T", "induced pluripotent stem cell",
        "FT819", "Fate Therapeutics", "universal donor",
        "master cell bank", "allogeneic iPSC",
    ],

    # --- Automated manufacturing ---
    "automated manufacturing": [
        "CliniMACS Prodigy", "automated manufacturing", "closed-system",
        "Cocoon platform", "Lonza Cocoon", "Miltenyi Prodigy", "Sepax",
        "decentralized manufacturing", "point-of-care manufacturing",
        "bedside manufacturing", "closed system", "point-of-care automated",
    ],

    # --- Additional manufacturing expansions ---
    "vector_production": [
        "viral vector", "lentiviral vector", "retroviral vector",
        "vector titer", "GMP vector", "vector manufacturing",
        "TU/mL", "transduction unit", "vector copy number", "VCN",
        "RCL testing",
    ],
    "quality_control": [
        "quality control", "release testing", "CAR expression",
        "potency assay", "sterility test", "mycoplasma", "endotoxin",
        "USP 71", "BacT/ALERT", "rapid sterility",
    ],
    "supply_chain": [
        "supply chain", "cold chain", "cryoshipper",
        "chain of identity", "chain of custody",
        "vein to vein logistics", "thaw to infuse",
        "scheduling coordination",
    ],
    "cost_economics": [
        "manufacturing cost", "COGS", "cost per dose",
        "reimbursement", "value-based pricing",
        "outcomes-based contract", "indication-based pricing",
    ],
    "ipsc_manufacturing": [
        "iPSC derived", "induced pluripotent stem cell",
        "off-the-shelf CAR-T", "master cell bank", "FT819",
        "universal donor", "iPSC differentiation",
    ],
    "activation_media": [
        "T-cell activation", "IL-7 IL-15", "IL-2 expansion",
        "TexMACS", "OpTmizer", "X-VIVO",
        "anti-CD3 CD28", "TransAct", "Dynabeads",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 5. MECHANISM_EXPANSION
#    Maps mechanism / biology terms to related concepts
# ═══════════════════════════════════════════════════════════════════════

MECHANISM_EXPANSION: Dict[str, List[str]] = {
    # --- Resistance / Antigen escape ---
    "resistance": [
        "antigen loss", "antigen escape", "antigen-negative relapse",
        "lineage switch", "myeloid switch",
        "trogocytosis", "antigen masking", "epitope masking",
        "CD19-negative", "BCMA-negative",
        "bispecific CAR", "dual-targeting",
        "tandem CAR", "OR-gate",
    ],
    "antigen escape": [
        "antigen escape", "antigen loss", "antigen-negative relapse",
        "CD19-negative", "lineage switch",
        "bispecific", "dual-targeting",
        "alternative splicing", "truncated antigen",
    ],
    "antigen loss": [
        "antigen loss", "antigen escape", "antigen-negative relapse",
        "lineage switch", "bispecific", "dual-targeting",
    ],
    "lineage switch": [
        "lineage switch", "myeloid switch", "phenotypic switch",
        "mixed phenotype acute leukemia", "MPAL",
        "antigen escape", "CD19-negative AML",
    ],
    "trogocytosis": [
        "trogocytosis", "antigen transfer", "antigen stripping",
        "CAR-mediated trogocytosis", "fratricide",
        "antigen density reduction",
    ],

    # --- T-cell exhaustion ---
    "exhaustion": [
        "T-cell exhaustion", "CAR-T exhaustion",
        "PD-1", "PDCD1", "LAG-3", "LAG3",
        "TIM-3", "HAVCR2", "TIGIT",
        "TOX", "NR4A", "transcription factor",
        "tonic signaling", "chronic stimulation",
        "epigenetic remodeling", "chromatin accessibility",
        "terminal differentiation",
        "checkpoint blockade", "PD-1 knockout",
        "anti-PD-1", "ipilimumab",
    ],
    "pd-1": [
        "PD-1", "PDCD1", "programmed death-1",
        "checkpoint inhibitor", "pembrolizumab", "nivolumab",
        "PD-L1", "CD274",
        "T-cell exhaustion", "checkpoint blockade",
    ],
    "tonic signaling": [
        "tonic signaling", "ligand-independent signaling",
        "scFv clustering", "antigen-independent activation",
        "T-cell exhaustion", "premature differentiation",
        "FMC63", "4-1BB amelioration",
    ],

    # --- Persistence ---
    "persistence": [
        "T-cell persistence", "CAR-T persistence",
        "long-term remission", "B-cell aplasia duration",
        "memory T cells", "Tcm", "central memory",
        "stem cell memory", "Tscm",
        "4-1BB", "CD28 costimulation",
        "IL-7", "IL-15", "homeostatic cytokines",
        "qPCR transgene", "flow CAR detection",
    ],
    "memory": [
        "memory T cells", "central memory", "Tcm",
        "effector memory", "Tem",
        "stem cell memory", "Tscm",
        "naive T cells", "Tn",
        "T-cell differentiation",
        "CD62L", "CCR7", "CD45RA", "CD45RO",
        "long-lived persistence",
    ],

    # --- Costimulation ---
    "costimulation": [
        "CD28", "4-1BB", "CD137", "TNFRSF9",
        "OX40", "CD134", "ICOS", "CD278",
        "costimulatory domain", "second generation CAR",
        "third generation CAR",
        "CD28 versus 4-1BB", "signaling kinetics",
        "PI3K", "NF-kB", "TRAF",
    ],
    "4-1bb": [
        "4-1BB", "CD137", "TNFRSF9",
        "costimulatory domain", "TRAF1", "TRAF2",
        "NF-kB signaling", "persistence advantage",
        "oxidative metabolism", "central memory",
        "tisagenlecleucel", "Kymriah",
    ],
    "cd28 costimulation": [
        "CD28", "costimulatory domain",
        "PI3K signaling", "rapid expansion",
        "effector function", "glycolytic metabolism",
        "axicabtagene ciloleucel", "Yescarta",
    ],

    # --- Signaling ---
    "signaling": [
        "CD3-zeta", "CD247", "ITAM",
        "immunoreceptor tyrosine-based activation motif",
        "ZAP-70", "LAT", "SLP-76",
        "phosphorylation", "signal transduction",
        "proximal signaling", "distal signaling",
        "PI3K", "AKT", "mTOR", "NF-kB", "NFAT", "AP-1",
    ],
    "cd3-zeta": [
        "CD3-zeta", "CD247", "ITAM",
        "signal 1", "activation domain",
        "first-generation CAR",
    ],

    # --- Trafficking / Infiltration ---
    "trafficking": [
        "T-cell trafficking", "tumor infiltration",
        "homing", "chemokine receptor",
        "CXCR1", "CXCR2", "CXCR4", "CCR2", "CCR4",
        "IL-8 receptor", "tumor homing",
        "extravasation", "adhesion molecules",
        "integrin", "selectin",
    ],
    "tumor microenvironment": [
        "tumor microenvironment", "TME",
        "immunosuppressive", "regulatory T cells", "Tregs",
        "myeloid-derived suppressor cells", "MDSC",
        "tumor-associated macrophages", "TAM",
        "TGF-beta", "IL-10", "PGE2", "IDO",
        "hypoxia", "HIF-1alpha",
        "extracellular matrix", "fibrotic stroma",
        "checkpoint ligands", "PD-L1",
    ],
    "tme": [
        "tumor microenvironment", "TME",
        "immunosuppressive", "Tregs", "MDSC",
        "TGF-beta", "hypoxia",
    ],

    # --- Cytokine biology ---
    "cytokine": [
        "cytokine", "IFN-gamma", "TNF-alpha", "IL-2",
        "IL-6", "GM-CSF", "IL-10", "IL-17",
        "granzyme B", "perforin",
        "cytolytic activity", "degranulation",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 6. CONSTRUCT_EXPANSION
#    Maps CAR design / engineering terms to related concepts
# ═══════════════════════════════════════════════════════════════════════

CONSTRUCT_EXPANSION: Dict[str, List[str]] = {
    # --- scFv / binding domain ---
    "scfv": [
        "scFv", "single-chain variable fragment",
        "VH", "VL", "variable heavy", "variable light",
        "linker", "G4S linker",
        "humanized", "fully human",
        "VHH", "nanobody", "camelid",
        "binding affinity", "Kd",
        "FMC63", "m971",
    ],
    "nanobody": [
        "nanobody", "VHH", "single-domain antibody",
        "camelid", "llama-derived",
        "small format", "biparatopic",
    ],

    # --- Hinge / spacer ---
    "hinge": [
        "hinge", "spacer", "extracellular spacer",
        "IgG1 hinge", "IgG4 hinge",
        "CD8-alpha hinge", "CD28 hinge",
        "spacer length", "immunological synapse distance",
        "Fc receptor binding", "FcgammaR",
    ],

    # --- Transmembrane domain ---
    "transmembrane": [
        "transmembrane domain", "TM domain",
        "CD8-alpha TM", "CD28 TM",
        "membrane anchoring", "dimerization",
        "signal transduction",
    ],

    # --- Armored / 4th gen ---
    "armored": [
        "armored CAR", "fourth generation", "4th generation",
        "TRUCK", "T cells redirected for universal cytokine killing",
        "IL-12 secreting", "IL-15 secreting", "IL-18 secreting",
        "IL-21 armored",
        "dominant-negative TGF-beta receptor", "dnTGFbRII",
        "PD-1 dominant negative", "switch receptor",
        "constitutive cytokine",
    ],
    "4th generation": [
        "fourth generation CAR", "4th generation",
        "armored CAR", "TRUCK",
        "cytokine-secreting CAR",
        "IL-12", "IL-15", "IL-18",
    ],
    "truck": [
        "TRUCK", "T cells redirected for universal cytokine killing",
        "armored CAR", "IL-12 secreting",
        "cytokine payload", "transgenic cytokine",
    ],

    # --- Bispecific / multi-target ---
    "bispecific": [
        "bispecific CAR", "dual-targeting CAR",
        "tandem CAR", "bivalent CAR",
        "OR-gate logic", "AND-gate logic",
        "CD19/CD22", "BCMA/CD38", "BCMA/GPRC5D",
        "antigen escape prevention",
        "loop CAR", "split CAR",
    ],
    "tandem": [
        "tandem CAR", "bispecific", "bivalent",
        "OR-gate", "two binding domains",
        "CD19/CD22 tandem",
    ],
    "logic gate": [
        "logic gate", "AND-gate", "OR-gate",
        "NOT-gate", "synNotch",
        "synthetic biology", "Boolean logic",
        "split CAR", "if-then circuit",
        "on-switch", "safety circuit",
    ],

    # --- Universal / allogeneic ---
    "universal": [
        "universal CAR-T", "allogeneic CAR-T",
        "off-the-shelf", "donor-derived",
        "TRAC knockout", "TCR knockout",
        "beta-2-microglobulin knockout", "B2M knockout",
        "HLA-negative", "MHC class I knockout",
        "CRISPR-Cas9", "TALEN", "zinc finger nuclease",
        "gene editing", "base editing",
        "NK cell evasion", "HLA-E overexpression",
        "CD52 knockout", "alemtuzumab-resistant",
        "UCART19", "ALLO-501",
    ],
    "allogeneic": [
        "allogeneic", "off-the-shelf", "donor-derived",
        "universal CAR-T",
        "TRAC knockout", "B2M knockout", "HLA knockout",
        "GvHD prevention", "gene editing",
        "CRISPR", "TALEN",
    ],
    "gene editing": [
        "gene editing", "CRISPR-Cas9", "CRISPR",
        "TALEN", "zinc finger nuclease", "ZFN",
        "base editing", "prime editing",
        "TRAC knockout", "B2M knockout", "PD-1 knockout",
        "homology-directed repair", "HDR",
        "NHEJ", "guide RNA", "sgRNA",
        "off-target effects", "translocation",
    ],

    # --- Safety switches ---
    "safety switch": [
        "safety switch", "kill switch", "suicide gene",
        "iCasp9", "inducible caspase 9",
        "AP1903", "rimiducid",
        "truncated EGFR", "EGFRt", "cetuximab",
        "CD20 safety switch", "rituximab",
        "HSV-TK", "herpes simplex thymidine kinase",
        "ganciclovir",
    ],
    "icasp9": [
        "iCasp9", "inducible caspase 9",
        "AP1903", "rimiducid", "dimerizer",
        "safety switch", "suicide gene",
    ],

    # --- CAR generations overview ---
    "first generation": [
        "first generation CAR", "1st generation",
        "CD3-zeta only", "no costimulation",
        "limited persistence", "historical",
    ],
    "second generation": [
        "second generation CAR", "2nd generation",
        "CD28 or 4-1BB costimulation",
        "CD3-zeta plus costimulatory",
        "FDA-approved products",
        "tisagenlecleucel", "axicabtagene ciloleucel",
    ],
    "third generation": [
        "third generation CAR", "3rd generation",
        "dual costimulatory domains",
        "CD28 plus 4-1BB", "CD28 plus OX40",
        "enhanced signaling",
    ],

    # --- Next-gen platforms ---
    "synnotch": [
        "synNotch", "synthetic Notch",
        "AND-gate", "logic-gated",
        "priming receptor", "response element",
        "tissue-specific", "conditional activation",
    ],
    "adapter car": [
        "adapter CAR", "switchable CAR",
        "FITC adapter", "leucine zipper",
        "SpyCatcher/SpyTag", "universal adapter",
        "dose-titratable", "modular targeting",
    ],

    # --- Additional construct expansions ---
    "fratricide_resistant": [
        "fratricide resistant", "CD7 knockout", "CD5 knockout",
        "T-cell antigen targeting", "self-antigen CAR",
        "gene editing fratricide",
    ],
    "multi_antigen": [
        "multi-antigen", "sequential targeting", "combinatorial antigen",
        "AND gate", "OR gate", "Boolean logic CAR", "multi-specific",
    ],
    "car_nk": [
        "CAR-NK", "natural killer cell", "NK cell therapy",
        "off-the-shelf NK", "cord blood NK", "iPSC-NK",
        "FT596", "NK-92",
    ],
    "tcr_therapy": [
        "TCR therapy", "TCR-T", "T-cell receptor",
        "NY-ESO-1", "MAGE-A4", "TCR-redirected",
        "peptide-MHC", "engineered TCR",
    ],
    "next_gen_costim": [
        "next generation costimulation", "CD28 4-1BB hybrid",
        "ICOS", "OX40", "CD27", "MyD88/CD40",
        "toll-like receptor signaling",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 7. SAFETY_EXPANSION
# ═══════════════════════════════════════════════════════════════════════

SAFETY_EXPANSION: Dict[str, List[str]] = {
    "adverse event": [
        "adverse event", "AE", "safety signal", "toxicity", "side effect",
        "post-market", "pharmacovigilance", "FAERS", "label update",
    ],
    "pharmacovigilance": [
        "pharmacovigilance", "post-market surveillance", "FAERS", "safety monitoring",
        "adverse event reporting", "label change", "boxed warning", "REMS",
    ],
    "rems": [
        "REMS", "risk evaluation", "mitigation strategy", "certification",
        "treatment center", "safety protocol", "CRS management", "ICANS management",
    ],
    "secondary malignancy": [
        "secondary malignancy", "T-cell lymphoma", "MDS", "AML",
        "insertional mutagenesis", "oncogenesis", "long-term safety",
        "FDA boxed warning", "15-year follow-up",
    ],
    "long-term safety": [
        "long-term safety", "prolonged cytopenia", "B-cell aplasia",
        "hypogammaglobulinemia", "IVIG", "immunoglobulin replacement",
        "late-onset", "delayed toxicity", "chronic GVHD",
    ],
    "cytopenia": [
        "cytopenia", "neutropenia", "thrombocytopenia", "anemia",
        "pancytopenia", "prolonged cytopenia", "bone marrow suppression",
        "G-CSF", "platelet transfusion",
    ],
    "infection": [
        "infection", "opportunistic infection", "bacterial", "fungal",
        "viral reactivation", "CMV", "HHV-6", "aspergillus",
        "hypogammaglobulinemia", "immunodeficiency",
    ],
    "cardiac": [
        "cardiac toxicity", "cardiomyopathy", "arrhythmia", "troponin",
        "heart failure", "cardiac arrest", "myocarditis",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 8. BIOMARKER_EXPANSION
# ═══════════════════════════════════════════════════════════════════════

BIOMARKER_EXPANSION: Dict[str, List[str]] = {
    "biomarker": [
        "biomarker", "predictive marker", "prognostic marker",
        "pharmacodynamic", "response prediction", "outcome prediction",
    ],
    "ferritin": [
        "ferritin", "serum ferritin", "CRS prediction",
        "inflammatory marker", "iron storage", "hyperferritinemia",
    ],
    "crp": [
        "CRP", "C-reactive protein", "inflammatory biomarker",
        "acute phase", "CRS monitoring",
    ],
    "il-6": [
        "IL-6", "interleukin-6", "cytokine storm", "tocilizumab target",
        "CRS biomarker", "inflammatory cytokine",
    ],
    "mrd": [
        "MRD", "minimal residual disease", "flow cytometry MRD",
        "PCR MRD", "measurable residual disease", "MRD negative",
        "deep response", "molecular remission",
        "flow MRD", "NGS MRD",
        "10^-4 sensitivity", "10^-5 sensitivity", "10^-6 sensitivity",
    ],
    "car expansion": [
        "CAR-T expansion", "CAR expansion", "peak expansion", "Cmax", "Tmax", "AUC",
        "transgene copies", "persistence",
        "CAR copies per microgram DNA", "expansion kinetics",
        "qPCR", "in vivo proliferation", "pharmacokinetics",
    ],
    "exhaustion marker": [
        "exhaustion", "PD-1", "LAG-3", "TIM-3", "TOX", "NR4A",
        "T-cell dysfunction", "checkpoint", "terminal differentiation",
        "epigenetic exhaustion",
    ],
    "ldh": [
        "LDH", "lactate dehydrogenase", "tumor burden marker",
        "prognostic factor", "metabolic marker",
    ],
    "ctdna": [
        "ctDNA", "circulating tumor DNA", "liquid biopsy",
        "cell-free DNA", "molecular response", "genomic profiling",
        "cfDNA dynamics",
    ],
    "sbcma": [
        "sBCMA", "soluble BCMA", "BCMA shedding",
        "gamma-secretase", "decoy antigen", "BCMA resistance",
    ],
    "tcm": [
        "Tcm", "central memory", "CD45RA-", "CCR7+", "CD62L+",
        "T-cell fitness", "memory phenotype", "naive T-cell",
    ],
    "cd4 cd8": [
        "CD4:CD8", "CD4/CD8 ratio", "T-cell composition",
        "defined composition", "product phenotype",
    ],

    # --- TOX exhaustion ---
    "tox exhaustion": [
        "TOX", "TOX transcription factor", "thymocyte selection",
        "TOX exhaustion", "TOX+", "terminal exhaustion marker",
    ],

    # --- NR4A ---
    "nr4a": [
        "NR4A", "NR4A1", "NR4A2", "NR4A3", "NR4A family",
        "NR4A knockout", "exhaustion transcription factor", "Nur77",
    ],

    # --- Antigen density ---
    "antigen density": [
        "antigen density", "target expression level", "QIFIKIT",
        "molecules per cell", "antigen-low escape", "antigen threshold",
    ],

    # --- Tumor burden ---
    "tumor burden": [
        "tumor burden", "metabolic tumor volume", "MTV", "bulky disease",
        "high tumor burden", "baseline disease burden", "marrow blast percentage",
    ],

    # --- D-dimer ---
    "d-dimer": [
        "D-dimer", "fibrin degradation products", "FDP",
        "coagulation marker", "DIC marker", "fibrinolysis",
    ],

    # --- Troponin ---
    "troponin": [
        "troponin", "hs-cTnI", "hs-cTnT", "cardiac troponin",
        "myocardial injury marker", "cardiac biomarker",
    ],

    # --- Additional biomarker expansions ---
    "t_cell_phenotype": [
        "T-cell phenotype", "central memory", "Tcm",
        "stem cell memory", "Tscm", "effector memory", "Tem",
        "terminally differentiated", "Ttd", "naive T-cell", "Tn",
    ],
    "cytokine_profile": [
        "cytokine profile", "cytokine storm", "IL-6 level",
        "IFN-gamma", "IL-1 beta", "TNF-alpha", "IL-10",
        "IL-2 receptor", "cytokine kinetics",
    ],
    "ang2_ang1": [
        "angiopoietin", "Ang-2", "Ang-1", "Ang-2/Ang-1 ratio",
        "endothelial activation", "ICANS biomarker",
        "vascular permeability",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 9. REGULATORY_EXPANSION
# ═══════════════════════════════════════════════════════════════════════

REGULATORY_EXPANSION: Dict[str, List[str]] = {
    "fda": [
        "FDA", "Food and Drug Administration", "regulatory",
        "BLA", "biologics license", "approval", "label",
    ],
    "bla": [
        "BLA", "biologics license application", "regulatory submission",
        "NDA", "marketing authorization", "approval pathway",
    ],
    "breakthrough": [
        "breakthrough therapy", "BTD", "expedited program",
        "accelerated development", "FDA designation",
    ],
    "rmat": [
        "RMAT", "regenerative medicine advanced therapy",
        "cell therapy designation", "expedited approval",
    ],
    "accelerated approval": [
        "accelerated approval", "surrogate endpoint",
        "confirmatory trial", "conditional approval", "full approval",
    ],
    "ema": [
        "EMA", "European Medicines Agency", "CHMP",
        "marketing authorization", "conditional approval", "EU approval",
        "conditional marketing authorization", "PRIME designation",
        "ATMP classification", "advanced therapy",
    ],
    "label update": [
        "label update", "prescribing information", "boxed warning",
        "black box warning", "safety communication", "Dear Doctor letter",
    ],
    "post-marketing": [
        "post-marketing", "Phase 4", "post-approval",
        "registry study", "long-term follow-up", "PMR", "PMC",
    ],

    # --- Additional regulatory expansions ---
    "pmda": [
        "PMDA", "Japan regulatory", "Japanese approval",
        "Pharmaceuticals and Medical Devices Agency",
    ],
    "health_canada": [
        "Health Canada", "Canadian approval", "NOC",
        "Notice of Compliance",
    ],
    "ema_specifics": [
        "EMA", "European Medicines Agency", "CHMP",
        "conditional marketing authorization", "PRIME designation",
        "ATMP classification", "advanced therapy",
    ],
    "secondary_malignancy_regulatory": [
        "secondary malignancy", "T-cell lymphoma", "boxed warning",
        "FDA black box", "class-wide label", "REMS update",
        "15-year follow-up", "insertional mutagenesis",
    ],
    "real_world_regulatory": [
        "real world evidence regulatory", "post-marketing surveillance",
        "PASS", "PAES", "registry mandate",
        "long-term follow-up study", "15-year LTFU",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 10. SEQUENCE_EXPANSION
# ═══════════════════════════════════════════════════════════════════════

SEQUENCE_EXPANSION: Dict[str, List[str]] = {
    "scfv": [
        "scFv", "single-chain variable fragment", "antibody fragment",
        "VH", "VL", "linker", "binding domain",
    ],
    "cdr": [
        "CDR", "complementarity determining region", "CDR3",
        "hypervariable region", "antigen binding loop", "paratope",
    ],
    "binding affinity": [
        "binding affinity", "Kd", "dissociation constant",
        "kon", "koff", "SPR", "BLI", "affinity maturation",
    ],
    "humanization": [
        "humanization", "humanized antibody", "CDR grafting",
        "framework", "deimmunization", "anti-drug antibody", "ADA",
        "immunogenicity",
    ],
    "fmc63": [
        "FMC63", "anti-CD19 scFv", "murine origin",
        "tisagenlecleucel binder", "axicabtagene binder",
    ],
    "nanobody": [
        "nanobody", "VHH", "single-domain antibody", "camelid",
        "llama antibody", "sdAb", "LCAR-B38M",
    ],
    "darpin": [
        "DARPin", "designed ankyrin repeat", "non-antibody scaffold",
        "Centyrin", "fibronectin", "alternative binding domain",
    ],
    "bispecific car": [
        "bispecific CAR", "tandem CAR", "dual-targeting",
        "bicistronic", "loop CAR", "CD19/CD22", "OR-gate",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 11. REALWORLD_EXPANSION
# ═══════════════════════════════════════════════════════════════════════

REALWORLD_EXPANSION: Dict[str, List[str]] = {
    "real-world": [
        "real-world", "RWE", "real-world evidence", "real-world data",
        "clinical practice", "commercial experience", "post-approval",
    ],
    "cibmtr": [
        "CIBMTR", "Center for International Blood and Marrow Transplant Research",
        "transplant registry", "national registry",
    ],
    "registry": [
        "registry", "CIBMTR", "EBMT", "DESCAR-T", "observational study",
        "national registry", "multi-center registry",
    ],
    "community": [
        "community center", "community practice", "non-academic",
        "community oncology", "access to care", "referral pattern",
    ],
    "academic": [
        "academic center", "academic medical center", "tertiary center",
        "CAR-T center of excellence", "FACT accredited",
    ],
    "elderly": [
        "elderly", "older adults", "geriatric", "age ≥65",
        "age ≥70", "frailty", "comorbidities", "fitness",
    ],
    "bridging therapy": [
        "bridging therapy", "bridging chemotherapy", "bridging radiation",
        "pre-CAR-T", "disease control", "tumor debulking",
    ],
    "disparities": [
        "disparities", "health disparities", "racial disparities", "ethnic disparities",
        "racial", "ethnic", "socioeconomic",
        "insurance barriers", "financial toxicity", "access equity",
        "underrepresented populations",
        "access", "underserved", "minority", "equity",
    ],
    "resource utilization": [
        "resource utilization", "ICU admission", "readmission",
        "length of stay", "cost", "healthcare economics", "HCRU",
    ],
    "long-term follow-up": [
        "long-term follow-up", "durability", "late relapse",
        "5-year", "3-year", "sustained remission", "cure",
    ],

    # --- Additional real-world expansions ---
    "outpatient": [
        "outpatient CAR-T", "ambulatory", "home monitoring",
        "remote monitoring", "telehealth post CAR-T", "early discharge",
    ],
    "cost_rwe": [
        "total cost of care", "hospital cost", "ICU cost",
        "readmission cost", "QALY", "cost-effectiveness",
        "ICER", "budget impact", "payer perspective",
    ],
    "caregiver": [
        "caregiver", "caregiver burden", "quality of life",
        "patient-reported outcomes", "PRO", "FACT-Lym",
        "EQ-5D", "long-term survivorship",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 12. IMMUNOGENICITY_EXPANSION
#     Maps HLA, ADA, and immunogenicity terms to related concepts
# ═══════════════════════════════════════════════════════════════════════

IMMUNOGENICITY_EXPANSION: Dict[str, List[str]] = {
    "hla": [
        "HLA", "human leukocyte antigen", "MHC", "major histocompatibility complex",
        "HLA-DRB1", "HLA-A*02:01", "antigen presentation", "T-cell epitope",
    ],
    "immunogenicity": [
        "immunogenicity", "immunogenic", "anti-drug antibody", "ADA",
        "neutralizing antibody", "binding antibody", "titer",
        "pre-existing immunity", "HAMA",
    ],
    "ada": [
        "anti-drug antibody", "ADA", "neutralizing antibody", "NAb",
        "binding antibody", "immunogenicity testing", "tiered assay",
        "drug tolerance", "ADA incidence",
    ],
    "humanization": [
        "humanization", "humanized", "CDR grafting", "framework selection",
        "back-mutation", "VH3-23", "VK1-39", "deimmunization",
        "fully human", "phage display",
    ],
    "deimmunization": [
        "deimmunization", "deimmunized", "T-cell epitope removal",
        "EpiMatrix", "NetMHCIIpan", "epitope prediction",
        "framework shuffling", "germline humanization",
    ],
    "anti-drug antibody": [
        "anti-drug antibody", "ADA", "HAMA", "anti-murine antibody",
        "immunogenicity", "CAR-T persistence", "neutralizing",
        "infusion reaction", "anaphylaxis",
    ],
    "cdr grafting": [
        "CDR grafting", "complementarity-determining region", "framework region",
        "humanization", "VH germline", "VL germline", "Kabat numbering",
    ],
    "t-cell epitope": [
        "T-cell epitope", "MHC-II", "HLA-DRB1", "CD4 T helper",
        "ELISpot", "IFN-gamma", "DC-T cell assay",
        "in silico prediction", "NetMHCIIpan",
    ],
    "elispot": [
        "ELISpot", "enzyme-linked immunospot", "IFN-gamma spot",
        "T-cell response", "immunogenicity testing", "in vitro assay",
    ],
    "framework shuffling": [
        "framework shuffling", "framework selection", "germline framework",
        "VH3 family", "VK1 family", "CDR loop grafting",
        "thermal stability", "aggregation resistance",
    ],
    "hama": [
        "HAMA", "human anti-mouse antibody", "anti-murine",
        "murine scFv", "FMC63", "pre-existing immunity",
        "cross-reactivity", "clearance",
    ],
    "netmhciipan": [
        "NetMHCIIpan", "MHC-II prediction", "epitope prediction",
        "EpiMatrix", "IEDB", "immunoinformatics", "in silico immunogenicity",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# ALL EXPANSION MAPS — ordered list for iteration
# ═══════════════════════════════════════════════════════════════════════

ALL_EXPANSION_MAPS: List[tuple] = [
    ("Target Antigen", TARGET_ANTIGEN_EXPANSION),
    ("Disease", DISEASE_EXPANSION),
    ("Toxicity", TOXICITY_EXPANSION),
    ("Manufacturing", MANUFACTURING_EXPANSION),
    ("Mechanism", MECHANISM_EXPANSION),
    ("Construct", CONSTRUCT_EXPANSION),
    ("Safety", SAFETY_EXPANSION),
    ("Biomarker", BIOMARKER_EXPANSION),
    ("Regulatory", REGULATORY_EXPANSION),
    ("Sequence", SEQUENCE_EXPANSION),
    ("RealWorld", REALWORLD_EXPANSION),
    ("Immunogenicity", IMMUNOGENICITY_EXPANSION),
]


# ═══════════════════════════════════════════════════════════════════════
# EXPANSION FUNCTION
# ═══════════════════════════════════════════════════════════════════════


def expand_query(query: str) -> List[str]:
    """Extract expansion terms from a user query.

    Scans the query for keywords matching any expansion map,
    returns a deduplicated list of related terms to broaden
    the search across the 5 CAR-T Milvus collections.

    Args:
        query: Raw user question (e.g., "What causes CRS after
               CD19 CAR-T therapy?")

    Returns:
        Deduplicated list of expansion terms.  Empty list if no
        keywords matched.

    Example::

        >>> expand_query("Why do patients relapse with CD19-negative disease?")
        ['CD19', 'B-ALL', 'DLBCL', ..., 'antigen loss', 'antigen escape', ...]
    """
    query_lower = query.lower()
    matched_terms: Set[str] = set()

    for category, mapping in ALL_EXPANSION_MAPS:
        for keyword, terms in mapping.items():
            if keyword in query_lower:
                matched_terms.update(terms)
                logger.info(
                    f"CAR-T expansion [{category}]: '{keyword}' -> "
                    f"{len(terms)} terms"
                )

    result = sorted(matched_terms)
    if result:
        logger.info(
            f"Query expansion produced {len(result)} unique terms "
            f"from query: {query[:80]}..."
        )
    return result


def expand_query_by_category(query: str) -> Dict[str, List[str]]:
    """Like expand_query but returns terms grouped by expansion category.

    Useful when different collections should weight different
    categories (e.g., manufacturing terms are most relevant to
    the cart_manufacturing collection).

    Args:
        query: Raw user question.

    Returns:
        Dict mapping category name to list of matched terms.
        Only categories with at least one match are included.

    Example::

        >>> expand_query_by_category("transduction efficiency for CD19 CAR")
        {
            'Target Antigen': ['CD19', 'B-ALL', ...],
            'Manufacturing': ['transduction efficiency', 'lentiviral vector', ...],
        }
    """
    query_lower = query.lower()
    categories: Dict[str, Set[str]] = {}

    for category, mapping in ALL_EXPANSION_MAPS:
        for keyword, terms in mapping.items():
            if keyword in query_lower:
                if category not in categories:
                    categories[category] = set()
                categories[category].update(terms)
                logger.debug(
                    f"CAR-T expansion [{category}]: '{keyword}' matched"
                )

    # Convert sets to sorted lists
    return {cat: sorted(terms) for cat, terms in categories.items()}


def get_expansion_stats() -> Dict[str, int]:
    """Return the number of keywords and total terms per expansion map.

    Useful for logging / health checks.

    Returns:
        Dict with category names as keys, keyword counts as values.
    """
    stats: Dict[str, int] = {}
    for category, mapping in ALL_EXPANSION_MAPS:
        total_terms = sum(len(v) for v in mapping.values())
        stats[category] = {
            "keywords": len(mapping),
            "total_terms": total_terms,
        }
    return stats
