"""
Oncology-Specific Query Expansion Maps
=======================================

Provides domain-aware query expansion for precision oncology searches.
Maps oncology keywords (cancer types, genes, therapies, biomarkers, pathways,
resistance mechanisms, clinical terms, trials, immunotherapy, surgery/radiation,
toxicity, and genomics) to lists of semantically related terms for improved
RAG retrieval recall.

Author: Adam Jones
Date:   February 2026
"""

from typing import List


# ---------------------------------------------------------------------------
# 1. Cancer type expansions
# ---------------------------------------------------------------------------
CANCER_TYPE_EXPANSIONS = {
    "NSCLC": [
        "non-small cell lung cancer", "lung adenocarcinoma",
        "lung squamous cell", "EGFR-mutant lung", "ALK-positive lung",
    ],
    "SCLC": [
        "small cell lung cancer", "neuroendocrine lung",
        "extensive-stage SCLC", "limited-stage SCLC", "ES-SCLC",
    ],
    "CRC": [
        "colorectal cancer", "colon cancer", "rectal cancer",
        "metastatic colorectal", "mCRC", "left-sided colon",
    ],
    "HCC": [
        "hepatocellular carcinoma", "liver cancer",
        "hepatobiliary", "fibrolamellar HCC", "Child-Pugh",
    ],
    "TNBC": [
        "triple-negative breast cancer", "basal-like breast",
        "ER-negative PR-negative HER2-negative", "BRCA-mutant breast", "TNBC subtype",
    ],
    "AML": [
        "acute myeloid leukemia", "FLT3-mutant AML",
        "NPM1-mutant AML", "IDH-mutant AML", "secondary AML",
    ],
    "melanoma": [
        "cutaneous melanoma", "BRAF-mutant melanoma", "uveal melanoma",
        "acral melanoma", "mucosal melanoma", "metastatic melanoma",
    ],
    "GBM": [
        "glioblastoma multiforme", "high-grade glioma",
        "IDH-wildtype glioblastoma", "MGMT-methylated GBM", "brain tumor",
    ],
    "RCC": [
        "renal cell carcinoma", "clear cell RCC", "papillary RCC",
        "chromophobe RCC", "kidney cancer", "metastatic RCC",
    ],
    "PDAC": [
        "pancreatic ductal adenocarcinoma", "pancreatic cancer",
        "KRAS-mutant pancreatic", "metastatic pancreatic",
        "borderline resectable pancreatic",
    ],
    "cholangiocarcinoma": [
        "bile duct cancer", "intrahepatic cholangiocarcinoma",
        "extrahepatic cholangiocarcinoma", "biliary tract cancer",
        "FGFR2-fusion cholangiocarcinoma", "IDH1-mutant bile duct",
    ],
    "endometrial": [
        "endometrial cancer", "uterine cancer", "endometrial carcinoma",
        "POLE-ultramutated endometrial", "serous endometrial", "MSI-H endometrial",
    ],
    "sarcoma": [
        "soft tissue sarcoma", "bone sarcoma", "leiomyosarcoma",
        "liposarcoma", "undifferentiated pleomorphic sarcoma", "Ewing sarcoma",
    ],
    "thyroid": [
        "thyroid cancer", "papillary thyroid carcinoma",
        "medullary thyroid cancer", "anaplastic thyroid", "differentiated thyroid",
    ],
    "mesothelioma": [
        "malignant pleural mesothelioma", "peritoneal mesothelioma",
        "BAP1-altered mesothelioma", "asbestos-related cancer",
        "mesothelioma immunotherapy",
    ],
    "HNSCC": [
        "head and neck squamous cell carcinoma", "oropharyngeal cancer",
        "HPV-positive HNSCC", "HPV-negative HNSCC",
        "laryngeal cancer", "oral cavity cancer",
    ],
}

# ---------------------------------------------------------------------------
# 2. Gene expansions
# ---------------------------------------------------------------------------
GENE_EXPANSIONS = {
    "EGFR": [
        "epidermal growth factor receptor", "EGFR L858R",
        "EGFR exon 19 deletion", "EGFR T790M", "EGFR C797S",
    ],
    "KRAS": [
        "Kirsten rat sarcoma", "KRAS G12C", "KRAS G12D",
        "KRAS G12V", "KRAS amplification",
    ],
    "BRAF": [
        "B-Raf proto-oncogene", "BRAF V600E", "BRAF V600K",
        "BRAF class I mutation", "BRAF fusion",
    ],
    "ALK": [
        "anaplastic lymphoma kinase", "EML4-ALK",
        "ALK rearrangement", "ALK fusion", "ALK amplification",
    ],
    "TP53": [
        "tumor protein p53", "p53 loss-of-function",
        "TP53 missense mutation", "TP53 deletion", "Li-Fraumeni",
    ],
    "BRCA1": [
        "breast cancer gene 1", "BRCA1 germline mutation",
        "BRCA1 somatic mutation", "homologous recombination deficient",
        "PARP inhibitor sensitivity",
    ],
    "BRCA2": [
        "breast cancer gene 2", "BRCA2 germline mutation",
        "BRCA2 frameshift", "Fanconi anemia pathway",
        "HRD-positive", "BRCA2 reversion mutation",
    ],
    "PIK3CA": [
        "phosphatidylinositol 3-kinase catalytic subunit alpha",
        "PIK3CA H1047R", "PIK3CA E545K", "PI3K pathway",
        "PI3K inhibitor sensitivity",
    ],
    "MET": [
        "mesenchymal-epithelial transition factor", "MET exon 14 skipping",
        "MET amplification", "MET overexpression", "c-MET",
    ],
    "RET": [
        "rearranged during transfection", "RET fusion",
        "RET M918T", "RET inhibitor", "KIF5B-RET",
    ],
    "NRG1": [
        "neuregulin 1", "NRG1 fusion", "CD74-NRG1",
        "HER2-HER3 signaling", "zenocutuzumab target",
    ],
    "TROP2": [
        "trophoblast cell surface antigen 2", "TACSTD2",
        "TROP2 overexpression", "sacituzumab govitecan target",
        "TROP2-directed ADC",
    ],
    "DLL3": [
        "delta-like ligand 3", "DLL3 expression",
        "Notch inhibitory ligand", "tarlatamab target", "SCLC biomarker",
    ],
    "EZH2": [
        "enhancer of zeste homolog 2", "EZH2 gain-of-function",
        "PRC2 catalytic subunit", "H3K27 trimethylation",
        "tazemetostat target",
    ],
    "CDH1": [
        "E-cadherin", "CDH1 germline mutation",
        "hereditary diffuse gastric cancer", "lobular breast cancer gene",
        "CDH1 loss",
    ],
    "MAP2K1": [
        "MEK1", "MAP2K1 mutation", "MEK inhibitor target",
        "MAPK pathway kinase", "selumetinib target",
    ],
}

# ---------------------------------------------------------------------------
# 3. Therapy expansions
# ---------------------------------------------------------------------------
THERAPY_EXPANSIONS = {
    "osimertinib": [
        "Tagrisso", "EGFR TKI", "third-generation EGFR inhibitor",
        "FLAURA", "FLAURA2",
    ],
    "sotorasib": [
        "Lumakras", "KRAS G12C inhibitor", "AMG 510",
        "CodeBreaK 200", "covalent KRAS inhibitor",
    ],
    "pembrolizumab": [
        "Keytruda", "anti-PD-1", "checkpoint inhibitor",
        "KEYNOTE trials", "PD-1 blockade",
    ],
    "trastuzumab": [
        "Herceptin", "anti-HER2", "HER2-targeted therapy",
        "trastuzumab deruxtecan", "T-DXd",
    ],
    "olaparib": [
        "Lynparza", "PARP inhibitor", "PARPi",
        "synthetic lethality", "BRCA-directed therapy",
    ],
    "nivolumab": [
        "Opdivo", "anti-PD-1", "immune checkpoint inhibitor",
        "CheckMate trials", "PD-1 antibody",
    ],
    "dabrafenib": [
        "Tafinlar", "BRAF inhibitor", "BRAF V600E inhibitor",
        "dabrafenib-trametinib combination", "COMBI-d",
    ],
    "venetoclax": [
        "Venclexta", "BCL-2 inhibitor", "BH3 mimetic",
        "venetoclax-obinutuzumab", "VIALE-A",
    ],
    "lorlatinib": [
        "Lorbrena", "third-generation ALK inhibitor",
        "ALK TKI", "CROWN trial", "ALK-positive NSCLC treatment",
    ],
    "encorafenib": [
        "Braftovi", "BRAF inhibitor", "encorafenib-cetuximab",
        "BEACON CRC", "RAF inhibitor",
    ],
    "sacituzumab": [
        "Trodelvy", "TROP2 ADC", "sacituzumab govitecan",
        "ASCENT trial", "SN-38 payload",
    ],
    "enfortumab": [
        "Padcev", "enfortumab vedotin", "Nectin-4 ADC",
        "EV-302", "MMAE payload",
    ],
    "tarlatamab": [
        "Imdelltra", "DLL3 BiTE", "bispecific T-cell engager",
        "DeLLphi-301", "DLL3-CD3 bispecific",
    ],
    "zolbetuximab": [
        "Vyloy", "anti-CLDN18.2", "Claudin 18.2 antibody",
        "SPOTLIGHT trial", "GLOW trial",
    ],
    "tazemetostat": [
        "Tazverik", "EZH2 inhibitor", "PRC2 inhibitor",
        "H3K27me3 inhibitor", "epigenetic therapy",
    ],
    "capivasertib": [
        "Truqap", "AKT inhibitor", "pan-AKT kinase inhibitor",
        "CAPItello-291", "PI3K pathway inhibitor",
    ],
}

# ---------------------------------------------------------------------------
# 4. Biomarker expansions
# ---------------------------------------------------------------------------
BIOMARKER_EXPANSIONS = {
    "MSI-H": [
        "microsatellite instability high", "mismatch repair deficient",
        "dMMR", "MLH1 loss", "immunotherapy biomarker",
    ],
    "TMB": [
        "tumor mutational burden", "TMB-high",
        "mutations per megabase", "hypermutated", "neoantigen load",
    ],
    "PD-L1": [
        "programmed death-ligand 1", "TPS score",
        "CPS score", "PD-L1 IHC 22C3", "PD-L1 expression",
    ],
    "HRD": [
        "homologous recombination deficiency", "HRD score",
        "genomic instability score", "LOH", "BRCA-ness",
    ],
    "HER2": [
        "human epidermal growth factor receptor 2", "ERBB2",
        "HER2-positive", "HER2-low", "HER2 amplification",
    ],
    "ALK rearrangement": [
        "ALK fusion", "EML4-ALK", "ALK FISH-positive",
        "ALK IHC", "ALK break-apart",
    ],
    "NTRK fusion": [
        "neurotrophic receptor tyrosine kinase fusion",
        "pan-TRK IHC", "NTRK1 fusion", "NTRK2 fusion", "NTRK3 fusion",
        "larotrectinib biomarker",
    ],
    "ctDNA": [
        "circulating tumor DNA", "liquid biopsy",
        "cell-free DNA", "cfDNA", "minimal residual disease",
        "molecular residual disease",
    ],
    "FGFR alteration": [
        "fibroblast growth factor receptor", "FGFR2 fusion",
        "FGFR3 mutation", "FGFR amplification",
        "FGFR inhibitor sensitivity",
    ],
    "ROS1": [
        "ROS proto-oncogene 1", "ROS1 fusion",
        "ROS1 rearrangement", "CD74-ROS1", "crizotinib biomarker",
    ],
    "FRalpha": [
        "folate receptor alpha", "FOLR1 expression",
        "mirvetuximab soravtansine biomarker", "FRα IHC",
        "platinum-resistant ovarian biomarker",
    ],
    "Nectin-4": [
        "Nectin-4 expression", "PVRL4", "enfortumab vedotin biomarker",
        "urothelial carcinoma target", "cell adhesion molecule",
    ],
    "Claudin 18.2": [
        "CLDN18.2 expression", "CLDN18.2 IHC",
        "zolbetuximab companion diagnostic", "gastric cancer biomarker",
        "tight junction protein",
    ],
    "GIS": [
        "genomic instability score", "Myriad myChoice",
        "LOH-TAI-LST composite", "HRD score",
        "PARP inhibitor predictive biomarker",
    ],
}

# ---------------------------------------------------------------------------
# 5. Pathway expansions
# ---------------------------------------------------------------------------
PATHWAY_EXPANSIONS = {
    "MAPK": [
        "RAS-RAF-MEK-ERK", "MAPK cascade", "MAPK pathway",
        "ERK signaling", "RAS signaling",
    ],
    "PI3K/AKT/mTOR": [
        "PI3K pathway", "AKT signaling", "mTOR signaling",
        "PTEN loss", "PI3K-AKT-mTOR axis",
    ],
    "Wnt": [
        "Wnt/beta-catenin", "canonical Wnt pathway",
        "APC mutation", "CTNNB1 signaling", "Wnt target genes",
    ],
    "Hedgehog": [
        "Hedgehog signaling pathway", "SHH pathway",
        "SMO inhibitor", "GLI activation", "Patched receptor",
    ],
    "Notch": [
        "Notch signaling", "Notch1 activation",
        "gamma-secretase", "DLL3 ligand", "Notch pathway crosstalk",
    ],
    "JAK/STAT": [
        "Janus kinase", "STAT signaling",
        "JAK2 V617F", "cytokine signaling", "JAK inhibitor",
    ],
    "DNA damage repair": [
        "DDR pathway", "homologous recombination",
        "non-homologous end joining", "base excision repair",
        "mismatch repair", "PARP dependency",
    ],
    "angiogenesis": [
        "VEGF signaling", "VEGFR pathway",
        "HIF-1 alpha", "anti-angiogenic", "tumor vasculature",
    ],
    "cell cycle": [
        "CDK4/6", "RB pathway", "cyclin D1",
        "G1/S checkpoint", "cell cycle arrest",
    ],
    "apoptosis": [
        "programmed cell death", "BCL-2 family",
        "caspase activation", "BH3 mimetic", "intrinsic apoptosis pathway",
    ],
}

# ---------------------------------------------------------------------------
# 6. Resistance expansions
# ---------------------------------------------------------------------------
RESISTANCE_EXPANSIONS = {
    "resistance": [
        "acquired resistance", "secondary mutations",
        "bypass pathway", "treatment failure", "dose escalation",
    ],
    "EGFR resistance": [
        "T790M gatekeeper mutation", "C797S mutation",
        "MET amplification bypass", "HER2 amplification bypass",
        "small cell transformation",
    ],
    "ALK resistance": [
        "ALK solvent-front mutation", "ALK G1202R",
        "ALK compound mutations", "ALK L1196M gatekeeper",
        "ALK resistance cascade",
    ],
    "BRAF resistance": [
        "MAPK reactivation", "NRAS mutation",
        "MEK1 mutation", "BRAF amplification",
        "ERK pathway feedback",
    ],
    "immunotherapy resistance": [
        "immune escape", "antigen presentation loss",
        "beta-2-microglobulin loss", "JAK1/2 mutation",
        "T cell exhaustion", "cold tumor",
    ],
    "chemotherapy resistance": [
        "multidrug resistance", "P-glycoprotein",
        "DNA repair upregulation", "efflux pump",
        "apoptosis evasion",
    ],
    "endocrine resistance": [
        "ESR1 mutation", "ER-independent growth",
        "CDK4/6 inhibitor resistance", "PI3K activation",
        "hormone receptor loss",
    ],
    "PARP inhibitor resistance": [
        "BRCA reversion mutation", "53BP1 loss",
        "HR restoration", "drug efflux", "PARP1 mutation",
    ],
    "anti-HER2 resistance": [
        "HER2 truncation", "p95-HER2",
        "PI3K pathway activation", "HER3 upregulation",
        "trastuzumab resistance",
    ],
    "targeted therapy resistance": [
        "on-target resistance", "off-target resistance",
        "lineage plasticity", "phenotypic switching",
        "tumor heterogeneity",
    ],
    "ADC resistance": [
        "antibody-drug conjugate resistance", "antigen downregulation",
        "payload efflux", "linker cleavage failure",
        "bystander effect loss",
    ],
    "bispecific resistance": [
        "bispecific antibody resistance", "target antigen loss",
        "T cell exhaustion", "cytokine release avoidance",
        "antigen escape",
    ],
    "PI3K inhibitor resistance": [
        "PI3K pathway reactivation", "insulin feedback loop",
        "PTEN restoration failure", "mTOR bypass activation",
        "AKT-independent signaling",
    ],
}

# ---------------------------------------------------------------------------
# 7. Clinical expansions
# ---------------------------------------------------------------------------
CLINICAL_EXPANSIONS = {
    "response": [
        "overall response rate", "ORR",
        "complete response", "partial response", "RECIST criteria",
    ],
    "survival": [
        "overall survival", "OS", "progression-free survival",
        "PFS", "disease-free survival", "DFS",
    ],
    "staging": [
        "TNM staging", "AJCC staging",
        "stage IV", "metastatic", "locally advanced",
    ],
    "performance status": [
        "ECOG PS", "Karnofsky performance",
        "ECOG 0-1", "functional status", "patient fitness",
    ],
    "progression": [
        "disease progression", "progressive disease",
        "RECIST progression", "radiographic progression",
        "clinical deterioration",
    ],
    "remission": [
        "complete remission", "partial remission",
        "minimal residual disease negative", "MRD-negative",
        "deep molecular response",
    ],
    "relapse": [
        "disease relapse", "recurrence", "refractory disease",
        "relapsed/refractory", "second-line therapy",
    ],
    "prognosis": [
        "prognostic factor", "risk stratification",
        "favorable prognosis", "poor prognosis", "survival prediction",
    ],
    "comorbidity": [
        "comorbid conditions", "Charlson comorbidity index",
        "organ dysfunction", "renal impairment",
        "hepatic impairment", "cardiac comorbidity",
    ],
    "molecular profiling": [
        "comprehensive genomic profiling", "CGP",
        "FoundationOne", "Tempus xT", "MSK-IMPACT",
        "tumor molecular testing",
    ],
}

# ---------------------------------------------------------------------------
# 8. Trial expansions
# ---------------------------------------------------------------------------
TRIAL_EXPANSIONS = {
    "clinical trial": [
        "phase III", "randomized controlled",
        "progression-free survival", "overall survival", "enrollment",
    ],
    "phase I": [
        "dose escalation", "first-in-human",
        "maximum tolerated dose", "dose-finding", "safety run-in",
    ],
    "phase II": [
        "signal-finding", "expansion cohort",
        "single-arm trial", "Simon two-stage",
        "objective response rate endpoint",
    ],
    "phase III": [
        "pivotal trial", "registration trial",
        "randomized controlled trial", "primary endpoint PFS",
        "superiority design",
    ],
    "basket trial": [
        "histology-agnostic", "biomarker-selected",
        "tissue-agnostic", "tumor-agnostic", "master protocol",
    ],
    "umbrella trial": [
        "biomarker-driven", "multi-arm trial",
        "adaptive platform", "molecular cohort",
        "precision medicine trial",
    ],
    "real-world evidence": [
        "RWE", "real-world data", "retrospective cohort",
        "electronic health records", "claims data analysis",
    ],
    "endpoint": [
        "primary endpoint", "secondary endpoint",
        "surrogate endpoint", "composite endpoint",
        "time-to-event endpoint",
    ],
    "eligibility": [
        "inclusion criteria", "exclusion criteria",
        "prior lines of therapy", "measurable disease",
        "brain metastases allowed",
    ],
    "randomization": [
        "stratification", "1:1 randomization",
        "crossover allowed", "intent-to-treat",
        "per-protocol population",
    ],
}

# ---------------------------------------------------------------------------
# 9. Immunotherapy expansions
# ---------------------------------------------------------------------------
IMMUNOTHERAPY_EXPANSIONS = {
    "immunotherapy": [
        "checkpoint inhibitor", "anti-PD-1", "anti-PD-L1",
        "anti-CTLA-4", "immune checkpoint blockade",
    ],
    "CAR-T": [
        "chimeric antigen receptor T-cell", "CAR-T cell therapy",
        "adoptive cell therapy", "autologous CAR-T",
        "CD19-directed CAR-T",
    ],
    "bispecific": [
        "bispecific antibody", "bispecific T-cell engager",
        "BiTE", "dual-target immunotherapy",
        "T-cell redirecting antibody",
    ],
    "TIL therapy": [
        "tumor-infiltrating lymphocyte", "lifileucel",
        "adoptive TIL transfer", "ex vivo TIL expansion",
        "TIL-based immunotherapy",
    ],
    "cancer vaccine": [
        "neoantigen vaccine", "mRNA cancer vaccine",
        "dendritic cell vaccine", "personalized vaccine",
        "tumor-associated antigen",
    ],
    "combination immunotherapy": [
        "ipilimumab-nivolumab", "dual checkpoint blockade",
        "PD-1 plus CTLA-4", "chemo-immunotherapy",
        "IO combination",
    ],
    "immune microenvironment": [
        "tumor immune microenvironment", "TIME",
        "hot tumor", "cold tumor", "immune desert",
        "immune-inflamed",
    ],
    "cytokine therapy": [
        "interleukin-2", "IL-2", "interferon-alpha",
        "high-dose IL-2", "cytokine storm",
    ],
    "ADC": [
        "antibody-drug conjugate", "trastuzumab deruxtecan",
        "sacituzumab govitecan", "enfortumab vedotin",
        "payload-linker technology",
    ],
    "oncolytic virus": [
        "oncolytic virotherapy", "T-VEC", "talimogene laherparepvec",
        "viral immunotherapy", "tumor-selective replication",
    ],
    "anti-CTLA-4": [
        "ipilimumab", "tremelimumab", "CTLA-4 blockade",
        "dual checkpoint blockade", "Treg depletion",
    ],
    "TIGIT": [
        "T cell immunoreceptor with Ig and ITIM domains",
        "tiragolumab", "anti-TIGIT antibody",
        "TIGIT checkpoint", "CITYSCAPE trial",
    ],
    "LAG-3": [
        "lymphocyte activation gene 3", "relatlimab",
        "LAG-3 checkpoint inhibitor", "Opdualag",
        "relatlimab-nivolumab combination",
    ],
}

# ---------------------------------------------------------------------------
# 10. Surgery / radiation expansions
# ---------------------------------------------------------------------------
SURGERY_RADIATION_EXPANSIONS = {
    "neoadjuvant": [
        "preoperative", "before surgery", "downstaging",
        "pathologic complete response", "pCR",
    ],
    "adjuvant": [
        "postoperative", "after surgery", "adjuvant chemotherapy",
        "adjuvant immunotherapy", "risk reduction",
    ],
    "radiation": [
        "radiotherapy", "external beam radiation", "EBRT",
        "intensity-modulated radiation therapy", "IMRT",
    ],
    "stereotactic": [
        "SBRT", "stereotactic body radiotherapy",
        "SRS", "stereotactic radiosurgery", "CyberKnife",
    ],
    "proton therapy": [
        "proton beam therapy", "particle therapy",
        "Bragg peak", "pencil beam scanning",
        "proton vs photon",
    ],
    "surgery": [
        "surgical resection", "curative resection",
        "R0 resection", "margin status",
        "minimally invasive surgery",
    ],
    "chemoradiation": [
        "concurrent chemoradiation", "chemoRT",
        "radiosensitization", "definitive chemoradiation",
        "sequential chemoradiation",
    ],
    "brachytherapy": [
        "internal radiation", "interstitial brachytherapy",
        "intracavitary brachytherapy", "high-dose rate",
        "HDR brachytherapy",
    ],
    "ablation": [
        "radiofrequency ablation", "RFA", "microwave ablation",
        "cryoablation", "thermal ablation",
    ],
    "perioperative": [
        "perioperative therapy", "perioperative chemotherapy",
        "perioperative immunotherapy", "sandwich approach",
        "total neoadjuvant therapy", "TNT",
    ],
}

# ---------------------------------------------------------------------------
# 11. Toxicity expansions
# ---------------------------------------------------------------------------
TOXICITY_EXPANSIONS = {
    "side effects": [
        "adverse events", "dose-limiting toxicity",
        "grade 3-4", "treatment discontinuation",
        "immune-related adverse events",
    ],
    "irAE": [
        "immune-related adverse event", "checkpoint toxicity",
        "autoimmune toxicity", "irAE management",
        "steroid taper", "immunosuppression",
    ],
    "hepatotoxicity": [
        "liver toxicity", "transaminase elevation",
        "drug-induced liver injury", "ALT/AST elevation",
        "hepatic failure",
    ],
    "cardiotoxicity": [
        "cardiac toxicity", "QTc prolongation",
        "left ventricular dysfunction", "ejection fraction decrease",
        "myocarditis",
    ],
    "nephrotoxicity": [
        "renal toxicity", "creatinine elevation",
        "acute kidney injury", "proteinuria",
        "renal dose adjustment",
    ],
    "dermatologic toxicity": [
        "rash", "hand-foot syndrome", "acneiform rash",
        "Stevens-Johnson syndrome", "skin toxicity",
    ],
    "neurotoxicity": [
        "peripheral neuropathy", "CNS toxicity",
        "immune effector cell-associated neurotoxicity syndrome",
        "ICANS", "encephalopathy",
    ],
    "myelosuppression": [
        "neutropenia", "thrombocytopenia", "anemia",
        "febrile neutropenia", "bone marrow suppression",
        "G-CSF support",
    ],
    "GI toxicity": [
        "diarrhea", "colitis", "nausea and vomiting",
        "mucositis", "immune-mediated colitis",
    ],
    "pneumonitis": [
        "drug-induced pneumonitis", "interstitial lung disease",
        "immune-mediated pneumonitis", "radiation pneumonitis",
        "checkpoint pneumonitis",
    ],
}

# ---------------------------------------------------------------------------
# 12. Genomics expansions
# ---------------------------------------------------------------------------
GENOMICS_EXPANSIONS = {
    "sequencing": [
        "next-generation sequencing", "NGS",
        "whole exome", "targeted panel", "liquid biopsy", "ctDNA",
    ],
    "WGS": [
        "whole genome sequencing", "germline WGS",
        "tumor-normal WGS", "structural variant detection",
        "copy number analysis",
    ],
    "RNA-seq": [
        "RNA sequencing", "transcriptome profiling",
        "gene expression", "fusion detection",
        "differential expression",
    ],
    "single-cell": [
        "single-cell sequencing", "scRNA-seq",
        "single-cell ATAC-seq", "tumor heterogeneity",
        "clonal architecture",
    ],
    "variant calling": [
        "somatic variant calling", "germline variant calling",
        "variant allele frequency", "VAF",
        "tumor purity", "clonality",
    ],
    "copy number": [
        "copy number alteration", "CNA",
        "gene amplification", "homozygous deletion",
        "chromosomal instability",
    ],
    "epigenomics": [
        "DNA methylation", "histone modification",
        "ATAC-seq", "chromatin accessibility",
        "epigenetic silencing",
    ],
    "proteomics": [
        "mass spectrometry", "protein expression",
        "phosphoproteomics", "RPPA", "IHC quantification",
    ],
    "spatial omics": [
        "spatial transcriptomics", "Visium", "MERFISH",
        "spatial proteomics", "multiplexed imaging",
    ],
    "multi-omics": [
        "multi-omic integration", "genomics-transcriptomics",
        "integrated molecular profiling", "TCGA",
        "pan-cancer atlas",
    ],
    "methylation profiling": [
        "DNA methylation array", "Illumina EPIC array",
        "tumor classification by methylation", "CpG island methylation",
        "methylation-based CNS tumor classifier",
    ],
    "fragmentomics": [
        "cfDNA fragment size analysis", "nucleosome positioning",
        "DELFI fragmentomics", "cell-free DNA fragmentation",
        "fragment-based cancer detection",
    ],
}

# ---------------------------------------------------------------------------
# All expansion maps collected for iteration
# ---------------------------------------------------------------------------
ALL_EXPANSION_MAPS = [
    CANCER_TYPE_EXPANSIONS,
    GENE_EXPANSIONS,
    THERAPY_EXPANSIONS,
    BIOMARKER_EXPANSIONS,
    PATHWAY_EXPANSIONS,
    RESISTANCE_EXPANSIONS,
    CLINICAL_EXPANSIONS,
    TRIAL_EXPANSIONS,
    IMMUNOTHERAPY_EXPANSIONS,
    SURGERY_RADIATION_EXPANSIONS,
    TOXICITY_EXPANSIONS,
    GENOMICS_EXPANSIONS,
]


def expand_query(query: str) -> List[str]:
    """Return up to 10 expansion terms for *query* by checking all maps.

    The function performs a case-insensitive substring match of every key in
    each expansion dictionary against the incoming query.  Matched expansion
    terms are accumulated (preserving insertion order, no duplicates) and the
    first 10 unique terms are returned.

    Parameters
    ----------
    query : str
        Free-text oncology question or keyword string.

    Returns
    -------
    List[str]
        Up to 10 semantically related expansion terms drawn from the 12
        domain-specific maps.
    """
    query_lower = query.lower()
    seen: set = set()
    expansions: List[str] = []

    for exp_map in ALL_EXPANSION_MAPS:
        for key, terms in exp_map.items():
            if key.lower() in query_lower:
                for term in terms:
                    term_lower = term.lower()
                    if term_lower not in seen:
                        seen.add(term_lower)
                        expansions.append(term)
                        if len(expansions) >= 10:
                            return expansions

    return expansions
