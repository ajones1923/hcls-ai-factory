"""Clinical Trial Intelligence Agent — Query Expansion System.

Maps clinical trial terminology to related terms, enabling comprehensive
retrieval across trial databases.  Each expansion map covers a domain
(therapeutic areas, phases, drugs, biomarkers, endpoints, regulatory terms,
trial designs, populations, safety, geography, and procedures) and widens
the semantic net for search queries.

The QueryExpander class orchestrates alias resolution, entity detection,
workflow-aware term boosting, and MeSH-style term expansion.

Author: Adam Jones
Date: March 2026
"""

import re
from typing import Dict, List, Optional

from src.models import TrialWorkflowType


# ═══════════════════════════════════════════════════════════════════════
# 1. ENTITY ALIASES (100+ mappings)
# ═══════════════════════════════════════════════════════════════════════

ENTITY_ALIASES: Dict[str, str] = {
    # --- Cancer types ---
    "NSCLC": "non-small cell lung cancer",
    "SCLC": "small cell lung cancer",
    "mCRC": "metastatic colorectal cancer",
    "CRC": "colorectal cancer",
    "HCC": "hepatocellular carcinoma",
    "RCC": "renal cell carcinoma",
    "TNBC": "triple-negative breast cancer",
    "DLBCL": "diffuse large B-cell lymphoma",
    "CLL": "chronic lymphocytic leukemia",
    "AML": "acute myeloid leukemia",
    "ALL": "acute lymphoblastic leukemia",
    "MDS": "myelodysplastic syndromes",
    "MM": "multiple myeloma",
    "GIST": "gastrointestinal stromal tumor",
    "GBM": "glioblastoma multiforme",
    "NPC": "nasopharyngeal carcinoma",
    "GEJ": "gastroesophageal junction cancer",
    "HNSCC": "head and neck squamous cell carcinoma",
    # --- Cardiovascular ---
    "ACS": "acute coronary syndrome",
    "MI": "myocardial infarction",
    "HFrEF": "heart failure with reduced ejection fraction",
    "HFpEF": "heart failure with preserved ejection fraction",
    "AF": "atrial fibrillation",
    "PAH": "pulmonary arterial hypertension",
    "ASCVD": "atherosclerotic cardiovascular disease",
    "CVD": "cardiovascular disease",
    # --- Neuroscience ---
    "AD": "Alzheimer's disease",
    "PD": "Parkinson's disease",
    "MS": "multiple sclerosis",
    "MDD": "major depressive disorder",
    "ALS": "amyotrophic lateral sclerosis",
    "TRD": "treatment-resistant depression",
    # --- Metabolic ---
    "T2DM": "type 2 diabetes mellitus",
    "T1DM": "type 1 diabetes mellitus",
    "NASH": "nonalcoholic steatohepatitis",
    "MASH": "metabolic dysfunction-associated steatohepatitis",
    "NAFLD": "nonalcoholic fatty liver disease",
    "MASLD": "metabolic dysfunction-associated steatotic liver disease",
    "CKD": "chronic kidney disease",
    "ESKD": "end-stage kidney disease",
    # --- Autoimmune / Inflammatory ---
    "RA": "rheumatoid arthritis",
    "SLE": "systemic lupus erythematosus",
    "IBD": "inflammatory bowel disease",
    "UC": "ulcerative colitis",
    "CD": "Crohn's disease",
    "AS": "ankylosing spondylitis",
    "PsA": "psoriatic arthritis",
    "JIA": "juvenile idiopathic arthritis",
    "SSc": "systemic sclerosis",
    "IPF": "idiopathic pulmonary fibrosis",
    "ITP": "immune thrombocytopenia",
    # --- Respiratory ---
    "COPD": "chronic obstructive pulmonary disease",
    "RSV": "respiratory syncytial virus",
    "CF": "cystic fibrosis",
    # --- Drug brand → generic ---
    "Keytruda": "pembrolizumab",
    "Opdivo": "nivolumab",
    "Tecentriq": "atezolizumab",
    "Imfinzi": "durvalumab",
    "Bavencio": "avelumab",
    "Herceptin": "trastuzumab",
    "Avastin": "bevacizumab",
    "Revlimid": "lenalidomide",
    "Entresto": "sacubitril/valsartan",
    "Jardiance": "empagliflozin",
    "Farxiga": "dapagliflozin",
    "Ozempic": "semaglutide",
    "Wegovy": "semaglutide",
    "Mounjaro": "tirzepatide",
    "Zepbound": "tirzepatide",
    "Humira": "adalimumab",
    "Stelara": "ustekinumab",
    "Dupixent": "dupilumab",
    "Xarelto": "rivaroxaban",
    "Eliquis": "apixaban",
    "Tagrisso": "osimertinib",
    "Enhertu": "trastuzumab deruxtecan",
    "Lynparza": "olaparib",
    "Kisqali": "ribociclib",
    "Ibrance": "palbociclib",
    "Verzenio": "abemaciclib",
    "Kymriah": "tisagenlecleucel",
    "Yescarta": "axicabtagene ciloleucel",
    "Breyanzi": "lisocabtagene maraleucel",
    "Abecma": "idecabtagene vicleucel",
    "Carvykti": "ciltacabtagene autoleucel",
    "Zolgensma": "onasemnogene abeparvovec",
    # --- Biomarker abbreviations ---
    "TMB": "tumor mutational burden",
    "MSI": "microsatellite instability",
    "MSI-H": "microsatellite instability high",
    "dMMR": "deficient mismatch repair",
    "pMMR": "proficient mismatch repair",
    "CPS": "combined positive score (PD-L1)",
    "TPS": "tumor proportion score (PD-L1)",
    "MRD": "minimal residual disease",
    "ctDNA": "circulating tumor DNA",
    "cfDNA": "cell-free DNA",
    "NfL": "neurofilament light chain",
    "FeNO": "fractional exhaled nitric oxide",
    # --- Trial terminology ---
    "RCT": "randomized controlled trial",
    "SAE": "serious adverse event",
    "AE": "adverse event",
    "DSMB": "Data Safety Monitoring Board",
    "DMC": "Data Monitoring Committee",
    "IRB": "Institutional Review Board",
    "IEC": "Independent Ethics Committee",
    "ICF": "informed consent form",
    "CRO": "contract research organization",
    "EDC": "electronic data capture",
    "IRT": "interactive response technology",
    "IWRS": "interactive web response system",
    "eCRF": "electronic case report form",
    "SOC": "standard of care",
    "BSC": "best supportive care",
    "RP2D": "recommended phase 2 dose",
    "MTD": "maximum tolerated dose",
    "DLT": "dose-limiting toxicity",
    "PDAC": "pancreatic ductal adenocarcinoma",
    "MAFLD": "metabolic dysfunction-associated fatty liver disease",
    "PNH": "paroxysmal nocturnal hemoglobinuria",
    "C3G": "C3 glomerulopathy",
    # --- Drug brand → generic (expanded) ---
    "Padcev": "enfortumab vedotin",
    "Trodelvy": "sacituzumab govitecan",
    "Rezdiffra": "resmetirom",
    "Tzield": "teplizumab",
    "Reblozyl": "luspatercept",
    "Truqap": "capivasertib",
    # --- Drug modality abbreviations ---
    "ADC": "antibody-drug conjugate",
    "bispecific": "bispecific antibody",
    "TCE": "T-cell engager",
    # --- Endpoint abbreviations ---
    "ORR": "objective response rate",
    "mPFS": "median progression-free survival",
    "mOS": "median overall survival",
    "DOR": "duration of response",
    "DCR": "disease control rate",
    "CBR": "clinical benefit rate",
    # --- Safety abbreviations ---
    "DILI": "drug-induced liver injury",
    "irAE": "immune-related adverse event",
    "CRS": "cytokine release syndrome",
    # --- Review committee abbreviations ---
    "BICR": "blinded independent central review",
    "IRC": "independent review committee",
}


# ═══════════════════════════════════════════════════════════════════════
# 2. THERAPEUTIC AREA MAP
# ═══════════════════════════════════════════════════════════════════════

THERAPEUTIC_AREA_MAP: Dict[str, List[str]] = {
    "oncology": [
        "cancer", "tumor", "tumour", "neoplasm", "malignancy", "carcinoma",
        "sarcoma", "lymphoma", "leukemia", "melanoma", "immuno-oncology",
        "checkpoint inhibitor", "PD-1", "PD-L1", "CTLA-4", "CAR-T",
        "ADC", "antibody-drug conjugate", "targeted therapy", "chemotherapy",
    ],
    "cardiovascular": [
        "heart", "cardiac", "vascular", "coronary", "arrhythmia",
        "heart failure", "hypertension", "atherosclerosis", "stroke",
        "thrombosis", "anticoagulant", "statin", "SGLT2i", "ARNI",
    ],
    "neuroscience": [
        "neurological", "neurodegenerative", "CNS", "brain", "cognitive",
        "dementia", "Alzheimer", "Parkinson", "epilepsy", "migraine",
        "depression", "schizophrenia", "anxiety", "neuropathy", "ALS",
    ],
    "immunology": [
        "autoimmune", "inflammatory", "immune-mediated", "rheumatology",
        "biologic", "TNF", "IL-17", "IL-23", "JAK inhibitor",
        "anti-inflammatory", "immunosuppressant", "rheumatoid", "lupus",
    ],
    "infectious_disease": [
        "infection", "antibiotic", "antiviral", "antifungal", "vaccine",
        "antimicrobial", "pathogen", "viral", "bacterial", "pandemic",
        "HIV", "hepatitis", "tuberculosis", "malaria", "resistance",
    ],
    "rare_diseases": [
        "orphan drug", "rare disease", "genetic disorder", "gene therapy",
        "enzyme replacement", "substrate reduction", "ultra-rare",
        "inborn error", "lysosomal storage", "neuromuscular", "SMA", "DMD",
    ],
    "metabolic": [
        "diabetes", "obesity", "metabolic syndrome", "HbA1c", "insulin",
        "GLP-1", "GIP", "incretin", "SGLT2", "lipid", "cholesterol",
        "fatty liver", "NASH", "MASH", "weight loss", "BMI",
    ],
    "respiratory": [
        "pulmonary", "lung", "airway", "asthma", "COPD", "fibrosis",
        "bronchitis", "emphysema", "FEV1", "spirometry", "inhaler",
        "bronchodilator", "eosinophilic", "cystic fibrosis",
    ],
    "hematology": [
        "blood", "hematologic", "coagulation", "hemophilia", "anemia",
        "thrombocytopenia", "myeloma", "leukemia", "lymphoma",
        "transplant", "CAR-T", "gene therapy", "sickle cell",
    ],
    "gastroenterology": [
        "gastrointestinal", "GI", "liver", "hepatic", "intestinal",
        "colitis", "Crohn", "IBD", "endoscopy", "GERD", "celiac",
        "cirrhosis", "pancreatitis", "microbiome",
    ],
    "dermatology": [
        "skin", "dermatologic", "cutaneous", "psoriasis", "eczema",
        "atopic dermatitis", "alopecia", "vitiligo", "acne",
        "PASI", "EASI", "IGA", "topical", "biologic",
    ],
    "ophthalmology": [
        "eye", "ocular", "retinal", "ophthalmic", "macular degeneration",
        "glaucoma", "diabetic retinopathy", "dry eye", "anti-VEGF",
        "intravitreal", "visual acuity", "IOP",
    ],
    "gene_cell_therapy": [
        "gene therapy", "cell therapy", "ATMP", "CAR-T", "AAV",
        "lentiviral", "CRISPR", "gene editing", "ex vivo", "in vivo",
        "transgene", "vector", "chimeric antigen receptor",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 3. PHASE MAP
# ═══════════════════════════════════════════════════════════════════════

PHASE_MAP: Dict[str, List[str]] = {
    "phase 1": ["phase I", "first-in-human", "FIH", "dose escalation",
                "dose finding", "MTD", "DLT", "PK study", "safety study"],
    "phase 2": ["phase II", "proof of concept", "POC", "dose ranging",
                "phase 2a", "phase 2b", "signal finding", "exploratory efficacy"],
    "phase 3": ["phase III", "pivotal", "confirmatory", "registration trial",
                "registration study", "pivotal trial", "adequately controlled"],
    "phase 4": ["phase IV", "post-marketing", "post-approval", "PMR", "PMC",
                "real-world evidence", "RWE", "observational", "registry"],
    "phase 1/2": ["phase I/II", "combined phase", "seamless", "first-in-human to POC"],
    "phase 2/3": ["phase II/III", "seamless phase 2/3", "adaptive seamless",
                  "combined pivotal"],
    "preclinical": ["nonclinical", "IND-enabling", "GLP toxicology",
                    "animal study", "in vitro", "in vivo preclinical"],
}


# ═══════════════════════════════════════════════════════════════════════
# 4. DRUG SYNONYM MAP (30+ drugs)
# ═══════════════════════════════════════════════════════════════════════

DRUG_SYNONYM_MAP: Dict[str, List[str]] = {
    "pembrolizumab": ["Keytruda", "MK-3475", "lambrolizumab", "anti-PD-1"],
    "nivolumab": ["Opdivo", "BMS-936558", "MDX-1106", "anti-PD-1"],
    "atezolizumab": ["Tecentriq", "MPDL3280A", "anti-PD-L1"],
    "durvalumab": ["Imfinzi", "MEDI4736", "anti-PD-L1"],
    "ipilimumab": ["Yervoy", "MDX-010", "BMS-734016", "anti-CTLA-4"],
    "trastuzumab": ["Herceptin", "anti-HER2", "rhuMAb HER2"],
    "trastuzumab_deruxtecan": ["Enhertu", "T-DXd", "DS-8201", "fam-trastuzumab deruxtecan"],
    "bevacizumab": ["Avastin", "anti-VEGF", "rhuMAb-VEGF"],
    "osimertinib": ["Tagrisso", "AZD9291", "EGFR TKI", "third-generation EGFR inhibitor"],
    "sotorasib": ["Lumakras", "AMG 510", "KRAS G12C inhibitor"],
    "adagrasib": ["Krazati", "MRTX849", "KRAS G12C inhibitor"],
    "olaparib": ["Lynparza", "AZD2281", "PARP inhibitor"],
    "sacubitril_valsartan": ["Entresto", "LCZ696", "ARNI", "angiotensin receptor neprilysin inhibitor"],
    "empagliflozin": ["Jardiance", "BI 10773", "SGLT2 inhibitor"],
    "dapagliflozin": ["Farxiga", "Forxiga", "BMS-512148", "SGLT2 inhibitor"],
    "semaglutide": ["Ozempic", "Wegovy", "Rybelsus", "GLP-1 receptor agonist", "GLP-1 RA"],
    "tirzepatide": ["Mounjaro", "Zepbound", "LY3298176", "GLP-1/GIP dual agonist"],
    "adalimumab": ["Humira", "D2E7", "anti-TNF", "TNF inhibitor"],
    "ustekinumab": ["Stelara", "CNTO 1275", "anti-IL-12/23"],
    "dupilumab": ["Dupixent", "REGN668", "SAR231893", "anti-IL-4Ralpha"],
    "risankizumab": ["Skyrizi", "BI 655066", "anti-IL-23"],
    "secukinumab": ["Cosentyx", "AIN457", "anti-IL-17A"],
    "abemaciclib": ["Verzenio", "LY2835219", "CDK4/6 inhibitor"],
    "palbociclib": ["Ibrance", "PD-0332991", "CDK4/6 inhibitor"],
    "ribociclib": ["Kisqali", "LEE011", "CDK4/6 inhibitor"],
    "lecanemab": ["Leqembi", "BAN2401", "anti-amyloid", "anti-Abeta protofibril"],
    "donanemab": ["Kisunla", "LY3002813", "anti-amyloid", "anti-N3pG Abeta"],
    "vericiguat": ["Verquvo", "BAY 1021189", "soluble guanylate cyclase stimulator"],
    "evolocumab": ["Repatha", "AMG 145", "PCSK9 inhibitor"],
    "inclisiran": ["Leqvio", "ALN-PCSSC", "PCSK9 siRNA"],
    "onasemnogene_abeparvovec": ["Zolgensma", "AVXS-101", "SMN1 gene therapy"],
    "tisagenlecleucel": ["Kymriah", "CTL019", "anti-CD19 CAR-T"],
    "axicabtagene_ciloleucel": ["Yescarta", "KTE-C19", "anti-CD19 CAR-T"],
}


# ═══════════════════════════════════════════════════════════════════════
# 5. BIOMARKER MAP (20+ biomarkers)
# ═══════════════════════════════════════════════════════════════════════

BIOMARKER_MAP: Dict[str, List[str]] = {
    "PD-L1": ["programmed death-ligand 1", "CD274", "PD-L1 expression",
              "TPS", "CPS", "IC score", "SP263", "22C3", "28-8", "SP142"],
    "TMB": ["tumor mutational burden", "mutational load", "TMB-high",
            "TMB-H", "somatic mutations per megabase"],
    "MSI": ["microsatellite instability", "MSI-H", "MSI-high", "MSI-L",
            "MSS", "microsatellite stable", "Lynch syndrome"],
    "dMMR": ["deficient mismatch repair", "mismatch repair deficient",
             "MLH1", "MSH2", "MSH6", "PMS2", "MMR"],
    "HER2": ["ERBB2", "human epidermal growth factor receptor 2",
             "HER2-positive", "HER2-negative", "HER2-low", "HER2 IHC",
             "HER2 FISH", "HER2 amplification"],
    "EGFR": ["epidermal growth factor receptor", "EGFR mutation",
             "exon 19 deletion", "L858R", "T790M", "C797S", "exon 20 insertion"],
    "BRCA": ["BRCA1", "BRCA2", "BRCA mutation", "homologous recombination deficiency",
             "HRD", "germline BRCA", "somatic BRCA"],
    "ALK": ["anaplastic lymphoma kinase", "ALK rearrangement", "ALK fusion",
            "ALK-positive", "EML4-ALK"],
    "KRAS": ["KRAS G12C", "KRAS G12D", "KRAS G12V", "KRAS mutation",
             "Kirsten rat sarcoma", "RAS mutation"],
    "BRAF": ["BRAF V600E", "BRAF V600K", "BRAF mutation", "BRAF inhibitor"],
    "NTRK": ["NTRK fusion", "neurotrophic tyrosine receptor kinase",
             "NTRK1", "NTRK2", "NTRK3", "TRK fusion"],
    "RET": ["RET fusion", "RET mutation", "RET rearrangement", "RET alteration"],
    "MET": ["MET exon 14 skipping", "MET amplification", "c-MET",
            "hepatocyte growth factor receptor"],
    "FGFR": ["FGFR alteration", "FGFR2 fusion", "FGFR3 mutation",
             "fibroblast growth factor receptor"],
    "PIK3CA": ["PIK3CA mutation", "PI3K pathway", "PI3Kalpha"],
    "NT-proBNP": ["N-terminal pro-brain natriuretic peptide", "BNP",
                  "brain natriuretic peptide", "natriuretic peptide"],
    "troponin": ["troponin I", "troponin T", "high-sensitivity troponin",
                 "hs-cTnI", "hs-cTnT", "cardiac troponin"],
    "HbA1c": ["glycated hemoglobin", "hemoglobin A1c", "A1C",
              "glycosylated hemoglobin"],
    "ctDNA": ["circulating tumor DNA", "cell-free DNA", "cfDNA",
              "liquid biopsy", "molecular residual disease"],
    "MRD": ["minimal residual disease", "measurable residual disease",
            "MRD-negative", "MRD-positive", "flow cytometry MRD", "NGS MRD"],
    "CRP": ["C-reactive protein", "hsCRP", "high-sensitivity CRP",
            "inflammatory marker"],
    "NfL": ["neurofilament light chain", "neurofilament", "NfL serum",
            "NfL CSF", "neurodegeneration biomarker"],
}


# ═══════════════════════════════════════════════════════════════════════
# 6. ENDPOINT MAP
# ═══════════════════════════════════════════════════════════════════════

ENDPOINT_MAP: Dict[str, List[str]] = {
    "OS": ["overall survival", "survival", "mortality", "death",
           "median survival", "survival rate", "OS benefit"],
    "PFS": ["progression-free survival", "time to progression",
            "TTP", "disease progression", "radiographic PFS", "rPFS"],
    "ORR": ["objective response rate", "overall response rate",
            "response rate", "tumor response", "RECIST response"],
    "DFS": ["disease-free survival", "relapse-free survival",
            "recurrence-free survival", "RFS"],
    "CR": ["complete response", "complete remission", "CRR",
           "complete response rate", "pathological complete response", "pCR"],
    "PR": ["partial response", "partial remission", "tumor shrinkage"],
    "DOR": ["duration of response", "response duration",
            "time to response", "TTR", "median DOR"],
    "EFS": ["event-free survival", "composite event endpoint"],
    "MACE": ["major adverse cardiovascular events", "cardiovascular death",
             "myocardial infarction", "stroke", "3-point MACE", "4-point MACE"],
    "QoL": ["quality of life", "health-related quality of life", "HRQoL",
            "EQ-5D", "SF-36", "FACT", "EORTC QLQ"],
    "TTE": ["time to event", "Kaplan-Meier", "survival analysis",
            "hazard ratio", "HR", "log-rank test"],
    "PRO": ["patient-reported outcome", "ePRO", "patient-reported",
            "PGIC", "PGI-C", "NRS", "VAS", "PROMIS"],
    "iDFS": ["invasive disease-free survival", "adjuvant endpoint"],
    "MFS": ["metastasis-free survival", "distant metastasis-free survival"],
    "TTD": ["time to deterioration", "symptom deterioration",
            "functional deterioration"],
}


# ═══════════════════════════════════════════════════════════════════════
# 7. REGULATORY MAP
# ═══════════════════════════════════════════════════════════════════════

REGULATORY_MAP: Dict[str, List[str]] = {
    "IND": ["investigational new drug", "IND application", "IND submission",
            "pre-IND meeting", "IND safety report", "clinical hold"],
    "NDA": ["new drug application", "NDA submission", "NDA filing",
            "NDA approval", "505(b)(1)", "505(b)(2)"],
    "BLA": ["biologics license application", "BLA submission", "BLA approval",
            "biologics", "biological product"],
    "MAA": ["marketing authorisation application", "MAA submission",
            "EMA approval", "centralised procedure", "EU approval"],
    "CRL": ["complete response letter", "refuse to file", "RTF",
            "approvable letter", "not approvable letter"],
    "REMS": ["risk evaluation and mitigation strategy", "risk management",
             "REMS program", "medication guide"],
    "ETASU": ["elements to assure safe use", "restricted distribution",
              "prescriber certification", "patient registry"],
    "BTD": ["breakthrough therapy designation", "breakthrough therapy",
            "breakthrough designation", "expedited development"],
    "FTD": ["fast track designation", "fast track", "rolling review",
            "rolling submission"],
    "AA": ["accelerated approval", "conditional approval", "surrogate endpoint",
           "confirmatory trial", "post-marketing commitment"],
    "PR": ["priority review", "priority review voucher", "PRV",
           "8-month review", "6-month review"],
    "RMAT": ["regenerative medicine advanced therapy", "RMAT designation",
             "cell therapy designation", "gene therapy designation"],
    "PRIME": ["PRIME designation", "priority medicines", "EMA PRIME",
              "enhanced scientific support"],
    "SPA": ["special protocol assessment", "SPA agreement", "study design agreement"],
    "PDUFA": ["Prescription Drug User Fee Act", "PDUFA date", "action date",
              "user fee", "review timeline"],
    "AdCom": ["advisory committee", "FDA advisory committee", "ODAC",
              "advisory panel", "expert panel review"],
    "PAI": ["pre-approval inspection", "GMP inspection", "facility inspection"],
    "RTOR": ["real-time oncology review", "RTOR program", "expedited oncology review"],
    "Orbis": ["Project Orbis", "simultaneous review", "multi-country review",
              "collaborative review"],
}


# ═══════════════════════════════════════════════════════════════════════
# 8. DESIGN MAP
# ═══════════════════════════════════════════════════════════════════════

DESIGN_MAP: Dict[str, List[str]] = {
    "adaptive": ["adaptive design", "interim analysis", "sample size reestimation",
                 "response adaptive randomization", "Bayesian adaptive",
                 "group sequential", "alpha spending"],
    "basket": ["basket trial", "histology-independent", "tumor-agnostic",
               "biomarker-selected", "master protocol"],
    "umbrella": ["umbrella trial", "biomarker-stratified", "multiple biomarker arms",
                 "master protocol", "single disease multiple treatments"],
    "platform": ["platform trial", "perpetual trial", "shared control",
                 "multi-arm multi-stage", "MAMS", "master protocol",
                 "RECOVERY design", "I-SPY design"],
    "crossover": ["crossover design", "cross-over", "AB/BA design",
                  "washout period", "carryover effect", "period effect"],
    "non-inferiority": ["non-inferiority trial", "NI trial", "NI margin",
                        "non-inferiority margin", "delta", "active comparator"],
    "superiority": ["superiority trial", "superiority test", "superiority margin",
                    "placebo-controlled", "active-controlled superiority"],
    "equivalence": ["equivalence trial", "bioequivalence", "therapeutic equivalence",
                    "equivalence margin", "two one-sided tests", "TOST"],
    "factorial": ["factorial design", "2x2 factorial", "multi-factorial",
                  "interaction test", "combination assessment"],
    "enrichment": ["enrichment design", "biomarker enrichment", "predictive enrichment",
                   "prognostic enrichment", "run-in period enrichment"],
    "decentralized": ["DCT", "decentralized clinical trial", "virtual trial",
                      "hybrid trial", "remote trial", "site-less trial"],
    "real_world": ["real-world study", "RWE", "real-world evidence",
                   "pragmatic trial", "EHR-based", "claims-based",
                   "observational study", "registry study"],
    "single_arm": ["single-arm trial", "single-arm study", "uncontrolled",
                   "historical control", "external control", "synthetic control"],
    "seamless": ["seamless design", "phase 2/3", "combined phase",
                 "adaptive seamless", "inferentially seamless"],
}


# ═══════════════════════════════════════════════════════════════════════
# 9. POPULATION MAP
# ═══════════════════════════════════════════════════════════════════════

POPULATION_MAP: Dict[str, List[str]] = {
    "pediatric": ["children", "adolescent", "infant", "neonatal", "juvenile",
                  "pediatric population", "PREA", "pediatric study plan",
                  "pediatric investigation plan", "PIP", "age 0-17"],
    "geriatric": ["elderly", "older adults", "age 65+", "age 75+",
                  "geriatric population", "frailty", "polypharmacy",
                  "ICH E7", "age-related"],
    "treatment_naive": ["treatment-naive", "first-line", "1L", "untreated",
                        "newly diagnosed", "de novo", "frontline"],
    "treatment_experienced": ["previously treated", "prior therapy", "relapsed",
                              "second-line", "2L", "third-line", "3L",
                              "later line", "pretreated"],
    "refractory": ["treatment-refractory", "resistant", "non-responder",
                   "primary resistance", "acquired resistance",
                   "relapsed/refractory", "r/r"],
    "biomarker_positive": ["biomarker-selected", "mutation-positive",
                           "biomarker-enriched", "molecularly selected",
                           "companion diagnostic required"],
    "biomarker_negative": ["biomarker-negative", "wild-type", "WT",
                           "mutation-negative", "all-comers"],
    "pregnant": ["pregnancy", "maternal", "gestational", "lactation",
                 "reproductive toxicology", "DART study"],
    "rare_population": ["ultra-rare", "orphan", "small population",
                        "n-of-1", "limited population"],
    "diverse": ["diversity", "health equity", "underrepresented",
                "minority enrollment", "racial disparity",
                "FDA diversity action plan", "FDORA"],
}


# ═══════════════════════════════════════════════════════════════════════
# 10. SAFETY MAP
# ═══════════════════════════════════════════════════════════════════════

SAFETY_MAP: Dict[str, List[str]] = {
    "SUSAR": ["suspected unexpected serious adverse reaction",
              "unexpected SAR", "expedited safety report", "IND safety report"],
    "SAE": ["serious adverse event", "serious adverse experience",
            "hospitalization", "life-threatening", "death", "disability"],
    "DSMB": ["data safety monitoring board", "data monitoring committee",
             "DMC", "independent monitoring", "interim safety review",
             "DSMB charter", "unblinded review"],
    "CIOMS": ["Council for International Organizations of Medical Sciences",
              "CIOMS form", "CIOMS I", "expedited report", "individual case safety report"],
    "MedDRA": ["Medical Dictionary for Regulatory Activities",
               "adverse event coding", "preferred term", "PT",
               "system organ class", "SOC", "HLGT", "HLT", "LLT"],
    "SOC": ["system organ class", "MedDRA SOC", "organ class",
            "body system classification"],
    "PT": ["preferred term", "MedDRA preferred term", "adverse event term",
           "coded term"],
    "REMS": ["risk evaluation and mitigation strategy", "risk management plan",
             "RMP", "risk minimization", "ETASU"],
    "AESI": ["adverse event of special interest", "targeted adverse event",
             "protocol-defined AE", "AE of clinical interest"],
    "CRS": ["cytokine release syndrome", "cytokine storm", "CAR-T toxicity",
            "immune-related toxicity", "irAE"],
    "DLT": ["dose-limiting toxicity", "dose escalation", "maximum tolerated dose",
            "MTD", "3+3 design", "BOIN", "CRM"],
    "QTc": ["QT prolongation", "QTc interval", "cardiac safety", "TQT study",
            "thorough QT", "ICH E14", "arrhythmia risk"],
    "DILI": ["drug-induced liver injury", "hepatotoxicity", "Hy's Law",
             "ALT elevation", "AST elevation", "liver safety"],
    "irAE": ["immune-related adverse event", "immune-mediated toxicity",
             "checkpoint toxicity", "colitis", "pneumonitis", "hepatitis",
             "thyroiditis", "hypophysitis"],
}


# ═══════════════════════════════════════════════════════════════════════
# 11. GEOGRAPHIC MAP
# ═══════════════════════════════════════════════════════════════════════

GEOGRAPHIC_MAP: Dict[str, List[str]] = {
    "north_america": ["United States", "USA", "US", "Canada",
                      "FDA", "Health Canada", "North American sites"],
    "europe": ["EU", "European Union", "EMA", "MHRA", "UK",
               "Germany", "France", "Spain", "Italy", "Netherlands",
               "EU member states", "European sites"],
    "asia_pacific": ["Japan", "PMDA", "China", "NMPA", "South Korea", "MFDS",
                     "Australia", "TGA", "India", "CDSCO",
                     "Asia-Pacific", "APAC"],
    "latin_america": ["Brazil", "ANVISA", "Mexico", "COFEPRIS",
                      "Argentina", "ANMAT", "Colombia", "Chile",
                      "Latin America", "LATAM"],
    "middle_east_africa": ["South Africa", "SAHPRA", "Israel",
                           "Saudi Arabia", "SFDA", "UAE",
                           "Middle East", "Africa", "MEA"],
    "global": ["global trial", "multinational", "international",
               "multi-regional", "MRCT", "ICH E17",
               "multi-country", "worldwide"],
}


# ═══════════════════════════════════════════════════════════════════════
# 12. PROCEDURE MAP
# ═══════════════════════════════════════════════════════════════════════

PROCEDURE_MAP: Dict[str, List[str]] = {
    "biopsy": ["tissue biopsy", "core biopsy", "fine needle aspiration",
               "FNA", "excisional biopsy", "liquid biopsy", "tumor biopsy"],
    "imaging": ["CT scan", "MRI", "PET scan", "PET-CT", "ultrasound",
                "X-ray", "mammography", "RECIST", "iRECIST",
                "tumor measurement", "radiographic assessment"],
    "endoscopy": ["colonoscopy", "upper endoscopy", "EGD",
                  "bronchoscopy", "cystoscopy", "endoscopic assessment"],
    "surgery": ["resection", "surgical excision", "debulking",
                "transplant", "surgical intervention"],
    "infusion": ["IV infusion", "intravenous administration",
                 "subcutaneous injection", "intramuscular injection",
                 "drug administration", "infusion reaction"],
    "blood_draw": ["phlebotomy", "venipuncture", "blood sample",
                   "serum collection", "plasma collection",
                   "PK sampling", "pharmacokinetic sampling"],
    "genetic_testing": ["next-generation sequencing", "NGS", "PCR",
                        "FISH", "IHC", "immunohistochemistry",
                        "companion diagnostic", "CDx", "molecular profiling",
                        "comprehensive genomic profiling", "CGP"],
    "ecg": ["electrocardiogram", "ECG", "EKG", "12-lead ECG",
            "Holter monitor", "cardiac monitoring", "QTc assessment"],
    "spirometry": ["pulmonary function test", "PFT", "FEV1", "FVC",
                   "DLCO", "lung function test", "bronchodilator reversibility"],
    "vaccination": ["immunization", "vaccine administration",
                    "booster dose", "prime-boost", "seroconversion",
                    "neutralizing antibody", "titer"],
}


# ═══════════════════════════════════════════════════════════════════════
# 13. WORKFLOW TERMS
# ═══════════════════════════════════════════════════════════════════════

_WORKFLOW_TERMS: Dict[TrialWorkflowType, List[str]] = {
    TrialWorkflowType.PATIENT_MATCHING: [
        "clinical trial", "study", "NCT", "ClinicalTrials.gov",
        "enrollment", "recruiting", "active", "completed",
        "inclusion criteria", "exclusion criteria", "eligibility",
    ],
    TrialWorkflowType.PROTOCOL_DESIGN: [
        "protocol", "study design", "randomization", "blinding",
        "stratification", "sample size", "power calculation",
        "inclusion criteria", "exclusion criteria", "endpoints",
        "statistical analysis plan", "SAP",
    ],
    TrialWorkflowType.ENDPOINT_STRATEGY: [
        "endpoint", "primary endpoint", "secondary endpoint",
        "surrogate endpoint", "composite endpoint", "multiplicity",
        "alpha spending", "hierarchical testing", "estimand",
        "ITT", "per-protocol", "sensitivity analysis",
    ],
    TrialWorkflowType.SAFETY_SIGNAL: [
        "safety", "adverse event", "SAE", "SUSAR", "DSMB",
        "safety signal", "disproportionality", "PRR", "ROR",
        "MedDRA", "CIOMS", "pharmacovigilance", "PBRER",
        "benefit-risk", "REMS",
    ],
    TrialWorkflowType.REGULATORY_STRATEGY: [
        "NDA", "BLA", "MAA", "IND", "CTA", "CTD",
        "regulatory submission", "approval", "label",
        "advisory committee", "PDUFA", "complete response letter",
        "FDA", "EMA", "PMDA", "Health Canada",
    ],
    TrialWorkflowType.SITE_SELECTION: [
        "site selection", "investigator", "principal investigator",
        "site feasibility", "enrollment rate", "screen failure",
        "geographic distribution", "site activation",
        "site initiation visit", "SIV",
    ],
    TrialWorkflowType.RECRUITMENT_OPTIMIZATION: [
        "recruitment", "enrollment", "screening", "informed consent",
        "patient identification", "referral", "retention",
        "diversity", "inclusion", "decentralized",
        "digital recruitment", "social media recruitment",
    ],
    TrialWorkflowType.COMPETITIVE_INTELLIGENCE: [
        "competitor", "competitive landscape", "pipeline",
        "market share", "differentiation", "head-to-head",
        "standard of care", "comparator", "market access",
        "pricing", "reimbursement", "HTA",
    ],
    TrialWorkflowType.BIOMARKER_STRATEGY: [
        "biomarker", "companion diagnostic", "CDx",
        "predictive biomarker", "prognostic biomarker",
        "surrogate biomarker", "enrichment", "stratification",
        "molecular profiling", "precision medicine", "targeted therapy",
    ],
    TrialWorkflowType.ADAPTIVE_DESIGN: [
        "adaptive design", "interim analysis", "futility",
        "sample size reestimation", "response adaptive",
        "Bayesian", "group sequential", "platform trial",
        "master protocol", "seamless design",
    ],
    TrialWorkflowType.GENERAL: [
        "clinical trial", "clinical study", "clinical research",
        "investigational", "experimental", "therapeutic",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 14. ENTITY CATEGORIES (for entity detection)
# ═══════════════════════════════════════════════════════════════════════

_ENTITY_CATEGORIES: Dict[str, List[str]] = {
    "therapeutic_areas": [
        "oncology", "cardiovascular", "neuroscience", "immunology",
        "infectious disease", "rare disease", "metabolic", "respiratory",
        "hematology", "gastroenterology", "dermatology", "ophthalmology",
        "gene therapy", "cell therapy",
    ],
    "trial_phases": [
        "phase 1", "phase 2", "phase 3", "phase 4",
        "phase I", "phase II", "phase III", "phase IV",
        "pivotal", "preclinical", "first-in-human",
    ],
    "endpoints": [
        "overall survival", "progression-free survival", "response rate",
        "disease-free survival", "complete response", "partial response",
        "duration of response", "quality of life", "MACE",
        "patient-reported outcome",
    ],
    "regulatory": [
        "FDA", "EMA", "PMDA", "Health Canada", "TGA", "MHRA",
        "IND", "NDA", "BLA", "MAA", "breakthrough therapy",
        "accelerated approval", "fast track", "priority review",
    ],
    "designs": [
        "adaptive", "basket trial", "umbrella trial", "platform trial",
        "crossover", "non-inferiority", "randomized", "double-blind",
        "open-label", "single-arm", "decentralized",
    ],
    "biomarkers": [
        "PD-L1", "TMB", "MSI", "HER2", "EGFR", "BRCA", "ALK",
        "KRAS", "BRAF", "NTRK", "RET", "MET", "ctDNA", "MRD",
    ],
    "drugs": [
        "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab",
        "trastuzumab", "bevacizumab", "osimertinib", "olaparib",
        "empagliflozin", "dapagliflozin", "semaglutide", "tirzepatide",
        "adalimumab", "dupilumab", "lecanemab",
    ],
    "safety_terms": [
        "adverse event", "serious adverse event", "DSMB",
        "dose-limiting toxicity", "safety signal", "REMS",
        "pharmacovigilance", "MedDRA",
    ],
    "populations": [
        "pediatric", "geriatric", "treatment-naive", "refractory",
        "biomarker-positive", "rare population",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 15. TRIAL DESIGN ADVANCED MAP
# ═══════════════════════════════════════════════════════════════════════

TRIAL_DESIGN_ADVANCED_MAP: Dict[str, List[str]] = {
    "bayesian_adaptive": [
        "Bayesian adaptive design", "posterior probability", "predictive probability",
        "Bayesian interim analysis", "Bayesian decision rule", "prior distribution",
    ],
    "response_adaptive_randomization": [
        "RAR", "response adaptive", "outcome adaptive randomization",
        "play-the-winner", "adaptive allocation", "urn model",
    ],
    "biomarker_stratified": [
        "biomarker-stratified design", "stratified medicine", "marker-stratified",
        "predictive biomarker design", "biomarker-guided randomization",
    ],
    "dose_escalation_3plus3": [
        "3+3 design", "rule-based dose escalation", "traditional dose escalation",
        "three plus three", "modified 3+3",
    ],
    "boin_design": [
        "BOIN", "Bayesian optimal interval", "interval-based dose finding",
        "BOIN12", "BOIN-ET",
    ],
    "crm_design": [
        "CRM", "continual reassessment method", "model-based dose escalation",
        "modified CRM", "TITE-CRM", "time-to-event CRM",
    ],
    "keyboard_design": [
        "keyboard design", "key-based dose finding", "mTPI-2",
        "modified toxicity probability interval",
    ],
    "simon_two_stage": [
        "Simon two-stage", "Simon's two-stage design", "optimal two-stage",
        "minimax two-stage", "single-arm phase II",
    ],
    "single_arm_with_historical_control": [
        "single-arm with historical control", "external comparator",
        "historical comparison", "propensity score-matched historical",
    ],
    "external_control_arm": [
        "external control arm", "ECA", "real-world data control",
        "registry-based control", "natural history study control",
    ],
    "synthetic_control": [
        "synthetic control arm", "digital twin control", "in-silico control",
        "synthetic comparator", "model-based control",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 16. PATIENT SUBGROUP MAP
# ═══════════════════════════════════════════════════════════════════════

PATIENT_SUBGROUP_MAP: Dict[str, List[str]] = {
    "treatment_naive": [
        "treatment-naive", "untreated", "no prior therapy", "newly diagnosed",
        "de novo", "chemo-naive", "frontline",
    ],
    "first_line": [
        "first-line", "1L", "initial therapy", "frontline treatment",
        "first-line setting",
    ],
    "second_line": [
        "second-line", "2L", "prior line of therapy", "one prior regimen",
        "second-line setting",
    ],
    "third_line_plus": [
        "third-line", "3L", "3L+", "heavily pretreated", "multiple prior lines",
        "later-line", "fourth-line", "4L", "fifth-line",
    ],
    "maintenance": [
        "maintenance therapy", "maintenance treatment", "continuation therapy",
        "switch maintenance", "continuation maintenance",
    ],
    "adjuvant": [
        "adjuvant therapy", "adjuvant treatment", "post-surgical",
        "adjuvant setting", "adjuvant chemotherapy",
    ],
    "neoadjuvant": [
        "neoadjuvant therapy", "pre-surgical", "induction therapy",
        "neoadjuvant setting", "preoperative treatment",
    ],
    "perioperative": [
        "perioperative therapy", "peri-operative", "pre- and post-surgical",
        "perioperative setting", "surgical window",
    ],
    "metastatic": [
        "metastatic disease", "stage IV", "distant metastasis",
        "metastatic setting", "advanced/metastatic",
    ],
    "locally_advanced": [
        "locally advanced", "unresectable", "stage III",
        "locally advanced unresectable", "non-metastatic",
    ],
    "early_stage": [
        "early-stage", "stage I", "stage II", "operable",
        "resectable", "curable intent",
    ],
    "recurrent": [
        "recurrent disease", "relapsed", "disease recurrence",
        "recurrent/metastatic", "locoregional recurrence",
    ],
    "measurable_disease": [
        "measurable disease", "RECIST measurable", "target lesion",
        "measurable by RECIST 1.1", "bidimensionally measurable",
    ],
    "evaluable_disease": [
        "evaluable disease", "non-measurable disease", "evaluable only",
        "bone-only disease", "non-target lesion",
    ],
    "minimal_residual": [
        "minimal residual disease", "MRD", "MRD-negative", "MRD-positive",
        "molecular residual disease", "measurable residual disease",
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 17. COMPARATIVE PATTERNS
# ═══════════════════════════════════════════════════════════════════════

COMPARATIVE_PATTERNS: List[str] = [
    r"\bvs\.?\b",
    r"\bversus\b",
    r"\bcompared\s+(?:to|with)\b",
    r"\bhead[\s-]to[\s-]head\b",
    r"\bnon[\s-]?inferior(?:ity)?\b",
    r"\bsuperior(?:ity)?\b",
    r"\bequivalen(?:t|ce)\b",
    r"\bbetter\s+than\b",
    r"\bworse\s+than\b",
    r"\bdifference\s+between\b",
]


# ═══════════════════════════════════════════════════════════════════════
# QUERY EXPANDER CLASS
# ═══════════════════════════════════════════════════════════════════════

class QueryExpander:
    """Expand clinical trial queries with related domain terms.

    Orchestrates 12 domain-specific expansion maps, entity alias
    resolution, comparative-query detection, and workflow-aware
    term boosting to maximise recall across trial databases.
    """

    def __init__(self) -> None:
        self.expansion_maps: List[Dict[str, List[str]]] = [
            THERAPEUTIC_AREA_MAP,
            PHASE_MAP,
            DRUG_SYNONYM_MAP,
            BIOMARKER_MAP,
            ENDPOINT_MAP,
            REGULATORY_MAP,
            DESIGN_MAP,
            POPULATION_MAP,
            SAFETY_MAP,
            GEOGRAPHIC_MAP,
            PROCEDURE_MAP,
            TRIAL_DESIGN_ADVANCED_MAP,
            PATIENT_SUBGROUP_MAP,
        ]
        self.entity_aliases: Dict[str, str] = ENTITY_ALIASES
        self._comparative_re = re.compile(
            "|".join(f"({p})" for p in COMPARATIVE_PATTERNS),
            re.IGNORECASE,
        )

    # ───────────────────────────────────────────────────────────────
    # Public API
    # ───────────────────────────────────────────────────────────────

    def expand(
        self,
        query: str,
        workflow: Optional[TrialWorkflowType] = None,
    ) -> dict:
        """Expand a query with related clinical trial terms.

        Parameters
        ----------
        query : str
            The raw user query.
        workflow : TrialWorkflowType, optional
            If provided, additional workflow-specific terms are appended.

        Returns
        -------
        dict
            ``original``            - the raw query text
            ``expanded_terms``      - list of additional search terms
            ``detected_entities``   - dict of categorised entities found
            ``is_comparative``      - whether this is a comparison query
            ``comparison``          - parsed comparison details (or None)
            ``workflow_hint``       - workflow type used for expansion
        """
        query_lower = query.lower().strip()

        # Collect expanded terms from all maps
        expanded_terms: List[str] = []
        for exp_map in self.expansion_maps:
            for trigger, terms in exp_map.items():
                if self._trigger_matches(trigger, query_lower):
                    for term in terms:
                        if term.lower() not in query_lower and term not in expanded_terms:
                            expanded_terms.append(term)

        # Add workflow-specific terms
        if workflow is not None:
            workflow_terms = self.get_workflow_terms(workflow)
            for term in workflow_terms:
                if term.lower() not in query_lower and term not in expanded_terms:
                    expanded_terms.append(term)

        # Detect entities
        detected_entities = self.detect_entities(query)

        # Detect comparative queries
        comparison = self._detect_comparative(query)
        is_comparative = comparison is not None

        return {
            "original": query,
            "expanded_terms": expanded_terms,
            "detected_entities": detected_entities,
            "is_comparative": is_comparative,
            "comparison": comparison,
            "workflow_hint": workflow.value if workflow else None,
        }

    def detect_entities(self, text: str) -> dict:
        """Detect clinical trial entities in text by category.

        Parameters
        ----------
        text : str
            Input text to scan.

        Returns
        -------
        dict
            Keys are category names; values are lists of matched entity strings.
        """
        detected: Dict[str, List[str]] = {}

        for category, entities in _ENTITY_CATEGORIES.items():
            matches: List[str] = []
            for entity in entities:
                pattern = re.compile(
                    r"\b" + re.escape(entity) + r"\b", re.IGNORECASE
                )
                if pattern.search(text):
                    matches.append(entity)
            if matches:
                detected[category] = matches

        # Also detect aliases that appear directly
        alias_matches: List[str] = []
        for abbr in self.entity_aliases:
            pattern = re.compile(r"\b" + re.escape(abbr) + r"\b")
            if pattern.search(text):
                canonical = self.entity_aliases[abbr]
                entry = f"{abbr} ({canonical})"
                if entry not in alias_matches:
                    alias_matches.append(entry)

        if alias_matches:
            detected["resolved_aliases"] = alias_matches

        return detected

    def get_workflow_terms(self, workflow: TrialWorkflowType) -> List[str]:
        """Return additional search terms for a specific workflow.

        Parameters
        ----------
        workflow : TrialWorkflowType
            The clinical trial workflow.

        Returns
        -------
        list of str
            Terms relevant to the workflow.
        """
        return list(_WORKFLOW_TERMS.get(workflow, []))

    def expand_mesh_terms(self, term: str) -> dict:
        """Expand a term with broader, narrower, and related MeSH-style terms.

        Parameters
        ----------
        term : str
            The term to expand.

        Returns
        -------
        dict
            ``broader``  - broader/parent terms
            ``narrower`` - narrower/child terms
            ``related``  - related terms at same level
        """
        term_lower = term.lower().strip()

        broader: List[str] = []
        narrower: List[str] = []
        related: List[str] = []

        # Check therapeutic area maps for hierarchy
        for area, terms in THERAPEUTIC_AREA_MAP.items():
            terms_lower = [t.lower() for t in terms]
            if term_lower == area:
                narrower.extend(terms)
            elif term_lower in terms_lower:
                broader.append(area)
                related.extend(
                    t for t in terms
                    if t.lower() != term_lower
                )

        # Check biomarker maps
        for biomarker, terms in BIOMARKER_MAP.items():
            terms_lower = [t.lower() for t in terms]
            if term_lower == biomarker.lower():
                narrower.extend(terms)
                broader.append("biomarker")
            elif term_lower in terms_lower:
                broader.append(biomarker)
                related.extend(
                    t for t in terms
                    if t.lower() != term_lower
                )

        # Check drug synonym maps
        for drug, synonyms in DRUG_SYNONYM_MAP.items():
            synonyms_lower = [s.lower() for s in synonyms]
            if term_lower == drug.lower():
                related.extend(synonyms)
            elif term_lower in synonyms_lower:
                broader.append(drug)
                related.extend(
                    s for s in synonyms
                    if s.lower() != term_lower
                )

        # Check endpoint maps
        for endpoint, terms in ENDPOINT_MAP.items():
            terms_lower = [t.lower() for t in terms]
            if term_lower == endpoint.lower():
                narrower.extend(terms)
                broader.append("endpoint")
            elif term_lower in terms_lower:
                broader.append(endpoint)
                related.extend(
                    t for t in terms
                    if t.lower() != term_lower
                )

        # Deduplicate
        broader = list(dict.fromkeys(broader))
        narrower = list(dict.fromkeys(narrower))
        related = list(dict.fromkeys(related))

        return {
            "broader": broader,
            "narrower": narrower,
            "related": related,
        }

    # ───────────────────────────────────────────────────────────────
    # Internal helpers
    # ───────────────────────────────────────────────────────────────

    @staticmethod
    def _trigger_matches(trigger: str, query_lower: str) -> bool:
        """Check if a trigger term appears as a whole word in the query."""
        pattern = re.compile(r"\b" + re.escape(trigger.lower()) + r"\b")
        return bool(pattern.search(query_lower))

    def _detect_comparative(self, text: str) -> Optional[dict]:
        """Detect and parse comparative queries."""
        match = self._comparative_re.search(text)
        if match is None:
            return None

        result: dict = {
            "pattern_matched": match.group(0),
            "raw_sides": None,
        }

        vs_pattern = re.compile(
            r"(.+?)\s+(?:vs\.?|versus|compared\s+to|compared\s+with)\s+(.+)",
            re.IGNORECASE,
        )
        vs_match = vs_pattern.search(text)
        if vs_match:
            result["raw_sides"] = (
                vs_match.group(1).strip(),
                vs_match.group(2).strip(),
            )

        return result
