"""Domain-specific query expansion for the Pharmacogenomics Intelligence Agent.

Provides 14 expansion maps that enrich user queries with pharmacogenomic
synonyms, related genes, drug names, phenotypes, and clinical concepts.
This dramatically improves recall when searching across the 15 Milvus
collections, because users may phrase queries using brand names, gene
aliases, phenotype shorthand, or disease-oriented language rather than
the precise terminology stored in the vector database.

Follows the same expansion-map pattern used by:
  - cart_intelligence_agent/src/query_expansion.py
  - biomarker_intelligence_agent/src/query_expansion.py

Usage:
    from src.query_expansion import expand_query, expand_query_by_category

    terms = expand_query("warfarin dosing for CYP2C9 poor metabolizer")
    categorized = expand_query_by_category("carbamazepine SJS Asian")
"""

from typing import Dict, List, Tuple

# ===================================================================
# 1. DRUG_EXPANSION  (40+ entries)
# Maps drug names (generic / brand) to related PGx terms.
# ===================================================================

DRUG_EXPANSION: Dict[str, List[str]] = {
    # --- Anticoagulants ---
    "warfarin": [
        "coumadin", "anticoagulant", "INR", "blood thinner", "CYP2C9",
        "VKORC1", "vitamin K antagonist", "bleeding risk", "dose adjustment",
        "CYP4F2", "thromboembolism", "atrial fibrillation", "DVT",
        "international normalized ratio", "prothrombin time",
    ],
    # --- Opioids ---
    "codeine": [
        "opioid", "pain", "CYP2D6", "morphine", "prodrug", "ultra-rapid",
        "poor metabolizer", "respiratory depression", "analgesic",
        "tylenol 3", "acetaminophen with codeine", "tramadol alternative",
    ],
    "tramadol": [
        "opioid", "pain", "CYP2D6", "O-desmethyltramadol", "analgesic",
        "prodrug", "seizure risk", "serotonin syndrome", "ultra-rapid metabolizer",
    ],
    "oxycodone": [
        "opioid", "pain", "CYP2D6", "CYP3A4", "oxymorphone", "analgesic",
        "controlled substance", "oxycontin", "percocet",
    ],
    "hydrocodone": [
        "opioid", "pain", "CYP2D6", "hydromorphone", "vicodin",
        "norco", "analgesic", "prodrug",
    ],
    # --- Antiplatelets ---
    "clopidogrel": [
        "plavix", "antiplatelet", "stent", "CYP2C19", "prasugrel",
        "ticagrelor", "stent thrombosis", "ACS", "PCI", "acute coronary syndrome",
        "percutaneous coronary intervention", "dual antiplatelet therapy", "DAPT",
        "loss of function", "poor metabolizer",
    ],
    # --- Statins ---
    "simvastatin": [
        "zocor", "statin", "SLCO1B1", "myopathy", "rhabdomyolysis",
        "cholesterol", "pravastatin", "rosuvastatin", "HMG-CoA reductase",
        "lipid-lowering", "muscle toxicity", "CPK elevation",
    ],
    "atorvastatin": [
        "lipitor", "statin", "SLCO1B1", "CYP3A4", "cholesterol",
        "HMG-CoA reductase", "myopathy", "lipid-lowering",
    ],
    "rosuvastatin": [
        "crestor", "statin", "SLCO1B1", "ABCG2", "cholesterol",
        "HMG-CoA reductase", "myopathy", "lipid-lowering",
    ],
    "pravastatin": [
        "pravachol", "statin", "SLCO1B1", "cholesterol", "HMG-CoA reductase",
        "hydrophilic statin", "lipid-lowering",
    ],
    # --- Breast cancer ---
    "tamoxifen": [
        "breast cancer", "CYP2D6", "endoxifen", "ER-positive", "SERM",
        "selective estrogen receptor modulator", "nolvadex", "adjuvant therapy",
        "poor metabolizer", "aromatase inhibitor alternative",
    ],
    # --- Chemotherapy ---
    "5-fluorouracil": [
        "5-FU", "capecitabine", "DPYD", "DPD deficiency", "mucositis",
        "neutropenia", "chemotherapy toxicity", "fluoropyrimidine",
        "dihydropyrimidine dehydrogenase", "dose reduction", "xeloda",
        "colorectal cancer", "gastrointestinal cancer",
    ],
    "capecitabine": [
        "xeloda", "5-FU prodrug", "DPYD", "DPD deficiency", "fluoropyrimidine",
        "oral chemotherapy", "breast cancer", "colorectal cancer",
        "mucositis", "neutropenia", "hand-foot syndrome",
    ],
    "irinotecan": [
        "camptosar", "UGT1A1", "SN-38", "neutropenia", "diarrhea",
        "Gilbert syndrome", "UGT1A1*28", "topoisomerase inhibitor",
        "colorectal cancer", "dose reduction",
    ],
    "mercaptopurine": [
        "6-MP", "purinethol", "TPMT", "NUDT15", "thiopurine",
        "leukemia", "ALL", "myelosuppression", "azathioprine related",
        "immunosuppressant", "dose reduction",
    ],
    "thioguanine": [
        "6-TG", "tabloid", "TPMT", "NUDT15", "thiopurine",
        "leukemia", "myelosuppression",
    ],
    "cisplatin": [
        "platinol", "TPMT", "ototoxicity", "nephrotoxicity",
        "platinum-based", "chemotherapy", "hearing loss",
    ],
    # --- Antivirals ---
    "abacavir": [
        "HIV", "HLA-B*57:01", "hypersensitivity", "antiretroviral",
        "ziagen", "nucleoside reverse transcriptase inhibitor", "NRTI",
        "mandatory screening", "abacavir hypersensitivity reaction",
    ],
    # --- Anticonvulsants ---
    "carbamazepine": [
        "tegretol", "epilepsy", "HLA-B*15:02", "HLA-A*31:01", "SJS", "TEN",
        "DRESS", "Stevens-Johnson syndrome", "toxic epidermal necrolysis",
        "anticonvulsant", "mood stabilizer", "trigeminal neuralgia",
        "Southeast Asian", "screening",
    ],
    "phenytoin": [
        "dilantin", "epilepsy", "CYP2C9", "HLA-B*15:02", "SJS", "TEN",
        "anticonvulsant", "zero-order kinetics", "narrow therapeutic index",
        "dose adjustment", "Stevens-Johnson syndrome",
    ],
    "oxcarbazepine": [
        "trileptal", "epilepsy", "HLA-B*15:02", "SJS", "TEN",
        "anticonvulsant", "carbamazepine alternative",
    ],
    "lamotrigine": [
        "lamictal", "epilepsy", "bipolar", "SJS", "TEN",
        "anticonvulsant", "mood stabilizer", "HLA-B*15:02",
    ],
    # --- Immunosuppressants ---
    "tacrolimus": [
        "immunosuppressant", "transplant", "CYP3A5", "organ rejection",
        "calcineurin inhibitor", "prograf", "kidney transplant",
        "liver transplant", "trough level", "therapeutic drug monitoring",
        "CYP3A4", "narrow therapeutic index",
    ],
    "azathioprine": [
        "imuran", "TPMT", "NUDT15", "thiopurine", "immunosuppressant",
        "autoimmune", "IBD", "myelosuppression", "dose reduction",
        "rheumatoid arthritis", "Crohn disease",
    ],
    "mycophenolate": [
        "cellcept", "myfortic", "UGT1A9", "immunosuppressant", "transplant",
        "organ rejection", "mycophenolic acid", "MPA",
    ],
    # --- Antidepressants (SSRIs) ---
    "fluoxetine": [
        "prozac", "SSRI", "CYP2D6", "antidepressant", "depression",
        "serotonin", "CYP2D6 inhibitor", "strong inhibitor",
        "phenoconversion", "drug interaction",
    ],
    "paroxetine": [
        "paxil", "SSRI", "CYP2D6", "antidepressant", "depression",
        "anxiety", "strong CYP2D6 inhibitor", "phenoconversion",
        "serotonin", "dose adjustment",
    ],
    "sertraline": [
        "zoloft", "SSRI", "CYP2C19", "CYP2B6", "antidepressant",
        "depression", "anxiety", "serotonin", "dose adjustment",
    ],
    "escitalopram": [
        "lexapro", "SSRI", "CYP2C19", "antidepressant", "depression",
        "anxiety", "serotonin", "citalopram", "dose adjustment",
        "poor metabolizer", "ultra-rapid metabolizer",
    ],
    "citalopram": [
        "celexa", "SSRI", "CYP2C19", "antidepressant", "depression",
        "QT prolongation", "serotonin", "dose adjustment",
    ],
    "fluvoxamine": [
        "luvox", "SSRI", "CYP2D6", "CYP1A2", "antidepressant",
        "OCD", "strong CYP1A2 inhibitor", "drug interaction",
    ],
    # --- Antidepressants (TCAs) ---
    "amitriptyline": [
        "elavil", "TCA", "tricyclic antidepressant", "CYP2D6", "CYP2C19",
        "nortriptyline", "depression", "neuropathic pain", "migraine",
        "dose adjustment", "poor metabolizer", "ultra-rapid metabolizer",
    ],
    "nortriptyline": [
        "pamelor", "TCA", "tricyclic antidepressant", "CYP2D6",
        "depression", "neuropathic pain", "dose adjustment",
        "amitriptyline metabolite", "therapeutic drug monitoring",
    ],
    "imipramine": [
        "tofranil", "TCA", "tricyclic antidepressant", "CYP2D6", "CYP2C19",
        "depression", "enuresis", "desipramine metabolite",
    ],
    "desipramine": [
        "norpramin", "TCA", "tricyclic antidepressant", "CYP2D6",
        "depression", "dose adjustment", "imipramine metabolite",
    ],
    "clomipramine": [
        "anafranil", "TCA", "tricyclic antidepressant", "CYP2D6", "CYP2C19",
        "OCD", "depression", "desmethylclomipramine",
    ],
    "doxepin": [
        "sinequan", "TCA", "tricyclic antidepressant", "CYP2D6", "CYP2C19",
        "depression", "insomnia", "antihistamine",
    ],
    # --- Antidepressants (other) ---
    "venlafaxine": [
        "effexor", "SNRI", "CYP2D6", "antidepressant", "depression",
        "O-desmethylvenlafaxine", "desvenlafaxine", "anxiety",
        "dose adjustment", "serotonin norepinephrine",
    ],
    "vortioxetine": [
        "trintellix", "CYP2D6", "antidepressant", "multimodal",
        "serotonin modulator", "depression", "dose adjustment",
    ],
    # --- Antipsychotics ---
    "aripiprazole": [
        "abilify", "antipsychotic", "CYP2D6", "CYP3A4", "schizophrenia",
        "bipolar", "dose adjustment", "partial dopamine agonist",
        "poor metabolizer",
    ],
    "clozapine": [
        "clozaril", "antipsychotic", "CYP1A2", "CYP2D6",
        "schizophrenia", "agranulocytosis", "treatment-resistant",
        "REMS", "ANC monitoring", "metabolic syndrome",
    ],
    "haloperidol": [
        "haldol", "antipsychotic", "CYP2D6", "CYP3A4",
        "schizophrenia", "QT prolongation", "EPS",
        "extrapyramidal symptoms", "dose adjustment",
    ],
    "risperidone": [
        "risperdal", "antipsychotic", "CYP2D6", "paliperidone",
        "schizophrenia", "bipolar", "dose adjustment",
    ],
    "pimozide": [
        "orap", "antipsychotic", "CYP2D6", "CYP3A4",
        "Tourette syndrome", "QT prolongation", "contraindicated",
    ],
    # --- Proton pump inhibitors ---
    "omeprazole": [
        "prilosec", "proton pump inhibitor", "PPI", "CYP2C19",
        "GERD", "acid reflux", "Helicobacter pylori", "dose adjustment",
        "poor metabolizer", "ultra-rapid metabolizer",
    ],
    "lansoprazole": [
        "prevacid", "proton pump inhibitor", "PPI", "CYP2C19",
        "GERD", "acid reflux", "dose adjustment",
    ],
    "pantoprazole": [
        "protonix", "proton pump inhibitor", "PPI", "CYP2C19",
        "GERD", "acid reflux",
    ],
    "dexlansoprazole": [
        "dexilant", "proton pump inhibitor", "PPI", "CYP2C19",
        "GERD", "acid reflux",
    ],
    # --- Anti-gout ---
    "allopurinol": [
        "zyloprim", "HLA-B*58:01", "SJS", "TEN", "gout",
        "hyperuricemia", "xanthine oxidase inhibitor", "uric acid",
        "hypersensitivity", "DRESS", "febuxostat alternative",
    ],
    # --- Cardiovascular ---
    "metoprolol": [
        "lopressor", "toprol", "beta blocker", "CYP2D6",
        "heart failure", "hypertension", "atrial fibrillation",
        "dose adjustment", "poor metabolizer",
    ],
    "propranolol": [
        "inderal", "beta blocker", "CYP2D6", "hypertension",
        "migraine", "tremor", "anxiety",
    ],
    "propafenone": [
        "rythmol", "antiarrhythmic", "CYP2D6", "atrial fibrillation",
        "dose adjustment", "poor metabolizer",
    ],
    "ivacaftor": [
        "kalydeco", "CFTR", "cystic fibrosis", "G551D", "CFTR modulator",
        "potentiator", "sweat chloride",
    ],
    # --- Antifungals ---
    "voriconazole": [
        "vfend", "antifungal", "CYP2C19", "azole", "aspergillosis",
        "invasive fungal infection", "therapeutic drug monitoring",
        "dose adjustment", "poor metabolizer", "ultra-rapid metabolizer",
    ],
    # --- Antiemetics ---
    "ondansetron": [
        "zofran", "antiemetic", "CYP2D6", "serotonin 5-HT3 antagonist",
        "nausea", "vomiting", "chemotherapy-induced", "ultra-rapid metabolizer",
    ],
    # --- Autoimmune / Anti-inflammatory ---
    "sulfasalazine": [
        "azulfidine", "NAT2", "rheumatoid arthritis", "ulcerative colitis",
        "IBD", "slow acetylator", "rapid acetylator",
    ],
    # --- Atomoxetine ---
    "atomoxetine": [
        "strattera", "CYP2D6", "ADHD", "attention deficit",
        "norepinephrine reuptake inhibitor", "dose adjustment",
        "poor metabolizer",
    ],
    # --- Eliglustat ---
    "eliglustat": [
        "cerdelga", "CYP2D6", "Gaucher disease", "glucosylceramide synthase",
        "dose adjustment", "poor metabolizer", "contraindicated",
    ],
    # --- Rasburicase ---
    "rasburicase": [
        "elitek", "G6PD", "glucose-6-phosphate dehydrogenase",
        "tumor lysis syndrome", "hemolytic anemia", "contraindicated",
        "uric acid", "methemoglobinemia",
    ],
    # --- Primaquine ---
    "primaquine": [
        "G6PD", "glucose-6-phosphate dehydrogenase", "malaria",
        "hemolytic anemia", "Plasmodium vivax", "screening required",
    ],
    # --- Dapsone ---
    "dapsone": [
        "G6PD", "leprosy", "dermatitis herpetiformis", "PCP prophylaxis",
        "hemolytic anemia", "methemoglobinemia",
    ],
    # --- Trastuzumab ---
    "trastuzumab": [
        "herceptin", "HER2", "ERBB2", "breast cancer", "gastric cancer",
        "monoclonal antibody", "HER2-positive", "companion diagnostic",
    ],
    # --- Cetuximab ---
    "cetuximab": [
        "erbitux", "EGFR", "KRAS", "RAS", "colorectal cancer",
        "KRAS wild-type", "monoclonal antibody", "companion diagnostic",
    ],
}

# ===================================================================
# 2. GENE_EXPANSION  (25 entries)
# Maps pharmacogene symbols to related terms.
# ===================================================================

GENE_EXPANSION: Dict[str, List[str]] = {
    "CYP2D6": [
        "cytochrome P450 2D6", "debrisoquine hydroxylase", "22q13.2",
        "opioids", "antidepressants", "tamoxifen", "poor metabolizer",
        "ultra-rapid metabolizer", "intermediate metabolizer",
        "gene duplication", "star alleles", "activity score",
    ],
    "CYP2C19": [
        "cytochrome P450 2C19", "S-mephenytoin hydroxylase", "10q23.33",
        "clopidogrel", "proton pump inhibitor", "voriconazole", "SSRI",
        "poor metabolizer", "ultra-rapid metabolizer", "loss of function",
        "*2", "*3", "*17",
    ],
    "CYP2C9": [
        "cytochrome P450 2C9", "10q23.33", "warfarin", "phenytoin",
        "NSAIDs", "poor metabolizer", "intermediate metabolizer",
        "*2", "*3", "vitamin K antagonist",
    ],
    "CYP3A5": [
        "cytochrome P450 3A5", "7q22.1", "tacrolimus", "transplant",
        "expresser", "non-expresser", "*1", "*3",
        "immunosuppressant metabolism",
    ],
    "CYP3A4": [
        "cytochrome P450 3A4", "7q22.1", "drug metabolism",
        "most abundant CYP", "CYP3A5", "midazolam", "cyclosporine",
        "enzyme induction", "enzyme inhibition",
    ],
    "CYP1A2": [
        "cytochrome P450 1A2", "15q24.1", "clozapine", "caffeine",
        "theophylline", "fluvoxamine", "smoking induction",
        "polycyclic aromatic hydrocarbons",
    ],
    "CYP2B6": [
        "cytochrome P450 2B6", "19q13.2", "efavirenz", "methadone",
        "cyclophosphamide", "bupropion", "ketamine",
    ],
    "CYP4F2": [
        "cytochrome P450 4F2", "19p13.12", "warfarin",
        "vitamin K oxidase", "VKORC1 interaction",
    ],
    "VKORC1": [
        "vitamin K epoxide reductase complex subunit 1", "16p11.2",
        "warfarin", "anticoagulant", "vitamin K cycle",
        "-1639G>A", "rs9923231", "dose sensitivity",
    ],
    "SLCO1B1": [
        "solute carrier organic anion transporter 1B1", "OATP1B1",
        "12p12.1", "simvastatin", "statin", "myopathy",
        "*5", "*15", "rs4149056", "hepatic uptake transporter",
    ],
    "ABCB1": [
        "ATP binding cassette subfamily B member 1", "P-glycoprotein",
        "MDR1", "7q21.12", "drug efflux", "blood-brain barrier",
        "multidrug resistance", "digoxin", "C3435T",
    ],
    "ABCG2": [
        "ATP binding cassette subfamily G member 2", "BCRP",
        "breast cancer resistance protein", "4q22.1",
        "rosuvastatin", "drug efflux",
    ],
    "DPYD": [
        "dihydropyrimidine dehydrogenase", "1p21.3",
        "5-fluorouracil", "capecitabine", "DPD deficiency",
        "fluoropyrimidine toxicity", "*2A", "c.1905+1G>A",
        "pre-treatment screening",
    ],
    "TPMT": [
        "thiopurine S-methyltransferase", "6p22.3",
        "azathioprine", "mercaptopurine", "thioguanine",
        "thiopurine", "myelosuppression", "*2", "*3A", "*3C",
    ],
    "NUDT15": [
        "nudix hydrolase 15", "13q14.2", "thiopurine",
        "azathioprine", "mercaptopurine",
        "East Asian", "myelosuppression", "R139C",
    ],
    "UGT1A1": [
        "UDP glucuronosyltransferase 1A1", "2q37.1",
        "irinotecan", "atazanavir", "Gilbert syndrome",
        "*28", "*6", "bilirubin", "glucuronidation",
    ],
    "G6PD": [
        "glucose-6-phosphate dehydrogenase", "Xq28",
        "rasburicase", "primaquine", "dapsone",
        "hemolytic anemia", "X-linked", "favism",
        "oxidative stress", "methemoglobinemia",
    ],
    "NAT2": [
        "N-acetyltransferase 2", "8p22", "isoniazid",
        "sulfasalazine", "hydralazine", "slow acetylator",
        "rapid acetylator", "drug-induced lupus",
    ],
    "IFNL3": [
        "interferon lambda 3", "IL28B", "19q13.2",
        "hepatitis C", "peginterferon", "ribavirin",
        "sustained virologic response", "rs12979860",
    ],
    "CYP2C8": [
        "cytochrome P450 2C8", "10q23.33", "paclitaxel",
        "amodiaquine", "repaglinide", "thiazolidinediones",
    ],
    "RYR1": [
        "ryanodine receptor 1", "19q13.2",
        "malignant hyperthermia", "volatile anesthetics",
        "succinylcholine", "dantrolene",
    ],
    "CACNA1S": [
        "calcium voltage-gated channel subunit alpha1 S", "1q32.1",
        "malignant hyperthermia", "volatile anesthetics",
    ],
    "HLA-A": [
        "human leukocyte antigen A", "6p21.3", "MHC class I",
        "carbamazepine", "HLA-A*31:01", "DRESS",
        "drug hypersensitivity", "immune-mediated",
    ],
    "HLA-B": [
        "human leukocyte antigen B", "6p21.3", "MHC class I",
        "abacavir", "carbamazepine", "allopurinol",
        "HLA-B*57:01", "HLA-B*15:02", "HLA-B*58:01",
        "SJS", "TEN", "hypersensitivity",
    ],
    "CFTR": [
        "cystic fibrosis transmembrane conductance regulator", "7q31.2",
        "ivacaftor", "lumacaftor", "tezacaftor", "elexacaftor",
        "cystic fibrosis", "F508del", "G551D", "CFTR modulator",
    ],
}

# ===================================================================
# 3. PHENOTYPE_EXPANSION  (5 entries)
# Maps metabolizer phenotype terms to related concepts.
# ===================================================================

PHENOTYPE_EXPANSION: Dict[str, List[str]] = {
    "poor metabolizer": [
        "PM", "no function", "reduced metabolism", "drug accumulation",
        "toxicity risk", "dose reduction", "two no-function alleles",
        "activity score 0", "null alleles", "adverse drug reaction",
    ],
    "intermediate metabolizer": [
        "IM", "decreased function", "reduced activity",
        "one functional allele", "moderate dose reduction",
        "activity score 0.5-1.0", "partial metabolism",
    ],
    "normal metabolizer": [
        "NM", "extensive metabolizer", "EM", "normal function",
        "two functional alleles", "standard dosing",
        "activity score 1.0-2.0", "wild-type",
    ],
    "ultra-rapid metabolizer": [
        "UM", "increased activity", "rapid activation",
        "subtherapeutic levels", "prodrug toxicity",
        "gene duplication", "increased function",
        "activity score >2.0", "treatment failure",
    ],
    "rapid metabolizer": [
        "RM", "increased function", "one increased-function allele",
        "activity score 1.5-2.0", "CYP2C19*17",
    ],
}

# ===================================================================
# 4. TOXICITY_EXPANSION  (15 entries)
# Maps adverse drug reaction types to PGx-related terms.
# ===================================================================

TOXICITY_EXPANSION: Dict[str, List[str]] = {
    "myopathy": [
        "muscle pain", "rhabdomyolysis", "CPK elevation", "CK elevation",
        "statin", "SLCO1B1", "muscle weakness", "myalgia",
        "creatine kinase", "dose-dependent",
    ],
    "SJS/TEN": [
        "Stevens-Johnson syndrome", "toxic epidermal necrolysis",
        "carbamazepine", "HLA-B*15:02", "phenytoin", "allopurinol",
        "HLA-B*58:01", "skin blistering", "mucosal erosion",
        "life-threatening", "dermatologic emergency",
    ],
    "serotonin syndrome": [
        "serotonergic", "SSRI", "CYP2D6", "drug interaction",
        "hyperthermia", "myoclonus", "agitation", "diaphoresis",
        "serotonin toxicity", "linezolid", "MAO inhibitor",
    ],
    "DRESS": [
        "drug reaction with eosinophilia and systemic symptoms",
        "HLA-A*31:01", "carbamazepine", "allopurinol",
        "eosinophilia", "organ involvement", "hepatitis",
        "rash", "fever", "lymphadenopathy",
    ],
    "neutropenia": [
        "low neutrophils", "ANC", "absolute neutrophil count",
        "bone marrow suppression", "myelosuppression",
        "TPMT", "NUDT15", "DPYD", "clozapine",
        "febrile neutropenia", "infection risk",
    ],
    "hepatotoxicity": [
        "liver injury", "DILI", "drug-induced liver injury",
        "ALT elevation", "AST elevation", "jaundice",
        "hepatic failure", "NAT2", "isoniazid",
    ],
    "nephrotoxicity": [
        "kidney injury", "renal toxicity", "acute kidney injury",
        "creatinine elevation", "GFR decline",
        "cisplatin", "aminoglycosides",
    ],
    "ototoxicity": [
        "hearing loss", "tinnitus", "cisplatin", "aminoglycosides",
        "TPMT", "MT-RNR1", "cochlear damage",
    ],
    "QT prolongation": [
        "long QT", "torsades de pointes", "cardiac arrhythmia",
        "CYP2D6", "citalopram", "haloperidol", "pimozide",
        "sudden cardiac death", "ECG monitoring",
    ],
    "bleeding": [
        "hemorrhage", "warfarin", "CYP2C9", "VKORC1",
        "INR supratherapeutic", "anticoagulant toxicity",
        "gastrointestinal bleed", "intracranial hemorrhage",
    ],
    "hemolytic anemia": [
        "G6PD", "glucose-6-phosphate dehydrogenase deficiency",
        "rasburicase", "primaquine", "dapsone",
        "oxidative stress", "hemolysis", "methemoglobinemia",
    ],
    "hypersensitivity": [
        "drug allergy", "anaphylaxis", "HLA", "immune-mediated",
        "abacavir", "HLA-B*57:01", "skin reaction",
        "maculopapular rash", "fever",
    ],
    "myelosuppression": [
        "bone marrow suppression", "pancytopenia", "leukopenia",
        "thrombocytopenia", "TPMT", "NUDT15", "thiopurine",
        "dose-limiting toxicity",
    ],
    "mucositis": [
        "oral mucositis", "stomatitis", "DPYD", "5-fluorouracil",
        "chemotherapy", "mouth sores", "gastrointestinal",
    ],
    "malignant hyperthermia": [
        "MH", "RYR1", "CACNA1S", "volatile anesthetics",
        "succinylcholine", "dantrolene", "muscle rigidity",
        "hyperthermia", "hyperkalemia", "anesthesia emergency",
    ],
}

# ===================================================================
# 5. GUIDELINE_EXPANSION  (5 entries)
# Maps pharmacogenomics guideline organizations to related terms.
# ===================================================================

GUIDELINE_EXPANSION: Dict[str, List[str]] = {
    "CPIC": [
        "Clinical Pharmacogenetics Implementation Consortium",
        "CPIC guideline", "level A", "level B",
        "gene-drug pair", "clinical recommendation",
        "dosing guideline", "pharmacogenomics guideline",
        "cpicpgx.org", "activity score",
    ],
    "DPWG": [
        "Dutch Pharmacogenetics Working Group", "KNMP",
        "Royal Dutch Pharmacists Association",
        "European pharmacogenomics", "dosing recommendation",
        "gene-drug interaction", "Dutch guideline",
    ],
    "FDA": [
        "Food and Drug Administration", "drug label",
        "boxed warning", "pharmacogenomic biomarker",
        "table of pharmacogenomic biomarkers",
        "companion diagnostic", "required testing",
        "FDA labeling", "drug approval",
    ],
    "PharmGKB": [
        "Pharmacogenomics Knowledge Base", "clinical annotation",
        "variant annotation", "drug label annotation",
        "pathway", "VIP gene", "pharmgkb.org",
        "Stanford", "evidence level",
    ],
    "PharmVar": [
        "Pharmacogene Variation Consortium", "star allele nomenclature",
        "allele definition", "haplotype", "CYP allele",
        "reference sequence", "suballele", "pharmvar.org",
    ],
}

# ===================================================================
# 6. DOSING_EXPANSION  (10 entries)
# Maps dosing-related terms to PGx concepts.
# ===================================================================

DOSING_EXPANSION: Dict[str, List[str]] = {
    "warfarin dosing": [
        "IWPC", "international warfarin pharmacogenetics consortium",
        "INR", "CYP2C9", "VKORC1", "CYP4F2",
        "weekly dose", "anticoagulation", "dose algorithm",
        "initiation dose", "maintenance dose", "genotype-guided",
    ],
    "tacrolimus dosing": [
        "CYP3A5", "trough level", "transplant", "calcineurin inhibitor",
        "therapeutic drug monitoring", "target range",
        "expresser", "non-expresser", "dose titration",
    ],
    "clopidogrel dosing": [
        "CYP2C19", "loading dose", "maintenance dose",
        "antiplatelet", "prasugrel alternative", "ticagrelor alternative",
        "stent thrombosis prevention", "loss of function",
    ],
    "opioid dosing": [
        "CYP2D6", "codeine", "tramadol", "hydrocodone",
        "morphine equivalent", "equianalgesic",
        "ultra-rapid metabolizer", "poor metabolizer",
        "avoid in UM", "alternative analgesic",
    ],
    "antidepressant dosing": [
        "CYP2D6", "CYP2C19", "SSRI", "TCA", "dose adjustment",
        "titration", "starting dose", "maximum dose",
        "poor metabolizer", "ultra-rapid metabolizer",
        "serotonin syndrome risk",
    ],
    "thiopurine dosing": [
        "TPMT", "NUDT15", "azathioprine", "mercaptopurine",
        "6-MP", "dose reduction", "myelosuppression",
        "start low", "activity score", "genotype-guided",
    ],
    "fluoropyrimidine dosing": [
        "DPYD", "5-fluorouracil", "capecitabine", "DPD deficiency",
        "dose reduction", "contraindicated", "pre-treatment screening",
        "DPYD activity score", "genotype-guided dosing",
    ],
    "statin dosing": [
        "SLCO1B1", "simvastatin", "myopathy risk", "dose limitation",
        "20 mg maximum", "alternative statin", "pravastatin",
        "rosuvastatin", "ABCG2",
    ],
    "PPI dosing": [
        "CYP2C19", "omeprazole", "lansoprazole", "pantoprazole",
        "proton pump inhibitor", "Helicobacter pylori",
        "poor metabolizer", "ultra-rapid metabolizer",
        "increase dose for RM/UM", "eradication therapy",
    ],
    "voriconazole dosing": [
        "CYP2C19", "antifungal", "therapeutic drug monitoring",
        "trough level", "poor metabolizer", "ultra-rapid metabolizer",
        "aspergillosis", "invasive fungal infection",
    ],
}

# ===================================================================
# 7. HLA_EXPANSION  (12 entries)
# Maps HLA alleles to associated drugs, reactions, and populations.
# ===================================================================

HLA_EXPANSION: Dict[str, List[str]] = {
    "HLA-B*57:01": [
        "abacavir", "hypersensitivity", "HIV", "mandatory screening",
        "abacavir hypersensitivity reaction", "fever", "rash",
        "gastrointestinal", "respiratory", "CPIC level A",
    ],
    "HLA-B*15:02": [
        "carbamazepine", "SJS", "TEN", "Southeast Asian",
        "screening", "phenytoin", "oxcarbazepine",
        "Stevens-Johnson syndrome", "toxic epidermal necrolysis",
        "Han Chinese", "Thai", "Filipino", "Malaysian",
    ],
    "HLA-B*58:01": [
        "allopurinol", "SJS", "TEN", "DRESS",
        "gout", "hyperuricemia", "febuxostat alternative",
        "African American", "Southeast Asian", "Korean",
    ],
    "HLA-A*31:01": [
        "carbamazepine", "DRESS", "SJS", "TEN",
        "European", "Japanese", "drug reaction with eosinophilia",
        "anticonvulsant hypersensitivity",
    ],
    "HLA-B*13:01": [
        "dapsone", "hypersensitivity", "leprosy",
        "Southeast Asian", "DRESS",
    ],
    "HLA-B*35:05": [
        "nevirapine", "hypersensitivity", "HIV",
        "skin rash", "Southeast Asian",
    ],
    "HLA-A*02:01": [
        "thalidomide", "erythema nodosum leprosum",
        "multiple myeloma", "immunomodulatory",
    ],
    "HLA-DRB1*07:01": [
        "lapatinib", "hepatotoxicity", "breast cancer",
        "HER2-positive", "liver injury",
    ],
    "HLA-B*44:03": [
        "methazolamide", "SJS", "TEN",
        "carbonic anhydrase inhibitor", "glaucoma",
    ],
    "HLA-A*33:03": [
        "ticlopidine", "hepatotoxicity",
        "antiplatelet", "liver injury",
    ],
    "HLA-DQA1*02:01": [
        "lapatinib", "hepatotoxicity",
        "breast cancer", "tyrosine kinase inhibitor",
    ],
    "HLA-B*53:01": [
        "abacavir", "hypersensitivity risk",
        "African descent", "HIV",
    ],
}

# ===================================================================
# 8. INTERACTION_EXPANSION  (10 entries)
# Maps drug interaction terms to PGx concepts.
# ===================================================================

INTERACTION_EXPANSION: Dict[str, List[str]] = {
    "phenoconversion": [
        "drug-drug-gene interaction", "CYP inhibitor", "CYP inducer",
        "effective phenotype", "genotype-predicted phenotype",
        "inhibitor reclassification", "adjusted phenotype",
        "functional phenotype", "clinical phenotype",
    ],
    "polypharmacy": [
        "multiple medications", "drug interaction", "medication review",
        "elderly", "comorbidity", "deprescribing",
        "pharmacogenomics-guided", "medication reconciliation",
        "drug burden index",
    ],
    "CYP inhibition": [
        "enzyme inhibition", "competitive inhibition",
        "mechanism-based inhibition", "time-dependent inhibition",
        "strong inhibitor", "moderate inhibitor", "weak inhibitor",
        "fluoxetine", "paroxetine", "quinidine", "bupropion",
    ],
    "CYP induction": [
        "enzyme induction", "increased metabolism", "rifampin",
        "carbamazepine", "phenobarbital", "St John's wort",
        "CYP3A4 induction", "pregnane X receptor", "PXR",
    ],
    "drug-drug interaction": [
        "DDI", "concomitant medication", "comedication",
        "interacting drug", "substrate", "inhibitor", "inducer",
        "contraindicated combination", "dose adjustment required",
    ],
    "drug-gene interaction": [
        "DGI", "pharmacogenomic interaction", "gene-drug pair",
        "genotype-guided therapy", "actionable gene-drug",
        "CPIC level A", "clinical annotation",
    ],
    "drug-drug-gene interaction": [
        "DDGI", "phenoconversion", "triple interaction",
        "combined effect", "CYP inhibitor plus genotype",
        "clinical complexity",
    ],
    "substrate": [
        "metabolized by", "CYP substrate", "enzyme substrate",
        "biotransformation", "clearance pathway", "metabolic pathway",
    ],
    "prodrug activation": [
        "prodrug", "bioactivation", "active metabolite",
        "codeine to morphine", "clopidogrel activation",
        "CYP2D6", "CYP2C19", "tamoxifen to endoxifen",
    ],
    "narrow therapeutic index": [
        "NTI", "narrow therapeutic window", "small safety margin",
        "warfarin", "phenytoin", "tacrolimus", "theophylline",
        "therapeutic drug monitoring", "dose-critical",
    ],
}

# ===================================================================
# 9. POPULATION_EXPANSION  (8 entries)
# Maps population/ethnicity terms to PGx-relevant allele frequencies.
# ===================================================================

POPULATION_EXPANSION: Dict[str, List[str]] = {
    "African American": [
        "CYP2D6*17", "CYP3A5*1", "warfarin dose", "allele frequency",
        "African descent", "Black", "CYP2D6*29", "CYP2C9*5", "CYP2C9*8",
        "CYP2C9*11", "lower warfarin dose", "HLA-B*58:01 prevalence",
    ],
    "East Asian": [
        "CYP2D6*10", "CYP2C19*2", "CYP2C19*3", "HLA-B*15:02",
        "Han Chinese", "Japanese", "Korean",
        "NUDT15 R139C", "UGT1A1*6", "allele frequency",
        "carbamazepine screening", "thiopurine sensitivity",
    ],
    "European": [
        "Caucasian", "White", "CYP2D6*4", "CYP2D6*5", "CYP2D6*41",
        "CYP2C9*2", "CYP2C9*3", "DPYD*2A", "SLCO1B1*5",
        "HLA-A*31:01", "allele frequency",
    ],
    "South Asian": [
        "Indian", "Pakistani", "Bangladeshi", "Sri Lankan",
        "CYP2C19*2", "CYP2D6*10", "HLA-B*15:02",
        "allele frequency", "pharmacogenomic diversity",
    ],
    "Latino": [
        "Hispanic", "Latin American", "Mexican American",
        "CYP2D6*4", "CYP2C9*2", "CYP2C19*2",
        "allele frequency", "admixture", "population stratification",
    ],
    "Middle Eastern": [
        "Arab", "North African", "CYP2D6*10",
        "CYP2C19*2", "allele frequency",
        "underrepresented", "pharmacogenomic data gap",
    ],
    "Oceanian": [
        "Pacific Islander", "Polynesian", "Melanesian",
        "CYP2D6*10", "allele frequency",
        "underrepresented", "pharmacogenomic data gap",
    ],
    "Ashkenazi Jewish": [
        "CYP2C19*2", "BRCA1", "BRCA2",
        "founder effect", "allele frequency",
        "genetic testing", "population genetics",
    ],
}

# ===================================================================
# 10. DISEASE_EXPANSION  (15 entries)
# Maps diseases to PGx-relevant drugs and genes.
# ===================================================================

DISEASE_EXPANSION: Dict[str, List[str]] = {
    "depression": [
        "SSRI", "SNRI", "TCA", "CYP2D6", "CYP2C19",
        "antidepressant selection", "fluoxetine", "sertraline",
        "escitalopram", "amitriptyline", "venlafaxine",
        "treatment-resistant", "major depressive disorder",
    ],
    "atrial fibrillation": [
        "warfarin", "CYP2C9", "VKORC1", "anticoagulation",
        "stroke prevention", "INR", "DOAC",
        "propafenone", "CYP2D6", "flecainide",
    ],
    "acute coronary syndrome": [
        "clopidogrel", "CYP2C19", "antiplatelet", "PCI",
        "stent thrombosis", "prasugrel", "ticagrelor",
        "DAPT", "MACE", "myocardial infarction",
    ],
    "epilepsy": [
        "carbamazepine", "phenytoin", "HLA-B*15:02", "CYP2C9",
        "anticonvulsant", "lamotrigine", "oxcarbazepine",
        "SJS risk", "seizure", "antiepileptic drug",
    ],
    "breast cancer": [
        "tamoxifen", "CYP2D6", "endoxifen", "ER-positive",
        "aromatase inhibitor", "trastuzumab", "HER2",
        "adjuvant therapy", "SERM",
    ],
    "colorectal cancer": [
        "5-fluorouracil", "capecitabine", "DPYD", "irinotecan",
        "UGT1A1", "cetuximab", "KRAS", "RAS testing",
        "fluoropyrimidine", "chemotherapy",
    ],
    "HIV": [
        "abacavir", "HLA-B*57:01", "efavirenz", "CYP2B6",
        "antiretroviral", "NRTI", "NNRTI",
        "screening", "hypersensitivity",
    ],
    "organ transplant": [
        "tacrolimus", "CYP3A5", "mycophenolate", "cyclosporine",
        "immunosuppressant", "rejection", "kidney transplant",
        "liver transplant", "therapeutic drug monitoring",
    ],
    "gout": [
        "allopurinol", "HLA-B*58:01", "febuxostat",
        "hyperuricemia", "uric acid", "xanthine oxidase",
        "SJS risk", "DRESS risk",
    ],
    "IBD": [
        "inflammatory bowel disease", "Crohn disease", "ulcerative colitis",
        "azathioprine", "mercaptopurine", "TPMT", "NUDT15",
        "thiopurine", "myelosuppression", "sulfasalazine", "NAT2",
    ],
    "pain": [
        "opioid", "codeine", "tramadol", "CYP2D6",
        "analgesic", "neuropathic pain", "acute pain",
        "chronic pain", "oxycodone", "hydrocodone",
        "ultra-rapid metabolizer", "poor metabolizer",
    ],
    "schizophrenia": [
        "clozapine", "CYP1A2", "aripiprazole", "CYP2D6",
        "haloperidol", "risperidone", "antipsychotic",
        "treatment-resistant", "metabolic side effects",
    ],
    "GERD": [
        "gastroesophageal reflux disease", "omeprazole", "lansoprazole",
        "CYP2C19", "proton pump inhibitor", "PPI",
        "Helicobacter pylori", "acid reflux",
    ],
    "cystic fibrosis": [
        "CFTR", "ivacaftor", "lumacaftor", "elexacaftor", "tezacaftor",
        "F508del", "G551D", "CFTR modulator", "sweat chloride",
        "companion diagnostic", "precision medicine",
    ],
    "malaria": [
        "primaquine", "G6PD", "chloroquine", "hemolytic anemia",
        "Plasmodium vivax", "radical cure", "screening",
    ],
}

# ===================================================================
# 11. SAFETY_EXPANSION  (8 entries)
# Maps safety terms to PGx concepts.
# ===================================================================

SAFETY_EXPANSION: Dict[str, List[str]] = {
    "adverse drug reaction": [
        "ADR", "pharmacogenomics", "genetic testing", "preventable",
        "side effect", "drug safety", "pharmacovigilance",
        "post-market surveillance", "serious adverse event", "SAE",
    ],
    "drug allergy": [
        "HLA", "hypersensitivity", "immune-mediated",
        "SJS", "TEN", "DRESS", "anaphylaxis",
        "drug-induced", "cross-reactivity",
    ],
    "black box warning": [
        "boxed warning", "FDA", "serious risk", "death",
        "life-threatening", "clopidogrel CYP2C19",
        "codeine pediatric", "pharmacogenomic biomarker",
    ],
    "medication error": [
        "prescribing error", "dosing error", "wrong dose",
        "pharmacogenomics", "clinical decision support",
        "preventable adverse event", "patient safety",
    ],
    "drug recall": [
        "safety signal", "post-market", "FDA",
        "pharmacovigilance", "risk-benefit", "REMS",
    ],
    "pediatric safety": [
        "children", "neonatal", "codeine", "CYP2D6",
        "ultra-rapid metabolizer", "respiratory depression",
        "FDA contraindication", "tonsillectomy", "adenoidectomy",
    ],
    "geriatric safety": [
        "elderly", "polypharmacy", "renal function",
        "hepatic function", "age-related metabolism",
        "drug accumulation", "falls risk", "deprescribing",
    ],
    "pregnancy safety": [
        "teratogenicity", "pregnancy category", "lactation",
        "maternal pharmacogenomics", "fetal exposure",
        "dose adjustment in pregnancy", "CYP enzyme changes",
    ],
}

# ===================================================================
# 12. CLINICAL_EXPANSION  (8 entries)
# Maps clinical scenario terms to PGx concepts.
# ===================================================================

CLINICAL_EXPANSION: Dict[str, List[str]] = {
    "pre-emptive testing": [
        "PGx panel", "pharmacogenomic", "genotype", "lifetime",
        "pre-prescription", "prospective testing", "multi-gene panel",
        "preemptive pharmacogenomics", "anticipatory guidance",
        "once-in-a-lifetime test", "combinatorial PGx",
    ],
    "therapeutic drug monitoring": [
        "TDM", "drug levels", "dose adjustment", "trough level",
        "peak level", "steady state", "pharmacokinetics",
        "tacrolimus", "voriconazole", "phenytoin",
    ],
    "genotype-guided therapy": [
        "precision medicine", "personalized medicine",
        "pharmacogenomics-informed", "individualized dosing",
        "gene-guided prescribing", "tailored therapy",
        "right drug right dose", "companion diagnostic",
    ],
    "reactive testing": [
        "post-adverse event", "single-gene test",
        "indication-based testing", "triggered by ADR",
        "after treatment failure", "diagnostic testing",
    ],
    "clinical validity": [
        "analytical validity", "clinical utility",
        "positive predictive value", "negative predictive value",
        "sensitivity", "specificity", "evidence-based",
    ],
    "point-of-care testing": [
        "rapid genotyping", "bedside testing", "POCT",
        "turnaround time", "CYP2C19 rapid test",
        "clopidogrel", "cardiac catheterization", "Spartan RX",
    ],
    "pharmacogenomic consultation": [
        "PGx consult", "clinical pharmacist", "genetic counselor",
        "interpretation", "medication review",
        "genotype-phenotype correlation", "clinical recommendation",
    ],
    "treatment failure": [
        "non-response", "therapeutic failure", "inefficacy",
        "subtherapeutic", "ultra-rapid metabolizer",
        "prodrug activation", "CYP2C19 clopidogrel",
        "CYP2D6 codeine", "treatment resistance",
    ],
}

# ===================================================================
# 13. IMPLEMENTATION_EXPANSION  (6 entries)
# Maps PGx implementation terms to related concepts.
# ===================================================================

IMPLEMENTATION_EXPANSION: Dict[str, List[str]] = {
    "clinical decision support": [
        "CDS", "EHR integration", "alert", "prescribing",
        "electronic health record", "best practice alert",
        "interruptive alert", "informational alert",
        "CPIC guideline integration", "clinical workflow",
    ],
    "IGNITE": [
        "implementation network", "health system", "adoption",
        "Implementing Genomics in Practice",
        "NIH", "pragmatic clinical trial", "dissemination",
    ],
    "eMERGE": [
        "electronic Medical Records and Genomics",
        "biobank", "genotype-phenotype", "EHR-linked",
        "return of results", "clinical sequencing",
    ],
    "pharmacist role": [
        "clinical pharmacist", "pharmacogenomics service",
        "medication therapy management", "MTM",
        "drug interaction review", "dose optimization",
        "pharmacy practice", "PGx champion",
    ],
    "cost-effectiveness": [
        "economic evaluation", "cost-benefit analysis",
        "QALY", "ICER", "incremental cost-effectiveness ratio",
        "value-based", "reimbursement", "payer coverage",
        "return on investment", "health economics",
    ],
    "education and training": [
        "pharmacogenomics education", "curriculum",
        "competency", "workforce development",
        "continuing education", "CME", "pharmacy school",
        "medical school", "genetic literacy",
    ],
}

# ===================================================================
# 14. REGULATORY_EXPANSION  (5 entries)
# Maps regulatory terms to PGx-related concepts.
# ===================================================================

REGULATORY_EXPANSION: Dict[str, List[str]] = {
    "FDA labeling": [
        "boxed warning", "pharmacogenomic biomarker", "drug label",
        "table of pharmacogenomic biomarkers in drug labeling",
        "required testing", "recommended testing", "actionable PGx",
        "informative PGx", "label change",
    ],
    "CLIA": [
        "Clinical Laboratory Improvement Amendments",
        "clinical laboratory", "genetic testing", "certification",
        "laboratory-developed test", "LDT", "quality assurance",
        "proficiency testing", "accreditation",
    ],
    "companion diagnostic": [
        "CDx", "FDA-approved diagnostic", "co-development",
        "HER2 testing", "KRAS testing", "CFTR testing",
        "trastuzumab", "ivacaftor", "biomarker-driven therapy",
    ],
    "genetic non-discrimination": [
        "GINA", "Genetic Information Nondiscrimination Act",
        "insurance coverage", "employment protection",
        "privacy", "genetic privacy", "patient consent",
        "informed consent", "ethical considerations",
    ],
    "direct-to-consumer": [
        "DTC", "consumer genomics", "23andMe", "pharmacogenomics panel",
        "clinical confirmation", "clinical-grade testing",
        "limitations", "FDA authorization",
    ],
}


# ===================================================================
# ALL_EXPANSION_MAPS: ordered list of (name, dict) tuples
# ===================================================================

ALL_EXPANSION_MAPS: List[Tuple[str, Dict[str, List[str]]]] = [
    ("DRUG_EXPANSION", DRUG_EXPANSION),
    ("GENE_EXPANSION", GENE_EXPANSION),
    ("PHENOTYPE_EXPANSION", PHENOTYPE_EXPANSION),
    ("TOXICITY_EXPANSION", TOXICITY_EXPANSION),
    ("GUIDELINE_EXPANSION", GUIDELINE_EXPANSION),
    ("DOSING_EXPANSION", DOSING_EXPANSION),
    ("HLA_EXPANSION", HLA_EXPANSION),
    ("INTERACTION_EXPANSION", INTERACTION_EXPANSION),
    ("POPULATION_EXPANSION", POPULATION_EXPANSION),
    ("DISEASE_EXPANSION", DISEASE_EXPANSION),
    ("SAFETY_EXPANSION", SAFETY_EXPANSION),
    ("CLINICAL_EXPANSION", CLINICAL_EXPANSION),
    ("IMPLEMENTATION_EXPANSION", IMPLEMENTATION_EXPANSION),
    ("REGULATORY_EXPANSION", REGULATORY_EXPANSION),
]


# ===================================================================
# FUNCTIONS
# ===================================================================


def expand_query(query: str) -> List[str]:
    """Expand a pharmacogenomics query with domain-specific synonyms.

    Scans the query (case-insensitive) against all 14 expansion maps
    and returns a deduplicated, sorted list of expansion terms that
    matched any key in any map.

    Args:
        query: Natural-language query string.

    Returns:
        Sorted list of unique expansion terms found across all maps.
        Returns an empty list if no keys match.

    Example:
        >>> terms = expand_query("warfarin dosing for CYP2C9 poor metabolizer")
        >>> "VKORC1" in terms
        True
        >>> "dose reduction" in terms
        True
    """
    query_lower = query.lower()
    expanded: set[str] = set()

    for _name, expansion_map in ALL_EXPANSION_MAPS:
        for key, terms in expansion_map.items():
            if key.lower() in query_lower:
                expanded.update(terms)

    return sorted(expanded)


def expand_query_by_category(query: str) -> Dict[str, List[str]]:
    """Expand a query and group results by expansion-map category.

    Args:
        query: Natural-language query string.

    Returns:
        Dictionary mapping category names (e.g., ``"DRUG_EXPANSION"``)
        to lists of matched expansion terms.  Categories with no
        matches are omitted.

    Example:
        >>> result = expand_query_by_category("carbamazepine SJS Asian")
        >>> "DRUG_EXPANSION" in result
        True
        >>> "HLA_EXPANSION" in result
        True
    """
    query_lower = query.lower()
    categorized: Dict[str, List[str]] = {}

    for name, expansion_map in ALL_EXPANSION_MAPS:
        matched_terms: set[str] = set()
        for key, terms in expansion_map.items():
            if key.lower() in query_lower:
                matched_terms.update(terms)
        if matched_terms:
            categorized[name] = sorted(matched_terms)

    return categorized


def get_expansion_stats() -> Dict[str, object]:
    """Return statistics about the expansion maps.

    Returns:
        Dictionary with keys:
            - ``total_maps``: Number of expansion maps (14).
            - ``total_keys``: Total number of keys across all maps.
            - ``total_terms``: Total number of expansion terms (with duplicates).
            - ``unique_terms``: Number of globally unique expansion terms.
            - ``per_map``: Dict mapping each map name to its key count and
              term count.

    Example:
        >>> stats = get_expansion_stats()
        >>> stats["total_maps"]
        14
    """
    total_keys = 0
    total_terms = 0
    all_terms: set[str] = set()
    per_map: Dict[str, Dict[str, int]] = {}

    for name, expansion_map in ALL_EXPANSION_MAPS:
        map_keys = len(expansion_map)
        map_terms = sum(len(v) for v in expansion_map.values())
        total_keys += map_keys
        total_terms += map_terms
        for terms in expansion_map.values():
            all_terms.update(terms)
        per_map[name] = {"keys": map_keys, "terms": map_terms}

    return {
        "total_maps": len(ALL_EXPANSION_MAPS),
        "total_keys": total_keys,
        "total_terms": total_terms,
        "unique_terms": len(all_terms),
        "per_map": per_map,
    }
