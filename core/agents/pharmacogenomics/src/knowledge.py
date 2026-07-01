"""Pharmacogenomics Intelligence Agent — Knowledge Graph.

Extends the Clinker pattern from core/engines/precision-intelligence/src/knowledge.py,
adapted for the pharmacogenomics domain. Contains:

1. PHARMACOGENES: 25 gene entries with clinical pharmacogenomic data
2. METABOLIZER_PHENOTYPES: 5 phenotype classifications
3. DRUG_CATEGORIES: 12 therapeutic categories with member drugs
4. CYP_INHIBITORS: 4 enzymes with strong/moderate/weak inhibitors
5. CYP_INDUCERS: 3 enzymes with strong/moderate inducers
6. HLA_DRUG_ASSOCIATIONS: 12 HLA-drug hypersensitivity pairs
7. DRUG_ALTERNATIVES: 30+ gene-phenotype-specific drug substitutions
8. ACTIVITY_SCORE_TABLES: CYP2D6 and DPYD activity score mappings
9. ENTITY_ALIASES: 80+ aliases for drugs, genes, and phenotypes

Author: Adam Jones
Date: March 2026
"""

from typing import Any, Dict, List, Optional


# =============================================================================
# 1. PHARMACOGENES — 25 pharmacogene entries with clinical data
# =============================================================================

PHARMACOGENES: Dict[str, Dict[str, Any]] = {
    "CYP2D6": {
        "full_name": "Cytochrome P450 2D6",
        "chromosome": "22q13.2",
        "function": "Phase I oxidative metabolism of ~25% of clinically used drugs",
        "substrates_count": 120,
        "percent_drugs_metabolized": 25,
        "star_alleles_defined": 150,
        "key_variants": [
            "*1 (normal function)", "*2 (normal function)", "*3 (no function — frameshift)",
            "*4 (no function — splicing defect)", "*5 (gene deletion)",
            "*6 (no function — frameshift)", "*9 (decreased function)",
            "*10 (decreased function)", "*17 (decreased function)",
            "*29 (decreased function)", "*41 (decreased function)",
            "*1xN/*2xN (gene duplication — ultra-rapid)",
        ],
        "structural_variation": True,
        "complexity_level": "very high — gene deletions, duplications, tandems, hybrids with CYP2D7",
        "cpic_guidelines": [
            "codeine", "tramadol", "tamoxifen", "ondansetron",
            "tropisetron", "atomoxetine", "venlafaxine", "eliglustat",
        ],
    },
    "CYP2C19": {
        "full_name": "Cytochrome P450 2C19",
        "chromosome": "10q23.33",
        "function": "Phase I metabolism of PPIs, antidepressants, antiplatelet agents",
        "substrates_count": 80,
        "percent_drugs_metabolized": 10,
        "star_alleles_defined": 38,
        "key_variants": [
            "*1 (normal function)", "*2 (no function — splicing defect, c.681G>A)",
            "*3 (no function — premature stop, c.636G>A)",
            "*4 (no function)", "*17 (increased function — c.-806C>T)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate",
        "cpic_guidelines": [
            "clopidogrel", "voriconazole", "escitalopram", "citalopram",
            "sertraline", "amitriptyline", "clomipramine", "doxepin",
            "imipramine", "trimipramine",
        ],
    },
    "CYP2C9": {
        "full_name": "Cytochrome P450 2C9",
        "chromosome": "10q23.33",
        "function": "Phase I metabolism of warfarin, NSAIDs, sulfonylureas, phenytoin",
        "substrates_count": 100,
        "percent_drugs_metabolized": 15,
        "star_alleles_defined": 75,
        "key_variants": [
            "*1 (normal function)", "*2 (decreased function — R144C, rs1799853)",
            "*3 (decreased function — I359L, rs1057910)",
            "*5 (decreased function)", "*6 (no function)",
            "*8 (decreased function)", "*11 (decreased function)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate",
        "cpic_guidelines": [
            "warfarin", "phenytoin", "celecoxib", "flurbiprofen",
            "ibuprofen", "lornoxicam", "meloxicam", "piroxicam",
            "tenoxicam", "siponimod",
        ],
    },
    "CYP3A4": {
        "full_name": "Cytochrome P450 3A4",
        "chromosome": "7q22.1",
        "function": "Most abundant hepatic CYP; Phase I metabolism of ~30% of all drugs",
        "substrates_count": 250,
        "percent_drugs_metabolized": 30,
        "star_alleles_defined": 34,
        "key_variants": [
            "*1 (normal function)", "*1B (possible increased expression)",
            "*22 (decreased function — intron 6 variant, rs35599367)",
        ],
        "structural_variation": False,
        "complexity_level": "low — few clinically actionable variants; phenotype driven by DDI",
        "cpic_guidelines": [
            "tacrolimus",
        ],
    },
    "CYP3A5": {
        "full_name": "Cytochrome P450 3A5",
        "chromosome": "7q22.1",
        "function": "Phase I metabolism; clinically relevant for tacrolimus dosing",
        "substrates_count": 50,
        "percent_drugs_metabolized": 5,
        "star_alleles_defined": 11,
        "key_variants": [
            "*1 (normal function — expresser allele)",
            "*3 (no function — splice defect, rs776746, most common loss allele)",
            "*6 (no function)", "*7 (no function)",
        ],
        "structural_variation": False,
        "complexity_level": "low",
        "cpic_guidelines": [
            "tacrolimus",
        ],
    },
    "CYP1A2": {
        "full_name": "Cytochrome P450 1A2",
        "chromosome": "15q24.1",
        "function": "Phase I metabolism of caffeine, theophylline, clozapine, olanzapine",
        "substrates_count": 60,
        "percent_drugs_metabolized": 8,
        "star_alleles_defined": 21,
        "key_variants": [
            "*1 (normal function)", "*1C (decreased inducibility)",
            "*1F (increased inducibility — rs762551 A>C)",
            "*1K (decreased function)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — phenotype heavily influenced by inducers (smoking, diet)",
        "cpic_guidelines": [],
    },
    "CYP2B6": {
        "full_name": "Cytochrome P450 2B6",
        "chromosome": "19q13.2",
        "function": "Phase I metabolism of efavirenz, methadone, bupropion, cyclophosphamide",
        "substrates_count": 40,
        "percent_drugs_metabolized": 4,
        "star_alleles_defined": 38,
        "key_variants": [
            "*1 (normal function)", "*4 (increased function)",
            "*6 (decreased function — Q172H + K262R, rs3745274 + rs2279343)",
            "*9 (decreased function)", "*18 (no function)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate",
        "cpic_guidelines": [
            "efavirenz",
        ],
    },
    "CYP4F2": {
        "full_name": "Cytochrome P450 4F2",
        "chromosome": "19p13.12",
        "function": "Vitamin K1 oxidase; modulates warfarin and phylloquinone metabolism",
        "substrates_count": 10,
        "percent_drugs_metabolized": 1,
        "star_alleles_defined": 4,
        "key_variants": [
            "*1 (normal function)",
            "*3 (decreased function — V433M, rs2108622)",
        ],
        "structural_variation": False,
        "complexity_level": "low",
        "cpic_guidelines": [
            "warfarin",
        ],
    },
    "VKORC1": {
        "full_name": "Vitamin K Epoxide Reductase Complex Subunit 1",
        "chromosome": "16p11.2",
        "function": "Target of warfarin; recycles vitamin K for clotting factor activation",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "-1639G>A (rs9923231) — increased sensitivity to warfarin",
            "D36Y — warfarin resistance",
            "V66M — warfarin resistance",
            "R58G — incomplete warfarin resistance",
        ],
        "structural_variation": False,
        "complexity_level": "low — key SNP -1639G>A explains most variability",
        "cpic_guidelines": [
            "warfarin",
        ],
    },
    "SLCO1B1": {
        "full_name": "Solute Carrier Organic Anion Transporter Family Member 1B1",
        "chromosome": "12p12.2",
        "function": "Hepatic uptake transporter for statins, bilirubin, methotrexate",
        "substrates_count": 30,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": 46,
        "key_variants": [
            "*1 (normal function)",
            "*5 (decreased function — V174A, rs4149056, c.521T>C)",
            "*14 (decreased function)",
            "*15 (*1B + *5 haplotype — decreased function)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate",
        "cpic_guidelines": [
            "simvastatin", "atorvastatin", "rosuvastatin", "pravastatin",
            "lovastatin", "fluvastatin", "pitavastatin",
        ],
    },
    "DPYD": {
        "full_name": "Dihydropyrimidine Dehydrogenase",
        "chromosome": "1p21.3",
        "function": "Rate-limiting enzyme in pyrimidine catabolism; inactivates 5-fluorouracil",
        "substrates_count": 5,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": 25,
        "key_variants": [
            "*1 (normal function)",
            "*2A (no function — IVS14+1G>A, rs3918290, exon 14 skipping)",
            "*13 (no function — T1679G, rs55886062, I560S)",
            "c.2846A>T (decreased function — rs67376798, D949V)",
            "c.1129-5923C>G (decreased function — rs75017182, HapB3)",
        ],
        "structural_variation": False,
        "complexity_level": "high — 80% of toxicity from 4 key variants",
        "cpic_guidelines": [
            "fluorouracil", "capecitabine", "tegafur",
        ],
    },
    "TPMT": {
        "full_name": "Thiopurine S-Methyltransferase",
        "chromosome": "6p22.3",
        "function": "Methylation (inactivation) of thiopurine drugs; prevents toxic metabolite accumulation",
        "substrates_count": 5,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": 44,
        "key_variants": [
            "*1 (normal function)",
            "*2 (no function — A80P, rs1800462)",
            "*3A (no function — A154T + Y240C, rs1800460 + rs1142345)",
            "*3B (no function — A154T, rs1800460)",
            "*3C (no function — Y240C, rs1142345)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — 3 alleles account for >95% of low activity",
        "cpic_guidelines": [
            "azathioprine", "mercaptopurine", "thioguanine",
        ],
    },
    "NUDT15": {
        "full_name": "Nudix Hydrolase 15",
        "chromosome": "13q14.11",
        "function": "Dephosphorylation of thiopurine active metabolites (TGTP/TdGTP)",
        "substrates_count": 3,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": 18,
        "key_variants": [
            "*1 (normal function)",
            "*2 (no function — p.R139C, rs116855232)",
            "*3 (no function — p.R139H, rs147390019)",
            "*5 (uncertain function)",
            "*6 (no function)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — especially important in East Asian populations",
        "cpic_guidelines": [
            "azathioprine", "mercaptopurine", "thioguanine",
        ],
    },
    "UGT1A1": {
        "full_name": "UDP-Glucuronosyltransferase 1A1",
        "chromosome": "2q37.1",
        "function": "Phase II glucuronidation of bilirubin, SN-38 (irinotecan active metabolite)",
        "substrates_count": 15,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": 113,
        "key_variants": [
            "*1 (normal function — 6 TA repeats in promoter)",
            "*28 (decreased function — 7 TA repeats, rs8175347, Gilbert syndrome)",
            "*6 (decreased function — G71R, rs4148323, common in East Asians)",
            "*36 (increased function — 5 TA repeats)",
            "*37 (decreased function — 8 TA repeats)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — TA repeat polymorphism in promoter",
        "cpic_guidelines": [
            "irinotecan", "atazanavir",
        ],
    },
    "G6PD": {
        "full_name": "Glucose-6-Phosphate Dehydrogenase",
        "chromosome": "Xq28",
        "function": "Rate-limiting enzyme in pentose phosphate pathway; protects RBCs from oxidative damage",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "B (normal — reference allele)",
            "A+ (normal — Val68Met, Asn126Asp)",
            "A- (deficient — Val68Met, Asn126Asp, Leu323Pro, rs1050828+rs1050829)",
            "Mediterranean (deficient — Ser188Phe, rs5030868)",
            "Canton (deficient — Arg459Leu)",
            "Mahidol (deficient — Gly163Ser)",
        ],
        "structural_variation": False,
        "complexity_level": "high — X-linked; >200 known variants; WHO classification system",
        "cpic_guidelines": [
            "rasburicase", "dapsone", "primaquine", "tafenoquine",
        ],
    },
    "NAT2": {
        "full_name": "N-Acetyltransferase 2",
        "chromosome": "8p22",
        "function": "Phase II acetylation of isoniazid, hydralazine, sulfonamides, procainamide",
        "substrates_count": 15,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": 88,
        "key_variants": [
            "*4 (rapid acetylator — reference)",
            "*5 (slow acetylator — I114T, rs1801280)",
            "*6 (slow acetylator — R197Q, rs1799930)",
            "*7 (slow acetylator — G286E, rs1799931)",
            "*13 (rapid acetylator)",
            "*14 (slow acetylator — R64Q, rs1801279)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — bimodal/trimodal acetylator phenotype",
        "cpic_guidelines": [
            "isoniazid",
        ],
    },
    "ABCB1": {
        "full_name": "ATP-Binding Cassette Subfamily B Member 1 (P-glycoprotein/MDR1)",
        "chromosome": "7q21.12",
        "function": "Efflux transporter at blood-brain barrier, gut, liver, kidney; limits drug absorption and CNS penetration",
        "substrates_count": 100,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "C3435T (rs1045642 — synonymous but affects mRNA stability/folding)",
            "G2677T/A (rs2032582 — A893S/T, nonsynonymous)",
            "C1236T (rs1128503 — synonymous, LD with 3435T)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — haplotype effects; inconsistent clinical data",
        "cpic_guidelines": [],
    },
    "CYP2C8": {
        "full_name": "Cytochrome P450 2C8",
        "chromosome": "10q23.33",
        "function": "Phase I metabolism of paclitaxel, amodiaquine, repaglinide, rosiglitazone",
        "substrates_count": 30,
        "percent_drugs_metabolized": 5,
        "star_alleles_defined": 14,
        "key_variants": [
            "*1 (normal function)",
            "*2 (decreased function — I269F, rs11572103)",
            "*3 (decreased function — R139K, rs10509681)",
            "*4 (decreased function — I264M, rs1058930)",
        ],
        "structural_variation": False,
        "complexity_level": "low",
        "cpic_guidelines": [],
    },
    "IFNL3": {
        "full_name": "Interferon Lambda 3 (IL28B)",
        "chromosome": "19q13.2",
        "function": "Innate immune response to viral infection; predicts HCV treatment response",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "rs12979860 C/T (CC = favorable response to interferon-based HCV therapy)",
            "rs8099917 T/G (TT = favorable response)",
        ],
        "structural_variation": False,
        "complexity_level": "low — single SNP drives most predictive value",
        "cpic_guidelines": [
            "peginterferon alfa-2a", "peginterferon alfa-2b",
        ],
    },
    "RYR1": {
        "full_name": "Ryanodine Receptor 1",
        "chromosome": "19q13.2",
        "function": "Skeletal muscle calcium release channel; mutations cause malignant hyperthermia susceptibility",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "c.1021G>A (p.G341R — MH causative)",
            "c.1840C>T (p.R614C — MH causative, most common)",
            "c.7300G>A (p.G2434R — MH causative)",
            "c.7354C>T (p.R2452W — MH causative)",
            "c.6502G>A (p.V2168M — MH causative)",
        ],
        "structural_variation": False,
        "complexity_level": "high — large gene, >400 variants, incomplete penetrance",
        "cpic_guidelines": [
            "volatile anesthetics (desflurane, sevoflurane, isoflurane, halothane)",
            "succinylcholine",
        ],
    },
    "CACNA1S": {
        "full_name": "Calcium Voltage-Gated Channel Subunit Alpha1 S",
        "chromosome": "1q32.1",
        "function": "L-type calcium channel in skeletal muscle; mutations cause malignant hyperthermia susceptibility",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "c.520C>T (p.R174W — MH causative)",
            "c.3257G>A (p.R1086H — MH susceptible)",
        ],
        "structural_variation": False,
        "complexity_level": "moderate — fewer known variants than RYR1",
        "cpic_guidelines": [
            "volatile anesthetics", "succinylcholine",
        ],
    },
    "CFTR": {
        "full_name": "Cystic Fibrosis Transmembrane Conductance Regulator",
        "chromosome": "7q31.2",
        "function": "Chloride/bicarbonate channel; mutations cause cystic fibrosis and guide CFTR modulator therapy",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "F508del (most common — protein misfolding, class II)",
            "G551D (gating defect, class III — ivacaftor responsive)",
            "G542X (nonsense, class I)", "W1282X (nonsense, class I)",
            "N1303K (class II)", "R117H (variable severity, class IV)",
        ],
        "structural_variation": False,
        "complexity_level": "high — >2000 variants identified, 6 functional classes",
        "cpic_guidelines": [
            "ivacaftor", "lumacaftor/ivacaftor", "tezacaftor/ivacaftor",
            "elexacaftor/tezacaftor/ivacaftor",
        ],
    },
    "F5": {
        "full_name": "Coagulation Factor V",
        "chromosome": "1q24.2",
        "function": "Coagulation cascade cofactor; Factor V Leiden increases VTE risk, affects anticoagulant selection",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "Factor V Leiden (R506Q, rs6025 — activated protein C resistance)",
            "Factor V Cambridge (R306T)", "Factor V Hong Kong (R306G)",
        ],
        "structural_variation": False,
        "complexity_level": "low — single major variant (Leiden)",
        "cpic_guidelines": [
            "hormonal contraceptives",
        ],
    },
    "MTHFR": {
        "full_name": "Methylenetetrahydrofolate Reductase",
        "chromosome": "1p36.22",
        "function": "Folate metabolism; converts 5,10-methyleneTHF to 5-methylTHF for homocysteine remethylation",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "C677T (A222V, rs1801133 — thermolabile variant, 30% reduced activity per T allele)",
            "A1298C (E429A, rs1801131 — mildly reduced activity)",
        ],
        "structural_variation": False,
        "complexity_level": "low — two common SNPs, but clinical actionability debated",
        "cpic_guidelines": [],
    },
    "HLA": {
        "full_name": "Human Leukocyte Antigen Complex (umbrella)",
        "chromosome": "6p21.3",
        "function": "Immune recognition; specific alleles predict severe drug hypersensitivity reactions (SJS/TEN, DRESS, DILI)",
        "substrates_count": None,
        "percent_drugs_metabolized": None,
        "star_alleles_defined": None,
        "key_variants": [
            "HLA-B*57:01 (abacavir hypersensitivity)",
            "HLA-B*15:02 (carbamazepine/phenytoin SJS/TEN in SE Asians)",
            "HLA-B*58:01 (allopurinol SJS/TEN)",
            "HLA-A*31:01 (carbamazepine DRESS/SJS in Europeans/Japanese)",
            "HLA-B*15:11 (carbamazepine SJS in Japanese/Korean)",
            "HLA-B*57:01 (flucloxacillin DILI)",
        ],
        "structural_variation": False,
        "complexity_level": "very high — most polymorphic region in human genome; >30,000 alleles across loci",
        "cpic_guidelines": [
            "abacavir", "carbamazepine", "oxcarbazepine", "phenytoin",
            "allopurinol",
        ],
    },
}


# =============================================================================
# 2. METABOLIZER_PHENOTYPES — 5 standard classifications
# =============================================================================

METABOLIZER_PHENOTYPES: Dict[str, Dict[str, str]] = {
    "ultra_rapid": {
        "abbreviation": "UM",
        "clinical_meaning": (
            "Metabolizes the substrate drug much faster than normal. For prodrugs "
            "(e.g., codeine → morphine), this leads to dangerously high active "
            "metabolite concentrations. For active drugs, this causes therapeutic "
            "failure due to rapid clearance."
        ),
        "risk": (
            "Prodrugs: toxicity/overdose from excessive active metabolite formation. "
            "Active drugs: subtherapeutic levels, treatment failure. "
            "Example: CYP2D6 UM + codeine → morphine toxicity, respiratory depression, death."
        ),
    },
    "rapid": {
        "abbreviation": "RM",
        "clinical_meaning": (
            "Metabolizes the substrate drug faster than normal but less than ultra-rapid. "
            "Assigned when activity score falls between normal and ultra-rapid thresholds. "
            "Clinical significance varies by drug."
        ),
        "risk": (
            "Prodrugs: somewhat increased active metabolite exposure. "
            "Active drugs: may need dose increases. "
            "Clinical impact generally milder than UM but monitoring recommended."
        ),
    },
    "normal": {
        "abbreviation": "NM",
        "clinical_meaning": (
            "Metabolizes the substrate drug at the expected rate. Standard dosing "
            "recommendations apply. Most individuals fall into this category."
        ),
        "risk": (
            "Standard risk. Drug dosing per label. No pharmacogenomic dose adjustment needed."
        ),
    },
    "intermediate": {
        "abbreviation": "IM",
        "clinical_meaning": (
            "Metabolizes the substrate drug more slowly than normal. One functional "
            "and one non-functional allele, or two decreased-function alleles. "
            "Drug clearance is reduced, leading to higher plasma concentrations."
        ),
        "risk": (
            "Active drugs: increased plasma levels, higher risk of ADRs at standard doses. "
            "Prodrugs: reduced conversion to active metabolite, possible decreased efficacy. "
            "Dose reduction or alternative drug often recommended."
        ),
    },
    "poor": {
        "abbreviation": "PM",
        "clinical_meaning": (
            "Minimal to no functional enzyme activity. Two non-functional alleles. "
            "Drug clearance dramatically reduced for substrates. "
            "Markedly increased risk of ADRs at standard doses."
        ),
        "risk": (
            "Active drugs: significantly elevated plasma levels, high ADR risk, "
            "contraindication in some cases. Prodrugs: negligible conversion to active "
            "metabolite, therapeutic failure. "
            "Alternative drug or major dose reduction typically required."
        ),
    },
}


# =============================================================================
# 3. DRUG_CATEGORIES — 12 therapeutic categories
# =============================================================================

DRUG_CATEGORIES: Dict[str, Dict[str, Any]] = {
    "opioids": {
        "description": "Opioid analgesics with pharmacogenomic relevance (CYP2D6-mediated metabolism)",
        "drugs": [
            "codeine", "tramadol", "hydrocodone", "oxycodone",
            "morphine", "fentanyl", "methadone", "buprenorphine",
        ],
        "primary_genes": ["CYP2D6", "CYP3A4", "CYP2B6"],
    },
    "anticoagulants": {
        "description": "Anticoagulant and antiplatelet agents with pharmacogenomic dosing",
        "drugs": [
            "warfarin", "clopidogrel", "heparin", "enoxaparin",
            "rivaroxaban", "apixaban", "dabigatran", "edoxaban",
        ],
        "primary_genes": ["CYP2C9", "VKORC1", "CYP4F2", "CYP2C19"],
    },
    "antidepressants": {
        "description": "SSRIs, SNRIs, and tricyclics with CYP2D6/CYP2C19-guided dosing",
        "drugs": [
            "fluoxetine", "paroxetine", "sertraline", "citalopram",
            "escitalopram", "venlafaxine", "duloxetine", "amitriptyline",
            "nortriptyline", "clomipramine", "doxepin", "imipramine",
            "desipramine", "trimipramine", "fluvoxamine", "mirtazapine",
        ],
        "primary_genes": ["CYP2D6", "CYP2C19"],
    },
    "antipsychotics": {
        "description": "Antipsychotic medications metabolized by CYP2D6 and CYP1A2",
        "drugs": [
            "haloperidol", "pimozide", "aripiprazole", "brexpiprazole",
            "clozapine", "olanzapine", "risperidone", "quetiapine",
            "iloperidone", "thioridazine",
        ],
        "primary_genes": ["CYP2D6", "CYP1A2", "CYP3A4"],
    },
    "statins": {
        "description": "HMG-CoA reductase inhibitors with SLCO1B1 myopathy risk",
        "drugs": [
            "simvastatin", "atorvastatin", "rosuvastatin", "pravastatin",
            "lovastatin", "fluvastatin", "pitavastatin",
        ],
        "primary_genes": ["SLCO1B1"],
    },
    "chemotherapy": {
        "description": "Chemotherapy agents with pharmacogenomic toxicity prediction",
        "drugs": [
            "fluorouracil", "capecitabine", "tegafur", "irinotecan",
            "mercaptopurine", "azathioprine", "thioguanine",
            "tamoxifen", "cyclophosphamide", "cisplatin", "carboplatin",
        ],
        "primary_genes": ["DPYD", "UGT1A1", "TPMT", "NUDT15", "CYP2D6"],
    },
    "anticonvulsants": {
        "description": "Antiepileptic drugs with HLA and CYP-based safety/dosing",
        "drugs": [
            "carbamazepine", "oxcarbazepine", "phenytoin", "fosphenytoin",
            "lamotrigine", "valproic acid", "levetiracetam", "phenobarbital",
        ],
        "primary_genes": ["HLA-B", "HLA-A", "CYP2C9", "CYP2C19"],
    },
    "antivirals": {
        "description": "Antiviral agents with pharmacogenomic relevance",
        "drugs": [
            "abacavir", "efavirenz", "nevirapine", "atazanavir",
            "peginterferon alfa-2a", "peginterferon alfa-2b", "ribavirin",
        ],
        "primary_genes": ["HLA-B", "CYP2B6", "UGT1A1", "IFNL3"],
    },
    "immunosuppressants": {
        "description": "Immunosuppressant drugs with CYP3A5/TPMT-guided dosing",
        "drugs": [
            "tacrolimus", "cyclosporine", "mycophenolate", "sirolimus",
            "everolimus", "azathioprine",
        ],
        "primary_genes": ["CYP3A5", "CYP3A4", "TPMT", "NUDT15"],
    },
    "cardiovascular": {
        "description": "Cardiovascular drugs beyond anticoagulants with PGx relevance",
        "drugs": [
            "metoprolol", "carvedilol", "propranolol", "flecainide",
            "propafenone", "hydralazine", "isosorbide dinitrate", "ivabradine",
        ],
        "primary_genes": ["CYP2D6", "NAT2"],
    },
    "proton_pump_inhibitors": {
        "description": "PPIs with CYP2C19-dependent efficacy and dosing",
        "drugs": [
            "omeprazole", "esomeprazole", "lansoprazole", "dexlansoprazole",
            "pantoprazole", "rabeprazole",
        ],
        "primary_genes": ["CYP2C19"],
    },
    "anti_gout": {
        "description": "Gout therapies with HLA-based safety screening",
        "drugs": [
            "allopurinol", "febuxostat", "colchicine", "probenecid",
            "rasburicase",
        ],
        "primary_genes": ["HLA-B", "G6PD"],
    },
}


# =============================================================================
# 4. CYP_INHIBITORS — Drug-drug interaction knowledge
# =============================================================================

CYP_INHIBITORS: Dict[str, Dict[str, List[str]]] = {
    "CYP2D6": {
        "strong": [
            "fluoxetine", "paroxetine", "bupropion", "quinidine",
            "cinacalcet", "terbinafine",
        ],
        "moderate": [
            "duloxetine", "sertraline", "abiraterone", "mirabegron",
            "diphenhydramine",
        ],
        "weak": [
            "amiodarone", "cimetidine", "escitalopram", "haloperidol",
            "ranitidine", "venlafaxine", "celecoxib",
        ],
    },
    "CYP3A4": {
        "strong": [
            "ketoconazole", "itraconazole", "voriconazole", "posaconazole",
            "clarithromycin", "telithromycin", "ritonavir", "cobicistat",
            "nefazodone", "idelalisib", "nelfinavir",
        ],
        "moderate": [
            "erythromycin", "fluconazole", "diltiazem", "verapamil",
            "aprepitant", "crizotinib", "cyclosporine", "dronedarone",
            "grapefruit juice",
        ],
        "weak": [
            "alprazolam", "amiodarone", "amlodipine", "cilostazol",
            "cimetidine", "ranitidine", "ticagrelor",
        ],
    },
    "CYP2C19": {
        "strong": [
            "fluconazole", "fluoxetine", "fluvoxamine", "ticlopidine",
        ],
        "moderate": [
            "esomeprazole", "omeprazole", "voriconazole", "moclobemide",
            "felbamate",
        ],
        "weak": [
            "armodafinil", "cimetidine", "etravirine",
        ],
    },
    "CYP1A2": {
        "strong": [
            "fluvoxamine", "ciprofloxacin", "enoxacin",
        ],
        "moderate": [
            "mexiletine", "oral contraceptives", "vemurafenib", "methoxsalen",
        ],
        "weak": [
            "acyclovir", "allopurinol", "caffeine", "cimetidine",
            "echinacea", "famotidine", "zileuton",
        ],
    },
}


# =============================================================================
# 5. CYP_INDUCERS — Drug-drug interaction knowledge
# =============================================================================

CYP_INDUCERS: Dict[str, Dict[str, List[str]]] = {
    "CYP3A4": {
        "strong": [
            "rifampin", "carbamazepine", "phenytoin", "phenobarbital",
            "St. John's wort", "enzalutamide", "mitotane", "apalutamide",
        ],
        "moderate": [
            "efavirenz", "etravirine", "bosentan", "modafinil",
            "nafcillin", "dabrafenib",
        ],
    },
    "CYP1A2": {
        "strong": [
            "smoking (tobacco)", "charcoal-broiled foods", "rifampin",
            "phenytoin", "carbamazepine",
        ],
        "moderate": [
            "phenobarbital", "modafinil", "broccoli/cruciferous vegetables",
            "omeprazole (at high doses)", "insulin",
        ],
    },
    "CYP2C19": {
        "strong": [
            "rifampin", "St. John's wort",
        ],
        "moderate": [
            "carbamazepine", "phenytoin", "ritonavir", "efavirenz",
            "enzalutamide", "dabrafenib",
        ],
    },
}


# =============================================================================
# 6. HLA_DRUG_ASSOCIATIONS — 12 HLA-drug hypersensitivity pairs
# =============================================================================

HLA_DRUG_ASSOCIATIONS: Dict[str, Dict[str, Any]] = {
    "HLA-B*57:01_abacavir": {
        "hla_allele": "HLA-B*57:01",
        "drug": "abacavir",
        "reaction_type": "Abacavir Hypersensitivity Syndrome (AHS)",
        "risk_if_positive": "~50% risk of hypersensitivity reaction; can be fatal on rechallenge",
        "severity": "potentially fatal",
        "screening_mandatory": True,
        "recommendation": "Do NOT prescribe abacavir if HLA-B*57:01 positive. Use alternative NRTI (tenofovir).",
        "prevalence_by_population": {
            "European": "6-8%", "African": "1-2%", "East_Asian": "<1%",
            "South_Asian": "4-5%", "Hispanic": "2-3%",
        },
    },
    "HLA-B*15:02_carbamazepine": {
        "hla_allele": "HLA-B*15:02",
        "drug": "carbamazepine",
        "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)",
        "risk_if_positive": "~5-10% risk of SJS/TEN; mortality 10-30% for TEN",
        "severity": "life-threatening",
        "screening_mandatory": True,
        "recommendation": (
            "Do NOT prescribe carbamazepine if HLA-B*15:02 positive. "
            "Test all patients of Southeast Asian ancestry before initiating. "
            "Use alternative: levetiracetam, lamotrigine (with caution), valproic acid."
        ),
        "prevalence_by_population": {
            "Southeast_Asian": "10-15%", "East_Asian": "2-5%", "South_Asian": "2-4%",
            "European": "<1%", "African": "<1%",
        },
    },
    "HLA-B*15:02_oxcarbazepine": {
        "hla_allele": "HLA-B*15:02",
        "drug": "oxcarbazepine",
        "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)",
        "risk_if_positive": "Increased risk of SJS/TEN, lower than carbamazepine",
        "severity": "life-threatening",
        "screening_mandatory": True,
        "recommendation": (
            "Do NOT prescribe oxcarbazepine if HLA-B*15:02 positive. "
            "Same population screening recommendations as carbamazepine."
        ),
        "prevalence_by_population": {
            "Southeast_Asian": "10-15%", "East_Asian": "2-5%", "South_Asian": "2-4%",
            "European": "<1%", "African": "<1%",
        },
    },
    "HLA-B*15:02_phenytoin": {
        "hla_allele": "HLA-B*15:02",
        "drug": "phenytoin",
        "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN)",
        "risk_if_positive": "Increased risk of SJS/TEN; odds ratio ~4-5 in Southeast Asian populations",
        "severity": "life-threatening",
        "screening_mandatory": True,
        "recommendation": (
            "Consider avoiding phenytoin if HLA-B*15:02 positive. "
            "FDA black box warning for SE Asian descent. "
            "Alternatives: levetiracetam, valproic acid."
        ),
        "prevalence_by_population": {
            "Southeast_Asian": "10-15%", "East_Asian": "2-5%", "South_Asian": "2-4%",
            "European": "<1%", "African": "<1%",
        },
    },
    "HLA-A*31:01_carbamazepine": {
        "hla_allele": "HLA-A*31:01",
        "drug": "carbamazepine",
        "reaction_type": "Drug Reaction with Eosinophilia and Systemic Symptoms (DRESS) / SJS/maculopapular exanthema",
        "risk_if_positive": "~5-26% risk of hypersensitivity (DRESS > SJS); variable across populations",
        "severity": "serious to life-threatening",
        "screening_mandatory": True,
        "recommendation": (
            "Consider avoiding carbamazepine if HLA-A*31:01 positive. "
            "Primarily relevant in European and Japanese populations. "
            "CPIC Level A recommendation."
        ),
        "prevalence_by_population": {
            "European": "2-5%", "Japanese": "8-10%", "Korean": "7-8%",
            "African": "1-3%", "East_Asian": "5-8%",
        },
    },
    "HLA-B*58:01_allopurinol": {
        "hla_allele": "HLA-B*58:01",
        "drug": "allopurinol",
        "reaction_type": "Stevens-Johnson Syndrome / Toxic Epidermal Necrolysis (SJS/TEN) / DRESS",
        "risk_if_positive": "~2-5% risk of severe cutaneous adverse reaction; OR >80 in some studies",
        "severity": "life-threatening",
        "screening_mandatory": True,
        "recommendation": (
            "Do NOT prescribe allopurinol if HLA-B*58:01 positive. "
            "Use febuxostat or probenecid as alternative. "
            "CPIC, ACMG, ACR all recommend pre-treatment testing."
        ),
        "prevalence_by_population": {
            "Southeast_Asian": "6-8%", "East_Asian": "6-10%", "Korean": "12%",
            "African_American": "3-4%", "European": "1-2%",
        },
    },
    "HLA-B*57:01_flucloxacillin": {
        "hla_allele": "HLA-B*57:01",
        "drug": "flucloxacillin",
        "reaction_type": "Drug-Induced Liver Injury (DILI)",
        "risk_if_positive": "~1 in 500-1000 HLA-B*57:01 carriers; OR ~80",
        "severity": "serious — cholestatic hepatitis, can be prolonged",
        "screening_mandatory": False,
        "recommendation": (
            "Consider HLA-B*57:01 testing before flucloxacillin use. "
            "Low absolute risk limits universal screening cost-effectiveness. "
            "Monitor LFTs if prescribed."
        ),
        "prevalence_by_population": {
            "European": "6-8%", "African": "1-2%", "East_Asian": "<1%",
        },
    },
    "HLA-DRB1*07:01_lapatinib": {
        "hla_allele": "HLA-DRB1*07:01",
        "drug": "lapatinib",
        "reaction_type": "Drug-Induced Liver Injury (DILI) — hepatotoxicity",
        "risk_if_positive": "~2% ALT elevation risk in carriers (vs 0.3% non-carriers)",
        "severity": "moderate — hepatotoxicity, usually reversible",
        "screening_mandatory": False,
        "recommendation": (
            "HLA-DRB1*07:01 testing may be considered before lapatinib. "
            "Enhanced LFT monitoring in carriers. Not yet a CPIC guideline."
        ),
        "prevalence_by_population": {
            "European": "25-30%", "African": "10-15%", "East_Asian": "5-10%",
        },
    },
    "HLA-B*15:11_carbamazepine": {
        "hla_allele": "HLA-B*15:11",
        "drug": "carbamazepine",
        "reaction_type": "Stevens-Johnson Syndrome (SJS)",
        "risk_if_positive": "Increased SJS risk in Japanese and Korean populations",
        "severity": "life-threatening",
        "screening_mandatory": False,
        "recommendation": (
            "Consider testing for HLA-B*15:11 in Japanese/Korean patients before "
            "carbamazepine initiation, in addition to HLA-B*15:02 and HLA-A*31:01."
        ),
        "prevalence_by_population": {
            "Japanese": "1-2%", "Korean": "1-2%", "European": "<0.5%",
        },
    },
    "HLA-A*02:01_thalidomide": {
        "hla_allele": "HLA-A*02:01",
        "drug": "thalidomide",
        "reaction_type": "Improved antitumor response / erythema nodosum leprosum",
        "risk_if_positive": "Better response to thalidomide-based regimens in myeloma",
        "severity": "pharmacogenomic efficacy marker (not toxicity)",
        "screening_mandatory": False,
        "recommendation": (
            "HLA-A*02:01 may predict improved outcomes with thalidomide-based "
            "therapy. Research context; not yet actionable clinically."
        ),
        "prevalence_by_population": {
            "European": "25-50%", "African": "10-20%", "East_Asian": "5-10%",
        },
    },
    "HLA-B*35:05_nevirapine": {
        "hla_allele": "HLA-B*35:05",
        "drug": "nevirapine",
        "reaction_type": "Hepatotoxicity and skin rash (DRESS-like)",
        "risk_if_positive": "Increased risk of nevirapine-associated hypersensitivity",
        "severity": "serious — hepatotoxicity, rash; can progress to DRESS",
        "screening_mandatory": False,
        "recommendation": (
            "Consider HLA-B*35:05 testing in populations with high allele frequency. "
            "Alternative NNRTIs available (efavirenz, rilpivirine, doravirine). "
            "CD4 count >250 in women / >400 in men is an independent risk factor."
        ),
        "prevalence_by_population": {
            "European": "10-15%", "South_Asian": "8-12%", "African": "5-10%",
            "East_Asian": "5-8%",
        },
    },
    "HLA-DPB1*03:01_AERD": {
        "hla_allele": "HLA-DPB1*03:01",
        "drug": "aspirin",
        "reaction_type": "Aspirin-Exacerbated Respiratory Disease (AERD / Samter's triad)",
        "risk_if_positive": "Increased risk of AERD in susceptible individuals with asthma/nasal polyps",
        "severity": "moderate — bronchospasm, rhinosinusitis, nasal polyps",
        "screening_mandatory": False,
        "recommendation": (
            "HLA-DPB1*03:01 has been associated with AERD susceptibility. "
            "Aspirin desensitization may be considered. "
            "COX-2 selective inhibitors (celecoxib) as alternatives."
        ),
        "prevalence_by_population": {
            "European": "10-15%", "East_Asian": "5-10%", "African": "5-8%",
        },
    },
}


# =============================================================================
# 7. DRUG_ALTERNATIVES — Gene-phenotype-specific substitutions
# =============================================================================

DRUG_ALTERNATIVES: Dict[str, Dict[str, Any]] = {
    # --- Codeine / CYP2D6 ---
    "codeine_CYP2D6_PM": {
        "drug": "codeine",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 PM cannot convert codeine to morphine; analgesic failure expected.",
        "alternatives": [
            {"drug": "morphine", "rationale": "Direct mu-opioid agonist, no CYP2D6 activation required"},
            {"drug": "oxycodone", "rationale": "Partially CYP2D6-dependent but retains analgesic effect in PMs"},
            {"drug": "hydromorphone", "rationale": "No CYP2D6 dependence for activation"},
            {"drug": "non-opioid analgesics", "rationale": "NSAIDs, acetaminophen, gabapentinoids for mild-moderate pain"},
        ],
    },
    "codeine_CYP2D6_UM": {
        "drug": "codeine",
        "gene": "CYP2D6",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": (
            "CYP2D6 UM converts codeine to morphine too rapidly; risk of respiratory depression "
            "and death (FDA black box). Especially dangerous in children and breastfeeding mothers."
        ),
        "alternatives": [
            {"drug": "morphine (low dose, monitored)", "rationale": "Avoid prodrug conversion variability; titrate carefully"},
            {"drug": "hydromorphone", "rationale": "Not a CYP2D6 prodrug"},
            {"drug": "non-opioid analgesics", "rationale": "Preferred if pain severity allows"},
        ],
    },
    "tramadol_CYP2D6_PM": {
        "drug": "tramadol",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 converts tramadol to O-desmethyltramadol (active). PMs have reduced efficacy.",
        "alternatives": [
            {"drug": "morphine", "rationale": "No CYP2D6 activation required"},
            {"drug": "hydromorphone", "rationale": "Direct-acting opioid"},
            {"drug": "tapentadol", "rationale": "Dual mechanism; not CYP2D6 dependent"},
        ],
    },
    "tramadol_CYP2D6_UM": {
        "drug": "tramadol",
        "gene": "CYP2D6",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": "Rapid conversion to active metabolite; respiratory depression risk similar to codeine UM.",
        "alternatives": [
            {"drug": "non-opioid analgesics", "rationale": "Preferred if pain allows"},
            {"drug": "tapentadol", "rationale": "Not CYP2D6 dependent for activation"},
        ],
    },
    # --- Clopidogrel / CYP2C19 ---
    "clopidogrel_CYP2C19_PM": {
        "drug": "clopidogrel",
        "gene": "CYP2C19",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "Clopidogrel is a prodrug requiring CYP2C19 activation. PMs have dramatically reduced "
            "active metabolite, leading to insufficient platelet inhibition and increased MACE/stent thrombosis."
        ),
        "alternatives": [
            {"drug": "prasugrel", "rationale": "Alternative P2Y12 inhibitor; not CYP2C19 dependent. Contraindicated if prior stroke/TIA or >=75y or <60kg."},
            {"drug": "ticagrelor", "rationale": "Direct-acting P2Y12 inhibitor; no CYP2C19 activation needed. Twice-daily dosing, dyspnea side effect."},
        ],
    },
    "clopidogrel_CYP2C19_IM": {
        "drug": "clopidogrel",
        "gene": "CYP2C19",
        "phenotype": "Intermediate Metabolizer",
        "rationale": "Reduced CYP2C19 activity impairs clopidogrel activation; increased cardiovascular event risk.",
        "alternatives": [
            {"drug": "prasugrel", "rationale": "Not CYP2C19 dependent. Consider bleeding risk."},
            {"drug": "ticagrelor", "rationale": "Direct-acting; preferred in ACS if no contraindications."},
            {"drug": "clopidogrel (double dose)", "rationale": "Some evidence for 150mg maintenance in IMs, but less robust than switching."},
        ],
    },
    # --- Simvastatin / SLCO1B1 ---
    "simvastatin_SLCO1B1_decreased": {
        "drug": "simvastatin",
        "gene": "SLCO1B1",
        "phenotype": "Decreased Function",
        "rationale": "SLCO1B1 decreased function (*1/*5) increases simvastatin exposure 4.5-fold, elevating myopathy risk",
        "alternatives": [
            {"drug": "rosuvastatin", "rationale": "Less SLCO1B1-dependent transport"},
            {"drug": "pravastatin", "rationale": "Hydrophilic statin, less myopathy risk"},
            {"drug": "fluvastatin", "rationale": "Minimal SLCO1B1 transport"},
        ],
    },
    "simvastatin_SLCO1B1_poor": {
        "drug": "simvastatin",
        "gene": "SLCO1B1",
        "phenotype": "Poor Function (521 CC homozygote)",
        "rationale": (
            "SLCO1B1*5 homozygotes (521 CC) have ~17-fold increased risk of simvastatin myopathy. "
            "Impaired hepatic uptake leads to elevated systemic statin exposure."
        ),
        "alternatives": [
            {"drug": "rosuvastatin", "rationale": "Lower myopathy risk; less SLCO1B1-dependent"},
            {"drug": "pravastatin", "rationale": "Hydrophilic; lower myopathy risk"},
            {"drug": "fluvastatin", "rationale": "Minimal SLCO1B1 dependence"},
            {"drug": "simvastatin (low dose <=20mg)", "rationale": "If simvastatin preferred, limit to 20mg/day max. Avoid 80mg."},
        ],
    },
    "simvastatin_SLCO1B1_intermediate": {
        "drug": "simvastatin",
        "gene": "SLCO1B1",
        "phenotype": "Intermediate Function (521 TC heterozygote)",
        "rationale": "SLCO1B1*5 heterozygotes have ~3-fold increased myopathy risk with simvastatin.",
        "alternatives": [
            {"drug": "simvastatin (max 40mg)", "rationale": "Dose-limit; avoid 80mg"},
            {"drug": "rosuvastatin", "rationale": "Lower myopathy risk"},
            {"drug": "pravastatin", "rationale": "Hydrophilic alternative"},
        ],
    },
    # --- Warfarin / CYP2C9 ---
    "warfarin_CYP2C9_PM": {
        "drug": "warfarin",
        "gene": "CYP2C9",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2C9 PMs require significantly lower warfarin doses (often 50-80% reduction). "
            "Combined with VKORC1 -1639 AA, some patients need <1mg/day. "
            "Elevated bleeding risk at standard doses."
        ),
        "alternatives": [
            {"drug": "warfarin (reduced dose per algorithm)", "rationale": "Use FDA-table or IWPC algorithm incorporating CYP2C9, VKORC1, CYP4F2, clinical factors"},
            {"drug": "apixaban", "rationale": "DOAC; no CYP2C9/VKORC1 dependence; fixed dosing; lower ICH risk"},
            {"drug": "rivaroxaban", "rationale": "DOAC alternative; CYP3A4/P-gp substrate but no CYP2C9 issue"},
            {"drug": "edoxaban", "rationale": "DOAC alternative with fixed dosing"},
            {"drug": "dabigatran", "rationale": "Direct thrombin inhibitor; renal clearance, idarucizumab reversal available"},
        ],
    },
    "warfarin_VKORC1_sensitive": {
        "drug": "warfarin",
        "gene": "VKORC1",
        "phenotype": "Sensitive (-1639 AA genotype)",
        "rationale": "VKORC1 -1639 AA homozygotes require ~50% lower warfarin dose than GG carriers.",
        "alternatives": [
            {"drug": "warfarin (reduced dose per algorithm)", "rationale": "Incorporate VKORC1 genotype into dosing algorithm (Gage or IWPC)"},
            {"drug": "apixaban", "rationale": "DOAC; no VKORC1 dependence"},
            {"drug": "rivaroxaban", "rationale": "DOAC; fixed dosing"},
        ],
    },
    # --- 5-FU / DPYD ---
    "fluorouracil_DPYD_PM": {
        "drug": "fluorouracil (5-FU)",
        "gene": "DPYD",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "DPYD PMs (activity score 0-0.5) cannot adequately catabolize 5-FU. "
            "Fatal toxicity (severe mucositis, myelosuppression, neurotoxicity) reported in ~10-20% "
            "of untested DPYD-deficient patients receiving standard-dose 5-FU/capecitabine."
        ),
        "alternatives": [
            {"drug": "5-FU at 25-50% dose with monitoring (AS 0.5)", "rationale": "CPIC: 50% reduction for AS 1.0; avoid or 75% reduction for AS 0.5"},
            {"drug": "AVOID 5-FU/capecitabine entirely (AS 0)", "rationale": "Complete DPD deficiency — any dose can be fatal. Use non-fluoropyrimidine regimen."},
            {"drug": "raltitrexed (Tomudex)", "rationale": "Thymidylate synthase inhibitor; not DPD-dependent"},
            {"drug": "oxaliplatin-based regimen (without 5-FU)", "rationale": "When fluoropyrimidine cannot be used"},
        ],
    },
    "fluorouracil_DPYD_IM": {
        "drug": "fluorouracil (5-FU)",
        "gene": "DPYD",
        "phenotype": "Intermediate Metabolizer",
        "rationale": "DPYD IMs (activity score 1.0-1.5) have increased toxicity risk; dose reduction recommended.",
        "alternatives": [
            {"drug": "5-FU at 50% dose (AS 1.0)", "rationale": "CPIC recommends 50% reduction for AS 1.0, titrate up based on tolerance"},
            {"drug": "5-FU at 75% dose (AS 1.5)", "rationale": "CPIC recommends 25% reduction for AS 1.5, titrate based on tolerance"},
            {"drug": "capecitabine at same proportional reduction", "rationale": "Same DPYD pathway as 5-FU; same dose reduction principles apply"},
        ],
    },
    "capecitabine_DPYD_PM": {
        "drug": "capecitabine",
        "gene": "DPYD",
        "phenotype": "Poor Metabolizer",
        "rationale": "Capecitabine is a 5-FU prodrug; same DPD deficiency toxicity risk applies.",
        "alternatives": [
            {"drug": "AVOID capecitabine entirely (AS 0)", "rationale": "Complete DPD deficiency is a contraindication"},
            {"drug": "capecitabine at 25-50% dose (AS 0.5)", "rationale": "With intensive monitoring and growth factor support"},
            {"drug": "raltitrexed", "rationale": "Non-fluoropyrimidine TS inhibitor"},
        ],
    },
    # --- Carbamazepine / HLA ---
    "carbamazepine_HLA_A3101": {
        "drug": "carbamazepine",
        "gene": "HLA-A",
        "phenotype": "HLA-A*31:01 Positive",
        "rationale": "HLA-A*31:01 carriers have 5-26% risk of DRESS, SJS, or MPE",
        "alternatives": [
            {"drug": "levetiracetam", "rationale": "No HLA association"},
            {"drug": "valproic acid", "rationale": "No DRESS risk"},
            {"drug": "lacosamide", "rationale": "No HLA association"},
        ],
    },
    "carbamazepine_HLA_B1502": {
        "drug": "carbamazepine",
        "gene": "HLA-B",
        "phenotype": "HLA-B*15:02 Positive",
        "rationale": "HLA-B*15:02 carriers have 5-10% risk of fatal SJS/TEN with carbamazepine. FDA Boxed Warning",
        "alternatives": [
            {"drug": "levetiracetam", "rationale": "No HLA association, well tolerated"},
            {"drug": "valproic acid", "rationale": "Different mechanism, no SJS risk"},
            {"drug": "lacosamide", "rationale": "Newer agent, no HLA association"},
            {"drug": "lamotrigine", "rationale": "Lower SJS risk if HLA-B*15:02 negative; slow titration required"},
        ],
    },
    "carbamazepine_HLA_B1502_positive": {
        "drug": "carbamazepine",
        "gene": "HLA-B*15:02",
        "phenotype": "Positive",
        "rationale": "HLA-B*15:02 carriers have ~5-10% risk of SJS/TEN with carbamazepine. Contraindicated.",
        "alternatives": [
            {"drug": "levetiracetam", "rationale": "No HLA-mediated hypersensitivity risk; renal clearance"},
            {"drug": "valproic acid", "rationale": "Different mechanism; no HLA-B*15:02 association. Monitor LFTs."},
            {"drug": "lamotrigine", "rationale": "Lower SJS risk but still possible; slow titration required. Consider HLA testing."},
            {"drug": "lacosamide", "rationale": "Newer AED; no known HLA association"},
        ],
    },
    "carbamazepine_HLA_A3101_positive": {
        "drug": "carbamazepine",
        "gene": "HLA-A*31:01",
        "phenotype": "Positive",
        "rationale": "HLA-A*31:01 carriers have increased risk of DRESS/maculopapular exanthema with carbamazepine.",
        "alternatives": [
            {"drug": "levetiracetam", "rationale": "No HLA-mediated risk"},
            {"drug": "valproic acid", "rationale": "No HLA-A*31:01 association"},
            {"drug": "lamotrigine", "rationale": "Use with slow titration"},
            {"drug": "oxcarbazepine", "rationale": "Cross-reactivity possible but lower; use with caution and monitoring"},
        ],
    },
    # --- Abacavir / HLA ---
    "abacavir_HLA_B5701_positive": {
        "drug": "abacavir",
        "gene": "HLA-B*57:01",
        "phenotype": "Positive",
        "rationale": (
            "HLA-B*57:01 positive patients have ~50% risk of abacavir hypersensitivity syndrome. "
            "Rechallenge can be fatal. CPIC Level A — mandatory screening."
        ),
        "alternatives": [
            {"drug": "tenofovir disoproxil fumarate (TDF)", "rationale": "Standard NRTI backbone alternative; monitor renal function"},
            {"drug": "tenofovir alafenamide (TAF)", "rationale": "Improved renal/bone safety vs TDF; preferred in many guidelines"},
            {"drug": "zidovudine", "rationale": "Older NRTI; more side effects but no HLA-B*57:01 association"},
        ],
    },
    # --- Tamoxifen / CYP2D6 ---
    "tamoxifen_CYP2D6_PM": {
        "drug": "tamoxifen",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "Tamoxifen requires CYP2D6 conversion to endoxifen (active metabolite). "
            "PMs have 75% lower endoxifen levels and may have increased breast cancer recurrence risk."
        ),
        "alternatives": [
            {"drug": "aromatase inhibitor (anastrozole, letrozole, exemestane)", "rationale": "Not CYP2D6 dependent; preferred in postmenopausal ER+ breast cancer"},
            {"drug": "tamoxifen with CYP2D6 inhibitor avoidance", "rationale": "Ensure no concomitant strong CYP2D6 inhibitors (fluoxetine, paroxetine)"},
            {"drug": "toremifene", "rationale": "Similar mechanism to tamoxifen but less CYP2D6 dependent; limited data"},
        ],
    },
    "tamoxifen_CYP2D6_IM": {
        "drug": "tamoxifen",
        "gene": "CYP2D6",
        "phenotype": "Intermediate Metabolizer",
        "rationale": "Reduced endoxifen formation in CYP2D6 IMs; consider alternatives or dose increase.",
        "alternatives": [
            {"drug": "aromatase inhibitor", "rationale": "Not CYP2D6 dependent; if postmenopausal"},
            {"drug": "tamoxifen (40mg dose increase)", "rationale": "Some evidence 40mg compensates for reduced CYP2D6 activity in IMs"},
        ],
    },
    # --- Azathioprine / TPMT ---
    "azathioprine_TPMT_PM": {
        "drug": "azathioprine",
        "gene": "TPMT",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "TPMT PMs accumulate toxic thioguanine nucleotides (TGNs). Standard doses cause "
            "life-threatening myelosuppression (pancytopenia). Dose must be reduced 90% or drug avoided."
        ),
        "alternatives": [
            {"drug": "azathioprine at 10% standard dose (3x/week)", "rationale": "CPIC: start at drastically reduced dose; monitor CBC weekly for 8 weeks"},
            {"drug": "mycophenolate mofetil", "rationale": "Alternative immunosuppressant; not TPMT dependent"},
            {"drug": "methotrexate", "rationale": "Alternative for autoimmune conditions; different metabolism"},
        ],
    },
    "azathioprine_TPMT_IM": {
        "drug": "azathioprine",
        "gene": "TPMT",
        "phenotype": "Intermediate Metabolizer",
        "rationale": "TPMT IMs have moderately elevated TGN levels; dose reduction to 30-80% recommended.",
        "alternatives": [
            {"drug": "azathioprine at 50-80% dose", "rationale": "CPIC: reduce starting dose; titrate based on TGN levels and CBC"},
            {"drug": "mycophenolate mofetil", "rationale": "Non-TPMT alternative immunosuppressant"},
        ],
    },
    "azathioprine_NUDT15_PM": {
        "drug": "azathioprine",
        "gene": "NUDT15",
        "phenotype": "Poor Metabolizer",
        "rationale": "NUDT15 PMs (especially *2/*2, *3/*3) are at extreme myelosuppression risk with thiopurines.",
        "alternatives": [
            {"drug": "azathioprine at 10% dose", "rationale": "CPIC: drastic dose reduction; consider avoiding entirely"},
            {"drug": "mycophenolate mofetil", "rationale": "Non-thiopurine immunosuppressant"},
        ],
    },
    # --- Tacrolimus / CYP3A5 ---
    "tacrolimus_CYP3A5_expresser": {
        "drug": "tacrolimus",
        "gene": "CYP3A5",
        "phenotype": "Expresser (*1/*1 or *1/*3)",
        "rationale": (
            "CYP3A5 expressers metabolize tacrolimus faster, requiring higher starting doses "
            "(typically 1.5-2x standard) to achieve target trough levels."
        ),
        "alternatives": [
            {"drug": "tacrolimus 0.3 mg/kg/day", "rationale": "CPIC: increase starting dose 1.5-2x for CYP3A5 expressers"},
            {"drug": "cyclosporine", "rationale": "Alternative calcineurin inhibitor; not CYP3A5 dependent for dosing"},
        ],
    },
    "tacrolimus_CYP3A5_non_expresser": {
        "drug": "tacrolimus",
        "gene": "CYP3A5",
        "phenotype": "Non-expresser (*3/*3)",
        "rationale": "CYP3A5 non-expressers metabolize tacrolimus at standard rate. Use standard starting dose.",
        "alternatives": [
            {"drug": "tacrolimus 0.15 mg/kg/day", "rationale": "Standard starting dose per CPIC for non-expressers"},
        ],
    },
    # --- Allopurinol / HLA ---
    "allopurinol_HLA_B5801_positive": {
        "drug": "allopurinol",
        "gene": "HLA-B*58:01",
        "phenotype": "Positive",
        "rationale": "HLA-B*58:01 carriers have >80-fold increased risk of allopurinol-induced SJS/TEN/DRESS.",
        "alternatives": [
            {"drug": "febuxostat", "rationale": "Alternative xanthine oxidase inhibitor; no HLA-B*58:01 association. FDA CV boxed warning (CARES trial)."},
            {"drug": "probenecid", "rationale": "Uricosuric agent; different mechanism. Requires adequate renal function."},
            {"drug": "pegloticase", "rationale": "Recombinant uricase for refractory gout; IV infusion"},
        ],
    },
    # --- Irinotecan / UGT1A1 ---
    "irinotecan_UGT1A1_PM": {
        "drug": "irinotecan",
        "gene": "UGT1A1",
        "phenotype": "Poor Metabolizer (*28/*28 homozygote)",
        "rationale": (
            "UGT1A1*28 homozygotes have reduced glucuronidation of SN-38 (active metabolite). "
            "Increased risk of severe neutropenia and diarrhea at standard doses."
        ),
        "alternatives": [
            {"drug": "irinotecan at reduced dose (30% reduction)", "rationale": "FDA label: reduce starting dose by >=1 level for UGT1A1*28 homozygotes"},
            {"drug": "oxaliplatin-based regimen", "rationale": "FOLFOX instead of FOLFIRI if clinically appropriate"},
        ],
    },
    # --- Atomoxetine / CYP2D6 ---
    "atomoxetine_CYP2D6_PM": {
        "drug": "atomoxetine",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 PMs have ~10-fold higher atomoxetine AUC; increased side effects (insomnia, tachycardia).",
        "alternatives": [
            {"drug": "atomoxetine at 50% dose (0.5mg/kg)", "rationale": "CPIC/FDA: reduce starting dose; slower titration"},
            {"drug": "methylphenidate", "rationale": "Alternative ADHD medication; not CYP2D6 dependent"},
            {"drug": "lisdexamfetamine", "rationale": "Alternative ADHD medication; prodrug not CYP2D6 dependent"},
        ],
    },
    # --- Voriconazole / CYP2C19 ---
    "voriconazole_CYP2C19_PM": {
        "drug": "voriconazole",
        "gene": "CYP2C19",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2C19 PMs have ~4-fold higher voriconazole exposure; hepatotoxicity and neurotoxicity risk.",
        "alternatives": [
            {"drug": "voriconazole at 50% dose", "rationale": "With therapeutic drug monitoring (target trough 1-5 mcg/mL)"},
            {"drug": "posaconazole", "rationale": "Broad-spectrum azole; not CYP2C19 dependent"},
            {"drug": "isavuconazole", "rationale": "Alternative azole with less CYP2C19 dependence; fewer DDIs"},
        ],
    },
    "voriconazole_CYP2C19_UM": {
        "drug": "voriconazole",
        "gene": "CYP2C19",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": "CYP2C19 UMs have subtherapeutic voriconazole levels; treatment failure risk.",
        "alternatives": [
            {"drug": "voriconazole at higher dose with TDM", "rationale": "Increase dose and monitor troughs closely (target >1 mcg/mL)"},
            {"drug": "posaconazole", "rationale": "Not CYP2C19 dependent; good for mold coverage"},
            {"drug": "isavuconazole", "rationale": "Alternative broad-spectrum azole"},
        ],
    },
    # --- Ondansetron / CYP2D6 ---
    "ondansetron_CYP2D6_UM": {
        "drug": "ondansetron",
        "gene": "CYP2D6",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": "CYP2D6 UMs may have reduced ondansetron efficacy due to rapid clearance.",
        "alternatives": [
            {"drug": "granisetron", "rationale": "5-HT3 antagonist; less CYP2D6 dependent"},
            {"drug": "palonosetron", "rationale": "Longer-acting 5-HT3 antagonist; different metabolism"},
        ],
    },
    # --- Doxepin / CYP2D6 ---
    "doxepin_CYP2D6_PM": {
        "drug": "doxepin",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 PM causes elevated doxepin levels; increased sedation, anticholinergic effects, and QT risk",
        "alternatives": [
            {"drug": "mirtazapine", "rationale": "Sedating antidepressant, different pathway"},
            {"drug": "trazodone", "rationale": "Sedating, CYP3A4 substrate"},
            {"drug": "sertraline", "rationale": "SSRI alternative"},
        ],
    },
    # --- Efavirenz / CYP2B6 ---
    "efavirenz_CYP2B6_PM": {
        "drug": "efavirenz",
        "gene": "CYP2B6",
        "phenotype": "Poor Metabolizer (*6/*6)",
        "rationale": (
            "CYP2B6*6 homozygotes have 3-4x higher efavirenz levels. CNS side effects "
            "(vivid dreams, dizziness, depression) are more severe and persistent."
        ),
        "alternatives": [
            {"drug": "efavirenz at reduced dose (400mg or 200mg)", "rationale": "ENCORE1 trial supports 400mg; some CYP2B6 PMs may need 200mg"},
            {"drug": "dolutegravir", "rationale": "INSTI-based regimen; no CYP2B6 dependence; preferred 1st line"},
            {"drug": "rilpivirine", "rationale": "NNRTI alternative; not CYP2B6 dependent; requires food and VL<100k"},
        ],
    },
    # --- Escitalopram / CYP2C19 ---
    "escitalopram_CYP2C19_PM": {
        "drug": "escitalopram",
        "gene": "CYP2C19",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2C19 PM leads to increased escitalopram exposure; max 10mg/day due to QT risk",
        "alternatives": [
            {"drug": "sertraline", "rationale": "Alternative SSRI; dose adjust for CYP2C19 PM"},
            {"drug": "mirtazapine", "rationale": "NaSSA, minimal CYP2C19 dependence"},
            {"drug": "bupropion", "rationale": "NDRI, CYP2B6 substrate"},
        ],
    },
    # --- Fluvoxamine / CYP2D6 ---
    "fluvoxamine_CYP2D6_PM": {
        "drug": "fluvoxamine",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 PM may increase fluvoxamine side effects via minor metabolic pathway",
        "alternatives": [
            {"drug": "sertraline", "rationale": "Alternative SSRI with different CYP profile"},
            {"drug": "escitalopram", "rationale": "CYP2C19 primary; adjust if also CYP2C19 PM"},
        ],
    },
    # --- Phenytoin / CYP2C9 ---
    "phenytoin_CYP2C9_PM": {
        "drug": "phenytoin",
        "gene": "CYP2C9",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2C9 PMs have markedly reduced phenytoin clearance; toxicity at standard doses (nystagmus, ataxia, seizure).",
        "alternatives": [
            {"drug": "phenytoin at 50% dose", "rationale": "CPIC: reduce maintenance dose by >=50% for CYP2C9 PM; use TDM"},
            {"drug": "levetiracetam", "rationale": "Not CYP2C9 dependent; renal clearance"},
            {"drug": "valproic acid", "rationale": "Alternative AED; different metabolism"},
            {"drug": "lacosamide", "rationale": "Newer AED; minimal CYP dependence"},
        ],
    },
    # --- Phenytoin / HLA ---
    "phenytoin_HLA_B1502": {
        "drug": "phenytoin",
        "gene": "HLA-B",
        "phenotype": "HLA-B*15:02 Positive",
        "rationale": "HLA-B*15:02 carriers have elevated SJS/TEN risk with phenytoin, though lower than carbamazepine",
        "alternatives": [
            {"drug": "levetiracetam", "rationale": "No HLA-associated SJS risk"},
            {"drug": "valproic acid", "rationale": "No SJS risk"},
            {"drug": "lacosamide", "rationale": "No HLA association"},
        ],
    },
    # --- Mercaptopurine / TPMT+NUDT15 ---
    "mercaptopurine_TPMT_PM": {
        "drug": "mercaptopurine",
        "gene": "TPMT",
        "phenotype": "Poor Metabolizer",
        "rationale": "Same TPMT pathway as azathioprine. Drastic dose reduction or avoidance required.",
        "alternatives": [
            {"drug": "mercaptopurine at 10% dose (3x/week)", "rationale": "CPIC: reduce to 10% standard dose; monitor CBC weekly"},
            {"drug": "alternative ALL maintenance regimen", "rationale": "Consult oncology for non-thiopurine options"},
        ],
    },
    "thioguanine_TPMT_PM": {
        "drug": "thioguanine",
        "gene": "TPMT",
        "phenotype": "Poor Metabolizer",
        "rationale": "Thioguanine has complex metabolism; TPMT PMs at increased myelosuppression risk.",
        "alternatives": [
            {"drug": "thioguanine at reduced dose", "rationale": "CPIC: reduce dose and monitor; TGN levels if available"},
        ],
    },
    # --- Clozapine / CYP1A2 ---
    "clomipramine_CYP2D6_PM": {
        "drug": "clomipramine",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 PM causes elevated clomipramine and desmethylclomipramine levels; increased anticholinergic effects and seizure risk",
        "alternatives": [
            {"drug": "sertraline", "rationale": "SSRI, safer CYP profile"},
            {"drug": "fluoxetine", "rationale": "SSRI, note CYP2D6 inhibition"},
            {"drug": "fluvoxamine", "rationale": "SSRI, effective for OCD"},
        ],
    },
    "clozapine_CYP1A2_PM": {
        "drug": "clozapine",
        "gene": "CYP1A2",
        "phenotype": "Poor Metabolizer / Non-Inducer (non-smoker)",
        "rationale": (
            "CYP1A2 is the primary enzyme for clozapine metabolism. Non-smokers and "
            "CYP1A2 poor metabolizers have significantly higher clozapine levels, "
            "increasing risk of seizures, sedation, and metabolic syndrome."
        ),
        "alternatives": [
            {"drug": "clozapine at reduced dose with TDM", "rationale": "Target trough 350-600 ng/mL; reduce dose 30-50% in non-smokers vs smoker dose"},
            {"drug": "olanzapine", "rationale": "Alternative atypical antipsychotic; also CYP1A2 substrate but better tolerability window"},
            {"drug": "quetiapine", "rationale": "CYP3A4 substrate; no CYP1A2 dependence"},
        ],
    },
    "clozapine_CYP1A2_smoking_cessation": {
        "drug": "clozapine",
        "gene": "CYP1A2",
        "phenotype": "Smoking cessation (CYP1A2 de-induction)",
        "rationale": (
            "Smoking induces CYP1A2 activity 1.5-2x. Upon cessation, CYP1A2 returns to baseline "
            "within 3-7 days, causing clozapine levels to rise 50-70%. Risk of toxicity."
        ),
        "alternatives": [
            {"drug": "clozapine dose reduction 30-50% upon cessation", "rationale": "Reduce dose within 2-4 days of smoking cessation; monitor levels"},
        ],
    },
    # --- Aripiprazole / CYP2D6 ---
    "aripiprazole_CYP2D6_PM": {
        "drug": "aripiprazole",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2D6 PMs have ~80% higher aripiprazole exposure. Increased risk of akathisia, EPS.",
        "alternatives": [
            {"drug": "aripiprazole at 50-67% dose", "rationale": "FDA label: reduce to 50-67% of normal dose for CYP2D6 PMs"},
            {"drug": "quetiapine", "rationale": "CYP3A4 substrate; no CYP2D6 dependence"},
            {"drug": "ziprasidone", "rationale": "Aldehyde oxidase metabolism; minimal CYP2D6"},
        ],
    },
    # --- Metoprolol / CYP2D6 ---
    "metoprolol_CYP2D6_PM": {
        "drug": "metoprolol",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2D6 PMs have 5-fold higher metoprolol AUC. Risk of excessive "
            "bradycardia, hypotension, and heart failure worsening."
        ),
        "alternatives": [
            {"drug": "metoprolol at 25-50% dose", "rationale": "Titrate slowly; monitor HR and BP closely"},
            {"drug": "bisoprolol", "rationale": "50% renal / 50% hepatic clearance; less CYP2D6 dependent"},
            {"drug": "atenolol", "rationale": "Renally eliminated; no CYP2D6 dependence"},
            {"drug": "carvedilol", "rationale": "Multiple CYP pathways (2D6, 2C9, 1A2); less impact from single CYP2D6 PM"},
        ],
    },
    "metoprolol_CYP2D6_UM": {
        "drug": "metoprolol",
        "gene": "CYP2D6",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": "CYP2D6 UMs clear metoprolol rapidly; subtherapeutic levels and inadequate heart rate control.",
        "alternatives": [
            {"drug": "metoprolol at higher dose", "rationale": "May need 2-3x standard dose; guided by HR/BP response"},
            {"drug": "bisoprolol", "rationale": "Less CYP2D6 dependent; more predictable pharmacokinetics"},
            {"drug": "atenolol", "rationale": "No CYP2D6 dependence"},
        ],
    },
    # --- Rasburicase / G6PD ---
    "rasburicase_G6PD_deficient": {
        "drug": "rasburicase",
        "gene": "G6PD",
        "phenotype": "Deficient",
        "rationale": (
            "Rasburicase generates hydrogen peroxide as a byproduct. In G6PD-deficient patients, "
            "RBCs cannot neutralize H2O2 via glutathione pathway, causing hemolytic anemia. "
            "CONTRAINDICATED in G6PD deficiency — FDA black box warning."
        ),
        "alternatives": [
            {"drug": "allopurinol", "rationale": "Xanthine oxidase inhibitor; no oxidative hemolysis risk. Slower onset."},
            {"drug": "febuxostat", "rationale": "Alternative XOI for tumor lysis prophylaxis if allopurinol not tolerated"},
            {"drug": "aggressive IV hydration + monitoring", "rationale": "Non-pharmacologic TLS prevention strategy"},
        ],
    },
    # --- Isoniazid / NAT2 ---
    "isoniazid_NAT2_slow": {
        "drug": "isoniazid",
        "gene": "NAT2",
        "phenotype": "Slow Acetylator",
        "rationale": (
            "NAT2 slow acetylators have 5-6x higher isoniazid AUC. Increased risk of "
            "peripheral neuropathy (40-50% vs 5% in rapid acetylators) and hepatotoxicity. "
            "Pyridoxine (vitamin B6) supplementation essential."
        ),
        "alternatives": [
            {"drug": "isoniazid at reduced dose or frequency", "rationale": "Consider 5mg/kg 3x/week instead of daily; always add pyridoxine 25-50mg/day"},
            {"drug": "rifampin-based regimen (4R)", "rationale": "4 months rifampin alone for LTBI; avoids isoniazid entirely"},
            {"drug": "rifapentine + isoniazid (3HP)", "rationale": "Once-weekly for 12 weeks; may still have NAT2 effect but shorter exposure"},
        ],
    },
    # --- Flecainide / CYP2D6 ---
    "flecainide_CYP2D6_PM": {
        "drug": "flecainide",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2D6 PMs have ~2-fold higher flecainide exposure. Increased risk of "
            "proarrhythmia (QRS widening, ventricular tachycardia) at standard doses."
        ),
        "alternatives": [
            {"drug": "flecainide at 50% dose with ECG monitoring", "rationale": "Reduce dose; target trough 200-800 ng/mL; monitor QRS duration"},
            {"drug": "propafenone", "rationale": "Also CYP2D6 substrate but may have better therapeutic window"},
            {"drug": "amiodarone", "rationale": "Multi-channel blocker; not CYP2D6 dependent; more side effects"},
            {"drug": "sotalol", "rationale": "Renally cleared; no CYP2D6 dependence"},
        ],
    },
    # --- Nortriptyline / CYP2D6 ---
    "nortriptyline_CYP2D6_PM": {
        "drug": "nortriptyline",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2D6 PMs have markedly elevated nortriptyline levels. Increased risk of "
            "QTc prolongation, anticholinergic toxicity, and sedation."
        ),
        "alternatives": [
            {"drug": "nortriptyline at 50% dose", "rationale": "CPIC: 50% dose reduction; target plasma level 50-150 ng/mL"},
            {"drug": "sertraline", "rationale": "SSRI alternative; CYP2C19 pathway, less CYP2D6 impact"},
            {"drug": "escitalopram", "rationale": "SSRI alternative; primarily CYP2C19"},
        ],
    },
    "nortriptyline_CYP2D6_UM": {
        "drug": "nortriptyline",
        "gene": "CYP2D6",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": "CYP2D6 UMs have subtherapeutic nortriptyline levels; treatment failure likely.",
        "alternatives": [
            {"drug": "nortriptyline at higher dose with TDM", "rationale": "Increase dose; monitor levels; may not achieve therapeutic range"},
            {"drug": "sertraline", "rationale": "Not CYP2D6 dependent; alternative antidepressant"},
            {"drug": "venlafaxine", "rationale": "SNRI with CYP2D6 metabolism but retains efficacy as both parent and metabolite are active"},
        ],
    },
    # --- Amitriptyline / CYP2D6+CYP2C19 ---
    "amitriptyline_CYP2C19_PM": {
        "drug": "amitriptyline",
        "gene": "CYP2C19",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2C19 PMs have reduced demethylation of amitriptyline to nortriptyline. "
            "Higher parent drug (amitriptyline) levels with more sedation and anticholinergic effects."
        ),
        "alternatives": [
            {"drug": "amitriptyline at 50% dose", "rationale": "CPIC: 50% dose reduction for CYP2C19 PM"},
            {"drug": "nortriptyline", "rationale": "Active metabolite; bypasses CYP2C19 step (but CYP2D6 dependent)"},
            {"drug": "SSRI (sertraline, escitalopram)", "rationale": "Alternative antidepressant class"},
        ],
    },
    "amitriptyline_CYP2D6_PM": {
        "drug": "amitriptyline",
        "gene": "CYP2D6",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2D6 PMs have reduced hydroxylation of both amitriptyline and nortriptyline. "
            "Elevated total TCA levels; increased toxicity risk."
        ),
        "alternatives": [
            {"drug": "amitriptyline at 50% dose", "rationale": "CPIC: 50% dose reduction for CYP2D6 PM; monitor TCA levels"},
            {"drug": "sertraline", "rationale": "SSRI; not CYP2D6 dependent for efficacy"},
            {"drug": "escitalopram", "rationale": "SSRI; primarily CYP2C19 metabolism"},
        ],
    },
    # --- Citalopram / CYP2C19 ---
    "citalopram_CYP2C19_PM": {
        "drug": "citalopram",
        "gene": "CYP2C19",
        "phenotype": "Poor Metabolizer",
        "rationale": (
            "CYP2C19 PMs have ~2-fold higher citalopram exposure. FDA max dose is 20mg/day "
            "for CYP2C19 PMs due to QTc prolongation risk."
        ),
        "alternatives": [
            {"drug": "citalopram max 20mg/day", "rationale": "FDA label: do not exceed 20mg/day in CYP2C19 PMs"},
            {"drug": "sertraline", "rationale": "Less CYP2C19 dependent; alternative SSRI"},
            {"drug": "fluoxetine", "rationale": "CYP2D6 substrate; alternative SSRI if CYP2D6 normal"},
        ],
    },
    "citalopram_CYP2C19_UM": {
        "drug": "citalopram",
        "gene": "CYP2C19",
        "phenotype": "Ultra-rapid Metabolizer",
        "rationale": "CYP2C19 UMs rapidly clear citalopram; may need higher dose or alternative.",
        "alternatives": [
            {"drug": "citalopram at higher dose", "rationale": "May need doses above standard starting dose; monitor response"},
            {"drug": "sertraline", "rationale": "Less CYP2C19 dependent SSRI"},
            {"drug": "fluoxetine", "rationale": "CYP2D6 substrate; not CYP2C19 dependent"},
        ],
    },
    # --- Sertraline / CYP2C19 ---
    "sertraline_CYP2C19_PM": {
        "drug": "sertraline",
        "gene": "CYP2C19",
        "phenotype": "Poor Metabolizer",
        "rationale": "CYP2C19 PM may increase sertraline exposure; consider dose reduction",
        "alternatives": [
            {"drug": "bupropion", "rationale": "CYP2B6 substrate, no CYP2C19 dependence"},
            {"drug": "mirtazapine", "rationale": "NaSSA, minimal CYP2C19 involvement"},
        ],
    },
    # --- Succinylcholine / RYR1/CACNA1S ---
    "succinylcholine_RYR1_MH_susceptible": {
        "drug": "succinylcholine",
        "gene": "RYR1",
        "phenotype": "Malignant Hyperthermia Susceptible",
        "rationale": (
            "RYR1 pathogenic variant carriers are at high risk of malignant hyperthermia "
            "when exposed to succinylcholine or volatile anesthetics. MH mortality ~5-10% "
            "even with dantrolene rescue."
        ),
        "alternatives": [
            {"drug": "rocuronium", "rationale": "Non-depolarizing neuromuscular blocker; no MH trigger. Sugammadex reversal available."},
            {"drug": "vecuronium", "rationale": "Non-depolarizing agent; safe in MH-susceptible patients"},
            {"drug": "cisatracurium", "rationale": "Non-depolarizing; Hofmann elimination; no MH risk"},
        ],
    },
    # --- Volatile anesthetics / RYR1 ---
    "volatile_anesthetics_RYR1_MH_susceptible": {
        "drug": "volatile anesthetics (desflurane, sevoflurane, isoflurane)",
        "gene": "RYR1",
        "phenotype": "Malignant Hyperthermia Susceptible",
        "rationale": (
            "All volatile anesthetics trigger MH in susceptible individuals via uncontrolled "
            "RYR1 calcium release. Presents as hyperthermia, muscle rigidity, rhabdomyolysis, "
            "hyperkalemia, metabolic acidosis."
        ),
        "alternatives": [
            {"drug": "total intravenous anesthesia (TIVA)", "rationale": "Propofol + remifentanil; no MH trigger"},
            {"drug": "nitrous oxide (safe)", "rationale": "Not an MH trigger; can supplement TIVA"},
            {"drug": "regional anesthesia", "rationale": "Neuraxial or peripheral nerve blocks when surgically appropriate"},
        ],
    },
    # --- Primaquine / G6PD ---
    "primaquine_G6PD_deficient": {
        "drug": "primaquine",
        "gene": "G6PD",
        "phenotype": "Deficient",
        "rationale": (
            "Primaquine causes oxidative stress leading to hemolytic anemia in G6PD-deficient "
            "patients. Severity depends on G6PD variant class (Class I-III)."
        ),
        "alternatives": [
            {"drug": "tafenoquine (with G6PD testing)", "rationale": "Also requires G6PD testing; longer half-life means prolonged hemolysis risk if deficient"},
            {"drug": "chloroquine (for P. vivax suppression)", "rationale": "Does not cause hemolysis but does not achieve radical cure (hypnozoite eradication)"},
            {"drug": "primaquine 0.75mg/kg weekly x 8 weeks", "rationale": "WHO: weekly dosing may be tolerated in mild G6PD deficiency (Class III/A-)"},
        ],
    },
    # --- Hormonal contraceptives / F5 ---
    "hormonal_contraceptives_F5_Leiden": {
        "drug": "combined hormonal contraceptives",
        "gene": "F5",
        "phenotype": "Factor V Leiden carrier",
        "rationale": (
            "Factor V Leiden heterozygotes have 3-8x baseline VTE risk. Combined oral "
            "contraceptives add 3-4x multiplicative risk, yielding ~15-30x VTE risk. "
            "Homozygotes have even higher risk."
        ),
        "alternatives": [
            {"drug": "progestin-only pill", "rationale": "Minimal VTE risk increase; effective contraception"},
            {"drug": "levonorgestrel IUD", "rationale": "Local progestin; negligible systemic VTE effect"},
            {"drug": "copper IUD", "rationale": "Non-hormonal; no VTE risk"},
            {"drug": "etonogestrel implant", "rationale": "Progestin-only; low systemic VTE impact"},
        ],
    },
}


# =============================================================================
# 8. ACTIVITY_SCORE_TABLES — CPIC-based activity score mappings
# =============================================================================

ACTIVITY_SCORE_TABLES: Dict[str, Dict[str, Any]] = {
    "CYP2D6": {
        "description": "CYP2D6 allele activity values per CPIC 2019 guidelines",
        "allele_scores": {
            "*1": 1.0, "*2": 1.0, "*35": 1.0, "*9": 0.5, "*10": 0.25,
            "*17": 0.5, "*29": 0.5, "*41": 0.5,
            "*3": 0, "*4": 0, "*5": 0, "*6": 0, "*7": 0, "*8": 0,
            "*11": 0, "*12": 0, "*14": 0, "*15": 0, "*36": 0,
        },
        "phenotype_thresholds": {
            "Ultra-rapid Metabolizer": ">2.25 (includes gene duplication of functional allele)",
            "Rapid Metabolizer": "Not used for CYP2D6 — skip from NM to UM",
            "Normal Metabolizer": "1.25-2.25",
            "Intermediate Metabolizer": "0.25-1.0",
            "Poor Metabolizer": "0",
        },
        "duplication_note": (
            "Gene duplications (*1xN, *2xN) add the allele score for each copy. "
            "E.g., *1/*2x2 = 1.0 + 2.0 = 3.0 (UM). *41x2/*41 = 1.0 + 0.5 = 1.5 (NM)."
        ),
    },
    "DPYD": {
        "description": "DPYD activity score per CPIC 2018 fluoropyrimidine guideline",
        "allele_scores": {
            "reference (normal)": 1.0,
            "c.1905+1G>A (*2A, rs3918290)": 0,
            "c.1679T>G (*13, rs55886062)": 0,
            "c.2846A>T (rs67376798)": 0.5,
            "c.1129-5923C>G (HapB3, rs75017182)": 0.5,
            "c.1601G>A (*4, rs1801158)": 0.5,
            "c.1236G>A (*5, rs56038477)": 0.5,
        },
        "phenotype_thresholds": {
            "Normal Metabolizer": "2.0 (two normal alleles)",
            "Intermediate Metabolizer (reduced risk)": "1.5 (one normal + one decreased)",
            "Intermediate Metabolizer (increased risk)": "1.0 (two decreased or one normal + one no function)",
            "Poor Metabolizer": "0-0.5 (two no-function, or one no-function + one decreased)",
        },
        "dosing_recommendations": {
            "AS 2.0": "100% standard dose",
            "AS 1.5": "75% dose (25% reduction); titrate based on toxicity/response",
            "AS 1.0": "50% dose (50% reduction); titrate based on toxicity/response",
            "AS 0.5": "25-50% dose with intensive monitoring; consider alternative",
            "AS 0": "AVOID fluoropyrimidines entirely; use alternative regimen",
        },
    },
}


# =============================================================================
# 8b. POPULATION_FREQUENCIES — Allele frequencies across global populations
# =============================================================================

POPULATION_FREQUENCIES: Dict[str, Dict[str, Dict[str, str]]] = {
    "CYP2D6": {
        "*4 (no function)": {
            "European": "20-25%", "African": "2-6%", "East_Asian": "1%",
            "South_Asian": "7-10%", "Middle_Eastern": "5-8%",
        },
        "*5 (gene deletion)": {
            "European": "2-7%", "African": "4-6%", "East_Asian": "6%",
            "South_Asian": "3-5%",
        },
        "*10 (decreased function)": {
            "European": "1-2%", "African": "3-5%", "East_Asian": "35-50%",
            "South_Asian": "15-25%",
        },
        "*17 (decreased function)": {
            "European": "<1%", "African": "20-35%", "East_Asian": "<1%",
        },
        "*41 (decreased function)": {
            "European": "7-10%", "African": "5-8%", "Middle_Eastern": "15-20%",
            "East_Asian": "2-5%",
        },
        "gene duplication (*1xN, *2xN)": {
            "European": "1-5%", "African": "2-10%", "East_Asian": "0-2%",
            "Middle_Eastern": "10-15%", "North_African": "10-29%",
            "Oceanian": "5-10%",
        },
        "PM phenotype (overall)": {
            "European": "5-10%", "African": "1-3%", "East_Asian": "<1%",
            "Middle_Eastern": "1-3%",
        },
        "UM phenotype (overall)": {
            "European": "1-5%", "African": "2-10%", "Middle_Eastern": "10-16%",
            "North_African": "10-29%", "East_Asian": "<1%",
        },
    },
    "CYP2C19": {
        "*2 (no function)": {
            "European": "12-15%", "African": "15-18%", "East_Asian": "25-35%",
            "South_Asian": "30-40%", "Oceanian": "40-70%",
        },
        "*3 (no function)": {
            "European": "<1%", "African": "<1%", "East_Asian": "5-10%",
            "South_Asian": "2-5%",
        },
        "*17 (increased function)": {
            "European": "20-25%", "African": "15-20%", "East_Asian": "<5%",
            "Middle_Eastern": "15-20%",
        },
        "PM phenotype (overall)": {
            "European": "2-5%", "African": "3-5%", "East_Asian": "12-23%",
            "South_Asian": "10-15%", "Oceanian": "40-70%",
        },
        "UM phenotype (overall)": {
            "European": "5-10%", "African": "4-8%", "East_Asian": "<1%",
            "Middle_Eastern": "5-10%",
        },
    },
    "CYP2C9": {
        "*2 (decreased function)": {
            "European": "10-15%", "African": "1-3%", "East_Asian": "<1%",
            "South_Asian": "4-8%",
        },
        "*3 (decreased function)": {
            "European": "5-8%", "African": "<1%", "East_Asian": "2-5%",
            "South_Asian": "5-10%",
        },
        "PM phenotype (overall)": {
            "European": "1-3%", "African": "<0.5%", "East_Asian": "<0.5%",
        },
    },
    "DPYD": {
        "*2A (no function)": {
            "European": "1-2%", "African": "0.1%", "East_Asian": "<0.1%",
        },
        "c.2846A>T (decreased function)": {
            "European": "1-1.5%", "African": "<0.5%", "East_Asian": "<0.5%",
        },
        "HapB3 (decreased function)": {
            "European": "2-5%", "African": "<1%", "East_Asian": "<1%",
        },
        "DPD deficiency (any)": {
            "European": "3-8%", "African": "1-3%", "East_Asian": "1-2%",
        },
    },
    "TPMT": {
        "*3A (no function)": {
            "European": "3-5%", "African": "<1%", "East_Asian": "<1%",
        },
        "*3C (no function)": {
            "European": "<1%", "African": "5-8%", "East_Asian": "2-5%",
        },
        "PM phenotype (overall)": {
            "European": "0.3%", "African": "0.2%", "East_Asian": "0.1%",
        },
        "IM phenotype (overall)": {
            "European": "10%", "African": "5-7%", "East_Asian": "3-5%",
        },
    },
    "NUDT15": {
        "*2 (no function, R139C)": {
            "European": "<1%", "African": "<0.5%", "East_Asian": "5-10%",
            "South_Asian": "3-7%", "Hispanic": "2-5%",
        },
        "*3 (no function, R139H)": {
            "European": "<0.5%", "East_Asian": "2-4%",
        },
    },
    "SLCO1B1": {
        "*5 (521T>C, decreased function)": {
            "European": "15-20%", "African": "2-5%", "East_Asian": "10-15%",
            "South_Asian": "8-12%",
        },
    },
    "CYP3A5": {
        "*3 (non-expresser)": {
            "European": "90-95%", "African": "30-50%", "East_Asian": "60-75%",
            "South_Asian": "60-70%", "Hispanic": "70-80%",
        },
        "*1 (expresser)": {
            "European": "5-10%", "African": "50-70%", "East_Asian": "25-40%",
            "South_Asian": "30-40%",
        },
    },
    "UGT1A1": {
        "*28 (7 TA repeats)": {
            "European": "25-35%", "African": "40-55%", "East_Asian": "10-15%",
        },
        "*6 (G71R)": {
            "European": "<1%", "East_Asian": "15-25%", "South_Asian": "5-10%",
        },
    },
    "NAT2": {
        "slow acetylator phenotype": {
            "European": "50-60%", "African": "40-50%", "East_Asian": "10-25%",
            "South_Asian": "30-40%", "Middle_Eastern": "55-65%",
        },
    },
    "G6PD": {
        "deficiency (any variant)": {
            "African_Male": "10-15%", "Mediterranean_Male": "5-10%",
            "Southeast_Asian_Male": "5-15%", "Middle_Eastern_Male": "5-10%",
            "European_Male": "1-2%", "East_Asian_Male": "2-5%",
        },
    },
    "VKORC1": {
        "-1639 A allele (warfarin sensitive)": {
            "European": "35-40%", "African": "10-15%", "East_Asian": "85-95%",
            "South_Asian": "45-55%",
        },
    },
}


# =============================================================================
# 8c. CLINICAL_DECISION_NOTES — Key clinical pearls and caveats
# =============================================================================

CLINICAL_DECISION_NOTES: Dict[str, str] = {
    "phenoconversion": (
        "Phenoconversion occurs when a concomitant CYP inhibitor converts a genotypic "
        "Normal Metabolizer into a phenotypic Poor Metabolizer. Example: a CYP2D6 NM "
        "taking paroxetine (strong CYP2D6 inhibitor) effectively becomes a CYP2D6 PM "
        "for codeine metabolism. Always assess both genotype AND concomitant medications."
    ),
    "gene_gene_interactions": (
        "Patients may carry actionable variants in multiple pharmacogenes simultaneously. "
        "Example: CYP2C9 PM + VKORC1 -1639 AA requires much lower warfarin doses than "
        "either alone. CYP2D6 PM + CYP2C19 PM narrows antidepressant choices significantly. "
        "Always evaluate the complete pharmacogenomic panel, not single genes."
    ),
    "prodrug_vs_active_drug": (
        "The clinical impact of a PM or UM phenotype reverses depending on whether "
        "the drug is a prodrug or active drug. PM + prodrug = reduced efficacy (cannot "
        "activate). PM + active drug = toxicity risk (cannot clear). UM + prodrug = "
        "toxicity (excessive activation). UM + active drug = reduced efficacy (rapid clearance)."
    ),
    "pediatric_considerations": (
        "Pharmacogenomic effects may be amplified in pediatric patients due to developmental "
        "differences in CYP expression. CYP2D6 activity increases from birth to age 10. "
        "CYP3A4 is relatively high in infancy. Codeine is contraindicated in children <12 "
        "and in post-tonsillectomy patients <18 regardless of CYP2D6 status (FDA)."
    ),
    "pregnancy_considerations": (
        "CYP2D6 activity increases during pregnancy (especially 2nd/3rd trimester) due to "
        "hormonal induction. A genotypic IM may phenotypically normalize. CYP1A2 activity "
        "decreases. Dose adjustments may be needed for affected substrates during pregnancy "
        "and again postpartum when enzyme activity returns to baseline."
    ),
    "ethnicity_and_allele_frequency": (
        "Pharmacogene allele frequencies vary dramatically by population. CYP2D6*10 is "
        "common in East Asians (35-50%) but rare in Europeans (1-2%). CYP2D6*17 is "
        "common in Africans (20-35%) but rare elsewhere. CYP2C19*2 frequency ranges "
        "from 12% (European) to 70% (Oceanian/Melanesian). Always consider patient "
        "ancestry when interpreting results and selecting pre-emptive panels."
    ),
    "preemptive_vs_reactive_testing": (
        "Preemptive testing (panel-based, before prescribing) is more cost-effective than "
        "reactive testing (single gene, after adverse event). The PREPARE trial (Swen et al., "
        "Lancet 2023) demonstrated 30% reduction in ADRs with preemptive 12-gene PGx panel. "
        "Major health systems (St. Jude, Vanderbilt, Mayo) have implemented preemptive programs."
    ),
    "star_allele_limitations": (
        "Star allele nomenclature has limitations: (1) Not all variants are captured by "
        "standard genotyping panels; rare alleles may be missed, defaulting to *1 (reference). "
        "(2) Structural variants (CYP2D6 deletions, duplications, hybrids) require copy number "
        "analysis. (3) Novel variants may not have assigned star alleles. (4) Phasing "
        "(determining which alleles are on which chromosome) affects diplotype assignment."
    ),
    "therapeutic_drug_monitoring": (
        "TDM complements PGx testing for narrow therapeutic index drugs. Use TDM to: "
        "(1) confirm predicted phenotype effect, (2) adjust doses in phenoconversion, "
        "(3) account for non-genetic factors (organ function, age, weight), "
        "(4) guide dose titration after initial PGx-guided starting dose. "
        "Key TDM drugs: tacrolimus, voriconazole, clozapine, tricyclics, aminoglycosides."
    ),
    "cpic_vs_dpwg_differences": (
        "CPIC (Clinical Pharmacogenetics Implementation Consortium) and DPWG (Dutch "
        "Pharmacogenetics Working Group) sometimes differ in recommendations. Key differences: "
        "(1) DPWG may be more conservative for some CYP2D6 IM scenarios, "
        "(2) CPIC provides activity scores while DPWG uses categorical phenotypes, "
        "(3) Both are evidence-based but use different grading systems. "
        "When in doubt, follow the more conservative recommendation."
    ),
}


# =============================================================================
# 9. ENTITY_ALIASES — 80+ aliases for entity resolution
# =============================================================================

ENTITY_ALIASES: Dict[str, Dict[str, str]] = {
    # --- Drug brand → generic ---
    "COUMADIN": {"type": "drug", "canonical": "warfarin"},
    "JANTOVEN": {"type": "drug", "canonical": "warfarin"},
    "PLAVIX": {"type": "drug", "canonical": "clopidogrel"},
    "LIPITOR": {"type": "drug", "canonical": "atorvastatin"},
    "ZOCOR": {"type": "drug", "canonical": "simvastatin"},
    "CRESTOR": {"type": "drug", "canonical": "rosuvastatin"},
    "PRAVACHOL": {"type": "drug", "canonical": "pravastatin"},
    "LESCOL": {"type": "drug", "canonical": "fluvastatin"},
    "LIVALO": {"type": "drug", "canonical": "pitavastatin"},
    "MEVACOR": {"type": "drug", "canonical": "lovastatin"},
    "PROZAC": {"type": "drug", "canonical": "fluoxetine"},
    "PAXIL": {"type": "drug", "canonical": "paroxetine"},
    "ZOLOFT": {"type": "drug", "canonical": "sertraline"},
    "CELEXA": {"type": "drug", "canonical": "citalopram"},
    "LEXAPRO": {"type": "drug", "canonical": "escitalopram"},
    "EFFEXOR": {"type": "drug", "canonical": "venlafaxine"},
    "CYMBALTA": {"type": "drug", "canonical": "duloxetine"},
    "ELAVIL": {"type": "drug", "canonical": "amitriptyline"},
    "PAMELOR": {"type": "drug", "canonical": "nortriptyline"},
    "ANAFRANIL": {"type": "drug", "canonical": "clomipramine"},
    "TEGRETOL": {"type": "drug", "canonical": "carbamazepine"},
    "TRILEPTAL": {"type": "drug", "canonical": "oxcarbazepine"},
    "DILANTIN": {"type": "drug", "canonical": "phenytoin"},
    "LAMICTAL": {"type": "drug", "canonical": "lamotrigine"},
    "DEPAKOTE": {"type": "drug", "canonical": "valproic acid"},
    "KEPPRA": {"type": "drug", "canonical": "levetiracetam"},
    "ZIAGEN": {"type": "drug", "canonical": "abacavir"},
    "SUSTIVA": {"type": "drug", "canonical": "efavirenz"},
    "VIRAMUNE": {"type": "drug", "canonical": "nevirapine"},
    "REYATAZ": {"type": "drug", "canonical": "atazanavir"},
    "PROGRAF": {"type": "drug", "canonical": "tacrolimus"},
    "NEORAL": {"type": "drug", "canonical": "cyclosporine"},
    "SANDIMMUNE": {"type": "drug", "canonical": "cyclosporine"},
    "IMURAN": {"type": "drug", "canonical": "azathioprine"},
    "PURINETHOL": {"type": "drug", "canonical": "mercaptopurine"},
    "TABLOID": {"type": "drug", "canonical": "thioguanine"},
    "NOLVADEX": {"type": "drug", "canonical": "tamoxifen"},
    "XELODA": {"type": "drug", "canonical": "capecitabine"},
    "ADRUCIL": {"type": "drug", "canonical": "fluorouracil"},
    "CAMPTOSAR": {"type": "drug", "canonical": "irinotecan"},
    "ZYLOPRIM": {"type": "drug", "canonical": "allopurinol"},
    "ULORIC": {"type": "drug", "canonical": "febuxostat"},
    "PRILOSEC": {"type": "drug", "canonical": "omeprazole"},
    "NEXIUM": {"type": "drug", "canonical": "esomeprazole"},
    "PREVACID": {"type": "drug", "canonical": "lansoprazole"},
    "PROTONIX": {"type": "drug", "canonical": "pantoprazole"},
    "ACIPHEX": {"type": "drug", "canonical": "rabeprazole"},
    "ELIQUIS": {"type": "drug", "canonical": "apixaban"},
    "XARELTO": {"type": "drug", "canonical": "rivaroxaban"},
    "PRADAXA": {"type": "drug", "canonical": "dabigatran"},
    "SAVAYSA": {"type": "drug", "canonical": "edoxaban"},
    "STRATTERA": {"type": "drug", "canonical": "atomoxetine"},
    "VFEND": {"type": "drug", "canonical": "voriconazole"},
    "WELLBUTRIN": {"type": "drug", "canonical": "bupropion"},
    "ABILIFY": {"type": "drug", "canonical": "aripiprazole"},
    "CLOZARIL": {"type": "drug", "canonical": "clozapine"},
    "ZYPREXA": {"type": "drug", "canonical": "olanzapine"},
    "RISPERDAL": {"type": "drug", "canonical": "risperidone"},
    "SEROQUEL": {"type": "drug", "canonical": "quetiapine"},
    "HALDOL": {"type": "drug", "canonical": "haloperidol"},
    "LOPRESSOR": {"type": "drug", "canonical": "metoprolol"},
    "TOPROL": {"type": "drug", "canonical": "metoprolol"},
    "COREG": {"type": "drug", "canonical": "carvedilol"},
    "INDERAL": {"type": "drug", "canonical": "propranolol"},
    "TRIKAFTA": {"type": "drug", "canonical": "elexacaftor/tezacaftor/ivacaftor"},
    "KALYDECO": {"type": "drug", "canonical": "ivacaftor"},
    "ORKAMBI": {"type": "drug", "canonical": "lumacaftor/ivacaftor"},
    "SYMDEKO": {"type": "drug", "canonical": "tezacaftor/ivacaftor"},
    # --- Gene aliases ---
    "CYTOCHROME P450 2D6": {"type": "gene", "canonical": "CYP2D6"},
    "CYTOCHROME P450 2C19": {"type": "gene", "canonical": "CYP2C19"},
    "CYTOCHROME P450 2C9": {"type": "gene", "canonical": "CYP2C9"},
    "CYTOCHROME P450 3A4": {"type": "gene", "canonical": "CYP3A4"},
    "CYTOCHROME P450 3A5": {"type": "gene", "canonical": "CYP3A5"},
    "CYTOCHROME P450 1A2": {"type": "gene", "canonical": "CYP1A2"},
    "CYTOCHROME P450 2B6": {"type": "gene", "canonical": "CYP2B6"},
    "P-GLYCOPROTEIN": {"type": "gene", "canonical": "ABCB1"},
    "P-GP": {"type": "gene", "canonical": "ABCB1"},
    "PGP": {"type": "gene", "canonical": "ABCB1"},
    "MDR1": {"type": "gene", "canonical": "ABCB1"},
    "MULTIDRUG RESISTANCE 1": {"type": "gene", "canonical": "ABCB1"},
    "DPD": {"type": "gene", "canonical": "DPYD"},
    "DIHYDROPYRIMIDINE DEHYDROGENASE": {"type": "gene", "canonical": "DPYD"},
    "OATP1B1": {"type": "gene", "canonical": "SLCO1B1"},
    "IL28B": {"type": "gene", "canonical": "IFNL3"},
    "INTERFERON LAMBDA 3": {"type": "gene", "canonical": "IFNL3"},
    "FACTOR V": {"type": "gene", "canonical": "F5"},
    "FACTOR V LEIDEN": {"type": "gene", "canonical": "F5"},
    "RYANODINE RECEPTOR": {"type": "gene", "canonical": "RYR1"},
    "G6PD DEFICIENCY": {"type": "gene", "canonical": "G6PD"},
    "GLUCOSE-6-PHOSPHATE DEHYDROGENASE": {"type": "gene", "canonical": "G6PD"},
    "NAT-2": {"type": "gene", "canonical": "NAT2"},
    "N-ACETYLTRANSFERASE 2": {"type": "gene", "canonical": "NAT2"},
    "MTHFR GENE": {"type": "gene", "canonical": "MTHFR"},
    "METHYLENETETRAHYDROFOLATE REDUCTASE": {"type": "gene", "canonical": "MTHFR"},
    "UGT": {"type": "gene", "canonical": "UGT1A1"},
    "GILBERT SYNDROME GENE": {"type": "gene", "canonical": "UGT1A1"},
    # --- Phenotype aliases ---
    "PM": {"type": "phenotype", "canonical": "poor"},
    "POOR METABOLIZER": {"type": "phenotype", "canonical": "poor"},
    "IM": {"type": "phenotype", "canonical": "intermediate"},
    "INTERMEDIATE METABOLIZER": {"type": "phenotype", "canonical": "intermediate"},
    "NM": {"type": "phenotype", "canonical": "normal"},
    "EM": {"type": "phenotype", "canonical": "normal"},
    "NORMAL METABOLIZER": {"type": "phenotype", "canonical": "normal"},
    "EXTENSIVE METABOLIZER": {"type": "phenotype", "canonical": "normal"},
    "RM": {"type": "phenotype", "canonical": "rapid"},
    "RAPID METABOLIZER": {"type": "phenotype", "canonical": "rapid"},
    "UM": {"type": "phenotype", "canonical": "ultra_rapid"},
    "ULTRA-RAPID METABOLIZER": {"type": "phenotype", "canonical": "ultra_rapid"},
    "ULTRARAPID METABOLIZER": {"type": "phenotype", "canonical": "ultra_rapid"},
    "ULTRA RAPID METABOLIZER": {"type": "phenotype", "canonical": "ultra_rapid"},
    # --- Phenotype descriptors ---
    "SLOW METABOLIZER": {"type": "phenotype", "canonical": "poor"},
    "FAST METABOLIZER": {"type": "phenotype", "canonical": "ultra_rapid"},
    "SLOW ACETYLATOR": {"type": "phenotype", "canonical": "poor"},
    "RAPID ACETYLATOR": {"type": "phenotype", "canonical": "ultra_rapid"},
    "EXPRESSER": {"type": "phenotype", "canonical": "normal"},
    "NON-EXPRESSER": {"type": "phenotype", "canonical": "poor"},
}


# =============================================================================
# PEDIATRIC PHARMACOGENOMICS
# =============================================================================

PEDIATRIC_PGX: Dict[str, Dict[str, Any]] = {
    "TPMT_6MP": {
        "gene": "TPMT",
        "drug": "6-mercaptopurine",
        "indication": "Pediatric acute lymphoblastic leukemia (ALL) maintenance therapy",
        "dosing_by_phenotype": {
            "poor_metabolizer_TPMT_3A_3A": "Start at 10% of standard dose",
            "intermediate_TPMT_1_3A": "Start at 50% of standard dose",
            "normal_metabolizer": "Standard dosing with monitoring",
        },
        "monitoring": "Weekly CBC during dose titration; TGN levels for dose optimization",
        "guideline": "CPIC Guideline for Thiopurines and TPMT/NUDT15",
    },
    "NUDT15_6MP": {
        "gene": "NUDT15",
        "drug": "6-mercaptopurine",
        "indication": "Pediatric ALL maintenance therapy",
        "dosing_by_phenotype": {
            "poor_metabolizer_NUDT15_3_3": "Start at 10% of standard dose",
            "intermediate_NUDT15_1_3": "Start at 25-50% of standard dose",
        },
        "population_note": "Higher prevalence of risk alleles in East Asian populations "
                          "(NUDT15 *3 allele frequency ~10% in East Asian vs <1% in European)",
        "guideline": "CPIC Guideline for Thiopurines and TPMT/NUDT15",
    },
    "UGT1A1_irinotecan": {
        "gene": "UGT1A1",
        "drug": "Irinotecan",
        "indication": "Pediatric solid tumors (rhabdomyosarcoma, neuroblastoma, Ewing sarcoma)",
        "dosing_by_genotype": {
            "UGT1A1_28_28": "30% dose reduction from standard",
        },
        "key_toxicity": "Severe diarrhea and neutropenia; monitor closely in *28/*28",
        "guideline": "CPIC/DPWG Guideline for UGT1A1 and Irinotecan",
    },
    "CYP3A5_vincristine": {
        "gene": "CYP3A4/CYP3A5",
        "drug": "Vincristine",
        "indication": "Pediatric ALL, lymphoma, solid tumors",
        "dosing_by_genotype": {
            "CYP3A5_1_1_rapid_clearance": "May need dose escalation due to faster metabolism",
            "CYP3A5_3_3_standard": "Standard dosing",
        },
        "clinical_note": "CYP3A5 expressors (*1/*1) may have reduced vincristine exposure "
                        "and potentially less neurotoxicity but also less efficacy",
    },
    "DPYD_fluoropyrimidines": {
        "gene": "DPYD",
        "drug": "5-FU / capecitabine",
        "indication": "Rare in pediatrics; used in some solid tumors (hepatoblastoma, NPC)",
        "dosing_by_phenotype": {
            "poor_metabolizer": "Contraindicated (DPD deficiency → life-threatening toxicity)",
            "intermediate_metabolizer": "50% dose reduction",
        },
        "guideline": "CPIC Guideline for Fluoropyrimidines and DPYD",
    },
    "MTHFR_methotrexate": {
        "gene": "MTHFR",
        "drug": "Methotrexate",
        "indication": "Pediatric ALL (high-dose MTX protocol)",
        "variant": "C677T (rs1801133)",
        "clinical_impact": "Reduced folate metabolism → increased MTX toxicity risk "
                          "(mucositis, myelosuppression, leukoencephalopathy)",
        "management": "Enhanced leucovorin rescue; close monitoring of MTX clearance",
    },
    "asparaginase_immunogenicity": {
        "gene": "No established PGx markers",
        "drug": "Asparaginase (E. coli, Erwinia, PEG-asparaginase)",
        "indication": "Pediatric ALL induction and consolidation",
        "clinical_note": "Anti-asparaginase antibodies affect efficacy; silent inactivation "
                        "in up to 30% of patients. Monitor trough asparaginase activity levels. "
                        "Switch to Erwinia or recombinant formulation if antibodies detected.",
    },
}


# =============================================================================
# PUBLIC FUNCTIONS — Knowledge graph query interface
# =============================================================================

def get_gene_context(gene: str) -> str:
    """Return formatted pharmacogene knowledge for a given gene.

    Args:
        gene: Gene symbol (e.g. 'CYP2D6', 'DPYD').

    Returns:
        Formatted string with gene knowledge, or empty string if not found.
    """
    key = gene.strip().upper()
    # Check alias table
    if key in ENTITY_ALIASES and ENTITY_ALIASES[key].get("type") == "gene":
        key = ENTITY_ALIASES[key]["canonical"]

    data = PHARMACOGENES.get(key)
    if not data:
        # Try partial match
        for k in PHARMACOGENES:
            if key in k.upper():
                data = PHARMACOGENES[k]
                key = k
                break
    if not data:
        return ""

    lines = [
        f"## Pharmacogene: {key} ({data['full_name']})",
        f"- Chromosome: {data['chromosome']}",
        f"- Function: {data['function']}",
    ]
    if data.get("substrates_count"):
        lines.append(f"- Known substrates: ~{data['substrates_count']}")
    if data.get("percent_drugs_metabolized"):
        lines.append(f"- Percent of drugs metabolized: ~{data['percent_drugs_metabolized']}%")
    if data.get("star_alleles_defined"):
        lines.append(f"- Star alleles defined: {data['star_alleles_defined']}")

    lines.append(f"- Key variants: {', '.join(data['key_variants'])}")
    lines.append(f"- Structural variation: {'Yes' if data['structural_variation'] else 'No'}")
    lines.append(f"- Complexity: {data['complexity_level']}")

    if data.get("cpic_guidelines"):
        lines.append(f"- CPIC guidelines exist for: {', '.join(data['cpic_guidelines'])}")

    # Add activity score info if available
    if key in ACTIVITY_SCORE_TABLES:
        ast = ACTIVITY_SCORE_TABLES[key]
        lines.append(f"\n### Activity Score System ({key})")
        lines.append(f"- {ast['description']}")
        for pheno, threshold in ast["phenotype_thresholds"].items():
            lines.append(f"  - {pheno}: {threshold}")

    return "\n".join(lines)


def get_phenotype_context(phenotype: str) -> str:
    """Return formatted metabolizer phenotype knowledge.

    Args:
        phenotype: Phenotype name or abbreviation (e.g. 'PM', 'poor', 'Ultra-rapid').

    Returns:
        Formatted string with phenotype information.
    """
    key = phenotype.strip().upper()
    # Check alias table
    if key in ENTITY_ALIASES and ENTITY_ALIASES[key].get("type") == "phenotype":
        key = ENTITY_ALIASES[key]["canonical"]
    else:
        key = phenotype.strip().lower().replace("-", "_").replace(" ", "_")

    data = METABOLIZER_PHENOTYPES.get(key)
    if not data:
        for k, v in METABOLIZER_PHENOTYPES.items():
            if key.upper() == v["abbreviation"] or key.lower() in k:
                data = v
                key = k
                break
    if not data:
        return ""

    lines = [
        f"## Metabolizer Phenotype: {key.replace('_', ' ').title()} ({data['abbreviation']})",
        f"- Clinical meaning: {data['clinical_meaning']}",
        f"- Risk profile: {data['risk']}",
    ]
    return "\n".join(lines)


def get_drug_category_context(category: str) -> str:
    """Return formatted drug category knowledge.

    Args:
        category: Category name (e.g. 'opioids', 'statins', 'antidepressants').

    Returns:
        Formatted string with drug category information.
    """
    key = category.strip().lower().replace(" ", "_").replace("-", "_")
    data = DRUG_CATEGORIES.get(key)
    if not data:
        for k in DRUG_CATEGORIES:
            if key in k.lower() or k.lower() in key:
                data = DRUG_CATEGORIES[k]
                key = k
                break
    if not data:
        return ""

    lines = [
        f"## Drug Category: {key.replace('_', ' ').title()}",
        f"- Description: {data['description']}",
        f"- Drugs: {', '.join(data['drugs'])}",
        f"- Primary pharmacogenes: {', '.join(data['primary_genes'])}",
    ]
    return "\n".join(lines)


def get_hla_context(hla_allele: str) -> str:
    """Return formatted HLA-drug association knowledge.

    Args:
        hla_allele: HLA allele (e.g. 'HLA-B*57:01') or drug name.

    Returns:
        Formatted string with all HLA-drug associations matching the query.
    """
    query = hla_allele.strip().upper()
    matches = []

    for key, data in HLA_DRUG_ASSOCIATIONS.items():
        if (query in data["hla_allele"].upper() or
                query in data["drug"].upper() or
                query in key.upper()):
            matches.append(data)

    if not matches:
        return ""

    lines = []
    for data in matches:
        lines.append(f"## HLA-Drug Association: {data['hla_allele']} / {data['drug']}")
        lines.append(f"- Reaction type: {data['reaction_type']}")
        lines.append(f"- Risk if positive: {data['risk_if_positive']}")
        lines.append(f"- Severity: {data['severity']}")
        lines.append(f"- Screening mandatory: {'Yes' if data['screening_mandatory'] else 'No'}")
        lines.append(f"- Recommendation: {data['recommendation']}")
        prev_parts = [f"{pop}: {freq}" for pop, freq in data["prevalence_by_population"].items()]
        lines.append(f"- Prevalence: {'; '.join(prev_parts)}")
        lines.append("")

    return "\n".join(lines)


def get_inhibitor_context(enzyme: str) -> str:
    """Return formatted CYP inhibitor/inducer knowledge for an enzyme.

    Args:
        enzyme: Enzyme name (e.g. 'CYP2D6', 'CYP3A4').

    Returns:
        Formatted string with inhibitor and inducer lists.
    """
    key = enzyme.strip().upper()
    lines = []

    if key in CYP_INHIBITORS:
        inh = CYP_INHIBITORS[key]
        lines.append(f"## {key} Inhibitors")
        lines.append(f"- Strong: {', '.join(inh['strong'])}")
        lines.append(f"- Moderate: {', '.join(inh['moderate'])}")
        if "weak" in inh:
            lines.append(f"- Weak: {', '.join(inh['weak'])}")

    if key in CYP_INDUCERS:
        ind = CYP_INDUCERS[key]
        lines.append(f"\n## {key} Inducers")
        lines.append(f"- Strong: {', '.join(ind['strong'])}")
        lines.append(f"- Moderate: {', '.join(ind['moderate'])}")

    return "\n".join(lines) if lines else ""


def get_alternative_drugs(drug: str, gene: str, phenotype: str) -> str:
    """Return formatted drug alternative recommendations.

    Args:
        drug: Drug name (generic or brand).
        gene: Gene symbol.
        phenotype: Metabolizer phenotype.

    Returns:
        Formatted string with alternative drug recommendations.
    """
    # Resolve brand names
    drug_clean = drug.strip().upper()
    if drug_clean in ENTITY_ALIASES and ENTITY_ALIASES[drug_clean].get("type") == "drug":
        drug_clean = ENTITY_ALIASES[drug_clean]["canonical"]
    else:
        drug_clean = drug.strip().lower()

    gene_clean = gene.strip().upper()
    if gene_clean in ENTITY_ALIASES and ENTITY_ALIASES[gene_clean].get("type") == "gene":
        gene_clean = ENTITY_ALIASES[gene_clean]["canonical"]

    # Search for matching alternatives
    matches = []
    for key, data in DRUG_ALTERNATIVES.items():
        key_lower = key.lower()
        if (drug_clean.lower() in key_lower and
                gene_clean.lower() in key_lower):
            matches.append(data)

    if not matches:
        # Broader search — drug name only
        for key, data in DRUG_ALTERNATIVES.items():
            if drug_clean.lower() in data["drug"].lower():
                matches.append(data)

    if not matches:
        return ""

    lines = []
    for data in matches:
        lines.append(f"## Alternative Drugs: {data['drug']} ({data['gene']} {data['phenotype']})")
        lines.append(f"- Rationale: {data['rationale']}")
        lines.append("- Alternatives:")
        for alt in data["alternatives"]:
            lines.append(f"  - **{alt['drug']}**: {alt['rationale']}")
        lines.append("")

    return "\n".join(lines)


def get_drug_context(drug: str) -> str:
    """Return formatted knowledge context for a specific drug.

    Searches across drug categories and alternative mappings to find
    the drug and return its pharmacogenomic context.

    Args:
        drug: Drug name (generic or brand).

    Returns:
        Formatted string with drug PGx context.
    """
    drug_clean = drug.strip().lower()

    # Resolve brand names via aliases
    drug_upper = drug.strip().upper()
    if drug_upper in ENTITY_ALIASES and ENTITY_ALIASES[drug_upper].get("type") == "drug":
        drug_clean = ENTITY_ALIASES[drug_upper]["canonical"].lower()

    lines = []

    # Find in drug categories
    for cat_key, cat_data in DRUG_CATEGORIES.items():
        if drug_clean in [d.lower() for d in cat_data["drugs"]]:
            lines.append(f"## Drug: {drug_clean.title()}")
            lines.append(f"- Category: {cat_key.replace('_', ' ').title()}")
            lines.append(f"- Primary pharmacogenes: {', '.join(cat_data['primary_genes'])}")
            break

    # Find in alternative drug mappings
    for key, data in DRUG_ALTERNATIVES.items():
        if drug_clean == data["drug"].lower():
            lines.append(f"- Gene: {data['gene']}")
            lines.append(f"- Affected phenotype: {data['phenotype']}")
            lines.append(f"- Rationale: {data['rationale']}")
            alts = ", ".join(a["drug"] for a in data["alternatives"])
            lines.append(f"- Alternatives: {alts}")
            break

    return "\n".join(lines)


def get_all_context_for_query(query: str) -> str:
    """Build comprehensive knowledge context for a free-text query.

    Scans the query for gene names, drug names, phenotypes, and HLA alleles,
    then aggregates all relevant knowledge graph entries.

    Args:
        query: Free-text clinical question.

    Returns:
        Aggregated knowledge context string.
    """
    query_upper = query.upper()
    sections: List[str] = []

    # 1. Check for pharmacogenes
    for gene in PHARMACOGENES:
        if gene.upper() in query_upper:
            ctx = get_gene_context(gene)
            if ctx:
                sections.append(ctx)

    # 2. Check for gene aliases
    for alias, info in ENTITY_ALIASES.items():
        if info.get("type") == "gene" and alias in query_upper:
            ctx = get_gene_context(info["canonical"])
            if ctx and ctx not in sections:
                sections.append(ctx)

    # 3. Check for metabolizer phenotypes
    for pheno in METABOLIZER_PHENOTYPES:
        abbrev = METABOLIZER_PHENOTYPES[pheno]["abbreviation"]
        pheno_words = pheno.replace("_", " ")
        if abbrev in query_upper or pheno_words.upper() in query_upper:
            ctx = get_phenotype_context(pheno)
            if ctx:
                sections.append(ctx)

    # 4. Check for drug categories
    for cat in DRUG_CATEGORIES:
        cat_words = cat.replace("_", " ")
        if cat_words.upper() in query_upper:
            ctx = get_drug_category_context(cat)
            if ctx:
                sections.append(ctx)

    # 5. Check for HLA alleles
    hla_patterns = [
        "HLA-B*57:01", "HLA-B*15:02", "HLA-B*58:01", "HLA-A*31:01",
        "HLA-B*15:11", "HLA-A*02:01", "HLA-B*35:05", "HLA-DRB1*07:01",
        "HLA-DPB1*03:01",
    ]
    for hla in hla_patterns:
        if hla.upper().replace(":", "").replace("*", "") in query_upper.replace(":", "").replace("*", ""):
            ctx = get_hla_context(hla)
            if ctx and ctx not in sections:
                sections.append(ctx)

    # 6. Check for specific drugs by brand or generic name
    for alias, info in ENTITY_ALIASES.items():
        if info.get("type") == "drug" and alias in query_upper:
            generic = info["canonical"]
            # Try to find alternatives
            for key in DRUG_ALTERNATIVES:
                if generic.lower() in key.lower():
                    alt_ctx = get_alternative_drugs(generic, "", "")
                    if alt_ctx and alt_ctx not in sections:
                        sections.append(alt_ctx)
                    break

    # 7. Check for drug generic names in alternatives
    for key, data in DRUG_ALTERNATIVES.items():
        if data["drug"].upper() in query_upper:
            ctx = get_alternative_drugs(data["drug"], data["gene"], data["phenotype"])
            if ctx and ctx not in sections:
                sections.append(ctx)

    # 8. Check for CYP inhibitor/inducer queries
    inhibitor_keywords = ["INHIBITOR", "INDUCER", "INTERACTION", "DDI", "DRUG INTERACTION"]
    if any(kw in query_upper for kw in inhibitor_keywords):
        for enzyme in CYP_INHIBITORS:
            if enzyme in query_upper:
                ctx = get_inhibitor_context(enzyme)
                if ctx and ctx not in sections:
                    sections.append(ctx)

    return "\n\n---\n\n".join(sections) if sections else ""


def get_knowledge_stats() -> Dict[str, int]:
    """Return counts of all knowledge graph entries.

    Returns:
        Dictionary with counts for each knowledge category.
    """
    return {
        "pharmacogenes": len(PHARMACOGENES),
        "metabolizer_phenotypes": len(METABOLIZER_PHENOTYPES),
        "drug_categories": len(DRUG_CATEGORIES),
        "drugs_tracked": sum(len(v["drugs"]) for v in DRUG_CATEGORIES.values()),
        "cyp_inhibitor_enzymes": len(CYP_INHIBITORS),
        "cyp_inducer_enzymes": len(CYP_INDUCERS),
        "hla_drug_associations": len(HLA_DRUG_ASSOCIATIONS),
        "drug_alternative_mappings": len(DRUG_ALTERNATIVES),
        "activity_score_tables": len(ACTIVITY_SCORE_TABLES),
        "entity_aliases": len(ENTITY_ALIASES),
    }


def resolve_comparison_entity(text: str) -> Optional[Dict[str, str]]:
    """Resolve a raw text string to a known pharmacogenomics entity.

    Checks gene names, drug aliases, phenotypes, and HLA alleles
    in priority order.

    Args:
        text: Raw text to resolve (e.g. 'Plavix', 'CYP2D6', 'PM').

    Returns:
        Dict with 'type' and 'canonical' keys, or None if not recognized.
    """
    cleaned = text.strip().upper()

    # 1. Exact gene match
    for key in PHARMACOGENES:
        if key.upper() == cleaned:
            return {"type": "gene", "canonical": key}

    # 2. Entity alias table
    if cleaned in ENTITY_ALIASES:
        return dict(ENTITY_ALIASES[cleaned])

    # 3. HLA alleles
    for key in HLA_DRUG_ASSOCIATIONS:
        data = HLA_DRUG_ASSOCIATIONS[key]
        hla_norm = data["hla_allele"].upper().replace(":", "").replace("*", "")
        cleaned_norm = cleaned.replace(":", "").replace("*", "").replace("-", "")
        if hla_norm.replace("-", "") == cleaned_norm:
            return {"type": "hla", "canonical": data["hla_allele"]}

    # 4. Drug categories
    for key in DRUG_CATEGORIES:
        if cleaned.lower().replace(" ", "_") == key:
            return {"type": "drug_category", "canonical": key}

    # 5. Phenotypes (substring)
    for key in METABOLIZER_PHENOTYPES:
        abbrev = METABOLIZER_PHENOTYPES[key]["abbreviation"]
        if cleaned == abbrev or cleaned.lower().replace(" ", "_") == key:
            return {"type": "phenotype", "canonical": key}

    # 6. Drug generic names in alternatives
    for alt_key, data in DRUG_ALTERNATIVES.items():
        if cleaned.lower() == data["drug"].lower():
            return {"type": "drug", "canonical": data["drug"]}

    return None


def get_comparison_context(entity_a: Dict[str, str], entity_b: Dict[str, str]) -> str:
    """Build side-by-side knowledge graph context for two entities.

    Reuses existing get_gene_context / get_hla_context /
    get_phenotype_context / get_drug_category_context depending on type.

    Args:
        entity_a: Resolved entity dict with 'type' and 'canonical'.
        entity_b: Resolved entity dict with 'type' and 'canonical'.

    Returns:
        Formatted comparison context string with both entities' data.
    """
    def _get_entity_context(entity: Dict[str, str]) -> str:
        etype = entity["type"]
        canonical = entity["canonical"]
        if etype == "gene":
            return get_gene_context(canonical)
        elif etype == "drug":
            # Return all alternative mappings for this drug
            ctx_parts = []
            for key, data in DRUG_ALTERNATIVES.items():
                if canonical.lower() in data["drug"].lower():
                    ctx_parts.append(
                        get_alternative_drugs(data["drug"], data["gene"], data["phenotype"])
                    )
            return "\n".join(ctx_parts) if ctx_parts else f"Drug: {canonical}"
        elif etype == "hla":
            return get_hla_context(canonical)
        elif etype == "phenotype":
            return get_phenotype_context(canonical)
        elif etype == "drug_category":
            return get_drug_category_context(canonical)
        return ""

    sections = []
    ctx_a = _get_entity_context(entity_a)
    ctx_b = _get_entity_context(entity_b)

    if ctx_a:
        sections.append(f"### {entity_a['canonical']}\n{ctx_a}")
    if ctx_b:
        sections.append(f"### {entity_b['canonical']}\n{ctx_b}")

    return "\n\n---\n\n".join(sections)
