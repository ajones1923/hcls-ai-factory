"""Oncology Intelligence knowledge graph — targets, therapies, resistance, pathways, biomarkers.

Provides structured domain knowledge for augmenting RAG queries with
curated oncology intelligence. Each helper function returns formatted
context strings that are injected into LLM prompts alongside retrieved
evidence from Milvus collections.

Follows the same pattern as:
  cart_intelligence_agent/src/knowledge.py

Author: Adam Jones
Date: February 2026
"""

from typing import Dict, List, Optional


# ═══════════════════════════════════════════════════════════════════════
# 1. ACTIONABLE TARGETS (~40 genes)
# ═══════════════════════════════════════════════════════════════════════

ACTIONABLE_TARGETS: Dict[str, Dict] = {
    "BRAF": {
        "gene": "BRAF",
        "full_name": "B-Raf Proto-Oncogene",
        "cancer_types": ["melanoma", "NSCLC", "colorectal", "thyroid", "hairy cell leukemia"],
        "key_variants": ["V600E", "V600K", "V600D", "V600R", "class II", "class III"],
        "targeted_therapies": ["vemurafenib", "dabrafenib", "encorafenib"],
        "combination_therapies": ["dabrafenib + trametinib", "encorafenib + binimetinib", "encorafenib + cetuximab"],
        "resistance_mutations": ["MEK1/2 mutations", "NRAS activation", "BRAF amplification", "MAP2K1 C121S"],
        "pathway": "MAPK",
        "evidence_level": "A",
        "description": "BRAF V600E is the most common BRAF mutation, found in ~50% of melanomas. "
                       "FDA-approved targeted therapies include BRAF+MEK inhibitor combinations. "
                       "In CRC, BRAF V600E requires triplet therapy (BRAF+MEK+EGFR inhibition).",
    },
    "EGFR": {
        "gene": "EGFR",
        "full_name": "Epidermal Growth Factor Receptor",
        "cancer_types": ["NSCLC", "head and neck", "colorectal", "glioblastoma"],
        "key_variants": ["L858R", "exon 19 deletion", "T790M", "C797S", "exon 20 insertion", "S768I", "L861Q", "G719X"],
        "targeted_therapies": ["osimertinib", "erlotinib", "gefitinib", "afatinib", "dacomitinib", "amivantamab"],
        "combination_therapies": ["osimertinib + chemotherapy", "amivantamab + lazertinib"],
        "resistance_mutations": ["T790M", "C797S", "MET amplification", "HER2 amplification",
                                  "small cell transformation", "BRAF V600E", "PIK3CA mutations"],
        "pathway": "MAPK",
        "evidence_level": "A",
        "description": "EGFR mutations occur in ~15-20% of NSCLC (higher in never-smokers, Asian populations). "
                       "Osimertinib is standard first-line for EGFR-mutant NSCLC. T790M is the most common "
                       "acquired resistance mechanism to 1st/2nd-gen TKIs.",
    },
    "ALK": {
        "gene": "ALK",
        "full_name": "Anaplastic Lymphoma Kinase",
        "cancer_types": ["NSCLC", "anaplastic large cell lymphoma", "neuroblastoma"],
        "key_variants": ["EML4-ALK fusion", "ALK amplification", "F1174L", "R1275Q"],
        "targeted_therapies": ["alectinib", "lorlatinib", "brigatinib", "ceritinib", "crizotinib"],
        "combination_therapies": [],
        "resistance_mutations": ["G1202R", "L1196M", "I1171T/N/S", "compound mutations",
                                  "MET amplification", "EGFR activation"],
        "pathway": "MAPK/PI3K",
        "evidence_level": "A",
        "description": "ALK fusions occur in ~5% of NSCLC. Alectinib is preferred first-line. "
                       "Lorlatinib is effective against most single ALK resistance mutations "
                       "including G1202R. Sequential ALK TKI therapy can extend disease control.",
    },
    "ROS1": {
        "gene": "ROS1",
        "full_name": "ROS Proto-Oncogene 1",
        "cancer_types": ["NSCLC", "cholangiocarcinoma", "glioblastoma"],
        "key_variants": ["CD74-ROS1 fusion", "SLC34A2-ROS1 fusion", "SDC4-ROS1 fusion"],
        "targeted_therapies": ["crizotinib", "entrectinib", "lorlatinib", "repotrectinib"],
        "combination_therapies": [],
        "resistance_mutations": ["G2032R", "D2033N", "L2026M"],
        "pathway": "MAPK/PI3K",
        "evidence_level": "A",
        "description": "ROS1 fusions occur in ~1-2% of NSCLC. Entrectinib and crizotinib are "
                       "first-line options. Repotrectinib shows activity against G2032R resistance.",
    },
    "KRAS": {
        "gene": "KRAS",
        "full_name": "KRAS Proto-Oncogene GTPase",
        "cancer_types": ["NSCLC", "colorectal", "pancreatic"],
        "key_variants": ["G12C", "G12D", "G12V", "G13D", "Q61H"],
        "targeted_therapies": ["sotorasib", "adagrasib"],
        "combination_therapies": ["sotorasib + panitumumab", "adagrasib + cetuximab"],
        "resistance_mutations": ["Y96D", "R68S", "H95D/Q/R", "secondary KRAS mutations",
                                  "MET amplification", "EGFR activation", "RAS-MAPK reactivation"],
        "pathway": "MAPK",
        "evidence_level": "A",
        "description": "KRAS G12C occurs in ~13% of NSCLC and ~3% of CRC. Sotorasib and adagrasib "
                       "are FDA-approved for KRAS G12C-mutant NSCLC. Other KRAS alleles (G12D) are "
                       "targets of emerging therapies including MRTX1133.",
    },
    "HER2": {
        "gene": "ERBB2",
        "full_name": "Human Epidermal Growth Factor Receptor 2",
        "cancer_types": ["breast", "gastric", "NSCLC", "colorectal", "endometrial"],
        "key_variants": ["amplification", "overexpression", "S310F", "exon 20 insertions", "V777L"],
        "targeted_therapies": ["trastuzumab", "pertuzumab", "trastuzumab deruxtecan", "tucatinib",
                                "trastuzumab emtansine", "margetuximab"],
        "combination_therapies": ["trastuzumab + pertuzumab + docetaxel", "tucatinib + trastuzumab + capecitabine"],
        "resistance_mutations": ["PIK3CA mutations", "PTEN loss", "HER2 truncation (p95)", "NRG1 activation"],
        "pathway": "PI3K/MAPK",
        "evidence_level": "A",
        "description": "HER2 amplification defines a breast cancer subtype (~20%). Trastuzumab deruxtecan "
                       "(T-DXd/Enhertu) has shown activity across HER2-amplified and HER2-low cancers. "
                       "HER2 mutations in NSCLC are targetable with T-DXd.",
    },
    "NTRK": {
        "gene": "NTRK1/2/3",
        "full_name": "Neurotrophic Tyrosine Receptor Kinase",
        "cancer_types": ["tissue-agnostic", "secretory breast", "infantile fibrosarcoma", "thyroid",
                          "salivary gland", "NSCLC", "CRC"],
        "key_variants": ["ETV6-NTRK3 fusion", "TPM3-NTRK1 fusion", "LMNA-NTRK1 fusion"],
        "targeted_therapies": ["larotrectinib", "entrectinib", "repotrectinib", "selitrectinib"],
        "combination_therapies": [],
        "resistance_mutations": ["G595R (NTRK1)", "G623R (NTRK3)", "F589L", "xDFG mutations"],
        "pathway": "MAPK/PI3K",
        "evidence_level": "A",
        "description": "NTRK fusions are rare (<1% pan-cancer) but highly actionable. Larotrectinib "
                       "and entrectinib have tissue-agnostic FDA approvals. ORR >75% across tumor types.",
    },
    "RET": {
        "gene": "RET",
        "full_name": "Rearranged during Transfection",
        "cancer_types": ["NSCLC", "medullary thyroid", "papillary thyroid"],
        "key_variants": ["KIF5B-RET fusion", "CCDC6-RET fusion", "M918T", "C634R"],
        "targeted_therapies": ["selpercatinib", "pralsetinib"],
        "combination_therapies": [],
        "resistance_mutations": ["G810R/S/C", "V804M", "Y806C/H"],
        "pathway": "MAPK/PI3K",
        "evidence_level": "A",
        "description": "RET fusions occur in ~1-2% of NSCLC and RET mutations in >90% of medullary "
                       "thyroid cancer. Selpercatinib is the preferred selective RET inhibitor.",
    },
    "MET": {
        "gene": "MET",
        "full_name": "MET Proto-Oncogene (Hepatocyte Growth Factor Receptor)",
        "cancer_types": ["NSCLC", "renal cell", "hepatocellular", "gastric"],
        "key_variants": ["exon 14 skipping", "amplification", "Y1003 mutations"],
        "targeted_therapies": ["capmatinib", "tepotinib", "crizotinib"],
        "combination_therapies": ["osimertinib + savolitinib"],
        "resistance_mutations": ["D1228N/V", "Y1230H/C", "KRAS amplification"],
        "pathway": "MAPK/PI3K",
        "evidence_level": "A",
        "description": "MET exon 14 skipping occurs in ~3-4% of NSCLC. MET amplification is a key "
                       "resistance mechanism to EGFR TKIs. Capmatinib and tepotinib are FDA-approved.",
    },
    "FGFR": {
        "gene": "FGFR1/2/3",
        "full_name": "Fibroblast Growth Factor Receptor",
        "cancer_types": ["bladder/urothelial", "cholangiocarcinoma", "endometrial", "gastric"],
        "key_variants": ["FGFR2 fusions", "FGFR3 S249C", "FGFR3 R248C", "FGFR2 amplification"],
        "targeted_therapies": ["erdafitinib", "futibatinib", "pemigatinib", "infigratinib"],
        "combination_therapies": [],
        "resistance_mutations": ["V565I (gatekeeper)", "N550H/K", "FGFR2 V564F"],
        "pathway": "MAPK/PI3K",
        "evidence_level": "A",
        "description": "FGFR alterations are common in bladder cancer (~20%) and cholangiocarcinoma "
                       "(~15% FGFR2 fusions). Erdafitinib is approved for FGFR-altered urothelial cancer.",
    },
    "PIK3CA": {
        "gene": "PIK3CA",
        "full_name": "Phosphatidylinositol-4,5-Bisphosphate 3-Kinase Catalytic Subunit Alpha",
        "cancer_types": ["breast", "endometrial", "cervical", "head and neck", "colorectal"],
        "key_variants": ["H1047R", "E545K", "E542K", "C420R", "N345K"],
        "targeted_therapies": ["alpelisib", "inavolisib"],
        "combination_therapies": ["alpelisib + fulvestrant", "inavolisib + palbociclib + fulvestrant"],
        "resistance_mutations": ["PTEN loss", "AKT1 E17K", "mTOR mutations"],
        "pathway": "PI3K/AKT/mTOR",
        "evidence_level": "A",
        "description": "PIK3CA mutations occur in ~40% of HR+/HER2- breast cancer. Alpelisib + fulvestrant "
                       "is approved for PIK3CA-mutant advanced breast cancer after progression on endocrine therapy.",
    },
    "IDH1": {
        "gene": "IDH1",
        "full_name": "Isocitrate Dehydrogenase 1",
        "cancer_types": ["glioma", "AML", "cholangiocarcinoma"],
        "key_variants": ["R132H", "R132C", "R132G", "R132S"],
        "targeted_therapies": ["ivosidenib", "vorasidenib"],
        "combination_therapies": ["ivosidenib + azacitidine"],
        "resistance_mutations": ["second-site IDH1 mutations", "RTK pathway activation"],
        "pathway": "metabolic",
        "evidence_level": "A",
        "description": "IDH1 R132H is the most common mutation. Ivosidenib is approved for IDH1-mutant "
                       "AML and cholangiocarcinoma. Vorasidenib (brain-penetrant dual IDH1/2) approved for "
                       "low-grade glioma.",
    },
    "IDH2": {
        "gene": "IDH2",
        "full_name": "Isocitrate Dehydrogenase 2",
        "cancer_types": ["AML", "angioimmunoblastic T-cell lymphoma"],
        "key_variants": ["R140Q", "R172K"],
        "targeted_therapies": ["enasidenib"],
        "combination_therapies": ["enasidenib + azacitidine"],
        "resistance_mutations": ["second-site IDH2 mutations", "RAS pathway activation"],
        "pathway": "metabolic",
        "evidence_level": "A",
        "description": "IDH2 mutations occur in ~10-15% of AML. Enasidenib is approved for relapsed/refractory "
                       "IDH2-mutant AML. R140Q is more common than R172K.",
    },
    "BRCA1": {
        "gene": "BRCA1",
        "full_name": "BRCA1 DNA Repair Associated",
        "cancer_types": ["breast", "ovarian", "prostate", "pancreatic"],
        "key_variants": ["various truncating/missense mutations", "185delAG", "5382insC",
                          "large rearrangements", "promoter methylation (somatic)"],
        "targeted_therapies": ["olaparib", "rucaparib", "niraparib", "talazoparib"],
        "combination_therapies": ["olaparib + bevacizumab", "talazoparib + enzalutamide"],
        "resistance_mutations": ["reversion mutations", "53BP1 loss", "BRCA1 promoter demethylation",
                                  "RAD51 upregulation", "drug efflux (ABCB1)"],
        "pathway": "DDR",
        "evidence_level": "A",
        "description": "BRCA1 mutations confer homologous recombination deficiency (HRD), creating synthetic "
                       "lethality with PARP inhibitors. Germline testing recommended for breast, ovarian, "
                       "prostate, and pancreatic cancer.",
    },
    "BRCA2": {
        "gene": "BRCA2",
        "full_name": "BRCA2 DNA Repair Associated",
        "cancer_types": ["breast", "ovarian", "prostate", "pancreatic"],
        "key_variants": ["various truncating/missense mutations", "6174delT"],
        "targeted_therapies": ["olaparib", "rucaparib", "niraparib", "talazoparib"],
        "combination_therapies": ["olaparib + abiraterone"],
        "resistance_mutations": ["reversion mutations", "RAD51C/D mutations", "drug efflux"],
        "pathway": "DDR",
        "evidence_level": "A",
        "description": "BRCA2 mutations create HRD and PARP inhibitor sensitivity. More common than BRCA1 "
                       "in prostate and pancreatic cancers. Platinum sensitivity often correlates.",
    },
    "TP53": {
        "gene": "TP53",
        "full_name": "Tumor Protein P53",
        "cancer_types": ["pan-cancer"],
        "key_variants": ["R175H", "R248W", "R273H", "G245S", "R249S", "Y220C", "hotspot missense"],
        "targeted_therapies": [],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "cell_cycle/apoptosis",
        "evidence_level": "D",
        "description": "TP53 is the most frequently mutated gene in cancer (~50% of all cancers). "
                       "Currently no direct targeted therapies are FDA-approved, though Y220C-targeting "
                       "compounds are in development. TP53 mutations are prognostic in most cancer types.",
    },
    "PTEN": {
        "gene": "PTEN",
        "full_name": "Phosphatase and Tensin Homolog",
        "cancer_types": ["prostate", "endometrial", "glioblastoma", "breast"],
        "key_variants": ["deletion", "truncating mutations", "promoter methylation"],
        "targeted_therapies": [],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "PI3K/AKT/mTOR",
        "evidence_level": "C",
        "description": "PTEN loss activates PI3K/AKT signaling. No direct PTEN-restoring therapies exist, "
                       "but PTEN loss may predict response to PI3K/AKT/mTOR inhibitors. Common in prostate "
                       "and endometrial cancers.",
    },
    "CDKN2A": {
        "gene": "CDKN2A",
        "full_name": "Cyclin Dependent Kinase Inhibitor 2A (p16/p14ARF)",
        "cancer_types": ["melanoma", "NSCLC", "pancreatic", "mesothelioma", "glioma"],
        "key_variants": ["homozygous deletion", "promoter methylation", "missense mutations"],
        "targeted_therapies": ["palbociclib", "ribociclib", "abemaciclib"],
        "combination_therapies": ["CDK4/6i + endocrine therapy"],
        "resistance_mutations": ["RB1 loss", "CDK6 amplification", "CCNE1 amplification"],
        "pathway": "cell_cycle",
        "evidence_level": "B",
        "description": "CDKN2A loss removes the p16-mediated brake on CDK4/6, making tumors potentially "
                       "sensitive to CDK4/6 inhibitors. Homozygous deletion is a poor prognostic marker "
                       "in mesothelioma and glioma.",
    },
    "STK11": {
        "gene": "STK11",
        "full_name": "Serine/Threonine Kinase 11 (LKB1)",
        "cancer_types": ["NSCLC", "cervical"],
        "key_variants": ["truncating mutations", "deletions"],
        "targeted_therapies": [],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "AMPK/mTOR",
        "evidence_level": "C",
        "description": "STK11/LKB1 loss in NSCLC (~15-20%) is associated with primary resistance to "
                       "PD-1/PD-L1 immunotherapy, especially in KRAS co-mutated tumors. No direct targeted "
                       "therapy, but important for treatment selection.",
    },
    "ESR1": {
        "gene": "ESR1",
        "full_name": "Estrogen Receptor 1",
        "cancer_types": ["breast"],
        "key_variants": ["Y537S", "D538G", "E380Q", "L536R", "ligand-binding domain mutations"],
        "targeted_therapies": ["elacestrant", "fulvestrant"],
        "combination_therapies": ["elacestrant monotherapy"],
        "resistance_mutations": ["compound ESR1 mutations", "RB1 loss", "MYC amplification"],
        "pathway": "estrogen_signaling",
        "evidence_level": "A",
        "description": "ESR1 mutations (~30-40% in metastatic HR+ breast cancer) confer resistance to "
                       "aromatase inhibitors. Elacestrant (oral SERD) is FDA-approved for ESR1-mutant "
                       "advanced breast cancer after prior endocrine therapy + CDK4/6i.",
    },
    "APC": {
        "gene": "APC",
        "full_name": "Adenomatous Polyposis Coli",
        "cancer_types": ["colorectal"],
        "key_variants": ["truncating mutations", "R1450*", "mutation cluster region"],
        "targeted_therapies": [],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "WNT",
        "evidence_level": "D",
        "description": "APC mutations initiate >80% of colorectal cancers through WNT pathway activation. "
                       "No direct targeted therapies, but foundational for CRC pathogenesis understanding. "
                       "Germline APC mutations cause familial adenomatous polyposis (FAP).",
    },
    "PALB2": {
        "gene": "PALB2",
        "full_name": "Partner and Localizer of BRCA2",
        "cancer_types": ["breast", "pancreatic", "ovarian"],
        "key_variants": ["truncating mutations"],
        "targeted_therapies": ["olaparib", "rucaparib"],
        "combination_therapies": [],
        "resistance_mutations": ["reversion mutations"],
        "pathway": "DDR",
        "evidence_level": "B",
        "description": "PALB2 mutations confer HRD similar to BRCA1/2. Classified as a high-penetrance "
                       "breast cancer gene. PARP inhibitor sensitivity demonstrated in clinical trials.",
    },
    "ATM": {
        "gene": "ATM",
        "full_name": "Ataxia Telangiectasia Mutated",
        "cancer_types": ["prostate", "pancreatic", "breast", "CLL"],
        "key_variants": ["truncating mutations", "missense variants (VUS common)"],
        "targeted_therapies": ["olaparib"],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "DDR",
        "evidence_level": "C",
        "description": "ATM loss impairs DNA damage response. Evidence for PARP inhibitor benefit is less "
                       "robust than BRCA1/2 — primarily studied in mCRPC (PROfound trial ATM subgroup).",
    },
    "POLE": {
        "gene": "POLE",
        "full_name": "DNA Polymerase Epsilon Catalytic Subunit",
        "cancer_types": ["endometrial", "colorectal"],
        "key_variants": ["P286R", "V411L", "S459F", "exonuclease domain mutations"],
        "targeted_therapies": [],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "DDR",
        "evidence_level": "B",
        "description": "POLE exonuclease domain mutations create ultramutated tumors (>100 mut/Mb) with "
                       "excellent prognosis and high immunotherapy response rates. Important for treatment "
                       "de-escalation decisions in endometrial cancer.",
    },
    "ARID1A": {
        "gene": "ARID1A",
        "full_name": "AT-Rich Interaction Domain 1A",
        "cancer_types": ["ovarian (clear cell)", "gastric", "endometrial", "cholangiocarcinoma"],
        "key_variants": ["truncating mutations", "frameshifts"],
        "targeted_therapies": [],
        "combination_therapies": [],
        "resistance_mutations": [],
        "pathway": "chromatin_remodeling",
        "evidence_level": "C",
        "description": "ARID1A is a SWI/SNF chromatin remodeling subunit mutated in ~50% of ovarian clear "
                       "cell carcinoma. Loss may confer sensitivity to EZH2 inhibitors and ATR inhibitors. "
                       "Emerging therapeutic target.",
    },
    "MSI_H": {
        "gene": "MSI-H/dMMR",
        "full_name": "Microsatellite Instability High / Mismatch Repair Deficient",
        "cancer_types": ["tissue-agnostic", "colorectal", "endometrial", "gastric"],
        "key_variants": ["MLH1 promoter methylation", "MLH1/MSH2/MSH6/PMS2 loss"],
        "targeted_therapies": ["pembrolizumab", "nivolumab", "dostarlimab"],
        "combination_therapies": ["nivolumab + ipilimumab"],
        "resistance_mutations": ["B2M loss", "JAK1/2 mutations", "STK11 co-mutation"],
        "pathway": "DDR/immune",
        "evidence_level": "A",
        "description": "MSI-H/dMMR is a tissue-agnostic biomarker for immunotherapy (pembrolizumab, "
                       "dostarlimab). KEYNOTE-177 established first-line pembrolizumab for MSI-H CRC. "
                       "CheckMate-8HW showed nivo+ipi superiority.",
    },
    "NRG1": {
        "gene": "NRG1",
        "full_name": "Neuregulin 1 Fusions",
        "cancer_types": ["NSCLC", "pancreatic", "breast", "cholangiocarcinoma"],
        "key_variants": ["CD74-NRG1 fusion", "SLC3A2-NRG1 fusion", "ATP1B1-NRG1 fusion"],
        "targeted_therapies": ["afatinib", "zenocutuzumab"],
        "combination_therapies": ["afatinib + cetuximab"],
        "resistance_mutations": ["secondary NRG1 fusions", "HER3 bypass"],
        "pathway": "HER2/HER3",
        "evidence_level": "B",
        "description": "NRG1 fusions drive HER2-HER3 heterodimer signaling. Rare (~0.2% pan-cancer) "
                       "but highly actionable. Zenocutuzumab (anti-HER2/HER3 bispecific) showed 35% "
                       "ORR in eNRGy trial.",
    },
    "KRASG12D": {
        "gene": "KRAS G12D",
        "full_name": "KRAS G12D Allele-Specific",
        "cancer_types": ["pancreatic", "CRC", "NSCLC"],
        "key_variants": ["G12D"],
        "targeted_therapies": ["MRTX1133 (investigational)", "divarasib (cross-allele, investigational)"],
        "combination_therapies": [],
        "resistance_mutations": ["adaptive MAPK reactivation", "RTK bypass"],
        "pathway": "MAPK",
        "evidence_level": "D",
        "description": "KRAS G12D is the most common KRAS mutation in pancreatic cancer (~35%). "
                       "MRTX1133 is a first-in-class selective G12D inhibitor in Phase 1/2.",
    },
    "TROP2": {
        "gene": "TROP2",
        "full_name": "Trophoblast Cell Surface Antigen 2",
        "cancer_types": ["TNBC", "urothelial", "NSCLC", "endometrial", "HR+/HER2- breast"],
        "key_variants": ["TROP2 overexpression"],
        "targeted_therapies": ["sacituzumab govitecan", "datopotamab deruxtecan"],
        "combination_therapies": ["sacituzumab + pembrolizumab"],
        "resistance_mutations": ["TROP2 downregulation", "TOP1 mutations", "drug efflux (ABCG2)"],
        "pathway": "signaling/ADC target",
        "evidence_level": "A",
        "description": "TROP2 is a transmembrane glycoprotein overexpressed in >80% of TNBC. "
                       "Sacituzumab govitecan FDA-approved for mTNBC, HR+/HER2- breast, and mUC.",
    },
    "NECTIN4": {
        "gene": "NECTIN4",
        "full_name": "Nectin Cell Adhesion Molecule 4",
        "cancer_types": ["urothelial", "breast", "NSCLC", "gastric"],
        "key_variants": ["Nectin-4 overexpression"],
        "targeted_therapies": ["enfortumab vedotin"],
        "combination_therapies": ["enfortumab vedotin + pembrolizumab"],
        "resistance_mutations": ["Nectin-4 downregulation", "drug efflux", "EMT"],
        "pathway": "cell adhesion/ADC target",
        "evidence_level": "A",
        "description": "Nectin-4 is expressed in >97% of urothelial carcinoma. Enfortumab vedotin + "
                       "pembrolizumab is first-line for cisplatin-ineligible mUC (EV-302/KEYNOTE-A39).",
    },
    "FOLR1": {
        "gene": "FOLR1",
        "full_name": "Folate Receptor Alpha",
        "cancer_types": ["ovarian (high-grade serous)", "NSCLC", "endometrial"],
        "key_variants": ["FRα overexpression (IHC PS2+, ≥75% cells)"],
        "targeted_therapies": ["mirvetuximab soravtansine"],
        "combination_therapies": ["mirvetuximab + bevacizumab (investigating)"],
        "resistance_mutations": ["FRα loss", "tubulin mutations"],
        "pathway": "folate transport/ADC target",
        "evidence_level": "A",
        "description": "FRα is overexpressed in ~40% of platinum-resistant ovarian cancer. Mirvetuximab "
                       "soravtansine FDA-approved for FRα-high, platinum-resistant ovarian cancer "
                       "(SORAYA, MIRASOL trials).",
    },
    "CLDN18_2": {
        "gene": "CLDN18.2",
        "full_name": "Claudin 18.2",
        "cancer_types": ["gastric", "gastroesophageal junction", "pancreatic"],
        "key_variants": ["CLDN18.2 expression (IHC ≥2+ in ≥75% cells)"],
        "targeted_therapies": ["zolbetuximab"],
        "combination_therapies": ["zolbetuximab + CAPOX", "zolbetuximab + mFOLFOX6"],
        "resistance_mutations": ["CLDN18.2 downregulation", "isoform switching (CLDN18.1)"],
        "pathway": "tight junction/targeted",
        "evidence_level": "A",
        "description": "CLDN18.2 is a tight junction protein expressed in ~38% of gastric cancers. "
                       "Zolbetuximab + chemo improved OS in SPOTLIGHT and GLOW trials. First targeted "
                       "therapy for CLDN18.2+ gastric cancer.",
    },
    "DLL3": {
        "gene": "DLL3",
        "full_name": "Delta-Like Ligand 3",
        "cancer_types": ["SCLC", "large cell neuroendocrine", "neuroendocrine tumors"],
        "key_variants": ["DLL3 surface expression"],
        "targeted_therapies": ["tarlatamab", "rovalpituzumab tesirine (withdrawn)"],
        "combination_therapies": ["tarlatamab monotherapy or + chemo"],
        "resistance_mutations": ["DLL3 loss", "immune evasion", "neuroendocrine-to-non-NE transition"],
        "pathway": "Notch (inhibitory ligand)",
        "evidence_level": "B",
        "description": "DLL3 is expressed on >80% of SCLC tumors but absent from normal tissue. "
                       "Tarlatamab (DLL3xCD3 BiTE) showed 40% ORR in DeLLphi-301 trial for relapsed "
                       "SCLC. FDA accelerated approval 2024.",
    },
    "BARD1": {
        "gene": "BARD1",
        "full_name": "BRCA1-Associated RING Domain 1",
        "cancer_types": ["ovarian", "breast"],
        "key_variants": ["truncating mutations", "splice site variants"],
        "targeted_therapies": ["olaparib (by HRD inference)", "rucaparib"],
        "combination_therapies": ["PARP + platinum"],
        "resistance_mutations": ["HR restoration", "reversion mutations"],
        "pathway": "DDR",
        "evidence_level": "C",
        "description": "BARD1 forms obligate heterodimer with BRCA1 for HR repair. BARD1 pathogenic "
                       "variants confer moderate-to-high cancer risk and potential PARPi sensitivity.",
    },
    "RAD51C": {
        "gene": "RAD51C",
        "full_name": "RAD51 Paralog C",
        "cancer_types": ["ovarian", "breast"],
        "key_variants": ["truncating mutations", "large deletions"],
        "targeted_therapies": ["olaparib", "rucaparib", "niraparib"],
        "combination_therapies": ["PARP + platinum"],
        "resistance_mutations": ["reversion mutations", "RAD51 foci restoration"],
        "pathway": "DDR",
        "evidence_level": "B",
        "description": "RAD51C mutations confer HRD and PARPi sensitivity comparable to BRCA2. "
                       "Included in most HRD gene panels. ~1% of ovarian cancers.",
    },
    "CDH1": {
        "gene": "CDH1",
        "full_name": "E-Cadherin",
        "cancer_types": ["lobular breast", "hereditary diffuse gastric"],
        "key_variants": ["germline truncating mutations", "splice site", "promoter methylation (somatic)"],
        "targeted_therapies": ["CDK4/6 inhibitors (lobular breast responds well)"],
        "combination_therapies": ["CDK4/6i + endocrine therapy"],
        "resistance_mutations": ["EMT completion", "CDK6 amplification"],
        "pathway": "cell adhesion/WNT",
        "evidence_level": "B",
        "description": "Germline CDH1 mutations cause hereditary diffuse gastric cancer syndrome (HDGC) "
                       "with >80% lifetime gastric cancer risk. Somatic CDH1 loss defines lobular breast "
                       "cancer. Prophylactic gastrectomy recommended for carriers.",
    },
    "MAP2K1": {
        "gene": "MAP2K1",
        "full_name": "MEK1 (Mitogen-Activated Protein Kinase Kinase 1)",
        "cancer_types": ["melanoma (BRAF-WT)", "NSCLC", "histiocytic neoplasms", "low-grade serous ovarian"],
        "key_variants": ["K57N", "C121S", "Q56P", "E203K"],
        "targeted_therapies": ["trametinib", "cobimetinib", "binimetinib", "selumetinib"],
        "combination_therapies": ["trametinib + dabrafenib (for BRAF-mutant context)",
                                   "selumetinib monotherapy (NF1 plexiform neurofibromas)"],
        "resistance_mutations": ["MEK2 mutations", "ERK amplification", "PI3K bypass"],
        "pathway": "MAPK",
        "evidence_level": "A",
        "description": "MAP2K1 (MEK1) mutations activate MAPK signaling. MEK inhibitors are approved in "
                       "combination with BRAF inhibitors for BRAF-mutant cancers and as monotherapy for "
                       "NF1 plexiform neurofibromas (selumetinib).",
    },
    "SMARCB1": {
        "gene": "SMARCB1",
        "full_name": "SMARCB1/INI1 (SWI/SNF Related Matrix Associated Actin Dependent Regulator)",
        "cancer_types": ["rhabdoid tumor", "epithelioid sarcoma", "meningioma"],
        "key_variants": ["biallelic inactivation", "homozygous deletion", "truncating mutations"],
        "targeted_therapies": ["tazemetostat"],
        "combination_therapies": [],
        "resistance_mutations": ["SWI/SNF complex rewiring"],
        "pathway": "chromatin_remodeling",
        "evidence_level": "A",
        "description": "SMARCB1 (INI1) loss is pathognomonic for malignant rhabdoid tumors and epithelioid "
                       "sarcomas. Tazemetostat (EZH2 inhibitor) is FDA-approved for INI1-negative "
                       "epithelioid sarcoma.",
    },
    "EZH2": {
        "gene": "EZH2",
        "full_name": "Enhancer of Zeste Homolog 2",
        "cancer_types": ["follicular lymphoma", "epithelioid sarcoma", "DLBCL"],
        "key_variants": ["Y646N", "Y646F", "Y646H", "A677G", "A687V (gain-of-function)"],
        "targeted_therapies": ["tazemetostat"],
        "combination_therapies": ["tazemetostat + R-CHOP (investigating)"],
        "resistance_mutations": ["SWI/SNF mutations", "secondary EZH2 mutations"],
        "pathway": "chromatin_remodeling",
        "evidence_level": "A",
        "description": "EZH2 gain-of-function mutations occur in ~20% of follicular lymphoma. "
                       "Tazemetostat FDA-approved for EZH2-mutant FL and epithelioid sarcoma. EZH2 is "
                       "the catalytic subunit of PRC2.",
    },
    "CTLA4": {
        "gene": "CTLA-4",
        "full_name": "Cytotoxic T-Lymphocyte Associated Protein 4",
        "cancer_types": ["melanoma", "RCC", "MSI-H (tissue-agnostic)", "HCC", "NSCLC"],
        "key_variants": ["N/A (drug target, not mutated)"],
        "targeted_therapies": ["ipilimumab", "tremelimumab"],
        "combination_therapies": ["ipilimumab + nivolumab", "tremelimumab + durvalumab"],
        "resistance_mutations": ["Treg expansion", "compensatory immune checkpoints (TIGIT, LAG-3)", "cold TME"],
        "pathway": "immune checkpoint",
        "evidence_level": "A",
        "description": "CTLA-4 blockade was the first FDA-approved checkpoint inhibitor (ipilimumab, 2011). "
                       "Ipilimumab + nivolumab dual checkpoint blockade is standard in melanoma, RCC, "
                       "and MSI-H CRC.",
    },
    # ── Pediatric Oncology Targets ──────────────────────────────────────
    "ETV6_RUNX1": {
        "gene": "ETV6-RUNX1",
        "full_name": "ETV6-RUNX1 Fusion (t(12;21)) — Pediatric B-ALL",
        "cancer_types": ["pediatric B-ALL"],
        "key_variants": ["ETV6-RUNX1 fusion", "t(12;21)(p13;q22)"],
        "targeted_therapies": ["vincristine", "daunorubicin", "methotrexate", "L-asparaginase",
                                "6-mercaptopurine"],
        "combination_therapies": ["COG AALL0434 regimen", "BFM-based chemotherapy"],
        "resistance_mutations": ["CREBBP mutations", "NR3C1 (glucocorticoid resistance)",
                                  "NT5C2 mutations (thiopurine resistance)"],
        "pathway": "transcription_regulation",
        "evidence_level": "A",
        "description": "ETV6-RUNX1 is the most common fusion in pediatric B-ALL (~25%). "
                       "Associated with favorable prognosis (>90% EFS). Standard treatment "
                       "uses risk-adapted chemotherapy. CD19/CD22 are therapeutic targets "
                       "for relapsed disease.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "1-10 years (peak incidence)",
            "frequency": "~25% of pediatric B-ALL",
            "prognosis": "favorable (>90% 5-year EFS)",
            "targets": ["CD19", "CD22", "CRLF2"],
        },
    },
    "BCR_ABL1": {
        "gene": "BCR-ABL1",
        "full_name": "BCR-ABL1 Fusion (Philadelphia Chromosome) — Pediatric ALL",
        "cancer_types": ["pediatric B-ALL", "CML", "Ph+ ALL"],
        "key_variants": ["BCR-ABL1 p190 fusion", "BCR-ABL1 p210 fusion", "t(9;22)(q34;q11)"],
        "targeted_therapies": ["imatinib", "dasatinib", "ponatinib"],
        "combination_therapies": ["dasatinib + chemotherapy (COG AALL0434)",
                                   "imatinib + hyper-CVAD"],
        "resistance_mutations": ["T315I (gatekeeper)", "E255K", "Y253H", "compound mutations"],
        "pathway": "ABL kinase/MAPK/STAT5",
        "evidence_level": "A",
        "description": "BCR-ABL1 (Ph+) occurs in ~3-5% of pediatric ALL but ~25% of adult ALL. "
                       "Addition of TKIs (dasatinib/imatinib) to chemotherapy has improved "
                       "outcomes dramatically. T315I is the most common resistance mutation.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "all pediatric ages (increases with age)",
            "frequency": "~3-5% of pediatric ALL",
            "prognosis": "historically poor, improved with TKI + chemo (>70% EFS)",
            "targets": ["CD19", "CD22"],
        },
    },
    "NOTCH1": {
        "gene": "NOTCH1",
        "full_name": "Notch Receptor 1 — Pediatric T-ALL",
        "cancer_types": ["pediatric T-ALL", "T-ALL", "CLL"],
        "key_variants": ["activating mutations (HD domain)", "PEST domain truncations",
                          "juxtamembrane expansion mutations"],
        "targeted_therapies": [],
        "combination_therapies": ["COG T-ALL chemotherapy regimens",
                                   "gamma-secretase inhibitors (investigational)"],
        "resistance_mutations": [],
        "pathway": "Notch",
        "evidence_level": "B",
        "description": "NOTCH1 activating mutations occur in >60% of T-ALL cases. "
                       "Paradoxically, NOTCH1 mutations are associated with favorable prognosis "
                       "in pediatric T-ALL. CDKN2A deletions are frequent co-occurring events. "
                       "Gamma-secretase inhibitors are under investigation.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "childhood/adolescent",
            "frequency": ">60% of T-ALL",
            "prognosis": "favorable (NOTCH1-mutant)",
            "co_alterations": ["CDKN2A deletions", "FBXW7 mutations"],
        },
    },
    "MYCN": {
        "gene": "MYCN",
        "full_name": "MYCN Proto-Oncogene (N-Myc) — Neuroblastoma",
        "cancer_types": ["neuroblastoma", "retinoblastoma", "medulloblastoma",
                          "rhabdomyosarcoma", "Wilms tumor"],
        "key_variants": ["amplification (>10 copies)", "gain"],
        "targeted_therapies": ["dinutuximab (anti-GD2)", "crizotinib (ALK co-mutation)",
                                "isotretinoin (differentiation)"],
        "combination_therapies": ["dinutuximab + GM-CSF + IL-2 + isotretinoin",
                                   "high-dose chemo + autologous SCT"],
        "resistance_mutations": ["ALK mutations", "RAS-MAPK activation", "telomere maintenance"],
        "pathway": "MYC transcription",
        "evidence_level": "A",
        "description": "MYCN amplification defines high-risk neuroblastoma (~20% of cases) and "
                       "is the single most important prognostic factor. MYCN-amplified tumors "
                       "require intensive multimodal therapy including high-dose chemotherapy, "
                       "surgery, radiation, autologous SCT, and anti-GD2 immunotherapy.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "median diagnosis 17 months",
            "frequency": "~20% of neuroblastoma",
            "prognosis": "poor (MYCN amp = high risk, <50% 5yr OS)",
            "risk_stratification": "MYCN amplification = automatic high-risk classification",
            "co_alterations": ["1p deletion", "11q deletion", "ALK mutations"],
            "targets": ["GD2 (dinutuximab)", "ALK (crizotinib)"],
        },
    },
    "PTCH1": {
        "gene": "PTCH1",
        "full_name": "Patched 1 — Medulloblastoma SHH Subtype",
        "cancer_types": ["medulloblastoma (SHH)", "basal cell carcinoma", "Gorlin syndrome"],
        "key_variants": ["loss-of-function mutations", "deletions", "splice site mutations"],
        "targeted_therapies": ["vismodegib", "sonidegib"],
        "combination_therapies": [],
        "resistance_mutations": ["SMO mutations (D473H)", "GLI2 amplification", "PI3K activation"],
        "pathway": "Hedgehog (SHH)",
        "evidence_level": "B",
        "description": "PTCH1 loss-of-function activates the SHH pathway constitutively. In "
                       "pediatric medulloblastoma, the SHH subtype (~30%) is driven by PTCH1, "
                       "SMO, or SUFU mutations. Vismodegib shows activity in adult SHH-MB but "
                       "is contraindicated in skeletally immature children due to growth plate "
                       "effects. Medulloblastoma subtypes: WNT (best prognosis), SHH, Group 3 "
                       "(worst prognosis), Group 4.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "bimodal: infants (<3yr) and adolescents",
            "frequency": "~30% of medulloblastoma",
            "subtypes": ["WNT", "SHH", "Group 3", "Group 4"],
            "shh_genes": ["PTCH1", "SMO", "SUFU", "GLI1", "GLI2"],
            "caution": "Vismodegib contraindicated in skeletally immature patients",
        },
    },
    "PAX3_FOXO1": {
        "gene": "PAX3-FOXO1",
        "full_name": "PAX3-FOXO1 Fusion — Alveolar Rhabdomyosarcoma",
        "cancer_types": ["alveolar rhabdomyosarcoma"],
        "key_variants": ["PAX3-FOXO1 fusion", "t(2;13)(q35;q14)", "PAX7-FOXO1 (less common)"],
        "targeted_therapies": [],
        "combination_therapies": ["VAC (vincristine + actinomycin D + cyclophosphamide)",
                                   "VI (vincristine + irinotecan)"],
        "resistance_mutations": [],
        "pathway": "transcription_regulation",
        "evidence_level": "A",
        "description": "PAX3-FOXO1 fusion defines alveolar rhabdomyosarcoma and confers poor "
                       "prognosis compared to fusion-negative (embryonal) RMS. Embryonal RMS "
                       "is driven by RAS pathway mutations. Treatment is multimodal: surgery, "
                       "radiation, and combination chemotherapy (VAC backbone).",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "adolescents (alveolar); young children (embryonal)",
            "frequency": "~20% of pediatric soft tissue sarcomas",
            "embryonal_drivers": ["RAS pathway (NRAS, KRAS, HRAS)", "FGFR4", "PIK3CA"],
            "prognosis": "fusion-positive: poor; fusion-negative: favorable",
        },
    },
    "EWSR1_FLI1": {
        "gene": "EWSR1-FLI1",
        "full_name": "EWSR1-FLI1 Fusion — Ewing Sarcoma",
        "cancer_types": ["Ewing sarcoma", "primitive neuroectodermal tumor"],
        "key_variants": ["EWSR1-FLI1 fusion", "t(11;22)(q24;q12)", "EWSR1-ERG (variant)"],
        "targeted_therapies": [],
        "combination_therapies": ["VDC/IE (vincristine + doxorubicin + cyclophosphamide "
                                   "alternating with ifosfamide + etoposide)"],
        "resistance_mutations": ["STAG2 loss", "TP53 mutations", "CDKN2A deletion"],
        "pathway": "transcription_regulation",
        "evidence_level": "A",
        "description": "EWSR1-FLI1 fusion is pathognomonic for Ewing sarcoma (~85% of cases). "
                       "Diagnosis requires FISH or RT-PCR confirmation. Standard treatment is "
                       "compressed VDC/IE chemotherapy with local control (surgery and/or "
                       "radiation). 5-year OS ~70% for localized disease.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "10-20 years (peak incidence)",
            "frequency": "~200 cases/year in US",
            "prognosis": "localized ~70% 5yr OS; metastatic ~30% 5yr OS",
        },
    },
    "WT1": {
        "gene": "WT1",
        "full_name": "Wilms Tumor 1 — Nephroblastoma",
        "cancer_types": ["Wilms tumor (nephroblastoma)", "AML"],
        "key_variants": ["WT1 truncating mutations", "WT1 deletions", "11p13 deletion"],
        "targeted_therapies": [],
        "combination_therapies": ["vincristine + actinomycin D (favorable histology)",
                                   "vincristine + actinomycin D + doxorubicin (higher stage)"],
        "resistance_mutations": ["TP53 mutations (anaplastic histology)"],
        "pathway": "transcription_regulation",
        "evidence_level": "B",
        "description": "WT1 mutations occur in ~10-20% of Wilms tumors. Other drivers include "
                       "WTX (AMER1) and CTNNB1 (beta-catenin). Wilms tumor has >90% overall "
                       "survival with combined surgery and chemotherapy. Anaplastic histology "
                       "(TP53 mutations) confers worse prognosis.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "median diagnosis 3-4 years",
            "frequency": "~500 cases/year in US",
            "co_alterations": ["WTX (AMER1)", "CTNNB1", "TP53 (anaplastic)"],
            "syndromes": ["WAGR (WT1 deletion)", "Denys-Drash", "Beckwith-Wiedemann"],
            "prognosis": "favorable histology >90% OS; anaplastic ~50% OS",
        },
    },
    "CTNNB1_PED": {
        "gene": "CTNNB1",
        "full_name": "Beta-Catenin — Hepatoblastoma & Pediatric Tumors",
        "cancer_types": ["hepatoblastoma", "Wilms tumor", "medulloblastoma (WNT subtype)",
                          "desmoid tumors", "pilomatricoma"],
        "key_variants": ["exon 3 activating mutations", "S45del", "T41A", "S33C"],
        "targeted_therapies": [],
        "combination_therapies": ["cisplatin-based therapy (hepatoblastoma: PLADO regimen)",
                                   "cisplatin + doxorubicin (C5V)"],
        "resistance_mutations": [],
        "pathway": "WNT",
        "evidence_level": "B",
        "description": "CTNNB1 activating mutations drive hepatoblastoma (~80% of cases) and "
                       "WNT-subtype medulloblastoma. Hepatoblastoma is treated with cisplatin-based "
                       "chemotherapy and surgical resection. Liver transplant for unresectable disease. "
                       "CTNNB1-mutant WNT medulloblastoma has the best prognosis of all MB subtypes.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "hepatoblastoma: <3 years; WNT-MB: older children",
            "frequency": "hepatoblastoma ~80%; WNT-MB ~10% of all MB",
            "prognosis": "hepatoblastoma: >80% OS; WNT-MB: >95% OS",
        },
    },
    "RB1": {
        "gene": "RB1",
        "full_name": "Retinoblastoma 1 — Retinoblastoma & Cell Cycle Regulation",
        "cancer_types": ["retinoblastoma", "osteosarcoma", "SCLC", "bladder"],
        "key_variants": ["biallelic loss", "truncating mutations", "deletions",
                          "promoter methylation"],
        "targeted_therapies": [],
        "combination_therapies": ["intra-arterial melphalan (retinoblastoma)",
                                   "intravitreal melphalan", "systemic chemoreduction (VEC)"],
        "resistance_mutations": [],
        "pathway": "cell_cycle",
        "evidence_level": "B",
        "description": "RB1 biallelic loss is the hallmark of retinoblastoma (Knudson two-hit "
                       "hypothesis). Hereditary retinoblastoma (~40%) involves germline RB1 "
                       "mutation + somatic second hit. Treatment includes focal therapy, "
                       "chemoreduction, intra-arterial chemotherapy, and enucleation. "
                       "Survivors have increased risk of secondary malignancies (osteosarcoma).",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "median diagnosis 18 months (bilateral), 24 months (unilateral)",
            "frequency": "~300 cases/year in US",
            "hereditary": "~40% germline RB1 mutation (bilateral, earlier onset)",
            "secondary_cancer_risk": "osteosarcoma, soft tissue sarcoma, melanoma",
            "prognosis": ">95% OS in developed countries",
        },
    },
    "KMT2A": {
        "gene": "KMT2A",
        "full_name": "KMT2A (MLL) Rearrangements — Pediatric AML/ALL",
        "cancer_types": ["pediatric AML", "infant ALL", "therapy-related AML"],
        "key_variants": ["KMT2A-MLLT3 (t(9;11))", "KMT2A-AFF1 (t(4;11))", "KMT2A-MLLT1 (t(11;19))",
                          "KMT2A-MLLT10 (t(10;11))", "NUP98 fusions (NUP98-NSD1, NUP98-KDM5A)"],
        "targeted_therapies": ["menin inhibitors (revumenib, ziftomenib — investigational)"],
        "combination_therapies": ["intensive AML chemotherapy (ADE: cytarabine + daunorubicin + etoposide)"],
        "resistance_mutations": ["menin inhibitor resistance mutations"],
        "pathway": "chromatin_remodeling/transcription",
        "evidence_level": "B",
        "description": "KMT2A (MLL) rearrangements occur in ~20% of pediatric AML and ~80% of "
                       "infant ALL. KMT2A-AFF1 in infant ALL confers very poor prognosis. "
                       "NUP98 fusions are another high-risk pediatric AML subtype. Menin "
                       "inhibitors (revumenib) are emerging as targeted therapy for KMT2A-r "
                       "and NUP98-r leukemias. FLT3-ITD co-occurs in ~15% of pediatric AML "
                       "and is targetable with midostaurin/gilteritinib.",
        "pediatric": True,
        "pediatric_context": {
            "age_group": "infants (ALL); all ages (AML)",
            "frequency": "~20% of pediatric AML; ~80% of infant ALL",
            "prognosis": "KMT2A-r infant ALL: <40% EFS; KMT2A-r AML: ~50% EFS",
            "co_alterations": ["FLT3-ITD (~15% pediatric AML)", "NUP98 fusions"],
            "emerging_therapy": "menin inhibitors (revumenib, ziftomenib)",
        },
    },
}


# ═══════════════════════════════════════════════════════════════════════
# 2. THERAPY MAP (~30 drugs)
# ═══════════════════════════════════════════════════════════════════════

THERAPY_MAP: Dict[str, Dict] = {
    "osimertinib": {
        "drug_name": "osimertinib",
        "brand_name": "Tagrisso",
        "category": "targeted",
        "targets": ["EGFR"],
        "approved_indications": ["EGFR-mutant NSCLC (1L)", "EGFR T790M+ NSCLC (2L+)", "adjuvant EGFR+ NSCLC"],
        "mechanism": "3rd-generation EGFR TKI, irreversible binding to C797, mutant-selective",
        "key_trials": ["FLAURA", "FLAURA2", "AURA3", "ADAURA"],
        "resistance_mechanisms": ["C797S", "MET amplification", "small cell transformation", "HER2 amp"],
        "evidence_level": "A",
    },
    "pembrolizumab": {
        "drug_name": "pembrolizumab",
        "brand_name": "Keytruda",
        "category": "immunotherapy",
        "targets": ["PD-1"],
        "approved_indications": ["NSCLC (PD-L1≥50% 1L mono, PD-L1≥1% + chemo)", "melanoma",
                                  "MSI-H/dMMR (tissue-agnostic)", "TMB-H (≥10 mut/Mb)",
                                  "head and neck", "cervical", "gastric", "TNBC", "RCC"],
        "mechanism": "Anti-PD-1 monoclonal antibody, blocks PD-1/PD-L1 interaction to restore T-cell activity",
        "key_trials": ["KEYNOTE-024", "KEYNOTE-189", "KEYNOTE-158", "KEYNOTE-177", "KEYNOTE-826"],
        "resistance_mechanisms": ["B2M loss", "JAK1/2 loss", "STK11 co-mutation", "WNT/β-catenin activation"],
        "evidence_level": "A",
    },
    "nivolumab": {
        "drug_name": "nivolumab",
        "brand_name": "Opdivo",
        "category": "immunotherapy",
        "targets": ["PD-1"],
        "approved_indications": ["melanoma", "NSCLC", "RCC", "Hodgkin lymphoma", "head and neck",
                                  "urothelial", "MSI-H CRC", "HCC", "gastric"],
        "mechanism": "Anti-PD-1 monoclonal antibody",
        "key_trials": ["CheckMate-227", "CheckMate-9LA", "CheckMate-274", "CheckMate-8HW"],
        "resistance_mechanisms": ["B2M loss", "JAK1/2 mutations", "PTEN loss"],
        "evidence_level": "A",
    },
    "vemurafenib": {
        "drug_name": "vemurafenib",
        "brand_name": "Zelboraf",
        "category": "targeted",
        "targets": ["BRAF V600"],
        "approved_indications": ["BRAF V600E melanoma", "BRAF V600E Erdheim-Chester disease"],
        "mechanism": "Selective BRAF V600 kinase inhibitor",
        "key_trials": ["BRIM-3"],
        "resistance_mechanisms": ["MEK1/2 mutations", "NRAS activation", "BRAF amplification", "COT overexpression"],
        "evidence_level": "A",
    },
    "dabrafenib": {
        "drug_name": "dabrafenib",
        "brand_name": "Tafinlar",
        "category": "targeted",
        "targets": ["BRAF V600"],
        "approved_indications": ["BRAF V600E/K melanoma", "BRAF V600E NSCLC", "BRAF V600E ATC"],
        "mechanism": "Selective BRAF kinase inhibitor, typically combined with trametinib (MEK inhibitor)",
        "key_trials": ["COMBI-d", "COMBI-v", "BRF113928"],
        "resistance_mechanisms": ["MEK mutations", "RAS activation", "BRAF amplification"],
        "evidence_level": "A",
    },
    "sotorasib": {
        "drug_name": "sotorasib",
        "brand_name": "Lumakras",
        "category": "targeted",
        "targets": ["KRAS G12C"],
        "approved_indications": ["KRAS G12C NSCLC (2L+)"],
        "mechanism": "Irreversible, selective KRAS G12C inhibitor — locks KRAS in inactive GDP-bound state",
        "key_trials": ["CodeBreaK 100", "CodeBreaK 200"],
        "resistance_mechanisms": ["Y96D", "R68S", "H95 mutations", "MET amplification",
                                   "secondary KRAS mutations", "RAS-MAPK reactivation"],
        "evidence_level": "A",
    },
    "adagrasib": {
        "drug_name": "adagrasib",
        "brand_name": "Krazati",
        "category": "targeted",
        "targets": ["KRAS G12C"],
        "approved_indications": ["KRAS G12C NSCLC (2L+)"],
        "mechanism": "Irreversible KRAS G12C inhibitor with long half-life and CNS penetration",
        "key_trials": ["KRYSTAL-1", "KRYSTAL-7"],
        "resistance_mechanisms": ["secondary KRAS mutations", "MET amplification", "MAPK reactivation"],
        "evidence_level": "A",
    },
    "larotrectinib": {
        "drug_name": "larotrectinib",
        "brand_name": "Vitrakvi",
        "category": "targeted",
        "targets": ["NTRK1/2/3"],
        "approved_indications": ["NTRK fusion-positive solid tumors (tissue-agnostic)"],
        "mechanism": "Selective pan-TRK inhibitor",
        "key_trials": ["NAVIGATE", "SCOUT", "LOXO-TRK-14001"],
        "resistance_mechanisms": ["solvent front mutations (G595R, G623R)", "xDFG mutations"],
        "evidence_level": "A",
    },
    "entrectinib": {
        "drug_name": "entrectinib",
        "brand_name": "Rozlytrek",
        "category": "targeted",
        "targets": ["NTRK1/2/3", "ROS1", "ALK"],
        "approved_indications": ["NTRK fusion-positive solid tumors", "ROS1+ NSCLC"],
        "mechanism": "Pan-TRK, ROS1, and ALK inhibitor with CNS penetration",
        "key_trials": ["STARTRK-2", "STARTRK-1"],
        "resistance_mechanisms": ["solvent front mutations", "xDFG mutations", "G2032R (ROS1)"],
        "evidence_level": "A",
    },
    "selpercatinib": {
        "drug_name": "selpercatinib",
        "brand_name": "Retevmo",
        "category": "targeted",
        "targets": ["RET"],
        "approved_indications": ["RET fusion-positive NSCLC", "RET-mutant medullary thyroid",
                                  "RET fusion-positive thyroid (age≥12)"],
        "mechanism": "Highly selective RET kinase inhibitor designed to minimize off-target VEGFR activity",
        "key_trials": ["LIBRETTO-001", "LIBRETTO-431"],
        "resistance_mechanisms": ["G810R/S/C", "V804M", "Y806C/H"],
        "evidence_level": "A",
    },
    "capmatinib": {
        "drug_name": "capmatinib",
        "brand_name": "Tabrecta",
        "category": "targeted",
        "targets": ["MET"],
        "approved_indications": ["MET exon 14 skipping NSCLC"],
        "mechanism": "Selective MET kinase inhibitor",
        "key_trials": ["GEOMETRY mono-1"],
        "resistance_mechanisms": ["D1228N/V", "Y1230H/C", "KRAS amplification"],
        "evidence_level": "A",
    },
    "erdafitinib": {
        "drug_name": "erdafitinib",
        "brand_name": "Balversa",
        "category": "targeted",
        "targets": ["FGFR1/2/3/4"],
        "approved_indications": ["FGFR-altered urothelial carcinoma"],
        "mechanism": "Pan-FGFR tyrosine kinase inhibitor",
        "key_trials": ["THOR"],
        "resistance_mechanisms": ["gatekeeper mutations (V565I)", "FGFR2 secondary mutations"],
        "evidence_level": "A",
    },
    "alpelisib": {
        "drug_name": "alpelisib",
        "brand_name": "Piqray",
        "category": "targeted",
        "targets": ["PI3K alpha (PIK3CA)"],
        "approved_indications": ["PIK3CA-mutant HR+/HER2- advanced breast cancer (with fulvestrant)"],
        "mechanism": "PI3Kα-selective inhibitor",
        "key_trials": ["SOLAR-1"],
        "resistance_mechanisms": ["PTEN loss", "AKT1 E17K", "mTOR mutations"],
        "evidence_level": "A",
    },
    "olaparib": {
        "drug_name": "olaparib",
        "brand_name": "Lynparza",
        "category": "targeted",
        "targets": ["PARP1/2"],
        "approved_indications": ["BRCA-mutant breast cancer", "BRCA-mutant ovarian cancer",
                                  "BRCA-mutant prostate cancer", "BRCA-mutant pancreatic cancer",
                                  "HRD+ ovarian cancer (with bevacizumab)"],
        "mechanism": "PARP1/2 inhibitor — traps PARP on DNA, synthetic lethality with HRD",
        "key_trials": ["OlympiAD", "SOLO-1", "SOLO-2", "PROfound", "POLO", "PAOLA-1"],
        "resistance_mechanisms": ["BRCA reversion mutations", "53BP1 loss", "RAD51 upregulation",
                                   "drug efflux (ABCB1)"],
        "evidence_level": "A",
    },
    "ivosidenib": {
        "drug_name": "ivosidenib",
        "brand_name": "Tibsovo",
        "category": "targeted",
        "targets": ["IDH1"],
        "approved_indications": ["IDH1-mutant AML (newly diagnosed, R/R)",
                                  "IDH1-mutant cholangiocarcinoma"],
        "mechanism": "IDH1 inhibitor — blocks mutant IDH1 enzyme, reduces 2-HG oncometabolite",
        "key_trials": ["ClarIDHy", "AGILE"],
        "resistance_mechanisms": ["second-site IDH1 mutations", "RTK pathway activation"],
        "evidence_level": "A",
    },
    "trastuzumab_deruxtecan": {
        "drug_name": "trastuzumab deruxtecan",
        "brand_name": "Enhertu",
        "category": "targeted",
        "targets": ["HER2"],
        "approved_indications": ["HER2+ metastatic breast cancer (2L+)", "HER2-low breast cancer",
                                  "HER2+ gastric cancer", "HER2-mutant NSCLC"],
        "mechanism": "HER2-directed antibody-drug conjugate (ADC) with topoisomerase I inhibitor payload; "
                     "bystander effect enables activity in HER2-low tumors",
        "key_trials": ["DESTINY-Breast03", "DESTINY-Breast04", "DESTINY-Lung01", "DESTINY-Gastric01"],
        "resistance_mechanisms": ["HER2 downregulation", "drug efflux", "SLX4-dependent DNA repair"],
        "evidence_level": "A",
    },
    "palbociclib": {
        "drug_name": "palbociclib",
        "brand_name": "Ibrance",
        "category": "targeted",
        "targets": ["CDK4/6"],
        "approved_indications": ["HR+/HER2- advanced breast cancer (with endocrine therapy)"],
        "mechanism": "Selective CDK4/6 inhibitor — blocks RB phosphorylation, induces G1 arrest",
        "key_trials": ["PALOMA-2", "PALOMA-3"],
        "resistance_mechanisms": ["RB1 loss", "CDK6 amplification", "CCNE1 amplification",
                                   "FGFR pathway activation"],
        "evidence_level": "A",
    },
    "lorlatinib": {
        "drug_name": "lorlatinib",
        "brand_name": "Lorbrena",
        "category": "targeted",
        "targets": ["ALK", "ROS1"],
        "approved_indications": ["ALK+ NSCLC (1L and 2L+)"],
        "mechanism": "3rd-generation ALK TKI — macrocyclic structure with CNS penetration, "
                     "active against most ALK resistance mutations including G1202R",
        "key_trials": ["CROWN"],
        "resistance_mechanisms": ["compound ALK mutations", "ALK-independent bypass pathways"],
        "evidence_level": "A",
    },
    "alectinib": {
        "drug_name": "alectinib",
        "brand_name": "Alecensa",
        "category": "targeted",
        "targets": ["ALK"],
        "approved_indications": ["ALK+ NSCLC (1L and 2L+)"],
        "mechanism": "2nd-generation ALK inhibitor with CNS penetration",
        "key_trials": ["ALEX", "J-ALEX", "ALUR"],
        "resistance_mechanisms": ["G1202R", "I1171T/N/S", "V1180L", "compound mutations"],
        "evidence_level": "A",
    },
    "dostarlimab": {
        "drug_name": "dostarlimab",
        "brand_name": "Jemperli",
        "category": "immunotherapy",
        "targets": ["PD-1"],
        "approved_indications": ["dMMR recurrent/advanced endometrial cancer",
                                  "dMMR recurrent/advanced solid tumors"],
        "mechanism": "Anti-PD-1 monoclonal antibody",
        "key_trials": ["GARNET"],
        "resistance_mechanisms": ["B2M loss", "JAK1/2 mutations"],
        "evidence_level": "A",
    },
    "sacituzumab_govitecan": {
        "drug_name": "sacituzumab govitecan",
        "brand_name": "Trodelvy",
        "category": "targeted (ADC)",
        "targets": ["TROP2"],
        "approved_indications": ["mTNBC (2L+)", "HR+/HER2- mBC (2L+)", "mUC (2L+)"],
        "mechanism": "TROP2-directed ADC with SN-38 (topoisomerase I inhibitor) payload; "
                     "moderate DAR of ~7.6",
        "key_trials": ["ASCENT", "TROPiCS-02", "TROPHY-U-01"],
        "resistance_mechanisms": ["TROP2 downregulation", "TOP1 mutations", "drug efflux"],
        "evidence_level": "A",
    },
    "enfortumab_vedotin": {
        "drug_name": "enfortumab vedotin",
        "brand_name": "Padcev",
        "category": "targeted (ADC)",
        "targets": ["Nectin-4"],
        "approved_indications": ["locally advanced/mUC (1L with pembro, 2L+ mono)"],
        "mechanism": "Nectin-4-directed ADC with MMAE payload; disrupts microtubule assembly",
        "key_trials": ["EV-301", "EV-302/KEYNOTE-A39"],
        "resistance_mechanisms": ["Nectin-4 downregulation", "drug efflux"],
        "evidence_level": "A",
    },
    "mirvetuximab_soravtansine": {
        "drug_name": "mirvetuximab soravtansine",
        "brand_name": "Elahere",
        "category": "targeted (ADC)",
        "targets": ["FRα"],
        "approved_indications": ["FRα-high platinum-resistant ovarian cancer"],
        "mechanism": "FRα-directed ADC with DM4 (maytansinoid) payload",
        "key_trials": ["SORAYA", "MIRASOL"],
        "resistance_mechanisms": ["FRα loss", "tubulin mutations"],
        "evidence_level": "A",
    },
    "zolbetuximab": {
        "drug_name": "zolbetuximab",
        "brand_name": "Vyloy",
        "category": "targeted",
        "targets": ["CLDN18.2"],
        "approved_indications": ["CLDN18.2+ HER2- gastric/GEJ adenocarcinoma (1L + chemo)"],
        "mechanism": "Anti-CLDN18.2 chimeric IgG1 mAb, ADCC and CDC",
        "key_trials": ["SPOTLIGHT", "GLOW"],
        "resistance_mechanisms": ["CLDN18.2 downregulation", "isoform switching"],
        "evidence_level": "A",
    },
    "tarlatamab": {
        "drug_name": "tarlatamab",
        "brand_name": "Imdelltra",
        "category": "immunotherapy (bispecific)",
        "targets": ["DLL3", "CD3"],
        "approved_indications": ["extensive-stage SCLC (2L+)"],
        "mechanism": "DLL3xCD3 bispecific T-cell engager (BiTE) — redirects T cells to "
                     "DLL3-expressing tumor cells",
        "key_trials": ["DeLLphi-301"],
        "resistance_mechanisms": ["DLL3 loss", "T cell exhaustion"],
        "evidence_level": "B",
    },
    "tazemetostat": {
        "drug_name": "tazemetostat",
        "brand_name": "Tazverik",
        "category": "targeted",
        "targets": ["EZH2"],
        "approved_indications": ["EZH2-mutant relapsed/refractory FL",
                                  "INI1-negative epithelioid sarcoma"],
        "mechanism": "Selective EZH2 inhibitor — blocks H3K27 trimethylation, derepresses "
                     "tumor suppressors",
        "key_trials": ["NCT01897571 (E7438-G000-101)"],
        "resistance_mechanisms": ["SWI/SNF rewiring", "secondary EZH2 mutations"],
        "evidence_level": "A",
    },
    "ipilimumab": {
        "drug_name": "ipilimumab",
        "brand_name": "Yervoy",
        "category": "immunotherapy",
        "targets": ["CTLA-4"],
        "approved_indications": ["melanoma", "RCC (+ nivolumab)", "CRC MSI-H (+ nivolumab)",
                                  "HCC (+ nivolumab)", "NSCLC (+ nivolumab + chemo)",
                                  "mesothelioma (+ nivolumab)"],
        "mechanism": "Anti-CTLA-4 fully human IgG1 mAb, augments T cell activation/proliferation "
                     "by blocking CTLA-4 inhibitory signaling",
        "key_trials": ["CheckMate-067", "CheckMate-214", "CheckMate-142",
                        "CheckMate-8HW", "CheckMate-743"],
        "resistance_mechanisms": ["Treg rebound", "cold TME", "compensatory checkpoints"],
        "evidence_level": "A",
    },
    "capivasertib": {
        "drug_name": "capivasertib",
        "brand_name": "Truqap",
        "category": "targeted",
        "targets": ["AKT1/2/3"],
        "approved_indications": ["PIK3CA/AKT1/PTEN-altered HR+/HER2- advanced breast cancer "
                                  "(+ fulvestrant)"],
        "mechanism": "Selective pan-AKT kinase inhibitor",
        "key_trials": ["CAPItello-291"],
        "resistance_mechanisms": ["mTOR mutations", "MAPK bypass", "AKT3 upregulation"],
        "evidence_level": "A",
    },
    "elacestrant": {
        "drug_name": "elacestrant",
        "brand_name": "Orserdu",
        "category": "targeted (hormonal)",
        "targets": ["ESR1"],
        "approved_indications": ["ESR1-mutant ER+/HER2- advanced breast cancer (2L+)"],
        "mechanism": "Oral selective estrogen receptor degrader (SERD) — degrades both "
                     "wild-type and mutant ER",
        "key_trials": ["EMERALD"],
        "resistance_mechanisms": ["compound ESR1 mutations", "CDK4/6 pathway bypass"],
        "evidence_level": "A",
    },
    "repotrectinib": {
        "drug_name": "repotrectinib",
        "brand_name": "Augtyro",
        "category": "targeted",
        "targets": ["ROS1", "NTRK1/2/3", "ALK"],
        "approved_indications": ["ROS1+ NSCLC (1L and previously treated)"],
        "mechanism": "Next-gen macrocyclic TKI designed to overcome ROS1 G2032R and NTRK "
                     "solvent-front resistance mutations",
        "key_trials": ["TRIDENT-1"],
        "resistance_mechanisms": ["compound ROS1 mutations", "bypass pathway activation"],
        "evidence_level": "A",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# 3. RESISTANCE MAP (~20 entries)
# ═══════════════════════════════════════════════════════════════════════

RESISTANCE_MAP: Dict[str, List[Dict]] = {
    "EGFR TKIs (1st/2nd gen)": [
        {"mutation": "T790M", "gene": "EGFR", "frequency": "~50-60%",
         "bypass": None, "next_line": ["osimertinib"]},
        {"mutation": "MET amplification", "gene": "MET", "frequency": "~5-10%",
         "bypass": "PI3K/AKT", "next_line": ["osimertinib + MET inhibitor"]},
        {"mutation": "HER2 amplification", "gene": "ERBB2", "frequency": "~5%",
         "bypass": "HER2 signaling", "next_line": ["T-DXd"]},
        {"mutation": "Small cell transformation", "gene": "RB1/TP53 loss", "frequency": "~3-10%",
         "bypass": "lineage switch", "next_line": ["platinum + etoposide"]},
    ],
    "osimertinib": [
        {"mutation": "C797S", "gene": "EGFR", "frequency": "~10-15%",
         "bypass": None, "next_line": ["combination EGFR strategies", "amivantamab + lazertinib"]},
        {"mutation": "MET amplification", "gene": "MET", "frequency": "~15-20%",
         "bypass": "PI3K/AKT", "next_line": ["osimertinib + savolitinib"]},
        {"mutation": "BRAF V600E", "gene": "BRAF", "frequency": "~3%",
         "bypass": "MAPK", "next_line": ["BRAF+MEK inhibitors"]},
    ],
    "BRAF V600 inhibitors": [
        {"mutation": "MEK1/2 mutations", "gene": "MAP2K1/2", "frequency": "~5-10%",
         "bypass": "MAPK reactivation", "next_line": ["MEK inhibitor switch"]},
        {"mutation": "NRAS mutations", "gene": "NRAS", "frequency": "~10-20%",
         "bypass": "MAPK reactivation", "next_line": ["immunotherapy"]},
        {"mutation": "BRAF amplification", "gene": "BRAF", "frequency": "~10%",
         "bypass": "MAPK overdrive", "next_line": ["dose escalation", "immunotherapy"]},
    ],
    "ALK TKIs": [
        {"mutation": "G1202R", "gene": "ALK", "frequency": "~20-30% (post-crizotinib)",
         "bypass": None, "next_line": ["lorlatinib"]},
        {"mutation": "compound mutations", "gene": "ALK", "frequency": "~10-15% (post-lorlatinib)",
         "bypass": None, "next_line": ["chemotherapy", "clinical trials"]},
        {"mutation": "MET amplification", "gene": "MET", "frequency": "~5%",
         "bypass": "PI3K", "next_line": ["ALK TKI + MET inhibitor"]},
    ],
    "KRAS G12C inhibitors": [
        {"mutation": "Y96D", "gene": "KRAS", "frequency": "~5-10%",
         "bypass": None, "next_line": ["next-gen KRAS G12C inhibitors"]},
        {"mutation": "Secondary KRAS mutations", "gene": "KRAS", "frequency": "~10%",
         "bypass": "RAS-MAPK", "next_line": ["SOS1 inhibitors", "SHP2 inhibitors"]},
        {"mutation": "MET amplification", "gene": "MET", "frequency": "~10%",
         "bypass": "RTK bypass", "next_line": ["KRAS G12Ci + MET inhibitor"]},
    ],
    "PARP inhibitors": [
        {"mutation": "BRCA reversion mutations", "gene": "BRCA1/2", "frequency": "~20-30%",
         "bypass": "HR restoration", "next_line": ["platinum rechallenge", "immunotherapy"]},
        {"mutation": "53BP1/RIF1 loss", "gene": "TP53BP1", "frequency": "~5%",
         "bypass": "NHEJ restoration", "next_line": ["immunotherapy"]},
        {"mutation": "Drug efflux (ABCB1)", "gene": "ABCB1", "frequency": "~5%",
         "bypass": "drug clearance", "next_line": ["switch PARP inhibitor"]},
    ],
    "anti-PD-1/PD-L1": [
        {"mutation": "B2M loss", "gene": "B2M", "frequency": "~5-10%",
         "bypass": "MHC-I loss", "next_line": ["CTLA-4 combination", "cellular therapy"]},
        {"mutation": "JAK1/2 mutations", "gene": "JAK1/JAK2", "frequency": "~5%",
         "bypass": "IFN-gamma signaling loss", "next_line": ["cellular therapy"]},
        {"mutation": "STK11 co-mutation", "gene": "STK11", "frequency": "~15-20% (NSCLC)",
         "bypass": "cold microenvironment", "next_line": ["combination strategies"]},
    ],
    "CDK4/6 inhibitors": [
        {"mutation": "RB1 loss", "gene": "RB1", "frequency": "~5-10%",
         "bypass": "target loss", "next_line": ["chemotherapy", "PI3K inhibitors"]},
        {"mutation": "CCNE1 amplification", "gene": "CCNE1", "frequency": "~10%",
         "bypass": "CDK2-driven G1/S", "next_line": ["chemotherapy"]},
        {"mutation": "CDK6 amplification", "gene": "CDK6", "frequency": "~5%",
         "bypass": "target amplification", "next_line": ["dose escalation"]},
    ],
    "HER2 ADC (T-DXd)": [
        {"mutation": "HER2 downregulation", "gene": "ERBB2", "frequency": "~15-20%",
         "bypass": "antigen loss", "next_line": ["switch ADC target or chemo"]},
        {"mutation": "SLX4-dependent DNA repair", "gene": "SLX4", "frequency": "~5-10%",
         "bypass": "enhanced DNA repair", "next_line": ["PARP combination"]},
        {"mutation": "Drug efflux (ABCB1)", "gene": "ABCB1", "frequency": "~5%",
         "bypass": "drug clearance", "next_line": ["alternative ADC or payload"]},
    ],
    "anti-CTLA-4": [
        {"mutation": "Treg rebound expansion", "gene": "FOXP3+ Tregs", "frequency": "~20%",
         "bypass": "immunosuppressive TME", "next_line": ["anti-PD-1 combination", "Treg depletion"]},
        {"mutation": "ICOS upregulation", "gene": "ICOS", "frequency": "~10%",
         "bypass": "alternative co-stimulation", "next_line": ["anti-ICOS agents"]},
        {"mutation": "Cold tumor microenvironment", "gene": "multiple", "frequency": "~30%",
         "bypass": "insufficient T cell infiltration", "next_line": ["radiation + IO", "oncolytic virus"]},
    ],
    "PI3K inhibitors": [
        {"mutation": "Insulin-feedback hyperglycemia", "gene": "insulin/IGF1R", "frequency": "~60% (clinical)",
         "bypass": "insulin receptor activation", "next_line": ["dietary management + metformin"]},
        {"mutation": "KRAS activation", "gene": "KRAS", "frequency": "~5-10%",
         "bypass": "MAPK bypass", "next_line": ["PI3Ki + MEK inhibitor"]},
        {"mutation": "mTOR pathway bypass", "gene": "MTOR", "frequency": "~5%",
         "bypass": "downstream activation", "next_line": ["PI3Ki + mTOR inhibitor"]},
    ],
    "FGFR inhibitors": [
        {"mutation": "V565I gatekeeper mutation", "gene": "FGFR2", "frequency": "~10-15%",
         "bypass": "reduced drug binding", "next_line": ["next-gen FGFR inhibitor"]},
        {"mutation": "Polyclonal FGFR2 secondaries", "gene": "FGFR2", "frequency": "~20%",
         "bypass": "multiple resistance clones", "next_line": ["combination strategies"]},
        {"mutation": "MAPK reactivation", "gene": "KRAS/BRAF", "frequency": "~10%",
         "bypass": "bypass pathway", "next_line": ["FGFRi + MEK inhibitor"]},
    ],
}


# ═══════════════════════════════════════════════════════════════════════
# 4. PATHWAY MAP (~10 pathways)
# ═══════════════════════════════════════════════════════════════════════

PATHWAY_MAP: Dict[str, Dict] = {
    "MAPK": {
        "pathway_name": "MAPK (RAS-RAF-MEK-ERK)",
        "key_genes": ["KRAS", "NRAS", "HRAS", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "NF1"],
        "therapeutic_targets": ["BRAF V600 (vemurafenib/dabrafenib)", "MEK1/2 (trametinib/cobimetinib)",
                                 "KRAS G12C (sotorasib/adagrasib)", "ERK1/2 (ulixertinib)"],
        "cross_talk": ["PI3K/AKT (feedback activation upon MAPK inhibition)", "WNT", "cell cycle (cyclin D1)"],
        "clinical_relevance": "Most frequently altered oncogenic pathway. MAPK reactivation is the dominant "
                              "resistance mechanism to targeted therapies in this pathway.",
    },
    "PI3K_AKT_mTOR": {
        "pathway_name": "PI3K/AKT/mTOR",
        "key_genes": ["PIK3CA", "PIK3R1", "PTEN", "AKT1", "AKT2", "MTOR", "TSC1", "TSC2", "STK11"],
        "therapeutic_targets": ["PI3Kα (alpelisib)", "AKT (capivasertib)", "mTOR (everolimus/temsirolimus)",
                                 "dual PI3K/mTOR (future)"],
        "cross_talk": ["MAPK (reciprocal feedback)", "estrogen receptor signaling", "cell cycle"],
        "clinical_relevance": "Activated in >50% of cancers through PIK3CA mutation, PTEN loss, or AKT activation. "
                              "Critical for treatment resistance across tumor types.",
    },
    "DDR": {
        "pathway_name": "DNA Damage Response",
        "key_genes": ["BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "PALB2", "RAD51",
                       "POLE", "MSH2", "MLH1", "MSH6", "PMS2"],
        "therapeutic_targets": ["PARP1/2 (olaparib/rucaparib/niraparib)", "ATR (berzosertib)",
                                 "WEE1 (adavosertib)", "CHK1 (prexasertib)"],
        "cross_talk": ["cell cycle checkpoints", "immune response (cGAS-STING)", "apoptosis"],
        "clinical_relevance": "HRD creates synthetic lethality with PARP inhibitors. MSI-H/dMMR predicts "
                              "immunotherapy response. Foundation of precision oncology biomarker-guided therapy.",
    },
    "cell_cycle": {
        "pathway_name": "Cell Cycle (CDK-RB-E2F)",
        "key_genes": ["CDK4", "CDK6", "CDKN2A", "CCND1", "CCNE1", "RB1", "CDK2", "TP53", "MDM2"],
        "therapeutic_targets": ["CDK4/6 (palbociclib/ribociclib/abemaciclib)", "CDK2 (emerging)",
                                 "MDM2 (nutlin derivatives)"],
        "cross_talk": ["estrogen receptor (ER drives CCND1)", "MAPK (drives CCND1)",
                       "PI3K/AKT (drives cell cycle entry)"],
        "clinical_relevance": "CDK4/6 inhibitors transformed HR+ breast cancer treatment. RB1 loss is a "
                              "key resistance mechanism. CCNE1 amplification predicts CDK4/6i resistance.",
    },
    "WNT": {
        "pathway_name": "WNT/β-catenin",
        "key_genes": ["APC", "CTNNB1", "AXIN1", "AXIN2", "RNF43", "RSPO2", "RSPO3"],
        "therapeutic_targets": ["porcupine inhibitors (emerging)", "tankyrase inhibitors (emerging)"],
        "cross_talk": ["MAPK", "PI3K", "Notch"],
        "clinical_relevance": "WNT activation (APC loss) initiates >80% of CRC. β-catenin activation "
                              "may confer immunotherapy resistance. Limited druggable targets currently.",
    },
    "JAK_STAT": {
        "pathway_name": "JAK-STAT Signaling",
        "key_genes": ["JAK1", "JAK2", "JAK3", "STAT3", "STAT5", "SOCS1"],
        "therapeutic_targets": ["JAK1/2 (ruxolitinib)", "JAK2 (fedratinib)", "STAT3 (emerging)"],
        "cross_talk": ["immune response", "MAPK", "PI3K"],
        "clinical_relevance": "JAK2 V617F defines myeloproliferative neoplasms. JAK1/2 loss-of-function "
                              "mutations confer immunotherapy resistance through impaired IFN-gamma signaling.",
    },
    "apoptosis": {
        "pathway_name": "Apoptosis (BCL-2 Family)",
        "key_genes": ["BCL2", "BCL2L1", "MCL1", "BAX", "BAK1", "BIM", "TP53"],
        "therapeutic_targets": ["BCL-2 (venetoclax)", "MCL-1 (emerging)", "BCL-XL (emerging)"],
        "cross_talk": ["TP53/MDM2", "cell cycle", "DDR"],
        "clinical_relevance": "Venetoclax + azacitidine is standard for unfit AML. BCL-2 dependence in CLL "
                              "makes venetoclax-based regimens highly effective.",
    },
    "angiogenesis": {
        "pathway_name": "Angiogenesis (VEGF/VEGFR)",
        "key_genes": ["VEGFA", "KDR", "FLT1", "FLT4", "PDGFRA", "PDGFRB"],
        "therapeutic_targets": ["VEGF-A (bevacizumab)", "VEGFR (sunitinib/pazopanib/cabozantinib/lenvatinib)",
                                 "multi-kinase (sorafenib)"],
        "cross_talk": ["HIF pathway", "MAPK", "PI3K"],
        "clinical_relevance": "Anti-angiogenics are backbone therapy in RCC, HCC, CRC, ovarian, and cervical cancer. "
                              "Combination with immunotherapy (lenvatinib+pembrolizumab) is paradigm-shifting.",
    },
    "Notch": {
        "pathway_name": "Notch Signaling",
        "key_genes": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL3", "JAG1"],
        "therapeutic_targets": ["DLL3 (tarlatamab — BiTE)", "gamma-secretase inhibitors (emerging)"],
        "cross_talk": ["WNT", "Hedgehog", "MYC"],
        "clinical_relevance": "NOTCH1 mutations occur in T-ALL, CLL. DLL3 is targetable in SCLC with "
                              "tarlatamab (bispecific T-cell engager). Notch is context-dependent — oncogenic "
                              "in some cancers, tumor-suppressive in others.",
    },
    "Hedgehog": {
        "pathway_name": "Hedgehog (HH) Signaling",
        "key_genes": ["PTCH1", "SMO", "GLI1", "GLI2", "SUFU"],
        "therapeutic_targets": ["SMO (vismodegib/sonidegib)"],
        "cross_talk": ["WNT", "Notch", "PI3K"],
        "clinical_relevance": "Vismodegib and sonidegib are approved for advanced basal cell carcinoma. "
                              "PTCH1 loss-of-function mutations activate pathway constitutively.",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# 5. BIOMARKER PANELS (~15 entries)
# ═══════════════════════════════════════════════════════════════════════

BIOMARKER_PANELS: Dict[str, Dict] = {
    "TMB-H": {
        "name": "Tumor Mutational Burden — High",
        "type": "predictive",
        "testing_method": "NGS panel (F1CDx, MSK-IMPACT) or WES",
        "clinical_cutoff": "≥10 mutations/Mb",
        "predictive_for": ["pembrolizumab (tissue-agnostic)"],
        "cancer_types": ["tissue-agnostic"],
        "evidence_level": "A",
        "description": "TMB-H (≥10 mut/Mb) is an FDA-approved tissue-agnostic biomarker for pembrolizumab. "
                       "Higher TMB generally correlates with neoantigen load and immunotherapy benefit.",
    },
    "MSI-H": {
        "name": "Microsatellite Instability — High / dMMR",
        "type": "predictive",
        "testing_method": "IHC (MMR proteins), PCR (microsatellite markers), NGS",
        "clinical_cutoff": "Loss of MLH1/MSH2/MSH6/PMS2 by IHC or MSI by PCR/NGS",
        "predictive_for": ["pembrolizumab", "nivolumab", "dostarlimab", "nivolumab + ipilimumab"],
        "cancer_types": ["tissue-agnostic", "colorectal", "endometrial", "gastric"],
        "evidence_level": "A",
        "description": "MSI-H/dMMR is the strongest predictive biomarker for immunotherapy response. "
                       "First tissue-agnostic FDA approval. Present in ~15% of CRC, ~25% of endometrial.",
    },
    "PD-L1": {
        "name": "PD-L1 Expression (TPS/CPS)",
        "type": "predictive",
        "testing_method": "IHC (22C3, 28-8, SP142, SP263 antibodies)",
        "clinical_cutoff": "TPS≥50% (1L NSCLC mono), TPS≥1% (1L NSCLC + chemo), CPS≥1-10 (various)",
        "predictive_for": ["pembrolizumab", "atezolizumab", "durvalumab"],
        "cancer_types": ["NSCLC", "head and neck", "gastric", "cervical", "TNBC", "urothelial"],
        "evidence_level": "A",
        "description": "PD-L1 TPS≥50% identifies NSCLC patients who benefit from pembrolizumab monotherapy. "
                       "Scoring systems (TPS vs CPS) and cutoffs vary by indication and assay.",
    },
    "HRD": {
        "name": "Homologous Recombination Deficiency",
        "type": "predictive",
        "testing_method": "Genomic scar assays (Myriad myChoice HRD, F1CDx LOH)",
        "clinical_cutoff": "myChoice HRD score ≥42, F1CDx LOH ≥16%",
        "predictive_for": ["PARP inhibitors (olaparib, niraparib)"],
        "cancer_types": ["ovarian", "breast", "prostate", "pancreatic"],
        "evidence_level": "A",
        "description": "HRD status predicts PARP inhibitor benefit beyond BRCA mutations. Composite scores "
                       "incorporate LOH, TAI, and LST. Most validated in ovarian cancer (PAOLA-1, PRIMA).",
    },
    "BRCA_status": {
        "name": "BRCA1/2 Mutation Status",
        "type": "predictive",
        "testing_method": "NGS panel, germline and somatic testing",
        "clinical_cutoff": "Pathogenic/likely pathogenic variants",
        "predictive_for": ["PARP inhibitors", "platinum chemotherapy"],
        "cancer_types": ["breast", "ovarian", "prostate", "pancreatic"],
        "evidence_level": "A",
        "description": "BRCA1/2 mutations are the strongest predictors of PARP inhibitor response. "
                       "Both germline and somatic mutations confer benefit. Platinum sensitivity often correlates.",
    },
    "EGFR_mutation": {
        "name": "EGFR Activating Mutations",
        "type": "predictive",
        "testing_method": "NGS, PCR (cobas), liquid biopsy (ctDNA)",
        "clinical_cutoff": "Sensitizing mutations (L858R, exon 19 del, others)",
        "predictive_for": ["osimertinib", "erlotinib", "gefitinib", "afatinib"],
        "cancer_types": ["NSCLC"],
        "evidence_level": "A",
        "description": "EGFR mutations are the most common actionable target in NSCLC (~15-20% overall, "
                       "~50% in Asian never-smokers). Testing at diagnosis is mandatory per NCCN guidelines.",
    },
    "ALK_fusion": {
        "name": "ALK Rearrangement",
        "type": "predictive",
        "testing_method": "FISH, IHC (D5F3), NGS",
        "clinical_cutoff": "Presence of ALK fusion",
        "predictive_for": ["alectinib", "lorlatinib", "brigatinib", "crizotinib"],
        "cancer_types": ["NSCLC"],
        "evidence_level": "A",
        "description": "ALK fusions occur in ~5% of NSCLC, typically in younger never/light-smokers. "
                       "5-year OS >60% with sequential ALK TKI therapy.",
    },
    "ctDNA_MRD": {
        "name": "Circulating Tumor DNA / Minimal Residual Disease",
        "type": "monitoring",
        "testing_method": "Liquid biopsy (Guardant, Signatera, FoundationOne Liquid CDx)",
        "clinical_cutoff": "ctDNA detected vs not detected (binary); variant allele frequency",
        "predictive_for": ["recurrence risk", "treatment response monitoring"],
        "cancer_types": ["CRC", "breast", "NSCLC", "bladder"],
        "evidence_level": "B",
        "description": "ctDNA MRD detection post-surgery identifies patients at high recurrence risk. "
                       "ctDNA clearance during treatment correlates with response. Actively studied for "
                       "therapy escalation/de-escalation decisions (CIRCULATE, DYNAMIC trials).",
    },
    "tumor_fraction": {
        "name": "Tumor Fraction (ctDNA)",
        "type": "monitoring",
        "testing_method": "Liquid biopsy, whole-genome sequencing of cfDNA",
        "clinical_cutoff": "Variable; typically >1% for reliable variant calling",
        "predictive_for": ["treatment response", "tumor burden estimation"],
        "cancer_types": ["pan-cancer"],
        "evidence_level": "C",
        "description": "Tumor fraction measures the proportion of cfDNA derived from tumor cells. "
                       "Low tumor fraction can lead to false negatives in liquid biopsy testing.",
    },
    "BRAF_V600": {
        "name": "BRAF V600 Mutation",
        "type": "predictive",
        "testing_method": "NGS, PCR (cobas), IHC (VE1)",
        "clinical_cutoff": "BRAF V600E/K/D/R detected",
        "predictive_for": ["dabrafenib+trametinib", "vemurafenib+cobimetinib",
                            "encorafenib+binimetinib", "encorafenib+cetuximab (CRC)"],
        "cancer_types": ["melanoma", "NSCLC", "CRC", "thyroid", "hairy cell leukemia"],
        "evidence_level": "A",
        "description": "BRAF V600E is highly actionable across tumor types. Treatment approach differs: "
                       "melanoma uses BRAF+MEK, CRC requires BRAF+MEK+anti-EGFR triplet.",
    },
    "KRAS_G12C": {
        "name": "KRAS G12C Mutation",
        "type": "predictive",
        "testing_method": "NGS, PCR",
        "clinical_cutoff": "KRAS G12C detected",
        "predictive_for": ["sotorasib", "adagrasib"],
        "cancer_types": ["NSCLC", "CRC"],
        "evidence_level": "A",
        "description": "KRAS G12C occurs in ~13% of NSCLC and ~3% of CRC. First directly targeted "
                       "KRAS mutations in clinical practice. sotorasib and adagrasib FDA-approved for NSCLC.",
    },
    "RET_fusion": {
        "name": "RET Fusion/Mutation",
        "type": "predictive",
        "testing_method": "NGS (DNA and RNA-based), FISH",
        "clinical_cutoff": "RET fusion or activating mutation detected",
        "predictive_for": ["selpercatinib", "pralsetinib"],
        "cancer_types": ["NSCLC", "thyroid"],
        "evidence_level": "A",
        "description": "RET fusions occur in ~1-2% of NSCLC. RNA-based NGS recommended as some fusions "
                       "may be missed by DNA-only testing.",
    },
    "NTRK_fusion": {
        "name": "NTRK Fusion",
        "type": "predictive",
        "testing_method": "NGS (RNA-based preferred), IHC (screening), FISH",
        "clinical_cutoff": "NTRK1/2/3 fusion detected",
        "predictive_for": ["larotrectinib", "entrectinib"],
        "cancer_types": ["tissue-agnostic"],
        "evidence_level": "A",
        "description": "NTRK fusions are rare (<1% pan-cancer) but highly actionable with ORR >75%. "
                       "Tissue-agnostic FDA approval. RNA-based testing increases sensitivity.",
    },
    "HER2_amp": {
        "name": "HER2 Amplification/Overexpression",
        "type": "predictive",
        "testing_method": "IHC (0/1+/2+/3+), FISH, NGS",
        "clinical_cutoff": "IHC 3+ or IHC 2+/FISH-amplified; HER2-low = IHC 1+ or IHC 2+/FISH-negative",
        "predictive_for": ["trastuzumab", "pertuzumab", "T-DXd", "tucatinib"],
        "cancer_types": ["breast", "gastric", "CRC", "NSCLC", "endometrial"],
        "evidence_level": "A",
        "description": "HER2 testing defines a major breast cancer subtype (~20%). The HER2-low paradigm "
                       "(IHC 1+/2+, FISH-negative) emerged with T-DXd (DESTINY-Breast04).",
    },
    "FRalpha": {
        "name": "Folate Receptor Alpha",
        "type": "therapeutic_selection",
        "testing_method": "IHC (LDT, Ventana FOLR1-2.1 antibody)",
        "clinical_cutoff": "PS2+ (≥75% membrane staining at any intensity)",
        "predictive_for": ["mirvetuximab soravtansine"],
        "cancer_types": ["ovarian"],
        "evidence_level": "A",
        "description": "FRα expression identifies platinum-resistant ovarian cancer patients likely to "
                       "benefit from mirvetuximab soravtansine. ~40% of platinum-resistant ovarian is "
                       "FRα-high.",
    },
    "Nectin4": {
        "name": "Nectin-4 Expression",
        "type": "therapeutic_selection",
        "testing_method": "IHC",
        "clinical_cutoff": "Broadly expressed in urothelial (>97%), no cutoff required for enfortumab vedotin",
        "predictive_for": ["enfortumab vedotin"],
        "cancer_types": ["urothelial", "bladder"],
        "evidence_level": "A",
        "description": "Nectin-4 is near-universally expressed in urothelial carcinoma, making it an "
                       "ideal ADC target. No companion diagnostic required for enfortumab vedotin.",
    },
    "CLDN18_2": {
        "name": "Claudin 18.2 Expression",
        "type": "predictive",
        "testing_method": "IHC (VENTANA CLDN18 [43-14A] RxDx assay)",
        "clinical_cutoff": "≥2+ in ≥75% of tumor cells",
        "predictive_for": ["zolbetuximab"],
        "cancer_types": ["gastric", "GEJ"],
        "evidence_level": "A",
        "description": "CLDN18.2 IHC is the companion diagnostic for zolbetuximab. ~38% of gastric "
                       "cancers are CLDN18.2-positive at the approved cutoff.",
    },
    "MRD_heme": {
        "name": "Minimal Residual Disease (Hematologic)",
        "type": "monitoring",
        "testing_method": "multiparameter flow cytometry (MFC), NGS (ClonoSEQ)",
        "clinical_cutoff": "MRD-negative (<10⁻⁴ or <10⁻⁵)",
        "predictive_for": ["treatment de-escalation", "prognosis"],
        "cancer_types": ["ALL", "multiple_myeloma", "AML", "CLL"],
        "evidence_level": "B",
        "description": "MRD negativity in heme malignancies is a strong prognostic marker. ClonoSEQ "
                       "(NGS-based) is FDA-cleared for ALL, CLL, and myeloma MRD assessment.",
    },
    "tTMB": {
        "name": "Tissue vs Blood TMB",
        "type": "predictive",
        "testing_method": "tissue NGS (F1CDx) or blood NGS (F1Liquid CDx, Guardant360)",
        "clinical_cutoff": "≥10 mut/Mb (tissue); ≥16 mut/Mb (blood, varies by platform)",
        "predictive_for": ["pembrolizumab (tissue-agnostic)", "immunotherapy broadly"],
        "cancer_types": ["tissue-agnostic"],
        "evidence_level": "B",
        "description": "Blood-based TMB (bTMB) is emerging as a non-invasive alternative to tissue TMB. "
                       "Concordance with tTMB varies (~65-80%). F1Liquid CDx bTMB≥16 correlates with "
                       "clinical benefit.",
    },
    "GIS": {
        "name": "Genomic Instability Score",
        "type": "predictive",
        "testing_method": "Myriad myChoice CDx (LOH + TAI + LST composite)",
        "clinical_cutoff": "GIS ≥42",
        "predictive_for": ["PARP inhibitors (olaparib + bevacizumab in ovarian)"],
        "cancer_types": ["ovarian", "breast"],
        "evidence_level": "A",
        "description": "GIS captures HRD beyond BRCA mutations by measuring genomic scarring. "
                       "myChoice CDx GIS ≥42 (or BRCA-mutant) identifies ovarian cancer patients "
                       "benefiting from olaparib maintenance (PAOLA-1, PRIMA trials).",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# 6. ENTITY ALIASES
# ═══════════════════════════════════════════════════════════════════════

ENTITY_ALIASES: Dict[str, str] = {
    # Drug brand → generic
    "keytruda": "pembrolizumab", "opdivo": "nivolumab", "tagrisso": "osimertinib",
    "zelboraf": "vemurafenib", "tafinlar": "dabrafenib", "mekinist": "trametinib",
    "lumakras": "sotorasib", "krazati": "adagrasib", "vitrakvi": "larotrectinib",
    "rozlytrek": "entrectinib", "retevmo": "selpercatinib", "gavreto": "pralsetinib",
    "tabrecta": "capmatinib", "tepmetko": "tepotinib", "balversa": "erdafitinib",
    "piqray": "alpelisib", "lynparza": "olaparib", "rubraca": "rucaparib",
    "zejula": "niraparib", "talzenna": "talazoparib", "tibsovo": "ivosidenib",
    "idhifa": "enasidenib", "enhertu": "trastuzumab_deruxtecan", "t-dxd": "trastuzumab_deruxtecan",
    "ibrance": "palbociclib", "kisqali": "ribociclib", "verzenio": "abemaciclib",
    "lorbrena": "lorlatinib", "alecensa": "alectinib", "xalkori": "crizotinib",
    "jemperli": "dostarlimab",
    # Gene aliases
    "her2": "HER2", "erbb2": "HER2", "egfrviii": "EGFR", "lkb1": "STK11",
    "p16": "CDKN2A", "p14arf": "CDKN2A", "ntrk": "NTRK",
    # Cancer type aliases
    "non-small cell": "NSCLC", "non small cell": "NSCLC", "lung adenocarcinoma": "NSCLC",
    "lung squamous": "NSCLC", "small cell lung": "SCLC", "triple negative": "TNBC",
    "crc": "colorectal", "colon cancer": "colorectal", "rectal cancer": "colorectal",
    "rcc": "renal cell", "hcc": "hepatocellular",
    # Biomarker aliases
    "msi high": "MSI-H", "msi-high": "MSI-H", "microsatellite instability": "MSI-H",
    "mmr deficient": "MSI-H", "dmmr": "MSI-H", "tmb high": "TMB-H", "tmb-high": "TMB-H",
    "pd-l1": "PD-L1", "pdl1": "PD-L1", "hrd": "HRD",
    # New drug brand → generic
    "trodelvy": "sacituzumab_govitecan", "padcev": "enfortumab_vedotin",
    "elahere": "mirvetuximab_soravtansine", "vyloy": "zolbetuximab",
    "imdelltra": "tarlatamab", "tazverik": "tazemetostat",
    "yervoy": "ipilimumab", "truqap": "capivasertib",
    "orserdu": "elacestrant", "augtyro": "repotrectinib",
    # New gene aliases
    "claudin 18.2": "CLDN18_2", "claudin18.2": "CLDN18_2", "cldn18": "CLDN18_2",
    "nectin-4": "NECTIN4", "nectin4": "NECTIN4",
    "folate receptor": "FOLR1", "fra": "FOLR1",
    "dll3": "DLL3", "trop-2": "TROP2", "trop2": "TROP2",
    "e-cadherin": "CDH1", "mek1": "MAP2K1",
    "ini1": "SMARCB1", "snf5": "SMARCB1",
    # Pediatric cancer aliases
    "all": "ETV6_RUNX1", "b-all": "ETV6_RUNX1", "b-cell all": "ETV6_RUNX1",
    "t-all": "NOTCH1", "t-cell all": "NOTCH1",
    "neuroblastoma": "MYCN", "nbl": "MYCN",
    "medulloblastoma": "PTCH1", "mb": "PTCH1",
    "rhabdomyosarcoma": "PAX3_FOXO1", "rms": "PAX3_FOXO1",
    "ewing sarcoma": "EWSR1_FLI1", "ewing": "EWSR1_FLI1",
    "wilms tumor": "WT1", "nephroblastoma": "WT1",
    "hepatoblastoma": "CTNNB1_PED", "hbl": "CTNNB1_PED",
    "retinoblastoma": "RB1", "rb": "RB1",
    "infant all": "KMT2A", "mll": "KMT2A", "mll rearranged": "KMT2A",
    "philadelphia": "BCR_ABL1", "ph+": "BCR_ABL1", "ph+ all": "BCR_ABL1",
}


# ═══════════════════════════════════════════════════════════════════════
# 7. PEDIATRIC DOSING PRINCIPLES
# ═══════════════════════════════════════════════════════════════════════

PEDIATRIC_DOSING: Dict[str, Dict] = {
    "general_principles": {
        "name": "Pediatric Dosing — General Principles",
        "description": "Pediatric oncology dosing is BSA-based (mg/m²), not fixed adult "
                       "doses. Weight-based dosing (mg/kg) used for infants <10 kg. "
                       "Pharmacokinetics differ from adults: higher drug clearance relative "
                       "to body weight, different volume of distribution, immature hepatic "
                       "metabolism in neonates/infants.",
        "key_rules": [
            "BSA-based dosing (mg/m²) for most chemotherapy agents",
            "Weight-based (mg/kg) for infants <10 kg or BSA <0.6 m²",
            "Never exceed adult maximum doses regardless of BSA",
            "Dose reductions for organ dysfunction (hepatic, renal)",
            "Age-specific pharmacokinetic considerations",
        ],
    },
    "anthracycline_limits": {
        "name": "Cumulative Anthracycline Limits",
        "description": "Anthracyclines (doxorubicin, daunorubicin, epirubicin, idarubicin) "
                       "cause cumulative cardiotoxicity. Pediatric patients are especially "
                       "vulnerable due to developing myocardium and long life expectancy "
                       "post-treatment.",
        "key_rules": [
            "Cumulative anthracycline limit: 450 mg/m² lifetime (doxorubicin equivalent)",
            "Dexrazoxane cardioprotection recommended for cumulative dose >300 mg/m²",
            "Serial echocardiography (shortening fraction) monitoring required",
            "Late cardiotoxicity can manifest decades after treatment",
            "Conversion: daunorubicin 1:1, epirubicin 0.67:1, idarubicin 5:1 (to doxorubicin equivalent)",
        ],
    },
    "methotrexate_intrathecal": {
        "name": "Methotrexate — Age-Adjusted Intrathecal Dosing",
        "description": "Intrathecal (IT) methotrexate dosing in pediatric ALL is age-adjusted "
                       "(not BSA-based) because CSF volume reaches adult levels by age 3 "
                       "while body surface area continues to grow.",
        "key_rules": [
            "Age 1-2 years: 8 mg IT methotrexate",
            "Age 2-3 years: 10 mg IT methotrexate",
            "Age ≥3 years: 12 mg IT methotrexate",
            "Maximum IT methotrexate dose: 15 mg (adolescents/adults)",
            "Triple intrathecal therapy: methotrexate + cytarabine + hydrocortisone",
            "Leucovorin rescue for high-dose systemic methotrexate (≥1 g/m²)",
        ],
    },
    "vincristine_cap": {
        "name": "Vincristine — Absolute Dose Cap",
        "description": "Vincristine has an absolute dose cap of 2 mg regardless of BSA to "
                       "prevent severe neurotoxicity. This is one of the most critical "
                       "safety rules in pediatric oncology.",
        "key_rules": [
            "Absolute maximum: 2 mg total dose regardless of BSA",
            "Standard dose: 1.5 mg/m² (capped at 2 mg)",
            "Peripheral neuropathy is dose-limiting toxicity",
            "Dose reduction for hepatic dysfunction (bilirubin >3 mg/dL → 50% reduction)",
            "NEVER administer intrathecally (universally fatal)",
        ],
    },
    "asparaginase": {
        "name": "L-Asparaginase — Pediatric Considerations",
        "description": "L-asparaginase is a cornerstone of ALL therapy. Pegylated formulation "
                       "(pegaspargase) is standard. Hypersensitivity, pancreatitis, and "
                       "thrombosis are key toxicities.",
        "key_rules": [
            "Pegaspargase: 2,500 IU/m² IM or IV (standard)",
            "Monitor asparaginase activity levels (therapeutic >0.1 IU/mL)",
            "Switch to Erwinia asparaginase (now recombinant: Rylaze) for hypersensitivity",
            "Monitor for pancreatitis, thrombosis (antithrombin III), hepatotoxicity",
            "Silent inactivation: check asparaginase activity if subtherapeutic response",
        ],
    },
}


def get_pediatric_dosing_context(agent: str) -> str:
    """Return formatted pediatric dosing context for a given drug or principle.

    Args:
        agent: Drug name or dosing topic (e.g. 'vincristine', 'anthracycline',
               'methotrexate', 'general').

    Returns:
        Formatted string with pediatric dosing guidelines.
    """
    agent_lower = agent.lower().strip()

    # Direct key match
    if agent_lower in PEDIATRIC_DOSING:
        entry = PEDIATRIC_DOSING[agent_lower]
    else:
        # Fuzzy match
        entry = None
        for key, val in PEDIATRIC_DOSING.items():
            if agent_lower in key or agent_lower in val["name"].lower():
                entry = val
                break
    if not entry:
        return ""

    lines = [
        f"Pediatric Dosing: {entry['name']}",
        f"  {entry['description']}",
        "  Key Rules:",
    ]
    for rule in entry["key_rules"]:
        lines.append(f"    - {rule}")
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════
# 8. HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════


def classify_variant_actionability(gene: str, variant: str) -> str:
    """Classify a variant's actionability against known ACTIONABLE_TARGETS.

    This is the shared implementation used by case_manager and cross_modal
    to avoid duplicating classification logic.

    Args:
        gene: Gene symbol (e.g. "EGFR", "BRAF").
        variant: Variant description (e.g. "L858R", "V600E").

    Returns:
        Evidence level string ("A", "B", "C", "D") or "VUS" if not recognized.
    """
    gene_upper = gene.upper().strip()
    if gene_upper not in ACTIONABLE_TARGETS:
        return "VUS"

    target_info = ACTIONABLE_TARGETS[gene_upper]

    # Check key_variants list for a match
    key_variants = target_info.get("key_variants", [])
    variant_upper = variant.upper().strip()
    for known_var in key_variants:
        if known_var.upper() in variant_upper or variant_upper in known_var.upper():
            return target_info.get("evidence_level", "C")

    # Also support legacy "variants" dict format
    known_variants = target_info.get("variants", {})
    if isinstance(known_variants, dict):
        for known_var, details in known_variants.items():
            if known_var.upper() in variant_upper:
                ev = details.get("evidence_level", "C") if isinstance(details, dict) else "C"
                return ev

    # Gene-level actionability (some genes are actionable regardless of variant)
    if target_info.get("gene_level_actionable", False):
        return target_info.get("default_evidence_level", "C")

    return "VUS"


def get_target_context(gene: str) -> str:
    """Return formatted knowledge context for an actionable target gene."""
    target = ACTIONABLE_TARGETS.get(gene.upper()) or ACTIONABLE_TARGETS.get(gene)
    if not target:
        alias = ENTITY_ALIASES.get(gene.lower(), "").upper()
        target = ACTIONABLE_TARGETS.get(alias)
    if not target:
        return ""

    lines = [
        f"Target: {target['gene']} ({target['full_name']})",
        f"  Cancer types: {', '.join(target['cancer_types'])}",
        f"  Key variants: {', '.join(target['key_variants'][:6])}",
        f"  Targeted therapies: {', '.join(target['targeted_therapies'][:5]) or 'None approved'}",
        f"  Resistance: {', '.join(target['resistance_mutations'][:4]) or 'Unknown'}",
        f"  Pathway: {target['pathway']}",
        f"  Evidence level: {target['evidence_level']}",
        f"  {target['description']}",
    ]
    return "\n".join(lines)


def get_therapy_context(drug: str) -> str:
    """Return formatted knowledge context for a therapy."""
    key = drug.lower().replace(" ", "_").replace("-", "_")
    therapy = THERAPY_MAP.get(key)
    if not therapy:
        resolved = ENTITY_ALIASES.get(drug.lower(), "")
        therapy = THERAPY_MAP.get(resolved)
    if not therapy:
        return ""

    lines = [
        f"Therapy: {therapy['drug_name']} ({therapy['brand_name']})",
        f"  Category: {therapy['category']}",
        f"  Targets: {', '.join(therapy['targets'])}",
        f"  Approved indications: {'; '.join(therapy['approved_indications'][:4])}",
        f"  Mechanism: {therapy['mechanism']}",
        f"  Key trials: {', '.join(therapy['key_trials'][:4])}",
        f"  Resistance: {', '.join(therapy['resistance_mechanisms'][:3])}",
    ]
    return "\n".join(lines)


def get_resistance_context(drug: str) -> str:
    """Return resistance mechanisms for a therapy class."""
    mechanisms = RESISTANCE_MAP.get(drug)
    if not mechanisms:
        for key in RESISTANCE_MAP:
            if drug.lower() in key.lower():
                mechanisms = RESISTANCE_MAP[key]
                break
    if not mechanisms:
        return ""

    lines = [f"Resistance mechanisms for {drug}:"]
    for m in mechanisms:
        line = f"  - {m['mutation']} ({m['gene']}): ~{m['frequency']}"
        if m.get("bypass"):
            line += f" | bypass: {m['bypass']}"
        if m.get("next_line"):
            line += f" → {', '.join(m['next_line'])}"
        lines.append(line)
    return "\n".join(lines)


def get_pathway_context(pathway: str) -> str:
    """Return formatted context for an oncogenic pathway."""
    pw = PATHWAY_MAP.get(pathway.upper()) or PATHWAY_MAP.get(pathway)
    if not pw:
        for key, val in PATHWAY_MAP.items():
            if pathway.lower() in key.lower() or pathway.lower() in val["pathway_name"].lower():
                pw = val
                break
    if not pw:
        return ""

    lines = [
        f"Pathway: {pw['pathway_name']}",
        f"  Key genes: {', '.join(pw['key_genes'][:8])}",
        f"  Therapeutic targets: {'; '.join(pw['therapeutic_targets'][:4])}",
        f"  Cross-talk: {', '.join(pw['cross_talk'][:3])}",
        f"  {pw['clinical_relevance']}",
    ]
    return "\n".join(lines)


def get_biomarker_context(biomarker: str) -> str:
    """Return formatted context for a biomarker panel."""
    bm = BIOMARKER_PANELS.get(biomarker)
    if not bm:
        alias = ENTITY_ALIASES.get(biomarker.lower(), "")
        bm = BIOMARKER_PANELS.get(alias)
    if not bm:
        for key, val in BIOMARKER_PANELS.items():
            if biomarker.lower() in key.lower() or biomarker.lower() in val["name"].lower():
                bm = val
                break
    if not bm:
        return ""

    lines = [
        f"Biomarker: {bm['name']}",
        f"  Type: {bm['type']}",
        f"  Testing: {bm['testing_method']}",
        f"  Cutoff: {bm['clinical_cutoff']}",
        f"  Predictive for: {', '.join(bm['predictive_for'][:3])}",
        f"  Cancer types: {', '.join(bm['cancer_types'][:4])}",
        f"  Evidence level: {bm['evidence_level']}",
        f"  {bm['description']}",
    ]
    return "\n".join(lines)


def resolve_comparison_entity(raw: str) -> Optional[Dict]:
    """Resolve a raw entity string to a canonical entity for comparison.

    Searches targets, therapies, pathways, and biomarkers in that order.

    Returns:
        Dict with 'canonical', 'type', 'target' (gene if applicable), or None.
    """
    normalized = raw.strip()
    upper = normalized.upper()

    # Check entity aliases first
    alias = ENTITY_ALIASES.get(normalized.lower())
    if alias:
        normalized = alias
        upper = normalized.upper()

    # Check targets
    if upper in ACTIONABLE_TARGETS:
        t = ACTIONABLE_TARGETS[upper]
        return {"canonical": t["gene"], "type": "target", "target": t["gene"], "data": t}

    # Check therapies
    key = normalized.lower().replace(" ", "_").replace("-", "_")
    if key in THERAPY_MAP:
        t = THERAPY_MAP[key]
        return {"canonical": t["drug_name"], "type": "therapy", "target": None, "data": t}

    # Check pathways
    for pw_key, pw in PATHWAY_MAP.items():
        if upper == pw_key or normalized.lower() in pw["pathway_name"].lower():
            return {"canonical": pw["pathway_name"], "type": "pathway", "target": None, "data": pw}

    # Check biomarkers
    for bm_key, bm in BIOMARKER_PANELS.items():
        if upper == bm_key.upper() or normalized.lower() in bm["name"].lower():
            return {"canonical": bm["name"], "type": "biomarker", "target": None, "data": bm}

    return None


def get_comparison_context(entity_a: Dict, entity_b: Dict) -> str:
    """Generate formatted comparison context between two entities."""
    lines = [f"Comparison: {entity_a['canonical']} vs {entity_b['canonical']}"]

    if entity_a["type"] == "target" and entity_b["type"] == "target":
        a, b = entity_a["data"], entity_b["data"]
        lines.extend([
            f"\n{a['gene']}:",
            f"  Cancer types: {', '.join(a['cancer_types'][:4])}",
            f"  Therapies: {', '.join(a['targeted_therapies'][:3])}",
            f"  Pathway: {a['pathway']}",
            f"\n{b['gene']}:",
            f"  Cancer types: {', '.join(b['cancer_types'][:4])}",
            f"  Therapies: {', '.join(b['targeted_therapies'][:3])}",
            f"  Pathway: {b['pathway']}",
        ])
    elif entity_a["type"] == "therapy" and entity_b["type"] == "therapy":
        a, b = entity_a["data"], entity_b["data"]
        lines.extend([
            f"\n{a['drug_name']} ({a['brand_name']}):",
            f"  Targets: {', '.join(a['targets'])}",
            f"  Category: {a['category']}",
            f"  Indications: {'; '.join(a['approved_indications'][:3])}",
            f"\n{b['drug_name']} ({b['brand_name']}):",
            f"  Targets: {', '.join(b['targets'])}",
            f"  Category: {b['category']}",
            f"  Indications: {'; '.join(b['approved_indications'][:3])}",
        ])
    else:
        for entity in [entity_a, entity_b]:
            if entity["type"] == "target":
                lines.append(f"\n{get_target_context(entity['canonical'])}")
            elif entity["type"] == "therapy":
                lines.append(f"\n{get_therapy_context(entity['canonical'])}")
            elif entity["type"] == "pathway":
                lines.append(f"\n{get_pathway_context(entity['canonical'])}")
            elif entity["type"] == "biomarker":
                lines.append(f"\n{get_biomarker_context(entity['canonical'])}")

    return "\n".join(lines)
