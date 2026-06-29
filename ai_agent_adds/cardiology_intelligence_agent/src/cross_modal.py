"""Cross-modal genomic triggers for the Cardiology Intelligence Agent.

Links cardiac imaging findings, ECG abnormalities, and clinical data to
genomic workup recommendations by querying the shared genomic_evidence
collection (~3.5 M variants indexed via BGE-small-en-v1.5 embeddings in
Milvus).

When a clinical finding meets predefined threshold criteria -- for
example unexplained left ventricular hypertrophy >= 15 mm, or QTc > 480
ms -- the engine emits a CrossModalTrigger that names the relevant gene
panel and builds a vector-search query against genomic_evidence.  This
bridges bedside cardiology with the HCLS AI Factory's genomic pipeline
running on the same DGX Spark, enabling clinicians to move from imaging
to molecular insight in seconds rather than weeks.

Trigger categories:
    - Cardiomyopathies (HCM, DCM, ACM, LVNC, Fabry, ATTR amyloid)
    - Channelopathies (LQTS, Brugada, CPVT)
    - Vascular / lipid (premature CAD / FH, heritable aortopathy)
    - Family screening (first-degree relative with SCD < 40)

Each trigger carries an urgency level (critical / high / moderate / low)
that drives clinical workflow prioritisation and cascade logic.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import hashlib
import logging
import time
from collections import defaultdict
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from src.models import CrossModalTrigger, SeverityLevel

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

GENOMIC_COLLECTION = "genomic_evidence"
TOP_K_PER_QUERY = 8
SCORE_THRESHOLD = 0.38
EMBEDDING_DIM = 384  # BGE-small-en-v1.5

# Urgency string -> SeverityLevel mapping for the trigger map
_URGENCY_MAP: Dict[str, SeverityLevel] = {
    "critical": SeverityLevel.CRITICAL,
    "high": SeverityLevel.HIGH,
    "moderate": SeverityLevel.MODERATE,
    "low": SeverityLevel.LOW,
}

# Ordered urgency levels for sorting (most urgent first)
_URGENCY_PRIORITY: Dict[SeverityLevel, int] = {
    SeverityLevel.CRITICAL: 0,
    SeverityLevel.HIGH: 1,
    SeverityLevel.MODERATE: 2,
    SeverityLevel.LOW: 3,
    SeverityLevel.INFORMATIONAL: 4,
}


# ═══════════════════════════════════════════════════════════════════════
# GENOMIC TRIGGER MAP
# ═══════════════════════════════════════════════════════════════════════
# Comprehensive mapping of clinical cardiac findings to gene panels.
# Each entry specifies the conditions under suspicion, the gene panel
# to order / query, the clinical criteria that must be met, urgency,
# and (optionally) the supporting guideline.

GENOMIC_TRIGGER_MAP: Dict[str, Dict[str, Any]] = {
    # ── Cardiomyopathies ──────────────────────────────────────────────
    "unexplained_lvh": {
        "conditions": ["Hypertrophic Cardiomyopathy"],
        "gene_panel": [
            "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1",
            "ACTC1", "MYL2", "MYL3", "GLA", "PRKAG2", "LAMP2",
        ],
        "criteria": (
            "LV wall thickness >=15mm (or >=13mm with family history) "
            "without identifiable cause"
        ),
        "urgency": "high",
        "guideline": "ACC/AHA 2024 HCM Guidelines",
    },
    "unexplained_dcm": {
        "conditions": ["Dilated Cardiomyopathy"],
        "gene_panel": [
            "TTN", "LMNA", "RBM20", "MYH7", "TNNT2",
            "DSP", "FLNC", "BAG3", "SCN5A", "PLN",
        ],
        "criteria": "LVEF <45% with LV dilation, age <60, no ischemic etiology",
        "urgency": "high",
        "guideline": "ACC/AHA 2022 HF Guidelines",
    },
    "arrhythmogenic_cm": {
        "conditions": ["Arrhythmogenic Cardiomyopathy"],
        "gene_panel": [
            "PKP2", "DSP", "DSG2", "DSC2", "JUP",
            "TMEM43", "PLN", "FLNC", "DES", "LMNA",
        ],
        "criteria": (
            "RV dilation/dysfunction with ventricular arrhythmias, "
            "fibro-fatty replacement on CMR"
        ),
        "urgency": "high",
    },
    "lvnc": {
        "conditions": ["Left Ventricular Non-Compaction"],
        "gene_panel": ["MYH7", "MYBPC3", "TTN", "LMNA", "TAZ", "ACTC1"],
        "criteria": "Non-compacted to compacted ratio >2.3 (diastole) on CMR",
        "urgency": "moderate",
    },
    "cardiac_amyloid": {
        "conditions": ["Cardiac Amyloidosis (ATTR)", "Cardiac Amyloidosis (AL)"],
        "gene_panel": ["TTR"],
        "criteria": (
            "Unexplained LVH + diastolic dysfunction + low voltage ECG, "
            "or positive Tc-99m PYP/DPD scan"
        ),
        "urgency": "high",
    },
    "fabry_disease": {
        "conditions": ["Fabry Disease"],
        "gene_panel": ["GLA"],
        "criteria": (
            "Unexplained LVH + renal dysfunction + neuropathy, "
            "or LVH in young patient"
        ),
        "urgency": "moderate",
    },

    # ── Channelopathies ───────────────────────────────────────────────
    "long_qt": {
        "conditions": ["Long QT Syndrome"],
        "gene_panel": [
            "KCNQ1", "KCNH2", "SCN5A", "KCNJ2",
            "CALM1", "CALM2", "CALM3", "TRDN", "ANK2",
        ],
        "criteria": "QTc >480ms or >460ms with syncope/family hx SCD",
        "urgency": "critical",
    },
    "brugada_pattern": {
        "conditions": ["Brugada Syndrome"],
        "gene_panel": [
            "SCN5A", "CACNA1C", "CACNB2", "SCN1B",
            "SCN2B", "SCN3B", "GPD1L", "HEY2",
        ],
        "criteria": (
            "Type 1 Brugada ECG pattern "
            "(coved ST elevation >=2mm in V1-V3)"
        ),
        "urgency": "critical",
    },
    "cpvt_suspected": {
        "conditions": ["Catecholaminergic Polymorphic VT"],
        "gene_panel": ["RYR2", "CASQ2", "TRDN", "CALM1", "TECRL"],
        "criteria": (
            "Exercise or catecholamine-induced bidirectional/polymorphic VT "
            "in structurally normal heart"
        ),
        "urgency": "critical",
    },

    # ── Vascular / Lipid ─────────────────────────────────────────────
    "premature_cad": {
        "conditions": ["Familial Hypercholesterolemia"],
        "gene_panel": ["LDLR", "PCSK9", "APOB", "APOE", "LDLRAP1"],
        "criteria": "CAD age <55 (male) or <65 (female), or LDL >=190 mg/dL",
        "urgency": "moderate",
    },
    "aortic_dilation": {
        "conditions": [
            "Heritable Thoracic Aortic Disease",
            "Marfan Syndrome",
            "Loeys-Dietz Syndrome",
        ],
        "gene_panel": [
            "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2",
            "MYH11", "COL3A1", "PRKG1", "LOX",
        ],
        "criteria": (
            "Aortic root >=4.0 cm (or Z-score >=2) at age <50 "
            "or with systemic features"
        ),
        "urgency": "high",
    },

    # ── Infiltrative / Metabolic ─────────────────────────────────────
    "cardiac_sarcoidosis": {
        "conditions": ["cardiac sarcoidosis", "heart block", "VT with sarcoidosis"],
        "gene_panel": ["TNF", "IL6R", "BTNL2", "HLA-DRB1"],
        "criteria": (
            "Biopsy-confirmed non-caseating granulomas OR clinical diagnosis "
            "(HRS 2014 criteria). CMR: patchy mid-wall LGE + FDG PET: focal uptake."
        ),
        "urgency": "high",
        "guideline": "2014 HRS Expert Consensus on Arrhythmias in Sarcoidosis",
    },
    "hemochromatosis": {
        "conditions": ["hereditary hemochromatosis", "iron overload cardiomyopathy"],
        "gene_panel": ["HFE", "TFR2", "HJV", "HAMP", "SLC40A1"],
        "criteria": (
            "CMR T2* <20ms (myocardial iron overload) or serum ferritin >1000 "
            "with TSAT >45%"
        ),
        "urgency": "high",
    },

    # ── Pregnancy-Associated ──────────────────────────────────────────
    "peripartum_cardiomyopathy": {
        "conditions": ["peripartum cardiomyopathy", "postpartum heart failure"],
        "gene_panel": [
            "TTN", "MYBPC3", "MYH7", "BAG3", "LMNA", "RBM20", "FLNC",
        ],
        "criteria": (
            "LVEF <45%, no other cause identified, onset in last month of "
            "pregnancy or within 5 months postpartum"
        ),
        "urgency": "high",
        "guideline": "ESC 2018 Guidelines on Cardiovascular Disease in Pregnancy",
    },

    # ── Structural / Morphological ────────────────────────────────────
    "lvnc_phenotype": {
        "conditions": ["left ventricular non-compaction", "LVNC"],
        "gene_panel": [
            "MYH7", "MYBPC3", "TNNT2", "ACTC1", "MIB1",
            "DTNA", "LDB3", "TAZ",
        ],
        "criteria": (
            "NC/C ratio >2.3 (Jenni) or >2.3 in diastole (Petersen), "
            "fractal dimension analysis. TAZ (Barth syndrome) if male."
        ),
        "urgency": "moderate",
        "guideline": "2023 ESC Cardiomyopathy Guidelines",
    },
    "bicuspid_aortopathy": {
        "conditions": [
            "bicuspid aortic valve", "bicuspid aortopathy",
            "ascending aortic aneurysm",
        ],
        "gene_panel": ["NOTCH1", "SMAD6", "GATA5", "GATA4", "NKX2-5", "ROBO4"],
        "criteria": (
            "BAV + aorta >=4.0cm or growth rate >=3mm/yr or family history "
            "of aortic dissection. Urgency escalates to high if >5.0cm or "
            "rapid growth."
        ),
        "urgency": "moderate",
        "guideline": "2022 ACC/AHA Aortic Disease Guidelines",
    },

    # ── Channelopathies (additional) ──────────────────────────────────
    "catecholaminergic_polymorphic_vt": {
        "conditions": [
            "CPVT", "exercise-induced VT",
            "catecholaminergic polymorphic VT",
        ],
        "gene_panel": [
            "RYR2", "CASQ2", "TRDN", "CALM1", "CALM2", "CALM3", "TECRL",
        ],
        "criteria": (
            "Exercise or isoproterenol provoked bidirectional VT, age <40, "
            "structurally normal heart"
        ),
        "urgency": "critical",
        "guideline": "2022 ESC Ventricular Arrhythmia Guidelines",
    },

    # ── Family Screening ─────────────────────────────────────────────
    "scd_family_history": {
        "conditions": ["Sudden Cardiac Death - Family Screening"],
        "gene_panel": [
            "SCN5A", "KCNQ1", "KCNH2", "MYH7", "MYBPC3",
            "PKP2", "DSP", "RYR2", "LMNA", "TNNT2", "TTN",
        ],
        "criteria": (
            "First-degree relative with SCD age <40 "
            "or diagnosed inherited cardiac condition"
        ),
        "urgency": "high",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# FAMILY SCREENING PROTOCOLS
# ═══════════════════════════════════════════════════════════════════════
# Mapping from condition name to the recommended family screening
# protocol: who to screen, what modalities, and interval.

FAMILY_SCREENING_PROTOCOLS: Dict[str, Dict[str, Any]] = {
    "Hypertrophic Cardiomyopathy": {
        "who": "All first-degree relatives",
        "modalities": ["Echocardiography", "12-lead ECG", "Genetic testing"],
        "initial_screening_age": 10,
        "repeat_interval_years": 3,
        "repeat_until_age": 60,
        "notes": (
            "If pathogenic variant identified, cascade genetic testing "
            "replaces serial imaging for genotype-negative relatives."
        ),
        "guideline": "ACC/AHA 2024 HCM Guidelines",
    },
    "Dilated Cardiomyopathy": {
        "who": "All first-degree relatives",
        "modalities": ["Echocardiography", "12-lead ECG"],
        "initial_screening_age": 10,
        "repeat_interval_years": 3,
        "repeat_until_age": 60,
        "notes": (
            "If LMNA variant present, lower threshold for ICD discussion "
            "in relatives who carry the variant."
        ),
        "guideline": "ACC/AHA 2022 HF Guidelines",
    },
    "Arrhythmogenic Cardiomyopathy": {
        "who": "All first-degree relatives",
        "modalities": [
            "12-lead ECG", "Signal-averaged ECG", "Holter monitor",
            "Echocardiography", "Cardiac MRI",
        ],
        "initial_screening_age": 10,
        "repeat_interval_years": 2,
        "repeat_until_age": 60,
        "notes": (
            "Exercise restriction counselling for all variant carriers, "
            "even if phenotype-negative."
        ),
        "guideline": "2019 HRS Expert Consensus on ACM",
    },
    "Long QT Syndrome": {
        "who": "All first-degree relatives",
        "modalities": ["12-lead ECG with QTc measurement"],
        "initial_screening_age": 0,
        "repeat_interval_years": 1,
        "repeat_until_age": 40,
        "notes": (
            "Immediate beta-blocker consideration for genotype-positive "
            "relatives.  Avoid QT-prolonging medications (CredibleMeds list)."
        ),
        "guideline": "2017 AHA/ACC/HRS Ventricular Arrhythmia Guidelines",
    },
    "Brugada Syndrome": {
        "who": "All first-degree relatives",
        "modalities": [
            "12-lead ECG", "Ajmaline/Procainamide provocation test",
        ],
        "initial_screening_age": 16,
        "repeat_interval_years": 2,
        "repeat_until_age": 50,
        "notes": (
            "Fever management protocol for all variant carriers.  "
            "Avoid class I antiarrhythmics and Brugada-provocative drugs."
        ),
        "guideline": "2022 ESC Ventricular Arrhythmia Guidelines",
    },
    "Catecholaminergic Polymorphic VT": {
        "who": "All first-degree relatives",
        "modalities": ["Exercise stress test", "12-lead ECG"],
        "initial_screening_age": 5,
        "repeat_interval_years": 2,
        "repeat_until_age": 40,
        "notes": (
            "Mandatory exercise restriction counselling.  Cascade genetic "
            "testing strongly recommended given high penetrance of RYR2 variants."
        ),
        "guideline": "2017 AHA/ACC/HRS Ventricular Arrhythmia Guidelines",
    },
    "Familial Hypercholesterolemia": {
        "who": "All first-degree relatives",
        "modalities": ["Fasting lipid panel", "Genetic testing"],
        "initial_screening_age": 2,
        "repeat_interval_years": 5,
        "repeat_until_age": None,
        "notes": (
            "Cascade testing is the most cost-effective approach.  "
            "Initiate statin therapy in children age >=10 if FH confirmed."
        ),
        "guideline": "2018 AHA/ACC Cholesterol Guidelines",
    },
    "Heritable Thoracic Aortic Disease": {
        "who": "All first-degree relatives",
        "modalities": [
            "Echocardiography (aortic root)", "Genetic testing",
            "Ophthalmology evaluation (if Marfan suspected)",
        ],
        "initial_screening_age": 0,
        "repeat_interval_years": 2,
        "repeat_until_age": None,
        "notes": (
            "More frequent imaging (annually) if aortic dimension at upper "
            "limit of normal.  Blood pressure management critical."
        ),
        "guideline": "2022 ACC/AHA Aortic Disease Guidelines",
    },
    "Marfan Syndrome": {
        "who": "All first-degree relatives",
        "modalities": [
            "Echocardiography (aortic root)", "Genetic testing (FBN1)",
            "Ophthalmology evaluation", "Skeletal survey",
        ],
        "initial_screening_age": 0,
        "repeat_interval_years": 1,
        "repeat_until_age": None,
        "notes": (
            "Annual aortic root imaging mandatory.  Surgical threshold "
            "is 5.0 cm (or 4.5 cm with risk factors / rapid growth)."
        ),
        "guideline": "2022 ACC/AHA Aortic Disease Guidelines",
    },
    "Loeys-Dietz Syndrome": {
        "who": "All first-degree relatives",
        "modalities": [
            "Echocardiography", "CTA/MRA of entire aorta",
            "Genetic testing (TGFBR1/TGFBR2/SMAD3)",
        ],
        "initial_screening_age": 0,
        "repeat_interval_years": 1,
        "repeat_until_age": None,
        "notes": (
            "Head-to-pelvis vascular imaging recommended.  "
            "Lower surgical threshold (4.0-4.2 cm) for aortic root."
        ),
        "guideline": "2022 ACC/AHA Aortic Disease Guidelines",
    },
    "Cardiac Amyloidosis (ATTR)": {
        "who": "All first-degree relatives if hereditary TTR variant identified",
        "modalities": [
            "Genetic testing for TTR variant", "Echocardiography",
            "Tc-99m PYP/DPD scintigraphy",
        ],
        "initial_screening_age": 30,
        "repeat_interval_years": 3,
        "repeat_until_age": None,
        "notes": (
            "Early detection enables tafamidis initiation before NYHA III.  "
            "Neurological evaluation for polyneuropathy in TTR carriers."
        ),
        "guideline": "2023 ACC Expert Consensus on Cardiac Amyloidosis",
    },
    "Fabry Disease": {
        "who": "All first-degree relatives (X-linked: both male and female)",
        "modalities": [
            "GLA enzyme activity (males)", "Genetic testing (females)",
            "Echocardiography", "Renal function panel",
        ],
        "initial_screening_age": 0,
        "repeat_interval_years": 2,
        "repeat_until_age": None,
        "notes": (
            "Enzyme replacement therapy or oral chaperone should be "
            "considered early in confirmed carriers.  Female heterozygotes "
            "can be significantly affected."
        ),
        "guideline": "2020 ESC Cardiomyopathy Guidelines",
    },
    "Sudden Cardiac Death - Family Screening": {
        "who": "All first-degree relatives of SCD victim",
        "modalities": [
            "12-lead ECG", "Echocardiography", "Exercise stress test",
            "Cardiac MRI (if abnormality detected)", "Genetic testing",
        ],
        "initial_screening_age": 5,
        "repeat_interval_years": 2,
        "repeat_until_age": 50,
        "notes": (
            "Post-mortem genetic testing (molecular autopsy) of the SCD "
            "victim is strongly recommended to guide family screening."
        ),
        "guideline": "2022 ESC Ventricular Arrhythmia Guidelines",
    },
    "Left Ventricular Non-Compaction": {
        "who": "All first-degree relatives",
        "modalities": ["Echocardiography", "Cardiac MRI", "12-lead ECG"],
        "initial_screening_age": 10,
        "repeat_interval_years": 3,
        "repeat_until_age": 60,
        "notes": (
            "LVNC diagnosis can be challenging; cardiac MRI with NC/C "
            "ratio measurement is preferred.  Overlap with DCM and HCM "
            "phenotypes is common."
        ),
        "guideline": "2023 ESC Cardiomyopathy Guidelines",
    },
    "Peripartum Cardiomyopathy": {
        "who": "Female first-degree relatives (screen during and after pregnancy)",
        "modalities": [
            "Echocardiography", "12-lead ECG",
            "BNP/NT-proBNP monitoring",
        ],
        "initial_screening_age": 16,
        "repeat_interval_years": None,
        "repeat_until_age": None,
        "notes": (
            "Screen female first-degree relatives during third trimester and "
            "at 6 weeks postpartum.  Genetic overlap with DCM (TTN truncating "
            "variants in ~15% of PPCM).  Counsel regarding recurrence risk "
            "in subsequent pregnancies if LVEF does not fully recover."
        ),
        "guideline": "ESC 2018 Guidelines on Cardiovascular Disease in Pregnancy",
    },
    "Bicuspid Aortopathy": {
        "who": "All first-degree relatives",
        "modalities": [
            "Echocardiography (aortic valve and ascending aorta)",
            "CTA or MRA of thoracic aorta (if echo abnormal)",
        ],
        "initial_screening_age": 10,
        "repeat_interval_years": 5,
        "repeat_until_age": None,
        "notes": (
            "BAV is present in ~1-2% of general population; heritability "
            "is high (~90%).  Screen for BAV and aortic dilation.  More "
            "frequent imaging (annually) if aorta >4.0 cm.  Surgical "
            "threshold 5.5 cm (or 5.0 cm with risk factors)."
        ),
        "guideline": "2022 ACC/AHA Aortic Disease Guidelines",
    },
}


# ═══════════════════════════════════════════════════════════════════════
# GENE PANEL COSTS (approximate, USD, for patient counseling)
# ═══════════════════════════════════════════════════════════════════════
# These are representative 2025-2026 commercial panel prices used to
# support informed consent discussions.  Actual costs depend on
# laboratory, insurance, and region.

GENE_PANEL_COSTS: Dict[str, Dict[str, Any]] = {
    "HCM_panel": {
        "genes": 11,
        "representative_genes": ["MYH7", "MYBPC3", "TNNT2", "TNNI3"],
        "estimated_cost_usd": 1500,
        "turnaround_days": 21,
        "insurance_coverage": "Typically covered with clinical indication",
        "notes": "Many labs offer comprehensive cardiomyopathy panels at similar cost.",
    },
    "DCM_panel": {
        "genes": 10,
        "representative_genes": ["TTN", "LMNA", "RBM20", "MYH7"],
        "estimated_cost_usd": 1500,
        "turnaround_days": 21,
        "insurance_coverage": "Typically covered with clinical indication",
        "notes": (
            "TTN truncating variants are the most common finding "
            "(~20-25% of familial DCM)."
        ),
    },
    "ACM_panel": {
        "genes": 10,
        "representative_genes": ["PKP2", "DSP", "DSG2", "DSC2"],
        "estimated_cost_usd": 1500,
        "turnaround_days": 21,
        "insurance_coverage": "Typically covered with clinical indication",
        "notes": "PKP2 is the most common gene, accounting for ~40% of positive tests.",
    },
    "channelopathy_panel": {
        "genes": 15,
        "representative_genes": ["KCNQ1", "KCNH2", "SCN5A", "RYR2"],
        "estimated_cost_usd": 1800,
        "turnaround_days": 28,
        "insurance_coverage": "Typically covered with clinical indication",
        "notes": (
            "Combined panel covering LQTS, Brugada, and CPVT.  "
            "Yield is highest when phenotype is well characterised."
        ),
    },
    "FH_panel": {
        "genes": 5,
        "representative_genes": ["LDLR", "PCSK9", "APOB"],
        "estimated_cost_usd": 500,
        "turnaround_days": 14,
        "insurance_coverage": "Covered when LDL >=190 or clinical FH criteria met",
        "notes": "Cascade testing in relatives is highly cost-effective.",
    },
    "aortopathy_panel": {
        "genes": 9,
        "representative_genes": ["FBN1", "TGFBR1", "TGFBR2", "COL3A1"],
        "estimated_cost_usd": 2000,
        "turnaround_days": 28,
        "insurance_coverage": "Typically covered with clinical indication",
        "notes": (
            "FBN1 testing alone may miss Loeys-Dietz and vascular "
            "Ehlers-Danlos.  Comprehensive panel recommended."
        ),
    },
    "TTR_single_gene": {
        "genes": 1,
        "representative_genes": ["TTR"],
        "estimated_cost_usd": 250,
        "turnaround_days": 7,
        "insurance_coverage": "Covered with positive PYP/DPD scan or clinical suspicion",
        "notes": (
            "Val122Ile (V142I) variant present in ~3.4% of African Americans.  "
            "Rapid turnaround supports timely tafamidis initiation."
        ),
    },
    "GLA_single_gene": {
        "genes": 1,
        "representative_genes": ["GLA"],
        "estimated_cost_usd": 300,
        "turnaround_days": 10,
        "insurance_coverage": "Covered with clinical suspicion",
        "notes": (
            "Enzyme activity assay (alpha-galactosidase A) is the first-line "
            "test in males; genetic testing is required to confirm in females."
        ),
    },
    "comprehensive_SCD_panel": {
        "genes": 30,
        "representative_genes": [
            "SCN5A", "KCNQ1", "KCNH2", "MYH7", "MYBPC3",
            "PKP2", "RYR2", "LMNA", "TTN",
        ],
        "estimated_cost_usd": 2500,
        "turnaround_days": 35,
        "insurance_coverage": "Usually covered for SCD family screening",
        "notes": (
            "Broadest panel for unexplained SCD in young individuals.  "
            "Molecular autopsy of decedent tissue is preferred if available."
        ),
    },
    "LVNC_panel": {
        "genes": 6,
        "representative_genes": ["MYH7", "MYBPC3", "TTN", "LMNA"],
        "estimated_cost_usd": 1500,
        "turnaround_days": 21,
        "insurance_coverage": "Typically covered with clinical indication",
        "notes": (
            "Significant genetic overlap with HCM and DCM panels; "
            "some labs offer a combined cardiomyopathy super-panel."
        ),
    },
}


# ═══════════════════════════════════════════════════════════════════════
# PHENOTYPE -> GENE MAP  (reverse index)
# ═══════════════════════════════════════════════════════════════════════
# Built automatically from GENOMIC_TRIGGER_MAP so it stays in sync.
# Maps each gene to the set of phenotypic findings it is associated
# with, enabling a gene-first lookup (e.g., "patient has LMNA variant
# -- what phenotypes should we watch for?").


def _build_phenotype_gene_map() -> Dict[str, List[str]]:
    """Build a reverse mapping from gene symbol to associated phenotypes.

    Returns a dict where keys are gene symbols (e.g. ``"MYH7"``) and
    values are sorted lists of phenotype/finding keys from the trigger
    map (e.g. ``["unexplained_dcm", "unexplained_lvh"]``).
    """
    gene_to_phenotypes: Dict[str, set] = defaultdict(set)
    for finding_key, entry in GENOMIC_TRIGGER_MAP.items():
        for gene in entry.get("gene_panel", []):
            gene_to_phenotypes[gene].add(finding_key)
    return {
        gene: sorted(phenotypes)
        for gene, phenotypes in sorted(gene_to_phenotypes.items())
    }


PHENOTYPE_GENE_MAP: Dict[str, List[str]] = _build_phenotype_gene_map()


# ═══════════════════════════════════════════════════════════════════════
# INTERNAL TYPES
# ═══════════════════════════════════════════════════════════════════════
# CrossModalTrigger from models.py is lean (trigger_source, finding,
# gene_panel, conditions, rationale).  Internally we track additional
# metadata in _TriggerContext to drive cascade logic, reporting, and
# genomic query construction without bloating the Pydantic model.


class _TriggerContext:
    """Internal bookkeeping attached to each trigger during evaluation.

    Not serialised or returned to callers.  Used by the engine to track
    urgency, criteria, guideline, cascade state, and confidence.
    """

    __slots__ = (
        "finding_key", "urgency", "trigger_type", "criteria",
        "guideline", "confidence", "cascade_children", "extra",
    )

    def __init__(
        self,
        finding_key: str,
        urgency: SeverityLevel = SeverityLevel.MODERATE,
        trigger_type: str = "imaging",
        criteria: str = "",
        guideline: str = "",
        confidence: float = 1.0,
        extra: str = "",
    ) -> None:
        self.finding_key = finding_key
        self.urgency = urgency
        self.trigger_type = trigger_type
        self.criteria = criteria
        self.guideline = guideline
        self.confidence = confidence
        self.cascade_children: List[str] = []
        self.extra = extra


# ═══════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════


def _make_trigger_id(finding_key: str, extra: str = "") -> str:
    """Generate a deterministic short ID for a trigger instance.

    Combines the finding key with optional extra context (e.g.,
    measurement value) to produce a hex digest suitable for
    deduplication.
    """
    raw = f"{finding_key}:{extra}:{time.time()}"
    return hashlib.md5(raw.encode()).hexdigest()[:12]


def _severity_from_urgency(urgency_str: str) -> SeverityLevel:
    """Convert a trigger-map urgency string to a SeverityLevel enum."""
    return _URGENCY_MAP.get(urgency_str, SeverityLevel.MODERATE)


def _build_rationale(
    entry: Dict[str, Any],
    trigger_type: str,
    extra: str,
    confidence: float,
) -> str:
    """Compose the rationale string for a CrossModalTrigger.

    Packs the criteria, guideline, urgency, trigger type, and confidence
    into a single human-readable string so that downstream consumers
    (and the LLM synthesis layer) have full context.
    """
    parts: List[str] = []
    criteria = entry.get("criteria", "")
    if criteria:
        parts.append(f"Criteria: {criteria}")
    guideline = entry.get("guideline", "")
    if guideline:
        parts.append(f"Guideline: {guideline}")
    urgency = entry.get("urgency", "moderate")
    parts.append(f"Urgency: {urgency}")
    parts.append(f"Source: {trigger_type}")
    if extra:
        parts.append(f"Detail: {extra}")
    if confidence < 1.0:
        parts.append(f"Confidence: {confidence:.0%}")
    return " | ".join(parts)


def prioritize_triggers(
    triggers: List[Tuple[CrossModalTrigger, _TriggerContext]],
) -> List[Tuple[CrossModalTrigger, _TriggerContext]]:
    """Sort trigger pairs by clinical urgency (most urgent first).

    Ordering:
        1. critical  (immediate action -- channelopathies, SCD risk)
        2. high      (urgent workup within days)
        3. moderate  (workup within weeks)
        4. low       (elective / informational)

    Within the same urgency tier, triggers are stable-sorted in
    discovery order.

    Parameters
    ----------
    triggers:
        Unsorted list of ``(CrossModalTrigger, _TriggerContext)`` pairs.

    Returns
    -------
    A new list sorted in descending urgency order.
    """
    return sorted(
        triggers,
        key=lambda pair: _URGENCY_PRIORITY.get(pair[1].urgency, 99),
    )


def prioritize_trigger_models(
    triggers: List[CrossModalTrigger],
) -> List[CrossModalTrigger]:
    """Sort standalone CrossModalTrigger models by urgency extracted from rationale.

    This is a convenience for callers who only have the Pydantic models
    (without the internal ``_TriggerContext``).  It parses the urgency
    from the rationale string.
    """
    def _extract_urgency_rank(t: CrossModalTrigger) -> int:
        rationale_lower = t.rationale.lower()
        if "urgency: critical" in rationale_lower:
            return 0
        if "urgency: high" in rationale_lower:
            return 1
        if "urgency: moderate" in rationale_lower:
            return 2
        if "urgency: low" in rationale_lower:
            return 3
        return 4

    return sorted(triggers, key=_extract_urgency_rank)


def genes_for_phenotype(finding_key: str) -> List[str]:
    """Return the gene panel for a given clinical finding key.

    Convenience wrapper around GENOMIC_TRIGGER_MAP for quick lookups.
    Returns an empty list if the finding key is not recognised.
    """
    entry = GENOMIC_TRIGGER_MAP.get(finding_key, {})
    return list(entry.get("gene_panel", []))


def phenotypes_for_gene(gene: str) -> List[str]:
    """Return all clinical finding keys associated with a gene.

    Uses the pre-built PHENOTYPE_GENE_MAP reverse index.
    """
    return list(PHENOTYPE_GENE_MAP.get(gene, []))


def estimate_panel_cost(trigger: CrossModalTrigger) -> Optional[Dict[str, Any]]:
    """Look up the approximate panel cost for a given trigger.

    Matches the trigger's gene panel against GENE_PANEL_COSTS entries
    by checking for overlap with the representative genes.  Returns
    the best-matching panel cost entry, or None if no match is found.
    """
    trigger_genes = set(trigger.gene_panel)
    best_match: Optional[str] = None
    best_overlap = 0

    for panel_name, panel_info in GENE_PANEL_COSTS.items():
        rep_genes = set(panel_info.get("representative_genes", []))
        overlap = len(trigger_genes & rep_genes)
        if overlap > best_overlap:
            best_overlap = overlap
            best_match = panel_name

    if best_match is not None:
        return {"panel_name": best_match, **GENE_PANEL_COSTS[best_match]}
    return None


# ═══════════════════════════════════════════════════════════════════════
# CROSS-MODAL ENGINE
# ═══════════════════════════════════════════════════════════════════════


class CrossModalEngine:
    """Engine that evaluates clinical findings against genomic trigger rules.

    The engine is stateless with respect to patient data: each call to
    ``evaluate_triggers``, ``check_imaging_triggers``, etc. takes the
    relevant clinical observations as arguments and returns a list of
    ``CrossModalTrigger`` instances.  The caller (typically the
    CardioAgent orchestration layer) is responsible for aggregating
    triggers across modalities and forwarding them to the genomic
    query pipeline.

    Attributes
    ----------
    trigger_map : dict
        Reference to ``GENOMIC_TRIGGER_MAP``.
    cascade_enabled : bool
        When True (default), ``get_cascade_triggers`` is invoked
        automatically at the end of ``evaluate_triggers``.
    """

    def __init__(self, *, cascade_enabled: bool = True) -> None:
        self.trigger_map: Dict[str, Dict[str, Any]] = GENOMIC_TRIGGER_MAP
        self.cascade_enabled = cascade_enabled
        # Internal context store, keyed by finding key, populated
        # during a single evaluate_triggers call.
        self._contexts: Dict[str, _TriggerContext] = {}
        logger.info(
            "CrossModalEngine initialised with %d trigger rules (cascade=%s)",
            len(self.trigger_map),
            self.cascade_enabled,
        )

    # ── Primary evaluation entry point ────────────────────────────────

    def evaluate_triggers(
        self,
        findings: Dict[str, Any],
    ) -> List[CrossModalTrigger]:
        """Evaluate a bundle of clinical findings against all trigger rules.

        ``findings`` is a flexible dict whose keys may include:

        Imaging fields:
            ``modality`` (str), ``measurements`` (dict),
            ``imaging_findings`` (list[str])

        ECG fields:
            ``qtc`` (int, ms), ``rhythm`` (str), ``morphology`` (str)

        Clinical fields:
            ``age`` (int), ``sex`` (str), ``ldl`` (float, mg/dL),
            ``conditions`` (list[str])

        Family history:
            ``family_hx`` (dict) -- keys are condition names, values
            are booleans or detail strings.

        Returns a deduplicated, priority-sorted list of triggers.
        """
        pairs: List[Tuple[CrossModalTrigger, _TriggerContext]] = []
        seen_keys: set = set()
        self._contexts = {}

        # 1. Imaging triggers
        if any(k in findings for k in ("modality", "measurements", "imaging_findings")):
            for trigger, ctx in self._check_imaging_internal(
                modality=findings.get("modality", ""),
                measurements=findings.get("measurements", {}),
                findings_list=findings.get("imaging_findings", []),
            ):
                if ctx.finding_key not in seen_keys:
                    pairs.append((trigger, ctx))
                    seen_keys.add(ctx.finding_key)
                    self._contexts[ctx.finding_key] = ctx

        # 2. ECG triggers
        if any(k in findings for k in ("qtc", "rhythm", "morphology")):
            for trigger, ctx in self._check_ecg_internal(
                qtc=findings.get("qtc", 0),
                rhythm=findings.get("rhythm", ""),
                morphology=findings.get("morphology", ""),
                family_hx=findings.get("family_hx", {}),
            ):
                if ctx.finding_key not in seen_keys:
                    pairs.append((trigger, ctx))
                    seen_keys.add(ctx.finding_key)
                    self._contexts[ctx.finding_key] = ctx

        # 3. Clinical triggers
        if any(k in findings for k in ("age", "ldl", "conditions")):
            for trigger, ctx in self._check_clinical_internal(
                age=findings.get("age", 0),
                sex=findings.get("sex", "unknown"),
                ldl=findings.get("ldl", 0.0),
                family_hx=findings.get("family_hx", {}),
                conditions=findings.get("conditions", []),
            ):
                if ctx.finding_key not in seen_keys:
                    pairs.append((trigger, ctx))
                    seen_keys.add(ctx.finding_key)
                    self._contexts[ctx.finding_key] = ctx

        # 4. Cascade triggers
        if self.cascade_enabled and pairs:
            for trigger, ctx in self._get_cascade_internal(pairs):
                if ctx.finding_key not in seen_keys:
                    pairs.append((trigger, ctx))
                    seen_keys.add(ctx.finding_key)
                    self._contexts[ctx.finding_key] = ctx

        # 5. Sort by urgency
        pairs = prioritize_triggers(pairs)

        result = [t for t, _ in pairs]
        logger.info(
            "evaluate_triggers produced %d triggers from %d finding keys",
            len(result),
            len(findings),
        )
        return result

    # ── Public convenience methods ────────────────────────────────────

    def check_imaging_triggers(
        self,
        modality: str,
        measurements: Dict[str, Any],
        findings: List[str],
    ) -> List[CrossModalTrigger]:
        """Check imaging findings for genomic triggers.

        Evaluates echocardiographic, CT, and MRI measurements against
        the thresholds encoded in the trigger map.

        Parameters
        ----------
        modality:
            Imaging modality (see ``ImagingModality`` enum values).
        measurements:
            Quantitative measurements, e.g.
            ``{"lv_wall_thickness_mm": 16, "lvef_percent": 38,
              "aortic_root_cm": 4.2, "nc_c_ratio": 2.5}``.
        findings:
            Free-text or coded findings list, e.g.
            ``["fibro-fatty replacement", "mid-wall LGE"]``.

        Returns
        -------
        List of CrossModalTrigger instances for any matched rules.
        """
        return [
            t for t, _ in self._check_imaging_internal(
                modality, measurements, findings,
            )
        ]

    def check_ecg_triggers(
        self,
        qtc: int,
        rhythm: str,
        morphology: str,
        family_hx: Dict[str, Any],
    ) -> List[CrossModalTrigger]:
        """Check ECG findings for channelopathy and arrhythmia triggers.

        Parameters
        ----------
        qtc:
            Corrected QT interval in milliseconds.
        rhythm:
            Rhythm description (e.g. ``"sinus"``, ``"polymorphic_vt"``).
        morphology:
            ECG morphology description (e.g. ``"brugada_type1"``).
        family_hx:
            Family history dict.  Relevant keys include
            ``"scd"`` (sudden cardiac death in relative < 40),
            ``"long_qt"``, ``"brugada"``.

        Returns
        -------
        List of CrossModalTrigger instances for any matched rules.
        """
        return [
            t for t, _ in self._check_ecg_internal(
                qtc, rhythm, morphology, family_hx,
            )
        ]

    def check_clinical_triggers(
        self,
        age: int,
        ldl: float = 0.0,
        family_hx: Optional[Dict[str, Any]] = None,
        conditions: Optional[List[str]] = None,
        sex: str = "unknown",
    ) -> List[CrossModalTrigger]:
        """Check clinical data for genetic testing triggers.

        Parameters
        ----------
        age:
            Patient age in years.
        ldl:
            LDL cholesterol in mg/dL.
        family_hx:
            Family history dict.
        conditions:
            List of active problem-list entries or ICD codes.
        sex:
            ``"male"`` or ``"female"``.

        Returns
        -------
        List of CrossModalTrigger instances for any matched rules.
        """
        return [
            t for t, _ in self._check_clinical_internal(
                age, sex, ldl, family_hx or {}, conditions or [],
            )
        ]

    def build_genomic_query(
        self,
        trigger: CrossModalTrigger,
    ) -> Dict[str, Any]:
        """Build a query payload for the genomic_evidence Milvus collection.

        The returned dict is structured for use with the HCLS AI Factory
        RAG pipeline.  It includes:

        - ``collection``: target collection name (``genomic_evidence``).
        - ``query_texts``: list of natural-language query strings that
          will be embedded via BGE-small-en-v1.5.
        - ``gene_filter``: Milvus expression to restrict results to
          relevant genes.
        - ``top_k``: number of results per query.
        - ``score_threshold``: minimum similarity score.
        - ``metadata``: additional context passed through to the result.

        Parameters
        ----------
        trigger:
            The ``CrossModalTrigger`` to build a query for.

        Returns
        -------
        Query dict ready for submission to the RAG pipeline.
        """
        genes = trigger.gene_panel
        conditions = trigger.conditions

        # Build natural-language query strings
        query_texts: List[str] = []
        for condition in conditions:
            query_texts.append(
                f"{condition} genetic variants pathogenic likely pathogenic"
            )
        for gene in genes[:5]:  # Limit to top 5 genes to control query volume
            query_texts.append(
                f"{gene} gene cardiac variant clinical significance"
            )
        if not query_texts:
            query_texts.append(
                f"cardiac genetic testing {trigger.finding}"
            )

        # Build Milvus filter expression for gene symbols
        if genes:
            gene_list_str = ", ".join(f'"{g}"' for g in genes)
            gene_filter = f"gene_symbol in [{gene_list_str}]"
        else:
            gene_filter = ""

        return {
            "collection": GENOMIC_COLLECTION,
            "query_texts": query_texts,
            "gene_filter": gene_filter,
            "top_k": TOP_K_PER_QUERY,
            "score_threshold": SCORE_THRESHOLD,
            "embedding_dim": EMBEDDING_DIM,
            "metadata": {
                "trigger_source": trigger.trigger_source,
                "finding": trigger.finding,
                "conditions": conditions,
            },
        }

    def format_trigger_report(
        self,
        triggers: List[CrossModalTrigger],
    ) -> str:
        """Format a list of triggers into a human-readable clinical report.

        The report is structured for inclusion in a clinical decision
        support note.  Triggers are grouped by urgency, with critical
        items highlighted.

        Parameters
        ----------
        triggers:
            List of triggers (will be re-sorted internally).

        Returns
        -------
        Multi-line string with formatted trigger report.
        """
        if not triggers:
            return "No cross-modal genomic triggers identified."

        sorted_triggers = prioritize_trigger_models(triggers)
        lines: List[str] = [
            "=" * 72,
            "CROSS-MODAL GENOMIC TRIGGER REPORT",
            f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
            f"Total triggers: {len(sorted_triggers)}",
            "=" * 72,
            "",
        ]

        # Group by urgency (parsed from rationale)
        urgency_groups: Dict[str, List[CrossModalTrigger]] = defaultdict(list)
        for t in sorted_triggers:
            rationale_lower = t.rationale.lower()
            if "urgency: critical" in rationale_lower:
                urgency_groups["critical"].append(t)
            elif "urgency: high" in rationale_lower:
                urgency_groups["high"].append(t)
            elif "urgency: moderate" in rationale_lower:
                urgency_groups["moderate"].append(t)
            elif "urgency: low" in rationale_lower:
                urgency_groups["low"].append(t)
            else:
                urgency_groups["moderate"].append(t)

        urgency_order = ["critical", "high", "moderate", "low"]
        urgency_labels = {
            "critical": "CRITICAL -- Immediate Action Required",
            "high": "HIGH PRIORITY -- Urgent Workup Recommended",
            "moderate": "MODERATE -- Workup Within Weeks",
            "low": "LOW -- Elective / Informational",
        }

        for urgency_key in urgency_order:
            group = urgency_groups.get(urgency_key, [])
            if not group:
                continue

            label = urgency_labels.get(urgency_key, urgency_key.upper())
            lines.append(f"--- {label} ---")
            lines.append("")

            for i, trigger in enumerate(group, 1):
                finding_display = trigger.finding.replace("_", " ").title()
                lines.append(f"  {i}. {finding_display}")
                lines.append(f"     Suspected: {', '.join(trigger.conditions)}")
                lines.append(f"     Gene panel: {', '.join(trigger.gene_panel)}")
                lines.append(f"     Source: {trigger.trigger_source}")

                # Parse criteria and guideline from rationale
                for segment in trigger.rationale.split(" | "):
                    if segment.startswith("Criteria:"):
                        lines.append(f"     {segment}")
                    elif segment.startswith("Guideline:"):
                        lines.append(f"     {segment}")

                # Panel cost estimate
                cost_info = estimate_panel_cost(trigger)
                if cost_info:
                    lines.append(
                        f"     Est. cost: ${cost_info['estimated_cost_usd']:,} "
                        f"({cost_info['turnaround_days']}d turnaround)"
                    )

                # Family screening
                for condition in trigger.conditions:
                    protocol = FAMILY_SCREENING_PROTOCOLS.get(condition)
                    if protocol:
                        lines.append(
                            f"     Family screening: {protocol['who']}, "
                            f"starting age {protocol['initial_screening_age']}"
                        )
                        break

                lines.append("")

        # Summary footer
        gene_set: set = set()
        for t in sorted_triggers:
            gene_set.update(t.gene_panel)

        lines.append("-" * 72)
        lines.append(f"Total unique genes across all panels: {len(gene_set)}")
        lines.append(f"Genes: {', '.join(sorted(gene_set))}")
        lines.append("")
        lines.append(
            "Note: Genetic testing recommendations are based on clinical "
            "criteria and published guidelines.  Final testing decisions "
            "should incorporate genetic counseling and patient preferences."
        )
        lines.append("=" * 72)

        return "\n".join(lines)

    def get_cascade_triggers(
        self,
        initial_triggers: List[CrossModalTrigger],
    ) -> List[CrossModalTrigger]:
        """Check if initial triggers cascade into additional testing.

        Public wrapper that works with standalone CrossModalTrigger
        instances (without internal context).  Reconstructs contexts
        from the trigger's finding key and rationale.

        Parameters
        ----------
        initial_triggers:
            Triggers already identified by the primary evaluation.

        Returns
        -------
        Additional triggers generated by cascade logic (may be empty).
        """
        pairs: List[Tuple[CrossModalTrigger, _TriggerContext]] = []
        for t in initial_triggers:
            ctx = _TriggerContext(
                finding_key=t.finding,
                urgency=_severity_from_urgency(
                    self._parse_urgency_from_rationale(t.rationale)
                ),
                trigger_type=t.trigger_source,
            )
            pairs.append((t, ctx))

        cascade_pairs = self._get_cascade_internal(pairs)
        return [t for t, _ in cascade_pairs]

    # ── Batch evaluation ──────────────────────────────────────────────

    def evaluate_batch(
        self,
        patients: List[Dict[str, Any]],
    ) -> Dict[str, List[CrossModalTrigger]]:
        """Evaluate triggers for a batch of patients.

        Parameters
        ----------
        patients:
            List of dicts, each with a ``"patient_id"`` key and the
            same finding keys accepted by ``evaluate_triggers``.

        Returns
        -------
        Dict mapping patient_id to their trigger list.
        """
        results: Dict[str, List[CrossModalTrigger]] = {}
        for patient in patients:
            pid = patient.get("patient_id", f"unknown_{id(patient)}")
            try:
                triggers = self.evaluate_triggers(patient)
                results[pid] = triggers
            except Exception:
                logger.exception("Error evaluating triggers for patient %s", pid)
                results[pid] = []
        return results

    # ── Summary statistics ────────────────────────────────────────────

    def summarize_triggers(
        self,
        triggers: List[CrossModalTrigger],
    ) -> Dict[str, Any]:
        """Generate summary statistics for a set of triggers.

        Useful for dashboard displays and aggregate reporting.

        Returns a dict with keys: ``total``, ``by_urgency``,
        ``by_source``, ``unique_genes``, ``unique_conditions``,
        ``estimated_total_cost``.
        """
        by_urgency: Dict[str, int] = defaultdict(int)
        by_source: Dict[str, int] = defaultdict(int)
        all_genes: set = set()
        all_conditions: set = set()
        total_cost = 0

        for t in triggers:
            # Parse urgency from rationale
            urgency_str = self._parse_urgency_from_rationale(t.rationale)
            by_urgency[urgency_str] += 1
            by_source[t.trigger_source] += 1
            all_genes.update(t.gene_panel)
            all_conditions.update(t.conditions)
            cost_info = estimate_panel_cost(t)
            if cost_info:
                total_cost += cost_info["estimated_cost_usd"]

        return {
            "total": len(triggers),
            "by_urgency": dict(by_urgency),
            "by_source": dict(by_source),
            "unique_genes": sorted(all_genes),
            "unique_conditions": sorted(all_conditions),
            "estimated_total_cost_usd": total_cost,
        }

    # ═══════════════════════════════════════════════════════════════════
    # INTERNAL METHODS
    # ═══════════════════════════════════════════════════════════════════

    def _check_imaging_internal(
        self,
        modality: str,
        measurements: Dict[str, Any],
        findings_list: List[str],
    ) -> List[Tuple[CrossModalTrigger, _TriggerContext]]:
        """Internal imaging trigger evaluation returning (trigger, context) pairs."""
        pairs: List[Tuple[CrossModalTrigger, _TriggerContext]] = []

        lv_wall = measurements.get("lv_wall_thickness_mm", 0)
        lvef = measurements.get("lvef_percent", 100)
        lv_dilated = measurements.get("lv_dilated", False)
        aortic_root = measurements.get("aortic_root_cm", 0)
        nc_c_ratio = measurements.get("nc_c_ratio", 0)
        rv_dilated = measurements.get("rv_dilated", False)
        rv_dysfunction = measurements.get("rv_dysfunction", False)
        pyp_positive = measurements.get("pyp_dpd_positive", False)
        low_voltage = measurements.get("low_voltage_ecg", False)
        age = measurements.get("patient_age", 999)
        family_hx = measurements.get("family_hx_cm", False)
        diastolic_dysfunction = measurements.get("diastolic_dysfunction", False)
        renal_dysfunction = measurements.get("renal_dysfunction", False)
        neuropathy = measurements.get("neuropathy", False)

        findings_lower = [f.lower() for f in findings_list]

        # ── Unexplained LVH ──
        threshold = 13 if family_hx else 15
        if lv_wall >= threshold:
            pairs.append(self._create_trigger_pair(
                "unexplained_lvh",
                trigger_type="imaging",
                extra=f"wall={lv_wall}mm",
            ))

        # ── Unexplained DCM ──
        if lvef < 45 and lv_dilated and age < 60:
            pairs.append(self._create_trigger_pair(
                "unexplained_dcm",
                trigger_type="imaging",
                extra=f"ef={lvef}%",
            ))

        # ── Arrhythmogenic CM ──
        has_vent_arrhythmia = any(
            term in " ".join(findings_lower)
            for term in [
                "ventricular arrhythmia", "ventricular tachycardia", "vt",
            ]
        )
        has_fibro_fatty = any(
            "fibro-fatty" in f or "fibro fatty" in f for f in findings_lower
        )
        if (rv_dilated or rv_dysfunction) and (has_vent_arrhythmia or has_fibro_fatty):
            pairs.append(self._create_trigger_pair(
                "arrhythmogenic_cm",
                trigger_type="imaging",
                extra="rv_abnl+arrhythmia",
            ))

        # ── LVNC ──
        if nc_c_ratio > 2.3:
            pairs.append(self._create_trigger_pair(
                "lvnc",
                trigger_type="imaging",
                extra=f"nc_c={nc_c_ratio}",
            ))

        # ── Aortic dilation ──
        if aortic_root >= 4.0 and age < 50:
            pairs.append(self._create_trigger_pair(
                "aortic_dilation",
                trigger_type="imaging",
                extra=f"root={aortic_root}cm",
            ))

        # ── Cardiac amyloid ──
        if pyp_positive or (lv_wall >= 13 and diastolic_dysfunction and low_voltage):
            pairs.append(self._create_trigger_pair(
                "cardiac_amyloid",
                trigger_type="imaging",
                extra="pyp" if pyp_positive else "lvh+dd+lv",
            ))

        # ── Fabry disease ──
        if lv_wall >= 13 and age < 50 and (renal_dysfunction or neuropathy):
            pairs.append(self._create_trigger_pair(
                "fabry_disease",
                trigger_type="imaging",
                extra=f"wall={lv_wall}mm,age={age}",
            ))

        logger.debug(
            "check_imaging_triggers: modality=%s -> %d triggers",
            modality, len(pairs),
        )
        return pairs

    def _check_ecg_internal(
        self,
        qtc: int,
        rhythm: str,
        morphology: str,
        family_hx: Dict[str, Any],
    ) -> List[Tuple[CrossModalTrigger, _TriggerContext]]:
        """Internal ECG trigger evaluation returning (trigger, context) pairs."""
        pairs: List[Tuple[CrossModalTrigger, _TriggerContext]] = []
        has_syncope = family_hx.get("syncope", False)
        has_scd_hx = family_hx.get("scd", False)
        morphology_lower = morphology.lower() if morphology else ""
        rhythm_lower = rhythm.lower() if rhythm else ""

        # ── Long QT ──
        if qtc > 480 or (qtc > 460 and (has_syncope or has_scd_hx)):
            pairs.append(self._create_trigger_pair(
                "long_qt",
                trigger_type="ecg",
                extra=f"qtc={qtc}ms",
            ))

        # ── Brugada ──
        brugada_keywords = [
            "brugada", "coved_st", "coved st", "type_1", "type 1",
        ]
        if any(kw in morphology_lower for kw in brugada_keywords):
            pairs.append(self._create_trigger_pair(
                "brugada_pattern",
                trigger_type="ecg",
                extra=f"morph={morphology}",
            ))

        # ── CPVT ──
        cpvt_rhythms = [
            "bidirectional_vt", "bidirectional vt",
            "polymorphic_vt", "polymorphic vt",
            "catecholaminergic",
        ]
        if any(r in rhythm_lower for r in cpvt_rhythms):
            pairs.append(self._create_trigger_pair(
                "cpvt_suspected",
                trigger_type="ecg",
                extra=f"rhythm={rhythm}",
            ))

        # ── Family history SCD (ECG context) ──
        if has_scd_hx:
            existing_findings = {ctx.finding_key for _, ctx in pairs}
            if "scd_family_history" not in existing_findings:
                pairs.append(self._create_trigger_pair(
                    "scd_family_history",
                    trigger_type="ecg",
                    extra="ecg_screen",
                ))

        logger.debug(
            "check_ecg_triggers: qtc=%d, rhythm=%s -> %d triggers",
            qtc, rhythm, len(pairs),
        )
        return pairs

    def _check_clinical_internal(
        self,
        age: int,
        sex: str,
        ldl: float,
        family_hx: Dict[str, Any],
        conditions: List[str],
    ) -> List[Tuple[CrossModalTrigger, _TriggerContext]]:
        """Internal clinical trigger evaluation returning (trigger, context) pairs."""
        pairs: List[Tuple[CrossModalTrigger, _TriggerContext]] = []
        conditions_lower = [c.lower() for c in conditions]
        sex_lower = sex.lower() if sex else "unknown"

        # ── Premature CAD / Familial Hypercholesterolemia ──
        premature_cad_threshold = 55 if sex_lower == "male" else 65
        has_cad = any(
            term in " ".join(conditions_lower)
            for term in [
                "cad", "coronary artery disease",
                "myocardial infarction", "mi", "stemi", "nstemi",
            ]
        )
        if (has_cad and age < premature_cad_threshold) or ldl >= 190:
            pairs.append(self._create_trigger_pair(
                "premature_cad",
                trigger_type="clinical",
                extra=f"age={age},ldl={ldl}",
            ))

        # ── Aortic dilation (clinical suspicion) ──
        marfan_features = any(
            term in " ".join(conditions_lower)
            for term in [
                "marfan", "aortic dilation", "aortic aneurysm",
                "loeys-dietz", "ehlers-danlos",
            ]
        )
        if marfan_features and age < 50:
            pairs.append(self._create_trigger_pair(
                "aortic_dilation",
                trigger_type="clinical",
                extra=f"clinical_suspicion,age={age}",
            ))

        # ── SCD family history ──
        if family_hx.get("scd", False) or family_hx.get("sudden_death", False):
            scd_age = family_hx.get(
                "scd_age", family_hx.get("sudden_death_age", 999)
            )
            if scd_age < 40:
                pairs.append(self._create_trigger_pair(
                    "scd_family_history",
                    trigger_type="clinical",
                    extra=f"scd_age={scd_age}",
                ))

        # ── Family history of cardiomyopathy ──
        if family_hx.get("cardiomyopathy", False):
            cm_type = family_hx.get("cardiomyopathy_type", "unspecified")
            mapping = {
                "hcm": "unexplained_lvh",
                "dcm": "unexplained_dcm",
                "acm": "arrhythmogenic_cm",
                "arvc": "arrhythmogenic_cm",
                "lvnc": "lvnc",
            }
            finding_key = mapping.get(cm_type.lower(), "scd_family_history")
            pairs.append(self._create_trigger_pair(
                finding_key,
                trigger_type="family_history",
                extra=f"fhx_cm={cm_type}",
            ))

        # ── Family history of channelopathy ──
        if family_hx.get("long_qt", False):
            pairs.append(self._create_trigger_pair(
                "long_qt",
                trigger_type="family_history",
                extra="fhx_lqts",
            ))
        if family_hx.get("brugada", False):
            pairs.append(self._create_trigger_pair(
                "brugada_pattern",
                trigger_type="family_history",
                extra="fhx_brugada",
            ))

        # ── Fabry disease (clinical suspicion) ──
        fabry_features = any(
            term in " ".join(conditions_lower)
            for term in ["fabry", "alpha-galactosidase", "angiokeratoma"]
        )
        if fabry_features:
            pairs.append(self._create_trigger_pair(
                "fabry_disease",
                trigger_type="clinical",
                extra="clinical_suspicion",
            ))

        # ── Cardiac amyloidosis (clinical suspicion) ──
        amyloid_features = any(
            term in " ".join(conditions_lower)
            for term in ["amyloid", "attr", "transthyretin", "al amyloid"]
        )
        if amyloid_features:
            pairs.append(self._create_trigger_pair(
                "cardiac_amyloid",
                trigger_type="clinical",
                extra="clinical_suspicion",
            ))

        logger.debug(
            "check_clinical_triggers: age=%d, sex=%s, ldl=%.1f -> %d triggers",
            age, sex_lower, ldl, len(pairs),
        )
        return pairs

    def _get_cascade_internal(
        self,
        initial_pairs: List[Tuple[CrossModalTrigger, _TriggerContext]],
    ) -> List[Tuple[CrossModalTrigger, _TriggerContext]]:
        """Internal cascade logic returning (trigger, context) pairs.

        Cascade rules encode clinical logic where finding one condition
        warrants screening for another:

        - HCM (unexplained_lvh) -> also screen for Fabry (GLA) and
          cardiac amyloid (TTR) if not already triggered.
        - DCM with conduction disease (LMNA) -> screen for ACM
          (desmosomal genes overlap).
        - Any channelopathy -> trigger SCD family screening.
        - LVNC -> check HCM and DCM (phenotypic overlap).
        - ACM -> also screen LMNA / DCM overlap.
        """
        cascade: List[Tuple[CrossModalTrigger, _TriggerContext]] = []
        existing_keys = {ctx.finding_key for _, ctx in initial_pairs}

        for _, ctx in initial_pairs:
            finding = ctx.finding_key

            # HCM -> also check Fabry and ATTR amyloid
            if finding == "unexplained_lvh":
                if "fabry_disease" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "fabry_disease",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.6,
                    )
                    cascade.append(pair)
                    existing_keys.add("fabry_disease")

                if "cardiac_amyloid" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "cardiac_amyloid",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.5,
                    )
                    cascade.append(pair)
                    existing_keys.add("cardiac_amyloid")

            # DCM -> also check ACM (desmosomal overlap)
            if finding == "unexplained_dcm":
                if "arrhythmogenic_cm" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "arrhythmogenic_cm",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.5,
                    )
                    cascade.append(pair)
                    existing_keys.add("arrhythmogenic_cm")

            # Any channelopathy -> SCD family screening
            if finding in ("long_qt", "brugada_pattern", "cpvt_suspected"):
                if "scd_family_history" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "scd_family_history",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.8,
                    )
                    cascade.append(pair)
                    existing_keys.add("scd_family_history")

            # ACM -> also screen LMNA (overlap with DCM+conduction disease)
            if finding == "arrhythmogenic_cm":
                if "unexplained_dcm" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "unexplained_dcm",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.4,
                    )
                    cascade.append(pair)
                    existing_keys.add("unexplained_dcm")

            # LVNC -> also check HCM and DCM (phenotypic overlap)
            if finding == "lvnc":
                if "unexplained_lvh" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "unexplained_lvh",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.4,
                    )
                    cascade.append(pair)
                    existing_keys.add("unexplained_lvh")

                if "unexplained_dcm" not in existing_keys:
                    pair = self._create_trigger_pair(
                        "unexplained_dcm",
                        trigger_type="cascade",
                        extra=f"from_{finding}",
                        confidence=0.4,
                    )
                    cascade.append(pair)
                    existing_keys.add("unexplained_dcm")

        if cascade:
            logger.info(
                "Cascade logic generated %d additional triggers from %d initial",
                len(cascade),
                len(initial_pairs),
            )

        return cascade

    # ── Trigger construction ──────────────────────────────────────────

    def _create_trigger_pair(
        self,
        finding_key: str,
        *,
        trigger_type: str = "imaging",
        extra: str = "",
        confidence: float = 1.0,
    ) -> Tuple[CrossModalTrigger, _TriggerContext]:
        """Create a (CrossModalTrigger, _TriggerContext) pair from a trigger map entry.

        Parameters
        ----------
        finding_key:
            Key into ``GENOMIC_TRIGGER_MAP`` (e.g. ``"unexplained_lvh"``).
        trigger_type:
            Source category (``"imaging"``, ``"ecg"``, ``"clinical"``,
            ``"family_history"``, ``"cascade"``).
        extra:
            Additional context string appended to the trigger ID.
        confidence:
            Confidence level (0.0 -- 1.0).  Cascade triggers typically
            have reduced confidence.

        Returns
        -------
        A ``(CrossModalTrigger, _TriggerContext)`` tuple.
        """
        entry = self.trigger_map.get(finding_key, {})
        urgency_str = entry.get("urgency", "moderate")

        rationale = _build_rationale(entry, trigger_type, extra, confidence)

        trigger = CrossModalTrigger(
            trigger_source=trigger_type,
            finding=finding_key,
            gene_panel=list(entry.get("gene_panel", [])),
            conditions=list(entry.get("conditions", [])),
            rationale=rationale,
        )

        ctx = _TriggerContext(
            finding_key=finding_key,
            urgency=_severity_from_urgency(urgency_str),
            trigger_type=trigger_type,
            criteria=entry.get("criteria", ""),
            guideline=entry.get("guideline", ""),
            confidence=confidence,
            extra=extra,
        )

        return trigger, ctx

    @staticmethod
    def _parse_urgency_from_rationale(rationale: str) -> str:
        """Extract the urgency string from a rationale field."""
        for segment in rationale.split(" | "):
            if segment.lower().startswith("urgency:"):
                return segment.split(":", 1)[1].strip().lower()
        return "moderate"


# ═══════════════════════════════════════════════════════════════════════
# MODULE-LEVEL CONVENIENCE API
# ═══════════════════════════════════════════════════════════════════════

# Singleton engine for simple use-cases (import and call directly
# without explicit instantiation).
_default_engine: Optional[CrossModalEngine] = None


def get_engine() -> CrossModalEngine:
    """Return (and lazily create) the module-level CrossModalEngine."""
    global _default_engine
    if _default_engine is None:
        _default_engine = CrossModalEngine()
    return _default_engine


def evaluate(findings: Dict[str, Any]) -> List[CrossModalTrigger]:
    """Convenience wrapper: evaluate findings using the default engine."""
    return get_engine().evaluate_triggers(findings)


def format_report(triggers: List[CrossModalTrigger]) -> str:
    """Convenience wrapper: format triggers using the default engine."""
    return get_engine().format_trigger_report(triggers)
