"""Unified demo patient data for cross-agent demonstration.

Provides a consistent VCP/FTD demo patient (DEMO-VCP-001) that flows
through all 5 intelligence agents with domain-specific data.
"""

from __future__ import annotations

# ─── Core Patient Identity ───────────────────────────────────────
DEMO_PATIENT_ID = "DEMO-VCP-001"
DEMO_PATIENT_NAME = "Demo Patient (VCP p.R155H)"
DEMO_PATIENT_AGE = 52
DEMO_PATIENT_SEX = "M"

# ─── Biomarker Agent Data ────────────────────────────────────────
DEMO_BIOMARKERS = {
    "albumin": 3.8,
    "creatinine": 1.1,
    "glucose": 108.0,
    "c_reactive_protein": 4.2,
    "lymphocyte_pct": 28.0,
    "mcv": 91.0,
    "rdw": 14.1,
    "alkaline_phosphatase": 72.0,
    "wbc": 6.8,
    "platelets": 245.0,
    "alt": 32.0,
    "ast": 28.0,
    "total_cholesterol": 218.0,
    "ldl": 142.0,
    "hdl": 48.0,
    "triglycerides": 165.0,
    "hs_crp": 3.8,
    "d_dimer": 0.42,
    "troponin_i": 0.01,
    "nt_probnp": 125.0,
    "hba1c": 5.9,
    "tsh": 2.1,
    "ferritin": 180.0,
    "vitamin_d": 28.0,
    "homocysteine": 14.2,
}

DEMO_GENOTYPES = {
    "rs1801133": "CT",   # MTHFR C677T heterozygous
    "rs4149056": "TC",   # SLCO1B1 - intermediate function
    "rs9923231": "CT",   # VKORC1 - intermediate sensitivity
    "rs1799853": "CT",   # CYP2C9 *2
}

DEMO_STAR_ALLELES = {
    "CYP2D6": "*4/*4",     # Poor metabolizer
    "CYP2C19": "*1/*2",    # Intermediate metabolizer
    "CYP2C9": "*1/*2",     # Intermediate metabolizer
    "SLCO1B1": "*1/*5",    # Intermediate function
    "VKORC1": "AG",        # Intermediate warfarin sensitivity
    "DPYD": "*1/*1",       # Normal
    "CYP3A5": "*3/*3",     # Poor expresser
    "TPMT": "*1/*1",       # Normal
}

# ─── Oncology Agent Data ─────────────────────────────────────────
DEMO_ONCOLOGY = {
    "cancer_type": "NSCLC",
    "stage": "IIIA",
    "variants": [
        {"gene": "VCP", "variant": "p.R155H", "vaf": 0.35, "type": "missense"},
        {"gene": "TP53", "variant": "p.R248W", "vaf": 0.22, "type": "missense"},
        {"gene": "EGFR", "variant": "L858R", "vaf": 0.41, "type": "missense"},
    ],
    "biomarkers": {
        "tmb": 8.2,
        "msi_status": "MSS",
        "pd_l1_tps": 45,
        "hrd_score": 12,
    },
    "prior_therapies": ["Carboplatin/Pemetrexed"],
}

# ─── CAR-T Agent Data ────────────────────────────────────────────
DEMO_CART = {
    "target_antigen": "CD19",
    "indication": "B-cell ALL relapsed/refractory",
    "question": "What are the latest CAR-T constructs targeting CD19 with reduced CRS risk for relapsed B-cell ALL?",
}

# ─── Imaging Agent Data ─────────────────────────────────────────
DEMO_IMAGING = {
    "modality": "MRI",
    "body_region": "brain",
    "question": "Analyze brain MRI findings in a 52-year-old male with VCP mutation and suspected frontotemporal dementia. What patterns of cortical atrophy are expected?",
}

# ─── Autoimmune Agent Data ───────────────────────────────────────
DEMO_AUTOIMMUNE = {
    "primary_diagnosis": "Inclusion Body Myopathy (VCP-related)",
    "hla_alleles": ["B*08:01", "DRB1*03:01"],
    "autoantibodies": {
        "anti_nuclear": {"titer": "1:160", "pattern": "speckled"},
        "anti_jo1": {"result": "negative"},
        "crp": 4.2,
    },
}


def get_demo_patient_summary() -> str:
    """Return a one-line summary of the demo patient."""
    return (
        f"{DEMO_PATIENT_ID}: {DEMO_PATIENT_AGE}yo {DEMO_PATIENT_SEX}, "
        f"VCP p.R155H, CYP2D6 *4/*4 (PM), pre-diabetic trajectory"
    )
