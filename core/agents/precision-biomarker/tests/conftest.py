"""Shared pytest fixtures for Biomarker Intelligence Agent test suite.

Provides mock embedder, LLM client, collection manager, and sample
patient/biomarker data so that tests run without Milvus or external services.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so ``from src.`` imports work
# regardless of how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.models import (  # noqa: E402
    CrossCollectionResult,
    PatientProfile,
    SearchHit,
)

# =====================================================================
# MOCK EMBEDDER
# =====================================================================


@pytest.fixture
def mock_embedder():
    """Return a mock embedder that produces 384-dim zero vectors."""
    embedder = MagicMock()
    embedder.embed_text.return_value = [0.0] * 384
    embedder.encode.return_value = [[0.0] * 384]
    return embedder


# =====================================================================
# MOCK LLM CLIENT
# =====================================================================


@pytest.fixture
def mock_llm_client():
    """Return a mock LLM client that always responds with 'Mock response'."""
    client = MagicMock()
    client.generate.return_value = "Mock response"
    client.generate_stream.return_value = iter(["Mock ", "response"])
    return client


# =====================================================================
# MOCK COLLECTION MANAGER
# =====================================================================


@pytest.fixture
def mock_collection_manager():
    """Return a MagicMock collection manager with sane defaults.

    - search()      -> empty list
    - search_all()  -> empty dict of lists for all 14 collections
    - get_collection_stats() -> dummy counts for all 14 collections
    - connect() / disconnect() -> no-ops
    """
    manager = MagicMock()

    manager.search.return_value = []

    collection_names = [
        "biomarker_reference",
        "biomarker_genetic_variants",
        "biomarker_pgx_rules",
        "biomarker_disease_trajectories",
        "biomarker_clinical_evidence",
        "biomarker_nutrition",
        "biomarker_drug_interactions",
        "biomarker_aging_markers",
        "biomarker_genotype_adjustments",
        "biomarker_monitoring",
        "biomarker_critical_values",
        "biomarker_discordance_rules",
        "biomarker_aj_carrier_screening",
        "genomic_evidence",
    ]
    manager.search_all.return_value = {name: [] for name in collection_names}

    manager.get_collection_stats.return_value = {
        name: 42 for name in collection_names
    }

    manager.connect.return_value = None
    manager.disconnect.return_value = None

    return manager


# =====================================================================
# SAMPLE PATIENT PROFILE
# =====================================================================


@pytest.fixture
def sample_patient_profile():
    """Return a complete PatientProfile with biomarkers, genotypes, star alleles."""
    return PatientProfile(
        patient_id="HG002",
        age=45,
        sex="M",
        biomarkers={
            "albumin": 4.2,
            "creatinine": 0.9,
            "glucose": 105.0,
            "hs_crp": 2.1,
            "lymphocyte_pct": 28.0,
            "mcv": 92.0,
            "rdw": 13.5,
            "alkaline_phosphatase": 72.0,
            "wbc": 6.8,
            "hba1c": 5.9,
            "fasting_glucose": 105.0,
            "fasting_insulin": 12.0,
            "ldl_c": 145.0,
            "hdl_c": 42.0,
            "triglycerides": 180.0,
            "lpa": 65.0,
            "tsh": 3.2,
            "free_t4": 1.1,
            "free_t3": 2.4,
            "alt": 38.0,
            "ast": 32.0,
            "ferritin": 280.0,
            "transferrin_saturation": 35.0,
            "vitamin_d": 22.0,
            "vitamin_b12": 350.0,
            "folate": 8.0,
            "omega3_index": 4.5,
            "magnesium": 1.9,
        },
        genotypes={
            "TCF7L2_rs7903146": "CT",
            "PNPLA3_rs738409": "CG",
            "DIO2_rs225014": "GA",
            "APOE": "E3/E4",
            "HFE_rs1800562": "GG",
            "MTHFR_rs1801133": "CT",
            "FADS1_rs174546": "CT",
        },
        star_alleles={
            "CYP2D6": "*1/*4",
            "CYP2C19": "*1/*2",
            "SLCO1B1": "*1/*5",
            "VKORC1": "A/G",
            "MTHFR": "CT",
            "TPMT": "*1/*1",
            "CYP2C9": "*1/*1",
        },
    )


# =====================================================================
# SAMPLE BIOMARKERS
# =====================================================================


@pytest.fixture
def sample_biomarkers():
    """Return a dict of common biomarker values for a 45-year-old male."""
    return {
        "albumin": 4.2,
        "creatinine": 0.9,
        "glucose": 95.0,
        "hs_crp": 1.5,
        "lymphocyte_pct": 30.0,
        "mcv": 89.0,
        "rdw": 13.0,
        "alkaline_phosphatase": 65.0,
        "wbc": 6.0,
        "hba1c": 5.4,
        "fasting_glucose": 95.0,
        "fasting_insulin": 8.0,
        "ldl_c": 120.0,
        "hdl_c": 55.0,
        "triglycerides": 130.0,
        "tsh": 2.5,
        "free_t4": 1.2,
        "free_t3": 3.0,
        "alt": 25.0,
        "ast": 22.0,
        "ferritin": 150.0,
        "transferrin_saturation": 30.0,
        "vitamin_d": 35.0,
        "vitamin_b12": 500.0,
        "folate": 12.0,
        "omega3_index": 6.5,
        "magnesium": 2.1,
    }


# =====================================================================
# SAMPLE GENOTYPES
# =====================================================================


@pytest.fixture
def sample_genotypes():
    """Return a dict of common genotypes (rsID -> genotype string)."""
    return {
        "TCF7L2_rs7903146": "CT",
        "PNPLA3_rs738409": "CC",
        "DIO2_rs225014": "GG",
        "APOE": "E3/E3",
        "HFE_rs1800562": "GG",
        "HFE_rs1799945": "CC",
        "MTHFR_rs1801133": "CC",
        "FADS1_rs174546": "TT",
        "FADS2_rs1535": "GG",
        "PPARG_rs1801282": "CC",
        "SLC30A8_rs13266634": "CC",
        "KCNJ11_rs5219": "CC",
        "PCSK9_rs11591147": "GG",
    }


# =====================================================================
# SAMPLE STAR ALLELES
# =====================================================================


@pytest.fixture
def sample_star_alleles():
    """Return a dict of star alleles for PGx testing."""
    return {
        "CYP2D6": "*1/*1",
        "CYP2C19": "*1/*1",
        "SLCO1B1": "*1/*1",
        "VKORC1": "G/G",
        "MTHFR": "CC",
        "TPMT": "*1/*1",
        "CYP2C9": "*1/*1",
    }


# =====================================================================
# SAMPLE SEARCH DATA
# =====================================================================


@pytest.fixture
def sample_search_hits():
    """Return a list of 5 SearchHit objects spanning different collections."""
    return [
        SearchHit(
            collection="biomarker_reference",
            id="ref-hba1c",
            score=0.92,
            text="HbA1c (glycated hemoglobin) reflects average blood glucose over 2-3 months.",
            metadata={"name": "HbA1c", "unit": "%", "category": "Diabetes"},
        ),
        SearchHit(
            collection="biomarker_genetic_variants",
            id="var-mthfr-c677t",
            score=0.87,
            text="MTHFR C677T reduces enzyme activity by 30-70%, affecting folate metabolism.",
            metadata={"gene": "MTHFR", "rs_id": "rs1801133"},
        ),
        SearchHit(
            collection="biomarker_pgx_rules",
            id="pgx-cyp2d6-codeine",
            score=0.83,
            text="CYP2D6 poor metabolizers should avoid codeine due to lack of analgesic effect.",
            metadata={"gene": "CYP2D6", "drug": "Codeine"},
        ),
        SearchHit(
            collection="biomarker_clinical_evidence",
            id="evi-phenoage-001",
            score=0.78,
            text="PhenoAge uses 9 blood biomarkers to estimate biological age (Levine 2018).",
            metadata={"pmid": "29676998", "year": 2018},
        ),
        SearchHit(
            collection="biomarker_disease_trajectories",
            id="traj-diabetes-pre",
            score=0.71,
            text="Pre-diabetes stage: HbA1c 5.7-6.4%, fasting glucose 100-125 mg/dL.",
            metadata={"disease": "diabetes", "stage": "pre_diabetic"},
        ),
    ]


@pytest.fixture
def sample_evidence(sample_search_hits):
    """Return a CrossCollectionResult populated with 5 sample hits."""
    return CrossCollectionResult(
        query="What biomarkers indicate pre-diabetes risk?",
        hits=sample_search_hits,
        knowledge_context="## Diabetes Biomarkers\n- HbA1c, fasting glucose, HOMA-IR",
        total_collections_searched=10,
        search_time_ms=42.5,
    )
