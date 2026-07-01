"""Shared pytest fixtures for the Pharmacogenomics Intelligence Agent test suite.

Provides mock embedder, LLM client, collection manager, sample gene data,
sample patient profiles, and mock settings so tests run without Milvus,
Anthropic, or SentenceTransformer dependencies.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so ``from src.…`` / ``from config.…``
# imports work regardless of how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.models import CrossCollectionResult, SearchHit, PGxAlert, AlertLevel  # noqa: E402


# ═══════════════════════════════════════════════════════════════════════
# MOCK EMBEDDER
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def mock_embedder():
    """Return a mock embedder that produces 384-dim zero vectors."""
    embedder = MagicMock()
    embedder.embed_text.return_value = [0.0] * 384
    return embedder


# ═══════════════════════════════════════════════════════════════════════
# MOCK LLM CLIENT
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def mock_llm_client():
    """Return a mock LLM client that always responds with 'Mock response'."""
    client = MagicMock()
    client.generate.return_value = "Mock response"
    client.generate_stream.return_value = iter(["Mock ", "response"])
    return client


# ═══════════════════════════════════════════════════════════════════════
# MOCK COLLECTION MANAGER
# ═══════════════════════════════════════════════════════════════════════

COLLECTION_NAMES = [
    "pgx_gene_reference", "pgx_drug_guidelines", "pgx_drug_interactions",
    "pgx_hla_hypersensitivity", "pgx_phenoconversion", "pgx_dosing_algorithms",
    "pgx_clinical_evidence", "pgx_population_data", "pgx_clinical_trials",
    "pgx_fda_labels", "pgx_drug_alternatives", "pgx_patient_profiles",
    "pgx_implementation", "pgx_education", "genomic_evidence",
]


@pytest.fixture
def mock_collection_manager():
    """Return a MagicMock collection manager with sane defaults."""
    manager = MagicMock()
    manager.search.return_value = []
    manager.search_all.return_value = {name: [] for name in COLLECTION_NAMES}
    manager.get_collection_stats.return_value = {name: 100 for name in COLLECTION_NAMES}
    return manager


# ═══════════════════════════════════════════════════════════════════════
# SAMPLE DATA FIXTURES
# ═══════════════════════════════════════════════════════════════════════

@pytest.fixture
def sample_search_hits():
    """Return a list of sample SearchHit objects for testing."""
    return [
        SearchHit(collection="GeneReference", id="CYP2D6_star4", score=0.92,
                  text="CYP2D6 *4 no function splicing defect"),
        SearchHit(collection="DrugGuideline", id="cpic_codeine_2d6", score=0.88,
                  text="CPIC codeine CYP2D6 poor metabolizer avoid"),
        SearchHit(collection="Evidence", id="12345678", score=0.75,
                  text="Clinical evidence for CYP2D6 codeine interaction"),
    ]


@pytest.fixture
def sample_evidence(sample_search_hits):
    """Return a sample CrossCollectionResult."""
    return CrossCollectionResult(
        query="Is codeine safe for CYP2D6 poor metabolizer?",
        hits=sample_search_hits,
        knowledge_context="## CYP2D6 context",
        total_collections_searched=15,
        search_time_ms=42.5,
    )


@pytest.fixture
def sample_alerts():
    """Return sample PGxAlert instances."""
    return [
        PGxAlert(
            alert_level=AlertLevel.CRITICAL,
            gene="CYP2D6",
            drug="codeine",
            phenotype="Poor Metabolizer",
            recommendation="AVOID codeine. Insufficient conversion to morphine.",
        ),
        PGxAlert(
            alert_level=AlertLevel.WARNING,
            gene="CYP2C19",
            drug="clopidogrel",
            phenotype="Intermediate Metabolizer",
            recommendation="Consider alternative antiplatelet.",
        ),
    ]


@pytest.fixture
def sample_patient_profile():
    """Return a sample patient phenotype profile from PhenotypeTranslator."""
    return {
        "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1", "activity_score": 2.0},
        "CYP2C19": {"phenotype": "Intermediate Metabolizer", "diplotype": "*1/*2", "activity_score": 1.0},
        "CYP2C9": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1", "activity_score": 2.0},
        "CYP3A5": {"phenotype": "CYP3A5 Non-expresser", "diplotype": "*3/*3", "activity_score": 0.0},
        "DPYD": {"phenotype": "DPD Normal Activity", "diplotype": "*1/*1", "activity_score": 2.0},
    }
