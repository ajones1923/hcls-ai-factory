"""Shared pytest fixtures for CAR-T Intelligence Agent test suite.

Provides mock embedder, LLM client, collection manager, and sample
search results so that tests run without Milvus or external services.

Author: Adam Jones
Date: February 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so ``from src.…`` imports work
# regardless of how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.models import CrossCollectionResult, SearchHit  # noqa: E402


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


@pytest.fixture
def mock_collection_manager():
    """Return a MagicMock collection manager with sane defaults.

    - search()      -> empty list
    - search_all()  -> empty dict of lists for all 10 collections
    - get_collection_stats() -> dummy counts for all 10 collections
    - connect() / disconnect() -> no-ops
    """
    manager = MagicMock()

    manager.search.return_value = []

    collection_names = [
        "cart_literature",
        "cart_trials",
        "cart_constructs",
        "cart_assays",
        "cart_manufacturing",
        "cart_safety",
        "cart_biomarkers",
        "cart_regulatory",
        "cart_sequences",
        "cart_realworld",
    ]
    manager.search_all.return_value = {name: [] for name in collection_names}

    manager.get_collection_stats.return_value = {
        name: 42 for name in collection_names
    }

    manager.connect.return_value = None
    manager.disconnect.return_value = None

    return manager


# ═══════════════════════════════════════════════════════════════════════
# SAMPLE SEARCH DATA
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def sample_search_hits():
    """Return a list of 5 SearchHit objects spanning different collections."""
    return [
        SearchHit(
            collection="Literature",
            id="12345678",
            score=0.92,
            text="CD19 CAR-T therapy achieves high response rates in B-ALL.",
            metadata={"title": "CD19 CAR-T in B-ALL", "year": 2023, "target_antigen": "CD19"},
        ),
        SearchHit(
            collection="Trial",
            id="NCT03958656",
            score=0.87,
            text="Phase 2 study of tisagenlecleucel in pediatric B-ALL.",
            metadata={"phase": "Phase 2", "status": "Completed", "sponsor": "Novartis"},
        ),
        SearchHit(
            collection="Construct",
            id="construct-kymriah",
            score=0.83,
            text="Tisagenlecleucel: 4-1BB costimulatory domain, FMC63 scFv.",
            metadata={"name": "Kymriah", "generation": "2nd", "costimulatory_domain": "4-1BB"},
        ),
        SearchHit(
            collection="Safety",
            id="safety-crs-001",
            score=0.78,
            text="Grade 3+ CRS occurred in 22% of patients receiving tisagenlecleucel.",
            metadata={"product": "Kymriah", "event_type": "CRS", "severity_grade": "Grade 3-4"},
        ),
        SearchHit(
            collection="Manufacturing",
            id="mfg-lenti-001",
            score=0.71,
            text="Lentiviral transduction efficiency of 45% with MOI 5.",
            metadata={"process_step": "transduction", "parameter": "MOI", "parameter_value": "5"},
        ),
    ]


@pytest.fixture
def sample_evidence(sample_search_hits):
    """Return a CrossCollectionResult populated with 5 sample hits."""
    return CrossCollectionResult(
        query="What is the efficacy of CD19 CAR-T therapy?",
        hits=sample_search_hits,
        knowledge_context="## Target Antigen: CD19\n- **Protein:** B-Lymphocyte Antigen CD19",
        total_collections_searched=10,
        search_time_ms=42.5,
    )
