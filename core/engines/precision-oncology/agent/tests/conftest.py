"""
Shared pytest fixtures for the Oncology Intelligence Agent test suite.
=================================================================
Provides mock objects for Milvus, embeddings, LLM clients, and
reusable sample data following the CAR-T Intelligence Agent pattern.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

# ---------------------------------------------------------------------------
# Ensure src is importable
# ---------------------------------------------------------------------------
_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.models import CrossCollectionResult, SearchHit


# ---------------------------------------------------------------------------
# Collection name constants (all 11 collections)
# ---------------------------------------------------------------------------
ALL_COLLECTION_NAMES = [
    "onco_literature",
    "onco_trials",
    "onco_variants",
    "onco_biomarkers",
    "onco_therapies",
    "onco_pathways",
    "onco_guidelines",
    "onco_resistance",
    "onco_outcomes",
    "onco_cases",
    "genomic_evidence",
]


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def mock_embedder():
    """Mock embedding model that returns 384-dimensional zero vectors."""
    embedder = MagicMock()
    embedder.embed_text.return_value = [0.0] * 384
    embedder.encode.return_value = [0.0] * 384
    embedder.embed.return_value = [0.0] * 384
    return embedder


@pytest.fixture
def mock_llm_client():
    """Mock LLM client that returns a fixed response string."""
    client = MagicMock()
    client.chat.return_value = "Mock response"
    client.chat_stream.return_value = iter(["Mock ", "response"])
    return client


@pytest.fixture
def mock_collection_manager():
    """Mock Milvus collection manager with all 11 collection names.

    - search() returns an empty list
    - search_all() returns a dict mapping each collection to an empty list
    - get_collection_stats() returns 42 entities for each collection
    - query() returns an empty list
    - insert() / insert_batch() returns 0
    """
    manager = MagicMock()

    # search returns empty list
    manager.search.return_value = []

    # search_all returns dict with empty lists per collection
    manager.search_all.return_value = {name: [] for name in ALL_COLLECTION_NAMES}

    # get_collection_stats returns 42 entities for each collection
    def _stats(name, **kwargs):
        return {
            "name": name,
            "num_entities": 42,
            "fields": ["id", "embedding", "text_chunk"],
        }

    manager.get_collection_stats.side_effect = _stats

    # query returns empty list
    manager.query.return_value = []

    # insert / insert_batch returns 0
    manager.insert.return_value = 0
    manager.insert_batch.return_value = 0

    return manager


@pytest.fixture
def sample_search_hits():
    """Five SearchHit objects spanning different oncology collections."""
    return [
        SearchHit(
            collection="onco_literature",
            id="PMID:33096080",
            score=0.92,
            text="BRAF V600E is the most common activating mutation in melanoma, "
                 "occurring in approximately 50% of cases. Combination BRAF+MEK "
                 "inhibitor therapy significantly improves progression-free survival.",
            metadata={"gene": "BRAF", "variant": "V600E", "cancer_type": "melanoma"},
        ),
        SearchHit(
            collection="onco_trials",
            id="NCT02628067",
            score=0.88,
            text="KEYNOTE-158: Phase II basket trial of pembrolizumab in patients "
                 "with advanced solid tumors harboring MSI-H or TMB-H biomarkers.",
            metadata={"phase": "Phase 2", "status": "Active, not recruiting"},
        ),
        SearchHit(
            collection="onco_variants",
            id="CIViC:12",
            score=0.85,
            text="EGFR L858R is a sensitising mutation in exon 21 of EGFR, "
                 "conferring responsiveness to EGFR tyrosine kinase inhibitors "
                 "including osimertinib, erlotinib, gefitinib, and afatinib.",
            metadata={"gene": "EGFR", "variant_type": "snv", "evidence_level": "A"},
        ),
        SearchHit(
            collection="onco_therapies",
            id="THERAPY:osimertinib",
            score=0.81,
            text="Osimertinib (Tagrisso) is a 3rd-generation EGFR TKI approved for "
                 "first-line treatment of EGFR-mutant NSCLC (FLAURA trial).",
            metadata={"drug_name": "osimertinib", "category": "targeted"},
        ),
        SearchHit(
            collection="onco_guidelines",
            id="NCCN:NSCLC:2025.2",
            score=0.78,
            text="NCCN NSCLC Guidelines v2.2025 recommend molecular testing for "
                 "EGFR, ALK, ROS1, BRAF, KRAS G12C, NTRK, MET, RET, HER2 in "
                 "non-squamous and never/light-smoker squamous NSCLC.",
            metadata={"org": "NCCN", "cancer_type": "nsclc", "year": 2025},
        ),
    ]


@pytest.fixture
def sample_evidence(sample_search_hits):
    """CrossCollectionResult populated with sample oncology hits."""
    return CrossCollectionResult(
        query="BRAF V600E targeted therapy melanoma",
        hits=sample_search_hits,
        total_collections_searched=5,
        search_time_ms=42.7,
    )


@pytest.fixture
def sample_settings():
    """Create an OncoSettings instance with test defaults.

    Uses environment override to avoid touching real config files.
    """
    try:
        from config.settings import OncoSettings
        return OncoSettings(
            MILVUS_HOST="localhost",
            MILVUS_PORT=19530,
            EMBEDDING_DIM=384,
            TOP_K=5,
            SCORE_THRESHOLD=0.4,
            LLM_PROVIDER="anthropic",
            LLM_MODEL="claude-sonnet-4-20250514",
        )
    except Exception:
        # Fallback: return a MagicMock with the expected attributes
        s = MagicMock()
        s.MILVUS_HOST = "localhost"
        s.MILVUS_PORT = 19530
        s.EMBEDDING_DIM = 384
        s.TOP_K = 5
        s.SCORE_THRESHOLD = 0.4
        s.LLM_PROVIDER = "anthropic"
        s.LLM_MODEL = "claude-sonnet-4-20250514"
        s.COLLECTION_LITERATURE = "onco_literature"
        s.COLLECTION_TRIALS = "onco_trials"
        s.COLLECTION_VARIANTS = "onco_variants"
        s.COLLECTION_BIOMARKERS = "onco_biomarkers"
        s.COLLECTION_THERAPIES = "onco_therapies"
        s.COLLECTION_PATHWAYS = "onco_pathways"
        s.COLLECTION_GUIDELINES = "onco_guidelines"
        s.COLLECTION_RESISTANCE = "onco_resistance"
        s.COLLECTION_OUTCOMES = "onco_outcomes"
        s.COLLECTION_CASES = "onco_cases"
        s.COLLECTION_GENOMIC = "genomic_evidence"
        s.WEIGHT_VARIANTS = 0.18
        s.WEIGHT_LITERATURE = 0.16
        s.WEIGHT_THERAPIES = 0.14
        s.WEIGHT_GUIDELINES = 0.12
        s.WEIGHT_TRIALS = 0.10
        s.WEIGHT_BIOMARKERS = 0.08
        s.WEIGHT_RESISTANCE = 0.07
        s.WEIGHT_PATHWAYS = 0.06
        s.WEIGHT_OUTCOMES = 0.04
        s.WEIGHT_CASES = 0.02
        s.WEIGHT_GENOMIC = 0.03
        return s
