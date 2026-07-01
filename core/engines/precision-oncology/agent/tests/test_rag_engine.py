"""
Tests for the OncoRAGEngine multi-collection RAG engine.
=========================================================
Validates COLLECTION_CONFIG, weight sums, system prompt content,
prompt construction, and citation formatting.

All external dependencies (Milvus, LLM, embeddings) are mocked.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

# Import real models first so we can wire them into the agent.src.models alias
import src.models as _real_models

# We must mock the settings import before importing rag_engine
# because it references settings.WEIGHT_* at module level.
_mock_settings = MagicMock()
_mock_settings.WEIGHT_VARIANTS = 0.18
_mock_settings.WEIGHT_LITERATURE = 0.16
_mock_settings.WEIGHT_THERAPIES = 0.14
_mock_settings.WEIGHT_GUIDELINES = 0.12
_mock_settings.WEIGHT_TRIALS = 0.10
_mock_settings.WEIGHT_BIOMARKERS = 0.08
_mock_settings.WEIGHT_RESISTANCE = 0.07
_mock_settings.WEIGHT_PATHWAYS = 0.06
_mock_settings.WEIGHT_OUTCOMES = 0.04
_mock_settings.WEIGHT_CASES = 0.02
_mock_settings.WEIGHT_GENOMIC = 0.03

# Mock both config.settings AND agent.src.models (rag_engine uses agent.src.models)
_agent_mock = MagicMock()
_agent_src_mock = MagicMock()

with patch.dict("sys.modules", {
    "config.settings": MagicMock(settings=_mock_settings),
    "config": MagicMock(),
    "agent": _agent_mock,
    "agent.src": _agent_src_mock,
    "agent.src.models": _real_models,
}):
    sys.modules["config.settings"].settings = _mock_settings
    from src.rag_engine import (
        COLLECTION_CONFIG,
        ONCO_SYSTEM_PROMPT,
        OncoRAGEngine,
    )


# ═══════════════════════════════════════════════════════════════════════════
# COLLECTION_CONFIG Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCollectionConfig:
    """Verify COLLECTION_CONFIG has all collections with proper weights."""

    EXPECTED_COLLECTIONS = [
        "onco_variants",
        "onco_literature",
        "onco_therapies",
        "onco_guidelines",
        "onco_trials",
        "onco_biomarkers",
        "onco_resistance",
        "onco_pathways",
        "onco_outcomes",
        "onco_cases",
        "genomic_evidence",
    ]

    def test_all_collections_present(self):
        """COLLECTION_CONFIG should have all 11 collections."""
        for name in self.EXPECTED_COLLECTIONS:
            assert name in COLLECTION_CONFIG, f"Missing collection: {name}"

    def test_collection_count(self):
        """COLLECTION_CONFIG should have exactly 11 entries."""
        assert len(COLLECTION_CONFIG) == 11

    def test_each_collection_has_weight(self):
        """Every collection entry should have a 'weight' key."""
        for name, cfg in COLLECTION_CONFIG.items():
            assert "weight" in cfg, f"{name} missing 'weight' key"

    def test_each_collection_has_label(self):
        """Every collection entry should have a 'label' key."""
        for name, cfg in COLLECTION_CONFIG.items():
            assert "label" in cfg, f"{name} missing 'label' key"

    def test_weights_sum_approximately_one(self):
        """Collection weights should sum to approximately 1.0."""
        total = sum(cfg["weight"] for cfg in COLLECTION_CONFIG.values())
        assert abs(total - 1.0) < 0.05, (
            f"Weights sum to {total:.3f}, expected approximately 1.0"
        )

    def test_individual_weights_positive(self):
        """Each weight should be positive."""
        for name, cfg in COLLECTION_CONFIG.items():
            assert cfg["weight"] > 0, f"{name} weight should be positive"

    def test_variants_highest_weight(self):
        """onco_variants should have the highest weight."""
        variant_weight = COLLECTION_CONFIG["onco_variants"]["weight"]
        for name, cfg in COLLECTION_CONFIG.items():
            if name != "onco_variants":
                assert variant_weight >= cfg["weight"], (
                    f"onco_variants weight ({variant_weight}) should be >= {name} weight ({cfg['weight']})"
                )


# ═══════════════════════════════════════════════════════════════════════════
# System Prompt Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestSystemPrompt:
    """Verify the system prompt contains oncology-specific content."""

    def test_prompt_not_empty(self):
        assert len(ONCO_SYSTEM_PROMPT) > 100

    def test_contains_precision_oncology(self):
        assert "Oncology Intelligence Agent" in ONCO_SYSTEM_PROMPT

    def test_contains_molecular_profiling(self):
        assert "Molecular profiling" in ONCO_SYSTEM_PROMPT

    def test_contains_variant_interpretation(self):
        assert "Variant interpretation" in ONCO_SYSTEM_PROMPT

    def test_contains_therapy_selection(self):
        assert "Therapy selection" in ONCO_SYSTEM_PROMPT

    def test_contains_clinical_trial_matching(self):
        assert "Clinical trial matching" in ONCO_SYSTEM_PROMPT

    def test_contains_resistance_mechanisms(self):
        assert "Resistance mechanisms" in ONCO_SYSTEM_PROMPT

    def test_contains_biomarker_assessment(self):
        assert "Biomarker assessment" in ONCO_SYSTEM_PROMPT

    def test_contains_cite_evidence(self):
        assert "Cite evidence" in ONCO_SYSTEM_PROMPT

    def test_contains_nccn_reference(self):
        assert "NCCN" in ONCO_SYSTEM_PROMPT

    def test_contains_esmo_reference(self):
        assert "ESMO" in ONCO_SYSTEM_PROMPT

    def test_contains_uncertainty_acknowledgement(self):
        assert "uncertainty" in ONCO_SYSTEM_PROMPT.lower()


# ═══════════════════════════════════════════════════════════════════════════
# OncoRAGEngine Unit Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestOncoRAGEngineInit:
    """Test engine initialization with mock dependencies."""

    def test_initialization(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = OncoRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )
        assert engine.collection_manager is mock_collection_manager
        assert engine.embedder is mock_embedder
        assert engine.llm_client is mock_llm_client
        assert engine.knowledge is None
        assert engine.query_expander is None

    def test_initialization_with_knowledge(
        self, mock_collection_manager, mock_embedder, mock_llm_client
    ):
        mock_knowledge = MagicMock()
        engine = OncoRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            knowledge=mock_knowledge,
        )
        assert engine.knowledge is mock_knowledge

    def test_initialization_with_query_expander(
        self, mock_collection_manager, mock_embedder, mock_llm_client
    ):
        expander = MagicMock(return_value=["expanded term"])
        engine = OncoRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            query_expander=expander,
        )
        assert engine.query_expander is expander


class TestBuildPrompt:
    """Test _build_prompt includes evidence and query."""

    @pytest.fixture
    def engine(self, mock_collection_manager, mock_embedder, mock_llm_client):
        return OncoRAGEngine(
            collection_manager=mock_collection_manager,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

    def test_prompt_contains_question(self, engine):
        question = "What is the role of BRAF V600E in melanoma?"
        prompt = engine._build_prompt(question, [])
        assert question in prompt

    def test_prompt_contains_evidence_header(self, engine):
        prompt = engine._build_prompt("test question", [])
        assert "Retrieved Evidence" in prompt

    def test_prompt_contains_question_header(self, engine):
        prompt = engine._build_prompt("test question", [])
        assert "=== Question ===" in prompt

    def test_prompt_contains_instruction(self, engine):
        prompt = engine._build_prompt("test question", [])
        assert "well-cited answer" in prompt

    def test_prompt_with_evidence(self, engine, sample_search_hits):
        prompt = engine._build_prompt("EGFR therapy", sample_search_hits)
        # Should contain text from the evidence hits
        assert "BRAF V600E" in prompt or "EGFR" in prompt
        assert "score" in prompt


class TestScoreRelevance:
    """Test _score_relevance static method."""

    def test_high_relevance(self):
        assert OncoRAGEngine._score_relevance(0.90) == "high"

    def test_high_boundary(self):
        assert OncoRAGEngine._score_relevance(0.85) == "high"

    def test_medium_relevance(self):
        assert OncoRAGEngine._score_relevance(0.70) == "medium"

    def test_medium_boundary(self):
        assert OncoRAGEngine._score_relevance(0.65) == "medium"

    def test_low_relevance(self):
        assert OncoRAGEngine._score_relevance(0.40) == "low"

    def test_zero_score(self):
        assert OncoRAGEngine._score_relevance(0.0) == "low"


class TestFormatCitation:
    """Test _format_citation static method."""

    def test_pubmed_citation(self):
        citation = OncoRAGEngine._format_citation("onco_literature", "PMID:33096080")
        assert "pubmed.ncbi.nlm.nih.gov/33096080" in citation
        assert "PubMed" in citation

    def test_nct_citation(self):
        citation = OncoRAGEngine._format_citation("onco_trials", "NCT02628067")
        assert "clinicaltrials.gov" in citation
        assert "NCT02628067" in citation

    def test_generic_citation(self):
        citation = OncoRAGEngine._format_citation("onco_variants", "CIViC:12")
        assert "CIViC:12" in citation


class TestIsComparative:
    """Test _is_comparative static method."""

    def test_vs_query(self):
        assert OncoRAGEngine._is_comparative("osimertinib vs erlotinib") is True

    def test_versus_query(self):
        assert OncoRAGEngine._is_comparative("dabrafenib versus vemurafenib") is True

    def test_compare_query(self):
        assert OncoRAGEngine._is_comparative("compare BRAF and KRAS") is True

    def test_difference_query(self):
        assert OncoRAGEngine._is_comparative("difference between EGFR and ALK") is True

    def test_non_comparative_query(self):
        assert OncoRAGEngine._is_comparative("EGFR mutations in NSCLC") is False
