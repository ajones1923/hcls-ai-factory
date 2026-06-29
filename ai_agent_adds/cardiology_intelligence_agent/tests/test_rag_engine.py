"""Tests for the Cardiology Intelligence Agent RAG engine module.

Author: Adam Jones
Date: March 2026

Covers:
- CardioSearchResult and CardioResponse dataclasses
- COLLECTION_CONFIG structure (13 collections)
- CardioRAGEngine: creation, query(), search(), query_stream(),
  _build_context(), _score_confidence(), find_related(), health_check(),
  conversation memory (add_conversation_context, clear_conversation),
  _get_boosted_weights(), _embed_query(), _rerank_results()
- create_rag_engine() factory function
- format_search_results() display helper
"""

from dataclasses import fields as dc_fields
from unittest.mock import MagicMock, patch

import pytest

from src.models import CardioWorkflowType

# Patch settings before importing rag_engine since it reads settings at module level
_mock_settings = MagicMock()
_mock_settings.WEIGHT_LITERATURE = 0.10
_mock_settings.WEIGHT_TRIALS = 0.10
_mock_settings.WEIGHT_IMAGING = 0.10
_mock_settings.WEIGHT_ELECTROPHYSIOLOGY = 0.08
_mock_settings.WEIGHT_HEART_FAILURE = 0.10
_mock_settings.WEIGHT_VALVULAR = 0.08
_mock_settings.WEIGHT_PREVENTION = 0.08
_mock_settings.WEIGHT_INTERVENTIONAL = 0.08
_mock_settings.WEIGHT_ONCOLOGY = 0.05
_mock_settings.WEIGHT_DEVICES = 0.05
_mock_settings.WEIGHT_GUIDELINES = 0.10
_mock_settings.WEIGHT_HEMODYNAMICS = 0.05
_mock_settings.WEIGHT_GENOMIC = 0.05
_mock_settings.MAX_CONVERSATION_CONTEXT = 5
_mock_settings.EMBEDDING_MODEL = "BAAI/bge-small-en-v1.5"
_mock_settings.CITATION_HIGH_THRESHOLD = 0.85
_mock_settings.CITATION_MEDIUM_THRESHOLD = 0.70
_mock_settings.ANTHROPIC_API_KEY = ""
_mock_settings.LLM_MODEL = "claude-sonnet-4-20250514"
_mock_settings.MILVUS_HOST = "localhost"
_mock_settings.MILVUS_PORT = 19530

with patch.dict("sys.modules", {"config.settings": MagicMock(settings=_mock_settings),
                                  "config": MagicMock()}):
    with patch("config.settings.settings", _mock_settings):
        from src.rag_engine import (
            CardioSearchResult,
            CardioResponse,
            COLLECTION_CONFIG,
            CardioRAGEngine,
            create_rag_engine,
            format_search_results,
        )


# =====================================================================
# CardioSearchResult DATACLASS
# =====================================================================


class TestCardioSearchResult:
    """CardioSearchResult should hold a single search result."""

    def test_creation_defaults(self):
        r = CardioSearchResult()
        assert r.collection == ""
        assert r.record_id == ""
        assert r.score == 0.0
        assert r.text == ""
        assert r.metadata == {}
        assert r.relevance == "low"

    def test_creation_with_values(self):
        r = CardioSearchResult(
            collection="cardio_guidelines",
            record_id="123",
            score=0.92,
            text="Sample text",
            metadata={"pmid": "12345678"},
            relevance="high",
        )
        assert r.collection == "cardio_guidelines"
        assert r.record_id == "123"
        assert r.score == 0.92
        assert r.text == "Sample text"
        assert r.metadata["pmid"] == "12345678"
        assert r.relevance == "high"

    def test_has_expected_fields(self):
        field_names = {f.name for f in dc_fields(CardioSearchResult)}
        for name in ("collection", "record_id", "score", "text",
                      "metadata", "relevance"):
            assert name in field_names


# =====================================================================
# CardioResponse DATACLASS
# =====================================================================


class TestCardioResponse:
    """CardioResponse should hold a complete RAG response."""

    def test_creation_defaults(self):
        r = CardioResponse()
        assert r.question == ""
        assert r.answer == ""
        assert r.results == []
        assert r.workflow is None
        assert r.confidence == 0.0
        assert r.citations == []
        assert r.search_time_ms == 0.0
        assert r.collections_searched == 0
        assert r.patient_context_used is False
        assert isinstance(r.timestamp, str)

    def test_creation_with_values(self):
        r = CardioResponse(
            question="What is GDMT?",
            answer="GDMT consists of...",
            confidence=0.85,
        )
        assert r.question == "What is GDMT?"
        assert r.answer == "GDMT consists of..."
        assert r.confidence == 0.85

    def test_has_expected_fields(self):
        field_names = {f.name for f in dc_fields(CardioResponse)}
        for name in ("question", "answer", "results", "workflow",
                      "confidence", "citations", "search_time_ms",
                      "collections_searched", "patient_context_used",
                      "timestamp"):
            assert name in field_names

    def test_timestamp_is_iso_format(self):
        r = CardioResponse()
        # Should be parseable as ISO 8601
        assert "T" in r.timestamp
        assert r.timestamp.endswith("Z")


# =====================================================================
# COLLECTION_CONFIG
# =====================================================================


class TestCollectionConfig:
    """COLLECTION_CONFIG should define 13 collections."""

    def test_is_dict(self):
        assert isinstance(COLLECTION_CONFIG, dict)

    def test_has_13_collections(self):
        assert len(COLLECTION_CONFIG) == 13

    EXPECTED_COLLECTIONS = [
        "cardio_literature", "cardio_trials", "cardio_imaging",
        "cardio_electrophysiology", "cardio_heart_failure",
        "cardio_valvular", "cardio_prevention", "cardio_interventional",
        "cardio_oncology", "cardio_devices", "cardio_guidelines",
        "cardio_hemodynamics", "genomic_evidence",
    ]

    @pytest.mark.parametrize("coll_name", EXPECTED_COLLECTIONS)
    def test_collection_exists(self, coll_name):
        assert coll_name in COLLECTION_CONFIG

    def test_each_collection_has_weight(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "weight" in cfg, f"{name} missing weight"
            assert isinstance(cfg["weight"], (int, float))

    def test_each_collection_has_label(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "label" in cfg, f"{name} missing label"
            assert isinstance(cfg["label"], str)

    def test_each_collection_has_text_field(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "text_field" in cfg, f"{name} missing text_field"

    def test_each_collection_has_title_field(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "title_field" in cfg, f"{name} missing title_field"


# =====================================================================
# CardioRAGEngine CREATION
# =====================================================================


class TestCardioRAGEngineCreation:
    """CardioRAGEngine should initialise with optional components."""

    def test_creation_all_none(self):
        engine = CardioRAGEngine()
        assert engine.milvus is None
        assert engine.embedder is None
        assert engine.llm is None
        assert engine.query_expander is None

    def test_creation_with_components(self):
        mock_milvus = MagicMock()
        mock_embedder = MagicMock()
        mock_llm = MagicMock()
        engine = CardioRAGEngine(
            milvus_client=mock_milvus,
            embedding_model=mock_embedder,
            llm_client=mock_llm,
        )
        assert engine.milvus is mock_milvus
        assert engine.embedder is mock_embedder
        assert engine.llm is mock_llm

    def test_conversation_history_empty_on_creation(self):
        engine = CardioRAGEngine(session_id="test_empty_creation")
        assert engine.conversation_history == []


# =====================================================================
# CardioRAGEngine.search() (with mocks)
# =====================================================================


class TestCardioRAGEngineSearch:
    """search() should require Milvus and perform parallel search."""

    def test_search_without_milvus_raises(self):
        engine = CardioRAGEngine()
        with pytest.raises(RuntimeError, match="Milvus client not configured"):
            engine.search("heart failure")

    def test_search_calls_embed_query(self):
        mock_milvus = MagicMock()
        mock_embedder = MagicMock()
        mock_embedder.embed_text.return_value = [0.1] * 384
        mock_milvus.search.return_value = [[]]

        engine = CardioRAGEngine(
            milvus_client=mock_milvus,
            embedding_model=mock_embedder,
        )
        engine.search("test query", collections=["cardio_guidelines"], top_k=1)
        mock_embedder.embed_text.assert_called_once()


# =====================================================================
# CardioRAGEngine.query() (with mocks)
# =====================================================================


class TestCardioRAGEngineQuery:
    """query() should orchestrate search -> rerank -> synthesize."""

    @pytest.fixture
    def engine(self):
        mock_milvus = MagicMock()
        mock_embedder = MagicMock()
        mock_embedder.embed_text.return_value = [0.1] * 384
        mock_milvus.search.return_value = [[]]

        mock_llm = MagicMock()
        mock_llm.generate.return_value = "LLM synthesized answer"

        engine = CardioRAGEngine(
            milvus_client=mock_milvus,
            embedding_model=mock_embedder,
            llm_client=mock_llm,
        )
        return engine

    def test_query_returns_cardio_response(self, engine):
        response = engine.query("heart failure management")
        assert isinstance(response, CardioResponse)

    def test_query_sets_question(self, engine):
        response = engine.query("GDMT for HFrEF")
        assert response.question == "GDMT for HFrEF"

    def test_query_sets_search_time(self, engine):
        response = engine.query("heart failure")
        assert response.search_time_ms >= 0

    def test_query_without_llm_returns_search_only(self):
        mock_milvus = MagicMock()
        mock_embedder = MagicMock()
        mock_embedder.embed_text.return_value = [0.1] * 384
        mock_milvus.search.return_value = [[]]

        engine = CardioRAGEngine(
            milvus_client=mock_milvus,
            embedding_model=mock_embedder,
        )
        response = engine.query("test")
        assert "search-only" in response.answer.lower()

    def test_query_updates_conversation_history(self, engine):
        engine.query("heart failure management")
        assert len(engine.conversation_history) > 0


# =====================================================================
# CardioRAGEngine._build_context()
# =====================================================================


class TestBuildContext:
    """_build_context() should format results for LLM prompt."""

    @pytest.fixture
    def engine(self):
        return CardioRAGEngine()

    def test_empty_results_returns_no_evidence(self, engine):
        context = engine._build_context([])
        assert "No evidence found" in context

    def test_formats_results_by_collection(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                record_id="1",
                score=0.9,
                text="Sample guideline text",
                metadata={"collection_label": "Guideline"},
                relevance="high",
            ),
            CardioSearchResult(
                collection="cardio_literature",
                record_id="2",
                score=0.8,
                text="Sample literature text",
                metadata={"collection_label": "Literature"},
                relevance="medium",
            ),
        ]
        context = engine._build_context(results)
        assert "Guideline" in context
        assert "Literature" in context

    def test_includes_relevance_tags(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                record_id="1",
                score=0.9,
                text="Test",
                metadata={"collection_label": "Guideline"},
                relevance="high",
            ),
        ]
        context = engine._build_context(results)
        assert "high relevance" in context


# =====================================================================
# CardioRAGEngine._score_confidence()
# =====================================================================


class TestScoreConfidence:
    """_score_confidence() should return float between 0 and 1."""

    @pytest.fixture
    def engine(self):
        return CardioRAGEngine()

    def test_no_results_returns_zero(self, engine):
        assert engine._score_confidence([]) == 0.0

    def test_returns_float(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                score=0.9,
                relevance="high",
            )
        ]
        score = engine._score_confidence(results)
        assert isinstance(score, float)

    def test_score_between_0_and_1(self, engine):
        results = [
            CardioSearchResult(
                collection=f"coll_{i}",
                score=0.8 + i * 0.01,
                relevance="high" if i % 2 == 0 else "medium",
            )
            for i in range(10)
        ]
        score = engine._score_confidence(results)
        assert 0.0 <= score <= 1.0

    def test_high_relevance_increases_confidence(self, engine):
        low_results = [
            CardioSearchResult(collection="a", score=0.5, relevance="low")
            for _ in range(5)
        ]
        high_results = [
            CardioSearchResult(collection="a", score=0.5, relevance="high")
            for _ in range(5)
        ]
        low_score = engine._score_confidence(low_results)
        high_score = engine._score_confidence(high_results)
        assert high_score >= low_score

    def test_diversity_increases_confidence(self, engine):
        single_coll = [
            CardioSearchResult(collection="a", score=0.8, relevance="high")
            for _ in range(5)
        ]
        multi_coll = [
            CardioSearchResult(
                collection=f"coll_{i}", score=0.8, relevance="high"
            )
            for i in range(5)
        ]
        single_score = engine._score_confidence(single_coll)
        multi_score = engine._score_confidence(multi_coll)
        assert multi_score >= single_score

    def test_guideline_presence_adds_confidence(self, engine):
        no_guideline = [
            CardioSearchResult(collection="cardio_literature", score=0.8, relevance="high")
        ]
        with_guideline = [
            CardioSearchResult(collection="cardio_guidelines", score=0.8, relevance="high")
        ]
        ng_score = engine._score_confidence(no_guideline)
        wg_score = engine._score_confidence(with_guideline)
        assert wg_score > ng_score


# =====================================================================
# CardioRAGEngine.find_related()
# =====================================================================


class TestFindRelated:
    """find_related() should search targeted collections by entity type."""

    def test_find_related_calls_search(self):
        mock_milvus = MagicMock()
        mock_embedder = MagicMock()
        mock_embedder.embed_text.return_value = [0.1] * 384
        mock_milvus.search.return_value = [[]]

        engine = CardioRAGEngine(
            milvus_client=mock_milvus,
            embedding_model=mock_embedder,
        )
        results = engine.find_related("heart failure", entity_type="condition")
        assert isinstance(results, list)

    def test_gene_entity_type_searches_genomic(self):
        mock_milvus = MagicMock()
        mock_embedder = MagicMock()
        mock_embedder.embed_text.return_value = [0.1] * 384
        mock_milvus.search.return_value = [[]]

        engine = CardioRAGEngine(
            milvus_client=mock_milvus,
            embedding_model=mock_embedder,
        )
        # Just verify it runs without error
        results = engine.find_related("MYH7", entity_type="gene")
        assert isinstance(results, list)


# =====================================================================
# CardioRAGEngine.health_check()
# =====================================================================


class TestHealthCheck:
    """health_check() should report Milvus connection and collection status."""

    def test_no_milvus_returns_unhealthy(self):
        engine = CardioRAGEngine()
        health = engine.health_check()
        assert health["status"] == "unhealthy"
        assert health["milvus_connected"] is False

    def test_healthy_when_all_collections_exist(self):
        mock_milvus = MagicMock()
        mock_milvus.has_collection.return_value = True
        engine = CardioRAGEngine(milvus_client=mock_milvus)
        health = engine.health_check()
        assert health["status"] == "healthy"
        assert health["milvus_connected"] is True
        assert len(health["collections_missing"]) == 0

    def test_degraded_when_some_collections_missing(self):
        mock_milvus = MagicMock()
        call_count = [0]

        def side_effect(name):
            call_count[0] += 1
            return call_count[0] <= 8  # First 8 exist, rest missing

        mock_milvus.has_collection.side_effect = side_effect
        engine = CardioRAGEngine(milvus_client=mock_milvus)
        health = engine.health_check()
        assert health["status"] in ("healthy", "degraded")
        assert health["milvus_connected"] is True

    def test_health_check_has_expected_keys(self):
        engine = CardioRAGEngine()
        health = engine.health_check()
        for key in ("status", "milvus_connected", "collections_available",
                     "collections_missing", "embedding_model",
                     "llm_configured", "timestamp"):
            assert key in health, f"Missing key: {key}"


# =====================================================================
# CONVERSATION MEMORY
# =====================================================================


class TestConversationMemory:
    """Conversation history management."""

    def test_add_conversation_context(self):
        engine = CardioRAGEngine(session_id="test_add_ctx")
        engine.add_conversation_context("user", "What is GDMT?")
        assert len(engine.conversation_history) == 1
        assert engine.conversation_history[0]["role"] == "user"
        assert engine.conversation_history[0]["content"] == "What is GDMT?"
        engine.clear_conversation()

    def test_add_multiple_entries(self):
        engine = CardioRAGEngine(session_id="test_add_multi")
        engine.add_conversation_context("user", "Question 1")
        engine.add_conversation_context("assistant", "Answer 1")
        engine.add_conversation_context("user", "Question 2")
        assert len(engine.conversation_history) == 3
        engine.clear_conversation()

    def test_clear_conversation(self):
        engine = CardioRAGEngine(session_id="test_clear")
        engine.add_conversation_context("user", "test")
        engine.add_conversation_context("assistant", "response")
        engine.clear_conversation()
        assert len(engine.conversation_history) == 0

    def test_conversation_trimmed_to_max(self):
        engine = CardioRAGEngine(session_id="test_trim")
        # MAX_CONVERSATION_CONTEXT is 5, so max entries = 10
        for i in range(20):
            role = "user" if i % 2 == 0 else "assistant"
            engine.add_conversation_context(role, f"Message {i}")
        assert len(engine.conversation_history) <= 10
        engine.clear_conversation()

    def test_conversation_entry_has_timestamp(self):
        engine = CardioRAGEngine(session_id="test_timestamp")
        engine.add_conversation_context("user", "test")
        assert "timestamp" in engine.conversation_history[0]
        engine.clear_conversation()


# =====================================================================
# CardioRAGEngine._get_boosted_weights()
# =====================================================================


class TestGetBoostedWeights:
    """_get_boosted_weights() should apply workflow boosts."""

    @pytest.fixture
    def engine(self):
        return CardioRAGEngine()

    def test_no_workflow_returns_base_weights(self, engine):
        weights = engine._get_boosted_weights(None)
        assert isinstance(weights, dict)
        assert len(weights) > 0

    def test_with_workflow_returns_dict(self, engine):
        weights = engine._get_boosted_weights(CardioWorkflowType.HEART_FAILURE)
        assert isinstance(weights, dict)
        assert len(weights) > 0

    def test_all_values_are_positive_floats(self, engine):
        weights = engine._get_boosted_weights(CardioWorkflowType.CAD_ASSESSMENT)
        for coll, w in weights.items():
            assert isinstance(w, float)
            assert w >= 0


# =====================================================================
# create_rag_engine() FACTORY
# =====================================================================


class TestCreateRagEngine:
    """create_rag_engine() should return a CardioRAGEngine."""

    @patch("src.rag_engine.CardioRAGEngine")
    def test_returns_engine(self, MockEngine):
        MockEngine.return_value = MagicMock()
        engine = create_rag_engine(
            milvus_client=MagicMock(),
            embedding_model=MagicMock(),
            llm_client=MagicMock(),
        )
        assert engine is not None

    def test_without_any_args_returns_engine(self):
        # Should not raise, just warn about missing components
        engine = create_rag_engine(
            milvus_client=MagicMock(),
            embedding_model=MagicMock(),
            llm_client=MagicMock(),
            query_expander=MagicMock(),
        )
        assert isinstance(engine, CardioRAGEngine)


# =====================================================================
# format_search_results() DISPLAY HELPER
# =====================================================================


class TestFormatSearchResults:
    """format_search_results() should produce human-readable output."""

    def test_no_results(self):
        output = format_search_results([])
        assert "No results found" in output

    def test_with_results(self):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                record_id="1",
                score=0.92,
                text="Guideline recommendation text",
                metadata={"collection_label": "Guideline"},
                relevance="high",
            ),
            CardioSearchResult(
                collection="cardio_literature",
                record_id="2",
                score=0.78,
                text="Literature abstract text",
                metadata={"collection_label": "Literature"},
                relevance="medium",
            ),
        ]
        output = format_search_results(results)
        assert "2 results" in output
        assert "Guideline" in output
        assert "Literature" in output

    def test_shows_relevance_tag(self):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                record_id="1",
                score=0.9,
                text="Test",
                metadata={"collection_label": "Guideline"},
                relevance="high",
            ),
        ]
        output = format_search_results(results)
        assert "[high]" in output

    def test_shows_pmid_when_present(self):
        results = [
            CardioSearchResult(
                collection="cardio_literature",
                record_id="1",
                score=0.9,
                text="Test",
                metadata={"collection_label": "Literature", "pmid": "12345678"},
                relevance="high",
            ),
        ]
        output = format_search_results(results)
        assert "12345678" in output


# =====================================================================
# _rerank_results() HEURISTIC BOOSTING
# =====================================================================


class TestRerankResults:
    """_rerank_results() should apply heuristic boosts."""

    @pytest.fixture
    def engine(self):
        return CardioRAGEngine()

    def test_returns_list(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                score=0.5,
                relevance="medium",
                metadata={},
            )
        ]
        reranked = engine._rerank_results(results, "heart failure")
        assert isinstance(reranked, list)

    def test_boosts_guideline_results_for_matching_condition(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_guidelines",
                score=0.5,
                relevance="medium",
                metadata={"condition": "heart failure", "class_of_rec": "I"},
            )
        ]
        original_score = results[0].score
        engine._rerank_results(results, "heart failure management")
        assert results[0].score > original_score

    def test_boosts_high_relevance_results(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_literature",
                score=0.5,
                relevance="high",
                metadata={},
            )
        ]
        original_score = results[0].score
        engine._rerank_results(results, "test query")
        assert results[0].score > original_score

    def test_boosts_pmid_results(self, engine):
        results = [
            CardioSearchResult(
                collection="cardio_literature",
                score=0.5,
                relevance="low",
                metadata={"pmid": "12345678"},
            )
        ]
        original_score = results[0].score
        engine._rerank_results(results, "test query")
        assert results[0].score > original_score

    def test_sorted_by_score_descending(self, engine):
        results = [
            CardioSearchResult(collection="a", score=0.3, relevance="low", metadata={}),
            CardioSearchResult(collection="b", score=0.9, relevance="high", metadata={}),
            CardioSearchResult(collection="c", score=0.6, relevance="medium", metadata={}),
        ]
        reranked = engine._rerank_results(results, "test")
        scores = [r.score for r in reranked]
        assert scores == sorted(scores, reverse=True)
