"""Tests for Imaging Intelligence Agent RAG engine.

Validates ImagingRAGEngine initialization, retrieve(), query(),
retrieve_comparative(), find_related(), prompt building, and
knowledge context augmentation -- all using mock dependencies.

Author: Adam Jones
Date: February 2026
"""

import pytest
from unittest.mock import MagicMock, patch, call

from src.models import (
    AgentResponse,
    ComparativeResult,
    CrossCollectionResult,
    SearchHit,
)
from src.rag_engine import (
    COLLECTION_CONFIG,
    SYSTEM_PROMPT,
    ImagingRAGEngine,
)


# ===================================================================
# INITIALIZATION
# ===================================================================


class TestImagingRAGEngineInit:
    """Tests for ImagingRAGEngine initialization."""

    def test_init_stores_dependencies(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        assert engine.collection_manager is mock_collection_manager
        assert engine.embedder is mock_embedder
        assert engine.llm_client is mock_llm_client
        assert engine.nim_manager is None

    def test_init_with_nim_manager(self, mock_collection_manager, mock_embedder, mock_llm_client, mock_nim_services):
        engine = ImagingRAGEngine(
            mock_collection_manager, mock_embedder, mock_llm_client,
            nim_service_manager=mock_nim_services,
        )
        assert engine.nim_manager is mock_nim_services

    def test_system_prompt_is_set(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        assert engine.system_prompt == SYSTEM_PROMPT
        assert "medical imaging" in engine.system_prompt.lower()


# ===================================================================
# COLLECTION CONFIG
# ===================================================================


class TestCollectionConfig:
    """Tests for the COLLECTION_CONFIG constant."""

    def test_has_11_collections(self):
        assert len(COLLECTION_CONFIG) == 11

    def test_all_collections_have_weight(self):
        for name, config in COLLECTION_CONFIG.items():
            assert "weight" in config, f"{name} missing 'weight'"
            assert 0.0 < config["weight"] <= 1.0

    def test_all_collections_have_label(self):
        for name, config in COLLECTION_CONFIG.items():
            assert "label" in config, f"{name} missing 'label'"

    def test_weights_sum_approximately_to_one(self):
        total = sum(cfg["weight"] for cfg in COLLECTION_CONFIG.values())
        assert 0.9 <= total <= 1.1, f"Weights sum to {total}, expected ~1.0"

    def test_imaging_literature_weight(self):
        assert COLLECTION_CONFIG["imaging_literature"]["weight"] == 0.18

    def test_genomic_evidence_included(self):
        assert "genomic_evidence" in COLLECTION_CONFIG


# ===================================================================
# _embed_query
# ===================================================================


class TestEmbedQuery:
    """Tests for ImagingRAGEngine._embed_query()."""

    def test_embed_query_calls_embedder(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine._embed_query("test query")
        mock_embedder.encode.assert_called_once()
        assert isinstance(result, list)

    def test_embed_query_returns_list(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine._embed_query("hemorrhage CT")
        assert isinstance(result, list)
        assert len(result) == 384


# ===================================================================
# _get_knowledge_context
# ===================================================================


class TestGetKnowledgeContext:
    """Tests for ImagingRAGEngine._get_knowledge_context()."""

    def test_returns_pathology_context_for_hemorrhage(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        ctx = engine._get_knowledge_context("intracranial hemorrhage on CT")
        assert len(ctx) > 0
        assert "Intracranial Hemorrhage" in ctx

    def test_returns_modality_context_for_ct(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        ctx = engine._get_knowledge_context("ct scan findings")
        assert len(ctx) > 0

    def test_returns_anatomy_context_for_brain(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        ctx = engine._get_knowledge_context("brain mri lesion")
        assert len(ctx) > 0

    def test_returns_empty_for_unknown_query(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        ctx = engine._get_knowledge_context("xyzzy foobar qwerty")
        assert ctx == ""


# ===================================================================
# _is_comparative
# ===================================================================


class TestIsComparative:
    """Tests for ImagingRAGEngine._is_comparative()."""

    def test_vs_is_comparative(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        assert engine._is_comparative("CT vs MRI for hemorrhage") is True

    def test_versus_is_comparative(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        assert engine._is_comparative("CT versus MRI") is True

    def test_compared_to_is_comparative(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        assert engine._is_comparative("VISTA-3D compared to nnU-Net") is True

    def test_normal_query_is_not_comparative(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        assert engine._is_comparative("What is VISTA-3D?") is False


# ===================================================================
# _build_prompt
# ===================================================================


class TestBuildPrompt:
    """Tests for ImagingRAGEngine._build_prompt()."""

    def test_returns_list_of_messages(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        evidence = CrossCollectionResult(query="test", hits=[])
        messages = engine._build_prompt("test question", evidence)
        assert isinstance(messages, list)
        assert len(messages) == 2

    def test_first_message_is_system(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        evidence = CrossCollectionResult(query="test", hits=[])
        messages = engine._build_prompt("test question", evidence)
        assert messages[0]["role"] == "system"
        assert SYSTEM_PROMPT in messages[0]["content"]

    def test_second_message_is_user(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        evidence = CrossCollectionResult(query="test", hits=[])
        messages = engine._build_prompt("what is VISTA-3D?", evidence)
        assert messages[1]["role"] == "user"
        assert "what is VISTA-3D?" in messages[1]["content"]

    def test_includes_evidence_text(self, mock_collection_manager, mock_embedder, mock_llm_client, sample_search_hits):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        evidence = CrossCollectionResult(
            query="test",
            hits=sample_search_hits,
            total_collections_searched=5,
        )
        messages = engine._build_prompt("hemorrhage AI", evidence)
        user_content = messages[1]["content"]
        assert "Retrieved Evidence" in user_content
        assert "5 results" in user_content

    def test_includes_knowledge_context(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        evidence = CrossCollectionResult(
            query="test", hits=[],
            knowledge_context="## Pathology: Hemorrhage\nImportant info.",
        )
        messages = engine._build_prompt("test", evidence)
        user_content = messages[1]["content"]
        assert "Domain Knowledge" in user_content
        assert "Hemorrhage" in user_content

    def test_includes_conversation_context(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        evidence = CrossCollectionResult(query="test", hits=[])
        messages = engine._build_prompt(
            "follow-up question", evidence,
            conversation_context="Previous Q: What is VISTA-3D?",
        )
        user_content = messages[1]["content"]
        assert "Conversation Context" in user_content
        assert "Previous Q" in user_content


# ===================================================================
# retrieve()
# ===================================================================


class TestRetrieve:
    """Tests for ImagingRAGEngine.retrieve()."""

    def test_retrieve_calls_embedder(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("hemorrhage detection")
        mock_embedder.encode.assert_called_once()

    def test_retrieve_searches_all_collections(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("brain tumor segmentation")
        assert mock_collection_manager.search.call_count == 11  # all collections

    def test_retrieve_returns_cross_collection_result(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("lung nodule")
        assert isinstance(result, CrossCollectionResult)
        assert result.query == "lung nodule"

    def test_retrieve_respects_collections_filter(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve(
            "test",
            collections_filter=["imaging_literature", "imaging_trials"],
        )
        assert mock_collection_manager.search.call_count == 2

    def test_retrieve_includes_knowledge_context(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("intracranial hemorrhage detection")
        # Should find hemorrhage in knowledge and add context
        assert len(result.knowledge_context) > 0

    def test_retrieve_reports_collections_searched(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("test")
        assert result.total_collections_searched == 11

    def test_retrieve_has_positive_search_time(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("test")
        assert result.search_time_ms >= 0

    def test_retrieve_handles_search_failure(self, mock_collection_manager, mock_embedder, mock_llm_client):
        mock_collection_manager.search.side_effect = Exception("Connection failed")
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("test")
        # Should not crash; returns empty results
        assert isinstance(result, CrossCollectionResult)
        assert result.hit_count == 0

    def test_retrieve_hits_are_sorted_by_score(self, mock_collection_manager, mock_embedder, mock_llm_client):
        mock_collection_manager.search.return_value = [
            {"id": "low", "score": 0.3, "text_chunk": "low score"},
            {"id": "high", "score": 0.9, "text_chunk": "high score"},
        ]
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve("test", collections_filter=["imaging_literature"])
        if result.hit_count >= 2:
            assert result.hits[0].score >= result.hits[1].score


# ===================================================================
# query()
# ===================================================================


class TestQuery:
    """Tests for ImagingRAGEngine.query() (full RAG pipeline)."""

    def test_query_returns_string(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        answer = engine.query("What is VISTA-3D?")
        assert isinstance(answer, str)
        assert len(answer) > 0

    def test_query_calls_llm_generate(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        engine.query("hemorrhage detection models")
        mock_llm_client.generate.assert_called_once()

    def test_query_passes_messages_to_llm(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        engine.query("lung nodule classification")
        args, _ = mock_llm_client.generate.call_args
        messages = args[0]
        assert isinstance(messages, list)
        assert messages[0]["role"] == "system"

    def test_query_with_conversation_context(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        engine.query(
            "What about sensitivity?",
            conversation_context="Q: What is VISTA-3D?",
        )
        mock_llm_client.generate.assert_called_once()


# ===================================================================
# retrieve_comparative()
# ===================================================================


class TestRetrieveComparative:
    """Tests for ImagingRAGEngine.retrieve_comparative()."""

    def test_returns_comparative_result(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve_comparative("CT vs MRI for hemorrhage")
        assert isinstance(result, ComparativeResult)

    def test_extracts_entities(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve_comparative("CT vs MRI for brain imaging")
        assert result.entity_a != ""
        assert result.entity_b != ""

    def test_searches_both_entities(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve_comparative("VISTA-3D versus nnU-Net")
        # Two retrieve calls -> 2 x 11 collections = 22 searches
        assert mock_collection_manager.search.call_count == 22

    def test_includes_comparison_context(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve_comparative("ct vs mri")
        # ct and mri are both resolvable entities
        assert len(result.comparison_context) > 0

    def test_total_hits_from_both_sides(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.retrieve_comparative("CT vs MRI")
        assert result.total_hits == result.evidence_a.hit_count + result.evidence_b.hit_count


# ===================================================================
# _parse_comparison_entities
# ===================================================================


class TestParseComparisonEntities:
    """Tests for ImagingRAGEngine._parse_comparison_entities()."""

    def test_parses_vs(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        a, b = engine._parse_comparison_entities("CT vs MRI")
        assert a == "ct"
        assert b == "mri"

    def test_parses_versus(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        a, b = engine._parse_comparison_entities("VISTA-3D versus nnU-Net")
        assert "vista" in a
        assert "nnu-net" in b

    def test_parses_compared_to(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        a, b = engine._parse_comparison_entities("DenseNet compared to ResNet")
        assert "densenet" in a
        assert "resnet" in b

    def test_no_separator_returns_full_query(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        a, b = engine._parse_comparison_entities("What is VISTA-3D?")
        assert a == "What is VISTA-3D?"
        assert b == ""


# ===================================================================
# find_related()
# ===================================================================


class TestFindRelated:
    """Tests for ImagingRAGEngine.find_related()."""

    def test_returns_dict(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.find_related("VISTA-3D")
        assert isinstance(result, dict)

    def test_searches_all_collections(self, mock_collection_manager, mock_embedder, mock_llm_client):
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        engine.find_related("hemorrhage")
        assert mock_collection_manager.search.call_count == 11

    def test_handles_search_errors(self, mock_collection_manager, mock_embedder, mock_llm_client):
        mock_collection_manager.search.side_effect = Exception("fail")
        engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
        result = engine.find_related("test")
        assert isinstance(result, dict)
        assert len(result) == 0
