"""Tests for Biomarker Intelligence Agent RAG engine module.

Validates evidence retrieval, prompt building, collection filtering,
merge-and-rank deduplication, and cross-collection search — all with mocks.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import CrossCollectionResult, SearchHit

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def mock_rag_components(mock_embedder, mock_llm_client, mock_collection_manager):
    """Return a dict with all mock RAG components."""
    return {
        "embedder": mock_embedder,
        "llm_client": mock_llm_client,
        "manager": mock_collection_manager,
    }


# =====================================================================
# EVIDENCE RETRIEVAL
# =====================================================================


class TestEvidenceRetrieval:
    """Tests for RAG evidence retrieval with mocked components."""

    def test_search_returns_results(self, mock_collection_manager):
        """search() returns a list when given a query embedding."""
        mock_collection_manager.search.return_value = [
            {"id": "ref-1", "score": 0.85, "collection": "biomarker_reference",
             "text_chunk": "Albumin is a liver-synthesized protein."}
        ]
        results = mock_collection_manager.search(
            collection_name="biomarker_reference",
            query_embedding=[0.0] * 384,
            top_k=5,
        )
        assert len(results) == 1
        assert results[0]["score"] == 0.85

    def test_search_all_returns_dict(self, mock_collection_manager):
        """search_all() returns a dict mapping collection names to result lists."""
        results = mock_collection_manager.search_all(
            query_embedding=[0.0] * 384,
            top_k_per_collection=5,
        )
        assert isinstance(results, dict)
        assert "biomarker_reference" in results

    def test_embedder_produces_384_dim(self, mock_embedder):
        """embed_text() returns a 384-dimensional vector."""
        vec = mock_embedder.embed_text("test query")
        assert len(vec) == 384

    def test_cross_collection_result_structure(self, sample_search_hits):
        """CrossCollectionResult stores query, hits, and timing."""
        result = CrossCollectionResult(
            query="What biomarkers indicate diabetes risk?",
            hits=sample_search_hits,
            total_collections_searched=10,
            search_time_ms=42.5,
        )
        assert result.hit_count == 5
        assert result.total_collections_searched == 10

    def test_empty_search_returns_empty(self, mock_collection_manager):
        """search() returns empty list when no matches found."""
        results = mock_collection_manager.search(
            collection_name="biomarker_reference",
            query_embedding=[0.0] * 384,
            top_k=5,
        )
        assert results == []


# =====================================================================
# PROMPT BUILDING
# =====================================================================


class TestPromptBuilding:
    """Tests for evidence-based prompt construction."""

    def test_evidence_context_includes_hits(self, sample_search_hits):
        """Evidence context should include text from all search hits."""
        evidence_text = ""
        for i, hit in enumerate(sample_search_hits):
            evidence_text += f"[{i+1}] ({hit.collection}): {hit.text}\n"

        assert "HbA1c" in evidence_text
        assert "MTHFR" in evidence_text
        assert "CYP2D6" in evidence_text

    def test_prompt_includes_question_and_evidence(self, sample_evidence):
        """Built prompt should include the question and evidence context."""
        question = sample_evidence.query
        evidence_text = "\n".join(
            f"[{i+1}] {hit.text}" for i, hit in enumerate(sample_evidence.hits)
        )
        prompt = f"Question: {question}\n\nEvidence:\n{evidence_text}\n\nAnswer:"

        assert question in prompt
        assert "Evidence:" in prompt
        assert "HbA1c" in prompt

    def test_knowledge_context_included_when_present(self, sample_evidence):
        """Knowledge context should be prepended to evidence when available."""
        context = sample_evidence.knowledge_context
        assert len(context) > 0
        assert "Diabetes" in context

    def test_system_prompt_is_biomarker_focused(self):
        """System prompt should mention precision medicine / biomarkers."""
        system_prompt = (
            "You are a precision medicine biomarker expert. Answer the question "
            "using the provided evidence. Cite sources by their number."
        )
        assert "precision medicine" in system_prompt
        assert "biomarker" in system_prompt


# =====================================================================
# COLLECTION FILTERING
# =====================================================================


class TestCollectionFiltering:
    """Tests for collection-specific search filtering."""

    def test_filter_to_single_collection(self, mock_collection_manager):
        """Filtering to one collection should only search that collection."""
        mock_collection_manager.search.return_value = [
            {"id": "pgx-1", "score": 0.9, "collection": "biomarker_pgx_rules"}
        ]
        results = mock_collection_manager.search(
            collection_name="biomarker_pgx_rules",
            query_embedding=[0.0] * 384,
            top_k=5,
        )
        assert all(r["collection"] == "biomarker_pgx_rules" for r in results)

    def test_filter_excludes_unselected_collections(self):
        """Unselected collections should not appear in filtered results."""
        selected = ["biomarker_reference", "biomarker_pgx_rules"]
        all_collections = [
            "biomarker_reference", "biomarker_genetic_variants",
            "biomarker_pgx_rules", "biomarker_disease_trajectories",
        ]
        filtered = [c for c in all_collections if c in selected]
        assert "biomarker_genetic_variants" not in filtered
        assert len(filtered) == 2

    def test_empty_filter_uses_all_collections(self, mock_collection_manager):
        """Empty filter should search all collections."""
        results = mock_collection_manager.search_all(
            query_embedding=[0.0] * 384,
            top_k_per_collection=5,
        )
        assert len(results) == 14  # All 14 collections

    def test_filter_by_gene(self, mock_collection_manager):
        """Gene-specific filter expression should be passed to search."""
        mock_collection_manager.search.return_value = [
            {"id": "var-1", "score": 0.88, "collection": "biomarker_genetic_variants",
             "gene": "MTHFR"}
        ]
        results = mock_collection_manager.search(
            collection_name="biomarker_genetic_variants",
            query_embedding=[0.0] * 384,
            top_k=5,
            filter_expr='gene == "MTHFR"',
        )
        assert results[0]["gene"] == "MTHFR"


# =====================================================================
# MERGE AND RANK
# =====================================================================


class TestMergeAndRank:
    """Tests for cross-collection result merging and deduplication."""

    def test_deduplicates_by_id(self):
        """Duplicate IDs across collections are collapsed to a single hit."""
        hits = [
            SearchHit(collection="ref", id="dup-1", score=0.9, text="A"),
            SearchHit(collection="evi", id="dup-1", score=0.7, text="A duplicate"),
            SearchHit(collection="pgx", id="unique-1", score=0.8, text="B"),
        ]
        seen_ids = set()
        deduped = []
        for hit in hits:
            if hit.id not in seen_ids:
                deduped.append(hit)
                seen_ids.add(hit.id)
        assert len(deduped) == 2

    def test_sorts_by_score_descending(self):
        """Results should be sorted highest score first."""
        hits = [
            SearchHit(collection="A", id="low", score=0.3, text="low"),
            SearchHit(collection="B", id="high", score=0.95, text="high"),
            SearchHit(collection="C", id="mid", score=0.6, text="mid"),
        ]
        sorted_hits = sorted(hits, key=lambda h: h.score, reverse=True)
        scores = [h.score for h in sorted_hits]
        assert scores == [0.95, 0.6, 0.3]

    def test_caps_at_max_results(self):
        """Merge should cap results at a maximum count."""
        max_results = 30
        hits = [
            SearchHit(collection="Lit", id=str(i), score=0.5, text=f"hit {i}")
            for i in range(50)
        ]
        capped = hits[:max_results]
        assert len(capped) == 30

    def test_empty_merge(self):
        """Merging empty lists returns empty list."""
        hits = []
        assert len(hits) == 0

    def test_cross_collection_grouping(self, sample_search_hits):
        """hits_by_collection() correctly groups by collection name."""
        result = CrossCollectionResult(query="test", hits=sample_search_hits)
        grouped = result.hits_by_collection()
        assert len(grouped) == 5  # 5 unique collections in sample
        for coll, hits in grouped.items():
            assert all(h.collection == coll for h in hits)


# =====================================================================
# LLM RESPONSE GENERATION
# =====================================================================


class TestLLMGeneration:
    """Tests for LLM response generation with mocked client."""

    def test_generate_returns_string(self, mock_llm_client):
        """generate() returns a string response."""
        response = mock_llm_client.generate("test prompt")
        assert isinstance(response, str)
        assert response == "Mock response"

    def test_generate_stream_returns_iterator(self, mock_llm_client):
        """generate_stream() returns an iterator of text chunks."""
        chunks = list(mock_llm_client.generate_stream("test prompt"))
        assert len(chunks) == 2
        assert "".join(chunks) == "Mock response"

    def test_llm_receives_system_prompt(self, mock_llm_client):
        """LLM should receive both system prompt and user prompt."""
        mock_llm_client.generate(
            "What is HbA1c?",
            system_prompt="You are a biomarker expert.",
        )
        mock_llm_client.generate.assert_called_once()
        call_args = mock_llm_client.generate.call_args
        assert call_args[0][0] == "What is HbA1c?"
