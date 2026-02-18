"""
Tests for RAG Engine module.
"""
import sys
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.rag_engine import (
    GENOMICS_RAG_SYSTEM_PROMPT,
    NEURODEGENERATION_GENES,
    PHARMACOGENOMIC_GENES,
    RAGEngine,
)


class TestRAGEngineQueryExpansion:
    """Tests for query expansion functionality."""

    def test_neurodegeneration_expansion_ftd(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test FTD-related queries expand to relevant genes."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("What genes are associated with frontotemporal dementia?")

        assert "VCP" in genes
        assert "C9orf72" in genes
        assert "GRN" in genes
        assert "MAPT" in genes

    def test_neurodegeneration_expansion_vcp(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test VCP queries expand correctly."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("Tell me about VCP mutations")

        assert "VCP" in genes

    def test_neurodegeneration_expansion_p97(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test p97 (VCP alternate name) queries expand to VCP."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("What about p97 protein dysfunction?")

        assert "VCP" in genes

    def test_pharmacogenomic_expansion_ppi(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test PPI drug queries expand to CYP genes."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("How will I respond to omeprazole?")

        assert "CYP2C19" in genes

    def test_pharmacogenomic_expansion_warfarin(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test warfarin queries expand correctly."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("What about warfarin dosing?")

        assert "CYP2C9" in genes
        assert "VKORC1" in genes

    def test_combined_expansion(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test queries that match multiple categories."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("dementia drug metabolism")

        # Should include both neurodegeneration and pharmacogenomic genes
        assert "VCP" in genes  # neurodegeneration
        assert "CYP2D6" in genes  # drug metabolism

    def test_no_expansion_for_unrelated(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test unrelated queries don't expand."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        genes = engine._get_expanded_genes("What is DNA?")

        assert len(genes) == 0


class TestRAGEngineRetrieval:
    """Tests for evidence retrieval."""

    def test_retrieve_basic(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test basic retrieval."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        results = engine.retrieve("Tell me about variants in VCP")

        # Should have called embedder and milvus
        mock_embedder.embed_query.assert_called_once()
        mock_milvus_client.search.assert_called_once()

        assert len(results) >= 0

    def test_retrieve_with_filter(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test retrieval with filter expression."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        results = engine.retrieve(
            "What variants are on chromosome 9?",
            filter_expr='chrom == "chr9"'
        )

        # Should pass filter to milvus
        call_args = mock_milvus_client.search.call_args
        assert call_args[1]["filter_expr"] == 'chrom == "chr9"'

    def test_retrieve_custom_top_k(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test retrieval with custom top_k."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            top_k=5,
        )

        results = engine.retrieve("variants", top_k=20)

        call_args = mock_milvus_client.search.call_args
        assert call_args[1]["top_k"] == 20


class TestRAGEngineQuery:
    """Tests for RAG query (retrieve + generate)."""

    def test_query_basic(self, mock_milvus_client, mock_embedder, mock_llm_client, sample_evidence_results):
        """Test basic query."""
        mock_milvus_client.search.return_value = sample_evidence_results

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        result = engine.query("What VCP variants were found?")

        assert "answer" in result
        assert "query" in result
        assert result["query"] == "What VCP variants were found?"

        # LLM should have been called
        mock_llm_client.generate.assert_called_once()

    def test_query_includes_evidence(self, mock_milvus_client, mock_embedder, mock_llm_client, sample_evidence_results):
        """Test query returns evidence when requested."""
        mock_milvus_client.search.return_value = sample_evidence_results

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        result = engine.query("What VCP variants were found?", include_evidence=True)

        assert "evidence" in result
        assert "evidence_count" in result
        assert result["evidence_count"] == len(sample_evidence_results)

    def test_query_without_evidence(self, mock_milvus_client, mock_embedder, mock_llm_client, sample_evidence_results):
        """Test query can exclude evidence from response."""
        mock_milvus_client.search.return_value = sample_evidence_results

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        result = engine.query("What VCP variants were found?", include_evidence=False)

        assert "evidence" not in result


class TestRAGEngineFormatting:
    """Tests for context formatting."""

    def test_format_evidence_context(self, mock_milvus_client, mock_embedder, mock_llm_client, sample_evidence_results):
        """Test evidence formatting for LLM context."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        context = engine._format_evidence_context(sample_evidence_results)

        assert "Evidence 1" in context
        assert "Evidence 2" in context
        assert "VCP" in context
        assert "similarity:" in context

    def test_format_empty_evidence(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test formatting when no evidence found."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        context = engine._format_evidence_context([])

        assert "No relevant evidence" in context

    def test_build_prompt(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test prompt building."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        prompt = engine._build_prompt(
            query="What VCP variants exist?",
            context="## Evidence\nVariant at chr9..."
        )

        assert "What VCP variants exist?" in prompt
        assert "## Evidence" in prompt
        assert "Instructions" in prompt


class TestRAGEngineGeneSearch:
    """Tests for gene-specific search methods."""

    def test_search_gene(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test searching by gene."""
        mock_milvus_client.search_by_gene.return_value = [
            {"id": "v1", "gene": "VCP"},
            {"id": "v2", "gene": "VCP"},
        ]

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        results = engine.search_gene("VCP")

        mock_milvus_client.search_by_gene.assert_called_with("VCP")
        assert len(results) == 2

    def test_search_region(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test searching by genomic region."""
        mock_milvus_client.search_by_region.return_value = [
            {"id": "v1", "chrom": "chr9", "pos": 123456},
        ]

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        results = engine.search_region("chr9", 100000, 200000)

        mock_milvus_client.search_by_region.assert_called_with("chr9", 100000, 200000)


class TestRAGEngineStreaming:
    """Tests for streaming functionality."""

    def test_query_stream(self, mock_milvus_client, mock_embedder, mock_llm_client, sample_evidence_results):
        """Test streaming query."""
        mock_milvus_client.search.return_value = sample_evidence_results

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        stream_results = list(engine.query_stream("What about VCP?"))

        # Should yield evidence first, then tokens, then done
        assert stream_results[0]["type"] == "evidence"
        assert stream_results[0]["content"] == sample_evidence_results

        # Should have token events
        token_events = [r for r in stream_results if r["type"] == "token"]
        assert len(token_events) > 0

        # Should end with done
        assert stream_results[-1]["type"] == "done"


class TestRAGEngineSystemPrompt:
    """Tests for system prompt configuration."""

    def test_default_system_prompt(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test default system prompt is used."""
        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
        )

        assert engine.system_prompt == GENOMICS_RAG_SYSTEM_PROMPT

    def test_custom_system_prompt(self, mock_milvus_client, mock_embedder, mock_llm_client):
        """Test custom system prompt can be provided."""
        custom_prompt = "You are a custom genomics assistant."

        engine = RAGEngine(
            milvus_client=mock_milvus_client,
            embedder=mock_embedder,
            llm_client=mock_llm_client,
            system_prompt=custom_prompt,
        )

        assert engine.system_prompt == custom_prompt

    def test_system_prompt_contains_vcp_guidance(self):
        """Test default prompt has VCP/FTD guidance."""
        assert "VCP" in GENOMICS_RAG_SYSTEM_PROMPT
        assert "p97" in GENOMICS_RAG_SYSTEM_PROMPT
        assert "frontotemporal" in GENOMICS_RAG_SYSTEM_PROMPT.lower()
        assert "IBMPFD" in GENOMICS_RAG_SYSTEM_PROMPT

    def test_system_prompt_contains_pharmacogenomics(self):
        """Test default prompt has pharmacogenomic guidance."""
        assert "CYP2C19" in GENOMICS_RAG_SYSTEM_PROMPT
        assert "CYP2D6" in GENOMICS_RAG_SYSTEM_PROMPT
        assert "pharmacogenomic" in GENOMICS_RAG_SYSTEM_PROMPT.lower()


class TestQueryExpansionMappings:
    """Tests for the gene expansion mappings."""

    def test_neurodegeneration_genes_complete(self):
        """Test neurodegeneration gene mappings are complete."""
        # FTD genes
        assert "VCP" in NEURODEGENERATION_GENES.get("ftd", [])
        assert "C9orf72" in NEURODEGENERATION_GENES.get("ftd", [])
        assert "GRN" in NEURODEGENERATION_GENES.get("ftd", [])
        assert "MAPT" in NEURODEGENERATION_GENES.get("ftd", [])

        # VCP aliases
        assert "VCP" in NEURODEGENERATION_GENES.get("vcp", [])
        assert "VCP" in NEURODEGENERATION_GENES.get("p97", [])
        assert "VCP" in NEURODEGENERATION_GENES.get("ibmpfd", [])

    def test_pharmacogenomic_genes_complete(self):
        """Test pharmacogenomic gene mappings are complete."""
        # PPIs
        assert "CYP2C19" in PHARMACOGENOMIC_GENES.get("ppi", [])
        assert "CYP2C19" in PHARMACOGENOMIC_GENES.get("omeprazole", [])

        # Warfarin
        assert "CYP2C9" in PHARMACOGENOMIC_GENES.get("warfarin", [])
        assert "VKORC1" in PHARMACOGENOMIC_GENES.get("warfarin", [])

        # Codeine
        assert "CYP2D6" in PHARMACOGENOMIC_GENES.get("codeine", [])
