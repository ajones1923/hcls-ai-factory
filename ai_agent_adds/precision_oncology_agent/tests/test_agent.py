"""
Tests for the OncoIntelligenceAgent.
=====================================
Validates SearchPlan creation, search planning, evidence evaluation,
and agent initialization with mock dependencies.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.agent import (
    KNOWN_CANCER_TYPES,
    KNOWN_GENES,
    OncoIntelligenceAgent,
    SearchPlan,
    _CANCER_ALIASES,
)


# ═══════════════════════════════════════════════════════════════════════════
# SearchPlan Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestSearchPlan:
    """Test SearchPlan dataclass creation with oncology-specific plans."""

    def test_basic_creation(self):
        plan = SearchPlan(question="What is BRAF V600E?")
        assert plan.question == "What is BRAF V600E?"
        assert plan.search_strategy == "broad"
        assert plan.identified_topics == []
        assert plan.target_genes == []
        assert plan.relevant_cancer_types == []
        assert plan.sub_questions == []

    def test_targeted_plan(self):
        plan = SearchPlan(
            question="EGFR L858R treatment in NSCLC",
            target_genes=["EGFR"],
            relevant_cancer_types=["NSCLC"],
            search_strategy="targeted",
            identified_topics=["targeted therapy"],
        )
        assert plan.search_strategy == "targeted"
        assert "EGFR" in plan.target_genes
        assert "NSCLC" in plan.relevant_cancer_types

    def test_comparative_plan(self):
        plan = SearchPlan(
            question="osimertinib vs erlotinib",
            search_strategy="comparative",
            target_genes=["EGFR"],
            identified_topics=["targeted therapy"],
        )
        assert plan.search_strategy == "comparative"

    def test_multi_gene_plan(self):
        plan = SearchPlan(
            question="EGFR and ALK testing in NSCLC",
            target_genes=["EGFR", "ALK"],
            relevant_cancer_types=["NSCLC"],
            search_strategy="broad",
        )
        assert len(plan.target_genes) == 2
        assert "ALK" in plan.target_genes


# ═══════════════════════════════════════════════════════════════════════════
# Known Vocabularies Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestKnownGenes:
    """Test KNOWN_GENES set has expected oncology genes."""

    @pytest.mark.parametrize("gene", [
        "BRAF", "EGFR", "ALK", "ROS1", "KRAS", "HER2",
        "NTRK", "RET", "MET", "FGFR", "PIK3CA", "IDH1",
        "IDH2", "BRCA1", "BRCA2", "TP53", "PTEN", "CDKN2A",
        "STK11", "ESR1",
    ])
    def test_gene_in_known_set(self, gene):
        assert gene in KNOWN_GENES

    def test_known_genes_count(self):
        assert len(KNOWN_GENES) >= 20


class TestKnownCancerTypes:
    """Test KNOWN_CANCER_TYPES set has expected cancer types."""

    @pytest.mark.parametrize("cancer_type", [
        "NSCLC", "BREAST", "MELANOMA", "COLORECTAL", "PANCREATIC",
        "OVARIAN", "PROSTATE", "GLIOBLASTOMA", "AML", "CML",
    ])
    def test_cancer_type_in_known_set(self, cancer_type):
        assert cancer_type in KNOWN_CANCER_TYPES


class TestCancerAliases:
    """Test _CANCER_ALIASES resolve correctly."""

    @pytest.mark.parametrize("alias,expected", [
        ("lung", "NSCLC"),
        ("colon", "COLORECTAL"),
        ("crc", "COLORECTAL"),
        ("melanoma", "MELANOMA"),
        ("kidney", "RENAL"),
        ("liver", "HEPATOCELLULAR"),
        ("gbm", "GLIOBLASTOMA"),
        ("stomach", "GASTRIC"),
    ])
    def test_alias_resolution(self, alias, expected):
        assert _CANCER_ALIASES[alias] == expected


# ═══════════════════════════════════════════════════════════════════════════
# OncoIntelligenceAgent Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestAgentInit:
    """Test agent initialization with mock dependencies."""

    def test_basic_init(self):
        mock_rag = MagicMock()
        agent = OncoIntelligenceAgent(rag_engine=mock_rag)
        assert agent.rag_engine is mock_rag

    def test_agent_constants(self):
        mock_rag = MagicMock()
        agent = OncoIntelligenceAgent(rag_engine=mock_rag)
        assert agent.MAX_RETRIES >= 1
        assert agent.MIN_SUFFICIENT_HITS >= 1
        assert agent.MIN_COLLECTIONS_FOR_SUFFICIENT >= 1


class TestSearchPlanGeneration:
    """Test search_plan() returns valid SearchPlan for oncology queries."""

    @pytest.fixture
    def agent(self):
        return OncoIntelligenceAgent(rag_engine=MagicMock())

    def test_plan_for_single_gene_query(self, agent):
        plan = agent.search_plan("EGFR mutations in NSCLC")
        assert "EGFR" in plan.target_genes
        assert "NSCLC" in plan.relevant_cancer_types
        assert plan.search_strategy == "targeted"

    def test_plan_for_broad_query(self, agent):
        plan = agent.search_plan("What are the latest advances in immunotherapy?")
        assert plan.search_strategy == "broad"
        assert "immunotherapy response" in plan.identified_topics

    def test_plan_for_comparative_query(self, agent):
        plan = agent.search_plan("compare EGFR vs ALK in NSCLC")
        assert plan.search_strategy == "comparative"
        assert "EGFR" in plan.target_genes
        assert "ALK" in plan.target_genes

    def test_plan_identifies_resistance_topic(self, agent):
        plan = agent.search_plan("EGFR resistance mechanisms in NSCLC")
        assert "therapeutic resistance" in plan.identified_topics

    def test_plan_identifies_biomarker_topic(self, agent):
        plan = agent.search_plan("biomarker testing for NSCLC")
        assert "biomarker identification" in plan.identified_topics

    def test_plan_identifies_clinical_trial_topic(self, agent):
        plan = agent.search_plan("clinical trial for BRAF melanoma")
        assert "clinical trials" in plan.identified_topics

    def test_plan_identifies_tmb_topic(self, agent):
        plan = agent.search_plan("tumor mutational burden and immunotherapy")
        assert "TMB" in plan.identified_topics

    def test_plan_identifies_msi_topic(self, agent):
        plan = agent.search_plan("microsatellite instability in colorectal cancer")
        assert "MSI / microsatellite instability" in plan.identified_topics

    def test_plan_identifies_ctdna_topic(self, agent):
        plan = agent.search_plan("ctdna monitoring in breast cancer")
        assert "liquid biopsy / ctDNA" in plan.identified_topics

    def test_plan_cancer_alias_resolution(self, agent):
        plan = agent.search_plan("EGFR in lung cancer")
        assert "NSCLC" in plan.relevant_cancer_types

    def test_plan_decomposition_multi_gene(self, agent):
        plan = agent.search_plan("EGFR and BRAF in NSCLC")
        assert len(plan.sub_questions) >= 2

    def test_plan_returns_searchplan_type(self, agent):
        plan = agent.search_plan("BRAF therapy")
        assert isinstance(plan, SearchPlan)


class TestEvaluateEvidence:
    """Test evaluate_evidence filters by score threshold."""

    @pytest.fixture
    def agent(self):
        return OncoIntelligenceAgent(rag_engine=MagicMock())

    def test_empty_evidence_is_insufficient(self, agent):
        assert agent.evaluate_evidence([]) == "insufficient"

    def test_few_hits_is_partial(self, agent):
        # Create a single evidence item without a collection attribute
        mock_evidence = MagicMock()
        mock_evidence.collection = "onco_variants"
        result = agent.evaluate_evidence([mock_evidence])
        assert result in ("partial", "insufficient")

    def test_sufficient_evidence(self, agent):
        # Create enough evidence items from multiple collections
        evidence_items = []
        for i, col in enumerate(["onco_variants", "onco_literature", "onco_therapies"]):
            item = MagicMock()
            item.collection = col
            evidence_items.append(item)
        # Need at least MIN_SUFFICIENT_HITS and MIN_COLLECTIONS_FOR_SUFFICIENT
        result = agent.evaluate_evidence(evidence_items)
        assert result in ("sufficient", "partial")


class TestFallbackQueries:
    """Test _generate_fallback_queries produces broader queries."""

    @pytest.fixture
    def agent(self):
        return OncoIntelligenceAgent(rag_engine=MagicMock())

    def test_fallback_with_genes(self, agent):
        plan = SearchPlan(
            question="EGFR therapy",
            target_genes=["EGFR"],
        )
        fallbacks = agent._generate_fallback_queries(plan)
        assert len(fallbacks) >= 1
        assert any("EGFR" in fb for fb in fallbacks)

    def test_fallback_with_cancer_types(self, agent):
        plan = SearchPlan(
            question="melanoma treatment",
            relevant_cancer_types=["MELANOMA"],
        )
        fallbacks = agent._generate_fallback_queries(plan)
        assert any("MELANOMA" in fb for fb in fallbacks)

    def test_fallback_generic(self, agent):
        plan = SearchPlan(question="latest oncology advances")
        fallbacks = agent._generate_fallback_queries(plan)
        assert len(fallbacks) >= 1
        assert any("precision oncology" in fb.lower() for fb in fallbacks)
