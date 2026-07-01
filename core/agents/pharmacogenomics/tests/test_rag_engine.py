"""Tests for src/rag_engine.py — PGxRAGEngine.

Tests system prompt, collection config, workflow classification,
boosted weight computation, and known gene set.

Author: Adam Jones
Date: March 2026
"""

import pytest
from unittest.mock import MagicMock

from src.rag_engine import (
    PGX_SYSTEM_PROMPT,
    COLLECTION_CONFIG,
    WORKFLOW_COLLECTION_BOOST,
    _KNOWN_GENES,
    PGxRAGEngine,
)
from src.models import PGxWorkflowType


# ═══════════════════════════════════════════════════════════════════════
# SYSTEM PROMPT
# ═══════════════════════════════════════════════════════════════════════

class TestSystemPrompt:
    def test_not_empty(self):
        assert len(PGX_SYSTEM_PROMPT) > 500

    def test_mentions_pharmacogenomics(self):
        assert "pharmacogenomics" in PGX_SYSTEM_PROMPT.lower()

    def test_mentions_star_allele(self):
        assert "star allele" in PGX_SYSTEM_PROMPT.lower()

    def test_mentions_cpic(self):
        assert "CPIC" in PGX_SYSTEM_PROMPT

    def test_mentions_hla(self):
        assert "HLA" in PGX_SYSTEM_PROMPT

    def test_mentions_phenoconversion(self):
        assert "phenoconversion" in PGX_SYSTEM_PROMPT.lower()

    def test_mentions_dosing(self):
        assert "dosing" in PGX_SYSTEM_PROMPT.lower()


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION CONFIG
# ═══════════════════════════════════════════════════════════════════════

class TestCollectionConfig:
    def test_has_15_entries(self):
        assert len(COLLECTION_CONFIG) == 15

    @pytest.mark.parametrize("name", [
        "pgx_gene_reference", "pgx_drug_guidelines", "pgx_drug_interactions",
        "pgx_hla_hypersensitivity", "pgx_phenoconversion", "pgx_dosing_algorithms",
        "pgx_clinical_evidence", "pgx_population_data", "pgx_clinical_trials",
        "pgx_fda_labels", "pgx_drug_alternatives", "pgx_patient_profiles",
        "pgx_implementation", "pgx_education", "genomic_evidence",
    ])
    def test_collection_present(self, name):
        assert name in COLLECTION_CONFIG

    def test_each_has_weight(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "weight" in cfg
            assert isinstance(cfg["weight"], float)
            assert cfg["weight"] > 0

    def test_each_has_label(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "label" in cfg
            assert isinstance(cfg["label"], str)

    def test_weights_sum_approximately_one(self):
        total = sum(cfg["weight"] for cfg in COLLECTION_CONFIG.values())
        assert abs(total - 1.0) < 0.05

    def test_has_gene_and_drug_field_flags(self):
        for name, cfg in COLLECTION_CONFIG.items():
            assert "has_gene_field" in cfg
            assert "has_drug_field" in cfg

    def test_drug_guidelines_has_both_fields(self):
        cfg = COLLECTION_CONFIG["pgx_drug_guidelines"]
        assert cfg["has_gene_field"] is True
        assert cfg["has_drug_field"] is True


# ═══════════════════════════════════════════════════════════════════════
# WORKFLOW COLLECTION BOOST
# ═══════════════════════════════════════════════════════════════════════

class TestWorkflowCollectionBoost:
    def test_has_6_workflow_types(self):
        assert len(WORKFLOW_COLLECTION_BOOST) == 6

    @pytest.mark.parametrize("wf", [
        PGxWorkflowType.GENE_QUERY,
        PGxWorkflowType.DRUG_QUERY,
        PGxWorkflowType.PROFILE_QUERY,
        PGxWorkflowType.INTERACTION_QUERY,
        PGxWorkflowType.DOSING_QUERY,
        PGxWorkflowType.HLA_SCREEN,
    ])
    def test_workflow_present(self, wf):
        assert wf in WORKFLOW_COLLECTION_BOOST

    def test_gene_query_boosts_gene_reference(self):
        collections = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.GENE_QUERY]
        assert "pgx_gene_reference" in collections

    def test_drug_query_boosts_guidelines(self):
        collections = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.DRUG_QUERY]
        assert "pgx_drug_guidelines" in collections

    def test_hla_screen_boosts_hypersensitivity(self):
        collections = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.HLA_SCREEN]
        assert "pgx_hla_hypersensitivity" in collections

    def test_interaction_boosts_phenoconversion(self):
        collections = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.INTERACTION_QUERY]
        assert "pgx_phenoconversion" in collections

    def test_dosing_boosts_algorithms(self):
        collections = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.DOSING_QUERY]
        assert "pgx_dosing_algorithms" in collections


# ═══════════════════════════════════════════════════════════════════════
# KNOWN GENES
# ═══════════════════════════════════════════════════════════════════════

class TestKnownGenes:
    def test_not_empty(self):
        assert len(_KNOWN_GENES) > 20

    @pytest.mark.parametrize("gene", [
        "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "DPYD",
        "TPMT", "NUDT15", "UGT1A1", "SLCO1B1", "VKORC1",
        "G6PD", "HLA-A", "HLA-B",
    ])
    def test_gene_present(self, gene):
        assert gene in _KNOWN_GENES


# ═══════════════════════════════════════════════════════════════════════
# PGxRAGEngine
# ═══════════════════════════════════════════════════════════════════════

class TestPGxRAGEngine:
    @pytest.fixture
    def engine(self):
        mock_cm = MagicMock()
        mock_embedder = MagicMock()
        mock_llm = MagicMock()
        return PGxRAGEngine(mock_cm, mock_embedder, mock_llm)

    def test_init(self, engine):
        assert engine.collections is not None
        assert engine.embedder is not None
        assert engine.llm is not None
        assert engine.knowledge is None
        assert engine.expander is None

    def test_init_with_optional(self):
        mock_cm = MagicMock()
        mock_emb = MagicMock()
        mock_llm = MagicMock()
        mock_know = MagicMock()
        mock_exp = MagicMock()
        engine = PGxRAGEngine(mock_cm, mock_emb, mock_llm,
                              knowledge=mock_know, query_expander=mock_exp)
        assert engine.knowledge is mock_know
        assert engine.expander is mock_exp

    # -- _classify_workflow --

    def test_classify_gene_query(self, engine):
        wfs = engine._classify_workflow("Tell me about CYP2D6")
        assert PGxWorkflowType.GENE_QUERY in wfs

    def test_classify_hla_screen(self, engine):
        wfs = engine._classify_workflow("Should I screen for HLA-B*57:01 before abacavir?")
        assert PGxWorkflowType.HLA_SCREEN in wfs

    def test_classify_dosing(self, engine):
        wfs = engine._classify_workflow("What dose of warfarin should I use?")
        assert PGxWorkflowType.DOSING_QUERY in wfs

    def test_classify_interaction(self, engine):
        wfs = engine._classify_workflow("What is the phenoconversion risk with fluoxetine?")
        assert PGxWorkflowType.INTERACTION_QUERY in wfs

    def test_classify_profile(self, engine):
        wfs = engine._classify_workflow("Give me a comprehensive patient profile panel")
        assert PGxWorkflowType.PROFILE_QUERY in wfs

    def test_classify_star_allele_pattern(self, engine):
        wfs = engine._classify_workflow("What does *4/*4 mean?")
        assert PGxWorkflowType.GENE_QUERY in wfs

    def test_classify_default_when_empty(self, engine):
        wfs = engine._classify_workflow("hello world")
        assert len(wfs) >= 1  # Should default to PROFILE_QUERY

    # -- _compute_boosted_weights --

    def test_boosted_weights_sum_to_one(self, engine):
        wfs = [PGxWorkflowType.GENE_QUERY]
        weights = engine._compute_boosted_weights(wfs)
        total = sum(weights.values())
        assert abs(total - 1.0) < 0.001

    def test_boosted_weights_increase_target(self, engine):
        base_weights = engine._compute_boosted_weights([])
        boosted_weights = engine._compute_boosted_weights([PGxWorkflowType.HLA_SCREEN])
        # HLA collection should have higher relative weight
        assert boosted_weights.get("pgx_hla_hypersensitivity", 0) >= \
               base_weights.get("pgx_hla_hypersensitivity", 0)

    def test_boosted_weights_has_all_collections(self, engine):
        weights = engine._compute_boosted_weights([PGxWorkflowType.GENE_QUERY])
        assert len(weights) == 15
