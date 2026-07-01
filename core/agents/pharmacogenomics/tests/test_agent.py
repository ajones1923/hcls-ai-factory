"""Tests for src/agent.py — PGxIntelligenceAgent.

Tests PHARMACOGENES, HLA_ALLELES, DRUG_CATEGORIES knowledge dicts,
SearchPlan dataclass, search_plan method, and workflow boost mapping.

Author: Adam Jones
Date: March 2026
"""

import pytest
from unittest.mock import MagicMock

from src.agent import (
    PHARMACOGENES,
    HLA_ALLELES,
    DRUG_CATEGORIES,
    PGX_SYSTEM_PROMPT,
    WORKFLOW_COLLECTION_BOOST,
    SearchPlan,
    PGxIntelligenceAgent,
)
from src.models import PGxWorkflowType


# ═══════════════════════════════════════════════════════════════════════
# PHARMACOGENES
# ═══════════════════════════════════════════════════════════════════════

class TestPharmacogenes:
    def test_has_25_genes(self):
        assert len(PHARMACOGENES) == 25

    @pytest.mark.parametrize("gene", [
        "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "CYP2B6",
        "DPYD", "TPMT", "NUDT15", "UGT1A1", "SLCO1B1",
        "VKORC1", "CYP4F2", "G6PD", "IFNL3", "RYR1", "CACNA1S",
        "ABCB1", "CFTR", "CYP1A2", "CYP2C8", "CYP3A4",
        "F5", "HLA", "MTHFR", "NAT2",
    ])
    def test_gene_present(self, gene):
        assert gene in PHARMACOGENES

    def test_each_has_cpic_guidelines(self):
        for gene, info in PHARMACOGENES.items():
            assert "cpic_guidelines" in info
            assert isinstance(info["cpic_guidelines"], list)

    def test_most_have_guidelines(self):
        with_guidelines = [g for g, i in PHARMACOGENES.items() if len(i["cpic_guidelines"]) > 0]
        assert len(with_guidelines) >= 20

    def test_each_has_function(self):
        for gene, info in PHARMACOGENES.items():
            assert "function" in info
            assert len(info["function"]) > 10

    def test_cyp2d6_guidelines_include_codeine(self):
        assert "codeine" in PHARMACOGENES["CYP2D6"]["cpic_guidelines"]

    def test_cyp2c19_guidelines_include_clopidogrel(self):
        assert "clopidogrel" in PHARMACOGENES["CYP2C19"]["cpic_guidelines"]

    def test_dpyd_guidelines_include_fluorouracil(self):
        assert "fluorouracil" in PHARMACOGENES["DPYD"]["cpic_guidelines"]

    def test_slco1b1_guidelines_include_simvastatin(self):
        assert "simvastatin" in PHARMACOGENES["SLCO1B1"]["cpic_guidelines"]


# ═══════════════════════════════════════════════════════════════════════
# HLA_ALLELES
# ═══════════════════════════════════════════════════════════════════════

class TestHLAAlleles:
    def test_has_9_alleles(self):
        assert len(HLA_ALLELES) == 9

    @pytest.mark.parametrize("allele", [
        "HLA-B*57:01", "HLA-B*58:01", "HLA-B*15:02",
        "HLA-A*31:01", "HLA-A*02:01", "HLA-B*15:11",
        "HLA-B*35:05", "HLA-DPB1*03:01", "HLA-DRB1*07:01",
    ])
    def test_allele_present(self, allele):
        assert allele in HLA_ALLELES

    def test_abacavir_b5701(self):
        info = HLA_ALLELES["HLA-B*57:01"]
        assert "abacavir" in info["drugs"]

    def test_carbamazepine_b1502(self):
        info = HLA_ALLELES["HLA-B*15:02"]
        assert "carbamazepine" in info["drugs"]
        assert any("SJS" in rt for rt in info["reaction_types"])

    def test_each_has_drugs(self):
        for allele, info in HLA_ALLELES.items():
            assert "drugs" in info
            assert len(info["drugs"]) > 0


# ═══════════════════════════════════════════════════════════════════════
# DRUG_CATEGORIES
# ═══════════════════════════════════════════════════════════════════════

class TestDrugCategories:
    def test_has_12_categories(self):
        assert len(DRUG_CATEGORIES) == 12

    @pytest.mark.parametrize("cat", [
        "opioids", "anticoagulants", "antidepressants", "antipsychotics",
        "statins", "chemotherapy", "anticonvulsants", "antivirals",
        "immunosuppressants", "cardiovascular", "proton_pump_inhibitors",
        "anti_gout",
    ])
    def test_category_present(self, cat):
        assert cat in DRUG_CATEGORIES

    def test_opioids_include_codeine(self):
        assert "codeine" in DRUG_CATEGORIES["opioids"]["drugs"]

    def test_anticoagulants_include_warfarin(self):
        assert "warfarin" in DRUG_CATEGORIES["anticoagulants"]["drugs"]

    def test_statins_include_simvastatin(self):
        assert "simvastatin" in DRUG_CATEGORIES["statins"]["drugs"]


# ═══════════════════════════════════════════════════════════════════════
# PGX_SYSTEM_PROMPT
# ═══════════════════════════════════════════════════════════════════════

class TestAgentSystemPrompt:
    def test_not_empty(self):
        assert len(PGX_SYSTEM_PROMPT) > 200

    def test_mentions_critical_alerts(self):
        assert "CRITICAL" in PGX_SYSTEM_PROMPT

    def test_mentions_cpic(self):
        assert "CPIC" in PGX_SYSTEM_PROMPT

    def test_mentions_pmid(self):
        assert "PMID" in PGX_SYSTEM_PROMPT


# ═══════════════════════════════════════════════════════════════════════
# WORKFLOW_COLLECTION_BOOST
# ═══════════════════════════════════════════════════════════════════════

class TestAgentWorkflowBoost:
    def test_has_6_workflows(self):
        assert len(WORKFLOW_COLLECTION_BOOST) == 6

    def test_gene_query_boosts_gene_reference(self):
        boosts = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.GENE_QUERY]
        assert "pgx_gene_reference" in boosts
        assert boosts["pgx_gene_reference"] >= 1.5

    def test_drug_query_boosts_guidelines(self):
        boosts = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.DRUG_QUERY]
        assert "pgx_drug_guidelines" in boosts
        assert boosts["pgx_drug_guidelines"] >= 1.5

    def test_hla_screen_highest_boost(self):
        boosts = WORKFLOW_COLLECTION_BOOST[PGxWorkflowType.HLA_SCREEN]
        assert boosts["pgx_hla_hypersensitivity"] >= 2.0

    def test_boost_values_are_floats(self):
        for wf, boosts in WORKFLOW_COLLECTION_BOOST.items():
            for coll, val in boosts.items():
                assert isinstance(val, float)
                assert val >= 1.0


# ═══════════════════════════════════════════════════════════════════════
# SearchPlan
# ═══════════════════════════════════════════════════════════════════════

class TestSearchPlan:
    def test_creation(self):
        plan = SearchPlan(question="test question")
        assert plan.question == "test question"
        assert plan.genes == []
        assert plan.drugs == []
        assert plan.relevant_workflows == []
        assert plan.search_strategy == "broad"
        assert plan.sub_questions == []

    def test_with_genes_and_drugs(self):
        plan = SearchPlan(
            question="q",
            genes=["CYP2D6"],
            drugs=["codeine"],
            relevant_workflows=[PGxWorkflowType.GENE_QUERY],
        )
        assert plan.genes == ["CYP2D6"]
        assert plan.drugs == ["codeine"]


# ═══════════════════════════════════════════════════════════════════════
# PGxIntelligenceAgent
# ═══════════════════════════════════════════════════════════════════════

class TestPGxIntelligenceAgent:
    @pytest.fixture
    def agent(self):
        mock_rag = MagicMock()
        return PGxIntelligenceAgent(mock_rag)

    def test_init(self, agent):
        assert agent.rag is not None

    # -- search_plan --

    def test_search_plan_identifies_cyp2d6(self, agent):
        plan = agent.search_plan("Is codeine safe for a CYP2D6 ultra-rapid metabolizer?")
        assert "CYP2D6" in plan.genes

    def test_search_plan_identifies_drug(self, agent):
        plan = agent.search_plan("Is codeine safe for a CYP2D6 ultra-rapid metabolizer?")
        assert "codeine" in plan.drugs

    def test_search_plan_drug_safety_subqs(self, agent):
        plan = agent.search_plan("Is codeine safe for CYP2D6 UM?")
        assert len(plan.sub_questions) > 0

    def test_search_plan_dosing_subqs(self, agent):
        plan = agent.search_plan("What dose of warfarin for CYP2C9 *1/*3?")
        assert len(plan.sub_questions) > 0
        assert PGxWorkflowType.DOSING_QUERY in plan.relevant_workflows

    def test_search_plan_hla_workflow(self, agent):
        plan = agent.search_plan("Should I screen for HLA-B*57:01 before abacavir?")
        assert PGxWorkflowType.HLA_SCREEN in plan.relevant_workflows

    def test_search_plan_interaction_subqs(self, agent):
        plan = agent.search_plan("What drug-drug interaction does fluoxetine cause?")
        assert PGxWorkflowType.INTERACTION_QUERY in plan.relevant_workflows

    def test_search_plan_comparative_strategy(self, agent):
        plan = agent.search_plan("Compare clopidogrel vs ticagrelor in CYP2C19 PM")
        assert plan.search_strategy == "comparative"

    def test_search_plan_default_workflow(self, agent):
        plan = agent.search_plan("tell me about genetics")
        assert len(plan.relevant_workflows) >= 1

    def test_search_plan_gene_interpretation(self, agent):
        plan = agent.search_plan("What does CYP2D6 *4/*4 mean for my patient?")
        assert len(plan.sub_questions) > 0

    def test_search_plan_population_strategy(self, agent):
        plan = agent.search_plan("How do CYP2D6 allele frequencies differ across populations?")
        assert plan.search_strategy == "clinical"

    def test_search_plan_targeted_strategy(self, agent):
        plan = agent.search_plan("CYP2D6 codeine interaction")
        assert plan.search_strategy == "targeted"

    def test_search_plan_hla_allele_in_genes(self, agent):
        plan = agent.search_plan("screen for HLA-B*57:01")
        assert any("HLA" in g for g in plan.genes)

    def test_search_plan_topics(self, agent):
        plan = agent.search_plan("codeine dosing for CYP2D6")
        assert len(plan.identified_topics) > 0
