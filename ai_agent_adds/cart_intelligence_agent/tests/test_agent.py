"""Tests for CAR-T Intelligence Agent autonomous reasoning module.

Validates SearchPlan creation, search planning logic (target identification,
comparative detection, sub-question generation), and evidence evaluation.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.agent import CARTIntelligenceAgent, SearchPlan
from src.models import CARTStage, CrossCollectionResult, SearchHit


# ═══════════════════════════════════════════════════════════════════════
# FIXTURES
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture
def mock_rag_engine(mock_embedder, mock_llm_client, mock_collection_manager):
    """Return a mock CARTRAGEngine for agent tests."""
    from unittest.mock import MagicMock

    engine = MagicMock()
    engine.retrieve.return_value = CrossCollectionResult(
        query="test",
        hits=[],
        total_collections_searched=10,
        search_time_ms=10.0,
    )
    engine.query.return_value = "Mock answer from LLM."
    return engine


@pytest.fixture
def agent(mock_rag_engine):
    """Return a CARTIntelligenceAgent with a mock RAG engine."""
    return CARTIntelligenceAgent(mock_rag_engine)


# ═══════════════════════════════════════════════════════════════════════
# SEARCH PLAN CREATION
# ═══════════════════════════════════════════════════════════════════════


class TestSearchPlan:
    """Tests for the SearchPlan dataclass."""

    def test_create_default_plan(self):
        """SearchPlan can be created with just a question."""
        plan = SearchPlan(question="What causes CRS?")
        assert plan.question == "What causes CRS?"
        assert plan.identified_topics == []
        assert plan.target_antigens == []
        assert plan.relevant_stages == []
        assert plan.search_strategy == "broad"
        assert plan.sub_questions == []

    def test_create_with_all_fields(self):
        """SearchPlan accepts all fields."""
        plan = SearchPlan(
            question="Compare CD19 vs BCMA",
            identified_topics=["efficacy", "toxicity"],
            target_antigens=["CD19", "BCMA"],
            relevant_stages=[CARTStage.CLINICAL],
            search_strategy="comparative",
            sub_questions=["What are CD19 advantages?"],
        )
        assert plan.search_strategy == "comparative"
        assert len(plan.target_antigens) == 2
        assert len(plan.sub_questions) == 1


# ═══════════════════════════════════════════════════════════════════════
# SEARCH PLAN — TARGET IDENTIFICATION
# ═══════════════════════════════════════════════════════════════════════


class TestSearchPlanTargetIdentification:
    """Tests for search_plan() target antigen identification."""

    def test_identifies_cd19(self, agent):
        """search_plan detects CD19 in the question."""
        plan = agent.search_plan("What is the efficacy of CD19 CAR-T therapy?")
        assert "CD19" in plan.target_antigens

    def test_identifies_bcma(self, agent):
        """search_plan detects BCMA in the question."""
        plan = agent.search_plan("BCMA CAR-T resistance mechanisms")
        assert "BCMA" in plan.target_antigens

    def test_identifies_multiple_targets(self, agent):
        """search_plan detects multiple target antigens."""
        plan = agent.search_plan("Compare CD19 and BCMA CAR-T therapies")
        assert "CD19" in plan.target_antigens
        assert "BCMA" in plan.target_antigens

    @pytest.mark.parametrize(
        "antigen",
        ["CD22", "CD20", "CD30", "CD33", "HER2", "GD2", "GPC3", "PSMA", "ROR1"],
    )
    def test_identifies_various_antigens(self, agent, antigen):
        """search_plan detects a wide range of target antigens."""
        plan = agent.search_plan(f"What trials target {antigen}?")
        assert antigen in plan.target_antigens

    def test_no_antigens_for_generic_query(self, agent):
        """search_plan returns no antigens for a generic manufacturing question."""
        plan = agent.search_plan("What is the optimal expansion protocol?")
        # May detect 0 or just know about generic terms
        # This query has no specific antigen mention
        assert len(plan.target_antigens) == 0


# ═══════════════════════════════════════════════════════════════════════
# SEARCH PLAN — COMPARATIVE DETECTION
# ═══════════════════════════════════════════════════════════════════════


class TestSearchPlanComparative:
    """Tests for search_plan() comparative strategy detection."""

    @pytest.mark.parametrize(
        "question",
        [
            "Compare CD19 and BCMA",
            "CD28 vs 4-1BB costimulatory domains",
            "Kymriah versus Yescarta for DLBCL",
        ],
    )
    def test_detects_comparative_queries(self, agent, question):
        """search_plan sets strategy to 'comparative' for comparison questions."""
        plan = agent.search_plan(question)
        assert plan.search_strategy == "comparative"

    def test_non_comparative_is_broad_or_targeted(self, agent):
        """Non-comparative queries get 'broad' or 'targeted' strategy."""
        plan = agent.search_plan("What is the efficacy of CD19 CAR-T?")
        assert plan.search_strategy in ("broad", "targeted")


# ═══════════════════════════════════════════════════════════════════════
# SEARCH PLAN — STAGE IDENTIFICATION
# ═══════════════════════════════════════════════════════════════════════


class TestSearchPlanStages:
    """Tests for search_plan() development stage identification."""

    def test_clinical_stage_for_trial_query(self, agent):
        """A query about trials identifies the CLINICAL stage."""
        plan = agent.search_plan("What clinical trial results exist for CD19?")
        assert CARTStage.CLINICAL in plan.relevant_stages

    def test_car_design_stage(self, agent):
        """A query about constructs identifies the CAR_DESIGN stage."""
        plan = agent.search_plan("What scFv domains are used in CD19 CARs?")
        assert CARTStage.CAR_DESIGN in plan.relevant_stages

    def test_vector_eng_stage(self, agent):
        """A query about vectors identifies the VECTOR_ENG stage."""
        plan = agent.search_plan("Lentiviral vector production for CAR-T")
        assert CARTStage.VECTOR_ENG in plan.relevant_stages

    def test_testing_stage(self, agent):
        """A query about assays identifies the TESTING stage."""
        plan = agent.search_plan("In vitro cytotoxicity assay results")
        assert CARTStage.TESTING in plan.relevant_stages

    def test_target_id_stage(self, agent):
        """A query about antigen expression identifies the TARGET_ID stage."""
        plan = agent.search_plan("Target antigen expression profiling")
        assert CARTStage.TARGET_ID in plan.relevant_stages


# ═══════════════════════════════════════════════════════════════════════
# SEARCH PLAN — SUB-QUESTION GENERATION
# ═══════════════════════════════════════════════════════════════════════


class TestSearchPlanSubQuestions:
    """Tests for search_plan() sub-question decomposition."""

    def test_why_fail_generates_sub_questions(self, agent):
        """'Why...fail' questions generate resistance/manufacturing/patient sub-queries."""
        plan = agent.search_plan("Why do CD19 CAR-T therapies fail in some patients?")
        assert len(plan.sub_questions) >= 2
        # Check sub-questions cover different angles
        sub_text = " ".join(plan.sub_questions).lower()
        assert "resistance" in sub_text
        assert "manufacturing" in sub_text

    def test_why_fail_includes_antigen_in_sub_questions(self, agent):
        """When an antigen is detected, sub-questions mention it."""
        plan = agent.search_plan("Why does BCMA CAR-T fail?")
        assert len(plan.sub_questions) >= 2
        sub_text = " ".join(plan.sub_questions)
        assert "BCMA" in sub_text

    def test_non_why_fail_has_no_resistance_subs(self, agent):
        """A straightforward efficacy question does not generate failure sub-queries."""
        plan = agent.search_plan("What is the ORR for CD19 CAR-T?")
        # May have sub-questions for comparison, but not failure-specific
        sub_text = " ".join(plan.sub_questions).lower()
        assert "resistance" not in sub_text or len(plan.sub_questions) == 0

    @pytest.mark.parametrize(
        "question,expected_keyword",
        [
            ("How does CD19 CAR-T mechanism work?", "mechanism"),
            ("What biomarkers predict CAR-T response?", "biomarker"),
            ("What are the CAR-T production challenges?", "manufacturing"),
            ("What is the toxicity profile of BCMA CAR-T?", "toxicity"),
            ("What are CAR-T cost and access barriers?", "pricing"),
            ("What is the FDA approval pathway for CAR-T?", "approval"),
        ],
    )
    def test_expanded_patterns_generate_sub_questions(
        self, agent, question, expected_keyword
    ):
        """Expanded question patterns generate relevant sub-questions."""
        plan = agent.search_plan(question)
        assert len(plan.sub_questions) >= 2
        sub_text = " ".join(plan.sub_questions).lower()
        assert expected_keyword in sub_text

    def test_mechanism_pattern_generates_three_sub_questions(self, agent):
        """MECHANISM pattern generates sub-questions about mechanism, evidence, resistance."""
        plan = agent.search_plan("How does CD19 CAR-T mechanism of action work?")
        assert len(plan.sub_questions) == 3
        sub_text = " ".join(plan.sub_questions).lower()
        assert "molecular mechanism" in sub_text
        assert "clinical evidence" in sub_text
        assert "resistance" in sub_text

    def test_safety_pattern_generates_four_sub_questions(self, agent):
        """SAFETY pattern generates sub-questions about incidence, management, risk, biomarkers."""
        plan = agent.search_plan("What is the safety profile of CAR-T toxicity?")
        assert len(plan.sub_questions) == 4
        sub_text = " ".join(plan.sub_questions).lower()
        assert "incidence" in sub_text
        assert "management" in sub_text
        assert "risk factors" in sub_text

    def test_regulatory_pattern(self, agent):
        """REGULATORY / FDA patterns generate approval-related sub-questions."""
        plan = agent.search_plan("What is the FDA regulatory pathway for CAR-T?")
        assert len(plan.sub_questions) == 3
        sub_text = " ".join(plan.sub_questions).lower()
        assert "approval pathway" in sub_text or "regulatory" in sub_text
        assert "precedent" in sub_text

    def test_production_cmc_pattern(self, agent):
        """PRODUCTION / CMC patterns generate manufacturing sub-questions."""
        plan = agent.search_plan("What are the CMC requirements for CAR-T production?")
        assert len(plan.sub_questions) == 3
        sub_text = " ".join(plan.sub_questions).lower()
        assert "process parameters" in sub_text or "manufacturing" in sub_text
        assert "failure modes" in sub_text
        assert "cost" in sub_text

    def test_biomarker_predict_pattern(self, agent):
        """PREDICT / BIOMARKER patterns generate biomarker sub-questions."""
        plan = agent.search_plan("Can we predict CAR-T response with biomarkers?")
        assert len(plan.sub_questions) == 3
        sub_text = " ".join(plan.sub_questions).lower()
        assert "biomarker" in sub_text
        assert "cutoff" in sub_text
        assert "validation" in sub_text

    def test_cost_access_disparity_pattern(self, agent):
        """COST / ACCESS / DISPARITY patterns generate equity sub-questions."""
        plan = agent.search_plan("What are the health disparity issues with CAR-T access?")
        assert len(plan.sub_questions) == 3
        sub_text = " ".join(plan.sub_questions).lower()
        assert "pricing" in sub_text
        assert "access" in sub_text
        assert "equity" in sub_text or "disparity" in sub_text


# ═══════════════════════════════════════════════════════════════════════
# EVIDENCE EVALUATION
# ═══════════════════════════════════════════════════════════════════════


class TestEvaluateEvidence:
    """Tests for evaluate_evidence() quality assessment."""

    def test_sufficient_evidence(self, agent):
        """Evidence from 3+ collections with 10+ hits is 'sufficient'."""
        evidence = CrossCollectionResult(
            query="test",
            hits=[
                SearchHit(collection="Literature", id=str(i), score=0.8, text="hit")
                for i in range(5)
            ]
            + [
                SearchHit(collection="Trial", id=f"NCT{i:08d}", score=0.7, text="hit")
                for i in range(3)
            ]
            + [
                SearchHit(collection="Construct", id=f"c{i}", score=0.6, text="hit")
                for i in range(3)
            ],
        )
        assert agent.evaluate_evidence(evidence) == "sufficient"

    def test_partial_evidence(self, agent):
        """Evidence from 2 collections with 5+ hits is 'partial'."""
        evidence = CrossCollectionResult(
            query="test",
            hits=[
                SearchHit(collection="Literature", id=str(i), score=0.8, text="hit")
                for i in range(3)
            ]
            + [
                SearchHit(collection="Trial", id=f"NCT{i:08d}", score=0.7, text="hit")
                for i in range(3)
            ],
        )
        assert agent.evaluate_evidence(evidence) == "partial"

    def test_insufficient_evidence_no_hits(self, agent):
        """Zero hits means 'insufficient' evidence."""
        evidence = CrossCollectionResult(query="test", hits=[])
        assert agent.evaluate_evidence(evidence) == "insufficient"

    def test_run_passes_stages_to_retrieve(self, agent, mock_rag_engine):
        """run() passes the search plan's relevant_stages to rag.retrieve()."""
        # Ensure llm.generate returns a string so AgentResponse validates
        mock_rag_engine.llm.generate.return_value = "Mock LLM answer"
        mock_rag_engine._build_prompt.return_value = "Mock prompt"
        agent.run("What clinical trial results exist for CD19 CAR-T?")
        # Verify retrieve was called with stages keyword argument
        call_kwargs = mock_rag_engine.retrieve.call_args
        assert "stages" in call_kwargs.kwargs or (
            len(call_kwargs.args) > 1  # positional fallback
        )
        # The stages should include CLINICAL since the question mentions "trial"
        stages_arg = call_kwargs.kwargs.get("stages", None)
        assert stages_arg is not None
        assert CARTStage.CLINICAL in stages_arg

    def test_insufficient_evidence_few_hits(self, agent):
        """Very few hits from one collection is 'insufficient'."""
        evidence = CrossCollectionResult(
            query="test",
            hits=[
                SearchHit(collection="Literature", id="1", score=0.5, text="hit"),
                SearchHit(collection="Literature", id="2", score=0.4, text="hit"),
            ],
        )
        assert agent.evaluate_evidence(evidence) == "insufficient"

    @pytest.mark.parametrize(
        "num_collections,num_hits,expected",
        [
            (0, 0, "insufficient"),
            (1, 3, "insufficient"),
            (2, 6, "partial"),
            (2, 8, "partial"),
            (3, 12, "sufficient"),
            (4, 16, "sufficient"),
        ],
    )
    def test_parametrized_quality_levels(
        self, agent, num_collections, num_hits, expected
    ):
        """Parametrized test of evidence quality classification."""
        collection_names = [
            "Literature", "Trial", "Construct", "Safety",
            "Manufacturing", "Biomarker",
        ]
        hits = []
        # Round-robin distribute hits across the specified number of collections
        for h_idx in range(num_hits):
            coll = collection_names[h_idx % max(num_collections, 1)]
            hits.append(
                SearchHit(
                    collection=coll,
                    id=f"{coll}-{h_idx}",
                    score=0.7,
                    text="evidence",
                )
            )

        evidence = CrossCollectionResult(query="test", hits=hits)
        result = agent.evaluate_evidence(evidence)
        assert result == expected
