"""Tests for the AIQ-powered agentic reasoning engine.

Validates the Plan -> Execute -> Reflect -> Refine -> Report cycle
using mocked dependencies (no Milvus, NIM, or AIQ required).

All tests work without the aiq package installed.

Author: Adam Jones
Date: April 2026
"""

import asyncio
from unittest.mock import MagicMock, patch

import pytest

from src.agentic.imaging_agent_aiq import (
    AgenticResult,
    ImagingAgenticEngine,
    ReasoningStep,
    _classify_query,
    _compute_confidence,
    _detect_contradictions,
    _identify_workflow,
)
from src.agentic.tools import IMAGING_TOOLS, ImagingTool
from src.models import (
    CrossCollectionResult,
    CrossModalResult,
    FindingSeverity,
    SearchHit,
    WorkflowResult,
    WorkflowStatus,
)


# ===================================================================
# HELPERS
# ===================================================================


def _run(coro):
    """Run an async coroutine synchronously for tests."""
    return asyncio.get_event_loop().run_until_complete(coro)


def _make_hits(count: int, score: float = 0.8) -> list:
    """Create a list of SearchHit objects."""
    return [
        SearchHit(
            collection="imaging_literature",
            id=f"hit-{i}",
            score=score,
            text=f"Evidence text item {i} about imaging AI.",
        )
        for i in range(count)
    ]


def _make_evidence(count: int, score: float = 0.8) -> CrossCollectionResult:
    """Create a CrossCollectionResult with N hits."""
    return CrossCollectionResult(
        query="test query",
        hits=_make_hits(count, score),
        knowledge_context="Test knowledge context",
        total_collections_searched=11,
        search_time_ms=25.0,
    )


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def mock_rag():
    """Return a mock RAG engine."""
    rag = MagicMock()
    rag.retrieve.return_value = _make_evidence(5)
    rag.retrieve_comparative.return_value = MagicMock(
        evidence_a=_make_evidence(3),
        evidence_b=_make_evidence(3),
        entity_a="ct",
        entity_b="mri",
        shared_evidence=[],
        comparison_context="CT vs MRI comparison",
        total_search_time_ms=50.0,
    )
    rag.query.return_value = (
        "Based on the evidence, AI models achieve high sensitivity "
        "for hemorrhage detection on CT imaging."
    )
    rag.collection_manager = MagicMock()
    rag.collection_manager.get_collection_stats.return_value = {
        "imaging_literature": 100,
        "imaging_trials": 50,
    }
    rag._is_comparative = MagicMock(return_value=False)
    return rag


@pytest.fixture
def mock_workflows():
    """Return a mock workflow registry."""
    mock_ct_head = MagicMock()
    mock_ct_head_instance = MagicMock()
    mock_ct_head_instance.run.return_value = WorkflowResult(
        workflow_name="ct_head_hemorrhage",
        status=WorkflowStatus.COMPLETED,
        findings=[{"category": "hemorrhage", "severity": "urgent"}],
        measurements={"volume_ml": 12.5},
        classification="urgent_hemorrhage",
        severity=FindingSeverity.URGENT,
        is_mock=True,
    )
    mock_ct_head.return_value = mock_ct_head_instance

    mock_lung = MagicMock()
    mock_lung_instance = MagicMock()
    mock_lung_instance.run.return_value = WorkflowResult(
        workflow_name="ct_chest_lung_nodule",
        status=WorkflowStatus.COMPLETED,
        classification="Lung-RADS 4A",
        severity=FindingSeverity.SIGNIFICANT,
        is_mock=True,
    )
    mock_lung.return_value = mock_lung_instance

    return {
        "ct_head_hemorrhage": mock_ct_head,
        "ct_chest_lung_nodule": mock_lung,
    }


@pytest.fixture
def mock_cross_modal():
    """Return a mock CrossModalTrigger."""
    trigger = MagicMock()
    trigger.enabled = True
    trigger._query_genomics.return_value = CrossModalResult(
        trigger_reason="Agentic refinement",
        genomic_context=["[score=0.85] EGFR mutation in lung cancer"],
        genomic_hit_count=3,
        query_count=1,
        enrichment_summary="Retrieved 3 genomic evidence records.",
    )
    trigger.evaluate.return_value = None
    return trigger


@pytest.fixture
def mock_cross_agent():
    """Return a mock cross-agent module."""
    return MagicMock()


@pytest.fixture
def mock_guardrails():
    """Return a mock ClinicalGuardrails."""
    guardrails = MagicMock()

    # Default: input passes, output passes
    input_result = MagicMock()
    input_result.blocked = False
    input_result.modified_text = None
    input_result.flags = []
    guardrails.check_input.return_value = input_result

    output_result = MagicMock()
    output_result.passed = True
    output_result.modified_text = None
    output_result.flags = []
    guardrails.check_output.return_value = output_result

    return guardrails


@pytest.fixture
def engine(mock_rag, mock_workflows, mock_cross_modal, mock_cross_agent):
    """Create an ImagingAgenticEngine with mock dependencies."""
    return ImagingAgenticEngine(
        rag_engine=mock_rag,
        workflow_registry=mock_workflows,
        cross_modal=mock_cross_modal,
        cross_agent=mock_cross_agent,
        max_refinement_rounds=2,
        confidence_threshold=0.7,
        min_evidence_count=3,
    )


@pytest.fixture
def engine_with_guardrails(
    mock_rag, mock_workflows, mock_cross_modal, mock_cross_agent, mock_guardrails
):
    """Create an ImagingAgenticEngine with guardrails."""
    return ImagingAgenticEngine(
        rag_engine=mock_rag,
        workflow_registry=mock_workflows,
        cross_modal=mock_cross_modal,
        cross_agent=mock_cross_agent,
        guardrails=mock_guardrails,
    )


# ===================================================================
# QUERY CLASSIFICATION
# ===================================================================


class TestQueryClassification:
    """Tests for query classification logic."""

    def test_evidence_query(self):
        assert _classify_query("What are the latest deep learning architectures for medical imaging?") == "evidence_query"

    def test_workflow_query_hemorrhage(self):
        assert _classify_query("Analyze this CT head for hemorrhage") == "workflow_query"

    def test_workflow_query_cxr(self):
        assert _classify_query("Triage this chest x-ray") == "workflow_query"

    def test_comparative_query_vs(self):
        assert _classify_query("Compare CT vs MRI for brain hemorrhage") == "comparative_query"

    def test_comparative_query_versus(self):
        assert _classify_query("CT versus MRI for hemorrhage detection") == "comparative_query"

    def test_cross_modal_query(self):
        assert _classify_query("What genomic mutations are associated with this finding?") == "cross_modal_query"

    def test_multi_step_query(self):
        assert _classify_query("Detect hemorrhage and check genomic variants") == "multi_step_query"

    def test_workflow_identification(self):
        assert _identify_workflow("Analyze CT head for hemorrhage") == "ct_head_hemorrhage"
        assert _identify_workflow("Evaluate lung nodule on CT") == "ct_chest_lung_nodule"
        assert _identify_workflow("CXR triage") == "cxr_rapid_findings"
        assert _identify_workflow("What is deep learning?") is None


# ===================================================================
# PLAN
# ===================================================================


class TestPlanEvidenceQuery:
    """test_plan_evidence_query — evidence search query plans evidence search."""

    def test_plan_evidence_query(self, engine):
        step, plan_data = engine._plan(
            "What are the latest deep learning architectures for medical imaging?"
        )
        assert step.step == "plan"
        assert "evidence_query" in step.action
        assert plan_data["search_collections"] is True
        assert plan_data["run_workflow"] is False
        assert plan_data["comparative"] is False

    def test_plan_workflow_query(self, engine):
        """test_plan_workflow_query — workflow query plans workflow execution."""
        step, plan_data = engine._plan(
            "Analyze this CT head for hemorrhage"
        )
        assert step.step == "plan"
        assert "workflow_query" in step.action
        assert plan_data["run_workflow"] is True
        assert plan_data["workflow_name"] == "ct_head_hemorrhage"

    def test_plan_comparative_query(self, engine):
        """test_plan_comparative_query — comparative query plans dual retrieval."""
        step, plan_data = engine._plan(
            "Compare CT vs MRI for brain hemorrhage"
        )
        assert step.step == "plan"
        assert "comparative_query" in step.action
        assert plan_data["comparative"] is True


# ===================================================================
# EXECUTE
# ===================================================================


class TestExecute:
    """test_execute_collects_results — execution step returns evidence."""

    def test_execute_collects_evidence(self, engine, mock_rag):
        step, evidence = _run(
            engine._execute(
                {"search_collections": True, "comparative": False,
                 "run_workflow": False, "cross_modal": False},
                "test query",
                {},
            )
        )
        assert step.step == "execute"
        assert evidence["total_hit_count"] == 5
        mock_rag.retrieve.assert_called_once()

    def test_execute_runs_workflow(self, engine, mock_workflows):
        step, evidence = _run(
            engine._execute(
                {"search_collections": True, "comparative": False,
                 "run_workflow": True, "workflow_name": "ct_head_hemorrhage",
                 "cross_modal": False},
                "Analyze CT head for hemorrhage",
                {},
            )
        )
        assert "ct_head_hemorrhage" in evidence["workflows_executed"]
        mock_workflows["ct_head_hemorrhage"].assert_called_once()

    def test_execute_comparative(self, engine, mock_rag):
        step, evidence = _run(
            engine._execute(
                {"search_collections": True, "comparative": True,
                 "run_workflow": False, "cross_modal": False},
                "CT vs MRI for hemorrhage",
                {},
            )
        )
        mock_rag.retrieve_comparative.assert_called_once()


# ===================================================================
# REFLECT
# ===================================================================


class TestReflect:
    """Tests for the reflection step."""

    def test_reflect_sufficient_evidence(self, engine):
        """test_reflect_sufficient_evidence — 5+ citations -> high confidence."""
        evidence = {
            "total_hit_count": 8,
            "top_scores": [0.9, 0.85, 0.82, 0.78, 0.75],
            "evidence_texts": [f"Evidence {i}" for i in range(8)],
        }
        step, gaps = engine._reflect(evidence)
        assert step.step == "reflect"
        assert step.confidence >= 0.7
        assert "insufficient_evidence" not in gaps

    def test_reflect_insufficient_evidence(self, engine):
        """test_reflect_insufficient_evidence — 1 citation -> low confidence."""
        evidence = {
            "total_hit_count": 1,
            "top_scores": [0.5],
            "evidence_texts": ["Single evidence item"],
        }
        step, gaps = engine._reflect(evidence)
        assert step.confidence < 0.7
        assert "insufficient_evidence" in gaps

    def test_reflect_contradictory_evidence(self, engine):
        """Contradiction detection lowers confidence and flags gap."""
        evidence = {
            "total_hit_count": 6,
            "top_scores": [0.9, 0.85, 0.8],
            "evidence_texts": [
                "CT is superior for hemorrhage detection",
                "CT shows higher sensitivity for acute bleeds",
                "CT is inferior to MRI for chronic hemorrhage",
                "CT is not recommended for posterior fossa",
            ],
        }
        step, gaps = engine._reflect(evidence)
        assert "contradictory_evidence" in gaps
        # Confidence should be lower than without contradiction
        assert step.confidence < 0.85

    def test_reflect_low_relevance_scores(self, engine):
        """Low similarity scores trigger the low_relevance gap."""
        evidence = {
            "total_hit_count": 5,
            "top_scores": [0.3, 0.25, 0.2],
            "evidence_texts": [f"Evidence {i}" for i in range(5)],
        }
        step, gaps = engine._reflect(evidence)
        assert "low_relevance_scores" in gaps


# ===================================================================
# REFINE
# ===================================================================


class TestRefine:
    """Tests for the refinement step."""

    def test_refine_expands_search(self, engine, mock_rag):
        """test_refine_expands_search — low confidence -> additional collections."""
        evidence = {
            "total_hit_count": 1,
            "top_scores": [0.5],
            "evidence_texts": ["Minimal evidence"],
            "cross_modal_triggered": False,
        }
        step, new_evidence = _run(
            engine._refine(
                ["insufficient_evidence"],
                evidence,
                "test query",
                {},
            )
        )
        assert step.step == "refine"
        assert new_evidence["total_hit_count"] > 0
        # Should have called retrieve with expanded top_k
        assert mock_rag.retrieve.call_count >= 1

    def test_refine_triggers_cross_modal(self, engine, mock_cross_modal):
        """test_refine_triggers_cross_modal — genomic enrichment on gap."""
        evidence = {
            "total_hit_count": 5,
            "top_scores": [0.6],
            "evidence_texts": ["Some evidence"],
            "cross_modal_triggered": False,
        }
        step, new_evidence = _run(
            engine._refine(
                ["no_genomic_context"],
                evidence,
                "genomic variant query",
                {},
            )
        )
        assert new_evidence["cross_modal_triggered"] is True
        mock_cross_modal._query_genomics.assert_called_once()


# ===================================================================
# FULL CYCLE
# ===================================================================


class TestFullCycle:
    """test_full_cycle — complete Plan -> Execute -> Reflect -> Refine -> Report."""

    def test_full_reason_cycle(self, engine):
        result = _run(
            engine.reason("What is the evidence for AI lung nodule detection?")
        )
        assert isinstance(result, AgenticResult)
        assert result.query == "What is the evidence for AI lung nodule detection?"
        assert len(result.answer) > 0
        assert len(result.reasoning_chain) >= 3  # plan, execute, reflect, report
        assert result.evidence_count > 0
        assert result.total_duration_ms > 0

    def test_full_cycle_workflow_query(self, engine):
        result = _run(
            engine.reason("Analyze this CT head for hemorrhage")
        )
        assert isinstance(result, AgenticResult)
        assert "ct_head_hemorrhage" in result.workflows_executed

    def test_full_cycle_comparative(self, engine):
        result = _run(
            engine.reason("Compare CT vs MRI for brain hemorrhage")
        )
        assert isinstance(result, AgenticResult)
        assert result.confidence >= 0.0


# ===================================================================
# REASONING CHAIN TRACKING
# ===================================================================


class TestReasoningChain:
    """test_reasoning_chain_tracked — all steps recorded in AgenticResult."""

    def test_chain_has_plan(self, engine):
        result = _run(engine.reason("Test query"))
        steps = [s.step for s in result.reasoning_chain]
        assert "plan" in steps

    def test_chain_has_execute(self, engine):
        result = _run(engine.reason("Test query"))
        steps = [s.step for s in result.reasoning_chain]
        assert "execute" in steps

    def test_chain_has_reflect(self, engine):
        result = _run(engine.reason("Test query"))
        steps = [s.step for s in result.reasoning_chain]
        assert "reflect" in steps

    def test_chain_has_report(self, engine):
        result = _run(engine.reason("Test query"))
        steps = [s.step for s in result.reasoning_chain]
        assert "report" in steps

    def test_chain_step_durations(self, engine):
        result = _run(engine.reason("Test query"))
        for step in result.reasoning_chain:
            assert step.duration_ms >= 0.0

    def test_chain_step_actions_nonempty(self, engine):
        result = _run(engine.reason("Test query"))
        for step in result.reasoning_chain:
            assert len(step.action) > 0


# ===================================================================
# AIQ FALLBACK
# ===================================================================


class TestAIQFallback:
    """test_aiq_fallback — works when aiq package not installed."""

    def test_aiq_not_available(self, engine):
        """Engine should work without AIQ Toolkit."""
        assert engine._aiq_available is False

    def test_reason_without_aiq(self, engine):
        """Full cycle should succeed without AIQ."""
        result = _run(engine.reason("Evidence for lung cancer detection?"))
        assert isinstance(result, AgenticResult)
        assert len(result.answer) > 0

    def test_aiq_detection(self):
        """_check_aiq returns False when aiq is not installed."""
        engine = ImagingAgenticEngine(rag_engine=MagicMock())
        assert engine._aiq_available is False


# ===================================================================
# GUARDRAILS INTEGRATION
# ===================================================================


class TestGuardrailsIntegration:
    """test_guardrails_integration — output checked by guardrails."""

    def test_input_blocked(self, engine_with_guardrails, mock_guardrails):
        """Blocked input returns immediately with zero confidence."""
        blocked_result = MagicMock()
        blocked_result.blocked = True
        blocked_result.flags = ["unsafe_request"]
        blocked_result.modified_text = None
        mock_guardrails.check_input.return_value = blocked_result

        result = _run(engine_with_guardrails.reason("prescribe medication"))
        assert result.confidence == 0.0
        assert "blocked" in result.answer.lower()

    def test_output_modified(self, engine_with_guardrails, mock_guardrails):
        """Output guardrails can modify the final answer."""
        output_result = MagicMock()
        output_result.modified_text = "Modified answer with disclaimer."
        output_result.flags = ["disclaimer_missing"]
        mock_guardrails.check_output.return_value = output_result

        result = _run(
            engine_with_guardrails.reason("Evidence for hemorrhage detection?")
        )
        assert "Modified answer with disclaimer." in result.answer

    def test_guardrails_check_output_called(
        self, engine_with_guardrails, mock_guardrails
    ):
        """check_output is called with the answer and evidence count."""
        _run(
            engine_with_guardrails.reason("Evidence for hemorrhage detection?")
        )
        mock_guardrails.check_output.assert_called_once()


# ===================================================================
# CONFIDENCE SCORING
# ===================================================================


class TestConfidenceScoring:
    """test_confidence_scoring — confidence 0.0-1.0 based on evidence quality."""

    def test_confidence_zero_no_evidence(self):
        c = _compute_confidence(
            evidence_count=0, top_scores=[], min_evidence=3
        )
        assert c == 0.0

    def test_confidence_high_good_evidence(self):
        c = _compute_confidence(
            evidence_count=10,
            top_scores=[0.95, 0.90, 0.88, 0.85, 0.80],
            min_evidence=3,
        )
        assert c >= 0.8

    def test_confidence_reduced_by_contradiction(self):
        c_no = _compute_confidence(
            evidence_count=5,
            top_scores=[0.9, 0.85],
            min_evidence=3,
            contradiction_detected=False,
        )
        c_yes = _compute_confidence(
            evidence_count=5,
            top_scores=[0.9, 0.85],
            min_evidence=3,
            contradiction_detected=True,
        )
        assert c_yes < c_no

    def test_confidence_bounded_0_1(self):
        c = _compute_confidence(
            evidence_count=100,
            top_scores=[1.0, 1.0, 1.0],
            min_evidence=1,
        )
        assert 0.0 <= c <= 1.0

    def test_confidence_partial_evidence(self):
        c = _compute_confidence(
            evidence_count=1,
            top_scores=[0.6],
            min_evidence=5,
        )
        assert 0.0 < c < 0.7


# ===================================================================
# MAX REFINEMENT ROUNDS
# ===================================================================


class TestMaxRefinementRounds:
    """test_max_refinement_rounds — stops after max rounds even if confidence low."""

    def test_stops_at_max_rounds(self, mock_rag, mock_cross_modal):
        """Even with persistently low confidence, refinement is bounded."""
        # Make RAG return very few results so confidence stays low
        mock_rag.retrieve.return_value = _make_evidence(1, score=0.3)

        engine = ImagingAgenticEngine(
            rag_engine=mock_rag,
            cross_modal=mock_cross_modal,
            max_refinement_rounds=2,
            confidence_threshold=0.99,  # Very high threshold
            min_evidence_count=10,  # High bar
        )

        result = _run(engine.reason("Obscure imaging question"))
        assert result.refinement_rounds <= 2

    def test_zero_refinement_when_confident(self, engine):
        """No refinement when initial confidence exceeds threshold."""
        result = _run(
            engine.reason("What is the evidence for AI lung nodule detection?")
        )
        # With 5 hits and 0.8 scores, confidence should be above 0.7
        # So no refinement should occur
        assert result.refinement_rounds == 0


# ===================================================================
# CONTRADICTION DETECTION
# ===================================================================


class TestContradictionDetection:
    """Tests for the lightweight contradiction detector."""

    def test_no_contradiction_uniform(self):
        texts = [
            "CT is preferred for acute hemorrhage",
            "CT recommended as first-line imaging",
        ]
        assert _detect_contradictions(texts) is False

    def test_contradiction_opposing(self):
        texts = [
            "CT is superior for hemorrhage detection",
            "CT is inferior for posterior fossa lesions",
        ]
        assert _detect_contradictions(texts) is True

    def test_no_contradiction_single(self):
        assert _detect_contradictions(["Single text"]) is False

    def test_no_contradiction_empty(self):
        assert _detect_contradictions([]) is False


# ===================================================================
# TOOL REGISTRY
# ===================================================================


class TestToolRegistry:
    """Verify tool definitions are complete and valid."""

    def test_tool_count(self):
        assert len(IMAGING_TOOLS) == 6

    def test_tool_names_unique(self):
        names = [t.name for t in IMAGING_TOOLS]
        assert len(names) == len(set(names))

    def test_tool_has_function(self):
        for tool in IMAGING_TOOLS:
            assert callable(tool.function)

    def test_tool_has_description(self):
        for tool in IMAGING_TOOLS:
            assert len(tool.description) > 10

    def test_expected_tools_present(self):
        names = {t.name for t in IMAGING_TOOLS}
        expected = {
            "search_imaging_evidence",
            "run_clinical_workflow",
            "query_genomic_evidence",
            "query_cross_agent",
            "extract_radiomics",
            "parse_radiology_report",
        }
        assert expected == names


# ===================================================================
# EVIDENCE MERGING
# ===================================================================


class TestEvidenceMerging:
    """Tests for evidence data merging during refinement."""

    def test_merge_adds_hit_counts(self, engine):
        target = {
            "total_hit_count": 5,
            "top_scores": [0.9],
            "evidence_texts": ["text1"],
            "cross_modal_triggered": False,
            "cross_agents_queried": [],
        }
        new = {
            "total_hit_count": 3,
            "top_scores": [0.8],
            "evidence_texts": ["text2"],
            "cross_modal_triggered": True,
            "cross_agents_queried": ["oncology"],
        }
        engine._merge_evidence(target, new)
        assert target["total_hit_count"] == 8
        assert len(target["top_scores"]) == 2
        assert len(target["evidence_texts"]) == 2
        assert target["cross_modal_triggered"] is True
        assert "oncology" in target["cross_agents_queried"]

    def test_merge_preserves_existing_flags(self, engine):
        target = {
            "total_hit_count": 5,
            "top_scores": [],
            "evidence_texts": [],
            "cross_modal_triggered": True,
            "cross_agents_queried": ["cardiology"],
        }
        new = {
            "total_hit_count": 0,
            "top_scores": [],
            "evidence_texts": [],
            "cross_modal_triggered": False,
            "cross_agents_queried": [],
        }
        engine._merge_evidence(target, new)
        assert target["cross_modal_triggered"] is True
        assert "cardiology" in target["cross_agents_queried"]
