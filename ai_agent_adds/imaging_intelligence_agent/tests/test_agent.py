"""Tests for Imaging Intelligence Agent orchestrator.

Validates ImagingIntelligenceAgent.search_plan(), run(),
evaluate_evidence(), generate_report(), and invoke_workflow()
using mock RAG engine, workflow registry, and NIM manager.

Author: Adam Jones
Date: February 2026
"""

import pytest
from unittest.mock import MagicMock, patch

from src.agent import ImagingIntelligenceAgent, SearchPlan
from src.models import (
    AgentQuery,
    AgentResponse,
    BodyRegion,
    CrossCollectionResult,
    FindingSeverity,
    ImagingModality,
    SearchHit,
    WorkflowResult,
    WorkflowStatus,
)
from src.rag_engine import ImagingRAGEngine


# ===================================================================
# FIXTURES
# ===================================================================


@pytest.fixture
def mock_rag_engine(mock_collection_manager, mock_embedder, mock_llm_client):
    """Create a mock-backed RAG engine for agent tests."""
    engine = ImagingRAGEngine(mock_collection_manager, mock_embedder, mock_llm_client)
    return engine


@pytest.fixture
def mock_workflow_registry():
    """Return a dict of mock workflow classes."""
    mock_ct_head = MagicMock()
    mock_ct_head_instance = MagicMock()
    mock_ct_head_instance.run.return_value = WorkflowResult(
        workflow_name="ct_head_hemorrhage",
        status=WorkflowStatus.COMPLETED,
        findings=[{"category": "hemorrhage", "severity": "urgent"}],
        measurements={"volume_ml": 12.5},
        severity=FindingSeverity.URGENT,
        is_mock=True,
    )
    mock_ct_head.return_value = mock_ct_head_instance

    mock_lung = MagicMock()
    mock_lung_instance = MagicMock()
    mock_lung_instance.run.return_value = WorkflowResult(
        workflow_name="ct_chest_lung_nodule",
        status=WorkflowStatus.COMPLETED,
        severity=FindingSeverity.SIGNIFICANT,
        is_mock=True,
    )
    mock_lung.return_value = mock_lung_instance

    return {
        "ct_head_hemorrhage": mock_ct_head,
        "ct_chest_lung_nodule": mock_lung,
    }


@pytest.fixture
def agent(mock_rag_engine, mock_workflow_registry, mock_nim_services):
    """Create an ImagingIntelligenceAgent with mock dependencies."""
    return ImagingIntelligenceAgent(
        rag_engine=mock_rag_engine,
        workflow_registry=mock_workflow_registry,
        nim_manager=mock_nim_services,
    )


# ===================================================================
# search_plan() — MODALITY IDENTIFICATION
# ===================================================================


class TestSearchPlanModalities:
    """Tests for search_plan() modality identification."""

    def test_identifies_ct(self, agent):
        plan = agent.search_plan("CT scan of the head for hemorrhage")
        assert "ct" in plan.modalities

    def test_identifies_mri(self, agent):
        plan = agent.search_plan("brain MRI with contrast")
        assert "mri" in plan.modalities

    def test_identifies_xray(self, agent):
        plan = agent.search_plan("chest x-ray shows pneumothorax")
        assert "xray" in plan.modalities

    def test_identifies_cxr(self, agent):
        plan = agent.search_plan("CXR rapid triage in ED")
        assert "cxr" in plan.modalities

    def test_identifies_ultrasound(self, agent):
        plan = agent.search_plan("renal ultrasound findings")
        assert "ultrasound" in plan.modalities

    def test_identifies_pet(self, agent):
        plan = agent.search_plan("FDG PET for oncology staging")
        assert "pet" in plan.modalities

    def test_identifies_mammography(self, agent):
        plan = agent.search_plan("screening mammography results")
        assert "mammography" in plan.modalities

    def test_no_modality_for_generic_query(self, agent):
        plan = agent.search_plan("What is deep learning?")
        assert len(plan.modalities) == 0


# ===================================================================
# search_plan() — BODY REGION IDENTIFICATION
# ===================================================================


class TestSearchPlanBodyRegions:
    """Tests for search_plan() body region identification."""

    def test_identifies_head(self, agent):
        plan = agent.search_plan("head CT for trauma")
        assert "head" in plan.body_regions

    def test_identifies_brain(self, agent):
        plan = agent.search_plan("brain tumor on MRI")
        assert "head" in plan.body_regions  # 'brain' maps to 'head' keywords

    def test_identifies_chest(self, agent):
        plan = agent.search_plan("chest CT for lung nodule")
        assert "chest" in plan.body_regions

    def test_identifies_abdomen(self, agent):
        plan = agent.search_plan("abdominal CT for liver lesion")
        assert "abdomen" in plan.body_regions

    def test_identifies_spine(self, agent):
        plan = agent.search_plan("lumbar spine MRI")
        assert "spine" in plan.body_regions

    def test_identifies_cardiac(self, agent):
        plan = agent.search_plan("cardiac MRI for ejection fraction")
        assert "cardiac" in plan.body_regions

    def test_identifies_breast(self, agent):
        plan = agent.search_plan("breast MRI findings")
        assert "breast" in plan.body_regions

    def test_no_region_for_generic_query(self, agent):
        plan = agent.search_plan("How does VISTA-3D work?")
        assert len(plan.body_regions) == 0


# ===================================================================
# search_plan() — SEARCH STRATEGY
# ===================================================================


class TestSearchPlanStrategy:
    """Tests for search_plan() strategy determination."""

    def test_comparative_strategy_for_vs(self, agent):
        plan = agent.search_plan("CT vs MRI for hemorrhage detection")
        assert plan.search_strategy == "comparative"

    def test_targeted_strategy_for_specific_query(self, agent):
        plan = agent.search_plan("CT head hemorrhage detection models")
        assert plan.search_strategy == "targeted"

    def test_broad_strategy_for_generic_query(self, agent):
        plan = agent.search_plan("What are the latest AI models?")
        assert plan.search_strategy == "broad"


# ===================================================================
# search_plan() — NIM RECOMMENDATIONS
# ===================================================================


class TestSearchPlanNIMs:
    """Tests for search_plan() NIM recommendations."""

    def test_recommends_vista3d_for_segmentation(self, agent):
        plan = agent.search_plan("organ segmentation from CT volume")
        assert "vista3d" in plan.recommended_nims

    def test_recommends_maisi_for_synthetic(self, agent):
        plan = agent.search_plan("generate synthetic CT training data")
        assert "maisi" in plan.recommended_nims

    def test_recommends_vilam3_for_report(self, agent):
        plan = agent.search_plan("generate a radiology report for this image")
        assert "vilam3" in plan.recommended_nims

    def test_no_nim_for_generic_query(self, agent):
        plan = agent.search_plan("What is Lung-RADS?")
        assert len(plan.recommended_nims) == 0


# ===================================================================
# evaluate_evidence()
# ===================================================================


class TestEvaluateEvidence:
    """Tests for evaluate_evidence() classification."""

    def test_sufficient_with_10_hits(self, agent):
        evidence = CrossCollectionResult(
            query="test",
            hits=[SearchHit(collection="lit", id=str(i), score=0.8, text="t") for i in range(10)],
        )
        assert agent.evaluate_evidence(evidence) == "sufficient"

    def test_sufficient_with_15_hits(self, agent):
        evidence = CrossCollectionResult(
            query="test",
            hits=[SearchHit(collection="lit", id=str(i), score=0.8, text="t") for i in range(15)],
        )
        assert agent.evaluate_evidence(evidence) == "sufficient"

    def test_partial_with_5_hits(self, agent):
        evidence = CrossCollectionResult(
            query="test",
            hits=[SearchHit(collection="lit", id=str(i), score=0.8, text="t") for i in range(5)],
        )
        assert agent.evaluate_evidence(evidence) == "partial"

    def test_partial_with_3_hits(self, agent):
        evidence = CrossCollectionResult(
            query="test",
            hits=[SearchHit(collection="lit", id=str(i), score=0.8, text="t") for i in range(3)],
        )
        assert agent.evaluate_evidence(evidence) == "partial"

    def test_insufficient_with_2_hits(self, agent):
        evidence = CrossCollectionResult(
            query="test",
            hits=[SearchHit(collection="lit", id=str(i), score=0.8, text="t") for i in range(2)],
        )
        assert agent.evaluate_evidence(evidence) == "insufficient"

    def test_insufficient_with_0_hits(self, agent):
        evidence = CrossCollectionResult(query="test", hits=[])
        assert agent.evaluate_evidence(evidence) == "insufficient"


# ===================================================================
# invoke_workflow()
# ===================================================================


class TestInvokeWorkflow:
    """Tests for invoke_workflow()."""

    def test_invokes_known_workflow(self, agent):
        result, cross_modal = agent.invoke_workflow("ct_head_hemorrhage", "/tmp/test.nii.gz")
        assert result is not None
        assert result.workflow_name == "ct_head_hemorrhage"

    def test_returns_none_for_unknown_workflow(self, agent):
        result, cross_modal = agent.invoke_workflow("nonexistent_workflow")
        assert result is None
        assert cross_modal is None

    def test_workflow_receives_input_path(self, agent, mock_workflow_registry):
        agent.invoke_workflow("ct_head_hemorrhage", "/data/ct_head.nii.gz")
        # The mock class is called to create an instance, then .run() is called
        mock_workflow_registry["ct_head_hemorrhage"].assert_called_once()


# ===================================================================
# run() — FULL PIPELINE
# ===================================================================


class TestAgentRun:
    """Tests for run() full agent pipeline."""

    def test_run_returns_agent_response(self, agent):
        query = AgentQuery(question="What are AI models for hemorrhage detection?")
        response = agent.run(query)
        assert isinstance(response, AgentResponse)
        assert response.question == query.question

    def test_run_includes_answer(self, agent):
        query = AgentQuery(question="Brain tumor segmentation models")
        response = agent.run(query)
        assert len(response.answer) > 0

    def test_run_includes_evidence(self, agent):
        query = AgentQuery(question="Lung nodule classification")
        response = agent.run(query)
        assert isinstance(response.evidence, CrossCollectionResult)

    def test_run_has_timestamp(self, agent):
        query = AgentQuery(question="Test question")
        response = agent.run(query)
        assert len(response.timestamp) > 0


# ===================================================================
# generate_report()
# ===================================================================


class TestGenerateReport:
    """Tests for generate_report()."""

    def test_report_is_markdown(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert isinstance(report, str)
        assert "# Imaging Intelligence Report" in report

    def test_report_includes_query(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert sample_agent_response.question in report

    def test_report_includes_answer(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert "3D U-Net" in report  # part of the answer

    def test_report_includes_evidence_summary(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert "Evidence Summary" in report
        assert "Total evidence items" in report

    def test_report_includes_nim_services(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert "NIM Services Used" in report
        assert "vista3d" in report

    def test_report_includes_workflow_results(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert "Workflow Results" in report
        assert "ct_head_hemorrhage" in report

    def test_report_includes_disclaimer(self, agent, sample_agent_response):
        report = agent.generate_report(sample_agent_response)
        assert "clinician review" in report.lower()

    def test_report_without_nim_services(self, agent, sample_evidence):
        response = AgentResponse(
            question="simple question",
            answer="simple answer",
            evidence=sample_evidence,
            nim_services_used=[],
            workflow_results=[],
        )
        report = agent.generate_report(response)
        assert "NIM Services Used" not in report

    def test_report_without_workflow_results(self, agent, sample_evidence):
        response = AgentResponse(
            question="simple question",
            answer="simple answer",
            evidence=sample_evidence,
            workflow_results=[],
        )
        report = agent.generate_report(response)
        assert "Workflow Results" not in report
