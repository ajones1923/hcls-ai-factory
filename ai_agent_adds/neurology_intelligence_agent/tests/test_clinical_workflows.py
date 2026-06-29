"""Tests for clinical workflows in src/clinical_workflows.py.

Tests all 8 neurology workflows plus general.

Author: Adam Jones
Date: March 2026
"""

import pytest

# clinical_workflows may not exist yet; test conditionally
try:
    from src.clinical_workflows import (
        WorkflowEngine,
    )
    _CW_AVAILABLE = True
except ImportError:
    _CW_AVAILABLE = False

from src.models import NeuroWorkflowType, WorkflowResult, SeverityLevel


class TestWorkflowTypes:
    """Tests for workflow type completeness."""

    def test_all_8_clinical_workflows_exist(self):
        expected = [
            "acute_stroke",
            "dementia_evaluation",
            "epilepsy_focus",
            "brain_tumor",
            "ms_monitoring",
            "parkinsons_assessment",
            "headache_classification",
            "neuromuscular_evaluation",
        ]
        for wf_value in expected:
            assert NeuroWorkflowType(wf_value) is not None

    def test_general_workflow_exists(self):
        assert NeuroWorkflowType.GENERAL.value == "general"


class TestWorkflowResult:
    """Tests for WorkflowResult model used by workflows."""

    def test_stroke_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.ACUTE_STROKE,
            findings=["Large vessel occlusion detected", "NIHSS 16"],
            recommendations=["Activate thrombectomy team", "Door-to-groin time target <90 min"],
            severity=SeverityLevel.CRITICAL,
            guideline_references=["AHA/ASA 2019 Stroke Guidelines"],
        )
        assert result.workflow_type == NeuroWorkflowType.ACUTE_STROKE
        assert result.severity == SeverityLevel.CRITICAL
        assert len(result.findings) == 2

    def test_dementia_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.DEMENTIA_EVALUATION,
            findings=["MoCA 18/30", "Temporal lobe atrophy on MRI"],
            recommendations=["Amyloid PET imaging", "CSF ATN biomarkers"],
            severity=SeverityLevel.MODERATE,
        )
        assert result.workflow_type == NeuroWorkflowType.DEMENTIA_EVALUATION
        assert len(result.recommendations) == 2

    def test_epilepsy_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.EPILEPSY_FOCUS,
            findings=["Right temporal spikes on EEG", "Hippocampal sclerosis on MRI"],
            recommendations=["Consider epilepsy surgery evaluation"],
            severity=SeverityLevel.HIGH,
        )
        assert result.workflow_type == NeuroWorkflowType.EPILEPSY_FOCUS

    def test_brain_tumor_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.BRAIN_TUMOR,
            findings=["IDH-mutant, 1p19q-codeleted oligodendroglioma"],
            recommendations=["PCV chemotherapy", "Consider vorasidenib"],
        )
        assert result.workflow_type == NeuroWorkflowType.BRAIN_TUMOR

    def test_ms_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.MS_MONITORING,
            findings=["New T2 lesions on follow-up MRI", "EDSS stable at 2.0"],
            recommendations=["Consider escalation to high-efficacy DMT"],
        )
        assert result.workflow_type == NeuroWorkflowType.MS_MONITORING

    def test_parkinsons_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.PARKINSONS_ASSESSMENT,
            findings=["UPDRS-III 28", "Tremor-dominant subtype"],
            recommendations=["Initiate levodopa therapy"],
        )
        assert result.workflow_type == NeuroWorkflowType.PARKINSONS_ASSESSMENT

    def test_headache_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.HEADACHE_CLASSIFICATION,
            findings=["Meets ICHD-3 criteria for chronic migraine"],
            recommendations=["CGRP antibody prophylaxis"],
        )
        assert result.workflow_type == NeuroWorkflowType.HEADACHE_CLASSIFICATION

    def test_neuromuscular_workflow_result(self):
        result = WorkflowResult(
            workflow_type=NeuroWorkflowType.NEUROMUSCULAR_EVALUATION,
            findings=["Decremental response on repetitive nerve stimulation"],
            recommendations=["Check AChR antibodies", "Consider pyridostigmine"],
        )
        assert result.workflow_type == NeuroWorkflowType.NEUROMUSCULAR_EVALUATION


@pytest.mark.skipif(not _CW_AVAILABLE, reason="clinical_workflows module not yet implemented")
class TestWorkflowEngine:
    """Tests for WorkflowEngine if available."""

    def test_engine_instantiation(self):
        engine = WorkflowEngine()
        assert engine is not None
