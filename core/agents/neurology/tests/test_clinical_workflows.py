"""Tests for clinical workflows in src/clinical_workflows.py.

Tests all 8 neurology workflows plus general.

Author: Adam Jones
Date: March 2026

CORRECTED (2026-11-03). This file was 125 lines of tests for a 2,577-line module and invoked
none of it. Its eight `TestWorkflowResult` tests -- one per workflow -- each built a
`WorkflowResult` with hard-coded findings and recommendations and then asserted that the enum
round-tripped: `test_stroke_workflow_result` typed "Activate thrombectomy team" into a model
and read it back. No `execute()` was ever called. The only live assertion in the file was
`WorkflowEngine() is not None`, and it sat under
`skipif(not _CW_AVAILABLE, reason="clinical_workflows module not yet implemented")` -- the
module has existed since March 2026, so the guard was inert and the `try/except ImportError`
around the import could only hide a real breakage.

That is why N3, N6, N7, N10, N11, N12, N13, N14, N16 and N19 all shipped inside workflows that
this file appeared to cover.

Each test now dispatches the workflow it names through `WorkflowEngine` and asserts on what
the workflow produced. The behaviour-specific reference cases for the defects live in
`test_stroke_triage_regressions.py` and `test_missing_input_regressions.py`.
"""
import pytest

from src.clinical_workflows import WorkflowEngine
from src.clinical_scales import NIHSSCalculator
from src.models import NeuroWorkflowType, SeverityLevel, WorkflowResult


@pytest.fixture
def engine():
    return WorkflowEngine()


def _complete_nihss(total):
    scores = {item: 0 for item in NIHSSCalculator.ITEMS}
    remaining = total
    for item, max_val in NIHSSCalculator.ITEMS.items():
        take = min(max_val, remaining)
        scores[item] = take
        remaining -= take
    assert remaining == 0
    return scores


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

    def test_every_clinical_workflow_type_is_registered_on_the_engine(self, engine):
        """The enum advertising a workflow is not the same as the engine serving it."""
        registered = set(engine.get_available_workflows())
        for wf in NeuroWorkflowType:
            if wf is NeuroWorkflowType.GENERAL:
                continue
            assert wf in registered, f"{wf.value} is not registered"


class TestWorkflowExecution:
    """Each workflow is dispatched and its own output is asserted on."""

    def test_stroke_workflow_runs_and_reports_the_scale_and_the_type(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.ACUTE_STROKE, {
            "onset_hours": 2.0,
            "nihss_scores": _complete_nihss(16),
            "large_vessel_occlusion": True,
            "prestroke_mrs": 0,
            "ct_affected_regions": ["M1"],
            "blood_pressure_systolic": 150,
            "blood_glucose": 110.0,
            "on_anticoagulation": False,
            "platelet_count": 250,
            "age": 68,
        })
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == NeuroWorkflowType.ACUTE_STROKE
        assert result.severity == SeverityLevel.CRITICAL
        assert any("NIHSS total: 16" in f for f in result.findings)
        assert any("THROMBECTOMY ELIGIBLE" in f for f in result.findings)
        assert any("neurointerventional" in r.lower() for r in result.recommendations)
        assert result.guideline_references

    def test_dementia_workflow_runs_and_stages_the_biomarkers(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.DEMENTIA_EVALUATION, {
            "moca_domain_scores": {
                "visuospatial": 3, "naming": 2, "attention": 4, "language": 2,
                "abstraction": 1, "delayed_recall": 1, "orientation": 5,
            },
            "education_years": 16,
            "atrophy_pattern": "hippocampal",
            "amyloid_status": "positive",
            "tau_status": "positive",
        })
        assert result.workflow_type == NeuroWorkflowType.DEMENTIA_EVALUATION
        assert any("MoCA: 18/30" in f for f in result.findings)
        assert any("A+T+" in f for f in result.findings)
        assert any("alzheimers" in f.lower() for f in result.findings)
        assert result.recommendations

    def test_epilepsy_workflow_runs_and_localises(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.EPILEPSY_FOCUS, {
            "seizure_type": "focal_impaired_awareness",
            "age_at_onset": 22,
            "eeg_findings": {"interictal_focus": "right temporal"},
            "mri_findings": {"mesial_temporal_sclerosis": True,
                             "lesion_location": "right temporal"},
        })
        assert result.workflow_type == NeuroWorkflowType.EPILEPSY_FOCUS
        assert any("temporal" in f.lower() for f in result.findings)
        assert result.recommendations

    def test_brain_tumor_workflow_runs_and_integrates_the_markers(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.BRAIN_TUMOR, {
            "histology": "oligodendroglioma",
            "molecular_markers": ["idh_mutant", "1p19q_codeleted"],
            "size_cm": 3.2,
            "kps": 90,
        })
        assert result.workflow_type == NeuroWorkflowType.BRAIN_TUMOR
        assert any("1p/19q" in f for f in result.findings)
        assert any("PCV" in r for r in result.recommendations)

    def test_ms_workflow_runs_and_reports_disease_activity(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.MS_MONITORING, {
            "phenotype": "rrms",
            "edss_current": 2.0,
            "edss_prior": 2.0,
            "new_t2_lesions": 3,
            "gad_enhancing": 1,
            "relapses_last_year": 1,
            "relapses_last_2yr": 2,
            "current_dmt": "dimethyl_fumarate",
            "dmt_category": "platform",
        })
        assert result.workflow_type == NeuroWorkflowType.MS_MONITORING
        assert any("NEDA-3 STATUS: NOT MET" in f for f in result.findings)
        assert any("high-efficacy" in r.lower() for r in result.recommendations)

    def test_parkinsons_workflow_runs_and_scores_the_motor_exam(self, engine):
        from src.clinical_scales import UPDRSCalculator

        result = engine.run_workflow(NeuroWorkflowType.PARKINSONS_ASSESSMENT, {
            "updrs_scores": {item: 1 for item in UPDRSCalculator.ITEMS},
            "hoehn_yahr_stage": 2.0,
            "symptom_duration_years": 3.0,
            "age": 58,
        })
        assert result.workflow_type == NeuroWorkflowType.PARKINSONS_ASSESSMENT
        assert any("UPDRS" in f for f in result.findings)
        # Levodopa response not recorded, so the workflow asks for the trial rather than
        # assuming one happened.
        assert any("levodopa trial" in r.lower() for r in result.recommendations)

    def test_headache_workflow_runs_and_classifies_per_ichd3(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.HEADACHE_CLASSIFICATION, {
            "headache_features": {
                "unilateral": True, "quality": "throbbing", "intensity": "severe",
                "worsened_by_activity": True, "nausea": True,
                "photophobia": True, "phonophobia": True,
                "duration_hours": 12.0, "frequency_per_month": 18,
                "migraine_days_per_month": 12, "duration_months": 8,
            },
            "num_failed_preventives": 2,
        })
        assert result.workflow_type == NeuroWorkflowType.HEADACHE_CLASSIFICATION
        assert any("chronic migraine" in f.lower() for f in result.findings)
        assert any("CGRP" in r for r in result.recommendations)

    def test_neuromuscular_workflow_runs_and_classifies_the_emg_pattern(self, engine):
        result = engine.run_workflow(NeuroWorkflowType.NEUROMUSCULAR_EVALUATION, {
            "emg_ncs": {"decremental_response": True, "jitter_increased": True},
            "weakness_pattern": "bulbar",
            "progression": "subacute",
        })
        assert result.workflow_type == NeuroWorkflowType.NEUROMUSCULAR_EVALUATION
        assert any("NEUROMUSCULAR JUNCTION" in f for f in result.findings)
        assert any("Myasthenia gravis" in f for f in result.findings)


class TestWorkflowEngine:
    """Tests for WorkflowEngine.

    The import is now unconditional: `try: ... except ImportError: _CW_AVAILABLE = False`
    around a module that has existed since March 2026 could only mask a real breakage, and the
    `skipif` it fed was never applied anyway.
    """

    def test_engine_instantiation(self, engine):
        assert engine is not None

    def test_detect_workflow_routes_a_stroke_question(self, engine):
        detected = engine.detect_workflow("68yo with acute stroke, NIHSS 14, tPA?")
        assert detected == NeuroWorkflowType.ACUTE_STROKE

    def test_an_unknown_workflow_type_raises(self, engine):
        with pytest.raises((KeyError, ValueError)):
            engine.run_workflow("not_a_workflow", {})
