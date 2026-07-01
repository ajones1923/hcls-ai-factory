"""Integration tests for clinical trial workflow execution.

Tests that each workflow actually runs end-to-end with sample data,
producing meaningful findings, recommendations, and scores -- not just
model instantiation.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.clinical_workflows import (
    EligibilityOptimizationWorkflow,
    PatientMatchingWorkflow,
    ProtocolDesignWorkflow,
    SafetySignalWorkflow,
    SiteSelectionWorkflow,
    WorkflowEngine,
)
from src.models import TrialWorkflowType, WorkflowResult


# ===================================================================
# TestProtocolDesignExecution
# ===================================================================
class TestProtocolDesignExecution:
    """Run ProtocolDesignWorkflow with realistic inputs."""

    def test_oncology_phase_iii_findings_nonempty(self):
        wf = ProtocolDesignWorkflow()
        result = wf.run({
            "indication": "oncology",
            "phase": "phase_iii",
            "comparator": "placebo",
            "target_population": "adult NSCLC patients",
            "mechanism_of_action": "targeted kinase inhibitor",
        })
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == TrialWorkflowType.PROTOCOL_DESIGN
        assert len(result.findings) > 0

    def test_complexity_score_returned(self):
        wf = ProtocolDesignWorkflow()
        result = wf.run({
            "indication": "cardiology",
            "phase": "phase_ii",
            "comparator": "standard_of_care",
            "design_type": "adaptive",
            "num_arms": 3,
        })
        # Complexity score is embedded in findings text
        complexity_findings = [
            f for f in result.findings if "complexity score" in f.lower()
        ]
        assert len(complexity_findings) >= 1, (
            "Expected at least one finding mentioning complexity score"
        )

    def test_recommendations_nonempty(self):
        wf = ProtocolDesignWorkflow()
        result = wf.run({
            "indication": "neurology",
            "phase": "phase_iii",
            "comparator": "placebo",
        })
        assert len(result.recommendations) > 0

    def test_biomarker_enrichment_bonus(self):
        wf = ProtocolDesignWorkflow()
        result = wf.run({
            "indication": "oncology",
            "phase": "phase_ii",
            "mechanism_of_action": "biomarker-driven targeted therapy",
        })
        enrichment_findings = [
            f for f in result.findings if "enrichment" in f.lower()
        ]
        assert len(enrichment_findings) >= 1


# ===================================================================
# TestPatientMatchingExecution
# ===================================================================
class TestPatientMatchingExecution:
    """Run PatientMatchingWorkflow with sample patient profiles."""

    def test_matching_with_trials(self):
        wf = PatientMatchingWorkflow()
        result = wf.run({
            "patient": {
                "age": 55,
                "sex": "male",
                "diagnosis": "non-small cell lung cancer",
                "biomarkers": ["EGFR+", "PD-L1 high"],
                "medications": ["pembrolizumab"],
                "comorbidities": ["hypertension"],
            },
            "trials": [
                {
                    "trial_id": "NCT-001",
                    "title": "Phase III NSCLC EGFR study",
                    "phase": "phase_iii",
                    "status": "recruiting",
                    "inclusion_criteria": [
                        "Age >= 18 years",
                        "Histologically confirmed non-small cell lung cancer",
                        "EGFR mutation positive",
                    ],
                    "exclusion_criteria": [
                        "Active autoimmune disease",
                        "Brain metastases",
                    ],
                },
                {
                    "trial_id": "NCT-002",
                    "title": "Phase II breast cancer study",
                    "phase": "phase_ii",
                    "status": "recruiting",
                    "inclusion_criteria": [
                        "Female only",
                        "HER2+ breast cancer",
                    ],
                    "exclusion_criteria": [],
                },
            ],
        })
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == TrialWorkflowType.PATIENT_MATCHING
        assert len(result.findings) > 0
        # Should have evaluated 2 trials
        eval_findings = [f for f in result.findings if "Evaluated" in f]
        assert len(eval_findings) >= 1

    def test_missing_patient_fields_graceful(self):
        wf = PatientMatchingWorkflow()
        result = wf.run({
            "patient": {},
            "trials": [
                {
                    "trial_id": "NCT-EMPTY",
                    "title": "Minimal trial",
                    "inclusion_criteria": ["Age >= 18"],
                    "exclusion_criteria": [],
                },
            ],
        })
        assert isinstance(result, WorkflowResult)
        assert len(result.findings) > 0

    def test_no_trials_returns_informational(self):
        wf = PatientMatchingWorkflow()
        result = wf.run({
            "patient": {"age": 30, "diagnosis": "asthma"},
            "trials": [],
        })
        # Preprocess may inject INPUT WARNING; core finding follows
        matching_findings = [
            f for f in result.findings
            if "no trials" in f.lower() or "no matching" in f.lower()
        ]
        assert len(matching_findings) >= 1
        assert result.confidence == 0.0


# ===================================================================
# TestSiteSelectionExecution
# ===================================================================
class TestSiteSelectionExecution:
    """Run SiteSelectionWorkflow with geographic site data."""

    def test_sites_scored_and_ranked(self):
        wf = SiteSelectionWorkflow()
        result = wf.run({
            "sites": [
                {
                    "site_id": "SITE-A",
                    "facility_name": "Memorial Hospital",
                    "city": "New York",
                    "country": "US",
                    "enrollment_rate": 5.0,
                    "screen_failure_rate": 0.20,
                    "investigator_h_index": 45,
                    "therapeutic_experience": 15,
                    "population_access": 2_000_000,
                    "diversity_index": 0.65,
                    "regulatory_readiness": 0.90,
                },
                {
                    "site_id": "SITE-B",
                    "facility_name": "City Clinic",
                    "city": "London",
                    "country": "UK",
                    "enrollment_rate": 2.0,
                    "screen_failure_rate": 0.40,
                    "investigator_h_index": 20,
                    "therapeutic_experience": 5,
                    "population_access": 500_000,
                    "diversity_index": 0.35,
                    "regulatory_readiness": 0.70,
                },
            ],
            "target_enrollment": 200,
            "therapeutic_area": "oncology",
        })
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == TrialWorkflowType.SITE_SELECTION
        assert len(result.findings) > 0
        # Should mention both sites evaluated
        assert any("2 sites" in f for f in result.findings)

    def test_empty_sites_returns_empty_rankings(self):
        wf = SiteSelectionWorkflow()
        result = wf.run({"sites": [], "target_enrollment": 100})
        # Preprocess may inject INPUT WARNING before core finding
        site_findings = [
            f for f in result.findings
            if "no" in f.lower() and ("site" in f.lower() or "candidate" in f.lower())
        ]
        assert len(site_findings) >= 1


# ===================================================================
# TestEligibilityOptimizationExecution
# ===================================================================
class TestEligibilityOptimizationExecution:
    """Run EligibilityOptimizationWorkflow with sample criteria."""

    def test_recommendations_for_restrictive_criteria(self):
        wf = EligibilityOptimizationWorkflow()
        result = wf.run({
            "eligibility_criteria": [
                {"text": "ECOG 0-1 performance status", "type": "inclusion"},
                {"text": "No prior therapy for advanced disease", "type": "inclusion"},
                {"text": "Hemoglobin >= 10 g/dL", "type": "inclusion"},
                {"text": "No brain metastases", "type": "exclusion"},
                {"text": "Creatinine clearance >= 60 mL/min", "type": "inclusion"},
            ],
            "indication": "NSCLC",
        })
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == TrialWorkflowType.ELIGIBILITY_OPTIMIZATION
        assert len(result.findings) > 0
        # Should identify restrictive criteria
        restrictive_findings = [
            f for f in result.findings if "restrictive" in f.lower()
        ]
        assert len(restrictive_findings) >= 1
        assert len(result.recommendations) > 0

    def test_empty_criteria_returns_informational(self):
        wf = EligibilityOptimizationWorkflow()
        result = wf.run({"eligibility_criteria": []})
        assert result.confidence == 0.0


# ===================================================================
# TestSafetySignalExecution
# ===================================================================
class TestSafetySignalExecution:
    """Run SafetySignalWorkflow with sample adverse event data."""

    def test_prr_ror_calculated(self):
        wf = SafetySignalWorkflow()
        result = wf.run({
            "adverse_event": "hepatotoxicity",
            "drug": "DrugX",
            "event_count_drug": 15,
            "total_drug": 200,
            "event_count_comparator": 3,
            "total_comparator": 200,
            "dechallenge_positive": True,
            "rechallenge_positive": False,
            "time_to_onset": "weeks",
            "literature_reports": 5,
        })
        assert isinstance(result, WorkflowResult)
        assert result.workflow_type == TrialWorkflowType.SAFETY_SIGNAL

        # PRR and ROR should be in findings
        prr_findings = [f for f in result.findings if "PRR" in f]
        ror_findings = [f for f in result.findings if "ROR" in f]
        assert len(prr_findings) >= 1, "Expected PRR in findings"
        assert len(ror_findings) >= 1, "Expected ROR in findings"

    def test_signal_detection_with_high_rate(self):
        wf = SafetySignalWorkflow()
        result = wf.run({
            "adverse_event": "rash",
            "drug": "DrugY",
            "event_count_drug": 20,
            "total_drug": 100,
            "event_count_comparator": 2,
            "total_comparator": 100,
            "time_to_onset": "days",
            "dechallenge_positive": True,
            "rechallenge_positive": True,
            "literature_reports": 15,
        })
        signal_findings = [f for f in result.findings if "SIGNAL DETECTED" in f]
        assert len(signal_findings) >= 1
        # Should have recommendations
        assert len(result.recommendations) > 0

    def test_causality_assessment_present(self):
        wf = SafetySignalWorkflow()
        result = wf.run({
            "adverse_event": "nausea",
            "drug": "DrugZ",
            "event_count_drug": 5,
            "total_drug": 100,
            "event_count_comparator": 4,
            "total_comparator": 100,
        })
        causality_findings = [
            f for f in result.findings if "causality" in f.lower()
        ]
        assert len(causality_findings) >= 1


# ===================================================================
# TestWorkflowEngineDispatch
# ===================================================================
class TestWorkflowEngineDispatch:
    """Test WorkflowEngine dispatch and detection."""

    def test_all_workflow_types_dispatch(self):
        engine = WorkflowEngine()
        expected_types = [
            TrialWorkflowType.PROTOCOL_DESIGN,
            TrialWorkflowType.PATIENT_MATCHING,
            TrialWorkflowType.SITE_SELECTION,
            TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
            TrialWorkflowType.ADAPTIVE_DESIGN,
            TrialWorkflowType.SAFETY_SIGNAL,
            TrialWorkflowType.REGULATORY_DOCS,
            TrialWorkflowType.COMPETITIVE_INTEL,
            TrialWorkflowType.DIVERSITY_ASSESSMENT,
            TrialWorkflowType.DECENTRALIZED_PLANNING,
        ]
        for wf_type in expected_types:
            result = engine.run_workflow(wf_type, {})
            assert isinstance(result, WorkflowResult), (
                f"Workflow {wf_type.value} did not return WorkflowResult"
            )
            assert result.workflow_type == wf_type

    def test_unknown_type_raises_error(self):
        engine = WorkflowEngine()
        with pytest.raises(ValueError, match="Unknown workflow type"):
            engine.run_workflow(TrialWorkflowType.GENERAL, {})

    def test_detect_workflow_protocol(self):
        engine = WorkflowEngine()
        detected = engine.detect_workflow("How should I design my trial protocol?")
        assert detected == TrialWorkflowType.PROTOCOL_DESIGN

    def test_detect_workflow_safety(self):
        engine = WorkflowEngine()
        detected = engine.detect_workflow("Analyze this adverse event safety signal")
        assert detected == TrialWorkflowType.SAFETY_SIGNAL

    def test_detect_workflow_general_fallback(self):
        engine = WorkflowEngine()
        detected = engine.detect_workflow("Tell me about something unrelated")
        assert detected == TrialWorkflowType.GENERAL

    def test_list_workflows_complete(self):
        engine = WorkflowEngine()
        workflows = engine.list_workflows()
        assert len(workflows) == 10
        assert "protocol_design" in workflows
        assert "safety_signal" in workflows
        assert "patient_matching" in workflows
