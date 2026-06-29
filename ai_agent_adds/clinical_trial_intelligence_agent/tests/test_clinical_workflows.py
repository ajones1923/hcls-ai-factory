"""Tests for clinical workflow result generation.

Tests all 10+1 workflow types with sample inputs, verifying that
WorkflowResult objects can be constructed and serialized correctly.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    TrialWorkflowType,
    SeverityLevel,
    WorkflowResult,
    OverallMatch,
    TrialPhase,
    TrialStatus,
    EligibilityAnalysis,
    SiteScore,
    SafetySignal,
    CompetitorProfile,
)


class TestProtocolDesignWorkflow:
    """Test protocol design workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PROTOCOL_DESIGN,
            findings=[
                "Protocol has 25 procedures across 15 visits",
                "Eligibility criteria count (20) is above median",
            ],
            recommendations=[
                "Consider reducing visit frequency in maintenance phase",
                "Simplify lab panel to essential biomarkers only",
            ],
            severity=SeverityLevel.MODERATE,
            confidence=0.82,
        )
        assert wr.workflow_type == TrialWorkflowType.PROTOCOL_DESIGN
        assert len(wr.findings) == 2

    def test_serialization(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PROTOCOL_DESIGN,
            findings=["Complex protocol identified"],
            confidence=0.7,
        )
        d = wr.model_dump()
        assert d["workflow_type"] == "protocol_design"


class TestPatientMatchingWorkflow:
    """Test patient matching workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            findings=[
                "3 matching Phase III trials found",
                "Best match: NCT02477436 (KEYNOTE-024) with 85% score",
            ],
            recommendations=[
                "Consider NCT02477436 for PD-L1 high NSCLC",
                "NCT03521986 for HER2+ breast cancer",
            ],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.9,
        )
        assert wr.workflow_type == TrialWorkflowType.PATIENT_MATCHING

    def test_with_match_objects(self):
        match = OverallMatch(
            trial_id="NCT02477436",
            trial_title="KEYNOTE-024",
            phase=TrialPhase.PHASE_III,
            status=TrialStatus.COMPLETED,
            inclusion_met=8,
            inclusion_total=10,
            exclusion_clear=5,
            exclusion_total=5,
            overall_score=0.85,
            confidence=0.9,
        )
        assert match.overall_score == 0.85


class TestSiteSelectionWorkflow:
    """Test site selection workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.SITE_SELECTION,
            findings=[
                "Top 5 sites identified by enrollment rate and diversity",
                "Mayo Clinic Rochester leads with 5.2 patients/month",
            ],
            recommendations=[
                "Prioritize high-diversity sites for enrollment equity",
            ],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.85,
        )
        assert wr.workflow_type == TrialWorkflowType.SITE_SELECTION

    def test_site_score_model(self):
        site = SiteScore(
            site_id="SITE001",
            facility_name="Mayo Clinic",
            city="Rochester",
            country="USA",
            enrollment_rate=5.2,
            screen_failure_rate=0.15,
            diversity_index=0.7,
            overall_score=0.88,
        )
        assert site.overall_score == 0.88


class TestEligibilityOptimizationWorkflow:
    """Test eligibility optimization workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.ELIGIBILITY_OPTIMIZATION,
            findings=[
                "3 restrictive criteria identified",
                "ECOG status requirement excludes 30% of target population",
            ],
            recommendations=[
                "Broaden ECOG to 0-2 based on competitor analysis",
            ],
            severity=SeverityLevel.MODERATE,
            confidence=0.78,
        )
        assert wr.workflow_type == TrialWorkflowType.ELIGIBILITY_OPTIMIZATION

    def test_eligibility_analysis_model(self):
        ea = EligibilityAnalysis(
            criterion="ECOG performance status 0-1",
            population_impact=0.3,
            scientific_justification_score=0.6,
            recommendation="Consider broadening to 0-2",
        )
        assert ea.population_impact == 0.3


class TestAdaptiveDesignWorkflow:
    """Test adaptive design workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.ADAPTIVE_DESIGN,
            findings=[
                "Bayesian adaptive design recommended for dose-finding",
                "Sample size re-estimation at interim analysis",
            ],
            recommendations=[
                "Implement response-adaptive randomization",
                "Plan interim analysis at 50% enrollment",
            ],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.75,
        )
        assert wr.workflow_type == TrialWorkflowType.ADAPTIVE_DESIGN


class TestSafetySignalWorkflow:
    """Test safety signal workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.SAFETY_SIGNAL,
            findings=[
                "Potential hepatotoxicity signal detected (PRR=2.3)",
                "Frequency: 5% in treatment vs 1% in control",
            ],
            recommendations=[
                "Implement enhanced liver function monitoring",
                "Consider DILI risk mitigation strategy",
            ],
            severity=SeverityLevel.HIGH,
            confidence=0.85,
        )
        assert wr.severity == SeverityLevel.HIGH

    def test_safety_signal_model(self):
        signal = SafetySignal(
            event_type="Hepatotoxicity",
            severity=SeverityLevel.HIGH,
            frequency=0.05,
            prr=2.3,
            ror=2.5,
            causality_assessment="probable",
        )
        assert signal.prr == 2.3


class TestRegulatoryDocsWorkflow:
    """Test regulatory documentation workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.REGULATORY_DOCS,
            findings=[
                "IND submission package generated",
                "Clinical pharmacology section requires additional PK data",
            ],
            recommendations=[
                "Complete PK bridging study before IND submission",
            ],
            guideline_references=[
                "ICH E6(R2) Good Clinical Practice",
                "FDA Guidance: IND Applications",
            ],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.8,
        )
        assert len(wr.guideline_references) == 2


class TestCompetitiveIntelWorkflow:
    """Test competitive intelligence workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.COMPETITIVE_INTEL,
            findings=[
                "5 competing Phase III trials in same indication",
                "Lead competitor at 75% enrollment",
            ],
            recommendations=[
                "Accelerate enrollment to maintain competitive position",
            ],
            severity=SeverityLevel.MODERATE,
            confidence=0.7,
        )
        assert wr.workflow_type == TrialWorkflowType.COMPETITIVE_INTEL

    def test_competitor_profile(self):
        comp = CompetitorProfile(
            trial_id="NCT99999999",
            sponsor="Competitor Corp",
            phase=TrialPhase.PHASE_III,
            indication="NSCLC",
            mechanism="PD-1 inhibitor",
            enrollment_progress=0.75,
            threat_level=SeverityLevel.HIGH,
        )
        assert comp.enrollment_progress == 0.75


class TestDiversityAssessmentWorkflow:
    """Test diversity assessment workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.DIVERSITY_ASSESSMENT,
            findings=[
                "Current enrollment is 70% Caucasian",
                "Hispanic representation is below target (5% vs 15% target)",
            ],
            recommendations=[
                "Add community health center sites",
                "Implement multilingual recruitment materials",
            ],
            severity=SeverityLevel.MODERATE,
            confidence=0.8,
        )
        assert wr.workflow_type == TrialWorkflowType.DIVERSITY_ASSESSMENT


class TestDecentralizedPlanningWorkflow:
    """Test decentralized trial planning workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.DECENTRALIZED_PLANNING,
            findings=[
                "70% of visits can be conducted remotely",
                "eConsent platform compatible with 95% of site systems",
            ],
            recommendations=[
                "Implement hybrid DCT model with home health visits",
                "Deploy wearable for continuous endpoint monitoring",
            ],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.85,
        )
        assert wr.workflow_type == TrialWorkflowType.DECENTRALIZED_PLANNING


class TestGeneralWorkflow:
    """Test general query workflow."""

    def test_creates_result(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.GENERAL,
            findings=["General information provided"],
            confidence=0.6,
        )
        assert wr.workflow_type == TrialWorkflowType.GENERAL


class TestWorkflowCrossAgentTriggers:
    """Test cross-agent trigger fields in WorkflowResult."""

    def test_with_triggers(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            findings=["Match found"],
            cross_agent_triggers=[
                "oncology:molecular_profiling",
                "pgx:cyp2d6_screening",
                "cardiology:qt_assessment",
            ],
            confidence=0.85,
        )
        assert len(wr.cross_agent_triggers) == 3
        assert "oncology:molecular_profiling" in wr.cross_agent_triggers

    def test_empty_triggers(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.GENERAL,
            findings=["Simple query"],
        )
        assert wr.cross_agent_triggers == []


class TestAllWorkflowTypes:
    """Verify all workflow types can be used in WorkflowResult."""

    @pytest.mark.parametrize("wtype", list(TrialWorkflowType))
    def test_workflow_type_in_result(self, wtype):
        wr = WorkflowResult(
            workflow_type=wtype,
            findings=[f"Test finding for {wtype.value}"],
            confidence=0.5,
        )
        assert wr.workflow_type == wtype
        d = wr.model_dump()
        assert d["workflow_type"] == wtype.value
