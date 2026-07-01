"""Tests for decision support components.

Covers:
  - Protocol complexity scoring
  - Eligibility analysis
  - Safety signal assessment
  - Confidence calibration

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    ProtocolComplexity,
    EligibilityAnalysis,
    SafetySignal,
    SeverityLevel,
    MatchScore,
    CriterionType,
    OverallMatch,
    TrialPhase,
    TrialStatus,
    WorkflowResult,
    TrialWorkflowType,
)


class TestConfidenceCalibration:
    """Test confidence score calibration via WorkflowResult."""

    def test_high_confidence(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            findings=["Strong match found"],
            confidence=0.95,
        )
        assert wr.confidence >= 0.9

    def test_low_confidence(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.GENERAL,
            findings=["Uncertain result"],
            confidence=0.2,
        )
        assert wr.confidence <= 0.3

    def test_confidence_bounds(self):
        with pytest.raises(Exception):
            WorkflowResult(
                workflow_type=TrialWorkflowType.GENERAL,
                confidence=-0.1,
            )
        with pytest.raises(Exception):
            WorkflowResult(
                workflow_type=TrialWorkflowType.GENERAL,
                confidence=1.5,
            )

    def test_zero_confidence(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.GENERAL,
            confidence=0.0,
        )
        assert wr.confidence == 0.0

    def test_full_confidence(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.GENERAL,
            confidence=1.0,
        )
        assert wr.confidence == 1.0


class TestComplexityScorer:
    """Test protocol complexity scoring."""

    def test_simple_protocol(self):
        pc = ProtocolComplexity(
            procedure_count=5,
            visit_count=4,
            endpoint_count=2,
            eligibility_criteria_count=8,
            complexity_score=0.2,
            percentile_rank=15.0,
        )
        assert pc.complexity_score < 0.3
        assert pc.percentile_rank < 25.0

    def test_complex_protocol(self):
        pc = ProtocolComplexity(
            procedure_count=40,
            visit_count=25,
            endpoint_count=8,
            eligibility_criteria_count=35,
            complexity_score=0.85,
            percentile_rank=90.0,
        )
        assert pc.complexity_score > 0.7
        assert pc.percentile_rank > 75.0

    def test_boundary_values(self):
        pc = ProtocolComplexity(
            procedure_count=0,
            visit_count=0,
            endpoint_count=0,
            eligibility_criteria_count=0,
            complexity_score=0.0,
            percentile_rank=0.0,
        )
        assert pc.complexity_score == 0.0

    def test_max_values(self):
        pc = ProtocolComplexity(
            procedure_count=100,
            visit_count=50,
            endpoint_count=20,
            eligibility_criteria_count=50,
            complexity_score=1.0,
            percentile_rank=100.0,
        )
        assert pc.complexity_score == 1.0
        assert pc.percentile_rank == 100.0


class TestEligibilityAnalyzer:
    """Test eligibility criterion analysis."""

    def test_restrictive_criterion(self):
        ea = EligibilityAnalysis(
            criterion="ECOG performance status 0-1 only",
            population_impact=0.35,
            scientific_justification_score=0.6,
            competitor_comparison="Competitors allow ECOG 0-2",
            recommendation="Consider broadening to 0-2",
        )
        assert ea.population_impact > 0.3
        assert ea.scientific_justification_score < 0.7

    def test_well_justified_criterion(self):
        ea = EligibilityAnalysis(
            criterion="Confirmed HER2-positive status by IHC/FISH",
            population_impact=0.5,
            scientific_justification_score=0.95,
            recommendation="Retain - essential for mechanism of action",
        )
        assert ea.scientific_justification_score > 0.9

    def test_impact_bounds(self):
        with pytest.raises(Exception):
            EligibilityAnalysis(
                criterion="x",
                population_impact=-0.1,
                scientific_justification_score=0.5,
            )
        with pytest.raises(Exception):
            EligibilityAnalysis(
                criterion="x",
                population_impact=1.1,
                scientific_justification_score=0.5,
            )

    def test_match_score_met(self):
        ms = MatchScore(
            criterion_text="Age >= 18",
            criterion_type=CriterionType.INCLUSION,
            met=True,
            confidence=0.99,
            evidence="Patient age is 55",
        )
        assert ms.met is True
        assert ms.confidence > 0.9

    def test_match_score_not_met(self):
        ms = MatchScore(
            criterion_text="No prior immunotherapy",
            criterion_type=CriterionType.EXCLUSION,
            met=False,
            confidence=0.85,
            evidence="Patient received pembrolizumab 6 months ago",
        )
        assert ms.met is False


class TestSafetySignalAssessment:
    """Test safety signal detection and assessment."""

    def test_critical_signal(self):
        signal = SafetySignal(
            event_type="Sudden cardiac death",
            severity=SeverityLevel.CRITICAL,
            frequency=0.01,
            prr=5.0,
            ror=4.8,
            causality_assessment="probable",
        )
        assert signal.severity == SeverityLevel.CRITICAL
        assert signal.prr > 2.0

    def test_low_signal(self):
        signal = SafetySignal(
            event_type="Headache",
            severity=SeverityLevel.LOW,
            frequency=0.15,
            prr=0.8,
            causality_assessment="unlikely",
        )
        assert signal.severity == SeverityLevel.LOW
        assert signal.prr < 1.0

    def test_signal_without_prr(self):
        signal = SafetySignal(
            event_type="Nausea",
            severity=SeverityLevel.LOW,
            frequency=0.1,
        )
        assert signal.prr is None
        assert signal.ror is None
        assert signal.causality_assessment == ""

    def test_all_severity_levels(self):
        for severity in SeverityLevel:
            signal = SafetySignal(
                event_type=f"Test event ({severity.value})",
                severity=severity,
                frequency=0.05,
            )
            assert signal.severity == severity


class TestOverallMatchScoring:
    """Test overall patient-trial match scoring."""

    def test_perfect_match(self):
        m = OverallMatch(
            trial_id="NCT001",
            trial_title="Perfect Match Trial",
            phase=TrialPhase.PHASE_III,
            status=TrialStatus.RECRUITING,
            inclusion_met=10,
            inclusion_total=10,
            exclusion_clear=5,
            exclusion_total=5,
            overall_score=1.0,
            confidence=0.95,
        )
        assert m.overall_score == 1.0
        assert m.inclusion_met == m.inclusion_total
        assert m.exclusion_clear == m.exclusion_total

    def test_partial_match(self):
        m = OverallMatch(
            trial_id="NCT002",
            trial_title="Partial Match Trial",
            phase=TrialPhase.PHASE_II,
            status=TrialStatus.RECRUITING,
            inclusion_met=6,
            inclusion_total=10,
            exclusion_clear=4,
            exclusion_total=5,
            overall_score=0.55,
            confidence=0.7,
        )
        assert m.overall_score < 0.6
        assert m.inclusion_met < m.inclusion_total

    def test_no_match(self):
        m = OverallMatch(
            trial_id="NCT003",
            trial_title="No Match Trial",
            phase=TrialPhase.PHASE_I,
            status=TrialStatus.RECRUITING,
            inclusion_met=1,
            inclusion_total=10,
            exclusion_clear=2,
            exclusion_total=5,
            overall_score=0.1,
            confidence=0.9,
        )
        assert m.overall_score < 0.2
