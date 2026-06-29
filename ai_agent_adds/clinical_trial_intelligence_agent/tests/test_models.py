"""Tests for all enums and Pydantic models in src/models.py.

Covers:
  - Enum member counts and values
  - Pydantic model validation (required fields, constraints)
  - Pydantic model serialization (model_dump)
  - Edge cases (boundary values, optional fields)

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    # Enums
    TrialWorkflowType,
    TrialPhase,
    TrialStatus,
    EvidenceLevel,
    CriterionType,
    EndpointType,
    RegulatoryAgency,
    DocumentType,
    SeverityLevel,
    TherapeuticArea,
    DCTComponent,
    # Pydantic models
    TrialQuery,
    TrialSearchResult,
    MatchScore,
    OverallMatch,
    PatientProfile,
    EligibilityAnalysis,
    SiteScore,
    SafetySignal,
    CompetitorProfile,
    ProtocolComplexity,
    WorkflowResult,
    TrialResponse,
    # Dataclass
    SearchPlan,
)


# ===================================================================
# ENUM TESTS
# ===================================================================


class TestTrialWorkflowType:
    """Test TrialWorkflowType enum."""

    def test_member_count(self):
        assert len(TrialWorkflowType) == 19

    def test_all_values(self):
        expected = {
            "protocol_design", "patient_matching", "site_selection",
            "eligibility_optimization", "eligibility_analysis",
            "endpoint_strategy", "adaptive_design",
            "safety_signal", "safety_monitoring",
            "regulatory_docs", "regulatory_strategy",
            "competitive_intel", "competitive_intelligence",
            "biomarker_strategy", "rwe_analysis", "recruitment_optimization",
            "diversity_assessment", "decentralized_planning", "general",
        }
        actual = {m.value for m in TrialWorkflowType}
        assert actual == expected

    def test_string_enum(self):
        assert TrialWorkflowType.PROTOCOL_DESIGN == "protocol_design"
        assert isinstance(TrialWorkflowType.GENERAL, str)


class TestTrialPhase:
    """Test TrialPhase enum."""

    def test_member_count(self):
        assert len(TrialPhase) == 7

    def test_values(self):
        expected = {
            "phase_i", "phase_i_ii", "phase_ii", "phase_ii_iii",
            "phase_iii", "phase_iv", "not_applicable",
        }
        assert {m.value for m in TrialPhase} == expected


class TestTrialStatus:
    """Test TrialStatus enum."""

    def test_member_count(self):
        assert len(TrialStatus) == 7

    def test_values(self):
        expected = {
            "recruiting", "active_not_recruiting", "completed",
            "terminated", "withdrawn", "suspended", "not_yet_recruiting",
        }
        assert {m.value for m in TrialStatus} == expected


class TestEvidenceLevel:
    """Test EvidenceLevel enum."""

    def test_member_count(self):
        assert len(EvidenceLevel) == 6

    def test_hierarchy(self):
        levels = [e.value for e in EvidenceLevel]
        assert "a1" in levels
        assert "e" in levels


class TestCriterionType:
    """Test CriterionType enum."""

    def test_member_count(self):
        assert len(CriterionType) == 2

    def test_values(self):
        assert CriterionType.INCLUSION.value == "inclusion"
        assert CriterionType.EXCLUSION.value == "exclusion"


class TestEndpointType:
    """Test EndpointType enum."""

    def test_member_count(self):
        assert len(EndpointType) == 4

    def test_values(self):
        expected = {"primary", "secondary", "exploratory", "safety"}
        assert {m.value for m in EndpointType} == expected


class TestRegulatoryAgency:
    """Test RegulatoryAgency enum."""

    def test_member_count(self):
        assert len(RegulatoryAgency) == 6

    def test_fda_present(self):
        assert RegulatoryAgency.FDA.value == "fda"


class TestDocumentType:
    """Test DocumentType enum."""

    def test_member_count(self):
        assert len(DocumentType) == 6

    def test_values(self):
        expected = {"ind", "csr", "briefing", "psp", "rmp", "dsur"}
        assert {m.value for m in DocumentType} == expected


class TestSeverityLevel:
    """Test SeverityLevel enum."""

    def test_member_count(self):
        assert len(SeverityLevel) == 5

    def test_values(self):
        expected = {"critical", "high", "moderate", "low", "informational"}
        assert {m.value for m in SeverityLevel} == expected


class TestTherapeuticArea:
    """Test TherapeuticArea enum."""

    def test_member_count(self):
        assert len(TherapeuticArea) == 13

    def test_oncology_present(self):
        assert TherapeuticArea.ONCOLOGY.value == "oncology"

    def test_other_present(self):
        assert TherapeuticArea.OTHER.value == "other"


class TestDCTComponent:
    """Test DCTComponent enum."""

    def test_member_count(self):
        assert len(DCTComponent) == 7

    def test_values(self):
        expected = {
            "econsent", "telemedicine", "home_health", "local_labs",
            "wearables", "epro_ecoa", "direct_to_patient",
        }
        assert {m.value for m in DCTComponent} == expected


# ===================================================================
# PYDANTIC MODEL TESTS
# ===================================================================


class TestTrialQuery:
    """Test TrialQuery model."""

    def test_minimal_valid(self):
        q = TrialQuery(question="What trials are available?")
        assert q.question == "What trials are available?"
        assert q.workflow_type is None
        assert q.top_k == 5

    def test_full_model(self):
        q = TrialQuery(
            question="Match patient to trials",
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            patient_context={"age": 55},
            top_k=10,
            include_guidelines=False,
        )
        assert q.workflow_type == TrialWorkflowType.PATIENT_MATCHING
        assert q.top_k == 10
        assert q.include_guidelines is False

    def test_empty_question_fails(self):
        with pytest.raises(Exception):
            TrialQuery(question="")

    def test_top_k_bounds(self):
        with pytest.raises(Exception):
            TrialQuery(question="test", top_k=0)
        with pytest.raises(Exception):
            TrialQuery(question="test", top_k=51)

    def test_serialization(self):
        q = TrialQuery(question="test query")
        d = q.model_dump()
        assert d["question"] == "test query"
        assert "workflow_type" in d


class TestTrialSearchResult:
    """Test TrialSearchResult model."""

    def test_valid(self):
        r = TrialSearchResult(
            collection="trial_protocols",
            content="KEYNOTE-024 results",
            score=0.92,
        )
        assert r.collection == "trial_protocols"
        assert r.score == 0.92
        assert r.metadata == {}

    def test_negative_score_fails(self):
        with pytest.raises(Exception):
            TrialSearchResult(collection="x", content="y", score=-0.1)


class TestMatchScore:
    """Test MatchScore model."""

    def test_valid(self):
        m = MatchScore(
            criterion_text="Age >= 18",
            criterion_type=CriterionType.INCLUSION,
            met=True,
            confidence=0.95,
        )
        assert m.met is True
        assert m.confidence == 0.95

    def test_confidence_bounds(self):
        with pytest.raises(Exception):
            MatchScore(
                criterion_text="x",
                criterion_type=CriterionType.INCLUSION,
                met=True,
                confidence=1.5,
            )


class TestOverallMatch:
    """Test OverallMatch model."""

    def test_valid(self):
        m = OverallMatch(
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
        assert m.trial_id == "NCT02477436"
        assert m.overall_score == 0.85

    def test_serialization(self):
        m = OverallMatch(
            trial_id="NCT001",
            trial_title="Test Trial",
            phase=TrialPhase.PHASE_II,
            status=TrialStatus.RECRUITING,
            inclusion_met=5,
            inclusion_total=5,
            exclusion_clear=3,
            exclusion_total=3,
            overall_score=1.0,
            confidence=0.8,
        )
        d = m.model_dump()
        assert d["phase"] == "phase_ii"
        assert d["status"] == "recruiting"


class TestPatientProfile:
    """Test PatientProfile model."""

    def test_minimal(self):
        p = PatientProfile()
        assert p.age is None
        assert p.biomarkers == []

    def test_full(self):
        p = PatientProfile(
            age=55,
            sex="male",
            diagnosis="NSCLC",
            biomarkers=["EGFR T790M", "PD-L1 80%"],
            medications=["osimertinib"],
            genomic_variants=["EGFR exon 19 deletion"],
            comorbidities=["hypertension"],
            geographic_location="Boston, MA",
        )
        assert p.age == 55
        assert len(p.biomarkers) == 2

    def test_age_bounds(self):
        with pytest.raises(Exception):
            PatientProfile(age=-1)
        with pytest.raises(Exception):
            PatientProfile(age=121)


class TestEligibilityAnalysis:
    """Test EligibilityAnalysis model."""

    def test_valid(self):
        ea = EligibilityAnalysis(
            criterion="ECOG performance status 0-1",
            population_impact=0.3,
            scientific_justification_score=0.8,
        )
        assert ea.population_impact == 0.3


class TestSiteScore:
    """Test SiteScore model."""

    def test_valid(self):
        s = SiteScore(
            site_id="SITE001",
            facility_name="Mayo Clinic",
            city="Rochester",
            country="USA",
            enrollment_rate=5.2,
            screen_failure_rate=0.15,
            diversity_index=0.7,
            overall_score=0.88,
        )
        assert s.overall_score == 0.88

    def test_rate_bounds(self):
        with pytest.raises(Exception):
            SiteScore(
                site_id="x", facility_name="x", city="x", country="x",
                enrollment_rate=-1, screen_failure_rate=0.1,
                diversity_index=0.5, overall_score=0.5,
            )


class TestSafetySignal:
    """Test SafetySignal model."""

    def test_valid(self):
        s = SafetySignal(
            event_type="Hepatotoxicity",
            severity=SeverityLevel.HIGH,
            frequency=0.05,
            prr=2.3,
            causality_assessment="probable",
        )
        assert s.event_type == "Hepatotoxicity"
        assert s.prr == 2.3


class TestCompetitorProfile:
    """Test CompetitorProfile model."""

    def test_valid(self):
        c = CompetitorProfile(
            trial_id="NCT99999999",
            sponsor="Competitor Corp",
            phase=TrialPhase.PHASE_III,
            indication="NSCLC",
            mechanism="PD-1 inhibitor",
            enrollment_progress=0.75,
            threat_level=SeverityLevel.HIGH,
        )
        assert c.enrollment_progress == 0.75
        assert c.threat_level == SeverityLevel.HIGH


class TestProtocolComplexity:
    """Test ProtocolComplexity model."""

    def test_valid(self):
        pc = ProtocolComplexity(
            procedure_count=25,
            visit_count=15,
            endpoint_count=4,
            eligibility_criteria_count=20,
            complexity_score=0.65,
            percentile_rank=72.0,
        )
        assert pc.complexity_score == 0.65

    def test_score_bounds(self):
        with pytest.raises(Exception):
            ProtocolComplexity(
                procedure_count=0, visit_count=0, endpoint_count=0,
                eligibility_criteria_count=0, complexity_score=1.5,
                percentile_rank=50.0,
            )


class TestWorkflowResult:
    """Test WorkflowResult model."""

    def test_valid(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.PATIENT_MATCHING,
            findings=["Found 3 matching trials"],
            recommendations=["Consider NCT02477436"],
            severity=SeverityLevel.INFORMATIONAL,
            confidence=0.9,
        )
        assert wr.workflow_type == TrialWorkflowType.PATIENT_MATCHING
        assert len(wr.findings) == 1

    def test_defaults(self):
        wr = WorkflowResult(workflow_type=TrialWorkflowType.GENERAL)
        assert wr.findings == []
        assert wr.recommendations == []
        assert wr.severity == SeverityLevel.INFORMATIONAL
        assert wr.confidence == 0.0

    def test_serialization(self):
        wr = WorkflowResult(
            workflow_type=TrialWorkflowType.SAFETY_SIGNAL,
            findings=["Signal detected"],
            confidence=0.75,
        )
        d = wr.model_dump()
        assert d["workflow_type"] == "safety_signal"
        assert d["confidence"] == 0.75


class TestTrialResponse:
    """Test TrialResponse model."""

    def test_valid(self):
        tr = TrialResponse(
            answer="Based on the analysis, 3 trials match.",
            confidence=0.85,
        )
        assert tr.answer.startswith("Based")
        assert tr.citations == []
        assert tr.matches == []

    def test_full(self):
        tr = TrialResponse(
            answer="Analysis complete.",
            citations=[{"collection": "trial_protocols", "score": 0.9}],
            workflow_results=[
                WorkflowResult(
                    workflow_type=TrialWorkflowType.PATIENT_MATCHING,
                    confidence=0.8,
                )
            ],
            confidence=0.8,
        )
        assert len(tr.workflow_results) == 1


class TestSearchPlan:
    """Test SearchPlan dataclass."""

    def test_defaults(self):
        sp = SearchPlan()
        assert sp.question == ""
        assert sp.therapeutic_areas == []
        assert sp.search_strategy == "broad"

    def test_full(self):
        sp = SearchPlan(
            question="Find NSCLC trials with EGFR mutations",
            therapeutic_areas=["oncology"],
            drugs=["osimertinib"],
            biomarkers=["EGFR T790M"],
            relevant_workflows=[TrialWorkflowType.PATIENT_MATCHING],
            search_strategy="targeted",
            sub_questions=["What is EGFR?"],
            identified_topics=["NSCLC", "EGFR"],
        )
        assert sp.search_strategy == "targeted"
        assert len(sp.relevant_workflows) == 1
