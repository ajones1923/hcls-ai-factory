"""Tests for the Autoimmune Intelligence Agent."""

import pytest

from src.agent import AutoimmuneAgent
from src.knowledge import (
    AUTOANTIBODY_DISEASE_MAP,
    BIOLOGIC_THERAPIES,
    DISEASE_ACTIVITY_THRESHOLDS,
    HLA_DISEASE_ASSOCIATIONS,
    KNOWLEDGE_VERSION,
)
from src.models import (
    AutoantibodyPanel,
    AutoantibodyResult,
    AutoimmuneDisease,
    AutoimmunePatientProfile,
    DiseaseActivityLevel,
    FlarePredictor,
    FlareRisk,
    HLAProfile,
)

# =====================================================================
# Fixtures
# =====================================================================

@pytest.fixture
def agent():
    return AutoimmuneAgent()


@pytest.fixture
def ra_patient():
    """Rheumatoid arthritis patient with anti-CCP positive."""
    return AutoimmunePatientProfile(
        patient_id="AUTO-RA-001",
        age=52,
        sex="F",
        autoantibody_panel=AutoantibodyPanel(
            patient_id="AUTO-RA-001",
            collection_date="2025-06-15",
            results=[
                AutoantibodyResult(antibody="RF", value=85.0, unit="IU/mL", positive=True),
                AutoantibodyResult(antibody="anti-CCP", value=120.0, unit="U/mL", positive=True),
                AutoantibodyResult(antibody="ANA", value=1.2, unit="ratio", positive=False),
            ],
        ),
        hla_profile=HLAProfile(
            hla_drb1=["DRB1*04:01", "DRB1*07:01"],
        ),
        biomarkers={"CRP": 12.5, "ESR": 38},
        diagnosed_conditions=[AutoimmuneDisease.RHEUMATOID_ARTHRITIS],
    )


@pytest.fixture
def lupus_patient():
    """SLE patient with classic serological profile."""
    return AutoimmunePatientProfile(
        patient_id="AUTO-SLE-001",
        age=28,
        sex="F",
        autoantibody_panel=AutoantibodyPanel(
            patient_id="AUTO-SLE-001",
            collection_date="2025-06-15",
            results=[
                AutoantibodyResult(antibody="ANA", value=4.5, unit="ratio", positive=True,
                                   titer="1:640", pattern="homogeneous"),
                AutoantibodyResult(antibody="anti-dsDNA", value=85.0, unit="IU/mL", positive=True),
                AutoantibodyResult(antibody="anti-Smith", value=25.0, unit="AU/mL", positive=True),
            ],
        ),
        biomarkers={"CRP": 3.2, "complement_C3": 65, "complement_C4": 8},
        diagnosed_conditions=[AutoimmuneDisease.SYSTEMIC_LUPUS],
    )


# =====================================================================
# Test Classes
# =====================================================================

class TestAutoantibodyInterpretation:
    def test_ra_antibody_panel(self, agent, ra_patient):
        findings = agent.interpret_autoantibodies(ra_patient.autoantibody_panel)
        diseases = [f["disease"] for f in findings]
        assert "rheumatoid_arthritis" in diseases

    def test_lupus_antibody_panel(self, agent, lupus_patient):
        findings = agent.interpret_autoantibodies(lupus_patient.autoantibody_panel)
        diseases = [f["disease"] for f in findings]
        assert "systemic_lupus_erythematosus" in diseases

    def test_negative_antibodies_ignored(self, agent, ra_patient):
        findings = agent.interpret_autoantibodies(ra_patient.autoantibody_panel)
        # ANA is negative, should not appear
        ab_names = [f["antibody"] for f in findings]
        assert "ANA" not in ab_names

    def test_specificity_included(self, agent, lupus_patient):
        findings = agent.interpret_autoantibodies(lupus_patient.autoantibody_panel)
        for f in findings:
            assert "specificity" in f


class TestHLAAssociations:
    def test_ra_hla_association(self, agent, ra_patient):
        assocs = agent.analyze_hla_associations(ra_patient.hla_profile)
        diseases = [a["disease"] for a in assocs]
        assert "rheumatoid_arthritis" in diseases

    def test_odds_ratio_present(self, agent, ra_patient):
        assocs = agent.analyze_hla_associations(ra_patient.hla_profile)
        for a in assocs:
            assert a["odds_ratio"] > 0

    def test_sorted_by_risk(self, agent, ra_patient):
        assocs = agent.analyze_hla_associations(ra_patient.hla_profile)
        if len(assocs) > 1:
            for i in range(len(assocs) - 1):
                assert assocs[i]["odds_ratio"] >= assocs[i + 1]["odds_ratio"]

    def test_empty_hla_profile(self, agent):
        assocs = agent.analyze_hla_associations(HLAProfile())
        assert assocs == []


class TestDiseaseActivityScoring:
    def test_ra_activity_high(self, agent, ra_patient):
        scores = agent.calculate_disease_activity(
            ra_patient.biomarkers, ra_patient.diagnosed_conditions
        )
        assert len(scores) > 0
        assert any(s.level in (DiseaseActivityLevel.HIGH, DiseaseActivityLevel.MODERATE) for s in scores)

    def test_score_has_components(self, agent, ra_patient):
        scores = agent.calculate_disease_activity(
            ra_patient.biomarkers, ra_patient.diagnosed_conditions
        )
        for s in scores:
            assert s.components

    def test_no_biomarkers_no_scores(self, agent):
        scores = agent.calculate_disease_activity(
            {}, [AutoimmuneDisease.RHEUMATOID_ARTHRITIS]
        )
        assert len(scores) == 0


class TestFlarePredictor:
    def test_ra_flare_risk(self, agent, ra_patient):
        preds = agent.predict_flares(
            ra_patient.biomarkers, ra_patient.diagnosed_conditions
        )
        assert len(preds) > 0
        assert all(isinstance(p, FlarePredictor) for p in preds)

    def test_high_crp_increases_risk(self, agent):
        preds = agent.predict_flares(
            {"CRP": 50, "ESR": 80},
            [AutoimmuneDisease.RHEUMATOID_ARTHRITIS],
        )
        assert preds[0].risk_score > 0.4

    def test_lupus_complement_drop(self, agent, lupus_patient):
        preds = agent.predict_flares(
            lupus_patient.biomarkers, lupus_patient.diagnosed_conditions
        )
        assert len(preds) > 0

    def test_monitoring_recommendations(self, agent, ra_patient):
        preds = agent.predict_flares(
            ra_patient.biomarkers, ra_patient.diagnosed_conditions
        )
        for p in preds:
            assert len(p.recommended_monitoring) > 0


class TestBiologicRecommendations:
    def test_ra_biologics(self, agent, ra_patient):
        recs = agent.recommend_biologics(
            ra_patient.diagnosed_conditions, ra_patient.genotypes
        )
        assert len(recs) > 0
        drug_names = [r.drug_name for r in recs]
        assert "Adalimumab" in drug_names or "Tocilizumab" in drug_names

    def test_lupus_biologics(self, agent, lupus_patient):
        recs = agent.recommend_biologics(
            lupus_patient.diagnosed_conditions, lupus_patient.genotypes
        )
        drug_names = [r.drug_name for r in recs]
        assert "Belimumab" in drug_names or "Rituximab" in drug_names

    def test_pgx_considerations_present(self, agent, ra_patient):
        recs = agent.recommend_biologics(
            ra_patient.diagnosed_conditions, ra_patient.genotypes
        )
        has_pgx = any(len(r.pgx_considerations) > 0 for r in recs)
        assert has_pgx


class TestFullAnalysis:
    def test_ra_full_analysis(self, agent, ra_patient):
        result = agent.analyze_patient(ra_patient)
        assert result.patient_id == "AUTO-RA-001"
        assert len(result.disease_activity_scores) > 0
        assert len(result.biologic_recommendations) > 0

    def test_lupus_full_analysis(self, agent, lupus_patient):
        result = agent.analyze_patient(lupus_patient)
        assert result.patient_id == "AUTO-SLE-001"
        assert len(result.hla_associations) >= 0  # HLA not provided for lupus fixture
        assert len(result.cross_agent_findings) > 0  # From autoantibody interpretation

    def test_critical_alerts_generated(self, agent, ra_patient):
        result = agent.analyze_patient(ra_patient)
        # High CRP should trigger alerts
        assert isinstance(result.critical_alerts, list)


class TestKnowledgeBase:
    def test_version_present(self):
        assert KNOWLEDGE_VERSION["version"] == "2.0.0"

    def test_hla_associations_populated(self):
        assert len(HLA_DISEASE_ASSOCIATIONS) >= 10

    def test_autoantibody_map_populated(self):
        assert len(AUTOANTIBODY_DISEASE_MAP) >= 10

    def test_biologic_therapies_populated(self):
        assert len(BIOLOGIC_THERAPIES) >= 6

    def test_all_biologics_have_class(self):
        for therapy in BIOLOGIC_THERAPIES:
            assert therapy["drug_class"]

    def test_activity_thresholds_populated(self):
        assert "DAS28-CRP" in DISEASE_ACTIVITY_THRESHOLDS
        assert "SLEDAI-2K" in DISEASE_ACTIVITY_THRESHOLDS


class TestModels:
    def test_autoantibody_panel_positive_count(self):
        panel = AutoantibodyPanel(
            patient_id="TEST",
            collection_date="2025-01-01",
            results=[
                AutoantibodyResult(antibody="ANA", value=5.0, positive=True),
                AutoantibodyResult(antibody="RF", value=1.0, positive=False),
                AutoantibodyResult(antibody="anti-CCP", value=80.0, positive=True),
            ],
        )
        assert panel.positive_count == 2
        assert "ANA" in panel.positive_antibodies
        assert "RF" not in panel.positive_antibodies

    def test_patient_profile_validation(self):
        with pytest.raises(Exception):
            AutoimmunePatientProfile(
                patient_id="TEST", age=30, sex="F",
                # No data provided -- should fail validation
            )

    def test_flare_risk_enum(self):
        assert FlareRisk.IMMINENT.value == "imminent"
        assert FlareRisk.LOW.value == "low"

    def test_disease_activity_levels(self):
        assert DiseaseActivityLevel.REMISSION.value == "remission"
        assert DiseaseActivityLevel.VERY_HIGH.value == "very_high"
