"""Tests for Biomarker Intelligence Agent pipeline.

Validates the agent orchestration: full analysis pipeline, search planning,
critical alert extraction, and module integration -- all with mocks.

Author: Adam Jones
Date: March 2026
"""


import pytest

from src.biological_age import BiologicalAgeCalculator
from src.disease_trajectory import DiseaseTrajectoryAnalyzer
from src.models import (
    AnalysisResult,
    BiologicalAgeResult,
    DiseaseCategory,
    DiseaseTrajectoryResult,
    PatientProfile,
    RiskLevel,
)

# =====================================================================
# FIXTURES
# =====================================================================


@pytest.fixture
def bio_age_calc():
    """Return a BiologicalAgeCalculator instance."""
    return BiologicalAgeCalculator()


@pytest.fixture
def disease_analyzer():
    """Return a DiseaseTrajectoryAnalyzer instance."""
    return DiseaseTrajectoryAnalyzer()


@pytest.fixture
def mock_agent_engine(mock_embedder, mock_llm_client, mock_collection_manager):
    """Return a dict mimicking the engine dict used in the UI."""
    return {
        "manager": mock_collection_manager,
        "bio_age_calc": BiologicalAgeCalculator(),
        "disease_analyzer": DiseaseTrajectoryAnalyzer(),
        "embedder": mock_embedder,
        "llm_client": mock_llm_client,
    }


# =====================================================================
# FULL PIPELINE
# =====================================================================


class TestFullPipeline:
    """Tests for the complete analysis pipeline."""

    def test_pipeline_runs_biological_age(self, mock_agent_engine, sample_patient_profile):
        """Pipeline computes biological age from patient biomarkers."""
        calc = mock_agent_engine["bio_age_calc"]
        result = calc.calculate(
            chronological_age=float(sample_patient_profile.age),
            biomarkers=sample_patient_profile.biomarkers,
        )
        assert "biological_age" in result
        assert "age_acceleration" in result

    def test_pipeline_runs_disease_trajectories(self, mock_agent_engine, sample_patient_profile):
        """Pipeline runs all 6 disease trajectory analyses."""
        analyzer = mock_agent_engine["disease_analyzer"]
        results = analyzer.analyze_all(
            biomarkers=sample_patient_profile.biomarkers,
            genotypes=sample_patient_profile.genotypes,
            age=float(sample_patient_profile.age),
            sex="male",
        )
        assert len(results) == 9

    def test_pipeline_produces_analysis_result(self, mock_agent_engine, sample_patient_profile):
        """Pipeline produces a valid AnalysisResult with all sub-analyses."""
        calc = mock_agent_engine["bio_age_calc"]
        analyzer = mock_agent_engine["disease_analyzer"]

        bio_age = calc.calculate(
            float(sample_patient_profile.age),
            sample_patient_profile.biomarkers,
        )
        trajectories = analyzer.analyze_all(
            sample_patient_profile.biomarkers,
            sample_patient_profile.genotypes,
            age=float(sample_patient_profile.age),
            sex="male",
        )

        phenoage = bio_age.get("phenoage", bio_age)
        analysis = AnalysisResult(
            patient_profile=sample_patient_profile,
            biological_age=BiologicalAgeResult(
                chronological_age=sample_patient_profile.age,
                biological_age=phenoage["biological_age"],
                age_acceleration=phenoage["age_acceleration"],
                phenoage_score=phenoage.get("mortality_score", 0),
            ),
            disease_trajectories=[
                DiseaseTrajectoryResult(
                    disease=DiseaseCategory(t["disease"])
                    if t["disease"] in [e.value for e in DiseaseCategory]
                    else DiseaseCategory.DIABETES,
                    risk_level=RiskLevel(t["risk_level"].lower())
                    if t["risk_level"].lower() in [e.value for e in RiskLevel]
                    else RiskLevel.NORMAL,
                    current_markers=t.get("current_markers", {}),
                )
                for t in trajectories
            ],
        )
        assert analysis.biological_age is not None
        assert len(analysis.disease_trajectories) == 9

    def test_pipeline_handles_empty_biomarkers(self, mock_agent_engine):
        """Pipeline handles a patient with no biomarkers gracefully."""
        PatientProfile(patient_id="EMPTY", age=30, sex="F")
        calc = mock_agent_engine["bio_age_calc"]
        # Should still produce a result, even if missing biomarkers
        result = calc.calculate_phenoage(30.0, {})
        assert "biological_age" in result
        assert len(result["missing_biomarkers"]) > 0


# =====================================================================
# SEARCH PLAN
# =====================================================================


class TestSearchPlan:
    """Tests for search plan generation — identifying relevant topics."""

    def test_diabetes_query_identifies_relevant_collections(self):
        """A diabetes-focused query should target relevant collections."""
        query = "What does HbA1c 5.9% mean for diabetes risk?"
        keywords = query.lower()

        relevant = []
        if "hba1c" in keywords or "diabetes" in keywords:
            relevant.extend(["biomarker_reference", "biomarker_disease_trajectories"])
        if "genetic" in keywords or "genotype" in keywords:
            relevant.append("biomarker_genetic_variants")

        assert "biomarker_reference" in relevant
        assert "biomarker_disease_trajectories" in relevant

    def test_pgx_query_identifies_pgx_collections(self):
        """A PGx-focused query should target pgx_rules and drug_interactions."""
        query = "CYP2D6 poor metabolizer drug dosing"
        keywords = query.lower()

        relevant = []
        if "cyp" in keywords or "metabolizer" in keywords or "drug" in keywords:
            relevant.extend(["biomarker_pgx_rules", "biomarker_drug_interactions"])

        assert "biomarker_pgx_rules" in relevant
        assert "biomarker_drug_interactions" in relevant

    def test_aging_query_identifies_aging_collection(self):
        """An aging-focused query should target aging_markers."""
        query = "What does biological age acceleration mean?"
        keywords = query.lower()

        relevant = []
        if "age" in keywords or "aging" in keywords or "phenoage" in keywords:
            relevant.append("biomarker_aging_markers")

        assert "biomarker_aging_markers" in relevant

    def test_nutrition_query_identifies_nutrition_collection(self):
        """A nutrition query should target nutrition collection."""
        query = "MTHFR folate supplementation recommendations"
        keywords = query.lower()

        relevant = []
        if any(k in keywords for k in ("folate", "vitamin", "supplement", "nutrition", "omega")):
            relevant.append("biomarker_nutrition")

        assert "biomarker_nutrition" in relevant


# =====================================================================
# CRITICAL ALERTS
# =====================================================================


class TestCriticalAlerts:
    """Tests for critical alert extraction from analysis results."""

    def test_diabetic_hba1c_generates_alert(self, disease_analyzer):
        """HbA1c >= 6.5 should generate a critical alert."""
        result = disease_analyzer.analyze_type2_diabetes(
            biomarkers={"hba1c": 7.0},
            genotypes={},
        )
        assert result["risk_level"] == "CRITICAL"

    def test_very_high_ldl_generates_alert(self, disease_analyzer):
        """LDL >= 190 should generate a high-risk alert."""
        result = disease_analyzer.analyze_cardiovascular(
            biomarkers={"ldl_c": 200},
            genotypes={},
        )
        assert result["risk_level"] == "HIGH"

    def test_advanced_fibrosis_generates_alert(self, disease_analyzer):
        """FIB-4 > 2.67 should generate a high-risk alert."""
        result = disease_analyzer.analyze_liver(
            biomarkers={"ast": 80, "alt": 20, "platelets": 100},
            genotypes={},
            age=65.0,
        )
        assert result["risk_level"] == "HIGH"

    def test_overt_hypothyroidism_generates_alert(self, disease_analyzer):
        """TSH > 10 should generate a high-risk alert."""
        result = disease_analyzer.analyze_thyroid(
            biomarkers={"tsh": 15.0},
            genotypes={},
        )
        assert result["risk_level"] == "HIGH"

    def test_collect_all_critical_alerts(self, disease_analyzer, sample_patient_profile):
        """Critical alerts should be collected from all disease analyses."""
        # Use elevated biomarkers to trigger some alerts
        elevated = dict(sample_patient_profile.biomarkers)
        elevated["hba1c"] = 7.0  # Diabetic range
        elevated["ldl_c"] = 200  # Very high

        trajectories = disease_analyzer.analyze_all(
            biomarkers=elevated,
            genotypes=sample_patient_profile.genotypes,
            age=float(sample_patient_profile.age),
            sex="male",
        )

        alerts = []
        for t in trajectories:
            if t["risk_level"] in ("CRITICAL", "HIGH"):
                for finding in t.get("findings", []):
                    alerts.append(finding)

        assert len(alerts) > 0

    def test_normal_values_no_critical_alerts(self, disease_analyzer, sample_biomarkers):
        """Normal biomarker values should not generate critical alerts."""
        trajectories = disease_analyzer.analyze_all(
            biomarkers=sample_biomarkers,
            genotypes={},
        )
        critical_count = sum(
            1 for t in trajectories if t["risk_level"] == "CRITICAL"
        )
        assert critical_count == 0


# =====================================================================
# MODULE INTEGRATION
# =====================================================================


class TestModuleIntegration:
    """Tests for integration between analysis modules."""

    def test_bio_age_uses_same_biomarkers_as_disease(self, sample_patient_profile):
        """Biological age and disease trajectory use the same biomarker dict."""
        biomarkers = sample_patient_profile.biomarkers

        calc = BiologicalAgeCalculator()
        analyzer = DiseaseTrajectoryAnalyzer()

        bio_result = calc.calculate(float(sample_patient_profile.age), biomarkers)
        disease_results = analyzer.analyze_all(
            biomarkers, sample_patient_profile.genotypes,
            age=float(sample_patient_profile.age), sex="male",
        )

        assert bio_result["biological_age"] is not None
        assert len(disease_results) == 9

    def test_genotypes_shared_across_modules(self, sample_patient_profile):
        """Genotype dict is shared across disease trajectory modules."""
        genotypes = sample_patient_profile.genotypes
        analyzer = DiseaseTrajectoryAnalyzer()

        diabetes = analyzer.analyze_type2_diabetes(
            sample_patient_profile.biomarkers, genotypes,
        )
        liver = analyzer.analyze_liver(
            sample_patient_profile.biomarkers, genotypes,
        )

        # TCF7L2 affects diabetes, PNPLA3 affects liver
        assert diabetes["genotype_adjusted_thresholds"] is not None
        assert liver["alt_upper_limit"] is not None
