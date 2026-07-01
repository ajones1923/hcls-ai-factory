"""Tests for Biomarker Intelligence Agent Pydantic data models.

Validates all 10 collection models, enum values, embedding text generation,
search result models, patient profile, and analysis result models.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.models import (
    AgentQuery,
    AgingMarker,
    # Patient / Analysis models
    AnalysisResult,
    BiologicalAgeResult,
    BiomarkerReference,
    ClinicalEvidence,
    # Enums
    ClockType,
    CPICLevel,
    # Result models
    CrossCollectionResult,
    DiseaseCategory,
    DiseaseTrajectory,
    DrugInteraction,
    GeneticVariant,
    GenotypeAdjustment,
    MetabolizerPhenotype,
    MonitoringProtocol,
    NutritionGuideline,
    PatientProfile,
    PGxRule,
    RiskLevel,
    SearchHit,
    Zygosity,
)

# =====================================================================
# COLLECTION MODEL CREATION
# =====================================================================


class TestBiomarkerReference:
    """Tests for the BiomarkerReference model."""

    def test_create_with_valid_data(self):
        """BiomarkerReference can be instantiated with all required fields."""
        ref = BiomarkerReference(
            id="ref-albumin",
            name="Albumin",
            text_chunk="Albumin is a liver-synthesized protein reflecting nutritional status.",
            unit="g/dL",
            category="CMP",
            ref_range_min=3.5,
            ref_range_max=5.5,
        )
        assert ref.id == "ref-albumin"
        assert ref.name == "Albumin"
        assert ref.unit == "g/dL"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() produces a non-empty string."""
        ref = BiomarkerReference(
            id="ref-crp",
            name="hs-CRP",
            text_chunk="High-sensitivity C-reactive protein measures systemic inflammation.",
            unit="mg/L",
            category="Inflammation",
            clinical_significance="Elevated hs-CRP indicates cardiovascular risk.",
        )
        text = ref.to_embedding_text()
        assert len(text) > 0
        assert "hs-CRP" in text
        assert "Inflammation" in text
        assert "cardiovascular" in text

    def test_defaults(self):
        """Default values are applied correctly."""
        ref = BiomarkerReference(
            id="ref-1",
            name="Test",
            text_chunk="Test chunk",
        )
        assert ref.unit == ""
        assert ref.category == ""
        assert ref.ref_range_min == 0.0
        assert ref.ref_range_max == 0.0
        assert ref.genetic_modifiers == ""


class TestGeneticVariant:
    """Tests for the GeneticVariant model."""

    def test_create_with_valid_data(self):
        """GeneticVariant can be instantiated with required fields."""
        var = GeneticVariant(
            id="var-mthfr-c677t",
            gene="MTHFR",
            rs_id="rs1801133",
            text_chunk="MTHFR C677T reduces enzyme activity by 30-70%.",
            risk_allele="T",
        )
        assert var.gene == "MTHFR"
        assert var.rs_id == "rs1801133"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes gene and mechanism."""
        var = GeneticVariant(
            id="var-tcf7l2",
            gene="TCF7L2",
            rs_id="rs7903146",
            text_chunk="TCF7L2 rs7903146 T allele increases T2D risk.",
            mechanism="Beta cell dysfunction via Wnt signaling.",
            disease_associations="Type 2 Diabetes",
        )
        text = var.to_embedding_text()
        assert "TCF7L2" in text
        assert "Beta cell" in text
        assert "Type 2 Diabetes" in text


class TestPGxRule:
    """Tests for the PGxRule model."""

    def test_create_with_valid_data(self):
        """PGxRule can be created with required fields."""
        rule = PGxRule(
            id="pgx-001",
            gene="CYP2D6",
            star_alleles="*4/*4",
            drug="Codeine",
            text_chunk="CYP2D6 poor metabolizers cannot convert codeine to morphine.",
        )
        assert rule.gene == "CYP2D6"
        assert rule.drug == "Codeine"
        assert rule.phenotype == MetabolizerPhenotype.NORMAL

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes gene, drug, and recommendation."""
        rule = PGxRule(
            id="pgx-002",
            gene="CYP2C19",
            star_alleles="*2/*2",
            drug="Clopidogrel",
            text_chunk="CYP2C19 poor metabolizers have reduced clopidogrel activation.",
            recommendation="Use prasugrel or ticagrelor instead.",
            phenotype=MetabolizerPhenotype.POOR,
            cpic_level=CPICLevel.LEVEL_1A,
        )
        text = rule.to_embedding_text()
        assert "CYP2C19" in text
        assert "Clopidogrel" in text
        assert "prasugrel" in text


class TestDiseaseTrajectory:
    """Tests for the DiseaseTrajectory model."""

    def test_create_with_valid_data(self):
        """DiseaseTrajectory can be created with required fields."""
        traj = DiseaseTrajectory(
            id="traj-001",
            text_chunk="Pre-diabetes stage with HbA1c 5.7-6.4%.",
            disease=DiseaseCategory.DIABETES,
            stage="pre_diabetic",
        )
        assert traj.disease == DiseaseCategory.DIABETES
        assert traj.stage == "pre_diabetic"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes disease and intervention."""
        traj = DiseaseTrajectory(
            id="traj-002",
            text_chunk="Early cardiovascular risk with elevated LDL.",
            disease=DiseaseCategory.CARDIOVASCULAR,
            stage="elevated_risk",
            intervention_window="Lifestyle modification + statin therapy",
            risk_reduction_pct=30.0,
        )
        text = traj.to_embedding_text()
        assert "cardiovascular" in text
        assert "Lifestyle" in text
        assert "30.0%" in text


class TestClinicalEvidence:
    """Tests for the ClinicalEvidence model."""

    def test_create_with_valid_data(self):
        """ClinicalEvidence can be created with required fields."""
        evi = ClinicalEvidence(
            id="evi-001",
            title="PhenoAge Biomarker Study",
            text_chunk="PhenoAge uses 9 routine blood biomarkers.",
        )
        assert evi.title == "PhenoAge Biomarker Study"
        assert evi.pmid == ""

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes title and finding."""
        evi = ClinicalEvidence(
            id="evi-002",
            title="Levine PhenoAge 2018",
            text_chunk="PhenoAge predicts all-cause mortality.",
            finding="Biological age acceleration correlates with mortality risk.",
            year=2018,
            disease_area="Aging",
        )
        text = evi.to_embedding_text()
        assert "Levine" in text
        assert "mortality" in text


class TestNutritionGuideline:
    """Tests for the NutritionGuideline model."""

    def test_create_with_valid_data(self):
        """NutritionGuideline can be created with required fields."""
        guide = NutritionGuideline(
            id="nut-001",
            nutrient="Folate",
            text_chunk="Folate supplementation for MTHFR carriers.",
        )
        assert guide.nutrient == "Folate"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes nutrient, form, and dose."""
        guide = NutritionGuideline(
            id="nut-002",
            nutrient="Omega-3",
            text_chunk="Omega-3 supplementation for FADS1 CC genotype.",
            genetic_context="FADS1 rs174546 CC",
            recommended_form="EPA/DHA fish oil",
            dose_range="2-4 g/day",
        )
        text = guide.to_embedding_text()
        assert "Omega-3" in text
        assert "FADS1" in text
        assert "EPA/DHA" in text


class TestDrugInteraction:
    """Tests for the DrugInteraction model."""

    def test_create_with_valid_data(self):
        """DrugInteraction can be created with required fields."""
        interact = DrugInteraction(
            id="dint-001",
            drug="Simvastatin",
            gene="SLCO1B1",
            text_chunk="SLCO1B1 variants reduce simvastatin clearance.",
        )
        assert interact.drug == "Simvastatin"
        assert interact.gene == "SLCO1B1"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes drug, gene, and alternative."""
        interact = DrugInteraction(
            id="dint-002",
            drug="Warfarin",
            gene="CYP2C9",
            text_chunk="CYP2C9 variants require warfarin dose adjustment.",
            interaction_type="substrate",
            severity="major",
            alternative="Use alternative anticoagulant or adjust dose",
        )
        text = interact.to_embedding_text()
        assert "Warfarin" in text
        assert "CYP2C9" in text
        assert "alternative" in text.lower()


class TestAgingMarker:
    """Tests for the AgingMarker model."""

    def test_create_with_valid_data(self):
        """AgingMarker can be created with required fields."""
        marker = AgingMarker(
            id="aging-001",
            marker_name="Albumin",
            text_chunk="Albumin declines with aging, contributing to PhenoAge.",
        )
        assert marker.marker_name == "Albumin"
        assert marker.clock_type == ClockType.PHENOAGE

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes marker name and coefficient."""
        marker = AgingMarker(
            id="aging-002",
            marker_name="RDW",
            text_chunk="RDW has the largest PhenoAge coefficient.",
            clock_type=ClockType.PHENOAGE,
            coefficient=0.3306,
            interpretation="Higher RDW indicates accelerated aging.",
        )
        text = marker.to_embedding_text()
        assert "RDW" in text
        assert "PhenoAge" in text


class TestGenotypeAdjustment:
    """Tests for the GenotypeAdjustment model."""

    def test_create_with_valid_data(self):
        """GenotypeAdjustment can be created with required fields."""
        adj = GenotypeAdjustment(
            id="adj-001",
            biomarker="ALT",
            gene="PNPLA3",
            text_chunk="PNPLA3 GG genotype lowers ALT upper limit.",
        )
        assert adj.biomarker == "ALT"
        assert adj.gene == "PNPLA3"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes biomarker, gene, and rationale."""
        adj = GenotypeAdjustment(
            id="adj-002",
            biomarker="Fasting Glucose",
            gene="TCF7L2",
            rs_id="rs7903146",
            text_chunk="TCF7L2 TT genotype warrants tighter glucose threshold.",
            rationale="TT carriers have increased beta-cell dysfunction.",
        )
        text = adj.to_embedding_text()
        assert "TCF7L2" in text
        assert "Fasting Glucose" in text
        assert "beta-cell" in text


class TestMonitoringProtocol:
    """Tests for the MonitoringProtocol model."""

    def test_create_with_valid_data(self):
        """MonitoringProtocol can be created with required fields."""
        proto = MonitoringProtocol(
            id="mon-001",
            condition="Type 2 Diabetes",
            text_chunk="Quarterly HbA1c monitoring for diabetic patients.",
        )
        assert proto.condition == "Type 2 Diabetes"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes condition and panel."""
        proto = MonitoringProtocol(
            id="mon-002",
            condition="Hypothyroidism",
            text_chunk="TSH monitoring for hypothyroid patients on levothyroxine.",
            biomarker_panel="TSH, free T4, free T3",
            frequency="every 6-8 weeks initially",
            trigger_values="TSH >4.0 mIU/L",
        )
        text = proto.to_embedding_text()
        assert "Hypothyroidism" in text
        assert "TSH" in text
        assert "6-8 weeks" in text


# =====================================================================
# PARAMETRIZED: ALL 10 MODELS PRODUCE NON-EMPTY EMBEDDING TEXT
# =====================================================================


@pytest.mark.parametrize(
    "model_cls,kwargs",
    [
        (BiomarkerReference, {"id": "1", "name": "N", "text_chunk": "C"}),
        (GeneticVariant, {"id": "2", "gene": "G", "rs_id": "rs1", "text_chunk": "C"}),
        (PGxRule, {"id": "3", "gene": "G", "star_alleles": "*1/*1", "drug": "D", "text_chunk": "C"}),
        (DiseaseTrajectory, {"id": "4", "text_chunk": "C"}),
        (ClinicalEvidence, {"id": "5", "title": "T", "text_chunk": "C"}),
        (NutritionGuideline, {"id": "6", "nutrient": "N", "text_chunk": "C"}),
        (DrugInteraction, {"id": "7", "drug": "D", "gene": "G", "text_chunk": "C"}),
        (AgingMarker, {"id": "8", "marker_name": "M", "text_chunk": "C"}),
        (GenotypeAdjustment, {"id": "9", "biomarker": "B", "gene": "G", "text_chunk": "C"}),
        (MonitoringProtocol, {"id": "10", "condition": "C", "text_chunk": "C"}),
    ],
    ids=[
        "BiomarkerReference",
        "GeneticVariant",
        "PGxRule",
        "DiseaseTrajectory",
        "ClinicalEvidence",
        "NutritionGuideline",
        "DrugInteraction",
        "AgingMarker",
        "GenotypeAdjustment",
        "MonitoringProtocol",
    ],
)
def test_all_models_embedding_text(model_cls, kwargs):
    """Every collection model's to_embedding_text() returns a non-empty string."""
    instance = model_cls(**kwargs)
    text = instance.to_embedding_text()
    assert isinstance(text, str)
    assert len(text) > 0


# =====================================================================
# ENUM VALUES
# =====================================================================


class TestEnums:
    """Verify enum values match expected domain constants."""

    def test_risk_level_values(self):
        """RiskLevel has 5 expected levels."""
        values = [r.value for r in RiskLevel]
        assert "critical" in values
        assert "high" in values
        assert "moderate" in values
        assert "low" in values
        assert "normal" in values

    def test_clock_type_values(self):
        """ClockType has PhenoAge and GrimAge."""
        assert ClockType.PHENOAGE.value == "PhenoAge"
        assert ClockType.GRIMAGE.value == "GrimAge"

    def test_disease_category_values(self):
        """DiseaseCategory has all 9 expected categories."""
        assert len(DiseaseCategory) == 9
        assert DiseaseCategory.DIABETES.value == "diabetes"
        assert DiseaseCategory.CARDIOVASCULAR.value == "cardiovascular"
        assert DiseaseCategory.LIVER.value == "liver"

    def test_metabolizer_phenotype_values(self):
        """MetabolizerPhenotype has 4 expected phenotypes."""
        assert len(MetabolizerPhenotype) == 4
        assert MetabolizerPhenotype.POOR.value == "poor"
        assert MetabolizerPhenotype.ULTRA_RAPID.value == "ultra_rapid"

    def test_cpic_level_values(self):
        """CPICLevel has 5 evidence levels."""
        assert CPICLevel.LEVEL_1A.value == "1A"
        assert CPICLevel.LEVEL_3.value == "3"

    def test_zygosity_values(self):
        """Zygosity has 3 classifications."""
        assert len(Zygosity) == 3
        assert Zygosity.HETEROZYGOUS.value == "heterozygous"
        assert Zygosity.HOMOZYGOUS_ALT.value == "homozygous_alt"


# =====================================================================
# SEARCH RESULT MODELS
# =====================================================================


class TestSearchHit:
    """Tests for SearchHit creation and fields."""

    def test_create_search_hit(self):
        """SearchHit stores collection, id, score, text, and metadata."""
        hit = SearchHit(
            collection="biomarker_reference",
            id="ref-albumin",
            score=0.85,
            text="Albumin reference range 3.5-5.5 g/dL.",
            metadata={"unit": "g/dL"},
        )
        assert hit.collection == "biomarker_reference"
        assert hit.score == 0.85
        assert hit.metadata["unit"] == "g/dL"

    def test_default_metadata_is_empty_dict(self):
        """SearchHit metadata defaults to an empty dict."""
        hit = SearchHit(collection="test", id="1", score=0.5, text="test")
        assert hit.metadata == {}


class TestCrossCollectionResult:
    """Tests for CrossCollectionResult creation and properties."""

    def test_create_cross_collection_result(self, sample_search_hits):
        """CrossCollectionResult stores query, hits, and metrics."""
        result = CrossCollectionResult(
            query="test query",
            hits=sample_search_hits,
            total_collections_searched=10,
            search_time_ms=50.0,
        )
        assert result.query == "test query"
        assert result.hit_count == 5
        assert result.total_collections_searched == 10

    def test_hits_by_collection(self, sample_search_hits):
        """hits_by_collection() groups results correctly."""
        result = CrossCollectionResult(query="test", hits=sample_search_hits)
        grouped = result.hits_by_collection()
        assert "biomarker_reference" in grouped
        assert "biomarker_genetic_variants" in grouped
        assert len(grouped["biomarker_reference"]) == 1

    def test_empty_result(self):
        """An empty CrossCollectionResult has hit_count == 0."""
        result = CrossCollectionResult(query="empty")
        assert result.hit_count == 0
        assert result.hits_by_collection() == {}


# =====================================================================
# PATIENT / ANALYSIS MODELS
# =====================================================================


class TestPatientProfile:
    """Tests for PatientProfile model."""

    def test_create_patient(self, sample_patient_profile):
        """PatientProfile can be created with all fields."""
        assert sample_patient_profile.patient_id == "HG002"
        assert sample_patient_profile.age == 45
        assert sample_patient_profile.sex == "M"
        assert "albumin" in sample_patient_profile.biomarkers
        assert "TCF7L2_rs7903146" in sample_patient_profile.genotypes
        assert "CYP2D6" in sample_patient_profile.star_alleles

    def test_create_minimal_patient(self):
        """PatientProfile can be created with just required fields."""
        patient = PatientProfile(patient_id="P001", age=30, sex="F")
        assert patient.biomarkers == {}
        assert patient.genotypes == {}
        assert patient.star_alleles == {}


class TestAnalysisResult:
    """Tests for AnalysisResult model."""

    def test_create_analysis_result(self, sample_patient_profile):
        """AnalysisResult can be created with a patient and biological age."""
        result = AnalysisResult(
            patient_profile=sample_patient_profile,
            biological_age=BiologicalAgeResult(
                chronological_age=45,
                biological_age=47.2,
                age_acceleration=2.2,
                phenoage_score=0.003,
            ),
        )
        assert result.biological_age.biological_age == 47.2
        assert result.disease_trajectories == []
        assert result.critical_alerts == []

    def test_timestamp_is_set(self, sample_patient_profile):
        """AnalysisResult has a timestamp set by default."""
        result = AnalysisResult(
            patient_profile=sample_patient_profile,
            biological_age=BiologicalAgeResult(
                chronological_age=45,
                biological_age=45.0,
                age_acceleration=0.0,
                phenoage_score=0.002,
            ),
        )
        assert result.timestamp is not None
        assert "T" in result.timestamp


class TestAgentQuery:
    """Tests for the AgentQuery input model."""

    def test_create_with_required_only(self):
        """AgentQuery requires only the question field."""
        query = AgentQuery(question="What does elevated HbA1c mean?")
        assert query.question == "What does elevated HbA1c mean?"
        assert query.include_genomic is True

    def test_create_with_patient(self, sample_patient_profile):
        """AgentQuery can include a patient profile."""
        query = AgentQuery(
            question="Analyze my biomarkers",
            patient_profile=sample_patient_profile,
        )
        assert query.patient_profile is not None
        assert query.patient_profile.age == 45
