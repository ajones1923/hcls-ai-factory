"""Tests for CAR-T Intelligence Agent Pydantic data models.

Validates all 10 collection models, enum values, embedding text generation,
and search result models.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.models import (
    # Enums
    AssayType,
    BiomarkerType,
    CARGeneration,
    CARTStage,
    EvidenceLevel,
    FDAStatus,
    ProcessStep,
    RegulatoryEvent,
    RWEStudyType,
    SafetyEventType,
    SourceType,
    TrialPhase,
    TrialStatus,
    # Collection models
    AssayResult,
    BiomarkerRecord,
    CARConstruct,
    CARTLiterature,
    ClinicalTrial,
    ManufacturingRecord,
    RealWorldRecord,
    RegulatoryRecord,
    SafetyRecord,
    SequenceRecord,
    # Result models
    AgentQuery,
    ComparativeResult,
    CrossCollectionResult,
    SearchHit,
)


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION MODEL CREATION
# ═══════════════════════════════════════════════════════════════════════


class TestCARTLiterature:
    """Tests for the CARTLiterature model."""

    def test_create_with_valid_data(self):
        """CARTLiterature can be instantiated with all required fields."""
        lit = CARTLiterature(
            id="12345678",
            title="CD19 CAR-T therapy in B-ALL",
            text_chunk="This study evaluated tisagenlecleucel in pediatric patients.",
            year=2023,
        )
        assert lit.id == "12345678"
        assert lit.year == 2023
        assert lit.source_type == SourceType.PUBMED

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() produces a non-empty string."""
        lit = CARTLiterature(
            id="99999999",
            title="BCMA CAR-T for multiple myeloma",
            text_chunk="Ciltacabtagene autoleucel achieved 97% ORR.",
            year=2024,
            target_antigen="BCMA",
            disease="Multiple Myeloma",
        )
        text = lit.to_embedding_text()
        assert len(text) > 0
        assert "BCMA" in text
        assert "Multiple Myeloma" in text

    def test_defaults(self):
        """Default values are applied correctly."""
        lit = CARTLiterature(
            id="1",
            title="Test",
            text_chunk="Test chunk",
            year=2025,
        )
        assert lit.source_type == SourceType.PUBMED
        assert lit.cart_stage == CARTStage.CLINICAL
        assert lit.target_antigen == ""
        assert lit.disease == ""


class TestClinicalTrial:
    """Tests for the ClinicalTrial model."""

    def test_create_with_valid_data(self):
        """ClinicalTrial can be instantiated with an NCT ID."""
        trial = ClinicalTrial(
            id="NCT03958656",
            title="A Study of Tisagenlecleucel in Pediatric B-ALL",
            text_summary="Phase 2 study evaluating efficacy and safety.",
        )
        assert trial.id == "NCT03958656"
        assert trial.phase == TrialPhase.NA

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes target and disease info."""
        trial = ClinicalTrial(
            id="NCT12345678",
            title="CD19 CAR-T Trial",
            text_summary="Evaluating CD19-directed CAR-T in DLBCL.",
            target_antigen="CD19",
            disease="DLBCL",
            outcome_summary="ORR 82%",
        )
        text = trial.to_embedding_text()
        assert "CD19" in text
        assert "DLBCL" in text
        assert "ORR 82%" in text


class TestCARConstruct:
    """Tests for the CARConstruct model."""

    def test_create_with_valid_data(self):
        """CARConstruct can be created with required fields."""
        construct = CARConstruct(
            id="construct-001",
            name="Kymriah (tisagenlecleucel)",
            text_summary="Second-generation CAR with 4-1BB costimulation.",
            target_antigen="CD19",
        )
        assert construct.name == "Kymriah (tisagenlecleucel)"
        assert construct.signaling_domain == "CD3-zeta"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes target and generation."""
        construct = CARConstruct(
            id="c-002",
            name="Test CAR",
            text_summary="Test summary",
            target_antigen="BCMA",
            costimulatory_domain="4-1BB",
            known_toxicities="CRS, ICANS",
        )
        text = construct.to_embedding_text()
        assert "BCMA" in text
        assert "4-1BB" in text
        assert "CRS" in text


class TestAssayResult:
    """Tests for the AssayResult model."""

    def test_create_with_valid_data(self):
        """AssayResult can be created with required fields."""
        assay = AssayResult(
            id="assay-001",
            text_summary="Cytotoxicity assay of FMC63 CAR-T against Nalm-6.",
        )
        assert assay.assay_type == AssayType.CYTOTOXICITY
        assert assay.metric_value == 0.0

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes cell line and metrics."""
        assay = AssayResult(
            id="assay-002",
            text_summary="51Cr release assay at 10:1 E:T ratio.",
            cell_line="Nalm-6",
            key_metric="% lysis",
            metric_value=85.3,
            outcome="success",
        )
        text = assay.to_embedding_text()
        assert "Nalm-6" in text
        assert "85.3" in text
        assert "success" in text


class TestManufacturingRecord:
    """Tests for the ManufacturingRecord model."""

    def test_create_with_valid_data(self):
        """ManufacturingRecord can be created with required fields."""
        mfg = ManufacturingRecord(
            id="mfg-001",
            text_summary="Lentiviral transduction at MOI 5.",
        )
        assert mfg.process_step == ProcessStep.TRANSDUCTION
        assert mfg.vector_type == "lentiviral"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes parameter and spec info."""
        mfg = ManufacturingRecord(
            id="mfg-002",
            text_summary="Expansion in G-Rex bioreactor for 10 days.",
            parameter="fold expansion",
            parameter_value="350x",
            target_spec=">100x",
            met_spec="yes",
        )
        text = mfg.to_embedding_text()
        assert "350x" in text
        assert ">100x" in text
        assert "yes" in text


class TestSafetyRecord:
    """Tests for the SafetyRecord model."""

    def test_create_with_valid_data(self):
        """SafetyRecord can be created with required fields."""
        safety = SafetyRecord(
            id="safety-001",
            text_summary="CRS grade 3 in 22% of patients.",
        )
        assert safety.event_type == SafetyEventType.CRS

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes product, event, and management."""
        safety = SafetyRecord(
            id="safety-002",
            text_summary="ICANS management with dexamethasone.",
            product="Yescarta",
            event_type=SafetyEventType.ICANS,
            management_protocol="Dexamethasone 10mg IV q6h",
        )
        text = safety.to_embedding_text()
        assert "Yescarta" in text
        assert "ICANS" in text
        assert "Dexamethasone" in text


class TestBiomarkerRecord:
    """Tests for the BiomarkerRecord model."""

    def test_create_with_valid_data(self):
        """BiomarkerRecord can be created with required fields."""
        bio = BiomarkerRecord(
            id="bio-001",
            text_summary="Ferritin >500 predicts severe CRS.",
        )
        assert bio.biomarker_type == BiomarkerType.PREDICTIVE

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes biomarker name and method."""
        bio = BiomarkerRecord(
            id="bio-002",
            text_summary="CRP as early CRS predictor.",
            biomarker_name="CRP",
            assay_method="Immunoturbidimetry",
            associated_outcome="CRS severity",
        )
        text = bio.to_embedding_text()
        assert "CRP" in text
        assert "Immunoturbidimetry" in text
        assert "CRS severity" in text


class TestRegulatoryRecord:
    """Tests for the RegulatoryRecord model."""

    def test_create_with_valid_data(self):
        """RegulatoryRecord can be created with required fields."""
        reg = RegulatoryRecord(
            id="reg-001",
            text_summary="FDA approval of Kymriah for pediatric B-ALL.",
        )
        assert reg.regulatory_event == RegulatoryEvent.BLA
        assert reg.agency == "FDA"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes product and indication."""
        reg = RegulatoryRecord(
            id="reg-002",
            text_summary="Breakthrough therapy designation for Carvykti.",
            product="Carvykti",
            regulatory_event=RegulatoryEvent.BREAKTHROUGH_THERAPY,
            indication="Multiple Myeloma",
        )
        text = reg.to_embedding_text()
        assert "Carvykti" in text
        assert "breakthrough_therapy" in text
        assert "Multiple Myeloma" in text


class TestSequenceRecord:
    """Tests for the SequenceRecord model."""

    def test_create_with_valid_data(self):
        """SequenceRecord can be created with required fields."""
        seq = SequenceRecord(
            id="seq-001",
            text_summary="FMC63 scFv binding to CD19 at 0.3 nM Kd.",
        )
        assert seq.species_origin == "murine"

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes clone and affinity data."""
        seq = SequenceRecord(
            id="seq-002",
            text_summary="Humanized anti-BCMA nanobody.",
            construct_name="LCAR-B38M",
            scfv_clone="LCAR-B38M VHH",
            binding_affinity_kd="0.18 nM",
            species_origin="humanized",
        )
        text = seq.to_embedding_text()
        assert "LCAR-B38M" in text
        assert "0.18 nM" in text
        assert "humanized" in text


class TestRealWorldRecord:
    """Tests for the RealWorldRecord model."""

    def test_create_with_valid_data(self):
        """RealWorldRecord can be created with required fields."""
        rwe = RealWorldRecord(
            id="rwe-001",
            text_summary="CIBMTR registry outcomes for axi-cel in DLBCL.",
        )
        assert rwe.study_type == RWEStudyType.RETROSPECTIVE

    def test_to_embedding_text_non_empty(self):
        """to_embedding_text() includes product and endpoint."""
        rwe = RealWorldRecord(
            id="rwe-002",
            text_summary="Community vs academic center outcomes.",
            product="Yescarta",
            data_source="CIBMTR",
            primary_endpoint="12-month PFS",
            outcome_value="42%",
            special_population="elderly (age >65)",
        )
        text = rwe.to_embedding_text()
        assert "Yescarta" in text
        assert "CIBMTR" in text
        assert "42%" in text
        assert "elderly" in text


# ═══════════════════════════════════════════════════════════════════════
# PARAMETRIZED: ALL 10 MODELS PRODUCE NON-EMPTY EMBEDDING TEXT
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize(
    "model_cls,kwargs",
    [
        (CARTLiterature, {"id": "1", "title": "T", "text_chunk": "C", "year": 2024}),
        (ClinicalTrial, {"id": "NCT00000001", "title": "T", "text_summary": "S"}),
        (CARConstruct, {"id": "c1", "name": "N", "text_summary": "S", "target_antigen": "CD19"}),
        (AssayResult, {"id": "a1", "text_summary": "S"}),
        (ManufacturingRecord, {"id": "m1", "text_summary": "S"}),
        (SafetyRecord, {"id": "s1", "text_summary": "S"}),
        (BiomarkerRecord, {"id": "b1", "text_summary": "S"}),
        (RegulatoryRecord, {"id": "r1", "text_summary": "S"}),
        (SequenceRecord, {"id": "sq1", "text_summary": "S"}),
        (RealWorldRecord, {"id": "rw1", "text_summary": "S"}),
    ],
    ids=[
        "CARTLiterature",
        "ClinicalTrial",
        "CARConstruct",
        "AssayResult",
        "ManufacturingRecord",
        "SafetyRecord",
        "BiomarkerRecord",
        "RegulatoryRecord",
        "SequenceRecord",
        "RealWorldRecord",
    ],
)
def test_all_models_embedding_text(model_cls, kwargs):
    """Every collection model's to_embedding_text() returns a non-empty string."""
    instance = model_cls(**kwargs)
    text = instance.to_embedding_text()
    assert isinstance(text, str)
    assert len(text) > 0


# ═══════════════════════════════════════════════════════════════════════
# ENUM VALUES
# ═══════════════════════════════════════════════════════════════════════


class TestEnums:
    """Verify enum values match expected domain constants."""

    def test_cart_stage_values(self):
        """CARTStage has the 5 expected development stages."""
        values = [s.value for s in CARTStage]
        assert "target_id" in values
        assert "car_design" in values
        assert "vector_eng" in values
        assert "testing" in values
        assert "clinical" in values

    def test_source_type_values(self):
        """SourceType has 5 publication source types."""
        assert len(SourceType) == 5
        assert SourceType.PUBMED.value == "pubmed"
        assert SourceType.PATENT.value == "patent"

    def test_trial_phase_values(self):
        """TrialPhase covers Early Phase 1 through Phase 4 plus N/A."""
        assert TrialPhase.EARLY_1.value == "Early Phase 1"
        assert TrialPhase.PHASE_3.value == "Phase 3"
        assert TrialPhase.NA.value == "N/A"

    def test_trial_status_values(self):
        """TrialStatus covers recruitment lifecycle statuses."""
        assert TrialStatus.RECRUITING.value == "Recruiting"
        assert TrialStatus.COMPLETED.value == "Completed"
        assert TrialStatus.TERMINATED.value == "Terminated"

    def test_car_generation_values(self):
        """CARGeneration covers 1st through 4th plus armored/universal."""
        assert len(CARGeneration) == 6
        assert CARGeneration.SECOND.value == "2nd"
        assert CARGeneration.ARMORED.value == "armored"

    def test_assay_type_values(self):
        """AssayType includes cytotoxicity, cytokine, and in_vivo."""
        assert AssayType.CYTOTOXICITY.value == "cytotoxicity"
        assert AssayType.IN_VIVO.value == "in_vivo"
        assert AssayType.EXHAUSTION.value == "exhaustion"

    def test_process_step_values(self):
        """ProcessStep covers manufacturing lifecycle."""
        assert ProcessStep.TRANSDUCTION.value == "transduction"
        assert ProcessStep.CRYOPRESERVATION.value == "cryo"

    def test_fda_status_values(self):
        """FDAStatus includes approved through discontinued."""
        assert FDAStatus.APPROVED.value == "approved"
        assert FDAStatus.DISCONTINUED.value == "discontinued"

    def test_safety_event_type_values(self):
        """SafetyEventType covers CRS, ICANS, and other key toxicities."""
        assert SafetyEventType.CRS.value == "CRS"
        assert SafetyEventType.ICANS.value == "ICANS"
        assert SafetyEventType.SECONDARY_MALIGNANCY.value == "secondary_malignancy"

    def test_biomarker_type_values(self):
        """BiomarkerType includes predictive, prognostic, and resistance."""
        assert BiomarkerType.PREDICTIVE.value == "predictive"
        assert BiomarkerType.RESISTANCE.value == "resistance"

    def test_evidence_level_values(self):
        """EvidenceLevel has 3 tiers."""
        assert len(EvidenceLevel) == 3
        assert EvidenceLevel.VALIDATED.value == "validated"

    def test_regulatory_event_values(self):
        """RegulatoryEvent covers BLA through complete response."""
        assert RegulatoryEvent.BLA.value == "BLA"
        assert RegulatoryEvent.RMAT.value == "RMAT"
        assert RegulatoryEvent.COMPLETE_RESPONSE.value == "complete_response"

    def test_rwe_study_type_values(self):
        """RWEStudyType covers retrospective, registry, and meta-analysis."""
        assert RWEStudyType.RETROSPECTIVE.value == "retrospective"
        assert RWEStudyType.META_ANALYSIS.value == "meta_analysis"


# ═══════════════════════════════════════════════════════════════════════
# SEARCH RESULT MODELS
# ═══════════════════════════════════════════════════════════════════════


class TestSearchHit:
    """Tests for SearchHit creation and fields."""

    def test_create_search_hit(self):
        """SearchHit stores collection, id, score, text, and metadata."""
        hit = SearchHit(
            collection="Literature",
            id="12345678",
            score=0.85,
            text="CAR-T therapy study results.",
            metadata={"year": 2023},
        )
        assert hit.collection == "Literature"
        assert hit.score == 0.85
        assert hit.metadata["year"] == 2023

    def test_default_metadata_is_empty_dict(self):
        """SearchHit metadata defaults to an empty dict."""
        hit = SearchHit(collection="Trial", id="NCT00000001", score=0.5, text="test")
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
        assert "Literature" in grouped
        assert "Trial" in grouped
        assert len(grouped["Literature"]) == 1

    def test_empty_result(self):
        """An empty CrossCollectionResult has hit_count == 0."""
        result = CrossCollectionResult(query="empty")
        assert result.hit_count == 0
        assert result.hits_by_collection() == {}


class TestComparativeResult:
    """Tests for ComparativeResult creation and properties."""

    def test_create_comparative_result(self):
        """ComparativeResult links two CrossCollectionResults."""
        ev_a = CrossCollectionResult(
            query="CD19",
            hits=[SearchHit(collection="Literature", id="1", score=0.9, text="A")],
        )
        ev_b = CrossCollectionResult(
            query="BCMA",
            hits=[
                SearchHit(collection="Literature", id="2", score=0.85, text="B"),
                SearchHit(collection="Trial", id="NCT00000001", score=0.8, text="C"),
            ],
        )
        comp = ComparativeResult(
            query="Compare CD19 vs BCMA",
            entity_a="CD19",
            entity_b="BCMA",
            evidence_a=ev_a,
            evidence_b=ev_b,
        )
        assert comp.total_hits == 3
        assert comp.entity_a == "CD19"
        assert comp.entity_b == "BCMA"


# ═══════════════════════════════════════════════════════════════════════
# AGENT QUERY MODEL
# ═══════════════════════════════════════════════════════════════════════


class TestAgentQuery:
    """Tests for the AgentQuery input model."""

    def test_create_with_required_only(self):
        """AgentQuery requires only the question field."""
        query = AgentQuery(question="What causes CRS?")
        assert query.question == "What causes CRS?"
        assert query.target_antigen is None
        assert query.cart_stage is None
        assert query.include_genomic is True

    def test_create_with_all_fields(self):
        """AgentQuery accepts all optional fields."""
        query = AgentQuery(
            question="CD19 CAR-T resistance",
            target_antigen="CD19",
            cart_stage=CARTStage.CLINICAL,
            include_genomic=False,
        )
        assert query.target_antigen == "CD19"
        assert query.cart_stage == CARTStage.CLINICAL
        assert query.include_genomic is False
