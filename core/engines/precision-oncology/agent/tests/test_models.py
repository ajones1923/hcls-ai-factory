"""
Tests for Oncology Intelligence Agent Pydantic models.
====================================================
Validates all enums, domain models, embedding text generation,
and agent I/O types defined in src/models.py.
"""

import sys
from datetime import datetime
from pathlib import Path

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.models import (
    AgentQuery,
    AgentResponse,
    BiomarkerType,
    CancerType,
    CaseSnapshot,
    CrossCollectionResult,
    EvidenceLevel,
    GuidelineOrg,
    MTBPacket,
    OncologyBiomarker,
    OncologyGuideline,
    OncologyTherapy,
    OncologyTrial,
    OncologyVariant,
    OutcomeRecord,
    PathwayName,
    ResistanceMechanism,
    ResponseCategory,
    SearchHit,
    SourceType,
    TherapyCategory,
    TrialPhase,
    TrialStatus,
    VariantType,
)


# ═══════════════════════════════════════════════════════════════════════════
# Enum Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCancerType:
    """Verify CancerType enum members."""

    def test_nsclc(self):
        assert CancerType.NSCLC.value == "nsclc"

    def test_breast(self):
        assert CancerType.BREAST.value == "breast"

    def test_melanoma(self):
        assert CancerType.MELANOMA.value == "melanoma"

    def test_colorectal(self):
        assert CancerType.COLORECTAL.value == "colorectal"

    def test_other(self):
        assert CancerType.OTHER.value == "other"

    def test_expected_members(self):
        names = {m.name for m in CancerType}
        expected = {
            "NSCLC", "SCLC", "BREAST", "COLORECTAL", "MELANOMA",
            "PANCREATIC", "OVARIAN", "PROSTATE", "RENAL", "BLADDER",
            "HEAD_NECK", "HEPATOCELLULAR", "GASTRIC", "GLIOBLASTOMA",
            "AML", "CML", "ALL", "CLL", "DLBCL", "MULTIPLE_MYELOMA", "OTHER",
        }
        assert expected.issubset(names)


class TestVariantType:
    """Verify VariantType enum members."""

    def test_snv(self):
        assert VariantType.SNV.value == "snv"

    def test_fusion(self):
        assert VariantType.FUSION.value == "fusion"

    def test_cnv_amp(self):
        assert VariantType.CNV_AMP.value == "cnv_amplification"

    def test_expected_members(self):
        names = {m.name for m in VariantType}
        expected = {"SNV", "INDEL", "CNV_AMP", "CNV_DEL", "FUSION", "REARRANGEMENT", "SV"}
        assert expected == names


class TestEvidenceLevel:
    """Verify EvidenceLevel enum members."""

    def test_level_a(self):
        assert EvidenceLevel.LEVEL_A.value == "A"

    def test_level_e(self):
        assert EvidenceLevel.LEVEL_E.value == "E"

    def test_expected_members(self):
        names = {m.name for m in EvidenceLevel}
        assert names == {"LEVEL_A", "LEVEL_B", "LEVEL_C", "LEVEL_D", "LEVEL_E"}


class TestTherapyCategory:
    """Verify TherapyCategory enum members."""

    def test_targeted(self):
        assert TherapyCategory.TARGETED.value == "targeted"

    def test_immunotherapy(self):
        assert TherapyCategory.IMMUNOTHERAPY.value == "immunotherapy"

    def test_expected_members(self):
        names = {m.name for m in TherapyCategory}
        expected = {
            "TARGETED", "IMMUNOTHERAPY", "CHEMOTHERAPY", "HORMONAL",
            "COMBINATION", "RADIOTHERAPY", "CELL_THERAPY",
            "ADC", "BISPECIFIC",
        }
        assert expected == names


class TestTrialPhase:
    """Verify TrialPhase enum members."""

    def test_phase_3(self):
        assert TrialPhase.PHASE_3.value == "Phase 3"

    def test_na(self):
        assert TrialPhase.NA.value == "N/A"

    def test_expected_members(self):
        names = {m.name for m in TrialPhase}
        expected = {
            "EARLY_PHASE_1", "PHASE_1", "PHASE_1_2", "PHASE_2",
            "PHASE_2_3", "PHASE_3", "PHASE_4", "NA",
        }
        assert expected == names


class TestTrialStatus:
    """Verify TrialStatus enum members."""

    def test_recruiting(self):
        assert TrialStatus.RECRUITING.value == "Recruiting"

    def test_completed(self):
        assert TrialStatus.COMPLETED.value == "Completed"

    def test_expected_members(self):
        names = {m.name for m in TrialStatus}
        expected = {
            "NOT_YET_RECRUITING", "RECRUITING", "ENROLLING_BY_INVITATION",
            "ACTIVE_NOT_RECRUITING", "SUSPENDED", "TERMINATED",
            "COMPLETED", "WITHDRAWN", "UNKNOWN",
        }
        assert expected == names


class TestResponseCategory:
    """Verify ResponseCategory enum members."""

    def test_cr(self):
        assert ResponseCategory.CR.value == "complete_response"

    def test_pd(self):
        assert ResponseCategory.PD.value == "progressive_disease"

    def test_expected_members(self):
        names = {m.name for m in ResponseCategory}
        assert names == {"CR", "PR", "SD", "PD", "NE"}


class TestBiomarkerType:
    """Verify BiomarkerType enum members."""

    def test_predictive(self):
        assert BiomarkerType.PREDICTIVE.value == "predictive"

    def test_resistance(self):
        assert BiomarkerType.RESISTANCE.value == "resistance"

    def test_expected_members(self):
        names = {m.name for m in BiomarkerType}
        expected = {
            "PREDICTIVE", "PROGNOSTIC", "DIAGNOSTIC",
            "MONITORING", "RESISTANCE", "PHARMACODYNAMIC",
            "SCREENING", "THERAPEUTIC_SELECTION",
        }
        assert expected == names


class TestPathwayName:
    """Verify PathwayName enum members."""

    def test_mapk(self):
        assert PathwayName.MAPK.value == "mapk"

    def test_pi3k(self):
        assert PathwayName.PI3K_AKT_MTOR.value == "pi3k_akt_mtor"

    def test_expected_members(self):
        names = {m.name for m in PathwayName}
        expected = {
            "MAPK", "PI3K_AKT_MTOR", "DDR", "CELL_CYCLE", "APOPTOSIS",
            "WNT", "NOTCH", "HEDGEHOG", "JAK_STAT", "ANGIOGENESIS",
            "HIPPO", "NF_KB", "TGF_BETA",
        }
        assert expected == names


class TestGuidelineOrg:
    """Verify GuidelineOrg enum members."""

    def test_nccn(self):
        assert GuidelineOrg.NCCN.value == "NCCN"

    def test_esmo(self):
        assert GuidelineOrg.ESMO.value == "ESMO"

    def test_expected_members(self):
        names = {m.name for m in GuidelineOrg}
        assert names == {"NCCN", "ESMO", "ASCO", "WHO", "CAP_AMP", "FDA", "EMA", "AACR"}


class TestSourceType:
    """Verify SourceType enum members."""

    def test_pubmed(self):
        assert SourceType.PUBMED.value == "pubmed"

    def test_expected_members(self):
        names = {m.name for m in SourceType}
        assert names == {"PUBMED", "PMC", "PREPRINT", "MANUAL"}


# ═══════════════════════════════════════════════════════════════════════════
# Domain Model Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestOncologyVariant:
    """Test OncologyVariant creation and to_embedding_text()."""

    @pytest.fixture
    def variant(self):
        return OncologyVariant(
            id="CIViC:12",
            gene="EGFR",
            variant_name="L858R",
            variant_type=VariantType.SNV,
            cancer_type=CancerType.NSCLC,
            evidence_level=EvidenceLevel.LEVEL_A,
            drugs=["osimertinib", "erlotinib"],
            text_summary="EGFR L858R is a sensitising mutation.",
            clinical_significance="Pathogenic",
        )

    def test_creation(self, variant):
        assert variant.gene == "EGFR"
        assert variant.variant_name == "L858R"
        assert variant.variant_type == VariantType.SNV
        assert variant.evidence_level == EvidenceLevel.LEVEL_A

    def test_embedding_text(self, variant):
        text = variant.to_embedding_text()
        assert "EGFR" in text
        assert "L858R" in text
        assert "snv" in text
        assert "A" in text
        assert "osimertinib" in text
        assert "nsclc" in text
        assert "Pathogenic" in text


class TestOncologyBiomarker:
    """Test OncologyBiomarker creation and to_embedding_text()."""

    @pytest.fixture
    def biomarker(self):
        return OncologyBiomarker(
            id="BM:TMB-H",
            name="Tumor Mutational Burden High",
            biomarker_type=BiomarkerType.PREDICTIVE,
            cancer_types=[CancerType.NSCLC, CancerType.MELANOMA],
            predictive_value="Pembrolizumab response",
            testing_method="NGS panel",
            clinical_cutoff=">=10 mut/Mb",
            text_summary="TMB-H predicts immunotherapy benefit.",
            evidence_level=EvidenceLevel.LEVEL_A,
        )

    def test_creation(self, biomarker):
        assert biomarker.name == "Tumor Mutational Burden High"
        assert biomarker.biomarker_type == BiomarkerType.PREDICTIVE
        assert CancerType.NSCLC in biomarker.cancer_types

    def test_embedding_text(self, biomarker):
        text = biomarker.to_embedding_text()
        assert "Tumor Mutational Burden" in text
        assert "predictive" in text
        assert "NGS panel" in text
        assert ">=10 mut/Mb" in text


class TestOncologyTherapy:
    """Test OncologyTherapy creation and to_embedding_text()."""

    @pytest.fixture
    def therapy(self):
        return OncologyTherapy(
            id="TX:osimertinib",
            drug_name="osimertinib",
            category=TherapyCategory.TARGETED,
            targets=["EGFR"],
            approved_indications=["EGFR-mutant NSCLC"],
            evidence_level=EvidenceLevel.LEVEL_A,
            text_summary="Third-generation EGFR TKI.",
            mechanism_of_action="Irreversible EGFR C797 binding",
        )

    def test_creation(self, therapy):
        assert therapy.drug_name == "osimertinib"
        assert therapy.category == TherapyCategory.TARGETED
        assert "EGFR" in therapy.targets

    def test_embedding_text(self, therapy):
        text = therapy.to_embedding_text()
        assert "osimertinib" in text
        assert "targeted" in text
        assert "EGFR" in text
        assert "Irreversible" in text


class TestOncologyTrial:
    """Test OncologyTrial creation and to_embedding_text()."""

    @pytest.fixture
    def trial(self):
        return OncologyTrial(
            id="NCT02628067",
            title="KEYNOTE-158",
            text_summary="Phase II basket trial of pembrolizumab in MSI-H tumors.",
            phase=TrialPhase.PHASE_2,
            status=TrialStatus.ACTIVE_NOT_RECRUITING,
            cancer_types=[CancerType.COLORECTAL, CancerType.GASTRIC],
            biomarker_criteria=["MSI-H", "TMB-H"],
            outcome_summary="ORR 34.3%",
        )

    def test_creation(self, trial):
        assert trial.id == "NCT02628067"
        assert trial.phase == TrialPhase.PHASE_2
        assert trial.status == TrialStatus.ACTIVE_NOT_RECRUITING
        assert CancerType.COLORECTAL in trial.cancer_types

    def test_embedding_text(self, trial):
        text = trial.to_embedding_text()
        assert "KEYNOTE-158" in text
        assert "Phase 2" in text
        assert "Active, not recruiting" in text
        assert "MSI-H" in text
        assert "ORR 34.3%" in text


class TestOncologyGuideline:
    """Test OncologyGuideline creation and to_embedding_text()."""

    @pytest.fixture
    def guideline(self):
        return OncologyGuideline(
            id="GL:NCCN-NSCLC-2025.2",
            org=GuidelineOrg.NCCN,
            cancer_type=CancerType.NSCLC,
            version="2.2025",
            year=2025,
            key_recommendations=[
                "Molecular testing for all non-squamous NSCLC",
                "Osimertinib first-line for EGFR-mutant NSCLC",
            ],
            text_summary="NCCN NSCLC guidelines version 2.2025.",
            evidence_level=EvidenceLevel.LEVEL_A,
        )

    def test_creation(self, guideline):
        assert guideline.org == GuidelineOrg.NCCN
        assert guideline.cancer_type == CancerType.NSCLC
        assert guideline.year == 2025
        assert len(guideline.key_recommendations) == 2

    def test_embedding_text(self, guideline):
        text = guideline.to_embedding_text()
        assert "NCCN" in text
        assert "nsclc" in text
        assert "2.2025" in text
        assert "2025" in text
        assert "Molecular testing" in text


class TestResistanceMechanism:
    """Test ResistanceMechanism creation and to_embedding_text()."""

    @pytest.fixture
    def resistance(self):
        return ResistanceMechanism(
            id="RM:EGFR-T790M",
            primary_therapy="erlotinib",
            gene="EGFR",
            mechanism="T790M gatekeeper mutation",
            bypass_pathway="MAPK reactivation",
            alternative_therapies=["osimertinib"],
            text_summary="T790M confers resistance to first-gen EGFR TKIs.",
        )

    def test_creation(self, resistance):
        assert resistance.gene == "EGFR"
        assert resistance.mechanism == "T790M gatekeeper mutation"
        assert "osimertinib" in resistance.alternative_therapies

    def test_embedding_text(self, resistance):
        text = resistance.to_embedding_text()
        assert "erlotinib" in text
        assert "EGFR" in text
        assert "T790M" in text
        assert "MAPK" in text
        assert "osimertinib" in text


class TestOutcomeRecord:
    """Test OutcomeRecord creation and to_embedding_text()."""

    @pytest.fixture
    def outcome(self):
        return OutcomeRecord(
            id="OUT:001",
            case_id="CASE:XYZ",
            therapy="osimertinib",
            cancer_type=CancerType.NSCLC,
            response=ResponseCategory.PR,
            duration_months=14.2,
            toxicities=["rash", "diarrhea"],
            biomarkers_at_baseline={"EGFR": "L858R", "PD-L1_TPS": "40"},
            text_summary="Partial response on osimertinib for 14 months.",
        )

    def test_creation(self, outcome):
        assert outcome.therapy == "osimertinib"
        assert outcome.response == ResponseCategory.PR
        assert outcome.duration_months == 14.2

    def test_embedding_text(self, outcome):
        text = outcome.to_embedding_text()
        assert "osimertinib" in text
        assert "nsclc" in text
        assert "partial_response" in text
        assert "14.2 months" in text
        assert "rash" in text
        assert "EGFR=L858R" in text


class TestCaseSnapshot:
    """Test CaseSnapshot creation."""

    def test_creation_with_case_id(self):
        case = CaseSnapshot(
            case_id="CASE:001",
            patient_id="PT-001",
            cancer_type="nsclc",
            stage="IV",
            variants=["EGFR L858R", "TP53 R175H"],
            biomarkers={"PD-L1_TPS": "80", "TMB": "12.5"},
            prior_therapies=["carboplatin + pemetrexed"],
            text_summary="Stage IV NSCLC with EGFR L858R.",
        )
        assert case.patient_id == "PT-001"
        assert case.cancer_type == "nsclc"
        assert case.stage == "IV"
        assert len(case.variants) == 2
        assert "PD-L1_TPS" in case.biomarkers

    def test_creation_with_dict_variants(self):
        """CaseSnapshot should accept variant dicts (used by case_manager)."""
        case = CaseSnapshot(
            case_id="CASE:002",
            patient_id="PT-002",
            cancer_type="NSCLC",
            stage="IV",
            variants=[{"gene": "EGFR", "variant": "L858R", "actionability": "A"}],
            biomarkers={"TMB": 12.5},
            text_summary="NSCLC case.",
        )
        assert len(case.variants) == 1
        assert case.variants[0]["gene"] == "EGFR"

    def test_embedding_text_with_string_variants(self):
        case = CaseSnapshot(
            case_id="CASE:003",
            patient_id="PT-003",
            cancer_type="nsclc",
            variants=["EGFR L858R"],
            text_summary="Test case.",
        )
        text = case.to_embedding_text()
        assert "EGFR L858R" in text

    def test_embedding_text_with_dict_variants(self):
        case = CaseSnapshot(
            case_id="CASE:004",
            patient_id="PT-004",
            cancer_type="nsclc",
            variants=[{"gene": "BRAF", "variant": "V600E"}],
            text_summary="Test case.",
        )
        text = case.to_embedding_text()
        assert "BRAF" in text

    def test_default_text_summary(self):
        case = CaseSnapshot(
            case_id="CASE:005",
            patient_id="PT-005",
            cancer_type="breast",
        )
        assert case.text_summary == ""


class TestMTBPacket:
    """Test MTBPacket creation."""

    def test_creation(self):
        packet = MTBPacket(
            case_id="CASE:001",
            patient_summary="Stage IV NSCLC with EGFR L858R.",
            variant_table=[{"gene": "EGFR", "variant": "L858R", "evidence_level": "A"}],
            evidence_table=[{"gene": "EGFR", "citations": []}],
            therapy_ranking=[{"drug_name": "osimertinib", "rank": 1}],
            trial_matches=[{"trial_id": "NCT00000001", "match_score": 0.9}],
            open_questions=["VUS in TP53"],
            citations=["PMID:12345678"],
        )
        assert packet.case_id == "CASE:001"
        assert len(packet.variant_table) == 1
        assert len(packet.therapy_ranking) == 1
        assert isinstance(packet.generated_at, datetime)


# ═══════════════════════════════════════════════════════════════════════════
# Search & Agent Model Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestSearchHit:
    """Test SearchHit creation."""

    def test_creation(self):
        hit = SearchHit(
            collection="onco_variants",
            id="CIViC:12",
            score=0.91,
            text="EGFR L858R is a sensitising mutation.",
            metadata={"gene": "EGFR"},
        )
        assert hit.collection == "onco_variants"
        assert hit.score == 0.91
        assert hit.metadata["gene"] == "EGFR"

    def test_default_metadata(self):
        hit = SearchHit(collection="onco_literature", id="x", score=0.5, text="test")
        assert hit.metadata == {}


class TestCrossCollectionResult:
    """Test CrossCollectionResult creation and helper methods."""

    def test_creation(self):
        result = CrossCollectionResult(
            query="BRAF V600E therapy",
            hits=[
                SearchHit(collection="onco_variants", id="v1", score=0.9, text="t1"),
                SearchHit(collection="onco_therapies", id="t1", score=0.8, text="t2"),
            ],
            total_collections_searched=2,
            search_time_ms=15.0,
        )
        assert result.hit_count == 2
        assert result.query == "BRAF V600E therapy"

    def test_hits_by_collection(self):
        result = CrossCollectionResult(
            query="test",
            hits=[
                SearchHit(collection="onco_variants", id="v1", score=0.9, text="a"),
                SearchHit(collection="onco_variants", id="v2", score=0.8, text="b"),
                SearchHit(collection="onco_therapies", id="t1", score=0.7, text="c"),
            ],
        )
        grouped = result.hits_by_collection()
        assert len(grouped["onco_variants"]) == 2
        assert len(grouped["onco_therapies"]) == 1


class TestAgentQuery:
    """Test AgentQuery creation."""

    def test_basic_query(self):
        q = AgentQuery(question="What is EGFR L858R?")
        assert q.question == "What is EGFR L858R?"
        assert q.include_genomic is True

    def test_with_filters(self):
        q = AgentQuery(
            question="BRAF therapy in melanoma",
            cancer_type=CancerType.MELANOMA,
            gene="BRAF",
        )
        assert q.cancer_type == CancerType.MELANOMA
        assert q.gene == "BRAF"


class TestAgentResponse:
    """Test AgentResponse creation."""

    def test_creation(self):
        evidence = CrossCollectionResult(query="test", hits=[])
        resp = AgentResponse(
            question="What is EGFR?",
            answer="EGFR is a receptor tyrosine kinase.",
            evidence=evidence,
            knowledge_used=["ACTIONABLE_TARGETS"],
        )
        assert resp.question == "What is EGFR?"
        assert resp.answer.startswith("EGFR")
        assert isinstance(resp.timestamp, datetime)
        assert "ACTIONABLE_TARGETS" in resp.knowledge_used
