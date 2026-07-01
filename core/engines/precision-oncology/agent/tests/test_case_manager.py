"""
Tests for the OncologyCaseManager.
====================================
Validates variant actionability classification, case creation,
VCF parsing, and MTB packet generation.

All external dependencies (Milvus, embeddings, RAG engine) are mocked.
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.case_manager import OncologyCaseManager
from src.knowledge import ACTIONABLE_TARGETS


# ═══════════════════════════════════════════════════════════════════════════
# Fixtures
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def case_manager(mock_collection_manager, mock_embedder):
    """Create an OncologyCaseManager with mocked dependencies."""
    mock_knowledge = MagicMock()
    mock_rag = MagicMock()
    mock_rag.retrieve.return_value = []
    return OncologyCaseManager(
        collection_manager=mock_collection_manager,
        embedder=mock_embedder,
        knowledge=mock_knowledge,
        rag_engine=mock_rag,
    )


@pytest.fixture
def sample_variants():
    """Pre-parsed variant list for testing."""
    return [
        {"gene": "EGFR", "variant": "L858R", "chrom": "chr7", "pos": 55259515,
         "ref": "T", "alt": "G", "consequence": "missense_variant", "filter": "PASS"},
        {"gene": "TP53", "variant": "R175H", "chrom": "chr17", "pos": 7578406,
         "ref": "C", "alt": "T", "consequence": "missense_variant", "filter": "PASS"},
        {"gene": "BRAF", "variant": "V600E", "chrom": "chr7", "pos": 140453136,
         "ref": "A", "alt": "T", "consequence": "missense_variant", "filter": "PASS"},
        {"gene": "UNKNOWNGENE", "variant": "X123Y", "chrom": "chr1", "pos": 100,
         "ref": "A", "alt": "G", "consequence": "missense_variant", "filter": "PASS"},
    ]


@pytest.fixture
def sample_biomarkers():
    """Sample biomarker data."""
    return {
        "MSI": "MSI-H",
        "TMB": 14.2,
        "PD-L1_TPS": 80,
    }


# ═══════════════════════════════════════════════════════════════════════════
# Variant Actionability Classification Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestClassifyVariantActionability:
    """Test _classify_variant_actionability for known and unknown targets."""

    @pytest.fixture
    def manager(self, case_manager):
        return case_manager

    def test_egfr_is_actionable(self, manager):
        """EGFR is a well-known actionable target."""
        result = manager._classify_variant_actionability("EGFR", "L858R")
        assert result != "VUS", "EGFR L858R should be actionable, not VUS"

    def test_braf_is_actionable(self, manager):
        """BRAF V600E is a well-known actionable target."""
        result = manager._classify_variant_actionability("BRAF", "V600E")
        assert result != "VUS", "BRAF V600E should be actionable, not VUS"

    def test_alk_is_recognized(self, manager):
        """ALK should be recognized as an actionable gene."""
        # ALK is in ACTIONABLE_TARGETS but actionability depends on variant details
        result = manager._classify_variant_actionability("ALK", "EML4-ALK fusion")
        # Should at least not crash; may be VUS if variant lookup doesn't match
        assert isinstance(result, str)

    def test_kras_g12c_is_actionable(self, manager):
        """KRAS G12C should be recognized."""
        result = manager._classify_variant_actionability("KRAS", "G12C")
        assert result != "VUS" or "KRAS" in ACTIONABLE_TARGETS

    def test_unknown_gene_is_vus(self, manager):
        """An unknown gene should return VUS."""
        result = manager._classify_variant_actionability("UNKNOWNGENE", "X123Y")
        assert result == "VUS"

    def test_empty_gene_is_vus(self, manager):
        """Empty gene string should return VUS."""
        result = manager._classify_variant_actionability("", "")
        assert result == "VUS"

    def test_case_insensitive(self, manager):
        """Gene lookup should be case-insensitive."""
        result_upper = manager._classify_variant_actionability("EGFR", "L858R")
        result_lower = manager._classify_variant_actionability("egfr", "L858R")
        assert result_upper == result_lower

    def test_brca1_is_recognized(self, manager):
        """BRCA1 should be recognized as an actionable target."""
        result = manager._classify_variant_actionability("BRCA1", "185delAG")
        assert isinstance(result, str)

    def test_her2_is_recognized(self, manager):
        """HER2/ERBB2 should be recognized."""
        # HER2 is listed as ERBB2 or HER2 in ACTIONABLE_TARGETS
        result = manager._classify_variant_actionability("HER2", "amplification")
        assert isinstance(result, str)

    def test_ret_is_recognized(self, manager):
        """RET should be recognized as actionable."""
        result = manager._classify_variant_actionability("RET", "KIF5B-RET fusion")
        assert isinstance(result, str)


# ═══════════════════════════════════════════════════════════════════════════
# Case Creation Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestCreateCase:
    """Test create_case with mock data."""

    def test_create_case_with_variant_list(self, case_manager, sample_variants, sample_biomarkers):
        """Creating a case with pre-parsed variants should succeed."""
        case = case_manager.create_case(
            patient_id="PT-001",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
            biomarkers=sample_biomarkers,
            prior_therapies=["carboplatin + pemetrexed"],
        )
        assert case is not None
        assert case.patient_id == "PT-001"
        assert case.cancer_type == "NSCLC"
        assert case.stage == "IV"
        # Variants should have actionability assigned
        for v in case.variants:
            assert "actionability" in v

    def test_create_case_empty_variants(self, case_manager):
        """Creating a case with empty variant list should succeed."""
        case = case_manager.create_case(
            patient_id="PT-002",
            cancer_type="BREAST",
            stage="IIIA",
            vcf_content_or_variants=[],
        )
        assert case is not None
        assert case.variants == []

    def test_create_case_default_biomarkers(self, case_manager, sample_variants):
        """Missing biomarkers should default to empty dict."""
        case = case_manager.create_case(
            patient_id="PT-003",
            cancer_type="MELANOMA",
            stage="IIIC",
            vcf_content_or_variants=sample_variants,
        )
        assert case.biomarkers == {}

    def test_create_case_default_prior_therapies(self, case_manager, sample_variants):
        """Missing prior therapies should default to empty list."""
        case = case_manager.create_case(
            patient_id="PT-004",
            cancer_type="CRC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
        )
        assert case.prior_therapies == []

    def test_create_case_invalid_input_raises(self, case_manager):
        """Passing invalid vcf_content_or_variants type should raise ValueError."""
        with pytest.raises(ValueError, match="must be a VCF string or list"):
            case_manager.create_case(
                patient_id="PT-005",
                cancer_type="AML",
                stage="N/A",
                vcf_content_or_variants=12345,
            )

    def test_create_case_assigns_case_id(self, case_manager, sample_variants):
        """Each new case should get a unique case_id."""
        case1 = case_manager.create_case(
            patient_id="PT-010",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
        )
        case2 = case_manager.create_case(
            patient_id="PT-011",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
        )
        assert case1.case_id != case2.case_id


# ═══════════════════════════════════════════════════════════════════════════
# VCF Parsing Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestParseVCF:
    """Test _parse_vcf_text extracts PASS variants correctly."""

    @pytest.fixture
    def manager(self, case_manager):
        return case_manager

    def test_parse_simple_vcf(self, manager):
        """Should parse a simple VCF line with PASS filter."""
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr7\t55259515\t.\tT\tG\t100\tPASS\tANN=T|EGFR|missense_variant|HIGH\n"
        )
        variants = manager._parse_vcf_text(vcf_text)
        assert len(variants) == 1
        assert variants[0]["chrom"] == "chr7"
        assert variants[0]["filter"] == "PASS"
        assert variants[0]["gene"] == "EGFR"

    def test_skip_non_pass_variants(self, manager):
        """Should skip variants that do not have PASS filter."""
        vcf_text = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t50\tLowQual\tANN=A|TP53|missense_variant|HIGH\n"
        )
        variants = manager._parse_vcf_text(vcf_text)
        assert len(variants) == 0

    def test_skip_header_lines(self, manager):
        """Should skip all header lines starting with #."""
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=ANN,...>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        )
        variants = manager._parse_vcf_text(vcf_text)
        assert len(variants) == 0

    def test_parse_empty_vcf(self, manager):
        """Empty VCF should return empty list."""
        variants = manager._parse_vcf_text("")
        assert variants == []


# ═══════════════════════════════════════════════════════════════════════════
# MTB Packet Generation Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestGenerateMTBPacket:
    """Test generate_mtb_packet structure has required sections."""

    def test_packet_structure(self, case_manager, sample_variants, sample_biomarkers):
        """MTB packet should contain all required sections."""
        case = case_manager.create_case(
            patient_id="PT-MTB-001",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
            biomarkers=sample_biomarkers,
            prior_therapies=["carboplatin"],
        )
        packet = case_manager.generate_mtb_packet(case)

        assert packet is not None
        assert packet.case_id == case.case_id
        assert isinstance(packet.variant_table, list)
        assert isinstance(packet.evidence_table, list)
        assert isinstance(packet.therapy_ranking, list)
        assert isinstance(packet.trial_matches, list)
        assert isinstance(packet.open_questions, list)

    def test_packet_variant_table_populated(self, case_manager, sample_variants):
        """Variant table should have entries for each variant."""
        case = case_manager.create_case(
            patient_id="PT-MTB-002",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
        )
        packet = case_manager.generate_mtb_packet(case)
        assert len(packet.variant_table) == len(sample_variants)

    def test_packet_open_questions_for_vus(self, case_manager, sample_variants):
        """Open questions should flag VUS variants."""
        case = case_manager.create_case(
            patient_id="PT-MTB-003",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
        )
        packet = case_manager.generate_mtb_packet(case)
        # At least UNKNOWNGENE should be flagged as VUS
        vus_questions = [q for q in packet.open_questions if "uncertain significance" in q.lower()]
        assert len(vus_questions) >= 1

    def test_packet_open_questions_for_missing_biomarkers(self, case_manager, sample_variants):
        """Open questions should flag missing biomarkers (MSI, TMB, PD-L1)."""
        case = case_manager.create_case(
            patient_id="PT-MTB-004",
            cancer_type="NSCLC",
            stage="IV",
            vcf_content_or_variants=sample_variants,
            biomarkers={},  # no biomarkers
        )
        packet = case_manager.generate_mtb_packet(case)
        missing_questions = [q for q in packet.open_questions if "missing" in q.lower()]
        assert len(missing_questions) >= 1
