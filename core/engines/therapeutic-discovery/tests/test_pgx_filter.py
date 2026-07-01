"""
Tests for Pharmacogenomic Drug Candidate Filter.

Validates PGx-based filtering and re-ranking of drug candidates
based on patient metabolizer phenotypes. These tests ensure the
closed-loop feedback from Biomarker Agent to Drug Discovery Pipeline
works correctly.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.pgx_filter import (
    CYP_SUBSTRATE_MAP,
    PHENOTYPE_RISK,
    MetabolismRisk,
    PGxDrugFilter,
    PGxFilterResult,
    PGxFilterSummary,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def poor_metabolizer_profile():
    """Patient who is a CYP2D6 poor metabolizer."""
    return {
        "gene_results": [
            {"gene": "CYP2D6", "phenotype": "Poor Metabolizer", "diplotype": "*4/*4"},
            {"gene": "CYP2C19", "phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        ]
    }


@pytest.fixture
def normal_metabolizer_profile():
    """Patient with all normal metabolizer phenotypes."""
    return {
        "gene_results": [
            {"gene": "CYP2D6", "phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
            {"gene": "CYP2C19", "phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
            {"gene": "CYP2C9", "phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        ]
    }


@pytest.fixture
def ultra_rapid_profile():
    """Patient who is a CYP2D6 ultra-rapid metabolizer."""
    return {
        "gene_results": [
            {"gene": "CYP2D6", "phenotype": "Ultra Rapid Metabolizer", "diplotype": "*1/*1xN"},
        ]
    }


@pytest.fixture
def multi_risk_profile():
    """Patient with multiple at-risk phenotypes."""
    return {
        "gene_results": [
            {"gene": "CYP2D6", "phenotype": "Poor Metabolizer", "diplotype": "*4/*4"},
            {"gene": "CYP2C19", "phenotype": "Ultra Rapid Metabolizer", "diplotype": "*17/*17"},
            {"gene": "CYP2C9", "phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        ]
    }


@pytest.fixture
def sample_candidates():
    """A set of drug candidates for testing."""
    return [
        {
            "name": "VCP-001",
            "smiles": "CC(C)C1=CC=CC=C1",
            "metabolism_enzymes": ["CYP2D6"],
        },
        {
            "name": "VCP-002",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP2C19"],
        },
        {
            "name": "VCP-003",
            "smiles": "CCCO",
            "metabolism_enzymes": [],  # Unknown metabolism
        },
        {
            "name": "VCP-004",
            "smiles": "C1CCCCC1",
            "metabolism_enzymes": ["CYP3A5"],
        },
        {
            "name": "VCP-005",
            "smiles": "C1=CC=CC=C1",
            "metabolism_enzymes": ["CYP2C9"],
        },
    ]


# ---------------------------------------------------------------------------
# Test MetabolismRisk enum
# ---------------------------------------------------------------------------

class TestMetabolismRisk:
    """Tests for the MetabolismRisk enum."""

    def test_risk_values(self):
        """Verify all expected risk levels exist."""
        assert MetabolismRisk.SAFE == "safe"
        assert MetabolismRisk.CAUTION == "caution"
        assert MetabolismRisk.CONTRAINDICATED == "contraindicated"
        assert MetabolismRisk.UNKNOWN == "unknown"

    def test_risk_is_string_enum(self):
        """MetabolismRisk should be usable as a string."""
        assert str(MetabolismRisk.SAFE) == "MetabolismRisk.SAFE"
        assert MetabolismRisk.SAFE.value == "safe"


# ---------------------------------------------------------------------------
# Test PGxDrugFilter initialization
# ---------------------------------------------------------------------------

class TestPGxDrugFilterInit:
    """Tests for PGxDrugFilter initialization."""

    def test_init_with_profile(self, poor_metabolizer_profile):
        """Filter initializes correctly with a PGx profile."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        assert f.pgx_profile == poor_metabolizer_profile
        assert "CYP2D6" in f._phenotype_map
        assert f._phenotype_map["CYP2D6"] == "poor_metabolizer"
        assert f._phenotype_map["CYP2C19"] == "normal_metabolizer"

    def test_init_without_profile(self):
        """Filter initializes correctly without a PGx profile."""
        f = PGxDrugFilter()
        assert f.pgx_profile == {}
        assert f._phenotype_map == {}

    def test_init_with_none(self):
        """Filter initializes correctly with None."""
        f = PGxDrugFilter(None)
        assert f.pgx_profile == {}
        assert f._phenotype_map == {}

    def test_init_with_empty_gene_results(self):
        """Filter handles empty gene_results list."""
        f = PGxDrugFilter({"gene_results": []})
        assert f._phenotype_map == {}

    def test_phenotype_normalization(self):
        """Phenotype strings are normalized to lowercase with underscores."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Ultra-Rapid Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        assert f._phenotype_map["CYP2D6"] == "ultra_rapid_metabolizer"


# ---------------------------------------------------------------------------
# Test poor metabolizer filtering
# ---------------------------------------------------------------------------

class TestPoorMetabolizerFiltering:
    """Tests for filtering with a poor metabolizer profile."""

    def test_poor_metabolizer_flags_cyp2d6_substrate(self, poor_metabolizer_profile):
        """A CYP2D6 poor metabolizer should flag CYP2D6 substrates."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidate = {
            "name": "TestDrug",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP2D6"],
        }
        result = f.assess_candidate(candidate, rank=1)

        assert result.metabolism_risk == MetabolismRisk.CONTRAINDICATED
        assert "CYP2D6" in result.affected_enzymes
        assert result.rank_adjustment == 100
        assert "AVOID" in result.recommendation

    def test_poor_metabolizer_safe_for_other_enzymes(self, poor_metabolizer_profile):
        """A CYP2D6 poor metabolizer should be safe for CYP3A5 substrates."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidate = {
            "name": "TestDrug",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP3A5"],
        }
        result = f.assess_candidate(candidate, rank=1)

        # CYP3A5 is not flagged (patient is only a CYP2D6 PM)
        assert result.metabolism_risk == MetabolismRisk.SAFE
        assert result.rank_adjustment == 0

    def test_poor_metabolizer_unknown_enzymes_conservative(self, poor_metabolizer_profile):
        """With no enzyme info, poor metabolizer should conservatively flag."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidate = {
            "name": "NovelDrug",
            "smiles": "CCCO",
            "metabolism_enzymes": [],  # Unknown metabolism
        }
        result = f.assess_candidate(candidate, rank=1)

        # Should conservatively flag due to CYP2D6 poor metabolizer
        assert result.metabolism_risk == MetabolismRisk.CONTRAINDICATED
        assert "CYP2D6" in result.affected_enzymes

    def test_poor_metabolizer_prodrug(self, poor_metabolizer_profile):
        """A poor metabolizer with a prodrug should be contraindicated."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidate = {
            "name": "Prodrug",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP2D6"],
            "is_prodrug": True,
        }
        result = f.assess_candidate(candidate, rank=1)

        assert result.metabolism_risk == MetabolismRisk.CONTRAINDICATED


# ---------------------------------------------------------------------------
# Test normal metabolizer (all safe)
# ---------------------------------------------------------------------------

class TestNormalMetabolizerFiltering:
    """Tests for filtering with a normal metabolizer profile."""

    def test_all_safe_with_normal_profile(self, normal_metabolizer_profile, sample_candidates):
        """All candidates should be safe for a normal metabolizer."""
        f = PGxDrugFilter(normal_metabolizer_profile)
        summary = f.filter_candidates(sample_candidates)

        assert summary.safe_count == summary.total_candidates
        assert summary.contraindicated_count == 0
        assert summary.caution_count == 0

    def test_normal_metabolizer_no_warnings(self, normal_metabolizer_profile):
        """Normal metabolizer should produce no warnings."""
        f = PGxDrugFilter(normal_metabolizer_profile)
        candidate = {
            "name": "SafeDrug",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP2D6"],
        }
        result = f.assess_candidate(candidate, rank=1)

        assert result.warnings == []
        assert result.rank_adjustment == 0

    def test_normal_metabolizer_ranking_preserved(self, normal_metabolizer_profile, sample_candidates):
        """Ranking should be preserved for normal metabolizer (no adjustments)."""
        f = PGxDrugFilter(normal_metabolizer_profile)
        summary = f.filter_candidates(sample_candidates)

        for result in summary.results:
            assert result.original_rank == result.adjusted_rank


# ---------------------------------------------------------------------------
# Test empty PGx profile (unknown)
# ---------------------------------------------------------------------------

class TestEmptyPGxProfile:
    """Tests for filtering with no PGx data."""

    def test_empty_profile_returns_unknown(self, sample_candidates):
        """All candidates should be UNKNOWN risk with no PGx data."""
        f = PGxDrugFilter()
        summary = f.filter_candidates(sample_candidates)

        assert summary.unknown_count == summary.total_candidates
        assert summary.safe_count == 0
        assert summary.contraindicated_count == 0

    def test_empty_profile_no_rank_change(self, sample_candidates):
        """Rankings should not change without PGx data."""
        f = PGxDrugFilter()
        summary = f.filter_candidates(sample_candidates)

        for result in summary.results:
            assert result.original_rank == result.adjusted_rank

    def test_empty_profile_recommendation(self):
        """Recommendation should mention PGx data not available."""
        f = PGxDrugFilter()
        result = f.assess_candidate({"name": "Drug"}, rank=1)

        assert "not available" in result.recommendation.lower()


# ---------------------------------------------------------------------------
# Test re-ranking with mixed risk levels
# ---------------------------------------------------------------------------

class TestMixedRiskReranking:
    """Tests for re-ranking with a mix of safe, caution, and contraindicated."""

    def test_mixed_reranking_order(self, poor_metabolizer_profile):
        """Contraindicated candidates should be demoted to the end."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidates = [
            {
                "name": "BadDrug",
                "smiles": "CCO",
                "metabolism_enzymes": ["CYP2D6"],  # Contraindicated (CYP2D6 PM)
            },
            {
                "name": "GoodDrug",
                "smiles": "CCCO",
                "metabolism_enzymes": ["CYP3A5"],  # Safe (not affected)
            },
            {
                "name": "OkDrug",
                "smiles": "CCCCO",
                "metabolism_enzymes": ["CYP2C19"],  # Safe (patient is normal for CYP2C19)
            },
        ]

        summary = f.filter_candidates(candidates)

        # Safe candidates should come first in adjusted ranking
        safe_adjusted = [r.adjusted_rank for r in summary.results if r.metabolism_risk == MetabolismRisk.SAFE]
        contra_adjusted = [r.adjusted_rank for r in summary.results if r.metabolism_risk == MetabolismRisk.CONTRAINDICATED]

        if safe_adjusted and contra_adjusted:
            assert max(safe_adjusted) < min(contra_adjusted)

    def test_contraindicated_demoted_by_100(self, poor_metabolizer_profile):
        """Contraindicated candidates get a rank adjustment of +100."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidate = {
            "name": "BadDrug",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP2D6"],
        }
        result = f.assess_candidate(candidate, rank=1)

        assert result.rank_adjustment == 100

    def test_caution_demoted_by_10(self):
        """Caution candidates get a rank adjustment of +10."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Ultra Rapid Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        # Ultra rapid metabolizer + substrate = CAUTION
        candidate = {
            "name": "CautionDrug",
            "smiles": "CCO",
            "metabolism_enzymes": ["CYP2D6"],
        }
        result = f.assess_candidate(candidate, rank=1)

        assert result.metabolism_risk == MetabolismRisk.CAUTION
        assert result.rank_adjustment == 10

    def test_reranking_assigns_sequential_ranks(self, poor_metabolizer_profile, sample_candidates):
        """After re-ranking, adjusted ranks should be sequential 1..N."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        summary = f.filter_candidates(sample_candidates)

        adjusted_ranks = sorted([r.adjusted_rank for r in summary.results])
        expected = list(range(1, len(sample_candidates) + 1))
        assert adjusted_ranks == expected


# ---------------------------------------------------------------------------
# Test specific CYP2D6 poor metabolizer scenario
# ---------------------------------------------------------------------------

class TestCYP2D6PoorMetabolizer:
    """Detailed tests for the CYP2D6 poor metabolizer scenario.

    CYP2D6 is one of the most clinically important pharmacogenes. A poor
    metabolizer cannot convert prodrugs (codeine->morphine) or clear
    active drugs, leading to toxicity.
    """

    def test_cyp2d6_pm_substrate_contraindicated(self):
        """CYP2D6 PM should flag substrates as contraindicated."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Poor Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        result = f.assess_candidate(
            {"name": "MetoprololAnalog", "metabolism_enzymes": ["CYP2D6"]},
            rank=1,
        )
        assert result.metabolism_risk == MetabolismRisk.CONTRAINDICATED
        assert len(result.warnings) > 0
        assert any("avoid" in w.lower() for w in result.warnings)

    def test_cyp2d6_pm_prodrug_contraindicated(self):
        """CYP2D6 PM should flag prodrugs as contraindicated (no activation)."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Poor Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        result = f.assess_candidate(
            {"name": "CodeineAnalog", "metabolism_enzymes": ["CYP2D6"], "is_prodrug": True},
            rank=1,
        )
        assert result.metabolism_risk == MetabolismRisk.CONTRAINDICATED

    def test_cyp2d6_pm_warnings_include_gene(self):
        """Warnings should mention the gene name."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Poor Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        result = f.assess_candidate(
            {"name": "Drug", "metabolism_enzymes": ["CYP2D6"]},
            rank=1,
        )
        assert any("CYP2D6" in w for w in result.warnings)

    def test_cyp2d6_pm_affected_enzymes(self):
        """Affected enzymes list should include CYP2D6."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Poor Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        result = f.assess_candidate(
            {"name": "Drug", "metabolism_enzymes": ["CYP2D6"]},
            rank=1,
        )
        assert "CYP2D6" in result.affected_enzymes

    def test_cyp2d6_pm_original_rank_preserved(self):
        """Original rank should be preserved in the result."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Poor Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        result = f.assess_candidate(
            {"name": "Drug", "metabolism_enzymes": ["CYP2D6"]},
            rank=5,
        )
        assert result.original_rank == 5
        assert result.adjusted_rank == 105  # 5 + 100 demotion


# ---------------------------------------------------------------------------
# Test contraindicated count in summary
# ---------------------------------------------------------------------------

class TestContraindicatedCount:
    """Tests verifying correct counting of contraindicated candidates."""

    def test_contraindicated_count_matches(self, poor_metabolizer_profile):
        """Contraindicated count should match the number of flagged candidates."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        candidates = [
            {"name": "Drug1", "metabolism_enzymes": ["CYP2D6"]},  # Contraindicated
            {"name": "Drug2", "metabolism_enzymes": ["CYP2D6"]},  # Contraindicated
            {"name": "Drug3", "metabolism_enzymes": ["CYP3A5"]},  # Safe
            {"name": "Drug4", "metabolism_enzymes": ["CYP2C19"]},  # Safe (normal for 2C19)
        ]
        summary = f.filter_candidates(candidates)

        assert summary.contraindicated_count == 2
        assert summary.safe_count == 2

    def test_zero_contraindicated_with_normal_profile(self, normal_metabolizer_profile):
        """Normal metabolizer should have zero contraindicated."""
        f = PGxDrugFilter(normal_metabolizer_profile)
        candidates = [
            {"name": "Drug1", "metabolism_enzymes": ["CYP2D6"]},
            {"name": "Drug2", "metabolism_enzymes": ["CYP2C19"]},
        ]
        summary = f.filter_candidates(candidates)

        assert summary.contraindicated_count == 0

    def test_all_contraindicated_with_multiple_risks(self):
        """If all candidates hit at-risk enzymes, all should be flagged."""
        profile = {
            "gene_results": [
                {"gene": "CYP2D6", "phenotype": "Poor Metabolizer"},
                {"gene": "CYP2C19", "phenotype": "Poor Metabolizer"},
            ]
        }
        f = PGxDrugFilter(profile)
        candidates = [
            {"name": "Drug1", "metabolism_enzymes": ["CYP2D6"]},
            {"name": "Drug2", "metabolism_enzymes": ["CYP2C19"]},
            {"name": "Drug3", "metabolism_enzymes": ["CYP2D6", "CYP2C19"]},
        ]
        summary = f.filter_candidates(candidates)

        assert summary.contraindicated_count == 3

    def test_counts_sum_to_total(self, poor_metabolizer_profile, sample_candidates):
        """Risk counts should sum to total candidates."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        summary = f.filter_candidates(sample_candidates)

        total = (
            summary.safe_count
            + summary.caution_count
            + summary.contraindicated_count
            + summary.unknown_count
        )
        assert total == summary.total_candidates

    def test_empty_candidates_list(self, poor_metabolizer_profile):
        """Filtering empty candidates should return zero counts."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        summary = f.filter_candidates([])

        assert summary.total_candidates == 0
        assert summary.safe_count == 0
        assert summary.contraindicated_count == 0
        assert summary.caution_count == 0
        assert summary.unknown_count == 0
        assert summary.results == []


# ---------------------------------------------------------------------------
# Test ultra-rapid metabolizer
# ---------------------------------------------------------------------------

class TestUltraRapidMetabolizer:
    """Tests for ultra-rapid metabolizer scenarios."""

    def test_ultra_rapid_substrate_is_caution(self, ultra_rapid_profile):
        """Ultra-rapid metabolizer with substrate should be CAUTION (fast clearance)."""
        f = PGxDrugFilter(ultra_rapid_profile)
        result = f.assess_candidate(
            {"name": "Drug", "metabolism_enzymes": ["CYP2D6"]},
            rank=1,
        )
        assert result.metabolism_risk == MetabolismRisk.CAUTION

    def test_ultra_rapid_prodrug_is_contraindicated(self, ultra_rapid_profile):
        """Ultra-rapid metabolizer with prodrug should be CONTRAINDICATED (toxic activation)."""
        f = PGxDrugFilter(ultra_rapid_profile)
        result = f.assess_candidate(
            {"name": "Prodrug", "metabolism_enzymes": ["CYP2D6"], "is_prodrug": True},
            rank=1,
        )
        assert result.metabolism_risk == MetabolismRisk.CONTRAINDICATED


# ---------------------------------------------------------------------------
# Test data constants
# ---------------------------------------------------------------------------

class TestDataConstants:
    """Tests for the CYP_SUBSTRATE_MAP and PHENOTYPE_RISK constants."""

    def test_cyp_substrate_map_has_expected_genes(self):
        """CYP_SUBSTRATE_MAP should contain clinically important genes."""
        expected_genes = ["CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5", "DPYD", "TPMT", "SLCO1B1"]
        for gene in expected_genes:
            assert gene in CYP_SUBSTRATE_MAP, f"Missing gene: {gene}"

    def test_cyp_substrate_map_has_substrates(self):
        """Each gene in the map should have at least one substrate."""
        for gene, info in CYP_SUBSTRATE_MAP.items():
            assert len(info["substrates"]) > 0, f"{gene} has no substrates"

    def test_phenotype_risk_has_expected_phenotypes(self):
        """PHENOTYPE_RISK should have all standard phenotypes."""
        expected = [
            "poor_metabolizer",
            "intermediate_metabolizer",
            "normal_metabolizer",
            "rapid_metabolizer",
            "ultra_rapid_metabolizer",
        ]
        for phenotype in expected:
            assert phenotype in PHENOTYPE_RISK, f"Missing phenotype: {phenotype}"

    def test_phenotype_risk_has_both_drug_types(self):
        """Each phenotype should have risk levels for both substrates and prodrugs."""
        for phenotype, risk_map in PHENOTYPE_RISK.items():
            assert "substrate" in risk_map, f"{phenotype} missing substrate risk"
            assert "prodrug" in risk_map, f"{phenotype} missing prodrug risk"


# ---------------------------------------------------------------------------
# Test PGxFilterResult dataclass
# ---------------------------------------------------------------------------

class TestPGxFilterResult:
    """Tests for the PGxFilterResult dataclass."""

    def test_default_values(self):
        """Default values should be sensible."""
        result = PGxFilterResult(
            candidate_id="test",
            original_rank=1,
            adjusted_rank=1,
            metabolism_risk=MetabolismRisk.UNKNOWN,
        )
        assert result.affected_enzymes == []
        assert result.warnings == []
        assert result.recommendation == ""
        assert result.rank_adjustment == 0

    def test_result_fields_populated(self, poor_metabolizer_profile):
        """All fields should be populated after assessment."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        result = f.assess_candidate(
            {"name": "TestDrug", "metabolism_enzymes": ["CYP2D6"]},
            rank=3,
        )
        assert result.candidate_id == "TestDrug"
        assert result.original_rank == 3
        assert result.recommendation != ""


# ---------------------------------------------------------------------------
# Test PGxFilterSummary dataclass
# ---------------------------------------------------------------------------

class TestPGxFilterSummary:
    """Tests for the PGxFilterSummary dataclass."""

    def test_summary_contains_profile(self, poor_metabolizer_profile):
        """Summary should include the patient PGx profile."""
        f = PGxDrugFilter(poor_metabolizer_profile)
        summary = f.filter_candidates([{"name": "Drug"}])

        assert "CYP2D6" in summary.patient_pgx_profile
        assert summary.patient_pgx_profile["CYP2D6"] == "poor_metabolizer"

    def test_summary_total_candidates(self, normal_metabolizer_profile, sample_candidates):
        """Summary total_candidates should match input length."""
        f = PGxDrugFilter(normal_metabolizer_profile)
        summary = f.filter_candidates(sample_candidates)

        assert summary.total_candidates == len(sample_candidates)

    def test_summary_results_length(self, normal_metabolizer_profile, sample_candidates):
        """Summary results should have one entry per candidate."""
        f = PGxDrugFilter(normal_metabolizer_profile)
        summary = f.filter_candidates(sample_candidates)

        assert len(summary.results) == len(sample_candidates)
