"""Tests for CAR-T Intelligence Agent query expansion module.

Validates keyword-based query expansion across all 11 expansion map
categories, categorized expansion, and expansion statistics.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.query_expansion import (
    ALL_EXPANSION_MAPS,
    BIOMARKER_EXPANSION,
    CONSTRUCT_EXPANSION,
    DISEASE_EXPANSION,
    MANUFACTURING_EXPANSION,
    MECHANISM_EXPANSION,
    REALWORLD_EXPANSION,
    REGULATORY_EXPANSION,
    SAFETY_EXPANSION,
    SEQUENCE_EXPANSION,
    TARGET_ANTIGEN_EXPANSION,
    TOXICITY_EXPANSION,
    expand_query,
    expand_query_by_category,
    get_expansion_stats,
)


# ═══════════════════════════════════════════════════════════════════════
# expand_query()
# ═══════════════════════════════════════════════════════════════════════


class TestExpandQuery:
    """Tests for the expand_query() function."""

    def test_cd19_cart_returns_non_empty(self):
        """A query about 'CD19 CAR-T' produces expansion terms."""
        terms = expand_query("What is the efficacy of cd19 CAR-T therapy?")
        assert len(terms) > 0
        # Should include target antigen expansion terms
        assert any("CD19" in t for t in terms)

    def test_crs_query_returns_toxicity_terms(self):
        """A query about CRS produces toxicity-related expansion terms."""
        terms = expand_query("How to manage crs after CAR-T infusion?")
        assert len(terms) > 0
        assert any("tocilizumab" in t.lower() for t in terms)

    def test_manufacturing_query(self):
        """A query about transduction produces manufacturing terms."""
        terms = expand_query("What affects lentiviral transduction efficiency?")
        assert len(terms) > 0
        assert any("VCN" in t or "MOI" in t for t in terms)

    def test_unrelated_query_returns_empty_or_minimal(self):
        """A completely unrelated query returns an empty or very small list."""
        terms = expand_query("quantum computing in semiconductor design")
        assert len(terms) <= 5  # Should be empty or very small

    def test_expansion_is_deduplicated(self):
        """Expansion results contain no duplicate terms."""
        terms = expand_query("cd19 resistance antigen escape mechanisms")
        assert len(terms) == len(set(terms))

    def test_expansion_is_sorted(self):
        """Expansion results are returned in sorted order."""
        terms = expand_query("bcma multiple myeloma crs")
        assert terms == sorted(terms)

    @pytest.mark.parametrize(
        "query_keyword",
        [
            "cd19",
            "bcma",
            "crs",
            "icans",
            "transduction",
            "expansion",
            "exhaustion",
            "resistance",
            "scfv",
            "fda",
            "real-world",
        ],
    )
    def test_various_keywords_produce_terms(self, query_keyword):
        """Queries containing domain keywords produce at least one expansion term."""
        terms = expand_query(f"Tell me about {query_keyword}")
        assert len(terms) > 0, f"No expansion terms for keyword: {query_keyword}"


# ═══════════════════════════════════════════════════════════════════════
# expand_query_by_category()
# ═══════════════════════════════════════════════════════════════════════


class TestExpandQueryByCategory:
    """Tests for the expand_query_by_category() function."""

    def test_returns_categorized_dict(self):
        """expand_query_by_category() returns a dict with category keys."""
        result = expand_query_by_category("cd19 transduction efficiency")
        assert isinstance(result, dict)
        assert len(result) > 0

    def test_category_keys_match_all_expansion_maps(self):
        """Category keys are a subset of ALL_EXPANSION_MAPS names."""
        valid_categories = {name for name, _ in ALL_EXPANSION_MAPS}
        result = expand_query_by_category("cd19 crs bcma resistance exhaustion")
        for cat in result.keys():
            assert cat in valid_categories, f"Unknown category: {cat}"

    def test_category_values_are_sorted_lists(self):
        """Each category value is a sorted list of strings."""
        result = expand_query_by_category("cd19 crs transduction")
        for cat, terms in result.items():
            assert isinstance(terms, list), f"{cat} value is not a list"
            assert terms == sorted(terms), f"{cat} terms are not sorted"

    def test_empty_for_unrelated_query(self):
        """An unrelated query returns an empty dict."""
        result = expand_query_by_category("quantum computing")
        assert len(result) == 0

    def test_multi_category_query(self):
        """A multi-domain query populates multiple categories."""
        result = expand_query_by_category(
            "cd19 crs resistance lentiviral scfv fda real-world ferritin"
        )
        assert len(result) >= 3, "Expected at least 3 categories matched"


# ═══════════════════════════════════════════════════════════════════════
# get_expansion_stats()
# ═══════════════════════════════════════════════════════════════════════


class TestGetExpansionStats:
    """Tests for the get_expansion_stats() function."""

    def test_returns_correct_number_of_categories(self):
        """get_expansion_stats() returns exactly 12 expansion categories."""
        stats = get_expansion_stats()
        assert len(stats) == 12

    def test_each_category_has_keywords_and_total_terms(self):
        """Each category in stats has 'keywords' and 'total_terms' counts."""
        stats = get_expansion_stats()
        for cat, data in stats.items():
            assert "keywords" in data, f"{cat} missing 'keywords'"
            assert "total_terms" in data, f"{cat} missing 'total_terms'"
            assert data["keywords"] > 0, f"{cat} has 0 keywords"
            assert data["total_terms"] > 0, f"{cat} has 0 total_terms"

    def test_category_names_match_all_expansion_maps(self):
        """Stats category names match ALL_EXPANSION_MAPS names."""
        stats = get_expansion_stats()
        expected_names = {name for name, _ in ALL_EXPANSION_MAPS}
        assert set(stats.keys()) == expected_names


# ═══════════════════════════════════════════════════════════════════════
# ALL_EXPANSION_MAPS structure
# ═══════════════════════════════════════════════════════════════════════


class TestAllExpansionMaps:
    """Tests for the ALL_EXPANSION_MAPS registry."""

    def test_has_12_entries(self):
        """ALL_EXPANSION_MAPS contains exactly 12 expansion map tuples."""
        assert len(ALL_EXPANSION_MAPS) == 12

    def test_expected_category_names(self):
        """ALL_EXPANSION_MAPS contains all expected category names."""
        names = [name for name, _ in ALL_EXPANSION_MAPS]
        expected = [
            "Target Antigen",
            "Disease",
            "Toxicity",
            "Manufacturing",
            "Mechanism",
            "Construct",
            "Safety",
            "Biomarker",
            "Regulatory",
            "Sequence",
            "RealWorld",
            "Immunogenicity",
        ]
        assert names == expected

    def test_each_entry_is_name_dict_tuple(self):
        """Each entry in ALL_EXPANSION_MAPS is a (str, dict) tuple."""
        for entry in ALL_EXPANSION_MAPS:
            assert isinstance(entry, tuple)
            assert len(entry) == 2
            name, mapping = entry
            assert isinstance(name, str)
            assert isinstance(mapping, dict)

    def test_individual_maps_are_non_empty(self):
        """Each individual expansion map dictionary is non-empty."""
        assert len(TARGET_ANTIGEN_EXPANSION) > 0
        assert len(DISEASE_EXPANSION) > 0
        assert len(TOXICITY_EXPANSION) > 0
        assert len(MANUFACTURING_EXPANSION) > 0
        assert len(MECHANISM_EXPANSION) > 0
        assert len(CONSTRUCT_EXPANSION) > 0
        assert len(SAFETY_EXPANSION) > 0
        assert len(BIOMARKER_EXPANSION) > 0
        assert len(REGULATORY_EXPANSION) > 0
        assert len(SEQUENCE_EXPANSION) > 0
        assert len(REALWORLD_EXPANSION) > 0
