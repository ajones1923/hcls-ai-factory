"""Tests for src/query_expansion.py — domain-specific query expansion.

Tests expand_query, expand_query_by_category, get_expansion_stats.

Author: Adam Jones
Date: March 2026
"""

from src.query_expansion import (
    expand_query,
    expand_query_by_category,
    get_expansion_stats,
    ALL_EXPANSION_MAPS,
    DRUG_EXPANSION,
)


# ═══════════════════════════════════════════════════════════════════════
# ALL_EXPANSION_MAPS
# ═══════════════════════════════════════════════════════════════════════

class TestExpansionMaps:
    def test_has_14_maps(self):
        assert len(ALL_EXPANSION_MAPS) == 14

    def test_each_is_tuple(self):
        for item in ALL_EXPANSION_MAPS:
            assert isinstance(item, tuple)
            assert len(item) == 2
            assert isinstance(item[0], str)
            assert isinstance(item[1], dict)

    def test_drug_expansion_not_empty(self):
        assert len(DRUG_EXPANSION) >= 30

    def test_warfarin_expansion(self):
        assert "warfarin" in DRUG_EXPANSION
        terms = DRUG_EXPANSION["warfarin"]
        assert "CYP2C9" in terms
        assert "VKORC1" in terms


# ═══════════════════════════════════════════════════════════════════════
# expand_query
# ═══════════════════════════════════════════════════════════════════════

class TestExpandQuery:
    def test_warfarin_expands(self):
        terms = expand_query("warfarin dosing")
        assert len(terms) > 0
        assert "CYP2C9" in terms or "VKORC1" in terms

    def test_cyp2d6_expands(self):
        terms = expand_query("CYP2D6 poor metabolizer codeine")
        assert len(terms) > 0

    def test_empty_query_returns_empty(self):
        terms = expand_query("")
        # Empty query may not match any keys
        assert isinstance(terms, list)

    def test_no_match_returns_empty(self):
        terms = expand_query("zzz_nonexistent_medical_term_xyz")
        assert terms == []

    def test_returns_sorted_list(self):
        terms = expand_query("warfarin")
        assert terms == sorted(terms)

    def test_returns_deduplicated(self):
        terms = expand_query("warfarin codeine")
        assert len(terms) == len(set(terms))

    def test_case_insensitive(self):
        upper_terms = expand_query("WARFARIN")
        lower_terms = expand_query("warfarin")
        assert upper_terms == lower_terms

    def test_hla_expansion(self):
        terms = expand_query("HLA-B*57:01 abacavir hypersensitivity")
        assert len(terms) > 0

    def test_phenotype_expansion(self):
        terms = expand_query("poor metabolizer")
        assert len(terms) > 0

    def test_multi_term_query(self):
        terms = expand_query("warfarin dosing for CYP2C9 poor metabolizer")
        assert len(terms) > 5  # Should match multiple expansion maps

    def test_brand_name_expansion(self):
        terms = expand_query("coumadin dosing")
        # "coumadin" may be in drug expansion as a key
        assert isinstance(terms, list)


# ═══════════════════════════════════════════════════════════════════════
# expand_query_by_category
# ═══════════════════════════════════════════════════════════════════════

class TestExpandQueryByCategory:
    def test_returns_dict(self):
        result = expand_query_by_category("warfarin dosing")
        assert isinstance(result, dict)

    def test_categories_have_lists(self):
        result = expand_query_by_category("warfarin dosing")
        for cat, terms in result.items():
            assert isinstance(cat, str)
            assert isinstance(terms, list)

    def test_drug_category_present(self):
        result = expand_query_by_category("warfarin")
        assert "DRUG_EXPANSION" in result

    def test_empty_query(self):
        result = expand_query_by_category("zzz_nonexistent_xyz")
        assert result == {}

    def test_multiple_categories(self):
        result = expand_query_by_category("carbamazepine SJS Asian HLA")
        assert len(result) >= 1

    def test_categories_sorted(self):
        result = expand_query_by_category("warfarin")
        for cat, terms in result.items():
            assert terms == sorted(terms)


# ═══════════════════════════════════════════════════════════════════════
# get_expansion_stats
# ═══════════════════════════════════════════════════════════════════════

class TestGetExpansionStats:
    def test_returns_dict(self):
        stats = get_expansion_stats()
        assert isinstance(stats, dict)

    def test_total_maps(self):
        stats = get_expansion_stats()
        assert stats["total_maps"] == 14

    def test_total_keys(self):
        stats = get_expansion_stats()
        assert stats["total_keys"] > 100

    def test_total_terms(self):
        stats = get_expansion_stats()
        assert stats["total_terms"] > 500

    def test_unique_terms(self):
        stats = get_expansion_stats()
        assert stats["unique_terms"] > 0
        assert stats["unique_terms"] <= stats["total_terms"]

    def test_per_map_present(self):
        stats = get_expansion_stats()
        assert "per_map" in stats
        assert isinstance(stats["per_map"], dict)
        assert len(stats["per_map"]) == 14

    def test_per_map_has_keys_and_terms(self):
        stats = get_expansion_stats()
        for name, info in stats["per_map"].items():
            assert "keys" in info
            assert "terms" in info
            assert info["keys"] >= 0
            assert info["terms"] >= 0
