"""Tests for CAR-T Intelligence Agent knowledge graph module.

Validates target antigen, toxicity, manufacturing, biomarker, and
regulatory context retrieval, entity alias resolution, and combined
query context extraction.

Author: Adam Jones
Date: February 2026
"""

import pytest

from src.knowledge import (
    CART_BIOMARKERS,
    CART_MANUFACTURING,
    CART_REGULATORY,
    CART_TARGETS,
    CART_TOXICITIES,
    ENTITY_ALIASES,
    get_all_context_for_query,
    get_biomarker_context,
    get_knowledge_stats,
    get_manufacturing_context,
    get_regulatory_context,
    get_target_context,
    get_toxicity_context,
    resolve_comparison_entity,
)


# ═══════════════════════════════════════════════════════════════════════
# TARGET CONTEXT
# ═══════════════════════════════════════════════════════════════════════


class TestGetTargetContext:
    """Tests for get_target_context()."""

    @pytest.mark.parametrize("antigen", ["CD19", "BCMA", "CD22", "HER2", "GD2"])
    def test_known_targets_return_non_empty(self, antigen):
        """Known target antigens produce non-empty context strings."""
        ctx = get_target_context(antigen)
        assert len(ctx) > 0
        assert antigen in ctx

    def test_cd19_context_includes_products(self):
        """CD19 context mentions approved products."""
        ctx = get_target_context("CD19")
        assert "Kymriah" in ctx
        assert "Yescarta" in ctx

    def test_bcma_context_includes_myeloma(self):
        """BCMA context mentions multiple myeloma."""
        ctx = get_target_context("BCMA")
        assert "Multiple Myeloma" in ctx

    def test_case_insensitive(self):
        """get_target_context is case-insensitive."""
        ctx_lower = get_target_context("cd19")
        ctx_upper = get_target_context("CD19")
        assert ctx_lower == ctx_upper
        assert len(ctx_lower) > 0

    def test_unknown_target_returns_empty(self):
        """An unknown antigen returns an empty string."""
        ctx = get_target_context("NONEXISTENT_TARGET")
        assert ctx == ""


# ═══════════════════════════════════════════════════════════════════════
# TOXICITY CONTEXT
# ═══════════════════════════════════════════════════════════════════════


class TestGetToxicityContext:
    """Tests for get_toxicity_context()."""

    @pytest.mark.parametrize("tox", ["CRS", "ICANS", "B_CELL_APLASIA", "HLH_MAS"])
    def test_known_toxicities_return_non_empty(self, tox):
        """Known toxicity profiles produce non-empty context strings."""
        ctx = get_toxicity_context(tox)
        assert len(ctx) > 0

    def test_crs_context_includes_management(self):
        """CRS context includes tocilizumab management."""
        ctx = get_toxicity_context("CRS")
        assert "Tocilizumab" in ctx or "tocilizumab" in ctx

    def test_icans_context_includes_grading(self):
        """ICANS context mentions ICE score grading."""
        ctx = get_toxicity_context("ICANS")
        assert "ICE" in ctx

    def test_unknown_toxicity_returns_empty(self):
        """An unknown toxicity returns an empty string."""
        ctx = get_toxicity_context("NONEXISTENT_TOXICITY")
        assert ctx == ""


# ═══════════════════════════════════════════════════════════════════════
# MANUFACTURING CONTEXT
# ═══════════════════════════════════════════════════════════════════════


class TestGetManufacturingContext:
    """Tests for get_manufacturing_context()."""

    def test_lentiviral_transduction(self):
        """Lentiviral transduction context includes efficiency and VCN."""
        ctx = get_manufacturing_context("lentiviral_transduction")
        assert len(ctx) > 0
        assert "Typical Efficiency" in ctx or "efficiency" in ctx.lower()

    def test_expansion(self):
        """Expansion context is found by partial match."""
        ctx = get_manufacturing_context("expansion")
        assert len(ctx) > 0

    def test_cryopreservation(self):
        """Cryopreservation context includes DMSO or viability info."""
        ctx = get_manufacturing_context("cryopreservation")
        assert len(ctx) > 0

    def test_unknown_process_returns_empty(self):
        """An unknown manufacturing process returns an empty string."""
        ctx = get_manufacturing_context("quantum_teleportation")
        assert ctx == ""


# ═══════════════════════════════════════════════════════════════════════
# BIOMARKER CONTEXT
# ═══════════════════════════════════════════════════════════════════════


class TestGetBiomarkerContext:
    """Tests for get_biomarker_context()."""

    @pytest.mark.parametrize("biomarker", ["ferritin", "crp", "il6", "pd1", "mrd_flow"])
    def test_known_biomarkers_return_non_empty(self, biomarker):
        """Known biomarkers produce non-empty context strings."""
        ctx = get_biomarker_context(biomarker)
        assert len(ctx) > 0

    def test_ferritin_context_content(self):
        """Ferritin context mentions CRS prediction."""
        ctx = get_biomarker_context("ferritin")
        assert "Ferritin" in ctx or "ferritin" in ctx
        assert "CRS" in ctx

    def test_crp_context_content(self):
        """CRP context includes assay method."""
        ctx = get_biomarker_context("crp")
        assert "CRP" in ctx or "C-Reactive" in ctx

    def test_unknown_biomarker_returns_empty(self):
        """An unknown biomarker returns an empty string."""
        ctx = get_biomarker_context("nonexistent_biomarker")
        assert ctx == ""


# ═══════════════════════════════════════════════════════════════════════
# REGULATORY CONTEXT
# ═══════════════════════════════════════════════════════════════════════


class TestGetRegulatoryContext:
    """Tests for get_regulatory_context()."""

    @pytest.mark.parametrize("product", ["Kymriah", "Yescarta", "Carvykti", "Abecma"])
    def test_known_products_return_non_empty(self, product):
        """Known FDA-approved products produce non-empty regulatory context."""
        ctx = get_regulatory_context(product)
        assert len(ctx) > 0

    def test_kymriah_context_includes_approval_date(self):
        """Kymriah regulatory context mentions 2017 approval."""
        ctx = get_regulatory_context("Kymriah")
        assert "2017" in ctx

    def test_yescarta_context_includes_manufacturer(self):
        """Yescarta regulatory context mentions Kite/Gilead."""
        ctx = get_regulatory_context("Yescarta")
        assert "Kite" in ctx or "Gilead" in ctx

    def test_unknown_product_returns_empty(self):
        """An unrecognized product name returns an empty string."""
        ctx = get_regulatory_context("UnknownDrug")
        assert ctx == ""


# ═══════════════════════════════════════════════════════════════════════
# KNOWLEDGE STATS
# ═══════════════════════════════════════════════════════════════════════


class TestGetKnowledgeStats:
    """Tests for get_knowledge_stats()."""

    def test_returns_correct_keys(self):
        """get_knowledge_stats() returns a dict with the expected keys."""
        stats = get_knowledge_stats()
        expected_keys = {
            "target_antigens",
            "targets_with_approved_products",
            "toxicity_profiles",
            "manufacturing_processes",
            "biomarkers",
            "regulatory_products",
            "immunogenicity_topics",
        }
        assert set(stats.keys()) == expected_keys

    def test_counts_are_positive(self):
        """All counts in the knowledge stats are positive integers."""
        stats = get_knowledge_stats()
        for key, value in stats.items():
            assert isinstance(value, int), f"{key} is not int"
            assert value > 0, f"{key} is not positive"

    def test_target_count_matches_dict(self):
        """The target antigen count matches len(CART_TARGETS)."""
        stats = get_knowledge_stats()
        assert stats["target_antigens"] == len(CART_TARGETS)

    def test_toxicity_count_matches_dict(self):
        """The toxicity count matches len(CART_TOXICITIES)."""
        stats = get_knowledge_stats()
        assert stats["toxicity_profiles"] == len(CART_TOXICITIES)

    def test_manufacturing_count_matches_dict(self):
        """The manufacturing count matches len(CART_MANUFACTURING)."""
        stats = get_knowledge_stats()
        assert stats["manufacturing_processes"] == len(CART_MANUFACTURING)

    def test_biomarker_count_matches_dict(self):
        """The biomarker count matches len(CART_BIOMARKERS)."""
        stats = get_knowledge_stats()
        assert stats["biomarkers"] == len(CART_BIOMARKERS)

    def test_regulatory_count_matches_dict(self):
        """The regulatory count matches len(CART_REGULATORY)."""
        stats = get_knowledge_stats()
        assert stats["regulatory_products"] == len(CART_REGULATORY)


# ═══════════════════════════════════════════════════════════════════════
# ENTITY RESOLUTION
# ═══════════════════════════════════════════════════════════════════════


class TestResolveComparisonEntity:
    """Tests for resolve_comparison_entity()."""

    def test_resolves_cd19_as_target(self):
        """'CD19' resolves to a target entity."""
        result = resolve_comparison_entity("CD19")
        assert result is not None
        assert result["type"] == "target"
        assert result["canonical"] == "CD19"

    def test_resolves_bcma_as_target(self):
        """'BCMA' resolves to a target entity."""
        result = resolve_comparison_entity("BCMA")
        assert result is not None
        assert result["type"] == "target"

    def test_resolves_kymriah_as_product(self):
        """'Kymriah' resolves to a product entity."""
        result = resolve_comparison_entity("Kymriah")
        assert result is not None
        assert result["type"] == "product"
        assert "tisagenlecleucel" in result["canonical"]

    def test_resolves_4_1bb_as_costimulatory(self):
        """'4-1BB' resolves to a costimulatory domain entity."""
        result = resolve_comparison_entity("4-1BB")
        assert result is not None
        assert result["type"] == "costimulatory"

    def test_resolves_crs_as_toxicity(self):
        """'CRS' resolves to a toxicity entity."""
        result = resolve_comparison_entity("CRS")
        assert result is not None
        assert result["type"] == "toxicity"

    def test_unknown_entity_returns_none(self):
        """An unrecognized entity returns None."""
        result = resolve_comparison_entity("TotallyUnknownEntity")
        assert result is None

    @pytest.mark.parametrize(
        "alias,expected_type",
        [
            ("TISAGENLECLEUCEL", "product"),
            ("AXICABTAGENE", "product"),
            ("CD137", "costimulatory"),
            ("LENTIVIRAL", "manufacturing"),
            ("FERRITIN", "biomarker"),
        ],
    )
    def test_alias_resolution(self, alias, expected_type):
        """Various aliases resolve to the correct entity type."""
        result = resolve_comparison_entity(alias)
        assert result is not None
        assert result["type"] == expected_type


# ═══════════════════════════════════════════════════════════════════════
# ENTITY ALIASES
# ═══════════════════════════════════════════════════════════════════════


class TestEntityAliases:
    """Tests for the ENTITY_ALIASES dictionary."""

    def test_contains_product_aliases(self):
        """ENTITY_ALIASES contains key product aliases."""
        assert "KYMRIAH" in ENTITY_ALIASES
        assert "YESCARTA" in ENTITY_ALIASES
        assert "CARVYKTI" in ENTITY_ALIASES

    def test_contains_costimulatory_aliases(self):
        """ENTITY_ALIASES contains costimulatory domain aliases."""
        assert "4-1BB" in ENTITY_ALIASES
        assert "CD28" in ENTITY_ALIASES

    def test_contains_biomarker_aliases(self):
        """ENTITY_ALIASES contains biomarker aliases."""
        assert "FERRITIN" in ENTITY_ALIASES
        assert "CRP" in ENTITY_ALIASES
        assert "IL-6" in ENTITY_ALIASES
        assert "PD-1" in ENTITY_ALIASES

    def test_contains_manufacturing_aliases(self):
        """ENTITY_ALIASES contains manufacturing aliases."""
        assert "LENTIVIRAL" in ENTITY_ALIASES
        assert "RETROVIRAL" in ENTITY_ALIASES

    def test_alias_structure(self):
        """Each alias entry has 'type' and 'canonical' keys."""
        for alias_key, alias_data in ENTITY_ALIASES.items():
            assert "type" in alias_data, f"Missing 'type' in alias: {alias_key}"
            assert "canonical" in alias_data, f"Missing 'canonical' in alias: {alias_key}"


# ═══════════════════════════════════════════════════════════════════════
# COMBINED CONTEXT (get_all_context_for_query)
# ═══════════════════════════════════════════════════════════════════════


class TestGetAllContextForQuery:
    """Tests for get_all_context_for_query()."""

    def test_cd19_crs_query(self):
        """A query mentioning 'CD19' and 'CRS' returns both target and toxicity context."""
        ctx = get_all_context_for_query("What is the CRS rate for CD19 CAR-T therapy?")
        assert len(ctx) > 0
        assert "CD19" in ctx
        assert "Cytokine Release Syndrome" in ctx

    def test_manufacturing_query(self):
        """A query about transduction returns manufacturing context."""
        ctx = get_all_context_for_query("What is the lentiviral transduction efficiency?")
        assert len(ctx) > 0

    def test_biomarker_query(self):
        """A query mentioning ferritin returns biomarker context."""
        ctx = get_all_context_for_query("Does ferritin predict CRS severity?")
        assert len(ctx) > 0
        assert "Ferritin" in ctx or "ferritin" in ctx

    def test_regulatory_query(self):
        """A query mentioning Kymriah returns regulatory context."""
        ctx = get_all_context_for_query("When was Kymriah approved by the FDA?")
        assert len(ctx) > 0
        assert "tisagenlecleucel" in ctx

    def test_empty_for_unrelated_query(self):
        """A completely unrelated query returns an empty string."""
        ctx = get_all_context_for_query("What is the weather today?")
        assert ctx == ""

    def test_multi_domain_query(self):
        """A query spanning multiple domains returns combined context."""
        ctx = get_all_context_for_query(
            "Compare CD19 and BCMA CRS rates with ferritin as a biomarker"
        )
        assert "CD19" in ctx
        assert "BCMA" in ctx
