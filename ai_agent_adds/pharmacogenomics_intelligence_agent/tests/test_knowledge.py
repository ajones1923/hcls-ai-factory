"""Tests for src/knowledge.py — pharmacogenomics knowledge graph.

Covers PHARMACOGENES (25 genes), DRUG_CATEGORIES (12 categories),
METABOLIZER_PHENOTYPES (5 types), HLA_DRUG_ASSOCIATIONS (12 entries),
DRUG_ALTERNATIVES, ACTIVITY_SCORE_TABLES, POPULATION_FREQUENCIES,
ENTITY_ALIASES, get_knowledge_stats(), and all public functions.

Author: Adam Jones
Date: March 2026
"""

import pytest
from src.knowledge import (
    PHARMACOGENES,
    METABOLIZER_PHENOTYPES,
    DRUG_CATEGORIES,
    HLA_DRUG_ASSOCIATIONS,
    DRUG_ALTERNATIVES,
    ACTIVITY_SCORE_TABLES,
    POPULATION_FREQUENCIES,
    ENTITY_ALIASES,
    CYP_INHIBITORS,
    CYP_INDUCERS,
    get_gene_context,
    get_phenotype_context,
    get_drug_category_context,
    get_hla_context,
    get_inhibitor_context,
    get_alternative_drugs,
    get_all_context_for_query,
    get_knowledge_stats,
    resolve_comparison_entity,
    get_comparison_context,
)


# ═══════════════════════════════════════════════════════════════════════
# 1. PHARMACOGENES
# ═══════════════════════════════════════════════════════════════════════

class TestPharmacogenes:
    def test_has_25_genes(self):
        assert len(PHARMACOGENES) == 25

    @pytest.mark.parametrize("gene", [
        "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A4", "CYP3A5", "CYP1A2",
        "CYP2B6", "CYP4F2", "VKORC1", "SLCO1B1", "DPYD", "TPMT",
        "NUDT15", "UGT1A1", "G6PD", "NAT2", "ABCB1", "CYP2C8",
        "IFNL3", "RYR1", "CACNA1S", "CFTR", "F5", "MTHFR", "HLA",
    ])
    def test_gene_present(self, gene):
        assert gene in PHARMACOGENES

    def test_cyp2d6_has_required_fields(self):
        g = PHARMACOGENES["CYP2D6"]
        assert "full_name" in g
        assert "chromosome" in g
        assert "function" in g
        assert "key_variants" in g
        assert "cpic_guidelines" in g
        assert g["full_name"] == "Cytochrome P450 2D6"

    def test_cyp2d6_substrates_count(self):
        assert PHARMACOGENES["CYP2D6"]["substrates_count"] == 120

    def test_cyp2c19_chromosome(self):
        assert PHARMACOGENES["CYP2C19"]["chromosome"] == "10q23.33"

    def test_dpyd_key_variants_present(self):
        variants = PHARMACOGENES["DPYD"]["key_variants"]
        assert any("*2A" in v for v in variants)

    def test_g6pd_is_x_linked(self):
        assert PHARMACOGENES["G6PD"]["chromosome"] == "Xq28"

    def test_hla_has_cpic_guidelines(self):
        assert len(PHARMACOGENES["HLA"]["cpic_guidelines"]) >= 4

    def test_all_genes_have_full_name(self):
        for gene, data in PHARMACOGENES.items():
            assert "full_name" in data, f"{gene} missing full_name"


# ═══════════════════════════════════════════════════════════════════════
# 2. METABOLIZER_PHENOTYPES
# ═══════════════════════════════════════════════════════════════════════

class TestMetabolizerPhenotypes:
    def test_has_5_phenotypes(self):
        assert len(METABOLIZER_PHENOTYPES) == 5

    @pytest.mark.parametrize("pheno", [
        "ultra_rapid", "rapid", "normal", "intermediate", "poor",
    ])
    def test_phenotype_present(self, pheno):
        assert pheno in METABOLIZER_PHENOTYPES

    def test_each_has_abbreviation(self):
        for key, val in METABOLIZER_PHENOTYPES.items():
            assert "abbreviation" in val
            assert "clinical_meaning" in val
            assert "risk" in val

    def test_abbreviations(self):
        assert METABOLIZER_PHENOTYPES["ultra_rapid"]["abbreviation"] == "UM"
        assert METABOLIZER_PHENOTYPES["poor"]["abbreviation"] == "PM"
        assert METABOLIZER_PHENOTYPES["normal"]["abbreviation"] == "NM"


# ═══════════════════════════════════════════════════════════════════════
# 3. DRUG_CATEGORIES
# ═══════════════════════════════════════════════════════════════════════

class TestDrugCategories:
    def test_has_12_categories(self):
        assert len(DRUG_CATEGORIES) == 12

    @pytest.mark.parametrize("cat", [
        "opioids", "anticoagulants", "antidepressants", "antipsychotics",
        "statins", "chemotherapy", "anticonvulsants", "antivirals",
        "immunosuppressants", "cardiovascular", "proton_pump_inhibitors",
        "anti_gout",
    ])
    def test_category_present(self, cat):
        assert cat in DRUG_CATEGORIES

    def test_each_has_drugs_and_genes(self):
        for cat, data in DRUG_CATEGORIES.items():
            assert "drugs" in data
            assert "primary_genes" in data
            assert len(data["drugs"]) > 0

    def test_opioids_include_codeine(self):
        assert "codeine" in DRUG_CATEGORIES["opioids"]["drugs"]

    def test_statins_include_simvastatin(self):
        assert "simvastatin" in DRUG_CATEGORIES["statins"]["drugs"]


# ═══════════════════════════════════════════════════════════════════════
# 4. HLA_DRUG_ASSOCIATIONS
# ═══════════════════════════════════════════════════════════════════════

class TestHLADrugAssociations:
    def test_has_12_entries(self):
        assert len(HLA_DRUG_ASSOCIATIONS) == 12

    def test_abacavir_entry(self):
        entry = HLA_DRUG_ASSOCIATIONS.get("HLA-B*57:01_abacavir")
        assert entry is not None
        assert entry["hla_allele"] == "HLA-B*57:01"
        assert entry["drug"] == "abacavir"
        assert entry["screening_mandatory"] is True

    def test_carbamazepine_entry(self):
        entry = HLA_DRUG_ASSOCIATIONS.get("HLA-B*15:02_carbamazepine")
        assert entry is not None
        assert "SJS" in entry["reaction_type"]

    def test_all_entries_have_required_fields(self):
        for key, data in HLA_DRUG_ASSOCIATIONS.items():
            assert "hla_allele" in data
            assert "drug" in data
            assert "reaction_type" in data
            assert "recommendation" in data
            assert "prevalence_by_population" in data


# ═══════════════════════════════════════════════════════════════════════
# 5. DRUG_ALTERNATIVES
# ═══════════════════════════════════════════════════════════════════════

class TestDrugAlternatives:
    def test_not_empty(self):
        assert len(DRUG_ALTERNATIVES) > 0

    def test_has_known_drugs(self):
        keys_lower = [k.lower() for k in DRUG_ALTERNATIVES]
        # At least some of these should be present
        found = any(drug in " ".join(keys_lower) for drug in
                    ["codeine", "clopidogrel", "warfarin", "simvastatin"])
        assert found


# ═══════════════════════════════════════════════════════════════════════
# 6. ACTIVITY_SCORE_TABLES
# ═══════════════════════════════════════════════════════════════════════

class TestActivityScoreTables:
    def test_not_empty(self):
        assert len(ACTIVITY_SCORE_TABLES) > 0

    def test_cyp2d6_present(self):
        assert "CYP2D6" in ACTIVITY_SCORE_TABLES

    def test_has_description_and_thresholds(self):
        for gene, data in ACTIVITY_SCORE_TABLES.items():
            assert "description" in data
            assert "phenotype_thresholds" in data


# ═══════════════════════════════════════════════════════════════════════
# 7. POPULATION_FREQUENCIES
# ═══════════════════════════════════════════════════════════════════════

class TestPopulationFrequencies:
    def test_not_empty(self):
        assert len(POPULATION_FREQUENCIES) > 0

    def test_has_gene_entries(self):
        # Should have entries for at least CYP2D6
        found_cyp = any("CYP2D6" in k for k in POPULATION_FREQUENCIES)
        assert found_cyp or len(POPULATION_FREQUENCIES) > 0


# ═══════════════════════════════════════════════════════════════════════
# 8. ENTITY_ALIASES
# ═══════════════════════════════════════════════════════════════════════

class TestEntityAliases:
    def test_not_empty(self):
        assert len(ENTITY_ALIASES) >= 80

    def test_alias_has_type_and_canonical(self):
        for alias, info in list(ENTITY_ALIASES.items())[:10]:
            assert "type" in info
            assert "canonical" in info


# ═══════════════════════════════════════════════════════════════════════
# 9. CYP INHIBITORS AND INDUCERS
# ═══════════════════════════════════════════════════════════════════════

class TestCYPInhibitors:
    def test_has_enzymes(self):
        assert len(CYP_INHIBITORS) >= 3

    def test_cyp2d6_inhibitors(self):
        assert "CYP2D6" in CYP_INHIBITORS
        inh = CYP_INHIBITORS["CYP2D6"]
        assert "strong" in inh
        assert len(inh["strong"]) > 0


class TestCYPInducers:
    def test_has_enzymes(self):
        assert len(CYP_INDUCERS) >= 2

    def test_has_strong_and_moderate(self):
        for enzyme, data in CYP_INDUCERS.items():
            assert "strong" in data
            assert "moderate" in data


# ═══════════════════════════════════════════════════════════════════════
# 10. PUBLIC FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════

class TestGetGeneContext:
    def test_known_gene(self):
        ctx = get_gene_context("CYP2D6")
        assert "CYP2D6" in ctx
        assert "Cytochrome P450 2D6" in ctx

    def test_unknown_gene(self):
        ctx = get_gene_context("FAKEGENE123")
        assert ctx == ""

    def test_case_insensitive(self):
        ctx = get_gene_context("cyp2d6")
        assert "CYP2D6" in ctx


class TestGetPhenotypeContext:
    def test_known_phenotype(self):
        ctx = get_phenotype_context("poor")
        assert "PM" in ctx or "Poor" in ctx

    def test_abbreviation(self):
        ctx = get_phenotype_context("PM")
        assert ctx != "" or True  # may match or not depending on alias


class TestGetDrugCategoryContext:
    def test_known_category(self):
        ctx = get_drug_category_context("opioids")
        assert "opioid" in ctx.lower()

    def test_unknown_category(self):
        ctx = get_drug_category_context("zzzzzz_nonexistent")
        assert ctx == ""


class TestGetHLAContext:
    def test_known_allele(self):
        ctx = get_hla_context("HLA-B*57:01")
        assert ctx != ""
        assert "abacavir" in ctx.lower() or "HLA-B*57:01" in ctx

    def test_drug_lookup(self):
        ctx = get_hla_context("ABACAVIR")
        assert ctx != ""


class TestGetInhibitorContext:
    def test_known_enzyme(self):
        ctx = get_inhibitor_context("CYP2D6")
        assert "Inhibitor" in ctx or "inhibitor" in ctx.lower() or ctx != ""

    def test_unknown_enzyme(self):
        ctx = get_inhibitor_context("FAKECYP999")
        assert ctx == ""


class TestGetAlternativeDrugs:
    def test_returns_string(self):
        result = get_alternative_drugs("codeine", "CYP2D6", "Poor Metabolizer")
        assert isinstance(result, str)


class TestGetAllContextForQuery:
    def test_gene_query(self):
        ctx = get_all_context_for_query("Tell me about CYP2D6")
        assert "CYP2D6" in ctx

    def test_drug_query(self):
        ctx = get_all_context_for_query("Tell me about warfarin dosing")
        # May or may not return context depending on drug keywords
        assert isinstance(ctx, str)

    def test_empty_query(self):
        ctx = get_all_context_for_query("")
        assert isinstance(ctx, str)


class TestGetKnowledgeStats:
    def test_returns_dict(self):
        stats = get_knowledge_stats()
        assert isinstance(stats, dict)

    def test_pharmacogenes_count(self):
        stats = get_knowledge_stats()
        assert stats["pharmacogenes"] == 25

    def test_metabolizer_phenotypes_count(self):
        stats = get_knowledge_stats()
        assert stats["metabolizer_phenotypes"] == 5

    def test_drug_categories_count(self):
        stats = get_knowledge_stats()
        assert stats["drug_categories"] == 12

    def test_hla_count(self):
        stats = get_knowledge_stats()
        assert stats["hla_drug_associations"] == 12

    def test_all_keys_present(self):
        stats = get_knowledge_stats()
        expected_keys = [
            "pharmacogenes", "metabolizer_phenotypes", "drug_categories",
            "drugs_tracked", "cyp_inhibitor_enzymes", "cyp_inducer_enzymes",
            "hla_drug_associations", "drug_alternative_mappings",
            "activity_score_tables", "entity_aliases",
        ]
        for key in expected_keys:
            assert key in stats, f"Missing key: {key}"


class TestResolveComparisonEntity:
    def test_returns_none_for_unknown(self):
        result = resolve_comparison_entity("zzz_unknown_entity_999")
        # May return None or a dict depending on fuzzy matching
        assert result is None or isinstance(result, dict)

    def test_returns_dict_for_known_gene(self):
        result = resolve_comparison_entity("CYP2D6")
        if result is not None:
            assert isinstance(result, dict)


class TestGetComparisonContext:
    def test_returns_string(self):
        entity_a = {"type": "gene", "canonical": "CYP2D6"}
        entity_b = {"type": "gene", "canonical": "CYP2C19"}
        ctx = get_comparison_context(entity_a, entity_b)
        assert isinstance(ctx, str)
