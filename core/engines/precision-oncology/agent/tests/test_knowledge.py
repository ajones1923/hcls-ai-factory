"""
Tests for the precision oncology knowledge graph.
==================================================
Validates ACTIONABLE_TARGETS, THERAPY_MAP, RESISTANCE_MAP,
PATHWAY_MAP, BIOMARKER_PANELS, ENTITY_ALIASES, and all
helper functions.
"""

import sys
from pathlib import Path

import pytest

_AGENT_ROOT = Path(__file__).resolve().parents[1]
if str(_AGENT_ROOT) not in sys.path:
    sys.path.insert(0, str(_AGENT_ROOT))

from src.knowledge import (
    ACTIONABLE_TARGETS,
    BIOMARKER_PANELS,
    ENTITY_ALIASES,
    PATHWAY_MAP,
    RESISTANCE_MAP,
    THERAPY_MAP,
    get_biomarker_context,
    get_pathway_context,
    get_resistance_context,
    get_target_context,
    get_therapy_context,
    resolve_comparison_entity,
)


# ═══════════════════════════════════════════════════════════════════════════
# ACTIONABLE_TARGETS Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestActionableTargets:
    """Test ACTIONABLE_TARGETS has expected genes."""

    @pytest.mark.parametrize("gene", [
        "BRAF", "EGFR", "ALK", "ROS1", "KRAS", "HER2",
        "NTRK", "RET", "MET", "FGFR", "PIK3CA", "IDH1",
        "IDH2", "BRCA1", "BRCA2", "TP53", "PTEN", "CDKN2A",
        "STK11", "ESR1",
    ])
    def test_gene_present(self, gene):
        assert gene in ACTIONABLE_TARGETS, f"Missing gene: {gene}"

    def test_minimum_gene_count(self):
        """Should have at least 20 actionable targets."""
        assert len(ACTIONABLE_TARGETS) >= 20

    @pytest.mark.parametrize("gene", ["BRAF", "EGFR", "ALK", "KRAS"])
    def test_gene_has_required_fields(self, gene):
        """Each gene entry should have required fields."""
        target = ACTIONABLE_TARGETS[gene]
        assert "gene" in target
        assert "full_name" in target
        assert "cancer_types" in target
        assert "key_variants" in target
        assert "targeted_therapies" in target
        assert "pathway" in target
        assert "evidence_level" in target
        assert "description" in target

    def test_braf_cancer_types(self):
        """BRAF should include melanoma in cancer types."""
        assert "melanoma" in ACTIONABLE_TARGETS["BRAF"]["cancer_types"]

    def test_egfr_cancer_types(self):
        """EGFR should include NSCLC in cancer types."""
        assert "NSCLC" in ACTIONABLE_TARGETS["EGFR"]["cancer_types"]

    def test_braf_has_therapies(self):
        """BRAF should have targeted therapies."""
        assert len(ACTIONABLE_TARGETS["BRAF"]["targeted_therapies"]) >= 1

    def test_egfr_has_osimertinib(self):
        """EGFR should include osimertinib in targeted therapies."""
        assert "osimertinib" in ACTIONABLE_TARGETS["EGFR"]["targeted_therapies"]

    def test_braf_evidence_level_a(self):
        """BRAF should have evidence level A."""
        assert ACTIONABLE_TARGETS["BRAF"]["evidence_level"] == "A"

    def test_braf_has_resistance_mutations(self):
        """BRAF should have known resistance mutations."""
        assert len(ACTIONABLE_TARGETS["BRAF"]["resistance_mutations"]) >= 1

    def test_brca1_pathway_is_ddr(self):
        """BRCA1 should be in the DDR pathway."""
        assert "DDR" in ACTIONABLE_TARGETS["BRCA1"]["pathway"]

    def test_msi_h_is_tissue_agnostic(self):
        """MSI_H should include tissue-agnostic in cancer types."""
        if "MSI_H" in ACTIONABLE_TARGETS:
            assert "tissue-agnostic" in ACTIONABLE_TARGETS["MSI_H"]["cancer_types"]


# ═══════════════════════════════════════════════════════════════════════════
# THERAPY_MAP Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestTherapyMap:
    """Test THERAPY_MAP has expected drugs."""

    @pytest.mark.parametrize("drug", [
        "osimertinib", "pembrolizumab", "nivolumab", "vemurafenib",
        "dabrafenib", "sotorasib", "adagrasib", "larotrectinib",
        "entrectinib", "selpercatinib", "capmatinib", "olaparib",
        "alpelisib", "ivosidenib",
    ])
    def test_drug_present(self, drug):
        assert drug in THERAPY_MAP, f"Missing drug: {drug}"

    def test_minimum_drug_count(self):
        """Should have at least 15 drugs."""
        assert len(THERAPY_MAP) >= 15

    @pytest.mark.parametrize("drug", ["osimertinib", "pembrolizumab", "dabrafenib"])
    def test_drug_has_required_fields(self, drug):
        """Each drug entry should have required fields."""
        therapy = THERAPY_MAP[drug]
        assert "drug_name" in therapy
        assert "brand_name" in therapy
        assert "category" in therapy
        assert "targets" in therapy
        assert "approved_indications" in therapy
        assert "mechanism" in therapy
        assert "key_trials" in therapy
        assert "evidence_level" in therapy

    def test_osimertinib_targets_egfr(self):
        """Osimertinib should target EGFR."""
        assert "EGFR" in THERAPY_MAP["osimertinib"]["targets"]

    def test_osimertinib_brand_is_tagrisso(self):
        """Osimertinib brand name should be Tagrisso."""
        assert THERAPY_MAP["osimertinib"]["brand_name"] == "Tagrisso"

    def test_pembrolizumab_is_immunotherapy(self):
        """Pembrolizumab should be categorized as immunotherapy."""
        assert THERAPY_MAP["pembrolizumab"]["category"] == "immunotherapy"

    def test_pembrolizumab_targets_pd1(self):
        """Pembrolizumab should target PD-1."""
        assert "PD-1" in THERAPY_MAP["pembrolizumab"]["targets"]


# ═══════════════════════════════════════════════════════════════════════════
# RESISTANCE_MAP Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestResistanceMap:
    """Test RESISTANCE_MAP has expected mechanisms."""

    @pytest.mark.parametrize("drug_class", [
        "osimertinib",
        "BRAF V600 inhibitors",
        "ALK TKIs",
        "KRAS G12C inhibitors",
        "PARP inhibitors",
        "anti-PD-1/PD-L1",
    ])
    def test_drug_class_present(self, drug_class):
        assert drug_class in RESISTANCE_MAP, f"Missing drug class: {drug_class}"

    def test_minimum_entry_count(self):
        """Should have at least 5 resistance entries."""
        assert len(RESISTANCE_MAP) >= 5

    def test_osimertinib_has_c797s(self):
        """Osimertinib resistance should include C797S."""
        mechanisms = RESISTANCE_MAP["osimertinib"]
        mutations = [m["mutation"] for m in mechanisms]
        assert "C797S" in mutations

    def test_each_entry_has_required_fields(self):
        """Each resistance mechanism entry should have required fields."""
        for drug_class, mechanisms in RESISTANCE_MAP.items():
            for m in mechanisms:
                assert "mutation" in m, f"{drug_class} entry missing 'mutation'"
                assert "gene" in m, f"{drug_class} entry missing 'gene'"
                assert "frequency" in m, f"{drug_class} entry missing 'frequency'"

    def test_parp_inhibitor_reversion(self):
        """PARP inhibitor resistance should include reversion mutations."""
        mechanisms = RESISTANCE_MAP["PARP inhibitors"]
        mutations = [m["mutation"] for m in mechanisms]
        assert any("reversion" in m.lower() for m in mutations)


# ═══════════════════════════════════════════════════════════════════════════
# PATHWAY_MAP Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestPathwayMap:
    """Test PATHWAY_MAP has expected pathways."""

    @pytest.mark.parametrize("pathway", [
        "MAPK", "PI3K_AKT_mTOR", "DDR", "cell_cycle",
        "WNT", "JAK_STAT", "apoptosis", "angiogenesis",
        "Notch", "Hedgehog",
    ])
    def test_pathway_present(self, pathway):
        assert pathway in PATHWAY_MAP, f"Missing pathway: {pathway}"

    def test_minimum_pathway_count(self):
        """Should have at least 8 pathways."""
        assert len(PATHWAY_MAP) >= 8

    @pytest.mark.parametrize("pathway", ["MAPK", "PI3K_AKT_mTOR", "DDR"])
    def test_pathway_has_required_fields(self, pathway):
        """Each pathway should have required fields."""
        pw = PATHWAY_MAP[pathway]
        assert "pathway_name" in pw
        assert "key_genes" in pw
        assert "therapeutic_targets" in pw
        assert "cross_talk" in pw
        assert "clinical_relevance" in pw

    def test_mapk_has_kras(self):
        """MAPK pathway should include KRAS."""
        assert "KRAS" in PATHWAY_MAP["MAPK"]["key_genes"]

    def test_mapk_has_braf(self):
        """MAPK pathway should include BRAF."""
        assert "BRAF" in PATHWAY_MAP["MAPK"]["key_genes"]

    def test_ddr_has_brca1(self):
        """DDR pathway should include BRCA1."""
        assert "BRCA1" in PATHWAY_MAP["DDR"]["key_genes"]

    def test_pi3k_has_pik3ca(self):
        """PI3K/AKT/mTOR pathway should include PIK3CA."""
        assert "PIK3CA" in PATHWAY_MAP["PI3K_AKT_mTOR"]["key_genes"]


# ═══════════════════════════════════════════════════════════════════════════
# BIOMARKER_PANELS Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestBiomarkerPanels:
    """Test BIOMARKER_PANELS has expected biomarkers."""

    @pytest.mark.parametrize("biomarker", [
        "TMB-H", "MSI-H", "PD-L1", "HRD", "BRCA_status",
        "EGFR_mutation", "ALK_fusion", "BRAF_V600",
        "KRAS_G12C", "RET_fusion", "NTRK_fusion", "HER2_amp",
    ])
    def test_biomarker_present(self, biomarker):
        assert biomarker in BIOMARKER_PANELS, f"Missing biomarker: {biomarker}"

    def test_minimum_biomarker_count(self):
        """Should have at least 10 biomarker panels."""
        assert len(BIOMARKER_PANELS) >= 10

    @pytest.mark.parametrize("biomarker", ["TMB-H", "MSI-H", "PD-L1"])
    def test_biomarker_has_required_fields(self, biomarker):
        """Each biomarker should have required fields."""
        bm = BIOMARKER_PANELS[biomarker]
        assert "name" in bm
        assert "type" in bm
        assert "testing_method" in bm
        assert "clinical_cutoff" in bm
        assert "evidence_level" in bm
        assert "description" in bm

    def test_tmb_is_predictive(self):
        """TMB-H should be a predictive biomarker."""
        assert BIOMARKER_PANELS["TMB-H"]["type"] == "predictive"

    def test_msi_h_evidence_level_a(self):
        """MSI-H should have evidence level A."""
        assert BIOMARKER_PANELS["MSI-H"]["evidence_level"] == "A"

    def test_pdl1_has_nsclc(self):
        """PD-L1 should include NSCLC in cancer types."""
        assert "NSCLC" in BIOMARKER_PANELS["PD-L1"]["cancer_types"]


# ═══════════════════════════════════════════════════════════════════════════
# ENTITY_ALIASES Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestEntityAliases:
    """Test ENTITY_ALIASES maps correctly."""

    @pytest.mark.parametrize("alias,expected", [
        ("keytruda", "pembrolizumab"),
        ("opdivo", "nivolumab"),
        ("tagrisso", "osimertinib"),
        ("zelboraf", "vemurafenib"),
        ("tafinlar", "dabrafenib"),
        ("lumakras", "sotorasib"),
        ("vitrakvi", "larotrectinib"),
        ("rozlytrek", "entrectinib"),
        ("retevmo", "selpercatinib"),
        ("lynparza", "olaparib"),
        ("enhertu", "trastuzumab_deruxtecan"),
        ("ibrance", "palbociclib"),
    ])
    def test_drug_alias(self, alias, expected):
        assert ENTITY_ALIASES[alias] == expected

    @pytest.mark.parametrize("alias,expected", [
        ("her2", "HER2"),
        ("erbb2", "HER2"),
        ("lkb1", "STK11"),
        ("p16", "CDKN2A"),
    ])
    def test_gene_alias(self, alias, expected):
        assert ENTITY_ALIASES[alias] == expected

    @pytest.mark.parametrize("alias,expected", [
        ("non-small cell", "NSCLC"),
        ("lung adenocarcinoma", "NSCLC"),
        ("crc", "colorectal"),
        ("rcc", "renal cell"),
        ("hcc", "hepatocellular"),
    ])
    def test_cancer_type_alias(self, alias, expected):
        assert ENTITY_ALIASES[alias] == expected

    @pytest.mark.parametrize("alias,expected", [
        ("msi high", "MSI-H"),
        ("dmmr", "MSI-H"),
        ("tmb high", "TMB-H"),
        ("pdl1", "PD-L1"),
    ])
    def test_biomarker_alias(self, alias, expected):
        assert ENTITY_ALIASES[alias] == expected

    def test_minimum_alias_count(self):
        """Should have a substantial number of aliases."""
        assert len(ENTITY_ALIASES) >= 30


# ═══════════════════════════════════════════════════════════════════════════
# Helper Function Tests: get_target_context
# ═══════════════════════════════════════════════════════════════════════════


class TestGetTargetContext:
    """Test get_target_context returns data for known targets."""

    def test_known_gene_returns_context(self):
        context = get_target_context("BRAF")
        assert len(context) > 0
        assert "BRAF" in context

    def test_egfr_context(self):
        context = get_target_context("EGFR")
        assert "EGFR" in context
        assert "osimertinib" in context or "Epidermal" in context

    def test_context_has_cancer_types(self):
        context = get_target_context("BRAF")
        assert "Cancer types:" in context or "melanoma" in context

    def test_context_has_therapies(self):
        context = get_target_context("BRAF")
        assert "Targeted therapies:" in context or "vemurafenib" in context

    def test_context_has_pathway(self):
        context = get_target_context("BRAF")
        assert "Pathway:" in context

    def test_context_has_evidence_level(self):
        context = get_target_context("BRAF")
        assert "Evidence level:" in context

    def test_unknown_gene_returns_empty(self):
        context = get_target_context("UNKNOWNGENE123")
        assert context == ""

    def test_case_insensitive(self):
        context_upper = get_target_context("BRAF")
        context_lower = get_target_context("braf")
        # Both should return non-empty context
        assert len(context_upper) > 0
        assert len(context_lower) > 0

    def test_alias_resolution(self):
        """Should resolve aliases like her2 -> HER2."""
        context = get_target_context("her2")
        assert len(context) > 0


# ═══════════════════════════════════════════════════════════════════════════
# Helper Function Tests: get_therapy_context
# ═══════════════════════════════════════════════════════════════════════════


class TestGetTherapyContext:
    """Test get_therapy_context returns data for known drugs."""

    def test_known_drug_returns_context(self):
        context = get_therapy_context("osimertinib")
        assert len(context) > 0
        assert "osimertinib" in context

    def test_pembrolizumab_context(self):
        context = get_therapy_context("pembrolizumab")
        assert "pembrolizumab" in context
        assert "Keytruda" in context

    def test_context_has_category(self):
        context = get_therapy_context("osimertinib")
        assert "Category:" in context

    def test_context_has_targets(self):
        context = get_therapy_context("osimertinib")
        assert "Targets:" in context or "EGFR" in context

    def test_context_has_mechanism(self):
        context = get_therapy_context("osimertinib")
        assert "Mechanism:" in context

    def test_unknown_drug_returns_empty(self):
        context = get_therapy_context("fake_drug_xyz")
        assert context == ""

    def test_alias_resolution(self):
        """Should resolve brand name aliases like tagrisso -> osimertinib."""
        context = get_therapy_context("tagrisso")
        assert len(context) > 0
        assert "osimertinib" in context


# ═══════════════════════════════════════════════════════════════════════════
# Helper Function Tests: get_resistance_context
# ═══════════════════════════════════════════════════════════════════════════


class TestGetResistanceContext:
    """Test get_resistance_context returns data for known mechanisms."""

    def test_osimertinib_resistance(self):
        context = get_resistance_context("osimertinib")
        assert len(context) > 0
        assert "C797S" in context

    def test_parp_inhibitor_resistance(self):
        context = get_resistance_context("PARP inhibitors")
        assert len(context) > 0
        assert "reversion" in context.lower()

    def test_context_has_mutations(self):
        context = get_resistance_context("osimertinib")
        assert "mutation" in context.lower() or "C797S" in context

    def test_unknown_drug_returns_empty(self):
        context = get_resistance_context("nonexistent_drug_xyz")
        assert context == ""

    def test_partial_name_match(self):
        """Should match via partial drug class name."""
        context = get_resistance_context("BRAF")
        # Should find "BRAF V600 inhibitors"
        assert len(context) > 0


# ═══════════════════════════════════════════════════════════════════════════
# Helper Function Tests: get_pathway_context
# ═══════════════════════════════════════════════════════════════════════════


class TestGetPathwayContext:
    """Test get_pathway_context returns data for known pathways."""

    def test_mapk_pathway(self):
        context = get_pathway_context("MAPK")
        assert len(context) > 0
        assert "MAPK" in context

    def test_pi3k_pathway(self):
        context = get_pathway_context("PI3K_AKT_mTOR")
        assert len(context) > 0
        assert "PI3K" in context

    def test_ddr_pathway(self):
        context = get_pathway_context("DDR")
        assert len(context) > 0

    def test_context_has_key_genes(self):
        context = get_pathway_context("MAPK")
        assert "Key genes:" in context or "KRAS" in context

    def test_context_has_therapeutic_targets(self):
        context = get_pathway_context("MAPK")
        assert "Therapeutic targets:" in context

    def test_context_has_clinical_relevance(self):
        context = get_pathway_context("MAPK")
        assert len(context) > 50  # Should have substantial content

    def test_unknown_pathway_returns_empty(self):
        context = get_pathway_context("nonexistent_pathway_xyz")
        assert context == ""

    def test_partial_name_match(self):
        """Should match via partial pathway name."""
        context = get_pathway_context("DNA Damage")
        assert len(context) > 0


# ═══════════════════════════════════════════════════════════════════════════
# Helper Function Tests: get_biomarker_context
# ═══════════════════════════════════════════════════════════════════════════


class TestGetBiomarkerContext:
    """Test get_biomarker_context returns data for known biomarkers."""

    def test_tmb_h_context(self):
        context = get_biomarker_context("TMB-H")
        assert len(context) > 0
        assert "TMB" in context or "Tumor Mutational" in context

    def test_msi_h_context(self):
        context = get_biomarker_context("MSI-H")
        assert len(context) > 0
        assert "MSI" in context or "Microsatellite" in context

    def test_pdl1_context(self):
        context = get_biomarker_context("PD-L1")
        assert len(context) > 0
        assert "PD-L1" in context

    def test_context_has_testing_method(self):
        context = get_biomarker_context("TMB-H")
        assert "Testing:" in context or "NGS" in context

    def test_context_has_cutoff(self):
        context = get_biomarker_context("TMB-H")
        assert "Cutoff:" in context or "10" in context

    def test_unknown_biomarker_returns_empty(self):
        context = get_biomarker_context("nonexistent_biomarker_xyz")
        assert context == ""


# ═══════════════════════════════════════════════════════════════════════════
# resolve_comparison_entity Tests
# ═══════════════════════════════════════════════════════════════════════════


class TestResolveComparisonEntity:
    """Test resolve_comparison_entity resolves aliases and entity types."""

    def test_resolve_gene_target(self):
        result = resolve_comparison_entity("BRAF")
        assert result is not None
        assert result["type"] == "target"
        assert result["canonical"] == "BRAF"

    def test_resolve_therapy(self):
        result = resolve_comparison_entity("osimertinib")
        assert result is not None
        assert result["type"] == "therapy"
        assert result["canonical"] == "osimertinib"

    def test_resolve_pathway(self):
        result = resolve_comparison_entity("MAPK")
        assert result is not None
        assert result["type"] == "pathway"

    def test_resolve_biomarker(self):
        result = resolve_comparison_entity("TMB-H")
        assert result is not None
        assert result["type"] == "biomarker"

    def test_resolve_brand_name_alias(self):
        """Brand name should resolve to generic drug."""
        result = resolve_comparison_entity("tagrisso")
        assert result is not None
        assert result["type"] == "therapy"
        assert result["canonical"] == "osimertinib"

    def test_resolve_gene_alias(self):
        """Gene alias should resolve to canonical gene."""
        result = resolve_comparison_entity("her2")
        assert result is not None
        assert result["type"] == "target"

    def test_unknown_entity_returns_none(self):
        result = resolve_comparison_entity("completely_unknown_entity_xyz")
        assert result is None

    def test_result_has_data(self):
        """Resolved entity should include the underlying data dict."""
        result = resolve_comparison_entity("EGFR")
        assert result is not None
        assert "data" in result
        assert isinstance(result["data"], dict)
