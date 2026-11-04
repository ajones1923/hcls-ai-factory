"""Tests for src/hla_screener.py — HLAScreener.

Tests all 12 HLA-drug associations, positive/negative screening,
severity levels, alternative recommendations.

Author: Adam Jones
Date: March 2026
"""

import pytest
from src.hla_screener import (
    HLAStatus,
    ReactionSeverity,
    HLAScreenResult,
    HLA_DRUG_ASSOCIATIONS,
    HLAScreener,
)


# ═══════════════════════════════════════════════════════════════════════
# ENUMS
# ═══════════════════════════════════════════════════════════════════════

class TestHLAStatusEnum:
    def test_values(self):
        assert HLAStatus.SAFE == "SAFE"
        assert HLAStatus.CONTRAINDICATED == "CONTRAINDICATED"
        assert HLAStatus.HIGH_RISK == "HIGH_RISK"
        assert HLAStatus.UNKNOWN == "UNKNOWN"


class TestReactionSeverityEnum:
    def test_values(self):
        assert ReactionSeverity.FATAL == "fatal"
        assert ReactionSeverity.SEVERE == "severe"
        assert ReactionSeverity.MODERATE == "moderate"
        assert ReactionSeverity.MILD == "mild"


# ═══════════════════════════════════════════════════════════════════════
# HLA_DRUG_ASSOCIATIONS KNOWLEDGE BASE
# ═══════════════════════════════════════════════════════════════════════

class TestHLADrugAssociationsKB:
    def test_abacavir_present(self):
        assert "abacavir" in HLA_DRUG_ASSOCIATIONS

    def test_carbamazepine_has_two_alleles(self):
        assert "carbamazepine" in HLA_DRUG_ASSOCIATIONS
        assert len(HLA_DRUG_ASSOCIATIONS["carbamazepine"]) == 2

    def test_nevirapine_has_two_alleles(self):
        assert "nevirapine" in HLA_DRUG_ASSOCIATIONS
        assert len(HLA_DRUG_ASSOCIATIONS["nevirapine"]) == 2

    @pytest.mark.parametrize("drug", [
        "abacavir", "carbamazepine", "oxcarbazepine", "phenytoin",
        "allopurinol", "flucloxacillin", "lamotrigine", "dapsone",
        "ticlopidine", "nevirapine",
    ])
    def test_drug_present(self, drug):
        assert drug in HLA_DRUG_ASSOCIATIONS

    def test_abacavir_hla_b57_01(self):
        assocs = HLA_DRUG_ASSOCIATIONS["abacavir"]
        assert assocs[0]["hla_allele"] == "HLA-B*57:01"
        assert assocs[0]["severity"] == ReactionSeverity.SEVERE


# ═══════════════════════════════════════════════════════════════════════
# HLAScreener
# ═══════════════════════════════════════════════════════════════════════

class TestHLAScreener:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.screener = HLAScreener()

    # -- screen_drug: positive carrier --

    def test_abacavir_positive(self):
        hla = {"HLA-B": ["*57:01", "*44:02"]}
        results = self.screener.screen_drug("abacavir", hla)
        assert len(results) >= 1
        assert results[0].status == HLAStatus.CONTRAINDICATED
        assert results[0].hla_allele == "HLA-B*57:01"
        assert len(results[0].alternatives) > 0

    def test_carbamazepine_hla_b15_02_positive(self):
        hla = {"HLA-B": ["*15:02", "*44:02"]}
        results = self.screener.screen_drug("carbamazepine", hla)
        assert any(r.hla_allele == "HLA-B*15:02" for r in results)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    def test_carbamazepine_hla_a31_01_positive(self):
        hla = {"HLA-A": ["*31:01", "*02:01"], "HLA-B": ["*44:02"]}
        results = self.screener.screen_drug("carbamazepine", hla)
        assert any(r.hla_allele == "HLA-A*31:01" for r in results)

    def test_allopurinol_positive(self):
        hla = {"HLA-B": ["*58:01"]}
        results = self.screener.screen_drug("allopurinol", hla)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    def test_dapsone_positive(self):
        hla = {"HLA-B": ["*13:01"]}
        results = self.screener.screen_drug("dapsone", hla)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    # -- screen_drug: negative carrier (safe) --

    def test_abacavir_negative(self):
        hla = {"HLA-B": ["*44:02", "*07:02"]}
        results = self.screener.screen_drug("abacavir", hla)
        assert results[0].status == HLAStatus.SAFE

    def test_carbamazepine_negative(self):
        hla = {"HLA-B": ["*07:02"], "HLA-A": ["*02:01"]}
        results = self.screener.screen_drug("carbamazepine", hla)
        assert all(r.status == HLAStatus.SAFE for r in results)

    # -- screen_drug: unknown drug --

    def test_unknown_drug(self):
        hla = {"HLA-B": ["*57:01"]}
        results = self.screener.screen_drug("metformin", hla)
        assert len(results) == 1
        assert results[0].status == HLAStatus.UNKNOWN

    # -- screen_drug: case insensitive --

    def test_case_insensitive(self):
        hla = {"HLA-B": ["*57:01"]}
        results = self.screener.screen_drug("ABACAVIR", hla)
        assert results[0].status == HLAStatus.CONTRAINDICATED

    # -- screen_all_drugs --

    def test_screen_all_drugs_b57_01(self):
        hla = {"HLA-B": ["*57:01"]}
        results = self.screener.screen_all_drugs(hla)
        drugs = {r.drug for r in results}
        assert "abacavir" in drugs
        assert "flucloxacillin" in drugs

    def test_screen_all_drugs_sorted_by_severity(self):
        hla = {"HLA-B": ["*57:01", "*15:02"]}
        results = self.screener.screen_all_drugs(hla)
        if len(results) >= 2:
            severity_order = {
                ReactionSeverity.FATAL: 0,
                ReactionSeverity.SEVERE: 1,
                ReactionSeverity.MODERATE: 2,
                ReactionSeverity.MILD: 3,
            }
            for i in range(len(results) - 1):
                assert severity_order.get(results[i].severity, 99) <= \
                       severity_order.get(results[i + 1].severity, 99)

    def test_screen_all_drugs_no_risk_alleles(self):
        """A clean panel now requires every locus in the KB to be typed.

        CORRECTED (pgx audit §2): this case previously passed typing at
        HLA-A and HLA-B only and asserted an empty panel, which pinned the
        defect — HLA-DRB1 and HLA-DQB1 were never typed, so the DRB1/DQB1
        drugs (lamotrigine, clozapine) were silently reported as clean.
        The intent of the test is unchanged ("no risk alleles ⇒ no
        actionable results"); it now supplies the typing that intent
        assumes.  An untyped locus is covered by
        TestUntypedHLAIsNotNegative below.
        """
        hla = {
            "HLA-A": ["*02:01"],
            "HLA-B": ["*07:02"],
            "HLA-DRB1": ["*15:01"],
            "HLA-DQB1": ["*06:02"],
        }
        results = self.screener.screen_all_drugs(hla)
        assert len(results) == 0  # No actionable results

    # -- allele matching --

    def test_prefix_matching(self):
        """Patient allele HLA-B*57:01:01 should match risk allele HLA-B*57:01."""
        hla = {"HLA-B": ["*57:01:01"]}
        results = self.screener.screen_drug("abacavir", hla)
        assert results[0].status == HLAStatus.CONTRAINDICATED

    # -- private helpers --

    def test_flatten_hla_typing(self):
        hla = {"HLA-B": ["*57:01", "*44:02"], "HLA-A": ["*31:01"]}
        flat = HLAScreener._flatten_hla_typing(hla)
        assert "HLA-B*57:01" in flat
        assert "HLA-B*44:02" in flat
        assert "HLA-A*31:01" in flat

    def test_flatten_without_star(self):
        hla = {"HLA-B": ["57:01"]}
        flat = HLAScreener._flatten_hla_typing(hla)
        assert "HLA-B*57:01" in flat

    def test_check_allele_match_exact(self):
        assert HLAScreener._check_allele_match("HLA-B*57:01", {"HLA-B*57:01"})

    def test_check_allele_match_prefix(self):
        assert HLAScreener._check_allele_match("HLA-B*57:01", {"HLA-B*57:01:01"})

    def test_check_allele_no_match(self):
        assert not HLAScreener._check_allele_match("HLA-B*57:01", {"HLA-B*44:02"})

    def test_determine_status_fatal(self):
        assert HLAScreener._determine_status(ReactionSeverity.FATAL) == HLAStatus.CONTRAINDICATED

    def test_determine_status_severe(self):
        assert HLAScreener._determine_status(ReactionSeverity.SEVERE) == HLAStatus.CONTRAINDICATED

    def test_determine_status_moderate(self):
        assert HLAScreener._determine_status(ReactionSeverity.MODERATE) == HLAStatus.HIGH_RISK

    # -- HLAScreenResult dataclass --

    def test_screen_result_creation(self):
        r = HLAScreenResult(drug="abacavir", hla_allele="HLA-B*57:01",
                            status=HLAStatus.CONTRAINDICATED,
                            reaction_type="AHS",
                            severity=ReactionSeverity.SEVERE,
                            recommendation="do not prescribe",
                            alternatives=["tenofovir"])
        assert r.drug == "abacavir"
        assert len(r.alternatives) == 1

    # -- evidence_level in results --

    def test_result_has_evidence_level(self):
        hla = {"HLA-B": ["*57:01"]}
        results = self.screener.screen_drug("abacavir", hla)
        assert results[0].evidence_level != ""

    # -- population_risk in results --

    def test_result_has_population_risk(self):
        hla = {"HLA-B": ["*57:01"]}
        results = self.screener.screen_drug("abacavir", hla)
        assert results[0].population_risk != ""


# ═══════════════════════════════════════════════════════════════════════
# NEW HLA-DRUG ASSOCIATIONS TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestUntypedHLAIsNotNegative:
    """Regression for pgx audit §2 — untyped must never read as negative.

    Primary source: CPIC's own recommendation table keys abacavir on
    HLA-B "*57:01 negative" / "*57:01 positive" and nothing else
    (api.cpicpgx.org/v1/recommendation, drugid RxNorm:190521, retrieved
    2026-11-03); the FDA abacavir boxed warning likewise requires a
    result before prescribing.  An absent test satisfies neither row.

    Before the fix, screen_drug('abacavir', {}) returned HLAStatus.SAFE
    with the text "abacavir is safe to prescribe based on HLA screening",
    and screen_all_drugs({}) returned [] — a clean panel — for a patient
    who had never been HLA typed.  get_hla_from_vcf returns {} for every
    ordinary WGS/WES VCF, so this was the default path.
    """

    @pytest.fixture(autouse=True)
    def setup(self):
        self.screener = HLAScreener()

    @pytest.mark.parametrize("drug", ["abacavir", "carbamazepine", "allopurinol"])
    def test_no_typing_is_unknown_not_safe(self, drug):
        results = self.screener.screen_drug(drug, {})
        assert len(results) == 1
        assert results[0].status == HLAStatus.UNKNOWN
        assert results[0].status != HLAStatus.SAFE
        assert "not a negative result" in results[0].recommendation.lower()
        assert "safe to prescribe" not in results[0].recommendation.lower()

    def test_untyped_result_names_the_locus(self):
        results = self.screener.screen_drug("abacavir", {})
        assert "HLA-B" in results[0].recommendation

    def test_empty_allele_list_is_also_untyped(self):
        """A locus key present with an empty list is still not a result."""
        results = self.screener.screen_drug("abacavir", {"HLA-B": []})
        assert results[0].status == HLAStatus.UNKNOWN

    def test_wrong_locus_typed_does_not_clear_the_drug(self):
        """HLA-A typing says nothing about HLA-B*57:01."""
        results = self.screener.screen_drug("abacavir", {"HLA-A": ["*02:01"]})
        assert results[0].status == HLAStatus.UNKNOWN

    def test_partial_typing_clears_only_the_typed_locus(self):
        """Carbamazepine needs HLA-B *and* HLA-A; HLA-B alone is not enough."""
        results = self.screener.screen_drug(
            "carbamazepine", {"HLA-B": ["*07:02"]})
        assert results[0].status == HLAStatus.UNKNOWN
        assert "HLA-A" in results[0].recommendation

    def test_typed_and_negative_is_still_safe(self):
        """The fix must not turn a real negative into an UNKNOWN."""
        results = self.screener.screen_drug(
            "abacavir", {"HLA-B": ["*44:02", "*07:02"]})
        assert results[0].status == HLAStatus.SAFE

    def test_untyped_panel_is_not_a_clean_panel(self):
        results = self.screener.screen_all_drugs({})
        assert len(results) > 0
        assert all(r.status == HLAStatus.UNKNOWN for r in results)
        assert "abacavir" in {r.drug for r in results}

    def test_untyped_panel_flags_every_drug_in_the_kb(self):
        results = self.screener.screen_all_drugs({})
        assert {r.drug for r in results} == set(HLA_DRUG_ASSOCIATIONS)

    def test_carrier_still_wins_over_untyped(self):
        """A positive result at one locus is reported, not masked."""
        results = self.screener.screen_drug("carbamazepine", {"HLA-B": ["*15:02"]})
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)
        assert all(r.status != HLAStatus.SAFE for r in results)

    def test_required_loci_derived_from_kb(self):
        assert HLAScreener._required_loci(
            HLA_DRUG_ASSOCIATIONS["carbamazepine"]) == ["HLA-B", "HLA-A"]
        assert HLAScreener._required_loci(
            HLA_DRUG_ASSOCIATIONS["clozapine"]) == ["HLA-DQB1"]

    def test_negative_result_carries_its_own_caveat(self):
        """A typed negative is not a guarantee; CPIC says so on every row."""
        results = self.screener.screen_drug("abacavir", {"HLA-B": ["*44:02"]})
        assert "does not eliminate the risk" in results[0].recommendation


class TestExpandedHLAAssociations:
    """Tests for the 5 new HLA-drug associations."""

    @pytest.mark.parametrize("drug", [
        "sulfasalazine", "methazolamide", "clozapine",
        "trimethoprim-sulfamethoxazole", "minocycline",
    ])
    def test_new_drug_present(self, drug):
        assert drug in HLA_DRUG_ASSOCIATIONS

    def test_sulfasalazine_hla_b13_01(self):
        assocs = HLA_DRUG_ASSOCIATIONS["sulfasalazine"]
        assert assocs[0]["hla_allele"] == "HLA-B*13:01"
        assert assocs[0]["severity"] == ReactionSeverity.SEVERE

    def test_methazolamide_fatal(self):
        assocs = HLA_DRUG_ASSOCIATIONS["methazolamide"]
        assert assocs[0]["severity"] == ReactionSeverity.FATAL

    def test_clozapine_agranulocytosis(self):
        assocs = HLA_DRUG_ASSOCIATIONS["clozapine"]
        assert "agranulocytosis" in assocs[0]["reaction_type"].lower()

    def test_total_drug_count(self):
        assert len(HLA_DRUG_ASSOCIATIONS) >= 15


class TestExpandedHLAScreening:
    """Tests for screening the new drugs."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.screener = HLAScreener()

    def test_sulfasalazine_positive(self):
        hla = {"HLA-B": ["*13:01"]}
        results = self.screener.screen_drug("sulfasalazine", hla)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    def test_methazolamide_positive(self):
        hla = {"HLA-B": ["*59:01"]}
        results = self.screener.screen_drug("methazolamide", hla)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    def test_clozapine_positive(self):
        hla = {"HLA-DQB1": ["*05:02"]}
        results = self.screener.screen_drug("clozapine", hla)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    def test_minocycline_positive(self):
        hla = {"HLA-B": ["*35:02"]}
        results = self.screener.screen_drug("minocycline", hla)
        assert any(r.status == HLAStatus.CONTRAINDICATED for r in results)

    def test_sulfasalazine_negative(self):
        hla = {"HLA-B": ["*44:02"]}
        results = self.screener.screen_drug("sulfasalazine", hla)
        assert results[0].status == HLAStatus.SAFE

    def test_screen_all_includes_new_drugs(self):
        """Panel screen with risk alleles should include new drug associations."""
        hla = {"HLA-B": ["*13:01", "*59:01"]}
        results = self.screener.screen_all_drugs(hla)
        drugs_found = {r.drug for r in results}
        assert "sulfasalazine" in drugs_found
        assert "methazolamide" in drugs_found
