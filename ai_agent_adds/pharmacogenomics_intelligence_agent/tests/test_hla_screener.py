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
        hla = {"HLA-B": ["*07:02"], "HLA-A": ["*02:01"]}
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
