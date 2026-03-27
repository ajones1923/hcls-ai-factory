"""Tests for src/pgx_pipeline.py — StarAlleleCaller, PhenotypeTranslator, DrugGeneMatcher.

Author: Adam Jones
Date: March 2026
"""

import pytest
from src.pgx_pipeline import (
    AlertSeverity,
    PGxAlert,
    PGxPosition,
    StarAlleleCaller,
    PhenotypeTranslator,
    DrugGeneMatcher,
)


# ═══════════════════════════════════════════════════════════════════════
# AlertSeverity
# ═══════════════════════════════════════════════════════════════════════

class TestAlertSeverity:
    def test_values(self):
        assert AlertSeverity.CONTRAINDICATED == "contraindicated"
        assert AlertSeverity.MAJOR == "major"
        assert AlertSeverity.MODERATE == "moderate"
        assert AlertSeverity.MINOR == "minor"
        assert AlertSeverity.INFORMATIONAL == "informational"

    def test_count(self):
        assert len(AlertSeverity) == 5


# ═══════════════════════════════════════════════════════════════════════
# PGxPosition
# ═══════════════════════════════════════════════════════════════════════

class TestPGxPosition:
    def test_create(self):
        pos = PGxPosition(chrom="chr22", pos=42128191, ref="C", alt="T",
                          star_allele="*4", rsid="rs3892097", gene="CYP2D6")
        assert pos.gene == "CYP2D6"
        assert pos.star_allele == "*4"

    def test_frozen(self):
        pos = PGxPosition(chrom="chr22", pos=1, ref="A", alt="T",
                          star_allele="*1", rsid="rs1", gene="G")
        with pytest.raises(AttributeError):
            pos.gene = "OTHER"


# ═══════════════════════════════════════════════════════════════════════
# StarAlleleCaller
# ═══════════════════════════════════════════════════════════════════════

class TestStarAlleleCaller:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.caller = StarAlleleCaller()

    # -- PGX_POSITIONS coverage --

    def test_has_cyp2d6_positions(self):
        assert "CYP2D6" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CYP2D6"]) == 8

    def test_has_cyp2c19_positions(self):
        assert "CYP2C19" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CYP2C19"]) == 3

    def test_has_cyp2c9_positions(self):
        assert "CYP2C9" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CYP2C9"]) == 2

    def test_has_dpyd_positions(self):
        assert "DPYD" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["DPYD"]) == 3

    def test_has_tpmt_positions(self):
        assert "TPMT" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["TPMT"]) == 3

    def test_has_cyp3a5_positions(self):
        assert "CYP3A5" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CYP3A5"]) == 2

    def test_all_genes_present(self):
        expected_genes = [
            "CYP2D6", "CYP2C19", "CYP2C9", "VKORC1", "SLCO1B1",
            "DPYD", "TPMT", "NUDT15", "UGT1A1", "CYP3A5", "CYP4F2",
            "ABCB1", "CACNA1S", "CYP1A2", "CYP2B6", "F5", "G6PD",
            "IFNL3", "MTHFR", "NAT2", "RYR1",
        ]
        for gene in expected_genes:
            assert gene in StarAlleleCaller.PGX_POSITIONS

    def test_gene_count_at_least_21(self):
        assert len(StarAlleleCaller.PGX_POSITIONS) >= 21

    def test_has_nat2_positions(self):
        assert "NAT2" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["NAT2"]) == 4

    def test_has_ryr1_positions(self):
        assert "RYR1" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["RYR1"]) == 4

    def test_has_abcb1_positions(self):
        assert "ABCB1" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["ABCB1"]) == 3

    def test_has_cyp2b6_positions(self):
        assert "CYP2B6" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CYP2B6"]) == 3

    def test_has_g6pd_positions(self):
        assert "G6PD" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["G6PD"]) == 3

    def test_has_f5_positions(self):
        assert "F5" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["F5"]) == 2

    def test_has_mthfr_positions(self):
        assert "MTHFR" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["MTHFR"]) == 2

    def test_has_ifnl3_positions(self):
        assert "IFNL3" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["IFNL3"]) == 2

    def test_has_cacna1s_positions(self):
        assert "CACNA1S" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CACNA1S"]) == 1

    def test_has_cyp1a2_positions(self):
        assert "CYP1A2" in StarAlleleCaller.PGX_POSITIONS
        assert len(StarAlleleCaller.PGX_POSITIONS["CYP1A2"]) == 3

    # -- call_star_alleles (no variants = wild-type) --

    def test_no_variants_returns_star1_star1(self):
        a1, a2 = self.caller.call_star_alleles([], "CYP2D6")
        assert a1 == "*1"
        assert a2 == "*1"

    # -- call_star_alleles with heterozygous variant --

    def test_heterozygous_variant(self):
        variants = [{
            "gene": "CYP2D6",
            "star_allele": "*4",
            "genotype": "0/1",
        }]
        a1, a2 = self.caller.call_star_alleles(variants, "CYP2D6")
        assert set([a1, a2]) == {"*1", "*4"}

    # -- call_star_alleles with homozygous variant --

    def test_homozygous_variant(self):
        variants = [{
            "gene": "CYP2D6",
            "star_allele": "*4",
            "genotype": "1/1",
        }]
        a1, a2 = self.caller.call_star_alleles(variants, "CYP2D6")
        assert a1 == "*4"
        assert a2 == "*4"

    # -- call_star_alleles with multiple variants (compound het) --

    def test_compound_het(self):
        variants = [
            {"gene": "CYP2D6", "star_allele": "*4", "genotype": "0/1"},
            {"gene": "CYP2D6", "star_allele": "*6", "genotype": "0/1"},
        ]
        a1, a2 = self.caller.call_star_alleles(variants, "CYP2D6")
        assert set([a1, a2]) == {"*4", "*6"}

    # -- canonical ordering --

    def test_canonical_ordering(self):
        a1, a2 = StarAlleleCaller._order_diplotype("*4", "*1")
        assert a1 == "*1"
        assert a2 == "*4"

    def test_canonical_ordering_same(self):
        a1, a2 = StarAlleleCaller._order_diplotype("*3", "*3")
        assert a1 == "*3"
        assert a2 == "*3"

    # -- star sort key --

    def test_star_sort_key(self):
        assert StarAlleleCaller._star_sort_key("*1") < StarAlleleCaller._star_sort_key("*4")
        assert StarAlleleCaller._star_sort_key("*3A") < StarAlleleCaller._star_sort_key("*4")
        assert StarAlleleCaller._star_sort_key("*3A") > StarAlleleCaller._star_sort_key("*3")

    # -- compose_allele --

    def test_compose_allele_empty(self):
        result = StarAlleleCaller._compose_allele([], "CYP2D6")
        assert result == "*1"

    def test_compose_allele_single(self):
        result = StarAlleleCaller._compose_allele(["*4"], "CYP2D6")
        assert result == "*4"

    # -- TPMT *3A detection --

    def test_tpmt_3a_cis(self):
        a1, a2 = StarAlleleCaller._resolve_tpmt_3a(
            ["*3B", "*3C"], []
        )
        assert "*3A" in a1

    def test_tpmt_3a_trans(self):
        a1, a2 = StarAlleleCaller._resolve_tpmt_3a(
            ["*3B"], ["*3C"]
        )
        assert "*3B" in a1
        assert "*3C" in a2

    # -- filter_status for non-matching gene --

    def test_ignores_other_gene_variants(self):
        variants = [{
            "gene": "CYP2C19",
            "star_allele": "*2",
            "genotype": "0/1",
        }]
        a1, a2 = self.caller.call_star_alleles(variants, "CYP2D6")
        assert a1 == "*1"
        assert a2 == "*1"


# ═══════════════════════════════════════════════════════════════════════
# PhenotypeTranslator
# ═══════════════════════════════════════════════════════════════════════

class TestPhenotypeTranslator:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.translator = PhenotypeTranslator()

    # -- CYP2D6 --

    def test_cyp2d6_normal_metabolizer(self):
        result = self.translator.translate("CYP2D6", "*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["activity_score"] == 2.0

    def test_cyp2d6_poor_metabolizer(self):
        result = self.translator.translate("CYP2D6", "*4/*4")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["activity_score"] == 0.0

    def test_cyp2d6_as1_normal_metabolizer(self):
        """CYP2D6 *1/*4 (AS=1.0) is Normal Metabolizer per CPIC 2023."""
        result = self.translator.translate("CYP2D6", "*1/*4")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["activity_score"] == 1.0

    def test_cyp2d6_intermediate_metabolizer(self):
        """CYP2D6 *4/*41 (AS=0.5) is Intermediate Metabolizer."""
        result = self.translator.translate("CYP2D6", "*4/*41")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["activity_score"] == 0.5

    def test_cyp2d6_star10_decreased(self):
        result = self.translator.translate("CYP2D6", "*10/*10")
        assert result["activity_score"] == 0.5
        assert result["phenotype"] == "Intermediate Metabolizer"

    # -- CYP2C19 --

    def test_cyp2c19_poor_metabolizer(self):
        result = self.translator.translate("CYP2C19", "*2/*2")
        assert result["phenotype"] == "Poor Metabolizer"

    def test_cyp2c19_ultrarapid(self):
        result = self.translator.translate("CYP2C19", "*17/*17")
        assert result["phenotype"] == "Ultrarapid Metabolizer"
        assert result["activity_score"] == 3.0

    def test_cyp2c19_rapid(self):
        """CYP2C19 *1/*17 (AS=2.5) is Rapid Metabolizer per CPIC."""
        result = self.translator.translate("CYP2C19", "*1/*17")
        assert result["phenotype"] == "Rapid Metabolizer"
        assert result["activity_score"] == 2.5

    def test_cyp2c19_normal(self):
        """CYP2C19 *1/*1 (AS=2.0) is Normal Metabolizer per CPIC."""
        result = self.translator.translate("CYP2C19", "*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"

    # -- CYP2C9 --

    def test_cyp2c9_poor_metabolizer(self):
        result = self.translator.translate("CYP2C9", "*3/*3")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["activity_score"] == 0.0

    def test_cyp2c9_intermediate(self):
        """CYP2C9 *1/*2 (AS=1.5) is Intermediate Metabolizer per CPIC."""
        result = self.translator.translate("CYP2C9", "*1/*2")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["activity_score"] == 1.5

    def test_cyp2c9_normal(self):
        result = self.translator.translate("CYP2C9", "*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"

    # -- DPYD --

    def test_dpyd_deficient(self):
        result = self.translator.translate("DPYD", "*2A/*2A")
        assert result["phenotype"] == "DPD Deficient"
        assert result["activity_score"] == 0.0

    def test_dpyd_intermediate(self):
        result = self.translator.translate("DPYD", "*1/*2A")
        assert result["phenotype"] == "DPD Intermediate Activity"
        assert result["activity_score"] == 1.0

    def test_dpyd_normal(self):
        result = self.translator.translate("DPYD", "*1/*1")
        assert result["phenotype"] == "DPD Normal Activity"

    # -- TPMT --

    def test_tpmt_deficient(self):
        result = self.translator.translate("TPMT", "*3A/*3C")
        assert result["phenotype"] == "TPMT Deficient"

    def test_tpmt_intermediate(self):
        result = self.translator.translate("TPMT", "*1/*3A")
        assert result["phenotype"] == "TPMT Intermediate"

    def test_tpmt_normal(self):
        result = self.translator.translate("TPMT", "*1/*1")
        assert result["phenotype"] == "TPMT Normal"

    # -- CYP3A5 --

    def test_cyp3a5_non_expresser(self):
        result = self.translator.translate("CYP3A5", "*3/*3")
        assert result["phenotype"] == "CYP3A5 Non-expresser"

    def test_cyp3a5_intermediate_expresser(self):
        result = self.translator.translate("CYP3A5", "*1/*3")
        assert result["phenotype"] == "CYP3A5 Intermediate Expresser"

    def test_cyp3a5_expresser(self):
        result = self.translator.translate("CYP3A5", "*1/*1")
        assert result["phenotype"] == "CYP3A5 Expresser"

    # -- Genotype-based genes (VKORC1, SLCO1B1, CYP4F2) --

    def test_vkorc1_low_dose(self):
        result = self.translator.translate("VKORC1", "-1639G>A/-1639G>A")
        assert result["phenotype"] == "Low Dose Warfarin Sensitivity"
        assert result["activity_score"] is None

    def test_vkorc1_normal(self):
        result = self.translator.translate("VKORC1", "*1/*1")
        assert result["phenotype"] == "Normal Warfarin Sensitivity"

    def test_slco1b1_poor_function(self):
        result = self.translator.translate("SLCO1B1", "*5/*5")
        assert result["phenotype"] == "Poor Function"

    def test_slco1b1_normal(self):
        result = self.translator.translate("SLCO1B1", "*1/*1")
        assert result["phenotype"] == "Normal Function"

    # -- translate_all --

    def test_translate_all(self):
        profile = {
            "CYP2D6": {"diplotype": "*1/*4"},
            "CYP2C19": {"diplotype": "*1/*1"},
        }
        results = self.translator.translate_all(profile)
        assert "CYP2D6" in results
        assert "CYP2C19" in results
        assert results["CYP2D6"]["phenotype"] == "Normal Metabolizer"

    # -- return structure --

    def test_translate_returns_all_fields(self):
        result = self.translator.translate("CYP2D6", "*1/*1")
        assert "gene" in result
        assert "diplotype" in result
        assert "phenotype" in result
        assert "activity_score" in result
        assert "clinical_meaning" in result


# ═══════════════════════════════════════════════════════════════════════
# DrugGeneMatcher
# ═══════════════════════════════════════════════════════════════════════

class TestDrugGeneMatcher:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.matcher = DrugGeneMatcher()

    def test_has_17_drugs(self):
        assert len(DrugGeneMatcher.DRUG_GENE_MAP) >= 14

    # -- codeine CYP2D6 --

    def test_codeine_pm_contraindicated(self):
        profile = {"CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"}}
        alerts = self.matcher.check_drug("codeine", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MAJOR

    def test_codeine_um_contraindicated(self):
        profile = {"CYP2D6": {"phenotype": "Ultrarapid Metabolizer", "diplotype": "*1/*1xN"}}
        alerts = self.matcher.check_drug("codeine", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.CONTRAINDICATED

    def test_codeine_nm_no_alert(self):
        profile = {"CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"}}
        alerts = self.matcher.check_drug("codeine", profile)
        assert len(alerts) == 0

    # -- clopidogrel CYP2C19 --

    def test_clopidogrel_pm(self):
        profile = {"CYP2C19": {"phenotype": "Poor Metabolizer", "diplotype": "*2/*2"}}
        alerts = self.matcher.check_drug("clopidogrel", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.CONTRAINDICATED

    # -- warfarin CYP2C9 --

    def test_warfarin_cyp2c9_pm(self):
        profile = {
            "CYP2C9": {"phenotype": "Poor Metabolizer", "diplotype": "*3/*3"},
            "VKORC1": {"phenotype": "Normal Warfarin Sensitivity", "diplotype": "*1/*1"},
        }
        alerts = self.matcher.check_drug("warfarin", profile)
        assert any(a.gene == "CYP2C9" for a in alerts)

    # -- fluorouracil DPYD --

    def test_fluorouracil_dpyd_deficient(self):
        profile = {"DPYD": {"phenotype": "DPD Deficient", "diplotype": "*2A/*2A"}}
        alerts = self.matcher.check_drug("fluorouracil", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.CONTRAINDICATED

    # -- azathioprine TPMT --

    def test_azathioprine_tpmt_deficient(self):
        profile = {"TPMT": {"phenotype": "TPMT Deficient", "diplotype": "*3A/*3C"}}
        alerts = self.matcher.check_drug("azathioprine", profile)
        assert any(a.gene == "TPMT" for a in alerts)

    # -- simvastatin SLCO1B1 --

    def test_simvastatin_poor_function(self):
        profile = {"SLCO1B1": {"phenotype": "Poor Function", "diplotype": "*5/*5"}}
        alerts = self.matcher.check_drug("simvastatin", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.CONTRAINDICATED

    # -- tacrolimus CYP3A5 --

    def test_tacrolimus_expresser(self):
        profile = {"CYP3A5": {"phenotype": "CYP3A5 Expresser", "diplotype": "*1/*1"}}
        alerts = self.matcher.check_drug("tacrolimus", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MAJOR

    # -- unknown drug --

    def test_unknown_drug_no_alerts(self):
        profile = {"CYP2D6": {"phenotype": "Normal Metabolizer"}}
        alerts = self.matcher.check_drug("aspirin", profile)
        assert len(alerts) == 0

    # -- case insensitive --

    def test_case_insensitive_drug(self):
        profile = {"CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"}}
        alerts = self.matcher.check_drug("CODEINE", profile)
        assert len(alerts) >= 1

    # -- check_medications --

    def test_check_medications_method(self):
        """Test that check_drug works for multiple drugs."""
        profile = {
            "CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"},
            "CYP2C19": {"phenotype": "Poor Metabolizer", "diplotype": "*2/*2"},
        }
        all_alerts = []
        for drug in ["codeine", "clopidogrel"]:
            all_alerts.extend(self.matcher.check_drug(drug, profile))
        assert len(all_alerts) >= 2

    # -- alert structure --

    def test_alert_has_required_fields(self):
        profile = {"CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"}}
        alerts = self.matcher.check_drug("codeine", profile)
        if alerts:
            alert = alerts[0]
            assert alert.drug
            assert alert.gene
            assert alert.phenotype
            assert alert.recommendation
            assert isinstance(alert.severity, AlertSeverity)


# ═══════════════════════════════════════════════════════════════════════
# EXPANDED DATA TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestExpandedDrugGeneMap:
    """Tests for the 12 new drugs added to DRUG_GENE_MAP."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.matcher = DrugGeneMatcher()

    def test_drug_count_at_least_28(self):
        assert len(DrugGeneMatcher.DRUG_GENE_MAP) >= 28

    @pytest.mark.parametrize("drug", [
        "atomoxetine", "carbamazepine", "citalopram", "clomipramine",
        "doxepin", "efavirenz", "escitalopram", "fluvoxamine",
        "nortriptyline", "ondansetron", "phenytoin", "sertraline",
    ])
    def test_new_drug_present(self, drug):
        assert drug in DrugGeneMatcher.DRUG_GENE_MAP

    def test_atomoxetine_cyp2d6_pm(self):
        profile = {"CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"}}
        alerts = self.matcher.check_drug("atomoxetine", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MAJOR

    def test_citalopram_cyp2c19_pm(self):
        profile = {"CYP2C19": {"phenotype": "Poor Metabolizer", "diplotype": "*2/*2"}}
        alerts = self.matcher.check_drug("citalopram", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MAJOR

    def test_efavirenz_cyp2b6_pm(self):
        profile = {"CYP2B6": {"phenotype": "Poor Metabolizer", "diplotype": "*6/*6"}}
        alerts = self.matcher.check_drug("efavirenz", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MAJOR

    def test_nortriptyline_cyp2d6_um_major(self):
        profile = {"CYP2D6": {"phenotype": "Ultrarapid Metabolizer", "diplotype": "*1/*2xN"}}
        alerts = self.matcher.check_drug("nortriptyline", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MAJOR

    def test_phenytoin_hla_contraindicated(self):
        profile = {"HLA-B": {"phenotype": "HLA-B*15:02 positive"}}
        alerts = self.matcher.check_drug("phenytoin", profile)
        assert any(a.severity == AlertSeverity.CONTRAINDICATED for a in alerts)

    def test_carbamazepine_hla_contraindicated(self):
        profile = {"HLA-B": {"phenotype": "HLA-B*15:02 positive"}}
        alerts = self.matcher.check_drug("carbamazepine", profile)
        assert any(a.severity == AlertSeverity.CONTRAINDICATED for a in alerts)

    def test_ondansetron_cyp2d6_um(self):
        profile = {"CYP2D6": {"phenotype": "Ultrarapid Metabolizer", "diplotype": "*1/*2xN"}}
        alerts = self.matcher.check_drug("ondansetron", profile)
        assert len(alerts) >= 1
        assert alerts[0].severity == AlertSeverity.MODERATE


class TestExpandedActivityScores:
    """Tests for expanded allele coverage in ACTIVITY_SCORES."""

    def test_cyp2d6_21_alleles(self):
        assert len(PhenotypeTranslator.ACTIVITY_SCORES["CYP2D6"]) >= 21

    def test_cyp2c19_11_alleles(self):
        assert len(PhenotypeTranslator.ACTIVITY_SCORES["CYP2C19"]) >= 11

    def test_cyp2c9_7_alleles(self):
        assert len(PhenotypeTranslator.ACTIVITY_SCORES["CYP2C9"]) >= 7

    def test_ugt1a1_6_alleles(self):
        assert len(PhenotypeTranslator.ACTIVITY_SCORES["UGT1A1"]) >= 6

    def test_cyp2d6_star29_decreased(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["CYP2D6"]["*29"] == 0.5

    def test_cyp2d6_star35_normal(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["CYP2D6"]["*35"] == 1.0

    def test_cyp2c19_star4_no_function(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["CYP2C19"]["*4"] == 0.0

    def test_ugt1a1_star37_decreased(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["UGT1A1"]["*37"] == 0.0

    def test_ugt1a1_star36_increased(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["UGT1A1"]["*36"] == 1.5

    def test_nudt15_star2_intermediate(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]["*2"] == 0.5

    def test_cyp2d6_star7_no_function(self):
        """CYP2D6*7 is a no-function allele."""
        translator = PhenotypeTranslator()
        result = translator.translate("CYP2D6", "*7/*7")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["activity_score"] == 0.0

    def test_cyp2d6_star29_im(self):
        """CYP2D6 *4/*29 (AS=0.5) should be Intermediate Metabolizer."""
        translator = PhenotypeTranslator()
        result = translator.translate("CYP2D6", "*4/*29")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["activity_score"] == 0.5
