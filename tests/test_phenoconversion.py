"""Tests for src/phenoconversion.py — PhenoconversionDetector.

Tests all 30+ inhibitors, 8 inducers, phenotype shift for each enzyme,
substrate identification, and combined medication review.

Author: Adam Jones
Date: March 2026
"""

import pytest
from src.phenoconversion import (
    InhibitorStrength,
    InducerStrength,
    CYP_INHIBITORS,
    CYP_INDUCERS,
    CYP_SUBSTRATES,
    INHIBITION_SHIFT,
    INDUCTION_SHIFT,
    PhenoconversionAlert,
    PhenoconversionDetector,
)


# ═══════════════════════════════════════════════════════════════════════
# KNOWLEDGE BASE TESTS
# ═══════════════════════════════════════════════════════════════════════

class TestCYPInhibitorKnowledgeBase:
    def test_has_30_plus_inhibitors(self):
        assert len(CYP_INHIBITORS) >= 30

    @pytest.mark.parametrize("drug,enzyme", [
        ("fluoxetine", "CYP2D6"),
        ("paroxetine", "CYP2D6"),
        ("bupropion", "CYP2D6"),
        ("quinidine", "CYP2D6"),
        ("terbinafine", "CYP2D6"),
        ("duloxetine", "CYP2D6"),
        ("fluvoxamine", "CYP2C19"),
        ("fluconazole", "CYP2C19"),
        ("omeprazole", "CYP2C19"),
        ("voriconazole", "CYP2C19"),
        ("amiodarone", "CYP2C9"),
        ("miconazole", "CYP2C9"),
        ("ketoconazole", "CYP3A5"),
        ("clarithromycin", "CYP3A5"),
        ("diltiazem", "CYP3A5"),
    ])
    def test_inhibitor_enzyme_pair(self, drug, enzyme):
        assert drug in CYP_INHIBITORS
        assert enzyme in CYP_INHIBITORS[drug]

    def test_fluoxetine_is_strong_2d6(self):
        assert CYP_INHIBITORS["fluoxetine"]["CYP2D6"] == InhibitorStrength.STRONG

    def test_duloxetine_is_moderate_2d6(self):
        assert CYP_INHIBITORS["duloxetine"]["CYP2D6"] == InhibitorStrength.MODERATE

    def test_celecoxib_is_weak_2d6(self):
        assert CYP_INHIBITORS["celecoxib"]["CYP2D6"] == InhibitorStrength.WEAK


class TestCYPInducerKnowledgeBase:
    def test_has_10_plus_inducers(self):
        assert len(CYP_INDUCERS) >= 10

    @pytest.mark.parametrize("drug", [
        "rifampin", "carbamazepine", "phenytoin", "phenobarbital",
        "st_johns_wort", "efavirenz", "modafinil", "dexamethasone",
    ])
    def test_inducer_present(self, drug):
        assert drug in CYP_INDUCERS

    def test_rifampin_strong_multiple(self):
        rif = CYP_INDUCERS["rifampin"]
        assert rif["CYP3A5"] == InducerStrength.STRONG
        assert rif["CYP2C19"] == InducerStrength.STRONG


class TestCYPSubstrates:
    def test_has_5_enzymes(self):
        assert len(CYP_SUBSTRATES) >= 5

    def test_cyp2d6_substrates(self):
        subs = CYP_SUBSTRATES["CYP2D6"]
        assert "codeine" in subs
        assert "tamoxifen" in subs
        assert "metoprolol" in subs

    def test_cyp2c19_substrates(self):
        subs = CYP_SUBSTRATES["CYP2C19"]
        assert "clopidogrel" in subs

    def test_cyp2c9_substrates(self):
        subs = CYP_SUBSTRATES["CYP2C9"]
        assert "warfarin" in subs

    def test_cyp3a5_substrates(self):
        subs = CYP_SUBSTRATES["CYP3A5"]
        assert "tacrolimus" in subs


class TestInhibitionShift:
    def test_nm_strong_becomes_pm(self):
        assert INHIBITION_SHIFT["Normal Metabolizer"][InhibitorStrength.STRONG] == "Poor Metabolizer"

    def test_nm_moderate_becomes_im(self):
        assert INHIBITION_SHIFT["Normal Metabolizer"][InhibitorStrength.MODERATE] == "Intermediate Metabolizer"

    def test_nm_weak_stays_nm(self):
        assert INHIBITION_SHIFT["Normal Metabolizer"][InhibitorStrength.WEAK] == "Normal Metabolizer"

    def test_im_strong_becomes_pm(self):
        assert INHIBITION_SHIFT["Intermediate Metabolizer"][InhibitorStrength.STRONG] == "Poor Metabolizer"

    def test_pm_stays_pm(self):
        assert INHIBITION_SHIFT["Poor Metabolizer"][InhibitorStrength.STRONG] == "Poor Metabolizer"

    def test_um_strong_becomes_nm(self):
        assert INHIBITION_SHIFT["Ultrarapid Metabolizer"][InhibitorStrength.STRONG] == "Normal Metabolizer"

    def test_cyp3a5_expresser_strong(self):
        assert INHIBITION_SHIFT["CYP3A5 Expresser"][InhibitorStrength.STRONG] == "CYP3A5 Non-expresser"


class TestInductionShift:
    def test_nm_strong_becomes_um(self):
        assert INDUCTION_SHIFT["Normal Metabolizer"][InducerStrength.STRONG] == "Ultrarapid Metabolizer"

    def test_pm_strong_becomes_im(self):
        assert INDUCTION_SHIFT["Poor Metabolizer"][InducerStrength.STRONG] == "Intermediate Metabolizer"

    def test_im_strong_becomes_nm(self):
        assert INDUCTION_SHIFT["Intermediate Metabolizer"][InducerStrength.STRONG] == "Normal Metabolizer"


# ═══════════════════════════════════════════════════════════════════════
# PhenoconversionDetector
# ═══════════════════════════════════════════════════════════════════════

class TestPhenoconversionDetector:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.detector = PhenoconversionDetector()

    # -- detect() --

    def test_detect_fluoxetine_shifts_cyp2d6(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine"])
        assert len(alerts) >= 1
        alert = alerts[0]
        assert alert.gene == "CYP2D6"
        assert alert.genetic_phenotype == "Normal Metabolizer"
        assert alert.effective_phenotype == "Poor Metabolizer"
        assert alert.precipitant_drug == "fluoxetine"
        assert alert.inhibitor_strength == "strong"

    def test_detect_no_shift_for_pm(self):
        profile = {
            "CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine"])
        assert len(alerts) == 0  # PM cannot shift further

    def test_detect_moderate_inhibitor_shifts_nm_to_im(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["duloxetine"])
        assert len(alerts) >= 1
        assert alerts[0].effective_phenotype == "Intermediate Metabolizer"

    def test_detect_weak_inhibitor_no_shift_nm(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["celecoxib"])
        assert len(alerts) == 0  # Weak doesn't shift NM

    def test_detect_inducer_rifampin(self):
        profile = {
            "CYP2C19": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["rifampin"])
        cyp2c19_alerts = [a for a in alerts if a.gene == "CYP2C19"]
        assert len(cyp2c19_alerts) >= 1
        assert cyp2c19_alerts[0].effective_phenotype == "Ultrarapid Metabolizer"

    def test_detect_affected_substrates(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine", "codeine"])
        assert len(alerts) >= 1
        assert "codeine" in alerts[0].affected_substrates

    def test_detect_no_matching_gene(self):
        profile = {
            "CYP3A5": {"phenotype": "CYP3A5 Non-expresser", "diplotype": "*3/*3"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine"])
        # fluoxetine inhibits CYP2D6, not CYP3A5
        cyp3a5_alerts = [a for a in alerts if a.gene == "CYP3A5"]
        assert len(cyp3a5_alerts) == 0

    def test_detect_multiple_inhibitors(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
            "CYP2C19": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine", "voriconazole"])
        genes_affected = {a.gene for a in alerts}
        assert "CYP2D6" in genes_affected
        assert "CYP2C19" in genes_affected

    # -- get_effective_phenotype() --

    def test_effective_phenotype_no_meds(self):
        result = self.detector.get_effective_phenotype(
            "CYP2D6", "Normal Metabolizer", []
        )
        assert result == "Normal Metabolizer"

    def test_effective_phenotype_strong_inhibitor(self):
        result = self.detector.get_effective_phenotype(
            "CYP2D6", "Normal Metabolizer", ["fluoxetine"]
        )
        assert result == "Poor Metabolizer"

    def test_effective_phenotype_inhibition_over_induction(self):
        # When both present, inhibition takes precedence
        result = self.detector.get_effective_phenotype(
            "CYP2C19", "Normal Metabolizer", ["voriconazole", "rifampin"]
        )
        # voriconazole is strong CYP2C19 inhibitor, rifampin is strong inducer
        assert result == "Poor Metabolizer"

    def test_effective_phenotype_strongest_wins(self):
        result = self.detector.get_effective_phenotype(
            "CYP2D6", "Normal Metabolizer", ["celecoxib", "fluoxetine"]
        )
        # fluoxetine (strong) should dominate celecoxib (weak)
        assert result == "Poor Metabolizer"

    # -- get_all_effective_phenotypes() --

    def test_all_effective_phenotypes(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
            "CYP2C19": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        result = self.detector.get_all_effective_phenotypes(profile, ["fluoxetine"])
        assert result["CYP2D6"]["phenoconverted"] is True
        assert result["CYP2D6"]["effective_phenotype"] == "Poor Metabolizer"
        assert result["CYP2C19"]["phenoconverted"] is False

    def test_all_effective_no_meds(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        result = self.detector.get_all_effective_phenotypes(profile, [])
        assert result["CYP2D6"]["phenoconverted"] is False

    # -- clinical significance --

    def test_significance_with_substrates(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine", "codeine"])
        assert "CLINICAL ACTION REQUIRED" in alerts[0].clinical_significance

    def test_significance_without_substrates(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["fluoxetine"])
        assert "monitor" in alerts[0].clinical_significance.lower()

    # -- CYP3A5 phenoconversion --

    def test_cyp3a5_expresser_with_ketoconazole(self):
        result = self.detector.get_effective_phenotype(
            "CYP3A5", "CYP3A5 Expresser", ["ketoconazole"]
        )
        assert result == "CYP3A5 Non-expresser"

    # -- case insensitivity --

    def test_case_insensitive_medications(self):
        profile = {
            "CYP2D6": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
        }
        alerts = self.detector.detect(profile, ["FLUOXETINE"])
        assert len(alerts) >= 1


# ═══════════════════════════════════════════════════════════════════════
# NEW DATA EXPANSION TESTS
# ═══════════════════════════════════════════════════════════════════════


class TestExpandedSubstrates:
    """Tests for CYP1A2 and CYP2B6 substrate lists."""

    def test_cyp1a2_substrates_present(self):
        assert "CYP1A2" in CYP_SUBSTRATES

    def test_cyp1a2_clozapine(self):
        assert "clozapine" in CYP_SUBSTRATES["CYP1A2"]

    def test_cyp1a2_theophylline(self):
        assert "theophylline" in CYP_SUBSTRATES["CYP1A2"]

    def test_cyp2b6_substrates_present(self):
        assert "CYP2B6" in CYP_SUBSTRATES

    def test_cyp2b6_efavirenz(self):
        assert "efavirenz" in CYP_SUBSTRATES["CYP2B6"]

    def test_cyp2b6_bupropion(self):
        assert "bupropion" in CYP_SUBSTRATES["CYP2B6"]

    def test_total_enzyme_count(self):
        assert len(CYP_SUBSTRATES) >= 7


class TestExpandedInhibitors:
    """Tests for new CYP2B6 and CYP1A2 inhibitors."""

    def test_ticlopidine_cyp2b6(self):
        assert "CYP2B6" in CYP_INHIBITORS["ticlopidine"]
        assert CYP_INHIBITORS["ticlopidine"]["CYP2B6"] == InhibitorStrength.STRONG

    def test_clopidogrel_cyp2b6(self):
        assert "clopidogrel" in CYP_INHIBITORS
        assert CYP_INHIBITORS["clopidogrel"]["CYP2B6"] == InhibitorStrength.MODERATE

    def test_ritonavir_multi_enzyme(self):
        assert "ritonavir" in CYP_INHIBITORS
        assert CYP_INHIBITORS["ritonavir"]["CYP2B6"] == InhibitorStrength.STRONG
        assert CYP_INHIBITORS["ritonavir"]["CYP3A4"] == InhibitorStrength.STRONG

    def test_ciprofloxacin_cyp1a2(self):
        assert "ciprofloxacin" in CYP_INHIBITORS
        assert CYP_INHIBITORS["ciprofloxacin"]["CYP1A2"] == InhibitorStrength.STRONG

    def test_cimetidine_merged(self):
        """Cimetidine should inhibit both CYP2D6 and CYP2C19."""
        assert "CYP2D6" in CYP_INHIBITORS["cimetidine"]
        assert "CYP2C19" in CYP_INHIBITORS["cimetidine"]


class TestExpandedInducers:
    """Tests for new CYP2B6 and CYP1A2 inducers."""

    def test_rifampin_cyp2b6(self):
        assert "CYP2B6" in CYP_INDUCERS["rifampin"]
        assert CYP_INDUCERS["rifampin"]["CYP2B6"] == InducerStrength.STRONG

    def test_tobacco_smoking(self):
        assert "tobacco_smoking" in CYP_INDUCERS
        assert CYP_INDUCERS["tobacco_smoking"]["CYP1A2"] == InducerStrength.STRONG

    def test_charbroiled_food(self):
        assert "charbroiled_food" in CYP_INDUCERS


class TestCYP1A2Phenoconversion:
    """Test phenoconversion detection for CYP1A2."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.detector = PhenoconversionDetector()

    def test_fluvoxamine_inhibits_cyp1a2(self):
        eff = self.detector.get_effective_phenotype(
            "CYP1A2", "Normal Metabolizer", ["fluvoxamine"],
        )
        assert eff == "Poor Metabolizer"

    def test_ciprofloxacin_inhibits_cyp1a2(self):
        eff = self.detector.get_effective_phenotype(
            "CYP1A2", "Normal Metabolizer", ["ciprofloxacin"],
        )
        assert eff == "Poor Metabolizer"


class TestCYP2B6Phenoconversion:
    """Test phenoconversion detection for CYP2B6."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.detector = PhenoconversionDetector()

    def test_ticlopidine_inhibits_cyp2b6(self):
        eff = self.detector.get_effective_phenotype(
            "CYP2B6", "Normal Metabolizer", ["ticlopidine"],
        )
        assert eff == "Poor Metabolizer"

    def test_rifampin_induces_cyp2b6(self):
        eff = self.detector.get_effective_phenotype(
            "CYP2B6", "Normal Metabolizer", ["rifampin"],
        )
        assert eff == "Ultrarapid Metabolizer"
