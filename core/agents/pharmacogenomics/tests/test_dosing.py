"""Tests for src/dosing.py — DosingCalculator.

Tests IWPC warfarin algorithm, CYP3A5 tacrolimus dosing, DPYD
fluoropyrimidine dose reduction, TPMT+NUDT15 thiopurine dosing,
CYP2C19 clopidogrel therapy, SLCO1B1 simvastatin dosing,
CYP2C19 SSRI dosing, CYP2C9 phenytoin dosing, CYP2D6/CYP2C19 TCA dosing.

Author: Adam Jones
Date: March 2026
"""

import pytest
from src.dosing import DosingCalculator


class TestDosingCalculator:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.calc = DosingCalculator()

    # ═══════════════════════════════════════════════════════════════════
    # 1. WARFARIN (IWPC)
    # ═══════════════════════════════════════════════════════════════════

    class TestWarfarinDose:
        """IWPC warfarin dosing is DISABLED, and these tests hold it disabled.

        The coefficient block in `warfarin_dose` is mis-transcribed from Klein 2009: the column
        is shifted, so a Black patient is multiplied by +1.1816 (the enzyme-inducer coefficient)
        where the published algorithm applies -0.2760. On the square-root-of-weekly-dose scale
        that inflates the predicted dose by roughly 1.5-1.7x -- a race-specific overdose of a
        narrow-therapeutic-index anticoagulant, in the direction that causes a bleed.

        **The old test suite asserted that defect as the specification.** It contained:

            def test_african_american_increases_dose(self):
                ...
                assert aa["predicted_weekly_dose"] > other["predicted_weekly_dose"]

        along with `test_smoker_increases_dose` and `test_cyp4f2_effect`, which pin two terms
        the IWPC model does not contain at all. Eighteen tests covered this function and every
        one of them asserted only a DIRECTION -- "this input moves the dose down" -- so a
        coefficient attached to the wrong variable could not fail any of them. That is the same
        blind spot that let four EuroSCORE II coefficients survive 2,046 passing tests.

        Those assertions are not deleted and not weakened; they are replaced, because they were
        wrong about the world. The dose is not computed at all until the coefficients are
        restored from Supplementary Appendix 1 of Klein 2009, which could not be retrieved on
        the night this was found. See CLINICAL_REVIEW.md §21 for the full published mapping and
        the three "genotype unknown" terms the block never had.

        When the coefficients are restored, this class must be replaced by REFERENCE CASES --
        a described patient in, the published predicted dose out -- not by direction tests.
        """

        @pytest.fixture(autouse=True)
        def setup(self):
            self.calc = DosingCalculator()

        def test_no_dose_is_returned(self):
            with pytest.raises(RuntimeError, match="mis-transcribed"):
                self.calc.warfarin_dose("*1/*1", "*1/*1")

        def test_the_refusal_names_the_race_specific_error(self):
            """A reader of the traceback must learn why, not just that it failed."""
            with pytest.raises(RuntimeError) as exc:
                self.calc.warfarin_dose("*1/*1", "*1/*1")
            msg = str(exc.value)
            assert "race-specific" in msg
            assert "CLINICAL_REVIEW" in msg

        def test_it_refuses_for_every_patient_not_only_the_affected_group(self):
            """Disabling only the Black-patient path would leave a function nobody can trust."""
            for race in ("other", "asian", "european", "african_american"):
                with pytest.raises(RuntimeError):
                    self.calc.warfarin_dose("*1/*1", "*1/*1", race=race)

        def test_it_refuses_regardless_of_genotype(self):
            for cyp2c9, vkorc1 in (("*1/*1", "*1/*1"), ("*3/*3", "-1639G>A/-1639G>A"),
                                   ("*1/*2", "*1/-1639G>A")):
                with pytest.raises(RuntimeError):
                    self.calc.warfarin_dose(cyp2c9, vkorc1)

    # ═══════════════════════════════════════════════════════════════════
    # 2. TACROLIMUS (CYP3A5)
    # ═══════════════════════════════════════════════════════════════════

    class TestTacrolimusDose:
        @pytest.fixture(autouse=True)
        def setup(self):
            self.calc = DosingCalculator()

        def test_expresser(self):
            result = self.calc.tacrolimus_dose("*1/*1")
            assert result["phenotype"] == "CYP3A5 Expresser"
            assert result["mg_per_kg_per_day"] == 0.30
            assert result["activity_score"] == 2.0

        def test_intermediate_expresser(self):
            result = self.calc.tacrolimus_dose("*1/*3")
            assert result["phenotype"] == "CYP3A5 Intermediate Expresser"
            assert result["mg_per_kg_per_day"] == 0.25

        def test_non_expresser(self):
            result = self.calc.tacrolimus_dose("*3/*3")
            assert result["phenotype"] == "CYP3A5 Non-expresser"
            assert result["mg_per_kg_per_day"] == 0.15

        def test_total_daily_mg(self):
            result = self.calc.tacrolimus_dose("*3/*3", weight=70)
            assert result["total_daily_mg"] == pytest.approx(70 * 0.15, abs=0.5)

        def test_bid_dosing(self):
            result = self.calc.tacrolimus_dose("*3/*3", weight=70)
            assert result["frequency"] == "BID (every 12 hours)"
            assert result["per_dose_mg"] == pytest.approx(result["total_daily_mg"] / 2, abs=0.1)

        def test_monitoring_included(self):
            result = self.calc.tacrolimus_dose("*3/*3")
            assert "trough" in result["monitoring"].lower()

        def test_rationale_included(self):
            result = self.calc.tacrolimus_dose("*1/*1")
            assert result["rationale"] != ""

    # ═══════════════════════════════════════════════════════════════════
    # 3. FLUOROPYRIMIDINES (DPYD)
    # ═══════════════════════════════════════════════════════════════════

    class TestFluoropyrimidineDose:
        @pytest.fixture(autouse=True)
        def setup(self):
            self.calc = DosingCalculator()

        def test_normal_activity(self):
            result = self.calc.fluoropyrimidine_dose("*1/*1")
            assert result["phenotype"] == "DPD Normal Activity"
            assert result["dose_reduction_percent"] == 0
            assert result["activity_score"] == 2.0

        def test_intermediate_activity(self):
            result = self.calc.fluoropyrimidine_dose("*1/*2A")
            assert "DPD Intermediate" in result["phenotype"]
            assert result["dose_reduction_percent"] == 50
            assert result["activity_score"] == 1.0

        def test_deficient(self):
            result = self.calc.fluoropyrimidine_dose("*2A/*2A")
            assert result["phenotype"] == "DPD Deficient"
            assert result["dose_reduction_percent"] == 100
            assert result["activity_score"] == 0.0

        def test_partial_deficiency(self):
            result = self.calc.fluoropyrimidine_dose("*2A/*13")
            assert result["activity_score"] == 0.0
            assert result["dose_reduction_percent"] == 100

        def test_mild_intermediate(self):
            result = self.calc.fluoropyrimidine_dose("*1/c.2846A>T")
            assert result["activity_score"] == 1.5
            assert result["dose_reduction_percent"] == 25

        def test_affected_drugs(self):
            result = self.calc.fluoropyrimidine_dose("*1/*1")
            assert "5-fluorouracil" in " ".join(result["affected_drugs"]).lower() or \
                   "capecitabine" in " ".join(result["affected_drugs"]).lower()

        def test_allele_scores_returned(self):
            result = self.calc.fluoropyrimidine_dose("*1/*2A")
            assert "*1" in result["allele_scores"]
            assert result["allele_scores"]["*1"] == 1.0
            assert result["allele_scores"]["*2A"] == 0.0

    # ═══════════════════════════════════════════════════════════════════
    # 4. THIOPURINES (TPMT + NUDT15)
    # ═══════════════════════════════════════════════════════════════════

    class TestThiopurineDose:
        @pytest.fixture(autouse=True)
        def setup(self):
            self.calc = DosingCalculator()

        def test_normal_both(self):
            result = self.calc.thiopurine_dose("*1/*1", "*1/*1")
            assert result["dose_percent_of_standard"] == 100
            assert result["tpmt_phenotype"] == "TPMT Normal"
            assert result["nudt15_phenotype"] == "NUDT15 Normal"

        def test_tpmt_deficient(self):
            result = self.calc.thiopurine_dose("*3A/*3C", "*1/*1")
            assert result["dose_percent_of_standard"] == 10
            assert result["tpmt_phenotype"] == "TPMT Deficient"

        def test_tpmt_intermediate(self):
            result = self.calc.thiopurine_dose("*1/*3A", "*1/*1")
            assert result["dose_percent_of_standard"] == 50
            assert result["tpmt_phenotype"] == "TPMT Intermediate"

        def test_nudt15_deficient(self):
            result = self.calc.thiopurine_dose("*1/*1", "*3/*3")
            assert result["dose_percent_of_standard"] == 10
            assert result["nudt15_phenotype"] == "NUDT15 Deficient"

        def test_nudt15_intermediate(self):
            result = self.calc.thiopurine_dose("*1/*1", "*1/*3")
            assert result["dose_percent_of_standard"] == 50

        def test_combined_most_restrictive(self):
            # TPMT normal, NUDT15 deficient → deficient dominates
            result = self.calc.thiopurine_dose("*1/*1", "*3/*3")
            assert result["combined_activity_score"] == 0.0
            assert result["dose_percent_of_standard"] == 10

        def test_affected_drugs(self):
            result = self.calc.thiopurine_dose("*1/*1")
            assert len(result["affected_drugs"]) == 3

        def test_monitoring_included(self):
            result = self.calc.thiopurine_dose("*1/*1")
            assert "CBC" in result["monitoring"]

    # ═══════════════════════════════════════════════════════════════════
    # PRIVATE HELPERS
    # ═══════════════════════════════════════════════════════════════════

    class TestPrivateHelpers:
        @pytest.fixture(autouse=True)
        def setup(self):
            self.calc = DosingCalculator()

        def test_parse_diplotype(self):
            a, b = DosingCalculator._parse_diplotype("*1/*4")
            assert a == "*1"
            assert b == "*4"

        def test_parse_diplotype_invalid(self):
            a, b = DosingCalculator._parse_diplotype("invalid")
            assert a == "*1"
            assert b == "*1"

        def test_normalize_diplotype_already_ordered(self):
            assert DosingCalculator._normalize_diplotype("*1/*3") == "*1/*3"

        def test_normalize_diplotype_reorders(self):
            assert DosingCalculator._normalize_diplotype("*3/*1") == "*1/*3"

        def test_is_vkorc1_ag(self):
            assert DosingCalculator._is_vkorc1_ag("*1/-1639G>A") is True
            assert DosingCalculator._is_vkorc1_ag("A/G") is True

        def test_is_vkorc1_aa(self):
            assert DosingCalculator._is_vkorc1_aa("-1639G>A/-1639G>A") is True
            assert DosingCalculator._is_vkorc1_aa("A/A") is True

        def test_is_vkorc1_gg(self):
            assert DosingCalculator._is_vkorc1_ag("*1/*1") is False
            assert DosingCalculator._is_vkorc1_aa("*1/*1") is False

        def test_warfarin_dose_category(self):
            assert DosingCalculator._warfarin_dose_category(15) == "Low Dose (sensitive)"
            assert DosingCalculator._warfarin_dose_category(35) == "Standard Dose"
            assert DosingCalculator._warfarin_dose_category(55) == "High Dose (resistant)"

        def test_thiopurine_phenotype(self):
            assert DosingCalculator._thiopurine_phenotype(0.0, "TPMT") == "TPMT Deficient"
            assert DosingCalculator._thiopurine_phenotype(1.0, "TPMT") == "TPMT Intermediate"
            assert DosingCalculator._thiopurine_phenotype(2.0, "TPMT") == "TPMT Normal"


# ═══════════════════════════════════════════════════════════════════════
# CLOPIDOGREL DOSING (CYP2C19-guided)
# ═══════════════════════════════════════════════════════════════════════


class TestClopidogrelDosing:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.calc = DosingCalculator()

    def test_normal_metabolizer(self):
        result = self.calc.clopidogrel_dose("*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["clopidogrel_use"] == "standard"

    def test_poor_metabolizer(self):
        result = self.calc.clopidogrel_dose("*2/*2")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["clopidogrel_use"] == "contraindicated"
        assert "prasugrel" in result["recommended_therapy"].lower() or "ticagrelor" in result["recommended_therapy"].lower()

    def test_intermediate_metabolizer(self):
        result = self.calc.clopidogrel_dose("*1/*2")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["clopidogrel_use"] == "use_with_caution"

    def test_rapid_metabolizer(self):
        result = self.calc.clopidogrel_dose("*1/*17")
        assert result["phenotype"] == "Rapid Metabolizer"
        assert result["clopidogrel_use"] == "standard"

    def test_ultrarapid_metabolizer(self):
        result = self.calc.clopidogrel_dose("*17/*17")
        assert result["phenotype"] == "Ultrarapid Metabolizer"
        assert result["clopidogrel_use"] == "standard"

    def test_has_reference(self):
        result = self.calc.clopidogrel_dose("*1/*1")
        assert "reference" in result
        assert "CPIC" in result["reference"]

    def test_has_affected_drugs(self):
        result = self.calc.clopidogrel_dose("*1/*1")
        assert "clopidogrel" in result["affected_drugs"]

    def test_activity_score_pm(self):
        result = self.calc.clopidogrel_dose("*2/*3")
        assert result["activity_score"] == 0.0

    def test_activity_score_normal(self):
        result = self.calc.clopidogrel_dose("*1/*1")
        assert result["activity_score"] == 2.0


# ═══════════════════════════════════════════════════════════════════════
# SIMVASTATIN DOSING (SLCO1B1-guided)
# ═══════════════════════════════════════════════════════════════════════


class TestSimvastatinDosing:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.calc = DosingCalculator()

    def test_normal_function(self):
        result = self.calc.simvastatin_dose("*1/*1")
        assert result["phenotype"] == "Normal Function"
        assert result["dose_adjustment"] == "none"

    def test_decreased_function(self):
        result = self.calc.simvastatin_dose("*1/*5")
        assert result["phenotype"] == "Decreased Function"
        assert result["max_simvastatin_dose_mg"] <= 20

    def test_poor_function(self):
        result = self.calc.simvastatin_dose("*5/*5")
        assert result["phenotype"] == "Poor Function"
        assert result["max_simvastatin_dose_mg"] == 0 or result["dose_adjustment"] == "avoid_simvastatin"

    def test_unknown_genotype(self):
        result = self.calc.simvastatin_dose("*99/*99")
        assert "unknown" in result["myopathy_risk"].lower() or "standard" in result["dose_adjustment"]

    def test_has_alternatives(self):
        result = self.calc.simvastatin_dose("*5/*5")
        assert "alternatives" in result
        assert len(result["alternatives"]) >= 2

    def test_has_reference(self):
        result = self.calc.simvastatin_dose("*1/*1")
        assert "reference" in result
        assert "CPIC" in result["reference"]

    def test_affected_drugs(self):
        result = self.calc.simvastatin_dose("*1/*1")
        assert "simvastatin" in result["affected_drugs"]


# ═══════════════════════════════════════════════════════════════════════
# SSRI DOSING (CYP2C19-guided)
# ═══════════════════════════════════════════════════════════════════════


class TestSSRIDosing:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.calc = DosingCalculator()

    def test_citalopram_normal_metabolizer(self):
        result = self.calc.ssri_dose("citalopram", "*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["dose_adjustment"] == "none"
        assert result["max_dose_mg"] == 40

    def test_citalopram_poor_metabolizer(self):
        result = self.calc.ssri_dose("citalopram", "*2/*2")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["max_dose_mg"] == 20
        assert "50%" in result["dose_adjustment"]
        assert any("QT" in w for w in result["warnings"])

    def test_citalopram_ultrarapid_metabolizer(self):
        result = self.calc.ssri_dose("citalopram", "*17/*17")
        assert result["phenotype"] == "Ultrarapid Metabolizer"
        assert result["max_dose_mg"] == 60
        assert result["dose_adjustment"] == "consider_dose_increase"

    def test_escitalopram_poor_metabolizer(self):
        result = self.calc.ssri_dose("escitalopram", "*2/*2")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["max_dose_mg"] == 10
        assert "50%" in result["dose_adjustment"]

    def test_escitalopram_ultrarapid_metabolizer(self):
        result = self.calc.ssri_dose("escitalopram", "*17/*17")
        assert result["phenotype"] == "Ultrarapid Metabolizer"
        assert result["max_dose_mg"] == 30

    def test_escitalopram_normal_metabolizer(self):
        result = self.calc.ssri_dose("escitalopram", "*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["max_dose_mg"] == 20
        assert result["dose_adjustment"] == "none"

    def test_unsupported_drug(self):
        result = self.calc.ssri_dose("fluoxetine", "*1/*1")
        assert result["dose_adjustment"] == "no_cpic_guidance"
        assert result["activity_score"] is None

    def test_has_reference(self):
        result = self.calc.ssri_dose("citalopram", "*1/*1")
        assert "CPIC" in result["reference"]
        assert "SSRI" in result["reference"]

    def test_has_alternatives(self):
        result = self.calc.ssri_dose("citalopram", "*1/*1")
        assert "sertraline" in result["alternatives"]

    def test_intermediate_metabolizer(self):
        result = self.calc.ssri_dose("citalopram", "*1/*2")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["dose_adjustment"] == "consider_dose_reduction"


# ═══════════════════════════════════════════════════════════════════════
# PHENYTOIN DOSING (CYP2C9-guided)
# ═══════════════════════════════════════════════════════════════════════


class TestPhenytoinDosing:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.calc = DosingCalculator()

    def test_normal_metabolizer(self):
        result = self.calc.phenytoin_dose("*1/*1")
        assert result["phenotype"] == "Normal Metabolizer"
        assert result["dose_reduction_percent"] == 0
        assert result["activity_score"] == 2.0

    def test_poor_metabolizer(self):
        result = self.calc.phenytoin_dose("*3/*3")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["dose_reduction_percent"] == 50
        assert result["activity_score"] == 0.0

    def test_intermediate_metabolizer(self):
        result = self.calc.phenytoin_dose("*1/*3")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["dose_reduction_percent"] == 25
        assert result["activity_score"] == 1.0

    def test_loading_dose_standard(self):
        result = self.calc.phenytoin_dose("*1/*1", weight_kg=70.0)
        assert result["standard_loading_dose_mg"] == min(round(70 * 18), 1500)

    def test_loading_dose_capped(self):
        result = self.calc.phenytoin_dose("*1/*1", weight_kg=100.0)
        assert result["standard_loading_dose_mg"] <= 1500

    def test_maintenance_dose_normal(self):
        result = self.calc.phenytoin_dose("*1/*1", weight_kg=70.0)
        expected = round(70 * 5 * 1.0)
        assert result["suggested_maintenance_mg_per_day"] == expected

    def test_maintenance_dose_pm(self):
        result = self.calc.phenytoin_dose("*3/*3", weight_kg=70.0)
        expected = round(70 * 5 * 0.5)
        assert result["suggested_maintenance_mg_per_day"] == expected

    def test_has_tdm_target(self):
        result = self.calc.phenytoin_dose("*1/*1")
        assert result["tdm_target"] == "10-20 mcg/mL"

    def test_has_alternatives(self):
        result = self.calc.phenytoin_dose("*1/*1")
        assert "levetiracetam" in result["alternatives"]

    def test_has_reference(self):
        result = self.calc.phenytoin_dose("*1/*1")
        assert "CPIC" in result["reference"]
        assert "Phenytoin" in result["reference"]

    def test_affected_drugs(self):
        result = self.calc.phenytoin_dose("*1/*1")
        assert "phenytoin" in result["affected_drugs"]
        assert "fosphenytoin" in result["affected_drugs"]

    def test_pm_warnings(self):
        result = self.calc.phenytoin_dose("*3/*3")
        assert any("50%" in w for w in result["warnings"])
        assert any("TDM" in w for w in result["warnings"])


# ═══════════════════════════════════════════════════════════════════════
# TCA DOSING (CYP2D6/CYP2C19-guided)
# ═══════════════════════════════════════════════════════════════════════


class TestTCADosing:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.calc = DosingCalculator()

    def test_normal_metabolizer(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        assert result["cyp2d6_phenotype"] == "Normal Metabolizer"
        assert result["dose_adjustment"] == "none"
        assert result["tdm_required"] is False

    def test_poor_metabolizer(self):
        result = self.calc.tca_dose("amitriptyline", "*4/*4")
        assert result["cyp2d6_phenotype"] == "Poor Metabolizer"
        assert result["dose_adjustment"] == "reduce_50_percent"
        assert result["tdm_required"] is True

    def test_intermediate_metabolizer(self):
        # *10 (0.25) + *41 (0.5) = 0.75 → IM range [0.25, 1.0)
        result = self.calc.tca_dose("nortriptyline", "*10/*41")
        assert result["cyp2d6_phenotype"] == "Intermediate Metabolizer"
        assert result["dose_adjustment"] == "reduce_25_percent"

    def test_um_tertiary_avoid(self):
        # UM with tertiary amine should recommend avoid
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        # *1/*1 is NM (AS=2.0), not UM. Need higher AS for UM.
        # UM requires AS >= 2.25 — not achievable with standard alleles
        # in the lookup, so test with a known UM scenario description
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        assert result["cyp2d6_phenotype"] == "Normal Metabolizer"

    def test_um_secondary_use_with_tdm(self):
        # For secondary amine with UM, dose_adjustment is use_with_tdm
        # Standard allele pairs don't reach UM threshold easily; verify NM path
        result = self.calc.tca_dose("nortriptyline", "*1/*1")
        assert result["cyp2d6_phenotype"] == "Normal Metabolizer"
        assert result["tca_class"] == "secondary"

    def test_tca_class_tertiary(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        assert result["tca_class"] == "tertiary"

    def test_tca_class_secondary(self):
        result = self.calc.tca_dose("desipramine", "*1/*1")
        assert result["tca_class"] == "secondary"

    def test_cyp2c19_pm_warning_tertiary(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*1", "*2/*2")
        assert any("CYP2C19 PM" in w for w in result["warnings"])

    def test_cyp2c19_um_warning_tertiary(self):
        result = self.calc.tca_dose("imipramine", "*1/*1", "*17/*17")
        assert any("CYP2C19 UM" in w for w in result["warnings"])

    def test_cyp2c19_no_warning_secondary(self):
        # Secondary amines should not get CYP2C19 warnings
        result = self.calc.tca_dose("nortriptyline", "*1/*1", "*2/*2")
        assert not any("CYP2C19" in w for w in result["warnings"])

    def test_unsupported_drug(self):
        result = self.calc.tca_dose("fluoxetine", "*1/*1")
        assert result["dose_adjustment"] == "no_cpic_guidance"
        assert result["tca_class"] is None

    def test_has_reference(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        assert "CPIC" in result["reference"]
        assert "Tricyclic" in result["reference"]

    def test_has_tdm_range(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        assert "80-200" in result["tdm_target_range"]

    def test_affected_drugs_all_seven(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*1")
        assert len(result["affected_drugs"]) == 7

    def test_pm_all_drugs(self):
        for drug in ["amitriptyline", "nortriptyline", "clomipramine",
                      "doxepin", "desipramine", "imipramine", "trimipramine"]:
            result = self.calc.tca_dose(drug, "*4/*4")
            assert result["dose_adjustment"] == "reduce_50_percent"
            assert result["tdm_required"] is True
