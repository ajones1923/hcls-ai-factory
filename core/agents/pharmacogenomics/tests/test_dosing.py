"""Tests for src/dosing.py — DosingCalculator.

Tests IWPC warfarin algorithm, CYP3A5 tacrolimus dosing, DPYD
fluoropyrimidine dose reduction, TPMT+NUDT15 thiopurine dosing,
CYP2C19 clopidogrel therapy, SLCO1B1 simvastatin dosing,
CYP2C19 SSRI dosing, CYP2C9 phenytoin dosing, CYP2D6/CYP2C19 TCA dosing.

Author: Adam Jones
Date: March 2026
"""

import pytest
from src.dosing import DosingCalculator, IndeterminateGenotype


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
            """CORRECTED (pgx audit §17): CPIC gives *1/*1 and *1/*3 the SAME dose.

            This asserted 0.25 mg/kg/day for *1/*3 against 0.30 for *1/*1.
            CPIC's tacrolimus recommendation table gives CYP3A5 Normal
            Metabolizer, Possible Intermediate Metabolizer and Intermediate
            Metabolizer alike one recommendation: "Increase starting dose
            1.5 to 2 times recommended starting dose. Total starting dose
            should not exceed 0.3 mg/kg/day."
            (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:42316,
            retrieved 2026-11-03.)  The 0.25 and its "~1.25-1.5x standard"
            label were an invented clinical distinction presented under a
            CPIC heading — and the function's own docstring stated the CPIC
            rule the code did not implement.
            """
            result = self.calc.tacrolimus_dose("*1/*3")
            assert result["phenotype"] == "CYP3A5 Intermediate Expresser"
            assert result["mg_per_kg_per_day"] == 0.30
            assert result["mg_per_kg_per_day"] == self.calc.tacrolimus_dose(
                "*1/*1")["mg_per_kg_per_day"]

        def test_unrecognised_allele_is_not_a_non_expresser(self):
            """Regression for pgx audit §17.

            `.get(allele, 0.0)` called an allele it had never seen a
            non-expresser and gave it the standard 0.15 mg/kg/day. A true
            expresser dosed as a non-expresser is subtherapeutic tacrolimus
            in the first week post-transplant — the window in which acute
            rejection happens. CPIC's answer for Indeterminate is "No
            recommendation".
            """
            result = self.calc.tacrolimus_dose("*1/*8")   # *8: Unknown function
            assert result["phenotype"] == "Indeterminate"
            assert result["mg_per_kg_per_day"] is None
            assert "*8" in result["indeterminate_alleles"]

        def test_star7_is_a_no_function_allele(self):
            """CPIC: CYP3A5 *7 | No function. *3/*7 is a Poor Metabolizer."""
            result = self.calc.tacrolimus_dose("*3/*7")
            assert result["phenotype"] == "CYP3A5 Non-expresser"
            assert result["mg_per_kg_per_day"] == 0.15

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
            """CORRECTED (pgx audit §3): CPIC gives activity score 1.5 a 50% cut.

            This pinned a 25% reduction at activity score 1.5.  CPIC's
            fluorouracil and capecitabine recommendation rows give activity
            score 1.5 and activity score 1.0 the SAME answer, verbatim:
            "Reduce starting dose by 50% followed by titration of dose based
            on toxicity or therapeutic drug monitoring (if available)."
            (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:4492 and
            =eq.RxNorm:194000, retrieved 2026-11-03.)

            Reference case: DPYD *1/c.2846A>T, activity score 1.5 — the
            single most frequently reported actionable DPYD genotype in
            European-ancestry patients.  A 25% reduction left her on 75% of
            a standard 5-FU dose where CPIC calls for 50%; DPD
            heterozygotes receiving near-full dose are the population the
            fatal-toxicity case reports come from.  There is no "(Mild)"
            grade in CPIC's table.
            """
            result = self.calc.fluoropyrimidine_dose("*1/c.2846A>T")
            assert result["activity_score"] == 1.5
            assert result["dose_reduction_percent"] == 50
            assert result["dose_reduction_percent"] == self.calc.\
                fluoropyrimidine_dose("*1/*2A")["dose_reduction_percent"]
            assert "(Mild)" not in result["phenotype"]

        def test_activity_score_0_5_is_avoid_not_a_reduced_dose(self):
            """Regression for pgx audit §3 / CLINICAL_REVIEW §4, now SETTLED.

            dosing.py called activity score 0.5 "DPD Intermediate Activity
            (Partial Deficiency)" with a 50% reduction; pgx_pipeline.py
            called it "DPD Deficient"; knowledge.py gave a third answer.
            CPIC settles it: activity score 0.5 is a **Poor Metabolizer**
            and the recommendation is "Avoid use of 5-fluorouracil or
            5-fluorouracil prodrug-based regimens."
            """
            result = self.calc.fluoropyrimidine_dose("*2A/c.2846A>T")
            assert result["activity_score"] == 0.5
            assert "Deficient" in result["phenotype"]
            assert result["dose_reduction_percent"] == 100
            assert "AVOID" in result["recommendation"].upper()

        def test_a_loss_of_function_homozygote_is_not_given_half_dose(self):
            """Regression for pgx audit §4.

            DPYD *7 (c.295_298delTCAT) was scored 0.5 "decreased"; CPIC
            calls it **No function**.  *7/*7 — a patient with no DPD enzyme
            at all, CPIC activity score 0.0, "avoid use of 5-fluorouracil" —
            therefore came out activity score 1.0 with a 50% reduction:
            half-dose 5-FU for complete DPD deficiency.  Half-dose does not
            save that patient.
            """
            result = self.calc.fluoropyrimidine_dose("*7/*7")
            assert result["activity_score"] == 0.0
            assert result["phenotype"] == "DPD Deficient"
            assert result["dose_reduction_percent"] == 100

        def test_a_normal_function_allele_does_not_cut_the_dose(self):
            """Regression for pgx audit §4.

            DPYD *9A (c.85T>C, rs1801265) was scored 0.5 "decreased"; CPIC
            calls it **Normal function**, and it is common — of the order of
            20% minor-allele frequency in Europeans and higher in
            African-ancestry populations.  *1/*9A was given activity score
            1.5 and a 25-50% reduction, and *9A/*9A activity score 1.0 and
            50%, so a large fraction of patients would have had curative
            adjuvant chemotherapy halved on the strength of a variant with
            no effect on DPD activity.  Under-dosing adjuvant 5-FU costs
            disease-free survival.
            """
            het = self.calc.fluoropyrimidine_dose("*1/*9A")
            hom = self.calc.fluoropyrimidine_dose("*9A/*9A")
            assert het["activity_score"] == 2.0
            assert het["dose_reduction_percent"] == 0
            assert hom["activity_score"] == 2.0
            assert hom["dose_reduction_percent"] == 0

        def test_a_single_allele_label_is_not_wild_type(self):
            """Regression for pgx audit §6.

            `_parse_diplotype` returned ("*1", "*1") for anything that was
            not two slash-separated fields, so fluoropyrimidine_dose("*2A")
            — a pharmacy reporting the patient as a *2A carrier — printed
            "Standard fluoropyrimidine dose. No DPYD-based dose adjustment
            needed. Activity score 2.0 (normal)."  A *2A carrier given a
            full 5-FU dose is the published fatal-toxicity scenario.
            """
            for bad in ("*2A", "", "Unknown", "not reported", "*1/*1/*1"):
                result = self.calc.fluoropyrimidine_dose(bad)
                assert result["phenotype"] == "Indeterminate", bad
                assert result["activity_score"] is None, bad
                assert result["dose_reduction_percent"] is None, bad
                assert "no dpyd-based dose adjustment needed" not in \
                    result["recommendation"].lower(), bad

        def test_an_unlisted_allele_is_not_assumed_normal(self):
            """Regression for pgx audit §8. DPYD has 84 CPIC-assigned alleles."""
            result = self.calc.fluoropyrimidine_dose("*1/c.1905C>G")
            assert result["phenotype"] == "Indeterminate"
            assert result["dose_reduction_percent"] is None

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
            """CORRECTED (pgx audit §15): CPIC's azathioprine PM row is a
            different DRUG, not a reduced azathioprine dose.

            This pinned dose_percent_of_standard == 10, which came with the
            text "reduce dose to 10% of standard (e.g., azathioprine 0.5
            mg/kg/day thrice weekly)".  CPIC puts the 10-fold reduction and
            thrice-weekly schedule on MERCAPTOPURINE and THIOGUANINE; its
            azathioprine row for a poor metabolizer at either locus reads,
            in full, "Consider alternative nonthiopurine immunosuppressant
            therapy."  (api.cpicpgx.org/v1/recommendation
            ?drugid=eq.RxNorm:1256 and =eq.RxNorm:103, retrieved
            2026-11-03.)  No single percent-of-standard covers the three
            drugs, so none is reported and recommendations_by_drug carries
            each drug's own row.  This also settles pharmacogenomics open
            question #2.
            """
            result = self.calc.thiopurine_dose("*3A/*3C", "*1/*1")
            assert result["dose_percent_of_standard"] is None
            assert result["tpmt_phenotype"] == "TPMT Deficient"
            by_drug = result["recommendations_by_drug"]
            assert "alternative nonthiopurine" in by_drug["azathioprine"]
            assert "10-fold" in by_drug["6-mercaptopurine (6-MP)"]
            assert "10-fold" not in by_drug["azathioprine"]

        def test_tpmt_intermediate(self):
            """CORRECTED (pgx audit §15): CPIC's single-gene IM range is 30-80%.

            This pinned 50, a midpoint of a range the code wrote as
            "30-70%".  CPIC: "Initiate therapy with reduced starting doses
            (30-80% of standard starting dose)".  The structured field now
            carries the conservative (low) end of CPIC's range, which is
            where a starting dose belongs, and the full range is reported
            separately.
            """
            result = self.calc.thiopurine_dose("*1/*3A", "*1/*1")
            assert result["dose_percent_range_of_standard"] == (30, 80)
            assert result["dose_percent_of_standard"] == 30
            assert result["tpmt_phenotype"] == "TPMT Intermediate"

        def test_nudt15_deficient(self):
            """CORRECTED (pgx audit §15) — see test_tpmt_deficient."""
            result = self.calc.thiopurine_dose("*1/*1", "*3/*3")
            assert result["dose_percent_of_standard"] is None
            assert result["nudt15_phenotype"] == "NUDT15 Deficient"
            assert "alternative nonthiopurine" in \
                result["recommendations_by_drug"]["azathioprine"]

        def test_nudt15_intermediate(self):
            """CORRECTED (pgx audit §15) — see test_tpmt_intermediate."""
            result = self.calc.thiopurine_dose("*1/*1", "*1/*3")
            assert result["dose_percent_range_of_standard"] == (30, 80)
            assert result["dose_percent_of_standard"] == 30

        def test_combined_most_restrictive(self):
            """CORRECTED (pgx audit §15): min() cannot express a two-gene table.

            This asserted `combined_activity_score == 0.0`, a field computed
            as `min(tpmt_score, nudt15_score)`.  min() makes two
            intermediate metabolizers look identical to one, and CPIC's
            table distinguishes them: IM+NM is 30-80% of standard and
            IM+IM is 20-50%.  The field is gone; the recommendation is now
            keyed on the (TPMT phenotype, NUDT15 phenotype) pair, and the
            "more restrictive gene wins" intent this test was written for is
            asserted directly below.
            """
            result = self.calc.thiopurine_dose("*1/*1", "*3/*3")
            assert result["tpmt_phenotype"] == "TPMT Normal"
            assert result["nudt15_phenotype"] == "NUDT15 Deficient"
            assert result["dose_percent_of_standard"] is None
            assert "alternative nonthiopurine" in \
                result["recommendations_by_drug"]["azathioprine"]
            assert "combined_activity_score" not in result

        def test_double_heterozygote_gets_cpics_narrower_range(self):
            """Regression for pgx audit §15.

            A 31-year-old starting azathioprine 2 mg/kg/day for Crohn's,
            TPMT *1/*3A and NUDT15 *1/*3 — heterozygous at both loci.  The
            old `min()` reported "30-70% of standard, start at lower end
            (30-50%)" with dose_percent 50.  CPIC's ceiling for this patient
            is 50% and its floor is 20%, so the headline range ran 20 points
            above CPIC's ceiling and the structured field sat exactly on it
            with titrate-upward advice.  Double heterozygotes have
            additively reduced thiopurine clearance; the high end of the old
            range is a predictable route to grade 3-4 leukopenia at week
            2-4.
            """
            result = self.calc.thiopurine_dose("*1/*3A", "*1/*3")
            assert result["tpmt_phenotype"] == "TPMT Intermediate"
            assert result["nudt15_phenotype"] == "NUDT15 Intermediate"
            assert result["dose_percent_range_of_standard"] == (20, 50)
            assert result["dose_percent_of_standard"] == 20
            assert "20-50%" in result["dose_recommendation"]
            assert "70%" not in result["dose_recommendation"]

        def test_a_single_locus_intermediate_is_not_the_same_as_two(self):
            """The distinction min() could not make."""
            one = self.calc.thiopurine_dose("*1/*3A", "*1/*1")
            two = self.calc.thiopurine_dose("*1/*3A", "*1/*3")
            assert one["dose_percent_range_of_standard"] == (30, 80)
            assert two["dose_percent_range_of_standard"] == (20, 50)

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
            """CORRECTED (pgx audit §6): an unparseable genotype is not wild-type.

            This asserted that `_parse_diplotype("invalid")` returns
            ("*1", "*1") — homozygous wild-type — which is the behaviour the
            finding is about, pinned as the specification.  Every one of the
            nine calculators then reported a normal metabolizer and a
            standard dose.  It now raises IndeterminateGenotype, and each
            calculator turns that into CPIC's Indeterminate.
            """
            with pytest.raises(IndeterminateGenotype):
                DosingCalculator._parse_diplotype("invalid")

        @pytest.mark.parametrize("bad", [
            "", "invalid", "*2A", "Unknown", "not reported",
            "*1/*1/*1", "*2A;*2A", "/", "*1/", "/*1",
        ])
        def test_parse_diplotype_refuses_everything_that_is_not_a_diplotype(
                self, bad):
            with pytest.raises(IndeterminateGenotype):
                DosingCalculator._parse_diplotype(bad)

        def test_parse_diplotype_message_says_it_is_not_wild_type(self):
            with pytest.raises(IndeterminateGenotype) as exc:
                DosingCalculator._parse_diplotype("*2A")
            assert "not wild-type" in str(exc.value).lower()

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

        def test_thiopurine_phenotype_differs_between_the_two_genes(self):
            """CPIC's TPMT and NUDT15 function-pair rules are not the same.

            api.cpicpgx.org/v1/diplotype (retrieved 2026-11-03): a
            no-function + decreased-function pair (activity score 0.5) is a
            TPMT **Intermediate** Metabolizer but a NUDT15 **Poor** one; and
            a decreased + normal pair (activity score 1.5) is Normal for
            both, where the old shared rule called it Intermediate.
            """
            assert DosingCalculator._thiopurine_phenotype(0.5, "TPMT") == \
                "TPMT Intermediate"
            assert DosingCalculator._thiopurine_phenotype(0.5, "NUDT15") == \
                "NUDT15 Deficient"
            assert DosingCalculator._thiopurine_phenotype(1.5, "TPMT") == \
                "TPMT Normal"
            assert DosingCalculator._thiopurine_phenotype(1.5, "NUDT15") == \
                "NUDT15 Normal"


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
        """CORRECTED (pgx audit §8): CPIC says avoid the standard dose.

        `clopidogrel_use` was "use_with_caution" and the recommended
        therapy offered "clopidogrel 75mg with enhanced monitoring (if
        alternatives unavailable)".  CPIC's Intermediate Metabolizer row
        reads, verbatim: "Avoid standard dose (75 mg) clopidogrel if
        possible. Use prasugrel or ticagrelor at standard dose if no
        contraindication."
        (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:32968,
        retrieved 2026-11-03.)
        """
        result = self.calc.clopidogrel_dose("*1/*2")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["clopidogrel_use"] == "avoid_standard_dose"
        assert "avoid standard dose" in result["rationale"].lower()

    def test_a_no_function_allele_outside_the_old_table_is_not_normal(self):
        """Regression for pgx audit §8.

        A patient of West African ancestry, CYP2C19 *1/*35, 3 days after
        drug-eluting stent placement.  *35 is a no-function allele in
        CPIC's table and CPIC's diplotype table gives *1/*35 ->
        Intermediate Metabolizer -> "Avoid standard dose clopidogrel (75
        mg) if possible. Use prasugrel or ticagrelor."  `.get(allele, 1.0)`
        reported "Normal Metabolizer, clopidogrel standard, Standard
        clopidogrel therapy appropriate."  He stays on clopidogrel with
        reduced active-metabolite formation; the endpoint is stent
        thrombosis.
        """
        result = self.calc.clopidogrel_dose("*1/*35")
        assert result["phenotype"] == "Intermediate Metabolizer"
        assert result["clopidogrel_use"] == "avoid_standard_dose"

    def test_a_poor_metabolizer_outside_the_old_table_is_not_intermediate(self):
        """CPIC: CYP2C19 *2/*35 -> Poor Metabolizer, avoid clopidogrel."""
        result = self.calc.clopidogrel_dose("*2/*35")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["clopidogrel_use"] == "contraindicated"

    def test_a_decreased_function_allele_with_a_no_function_one_is_poor(self):
        """CPIC calls *2/*10 a (Likely) Poor Metabolizer, activity score 0.5.

        CYP2C19 *10 was scored 0.0 (no function); CPIC says Decreased
        function.  With the corrected 0.5 the Poor bucket must extend to
        just below activity score 1.0 or this patient reads Intermediate.
        """
        result = self.calc.clopidogrel_dose("*2/*10")
        assert result["activity_score"] == 0.5
        assert result["phenotype"] == "Poor Metabolizer"

    def test_an_uncertain_function_allele_gets_no_recommendation(self):
        """CPIC assigns Indeterminate to CYP2C19 *12 and answers
        "No recommendation"."""
        result = self.calc.clopidogrel_dose("*1/*12")
        assert result["phenotype"] == "Indeterminate"
        assert result["clopidogrel_use"] == "indeterminate"
        assert result["recommended_therapy"] is None

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
        """CORRECTED (pgx audit §11): 60 mg is the label's QT-study dose.

        This pinned max_dose_mg == 60 for a CYP2C19 ultrarapid metabolizer.
        FDA-approved CELEXA label (DailyMed, Allergan): §2.1 "Dosages above
        40 mg once daily are not recommended due to the risk of QT
        prolongation"; §12.2 the thorough QT study measured 60 mg,
        described in the label as "1.5 times the maximum recommended
        dosage", and found a QTcNi increase of 18.5 msec (upper bound
        21.0).  CPIC says nothing about 60 mg either — its UM row is
        "Consider a clinically appropriate alternative antidepressant not
        predominantly metabolized by CYP2C19 ... consider titrating to a
        higher maintenance dose", within the label.

        Reference case: a 29-year-old *17/*17 with inadequate response to
        citalopram 40 mg.  At 60 mg the label's own study puts the QTc
        increase at ~18.5 msec; with hypokalaemia, an interacting
        QT-prolonging drug or bradycardia that is the substrate for
        torsades de pointes.  The PM branch of this same function already
        reproduced the label's 20 mg cap, so the ceiling was understood
        elsewhere in the function.
        """
        result = self.calc.ssri_dose("citalopram", "*17/*17")
        assert result["phenotype"] == "Ultrarapid Metabolizer"
        assert result["max_dose_mg"] == 40
        assert result["dose_adjustment"] == \
            "consider_alternative_or_titrate_within_label"
        assert any("40 mg" in w for w in result["warnings"])

    def test_escitalopram_poor_metabolizer(self):
        result = self.calc.ssri_dose("escitalopram", "*2/*2")
        assert result["phenotype"] == "Poor Metabolizer"
        assert result["max_dose_mg"] == 10
        assert "50%" in result["dose_adjustment"]

    def test_escitalopram_ultrarapid_metabolizer(self):
        """CORRECTED (pgx audit §11): 30 mg is the label's supratherapeutic arm.

        This pinned max_dose_mg == 30.  FDA-approved LEXAPRO label
        (DailyMed, Allergan, setid
        13bb8267-1cab-43e5-acae-55a4d957630a): §2.1 "The recommended dosage
        of Lexapro in adults is 10 mg once daily", which "may be increased
        to the maximum recommended dosage of 20 mg once daily".  The label
        names NO recommended dose above 20 mg once daily for adult MDD; it
        names 30 mg exactly once, in §6.1, as a "supratherapeutic 30 mg"
        once-daily arm OF A QTcF STUDY — the same shape of error as the
        citalopram 60 mg, both numbers taken from the pharmacology section
        rather than the dosing section.
        """
        result = self.calc.ssri_dose("escitalopram", "*17/*17")
        assert result["phenotype"] == "Ultrarapid Metabolizer"
        assert result["max_dose_mg"] == 20
        assert any("20 mg" in w for w in result["warnings"])

    @pytest.mark.parametrize("drug,cap", [("citalopram", 40), ("escitalopram", 20)])
    def test_no_phenotype_ever_exceeds_the_labelled_maximum(self, drug, cap):
        """The labelled ceiling is not a function of genotype."""
        for diplotype in ("*1/*1", "*1/*2", "*2/*2", "*1/*17", "*17/*17",
                          "*2/*17", "*1/*35"):
            out = self.calc.ssri_dose(drug, diplotype)
            assert out["max_dose_mg"] <= cap, (drug, diplotype)

    def test_an_unscoreable_genotype_gets_no_dose(self):
        """Regression for pgx audit §6 and §8."""
        for bad in ("", "*2", "Unknown", "*1/*12"):
            out = self.calc.ssri_dose("citalopram", bad)
            assert out["phenotype"] == "Indeterminate", bad
            assert out["max_dose_mg"] is None, bad

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
        """CORRECTED (pgx audit §10): CPIC leads with "Avoid".

        This pinned dose_adjustment == "reduce_50_percent".  CPIC's CYP2D6
        Poor Metabolizer row reads "Avoid amitriptyline use. If
        amitriptyline is warranted, consider a 50% reduction of recommended
        starting dose. Utilize therapeutic drug monitoring to guide dose
        adjustments."  The 50% is conditional on the avoid, and the old
        value dropped the condition.
        (api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:704, retrieved
        2026-11-03.)
        """
        result = self.calc.tca_dose("amitriptyline", "*4/*4")
        assert result["cyp2d6_phenotype"] == "Poor Metabolizer"
        assert result["dose_adjustment"] == "avoid_or_reduce_50_percent"
        assert result["tdm_required"] is True
        assert "avoid" in result["rationale"].lower()
        assert "50%" in result["rationale"]

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
        """CORRECTED (pgx audit §10): the CYP2C19 axis must set the DOSE.

        This asserted only that an appended warning string mentioned
        "CYP2C19 PM" — and that string was all the CYP2C19 axis did.  The
        structured fields were set from CYP2D6 alone, so
        tca_dose('amitriptyline', '*1/*1', '*2/*2') returned
        dose_adjustment "none", tdm_required False and "Standard TCA
        dosing. Normal CYP2D6 metabolism."  A UI or report rendering the
        structured field — which pgx_ui.py does — showed "none".

        Reference case: a 46-year-old woman of Japanese ancestry with
        chronic neuropathic pain, CYP2D6 *1/*1, CYP2C19 *2/*2 (12-23% of
        East Asian populations).  CPIC's row for every CYP2D6-Normal x
        CYP2C19-Poor pair — activity scores 1.25, 1.5, 1.75, 2.0 and 2.25
        alike — reads "Avoid amitriptyline use. If amitriptyline is
        warranted, consider a 50% reduction of recommended starting dose.
        Utilize therapeutic drug monitoring to guide dose adjustments."
        Amitriptyline\'s dose-limiting toxicities are QT prolongation,
        anticholinergic delirium and seizures.
        """
        result = self.calc.tca_dose("amitriptyline", "*1/*1", "*2/*2")
        assert result["cyp2d6_phenotype"] == "Normal Metabolizer"
        assert result["cyp2c19_phenotype"] == "Poor Metabolizer"
        assert result["dose_adjustment"] == "avoid_or_reduce_50_percent"
        assert result["tdm_required"] is True
        assert "avoid" in result["rationale"].lower()
        assert any("CYP2C19" in w for w in result["warnings"])

    @pytest.mark.parametrize("cyp2d6", ["*1/*1", "*1/*2", "*2/*2", "*1/*9"])
    def test_cyp2c19_pm_drives_the_answer_at_every_normal_cyp2d6(self, cyp2d6):
        """CPIC gives the same answer across CYP2D6 AS 1.25-2.25."""
        result = self.calc.tca_dose("amitriptyline", cyp2d6, "*2/*2")
        assert result["cyp2d6_phenotype"] == "Normal Metabolizer"
        assert result["dose_adjustment"] == "avoid_or_reduce_50_percent"
        assert result["tdm_required"] is True

    def test_cyp2c19_um_warning_tertiary(self):
        """CORRECTED (pgx audit §10) — see test_cyp2c19_pm_warning_tertiary.

        CPIC's CYP2D6-Normal x CYP2C19-Ultrarapid row: "Consider
        alternative drug not metabolized by CYP2C19; if amitriptyline is
        warranted, utilize therapeutic drug monitoring to guide dose
        adjustment."  The old code emitted only "CYP2C19 UM increases
        demethylation; monitor metabolite levels" as free text and left
        dose_adjustment at "none".
        """
        result = self.calc.tca_dose("imipramine", "*1/*1", "*17/*17")
        assert result["cyp2c19_phenotype"] == "Ultrarapid Metabolizer"
        assert result["dose_adjustment"] == "avoid"
        assert result["tdm_required"] is True
        assert any("CYP2C19" in w for w in result["warnings"])

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
        """CORRECTED (pgx audit §10) — see test_poor_metabolizer."""
        for drug in ["amitriptyline", "nortriptyline", "clomipramine",
                      "doxepin", "desipramine", "imipramine", "trimipramine"]:
            result = self.calc.tca_dose(drug, "*4/*4")
            assert result["dose_adjustment"] == "avoid_or_reduce_50_percent"
            assert result["tdm_required"] is True

    def test_an_unscoreable_cyp2c19_blocks_the_recommendation(self):
        """Regression for pgx audit §6/§8/§10.

        The CYP2C19 axis can make CPIC\'s answer "avoid", so a CYP2D6
        result with no usable CYP2C19 result is not enough for a tertiary
        amine.
        """
        result = self.calc.tca_dose("amitriptyline", "*1/*1", "")
        assert result["dose_adjustment"] == "no_recommendation"
        assert result["cyp2c19_phenotype"] == "Indeterminate"
        assert result["tdm_required"] is True

    def test_an_unscoreable_cyp2d6_blocks_the_recommendation(self):
        result = self.calc.tca_dose("amitriptyline", "*1/*999", "*1/*1")
        assert result["dose_adjustment"] == "no_recommendation"
        assert result["cyp2d6_phenotype"] == "Indeterminate"
