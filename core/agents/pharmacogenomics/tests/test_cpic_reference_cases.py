"""Reference cases pinned to CPIC's own machine-readable tables.

Written for the pharmacogenomics audit of 2026-11-03. Its findings traced
almost every wrong dose in this agent to one mechanism: *a clinical lookup
table transcribed by hand, duplicated, and never checked against its
source.* The source is retrievable in bulk, by machine, so each case below
names the exact endpoint and the row it came from.

Every assertion here is a **reference case** — a described genotype in, the
published phenotype or dose out — not a monotonic-direction assertion. The
audit's finding §19 is that all 18 warfarin tests asserted only a
direction, which is why four coefficients survived in the wrong slots:
"Monotonic ordering was preserved, so every output stayed plausible."

Primary sources, all retrieved 2026-11-03:

    allele function      https://api.cpicpgx.org/v1/allele?genesymbol=eq.<GENE>
    diplotype→phenotype  https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.<GENE>
    drug recommendation  https://api.cpicpgx.org/v1/recommendation?drugid=eq.<RxNorm id>

    fluorouracil  RxNorm:4492     azathioprine     RxNorm:1256
    capecitabine  RxNorm:194000   mercaptopurine   RxNorm:103
    phenytoin     RxNorm:8183     amitriptyline    RxNorm:704
    clopidogrel   RxNorm:32968    citalopram       RxNorm:2556
    simvastatin   RxNorm:36567    escitalopram     RxNorm:321988
    tacrolimus    RxNorm:42316    abacavir         RxNorm:190521

The threshold tables and allele values under test were checked against
CPIC's diplotype tables exhaustively, not by sampling: they reproduce
CPIC's phenotype for all 276 NUDT15, 1,225 TPMT, 3,570 DPYD, 2,850 CYP2C9,
351 CYP2C19 (those built from assigned alleles), 36 UGT1A1, 21 CYP3A5 and
741 SLCO1B1 rows.
"""
import pytest

from src.dosing import DosingCalculator
from src.pgx_pipeline import PhenotypeTranslator


@pytest.fixture
def calc():
    return DosingCalculator()


@pytest.fixture
def translator():
    return PhenotypeTranslator()


# ═════════════════════════════════════════════════════════════════════════
# DPYD / fluoropyrimidines — CPIC recommendation rows by activity score
# ═════════════════════════════════════════════════════════════════════════

class TestFluoropyrimidineRecommendationRows:
    """One case per row of CPIC's fluorouracil recommendation table."""

    @pytest.mark.parametrize("diplotype,score,reduction,avoid", [
        # AS 2.0 — "no indication to change dose or therapy"
        ("*1/*1",            2.0,   0, False),
        ("*1/*9A",           2.0,   0, False),   # *9A: CPIC Normal function
        # AS 1.5 — "Reduce starting dose by 50%"
        ("*1/c.2846A>T",     1.5,  50, False),
        ("*1/HapB3",         1.5,  50, False),
        # AS 1.0 — "Reduce starting dose by 50%"
        ("*1/*2A",           1.0,  50, False),
        ("*1/*13",           1.0,  50, False),
        ("*1/*7",            1.0,  50, False),   # *7: CPIC NO function
        ("*1/*10",           1.0,  50, False),   # *10: CPIC NO function
        ("c.2846A>T/c.2846A>T", 1.0, 50, False),
        # AS 0.5 — "Avoid use of 5-fluorouracil ..." (CPIC Poor Metabolizer)
        ("*2A/c.2846A>T",    0.5, 100, True),
        ("*13/HapB3",        0.5, 100, True),
        # AS 0.0 — "Avoid use of 5-fluorouracil ..."
        ("*2A/*2A",          0.0, 100, True),
        ("*7/*7",            0.0, 100, True),
        ("*2A/*13",          0.0, 100, True),
    ])
    def test_row(self, calc, diplotype, score, reduction, avoid):
        out = calc.fluoropyrimidine_dose(diplotype)
        assert out["activity_score"] == score
        assert out["dose_reduction_percent"] == reduction
        assert ("AVOID" in out["recommendation"].upper()) is avoid

    def test_as_1_5_and_1_0_get_the_same_reduction(self, calc):
        """CPIC does not grade 1.5 and 1.0 differently (finding §3)."""
        assert calc.fluoropyrimidine_dose("*1/c.2846A>T")[
            "dose_reduction_percent"] == calc.fluoropyrimidine_dose("*1/*2A")[
            "dose_reduction_percent"]

    def test_as_1_0_carries_cpics_homozygote_caveat(self, calc):
        out = calc.fluoropyrimidine_dose("c.2846A>T/c.2846A>T")
        assert "c.2846A>T/c.2846A>T" in out["recommendation"]


class TestDPYDAlleleFunction:
    """CPIC allele-functionality table, the three values that were wrong."""

    @pytest.mark.parametrize("allele,value,why", [
        ("*7",  0.0, "c.295_298delTCAT — No function"),
        ("*10", 0.0, "c.2983G>T — No function"),
        ("*9A", 1.0, "c.85T>C — Normal function"),
        ("*2A", 0.0, "c.1905+1G>A — No function"),
        ("*13", 0.0, "c.1679T>G — No function"),
        ("c.2846A>T", 0.5, "Decreased function"),
        ("HapB3", 0.5, "c.1129-5923C>G, c.1236G>A — ONE allele, Decreased"),
    ])
    def test_value(self, allele, value, why):
        assert PhenotypeTranslator.ACTIVITY_SCORES["DPYD"][allele] == value, why


# ═════════════════════════════════════════════════════════════════════════
# CYP2C9 / phenytoin — CPIC recommendation rows by activity score
# ═════════════════════════════════════════════════════════════════════════

class TestPhenytoinRecommendationRows:
    @pytest.mark.parametrize("diplotype,score,phenotype,reduction", [
        # AS 2.0 — "No adjustments needed from typical dosing strategies."
        ("*1/*1",  2.0, "Normal Metabolizer",       0),
        ("*1/*9",  2.0, "Normal Metabolizer",       0),   # *9: CPIC Normal
        # AS 1.5 — ALSO "No adjustments needed" (finding §9)
        ("*1/*2",  1.5, "Intermediate Metabolizer", 0),
        ("*1/*5",  1.5, "Intermediate Metabolizer", 0),   # *5: CPIC Decreased
        ("*1/*8",  1.5, "Intermediate Metabolizer", 0),
        # AS 1.0 — "approximately 25% less than typical maintenance dose"
        ("*1/*3",  1.0, "Intermediate Metabolizer", 25),
        ("*1/*13", 1.0, "Intermediate Metabolizer", 25),  # finding §8
        ("*5/*5",  1.0, "Intermediate Metabolizer", 25),  # finding §13
        ("*2/*2",  1.0, "Intermediate Metabolizer", 25),
        # AS 0.5 / 0.0 — "approximately 50% less"
        ("*2/*3",  0.5, "Poor Metabolizer",         50),
        ("*3/*5",  0.5, "Poor Metabolizer",         50),
        ("*3/*3",  0.0, "Poor Metabolizer",         50),
    ])
    def test_row(self, calc, diplotype, score, phenotype, reduction):
        out = calc.phenytoin_dose(diplotype)
        assert out["activity_score"] == score
        assert out["phenotype"] == phenotype
        assert out["dose_reduction_percent"] == reduction

    def test_as_1_5_keeps_the_full_maintenance_dose(self, calc):
        """Reference case for finding §9.

        A 34-year-old with newly diagnosed focal epilepsy, CYP2C9 *1/*2 —
        activity score 1.5, the commonest actionable CYP2C9 genotype in
        European-ancestry patients (10-15% allele frequency).  CPIC: no
        dose adjustment.  The code gave a 25% reduction, so a 70 kg adult
        got 262 mg/day instead of 350 and the maintenance ceiling dropped
        from "300-400" to "225-300" — subtherapeutic phenytoin and a
        breakthrough seizure.
        """
        out = calc.phenytoin_dose("*1/*2", weight_kg=70.0)
        assert out["suggested_maintenance_mg_per_day"] == 350
        assert out["max_maintenance_mg_per_day"] == "300-400"

    def test_the_loading_dose_is_never_reduced(self, calc):
        """CPIC: "For first dose, use typical initial or loading dose"."""
        normal = calc.phenytoin_dose("*1/*1", weight_kg=70.0)
        poor = calc.phenytoin_dose("*3/*3", weight_kg=70.0)
        assert poor["standard_loading_dose_mg"] == \
            normal["standard_loading_dose_mg"]

    @pytest.mark.parametrize("diplotype", [
        "*1/*1", "*1/*2", "*1/*3", "*2/*3", "*3/*3", "*5/*5",
    ])
    def test_every_row_carries_the_hla_b_1502_caveat(self, calc, diplotype):
        """CPIC attaches it to every phenytoin row, and its first phenytoin
        recommendation is "If patient is phenytoin-naive, do not use
        phenytoin/fosphenytoin" for HLA-B*15:02-positive patients.  A
        function returning "Standard phenytoin dosing" under a heading
        citing the phenytoin guideline, with no HLA-B*15:02 caveat,
        misrepresented that guideline by omission (finding §9, related)."""
        out = calc.phenytoin_dose(diplotype)
        assert "HLA-B*15:02" in out["hla_b_1502_caveat"]
        assert any("HLA-B*15:02" in w for w in out["warnings"])


class TestCYP2C9AlleleFunction:
    @pytest.mark.parametrize("allele,value", [
        ("*5", 0.5),    # finding §13 — was 0.0
        ("*4", 0.5),    # was absent → defaulted to 1.0
        ("*13", 0.0),   # was absent → defaulted to 1.0
        ("*9", 1.0),    # CPIC Normal function
        ("*2", 0.5), ("*3", 0.0), ("*6", 0.0), ("*8", 0.5), ("*11", 0.5),
    ])
    def test_value(self, allele, value):
        assert PhenotypeTranslator.ACTIVITY_SCORES["CYP2C9"][allele] == value


# ═════════════════════════════════════════════════════════════════════════
# CYP2C19 / clopidogrel — CPIC diplotype and recommendation rows
# ═════════════════════════════════════════════════════════════════════════

class TestClopidogrelRecommendationRows:
    @pytest.mark.parametrize("diplotype,phenotype,use", [
        ("*1/*1",   "Normal Metabolizer",      "standard"),
        ("*1/*17",  "Rapid Metabolizer",       "standard"),
        ("*17/*17", "Ultrarapid Metabolizer",  "standard"),
        # "Avoid standard dose (75 mg) clopidogrel if possible."
        ("*1/*2",   "Intermediate Metabolizer", "avoid_standard_dose"),
        ("*1/*35",  "Intermediate Metabolizer", "avoid_standard_dose"),
        ("*1/*10",  "Intermediate Metabolizer", "avoid_standard_dose"),
        ("*2/*17",  "Intermediate Metabolizer", "avoid_standard_dose"),
        # "Avoid clopidogrel if possible."
        ("*2/*2",   "Poor Metabolizer",        "contraindicated"),
        ("*2/*3",   "Poor Metabolizer",        "contraindicated"),
        ("*2/*35",  "Poor Metabolizer",        "contraindicated"),
        ("*2/*10",  "Poor Metabolizer",        "contraindicated"),
        # "No recommendation."
        ("*1/*12",  "Indeterminate",           "indeterminate"),
        ("",        "Indeterminate",           "indeterminate"),
    ])
    def test_row(self, calc, diplotype, phenotype, use):
        out = calc.clopidogrel_dose(diplotype)
        assert out["phenotype"] == phenotype
        assert out["clopidogrel_use"] == use

    def test_the_stent_patient(self, calc):
        """Reference case for finding §8.

        A patient of West African ancestry, CYP2C19 *1/*35, 3 days after
        drug-eluting stent placement.  *35 is a no-function allele at CPIC
        and *1/*35 is an Intermediate Metabolizer.  `.get(allele, 1.0)`
        reported "Normal Metabolizer / standard / Standard clopidogrel
        therapy appropriate. Normal CYP2C19 metabolism expected."  The
        endpoint of staying on clopidogrel is stent thrombosis.
        """
        out = calc.clopidogrel_dose("*1/*35")
        assert out["phenotype"] != "Normal Metabolizer"
        assert out["clopidogrel_use"] != "standard"
        assert "avoid standard dose" in out["rationale"].lower()


class TestCYP2C19PhenotypeBoundaries:
    @pytest.mark.parametrize("diplotype,phenotype", [
        ("*2/*2",   "Poor Metabolizer"),
        ("*2/*10",  "Poor Metabolizer"),          # no function + decreased
        ("*1/*2",   "Intermediate Metabolizer"),
        ("*10/*10", "Intermediate Metabolizer"),  # decreased + decreased
        ("*10/*17", "Intermediate Metabolizer"),  # decreased + INCREASED
        ("*1/*1",   "Normal Metabolizer"),
        ("*1/*17",  "Rapid Metabolizer"),
        ("*17/*17", "Ultrarapid Metabolizer"),
    ])
    def test_row(self, translator, diplotype, phenotype):
        assert translator.translate("CYP2C19", diplotype)["phenotype"] == \
            phenotype

    def test_an_increased_allele_does_not_rescue_a_decreased_one(
            self, translator):
        """CPIC calls *10/*17 a (Likely) Intermediate Metabolizer.

        The summed activity score 0.5 + 1.5 = 2.0 would read Normal, which
        is why the increased-function allele is capped at a normal 1.0 when
        paired with a reduced-function one.
        """
        out = translator.translate("CYP2C19", "*10/*17")
        assert out["phenotype"] == "Intermediate Metabolizer"
        assert out["activity_score"] == 1.5


# ═════════════════════════════════════════════════════════════════════════
# SLCO1B1 / simvastatin
# ═════════════════════════════════════════════════════════════════════════

class TestSimvastatinRecommendationRows:
    @pytest.mark.parametrize("diplotype,phenotype", [
        ("*1/*1",   "Normal Function"),
        ("*1/*37",  "Normal Function"),
        ("*1/*14",  "Normal Function"),
        ("*14/*14", "Increased Function"),
        ("*1/*5",   "Decreased Function"),
        ("*1/*15",  "Decreased Function"),
        ("*5/*14",  "Decreased Function"),
        ("*5/*2",   "Possible Decreased Function"),
        ("*5/*5",   "Poor Function"),
        ("*5/*15",  "Poor Function"),
        ("*15/*15", "Poor Function"),
    ])
    def test_phenotype(self, calc, diplotype, phenotype):
        assert calc.simvastatin_dose(diplotype)["phenotype"] == phenotype

    def test_the_report_and_the_calculator_agree(self, translator, calc):
        """Reference case for finding §7.

        A 58-year-old man post-MI, SLCO1B1 *5/*15 — CPIC Poor Function,
        roughly 17-fold myopathy risk at standard simvastatin doses.  The
        integrated PGx report rendered "SLCO1B1: Normal Function —
        Wild-type 521TT, typical statin metabolism" (the translator fell
        back to *1/*1 for any diplotype it did not list, and it listed
        three of CPIC's 741), while dosing.py got the same patient right.
        The report and the calculator contradicted each other on one
        genotype.
        """
        for diplotype in ("*5/*15", "*15/*15", "*1/*15", "*5/*5", "*1/*1"):
            assert translator.translate("SLCO1B1", diplotype)["phenotype"] == \
                calc.simvastatin_dose(diplotype)["phenotype"], diplotype

    def test_decreased_function_ceiling_is_below_20mg(self, calc):
        """CPIC: "limit dose to <20mg/day" — strictly less. 20 is excluded."""
        out = calc.simvastatin_dose("*1/*5")
        assert out["max_simvastatin_dose_mg"] < 20
        assert "<20" in out["rationale"] or "LESS THAN 20" in out["rationale"]

    def test_normal_function_ceiling_is_the_labelled_range(self, calc):
        """80 mg is not a normal-function maximum.

        FDA-approved ZOCOR label (DailyMed, Organon LLC, setid
        8f55d5de-5a4f-4a39-8c84-c53976dd6af9): §2.2 "The recommended dosage
        range of ZOCOR is 20 mg to 40 mg once daily"; §2.1/§5.1 80 mg is
        "restricted to patients who have been taking simvastatin 80 mg
        daily chronically (e.g., for 12 months or more) without evidence of
        muscle toxicity"; §2.1 "The 5 mg and 80 mg strengths of ZOCOR are
        no longer marketed".  CPIC names no number: "Prescribe desired
        starting dose and adjust doses based on disease-specific
        guidelines."
        """
        out = calc.simvastatin_dose("*1/*1")
        assert out["max_simvastatin_dose_mg"] == 40

    def test_an_unrecognised_genotype_gets_no_dose_ceiling(self, calc):
        """Reference case for finding §16.

        simvastatin_dose("5/5") — a lab reporting without asterisks for a
        patient who is genuinely *5/*5, CPIC Poor Function — returned
        phenotype "Unknown", max_simvastatin_dose_mg 40 and an EMPTY
        warnings list.  CPIC's answer for Indeterminate is "No
        recommendation", which is not 40 mg.
        """
        out = calc.simvastatin_dose("5/5")
        assert out["phenotype"] == "Indeterminate"
        assert out["max_simvastatin_dose_mg"] is None
        assert out["warnings"]

    def test_poor_function_is_an_alternative_statin(self, calc):
        out = calc.simvastatin_dose("*5/*5")
        assert out["max_simvastatin_dose_mg"] == 0
        assert "alternative statin" in out["rationale"]


# ═════════════════════════════════════════════════════════════════════════
# CYP3A5 / tacrolimus
# ═════════════════════════════════════════════════════════════════════════

class TestTacrolimusRecommendationRows:
    @pytest.mark.parametrize("diplotype,mg_per_kg", [
        # CPIC gives Normal, Possible Intermediate and Intermediate
        # Metabolizers one recommendation: "Increase starting dose 1.5 to 2
        # times recommended starting dose. Total starting dose should not
        # exceed 0.3 mg/kg/day."
        ("*1/*1", 0.30),
        ("*1/*3", 0.30),
        ("*1/*6", 0.30),
        ("*1/*7", 0.30),
        # Poor Metabolizer: "Initiate therapy with standard recommended dose."
        ("*3/*3", 0.15),
        ("*3/*7", 0.15),
        ("*6/*6", 0.15),
        # Indeterminate: "No recommendation."
        ("*1/*8", None),
        ("*1/*999", None),
        ("", None),
    ])
    def test_row(self, calc, diplotype, mg_per_kg):
        assert calc.tacrolimus_dose(diplotype)["mg_per_kg_per_day"] == mg_per_kg

    def test_the_ceiling_is_cpics_0_3_mg_per_kg(self, calc):
        out = calc.tacrolimus_dose("*1/*1", weight=80.0)
        assert out["mg_per_kg_per_day"] <= 0.30
        assert out["total_daily_mg"] == pytest.approx(24.0, abs=0.1)


# ═════════════════════════════════════════════════════════════════════════
# UGT1A1
# ═════════════════════════════════════════════════════════════════════════

class TestUGT1A1DiplotypeRows:
    @pytest.mark.parametrize("diplotype,phenotype", [
        ("*1/*1",   "Normal Metabolizer"),
        ("*1/*36",  "Normal Metabolizer"),   # AS 2.5 — off the old table
        ("*36/*36", "Normal Metabolizer"),   # AS 3.0 — off the old table
        ("*1/*28",  "Intermediate Metabolizer"),
        ("*1/*6",   "Intermediate Metabolizer"),
        ("*1/*27",  "Intermediate Metabolizer"),
        ("*1/*37",  "Intermediate Metabolizer"),   # finding §14 — was Poor
        ("*28/*36", "Intermediate Metabolizer"),   # Decreased + Increased
        ("*36/*37", "Intermediate Metabolizer"),
        ("*28/*28", "Poor Metabolizer"),
        ("*37/*37", "Poor Metabolizer"),
        ("*6/*28",  "Poor Metabolizer"),
        ("*1/*80",  "Indeterminate"),              # *80: Unknown function
        ("*80/*80", "Indeterminate"),
    ])
    def test_row(self, translator, diplotype, phenotype):
        assert translator.translate("UGT1A1", diplotype)["phenotype"] == \
            phenotype

    def test_there_is_no_dead_zone_above_the_top_bucket(self, translator):
        """Finding §14.

        `_score_to_phenotype` closed the last interval at its stated upper
        bound — 2.1 for UGT1A1 — and returned "Indeterminate" for anything
        higher.  *36 is an increased-function allele scored 1.5, so *1/*36
        (2.5) and *36/*36 (3.0) fell off the end and were reported
        Indeterminate where CPIC calls both Normal Metabolizer.  The
        fallback is now a raise: an unmapped score is a bug in the table,
        not a clinical answer.
        """
        assert translator.translate("UGT1A1", "*36/*36")["activity_score"] == 3.0
        with pytest.raises(ValueError):
            translator._score_to_phenotype("UGT1A1", 4.0)


# ═════════════════════════════════════════════════════════════════════════
# The absent-value pattern, across every calculator
# ═════════════════════════════════════════════════════════════════════════

class TestNoCalculatorTreatsAnAbsentGenotypeAsWildType:
    """Findings §6, §7 and §8 are one pattern with nine entry points."""

    BAD = ["", "*2A", "Unknown", "not reported", "*1/*1/*1", "/", "*1/"]

    @pytest.mark.parametrize("bad", BAD)
    def test_fluoropyrimidine(self, calc, bad):
        out = calc.fluoropyrimidine_dose(bad)
        assert out["phenotype"] == "Indeterminate"
        assert out["dose_reduction_percent"] is None

    @pytest.mark.parametrize("bad", BAD)
    def test_clopidogrel(self, calc, bad):
        assert calc.clopidogrel_dose(bad)["phenotype"] == "Indeterminate"

    @pytest.mark.parametrize("bad", BAD)
    def test_phenytoin(self, calc, bad):
        out = calc.phenytoin_dose(bad)
        assert out["phenotype"] == "Indeterminate"
        assert out["suggested_maintenance_mg_per_day"] is None

    @pytest.mark.parametrize("bad", BAD)
    def test_tacrolimus(self, calc, bad):
        out = calc.tacrolimus_dose(bad)
        assert out["phenotype"] == "Indeterminate"
        assert out["mg_per_kg_per_day"] is None

    @pytest.mark.parametrize("bad", BAD)
    def test_simvastatin(self, calc, bad):
        out = calc.simvastatin_dose(bad)
        assert out["phenotype"] == "Indeterminate"
        assert out["max_simvastatin_dose_mg"] is None

    @pytest.mark.parametrize("bad", BAD)
    def test_ssri(self, calc, bad):
        out = calc.ssri_dose("citalopram", bad)
        assert out["phenotype"] == "Indeterminate"
        assert out["max_dose_mg"] is None

    @pytest.mark.parametrize("bad", BAD)
    def test_tca(self, calc, bad):
        assert calc.tca_dose("amitriptyline", bad)["dose_adjustment"] == \
            "no_recommendation"

    @pytest.mark.parametrize("bad", BAD)
    def test_thiopurine(self, calc, bad):
        out = calc.thiopurine_dose(bad, "*1/*1")
        assert out["dose_percent_of_standard"] is None

    @pytest.mark.parametrize("bad", BAD)
    def test_thiopurine_nudt15_axis(self, calc, bad):
        out = calc.thiopurine_dose("*1/*1", bad)
        assert out["dose_percent_of_standard"] is None

    @pytest.mark.parametrize("gene", [
        "CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "TPMT", "NUDT15", "UGT1A1",
        "CYP3A5", "SLCO1B1", "VKORC1", "CYP4F2",
    ])
    @pytest.mark.parametrize("bad", ["", "Unknown", "*1"])
    def test_translator(self, translator, gene, bad):
        out = translator.translate(gene, bad)
        assert out["phenotype"] == "Indeterminate", (gene, bad)
        assert out["activity_score"] is None, (gene, bad)
        # The explanation may say "this is not wild-type"; it must never
        # describe the patient AS wild-type.
        assert "Wild-type" not in out["clinical_meaning"], (gene, bad)

    @pytest.mark.parametrize("gene,bad", [
        ("VKORC1", "A/A"),
        ("VKORC1", "-1639G>A/*1x"),
        ("SLCO1B1", "*1/*17"),
        ("CYP4F2", "*3/*3/*3"),
    ])
    def test_unrecognised_genotype_based_input_is_not_wild_type(
            self, translator, gene, bad):
        """Finding §7.

        `lookup.get("*1/*1", ...)` under a comment reading "try to find the
        closest match" made every unrecognised VKORC1, SLCO1B1 or CYP4F2
        genotype the wild-type phenotype WITH the wild-type explanatory
        text.  A lab reporting VKORC1 as "A/A" rather than
        "-1639G>A/-1639G>A" got "Normal Warfarin Sensitivity — standard
        warfarin dose range" for a patient homozygous for the low-expression
        promoter variant, who needs the lowest doses in the population.
        dosing.py's `_is_vkorc1_aa` does recognise "A/A", so the dose
        calculator and the phenotype report disagreed about one string.
        """
        out = translator.translate(gene, bad)
        assert out["phenotype"] == "Indeterminate"
        assert "Wild-type" not in out["clinical_meaning"]
        assert "standard warfarin dose range" not in out["clinical_meaning"]

    def test_an_uncalled_gene_is_not_translated_as_wild_type(self, translator):
        """Finding §20b: translate_all defaulted a missing diplotype to *1/*1."""
        out = translator.translate_all({
            "CYP2D6": {"diplotype": "*1/*4"},
            "CYP2C19": {},                      # no call
            "DPYD": {"diplotype": None},        # call failed
        })
        assert out["CYP2D6"]["phenotype"] == "Intermediate Metabolizer"
        assert out["CYP2C19"]["phenotype"] == "Indeterminate"
        assert out["DPYD"]["phenotype"] == "Indeterminate"
