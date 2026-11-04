"""NUDT15 and TPMT thiopurine phenotypes, pinned to CPIC Table 1.

The defect this file exists for: `PhenotypeTranslator.ACTIVITY_SCORES` scored NUDT15 `*2` at
**0.5**, a decreased-function value. CPIC Table 1 lists `*2` as a **no function** allele. So
`*2/*2` -- "two no function alleles", a CPIC **poor metabolizer** -- summed to 1.0 and reported
**"NUDT15 Intermediate"**.

A NUDT15 poor metabolizer reported as intermediate is told 30-70% of a standard thiopurine dose
where CPIC says avoid the drug, or 10% with intensive monitoring. What follows a full-ish dose
in a NUDT15 PM is life-threatening myelosuppression. `dosing.py` had `*2` right at 0.0 the whole
time -- two copies of one table, disagreeing, and the one a clinician reads was the wrong one.

Three further errors came out of reading the table rather than summarising it:

  *4, *5     scored 0.0 (no function) in dosing.py -- CPIC says UNCLEAR function
  *6         scored 1.0 (normal) in both -- CPIC says UNCLEAR function
  unknown    any unrecognised allele scored 1.0, i.e. reported as fully functional

CPIC assigns **Indeterminate** to a diplotype carrying an unclear-function allele. An activity
score cannot express that, so the translator and the dose calculator now both refuse to produce
one rather than inventing a number in the middle.

Source: Relling MV, Schwab M, Whirl-Carrillo M, Suarez-Kurtz G, Pui CH, Stein CM, Moyer AM,
Evans WE, Klein TE, Antillon-Klussmann FG, Caudle KE, Kato M, Yeoh AEJ, Schmiegelow K, Yang JJ.
"Clinical Pharmacogenetics Implementation Consortium Guideline for Thiopurine Dosing Based on
TPMT and NUDT15 Genotypes: 2018 Update." Clin Pharmacol Ther 2019;105(5):1095-1105, Table 1.
doi:10.1002/cpt.1304, PMID 30447069. https://pmc.ncbi.nlm.nih.gov/articles/PMC6576267/
"""
import pytest

from src.dosing import _NUDT15_ACTIVITY, _NUDT15_UNCERTAIN, DosingCalculator
from src.pgx_pipeline import PhenotypeTranslator


@pytest.fixture
def translator():
    return PhenotypeTranslator()


@pytest.fixture
def calc():
    return DosingCalculator()


class TestTheAlleleFunctionAssignments:
    """Table 1 assigns a function to three NUDT15 alleles. It does not assign one to the rest."""

    @pytest.mark.parametrize("allele", ["*2", "*3"])
    def test_the_no_function_alleles_score_zero(self, allele):
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"][allele] == 0.0

    def test_star1_is_the_only_functional_allele(self):
        scores = PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]
        assert scores["*1"] == 1.0
        assert [a for a, v in scores.items() if v > 0] == ["*1"]

    def test_star2_is_not_a_decreased_function_allele(self):
        """The 0.5 that made *2/*2 come out Intermediate."""
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]["*2"] != 0.5

    @pytest.mark.parametrize("allele", ["*4", "*5", "*6", "*7", "*8", "*9"])
    def test_the_unclear_function_alleles_carry_no_score_at_all(self, allele):
        """"these variants (*4 to *9) are designated as unclear function" -- Relling 2019."""
        assert allele not in PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]
        assert allele in PhenotypeTranslator.UNCERTAIN_FUNCTION_ALLELES["NUDT15"]
        assert allele not in _NUDT15_ACTIVITY
        assert allele in _NUDT15_UNCERTAIN


class TestThePhenotypesTable1Names:
    @pytest.mark.parametrize("diplotype,phenotype", [
        ("*1/*1", "NUDT15 Normal"),            # two normal function alleles
        ("*1/*2", "NUDT15 Intermediate"),      # one normal + one no function
        ("*1/*3", "NUDT15 Intermediate"),
        ("*2/*2", "NUDT15 Deficient"),         # two no function -- the fixed case
        ("*3/*3", "NUDT15 Deficient"),
        ("*2/*3", "NUDT15 Deficient"),
    ])
    def test_translator_phenotype(self, translator, diplotype, phenotype):
        assert translator.translate("NUDT15", diplotype)["phenotype"] == phenotype

    def test_two_no_function_alleles_are_not_intermediate(self, translator):
        """The regression this file is named for."""
        result = translator.translate("NUDT15", "*2/*2")
        assert result["activity_score"] == 0.0
        assert result["phenotype"] != "NUDT15 Intermediate"


class TestTheIndeterminateCases:
    @pytest.mark.parametrize("diplotype", ["*1/*4", "*1/*5", "*4/*5", "*5/*6", "*2/*5", "*1/*6"])
    def test_an_unclear_allele_yields_indeterminate_not_a_guess(self, translator, diplotype):
        result = translator.translate("NUDT15", diplotype)
        assert result["phenotype"] == "Indeterminate"
        assert result["activity_score"] is None

    def test_the_offending_allele_is_named(self, translator):
        result = translator.translate("NUDT15", "*1/*6")
        assert result["unscoreable_alleles"] == ["*6"]
        assert "uncertain function" in result["clinical_meaning"]

    def test_an_unrecognised_allele_is_not_assumed_normal(self, translator):
        """`.get(allele, 1.0)` used to report a never-seen allele as fully functional."""
        result = translator.translate("NUDT15", "*1/*999")
        assert result["phenotype"] == "Indeterminate"
        assert result["activity_score"] is None
        assert "Do not infer normal function" in result["clinical_meaning"]

    def test_star1_slash_star6_no_longer_reads_as_normal(self, translator):
        assert translator.translate("NUDT15", "*1/*6")["phenotype"] != "NUDT15 Normal"


class TestTheDoseRecommendationFollows:
    @pytest.mark.parametrize("diplotype,percent", [
        ("*1/*1", 100),
        ("*1/*2", 50),
        ("*2/*2", 10),
        ("*2/*3", 10),
    ])
    def test_dose_percent_of_standard(self, calc, diplotype, percent):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype=diplotype)
        assert out["dose_percent_of_standard"] == percent

    def test_a_poor_metabolizer_is_told_to_avoid_the_drug(self, calc):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype="*2/*2")
        assert out["dose_percent_of_standard"] == 10
        assert "AVOID" in out["dose_recommendation"]

    @pytest.mark.parametrize("diplotype", ["*1/*4", "*1/*6", "*2/*5"])
    def test_an_unclear_allele_gets_no_number(self, calc, diplotype):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype=diplotype)
        assert out["dose_percent_of_standard"] is None
        assert out["nudt15_phenotype"] == "NUDT15 Indeterminate"
        assert "uncertain function" in out["dose_recommendation"]

    def test_an_unrecognised_tpmt_allele_gets_no_number_either(self, calc):
        """TPMT names five alleles here and dozens in CPIC's table."""
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*31", nudt15_diplotype="*1/*1")
        assert out["dose_percent_of_standard"] is None
        assert out["tpmt_phenotype"] == "TPMT Indeterminate"
        assert "*31" in out["dose_recommendation"]

    def test_the_guideline_is_cited_in_the_indeterminate_result(self, calc):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype="*1/*4")
        assert "Relling" in out["guideline"]


class TestTheTwoCopiesOfTheTableAgree:
    """The root cause: one clinical table, two implementations, silently divergent."""

    def test_nudt15_scores_match_between_the_translator_and_the_dose_calculator(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"] == _NUDT15_ACTIVITY

    def test_the_uncertain_sets_match(self):
        assert set(PhenotypeTranslator.UNCERTAIN_FUNCTION_ALLELES["NUDT15"]) == set(
            _NUDT15_UNCERTAIN)
