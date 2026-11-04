"""NUDT15 and TPMT thiopurine phenotypes, pinned to CPIC's MAINTAINED allele table.

CORRECTED 2026-11-03 (pgx audit §5). This file was written for CLINICAL_REVIEW §3c and
sourced to Relling 2019 Table 1, which in 2018 designated NUDT15 `*4` to `*9` "unclear
function". It therefore asserted, as the specification, that those six alleles must carry
**no activity value at all** and that `*1/*4`, `*1/*6`, `*2/*5` and `*3/*6` must come back
`Indeterminate` with no dose.

**CPIC has since assigned them.** A guideline paper and the table the issuing body maintains
and republishes are different artefacts, and for allele function the maintained one governs.
CPIC's current NUDT15 allele-functionality table says:

    *4, *6, *7, *9, *10, *13, *14, *17, *18, *19   No function
    *5, *8                                          Decreased function
    *11, *12, *15, *16, *20, *21, *22, *23          Uncertain function

and its diplotype table follows from that without ambiguity:

    *1/*4  -> Intermediate Metabolizer      *2/*4  -> Poor Metabolizer
    *1/*6  -> Intermediate Metabolizer      *3/*6  -> Poor Metabolizer
    *1/*5  -> Normal Metabolizer            *6/*6  -> Poor Metabolizer
    *1/*8  -> Normal Metabolizer            *1/*11 -> Indeterminate

So the §3c fix re-created, through a different door, the exact harm it was written to
prevent: a 9-year-old of Han Chinese ancestry starting 6-mercaptopurine maintenance for ALL
with NUDT15 `*3/*6` -- a CPIC **Poor Metabolizer** -- was told "no genotype-guided dose can
be given ... dose on clinical grounds with the standard CBC monitoring", which for a
thiopurine in a child is a full weight-based dose and profound, prolonged myelosuppression.
`*1/*4` and `*1/*6`, CPIC Intermediate Metabolizers who should start at a reduced dose, got
the same answer. The caveat text was worse than silence: it was factually wrong as of
CPIC's current table ("has uncertain function in CPIC's guideline") and it read as a
confident justification for doing nothing.

Indeterminate is still the right answer -- for the eight alleles CPIC still calls uncertain,
which the old set did not contain and which the "not in the allele table" branch had been
catching correctly all along.

The cases this file still pins from the original §3c defect are unchanged and still pass:
`*2` is a no-function allele, `*2/*2` is a Poor Metabolizer and not an Intermediate one, and
an unrecognised allele is not assumed normal.

Primary sources
---------------
CPIC allele-functionality table:
    https://api.cpicpgx.org/v1/allele?genesymbol=eq.NUDT15  (retrieved 2026-11-03)
CPIC diplotype-to-phenotype table (276 rows, every one reproduced by the values and
thresholds under test):
    https://api.cpicpgx.org/v1/diplotype?genesymbol=eq.NUDT15  (retrieved 2026-11-03)
CPIC azathioprine / mercaptopurine recommendation tables (the two-gene dose rows):
    https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:1256   (azathioprine)
    https://api.cpicpgx.org/v1/recommendation?drugid=eq.RxNorm:103    (mercaptopurine)
Guideline paper, for the framework rather than the allele values: Relling MV, Schwab M,
Whirl-Carrillo M, et al. "Clinical Pharmacogenetics Implementation Consortium Guideline for
Thiopurine Dosing Based on TPMT and NUDT15 Genotypes: 2018 Update." Clin Pharmacol Ther
2019;105(5):1095-1105. doi:10.1002/cpt.1304, PMID 30447069.
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
    """CPIC's current table, allele by allele."""

    @pytest.mark.parametrize("allele", [
        "*2", "*3", "*4", "*6", "*7", "*9", "*10", "*13", "*14", "*17", "*18", "*19",
    ])
    def test_the_no_function_alleles_score_zero(self, allele):
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"][allele] == 0.0

    @pytest.mark.parametrize("allele", ["*5", "*8"])
    def test_the_decreased_function_alleles_score_a_half(self, allele):
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"][allele] == 0.5

    def test_star1_is_the_only_normal_function_allele(self):
        scores = PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]
        assert scores["*1"] == 1.0
        assert [a for a, v in scores.items() if v >= 1.0] == ["*1"]

    def test_star2_is_not_a_decreased_function_allele(self):
        """The 0.5 that made *2/*2 come out Intermediate (CLINICAL_REVIEW §3c)."""
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]["*2"] != 0.5

    @pytest.mark.parametrize("allele", ["*4", "*6", "*7", "*9"])
    def test_the_superseded_alleles_are_no_longer_treated_as_unscoreable(self, allele):
        """These four are the §5 regression: CPIC assigns them NO function."""
        assert allele not in PhenotypeTranslator.UNCERTAIN_FUNCTION_ALLELES["NUDT15"]
        assert allele not in _NUDT15_UNCERTAIN
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"][allele] == 0.0
        assert _NUDT15_ACTIVITY[allele] == 0.0

    @pytest.mark.parametrize("allele", [
        "*11", "*12", "*15", "*16", "*20", "*21", "*22", "*23",
    ])
    def test_the_uncertain_function_alleles_carry_no_score_at_all(self, allele):
        """CPIC's Indeterminate set -- and only CPIC's Indeterminate set."""
        assert allele not in PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"]
        assert allele in PhenotypeTranslator.UNCERTAIN_FUNCTION_ALLELES["NUDT15"]
        assert allele not in _NUDT15_ACTIVITY
        assert allele in _NUDT15_UNCERTAIN


class TestThePhenotypesCPICsDiplotypeTableNames:
    @pytest.mark.parametrize("diplotype,phenotype", [
        # unchanged by this correction
        ("*1/*1", "NUDT15 Normal"),            # normal  + normal
        ("*1/*2", "NUDT15 Intermediate"),      # normal  + no function
        ("*1/*3", "NUDT15 Intermediate"),
        ("*2/*2", "NUDT15 Deficient"),         # the §3c case
        ("*3/*3", "NUDT15 Deficient"),
        ("*2/*3", "NUDT15 Deficient"),
        # the §5 regression: previously all six were "Indeterminate"
        ("*1/*4", "NUDT15 Intermediate"),      # normal + no function
        ("*1/*6", "NUDT15 Intermediate"),
        ("*1/*5", "NUDT15 Normal"),            # normal + decreased
        ("*1/*8", "NUDT15 Normal"),
        ("*3/*6", "NUDT15 Deficient"),         # no function + no function
        ("*2/*4", "NUDT15 Deficient"),
        ("*6/*6", "NUDT15 Deficient"),
        ("*2/*5", "NUDT15 Deficient"),         # no function + decreased
        ("*5/*5", "NUDT15 Intermediate"),      # decreased + decreased
        ("*5/*8", "NUDT15 Intermediate"),
    ])
    def test_translator_phenotype(self, translator, diplotype, phenotype):
        assert translator.translate("NUDT15", diplotype)["phenotype"] == phenotype

    def test_two_no_function_alleles_are_not_intermediate(self, translator):
        """The §3c regression this file was originally named for."""
        result = translator.translate("NUDT15", "*2/*2")
        assert result["activity_score"] == 0.0
        assert result["phenotype"] != "NUDT15 Intermediate"

    def test_a_no_function_plus_decreased_pair_is_poor_not_intermediate(self, translator):
        """NUDT15 differs from TPMT here: CPIC calls *2/*5 a Poor Metabolizer."""
        assert translator.translate("NUDT15", "*2/*5")["phenotype"] == "NUDT15 Deficient"

    def test_a_decreased_plus_normal_pair_is_normal_not_intermediate(self, translator):
        """CPIC calls *1/*5 a Normal Metabolizer; AS 1.5 must not read Intermediate."""
        out = translator.translate("NUDT15", "*1/*5")
        assert out["activity_score"] == 1.5
        assert out["phenotype"] == "NUDT15 Normal"


class TestTheIndeterminateCases:
    @pytest.mark.parametrize("diplotype", ["*1/*11", "*1/*15", "*11/*12", "*2/*20"])
    def test_an_uncertain_allele_yields_indeterminate_not_a_guess(
            self, translator, diplotype):
        result = translator.translate("NUDT15", diplotype)
        assert result["phenotype"] == "Indeterminate"
        assert result["activity_score"] is None

    def test_the_offending_allele_is_named(self, translator):
        result = translator.translate("NUDT15", "*1/*11")
        assert result["unscoreable_alleles"] == ["*11"]
        assert "uncertain function" in result["clinical_meaning"]

    def test_an_unrecognised_allele_is_not_assumed_normal(self, translator):
        """`.get(allele, 1.0)` used to report a never-seen allele as fully functional."""
        result = translator.translate("NUDT15", "*1/*999")
        assert result["phenotype"] == "Indeterminate"
        assert result["activity_score"] is None
        assert "Do not infer normal function" in result["clinical_meaning"]

    def test_star1_slash_star6_no_longer_reads_as_indeterminate(self, translator):
        """*6 is no function at CPIC; *1/*6 is an Intermediate Metabolizer."""
        assert translator.translate("NUDT15", "*1/*6")["phenotype"] == "NUDT15 Intermediate"


class TestTheDoseRecommendationFollows:
    """CPIC's azathioprine and mercaptopurine recommendation rows, by phenotype pair."""

    @pytest.mark.parametrize("diplotype,percent", [
        ("*1/*1", 100),   # NM + NM: standard starting dose
        ("*1/*2", 30),    # NM + IM: 30-80% of standard, reported at the low end
        ("*1/*4", 30),    # the §5 regression -- previously no dose at all
        ("*1/*6", 30),
        ("*1/*5", 100),   # NM + NM (decreased + normal): standard
        ("*2/*2", None),  # NM + PM: alternative therapy, no azathioprine dose
        ("*2/*3", None),
        ("*3/*6", None),  # the §5 clinical scenario
    ])
    def test_dose_percent_of_standard(self, calc, diplotype, percent):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype=diplotype)
        assert out["dose_percent_of_standard"] == percent

    def test_a_poor_metabolizer_is_told_to_consider_an_alternative(self, calc):
        """CPIC azathioprine, NUDT15 PM: "Consider alternative nonthiopurine
        immunosuppressant therapy." -- not a reduced azathioprine dose."""
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype="*2/*2")
        assert out["dose_percent_of_standard"] is None
        assert "alternative nonthiopurine" in out["dose_recommendation"]

    def test_the_poor_metabolizer_who_used_to_be_sent_away_now_gets_an_answer(self, calc):
        """NUDT15 *3/*6: the 9-year-old on 6-MP in the §5 failure scenario."""
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype="*3/*6")
        assert out["nudt15_phenotype"] == "NUDT15 Deficient"
        assert "clinical grounds" not in out["dose_recommendation"]
        assert "alternative nonthiopurine" in out["dose_recommendation"]

    def test_mercaptopurine_poor_metabolizers_get_the_10_fold_reduction(self, calc):
        """CPIC puts the 10-fold/thrice-weekly figure on mercaptopurine and
        thioguanine, NOT on azathioprine."""
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype="*2/*2")
        by_drug = out["recommendations_by_drug"]
        assert "10-fold" in by_drug["6-mercaptopurine (6-MP)"]
        assert "10-fold" not in by_drug["azathioprine"]

    @pytest.mark.parametrize("diplotype", ["*1/*11", "*1/*15", "*2/*22"])
    def test_an_uncertain_allele_gets_no_number(self, calc, diplotype):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype=diplotype)
        assert out["dose_percent_of_standard"] is None
        assert out["nudt15_phenotype"] == "NUDT15 Indeterminate"
        assert "uncertain function" in out["dose_recommendation"]

    def test_an_unrecognised_tpmt_allele_gets_no_number_either(self, calc):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*77", nudt15_diplotype="*1/*1")
        assert out["dose_percent_of_standard"] is None
        assert out["tpmt_phenotype"] == "TPMT Indeterminate"
        assert "*77" in out["dose_recommendation"]

    def test_the_guideline_is_cited_in_the_indeterminate_result(self, calc):
        out = calc.thiopurine_dose(tpmt_diplotype="*1/*1", nudt15_diplotype="*1/*11")
        assert "Relling" in out["guideline"]


class TestTPMTIsNotNUDT15:
    """The two genes' function-pair rules differ, and CPIC says so."""

    def test_tpmt_no_function_plus_decreased_is_intermediate(self, translator):
        """CPIC TPMT *2/*8 -> Intermediate Metabolizer (NUDT15's *2/*5 -> Poor)."""
        assert translator.translate("TPMT", "*2/*8")["phenotype"] == "TPMT Intermediate"

    def test_tpmt_decreased_plus_normal_is_normal(self, translator):
        out = translator.translate("TPMT", "*1/*8")
        assert out["activity_score"] == 1.5
        assert out["phenotype"] == "TPMT Normal"


class TestTheTwoCopiesOfTheTableAgree:
    """The root cause: one clinical table, two implementations, silently divergent."""

    def test_nudt15_scores_match_between_the_translator_and_the_dose_calculator(self):
        assert PhenotypeTranslator.ACTIVITY_SCORES["NUDT15"] == _NUDT15_ACTIVITY

    def test_the_uncertain_sets_match(self):
        assert set(PhenotypeTranslator.UNCERTAIN_FUNCTION_ALLELES["NUDT15"]) == set(
            _NUDT15_UNCERTAIN)
