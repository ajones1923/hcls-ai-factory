"""CYP2D6 activity score to phenotype, pinned to the CPIC/DPWG consensus table.

CLINICAL_REVIEW §3 recorded a contradiction inside this one agent: `pgx_pipeline.py` put
activity score 1.0 in **Normal Metabolizer** while `knowledge.py`, two files away, defined NM as
1.25-2.25. The consensus paper settles it in its own words -- "reclassify AS of 1 from NM to IM"
-- and this file holds both halves of the repo to it.

Two boundaries were wrong, not one:

  AS 1.00  (*1/*4, *1/*5, *41/*41)   NM -> IM    standard dosing where CPIC asks for a reduction
  AS 2.25  (*1x2/*41)                UM -> NM    2.25 is the inclusive TOP of Normal

Source: Caudle KE, Sangkuhl K, Whirl-Carrillo M, Swen JJ, Haidar CE, Klein TE, Gammal RS,
Relling MV, Scott SA, Hertz DL, Guchelaar HJ, Gaedigk A. "Standardizing CYP2D6 Genotype to
Phenotype Translation: Consensus Recommendations from the Clinical Pharmacogenetics
Implementation Consortium and Dutch Pharmacogenetics Working Group."
Clin Transl Sci 2020;13(1):116-124, Table 3. doi:10.1111/cts.12692, PMID 31647186.
https://pmc.ncbi.nlm.nih.gov/articles/PMC6951851/
"""
import pytest

from src.dosing import DosingCalculator
from src.knowledge import ACTIVITY_SCORE_TABLES
from src.pgx_pipeline import PhenotypeTranslator

# Table 3, consensus column: PM = 0, IM = 0 < AS < 1.25, NM = 1.25 <= AS <= 2.25, UM > 2.25.
CONSENSUS = [
    (0.00, "Poor Metabolizer"),
    (0.25, "Intermediate Metabolizer"),
    (0.50, "Intermediate Metabolizer"),
    (0.75, "Intermediate Metabolizer"),
    (1.00, "Intermediate Metabolizer"),   # the reclassification
    (1.25, "Normal Metabolizer"),         # the floor of Normal
    (1.50, "Normal Metabolizer"),
    (2.00, "Normal Metabolizer"),
    (2.25, "Normal Metabolizer"),         # inclusive top of Normal
    (2.50, "Ultrarapid Metabolizer"),
    (3.00, "Ultrarapid Metabolizer"),
]


@pytest.fixture
def translator():
    return PhenotypeTranslator()


class TestTheTranslatorFollowsTheConsensusTable:
    @pytest.mark.parametrize("score,phenotype", CONSENSUS)
    def test_activity_score_maps_to_the_consensus_phenotype(self, translator, score,
                                                            phenotype):
        assert translator._score_to_phenotype("CYP2D6", score) == phenotype

    def test_as_one_is_intermediate_not_normal(self, translator):
        """The single assertion CLINICAL_REVIEW §3 was about."""
        assert translator._score_to_phenotype("CYP2D6", 1.0) == "Intermediate Metabolizer"
        assert translator._score_to_phenotype("CYP2D6", 1.0) != "Normal Metabolizer"

    def test_two_twenty_five_is_normal_not_ultrarapid(self, translator):
        assert translator._score_to_phenotype("CYP2D6", 2.25) == "Normal Metabolizer"


class TestTheDiplotypesTheConsensusNames:
    """Table 3 lists example diplotypes per phenotype; these are the ones this repo can score."""

    @pytest.mark.parametrize("diplotype,score,phenotype", [
        ("*4/*4",   0.0, "Poor Metabolizer"),
        ("*4/*41",  0.5, "Intermediate Metabolizer"),
        # This repo scores *41 at 0.5, where the CPIC allele-functionality table assigns 0.25.
        # The phenotype is Intermediate either way (AS 1.0 here, AS 0.5 there), so the
        # difference does not change this result -- it is recorded in CLINICAL_REVIEW §3a
        # rather than changed, because the allele table itself was not retrieved.
        ("*41/*41", 1.0, "Intermediate Metabolizer"),
        ("*1/*4",   1.0, "Intermediate Metabolizer"),
        ("*1/*1",   2.0, "Normal Metabolizer"),
    ])
    def test_diplotype(self, translator, diplotype, score, phenotype):
        result = translator.translate("CYP2D6", diplotype)
        assert result["activity_score"] == score, diplotype
        assert result["phenotype"] == phenotype, diplotype


class TestTheDosingPathAgrees:
    """The TCA calculator carried its own copy of these thresholds, with the same two errors."""

    @staticmethod
    def _tca_phenotype(diplotype: str) -> str:
        out = DosingCalculator().tca_dose(
            drug="amitriptyline", cyp2d6_diplotype=diplotype, cyp2c19_diplotype="*1/*1")
        return out["cyp2d6_phenotype"]

    @pytest.mark.parametrize("diplotype,phenotype", [
        ("*4/*4",  "Poor Metabolizer"),
        ("*1/*4",  "Intermediate Metabolizer"),   # AS 1.0 — was Normal
        ("*1/*1",  "Normal Metabolizer"),
    ])
    def test_tca_dosing_uses_the_consensus_phenotype(self, diplotype, phenotype):
        assert self._tca_phenotype(diplotype) == phenotype

    def test_an_as_one_patient_no_longer_gets_standard_dosing(self):
        """CPIC asks for a 25% reduction in CYP2D6 IM; the old boundary gave standard dosing."""
        assert self._tca_phenotype("*1/*4") == "Intermediate Metabolizer"


class TestTheKnowledgeBaseAndTheCodeNoLongerDisagree:
    """§3 was a contradiction between two files, so the test covers both sides."""

    def test_the_knowledge_base_states_the_consensus_ranges(self):
        thresholds = ACTIVITY_SCORE_TABLES["CYP2D6"]["phenotype_thresholds"]
        assert thresholds["Normal Metabolizer"].startswith("1.25")
        assert thresholds["Poor Metabolizer"] == "0"

    def test_the_documented_range_and_the_computed_phenotype_agree_at_the_boundary(
            self, translator):
        """knowledge.py says NM starts at 1.25; the translator must not start it earlier."""
        assert translator._score_to_phenotype("CYP2D6", 1.0) != "Normal Metabolizer"
        assert translator._score_to_phenotype("CYP2D6", 1.25) == "Normal Metabolizer"
