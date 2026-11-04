"""One clinical table, three copies: they must not disagree.

Allele activity values live in three places in this agent -- `PhenotypeTranslator.ACTIVITY_SCORES`
(pgx_pipeline.py), the `_<GENE>_ACTIVITY` dicts (dosing.py), and
`ACTIVITY_SCORE_TABLES[gene]["allele_scores"]` (knowledge.py). Nothing reconciled them, and twice
now they have drifted into a clinical defect:

  CYP2D6   both copies of PHENOTYPE_THRESHOLDS carried the same two wrong boundaries
  NUDT15   the translator scored *2 at 0.5 and the dose calculator at 0.0; the translator's
           value turned a CPIC poor metabolizer into "Intermediate" (CLINICAL_REVIEW §3c)

The fix for the duplication itself is a refactor -- one table, imported -- which is not something
to do in the same change as a clinical correction. This test is the guard in the meantime: it
fails the moment two copies disagree about an allele they both define, or one grows an allele the
other does not know.

Divergences that exist and are deliberate (or recorded, not yet settled) are listed in
KNOWN_ASYMMETRIES with the reason. Adding to that list is a decision, which is the point.
"""
import pytest

from src import dosing
from src.knowledge import ACTIVITY_SCORE_TABLES
from src.pgx_pipeline import PhenotypeTranslator

TRANSLATOR = PhenotypeTranslator.ACTIVITY_SCORES

# gene -> alleles one copy defines and another does not, with why that is acceptable today.
KNOWN_ASYMMETRIES = {
    # CLOSED 2026-11-03 (pgx audit §4, §17): the CYP3A5 and DPYD entries that stood here are
    # gone. Both said, in effect, "the CPIC table has not been retrieved". It has now --
    # api.cpicpgx.org/v1/allele?genesymbol=eq.CYP3A5 gives *7 | No function, and
    # ?genesymbol=eq.DPYD gives *7 and *10 No function and *9A Normal function -- so the
    # translator's copy was filled in from the source rather than left thinner, and the
    # allowlist entries describe nothing any more.
    #
    # knowledge.py's CYP2D6 list predates two alleles the translator scores.
    ("CYP2D6", "knowledge"): {"*40", "*45"},
}


def _dosing_table(gene: str):
    return getattr(dosing, f"_{gene}_ACTIVITY", None)


def _knowledge_table(gene: str):
    entry = ACTIVITY_SCORE_TABLES.get(gene)
    if not entry:
        return None
    scores = entry.get("allele_scores")
    # knowledge.py keys DPYD by c. nomenclature; comparing those to star alleles is noise.
    return scores if gene != "DPYD" else None


GENES_WITH_A_SECOND_COPY = sorted(
    g for g in TRANSLATOR if _dosing_table(g) is not None or _knowledge_table(g) is not None
)


class TestTheCopiesAgreeOnEveryAlleleTheyShare:
    """A value disagreement is the dangerous case: both tables answer, differently."""

    @pytest.mark.parametrize("gene", GENES_WITH_A_SECOND_COPY)
    def test_dosing_agrees_with_the_translator(self, gene):
        other = _dosing_table(gene)
        if other is None:
            pytest.skip(f"dosing.py has no table for {gene}")
        disagreements = {
            a: (TRANSLATOR[gene][a], other[a])
            for a in sorted(set(TRANSLATOR[gene]) & set(other))
            if TRANSLATOR[gene][a] != other[a]
        }
        assert not disagreements, (
            f"{gene}: pgx_pipeline.py and dosing.py disagree about "
            f"{ {a: f'{t} vs {d}' for a, (t, d) in disagreements.items()} }. One of them is "
            f"wrong, and whichever a clinician reads is the one that matters.")

    @pytest.mark.parametrize("gene", GENES_WITH_A_SECOND_COPY)
    def test_knowledge_agrees_with_the_translator(self, gene):
        other = _knowledge_table(gene)
        if other is None:
            pytest.skip(f"knowledge.py has no comparable table for {gene}")
        disagreements = {
            a: (TRANSLATOR[gene][a], other[a])
            for a in sorted(set(TRANSLATOR[gene]) & set(other))
            if TRANSLATOR[gene][a] != other[a]
        }
        assert not disagreements, f"{gene}: knowledge.py disagrees: {disagreements}"


class TestNoCopyQuietlyGrowsAnAlleleTheOthersLack:
    """An allele in one table and not another is how the next divergence starts."""

    @pytest.mark.parametrize("gene", GENES_WITH_A_SECOND_COPY)
    def test_dosing_has_no_unrecorded_extra_alleles(self, gene):
        other = _dosing_table(gene)
        if other is None:
            pytest.skip(f"dosing.py has no table for {gene}")
        extra = set(other) - set(TRANSLATOR[gene]) - KNOWN_ASYMMETRIES.get((gene, "dosing"), set())
        missing = set(TRANSLATOR[gene]) - set(other) - KNOWN_ASYMMETRIES.get(
            (gene, "dosing"), set())
        assert not extra, f"{gene}: dosing.py scores {sorted(extra)}, the translator does not"
        assert not missing, (
            f"{gene}: the translator scores {sorted(missing)}, dosing.py does not -- and "
            f"dosing.py would fall back to its default for them")

    @pytest.mark.parametrize("gene", GENES_WITH_A_SECOND_COPY)
    def test_knowledge_has_no_unrecorded_extra_alleles(self, gene):
        other = _knowledge_table(gene)
        if other is None:
            pytest.skip(f"knowledge.py has no comparable table for {gene}")
        allow = KNOWN_ASYMMETRIES.get((gene, "knowledge"), set())
        extra = set(other) - set(TRANSLATOR[gene]) - allow
        missing = set(TRANSLATOR[gene]) - set(other) - allow
        assert not extra, f"{gene}: knowledge.py scores {sorted(extra)}, the translator does not"
        assert not missing, (
            f"{gene}: the translator scores {sorted(missing)}, knowledge.py does not")


class TestTheAllowlistStaysHonest:
    def test_every_recorded_asymmetry_is_still_real(self):
        """An entry that no longer describes anything should be deleted, not left to rot."""
        stale = []
        for (gene, which), alleles in KNOWN_ASYMMETRIES.items():
            other = _dosing_table(gene) if which == "dosing" else _knowledge_table(gene)
            if other is None:
                continue
            both = set(TRANSLATOR.get(gene, {})) & set(other)
            for a in alleles:
                if a in both:
                    stale.append(f"{gene}/{which}/{a}")
        assert not stale, f"KNOWN_ASYMMETRIES lists resolved entries: {stale}"

    def test_no_gene_is_exempted_wholesale(self):
        """The allowlist names alleles, never a whole gene -- otherwise it hides everything."""
        for (gene, which), alleles in KNOWN_ASYMMETRIES.items():
            assert alleles, f"{gene}/{which} has an empty allowlist entry"
            assert set(alleles) != set(TRANSLATOR.get(gene, {})), (
                f"{gene}/{which} exempts the entire table")
