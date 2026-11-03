"""Therapy ranking must not borrow evidence from a different variant or a different tumour.

Three defects this pins:

1. `_identify_variant_therapies` handed an UNMATCHED variant the gene's top evidence level.
   EGFR carries evidence_level "A" and lists T790M, C797S and exon-20 insertion among its
   key_variants, so an osimertinib-resistant C797S case was offered osimertinib and erlotinib
   at level A — the level L858R earned.

2. Neither the variant nor the combination path read `cancer_type`, though both took it as an
   argument. A BRAF V600E melanoma case was offered "encorafenib + cetuximab" — BEACON,
   colorectal only — at level A.

3. Combination regimens were appended AFTER the resistance, contraindication and evidence
   steps with `resistance_flag=False` and `supporting_evidence=[]` hard-coded, so a patient
   who had already failed dabrafenib + trametinib was re-offered it unflagged.
"""
import pytest

from src.therapy_ranker import TherapyRanker


@pytest.fixture
def ranker():
    """A ranker with no vector backend; retrieval is not what these tests exercise.

    _fetch_supporting_evidence degrades to an empty list without a collection manager, which
    is the condition these tests want: gating must hold with or without retrieval.
    """
    return TherapyRanker(collection_manager=None, embedder=None, knowledge=None)


def _names(therapies):
    return [t["drug_name"].lower() for t in therapies]


def _find(therapies, drug):
    return next((t for t in therapies if t["drug_name"].lower() == drug.lower()), None)


class TestEvidenceIsNotBorrowedFromAnotherVariant:
    def test_unmatched_variant_does_not_get_level_a(self, ranker):
        """P848L is not among EGFR's key_variants, so it earns none of their evidence."""
        out = ranker._identify_variant_therapies("EGFR", "P848L", "NSCLC")
        assert out, "EGFR should still yield gene-level candidates"
        assert all(t["evidence_level"] != "A" for t in out), (
            "an unlisted variant must not inherit EGFR's level A"
        )

    def test_unmatched_variant_is_marked_as_gene_level_inference(self, ranker):
        out = ranker._identify_variant_therapies("EGFR", "P848L", "NSCLC")
        assert all(t["gene_level_inference"] for t in out)
        assert all(t["variant_matched"] is False for t in out)

    def test_matched_variant_keeps_its_evidence_level(self, ranker):
        out = ranker._identify_variant_therapies("EGFR", "L858R", "NSCLC")
        assert any(t["evidence_level"] == "A" for t in out)
        assert all(t["variant_matched"] for t in out)

    def test_resistance_mutation_is_flagged(self, ranker):
        """C797S is listed BOTH as a key variant and as a resistance mechanism.

        It therefore stays "matched" and keeps its evidence level, because which drug a
        resistance mutation defeats is drug-specific and is not encoded in this knowledge
        base. What must not happen is it passing through silently.
        """
        out = ranker._identify_variant_therapies("EGFR", "C797S", "NSCLC")
        assert all(t["variant_is_resistance_mutation"] for t in out)
        assert all("resistance mechanism" in t["review_note"] for t in out)


class TestCancerTypeGating:
    def test_off_indication_tumour_is_downgraded(self, ranker):
        """EGFR's indications do not include pancreatic cancer."""
        on = ranker._identify_variant_therapies("EGFR", "L858R", "NSCLC")
        off = ranker._identify_variant_therapies("EGFR", "L858R", "pancreatic")
        assert all(t["indicated_for_cancer_type"] for t in on)
        assert all(t["indicated_for_cancer_type"] is False for t in off)
        assert _find(off, "osimertinib")["evidence_level"] != _find(on, "osimertinib")["evidence_level"]

    def test_synonym_still_matches(self, ranker):
        """"CRC" must resolve to the "colorectal" spelling used in the knowledge base."""
        assert ranker._cancer_type_applies("CRC", ["colorectal"]) is True
        assert ranker._cancer_type_applies("metastatic colorectal carcinoma", ["colorectal"]) is True

    def test_tissue_agnostic_matches_anything(self, ranker):
        assert ranker._cancer_type_applies("anything", ["tissue-agnostic"]) is True

    def test_unstated_indication_is_not_treated_as_approved(self, ranker):
        assert ranker._cancer_type_applies("NSCLC", []) is False


class TestCombinationRegimens:
    def test_crc_only_combo_is_not_offered_for_melanoma(self, ranker):
        out = ranker.rank_therapies(
            cancer_type="melanoma",
            variants=[{"gene": "BRAF", "variant": "V600E"}],
            biomarkers={},
        )
        assert "encorafenib + cetuximab" not in _names(out), (
            "BEACON is colorectal-only; it must not appear for melanoma"
        )

    def test_melanoma_combo_is_still_offered_for_melanoma(self, ranker):
        out = ranker.rank_therapies(
            cancer_type="melanoma",
            variants=[{"gene": "BRAF", "variant": "V600E"}],
            biomarkers={},
        )
        assert "dabrafenib + trametinib" in _names(out)

    def test_previously_failed_combo_is_flagged_not_silently_reoffered(self, ranker):
        out = ranker.rank_therapies(
            cancer_type="melanoma",
            variants=[{"gene": "BRAF", "variant": "V600E"}],
            biomarkers={},
            prior_therapies=["dabrafenib", "trametinib"],
        )
        combo = _find(out, "dabrafenib + trametinib")
        if combo is not None:
            assert combo["resistance_flag"] or combo["contraindication_flag"], (
                "a combo whose components already failed must carry a flag"
            )

    def test_combos_reach_the_evidence_step(self, ranker):
        """Combos were appended after evidence retrieval, so the key was always absent."""
        out = ranker.rank_therapies(
            cancer_type="melanoma",
            variants=[{"gene": "BRAF", "variant": "V600E"}],
            biomarkers={},
        )
        for t in out:
            assert "supporting_evidence" in t
            assert "resistance_flag" in t


class TestBiomarkerTypeRobustness:
    @pytest.mark.parametrize("value", [None, True, 1])
    def test_non_string_biomarker_does_not_raise(self, ranker, value):
        """`biomarkers.get("MSI", "").upper()` raised AttributeError on a None or bool."""
        ranker.rank_therapies(
            cancer_type="colorectal",
            variants=[],
            biomarkers={"MSI": value, "BRCA": value, "PTEN": value, "NTRK": value},
        )
