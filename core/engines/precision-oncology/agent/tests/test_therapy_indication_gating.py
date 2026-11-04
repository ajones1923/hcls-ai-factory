"""A recommendation may not claim an evidence tier the drug's own label does not support.

Three defects this pins, all of which shipped with the suite green:

1. **A-8** `_identify_biomarker_therapies` took `cancer_type` and never read it, so the
   MSI, TMB, HRD/BRCA, PTEN, PD-L1 and NTRK blocks emitted unconditionally. A PD-L1
   TPS 80 PROSTATE case got "First-line monotherapy for PD-L1 TPS >= 50% NSCLC
   (KEYNOTE-024)" at evidence level A as its rank-1 recommendation; an MSI-H
   ENDOMETRIAL case got nivolumab at level A quoting a colorectal approval. The
   biomarker dicts also omitted `indicated_for_cancer_type`, so nothing downstream
   could warn. The fix had landed on two of the three identification paths.

2. **A-20** The surviving gate trusted `ACTIONABLE_TARGETS[gene]["cancer_types"]`, a
   per-GENE list of tumours where the gene matters, as if it were the drugs' approved
   indications. Being the wider list, it meant the gate never fired: EGFR L858R in
   GLIOBLASTOMA returned osimertinib, erlotinib, gefitinib, afatinib and dacomitinib at
   level A with `indicated_for_cancer_type=True`.

3. **A-21** Every `EVIDENCE_LEVEL_ORDER` lookup defaulted to 4, which is level E's rank,
   so an unrecognised level sorted as E rather than last and could win the dedupe
   against a genuine level-E entry.

Every expected indication below comes from the drug's own FDA label, retrieved from
DailyMed on 2026-11-03; the citation is on the `VERIFIED_DRUG_INDICATIONS` entry.
"""
import pytest

from src.therapy_ranker import (
    EVIDENCE_LEVEL_ORDER,
    TMB_H_LABEL_THRESHOLD,
    TMB_H_TRIAL_THRESHOLD,
    UNKNOWN_EVIDENCE_RANK,
    TherapyRanker,
    evidence_rank,
)


@pytest.fixture
def ranker():
    return TherapyRanker(collection_manager=None, embedder=None, knowledge=None)


def _find(therapies, drug):
    return next((t for t in therapies if t["drug_name"].lower() == drug.lower()), None)


class TestEgfrTkisAreNotLevelAInGlioblastoma:
    """A-20. EGFR TKI labels are NSCLC-only; EGFR's gene list names glioblastoma."""

    @pytest.mark.parametrize("drug", [
        "osimertinib", "erlotinib", "gefitinib", "afatinib", "dacomitinib", "amivantamab",
    ])
    def test_not_level_a_and_not_indicated(self, ranker, drug):
        out = ranker._identify_variant_therapies("EGFR", "L858R", "glioblastoma")
        t = _find(out, drug)
        assert t is not None, "the candidate should still be reported, just not as approved"
        assert t["indicated_for_cancer_type"] is False
        assert t["evidence_level"] != "A"

    def test_note_names_the_label(self, ranker):
        out = ranker._identify_variant_therapies("EGFR", "L858R", "glioblastoma")
        note = _find(out, "osimertinib")["review_note"]
        assert "not approved in glioblastoma" in note
        assert "TAGRISSO" in note

    def test_nsclc_keeps_level_a(self, ranker):
        """The gate must not cost the drugs their real indication."""
        out = ranker._identify_variant_therapies("EGFR", "L858R", "NSCLC")
        for drug in ("osimertinib", "gefitinib", "afatinib", "dacomitinib"):
            t = _find(out, drug)
            assert t["evidence_level"] == "A"
            assert t["indicated_for_cancer_type"] is True

    def test_colorectal_egfr_tkis_are_gated(self, ranker):
        """EGFR's gene list names colorectal, but for antibodies on RAS-wild-type status.

        No EGFR TKI has a colorectal indication; the gene list conflated two different
        uses of the word "EGFR".
        """
        out = ranker._identify_variant_therapies("EGFR", "L858R", "colorectal")
        assert all(t["indicated_for_cancer_type"] is False for t in out)
        assert all(t["evidence_level"] != "A" for t in out)


class TestPerDrugRatherThanPerGene:
    """The gate's unit is the (drug, driver) pair, not the gene."""

    def test_lorlatinib_has_no_ros1_indication(self, ranker):
        """LORBRENA is ALK-positive NSCLC only; ROS1's gene list carries lorlatinib."""
        out = ranker._identify_variant_therapies("ROS1", "ROS1 fusion", "NSCLC")
        lorla = _find(out, "lorlatinib")
        assert lorla["indicated_for_cancer_type"] is False
        assert "no approved indication" in lorla["review_note"].lower()
        # while the drugs that do have a ROS1 NSCLC indication keep it
        assert _find(out, "crizotinib")["indicated_for_cancer_type"] is True
        assert _find(out, "entrectinib")["indicated_for_cancer_type"] is True

    def test_same_drug_differs_by_driver(self, ranker):
        """entrectinib: ROS1 -> NSCLC only; NTRK -> any site."""
        ros1_in_prostate = _find(
            ranker._identify_variant_therapies("ROS1", "ROS1 fusion", "prostate"), "entrectinib"
        )
        ntrk_in_prostate = _find(
            ranker._identify_biomarker_therapies({"NTRK": "fusion"}, "prostate"), "entrectinib"
        )
        assert ros1_in_prostate["indicated_for_cancer_type"] is False
        assert ntrk_in_prostate["indicated_for_cancer_type"] is True
        assert ntrk_in_prostate["evidence_level"] == "A"

    def test_limitation_of_use_is_honoured(self, ranker):
        """TAFINLAR's tissue-agnostic indication explicitly excludes colorectal cancer."""
        out = ranker._identify_variant_therapies("BRAF", "V600E", "colorectal")
        dab = _find(out, "dabrafenib")
        assert dab["indicated_for_cancer_type"] is False
        assert "Limitation of Use" in dab["review_note"]
        # BRAFTOVI does have a BRAF V600E mCRC indication, so it must survive
        assert _find(out, "encorafenib")["indicated_for_cancer_type"] is True
        assert _find(out, "encorafenib")["evidence_level"] == "A"

    def test_unverified_drug_is_declared_not_asserted(self, ranker):
        """A drug with no retrieved label must say so rather than imply an approval."""
        out = ranker._identify_variant_therapies("MET", "exon 14 skipping", "NSCLC")
        tepotinib = _find(out, "tepotinib")
        assert tepotinib is not None
        assert tepotinib["indication_unverified"] is True
        assert "NOT been verified against a label" in tepotinib["review_note"]
        # capmatinib's label WAS read, so it is not flagged
        assert _find(out, "capmatinib")["indication_unverified"] is False


class TestBiomarkerPathReadsCancerType:
    """A-8. The third identification path ignored the tumour entirely."""

    def test_pdl1_tps_outside_nsclc_is_not_level_a(self, ranker):
        out = ranker._identify_biomarker_therapies(
            {"PD-L1_TPS": 80}, "prostate adenocarcinoma"
        )
        pembro = _find(out, "pembrolizumab")
        assert pembro["evidence_level"] != "A"
        assert pembro["indicated_for_cancer_type"] is False
        assert "Tumor Proportion Score" in pembro["review_note"]

    def test_pdl1_tps_in_nsclc_keeps_level_a(self, ranker):
        out = ranker._identify_biomarker_therapies({"PD-L1_TPS": 80}, "NSCLC")
        pembro = _find(out, "pembrolizumab")
        assert pembro["evidence_level"] == "A"
        assert pembro["indicated_for_cancer_type"] is True

    def test_pdl1_with_no_cancer_type_is_not_level_a(self, ranker):
        out = ranker._identify_biomarker_therapies({"PD-L1_TPS": 80}, "")
        assert _find(out, "pembrolizumab")["evidence_level"] != "A"

    def test_msi_h_nivolumab_is_colorectal_only(self, ranker):
        endometrial = ranker._identify_biomarker_therapies({"MSI": "MSI-H"}, "endometrial carcinoma")
        colorectal = ranker._identify_biomarker_therapies({"MSI": "MSI-H"}, "colorectal")
        assert _find(endometrial, "nivolumab")["indicated_for_cancer_type"] is False
        assert _find(endometrial, "nivolumab")["evidence_level"] != "A"
        assert _find(colorectal, "nivolumab")["evidence_level"] == "A"

    def test_msi_h_tissue_agnostic_drugs_keep_level_a_anywhere(self, ranker):
        out = ranker._identify_biomarker_therapies({"MSI": "MSI-H"}, "endometrial carcinoma")
        for drug in ("pembrolizumab", "dostarlimab"):
            assert _find(out, drug)["evidence_level"] == "A"
            assert _find(out, drug)["indicated_for_cancer_type"] is True

    def test_parp_inhibitors_are_gated_per_drug(self, ranker):
        """RUBRACA has no breast indication; TALZENNA has no ovarian indication."""
        ovarian = ranker._identify_biomarker_therapies({"BRCA": "BRCA2"}, "ovarian")
        assert _find(ovarian, "olaparib")["evidence_level"] == "A"
        assert _find(ovarian, "rucaparib")["evidence_level"] == "A"
        assert _find(ovarian, "niraparib")["evidence_level"] == "A"
        assert _find(ovarian, "talazoparib")["indicated_for_cancer_type"] is False

        breast = ranker._identify_biomarker_therapies({"BRCA": "BRCA1"}, "breast")
        assert _find(breast, "olaparib")["evidence_level"] == "A"
        assert _find(breast, "talazoparib")["evidence_level"] == "A"
        assert _find(breast, "rucaparib")["indicated_for_cancer_type"] is False
        assert _find(breast, "niraparib")["indicated_for_cancer_type"] is False

    def test_brca_in_glioblastoma_yields_no_level_a_parp(self, ranker):
        out = ranker._identify_biomarker_therapies({"BRCA": "BRCA2"}, "glioblastoma")
        parps = [t for t in out if t["category"] == "PARP inhibitor"]
        assert parps
        assert all(t["evidence_level"] != "A" for t in parps)
        assert all(t["indicated_for_cancer_type"] is False for t in parps)

    def test_pten_loss_has_no_alpelisib_indication(self, ranker):
        """PIQRAY selects on a PIK3CA mutation, not on PTEN status."""
        out = ranker._identify_biomarker_therapies({"PTEN": "LOSS"}, "prostate")
        alpelisib = _find(out, "alpelisib")
        assert alpelisib["indicated_for_cancer_type"] is False
        assert "PIK3CA" in alpelisib["review_note"]

    def test_every_biomarker_entry_carries_the_indication_flag(self, ranker):
        """The biomarker dicts omitted the key, so no consumer could warn."""
        out = ranker._identify_biomarker_therapies(
            {"MSI": "MSI-H", "TMB": 42.0, "BRCA": "BRCA2", "PTEN": "LOSS",
             "PD-L1_TPS": 80, "NTRK": "fusion"},
            "prostate adenocarcinoma",
        )
        assert out
        for t in out:
            assert "indicated_for_cancer_type" in t
            assert "approved_cancer_types" in t
            assert "indication_unverified" in t

    def test_full_rank_does_not_put_an_off_label_level_a_first(self, ranker):
        """The audit's case: mCRPC with PD-L1 TPS 80 and a somatic BRCA2 mutation."""
        out = ranker.rank_therapies(
            cancer_type="prostate adenocarcinoma",
            variants=[],
            biomarkers={"PD-L1_TPS": 80, "BRCA": "BRCA2"},
        )
        pembro = _find(out, "pembrolizumab")
        assert pembro["evidence_level"] != "A"
        assert pembro["indicated_for_cancer_type"] is False
        # olaparib and rucaparib DO have prostate indications and should lead
        assert _find(out, "olaparib")["evidence_level"] == "A"
        assert _find(out, "rucaparib")["evidence_level"] == "A"
        assert _find(out, "niraparib")["indicated_for_cancer_type"] is False


class TestTmbThreshold:
    """The label defines TMB-H as >=13 mut/Mb; the code fired at >=10 claiming approval."""

    def test_label_threshold_is_thirteen(self):
        assert TMB_H_LABEL_THRESHOLD == 13.0
        assert TMB_H_TRIAL_THRESHOLD == 10.0

    @pytest.mark.parametrize("tmb", [13.0, 15.0, 42.0])
    def test_at_or_above_label_threshold_is_level_a(self, ranker, tmb):
        out = ranker._identify_biomarker_therapies({"TMB": tmb}, "NSCLC")
        pembro = _find(out, "pembrolizumab")
        assert pembro["evidence_level"] == "A"
        assert ">=13 mut/Mb" in pembro["guideline_recommendation"]

    @pytest.mark.parametrize("tmb", [10.0, 11.5, 12.9])
    def test_between_trial_and_label_threshold_is_not_approved(self, ranker, tmb):
        out = ranker._identify_biomarker_therapies({"TMB": tmb}, "NSCLC")
        pembro = _find(out, "pembrolizumab")
        assert pembro is not None, "the candidate is still reported, not dropped"
        assert pembro["evidence_level"] != "A"
        assert "BELOW the labelled TMB-H definition" in pembro["guideline_recommendation"]

    def test_below_trial_threshold_emits_nothing(self, ranker):
        out = ranker._identify_biomarker_therapies({"TMB": 6.0}, "NSCLC")
        assert _find(out, "pembrolizumab") is None

    def test_atezolizumab_has_no_tmb_indication(self, ranker):
        out = ranker._identify_biomarker_therapies({"TMB": 42.0}, "NSCLC")
        atezo = _find(out, "atezolizumab")
        assert atezo["indicated_for_cancer_type"] is False
        assert "does not mention tumour" in atezo["review_note"]


class TestUnknownEvidenceLevelSortsLast:
    """A-21. Every lookup defaulted to 4 — level E's rank — not last."""

    def test_unknown_rank_is_weaker_than_vus(self):
        assert UNKNOWN_EVIDENCE_RANK > EVIDENCE_LEVEL_ORDER["VUS"]
        assert UNKNOWN_EVIDENCE_RANK > EVIDENCE_LEVEL_ORDER["E"]

    @pytest.mark.parametrize("level", [None, "", "I", "IA", "II", "R1", "unknown"])
    def test_unrecognised_levels_rank_last(self, level):
        assert evidence_rank(level) == UNKNOWN_EVIDENCE_RANK

    @pytest.mark.parametrize("level,rank", [
        ("A", 0), ("B", 1), ("C", 2), ("D", 3), ("E", 4), ("VUS", 5),
        ("a", 0), (" b ", 1), ("vus", 5),
    ])
    def test_recognised_levels_keep_their_rank(self, level, rank):
        assert evidence_rank(level) == rank

    def test_level_e_beats_an_unrecognised_level_in_the_ranking(self, ranker):
        ranked = ranker._assign_final_ranks([
            {"drug_name": "amp_tier_one", "evidence_level": "I"},
            {"drug_name": "genuine_level_e", "evidence_level": "E"},
        ])
        assert [t["drug_name"] for t in ranked] == ["genuine_level_e", "amp_tier_one"]

    def test_dedupe_prefers_level_e_over_an_unrecognised_level(self, ranker):
        """The dedupe compared ranks, so an "I" tied with E and could win on arrival order."""
        assert evidence_rank("E") < evidence_rank("I")


class TestExportSurfacesTheIndicationMismatch:
    """A flag nothing renders is not a warning.

    The markdown export's Notes column read `notes`/`rationale`, neither of which the
    ranker emits, so an off-label recommendation reached the tumour board showing only
    its evidence tier.
    """

    def test_markdown_marks_an_off_indication_therapy(self):
        from src.export import export_markdown

        md = export_markdown({
            "case_id": "T-1",
            "cancer_type": "glioblastoma",
            "therapies": [{
                "name": "osimertinib",
                "evidence_level": "C",
                "indicated_for_cancer_type": False,
                "review_note": "osimertinib is not approved in glioblastoma.",
            }],
        })
        assert "NOT INDICATED for this tumour type" in md
        assert "not approved in glioblastoma" in md

    def test_markdown_marks_an_unverified_indication(self):
        from src.export import export_markdown

        md = export_markdown({
            "case_id": "T-2",
            "cancer_type": "NSCLC",
            "therapies": [{
                "name": "tepotinib",
                "evidence_level": "A",
                "indicated_for_cancer_type": True,
                "indication_unverified": True,
            }],
        })
        assert "Indication unverified against a label" in md

    def test_fhir_medicationrequest_carries_the_refusal(self):
        from src.export import export_fhir_r4

        bundle = export_fhir_r4({
            "case_id": "T-3",
            "cancer_type": "prostate adenocarcinoma",
            "therapies": [{
                "name": "pembrolizumab",
                "evidence_level": "C",
                "indicated_for_cancer_type": False,
                "review_note": "TPS indications are NSCLC.",
            }],
        }, patient_id="PT-3")
        meds = [e["resource"] for e in bundle["entry"]
                if e["resource"]["resourceType"] == "MedicationRequest"]
        assert meds
        notes = " ".join(n["text"] for n in meds[0]["note"])
        assert "NOT within the approved indications" in notes
