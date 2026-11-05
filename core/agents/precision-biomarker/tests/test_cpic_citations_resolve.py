"""Every CPIC citation must point at the guideline it claims to.

Found by hand on 2026-11-04, reviewing precision-biomarker for check 4 (is a clinical constant
cited to the right primary source). The sweep scanner had reported this subject clean on classes
A, B and D -- correctly, because a wrong citation is not a pattern any of them looks for.

`CPIC_GUIDELINE_VERSIONS` is the module's provenance record, and `map_all()` surfaces it as
`guideline_versions`. Three of its fourteen entries cited papers with nothing to do with
pharmacogenomics, resolved against the NCBI E-utilities API:

    CYP2C19  34697867 -> abuse prevention for adults with intellectual disabilities (2022)
    UGT1A1   32233013 -> diabetes as a risk factor in COVID-19 (2020)
    CYP2C9 } 27441996 -> "Standardizing terms for clinical pharmacogenetic test results"
    VKORC1 }             -- a real CPIC paper, but a terminology consensus with no dosing in
                           it, cited by the two genes that drive warfarin mg/day output here.

A wrong citation in an audit trail is worse than no citation. A missing one is visibly missing;
this one hands a reviewer something that looks authoritative, and only a lookup shows it is a
COVID paper. Same family as the rest of this review: the machine-readable field said something
the prose beside it did not.

WHAT THIS TEST DOES AND DOES NOT ESTABLISH. It pins each PMID to a fragment of the title the
NCBI API returned on 2026-11-04, so a future transposed digit fails here instead of shipping.
It is offline and deterministic by design -- a test that called PubMed would be a network
dependency in the unit suite and would go red on an API outage rather than on a defect.

It does NOT verify that the dosing numbers in this module match those guidelines. That is the
rest of check 4, it needs a clinician reading both, and it is recorded as open in
CLINICAL_REVIEW.md. Nothing here should be read as evidence the doses are correct.
"""
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.pharmacogenomics import CPIC_GUIDELINE_VERSIONS  # noqa: E402

# PMID -> a distinctive fragment of the title NCBI returned for it on 2026-11-04.
# Verified with:
#   curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=<PMID>&retmode=json"
VERIFIED_2026_11_04 = {
    "33387367": "CYP2D6, OPRM1, and COMT Genotypes and Select Opioid Therapy",
    "35034351": "CYP2C19 Genotype and Clopidogrel Therapy: 2022 Update",
    "28198005": "Pharmacogenetics-Guided Warfarin Dosing: 2017 Update",
    "35152405": "SLCO1B1, ABCG2, and CYP2C9 genotypes and Statin-Associated",
    "30447069": "Thiopurine Dosing Based on TPMT and NUDT15 Genotypes",
    "29152729": "Dihydropyrimidine Dehydrogenase Genotype and",
    "24561393": "HLA-B Genotype and Abacavir Dosing: 2014 update",
    "26094938": "allopurinol dosing: 2015 update",
    "25801146": "CYP3A5 Genotype and Tacrolimus Dosing",
    "26417955": "UGT1A1 and Atazanavir Prescribing",
}

# The exact citations the hand review found wrong. Pinned by PMID so that re-introducing one
# fails loudly rather than quietly becoming provenance again.
REJECTED = {
    "34697867": "abuse prevention in adults with intellectual disabilities -- not a CPIC guideline",
    "32233013": "diabetes as a COVID-19 risk factor -- not a CPIC guideline",
    "27441996": "CPIC terminology standardization -- a real CPIC paper, but contains no dosing",
}


def entries_with_a_pmid():
    return {g: e for g, e in CPIC_GUIDELINE_VERSIONS.items() if e.get("pmid", "N/A") != "N/A"}


class TestNoKnownBadCitation:
    def test_the_three_wrong_pmids_are_gone(self):
        """The finding in one assertion."""
        used = {e["pmid"]: g for g, e in entries_with_a_pmid().items()}
        back = [f"  {used[p]} cites {p}: {why}" for p, why in REJECTED.items() if p in used]
        assert not back, "a citation known to be wrong is back in the provenance record:\n" + "\n".join(back)


class TestEveryCitationIsVerified:
    @pytest.mark.parametrize("gene", sorted(entries_with_a_pmid()))
    def test_pmid_was_resolved_against_ncbi(self, gene):
        pmid = CPIC_GUIDELINE_VERSIONS[gene]["pmid"]
        assert pmid in VERIFIED_2026_11_04, (
            f"{gene} cites PMID {pmid}, which nobody has resolved. Look it up against the NCBI "
            f"API, add it to VERIFIED_2026_11_04 with its real title, and only then ship it."
        )

    @pytest.mark.parametrize("gene", sorted(entries_with_a_pmid()))
    def test_the_recorded_title_matches_the_pmid(self, gene):
        """Guards the pairing, not just the PMID.

        Swapping two genes' PMIDs would leave both present and both verified; only checking the
        title the entry claims against the title that PMID really has catches it.
        """
        entry = CPIC_GUIDELINE_VERSIONS[gene]
        fragment = VERIFIED_2026_11_04[entry["pmid"]]
        claimed = re.sub(r"\s+", " ", entry["title"])
        assert fragment.lower() in claimed.lower(), (
            f"{gene}: entry says {claimed!r}, but PMID {entry['pmid']} is {fragment!r}"
        )


class TestTheRecordIsCompleteEnoughToBeProvenance:
    def test_every_gene_declares_what_its_citation_covers(self):
        for gene, entry in CPIC_GUIDELINE_VERSIONS.items():
            assert "covers" in entry, f"{gene} does not say which drugs its citation governs"
            assert "title" in entry, f"{gene} does not record the paper its PMID resolves to"

    def test_genes_without_a_guideline_say_so(self):
        """MTHFR and G6PD have no CPIC guideline. That must read as a stated absence.

        The original record used "N/A" in every field, which is indistinguishable from a gap
        nobody filled in.
        """
        for gene in ("MTHFR", "G6PD"):
            entry = CPIC_GUIDELINE_VERSIONS[gene]
            assert entry["level"] == "Informational"
            assert "no CPIC guideline" in entry["title"], (
                f"{gene} should state plainly that no guideline exists")

    def test_update_is_a_date_not_a_guess(self):
        for gene, entry in CPIC_GUIDELINE_VERSIONS.items():
            upd = entry["update"]
            if upd == "N/A":
                continue
            assert re.fullmatch(r"\d{4}-\d{2}", upd), f"{gene}: {upd!r} is not a YYYY-MM date"

    def test_multi_drug_genes_name_their_other_guidelines(self):
        """CPIC publishes per gene-drug pair; this module recommends across many drugs.

        CYP2D6 drives tamoxifen and SSRI recommendations as well as opioids, and CYP2C19 drives
        SSRIs as well as clopidogrel. One PMID cannot be provenance for all of that, so those
        two must carry the other guidelines they rely on.
        """
        for gene in ("CYP2D6", "CYP2C19"):
            also = CPIC_GUIDELINE_VERSIONS[gene].get("also", "")
            assert re.search(r"\d{8}", also), (
                f"{gene} recommends across drugs from more than one CPIC guideline; `also` must "
                f"name the others by PMID")
