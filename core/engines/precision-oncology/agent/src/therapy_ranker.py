"""
Oncology Intelligence Agent - Therapy Ranker
Evidence-based therapy ranking with citations and resistance awareness.

Ranks targeted therapies by evidence level, cross-references prior therapy
history to detect resistance patterns, and retrieves supporting evidence
from RAG collections. Designed for Molecular Tumor Board decision support.

Author: Adam Jones
Date: February 2026
License: Apache 2.0
"""

import logging
from typing import Any, Dict, List, Optional

from src.models import CaseSnapshot
from src.knowledge import (
    ACTIONABLE_TARGETS,
    THERAPY_MAP,
    RESISTANCE_MAP,
    BIOMARKER_PANELS,
)

logger = logging.getLogger(__name__)


# Evidence level ordering: lower value = stronger evidence
EVIDENCE_LEVEL_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "VUS": 5}

# An evidence level this module does not recognise must sort LAST, not into the
# middle of the ladder. Every lookup used to default to 4, which is level E's
# rank: a therapy whose level was missing, or written in another vocabulary
# (AMP/ASCO/CAP "I"/"II"/"IA", CIViC "R1"), ranked as though it were level E and
# could win the dedupe against a genuine level-E entry.
UNKNOWN_EVIDENCE_RANK = max(EVIDENCE_LEVEL_ORDER.values()) + 1

_warned_evidence_levels: set = set()


def evidence_rank(level: Optional[str]) -> int:
    """Rank an evidence level for sorting; an unrecognised level sorts last."""
    key = (level or "").strip().upper()
    if key in EVIDENCE_LEVEL_ORDER:
        return EVIDENCE_LEVEL_ORDER[key]
    if key not in _warned_evidence_levels:
        _warned_evidence_levels.add(key)
        logger.warning(
            "Unrecognised evidence level %r: ranking it last (weaker than E/VUS)",
            level,
        )
    return UNKNOWN_EVIDENCE_RANK


# ---------------------------------------------------------------------------
# Per-DRUG (and per-driver) approved indications, from the drug's own FDA label.
#
# WHY THIS EXISTS: `ACTIONABLE_TARGETS[gene]["cancer_types"]` is a per-GENE list
# — "tumours in which this gene matters" — and the gate was using it as
# "tumours in which these drugs are approved". The two are not the same list,
# and the gene list is the wider of the two, so the gate never fired.
# `ACTIONABLE_TARGETS["EGFR"]["cancer_types"]` is
# ['NSCLC', 'head and neck', 'colorectal', 'glioblastoma'] at evidence_level A,
# so an EGFR L858R GLIOBLASTOMA returned osimertinib, erlotinib, gefitinib,
# afatinib and dacomitinib at level A with indicated_for_cancer_type=True —
# none of which is approved in glioblastoma.
#
# Keys are (drug lower-cased, driver) where driver is the gene symbol or the
# biomarker name, because an indication is a (drug, biomarker, tumour) triple:
# entrectinib is approved in ROS1-positive NSCLC *only*, and in NTRK
# fusion-positive solid tumours of *any* site.
#
# Each entry carries the label it was read from. Every one below was retrieved
# on 2026-11-03 from DailyMed (dailymed.nlm.nih.gov/dailymed/services/v2/spls/
# <setid>.xml), section 1 INDICATIONS AND USAGE. Nothing is entered here from
# memory: a drug with no retrieved label is simply absent, and the gate reports
# that as unverified rather than guessing.
#
#   tumors:   tumour terms the label's INDICATIONS AND USAGE covers for this
#             driver; ["tissue-agnostic"] for a site-agnostic indication; []
#             when the label carries NO indication for this driver.
#   excluded: tumours the label explicitly carves out (Limitations of Use).
VERIFIED_DRUG_INDICATIONS: Dict[tuple, Dict[str, Any]] = {
    # --- EGFR ---------------------------------------------------------------
    ("osimertinib", "EGFR"): {
        "tumors": ["NSCLC"],
        "label": "TAGRISSO, AstraZeneca, DailyMed setid 5e81b4a7-b971-45e1-9c31-29cea8c87ce7 "
                 "(published 2026-09-25) §1.1-1.5: NSCLC only",
    },
    ("gefitinib", "EGFR"): {
        "tumors": ["NSCLC"],
        "label": "Gefitinib tablets, DailyMed setid 3dd135f0-5db1-4236-9756-04533b66dc9d §1: "
                 "metastatic NSCLC with exon 19 del / exon 21 L858R only",
    },
    ("afatinib", "EGFR"): {
        "tumors": ["NSCLC"],
        "label": "Afatinib tablets, DailyMed setid ec5b9e00-3eab-4dbd-9bb3-75865b1af485 "
                 "(published 2026-09-04) §1.1-1.2: NSCLC only",
    },
    ("dacomitinib", "EGFR"): {
        "tumors": ["NSCLC"],
        "label": "VIZIMPRO, Pfizer, DailyMed setid 4ab27d2f-e385-4e9c-b324-fa69c10b855a "
                 "(published 2026-06-25) §1: metastatic NSCLC only",
    },
    ("erlotinib", "EGFR"): {
        "tumors": ["NSCLC"],
        "excluded": [],
        "label": "TARCEVA, DailyMed setid 5783976b-0862-44e5-9b19-0f41d236d5c3 §1.1-1.2: NSCLC, "
                 "and pancreatic cancer with gemcitabine. The pancreatic indication is NOT "
                 "EGFR-mutation-selected, so it is not listed as an EGFR-driver indication here.",
    },
    ("amivantamab", "EGFR"): {
        "tumors": ["NSCLC"],
        "label": "RYBREVANT FASPRO, Janssen, DailyMed setid "
                 "9e58b045-e352-4f62-99bf-77fae1bebc69 (published 2026-10-30) §1.1-1.4: NSCLC only",
    },
    # --- ALK ----------------------------------------------------------------
    ("alectinib", "ALK"): {
        "tumors": ["NSCLC"],
        "label": "ALECENSA, Genentech, DailyMed setid (published 2026-08-27) §1.1-1.2: ALK+ NSCLC only",
    },
    ("brigatinib", "ALK"): {
        "tumors": ["NSCLC"],
        "label": "ALUNBRIG, Takeda, DailyMed (published 2026-01-19) §1: ALK+ metastatic NSCLC only",
    },
    ("ceritinib", "ALK"): {
        "tumors": ["NSCLC"],
        "label": "ZYKADIA, Novartis, DailyMed (published 2026-10-30) §1: ALK+ metastatic NSCLC only",
    },
    ("lorlatinib", "ALK"): {
        "tumors": ["NSCLC"],
        "label": "LORBRENA, Pfizer, DailyMed setid 2b34d62d-e02a-4af3-bc0d-1571dd4ee76d "
                 "(published 2026-06-29) §1: ALK+ metastatic NSCLC only",
    },
    ("crizotinib", "ALK"): {
        "tumors": ["NSCLC", "anaplastic large cell lymphoma", "inflammatory myofibroblastic tumor"],
        "label": "XALKORI, Pfizer, DailyMed setid 2a51b0de-47d6-455e-a94c-d2c737b04ff7 "
                 "(published 2026-11-02) §1.1-1.3",
    },
    # --- ROS1 ---------------------------------------------------------------
    ("crizotinib", "ROS1"): {
        "tumors": ["NSCLC"],
        "label": "XALKORI §1.1: ROS1-positive metastatic NSCLC (the ALCL and IMT indications are "
                 "ALK-positive, not ROS1)",
    },
    ("entrectinib", "ROS1"): {
        "tumors": ["NSCLC"],
        "label": "ROZLYTREK, Genentech, DailyMed setid c7c71b0c-2549-4495-86b6-c2807fa54908 "
                 "(published 2026-10-16) §1.1: ROS1-positive metastatic NSCLC",
    },
    ("repotrectinib", "ROS1"): {
        "tumors": ["NSCLC"],
        "label": "AUGTYRO, E.R. Squibb, DailyMed setid fb526827-40ba-4462-94cf-179ae3b0cb8a "
                 "(published 2026-08-11) §1.1: ROS1-positive NSCLC",
    },
    ("lorlatinib", "ROS1"): {
        "tumors": [],
        "label": "LORBRENA §1 carries an ALK-positive NSCLC indication ONLY — no ROS1 indication "
                 "at any site (DailyMed, published 2026-06-29)",
    },
    # --- BRAF ---------------------------------------------------------------
    ("vemurafenib", "BRAF"): {
        "tumors": ["melanoma", "Erdheim-Chester"],
        "label": "ZELBORAF, Genentech, DailyMed (published 2026-09-04) §1.1-1.2",
    },
    ("dabrafenib", "BRAF"): {
        "tumors": ["tissue-agnostic"],
        "excluded": ["colorectal"],
        "label": "TAFINLAR, Novartis, DailyMed (published 2026-05-21) §1: melanoma, NSCLC, "
                 "anaplastic thyroid, low-grade glioma, and (with trametinib) BRAF V600E solid "
                 "tumours; Limitation of Use excludes colorectal cancer for known intrinsic "
                 "resistance to BRAF inhibition",
    },
    ("encorafenib", "BRAF"): {
        "tumors": ["melanoma", "colorectal", "NSCLC"],
        "label": "BRAFTOVI, Array BioPharma, DailyMed (published 2026-03-02) §1.1-1.3",
    },
    # --- KRAS G12C ----------------------------------------------------------
    ("sotorasib", "KRAS"): {
        "tumors": ["NSCLC", "colorectal"],
        "label": "LUMAKRAS, Amgen, DailyMed (published 2025-01-30) §1.1-1.2: KRAS G12C NSCLC, and "
                 "KRAS G12C mCRC with panitumumab",
    },
    ("adagrasib", "KRAS"): {
        "tumors": ["NSCLC"],
        "label": "KRAZATI, Mirati, DailyMed (published 2026-10-23) §1: KRAS G12C locally advanced "
                 "or metastatic NSCLC; no colorectal indication in the retrieved label",
    },
    # --- HER2 ---------------------------------------------------------------
    ("trastuzumab deruxtecan", "HER2"): {
        "tumors": ["tissue-agnostic"],
        "label": "ENHERTU, Daiichi Sankyo, DailyMed (published 2026-05-28) §1: breast, NSCLC, "
                 "gastric, and HER2 IHC3+ solid tumours after prior systemic treatment",
    },
    ("tucatinib", "HER2"): {
        "tumors": ["breast", "colorectal"],
        "label": "TUKYSA, Seagen, DailyMed (published 2026-10-19) §1.1: HER2+ breast cancer, and "
                 "RAS wild-type HER2+ colorectal cancer with trastuzumab",
    },
    # --- NTRK (tissue-agnostic) ---------------------------------------------
    ("larotrectinib", "NTRK"): {
        "tumors": ["tissue-agnostic"],
        "label": "VITRAKVI, Bayer, DailyMed (published 2026-05-25) §1: NTRK fusion-positive solid "
                 "tumours, any site",
    },
    ("entrectinib", "NTRK"): {
        "tumors": ["tissue-agnostic"],
        "label": "ROZLYTREK §1.2: NTRK fusion-positive solid tumours, any site",
    },
    ("repotrectinib", "NTRK"): {
        "tumors": ["tissue-agnostic"],
        "label": "AUGTYRO §1.2: NTRK fusion-positive solid tumours, any site (12 years and older)",
    },
    # --- RET ----------------------------------------------------------------
    ("selpercatinib", "RET"): {
        "tumors": ["tissue-agnostic"],
        "label": "RETEVMO, Eli Lilly, DailyMed (published 2026-09-14) §1.1-1.4: RET fusion NSCLC, "
                 "RET-mutant medullary thyroid, RET fusion thyroid, and RET fusion solid tumours",
    },
    # --- MET ----------------------------------------------------------------
    ("capmatinib", "MET"): {
        "tumors": ["NSCLC"],
        "label": "TABRECTA, Novartis, DailyMed (published 2026-10-28) §1: MET exon 14 skipping "
                 "metastatic NSCLC only",
    },
    # --- PIK3CA / PTEN ------------------------------------------------------
    ("alpelisib", "PIK3CA"): {
        "tumors": ["breast"],
        "label": "PIQRAY, Novartis, DailyMed (published 2025-12-16) §1: PIK3CA-mutated HR+/HER2- "
                 "advanced or metastatic breast cancer with fulvestrant. (VIJOICE, the same "
                 "molecule, is indicated for PIK3CA-Related Overgrowth Spectrum, not cancer.)",
    },
    ("alpelisib", "PTEN"): {
        "tumors": [],
        "label": "PIQRAY §1 is selected on a PIK3CA mutation, not on PTEN status; no PTEN-loss "
                 "indication exists in the retrieved label (DailyMed, published 2025-12-16)",
    },
    # --- IDH1 ---------------------------------------------------------------
    ("ivosidenib", "IDH1"): {
        "tumors": ["AML", "myelodysplastic syndromes", "cholangiocarcinoma"],
        "label": "TIBSOVO, Servier, DailyMed (published 2025-12-22) §1.1-1.4",
    },
    ("vorasidenib", "IDH1"): {
        "tumors": ["astrocytoma", "oligodendroglioma", "glioma"],
        "label": "VORANIGO, Servier, DailyMed (published 2026-09-17) §1: Grade 2 astrocytoma or "
                 "oligodendroglioma with a susceptible IDH1/IDH2 mutation",
    },
    # --- PARP inhibitors, BRCA/HRD ------------------------------------------
    ("olaparib", "BRCA"): {
        "tumors": ["ovarian", "fallopian tube", "primary peritoneal", "breast", "prostate",
                   "pancreatic"],
        "label": "LYNPARZA, AstraZeneca, DailyMed (published 2026) §1.1-1.8",
    },
    ("rucaparib", "BRCA"): {
        "tumors": ["ovarian", "fallopian tube", "primary peritoneal", "prostate"],
        "label": "RUBRACA, DailyMed §1.1-1.2: BRCA-mutated recurrent ovarian cancer, and "
                 "BRCA-mutated mCRPC. No breast indication.",
    },
    ("niraparib", "BRCA"): {
        "tumors": ["ovarian", "fallopian tube", "primary peritoneal"],
        "label": "ZEJULA, DailyMed §1.1-1.2: advanced/recurrent epithelial ovarian, fallopian "
                 "tube or primary peritoneal cancer only",
    },
    ("talazoparib", "BRCA"): {
        "tumors": ["breast", "prostate"],
        "label": "TALZENNA, Pfizer, DailyMed §1.1-1.2: gBRCAm HER2-negative breast cancer, and "
                 "HRR-mutated metastatic prostate cancer with enzalutamide",
    },
    # --- Immunotherapy biomarkers -------------------------------------------
    ("pembrolizumab", "MSI-H"): {
        "tumors": ["tissue-agnostic"],
        "label": "KEYTRUDA, Merck, DailyMed setid 9333c79b-d487-4538-a9f0-71b91a02b287 "
                 "(published 2026-11-02): MSI-H/dMMR solid tumours, any site",
    },
    ("nivolumab", "MSI-H"): {
        "tumors": ["colorectal"],
        "label": "OPDIVO, Bristol-Myers Squibb, DailyMed (published 2026): the MSI-H/dMMR "
                 "indications are colorectal cancer only — there is no tissue-agnostic MSI-H "
                 "indication in the retrieved label",
    },
    ("dostarlimab", "MSI-H"): {
        "tumors": ["tissue-agnostic"],
        "label": "JEMPERLI, GSK, DailyMed (published 2026) §1.3: dMMR recurrent or advanced solid "
                 "tumours, any site (§1.1 endometrial and §1.2 MSI-H/dMMR rectal are "
                 "site-specific)",
    },
    ("pembrolizumab", "TMB-H"): {
        "tumors": ["tissue-agnostic"],
        "label": "KEYTRUDA setid 9333c79b-d487-4538-a9f0-71b91a02b287 (published 2026-11-02) "
                 "§1.16: TMB-H [>=13 mut/Mb] solid tumours, any site",
    },
    ("atezolizumab", "TMB-H"): {
        "tumors": [],
        "label": "TECENTRIQ, Genentech, DailyMed setid 6fa682c9-a312-4932-9831-f286908660ee "
                 "(published 2026-10-16): INDICATIONS AND USAGE does not mention tumour "
                 "mutational burden at all — there is no TMB-H indication at any site",
    },
    ("pembrolizumab", "PD-L1_TPS"): {
        "tumors": ["NSCLC"],
        "label": "KEYTRUDA setid 9333c79b-d487-4538-a9f0-71b91a02b287 (published 2026-11-02): the "
                 "Tumor Proportion Score (TPS) indications are NSCLC; other tumours are selected "
                 "on Combined Positive Score (CPS), not TPS",
    },
}

# Pembrolizumab's TMB-H indication threshold, from the label cited above. The
# code used >= 10 mut/Mb and said "FDA-approved for TMB-H (>=10 mut/Mb)".
TMB_H_LABEL_THRESHOLD = 13.0
# KEYNOTE-158's TMB-H cohort cut-off, and the figure BIOMARKER_PANELS still
# carries. Between the two the recommendation is reported, but not as approved.
TMB_H_TRIAL_THRESHOLD = 10.0


class TherapyRanker:
    """Evidence-based therapy ranking engine for precision oncology.

    Integrates variant-level, biomarker-level, and resistance-aware analysis
    to produce a prioritized list of candidate therapies with citations.

    Ranking strategy:
      1. Identify variant-driven therapies from ACTIONABLE_TARGETS
      2. Identify biomarker-driven therapies (MSI-H, TMB-H, HRD, etc.)
      3. Rank by evidence level (A > B > C > D > E)
      4. Flag resistance based on prior therapy history
      5. Flag contraindications (same drug class as prior failed therapy)
      6. Attach supporting evidence from onco_therapies and onco_literature
    """

    def __init__(self, collection_manager, embedder, knowledge):
        """Initialize the therapy ranker.

        Args:
            collection_manager: Milvus collection manager for vector storage.
            embedder: Embedding model wrapper (e.g., BGE-small-en-v1.5).
            knowledge: Knowledge module providing ACTIONABLE_TARGETS,
                THERAPY_MAP, RESISTANCE_MAP, BIOMARKER_PANELS.
        """
        self.collection_manager = collection_manager
        self.embedder = embedder
        self.knowledge = knowledge

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def rank_therapies(
        self,
        cancer_type: str,
        variants: List[Dict],
        biomarkers: Dict[str, Any],
        prior_therapies: Optional[List[str]] = None,
    ) -> List[Dict]:
        """Produce a ranked list of candidate therapies for the patient.

        Args:
            cancer_type: Cancer type (e.g. "NSCLC", "CRC", "Melanoma").
            variants: List of variant dicts, each with at least 'gene' and
                'variant' keys.
            biomarkers: Dict of biomarker results, e.g.
                {"MSI": "MSI-H", "TMB": 14.2, "PD-L1_TPS": 80, "HRD": True}.
            prior_therapies: Optional list of previously administered drugs.

        Returns:
            Ranked list of therapy dicts, each containing:
              rank, drug_name, brand_name, category, targets, evidence_level,
              supporting_evidence, resistance_flag, contraindication_flag,
              guideline_recommendation.
        """
        prior_therapies = prior_therapies or []
        prior_lower = [t.lower().strip() for t in prior_therapies]

        # Step 1: Identify variant-driven therapies
        variant_therapies = []
        for v in variants:
            gene = v.get("gene", "")
            variant_str = v.get("variant", v.get("hgvs", ""))
            vt = self._identify_variant_therapies(gene, variant_str, cancer_type)
            variant_therapies.extend(vt)

        logger.info("Identified %d variant-driven therapy candidates", len(variant_therapies))

        # Step 2: Identify biomarker-driven therapies
        biomarker_therapies = self._identify_biomarker_therapies(biomarkers, cancer_type)
        logger.info("Identified %d biomarker-driven therapy candidates", len(biomarker_therapies))

        # Combine and deduplicate by drug_name
        all_therapies = {}
        for t in variant_therapies + biomarker_therapies:
            drug = t.get("drug_name", "").lower()
            if drug and drug not in all_therapies:
                all_therapies[drug] = t
            elif drug and drug in all_therapies:
                # Keep the one with stronger evidence
                existing = all_therapies[drug]
                existing_level = evidence_rank(existing.get("evidence_level"))
                new_level = evidence_rank(t.get("evidence_level"))
                if new_level < existing_level:
                    all_therapies[drug] = t

        therapies = list(all_therapies.values())

        # Step 3: Rank by evidence level
        therapies.sort(key=lambda t: evidence_rank(t.get("evidence_level")))

        # Step 3.5: Identify combination regimens BEFORE the safety checks below, so that
        # resistance, contraindication and evidence lookup apply to them too. This used to sit
        # after step 6 as "step 6.5", which is why combos carried hard-coded False flags.
        for combo in self._identify_combo_therapies(therapies, cancer_type):
            if combo["drug_name"].lower() not in {t["drug_name"].lower() for t in therapies}:
                therapies.append(combo)

        # Step 4: Check for resistance
        for t in therapies:
            drug_name = t.get("drug_name", "")
            resistance_info = self._check_resistance(drug_name, prior_therapies)
            if resistance_info:
                t["resistance_flag"] = True
                t["resistance_detail"] = resistance_info
            else:
                t["resistance_flag"] = False
                t["resistance_detail"] = None

        # Step 5: Check for contraindications (same class as prior failed therapy)
        for t in therapies:
            t["contraindication_flag"] = self._check_contraindication(
                t, prior_lower
            )
            # _check_contraindication and _check_resistance both match on drug_name, which
            # for a combination is "dabrafenib + trametinib" and matches no single prior
            # therapy. Look inside the regimen so a patient who already failed a component
            # is not re-offered the combination unflagged.
            used_components = [
                c for c in t.get("components", []) if c.lower().strip() in prior_lower
            ]
            if used_components:
                t["contraindication_flag"] = True
                detail = (
                    "already received: " + ", ".join(sorted(used_components))
                )
                t["resistance_detail"] = (
                    f"{t.get('resistance_detail')}; {detail}"
                    if t.get("resistance_detail") else detail
                )

        # Step 6: Retrieve supporting evidence
        for t in therapies:
            drug_name = t.get("drug_name", "")
            evidence = self._fetch_supporting_evidence(drug_name, cancer_type)
            t["supporting_evidence"] = evidence

        # Assign final ranks (resistance/contraindicated therapies demoted)
        therapies = self._assign_final_ranks(therapies)

        logger.info("Returning %d ranked therapies for %s", len(therapies), cancer_type)
        return therapies

    def rank_for_case(self, case: CaseSnapshot) -> List[Dict]:
        """Convenience method to rank therapies for an existing CaseSnapshot.

        Args:
            case: A CaseSnapshot object.

        Returns:
            Ranked list of therapy dicts.
        """
        return self.rank_therapies(
            cancer_type=case.cancer_type,
            variants=case.variants,
            biomarkers=case.biomarkers or {},
            prior_therapies=case.prior_therapies or [],
        )

    # ------------------------------------------------------------------
    # Variant-driven therapy identification
    # ------------------------------------------------------------------

    # Synonyms seen in case data for the cancer_types spellings in knowledge.py. Matching is
    # substring-based on top of this, so "metastatic colorectal carcinoma" still matches
    # "colorectal".
    _CANCER_SYNONYMS = {
        "crc": "colorectal",
        "colorectal carcinoma": "colorectal",
        "nsclc": "nsclc",
        "non-small cell lung cancer": "nsclc",
        "lung adenocarcinoma": "nsclc",
        "rcc": "renal cell",
        "hcc": "hepatocellular",
        "gbm": "glioblastoma",
        "breast carcinoma": "breast",
    }

    @staticmethod
    def _downgrade(level: str, steps: int = 1) -> str:
        """Move an evidence level down `steps` tiers, bottoming out at the weakest."""
        ladder = ["A", "B", "C", "D", "E"]
        try:
            return ladder[min(len(ladder) - 1, ladder.index(level.upper()) + steps)]
        except ValueError:
            return "E"

    def _cancer_type_applies(self, cancer_type: str, allowed: List[str]) -> bool:
        """Is `cancer_type` among the indications `allowed`?

        Nothing gated on cancer type before this: _identify_variant_therapies and
        _identify_combo_therapies both took a cancer_type argument and never read it, so a
        BRAF V600E melanoma case was offered "encorafenib + cetuximab" (BEACON, colorectal
        only) at evidence level A, and a PD-L1 TPS>=50 case got the NSCLC first-line
        pembrolizumab line whatever the tumour was.

        An empty `allowed` means the source data does not state an indication, which is
        treated as unknown rather than approved.
        """
        if not allowed:
            return False
        if any("tissue-agnostic" in a.lower() for a in allowed):
            return True
        if not cancer_type:
            return False
        ct = cancer_type.lower().strip()
        ct = self._CANCER_SYNONYMS.get(ct, ct)
        for a in allowed:
            al = a.lower().strip()
            if al in ct or ct in al:
                return True
        return False

    def _drug_indication_status(
        self, drug_name: str, driver: str, cancer_type: str
    ) -> Dict[str, Any]:
        """Is THIS DRUG approved for THIS DRIVER in THIS TUMOUR, per its own label?

        The per-gene `cancer_types` list answers a different question ("is this gene
        relevant in this tumour"), and because it is the wider list, using it as the
        gate meant the gate never fired: EGFR TKIs came back at level A with
        indicated_for_cancer_type=True for glioblastoma, because EGFR's gene-level list
        names glioblastoma.

        Returns a dict with:
          indicated: True  - the drug's label covers this driver in this tumour
                     False - the label was retrieved and does NOT cover it
                     None  - no label has been retrieved for this (drug, driver);
                             the answer is unknown and must not be asserted
          approved_cancer_types: what the label says, when known
          label: the citation
          note: reviewer-facing explanation when indicated is not True
        """
        entry = VERIFIED_DRUG_INDICATIONS.get((drug_name.lower().strip(), driver.upper().strip()))
        if entry is None:
            entry = VERIFIED_DRUG_INDICATIONS.get((drug_name.lower().strip(), driver.strip()))
        if entry is None:
            return {
                "indicated": None,
                "approved_cancer_types": [],
                "label": "",
                "note": (
                    f"The approved indications for {drug_name} with a {driver} driver have NOT "
                    f"been verified against a label. The tumour list carried here is the "
                    f"gene-level list from ACTIONABLE_TARGETS, which states where the GENE is "
                    f"relevant, not where the DRUG is approved. Confirm against the FDA label "
                    f"(Drugs@FDA / DailyMed) or the NCCN guideline for this tumour before acting."
                ),
            }

        tumors = list(entry.get("tumors", []))
        excluded = list(entry.get("excluded", []))
        label = entry.get("label", "")

        if excluded and self._cancer_type_applies(cancer_type, excluded):
            return {
                "indicated": False,
                "approved_cancer_types": tumors,
                "label": label,
                "note": (
                    f"{drug_name} is explicitly NOT indicated in {cancer_type} for a {driver} "
                    f"driver: the label carries a Limitation of Use excluding it. Label: {label}"
                ),
            }

        if not tumors:
            return {
                "indicated": False,
                "approved_cancer_types": [],
                "label": label,
                "note": (
                    f"{drug_name} has NO approved indication for a {driver} driver at any tumour "
                    f"site. Label: {label}"
                ),
            }

        if self._cancer_type_applies(cancer_type, tumors):
            return {
                "indicated": True,
                "approved_cancer_types": tumors,
                "label": label,
                "note": "",
            }

        return {
            "indicated": False,
            "approved_cancer_types": tumors,
            "label": label,
            "note": (
                f"{drug_name} is not approved in {cancer_type or 'an unstated tumour type'} for a "
                f"{driver} driver; its labelled indications are {tumors}. Label: {label}"
            ),
        }

    def _identify_variant_therapies(
        self, gene: str, variant: str, cancer_type: str
    ) -> List[Dict]:
        """Check ACTIONABLE_TARGETS for therapies targeting a specific variant.

        Args:
            gene: Gene symbol (e.g. "EGFR").
            variant: Variant string (e.g. "L858R").
            cancer_type: Cancer type for context-specific recommendations.

        Returns:
            List of therapy dicts with drug_name, brand_name, category,
            targets, evidence_level, guideline_recommendation.
        """
        gene_upper = gene.upper().strip()
        if gene_upper not in ACTIONABLE_TARGETS:
            return []

        target_info = ACTIONABLE_TARGETS[gene_upper]
        therapies = []

        # Check if variant matches known actionable variants
        key_variants = target_info.get("key_variants", [])
        variant_matched = any(
            kv.upper() in variant.upper()
            for kv in key_variants
        ) if key_variants and variant else False

        # Get drug list from targeted_therapies or drugs
        drugs = target_info.get("targeted_therapies", target_info.get("drugs", []))
        gene_level_inference = False
        if variant_matched:
            evidence_level = target_info.get("evidence_level", "C")
        elif target_info.get("gene_level_actionable", False):
            evidence_level = target_info.get("default_evidence_level", "C")
            gene_level_inference = True
        elif key_variants:
            # The gene has known actionable variants and THIS variant is not one of them.
            # This branch used to hand back the gene's top evidence_level unchanged, so an
            # EGFR C797S or an exon-20 insertion was offered osimertinib and erlotinib at
            # level A — the level earned by L858R, not by this variant. Evidence for a
            # different variant is not evidence for this one: infer, but say so and drop a
            # tier.
            evidence_level = self._downgrade(target_info.get("evidence_level", "C"))
            gene_level_inference = True
        else:
            return []

        # A variant listed as a resistance mechanism for this gene is surfaced with a flag
        # rather than silently ranked. Which drug it defeats is drug-specific (T790M is
        # osimertinib's indication and 1st/2nd-gen TKI resistance), so this does not presume
        # to pick: it marks the variant so a reviewer sees it.
        resistance_variants = [r.upper() for r in target_info.get("resistance_mutations", [])]
        variant_is_resistance = bool(variant) and any(
            r in variant.upper() for r in resistance_variants
        )
        review_note = ""
        if variant_is_resistance:
            review_note = (
                f"{variant} is listed as a resistance mechanism for {gene_upper}. Whether it "
                f"defeats a given drug in this list is drug-specific and is NOT encoded here "
                f"(for example T790M is osimertinib's indication while conferring resistance "
                f"to earlier-generation TKIs). Confirm per drug before acting."
            )

        # The gene-level list says where the GENE is relevant. It is kept as a
        # fallback for drugs whose own label has not been read, but it is no longer
        # what decides `indicated_for_cancer_type` on its own.
        gene_list_applies = self._cancer_type_applies(
            cancer_type, target_info.get("cancer_types", [])
        )

        for drug in drugs:
            if isinstance(drug, str):
                drug_name = drug
                drug_details = THERAPY_MAP.get(drug_name.lower(), {})
            elif isinstance(drug, dict):
                drug_name = drug.get("name", "")
                drug_details = drug
            else:
                continue

            # Per-drug gate: the drug's own label decides whether it is indicated
            # for this driver in this tumour.
            status = self._drug_indication_status(drug_name, gene_upper, cancer_type)
            drug_level = evidence_level
            drug_note = review_note
            if status["indicated"] is True:
                indicated = True
                approved_types = status["approved_cancer_types"]
            elif status["indicated"] is False:
                # The label was retrieved and does not cover this tumour. That is a
                # stronger statement than "the gene list omits it", so drop two
                # tiers rather than one, and say why.
                indicated = False
                approved_types = status["approved_cancer_types"]
                drug_level = self._downgrade(evidence_level, steps=2)
                drug_note = f"{drug_note} {status['note']}".strip()
            else:
                # No label retrieved for this (drug, driver): unknown, not approved.
                # Behave as before — gate on the gene-level list — but mark the
                # recommendation as unverified so no reader mistakes it for a
                # label-confirmed indication.
                indicated = gene_list_applies
                approved_types = list(target_info.get("cancer_types", []))
                if not gene_list_applies:
                    drug_level = self._downgrade(evidence_level)
                drug_note = f"{drug_note} {status['note']}".strip()

            therapies.append({
                "drug_name": drug_name,
                "brand_name": drug_details.get("brand_name", ""),
                "category": drug_details.get("category", "targeted therapy"),
                "targets": [gene_upper],
                "evidence_level": drug_level,
                "guideline_recommendation": drug_details.get(
                    "guideline", f"Consider for {gene_upper}-mutant {cancer_type}"
                ),
                "source": "variant",
                "source_gene": gene_upper,
                "source_variant": variant,
                "variant_matched": variant_matched,
                # True when the recommendation rests on the gene rather than on this variant.
                "gene_level_inference": gene_level_inference,
                "variant_is_resistance_mutation": variant_is_resistance,
                "review_note": drug_note,
                "indicated_for_cancer_type": indicated,
                "approved_cancer_types": list(approved_types),
                # True when no FDA label has been read for this (drug, driver), so
                # `approved_cancer_types` is the gene-level list, not an indication.
                "indication_unverified": status["indicated"] is None,
                "indication_label": status["label"],
            })

        return therapies

    # ------------------------------------------------------------------
    # Biomarker-driven therapy identification
    # ------------------------------------------------------------------

    def _identify_biomarker_therapies(
        self, biomarkers: Dict[str, Any], cancer_type: str
    ) -> List[Dict]:
        """Check BIOMARKER_PANELS for therapies driven by biomarker status.

        Standard biomarker-therapy mappings:
          - MSI-H / dMMR -> pembrolizumab (Keytruda)
          - TMB-H (>=10 mut/Mb) -> pembrolizumab
          - HRD / BRCA -> PARP inhibitors (olaparib, rucaparib, niraparib)
          - PD-L1 TPS >=50% -> pembrolizumab first-line
          - NTRK fusion -> larotrectinib, entrectinib

        Args:
            biomarkers: Patient biomarker dict.
            cancer_type: Cancer type for context.

        Returns:
            List of therapy dicts.
        """
        therapies = []

        def _entry(
            drug_name: str,
            brand: str,
            category: str,
            targets: List[str],
            level: str,
            guideline: str,
            driver: str,
            source_biomarker: str,
            source: str = "biomarker",
        ) -> Dict[str, Any]:
            """Build a biomarker therapy dict with the per-drug indication gate applied.

            This function used to take `cancer_type` and never read it: the MSI, TMB,
            HRD/BRCA, PTEN, PD-L1 and NTRK blocks all emitted unconditionally, so a
            PD-L1 TPS 80 PROSTATE case got "First-line monotherapy for PD-L1 TPS >= 50%
            NSCLC (KEYNOTE-024)" at level A as its rank-1 recommendation, and an MSI-H
            ENDOMETRIAL case got nivolumab at level A quoting a colorectal approval.
            The biomarker dicts also omitted `indicated_for_cancer_type`, so no
            downstream consumer could warn.
            """
            status = self._drug_indication_status(drug_name, driver, cancer_type)
            out_level = level
            note = ""
            if status["indicated"] is False:
                out_level = self._downgrade(level, steps=2)
                note = status["note"]
            elif status["indicated"] is None:
                note = status["note"]
            return {
                "drug_name": drug_name,
                "brand_name": brand,
                "category": category,
                "targets": list(targets),
                "evidence_level": out_level,
                "guideline_recommendation": guideline,
                "source": source,
                "source_biomarker": source_biomarker,
                "indicated_for_cancer_type": status["indicated"],
                "approved_cancer_types": status["approved_cancer_types"],
                "indication_unverified": status["indicated"] is None,
                "indication_label": status["label"],
                "review_note": note,
            }

        # MSI-H -> pembrolizumab, nivolumab, dostarlimab
        msi = str(biomarkers.get("MSI") or "").upper()
        if msi in ("MSI-H", "MSI-HIGH", "DMMR"):
            msi_drugs = [
                ("pembrolizumab", "Keytruda", "A",
                 "FDA-approved for MSI-H/dMMR solid tumors regardless of tumor type (KEYNOTE-158, KEYNOTE-177)."),
                ("nivolumab", "Opdivo", "A",
                 "FDA-approved for MSI-H/dMMR metastatic colorectal cancer (CheckMate-142)."),
                ("dostarlimab", "Jemperli", "A",
                 "FDA-approved for dMMR solid tumors (GARNET trial)."),
            ]
            for drug_name, brand, level, guideline in msi_drugs:
                therapies.append(_entry(
                    drug_name, brand, "immunotherapy", ["MSI-H"], level, guideline,
                    driver="MSI-H", source_biomarker="MSI-H",
                ))

        # TMB-H -> pembrolizumab, atezolizumab
        #
        # The labelled definition of TMB-H is >= 13 mut/Mb (KEYTRUDA, DailyMed setid
        # 9333c79b-d487-4538-a9f0-71b91a02b287, published 2026-11-02, §1.16), not the
        # >= 10 this block used while claiming "FDA-approved for TMB-H (>=10 mut/Mb)".
        # KEYNOTE-158's cohort used >= 10, which is where that figure came from, so a
        # patient between the two is still reported — as a trial-supported candidate,
        # not as an approved indication.
        tmb = biomarkers.get("TMB")
        if tmb is not None and isinstance(tmb, (int, float)) and tmb >= TMB_H_TRIAL_THRESHOLD:
            meets_label = tmb >= TMB_H_LABEL_THRESHOLD
            if meets_label:
                pembro_level = "A"
                pembro_guideline = (
                    f"FDA-approved for TMB-H (>={TMB_H_LABEL_THRESHOLD:.0f} mut/Mb) solid tumors. "
                    f"Patient TMB: {tmb} mut/Mb (KEYNOTE-158)."
                )
            else:
                pembro_level = "B"
                pembro_guideline = (
                    f"Patient TMB {tmb} mut/Mb is BELOW the labelled TMB-H definition of "
                    f">={TMB_H_LABEL_THRESHOLD:.0f} mut/Mb (KEYTRUDA label, DailyMed, published "
                    f"2026-11-02, section 1.16), so this is NOT within the approved indication. "
                    f"KEYNOTE-158's TMB-H cohort used >={TMB_H_TRIAL_THRESHOLD:.0f} mut/Mb."
                )
            tmb_drugs = [
                ("pembrolizumab", "Keytruda", pembro_level, pembro_guideline),
                ("atezolizumab", "Tecentriq", "B",
                 f"Evidence supports activity in TMB-H tumors. Patient TMB: {tmb} mut/Mb (IMvigor211, POPLAR)."),
            ]
            for drug_name, brand, level, guideline in tmb_drugs:
                therapies.append(_entry(
                    drug_name, brand, "immunotherapy", ["TMB-H"], level, guideline,
                    driver="TMB-H", source_biomarker="TMB-H",
                ))

        # HRD / BRCA -> PARP inhibitors + platinum-based chemotherapy
        hrd = biomarkers.get("HRD")
        brca = str(biomarkers.get("BRCA") or "").upper()
        if hrd in (True, "positive", "Positive") or brca in ("BRCA1", "BRCA2", "POSITIVE"):
            parp_drugs = [
                ("olaparib", "Lynparza"),
                ("rucaparib", "Rubraca"),
                ("niraparib", "Zejula"),
                ("talazoparib", "Talzenna"),
            ]
            for drug_name, brand in parp_drugs:
                therapies.append(_entry(
                    drug_name, brand, "PARP inhibitor", ["HRD", "BRCA"],
                    "A" if brca else "B",
                    f"Consider {brand} ({drug_name}) for HRD-positive / BRCA-mutant tumors.",
                    driver="BRCA", source_biomarker="HRD/BRCA",
                ))
            # Platinum compounds are also effective in HRD/BRCA. This is a biology
            # claim from the literature, not a labelled biomarker indication: no
            # platinum label selects patients on BRCA or HRD status. It is left at
            # level B and marked unverified rather than being re-levelled here —
            # recorded as an open question in CLINICAL_REVIEW.md.
            for plat_name in ("carboplatin", "cisplatin"):
                therapies.append(_entry(
                    plat_name, "", "platinum chemotherapy", ["HRD", "BRCA"], "B",
                    "BRCA/HRD-positive tumors show enhanced platinum sensitivity.",
                    driver="BRCA", source_biomarker="HRD/BRCA",
                ))

        # PTEN loss -> consider mTOR/PI3K inhibitors
        pten = str(biomarkers.get("PTEN") or "").upper()
        if pten in ("LOSS", "DELETED", "MUTATION", "LOST"):
            therapies.append(_entry(
                "alpelisib", "Piqray", "PI3K inhibitor", ["PTEN", "PI3K"], "C",
                "Consider PI3K pathway inhibition for PTEN-loss tumors. "
                "PI3K/AKT/mTOR pathway activation may predict response.",
                driver="PTEN", source_biomarker="PTEN_loss",
            ))

        # PD-L1 TPS >= 50% -> pembrolizumab first-line
        pdl1 = biomarkers.get("PD-L1_TPS")
        if pdl1 is not None and isinstance(pdl1, (int, float)) and pdl1 >= 50:
            therapies.append(_entry(
                "pembrolizumab", "Keytruda", "immunotherapy", ["PD-L1"], "A",
                f"First-line monotherapy for PD-L1 TPS >= 50% NSCLC. "
                f"Patient PD-L1 TPS: {pdl1}% (KEYNOTE-024).",
                driver="PD-L1_TPS", source_biomarker="PD-L1_TPS",
            ))

        # NTRK fusion -> larotrectinib, entrectinib
        ntrk = str(biomarkers.get("NTRK") or "").upper()
        if ntrk in ("FUSION", "POSITIVE", "NTRK1", "NTRK2", "NTRK3"):
            therapies.append(_entry(
                "larotrectinib", "Vitrakvi", "targeted therapy", ["NTRK"], "A",
                "FDA-approved for NTRK fusion-positive solid tumors (NAVIGATE, SCOUT).",
                driver="NTRK", source_biomarker="NTRK",
            ))
            therapies.append(_entry(
                "entrectinib", "Rozlytrek", "targeted therapy", ["NTRK", "ROS1"], "A",
                "FDA-approved for NTRK fusion-positive solid tumors (STARTRK-2).",
                driver="NTRK", source_biomarker="NTRK",
            ))

        # Also check BIOMARKER_PANELS for any additional mappings
        for panel_key, panel_info in BIOMARKER_PANELS.items():
            marker_name = panel_info.get("marker", "")
            threshold = panel_info.get("threshold")
            patient_value = biomarkers.get(marker_name)

            if patient_value is None:
                continue

            # Check if patient meets threshold
            meets = False
            if threshold is None:
                meets = bool(patient_value)
            elif isinstance(patient_value, (int, float)) and isinstance(threshold, (int, float)):
                meets = patient_value >= threshold
            elif isinstance(patient_value, str):
                positive_values = panel_info.get("positive_values", [])
                meets = patient_value.upper() in [pv.upper() for pv in positive_values]

            if meets:
                # NOTE: as of 2026-11-03 no entry in BIOMARKER_PANELS carries a
                # `marker`, `threshold` or `drugs` key, so this loop never emits
                # anything (marker_name is "" and `biomarkers.get("")` is None, so the
                # `continue` above always fires). It is gated anyway, so that the day
                # those keys appear it does not reproduce the ungated-biomarker defect
                # this commit fixes. Recorded in CLINICAL_REVIEW.md.
                panel_drugs = panel_info.get("drugs", [])
                panel_cancer_types = panel_info.get("cancer_types", [])
                panel_applies = self._cancer_type_applies(cancer_type, panel_cancer_types)
                panel_level = panel_info.get("evidence_level", "C")
                if not panel_applies:
                    panel_level = self._downgrade(panel_level)
                for pd in panel_drugs:
                    drug_name = pd if isinstance(pd, str) else pd.get("name", "")
                    if drug_name and drug_name.lower() not in {t["drug_name"].lower() for t in therapies}:
                        therapies.append({
                            "drug_name": drug_name,
                            "brand_name": pd.get("brand_name", "") if isinstance(pd, dict) else "",
                            "category": "biomarker-directed",
                            "targets": [marker_name],
                            "evidence_level": panel_level,
                            "guideline_recommendation": panel_info.get("guideline", ""),
                            "source": "biomarker_panel",
                            "source_biomarker": panel_key,
                            "indicated_for_cancer_type": panel_applies,
                            "approved_cancer_types": list(panel_cancer_types),
                            "indication_unverified": True,
                            "review_note": (
                                f"{drug_name} comes from BIOMARKER_PANELS[{panel_key!r}], whose "
                                f"cancer_types list has not been verified against a label."
                            ),
                        })

        return therapies

    # ------------------------------------------------------------------
    # Resistance and contraindication checks
    # ------------------------------------------------------------------

    # Drug class groupings for same-mechanism resistance detection
    _DRUG_CLASS_GROUPS = {
        "egfr_tki_1g": ["erlotinib", "gefitinib"],
        "egfr_tki_2g": ["afatinib", "dacomitinib"],
        "egfr_tki_3g": ["osimertinib"],
        "alk_tki": ["crizotinib", "ceritinib", "alectinib", "brigatinib", "lorlatinib"],
        "braf_inhibitor": ["vemurafenib", "dabrafenib", "encorafenib"],
        "mek_inhibitor": ["trametinib", "cobimetinib", "binimetinib"],
        "anti_pd1": ["pembrolizumab", "nivolumab", "dostarlimab", "cemiplimab"],
        "anti_pdl1": ["atezolizumab", "durvalumab", "avelumab"],
        "parp_inhibitor": ["olaparib", "rucaparib", "niraparib", "talazoparib"],
        "kras_g12c": ["sotorasib", "adagrasib"],
    }

    def _check_resistance(
        self, drug: str, prior_therapies: List[str]
    ) -> Optional[Dict]:
        """Check if a candidate drug has known resistance from prior therapy.

        Uses RESISTANCE_MAP for mutation-level resistance, and _DRUG_CLASS_GROUPS
        for same-mechanism class-level resistance (e.g., prior erlotinib confers
        likely resistance to gefitinib since both are 1st-gen EGFR TKIs).

        Args:
            drug: Candidate drug name.
            prior_therapies: List of prior therapy names.

        Returns:
            Dict with resistance details if found, otherwise None.
        """
        drug_lower = drug.lower().strip()
        prior_lower = [p.lower().strip() for p in prior_therapies]

        # 1. Check RESISTANCE_MAP for mutation-level resistance
        resistance_entry = RESISTANCE_MAP.get(drug_lower)
        if resistance_entry:
            if isinstance(resistance_entry, list):
                trigger_mutations = [
                    m.get("mutation", "").lower()
                    for m in resistance_entry
                    if isinstance(m, dict) and m.get("mutation")
                ]
                overlapping = [
                    p for p in prior_lower
                    if any(trig in p or p in trig for trig in trigger_mutations)
                ]
                if overlapping:
                    return {
                        "drug": drug,
                        "prior_triggers": overlapping,
                        "mechanism": "; ".join(
                            m.get("mutation", "") for m in resistance_entry
                            if isinstance(m, dict)
                        ),
                        "alternatives": [
                            alt for m in resistance_entry
                            if isinstance(m, dict)
                            for alt in m.get("next_line", [])
                        ],
                    }
            elif isinstance(resistance_entry, dict):
                triggers = resistance_entry.get("resistance_triggers", [])
                trigger_lower = [t.lower() for t in triggers]
                overlapping = [
                    p for p in prior_lower
                    if any(trig in p or p in trig for trig in trigger_lower)
                ]
                if overlapping:
                    return {
                        "drug": drug,
                        "prior_triggers": overlapping,
                        "mechanism": resistance_entry.get("mechanism", "Unknown resistance mechanism"),
                        "alternatives": resistance_entry.get("alternatives", []),
                    }

        # 2. Check same-mechanism class-level resistance
        for class_name, class_drugs in self._DRUG_CLASS_GROUPS.items():
            if drug_lower in class_drugs:
                prior_same_class = [
                    p for p in prior_lower
                    if p in class_drugs and p != drug_lower
                ]
                if prior_same_class:
                    return {
                        "drug": drug,
                        "prior_triggers": prior_same_class,
                        "mechanism": (
                            f"Same drug class ({class_name.replace('_', ' ')}): "
                            f"prior {', '.join(prior_same_class)} suggests likely "
                            f"cross-resistance to {drug}"
                        ),
                        "alternatives": [],
                    }

        return None

    def _check_contraindication(self, therapy: Dict, prior_lower: List[str]) -> bool:
        """Check if therapy is in the same drug class as a prior failed therapy.

        A contraindication flag is raised when the patient has previously
        received a drug in the same therapeutic category/class.

        Args:
            therapy: Therapy dict with category and drug_name.
            prior_lower: Lowercased list of prior therapy names.

        Returns:
            True if contraindicated, False otherwise.
        """
        category = therapy.get("category", "").lower()
        drug_name = therapy.get("drug_name", "").lower()

        # Direct match: same drug was used before
        if drug_name in prior_lower:
            return True

        # Same class check via THERAPY_MAP
        for prior in prior_lower:
            prior_info = THERAPY_MAP.get(prior, {})
            prior_category = prior_info.get("category", "").lower()
            if prior_category and prior_category == category:
                # Same broad category might be a concern, but only flag if
                # the specific drug class matches
                prior_class = prior_info.get("drug_class", "").lower()
                therapy_class = THERAPY_MAP.get(drug_name, {}).get("drug_class", "").lower()
                if prior_class and therapy_class and prior_class == therapy_class:
                    return True

        return False

    # ------------------------------------------------------------------
    # Evidence retrieval
    # ------------------------------------------------------------------

    def _fetch_supporting_evidence(
        self, drug_name: str, cancer_type: str
    ) -> List[Dict]:
        """Search onco_therapies and onco_literature for supporting evidence.

        Args:
            drug_name: Drug name to search for.
            cancer_type: Cancer type for context.

        Returns:
            List of evidence dicts with source, text, and score.
        """
        query_text = f"{drug_name} {cancer_type} clinical evidence efficacy"
        evidence = []

        for collection in ("onco_therapies", "onco_literature"):
            try:
                embedding = self.embedder.embed(query_text)
                results = self.collection_manager.search(
                    collection_name=collection,
                    query_vector=embedding,
                    top_k=3,
                    output_fields=["text", "source", "title"],
                )
                for r in results:
                    evidence.append({
                        "collection": collection,
                        "source": r.get("source", r.get("title", "")),
                        "text": r.get("text", "")[:500],
                        "score": round(r.get("score", 0.0), 4),
                    })
            except Exception as exc:
                logger.warning(
                    "Evidence retrieval from %s failed for %s: %s",
                    collection, drug_name, exc,
                )

        return evidence

    # ------------------------------------------------------------------
    # Combination therapy identification
    # ------------------------------------------------------------------

    # Known FDA-approved combination regimens
    _COMBO_REGIMENS = [
        {
            "combo_name": "dabrafenib + trametinib",
            "cancer_types": ["melanoma", "NSCLC", "thyroid"],
            "components": ["dabrafenib", "trametinib"],
            "target": "BRAF",
            "evidence_level": "A",
            "guideline": "FDA-approved BRAF+MEK combination for BRAF V600E/K melanoma, NSCLC, thyroid (COMBI-d, COMBI-v).",
        },
        {
            "combo_name": "encorafenib + binimetinib",
            "cancer_types": ["melanoma"],
            "components": ["encorafenib", "binimetinib"],
            "target": "BRAF",
            "evidence_level": "A",
            "guideline": "FDA-approved BRAF+MEK combination for BRAF V600E melanoma (COLUMBUS).",
        },
        {
            "combo_name": "encorafenib + cetuximab",
            "cancer_types": ["colorectal"],
            "components": ["encorafenib", "cetuximab"],
            "target": "BRAF",
            "evidence_level": "A",
            "guideline": "FDA-approved for BRAF V600E metastatic CRC (BEACON).",
        },
        {
            "combo_name": "ipilimumab + nivolumab",
            "cancer_types": ["melanoma", "renal cell", "colorectal", "hepatocellular", "NSCLC"],
            "components": ["ipilimumab", "nivolumab"],
            "target": "PD-1/CTLA-4",
            "evidence_level": "A",
            "guideline": "FDA-approved dual checkpoint for melanoma, RCC, MSI-H CRC, HCC, NSCLC (CheckMate-067, -214, -142).",
        },
        {
            "combo_name": "lenvatinib + pembrolizumab",
            "cancer_types": ["endometrial", "renal cell"],
            "components": ["lenvatinib", "pembrolizumab"],
            "target": "VEGFR/PD-1",
            "evidence_level": "A",
            "guideline": "FDA-approved for endometrial carcinoma and RCC (KEYNOTE-775, CLEAR).",
        },
        {
            "combo_name": "trastuzumab + pertuzumab",
            "cancer_types": ["breast"],
            "components": ["trastuzumab", "pertuzumab"],
            "target": "HER2",
            "evidence_level": "A",
            "guideline": "FDA-approved dual HER2 blockade for HER2+ breast cancer (CLEOPATRA).",
        },
    ]

    def _identify_combo_therapies(
        self, existing_therapies: List[Dict], cancer_type: str
    ) -> List[Dict]:
        """Identify known combination regimens when component drugs are present.

        If a single-agent therapy is identified, check if it belongs to an
        FDA-approved combination regimen and suggest the combo.

        Args:
            existing_therapies: Already-identified single-agent therapies.
            cancer_type: Cancer type for context.

        Returns:
            List of combo therapy dicts to add.
        """
        existing_drugs = {t["drug_name"].lower() for t in existing_therapies}
        combos = []

        for regimen in self._COMBO_REGIMENS:
            components = regimen["components"]
            # If at least one component is already in the therapy list
            overlap = [c for c in components if c.lower() in existing_drugs]
            if not overlap:
                continue
            # cancer_types comes from each regimen's own guideline string; a regimen is not
            # offered for a tumour its approval does not cover.
            if not self._cancer_type_applies(cancer_type, regimen.get("cancer_types", [])):
                continue
            combos.append({
                "drug_name": regimen["combo_name"],
                "components": list(components),
                "brand_name": "",
                "category": "combination regimen",
                "targets": [regimen["target"]],
                "evidence_level": regimen["evidence_level"],
                "guideline_recommendation": regimen["guideline"],
                "source": "combination",
                "approved_cancer_types": list(regimen.get("cancer_types", [])),
                "indicated_for_cancer_type": True,
                # resistance_flag / contraindication_flag / supporting_evidence are no
                # longer hard-coded here: combos used to be appended AFTER steps 4-6, so
                # they shipped resistance_flag=False and an empty evidence list by
                # construction, and a patient who had already failed dabrafenib+trametinib
                # was re-offered it unflagged. rank_therapies() now identifies combos
                # before those checks so the same loops cover them.
            })

        return combos

    # ------------------------------------------------------------------
    # Final ranking
    # ------------------------------------------------------------------

    def _assign_final_ranks(self, therapies: List[Dict]) -> List[Dict]:
        """Assign final ranks, demoting resistance-flagged and contraindicated therapies.

        Therapies are first sorted by evidence level, then those with
        resistance or contraindication flags are moved to the bottom of their
        evidence tier.

        Args:
            therapies: List of therapy dicts (already sorted by evidence level).

        Returns:
            Re-ranked list with 'rank' field assigned.
        """
        # Partition into clean and flagged
        clean = [t for t in therapies if not t.get("resistance_flag") and not t.get("contraindication_flag")]
        flagged = [t for t in therapies if t.get("resistance_flag") or t.get("contraindication_flag")]

        # Sort each group by evidence level
        clean.sort(key=lambda t: evidence_rank(t.get("evidence_level")))
        flagged.sort(key=lambda t: evidence_rank(t.get("evidence_level")))

        # Combine: clean first, then flagged
        ranked = clean + flagged
        for i, t in enumerate(ranked):
            t["rank"] = i + 1

        return ranked
