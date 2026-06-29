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
                existing_level = EVIDENCE_LEVEL_ORDER.get(existing.get("evidence_level", "E"), 4)
                new_level = EVIDENCE_LEVEL_ORDER.get(t.get("evidence_level", "E"), 4)
                if new_level < existing_level:
                    all_therapies[drug] = t

        therapies = list(all_therapies.values())

        # Step 3: Rank by evidence level
        therapies.sort(key=lambda t: EVIDENCE_LEVEL_ORDER.get(t.get("evidence_level", "E"), 4))

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

        # Step 6: Retrieve supporting evidence
        for t in therapies:
            drug_name = t.get("drug_name", "")
            evidence = self._fetch_supporting_evidence(drug_name, cancer_type)
            t["supporting_evidence"] = evidence

        # Step 6.5: Identify combination therapies
        combo_therapies = self._identify_combo_therapies(therapies, cancer_type)
        if combo_therapies:
            for combo in combo_therapies:
                drug_key = combo["drug_name"].lower()
                if drug_key not in {t["drug_name"].lower() for t in therapies}:
                    therapies.append(combo)

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
        if variant_matched:
            evidence_level = target_info.get("evidence_level", "C")
        elif target_info.get("gene_level_actionable", False):
            evidence_level = target_info.get("default_evidence_level", "C")
        elif key_variants:
            # Gene has known variants but this specific variant wasn't matched
            evidence_level = target_info.get("evidence_level", "C")
        else:
            return []

        for drug in drugs:
            if isinstance(drug, str):
                drug_name = drug
                drug_details = THERAPY_MAP.get(drug_name.lower(), {})
            elif isinstance(drug, dict):
                drug_name = drug.get("name", "")
                drug_details = drug
            else:
                continue

            therapies.append({
                "drug_name": drug_name,
                "brand_name": drug_details.get("brand_name", ""),
                "category": drug_details.get("category", "targeted therapy"),
                "targets": [gene_upper],
                "evidence_level": evidence_level,
                "guideline_recommendation": drug_details.get(
                    "guideline", f"Consider for {gene_upper}-mutant {cancer_type}"
                ),
                "source": "variant",
                "source_gene": gene_upper,
                "source_variant": variant,
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

        # MSI-H -> pembrolizumab, nivolumab, dostarlimab
        msi = biomarkers.get("MSI", "").upper()
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
                therapies.append({
                    "drug_name": drug_name,
                    "brand_name": brand,
                    "category": "immunotherapy",
                    "targets": ["MSI-H"],
                    "evidence_level": level,
                    "guideline_recommendation": guideline,
                    "source": "biomarker",
                    "source_biomarker": "MSI-H",
                })

        # TMB-H -> pembrolizumab, atezolizumab
        tmb = biomarkers.get("TMB")
        if tmb is not None and isinstance(tmb, (int, float)) and tmb >= 10:
            tmb_drugs = [
                ("pembrolizumab", "Keytruda", "A",
                 f"FDA-approved for TMB-H (>=10 mut/Mb) solid tumors. Patient TMB: {tmb} mut/Mb (KEYNOTE-158)."),
                ("atezolizumab", "Tecentriq", "B",
                 f"Evidence supports activity in TMB-H tumors. Patient TMB: {tmb} mut/Mb (IMvigor211, POPLAR)."),
            ]
            for drug_name, brand, level, guideline in tmb_drugs:
                therapies.append({
                    "drug_name": drug_name,
                    "brand_name": brand,
                    "category": "immunotherapy",
                    "targets": ["TMB-H"],
                    "evidence_level": level,
                    "guideline_recommendation": guideline,
                    "source": "biomarker",
                    "source_biomarker": "TMB-H",
                })

        # HRD / BRCA -> PARP inhibitors + platinum-based chemotherapy
        hrd = biomarkers.get("HRD")
        brca = biomarkers.get("BRCA", "").upper()
        if hrd in (True, "positive", "Positive") or brca in ("BRCA1", "BRCA2", "POSITIVE"):
            parp_drugs = [
                ("olaparib", "Lynparza"),
                ("rucaparib", "Rubraca"),
                ("niraparib", "Zejula"),
                ("talazoparib", "Talzenna"),
            ]
            for drug_name, brand in parp_drugs:
                therapies.append({
                    "drug_name": drug_name,
                    "brand_name": brand,
                    "category": "PARP inhibitor",
                    "targets": ["HRD", "BRCA"],
                    "evidence_level": "A" if brca else "B",
                    "guideline_recommendation": (
                        f"Consider {brand} ({drug_name}) for HRD-positive / BRCA-mutant tumors."
                    ),
                    "source": "biomarker",
                    "source_biomarker": "HRD/BRCA",
                })
            # Platinum compounds are also effective in HRD/BRCA
            for plat_name in ("carboplatin", "cisplatin"):
                therapies.append({
                    "drug_name": plat_name,
                    "brand_name": "",
                    "category": "platinum chemotherapy",
                    "targets": ["HRD", "BRCA"],
                    "evidence_level": "B",
                    "guideline_recommendation": (
                        "BRCA/HRD-positive tumors show enhanced platinum sensitivity."
                    ),
                    "source": "biomarker",
                    "source_biomarker": "HRD/BRCA",
                })

        # PTEN loss -> consider mTOR/PI3K inhibitors
        pten = biomarkers.get("PTEN", "").upper()
        if pten in ("LOSS", "DELETED", "MUTATION", "LOST"):
            therapies.append({
                "drug_name": "alpelisib",
                "brand_name": "Piqray",
                "category": "PI3K inhibitor",
                "targets": ["PTEN", "PI3K"],
                "evidence_level": "C",
                "guideline_recommendation": (
                    "Consider PI3K pathway inhibition for PTEN-loss tumors. "
                    "PI3K/AKT/mTOR pathway activation may predict response."
                ),
                "source": "biomarker",
                "source_biomarker": "PTEN_loss",
            })

        # PD-L1 TPS >= 50% -> pembrolizumab first-line
        pdl1 = biomarkers.get("PD-L1_TPS")
        if pdl1 is not None and isinstance(pdl1, (int, float)) and pdl1 >= 50:
            therapies.append({
                "drug_name": "pembrolizumab",
                "brand_name": "Keytruda",
                "category": "immunotherapy",
                "targets": ["PD-L1"],
                "evidence_level": "A",
                "guideline_recommendation": (
                    f"First-line monotherapy for PD-L1 TPS >= 50% NSCLC. "
                    f"Patient PD-L1 TPS: {pdl1}% (KEYNOTE-024)."
                ),
                "source": "biomarker",
                "source_biomarker": "PD-L1_TPS",
            })

        # NTRK fusion -> larotrectinib, entrectinib
        ntrk = biomarkers.get("NTRK", "").upper()
        if ntrk in ("FUSION", "POSITIVE", "NTRK1", "NTRK2", "NTRK3"):
            therapies.append({
                "drug_name": "larotrectinib",
                "brand_name": "Vitrakvi",
                "category": "targeted therapy",
                "targets": ["NTRK"],
                "evidence_level": "A",
                "guideline_recommendation": (
                    "FDA-approved for NTRK fusion-positive solid tumors (NAVIGATE, SCOUT)."
                ),
                "source": "biomarker",
                "source_biomarker": "NTRK",
            })
            therapies.append({
                "drug_name": "entrectinib",
                "brand_name": "Rozlytrek",
                "category": "targeted therapy",
                "targets": ["NTRK", "ROS1"],
                "evidence_level": "A",
                "guideline_recommendation": (
                    "FDA-approved for NTRK fusion-positive solid tumors (STARTRK-2)."
                ),
                "source": "biomarker",
                "source_biomarker": "NTRK",
            })

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
                panel_drugs = panel_info.get("drugs", [])
                for pd in panel_drugs:
                    drug_name = pd if isinstance(pd, str) else pd.get("name", "")
                    if drug_name and drug_name.lower() not in {t["drug_name"].lower() for t in therapies}:
                        therapies.append({
                            "drug_name": drug_name,
                            "brand_name": pd.get("brand_name", "") if isinstance(pd, dict) else "",
                            "category": "biomarker-directed",
                            "targets": [marker_name],
                            "evidence_level": panel_info.get("evidence_level", "C"),
                            "guideline_recommendation": panel_info.get("guideline", ""),
                            "source": "biomarker_panel",
                            "source_biomarker": panel_key,
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
            "components": ["dabrafenib", "trametinib"],
            "target": "BRAF",
            "evidence_level": "A",
            "guideline": "FDA-approved BRAF+MEK combination for BRAF V600E/K melanoma, NSCLC, thyroid (COMBI-d, COMBI-v).",
        },
        {
            "combo_name": "encorafenib + binimetinib",
            "components": ["encorafenib", "binimetinib"],
            "target": "BRAF",
            "evidence_level": "A",
            "guideline": "FDA-approved BRAF+MEK combination for BRAF V600E melanoma (COLUMBUS).",
        },
        {
            "combo_name": "encorafenib + cetuximab",
            "components": ["encorafenib", "cetuximab"],
            "target": "BRAF",
            "evidence_level": "A",
            "guideline": "FDA-approved for BRAF V600E metastatic CRC (BEACON).",
        },
        {
            "combo_name": "ipilimumab + nivolumab",
            "components": ["ipilimumab", "nivolumab"],
            "target": "PD-1/CTLA-4",
            "evidence_level": "A",
            "guideline": "FDA-approved dual checkpoint for melanoma, RCC, MSI-H CRC, HCC, NSCLC (CheckMate-067, -214, -142).",
        },
        {
            "combo_name": "lenvatinib + pembrolizumab",
            "components": ["lenvatinib", "pembrolizumab"],
            "target": "VEGFR/PD-1",
            "evidence_level": "A",
            "guideline": "FDA-approved for endometrial carcinoma and RCC (KEYNOTE-775, CLEAR).",
        },
        {
            "combo_name": "trastuzumab + pertuzumab",
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
            if overlap:
                combos.append({
                    "drug_name": regimen["combo_name"],
                    "brand_name": "",
                    "category": "combination regimen",
                    "targets": [regimen["target"]],
                    "evidence_level": regimen["evidence_level"],
                    "guideline_recommendation": regimen["guideline"],
                    "source": "combination",
                    "resistance_flag": False,
                    "resistance_detail": None,
                    "contraindication_flag": False,
                    "supporting_evidence": [],
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
        clean.sort(key=lambda t: EVIDENCE_LEVEL_ORDER.get(t.get("evidence_level", "E"), 4))
        flagged.sort(key=lambda t: EVIDENCE_LEVEL_ORDER.get(t.get("evidence_level", "E"), 4))

        # Combine: clean first, then flagged
        ranked = clean + flagged
        for i, t in enumerate(ranked):
            t["rank"] = i + 1

        return ranked
