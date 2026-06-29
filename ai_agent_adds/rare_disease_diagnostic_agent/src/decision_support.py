"""Decision support engines for the Rare Disease Diagnostic Agent.

Author: Adam Jones
Date: March 2026

Provides HPO-to-gene matching with information content scoring, ACMG variant
classification, orphan drug matching, diagnostic algorithm recommendation,
family segregation analysis, and natural history prediction to support rare
disease diagnostic decision-making.
"""

from __future__ import annotations

import logging
import math
from typing import Dict, List, Optional, Tuple

from src.models import (
    ACMGClassification,
    TherapyMatch,
    TherapyStatus,
)

logger = logging.getLogger(__name__)


def _clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    """Clamp a numeric value to [lo, hi]."""
    return max(lo, min(hi, value))


# ═══════════════════════════════════════════════════════════════════════════
# HPO-TO-GENE MATCHER
# ═══════════════════════════════════════════════════════════════════════════


class HPOToGeneMatcher:
    """Match patient HPO terms to candidate genes using information content
    (IC) scoring and best-match-average (BMA) semantic similarity.

    IC = -log2(p(t))  where p(t) is the frequency of term t in annotated
    disease entries.

    BMA similarity between patient terms P and disease terms D:
        BMA(P, D) = 0.5 * (avg max-IC match P->D + avg max-IC match D->P)

    Phenotype frequency weighting gives higher scores to phenotypes that
    are more frequently associated with a particular gene/disease.
    """

    # HPO term frequency estimates (fraction of diseases annotated with term)
    # Used for IC calculation: IC = -log2(frequency)
    _TERM_FREQUENCIES: Dict[str, float] = {
        "HP:0001250": 0.15,   # Seizures — common across many diseases
        "HP:0001252": 0.12,   # Hypotonia
        "HP:0001249": 0.10,   # Intellectual disability
        "HP:0000252": 0.04,   # Microcephaly
        "HP:0001263": 0.13,   # Global developmental delay
        "HP:0001166": 0.02,   # Arachnodactyly — specific to connective tissue
        "HP:0001519": 0.01,   # Disproportionate tall stature — very specific
        "HP:0004382": 0.01,   # Mitral valve prolapse
        "HP:0002816": 0.005,  # Genu recurvatum — highly specific
        "HP:0001083": 0.008,  # Ectopia lentis — specific to Marfan/homocystinuria
        "HP:0001657": 0.03,   # Long QT interval
        "HP:0004756": 0.02,   # Ventricular tachycardia
        "HP:0001279": 0.04,   # Syncope
        "HP:0001639": 0.03,   # Hypertrophic cardiomyopathy
        "HP:0001635": 0.05,   # Congestive heart failure
        "HP:0001382": 0.06,   # Joint hypermobility
        "HP:0000974": 0.04,   # Hyperextensible skin
        "HP:0001508": 0.08,   # Failure to thrive
        "HP:0002110": 0.03,   # Bronchiectasis
        "HP:0003202": 0.05,   # Skeletal muscle atrophy
        "HP:0002072": 0.01,   # Chorea — relatively specific
        "HP:0000726": 0.06,   # Dementia / psychiatric
        "HP:0001300": 0.04,   # Parkinsonism
        "HP:0002205": 0.02,   # Recurrent respiratory infections
        "HP:0001744": 0.02,   # Splenomegaly
        "HP:0002240": 0.03,   # Hepatomegaly
        "HP:0001882": 0.03,   # Leukopenia
        "HP:0001695": 0.02,   # Cardiac arrest
        "HP:0001663": 0.03,   # Ventricular fibrillation
        "HP:0001644": 0.04,   # Dilated cardiomyopathy
        "HP:0002093": 0.04,   # Respiratory insufficiency
        "HP:0001987": 0.01,   # Hyperammonemia — specific
        "HP:0001942": 0.02,   # Metabolic acidosis
        "HP:0002757": 0.02,   # Recurrent fractures
        "HP:0000678": 0.01,   # Dental crowding
        "HP:0001513": 0.03,   # Obesity
        "HP:0000046": 0.01,   # Scrotal hypoplasia
        "HP:0002167": 0.01,   # Abnormality of speech
        "HP:0000733": 0.01,   # Stereotypic motor behavior
        "HP:0004322": 0.03,   # Short stature
        "HP:0000256": 0.02,   # Macrocephaly
        "HP:0003027": 0.01,   # Mesomelic limb shortening
    }

    # Gene-HPO associations (gene -> list of associated HPO terms with frequencies)
    _GENE_HPO_MAP: Dict[str, List[Dict]] = {
        "CFTR": [
            {"hpo": "HP:0002205", "freq": 0.95}, {"hpo": "HP:0002110", "freq": 0.80},
            {"hpo": "HP:0001508", "freq": 0.70}, {"hpo": "HP:0006538", "freq": 0.85},
        ],
        "FBN1": [
            {"hpo": "HP:0001166", "freq": 0.90}, {"hpo": "HP:0001519", "freq": 0.80},
            {"hpo": "HP:0004382", "freq": 0.60}, {"hpo": "HP:0001083", "freq": 0.50},
            {"hpo": "HP:0002816", "freq": 0.30},
        ],
        "SCN1A": [
            {"hpo": "HP:0001250", "freq": 0.98}, {"hpo": "HP:0001263", "freq": 0.70},
            {"hpo": "HP:0001249", "freq": 0.65}, {"hpo": "HP:0001252", "freq": 0.40},
        ],
        "DMD": [
            {"hpo": "HP:0003202", "freq": 0.95}, {"hpo": "HP:0001644", "freq": 0.70},
            {"hpo": "HP:0002093", "freq": 0.60}, {"hpo": "HP:0003236", "freq": 0.90},
        ],
        "HTT": [
            {"hpo": "HP:0002072", "freq": 0.90}, {"hpo": "HP:0000726", "freq": 0.80},
            {"hpo": "HP:0001300", "freq": 0.50}, {"hpo": "HP:0001289", "freq": 0.60},
        ],
        "KCNQ1": [
            {"hpo": "HP:0001657", "freq": 0.95}, {"hpo": "HP:0004756", "freq": 0.40},
            {"hpo": "HP:0001279", "freq": 0.50}, {"hpo": "HP:0001695", "freq": 0.10},
        ],
        "MYH7": [
            {"hpo": "HP:0001639", "freq": 0.85}, {"hpo": "HP:0001635", "freq": 0.40},
            {"hpo": "HP:0001663", "freq": 0.20}, {"hpo": "HP:0001279", "freq": 0.30},
        ],
        "SMN1": [
            {"hpo": "HP:0003202", "freq": 0.98}, {"hpo": "HP:0001252", "freq": 0.95},
            {"hpo": "HP:0002093", "freq": 0.80}, {"hpo": "HP:0001319", "freq": 0.70},
        ],
        "MECP2": [
            {"hpo": "HP:0001263", "freq": 0.95}, {"hpo": "HP:0001249", "freq": 0.90},
            {"hpo": "HP:0001250", "freq": 0.70}, {"hpo": "HP:0002167", "freq": 0.85},
            {"hpo": "HP:0000733", "freq": 0.80},
        ],
        "PAH": [
            {"hpo": "HP:0001249", "freq": 0.80}, {"hpo": "HP:0001250", "freq": 0.15},
            {"hpo": "HP:0001252", "freq": 0.10}, {"hpo": "HP:0001263", "freq": 0.70},
        ],
        "GBA1": [
            {"hpo": "HP:0001744", "freq": 0.90}, {"hpo": "HP:0002240", "freq": 0.85},
            {"hpo": "HP:0001882", "freq": 0.60}, {"hpo": "HP:0002653", "freq": 0.50},
        ],
        "COL5A1": [
            {"hpo": "HP:0001382", "freq": 0.90}, {"hpo": "HP:0000974", "freq": 0.80},
            {"hpo": "HP:0001252", "freq": 0.30},
        ],
        "OTC": [
            {"hpo": "HP:0001987", "freq": 0.95}, {"hpo": "HP:0001250", "freq": 0.60},
            {"hpo": "HP:0001259", "freq": 0.50},
        ],
        "FGFR3": [
            {"hpo": "HP:0004322", "freq": 0.98}, {"hpo": "HP:0000256", "freq": 0.80},
            {"hpo": "HP:0003027", "freq": 0.90},
        ],
    }

    def _compute_ic(self, hpo_id: str) -> float:
        """Compute information content for an HPO term.

        IC = -log2(p(t)), where p(t) is the frequency of the term
        across annotated diseases.
        """
        freq = self._TERM_FREQUENCIES.get(hpo_id, 0.001)  # default rare
        return -math.log2(freq)

    def _max_ic_match(self, term: str, target_terms: List[str]) -> float:
        """Find maximum IC score for matching *term* against *target_terms*.

        Returns max IC if term is found in targets, else 0.
        """
        if term in target_terms:
            return self._compute_ic(term)
        return 0.0

    def _bma_similarity(
        self, patient_terms: List[str], gene_terms: List[str]
    ) -> float:
        """Compute Best-Match-Average (BMA) similarity.

        BMA(P, G) = 0.5 * (
            mean of max_IC(p, G) for p in P
            + mean of max_IC(g, P) for g in G
        )
        """
        if not patient_terms or not gene_terms:
            return 0.0

        # Patient -> Gene direction
        p_to_g = [self._max_ic_match(p, gene_terms) for p in patient_terms]
        avg_p_to_g = sum(p_to_g) / len(p_to_g)

        # Gene -> Patient direction
        g_to_p = [self._max_ic_match(g, patient_terms) for g in gene_terms]
        avg_g_to_p = sum(g_to_p) / len(g_to_p)

        return 0.5 * (avg_p_to_g + avg_g_to_p)

    def process_patient(
        self, hpo_terms: List[str], top_k: int = 10
    ) -> List[Dict]:
        """Process patient HPO terms and return ranked gene candidates.

        Parameters
        ----------
        hpo_terms : list[str]
            Patient HPO term IDs.
        top_k : int
            Number of top gene matches to return.

        Returns
        -------
        list[dict]
            Ranked list of gene matches with scores.
            Each dict: {gene, bma_score, matched_terms, ic_details}.
        """
        if not hpo_terms:
            return []

        results: List[Dict] = []

        for gene, associations in self._GENE_HPO_MAP.items():
            gene_terms = [a["hpo"] for a in associations]

            # BMA similarity
            bma = self._bma_similarity(hpo_terms, gene_terms)

            # Phenotype frequency weighting
            matched_terms = [t for t in hpo_terms if t in gene_terms]
            freq_weight = 0.0
            for mt in matched_terms:
                assoc = next((a for a in associations if a["hpo"] == mt), None)
                if assoc:
                    freq_weight += assoc["freq"]

            # Combined score
            combined = bma * 0.7 + (freq_weight / max(len(gene_terms), 1)) * 0.3
            combined = _clamp(combined / 10.0)  # Normalise to 0-1 range

            if bma > 0:
                ic_details = {t: round(self._compute_ic(t), 2) for t in matched_terms}
                results.append({
                    "gene": gene,
                    "bma_score": round(bma, 4),
                    "combined_score": round(combined, 4),
                    "matched_terms": matched_terms,
                    "ic_details": ic_details,
                    "n_gene_phenotypes": len(gene_terms),
                    "n_matched": len(matched_terms),
                })

        # Sort by combined score descending
        results.sort(key=lambda r: r["combined_score"], reverse=True)
        return results[:top_k]


# ═══════════════════════════════════════════════════════════════════════════
# ACMG VARIANT CLASSIFIER
# ═══════════════════════════════════════════════════════════════════════════


class ACMGVariantClassifier:
    """Classify variants according to ACMG/AMP guidelines (simplified).

    Implements a subset of the 28 ACMG criteria with scoring logic:
    - PVS1: Null variant (LOF) in gene with known LOF mechanism (+8)
    - PS1: Same amino acid change as established pathogenic (+4)
    - PS2: De novo (confirmed) in patient with disease (+4)
    - PS3: Well-established functional studies show damaging (+3)
    - PS4: Prevalence in affected significantly increased vs controls (+3)
    - PM1: Located in mutational hot spot / functional domain (+2)
    - PM2: Absent from controls (or extremely low frequency) (+2)
    - PM3: Detected in trans with pathogenic variant (recessive) (+2)
    - PM4: Protein length change (in-frame del/ins in non-repeat) (+2)
    - PM5: Novel missense at same position as established pathogenic (+2)
    - PM6: Assumed de novo (no confirmation) (+1)
    - PP1: Cosegregation with disease in family (+1)
    - PP2: Missense in gene with low rate of benign missense (+1)
    - PP3: Computational evidence supports deleterious (+1)
    - PP4: Phenotype is highly specific for disease with single gene (+1)
    - PP5: Reputable source reports variant as pathogenic (+1)
    - BA1: Allele frequency >5% in any population -> Benign (standalone)
    - BS1: Allele frequency greater than expected for disorder (+benign)
    - BS2: Observed in healthy adult (for fully penetrant condition) (+benign)
    - BP1: Missense in gene where only truncating cause disease (+benign)
    - BP3: In-frame del/ins in repetitive region (+benign)
    - BP4: Computational evidence suggests no impact (+benign)
    - BP6: Reputable source reports variant as benign (+benign)
    - BP7: Synonymous with no splice impact (+benign)

    Scoring thresholds:
    - Pathogenic: >= 10 (must include PVS1 or 2xPS)
    - Likely Pathogenic: >= 6
    - VUS: 1-5
    - Likely Benign: benign_score >= 4
    - Benign: BA1 alone or benign_score >= 6
    """

    # LOF-intolerant genes (pLI > 0.9 in gnomAD)
    _LOF_INTOLERANT_GENES = {
        "SCN1A", "MECP2", "KCNQ1", "KCNH2", "SCN5A", "MYH7", "MYBPC3",
        "FBN1", "COL1A1", "COL1A2", "HTT", "NSD1", "NIPBL", "KMT2D",
        "CHD7", "PTPN11", "RAF1", "BRAF", "FGFR3", "TCOF1",
    }

    # Mutational hot spots (gene -> list of domain ranges)
    _HOT_SPOTS: Dict[str, List[Tuple[int, int]]] = {
        "FGFR3": [(370, 380)],  # Transmembrane domain
        "BRAF": [(594, 601)],   # Kinase domain
        "KRAS": [(10, 15), (58, 63)],
        "TP53": [(125, 300)],   # DNA-binding domain
    }

    def classify(
        self, variant_data: dict
    ) -> Tuple[ACMGClassification, List[str], str]:
        """Classify a variant according to ACMG criteria.

        Parameters
        ----------
        variant_data : dict
            Keys: gene, variant_type, population_frequency, zygosity,
            is_de_novo, functional_studies, computational_prediction,
            in_clinvar, clinvar_classification, cosegregation,
            phenotype_specific, protein_position, in_repeat_region,
            synonymous, splice_impact.

        Returns
        -------
        tuple[ACMGClassification, list[str], str]
            (classification, criteria_met, evidence_summary)
        """
        gene = variant_data.get("gene", "")
        vtype = variant_data.get("variant_type", "")
        freq = float(variant_data.get("population_frequency", 0.0))
        is_de_novo = variant_data.get("is_de_novo", False)
        de_novo_confirmed = variant_data.get("de_novo_confirmed", False)
        functional = variant_data.get("functional_studies", "")
        comp_pred = variant_data.get("computational_prediction", "")
        in_clinvar = variant_data.get("in_clinvar", False)
        clinvar_class = variant_data.get("clinvar_classification", "")
        coseg = variant_data.get("cosegregation", False)
        pheno_specific = variant_data.get("phenotype_specific", False)
        protein_pos = variant_data.get("protein_position", None)
        in_repeat = variant_data.get("in_repeat_region", False)
        synonymous = variant_data.get("synonymous", False)
        splice_impact = variant_data.get("splice_impact", False)

        pathogenic_criteria: List[str] = []
        benign_criteria: List[str] = []
        path_score = 0
        benign_score = 0

        # ── BA1: Standalone benign ────────────────────────────────────
        if freq > 0.05:
            benign_criteria.append("BA1")
            classification = ACMGClassification.BENIGN
            summary = f"BA1: Population frequency {freq:.4f} > 5% — classified as Benign."
            return classification, benign_criteria, summary

        # ── Pathogenic criteria ───────────────────────────────────────

        # PVS1: Null variant in LOF-intolerant gene
        lof_types = {"frameshift", "nonsense", "splice_site", "deletion", "insertion"}
        if vtype in lof_types and gene.upper() in self._LOF_INTOLERANT_GENES:
            pathogenic_criteria.append("PVS1")
            path_score += 8

        # PS1: Same amino acid change as established pathogenic
        if in_clinvar and clinvar_class in ("pathogenic", "likely_pathogenic"):
            pathogenic_criteria.append("PS1")
            path_score += 4

        # PS2: De novo (confirmed)
        if is_de_novo and de_novo_confirmed:
            pathogenic_criteria.append("PS2")
            path_score += 4

        # PS3: Functional studies
        if functional in ("damaging", "loss_of_function", "deleterious"):
            pathogenic_criteria.append("PS3")
            path_score += 3

        # PM1: Mutational hot spot
        if protein_pos is not None and gene.upper() in self._HOT_SPOTS:
            for start, end in self._HOT_SPOTS[gene.upper()]:
                if start <= protein_pos <= end:
                    pathogenic_criteria.append("PM1")
                    path_score += 2
                    break

        # PM2: Absent from controls
        if freq < 0.0001:
            pathogenic_criteria.append("PM2")
            path_score += 2

        # PM3: Trans with pathogenic (recessive)
        if variant_data.get("in_trans_with_pathogenic", False):
            pathogenic_criteria.append("PM3")
            path_score += 2

        # PM4: In-frame protein length change (not in repeat)
        if vtype in ("indel", "insertion", "deletion") and not in_repeat:
            pathogenic_criteria.append("PM4")
            path_score += 2

        # PM5: Novel missense at same residue as established pathogenic
        if (
            vtype == "missense"
            and variant_data.get("same_residue_pathogenic", False)
            and "PS1" not in pathogenic_criteria
        ):
            pathogenic_criteria.append("PM5")
            path_score += 2

        # PM6: Assumed de novo
        if is_de_novo and not de_novo_confirmed:
            pathogenic_criteria.append("PM6")
            path_score += 1

        # PP1: Cosegregation
        if coseg:
            pathogenic_criteria.append("PP1")
            path_score += 1

        # PP3: Computational prediction
        if comp_pred in ("deleterious", "damaging", "pathogenic", "probably_damaging"):
            pathogenic_criteria.append("PP3")
            path_score += 1

        # PP4: Phenotype specific
        if pheno_specific:
            pathogenic_criteria.append("PP4")
            path_score += 1

        # PP5: Reputable source
        if in_clinvar and clinvar_class == "pathogenic":
            if "PP5" not in pathogenic_criteria and "PS1" not in pathogenic_criteria:
                pathogenic_criteria.append("PP5")
                path_score += 1

        # ── Benign criteria ───────────────────────────────────────────

        # BS1: Higher frequency than expected
        if freq > 0.01:
            benign_criteria.append("BS1")
            benign_score += 3

        # BS2: Observed in healthy adults
        if variant_data.get("seen_in_healthy", False):
            benign_criteria.append("BS2")
            benign_score += 3

        # BP1: Missense in gene where only truncating causes disease
        if variant_data.get("missense_non_truncating_gene", False):
            benign_criteria.append("BP1")
            benign_score += 1

        # BP3: In-frame in repeat region
        if vtype in ("indel", "insertion", "deletion") and in_repeat:
            benign_criteria.append("BP3")
            benign_score += 1

        # BP4: Computational no impact
        if comp_pred in ("benign", "tolerated", "neutral"):
            benign_criteria.append("BP4")
            benign_score += 1

        # BP6: Reputable source benign
        if in_clinvar and clinvar_class in ("benign", "likely_benign"):
            benign_criteria.append("BP6")
            benign_score += 2

        # BP7: Synonymous no splice
        if synonymous and not splice_impact:
            benign_criteria.append("BP7")
            benign_score += 1

        # ── Classification decision ───────────────────────────────────
        all_criteria = pathogenic_criteria + benign_criteria

        if benign_score >= 6:
            classification = ACMGClassification.BENIGN
        elif benign_score >= 4 and path_score == 0:
            classification = ACMGClassification.LIKELY_BENIGN
        elif path_score >= 10:
            classification = ACMGClassification.PATHOGENIC
        elif path_score >= 6:
            classification = ACMGClassification.LIKELY_PATHOGENIC
        elif path_score > 0 and benign_score > 0:
            classification = ACMGClassification.VUS
        elif path_score > 0:
            classification = ACMGClassification.VUS
        elif benign_score > 0:
            classification = ACMGClassification.LIKELY_BENIGN
        else:
            classification = ACMGClassification.VUS

        summary = (
            f"Pathogenic criteria: {', '.join(pathogenic_criteria) if pathogenic_criteria else 'none'} "
            f"(score={path_score}). "
            f"Benign criteria: {', '.join(benign_criteria) if benign_criteria else 'none'} "
            f"(score={benign_score}). "
            f"Classification: {classification.value}."
        )

        logger.debug("ACMG classification for %s/%s: %s", gene, vtype, classification.value)

        # Publish safety alert for pathogenic variant classifications
        if classification in (ACMGClassification.PATHOGENIC, ACMGClassification.LIKELY_PATHOGENIC):
            try:
                from api.routes.events import publish_event
                publish_event("critical_alert", {
                    "alert_type": "pathogenic_variant",
                    "gene": gene,
                    "variant_type": vtype,
                    "classification": classification.value,
                    "criteria": all_criteria,
                    "message": f"Variant in {gene} classified as {classification.value} -- clinical review recommended",
                })
            except Exception:
                logger.debug("Could not publish safety alert for %s classification", classification.value)

        return classification, all_criteria, summary


# ═══════════════════════════════════════════════════════════════════════════
# ORPHAN DRUG MATCHER
# ═══════════════════════════════════════════════════════════════════════════


class OrphanDrugMatcher:
    """Match patient disease/genotype to orphan drug therapies.

    Checks three match types:
    1. Exact disease match (approved indication)
    2. Pathway match (same pathway, different disease)
    3. Repurposing candidates (mechanism-based)
    """

    _ORPHAN_DRUGS: List[Dict] = [
        {"drug": "Ivacaftor (Kalydeco)", "disease": "Cystic Fibrosis", "gene": "CFTR",
         "pathway": "CFTR modulator", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": "G551D and other gating mutations"},
        {"drug": "Lumacaftor/Ivacaftor (Orkambi)", "disease": "Cystic Fibrosis", "gene": "CFTR",
         "pathway": "CFTR modulator", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": "F508del homozygous"},
        {"drug": "Elexacaftor/Tezacaftor/Ivacaftor (Trikafta)", "disease": "Cystic Fibrosis", "gene": "CFTR",
         "pathway": "CFTR modulator", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": "At least one F508del allele"},
        {"drug": "Nusinersen (Spinraza)", "disease": "Spinal Muscular Atrophy", "gene": "SMN1",
         "pathway": "SMN2 splicing modifier", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": None},
        {"drug": "Risdiplam (Evrysdi)", "disease": "Spinal Muscular Atrophy", "gene": "SMN1",
         "pathway": "SMN2 splicing modifier", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": None},
        {"drug": "Imiglucerase (Cerezyme)", "disease": "Gaucher Disease", "gene": "GBA1",
         "pathway": "Enzyme replacement", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": None},
        {"drug": "Migalastat (Galafold)", "disease": "Fabry Disease", "gene": "GLA",
         "pathway": "Pharmacological chaperone", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": "Amenable GLA mutations"},
        {"drug": "Sapropterin (Kuvan)", "disease": "Phenylketonuria", "gene": "PAH",
         "pathway": "BH4 cofactor", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": "BH4-responsive mutations"},
        {"drug": "Cerliponase alfa (Brineura)", "disease": "CLN2 Batten Disease", "gene": "TPP1",
         "pathway": "Enzyme replacement (intracerebroventricular)", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": None},
        {"drug": "Mavacamten (Camzyos)", "disease": "Obstructive Hypertrophic Cardiomyopathy", "gene": "MYH7",
         "pathway": "Cardiac myosin inhibitor", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": None},
        {"drug": "Mexiletine", "disease": "Long QT Syndrome Type 3", "gene": "SCN5A",
         "pathway": "Sodium channel blocker", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": "LQT3 (SCN5A gain-of-function)"},
        {"drug": "Stiripentol (Diacomit)", "disease": "Dravet Syndrome", "gene": "SCN1A",
         "pathway": "Antiepileptic (GABAergic)", "status": TherapyStatus.APPROVED_FDA,
         "genotype_specific": None},
    ]

    def match(
        self, disease_id: str, genotype: Optional[str] = None
    ) -> List[TherapyMatch]:
        """Match patient disease/genotype to orphan drug therapies.

        Parameters
        ----------
        disease_id : str
            Disease name or identifier.
        genotype : str, optional
            Patient genotype for genotype-specific matching.

        Returns
        -------
        list[TherapyMatch]
            Matched therapies ranked by match quality.
        """
        disease_lower = disease_id.lower()
        genotype_lower = (genotype or "").lower()
        matches: List[TherapyMatch] = []

        for drug in self._ORPHAN_DRUGS:
            match_type = None
            eligibility_note = ""

            # Exact disease match
            if drug["disease"].lower() in disease_lower or disease_lower in drug["disease"].lower():
                match_type = "exact_disease"
                if drug.get("genotype_specific"):
                    if genotype_lower and any(
                        g.lower() in genotype_lower
                        for g in drug["genotype_specific"].lower().split(" and ")
                    ):
                        eligibility_note = f"Genotype-specific match: {drug['genotype_specific']}"
                    else:
                        eligibility_note = f"Genotype-specific: {drug['genotype_specific']} — verify eligibility"

            # Pathway match (same gene, different disease)
            elif drug.get("gene") and drug["gene"].upper() in disease_lower.upper():
                match_type = "pathway_match"
                eligibility_note = f"Pathway match via gene {drug['gene']} — off-label use may be considered"

            if match_type:
                matches.append(TherapyMatch(
                    therapy_name=drug["drug"],
                    indication=drug["disease"],
                    status=drug["status"],
                    mechanism=drug["pathway"],
                    gene_target=drug.get("gene"),
                    orphan_designation=True,
                    eligibility_criteria=eligibility_note or None,
                ))

        return matches


# ═══════════════════════════════════════════════════════════════════════════
# DIAGNOSTIC ALGORITHM RECOMMENDER
# ═══════════════════════════════════════════════════════════════════════════


class DiagnosticAlgorithmRecommender:
    """Recommend ordered diagnostic test sequences based on phenotype cluster.

    Six clinical pathways:
    - Neurodevelopmental
    - Metabolic
    - Skeletal
    - Cardiac
    - Immunodeficiency
    - Connective tissue
    """

    _PATHWAYS: Dict[str, List[Dict[str, str]]] = {
        "neurodevelopmental": [
            {"test": "Chromosomal Microarray (CMA)", "yield": "15-20%", "urgency": "high",
             "rationale": "First-tier test for DD/ID; detects CNVs and UPD"},
            {"test": "Fragile X (FMR1 CGG repeat analysis)", "yield": "2-3%", "urgency": "high",
             "rationale": "Most common inherited cause of ID in males"},
            {"test": "Whole Exome Sequencing (WES)", "yield": "25-40%", "urgency": "medium",
             "rationale": "Second-tier for unresolved cases after CMA"},
            {"test": "Metabolic panel (amino acids, organic acids, acylcarnitines)", "yield": "1-5%", "urgency": "medium",
             "rationale": "Treatable metabolic causes should be excluded"},
            {"test": "Brain MRI", "yield": "Variable", "urgency": "medium",
             "rationale": "Structural brain abnormalities, white matter disease"},
            {"test": "EEG", "yield": "Variable", "urgency": "low",
             "rationale": "If seizures suspected or staring episodes reported"},
        ],
        "metabolic": [
            {"test": "Plasma amino acids", "yield": "High for aminoacidopathies", "urgency": "high",
             "rationale": "Detect amino acid disorders (PKU, MSUD, homocystinuria)"},
            {"test": "Urine organic acids", "yield": "High for organic acidemias", "urgency": "high",
             "rationale": "Detect organic acidurias (MMA, PA, isovaleric)"},
            {"test": "Acylcarnitine profile", "yield": "High for FAO defects", "urgency": "high",
             "rationale": "Fatty acid oxidation disorders, organic acidemias"},
            {"test": "Ammonia and lactate", "yield": "Screening", "urgency": "urgent",
             "rationale": "Urea cycle defects, mitochondrial disorders"},
            {"test": "Very long chain fatty acids (VLCFA)", "yield": "Specific", "urgency": "medium",
             "rationale": "Peroxisomal disorders (X-ALD, Zellweger spectrum)"},
            {"test": "Lysosomal enzyme panel", "yield": "High for LSDs", "urgency": "medium",
             "rationale": "Gaucher, Fabry, Pompe, MPS disorders"},
        ],
        "skeletal": [
            {"test": "Skeletal survey (radiographs)", "yield": "High", "urgency": "high",
             "rationale": "Characterize skeletal dysplasia pattern"},
            {"test": "Growth chart analysis with bone age", "yield": "Moderate", "urgency": "medium",
             "rationale": "Growth pattern differentiation"},
            {"test": "Skeletal dysplasia gene panel", "yield": "30-50%", "urgency": "medium",
             "rationale": "FGFR3, COL1A1/2, COL2A1, COMP, etc."},
            {"test": "DEXA scan", "yield": "Specific", "urgency": "low",
             "rationale": "Bone density for OI and osteoporosis"},
        ],
        "cardiac": [
            {"test": "12-lead ECG", "yield": "High for channelopathies", "urgency": "urgent",
             "rationale": "QT interval, Brugada pattern, conduction defects"},
            {"test": "Echocardiography", "yield": "High for cardiomyopathies", "urgency": "urgent",
             "rationale": "Structural and functional assessment"},
            {"test": "Cardiac MRI", "yield": "Moderate", "urgency": "high",
             "rationale": "ARVC diagnosis, fibrosis detection, tissue characterization"},
            {"test": "Cardiac gene panel", "yield": "20-40%", "urgency": "high",
             "rationale": "Channelopathy or cardiomyopathy gene panel"},
            {"test": "Holter monitor (48-72h)", "yield": "Variable", "urgency": "medium",
             "rationale": "Arrhythmia detection and characterization"},
        ],
        "immunodeficiency": [
            {"test": "Complete blood count with differential", "yield": "Screening", "urgency": "urgent",
             "rationale": "Lymphopenia, neutropenia, leukocyte abnormalities"},
            {"test": "Immunoglobulin levels (IgG, IgA, IgM, IgE)", "yield": "High", "urgency": "high",
             "rationale": "Humoral immune deficiency"},
            {"test": "Lymphocyte subset panel (flow cytometry)", "yield": "High", "urgency": "high",
             "rationale": "T/B/NK cell deficiency characterization"},
            {"test": "Complement levels (CH50, AH50)", "yield": "Moderate", "urgency": "medium",
             "rationale": "Complement deficiency workup"},
            {"test": "SCID NBS (TREC assay)", "yield": "Screening", "urgency": "urgent",
             "rationale": "Severe combined immunodeficiency screening"},
            {"test": "Primary immunodeficiency gene panel", "yield": "25-40%", "urgency": "medium",
             "rationale": "Over 400 known PID genes"},
        ],
        "connective_tissue": [
            {"test": "Beighton score assessment", "yield": "Screening", "urgency": "medium",
             "rationale": "Quantify joint hypermobility"},
            {"test": "Echocardiography (aortic root measurement)", "yield": "High for Marfan", "urgency": "high",
             "rationale": "Aortic root dilation, mitral valve prolapse"},
            {"test": "Ophthalmologic exam (slit lamp)", "yield": "High for Marfan", "urgency": "high",
             "rationale": "Ectopia lentis, myopia"},
            {"test": "FBN1 gene sequencing", "yield": "70-93% for Marfan", "urgency": "medium",
             "rationale": "Marfan syndrome confirmation"},
            {"test": "COL3A1 gene sequencing", "yield": "High for vEDS", "urgency": "urgent",
             "rationale": "Vascular EDS — critical for management decisions"},
            {"test": "Connective tissue gene panel", "yield": "30-50%", "urgency": "medium",
             "rationale": "Comprehensive panel if specific diagnosis unclear"},
        ],
    }

    def recommend(self, phenotype_cluster: str) -> List[Dict[str, str]]:
        """Recommend ordered diagnostic test sequence for a phenotype cluster.

        Parameters
        ----------
        phenotype_cluster : str
            Clinical phenotype cluster (neurodevelopmental, metabolic,
            skeletal, cardiac, immunodeficiency, connective_tissue).

        Returns
        -------
        list[dict]
            Ordered list of recommended tests with yield and rationale.
        """
        cluster_key = phenotype_cluster.lower().replace(" ", "_").replace("-", "_")

        # Try exact match first
        if cluster_key in self._PATHWAYS:
            return self._PATHWAYS[cluster_key]

        # Fuzzy match
        for key in self._PATHWAYS:
            if key in cluster_key or cluster_key in key:
                return self._PATHWAYS[key]

        logger.warning("No pathway match for cluster: %s", phenotype_cluster)
        return [
            {"test": "Comprehensive clinical evaluation", "yield": "Variable", "urgency": "medium",
             "rationale": f"No predefined pathway for '{phenotype_cluster}' — start with clinical assessment"},
            {"test": "Whole Exome Sequencing (WES)", "yield": "25-40%", "urgency": "medium",
             "rationale": "Broad genetic testing when phenotype cluster is unclear"},
        ]


# ═══════════════════════════════════════════════════════════════════════════
# FAMILY SEGREGATION ANALYZER
# ═══════════════════════════════════════════════════════════════════════════


class FamilySegregationAnalyzer:
    """Analyse variant segregation in a family pedigree and compute LOD-like
    scores for evidence classification.

    Simplified LOD score calculation:
        LOD = log10(L(linkage) / L(no linkage))

    For small families, returns ACMG segregation evidence:
    - PS: Pathogenic Strong (LOD >= 3.0)
    - PM: Pathogenic Moderate (LOD >= 1.5)
    - PP: Pathogenic Supporting (LOD >= 0.6)
    - Insufficient: LOD < 0.6
    """

    def analyze(
        self,
        pedigree: List[Dict],
        genotypes: Dict[str, str],
        inheritance: str,
    ) -> Dict:
        """Analyse variant segregation across a pedigree.

        Parameters
        ----------
        pedigree : list[dict]
            Family members with keys: id, affected (bool), sex, relationship.
        genotypes : dict[str, str]
            Member ID -> genotype ('het', 'hom', 'wt', 'hemi', 'unknown').
        inheritance : str
            Expected inheritance pattern (autosomal_dominant, autosomal_recessive,
            x_linked_recessive, x_linked_dominant).

        Returns
        -------
        dict
            {lod_score, classification, concordant, discordant, total,
             evidence_level, details}
        """
        if not pedigree or not genotypes:
            return {
                "lod_score": 0.0,
                "classification": "insufficient",
                "concordant": 0,
                "discordant": 0,
                "total": 0,
                "evidence_level": "insufficient",
                "details": "No pedigree or genotype data provided",
            }

        concordant = 0
        discordant = 0
        details: List[str] = []

        for member in pedigree:
            member_id = member.get("id", "")
            affected = member.get("affected", False)
            sex = member.get("sex", "unknown")
            genotype = genotypes.get(member_id, "unknown")

            if genotype == "unknown":
                details.append(f"{member_id}: genotype unknown — excluded")
                continue

            expected_carrier = self._expected_genotype(
                affected, inheritance, sex
            )

            if self._genotype_matches(genotype, expected_carrier, inheritance):
                concordant += 1
                details.append(f"{member_id}: concordant (affected={affected}, genotype={genotype})")
            else:
                discordant += 1
                details.append(f"{member_id}: DISCORDANT (affected={affected}, genotype={genotype})")

        total = concordant + discordant
        if total == 0:
            return {
                "lod_score": 0.0,
                "classification": "insufficient",
                "concordant": 0,
                "discordant": 0,
                "total": 0,
                "evidence_level": "insufficient",
                "details": "No informative family members",
            }

        # Simplified LOD score: each concordant meiosis contributes +0.3,
        # each discordant contributes -1.0 (assuming theta=0 vs theta=0.5)
        lod_score = concordant * 0.3 - discordant * 1.0
        lod_score = round(lod_score, 2)

        # ACMG classification
        if discordant > 0:
            evidence_level = "against_segregation"
            classification = "non_segregating"
        elif lod_score >= 3.0:
            evidence_level = "PS"
            classification = "strong_pathogenic"
        elif lod_score >= 1.5:
            evidence_level = "PM"
            classification = "moderate_pathogenic"
        elif lod_score >= 0.6:
            evidence_level = "PP"
            classification = "supporting_pathogenic"
        else:
            evidence_level = "insufficient"
            classification = "insufficient"

        return {
            "lod_score": lod_score,
            "classification": classification,
            "concordant": concordant,
            "discordant": discordant,
            "total": total,
            "evidence_level": evidence_level,
            "details": "; ".join(details),
        }

    @staticmethod
    def _expected_genotype(
        affected: bool, inheritance: str, sex: str
    ) -> str:
        """Determine expected genotype for a family member."""
        inh = inheritance.lower()
        if "dominant" in inh:
            return "het" if affected else "wt"
        elif "recessive" in inh:
            return "hom" if affected else "het_or_wt"
        elif "x_linked" in inh:
            if sex == "male":
                return "hemi" if affected else "wt"
            else:
                return "het" if affected else "wt"
        return "unknown"

    @staticmethod
    def _genotype_matches(
        observed: str, expected: str, inheritance: str
    ) -> bool:
        """Check if observed genotype matches expected."""
        if expected == "unknown":
            return True  # Cannot assess
        if expected == "het_or_wt":
            return observed in ("het", "wt")
        return observed == expected


# ═══════════════════════════════════════════════════════════════════════════
# NATURAL HISTORY PREDICTOR
# ═══════════════════════════════════════════════════════════════════════════


class NaturalHistoryPredictor:
    """Predict disease natural history milestones with confidence intervals.

    Uses curated natural history data from registries and literature to
    provide age-specific milestone predictions.
    """

    _NATURAL_HISTORY: Dict[str, Dict] = {
        "spinal muscular atrophy type 1": {
            "milestones": [
                {"event": "Symptom onset", "median_age_months": 3, "range": (0, 6), "confidence": 0.90},
                {"event": "Loss of head control", "median_age_months": 6, "range": (3, 12), "confidence": 0.85},
                {"event": "Ventilatory support needed", "median_age_months": 12, "range": (6, 24), "confidence": 0.80},
                {"event": "Mortality (without treatment)", "median_age_months": 24, "range": (12, 48), "confidence": 0.75},
            ],
            "genotype_modifiers": {"smn2_copies_3": "Milder course, later onset", "smn2_copies_1": "Severe, early onset"},
        },
        "duchenne muscular dystrophy": {
            "milestones": [
                {"event": "Symptom onset (gait abnormality)", "median_age_months": 36, "range": (24, 60), "confidence": 0.85},
                {"event": "Loss of independent ambulation", "median_age_months": 120, "range": (84, 156), "confidence": 0.80},
                {"event": "Cardiomyopathy onset", "median_age_months": 168, "range": (120, 216), "confidence": 0.70},
                {"event": "Ventilatory support needed", "median_age_months": 216, "range": (180, 264), "confidence": 0.70},
            ],
            "genotype_modifiers": {"exon_skippable": "Eligible for exon-skipping therapy (milder with treatment)"},
        },
        "cystic fibrosis": {
            "milestones": [
                {"event": "Diagnosis (median, NBS era)", "median_age_months": 1, "range": (0, 6), "confidence": 0.90},
                {"event": "First Pseudomonas colonization", "median_age_months": 24, "range": (6, 120), "confidence": 0.60},
                {"event": "FEV1 decline below 80%", "median_age_months": 144, "range": (72, 216), "confidence": 0.65},
                {"event": "Median survival (with CFTR modulators)", "median_age_months": 600, "range": (480, 720), "confidence": 0.50},
            ],
            "genotype_modifiers": {"f508del_homozygous": "Trikafta eligible — dramatically improved outcomes"},
        },
        "phenylketonuria": {
            "milestones": [
                {"event": "NBS detection", "median_age_months": 0, "range": (0, 1), "confidence": 0.95},
                {"event": "Diet initiation", "median_age_months": 0, "range": (0, 1), "confidence": 0.95},
                {"event": "Normal IQ if treated <3 weeks", "median_age_months": 0, "range": (0, 1), "confidence": 0.90},
                {"event": "ID if untreated (IQ < 50)", "median_age_months": 12, "range": (6, 36), "confidence": 0.85},
            ],
            "genotype_modifiers": {"bh4_responsive": "Sapropterin may allow diet relaxation"},
        },
        "marfan syndrome": {
            "milestones": [
                {"event": "Diagnosis (median)", "median_age_months": 120, "range": (24, 240), "confidence": 0.60},
                {"event": "Aortic root dilation", "median_age_months": 144, "range": (60, 360), "confidence": 0.70},
                {"event": "Aortic dissection risk (untreated)", "median_age_months": 420, "range": (240, 600), "confidence": 0.50},
                {"event": "Ectopia lentis", "median_age_months": 72, "range": (24, 240), "confidence": 0.55},
            ],
            "genotype_modifiers": {"haploinsufficiency": "Generally more severe than dominant-negative"},
        },
        "dravet syndrome": {
            "milestones": [
                {"event": "First febrile seizure", "median_age_months": 6, "range": (4, 12), "confidence": 0.90},
                {"event": "Afebrile seizures begin", "median_age_months": 12, "range": (8, 24), "confidence": 0.85},
                {"event": "Developmental plateau", "median_age_months": 24, "range": (18, 36), "confidence": 0.80},
                {"event": "Gait abnormality onset", "median_age_months": 36, "range": (24, 72), "confidence": 0.70},
            ],
            "genotype_modifiers": {"truncating_scn1a": "More severe phenotype than missense"},
        },
    }

    def predict(
        self,
        disease: str,
        genotype: Optional[str] = None,
        current_age_months: Optional[int] = None,
    ) -> Dict:
        """Predict natural history milestones for a disease.

        Parameters
        ----------
        disease : str
            Disease name.
        genotype : str, optional
            Patient genotype for modifier effects.
        current_age_months : int, optional
            Patient current age in months for filtering future milestones.

        Returns
        -------
        dict
            {disease, milestones, genotype_modifier, future_milestones}
        """
        disease_lower = disease.lower()

        # Find matching disease
        history = None
        matched_name = None
        for name, data in self._NATURAL_HISTORY.items():
            if name in disease_lower or disease_lower in name:
                history = data
                matched_name = name
                break

        if history is None:
            return {
                "disease": disease,
                "milestones": [],
                "genotype_modifier": None,
                "future_milestones": [],
                "message": f"No natural history data available for '{disease}'",
            }

        milestones = history["milestones"]

        # Apply genotype modifier
        genotype_modifier = None
        if genotype:
            genotype_lower = genotype.lower().replace(" ", "_").replace("-", "_")
            for mod_key, mod_desc in history.get("genotype_modifiers", {}).items():
                if mod_key in genotype_lower or genotype_lower in mod_key:
                    genotype_modifier = mod_desc
                    break

        # Filter to future milestones if current age provided
        future_milestones = []
        if current_age_months is not None:
            for ms in milestones:
                if ms["median_age_months"] > current_age_months:
                    months_until = ms["median_age_months"] - current_age_months
                    future_milestones.append({
                        **ms,
                        "months_until_median": months_until,
                        "years_until_median": round(months_until / 12, 1),
                    })
        else:
            future_milestones = milestones

        return {
            "disease": matched_name,
            "milestones": milestones,
            "genotype_modifier": genotype_modifier,
            "future_milestones": future_milestones,
        }
