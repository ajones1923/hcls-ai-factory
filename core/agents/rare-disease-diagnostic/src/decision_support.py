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
import re
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


def _normalise_disease(name: str) -> str:
    """Normalise a disease name for *equality* comparison.

    Case, punctuation and whitespace are folded; nothing else. This exists so
    disease lookups can require a whole-name match instead of the bidirectional
    substring test that returned mavacamten for "cardiomyopathy" and the
    Duchenne natural history for "muscular dystrophy".
    """
    return re.sub(r"[^a-z0-9]+", " ", (name or "").lower()).strip()


# ═══════════════════════════════════════════════════════════════════════════
# ACMG/AMP TABLE 5 COMBINING RULES
# ═══════════════════════════════════════════════════════════════════════════
#
# Richards S, Aziz N, Bale S, et al. "Standards and guidelines for the
# interpretation of sequence variants: a joint consensus recommendation of the
# American College of Medical Genetics and Genomics and the Association for
# Molecular Pathology." Genet Med 2015;17(5):405-424. PMC4544753.
# Table 5, "Rules for combining criteria to classify sequence variants",
# retrieved 2026-11-03 from https://pmc.ncbi.nlm.nih.gov/articles/PMC4544753/
# and transcribed verbatim into the predicates below.
#
# These rule sets replace a point-scoring approximation (PVS1=8, PS=4/3,
# PM=2, PP=1, with Pathogenic at >=10) that disagreed with Table 5 in both
# directions -- it called PVS1 alone Likely pathogenic and PS3+PM2 a VUS.
# Table 5 is combinatorial over strength BUCKETS, not additive over weights,
# so the weights are gone rather than retuned.

#: Strength bucket for each ACMG/AMP criterion code (Richards 2015 Tables 3/4).
_ACMG_BUCKETS: Dict[str, str] = {
    # Pathogenic: Very Strong / Strong / Moderate / Supporting
    "PVS1": "PVS",
    "PS1": "PS", "PS2": "PS", "PS3": "PS", "PS4": "PS",
    "PM1": "PM", "PM2": "PM", "PM3": "PM", "PM4": "PM", "PM5": "PM",
    "PM6": "PM",
    "PP1": "PP", "PP2": "PP", "PP3": "PP", "PP4": "PP", "PP5": "PP",
    # Benign: Stand-alone / Strong / Supporting
    "BA1": "BA",
    "BS1": "BS", "BS2": "BS", "BS3": "BS", "BS4": "BS",
    "BP1": "BP", "BP2": "BP", "BP3": "BP", "BP4": "BP", "BP5": "BP",
    "BP6": "BP", "BP7": "BP",
}


def _bucket_counts(criteria: List[str]) -> Dict[str, int]:
    """Count criterion codes by Richards 2015 strength bucket."""
    counts = {b: 0 for b in ("PVS", "PS", "PM", "PP", "BA", "BS", "BP")}
    for code in criteria:
        bucket = _ACMG_BUCKETS.get(code)
        if bucket is not None:
            counts[bucket] += 1
    return counts


def _pathogenic_rule(pvs: int, ps: int, pm: int, pp: int) -> Optional[str]:
    """Return the Table 5 *Pathogenic* rule that fires, or None."""
    if pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2):
        return "1 Very Strong (PVS1) AND (>=1 Strong OR >=2 Moderate OR 1 Moderate+1 Supporting OR >=2 Supporting)"
    if ps >= 2:
        return ">=2 Strong (PS1-PS4)"
    if ps >= 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)):
        return "1 Strong AND (>=3 Moderate OR 2 Moderate+>=2 Supporting OR 1 Moderate+>=4 Supporting)"
    return None


def _likely_pathogenic_rule(pvs: int, ps: int, pm: int, pp: int) -> Optional[str]:
    """Return the Table 5 *Likely pathogenic* rule that fires, or None."""
    if pvs >= 1 and pm >= 1:
        return "1 Very Strong (PVS1) AND 1 Moderate"
    if ps >= 1 and 1 <= pm <= 2:
        return "1 Strong AND 1-2 Moderate"
    if ps >= 1 and pp >= 2:
        return "1 Strong AND >=2 Supporting"
    if pm >= 3:
        return ">=3 Moderate"
    if pm >= 2 and pp >= 2:
        return "2 Moderate AND >=2 Supporting"
    if pm >= 1 and pp >= 4:
        return "1 Moderate AND >=4 Supporting"
    return None


def _benign_rule(ba: int, bs: int) -> Optional[str]:
    """Return the Table 5 *Benign* rule that fires, or None."""
    if ba >= 1:
        return "1 Stand-Alone (BA1)"
    if bs >= 2:
        return ">=2 Strong (BS1-BS4)"
    return None


def _likely_benign_rule(bs: int, bp: int) -> Optional[str]:
    """Return the Table 5 *Likely benign* rule that fires, or None."""
    if bs >= 1 and bp >= 1:
        return "1 Strong (BS1-BS4) AND 1 Supporting (BP1-BP7)"
    if bp >= 2:
        return ">=2 Supporting (BP1-BP7)"
    return None


def combine_acmg_criteria(
    criteria: List[str],
) -> Tuple[ACMGClassification, str]:
    """Combine met ACMG/AMP criteria into a classification per Richards Table 5.

    Parameters
    ----------
    criteria : list[str]
        Criterion codes that were met, e.g. ``["PVS1", "PM2"]``.

    Returns
    -------
    tuple[ACMGClassification, str]
        ``(classification, rule_that_fired)``. The rule string makes the call
        auditable against Table 5.

    Notes
    -----
    Richards 2015 Table 5 footnote: *"Variants should be classified as
    Uncertain Significance if other criteria are unmet or the criteria for
    benign and pathogenic are contradictory"*. Both halves of that sentence
    are implemented here -- a variant satisfying a pathogenic rule set AND a
    benign rule set returns VUS rather than whichever side was tested first.
    """
    c = _bucket_counts(criteria)
    pvs, ps, pm, pp = c["PVS"], c["PS"], c["PM"], c["PP"]
    ba, bs, bp = c["BA"], c["BS"], c["BP"]

    path_rule = _pathogenic_rule(pvs, ps, pm, pp)
    lp_rule = _likely_pathogenic_rule(pvs, ps, pm, pp)
    ben_rule = _benign_rule(ba, bs)
    lb_rule = _likely_benign_rule(bs, bp)

    path_side = path_rule or lp_rule
    benign_side = ben_rule or lb_rule

    # Table 5 footnote: contradictory benign and pathogenic criteria -> VUS.
    if path_side and benign_side:
        return (
            ACMGClassification.VUS,
            f"contradictory evidence (pathogenic rule '{path_side}' and "
            f"benign rule '{benign_side}' both satisfied) -> Uncertain significance",
        )
    if path_rule:
        return ACMGClassification.PATHOGENIC, path_rule
    if lp_rule:
        return ACMGClassification.LIKELY_PATHOGENIC, lp_rule
    if ben_rule:
        return ACMGClassification.BENIGN, ben_rule
    if lb_rule:
        return ACMGClassification.LIKELY_BENIGN, lb_rule
    # Table 5 footnote: other criteria unmet -> Uncertain significance.
    return (
        ACMGClassification.VUS,
        "no Table 5 rule set satisfied -> Uncertain significance",
    )


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

    def _compute_ic(self, hpo_id: str) -> Optional[float]:
        """Compute information content for an HPO term.

        IC = -log2(p(t)), where p(t) is the frequency of the term
        across annotated diseases.

        Returns ``None`` for a term with no frequency on file. An unannotated
        term has *undefined* IC, not maximal IC: the previous ``0.001`` default
        scored any unknown string -- including ``"banana"`` and six terms from
        this class's own gene map -- at IC 9.97, nearly 4x the IC of Seizures
        (2.74), so the matcher's highest-weighted evidence was systematically
        its least-known evidence.
        """
        freq = self._TERM_FREQUENCIES.get(hpo_id)
        if freq is None:
            return None
        return -math.log2(freq)

    def _max_ic_match(self, term: str, target_terms: List[str]) -> float:
        """Find maximum IC score for matching *term* against *target_terms*.

        Returns the term's IC if it is found in the targets AND its frequency
        is on file, else 0.0 (an unscoreable term contributes nothing).
        """
        if term in target_terms:
            ic = self._compute_ic(term)
            if ic is not None:
                return ic
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
                ic_details: Dict[str, Optional[float]] = {}
                unscored: List[str] = []
                for t in matched_terms:
                    ic = self._compute_ic(t)
                    if ic is None:
                        ic_details[t] = None
                        unscored.append(t)
                    else:
                        ic_details[t] = round(ic, 2)
                results.append({
                    "gene": gene,
                    "bma_score": round(bma, 4),
                    "combined_score": round(combined, 4),
                    "matched_terms": matched_terms,
                    "ic_details": ic_details,
                    "unscored_terms": unscored,
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

    Criteria are combined by :func:`combine_acmg_criteria`, which implements
    Richards 2015 Table 5 verbatim. There are no criterion weights: a point
    system cannot reproduce Table 5 and the previous one did not.
    """

    #: Null variant types that may support PVS1 (Richards 2015 Table 3: "nonsense,
    #: frameshift, canonical +/-1 or 2 splice sites, initiation codon, single or
    #: multi-exon deletion"). Disjoint from _INFRAME_TYPES by construction.
    _NULL_TYPES = frozenset({
        "frameshift", "nonsense", "splice_site", "canonical_splice_site",
        "start_loss", "initiation_codon", "whole_exon_deletion",
        "multi_exon_deletion", "single_exon_deletion",
    })

    #: In-frame length changes that may support PM4/BP3 (Richards 2015 Table 3
    #: PM4: "Protein length changes due to in-frame deletions/insertions in a
    #: non-repeat region or stop-loss variants").
    _INFRAME_TYPES = frozenset({
        "inframe_deletion", "inframe_insertion", "inframe_indel", "stop_loss",
    })

    #: Variant-type labels that do not state whether the frame is preserved.
    #: PVS1 and PM4 are mutually exclusive by definition (a null variant is not
    #: an in-frame length change), so neither is applied to an ambiguous label
    #: and the ambiguity is reported.
    _FRAME_AMBIGUOUS_TYPES = frozenset({
        "deletion", "insertion", "indel", "duplication", "del", "ins",
    })

    #: ClinGen SVI "Updated Recommendation for the Benign Stand-Alone ACMG/AMP
    #: Criterion" -- Ghosh R, Harrison SM, Rehm HL, Plon SE, Biesecker LG,
    #: Hum Mutat 2018;39(11):1525-1530, PMC6188666, Table 1: variants to which
    #: BA1 must NOT be applied. Retrieved 2026-11-03 from
    #: https://pmc.ncbi.nlm.nih.gov/articles/PMC6188666/
    #: Keyed by (gene, HGVS c.) with the protein change recorded for the report.
    _BA1_EXCEPTIONS: Dict[Tuple[str, str], str] = {
        ("ACAD9", "C.-44_-41DUPTAAG"): "ACAD9 c.-44_-41dupTAAG",
        ("GJB2", "C.109G>A"): "GJB2 c.109G>A p.(Val37Ile)",
        ("HFE", "C.187C>G"): "HFE c.187C>G p.(His63Asp)",
        ("HFE", "C.845G>A"): "HFE c.845G>A p.(Cys282Tyr)",
        ("MEFV", "C.1105C>T"): "MEFV c.1105C>T p.(Pro369Ser)",
        ("MEFV", "C.1223G>A"): "MEFV c.1223G>A p.(Arg408Gln)",
        ("PIBF1", "C.1214G>A"): "PIBF1 c.1214G>A p.(Arg405Gln)",
        ("ACADS", "C.511C>T"): "ACADS c.511C>T p.(Arg171Trp)",
        ("BTD", "C.1330G>C"): "BTD c.1330G>C p.(Asp444His)",
    }

    #: Protein-level aliases for the same nine variants, so a caller supplying
    #: p. notation instead of c. notation still hits the exception list.
    _BA1_EXCEPTIONS_PROTEIN: Dict[Tuple[str, str], str] = {
        ("GJB2", "P.VAL37ILE"): "GJB2 c.109G>A p.(Val37Ile)",
        ("GJB2", "P.V37I"): "GJB2 c.109G>A p.(Val37Ile)",
        ("HFE", "P.HIS63ASP"): "HFE c.187C>G p.(His63Asp)",
        ("HFE", "P.H63D"): "HFE c.187C>G p.(His63Asp)",
        ("HFE", "P.CYS282TYR"): "HFE c.845G>A p.(Cys282Tyr)",
        ("HFE", "P.C282Y"): "HFE c.845G>A p.(Cys282Tyr)",
        ("MEFV", "P.PRO369SER"): "MEFV c.1105C>T p.(Pro369Ser)",
        ("MEFV", "P.ARG408GLN"): "MEFV c.1223G>A p.(Arg408Gln)",
        ("PIBF1", "P.ARG405GLN"): "PIBF1 c.1214G>A p.(Arg405Gln)",
        ("ACADS", "P.ARG171TRP"): "ACADS c.511C>T p.(Arg171Trp)",
        ("BTD", "P.ASP444HIS"): "BTD c.1330G>C p.(Asp444His)",
        ("BTD", "P.D444H"): "BTD c.1330G>C p.(Asp444His)",
    }

    #: ClinGen SVI BA1 specification: allele frequency >0.05 "in any general
    #: continental population dataset of at least 2,000 observed alleles"
    #: (Ghosh et al. 2018, PMC6188666).
    _BA1_MIN_OBSERVED_ALLELES = 2000

    # Mutational hot spots (gene -> list of domain ranges)
    _HOT_SPOTS: Dict[str, List[Tuple[int, int]]] = {
        "FGFR3": [(370, 380)],  # Transmembrane domain
        "BRAF": [(594, 601)],   # Kinase domain
        "KRAS": [(10, 15), (58, 63)],
        "TP53": [(125, 300)],   # DNA-binding domain
    }

    @classmethod
    def _ba1_exception(cls, variant_data: dict, gene: str) -> Optional[str]:
        """Return the ClinGen SVI BA1 exception label if this variant is one.

        Ghosh et al. 2018 (PMC6188666) Table 1 lists nine variants to which BA1
        must not be applied -- common, low-penetrance variants the frequency
        criterion was not designed for, including *HFE* c.845G>A p.(Cys282Tyr),
        the haemochromatosis variant that reaches ~6% in non-Finnish Europeans.
        """
        gene_u = (gene or "").upper()
        if not gene_u:
            return None
        for key in ("hgvs_c", "cdna_change", "variant_id", "hgvs"):
            value = variant_data.get(key)
            if value:
                hit = cls._BA1_EXCEPTIONS.get((gene_u, str(value).upper()))
                if hit:
                    return hit
        for key in ("hgvs_p", "protein_change", "amino_acid_change"):
            value = variant_data.get(key)
            if value:
                hit = cls._BA1_EXCEPTIONS_PROTEIN.get(
                    (gene_u, str(value).upper())
                )
                if hit:
                    return hit
        return None

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

            Evidence that ACMG requires explicitly and that cannot be inferred:

            ``lof_is_mechanism`` : bool
                Whether loss of function is a known disease mechanism for this
                gene. Required for PVS1 (Richards 2015 Table 3). Absent -> PVS1
                is not applied.
            ``same_amino_acid_as_known_pathogenic`` : bool
                Required for PS1. Presence in ClinVar is PP5, not PS1.
            ``allele_count_observed`` : int
                Observed alleles in the dataset reporting
                ``population_frequency``. Required for BA1 (ClinGen SVI
                minimum 2,000). Absent -> BA1 is not applied.
            ``hgvs_c`` / ``hgvs_p`` : str
                Used to check the ClinGen SVI BA1 exception list.

            ``population_frequency`` absent (``None``) means "not looked up",
            which grants neither PM2 nor BA1/BS1. It is not the same as 0.0.

        Returns
        -------
        tuple[ACMGClassification, list[str], str]
            (classification, criteria_met, evidence_summary)
        """
        gene = variant_data.get("gene", "")
        vtype = variant_data.get("variant_type", "")
        # A missing allele frequency is NOT 0.0. "The field was not supplied"
        # and "gnomAD was queried and the variant is absent" are different
        # facts; collapsing them granted PM2 -- the criterion whose purpose is
        # to prove the variant was looked for -- for not looking.
        raw_freq = variant_data.get("population_frequency", None)
        freq: Optional[float] = None if raw_freq is None else float(raw_freq)
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
        #: Criteria that were NOT applied, and why. These are reported so a
        #: reviewer can see which evidence was unavailable rather than inferring
        #: it was absent.
        not_applied: List[str] = []

        # ── BA1: Stand-alone benign, per the ClinGen SVI specification ────
        # Ghosh et al. 2018 (PMC6188666): "Allele frequency is >0.05 in any
        # general continental population dataset of at least 2,000 observed
        # alleles and found in a gene without a gene- or variant-specific BA1
        # modification." All three conditions are now required, and BA1 no
        # longer short-circuits the pathogenic evidence -- Table 5's
        # contradiction footnote resolves the conflict instead.
        ba1_exception = self._ba1_exception(variant_data, gene)
        observed_alleles = variant_data.get("allele_count_observed", None)
        if freq is None:
            not_applied.append(
                "BA1/BS1/PM2 not applied: population frequency not available"
            )
        elif freq > 0.05:
            if ba1_exception is not None:
                not_applied.append(
                    f"BA1 not applied: {ba1_exception} is on the ClinGen SVI "
                    "BA1 exception list (Ghosh et al. 2018, Table 1)"
                )
            elif observed_alleles is None:
                not_applied.append(
                    "BA1 not applied: the number of observed alleles in the "
                    "reporting dataset was not supplied, so the ClinGen SVI "
                    ">=2,000-allele requirement cannot be checked"
                )
            elif int(observed_alleles) < self._BA1_MIN_OBSERVED_ALLELES:
                not_applied.append(
                    f"BA1 not applied: only {int(observed_alleles)} observed "
                    f"alleles, below the ClinGen SVI minimum of "
                    f"{self._BA1_MIN_OBSERVED_ALLELES}"
                )
            else:
                benign_criteria.append("BA1")

        # ── Pathogenic criteria ───────────────────────────────────────

        # PVS1: Richards 2015 Table 3 -- "Null variant ... **in a gene where
        # loss of function (LOF) is a known mechanism of disease**". The gene
        # clause is the whole criterion, and it is a per-gene disease-mechanism
        # assertion, not a gene-level constraint metric: the previous pLI list
        # granted PVS1 to truncating FGFR3 (achondroplasia is gain-of-function),
        # HTT (CAG expansion), and BRAF/RAF1/PTPN11 (activating RASopathies).
        # No curated mechanism table was retrievable for this fix, so the
        # criterion now requires the caller to assert the mechanism and
        # otherwise refuses to apply PVS1. See CLINICAL_REVIEW.md.
        lof_is_mechanism = variant_data.get("lof_is_mechanism", None)
        if vtype in self._NULL_TYPES:
            if lof_is_mechanism is True:
                pathogenic_criteria.append("PVS1")
            elif lof_is_mechanism is False:
                not_applied.append(
                    f"PVS1 not applied: loss of function is not a known "
                    f"disease mechanism for {gene or 'this gene'}"
                )
            else:
                not_applied.append(
                    f"PVS1 not applied: whether loss of function is a known "
                    f"disease mechanism for {gene or 'this gene'} was not "
                    "asserted (Richards 2015 Table 3 requires it)"
                )
        elif vtype in self._FRAME_AMBIGUOUS_TYPES:
            not_applied.append(
                f"PVS1/PM4 not applied: variant type '{vtype}' does not state "
                "whether the reading frame is preserved; PVS1 (null) and PM4 "
                "(in-frame length change) are mutually exclusive and neither "
                "can be assigned"
            )

        # PS1: Richards 2015 Table 3 -- "Same amino acid change as a previously
        # established pathogenic variant regardless of nucleotide change". That
        # requires an amino-acid-level comparison, which only the caller can
        # supply. Presence in ClinVar is PP5 (Supporting), handled below; it was
        # previously granted PS1 (Strong), a 4x inflation of its weight.
        if variant_data.get("same_amino_acid_as_known_pathogenic", False):
            pathogenic_criteria.append("PS1")

        # PS2: De novo (confirmed)
        if is_de_novo and de_novo_confirmed:
            pathogenic_criteria.append("PS2")

        # PS3: Functional studies
        if functional in ("damaging", "loss_of_function", "deleterious"):
            pathogenic_criteria.append("PS3")

        # PM1: Mutational hot spot
        if protein_pos is not None and gene.upper() in self._HOT_SPOTS:
            for start, end in self._HOT_SPOTS[gene.upper()]:
                if start <= protein_pos <= end:
                    pathogenic_criteria.append("PM1")
                    break

        # PM2: Absent from controls. Only assignable when a frequency was
        # actually reported (see the `freq is None` branch above).
        if freq is not None and freq < 0.0001:
            pathogenic_criteria.append("PM2")

        # PM3: Trans with pathogenic (recessive)
        if variant_data.get("in_trans_with_pathogenic", False):
            pathogenic_criteria.append("PM3")

        # PM4: Richards 2015 Table 3 -- "Protein length changes due to
        # **in-frame** deletions/insertions in a non-repeat region or stop-loss
        # variants". In-frame only, so disjoint from the PVS1 null types above.
        if vtype in self._INFRAME_TYPES and not in_repeat:
            pathogenic_criteria.append("PM4")

        # PM5: Novel missense at same residue as established pathogenic
        if (
            vtype == "missense"
            and variant_data.get("same_residue_pathogenic", False)
            and "PS1" not in pathogenic_criteria
        ):
            pathogenic_criteria.append("PM5")

        # PM6: Assumed de novo
        if is_de_novo and not de_novo_confirmed:
            pathogenic_criteria.append("PM6")

        # PP1: Cosegregation
        if coseg:
            pathogenic_criteria.append("PP1")

        # PP3: Computational prediction
        if comp_pred in ("deleterious", "damaging", "pathogenic", "probably_damaging"):
            pathogenic_criteria.append("PP3")

        # PP4: Phenotype specific
        if pheno_specific:
            pathogenic_criteria.append("PP4")

        # PP5: Richards 2015 Table 3 -- "Reputable source recently reports
        # variant as pathogenic but the evidence is not available to the
        # laboratory to perform an independent evaluation." This is exactly what
        # a ClinVar assertion is, at Supporting strength.
        if in_clinvar and clinvar_class in ("pathogenic", "likely_pathogenic"):
            if "PS1" not in pathogenic_criteria:
                pathogenic_criteria.append("PP5")

        # ── Benign criteria ───────────────────────────────────────────

        # BS1: Higher frequency than expected
        if freq is not None and freq > 0.01:
            benign_criteria.append("BS1")

        # BS2: Observed in healthy adults
        if variant_data.get("seen_in_healthy", False):
            benign_criteria.append("BS2")

        # BP1: Missense in gene where only truncating causes disease
        if variant_data.get("missense_non_truncating_gene", False):
            benign_criteria.append("BP1")

        # BP3: Richards 2015 Table 3 -- "In-frame deletions/insertions in a
        # repetitive region without a known function". In-frame only, as PM4.
        if vtype in self._INFRAME_TYPES and in_repeat:
            benign_criteria.append("BP3")

        # BP4: Computational no impact
        if comp_pred in ("benign", "tolerated", "neutral"):
            benign_criteria.append("BP4")

        # BP6: Reputable source benign
        if in_clinvar and clinvar_class in ("benign", "likely_benign"):
            benign_criteria.append("BP6")

        # BP7: Synonymous no splice
        if synonymous and not splice_impact:
            benign_criteria.append("BP7")

        # ── Classification decision: Richards 2015 Table 5 ────────────
        all_criteria = pathogenic_criteria + benign_criteria
        classification, rule = combine_acmg_criteria(all_criteria)

        summary = (
            f"Pathogenic criteria: {', '.join(pathogenic_criteria) if pathogenic_criteria else 'none'}. "
            f"Benign criteria: {', '.join(benign_criteria) if benign_criteria else 'none'}. "
            f"Classification: {classification.value} "
            f"(Richards 2015 Table 5: {rule})."
        )
        if not_applied:
            summary += " Not applied: " + "; ".join(not_applied) + "."

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
        disease_norm = _normalise_disease(disease_id)
        genotype_lower = (genotype or "").lower()
        matches: List[TherapyMatch] = []

        for drug in self._ORPHAN_DRUGS:
            match_type = None
            eligibility_note = ""

            # Exact (normalised) disease-name match. A bidirectional substring
            # test dispensed the wrong drug: a query for "cardiomyopathy" is a
            # substring of "Obstructive Hypertrophic Cardiomyopathy" and
            # returned mavacamten -- a cardiac myosin inhibitor -- for a
            # dilated cardiomyopathy, where reducing contractility worsens
            # heart failure. Likewise "Long QT Syndrome" returned mexiletine,
            # catalogued only for LQT3, and "PAH" (pulmonary arterial
            # hypertension) returned sapropterin for phenylketonuria.
            if _normalise_disease(drug["disease"]) == disease_norm or disease_norm in {
                _normalise_disease(a) for a in drug.get("aliases", ())
            }:
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
            elif drug.get("gene") and re.search(
                r"\b" + re.escape(drug["gene"].upper()) + r"\b", disease_id.upper()
            ):
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

            if expected_carrier is None:
                # No expected genotype can be derived (inheritance mode not
                # recognised, or sex not recorded for an X-linked pedigree), so
                # this member contributes no informative meiosis. Jarvik &
                # Browning 2016 (PMC4908147) counts informative meioses only.
                details.append(
                    f"{member_id}: no expected genotype for inheritance "
                    f"'{inheritance}' — excluded (not an informative meiosis)"
                )
                continue

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
    ) -> Optional[str]:
        """Determine expected genotype for a family member.

        Returns ``None`` when the inheritance mode is not recognised, so the
        member contributes no informative meiosis.

        X-linked modes are tested FIRST. ``"x_linked_dominant"`` contains the
        substring ``"dominant"`` and ``"x_linked_recessive"`` contains
        ``"recessive"``, so testing the autosomal modes first made the
        hemizygosity branch unreachable for both: an affected hemizygous son in
        an X-linked dominant pedigree was scored DISCORDANT against an
        autosomal ``het`` expectation, and the family was reported as evidence
        *against* pathogenicity (the BS4 direction).
        """
        inh = (inheritance or "").lower()
        if "x_linked" in inh or "x-linked" in inh or inh in ("xld", "xlr", "xl"):
            if sex == "male":
                return "hemi" if affected else "wt"
            if sex == "female":
                # An affected female in an X-linked pedigree is heterozygous
                # (XLD) or homozygous (XLR); an unaffected female may be a
                # carrier, so she is not informative for either direction.
                if affected:
                    return "het_or_hom"
                return "het_or_wt"
            # Sex not recorded: hemizygosity cannot be evaluated.
            return None
        if "dominant" in inh:
            return "het" if affected else "wt"
        if "recessive" in inh:
            return "hom" if affected else "het_or_wt"
        return None

    @staticmethod
    def _genotype_matches(
        observed: str, expected: Optional[str], inheritance: str
    ) -> bool:
        """Check if observed genotype matches expected.

        ``expected is None`` is handled by the caller, which excludes the member
        rather than counting it. Previously an unrecognised inheritance mode
        returned the literal string ``"unknown"`` and this method returned
        ``True`` for it, so "cannot assess" was counted as a concordant meiosis
        worth +0.3 LOD -- ten relatives in a mitochondrial pedigree produced
        LOD 3.0 and "PS" strong segregation evidence with no genotype ever
        compared to an expectation.
        """
        if expected is None:
            return False
        if expected == "het_or_wt":
            return observed in ("het", "wt")
        if expected == "het_or_hom":
            return observed in ("het", "hom")
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
        disease_norm = _normalise_disease(disease)

        # Find matching disease by exact (normalised) name. The previous
        # bidirectional substring test counselled families to the wrong
        # trajectory: "muscular dystrophy" returned the DUCHENNE milestone
        # table (ambulation lost at a median of 120 months) for a child with
        # Becker muscular dystrophy, who typically walks into adulthood; and
        # "atrophy" returned SMA type 1.
        history = None
        matched_name = None
        for name, data in self._NATURAL_HISTORY.items():
            if _normalise_disease(name) == disease_norm:
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
