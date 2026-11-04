"""Behavioural reference tests for ACMG/AMP variant classification.

These tests exercise ``src.decision_support`` -- the module the agent's only
general-purpose variant classifier lives in. Before this file existed, the
3,181 lines of clinical logic in ``src/decision_support.py`` and
``src/clinical_workflows.py`` had **zero** behavioural tests: three test files
named after the clinical modules test the ingest parsers instead, and
``ACMGVariantClassifier.classify`` had no test of any kind.

Every assertion here is a *reference case* -- a described variant in, the
published classification out -- taken from:

* Richards S, Aziz N, Bale S, et al. "Standards and guidelines for the
  interpretation of sequence variants." Genet Med 2015;17(5):405-424.
  PMC4544753. Table 3 (pathogenic criteria), Table 4 (benign criteria),
  Table 5 (combining rules) and the Table 5 footnote.
* Ghosh R, Harrison SM, Rehm HL, Plon SE, Biesecker LG (ClinGen Sequence
  Variant Interpretation Working Group). "Updated recommendation for the
  benign stand-alone ACMG/AMP criterion." Hum Mutat 2018;39(11):1525-1530.
  PMC6188666. Specified BA1 and Table 1 (the BA1 exception list).

Both were retrieved 2026-11-03.
"""

import pytest

from src.decision_support import (
    ACMGVariantClassifier,
    FamilySegregationAnalyzer,
    HPOToGeneMatcher,
    NaturalHistoryPredictor,
    OrphanDrugMatcher,
    combine_acmg_criteria,
)
from src.models import ACMGClassification


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 1 — Richards 2015 Table 5: every rule set, one case each
# ═══════════════════════════════════════════════════════════════════════════

#: (criteria met, expected classification). One row per Table 5 rule set,
#: plus the two footnote cases. A point system cannot reproduce this table,
#: which is why these are asserted as exact labels and not as directions.
TABLE_5_RULE_SETS = [
    # --- Pathogenic -------------------------------------------------------
    # "1 Very Strong (PVS1) AND >=1 Strong (PS1-PS4)"
    (["PVS1", "PS3"], ACMGClassification.PATHOGENIC),
    # "1 Very Strong (PVS1) AND >=2 Moderate (PM1-PM6)"
    (["PVS1", "PM2", "PM4"], ACMGClassification.PATHOGENIC),
    # "1 Very Strong (PVS1) AND 1 Moderate AND 1 Supporting"
    (["PVS1", "PM2", "PP3"], ACMGClassification.PATHOGENIC),
    # "1 Very Strong (PVS1) AND >=2 Supporting (PP1-PP5)"
    (["PVS1", "PP1", "PP3"], ACMGClassification.PATHOGENIC),
    # ">=2 Strong (PS1-PS4)"
    (["PS2", "PS3"], ACMGClassification.PATHOGENIC),
    # "1 Strong AND >=3 Moderate"
    (["PS3", "PM1", "PM2", "PM3"], ACMGClassification.PATHOGENIC),
    # "1 Strong AND 2 Moderate AND >=2 Supporting"
    (["PS3", "PM1", "PM2", "PP1", "PP3"], ACMGClassification.PATHOGENIC),
    # "1 Strong AND 1 Moderate AND >=4 Supporting"
    (["PS3", "PM2", "PP1", "PP2", "PP3", "PP4"], ACMGClassification.PATHOGENIC),
    # --- Likely pathogenic ------------------------------------------------
    # "1 Very Strong (PVS1) AND 1 Moderate"
    (["PVS1", "PM2"], ACMGClassification.LIKELY_PATHOGENIC),
    # "1 Strong AND 1-2 Moderate"
    (["PS3", "PM2"], ACMGClassification.LIKELY_PATHOGENIC),
    (["PS3", "PM1", "PM2"], ACMGClassification.LIKELY_PATHOGENIC),
    # "1 Strong AND >=2 Supporting"
    (["PS3", "PP1", "PP3"], ACMGClassification.LIKELY_PATHOGENIC),
    # ">=3 Moderate"
    (["PM1", "PM2", "PM3"], ACMGClassification.LIKELY_PATHOGENIC),
    # "2 Moderate AND >=2 Supporting"
    (["PM1", "PM2", "PP1", "PP3"], ACMGClassification.LIKELY_PATHOGENIC),
    # "1 Moderate AND >=4 Supporting"
    (["PM2", "PP1", "PP2", "PP3", "PP4"], ACMGClassification.LIKELY_PATHOGENIC),
    # --- Benign -----------------------------------------------------------
    # "1 Stand-Alone (BA1)"
    (["BA1"], ACMGClassification.BENIGN),
    # ">=2 Strong (BS1-BS4)"
    (["BS1", "BS2"], ACMGClassification.BENIGN),
    (["BS3", "BS4"], ACMGClassification.BENIGN),
    # --- Likely benign ----------------------------------------------------
    # "1 Strong (BS1-BS4) AND 1 Supporting (BP1-BP7)"
    (["BS1", "BP7"], ACMGClassification.LIKELY_BENIGN),
    # ">=2 Supporting (BP1-BP7)"
    (["BP4", "BP7"], ACMGClassification.LIKELY_BENIGN),
]


@pytest.mark.parametrize("criteria,expected", TABLE_5_RULE_SETS)
def test_table_5_rule_sets(criteria, expected):
    """Each Richards 2015 Table 5 rule set yields exactly its published label."""
    got, rule = combine_acmg_criteria(criteria)
    assert got is expected, f"{criteria} -> {got.value}, expected {expected.value} ({rule})"


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 2 — Boundary cases on each side of every rule
# ═══════════════════════════════════════════════════════════════════════════

#: The seven rows the pre-fix point system got wrong, measured by the audit,
#: plus the boundaries either side of each Table 5 threshold. The point system
#: scored PVS1=8, PS=4/3, PM=2, PP=1 with Pathogenic at >=10, which both
#: over-called (PVS1 alone -> Likely pathogenic) and under-called
#: (PS3+PM2 -> VUS, a missed PKU diagnosis).
TABLE_5_BOUNDARIES = [
    # 1 Very Strong alone satisfies NO Table 5 rule set.
    (["PVS1"], ACMGClassification.VUS),
    (["PVS1", "PM2"], ACMGClassification.LIKELY_PATHOGENIC),
    (["PVS1", "PM2", "PM4"], ACMGClassification.PATHOGENIC),
    # 1 Strong alone is not enough; 2 Strong is Pathogenic whichever two.
    (["PS3"], ACMGClassification.VUS),
    (["PS1", "PS2"], ACMGClassification.PATHOGENIC),
    (["PS3", "PS4"], ACMGClassification.PATHOGENIC),
    # 1 Strong + 2 Moderate is Likely pathogenic; + 3 Moderate is Pathogenic.
    (["PS3", "PM1", "PM2"], ACMGClassification.LIKELY_PATHOGENIC),
    (["PS3", "PM1", "PM2", "PM3"], ACMGClassification.PATHOGENIC),
    # 2 Moderate alone is not enough; 3 Moderate is Likely pathogenic.
    (["PM1", "PM2"], ACMGClassification.VUS),
    (["PM1", "PM2", "PM3"], ACMGClassification.LIKELY_PATHOGENIC),
    # 1 Supporting benign alone is a VUS, not Likely benign.
    (["BP7"], ACMGClassification.VUS),
    (["BP4", "BP7"], ACMGClassification.LIKELY_BENIGN),
    # 1 Strong benign alone is a VUS, not Likely benign or Benign.
    (["BS1"], ACMGClassification.VUS),
    (["BS1", "BP7"], ACMGClassification.LIKELY_BENIGN),
    (["BS1", "BS2"], ACMGClassification.BENIGN),
    # 1 Supporting pathogenic alone is a VUS.
    (["PP3"], ACMGClassification.VUS),
    # Nothing met is a VUS, per the Table 5 footnote.
    ([], ACMGClassification.VUS),
]


@pytest.mark.parametrize("criteria,expected", TABLE_5_BOUNDARIES)
def test_table_5_boundaries(criteria, expected):
    """Either side of every Table 5 threshold lands on the published label."""
    got, rule = combine_acmg_criteria(criteria)
    assert got is expected, f"{criteria} -> {got.value}, expected {expected.value} ({rule})"


def test_pvs1_alone_is_not_likely_pathogenic():
    """PVS1 alone satisfies no Table 5 rule set, so it is a VUS.

    The point system scored PVS1 at 8 against a Likely-pathogenic threshold of
    6 and returned Likely pathogenic.
    """
    got, _ = combine_acmg_criteria(["PVS1"])
    assert got is ACMGClassification.VUS


def test_one_strong_plus_one_moderate_is_likely_pathogenic_not_vus():
    """PS3+PM2 is "1 Strong AND 1-2 Moderate" -> Likely pathogenic.

    The point system scored 3+2=5 against a Likely-pathogenic threshold of 6
    and returned VUS. Clinical consequence: a *PAH* missense with a published
    enzyme assay at 4% residual activity in a newborn with raised
    phenylalanine was a VUS, and a VUS does not trigger dietary phenylalanine
    restriction -- untreated PKU causes irreversible intellectual disability
    within weeks.
    """
    got, _ = combine_acmg_criteria(["PS3", "PM2"])
    assert got is ACMGClassification.LIKELY_PATHOGENIC


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 3 — The contradiction rule (Table 5 footnote)
# ═══════════════════════════════════════════════════════════════════════════


def test_contradictory_evidence_returns_vus():
    """Contradictory benign and pathogenic criteria -> Uncertain significance.

    Richards 2015 Table 5 footnote: "Variants should be classified as Uncertain
    Significance if other criteria are unmet or the criteria for benign and
    pathogenic are contradictory."

    The point system tested ``benign_score >= 6`` first and unconditionally, so
    three Strong pathogenic criteria including a functional enzyme assay were
    overwritten by founder-population frequency and the call came back
    ``benign`` -- enzyme replacement therapy for Gaucher disease not offered.
    """
    got, rule = combine_acmg_criteria(["PS1", "PS2", "PS3", "PP3", "BS1", "BS2"])
    assert got is ACMGClassification.VUS
    assert "contradictory" in rule


def test_contradiction_rule_reports_both_sides():
    """The returned rule names the pathogenic AND the benign rule that fired."""
    _, rule = combine_acmg_criteria(["PVS1", "PS1", "PM2", "BS1", "BS2"])
    assert "pathogenic rule" in rule and "benign rule" in rule


def test_pathogenic_does_not_silently_carry_benign_criteria():
    """A Pathogenic-satisfying set plus one Strong benign criterion is a VUS.

    Previously ``PVS1+PS1+PM2`` vs ``BS2`` returned ``pathogenic`` with ``BS2``
    sitting unremarked in the returned criteria list.
    """
    got, _ = combine_acmg_criteria(["PVS1", "PS1", "PM2", "BS2", "BS1"])
    assert got is ACMGClassification.VUS


def test_unknown_criterion_codes_are_ignored_not_counted():
    """A code outside Tables 3/4 contributes to no bucket."""
    got, _ = combine_acmg_criteria(["PVS1", "PM2", "NOT_A_CRITERION", "PX9"])
    assert got is ACMGClassification.LIKELY_PATHOGENIC


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 4 — Absent-field invariance: no criterion from a missing field
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def classifier():
    return ACMGVariantClassifier()


def test_empty_variant_data_yields_vus_with_no_criteria(classifier):
    """``classify({})`` returns VUS with an EMPTY criteria list.

    A missing ``population_frequency`` used to be read as ``0.0``, i.e. "absent
    from controls", so PM2 -- the criterion whose entire purpose is to prove the
    variant was *looked for* in population data -- was awarded for not looking.
    ``classify({})`` returned ``vus ['PM2']``.
    """
    cls, criteria, summary = classifier.classify({})
    assert cls is ACMGClassification.VUS
    assert criteria == []
    assert "population frequency not available" in summary


def test_missing_frequency_grants_no_pm2(classifier):
    """No frequency supplied -> no PM2, and the report says why."""
    cls, criteria, summary = classifier.classify({"gene": "ZZZZZ", "variant_type": "missense"})
    assert "PM2" not in criteria
    assert cls is ACMGClassification.VUS
    assert "population frequency not available" in summary


def test_reported_zero_frequency_does_grant_pm2(classifier):
    """An explicitly reported frequency of 0.0 is real evidence and grants PM2.

    "Not supplied" and "gnomAD was queried and the variant is absent" must be
    distinguishable, and the second one still earns the criterion.
    """
    _, criteria, _ = classifier.classify({
        "gene": "PAH", "variant_type": "missense", "population_frequency": 0.0,
    })
    assert "PM2" in criteria


def test_missing_frequency_grants_no_bs1_or_ba1(classifier):
    """No frequency supplied -> neither benign frequency criterion."""
    _, criteria, _ = classifier.classify({"gene": "GBA1", "variant_type": "missense"})
    assert "BS1" not in criteria
    assert "BA1" not in criteria


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 5 — Criterion disjointness
# ═══════════════════════════════════════════════════════════════════════════


def test_pvs1_and_pm4_are_never_both_granted(classifier):
    """A null variant is not an in-frame length change; they are exclusive.

    Richards 2015 Table 3: PVS1 covers "single or multi-exon deletion" as a
    **null** variant; PM4 covers "protein length changes due to **in-frame**
    deletions/insertions ... or stop-loss variants". "deletion" and "insertion"
    appeared in both tests, so one molecular event scored both criteria --
    ``{"gene": "FBN1", "variant_type": "deletion"}``, two dictionary keys,
    returned ``pathogenic ['PVS1','PM2','PM4']``.
    """
    for vtype in ("deletion", "insertion", "indel", "whole_exon_deletion",
                  "frameshift", "inframe_deletion", "inframe_insertion",
                  "stop_loss", "nonsense", "splice_site"):
        _, criteria, _ = classifier.classify({
            "gene": "FBN1", "variant_type": vtype,
            "lof_is_mechanism": True, "population_frequency": 0.0,
        })
        assert not ("PVS1" in criteria and "PM4" in criteria), vtype


def test_two_key_input_is_not_pathogenic(classifier):
    """A gene symbol and a variant-type string are not a Pathogenic call.

    ``classify({"gene": "FBN1", "variant_type": "deletion"})`` -- no position,
    no frequency, no phenotype, no mechanism assertion -- returned
    ``pathogenic`` and fired a cross-agent ``critical_alert``.
    """
    cls, criteria, _ = classifier.classify({"gene": "FBN1", "variant_type": "deletion"})
    assert cls is ACMGClassification.VUS
    assert criteria == []


def test_ambiguous_frame_grants_neither_pvs1_nor_pm4(classifier):
    """A variant type that does not state the frame earns neither criterion."""
    _, criteria, summary = classifier.classify({
        "gene": "COL1A1", "variant_type": "deletion",
        "lof_is_mechanism": True, "population_frequency": 0.0,
    })
    assert "PVS1" not in criteria and "PM4" not in criteria
    assert "reading frame is preserved" in summary


def test_pm4_and_bp3_are_never_both_granted(classifier):
    """PM4 is non-repeat, BP3 is repeat; they cannot both apply."""
    for in_repeat in (True, False):
        _, criteria, _ = classifier.classify({
            "gene": "X", "variant_type": "inframe_deletion",
            "in_repeat_region": in_repeat, "population_frequency": 0.0,
        })
        assert not ("PM4" in criteria and "BP3" in criteria)


def test_ps1_and_pp5_are_never_both_granted(classifier):
    """PS1 suppresses PP5, so a ClinVar assertion is not counted twice."""
    _, criteria, _ = classifier.classify({
        "gene": "X", "variant_type": "missense",
        "same_amino_acid_as_known_pathogenic": True,
        "in_clinvar": True, "clinvar_classification": "pathogenic",
        "population_frequency": 0.0,
    })
    assert "PS1" in criteria
    assert "PP5" not in criteria


def test_pp3_and_bp4_are_never_both_granted(classifier):
    """Computational evidence points one way only."""
    for pred in ("deleterious", "benign", "tolerated", "damaging"):
        _, criteria, _ = classifier.classify({
            "gene": "X", "variant_type": "missense",
            "computational_prediction": pred, "population_frequency": 0.0,
        })
        assert not ("PP3" in criteria and "BP4" in criteria), pred


def test_ps2_and_pm6_are_never_both_granted(classifier):
    """Confirmed de novo is PS2; assumed de novo is PM6. Not both."""
    for confirmed in (True, False):
        _, criteria, _ = classifier.classify({
            "gene": "X", "variant_type": "missense", "is_de_novo": True,
            "de_novo_confirmed": confirmed, "population_frequency": 0.0,
        })
        assert not ("PS2" in criteria and "PM6" in criteria)


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 6 — PVS1 is gated on disease MECHANISM, not gene constraint
# ═══════════════════════════════════════════════════════════════════════════

#: Genes whose disease mechanism is explicitly NOT loss of function, and which
#: the pre-fix hand-written "pLI > 0.9" list nonetheless granted PVS1 to:
#: FGFR3 (achondroplasia is gain-of-function), HTT (CAG expansion), BRAF /
#: RAF1 / PTPN11 (activating RASopathies), MYH7 (dominant-negative missense).
NON_LOF_MECHANISM_GENES = ["FGFR3", "HTT", "BRAF", "RAF1", "PTPN11", "MYH7"]


@pytest.mark.parametrize("gene", NON_LOF_MECHANISM_GENES)
def test_truncating_variant_in_non_lof_gene_is_not_pathogenic(classifier, gene):
    """A truncating variant in a gene where LOF is not the mechanism.

    Richards 2015 Table 3, PVS1: "Null variant ... **in a gene where loss of
    function (LOF) is a known mechanism of disease**".

    ``{"gene": "MYH7", "variant_type": "nonsense"}`` returned ``pathogenic
    ['PVS1','PM2']`` and fired a cross-agent critical alert, entering a family
    into cascade predictive testing on a variant class that is not an
    established mechanism for that gene. The *HTT* case sent at-risk relatives
    toward predictive testing for Huntington disease.
    """
    cls, criteria, _ = classifier.classify({
        "gene": gene, "variant_type": "nonsense", "population_frequency": 0.0,
    })
    assert "PVS1" not in criteria
    assert cls is not ACMGClassification.PATHOGENIC


def test_pvs1_requires_an_explicit_mechanism_assertion(classifier):
    """Without a mechanism assertion PVS1 is refused, and the report says so."""
    cls, criteria, summary = classifier.classify({
        "gene": "SCN1A", "variant_type": "nonsense", "population_frequency": 0.0,
    })
    assert "PVS1" not in criteria
    assert "was not asserted" in summary
    assert cls is not ACMGClassification.PATHOGENIC


def test_pvs1_is_applied_when_mechanism_is_asserted(classifier):
    """With LOF asserted, a null variant absent from controls is Likely path.

    Richards Table 5: 1 Very Strong + 1 Moderate = Likely pathogenic. Not
    Pathogenic -- the point system returned Pathogenic here, and "Pathogenic"
    closes a differential that a Likely-pathogenic call keeps open.
    """
    cls, criteria, _ = classifier.classify({
        "gene": "COL1A1", "variant_type": "frameshift",
        "lof_is_mechanism": True, "population_frequency": 0.0,
    })
    assert criteria == ["PVS1", "PM2"]
    assert cls is ACMGClassification.LIKELY_PATHOGENIC


def test_mechanism_explicitly_false_is_reported_differently(classifier):
    """"LOF is not the mechanism" and "nobody said" are distinguishable."""
    _, _, summary = classifier.classify({
        "gene": "FGFR3", "variant_type": "nonsense",
        "lof_is_mechanism": False, "population_frequency": 0.0,
    })
    assert "is not a known" in summary


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 7 — PS1 is PP5: a ClinVar assertion is Supporting, not Strong
# ═══════════════════════════════════════════════════════════════════════════


def test_clinvar_presence_alone_is_pp5_and_yields_vus(classifier):
    """A ClinVar Pathogenic assertion plus nothing else is a VUS.

    Richards 2015 Table 3: PS1 is "same amino acid change as a previously
    established pathogenic variant"; PP5 is "reputable source recently reports
    variant as pathogenic but the evidence is not available to the laboratory
    to perform an independent evaluation." The code granted PS1 -- Strong, 4x
    the weight -- for the latter, so a single-submitter ClinVar entry plus a
    free PM2 was laundered into a Likely pathogenic report.
    """
    cls, criteria, _ = classifier.classify({
        "in_clinvar": True, "clinvar_classification": "pathogenic",
    })
    assert criteria == ["PP5"]
    assert cls is ACMGClassification.VUS


def test_clinvar_plus_absent_from_controls_is_still_vus(classifier):
    """PP5 + PM2 is 1 Moderate + 1 Supporting: no Table 5 rule set."""
    cls, criteria, _ = classifier.classify({
        "in_clinvar": True, "clinvar_classification": "pathogenic",
        "population_frequency": 0.0,
    })
    assert sorted(criteria) == ["PM2", "PP5"]
    assert cls is ACMGClassification.VUS


def test_ps1_requires_amino_acid_level_evidence(classifier):
    """PS1 is granted only on an explicit amino-acid-level assertion."""
    _, criteria, _ = classifier.classify({
        "gene": "X", "variant_type": "missense",
        "same_amino_acid_as_known_pathogenic": True, "population_frequency": 0.0,
    })
    assert "PS1" in criteria


def test_likely_pathogenic_clinvar_is_not_previously_established(classifier):
    """A ClinVar "likely pathogenic" assertion is PP5, never PS1."""
    _, criteria, _ = classifier.classify({
        "in_clinvar": True, "clinvar_classification": "likely_pathogenic",
    })
    assert "PS1" not in criteria
    assert "PP5" in criteria


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 8 — The ClinGen SVI BA1 specification and exception list
# ═══════════════════════════════════════════════════════════════════════════

#: Ghosh et al. 2018 (PMC6188666) Table 1 -- the nine variants to which BA1
#: must not be applied. (gene, HGVS c.)
BA1_EXCEPTION_VARIANTS = [
    ("ACAD9", "c.-44_-41dupTAAG"),
    ("GJB2", "c.109G>A"),
    ("HFE", "c.187C>G"),
    ("HFE", "c.845G>A"),
    ("MEFV", "c.1105C>T"),
    ("MEFV", "c.1223G>A"),
    ("PIBF1", "c.1214G>A"),
    ("ACADS", "c.511C>T"),
    ("BTD", "c.1330G>C"),
]


@pytest.mark.parametrize("gene,hgvs_c", BA1_EXCEPTION_VARIANTS)
def test_ba1_exception_variants_do_not_get_ba1(classifier, gene, hgvs_c):
    """None of the nine ClinGen SVI Table 1 variants receives BA1.

    Specified BA1 (Ghosh et al. 2018): "Allele frequency is >0.05 in any
    general continental population dataset of at least 2,000 observed alleles
    **and found in a gene without a gene- or variant-specific BA1
    modification**." Table 1 is that modification list.
    """
    cls, criteria, summary = classifier.classify({
        "gene": gene, "hgvs_c": hgvs_c,
        "population_frequency": 0.08, "allele_count_observed": 120000,
    })
    assert "BA1" not in criteria
    assert cls is not ACMGClassification.BENIGN
    assert "exception list" in summary


def test_hfe_c282y_homozygote_is_not_benign(classifier):
    """A *HFE* C282Y homozygote with a ClinVar Pathogenic assertion.

    Clinical reference case: a 48-year-old man with fatigue, arthropathy,
    raised ferritin and transferrin saturation 72%, homozygous for *HFE*
    c.845G>A p.(Cys282Tyr) -- allele frequency ~6% in non-Finnish Europeans.
    The classifier returned **Benign** with the summary "Population frequency
    0.0637 > 5% -- classified as Benign" from an unconditional early return
    that never read the pathogenic evidence, so hereditary haemochromatosis
    was excluded on the report, therapeutic phlebotomy was not started, and
    iron overload progresses to cirrhosis and hepatocellular carcinoma.
    """
    cls, criteria, _ = classifier.classify({
        "gene": "HFE", "hgvs_c": "c.845G>A", "hgvs_p": "p.Cys282Tyr",
        "variant_type": "missense", "zygosity": "homozygous",
        "population_frequency": 0.0637, "allele_count_observed": 120000,
        "in_clinvar": True, "clinvar_classification": "pathogenic",
    })
    assert cls is not ACMGClassification.BENIGN
    assert "BA1" not in criteria


def test_ba1_exception_matches_on_protein_notation(classifier):
    """The exception list is reachable from p. notation as well as c."""
    _, criteria, _ = classifier.classify({
        "gene": "GJB2", "hgvs_p": "p.Val37Ile",
        "population_frequency": 0.09, "allele_count_observed": 50000,
    })
    assert "BA1" not in criteria


def test_ba1_requires_the_observed_allele_count(classifier):
    """Without an allele count the ClinGen >=2,000 threshold is uncheckable."""
    cls, criteria, summary = classifier.classify({
        "gene": "ZZZZZ", "population_frequency": 0.10,
    })
    assert "BA1" not in criteria
    assert cls is not ACMGClassification.BENIGN
    assert "observed alleles" in summary


def test_ba1_refused_below_2000_observed_alleles(classifier):
    """A 6% frequency from a 400-allele dataset does not reach BA1."""
    _, criteria, summary = classifier.classify({
        "gene": "ZZZZZ", "population_frequency": 0.06, "allele_count_observed": 400,
    })
    assert "BA1" not in criteria
    assert "below the ClinGen SVI minimum" in summary


def test_ba1_is_applied_when_the_specification_is_satisfied(classifier):
    """A genuinely common variant in an unmodified gene is Benign."""
    cls, criteria, _ = classifier.classify({
        "gene": "ZZZZZ", "population_frequency": 0.12, "allele_count_observed": 100000,
    })
    assert "BA1" in criteria
    assert cls is ACMGClassification.BENIGN


def test_ba1_does_not_short_circuit_pathogenic_evidence(classifier):
    """BA1 plus strong pathogenic evidence is a contradiction, so a VUS.

    The pre-fix code returned Benign from line 357 before any pathogenic
    criterion was read.
    """
    cls, criteria, _ = classifier.classify({
        "gene": "ZZZZZ", "variant_type": "missense",
        "population_frequency": 0.06, "allele_count_observed": 100000,
        "functional_studies": "damaging",
        "same_amino_acid_as_known_pathogenic": True,
    })
    assert "BA1" in criteria
    assert cls is ACMGClassification.VUS


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 9 — Family segregation analysis
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def segregation():
    return FamilySegregationAnalyzer()


def _ten_member_pedigree():
    return [{"id": f"m{i}", "affected": True, "sex": "female"} for i in range(10)]


@pytest.mark.parametrize("inheritance", ["mitochondrial", "", "unknown", "digenic"])
def test_unrecognised_inheritance_yields_no_segregation_evidence(segregation, inheritance):
    """An inheritance mode no rule understands produces NO evidence.

    ``_expected_genotype`` returned the literal string ``"unknown"`` for any
    unrecognised mode and ``_genotype_matches`` returned ``True`` for it, so
    "cannot assess" was counted as a concordant meiosis worth +0.3 LOD. Ten
    relatives in a mitochondrial pedigree measured ``lod_score 3.0,
    evidence_level "PS", concordant 10, discordant 0`` -- Pathogenic Strong
    cosegregation evidence manufactured from a pedigree in which no genotype
    was ever compared to an expectation. PS-strength segregation is exactly
    what tips a VUS to Likely pathogenic in a family.

    Jarvik GP & Browning BL, Am J Hum Genet 2016;98:1077-1081 (PMC4908147),
    Table 1: the LOD cutoffs are computed over *informative meioses*, and a
    relative whose expected genotype is unknown contributes none.
    """
    pedigree = _ten_member_pedigree()
    genotypes = {m["id"]: "het" for m in pedigree}
    result = segregation.analyze(pedigree, genotypes, inheritance)
    assert result["concordant"] == 0
    assert result["evidence_level"] == "insufficient"
    assert result["lod_score"] == 0.0


def test_x_linked_dominant_family_segregates(segregation):
    """An X-linked dominant pedigree is not reported as non-segregating.

    Clinical reference case: a three-generation X-linked hypophosphataemic
    rickets family (*PHEX*, XLD) -- affected mother ``het``, affected son
    ``hemi``, unaffected son ``wt``.

    ``"x_linked_dominant"`` contains the substring ``"dominant"``, so the
    autosomal branch consumed it and the hemizygosity branch was unreachable.
    The affected hemizygous son was scored DISCORDANT against an autosomal
    ``het`` expectation: ``lod_score -0.4``, ``evidence_level
    "against_segregation"``, ``classification "non_segregating"`` -- a
    correctly segregating X-linked variant reported as evidence *against*
    pathogenicity (the BS4 direction), so the family's diagnosis is lost.
    """
    pedigree = [
        {"id": "mother", "affected": True, "sex": "female"},
        {"id": "son1", "affected": True, "sex": "male"},
        {"id": "son2", "affected": False, "sex": "male"},
    ]
    genotypes = {"mother": "het", "son1": "hemi", "son2": "wt"}
    result = segregation.analyze(pedigree, genotypes, "x_linked_dominant")
    assert result["discordant"] == 0
    assert result["classification"] != "non_segregating"
    assert result["evidence_level"] != "against_segregation"


def test_x_linked_recessive_affected_male_is_concordant(segregation):
    """An affected male in an XLR pedigree is hemizygous, not homozygous."""
    pedigree = [
        {"id": "son", "affected": True, "sex": "male"},
        {"id": "mother", "affected": False, "sex": "female"},
    ]
    genotypes = {"son": "hemi", "mother": "het"}
    result = segregation.analyze(pedigree, genotypes, "x_linked_recessive")
    assert result["discordant"] == 0


def test_x_linked_without_recorded_sex_is_not_scored(segregation):
    """Hemizygosity cannot be evaluated without sex, so the member is excluded."""
    pedigree = [{"id": "p1", "affected": True, "sex": "unknown"}]
    result = segregation.analyze(pedigree, {"p1": "het"}, "x_linked_dominant")
    assert result["evidence_level"] == "insufficient"


def test_autosomal_dominant_still_works(segregation):
    """The autosomal branches are unchanged by the X-linked reordering."""
    pedigree = [
        {"id": "a", "affected": True, "sex": "male"},
        {"id": "b", "affected": False, "sex": "female"},
    ]
    result = segregation.analyze(
        pedigree, {"a": "het", "b": "wt"}, "autosomal_dominant"
    )
    assert result["concordant"] == 2
    assert result["discordant"] == 0


def test_autosomal_recessive_still_works(segregation):
    """Affected homozygote, unaffected heterozygous carrier: both concordant."""
    pedigree = [
        {"id": "a", "affected": True, "sex": "male"},
        {"id": "b", "affected": False, "sex": "female"},
    ]
    result = segregation.analyze(
        pedigree, {"a": "hom", "b": "het"}, "autosomal_recessive"
    )
    assert result["concordant"] == 2


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 10 — HPO information content
# ═══════════════════════════════════════════════════════════════════════════


def test_unknown_hpo_term_is_not_the_most_informative_term():
    """A term the system knows nothing about scores no IC at all.

    ``_TERM_FREQUENCIES.get(hpo_id, 0.001)`` gave any unknown string an IC of
    9.97 -- nearly 4x the IC of HP:0001250 Seizures (2.74) -- so the matcher's
    highest-weighted evidence was systematically its least-known evidence, and
    gene ranking was driven by which terms the author forgot to add to a
    43-entry dict. An unannotated term has undefined IC, not maximal IC.
    """
    matcher = HPOToGeneMatcher()
    assert matcher._compute_ic("banana") is None
    assert matcher._compute_ic("HP:9999999") is None
    assert matcher._compute_ic("HP:0006538") is None

    known = matcher._compute_ic("HP:0001250")
    assert known is not None and known > 0


def test_unknown_term_contributes_zero_to_similarity():
    """An unscoreable term adds nothing rather than dominating the score."""
    matcher = HPOToGeneMatcher()
    assert matcher._max_ic_match("HP:0006538", ["HP:0006538"]) == 0.0


def test_unscored_terms_are_reported():
    """Terms matched but not scoreable are named in the result."""
    matcher = HPOToGeneMatcher()
    results = matcher.process_patient(["HP:0001250"])
    assert all("unscored_terms" in r for r in results)


# ═══════════════════════════════════════════════════════════════════════════
# GROUP 11 — Disease-name matching must be exact, not substring
# ═══════════════════════════════════════════════════════════════════════════


def test_cardiomyopathy_query_does_not_dispense_mavacamten():
    """"cardiomyopathy" must not return a drug licensed for obstructive HCM.

    Clinical reference case: a 29-year-old with dilated cardiomyopathy and
    LVEF 28%. A bidirectional substring test ("cardiomyopathy" is a substring
    of "Obstructive Hypertrophic Cardiomyopathy") returned **mavacamten** with
    ``match_type "exact_disease"`` and ``status APPROVED_FDA``. Mavacamten is a
    cardiac myosin inhibitor; reducing contractility in a dilated
    cardiomyopathy at EF 28% worsens heart failure.
    """
    matcher = OrphanDrugMatcher()
    names = {m.therapy_name.lower() for m in matcher.match("cardiomyopathy")}
    assert not any("mavacamten" in n or "camzyos" in n for n in names)


def test_long_qt_query_does_not_dispense_lqt3_only_drug():
    """"Long QT Syndrome" must not return a drug catalogued for LQT3 only."""
    matcher = OrphanDrugMatcher()
    names = {m.therapy_name.lower() for m in matcher.match("Long QT Syndrome")}
    assert not any("mexiletine" in n for n in names)


def test_exact_disease_name_still_matches():
    """An exact catalogue name still returns its drug."""
    matcher = OrphanDrugMatcher()
    target = matcher._ORPHAN_DRUGS[0]
    names = {m.therapy_name for m in matcher.match(target["disease"])}
    assert target["drug"] in names


def test_muscular_dystrophy_query_does_not_return_duchenne_history():
    """"muscular dystrophy" must not return the Duchenne milestone table.

    Clinical reference case: a 7-year-old with **Becker** muscular dystrophy.
    The substring match returned the Duchenne natural history -- loss of
    independent ambulation at a median of 120 months, ventilatory support at
    216 months -- at ``confidence 0.80``, counselling the family to a Duchenne
    trajectory for a child who will typically walk into adulthood.
    """
    predictor = NaturalHistoryPredictor()
    result = predictor.predict("muscular dystrophy")
    assert result["milestones"] == []
    assert "No natural history data" in result.get("message", "")


def test_atrophy_query_does_not_return_sma_history():
    """A bare "atrophy" query must not return the SMA type 1 trajectory."""
    predictor = NaturalHistoryPredictor()
    result = predictor.predict("atrophy")
    assert result["milestones"] == []


def test_exact_natural_history_name_still_matches():
    """An exact catalogue name still returns its milestone table."""
    predictor = NaturalHistoryPredictor()
    name = next(iter(predictor._NATURAL_HISTORY))
    result = predictor.predict(name)
    assert result["milestones"]
