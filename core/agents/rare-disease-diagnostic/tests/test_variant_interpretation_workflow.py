"""Behavioural reference tests for the clinical workflows in src/clinical_workflows.py.

Before this file existed, ``src/clinical_workflows.py`` (2,118 lines, ten
clinical workflows, the Ghent criteria, a second ACMG classifier and the
metabolic emergency protocols) had no behavioural test: the file named after it
tests the ``DiagnosticWorkflowType`` enum.

Reference cases are drawn from:

* Richards S, et al. Genet Med 2015;17:405-424 (PMC4544753), Tables 3 and 5.
* Loeys BL, Dietz HC, Braverman AC, et al. "The revised Ghent nosology for the
  Marfan syndrome." J Med Genet 2010;47(7):476-485 (PMID 20591885) -- for the
  four no-family-history criteria the workflow implements. The paper itself is
  paywalled and was NOT retrieved for this branch, so the family-history
  criteria are deliberately unimplemented and the tests below assert that the
  workflow SAYS so rather than reporting "criteria not met".
"""

import pytest

from src.clinical_workflows import (
    ConnectiveTissueWorkflow,
    WESWGSInterpretationWorkflow,
)
from src.models import ACMGClassification


# ═══════════════════════════════════════════════════════════════════════════
# WES/WGS interpretation — there is now ONE ACMG classifier in this agent
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def wes():
    return WESWGSInterpretationWorkflow()


def _classification_of(result, variant_id="v1"):
    for vc in result.diagnostic_result.variants:
        if vc.variant_id == variant_id:
            return vc
    raise AssertionError(f"variant {variant_id} not in result")


def test_nonsense_variant_with_no_gene_is_not_pathogenic(wes):
    """``{"variant_type": "nonsense"}`` must not yield a Pathogenic call.

    This workflow carried its own classifier which appended PVS1 on variant
    type alone -- no gene consulted, not a mechanism list, not a constraint
    metric -- and combined it with a PM2 granted from a missing frequency field
    to reach Pathogenic.

    Clinical reference case: a 6-month-old with hypotonia and failure to
    thrive has exome sequencing and the pipeline hands this workflow a
    truncating *TTN* variant -- a variant class carried by roughly 1 in 500
    healthy adults. The agent returned **Pathogenic (PVS1, PM2)** at severity
    HIGH, so a benign truncating variant in a tolerant gene became the reported
    molecular diagnosis, the real cause was never pursued, and the parents were
    told they are carriers.

    Richards 2015 Table 3, PVS1: "Null variant ... **in a gene where loss of
    function (LOF) is a known mechanism of disease**."
    """
    result = wes.execute({"variants": [
        {"variant_id": "v1", "variant_type": "nonsense", "gene": "TTN"},
    ]})
    vc = _classification_of(result)
    assert vc.classification is not ACMGClassification.PATHOGENIC
    assert "PVS1" not in vc.acmg_criteria
    assert "PM2" not in vc.acmg_criteria


def test_variant_with_no_criteria_met_is_vus_not_likely_benign(wes):
    """Score 0 means "no evidence either way", which is a VUS.

    ``elif score <= 0: classification = LIKELY_BENIGN`` reported a variant with
    NO criteria met as a benign call, and the ranking step then sorted
    Likely-benign variants to the bottom where a top-N view drops them.

    Clinical reference case: a 3-year-old with global developmental delay and a
    novel missense variant in a gene with no ClinVar entry, no supplied
    in-silico prediction, and a phenotype the inline 21-disease table does not
    cover. The single most important output of an exome -- the un-annotated
    novel variant that needs functional follow-up and 12-month reanalysis --
    was filed as benign and hidden, and lost its reanalysis flag too, because
    the recommendation block offers reanalysis only for a VUS.

    Richards 2015 Table 5 footnote: unmet criteria -> Uncertain significance.
    Likely benign requires 1 Strong benign + 1 Supporting benign, or >=2
    Supporting benign -- and this workflow computes no benign criteria at all,
    so it could never legitimately reach Likely benign.
    """
    result = wes.execute({"variants": [
        {"variant_id": "v1", "gene": "ZZZZZ", "variant_type": "missense",
         "population_frequency": 0.001},
    ]})
    vc = _classification_of(result)
    assert vc.acmg_criteria == []
    assert vc.classification is ACMGClassification.VUS


def test_novel_variant_keeps_its_reanalysis_recommendation(wes):
    """A VUS earns the reanalysis recommendation that Likely benign loses."""
    result = wes.execute({"variants": [
        {"variant_id": "v1", "gene": "ZZZZZ", "variant_type": "missense",
         "population_frequency": 0.001},
    ]})
    assert any("reanalysis" in r for r in result.recommendations)


def test_missing_frequency_is_reported_not_treated_as_zero(wes):
    """A variant with no frequency is flagged, not silently filtered as 0.0."""
    result = wes.execute({"variants": [
        {"variant_id": "v1", "gene": "ZZZZZ", "variant_type": "missense"},
    ]})
    assert any("no population frequency" in f for f in result.findings)


def test_clinvar_presence_does_not_reach_likely_pathogenic(wes):
    """A ClinVar assertion is PP5 (Supporting), so it cannot carry a call alone.

    The workflow's own classifier granted PS1 -- Strong -- here.
    """
    result = wes.execute({"variants": [
        {"variant_id": "v1", "gene": "BRCA2", "variant_type": "missense",
         "in_clinvar": True, "clinvar_classification": "pathogenic",
         "population_frequency": 0.0},
    ]})
    vc = _classification_of(result)
    assert "PS1" not in vc.acmg_criteria
    assert vc.classification is ACMGClassification.VUS


def test_truncating_variant_with_lof_mechanism_is_likely_pathogenic(wes):
    """The legitimate case still classifies: PVS1 + PM2 = Likely pathogenic.

    Richards Table 5: "1 Very Strong (PVS1) AND 1 Moderate" = Likely
    pathogenic, NOT Pathogenic. A novel frameshift in *COL1A1* in a child with
    three long-bone fractures, absent from gnomAD, is Likely pathogenic --
    which keeps non-accidental injury and the other OI genes in play and
    mandates confirmatory work. The point system returned Pathogenic, which
    closes that differential.
    """
    result = wes.execute({"variants": [
        {"variant_id": "v1", "gene": "COL1A1", "variant_type": "frameshift",
         "lof_is_mechanism": True, "population_frequency": 0.0},
    ]})
    vc = _classification_of(result)
    assert vc.classification is ACMGClassification.LIKELY_PATHOGENIC
    assert sorted(vc.acmg_criteria) == ["PM2", "PVS1"]


def test_references_name_the_table_implemented(wes):
    """A guideline citation names the rule set the code actually implements.

    ``guideline_references`` previously carried a bare "ClinGen Sequence
    Variant Interpretation (SVI) Working Group recommendations" -- naming no
    recommendation -- beside a classifier that followed neither Table 5 nor any
    SVI recommendation. The citation is the artifact that makes a call
    reviewable; here it certified a rule set the code did not implement.
    """
    result = wes.execute({"variants": [
        {"variant_id": "v1", "gene": "X", "variant_type": "missense",
         "population_frequency": 0.0},
    ]})
    refs = " ".join(result.guideline_references)
    assert "Table 5" in refs
    assert "PMC4544753" in refs


# ═══════════════════════════════════════════════════════════════════════════
# Connective tissue — revised Ghent nosology
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def connective():
    return ConnectiveTissueWorkflow()


def test_unmeasured_aortic_root_is_not_asserted_as_normal(connective):
    """An unmeasured aortic root must never be reported as "Z=0".

    ``aortic_z = ghent.get("aortic_root_z", 0)`` made the findings string
    assert **"aortic Z=0"** -- a normal aortic root -- for a patient whose
    aorta had never been measured.
    """
    result = connective.execute({"ghent_criteria": {"ectopia_lentis": True}})
    findings = " ".join(result.findings)
    assert "Z=0" not in findings
    assert "not measured" in findings or "not provided" in findings


def test_unmeasured_aortic_root_prompts_measurement(connective):
    """The workflow asks for the measurement it needs."""
    result = connective.execute({"ghent_criteria": {"ectopia_lentis": True}})
    assert any("aortic root" in r.lower() for r in result.recommendations)


def test_family_history_is_not_silently_discarded(connective):
    """Family history of Marfan syndrome reaches the output.

    Clinical reference case: a 16-year-old with bilateral ectopia lentis whose
    father has genetically confirmed Marfan syndrome, aortic root not yet
    dilated. ``ghent.get("family_history", False)`` was a bare expression
    statement -- the value was computed and bound to nothing -- so every
    family-history criterion of the revised Ghent nosology was unimplemented
    and the code never said so. Measured output was "Ghent criteria NOT fully
    met (aortic Z=0, ...)" with the single recommendation "Consider FBN1 gene
    sequencing", and the Marfan management block -- annual aortic root imaging,
    beta-blocker/ARB therapy, ophthalmology referral, avoid contact sports and
    isometric exercise -- was never emitted. That is the canonical young-athlete
    aortic dissection scenario.

    The family-history criteria remain unimplemented on this branch because
    Loeys et al. 2010 could not be retrieved, so the test asserts the workflow
    DECLARES the gap and escalates rather than returning a bare "not met".
    """
    result = connective.execute({"ghent_criteria": {
        "ectopia_lentis": True, "family_history": True,
    }})
    findings = " ".join(result.findings)
    assert "family-history criteria" in findings
    assert "NOT evaluated" in findings
    assert any("family-history arm" in r for r in result.recommendations)
    assert any("Marfan" in t for t in result.cross_agent_triggers)


def test_no_family_history_ghent_criterion_still_fires(connective):
    """Aortic Z >= 2 plus ectopia lentis is still a Marfan diagnosis."""
    result = connective.execute({"ghent_criteria": {
        "aortic_root_z": 3.1, "ectopia_lentis": True,
    }})
    assert any("Ghent criteria MET" in f for f in result.findings)
    assert any("aortic root imaging" in r for r in result.recommendations)


def test_classical_eds_does_not_get_the_hypermobility_panel(connective):
    """An unrecognised EDS subtype must not silently take the mildest panel.

    Clinical reference case: a patient recorded as ``eds_subtype="classical
    EDS"``. ``_EDS_PANELS.get(subtype, hypermobility)`` returned the
    hypermobility panel -- *TNXB* only -- so the classical-EDS genes *COL5A1*
    and *COL5A2*, which are in the same table, were never ordered.
    """
    result = connective.execute({"eds_subtype": "classical EDS"})
    recs = " ".join(result.recommendations)
    assert "TNXB" not in recs or "COL5A1" in recs
    assert any("not a recognised key" in f for f in result.findings)


def test_unrecognised_eds_subtype_offers_the_full_differential(connective):
    """With no subtype match, the comprehensive panel is recommended."""
    result = connective.execute({"eds_subtype": "some unknown subtype"})
    recs = " ".join(result.recommendations)
    assert "COL5A1" in recs and "COL3A1" in recs


def test_recognised_eds_subtype_still_gets_its_panel(connective):
    """A known subtype key still selects its own panel."""
    workflow = ConnectiveTissueWorkflow()
    subtype = next(iter(workflow._EDS_PANELS))
    result = workflow.execute({"eds_subtype": subtype})
    assert any(subtype in r for r in result.recommendations)
