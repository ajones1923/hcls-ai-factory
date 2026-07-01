"""Tests for the TSC-Variant Curator: ACMG combinatorial rules + real recovery (PRD §3 FR-VC-*)."""
from __future__ import annotations

from src.agents.base import AgentContext
from src.agents.variant_curator.acmg import Criterion, classify
from src.agents.variant_curator.annotation import annotate, assign_criteria
from src.agents.variant_curator.runner import VariantCurator
from src.cohort.build import build_cohort
from src.cohort.loader import featured_map, load_patient


# ── ACMG combinatorial classifier (Richards 2015 Table 5) ──────────────────
def C(code, bucket):
    return Criterion(code, bucket)


def test_classify_pathogenic_pvs1_plus_strong():
    cls, _ = classify([C("PVS1", "PVS"), C("PS1", "PS"), C("PM2", "PM")])
    assert cls == "Pathogenic"


def test_classify_likely_pathogenic_strong_plus_moderate():
    # Patient A pattern: PVS1 downgraded to Strong (mosaic) + PM2 + PP4
    cls, rule = classify([C("PVS1", "PS"), C("PM2", "PM"), C("PP4", "PP")])
    assert cls == "Likely Pathogenic"
    assert "Strong" in rule


def test_classify_vus_moderate_plus_supporting():
    cls, _ = classify([C("PM2", "PM"), C("PP4", "PP")])
    assert cls == "Variant of Uncertain Significance"


def test_classify_pvs1_plus_moderate_is_likely_pathogenic():
    cls, _ = classify([C("PVS1", "PVS"), C("PM2", "PM")])
    assert cls == "Likely Pathogenic"


def test_classify_benign_and_conflict():
    assert classify([C("BA1", "BA")])[0] == "Benign"
    assert classify([C("BS1", "BS"), C("BP4", "BP")])[0] == "Likely Benign"
    # conflicting evidence resolves to VUS
    assert classify([C("PVS1", "PVS"), C("PM2", "PM"), C("BA1", "BA")])[0] == \
        "Variant of Uncertain Significance"


def test_annotation_assigns_pvs1_strong_for_mosaic_null():
    ann = annotate({"gene": "TSC2", "consequence": "frameshift", "hgvsc": "c.4255del"})
    crits = assign_criteria(ann, mosaic=True, phenotype_consistent=True)
    codes = {c.code: c.bucket for c in crits}
    assert codes["PVS1"] == "PS"      # downgraded for mosaic
    assert codes["PM2"] == "PM"
    assert codes["PP4"] == "PP"


# ── real recovery on the cohort ────────────────────────────────────────────
def test_patient_a_mosaic_recovery_likely_pathogenic(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    fmap = featured_map(tmp_path)
    cohort = load_patient(fmap["A"], tmp_path)
    # point the curator at the test cohort dir
    import config.settings as cs
    cs.settings.COHORT_DIR = tmp_path

    out = VariantCurator().run(AgentContext(patient_id=fmap["A"], cohort=cohort))
    p = out.data["primary"]
    assert p["gene"] == "TSC2"
    assert p["mosaic"] is True and p["recovered"] is True
    assert abs(p["vaf"] - 0.083) < 0.01
    assert p["acmg_classification"] == "Likely Pathogenic"
    assert {c["code"] for c in p["acmg_criteria"]} >= {"PVS1", "PM2", "PP4"}
    assert p["ddpcr_recommended"] is True
    assert out.data["review_status"] == "pending_molecular_geneticist"
    assert out.data["source"] == "vcf"


def test_patient_b_germline_pathogenic(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    fmap = featured_map(tmp_path)
    cohort = load_patient(fmap["B"], tmp_path)
    import config.settings as cs
    cs.settings.COHORT_DIR = tmp_path

    out = VariantCurator().run(AgentContext(patient_id=fmap["B"], cohort=cohort))
    p = out.data["primary"]
    assert p["gene"] == "TSC2"
    assert p["mosaic"] is False
    assert p["acmg_classification"] == "Pathogenic"   # known recurrent null + PM2


def test_nmi_patient_no_variant_recommends_tissue(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    import config.settings as cs
    cs.settings.COHORT_DIR = tmp_path
    from src.cohort.loader import load_manifest
    nmi = next(p for p in load_manifest(tmp_path)["patients"] if p["zygosity"] == "none")
    cohort = load_patient(nmi["patient_id"], tmp_path)
    out = VariantCurator().run(AgentContext(patient_id=nmi["patient_id"], cohort=cohort))
    assert out.data["primary"] is None
    assert "tissue sequencing" in out.data["note"]
