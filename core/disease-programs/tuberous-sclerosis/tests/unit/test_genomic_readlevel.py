"""Read-level genomic evidence + artifact rejection (FR-CG-2 / FR-VC-4)."""
from src.cohort.genomic import (
    count_filtered_artifacts, parse_variants, write_vcf,
)
from src.cohort.spec import build_roster


def _spec(tag):
    return next(s for s in build_roster() if s.featured == tag)


def test_vcf_carries_strand_resolved_read_level_fields(tmp_path):
    p = write_vcf(_spec("B"), tmp_path / "b.vcf", seed=1)
    recs = parse_variants(p)                      # PASS only
    assert len(recs) == 1
    r = recs[0]
    assert r["alt_reads"] > 0 and r["ref_reads"] > 0 and r["dp"] > 0
    assert 0.3 <= r["strand_balance"] <= 0.7      # a true variant is strand-balanced
    assert r["filter"] == "PASS"


def test_strand_biased_artifact_is_filtered_not_reported(tmp_path):
    p = write_vcf(_spec("B"), tmp_path / "b.vcf", seed=1)
    # the artifact exists in the file but is FILTER=strand_bias
    assert count_filtered_artifacts(p) == 1
    all_recs = parse_variants(p, include_filtered=True)
    art = [r for r in all_recs if r["filter"] == "strand_bias"]
    assert art and (art[0]["strand_balance"] >= 0.8 or art[0]["strand_balance"] <= 0.2)
    # the default reader (what the curator uses) never surfaces it
    assert all(r["filter"] == "PASS" for r in parse_variants(p))


def test_nmi_patient_has_no_pass_variant_only_artifact(tmp_path):
    nmi = next(s for s in build_roster() if s.gene is None)
    p = write_vcf(nmi, tmp_path / "nmi.vcf", seed=1)
    assert parse_variants(p) == []               # nothing reportable in blood — the diagnostic gap
    assert count_filtered_artifacts(p) == 1      # only the rejected artifact


def test_curator_reports_read_level_and_artifacts_filtered():
    from app._engine import featured, get_engine
    orch, _ = get_engine()
    vi = (orch.store.projection(featured()["B"]).get("variant_interp") or {}).get("primary") or {}
    assert vi.get("artifact_assessment") == "pass"
    assert vi.get("alt_reads") and vi.get("depth")
