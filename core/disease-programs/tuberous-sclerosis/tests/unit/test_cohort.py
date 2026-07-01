"""Tests for the synthetic cohort pipeline (PRD §3 FR-CG-*; NFR-R-3 determinism)."""
from __future__ import annotations

from src.cohort.build import build_cohort
from src.cohort.genomic import parse_variants
from src.cohort.loader import featured_map, load_manifest, load_patient


def test_composition_matches_epidemiology(tmp_path):
    m = build_cohort(out_dir=tmp_path, seed=42)
    assert m["n_patients"] == 50
    comp = m["composition"]
    assert comp["TSC2_germline"] == 30
    assert comp["TSC1_germline"] == 12
    assert comp["TSC2_mosaic"] == 5
    assert comp["TSC1_mosaic"] == 2
    assert comp["NMI_none"] == 1


def test_deterministic_same_seed_same_hash(tmp_path):
    a = build_cohort(out_dir=tmp_path / "a", seed=42)
    b = build_cohort(out_dir=tmp_path / "b", seed=42)
    assert a["cohort_hash"] == b["cohort_hash"]
    c = build_cohort(out_dir=tmp_path / "c", seed=43)
    assert c["cohort_hash"] != a["cohort_hash"]


def test_featured_patients_present_and_correct(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    fmap = featured_map(tmp_path)
    assert set(fmap) == {"A", "B", "C"}

    a = load_patient(fmap["A"], tmp_path)
    assert a["genomics"]["gene"] == "TSC2"
    assert a["genomics"]["tissue_vaf"] == 0.083          # mosaic recovery target
    assert "Likely Pathogenic" in a["genomics"]["expected_acmg"]

    b = load_patient(fmap["B"], tmp_path)
    assert b["sega_series_cm"] == [0.8, 1.1, 1.3]
    # TAND signal is now span-grounded in the notes (the real, offset-carrying source)
    academic_spans = [sp for note in b["notes"] for sp in note.get("spans", [])
                      if sp.get("kind") == "tand" and sp.get("cluster") == "academic"
                      and sp.get("polarity") == "present"]
    assert len(academic_spans) >= 3                      # flags TAND surveillance
    for note in b["notes"]:                              # every span verbatim-faithful to its text
        for sp in note.get("spans", []):
            assert note["text"][sp["start"]:sp["end"]] == sp["quote"]


def test_vcf_is_valid_and_carries_vaf(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    fmap = featured_map(tmp_path)
    recs = parse_variants(tmp_path / fmap["A"] / "variants.vcf")
    assert len(recs) == 1
    r = recs[0]
    assert r["gene"] == "TSC2"
    assert abs(r["af"] - 0.083) < 0.01                   # low-VAF mosaic in the VCF
    assert r["dp"] >= 200                                 # high coverage for detection


def test_mosaic_vafs_in_range_and_nmi_has_no_variant(tmp_path):
    build_cohort(out_dir=tmp_path, seed=42)
    man = load_manifest(tmp_path)
    mosaics = [p for p in man["patients"] if p["zygosity"] == "mosaic"]
    assert mosaics and all(0.04 <= p["vaf"] <= 0.12 for p in mosaics)
    nmi = next(p for p in man["patients"] if p["zygosity"] == "none")
    assert parse_variants(tmp_path / nmi["patient_id"] / "variants.vcf") == []
