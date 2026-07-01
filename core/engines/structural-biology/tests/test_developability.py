"""Tests for protein developability + design (B3/B4). Real computed metrics, no model."""
import pytest
from developability import DevelopabilityScorer, develop_flags, SequenceError

STABLE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR"


class TestScore:
    def test_real_metrics_present(self):
        s = DevelopabilityScorer().score(STABLE)
        m = s["metrics"]
        assert all(k in m for k in ("gravy", "instability_index", "isoelectric_point", "aromaticity"))
        assert s["verdict"] in ("clean", "caution", "high-risk")

    def test_flags_unstable(self):
        flags, verdict = develop_flags({"instability_index": 55, "gravy": 0.0, "isoelectric_point": 7.0})
        assert any("unstable" in f for f in flags)

    def test_flags_clean(self):
        flags, verdict = develop_flags({"instability_index": 30, "gravy": -0.2, "isoelectric_point": 7.0})
        assert flags == [] and verdict == "clean"

    def test_rejects_bad_sequence(self):
        with pytest.raises(SequenceError):
            DevelopabilityScorer().score("XZ123")


class TestOptimize:
    def test_proposes_improving_mutations(self):
        o = DevelopabilityScorer().optimize(STABLE, n=3)
        assert o["n_proposals"] >= 1
        # every proposal must actually lower the instability index
        for p in o["top"]:
            assert p["delta"] < 0 and p["instability_index"] < o["baseline_instability"]
        # sorted most-improving first
        assert o["top"] == sorted(o["top"], key=lambda x: x["delta"])
