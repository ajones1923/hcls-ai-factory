"""Tests for the Artifact envelope (PF-3)."""
import pytest

from hcls_common.artifact import (
    Artifact, Honesty, Maturity, Provenance,
    new_artifact, combine_honesty, weakest_maturity,
)
from hcls_common.capability_registry import ArtifactShape


# ── maturity ordering (load-bearing for the non-inflation law) ───────────────
class TestMaturityOrder:
    def test_ranks_strongest_to_weakest(self):
        assert Maturity.live.rank == 0
        assert Maturity.roadmap.rank == 4
        assert Maturity.live.rank < Maturity.preclinical.rank < Maturity.research_use.rank \
            < Maturity.hypothesis_only.rank < Maturity.roadmap.rank

    def test_weakest_picks_most_cautious(self):
        assert weakest_maturity([Maturity.live, Maturity.hypothesis_only, Maturity.preclinical]) \
            is Maturity.hypothesis_only
        assert weakest_maturity([]) is Maturity.live   # empty → live


# ── non-inflation combine (the "label = weakest leg" rule) ───────────────────
class TestCombineHonesty:
    def test_composite_takes_weakest_maturity_and_unions_labels(self):
        a = Honesty(Maturity.live, labels=["decision-support"], requires=["clinician-review"])
        b = Honesty(Maturity.hypothesis_only, labels=["rna_inferred"], requires=["tumor-board"])
        c = combine_honesty([a, b])
        assert c.maturity is Maturity.hypothesis_only            # weakest leg wins
        assert c.labels == ["decision-support", "rna_inferred"]  # union, sorted
        assert c.requires == ["clinician-review", "tumor-board"]

    def test_cannot_inflate_a_hypothesis_into_live(self):
        # a live-looking downstream may not be MORE confident than its weakest input
        combined = combine_honesty([Honesty(Maturity.live), Honesty(Maturity.hypothesis_only)])
        assert combined.maturity.rank >= Maturity.hypothesis_only.rank

    def test_empty_is_live(self):
        assert combine_honesty([]).maturity is Maturity.live


# ── envelope construction + round-trip ───────────────────────────────────────
class TestArtifact:
    def test_new_artifact_fills_id_and_ts(self):
        art = new_artifact(
            ArtifactShape.PGX_DIPLOTYPES, {"CYP2C19": "*2/*2"},
            producer_id="pgx-producer", maturity=Maturity.live, labels=["decision-support"],
        )
        assert art.id and art.provenance.ts                      # auto-filled
        assert art.shape is ArtifactShape.PGX_DIPLOTYPES
        assert art.honesty.maturity is Maturity.live
        assert art.provenance.producer_id == "pgx-producer"

    def test_round_trip(self):
        art = new_artifact(
            ArtifactShape.MTB_PACKET, {"tier": 1},
            producer_id="precision-oncology-engine",
            maturity=Maturity.preclinical, labels=["pediatric-caution"],
            requires=["tumor-board"], inputs=["a1", "a2"], patient_id="P0", run_id="R0",
            serving="native", id="fixed-id", ts="2026-07-06T00:00:00+00:00",
        )
        art2 = Artifact.from_dict(art.to_dict())
        assert art2 == art
        assert art2.provenance.inputs == ["a1", "a2"]
        assert art2.patient_id == "P0" and art2.run_id == "R0"

    def test_defaults_are_safe(self):
        prov = Provenance(producer_id="x")
        assert prov.serving == "native" and prov.inputs == []
        h = Honesty(Maturity.roadmap)
        assert h.labels == [] and h.requires == []
