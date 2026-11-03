"""Every row of ACMG/AMP 2015 Table 5, pinned.

Richards S et al., "Standards and guidelines for the interpretation of sequence variants",
Genet Med 2015;17(5):405-424, Table 5, gives the combining rules as an exact list of evidence
combinations. `acmg.classify` implements them, and a reading of it against the table matched --
but a reading is not a test, and these rules are the last step before a variant is called
Pathogenic for a patient.

Each case below is one row of Table 5, built from criterion buckets rather than from real
variants, so the combination under test is unambiguous. The near-miss cases matter as much as
the hits: one criterion short of a rule must come back VUS, not the next category down.
"""
import pytest

from src.agents.variant_curator.acmg import Criterion, classify


def crits(**buckets):
    """Build a criterion list from bucket counts, e.g. crits(PVS=1, PM=2)."""
    out = []
    for bucket, n in buckets.items():
        for i in range(n):
            out.append(Criterion(f"{bucket}{i + 1}", bucket, "synthetic"))
    return out


def call(**buckets):
    return classify(crits(**buckets))[0]


# =====================================================================
# Table 5 — Pathogenic
#   (i)  1 PVS1 AND (>=1 PS | >=2 PM | 1 PM + 1 PP | >=2 PP)
#   (ii) >=2 PS
#   (iii) 1 PS AND (>=3 PM | 2 PM + >=2 PP | 1 PM + >=4 PP)
# =====================================================================
class TestPathogenic:
    @pytest.mark.parametrize("buckets", [
        {"PVS": 1, "PS": 1},
        {"PVS": 1, "PM": 2},
        {"PVS": 1, "PM": 1, "PP": 1},
        {"PVS": 1, "PP": 2},
    ])
    def test_very_strong_plus_supporting_combination(self, buckets):
        assert call(**buckets) == "Pathogenic"

    def test_two_strong(self):
        assert call(PS=2) == "Pathogenic"

    @pytest.mark.parametrize("buckets", [
        {"PS": 1, "PM": 3},
        {"PS": 1, "PM": 2, "PP": 2},
        {"PS": 1, "PM": 1, "PP": 4},
    ])
    def test_one_strong_plus_moderate_combination(self, buckets):
        assert call(**buckets) == "Pathogenic"


# =====================================================================
# Table 5 — Likely pathogenic
#   1 PVS1 + 1 PM · 1 PS + 1-2 PM · 1 PS + >=2 PP · >=3 PM ·
#   2 PM + >=2 PP · 1 PM + >=4 PP
# =====================================================================
class TestLikelyPathogenic:
    @pytest.mark.parametrize("buckets", [
        {"PVS": 1, "PM": 1},
        {"PS": 1, "PM": 1},
        {"PS": 1, "PM": 2},
        {"PS": 1, "PP": 2},
        {"PM": 3},
        {"PM": 2, "PP": 2},
        {"PM": 1, "PP": 4},
    ])
    def test_each_likely_pathogenic_row(self, buckets):
        assert call(**buckets) == "Likely Pathogenic"

    def test_one_strong_and_three_moderate_is_pathogenic_not_likely(self):
        """The 1-2 PM row stops at two; three moves it up to Pathogenic."""
        assert call(PS=1, PM=2) == "Likely Pathogenic"
        assert call(PS=1, PM=3) == "Pathogenic"


# =====================================================================
# Table 5 — Benign and likely benign
# =====================================================================
class TestBenign:
    def test_ba1_is_stand_alone(self):
        assert call(BA=1) == "Benign"

    def test_two_strong_benign(self):
        assert call(BS=2) == "Benign"

    def test_one_strong_plus_one_supporting_is_likely_benign(self):
        assert call(BS=1, BP=1) == "Likely Benign"

    def test_two_supporting_is_likely_benign(self):
        assert call(BP=2) == "Likely Benign"

    def test_one_strong_benign_alone_is_not_enough(self):
        assert call(BS=1) == "Variant of Uncertain Significance"

    def test_one_supporting_benign_alone_is_not_enough(self):
        assert call(BP=1) == "Variant of Uncertain Significance"


# =====================================================================
# Near-misses: one criterion short must be VUS, not the next tier down
# =====================================================================
class TestNearMisses:
    @pytest.mark.parametrize("buckets", [
        {"PVS": 1},              # very strong alone
        {"PS": 1},               # one strong alone
        {"PM": 2},               # two moderate alone
        {"PM": 1, "PP": 3},      # one short of 1 PM + 4 PP
        {"PM": 2, "PP": 1},      # one short of 2 PM + 2 PP
        {"PP": 4},               # supporting only
        {},                      # nothing at all
    ])
    def test_insufficient_evidence_is_vus(self, buckets):
        assert call(**buckets) == "Variant of Uncertain Significance"

    def test_pvs1_alone_is_not_likely_pathogenic(self):
        """A null variant in a LOF gene with nothing else is VUS, not Likely Pathogenic."""
        assert call(PVS=1) == "Variant of Uncertain Significance"


# =====================================================================
# Conflicting evidence
# =====================================================================
class TestConflictingEvidence:
    def test_pathogenic_and_benign_evidence_together_is_vus(self):
        verdict, rule = classify(crits(PVS=1, PS=1, BP=2))
        assert verdict == "Variant of Uncertain Significance"
        assert "conflicting" in rule

    def test_likely_pathogenic_against_likely_benign_is_vus(self):
        verdict, rule = classify(crits(PM=3, BP=2))
        assert verdict == "Variant of Uncertain Significance"
        assert "conflicting" in rule

    def test_a_stand_alone_benign_still_conflicts_rather_than_winning(self):
        """BA1 is stand-alone for Benign, but pathogenic evidence beside it is a conflict a
        curator has to resolve, not something to resolve silently either way."""
        verdict, _ = classify(crits(PVS=1, PS=1, BA=1))
        assert verdict == "Variant of Uncertain Significance"


class TestTheRuleIsReported:
    def test_every_call_names_the_rule_that_fired(self):
        for buckets in ({"PS": 2}, {"PM": 3}, {"BA": 1}, {"BP": 2}, {}):
            verdict, rule = classify(crits(**buckets))
            assert rule, f"{buckets} -> {verdict} with no rule stated"
