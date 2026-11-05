"""Every engine and agent must declare an evidence tier, with a reason.

Written 2026-11-04. Of the 17 engines and agents, **15 were marked `live` with no maturity tier
at all.** `site_gen_matrix.py` badges a capability with its tier where one exists and otherwise
with its serving status, so those fifteen were badged a flat `live` on every generated surface --
the public site included -- while their own descriptions said things like "research-/trial-use
inputs, not routine diagnostics" and "drug-response predictions are hypotheses to test".

The project's own plain-English summary of itself was *more honest than its manifest*: the first
engines are solid and the rest have their foundations in place, which is an ordinary and
defensible position for an open research platform to be in. Nothing in the registry said so.
A reader who tried engine 6 expecting engine 1's depth would have been misled by a field, not by
a sentence -- and the field is what every badge, table and page is generated from.

This is the same defect this review keeps finding, one level up: the honest words were somewhere
no generator reads.

WHAT THIS GATE IS NOT. It does not judge whether a tier is *correct* -- no test can. It requires
that a tier and a basis exist, so the claim is explicit and attributable. `maturity_basis` is
what makes a tier checkable rather than merely asserted: it names the specific thing holding the
tier down ("three CPIC citations resolved to unrelated papers"), so a later reader can go and
disagree with it. The sibling tests guard the two directions a tier can lie:

  * test_verified_tier_is_earned.py   -- a tier claiming more than the description supports
  * test_caveats_reach_the_badge.py   -- a description caveating more than the tier admits
"""
import json
import pathlib

import pytest

REGISTRY = pathlib.Path(__file__).resolve().parents[1] / "capabilities.json"

SUBJECT_TYPES = {"engine", "agent"}
VALID_TIERS = {"verified", "preclinical", "research-use", "roadmap", "gated"}

# A basis shorter than this is a label, not a reason. The real ones run 150-300 characters;
# the bar is set low enough to permit a genuinely short true statement and high enough to
# reject "TBD", "see docs" or a repeat of the tier name.
MIN_BASIS = 40


def subjects():
    caps = json.loads(REGISTRY.read_text())["capabilities"]
    return [c for c in caps if c["type"] in SUBJECT_TYPES]


def test_the_roster_is_the_size_we_think_it_is():
    """If this drops, the rest of the file is checking a shorter list than it claims to."""
    found = subjects()
    assert len(found) >= 17, (
        f"expected at least 17 engines+agents (8 engines, 8 agents, 1 disease program), "
        f"found {len(found)}: {sorted(c['id'] for c in found)}")


@pytest.mark.parametrize("cap", subjects(), ids=lambda c: c["id"])
class TestEverySubjectIsTiered:
    def test_declares_a_tier(self, cap):
        """The finding in one assertion."""
        assert cap.get("maturity"), (
            f"{cap['id']} is marked {cap.get('status')} with no maturity tier. Every surface "
            f"generated from this manifest will badge it on serving status alone, which says "
            f"whether it runs and nothing about the evidence behind what it returns.")

    def test_the_tier_is_one_of_the_real_ones(self, cap):
        tier = cap.get("maturity")
        assert tier in VALID_TIERS, f"{cap['id']}: {tier!r} is not a MaturityTier value"

    def test_declares_why(self, cap):
        basis = (cap.get("maturity_basis") or "").strip()
        assert basis, (
            f"{cap['id']} declares maturity {cap.get('maturity')!r} with no basis. A tier on "
            f"its own is an assertion; a reader cannot tell whether someone judged the evidence "
            f"or someone needed the field filled in.")
        assert len(basis) >= MIN_BASIS, (
            f"{cap['id']}: basis is {len(basis)} chars ({basis!r}) -- too short to name what "
            f"holds the tier down")

    def test_the_basis_is_not_just_the_tier_restated(self, cap):
        """"research-use because it is for research use" is circular, not a basis."""
        basis = (cap.get("maturity_basis") or "").lower()
        tier = cap["maturity"].replace("-", " ")
        stripped = basis.replace(tier, "").replace(cap["maturity"], "")
        assert len(stripped) >= MIN_BASIS, (
            f"{cap['id']}: the basis says little beyond repeating the tier name")


class TestTheTiersSayWhatWeMeanThemToSay:
    def test_nothing_live_is_roadmap(self):
        """`roadmap` and `live` together is a contradiction a reader would catch instantly."""
        bad = [c["id"] for c in subjects()
               if c.get("maturity") == "roadmap" and c.get("status") == "live"]
        assert not bad, f"badged roadmap while serving: {bad}"

    def test_the_gradient_exists(self):
        """The point of tiering was to stop presenting 17 subjects as interchangeable.

        If every subject carried the same tier this file would pass while the manifest said
        nothing useful, so the distribution itself is asserted: more than one tier in use, and
        `verified` not handed out wholesale.
        """
        tiers = {c.get("maturity") for c in subjects()}
        assert len(tiers) >= 2, f"every subject is {tiers} -- that is not a gradient"
        verified = [c["id"] for c in subjects() if c.get("maturity") == "verified"]
        assert len(verified) <= 3, (
            f"{len(verified)} subjects claim `verified`: {verified}. That tier means proven "
            f"against real reference data with a recorded run, and it is the one a reader will "
            f"check first.")
