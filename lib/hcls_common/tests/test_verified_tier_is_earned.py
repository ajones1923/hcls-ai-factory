"""`maturity: verified` may not contradict the capability's own description.

Found by the claim-integrity review of `genomic-foundation`, 2026-11-04, check 6.

The maturity matrix defines `verified` as "additionally proven against real, recorded input
(e.g. on real HG002)", and it is the single most load-bearing word on the public site: it is
the badge rendered *instead of* the serving status, so it is the one claim a reader sees first.

`acmg-secondary-findings` declared `verified` while its own description said the panel is a
curated representation of ACMG SF v3.2, that ACMG has since published v3.3, and to "validate
panel vs official publication before clinical use". The prose was admirably honest and the
machine-readable tier the badge is drawn from contradicted it. Of the nine capabilities
declaring `verified`, it was the only one in that position.

That is the inverse of the TSC finding the same day, where the flagship declared no tier at
all and carried its caveats in prose. Both are the same underlying defect: the honest words
are somewhere a generator does not read.

`variant-store` is the worked example of the tier used correctly -- "Verified on real HG002
PASS SNVs: Ts/Tv~2.0 (human ~2.0-2.1)" is a measurement against recorded input, in the entry.
"""
import json
import pathlib

import pytest

REGISTRY = (pathlib.Path(__file__).resolve().parents[1] / "capabilities.json")

pytestmark = pytest.mark.skipif(not REGISTRY.exists(),
                                reason="run from a source checkout; capabilities.json not found")

#: Phrases that say, in the entry's own words, that the thing is not settled yet. A
#: capability may say any of these -- it may not say them AND claim to be verified.
_CONTRADICTS_VERIFIED = (
    "before clinical use",
    "not yet validated",
    "pending validation",
    "awaiting validation",
)


def capabilities() -> list[dict]:
    return json.loads(REGISTRY.read_text())["capabilities"]


def verified() -> list[dict]:
    return [c for c in capabilities() if c.get("maturity") == "verified"]


class TestVerifiedIsNotSelfContradictory:
    def test_no_verified_capability_asks_to_be_validated_first(self):
        bad = [(c["id"], p) for c in verified()
               for p in _CONTRADICTS_VERIFIED
               if p in c.get("description", "").lower()]
        assert not bad, (
            "these declare maturity=verified while their own description says they are not "
            f"settled: {bad}. `verified` is rendered on the public site INSTEAD of the "
            "serving status, so it is the first claim a reader sees.")

    def test_the_acmg_panel_is_not_badged_verified(self):
        """The specific entry that prompted this. Its panel is knowingly superseded."""
        acmg = next(c for c in capabilities() if c["id"] == "acmg-secondary-findings")
        assert acmg.get("maturity") != "verified", (
            "the ACMG SF panel is a curated representation of v3.2, which v3.3 supersedes, "
            "and the entry itself says to validate before clinical use")
        assert acmg.get("maturity") == "research-use"

    def test_a_verified_capability_is_live_and_not_mock(self):
        """The registry's own honesty rule, restated here so this file owns the tier."""
        for c in verified():
            assert c.get("status") == "live", f"{c['id']} is verified but not live"


class TestTheTierStillMeansSomething:
    def test_variant_store_keeps_it_and_says_why(self):
        """If this ever fails, the guard above has been satisfied by deleting the evidence."""
        vs = next(c for c in capabilities() if c["id"] == "variant-store")
        assert vs.get("maturity") == "verified"
        assert "HG002" in vs["description"], (
            "variant-store is the worked example of an earned `verified`: it names the real "
            "recorded input it was proven against")

    def test_at_least_one_clinical_capability_still_carries_a_tier(self):
        """Guard against the fix being applied by stripping tiers platform-wide."""
        tiered = [c for c in capabilities() if c.get("maturity")]
        assert len(tiered) >= 10, f"only {len(tiered)} capabilities declare a maturity tier"
