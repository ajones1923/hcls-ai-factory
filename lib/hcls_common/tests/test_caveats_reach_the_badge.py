"""A capability that caveats itself in prose must say so in the field the badge reads.

This is the structural finding of the first four subject reviews, which all failed the same
way: **the honest words were somewhere no generator reads.**

    TSC trajectory      honest module docstring   false API payload ("NHD mixed-effects")
    TSC tier            honest description        no maturity tier at all
    ACMG SF panel       honest description        maturity: verified
    imaging agent       honest description        no maturity tier at all

`site_gen_matrix.py` badges each capability with "the single most honest label": its evidence
tier where one is declared, **otherwise its serving status**. So a capability whose own
description says a modality is "planned/placeholder, not yet live" is badged a flat `live` on
the public site, and the caveat that exists is invisible to the one surface most readers see.

Four capabilities are in that position today (listed below). This test does not pick their
tiers -- which tier is right is a judgement with public consequences, and more than one is
defensible for most of them. It freezes the set, so a seventh cannot appear by accident while
the six are decided. The sibling `test_verified_tier_is_earned.py` guards the other direction,
where a declared tier claims more than the description supports.
"""
import json
import pathlib

import pytest

REGISTRY = (pathlib.Path(__file__).resolve().parents[1] / "capabilities.json")

pytestmark = pytest.mark.skipif(not REGISTRY.exists(),
                                reason="run from a source checkout; capabilities.json not found")

#: Phrases by which an entry tells the reader it is not simply "live".
CAVEATS = (
    "not yet live",
    "placeholder",
    "preclinical",
    "research-use",
    "research/augmentation",
    "never diagnostic",
    "before clinical use",
    "superseded",
    "not yet a request",
)

#: `live`, caveated in prose, no tier -- as measured on 2026-11-04 by the claim-integrity
#: review. Each is an open S4: the tier is the owner's call. The list may SHRINK freely; it
#: may not grow, which is the whole point of pinning it.
#:
#: A first draft of this list had six names, built with "validate" among the phrases above.
#: That matched any description telling a reader to validate something before clinical use,
#: which is good practice rather than a caveat on the capability's own evidence, and the
#: rot-guard below rejected the two it had swept in. The phrase list is deliberately narrow:
#: a loose one produces a long exemption list, and a long exemption list is how a guard stops
#: meaning anything.
#:
#: Recommended tiers, recorded so the decision has a starting point:
#:   imaging-intelligence-agent   research-use  (CT path is placeholder; MAISI never diagnostic)
#:   tuberous-sclerosis-engine    research-use  (flagship; gene-therapy arm preclinical)
#:   precision-biomarker-agent    research-use  (biological-age clocks are not clinical tests)
#:   genomics-engine              no tier needed -- its caveat is about serving SHAPE
#:                                ("not yet a request/response API"), not about evidence
ACCEPTED_UNTIERED = {
    "genomics-engine",
    "precision-biomarker-agent",
    "imaging-intelligence-agent",
    "tuberous-sclerosis-engine",
}


def capabilities() -> list[dict]:
    return json.loads(REGISTRY.read_text())["capabilities"]


def caveated_without_a_tier() -> dict[str, list[str]]:
    out = {}
    for c in capabilities():
        if c.get("status") != "live" or c.get("maturity"):
            continue
        hits = [p for p in CAVEATS if p in c.get("description", "").lower()]
        if hits:
            out[c["id"]] = hits
    return out


class TestTheSetCannotGrow:
    def test_no_new_capability_caveats_itself_without_declaring_a_tier(self):
        found = caveated_without_a_tier()
        new = {k: v for k, v in found.items() if k not in ACCEPTED_UNTIERED}
        assert not new, (
            "these are status=live with no maturity tier, and their own descriptions say "
            f"they are qualified: {new}. The public matrix badges a capability on its tier "
            "where one exists and on serving status otherwise, so the caveat never reaches "
            "the badge. Declare a tier, or add it to ACCEPTED_UNTIERED with a reason.")

    def test_the_accepted_list_does_not_outlive_the_problem(self):
        """A name that no longer qualifies must be removed, so the list cannot rot."""
        stale = ACCEPTED_UNTIERED - set(caveated_without_a_tier())
        assert not stale, (
            f"{stale} no longer need an exemption -- they declare a tier now, or their "
            "description no longer caveats them. Remove them from ACCEPTED_UNTIERED.")


class TestTheBadgeRuleItDependsOn:
    def test_the_matrix_still_falls_back_to_serving_status(self):
        """If this changes, the reasoning above needs rewriting rather than the test muting."""
        gen = REGISTRY.parents[2] / "scripts" / "site_gen_matrix.py"
        if not gen.exists():
            pytest.skip("site_gen_matrix.py not in this checkout")
        src = gen.read_text()
        assert "otherwise its" in src and "status" in src, (
            "site_gen_matrix.py no longer documents the tier-else-status fallback this "
            "guard is built on")
