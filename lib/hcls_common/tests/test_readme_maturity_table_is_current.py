"""The README's maturity table must still match the registry.

The table exists because the README is the surface most readers judge the project by, and a
hand-maintained table of 17 subjects is exactly the kind of thing that is true the day it is
written and wrong a month later. A stale honesty table is worse than none: it is an overclaim
with a timestamp on it.

So the registry is the source of truth and this test is the link. It checks membership, not
formatting -- the table may be reworded, regrouped or restyled freely, but a subject cannot
appear under a tier it does not hold, and no subject may be missing from it.
"""
import json
import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parents[3]
README = REPO / "README.md"
REGISTRY = REPO / "lib" / "hcls_common" / "capabilities.json"

BEGIN = "<!-- MATURITY-TABLE:BEGIN"
END = "<!-- MATURITY-TABLE:END -->"


def table_text() -> str:
    text = README.read_text()
    assert BEGIN in text and END in text, (
        "the README maturity table markers are gone. If the table was removed on purpose, "
        "delete this test in the same commit and say why; do not leave a guard pointing at "
        "nothing.")
    return text.split(BEGIN, 1)[1].split(END, 1)[0]


def subjects():
    caps = json.loads(REGISTRY.read_text())["capabilities"]
    return [c for c in caps if c["type"] in {"engine", "agent"}]


def rows():
    """tier -> the subject cell, for each row of the table."""
    out = {}
    for line in table_text().splitlines():
        m = re.match(r"\s*\|\s*\*\*`([a-z-]+)`\*\*\s*\|([^|]*)\|(.*)\|\s*$", line)
        if m:
            out[m.group(1)] = m.group(3)
    return out


class TestTheTableCoversTheRoster:
    def test_every_tier_in_use_has_a_row(self):
        in_use = {c["maturity"] for c in subjects()}
        assert in_use <= set(rows()), (
            f"tiers held by a subject but absent from the README table: "
            f"{sorted(in_use - set(rows()))}")

    def test_no_row_claims_a_tier_nobody_holds(self):
        in_use = {c["maturity"] for c in subjects()}
        extra = set(rows()) - in_use
        assert not extra, f"README table has rows for unused tiers: {sorted(extra)}"

    def test_every_subject_appears_under_its_own_tier(self):
        """The assertion that actually catches drift.

        Matched on a distinctive word from each subject's registry name rather than the full
        string, so the table can abbreviate ("CAR-T" for "CAR-T Intelligence Agent") without
        being wrong. The word has to be in the right row.
        """
        table = rows()
        misplaced = []
        for cap in subjects():
            # the most distinctive token in the name: longest word that is not a role noun
            noise = {"intelligence", "agent", "engine", "precision", "clinical", "analysis",
                     "foundations", "large-molecule", "/"}
            words = [w for w in re.split(r"[\s/]+", cap["name"]) if w.lower() not in noise]
            if words:
                token = max(words, key=len)
            else:
                # Every word is a role noun -- "Precision Intelligence Engine". Fall back to the
                # leading two words, which is what an abbreviated table cell would keep.
                token = " ".join(re.split(r"[\s/]+", cap["name"])[:2])
            cell = table.get(cap["maturity"], "")
            if token.lower() not in cell.lower():
                misplaced.append(
                    f"  {cap['id']} is {cap['maturity']} but {token!r} is not in that row")
        assert not misplaced, (
            "the README maturity table no longer matches the registry:\n" + "\n".join(misplaced)
            + "\n\nRegenerate the table from lib/hcls_common/capabilities.json.")

    def test_the_counts_in_the_prose_match(self):
        """The sentence under the table states the distribution; it drifts as easily as the table."""
        text = README.read_text()
        counts = {}
        for cap in subjects():
            counts[cap["maturity"]] = counts.get(cap["maturity"], 0) + 1
        spelled = {1: "one", 2: "two", 3: "three", 13: "thirteen", 14: "fourteen"}
        n = counts.get("research-use", 0)
        word = spelled.get(n, str(n))
        assert re.search(rf"{word}\s+are\s+\n?research-use", text) or f"{word} are" in text, (
            f"{n} subjects are research-use; the README prose should say so (looked for "
            f"{word!r}). Update the sentence under the table.")
