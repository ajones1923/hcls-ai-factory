#!/usr/bin/env python3
"""Generate each subject's maturity statement: the seven claim-integrity checks, measured.

A reader who wants to know whether this platform's claims hold should not have to take
anyone's word for it. That is the whole deliverable: one statement per subject, built from
measurements, ending in the commands that reproduce every number in it.

The rule that makes it worth reading
------------------------------------
**A check with no answer is reported as NOT ANSWERED, never as a pass.** Four of the seven
can be measured by machine and are filled here. Three (2, 4 and 6) require someone to read
the code, and this script will not pretend otherwise -- it emits them unanswered, with the
starting evidence a reviewer needs, until a human records a result.

That is deliberate. A generator that quietly printed "pass" for what it had not checked would
be the exact defect this program exists to find, committed by the tool meant to find it.

Relationship to the capability maturity matrix
----------------------------------------------
`site_gen_matrix.py` publishes what the registry **declares** -- status and evidence tier,
honesty by construction. This publishes what was **measured**. A subject can declare `live`
and measure unreachable, and the gap between those two pages is the thing worth seeing.

Deliberately NOT a CI gate
--------------------------
`--check` exists for local use and is not wired into CI. The output carries a generation date
and the short commit, both of which differ between a developer box and a runner building a
merge commit, so a gate on it would fail for reasons that have nothing to do with the claims
it describes -- and a gate that fails for the wrong reason is one people learn to switch off.

What does block a merge is the four checks underneath: `check_endpoints.py`,
`check_corpus_ratchet.py`, `check_test_depth.py` and `check_docs.py --all` all run with
`--enforce` in CI. A statement can go out of date; the measurement it quotes cannot regress
without a build going red.

Usage:
    .venv/bin/python scripts/gen_maturity.py                     # all 17 + the index
    .venv/bin/python scripts/gen_maturity.py cardiology
    .venv/bin/python scripts/gen_maturity.py --tests /tmp/subject_tests.json
    .venv/bin/python scripts/gen_maturity.py --check             # local staleness check
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "lib"))

from subject_facts import BASE, ROSTER, build  # noqa: E402

ENDPOINTS = ROOT / "docs" / "build" / "endpoint_snapshot.json"
ENDPOINT_BASELINE = ROOT / "docs" / "build" / "endpoint_baseline.json"
CORPUS = ROOT / "docs" / "build" / "corpus_snapshot.json"
CLINICAL_REVIEW = ROOT / "CLINICAL_REVIEW.md"

PASS, FAIL, UNANSWERED, NA = "PASS", "FAIL", "NOT ANSWERED", "n/a"

DISCLAIMER = ("> **Decision support for a qualified clinician — never autonomous diagnosis "
              "or prescribing.**")

#: Checks 2, 4 and 6 need a person. What each needs, so an unanswered row is actionable
#: rather than merely blank.
MANUAL = {
    2: ("No displayed value originates anywhere but a measurement",
        "Trace every field in every API response and dashboard panel to its source. "
        "A modeled value must carry an explicit `*_is_modeled` flag."),
    4: ("Every clinical or scientific constant is cited or recorded as open",
        "Each threshold, range and algorithm coefficient carries a primary source in a code "
        "comment, or a numbered `CLINICAL_REVIEW.md` section."),
    6: ("The declared maturity tier matches what the code does",
        "Read the declared tier against the code path it describes. Status and tier are "
        "orthogonal; both must hold."),
}


def _load(path: pathlib.Path) -> dict:
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------- check 1
def check_endpoints(facts: dict) -> tuple[str, str, str]:
    """Does every `live` capability of this subject answer where the registry says?"""
    snap = _load(ENDPOINTS)
    results = snap.get("results", {})
    if not results:
        return UNANSWERED, "", "No probe recorded. Run `scripts/probe_live_endpoints.py --write`."
    mine = [c for c in facts["capabilities"] if c.get("status") == "live" and c.get("endpoint")]
    if not mine:
        return NA, snap.get("measured_at", ""), "No `live` capability declares an endpoint."
    unprobed = [c["id"] for c in mine if c["id"] not in results]
    if unprobed:
        return UNANSWERED, snap.get("measured_at", ""), (
            f"Never probed: {', '.join(unprobed)}. An unprobed endpoint is unknown, not up.")
    down = [c for c in mine if results[c["id"]].get("reachable") is False]
    if down:
        known = set(_load(ENDPOINT_BASELINE).get("known_down", []))
        tag = " (accepted in the endpoint baseline)" if all(c["id"] in known for c in down) else ""
        return FAIL, snap.get("measured_at", ""), (
            "Advertised and serving nothing: "
            + ", ".join(f"`{c['id']}` at {c['endpoint']}" for c in down) + tag + ".")
    return PASS, snap.get("measured_at", ""), (
        f"{len(mine)} of {len(mine)} answered: "
        + ", ".join(f"{c['endpoint']}" for c in mine) + ".")


# --------------------------------------------------------------------------- check 3
def corpus_rows(slug: str) -> tuple[str, str, str, list[tuple[str, float, int, int]]]:
    """Retrieval coverage per workflow, from the committed corpus snapshot."""
    snap = _load(CORPUS)
    counts = {str(k): int(v) for k, v in snap.get("collections", {}).items()}
    if not counts:
        return UNANSWERED, "", "No corpus snapshot recorded.", []
    try:
        from check_corpus_weights import WEIGHT_TABLES, weight_maps

        from hcls_common.corpus_coverage import coverage
    except Exception as exc:
        return UNANSWERED, "", f"Weight tables unavailable: {exc}", []
    if slug not in WEIGHT_TABLES:
        return UNANSWERED, snap.get("measured_at", ""), (
            "No weight table is declared for this subject, so its retrieval is invisible to "
            "the corpus gate. That is a gap in the gate, not evidence of coverage."), []
    maps, err = weight_maps(slug)
    if err:
        return UNANSWERED, snap.get("measured_at", ""), f"Weight table did not import: {err}", []
    rows = []
    for workflow, weights in sorted(maps.items()):
        v = coverage(weights, {c: counts.get(c, 0) for c in weights})
        filled = sum(1 for c in weights if counts.get(c, 0) > 0)
        rows.append((workflow, v.coverage, filled, len(weights)))
    worst = min(r[1] for r in rows)
    verdict = PASS if worst >= 0.75 else FAIL
    note = (f"{sum(1 for r in rows if r[1] < 0.75)} of {len(rows)} workflows retrieve from "
            f"under 75% of the evidence they declare; worst {worst:.0%}.")
    return verdict, snap.get("measured_at", ""), note, rows


# --------------------------------------------------------------------------- check 5
def check_depth(facts: dict, tests: pathlib.Path | None) -> tuple[str, str, str]:
    if not tests or not tests.exists():
        return UNANSWERED, "", (
            "No test results supplied. Run `scripts/run_all_tests.py --json <path>` and pass "
            "it with `--tests`.")
    try:
        from check_test_depth import CLINICAL, FLOOR, measure
    except Exception as exc:
        return UNANSWERED, "", f"Depth checker unavailable: {exc}"
    key = (facts["kind"], facts["slug"])
    row = next((r for r in measure(tests) if (r["kind"], r["name"]) == key), None)
    if row is None:
        return UNANSWERED, "", "This subject is absent from the supplied test results."
    when = datetime.datetime.fromtimestamp(
        tests.stat().st_mtime, datetime.timezone.utc).replace(microsecond=0).isoformat()
    # `loc` here is TRACKED, non-test, non-vendored Python only, which is deliberately not the
    # same number as the whole-tree count in section 2. Say which, or a reader reconciling the
    # two concludes one of them is wrong.
    if key not in CLINICAL:
        return NA, when, (
            f"{row['tests']} tests over {row['loc']:,} source lines ({row['per100']}/100). "
            f"The floor applies to subjects that emit clinical output; this one does not.")
    verdict = PASS if row["per100"] >= FLOOR else FAIL
    return verdict, when, (
        f"{row['tests']} tests over {row['loc']:,} source lines (tracked, non-test) = "
        f"{row['per100']} per 100 (clinical floor {FLOOR}).")


# --------------------------------------------------------------------------- check 7
def check_docs(facts: dict) -> tuple[str, str, str]:
    """`check_docs.py --all` is platform-wide; report it scoped to this subject's files."""
    base = ROOT / facts["path"]
    try:
        tracked = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", f"{facts['path']}/*.md"],
            capture_output=True, text=True, check=True).stdout.split()
    except Exception:
        tracked = [str(p) for p in base.rglob("*.md")]
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "check_docs.py"), "--all"],
                       capture_output=True, text=True)
    when = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()
    if r.returncode != 0:
        return FAIL, when, (
            "`check_docs.py --all` fails platform-wide; this subject's "
            f"{len(tracked)} markdown files are inside that scope.")
    return PASS, when, (
        f"{len(tracked)} tracked markdown files, inside the platform-wide "
        "`check_docs.py --all` pass: every internal link resolves and every referenced "
        "script exists.")


# --------------------------------------------------------------------------- clinical review
def clinical_sections(facts: dict) -> list[str]:
    """CLINICAL_REVIEW.md sections that name this subject. A keyword match, not an index."""
    if not CLINICAL_REVIEW.exists():
        return []
    words = {facts["slug"], facts["slug"].replace("-", " "), facts["title"].lower()}
    out = []
    for line in CLINICAL_REVIEW.read_text().splitlines():
        m = re.match(r"^#{2,3} (.+)$", line)
        if m and any(w in m.group(1).lower() for w in words):
            out.append(m.group(1))
    return out


# --------------------------------------------------------------------------- render
def render(facts: dict, tests: pathlib.Path | None) -> str:
    slug, title, kind = facts["slug"], facts["title"], facts["kind"]
    c1, c1_when, c1_note = check_endpoints(facts)
    c3, c3_when, c3_note, c3_rows = corpus_rows(slug)
    c5, c5_when, c5_note = check_depth(facts, tests)
    c7, c7_when, c7_note = check_docs(facts)

    def when(s: str) -> str:
        return s[:10] if s else "—"

    rows = [
        (1, "Every `live` capability answers on its advertised port", c1, when(c1_when), c1_note),
        (2, MANUAL[2][0], UNANSWERED, "—", MANUAL[2][1]),
        (3, "Retrieval has evidence behind it, or says it does not", c3, when(c3_when), c3_note),
        (4, MANUAL[4][0], UNANSWERED, "—", MANUAL[4][1]),
        (5, "Test depth is proportionate to clinical output", c5, when(c5_when), c5_note),
        (6, MANUAL[6][0], UNANSWERED, "—", MANUAL[6][1]),
        (7, "Every documented command exists and runs", c7, when(c7_when), c7_note),
    ]
    answered = sum(1 for r in rows if r[2] != UNANSWERED)

    L = [f"# {title} — maturity statement", ""]
    L += [f"> Generated by `scripts/gen_maturity.py` at commit `{commit()}` on "
          f"{datetime.date.today().isoformat()}. Every figure below is a measurement; "
          f"section 7 is how you reproduce it.", "", DISCLAIMER, ""]
    L += [f"**{answered} of 7 checks answered.** A check with no answer is reported as "
          f"*not answered*, never as a pass.", ""]

    L += ["## 1. What this subject does", ""]
    if facts["capabilities"]:
        for c in facts["capabilities"]:
            L += [f"- **`{c['id']}`** — {c.get('description', '').strip()}", ""]
    else:
        L += ["This subject registers no capability. It is library or pipeline code.", ""]

    L += ["## 2. Status and maturity", "",
          "What the registry *declares*. Whether it holds is checks 1 and 6.", "",
          "| Capability | Status | Maturity tier | Endpoint |", "| --- | --- | --- | --- |"]
    if facts["capabilities"]:
        for c in facts["capabilities"]:
            L.append(f"| `{c['id']}` | {c.get('status', '—')} | "
                     f"{c.get('maturity', '*none declared*')} | "
                     f"{c.get('endpoint') or '—'} |")
    else:
        L.append("| *none registered* | — | — | — |")
    L += ["", f"Source tree: `{facts['path']}` · {facts['py_files']} Python files · "
              f"{facts['loc']:,} lines including tests · {facts['test_files']} test files. "
              f"Check 5 below counts non-test source only, so the two differ by design.", ""]

    L += ["## 3. The seven checks", "",
          "| # | Check | Result | Measured | Evidence |", "| --- | --- | --- | --- | --- |"]
    for n, name, verdict, date, note in rows:
        L.append(f"| {n} | {name} | **{verdict}** | {date} | {note} |")
    L += ["", "Checks 2, 4 and 6 cannot be answered by a script. They are answered by a "
              "person reading the code, and stay *not answered* until one has.", ""]

    L += ["## 4. What it is grounded in", ""]
    if c3_rows:
        L += ["Retrieval coverage per workflow, from the committed corpus snapshot. Coverage "
              "is the share of declared evidence weight that has any rows behind it.", "",
              "| Workflow | Coverage | Collections with rows |", "| --- | --- | --- |"]
        for workflow, cov, filled, total in sorted(c3_rows, key=lambda r: r[1]):
            L.append(f"| {workflow} | {cov:.0%} | {filled} of {total} |")
        L.append("")
    else:
        L += [c3_note, "",
              "A search of an empty collection succeeds and returns zero rows. Nothing in the "
              "retrieval path fails, so this number is the only thing that would notice.", ""]

    L += ["## 5. Known limitations", "",
          "<!-- REVIEW: capability boundaries, stated as boundaries rather than apologies. -->",
          "*Not yet reviewed.* This section is written during check 2 and check 6, and must "
          "state what this subject cannot do. See `CHANGELOG.md` for the platform-level list.",
          ""]

    L += ["## 6. Open clinical questions", ""]
    secs = clinical_sections(facts)
    if secs:
        L += [f"`CLINICAL_REVIEW.md` sections naming this subject ({len(secs)} matched by "
              "name — a keyword match, not an authoritative index; building that index is "
              "part of check 4):", ""]
        L += [f"- {s}" for s in secs] + [""]
    else:
        L += ["No `CLINICAL_REVIEW.md` section names this subject. That is not the same as "
              "having no unreviewed clinical values — check 4 is what establishes that, and "
              "it is not yet answered.", ""]

    L += ["## 7. How to verify this yourself", "",
          "Every number above comes from one of these. Run them from the repository root.", "",
          "```bash",
          "# check 1 — do the live capabilities answer? (needs the running instance)",
          ".venv/bin/python scripts/probe_live_endpoints.py",
          "",
          "# check 3 — retrieval coverage against the live corpus",
          ".venv/bin/python scripts/check_corpus_weights.py",
          "",
          "# check 5 — tests, and depth against the clinical floor",
          ".venv/bin/python scripts/run_all_tests.py --json /tmp/t.json",
          ".venv/bin/python scripts/check_test_depth.py /tmp/t.json",
          "",
          "# check 7 — every documented command exists, every internal link resolves",
          ".venv/bin/python scripts/check_docs.py --all",
          "",
          "# this page",
          f".venv/bin/python scripts/gen_maturity.py {slug} --tests /tmp/t.json",
          "```", ""]
    L += [f"Checks 2, 4 and 6 have no command. They are read, by a person, against "
          f"`{facts['path']}`.", ""]
    return "\n".join(L)


INDEX = ROOT / "docs" / "honesty" / "claim-integrity.md"
REPO_URL = "https://github.com/ajones1923/hcls-ai-factory/blob/main"


def render_index(rendered: list[tuple[dict, str]]) -> str:
    """One table over all 17, so a gap reads as a gap rather than as an absence of prose."""
    L = ["# Claim Integrity Review", "",
         "> Generated by `scripts/gen_maturity.py` at commit "
         f"`{commit()}` on {datetime.date.today().isoformat()}.", "",
         "The [Capability Maturity Matrix](maturity-matrix.md) shows what the registry "
         "**declares**. This shows what was **measured**, against seven checks applied "
         "identically to all 17 subjects. Where the two disagree, this page is the one with "
         "evidence behind it.", "",
         "**A check with no answer is shown as not answered, never as a pass.** Checks 2, 4 "
         "and 6 require a person reading the code and are unanswered until one has.", "",
         "| Subject | Answered | Pass | Fail | n/a | Statement |",
         "| --- | --- | --- | --- | --- | --- |"]
    tot_pass = tot_fail = tot_na = 0
    for facts, body in rendered:
        npass = body.count("**PASS**")
        nfail = body.count("**FAIL**")
        # n/a is an answer: the check was applied and does not apply to this subject. It is
        # counted separately so that pass + fail + n/a reconciles with "answered" on sight —
        # a column a reader has to derive by subtraction is a column that gets derived wrong.
        nna = len(re.findall(r"\| \*\*n/a\*\* \|", body))
        ans = npass + nfail + nna
        tot_pass += npass
        tot_fail += nfail
        tot_na += nna
        rel = f"{BASE[facts['kind']]}/{facts['slug']}/MATURITY.md"
        L.append(f"| {facts['title']} ({facts['kind']}) | {ans} of 7 | {npass} | {nfail} | "
                 f"{nna} | [MATURITY.md]({REPO_URL}/{rel}) |")
    n = len(rendered)
    tot_ans = tot_pass + tot_fail + tot_na
    L += ["", f"Across {n} subjects: **{tot_ans} of {n * 7} checks answered** — {tot_pass} "
              f"pass, {tot_fail} fail, {tot_na} not applicable. The {n * 7 - tot_ans} "
              f"unanswered are the review work itself, not a clean result.", "",
          "*n/a* means the check was applied and does not apply to that subject — a library "
          "surface with no served endpoint, or a subject that emits no clinical output. It "
          "is not a pass and not a gap.", "",
          "## What each check asks", "",
          "| # | Check | Answered by |", "| --- | --- | --- |",
          "| 1 | Every `live` capability answers on its advertised port | "
          "`scripts/probe_live_endpoints.py` |",
          f"| 2 | {MANUAL[2][0]} | a person |",
          "| 3 | Retrieval has evidence behind it, or says it does not | "
          "`scripts/check_corpus_weights.py` |",
          f"| 4 | {MANUAL[4][0]} | a person |",
          "| 5 | Test depth is proportionate to clinical output | "
          "`scripts/check_test_depth.py` |",
          f"| 6 | {MANUAL[6][0]} | a person |",
          "| 7 | Every documented command exists and runs | `scripts/check_docs.py --all` |",
          "",
          "Checks 1, 3, 5 and 7 are blocking CI gates, so a pass here cannot drift back "
          "unnoticed. Checks 2, 4 and 6 cannot be automated: a generator that printed "
          "*pass* for what it had not checked would be the exact defect this review exists "
          "to find.", ""]
    return "\n".join(L)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("subject", nargs="?", help="one slug (default: all 17 + the index)")
    ap.add_argument("--tests", type=pathlib.Path,
                    help="run_all_tests.py --json output, for check 5")
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any generated file differs from the tree "
                         "(local use: the output carries a generation date and commit, so "
                         "this is not a CI gate)")
    a = ap.parse_args()

    # build() walks the whole roster and owns the per-kind index that DEMO_KEY needs; the
    # cost of rebuilding all 17 is a repo scan, and reimplementing that counter here is how
    # two sources of the same fact start to disagree.
    everything = build()
    wanted = [f for f in everything if a.subject in (None, f["slug"])]
    if not wanted:
        print(f"  unknown subject '{a.subject}'. One of: "
              f"{', '.join(sorted({r[1] for r in ROSTER}))}")
        return 1

    stale, rendered = [], []
    for facts in wanted:
        kind, name = facts["kind"], facts["slug"]
        out = ROOT / BASE[kind] / name / "MATURITY.md"
        body = render(facts, a.tests)
        rendered.append((facts, body))
        if a.check:
            if not out.exists() or out.read_text() != body:
                stale.append(str(out.relative_to(ROOT)))
        else:
            out.write_text(body)
            print(f"  wrote {out.relative_to(ROOT)}")

    # The index only means anything over the whole roster; a single-subject run leaves the
    # existing one alone rather than rewriting it with one row.
    if a.subject is None:
        index = render_index(rendered)
        if a.check:
            if not INDEX.exists() or INDEX.read_text() != index:
                stale.append(str(INDEX.relative_to(ROOT)))
        else:
            INDEX.parent.mkdir(parents=True, exist_ok=True)
            INDEX.write_text(index)
            print(f"  wrote {INDEX.relative_to(ROOT)}")

    if a.check:
        if stale:
            print(f"  FAILURE  {len(stale)} maturity statement(s) are stale:")
            for s in stale:
                print(f"             {s}")
            print("           Regenerate: .venv/bin/python scripts/gen_maturity.py")
            return 1
        print(f"  OK — {len(wanted)} maturity statement(s) match the measurements.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
