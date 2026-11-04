#!/usr/bin/env python3
"""Sweep all 17 subjects for the defect classes the first five reviews actually found.

Reading five subjects produced a small, repeating taxonomy rather than seventeen unique
problems. This looks for those classes everywhere at once, so the worst remaining instance
surfaces now instead of in review order.

  A. A value presented as derived that is invented.
     TSC's trajectory prior reported `method: "population-prior shrinkage (NHD
     mixed-effects)"` while being fitted to 24 simulated patients generated from a hardcoded
     0.026 cm/month. Same shape as the GPU telemetry built from `hash(str(time.time()))`.

  B. A computed honesty verdict that reaches only the log.
     precision-oncology graded its retrieval `sufficient` / `partial` / `insufficient`, used
     it to decide whether to retry, then discarded it and synthesized anyway. Earlier,
     therapeutic-discovery's `assert_real` helpers sat with no caller.

  C. A declared maturity tier that contradicts the entry's own description, in either
     direction. Gated already, both ways, in lib/hcls_common/tests.

  D. A simulated number in a field named as if it were a measurement.
     clinical-imaging's mock VLM fills `confidence`; therapeutic-discovery renames its to
     `simulated_score` with a `score_type` no normalizer accepts.

WHAT THIS IS, AND IS NOT
------------------------
**A triage scanner, not a gate.** It is tuned for recall: it flags candidates to read, and
most flags will be fine on reading. It cannot prove anything is clean, and a subject with
zero flags has NOT passed check 2 -- it has not been looked at. Class C is the only one
decidable by machine, and it is already enforced as a test; it appears here for the report.

That distinction is the whole reason this file is not wired into CI. A scanner that said
PASS for what it had merely failed to pattern-match would be the defect this program exists
to find, one level up.

Usage:
    .venv/bin/python scripts/sweep_claim_classes.py              # the table
    .venv/bin/python scripts/sweep_claim_classes.py --detail     # every hit, with context
    .venv/bin/python scripts/sweep_claim_classes.py cardiology   # one subject
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from subject_facts import BASE, ROSTER  # noqa: E402

REGISTRY = ROOT / "lib" / "hcls_common" / "capabilities.json"
REVIEWED = ROOT / "docs" / "build" / "claim_review.json"

SKIP_PARTS = {".venv", "venv", "__pycache__", "node_modules", "site-packages",
              "vendor_rfdiffusion", "vendor_proteinmpnn", ".git", ".pytest_cache"}

#: Class A. `time.time()` is deliberately absent: in precision-oncology every occurrence was
#: a genuine elapsed-time measurement, and including it produced noise that buried the real
#: hits. A fabricated *duration* is a lesser problem than a fabricated clinical value.
_RANDOM_CALL = re.compile(
    r"\b(?:random\.(?:random|uniform|gauss|normal|randint|choice|sample|shuffle)"
    r"|np\.random\.\w+|numpy\.random\.\w+)\s*\(")

#: Class A, second form: a value keyed off a clock or an object hash rather than read.
_CLOCK_HASH = re.compile(r"hash\s*\(\s*str\s*\(\s*time\.time\s*\(\s*\)")

#: Class B. Status words a subject uses to grade its own output. If one of these is produced
#: in src/ and never appears on the response surface, the grade probably ends at the log.
_VERDICT_WORDS = ("insufficient", "unverified", "canned", "simulated", "degraded",
                  "not_grounded", "ungrounded", "no_evidence", "untestable")

#: Hits read and found to be fine. Kept visible rather than filtered out: a scanner that
#: silently drops what someone once dismissed stops being evidence of anything. Each entry
#: names the reason, and the list may shrink whenever the code changes under it.
REVIEWED_FALSE_POSITIVES = {
    ("precision-intelligence", "B"):
        "the verdict IS surfaced -- `_verify_answer` returns a dict attached as "
        "result['verify'], and a blocking honesty violation appends a research-use "
        "correction to the answer text. Flagged only because the literal word "
        '"unverified" never appears under api/ or app/.',
}

#: Class D. Field names that read as measurements. A simulated value in one of these, in a
#: file that also knows about mocking, is the imaging pattern.
_NEUTRAL_FIELDS = ("confidence", "score", "probability", "accuracy", "similarity",
                   "time_ms", "duration", "latency")

#: Where a subject's answers reach a caller. Class B asks whether a verdict ever gets here.
_SURFACE_DIRS = ("api", "app", "routes")
_SURFACE_FILES = ("models.py", "schemas.py", "main.py", "server.py")


def _py_files(base: pathlib.Path):
    if not base.is_dir():
        return []
    return [p for p in base.rglob("*.py")
            if not (set(p.parts) & SKIP_PARTS)
            and "test" not in p.name
            and "tests" not in p.parts]


def _is_comment(line: str) -> bool:
    return line.lstrip().startswith("#")


def _code_lines(text: str) -> dict[int, str]:
    """Line number -> source, with comments and string literals blanked out.

    Regex over raw source matched its own documentation: `genomic-foundation` was flagged for
    `hash(str(time.time()))` appearing in the COMMENT recording that the fabrication had been
    removed, and `clinical-imaging` for `np.random` inside quoted strings in a demo-guide
    generator. A scanner that reads prose as code reports defects that are descriptions of
    fixed defects, which is a particularly bad failure for this one.

    tokenize gives the real answer. A file that will not tokenize (syntax error, odd
    encoding) falls back to the raw lines: over-reporting is the safer direction here.
    """
    import io
    import tokenize

    lines = text.splitlines()
    out = {i: ln for i, ln in enumerate(lines, 1)}
    try:
        toks = list(tokenize.generate_tokens(io.StringIO(text).readline))
    except (tokenize.TokenError, SyntaxError, IndentationError):
        return out
    for tok in toks:
        if tok.type not in (tokenize.STRING, tokenize.COMMENT):
            continue
        for n in range(tok.start[0], tok.end[0] + 1):
            if n in out:
                # Blank the whole line: a call split across a string boundary is rare, and
                # losing it is better than matching the string's contents.
                out[n] = ""
    return out


def _mock_aware(text: str) -> bool:
    """Does this file know it might be producing simulated data?"""
    return any(k in text for k in ("is_mock", "_provenance", "MockClient", "mock_mode",
                                   "CANNED", "simulated", "Mock"))


def class_a(files) -> list[tuple[str, int, str]]:
    """Randomness or a clock producing a value, in code that is not evidently a mock."""
    out = []
    for p in files:
        try:
            text = p.read_text(errors="ignore")
        except OSError:
            continue
        aware = _mock_aware(text)
        # A seeded generator building a declared synthetic cohort is intentional and
        # labelled; still listed, but marked, so a reader can skip it quickly.
        synthetic = any(w in str(p).lower() for w in ("cohort", "synthetic", "generate_", "demo"))
        for n, line in _code_lines(text).items():
            if not line:
                continue
            if _RANDOM_CALL.search(line) or _CLOCK_HASH.search(line):
                tag = ("mock-aware file" if aware else "NOT mock-aware")
                if synthetic:
                    tag += ", synthetic-data path"
                out.append((str(p.relative_to(ROOT)), n, f"{line.strip()[:90]}  [{tag}]"))
    return out


def class_b(subject_dir: pathlib.Path, files) -> list[tuple[str, int, str]]:
    """A verdict vocabulary that exists in src/ and never reaches the response surface.

    The surface text has the DISCARDING PATTERN stripped out of it first. cart was a missed
    detection until 2026-11-04 for exactly that reason: `app/cart_ui.py` contained a
    copy-pasted `if quality == "insufficient" and plan.sub_questions:` -- the bug itself,
    duplicated into the UI -- and the word appearing there was counted as the verdict having
    been surfaced. A scanner that reads a second copy of the defect as the fix is worse than
    no scanner, so the comparison text excludes lines that are themselves the pattern.
    """
    surface = ""
    for p in _py_files(subject_dir):
        rel = p.relative_to(subject_dir)
        if any(d in rel.parts for d in _SURFACE_DIRS) or p.name in _SURFACE_FILES:
            try:
                for line in p.read_text(errors="ignore").splitlines():
                    stripped = line.strip()
                    # A comparison against the verdict is control flow, not a disclosure.
                    if re.match(r"(el)?if\s+\w+\s*==\s*[\"'](%s)[\"']"
                                % "|".join(_VERDICT_WORDS), stripped):
                        continue
                    surface += line + "\n"
            except OSError:
                pass
    out = []
    for p in files:
        rel = p.relative_to(subject_dir)
        if any(d in rel.parts for d in _SURFACE_DIRS) or p.name in _SURFACE_FILES:
            continue                                    # the surface itself
        try:
            text = p.read_text(errors="ignore")
        except OSError:
            continue
        # NOT blanked: class B is looking for a returned string LITERAL, so the literals
        # are the signal rather than the noise. Comments are still skipped.
        for n, line in enumerate(text.splitlines(), 1):
            if _is_comment(line) or "return" not in line:
                continue
            for w in _VERDICT_WORDS:
                if f'"{w}"' in line or f"'{w}'" in line:
                    if w not in surface:
                        out.append((str(p.relative_to(ROOT)), n,
                                    f'returns "{w}" — that word appears nowhere on the '
                                    f"API/app surface"))
                    break
    return out


def class_d(files) -> list[tuple[str, int, str]]:
    """A randomly produced value assigned to a measurement-sounding field."""
    out = []
    for p in files:
        try:
            text = p.read_text(errors="ignore")
        except OSError:
            continue
        if not _mock_aware(text):
            continue
        for n, line in _code_lines(text).items():
            if not line or not _RANDOM_CALL.search(line):
                continue
            lhs = line.split("=")[0] if "=" in line else line
            hit = next((f for f in _NEUTRAL_FIELDS if f in lhs.lower()), None)
            if hit:
                out.append((str(p.relative_to(ROOT)), n,
                            f"simulated value assigned to `{hit}` — a measurement-sounding "
                            f"name; compare therapeutic-discovery's `simulated_score`"))
    return out


def class_c() -> list[tuple[str, str]]:
    """Tier vs description, both directions. Decidable; also enforced as a test."""
    caps = json.loads(REGISTRY.read_text())["capabilities"]
    contradicts = ("before clinical use", "not yet validated", "pending validation")
    caveats = ("not yet live", "placeholder", "preclinical", "research-use",
               "research/augmentation", "never diagnostic", "superseded")
    out = []
    for c in caps:
        d = (c.get("description") or "").lower()
        if c.get("maturity") == "verified" and any(p in d for p in contradicts):
            out.append((c["id"], "declares `verified` and asks to be validated first"))
        if c.get("status") == "live" and not c.get("maturity") and any(p in d for p in caveats):
            out.append((c["id"], "live, caveats itself in prose, declares no tier"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("subject", nargs="?")
    ap.add_argument("--detail", action="store_true", help="print every hit")
    a = ap.parse_args()

    reviewed = set(json.loads(REVIEWED.read_text()).keys()) if REVIEWED.exists() else set()
    rows, details = [], []
    for kind, slug, title in ROSTER:
        if a.subject and a.subject != slug:
            continue
        d = ROOT / BASE[kind] / slug
        files = _py_files(d)
        A, B, D = class_a(files), class_b(d, files), class_d(files)
        rows.append((title, slug, len(files), len(A), len(B), len(D),
                     "reviewed" if slug in reviewed else ""))
        for name, hits in (("A", A), ("B", B), ("D", D)):
            for f, n, why in hits:
                details.append((slug, name, f, n, why))

    print("  Claim-class sweep — CANDIDATES TO READ, not verdicts.")
    print("  A zero here means nothing was pattern-matched; it does NOT mean the subject is")
    print("  clean, and it is not an answer to check 2.\n")
    print(f"  {'subject':<26} {'files':>5} {'A':>4} {'B':>4} {'D':>4}   state")
    print(f"  {'-'*26} {'-'*5} {'-'*4} {'-'*4} {'-'*4}   {'-'*8}")
    for title, slug, nf, na, nb, nd, state in rows:
        mark = "  <<" if (na or nb or nd) and not state else ""
        print(f"  {slug:<26} {nf:>5} {na:>4} {nb:>4} {nd:>4}   {state}{mark}")
    print(f"\n  A  randomness or a clock producing a value ({sum(r[3] for r in rows)} hits)")
    print(f"  B  a verdict that never reaches the response surface ({sum(r[4] for r in rows)})")
    print(f"  D  a simulated value in a measurement-named field ({sum(r[5] for r in rows)})")

    if REVIEWED_FALSE_POSITIVES:
        print(f"\n  Read and dismissed ({len(REVIEWED_FALSE_POSITIVES)}), still counted "
              f"above so the numbers stay honest:")
        for (slug, cls), why in sorted(REVIEWED_FALSE_POSITIVES.items()):
            print(f"       [{cls}] {slug}: {why}")

    c = class_c()
    print(f"\n  C  tier contradicts the entry's own description: {len(c)}")
    for cid, why in c:
        print(f"       {cid}: {why}")

    if a.detail:
        print("\n  --- detail ---")
        for slug, cls, f, n, why in details:
            print(f"  [{cls}] {slug}  {f}:{n}\n        {why}")
    else:
        print("\n  Re-run with --detail for every hit.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
