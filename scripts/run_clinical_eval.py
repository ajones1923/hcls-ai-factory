#!/usr/bin/env python3
"""Measure whether the agents' answers are CLINICALLY CORRECT, not merely returned.

The 17 subject suites (8,404 tests) prove the plumbing. They cannot tell you that an agent
said TPMT when the right answer was TPMT, or that it invented an NCT number. This asks a small
set of questions whose correct answers are externally verifiable, and scores the prose.

Scoring is deliberately crude and transparent -- substring presence for `expect`, absence for
`forbid` -- because a subtle scorer on six cases would be measuring the scorer. `forbid` is
what makes a case discriminating: it names the specific wrong answer a model plausibly gives.

    .venv/bin/python scripts/run_clinical_eval.py
    .venv/bin/python scripts/run_clinical_eval.py --json out.json
    .venv/bin/python scripts/run_clinical_eval.py pgx        # substring filter on id/agent
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
SPEC = ROOT / "demo" / "eval" / "clinical_questions.yaml"


def ask(case: dict, timeout: int = 240) -> tuple[str, str]:
    """-> (answer_text, error). Never raises."""
    body = json.dumps({case.get("field", "question"): case["question"]}).encode()
    url = f"http://localhost:{case['port']}{case['path']}"
    req = urllib.request.Request(url, data=body,
                                 headers={"Content-Type": "application/json"}, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            d = json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return "", f"HTTP {e.code}"
    except Exception as e:
        return "", type(e).__name__
    for k in ("answer", "response", "summary", "interpretation"):
        if isinstance(d.get(k), str) and d[k]:
            return d[k], ""
    return "", "no answer field in response"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("filter", nargs="?", default="")
    ap.add_argument("--json", dest="jsonout")
    a = ap.parse_args()

    try:
        import yaml
    except ImportError:
        print("pyyaml required: .venv/bin/pip install pyyaml"); return 1
    cases = yaml.safe_load(SPEC.read_text())["cases"]
    if a.filter:
        cases = [c for c in cases if a.filter in c["id"] or a.filter in c["agent"]]

    rows, passed = [], 0
    print(f"{'id':26s}{'agent':26s}{'verdict':9s}detail")
    for c in cases:
        ans, err = ask(c)
        low = ans.lower()
        missing = [t for t in c.get("expect", []) if t.lower() not in low]
        violated = [t for t in c.get("forbid", []) if t.lower() in low]
        if err:
            verdict, detail = "ERROR", err
        elif violated:
            verdict, detail = "WRONG", "asserted: " + "; ".join(violated)
        elif missing:
            verdict, detail = "MISS", "missing: " + ", ".join(missing)
        else:
            verdict, detail = "PASS", f"{len(ans)} chars"
            passed += 1
        print(f"  {c['id']:24s}{c['agent']:26s}{verdict:9s}{detail[:60]}")
        rows.append({**{k: c[k] for k in ("id", "agent")}, "verdict": verdict,
                     "detail": detail, "answer_chars": len(ans)})

    print(f"\n  {passed}/{len(cases)} clinically correct")
    if a.jsonout:
        json.dump(rows, open(a.jsonout, "w"), indent=1)
        print(f"  wrote {a.jsonout}")
    # A MISS is informative, not fatal; a WRONG or ERROR is a real failure.
    return 1 if any(r["verdict"] in ("WRONG", "ERROR") for r in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
