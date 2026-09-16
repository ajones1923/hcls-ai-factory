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


def present(term, text_low: str) -> bool:
    """True if `term` appears. A list term means "any of these equivalent spellings".

    Clinical facts have more than one correct rendering — HLA-B*27 and HLA-B27 are the same
    allele, GBA and GBA1 the same gene. An eval that fails on formatting measures formatting.
    """
    if isinstance(term, list):
        return any(str(x).lower() in text_low for x in term)
    return str(term).lower() in text_low


def wait_healthy(port: int, tries: int = 30, gap: int = 2) -> bool:
    """Block until `port` answers /health. Returns False if it never does."""
    import time
    for _ in range(tries):
        try:
            urllib.request.urlopen(f"http://localhost:{port}/health", timeout=3).read()
            return True
        except Exception:
            time.sleep(gap)
    return False


def ask(case: dict, timeout: int = 240) -> tuple[str, str]:
    """-> (answer_text, error). Never raises."""
    body = json.dumps({case.get("field", "question"): case["question"]}).encode()
    url = f"http://localhost:{case['port']}{case['path']}"
    import os
    headers = {"Content-Type": "application/json"}
    if os.getenv("HCLS_API_KEY"):           # gate is fail-closed once the key is set
        headers["X-API-Key"] = os.environ["HCLS_API_KEY"]
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
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

    # Wait for each service to answer /health first. Querying a just-restarted agent returns
    # its stub ("Search completed...") and scores as a MISS -- a false alarm, and false alarms
    # are how a check earns being ignored.
    import time
    for port in sorted({c["port"] for c in cases}):
        wait_healthy(port)

    rows, passed = [], 0
    print(f"{'id':26s}{'agent':26s}{'verdict':9s}detail")
    for c in cases:
        ans, err = ask(c)
        if err or any(not present(t, ans.lower()) for t in c.get("expect", [])):
            # One retry. These are LLM answers, so a single sample is noisy; a fact the agent
            # knows should survive a second ask. A case that fails twice is a real finding.
            #
            # A transport error usually means the supervisor cycled that service mid-run: it
            # restarts anything it finds down every 5 minutes, and a restart plus model load
            # outlasts any fixed sleep. Sleeping 2s and re-asking produced three false ERRORs
            # against pharmacogenomics on 2026-09-16 whose answers were in fact correct. Wait
            # for /health instead -- a false alarm is how a check earns being ignored.
            if err and not err.startswith("HTTP"):
                wait_healthy(c["port"])
            else:
                time.sleep(2)
            ans2, err2 = ask(c)
            if not err2 and len(ans2) > len(ans):
                ans, err = ans2, err2
        low = ans.lower()
        # An expect term may be a LIST of equivalent spellings ("any of"). Clinical facts have
        # more than one correct rendering -- HLA-B*27 and HLA-B27 are the same allele -- and an
        # eval that fails on formatting is measuring formatting.
        missing = [t if isinstance(t, str) else "|".join(t)
                   for t in c.get("expect", []) if not present(t, low)]
        violated = [t if isinstance(t, str) else "|".join(t)
                    for t in c.get("forbid", []) if present(t, low)]
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
