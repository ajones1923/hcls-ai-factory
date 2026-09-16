#!/usr/bin/env python3
"""Test-depth floor for subjects that emit clinical output (decision H-D5).

**Why a raw test count is the wrong instrument.** By absolute count the alarming subject is the
single-cell *engine* with 4 tests — but it is 146 lines of deterministic scanpy glue, and 4 tests
is proportionate. Meanwhile clinical-imaging looks well covered at 1,365 tests and is 65,000
lines, which is thinner than it appears. Normalising by source size inverts the ranking, so that
is what this measures:

    tests per 100 non-test source lines

**Why only clinical subjects.** A thin suite on a data-loading engine costs a broken build. A
thin suite on a subject that generates clinical prose costs a wrong answer that reads exactly
like a right one — the failure mode this platform keeps finding. The floor therefore applies to
subjects that emit clinical output, listed explicitly below rather than inferred, because
inferring it from "has a clinical eval case" is circular: the two worst faults found on
2026-09-16 (precision-biomarker, clinical-imaging) were both in subjects with no eval case.

**Why a ratchet, not just a floor.** Three clinical subjects sit at 1.0 against a median of 2.8.
Making 2.0 blocking today would either turn CI red on merge or force ~600 tests to be written in
a hurry, and tests written to satisfy a number are worth nothing. So both rules run:

  RATCHET (blocking)  no clinical subject may fall below its recorded baseline. Backsliding is
                      caught the day it happens, which is the part that actually decays.
  FLOOR   (blocking for NEW subjects, advisory for the three below it)
                      2.0 per 100 lines — below half the median is not a floor anyone has to
                      argue about. A new clinical subject starts at or above it.

**The baseline must be recorded in the environment that enforces it.** CI collects FEWER tests
than a developer box — 300 cart tests there against 415 here, 665 biomarker against 709 — because
optional and GPU-gated dependencies are absent, so those suites skip or fail to collect. A
baseline captured locally therefore fails instantly in CI, which is exactly what happened the
first time this ran. The committed baseline is the CI one; a local run measures at or above it and
passes, which is the right asymmetry.

Usage:
    .venv/bin/python scripts/run_all_tests.py --json /tmp/t.json
    .venv/bin/python scripts/check_test_depth.py /tmp/t.json           # report
    .venv/bin/python scripts/check_test_depth.py /tmp/t.json --enforce # CI gate
    .venv/bin/python scripts/check_test_depth.py /tmp/t.json --update-baseline
"""
from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
BASELINE = ROOT / "docs" / "build" / "test_depth_baseline.json"
BASE = {"engine": "core/engines", "agent": "core/agents", "program": "core/disease-programs"}

FLOOR = 2.0
# A subject may not drop more than this below its baseline before the ratchet trips. Test counts
# move by a test or two for reasons that are not regressions (a skip, a parametrise change); a
# tenth of a test per 100 lines is noise, a fifth is someone deleting a suite.
TOLERANCE = 0.15

# Subjects that emit CLINICAL output -- prose or a recommendation a clinician could act on.
# Listed explicitly, with the reason, because this is the decision H-D5 asked to have recorded.
CLINICAL: dict[tuple[str, str], str] = {
    ("agent", "cart"): "CAR-T construct and trial-matching prose",
    ("agent", "precision-biomarker"): "biomarker interpretation prose",
    ("agent", "pharmacogenomics"): "CPIC dosing recommendations",
    ("agent", "precision-autoimmune"): "autoimmune assessment prose",
    ("agent", "neurology"): "neurological assessment prose",
    ("agent", "clinical-trial"): "eligibility and trial-matching prose",
    ("agent", "rare-disease-diagnostic"): "differential-diagnosis prose",
    ("agent", "single-cell"): "cell-type interpretation prose",
    ("engine", "precision-oncology"): "molecular tumour board decision support",
    ("engine", "cardiology"): "cardiovascular decision support",
    ("engine", "clinical-imaging"): "radiology findings and CAD-RADS assessment",
    ("program", "tuberous-sclerosis"): "variant curation and therapeutic strategy for a patient",
}
# Deliberately OUT, and why -- a subject is exempt because of what it EMITS, not its size:
#   genomic-foundation      alignment/variant calling; correctness is the caller's, not prose
#   precision-intelligence  annotation plumbing beneath the agents
#   therapeutic-discovery   molecule generation and docking scores, not clinical advice
#   structural-biology      structure prediction
#   single-cell (engine)    deterministic scanpy annotation used BY the agent


def source_loc(kind: str, name: str) -> int:
    """Non-test, non-vendored python lines, counted over TRACKED files only.

    Walking the filesystem counts whatever happens to be lying in the working tree — a scratch
    script, a half-finished module — so the same commit measures differently on two machines and
    the ratchet moves for reasons that have nothing to do with tests. `git ls-files` is the same
    everywhere.
    """
    rel = f"{BASE[kind]}/{name}"
    out = subprocess.run(["git", "ls-files", "-z", f"{rel}/*.py"], cwd=ROOT,
                         capture_output=True, text=True).stdout.split("\0")
    total = 0
    for f in out:
        if not f or any(x in f for x in ("/venv/", "/tests/", "vendor_", "__pycache__")):
            continue
        fp = ROOT / f
        if fp.exists():
            total += sum(1 for _ in fp.open(errors="ignore"))
    return total


def measure(results_path: pathlib.Path) -> list[dict]:
    rows = json.loads(results_path.read_text())
    out = []
    for r in rows:
        loc = source_loc(r["kind"], r["name"])
        key = (r["kind"], r["name"])
        out.append({
            "kind": r["kind"], "name": r["name"], "tests": r["passed"], "loc": loc,
            "per100": round(r["passed"] / loc * 100, 2) if loc else 0.0,
            "clinical": key in CLINICAL,
            "emits": CLINICAL.get(key, ""),
        })
    out.sort(key=lambda x: x["per100"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="JSON from scripts/run_all_tests.py --json")
    ap.add_argument("--enforce", action="store_true", help="exit non-zero on a violation")
    ap.add_argument("--update-baseline", action="store_true")
    a = ap.parse_args()

    rows = measure(pathlib.Path(a.results))
    baseline = json.loads(BASELINE.read_text())["subjects"] if BASELINE.exists() else {}

    if a.update_baseline:
        BASELINE.parent.mkdir(parents=True, exist_ok=True)
        BASELINE.write_text(json.dumps({
            "_comment": ("Test-depth ratchet baseline (decision H-D5). Regenerate deliberately "
                         "with scripts/check_test_depth.py --update-baseline; a drop below these "
                         "numbers is a regression, not a new normal."),
            "floor_per_100_loc": FLOOR,
            "tolerance": TOLERANCE,
            "subjects": {f"{r['kind']}/{r['name']}": r["per100"] for r in rows if r["clinical"]},
        }, indent=1) + "\n")
        print(f"  wrote {BASELINE.relative_to(ROOT)}")
        return 0

    violations, below_floor = [], []
    print(f"  {'subject':26s}{'kind':8s}{'tests':>6s}{'LOC':>8s}{'per100':>8s}  emits")
    for r in rows:
        key = f"{r['kind']}/{r['name']}"
        mark = ""
        if r["clinical"]:
            base = baseline.get(key)
            if base is not None and r["per100"] < base - TOLERANCE:
                violations.append((key, r["per100"], base)); mark = "  << REGRESSED"
            elif r["per100"] < FLOOR:
                below_floor.append((key, r["per100"]))
                mark = "  << below floor" + ("" if base is not None else " (NEW SUBJECT)")
        print(f"  {r['name']:26s}{r['kind']:8s}{r['tests']:6d}{r['loc']:8d}"
              f"{r['per100']:8.1f}  {r['emits'][:34]}{mark}")

    clin = sorted(r["per100"] for r in rows if r["clinical"])
    if clin:
        print(f"\n  clinical subjects: n={len(clin)}  min={clin[0]:.1f}  "
              f"median={clin[len(clin)//2]:.1f}  max={clin[-1]:.1f}  floor={FLOOR}")

    new_below = [k for k, _ in below_floor if k not in baseline]
    for k, v, b in violations:
        print(f"  REGRESSION  {k}: {v} is below its baseline {b}")
    for k in new_below:
        print(f"  NEW SUBJECT BELOW FLOOR  {k} must start at or above {FLOOR}")
    if below_floor and not violations and not new_below:
        print("\n  Known-thin (advisory, recorded in HARDENING_WORKBOOK 3.4): "
              + ", ".join(f"{k} {v}" for k, v in below_floor))

    failed = bool(violations or new_below)
    if not failed:
        print("\n  OK — no clinical subject regressed, and no new one starts below the floor.")
    return 1 if (failed and a.enforce) else 0


if __name__ == "__main__":
    raise SystemExit(main())
