#!/usr/bin/env python3
"""Corpus-coverage ratchet: a workflow may not retrieve from less evidence than it used to.

v2.2.0 made every weighted subject report `retrieval_coverage` honestly, and the first honest
measurement was that **56 of 70 workflows retrieve from under 75% of the evidence they
declare**. Seeding fixes that once. This is what stops it drifting back.

Why a snapshot instead of a live query
--------------------------------------
`check_corpus_weights.py --enforce` measures the live instance, which is the right thing to run
on the box and impossible to run in CI: there is no Milvus on the runner. A gate that cannot run
where the merge happens is not a gate. So the division is:

  on the box   `check_corpus_weights.py --write-snapshot docs/build/corpus_snapshot.json`
  in CI        this script, recomputing coverage from that snapshot and the weight tables

The weight tables come from the subject source trees, so CI exercises the part that code changes
can break — a boost edited, a collection dropped from a map, a weight renormalized — while the
snapshot supplies the row counts it cannot measure.

**A stale snapshot is a failure, not a pass.** That is the whole reason for `--max-age-days`: a
snapshot frozen at the moment everything was green would otherwise certify a corpus that has
since been emptied, which is precisely the shape of defect this repo spent 2026-11-04 removing —
a check that reports success because it is measuring something that stopped changing.
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "scripts"))

DEFAULT_SNAPSHOT = ROOT / "docs" / "build" / "corpus_snapshot.json"
DEFAULT_BASELINE = ROOT / "docs" / "build" / "corpus_baseline.json"
MAX_AGE_DAYS = 45


def _load(path: pathlib.Path) -> dict:
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError as exc:
        print(f"  FAILURE  {path.name} is not readable JSON: {exc}")
        raise SystemExit(1)


def measure(counts: dict[str, int]) -> tuple[dict[str, float], list[str]]:
    """Per-workflow coverage from recorded counts, keyed "subject/workflow"."""
    from check_corpus_weights import WEIGHT_TABLES, weight_maps

    from hcls_common.corpus_coverage import coverage

    out: dict[str, float] = {}
    notes: list[str] = []
    for subject in sorted(WEIGHT_TABLES):
        maps, err = weight_maps(subject)
        if err:
            notes.append(f"{subject}: {err}")
            continue
        for workflow, weights in maps.items():
            v = coverage(weights, {c: counts.get(c, 0) for c in weights})
            out[f"{subject}/{workflow}"] = round(v.coverage, 4)
    return out, notes


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", type=pathlib.Path, default=DEFAULT_SNAPSHOT)
    ap.add_argument("--baseline", type=pathlib.Path, default=DEFAULT_BASELINE)
    ap.add_argument("--enforce", action="store_true", help="exit non-zero on a regression")
    ap.add_argument("--update-baseline", action="store_true",
                    help="record the current measurement as the new baseline")
    ap.add_argument("--max-age-days", type=int, default=MAX_AGE_DAYS)
    a = ap.parse_args()

    from hcls_common.corpus_coverage import compare_coverage

    snap = _load(a.snapshot)
    if not snap.get("collections"):
        print(f"  FAILURE  no corpus snapshot at {a.snapshot}.")
        print("           Run on the box, then commit the file:")
        print("           .venv/bin/python scripts/check_corpus_weights.py "
              "--write-snapshot docs/build/corpus_snapshot.json")
        return 1 if a.enforce else 0

    counts = {str(k): int(v) for k, v in snap["collections"].items()}
    measured_at = snap.get("measured_at", "")
    age_days = None
    if measured_at:
        try:
            when = datetime.datetime.fromisoformat(measured_at)
            if when.tzinfo is None:
                when = when.replace(tzinfo=datetime.timezone.utc)
            age_days = (datetime.datetime.now(datetime.timezone.utc) - when).days
        except ValueError:
            pass

    current, notes = measure(counts)
    if not current:
        print("  FAILURE  no weight table could be imported, so nothing was measured.")
        for n in notes:
            print(f"           ! {n}")
        return 1 if a.enforce else 0

    if a.update_baseline:
        a.baseline.write_text(json.dumps(
            {"recorded_from": measured_at, "coverage": current}, indent=1, sort_keys=True) + "\n")
        print(f"  baseline updated: {len(current)} workflows -> {a.baseline}")
        return 0

    baseline = _load(a.baseline).get("coverage", {})
    regressions, vanished = compare_coverage(baseline, current)

    total = len(current)
    below75 = sum(1 for v in current.values() if v < 0.75)
    mean = sum(current.values()) / total
    print(f"  {total} workflows measured from a snapshot of {len(counts)} collections "
          f"({sum(counts.values())} rows)")
    print(f"  mean coverage {mean:.0%}; {below75} of {total} below 75%")
    if age_days is not None:
        print(f"  snapshot taken {measured_at} ({age_days}d old)")

    failed = False

    if age_days is not None and age_days > a.max_age_days:
        print(f"\n  FAILURE  the snapshot is {age_days} days old (limit {a.max_age_days}).")
        print("           A frozen snapshot certifies a corpus nobody has looked at. Refresh")
        print("           it on the box and commit the result.")
        failed = True

    if vanished:
        print(f"\n  FAILURE  {len(vanished)} workflow(s) in the baseline are no longer measured:")
        for key in vanished[:10]:
            print(f"             {key}")
        print("           A workflow that stopped reporting looks exactly like one that never")
        print("           regressed. If it was removed on purpose, --update-baseline.")
        failed = True

    if regressions:
        print(f"\n  FAILURE  {len(regressions)} workflow(s) retrieve from less than before:")
        for r in regressions[:12]:
            print(f"             {r}")
        print("           Nothing in the retrieval path fails when a collection goes quiet:")
        print("           a search of an empty collection succeeds and returns zero rows.")
        failed = True

    for n in notes:
        print(f"  ! {n}")

    if failed:
        return 1 if a.enforce else 0

    if not baseline:
        print("\n  OK -- no baseline recorded yet. Run --update-baseline to pin this measurement.")
    else:
        print("\n  OK -- no workflow retrieves from less evidence than the baseline recorded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
