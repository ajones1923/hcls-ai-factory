#!/usr/bin/env python3
"""Endpoint ratchet: a capability registered `live` may not newly stop answering.

Check 1 of the claim-integrity rubric, as a merge gate. It reads the registry from the source
tree and the last probe from `docs/build/endpoint_snapshot.json`, and fails when:

  * a capability is `live` with an endpoint and the snapshot has no record of it -- a new or
    renamed capability has never been probed, and an unprobed endpoint is unknown, not up;
  * a capability was answering in the baseline and is not in the snapshot;
  * the snapshot is older than `--max-age-days`.

Why a ratchet rather than "all live capabilities must answer"
-------------------------------------------------------------
Because today nine of them do not, and a gate that fails on day one gets switched off on day
two. v2.2.0 is released and the box still runs September code, so nine capabilities advertise
endpoints the running instance does not serve. That is recorded in the baseline as a known,
dated debt with the reason -- not hidden, and not allowed to grow.

After the deploy, `--update-baseline` records the repaired state and this becomes what it
should be: every `live` capability answers, and any regression fails the build. The baseline
shrinking is the measure of progress; it must never grow without someone writing down why.

A capability coming back up does NOT fail -- it prints a notice asking for --update-baseline.
Failing a build because something got fixed teaches people to route around the gate.
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "lib" / "hcls_common" / "capabilities.json"
DEFAULT_SNAPSHOT = ROOT / "docs" / "build" / "endpoint_snapshot.json"
DEFAULT_BASELINE = ROOT / "docs" / "build" / "endpoint_baseline.json"
MAX_AGE_DAYS = 30


def _load(path: pathlib.Path) -> dict:
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError as exc:
        print(f"  FAILURE  {path.name} is not readable JSON: {exc}")
        raise SystemExit(1)


def live_ids(registry: pathlib.Path) -> set[str]:
    data = _load(registry)
    return {c["id"] for c in data.get("capabilities", [])
            if c.get("status") == "live" and c.get("endpoint")}


def snapshot_age_days(measured_at: str) -> int | None:
    try:
        when = datetime.datetime.fromisoformat(measured_at)
    except (TypeError, ValueError):
        return None
    if when.tzinfo is None:
        when = when.replace(tzinfo=datetime.timezone.utc)
    return (datetime.datetime.now(datetime.timezone.utc) - when).days


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--snapshot", type=pathlib.Path, default=DEFAULT_SNAPSHOT)
    ap.add_argument("--baseline", type=pathlib.Path, default=DEFAULT_BASELINE)
    ap.add_argument("--registry", type=pathlib.Path, default=REGISTRY)
    ap.add_argument("--enforce", action="store_true", help="exit non-zero on a failure")
    ap.add_argument("--update-baseline", action="store_true",
                    help="record the current measurement as the accepted state")
    ap.add_argument("--max-age-days", type=int, default=MAX_AGE_DAYS)
    ap.add_argument("--reason", default="",
                    help="why the known-down set is accepted, and what closes it")
    a = ap.parse_args()

    snap = _load(a.snapshot)
    # An ABSENT snapshot is the failure. A snapshot that recorded no results is not the same
    # thing -- it is a measurement that found nothing to probe, and collapsing the two would
    # report a real (if empty) probe as a missing one, then pass the unprobed check below
    # without ever running it.
    if not snap:
        print(f"  FAILURE  no endpoint snapshot at {a.snapshot}.")
        print("           Run on the box, then commit the file:")
        print("           .venv/bin/python scripts/probe_live_endpoints.py --write")
        return 1 if a.enforce else 0

    results = snap.get("results", {})
    declared = live_ids(a.registry)
    down = sorted(i for i in declared
                  if i in results and results[i].get("reachable") is False)
    up = sorted(i for i in declared
                if i in results and results[i].get("reachable") is True)
    unprobed = sorted(declared - set(results))

    if a.update_baseline:
        # A reason is required whenever anything is accepted. An outage recorded without one
        # becomes background noise within a release, and then nobody can say whether the list
        # is debt someone chose or debt that accumulated. Carried over when not re-stated.
        reason = a.reason or _load(a.baseline).get("reason", "")
        if down and not reason:
            print("  FAILURE  refusing to accept an outage with no reason recorded.")
            print("           Pass --reason \"...\" saying why these are down and what")
            print("           closes it. An unexplained baseline is indistinguishable from")
            print("           a gate somebody disabled.")
            return 1
        a.baseline.write_text(json.dumps({
            "recorded_from": snap.get("measured_at", ""),
            "note": "capabilities registered live that were NOT answering when recorded",
            "reason": reason,
            "known_down": down,
        }, indent=1, sort_keys=True) + "\n")
        print(f"  baseline updated: {len(down)} known-down -> {a.baseline}")
        return 0

    baseline = _load(a.baseline)
    known_down = set(baseline.get("known_down", []))
    age = snapshot_age_days(snap.get("measured_at", ""))

    print(f"  {len(up)} of {len(declared)} live capabilities answered "
          f"on {snap.get('host', 'an unnamed host')}")
    if age is not None:
        print(f"  probe taken {snap.get('measured_at')} ({age}d old)")
    if known_down:
        print(f"  {len(known_down)} known-down, accepted in the baseline")

    failed = False

    if age is not None and age > a.max_age_days:
        print(f"\n  FAILURE  the probe is {age} days old (limit {a.max_age_days}).")
        print("           A frozen probe certifies an instance nobody has looked at.")
        failed = True

    if unprobed:
        print(f"\n  FAILURE  {len(unprobed)} capability(ies) are registered `live` with an "
              f"endpoint and have never been probed:")
        for i in unprobed:
            print(f"             {i}")
        print("           An unprobed endpoint is unknown, not up. Re-probe on the box.")
        failed = True

    regressed = sorted(set(down) - known_down)
    if regressed:
        print(f"\n  FAILURE  {len(regressed)} capability(ies) advertise `live` and serve "
              f"nothing, and are not in the baseline:")
        for i in regressed:
            print(f"             {i}  {results[i].get('endpoint')}  "
                  f"{results[i].get('detail', '')}")
        print("           Either the service is down, or the registry claims more than the")
        print("           instance delivers. Both are the same thing to a visitor.")
        failed = True

    recovered = sorted(known_down - set(down))
    if recovered:
        print(f"\n  {len(recovered)} baselined capability(ies) now answer: "
              f"{', '.join(recovered)}")
        print("  Run --update-baseline so the debt record shrinks with the debt.")

    if failed:
        return 1 if a.enforce else 0

    if known_down:
        print(f"\n  OK -- no new capability stopped answering. {len(known_down)} still "
              f"advertised and unserved; see docs/build/endpoint_baseline.json.")
    else:
        print("\n  OK -- every capability registered `live` answers on its advertised port.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
