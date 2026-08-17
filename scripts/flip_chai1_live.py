#!/usr/bin/env python3
"""Flip Chai-1 from `planned` to `live` — but only once it actually answers.

Why this is a script and not a hand edit:

`validate_registry.py` checks that the manifest is well-formed. It does NOT check that a service
responds. So flipping `status` by hand produces a claim that nothing in the repo can catch — on a
public site whose honesty matrix is generated from this file. The registry's own rule ("a live
capability can never be mock-served") is meant to stop exactly that.

So this script probes the endpoint first and refuses to flip if it does not answer.

    python scripts/flip_chai1_live.py --endpoint localhost:8579        # verify, then flip
    python scripts/flip_chai1_live.py --endpoint localhost:8579 --dry-run
    python scripts/flip_chai1_live.py --force                          # flip WITHOUT probing

`--force` exists for the case where the service runs somewhere this machine cannot reach. It prints
a loud warning, because at that point the claim rests on your word rather than on a check.

After flipping, run:  PYTHONPATH=lib python scripts/validate_registry.py
"""
from __future__ import annotations
import argparse, json, pathlib, sys, urllib.error, urllib.request

REPO = pathlib.Path(__file__).resolve().parents[1]
REG = REPO / "lib" / "hcls_common" / "capabilities.json"
LEDGER = REPO / "docs" / "honesty" / "ledger.md"
CAP_ID = "chai1-structure"


def probe(endpoint: str, timeout: float = 5.0) -> tuple[bool, str]:
    """A live service should answer /health. Returns (ok, detail)."""
    url = endpoint if endpoint.startswith("http") else f"http://{endpoint}"
    for path in ("/health", "/"):
        try:
            with urllib.request.urlopen(url.rstrip("/") + path, timeout=timeout) as r:
                if 200 <= r.status < 400:
                    return True, f"{url}{path} -> {r.status}"
        except urllib.error.HTTPError as e:
            if 200 <= e.code < 500:            # answering at all is the signal
                return True, f"{url}{path} -> {e.code}"
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
    return False, f"{url} did not answer (/health, /) — {last}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--endpoint", help="host:port Chai-1 serves on, e.g. localhost:8579")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--force", action="store_true", help="skip the probe (unverified)")
    a = ap.parse_args()

    reg = json.loads(REG.read_text())
    cap = next((c for c in reg["capabilities"] if c["id"] == CAP_ID), None)
    if cap is None:
        print(f"ERROR: {CAP_ID} not in the registry"); return 1
    if cap["status"] == "live":
        print(f"{CAP_ID} is already live — nothing to do."); return 0

    if a.force:
        print("!! --force: flipping WITHOUT probing. The claim now rests on your word.\n")
    else:
        if not a.endpoint:
            print("ERROR: pass --endpoint host:port (or --force to skip the check)"); return 2
        ok, detail = probe(a.endpoint)
        print(f"probe: {detail}")
        if not ok:
            print("\nREFUSING TO FLIP — the service did not answer.\n"
                  "Chai-1 stays `planned`, which is the honest state until /fold really serves.")
            return 3

    cap["status"] = "live"
    if a.endpoint:
        cap["endpoint"] = a.endpoint
    cap["tags"] = [t for t in cap.get("tags", []) if t != "build-in-progress"]
    cap["description"] = cap["description"].replace(
        "STATUS: registered, BUILD IN PROGRESS", "STATUS: live")

    led = LEDGER.read_text()
    led_new = led.replace("Chai-1 is `planned`", "Chai-1 is `live`")

    if a.dry_run:
        print("\n--dry-run — nothing written. Would change:")
        print(f"  {REG.relative_to(REPO)}: status planned -> live"
              f"{', endpoint ' + a.endpoint if a.endpoint else ''}, drop build-in-progress tag")
        print(f"  {LEDGER.relative_to(REPO)}: {'ledger line updated' if led_new != led else 'NO ledger line matched — check by hand'}")
        return 0

    REG.write_text(json.dumps(reg, indent=2, ensure_ascii=False) + "\n")
    print(f"\nregistry: {CAP_ID} -> live")
    if led_new != led:
        LEDGER.write_text(led_new); print("ledger:   Chai-1 planned -> live")
    else:
        print("ledger:   !! no matching line — update docs/honesty/ledger.md by hand")

    print("\nNow run:  PYTHONPATH=lib python scripts/validate_registry.py")
    print("Chai-2 is untouched: it is gated by a partnership, not by this build.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
