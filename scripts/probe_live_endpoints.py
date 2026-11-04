#!/usr/bin/env python3
"""Probe every `live` capability on the running instance and record what answered.

Check 1 of the claim-integrity rubric: **a capability registered `live` must answer on the port
the registry advertises.** On 2026-11-03 nine did not -- :8126, :8507, :8521, :8523, :8526,
:8531, :8538, :8540 and :8544 were all closed while the registry advertised them as live
endpoints. Nothing in the repository could notice, because nothing in the repository looks at
the running machine.

Two tests already pin the *static* half of this: `test_supervisor_covers_live_ports.py` checks
that every `live` capability is supervised by `health-monitor.sh`, in both directions. Being
supervised is not the same as answering. A service can be listed, started, and crash-looping,
and the static check stays green -- which is exactly the "healthy means liveness, not progress"
failure this platform keeps finding in other forms.

Why a snapshot
--------------
The same division as the corpus ratchet, for the same reason: this measures the live instance,
which is right on the box and impossible in CI, because there is no instance on the runner.

  on the box   scripts/probe_live_endpoints.py --write docs/build/endpoint_snapshot.json
  in CI        scripts/check_endpoints.py --enforce

The registry is read fresh from the source tree in both places, so CI still catches what a code
change can break -- a capability flipped to `live`, an endpoint edited, a new port allocated --
while the snapshot supplies the one fact it cannot observe.

This probe only reads. It opens a TCP connection and, if that succeeds, issues one GET. It
starts nothing, stops nothing, and writes only the snapshot file.
"""
from __future__ import annotations

import argparse
import datetime
import json
import pathlib
import socket
import sys
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "lib" / "hcls_common" / "capabilities.json"
DEFAULT_SNAPSHOT = ROOT / "docs" / "build" / "endpoint_snapshot.json"

#: A port that is open but whose HTTP layer is unhappy is still a service that is up. The
#: distinction this check exists to make is "something is listening" vs "nothing is", so any
#: HTTP response counts, including 404 and 405 -- a POST-only route answering a GET with 405
#: has proved it is there.
CONNECT_TIMEOUT = 2.0
READ_TIMEOUT = 4.0


def live_capabilities(registry: pathlib.Path) -> list[dict]:
    """Every capability the registry calls `live` and gives an endpoint for.

    A `live` capability with no endpoint (honesty-gate, event-bus, workflow-composer ...) is a
    library surface, not a served one, and is out of scope here rather than a failure.
    """
    data = json.loads(registry.read_text())
    return [c for c in data.get("capabilities", [])
            if c.get("status") == "live" and c.get("endpoint")]


def split_endpoint(endpoint: str) -> tuple[str, int] | None:
    """"localhost:8126" -> ("localhost", 8126). Anything else is unprobeable, not failed."""
    host, _, port = str(endpoint).rpartition(":")
    if not host or not port.isdigit():
        return None
    return host, int(port)


def probe(host: str, port: int) -> dict:
    """Is something answering here? Returns a result dict, never raises."""
    try:
        with socket.create_connection((host, port), timeout=CONNECT_TIMEOUT):
            pass
    except OSError as exc:
        return {"reachable": False, "how": "tcp", "detail": str(exc)[:120]}

    # The port is open. Ask for a response so that a listening-but-dead process is told apart
    # from a serving one; any status answers that question.
    for path in ("/health", "/"):
        try:
            req = urllib.request.Request(f"http://{host}:{port}{path}", method="GET")
            with urllib.request.urlopen(req, timeout=READ_TIMEOUT) as r:
                return {"reachable": True, "how": "http", "status": r.status, "path": path}
        except urllib.error.HTTPError as exc:
            # 404/405 still proves a server is there and responding.
            return {"reachable": True, "how": "http", "status": exc.code, "path": path}
        except Exception:
            continue

    # Open socket, no HTTP. True for Milvus on :19530, which speaks gRPC.
    return {"reachable": True, "how": "tcp", "detail": "port open, no HTTP response"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", type=pathlib.Path, nargs="?", const=DEFAULT_SNAPSHOT,
                    help="write the snapshot (default path when the flag is bare)")
    ap.add_argument("--registry", type=pathlib.Path, default=REGISTRY)
    a = ap.parse_args()

    caps = live_capabilities(a.registry)
    if not caps:
        print(f"  FAILURE  no live capabilities with endpoints in {a.registry}")
        return 1

    results, unreachable = {}, []
    for cap in sorted(caps, key=lambda c: c["id"]):
        parts = split_endpoint(cap["endpoint"])
        if parts is None:
            results[cap["id"]] = {"endpoint": cap["endpoint"], "reachable": None,
                                  "detail": "endpoint is not host:port"}
            print(f"  {'?':<4} {cap['id']:<36} {cap['endpoint']}  (not host:port)")
            continue
        host, port = parts
        r = probe(host, port)
        r["endpoint"] = cap["endpoint"]
        results[cap["id"]] = r
        if r["reachable"]:
            mark, extra = "up", f"{r['how']} {r.get('status', '')}".strip()
        else:
            mark, extra = "DOWN", r.get("detail", "")
            unreachable.append(cap["id"])
        print(f"  {mark:<4} {cap['id']:<36} {cap['endpoint']:<22} {extra}")

    print(f"\n  {len(caps) - len(unreachable)} of {len(caps)} live capabilities answered")
    if unreachable:
        print(f"  {len(unreachable)} advertised as live and serving nothing: "
              f"{', '.join(unreachable)}")
        print("  A capability nothing serves is the one defect a visitor hits first.")

    if a.write:
        a.write.parent.mkdir(parents=True, exist_ok=True)
        a.write.write_text(json.dumps({
            "measured_at": datetime.datetime.now(datetime.timezone.utc)
                           .replace(microsecond=0).isoformat(),
            "host": socket.gethostname(),
            "results": results,
        }, indent=1, sort_keys=True) + "\n")
        print(f"\n  wrote {a.write}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
