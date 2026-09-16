#!/usr/bin/env python3
"""Load every Milvus collection into memory. Idempotent, cheap, safe to run on every tick.

Milvus does NOT load collections automatically after a restart. A seeded, indexed collection
answers `collection not loaded` until something calls load_collection() — which nothing did,
so a reboot silently reduced every agent to zero retrieval while all 32 services still
reported healthy. That is the same invisible-degradation shape as the missing .env, and it is
the most predictable way a reboot breaks this platform.

Only loads what is not already loaded, so the steady-state cost is one state query per
collection.

    .venv/bin/python scripts/load_collections.py           # load, print a summary
    .venv/bin/python scripts/load_collections.py --quiet   # only report problems
"""
from __future__ import annotations

import argparse
import os
import sys


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--uri", default=None)
    a = ap.parse_args()

    uri = a.uri or (f"http://{os.getenv('MILVUS_HOST', 'localhost')}:"
                    f"{os.getenv('MILVUS_PORT', '19530')}")
    try:
        from pymilvus import MilvusClient
        client = MilvusClient(uri=uri)
        names = sorted(client.list_collections())
    except Exception as exc:
        print(f"milvus unreachable at {uri}: {type(exc).__name__}", file=sys.stderr)
        return 1

    loaded = already = failed = 0
    for name in names:
        try:
            state = str(client.get_load_state(collection_name=name).get("state"))
            if state.endswith("Loaded"):
                already += 1
                continue
            client.load_collection(collection_name=name)
            loaded += 1
            if not a.quiet:
                print(f"  loaded {name}")
        except Exception as exc:
            failed += 1
            print(f"  FAILED {name}: {type(exc).__name__}", file=sys.stderr)

    if not a.quiet or loaded or failed:
        rows = sum(_rows(client, n) for n in names)
        print(f"collections: {len(names)} · loaded now {loaded} · already {already} "
              f"· failed {failed} · {rows:,} vectors")
    return 1 if failed else 0


def _rows(client, name: str) -> int:
    try:
        return int(client.get_collection_stats(name).get("row_count", 0))
    except Exception:
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
