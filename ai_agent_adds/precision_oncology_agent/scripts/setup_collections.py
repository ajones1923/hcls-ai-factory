#!/usr/bin/env python3
"""Create all Oncology Intelligence Milvus collections and optionally seed data.

Usage:
    python scripts/setup_collections.py [--drop-existing] [--seed]

Options:
    --drop-existing    Drop and recreate all collections
    --seed             Seed all collections with reference data after creation
    --host HOST        Milvus host (default: localhost)
    --port PORT        Milvus port (default: 19530)
"""

import argparse
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.collections import COLLECTION_SCHEMAS, OncoCollectionManager


SEED_SCRIPTS = [
    "seed_variants.py",
    "seed_biomarkers.py",
    "seed_therapies.py",
    "seed_pathways.py",
    "seed_guidelines.py",
    "seed_trials.py",
    "seed_resistance.py",
    "seed_outcomes.py",
    "seed_literature.py",
    "seed_cases.py",
    "seed_knowledge.py",
]


def main():
    parser = argparse.ArgumentParser(
        description="Setup Oncology Intelligence Milvus collections"
    )
    parser.add_argument(
        "--drop-existing",
        action="store_true",
        help="Drop and recreate all collections",
    )
    parser.add_argument(
        "--seed",
        action="store_true",
        help="Seed all collections with reference data after creation",
    )
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    args = parser.parse_args()

    print("=" * 60)
    print("Oncology Intelligence -- Collection Setup")
    print("=" * 60)

    # ── Connect to Milvus ──────────────────────────────────────────────
    kwargs = {}
    if args.host:
        kwargs["host"] = args.host
    if args.port:
        kwargs["port"] = args.port

    manager = OncoCollectionManager(**kwargs)
    print("\n[1/4] Connecting to Milvus...")
    manager.connect()
    print("  Connected successfully.")

    # ── Drop existing (if requested) ───────────────────────────────────
    if args.drop_existing:
        print("\n[2/4] Dropping existing collections...")
        for name in COLLECTION_SCHEMAS:
            print(f"  Dropping {name}...")
            manager.drop_collection(name)
        print("  All collections dropped.")
    else:
        print("\n[2/4] Skipping drop (use --drop-existing to recreate)")

    # ── Create all collections ─────────────────────────────────────────
    print(f"\n[3/4] Creating {len(COLLECTION_SCHEMAS)} collections...")
    collections = manager.create_all_collections()
    for name in sorted(collections.keys()):
        try:
            stats = manager.get_collection_stats(name)
            count = stats.get("num_entities", 0)
            fields = len(stats.get("fields", []))
            print(f"  {name:25s}  {count:>6,} records  ({fields} fields)")
        except Exception:
            print(f"  {name:25s}  created (no stats available)")

    # ── Seed data (if requested) ───────────────────────────────────────
    if args.seed:
        print("\n[4/4] Seeding collections with reference data...")
        manager.disconnect()  # Disconnect; seed scripts create their own connections

        scripts_dir = Path(__file__).resolve().parent
        for script_name in SEED_SCRIPTS:
            script_path = scripts_dir / script_name
            if not script_path.exists():
                print(f"\n  WARNING: Seed script not found: {script_name}")
                continue

            print(f"\n  Running {script_name}...")
            print(f"  {'─' * 50}")
            result = subprocess.run(
                [sys.executable, str(script_path)],
                capture_output=False,
            )
            if result.returncode != 0:
                print(f"  WARNING: {script_name} exited with code {result.returncode}")

        # Reconnect for final stats
        manager.connect()
    else:
        print("\n[4/4] Skipping seed (use --seed to populate reference data)")

    # ── Final stats ────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Final Collection Stats:")
    print(f"{'=' * 60}")
    total_records = 0
    for name in sorted(COLLECTION_SCHEMAS.keys()):
        try:
            stats = manager.get_collection_stats(name)
            count = stats.get("num_entities", 0)
            total_records += count
            print(f"  {name:25s}  {count:>6,} records")
        except Exception:
            print(f"  {name:25s}  (unavailable)")

    print(f"  {'─' * 40}")
    print(f"  {'TOTAL':25s}  {total_records:>6,} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
