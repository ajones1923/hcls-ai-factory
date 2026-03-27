#!/usr/bin/env python3
"""Seed the PGx knowledge graph into Milvus collections.

Reads all 14 reference JSON files from data/reference/, generates BGE-small-en-v1.5
embeddings for each record, and inserts them into the corresponding Milvus collections.

Assumes collections have already been created by setup_collections.py.

Usage:
    python scripts/seed_knowledge.py
    python scripts/seed_knowledge.py --host milvus-standalone --port 19530
    python scripts/seed_knowledge.py --batch-size 64
    python scripts/seed_knowledge.py --collections pgx_gene_reference pgx_drug_guidelines

Author: Adam Jones
Date: March 2026
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from src.collections import COLLECTION_SCHEMAS, PGxCollectionManager

# Import the seeding utilities from setup_collections
from scripts.setup_collections import (
    REFERENCE_FILE_MAP,
    seed_collection,
)


def main():
    parser = argparse.ArgumentParser(
        description="Seed PGx knowledge base from reference JSON files"
    )
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    parser.add_argument(
        "--batch-size",
        type=int,
        default=32,
        help="Number of records to embed and insert per batch (default: 32)",
    )
    parser.add_argument(
        "--collections",
        nargs="+",
        default=None,
        help="Specific collection(s) to seed (default: all)",
    )
    args = parser.parse_args()

    logger.info("=" * 65)
    logger.info("  PGx Intelligence Agent — Seed Knowledge Base")
    logger.info("=" * 65)

    # ── Load embedding model ──
    try:
        from sentence_transformers import SentenceTransformer

        logger.info("Loading embedding model: BAAI/bge-small-en-v1.5")
        model = SentenceTransformer("BAAI/bge-small-en-v1.5")
        logger.info("Embedding model loaded (384-dim)")
    except ImportError as e:
        logger.error(f"Cannot load embedding model: {e}")
        logger.error("Install: pip install sentence-transformers")
        sys.exit(1)

    # ── Connect to Milvus ──
    manager = PGxCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # ── Determine which collections to seed ──
    reference_dir = PROJECT_ROOT / "data" / "reference"
    if not reference_dir.exists():
        logger.error(f"Reference directory not found: {reference_dir}")
        manager.disconnect()
        sys.exit(1)

    target_map = REFERENCE_FILE_MAP.copy()
    if args.collections:
        target_map = {
            k: v for k, v in target_map.items() if k in args.collections
        }
        if not target_map:
            logger.error(f"None of the specified collections are valid: {args.collections}")
            logger.info(f"Valid collections: {list(REFERENCE_FILE_MAP.keys())}")
            manager.disconnect()
            sys.exit(1)

    # ── Seed each collection ──
    start = time.time()
    results: Dict[str, int] = {}
    files_found = 0
    files_missing = 0

    for collection_name, filename in target_map.items():
        filepath = reference_dir / filename
        if not filepath.exists():
            logger.warning(f"  Skipping {collection_name}: {filename} not found")
            results[collection_name] = 0
            files_missing += 1
            continue

        files_found += 1
        logger.info(f"  Loading {filename} -> {collection_name}")

        with open(filepath, "r") as f:
            records = json.load(f)

        if not isinstance(records, list):
            logger.error(f"  {filename}: expected JSON array, got {type(records).__name__}")
            results[collection_name] = 0
            continue

        logger.info(f"  {filename}: {len(records)} records to embed and insert")
        count = seed_collection(
            manager, collection_name, records, model, batch_size=args.batch_size
        )
        results[collection_name] = count
        logger.info(f"  {collection_name}: {count} records seeded")

    elapsed = time.time() - start

    # ── Summary ──
    logger.info("")
    logger.info("=" * 65)
    logger.info("  Seed Summary")
    logger.info("=" * 65)
    logger.info(f"  Reference files found:   {files_found}")
    logger.info(f"  Reference files missing: {files_missing}")
    logger.info(f"  Time elapsed:            {elapsed:.1f}s")
    logger.info("")

    total_seeded = 0
    for name, count in results.items():
        status = f"{count:>6,} records" if count > 0 else "  (skipped)"
        logger.info(f"  {name:<30s} {status}")
        total_seeded += count

    logger.info("")
    logger.info(f"  Total records seeded: {total_seeded:,}")
    logger.info("=" * 65)

    # ── Final stats ──
    stats = manager.get_collection_stats()
    logger.info("")
    logger.info("Final Milvus collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    manager.disconnect()
    logger.info("Seed complete.")


if __name__ == "__main__":
    main()
