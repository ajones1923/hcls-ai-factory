#!/usr/bin/env python3
"""Seed the imaging_devices Milvus collection from reference JSON data.

Loads device_seed_data.json, validates each record as an ImagingDevice
Pydantic model, generates BGE-small-en-v1.5 embeddings, and inserts
into the imaging_devices collection.

Usage:
    python scripts/seed_devices.py
    python scripts/seed_devices.py --host 10.0.0.1 --port 19530

Author: Adam Jones
Date: February 2026
"""

import argparse
import json
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger
from sentence_transformers import SentenceTransformer
from pymilvus import connections

from config.settings import settings
from src.collections import ImagingCollectionManager
from src.models import ImagingDevice

COLLECTION_NAME = "imaging_devices"
SEED_FILE = "device_seed_data.json"


def main():
    parser = argparse.ArgumentParser(
        description="Seed imaging_devices collection from reference data"
    )
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    parser.add_argument(
        "--batch-size", type=int, default=32, help="Embedding batch size"
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  Imaging Intelligence Agent — Seed Devices")
    print("=" * 65)

    # Connect to Milvus via collection manager
    manager = ImagingCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Load embedding model
    logger.info(f"Loading embedding model: {settings.EMBEDDING_MODEL}")
    embedder = SentenceTransformer(settings.EMBEDDING_MODEL)

    # Load seed data
    seed_file = PROJECT_ROOT / "data" / "reference" / SEED_FILE
    logger.info(f"Loading seed data from {seed_file}")
    with open(seed_file) as f:
        raw_data = json.load(f)
    logger.info(f"Loaded {len(raw_data)} raw records")

    # Parse and validate
    records = []
    for item in raw_data:
        try:
            record = ImagingDevice(**item)
            records.append(record)
        except Exception as e:
            logger.warning(f"Skipping invalid record: {e}")

    logger.info(f"Validated {len(records)} ImagingDevice records")

    if not records:
        logger.warning("No valid records to seed. Exiting.")
        manager.disconnect()
        return

    # Generate embeddings
    logger.info("Generating embeddings...")
    texts = [r.to_embedding_text() for r in records]
    embeddings = embedder.encode(texts, show_progress_bar=True, normalize_embeddings=True).tolist()

    # Build insert data
    insert_data = []
    for rec, emb in zip(records, embeddings):
        row = rec.model_dump()
        for k, v in row.items():
            if hasattr(v, "value"):
                row[k] = v.value
        row["embedding"] = emb
        insert_data.append(row)

    # Insert into collection
    count = manager.insert_batch(COLLECTION_NAME, insert_data, batch_size=args.batch_size)
    logger.info(f"Seeded {count} records into {COLLECTION_NAME}")

    # Show final stats
    stats = manager.get_collection_stats()
    logger.info(f"  {COLLECTION_NAME}: {stats.get(COLLECTION_NAME, 0):,} total records")

    manager.disconnect()

    print()
    print("=" * 65)
    print(f"  Seeded {count} devices into {COLLECTION_NAME}")
    print("=" * 65)


if __name__ == "__main__":
    main()
