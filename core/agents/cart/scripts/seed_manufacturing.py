#!/usr/bin/env python3
"""Seed the cart_manufacturing collection with curated manufacturing/CMC data.

Loads manufacturing process records from data/reference/manufacturing_seed_data.json
and ingests them into the cart_manufacturing Milvus collection.

Data covers the full manufacturing lifecycle:
  - Leukapheresis and starting material quality
  - Viral vector production (lentiviral, retroviral, transposon, mRNA)
  - Transduction efficiency and VCN
  - T-cell expansion protocols (IL-2, IL-7/IL-15, rapid)
  - Harvest, formulation, and dosing
  - Cryopreservation and shipping logistics
  - Release testing (sterility, potency, identity, RCL)
  - Manufacturing failure modes and cost analysis
  - Emerging platforms (POC, allogeneic, non-viral)

Usage:
    python3 scripts/seed_manufacturing.py

Author: Adam Jones
Date: February 2026
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer
from src.collections import CARTCollectionManager
from src.ingest.manufacturing_parser import ManufacturingIngestPipeline


class SimpleEmbedder:
    """Wrapper around SentenceTransformer with encode() API."""
    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "manufacturing_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("CAR-T Manufacturing Data Seeder")
    print("=" * 60)

    # Connect to Milvus
    print("\n[1/3] Connecting to Milvus...")
    manager = CARTCollectionManager()
    manager.connect()

    stats = manager.get_collection_stats()
    existing = stats.get("cart_manufacturing", 0)
    print(f"  cart_manufacturing currently has {existing} records")

    if existing > 0:
        print(f"  WARNING: Collection already has {existing} records.")
        print("  Proceeding — duplicates may be inserted.")

    # Load embedder
    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    # Run ingest pipeline
    print("\n[3/3] Ingesting manufacturing seed data...")
    pipeline = ManufacturingIngestPipeline(
        collection_manager=manager,
        embedder=embedder,
    )
    count = pipeline.run(data_file=str(seed_file))

    # Verify
    stats = manager.get_collection_stats()
    final = stats.get("cart_manufacturing", 0)

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} manufacturing records")
    print(f"  cart_manufacturing now has {final} records")
    print(f"  Seed file: {seed_file}")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
