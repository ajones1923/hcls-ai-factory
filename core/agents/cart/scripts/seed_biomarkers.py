#!/usr/bin/env python3
"""Seed the cart_biomarkers collection with curated biomarker data.

Usage:
    python3 scripts/seed_biomarkers.py
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer
from src.collections import CARTCollectionManager
from src.ingest.biomarker_parser import BiomarkerIngestPipeline


class SimpleEmbedder:
    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")
    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "biomarker_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("CAR-T Biomarker Data Seeder")
    print("=" * 60)

    print("\n[1/3] Connecting to Milvus...")
    manager = CARTCollectionManager()
    manager.connect()

    stats = manager.get_collection_stats()
    existing = stats.get("cart_biomarkers", 0)
    print(f"  cart_biomarkers currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Ingesting biomarker seed data...")
    pipeline = BiomarkerIngestPipeline(collection_manager=manager, embedder=embedder)
    count = pipeline.run(data_file=str(seed_file))

    stats = manager.get_collection_stats()
    final = stats.get("cart_biomarkers", 0)

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} biomarker records")
    print(f"  cart_biomarkers now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
