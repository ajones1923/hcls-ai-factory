#!/usr/bin/env python3
"""Seed the onco_resistance collection with curated resistance mechanism data.

Loads resistance_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_summary, and inserts records into the onco_resistance
Milvus collection.

Usage: python3 scripts/seed_resistance.py
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer

from src.collections import OncoCollectionManager


class SimpleEmbedder:
    """Thin wrapper around SentenceTransformer for pipeline compatibility."""

    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        return self.model.encode(texts).tolist()


def transform_resistance_records(seed_data, embedder):
    """Transform seed JSON records into onco_resistance schema with embeddings."""
    texts = [rec["text_summary"] for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        records.append({
            "id": rec.get("id", ""),
            "embedding": emb,
            "primary_therapy": rec.get("primary_therapy", ""),
            "gene": rec.get("gene", ""),
            "mechanism": rec.get("mechanism", ""),
            "bypass_pathway": rec.get("bypass_pathway", ""),
            "alternative_therapies": rec.get("alternative_therapies", ""),
            "text_summary": rec.get("text_summary", ""),
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "resistance_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Resistance Mechanism Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} resistance records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_resistance")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_resistance")
    print(f"  onco_resistance currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting resistance seed data...")
    records = transform_resistance_records(seed_data, embedder)
    count = manager.insert_batch("onco_resistance", records)

    try:
        stats = manager.get_collection_stats("onco_resistance")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} resistance records")
    print(f"  onco_resistance now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
