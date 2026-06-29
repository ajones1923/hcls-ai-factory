#!/usr/bin/env python3
"""Seed the onco_literature collection with curated literature data.

Loads literature_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_chunk, and inserts records into the onco_literature
Milvus collection.

Usage: python3 scripts/seed_literature.py
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


def transform_literature_records(seed_data, embedder):
    """Transform seed JSON records into onco_literature schema with embeddings."""
    texts = [rec.get("text_chunk", rec.get("text_summary", "")) for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        records.append({
            "id": str(rec["id"])[:100],
            "embedding": emb,
            "title": str(rec.get("title", ""))[:500],
            "text_chunk": str(rec.get("text_chunk", rec.get("text_summary", "")))[:3000],
            "source_type": str(rec.get("source_type", "pubmed"))[:20],
            "year": int(rec.get("year", 0)),
            "cancer_type": str(rec.get("cancer_type", ""))[:50],
            "gene": str(rec.get("gene", ""))[:50],
            "variant": str(rec.get("variant", ""))[:100],
            "keywords": str(rec.get("keywords", ""))[:1000] if isinstance(rec.get("keywords"), str) else ", ".join(rec.get("keywords", []))[:1000],
            "journal": str(rec.get("journal", ""))[:200],
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "literature_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Literature Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} literature records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_literature")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_literature")
    print(f"  onco_literature currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting literature seed data...")
    records = transform_literature_records(seed_data, embedder)
    count = manager.insert_batch("onco_literature", records)

    try:
        stats = manager.get_collection_stats("onco_literature")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} literature records")
    print(f"  onco_literature now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
