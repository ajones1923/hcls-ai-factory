#!/usr/bin/env python3
"""Seed the onco_guidelines collection with curated guideline data.

Loads guideline_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_summary, and inserts records into the onco_guidelines
Milvus collection.

Usage: python3 scripts/seed_guidelines.py
"""

import json
import re
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


def _extract_year(rec):
    """Extract year as int from last_updated (e.g. '2025-01') or version."""
    last_updated = rec.get("last_updated", "")
    if last_updated:
        match = re.search(r"(\d{4})", last_updated)
        if match:
            return int(match.group(1))
    version = rec.get("version", "")
    if version:
        match = re.search(r"(\d{4})", str(version))
        if match:
            return int(match.group(1))
    return 0


def transform_guideline_records(seed_data, embedder):
    """Transform seed JSON records into onco_guidelines schema with embeddings."""
    texts = [rec["text_summary"] for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        records.append({
            "id": rec.get("id", ""),
            "embedding": emb,
            "org": rec.get("organization", ""),
            "cancer_type": rec.get("cancer_type", ""),
            "version": rec.get("version", ""),
            "year": _extract_year(rec),
            "key_recommendations": rec.get("key_recommendations", ""),
            "text_summary": rec.get("text_summary", ""),
            "evidence_level": rec.get("evidence_level", ""),
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "guideline_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Guideline Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} guideline records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_guidelines")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_guidelines")
    print(f"  onco_guidelines currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting guideline seed data...")
    records = transform_guideline_records(seed_data, embedder)
    count = manager.insert_batch("onco_guidelines", records)

    try:
        stats = manager.get_collection_stats("onco_guidelines")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} guideline records")
    print(f"  onco_guidelines now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
