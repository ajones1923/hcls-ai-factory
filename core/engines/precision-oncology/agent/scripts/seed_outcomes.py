#!/usr/bin/env python3
"""Seed the onco_outcomes collection with curated treatment outcome data.

Loads outcome_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_summary, and inserts records into the onco_outcomes
Milvus collection.

Usage: python3 scripts/seed_outcomes.py
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


def _to_float(value):
    """Convert a value to float, handling None and strings."""
    if value is None:
        return 0.0
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0


def transform_outcome_records(seed_data, embedder):
    """Transform seed JSON records into onco_outcomes schema with embeddings."""
    texts = [rec["text_summary"] for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        records.append({
            "id": rec.get("id", ""),
            "embedding": emb,
            "case_id": rec.get("case_id", ""),
            "therapy": rec.get("therapy", ""),
            "cancer_type": rec.get("cancer_type", ""),
            "response": rec.get("response", ""),
            "duration_months": _to_float(rec.get("duration_months")),
            "toxicities": rec.get("toxicities", ""),
            "biomarkers_at_baseline": rec.get("biomarkers_at_baseline", ""),
            "text_summary": rec.get("text_summary", ""),
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "outcome_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Outcome Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} outcome records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_outcomes")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_outcomes")
    print(f"  onco_outcomes currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting outcome seed data...")
    records = transform_outcome_records(seed_data, embedder)
    count = manager.insert_batch("onco_outcomes", records)

    try:
        stats = manager.get_collection_stats("onco_outcomes")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} outcome records")
    print(f"  onco_outcomes now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
