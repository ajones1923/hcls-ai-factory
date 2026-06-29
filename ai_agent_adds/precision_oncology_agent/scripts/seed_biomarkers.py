#!/usr/bin/env python3
"""Seed the onco_biomarkers collection with curated biomarker data.

Loads biomarker_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_summary, and inserts records into the onco_biomarkers
Milvus collection.

Usage: python3 scripts/seed_biomarkers.py
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


def transform_biomarker_records(seed_data, embedder):
    """Transform seed JSON records into onco_biomarkers schema with embeddings."""
    texts = [rec["text_summary"] for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        records.append({
            "id": rec["id"],
            "embedding": emb,
            "name": rec["biomarker_name"],
            "biomarker_type": rec["biomarker_type"],
            "cancer_types": rec["cancer_types"],
            "predictive_value": rec["predictive_value"],
            "testing_method": rec["testing_method"],
            "clinical_cutoff": rec["clinical_cutoff"],
            "text_summary": rec["text_summary"],
            "evidence_level": "A",
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "biomarker_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Biomarker Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} biomarker records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_biomarkers")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_biomarkers")
    print(f"  onco_biomarkers currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting biomarker seed data...")
    records = transform_biomarker_records(seed_data, embedder)
    count = manager.insert_batch("onco_biomarkers", records)

    try:
        stats = manager.get_collection_stats("onco_biomarkers")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} biomarker records")
    print(f"  onco_biomarkers now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
