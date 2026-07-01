#!/usr/bin/env python3
"""Seed the onco_cases collection with curated patient case data.

Loads cases_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_summary, and inserts records into the onco_cases
Milvus collection.

Usage: python3 scripts/seed_cases.py
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


def transform_case_records(seed_data, embedder):
    """Transform seed JSON records into onco_cases schema with embeddings."""
    texts = [rec.get("text_summary", "") for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        variants = rec.get("variants", [])
        if isinstance(variants, list):
            variants = ", ".join(variants)

        biomarkers = rec.get("biomarkers", {})
        if isinstance(biomarkers, dict):
            biomarkers = ", ".join(f"{k}: {v}" for k, v in biomarkers.items())
        elif isinstance(biomarkers, list):
            biomarkers = ", ".join(biomarkers)

        prior_therapies = rec.get("prior_therapies", [])
        if isinstance(prior_therapies, list):
            prior_therapies = ", ".join(prior_therapies)

        records.append({
            "id": str(rec["id"])[:100],
            "embedding": emb,
            "patient_id": str(rec.get("patient_id", ""))[:100],
            "cancer_type": str(rec.get("cancer_type", ""))[:50],
            "stage": str(rec.get("stage", ""))[:20],
            "variants": str(variants)[:1000],
            "biomarkers": str(biomarkers)[:1000],
            "prior_therapies": str(prior_therapies)[:500],
            "text_summary": str(rec.get("text_summary", ""))[:3000],
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "cases_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Case Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} case records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_cases")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_cases")
    print(f"  onco_cases currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting case seed data...")
    records = transform_case_records(seed_data, embedder)
    count = manager.insert_batch("onco_cases", records)

    try:
        stats = manager.get_collection_stats("onco_cases")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} case records")
    print(f"  onco_cases now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
