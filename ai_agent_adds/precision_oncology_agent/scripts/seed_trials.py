#!/usr/bin/env python3
"""Seed the onco_trials collection with curated clinical trial data.

Loads trial_seed_data.json, generates BGE-small-en-v1.5 embeddings
from each text_summary, and inserts records into the onco_trials
Milvus collection.

Usage: python3 scripts/seed_trials.py
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


def _parse_enrollment(value):
    """Convert enrollment to int, handling strings like '1,200'."""
    if value is None:
        return 0
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    # String: strip commas and whitespace, then convert
    cleaned = re.sub(r"[,\s]", "", str(value))
    try:
        return int(cleaned)
    except (ValueError, TypeError):
        return 0


def _to_str(value):
    """Convert a value to string, joining lists with commas."""
    if value is None:
        return ""
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return str(value)


def transform_trial_records(seed_data, embedder):
    """Transform seed JSON records into onco_trials schema with embeddings."""
    texts = [rec["text_summary"] for rec in seed_data]
    embeddings = embedder.encode(texts)

    records = []
    for rec, emb in zip(seed_data, embeddings):
        # Use nct_id as id if id is missing, or vice versa
        rec_id = rec.get("id", "") or rec.get("nct_id", "")
        # Truncate to 20 chars for VARCHAR max 20
        rec_id = str(rec_id)[:20]

        records.append({
            "id": rec_id,
            "embedding": emb,
            "title": rec.get("title", ""),
            "text_summary": rec.get("text_summary", ""),
            "phase": _to_str(rec.get("phase", "")),
            "status": rec.get("status", ""),
            "sponsor": rec.get("sponsor", ""),
            "cancer_types": _to_str(rec.get("cancer_types", "")),
            "biomarker_criteria": rec.get("biomarker_criteria", ""),
            "enrollment": _parse_enrollment(rec.get("enrollment")),
            "start_year": 0,
            "outcome_summary": rec.get("outcome_summary", ""),
        })
    return records


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "trial_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Clinical Trial Data Seeder")
    print("=" * 60)

    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} trial records from seed file")

    print("\n[1/3] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_trials")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_trials")
    print(f"  onco_trials currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting trial seed data...")
    records = transform_trial_records(seed_data, embedder)
    count = manager.insert_batch("onco_trials", records)

    try:
        stats = manager.get_collection_stats("onco_trials")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} trial records")
    print(f"  onco_trials now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
