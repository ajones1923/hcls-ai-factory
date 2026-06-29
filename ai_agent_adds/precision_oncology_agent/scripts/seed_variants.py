#!/usr/bin/env python3
"""Seed the onco_variants collection with curated actionable variant data.

Loads variant_seed_data.json and uses CIViCIngestPipeline to parse and
embed each record into the onco_variants Milvus collection.

Usage: python3 scripts/seed_variants.py
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer

from src.collections import OncoCollectionManager
from src.ingest.civic_parser import CIViCIngestPipeline


class SimpleEmbedder:
    """Thin wrapper around SentenceTransformer for pipeline compatibility."""

    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "variant_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("Oncology Intelligence -- Variant Data Seeder")
    print("=" * 60)

    # Load seed data
    with open(seed_file, "r", encoding="utf-8") as f:
        seed_data = json.load(f)
    print(f"\n  Loaded {len(seed_data)} variant records from seed file")

    print("\n[1/4] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_variants")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_variants")
    print(f"  onco_variants currently has {existing} records")

    print("\n[2/4] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/4] Preparing seed data for CIViC pipeline format...")
    # Transform seed records into a format the CIViC parser can handle
    # We'll use the parse method's raw_data format, simulating CIViC-like records
    civic_records = []
    for record in seed_data:
        civic_records.append({
            "id": record["id"],
            "gene": record["gene"],
            "variant_name": record["variant"],
            "variant_type": record["variant_type"],
            "cancer_type": record["cancer_type"],
            "evidence_level": record["evidence_level"],
            "drugs": record["drugs"],
            "civic_id": record.get("civic_id", ""),
            "vrs_id": "",
            "text_summary": record["text_summary"],
            "text": record["text_summary"],
            "clinical_significance": record.get("clinical_significance", ""),
            "allele_frequency": 0.0,
            "source_type": "seed",
        })

    print(f"  Prepared {len(civic_records)} records for embedding")

    print("\n[4/4] Embedding and inserting variant seed data...")
    pipeline = CIViCIngestPipeline(collection_manager=manager, embedder=embedder)
    # Use embed_and_store directly with our pre-formatted records
    count = pipeline.embed_and_store(civic_records)

    try:
        stats = manager.get_collection_stats("onco_variants")
        final = stats.get("num_entities", 0)
    except Exception:
        final = count

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} variant records")
    print(f"  onco_variants now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
