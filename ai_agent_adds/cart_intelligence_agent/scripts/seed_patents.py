#!/usr/bin/env python3
"""Seed the cart_literature collection with curated patent data.

Loads patent records from data/reference/patent_seed_data.json, creates
CARTLiterature models with source_type=patent, embeds the text chunks,
and inserts into the cart_literature Milvus collection.

Usage:
    python3 scripts/seed_patents.py
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer
from src.collections import CARTCollectionManager
from src.models import CARTLiterature, CARTStage, SourceType


class SimpleEmbedder:
    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")
    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "patent_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("CAR-T Patent Data Seeder")
    print("=" * 60)

    # Load patent records
    with open(seed_file, "r") as f:
        raw_records = json.load(f)
    print(f"\nLoaded {len(raw_records)} patent records from {seed_file.name}")

    # Parse into CARTLiterature models
    records = []
    for data in raw_records:
        try:
            record = CARTLiterature(
                id=data["id"],
                title=data["title"],
                text_chunk=data["text_chunk"],
                source_type=SourceType.PATENT,
                year=data["year"],
                cart_stage=CARTStage(data["cart_stage"]),
                target_antigen=data.get("target_antigen", ""),
                disease=data.get("disease", ""),
                keywords=data.get("keywords", ""),
                journal=data.get("journal", ""),
            )
            records.append(record)
        except Exception as e:
            print(f"  WARNING: Failed to parse {data.get('id', '?')}: {e}")

    print(f"Parsed {len(records)} valid patent records")

    print("\n[1/3] Connecting to Milvus...")
    manager = CARTCollectionManager()
    manager.connect()

    stats = manager.get_collection_stats()
    existing = stats.get("cart_literature", 0)
    print(f"  cart_literature currently has {existing} records")

    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/3] Embedding and inserting patent records...")
    total_inserted = 0
    batch_size = 32

    for i in range(0, len(records), batch_size):
        batch = records[i : i + batch_size]
        texts = [r.to_embedding_text() for r in batch]
        embeddings = embedder.encode(texts)

        batch_dicts = []
        for record, embedding in zip(batch, embeddings):
            d = record.model_dump()
            d["embedding"] = embedding
            # Convert enums to string values
            for key, value in d.items():
                if hasattr(value, "value"):
                    d[key] = value.value
            batch_dicts.append(d)

        manager.insert_batch("cart_literature", batch_dicts)
        total_inserted += len(batch_dicts)
        print(f"  Inserted batch {i // batch_size + 1} ({total_inserted}/{len(records)})")

    stats = manager.get_collection_stats()
    final = stats.get("cart_literature", 0)

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {total_inserted} patent records")
    print(f"  cart_literature now has {final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
