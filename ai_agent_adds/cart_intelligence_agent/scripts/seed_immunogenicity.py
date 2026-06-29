#!/usr/bin/env python3
"""Seed immunogenicity data into cart_biomarkers and cart_sequences.

Loads HLA/ADA-focused biomarker records and deimmunized/humanized scFv
sequence records from their respective JSON seed files.

Usage:
    python3 scripts/seed_immunogenicity.py
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer
from src.collections import CARTCollectionManager
from src.ingest.biomarker_parser import BiomarkerIngestPipeline
from src.ingest.sequence_parser import SequenceIngestPipeline


class SimpleEmbedder:
    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")
    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    biomarker_file = PROJECT_ROOT / "data" / "reference" / "immunogenicity_biomarker_seed.json"
    sequence_file = PROJECT_ROOT / "data" / "reference" / "immunogenicity_sequence_seed.json"

    missing = []
    if not biomarker_file.exists():
        missing.append(str(biomarker_file))
    if not sequence_file.exists():
        missing.append(str(sequence_file))
    if missing:
        print(f"ERROR: Seed file(s) not found: {', '.join(missing)}")
        return 1

    print("=" * 60)
    print("CAR-T Immunogenicity Data Seeder")
    print("=" * 60)

    print("\n[1/4] Connecting to Milvus...")
    manager = CARTCollectionManager()
    manager.connect()

    stats = manager.get_collection_stats()
    bio_existing = stats.get("cart_biomarkers", 0)
    seq_existing = stats.get("cart_sequences", 0)
    print(f"  cart_biomarkers currently has {bio_existing} records")
    print(f"  cart_sequences currently has {seq_existing} records")

    print("\n[2/4] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/4] Ingesting immunogenicity biomarker data...")
    bio_pipeline = BiomarkerIngestPipeline(collection_manager=manager, embedder=embedder)
    bio_count = bio_pipeline.run(data_file=str(biomarker_file))

    print("\n[4/4] Ingesting immunogenicity sequence data...")
    seq_pipeline = SequenceIngestPipeline(collection_manager=manager, embedder=embedder)
    seq_count = seq_pipeline.run(data_file=str(sequence_file))

    stats = manager.get_collection_stats()
    bio_final = stats.get("cart_biomarkers", 0)
    seq_final = stats.get("cart_sequences", 0)

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {bio_count} immunogenicity biomarker records")
    print(f"  cart_biomarkers now has {bio_final} records")
    print(f"DONE: Inserted {seq_count} immunogenicity sequence records")
    print(f"  cart_sequences now has {seq_final} records")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
