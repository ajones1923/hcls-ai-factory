#!/usr/bin/env python3
"""Seed the cart_assays collection with curated assay data from published CAR-T papers.

Loads assay results from data/reference/assay_seed_data.json and ingests them
into the cart_assays Milvus collection via the AssayIngestPipeline.

Data sources include:
  - ELIANA (tisagenlecleucel, B-ALL) — Maude et al., NEJM 2018
  - ZUMA-1 (axicabtagene, DLBCL) — Neelapu et al., NEJM 2017
  - ZUMA-2 (brexucabtagene, MCL) — Wang et al., NEJM 2020
  - TRANSFORM (lisocabtagene, LBCL) — Kamdar et al., Lancet 2022
  - KarMMa (idecabtagene, MM) — Munshi et al., NEJM 2021
  - CARTITUDE-1 (ciltacabtagene, MM) — Berdeja et al., Lancet 2021
  - Preclinical characterization studies for all 6 FDA-approved products
  - Resistance mechanism studies (CD19 loss, BCMA loss, trogocytosis)
  - Comparative studies (CD28 vs 4-1BB, dual-targeting, armored CARs)

Usage:
    python3 scripts/seed_assays.py

Author: Adam Jones
Date: February 2026
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer
from src.collections import CARTCollectionManager
from src.ingest.assay_parser import AssayIngestPipeline


class SimpleEmbedder:
    """Wrapper around SentenceTransformer with encode() API."""
    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    seed_file = PROJECT_ROOT / "data" / "reference" / "assay_seed_data.json"
    if not seed_file.exists():
        print(f"ERROR: Seed file not found: {seed_file}")
        return 1

    print("=" * 60)
    print("CAR-T Assay Data Seeder")
    print("=" * 60)

    # Connect to Milvus
    print("\n[1/3] Connecting to Milvus...")
    manager = CARTCollectionManager()
    manager.connect()

    stats = manager.get_collection_stats()
    existing = stats.get("cart_assays", 0)
    print(f"  cart_assays currently has {existing} records")

    if existing > 0:
        print(f"  WARNING: Collection already has {existing} records.")
        print("  Proceeding — duplicates will be skipped by Milvus primary key.")

    # Load embedder
    print("\n[2/3] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    # Run ingest pipeline
    print("\n[3/3] Ingesting assay seed data...")
    pipeline = AssayIngestPipeline(
        collection_manager=manager,
        embedder=embedder,
    )
    count = pipeline.run(data_file=str(seed_file))

    # Verify
    stats = manager.get_collection_stats()
    final = stats.get("cart_assays", 0)

    print(f"\n{'=' * 60}")
    print(f"DONE: Inserted {count} assay records")
    print(f"  cart_assays now has {final} records")
    print(f"  Seed file: {seed_file}")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
