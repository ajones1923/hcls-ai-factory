#!/usr/bin/env python3
"""Bulk ingest CIViC (Clinical Interpretation of Variants in Cancer) data.

Fetches variant and evidence records from the CIViC public API, parses
them into onco_variants schema, generates embeddings, and stores in Milvus.

Usage:
    python scripts/ingest_civic.py [--max-results N] [--gene GENE]

Options:
    --max-results N    Maximum number of CIViC variants to fetch (default: 5000)
    --gene GENE        Filter by specific gene symbol (e.g., BRAF, EGFR)
"""

import argparse
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
    parser = argparse.ArgumentParser(
        description="Bulk ingest CIViC variant and evidence data"
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=5000,
        help="Maximum number of CIViC variants to fetch (default: 5000)",
    )
    parser.add_argument(
        "--gene",
        type=str,
        default=None,
        help="Filter by gene symbol (e.g., BRAF, EGFR)",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Oncology Intelligence -- CIViC Bulk Ingest")
    print("=" * 60)

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

    print(f"\n[3/4] Fetching CIViC data (max_results={args.max_results})...")
    if args.gene:
        print(f"  Filtering by gene: {args.gene}")

    pipeline = CIViCIngestPipeline(collection_manager=manager, embedder=embedder)

    print("\n[4/4] Running CIViC ingest pipeline...")
    kwargs = {"max_results": args.max_results}
    if args.gene:
        kwargs["query"] = args.gene

    count = pipeline.run(**kwargs)

    try:
        stats = manager.get_collection_stats("onco_variants")
        final = stats.get("num_entities", 0)
    except Exception:
        final = existing + count

    print(f"\n{'=' * 60}")
    print(f"DONE: Ingested {count} CIViC records")
    print(f"  onco_variants now has {final} records (was {existing})")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
