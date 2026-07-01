#!/usr/bin/env python3
"""Bulk ingest PubMed oncology literature into onco_literature.

Fetches precision oncology articles from PubMed via E-utilities API,
extracts cancer type and gene mentions, generates embeddings, and
stores in the onco_literature Milvus collection.

Usage:
    python scripts/ingest_pubmed.py [--max-results N] [--query QUERY]

Options:
    --max-results N    Maximum number of PubMed articles to fetch (default: 5000)
    --query QUERY      PubMed search query (default: precision oncology)
"""

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer

from src.collections import OncoCollectionManager
from src.ingest.literature_parser import PubMedIngestPipeline


class SimpleEmbedder:
    """Thin wrapper around SentenceTransformer for pipeline compatibility."""

    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    parser = argparse.ArgumentParser(
        description="Bulk ingest PubMed oncology literature"
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=5000,
        help="Maximum number of PubMed articles to fetch (default: 5000)",
    )
    parser.add_argument(
        "--query",
        type=str,
        default="precision oncology targeted therapy biomarker",
        help="PubMed search query",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Oncology Intelligence -- PubMed Bulk Ingest")
    print("=" * 60)

    print("\n[1/4] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_literature")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_literature")
    print(f"  onco_literature currently has {existing} records")

    print("\n[2/4] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/4] Preparing PubMed query...")
    print(f"  Query: {args.query}")
    print(f"  Max results: {args.max_results}")

    pipeline = PubMedIngestPipeline(collection_manager=manager, embedder=embedder)

    print("\n[4/4] Running PubMed ingest pipeline...")
    count = pipeline.run(query=args.query, max_results=args.max_results)

    try:
        stats = manager.get_collection_stats("onco_literature")
        final = stats.get("num_entities", 0)
    except Exception:
        final = existing + count

    print(f"\n{'=' * 60}")
    print(f"DONE: Ingested {count} PubMed articles")
    print(f"  onco_literature now has {final} records (was {existing})")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
