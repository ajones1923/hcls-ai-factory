#!/usr/bin/env python3
"""Bulk ingest ClinicalTrials.gov oncology trial data.

Fetches precision oncology clinical trials from the ClinicalTrials.gov
v2 API, extracts biomarker eligibility criteria, generates embeddings,
and stores in the onco_trials Milvus collection.

Usage:
    python scripts/ingest_clinical_trials.py [--max-results N] [--query QUERY]

Options:
    --max-results N    Maximum trials to fetch (default: 2000)
    --query QUERY      ClinicalTrials.gov search query
"""

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from sentence_transformers import SentenceTransformer

from src.collections import OncoCollectionManager
from src.ingest.clinical_trials_parser import ClinicalTrialsIngestPipeline


class SimpleEmbedder:
    """Thin wrapper around SentenceTransformer for pipeline compatibility."""

    def __init__(self):
        self.model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    def encode(self, texts):
        return self.model.encode(texts).tolist()


def main():
    parser = argparse.ArgumentParser(
        description="Bulk ingest ClinicalTrials.gov oncology trials"
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=2000,
        help="Maximum trials to fetch (default: 2000)",
    )
    parser.add_argument(
        "--query",
        type=str,
        default="precision oncology targeted therapy biomarker",
        help="ClinicalTrials.gov search query",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Oncology Intelligence -- ClinicalTrials.gov Bulk Ingest")
    print("=" * 60)

    print("\n[1/4] Connecting to Milvus...")
    manager = OncoCollectionManager()
    manager.connect()
    try:
        stats = manager.get_collection_stats("onco_trials")
        existing = stats.get("num_entities", 0)
    except Exception:
        existing = 0
        manager.create_collection("onco_trials")
    print(f"  onco_trials currently has {existing} records")

    print("\n[2/4] Loading BGE-small-en-v1.5 embedder...")
    embedder = SimpleEmbedder()

    print("\n[3/4] Preparing ClinicalTrials.gov query...")
    print(f"  Query: {args.query}")
    print(f"  Max results: {args.max_results}")

    pipeline = ClinicalTrialsIngestPipeline(
        collection_manager=manager, embedder=embedder
    )

    print("\n[4/4] Running ClinicalTrials.gov ingest pipeline...")
    count = pipeline.run(query=args.query, max_results=args.max_results)

    try:
        stats = manager.get_collection_stats("onco_trials")
        final = stats.get("num_entities", 0)
    except Exception:
        final = existing + count

    print(f"\n{'=' * 60}")
    print(f"DONE: Ingested {count} clinical trials")
    print(f"  onco_trials now has {final} records (was {existing})")
    print(f"{'=' * 60}")

    manager.disconnect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
