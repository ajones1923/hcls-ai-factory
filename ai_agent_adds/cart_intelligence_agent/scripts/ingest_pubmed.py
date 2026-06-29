#!/usr/bin/env python3
"""Ingest CAR-T literature from PubMed into the cart_literature collection.

Fetches abstracts via NCBI E-utilities, classifies by CAR-T development
stage, generates BGE-small embeddings, and stores in Milvus.

Usage:
    python scripts/ingest_pubmed.py --max-results 5000
    python scripts/ingest_pubmed.py --query "CD19 CAR-T" --max-results 500
    python scripts/ingest_pubmed.py --dry-run  # Preview without storing

Author: Adam Jones
Date: February 2026
"""

import argparse
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger


def main():
    parser = argparse.ArgumentParser(
        description="Ingest CAR-T literature from PubMed"
    )
    parser.add_argument(
        "--query",
        default=(
            '"chimeric antigen receptor"[Title/Abstract] OR '
            '"CAR-T"[Title/Abstract] OR '
            '"CAR T cell"[Title/Abstract]'
        ),
        help="PubMed search query",
    )
    parser.add_argument(
        "--max-results", type=int, default=5000,
        help="Maximum number of abstracts to fetch",
    )
    parser.add_argument(
        "--batch-size", type=int, default=64,
        help="Embedding and insert batch size",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Fetch and parse but don't store in Milvus",
    )
    parser.add_argument(
        "--host", default=None,
        help="Milvus host (default: localhost)",
    )
    parser.add_argument(
        "--port", type=int, default=None,
        help="Milvus port (default: 19530)",
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  CAR-T Intelligence Agent — PubMed Ingest")
    print("=" * 65)
    print(f"  Query: {args.query[:80]}...")
    print(f"  Max results: {args.max_results}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  Dry run: {args.dry_run}")
    print()

    start_time = time.time()

    # --- Step 1: Initialize PubMed client ---
    logger.info("Initializing PubMed client...")
    from src.utils.pubmed_client import PubMedClient

    client = PubMedClient()

    # --- Step 2: Search for PMIDs ---
    logger.info(f"Searching PubMed for up to {args.max_results} articles...")
    pmids = client.search(args.query, max_results=args.max_results)
    logger.info(f"Found {len(pmids)} PMIDs")

    if not pmids:
        logger.warning("No PMIDs found. Exiting.")
        return

    # --- Step 3: Fetch abstracts ---
    logger.info(f"Fetching abstracts for {len(pmids)} PMIDs...")
    articles = client.fetch_abstracts(pmids)
    logger.info(f"Fetched {len(articles)} article records")

    if not articles:
        logger.warning("No articles fetched. Exiting.")
        return

    # --- Step 4: Parse into CARTLiterature models ---
    logger.info("Parsing articles into CARTLiterature models...")
    from src.ingest.literature_parser import PubMedIngestPipeline
    from src.collections import CARTCollectionManager

    # We need a temporary pipeline just for parsing (no embedder needed yet)
    # Create a dummy embedder for the pipeline init
    class DummyEmbedder:
        def encode(self, texts):
            return [[0.0] * 384 for _ in texts]

    dummy_manager = None
    if not args.dry_run:
        dummy_manager = CARTCollectionManager(host=args.host, port=args.port)

    pipeline = PubMedIngestPipeline(
        collection_manager=dummy_manager,
        embedder=DummyEmbedder(),
        pubmed_client=client,
    )
    records = pipeline.parse(articles)
    logger.info(f"Parsed {len(records)} CARTLiterature records")

    # Show stage distribution
    from collections import Counter
    stage_counts = Counter(r.cart_stage.value for r in records)
    antigen_counts = Counter(r.target_antigen for r in records if r.target_antigen)
    logger.info(f"Stage distribution: {dict(stage_counts)}")
    logger.info(f"Top antigens: {dict(antigen_counts.most_common(10))}")

    if args.dry_run:
        logger.info("Dry run complete. No data stored.")
        elapsed = time.time() - start_time
        print(f"\n  Dry run completed in {elapsed:.1f}s")
        print(f"  {len(records)} records parsed (would be stored)")
        return

    # --- Step 5: Initialize embedder and Milvus ---
    logger.info("Loading BGE-small-en-v1.5 embedding model...")
    from sentence_transformers import SentenceTransformer

    model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    class SimpleEmbedder:
        def __init__(self, st_model):
            self._model = st_model

        def encode(self, texts):
            return self._model.encode(texts).tolist()

    embedder = SimpleEmbedder(model)

    # --- Step 6: Connect to Milvus and store ---
    logger.info("Connecting to Milvus...")
    manager = CARTCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Recreate pipeline with real embedder and manager
    pipeline = PubMedIngestPipeline(
        collection_manager=manager,
        embedder=embedder,
        pubmed_client=client,
    )

    logger.info(f"Embedding and storing {len(records)} records (batch_size={args.batch_size})...")
    count = pipeline.embed_and_store(records, "cart_literature", batch_size=args.batch_size)

    elapsed = time.time() - start_time
    logger.info(f"Ingest complete: {count} records stored in {elapsed:.1f}s")

    # Show final stats
    stats = manager.get_collection_stats()
    logger.info("Final collection stats:")
    for name, cnt in stats.items():
        logger.info(f"  {name}: {cnt:,} records")

    manager.disconnect()

    print()
    print("=" * 65)
    print("  PubMed ingest complete!")
    print(f"  Records stored: {count:,}")
    print(f"  Time elapsed: {elapsed:.1f}s")
    print("=" * 65)


if __name__ == "__main__":
    main()
