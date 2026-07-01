#!/usr/bin/env python3
"""Ingest imaging AI literature from PubMed into imaging_literature.

Fetches abstracts via NCBI E-utilities, classifies by imaging modality,
body region, and AI task, generates BGE-small embeddings, and stores
in the imaging_literature Milvus collection.

Usage:
    python scripts/ingest_pubmed.py
    python scripts/ingest_pubmed.py --max-results 500
    python scripts/ingest_pubmed.py --query "chest CT segmentation" --max-results 200
    python scripts/ingest_pubmed.py --dry-run

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
        description="Ingest imaging AI literature from PubMed"
    )
    parser.add_argument(
        "--query",
        default=(
            '"medical imaging artificial intelligence"[Title/Abstract] OR '
            '"radiology AI"[Title/Abstract] OR '
            '"deep learning imaging"[Title/Abstract] OR '
            '"CT segmentation"[Title/Abstract] OR '
            '"MRI AI detection"[Title/Abstract]'
        ),
        help="PubMed search query",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=5000,
        help="Maximum number of abstracts to fetch",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=64,
        help="Embedding and insert batch size",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Fetch and parse but don't store in Milvus",
    )
    parser.add_argument(
        "--host", default=None, help="Milvus host (default: localhost)"
    )
    parser.add_argument(
        "--port", type=int, default=None, help="Milvus port (default: 19530)"
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  Imaging Intelligence Agent — PubMed Ingest")
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

    # --- Step 4: Parse into ImagingLiterature models ---
    logger.info("Parsing articles into ImagingLiterature models...")
    from src.ingest.literature_parser import PubMedImagingIngestPipeline
    from src.collections import ImagingCollectionManager

    class DummyEmbedder:
        def encode(self, texts, **kwargs):
            return [[0.0] * 384 for _ in texts]

    dummy_manager = None
    if not args.dry_run:
        dummy_manager = ImagingCollectionManager(host=args.host, port=args.port)

    pipeline = PubMedImagingIngestPipeline(
        collection_manager=dummy_manager,
        embedder=DummyEmbedder(),
        pubmed_client=client,
    )
    records = pipeline.parse(articles)
    logger.info(f"Parsed {len(records)} ImagingLiterature records")

    # Show modality distribution
    from collections import Counter

    modality_counts = Counter(r.modality.value for r in records)
    region_counts = Counter(r.body_region.value for r in records)
    logger.info(f"Modality distribution: {dict(modality_counts)}")
    logger.info(f"Body region distribution: {dict(region_counts.most_common(10))}")

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

    # --- Step 6: Connect to Milvus and store ---
    logger.info("Connecting to Milvus...")
    manager = ImagingCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Recreate pipeline with real embedder and manager
    pipeline = PubMedImagingIngestPipeline(
        collection_manager=manager,
        embedder=model,
        pubmed_client=client,
    )

    logger.info(
        f"Embedding and storing {len(records)} records "
        f"(batch_size={args.batch_size})..."
    )
    count = pipeline.embed_and_store(
        records, "imaging_literature", batch_size=args.batch_size
    )

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
    print(f"  PubMed ingest complete!")
    print(f"  Records stored: {count:,}")
    print(f"  Time elapsed: {elapsed:.1f}s")
    print("=" * 65)


if __name__ == "__main__":
    main()
