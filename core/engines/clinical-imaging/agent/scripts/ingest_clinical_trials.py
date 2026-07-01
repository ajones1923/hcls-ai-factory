#!/usr/bin/env python3
"""Ingest imaging AI clinical trials from ClinicalTrials.gov.

Fetches trial records via the ClinicalTrials.gov API v2, parses JSON
into ImagingTrial models, generates BGE-small embeddings, and stores
in the imaging_trials Milvus collection.

Usage:
    python scripts/ingest_clinical_trials.py
    python scripts/ingest_clinical_trials.py --max-results 500
    python scripts/ingest_clinical_trials.py --dry-run

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
        description="Ingest imaging AI clinical trials from ClinicalTrials.gov"
    )
    parser.add_argument(
        "--condition",
        default="imaging AI",
        help="Condition search term",
    )
    parser.add_argument(
        "--intervention",
        default="artificial intelligence radiology",
        help="Intervention search term",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=500,
        help="Maximum number of trials to fetch",
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
    print("  Imaging Intelligence Agent — ClinicalTrials.gov Ingest")
    print("=" * 65)
    print(f"  Condition: {args.condition}")
    print(f"  Intervention: {args.intervention}")
    print(f"  Max results: {args.max_results}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  Dry run: {args.dry_run}")
    print()

    start_time = time.time()

    # --- Step 1: Initialize pipeline ---
    logger.info("Initializing clinical trials ingest pipeline...")
    from src.ingest.clinical_trials_parser import ImagingTrialsIngestPipeline
    from src.collections import ImagingCollectionManager

    class DummyEmbedder:
        def encode(self, texts, **kwargs):
            return [[0.0] * 384 for _ in texts]

    if args.dry_run:
        # Dry run: no Milvus connection needed
        pipeline = ImagingTrialsIngestPipeline(
            collection_manager=None,
            embedder=DummyEmbedder(),
        )

        logger.info("Fetching trials from ClinicalTrials.gov...")
        raw_data = pipeline.fetch(
            condition=args.condition,
            intervention=args.intervention,
            max_results=args.max_results,
        )
        records = pipeline.parse(raw_data)
        logger.info(f"Parsed {len(records)} ImagingTrial records")

        # Show distribution
        from collections import Counter

        phase_counts = Counter(r.phase.value for r in records)
        status_counts = Counter(r.status.value for r in records)
        modality_counts = Counter(r.modality.value for r in records)
        logger.info(f"Phase distribution: {dict(phase_counts)}")
        logger.info(f"Status distribution: {dict(status_counts)}")
        logger.info(f"Modality distribution: {dict(modality_counts)}")

        elapsed = time.time() - start_time
        print(f"\n  Dry run completed in {elapsed:.1f}s")
        print(f"  {len(records)} records parsed (would be stored)")
        return

    # --- Step 2: Initialize embedder and Milvus ---
    logger.info("Loading BGE-small-en-v1.5 embedding model...")
    from sentence_transformers import SentenceTransformer

    model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    logger.info("Connecting to Milvus...")
    manager = ImagingCollectionManager(host=args.host, port=args.port)
    manager.connect()

    pipeline = ImagingTrialsIngestPipeline(
        collection_manager=manager,
        embedder=model,
    )

    # --- Step 3: Run the full pipeline ---
    logger.info("Running ingest pipeline...")
    count = pipeline.run(
        condition=args.condition,
        intervention=args.intervention,
        max_results=args.max_results,
        batch_size=args.batch_size,
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
    print(f"  ClinicalTrials.gov ingest complete!")
    print(f"  Records stored: {count:,}")
    print(f"  Time elapsed: {elapsed:.1f}s")
    print("=" * 65)


if __name__ == "__main__":
    main()
