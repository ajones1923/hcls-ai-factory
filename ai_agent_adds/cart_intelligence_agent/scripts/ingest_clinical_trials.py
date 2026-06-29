#!/usr/bin/env python3
"""Ingest CAR-T clinical trials from ClinicalTrials.gov.

Fetches trial data via the ClinicalTrials.gov API v2, parses into
ClinicalTrial models, generates embeddings, and stores in Milvus.

Usage:
    python scripts/ingest_clinical_trials.py --max-results 1500
    python scripts/ingest_clinical_trials.py --condition "CAR-T CD19"
    python scripts/ingest_clinical_trials.py --dry-run

Author: Adam Jones
Date: February 2026
"""

import argparse
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger


def main():
    parser = argparse.ArgumentParser(
        description="Ingest CAR-T clinical trials from ClinicalTrials.gov"
    )
    parser.add_argument(
        "--max-results", type=int, default=1500,
        help="Maximum trials to fetch",
    )
    parser.add_argument(
        "--condition", type=str, default="CAR-T",
        help="Condition search term (default: CAR-T)",
    )
    parser.add_argument(
        "--intervention", type=str, default="chimeric antigen receptor",
        help="Intervention search term",
    )
    parser.add_argument(
        "--page-size", type=int, default=100,
        help="Results per API page (max 1000)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=64,
        help="Embedding and insert batch size",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Fetch and parse but don't store",
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
    print("  CAR-T Intelligence Agent — ClinicalTrials.gov Ingest")
    print("=" * 65)
    print(f"  Condition: {args.condition}")
    print(f"  Intervention: {args.intervention}")
    print(f"  Max results: {args.max_results}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  Dry run: {args.dry_run}")
    print()

    start_time = time.time()

    # --- Step 1: Initialize pipeline components ---
    from src.collections import CARTCollectionManager
    from src.ingest.clinical_trials_parser import ClinicalTrialsIngestPipeline

    # Dummy embedder for parsing phase
    class DummyEmbedder:
        def encode(self, texts):
            return [[0.0] * 384 for _ in texts]

    dummy_manager = None
    if not args.dry_run:
        dummy_manager = CARTCollectionManager(host=args.host, port=args.port)

    pipeline = ClinicalTrialsIngestPipeline(
        collection_manager=dummy_manager,
        embedder=DummyEmbedder(),
    )

    # --- Step 2: Fetch trials from ClinicalTrials.gov ---
    logger.info(f"Fetching up to {args.max_results} CAR-T trials from ClinicalTrials.gov...")
    raw_studies = pipeline.fetch(
        condition=args.condition,
        intervention=args.intervention,
        max_results=args.max_results,
        page_size=args.page_size,
    )
    logger.info(f"Fetched {len(raw_studies)} raw studies")

    if not raw_studies:
        logger.warning("No studies fetched. Exiting.")
        return

    # --- Step 3: Parse into ClinicalTrial models ---
    logger.info("Parsing studies into ClinicalTrial models...")
    records = pipeline.parse(raw_studies)
    logger.info(f"Parsed {len(records)} ClinicalTrial records")

    # Show statistics
    from collections import Counter
    phase_counts = Counter(r.phase.value for r in records)
    status_counts = Counter(r.status.value for r in records)
    antigen_counts = Counter(r.target_antigen for r in records if r.target_antigen)
    gen_counts = Counter(r.car_generation.value for r in records)

    logger.info(f"Phase distribution: {dict(phase_counts)}")
    logger.info(f"Status distribution: {dict(status_counts)}")
    logger.info(f"Top antigens: {dict(antigen_counts.most_common(10))}")
    logger.info(f"CAR generation: {dict(gen_counts)}")

    if args.dry_run:
        elapsed = time.time() - start_time
        logger.info("Dry run complete. No data stored.")
        print(f"\n  Dry run completed in {elapsed:.1f}s")
        print(f"  {len(records)} records parsed (would be stored)")
        return

    # --- Step 4: Initialize real embedder ---
    logger.info("Loading BGE-small-en-v1.5 embedding model...")
    from sentence_transformers import SentenceTransformer

    model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    class SimpleEmbedder:
        def __init__(self, st_model):
            self._model = st_model

        def encode(self, texts):
            return self._model.encode(texts).tolist()

    embedder = SimpleEmbedder(model)

    # --- Step 5: Connect to Milvus and store ---
    logger.info("Connecting to Milvus...")
    manager = CARTCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Recreate pipeline with real embedder and manager
    pipeline = ClinicalTrialsIngestPipeline(
        collection_manager=manager,
        embedder=embedder,
    )

    logger.info(f"Embedding and storing {len(records)} records (batch_size={args.batch_size})...")
    count = pipeline.embed_and_store(records, "cart_trials", batch_size=args.batch_size)

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
    print("  ClinicalTrials.gov ingest complete!")
    print(f"  Records stored: {count:,}")
    print(f"  Time elapsed: {elapsed:.1f}s")
    print("=" * 65)


if __name__ == "__main__":
    main()
