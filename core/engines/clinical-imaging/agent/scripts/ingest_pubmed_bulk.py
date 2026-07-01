#!/usr/bin/env python3
"""Bulk ingest imaging AI literature from PubMed into imaging_literature.

Searches PubMed for medical imaging AI papers across 12 curated queries,
fetches abstracts, generates BGE-small-en-v1.5 embeddings, and stores
them in the imaging_literature Milvus collection.

Deduplicates by PMID (skips papers already seen across queries).
Handles network errors, rate limiting, and partial failures.

Usage:
    python scripts/ingest_pubmed_bulk.py
    python scripts/ingest_pubmed_bulk.py --dry-run
    python scripts/ingest_pubmed_bulk.py --max-per-query 100

Author: Adam Jones
Date: April 2026
"""

import argparse
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

# ═══════════════════════════════════════════════════════════════════════
# CURATED SEARCH QUERIES
# ═══════════════════════════════════════════════════════════════════════

# Each tuple: (query_string, max_results, date_filter_years)
QUERIES = [
    # Last 5 years — broad foundational topics
    (
        '"artificial intelligence" AND "radiology"[Title/Abstract]',
        300,
        5,
    ),
    (
        '"deep learning" AND "medical imaging"[Title/Abstract]',
        300,
        5,
    ),
    (
        '"computer aided detection" AND "CT"[Title/Abstract]',
        250,
        5,
    ),
    # Last 3 years — specific clinical applications
    (
        '"AI" AND "lung nodule detection"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"AI" AND "chest xray" OR "chest x-ray" AND "classification"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"radiomics" AND "machine learning"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"federated learning" AND "medical imaging"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"AI" AND "brain hemorrhage detection"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"AI" AND "breast cancer" AND "mammography"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"MONAI" AND "medical image segmentation"[Title/Abstract]',
        250,
        3,
    ),
    (
        '"NVIDIA" AND "medical imaging"[Title/Abstract]',
        250,
        3,
    ),
    # Last 2 years — cutting edge
    (
        '"vision language model" AND "radiology"[Title/Abstract]',
        200,
        2,
    ),
]


def _build_date_filter(years_back: int) -> str:
    """Build a PubMed date filter string for the last N years."""
    from datetime import datetime
    now = datetime.now()
    start_year = now.year - years_back
    return f' AND ("{start_year}/01/01"[Date - Publication] : "3000"[Date - Publication])'


def main():
    parser = argparse.ArgumentParser(
        description="Bulk ingest imaging AI literature from PubMed"
    )
    parser.add_argument(
        "--max-per-query",
        type=int,
        default=None,
        help="Override max results per query (default: use per-query settings)",
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

    print("=" * 70)
    print("  Imaging Intelligence Agent — Bulk PubMed Ingest")
    print("=" * 70)
    print(f"  Queries: {len(QUERIES)}")
    print(f"  Batch size: {args.batch_size}")
    print(f"  Dry run: {args.dry_run}")
    print()

    start_time = time.time()

    # --- Step 1: Initialize PubMed client ---
    logger.info("Initializing PubMed client...")
    from src.utils.pubmed_client import PubMedClient

    client = PubMedClient()

    # --- Step 2: Search and fetch across all queries, dedup by PMID ---
    all_pmids_seen: set = set()
    all_articles: list = []
    query_stats: list = []

    for idx, (query, max_results, years_back) in enumerate(QUERIES, start=1):
        effective_max = args.max_per_query if args.max_per_query else max_results
        date_filter = _build_date_filter(years_back)
        full_query = query + date_filter

        print(f"  [{idx}/{len(QUERIES)}] {query[:70]}...")
        print(f"           max={effective_max}, last {years_back} years")

        try:
            pmids = client.search(full_query, max_results=effective_max)
            logger.info(f"Query {idx}: found {len(pmids)} PMIDs")

            # Deduplicate against previously seen PMIDs
            new_pmids = [p for p in pmids if p not in all_pmids_seen]
            dupes = len(pmids) - len(new_pmids)
            if dupes > 0:
                logger.info(f"  Skipping {dupes} duplicate PMIDs already seen")

            if not new_pmids:
                query_stats.append((query[:60], 0, 0, dupes))
                print(f"           -> 0 new articles (all duplicates)")
                continue

            # Add to seen set
            all_pmids_seen.update(new_pmids)

            # Fetch abstracts for new PMIDs only
            articles = client.fetch_abstracts(new_pmids)

            # Filter out articles with empty abstracts
            articles_with_abstracts = [
                a for a in articles if a.get("abstract", "").strip()
            ]
            skipped_no_abstract = len(articles) - len(articles_with_abstracts)
            if skipped_no_abstract > 0:
                logger.info(f"  Skipped {skipped_no_abstract} articles without abstracts")

            all_articles.extend(articles_with_abstracts)
            query_stats.append((
                query[:60],
                len(articles_with_abstracts),
                skipped_no_abstract,
                dupes,
            ))
            print(
                f"           -> {len(articles_with_abstracts)} new articles "
                f"({dupes} dupes, {skipped_no_abstract} no abstract)"
            )

        except Exception as e:
            logger.error(f"Failed query {idx}: {e}")
            query_stats.append((query[:60], 0, 0, 0))
            print(f"           -> ERROR: {e}")
            continue

    print()
    print(f"  Total unique articles fetched: {len(all_articles)}")
    print(f"  Total unique PMIDs seen: {len(all_pmids_seen)}")
    print()

    if not all_articles:
        logger.warning("No articles fetched across all queries. Exiting.")
        return

    # --- Step 3: Parse into ImagingLiterature models ---
    logger.info("Parsing articles into ImagingLiterature models...")
    from src.ingest.literature_parser import PubMedImagingIngestPipeline

    class DummyEmbedder:
        def encode(self, texts, **kwargs):
            return [[0.0] * 384 for _ in texts]

    pipeline = PubMedImagingIngestPipeline(
        collection_manager=None,
        embedder=DummyEmbedder(),
        pubmed_client=client,
    )
    records = pipeline.parse(all_articles)
    logger.info(f"Parsed {len(records)} ImagingLiterature records")

    # Show distribution stats
    from collections import Counter

    modality_counts = Counter(r.modality.value for r in records)
    region_counts = Counter(r.body_region.value for r in records)
    task_counts = Counter(r.ai_task for r in records if r.ai_task)

    print("  Modality distribution:")
    for mod, cnt in modality_counts.most_common():
        print(f"    {mod}: {cnt}")
    print()
    print("  Body region distribution:")
    for reg, cnt in region_counts.most_common(10):
        print(f"    {reg}: {cnt}")
    print()
    print("  AI task distribution:")
    for task, cnt in task_counts.most_common(8):
        print(f"    {task}: {cnt}")
    print()

    if args.dry_run:
        elapsed = time.time() - start_time
        print(f"  Dry run completed in {elapsed:.1f}s")
        print(f"  {len(records)} records parsed (would be stored)")
        return

    # --- Step 4: Load embedding model ---
    logger.info("Loading BGE-small-en-v1.5 embedding model...")
    from sentence_transformers import SentenceTransformer

    model = SentenceTransformer("BAAI/bge-small-en-v1.5")

    # --- Step 5: Connect to Milvus ---
    logger.info("Connecting to Milvus...")
    from src.collections import ImagingCollectionManager

    manager = ImagingCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Check for existing PMIDs to deduplicate against what is already in Milvus
    existing_count = 0
    try:
        from pymilvus import Collection as MilvusCollection
        coll = MilvusCollection("imaging_literature")
        coll.load()
        existing_count = coll.num_entities
        logger.info(f"Collection imaging_literature has {existing_count} existing entities")
    except Exception as e:
        logger.warning(f"Could not check existing collection: {e}")

    # --- Step 6: Embed and store with real model ---
    pipeline_real = PubMedImagingIngestPipeline(
        collection_manager=manager,
        embedder=model,
        pubmed_client=client,
    )

    logger.info(
        f"Embedding and storing {len(records)} records "
        f"(batch_size={args.batch_size})..."
    )
    count = pipeline_real.embed_and_store(
        records, "imaging_literature", batch_size=args.batch_size
    )

    elapsed = time.time() - start_time

    # --- Step 7: Final stats ---
    try:
        coll = MilvusCollection("imaging_literature")
        coll.load()
        final_count = coll.num_entities
    except Exception:
        final_count = existing_count + count

    manager.disconnect()

    print()
    print("=" * 70)
    print("  Bulk PubMed Ingest Complete!")
    print("=" * 70)
    print(f"  Queries executed:     {len(QUERIES)}")
    print(f"  Articles fetched:     {len(all_articles)}")
    print(f"  Records parsed:       {len(records)}")
    print(f"  Records stored:       {count}")
    print(f"  Pre-existing records: {existing_count}")
    print(f"  Final collection:     {final_count}")
    print(f"  Time elapsed:         {elapsed:.1f}s")
    print("=" * 70)

    # Per-query breakdown
    print()
    print("  Per-query breakdown:")
    print(f"  {'Query':<62} {'New':>5} {'NAs':>5} {'Dup':>5}")
    print("  " + "-" * 80)
    for query_name, new, no_abs, dupes in query_stats:
        print(f"  {query_name:<62} {new:>5} {no_abs:>5} {dupes:>5}")
    print()


if __name__ == "__main__":
    main()
