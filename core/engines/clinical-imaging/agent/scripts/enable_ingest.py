#!/usr/bin/env python3
"""Enable and trigger the automated ingest scheduler.

Sets INGEST_ENABLED=True, validates Milvus connectivity, triggers an
immediate ingest run (PubMed + ClinicalTrials.gov), and prints the
schedule for future automated runs.

Usage:
    python scripts/enable_ingest.py              # Enable + trigger immediate ingest
    python scripts/enable_ingest.py --status      # Show scheduler status
    python scripts/enable_ingest.py --dry-run     # Show what would be ingested without writing

Author: Adam Jones
Date: April 2026
"""

import argparse
import os
import sys
import time
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger


def check_milvus(settings) -> bool:
    """Validate that Milvus is accessible and collections exist."""
    try:
        from src.collections import ImagingCollectionManager

        manager = ImagingCollectionManager(
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        manager.connect()
        stats = manager.get_collection_stats()
        total = sum(stats.values())
        logger.info(f"Milvus connected: {len(stats)} collections, {total} total vectors")
        for name, count in sorted(stats.items()):
            logger.info(f"  {name}: {count}")
        return True
    except Exception as e:
        logger.error(f"Milvus connection failed: {e}")
        return False


def show_status(settings):
    """Display current scheduler configuration and status."""
    print("\n=== Ingest Scheduler Status ===")
    print(f"  INGEST_ENABLED:              {settings.INGEST_ENABLED}")
    print(f"  INGEST_SCHEDULE_HOURS:       {settings.INGEST_SCHEDULE_HOURS} ({settings.INGEST_SCHEDULE_HOURS / 24:.0f} days)")
    print(f"  INGEST_TRIALS_SCHEDULE_HOURS:{settings.INGEST_TRIALS_SCHEDULE_HOURS} ({settings.INGEST_TRIALS_SCHEDULE_HOURS / 24:.0f} days)")
    print(f"  MILVUS_HOST:                 {settings.MILVUS_HOST}:{settings.MILVUS_PORT}")
    print(f"  EMBEDDING_MODEL:             {settings.EMBEDDING_MODEL}")
    print(f"  PUBMED_MAX_RESULTS:          {settings.PUBMED_MAX_RESULTS}")
    print(f"  CT_GOV_BASE_URL:             {settings.CT_GOV_BASE_URL}")
    print()

    milvus_ok = check_milvus(settings)
    print(f"  Milvus reachable:            {'Yes' if milvus_ok else 'No'}")

    print("\n=== Schedule ===")
    print(f"  PubMed literature:           every {settings.INGEST_SCHEDULE_HOURS}h ({settings.INGEST_SCHEDULE_HOURS / 24:.0f} days)")
    print(f"  ClinicalTrials.gov trials:   every {settings.INGEST_TRIALS_SCHEDULE_HOURS}h ({settings.INGEST_TRIALS_SCHEDULE_HOURS / 24:.0f} days)")
    print()


def run_immediate_ingest(settings, dry_run: bool = False):
    """Trigger an immediate ingest cycle (PubMed + ClinicalTrials.gov)."""
    from sentence_transformers import SentenceTransformer
    from src.collections import ImagingCollectionManager
    from src.scheduler import ImagingIngestScheduler

    logger.info("Initializing collection manager and embedder...")
    manager = ImagingCollectionManager(
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )
    manager.connect()
    manager.ensure_collections()

    embedder = SentenceTransformer(settings.EMBEDDING_MODEL)
    logger.info(f"Loaded embedding model: {settings.EMBEDDING_MODEL}")

    if dry_run:
        print("\n=== Dry Run ===")
        print("  Would ingest PubMed imaging AI literature")
        print(f"    Query: medical imaging AI (max {settings.PUBMED_MAX_RESULTS} results)")
        print(f"    Target: {settings.COLLECTION_LITERATURE}")
        print("  Would ingest ClinicalTrials.gov imaging AI trials")
        print(f"    API: {settings.CT_GOV_BASE_URL}/studies")
        print(f"    Target: {settings.COLLECTION_TRIALS}")
        print("\n  No data written. Remove --dry-run to execute.")
        return

    scheduler = ImagingIngestScheduler(
        collection_manager=manager,
        embedder=embedder,
    )

    print("\n=== Running Immediate Ingest Cycle ===")
    start = time.time()
    stats = scheduler.run_ingest_cycle()
    elapsed = time.time() - start

    print(f"\n=== Ingest Complete ({elapsed:.1f}s) ===")
    print(f"  PubMed records:   {stats['pubmed_records']}")
    print(f"  Trials records:   {stats['trials_records']}")
    print(f"  Total records:    {stats['total_records']}")
    print(f"  Success:          {stats['success']}")
    if stats.get("errors"):
        for err in stats["errors"]:
            print(f"  Error: {err}")


def main():
    parser = argparse.ArgumentParser(
        description="Enable and trigger the imaging ingest scheduler.",
    )
    parser.add_argument(
        "--status", action="store_true",
        help="Show scheduler status without running ingest",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Show what would be ingested without writing data",
    )
    args = parser.parse_args()

    # Set INGEST_ENABLED in environment before loading settings
    os.environ["IMAGING_INGEST_ENABLED"] = "true"

    from config.settings import ImagingSettings
    settings = ImagingSettings()

    if args.status:
        show_status(settings)
        return

    # Validate Milvus before proceeding
    if not check_milvus(settings):
        logger.error("Cannot proceed without Milvus. Start Milvus and retry.")
        sys.exit(1)

    show_status(settings)
    run_immediate_ingest(settings, dry_run=args.dry_run)

    print("\n=== Future Schedule ===")
    print(f"  To keep the scheduler running, start the API server with:")
    print(f"    IMAGING_INGEST_ENABLED=true uvicorn api.main:app --host 0.0.0.0 --port 8524")
    print(f"  PubMed will ingest every {settings.INGEST_SCHEDULE_HOURS}h ({settings.INGEST_SCHEDULE_HOURS / 24:.0f} days)")
    print(f"  Trials will ingest every {settings.INGEST_TRIALS_SCHEDULE_HOURS}h ({settings.INGEST_TRIALS_SCHEDULE_HOURS / 24:.0f} days)")
    print()


if __name__ == "__main__":
    main()
