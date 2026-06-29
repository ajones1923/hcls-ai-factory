#!/usr/bin/env python3
# =============================================================================
# Cardiology Intelligence Agent - Setup Milvus Collections
# Author: Adam Jones
# Date: March 2026
# =============================================================================
"""Create and configure Milvus collections for the Cardiology Intelligence Agent.

This script initializes all 12 domain-specific vector collections used by the
Cardiology Intelligence Agent, including schema definition, IVF_FLAT indexing,
and optional seeding with curated knowledge base data.

Collections created:
    - cardio_conditions          Cardiovascular disease phenotypes and classifications
    - cardio_biomarkers          Cardiac biomarker reference ranges and clinical utility
    - cardio_drug_classes        Cardiac pharmacotherapy classes and mechanisms
    - cardio_genes               Cardiovascular genomics (genes, variants, inheritance)
    - cardio_imaging             Imaging protocols, measurements, and normal values
    - cardio_ecg                 ECG criteria and arrhythmia management
    - cardio_guidelines          ACC/AHA/ESC guideline recommendations
    - cardio_devices             FDA-cleared cardiac devices and implantables
    - cardio_hemodynamics        Hemodynamic parameters and cath-lab protocols
    - cardio_trials              Landmark cardiovascular clinical trials
    - cardio_pubmed              PubMed literature (populated via run_ingest.py)
    - cardio_clinical_trials     Active clinical trials (populated via run_ingest.py)

Usage:
    python scripts/setup_collections.py                  # create if not exists
    python scripts/setup_collections.py --drop-existing  # recreate from scratch
    python scripts/setup_collections.py --seed           # create and seed
    python scripts/setup_collections.py --list           # list current collections
"""

import argparse
import logging
import sys
import time
from pathlib import Path

# ---------------------------------------------------------------------------
# Path setup - ensure project root is importable
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pymilvus import (
    connections,
    utility,
    Collection,
    CollectionSchema,
)
from src.collections import (
    ALL_COLLECTIONS,
    INDEX_TYPE,
    METRIC_TYPE,
    NLIST,
)
from config.settings import settings

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _connect():
    """Establish Milvus connection using project settings."""
    logger.info(
        "Connecting to Milvus at %s:%s", settings.MILVUS_HOST, settings.MILVUS_PORT
    )
    connections.connect(
        alias="default",
        host=settings.MILVUS_HOST,
        port=settings.MILVUS_PORT,
    )


def _disconnect():
    """Gracefully disconnect from Milvus."""
    connections.disconnect("default")
    logger.info("Disconnected from Milvus.")


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------
def list_collections():
    """List all existing Milvus collections and their entity counts."""
    _connect()
    existing = utility.list_collections()
    cardio = [c for c in existing if c.startswith("cardio_")]
    logger.info("Found %d cardiology collections (of %d total):", len(cardio), len(existing))
    for name in sorted(cardio):
        coll = Collection(name)
        coll.flush()
        logger.info("  %-30s  %8d entities", name, coll.num_entities)
    _disconnect()


def create_collections(drop_existing: bool = False):
    """Create all 12 cardiology Milvus collections.

    Parameters
    ----------
    drop_existing : bool
        If True, drop and recreate collections that already exist.
        If False (default), skip collections that already exist.
    """
    _connect()
    created, skipped = 0, 0

    for coll_config in ALL_COLLECTIONS:
        name = coll_config.name

        # Handle existing collections
        if utility.has_collection(name):
            if drop_existing:
                logger.info("Dropping existing collection: %s", name)
                utility.drop_collection(name)
            else:
                logger.info("Collection %s already exists, skipping.", name)
                skipped += 1
                continue

        # Build schema from collection config
        schema = CollectionSchema(
            fields=coll_config.schema_fields,
            description=coll_config.description,
        )
        collection = Collection(name=name, schema=schema)

        # Create IVF_FLAT index on the embedding field
        index_params = {
            "index_type": INDEX_TYPE,
            "metric_type": METRIC_TYPE,
            "params": {"nlist": NLIST},
        }
        collection.create_index("embedding", index_params)
        logger.info("Created IVF_FLAT index on %s.embedding (nlist=%d)", name, NLIST)

        # Load collection into memory for immediate querying
        collection.load()
        created += 1
        logger.info("Created and loaded collection: %s", name)

    _disconnect()
    logger.info(
        "Setup complete. Created: %d | Skipped: %d | Total configured: %d",
        created,
        skipped,
        len(ALL_COLLECTIONS),
    )


def verify_collections():
    """Verify that all expected collections exist and are loaded."""
    _connect()
    missing = []
    for coll_config in ALL_COLLECTIONS:
        if not utility.has_collection(coll_config.name):
            missing.append(coll_config.name)
    _disconnect()

    if missing:
        logger.error("Missing collections: %s", ", ".join(missing))
        return False
    logger.info("All %d collections verified.", len(ALL_COLLECTIONS))
    return True


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Setup Cardiology Intelligence Agent Milvus collections",
    )
    parser.add_argument(
        "--drop-existing",
        action="store_true",
        help="Drop and recreate existing collections",
    )
    parser.add_argument(
        "--seed",
        action="store_true",
        help="Seed collections with curated knowledge base data after creation",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List existing cardiology collections and exit",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify all expected collections exist and exit",
    )
    args = parser.parse_args()

    if args.list:
        list_collections()
        return

    if args.verify:
        ok = verify_collections()
        sys.exit(0 if ok else 1)

    # Create collections
    start = time.time()
    create_collections(drop_existing=args.drop_existing)
    elapsed = time.time() - start
    logger.info("Collection setup took %.1fs", elapsed)

    # Optionally seed with curated data
    if args.seed:
        from scripts.seed_knowledge import seed_all

        seed_all()


if __name__ == "__main__":
    main()
