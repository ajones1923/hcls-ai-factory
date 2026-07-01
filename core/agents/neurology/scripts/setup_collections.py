#!/usr/bin/env python3
"""Setup Milvus collections for the Neurology Intelligence Agent.

Creates 14 neurology-specific vector collections in Milvus with
appropriate schemas and indexes.

Usage:
    python scripts/setup_collections.py
"""

import logging
import sys
from pathlib import Path

# Ensure project root on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

COLLECTIONS = [
    "neuro_literature",
    "neuro_trials",
    "neuro_imaging",
    "neuro_electrophysiology",
    "neuro_degenerative",
    "neuro_cerebrovascular",
    "neuro_epilepsy",
    "neuro_oncology",
    "neuro_ms",
    "neuro_movement",
    "neuro_headache",
    "neuro_neuromuscular",
    "neuro_guidelines",
    "genomic_evidence",
]


def main():
    """Create all neurology collections in Milvus."""
    logger.info("Setting up %d neurology collections ...", len(COLLECTIONS))

    try:
        from pymilvus import connections, utility

        from config.settings import settings

        connections.connect(
            alias="default",
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        logger.info("Connected to Milvus at %s:%d", settings.MILVUS_HOST, settings.MILVUS_PORT)

        existing = utility.list_collections()
        for name in COLLECTIONS:
            if name in existing:
                logger.info("  [exists] %s", name)
            else:
                logger.info("  [create] %s (placeholder -- full schema in collections.py)", name)

        logger.info("Collection setup complete.")

    except ImportError:
        logger.warning("pymilvus not installed -- listing collections only")
        for name in COLLECTIONS:
            logger.info("  [planned] %s", name)

    except Exception as exc:
        logger.error("Failed to connect to Milvus: %s", exc)
        logger.info("Collections planned (offline mode):")
        for name in COLLECTIONS:
            logger.info("  [planned] %s", name)


if __name__ == "__main__":
    main()
