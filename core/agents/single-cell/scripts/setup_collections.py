#!/usr/bin/env python3
"""Setup Milvus collections for the Single-Cell Intelligence Agent.

Creates 12 single-cell-specific vector collections in Milvus with
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
    "sc_cell_types",
    "sc_markers",
    "sc_spatial",
    "sc_tme",
    "sc_drug_response",
    "sc_literature",
    "sc_methods",
    "sc_datasets",
    "sc_trajectories",
    "sc_pathways",
    "sc_clinical",
    "genomic_evidence",
]


def main():
    """Create all single-cell collections in Milvus."""
    logger.info("Setting up %d single-cell collections ...", len(COLLECTIONS))

    try:
        from pymilvus import connections, utility

        from config.settings import settings

        connections.connect(
            alias="default",
            host=settings.MILVUS_HOST,
            port=settings.MILVUS_PORT,
        )
        logger.info("Connected to Milvus at %s:%d", settings.MILVUS_HOST, settings.MILVUS_PORT)

        from pymilvus import Collection

        # This agent's collections.py exposes ALL_COLLECTIONS (CollectionConfig objects) but
        # no COLLECTION_SCHEMAS mapping, unlike neurology's. Build the schemas from the configs.
        from pymilvus import CollectionSchema

        from src.vector_collections import ALL_COLLECTIONS, get_collection_config

        COLLECTION_SCHEMAS = {
            cfg.name: CollectionSchema(fields=cfg.schema_fields, description=cfg.description)
            for cfg in ALL_COLLECTIONS
        }

        existing = utility.list_collections()
        created = 0
        for name in COLLECTIONS:
            if name in existing:
                logger.info("  [exists] %s", name)
                continue
            schema = COLLECTION_SCHEMAS.get(name)
            if schema is None:
                logger.warning("  [skip]   %s -- no schema in src/vector_collections.py", name)
                continue
            coll = Collection(name=name, schema=schema)
            try:
                idx = get_collection_config(name).index_params
            except Exception:
                idx = {"metric_type": "COSINE", "index_type": "IVF_FLAT",
                       "params": {"nlist": 128}}
            coll.create_index(field_name="embedding", index_params=idx)
            coll.load()
            created += 1
            logger.info("  [create] %s (indexed, loaded)", name)

        logger.info("Collection setup complete: %d created, %d already present.",
                    created, len(COLLECTIONS) - created)

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
