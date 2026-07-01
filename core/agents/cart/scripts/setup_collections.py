#!/usr/bin/env python3
"""Create all CAR-T Milvus collections and seed with FDA-approved constructs.

Usage:
    python scripts/setup_collections.py [--drop-existing] [--seed-constructs]

Options:
    --drop-existing    Drop and recreate all collections
    --seed-constructs  Seed cart_constructs with 6 FDA-approved products
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from src.collections import CARTCollectionManager


def main():
    parser = argparse.ArgumentParser(description="Setup CAR-T Milvus collections")
    parser.add_argument("--drop-existing", action="store_true",
                       help="Drop and recreate all collections")
    parser.add_argument("--seed-constructs", action="store_true",
                       help="Seed cart_constructs with FDA-approved products")
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    args = parser.parse_args()

    # Connect to Milvus
    manager = CARTCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Create all collections
    logger.info("Creating all CAR-T collections...")
    manager.create_all_collections(drop_existing=args.drop_existing)

    # Show stats
    stats = manager.get_collection_stats()
    logger.info("Collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    # Seed FDA constructs if requested
    if args.seed_constructs:
        logger.info("Seeding FDA-approved constructs...")
        try:
            from sentence_transformers import SentenceTransformer

            model = SentenceTransformer("BAAI/bge-small-en-v1.5")

            from src.ingest.construct_parser import ConstructIngestPipeline

            class SimpleEmbedder:
                def __init__(self, model):
                    self._model = model

                def encode(self, texts):
                    return self._model.encode(texts).tolist()

            embedder = SimpleEmbedder(model)
            pipeline = ConstructIngestPipeline(manager, embedder)
            count = pipeline.run(include_fda_seed=True)
            logger.info(f"Seeded {count} FDA-approved constructs")
        except ImportError as e:
            logger.error(f"Cannot seed constructs - missing dependency: {e}")
        except Exception as e:
            logger.error(f"Error seeding constructs: {e}")

    # Final stats
    stats = manager.get_collection_stats()
    logger.info("Final collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    manager.disconnect()
    logger.info("Done!")


if __name__ == "__main__":
    main()
