#!/usr/bin/env python3
"""Create all Biomarker Intelligence Milvus collections.

Usage:
    python scripts/setup_collections.py [--drop-existing]

Options:
    --drop-existing    Drop and recreate all collections

Author: Adam Jones
Date: March 2026
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from src.collections import BiomarkerCollectionManager


def main():
    parser = argparse.ArgumentParser(description="Setup Biomarker Intelligence Milvus collections")
    parser.add_argument("--drop-existing", action="store_true",
                        help="Drop and recreate all collections")
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    args = parser.parse_args()

    # Connect to Milvus
    manager = BiomarkerCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Create all collections
    logger.info("Creating all Biomarker Intelligence collections...")
    manager.create_all_collections(drop_existing=args.drop_existing)

    # Show stats
    stats = manager.get_collection_stats()
    logger.info("Collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    total = sum(stats.values())
    logger.info(f"Total: {total:,} records across {len(stats)} collections")

    manager.disconnect()
    logger.info("Done!")


if __name__ == "__main__":
    main()
