#!/usr/bin/env python3
"""Create all 10 Imaging Intelligence Agent Milvus collections.

Usage:
    python scripts/setup_collections.py
    python scripts/setup_collections.py --drop-existing
    python scripts/setup_collections.py --host 10.0.0.1 --port 19530

Options:
    --drop-existing    Drop and recreate all collections
    --host             Milvus host (default: localhost)
    --port             Milvus port (default: 19530)

Author: Adam Jones
Date: February 2026
"""

import argparse
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from src.collections import ImagingCollectionManager


def main():
    parser = argparse.ArgumentParser(
        description="Create all Imaging Intelligence Agent Milvus collections"
    )
    parser.add_argument(
        "--drop-existing",
        action="store_true",
        help="Drop and recreate all collections",
    )
    parser.add_argument(
        "--host", default=None, help="Milvus host (default: localhost)"
    )
    parser.add_argument(
        "--port", type=int, default=None, help="Milvus port (default: 19530)"
    )
    args = parser.parse_args()

    print("=" * 65)
    print("  Imaging Intelligence Agent — Collection Setup")
    print("=" * 65)
    print(f"  Drop existing: {args.drop_existing}")
    print()

    # Connect to Milvus
    manager = ImagingCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Create all 10 imaging collections
    logger.info("Creating all imaging collections...")
    collections = manager.create_all_collections(drop_existing=args.drop_existing)
    logger.info(f"Created/verified {len(collections)} collections")

    # Show stats
    stats = manager.get_collection_stats()
    logger.info("Collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    created_count = len(collections)
    manager.disconnect()

    print()
    print("=" * 65)
    print(f"  Setup complete: {created_count} collections ready")
    print("=" * 65)


if __name__ == "__main__":
    main()
