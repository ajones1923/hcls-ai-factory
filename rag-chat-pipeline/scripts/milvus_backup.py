#!/usr/bin/env python3
"""
Milvus Backup/Restore for genomic_evidence collection.

Exports all non-vector fields to JSON. On restore, re-generates embeddings
from text_summary using BGE-small-en-v1.5 (same as original ingest).

Usage:
    python scripts/milvus_backup.py backup --output data/backups/genomic_evidence.json
    python scripts/milvus_backup.py restore --input data/backups/genomic_evidence.json
"""

import argparse
import json
import os
import sys
from pathlib import Path

# Add parent for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from loguru import logger

# Fields to export (all non-vector fields)
EXPORT_FIELDS = [
    "id", "chrom", "pos", "ref", "alt", "qual", "gene",
    "consequence", "impact", "genotype", "text_summary",
    "clinical_significance", "rsid", "disease_associations",
    "am_pathogenicity", "am_class",
]

COLLECTION_NAME = "genomic_evidence"


def backup_collection(output_path: str, batch_size: int = 10000):
    """Export all entities from genomic_evidence collection to JSON."""
    from pymilvus import connections, Collection, utility

    host = os.environ.get("MILVUS_HOST", "localhost")
    port = int(os.environ.get("MILVUS_PORT", "19530"))

    logger.info(f"Connecting to Milvus at {host}:{port}")
    connections.connect("default", host=host, port=port, timeout=10)

    if not utility.has_collection(COLLECTION_NAME):
        logger.error(f"Collection '{COLLECTION_NAME}' not found")
        connections.disconnect("default")
        sys.exit(1)

    collection = Collection(COLLECTION_NAME)
    collection.load()
    total = collection.num_entities
    logger.info(f"Collection has {total:,} entities")

    all_entities = []
    offset = 0
    while True:
        results = collection.query(
            expr="id >= 0",
            output_fields=EXPORT_FIELDS,
            limit=batch_size,
            offset=offset,
        )
        if not results:
            break
        all_entities.extend(results)
        offset += len(results)
        logger.info(f"  Exported {offset:,} / {total:,} entities")

    # Ensure output directory exists
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(all_entities, f, default=str)

    connections.disconnect("default")
    logger.info(f"Backed up {len(all_entities):,} entities to {output_path}")


def restore_collection(input_path: str, batch_size: int = 1000):
    """Re-ingest entities from backup file, re-generating embeddings."""
    from pymilvus import connections

    host = os.environ.get("MILVUS_HOST", "localhost")
    port = int(os.environ.get("MILVUS_PORT", "19530"))

    logger.info(f"Loading backup from {input_path}")
    with open(input_path) as f:
        entities = json.load(f)
    logger.info(f"Loaded {len(entities):,} entities")

    # Initialize embedder (same model as original ingest)
    try:
        from sentence_transformers import SentenceTransformer
    except ImportError:
        logger.error("sentence-transformers is required for restore. Install with: pip install sentence-transformers")
        sys.exit(1)

    model_name = os.environ.get("EMBEDDING_MODEL", "BAAI/bge-small-en-v1.5")
    logger.info(f"Loading embedding model: {model_name}")
    embedder = SentenceTransformer(model_name)

    logger.info(f"Connecting to Milvus at {host}:{port}")
    connections.connect("default", host=host, port=port, timeout=10)

    # Use existing MilvusClient to create/get collection
    from src.milvus_client import MilvusClient
    client = MilvusClient(host=host, port=port)
    collection = client.get_collection()

    # Insert in batches
    inserted = 0
    for i in range(0, len(entities), batch_size):
        batch = entities[i : i + batch_size]

        # Re-generate embeddings from text_summary
        texts = [e.get("text_summary", "") for e in batch]
        embeddings = embedder.encode(texts, normalize_embeddings=True).tolist()

        # Build insert data (exclude 'id' — Milvus auto-generates)
        insert_data = [
            [e.get("chrom", "") for e in batch],
            [int(e.get("pos", 0)) for e in batch],
            [e.get("ref", "") for e in batch],
            [e.get("alt", "") for e in batch],
            [float(e.get("qual", 0)) for e in batch],
            [e.get("gene", "") for e in batch],
            [e.get("consequence", "") for e in batch],
            [e.get("impact", "") for e in batch],
            [e.get("genotype", "") for e in batch],
            [e.get("text_summary", "") for e in batch],
            [e.get("clinical_significance", "") for e in batch],
            [e.get("rsid", "") for e in batch],
            [e.get("disease_associations", "") for e in batch],
            [float(e.get("am_pathogenicity", 0)) for e in batch],
            [e.get("am_class", "") for e in batch],
            embeddings,
        ]

        collection.insert(insert_data)
        inserted += len(batch)
        logger.info(f"  Inserted {inserted:,} / {len(entities):,} entities")

    collection.flush()
    connections.disconnect("default")
    logger.info(f"Restore complete: {inserted:,} entities inserted")


def main():
    parser = argparse.ArgumentParser(description="Milvus backup/restore for genomic_evidence")
    sub = parser.add_subparsers(dest="command", required=True)

    bp = sub.add_parser("backup", help="Export collection to JSON")
    bp.add_argument("--output", required=True, help="Output JSON file path")
    bp.add_argument("--batch-size", type=int, default=10000)

    rp = sub.add_parser("restore", help="Re-ingest from JSON backup")
    rp.add_argument("--input", required=True, help="Input JSON file path")
    rp.add_argument("--batch-size", type=int, default=1000)

    args = parser.parse_args()

    if args.command == "backup":
        backup_collection(args.output, args.batch_size)
    elif args.command == "restore":
        restore_collection(args.input, args.batch_size)


if __name__ == "__main__":
    main()
