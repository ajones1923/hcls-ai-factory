#!/usr/bin/env python3
"""Seed all 10 Biomarker Intelligence Milvus collections from JSON reference files.

Reads JSON seed files from data/reference/, embeds text_chunk fields with
BGE-small-en-v1.5, and inserts into the corresponding Milvus collections.

Usage:
    python scripts/seed_all.py [--host HOST] [--port PORT] [--batch-size 32]

Author: Adam Jones
Date: March 2026
"""

import argparse
import json
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger
from tqdm import tqdm

from src.collections import COLLECTION_MODELS, COLLECTION_SCHEMAS, BiomarkerCollectionManager

# Mapping of collection names to their JSON seed filenames
SEED_FILES = {
    "biomarker_reference": "biomarker_reference.json",
    "biomarker_genetic_variants": "biomarker_genetic_variants.json",
    "biomarker_pgx_rules": "biomarker_pgx_rules.json",
    "biomarker_disease_trajectories": "biomarker_disease_trajectories.json",
    "biomarker_clinical_evidence": "biomarker_clinical_evidence.json",
    "biomarker_nutrition": "biomarker_nutrition.json",
    "biomarker_drug_interactions": "biomarker_drug_interactions.json",
    "biomarker_aging_markers": "biomarker_aging_markers.json",
    "biomarker_genotype_adjustments": "biomarker_genotype_adjustments.json",
    "biomarker_monitoring": "biomarker_monitoring.json",
    "biomarker_critical_values": "biomarker_critical_values.json",
    "biomarker_discordance_rules": "biomarker_discordance_rules.json",
    "biomarker_aj_carrier_screening": "biomarker_aj_carrier_screening.json",
    "genomic_evidence": "biomarker_genomic_evidence.json",
}


def load_seed_file(reference_dir: Path, filename: str):
    """Load a JSON seed file and return the records list."""
    filepath = reference_dir / filename
    if not filepath.exists():
        logger.warning(f"Seed file not found: {filepath}")
        return []
    with open(filepath, "r") as f:
        data = json.load(f)
    # Support both list format and dict with "records" key
    if isinstance(data, list):
        return data
    elif isinstance(data, dict) and "records" in data:
        return data["records"]
    else:
        logger.warning(f"Unexpected format in {filepath}")
        return []


def embed_and_insert(
    manager: BiomarkerCollectionManager,
    embedder,
    collection_name: str,
    records: list,
    batch_size: int = 32,
) -> int:
    """Embed text_chunk fields and insert records into a Milvus collection.

    Args:
        manager: BiomarkerCollectionManager instance.
        embedder: SentenceTransformer model for encoding.
        collection_name: Target collection name.
        records: List of record dicts from the seed file.
        batch_size: Number of records to embed and insert at once.

    Returns:
        Total number of records inserted.
    """
    if not records:
        logger.info(f"  {collection_name}: no records to insert")
        return 0

    model_cls = COLLECTION_MODELS.get(collection_name)
    schema = COLLECTION_SCHEMAS.get(collection_name)
    if not schema:
        logger.error(f"  No schema found for {collection_name}")
        return 0

    # Get expected field names from schema (excluding 'embedding')
    field_names = [f.name for f in schema.fields if f.name != "embedding"]

    total_inserted = 0
    for i in tqdm(range(0, len(records), batch_size), desc=f"  {collection_name}"):
        batch = records[i : i + batch_size]

        # Pre-process: convert list fields to comma-separated strings for Milvus VARCHAR
        for rec in batch:
            for key, val in list(rec.items()):
                if isinstance(val, list):
                    rec[key] = ", ".join(str(v) for v in val)

        # Extract text for embedding
        texts = []
        for rec in batch:
            if model_cls:
                try:
                    model_instance = model_cls(**rec)
                    texts.append(model_instance.to_embedding_text())
                except Exception:
                    texts.append(rec.get("text_chunk", rec.get("text_summary", "")))
            else:
                texts.append(rec.get("text_chunk", rec.get("text_summary", "")))

        # Generate embeddings
        prefix = "Represent this sentence for searching relevant passages: "
        prefixed_texts = [prefix + t for t in texts]
        embeddings = embedder.encode(prefixed_texts).tolist()

        # Build insertion records with only fields present in the schema
        insert_records = []
        for j, rec in enumerate(batch):
            insert_rec = {"embedding": embeddings[j]}
            for field_name in field_names:
                if field_name in rec and rec[field_name] is not None:
                    insert_rec[field_name] = rec[field_name]
                elif field_name == "id":
                    insert_rec["id"] = rec.get("id", f"{collection_name}_{i + j}")
                else:
                    # Provide default values for missing fields
                    field_schema = next(
                        (f for f in schema.fields if f.name == field_name), None
                    )
                    if field_schema:
                        from pymilvus import DataType
                        if field_schema.dtype == DataType.VARCHAR:
                            insert_rec[field_name] = ""
                        elif field_schema.dtype == DataType.INT64:
                            insert_rec[field_name] = 0
                        elif field_schema.dtype == DataType.FLOAT:
                            insert_rec[field_name] = 0.0
                        else:
                            insert_rec[field_name] = ""
            insert_records.append(insert_rec)

        # Insert batch
        count = manager.insert_batch(collection_name, insert_records)
        total_inserted += count

    return total_inserted


def main():
    parser = argparse.ArgumentParser(description="Seed all Biomarker Intelligence collections")
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    parser.add_argument("--batch-size", type=int, default=32, help="Embedding batch size")
    args = parser.parse_args()

    reference_dir = PROJECT_ROOT / "data" / "reference"
    if not reference_dir.exists():
        logger.error(f"Reference directory not found: {reference_dir}")
        sys.exit(1)

    # Connect to Milvus
    manager = BiomarkerCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # Load embedding model
    logger.info("Loading BGE-small-en-v1.5 embedding model...")
    from sentence_transformers import SentenceTransformer
    embedder = SentenceTransformer("BAAI/bge-small-en-v1.5")
    logger.info("Embedding model loaded")

    # Seed each collection
    total_all = 0
    for collection_name, filename in SEED_FILES.items():
        logger.info(f"Seeding {collection_name} from {filename}...")
        records = load_seed_file(reference_dir, filename)
        if records:
            count = embed_and_insert(
                manager, embedder, collection_name, records,
                batch_size=args.batch_size,
            )
            total_all += count
            logger.info(f"  Inserted {count} records into {collection_name}")
        else:
            logger.warning(f"  No seed data found for {collection_name}")

    # Final stats
    logger.info("=" * 60)
    logger.info("Seeding complete!")
    stats = manager.get_collection_stats()
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")
    logger.info(f"Total inserted this run: {total_all:,}")

    manager.disconnect()


if __name__ == "__main__":
    main()
