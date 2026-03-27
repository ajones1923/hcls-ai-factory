#!/usr/bin/env python3
"""Create all 15 PGx Milvus collections with proper schemas.

Connects to Milvus, creates (or recreates) all pharmacogenomics collections
defined in src/collections.py, optionally drops existing data, and seeds
reference data from data/reference/*.json files.

Usage:
    python scripts/setup_collections.py
    python scripts/setup_collections.py --drop-existing
    python scripts/setup_collections.py --drop-existing --seed
    python scripts/setup_collections.py --host milvus-standalone --port 19530

Options:
    --drop-existing    Drop and recreate all collections (WARNING: deletes data)
    --seed             Load reference data from data/reference/*.json after creation
    --host HOST        Milvus server hostname (default: PGX_MILVUS_HOST or localhost)
    --port PORT        Milvus server port (default: PGX_MILVUS_PORT or 19530)

Author: Adam Jones
Date: March 2026
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from loguru import logger

from src.collections import COLLECTION_SCHEMAS, PGxCollectionManager


# ═══════════════════════════════════════════════════════════════════════
# REFERENCE DATA MAPPING — which JSON file seeds which collection
# ═══════════════════════════════════════════════════════════════════════

REFERENCE_FILE_MAP: Dict[str, str] = {
    "pgx_gene_reference": "gene_reference.json",
    "pgx_drug_guidelines": "drug_guidelines.json",
    "pgx_drug_interactions": "drug_interactions.json",
    "pgx_hla_hypersensitivity": "hla_hypersensitivity.json",
    "pgx_phenoconversion": "phenoconversion.json",
    "pgx_dosing_algorithms": "dosing_algorithms.json",
    "pgx_clinical_evidence": "clinical_evidence.json",
    "pgx_population_data": "population_data.json",
    "pgx_clinical_trials": "clinical_trials.json",
    "pgx_fda_labels": "fda_labels.json",
    "pgx_drug_alternatives": "drug_alternatives.json",
    "pgx_patient_profiles": "patient_profiles.json",
    "pgx_implementation": "implementation.json",
    "pgx_education": "education.json",
}


def _load_embedding_model():
    """Load the BGE-small-en-v1.5 sentence transformer."""
    from sentence_transformers import SentenceTransformer

    logger.info("Loading embedding model: BAAI/bge-small-en-v1.5")
    model = SentenceTransformer("BAAI/bge-small-en-v1.5")
    logger.info("Embedding model loaded (384-dim)")
    return model


def _build_text_for_embedding(record: Dict[str, Any]) -> str:
    """Build a text string from a record for embedding.

    Concatenates all string fields into a single passage suitable for
    BGE-small-en-v1.5 embedding. Skips 'id' and 'embedding' fields.
    """
    parts = []
    for key, value in record.items():
        if key in ("id", "embedding"):
            continue
        if isinstance(value, str) and value.strip():
            parts.append(f"{key}: {value}")
        elif isinstance(value, (int, float)):
            parts.append(f"{key}: {value}")
        elif isinstance(value, list):
            parts.append(f"{key}: {', '.join(str(v) for v in value)}")
    return ". ".join(parts)


def _get_collection_text_fields(collection_name: str) -> List[str]:
    """Return the non-embedding, non-id field names for a collection."""
    if collection_name not in COLLECTION_SCHEMAS:
        return []
    schema = COLLECTION_SCHEMAS[collection_name]
    return [
        field.name
        for field in schema.fields
        if field.name not in ("id", "embedding")
    ]


def _coerce_record(
    record: Dict[str, Any],
    collection_name: str,
) -> Dict[str, Any]:
    """Coerce a JSON record to match the Milvus schema field types.

    Ensures all expected fields are present with the correct type, filling
    missing fields with safe defaults (empty string, 0.0, etc.).
    """
    if collection_name not in COLLECTION_SCHEMAS:
        return record

    schema = COLLECTION_SCHEMAS[collection_name]
    coerced: Dict[str, Any] = {}

    for field in schema.fields:
        name = field.name
        if name == "embedding":
            continue  # handled separately

        value = record.get(name)

        # Determine type from dtype enum name
        dtype_name = field.dtype.name if hasattr(field.dtype, "name") else str(field.dtype)

        if "VARCHAR" in dtype_name:
            if value is None:
                coerced[name] = ""
            elif isinstance(value, list):
                coerced[name] = ", ".join(str(v) for v in value)
            else:
                coerced[name] = str(value)
                # Truncate to max_length if defined
                max_len = getattr(field, "max_length", None)
                if max_len and len(coerced[name].encode("utf-8")) > max_len:
                    encoded = coerced[name].encode("utf-8")[:max_len]
                    coerced[name] = encoded.decode("utf-8", errors="ignore")
        elif "FLOAT" in dtype_name or "DOUBLE" in dtype_name:
            coerced[name] = float(value) if value is not None else 0.0
        elif "INT" in dtype_name:
            coerced[name] = int(value) if value is not None else 0
        elif "BOOL" in dtype_name:
            coerced[name] = bool(value) if value is not None else False
        else:
            coerced[name] = value if value is not None else ""

    return coerced


def seed_collection(
    manager: PGxCollectionManager,
    collection_name: str,
    records: List[Dict[str, Any]],
    model,
    batch_size: int = 32,
) -> int:
    """Embed and insert records into a single collection.

    Args:
        manager: Connected PGxCollectionManager.
        collection_name: Target Milvus collection name.
        records: List of dicts, each representing one record.
        model: SentenceTransformer model for embedding.
        batch_size: Number of records to embed and insert at once.

    Returns:
        Number of records inserted.
    """
    if not records:
        logger.warning(f"  No records to seed for {collection_name}")
        return 0

    total_inserted = 0

    for i in range(0, len(records), batch_size):
        batch = records[i : i + batch_size]

        # Build embedding texts
        texts = [_build_text_for_embedding(r) for r in batch]
        embeddings = model.encode(texts).tolist()

        # Coerce and attach embeddings
        rows = []
        for record, emb in zip(batch, embeddings):
            coerced = _coerce_record(record, collection_name)
            coerced["embedding"] = emb
            rows.append(coerced)

        # Insert batch
        try:
            manager.insert_batch(collection_name, rows)
            total_inserted += len(rows)
        except Exception as e:
            logger.error(f"  Failed to insert batch {i // batch_size + 1} "
                         f"into {collection_name}: {e}")

    return total_inserted


def seed_from_reference_files(
    manager: PGxCollectionManager,
    reference_dir: Path,
    model,
) -> Dict[str, int]:
    """Load all reference JSON files and seed corresponding collections.

    Args:
        manager: Connected PGxCollectionManager.
        reference_dir: Path to data/reference/ directory.
        model: SentenceTransformer model for embedding.

    Returns:
        Dict mapping collection name -> number of records seeded.
    """
    results: Dict[str, int] = {}

    for collection_name, filename in REFERENCE_FILE_MAP.items():
        filepath = reference_dir / filename
        if not filepath.exists():
            logger.warning(f"  Reference file not found: {filepath}")
            results[collection_name] = 0
            continue

        logger.info(f"  Loading {filepath.name} -> {collection_name}")
        with open(filepath, "r") as f:
            records = json.load(f)

        if not isinstance(records, list):
            logger.error(f"  {filepath.name}: expected a JSON array, got {type(records).__name__}")
            results[collection_name] = 0
            continue

        count = seed_collection(manager, collection_name, records, model)
        results[collection_name] = count
        logger.info(f"  {collection_name}: seeded {count} records")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Create all 15 PGx Milvus collections"
    )
    parser.add_argument(
        "--drop-existing",
        action="store_true",
        help="Drop and recreate all collections (WARNING: deletes data)",
    )
    parser.add_argument(
        "--seed",
        action="store_true",
        help="Load reference data from data/reference/*.json after creation",
    )
    parser.add_argument("--host", default=None, help="Milvus host")
    parser.add_argument("--port", type=int, default=None, help="Milvus port")
    args = parser.parse_args()

    # ── Connect ──
    manager = PGxCollectionManager(host=args.host, port=args.port)
    manager.connect()

    # ── Create collections ──
    logger.info("=" * 65)
    logger.info("  PGx Intelligence Agent — Collection Setup")
    logger.info("=" * 65)

    start = time.time()

    logger.info(f"Creating all {len(COLLECTION_SCHEMAS)} PGx collections "
                f"(drop_existing={args.drop_existing})")
    manager.create_all_collections(drop_existing=args.drop_existing)

    elapsed = time.time() - start
    logger.info(f"Collections created in {elapsed:.1f}s")

    # ── Show initial stats ──
    stats = manager.get_collection_stats()
    logger.info("Collection stats after creation:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    # ── Seed reference data ──
    if args.seed:
        logger.info("")
        logger.info("Seeding reference data from data/reference/*.json ...")

        try:
            model = _load_embedding_model()
        except ImportError as e:
            logger.error(f"Cannot load embedding model: {e}")
            logger.error("Install: pip install sentence-transformers")
            manager.disconnect()
            sys.exit(1)

        reference_dir = PROJECT_ROOT / "data" / "reference"
        if not reference_dir.exists():
            logger.error(f"Reference directory not found: {reference_dir}")
            manager.disconnect()
            sys.exit(1)

        seed_results = seed_from_reference_files(manager, reference_dir, model)

        total_seeded = sum(seed_results.values())
        logger.info(f"Total records seeded: {total_seeded:,}")

    # ── Final stats ──
    stats = manager.get_collection_stats()
    logger.info("")
    logger.info("Final collection stats:")
    for name, count in stats.items():
        logger.info(f"  {name}: {count:,} records")

    manager.disconnect()
    logger.info("Setup complete.")


if __name__ == "__main__":
    main()
