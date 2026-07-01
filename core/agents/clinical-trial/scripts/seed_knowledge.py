"""Seed Milvus collections with curated knowledge for the Clinical Trial Intelligence Agent.

Populates the trial_protocols and trial_regulatory collections with
landmark trials and regulatory milestones from the ingest parsers.

Usage:
    python scripts/seed_knowledge.py [--host localhost] [--port 19530]

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List

# Ensure project root is on path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.ingest.clinicaltrials_parser import ClinicalTrialsParser, LANDMARK_TRIALS
from src.ingest.regulatory_parser import RegulatoryParser, REGULATORY_MILESTONES
from src.ingest.base import IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA COUNTS
# ===================================================================

EXPECTED_LANDMARK_TRIALS = len(LANDMARK_TRIALS)
EXPECTED_REGULATORY_MILESTONES = len(REGULATORY_MILESTONES)


# ===================================================================
# INSERT HELPER
# ===================================================================


def _insert_records(
    collection_name: str,
    records: List[IngestRecord],
    text_field: str = "text",
) -> int:
    """Generate embeddings and insert IngestRecord objects into a Milvus collection.

    Degrades gracefully: if pymilvus or sentence_transformers are not
    installed, or if Milvus is unreachable, logs a warning and returns
    the record count (as if it were a dry run).

    Parameters
    ----------
    collection_name : str
        Target Milvus collection name.
    records : list[IngestRecord]
        Records to insert.  Each must have a ``.text`` attribute.
    text_field : str
        Attribute name whose value is used to produce the embedding vector.

    Returns
    -------
    int
        Number of records inserted (or that would have been inserted on
        graceful degradation).
    """
    if not records:
        logger.info("No records to insert into '%s'.", collection_name)
        return 0

    # --- load optional dependencies ---
    try:
        from pymilvus import MilvusClient  # noqa: F811
    except ImportError:
        logger.warning(
            "pymilvus is not installed – skipping Milvus insert for %s "
            "(%d records would have been inserted).",
            collection_name,
            len(records),
        )
        return len(records)

    try:
        from sentence_transformers import SentenceTransformer
    except ImportError:
        logger.warning(
            "sentence_transformers is not installed – skipping Milvus "
            "insert for %s (%d records would have been inserted).",
            collection_name,
            len(records),
        )
        return len(records)

    # --- connect to Milvus ---
    try:
        from config.settings import settings

        client = MilvusClient(
            uri=f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}",
        )
    except Exception as exc:
        logger.warning(
            "Could not connect to Milvus for %s: %s – "
            "treating as dry run (%d records).",
            collection_name,
            exc,
            len(records),
        )
        return len(records)

    # --- generate embeddings ---
    try:
        model = SentenceTransformer("BAAI/bge-small-en-v1.5")
        texts = [getattr(r, text_field, "") for r in records]
        embeddings = model.encode(texts, show_progress_bar=False).tolist()
    except Exception as exc:
        logger.warning(
            "Embedding generation failed for %s: %s – "
            "treating as dry run (%d records).",
            collection_name,
            exc,
            len(records),
        )
        return len(records)

    # --- discover valid schema fields ---
    try:
        desc = client.describe_collection(collection_name)
        valid_fields = {f["name"] for f in desc.get("fields", [])}
    except Exception:
        valid_fields = None

    # --- insert ---
    try:
        insert_data = []
        for record, emb in zip(records, embeddings):
            entry = {
                "embedding": emb,
                "text_content": record.text[:8192],
            }
            # Include metadata fields that exist in the collection schema
            for k, v in record.metadata.items():
                val = str(v)[:1024] if isinstance(v, (list, dict)) else str(v)[:1024]
                entry[k] = val
            # Filter to only schema-valid fields
            if valid_fields is not None:
                entry = {k: v for k, v in entry.items() if k in valid_fields}
            insert_data.append(entry)
        client.insert(collection_name=collection_name, data=insert_data)
        client.flush(collection_name)
        logger.info(
            "Inserted %d records into Milvus collection '%s'.",
            len(insert_data),
            collection_name,
        )
        return len(insert_data)
    except Exception as exc:
        logger.warning(
            "Milvus insert failed for %s: %s",
            collection_name,
            exc,
        )
        return 0


# ===================================================================
# SEED FUNCTIONS
# ===================================================================


def seed_landmark_trials(dry_run: bool = False) -> int:
    """Seed the trial_protocols collection with landmark clinical trials.

    Args:
        dry_run: If True, count records without inserting.

    Returns:
        Number of records seeded.
    """
    logger.info("Seeding %d landmark trials ...", EXPECTED_LANDMARK_TRIALS)

    parser = ClinicalTrialsParser()
    raw_data = parser.seed_landmark_trials()
    records = parser.parse(raw_data)

    # Validate
    valid_records = [r for r in records if parser.validate_record(r)]
    logger.info(
        "Validated %d / %d landmark trial records.",
        len(valid_records), len(records),
    )

    if dry_run:
        logger.info("Dry run: would seed %d landmark trial records.", len(valid_records))
        return len(valid_records)

    return _insert_records("trial_protocols", valid_records, text_field="text")


def seed_regulatory_milestones(dry_run: bool = False) -> int:
    """Seed the trial_regulatory collection with FDA regulatory milestones.

    Args:
        dry_run: If True, count records without inserting.

    Returns:
        Number of records seeded.
    """
    logger.info("Seeding %d regulatory milestones ...", EXPECTED_REGULATORY_MILESTONES)

    parser = RegulatoryParser()
    records = parser.parse(REGULATORY_MILESTONES)

    # Validate
    valid_records = [r for r in records if parser.validate_record(r)]
    logger.info(
        "Validated %d / %d regulatory milestone records.",
        len(valid_records), len(records),
    )

    if dry_run:
        logger.info("Dry run: would seed %d regulatory milestone records.", len(valid_records))
        return len(valid_records)

    return _insert_records("trial_regulatory", valid_records, text_field="text")


def seed_all(dry_run: bool = False) -> Dict[str, int]:
    """Seed all knowledge collections.

    Args:
        dry_run: If True, count records without inserting.

    Returns:
        Dict mapping collection name to number of records seeded.
    """
    logger.info("Starting knowledge base seeding (dry_run=%s)...", dry_run)
    results = {}
    results["trial_protocols"] = seed_landmark_trials(dry_run=dry_run)
    results["trial_regulatory"] = seed_regulatory_milestones(dry_run=dry_run)
    logger.info("Seeding complete: %s", results)
    return results


# ===================================================================
# CLI ENTRY POINT
# ===================================================================


def main() -> None:
    """CLI entry point for seeding knowledge collections."""
    parser = argparse.ArgumentParser(
        description="Seed Clinical Trial Intelligence Agent knowledge collections"
    )
    parser.add_argument("--host", default="localhost", help="Milvus host")
    parser.add_argument("--port", type=int, default=19530, help="Milvus port")
    parser.add_argument("--dry-run", action="store_true", help="Log without inserting")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    seed_all(dry_run=args.dry_run)


if __name__ == "__main__":
    main()
