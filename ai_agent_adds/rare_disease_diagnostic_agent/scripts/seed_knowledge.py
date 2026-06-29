"""Seed Milvus collections with curated knowledge for the Rare Disease Diagnostic Agent.

Populates rd_diseases, rd_phenotypes, and rd_therapies collections with
seed data from the OMIM, HPO, Orphanet, and Gene Therapy parsers.

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

from src.ingest.omim_parser import OMIMParser, OMIM_DISEASES
from src.ingest.hpo_parser import HPOParser, HPO_TERMS
from src.ingest.orphanet_parser import OrphanetParser, ORPHANET_DISEASES
from src.ingest.gene_therapy_parser import GeneTherapyParser, ALL_GENE_THERAPIES
from src.ingest.base import IngestRecord

logger = logging.getLogger(__name__)


# ===================================================================
# SEED DATA COUNTS
# ===================================================================

EXPECTED_OMIM_DISEASES = len(OMIM_DISEASES)
EXPECTED_HPO_TERMS = len(HPO_TERMS)
EXPECTED_ORPHANET_DISEASES = len(ORPHANET_DISEASES)
EXPECTED_GENE_THERAPIES = len(ALL_GENE_THERAPIES)


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
        Number of records inserted (or that would have been inserted).
    """
    if not records:
        logger.info("No records to insert into '%s'.", collection_name)
        return 0

    try:
        from pymilvus import MilvusClient  # noqa: F811
    except ImportError:
        logger.warning(
            "pymilvus not installed -- dry run: %d records for '%s'",
            len(records), collection_name,
        )
        return len(records)

    try:
        from sentence_transformers import SentenceTransformer
    except ImportError:
        logger.warning(
            "sentence_transformers not installed -- dry run: %d records for '%s'",
            len(records), collection_name,
        )
        return len(records)

    try:
        model = SentenceTransformer("BAAI/bge-small-en-v1.5")
        texts = [getattr(r, text_field, r.text) for r in records]
        embeddings = model.encode(texts, show_progress_bar=False).tolist()

        client = MilvusClient(uri="http://localhost:19530")

        # Discover valid field names from the collection schema
        try:
            desc = client.describe_collection(collection_name)
            valid_fields = {f["name"] for f in desc.get("fields", [])}
        except Exception:
            valid_fields = None

        data = []
        for record, embedding in zip(records, embeddings):
            row = {
                "embedding": embedding,
                "text": record.text,
                "source": record.source,
            }
            # Only include metadata keys that match schema fields
            for k, v in record.metadata.items():
                if valid_fields is None or k in valid_fields:
                    row[k] = v
            # Remove any keys not in the schema
            if valid_fields is not None:
                row = {k: v for k, v in row.items() if k in valid_fields}
            data.append(row)

        client.insert(collection_name=collection_name, data=data)
        client.flush(collection_name)
        logger.info("Inserted %d records into '%s'.", len(data), collection_name)
        return len(data)

    except Exception as exc:
        logger.warning(
            "Milvus insert failed for '%s': %s -- dry run: %d records",
            collection_name, exc, len(records),
        )
        return 0


# ===================================================================
# SEED FUNCTIONS
# ===================================================================


def seed_omim() -> int:
    """Seed OMIM diseases into rd_diseases."""
    parser = OMIMParser()
    records, stats = parser.run()
    return _insert_records("rd_diseases", records)


def seed_hpo() -> int:
    """Seed HPO terms into rd_phenotypes."""
    parser = HPOParser()
    records, stats = parser.run()
    return _insert_records("rd_phenotypes", records)


def seed_orphanet() -> int:
    """Seed Orphanet diseases into rd_diseases."""
    parser = OrphanetParser()
    records, stats = parser.run()
    return _insert_records("rd_diseases", records)


def seed_gene_therapies() -> int:
    """Seed gene therapies into rd_therapies."""
    parser = GeneTherapyParser()
    records, stats = parser.run()
    return _insert_records("rd_therapies", records)


def seed_all() -> Dict[str, int]:
    """Seed all knowledge bases."""
    results = {
        "omim": seed_omim(),
        "hpo": seed_hpo(),
        "orphanet": seed_orphanet(),
        "gene_therapies": seed_gene_therapies(),
    }
    total = sum(results.values())
    logger.info("Seeded %d total records: %s", total, results)
    return results


# ===================================================================
# CLI
# ===================================================================


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Seed Rare Disease Diagnostic Agent knowledge base"
    )
    parser.add_argument("--host", default="localhost", help="Milvus host")
    parser.add_argument("--port", type=int, default=19530, help="Milvus port")
    parser.add_argument(
        "--source",
        choices=["omim", "hpo", "orphanet", "gene_therapy", "all"],
        default="all",
        help="Source to seed",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.source == "all":
        seed_all()
    elif args.source == "omim":
        seed_omim()
    elif args.source == "hpo":
        seed_hpo()
    elif args.source == "orphanet":
        seed_orphanet()
    elif args.source == "gene_therapy":
        seed_gene_therapies()


if __name__ == "__main__":
    main()
