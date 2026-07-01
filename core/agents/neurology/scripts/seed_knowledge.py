#!/usr/bin/env python3
"""Seed the neurology knowledge base with curated domain data.

Runs all three ingest parsers (PubMed, Neuroimaging, EEG) in seed mode
to populate the knowledge base with landmark papers, imaging protocols,
and EEG patterns.

Usage:
    python scripts/seed_knowledge.py
"""

import logging
import sys
from pathlib import Path
from typing import Any, List

# Ensure project root on sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from config.settings import settings

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


# ===================================================================
# INSERT HELPER
# ===================================================================


def _insert_records(
    collection_name: str,
    records: List[Any],
    text_field: str = "text",
) -> int:
    """Generate embeddings and insert records into a Milvus collection.

    Degrades gracefully: if pymilvus or sentence_transformers are not
    installed, or if Milvus is unreachable, logs a warning and returns
    the record count (as if it were a dry run).

    Parameters
    ----------
    collection_name : str
        Target Milvus collection name.
    records : list
        Records to insert.  Each must have a ``.text`` attribute or be a dict
        with a key matching *text_field*.
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
            "pymilvus is not installed -- skipping Milvus insert for %s "
            "(%d records would have been inserted).",
            collection_name,
            len(records),
        )
        return len(records)

    try:
        from sentence_transformers import SentenceTransformer
    except ImportError:
        logger.warning(
            "sentence-transformers is not installed -- skipping embedding "
            "generation for %s (%d records).",
            collection_name,
            len(records),
        )
        return len(records)

    # --- generate embeddings ---
    model = SentenceTransformer(settings.EMBEDDING_MODEL)
    texts = []
    for r in records:
        if isinstance(r, dict):
            texts.append(str(r.get(text_field, "")))
        elif hasattr(r, text_field):
            texts.append(str(getattr(r, text_field, "")))
        else:
            texts.append(str(r))

    embeddings = model.encode(texts, show_progress_bar=False).tolist()
    logger.info(
        "Generated %d embeddings for collection '%s'.",
        len(embeddings),
        collection_name,
    )

    # --- insert into Milvus ---
    try:
        client = MilvusClient(
            uri=f"http://{settings.MILVUS_HOST}:{settings.MILVUS_PORT}"
        )
        data_rows = []
        for i, rec in enumerate(records):
            row = dict(rec) if isinstance(rec, dict) else {"text": texts[i]}
            row["embedding"] = embeddings[i]
            data_rows.append(row)

        client.insert(collection_name=collection_name, data=data_rows)
        client.flush(collection_name)
        logger.info(
            "Inserted %d records into '%s'.", len(data_rows), collection_name
        )
        return len(data_rows)
    except Exception as exc:
        logger.warning(
            "Milvus insert failed for %s: %s",
            collection_name,
            exc,
        )
        return 0


def main():
    """Run all seed parsers and report results."""
    from src.ingest.pubmed_neuro_parser import PubMedNeuroParser
    from src.ingest.neuroimaging_parser import NeuroimagingParser
    from src.ingest.eeg_parser import EEGParser

    parsers = [
        ("PubMed Neurology", PubMedNeuroParser(), "neuro_literature"),
        ("Neuroimaging Protocols", NeuroimagingParser(), "neuro_imaging"),
        ("EEG Patterns", EEGParser(), "neuro_electrophysiology"),
    ]

    total_records = 0

    for name, parser, collection in parsers:
        logger.info("Running %s seed ingest ...", name)
        records, stats = parser.run()
        logger.info(
            "  %s: %d fetched, %d parsed, %d validated, %d errors (%.1fs)",
            name,
            stats.total_fetched,
            stats.total_parsed,
            stats.total_validated,
            stats.total_errors,
            stats.duration_seconds,
        )
        inserted = _insert_records(collection, records)
        total_records += inserted

    logger.info("Seed complete: %d total records inserted", total_records)


if __name__ == "__main__":
    main()
