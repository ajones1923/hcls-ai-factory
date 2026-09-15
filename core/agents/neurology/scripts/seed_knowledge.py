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
        # These collections are structured (neuro_electrophysiology declares test_type,
        # finding, pattern, ... and no generic text field), and dynamic fields are off. An
        # unexpected key therefore aborts the whole insert -- which is how a record carrying
        # `text` killed the batch and the seeder reported "0 total records inserted".
        # Project each row onto the collection's declared fields, and drop the ones that
        # carry no usable content rather than inserting a bare vector.
        try:
            _fields = client.describe_collection(collection_name).get("fields", [])
            valid = {fld["name"] for fld in _fields}
            # Dynamic fields are off and nothing is nullable, so every declared column must be
            # present on every row or the whole batch is rejected. The parsers legitimately do
            # not carry all of them (a PubMed record has no `sequence`), so supply a
            # type-appropriate empty rather than dropping real records.
            _required, _numeric, _integral = {}, set(), set()
            for fld in _fields:
                nm = fld["name"]
                if nm in ("id", "embedding") or fld.get("auto_id"):
                    continue
                # pymilvus returns a DataType enum whose str() is its NUMERIC CODE
                # ("10"), not its name -- so matching on str() silently classified every
                # float column as text. Use .name.
                _t = fld.get("type", "")
                t = getattr(_t, "name", str(_t)).upper()
                is_num = ("INT" in t or "FLOAT" in t or "DOUBLE" in t) and "VECTOR" not in t
                if is_num:
                    _numeric.add(nm)
                    if "INT" in t:
                        _integral.add(nm)
                _required[nm] = 0 if is_num else ""
        except Exception:
            valid, _required, _numeric, _integral = None, {}, set(), set()

        def _flatten(rec):
            """Row fields live INSIDE IngestRecord.metadata, not on the record itself.

            The parsers emit IngestRecord(text=..., metadata={pmid, title, ...}); the schema
            declares pmid/title/... as top-level columns. Reading the dataclass's own
            attributes yielded {text, metadata, collection_name, record_id, source}, which
            shares no column with the schema -- hence "0 records inserted" from a seeder that
            had just embedded 49 real publications.
            """
            if isinstance(rec, dict):
                base = dict(rec)
                meta = base.pop("metadata", None)
            else:
                base = {"text": getattr(rec, "text", ""),
                        "source": getattr(rec, "source", "")}
                meta = getattr(rec, "metadata", None)
            if isinstance(meta, str):
                try:
                    import ast as _ast
                    meta = _ast.literal_eval(meta)
                except Exception:
                    meta = None
            if isinstance(meta, dict):
                for k, v in meta.items():
                    base.setdefault(k, v)
            return base

        data_rows = []
        for i, rec in enumerate(records):
            row = _flatten(rec)
            if valid:
                dropped = {k for k in row if k not in valid}
                row = {k: v for k, v in row.items() if k in valid}
                if dropped and i == 0:
                    logger.info(
                        "  '%s': ignoring %d field(s) absent from the schema: %s",
                        collection_name, len(dropped), ", ".join(sorted(dropped)[:6]),
                    )
            if not row:
                continue
            # Coerce to the column's declared type. A blanket str() here silently broke
            # numeric columns ("{immune_score} field should be a float, but got a str").
            def _coerce(k, v):
                if k == "embedding":
                    return v
                if k in _numeric:
                    # int64 and float are different columns; coercing everything to float
                    # fails an int64 insert ("'float' object cannot be interpreted as an
                    # integer").
                    try:
                        return int(float(v)) if k in _integral else float(v)
                    except (TypeError, ValueError):
                        return 0 if k in _integral else 0.0
                return v if isinstance(v, (int, float, bool)) else str(v)[:4096]

            row = {k: _coerce(k, v) for k, v in row.items()}
            # Preserve the source text. The parsers and the schemas were designed
            # independently, so for some collections only one or two columns overlap and the
            # rest would be filled with empty defaults -- rows that embed correctly (the
            # vector is built from the real text) but return nothing readable. If the schema
            # offers a prose column and nothing has filled it, put the record's text there.
            _src_text = texts[i] if i < len(texts) else ""
            if _src_text and valid:
                for _cand in ("description", "text", "summary", "abstract",
                              "content", "finding", "clinical_correlation"):
                    if _cand in valid and not row.get(_cand):
                        row[_cand] = _src_text[:4096]
                        break
            for _k, _default in _required.items():
                row.setdefault(_k, _default)
            row["embedding"] = embeddings[i]
            data_rows.append(row)

        if not data_rows:
            logger.warning(
                "No rows for '%s' survived schema projection -- the seed records and the "
                "collection schema share no fields.", collection_name,
            )
            return 0

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
