"""Write IngestRecords to Milvus — one implementation, shared by every agent's ingest.

Several agents' `run_ingest.py` fetched, parsed, validated and then dropped everything on the
floor: they logged "60 records validated" and wrote nothing, so a working ClinicalTrials.gov or
PubMed ingest left the corpus exactly as it found it. That is why several agents sat at a few
dozen rows while their ingest scripts "worked".

This lives here rather than being copied per agent because the projection rules below were each
learned from a failed insert, and a copy that misses one fails silently — Milvus aborts the
whole batch on a single bad key, and callers log a warning and move on.

Projection rules, each from a real failure:

  * Dynamic fields are OFF and nothing is nullable, so every declared column must be present
    and no undeclared key may appear. Project onto the schema, then fill the remainder.
  * `str(DataType.FLOAT)` is its numeric CODE ("10"), not "FLOAT" — matching on str() silently
    classifies every float column as text. Use `.name`.
  * int64 and float are different: coercing both to float fails an int64 insert.
  * Only genuine ARRAY columns may receive a list. A list in a VARCHAR column aborts the batch.
  * BOOL is neither numeric nor text — defaulting it to "" fails the insert.
  * The primary key differs per collection (VARCHAR in some, int64 in others) and is not
    auto_id, so it must be supplied in the right type.
  * A content-derived id makes re-ingest idempotent via upsert instead of duplicating the
    corpus on every run.
"""
from __future__ import annotations

import hashlib
import logging
from collections import defaultdict
from typing import Any, Iterable

logger = logging.getLogger(__name__)

DEFAULT_EMBED_MODEL = "BAAI/bge-small-en-v1.5"
_TEXT_FIELDS = ("text", "text_chunk", "text_content", "text_summary", "description",
                "abstract", "content", "summary")


def _type_name(field: dict) -> str:
    t = field.get("type", "")
    return getattr(t, "name", str(t)).upper()


def _record_text(rec: Any) -> str:
    for attr in ("text", "content", "summary"):
        v = getattr(rec, attr, None) if not isinstance(rec, dict) else rec.get(attr)
        if isinstance(v, str) and v.strip():
            return v
    return str(rec)


def _record_meta(rec: Any) -> dict:
    m = getattr(rec, "metadata", None) if not isinstance(rec, dict) else rec.get("metadata")
    if isinstance(m, str):
        try:
            import ast
            m = ast.literal_eval(m)
        except Exception:
            m = None
    return m if isinstance(m, dict) else {}


def persist_records(records: Iterable[Any], default_collection: str, *,
                    milvus_uri: str = "http://localhost:19530",
                    embed_model: str = DEFAULT_EMBED_MODEL) -> int:
    """Embed and upsert records into Milvus. Returns rows written; never raises."""
    records = list(records or [])
    if not records:
        return 0
    try:
        from pymilvus import MilvusClient
        from sentence_transformers import SentenceTransformer
    except Exception as exc:                       # optional deps — degrade, don't crash
        logger.warning("Cannot persist (%s); install pymilvus + sentence-transformers", exc)
        return 0

    client = MilvusClient(uri=milvus_uri)
    model = SentenceTransformer(embed_model)

    groups: dict[str, list] = defaultdict(list)
    for r in records:
        coll = (getattr(r, "collection_name", None)
                or (r.get("collection_name") if isinstance(r, dict) else None)
                or default_collection)
        groups[coll].append(r)

    written = 0
    for coll, rows in groups.items():
        try:
            fields = client.describe_collection(coll).get("fields", [])
        except Exception:
            logger.warning("  %s: no such collection — run setup_collections.py first", coll)
            continue

        valid = {f["name"] for f in fields}
        numeric = {f["name"] for f in fields
                   if ("INT" in _type_name(f) or "FLOAT" in _type_name(f) or
                       "DOUBLE" in _type_name(f)) and "VECTOR" not in _type_name(f)}
        arrays = {f["name"] for f in fields if "ARRAY" in _type_name(f)}
        # BOOL is neither numeric nor text: defaulting it to "" fails the insert.
        bools = {f["name"] for f in fields if _type_name(f) == "BOOL"}
        text_field = next((c for c in _TEXT_FIELDS if c in valid), None)

        texts = [_record_text(r) for r in rows]
        vectors = model.encode(texts, show_progress_bar=False).tolist()

        payload, dropped_once = [], False
        for rec, text, vec in zip(rows, texts, vectors):
            row = {k: v for k, v in _record_meta(rec).items() if k in valid}
            if not dropped_once:
                extra = set(_record_meta(rec)) - valid
                if extra:
                    logger.info("  %s: ignoring %d field(s) absent from the schema: %s",
                                coll, len(extra), ", ".join(sorted(extra)[:6]))
                    dropped_once = True
            if text_field and not row.get(text_field):
                row[text_field] = text[:8192]

            for name in valid:                      # every declared column must be present
                if name in ("id", "embedding"):
                    continue
                row.setdefault(name, False if name in bools
                               else (0 if name in numeric else ""))

            for name, value in list(row.items()):   # coerce to the declared type
                if name in ("id", "embedding"):
                    continue
                if name in bools:
                    row[name] = (value if isinstance(value, bool)
                                 else str(value).strip().lower() in ("true", "1", "yes", "y"))
                elif name in numeric:
                    try:
                        row[name] = (int(float(value)) if "INT" in _type_name(
                            next(f for f in fields if f["name"] == name)) else float(value))
                    except (TypeError, ValueError):
                        row[name] = 0
                elif isinstance(value, list) and name not in arrays:
                    row[name] = ", ".join(str(x) for x in value)[:4096]
                elif not isinstance(value, (str, list)):
                    row[name] = str(value)[:4096]
                elif isinstance(value, str):
                    row[name] = value[:4096]

            if "id" in valid:
                digest = hashlib.sha1(f"{coll}:{text}".encode()).hexdigest()
                row["id"] = int(digest[:15], 16) if "id" in numeric else f"{coll}-{digest[:24]}"
            row["embedding"] = vec
            payload.append(row)

        try:
            client.upsert(collection_name=coll, data=payload)
            client.flush(coll)
            client.load_collection(collection_name=coll)   # else search says "not loaded"
            logger.info("  persisted %d rows into %s", len(payload), coll)
            written += len(payload)
        except Exception as exc:
            logger.warning("  %s: insert failed: %s", coll, str(exc)[:200])
    return written
