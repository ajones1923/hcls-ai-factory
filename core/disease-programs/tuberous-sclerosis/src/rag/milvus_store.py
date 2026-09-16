"""
Milvus vector-store backend (PRD §2.5.4). The production RAG store: the `tsc_literature`
collection, partitioned by source. Same interface as the in-memory store, so the retriever
is backend-agnostic. Requires a running Milvus (docker-compose service) and `pymilvus`.
Selected via TSC_USE_MILVUS=1; otherwise the in-memory cosine store is used.
"""
from __future__ import annotations

import numpy as np

from config.settings import settings
from src.rag.embeddings import embed

COLLECTION = "tsc_literature"


def _content_pk(item: dict) -> int:
    """A stable int64 key from the chunk's identity, so re-ingest is idempotent."""
    import hashlib
    ident = f"{item.get('id') or ''}|{item.get('source_uri', '')}|{item.get('text', '')}"
    return int.from_bytes(hashlib.blake2b(ident.encode(), digest_size=8).digest(), "big") >> 1


class MilvusVectorStore:  # pragma: no cover - needs a Milvus server
    def __init__(self, uri: str | None = None) -> None:
        from pymilvus import DataType, MilvusClient  # noqa: WPS433

        self.client = MilvusClient(uri=uri or settings.MILVUS_URI)
        self._dim = int(embed(["dimension probe"]).shape[1])
        self._legacy_auto_id = False
        if self.client.has_collection(COLLECTION):
            # A collection created before 2026-09-16 used auto_id, so every re-run of
            # scripts/load_rag.py APPENDED the whole seed corpus again -- 6 chunks became 12,
            # and duplicate passages skew retrieval while looking like a healthy corpus.
            # Not dropped automatically: by then it may hold real ingested literature.
            try:
                desc = self.client.describe_collection(COLLECTION)
                self._legacy_auto_id = bool(desc.get("auto_id"))
            except Exception:
                pass
        else:
            schema = self.client.create_schema(auto_id=False, enable_dynamic_field=True)
            schema.add_field("pk", DataType.INT64, is_primary=True)
            schema.add_field("vector", DataType.FLOAT_VECTOR, dim=self._dim)
            schema.add_field("text", DataType.VARCHAR, max_length=8192)
            schema.add_field("source_uri", DataType.VARCHAR, max_length=1024)
            schema.add_field("partition", DataType.VARCHAR, max_length=64)
            index = self.client.prepare_index_params()
            index.add_index("vector", metric_type="COSINE", index_type="IVF_FLAT", params={"nlist": 128})
            self.client.create_collection(COLLECTION, schema=schema, index_params=index)

    def upsert(self, items: list[dict]) -> None:
        vecs = embed([it["text"] for it in items])
        rows = [{
            "vector": vecs[i].tolist(), "text": it["text"],
            "source_uri": it.get("source_uri", ""), "partition": it.get("partition", ""),
            "pub_year": it.get("pub_year"), "section": it.get("section"),
            **({} if self._legacy_auto_id else {"pk": _content_pk(it)}),
        } for i, it in enumerate(items)]
        if self._legacy_auto_id:
            # Cannot dedupe without a deterministic key; preserve old behaviour and say so.
            import logging
            logging.getLogger(__name__).warning(
                "%s was created with auto_id: re-ingest APPENDS instead of replacing. Drop it "
                "and re-run scripts/load_rag.py to get idempotent upserts.", COLLECTION)
            self.client.insert(COLLECTION, rows)
        else:
            # Deterministic pk from the content, so re-running the loader replaces rather than
            # duplicates -- the same rule lib/hcls_common/ingest_persist.py follows.
            self.client.upsert(COLLECTION, rows)
        # Flush, or `num_entities` keeps reporting 0 while the rows sit in a growing segment.
        # The loader then prints "Ingested 6 chunks" and every corpus check, dashboard and
        # operator sees an EMPTY collection -- success reported, nothing visible. Found
        # 2026-09-16: tsc_literature read 0 entities until an explicit flush turned it into 6.
        try:
            self.client.flush(COLLECTION)
        except Exception:                      # older clients expose it only on the ORM handle
            try:
                from pymilvus import Collection
                Collection(COLLECTION).flush()
            except Exception:
                pass

    def search(self, query: str, k: int = 4, partition: str | None = None) -> list[dict]:
        q = embed([query])[0].tolist()
        flt = f'partition == "{partition}"' if partition else ""
        res = self.client.search(
            COLLECTION, data=[q], limit=k, filter=flt,
            output_fields=["text", "source_uri", "partition", "pub_year", "section"],
        )[0]
        return [{**h["entity"], "score": round(float(h["distance"]), 4)} for h in res]
