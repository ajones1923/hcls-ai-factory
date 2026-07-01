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


class MilvusVectorStore:  # pragma: no cover - needs a Milvus server
    def __init__(self, uri: str | None = None) -> None:
        from pymilvus import DataType, MilvusClient  # noqa: WPS433

        self.client = MilvusClient(uri=uri or settings.MILVUS_URI)
        self._dim = int(embed(["dimension probe"]).shape[1])
        if not self.client.has_collection(COLLECTION):
            schema = self.client.create_schema(auto_id=True, enable_dynamic_field=True)
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
        } for i, it in enumerate(items)]
        self.client.insert(COLLECTION, rows)

    def search(self, query: str, k: int = 4, partition: str | None = None) -> list[dict]:
        q = embed([query])[0].tolist()
        flt = f'partition == "{partition}"' if partition else ""
        res = self.client.search(
            COLLECTION, data=[q], limit=k, filter=flt,
            output_fields=["text", "source_uri", "partition", "pub_year", "section"],
        )[0]
        return [{**h["entity"], "score": round(float(h["distance"]), 4)} for h in res]
