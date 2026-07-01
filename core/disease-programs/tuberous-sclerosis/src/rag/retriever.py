"""
RAG store + retriever (PRD §2.5.4). An in-memory cosine store runs here with no Milvus
server; the Milvus backend is the prod target (same interface). Retrieval returns chunks
with their source_uri so downstream attribution is exact.
"""
from __future__ import annotations

import numpy as np

from src.rag.corpus import SEED_CORPUS
from src.rag.embeddings import embed, mode


class InMemoryVectorStore:
    def __init__(self) -> None:
        self._items: list[dict] = []
        self._mat: np.ndarray | None = None

    def upsert(self, items: list[dict]) -> None:
        self._items.extend(items)
        vecs = embed([it["text"] for it in items])
        self._mat = vecs if self._mat is None else np.vstack([self._mat, vecs])

    def search(self, query: str, k: int = 4, partition: str | None = None) -> list[dict]:
        if self._mat is None:
            return []
        q = embed([query])[0]
        sims = self._mat @ q
        idx = np.argsort(-sims)
        out = []
        for i in idx:
            it = self._items[int(i)]
            if partition and it.get("partition") != partition:
                continue
            out.append({**it, "score": round(float(sims[int(i)]), 4)})
            if len(out) >= k:
                break
        return out


_retriever: "TSCRetriever | None" = None


def make_vector_store():
    """Milvus in prod (TSC_USE_MILVUS=1); in-memory cosine store otherwise (the default)."""
    import os

    if os.environ.get("TSC_USE_MILVUS") == "1":
        from src.rag.milvus_store import MilvusVectorStore

        return MilvusVectorStore()
    return InMemoryVectorStore()


class TSCRetriever:
    def __init__(self, store=None) -> None:
        self.store = store or make_vector_store()
        if isinstance(self.store, InMemoryVectorStore) and not self.store._items:
            self.store.upsert(SEED_CORPUS)

    def retrieve(self, query: str, k: int = 4, partition: str | None = None) -> list[dict]:
        return self.store.search(query, k=k, partition=partition)

    @property
    def embedding_mode(self) -> str:
        return mode()


def get_retriever() -> TSCRetriever:
    global _retriever
    if _retriever is None:
        _retriever = TSCRetriever()
    return _retriever
