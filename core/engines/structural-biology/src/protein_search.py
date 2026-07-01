"""
Protein sequence search (B2) — a vector index of protein embeddings with BLAST-like
similarity search, backed by the factory's existing Milvus vector DB.

Two interchangeable backends:
  * ``MilvusBackend``   — real, uses the running Milvus (reuses the platform vector DB).
  * ``InMemoryBackend`` — exact cosine search; used in tests and as a zero-dep fallback.

Vectors come from ``ESM2Embedder`` (L2-normalized), so inner-product == cosine similarity.
"""
from __future__ import annotations

from typing import Any, Protocol, Sequence

from esm2_embeddings import ESM2Embedder, ESM2_DIM


class VectorBackend(Protocol):
    def ensure(self) -> None: ...
    def upsert(self, items: list[tuple[int, str, list[float]]]) -> None: ...
    def search(self, vec: Sequence[float], k: int) -> list[tuple[str, float]]: ...
    def count(self) -> int: ...


# --------------------------------------------------------------------------- #
class InMemoryBackend:
    """Exact cosine search over normalized vectors (dot product). Test/fallback backend."""

    def __init__(self, dim: int = ESM2_DIM) -> None:
        self.dim = dim
        self._items: list[tuple[int, str, list[float]]] = []

    def ensure(self) -> None:
        pass

    def upsert(self, items: list[tuple[int, str, list[float]]]) -> None:
        seen = {i for i, _, _ in self._items}
        for it in items:                       # dedup by id within the batch AND vs existing
            if it[0] not in seen:
                self._items.append(it)
                seen.add(it[0])

    def search(self, vec: Sequence[float], k: int) -> list[tuple[str, float]]:
        scored = [(name, sum(a * b for a, b in zip(vec, v))) for _, name, v in self._items]
        scored.sort(key=lambda t: -t[1])
        return scored[:k]

    def count(self) -> int:
        return len(self._items)


# --------------------------------------------------------------------------- #
class MilvusBackend:
    """Real Milvus-backed protein index (reuses the platform vector DB)."""

    def __init__(self, collection: str = "protein_sequences", dim: int = ESM2_DIM,
                 host: str = "127.0.0.1", port: str = "19530") -> None:
        self.collection_name = collection
        self.dim = dim
        self.host, self.port = host, port
        self._coll = None

    def ensure(self) -> None:
        from pymilvus import (Collection, CollectionSchema, DataType, FieldSchema,
                              connections, utility)
        connections.connect(host=self.host, port=self.port, timeout=10)
        if not utility.has_collection(self.collection_name):
            schema = CollectionSchema([
                FieldSchema("id", DataType.INT64, is_primary=True),
                FieldSchema("name", DataType.VARCHAR, max_length=256),
                FieldSchema("embedding", DataType.FLOAT_VECTOR, dim=self.dim),
            ], description="ESM-2 protein sequence embeddings")
            coll = Collection(self.collection_name, schema)
            coll.create_index("embedding", {"index_type": "FLAT", "metric_type": "IP"})
        self._coll = Collection(self.collection_name)
        self._coll.load()

    def upsert(self, items: list[tuple[int, str, list[float]]]) -> None:
        if self._coll is None:
            self.ensure()
        self._coll.insert([[i for i, _, _ in items],
                           [n for _, n, _ in items],
                           [v for _, _, v in items]])
        self._coll.flush()

    def search(self, vec: Sequence[float], k: int) -> list[tuple[str, float]]:
        if self._coll is None:
            self.ensure()
        res = self._coll.search([list(vec)], "embedding", {"metric_type": "IP"},
                                limit=k, output_fields=["name"])
        return [(h.entity.get("name"), float(h.score)) for h in res[0]]

    def count(self) -> int:
        return self._coll.num_entities if self._coll is not None else 0


# --------------------------------------------------------------------------- #
class ProteinSearchIndex:
    """Embed proteins and search them by sequence similarity."""

    def __init__(self, embedder: ESM2Embedder | None = None, backend: VectorBackend | None = None) -> None:
        self.embedder = embedder or ESM2Embedder()
        self.backend = backend or MilvusBackend()
        self._seqs: dict[str, str] = {}          # name->sequence, for the B3 SW re-rank

    def add(self, proteins: list[dict[str, Any]]) -> int:
        """proteins: [{"id": int, "name": str, "sequence": str}, ...]"""
        self.backend.ensure()
        items = [(p["id"], p["name"], self.embedder.embed(p["sequence"])) for p in proteins]
        self.backend.upsert(items)
        for p in proteins:
            self._seqs[p["name"]] = p["sequence"]
        return len(items)

    def search(self, sequence: str, top_k: int = 5, rerank: bool = False) -> list[dict[str, Any]]:
        """ANN cosine search; with rerank=True, recall a wider set then exact-SW re-rank (B3)."""
        recall = top_k * 4 if rerank else top_k
        qv = self.embedder.embed(sequence)
        hits = self.backend.search(qv, recall)
        results = [{"name": name, "cosine": round(score, 4)} for name, score in hits]
        if rerank:
            from protein_rerank import rerank as _rr
            results = _rr(sequence, results, self._seqs)[:top_k]
        return results

    def count(self) -> int:
        return self.backend.count()
