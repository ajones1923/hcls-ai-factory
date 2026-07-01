"""Tests for ESM-2 embeddings + protein search (B2). Model + Milvus are both injected out."""
import math

import pytest

from esm2_embeddings import ESM2Embedder, l2_normalize
from esmfold_service import SequenceError
from protein_search import InMemoryBackend, ProteinSearchIndex


# a tiny deterministic "model": maps sequences to known 4-d vectors
_TABLE = {
    "AAAA": [1.0, 0.0, 0.0, 0.0],
    "AAAC": [0.9, 0.1, 0.0, 0.0],   # close to AAAA
    "CCCC": [0.0, 0.0, 1.0, 0.0],   # orthogonal
}
def fake_embed(seq):
    return _TABLE.get(seq, [0.5, 0.5, 0.5, 0.5])

def embedder():
    return ESM2Embedder(dim=4, _embed_fn=fake_embed)


class TestEmbedder:
    def test_l2_normalize(self):
        v = l2_normalize([3.0, 4.0])
        assert math.isclose(v[0], 0.6) and math.isclose(v[1], 0.8)

    def test_embed_normalizes_and_validates(self):
        v = embedder().embed("CCCC")
        assert math.isclose(math.sqrt(sum(x * x for x in v)), 1.0, abs_tol=1e-6)

    def test_embed_rejects_bad_sequence(self):
        with pytest.raises(SequenceError):
            embedder().embed("XZ123")

    def test_embed_dim_mismatch_raises(self):
        bad = ESM2Embedder(dim=8, _embed_fn=fake_embed)   # model gives 4, expects 8
        with pytest.raises(ValueError, match="dim"):
            bad.embed("AAAA")

    def test_batch(self):
        assert len(embedder().embed_batch(["AAAA", "CCCC"])) == 2


class TestInMemoryBackend:
    def test_upsert_dedup_and_count(self):
        b = InMemoryBackend(dim=4)
        b.upsert([(1, "a", [1, 0, 0, 0]), (1, "a", [1, 0, 0, 0])])  # dup id
        assert b.count() == 1

    def test_cosine_ranking(self):
        b = InMemoryBackend(dim=4)
        b.upsert([(1, "x", [1, 0, 0, 0]), (2, "y", [0, 1, 0, 0])])
        hits = b.search([1, 0, 0, 0], k=2)
        assert hits[0][0] == "x" and hits[0][1] > hits[1][1]


class TestProteinSearchIndex:
    def test_add_and_self_retrieval(self):
        idx = ProteinSearchIndex(embedder(), InMemoryBackend(dim=4))
        n = idx.add([
            {"id": 1, "name": "ubq", "sequence": "AAAA"},
            {"id": 2, "name": "ubq_like", "sequence": "AAAC"},
            {"id": 3, "name": "other", "sequence": "CCCC"},
        ])
        assert n == 3 and idx.count() == 3
        hits = idx.search("AAAA", top_k=3)
        # exact match ranks first (cosine ~1.0), the near-variant second, orthogonal last
        assert hits[0]["name"] == "ubq" and hits[0]["cosine"] == pytest.approx(1.0, abs=1e-3)
        assert hits[1]["name"] == "ubq_like" and hits[1]["cosine"] > 0.9
        assert hits[2]["name"] == "other" and hits[2]["cosine"] == pytest.approx(0.0, abs=1e-6)
