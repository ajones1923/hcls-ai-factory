"""Tests for hcls_common.embedder — memory cache and base class logic."""

import numpy as np
import pytest

from hcls_common.embedder import BaseEmbedder, _EmbeddingMemoryCache


class TestEmbeddingMemoryCache:
    def test_put_and_get(self):
        cache = _EmbeddingMemoryCache(max_size=10)
        vec = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        cache.put("hello world", vec)
        result = cache.get("hello world")
        assert result is not None
        np.testing.assert_array_equal(result, vec)

    def test_miss_returns_none(self):
        cache = _EmbeddingMemoryCache(max_size=10)
        assert cache.get("nonexistent") is None

    def test_lru_eviction(self):
        cache = _EmbeddingMemoryCache(max_size=2)
        cache.put("a", np.array([1.0]))
        cache.put("b", np.array([2.0]))
        cache.put("c", np.array([3.0]))  # evicts "a"
        assert cache.get("a") is None
        assert cache.get("c") is not None

    def test_lru_access_preserves(self):
        cache = _EmbeddingMemoryCache(max_size=2)
        cache.put("a", np.array([1.0]))
        cache.put("b", np.array([2.0]))
        cache.get("a")  # move "a" to end
        cache.put("c", np.array([3.0]))  # should evict "b"
        assert cache.get("a") is not None
        assert cache.get("b") is None

    def test_get_many(self):
        cache = _EmbeddingMemoryCache(max_size=10)
        cache.put("alpha", np.array([1.0]))
        cache.put("beta", np.array([2.0]))
        found = cache.get_many(["alpha", "gamma", "beta"])
        assert 0 in found
        assert 1 not in found
        assert 2 in found

    def test_put_many(self):
        cache = _EmbeddingMemoryCache(max_size=10)
        texts = ["x", "y", "z"]
        vecs = np.array([[1.0], [2.0], [3.0]])
        cache.put_many(texts, vecs)
        assert cache.get("x") is not None
        assert cache.get("y") is not None
        assert cache.get("z") is not None

    def test_clear(self):
        cache = _EmbeddingMemoryCache(max_size=10)
        cache.put("a", np.array([1.0]))
        cache.clear()
        assert len(cache) == 0
        assert cache.get("a") is None

    def test_key_is_deterministic(self):
        key1 = _EmbeddingMemoryCache._key("test string")
        key2 = _EmbeddingMemoryCache._key("test string")
        assert key1 == key2

    def test_key_differs_for_different_text(self):
        key1 = _EmbeddingMemoryCache._key("alpha")
        key2 = _EmbeddingMemoryCache._key("beta")
        assert key1 != key2


class _DummyEmbedder(BaseEmbedder):
    """Concrete embedder for testing base class logic."""

    @property
    def provider_name(self) -> str:
        return "dummy"

    def _encode_batch(self, texts):
        # Return random vectors of the correct dimension
        return np.random.randn(len(texts), self.dimension).astype(np.float32)


class TestBaseEmbedder:
    def test_embed_text_returns_vector(self):
        embedder = _DummyEmbedder(dimension=8)
        vec = embedder.embed_text("test sentence")
        assert vec.shape == (8,)

    def test_embed_texts_returns_matrix(self):
        embedder = _DummyEmbedder(dimension=8)
        vecs = embedder.embed_texts(["hello", "world", "test"])
        assert vecs.shape == (3, 8)

    def test_embed_texts_empty_list(self):
        embedder = _DummyEmbedder(dimension=8)
        vecs = embedder.embed_texts([])
        assert vecs.shape == (0, 8)

    def test_embed_text_uses_cache(self):
        embedder = _DummyEmbedder(dimension=8)
        vec1 = embedder.embed_text("same text")
        vec2 = embedder.embed_text("same text")
        # Should be the same cached vector
        np.testing.assert_array_equal(vec1, vec2)

    def test_embed_query_adds_bge_prefix(self):
        embedder = _DummyEmbedder(
            model_name="BAAI/bge-small-en-v1.5", dimension=8
        )
        # Just verify it doesn't crash and returns correct shape
        vec = embedder.embed_query("what is BRCA1?")
        assert vec.shape == (8,)

    def test_clear_caches(self):
        embedder = _DummyEmbedder(dimension=8)
        embedder.embed_text("cached text")
        embedder.clear_caches()
        assert len(embedder._mem_cache) == 0
