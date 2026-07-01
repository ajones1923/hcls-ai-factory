"""Tests for hcls_common.milvus_client — validation helpers and TTLLRUCache."""

import time

import pytest

from hcls_common.milvus_client import (
    TTLLRUCache,
    UnifiedMilvusClient,
    sanitize_gene_name,
    sanitize_chromosome,
    validate_milvus_filter,
)


class TestSanitizeGeneName:
    def test_valid_gene(self):
        assert sanitize_gene_name("BRCA1") == "BRCA1"

    def test_strips_whitespace(self):
        assert sanitize_gene_name("  TP53  ") == "TP53"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            sanitize_gene_name("")

    def test_injection_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            sanitize_gene_name('BRCA1" OR 1=1 --')


class TestSanitizeChromosome:
    def test_valid(self):
        assert sanitize_chromosome("chr1") == "chr1"
        assert sanitize_chromosome("X") == "X"
        assert sanitize_chromosome("MT") == "MT"

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            sanitize_chromosome("")

    def test_invalid_raises(self):
        with pytest.raises(ValueError, match="Invalid"):
            sanitize_chromosome("chromosome_99")


class TestValidateMilvusFilter:
    def test_valid_expression(self):
        expr = 'gene == "EGFR" and impact == "HIGH"'
        assert validate_milvus_filter(expr) == expr

    def test_numeric_expression(self):
        expr = "pos >= 100 and pos <= 200"
        assert validate_milvus_filter(expr) == expr

    def test_in_expression(self):
        expr = 'gene in ["BRCA1", "BRCA2"]'
        assert validate_milvus_filter(expr) == expr

    def test_injection_raises(self):
        with pytest.raises(ValueError, match="disallowed"):
            validate_milvus_filter("gene == 'x'; system('rm -rf /')")

    def test_empty_raises(self):
        with pytest.raises(ValueError, match="empty"):
            validate_milvus_filter("")


class TestTTLLRUCache:
    def test_put_and_get(self):
        cache = TTLLRUCache(max_size=10, ttl=60.0)
        cache.put("k1", "value1")
        hit, val = cache.get("k1")
        assert hit is True
        assert val == "value1"

    def test_miss_returns_false(self):
        cache = TTLLRUCache(max_size=10, ttl=60.0)
        hit, val = cache.get("nonexistent")
        assert hit is False
        assert val is None

    def test_ttl_expiration(self):
        cache = TTLLRUCache(max_size=10, ttl=0.05)
        cache.put("k1", "value1")
        time.sleep(0.1)
        hit, val = cache.get("k1")
        assert hit is False

    def test_lru_eviction(self):
        cache = TTLLRUCache(max_size=2, ttl=60.0)
        cache.put("k1", "v1")
        cache.put("k2", "v2")
        cache.put("k3", "v3")  # should evict k1
        hit, _ = cache.get("k1")
        assert hit is False
        hit, val = cache.get("k3")
        assert hit is True
        assert val == "v3"

    def test_lru_access_order(self):
        cache = TTLLRUCache(max_size=2, ttl=60.0)
        cache.put("k1", "v1")
        cache.put("k2", "v2")
        # Access k1 to make it recently used
        cache.get("k1")
        cache.put("k3", "v3")  # should evict k2 (least recently used)
        hit_k1, _ = cache.get("k1")
        hit_k2, _ = cache.get("k2")
        assert hit_k1 is True
        assert hit_k2 is False

    def test_invalidate(self):
        cache = TTLLRUCache(max_size=10, ttl=60.0)
        cache.put("k1", "v1")
        cache.invalidate("k1")
        hit, _ = cache.get("k1")
        assert hit is False

    def test_clear(self):
        cache = TTLLRUCache(max_size=10, ttl=60.0)
        cache.put("k1", "v1")
        cache.put("k2", "v2")
        cache.clear()
        assert len(cache) == 0

    def test_make_key_deterministic(self):
        cache = TTLLRUCache()
        key1 = cache._make_key("foo", "bar", 42)
        key2 = cache._make_key("foo", "bar", 42)
        assert key1 == key2

    def test_make_key_different_args(self):
        cache = TTLLRUCache()
        key1 = cache._make_key("foo", "bar")
        key2 = cache._make_key("foo", "baz")
        assert key1 != key2


class TestUnifiedMilvusClientInit:
    def test_default_instantiation(self):
        client = UnifiedMilvusClient()
        assert client.host == "localhost"
        assert client.port == 19530
        assert client.collection_name == "genomic_evidence"
        assert client.embedding_dim == 384

    def test_custom_params(self):
        client = UnifiedMilvusClient(
            host="milvus-prod",
            port=19531,
            collection_name="test_collection",
            embedding_dim=768,
        )
        assert client.host == "milvus-prod"
        assert client.port == 19531
        assert client.collection_name == "test_collection"
        assert client.embedding_dim == 768

    def test_cache_initialized(self):
        client = UnifiedMilvusClient(cache_max_size=512, cache_ttl=120.0)
        assert client._cache is not None
        assert client._cache._max_size == 512
        assert client._cache._ttl == 120.0

    def test_default_output_fields(self):
        assert "gene" in UnifiedMilvusClient.DEFAULT_OUTPUT_FIELDS
        assert "chrom" in UnifiedMilvusClient.DEFAULT_OUTPUT_FIELDS
        assert "clinical_significance" in UnifiedMilvusClient.DEFAULT_OUTPUT_FIELDS
