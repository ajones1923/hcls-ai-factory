"""Tests for hcls_common.llm_client — JSON extraction, token budget, response cache, factory."""

import pytest

from hcls_common.llm_client import (
    BaseLLMClient,
    LLMClientFactory,
    _ResponseCache,
    _TokenBudget,
    extract_json,
)


class TestExtractJson:
    def test_bare_json_object(self):
        result = extract_json('{"key": "value"}')
        assert result == {"key": "value"}

    def test_bare_json_array(self):
        result = extract_json('[1, 2, 3]')
        assert result == [1, 2, 3]

    def test_markdown_fenced_json(self):
        text = '```json\n{"gene": "BRCA1", "impact": "HIGH"}\n```'
        result = extract_json(text)
        assert result == {"gene": "BRCA1", "impact": "HIGH"}

    def test_markdown_fenced_no_lang(self):
        text = '```\n{"gene": "TP53"}\n```'
        result = extract_json(text)
        assert result == {"gene": "TP53"}

    def test_trailing_comma_repair(self):
        text = '{"a": 1, "b": 2,}'
        result = extract_json(text)
        assert result == {"a": 1, "b": 2}

    def test_trailing_comma_array_repair(self):
        text = "[1, 2, 3,]"
        result = extract_json(text)
        assert result == [1, 2, 3]

    def test_json_embedded_in_text(self):
        text = 'Here is the result: {"answer": 42} as you can see.'
        result = extract_json(text)
        assert result == {"answer": 42}

    def test_invalid_json_raises(self):
        with pytest.raises(ValueError, match="Could not extract"):
            extract_json("This is not JSON at all.")

    def test_empty_string_raises(self):
        with pytest.raises(ValueError):
            extract_json("")

    def test_nested_json(self):
        text = '{"patient": {"id": "HG002", "variants": [1, 2]}}'
        result = extract_json(text)
        assert result["patient"]["id"] == "HG002"
        assert result["patient"]["variants"] == [1, 2]


class TestTokenBudget:
    def test_unlimited_budget(self):
        budget = _TokenBudget(daily_limit=0)
        assert budget.check() is True
        assert budget.remaining == -1

    def test_budget_tracking(self):
        budget = _TokenBudget(daily_limit=1000)
        budget.record(500)
        assert budget.check() is True
        assert budget.remaining == 500

    def test_budget_exhausted(self):
        budget = _TokenBudget(daily_limit=100)
        budget.record(100)
        assert budget.check() is False
        assert budget.remaining == 0

    def test_configure_changes_limit(self):
        budget = _TokenBudget(daily_limit=100)
        budget.configure(5000)
        assert budget._daily_limit == 5000

    def test_budget_over_records(self):
        budget = _TokenBudget(daily_limit=100)
        budget.record(150)
        assert budget.remaining == 0
        assert budget.check() is False


class TestResponseCache:
    def test_cache_miss(self):
        cache = _ResponseCache(max_size=10, ttl=60.0)
        result = cache.get("prompt", "system", "model", 0.7)
        assert result is None

    def test_cache_hit(self):
        cache = _ResponseCache(max_size=10, ttl=60.0)
        cache.put("prompt", "system", "model", 0.7, "response text")
        result = cache.get("prompt", "system", "model", 0.7)
        assert result == "response text"

    def test_different_params_miss(self):
        cache = _ResponseCache(max_size=10, ttl=60.0)
        cache.put("prompt", "system", "model", 0.7, "response text")
        # Different temperature
        result = cache.get("prompt", "system", "model", 0.3)
        assert result is None

    def test_disabled_when_ttl_zero(self):
        cache = _ResponseCache(max_size=10, ttl=0)
        cache.put("prompt", "system", "model", 0.7, "response")
        result = cache.get("prompt", "system", "model", 0.7)
        assert result is None

    def test_clear(self):
        cache = _ResponseCache(max_size=10, ttl=60.0)
        cache.put("prompt", "system", "model", 0.7, "response")
        cache.clear()
        result = cache.get("prompt", "system", "model", 0.7)
        assert result is None

    def test_lru_eviction(self):
        cache = _ResponseCache(max_size=2, ttl=60.0)
        cache.put("p1", "s", "m", 0.7, "r1")
        cache.put("p2", "s", "m", 0.7, "r2")
        cache.put("p3", "s", "m", 0.7, "r3")  # evicts p1
        assert cache.get("p1", "s", "m", 0.7) is None
        assert cache.get("p3", "s", "m", 0.7) == "r3"


class TestLLMClientFactory:
    def test_unknown_provider_raises(self):
        with pytest.raises(ValueError, match="Unknown LLM provider"):
            LLMClientFactory.create(provider="gpt4all")

    def test_anthropic_without_key_raises(self):
        """Anthropic client requires API key."""
        with pytest.raises((ValueError, ImportError)):
            LLMClientFactory.create(provider="anthropic", api_key=None)

    def test_openai_without_key_raises(self):
        """OpenAI client requires API key."""
        with pytest.raises((ValueError, ImportError)):
            LLMClientFactory.create(provider="openai", api_key=None)


class TestBaseLLMClient:
    def test_is_abstract(self):
        with pytest.raises(TypeError):
            BaseLLMClient()

    def test_provider_attribute(self):
        assert BaseLLMClient.provider == "unknown"
