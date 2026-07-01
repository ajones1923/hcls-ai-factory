"""Tests for hcls_common.config — HCLSSettings and get_settings()."""

import os
import pytest
from unittest.mock import patch

from hcls_common.config import HCLSSettings, get_settings


class TestHCLSSettingsDefaults:
    """Verify default values are sensible."""

    def test_default_service_name(self):
        s = HCLSSettings()
        assert s.service_name == "hcls-ai-factory"

    def test_default_milvus(self):
        s = HCLSSettings()
        assert s.milvus_host == "localhost"
        assert s.milvus_port == 19530
        assert s.milvus_collection == "genomic_evidence"

    def test_default_embedding(self):
        s = HCLSSettings()
        assert s.embedding_model == "BAAI/bge-small-en-v1.5"
        assert s.embedding_dimension == 384
        assert s.embedding_provider == "local_bge"

    def test_default_llm(self):
        s = HCLSSettings()
        assert s.llm_provider == "ollama"
        assert s.llm_temperature == 0.7
        assert s.llm_max_tokens == 4096

    def test_default_rag(self):
        s = HCLSSettings()
        assert s.rag_top_k == 10
        assert 0 <= s.rag_score_threshold <= 1.0

    def test_default_nim_urls(self):
        s = HCLSSettings()
        assert "8001" in s.nim_molmim_url
        assert "8002" in s.nim_diffdock_url


class TestHCLSSettingsValidation:
    """Verify field validators reject bad input."""

    def test_invalid_llm_provider(self):
        with pytest.raises(ValueError, match="llm_provider"):
            HCLSSettings(llm_provider="gpt4all")

    def test_valid_llm_providers(self):
        for provider in ("anthropic", "openai", "ollama", "vllm"):
            s = HCLSSettings(llm_provider=provider)
            assert s.llm_provider == provider

    def test_invalid_embedding_provider(self):
        with pytest.raises(ValueError, match="embedding_provider"):
            HCLSSettings(embedding_provider="cohere")

    def test_valid_embedding_providers(self):
        for provider in ("local_bge", "tei", "openai"):
            s = HCLSSettings(embedding_provider=provider)
            assert s.embedding_provider == provider

    def test_temperature_range(self):
        with pytest.raises(Exception):
            HCLSSettings(llm_temperature=3.0)

    def test_temperature_lower_bound(self):
        with pytest.raises(Exception):
            HCLSSettings(llm_temperature=-0.1)


class TestHCLSSettingsEnvOverride:
    """Verify environment variable overrides work."""

    def test_env_override_milvus_host(self):
        with patch.dict(os.environ, {"HCLS_MILVUS_HOST": "milvus-prod"}):
            s = HCLSSettings()
            assert s.milvus_host == "milvus-prod"

    def test_env_override_llm_model(self):
        with patch.dict(os.environ, {"HCLS_LLM_MODEL": "claude-sonnet-4-6"}):
            s = HCLSSettings()
            assert s.llm_model == "claude-sonnet-4-6"


class TestGetSettings:
    """Verify singleton caching."""

    def test_get_settings_returns_instance(self):
        get_settings.cache_clear()
        s = get_settings()
        assert isinstance(s, HCLSSettings)

    def test_get_settings_is_cached(self):
        get_settings.cache_clear()
        s1 = get_settings()
        s2 = get_settings()
        assert s1 is s2

    def test_cache_clear_reloads(self):
        get_settings.cache_clear()
        s1 = get_settings()
        get_settings.cache_clear()
        s2 = get_settings()
        # Different objects after cache clear
        assert s1 is not s2
