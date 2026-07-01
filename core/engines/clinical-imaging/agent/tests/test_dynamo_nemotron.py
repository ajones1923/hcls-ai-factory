"""Tests for NVIDIA Dynamo optimizer and Nemotron Nano client.

All tests work without dynamo or nemotron packages installed -- they
verify graceful degradation, mock behavior, and configuration defaults.

Author: Adam Jones
Date: April 2026
"""

from unittest.mock import patch

import pytest

from config.settings import ImagingSettings
from src.inference.dynamo_config import (
    MODEL_PRESETS,
    DynamoConfig,
    DynamoOptimizer,
    DynamoStatus,
)
from src.nim.base import BaseNIMClient
from src.nim.nemotron_nano_client import NemotronNanoClient


# ===================================================================
# DynamoOptimizer — Graceful Fallback
# ===================================================================


class TestDynamoGracefulFallback:
    """Dynamo works correctly when the dynamo package is not installed."""

    def test_dynamo_graceful_fallback(self):
        """Optimizer initializes without error when dynamo is missing."""
        with patch.dict("sys.modules", {"dynamo": None}):
            optimizer = DynamoOptimizer()
            assert optimizer._dynamo_available is False

    def test_dynamo_optimize_request_without_dynamo(self):
        """optimize_request returns empty dict when dynamo unavailable."""
        with patch.dict("sys.modules", {"dynamo": None}):
            optimizer = DynamoOptimizer()
            hints = optimizer.optimize_request("Analyze this CT scan", max_tokens=512)
            assert hints == {}


class TestDynamoStatusWhenUnavailable:
    """Status reporting when Dynamo is not available."""

    def test_dynamo_status_when_unavailable(self):
        """get_status returns available=False when dynamo not installed."""
        with patch.dict("sys.modules", {"dynamo": None}):
            optimizer = DynamoOptimizer()
            status = optimizer.get_status()
            assert isinstance(status, DynamoStatus)
            assert status.available is False
            assert status.scheduler == "none"
            assert status.active_requests == 0
            assert status.kv_cache_usage_gb == 0.0

    def test_dynamo_status_when_disabled(self):
        """get_status returns available=False when config disabled."""
        config = DynamoConfig(enabled=False)
        optimizer = DynamoOptimizer(config=config)
        status = optimizer.get_status()
        assert status.available is False


class TestDynamoConfigDefaults:
    """DynamoConfig has correct default values."""

    def test_dynamo_config_defaults(self):
        """Default config values match specification."""
        config = DynamoConfig()
        assert config.enabled is False
        assert config.prefill_batch_size == 4
        assert config.decode_batch_size == 8
        assert config.kv_cache_max_gb == 16.0
        assert config.scheduler == "round_robin"

    def test_dynamo_config_custom_values(self):
        """Config accepts custom values."""
        config = DynamoConfig(
            enabled=True,
            prefill_batch_size=16,
            decode_batch_size=32,
            kv_cache_max_gb=8.0,
            scheduler="load_balanced",
        )
        assert config.enabled is True
        assert config.prefill_batch_size == 16
        assert config.kv_cache_max_gb == 8.0


# ===================================================================
# DynamoOptimizer — KV Cache Estimation
# ===================================================================


class TestDynamoKVCacheEstimate:
    """KV cache size estimation for different models."""

    def test_dynamo_kv_cache_estimate_8b(self):
        """KV cache estimate for 8B model is in a reasonable range."""
        optimizer = DynamoOptimizer()
        # 4096 context length, 8B model
        kv_gb = optimizer.estimate_kv_cache_size(
            context_length=4096, model_params_b=8.0
        )
        # 8B model at 4K context should be well under 16 GB
        assert 0.0 < kv_gb < 16.0
        # Should be a meaningful amount (at least a few hundred MB)
        assert kv_gb > 0.01

    def test_dynamo_kv_cache_estimate_nano(self):
        """KV cache estimate for Nano is smaller than 8B."""
        optimizer = DynamoOptimizer()
        kv_8b = optimizer.estimate_kv_cache_size(
            context_length=4096, model_params_b=8.0
        )
        kv_nano = optimizer.estimate_kv_cache_size(
            context_length=4096, model_params_b=4.0
        )
        assert kv_nano < kv_8b

    def test_dynamo_kv_cache_scales_with_context(self):
        """KV cache grows with context length."""
        optimizer = DynamoOptimizer()
        kv_short = optimizer.estimate_kv_cache_size(
            context_length=1024, model_params_b=8.0
        )
        kv_long = optimizer.estimate_kv_cache_size(
            context_length=8192, model_params_b=8.0
        )
        assert kv_long > kv_short


# ===================================================================
# DynamoOptimizer — Recommended Configs
# ===================================================================


class TestDynamoRecommendedConfig:
    """Preset configurations for known models."""

    def test_dynamo_recommended_config_llama(self):
        """Llama-3-8B preset has correct values."""
        optimizer = DynamoOptimizer()
        config = optimizer.get_recommended_config("llama-3-8b")
        assert isinstance(config, DynamoConfig)
        assert config.prefill_batch_size == 4
        assert config.decode_batch_size == 8
        assert config.kv_cache_max_gb == 16.0
        assert config.scheduler == "prefill_priority"

    def test_dynamo_recommended_config_nano(self):
        """Nemotron Nano preset has correct values."""
        optimizer = DynamoOptimizer()
        config = optimizer.get_recommended_config("nemotron-nano")
        assert isinstance(config, DynamoConfig)
        assert config.prefill_batch_size == 8
        assert config.decode_batch_size == 16
        assert config.kv_cache_max_gb == 4.0
        assert config.scheduler == "round_robin"

    def test_dynamo_recommended_config_unknown(self):
        """Unknown model returns default config."""
        optimizer = DynamoOptimizer()
        config = optimizer.get_recommended_config("unknown-model-xyz")
        assert isinstance(config, DynamoConfig)
        # Should use defaults
        assert config.prefill_batch_size == 4
        assert config.decode_batch_size == 8


# ===================================================================
# NemotronNanoClient — Mock Generation
# ===================================================================


class TestNemotronNanoMockGenerate:
    """Mock generation returns clinical text."""

    def test_nemotron_nano_mock_generate(self):
        """Mock mode returns a non-empty clinical response."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=True
        )
        response = client.generate("What is the standard CT chest protocol?")
        assert isinstance(response, str)
        assert len(response) > 50
        assert "Nemotron Nano" in response

    def test_nemotron_nano_mock_generate_protocol(self):
        """Mock mode returns protocol-specific content for CT queries."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=True
        )
        response = client.generate("What is the protocol for CT chest?")
        assert "protocol" in response.lower() or "CT" in response

    def test_nemotron_nano_mock_generate_measurement(self):
        """Mock mode returns measurement-specific content."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=True
        )
        response = client.generate("What is the normal range for liver size?")
        assert "reference" in response.lower() or "range" in response.lower()


# ===================================================================
# NemotronNanoClient — Query Routing
# ===================================================================


class TestNemotronNanoQueryRouting:
    """Query classification for routine vs complex queries."""

    def test_nemotron_nano_routine_query(self):
        """Protocol lookup is classified as routine."""
        assert NemotronNanoClient.is_routine_query(
            "What is the protocol for CT chest?"
        ) is True

    def test_nemotron_nano_routine_query_define(self):
        """Definition query is classified as routine."""
        assert NemotronNanoClient.is_routine_query(
            "Define Hounsfield unit"
        ) is True

    def test_nemotron_nano_routine_query_device(self):
        """Device query is classified as routine."""
        assert NemotronNanoClient.is_routine_query(
            "List FDA cleared CT scanners"
        ) is True

    def test_nemotron_nano_complex_query(self):
        """Comparative query is classified as complex."""
        assert NemotronNanoClient.is_routine_query(
            "Compare CT vs MRI for hemorrhage detection"
        ) is False

    def test_nemotron_nano_complex_query_differential(self):
        """Differential diagnosis is classified as complex."""
        assert NemotronNanoClient.is_routine_query(
            "What is the differential diagnosis for a lung nodule?"
        ) is False

    def test_nemotron_nano_complex_query_genomic(self):
        """Genomic integration query is classified as complex."""
        assert NemotronNanoClient.is_routine_query(
            "Correlate the genomic findings with the MRI features"
        ) is False

    def test_nemotron_nano_unknown_defaults_complex(self):
        """Unrecognized queries default to complex (route to larger model)."""
        assert NemotronNanoClient.is_routine_query(
            "Tell me something interesting"
        ) is False


# ===================================================================
# NemotronNanoClient — Status & Inheritance
# ===================================================================


class TestNemotronNanoStatus:
    """Status reporting for Nemotron Nano client."""

    def test_nemotron_nano_status(self):
        """Returns 'mock' when service is unavailable but mock enabled."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=True
        )
        status = client.get_status()
        assert status == "mock"

    def test_nemotron_nano_status_unavailable(self):
        """Returns 'unavailable' when mock disabled and service down."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=False
        )
        status = client.get_status()
        assert status == "unavailable"


class TestNemotronNanoInheritance:
    """Verify proper BaseNIMClient inheritance."""

    def test_nemotron_nano_client_inherits_base(self):
        """NemotronNanoClient is a proper subclass of BaseNIMClient."""
        assert issubclass(NemotronNanoClient, BaseNIMClient)

    def test_nemotron_nano_has_required_methods(self):
        """Client implements all required abstract methods."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=True
        )
        # BaseNIMClient requires _mock_response
        assert hasattr(client, "_mock_response")
        assert callable(client._mock_response)
        # Should also have health_check from base
        assert hasattr(client, "health_check")
        assert hasattr(client, "is_available")
        assert hasattr(client, "get_status")

    def test_nemotron_nano_service_name(self):
        """Service name is set correctly."""
        client = NemotronNanoClient(
            base_url="http://localhost:8538", mock_enabled=True
        )
        assert client.service_name == "nemotron-nano"


# ===================================================================
# Settings — Dynamo & Nemotron Defaults
# ===================================================================


class TestSettingsDefaults:
    """Dynamo and Nemotron settings are disabled by default."""

    def test_settings_dynamo_defaults(self):
        """Dynamo settings default to disabled."""
        s = ImagingSettings()
        assert s.DYNAMO_ENABLED is False
        assert s.DYNAMO_PREFILL_BATCH == 4
        assert s.DYNAMO_DECODE_BATCH == 8
        assert s.DYNAMO_KV_CACHE_GB == 16.0

    def test_settings_nemotron_defaults(self):
        """Nemotron Nano settings default to disabled."""
        s = ImagingSettings()
        assert s.NIM_NEMOTRON_NANO_URL == "http://localhost:8538"
        assert s.NEMOTRON_NANO_ENABLED is False
        assert s.NEMOTRON_NANO_ROUTE_ROUTINE is True
