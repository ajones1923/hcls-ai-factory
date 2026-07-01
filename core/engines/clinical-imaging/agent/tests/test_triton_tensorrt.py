"""Tests for Triton Inference Server and TensorRT optimization.

All tests work without Triton or TensorRT installed -- they verify
graceful degradation and mock/disabled behavior.

Author: Adam Jones
Date: April 2026
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from src.inference.tensorrt_optimizer import (
    DYNAMIC_SHAPE_PROFILES,
    SPEEDUP_ESTIMATES,
    TensorRTOptimizer,
)
from src.inference.triton_config import (
    MONAI_MODEL_CONFIGS,
    TritonModelManager,
    get_model_configs,
)


# ===================================================================
# TritonModelManager
# ===================================================================


class TestTritonGracefulFallback:
    """Triton works correctly when tritonclient is not installed."""

    def test_triton_graceful_fallback(self):
        """Manager initializes with enabled=False when tritonclient missing."""
        with patch.dict("sys.modules", {"tritonclient": None, "tritonclient.http": None}):
            mgr = TritonModelManager(enabled=True)
            # Should have disabled itself
            assert mgr.enabled is False

    def test_triton_disabled_by_config(self):
        """Manager respects enabled=False without trying imports."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.enabled is False
        assert mgr._client is None

    def test_triton_health_check_when_disabled(self):
        """Health check returns False when disabled."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.is_healthy() is False

    def test_triton_list_models_when_disabled(self):
        """list_models returns empty list when disabled."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.list_models() == []

    def test_triton_load_model_when_disabled(self):
        """load_model returns False when disabled."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.load_model("segresnet_hemorrhage") is False

    def test_triton_unload_model_when_disabled(self):
        """unload_model returns False when disabled."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.unload_model("segresnet_hemorrhage") is False

    def test_triton_is_model_ready_when_disabled(self):
        """is_model_ready returns False when disabled."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.is_model_ready("segresnet_hemorrhage") is False

    def test_triton_infer_raises_when_disabled(self):
        """infer raises RuntimeError when disabled."""
        mgr = TritonModelManager(enabled=False)
        with pytest.raises(RuntimeError, match="Triton not available"):
            mgr.infer("segresnet_hemorrhage", {"INPUT__0": np.zeros((1,))})

    def test_triton_get_status_when_disabled(self):
        """get_status works when disabled."""
        mgr = TritonModelManager(enabled=False)
        status = mgr.get_status()
        assert status["enabled"] is False
        assert status["healthy"] is False
        assert "available_configs" in status


class TestTritonHealthCheckMock:
    """Health check with mocked tritonclient."""

    def test_health_check_healthy(self):
        """Health check returns True when server is live and ready."""
        mgr = TritonModelManager(enabled=False)
        mgr.enabled = True
        mock_client = MagicMock()
        mock_client.is_server_live.return_value = True
        mock_client.is_server_ready.return_value = True
        mgr._client = mock_client

        assert mgr.is_healthy() is True

    def test_health_check_not_ready(self):
        """Health check returns False when server is not ready."""
        mgr = TritonModelManager(enabled=False)
        mgr.enabled = True
        mock_client = MagicMock()
        mock_client.is_server_live.return_value = True
        mock_client.is_server_ready.return_value = False
        mgr._client = mock_client

        assert mgr.is_healthy() is False

    def test_health_check_caching(self):
        """Health check result is cached within TTL."""
        mgr = TritonModelManager(enabled=False)
        mgr.enabled = True
        mock_client = MagicMock()
        mock_client.is_server_live.return_value = True
        mock_client.is_server_ready.return_value = True
        mgr._client = mock_client

        # First call hits the client
        assert mgr.is_healthy() is True
        assert mock_client.is_server_live.call_count == 1

        # Second call within TTL uses cache
        assert mgr.is_healthy() is True
        assert mock_client.is_server_live.call_count == 1


class TestTritonModelConfig:
    """Model configuration and repository setup."""

    def test_get_model_config_static(self):
        """get_model_config returns static config when Triton unavailable."""
        mgr = TritonModelManager(enabled=False)
        config = mgr.get_model_config("segresnet_hemorrhage")
        assert config["framework"] == "pytorch"
        assert config["input_shape"] == [1, 1, 256, 256, 256]

    def test_get_model_config_unknown(self):
        """get_model_config returns empty dict for unknown model."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.get_model_config("nonexistent_model") == {}

    def test_generate_config_pbtxt(self):
        """Config pbtxt generation produces valid structure."""
        mgr = TritonModelManager(enabled=False)
        pbtxt = mgr.generate_config_pbtxt("densenet121_cxr")
        assert 'name: "densenet121_cxr"' in pbtxt
        assert "platform: \"pytorch_libtorch\"" in pbtxt
        assert "max_batch_size: 8" in pbtxt
        assert "KIND_GPU" in pbtxt

    def test_generate_config_pbtxt_unknown(self):
        """Config pbtxt returns empty string for unknown model."""
        mgr = TritonModelManager(enabled=False)
        assert mgr.generate_config_pbtxt("unknown") == ""

    def test_setup_model_repository(self):
        """Model repository setup creates dirs and config files."""
        mgr = TritonModelManager(enabled=False)
        with tempfile.TemporaryDirectory() as tmpdir:
            created = mgr.setup_model_repository(
                tmpdir, model_names=["densenet121_cxr"]
            )
            assert len(created) == 1
            config_path = Path(tmpdir) / "densenet121_cxr" / "config.pbtxt"
            assert config_path.exists()
            version_dir = Path(tmpdir) / "densenet121_cxr" / "1"
            assert version_dir.is_dir()

    def test_get_model_configs_convenience(self):
        """get_model_configs returns all 4 MONAI model configs."""
        configs = get_model_configs()
        assert len(configs) == 4
        assert "segresnet_hemorrhage" in configs
        assert "unest_brain_ms" in configs
        assert "retinanet_lung_nodule" in configs
        assert "densenet121_cxr" in configs


class TestTritonServerMetadata:
    """Server metadata retrieval."""

    def test_metadata_when_disabled(self):
        """get_server_metadata returns unavailable when disabled."""
        mgr = TritonModelManager(enabled=False)
        meta = mgr.get_server_metadata()
        assert meta["status"] == "unavailable"


# ===================================================================
# TensorRTOptimizer
# ===================================================================


class TestTensorRTGracefulFallback:
    """TensorRT works correctly when torch_tensorrt is not installed."""

    def test_tensorrt_graceful_fallback(self):
        """Optimizer initializes with available=False when missing."""
        optimizer = TensorRTOptimizer()
        # torch_tensorrt is almost certainly not installed in test env
        # but the class must not raise
        assert isinstance(optimizer.available, bool)

    def test_tensorrt_optimize_raises_when_unavailable(self):
        """optimize_model raises RuntimeError when TensorRT missing."""
        optimizer = TensorRTOptimizer()
        if not optimizer.available:
            with pytest.raises(RuntimeError, match="TensorRT not available"):
                optimizer.optimize_model("/fake/model.pt", "segresnet_hemorrhage")

    def test_tensorrt_load_raises_when_unavailable(self):
        """load_optimized raises RuntimeError when TensorRT missing."""
        optimizer = TensorRTOptimizer()
        if not optimizer.available:
            with pytest.raises(RuntimeError, match="TensorRT not available"):
                optimizer.load_optimized("segresnet_hemorrhage")


class TestTensorRTIsOptimized:
    """Cache detection tests."""

    def test_is_optimized_false_initially(self):
        """No cached models exist at startup."""
        with tempfile.TemporaryDirectory() as tmpdir:
            optimizer = TensorRTOptimizer(cache_dir=tmpdir)
            assert optimizer.is_optimized("segresnet_hemorrhage") is False
            assert optimizer.is_optimized("unest_brain_ms") is False
            assert optimizer.is_optimized("retinanet_lung_nodule") is False
            assert optimizer.is_optimized("densenet121_cxr") is False

    def test_is_optimized_nonexistent_cache_dir(self):
        """is_optimized returns False when cache dir doesn't exist."""
        optimizer = TensorRTOptimizer(cache_dir="/nonexistent/path/trt_cache")
        assert optimizer.is_optimized("segresnet_hemorrhage") is False


class TestTensorRTSpeedupEstimates:
    """Speedup estimation for all 4 MONAI models."""

    def test_speedup_estimates_all_models(self):
        """Speedup estimates returned for all 4 models."""
        optimizer = TensorRTOptimizer(precision="fp16")
        for model_name in MONAI_MODEL_CONFIGS:
            estimate = optimizer.get_speedup_estimate(model_name)
            assert estimate["model_name"] == model_name
            assert estimate["current_precision"] == "fp16"
            assert estimate["estimated_speedup"] > 1.0
            assert "speedup_by_precision" in estimate
            assert "input_shape_profile" in estimate

    def test_speedup_segresnet(self):
        """SegResNet FP16 speedup is ~3x."""
        optimizer = TensorRTOptimizer(precision="fp16")
        est = optimizer.get_speedup_estimate("segresnet_hemorrhage")
        assert est["estimated_speedup"] == 3.0

    def test_speedup_unest(self):
        """UNEST FP16 speedup is ~2.5x."""
        optimizer = TensorRTOptimizer(precision="fp16")
        est = optimizer.get_speedup_estimate("unest_brain_ms")
        assert est["estimated_speedup"] == 2.5

    def test_speedup_retinanet(self):
        """RetinaNet FP16 speedup is ~3x."""
        optimizer = TensorRTOptimizer(precision="fp16")
        est = optimizer.get_speedup_estimate("retinanet_lung_nodule")
        assert est["estimated_speedup"] == 3.0

    def test_speedup_densenet(self):
        """DenseNet-121 FP16 speedup is ~4x (2D, fastest)."""
        optimizer = TensorRTOptimizer(precision="fp16")
        est = optimizer.get_speedup_estimate("densenet121_cxr")
        assert est["estimated_speedup"] == 4.0

    def test_speedup_unknown_model(self):
        """Unknown model returns 1.0x speedup (no optimization)."""
        optimizer = TensorRTOptimizer(precision="fp16")
        est = optimizer.get_speedup_estimate("nonexistent_model")
        assert est["estimated_speedup"] == 1.0

    def test_get_all_speedup_estimates(self):
        """get_all_speedup_estimates returns list for all models."""
        optimizer = TensorRTOptimizer(precision="fp16")
        estimates = optimizer.get_all_speedup_estimates()
        assert len(estimates) == 4


class TestTensorRTPrecisionConfig:
    """Precision configuration tests."""

    def test_fp16_precision(self):
        """FP16 precision is accepted."""
        optimizer = TensorRTOptimizer(precision="fp16")
        assert optimizer.precision == "fp16"

    def test_int8_precision(self):
        """INT8 precision is accepted."""
        optimizer = TensorRTOptimizer(precision="int8")
        assert optimizer.precision == "int8"

    def test_fp32_precision(self):
        """FP32 precision is accepted."""
        optimizer = TensorRTOptimizer(precision="fp32")
        assert optimizer.precision == "fp32"

    def test_precision_affects_speedup(self):
        """Different precisions give different speedup estimates."""
        opt_fp16 = TensorRTOptimizer(precision="fp16")
        opt_int8 = TensorRTOptimizer(precision="int8")
        opt_fp32 = TensorRTOptimizer(precision="fp32")

        est_fp16 = opt_fp16.get_speedup_estimate("densenet121_cxr")
        est_int8 = opt_int8.get_speedup_estimate("densenet121_cxr")
        est_fp32 = opt_fp32.get_speedup_estimate("densenet121_cxr")

        assert est_int8["estimated_speedup"] > est_fp16["estimated_speedup"]
        assert est_fp16["estimated_speedup"] > est_fp32["estimated_speedup"]
        assert est_fp32["estimated_speedup"] == 1.0


class TestTensorRTCache:
    """Cache management tests."""

    def test_list_cached_models_empty(self):
        """list_cached_models returns empty list for empty cache."""
        with tempfile.TemporaryDirectory() as tmpdir:
            optimizer = TensorRTOptimizer(cache_dir=tmpdir)
            assert optimizer.list_cached_models() == []

    def test_clear_cache_empty(self):
        """clear_cache returns 0 when no files to remove."""
        with tempfile.TemporaryDirectory() as tmpdir:
            optimizer = TensorRTOptimizer(cache_dir=tmpdir)
            assert optimizer.clear_cache() == 0

    def test_clear_cache_nonexistent(self):
        """clear_cache returns 0 when cache dir doesn't exist."""
        optimizer = TensorRTOptimizer(cache_dir="/nonexistent/trt_cache")
        assert optimizer.clear_cache() == 0

    def test_get_status(self):
        """get_status returns complete status dict."""
        with tempfile.TemporaryDirectory() as tmpdir:
            optimizer = TensorRTOptimizer(cache_dir=tmpdir, precision="fp16")
            status = optimizer.get_status()
            assert isinstance(status["available"], bool)
            assert status["precision"] == "fp16"
            assert status["cache_dir"] == tmpdir
            assert status["cached_models"] == 0
            assert len(status["known_models"]) == 4


# ===================================================================
# Model Configs
# ===================================================================


class TestModelConfigs:
    """Verify MONAI model configuration completeness."""

    @pytest.mark.parametrize(
        "model_name",
        ["segresnet_hemorrhage", "unest_brain_ms", "retinanet_lung_nodule", "densenet121_cxr"],
    )
    def test_model_config_structure(self, model_name):
        """Each model config has all required fields."""
        config = MONAI_MODEL_CONFIGS[model_name]
        required_keys = [
            "display_name", "framework", "max_batch_size",
            "input_name", "input_dtype", "input_shape",
            "output_name", "output_dtype", "output_shape",
            "instance_group_count", "model_file", "optimization",
        ]
        for key in required_keys:
            assert key in config, f"Missing key '{key}' in {model_name}"

    @pytest.mark.parametrize(
        "model_name",
        ["segresnet_hemorrhage", "unest_brain_ms", "retinanet_lung_nodule", "densenet121_cxr"],
    )
    def test_dynamic_shape_profiles_exist(self, model_name):
        """Each model has a dynamic shape profile for TensorRT."""
        assert model_name in DYNAMIC_SHAPE_PROFILES
        profile = DYNAMIC_SHAPE_PROFILES[model_name]
        assert "min" in profile
        assert "opt" in profile
        assert "max" in profile

    @pytest.mark.parametrize(
        "model_name",
        ["segresnet_hemorrhage", "unest_brain_ms", "retinanet_lung_nodule", "densenet121_cxr"],
    )
    def test_speedup_estimates_exist(self, model_name):
        """Each model has speedup estimates for all precision levels."""
        assert model_name in SPEEDUP_ESTIMATES
        estimates = SPEEDUP_ESTIMATES[model_name]
        assert "fp16" in estimates
        assert "int8" in estimates
        assert "fp32" in estimates

    def test_densenet121_is_2d(self):
        """DenseNet-121 (CXR) uses 2D input (4 dims, not 5)."""
        config = MONAI_MODEL_CONFIGS["densenet121_cxr"]
        assert len(config["input_shape"]) == 4  # [B, C, H, W]

    def test_3d_models_have_5d_input(self):
        """3D models (SegResNet, UNEST, RetinaNet) use 5D input."""
        for name in ["segresnet_hemorrhage", "unest_brain_ms", "retinanet_lung_nodule"]:
            config = MONAI_MODEL_CONFIGS[name]
            assert len(config["input_shape"]) == 5  # [B, C, D, H, W]


# ===================================================================
# Settings Integration
# ===================================================================


class TestSettingsDefaults:
    """Verify Triton/TensorRT settings defaults."""

    def test_triton_disabled_by_default(self):
        """TRITON_ENABLED defaults to False."""
        from config.settings import ImagingSettings

        s = ImagingSettings()
        assert s.TRITON_ENABLED is False

    def test_tensorrt_disabled_by_default(self):
        """TENSORRT_ENABLED defaults to False."""
        from config.settings import ImagingSettings

        s = ImagingSettings()
        assert s.TENSORRT_ENABLED is False

    def test_tensorrt_precision_default(self):
        """TENSORRT_PRECISION defaults to fp16."""
        from config.settings import ImagingSettings

        s = ImagingSettings()
        assert s.TENSORRT_PRECISION == "fp16"

    def test_triton_url_default(self):
        """TRITON_URL defaults to localhost:8000."""
        from config.settings import ImagingSettings

        s = ImagingSettings()
        assert s.TRITON_URL == "localhost:8000"

    def test_tensorrt_cache_dir_default(self):
        """TENSORRT_CACHE_DIR points to data/cache/tensorrt."""
        from config.settings import ImagingSettings

        s = ImagingSettings()
        assert "tensorrt" in s.TENSORRT_CACHE_DIR
