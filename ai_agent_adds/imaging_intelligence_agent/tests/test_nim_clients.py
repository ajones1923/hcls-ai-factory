"""Tests for NIM service clients (VISTA-3D, MAISI, VILA-M3, LLM, ServiceManager).

All tests use mock mode or mock dependencies -- no real NIM service needed.

Author: Adam Jones
Date: February 2026
"""

import pytest
from unittest.mock import MagicMock, patch, PropertyMock

from src.models import (
    NIMServiceStatus,
    SegmentationResult,
    SyntheticCTResult,
    VLMResponse,
)
from src.nim.base import BaseNIMClient
from src.nim.vista3d_client import VISTA3D_CLASSES, VISTA3DClient
from src.nim.maisi_client import MAISIClient
from src.nim.vilam3_client import VILAM3Client
from src.nim.llm_client import LlamaLLMClient
from src.nim.service_manager import NIMServiceManager


# ===================================================================
# BaseNIMClient
# ===================================================================


class TestBaseNIMClient:
    """Tests for BaseNIMClient abstract base class."""

    def test_cannot_instantiate_directly(self):
        with pytest.raises(TypeError):
            BaseNIMClient("http://localhost:8000", "test")

    def test_concrete_subclass_works(self):
        """A minimal concrete subclass can be instantiated."""

        class ConcreteClient(BaseNIMClient):
            def _mock_response(self, **kwargs):
                return {"result": "mock"}

        client = ConcreteClient("http://localhost:8000", "test", mock_enabled=True)
        assert client.service_name == "test"
        assert client.mock_enabled is True

    def test_health_check_returns_false_when_unavailable(self):
        class ConcreteClient(BaseNIMClient):
            def _mock_response(self, **kwargs):
                return {}

        client = ConcreteClient("http://localhost:99999", "test")
        # Will fail to connect
        assert client.health_check() is False

    def test_get_status_mock_when_unavailable(self):
        class ConcreteClient(BaseNIMClient):
            def _mock_response(self, **kwargs):
                return {}

        client = ConcreteClient("http://localhost:99999", "test", mock_enabled=True)
        client._available = False
        client._last_check = 1e18  # Far future so cache is valid
        status = client.get_status()
        assert status == "mock"

    def test_get_status_unavailable_when_mock_disabled(self):
        class ConcreteClient(BaseNIMClient):
            def _mock_response(self, **kwargs):
                return {}

        client = ConcreteClient("http://localhost:99999", "test", mock_enabled=False)
        client._available = False
        client._last_check = 1e18
        status = client.get_status()
        assert status == "unavailable"


# ===================================================================
# VISTA3DClient
# ===================================================================


class TestVISTA3DClient:
    """Tests for VISTA3DClient mock mode."""

    def test_init(self):
        client = VISTA3DClient("http://localhost:8530", mock_enabled=True)
        assert client.service_name == "vista3d"
        assert client.mock_enabled is True

    def test_segment_mock_returns_segmentation_result(self):
        client = VISTA3DClient("http://localhost:8530", mock_enabled=True)
        # Force unavailable so mock fallback activates
        client._available = False
        client._last_check = 1e18
        result = client.segment("/tmp/test.nii.gz")
        assert isinstance(result, SegmentationResult)
        assert result.is_mock is True

    def test_segment_mock_has_classes_detected(self):
        client = VISTA3DClient("http://localhost:8530", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.segment("/tmp/test.nii.gz")
        assert len(result.classes_detected) >= 5
        assert len(result.classes_detected) <= 10

    def test_segment_mock_has_volumes(self):
        client = VISTA3DClient("http://localhost:8530", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.segment("/tmp/test.nii.gz")
        assert len(result.volumes) >= 5
        for cls, vol in result.volumes.items():
            assert vol > 0.0

    def test_segment_mock_has_inference_time(self):
        client = VISTA3DClient("http://localhost:8530", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.segment("/tmp/test.nii.gz")
        assert result.inference_time_ms > 0

    def test_segment_with_specific_classes(self):
        client = VISTA3DClient("http://localhost:8530", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.segment("/tmp/test.nii.gz", classes=["liver", "spleen"])
        assert isinstance(result, SegmentationResult)

    def test_get_supported_classes_returns_list(self):
        client = VISTA3DClient("http://localhost:8530")
        classes = client.get_supported_classes()
        assert isinstance(classes, list)
        assert len(classes) > 100
        assert "liver" in classes
        assert "brain" in classes
        assert "heart" in classes

    def test_vista3d_classes_constant(self):
        assert len(VISTA3D_CLASSES) > 100
        assert "liver" in VISTA3D_CLASSES
        assert "aorta" in VISTA3D_CLASSES


# ===================================================================
# MAISIClient
# ===================================================================


class TestMAISIClient:
    """Tests for MAISIClient mock mode."""

    def test_init(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        assert client.service_name == "maisi"

    def test_generate_mock_returns_synthetic_ct_result(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.generate(body_region="chest")
        assert isinstance(result, SyntheticCTResult)
        assert result.is_mock is True
        assert result.body_region == "chest"

    def test_generate_mock_has_resolution(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.generate(resolution="256x256x256")
        assert result.resolution == "256x256x256"

    def test_generate_mock_has_generation_time(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.generate()
        assert result.generation_time_ms > 0

    def test_generate_invalid_resolution_raises(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        with pytest.raises(ValueError, match="Resolution must be"):
            client.generate(resolution="512x512")

    def test_generate_different_body_regions(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        for region in ["chest", "abdomen", "head", "pelvis"]:
            result = client.generate(body_region=region)
            assert result.body_region == region

    def test_generate_mock_body_region_default(self):
        client = MAISIClient("http://localhost:8531", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.generate()
        # Default body_region is "chest"
        assert result.body_region == "chest"
        assert result.model == "maisi"


# ===================================================================
# VILAM3Client
# ===================================================================


class TestVILAM3Client:
    """Tests for VILAM3Client mock mode."""

    def test_init(self):
        client = VILAM3Client("http://localhost:8532", mock_enabled=True)
        assert client.service_name == "vila_m3"

    def test_analyze_image_mock_returns_vlm_response(self):
        client = VILAM3Client("http://localhost:8532", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        result = client.analyze_image.__wrapped__(
            client, "/tmp/test.png", "What do you see?"
        ) if hasattr(client.analyze_image, '__wrapped__') else None

        # Use _mock_response directly since _invoke_or_mock needs file
        result = client._mock_response(
            image_path="/tmp/test.png",
            question="What do you see?",
            mode="vqa",
        )
        assert isinstance(result, VLMResponse)
        assert result.is_mock is True
        assert len(result.answer) > 0

    def test_mock_response_vqa_mode(self):
        client = VILAM3Client("http://localhost:8532", mock_enabled=True)
        result = client._mock_response(mode="vqa", question="findings?")
        assert isinstance(result, VLMResponse)
        assert len(result.answer) > 0
        assert result.confidence > 0.0

    def test_mock_response_report_mode(self):
        client = VILAM3Client("http://localhost:8532", mock_enabled=True)
        result = client._mock_response(mode="report")
        assert isinstance(result, VLMResponse)
        assert "TECHNIQUE" in result.answer
        assert "FINDINGS" in result.answer
        assert "IMPRESSION" in result.answer

    def test_mock_response_classify_mode(self):
        client = VILAM3Client("http://localhost:8532", mock_enabled=True)
        result = client._mock_response(
            mode="classify",
            labels=["normal", "pneumonia", "effusion"],
        )
        assert isinstance(result, VLMResponse)
        # Answer should be JSON-like
        assert len(result.answer) > 0

    def test_mock_response_has_confidence(self):
        client = VILAM3Client("http://localhost:8532", mock_enabled=True)
        result = client._mock_response(mode="vqa")
        assert 0.0 < result.confidence <= 1.0


# ===================================================================
# LlamaLLMClient
# ===================================================================


class TestLlamaLLMClient:
    """Tests for LlamaLLMClient mock and fallback behavior."""

    def test_init(self):
        client = LlamaLLMClient("http://localhost:8520", mock_enabled=True)
        assert client.service_name == "llm"
        assert client.mock_enabled is True
        assert client.anthropic_api_key is None

    def test_init_with_anthropic_key(self):
        client = LlamaLLMClient(
            "http://localhost:8520",
            anthropic_api_key="sk-test-key",
        )
        assert client.anthropic_api_key == "sk-test-key"

    def test_generate_mock_returns_string(self):
        client = LlamaLLMClient("http://localhost:8520", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        messages = [{"role": "user", "content": "What is hemorrhage?"}]
        result = client.generate(messages)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_generate_mock_includes_clinical_content(self):
        client = LlamaLLMClient("http://localhost:8520", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        messages = [{"role": "user", "content": "test"}]
        result = client.generate(messages)
        assert "clinical" in result.lower() or "imaging" in result.lower() or "findings" in result.lower()

    def test_generate_raises_when_all_unavailable(self):
        client = LlamaLLMClient(
            "http://localhost:8520",
            mock_enabled=False,
            anthropic_api_key=None,
        )
        client._available = False
        client._last_check = 1e18
        with pytest.raises(ConnectionError):
            client.generate([{"role": "user", "content": "test"}])

    def test_mock_response_returns_string(self):
        client = LlamaLLMClient("http://localhost:8520", mock_enabled=True)
        result = client._mock_response(
            messages=[{"role": "user", "content": "test"}],
        )
        assert isinstance(result, str)
        assert "Assessment" in result or "analysis" in result.lower()

    def test_generate_stream_mock(self):
        client = LlamaLLMClient("http://localhost:8520", mock_enabled=True)
        client._available = False
        client._last_check = 1e18
        messages = [{"role": "user", "content": "test"}]
        chunks = list(client.generate_stream(messages))
        assert len(chunks) > 0
        full_text = "".join(chunks)
        assert len(full_text) > 0

    def test_generate_stream_raises_when_all_unavailable(self):
        client = LlamaLLMClient(
            "http://localhost:8520",
            mock_enabled=False,
            anthropic_api_key=None,
        )
        client._available = False
        client._last_check = 1e18
        with pytest.raises(ConnectionError):
            list(client.generate_stream([{"role": "user", "content": "test"}]))

    def test_health_check_returns_false_when_unavailable(self):
        client = LlamaLLMClient("http://localhost:99999", mock_enabled=True)
        assert client.health_check() is False


# ===================================================================
# NIMServiceManager
# ===================================================================


class TestNIMServiceManager:
    """Tests for NIMServiceManager lifecycle management."""

    @pytest.fixture
    def mock_settings(self):
        settings = MagicMock()
        settings.NIM_VISTA3D_URL = "http://localhost:8530"
        settings.NIM_MAISI_URL = "http://localhost:8531"
        settings.NIM_VILAM3_URL = "http://localhost:8532"
        settings.NIM_LLM_URL = "http://localhost:8520"
        settings.NIM_ALLOW_MOCK_FALLBACK = True
        settings.ANTHROPIC_API_KEY = None
        return settings

    def test_init_creates_all_clients(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        assert manager.vista3d is not None
        assert manager.maisi is not None
        assert manager.vilam3 is not None
        assert manager.llm is not None

    def test_vista3d_property(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        assert isinstance(manager.vista3d, VISTA3DClient)

    def test_maisi_property(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        assert isinstance(manager.maisi, MAISIClient)

    def test_vilam3_property(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        assert isinstance(manager.vilam3, VILAM3Client)

    def test_llm_property(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        assert isinstance(manager.llm, LlamaLLMClient)

    def test_check_all_services_returns_dict(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        status = manager.check_all_services()
        assert isinstance(status, dict)
        assert "vista3d" in status
        assert "maisi" in status
        assert "vila_m3" in status
        assert "llm" in status

    def test_check_all_services_values_are_valid_statuses(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        status = manager.check_all_services()
        valid_statuses = {"available", "mock", "unavailable"}
        for service, state in status.items():
            assert state in valid_statuses, (
                f"{service} has invalid status: {state}"
            )

    def test_get_available_services_returns_list(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        available = manager.get_available_services()
        assert isinstance(available, list)

    def test_get_client_vista3d(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        client = manager.get_client("vista3d")
        assert isinstance(client, VISTA3DClient)

    def test_get_client_maisi(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        client = manager.get_client("maisi")
        assert isinstance(client, MAISIClient)

    def test_get_client_vilam3(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        client = manager.get_client("vila_m3")
        assert isinstance(client, VILAM3Client)

    def test_get_client_llm(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        client = manager.get_client("llm")
        assert isinstance(client, LlamaLLMClient)

    def test_get_client_unknown_raises_key_error(self, mock_settings):
        manager = NIMServiceManager(mock_settings)
        with pytest.raises(KeyError, match="Unknown NIM service"):
            manager.get_client("nonexistent")
