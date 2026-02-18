"""
Tests for LLM Client module.
"""
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.llm_client import (
    AnthropicClient,
    BaseLLMClient,
    LLMClient,
    OllamaClient,
    OpenAIClient,
    VLLMClient,
)


def _mock_anthropic():
    """Create a mock anthropic module and inject it into sys.modules."""
    return patch.dict("sys.modules", {"anthropic": MagicMock()})


def _mock_openai():
    """Create a mock openai module and inject it into sys.modules."""
    return patch.dict("sys.modules", {"openai": MagicMock()})


class TestLLMClientFactory:
    """Tests for LLMClient factory."""

    def test_create_anthropic_client(self):
        """Test creating Anthropic client via factory."""
        with patch.dict(os.environ, {"ANTHROPIC_API_KEY": "test-key"}):  # noqa: SIM117
            with _mock_anthropic():
                client = LLMClient.create(provider="anthropic")
                assert isinstance(client, AnthropicClient)

    def test_create_openai_client(self):
        """Test creating OpenAI client via factory."""
        with patch.dict(os.environ, {"OPENAI_API_KEY": "test-key"}):  # noqa: SIM117
            with _mock_openai():
                client = LLMClient.create(provider="openai")
                assert isinstance(client, OpenAIClient)

    def test_create_ollama_client(self):
        """Test creating Ollama client via factory."""
        with _mock_openai():
            client = LLMClient.create(provider="ollama")
            assert isinstance(client, OllamaClient)

    def test_create_vllm_client(self):
        """Test creating vLLM client via factory."""
        with _mock_openai():
            client = LLMClient.create(provider="vllm")
            assert isinstance(client, VLLMClient)

    def test_create_unknown_provider(self):
        """Test error for unknown provider."""
        with pytest.raises(ValueError, match="Unknown LLM provider"):
            LLMClient.create(provider="unknown")

    def test_create_with_custom_model(self):
        """Test creating client with custom model."""
        with patch.dict(os.environ, {"ANTHROPIC_API_KEY": "test-key"}):  # noqa: SIM117
            with _mock_anthropic():
                client = LLMClient.create(
                    provider="anthropic",
                    model="claude-3-opus-20240229"
                )
                assert client.model == "claude-3-opus-20240229"


class TestAnthropicClient:
    """Tests for Anthropic client."""

    def test_init_with_api_key(self):
        """Test initialization with API key."""
        mock_anthropic = MagicMock()
        with patch.dict("sys.modules", {"anthropic": mock_anthropic}):
            client = AnthropicClient(api_key="test-key")
            assert client.api_key == "test-key"
            mock_anthropic.Anthropic.assert_called_with(api_key="test-key")

    def test_init_from_env(self):
        """Test initialization from environment variable."""
        with patch.dict(os.environ, {"ANTHROPIC_API_KEY": "env-key"}):  # noqa: SIM117
            with _mock_anthropic():
                client = AnthropicClient()
                assert client.api_key == "env-key"

    def test_init_no_key_raises(self):
        """Test error when no API key available."""
        with patch.dict(os.environ, {}, clear=True):
            # Remove ANTHROPIC_API_KEY if it exists
            os.environ.pop("ANTHROPIC_API_KEY", None)
            with _mock_anthropic():  # noqa: SIM117
                with pytest.raises(ValueError, match="ANTHROPIC_API_KEY"):
                    AnthropicClient()

    def test_generate(self):
        """Test generate method."""
        mock_anthropic = MagicMock()
        mock_message = Mock()
        mock_message.content = [Mock(text="Generated response")]
        mock_anthropic.Anthropic.return_value.messages.create.return_value = mock_message

        with patch.dict("sys.modules", {"anthropic": mock_anthropic}):
            client = AnthropicClient(api_key="test-key")
            response = client.generate("Test prompt", system_prompt="System")

            assert response == "Generated response"
            mock_anthropic.Anthropic.return_value.messages.create.assert_called_once()

    def test_generate_stream(self):
        """Test streaming generation."""
        mock_anthropic = MagicMock()
        # Mock the stream context manager
        mock_stream = Mock()
        mock_stream.__enter__ = Mock(return_value=mock_stream)
        mock_stream.__exit__ = Mock(return_value=False)
        mock_stream.text_stream = iter(["Hello", " ", "world"])
        mock_anthropic.Anthropic.return_value.messages.stream.return_value = mock_stream

        with patch.dict("sys.modules", {"anthropic": mock_anthropic}):
            client = AnthropicClient(api_key="test-key")
            tokens = list(client.generate_stream("Test prompt"))

            assert tokens == ["Hello", " ", "world"]


class TestOllamaClient:
    """Tests for Ollama client."""

    def test_init_default_host(self):
        """Test default host configuration."""
        with _mock_openai():
            client = OllamaClient()
            assert client.base_url == "http://localhost:11434/v1"

    def test_init_custom_host(self):
        """Test custom host configuration."""
        with _mock_openai():
            client = OllamaClient(host="http://gpu-server:11434")
            assert client.base_url == "http://gpu-server:11434/v1"

    def test_generate(self):
        """Test generate method."""
        mock_openai = MagicMock()
        mock_response = Mock()
        mock_response.choices = [Mock(message=Mock(content="Ollama response"))]
        mock_openai.OpenAI.return_value.chat.completions.create.return_value = mock_response

        with patch.dict("sys.modules", {"openai": mock_openai}):
            client = OllamaClient()
            response = client.generate("Test prompt")

            assert response == "Ollama response"


class TestVLLMClient:
    """Tests for vLLM client."""

    def test_init_default_config(self):
        """Test default configuration."""
        with _mock_openai():
            client = VLLMClient()
            assert client.base_url == "http://localhost:8080/v1"
            assert client.model == "meta-llama/Llama-3.1-8B-Instruct"

    def test_init_custom_config(self):
        """Test custom configuration."""
        with _mock_openai():
            client = VLLMClient(
                host="gpu-server",
                port=8000,
                model="meta-llama/Llama-3.1-70B-Instruct"
            )
            assert client.base_url == "http://gpu-server:8000/v1"
            assert client.model == "meta-llama/Llama-3.1-70B-Instruct"


class TestBaseLLMClient:
    """Tests for abstract base class."""

    def test_abstract_methods(self):
        """Test that BaseLLMClient cannot be instantiated."""
        with pytest.raises(TypeError):
            BaseLLMClient()

    def test_subclass_must_implement(self):
        """Test subclass must implement abstract methods."""
        class IncompleteClient(BaseLLMClient):
            pass

        with pytest.raises(TypeError):
            IncompleteClient()


class TestLLMClientIntegration:
    """Integration-style tests for LLM clients."""

    def test_client_with_temperature(self):
        """Test temperature parameter is passed."""
        mock_anthropic = MagicMock()
        mock_message = Mock()
        mock_message.content = [Mock(text="Response")]
        mock_anthropic.Anthropic.return_value.messages.create.return_value = mock_message

        with patch.dict("sys.modules", {"anthropic": mock_anthropic}):
            client = AnthropicClient(api_key="test-key")
            client.generate("Test", temperature=0.5)

            call_kwargs = mock_anthropic.Anthropic.return_value.messages.create.call_args[1]
            assert call_kwargs["temperature"] == 0.5

    def test_client_with_max_tokens(self):
        """Test max_tokens parameter is passed."""
        mock_anthropic = MagicMock()
        mock_message = Mock()
        mock_message.content = [Mock(text="Response")]
        mock_anthropic.Anthropic.return_value.messages.create.return_value = mock_message

        with patch.dict("sys.modules", {"anthropic": mock_anthropic}):
            client = AnthropicClient(api_key="test-key")
            client.generate("Test", max_tokens=2048)

            call_kwargs = mock_anthropic.Anthropic.return_value.messages.create.call_args[1]
            assert call_kwargs["max_tokens"] == 2048
