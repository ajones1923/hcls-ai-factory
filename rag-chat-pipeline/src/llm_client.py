"""
LLM Client - Interface for language model interactions.

Supports multiple providers:
- Anthropic (Claude)
- OpenAI (GPT)
- vLLM (local models)
"""
from typing import List, Dict, Optional, Generator
from abc import ABC, abstractmethod
import os
from loguru import logger


class BaseLLMClient(ABC):
    """Abstract base class for LLM clients."""

    @abstractmethod
    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        """Generate a response from the LLM."""
        pass

    @abstractmethod
    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        """Generate a streaming response from the LLM."""
        pass


class AnthropicClient(BaseLLMClient):
    """Client for Anthropic Claude models."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "claude-sonnet-4-20250514",
    ):
        try:
            import anthropic
        except ImportError:
            raise ImportError("anthropic package required. Install with: pip install anthropic")

        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not set")

        self.model = model
        self.client = anthropic.Anthropic(api_key=self.api_key)
        logger.info(f"Initialized Anthropic client with model: {model}")

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        message = self.client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_prompt or "",
            messages=[
                {"role": "user", "content": prompt}
            ]
        )
        return message.content[0].text

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        with self.client.messages.stream(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_prompt or "",
            messages=[
                {"role": "user", "content": prompt}
            ]
        ) as stream:
            for text in stream.text_stream:
                yield text


class OpenAIClient(BaseLLMClient):
    """Client for OpenAI GPT models."""

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "gpt-4-turbo-preview",
    ):
        try:
            import openai
        except ImportError:
            raise ImportError("openai package required. Install with: pip install openai")

        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError("OPENAI_API_KEY not set")

        self.model = model
        self.client = openai.OpenAI(api_key=self.api_key)
        logger.info(f"Initialized OpenAI client with model: {model}")

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        response = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
        )
        return response.choices[0].message.content

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        stream = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
            stream=True,
        )

        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


class OllamaClient(BaseLLMClient):
    """Client for Ollama server (OpenAI-compatible API)."""

    def __init__(
        self,
        host: str = None,
        model: str = "llama3.1:70b",
    ):
        try:
            import openai
        except ImportError:
            raise ImportError("openai package required. Install with: pip install openai")

        host = host or os.getenv("OLLAMA_HOST", "http://localhost:11434")
        # Ollama uses /v1 endpoint for OpenAI compatibility
        self.base_url = f"{host}/v1"
        self.model = model
        self.client = openai.OpenAI(
            base_url=self.base_url,
            api_key="ollama",  # Ollama doesn't require API key but needs non-empty value
        )
        logger.info(f"Initialized Ollama client at {host} with model: {model}")

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        response = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
        )
        return response.choices[0].message.content

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        stream = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
            stream=True,
        )

        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


class VLLMClient(BaseLLMClient):
    """Client for local vLLM server (OpenAI-compatible API)."""

    def __init__(
        self,
        host: str = None,
        port: int = None,
        model: str = "meta-llama/Llama-3.1-8B-Instruct",
    ):
        try:
            import openai
        except ImportError:
            raise ImportError("openai package required. Install with: pip install openai")

        host = host or os.getenv("VLLM_HOST", "localhost")
        port = port or int(os.getenv("VLLM_PORT", "8080"))
        self.base_url = f"http://{host}:{port}/v1"
        self.model = model
        self.client = openai.OpenAI(
            base_url=self.base_url,
            api_key="not-needed",  # vLLM doesn't require API key
        )
        logger.info(f"Initialized vLLM client at {self.base_url} with model: {model}")

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        response = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
        )
        return response.choices[0].message.content

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        messages = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        stream = self.client.chat.completions.create(
            model=self.model,
            messages=messages,
            max_tokens=max_tokens,
            temperature=temperature,
            stream=True,
        )

        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


class LLMClient:
    """
    Factory class to create appropriate LLM client based on provider.
    """

    @staticmethod
    def create(
        provider: str = "anthropic",
        model: Optional[str] = None,
        api_key: Optional[str] = None,
        **kwargs
    ) -> BaseLLMClient:
        """
        Create an LLM client.

        Args:
            provider: One of "anthropic", "openai", or "vllm"
            model: Model name (provider-specific)
            api_key: API key (not needed for vllm)
            **kwargs: Additional provider-specific arguments

        Returns:
            LLM client instance
        """
        provider = provider.lower()

        if provider == "anthropic":
            return AnthropicClient(
                api_key=api_key,
                model=model or "claude-sonnet-4-20250514",
            )
        elif provider == "openai":
            return OpenAIClient(
                api_key=api_key,
                model=model or "gpt-4-turbo-preview",
            )
        elif provider == "ollama":
            return OllamaClient(
                host=kwargs.get("host"),
                model=model or os.getenv("LLM_MODEL", "llama3.1:70b"),
            )
        elif provider == "vllm":
            return VLLMClient(
                host=kwargs.get("host"),
                port=kwargs.get("port"),
                model=model or "meta-llama/Llama-3.1-8B-Instruct",
            )
        else:
            raise ValueError(f"Unknown LLM provider: {provider}. Supported: anthropic, openai, ollama, vllm")
