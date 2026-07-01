"""
Unified LLM client for the HCLS AI Factory.

Consolidates and extends:
  - core/engines/precision-intelligence/src/llm_client.py  (BaseLLMClient ABC, AnthropicClient,
    OpenAIClient, OllamaClient, VLLMClient, LLMClient factory)

Enhancements over the original:
  - Exponential backoff (2s-32s) on 429 / 529 rate-limit responses
  - Daily token budget tracker
  - Response cache with TTL + LRU eviction
  - JSON extraction with repair (tolerates markdown fences, trailing commas)
  - Prometheus metrics (api_calls by provider/status, latency, input/output tokens)
  - Singleton ``get_llm_client()``
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import threading
import time
from abc import ABC, abstractmethod
from collections import OrderedDict
from datetime import date
from functools import lru_cache
from typing import Any, Dict, Generator, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Prometheus metrics (optional)
# ---------------------------------------------------------------------------

try:
    from prometheus_client import Counter, Histogram

    _api_calls = Counter(
        "hcls_llm_api_calls_total",
        "LLM API calls",
        ["provider", "status"],  # status: success | error | rate_limited
    )
    _latency = Histogram(
        "hcls_llm_latency_seconds",
        "LLM generation latency",
        ["provider"],
        buckets=(0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, 60.0),
    )
    _input_tokens = Counter(
        "hcls_llm_input_tokens_total",
        "Total input tokens",
        ["provider"],
    )
    _output_tokens = Counter(
        "hcls_llm_output_tokens_total",
        "Total output tokens",
        ["provider"],
    )
    _HAS_PROMETHEUS = True
except ImportError:
    _HAS_PROMETHEUS = False


# ---------------------------------------------------------------------------
# Daily token budget tracker
# ---------------------------------------------------------------------------

class _TokenBudget:
    """Thread-safe daily token budget tracker."""

    def __init__(self, daily_limit: int = 0) -> None:
        self._daily_limit = daily_limit  # 0 = unlimited
        self._lock = threading.Lock()
        self._date: str = ""
        self._used: int = 0

    def configure(self, daily_limit: int) -> None:
        self._daily_limit = daily_limit

    def record(self, tokens: int) -> None:
        with self._lock:
            today = date.today().isoformat()
            if today != self._date:
                self._date = today
                self._used = 0
            self._used += tokens

    def check(self) -> bool:
        """Return True if budget allows more tokens."""
        if self._daily_limit <= 0:
            return True
        with self._lock:
            today = date.today().isoformat()
            if today != self._date:
                return True
            return self._used < self._daily_limit

    @property
    def remaining(self) -> int:
        if self._daily_limit <= 0:
            return -1  # unlimited
        with self._lock:
            today = date.today().isoformat()
            if today != self._date:
                return self._daily_limit
            return max(0, self._daily_limit - self._used)


_budget = _TokenBudget()


# ---------------------------------------------------------------------------
# Response cache (TTL + LRU)
# ---------------------------------------------------------------------------

class _ResponseCache:
    """Thread-safe TTL + LRU cache for LLM responses."""

    def __init__(self, max_size: int = 256, ttl: float = 300.0) -> None:
        self._max_size = max_size
        self._ttl = ttl
        self._lock = threading.Lock()
        self._store: OrderedDict[str, Tuple[float, str]] = OrderedDict()

    def configure(self, max_size: int, ttl: float) -> None:
        self._max_size = max_size
        self._ttl = ttl

    @staticmethod
    def _key(prompt: str, system: str, model: str, temperature: float) -> str:
        raw = f"{model}|{temperature}|{system}|{prompt}"
        return hashlib.sha256(raw.encode("utf-8")).hexdigest()

    def get(
        self,
        prompt: str,
        system: str,
        model: str,
        temperature: float,
    ) -> Optional[str]:
        if self._ttl <= 0:
            return None
        key = self._key(prompt, system, model, temperature)
        with self._lock:
            entry = self._store.get(key)
            if entry is None:
                return None
            ts, value = entry
            if time.monotonic() - ts > self._ttl:
                del self._store[key]
                return None
            self._store.move_to_end(key)
            return value

    def put(
        self,
        prompt: str,
        system: str,
        model: str,
        temperature: float,
        response: str,
    ) -> None:
        if self._ttl <= 0:
            return
        key = self._key(prompt, system, model, temperature)
        with self._lock:
            if key in self._store:
                self._store.move_to_end(key)
            self._store[key] = (time.monotonic(), response)
            while len(self._store) > self._max_size:
                self._store.popitem(last=False)

    def clear(self) -> None:
        with self._lock:
            self._store.clear()


_resp_cache = _ResponseCache()


# ---------------------------------------------------------------------------
# JSON extraction with repair
# ---------------------------------------------------------------------------

def extract_json(text: str) -> Any:
    """
    Extract JSON from LLM output.  Handles:
      - Bare JSON objects / arrays
      - Markdown fenced code blocks (```json ... ```)
      - Trailing commas before closing brackets
      - Single quotes -> double quotes
    Returns the parsed Python object or raises ``ValueError``.
    """
    # Strip markdown fences
    fenced = re.search(r"```(?:json)?\s*([\s\S]*?)```", text)
    if fenced:
        text = fenced.group(1)

    text = text.strip()

    # Attempt direct parse first
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # Find first { or [ and last } or ]
    start = None
    end = None
    for i, c in enumerate(text):
        if c in "{[":
            start = i
            break
    if start is not None:
        closer = "}" if text[start] == "{" else "]"
        for i in range(len(text) - 1, start - 1, -1):
            if text[i] == closer:
                end = i + 1
                break

    if start is not None and end is not None:
        candidate = text[start:end]
        # Repair trailing commas:  ,\s*} -> }  and  ,\s*] -> ]
        candidate = re.sub(r",\s*}", "}", candidate)
        candidate = re.sub(r",\s*\]", "]", candidate)

        try:
            return json.loads(candidate)
        except json.JSONDecodeError:
            pass

        # Try replacing single quotes with double quotes
        try:
            return json.loads(candidate.replace("'", '"'))
        except json.JSONDecodeError:
            pass

    raise ValueError(f"Could not extract valid JSON from LLM response: {text[:200]}")


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------

class BaseLLMClient(ABC):
    """Abstract base class for LLM provider clients."""

    provider: str = "unknown"

    @abstractmethod
    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        """Generate a single response."""
        ...

    @abstractmethod
    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        """Generate a streaming response (yields text chunks)."""
        ...

    def generate_json(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.3,
    ) -> Any:
        """Generate and parse JSON (with repair)."""
        raw = self.generate(
            prompt=prompt,
            system_prompt=system_prompt,
            max_tokens=max_tokens,
            temperature=temperature,
        )
        return extract_json(raw)


# ---------------------------------------------------------------------------
# Retry / backoff mixin
# ---------------------------------------------------------------------------

class _RetryMixin:
    """
    Adds exponential backoff on 429 (Too Many Requests) and 529
    (Overloaded) responses.
    """

    _retry_base: float = 2.0
    _retry_max: float = 32.0
    _max_retries: int = 6

    def _should_retry(self, exc: Exception) -> bool:
        """Check if the exception indicates a retryable rate-limit error."""
        status = getattr(exc, "status_code", None) or getattr(
            getattr(exc, "response", None), "status_code", None
        )
        return status in (429, 529)

    def _generate_with_retry(
        self,
        generate_fn: Any,
        prompt: str,
        system_prompt: Optional[str],
        max_tokens: int,
        temperature: float,
    ) -> str:
        provider = getattr(self, "provider", "unknown")
        last_exc: Optional[Exception] = None

        for attempt in range(self._max_retries):
            # Budget check
            if not _budget.check():
                raise RuntimeError(
                    f"Daily token budget exhausted ({_budget.remaining} remaining)"
                )

            # Cache check
            cached = _resp_cache.get(
                prompt,
                system_prompt or "",
                getattr(self, "model", ""),
                temperature,
            )
            if cached is not None:
                return cached

            t0 = time.monotonic()
            try:
                result = generate_fn(prompt, system_prompt, max_tokens, temperature)
                elapsed = time.monotonic() - t0

                if _HAS_PROMETHEUS:
                    _api_calls.labels(provider=provider, status="success").inc()
                    _latency.labels(provider=provider).observe(elapsed)

                # Cache the response
                _resp_cache.put(
                    prompt,
                    system_prompt or "",
                    getattr(self, "model", ""),
                    temperature,
                    result,
                )

                return result

            except Exception as exc:
                last_exc = exc
                if self._should_retry(exc) and attempt < self._max_retries - 1:
                    delay = min(
                        self._retry_base * (2 ** attempt), self._retry_max
                    )
                    if _HAS_PROMETHEUS:
                        _api_calls.labels(
                            provider=provider, status="rate_limited"
                        ).inc()
                    logger.warning(
                        "LLM %s rate-limited (attempt %d/%d), backing off %.1fs: %s",
                        provider,
                        attempt + 1,
                        self._max_retries,
                        delay,
                        exc,
                    )
                    time.sleep(delay)
                else:
                    if _HAS_PROMETHEUS:
                        _api_calls.labels(provider=provider, status="error").inc()
                    raise

        raise RuntimeError(
            f"LLM {provider} failed after {self._max_retries} retries: {last_exc}"
        )


# ---------------------------------------------------------------------------
# AnthropicClient
# ---------------------------------------------------------------------------

class AnthropicClient(_RetryMixin, BaseLLMClient):
    """Client for Anthropic Claude models."""

    provider = "anthropic"

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "claude-sonnet-4-20250514",
    ) -> None:
        try:
            import anthropic
        except ImportError as exc:
            raise ImportError(
                "anthropic package required. Install with: pip install anthropic"
            ) from exc

        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not set")

        self.model = model
        self._client = anthropic.Anthropic(api_key=self.api_key)
        logger.info("Initialized Anthropic client: model=%s", model)

    def _raw_generate(
        self,
        prompt: str,
        system_prompt: Optional[str],
        max_tokens: int,
        temperature: float,
    ) -> str:
        message = self._client.messages.create(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_prompt or "",
            messages=[{"role": "user", "content": prompt}],
        )
        # Record token usage
        input_tok = getattr(message.usage, "input_tokens", 0)
        output_tok = getattr(message.usage, "output_tokens", 0)
        _budget.record(input_tok + output_tok)
        if _HAS_PROMETHEUS:
            _input_tokens.labels(provider=self.provider).inc(input_tok)
            _output_tokens.labels(provider=self.provider).inc(output_tok)
        return message.content[0].text

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        return self._generate_with_retry(
            self._raw_generate, prompt, system_prompt, max_tokens, temperature
        )

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        with self._client.messages.stream(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            system=system_prompt or "",
            messages=[{"role": "user", "content": prompt}],
        ) as stream:
            for text in stream.text_stream:
                yield text


# ---------------------------------------------------------------------------
# OpenAIClient
# ---------------------------------------------------------------------------

class OpenAIClient(_RetryMixin, BaseLLMClient):
    """Client for OpenAI GPT models."""

    provider = "openai"

    def __init__(
        self,
        api_key: Optional[str] = None,
        model: str = "gpt-4-turbo-preview",
    ) -> None:
        try:
            import openai
        except ImportError as exc:
            raise ImportError(
                "openai package required. Install with: pip install openai"
            ) from exc

        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        if not self.api_key:
            raise ValueError("OPENAI_API_KEY not set")

        self.model = model
        self._client = openai.OpenAI(api_key=self.api_key)
        logger.info("Initialized OpenAI client: model=%s", model)

    def _build_messages(
        self, prompt: str, system_prompt: Optional[str]
    ) -> List[Dict[str, str]]:
        msgs: List[Dict[str, str]] = []
        if system_prompt:
            msgs.append({"role": "system", "content": system_prompt})
        msgs.append({"role": "user", "content": prompt})
        return msgs

    def _raw_generate(
        self,
        prompt: str,
        system_prompt: Optional[str],
        max_tokens: int,
        temperature: float,
    ) -> str:
        response = self._client.chat.completions.create(
            model=self.model,
            messages=self._build_messages(prompt, system_prompt),
            max_tokens=max_tokens,
            temperature=temperature,
        )
        usage = response.usage
        if usage:
            _budget.record((usage.prompt_tokens or 0) + (usage.completion_tokens or 0))
            if _HAS_PROMETHEUS:
                _input_tokens.labels(provider=self.provider).inc(
                    usage.prompt_tokens or 0
                )
                _output_tokens.labels(provider=self.provider).inc(
                    usage.completion_tokens or 0
                )
        return response.choices[0].message.content

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        return self._generate_with_retry(
            self._raw_generate, prompt, system_prompt, max_tokens, temperature
        )

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        stream = self._client.chat.completions.create(
            model=self.model,
            messages=self._build_messages(prompt, system_prompt),
            max_tokens=max_tokens,
            temperature=temperature,
            stream=True,
        )
        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


# ---------------------------------------------------------------------------
# OllamaClient
# ---------------------------------------------------------------------------

class OllamaClient(_RetryMixin, BaseLLMClient):
    """Client for Ollama (OpenAI-compatible API at /v1)."""

    provider = "ollama"

    def __init__(
        self,
        host: Optional[str] = None,
        model: str = "llama3.1:70b",
    ) -> None:
        try:
            import openai
        except ImportError as exc:
            raise ImportError(
                "openai package required. Install with: pip install openai"
            ) from exc

        host = host or os.getenv("OLLAMA_HOST", "http://localhost:11434")
        self.base_url = f"{host}/v1"
        self.model = model
        self._client = openai.OpenAI(
            base_url=self.base_url,
            api_key="ollama",
        )
        logger.info("Initialized Ollama client: host=%s  model=%s", host, model)

    def _build_messages(
        self, prompt: str, system_prompt: Optional[str]
    ) -> List[Dict[str, str]]:
        msgs: List[Dict[str, str]] = []
        if system_prompt:
            msgs.append({"role": "system", "content": system_prompt})
        msgs.append({"role": "user", "content": prompt})
        return msgs

    def _raw_generate(
        self,
        prompt: str,
        system_prompt: Optional[str],
        max_tokens: int,
        temperature: float,
    ) -> str:
        response = self._client.chat.completions.create(
            model=self.model,
            messages=self._build_messages(prompt, system_prompt),
            max_tokens=max_tokens,
            temperature=temperature,
        )
        return response.choices[0].message.content

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        return self._generate_with_retry(
            self._raw_generate, prompt, system_prompt, max_tokens, temperature
        )

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        stream = self._client.chat.completions.create(
            model=self.model,
            messages=self._build_messages(prompt, system_prompt),
            max_tokens=max_tokens,
            temperature=temperature,
            stream=True,
        )
        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


# ---------------------------------------------------------------------------
# VLLMClient
# ---------------------------------------------------------------------------

class VLLMClient(_RetryMixin, BaseLLMClient):
    """Client for local vLLM server (OpenAI-compatible API)."""

    provider = "vllm"

    def __init__(
        self,
        host: Optional[str] = None,
        port: Optional[int] = None,
        model: str = "meta-llama/Llama-3.1-8B-Instruct",
    ) -> None:
        try:
            import openai
        except ImportError as exc:
            raise ImportError(
                "openai package required. Install with: pip install openai"
            ) from exc

        host = host or os.getenv("VLLM_HOST", "localhost")
        port = port or int(os.getenv("VLLM_PORT", "8080"))
        self.base_url = f"http://{host}:{port}/v1"
        self.model = model
        self._client = openai.OpenAI(
            base_url=self.base_url,
            api_key="not-needed",
        )
        logger.info(
            "Initialized vLLM client: %s  model=%s", self.base_url, model
        )

    def _build_messages(
        self, prompt: str, system_prompt: Optional[str]
    ) -> List[Dict[str, str]]:
        msgs: List[Dict[str, str]] = []
        if system_prompt:
            msgs.append({"role": "system", "content": system_prompt})
        msgs.append({"role": "user", "content": prompt})
        return msgs

    def _raw_generate(
        self,
        prompt: str,
        system_prompt: Optional[str],
        max_tokens: int,
        temperature: float,
    ) -> str:
        response = self._client.chat.completions.create(
            model=self.model,
            messages=self._build_messages(prompt, system_prompt),
            max_tokens=max_tokens,
            temperature=temperature,
        )
        return response.choices[0].message.content

    def generate(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> str:
        return self._generate_with_retry(
            self._raw_generate, prompt, system_prompt, max_tokens, temperature
        )

    def generate_stream(
        self,
        prompt: str,
        system_prompt: Optional[str] = None,
        max_tokens: int = 1024,
        temperature: float = 0.7,
    ) -> Generator[str, None, None]:
        stream = self._client.chat.completions.create(
            model=self.model,
            messages=self._build_messages(prompt, system_prompt),
            max_tokens=max_tokens,
            temperature=temperature,
            stream=True,
        )
        for chunk in stream:
            if chunk.choices[0].delta.content:
                yield chunk.choices[0].delta.content


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

class LLMClientFactory:
    """
    Factory to create the appropriate LLM client from provider name.
    """

    @staticmethod
    def create(
        provider: str = "anthropic",
        model: Optional[str] = None,
        api_key: Optional[str] = None,
        **kwargs: Any,
    ) -> BaseLLMClient:
        """
        Create an LLM client.

        Parameters
        ----------
        provider : str
            One of ``anthropic``, ``openai``, ``ollama``, ``vllm``.
        model : str or None
            Provider-specific model identifier.
        api_key : str or None
            API key (not needed for Ollama / vLLM).
        **kwargs
            Forwarded to the provider constructor (host, port, etc.).
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
            raise ValueError(
                f"Unknown LLM provider: {provider!r}. "
                "Supported: anthropic, openai, ollama, vllm"
            )


# ---------------------------------------------------------------------------
# Singleton accessor
# ---------------------------------------------------------------------------

_llm_instance: Optional[BaseLLMClient] = None
_llm_lock = threading.Lock()


def get_llm_client(settings: Any = None) -> BaseLLMClient:
    """
    Return a singleton LLM client based on ``HCLSSettings``.

    Pass ``settings=None`` to auto-load via ``get_settings()``.
    """
    global _llm_instance
    if _llm_instance is not None:
        return _llm_instance

    with _llm_lock:
        if _llm_instance is not None:
            return _llm_instance

        if settings is None:
            from hcls_common.config import get_settings

            settings = get_settings()

        provider = getattr(settings, "llm_provider", "ollama")
        model = getattr(settings, "llm_model", None)
        api_key = getattr(settings, "llm_api_key", None)
        daily_budget = getattr(settings, "llm_daily_token_budget", 0)
        cache_ttl = getattr(settings, "llm_response_cache_ttl", 300)
        cache_max = getattr(settings, "llm_response_cache_max", 256)
        retry_base = getattr(settings, "llm_retry_base", 2.0)
        retry_max = getattr(settings, "llm_retry_max", 32.0)

        # Configure globals
        _budget.configure(daily_budget)
        _resp_cache.configure(cache_max, cache_ttl)

        client = LLMClientFactory.create(
            provider=provider,
            model=model,
            api_key=api_key,
            host=getattr(settings, "ollama_host", None)
            if provider == "ollama"
            else getattr(settings, "vllm_host", None),
            port=getattr(settings, "vllm_port", None)
            if provider == "vllm"
            else None,
        )

        # Set retry params
        if isinstance(client, _RetryMixin):
            client._retry_base = retry_base
            client._retry_max = retry_max

        _llm_instance = client
        return _llm_instance


def reset_llm_client() -> None:
    """Clear the singleton (useful in tests)."""
    global _llm_instance
    with _llm_lock:
        _llm_instance = None
    _resp_cache.clear()
