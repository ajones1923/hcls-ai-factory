"""LLM NIM client with Llama-3 NIM primary, cloud NIM, and Anthropic Claude fallback.

Provides text generation for clinical reasoning, report synthesis,
and multi-turn conversation with automatic fallback chain:
  local NIM Llama-3 -> NVIDIA cloud NIM -> Anthropic Claude -> mock
"""

from typing import Any, Dict, Generator, List, Optional

from loguru import logger

from .base import BaseNIMClient


class LlamaLLMClient(BaseNIMClient):
    """Client for LLM inference via NIM (Llama-3) with cloud and Anthropic fallback.

    Priority order:
      1. Local NIM Llama-3 endpoint (OpenAI-compatible API)
      2. NVIDIA Cloud NIM (Llama-3.1-8B via integrate.api.nvidia.com)
      3. Anthropic Claude API (if anthropic_api_key provided)
      4. Mock response (if mock_enabled)

    Uses the OpenAI Python client for NIM communication since NVIDIA
    NIM exposes an OpenAI-compatible /v1/chat/completions endpoint.
    """

    def __init__(
        self,
        base_url: str,
        mock_enabled: bool = True,
        anthropic_api_key: Optional[str] = None,
        nvidia_api_key: Optional[str] = None,
        cloud_url: str = "https://integrate.api.nvidia.com/v1",
        cloud_llm_model: str = "meta/llama-3.1-8b-instruct",
        local_llm_model: str = "meta/llama3-70b-instruct",
    ):
        super().__init__(base_url, service_name="llm", mock_enabled=mock_enabled)
        self.anthropic_api_key = anthropic_api_key
        self.nvidia_api_key = nvidia_api_key
        self.cloud_url = cloud_url.rstrip("/")
        self.cloud_llm_model = cloud_llm_model
        self.local_llm_model = local_llm_model
        self._openai_client = None
        self._cloud_client = None
        self._anthropic_client = None
        self._cloud_available: Optional[bool] = None

    def health_check(self) -> bool:
        """Check NIM LLM availability via /v1/models endpoint (OpenAI-compat)."""
        try:
            import requests
            resp = requests.get(f"{self.base_url}/v1/models", timeout=5)
            return resp.status_code == 200
        except Exception:
            return False

    def cloud_health_check(self) -> bool:
        """Check NVIDIA cloud NIM availability."""
        if not self.nvidia_api_key:
            return False
        try:
            import requests
            resp = requests.get(
                f"{self.cloud_url}/models",
                headers={"Authorization": f"Bearer {self.nvidia_api_key}"},
                timeout=10,
            )
            return resp.status_code == 200
        except Exception:
            return False

    def is_cloud_available(self) -> bool:
        """Check if NVIDIA cloud NIM is available (cached)."""
        import time
        now = time.time()
        if self._cloud_available is None or (now - self._last_check) > self._check_interval:
            self._cloud_available = self.cloud_health_check()
            if self._cloud_available:
                logger.info(f"NVIDIA Cloud NIM LLM available at {self.cloud_url}")
            else:
                logger.debug("NVIDIA Cloud NIM LLM not available")
        return self._cloud_available

    def _get_openai_client(self):
        """Lazy-initialize the OpenAI client for local NIM."""
        if self._openai_client is None:
            try:
                from openai import OpenAI
                self._openai_client = OpenAI(
                    base_url=f"{self.base_url}/v1",
                    api_key="not-needed",  # NIM does not require an API key
                )
            except ImportError:
                logger.warning("openai package not installed; NIM LLM unavailable")
                return None
        return self._openai_client

    def _get_cloud_client(self):
        """Lazy-initialize the OpenAI client for NVIDIA cloud NIM."""
        if self._cloud_client is None and self.nvidia_api_key:
            try:
                from openai import OpenAI
                self._cloud_client = OpenAI(
                    base_url=self.cloud_url,
                    api_key=self.nvidia_api_key,
                )
            except ImportError:
                logger.warning("openai package not installed; cloud NIM unavailable")
                return None
        return self._cloud_client

    def _get_anthropic_client(self):
        """Lazy-initialize the Anthropic client."""
        if self._anthropic_client is None and self.anthropic_api_key:
            try:
                import anthropic
                self._anthropic_client = anthropic.Anthropic(
                    api_key=self.anthropic_api_key,
                )
            except ImportError:
                logger.warning(
                    "anthropic package not installed; Claude fallback unavailable"
                )
                return None
        return self._anthropic_client

    def _cloud_generate(
        self,
        messages: List[Dict],
        temperature: float = 0.1,
        max_tokens: int = 2048,
    ) -> Optional[str]:
        """Generate text via NVIDIA cloud NIM endpoint.

        Uses OpenAI-compatible API at integrate.api.nvidia.com.

        Args:
            messages: List of message dicts with "role" and "content" keys.
            temperature: Sampling temperature.
            max_tokens: Maximum response tokens.

        Returns:
            Generated text string, or None if cloud generation fails.
        """
        client = self._get_cloud_client()
        if client is None:
            return None

        try:
            response = client.chat.completions.create(
                model=self.cloud_llm_model,
                messages=messages,
                temperature=temperature,
                max_tokens=max_tokens,
            )
            text = response.choices[0].message.content
            logger.info(
                f"Cloud NIM LLM ({self.cloud_llm_model}) generated {len(text)} chars "
                f"({response.usage.total_tokens} tokens)"
            )
            return text
        except Exception as e:
            logger.error(f"Cloud NIM LLM generation failed: {e}")
            return None

    def _cloud_generate_stream(
        self,
        messages: List[Dict],
        temperature: float = 0.1,
    ) -> Optional[Generator[str, None, None]]:
        """Stream text via NVIDIA cloud NIM endpoint.

        Args:
            messages: List of message dicts.
            temperature: Sampling temperature.

        Returns:
            Generator yielding text chunks, or None if streaming fails.
        """
        client = self._get_cloud_client()
        if client is None:
            return None

        try:
            stream = client.chat.completions.create(
                model=self.cloud_llm_model,
                messages=messages,
                temperature=temperature,
                stream=True,
            )
            return stream
        except Exception as e:
            logger.error(f"Cloud NIM LLM streaming failed: {e}")
            return None

    def generate(
        self,
        messages: List[Dict],
        temperature: float = 0.1,
        max_tokens: int = 2048,
    ) -> str:
        """Generate a text response from the LLM.

        Attempts local NIM Llama-3 first, then NVIDIA cloud NIM,
        then Anthropic Claude, then mock if all else fails.

        Args:
            messages: List of message dicts with "role" and "content" keys.
                Example: [{"role": "user", "content": "Describe this finding."}]
            temperature: Sampling temperature (0.0 = deterministic).
            max_tokens: Maximum response tokens.

        Returns:
            Generated text response string.
        """
        # Attempt 1: Local NIM Llama-3
        if self.is_available():
            try:
                client = self._get_openai_client()
                if client is not None:
                    response = client.chat.completions.create(
                        model=self.local_llm_model,
                        messages=messages,
                        temperature=temperature,
                        max_tokens=max_tokens,
                    )
                    text = response.choices[0].message.content
                    logger.info(
                        f"NIM Llama-3 generated {len(text)} chars "
                        f"({response.usage.total_tokens} tokens)"
                    )
                    return text
            except Exception as e:
                logger.error(f"NIM LLM generation failed: {e}")

        # Attempt 2: NVIDIA Cloud NIM
        if self.nvidia_api_key:
            cloud_result = self._cloud_generate(messages, temperature, max_tokens)
            if cloud_result is not None:
                return cloud_result

        # Attempt 3: Anthropic Claude
        if self.anthropic_api_key:
            try:
                client = self._get_anthropic_client()
                if client is not None:
                    logger.info("Falling back to Anthropic Claude API")

                    # Convert messages: extract system message if present
                    system_msg = ""
                    chat_messages = []
                    for msg in messages:
                        if msg["role"] == "system":
                            system_msg = msg["content"]
                        else:
                            chat_messages.append(msg)

                    # Ensure at least one user message
                    if not chat_messages:
                        chat_messages = [{"role": "user", "content": "Hello."}]

                    kwargs: Dict[str, Any] = {
                        "model": "claude-sonnet-4-20250514",
                        "max_tokens": max_tokens,
                        "messages": chat_messages,
                    }
                    if system_msg:
                        kwargs["system"] = system_msg

                    response = client.messages.create(**kwargs)
                    text = response.content[0].text
                    logger.info(
                        f"Claude generated {len(text)} chars "
                        f"(input={response.usage.input_tokens}, "
                        f"output={response.usage.output_tokens})"
                    )
                    return text
            except Exception as e:
                logger.error(f"Anthropic Claude generation failed: {e}")

        # Attempt 4: Mock
        if self.mock_enabled:
            logger.info("Using mock LLM response (no NIM, cloud, or Claude available)")
            return self._mock_response(messages=messages)

        raise ConnectionError(
            "LLM unavailable: local NIM unreachable, cloud NIM failed, "
            "no Anthropic API key, and mock disabled"
        )

    def generate_stream(
        self,
        messages: List[Dict],
        temperature: float = 0.1,
    ) -> Generator[str, None, None]:
        """Stream a text response from the LLM token by token.

        Attempts local NIM Llama-3 streaming first, then NVIDIA cloud NIM
        streaming, then Anthropic Claude streaming, then yields mock
        response in chunks.

        Args:
            messages: List of message dicts with "role" and "content" keys.
            temperature: Sampling temperature.

        Yields:
            Text chunks as they are generated.
        """
        # Attempt 1: Local NIM Llama-3 streaming
        if self.is_available():
            try:
                client = self._get_openai_client()
                if client is not None:
                    stream = client.chat.completions.create(
                        model=self.local_llm_model,
                        messages=messages,
                        temperature=temperature,
                        stream=True,
                    )
                    for chunk in stream:
                        if chunk.choices and chunk.choices[0].delta.content:
                            yield chunk.choices[0].delta.content
                    return
            except Exception as e:
                logger.error(f"NIM LLM streaming failed: {e}")

        # Attempt 2: NVIDIA Cloud NIM streaming
        if self.nvidia_api_key:
            try:
                stream = self._cloud_generate_stream(messages, temperature)
                if stream is not None:
                    logger.info(f"Streaming via NVIDIA Cloud NIM ({self.cloud_llm_model})")
                    for chunk in stream:
                        if chunk.choices and chunk.choices[0].delta.content:
                            yield chunk.choices[0].delta.content
                    return
            except Exception as e:
                logger.error(f"Cloud NIM LLM streaming failed: {e}")

        # Attempt 3: Anthropic Claude streaming
        if self.anthropic_api_key:
            try:
                client = self._get_anthropic_client()
                if client is not None:
                    logger.info("Streaming via Anthropic Claude API")

                    system_msg = ""
                    chat_messages = []
                    for msg in messages:
                        if msg["role"] == "system":
                            system_msg = msg["content"]
                        else:
                            chat_messages.append(msg)

                    if not chat_messages:
                        chat_messages = [{"role": "user", "content": "Hello."}]

                    kwargs: Dict[str, Any] = {
                        "model": "claude-sonnet-4-20250514",
                        "max_tokens": 2048,
                        "messages": chat_messages,
                    }
                    if system_msg:
                        kwargs["system"] = system_msg

                    with client.messages.stream(**kwargs) as stream:
                        for text in stream.text_stream:
                            yield text
                    return
            except Exception as e:
                logger.error(f"Anthropic Claude streaming failed: {e}")

        # Attempt 4: Mock streaming
        if self.mock_enabled:
            logger.info("Using mock LLM streaming (no NIM, cloud, or Claude available)")
            mock_text = self._mock_response(messages=messages)
            # Yield word by word to simulate streaming
            words = mock_text.split(" ")
            for i, word in enumerate(words):
                yield word + (" " if i < len(words) - 1 else "")
            return

        raise ConnectionError(
            "LLM unavailable: local NIM unreachable, cloud NIM failed, "
            "no Anthropic API key, and mock disabled"
        )

    def get_status(self) -> str:
        """Return service status string.

        Returns "available" for local NIM, "cloud" for NVIDIA cloud NIM,
        "mock" for mock mode, or "unavailable".
        """
        if self.is_available():
            return "available"
        if self.nvidia_api_key and self.is_cloud_available():
            return "cloud"
        if self.anthropic_api_key:
            return "anthropic"
        if self.mock_enabled:
            return "mock"
        return "unavailable"

    def _mock_response(self, **kwargs) -> str:
        """Return a template clinical response string.

        Provides a plausible clinical AI assistant response for
        demonstration and testing purposes.
        """
        messages = kwargs.get("messages", [])

        # Extract last user message for context
        last_user_msg = ""
        for msg in reversed(messages):
            if msg.get("role") == "user":
                last_user_msg = msg.get("content", "")
                break

        return (
            "Based on the available clinical information and imaging data, "
            "here is my analysis:\n\n"
            "**Assessment:**\n"
            "The imaging findings are consistent with normal anatomical "
            "structures without evidence of acute pathology. The visualized "
            "organs demonstrate expected morphology and attenuation patterns.\n\n"
            "**Key Observations:**\n"
            "1. No focal lesions or mass effects identified\n"
            "2. Vascular structures appear patent and normal in caliber\n"
            "3. No evidence of lymphadenopathy\n"
            "4. Osseous structures are intact\n\n"
            "**Recommendations:**\n"
            "- Clinical correlation is recommended\n"
            "- Follow-up imaging per standard protocol if clinically indicated\n"
            "- Integration with laboratory findings may provide additional "
            "diagnostic clarity\n\n"
            "*Note: This is an AI-generated analysis for clinical decision "
            "support. All findings should be verified by a qualified radiologist.*"
        )
