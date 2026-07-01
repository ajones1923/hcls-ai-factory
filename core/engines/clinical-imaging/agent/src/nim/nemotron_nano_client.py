"""Nemotron Nano client for efficient edge LLM inference.

Nemotron Nano is specifically designed for edge deployment with
~4 GB memory footprint and faster inference for routine queries.

Slots into the LLM fallback chain:
  Claude Sonnet (cloud, highest quality)
    -> Llama-3 8B (local, 16 GB)
    -> Nemotron Nano (local, ~4 GB, fastest)  <- NEW
    -> Mock mode

Apache 2.0 Licensed.
"""

import re
from typing import Any, Dict, List, Optional

from loguru import logger

from .base import BaseNIMClient


# ── Query classification patterns ──────────────────────────────────

ROUTINE_PATTERNS: List[str] = [
    "protocol",
    "what is the",
    "what is a",
    "define",
    "definition",
    "measurement",
    "normal range",
    "reference range",
    "device",
    "fda cleared",
    "fda approved",
    "how many",
    "list",
    "standard",
    "describe",
    "name the",
    "what does",
    "abbreviation",
    "units",
    "frequency",
]

COMPLEX_PATTERNS: List[str] = [
    "differential",
    "compare",
    " vs ",
    "versus",
    "why",
    "analyze",
    "analysis",
    "cross-modal",
    "genomic",
    "treatment plan",
    "prognosis",
    "risk assessment",
    "risk stratification",
    "correlate",
    "integrate",
    "multi-step",
    "reasoning",
    "contraindicated",
    "recommend treatment",
    "staging",
]


class NemotronNanoClient(BaseNIMClient):
    """Client for Nemotron Nano edge LLM.

    Optimized for routine clinical queries that don't need
    complex reasoning -- evidence retrieval summaries,
    protocol lookups, measurement descriptions.

    Uses the OpenAI-compatible API exposed by NIM.
    """

    def __init__(
        self,
        base_url: str,
        mock_enabled: bool = True,
        model_name: str = "nvidia/nemotron-nano",
    ):
        super().__init__(
            base_url, service_name="nemotron-nano", mock_enabled=mock_enabled
        )
        self.model_name = model_name
        self._openai_client = None

    def health_check(self) -> bool:
        """Check Nemotron Nano NIM availability via /v1/models."""
        try:
            import requests

            resp = requests.get(f"{self.base_url}/v1/models", timeout=5)
            return resp.status_code == 200
        except Exception:
            return False

    def _get_openai_client(self):
        """Lazy-initialize the OpenAI client for Nemotron Nano NIM."""
        if self._openai_client is None:
            try:
                from openai import OpenAI

                self._openai_client = OpenAI(
                    base_url=f"{self.base_url}/v1",
                    api_key="not-needed",
                )
            except ImportError:
                logger.warning(
                    "openai package not installed; Nemotron Nano NIM unavailable"
                )
                return None
        return self._openai_client

    def generate(
        self,
        prompt: str,
        max_tokens: int = 1024,
        temperature: float = 0.3,
        system_prompt: Optional[str] = None,
    ) -> str:
        """Generate text response via Nemotron Nano.

        Uses OpenAI-compatible /v1/chat/completions endpoint.

        Args:
            prompt: User query text.
            max_tokens: Maximum response tokens.
            temperature: Sampling temperature (lower = more deterministic).
            system_prompt: Optional system prompt for context.

        Returns:
            Generated text string.
        """
        messages: List[Dict[str, str]] = []
        if system_prompt:
            messages.append({"role": "system", "content": system_prompt})
        messages.append({"role": "user", "content": prompt})

        # Attempt live NIM
        if self.is_available():
            try:
                client = self._get_openai_client()
                if client is not None:
                    response = client.chat.completions.create(
                        model=self.model_name,
                        messages=messages,
                        temperature=temperature,
                        max_tokens=max_tokens,
                    )
                    text = response.choices[0].message.content
                    logger.info(
                        f"Nemotron Nano generated {len(text)} chars "
                        f"({response.usage.total_tokens} tokens)"
                    )
                    return text
            except Exception as e:
                logger.error(f"Nemotron Nano generation failed: {e}")

        # Fall back to mock
        if self.mock_enabled:
            logger.info("Using mock response for Nemotron Nano (service unavailable)")
            return self._mock_response(prompt=prompt)

        raise ConnectionError(
            "Nemotron Nano unavailable and mock disabled"
        )

    def generate_chat(
        self,
        messages: List[Dict[str, str]],
        max_tokens: int = 1024,
        temperature: float = 0.3,
    ) -> str:
        """Generate a response from a full message list.

        Args:
            messages: List of message dicts with "role" and "content" keys.
            max_tokens: Maximum response tokens.
            temperature: Sampling temperature.

        Returns:
            Generated text string.
        """
        if self.is_available():
            try:
                client = self._get_openai_client()
                if client is not None:
                    response = client.chat.completions.create(
                        model=self.model_name,
                        messages=messages,
                        temperature=temperature,
                        max_tokens=max_tokens,
                    )
                    text = response.choices[0].message.content
                    logger.info(
                        f"Nemotron Nano chat generated {len(text)} chars"
                    )
                    return text
            except Exception as e:
                logger.error(f"Nemotron Nano chat generation failed: {e}")

        if self.mock_enabled:
            logger.info("Using mock response for Nemotron Nano chat")
            last_user_msg = ""
            for msg in reversed(messages):
                if msg.get("role") == "user":
                    last_user_msg = msg.get("content", "")
                    break
            return self._mock_response(prompt=last_user_msg)

        raise ConnectionError(
            "Nemotron Nano unavailable and mock disabled"
        )

    @staticmethod
    def is_routine_query(query: str) -> bool:
        """Determine if a query is routine (suitable for Nano) or complex.

        Routine queries: protocol lookups, measurement descriptions,
        device queries, simple evidence retrieval, definitions.

        Complex queries: differential diagnosis, multi-step reasoning,
        cross-modal analysis, comparative queries, treatment planning.

        Args:
            query: The user query string.

        Returns:
            True if the query is routine and can be handled by Nano,
            False if it requires a larger model.
        """
        query_lower = query.lower()

        # Check complex patterns first (they override routine)
        for pattern in COMPLEX_PATTERNS:
            if pattern in query_lower:
                return False

        # Check routine patterns
        for pattern in ROUTINE_PATTERNS:
            if pattern in query_lower:
                return True

        # Default: treat as complex (route to larger model)
        return False

    def _mock_response(self, **kwargs) -> str:
        """Return a mock clinical response for demo/testing.

        Provides a concise clinical response appropriate for the
        routine query types that Nemotron Nano handles.
        """
        prompt = kwargs.get("prompt", "")
        prompt_lower = prompt.lower()

        if "protocol" in prompt_lower or "ct" in prompt_lower:
            return (
                "Standard CT chest protocol:\n\n"
                "- Patient position: Supine, arms above head\n"
                "- Scan range: Lung apices to adrenal glands\n"
                "- Slice thickness: 1.25 mm (axial), 2.5 mm (coronal/sagittal)\n"
                "- kVp: 120 (standard), 100 (low dose)\n"
                "- mAs: Auto-modulated (noise index 12-15)\n"
                "- Contrast: 80-100 mL iodinated, 2.5-3.0 mL/s\n"
                "- Reconstruction: Lung (B60f) and Soft Tissue (B30f) kernels\n\n"
                "*Generated by Nemotron Nano for routine protocol queries.*"
            )

        if "normal range" in prompt_lower or "measurement" in prompt_lower:
            return (
                "Normal reference ranges for the requested measurement "
                "are within established clinical guidelines. Please specify "
                "the exact measurement or anatomical structure for detailed "
                "reference values.\n\n"
                "*Generated by Nemotron Nano for routine reference queries.*"
            )

        return (
            "Based on the available clinical reference data:\n\n"
            "The requested information has been retrieved from standard "
            "clinical protocols and guidelines. The findings are consistent "
            "with established medical practice.\n\n"
            "Key points:\n"
            "1. Standard protocols apply for this clinical scenario\n"
            "2. No special considerations identified\n"
            "3. Refer to institutional guidelines for site-specific details\n\n"
            "*Generated by Nemotron Nano for routine clinical queries.*"
        )
