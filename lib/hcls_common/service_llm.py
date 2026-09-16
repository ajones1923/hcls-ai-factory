"""One LLM client for every service that generates clinical prose.

Eight services each carried their own `_LLMClient` class, byte-for-byte alike apart from which
`settings.LLM_MODEL` they read. That is not a tidiness problem — it is why a single API change
cost eight edits and got seven of them:

  * `temperature` was REMOVED on the current models (Sonnet 5, Opus 5, the 4.6+ family). Sending
    it returns 400 "`temperature` is deprecated for this model." Every agent sent it, every
    agent's synthesis 400'd, and each route caught the error and fell back to a stub answer —
    so the services looked healthy and returned "Search completed. See evidence passages below."
    The shared `hcls_common.llm_client` was fixed first; these eight copies were missed, and the
    clinical eval is what caught it.

  * Model ids drift. 37 references to three dead models were live in this repo at once.

Both belong in one place. Per-service behaviour that genuinely differs — the default system prompt, the model choice —
stays a parameter. Five of the eight baked in their own system prompt; that is the only real
variation among them.
"""
from __future__ import annotations

import logging
from typing import Any, Iterator
from hcls_common.llm_text import first_text

logger = logging.getLogger(__name__)

# Models on which temperature / top_p / top_k are removed and return 400 if sent.
_NO_SAMPLING_PREFIXES = (
    "claude-opus-5", "claude-opus-4-8", "claude-opus-4-7", "claude-opus-4-6",
    "claude-sonnet-5", "claude-sonnet-4-6", "claude-fable-5", "claude-mythos-5",
)


def rejects_sampling(model: str) -> bool:
    """True when `model` rejects temperature/top_p/top_k."""
    m = (model or "").lower()
    return any(m.startswith(p) for p in _NO_SAMPLING_PREFIXES)


class ServiceLLM:
    """Thin wrapper over the Anthropic client with the sampling-parameter rule applied.

    Deliberately keeps the `temperature` argument in its own signature so callers written
    against the old shape keep working; it is simply not forwarded to models that reject it.
    """

    def __init__(self, model: str, *, client: Any = None, default_system: str = ""):
        self.model = model
        self.default_system = default_system      # the real per-service difference
        if client is None:                        # only the real path needs the SDK,
            import anthropic                      # so an injected client tests without it
            client = anthropic.Anthropic()
        self.client = client

    def _kwargs(self, prompt: str, system_prompt: str, max_tokens: int,
                temperature: float | None) -> dict:
        kw = {
            "model": self.model,
            "max_tokens": max_tokens,
            "system": system_prompt or self.default_system or "",
            "messages": [{"role": "user", "content": prompt}],
        }
        if temperature is not None and not rejects_sampling(self.model):
            kw["temperature"] = temperature
        return kw

    def generate(self, prompt: str, system_prompt: str = "",
                 max_tokens: int = 2048, temperature: float | None = 0.7) -> str:
        msg = self.client.messages.create(
            **self._kwargs(prompt, system_prompt, max_tokens, temperature))
        return first_text(msg)

    def generate_stream(self, prompt: str, system_prompt: str = "",
                        max_tokens: int = 2048,
                        temperature: float | None = 0.7) -> Iterator[str]:
        with self.client.messages.stream(
                **self._kwargs(prompt, system_prompt, max_tokens, temperature)) as stream:
            yield from stream.text_stream


def build_service_llm(model: str, *, service: str = "",
                      default_system: str = "") -> ServiceLLM | None:
    """Return a client, or None when the SDK or a credential is missing.

    Returns None rather than raising so a service still starts and degrades to retrieval-only —
    the posture every agent already had. The log line is the signal that synthesis is off; a
    silent None is how "healthy but answering with a stub" happened in the first place.
    """
    try:
        llm = ServiceLLM(model, default_system=default_system)
    except Exception as exc:
        logger.warning("LLM client unavailable for %s (%s): synthesis disabled, "
                       "retrieval still works", service or "service", exc)
        return None
    logger.info("Anthropic LLM client initialized (model=%s%s)", model,
                ", sampling params omitted" if rejects_sampling(model) else "")
    return llm
