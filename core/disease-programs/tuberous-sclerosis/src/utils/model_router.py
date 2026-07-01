"""
Tiered model router (PRD §2.2; master paper §16).

Reads config/model_policy.yaml and routes an agent step to its tier. When a real
Anthropic key is configured it calls the API; otherwise (offline / no key / tests)
it returns a deterministic, watermarked stub so the engine runs end-to-end without
network. Either way it returns the text plus the provenance fields the event log
records. Deterministic and classical "tiers" are handled by the agents directly;
this router covers the LLM tiers (haiku/sonnet/opus/local).
"""
from __future__ import annotations

import json
import time

import yaml

from config.settings import settings
from src.utils.cost import estimate_tokens, get_ledger
from src.utils.provenance import Provenance, stable_hash

_LLM_TIERS = {"haiku", "sonnet", "opus", "local"}


class ModelRouter:
    def __init__(self, policy_path=None) -> None:
        path = policy_path or settings.MODEL_POLICY_PATH
        with open(path) as f:
            self._policy = yaml.safe_load(f)
        self._models = self._policy["models"]
        self._client = None
        self._online = bool(settings.ANTHROPIC_API_KEY) and not settings.OFFLINE
        if self._online:
            try:
                import anthropic  # noqa: WPS433

                self._client = anthropic.Anthropic(api_key=settings.ANTHROPIC_API_KEY)
            except Exception:  # pragma: no cover - falls back to stub
                self._online = False
        get_ledger().set_cap(getattr(settings, "MAX_RUN_USD", 0.0))

    def tier_for(self, agent: str, step: str) -> str:
        return self._policy["agents"][agent][step]

    @staticmethod
    def _meter(tier: str, sent: str, received: str, resp, billed: bool) -> None:
        """Record token + cost usage for one call. Uses real API usage when available,
        a ~4-chars/token estimate otherwise. Offline stub calls are counted but never billed."""
        if resp is not None and getattr(resp, "usage", None) is not None:
            in_tok = int(getattr(resp.usage, "input_tokens", 0) or 0)
            out_tok = int(getattr(resp.usage, "output_tokens", 0) or 0)
        else:
            in_tok, out_tok = estimate_tokens(sent), estimate_tokens(received)
        get_ledger().record(tier, in_tok, out_tok, billed=billed)

    def call(
        self,
        agent: str,
        step: str,
        system: str,
        prompt: str,
        prompt_template_version: str = "v0",
        rag_source_uris: list[str] | None = None,
        max_tokens: int = 2048,
    ) -> tuple[str, Provenance]:
        """Run one LLM step. Returns (text, provenance)."""
        tier = self.tier_for(agent, step)
        if tier not in _LLM_TIERS:
            raise ValueError(
                f"{agent}.{step} is tier '{tier}' (deterministic/classical); "
                "it must not be routed through the LLM model router."
            )
        model_id = self._models[tier]
        input_hash = stable_hash({"system": system, "prompt": prompt})
        t0 = time.perf_counter()
        if self._online and self._client is not None:
            resp = self._client.messages.create(
                model=model_id,
                max_tokens=max_tokens,
                # Prompt caching (COST-12): the per-agent system prompt is reused across
                # all 50 cohort patients, so caching it cuts cached-input cost ~90% on a
                # full run. GA (no beta header); no-ops silently below the model minimum.
                system=[{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}],
                messages=[{"role": "user", "content": prompt}],
            )
            text = "".join(b.text for b in resp.content if getattr(b, "type", None) == "text")
            model_version = getattr(resp, "model", model_id)
            self._meter(tier, system + prompt, text, resp, billed=True)
        else:
            text = self._stub(agent, step, tier)
            model_version = f"{model_id} (OFFLINE STUB)"
            self._meter(tier, system + prompt, text, None, billed=False)
        latency_ms = (time.perf_counter() - t0) * 1000.0
        prov = Provenance(
            agent=agent,
            agent_version=settings.AGENT_VERSION,
            step=step,
            tier=tier,
            model_id=model_id,
            model_version=model_version,
            prompt_template_version=prompt_template_version,
            rag_source_uris=rag_source_uris or [],
            input_hash=input_hash,
            latency_ms=round(latency_ms, 2),
        )
        return text, prov

    @property
    def online(self) -> bool:
        return self._online and self._client is not None

    def call_structured(
        self,
        agent: str,
        step: str,
        system: str,
        prompt: str,
        schema_hint: dict,
        prompt_template_version: str = "v0",
        rag_source_uris: list[str] | None = None,
        max_tokens: int = 2048,
    ) -> tuple[dict | None, Provenance]:
        """Run an LLM step that must return JSON. Returns (parsed_dict_or_None, provenance).
        Offline (no key) returns None so the agent keeps its deterministic fallback."""
        tier = self.tier_for(agent, step)
        if tier not in _LLM_TIERS:
            raise ValueError(f"{agent}.{step} is tier '{tier}'; not an LLM step.")
        model_id = self._models[tier]
        input_hash = stable_hash({"system": system, "prompt": prompt})
        t0 = time.perf_counter()
        parsed: dict | None = None
        if self.online:
            full_system = (
                system + "\n\nReturn ONLY a single JSON object of exactly this shape "
                "(no prose, no code fence):\n" + json.dumps(schema_hint)
            )
            resp = self._client.messages.create(
                model=model_id, max_tokens=max_tokens,
                system=[{"type": "text", "text": full_system,            # COST-12 prompt caching
                         "cache_control": {"type": "ephemeral"}}],
                messages=[{"role": "user", "content": prompt}],
            )
            text = "".join(b.text for b in resp.content if getattr(b, "type", None) == "text")
            parsed = _extract_json(text)
            model_version = getattr(resp, "model", model_id)
            self._meter(tier, full_system + prompt, text, resp, billed=True)
        else:
            model_version = f"{model_id} (OFFLINE STUB)"
            self._meter(tier, system + prompt, "", None, billed=False)
        latency_ms = (time.perf_counter() - t0) * 1000.0
        prov = Provenance(
            agent=agent, agent_version=settings.AGENT_VERSION, step=step, tier=tier,
            model_id=model_id, model_version=model_version,
            prompt_template_version=prompt_template_version,
            rag_source_uris=rag_source_uris or [], input_hash=input_hash,
            latency_ms=round(latency_ms, 2),
        )
        return parsed, prov

    @staticmethod
    def _stub(agent: str, step: str, tier: str) -> str:
        return (
            f"[OFFLINE STUB — {tier} tier] {agent}.{step}: deterministic placeholder "
            "output for the walking-skeleton build. Real reasoning is produced when "
            "TSC_ANTHROPIC_API_KEY is configured. SYNTHETIC demonstration data."
        )


def _extract_json(text: str) -> dict | None:
    """Best-effort: pull the first JSON object out of a model response."""
    if not text:
        return None
    t = text.strip()
    if t.startswith("```"):
        t = t.strip("`")
    i, j = t.find("{"), t.rfind("}")
    if i == -1 or j == -1 or j < i:
        return None
    try:
        obj = json.loads(t[i:j + 1])
        return obj if isinstance(obj, dict) else None
    except Exception:
        return None


_router: ModelRouter | None = None


def get_router() -> ModelRouter:
    global _router
    if _router is None:
        _router = ModelRouter()
    return _router
