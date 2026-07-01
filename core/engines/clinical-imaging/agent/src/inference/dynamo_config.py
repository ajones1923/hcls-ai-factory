"""NVIDIA Dynamo inference optimization for LLM serving.

Provides disaggregated prefill/decode scheduling and LLM-aware
request routing for optimized LLM serving on DGX Spark.

Dynamo is NVIDIA's open-source framework for distributed and
optimized LLM inference. On DGX Spark (single GPU), it provides
KV cache management and prefill/decode disaggregation for better
GPU utilization during multi-turn medical conversations.

Apache 2.0 Licensed.
"""

from typing import Any, Dict, Optional

from loguru import logger
from pydantic import BaseModel


# ── Preset configurations per model ────────────────────────────────

MODEL_PRESETS: Dict[str, Dict[str, Any]] = {
    "llama-3-8b": {
        "prefill_batch_size": 4,
        "decode_batch_size": 8,
        "kv_cache_max_gb": 16.0,
        "scheduler": "prefill_priority",
    },
    "nemotron-nano": {
        "prefill_batch_size": 8,
        "decode_batch_size": 16,
        "kv_cache_max_gb": 4.0,
        "scheduler": "round_robin",
    },
}

# KV cache bytes per token per layer per head (fp16)
# Formula: 2 (K+V) * num_layers * num_kv_heads * head_dim * 2 (fp16 bytes)
_KV_CACHE_PARAMS: Dict[str, Dict[str, int]] = {
    "llama-3-8b": {
        "num_layers": 32,
        "num_kv_heads": 8,
        "head_dim": 128,
    },
    "nemotron-nano": {
        "num_layers": 24,
        "num_kv_heads": 8,
        "head_dim": 64,
    },
}

VALID_SCHEDULERS = ("round_robin", "load_balanced", "prefill_priority")


class DynamoConfig(BaseModel):
    """Configuration for NVIDIA Dynamo LLM inference optimization."""

    enabled: bool = False
    prefill_batch_size: int = 4
    decode_batch_size: int = 8
    kv_cache_max_gb: float = 16.0
    scheduler: str = "round_robin"


class DynamoStatus(BaseModel):
    """Runtime status of the Dynamo optimizer."""

    available: bool = False
    scheduler: str = "none"
    active_requests: int = 0
    kv_cache_usage_gb: float = 0.0
    prefill_queue_depth: int = 0
    decode_queue_depth: int = 0
    avg_ttft_ms: float = 0.0
    avg_tps: float = 0.0


class DynamoOptimizer:
    """NVIDIA Dynamo LLM inference optimizer.

    On DGX Spark (single GPU), provides:
    - Prefill/decode disaggregation for better GPU utilization
    - KV cache management for multi-turn medical conversations
    - Request scheduling optimized for RAG synthesis workloads

    Falls back to standard inference when Dynamo is not available.
    """

    def __init__(self, config: Optional[DynamoConfig] = None):
        self.config = config or DynamoConfig()
        self._dynamo_available = self._check_dynamo()
        if self._dynamo_available:
            logger.info(
                f"NVIDIA Dynamo available — scheduler={self.config.scheduler}, "
                f"kv_cache_max={self.config.kv_cache_max_gb:.1f} GB"
            )
        else:
            logger.info(
                "NVIDIA Dynamo not available — using standard inference. "
                "Install with: pip install nvidia-dynamo"
            )

    def _check_dynamo(self) -> bool:
        """Check whether the dynamo runtime is importable."""
        try:
            import dynamo  # noqa: F401

            return True
        except ImportError:
            return False

    def get_status(self) -> DynamoStatus:
        """Return current Dynamo optimizer status.

        When Dynamo is unavailable, returns a status object with
        available=False and zeroed metrics.
        """
        if not self._dynamo_available or not self.config.enabled:
            return DynamoStatus(
                available=False,
                scheduler="none",
            )

        # When Dynamo is available, query runtime metrics
        try:
            import dynamo  # noqa: F401

            # In a real deployment, these would come from the Dynamo
            # runtime API. For now, return the configured state.
            return DynamoStatus(
                available=True,
                scheduler=self.config.scheduler,
                active_requests=0,
                kv_cache_usage_gb=0.0,
                prefill_queue_depth=0,
                decode_queue_depth=0,
                avg_ttft_ms=0.0,
                avg_tps=0.0,
            )
        except Exception as e:
            logger.error(f"Failed to query Dynamo status: {e}")
            return DynamoStatus(available=False, scheduler="none")

    def optimize_request(
        self, prompt: str, max_tokens: int
    ) -> Dict[str, Any]:
        """Return optimization hints for a given request.

        These hints can be passed to the serving backend to influence
        scheduling, batching, and KV cache allocation.

        Args:
            prompt: The input prompt text.
            max_tokens: Maximum tokens to generate.

        Returns:
            Dict of optimization hints. When Dynamo is unavailable,
            returns an empty dict (standard inference path).
        """
        if not self._dynamo_available or not self.config.enabled:
            return {}

        prompt_tokens = len(prompt) // 4  # rough estimate
        estimated_kv_gb = self.estimate_kv_cache_size(
            context_length=prompt_tokens + max_tokens,
            model_params_b=8.0,
        )

        # Choose scheduling priority based on prompt length
        if prompt_tokens > 2048:
            priority = "prefill_priority"
        else:
            priority = self.config.scheduler

        hints: Dict[str, Any] = {
            "scheduler": priority,
            "prefill_batch_size": self.config.prefill_batch_size,
            "decode_batch_size": self.config.decode_batch_size,
            "estimated_kv_cache_gb": round(estimated_kv_gb, 3),
            "prompt_tokens_est": prompt_tokens,
            "max_tokens": max_tokens,
            "disaggregate_prefill": True,
        }

        logger.debug(
            f"Dynamo optimization hints: scheduler={priority}, "
            f"kv_cache_est={estimated_kv_gb:.3f} GB, "
            f"prompt_tokens~{prompt_tokens}"
        )
        return hints

    def estimate_kv_cache_size(
        self, context_length: int, model_params_b: float = 8.0
    ) -> float:
        """Estimate KV cache memory in GB for a given context length.

        Uses model-specific parameters when available, otherwise
        falls back to a formula based on total parameter count.

        Args:
            context_length: Total sequence length (prompt + generation).
            model_params_b: Model size in billions of parameters.

        Returns:
            Estimated KV cache size in GB.
        """
        # Try to find a matching preset for precise calculation
        for model_key, params in _KV_CACHE_PARAMS.items():
            preset = MODEL_PRESETS.get(model_key, {})
            preset_params_b = {
                "llama-3-8b": 8.0,
                "nemotron-nano": 4.0,
            }.get(model_key, 0.0)

            if abs(model_params_b - preset_params_b) < 0.5:
                # Precise calculation: 2 * num_layers * num_kv_heads * head_dim * 2 bytes * context_length
                bytes_per_token = (
                    2
                    * params["num_layers"]
                    * params["num_kv_heads"]
                    * params["head_dim"]
                    * 2  # fp16
                )
                total_bytes = bytes_per_token * context_length
                return total_bytes / (1024**3)

        # Generic fallback: ~0.5 GB per billion params per 1K context tokens
        gb = (model_params_b * 0.5 * context_length) / 1000.0
        return gb

    def get_recommended_config(self, model_name: str) -> DynamoConfig:
        """Return a recommended DynamoConfig for a known model.

        Args:
            model_name: Model identifier. Supported presets:
                - "llama-3-8b": Llama-3 8B (16 GB KV cache)
                - "nemotron-nano": Nemotron Nano (~4 GB KV cache)

        Returns:
            DynamoConfig with preset values, or default config if
            the model is not recognized.
        """
        preset = MODEL_PRESETS.get(model_name)
        if preset is None:
            logger.warning(
                f"No Dynamo preset for model '{model_name}'; "
                f"using defaults. Known presets: {list(MODEL_PRESETS.keys())}"
            )
            return DynamoConfig(enabled=self._dynamo_available)

        return DynamoConfig(
            enabled=self._dynamo_available,
            prefill_batch_size=preset["prefill_batch_size"],
            decode_batch_size=preset["decode_batch_size"],
            kv_cache_max_gb=preset["kv_cache_max_gb"],
            scheduler=preset["scheduler"],
        )
