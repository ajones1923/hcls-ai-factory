"""
Centralized configuration for the HCLS AI Factory precision medicine platform.

Consolidates settings from:
  - core/engines/precision-intelligence/config/settings.py   (Milvus, embedding, LLM, RAG params)
  - core/engines/therapeutic-discovery/src/nim_clients.py  (NIM service configs)
  - core/agents/cart/config/settings.py  (CAR-T agent)

All values can be overridden via environment variables (uppercased, prefixed with
HCLS_ for the main block, or via the .env file).  A cached singleton is exposed
via ``get_settings()``.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path
from typing import List, Optional

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Settings model
# ---------------------------------------------------------------------------

class HCLSSettings(BaseSettings):
    """
    Unified settings for every HCLS AI Factory component.

    Categories:
        service       -- service identity
        database      -- general DB knobs
        milvus        -- vector store connection
        embedding     -- model, dimension, batching
        llm           -- provider, model, key, generation params
        nim_services  -- MolMIM / DiffDock endpoints
        pipeline      -- worker / batch tuning
        observability -- Prometheus, OTLP
        security      -- rate limits, token budgets
    """

    model_config = SettingsConfigDict(
        env_prefix="HCLS_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
        case_sensitive=False,
    )

    # ── Service Identity ──────────────────────────────────────────────────
    service_name: str = Field(
        default="hcls-ai-factory",
        description="Logical service name used in tracing spans and metrics labels",
    )

    # ── Milvus ────────────────────────────────────────────────────────────
    milvus_host: str = Field(default="localhost", description="Milvus gRPC host")
    milvus_port: int = Field(default=19530, description="Milvus gRPC port")
    milvus_collection: str = Field(
        default="genomic_evidence",
        description="Default Milvus collection for genomic variant evidence",
    )
    milvus_pool_size: int = Field(
        default=4, description="Max simultaneous Milvus connections"
    )
    milvus_timeout: int = Field(
        default=30, description="Milvus operation timeout in seconds"
    )

    # ── Embedding ─────────────────────────────────────────────────────────
    embedding_model: str = Field(
        default="BAAI/bge-small-en-v1.5",
        description="HuggingFace model ID for local sentence-transformers embedder",
    )
    embedding_dimension: int = Field(default=384, description="Vector dimension")
    embedding_batch_size: int = Field(
        default=32, description="Batch size for encoding"
    )
    embedding_provider: str = Field(
        default="local_bge",
        description="Embedding backend: local_bge | tei | openai",
    )
    embedding_device: Optional[str] = Field(
        default=None,
        description="PyTorch device for local embedder (cpu, cuda, cuda:0, ...)",
    )
    tei_endpoint: Optional[str] = Field(
        default=None,
        description="Text Embeddings Inference HTTP endpoint (e.g. http://localhost:8081)",
    )
    embedding_cache_dir: Optional[Path] = Field(
        default=None,
        description="Directory for disk-cached embeddings (variant IDs -> vectors)",
    )

    # ── LLM ───────────────────────────────────────────────────────────────
    llm_provider: str = Field(
        default="ollama",
        description="LLM backend: anthropic | openai | ollama | vllm",
    )
    llm_model: str = Field(
        default="llama3.1:70b",
        description="Model identifier (provider-specific)",
    )
    llm_api_key: Optional[str] = Field(
        default=None,
        description="API key for Anthropic / OpenAI (also reads ANTHROPIC_API_KEY / OPENAI_API_KEY)",
    )
    llm_max_tokens: int = Field(default=4096, description="Max output tokens")
    llm_temperature: float = Field(default=0.7, ge=0.0, le=2.0)

    # Ollama-specific
    ollama_host: str = Field(default="http://localhost:11434")

    # vLLM-specific
    vllm_host: str = Field(default="localhost")
    vllm_port: int = Field(default=8080)
    vllm_model: str = Field(default="meta-llama/Llama-3.1-70B-Instruct")

    # LLM resilience
    llm_retry_base: float = Field(
        default=2.0,
        description="Base delay (seconds) for exponential backoff on 429/529",
    )
    llm_retry_max: float = Field(
        default=32.0, description="Max backoff delay in seconds"
    )
    llm_daily_token_budget: int = Field(
        default=0,
        description="Daily token budget (0 = unlimited)",
    )
    llm_response_cache_ttl: int = Field(
        default=300,
        description="Response cache TTL in seconds (0 = disabled)",
    )
    llm_response_cache_max: int = Field(
        default=256,
        description="Maximum entries in LLM response cache",
    )

    # ── NIM Services (Drug Discovery) ─────────────────────────────────────
    nim_molmim_host: str = Field(default="localhost")
    nim_molmim_port: int = Field(default=8001)
    nim_diffdock_host: str = Field(default="localhost")
    nim_diffdock_port: int = Field(default=8002)
    nim_api_version: str = Field(default="v1")
    nim_timeout: int = Field(default=300, description="NIM request timeout (s)")
    nim_max_retries: int = Field(default=5)

    # ── RAG Parameters ────────────────────────────────────────────────────
    rag_top_k: int = Field(default=10, description="Evidence items to retrieve")
    rag_score_threshold: float = Field(
        default=0.5, ge=0.0, le=1.0, description="Minimum cosine similarity"
    )

    # ── Pipeline Tuning ───────────────────────────────────────────────────
    pipeline_workers: int = Field(default=4, ge=1, description="Parallel workers")
    pipeline_batch_size: int = Field(
        default=1000, ge=1, description="Ingestion batch size"
    )

    # ── Observability ─────────────────────────────────────────────────────
    metrics_enabled: bool = Field(default=True, description="Enable Prometheus metrics")
    metrics_port: int = Field(
        default=9099, description="Prometheus metrics endpoint port"
    )
    otlp_endpoint: Optional[str] = Field(
        default=None,
        description="OpenTelemetry collector gRPC endpoint (e.g. http://localhost:4317)",
    )

    # ── Security ──────────────────────────────────────────────────────────
    rate_limit_max_requests: int = Field(
        default=100, description="Max requests per rate-limit window"
    )
    rate_limit_window_seconds: int = Field(
        default=60, description="Rate-limit sliding window in seconds"
    )

    # ── Quality Filters (Genomics) ────────────────────────────────────────
    min_variant_qual: float = Field(default=20.0, description="Min QUAL score")
    include_chromosomes: List[str] = Field(
        default=[
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
            "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
            "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
            "chrX", "chrY",
        ],
        description="Chromosomes to include during ingestion",
    )

    # ── Validators ────────────────────────────────────────────────────────

    @field_validator("llm_provider")
    @classmethod
    def _validate_llm_provider(cls, v: str) -> str:
        allowed = {"anthropic", "openai", "ollama", "vllm"}
        v = v.lower()
        if v not in allowed:
            raise ValueError(
                f"llm_provider must be one of {allowed}, got {v!r}"
            )
        return v

    @field_validator("embedding_provider")
    @classmethod
    def _validate_embedding_provider(cls, v: str) -> str:
        allowed = {"local_bge", "tei", "openai"}
        v = v.lower()
        if v not in allowed:
            raise ValueError(
                f"embedding_provider must be one of {allowed}, got {v!r}"
            )
        return v

    # ── Convenience properties ────────────────────────────────────────────

    @property
    def nim_molmim_url(self) -> str:
        """Full MolMIM base URL."""
        return f"http://{self.nim_molmim_host}:{self.nim_molmim_port}/{self.nim_api_version}"

    @property
    def nim_diffdock_url(self) -> str:
        """Full DiffDock base URL."""
        return f"http://{self.nim_diffdock_host}:{self.nim_diffdock_port}/{self.nim_api_version}"


# ---------------------------------------------------------------------------
# Singleton accessor
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def get_settings() -> HCLSSettings:
    """
    Return a cached singleton ``HCLSSettings`` instance.

    The first call reads from env / .env; subsequent calls return the same
    object.  To force a reload (e.g. in tests) call
    ``get_settings.cache_clear()`` first.
    """
    settings = HCLSSettings()
    logger.info(
        "HCLSSettings loaded  service=%s  milvus=%s:%d  llm=%s/%s  embed=%s",
        settings.service_name,
        settings.milvus_host,
        settings.milvus_port,
        settings.llm_provider,
        settings.llm_model,
        settings.embedding_model,
    )
    return settings
