"""Cardiology Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as the PGx agent.
"""

import logging
import os
from pathlib import Path
from typing import List, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


class CardioSettings(BaseSettings):
    """Configuration for Cardiology Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── RAG Pipeline (reuse existing) ──
    RAG_PIPELINE_ROOT: Path = Path(
        os.environ.get("CARDIO_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")
    )

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (12 cardiology-specific + 1 shared genomic_evidence)
    COLLECTION_LITERATURE: str = "cardio_literature"
    COLLECTION_TRIALS: str = "cardio_trials"
    COLLECTION_IMAGING: str = "cardio_imaging"
    COLLECTION_ELECTROPHYSIOLOGY: str = "cardio_electrophysiology"
    COLLECTION_HEART_FAILURE: str = "cardio_heart_failure"
    COLLECTION_VALVULAR: str = "cardio_valvular"
    COLLECTION_PREVENTION: str = "cardio_prevention"
    COLLECTION_INTERVENTIONAL: str = "cardio_interventional"
    COLLECTION_ONCOLOGY: str = "cardio_oncology"
    COLLECTION_DEVICES: str = "cardio_devices"
    COLLECTION_GUIDELINES: str = "cardio_guidelines"
    COLLECTION_HEMODYNAMICS: str = "cardio_hemodynamics"
    COLLECTION_GENOMIC: str = "genomic_evidence"  # Existing shared collection

    # ── Embeddings ──
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIMENSION: int = 384
    EMBEDDING_BATCH_SIZE: int = 32

    # ── LLM ──
    LLM_PROVIDER: str = "anthropic"
    LLM_MODEL: str = "claude-sonnet-4-6"
    ANTHROPIC_API_KEY: Optional[str] = None

    # ── RAG Search ──
    TOP_K_PER_COLLECTION: int = 5
    SCORE_THRESHOLD: float = 0.4

    # Collection search weights (must sum to ~1.0)
    WEIGHT_LITERATURE: float = 0.10
    WEIGHT_TRIALS: float = 0.08
    WEIGHT_IMAGING: float = 0.10
    WEIGHT_ELECTROPHYSIOLOGY: float = 0.08
    WEIGHT_HEART_FAILURE: float = 0.10
    WEIGHT_VALVULAR: float = 0.08
    WEIGHT_PREVENTION: float = 0.10
    WEIGHT_INTERVENTIONAL: float = 0.07
    WEIGHT_ONCOLOGY: float = 0.06
    WEIGHT_DEVICES: float = 0.04
    WEIGHT_GUIDELINES: float = 0.10
    WEIGHT_HEMODYNAMICS: float = 0.06
    WEIGHT_GENOMIC: float = 0.03

    # ── PubMed ──
    NCBI_API_KEY: Optional[str] = None
    PUBMED_MAX_RESULTS: int = 5000

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8126

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8536

    # ── Prometheus Metrics ──
    METRICS_ENABLED: bool = True

    # ── Scheduler ──
    INGEST_SCHEDULE_HOURS: int = 168  # Weekly (7 * 24)
    INGEST_ENABLED: bool = False

    # ── Conversation Memory ──
    MAX_CONVERSATION_CONTEXT: int = 3

    # ── Citation Scoring ──
    CITATION_HIGH_THRESHOLD: float = 0.75
    CITATION_MEDIUM_THRESHOLD: float = 0.60

    # ── Cross-Agent Integration ──
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    NEUROLOGY_AGENT_URL: str = "http://localhost:8528"
    IMAGING_AGENT_URL: str = "http://localhost:8524"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── Authentication ──
    API_KEY: str = ""  # Empty = no auth required

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8126,http://localhost:8536"

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="CARDIO_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )

    # ── Startup Validation ──

    def validate(self) -> List[str]:
        """Return a list of configuration warnings/errors (never raises)."""
        issues: List[str] = []

        if not self.MILVUS_HOST or not self.MILVUS_HOST.strip():
            issues.append("MILVUS_HOST is empty — Milvus connections will fail.")
        if not (1 <= self.MILVUS_PORT <= 65535):
            issues.append(
                f"MILVUS_PORT={self.MILVUS_PORT} is outside valid range (1-65535)."
            )

        if not self.ANTHROPIC_API_KEY:
            issues.append(
                "ANTHROPIC_API_KEY is not set — LLM features disabled, "
                "search-only mode available."
            )

        if not self.EMBEDDING_MODEL or not self.EMBEDDING_MODEL.strip():
            issues.append("EMBEDDING_MODEL is empty — embedding pipeline will fail.")

        for name, port in [("API_PORT", self.API_PORT), ("STREAMLIT_PORT", self.STREAMLIT_PORT)]:
            if not (1024 <= port <= 65535):
                issues.append(
                    f"{name}={port} is outside valid range (1024-65535)."
                )
        if self.API_PORT == self.STREAMLIT_PORT:
            issues.append(
                f"API_PORT and STREAMLIT_PORT are both {self.API_PORT} — port conflict."
            )

        weight_attrs = [
            attr for attr in dir(self)
            if attr.startswith("WEIGHT_") and isinstance(getattr(self, attr), float)
        ]
        weights = []
        for attr in weight_attrs:
            val = getattr(self, attr)
            if val < 0:
                issues.append(f"{attr}={val} is negative — weights must be >= 0.")
            weights.append(val)
        if weights:
            total = sum(weights)
            if abs(total - 1.0) > 0.05:
                issues.append(
                    f"Collection weights sum to {total:.4f}, expected ~1.0 "
                    f"(tolerance 0.05)."
                )

        if self.RAG_PIPELINE_ROOT and not self.RAG_PIPELINE_ROOT.is_dir():
            issues.append(
                f"RAG_PIPELINE_ROOT={self.RAG_PIPELINE_ROOT} does not exist or "
                f"is not a directory."
            )

        return issues

    def validate_or_warn(self) -> None:
        """Run validate() and log each issue as a warning."""
        for issue in self.validate():
            logger.warning("Cardio config: %s", issue)


settings = CardioSettings()
settings.validate_or_warn()
