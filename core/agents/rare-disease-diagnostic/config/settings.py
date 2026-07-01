"""Rare Disease Diagnostic Agent configuration.

Follows the same Pydantic BaseSettings pattern as the Clinical Trial agent.

Author: Adam Jones
Date: March 2026
"""

import logging
from pathlib import Path
from typing import List, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


class RareDiseaseSettings(BaseSettings):
    """Configuration for the Rare Disease Diagnostic Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (14 rare-disease-specific collections)
    COLLECTION_PHENOTYPES: str = "rd_phenotypes"
    COLLECTION_DISEASES: str = "rd_diseases"
    COLLECTION_GENES: str = "rd_genes"
    COLLECTION_VARIANTS: str = "rd_variants"
    COLLECTION_LITERATURE: str = "rd_literature"
    COLLECTION_TRIALS: str = "rd_trials"
    COLLECTION_THERAPIES: str = "rd_therapies"
    COLLECTION_CASE_REPORTS: str = "rd_case_reports"
    COLLECTION_GUIDELINES: str = "rd_guidelines"
    COLLECTION_PATHWAYS: str = "rd_pathways"
    COLLECTION_REGISTRIES: str = "rd_registries"
    COLLECTION_NATURAL_HISTORY: str = "rd_natural_history"
    COLLECTION_NEWBORN_SCREENING: str = "rd_newborn_screening"
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
    SCORE_THRESHOLD: float = 0.4

    # Per-collection TOP_K defaults
    TOP_K_PHENOTYPES: int = 50
    TOP_K_DISEASES: int = 30
    TOP_K_GENES: int = 30
    TOP_K_VARIANTS: int = 100
    TOP_K_LITERATURE: int = 20
    TOP_K_CASE_REPORTS: int = 15
    TOP_K_THERAPIES: int = 10
    TOP_K_GUIDELINES: int = 10
    TOP_K_TRIALS: int = 10
    TOP_K_PATHWAYS: int = 10
    TOP_K_REGISTRIES: int = 10
    TOP_K_NATURAL_HISTORY: int = 10
    TOP_K_NEWBORN_SCREENING: int = 10
    TOP_K_GENOMIC: int = 20

    # Collection search weights (must sum to ~1.0)
    WEIGHT_PHENOTYPES: float = 0.12
    WEIGHT_DISEASES: float = 0.11
    WEIGHT_GENES: float = 0.10
    WEIGHT_VARIANTS: float = 0.10
    WEIGHT_LITERATURE: float = 0.08
    WEIGHT_CASE_REPORTS: float = 0.07
    WEIGHT_THERAPIES: float = 0.07
    WEIGHT_TRIALS: float = 0.06
    WEIGHT_GUIDELINES: float = 0.06
    WEIGHT_PATHWAYS: float = 0.06
    WEIGHT_NATURAL_HISTORY: float = 0.05
    WEIGHT_REGISTRIES: float = 0.04
    WEIGHT_NEWBORN_SCREENING: float = 0.05
    WEIGHT_GENOMIC: float = 0.03

    # ── External APIs ──
    ORPHANET_API_KEY: Optional[str] = None
    NCBI_API_KEY: Optional[str] = None

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8134

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8544

    # ── Prometheus Metrics ──
    METRICS_ENABLED: bool = True

    # ── Scheduler ──
    INGEST_SCHEDULE_HOURS: int = 24
    INGEST_ENABLED: bool = False

    # ── Conversation Memory ──
    MAX_CONVERSATION_CONTEXT: int = 3

    # ── Citation Scoring ──
    CITATION_HIGH_THRESHOLD: float = 0.75
    CITATION_MEDIUM_THRESHOLD: float = 0.60

    # ── Authentication ──
    API_KEY: str = ""  # Empty = no auth required

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8134,http://localhost:8544"

    # ── Cross-Agent Integration ──
    GENOMICS_AGENT_URL: str = "http://localhost:8527"
    PGX_AGENT_URL: str = "http://localhost:8107"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    IMAGING_AGENT_URL: str = "http://localhost:8524"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="RD_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )

    # ── Startup Validation ──

    def validate(self) -> List[str]:
        """Return a list of configuration warnings/errors (never raises)."""
        issues: List[str] = []

        if not self.MILVUS_HOST or not self.MILVUS_HOST.strip():
            issues.append("MILVUS_HOST is empty -- Milvus connections will fail.")
        if not (1 <= self.MILVUS_PORT <= 65535):
            issues.append(
                f"MILVUS_PORT={self.MILVUS_PORT} is outside valid range (1-65535)."
            )

        if not self.ANTHROPIC_API_KEY:
            issues.append(
                "ANTHROPIC_API_KEY is not set -- LLM features disabled, "
                "search-only mode available."
            )

        if not self.EMBEDDING_MODEL or not self.EMBEDDING_MODEL.strip():
            issues.append("EMBEDDING_MODEL is empty -- embedding pipeline will fail.")

        for name, port in [("API_PORT", self.API_PORT), ("STREAMLIT_PORT", self.STREAMLIT_PORT)]:
            if not (1024 <= port <= 65535):
                issues.append(
                    f"{name}={port} is outside valid range (1024-65535)."
                )
        if self.API_PORT == self.STREAMLIT_PORT:
            issues.append(
                f"API_PORT and STREAMLIT_PORT are both {self.API_PORT} -- port conflict."
            )

        weight_attrs = [
            attr for attr in dir(self)
            if attr.startswith("WEIGHT_") and isinstance(getattr(self, attr), float)
        ]
        weights = []
        for attr in weight_attrs:
            val = getattr(self, attr)
            if val < 0:
                issues.append(f"{attr}={val} is negative -- weights must be >= 0.")
            weights.append(val)
        if weights:
            total = sum(weights)
            if abs(total - 1.0) > 0.05:
                issues.append(
                    f"Collection weights sum to {total:.4f}, expected ~1.0 "
                    f"(tolerance 0.05)."
                )

        return issues

    def validate_or_warn(self) -> None:
        """Run validate() and log each issue as a warning."""
        for issue in self.validate():
            logger.warning("RareDisease config: %s", issue)


settings = RareDiseaseSettings()
settings.validate_or_warn()
