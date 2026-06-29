"""Clinical Trial Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as the Cardiology agent.
"""

import logging
from pathlib import Path
from typing import List, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


class TrialSettings(BaseSettings):
    """Configuration for Clinical Trial Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (14 trial-specific collections)
    COLLECTION_PROTOCOLS: str = "trial_protocols"
    COLLECTION_ELIGIBILITY: str = "trial_eligibility"
    COLLECTION_ENDPOINTS: str = "trial_endpoints"
    COLLECTION_SITES: str = "trial_sites"
    COLLECTION_INVESTIGATORS: str = "trial_investigators"
    COLLECTION_RESULTS: str = "trial_results"
    COLLECTION_REGULATORY: str = "trial_regulatory"
    COLLECTION_LITERATURE: str = "trial_literature"
    COLLECTION_BIOMARKERS: str = "trial_biomarkers"
    COLLECTION_SAFETY: str = "trial_safety"
    COLLECTION_RWE: str = "trial_rwe"
    COLLECTION_ADAPTIVE: str = "trial_adaptive"
    COLLECTION_GUIDELINES: str = "trial_guidelines"
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
    WEIGHT_PROTOCOLS: float = 0.10
    WEIGHT_ELIGIBILITY: float = 0.09
    WEIGHT_ENDPOINTS: float = 0.08
    WEIGHT_SITES: float = 0.07
    WEIGHT_INVESTIGATORS: float = 0.05
    WEIGHT_RESULTS: float = 0.09
    WEIGHT_REGULATORY: float = 0.07
    WEIGHT_LITERATURE: float = 0.08
    WEIGHT_BIOMARKERS: float = 0.07
    WEIGHT_SAFETY: float = 0.08
    WEIGHT_RWE: float = 0.06
    WEIGHT_ADAPTIVE: float = 0.05
    WEIGHT_GUIDELINES: float = 0.08
    WEIGHT_GENOMIC: float = 0.03

    # ── External APIs ──
    CLINICALTRIALS_API_KEY: Optional[str] = None
    NCBI_API_KEY: Optional[str] = None

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8538

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8128

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
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8538,http://localhost:8128"

    # ── Cross-Agent Integration ──
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    PGX_AGENT_URL: str = "http://localhost:8107"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    RARE_DISEASE_AGENT_URL: str = "http://localhost:8134"
    NEUROLOGY_AGENT_URL: str = "http://localhost:8528"
    SINGLE_CELL_AGENT_URL: str = "http://localhost:8540"
    IMAGING_AGENT_URL: str = "http://localhost:8524"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="TRIAL_",
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

        return issues

    def validate_or_warn(self) -> None:
        """Run validate() and log each issue as a warning."""
        for issue in self.validate():
            logger.warning("Trial config: %s", issue)


settings = TrialSettings()
settings.validate_or_warn()
