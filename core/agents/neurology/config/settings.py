"""Neurology Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as the Rare Disease Diagnostic Agent
and Clinical Trial Intelligence Agent.

Author: Adam Jones
Date: March 2026
"""

import logging
from pathlib import Path
from typing import List, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


class NeuroSettings(BaseSettings):
    """Configuration for the Neurology Intelligence Agent."""

    # -- Paths --
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # -- Milvus --
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (14 neurology-specific collections)
    COLLECTION_LITERATURE: str = "neuro_literature"
    COLLECTION_TRIALS: str = "neuro_trials"
    COLLECTION_IMAGING: str = "neuro_imaging"
    COLLECTION_ELECTROPHYSIOLOGY: str = "neuro_electrophysiology"
    COLLECTION_DEGENERATIVE: str = "neuro_degenerative"
    COLLECTION_CEREBROVASCULAR: str = "neuro_cerebrovascular"
    COLLECTION_EPILEPSY: str = "neuro_epilepsy"
    COLLECTION_ONCOLOGY: str = "neuro_oncology"
    COLLECTION_MS: str = "neuro_ms"
    COLLECTION_MOVEMENT: str = "neuro_movement"
    COLLECTION_HEADACHE: str = "neuro_headache"
    COLLECTION_NEUROMUSCULAR: str = "neuro_neuromuscular"
    COLLECTION_GUIDELINES: str = "neuro_guidelines"
    COLLECTION_GENOMIC: str = "genomic_evidence"

    # -- Embeddings --
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIMENSION: int = 384
    EMBEDDING_BATCH_SIZE: int = 32

    # -- LLM --
    LLM_PROVIDER: str = "anthropic"
    LLM_MODEL: str = "claude-sonnet-4-6"
    ANTHROPIC_API_KEY: Optional[str] = None

    # -- RAG Search --
    SCORE_THRESHOLD: float = 0.4

    # Per-collection TOP_K defaults
    TOP_K_LITERATURE: int = 20
    TOP_K_TRIALS: int = 15
    TOP_K_IMAGING: int = 30
    TOP_K_ELECTROPHYSIOLOGY: int = 25
    TOP_K_DEGENERATIVE: int = 30
    TOP_K_CEREBROVASCULAR: int = 30
    TOP_K_EPILEPSY: int = 25
    TOP_K_ONCOLOGY: int = 20
    TOP_K_MS: int = 25
    TOP_K_MOVEMENT: int = 25
    TOP_K_HEADACHE: int = 20
    TOP_K_NEUROMUSCULAR: int = 20
    TOP_K_GUIDELINES: int = 15
    TOP_K_GENOMIC: int = 20

    # Collection search weights (must sum to ~1.0)
    WEIGHT_LITERATURE: float = 0.08
    WEIGHT_TRIALS: float = 0.06
    WEIGHT_IMAGING: float = 0.09
    WEIGHT_ELECTROPHYSIOLOGY: float = 0.07
    WEIGHT_DEGENERATIVE: float = 0.09
    WEIGHT_CEREBROVASCULAR: float = 0.09
    WEIGHT_EPILEPSY: float = 0.08
    WEIGHT_ONCOLOGY: float = 0.06
    WEIGHT_MS: float = 0.07
    WEIGHT_MOVEMENT: float = 0.07
    WEIGHT_HEADACHE: float = 0.06
    WEIGHT_NEUROMUSCULAR: float = 0.06
    WEIGHT_GUIDELINES: float = 0.07
    WEIGHT_GENOMIC: float = 0.05

    # -- External APIs --
    CLINICALTRIALS_API_KEY: Optional[str] = None
    NCBI_API_KEY: Optional[str] = None

    # -- API Server --
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8528

    # -- Streamlit --
    STREAMLIT_PORT: int = 8534

    # -- Prometheus Metrics --
    METRICS_ENABLED: bool = True

    # -- Scheduler --
    INGEST_SCHEDULE_HOURS: int = 24
    INGEST_ENABLED: bool = False

    # -- Conversation Memory --
    MAX_CONVERSATION_CONTEXT: int = 3

    # -- Citation Scoring --
    CITATION_HIGH_THRESHOLD: float = 0.75
    CITATION_MEDIUM_THRESHOLD: float = 0.60

    # -- Authentication --
    API_KEY: str = ""  # Empty = no auth required

    # -- CORS --
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8528,http://localhost:8529"

    # -- Cross-Agent Integration --
    GENOMICS_AGENT_URL: str = "http://localhost:8527"
    IMAGING_AGENT_URL: str = "http://localhost:8524"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    RARE_DISEASE_AGENT_URL: str = "http://localhost:8134"
    CROSS_AGENT_TIMEOUT: int = 30

    # -- Request Limits --
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="NEURO_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )

    # -- Startup Validation --

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
            logger.warning("Neuro config: %s", issue)


settings = NeuroSettings()
settings.validate_or_warn()
