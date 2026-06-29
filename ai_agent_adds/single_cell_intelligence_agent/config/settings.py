"""Single-Cell Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as the Rare Disease agent.

Author: Adam Jones
Date: March 2026
"""

import logging
from pathlib import Path
from typing import List, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


class SingleCellSettings(BaseSettings):
    """Configuration for the Single-Cell Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (12 single-cell-specific collections)
    COLLECTION_CELL_TYPES: str = "sc_cell_types"
    COLLECTION_MARKERS: str = "sc_markers"
    COLLECTION_SPATIAL: str = "sc_spatial"
    COLLECTION_TME: str = "sc_tme"
    COLLECTION_DRUG_RESPONSE: str = "sc_drug_response"
    COLLECTION_LITERATURE: str = "sc_literature"
    COLLECTION_METHODS: str = "sc_methods"
    COLLECTION_DATASETS: str = "sc_datasets"
    COLLECTION_TRAJECTORIES: str = "sc_trajectories"
    COLLECTION_PATHWAYS: str = "sc_pathways"
    COLLECTION_CLINICAL: str = "sc_clinical"
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
    TOP_K_CELL_TYPES: int = 50
    TOP_K_MARKERS: int = 40
    TOP_K_SPATIAL: int = 30
    TOP_K_TME: int = 30
    TOP_K_DRUG_RESPONSE: int = 20
    TOP_K_LITERATURE: int = 20
    TOP_K_METHODS: int = 15
    TOP_K_DATASETS: int = 15
    TOP_K_TRAJECTORIES: int = 20
    TOP_K_PATHWAYS: int = 20
    TOP_K_CLINICAL: int = 15
    TOP_K_GENOMIC: int = 20

    # Collection search weights (must sum to ~1.0)
    WEIGHT_CELL_TYPES: float = 0.14
    WEIGHT_MARKERS: float = 0.12
    WEIGHT_SPATIAL: float = 0.10
    WEIGHT_TME: float = 0.10
    WEIGHT_DRUG_RESPONSE: float = 0.09
    WEIGHT_LITERATURE: float = 0.08
    WEIGHT_METHODS: float = 0.07
    WEIGHT_DATASETS: float = 0.06
    WEIGHT_TRAJECTORIES: float = 0.07
    WEIGHT_PATHWAYS: float = 0.07
    WEIGHT_CLINICAL: float = 0.07
    WEIGHT_GENOMIC: float = 0.03

    # ── External APIs ──
    CELLXGENE_API_URL: Optional[str] = None
    NCBI_API_KEY: Optional[str] = None

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8540

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8130

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
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8130,http://localhost:8540"

    # ── Cross-Agent Integration ──
    GENOMICS_AGENT_URL: str = "http://localhost:8527"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── GPU / RAPIDS ──
    GPU_MEMORY_LIMIT_GB: int = 120

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="SC_",
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
            logger.warning("SingleCell config: %s", issue)


settings = SingleCellSettings()
settings.validate_or_warn()
