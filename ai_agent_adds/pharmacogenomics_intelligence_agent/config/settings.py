"""Pharmacogenomics Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as rag-chat-pipeline/config/settings.py.
"""

import logging
import os
from pathlib import Path
from typing import List, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict

logger = logging.getLogger(__name__)


class PGxSettings(BaseSettings):
    """Configuration for Pharmacogenomics Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── RAG Pipeline (reuse existing) ──
    RAG_PIPELINE_ROOT: Path = Path(
        os.environ.get("PGX_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")
    )

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (14 PGx-specific + 1 shared genomic_evidence)
    COLLECTION_GENE_REFERENCE: str = "pgx_gene_reference"
    COLLECTION_DRUG_GUIDELINES: str = "pgx_drug_guidelines"
    COLLECTION_DRUG_INTERACTIONS: str = "pgx_drug_interactions"
    COLLECTION_HLA_HYPERSENSITIVITY: str = "pgx_hla_hypersensitivity"
    COLLECTION_PHENOCONVERSION: str = "pgx_phenoconversion"
    COLLECTION_DOSING_ALGORITHMS: str = "pgx_dosing_algorithms"
    COLLECTION_CLINICAL_EVIDENCE: str = "pgx_clinical_evidence"
    COLLECTION_POPULATION_DATA: str = "pgx_population_data"
    COLLECTION_CLINICAL_TRIALS: str = "pgx_clinical_trials"
    COLLECTION_FDA_LABELS: str = "pgx_fda_labels"
    COLLECTION_DRUG_ALTERNATIVES: str = "pgx_drug_alternatives"
    COLLECTION_PATIENT_PROFILES: str = "pgx_patient_profiles"
    COLLECTION_IMPLEMENTATION: str = "pgx_implementation"
    COLLECTION_EDUCATION: str = "pgx_education"
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
    WEIGHT_GENE_REFERENCE: float = 0.10
    WEIGHT_DRUG_GUIDELINES: float = 0.14
    WEIGHT_DRUG_INTERACTIONS: float = 0.12
    WEIGHT_HLA_HYPERSENSITIVITY: float = 0.10
    WEIGHT_PHENOCONVERSION: float = 0.08
    WEIGHT_DOSING_ALGORITHMS: float = 0.07
    WEIGHT_CLINICAL_EVIDENCE: float = 0.08
    WEIGHT_POPULATION_DATA: float = 0.06
    WEIGHT_CLINICAL_TRIALS: float = 0.04
    WEIGHT_FDA_LABELS: float = 0.06
    WEIGHT_DRUG_ALTERNATIVES: float = 0.05
    WEIGHT_PATIENT_PROFILES: float = 0.03
    WEIGHT_IMPLEMENTATION: float = 0.02
    WEIGHT_EDUCATION: float = 0.02
    WEIGHT_GENOMIC: float = 0.03

    # ── PubMed ──
    NCBI_API_KEY: Optional[str] = None  # Optional, increases rate limit
    PUBMED_MAX_RESULTS: int = 5000

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8107
    API_KEY: str = ""

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8507

    # ── Prometheus Metrics ──
    METRICS_ENABLED: bool = True

    # ── Scheduler ──
    INGEST_SCHEDULE_HOURS: int = 168  # Weekly (7 * 24)
    INGEST_ENABLED: bool = False

    # ── Conversation Memory ──
    MAX_CONVERSATION_CONTEXT: int = 3  # Number of prior exchanges to inject

    # ── Citation Scoring ──
    CITATION_HIGH_THRESHOLD: float = 0.75
    CITATION_MEDIUM_THRESHOLD: float = 0.60

    # ── Cross-Agent Integration ──
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    NEUROLOGY_AGENT_URL: str = "http://localhost:8528"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8107,http://localhost:8507"

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="PGX_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )

    # ── Startup Validation ──

    def validate(self) -> List[str]:
        """Return a list of configuration warnings/errors (never raises)."""
        issues: List[str] = []

        # 1. Milvus connection info
        if not self.MILVUS_HOST or not self.MILVUS_HOST.strip():
            issues.append("MILVUS_HOST is empty — Milvus connections will fail.")
        if not (1 <= self.MILVUS_PORT <= 65535):
            issues.append(
                f"MILVUS_PORT={self.MILVUS_PORT} is outside valid range (1-65535)."
            )

        # 2. API key availability
        if not self.ANTHROPIC_API_KEY:
            issues.append(
                "ANTHROPIC_API_KEY is not set — LLM features disabled, "
                "search-only mode available."
            )

        # 3. Embedding model
        if not self.EMBEDDING_MODEL or not self.EMBEDDING_MODEL.strip():
            issues.append("EMBEDDING_MODEL is empty — embedding pipeline will fail.")

        # 4. Port conflicts
        for name, port in [("API_PORT", self.API_PORT), ("STREAMLIT_PORT", self.STREAMLIT_PORT)]:
            if not (1024 <= port <= 65535):
                issues.append(
                    f"{name}={port} is outside valid range (1024-65535)."
                )
        if self.API_PORT == self.STREAMLIT_PORT:
            issues.append(
                f"API_PORT and STREAMLIT_PORT are both {self.API_PORT} — port conflict."
            )

        # 5. Collection weights
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

        # 6. RAG pipeline root
        if self.RAG_PIPELINE_ROOT and not self.RAG_PIPELINE_ROOT.is_dir():
            issues.append(
                f"RAG_PIPELINE_ROOT={self.RAG_PIPELINE_ROOT} does not exist or "
                f"is not a directory."
            )

        return issues

    def validate_or_warn(self) -> None:
        """Run validate() and log each issue as a warning."""
        for issue in self.validate():
            logger.warning("PGx config: %s", issue)


settings = PGxSettings()
settings.validate_or_warn()
