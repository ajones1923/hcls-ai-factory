"""Biomarker Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as cart_intelligence_agent/config/settings.py.
"""

import os
from pathlib import Path
from typing import Optional

from pydantic import model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class PrecisionBiomarkerSettings(BaseSettings):
    """Configuration for Biomarker Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── RAG Pipeline (reuse existing) ──
    RAG_PIPELINE_ROOT: Path = Path(
        os.environ.get("BIOMARKER_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")
    )

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (13 biomarker-specific + 1 read-only genomic)
    COLLECTION_BIOMARKER_REF: str = "biomarker_reference"
    COLLECTION_GENETIC_VARIANTS: str = "biomarker_genetic_variants"
    COLLECTION_PGX_RULES: str = "biomarker_pgx_rules"
    COLLECTION_DISEASE_TRAJECTORIES: str = "biomarker_disease_trajectories"
    COLLECTION_CLINICAL_EVIDENCE: str = "biomarker_clinical_evidence"
    COLLECTION_NUTRITION: str = "biomarker_nutrition"
    COLLECTION_DRUG_INTERACTIONS: str = "biomarker_drug_interactions"
    COLLECTION_AGING_MARKERS: str = "biomarker_aging_markers"
    COLLECTION_GENOTYPE_ADJUSTMENTS: str = "biomarker_genotype_adjustments"
    COLLECTION_MONITORING: str = "biomarker_monitoring"
    COLLECTION_CRITICAL_VALUES: str = "biomarker_critical_values"
    COLLECTION_DISCORDANCE_RULES: str = "biomarker_discordance_rules"
    COLLECTION_AJ_CARRIER_SCREENING: str = "biomarker_aj_carrier_screening"
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
    WEIGHT_BIOMARKER_REF: float = 0.12
    WEIGHT_GENETIC_VARIANTS: float = 0.11
    WEIGHT_PGX_RULES: float = 0.10
    WEIGHT_DISEASE_TRAJECTORIES: float = 0.10
    WEIGHT_CLINICAL_EVIDENCE: float = 0.09
    WEIGHT_NUTRITION: float = 0.05
    WEIGHT_DRUG_INTERACTIONS: float = 0.07
    WEIGHT_AGING_MARKERS: float = 0.07
    WEIGHT_GENOTYPE_ADJUSTMENTS: float = 0.05
    WEIGHT_MONITORING: float = 0.05
    WEIGHT_CRITICAL_VALUES: float = 0.04
    WEIGHT_DISCORDANCE_RULES: float = 0.04
    WEIGHT_AJ_CARRIER_SCREENING: float = 0.03
    WEIGHT_GENOMIC_EVIDENCE: float = 0.08

    # ── Timeouts ──
    REQUEST_TIMEOUT_SECONDS: int = 60
    MILVUS_TIMEOUT_SECONDS: int = 10
    LLM_MAX_RETRIES: int = 3

    # ── Port validation ──
    @model_validator(mode="after")
    def _validate_settings(self) -> "PrecisionBiomarkerSettings":
        weights = [
            self.WEIGHT_BIOMARKER_REF, self.WEIGHT_GENETIC_VARIANTS,
            self.WEIGHT_PGX_RULES, self.WEIGHT_DISEASE_TRAJECTORIES,
            self.WEIGHT_CLINICAL_EVIDENCE, self.WEIGHT_NUTRITION,
            self.WEIGHT_DRUG_INTERACTIONS, self.WEIGHT_AGING_MARKERS,
            self.WEIGHT_GENOTYPE_ADJUSTMENTS, self.WEIGHT_MONITORING,
            self.WEIGHT_CRITICAL_VALUES, self.WEIGHT_DISCORDANCE_RULES,
            self.WEIGHT_AJ_CARRIER_SCREENING,
            self.WEIGHT_GENOMIC_EVIDENCE,
        ]
        total = sum(weights)
        if abs(total - 1.0) > 0.05:
            from loguru import logger
            logger.warning(f"Collection search weights sum to {total:.2f}, expected ~1.0")
        return self

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8529

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8533

    # ── Prometheus Metrics ──
    METRICS_ENABLED: bool = True

    # ── Conversation Memory ──
    MAX_CONVERSATION_CONTEXT: int = 3

    # ── Citation Scoring ──
    CITATION_HIGH_THRESHOLD: float = 0.75
    CITATION_MEDIUM_THRESHOLD: float = 0.60

    # ── Cross-Agent Integration ──
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    CART_AGENT_URL: str = "http://localhost:8522"
    PGX_AGENT_URL: str = "http://localhost:8107"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8528,http://localhost:8529"

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    # ── Authentication ──
    # WARNING: Empty API_KEY disables authentication entirely. In production
    # deployments handling patient health data (PHI/PII), always set
    # BIOMARKER_API_KEY to a strong secret to enforce request authentication
    # and comply with HIPAA / data-security requirements.
    API_KEY: str = ""  # Optional API key; empty = no auth required

    model_config = SettingsConfigDict(
        env_prefix="BIOMARKER_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )


settings = PrecisionBiomarkerSettings()
