"""CAR-T Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as rag-chat-pipeline/config/settings.py.
"""

import os
from pathlib import Path
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


class CARTSettings(BaseSettings):
    """Configuration for CAR-T Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── RAG Pipeline (reuse existing) ──
    RAG_PIPELINE_ROOT: Path = Path(
        os.environ.get("CART_RAG_PIPELINE_ROOT", "/app/rag-chat-pipeline")
    )

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names
    COLLECTION_LITERATURE: str = "cart_literature"
    COLLECTION_TRIALS: str = "cart_trials"
    COLLECTION_CONSTRUCTS: str = "cart_constructs"
    COLLECTION_ASSAYS: str = "cart_assays"
    COLLECTION_MANUFACTURING: str = "cart_manufacturing"
    COLLECTION_GENOMIC: str = "genomic_evidence"  # Existing collection
    COLLECTION_SAFETY: str = "cart_safety"
    COLLECTION_BIOMARKERS: str = "cart_biomarkers"
    COLLECTION_REGULATORY: str = "cart_regulatory"
    COLLECTION_SEQUENCES: str = "cart_sequences"
    COLLECTION_REALWORLD: str = "cart_realworld"

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
    WEIGHT_LITERATURE: float = 0.20
    WEIGHT_TRIALS: float = 0.16
    WEIGHT_CONSTRUCTS: float = 0.10
    WEIGHT_ASSAYS: float = 0.09
    WEIGHT_MANUFACTURING: float = 0.07
    WEIGHT_SAFETY: float = 0.08
    WEIGHT_BIOMARKERS: float = 0.08
    WEIGHT_REGULATORY: float = 0.06
    WEIGHT_SEQUENCES: float = 0.06
    WEIGHT_REALWORLD: float = 0.07
    WEIGHT_GENOMIC: float = 0.04

    # ── PubMed ──
    NCBI_API_KEY: Optional[str] = None  # Optional, increases rate limit
    PUBMED_MAX_RESULTS: int = 5000

    # ── ClinicalTrials.gov ──
    CT_GOV_BASE_URL: str = "https://clinicaltrials.gov/api/v2"

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8522
    API_KEY: str = ""

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8521

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
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    SINGLE_CELL_AGENT_URL: str = "http://localhost:8540"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8521,http://localhost:8522"

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    model_config = SettingsConfigDict(
        env_prefix="CART_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )


settings = CARTSettings()
