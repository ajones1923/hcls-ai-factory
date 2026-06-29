"""
Oncology Intelligence Agent - Configuration Settings
======================================================
Pydantic BaseSettings for the Oncology Intelligence RAG agent.
All values can be overridden via environment variables with the ONCO_ prefix.
"""

from pathlib import Path
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


class OncoSettings(BaseSettings):
    """Central configuration for the Oncology Intelligence Agent."""

    model_config = SettingsConfigDict(
        env_prefix="ONCO_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
    )

    # ── Paths ────────────────────────────────────────────────────────────
    PROJECT_ROOT: Path = Path(__file__).resolve().parents[2]
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = PROJECT_ROOT / "cache"
    REFERENCE_DIR: Path = PROJECT_ROOT / "reference"
    RAG_PIPELINE_ROOT: Path = Path(__file__).resolve().parents[4] / "rag-chat-pipeline"

    # ── Milvus ───────────────────────────────────────────────────────────
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # ── Collection Names ─────────────────────────────────────────────────
    COLLECTION_LITERATURE: str = "onco_literature"
    COLLECTION_TRIALS: str = "onco_trials"
    COLLECTION_VARIANTS: str = "onco_variants"
    COLLECTION_BIOMARKERS: str = "onco_biomarkers"
    COLLECTION_THERAPIES: str = "onco_therapies"
    COLLECTION_PATHWAYS: str = "onco_pathways"
    COLLECTION_GUIDELINES: str = "onco_guidelines"
    COLLECTION_RESISTANCE: str = "onco_resistance"
    COLLECTION_OUTCOMES: str = "onco_outcomes"
    COLLECTION_CASES: str = "onco_cases"
    COLLECTION_GENOMIC: str = "genomic_evidence"

    # ── Embeddings ───────────────────────────────────────────────────────
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIM: int = 384
    EMBEDDING_BATCH_SIZE: int = 32

    # ── LLM ──────────────────────────────────────────────────────────────
    LLM_PROVIDER: str = "anthropic"
    LLM_MODEL: str = "claude-sonnet-4-20250514"
    ANTHROPIC_API_KEY: Optional[str] = None

    # ── RAG Search ───────────────────────────────────────────────────────
    TOP_K: int = 5
    SCORE_THRESHOLD: float = 0.4

    # ── Collection Weights (sum ~1.0) ────────────────────────────────────
    WEIGHT_VARIANTS: float = 0.18
    WEIGHT_LITERATURE: float = 0.16
    WEIGHT_THERAPIES: float = 0.14
    WEIGHT_GUIDELINES: float = 0.12
    WEIGHT_TRIALS: float = 0.10
    WEIGHT_BIOMARKERS: float = 0.08
    WEIGHT_RESISTANCE: float = 0.07
    WEIGHT_PATHWAYS: float = 0.06
    WEIGHT_OUTCOMES: float = 0.04
    WEIGHT_CASES: float = 0.02
    WEIGHT_GENOMIC: float = 0.03

    # ── PubMed ───────────────────────────────────────────────────────────
    NCBI_API_KEY: Optional[str] = None
    PUBMED_MAX_RESULTS: int = 5000

    # ── ClinicalTrials.gov ───────────────────────────────────────────────
    CT_GOV_BASE_URL: str = "https://clinicaltrials.gov/api/v2"

    # ── CIViC ────────────────────────────────────────────────────────────
    CIVIC_BASE_URL: str = "https://civicdb.org/api"

    # ── API ──────────────────────────────────────────────────────────────
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8527
    API_KEY: str = ""

    # ── Streamlit ────────────────────────────────────────────────────────
    STREAMLIT_PORT: int = 8526

    # ── Metrics ──────────────────────────────────────────────────────────
    METRICS_ENABLED: bool = True

    # ── Scheduler ────────────────────────────────────────────────────────
    SCHEDULER_INTERVAL: str = "168h"

    # ── Conversation Memory ──────────────────────────────────────────────
    CONVERSATION_MEMORY_DEPTH: int = 3

    # ── Citation Thresholds ──────────────────────────────────────────────
    CITATION_STRONG_THRESHOLD: float = 0.75
    CITATION_MODERATE_THRESHOLD: float = 0.60

    # ── Cross-Modal ──────────────────────────────────────────────────────
    CROSS_MODAL_ENABLED: bool = True
    CROSS_MODAL_THRESHOLD: float = 0.40
    GENOMIC_TOP_K: int = 5
    IMAGING_TOP_K: int = 5

    # ── Cross-Agent URLs ─────────────────────────────────────────────
    CART_AGENT_URL: str = "http://localhost:8522"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    NEUROLOGY_AGENT_URL: str = "http://localhost:8528"
    PGX_AGENT_URL: str = "http://localhost:8107"
    IMAGING_AGENT_URL: str = "http://localhost:8524"
    SINGLE_CELL_AGENT_URL: str = "http://localhost:8540"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── Trial Matching ────────────────────────────────────────────────────
    TRIAL_WEIGHT_BIOMARKER: float = 0.40
    TRIAL_WEIGHT_SEMANTIC: float = 0.25
    TRIAL_WEIGHT_PHASE: float = 0.20
    TRIAL_WEIGHT_STATUS: float = 0.15

    # ── Agent ─────────────────────────────────────────────────────────────
    MIN_SUFFICIENT_HITS: int = 3
    MIN_COLLECTIONS_FOR_SUFFICIENT: int = 2
    MIN_SIMILARITY_SCORE: float = 0.30

    # ── API ───────────────────────────────────────────────────────────────
    API_BASE_URL: str = "http://localhost:8527"
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8526,http://localhost:8527"
    MAX_REQUEST_SIZE_MB: int = 10

    # ── PDF Branding ──────────────────────────────────────────────────────
    PDF_BRAND_COLOR_R: int = 118
    PDF_BRAND_COLOR_G: int = 185
    PDF_BRAND_COLOR_B: int = 0


# ── Singleton ────────────────────────────────────────────────────────────
settings = OncoSettings()
