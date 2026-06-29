"""Imaging Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as rag-chat-pipeline/config/settings.py.
"""

from pathlib import Path
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


class ImagingSettings(BaseSettings):
    """Configuration for Imaging Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── RAG Pipeline (reuse existing) ──
    RAG_PIPELINE_ROOT: Path = Path(
        "/home/adam/projects/hcls-ai-factory/rag-chat-pipeline"
    )

    # ── Milvus ──
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # Collection names (10 imaging-specific)
    COLLECTION_LITERATURE: str = "imaging_literature"
    COLLECTION_TRIALS: str = "imaging_trials"
    COLLECTION_FINDINGS: str = "imaging_findings"
    COLLECTION_PROTOCOLS: str = "imaging_protocols"
    COLLECTION_DEVICES: str = "imaging_devices"
    COLLECTION_ANATOMY: str = "imaging_anatomy"
    COLLECTION_BENCHMARKS: str = "imaging_benchmarks"
    COLLECTION_GUIDELINES: str = "imaging_guidelines"
    COLLECTION_REPORT_TEMPLATES: str = "imaging_report_templates"
    COLLECTION_DATASETS: str = "imaging_datasets"
    # Read-only cross-agent collection
    COLLECTION_GENOMIC: str = "genomic_evidence"  # Existing collection

    # ── Embeddings ──
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIMENSION: int = 384
    EMBEDDING_BATCH_SIZE: int = 32

    # ── LLM ──
    LLM_PROVIDER: str = "anthropic"
    LLM_MODEL: str = "claude-sonnet-4-20250514"
    ANTHROPIC_API_KEY: Optional[str] = None

    # ── NIM Configuration (active in Phase 1) ──
    NIM_LLM_URL: str = "http://localhost:8520/v1"
    NIM_VISTA3D_URL: str = "http://localhost:8530"
    NIM_MAISI_URL: str = "http://localhost:8531"
    NIM_VILAM3_URL: str = "http://localhost:8532"
    NIM_MODE: str = "local"  # "local", "cloud", or "mock"
    NIM_ALLOW_MOCK_FALLBACK: bool = True
    NGC_API_KEY: Optional[str] = None

    # ── NVIDIA Cloud NIM Endpoints ──
    NVIDIA_API_KEY: Optional[str] = None
    NIM_CLOUD_URL: str = "https://integrate.api.nvidia.com/v1"
    NIM_CLOUD_LLM_MODEL: str = "meta/llama-3.1-8b-instruct"
    NIM_CLOUD_VLM_MODEL: str = "meta/llama-3.2-11b-vision-instruct"

    # ── RAG Search ──
    TOP_K_PER_COLLECTION: int = 5
    SCORE_THRESHOLD: float = 0.4

    # Collection search weights (must sum to ~1.0)
    WEIGHT_LITERATURE: float = 0.18
    WEIGHT_TRIALS: float = 0.12
    WEIGHT_FINDINGS: float = 0.15
    WEIGHT_PROTOCOLS: float = 0.08
    WEIGHT_DEVICES: float = 0.08
    WEIGHT_ANATOMY: float = 0.06
    WEIGHT_BENCHMARKS: float = 0.08
    WEIGHT_GUIDELINES: float = 0.10
    WEIGHT_REPORT_TEMPLATES: float = 0.05
    WEIGHT_DATASETS: float = 0.06
    WEIGHT_GENOMIC: float = 0.04

    # ── PubMed ──
    NCBI_API_KEY: Optional[str] = None  # Optional, increases rate limit
    PUBMED_MAX_RESULTS: int = 5000

    # ── ClinicalTrials.gov ──
    CT_GOV_BASE_URL: str = "https://clinicaltrials.gov/api/v2"

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8524

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8525

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

    # ── Orthanc DICOM Server ──
    ORTHANC_URL: str = "http://localhost:8042"
    ORTHANC_USERNAME: str = "admin"
    ORTHANC_PASSWORD: str = ""
    DICOM_AUTO_INGEST: bool = False
    DICOM_WATCH_INTERVAL: int = 5  # seconds between Orthanc /changes polls

    # ── Cross-Modal ──
    CROSS_MODAL_ENABLED: bool = False

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8524,http://localhost:8525"

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 10

    # ── Local NIM LLM Model Name ──
    NIM_LOCAL_LLM_MODEL: str = "meta/llama3-70b-instruct"

    # ── Legacy alias (kept for backward compat) ──
    DICOM_SERVER_URL: str = "http://localhost:8042"

    model_config = SettingsConfigDict(
        env_prefix="IMAGING_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )


settings = ImagingSettings()
