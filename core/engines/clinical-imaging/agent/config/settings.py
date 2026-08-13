"""Imaging Intelligence Agent configuration.

Follows the same Pydantic BaseSettings pattern as core/engines/precision-intelligence/config/settings.py.
"""

import os
from enum import Enum
from pathlib import Path
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


class DeploymentTier(str, Enum):
    """Deployment tier controls which NIM backends and storage are available.

    - community:   20 free NVIDIA technologies (Apache 2.0 / BSD / MIT / Open Model)
    - enterprise:  28 total (adds AI Enterprise NIMs + the AI platform)
    - research:    adds noncommercial models (NV-Reason-CXR-3B, NV-Segment-CTMR)
    """

    COMMUNITY = "community"
    ENTERPRISE = "enterprise"
    RESEARCH = "research"


class ImagingSettings(BaseSettings):
    """Configuration for Imaging Intelligence Agent."""

    # ── Paths ──
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"
    REFERENCE_DIR: Path = DATA_DIR / "reference"

    # ── RAG Pipeline (reuse existing) ──
    RAG_PIPELINE_ROOT: Path = Path(
        os.environ.get("IMAGING_RAG_PIPELINE_ROOT", "/app/core/engines/precision-intelligence")
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

    # ── Deployment Tier ──
    DEPLOYMENT_TIER: str = "community"  # "community", "enterprise", "research"

    # ── NIM Configuration (active in Phase 1) ──
    NIM_LLM_URL: str = "http://localhost:8520/v1"
    NIM_VISTA3D_URL: str = "http://localhost:8530"
    NIM_MAISI_URL: str = "http://localhost:8531"
    NIM_VILAM3_URL: str = "http://localhost:8532"
    NIM_NV_SEGMENT_CT_URL: str = "http://localhost:8534"

    # ── NV-Generate ──
    NIM_NV_GENERATE_CT_URL: str = "http://localhost:8539"
    NIM_NV_GENERATE_MR_URL: str = "http://localhost:8540"

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
    WEIGHT_TRIALS: float = 0.10
    WEIGHT_FINDINGS: float = 0.13
    WEIGHT_PROTOCOLS: float = 0.06
    WEIGHT_DEVICES: float = 0.08
    WEIGHT_ANATOMY: float = 0.06
    WEIGHT_BENCHMARKS: float = 0.08
    WEIGHT_GUIDELINES: float = 0.10
    WEIGHT_REPORT_TEMPLATES: float = 0.05
    WEIGHT_DATASETS: float = 0.02
    WEIGHT_GENOMIC: float = 0.04

    # ── Analytics ──
    ANALYTICS_ENABLED: bool = True
    ANALYTICS_DEMO_STUDIES: int = 500

    # ── Real-Time Streaming ──
    STREAMING_ENABLED: bool = False  # Off by default
    STREAMING_DEFAULT_FPS: int = 30
    STREAMING_OUTPUT_DIR: str = str(DATA_DIR / "streaming_output")

    # ── GPU Preprocessing ──
    DALI_ENABLED: bool = True  # Auto-detected at runtime
    CUCIM_ENABLED: bool = True  # Auto-detected at runtime
    DEFAULT_CT_WINDOW: str = "soft_tissue"
    DEFAULT_TARGET_SPACING: float = 1.0  # mm isotropic

    # ── Radiomics ──
    COLLECTION_RADIOMICS: str = "imaging_radiomics"
    WEIGHT_RADIOMICS: float = 0.04
    RADIOMICS_FEATURE_CLASSES: str = "firstorder,shape,glcm"
    RADIOMICS_FILTERS: str = ""
    RADIOMICS_GPU_ENABLED: bool = True

    # ── Radiology Reports ──
    COLLECTION_REPORTS: str = "imaging_reports"
    WEIGHT_REPORTS: float = 0.06
    REPORT_EXTRACT_WITH_LLM: bool = False

    # ── PubMed ──
    NCBI_API_KEY: Optional[str] = None  # Optional, increases rate limit
    PUBMED_MAX_RESULTS: int = 5000

    # ── ClinicalTrials.gov ──
    CT_GOV_BASE_URL: str = "https://clinicaltrials.gov/api/v2"

    # ── API Server ──
    API_HOST: str = "0.0.0.0"
    API_PORT: int = 8524
    API_KEY: str = ""

    # ── Streamlit ──
    STREAMLIT_PORT: int = 8525

    # ── Prometheus Metrics ──
    METRICS_ENABLED: bool = True

    # ── Scheduler ──
    INGEST_SCHEDULE_HOURS: int = 168  # PubMed: Weekly (7 * 24)
    INGEST_TRIALS_SCHEDULE_HOURS: int = 720  # ClinicalTrials.gov: Monthly (30 * 24)
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

    # ── OHIF Viewer ──
    OHIF_URL: str = "http://localhost:8526"

    # ── MONAI Label ──
    MONAI_LABEL_URL: str = "http://localhost:8527"
    MONAI_LABEL_ENABLED: bool = True
    MONAI_LABEL_APP_DIR: str = str(DATA_DIR / "monai_label_apps")
    MONAI_LABEL_STUDIES_DIR: str = str(DATA_DIR / "monai_label_studies")

    # ── Internal API URL (used by Streamlit to reach the FastAPI server) ──
    API_BASE_URL: str = "http://localhost:8524"

    # ── Preview Generation ──
    PREVIEW_CACHE_DIR: str = str(DATA_DIR / "cache" / "previews")
    PREVIEW_DEFAULT_FPS: int = 8
    PREVIEW_DEFAULT_FORMAT: str = "mp4"
    PREVIEW_MAX_FRAMES: int = 200

    # ── Guardrails ──
    GUARDRAILS_ENABLED: bool = True
    GUARDRAILS_CONFIG_DIR: str = str(PROJECT_ROOT / "config" / "guardrails")
    GUARDRAILS_MODE: str = "lite"  # "full" (NeMo Guardrails) or "lite" (regex-based)

    # ── Cross-Agent Integration ──
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    NEUROLOGY_AGENT_URL: str = "http://localhost:8528"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── Agentic Reasoning ──
    AGENTIC_ENABLED: bool = True
    AGENTIC_MAX_REFINEMENT_ROUNDS: int = 2
    AGENTIC_CONFIDENCE_THRESHOLD: float = 0.7
    AGENTIC_MIN_EVIDENCE_COUNT: int = 3

    # ── Cross-Modal ──
    CROSS_MODAL_ENABLED: bool = False

    # ── CORS ──
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8524,http://localhost:8525,http://localhost:8550,http://192.168.68.107:8550,http://192.168.68.107:8525,http://192.168.68.107:3001,http://192.168.68.107:8524,http://192.168.68.107:3000,http://192.168.68.107:3002,http://192.168.68.107:3003,http://192.168.68.107:5173,http://192.168.68.107:5174,http://localhost:3000,http://localhost:3001,http://localhost:3002,http://localhost:3003,http://localhost:5173,http://localhost:5174"

    # ── Request Limits ──
    MAX_REQUEST_SIZE_MB: int = 500  # Increased for CT series uploads

    # ── Local NIM LLM Model Name ──
    NIM_LOCAL_LLM_MODEL: str = "meta/llama3-70b-instruct"

    # ── Dynamo ──
    DYNAMO_ENABLED: bool = False
    DYNAMO_PREFILL_BATCH: int = 4
    DYNAMO_DECODE_BATCH: int = 8
    DYNAMO_KV_CACHE_GB: float = 16.0

    # ── Nemotron Nano ──
    NIM_NEMOTRON_NANO_URL: str = "http://localhost:8538"
    NEMOTRON_NANO_ENABLED: bool = False
    NEMOTRON_NANO_ROUTE_ROUTINE: bool = True  # Auto-route routine queries to Nano

    # ── Triton Inference Server ──
    TRITON_URL: str = "localhost:8000"
    TRITON_ENABLED: bool = False  # Off by default, opt-in
    TRITON_GRPC_URL: str = "localhost:8001"

    # ── TensorRT ──
    TENSORRT_ENABLED: bool = False  # Off by default, opt-in
    TENSORRT_PRECISION: str = "fp16"  # "fp32", "fp16", "int8"
    TENSORRT_CACHE_DIR: str = str(DATA_DIR / "cache" / "tensorrt")

    # ── Legacy alias (kept for backward compat) ──
    DICOM_SERVER_URL: str = "http://localhost:8042"

    # ── Tier-aware properties ──
    @property
    def segmentation_backend(self) -> str:
        """Return the segmentation NIM backend for the active deployment tier."""
        if self.DEPLOYMENT_TIER == "community":
            return "nv-segment-ct"
        elif self.DEPLOYMENT_TIER == "enterprise":
            return "vista3d-nim"
        return "nv-segment-ctmr"

    @property
    def cxr_backend(self) -> str:
        """Return the chest X-ray classification backend for the active tier."""
        if self.DEPLOYMENT_TIER in ("community", "enterprise"):
            return "densenet121"
        return "nv-reason-cxr-3b"

    @property
    def embedding_backend(self) -> str:
        """Return the embedding model backend for the active tier."""
        if self.DEPLOYMENT_TIER == "enterprise":
            return "nemo-retriever"
        return "bge-small-en-v1.5"

    @property
    def storage_backend(self) -> str:
        """Return the storage backend for the active tier."""
        if self.DEPLOYMENT_TIER == "enterprise":
            return "vast-aios"
        return "local"

    model_config = SettingsConfigDict(
        env_prefix="IMAGING_",
        case_sensitive=False,
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )


settings = ImagingSettings()
