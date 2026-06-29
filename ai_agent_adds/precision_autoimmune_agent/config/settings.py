"""
Precision Autoimmune Agent — Settings

Pydantic BaseSettings with AUTO_ env prefix.
All values can be overridden via environment variables or .env file.
"""

from pathlib import Path
from typing import Dict

from pydantic import model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict


class AutoimmuneSettings(BaseSettings):
    """Central configuration for the Precision Autoimmune Agent."""

    model_config = SettingsConfigDict(
        env_prefix="AUTO_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # ── Paths ────────────────────────────────────────────────────────────
    PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

    @property
    def DATA_DIR(self) -> Path:
        return self.PROJECT_ROOT / "data"

    @property
    def CACHE_DIR(self) -> Path:
        return self.DATA_DIR / "cache"

    @property
    def REFERENCE_DIR(self) -> Path:
        return self.DATA_DIR / "reference"

    @property
    def DEMO_DATA_DIR(self) -> Path:
        return self.PROJECT_ROOT / "demo_data"

    # ── Milvus ───────────────────────────────────────────────────────────
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530

    # ── Collection names (14 autoimmune-domain collections) ──────────────
    COLL_CLINICAL_DOCUMENTS: str = "autoimmune_clinical_documents"
    COLL_PATIENT_LABS: str = "autoimmune_patient_labs"
    COLL_AUTOANTIBODY_PANELS: str = "autoimmune_autoantibody_panels"
    COLL_HLA_ASSOCIATIONS: str = "autoimmune_hla_associations"
    COLL_DISEASE_CRITERIA: str = "autoimmune_disease_criteria"
    COLL_DISEASE_ACTIVITY: str = "autoimmune_disease_activity"
    COLL_FLARE_PATTERNS: str = "autoimmune_flare_patterns"
    COLL_BIOLOGIC_THERAPIES: str = "autoimmune_biologic_therapies"
    COLL_PGX_RULES: str = "autoimmune_pgx_rules"
    COLL_CLINICAL_TRIALS: str = "autoimmune_clinical_trials"
    COLL_LITERATURE: str = "autoimmune_literature"
    COLL_PATIENT_TIMELINES: str = "autoimmune_patient_timelines"
    COLL_CROSS_DISEASE: str = "autoimmune_cross_disease"
    COLL_GENOMIC_EVIDENCE: str = "genomic_evidence"  # shared read-only

    # ── Embedding ────────────────────────────────────────────────────────
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIM: int = 384
    EMBEDDING_BATCH_SIZE: int = 32
    BGE_INSTRUCTION: str = "Represent this sentence for searching relevant passages: "

    # ── LLM ──────────────────────────────────────────────────────────────
    ANTHROPIC_API_KEY: str = ""
    LLM_MODEL: str = "claude-sonnet-4-20250514"
    LLM_MAX_TOKENS: int = 4096
    LLM_TEMPERATURE: float = 0.2

    # ── RAG parameters ───────────────────────────────────────────────────
    TOP_K_PER_COLLECTION: int = 5
    SCORE_THRESHOLD: float = 0.40
    MAX_EVIDENCE_ITEMS: int = 30
    CONVERSATION_MEMORY_SIZE: int = 3

    # ── Collection weights (sum ≈ 1.0) ───────────────────────────────────
    WEIGHT_CLINICAL_DOCUMENTS: float = 0.18
    WEIGHT_PATIENT_LABS: float = 0.14
    WEIGHT_AUTOANTIBODY_PANELS: float = 0.12
    WEIGHT_HLA_ASSOCIATIONS: float = 0.08
    WEIGHT_DISEASE_CRITERIA: float = 0.08
    WEIGHT_DISEASE_ACTIVITY: float = 0.07
    WEIGHT_FLARE_PATTERNS: float = 0.06
    WEIGHT_BIOLOGIC_THERAPIES: float = 0.06
    WEIGHT_PGX_RULES: float = 0.04
    WEIGHT_CLINICAL_TRIALS: float = 0.05
    WEIGHT_LITERATURE: float = 0.05
    WEIGHT_PATIENT_TIMELINES: float = 0.03
    WEIGHT_CROSS_DISEASE: float = 0.02
    WEIGHT_GENOMIC_EVIDENCE: float = 0.02

    # ── Ports ────────────────────────────────────────────────────────────
    STREAMLIT_PORT: int = 8531
    API_PORT: int = 8532

    # ── Authentication ───────────────────────────────────────────────────
    API_KEY: str = ""  # empty = no auth required
    CORS_ORIGINS: str = "http://localhost:8080,http://localhost:8531"
    MAX_REQUEST_SIZE_MB: int = 50  # PDF uploads

    # ── Document processing ──────────────────────────────────────────────
    MAX_CHUNK_SIZE: int = 2500
    CHUNK_OVERLAP: int = 200
    PDF_DPI: int = 200

    # ── Relevance thresholds ─────────────────────────────────────────────
    CITATION_HIGH: float = 0.80
    CITATION_MEDIUM: float = 0.60

    # ── Flare risk thresholds ──────────────────────────────────────────────
    FLARE_RISK_IMMINENT: float = 0.8
    FLARE_RISK_HIGH: float = 0.6
    FLARE_RISK_MODERATE: float = 0.4

    # ── Evidence display ───────────────────────────────────────────────────
    MAX_EVIDENCE_TEXT_LENGTH: int = 1500
    MAX_KNOWLEDGE_CONTEXT_ITEMS: int = 25

    # ── Streaming ──────────────────────────────────────────────────────────
    STREAMING_ENABLED: bool = True

    # ── Cross-Agent Integration ──────────────────────────────────────────
    ONCOLOGY_AGENT_URL: str = "http://localhost:8527"
    CARDIOLOGY_AGENT_URL: str = "http://localhost:8126"
    NEUROLOGY_AGENT_URL: str = "http://localhost:8528"
    PGX_AGENT_URL: str = "http://localhost:8107"
    BIOMARKER_AGENT_URL: str = "http://localhost:8529"
    TRIAL_AGENT_URL: str = "http://localhost:8538"
    CROSS_AGENT_TIMEOUT: int = 30

    # ── Timeouts ─────────────────────────────────────────────────────────
    REQUEST_TIMEOUT_SECONDS: int = 60
    MILVUS_TIMEOUT_SECONDS: int = 10
    LLM_MAX_RETRIES: int = 3

    # ── Logging ──────────────────────────────────────────────────────────
    LOG_LEVEL: str = "INFO"
    LOG_DIR: str = ""  # empty = <PROJECT_ROOT>/logs

    # ── Metrics ──────────────────────────────────────────────────────────
    METRICS_ENABLED: bool = True

    @model_validator(mode="after")
    def _validate_weights(self):
        """Warn if collection weights don't sum to ~1.0."""
        import logging
        total = (
            self.WEIGHT_CLINICAL_DOCUMENTS + self.WEIGHT_PATIENT_LABS +
            self.WEIGHT_AUTOANTIBODY_PANELS + self.WEIGHT_HLA_ASSOCIATIONS +
            self.WEIGHT_DISEASE_CRITERIA + self.WEIGHT_DISEASE_ACTIVITY +
            self.WEIGHT_FLARE_PATTERNS + self.WEIGHT_BIOLOGIC_THERAPIES +
            self.WEIGHT_PGX_RULES + self.WEIGHT_CLINICAL_TRIALS +
            self.WEIGHT_LITERATURE + self.WEIGHT_PATIENT_TIMELINES +
            self.WEIGHT_CROSS_DISEASE + self.WEIGHT_GENOMIC_EVIDENCE
        )
        if abs(total - 1.0) > 0.05:
            logging.getLogger(__name__).warning(
                f"Collection weights sum to {total:.3f}, expected ~1.0"
            )
        return self

    @property
    def collection_config(self) -> Dict[str, dict]:
        """Return collection name → {weight, label, name} mapping."""
        return {
            self.COLL_CLINICAL_DOCUMENTS: {
                "weight": self.WEIGHT_CLINICAL_DOCUMENTS,
                "label": "Clinical Document",
                "name": self.COLL_CLINICAL_DOCUMENTS,
            },
            self.COLL_PATIENT_LABS: {
                "weight": self.WEIGHT_PATIENT_LABS,
                "label": "Lab Result",
                "name": self.COLL_PATIENT_LABS,
            },
            self.COLL_AUTOANTIBODY_PANELS: {
                "weight": self.WEIGHT_AUTOANTIBODY_PANELS,
                "label": "Autoantibody",
                "name": self.COLL_AUTOANTIBODY_PANELS,
            },
            self.COLL_HLA_ASSOCIATIONS: {
                "weight": self.WEIGHT_HLA_ASSOCIATIONS,
                "label": "HLA Association",
                "name": self.COLL_HLA_ASSOCIATIONS,
            },
            self.COLL_DISEASE_CRITERIA: {
                "weight": self.WEIGHT_DISEASE_CRITERIA,
                "label": "Classification Criteria",
                "name": self.COLL_DISEASE_CRITERIA,
            },
            self.COLL_DISEASE_ACTIVITY: {
                "weight": self.WEIGHT_DISEASE_ACTIVITY,
                "label": "Disease Activity",
                "name": self.COLL_DISEASE_ACTIVITY,
            },
            self.COLL_FLARE_PATTERNS: {
                "weight": self.WEIGHT_FLARE_PATTERNS,
                "label": "Flare Pattern",
                "name": self.COLL_FLARE_PATTERNS,
            },
            self.COLL_BIOLOGIC_THERAPIES: {
                "weight": self.WEIGHT_BIOLOGIC_THERAPIES,
                "label": "Biologic Therapy",
                "name": self.COLL_BIOLOGIC_THERAPIES,
            },
            self.COLL_PGX_RULES: {
                "weight": self.WEIGHT_PGX_RULES,
                "label": "PGx Rule",
                "name": self.COLL_PGX_RULES,
            },
            self.COLL_CLINICAL_TRIALS: {
                "weight": self.WEIGHT_CLINICAL_TRIALS,
                "label": "Clinical Trial",
                "name": self.COLL_CLINICAL_TRIALS,
            },
            self.COLL_LITERATURE: {
                "weight": self.WEIGHT_LITERATURE,
                "label": "Literature",
                "name": self.COLL_LITERATURE,
                "year_field": "year",
            },
            self.COLL_PATIENT_TIMELINES: {
                "weight": self.WEIGHT_PATIENT_TIMELINES,
                "label": "Timeline",
                "name": self.COLL_PATIENT_TIMELINES,
            },
            self.COLL_CROSS_DISEASE: {
                "weight": self.WEIGHT_CROSS_DISEASE,
                "label": "Cross-Disease",
                "name": self.COLL_CROSS_DISEASE,
            },
            self.COLL_GENOMIC_EVIDENCE: {
                "weight": self.WEIGHT_GENOMIC_EVIDENCE,
                "label": "Genomic Evidence",
                "name": self.COLL_GENOMIC_EVIDENCE,
            },
        }

    @property
    def all_collection_names(self) -> list:
        """Return list of all collection names."""
        return list(self.collection_config.keys())


settings = AutoimmuneSettings()
