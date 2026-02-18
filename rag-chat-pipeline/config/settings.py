"""
RAG Chat Pipeline Configuration
"""
import os
from pathlib import Path
from typing import Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore"
    )

    # Project Paths
    PROJECT_ROOT: Path = Path(__file__).parent.parent
    DATA_DIR: Path = PROJECT_ROOT / "data"
    CACHE_DIR: Path = DATA_DIR / "cache"

    # HCLS AI Factory Root (parent of all pipelines)
    HCLS_ROOT: Path = PROJECT_ROOT.parent

    # Genomics Pipeline Integration
    # Can be overridden via GENOMICS_PIPELINE_DIR environment variable
    GENOMICS_PIPELINE_DIR: Path = Path(
        os.environ.get("GENOMICS_PIPELINE_DIR", str(PROJECT_ROOT.parent / "genomics-pipeline"))
    )
    VCF_INPUT_PATH: Path = Path(
        os.environ.get("VCF_INPUT_PATH", str(GENOMICS_PIPELINE_DIR / "data/output/HG002.genome.vcf.gz"))
    )

    # Milvus Configuration
    MILVUS_HOST: str = "localhost"
    MILVUS_PORT: int = 19530
    MILVUS_COLLECTION: str = "genomic_evidence"

    # Embedding Configuration
    EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
    EMBEDDING_DIMENSION: int = 384
    EMBEDDING_BATCH_SIZE: int = 32

    # VEP Annotation
    VEP_CACHE_DIR: Path = DATA_DIR / "vep_cache"
    VEP_SPECIES: str = "homo_sapiens"
    VEP_ASSEMBLY: str = "GRCh38"
    USE_VEP_ANNOTATION: bool = True

    # Gene Annotation (fallback if VEP not available)
    GENE_BED_FILE: Path | None = None  # Path to gene regions BED file

    # LLM Configuration
    LLM_PROVIDER: str = "ollama"  # "openai", "anthropic", "vllm", or "ollama"
    LLM_MODEL: str = "llama3.1:70b"  # 70B for DGX Spark
    LLM_API_KEY: str | None = None  # Set via ANTHROPIC_API_KEY or OPENAI_API_KEY env var
    LLM_MAX_TOKENS: int = 4096
    LLM_TEMPERATURE: float = 0.7

    # Ollama Configuration (default for ARM64 systems like DGX Spark)
    OLLAMA_HOST: str = "http://localhost:11434"

    # vLLM Local Server (alternative for x86 systems)
    VLLM_HOST: str = "localhost"
    VLLM_PORT: int = 8080
    VLLM_MODEL: str = "meta-llama/Llama-3.1-70B-Instruct"

    # RAG Configuration
    RAG_TOP_K: int = 10  # Number of evidence items to retrieve
    RAG_SCORE_THRESHOLD: float = 0.5  # Minimum similarity score

    # Streamlit Configuration
    STREAMLIT_PORT: int = 8501
    STREAMLIT_TITLE: str = "Genomic Evidence RAG Chat"

    # Ingestion Configuration
    INGESTION_BATCH_SIZE: int = 1000
    MAX_VARIANTS_TO_INGEST: int | None = None  # None = all variants

    # Quality Filters
    MIN_VARIANT_QUAL: float = 20.0
    INCLUDE_CHROMOSOMES: list = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
        "chrX", "chrY"
    ]


# Global settings instance
settings = Settings()
