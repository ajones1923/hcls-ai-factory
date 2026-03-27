"""Tests for config/settings.py — PGxSettings.

Tests default values, env_prefix, collection names, weight sum,
port assignments, and embedding configuration.

Author: Adam Jones
Date: March 2026
"""

import pytest
from config.settings import PGxSettings, settings


# ═══════════════════════════════════════════════════════════════════════
# SINGLETON
# ═══════════════════════════════════════════════════════════════════════

class TestSettingsSingleton:
    def test_settings_is_instance(self):
        assert isinstance(settings, PGxSettings)


# ═══════════════════════════════════════════════════════════════════════
# ENV PREFIX
# ═══════════════════════════════════════════════════════════════════════

class TestEnvPrefix:
    def test_env_prefix(self):
        assert settings.model_config["env_prefix"] == "PGX_"


# ═══════════════════════════════════════════════════════════════════════
# COLLECTION NAMES
# ═══════════════════════════════════════════════════════════════════════

class TestCollectionNames:
    def test_gene_reference(self):
        assert settings.COLLECTION_GENE_REFERENCE == "pgx_gene_reference"

    def test_drug_guidelines(self):
        assert settings.COLLECTION_DRUG_GUIDELINES == "pgx_drug_guidelines"

    def test_drug_interactions(self):
        assert settings.COLLECTION_DRUG_INTERACTIONS == "pgx_drug_interactions"

    def test_hla_hypersensitivity(self):
        assert settings.COLLECTION_HLA_HYPERSENSITIVITY == "pgx_hla_hypersensitivity"

    def test_phenoconversion(self):
        assert settings.COLLECTION_PHENOCONVERSION == "pgx_phenoconversion"

    def test_dosing_algorithms(self):
        assert settings.COLLECTION_DOSING_ALGORITHMS == "pgx_dosing_algorithms"

    def test_clinical_evidence(self):
        assert settings.COLLECTION_CLINICAL_EVIDENCE == "pgx_clinical_evidence"

    def test_population_data(self):
        assert settings.COLLECTION_POPULATION_DATA == "pgx_population_data"

    def test_clinical_trials(self):
        assert settings.COLLECTION_CLINICAL_TRIALS == "pgx_clinical_trials"

    def test_fda_labels(self):
        assert settings.COLLECTION_FDA_LABELS == "pgx_fda_labels"

    def test_drug_alternatives(self):
        assert settings.COLLECTION_DRUG_ALTERNATIVES == "pgx_drug_alternatives"

    def test_patient_profiles(self):
        assert settings.COLLECTION_PATIENT_PROFILES == "pgx_patient_profiles"

    def test_implementation(self):
        assert settings.COLLECTION_IMPLEMENTATION == "pgx_implementation"

    def test_education(self):
        assert settings.COLLECTION_EDUCATION == "pgx_education"

    def test_genomic(self):
        assert settings.COLLECTION_GENOMIC == "genomic_evidence"


# ═══════════════════════════════════════════════════════════════════════
# WEIGHTS
# ═══════════════════════════════════════════════════════════════════════

class TestWeights:
    def test_weights_sum_approximately_one(self):
        total = (
            settings.WEIGHT_GENE_REFERENCE
            + settings.WEIGHT_DRUG_GUIDELINES
            + settings.WEIGHT_DRUG_INTERACTIONS
            + settings.WEIGHT_HLA_HYPERSENSITIVITY
            + settings.WEIGHT_PHENOCONVERSION
            + settings.WEIGHT_DOSING_ALGORITHMS
            + settings.WEIGHT_CLINICAL_EVIDENCE
            + settings.WEIGHT_POPULATION_DATA
            + settings.WEIGHT_CLINICAL_TRIALS
            + settings.WEIGHT_FDA_LABELS
            + settings.WEIGHT_DRUG_ALTERNATIVES
            + settings.WEIGHT_PATIENT_PROFILES
            + settings.WEIGHT_IMPLEMENTATION
            + settings.WEIGHT_EDUCATION
            + settings.WEIGHT_GENOMIC
        )
        assert abs(total - 1.0) < 0.02

    def test_drug_guidelines_highest(self):
        assert settings.WEIGHT_DRUG_GUIDELINES >= 0.10

    def test_each_weight_positive(self):
        weights = [
            settings.WEIGHT_GENE_REFERENCE,
            settings.WEIGHT_DRUG_GUIDELINES,
            settings.WEIGHT_DRUG_INTERACTIONS,
            settings.WEIGHT_HLA_HYPERSENSITIVITY,
            settings.WEIGHT_PHENOCONVERSION,
            settings.WEIGHT_DOSING_ALGORITHMS,
            settings.WEIGHT_CLINICAL_EVIDENCE,
            settings.WEIGHT_POPULATION_DATA,
            settings.WEIGHT_CLINICAL_TRIALS,
            settings.WEIGHT_FDA_LABELS,
            settings.WEIGHT_DRUG_ALTERNATIVES,
            settings.WEIGHT_PATIENT_PROFILES,
            settings.WEIGHT_IMPLEMENTATION,
            settings.WEIGHT_EDUCATION,
            settings.WEIGHT_GENOMIC,
        ]
        for w in weights:
            assert w > 0


# ═══════════════════════════════════════════════════════════════════════
# PORTS
# ═══════════════════════════════════════════════════════════════════════

class TestPorts:
    def test_api_port(self):
        assert settings.API_PORT == 8107

    def test_streamlit_port(self):
        assert settings.STREAMLIT_PORT == 8507


# ═══════════════════════════════════════════════════════════════════════
# EMBEDDINGS
# ═══════════════════════════════════════════════════════════════════════

class TestEmbeddings:
    def test_embedding_dimension(self):
        assert settings.EMBEDDING_DIMENSION == 384

    def test_embedding_model(self):
        assert "bge-small-en" in settings.EMBEDDING_MODEL

    def test_batch_size(self):
        assert settings.EMBEDDING_BATCH_SIZE == 32


# ═══════════════════════════════════════════════════════════════════════
# MILVUS
# ═══════════════════════════════════════════════════════════════════════

class TestMilvus:
    def test_milvus_host(self):
        assert settings.MILVUS_HOST == "localhost"

    def test_milvus_port(self):
        assert settings.MILVUS_PORT == 19530


# ═══════════════════════════════════════════════════════════════════════
# RAG SEARCH
# ═══════════════════════════════════════════════════════════════════════

class TestRAGSearch:
    def test_top_k(self):
        assert settings.TOP_K_PER_COLLECTION == 5

    def test_score_threshold(self):
        assert settings.SCORE_THRESHOLD == 0.4


# ═══════════════════════════════════════════════════════════════════════
# SCHEDULER
# ═══════════════════════════════════════════════════════════════════════

class TestSchedulerSettings:
    def test_ingest_schedule_hours(self):
        assert settings.INGEST_SCHEDULE_HOURS == 168

    def test_ingest_disabled_by_default(self):
        assert settings.INGEST_ENABLED is False


# ═══════════════════════════════════════════════════════════════════════
# OTHER DEFAULTS
# ═══════════════════════════════════════════════════════════════════════

class TestOtherDefaults:
    def test_llm_provider(self):
        assert settings.LLM_PROVIDER == "anthropic"

    def test_max_conversation_context(self):
        assert settings.MAX_CONVERSATION_CONTEXT == 3

    def test_citation_thresholds(self):
        assert settings.CITATION_HIGH_THRESHOLD == 0.75
        assert settings.CITATION_MEDIUM_THRESHOLD == 0.60


# ═══════════════════════════════════════════════════════════════════════
# STARTUP VALIDATION
# ═══════════════════════════════════════════════════════════════════════

class TestValidation:
    """Tests for PGxSettings.validate() startup validation."""

    def test_default_settings_no_errors_except_optional(self):
        """Default settings should only warn about missing API key and RAG path."""
        s = PGxSettings()
        issues = s.validate()
        # Filter out expected warnings (API key not set, RAG path missing)
        unexpected = [
            i for i in issues
            if "ANTHROPIC_API_KEY" not in i and "RAG_PIPELINE_ROOT" not in i
        ]
        assert unexpected == [], f"Unexpected validation issues: {unexpected}"

    def test_invalid_milvus_port_zero(self):
        s = PGxSettings(MILVUS_PORT=0)
        issues = s.validate()
        assert any("MILVUS_PORT" in i for i in issues)

    def test_invalid_milvus_port_too_high(self):
        s = PGxSettings(MILVUS_PORT=70000)
        issues = s.validate()
        assert any("MILVUS_PORT" in i for i in issues)

    def test_empty_milvus_host(self):
        s = PGxSettings(MILVUS_HOST="")
        issues = s.validate()
        assert any("MILVUS_HOST" in i for i in issues)

    def test_negative_weight_produces_warning(self):
        s = PGxSettings(WEIGHT_GENE_REFERENCE=-0.5)
        issues = s.validate()
        assert any("WEIGHT_GENE_REFERENCE" in i and "negative" in i for i in issues)

    def test_port_conflict_produces_warning(self):
        s = PGxSettings(API_PORT=8000, STREAMLIT_PORT=8000)
        issues = s.validate()
        assert any("port conflict" in i.lower() for i in issues)

    def test_api_port_below_1024_produces_warning(self):
        s = PGxSettings(API_PORT=80)
        issues = s.validate()
        assert any("API_PORT" in i for i in issues)

    def test_empty_embedding_model_produces_warning(self):
        s = PGxSettings(EMBEDDING_MODEL="")
        issues = s.validate()
        assert any("EMBEDDING_MODEL" in i for i in issues)

    def test_validate_never_raises(self):
        """Validation must never raise, even with extreme values."""
        s = PGxSettings(
            MILVUS_HOST="",
            MILVUS_PORT=-1,
            EMBEDDING_MODEL="",
            API_PORT=0,
            STREAMLIT_PORT=0,
            WEIGHT_GENE_REFERENCE=-99.0,
        )
        # Should return issues, not raise
        issues = s.validate()
        assert len(issues) > 0

    def test_validate_or_warn_logs(self, caplog):
        """validate_or_warn() should log each issue as a warning."""
        import logging

        s = PGxSettings(MILVUS_HOST="")
        with caplog.at_level(logging.WARNING, logger="config.settings"):
            s.validate_or_warn()
        assert any("MILVUS_HOST" in r.message for r in caplog.records)
