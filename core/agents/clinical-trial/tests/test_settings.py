"""Tests for configuration settings in config/settings.py.

Covers:
  - Default values
  - Validation logic
  - Environment variable prefix

Author: Adam Jones
Date: March 2026
"""


from config.settings import TrialSettings


class TestSettingsDefaults:
    """Test that default settings have sensible values."""

    def test_milvus_defaults(self):
        s = TrialSettings()
        assert s.MILVUS_HOST == "localhost"
        assert s.MILVUS_PORT == 19530

    def test_embedding_defaults(self):
        s = TrialSettings()
        assert s.EMBEDDING_MODEL == "BAAI/bge-small-en-v1.5"
        assert s.EMBEDDING_DIMENSION == 384
        assert s.EMBEDDING_BATCH_SIZE == 32

    def test_llm_defaults(self):
        s = TrialSettings()
        assert s.LLM_PROVIDER == "anthropic"
        assert "claude" in s.LLM_MODEL

    def test_api_defaults(self):
        s = TrialSettings()
        assert s.API_HOST == "0.0.0.0"
        assert s.API_PORT == 8538

    def test_streamlit_defaults(self):
        s = TrialSettings()
        assert s.STREAMLIT_PORT == 8128

    def test_scheduler_defaults(self):
        s = TrialSettings()
        assert s.INGEST_SCHEDULE_HOURS == 24
        assert s.INGEST_ENABLED is False

    def test_search_defaults(self):
        s = TrialSettings()
        assert s.TOP_K_PER_COLLECTION == 5
        assert s.SCORE_THRESHOLD == 0.4

    def test_conversation_defaults(self):
        s = TrialSettings()
        assert s.MAX_CONVERSATION_CONTEXT == 3

    def test_citation_defaults(self):
        s = TrialSettings()
        assert s.CITATION_HIGH_THRESHOLD == 0.75
        assert s.CITATION_MEDIUM_THRESHOLD == 0.60


class TestSettingsCollectionNames:
    """Test collection name settings."""

    def test_protocol_collection(self):
        s = TrialSettings()
        assert s.COLLECTION_PROTOCOLS == "trial_protocols"

    def test_all_14_collections(self):
        s = TrialSettings()
        collection_attrs = [
            attr for attr in dir(s)
            if attr.startswith("COLLECTION_")
        ]
        assert len(collection_attrs) == 14


class TestSettingsWeights:
    """Test collection weight settings."""

    def test_weight_sum(self):
        s = TrialSettings()
        weight_attrs = [
            attr for attr in dir(s)
            if attr.startswith("WEIGHT_") and isinstance(getattr(s, attr), float)
        ]
        total = sum(getattr(s, attr) for attr in weight_attrs)
        assert abs(total - 1.0) <= 0.05, (
            f"Settings weights sum to {total:.4f}, expected ~1.0"
        )

    def test_all_weights_positive(self):
        s = TrialSettings()
        weight_attrs = [
            attr for attr in dir(s)
            if attr.startswith("WEIGHT_") and isinstance(getattr(s, attr), float)
        ]
        for attr in weight_attrs:
            val = getattr(s, attr)
            assert val >= 0, f"{attr}={val} is negative"


class TestSettingsValidation:
    """Test the validate() method."""

    def test_default_validation_warnings(self):
        s = TrialSettings()
        issues = s.validate()
        # Should warn about missing ANTHROPIC_API_KEY at minimum
        api_key_warning = any("ANTHROPIC_API_KEY" in i for i in issues)
        assert api_key_warning, "Should warn about missing ANTHROPIC_API_KEY"

    def test_valid_milvus_config(self):
        s = TrialSettings()
        issues = s.validate()
        milvus_errors = [i for i in issues if "MILVUS_HOST" in i]
        assert len(milvus_errors) == 0

    def test_port_conflict_detection(self):
        s = TrialSettings(API_PORT=8538, STREAMLIT_PORT=8538)
        issues = s.validate()
        conflict_warning = any("port conflict" in i.lower() for i in issues)
        assert conflict_warning


class TestSettingsEnvPrefix:
    """Test environment variable prefix."""

    def test_env_prefix(self):
        assert TrialSettings.model_config.get("env_prefix") == "TRIAL_"
