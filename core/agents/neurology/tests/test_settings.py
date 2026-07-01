"""Tests for configuration settings in config/settings.py.

Author: Adam Jones
Date: March 2026
"""


from config.settings import NeuroSettings, settings


class TestNeuroSettings:
    """Tests for NeuroSettings configuration."""

    def test_settings_instance_exists(self):
        assert settings is not None

    def test_is_neuro_settings(self):
        assert isinstance(settings, NeuroSettings)

    def test_default_milvus_host(self):
        assert settings.MILVUS_HOST == "localhost"

    def test_default_milvus_port(self):
        assert settings.MILVUS_PORT == 19530

    def test_embedding_model(self):
        assert settings.EMBEDDING_MODEL == "BAAI/bge-small-en-v1.5"

    def test_embedding_dimension(self):
        assert settings.EMBEDDING_DIMENSION == 384

    def test_api_port(self):
        assert settings.API_PORT == 8528

    def test_streamlit_port(self):
        assert settings.STREAMLIT_PORT == 8529

    def test_ports_differ(self):
        assert settings.API_PORT != settings.STREAMLIT_PORT

    def test_score_threshold(self):
        assert 0.0 <= settings.SCORE_THRESHOLD <= 1.0

    def test_collection_names(self):
        assert settings.COLLECTION_LITERATURE == "neuro_literature"
        assert settings.COLLECTION_IMAGING == "neuro_imaging"
        assert settings.COLLECTION_EPILEPSY == "neuro_epilepsy"

    def test_cross_agent_urls(self):
        assert settings.CARDIOLOGY_AGENT_URL.startswith("http")
        assert settings.RARE_DISEASE_AGENT_URL.startswith("http")

    def test_cross_agent_timeout(self):
        assert settings.CROSS_AGENT_TIMEOUT > 0


class TestSettingsValidation:
    """Tests for settings.validate() method."""

    def test_validate_returns_list(self):
        issues = settings.validate()
        assert isinstance(issues, list)

    def test_validate_or_warn_runs(self):
        # Should not raise
        settings.validate_or_warn()

    def test_weight_sum_close_to_one(self):
        weight_attrs = [
            attr for attr in dir(settings)
            if attr.startswith("WEIGHT_") and isinstance(getattr(settings, attr), float)
        ]
        weights = [getattr(settings, attr) for attr in weight_attrs]
        total = sum(weights)
        assert abs(total - 1.0) < 0.05, f"Weights sum to {total}"

    def test_all_top_k_positive(self):
        top_k_attrs = [
            attr for attr in dir(settings)
            if attr.startswith("TOP_K_") and isinstance(getattr(settings, attr), int)
        ]
        for attr in top_k_attrs:
            val = getattr(settings, attr)
            assert val > 0, f"{attr}={val} should be positive"

    def test_env_prefix(self):
        assert settings.model_config.get("env_prefix") == "NEURO_"
