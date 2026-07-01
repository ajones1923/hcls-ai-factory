"""Tests for configuration settings in config/settings.py.

Author: Adam Jones
Date: March 2026
"""


from config.settings import SingleCellSettings, settings


class TestSingleCellSettings:
    """Tests for SingleCellSettings configuration."""

    def test_settings_instance_exists(self):
        assert settings is not None

    def test_is_single_cell_settings(self):
        assert isinstance(settings, SingleCellSettings)

    def test_default_milvus_host(self):
        assert settings.MILVUS_HOST == "localhost"

    def test_default_milvus_port(self):
        assert settings.MILVUS_PORT == 19530

    def test_embedding_model(self):
        assert settings.EMBEDDING_MODEL == "BAAI/bge-small-en-v1.5"

    def test_embedding_dimension(self):
        assert settings.EMBEDDING_DIMENSION == 384

    def test_api_port(self):
        assert settings.API_PORT == 8540

    def test_streamlit_port(self):
        assert settings.STREAMLIT_PORT == 8130

    def test_ports_differ(self):
        assert settings.API_PORT != settings.STREAMLIT_PORT

    def test_score_threshold(self):
        assert 0.0 <= settings.SCORE_THRESHOLD <= 1.0

    def test_collection_names(self):
        assert settings.COLLECTION_CELL_TYPES == "sc_cell_types"
        assert settings.COLLECTION_MARKERS == "sc_markers"
        assert settings.COLLECTION_TME == "sc_tme"

    def test_cross_agent_urls(self):
        assert settings.ONCOLOGY_AGENT_URL.startswith("http")
        assert settings.BIOMARKER_AGENT_URL.startswith("http")

    def test_cross_agent_timeout(self):
        assert settings.CROSS_AGENT_TIMEOUT > 0


class TestSettingsValidation:
    """Tests for settings.validate() method."""

    def test_validate_returns_list(self):
        issues = settings.validate()
        assert isinstance(issues, list)

    def test_validate_or_warn_runs(self):
        """validate_or_warn should not raise."""
        settings.validate_or_warn()

    def test_weights_sum_approximately_one(self):
        weight_attrs = [
            attr for attr in dir(settings)
            if attr.startswith("WEIGHT_") and isinstance(getattr(settings, attr), float)
        ]
        weights = [getattr(settings, attr) for attr in weight_attrs]
        if weights:
            total = sum(weights)
            assert abs(total - 1.0) < 0.06, (
                f"Collection weights sum to {total}, expected ~1.0"
            )

    def test_valid_port_ranges(self):
        assert 1024 <= settings.API_PORT <= 65535
        assert 1024 <= settings.STREAMLIT_PORT <= 65535

    def test_embedding_batch_size_positive(self):
        assert settings.EMBEDDING_BATCH_SIZE > 0
