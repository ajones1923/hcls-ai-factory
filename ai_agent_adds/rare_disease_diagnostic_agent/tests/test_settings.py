"""Tests for config/settings.py.

Covers:
  - Default values
  - Weight validation
  - Port validation

Author: Adam Jones
Date: March 2026
"""


from config.settings import RareDiseaseSettings


class TestRareDiseaseSettings:
    """Tests for RareDiseaseSettings configuration."""

    def test_default_milvus_host(self):
        s = RareDiseaseSettings()
        assert s.MILVUS_HOST == "localhost"

    def test_default_milvus_port(self):
        s = RareDiseaseSettings()
        assert s.MILVUS_PORT == 19530

    def test_default_api_port(self):
        s = RareDiseaseSettings()
        assert s.API_PORT == 8134

    def test_default_streamlit_port(self):
        s = RareDiseaseSettings()
        assert s.STREAMLIT_PORT == 8544

    def test_default_embedding_model(self):
        s = RareDiseaseSettings()
        assert s.EMBEDDING_MODEL == "BAAI/bge-small-en-v1.5"

    def test_default_embedding_dimension(self):
        s = RareDiseaseSettings()
        assert s.EMBEDDING_DIMENSION == 384

    def test_14_collections_defined(self):
        s = RareDiseaseSettings()
        collection_attrs = [
            attr for attr in dir(s)
            if attr.startswith("COLLECTION_") and not attr.startswith("__")
        ]
        assert len(collection_attrs) == 14

    def test_weight_sum_approximately_one(self):
        s = RareDiseaseSettings()
        weight_attrs = [
            attr for attr in dir(s)
            if attr.startswith("WEIGHT_") and isinstance(getattr(s, attr), float)
        ]
        total = sum(getattr(s, attr) for attr in weight_attrs)
        assert abs(total - 1.0) < 0.05, f"Weight sum is {total}"

    def test_all_weights_non_negative(self):
        s = RareDiseaseSettings()
        weight_attrs = [
            attr for attr in dir(s)
            if attr.startswith("WEIGHT_") and isinstance(getattr(s, attr), float)
        ]
        for attr in weight_attrs:
            val = getattr(s, attr)
            assert val >= 0, f"{attr} = {val} is negative"

    def test_validate_no_critical_issues_with_defaults(self):
        s = RareDiseaseSettings()
        issues = s.validate()
        # Only expected issue is missing ANTHROPIC_API_KEY
        critical_issues = [i for i in issues if "MILVUS" in i or "PORT" in i]
        assert len(critical_issues) == 0

    def test_validate_warns_missing_api_key(self):
        s = RareDiseaseSettings()
        issues = s.validate()
        api_key_issues = [i for i in issues if "ANTHROPIC_API_KEY" in i]
        assert len(api_key_issues) == 1

    def test_validate_port_conflict(self):
        s = RareDiseaseSettings(API_PORT=8134, STREAMLIT_PORT=8134)
        issues = s.validate()
        conflict = [i for i in issues if "port conflict" in i.lower()]
        assert len(conflict) == 1

    def test_validate_invalid_milvus_port(self):
        s = RareDiseaseSettings(MILVUS_PORT=0)
        issues = s.validate()
        port_issues = [i for i in issues if "MILVUS_PORT" in i]
        assert len(port_issues) == 1

    def test_default_score_threshold(self):
        s = RareDiseaseSettings()
        assert s.SCORE_THRESHOLD == 0.4

    def test_default_cross_agent_timeout(self):
        s = RareDiseaseSettings()
        assert s.CROSS_AGENT_TIMEOUT == 30

    def test_default_ingest_not_enabled(self):
        s = RareDiseaseSettings()
        assert s.INGEST_ENABLED is False

    def test_default_metrics_enabled(self):
        s = RareDiseaseSettings()
        assert s.METRICS_ENABLED is True

    def test_env_prefix(self):
        assert RareDiseaseSettings.model_config["env_prefix"] == "RD_"
