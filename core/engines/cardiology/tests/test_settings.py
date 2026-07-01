"""Tests for the Cardiology Intelligence Agent configuration.

Validates CardioSettings creation, default values for ports, Milvus
configuration, collection names, embedding dimension, LLM model,
validation logic, weight constraints, and environment prefix.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import unittest

from config.settings import CardioSettings


# =====================================================================
# CardioSettings creation
# =====================================================================


class TestCardioSettingsCreation(unittest.TestCase):
    """Tests for CardioSettings instantiation."""

    def test_create_settings(self):
        settings = CardioSettings()
        self.assertIsInstance(settings, CardioSettings)

    def test_settings_is_base_settings(self):
        from pydantic_settings import BaseSettings
        self.assertTrue(issubclass(CardioSettings, BaseSettings))


# =====================================================================
# Default port values
# =====================================================================


class TestDefaultPorts(unittest.TestCase):
    """Tests for default port configuration."""

    def test_api_port_default(self):
        settings = CardioSettings()
        self.assertEqual(settings.API_PORT, 8126)

    def test_streamlit_port_default(self):
        settings = CardioSettings()
        self.assertEqual(settings.STREAMLIT_PORT, 8536)

    def test_milvus_port_default(self):
        settings = CardioSettings()
        self.assertEqual(settings.MILVUS_PORT, 19530)

    def test_api_port_in_valid_range(self):
        settings = CardioSettings()
        self.assertGreaterEqual(settings.API_PORT, 1024)
        self.assertLessEqual(settings.API_PORT, 65535)

    def test_streamlit_port_in_valid_range(self):
        settings = CardioSettings()
        self.assertGreaterEqual(settings.STREAMLIT_PORT, 1024)
        self.assertLessEqual(settings.STREAMLIT_PORT, 65535)

    def test_ports_are_different(self):
        settings = CardioSettings()
        self.assertNotEqual(settings.API_PORT, settings.STREAMLIT_PORT)


# =====================================================================
# Milvus configuration
# =====================================================================


class TestMilvusConfig(unittest.TestCase):
    """Tests for Milvus connection defaults."""

    def test_milvus_host_default(self):
        settings = CardioSettings()
        self.assertEqual(settings.MILVUS_HOST, "localhost")

    def test_milvus_host_not_empty(self):
        settings = CardioSettings()
        self.assertTrue(len(settings.MILVUS_HOST) > 0)


# =====================================================================
# Collection names
# =====================================================================


class TestCollectionNames(unittest.TestCase):
    """Tests for Milvus collection name defaults."""

    def test_collection_literature(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_LITERATURE, "cardio_literature")

    def test_collection_trials(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_TRIALS, "cardio_trials")

    def test_collection_imaging(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_IMAGING, "cardio_imaging")

    def test_collection_electrophysiology(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_ELECTROPHYSIOLOGY, "cardio_electrophysiology")

    def test_collection_heart_failure(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_HEART_FAILURE, "cardio_heart_failure")

    def test_collection_valvular(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_VALVULAR, "cardio_valvular")

    def test_collection_prevention(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_PREVENTION, "cardio_prevention")

    def test_collection_interventional(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_INTERVENTIONAL, "cardio_interventional")

    def test_collection_oncology(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_ONCOLOGY, "cardio_oncology")

    def test_collection_devices(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_DEVICES, "cardio_devices")

    def test_collection_guidelines(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_GUIDELINES, "cardio_guidelines")

    def test_collection_hemodynamics(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_HEMODYNAMICS, "cardio_hemodynamics")

    def test_collection_genomic(self):
        settings = CardioSettings()
        self.assertEqual(settings.COLLECTION_GENOMIC, "genomic_evidence")


# =====================================================================
# Embedding configuration
# =====================================================================


class TestEmbeddingConfig(unittest.TestCase):
    """Tests for embedding model defaults."""

    def test_embedding_dimension(self):
        settings = CardioSettings()
        self.assertEqual(settings.EMBEDDING_DIMENSION, 384)

    def test_embedding_model(self):
        settings = CardioSettings()
        self.assertEqual(settings.EMBEDDING_MODEL, "BAAI/bge-small-en-v1.5")

    def test_embedding_batch_size(self):
        settings = CardioSettings()
        self.assertEqual(settings.EMBEDDING_BATCH_SIZE, 32)


# =====================================================================
# LLM configuration
# =====================================================================


class TestLLMConfig(unittest.TestCase):
    """Tests for LLM configuration defaults."""

    def test_llm_model(self):
        settings = CardioSettings()
        self.assertEqual(settings.LLM_MODEL, "claude-sonnet-4-6")

    def test_llm_provider(self):
        settings = CardioSettings()
        self.assertEqual(settings.LLM_PROVIDER, "anthropic")


# =====================================================================
# validate() method
# =====================================================================


class TestValidateMethod(unittest.TestCase):
    """Tests for CardioSettings.validate() method."""

    def test_validate_returns_list(self):
        settings = CardioSettings()
        result = settings.validate()
        self.assertIsInstance(result, list)

    def test_validate_all_entries_are_strings(self):
        settings = CardioSettings()
        result = settings.validate()
        for item in result:
            self.assertIsInstance(item, str)

    def test_validate_catches_port_conflict(self):
        settings = CardioSettings()
        settings.API_PORT = 8000
        settings.STREAMLIT_PORT = 8000
        issues = settings.validate()
        port_issues = [i for i in issues if "port conflict" in i.lower() or "both" in i.lower()]
        self.assertGreater(len(port_issues), 0)

    def test_validate_catches_empty_milvus_host(self):
        settings = CardioSettings()
        settings.MILVUS_HOST = ""
        issues = settings.validate()
        host_issues = [i for i in issues if "MILVUS_HOST" in i]
        self.assertGreater(len(host_issues), 0)

    def test_validate_catches_invalid_milvus_port(self):
        settings = CardioSettings()
        settings.MILVUS_PORT = 0
        issues = settings.validate()
        port_issues = [i for i in issues if "MILVUS_PORT" in i]
        self.assertGreater(len(port_issues), 0)

    def test_validate_catches_missing_api_key(self):
        settings = CardioSettings()
        settings.ANTHROPIC_API_KEY = None
        issues = settings.validate()
        key_issues = [i for i in issues if "ANTHROPIC_API_KEY" in i]
        self.assertGreater(len(key_issues), 0)

    def test_validate_catches_empty_embedding_model(self):
        settings = CardioSettings()
        settings.EMBEDDING_MODEL = ""
        issues = settings.validate()
        model_issues = [i for i in issues if "EMBEDDING_MODEL" in i]
        self.assertGreater(len(model_issues), 0)

    def test_validate_catches_port_below_range(self):
        settings = CardioSettings()
        settings.API_PORT = 80
        issues = settings.validate()
        port_issues = [i for i in issues if "API_PORT" in i]
        self.assertGreater(len(port_issues), 0)


# =====================================================================
# Weights sum to ~1.0
# =====================================================================


class TestWeightsSum(unittest.TestCase):
    """Tests for collection search weight constraints."""

    def test_weights_sum_approximately_one(self):
        settings = CardioSettings()
        weights = [
            settings.WEIGHT_LITERATURE,
            settings.WEIGHT_TRIALS,
            settings.WEIGHT_IMAGING,
            settings.WEIGHT_ELECTROPHYSIOLOGY,
            settings.WEIGHT_HEART_FAILURE,
            settings.WEIGHT_VALVULAR,
            settings.WEIGHT_PREVENTION,
            settings.WEIGHT_INTERVENTIONAL,
            settings.WEIGHT_ONCOLOGY,
            settings.WEIGHT_DEVICES,
            settings.WEIGHT_GUIDELINES,
            settings.WEIGHT_HEMODYNAMICS,
            settings.WEIGHT_GENOMIC,
        ]
        total = sum(weights)
        self.assertAlmostEqual(total, 1.0, places=1)

    def test_all_weights_non_negative(self):
        settings = CardioSettings()
        weight_attrs = [
            attr for attr in dir(settings)
            if attr.startswith("WEIGHT_") and isinstance(getattr(settings, attr), float)
        ]
        for attr in weight_attrs:
            val = getattr(settings, attr)
            self.assertGreaterEqual(val, 0.0, f"{attr} is negative: {val}")

    def test_weights_within_valid_range(self):
        settings = CardioSettings()
        weight_attrs = [
            attr for attr in dir(settings)
            if attr.startswith("WEIGHT_") and isinstance(getattr(settings, attr), float)
        ]
        for attr in weight_attrs:
            val = getattr(settings, attr)
            self.assertLessEqual(val, 1.0, f"{attr} exceeds 1.0: {val}")


# =====================================================================
# env_prefix
# =====================================================================


class TestEnvPrefix(unittest.TestCase):
    """Tests for Pydantic Settings environment prefix."""

    def test_env_prefix_is_cardio(self):
        config = CardioSettings.model_config
        self.assertEqual(config.get("env_prefix"), "CARDIO_")


# =====================================================================
# Other defaults
# =====================================================================


class TestOtherDefaults(unittest.TestCase):
    """Tests for miscellaneous default values."""

    def test_top_k_per_collection(self):
        settings = CardioSettings()
        self.assertEqual(settings.TOP_K_PER_COLLECTION, 5)

    def test_score_threshold(self):
        settings = CardioSettings()
        self.assertEqual(settings.SCORE_THRESHOLD, 0.4)

    def test_metrics_enabled(self):
        settings = CardioSettings()
        self.assertTrue(settings.METRICS_ENABLED)

    def test_ingest_disabled_by_default(self):
        settings = CardioSettings()
        self.assertFalse(settings.INGEST_ENABLED)

    def test_ingest_schedule_hours(self):
        settings = CardioSettings()
        self.assertEqual(settings.INGEST_SCHEDULE_HOURS, 168)

    def test_max_conversation_context(self):
        settings = CardioSettings()
        self.assertEqual(settings.MAX_CONVERSATION_CONTEXT, 3)

    def test_api_host(self):
        settings = CardioSettings()
        self.assertEqual(settings.API_HOST, "0.0.0.0")

    def test_cors_origins_not_empty(self):
        settings = CardioSettings()
        self.assertTrue(len(settings.CORS_ORIGINS) > 0)

    def test_validate_or_warn_exists(self):
        settings = CardioSettings()
        self.assertTrue(hasattr(settings, "validate_or_warn"))
        self.assertTrue(callable(settings.validate_or_warn))

    def test_max_request_size_mb(self):
        settings = CardioSettings()
        self.assertEqual(settings.MAX_REQUEST_SIZE_MB, 10)


if __name__ == "__main__":
    unittest.main()
