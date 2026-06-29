"""Unit tests for Cardiology Intelligence Agent Milvus collection schemas.

Tests ALL_COLLECTIONS, WORKFLOW_COLLECTION_WEIGHTS, CollectionConfig,
and helper functions in src/collections.py.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.collections import (
    ALL_COLLECTIONS,
    CARDIO_ONCOLOGY_CONFIG,
    COLLECTION_NAMES,
    COLLECTION_SCHEMAS,
    CollectionConfig,
    DEVICES_CONFIG,
    ELECTROPHYSIOLOGY_CONFIG,
    GUIDELINES_CONFIG,
    HEART_FAILURE_CONFIG,
    HEMODYNAMICS_CONFIG,
    IMAGING_CONFIG,
    INTERVENTIONAL_CONFIG,
    LITERATURE_CONFIG,
    PREVENTION_CONFIG,
    TRIALS_CONFIG,
    VALVULAR_CONFIG,
    WORKFLOW_COLLECTION_WEIGHTS,
    get_all_collection_names,
    get_collection_config,
    get_search_weights,
)
from src.models import CardioWorkflowType


# ===================================================================
# ALL_COLLECTIONS (12 entries)
# ===================================================================


class TestAllCollections:
    """Tests for the ALL_COLLECTIONS list."""

    def test_total_count(self):
        assert len(ALL_COLLECTIONS) == 12

    def test_is_list(self):
        assert isinstance(ALL_COLLECTIONS, list)

    def test_all_entries_are_collection_config(self):
        for cfg in ALL_COLLECTIONS:
            assert isinstance(cfg, CollectionConfig)

    @pytest.mark.parametrize(
        "expected_name",
        [
            "cardio_literature",
            "cardio_trials",
            "cardio_imaging",
            "cardio_electrophysiology",
            "cardio_heart_failure",
            "cardio_valvular",
            "cardio_prevention",
            "cardio_interventional",
            "cardio_oncology",
            "cardio_devices",
            "cardio_guidelines",
            "cardio_hemodynamics",
        ],
    )
    def test_collection_name_exists(self, expected_name):
        names = [cfg.name for cfg in ALL_COLLECTIONS]
        assert expected_name in names, f"'{expected_name}' not found in ALL_COLLECTIONS"


# ===================================================================
# CollectionConfig attributes
# ===================================================================


class TestCollectionConfig:
    """Tests for CollectionConfig dataclass attributes."""

    def test_has_name(self):
        for cfg in ALL_COLLECTIONS:
            assert hasattr(cfg, "name")
            assert isinstance(cfg.name, str)
            assert len(cfg.name) > 0

    def test_has_description(self):
        for cfg in ALL_COLLECTIONS:
            assert hasattr(cfg, "description")
            assert isinstance(cfg.description, str)
            assert len(cfg.description) > 0

    def test_has_schema_fields(self):
        for cfg in ALL_COLLECTIONS:
            assert hasattr(cfg, "schema_fields")
            assert isinstance(cfg.schema_fields, list)
            assert len(cfg.schema_fields) > 0

    def test_has_estimated_records(self):
        for cfg in ALL_COLLECTIONS:
            assert hasattr(cfg, "estimated_records")
            assert isinstance(cfg.estimated_records, int)
            assert cfg.estimated_records >= 0

    def test_has_search_weight(self):
        for cfg in ALL_COLLECTIONS:
            assert hasattr(cfg, "search_weight")
            assert isinstance(cfg.search_weight, float)
            assert cfg.search_weight > 0.0

    def test_has_index_params(self):
        for cfg in ALL_COLLECTIONS:
            assert hasattr(cfg, "index_params")
            assert "metric_type" in cfg.index_params
            assert cfg.index_params["metric_type"] == "COSINE"


# ===================================================================
# Individual collection config tests
# ===================================================================


class TestIndividualConfigs:
    """Tests for each individual collection config."""

    def test_literature_config_name(self):
        assert LITERATURE_CONFIG.name == "cardio_literature"

    def test_literature_estimated_records(self):
        assert LITERATURE_CONFIG.estimated_records == 3000

    def test_trials_config_name(self):
        assert TRIALS_CONFIG.name == "cardio_trials"

    def test_trials_estimated_records(self):
        assert TRIALS_CONFIG.estimated_records == 500

    def test_imaging_config_name(self):
        assert IMAGING_CONFIG.name == "cardio_imaging"

    def test_imaging_estimated_records(self):
        assert IMAGING_CONFIG.estimated_records == 200

    def test_electrophysiology_config_name(self):
        assert ELECTROPHYSIOLOGY_CONFIG.name == "cardio_electrophysiology"

    def test_heart_failure_config_name(self):
        assert HEART_FAILURE_CONFIG.name == "cardio_heart_failure"

    def test_valvular_config_name(self):
        assert VALVULAR_CONFIG.name == "cardio_valvular"

    def test_prevention_config_name(self):
        assert PREVENTION_CONFIG.name == "cardio_prevention"

    def test_interventional_config_name(self):
        assert INTERVENTIONAL_CONFIG.name == "cardio_interventional"

    def test_cardio_oncology_config_name(self):
        assert CARDIO_ONCOLOGY_CONFIG.name == "cardio_oncology"

    def test_devices_config_name(self):
        assert DEVICES_CONFIG.name == "cardio_devices"

    def test_guidelines_config_name(self):
        assert GUIDELINES_CONFIG.name == "cardio_guidelines"

    def test_hemodynamics_config_name(self):
        assert HEMODYNAMICS_CONFIG.name == "cardio_hemodynamics"

    def test_literature_search_weight(self):
        assert LITERATURE_CONFIG.search_weight == 0.10

    def test_guidelines_search_weight(self):
        assert GUIDELINES_CONFIG.search_weight == 0.10

    def test_devices_search_weight(self):
        assert DEVICES_CONFIG.search_weight == 0.04


# ===================================================================
# COLLECTION_NAMES mapping
# ===================================================================


class TestCollectionNames:
    """Tests for the COLLECTION_NAMES alias mapping."""

    def test_is_dict(self):
        assert isinstance(COLLECTION_NAMES, dict)

    def test_count(self):
        assert len(COLLECTION_NAMES) == 12

    @pytest.mark.parametrize(
        "alias, full_name",
        [
            ("literature", "cardio_literature"),
            ("trials", "cardio_trials"),
            ("imaging", "cardio_imaging"),
            ("electrophysiology", "cardio_electrophysiology"),
            ("heart_failure", "cardio_heart_failure"),
            ("valvular", "cardio_valvular"),
            ("prevention", "cardio_prevention"),
            ("interventional", "cardio_interventional"),
            ("cardio_oncology", "cardio_oncology"),
            ("devices", "cardio_devices"),
            ("guidelines", "cardio_guidelines"),
            ("hemodynamics", "cardio_hemodynamics"),
        ],
    )
    def test_alias_mapping(self, alias, full_name):
        assert COLLECTION_NAMES[alias] == full_name


# ===================================================================
# COLLECTION_SCHEMAS
# ===================================================================


class TestCollectionSchemas:
    """Tests for the COLLECTION_SCHEMAS dict."""

    def test_is_dict(self):
        assert isinstance(COLLECTION_SCHEMAS, dict)

    def test_count(self):
        assert len(COLLECTION_SCHEMAS) == 12

    def test_keys_match_collection_names(self):
        expected = {cfg.name for cfg in ALL_COLLECTIONS}
        assert set(COLLECTION_SCHEMAS.keys()) == expected

    def test_schemas_have_description(self):
        for name, schema in COLLECTION_SCHEMAS.items():
            assert schema.description, f"Schema '{name}' has empty description"


# ===================================================================
# WORKFLOW_COLLECTION_WEIGHTS
# ===================================================================


class TestWorkflowCollectionWeights:
    """Tests for WORKFLOW_COLLECTION_WEIGHTS."""

    def test_is_dict(self):
        assert isinstance(WORKFLOW_COLLECTION_WEIGHTS, dict)

    def test_has_workflows(self):
        assert len(WORKFLOW_COLLECTION_WEIGHTS) >= 8

    def test_all_keys_are_workflow_types(self):
        for key in WORKFLOW_COLLECTION_WEIGHTS:
            assert isinstance(key, CardioWorkflowType)

    def test_all_values_are_dicts(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            assert isinstance(weights, dict), f"Weights for {wf} is not a dict"

    def test_all_weights_are_floats(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            for coll, weight in weights.items():
                assert isinstance(weight, float), (
                    f"Weight for {wf}/{coll} is not a float"
                )

    def test_all_weights_positive(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            for coll, weight in weights.items():
                assert weight > 0.0, f"Weight for {wf}/{coll} is not positive"

    def test_weights_sum_approximately_to_one(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 0.05, (
                f"Weights for {wf} sum to {total}, expected ~1.0"
            )

    def test_each_workflow_has_12_collections(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            assert len(weights) == 12, (
                f"Workflow {wf} has {len(weights)} weights, expected 12"
            )

    def test_heart_failure_workflow_dominant_collection(self):
        if CardioWorkflowType.HEART_FAILURE in WORKFLOW_COLLECTION_WEIGHTS:
            weights = WORKFLOW_COLLECTION_WEIGHTS[CardioWorkflowType.HEART_FAILURE]
            assert weights.get("cardio_heart_failure", 0) >= 0.20

    def test_arrhythmia_workflow_dominant_collection(self):
        if CardioWorkflowType.ARRHYTHMIA in WORKFLOW_COLLECTION_WEIGHTS:
            weights = WORKFLOW_COLLECTION_WEIGHTS[CardioWorkflowType.ARRHYTHMIA]
            assert weights.get("cardio_electrophysiology", 0) >= 0.20

    def test_cardio_oncology_workflow_dominant_collection(self):
        if CardioWorkflowType.CARDIO_ONCOLOGY in WORKFLOW_COLLECTION_WEIGHTS:
            weights = WORKFLOW_COLLECTION_WEIGHTS[CardioWorkflowType.CARDIO_ONCOLOGY]
            assert weights.get("cardio_oncology", 0) >= 0.20


# ===================================================================
# get_collection_config()
# ===================================================================


class TestGetCollectionConfig:
    """Tests for the get_collection_config() helper."""

    @pytest.mark.parametrize(
        "name",
        [
            "cardio_literature", "cardio_trials", "cardio_imaging",
            "cardio_electrophysiology", "cardio_heart_failure",
            "cardio_valvular", "cardio_prevention", "cardio_interventional",
            "cardio_oncology", "cardio_devices", "cardio_guidelines",
            "cardio_hemodynamics",
        ],
    )
    def test_lookup_by_full_name(self, name):
        cfg = get_collection_config(name)
        assert cfg.name == name

    @pytest.mark.parametrize(
        "alias, expected_name",
        [
            ("literature", "cardio_literature"),
            ("trials", "cardio_trials"),
            ("imaging", "cardio_imaging"),
            ("guidelines", "cardio_guidelines"),
            ("hemodynamics", "cardio_hemodynamics"),
        ],
    )
    def test_lookup_by_alias(self, alias, expected_name):
        cfg = get_collection_config(alias)
        assert cfg.name == expected_name

    def test_unknown_name_raises(self):
        with pytest.raises(ValueError, match="Unknown collection"):
            get_collection_config("nonexistent_collection")

    def test_returns_collection_config(self):
        cfg = get_collection_config("cardio_literature")
        assert isinstance(cfg, CollectionConfig)


# ===================================================================
# get_all_collection_names()
# ===================================================================


class TestGetAllCollectionNames:
    """Tests for the get_all_collection_names() helper."""

    def test_returns_list(self):
        names = get_all_collection_names()
        assert isinstance(names, list)

    def test_returns_12_names(self):
        names = get_all_collection_names()
        assert len(names) == 12

    def test_all_names_are_strings(self):
        for name in get_all_collection_names():
            assert isinstance(name, str)

    def test_contains_expected_names(self):
        names = get_all_collection_names()
        assert "cardio_literature" in names
        assert "cardio_guidelines" in names
        assert "cardio_hemodynamics" in names


# ===================================================================
# get_search_weights()
# ===================================================================


class TestGetSearchWeights:
    """Tests for the get_search_weights() helper."""

    def test_returns_dict(self):
        weights = get_search_weights()
        assert isinstance(weights, dict)

    def test_default_weights_have_12_entries(self):
        weights = get_search_weights()
        assert len(weights) == 12

    def test_default_weights_all_positive(self):
        weights = get_search_weights()
        for name, w in weights.items():
            assert w > 0.0, f"Default weight for '{name}' is not positive"

    def test_none_workflow_returns_defaults(self):
        weights = get_search_weights(None)
        assert len(weights) == 12

    def test_workflow_weights_override(self):
        """When a valid workflow is passed, weights should come from WORKFLOW_COLLECTION_WEIGHTS."""
        for wf_type in WORKFLOW_COLLECTION_WEIGHTS:
            weights = get_search_weights(wf_type)
            assert len(weights) == 12
            expected = WORKFLOW_COLLECTION_WEIGHTS[wf_type]
            for coll, w in expected.items():
                assert weights[coll] == w

    def test_returns_copy_not_reference(self):
        w1 = get_search_weights()
        w2 = get_search_weights()
        w1["cardio_literature"] = 999.0
        assert w2.get("cardio_literature") != 999.0
