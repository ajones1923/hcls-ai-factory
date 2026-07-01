"""Tests for Milvus collection configurations in src/collections.py.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.collections import (
    ALL_COLLECTIONS,
    COLLECTION_NAMES,
    COLLECTION_SCHEMAS,
    EMBEDDING_DIM,
    WORKFLOW_COLLECTION_WEIGHTS,
    get_all_collection_names,
    get_collection_config,
    get_search_weights,
)
from src.models import NeuroWorkflowType


class TestCollectionConstants:
    """Tests for collection-level constants."""

    def test_embedding_dim(self):
        assert EMBEDDING_DIM == 384

    def test_total_collection_count(self):
        assert len(ALL_COLLECTIONS) == 14

    def test_collection_names_count(self):
        assert len(COLLECTION_NAMES) == 14

    def test_schema_count(self):
        assert len(COLLECTION_SCHEMAS) == 14


class TestCollectionConfigs:
    """Tests for individual collection configurations."""

    def test_all_have_names(self):
        for cfg in ALL_COLLECTIONS:
            assert cfg.name, f"Collection missing name: {cfg}"

    def test_all_have_descriptions(self):
        for cfg in ALL_COLLECTIONS:
            assert cfg.description, f"Collection {cfg.name} missing description"

    def test_all_have_search_weights(self):
        for cfg in ALL_COLLECTIONS:
            assert cfg.search_weight >= 0.0
            assert cfg.search_weight <= 1.0

    def test_unique_names(self):
        names = [cfg.name for cfg in ALL_COLLECTIONS]
        assert len(names) == len(set(names))

    def test_expected_collections_present(self):
        names = get_all_collection_names()
        expected = [
            "neuro_literature", "neuro_trials", "neuro_imaging",
            "neuro_electrophysiology", "neuro_degenerative",
            "neuro_cerebrovascular", "neuro_epilepsy", "neuro_oncology",
            "neuro_ms", "neuro_movement", "neuro_headache",
            "neuro_neuromuscular", "neuro_guidelines", "genomic_evidence",
        ]
        for exp in expected:
            assert exp in names, f"Missing collection: {exp}"


class TestGetCollectionConfig:
    """Tests for get_collection_config lookup."""

    def test_lookup_by_full_name(self):
        cfg = get_collection_config("neuro_literature")
        assert cfg.name == "neuro_literature"

    def test_lookup_by_alias(self):
        cfg = get_collection_config("literature")
        assert cfg.name == "neuro_literature"

    def test_invalid_name_raises(self):
        with pytest.raises(ValueError):
            get_collection_config("nonexistent_collection")


class TestSearchWeights:
    """Tests for search weight retrieval."""

    def test_default_weights(self):
        weights = get_search_weights()
        assert len(weights) == 14

    def test_workflow_weights(self):
        for wf in NeuroWorkflowType:
            if wf in WORKFLOW_COLLECTION_WEIGHTS:
                weights = get_search_weights(wf)
                assert len(weights) == 14
                total = sum(weights.values())
                assert abs(total - 1.0) < 0.05, f"Weights for {wf} sum to {total}"

    def test_general_workflow_has_all_collections(self):
        weights = get_search_weights(NeuroWorkflowType.GENERAL)
        names = get_all_collection_names()
        for name in names:
            assert name in weights
