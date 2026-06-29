"""Tests for Milvus collection configurations in src/collections.py.

Author: Adam Jones
Date: March 2026
"""

import pytest

from src.collections import (
    ALL_COLLECTIONS,
    COLLECTION_NAMES,
    EMBEDDING_DIM,
    WORKFLOW_COLLECTION_WEIGHTS,
    get_all_collection_names,
    get_collection_config,
)
from src.models import SCWorkflowType


class TestCollectionConstants:
    """Tests for collection-level constants."""

    def test_embedding_dim(self):
        assert EMBEDDING_DIM == 384

    def test_total_collection_count(self):
        assert len(ALL_COLLECTIONS) == 12

    def test_collection_names_count(self):
        assert len(COLLECTION_NAMES) == 12


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
            "sc_cell_types",
            "sc_markers",
            "sc_spatial",
            "sc_tme",
            "sc_drug_response",
            "sc_literature",
            "sc_methods",
            "sc_datasets",
            "sc_trajectories",
            "sc_pathways",
            "sc_clinical",
            "genomic_evidence",
        ]
        for name in expected:
            assert name in names, f"Missing collection: {name}"


class TestCollectionLookup:
    """Tests for collection lookup functions."""

    def test_get_collection_config(self):
        cfg = get_collection_config("sc_cell_types")
        assert cfg is not None
        assert cfg.name == "sc_cell_types"

    def test_get_collection_config_missing(self):
        with pytest.raises(KeyError):
            get_collection_config("nonexistent_collection")

    def test_get_all_collection_names(self):
        names = get_all_collection_names()
        assert len(names) == 12
        assert "sc_cell_types" in names
        assert "sc_tme" in names


class TestWorkflowCollectionWeights:
    """Tests for workflow-specific collection weights."""

    def test_all_workflows_have_weights(self):
        for wf in SCWorkflowType:
            assert wf in WORKFLOW_COLLECTION_WEIGHTS, (
                f"Workflow {wf.value} missing from WORKFLOW_COLLECTION_WEIGHTS"
            )

    def test_weights_reference_valid_collections(self):
        all_names = set(get_all_collection_names())
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            for coll_name in weights:
                assert coll_name in all_names, (
                    f"Workflow {wf.value} references unknown collection: {coll_name}"
                )

    def test_weight_sums(self):
        for wf, weights in WORKFLOW_COLLECTION_WEIGHTS.items():
            total = sum(weights.values())
            assert abs(total - 1.0) < 0.05, (
                f"Workflow {wf.value} weights sum to {total}, expected ~1.0"
            )
