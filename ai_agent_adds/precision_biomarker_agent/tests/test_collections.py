"""Tests for Biomarker Intelligence Agent collection management module.

Validates collection schemas, models, and BiomarkerCollectionManager
operations — all with mocked pymilvus to avoid requiring a running
Milvus instance.

Author: Adam Jones
Date: March 2026
"""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so ``from src.`` imports work
# regardless of how pytest is invoked.
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.collections import (
    COLLECTION_MODELS,
    COLLECTION_SCHEMAS,
    EMBEDDING_DIM,
    BiomarkerCollectionManager,
)

# =====================================================================
# EXPECTED COLLECTION NAMES
# =====================================================================

ALL_14_COLLECTIONS = [
    "biomarker_reference",
    "biomarker_genetic_variants",
    "biomarker_pgx_rules",
    "biomarker_disease_trajectories",
    "biomarker_clinical_evidence",
    "biomarker_nutrition",
    "biomarker_drug_interactions",
    "biomarker_aging_markers",
    "biomarker_genotype_adjustments",
    "biomarker_monitoring",
    "biomarker_critical_values",
    "biomarker_discordance_rules",
    "biomarker_aj_carrier_screening",
    "genomic_evidence",
]


# =====================================================================
# COLLECTION SCHEMA TESTS
# =====================================================================


class TestCollectionSchemas:
    """Validate that all 14 collection schemas are defined correctly."""

    def test_all_14_schemas_defined(self):
        """COLLECTION_SCHEMAS contains exactly 14 entries."""
        assert len(COLLECTION_SCHEMAS) == 14

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_schema_exists(self, name):
        """Each expected collection name is present in COLLECTION_SCHEMAS."""
        assert name in COLLECTION_SCHEMAS

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_schema_has_id_field(self, name):
        """Each schema has an 'id' field."""
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f.name for f in schema.fields]
        assert "id" in field_names

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_schema_has_embedding_field(self, name):
        """Each schema has an 'embedding' field."""
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f.name for f in schema.fields]
        assert "embedding" in field_names

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_schema_has_text_field(self, name):
        """Each schema has a text content field ('text_chunk' or 'text_summary')."""
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f.name for f in schema.fields]
        assert "text_chunk" in field_names or "text_summary" in field_names

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_schema_embedding_dim(self, name):
        """Each schema's embedding field has the correct dimension."""
        schema = COLLECTION_SCHEMAS[name]
        for field in schema.fields:
            if field.name == "embedding":
                assert field.params.get("dim", None) == EMBEDDING_DIM
                break

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_schema_has_primary_key(self, name):
        """Each schema has exactly one primary key field."""
        schema = COLLECTION_SCHEMAS[name]
        pk_fields = [f for f in schema.fields if f.is_primary]
        assert len(pk_fields) == 1


# =====================================================================
# COLLECTION MODEL TESTS
# =====================================================================


class TestCollectionModels:
    """Validate that all 14 collection models are defined correctly."""

    def test_all_14_models_defined(self):
        """COLLECTION_MODELS contains exactly 14 entries."""
        assert len(COLLECTION_MODELS) == 14

    @pytest.mark.parametrize("name", ALL_14_COLLECTIONS)
    def test_model_exists(self, name):
        """Each expected collection name is present in COLLECTION_MODELS."""
        assert name in COLLECTION_MODELS

    def test_schemas_and_models_have_matching_keys(self):
        """COLLECTION_SCHEMAS and COLLECTION_MODELS have identical key sets."""
        assert set(COLLECTION_SCHEMAS.keys()) == set(COLLECTION_MODELS.keys())

    def test_genomic_evidence_model_is_none(self):
        """genomic_evidence is read-only and mapped to None."""
        assert COLLECTION_MODELS["genomic_evidence"] is None

    @pytest.mark.parametrize(
        "name",
        [n for n in ALL_14_COLLECTIONS if n != "genomic_evidence"],
    )
    def test_non_genomic_models_are_not_none(self, name):
        """All non-genomic collections have a non-None model class."""
        assert COLLECTION_MODELS[name] is not None


# =====================================================================
# BIOMARKER COLLECTION MANAGER TESTS
# =====================================================================


class TestBiomarkerCollectionManagerInit:
    """Test BiomarkerCollectionManager initialisation with mocked pymilvus."""

    def test_default_init(self):
        """Manager initialises with default host and port."""
        manager = BiomarkerCollectionManager()
        assert manager.host is not None
        assert manager.port is not None
        assert manager.embedding_dim == EMBEDDING_DIM

    def test_custom_host_port(self):
        """Manager accepts custom host and port."""
        manager = BiomarkerCollectionManager(host="10.0.0.1", port=9999)
        assert manager.host == "10.0.0.1"
        assert manager.port == 9999

    @patch("src.collections.connections")
    def test_connect_calls_pymilvus(self, mock_connections):
        """connect() calls pymilvus connections.connect."""
        manager = BiomarkerCollectionManager()
        manager.connect()
        mock_connections.connect.assert_called_once()

    @patch("src.collections.connections")
    def test_disconnect_calls_pymilvus(self, mock_connections):
        """disconnect() calls pymilvus connections.disconnect."""
        manager = BiomarkerCollectionManager()
        manager.disconnect()
        mock_connections.disconnect.assert_called_once_with("default")


class TestGetCollectionStats:
    """Test get_collection_stats with mocked pymilvus."""

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    def test_returns_dict_with_counts(self, mock_utility, mock_collection_cls):
        """get_collection_stats returns a dict mapping names to entity counts."""
        mock_utility.has_collection.return_value = True
        mock_coll_instance = MagicMock()
        mock_coll_instance.num_entities = 100
        mock_collection_cls.return_value = mock_coll_instance

        manager = BiomarkerCollectionManager()
        stats = manager.get_collection_stats()

        assert isinstance(stats, dict)
        assert len(stats) == 14
        for name in ALL_14_COLLECTIONS:
            assert name in stats
            assert stats[name] == 100

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    def test_missing_collection_shows_zero(self, mock_utility, mock_collection_cls):
        """Collections that do not exist show 0 count."""
        mock_utility.has_collection.return_value = False

        manager = BiomarkerCollectionManager()
        stats = manager.get_collection_stats()

        assert isinstance(stats, dict)
        for name in ALL_14_COLLECTIONS:
            assert stats[name] == 0


class TestSearchAll:
    """Test search_all with mocked pymilvus."""

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    def test_search_all_returns_dict_for_all_collections(
        self, mock_utility, mock_collection_cls
    ):
        """search_all returns results dict with an entry per collection."""
        mock_utility.has_collection.return_value = True
        mock_coll_instance = MagicMock()
        # pymilvus search returns a list of lists of hits
        mock_coll_instance.search.return_value = []
        mock_collection_cls.return_value = mock_coll_instance

        manager = BiomarkerCollectionManager()
        # Pre-populate _collections so get_collection doesn't try to create
        for name in ALL_14_COLLECTIONS:
            manager._collections[name] = mock_coll_instance

        query_vec = [0.0] * EMBEDDING_DIM
        results = manager.search_all(query_embedding=query_vec, top_k_per_collection=3)

        assert isinstance(results, dict)
        for name in ALL_14_COLLECTIONS:
            assert name in results


class TestInsertBatch:
    """Test insert_batch with mocked pymilvus."""

    def test_insert_batch_empty_records(self):
        """insert_batch with an empty list returns 0 without error."""
        manager = BiomarkerCollectionManager()
        mock_coll = MagicMock()
        mock_result = MagicMock()
        mock_result.insert_count = 0
        mock_coll.insert.return_value = mock_result
        manager._collections["biomarker_reference"] = mock_coll

        count = manager.insert_batch("biomarker_reference", [])
        assert count == 0
        mock_coll.insert.assert_called_once_with([])
        mock_coll.flush.assert_called_once()

    def test_insert_batch_with_records(self):
        """insert_batch inserts records and returns the count."""
        manager = BiomarkerCollectionManager()
        mock_coll = MagicMock()
        mock_result = MagicMock()
        mock_result.insert_count = 3
        mock_coll.insert.return_value = mock_result
        manager._collections["biomarker_reference"] = mock_coll

        records = [{"id": f"r{i}", "embedding": [0.0] * 384} for i in range(3)]
        count = manager.insert_batch("biomarker_reference", records)
        assert count == 3

    @patch("src.collections.utility")
    def test_insert_batch_unknown_collection_raises(self, mock_utility):
        """insert_batch raises ValueError for an unknown collection name."""
        mock_utility.has_collection.return_value = False
        manager = BiomarkerCollectionManager()
        with pytest.raises(ValueError, match="Unknown collection"):
            manager.insert_batch("nonexistent_collection", [])
