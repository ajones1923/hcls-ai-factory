"""
Tests for src/collections.py — AutoimmuneCollectionManager

All external dependencies (pymilvus) are mocked so tests run without Milvus.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest
from pymilvus import DataType, MilvusException

# We need to mock pymilvus before importing the module under test in some cases,
# but since collections.py imports pymilvus at module level and defines schemas
# using real pymilvus types, we import normally and mock at call sites.
from src.collections import (
    _DIM,
    COLLECTION_SCHEMAS,
    INDEX_PARAMS,
    SEARCH_PARAMS,
    AutoimmuneCollectionManager,
    _embedding,
    _float,
    _int,
    _pk,
    _varchar,
)

# ── Schema / helper tests ────────────────────────────────────────────────


class TestSchemaHelpers:
    """Test the field-schema helper functions."""

    def test_pk_defaults(self):
        f = _pk()
        assert f.name == "id"
        assert f.dtype == DataType.VARCHAR
        assert f.is_primary is True
        assert f.max_length == 128

    def test_pk_custom_length(self):
        f = _pk(max_len=256)
        assert f.max_length == 256

    def test_embedding_defaults(self):
        f = _embedding()
        assert f.name == "embedding"
        assert f.dtype == DataType.FLOAT_VECTOR
        assert f.dim == 384

    def test_embedding_custom_dim(self):
        f = _embedding(dim=768)
        assert f.dim == 768

    def test_varchar(self):
        f = _varchar("test_field", max_length=512)
        assert f.name == "test_field"
        assert f.dtype == DataType.VARCHAR
        assert f.max_length == 512

    def test_varchar_default_length(self):
        f = _varchar("content")
        assert f.max_length == 3000

    def test_int_field(self):
        f = _int("page_number")
        assert f.name == "page_number"
        assert f.dtype == DataType.INT64

    def test_float_field(self):
        f = _float("score")
        assert f.name == "score"
        assert f.dtype == DataType.FLOAT


class TestCollectionSchemas:
    """Verify that all 14 schemas are registered with correct structure."""

    EXPECTED_COLLECTIONS = [
        "autoimmune_clinical_documents",
        "autoimmune_patient_labs",
        "autoimmune_autoantibody_panels",
        "autoimmune_hla_associations",
        "autoimmune_disease_criteria",
        "autoimmune_disease_activity",
        "autoimmune_flare_patterns",
        "autoimmune_biologic_therapies",
        "autoimmune_pgx_rules",
        "autoimmune_clinical_trials",
        "autoimmune_literature",
        "autoimmune_patient_timelines",
        "autoimmune_cross_disease",
    ]

    def test_all_13_autoimmune_schemas_registered(self):
        """13 autoimmune schemas (genomic_evidence is not registered here)."""
        for name in self.EXPECTED_COLLECTIONS:
            assert name in COLLECTION_SCHEMAS, f"Missing schema: {name}"

    def test_schema_count(self):
        assert len(COLLECTION_SCHEMAS) == 13

    @pytest.mark.parametrize("name", EXPECTED_COLLECTIONS)
    def test_each_schema_has_id_and_embedding(self, name):
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f.name for f in schema.fields]
        assert "id" in field_names
        assert "embedding" in field_names

    @pytest.mark.parametrize("name", EXPECTED_COLLECTIONS)
    def test_each_schema_has_text_chunk(self, name):
        schema = COLLECTION_SCHEMAS[name]
        field_names = [f.name for f in schema.fields]
        assert "text_chunk" in field_names

    def test_clinical_documents_fields(self):
        schema = COLLECTION_SCHEMAS["autoimmune_clinical_documents"]
        field_names = [f.name for f in schema.fields]
        assert "patient_id" in field_names
        assert "doc_type" in field_names
        assert "specialty" in field_names
        assert "page_number" in field_names
        assert "chunk_index" in field_names

    def test_patient_labs_has_value_field(self):
        schema = COLLECTION_SCHEMAS["autoimmune_patient_labs"]
        field_map = {f.name: f for f in schema.fields}
        assert "value" in field_map
        assert field_map["value"].dtype == DataType.FLOAT

    def test_hla_associations_has_odds_ratio(self):
        schema = COLLECTION_SCHEMAS["autoimmune_hla_associations"]
        field_map = {f.name: f for f in schema.fields}
        assert "odds_ratio" in field_map
        assert field_map["odds_ratio"].dtype == DataType.FLOAT

    def test_clinical_trials_pk_length_256(self):
        schema = COLLECTION_SCHEMAS["autoimmune_clinical_trials"]
        pk = next(f for f in schema.fields if f.is_primary)
        assert pk.max_length == 256

    def test_literature_pk_length_256(self):
        schema = COLLECTION_SCHEMAS["autoimmune_literature"]
        pk = next(f for f in schema.fields if f.is_primary)
        assert pk.max_length == 256

    def test_cross_disease_has_co_occurrence_rate(self):
        schema = COLLECTION_SCHEMAS["autoimmune_cross_disease"]
        field_map = {f.name: f for f in schema.fields}
        assert "co_occurrence_rate" in field_map
        assert field_map["co_occurrence_rate"].dtype == DataType.FLOAT

    def test_embedding_dimension(self):
        """All embedding fields use 384 dimensions (BGE-small-en-v1.5)."""
        for name, schema in COLLECTION_SCHEMAS.items():
            emb = next(f for f in schema.fields if f.name == "embedding")
            assert emb.dim == _DIM, f"{name} embedding dim={emb.dim}, expected {_DIM}"


class TestIndexAndSearchParams:
    """Verify index/search parameter constants."""

    def test_index_params_metric(self):
        assert INDEX_PARAMS["metric_type"] == "COSINE"

    def test_index_params_type(self):
        assert INDEX_PARAMS["index_type"] == "IVF_FLAT"

    def test_index_params_nlist(self):
        assert INDEX_PARAMS["params"]["nlist"] == 1024

    def test_search_params_metric(self):
        assert SEARCH_PARAMS["metric_type"] == "COSINE"

    def test_search_params_nprobe(self):
        assert SEARCH_PARAMS["params"]["nprobe"] == 16


# ── CollectionManager tests ──────────────────────────────────────────────


class TestAutoimmuneCollectionManagerInit:
    """Test manager initialization and configuration."""

    def test_default_init(self):
        mgr = AutoimmuneCollectionManager()
        assert mgr.host == "localhost"
        assert mgr.port == 19530
        assert mgr.embedding_dim == 384
        assert mgr._alias == "autoimmune_agent"
        assert mgr._connected is False

    def test_custom_init(self):
        mgr = AutoimmuneCollectionManager(host="10.0.0.1", port=19531, embedding_dim=768)
        assert mgr.host == "10.0.0.1"
        assert mgr.port == 19531
        assert mgr.embedding_dim == 768

    def test_is_connected_property(self):
        mgr = AutoimmuneCollectionManager()
        assert mgr.is_connected is False
        mgr._connected = True
        assert mgr.is_connected is True


class TestCollectionManagerConnect:
    """Test connect/disconnect with mocked pymilvus."""

    @patch("src.collections.connections")
    def test_connect_success(self, mock_connections):
        mgr = AutoimmuneCollectionManager(host="testhost", port=19530)
        mgr.connect()
        mock_connections.connect.assert_called_once_with(
            alias="autoimmune_agent", host="testhost", port=19530,
        )
        assert mgr._connected is True

    @patch("src.collections.connections")
    def test_connect_already_connected_is_noop(self, mock_connections):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.connect()
        mock_connections.connect.assert_not_called()

    @patch("src.collections.connections")
    def test_connect_failure_raises(self, mock_connections):
        mock_connections.connect.side_effect = MilvusException(message="Connection refused")
        mgr = AutoimmuneCollectionManager()
        with pytest.raises(MilvusException):
            mgr.connect()
        assert mgr._connected is False

    @patch("src.collections.connections")
    def test_disconnect(self, mock_connections):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.disconnect()
        mock_connections.disconnect.assert_called_once_with("autoimmune_agent")
        assert mgr._connected is False

    @patch("src.collections.connections")
    def test_disconnect_when_not_connected_is_noop(self, mock_connections):
        mgr = AutoimmuneCollectionManager()
        mgr.disconnect()
        mock_connections.disconnect.assert_not_called()

    @patch("src.collections.connections")
    def test_ensure_connected_calls_connect(self, mock_connections):
        mgr = AutoimmuneCollectionManager()
        mgr._ensure_connected()
        mock_connections.connect.assert_called_once()
        assert mgr._connected is True

    @patch("src.collections.connections")
    def test_ensure_connected_skips_if_connected(self, mock_connections):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr._ensure_connected()
        mock_connections.connect.assert_not_called()


class TestCollectionManagerListCollections:
    """Test list_collections filtering."""

    @patch("src.collections.connections")
    @patch("src.collections.utility")
    def test_list_filters_autoimmune_and_genomic(self, mock_utility, mock_conn):
        mock_utility.list_collections.return_value = [
            "autoimmune_patient_labs",
            "genomic_evidence",
            "biomarker_cancer_variants",
            "autoimmune_hla_associations",
            "some_other_collection",
        ]
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        result = mgr.list_collections()
        assert "autoimmune_patient_labs" in result
        assert "autoimmune_hla_associations" in result
        assert "genomic_evidence" in result
        assert "biomarker_cancer_variants" not in result
        assert "some_other_collection" not in result

    @patch("src.collections.connections")
    @patch("src.collections.utility")
    def test_list_empty(self, mock_utility, mock_conn):
        mock_utility.list_collections.return_value = []
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        assert mgr.list_collections() == []


class TestCollectionManagerCreateCollection:
    """Test create_collection with mocked pymilvus."""

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_create_new_collection(self, mock_conn, mock_utility, mock_coll_cls):
        mock_utility.has_collection.return_value = False
        mock_coll = MagicMock()
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        result = mgr.create_collection("autoimmune_patient_labs")

        mock_coll.create_index.assert_called_once_with("embedding", INDEX_PARAMS)
        mock_coll.load.assert_called_once()
        assert result is mock_coll

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_create_existing_collection_no_drop(self, mock_conn, mock_utility, mock_coll_cls):
        mock_utility.has_collection.return_value = True
        mock_coll = MagicMock()
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.create_collection("autoimmune_patient_labs", drop_existing=False)

        mock_utility.drop_collection.assert_not_called()
        mock_coll.load.assert_called_once()

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_create_existing_collection_with_drop(self, mock_conn, mock_utility, mock_coll_cls):
        mock_utility.has_collection.return_value = True
        mock_coll = MagicMock()
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.create_collection("autoimmune_patient_labs", drop_existing=True)

        mock_utility.drop_collection.assert_called_once_with(
            "autoimmune_patient_labs", using="autoimmune_agent",
        )

    def test_create_unknown_collection_raises(self):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        with pytest.raises(ValueError, match="Unknown collection"):
            mgr.create_collection("nonexistent_collection")


class TestCollectionManagerInsert:
    """Test insert with validation and edge cases."""

    def _make_manager(self):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        return mgr

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_insert_empty_records_returns_zero(self, mock_conn, mock_coll_cls):
        mgr = self._make_manager()
        assert mgr.insert("autoimmune_patient_labs", []) == 0

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_insert_skips_records_without_embedding(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        mock_schema_field_id = MagicMock()
        mock_schema_field_id.name = "id"
        mock_schema_field_id.dtype = DataType.VARCHAR

        mock_schema_field_emb = MagicMock()
        mock_schema_field_emb.name = "embedding"
        mock_schema_field_emb.dtype = DataType.FLOAT_VECTOR

        mock_schema_field_text = MagicMock()
        mock_schema_field_text.name = "text_chunk"
        mock_schema_field_text.dtype = DataType.VARCHAR

        mock_coll.schema.fields = [mock_schema_field_id, mock_schema_field_emb, mock_schema_field_text]
        mock_coll_cls.return_value = mock_coll

        mgr = self._make_manager()
        # Record with no embedding
        records = [{"id": "r1", "text_chunk": "hello"}]
        result = mgr.insert("autoimmune_patient_labs", records)
        assert result == 0

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_insert_skips_wrong_dimension_embedding(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        mock_schema_field_id = MagicMock()
        mock_schema_field_id.name = "id"
        mock_schema_field_id.dtype = DataType.VARCHAR

        mock_schema_field_emb = MagicMock()
        mock_schema_field_emb.name = "embedding"
        mock_schema_field_emb.dtype = DataType.FLOAT_VECTOR

        mock_coll.schema.fields = [mock_schema_field_id, mock_schema_field_emb]
        mock_coll_cls.return_value = mock_coll

        mgr = self._make_manager()
        # Wrong dimension (10 instead of 384)
        records = [{"id": "r1", "embedding": [0.1] * 10}]
        result = mgr.insert("autoimmune_patient_labs", records)
        assert result == 0

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_insert_valid_records(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        # Minimal schema: id, embedding, text_chunk
        fields = []
        for fname, dtype in [("id", DataType.VARCHAR), ("embedding", DataType.FLOAT_VECTOR), ("text_chunk", DataType.VARCHAR)]:
            f = MagicMock()
            f.name = fname
            f.dtype = dtype
            f.max_length = 3000 if dtype == DataType.VARCHAR else None
            fields.append(f)
        mock_coll.schema.fields = fields

        mock_insert_result = MagicMock()
        mock_insert_result.insert_count = 2
        mock_coll.insert.return_value = mock_insert_result
        mock_coll_cls.return_value = mock_coll

        mgr = self._make_manager()
        records = [
            {"id": "r1", "embedding": [0.1] * 384, "text_chunk": "Record 1"},
            {"id": "r2", "embedding": [0.2] * 384, "text_chunk": "Record 2"},
        ]
        result = mgr.insert("autoimmune_patient_labs", records)
        assert result == 2
        mock_coll.flush.assert_called_once()

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_insert_milvus_exception_raises(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        fields = []
        for fname, dtype in [("id", DataType.VARCHAR), ("embedding", DataType.FLOAT_VECTOR)]:
            f = MagicMock()
            f.name = fname
            f.dtype = dtype
            f.max_length = 128 if dtype == DataType.VARCHAR else None
            fields.append(f)
        mock_coll.schema.fields = fields
        mock_coll.insert.side_effect = MilvusException(message="Insert failed")
        mock_coll_cls.return_value = mock_coll

        mgr = self._make_manager()
        records = [{"id": "r1", "embedding": [0.1] * 384}]
        with pytest.raises(MilvusException):
            mgr.insert("autoimmune_patient_labs", records)

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_insert_truncates_long_strings(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        f_id = MagicMock()
        f_id.name = "id"
        f_id.dtype = DataType.VARCHAR
        f_id.max_length = 128

        f_emb = MagicMock()
        f_emb.name = "embedding"
        f_emb.dtype = DataType.FLOAT_VECTOR
        f_emb.max_length = None

        f_text = MagicMock()
        f_text.name = "text_chunk"
        f_text.dtype = DataType.VARCHAR
        f_text.max_length = 10  # Very short for testing truncation

        mock_coll.schema.fields = [f_id, f_emb, f_text]

        mock_result = MagicMock()
        mock_result.insert_count = 1
        mock_coll.insert.return_value = mock_result
        mock_coll_cls.return_value = mock_coll

        mgr = self._make_manager()
        records = [{"id": "r1", "embedding": [0.1] * 384, "text_chunk": "A" * 100}]
        mgr.insert("autoimmune_patient_labs", records)

        # Verify the data passed to insert had the text truncated
        call_args = mock_coll.insert.call_args[0][0]
        text_data = call_args[2]  # third field is text_chunk
        assert len(text_data[0]) == 10


class TestCollectionManagerSearch:
    """Test search and search_all with mocked Milvus."""

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_search_returns_hits(self, mock_conn, mock_coll_cls):
        mock_hit = MagicMock()
        mock_hit.id = "hit1"
        mock_hit.score = 0.95
        mock_hit.entity.get.side_effect = lambda field, default="": {
            "text_chunk": "Sample text",
        }.get(field, default)

        mock_coll = MagicMock()
        mock_coll.schema.fields = [
            MagicMock(name="id"), MagicMock(name="embedding"), MagicMock(name="text_chunk"),
        ]
        mock_coll.search.return_value = [[mock_hit]]
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        results = mgr.search(
            "autoimmune_patient_labs",
            query_embedding=[0.1] * 384,
            top_k=5,
        )
        assert len(results) == 1
        assert results[0]["id"] == "hit1"
        assert results[0]["score"] == 0.95

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_search_with_filter_expr(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        mock_coll.schema.fields = [MagicMock(name="id"), MagicMock(name="embedding")]
        mock_coll.search.return_value = [[]]
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.search(
            "autoimmune_patient_labs",
            query_embedding=[0.1] * 384,
            filter_expr='patient_id == "P001"',
        )
        call_kwargs = mock_coll.search.call_args[1]
        assert call_kwargs["expr"] == 'patient_id == "P001"'

    @patch("src.collections.Collection")
    @patch("src.collections.connections")
    def test_search_exception_returns_empty(self, mock_conn, mock_coll_cls):
        mock_coll = MagicMock()
        mock_coll.schema.fields = [MagicMock(name="id"), MagicMock(name="embedding")]
        mock_coll.search.side_effect = Exception("Search failed")
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        results = mgr.search("autoimmune_patient_labs", [0.1] * 384)
        assert results == []

    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_search_all_parallel(self, mock_conn, mock_utility):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True

        mock_utility.list_collections.return_value = ["autoimmune_patient_labs", "autoimmune_hla_associations"]

        # Mock search method to return results per collection
        def fake_search(coll_name, emb, top_k, filter_expr=None, output_fields=None):
            return [{"id": f"{coll_name}_1", "score": 0.85, "text_chunk": "text"}]

        mgr.search = MagicMock(side_effect=fake_search)

        results = mgr.search_all(
            query_embedding=[0.1] * 384,
            top_k_per_collection=3,
        )
        assert len(results) == 2
        assert "autoimmune_patient_labs" in results
        assert "autoimmune_hla_associations" in results

    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_search_all_with_score_threshold(self, mock_conn, mock_utility):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mock_utility.list_collections.return_value = ["autoimmune_patient_labs"]

        def fake_search(coll_name, emb, top_k, filter_expr=None, output_fields=None):
            return [
                {"id": "high", "score": 0.9},
                {"id": "low", "score": 0.3},
            ]

        mgr.search = MagicMock(side_effect=fake_search)

        results = mgr.search_all(
            query_embedding=[0.1] * 384,
            score_threshold=0.5,
        )
        # Only the high-score hit should survive
        hits = results.get("autoimmune_patient_labs", [])
        assert len(hits) == 1
        assert hits[0]["id"] == "high"

    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_search_all_specific_collections(self, mock_conn, mock_utility):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True

        mgr.search = MagicMock(return_value=[{"id": "x", "score": 0.8}])
        results = mgr.search_all(
            query_embedding=[0.1] * 384,
            collections=["autoimmune_hla_associations"],
        )
        # Should only search the specified collection
        mgr.search.assert_called_once()
        assert "autoimmune_hla_associations" in results


class TestCollectionManagerStats:
    """Test get_collection_stats and get_collection_count."""

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_get_collection_count(self, mock_conn, mock_utility, mock_coll_cls):
        mock_coll = MagicMock()
        mock_coll.num_entities = 42
        mock_coll_cls.return_value = mock_coll

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        assert mgr.get_collection_count("autoimmune_patient_labs") == 42

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_get_collection_count_error_returns_zero(self, mock_conn, mock_utility, mock_coll_cls):
        mock_coll_cls.side_effect = Exception("not found")

        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        assert mgr.get_collection_count("nonexistent") == 0

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_get_collection_stats(self, mock_conn, mock_utility, mock_coll_cls):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.list_collections = MagicMock(return_value=["autoimmune_patient_labs"])

        mock_coll = MagicMock()
        mock_coll.num_entities = 100
        mock_coll_cls.return_value = mock_coll

        stats = mgr.get_collection_stats()
        assert stats["autoimmune_patient_labs"] == 100

    @patch("src.collections.Collection")
    @patch("src.collections.utility")
    @patch("src.collections.connections")
    def test_get_collection_stats_error_returns_minus_one(self, mock_conn, mock_utility, mock_coll_cls):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True
        mgr.list_collections = MagicMock(return_value=["broken_coll"])

        mock_coll_cls.side_effect = Exception("fail")

        stats = mgr.get_collection_stats()
        assert stats["broken_coll"] == -1


class TestInsertBatch:
    """Test insert_batch divides records into batches."""

    @patch("src.collections.connections")
    def test_insert_batch_calls_insert_multiple_times(self, mock_conn):
        mgr = AutoimmuneCollectionManager()
        mgr._connected = True

        call_count = 0
        def fake_insert(coll_name, records):
            nonlocal call_count
            call_count += 1
            return len(records)

        mgr.insert = MagicMock(side_effect=fake_insert)
        records = [{"id": f"r{i}", "embedding": [0.1] * 384} for i in range(250)]
        total = mgr.insert_batch("autoimmune_patient_labs", records, batch_size=100)

        assert total == 250
        assert mgr.insert.call_count == 3  # 100 + 100 + 50
