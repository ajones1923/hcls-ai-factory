"""Projection rules for writing records into a collection whose schema you did not choose.

Every rule here was learned from an insert that failed, or -- worse -- from one that was
reported as succeeding. The auto_id group is the 2026-11-04 case: all twelve `cardio_*`
collections declare `auto_id=True`, this module supplied a content-derived primary key anyway,
and Milvus rejected every row. The seeder on top of it then reported "Total records: 392 |
Failed: 0", because its own error paths returned the record count.
"""
import sys
import types
from unittest.mock import MagicMock

import pytest

from hcls_common.ingest_persist import persist_records


def _schema(fields, auto_id=False):
    return {"auto_id": auto_id, "fields": fields}


def _f(name, type_name, **kw):
    """A describe_collection() field. `type` carries .name, as pymilvus' DataType does."""
    return {"name": name, "type": types.SimpleNamespace(name=type_name), **kw}


CARDIO_IMAGING = _schema([
    _f("id", "INT64", is_primary=True, auto_id=True),
    _f("embedding", "FLOAT_VECTOR"),
    _f("modality", "VARCHAR", params={"max_length": 64}),
    _f("finding", "VARCHAR", params={"max_length": 512}),
    _f("cross_modal_trigger", "BOOL"),
], auto_id=True)

LITERATURE = _schema([
    _f("id", "VARCHAR", is_primary=True, params={"max_length": 256}),
    _f("embedding", "FLOAT_VECTOR"),
    _f("text_chunk", "VARCHAR", params={"max_length": 3000}),
    _f("year", "INT64"),
])


@pytest.fixture
def milvus(monkeypatch):
    """Stand in for pymilvus + sentence-transformers, and record what gets written."""
    client = MagicMock()
    pymilvus = types.ModuleType("pymilvus")
    pymilvus.MilvusClient = MagicMock(return_value=client)

    st = types.ModuleType("sentence_transformers")
    encoder = MagicMock()
    encoder.encode.side_effect = lambda texts, **kw: types.SimpleNamespace(
        tolist=lambda: [[0.0] * 384 for _ in texts])
    st.SentenceTransformer = MagicMock(return_value=encoder)

    monkeypatch.setitem(sys.modules, "pymilvus", pymilvus)
    monkeypatch.setitem(sys.modules, "sentence_transformers", st)
    return client


def _record(text="echo — LV systolic function", **meta):
    return types.SimpleNamespace(text=text, metadata=meta, collection="cardio_imaging",
                                 source="curated", source_id=None)


def _written(client):
    """The single batch handed to Milvus, whichever call was used."""
    call = client.insert.call_args or client.upsert.call_args
    assert call is not None, "nothing was written"
    return call.kwargs["data"]


class TestAutoIdCollections:
    def test_the_primary_key_is_not_supplied(self, milvus):
        """Milvus rejects a supplied pk on an auto_id collection -- the whole batch."""
        milvus.describe_collection.return_value = CARDIO_IMAGING
        persist_records([_record(modality="echo")], "cardio_imaging")
        assert "id" not in _written(milvus)[0]

    def test_it_inserts_rather_than_upserts(self, milvus):
        """`upsert` needs a key to match on. There isn't one, so it cannot be used."""
        milvus.describe_collection.return_value = CARDIO_IMAGING
        persist_records([_record(modality="echo")], "cardio_imaging")
        assert milvus.insert.called
        assert not milvus.upsert.called

    def test_an_id_in_the_record_metadata_is_dropped_too(self, milvus):
        """A parser that volunteers an id must not be able to abort the batch."""
        milvus.describe_collection.return_value = CARDIO_IMAGING
        persist_records([_record(modality="echo", id=7)], "cardio_imaging")
        assert "id" not in _written(milvus)[0]

    def test_the_rest_of_the_projection_still_applies(self, milvus):
        """Declared columns present, BOOL left boolean, counts returned."""
        milvus.describe_collection.return_value = CARDIO_IMAGING
        n = persist_records([_record(modality="echo")], "cardio_imaging")
        row = _written(milvus)[0]
        assert n == 1
        assert row["modality"] == "echo"
        assert row["finding"] == ""                  # declared, absent, not nullable
        assert row["cross_modal_trigger"] is False   # BOOL, not ""
        assert len(row["embedding"]) == 384


class TestKeyedCollections:
    def test_a_content_derived_key_is_supplied_and_upserted(self, milvus):
        """Where the schema expects a key, re-ingest stays idempotent."""
        milvus.describe_collection.return_value = LITERATURE
        persist_records([_record(text="an abstract", year=2024)], "autoimmune_literature")
        row = _written(milvus)[0]
        assert row["id"].startswith("autoimmune_literature-")
        assert milvus.upsert.called and not milvus.insert.called

    def test_the_same_text_produces_the_same_key(self, milvus):
        milvus.describe_collection.return_value = LITERATURE
        keys = []
        for _ in range(2):
            milvus.reset_mock()
            persist_records([_record(text="an abstract")], "autoimmune_literature")
            keys.append(_written(milvus)[0]["id"])
        assert keys[0] == keys[1]


class TestItRefusesToClaimSuccess:
    def test_a_missing_collection_writes_nothing_and_returns_zero(self, milvus):
        milvus.describe_collection.side_effect = RuntimeError("collection not found")
        assert persist_records([_record()], "cardio_conditions") == 0
        assert not milvus.insert.called and not milvus.upsert.called

    def test_a_failed_insert_returns_zero(self, milvus):
        milvus.describe_collection.return_value = CARDIO_IMAGING
        milvus.insert.side_effect = RuntimeError("length of varchar field exceeds max length")
        assert persist_records([_record()], "cardio_imaging") == 0
