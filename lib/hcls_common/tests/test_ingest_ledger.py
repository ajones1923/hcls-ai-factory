"""The provenance ledger, and the model conflict it exists to catch.

A mixed-embedding collection is this platform's signature failure mode: nothing raises, the
insert succeeds, and retrieval returns confident nonsense. Milvus validates vector *dimension*,
so a 768-dim model into a 384-dim collection fails loudly — but two different 384-dim models
insert equally happily into one collection, and a nearest-neighbour search across two embedding
spaces has no error to report. Dimension checking cannot see it; a recorded model can.
"""
import json

import pytest

from hcls_common.ingest_ledger import (
    IngestEvent,
    append,
    history,
    model_conflict,
    models_used,
    read,
)


@pytest.fixture
def ledger(tmp_path):
    return tmp_path / "prov.jsonl"


def _ev(collection="rd_literature", model="BAAI/bge-small-en-v1.5", written=10, **kw):
    return IngestEvent(collection=collection, source="pubmed", query="q",
                       embed_model=model, fetched=written, written=written, **kw)


class TestRoundTrip:
    def test_a_missing_ledger_is_empty_not_an_error(self, ledger):
        assert read(ledger) == []

    def test_append_then_read(self, ledger):
        append(_ev(), ledger)
        append(_ev(collection="trial_results"), ledger)
        got = read(ledger)
        assert [e.collection for e in got] == ["rd_literature", "trial_results"]
        assert got[0].source == "pubmed"

    def test_it_creates_its_parent_directory(self, tmp_path):
        nested = tmp_path / "docs" / "build" / "prov.jsonl"
        append(_ev(), nested)
        assert nested.exists() and len(read(nested)) == 1

    def test_every_event_carries_a_timestamp(self, ledger):
        append(_ev(), ledger)
        assert read(ledger)[0].at.endswith("+00:00")


class TestItSurvivesADamagedFile:
    def test_one_corrupt_line_does_not_hide_the_rest(self, ledger):
        append(_ev(collection="a"), ledger)
        with ledger.open("a") as fh:
            fh.write("{not json at all\n")
        append(_ev(collection="b"), ledger)
        assert [e.collection for e in read(ledger)] == ["a", "b"]

    def test_an_unknown_field_from_a_future_version_is_ignored(self, ledger):
        """Forward compatibility: a newer writer's extra key must not break an older reader."""
        with ledger.open("a") as fh:
            fh.write(json.dumps({
                "collection": "a", "source": "pubmed", "query": "q",
                "embed_model": "m", "fetched": 1, "written": 1,
                "something_added_later": True,
            }) + "\n")
        assert [e.collection for e in read(ledger)] == ["a"]

    def test_a_line_missing_a_required_field_is_skipped(self, ledger):
        with ledger.open("a") as fh:
            fh.write(json.dumps({"collection": "a"}) + "\n")
        append(_ev(collection="good"), ledger)
        assert [e.collection for e in read(ledger)] == ["good"]


class TestModelConflict:
    def test_no_history_is_no_conflict(self, ledger):
        assert model_conflict("rd_literature", "BAAI/bge-small-en-v1.5", ledger) is None

    def test_the_same_model_again_is_no_conflict(self, ledger):
        append(_ev(), ledger)
        assert model_conflict("rd_literature", "BAAI/bge-small-en-v1.5", ledger) is None

    def test_a_different_model_is_a_conflict(self, ledger):
        append(_ev(), ledger)
        msg = model_conflict("rd_literature", "sentence-transformers/all-MiniLM-L6-v2", ledger)
        assert msg is not None
        assert "mix embedding spaces" in msg
        assert "rd_literature" in msg

    def test_conflict_names_the_model_already_present(self, ledger):
        append(_ev(model="model-A"), ledger)
        msg = model_conflict("rd_literature", "model-B", ledger)
        assert "model-A" in msg

    def test_another_collection_s_model_is_irrelevant(self, ledger):
        """Collections are independent embedding spaces; only the target one matters."""
        append(_ev(collection="other", model="model-A"), ledger)
        assert model_conflict("rd_literature", "model-B", ledger) is None

    def test_a_run_that_wrote_nothing_does_not_pin_the_model(self, ledger):
        """A failed ingest must not lock a collection to a model it never actually stored."""
        append(_ev(model="model-A", written=0), ledger)
        assert models_used("rd_literature", ledger) == []
        assert model_conflict("rd_literature", "model-B", ledger) is None


class TestOperatorSummary:
    def test_history_totals_rows_per_collection(self, ledger):
        append(_ev(collection="a", written=10), ledger)
        append(_ev(collection="a", written=5), ledger)
        append(_ev(collection="b", written=7), ledger)
        assert list(history(ledger)) == [
            ("a", 15, "BAAI/bge-small-en-v1.5"),
            ("b", 7, "BAAI/bge-small-en-v1.5"),
        ]

    def test_the_ledger_holds_no_record_text(self, ledger):
        """Abstracts stay on the box under NCBI's terms; the ledger is metadata only."""
        append(_ev(), ledger)
        fields = set(json.loads(ledger.read_text().splitlines()[0]))
        assert not fields & {"abstract", "text", "title", "content"}
