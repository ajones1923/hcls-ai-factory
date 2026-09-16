"""Tests for the one shared Milvus similarity search.

Every clinical answer in the factory passes through this function, once per collection searched,
and five subjects carried it in copy. Its failure mode is an empty list — indistinguishable from
a corpus with no match — which is exactly why it needs tests rather than eyeballs.
"""
from hcls_common.vector_search import search_collection


class _Client:
    """Records the kwargs it was called with and returns a canned result set."""
    def __init__(self, results=None, raises=None):
        self.results, self.raises, self.last = results, raises, None

    def search(self, **kw):
        self.last = kw
        if self.raises:
            raise self.raises
        return self.results


class _Hit:                      # the ORM API's shape
    def __init__(self, id, score, entity):
        self.id, self.score, self.entity = id, score, entity


class TestRequestShape:
    def test_uses_search_params_not_the_orm_param(self):
        """MilvusClient takes `search_params`; passing the ORM's `param` raises inside the try
        and returns [] — a silent empty result set rather than an error."""
        c = _Client([[]])
        search_collection(c, "trials", [0.1, 0.2], 5)
        assert "search_params" in c.last and "param" not in c.last
        assert c.last["search_params"]["metric_type"] == "COSINE"
        assert c.last["collection_name"] == "trials"
        assert c.last["limit"] == 5
        assert c.last["data"] == [[0.1, 0.2]]

    def test_filter_is_sent_only_when_given(self):
        c = _Client([[]])
        search_collection(c, "trials", [0.1], 5)
        assert "filter" not in c.last
        search_collection(c, "trials", [0.1], 5, 'phase == "Phase 3"')
        assert c.last["filter"] == 'phase == "Phase 3"'

    def test_metric_and_nprobe_are_overridable(self):
        c = _Client([[]])
        search_collection(c, "t", [0.1], 5, metric_type="L2", nprobe=64)
        assert c.last["search_params"] == {"metric_type": "L2", "params": {"nprobe": 64}}


class TestMilvusClientShape:
    """Plain dicts: {id, distance, entity}."""

    def test_flattens_dict_hits(self):
        c = _Client([[{"id": 42, "distance": 0.91,
                       "entity": {"title": "ELIANA", "nct_id": "NCT02435849",
                                  "embedding": [0.0] * 384}}]])
        (r,) = search_collection(c, "trials", [0.1], 1)
        assert r["id"] == "42" and r["score"] == 0.91
        assert r["title"] == "ELIANA" and r["nct_id"] == "NCT02435849"

    def test_embedding_is_never_returned(self):
        c = _Client([[{"id": 1, "distance": 0.5, "entity": {"embedding": [0.0] * 384, "t": "x"}}]])
        (r,) = search_collection(c, "trials", [0.1], 1)
        assert "embedding" not in r and "embedding" not in r["metadata"]

    def test_metadata_excludes_the_reserved_keys(self):
        c = _Client([[{"id": 1, "distance": 0.5, "entity": {"title": "T"}}]])
        (r,) = search_collection(c, "trials", [0.1], 1)
        assert r["metadata"] == {"title": "T"}

    def test_score_falls_back_to_score_key_and_none(self):
        assert search_collection(_Client([[{"id": 1, "score": 0.7, "entity": {}}]]),
                                 "c", [0.1], 1)[0]["score"] == 0.7
        assert search_collection(_Client([[{"id": 1, "distance": None, "entity": {}}]]),
                                 "c", [0.1], 1)[0]["score"] == 0.0

    def test_missing_entity_is_tolerated(self):
        (r,) = search_collection(_Client([[{"id": 1, "distance": 0.5}]]), "c", [0.1], 1)
        assert r["id"] == "1"


class TestOrmShape:
    def test_flattens_hit_objects_with_dict_entity(self):
        c = _Client([[_Hit(7, 0.88, {"title": "ZUMA-1", "embedding": [0.0]})]])
        (r,) = search_collection(c, "trials", [0.1], 1)
        assert r["id"] == "7" and r["score"] == 0.88 and r["title"] == "ZUMA-1"
        assert "embedding" not in r

    def test_flattens_hit_objects_with_fields_entity(self):
        class _Ent:
            fields = {"title": "JULIET", "embedding": [0.0]}
        c = _Client([[_Hit(9, 0.5, _Ent())]])
        (r,) = search_collection(c, "trials", [0.1], 1)
        assert r["title"] == "JULIET" and "embedding" not in r


class TestDegradation:
    """A retrieval failure must degrade the answer, never break the request."""

    def test_exception_returns_empty_and_warns(self, caplog):
        c = _Client(raises=RuntimeError("collection not loaded"))
        with caplog.at_level("WARNING"):
            assert search_collection(c, "trials", [0.1], 5) == []
        assert "trials" in caplog.text

    def test_empty_and_none_results_are_empty_lists(self):
        assert search_collection(_Client([]), "c", [0.1], 5) == []
        assert search_collection(_Client(None), "c", [0.1], 5) == []
        assert search_collection(_Client([[]]), "c", [0.1], 5) == []
