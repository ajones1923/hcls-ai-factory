"""Retrieval coverage must be reported, not implied.

This engine is the extreme case. All twelve `cardio_*` collections were ABSENT
from Milvus on 2026-11-03 -- never created, not merely empty -- so retrieval could
only ever reach the shared `genomic_evidence`: **2% of the declared weight**, with
`collections_searched` reporting 13. `scripts/check_corpus.py` could not show this
at all, because it iterates collections that exist.

The risk calculators are pure functions over their inputs and do not touch Milvus;
this is about the retrieved half of the engine's output only.

A search of an empty or missing collection succeeds and returns zero rows. Nothing in the
retrieval path can tell that apart from a search that found nothing relevant, so the response has
to carry the difference itself. See `hcls_common.corpus_coverage` and
docs/build/EMPTY_COLLECTIONS.md.
"""
import pytest

from src.agent import CardioWorkflowType
from src.rag_engine import CardioRAGEngine

# Measured 2026-11-03: every cardio_* collection was absent from Milvus. The shared
# genomic_evidence collection, which this engine also searches, held 60 rows.
POPULATED = {"genomic_evidence"}


class FakeEmbedder:
    def embed_text(self, text):
        return [0.1] * 384


@pytest.fixture
def engine(monkeypatch):
    eng = CardioRAGEngine(milvus_client=object(), embedding_model=FakeEmbedder(),
                   llm_client=None, session_id="coverage-test")

    def fake_search(collection_name, query_vector, top_k, filter_expr=None):
        if collection_name not in eng._populated:
            return []
        return [{"id": f"{collection_name}-1", "score": 0.8,
                 "text": f"a row from {collection_name}"}]

    eng._populated = set(POPULATED)
    monkeypatch.setattr(eng, "_search_collection", fake_search)
    return eng


WORKFLOW = CardioWorkflowType.HEART_FAILURE


class TestThePartialCorpusIsReported:
    def test_the_verdict_is_attached(self, engine):
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        assert r.retrieval_coverage is not None
        assert r.retrieval_coverage["degraded"] is True

    def test_contributing_is_reported_alongside_searched(self, engine):
        """The number that moves with the corpus, next to the one that does not."""
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        assert r.collections_searched > r.collections_contributing
        assert r.collections_contributing == len(POPULATED & set(
            engine._get_boosted_weights(WORKFLOW)))

    def test_the_silent_sources_are_named_heaviest_first(self, engine):
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        silent = r.retrieval_coverage["silent_collections"]
        assert silent
        weights = engine._get_boosted_weights(WORKFLOW)
        assert weights[silent[0]] == max(weights[c] for c in silent)

    def test_the_answer_itself_discloses_it(self, engine):
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        assert "Partial corpus" in r.answer
        assert "about the corpus, not about the question" in r.answer

    def test_the_measured_two_percent(self, engine):
        """One shared collection out of thirteen declared sources."""
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        assert r.retrieval_coverage["coverage"] == pytest.approx(0.02, abs=0.015)
        assert r.collections_contributing == 1
        assert r.collections_searched >= 12

    def test_every_cardio_collection_is_named_silent(self, engine):
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        silent = r.retrieval_coverage["silent_collections"]
        assert all(c.startswith("cardio_") for c in silent)
        assert len(silent) == 12


class TestAFullCorpusStaysQuiet:
    def test_no_disclosure_when_everything_answers(self, engine):
        engine._populated = set(engine._get_boosted_weights(WORKFLOW))
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        assert r.retrieval_coverage["coverage"] == pytest.approx(1.0)
        assert r.retrieval_coverage["degraded"] is False
        assert r.retrieval_coverage["silent_collections"] == []
        assert "Partial corpus" not in r.answer


class TestTheEmptyCorpusCase:
    def test_nothing_anywhere_reports_zero_coverage(self, engine):
        engine._populated = set()
        r = engine.query("management of HFrEF with EF 28%", workflow=WORKFLOW)
        assert r.retrieval_coverage["coverage"] == 0.0
        assert r.collections_contributing == 0
        assert "Partial corpus" in r.answer
