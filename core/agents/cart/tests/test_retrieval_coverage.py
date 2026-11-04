"""Retrieval coverage must be reported, not implied.

Measured 2026-11-03: two of this agent's ten collections were empty — `cart_constructs` and
`cart_trials`, which are the two a construct-design question leans on. Coverage for the
car_design stage ran 38%.

`total_collections_searched` counted the question asked, not the yield: a search of an empty
collection succeeds and returns zero rows, and nothing in the retrieval path can tell that apart
from a search that found nothing relevant. See `hcls_common.corpus_coverage` and
docs/build/EMPTY_COLLECTIONS.md.
"""
import pytest

from src.models import AgentQuery, SearchHit
from src.rag_engine import COLLECTION_CONFIG, CARTRAGEngine

# Measured 2026-11-03: cart_constructs and cart_trials held nothing; the rest did.
EMPTY = {"cart_constructs", "cart_trials"}
POPULATED = set(COLLECTION_CONFIG) - EMPTY


class FakeEmbedder:
    def embed_text(self, text):
        return [0.1] * 384

    def embed(self, text):
        return [0.1] * 384


class FakeLLM:
    def generate(self, prompt, system_prompt=None, **kwargs):
        return "A CD19 construct with a 4-1BB costimulatory domain is the usual starting point."


@pytest.fixture
def engine(monkeypatch):
    eng = CARTRAGEngine(collection_manager=object(), embedder=FakeEmbedder(),
                        llm_client=FakeLLM())

    def fake_search_all(query_embedding, collections, top_k, filter_exprs,
                        weight_overrides=None):
        return [
            SearchHit(id=f"{c}-1", collection=c, score=0.8, text=f"a row from {c}",
                      metadata={})
            for c in collections if c in eng._populated
        ]

    eng._populated = set(POPULATED)
    monkeypatch.setattr(eng, "_search_all_collections", fake_search_all)
    monkeypatch.setattr(eng, "_build_prompt", lambda q, e: "prompt")
    return eng


class TestThePartialCorpusIsReported:
    def test_the_verdict_rides_on_the_evidence(self, engine):
        ev = engine.retrieve(AgentQuery(question="which construct for CD19?"))
        assert ev.retrieval_coverage is not None
        assert ev.retrieval_coverage["degraded"] is True

    def test_contributing_sits_next_to_searched(self, engine):
        ev = engine.retrieve(AgentQuery(question="which construct for CD19?"))
        assert ev.total_collections_searched == len(COLLECTION_CONFIG)
        assert ev.collections_contributing == len(POPULATED)
        assert ev.collections_contributing < ev.total_collections_searched

    def test_the_two_empty_collections_are_named(self, engine):
        ev = engine.retrieve(AgentQuery(question="which construct for CD19?"))
        assert set(ev.retrieval_coverage["silent_collections"]) == EMPTY

    def test_the_answer_discloses_it(self, engine):
        """The reader who sees only the prose must still be told."""
        answer = engine.query("which construct for CD19?")
        assert "Partial corpus" in answer
        assert "about the corpus, not about the question" in answer
        assert "cart_constructs" in answer

    def test_the_answer_itself_is_still_there(self, engine):
        """A disclosure appends; it does not replace."""
        answer = engine.query("which construct for CD19?")
        assert answer.startswith("A CD19 construct")


class TestAFullCorpusStaysQuiet:
    def test_no_disclosure_when_everything_answers(self, engine):
        engine._populated = set(COLLECTION_CONFIG)
        ev = engine.retrieve(AgentQuery(question="which construct for CD19?"))
        assert ev.retrieval_coverage["coverage"] == pytest.approx(1.0)
        assert ev.retrieval_coverage["degraded"] is False
        assert "Partial corpus" not in engine.query("which construct for CD19?")


class TestTheEmptyCorpusCase:
    def test_nothing_anywhere_reports_zero_coverage(self, engine):
        engine._populated = set()
        ev = engine.retrieve(AgentQuery(question="which construct for CD19?"))
        assert ev.retrieval_coverage["coverage"] == 0.0
        assert ev.collections_contributing == 0
        assert "Partial corpus" in engine.query("which construct for CD19?")
