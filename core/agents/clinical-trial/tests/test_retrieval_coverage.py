"""The response must say when the corpus, not the question, was the limit.

Measured 2026-11-03: ten of this agent's fourteen collections were empty, including
`trial_eligibility`, which PATIENT_MATCHING weights at 0.25 and SITE_SELECTION's `trial_sites`
at 0.25. 73% of the declared retrieval weight for patient matching pointed at nothing.

Nothing in the response said so. `collections_searched` read 14 either way, because a search of
an empty collection succeeds and returns zero rows -- indistinguishable, at every layer, from a
search that found nothing relevant. These tests hold the line on reporting the difference.
"""
import pytest

from src.agent import TrialWorkflowType
from src.rag_engine import TrialRAGEngine

# The four that held anything on the day this was written (156/150/22/60 rows).
POPULATED = {"trial_protocols", "trial_literature", "trial_regulatory", "genomic_evidence"}


class FakeEmbedder:
    def embed_text(self, text):
        return [0.1] * 384


@pytest.fixture
def engine(monkeypatch, tmp_path):
    """An engine whose Milvus returns rows for POPULATED and nothing for the rest."""
    eng = TrialRAGEngine(milvus_client=object(), embedding_model=FakeEmbedder(),
                         llm_client=None, session_id="coverage-test")

    def fake_search(collection_name, query_vector, top_k, filter_expr=None):
        if collection_name not in eng._populated:
            return []
        return [{"id": f"{collection_name}-1", "score": 0.8,
                 "text": f"a row from {collection_name}"}]

    eng._populated = set(POPULATED)
    monkeypatch.setattr(eng, "_search_collection", fake_search)
    return eng


class TestThePartialCorpusIsReported:
    def test_the_verdict_is_attached_to_the_response(self, engine):
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.retrieval_coverage is not None
        assert r.retrieval_coverage["degraded"] is True

    def test_coverage_is_the_measured_29_percent(self, engine):
        """28.8% through the engine's own weights.

        Note which weights those are. This agent carries two weight tables:
        `WORKFLOW_COLLECTION_WEIGHTS` in vector_collections.py, which reads as the declared
        budget and gives 27%, and `WORKFLOW_COLLECTION_BOOST` in agent.py, which is what
        `_get_boosted_weights` actually multiplies against the base weights. Only the second one
        reaches retrieval -- `get_search_weights()`, the accessor for the first, has no caller in
        this agent. The two agree to within two points here, so the finding does not turn on
        which is authoritative, but this test follows the code path that runs.
        """
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.retrieval_coverage["coverage"] == pytest.approx(0.288, abs=0.01)

    def test_contributing_is_reported_alongside_searched(self, engine):
        """The number that moves with the corpus, next to the one that does not."""
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.collections_searched == 14
        assert r.collections_contributing == 4

    def test_the_heaviest_empty_source_is_named(self, engine):
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.retrieval_coverage["silent_collections"][0] == "trial_eligibility"

    def test_the_answer_itself_discloses_it(self, engine):
        """A reader who sees only the prose must still be told."""
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert "Partial corpus" in r.answer
        assert "trial_eligibility" in r.answer
        assert "10 of 14" in r.answer

    def test_site_selection_is_hit_too(self, engine):
        """trial_sites at 0.25 and trial_investigators at 0.18 were both empty."""
        r = engine.query("where should we run this?",
                         workflow=TrialWorkflowType.SITE_SELECTION)
        assert r.retrieval_coverage["degraded"] is True
        assert "trial_sites" in r.retrieval_coverage["silent_collections"]
        assert "trial_investigators" in r.retrieval_coverage["silent_collections"]


class TestAFullCorpusStaysQuiet:
    def test_no_disclosure_when_everything_answers(self, engine):
        engine._populated = set(engine._get_boosted_weights(
            TrialWorkflowType.PATIENT_MATCHING).keys())
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.retrieval_coverage["coverage"] == pytest.approx(1.0)
        assert r.retrieval_coverage["degraded"] is False
        assert "Partial corpus" not in r.answer
        assert r.retrieval_coverage["silent_collections"] == []

    def test_one_trivial_gap_is_recorded_but_not_announced(self, engine):
        """trial_adaptive is 0.01 of the patient-matching budget."""
        weights = engine._get_boosted_weights(TrialWorkflowType.PATIENT_MATCHING)
        engine._populated = {c for c in weights if c != "trial_adaptive"}
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.retrieval_coverage["silent_collections"] == ["trial_adaptive"]
        assert r.retrieval_coverage["degraded"] is False
        assert "Partial corpus" not in r.answer


class TestTheOperatorScriptComputesTheSameWeights:
    """`scripts/check_corpus_weights.py` reports coverage per workflow across every subject, and
    to do that it reimplements `_get_boosted_weights` in a subprocess -- the agent modules cannot
    all be imported into one interpreter, since several shadow each other's `src` package.

    A reimplementation drifts. This pins the two together on the subject where the arithmetic
    was first checked: read the boost table as if it were a weight map (the obvious mistake) and
    patient matching reads 36% instead of 29%, because every collection the boost table does not
    name silently drops out of the denominator.
    """

    @staticmethod
    def _as_the_script_does(boosts, base):
        boosted = {n: w * float(boosts.get(n, 1.0)) for n, w in base.items()}
        total = sum(boosted.values())
        return {n: w / total for n, w in boosted.items()} if total > 0 else boosted

    def test_the_script_arithmetic_matches_the_engine(self, engine):
        from src.agent import WORKFLOW_COLLECTION_BOOST
        from src.rag_engine import COLLECTION_CONFIG

        base = {n: cfg.get("weight", 0.05) for n, cfg in COLLECTION_CONFIG.items()}
        for workflow, boosts in WORKFLOW_COLLECTION_BOOST.items():
            mine = self._as_the_script_does(boosts, base)
            theirs = engine._get_boosted_weights(workflow)
            assert mine.keys() == theirs.keys(), workflow
            for name in mine:
                assert mine[name] == pytest.approx(theirs[name]), f"{workflow} / {name}"


class TestTheEmptyCorpusCase:
    def test_nothing_anywhere_reports_zero_coverage(self, engine):
        engine._populated = set()
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert r.retrieval_coverage["coverage"] == 0.0
        assert r.collections_contributing == 0
        assert "Partial corpus" in r.answer

    def test_the_disclosure_does_not_blame_the_question(self, engine):
        """An empty corpus must not read as "no trials match this patient"."""
        engine._populated = set()
        r = engine.query("which trials fit this patient?",
                         workflow=TrialWorkflowType.PATIENT_MATCHING)
        assert "about the corpus, not about the question" in r.answer
