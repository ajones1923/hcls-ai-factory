"""With nothing retrieved, the prompt must not ask for citations that cannot exist.

G1. `NeuroRAGEngine._synthesize_response` built one prompt for every case. `_build_context`
returns the string "No evidence found in the knowledge base." for an empty result set, and the
instructions then asked for an answer "grounded in the retrieved evidence above" and spelled
out the citation formats to use -- PMID links, NCT numbers, [Collection:record-id]. Given a
format and nothing to fill it with, the model fills it from its own weights, which is how an
ungrounded answer comes back reading as sourced. Against the measured corpus (4,764 vectors
across 114 collections, 40 of them empty) a neurology query returns 0 rows routinely, and the
reproduction handed the model the tPA question with an empty evidence block, an instruction to
ground in it, and two worked examples of the citation format.

The same defect, and the same fix, is in the cart, pharmacogenomics, precision-biomarker,
precision-autoimmune and single-cell agents. Neurology was the one that was missed: its
grounding work went into api/routes/neuro_clinical.py and src/agent.py and left the RAG engine
on the old prompt.
"""
import pytest

from src.rag_engine import NeuroRAGEngine, NeuroSearchResult


class _CapturingLLM:
    """Stands in for the LLM client and keeps the prompt it was handed."""

    def __init__(self):
        self.prompt = None

    def generate(self, prompt, system_prompt=None, max_tokens=None, temperature=None):
        self.prompt = prompt
        return "an answer"


@pytest.fixture
def engine():
    """A real engine object without its constructor's Milvus/embedder dependencies."""
    eng = NeuroRAGEngine.__new__(NeuroRAGEngine)
    eng.llm = _CapturingLLM()
    eng._conversation_history = []
    return eng


_QUESTION = "68yo acute stroke, NIHSS 14, onset 2 hours -- is tPA indicated?"


def _prompt(engine, results):
    engine._synthesize_response(_QUESTION, results)
    return engine.llm.prompt


def _one_hit():
    return [
        NeuroSearchResult(
            collection="neuro_cerebrovascular",
            record_id="rec-1",
            score=0.88,
            text="Alteplase is indicated within 4.5 hours of established onset.",
            relevance="high",
        )
    ]


class TestNoEvidence:
    def test_prompt_does_not_ask_for_grounding_that_does_not_exist(self, engine):
        p = _prompt(engine, [])
        assert "grounded in the retrieved evidence above" not in p

    def test_prompt_states_there_is_no_evidence(self, engine):
        p = _prompt(engine, [])
        assert "NO RETRIEVED EVIDENCE" in p
        assert "no source material was available" in p

    def test_prompt_forbids_fabricated_citations(self, engine):
        p = _prompt(engine, [])
        assert "Do NOT cite" in p
        assert "do not fabricate" in p

    def test_no_citation_format_is_offered_to_fill(self, engine):
        """The example links were the mechanism: a template invites a plausible id."""
        p = _prompt(engine, [])
        assert "pubmed.ncbi.nlm.nih.gov/12345678" not in p
        assert "NCT01234567" not in p
        assert "[Collection:record-id]" not in p

    def test_the_old_placeholder_is_not_passed_off_as_evidence(self, engine):
        p = _prompt(engine, [])
        assert "## Retrieved Evidence" not in p
        assert "No evidence found in the knowledge base." not in p

    def test_no_treatment_verdict_is_requested_on_no_evidence(self, engine):
        """The reproduction was a tPA question. Nothing retrieved means no verdict."""
        p = _prompt(engine, [])
        assert "thrombolysis" in p
        assert "what data would be needed" in p

    def test_no_scale_score_is_requested_on_no_evidence(self, engine):
        p = _prompt(engine, [])
        assert "do not report a clinical scale score" in p.lower()


class TestWithEvidence:
    def test_the_grounded_prompt_is_unchanged(self, engine):
        p = _prompt(engine, _one_hit())
        assert "## Retrieved Evidence" in p
        assert "grounded in the retrieved evidence above" in p
        assert "NO RETRIEVED EVIDENCE" not in p

    def test_the_citation_format_is_still_offered_when_there_is_something_to_cite(self, engine):
        p = _prompt(engine, _one_hit())
        assert "pubmed.ncbi.nlm.nih.gov/12345678" in p

    def test_the_question_is_still_asked(self, engine):
        assert _QUESTION in _prompt(engine, _one_hit())
