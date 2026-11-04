"""With nothing retrieved, the prompt must not ask for citations that cannot exist.

`_synthesize_response` built one prompt for every case. `_build_context` returns the string
"No evidence found in the knowledge base." for an empty result set, and the instructions then
asked for an analysis "grounded in the retrieved evidence above" and spelled out the citation
formats to use — PMID links, NCT numbers, Cell Ontology ids, [Collection:record-id]. Given a
format and nothing to fill it with, the model fills it from its own weights, which is how an
ungrounded answer comes back reading as sourced.

The same defect, and the same fix, is in the cart, pharmacogenomics, precision-biomarker and
precision-autoimmune agents.
"""
import pytest

from src.rag_engine import SCSearchResult, SingleCellRAGEngine


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
    eng = SingleCellRAGEngine.__new__(SingleCellRAGEngine)
    eng.llm = _CapturingLLM()
    eng._conversation_history = []
    return eng


def _prompt(engine, results):
    engine._synthesize_response("Which cell types express CD19?", results)
    return engine.llm.prompt


def _one_hit():
    hit = SCSearchResult(collection="sc_cell_types", record_id=1, score=0.88,
                         text="CD19 marks B lineage cells.")
    hit.relevance = "high"
    return [hit]


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
        assert "CL_0000084" not in p

    def test_the_old_placeholder_is_not_passed_off_as_evidence(self, engine):
        p = _prompt(engine, [])
        assert "## Retrieved Evidence" not in p
        assert "No evidence found in the knowledge base." not in p

    def test_no_scored_annotation_is_requested_on_no_evidence(self, engine):
        p = _prompt(engine, [])
        assert "what data would be needed" in p


class TestBareLLMFallbackIsLabelled:
    """`/v1/sc/query` falls back to the bare LLM when the engine is unavailable.

    It returned `confidence=0.3` with `evidence=[]` — a number describing nothing, since
    there was no retrieval behind it to be 30 per cent confident in — and the caller could
    not tell that response apart from a grounded one.
    """

    @pytest.fixture
    def client(self, monkeypatch):
        from fastapi.testclient import TestClient

        import api.routes.sc_clinical as route_mod
        from api.main import app

        monkeypatch.setattr(route_mod, "_get_engine", lambda req: None)

        async def _fallback(req, prompt, system_prompt=""):
            _fallback.system_prompt = system_prompt
            return "general background"

        monkeypatch.setattr(route_mod, "_llm_fallback", _fallback)
        self._fallback = _fallback
        return TestClient(app)

    def test_confidence_is_zero_and_grounded_is_false(self, client):
        body = client.post("/v1/sc/query", json={"question": "What marks B cells?"}).json()
        assert body["confidence"] == 0.0
        assert body["grounded"] is False
        assert body["evidence"] == []

    def test_the_system_prompt_does_not_ask_for_evidence_based_answers(self, client):
        client.post("/v1/sc/query", json={"question": "What marks B cells?"})
        sp = self._fallback.system_prompt
        assert "no retrieved evidence" in sp.lower()
        assert "Do not cite" in sp


class TestWithEvidence:
    def test_the_grounded_prompt_is_unchanged(self, engine):
        p = _prompt(engine, _one_hit())
        assert "## Retrieved Evidence" in p
        assert "grounded in the retrieved evidence above" in p
        assert "NO RETRIEVED EVIDENCE" not in p

    def test_the_question_is_still_asked(self, engine):
        p = _prompt(engine, _one_hit())
        assert "Which cell types express CD19?" in p
