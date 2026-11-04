"""With nothing retrieved, the prompt must not ask for citations that cannot exist.

`_build_messages` built one user turn for every case. With no hits the evidence block is
empty, so the prompt read "**Retrieved Evidence (0 items ...):**" followed by nothing, and
then still instructed: "Cite sources using the [Label:ID] format above." There is no format
above and nothing to cite, so the model supplies both from its own weights — which is how an
ungrounded answer comes back looking sourced.

The same defect, and the same fix, is in the cart, pharmacogenomics, precision-biomarker and
single-cell agents.
"""
import pytest

from src.rag_engine import AutoimmuneRAGEngine, CrossCollectionResult, SearchHit


@pytest.fixture
def engine():
    """A real engine object without its constructor's Milvus/embedder dependencies."""
    eng = AutoimmuneRAGEngine.__new__(AutoimmuneRAGEngine)
    import threading
    eng._conversation_lock = threading.Lock()
    eng._conversation_history = []

    class _Settings:
        collection_config = {"ai_hla_associations": {"label": "HLA"}}

    eng.settings = _Settings()
    return eng


def _user_turn(engine, hits):
    result = CrossCollectionResult(query="q", hits=hits, total_collections_searched=6,
                                  search_time_ms=12.0)
    messages = engine._build_messages("Is this patient a candidate for rituximab?", result)
    assert messages[-1]["role"] == "user"
    return messages[-1]["content"]


def _one_hit():
    return [SearchHit(collection="ai_hla_associations", id="7", score=0.9,
                      text="HLA-DRB1*04:01 is associated with rheumatoid arthritis.",
                      metadata={"pmid": "12345678"}, relevance="high")]


class TestNoEvidence:
    def test_prompt_does_not_ask_for_citations(self, engine):
        turn = _user_turn(engine, [])
        assert "Cite sources using the [Label:ID] format above" not in turn

    def test_prompt_states_there_is_no_evidence(self, engine):
        turn = _user_turn(engine, [])
        assert "NO RETRIEVED EVIDENCE" in turn
        assert "no source material was available" in turn

    def test_prompt_forbids_fabricated_identifiers(self, engine):
        turn = _user_turn(engine, [])
        assert "Do NOT cite" in turn
        assert "do not fabricate" in turn

    def test_no_empty_evidence_header_is_passed_off_as_evidence(self, engine):
        """"Retrieved Evidence (0 items ...)" followed by nothing read like a retrieval."""
        turn = _user_turn(engine, [])
        assert "Retrieved Evidence (0 items" not in turn

    def test_no_therapy_recommendation_is_requested_on_no_evidence(self, engine):
        turn = _user_turn(engine, [])
        assert "what evidence would be needed" in turn

    def test_the_question_is_still_asked(self, engine):
        assert "rituximab" in _user_turn(engine, [])


class TestWithEvidence:
    def test_the_grounded_prompt_is_unchanged(self, engine):
        turn = _user_turn(engine, _one_hit())
        assert "Cite sources using the [Label:ID] format above" in turn
        assert "Retrieved Evidence (1 items" in turn
        assert "NO RETRIEVED EVIDENCE" not in turn

    def test_the_passage_is_in_the_prompt(self, engine):
        assert "HLA-DRB1*04:01" in _user_turn(engine, _one_hit())


class TestHistoryIsKeptInBothBranches:
    def test_history_precedes_the_ungrounded_turn_too(self, engine):
        """The ungrounded branch returns early; it must not drop conversation history."""
        engine._conversation_history = [{"question": "earlier?", "answer": "earlier answer"}]
        result = CrossCollectionResult(query="q", hits=[], total_collections_searched=6)
        messages = engine._build_messages("now?", result)
        assert len(messages) == 3
        assert messages[0]["content"] == "earlier?"
        assert messages[1]["role"] == "assistant"
        assert "NO RETRIEVED EVIDENCE" in messages[2]["content"]
