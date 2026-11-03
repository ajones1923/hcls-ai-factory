"""With nothing retrieved, the prompt must not ask for citations that cannot exist.

`_build_prompt` injected the literal "No evidence found." as the evidence block and then kept
its instruction to "provide a comprehensive answer grounded in the evidence above. Cite
sources using the clickable markdown links provided in each evidence item." There were no
items and no links, so the instruction asked the model to cite things it did not have.

The same defect, and the same fix, is in the cart and precision-biomarker agents.
"""
from types import SimpleNamespace

from src.rag_engine import PGxRAGEngine


def _prompt(engine, sections_present: bool):
    """Build a prompt with and without retrieval, going through the real code path."""
    hit = SimpleNamespace(
        collection="pgx_drug_guidelines", id="1", score=0.91,
        text="CYP2C19 poor metabolisers should avoid clopidogrel.", metadata={},
    )
    by_coll = {"pgx_drug_guidelines": [hit]} if sections_present else {}
    # hits_by_collection() is a method on CrossCollectionResult, not an attribute.
    evidence = SimpleNamespace(
        hits_by_collection=lambda: by_coll,
        knowledge_context="",
    )
    return engine._build_prompt("Is clopidogrel appropriate?", evidence)


class TestNoEvidence:
    def test_prompt_does_not_ask_for_citations(self):
        p = _prompt(PGxRAGEngine.__new__(PGxRAGEngine), sections_present=False)
        assert "Cite sources" not in p
        assert "grounded in the evidence above" not in p

    def test_prompt_states_there_is_no_evidence(self):
        p = _prompt(PGxRAGEngine.__new__(PGxRAGEngine), sections_present=False)
        assert "NO RETRIEVED EVIDENCE" in p
        assert "no source material was available" in p

    def test_prompt_forbids_fabricated_citations(self):
        p = _prompt(PGxRAGEngine.__new__(PGxRAGEngine), sections_present=False)
        assert "Do NOT cite" in p
        assert "do not fabricate" in p

    def test_the_old_placeholder_is_gone(self):
        """"No evidence found." used to be passed off as the evidence block."""
        p = _prompt(PGxRAGEngine.__new__(PGxRAGEngine), sections_present=False)
        assert "## Retrieved Evidence" not in p


class TestWithEvidence:
    def test_grounded_prompt_is_unchanged(self):
        p = _prompt(PGxRAGEngine.__new__(PGxRAGEngine), sections_present=True)
        assert "## Retrieved Evidence" in p
        assert "Cite sources" in p
        assert "NO RETRIEVED EVIDENCE" not in p

    def test_the_question_is_still_asked(self):
        p = _prompt(PGxRAGEngine.__new__(PGxRAGEngine), sections_present=True)
        assert "Is clopidogrel appropriate?" in p
