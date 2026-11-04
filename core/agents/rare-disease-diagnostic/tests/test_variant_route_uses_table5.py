"""`/variants/interpret` must classify with Table 5, not with the model's prose.

The deterministic ACMG/AMP classifier in `src/decision_support.py` was corrected to the
Richards 2015 Table 5 combining rules -- including the footnote that contradictory
pathogenic and benign evidence yields VUS -- and then **no HTTP caller reached it**. This
route asked the model for a classification and reported the label it stated, so the
corrected rules governed nothing a client could call. Fixing the classifier without
wiring it is a fix nobody can use, which is why these tests go through the HTTP layer
rather than calling `combine_acmg_criteria` directly (it has its own unit tests in
test_decision_support.py).

The division of labour under test: the model proposes CRITERIA from retrieved evidence,
Table 5 decides the CLASS.
"""
import os

import pytest

os.environ.setdefault("ANTHROPIC_API_KEY", "test-key-placeholder")

from fastapi.testclient import TestClient

from api.main import app


class _FakeLLM:
    """Returns one canned interpretation, and records the prompt it was given."""

    def __init__(self, answer: str):
        self.answer = answer
        self.prompts: list[str] = []

    def generate(self, prompt: str, max_tokens: int = 0) -> str:
        self.prompts.append(prompt)
        return self.answer


class _FakeEngine:
    """Non-empty retrieval, so the route does not take its no-evidence refusal path."""

    def search(self, query, top_k=10):
        return [{"content": "ClinVar: two submitters report this variant as pathogenic."}]


@pytest.fixture
def client():
    """A client that does not run the app lifespan.

    Startup builds a RAG engine and an LLM client from settings; this test supplies both
    directly, and running startup would also make the test depend on services that are
    not part of what it is checking.
    """
    return TestClient(app)


def _interpret(client, answer, **overrides):
    llm = _FakeLLM(answer)
    app.state.llm_client = llm
    app.state.engine = _FakeEngine()
    payload = {"gene": "BRCA1", "variant": "c.68_69delAG"}
    payload.update(overrides)
    response = client.post("/v1/diagnostic/variants/interpret", json=payload)
    assert response.status_code == 200, response.text
    return response.json()


class TestTable5GovernsTheClassification:
    def test_cited_criteria_decide_the_class(self, client):
        """PVS1 + PM2 is Table 5 'Likely pathogenic (i)', whatever the prose says."""
        data = _interpret(
            client,
            "PVS1 is met: the deletion causes a frameshift. PM2 is met: absent from gnomAD.\n"
            "Final classification: Likely Pathogenic",
        )
        assert data["classification"] == "Likely Pathogenic"
        assert data["classification_rule"], "the firing rule must be reported"
        assert "Table 5" in data["classification_source"]

    def test_table5_overrides_a_disagreeing_label(self, client):
        """The model's stated label loses to the rules applied to its own criteria.

        Here it cites strong pathogenic evidence and then states Benign. Reporting Benign
        would be the most dangerous available answer, and it is exactly what the route
        used to do.
        """
        data = _interpret(
            client,
            "PVS1 is met: nonsense variant in a gene where loss of function causes disease. "
            "PS2 is met: confirmed de novo. PM2 is met: absent from population databases.\n"
            "Final classification: Benign",
        )
        assert data["classification"] == "Pathogenic", data["notes"]
        joined = " ".join(data["notes"])
        assert "stated 'Benign'" in joined, joined
        assert "Table 5" in joined

    def test_contradictory_evidence_is_vus_not_a_side(self, client):
        """Richards Table 5 footnote: contradictory benign and pathogenic evidence -> VUS.

        Both sides must actually fire for the footnote to apply. PVS1 + PS2 + PM2 satisfies
        Pathogenic (i)(a), and BA1 is stand-alone Benign, so the two rule sets collide.
        """
        data = _interpret(
            client,
            "PVS1 is met: frameshift in a loss-of-function gene. PS2 is met: confirmed de "
            "novo. PM2 is met: absent from controls. "
            "BA1 is met: allele frequency is 7% in gnomAD African/African American.\n"
            "Final classification: Pathogenic",
        )
        assert data["classification"] == "VUS", data["notes"]
        assert "contradictor" in (data["classification_rule"] or "").lower()

    def test_pvs1_alone_against_ba1_is_benign_not_contradictory(self, client):
        """A surprising result worth pinning, because it looks like a bug and is not.

        PVS1 on its own satisfies no Table 5 pathogenic rule set -- Likely pathogenic (i) is
        PVS1 AND 1 PM, not PVS1 alone -- while BA1 is stand-alone Benign. So there is no
        collision to detect and Benign is the correct combination. I asserted VUS here
        first; Table 5 was right and the expectation was wrong.

        Note this is the combining rule in isolation. `ACMGVariantClassifier._ba1_exception`
        holds the Ghosh 2018 carve-out for genes where a high population frequency does not
        clear a variant, and CLINICAL_REVIEW records that scope.
        """
        data = _interpret(
            client,
            "PVS1 is met: frameshift in a loss-of-function gene. "
            "BA1 is met: allele frequency is 7% in gnomAD African/African American.\n"
            "Final classification: Pathogenic",
        )
        assert data["classification"] == "Benign", data["notes"]
        assert "contradictor" not in (data["classification_rule"] or "").lower()

    def test_criteria_marked_not_met_do_not_count(self, client):
        """Only criteria the interpretation says are MET may enter the combination."""
        data = _interpret(
            client,
            "PVS1 is not met: this is a missense variant. PM2 is met: absent from gnomAD.\n"
            "Final classification: Pathogenic",
        )
        # PM2 alone satisfies no Table 5 rule set.
        assert data["classification"] == "VUS", data["notes"]
        assert "PVS1" not in (data["classification_rule"] or "")


class TestAClassWithoutCriteriaIsNotAClassification:
    def test_label_with_no_cited_criteria_is_demoted(self, client):
        """ACMG defines the class AS the combination of criteria.

        A stated label with nothing to combine cannot be audited, so it is not reported as
        a classification. The converse rule already held: criteria are never derived from
        the classification.
        """
        data = _interpret(
            client,
            "This variant is concerning given the family history.\n"
            "Final classification: Pathogenic",
        )
        assert data["classification"] == "Unclassified", data["notes"]
        assert data["classification_rule"] is None
        joined = " ".join(data["notes"])
        assert "cited no ACMG criteria" in joined, joined
        # The model's claim is preserved rather than discarded.
        assert "'Pathogenic'" in joined

    def test_no_retrieved_evidence_still_refuses(self, client):
        """The existing refusal path must survive this change."""
        app.state.llm_client = _FakeLLM("Final classification: Pathogenic")

        class _Empty:
            def search(self, query, top_k=10):
                return []

        app.state.engine = _Empty()
        response = client.post(
            "/v1/diagnostic/variants/interpret",
            json={"gene": "BRCA1", "variant": "c.68_69delAG"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["classification"] == "Unclassified"
        assert any("No evidence was retrieved" in n for n in data["notes"])
