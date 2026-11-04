"""An answer with no retrieved evidence must not be presented as a grounded one.

`/v1/trial/query` used to start `confidence` at 0.5 and add 0.05 per hit, so a question that
retrieved nothing still returned 0.5 — a coin flip reported as moderate confidence — while the
prompt asked the model to "provide a comprehensive answer citing specific evidence" over an
empty context. `guidelines_cited` was then filled by looking for guideline names in the
model's own prose, so a sentence mentioning ICH E6 produced a citation of ICH E6 that nothing
had retrieved.
"""
from contextlib import asynccontextmanager

import pytest


@pytest.fixture
def client():
    from fastapi.testclient import TestClient

    from api.main import _rate_limit_store, app

    @asynccontextmanager
    async def _noop_lifespan(_app):
        yield

    _rate_limit_store.clear()
    app.router.lifespan_context = _noop_lifespan
    app.state.engine = None
    app.state.manager = None
    app.state.workflow_engine = None
    app.state.llm_client = None
    app.state.metrics = {"requests_total": 0, "query_requests_total": 0}
    with TestClient(app) as c:
        yield c


class _EmptyEngine:
    """Retrieval that is reachable and simply finds nothing."""

    def search(self, *_args, **_kwargs):
        return []


class _HitEngine:
    def __init__(self, text):
        self._text = text

    def search(self, *_args, **_kwargs):
        return [{"text": self._text, "collection": "trial_guidelines", "score": 0.8, "metadata": {}}]


class _StubLLM:
    def __init__(self, reply):
        self.reply = reply
        self.prompts = []

    def generate(self, prompt, **_kwargs):
        self.prompts.append(prompt)
        return self.reply


def _ask(client, question="What endpoint suits a phase III oncology trial?"):
    return client.post("/v1/trial/query", json={"question": question})


class TestNoEvidence:
    def test_confidence_is_zero_not_half(self, client):
        client.app.state.engine = _EmptyEngine()
        client.app.state.llm_client = _StubLLM("Generally, overall survival is preferred.")
        data = _ask(client).json()
        assert data["confidence"] == 0.0

    def test_response_is_flagged_ungrounded(self, client):
        client.app.state.engine = _EmptyEngine()
        client.app.state.llm_client = _StubLLM("Generally, overall survival is preferred.")
        assert _ask(client).json()["grounded"] is False

    def test_answer_says_it_is_ungrounded(self, client):
        client.app.state.engine = _EmptyEngine()
        client.app.state.llm_client = _StubLLM("Generally, overall survival is preferred.")
        assert "no evidence was retrieved" in _ask(client).json()["answer"].lower()

    def test_prompt_does_not_demand_citations_it_cannot_have(self, client):
        client.app.state.engine = _EmptyEngine()
        llm = _StubLLM("Generally, overall survival is preferred.")
        client.app.state.llm_client = llm
        _ask(client)
        prompt = llm.prompts[0].lower()
        assert "citing specific evidence" not in prompt
        assert "no evidence was retrieved" in prompt

    def test_guideline_named_only_by_the_model_is_not_cited(self, client):
        """The model writing "ICH E6" is not a citation when nothing retrieved it."""
        client.app.state.engine = _EmptyEngine()
        client.app.state.llm_client = _StubLLM("Per ICH E6 and 21 CFR, use overall survival.")
        assert _ask(client).json()["guidelines_cited"] == []


class TestWithEvidence:
    def test_grounded_answer_keeps_its_confidence(self, client):
        client.app.state.engine = _HitEngine("ICH E6 requires a prespecified primary endpoint.")
        client.app.state.llm_client = _StubLLM("Per ICH E6, prespecify the primary endpoint.")
        data = _ask(client).json()
        assert data["grounded"] is True
        assert data["confidence"] > 0.0

    def test_guideline_present_in_evidence_is_cited(self, client):
        client.app.state.engine = _HitEngine("ICH E6 requires a prespecified primary endpoint.")
        client.app.state.llm_client = _StubLLM("Per ICH E6, prespecify the primary endpoint.")
        assert "ICH E6" in _ask(client).json()["guidelines_cited"]

    def test_guideline_in_prose_but_not_in_evidence_is_not_cited(self, client):
        client.app.state.engine = _HitEngine("ICH E6 requires a prespecified primary endpoint.")
        client.app.state.llm_client = _StubLLM("Per ICH E9 the estimand must be defined.")
        assert "ICH E9" not in _ask(client).json()["guidelines_cited"]
