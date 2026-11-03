"""`/query/stream` served clinical prose with no gate at all.

The output-honesty gate is response middleware and skips `text/event-stream`, so this route --
which streams generated clinical answers token by token -- had no overclaim register and no
decision-support disclaimer, while the non-streaming `/query` beside it was gated.

These tests go through the real route and pin the wire format, because the fix is only useful
if a client can act on it: the terminal event has to arrive, say whether the answer was
withheld, and carry the replacement text.
"""
import json

import pytest
from fastapi.testclient import TestClient

CURE = "This platform cures the patient's disease."
CLINICAL = "Dosing for this patient follows the label. " * 8


def client_streaming(tokens):
    """The real app and the real route, with only the engine's token source stubbed."""
    from api.main import app

    class FakeEngine:
        llm = object()
        embedder = object()

        def query_stream(self, **kwargs):
            for t in tokens:
                yield {"type": "token", "content": t}
            yield {"type": "done"}

    app.state.engine = FakeEngine()
    return TestClient(app)


def events(response):
    out = []
    for line in response.text.splitlines():
        if line.startswith("data: "):
            body = line[6:]
            if body == "[DONE]":
                continue
            try:
                out.append(json.loads(body))
            except json.JSONDecodeError:
                pass
    return out


def post(client, question="What does this panel show?"):
    return client.post("/v1/query/stream", json={"question": question})


@pytest.fixture
def overclaiming():
    return client_streaming([CURE])


@pytest.fixture
def clinical():
    return client_streaming([CLINICAL])


class TestTheTerminalVerdictReachesTheClient:
    def test_an_honesty_event_is_emitted(self, clinical):
        kinds = [e.get("type") for e in events(post(clinical))]
        assert "honesty" in kinds

    def test_it_arrives_before_done(self, clinical):
        kinds = [e.get("type") for e in events(post(clinical)) if e.get("type") in ("honesty", "done")]
        assert kinds == ["honesty", "done"]

    def test_a_naive_concatenating_client_still_ends_up_with_the_gated_text(self, clinical):
        """The gate's delta goes out through the token channel, so a client that ignores the
        typed honesty event is not left holding the ungated prose."""
        toks = [e["content"] for e in events(post(clinical)) if e.get("type") == "token"]
        joined = "".join(toks)
        assert joined.startswith(CLINICAL)          # the answer itself still streamed
        assert "decision support" in joined.lower()  # and the verdict arrived with it


class TestAnOverclaimIsReported:
    def test_the_verdict_says_withheld(self, overclaiming):
        honesty = next(e for e in events(post(overclaiming)) if e.get("type") == "honesty")
        assert honesty["withheld"] is True

    def test_the_replacement_text_is_supplied(self, overclaiming):
        honesty = next(e for e in events(post(overclaiming)) if e.get("type") == "honesty")
        assert "Withheld by the output-honesty gate" in honesty["text"]
        assert "cures the patient" not in honesty["text"]

    def test_the_reason_is_supplied(self, overclaiming):
        honesty = next(e for e in events(post(overclaiming)) if e.get("type") == "honesty")
        assert any(f["severity"] == "block" for f in honesty["findings"])


class TestClinicalProseIsDisclaimed:
    def test_the_verdict_text_carries_the_disclaimer(self, clinical):
        honesty = next(e for e in events(post(clinical)) if e.get("type") == "honesty")
        assert "decision support" in honesty["text"].lower()
        assert honesty["withheld"] is False


class TestTheResponseSaysItWasGated:
    def test_the_header_claims_in_stream_gating_not_not_scanned(self, clinical):
        r = post(clinical)
        assert "in-stream" in r.headers["X-HCLS-Honesty"]
        assert "not-scanned" not in r.headers["X-HCLS-Honesty"]
