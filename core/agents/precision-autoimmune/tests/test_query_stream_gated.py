"""`/query/stream` served clinical prose with no gate at all.

The output-honesty gate is response middleware and skips `text/event-stream`, so this route --
which streams generated clinical answers token by token -- had no overclaim register and no
decision-support disclaimer, while the non-streaming `/query` beside it was gated.

Same defect and same fix as the precision-biomarker agent. The wire format differs (this route
emits `{"text": ...}` and a literal `[DONE]`), so the terminal verdict is pinned here too.
"""
import json

import pytest
from fastapi.testclient import TestClient

CURE = "This platform cures the patient's disease."
CLINICAL = "Therapy for this patient follows the label. " * 8


def client_streaming(tokens):
    import api.main as main_mod

    class FakeRag:
        def query_stream(self, **kwargs):
            yield from tokens

    main_mod._state["rag"] = FakeRag()
    return TestClient(main_mod.app)


def payloads(response):
    out = []
    for line in response.text.splitlines():
        if line.startswith("data: ") and line[6:] != "[DONE]":
            try:
                out.append(json.loads(line[6:]))
            except json.JSONDecodeError:
                pass
    return out


def post(client):
    return client.post("/query/stream", json={"question": "Is rituximab appropriate?"})


@pytest.fixture
def overclaiming():
    return client_streaming([CURE])


@pytest.fixture
def clinical():
    return client_streaming([CLINICAL])


class TestTheTerminalVerdictReachesTheClient:
    def test_an_honesty_payload_is_emitted(self, clinical):
        assert any("honesty" in p for p in payloads(post(clinical)))

    def test_a_naive_concatenating_client_still_ends_up_with_the_gated_text(self, clinical):
        """This route's contract is "concatenate every `text`". The verdict's delta rides that
        same channel, so a client that ignores `honesty` is not left holding ungated prose."""
        toks = [p["text"] for p in payloads(post(clinical)) if "text" in p]
        joined = "".join(toks)
        assert joined.startswith(CLINICAL)
        assert "decision support" in joined.lower()

    def test_done_still_terminates_the_stream(self, clinical):
        assert post(clinical).text.rstrip().endswith("data: [DONE]")


class TestAnOverclaimIsReported:
    def test_the_verdict_says_withheld(self, overclaiming):
        honesty = next(p["honesty"] for p in payloads(post(overclaiming)) if "honesty" in p)
        assert honesty["withheld"] is True
        assert "cures the patient" not in honesty["text"]

    def test_the_reason_is_supplied(self, overclaiming):
        honesty = next(p["honesty"] for p in payloads(post(overclaiming)) if "honesty" in p)
        assert any(f["severity"] == "block" for f in honesty["findings"])


class TestClinicalProseIsDisclaimed:
    def test_the_verdict_text_carries_the_disclaimer(self, clinical):
        honesty = next(p["honesty"] for p in payloads(post(clinical)) if "honesty" in p)
        assert "decision support" in honesty["text"].lower()


class TestTheResponseSaysItWasGated:
    def test_the_header_claims_in_stream_gating(self, clinical):
        hdr = post(clinical).headers["X-HCLS-Honesty"]
        assert "in-stream" in hdr and "not-scanned" not in hdr
