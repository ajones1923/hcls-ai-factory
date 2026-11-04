"""What the output-honesty gate can actually see.

The gate was documented as scanning "every JSON answer", and it selected one key:

    key = next((k for k in _ANSWER_KEYS
                if isinstance(payload.get(k), str) and len(payload[k]) > 200), None)

Four ways past it, all exercised below:
  * a string of 200 characters or fewer -- "This platform cures the patient's disease." is 41;
  * prose nested under another key, which is the common shape in this repo;
  * a second answer key in the same payload, once the first one matched;
  * a list payload, which was not looked at at all.

Two more things this file pins: the gate must report itself in `X-HCLS-Governed` (it never
did -- the contextvar it read is set in a different middleware context), and it must never
rewrite retrieved evidence, which is the one promise its docstring makes about content.
"""
from fastapi import FastAPI
from fastapi.testclient import TestClient
from starlette.responses import StreamingResponse

from hcls_common.api_gate import install_governance

CURE = "This platform cures the patient's disease."          # 41 chars, block severity
LONG_CLINICAL = "Dosing for this patient follows the label. " * 8


def client():
    app = FastAPI()
    install_governance(app, service="scopetest")

    @app.get("/short-overclaim")
    def short_overclaim():
        return {"answer": CURE}

    @app.get("/nested")
    def nested():
        return {"result": {"detail": {"interpretation": LONG_CLINICAL}}}

    @app.get("/two-answers")
    def two_answers():
        return {"answer": LONG_CLINICAL, "summary": CURE}

    @app.get("/list-payload")
    def list_payload():
        return [{"answer": LONG_CLINICAL}, {"answer": CURE}]

    @app.get("/new-key")
    def new_key():
        return {"recommendation": CURE}

    @app.get("/evidence-only")
    def evidence_only():
        # A retrieved passage that happens to contain a cure claim: it is a quotation, and the
        # clinician checks the answer against it. Rewriting it would destroy the check.
        return {"answer": "See the retrieved passage.",
                "citations": [{"summary": CURE, "url": "https://example.org/x"}]}

    @app.get("/no-answer-key")
    def no_answer_key():
        return {"count": 3, "collections": ["a", "b"]}

    @app.get("/stream")
    def stream():
        def gen():
            yield b'{"answer": "' + CURE.encode() + b'"}'
        return StreamingResponse(gen(), media_type="text/event-stream")

    return TestClient(app)


def _sev(body):
    return [f["severity"] for f in body.get("honesty", [])]


class TestLengthIsNotAShield:
    def test_a_short_block_claim_is_withheld(self):
        body = client().get("/short-overclaim").json()
        assert body["withheld"] is True
        assert "Withheld by the output-honesty gate" in body["answer"]
        assert CURE not in body["answer"]

    def test_it_is_reported_as_block_severity(self):
        assert "block" in _sev(client().get("/short-overclaim").json())


class TestDepth:
    def test_nested_prose_is_scanned(self):
        body = client().get("/nested").json()
        assert "decision support" in body["result"]["detail"]["interpretation"].lower()

    def test_the_finding_names_the_path(self):
        body = client().get("/nested").json()
        assert body["honesty"][0]["path"] == "result.detail.interpretation"


class TestEveryAnswerNotJustTheFirst:
    def test_both_answer_keys_are_handled(self):
        body = client().get("/two-answers").json()
        assert "decision support" in body["answer"].lower()      # disclaimed
        assert body["withheld"] is True                           # and the other withheld
        assert CURE not in body["summary"]

    def test_the_header_counts_both(self):
        assert "scanned=2" in client().get("/two-answers").headers["X-HCLS-Honesty"]


class TestListPayload:
    def test_a_list_body_is_scanned(self):
        body = client().get("/list-payload").json()
        assert "decision support" in body[0]["answer"].lower()
        assert CURE not in body[1]["answer"]

    def test_findings_cannot_be_attached_but_the_header_reports(self):
        """There is no object to hang `honesty` on, so the header carries the verdict."""
        r = client().get("/list-payload")
        assert isinstance(r.json(), list)
        assert "block=1" in r.headers["X-HCLS-Honesty"]


class TestKeyCoverage:
    def test_recommendation_is_an_answer_key(self):
        assert client().get("/new-key").json()["withheld"] is True


class TestRetrievedEvidenceIsNeverRewritten:
    def test_a_passage_is_returned_verbatim(self):
        body = client().get("/evidence-only").json()
        assert body["citations"][0]["summary"] == CURE
        assert "withheld" not in body

    def test_and_the_answer_beside_it_was_still_scanned(self):
        assert "output-honesty" in client().get("/evidence-only").headers["X-HCLS-Governed"]


class TestHeaders:
    def test_the_gate_reports_itself(self):
        """It never did: the contextvar was read in a context that had not seen the .set()."""
        r = client().get("/nested")
        assert "output-honesty" in r.headers["X-HCLS-Governed"]

    def test_no_header_when_there_was_nothing_to_scan(self):
        r = client().get("/no-answer-key")
        assert "X-HCLS-Governed" not in r.headers
        assert "X-HCLS-Honesty" not in r.headers

    def test_a_clean_scan_says_so(self):
        r = client().get("/evidence-only")
        assert r.headers["X-HCLS-Honesty"].startswith("scanned=1")
        assert "block" not in r.headers["X-HCLS-Honesty"]


class TestStreaming:
    def test_a_stream_is_marked_unscanned_rather_than_passed_off_as_gated(self):
        r = client().get("/stream")
        assert "not-scanned" in r.headers["X-HCLS-Honesty"]
        assert "X-HCLS-Governed" not in r.headers
        assert CURE in r.text      # unchanged, and honestly labelled as unchanged
