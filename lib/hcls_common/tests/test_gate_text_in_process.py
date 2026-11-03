"""The output gate is middleware, and five UIs never cross the HTTP boundary.

`install_output_honesty` is documented as automatic since 2026-09-15. It is a FastAPI response
middleware, so it only sees answers that go out over HTTP. The cart, pharmacogenomics,
precision-biomarker and precision-autoimmune UIs import their agent's RAG engine and call it
in-process, and the precision-intelligence chat UI calls `engine.query_stream` the same way --
so the gate never ran on the surface a clinician actually reads.

`gate_text` makes the same decisions outside HTTP, and these tests pin that the two agree.
"""
import pytest

from hcls_common.api_gate import DISCLAIMER, gate_text

CURE = "This platform cures the patient's disease."
CLINICAL = "Dosing for this patient follows the label. " * 8


class TestBlockSeverityIsWithheld:
    def test_a_cure_claim_does_not_reach_the_screen(self):
        shown, findings = gate_text(CURE)
        assert CURE not in shown
        assert "Withheld by the output-honesty gate" in shown
        assert any(f["severity"] == "block" for f in findings)

    def test_the_reason_is_shown_not_just_the_refusal(self):
        shown, _ = gate_text(CURE)
        assert "- " in shown          # the finding messages are listed

    def test_the_withheld_text_still_carries_the_disclaimer(self):
        shown, _ = gate_text(CURE)
        assert "decision support" in shown.lower()

    def test_enforcement_can_be_turned_off(self, monkeypatch):
        monkeypatch.setenv("HCLS_HONESTY_ENFORCE", "0")
        shown, findings = gate_text(CURE)
        assert CURE in shown                       # annotated, not withheld
        assert any(f["severity"] == "block" for f in findings)


class TestClinicalProseGainsTheDisclaimer:
    def test_a_long_clinical_answer_is_disclaimed(self):
        shown, findings = gate_text(CLINICAL)
        assert shown.endswith(DISCLAIMER)
        assert findings

    def test_an_answer_that_already_says_it_is_not_disclaimed_twice(self):
        once, _ = gate_text(CLINICAL)
        twice, _ = gate_text(once)
        assert twice == once


class TestItLeavesAloneWhatItShould:
    def test_a_short_status_string_is_untouched(self):
        assert gate_text("ok") == ("ok", [])

    @pytest.mark.parametrize("value", ["", None, 0, [], {"a": 1}])
    def test_non_text_passes_through(self, value):
        shown, findings = gate_text(value)
        assert shown == value and findings == []

    def test_a_retrieved_passage_listing_is_not_rewritten(self):
        """Evidence-only output is quoted source material, same rule as the middleware."""
        passages = "1. [pgx_guidelines] CYP2C19 poor metabolisers should avoid clopidogrel.\n"
        shown, _ = gate_text(passages)
        assert passages in shown


class TestItAgreesWithTheMiddleware:
    def test_the_same_text_gets_the_same_verdict_over_http(self):
        from fastapi import FastAPI
        from fastapi.testclient import TestClient

        from hcls_common.api_gate import install_governance

        app = FastAPI()
        install_governance(app, service="agreement")

        @app.get("/a")
        def a():
            return {"answer": CURE}

        @app.get("/b")
        def b():
            return {"answer": CLINICAL}

        client = TestClient(app)
        over_http = client.get("/a").json()["answer"]
        in_process = gate_text(CURE)[0]
        # Same verdict and the same reason, from one shared builder. The HTTP version adds one
        # sentence that is only true there -- that the retrieved evidence is still returned
        # alongside -- because the payload carries it and a bare string does not.
        for text in (over_http, in_process):
            assert text.startswith("**Withheld by the output-honesty gate.**")
            assert CURE not in text
            assert "HCLS_HONESTY_ENFORCE=0" in text
        assert "retrieved evidence is unchanged" in over_http
        assert "retrieved evidence is unchanged" not in in_process

        # Disclaiming, where there is no such difference, must match exactly.
        assert client.get("/b").json()["answer"] == gate_text(CLINICAL)[0]
