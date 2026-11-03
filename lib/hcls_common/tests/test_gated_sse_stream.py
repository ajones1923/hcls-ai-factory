"""A streamed clinical answer was the last ungated HTTP path in the factory.

`install_output_honesty` skips `text/event-stream` -- it cannot buffer a stream without
destroying the reason for streaming -- and `POST /query/stream` on the precision-biomarker and
precision-autoimmune agents streams generated clinical prose. So those two routes served
clinical answers with no overclaim register and no decision-support disclaimer, while the
non-streaming `/query` on the same services was gated.

`gate_token_stream` judges the COMPLETED answer and emits a terminal verdict. It is honestly
weaker than the non-streaming path: the tokens are already on the wire, so the verdict is a
correction after the fact and the client has to honour it. These tests pin that the verdict is
always emitted, that it agrees with the non-streaming gate, and that the tokens still flow
unbuffered while it happens.
"""
import pytest

from hcls_common.api_gate import gate_text, gate_token_stream

CURE_TOKENS = ["This platform ", "cures the patient's ", "disease."]
CLINICAL_TOKENS = ["Dosing for this patient follows the label. "] * 8


def run(tokens):
    tok, verdicts = [], []
    for kind, payload in gate_token_stream(iter(tokens)):
        (tok if kind == "token" else verdicts).append(payload)
    assert len(verdicts) == 1, "exactly one terminal verdict"
    return tok, verdicts[0]


class TestTheVerdictIsAlwaysEmitted:
    def test_a_clean_stream_still_ends_with_a_verdict(self):
        _, verdict = run(["Retrieved ", "3 references."])
        assert verdict["withheld"] is False

    def test_an_empty_stream_ends_with_a_verdict(self):
        tokens, verdict = run([])
        assert tokens == []
        assert verdict["withheld"] is False

    def test_the_verdict_comes_last(self):
        seq = [kind for kind, _ in gate_token_stream(iter(CURE_TOKENS))]
        assert seq[-1] == "verdict"
        assert seq[:-1] == ["token"] * len(CURE_TOKENS)


class TestBlockSeverityIsReported:
    def test_a_cure_claim_is_marked_withheld(self):
        _, verdict = run(CURE_TOKENS)
        assert verdict["withheld"] is True

    def test_the_replacement_text_does_not_contain_the_claim(self):
        _, verdict = run(CURE_TOKENS)
        assert "cures the patient" not in verdict["text"]
        assert "Withheld by the output-honesty gate" in verdict["text"]

    def test_the_findings_are_carried_so_a_client_can_show_why(self):
        _, verdict = run(CURE_TOKENS)
        assert any(f["severity"] == "block" for f in verdict["findings"])

    def test_the_tokens_were_still_delivered(self):
        """Honest about the limitation: a stream cannot be un-sent, only corrected."""
        tokens, _ = run(CURE_TOKENS)
        assert "".join(tokens) == "".join(CURE_TOKENS)

    def test_enforcement_off_does_not_mark_it_withheld(self, monkeypatch):
        monkeypatch.setenv("HCLS_HONESTY_ENFORCE", "0")
        _, verdict = run(CURE_TOKENS)
        assert verdict["withheld"] is False
        assert "cures the patient" in verdict["text"]


class TestTheDeltaLetsANaiveClientSeeTheVerdict:
    """Both routes' existing contract is "concatenate the tokens". A client written against
    that contract ignores an unfamiliar verdict event, so the verdict also rides the token
    channel as a delta -- otherwise the fix only protects clients updated for it."""

    def test_a_clean_answer_has_an_empty_delta(self):
        _, verdict = run(["Retrieved ", "3 references."])
        assert verdict["delta"] == ""

    def test_an_appended_disclaimer_is_the_delta(self):
        _, verdict = run(CLINICAL_TOKENS)
        assert verdict["delta"]
        assert "".join(CLINICAL_TOKENS) + verdict["delta"] == verdict["text"]

    def test_a_withheld_answer_appends_the_notice(self):
        """The tokens cannot be recalled, so the notice is appended rather than substituted."""
        _, verdict = run(CURE_TOKENS)
        assert "Withheld by the output-honesty gate" in verdict["delta"]

    def test_concatenating_tokens_plus_delta_never_loses_the_verdict(self):
        for tokens in (CURE_TOKENS, CLINICAL_TOKENS, ["ok"]):
            toks, verdict = run(tokens)
            naive = "".join(toks) + verdict["delta"]
            if verdict["findings"]:
                assert "decision support" in naive.lower() or verdict["withheld"]


class TestItAgreesWithTheNonStreamingGate:
    @pytest.mark.parametrize("tokens", [CURE_TOKENS, CLINICAL_TOKENS, ["ok"], []])
    def test_the_verdict_text_matches_gate_text(self, tokens):
        _, verdict = run(tokens)
        assert verdict["text"] == gate_text("".join(tokens))[0]


class TestClinicalProseGainsTheDisclaimer:
    def test_the_verdict_text_carries_the_disclaimer(self):
        _, verdict = run(CLINICAL_TOKENS)
        assert "decision support" in verdict["text"].lower()
        assert verdict["withheld"] is False


class TestStreamingIsNotBuffered:
    def test_tokens_are_yielded_before_the_source_is_exhausted(self):
        """Buffering the whole answer first would defeat the point of the route."""
        produced, consumed = [], []

        def source():
            for t in ("a", "b", "c"):
                produced.append(t)
                yield t

        for kind, payload in gate_token_stream(source()):
            if kind == "token":
                consumed.append(payload)
                # at this point the source must not have run ahead to the end
                assert len(produced) == len(consumed)

    def test_non_string_tokens_do_not_break_the_stream(self):
        tokens, verdict = run([1, 2])
        assert tokens == ["1", "2"]
        assert verdict["text"] == "12"
