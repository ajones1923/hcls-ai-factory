"""The verification gate must fail closed, and must read the whole text.

Two bypasses, both in the layer whose entire job is to not be fooled.

**It stamped "verified" during an outage.** `verify_claims` swallowed every exception from the
model and returned empty verdicts. Empty verdicts mean zero refuted and zero unsupported, and
`verify_text` reads those zeros as "nothing wrong with this text" -- so any LLM outage produced
`status: "verified"`, `llm_checked: True` on prose no adversary had ever seen. "Checked and
clean" and "could not be checked" are different results and only one of them may publish.

**It only ever looked at the first match.** `honesty_check` used `rx.search`, decided the rule on
that one occurrence, and moved on. The scoped block-rules downgrade to a warning when the
sentence reads as being about a third party -- so one third-party sentence placed ahead of your
own overclaim laundered the real claim, and the answer shipped with a warning instead of being
withheld.
"""
import pytest

from hcls_common.verify_gate import (
    ClaimCheckUnavailable,
    extract_claims,
    honesty_check,
    is_publishable,
    verify_claims,
    verify_text,
)


class ExplodingLLM:
    """Every call fails, as in an API outage, a bad key, or a network partition."""

    def __init__(self, exc=RuntimeError("upstream 503")):
        self.exc = exc
        self.calls = 0

    def generate_json(self, *a, **k):
        self.calls += 1
        raise self.exc


class LLMThatExtractsThenDies:
    """Extraction succeeds; the adversarial pass fails. The dangerous half-way state."""

    def __init__(self):
        self.calls = 0

    def generate_json(self, *a, **k):
        self.calls += 1
        if self.calls == 1:
            return {"claims": ["Everolimus reduced SEGA volume by 50% in this patient."]}
        raise RuntimeError("upstream 503")


class LLMReturningJunk:
    def __init__(self, payload):
        self.payload = payload
        self.calls = 0

    def generate_json(self, *a, **k):
        self.calls += 1
        if self.calls == 1:
            return {"claims": ["A checkable claim."]}
        return self.payload


CLINICAL = ("The patient's variant is pathogenic and the recommended dose is 10 mg daily. "
            "Research use only; decision support for a qualified clinician.")


class TestAnOutageCannotProduceAVerifiedVerdict:
    def test_status_is_not_verified_when_the_check_could_not_run(self):
        v = verify_text(CLINICAL, evidence=[{"text": "some evidence"}], llm=ExplodingLLM())
        assert v["status"] != "verified"
        assert v["status"] == "draft"

    def test_it_is_not_publishable(self):
        v = verify_text(CLINICAL, evidence=[{"text": "e"}], llm=ExplodingLLM())
        assert is_publishable(v) is False

    def test_llm_checked_does_not_claim_a_check_that_never_happened(self):
        v = verify_text(CLINICAL, evidence=[{"text": "e"}], llm=ExplodingLLM())
        assert v["summary"]["llm_checked"] is False
        assert v["summary"]["claim_check_unavailable"] is True

    def test_the_reason_is_reported(self):
        v = verify_text(CLINICAL, evidence=[{"text": "e"}], llm=ExplodingLLM())
        assert "503" in (v["summary"]["claim_check_error"] or "")

    def test_a_failure_halfway_through_also_fails_closed(self):
        """Claims extracted, then the adversarial pass dies: the claims are UNCHECKED."""
        v = verify_text(CLINICAL, evidence=[{"text": "e"}], llm=LLMThatExtractsThenDies())
        assert v["status"] == "draft"
        assert v["summary"]["llm_checked"] is False
        assert v["flagged"] and v["flagged"][0]["verdict"] == "unchecked"

    def test_extract_claims_raises_rather_than_reporting_no_claims(self):
        with pytest.raises(ClaimCheckUnavailable):
            extract_claims("some text", ExplodingLLM())

    def test_verify_claims_reports_unavailable_not_clean(self):
        r = verify_claims(CLINICAL, [{"text": "e"}], ExplodingLLM())
        assert r["available"] is False
        assert r["error"]


class TestARealCheckStillWorks:
    def test_text_with_no_checkable_claims_is_available_and_can_verify(self):
        """"Nothing to check" is a genuine clean result and must stay distinguishable."""
        class NoClaims:
            def generate_json(self, *a, **k):
                return {"claims": []}

        r = verify_claims("Research use only; decision support.", None, NoClaims())
        assert r["available"] is True
        assert r["error"] is None

    def test_a_supported_claim_verifies(self):
        class Supportive:
            def __init__(self):
                self.calls = 0

            def generate_json(self, *a, **k):
                self.calls += 1
                if self.calls == 1:
                    return {"claims": ["The dose is 10 mg."]}
                return {"verdicts": [{"claim_index": 1, "verdict": "supported", "why": "p.2"}]}

        v = verify_text(CLINICAL, evidence=[{"text": "The dose is 10 mg."}], llm=Supportive())
        assert v["status"] == "verified"
        assert v["summary"]["llm_checked"] is True
        assert v["summary"]["claim_check_unavailable"] is False

    def test_no_llm_means_no_claim_check_was_requested(self):
        v = verify_text("Research use only; decision support.", llm=None)
        assert v["summary"]["llm_checked"] is False
        assert v["summary"]["claim_check_unavailable"] is False


class TestAMalformedModelResponseIsUncheckedNotCrashed:
    @pytest.mark.parametrize("payload", [
        {"verdicts": "none"},                                  # not a list
        {"verdicts": ["not a dict"]},                          # entries not dicts
        {"verdicts": [{"claim_index": "two", "verdict": "supported"}]},   # unparseable index
        {},                                                    # no verdicts key
    ])
    def test_it_does_not_raise(self, payload):
        v = verify_text(CLINICAL, evidence=[{"text": "e"}], llm=LLMReturningJunk(payload))
        assert v["status"] in {"draft", "blocked"}
        assert v["status"] != "verified"


class TestEveryOccurrenceIsScanned:
    SELF_CLAIM = "This platform is FDA cleared for diagnostic use."

    def test_a_self_referential_block_claim_is_blocked_on_its_own(self):
        out = honesty_check(self.SELF_CLAIM)
        assert any(v["severity"] == "block" for v in out), out

    def test_a_third_party_sentence_first_does_not_launder_it(self):
        """The bypass: the first match scopes to a third party, the real claim comes later."""
        text = ("Everolimus is FDA approved for SEGA in tuberous sclerosis. "
                + self.SELF_CLAIM)
        out = honesty_check(text)
        assert any(v["severity"] == "block" for v in out), (
            f"the self-referential claim was laundered to {[v['severity'] for v in out]}")

    def test_the_excerpt_names_the_sentence_that_blocked(self):
        text = ("Everolimus is FDA approved for SEGA in tuberous sclerosis. "
                + self.SELF_CLAIM)
        blocking = [v for v in honesty_check(text) if v["severity"] == "block"]
        assert blocking
        assert "platform" in blocking[0]["excerpt"] or "This" in blocking[0]["excerpt"]

    def test_a_third_party_statement_alone_still_only_warns(self):
        """The scoping exists for a reason and must survive the fix."""
        out = honesty_check("Everolimus is FDA approved for SEGA in tuberous sclerosis.")
        assert not any(v["severity"] == "block" for v in out)

    def test_verify_text_blocks_the_laundered_case_end_to_end(self):
        text = ("Everolimus is FDA approved for SEGA in tuberous sclerosis. "
                + self.SELF_CLAIM)
        v = verify_text(text, llm=None)
        assert v["status"] == "blocked"
        assert is_publishable(v) is False
