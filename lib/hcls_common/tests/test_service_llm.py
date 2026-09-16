"""Tests for the one shared service LLM client (Phase 2 de-duplication).

These exist because the rule they cover was previously duplicated in nine places and got
applied correctly in one of them. The sampling-parameter rule is now testable in CI instead
of discoverable only by a 400 in production that degrades to a stub answer.
"""
import pytest

from hcls_common.service_llm import ServiceLLM, build_service_llm, rejects_sampling


class _FakeMessages:
    def __init__(self):
        self.last = None

    def create(self, **kw):
        self.last = kw
        return type("Msg", (), {"content": [type("B", (), {"text": "answer"})()]})()


class _FakeClient:
    def __init__(self):
        self.messages = _FakeMessages()


class TestSamplingRule:
    @pytest.mark.parametrize("model", [
        "claude-opus-5", "claude-sonnet-5", "claude-opus-4-8", "claude-opus-4-7",
        "claude-opus-4-6", "claude-sonnet-4-6", "claude-fable-5", "CLAUDE-OPUS-5",
    ])
    def test_current_models_reject_sampling(self, model):
        assert rejects_sampling(model)

    @pytest.mark.parametrize("model", [
        "claude-3-5-sonnet-20241022", "claude-haiku-4-5-20251001", "", None,
    ])
    def test_older_models_still_accept_it(self, model):
        assert not rejects_sampling(model)


class TestKwargs:
    def test_temperature_omitted_on_current_models(self):
        c = _FakeClient()
        ServiceLLM("claude-opus-5", client=c).generate("q", temperature=0.7)
        assert "temperature" not in c.messages.last, (
            "sending temperature to a current model returns 400 and the caller "
            "silently degrades to a stub answer")

    def test_temperature_forwarded_where_accepted(self):
        c = _FakeClient()
        ServiceLLM("claude-haiku-4-5-20251001", client=c).generate("q", temperature=0.2)
        assert c.messages.last["temperature"] == 0.2

    def test_explicit_none_never_sends_it(self):
        c = _FakeClient()
        ServiceLLM("claude-haiku-4-5-20251001", client=c).generate("q", temperature=None)
        assert "temperature" not in c.messages.last


class TestSystemPrompt:
    """The default system prompt is the only genuine difference among the eight services."""

    def test_default_system_used_when_caller_gives_none(self):
        c = _FakeClient()
        ServiceLLM("claude-opus-5", client=c, default_system="You are a CAR-T agent.").generate("q")
        assert c.messages.last["system"] == "You are a CAR-T agent."

    def test_per_call_system_prompt_wins(self):
        c = _FakeClient()
        ServiceLLM("claude-opus-5", client=c, default_system="default").generate("q", "per-call")
        assert c.messages.last["system"] == "per-call"

    def test_system_is_always_a_string(self):
        c = _FakeClient()
        ServiceLLM("claude-opus-5", client=c).generate("q")
        assert c.messages.last["system"] == ""


class TestGenerate:
    def test_returns_first_text_block(self):
        assert ServiceLLM("claude-opus-5", client=_FakeClient()).generate("q") == "answer"

    def test_model_and_prompt_are_passed_through(self):
        c = _FakeClient()
        ServiceLLM("claude-opus-5", client=c).generate("what is CD19?", max_tokens=99)
        assert c.messages.last["model"] == "claude-opus-5"
        assert c.messages.last["max_tokens"] == 99
        assert c.messages.last["messages"] == [{"role": "user", "content": "what is CD19?"}]


class TestBuildServiceLLM:
    def test_missing_sdk_or_credential_degrades_to_none(self, monkeypatch):
        """A service must still start and serve retrieval when synthesis is unavailable."""
        import hcls_common.service_llm as mod

        def _boom(*a, **k):
            raise RuntimeError("no API key")

        monkeypatch.setattr(mod, "ServiceLLM", _boom)
        assert build_service_llm("claude-opus-5", service="cart") is None

    def test_failure_is_logged_not_silent(self, monkeypatch, caplog):
        """A silent None is how 'healthy but answering with a stub' happened."""
        import hcls_common.service_llm as mod

        def _boom(*a, **k):
            raise RuntimeError("no API key")

        monkeypatch.setattr(mod, "ServiceLLM", _boom)
        with caplog.at_level("WARNING"):
            build_service_llm("claude-opus-5", service="cart")
        assert "cart" in caplog.text
