"""Tests for answer-text extraction.

`content[0].text` was taken at fifteen call sites. With extended thinking the first block is a
ThinkingBlock with no `.text`, and in clinical-imaging the resulting AttributeError was caught by
a fallback chain that served a MOCK radiology report from a live clinical endpoint.
"""
from hcls_common.llm_text import first_text


class _Block:
    def __init__(self, type=None, text=None):
        if type is not None:
            self.type = type
        if text is not None:
            self.text = text


class _Msg:
    def __init__(self, content):
        self.content = content


class TestThinkingBlocks:
    def test_thinking_block_first_is_skipped(self):
        """The exact production failure: ThinkingBlock has no .text at all."""
        msg = _Msg([_Block(type="thinking"), _Block(type="text", text="CAD-RADS 4A")])
        assert first_text(msg) == "CAD-RADS 4A"

    def test_redacted_thinking_is_skipped(self):
        msg = _Msg([_Block(type="redacted_thinking"), _Block(type="text", text="answer")])
        assert first_text(msg) == "answer"

    def test_tool_use_block_is_skipped(self):
        msg = _Msg([_Block(type="tool_use"), _Block(type="text", text="answer")])
        assert first_text(msg) == "answer"

    def test_first_text_block_wins_over_later_ones(self):
        msg = _Msg([_Block(type="thinking"), _Block(type="text", text="first"),
                    _Block(type="text", text="second")])
        assert first_text(msg) == "first"


class TestOrdinaryShapes:
    def test_single_text_block(self):
        assert first_text(_Msg([_Block(type="text", text="hello")])) == "hello"

    def test_block_without_a_type_but_with_text(self):
        """Most test doubles in this repo look like this."""
        assert first_text(_Msg([_Block(text="hello")])) == "hello"

    def test_dict_blocks(self):
        assert first_text(_Msg([{"type": "thinking"}, {"type": "text", "text": "hi"}])) == "hi"

    def test_a_raw_content_list_is_accepted(self):
        assert first_text([_Block(type="text", text="hi")]) == "hi"

    def test_a_plain_string_passes_through(self):
        assert first_text("hi") == "hi"


class TestNeverRaises:
    """A generation path that raises here degrades to its fallback — the failure this prevents."""

    def test_empty_content(self):
        assert first_text(_Msg([])) == ""

    def test_no_text_anywhere(self):
        assert first_text(_Msg([_Block(type="thinking"), _Block(type="tool_use")])) == ""

    def test_unexpected_shape(self):
        assert first_text(None) == ""
        assert first_text(_Msg(42)) == ""

    def test_default_is_returned(self):
        assert first_text(_Msg([]), default="(no answer)") == "(no answer)"
