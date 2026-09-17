# Copyright 2026 Adam Jones
# SPDX-License-Identifier: Apache-2.0
# Part of the HCLS AI Factory — https://github.com/ajones1923/hcls-ai-factory
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
"""Extract the answer text from an Anthropic message.

`response.content[0].text` is the obvious line and it is wrong. The content list is a sequence
of BLOCKS, and the first one is not necessarily the answer: with extended thinking enabled the
first block is a `ThinkingBlock`, which has no `.text` at all. The result is:

    AttributeError: 'ThinkingBlock' object has no attribute 'text'

That is not hypothetical. The clinical-imaging engine hit exactly this on 2026-09-16, and because
its generation path catches the exception and falls through to a MOCK response, the live endpoint
answered clinical questions with a canned "normal study" radiology report — 200 OK, fluent prose,
no relation to the question. The output-honesty gate is what surfaced it, by refuting a claim in
the mock that the retrieved evidence contradicted.

Fifteen call sites in this repo took `content[0].text`. One line, one shared implementation, so
the next block-type addition does not silently degrade a clinical answer.
"""
from __future__ import annotations

from typing import Any


def first_text(message: Any, default: str = "") -> str:
    """Return the first text block's text from an Anthropic message.

    Skips thinking, redacted-thinking, tool-use and any other non-text block. Accepts a raw
    content list as well as a message object, and never raises on an unexpected shape — a
    generation path that raises here degrades to whatever its fallback is, which is precisely
    the failure this exists to stop.
    """
    content = getattr(message, "content", message)
    if isinstance(content, str):
        return content
    if not isinstance(content, (list, tuple)):
        return default
    for block in content:
        if getattr(block, "type", None) == "text" and isinstance(getattr(block, "text", None), str):
            return block.text
        if isinstance(block, dict) and block.get("type") == "text" and isinstance(block.get("text"), str):
            return block["text"]
    # No typed text block: accept a bare `.text` attribute, which is what the older
    # single-block responses and most test doubles look like.
    for block in content:
        t = getattr(block, "text", None)
        if isinstance(t, str):
            return t
        if isinstance(block, dict) and isinstance(block.get("text"), str):
            return block["text"]
    return default
