"""A new route that streams generated prose must not quietly skip the gate.

The two that existed (`/query/stream` on the precision-biomarker and precision-autoimmune
agents) went unnoticed for months because the response middleware skips `text/event-stream` and
says so in a header nobody reads. This scans the API sources for the combination that defines
the defect -- an SSE response in a module that also pulls tokens from an LLM -- and requires
that module to use `gate_token_stream`.

It is a source scan, which is crude. The alternative is noticing the next one by accident.
"""
import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parents[3]

pytestmark = pytest.mark.skipif(not (REPO / "core").is_dir(),
                                reason="run from a source checkout")

# Pulling tokens from a model, as opposed to emitting status events or file bytes.
TOKEN_SOURCES = ("query_stream", "generate_stream", "text_stream", "chat_stream")
SSE = ("text/event-stream", "application/x-ndjson")


def api_modules():
    for base in (REPO / "core").glob("*/*/api"):
        yield from base.rglob("*.py")
    for base in (REPO / "core").glob("*/*/*/api"):
        yield from base.rglob("*.py")


def prose_streaming_modules():
    out = []
    for path in sorted(set(api_modules())):
        src = path.read_text(errors="replace")
        if any(s in src for s in SSE) and any(t in src for t in TOKEN_SOURCES):
            out.append((path.relative_to(REPO).as_posix(), src))
    return out


def test_the_scan_finds_the_two_known_routes():
    """If this breaks, the scan stopped matching and the guard below is vacuous."""
    found = {p for p, _ in prose_streaming_modules()}
    assert "core/agents/precision-biomarker/api/routes/analysis.py" in found
    assert "core/agents/precision-autoimmune/api/main.py" in found


def test_every_prose_stream_is_gated():
    ungated = [p for p, src in prose_streaming_modules() if "gate_token_stream" not in src]
    assert not ungated, (
        "these modules stream generated prose over SSE without gate_token_stream, so the "
        "answer reaches the caller with no overclaim register and no disclaimer: " + repr(ungated)
    )


def test_every_prose_stream_declares_in_stream_gating_in_its_header():
    """The middleware marks streams `not-scanned`; a gated one must correct that claim."""
    missing = [p for p, src in prose_streaming_modules() if "in-stream" not in src]
    assert not missing, f"gated prose streams not advertising it in X-HCLS-Honesty: {missing}"


def test_the_middleware_does_not_overwrite_a_route_header():
    """A route's own X-HCLS-Honesty must survive; otherwise the claim above is a lie."""
    gate_src = (REPO / "lib" / "hcls_common" / "api_gate.py").read_text()
    # Anchor on the code branch, not on the docstring that also names the media type.
    block = re.search(r'if any\(s in ctype for s in \([\s\S]{0,500}', gate_src)
    assert block, "the streaming branch moved; re-point this test"
    assert "setdefault" in block.group(0), (
        "the streaming branch assigns X-HCLS-Honesty directly, which would overwrite the "
        "header a self-gating route sets"
    )
