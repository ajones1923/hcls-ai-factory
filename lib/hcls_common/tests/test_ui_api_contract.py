"""Every UI call to a governed API must carry a credential.

On 2026-09-16 six Streamlit UIs were returning 401 on every query. They called their own agent's
API with `headers={"Content-Type": "application/json"}` and nothing else, and the moment the gate
went fail-closed each of them stopped working. Nothing noticed, because the entire verification
estate — 8,191 tests, 27 clinical eval cases, 17 demos — talks to API ports:

    ports exercised by all 17 demo transcripts:  8127 8508 8522 8524 8527 8529 8532 8536 8539 8541 8545
    UI ports exercised by anything:              none

This test closes that specific hole without needing a browser or a running service: it reads each
UI with the AST and asserts that any request aimed at the agent API passes headers. It is a
contract check, not a substitute for driving the UI — but it makes the regression that actually
happened impossible to reintroduce silently.
"""
from __future__ import annotations

import ast
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[3]
UI_FILES = sorted(ROOT.glob("core/*/*/app/*_ui.py")) + sorted(ROOT.glob("core/*/*/*/app/*_ui.py"))


def _api_calls(tree: ast.AST) -> list[ast.Call]:
    """requests.get/post/put calls whose URL mentions the agent API base."""
    out = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        f = node.func
        if not (isinstance(f, ast.Attribute) and f.attr in {"get", "post", "put", "patch"}):
            continue
        if not (isinstance(f.value, ast.Name) and f.value.id == "requests"):
            continue
        url = ast.dump(node.args[0]) if node.args else ""
        if "API_BASE" in url:            # only the governed agent API; Orthanc etc. carry their own auth
            out.append(node)
    return out


def test_ui_files_are_discovered():
    assert UI_FILES, "no *_ui.py found — the glob is wrong and this test would vacuously pass"


@pytest.mark.parametrize("path", UI_FILES, ids=lambda p: p.name)
def test_every_api_call_carries_a_credential(path: pathlib.Path):
    """Per CALL SITE, not per file.

    A first version of this test asserted only that a `headers=` kwarg existed and that the file
    mentioned auth_headers somewhere. Reintroducing the real bug in ONE of two call sites passed
    it — `headers={"Content-Type": "application/json"}` satisfies both, and the other call site
    kept the word `auth_headers` in the file. A test that passes while the defect is present is
    worse than no test, so this checks the headers expression of each individual call.
    """
    tree = ast.parse(path.read_text(errors="ignore"))
    calls = _api_calls(tree)
    if not calls:
        pytest.skip("does not call the agent API")
    for call in calls:
        headers = next((k.value for k in call.keywords if k.arg == "headers"), None)
        assert headers is not None, (
            f"{path.name}:{call.lineno} calls the agent API with no headers at all — the gate is "
            f"fail-closed, so this is a 401.")
        uses_auth = any(
            isinstance(n, ast.Call) and (
                (isinstance(n.func, ast.Name) and n.func.id == "auth_headers")
                or (isinstance(n.func, ast.Attribute) and n.func.attr == "auth_headers"))
            for n in ast.walk(headers))
        assert uses_auth, (
            f"{path.name}:{call.lineno} passes headers that never reach auth_headers(). "
            f"That is a 401 on every query, and the UI will report nothing until someone asks it "
            f"a question. Use hcls_common.api_client.auth_headers().")
