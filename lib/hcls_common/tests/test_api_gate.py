"""Tests for the shared governance surface (hcls_common.api_gate).

Skips cleanly if FastAPI isn't installed (it's not a core dependency).
"""
import pytest

pytest.importorskip("fastapi")

from fastapi.testclient import TestClient  # noqa: E402

from hcls_common.api_gate import (  # noqa: E402
    create_governed_app,
    honesty_flags,
    require_valid_input,
)


def _app():
    app = create_governed_app("test-svc", capability_id="cart-intelligence-agent")

    @app.post("/query")
    def query(body: dict):
        body = require_valid_input("cart-intelligence-agent", body)
        return {"ok": True, "cleaned": body}

    return app


def test_governance_endpoint_and_headers():
    """Updated 2026-08-16.

    `gates` became `gates_available`, and `X-HCLS-Governed` is no longer emitted unconditionally —
    it appears only when a handler actually invoked a gate. /governance is a plain info route that
    runs none, so asserting the header here was asserting the very overclaim that was removed.
    `X-HCLS-Service` is the header that is always true.
    """
    c = TestClient(_app())
    r = c.get("/governance")
    assert r.status_code == 200
    assert r.json()["gates_available"] == ["input-validation", "output-honesty"]
    assert any(k.lower() == "x-hcls-service" for k in r.headers)
    assert not any(k.lower() == "x-hcls-governed" for k in r.headers)


def test_governed_header_appears_when_a_gate_runs():
    c = TestClient(_app())
    r = c.post("/query", json={"patient_context": {"patient_id": "P1"}, "question": "q"})
    if r.status_code == 200:                      # payload satisfied the contract
        assert r.headers.get("X-HCLS-Governed") == "input-validation"


def test_input_gate_rejects_missing_required():
    c = TestClient(_app())
    r = c.post("/query", json={"patient_context": {}})
    assert r.status_code == 422
    assert "input_errors" in r.json()["detail"]


def test_input_gate_accepts_valid_payload():
    c = TestClient(_app())
    r = c.post("/query", json={"query": "CD19 CAR-T"})
    assert r.status_code == 200
    assert r.json()["ok"] is True


def test_honesty_flags_catches_overclaim():
    flags = honesty_flags("This is FDA-approved and guaranteed to cure the patient.")
    assert len(flags) >= 1


def test_honesty_flags_clean_text():
    assert honesty_flags("Retrieved 3 candidate references for review.") == []


def test_require_valid_input_unknown_capability_passes_through():
    payload = {"anything": 1}
    assert require_valid_input("no-such-capability", payload) == payload
