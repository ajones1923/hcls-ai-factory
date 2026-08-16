"""The governance header must not claim a gate that did not run.

Before the 2026-08-16 fix, `X-HCLS-Governed` was emitted on every response regardless of whether
any gate executed. These tests fail if that regresses.
"""
from fastapi import FastAPI
from fastapi.testclient import TestClient
from hcls_common.api_gate import install_governance, honesty_flags


def build():
    app = FastAPI()
    install_governance(app, service="testsvc", capability_id="cart-intelligence-agent")

    @app.get("/ungated")
    def ungated():
        return {"answer": "no gate ran here"}

    @app.get("/gated")
    def gated():
        honesty_flags("This is decision support for a qualified clinician.")
        return {"answer": "a gate ran"}

    return TestClient(app)


def test_no_governed_header_when_no_gate_ran():
    r = build().get("/ungated")
    assert r.status_code == 200
    assert "X-HCLS-Governed" not in r.headers        # the whole point
    assert r.headers["X-HCLS-Service"] == "testsvc"  # true, and still useful


def test_governed_header_lists_the_gate_that_ran():
    r = build().get("/gated")
    assert r.headers.get("X-HCLS-Governed") == "output-honesty"


def test_header_does_not_leak_between_requests():
    c = build()
    c.get("/gated")
    assert "X-HCLS-Governed" not in c.get("/ungated").headers


def test_governance_endpoint_states_gates_are_opt_in():
    body = build().get("/governance").json()
    assert "gates_are_opt_in" in body
    assert "does NOT gate" in body["gates_are_opt_in"]
    assert "auth" in body
