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


def test_governance_endpoint_describes_both_gates_accurately():
    """The OUTPUT gate stopped being opt-in on 2026-09-15; the input gate did not.

    The endpoint must not overstate either. Claiming automatic input validation would be the
    same class of error the X-HCLS-Governed header fix above exists to prevent.
    """
    body = build().get("/governance").json()
    assert "AUTOMATIC" in body["output_gate"]
    assert "must still call require_valid_input" in body["input_gate_is_opt_in"]
    assert "auth" in body


# ── the output gate is automatic (2026-09-15) ────────────────────────────────
def build_clinical():
    app = FastAPI()
    install_governance(app, service="testsvc", capability_id="cart-intelligence-agent")

    @app.get("/answer")
    def answer():
        # must contain a _CLINICAL_TRIGGER word (patient / therapy / dose / variant ...)
        # or the missing-disclaimer rule correctly does not fire.
        return {"answer": "Tisagenlecleucel therapy in a paediatric patient achieved an 81% "
                          "overall remission rate. " + "Dosing follows the label. " * 8}

    @app.get("/overclaim")
    def overclaim():
        return {"answer": "This therapy is clinically proven to cure the disease and is "
                          "FDA-approved for all patients. " + "It carries zero risk. " * 10}

    @app.get("/short")
    def short():
        return {"answer": "ok"}          # not clinical prose

    return TestClient(app)


def test_clinical_answer_gains_the_decision_support_disclaimer():
    body = build_clinical().get("/answer").json()
    assert "decision support for a qualified clinician" in body["answer"].lower()
    assert body["honesty"], "missing-disclaimer finding should be reported"


def test_short_answers_are_left_alone():
    """A status string is not clinical prose; appending a disclaimer would be noise."""
    body = build_clinical().get("/short").json()
    assert body["answer"] == "ok" and "honesty" not in body


def test_overclaims_are_surfaced_not_hidden():
    body = build_clinical().get("/overclaim").json()
    msgs = " ".join(f["message"] for f in body["honesty"])
    assert any(w in msgs for w in ("cure", "Cure", "Regulatory", "overclaim", "Efficacy"))
    assert any(f["severity"] == "block" for f in body["honesty"])
