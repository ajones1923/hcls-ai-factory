"""Negative-tested: a gate that never rejects is worse than no gate."""
import os, pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from hcls_common.api_auth import install_api_key_auth, resolve_key, auth_status

def build(**env):
    for k in list(os.environ):
        if k.startswith("HCLS_API_KEY"): del os.environ[k]
    os.environ.update(env)
    app = FastAPI()
    @app.get("/health")
    def h(): return {"ok": True}
    @app.get("/interpret")
    def i(): return {"variant": "VCP", "note": "decision support"}
    installed = install_api_key_auth(app, service="cart")
    return TestClient(app), installed

def test_disabled_by_default():
    c, installed = build()
    assert installed is False
    assert c.get("/interpret").status_code == 200      # unchanged posture

def test_fail_closed_when_key_set():
    c, installed = build(HCLS_API_KEY="s3cret")
    assert installed is True
    assert c.get("/interpret").status_code == 401      # no key -> rejected
    assert c.get("/interpret", headers={"X-API-Key": "wrong"}).status_code == 401
    assert c.get("/interpret", headers={"X-API-Key": "s3cret"}).status_code == 200

def test_health_stays_public():
    c, _ = build(HCLS_API_KEY="s3cret")
    assert c.get("/health").status_code == 200         # probes must not fail closed

def test_per_service_key_overrides_shared():
    c, _ = build(HCLS_API_KEY="shared", HCLS_API_KEY_CART="cartkey")
    assert c.get("/interpret", headers={"X-API-Key": "shared"}).status_code == 401
    assert c.get("/interpret", headers={"X-API-Key": "cartkey"}).status_code == 200

def test_status_never_leaks_the_key():
    os.environ["HCLS_API_KEY"] = "topsecret"
    s = auth_status("cart")
    assert s["api_key_required"] is True
    assert "topsecret" not in repr(s)

# ── S8 regression (2026-11-03 audit) ─────────────────────────────────────────────────────
# compare_digest on str raises TypeError unless both sides are ASCII. A non-ASCII header
# was therefore a 500 rather than a 401, and a non-ASCII HCLS_API_KEY took every gated
# route down while /health stayed 200 and auth_status() still said "fail-closed".

def test_non_ascii_header_is_rejected_not_an_error():
    """A non-ASCII byte in the header used to raise TypeError inside the middleware."""
    c, _ = build(HCLS_API_KEY="s3cret")
    # Sent as raw bytes, which is how a non-ASCII header actually arrives; Starlette
    # latin-1-decodes it, so the middleware sees a non-ASCII str.
    r = c.get("/interpret", headers={"X-API-Key": "pässwörd".encode("utf-8")})
    assert r.status_code == 401


def test_non_ascii_configured_key_still_gates():
    """A passphrase with an accent in it used to 500 every gated route on the service.

    Every request now gets a verdict — 401 here — rather than an unhandled exception. The
    matching half is asserted on keys_match directly, because what a non-ASCII key has to
    be matched against depends on how the client encoded the header, and that ambiguity is
    itself a reason to tell operators to use an ASCII key.
    """
    c, _ = build(HCLS_API_KEY="pässwörd-ünicode")
    assert c.get("/interpret").status_code == 401
    assert c.get("/interpret", headers={"X-API-Key": "wrong"}).status_code == 401
    assert c.get("/interpret",
                 headers={"X-API-Key": "pässwörd-ünicode".encode("utf-8")}).status_code == 401


def test_non_ascii_key_does_not_take_the_service_down():
    """The shape of the outage: every gated route 500s while /health stays 200."""
    c, _ = build(HCLS_API_KEY="pässwörd")
    assert c.get("/health").status_code == 200
    assert c.get("/interpret", headers={"X-API-Key": "wrong"}).status_code == 401


def test_keys_match_is_constant_time_and_total():
    import inspect
    from hcls_common.api_auth import keys_match
    assert "compare_digest" in inspect.getsource(keys_match)
    # Never raises, whatever it is handed.
    assert keys_match("", "k") is False
    assert keys_match("k", "") is False
    assert keys_match("ö", "ö") is True
    assert keys_match("ö", "o") is False
    # A lone surrogate is what latin-1/utf-8 round-tripping can leave behind.
    assert keys_match("\udcff", "k") is False
