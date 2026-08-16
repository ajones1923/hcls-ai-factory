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
