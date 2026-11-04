"""Security regression tests for the TSC clinical API (2026-11-03 audit, S6 and S7).

Two defects, both on an API that is open by default:

* CORS was `allow_origins=["*"]`, so any page a clinician happened to visit could read
  every cohort, agent, provenance and briefing route from their machine, with their
  browser as the proxy. Harm is bounded today only because the cohort is synthetic.
* The API-key check compared with `!=`, which leaks the key's length and matching prefix
  through response timing.

These tests assert the refusals, and that the wildcard cannot come back by default.
"""
import importlib

from fastapi.testclient import TestClient

import api.main as main_mod
from config.settings import settings


def _origins_of(app):
    for mw in app.user_middleware:
        if mw.cls.__name__ == "CORSMiddleware":
            return mw.kwargs.get("allow_origins"), mw.kwargs
    return None, {}


def test_cors_default_is_not_a_wildcard():
    app = main_mod.create_app()
    origins, _ = _origins_of(app)
    assert origins is not None, "CORS middleware is not installed"
    assert "*" not in origins


def test_cors_default_is_a_loopback_allowlist():
    app = main_mod.create_app()
    origins, _ = _origins_of(app)
    assert origins
    assert all(o.startswith("http://localhost") or o.startswith("http://127.0.0.1")
               for o in origins), origins


def test_cors_credentials_are_not_enabled():
    """A wildcard origin with credentials is the cookie-theft shape; keep them apart."""
    app = main_mod.create_app()
    _, kwargs = _origins_of(app)
    assert kwargs.get("allow_credentials") is False


def test_cors_headers_are_not_a_wildcard():
    app = main_mod.create_app()
    _, kwargs = _origins_of(app)
    assert "*" not in kwargs.get("allow_headers", [])
    assert "*" not in kwargs.get("allow_methods", [])


def test_foreign_origin_is_not_echoed_back():
    app = main_mod.create_app()
    with TestClient(app) as c:
        r = c.get("/health", headers={"Origin": "https://evil.example"})
    assert r.headers.get("access-control-allow-origin") != "https://evil.example"


def test_allowlisted_origin_is_still_served():
    app = main_mod.create_app()
    with TestClient(app) as c:
        r = c.get("/health", headers={"Origin": "http://localhost:8080"})
    assert r.status_code == 200
    assert r.headers.get("access-control-allow-origin") == "http://localhost:8080"


def test_cors_allowlist_is_configurable(monkeypatch):
    monkeypatch.setattr(settings, "CORS_ORIGINS", "https://tsc.example.org")
    app = main_mod.create_app()
    origins, _ = _origins_of(app)
    assert origins == ["https://tsc.example.org"]


def test_api_key_comparison_is_constant_time(monkeypatch):
    """Guard against a regression back to `!=`."""
    import inspect
    source = inspect.getsource(main_mod.create_app)
    assert "keys_match" in source
    assert 'request.headers.get("X-API-Key") != settings.API_KEY' not in source


def test_gated_api_refuses_wrong_and_missing_key(monkeypatch):
    monkeypatch.setattr(settings, "API_KEY", "tsc-test-key")
    app = main_mod.create_app()
    with TestClient(app) as c:
        assert c.get("/cohort").status_code == 401
        assert c.get("/cohort", headers={"X-API-Key": "wrong"}).status_code == 401
        assert c.get("/cohort", headers={"X-API-Key": "tsc-test-key"}).status_code == 200
        # Probes must not fail closed with the service.
        assert c.get("/health").status_code == 200


def test_gated_api_does_not_500_on_a_non_ascii_header(monkeypatch):
    """compare_digest on str raises on non-ASCII; the bytes comparison returns a verdict."""
    monkeypatch.setattr(settings, "API_KEY", "tsc-test-key")
    app = main_mod.create_app()
    with TestClient(app) as c:
        r = c.get("/cohort", headers={"X-API-Key": "pässwörd".encode("utf-8")})
    assert r.status_code == 401


def test_settings_module_imports_clean():
    """The allowlist lives in settings, so a typo there must not be silent."""
    import config.settings as s
    importlib.reload(s)
    assert s.settings.cors_origins()
