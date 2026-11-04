"""Security regression tests for this service's API-key gate (2026-11-03 audit, S7).

The gate had two faults, both of which this file pins down:

* it compared the supplied key with `!=`, which leaks the key's length and its matching
  prefix through response timing (and, had it been compare_digest on str, would have raised
  on a non-ASCII header);
* it read only this service's own API_KEY setting, while the shared gate installed a few
  lines above reads HCLS_API_KEY, and the two carried different open-path sets. So which
  endpoints were credentialed depended on which setting the operator happened to use, and
  /governance reported on only one of the two.

The gate now resolves one key (shared first, legacy setting as fallback), compares it in
constant time, and leaves open the same paths the shared gate does.
"""
import inspect

import pytest
from fastapi.testclient import TestClient

import api.main as main_mod
from config.settings import settings


GATED_PATH = "/__gate_probe__"


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for var in ("HCLS_API_KEY", "HCLS_API_KEY_RARE_DISEASE_DIAGNOSTIC"):
        monkeypatch.delenv(var, raising=False)
    original = getattr(settings, "API_KEY", None)
    # These services rate-limit per client IP in a per-process dict. The whole subject suite
    # shares that process, so a test asserting a 401 must not be answered with a 429 left
    # over from earlier tests.
    store = getattr(main_mod, "_rate_limit_store", None)
    if store is not None:
        store.clear()
    yield
    settings.API_KEY = original


def test_comparison_is_constant_time():
    """Guard against a regression back to `!=`."""
    source = inspect.getsource(main_mod)
    assert "_api_keys_match(" in source
    assert "provided != api_key" not in source
    assert "provided_key != api_key" not in source
    assert "compare_digest" in inspect.getsource(main_mod._api_keys_match)


def test_match_never_raises_and_rejects_mismatches():
    assert main_mod._api_keys_match("k", "k") is True
    assert main_mod._api_keys_match("k", "K") is False
    assert main_mod._api_keys_match("", "k") is False
    assert main_mod._api_keys_match("k", "") is False
    # A non-ASCII header used to be a TypeError inside the middleware.
    assert main_mod._api_keys_match("p\u00e4ssw\u00f6rd", "p\u00e4ssw\u00f6rd") is True
    assert main_mod._api_keys_match("p\u00e4ssw\u00f6rd", "other") is False


def test_shared_key_setting_enables_this_gate(monkeypatch):
    """HCLS_API_KEY used to enable only the shared gate, leaving this one off."""
    monkeypatch.setenv("HCLS_API_KEY", "shared-secret")
    settings.API_KEY = None
    assert main_mod._resolved_api_key() == "shared-secret"


def test_per_service_shared_key_wins(monkeypatch):
    monkeypatch.setenv("HCLS_API_KEY", "shared-secret")
    monkeypatch.setenv("HCLS_API_KEY_RARE_DISEASE_DIAGNOSTIC", "service-secret")
    assert main_mod._resolved_api_key() == "service-secret"


def test_legacy_setting_still_enables_the_gate(monkeypatch):
    """Operators who configured only the per-service setting must not lose their gate."""
    settings.API_KEY = "legacy-secret"
    assert main_mod._resolved_api_key() == "legacy-secret"


def test_both_gates_leave_the_same_paths_open():
    from hcls_common.api_auth import DEFAULT_OPEN_PATHS
    open_paths = main_mod._auth_open_paths({"/health"})
    for path in DEFAULT_OPEN_PATHS:
        assert path in open_paths, path


def test_gated_route_refuses_missing_and_wrong_key(monkeypatch):
    monkeypatch.setenv("HCLS_API_KEY", "shared-secret")
    with TestClient(main_mod.app) as c:
        assert c.get(GATED_PATH).status_code == 401
        assert c.get(GATED_PATH, headers={"X-API-Key": "wrong"}).status_code == 401
        ok = c.get(GATED_PATH, headers={"X-API-Key": "shared-secret"})
        assert ok.status_code != 401


def test_probes_stay_open_when_the_gate_is_on(monkeypatch):
    """A gate that fails closed with the service takes the whole fleet down."""
    monkeypatch.setenv("HCLS_API_KEY", "shared-secret")
    with TestClient(main_mod.app) as c:
        assert c.get("/health").status_code == 200


def test_non_ascii_header_is_401_not_500(monkeypatch):
    monkeypatch.setenv("HCLS_API_KEY", "shared-secret")
    with TestClient(main_mod.app) as c:
        r = c.get(GATED_PATH, headers={"X-API-Key": "p\u00e4ssw\u00f6rd".encode("utf-8")})
    assert r.status_code == 401
