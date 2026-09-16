"""Tests for first-party API call headers.

Six UIs called their own agent's API with no credential after the gate went fail-closed. Every
query returned 401 while the UI still loaded and reported nothing — the failure was one layer in.
"""
from hcls_common.api_client import auth_headers


class TestAuthHeaders:
    def test_key_is_attached(self, monkeypatch):
        monkeypatch.setenv("HCLS_API_KEY", "shared-key")
        assert auth_headers()["X-API-Key"] == "shared-key"

    def test_extra_headers_are_preserved(self, monkeypatch):
        monkeypatch.setenv("HCLS_API_KEY", "k")
        h = auth_headers({"Content-Type": "application/json"})
        assert h["Content-Type"] == "application/json" and h["X-API-Key"] == "k"

    def test_per_service_key_wins(self, monkeypatch):
        monkeypatch.setenv("HCLS_API_KEY", "shared")
        monkeypatch.setenv("HCLS_API_KEY_NEUROLOGY", "per-service")
        assert auth_headers(service="neurology")["X-API-Key"] == "per-service"

    def test_service_slug_normalises_hyphens(self, monkeypatch):
        monkeypatch.delenv("HCLS_API_KEY", raising=False)
        monkeypatch.setenv("HCLS_API_KEY_CLINICAL_TRIAL", "ct")
        assert auth_headers(service="clinical-trial")["X-API-Key"] == "ct"

    def test_no_key_configured_leaves_headers_alone(self, monkeypatch):
        monkeypatch.delenv("HCLS_API_KEY", raising=False)
        h = auth_headers({"Content-Type": "application/json"})
        assert h == {"Content-Type": "application/json"}, (
            "a deployment with the gate disabled must keep working")

    def test_caller_dict_is_not_mutated(self, monkeypatch):
        monkeypatch.setenv("HCLS_API_KEY", "k")
        base = {"Content-Type": "application/json"}
        auth_headers(base)
        assert "X-API-Key" not in base
