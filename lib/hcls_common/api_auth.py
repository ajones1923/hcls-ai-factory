"""Shared API-key gate for every service that emits clinical output.

Audit finding (2026-08-15): **1 of 12 FastAPI entrypoints enforced authentication.** The eleven
others — every intelligence agent plus the cardiology, imaging and oncology engines — exposed
clinical decision-support endpoints with no auth dependency at all. They accept patient context and
return variant interpretations, therapy suggestions and trial matches.

This generalises the pattern `core/disease-programs/tuberous-sclerosis` already used (its "P1-7"
gate), so there is ONE implementation rather than twelve.

Posture, deliberately unchanged from TSC's:

    key unset  -> open. The trusted-LAN / synthetic-demo posture the platform ships with.
    key set    -> FAIL CLOSED. Every route except the open paths requires a matching X-API-Key.

Turning it on is a one-line env change per service and needs no code edit:

    HCLS_API_KEY=<shared secret>            # applies to every service
    HCLS_API_KEY_CART=<per-service secret>  # overrides it for the cart agent

Health, docs and the OpenAPI schema stay public so probes, load balancers and the landing page keep
working when the gate is on.
"""
from __future__ import annotations

import hmac
import os
from typing import Iterable

# Probes and API discovery must not require a credential, or health checks fail closed with the
# service and every orchestrator marks it unhealthy.
DEFAULT_OPEN_PATHS: tuple[str, ...] = (
    "/health", "/healthz", "/", "/docs", "/openapi.json", "/redoc", "/governance", "/metrics",
)


def resolve_key(service: str = "", explicit: str | None = None) -> str | None:
    """Per-service key wins over the shared one. Returns None when auth is disabled."""
    if explicit:
        return explicit
    if service:
        slug = service.upper().replace("-", "_").replace(" ", "_")
        per = os.environ.get(f"HCLS_API_KEY_{slug}")
        if per:
            return per
    return os.environ.get("HCLS_API_KEY") or None


def install_api_key_auth(app, *, service: str = "", api_key: str | None = None,
                         open_paths: Iterable[str] = DEFAULT_OPEN_PATHS):
    """Attach the gate. No-op when no key is configured, so existing deployments do not break.

    Returns True if the gate was installed, False if auth is disabled.
    """
    key = resolve_key(service, api_key)
    if not key:
        return False

    from fastapi import Request
    from fastapi.responses import JSONResponse

    allowed = tuple(open_paths)

    @app.middleware("http")
    async def _require_api_key(request: "Request", call_next):
        if request.method == "OPTIONS" or request.url.path in allowed:
            return await call_next(request)
        supplied = request.headers.get("X-API-Key") or ""
        # constant-time: a plain == leaks key length and prefix through timing
        if not hmac.compare_digest(supplied, key):
            return JSONResponse(
                {"detail": "missing or invalid X-API-Key",
                 "service": service or "hcls-ai-factory"},
                status_code=401,
            )
        return await call_next(request)

    return True


def auth_status(service: str = "") -> dict:
    """For /health and /governance — reports whether the gate is on WITHOUT leaking the key."""
    key = resolve_key(service)
    return {
        "api_key_required": bool(key),
        "posture": "fail-closed" if key else "open (trusted network / demo)",
        "header": "X-API-Key" if key else None,
    }
