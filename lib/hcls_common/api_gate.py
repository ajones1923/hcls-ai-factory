"""Shared governance surface for factory services — the input-validation and
output-honesty gates, wired into any FastAPI service in one line.

The gates used to fire only on the MCP/composer path, so services that talked
to their own FastAPI routes bypassed governance entirely. This module makes the
gates a drop-in every service (and every new engine/agent) inherits.

Retrofit an existing app:

    from hcls_common.api_gate import install_governance
    install_governance(app, service="cart", capability_id="cart-intelligence-agent")

or build a pre-governed app (preferred for new services):

    from hcls_common.api_gate import create_governed_app
    app = create_governed_app("cart", capability_id="cart-intelligence-agent")

Then, in the endpoints:

    from hcls_common.api_gate import require_valid_input, honesty_flags
    payload = require_valid_input("cart-intelligence-agent", payload)   # 422 on bad input
    flags = honesty_flags(answer_text)                                  # overclaim scan (no LLM)

FastAPI is imported lazily so `import hcls_common` never requires web deps.
"""
from __future__ import annotations

import time
import uuid
from typing import Any


def create_governed_app(service: str, *, capability_id: str | None = None, **fastapi_kwargs):
    """Return a FastAPI app with the governance middleware + /governance already installed."""
    from fastapi import FastAPI

    app = FastAPI(**fastapi_kwargs)
    install_governance(app, service=service, capability_id=capability_id)
    return app


def install_governance(app, *, service: str = "", capability_id: str | None = None):
    """Attach the governance middleware (request id + timing + governed header) and a
    /governance info endpoint to an existing FastAPI app. Idempotent-ish; safe to call once."""
    from fastapi import Request

    @app.middleware("http")
    async def _governance(request: "Request", call_next):
        rid = request.headers.get("x-request-id") or uuid.uuid4().hex[:16]
        t0 = time.monotonic()
        response = await call_next(request)
        response.headers["X-Request-ID"] = rid
        response.headers["X-HCLS-Governed"] = service or "1"
        response.headers["X-HCLS-Duration-ms"] = f"{(time.monotonic() - t0) * 1000:.1f}"
        return response

    @app.get("/governance", tags=["status"])
    def _governance_info():
        return {
            "service": service,
            "capability_id": capability_id,
            "gates": ["input-validation", "output-honesty"],
            "how": {
                "input": "call require_valid_input(capability_id, payload) in POST handlers",
                "output": "call honesty_flags(text) (deterministic) or assert_publishable(text, llm=...)",
            },
        }

    return app


def require_valid_input(capability_id: str, payload: dict[str, Any] | None) -> dict[str, Any]:
    """Validate a request payload against a capability's input contract.

    Applies defaults, clamps out-of-range numerics (logged as WARN), and rejects
    ERROR issues (missing-required / enum-violation) with HTTP 422. Unknown
    capability ids pass through unchanged (never hard-fail a running service)."""
    from fastapi import HTTPException

    from hcls_common.capability_registry import get_registry, validate_inputs

    try:
        cap = get_registry().get(capability_id)
    except Exception:
        return dict(payload or {})
    cleaned, issues = validate_inputs(cap, payload)
    errors = [i for i in issues if i.startswith("ERROR")]
    if errors:
        raise HTTPException(status_code=422, detail={"input_errors": errors})
    return cleaned


def honesty_flags(text: str) -> list[dict]:
    """Deterministic overclaim / missing-disclaimer scan (never calls an LLM)."""
    from hcls_common.verify_gate import honesty_check

    return honesty_check(text)


def assert_publishable(text: str, *, evidence: Any = None, llm: Any = None) -> dict:
    """Full verify + honesty gate; raise HTTP 422 if the text is not publishable.

    Use for send-ready clinical text. With no `llm`, runs the deterministic layer only."""
    from fastapi import HTTPException

    from hcls_common.verify_gate import is_publishable, verify_text

    verdict = verify_text(text, evidence=evidence, llm=llm)
    if not is_publishable(verdict):
        raise HTTPException(status_code=422, detail={"verify": verdict})
    return verdict
