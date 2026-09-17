# Copyright 2026 Adam Jones
# SPDX-License-Identifier: Apache-2.0
# Part of the HCLS AI Factory — https://github.com/ajones1923/hcls-ai-factory
# Licensed under the Apache License, Version 2.0. See LICENSE and NOTICE.
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

import contextvars
import time
import uuid
from typing import Any

# Which gates actually executed on THIS request. The middleware cannot know -- the gates are
# called from inside handlers -- so they record here and the response header reports the truth.
# Holds a MUTABLE set. The middleware must not rebind this var from inside a handler:
# Starlette runs sync endpoints in a threadpool with a copied context, so a `.set()` made in
# the handler is invisible to the middleware that awaited it. Mutating one shared set works
# because both see the same object. A test asserts the header actually appears.
_GATES_RUN: contextvars.ContextVar[set] = contextvars.ContextVar("_hcls_gates_run", default=None)


def _mark_gate(name: str) -> None:
    bucket = _GATES_RUN.get()
    if bucket is not None:
        bucket.add(name)


def create_governed_app(service: str, *, capability_id: str | None = None, **fastapi_kwargs):
    """Return a FastAPI app with the governance middleware + /governance already installed."""
    from fastapi import FastAPI

    app = FastAPI(**fastapi_kwargs)
    install_governance(app, service=service, capability_id=capability_id)
    return app


def _auth_status(service: str) -> dict:
    try:
        from hcls_common.api_auth import auth_status
        return auth_status(service)
    except Exception:
        return {"api_key_required": None}


def install_governance(app, *, service: str = "", capability_id: str | None = None):
    """Attach the governance middleware (request id + timing + governed header) and a
    /governance info endpoint to an existing FastAPI app. Idempotent-ish; safe to call once."""
    from fastapi import Request

    @app.middleware("http")
    async def _governance(request: "Request", call_next):
        rid = request.headers.get("x-request-id") or uuid.uuid4().hex[:16]
        gates: set[str] = set()
        _GATES_RUN.set(gates)
        t0 = time.monotonic()
        response = await call_next(request)
        response.headers["X-Request-ID"] = rid
        response.headers["X-HCLS-Service"] = service or "hcls-ai-factory"
        # AUDIT FIX (2026-08-16): this used to emit `X-HCLS-Governed: <service>` on EVERY request,
        # including requests where no gate ran -- this middleware only adds a request id and timing.
        # A header asserting governance on an ungoverned request is worse than no header on a
        # project whose credibility rests on honesty by construction. Now it reports only the gates
        # that actually executed, and is omitted entirely when none did.
        ran = sorted(gates)
        if ran:
            response.headers["X-HCLS-Governed"] = ",".join(ran)
        response.headers["X-HCLS-Duration-ms"] = f"{(time.monotonic() - t0) * 1000:.1f}"
        return response

    # The output gate is no longer opt-in: see install_output_honesty for why.
    install_output_honesty(app, service=service)

    @app.get("/governance", tags=["status"])
    def _governance_info():
        return {
            "service": service,
            "capability_id": capability_id,
            "gates_available": ["input-validation", "output-honesty"],
            "output_gate": ("AUTOMATIC since 2026-09-15: every JSON answer over 200 chars is "
                            "scanned by the deterministic honesty register and carries the "
                            "decision-support disclaimer. Block-severity overclaims are WITHHELD "
                            "by default; set HCLS_HONESTY_ENFORCE=0 to annotate instead."),
            "input_gate_is_opt_in": ("a handler must still call require_valid_input(); "
                                     "X-HCLS-Governed lists only what actually ran."),
            "auth": _auth_status(service),
            "how": {
                "input": "call require_valid_input(capability_id, payload) in POST handlers",
                "output": "call honesty_flags(text) (deterministic) or assert_publishable(text, llm=...)",
            },
        }

    return app


# ── output-honesty response gate ─────────────────────────────────────────────
DISCLAIMER = (
    "\n\n---\n*Decision support for a qualified clinician — not a diagnosis, not a "
    "prescription, and not a substitute for clinical judgement. Generated from retrieved "
    "evidence; verify every claim against the cited primary source before acting on it.*"
)

# Response keys that carry prose a clinician might act on.
_ANSWER_KEYS = ("answer", "response", "summary", "interpretation", "narrative", "brief")



# Keys a service uses to return the passages an answer was built from.
_EVIDENCE_KEYS = ("citations", "sources", "evidence", "results", "matches")

_VERIFY_LLM: Any = None
_VERIFY_LLM_TRIED = False


def _verify_llm():
    """Lazy singleton LLM client for adversarial claim checking (None if no key).

    Built once per process: the middleware runs per request and constructing a client each
    time would add a connection setup to every answer.
    """
    global _VERIFY_LLM, _VERIFY_LLM_TRIED
    if _VERIFY_LLM_TRIED:
        return _VERIFY_LLM
    _VERIFY_LLM_TRIED = True
    import os
    if not (os.getenv("ANTHROPIC_API_KEY") or "").strip():
        return None
    try:
        from hcls_common.llm_client import AnthropicClient
        _VERIFY_LLM = AnthropicClient()
    except Exception:
        _VERIFY_LLM = None
    return _VERIFY_LLM


def _evidence_from(payload: dict) -> list:
    for k in _EVIDENCE_KEYS:
        v = payload.get(k)
        if isinstance(v, list) and v:
            return v
    return []


def install_output_honesty(app, *, service: str = ""):
    """Run the deterministic honesty gate over every JSON answer this service returns.

    Why this is middleware and not a line in each handler: the gates were opt-in, and the
    2026-08-15 audit found them called by 1 of 12 entrypoints. That was survivable while the
    agents returned retrieved passages. It stopped being survivable the moment an API key was
    present, because the same endpoints now return generated clinical prose -- dosing, response
    rates, management guidance -- and a 6,400-character answer was going out with no statement
    that it is decision support. On a project whose thesis is honesty by construction, that is
    the one failure that cannot be allowed to depend on twelve handlers remembering.

    Appends the decision-support disclaimer when clinical prose lacks one, attaches any
    overclaim findings to the payload under `honesty`, and reports them in X-HCLS-Honesty.
    Content is annotated, never silently rewritten -- except that a `block`-severity finding
    (a cure claim, diagnostic certainty, a clearance claim about THIS platform) replaces the
    generated prose, because those must not ship at all. Enforcing by default since
    2026-09-15; HCLS_HONESTY_ENFORCE=0 annotates instead. Retrieved evidence is never
    withheld -- only the generated text.
    """
    import json as _json
    import os as _os

    from fastapi import Request
    from starlette.responses import Response

    @app.middleware("http")
    async def _output_honesty(request: "Request", call_next):
        response = await call_next(request)
        ctype = response.headers.get("content-type", "")
        if "application/json" not in ctype or response.status_code >= 400:
            return response
        body = b""
        async for chunk in response.body_iterator:
            body += chunk
        try:
            payload = _json.loads(body)
        except Exception:
            return Response(content=body, status_code=response.status_code,
                            headers=dict(response.headers), media_type=response.media_type)

        if isinstance(payload, dict):
            key = next((k for k in _ANSWER_KEYS
                        if isinstance(payload.get(k), str) and len(payload[k]) > 200), None)
            if key:
                text = payload[key]
                try:
                    findings = honesty_flags(text)          # marks the gate as having run
                except Exception:
                    findings = []
                # ── adversarial layer ────────────────────────────────────────
                # The deterministic register above catches PHRASING ("cures", "zero risk").
                # It cannot catch a fabricated statistic or an invented trial number, which is
                # the actual failure mode of generated clinical prose. This extracts atomic
                # claims and asks a model to REFUTE each against the passages the answer was
                # built from -- refutation, not confirmation, is what catches invented numbers.
                # Off when there is no evidence to check against, or no key, or
                # HCLS_VERIFY_CLAIMS=0. It must never break generation.
                refuted: list = []
                if _os.getenv("HCLS_VERIFY_CLAIMS", "1") != "0":
                    evidence = _evidence_from(payload)
                    llm = _verify_llm() if evidence else None
                    if llm is not None:
                        try:
                            from hcls_common.verify_gate import verify_claims
                            rep = verify_claims(text, evidence, llm)
                            _mark_gate("claim-verification")
                            payload["verification"] = {
                                "claims": len(rep["claims"]),
                                "supported": rep["supported"],
                                "refuted": rep["refuted"],
                                "unsupported": rep["unsupported"],
                                "flagged": rep["flagged"][:8],
                            }
                            refuted = [x for x in rep["flagged"] if x.get("verdict") == "refuted"]
                        except Exception as exc:          # never break the answer
                            payload["verification"] = {"error": f"{type(exc).__name__}"}

                blocking = [f for f in findings if f.get("severity") == "block"]
                # ENFORCING BY DEFAULT since 2026-09-15 (Adam's decision). Set
                # HCLS_HONESTY_ENFORCE=0 to annotate instead of withhold.
                #
                # Safe to default on only because the regulatory/efficacy rules became
                # subject-aware first: "this platform is FDA-approved" blocks, while
                # "tisagenlecleucel is FDA-approved" degrades to a warning. Without that,
                # enforcement would have withheld correct clinical answers about approved
                # drugs -- worse than the problem it solves. Rules that are unsafe whatever
                # the subject (a cure claim, absolute certainty, zero risk) always block.
                if (blocking or refuted) and _os.getenv("HCLS_HONESTY_ENFORCE", "1") != "0":
                    reasons = [f.get("message", "") for f in blocking]
                    reasons += [f"Claim contradicted by the cited evidence: "
                                f"{r.get('claim', '')[:140]} — {r.get('why', '')[:140]}"
                                for r in refuted]
                    payload[key] = (
                        "**Withheld by the output-honesty gate.**\n\n"
                        + "\n".join(f"- {m}" for m in reasons)
                        + "\n\nThe generated text made a claim this platform must not publish. "
                          "The retrieved evidence is unchanged and still returned; only the "
                          "generated prose was withheld. Set HCLS_HONESTY_ENFORCE=0 to receive "
                          "it annotated instead." + DISCLAIMER)
                    payload["withheld"] = True
                elif not _has_disclaimer(text):
                    payload[key] = text + DISCLAIMER
                if findings:
                    payload["honesty"] = findings
                body = _json.dumps(payload).encode()

        headers = dict(response.headers)
        headers.pop("content-length", None)
        out = Response(content=body, status_code=response.status_code,
                       headers=headers, media_type="application/json")
        ran = sorted(_GATES_RUN.get(set()) or set())
        if ran:
            out.headers["X-HCLS-Governed"] = ",".join(ran)
        return out

    return app


def _has_disclaimer(text: str) -> bool:
    low = text.lower()
    return any(p in low for p in ("decision support", "not a diagnosis", "qualified clinician",
                                  "research use", "not a substitute for clinical"))


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
    _mark_gate("input-validation")
    cleaned, issues = validate_inputs(cap, payload)
    errors = [i for i in issues if i.startswith("ERROR")]
    if errors:
        raise HTTPException(status_code=422, detail={"input_errors": errors})
    return cleaned


def honesty_flags(text: str) -> list[dict]:
    """Deterministic overclaim / missing-disclaimer scan (never calls an LLM)."""
    from hcls_common.verify_gate import honesty_check

    _mark_gate("output-honesty")
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
