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
# AUDIT FIX (2026-11-03): the first six were the whole list, so an answer returned under
# `recommendation` or `rationale` -- which several routes do -- was never scanned.
_ANSWER_KEYS = (
    "answer", "response", "summary", "interpretation", "narrative", "brief",
    "recommendation", "rationale", "assessment", "explanation", "impression",
    "conclusion", "advice", "guidance", "report", "discussion", "findings",
    "answer_text", "clinical_interpretation", "plan",
)

# Keys a service uses to return the passages an answer was built from.
_EVIDENCE_KEYS = ("citations", "sources", "evidence", "results", "matches")

# Keys this gate writes into the payload itself. Never scan or rewrite our own annotations.
_GATE_KEYS = ("honesty", "verification", "withheld")

# The walk is bounded so a pathological payload cannot turn one response into a long scan.
_MAX_WALK_DEPTH = 6
_MAX_WALK_NODES = 2000


def _answer_slots(payload: Any) -> list[tuple]:
    """Every answer-key string in the payload, at any depth.

    Returns `(container, key, path, text)` so the caller can rewrite in place.

    AUDIT FIX (2026-11-03): this used to be `next((k for k in _ANSWER_KEYS ...))` over the
    top-level dict only, which left three ways past the gate: prose nested under another key
    (`{"result": {"interpretation": ...}}`), a second answer key in the same payload once the
    first one matched, and a list payload, which was not looked at at all.

    Subtrees under `_EVIDENCE_KEYS` are skipped. Retrieved evidence is quoted source material;
    the gate annotates generated prose and must never rewrite a passage a clinician will check
    the answer against.
    """
    out: list[tuple] = []
    budget = [0]

    def walk(node, path: str, depth: int) -> None:
        if depth > _MAX_WALK_DEPTH or budget[0] > _MAX_WALK_NODES:
            return
        if isinstance(node, dict):
            for k, v in node.items():
                budget[0] += 1
                if k in _EVIDENCE_KEYS or k in _GATE_KEYS:
                    continue
                p = f"{path}.{k}" if path else k
                if isinstance(v, str):
                    if k in _ANSWER_KEYS:
                        out.append((node, k, p, v))
                elif isinstance(v, (dict, list)):
                    walk(v, p, depth + 1)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                budget[0] += 1
                if isinstance(v, (dict, list)):
                    walk(v, f"{path}[{i}]", depth + 1)

    walk(payload, "", 0)
    return out


def _needs_disclaimer(text: str, findings: list) -> bool:
    """Whether to append the decision-support footer.

    Two triggers, deliberately: the register's own missing-disclaimer finding (which keys on
    clinical trigger words), and the pre-existing length rule, kept so that nothing which used
    to be disclaimed stops being disclaimed.
    """
    if _has_disclaimer(text):
        return False
    if len(text) > 200:
        return True
    return any("disclaimer" in (f.get("message") or "").lower() for f in findings)

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

    Scope, after the 2026-11-03 audit: every string under an `_ANSWER_KEYS` key at any depth
    in a JSON body, of any length, in a dict or a list payload. Subtrees under
    `_EVIDENCE_KEYS` are skipped so quoted source material is never rewritten.

    What it still does NOT cover, so nobody has to re-derive it: a streamed body
    (`text/event-stream`, `application/x-ndjson`) is passed through unscanned and marked
    `X-HCLS-Honesty: not-scanned`; prose returned under a key outside `_ANSWER_KEYS` is
    invisible to it; and a module that calls its own generator in-process, without crossing
    the HTTP boundary, bypasses middleware entirely -- those handlers must still call
    `honesty_flags()` or `assert_publishable()` themselves.
    """
    import json as _json
    import os as _os

    from fastapi import Request
    from starlette.responses import Response

    @app.middleware("http")
    async def _output_honesty(request: "Request", call_next):
        response = await call_next(request)
        ctype = response.headers.get("content-type", "")
        # A streamed answer cannot be scanned without buffering it, and buffering it would
        # defeat the point of streaming. The gate does not run, and says so rather than
        # letting the response look like every other gated one.
        if any(s in ctype for s in ("text/event-stream", "application/x-ndjson")):
            response.headers["X-HCLS-Honesty"] = (
                "not-scanned; streamed response (the gate cannot buffer a stream)")
            return response
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

        slots = _answer_slots(payload)
        # At most one claim-verification call per response: it costs an LLM round trip, so it
        # runs against the longest answer, which is the generated prose in every payload shape
        # seen in this repo. The deterministic register runs over all of them.
        primary = max(range(len(slots)), key=lambda i: len(slots[i][3])) if slots else -1
        all_findings: list[dict] = []
        changed = False
        withheld = 0
        verified = False

        for idx, (container, key, path, text) in enumerate(slots):
            try:
                findings = honesty_flags(text)              # marks the gate as having run
            except Exception:
                findings = []
            # ── adversarial layer ────────────────────────────────────────────
            # The deterministic register above catches PHRASING ("cures", "zero risk").
            # It cannot catch a fabricated statistic or an invented trial number, which is
            # the actual failure mode of generated clinical prose. This extracts atomic
            # claims and asks a model to REFUTE each against the passages the answer was
            # built from -- refutation, not confirmation, is what catches invented numbers.
            # Off when there is no evidence to check against, or no answer, or
            # HCLS_VERIFY_CLAIMS=0. It must never break generation.
            refuted: list = []
            if idx == primary and _os.getenv("HCLS_VERIFY_CLAIMS", "1") != "0":
                evidence = _evidence_from(payload) if isinstance(payload, dict) else []
                if not evidence and isinstance(container, dict):
                    evidence = _evidence_from(container)
                llm = _verify_llm() if evidence else None
                if llm is not None:
                    report: dict = {}
                    try:
                        from hcls_common.verify_gate import verify_claims
                        rep = verify_claims(text, evidence, llm)
                        _mark_gate("claim-verification")
                        verified = True
                        report = {
                            "claims": len(rep["claims"]),
                            "supported": rep["supported"],
                            "refuted": rep["refuted"],
                            "unsupported": rep["unsupported"],
                            "flagged": rep["flagged"][:8],
                        }
                        refuted = [x for x in rep["flagged"] if x.get("verdict") == "refuted"]
                    except Exception as exc:              # never break the answer
                        report = {"error": f"{type(exc).__name__}"}
                    if isinstance(payload, dict):
                        payload["verification"] = report

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
            #
            # Length is NOT a condition here. It used to be -- only strings over 200
            # characters were looked at -- and "This platform cures the patient's disease."
            # is 41 characters.
            if (blocking or refuted) and _os.getenv("HCLS_HONESTY_ENFORCE", "1") != "0":
                reasons = [f.get("message", "") for f in blocking]
                reasons += [f"Claim contradicted by the cited evidence: "
                            f"{r.get('claim', '')[:140]} — {r.get('why', '')[:140]}"
                            for r in refuted]
                container[key] = withheld_message(reasons, evidence_returned=True)
                withheld += 1
                changed = True
            elif _needs_disclaimer(text, findings):
                container[key] = text + DISCLAIMER
                changed = True
            for f in findings:
                all_findings.append({**f, "path": path})

        if all_findings and isinstance(payload, dict):
            payload["honesty"] = all_findings
            changed = True
        if withheld and isinstance(payload, dict):
            payload["withheld"] = True
            changed = True
        if changed:
            body = _json.dumps(payload).encode()

        headers = dict(response.headers)
        headers.pop("content-length", None)
        out = Response(content=body, status_code=response.status_code,
                       headers=headers, media_type="application/json")
        # AUDIT FIX (2026-11-03): this read the gates out of the contextvar, which the inner
        # middleware sets in a different context -- so the automatic output gate, the one gate
        # that always runs, never appeared in X-HCLS-Governed. This middleware knows what it
        # ran, so it unions its own gates with whatever the inner middleware reported.
        ran = {g for g in (headers.get("x-hcls-governed") or "").split(",") if g}
        if slots:
            ran.add("output-honesty")
        if verified:
            ran.add("claim-verification")
        if ran:
            out.headers["X-HCLS-Governed"] = ",".join(sorted(ran))
        # The docstring promised this header since 2026-09-15 and nothing ever set it.
        if slots:
            sev = [f.get("severity") for f in all_findings]
            out.headers["X-HCLS-Honesty"] = "; ".join(
                [f"scanned={len(slots)}"]
                + ([f"block={sev.count('block')}"] if "block" in sev else [])
                + ([f"warn={sev.count('warn')}"] if "warn" in sev else [])
                + ([f"withheld={withheld}"] if withheld else []))
        return out

    return app


def _has_disclaimer(text: str) -> bool:
    low = text.lower()
    return any(p in low for p in ("decision support", "not a diagnosis", "qualified clinician",
                                  "research use", "not a substitute for clinical"))


def withheld_message(reasons: list[str], *, evidence_returned: bool) -> str:
    """The text that replaces withheld prose.

    One builder for both callers -- the response middleware and `gate_text` -- so the wording
    cannot drift between the HTTP surface and the in-process UI surface. `evidence_returned`
    adds the sentence about retrieved passages, which is only true where a payload carries
    them alongside the answer.
    """
    body = "\n".join(f"- {m}" for m in reasons)
    tail = ("The generated text made a claim this platform must not publish. ")
    if evidence_returned:
        tail += ("The retrieved evidence is unchanged and still returned; only the generated "
                 "prose was withheld. ")
    tail += "Set HCLS_HONESTY_ENFORCE=0 to receive it annotated instead."
    return f"**Withheld by the output-honesty gate.**\n\n{body}\n\n{tail}{DISCLAIMER}"


def gate_text(text: str) -> tuple[str, list[dict]]:
    """Apply the output-honesty gate to one piece of generated prose, outside HTTP.

    Returns `(text_to_show, findings)`, with the same decisions the response middleware makes:
    a `block`-severity finding replaces the prose, anything else gets the decision-support
    disclaimer appended when it is missing.

    This exists because `install_output_honesty` is middleware, and a Streamlit UI that imports
    its agent's RAG engine and calls it in-process never crosses the HTTP boundary -- eight of
    the fourteen UIs in this repo work that way, so the gate that is "automatic since
    2026-09-15" never ran on the surface a clinician actually reads. Those UIs call this.

    `HCLS_HONESTY_ENFORCE=0` annotates instead of withholding, as it does in the middleware.
    """
    import os as _os

    if not isinstance(text, str) or not text:
        return text, []
    try:
        findings = honesty_flags(text)
    except Exception:
        return text, []

    blocking = [f for f in findings if f.get("severity") == "block"]
    if blocking and _os.getenv("HCLS_HONESTY_ENFORCE", "1") != "0":
        reasons = [f.get("message", "") for f in blocking]
        return withheld_message(reasons, evidence_returned=False), findings
    if _needs_disclaimer(text, findings):
        return text + DISCLAIMER, findings
    return text, findings


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
