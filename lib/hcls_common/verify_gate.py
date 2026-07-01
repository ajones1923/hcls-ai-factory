"""
Verify-gate / honesty register for the HCLS AI Factory (ORCH-11, 2026-06-29).

Ports the marketing/Paai adversarial verify-gate to HCLS as a reusable, non-bypassable
check on generated artifacts (RAG answers, target hypotheses, drug-discovery reports,
agent outputs). Two layers:

  1. honesty_check(text)  — DETERMINISTIC clinical/genomics overclaim register. No LLM,
     always runs. Catches the failure modes that actually matter on a clinical-adjacent
     platform: claiming regulatory clearance, diagnostic certainty, cures/guarantees, or
     issuing treatment directives — and flags missing research-use disclaimers.

  2. verify_claims(text, evidence, llm) — ADVERSARIAL grounding. A fast model extracts
     atomic claims; a high-stakes model is instructed to *refute* each against the cited
     evidence (refutation, not confirmation, is what catches hallucinated stats/quotes).
     Optional — degrades to honesty-only when no llm is supplied.

verify_text() combines both into one verdict; is_publishable() is the gate. The honesty
layer is the load-bearing, always-available piece (deterministic = predictable, the HCLS
honesty ethos); the LLM layer is the enhancement when a client is available.

The `llm` argument is any object exposing
    generate_json(prompt: str, system_prompt: str = ..., max_tokens: int = ...) -> dict
(e.g. an hcls_common LLM client) — keeping this module decoupled from any one provider.
"""
from __future__ import annotations

import re
from typing import Any, Optional

# ── Layer 1: deterministic clinical-honesty register ────────────────────────
# (pattern, severity, message). severity: "block" = must not publish; "warn" = review.
_OVERCLAIM_RULES: list[tuple[str, str, str]] = [
    (r"\bFDA[- ]?(approved|cleared)\b", "block", "Regulatory-clearance claim (FDA) — the platform is research-use, not FDA-cleared."),
    (r"\b(CE[- ]?marked|510\(k\) cleared)\b", "block", "Regulatory-clearance claim — not cleared."),
    (r"\b(definitive|confirmed|certain|conclusive)\s+(diagnos[ie]s|diagnosis)\b", "block", "Diagnostic-certainty overclaim — outputs are decision support, not a diagnosis."),
    (r"\b(confirms?|establishes?)\s+(the\s+)?diagnosis\b", "block", "Diagnostic-certainty overclaim."),
    (r"\b(cure[sd]?|will cure|proven to cure|eradicat\w+)\b", "block", "Cure claim — unsupported and unsafe to assert."),
    (r"\b(100%|guarantee[sd]?|guaranteed|certain to|always works|never fails)\b", "block", "Absolute-certainty overclaim."),
    (r"\b(you should (take|stop|start)|discontinue your|prescribe[sd]?)\b", "warn", "Treatment-directive language — outputs inform clinicians, they don't direct patients."),
    (r"\b(safe for all|no side effects|zero risk|completely safe)\b", "block", "Safety overclaim."),
]
_OVERCLAIM = [(re.compile(p, re.I), sev, msg) for (p, sev, msg) in _OVERCLAIM_RULES]

# A clinical-ish claim with none of these caveats present → missing-disclaimer warning.
_DISCLAIMER_HINTS = re.compile(
    r"research[- ]use|decision support|not (a )?(diagnos|clinical|substitute)|consult|"
    r"clinician review|investigational|for informational|synthetic", re.I)
_CLINICAL_TRIGGER = re.compile(
    r"\b(diagnos\w+|treat\w+|therap\w+|patient|prognos\w+|pathogenic|variant|dose|dosing)\b", re.I)


def honesty_check(text: str) -> list[dict]:
    """Deterministic scan for overclaims + missing disclaimers. Returns a list of
    violations: {severity, message, excerpt}. Never calls an LLM."""
    text = text or ""
    out: list[dict] = []
    for rx, sev, msg in _OVERCLAIM:
        m = rx.search(text)
        if m:
            s = max(0, m.start() - 30)
            out.append({"severity": sev, "message": msg,
                        "excerpt": text[s:m.end() + 30].strip()})
    if _CLINICAL_TRIGGER.search(text) and not _DISCLAIMER_HINTS.search(text):
        out.append({"severity": "warn",
                    "message": "Clinical content without a research-use / decision-support / clinician-review disclaimer.",
                    "excerpt": ""})
    return out


# ── Layer 2: adversarial LLM claim verification (optional) ───────────────────
_EXTRACT_SYS = ("You extract atomic, checkable factual claims from a generated answer. "
                "Return ONLY JSON: {\"claims\": [\"...\", ...]}. Max 20 claims; skip hedged/opinion text.")
_REFUTE_SYS = ("You are an adversarial fact-checker. For each claim, try to REFUTE it using "
               "ONLY the provided evidence. Verdict 'supported' only if the evidence directly "
               "supports it; 'refuted' if evidence contradicts it; 'unsupported' if evidence is "
               "silent. Default to 'unsupported' when unsure. Return ONLY JSON: "
               "{\"verdicts\": [{\"claim_index\": int, \"verdict\": \"supported|refuted|unsupported\", \"why\": \"...\"}]}.")


def extract_claims(text: str, llm: Any, max_claims: int = 20) -> list[str]:
    try:
        parsed = llm.generate_json(f"ANSWER:\n{text}", system_prompt=_EXTRACT_SYS, max_tokens=1024)
        claims = [c for c in (parsed or {}).get("claims", []) if isinstance(c, str) and c.strip()]
        return claims[:max_claims]
    except Exception:
        return []


def _evidence_block(evidence: Optional[list]) -> str:
    if not evidence:
        return "(no evidence provided)"
    lines = []
    for i, e in enumerate(evidence, 1):
        txt = e.get("text", e) if isinstance(e, dict) else e
        lines.append(f"[{i}] {str(txt)[:600]}")
    return "\n".join(lines)


def verify_claims(text: str, evidence: Optional[list], llm: Any) -> dict:
    """Extract claims and adversarially check each against evidence. Returns
    {claims, verdicts, supported, refuted, unsupported, flagged}."""
    claims = extract_claims(text, llm)
    if not claims:
        return {"claims": [], "verdicts": [], "supported": 0, "refuted": 0,
                "unsupported": 0, "flagged": []}
    numbered = "\n".join(f"{i+1}. {c}" for i, c in enumerate(claims))
    prompt = (f"EVIDENCE:\n{_evidence_block(evidence)}\n\nCLAIMS:\n{numbered}")
    try:
        parsed = llm.generate_json(prompt, system_prompt=_REFUTE_SYS, max_tokens=2048) or {}
        verdicts = parsed.get("verdicts", [])
    except Exception:
        verdicts = []
    sup = ref = uns = 0
    flagged = []
    for v in verdicts:
        verdict = (v.get("verdict") or "unsupported").lower()
        idx = int(v.get("claim_index", 0)) - 1
        claim = claims[idx] if 0 <= idx < len(claims) else ""
        if verdict == "supported":
            sup += 1
        elif verdict == "refuted":
            ref += 1
            flagged.append({"claim": claim, "verdict": "refuted", "why": v.get("why", "")})
        else:
            uns += 1
            flagged.append({"claim": claim, "verdict": "unsupported", "why": v.get("why", "")})
    return {"claims": claims, "verdicts": verdicts, "supported": sup,
            "refuted": ref, "unsupported": uns, "flagged": flagged}


# ── Combined gate ────────────────────────────────────────────────────────────
def verify_text(text: str, evidence: Optional[list] = None, llm: Any = None,
                *, unsupported_tolerance: int = 0) -> dict:
    """Run the honesty register (always) + adversarial verification (if llm given).
    Returns a verdict with `status` in {verified, draft, blocked}."""
    violations = honesty_check(text)
    blocked = any(v["severity"] == "block" for v in violations)
    warns = sum(1 for v in violations if v["severity"] == "warn")

    claim_report = verify_claims(text, evidence, llm) if llm is not None else None
    refuted = claim_report["refuted"] if claim_report else 0
    unsupported = claim_report["unsupported"] if claim_report else 0

    if blocked or refuted > 0:
        status = "blocked"
    elif warns > 0 or unsupported > unsupported_tolerance:
        status = "draft"
    else:
        status = "verified"

    return {
        "status": status,
        "honesty_violations": violations,
        "claims": claim_report["claims"] if claim_report else [],
        "flagged": claim_report["flagged"] if claim_report else [],
        "summary": {
            "blocked": blocked, "warnings": warns,
            "supported": claim_report["supported"] if claim_report else 0,
            "refuted": refuted, "unsupported": unsupported,
            "llm_checked": claim_report is not None,
        },
    }


def is_publishable(verdict: dict) -> bool:
    """The gate: only a clean 'verified' verdict may publish."""
    return verdict.get("status") == "verified"
