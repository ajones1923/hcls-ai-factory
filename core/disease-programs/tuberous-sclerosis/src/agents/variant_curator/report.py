"""
AI-labeled DRAFT molecular-genetics report (PRD §3 FR-VC-5; master paper §8).

Renders the Variant Curator's interpretation into a clinician-readable draft report and
writes it under tsc-reports/. The report is explicitly a DRAFT pending board-certified
molecular-geneticist sign-off — never an autonomous result. A lightweight sign-off
registry records the human gate (signer + decision) without adding an event type.
"""
from __future__ import annotations

import json
import threading
from pathlib import Path

from config.settings import settings

_REPORTS = Path(settings.DATA_DIR) / "tsc-reports"
_SIGNOFFS = Path(settings.STATE_DIR) / "signoffs.json"


def render_draft(patient_id: str, variant_interp: dict) -> str:
    """Markdown draft report from the curator's projection (variant_interp dict)."""
    p = (variant_interp or {}).get("primary") or {}
    lines = [
        f"# DRAFT Molecular Genetics Report — {patient_id}",
        "",
        f"> **{settings.WATERMARK}**  ",
        "> **STATUS: DRAFT — AI-generated decision support. Pending board-certified "
        "molecular-geneticist review and sign-off. Not a final clinical result.**",
        "",
        "## Result",
    ]
    if not p:
        lines += [
            "No reportable variant on the available sample (e.g., negative blood test).",
            "",
            (variant_interp or {}).get("note", ""),
        ]
    else:
        lines += [
            f"- **Gene:** {p.get('gene')}",
            f"- **Variant:** {p.get('hgvsc') or ''} {p.get('hgvsp') or ''}".rstrip(),
            f"- **Consequence:** {p.get('consequence')}",
            f"- **Classification (ACMG-AMP):** **{p.get('acmg_classification')}** "
            f"({p.get('acmg_rule')})",
            f"- **Zygosity / VAF:** {'mosaic' if p.get('mosaic') else 'germline'} · "
            f"VAF {p.get('vaf')}" + ("  · **recovered from affected tissue**" if p.get('recovered') else ""),
            f"- **Read-level QC:** depth {p.get('depth')}, alt reads {p.get('alt_reads')}, "
            f"strand balance {p.get('strand_balance')} → {p.get('artifact_assessment')}",
            f"- **Orthogonal validation:** {'ddPCR recommended' if p.get('ddpcr_recommended') else 'not indicated'}",
            f"- **Candidate artifacts rejected (strand-bias filter):** {variant_interp.get('artifacts_filtered', 0)}",
            "",
            "## ACMG-AMP criteria",
        ]
        for c in p.get("acmg_criteria", []):
            lines.append(f"- **{c.get('code')}** ({c.get('bucket')}): {c.get('rationale')}")
        nar = (variant_interp or {}).get("synthesis_narrative") or {}
        if nar.get("summary"):
            lines += ["", "## Interpretive narrative (AI-generated)", ""]
            for pc in nar.get("per_criterion", []):
                lines.append(f"- **{pc.get('code')}** — {pc.get('reasoning')}")
            lines += ["", f"_{nar['summary']}_"]
    lines += [
        "",
        "## Sign-off",
        "- [ ] Reviewed by board-certified molecular geneticist",
        "- Reviewer: __________________________   Date: __________",
        "",
        "_This draft is decision support and requires clinician review. SYNTHETIC demonstration data._",
    ]
    return "\n".join(lines)


def write_draft(patient_id: str, variant_interp: dict) -> Path:
    out = _REPORTS / patient_id / "variant_draft.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(render_draft(patient_id, variant_interp), encoding="utf-8")
    return out


# ── sign-off registry (the human gate) ──────────────────────────────────────
_lock = threading.Lock()


def _load() -> dict:
    try:
        return json.loads(_SIGNOFFS.read_text())
    except Exception:
        return {}


def record_signoff(patient_id: str, signer: str, decision: str, comment: str = "") -> dict:
    """Record a molecular-geneticist sign-off (approve | reject | amend)."""
    entry = {"patient_id": patient_id, "signer": signer, "decision": decision,
             "comment": comment, "review_status": f"signed_off:{decision}"}
    with _lock:
        data = _load()
        data[patient_id] = entry
        _SIGNOFFS.parent.mkdir(parents=True, exist_ok=True)
        _SIGNOFFS.write_text(json.dumps(data, indent=2))
    return entry


def get_signoff(patient_id: str) -> dict | None:
    return _load().get(patient_id)
