"""
Agent 4 — TAND Surveillance Agent (PRD §3 FR-TS-*; master paper §11).

Surfaces under-recognized TSC-Associated Neuropsychiatric Disorders from the longitudinal
record. Per-note discourse analysis over the six TAND clusters applying the
Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy; a deterministic
aggregation/scoring layer; an Opus briefing. Surfaced as pre-visit BRIEFING MATERIAL —
never interruptive alerts, never a diagnosis. The discourse-marker detection and scoring
are real offline; the Sonnet per-note pass and Opus briefing go live when keyed.
"""
from __future__ import annotations

from src.agents.base import Agent, AgentContext, AgentOutput
from src.agents.tand_surveillance.scoring import score
from src.agents.tand_surveillance.taxonomy import TAXONOMY_VERSION
from src.utils.model_router import get_router
from src.utils.prompts import load_prompt


class TANDSurveillance(Agent):
    name = "tand_surveillance"

    def run(self, ctx: AgentContext) -> AgentOutput:
        router = get_router()
        signals = ctx.cohort.get("tand_signals", [])
        notes = ctx.cohort.get("notes", [])

        # Sonnet per-note discourse analysis (real prompt/contract; live when keyed)
        provenance = []
        extra_signals = []
        for note in notes:
            parsed, prov = router.call_structured(
                self.name, "discourse_analysis",
                system=load_prompt("tand_surveillance", "discourse_analysis"),
                prompt=f"[{note.get('specialty')}] {note.get('text')}",
                schema_hint={"signals": [{"cluster": "academic", "markers": ["hedging"]}]},
                prompt_template_version=f"discourse_analysis/{TAXONOMY_VERSION}",
            )
            provenance.append(prov)
            for s in (parsed or {}).get("signals", []):
                if s.get("cluster"):
                    extra_signals.append({"cluster": s["cluster"], "source": f"llm:{note.get('date')}"})

        # deterministic aggregation/scoring over structured + (when keyed) LLM-derived signals
        scored = score(signals + extra_signals, notes)
        flagged = scored["flagged_clusters"]

        # Opus briefing summary (real prompt; live when keyed)
        _, prov_brief = router.call(
            self.name, "briefing_summary",
            system="Summarize aggregated TAND signals as clinician briefing material "
                   "(not an alert, not a diagnosis), with explicit confidence.",
            prompt=str(scored), prompt_template_version="briefing_summary/v0",
        )
        provenance.append(prov_brief)

        recommendation = (
            f"Pattern across {', '.join(flagged)} cluster(s) — consider TAND-L re-administration / "
            "developmental-behavioral referral. (Briefing item for clinician review, not an alert.)"
            if flagged else "No aggregated pattern crossing the briefing threshold."
        )
        data = {
            "clusters": scored["clusters"],
            "flagged_clusters": flagged,
            "discourse_highlights": scored["highlights"],
            "taxonomy_version": TAXONOMY_VERSION,
            "surfaced_as": "pre_visit_briefing",   # never an alert, never a diagnosis
            "recommendation": recommendation,
            "_build_status": "real discourse-marker detection + scoring; Sonnet/Opus live when keyed",
        }
        return self.ok(self.name, data, provenance)
