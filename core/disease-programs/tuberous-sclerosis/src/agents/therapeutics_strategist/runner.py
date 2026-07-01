"""
Agent 5 — TSC-Therapeutics Strategist (PRD §3 FR-TR-*; master paper §12). The most senior
agent; Opus-class is non-negotiable. Integrates all four prior agents (from the projection)
+ medications/adverse events + RAG over the TSC literature + a ClinicalTrials.gov snapshot
into a SIX-SECTION options brief, every claim source-attributed. Decision support, never a
treatment recommendation. RAG retrieval and trial matching are real here; the Opus narrative
goes live when keyed (offline, sections are assembled from real upstream data + retrieved sources).
"""
from __future__ import annotations

from src.agents.base import Agent, AgentContext, AgentOutput
from src.agents.therapeutics_strategist.trials import match_trials
from src.rag.retriever import get_retriever
from src.utils.model_router import get_router
from src.utils.prompts import load_prompt

SECTION_KEYS = [
    "current_therapy", "optimization", "combination",
    "trial_matching", "emerging_evidence", "open_questions",
]


def _profile(ctx: AgentContext) -> dict:
    proj = ctx.projection
    vi = (proj.get("variant_interp") or {}).get("primary") or {}
    demo = ctx.cohort.get("demographics", {})
    cls = vi.get("acmg_classification")
    seizures = (ctx.cohort.get("demographics", {}).get("severity") in ("moderate-severe", "severe"))
    return {
        "molecular_dx_confirmed": cls in ("Pathogenic", "Likely Pathogenic"),
        "age": demo.get("age"),
        "refractory_epilepsy": ctx.cohort.get("refractory_epilepsy",
                                              True if "refractory" in str(ctx.cohort).lower() else None),
        "aml_max_cm": ctx.cohort.get("aml_max_cm"),
        "classification": cls,
        "trajectory": proj.get("trajectory"),
        "tand": proj.get("tand_briefing"),
        "therapy": ctx.cohort.get("therapy") or (ctx.cohort.get("trial_matches") and "see record"),
    }


class TherapeuticsStrategist(Agent):
    name = "therapeutics_strategist"

    def run(self, ctx: AgentContext) -> AgentOutput:
        prof = _profile(ctx)
        retriever = get_retriever()

        # real RAG retrieval (emerging evidence) with exact source attribution
        evidence = retriever.retrieve(
            f"TSC therapeutics mTOR inhibitor trials {prof.get('classification')}", k=3
        )
        sources = [{"source_uri": e["source_uri"], "pub_year": e["pub_year"],
                    "section": e["section"], "score": e["score"]} for e in evidence]

        # real trial matching
        trial_matches = match_trials(prof)

        # six-section brief assembled from real upstream data + retrieved sources
        sections = {
            "current_therapy": f"Current regimen: {prof.get('therapy') or 'not documented'}. "
                               f"Molecular diagnosis: {prof.get('classification') or 'unconfirmed'}.",
            "optimization": "Within-class options (mTOR inhibitor dosing, anti-seizure regimen) "
                            "to weigh against tolerability; see cited evidence.",
            "combination": "Cross-class considerations (e.g., adjunct anti-seizure therapy; "
                           "AML-directed intervention if approaching threshold).",
            "trial_matching": "; ".join(f"{t['nct']}: {t['match']}" for t in trial_matches),
            "emerging_evidence": " | ".join(f"{e['text'][:120]}... [{e['source_uri']}]" for e in evidence),
            "open_questions": "Items requiring clinician judgment and family discussion "
                              "(verify stable-regimen duration; confirm surveillance cadence).",
        }

        # Opus synthesis narrative (real prompt; live when keyed)
        router = get_router()
        parsed, prov = router.call_structured(
            self.name, "options_brief",
            system=load_prompt("therapeutics_strategist", "options_brief"),
            prompt=f"profile={prof}\ntrials={trial_matches}\nsources={sources}",
            schema_hint={k: "..." for k in SECTION_KEYS},
            prompt_template_version="options_brief/v0",
            rag_source_uris=[e["source_uri"] for e in evidence],
        )
        if parsed:   # real model prose grounded in the retrieved sources, when keyed
            for k in SECTION_KEYS:
                if parsed.get(k):
                    sections[k] = parsed[k]

        data = {
            "sections": sections,
            "trial_matches": trial_matches,
            "sources": sources,
            "embedding_mode": retriever.embedding_mode,
            "framing": "decision-support, not treatment recommendation",
            "_build_status": "real RAG + trial matching; Opus narrative live when keyed",
        }
        return self.ok(self.name, data, [prov])
