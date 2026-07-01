"""
Agent 2 — TSC-Phenome Mapper (PRD §3 FR-PM-*; master paper §9). Runs FIRST.

Builds the longitudinal HPO profile every other agent reads: structured phenotypes
(real, normalized here), per-note HPO extraction (Sonnet — real prompt/contract, live
when keyed), an ITSC surveillance-gap report (real, deterministic), and a discordance
log. The structured + surveillance paths produce real output offline; per-note NLP
enriches the profile when TSC_ANTHROPIC_API_KEY is set.
"""
from __future__ import annotations

from src.agents.base import Agent, AgentContext, AgentOutput
from src.agents.phenome_mapper.itsc import analyze
from src.agents.phenome_mapper.schema import (
    HPOAssertion, PhenomeProfile, SurveillanceItem,
)
from src.utils.hpo import get_ontology
from src.utils.model_router import get_router
from src.utils.prompts import load_prompt


class PhenomeMapper(Agent):
    name = "phenome_mapper"

    def run(self, ctx: AgentContext) -> AgentOutput:
        router = get_router()
        onto = get_ontology()
        provenance = []
        terms: list[HPOAssertion] = []
        seen_ids: set[str] = set()
        mentions: dict[str, dict] = {}        # canonical id -> {label, present:[src], absent:[src]}
        dropped = 0
        bad_spans = 0

        def consider(hpo_id, label, onset, span, source, polarity="present", temporality="current"):
            """Ground one mention against the ontology; record its polarity for discordance
            detection; admit it as a term only if asserted PRESENT (negated mentions are
            logged, never added). Drops codes that are not real HPO terms."""
            nonlocal dropped
            r = onto.resolve(hpo_id, label)
            if r["status"] == "unknown":
                dropped += 1
                return
            cid = r["hpo_id"]
            m = mentions.setdefault(cid, {"label": r["label"], "present": [], "absent": []})
            m["absent" if polarity == "absent" else "present"].append(source)
            if polarity == "absent":
                return
            if cid in seen_ids:
                return
            seen_ids.add(cid)
            terms.append(HPOAssertion(
                hpo_id=cid, label=r["label"], onset=onset, evidence_span=span,
                source=source, validation=r["status"], polarity="present", temporality=temporality,
            ))

        # 1) structured phenotypes -> ontology-grounded HPO assertions (real)
        for p in ctx.cohort.get("phenotypes", []):
            consider(p.get("hpo_id", ""), p.get("label", ""), p.get("onset"),
                     p.get("evidence_span"), "structured",
                     temporality="onset" if p.get("onset") in ("infancy", "childhood") else "current")

        # 2a) per-note span extraction (REAL offline path): every HPO span is verbatim-validated
        #     against the note text, polarity-checked (negated mentions logged not admitted), and
        #     ontology-grounded — span-grounded NLP without requiring a model.
        for note in ctx.cohort.get("notes", []):
            text = note.get("text", "")
            src = f"{note.get('specialty')} {note.get('date')}"
            for sp in note.get("spans", []):
                if sp.get("kind") != "hpo":
                    continue
                if text[sp["start"]:sp["end"]] != sp.get("quote"):    # verbatim-span integrity
                    bad_spans += 1
                    continue
                consider(sp.get("hpo_id"), sp.get("label"), None, sp.get("quote"), src,
                         polarity=sp.get("polarity", "present"),
                         temporality=sp.get("temporality", "current"))

        # 2b) per-note extraction (Sonnet). Offline -> no assertions; when keyed, each parsed
        #     assertion is validated the same way before merge.
        for note in ctx.cohort.get("notes", []):
            parsed, prov = router.call_structured(
                self.name, "note_extraction",
                system=load_prompt("phenome_mapper", "note_extraction"),
                prompt=f"[{note.get('specialty')} {note.get('date')}] {note.get('text')}",
                schema_hint={"assertions": [{"hpo_id": "HP:0000000", "label": "",
                                             "onset": None, "evidence_span": ""}]},
                prompt_template_version="note_extraction/v0",
            )
            provenance.append(prov)
            for a in (parsed or {}).get("assertions", []):
                consider(a.get("hpo_id"), a.get("label"), a.get("onset"),
                         a.get("evidence_span"), "note", polarity=a.get("polarity", "present"))

        # 2c) discordance log (FR-PM-3): a concept asserted both present and absent across the
        #     record (e.g. infantile spasms active at onset, resolved now).
        discordances = [
            {"hpo_id": cid, "label": m["label"], "type": "present/absent conflict",
             "present_sources": m["present"], "absent_sources": m["absent"],
             "note": "feature present in part of the record and explicitly negated elsewhere "
                     "— likely resolved/longitudinal, surfaced for clinician reconciliation"}
            for cid, m in mentions.items() if m["present"] and m["absent"]
        ]

        # 3) ITSC surveillance-gap analysis (real, deterministic)
        labels = [t.label for t in terms]
        gaps = [SurveillanceItem(**g) for g in analyze(labels, ctx.cohort.get("surveillance_gaps", []))]

        profile = PhenomeProfile(
            hpo_terms=terms, surveillance_gaps=gaps, discordances=discordances, n_terms=len(terms),
            ontology_release=(f"HPO release loaded · {onto.n_terms} terms" if onto.available
                              else "unavailable (terms passed through unverified)"),
            n_dropped_unverified=dropped,
        )
        if not provenance:  # ensure at least one provenance record
            _, prov = router.call(
                self.name, "note_extraction", system="(no notes)", prompt="(none)",
                prompt_template_version="note_extraction/v0",
            )
            provenance.append(prov)
        return self.ok(self.name, profile.model_dump(), provenance)
