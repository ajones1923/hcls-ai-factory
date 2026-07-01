"""
The TSC-Orchestrator (PRD §2.6): a DETERMINISTIC event router — not an LLM, and not
a service the surfaces call directly. It dispatches agents in dependency order,
appends every output/decision to the event log, materializes the projection, and
assembles surfaces on demand from projections (never from live agent calls). A failed
agent writes `agent_failed` and leaves the surface section "pending" — never a silent gap.

(LangGraph is the intended production runtime per the PRD; the deterministic logic is
implemented here directly so the engine runs with no extra dependency. The graph
definition is a thin wrapper to add later.)
"""
from __future__ import annotations

from config.settings import settings
from src.agents.base import Agent, AgentContext
from src.orchestrator.ephemeral import Ephemeral, make_ephemeral
from src.orchestrator.events import Event, EventType
from src.orchestrator.policies import default_agents
from src.orchestrator.state import EventStore, make_store
from src.utils.provenance import Provenance

# agents whose only new input is a clinical note — the incremental-update set (FR-OR-3)
_NOTE_AGENTS = ("phenome_mapper", "tand_surveillance")
_DEFAULT_CLINICIAN = "tsc-clinic"


class Orchestrator:
    def __init__(self, store: EventStore | None = None, agents: list[Agent] | None = None,
                 ephemeral: Ephemeral | None = None) -> None:
        self.store = store or make_store()
        self.agents = agents or default_agents()
        self.ephemeral = ephemeral or make_ephemeral()

    # ── ingest ──────────────────────────────────────────────────────────
    def enroll(self, patient_id: str, cohort: dict | None = None) -> None:
        cohort = cohort or {}
        self.store.append(Event(patient_id=patient_id, event_type=EventType.patient_enrolled,
                                payload={"cohort_keys": sorted(cohort.keys())}))
        # In the demo the genomic substrate is the synthetic BAM/VCF, ready at enrollment.
        self.store.append(Event(patient_id=patient_id, event_type=EventType.genomic_substrate_ready,
                                payload={"source": "synthetic", "ref": f"tsc-cohort/{patient_id}"}))
        self.dispatch(patient_id, cohort)

    # ── dependency-ordered dispatch (PRD §2.6) ──────────────────────────
    def dispatch(self, patient_id: str, cohort: dict) -> None:
        for agent in self.agents:
            self.run_agent(patient_id, agent, cohort)

    def ingest_note(self, patient_id: str, note: dict, cohort: dict) -> list[str]:
        """Incremental update (FR-OR-3): a new clinical note re-runs ONLY the note-dependent
        agents (Phenome Mapper, then TAND) — not variant calling, trajectory, or therapeutics —
        and re-materializes the projection. Returns the agents re-run, in dependency order."""
        cohort.setdefault("notes", []).append(note)
        ran = []
        for agent in self.agents:                       # self.agents is already dependency-ordered
            if agent.name in _NOTE_AGENTS:
                self.run_agent(patient_id, agent, cohort)
                ran.append(agent.name)
        return ran

    def deliver_alert(self, clinician_id: str = _DEFAULT_CLINICIAN) -> dict:
        """Record one alert delivered to a clinician in the sliding weekly window (FR-SF-3)."""
        count = self.ephemeral.alert_window_incr(clinician_id)
        budget = settings.ALERTS_PER_CLINICIAN_PER_WEEK_MAX
        return {"clinician_id": clinician_id, "week_count": count, "weekly_budget": budget,
                "over_budget": count > budget, "recalibrate": count > budget}

    def run_agent(self, patient_id: str, agent: Agent, cohort: dict) -> None:
        """Run one agent against the current projection and append its output (or failure).
        Shared by the plain dispatcher and the LangGraph runtime."""
        if not self._deps_met(patient_id, agent):
            self._fail(patient_id, agent.name, "unmet dependencies")
            return
        ctx = AgentContext(patient_id=patient_id,
                           projection=self.store.projection(patient_id),
                           cohort=cohort)
        try:
            out = agent.run(ctx)
        except Exception as exc:  # conservative failure handling
            self._fail(patient_id, agent.name, repr(exc))
            return
        self.store.append(Event(
            patient_id=patient_id, event_type=agent.emits,
            payload=out.data, provenance={"records": out.provenance},
        ))

    def _deps_met(self, patient_id: str, agent: Agent) -> bool:
        have = {e.event_type for e in self.store.events_for(patient_id)}
        return all(dep in have for dep in agent.depends_on)

    def _fail(self, patient_id: str, agent_name: str, error: str) -> None:
        self.store.append(Event(
            patient_id=patient_id, event_type=EventType.agent_failed,
            payload={"agent": agent_name, "error": error},
            provenance=Provenance(agent=agent_name, step="dispatch").to_dict(),
        ))

    # ── demand-driven surface assembly (PRD §2.6) ───────────────────────
    def assemble_surface(self, patient_id: str, kind: str) -> dict:
        self.store.append(Event(patient_id=patient_id, event_type=EventType.surface_requested,
                                payload={"kind": kind}))
        proj = self.store.projection(patient_id)
        builder = {"briefing": self._briefing, "dashboard": self._dashboard, "alerts": self._alerts}[kind]
        payload = builder(patient_id, proj)
        payload["watermark"] = "SYNTHETIC"
        self.store.append(Event(patient_id=patient_id, event_type=EventType.surface_assembled,
                                payload={"kind": kind, "staleness": proj["staleness"]}))
        return payload

    # surface builders read projections only
    def _briefing(self, patient_id: str, proj: dict) -> dict:
        vi = (proj.get("variant_interp") or {}).get("primary") or {}
        traj = (proj.get("trajectory") or {}).get("lesions", [])
        tand = proj.get("tand_briefing") or {}
        actions = []
        if vi.get("recovered"):
            actions.append(
                f"Somatic mosaic {vi.get('gene')} variant recovered in tissue "
                f"(VAF {vi.get('vaf')}) — {vi.get('acmg_classification')}; "
                "pending molecular-geneticist review. (Blood testing was negative.)"
            )
        for les in traj:
            if les.get("crosses_in_12_18mo_window"):
                actions.append(f"{les['lesion']} forecast approaches discussion threshold "
                               f"in ~{les['months_to_threshold']} months — plan discussion.")
        if tand.get("flagged_clusters"):
            actions.append(f"TAND signals aggregating ({', '.join(tand['flagged_clusters'])}) — "
                           "consider formal evaluation. (Briefing item, not an alert.)")
        return {
            "kind": "briefing", "patient_id": patient_id,
            "header": {"genotype": vi.get("gene"), "variant": vi.get("hgvsc"),
                       "classification": vi.get("acmg_classification")},
            "action_items": actions[:3],
            "staleness": proj["staleness"],
        }

    def _dashboard(self, patient_id: str, proj: dict) -> dict:
        return {
            "kind": "dashboard", "patient_id": patient_id,
            "quadrants": {
                "variant_interpretation": proj.get("variant_interp"),
                "hpo_timeline": proj.get("hpo_profile"),
                "trajectory": proj.get("trajectory"),
                "tand_and_therapeutics": {
                    "tand": proj.get("tand_briefing"),
                    "therapeutics": proj.get("therapeutics"),
                },
            },
            "staleness": proj["staleness"],
        }

    def _alerts(self, patient_id: str, proj: dict) -> dict:
        alerts = []
        for les in (proj.get("trajectory") or {}).get("lesions", []):
            if les.get("crosses_in_12_18mo_window"):
                alerts.append({
                    "category": "trajectory_threshold_crossing",
                    "summary": f"{les['lesion']} forecast crosses threshold in "
                               f"~{les['months_to_threshold']} months",
                    "source_section": "trajectory", "dismissable": True,
                })
        # alert discipline (PRD §2.4 / FR-SF-3): track the per-clinician sliding weekly window
        budget = settings.ALERTS_PER_CLINICIAN_PER_WEEK_MAX
        week_count = self.ephemeral.alert_window_count(_DEFAULT_CLINICIAN)
        return {
            "kind": "alerts", "patient_id": patient_id, "alerts": alerts,
            "clinician_id": _DEFAULT_CLINICIAN, "weekly_budget": budget,
            "clinician_week_count": week_count,
            "over_budget": week_count > budget,
            "recalibrate": week_count > budget,
            "discipline_note": "Surfaced only with clear action + source; "
                               f"recalibrate if >{budget}/clinician/week.",
            "staleness": proj["staleness"],
        }
