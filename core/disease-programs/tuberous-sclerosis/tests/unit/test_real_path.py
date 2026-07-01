"""
Verifies the 'real-reasoning' path (PRD §3): when the model is online, agents parse and
merge structured output. Uses a mocked router — no API key, no spend — so the wiring that
activates the moment TSC_ANTHROPIC_API_KEY is set is itself covered by tests.
"""
from __future__ import annotations

from src.agents.base import AgentContext
from src.agents.phenome_mapper.runner import PhenomeMapper
from src.agents.tand_surveillance.runner import TANDSurveillance
from src.agents.therapeutics_strategist.runner import SECTION_KEYS, TherapeuticsStrategist
from src.agents.variant_curator.runner import VariantCurator
from src.utils.hpo import get_ontology
from src.utils.provenance import Provenance

_CANNED = {
    "note_extraction": {"assertions": [
        {"hpo_id": "HP:0009721", "label": "shagreen patch",   # real term (wrong-case label) -> grounded & merged
         "onset": "childhood", "evidence_span": "note span"},
        {"hpo_id": "HP:9999", "label": "llm-extracted finding",  # not in the ontology -> dropped, not trusted
         "onset": "childhood", "evidence_span": "note span"}]},
    "discourse_analysis": {"signals": [{"cluster": "psychiatric", "markers": ["hedging"]}]},
    "acmg_synthesis": {"per_criterion": [{"code": "PVS1", "reasoning": "null variant in LOF gene"}],
                       "summary": "Likely Pathogenic; recommend ddPCR"},
    "options_brief": {k: f"LLM section {k}" for k in SECTION_KEYS},
}


class FakeRouter:
    online = True

    def call_structured(self, agent, step, system, prompt, schema_hint, **kw):
        return _CANNED.get(step), Provenance(agent=agent, step=step, tier="sonnet", latency_ms=1.0)

    def call(self, agent, step, system, prompt, **kw):
        return "fake narrative", Provenance(agent=agent, step=step, tier="sonnet", latency_ms=1.0)


def _patch(monkeypatch, module):
    monkeypatch.setattr(f"src.agents.{module}.runner.get_router", lambda: FakeRouter())


def test_phenome_merges_llm_assertions(monkeypatch):
    _patch(monkeypatch, "phenome_mapper")
    cohort = {"phenotypes": [{"hpo_id": "HP:0002518", "label": "cortical tuber"}],
              "notes": [{"specialty": "neuro", "date": "2025", "text": "note text"}]}
    out = PhenomeMapper().run(AgentContext(patient_id="P", cohort=cohort))
    ids = {t["hpo_id"]: t for t in out.data["hpo_terms"]}
    # a real model-extracted term is grounded (label corrected) and merged as a note assertion;
    # the non-ontology code is dropped rather than trusted, with the count surfaced.
    if get_ontology().available:
        assert "HP:0009721" in ids and ids["HP:0009721"]["source"] == "note"
        assert ids["HP:0009721"]["label"] == "Shagreen patch"  # official label, not the model's casing
        assert "HP:9999" not in ids
        assert out.data["n_dropped_unverified"] >= 1
    else:
        assert "HP:0009721" in ids and ids["HP:0009721"]["source"] == "note"


def test_tand_incorporates_llm_signals(monkeypatch):
    _patch(monkeypatch, "tand_surveillance")
    cohort = {"tand_signals": [], "notes": [{"specialty": "pc", "date": "2025", "text": "t"}]}
    out = TANDSurveillance().run(AgentContext(patient_id="P", cohort=cohort))
    psy = out.data["clusters"]["psychiatric"]
    assert any(str(s).startswith("llm:") for s in psy["sources"])     # LLM signal folded in


def test_variant_curator_attaches_synthesis_narrative(monkeypatch):
    _patch(monkeypatch, "variant_curator")
    cohort = {"genomics": {"gene": "TSC2", "variant": "c.3037C>T p.Arg1013Ter"},
              "phenotypes": [{"hpo_id": "HP:0002518", "label": "cortical tuber"}]}
    out = VariantCurator().run(AgentContext(patient_id="P", cohort=cohort))
    assert out.data["synthesis_narrative"]["summary"].startswith("Likely Pathogenic")


def test_therapeutics_uses_llm_sections(monkeypatch):
    _patch(monkeypatch, "therapeutics_strategist")
    ctx = AgentContext(patient_id="P",
                       projection={"variant_interp": {"primary": {"acmg_classification": "Pathogenic"}}},
                       cohort={"demographics": {"age": 18}})
    out = TherapeuticsStrategist().run(ctx)
    assert out.data["sections"]["current_therapy"] == "LLM section current_therapy"
