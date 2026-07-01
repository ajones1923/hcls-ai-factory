# Neurology Agent — Disease Coverage & Agent-vs-Engine Decisions

*HCLS AI Factory · architecture note · 2026-06-29*

## Question

Does the Neurology Agent have room for everything related to **Parkinson's Disease, Alzheimer's
Disease, Dementia, ALS, and Multiple Sclerosis** — or do any of these need to become a separate
Engine or Agent?

## Answer (short)

**All five are already covered by the Neurology Agent, and none *needs* to be a separate service.**
The agent is built as a broad, domain-organized neurology platform with an extensible knowledge
layer. A disease graduates to its own **Engine** only as a deliberate *precision-medicine beachhead*
(the TSC pattern) — never out of architectural necessity.

---

## 1. Current coverage (already in the agent)

| Disease | Modeled as | Evidence in code |
|---|---|---|
| **Parkinson's** | Movement-disorders domain | MDS-UPDRS Part III + Hoehn-Yahr scales · Parkinson's assessment workflow · `neuro_movement` collection |
| **Alzheimer's** | Neurodegenerative domain | A/T/N (amyloid / tau / neurodegeneration) framework · MoCA · dementia-evaluation workflow |
| **Dementia** | Neurodegenerative domain | AD + FTD + Lewy body + prion as key conditions · reversible-dementia workup |
| **ALS** | Neuromuscular domain | ALSFRS-R scale · `neuro_neuromuscular` collection |
| **Multiple Sclerosis** | Inflammatory-demyelinating domain | EDSS + McDonald 2017 criteria · DMT guidelines · dedicated `neuro_ms` collection |

The agent's stated capabilities back this up: **8 clinical workflows** (acute stroke triage,
dementia evaluation, epilepsy classification, brain-tumor grading, MS monitoring, Parkinson's
assessment, headache classification, neuromuscular evaluation) and **10 clinical-scale calculators**
(NIHSS, GCS, MoCA, MDS-UPDRS III, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn-Yahr).

## 2. Why it has room: the architecture is domain-extensible

`knowledge.py` (~1,300 lines) organizes neurology into **10 disease domains / 15 neurodegenerative
diseases / 38+ genes**, each domain carrying its own knowledge entries, guidelines, scales, and a
**dedicated vector collection**:

```
neuro_cerebrovascular   neuro_degenerative   neuro_movement      neuro_ms
neuro_neuromuscular      neuro_epilepsy       neuro_headache      neuro_oncology
neuro_electrophysiology  neuro_imaging        neuro_guidelines    neuro_literature  neuro_trials
```

Adding depth to any disease = **extending its domain + collection**, not standing up a new service.
That is the right shape for a clinical-decision-support agent: breadth across all of neurology.

---

## 3. The decision rule: breadth (Agent) vs depth (Engine)

| | **Agent** (e.g. Neurology) | **Engine** (e.g. TSC) |
|---|---|---|
| Purpose | Clinical decision support across a whole specialty | End-to-end precision medicine for **one** disease |
| Scope | Broad — many diseases | Deep — genomics → variant → target → drug |
| Shape | One RAG + workflows + scales service | Orchestrator + multiple sub-agents + event store |
| Example | Neurology: 10 domains, 8 workflows | TSC: 5 sub-agents (variant curator, trajectory modeler, therapeutics strategist, phenome mapper, TAND surveillance) |
| When to create | Cover a specialty | Make a **strategic beachhead** on a single disease |

A disease becomes an Engine when you want to do for it what the **Tuberous Sclerosis Engine** does
for TSC — a deep, multi-agent, genomics-driven pipeline — *and* there's a business reason (a
beachhead, a partner, a flagship). It is a strategy decision, not a coverage gap.

---

## 4. Per-disease recommendation

- **Dementia → keep in the agent; do not separate.** It is an umbrella over AD / FTD / Lewy body /
  vascular, all already in the neurodegenerative domain. Note: the platform's **flagship demo *is* a
  dementia** — frontotemporal dementia via the *VCP* gene + CB-5083 — already running end-to-end
  across the three core engines (Genomics → Precision Intelligence → Therapeutic Discovery). Dementia
  already has the deepest cross-platform representation in the factory.

- **Alzheimer's, Parkinson's, ALS → keep in the agent now; candidates for a future Engine.** These
  three are the strongest candidates *if* you decide on a deep precision-medicine flagship, because
  each has rich, actionable genetics and an active drug pipeline:
  - Parkinson's — *LRRK2 / GBA / SNCA*
  - Alzheimer's — *APOE / APP / PSEN1 / PSEN2*
  - ALS — *C9orf72 / SOD1 / TARDBP / FUS* (antisense era: tofersen)

  Promote one to `engines/<disease>/` only when it earns beachhead status; otherwise the agent
  covers it.

- **Multiple Sclerosis → keep in Neurology, but cross-link to the Precision Autoimmune Agent.** MS is
  the one cross-cutting case: clinically neurology (EDSS, DMTs) but mechanistically
  autoimmune/demyelinating. Rather than duplicate the immunology, share MS knowledge between the two
  agents over the cross-agent event bus.

---

## 5. Bottom line

Keep all five diseases in the **Neurology Agent** — it already models them and has room to deepen any
of them. Reserve **Engine** status for a genuine end-to-end precision-medicine play on a *single*
disease (most likely **Alzheimer's, Parkinson's, or ALS**), following the TSC pattern. This keeps the
"Eight Agents" tier intact and makes "Engine" mean real depth.

### Quick reference

| Disease | Today | Could become its own Engine? |
|---|---|---|
| Parkinson's | Neurology Agent (movement domain) | Yes — strong precision-medicine candidate (beachhead) |
| Alzheimer's | Neurology Agent (neurodegenerative) | Yes — strong precision-medicine candidate (beachhead) |
| Dementia | Neurology Agent (umbrella) + FTD flagship demo | No — umbrella, not a single disease |
| ALS | Neurology Agent (neuromuscular) | Yes — strong precision-medicine candidate (beachhead) |
| Multiple Sclerosis | Neurology Agent (demyelinating) | No — keep as Agent; cross-link to Precision Autoimmune |
