# TSC Intelligence Engine — Product Requirements Document (v0.1)

### HCLS AI Factory, Engine 7 · Cincinnati Children's Demonstration Build

**Owner:** Adam M. Jones (architect, sole engineer for the v0.1 build)
**Target:** A working, end-to-end, 5-agent + orchestrator + 3-surface demonstration on a single NVIDIA DGX Spark (RunPod for extra GPU as needed), against a synthetic 50-patient TSC cohort.
**License:** Apache 2.0 · **Companion to:** *The TSC Intelligence Engine* (master research volume) · 2026

> **v0.1 = synthetic-data demonstration build.** No real patient data, no IRB, no institutional commitment. Epic/Clarity/Caboodle, biobank LIMS, imaging-AI/PACS, IAM, and institutional audit/SIEM are specified for the institutional path but are **out of scope for v0.1** and explicitly stubbed/simulated. Acceptance is defined by the criteria in §6.

---

## Implementation Status — v0.1 Built & Verified (June 2026)

This PRD has been **executed**. A working v0.1 is implemented and running on the NVIDIA DGX Spark.

- **Acceptance criteria AC-1 through AC-7 (the buildable set) are met.** The synthetic 50-patient cohort is generated and version-controlled; all five agents are operational against it; the orchestrator routes the 13 event types and maintains projections; the three clinician surfaces render; the audit/provenance trail is complete and queryable; the demo runbook executes end-to-end; the repository (code, config, prompts, cohort pipeline, this PRD) is complete. AC-8 (clinical review of sampled notes) and AC-9 (delivery at Cincinnati Children's) are inherently the institution's and remain.
- **~63 Python files, ~2,870 LOC; 41 automated tests pass in under one second.** Run: `python3 -m pytest tests/unit tests/eval -q` and `python3 scripts/dry_run_demo.py`.
- **Deterministic/classical FRs are fully real** — the ACMG-AMP combinatorial classifier (FR-VC-5/6), the Gaussian-process trajectory intervals (FR-TM-5), the Marshall-Hagedorn discourse detection (FR-TS-3), the ITSC surveillance-gap analyzer (FR-PM-7), and the trial eligibility matcher (FR-TR-3).
- **The LLM FRs reason for real when keyed** — verified live: the Variant Curator produced a genuine `claude-opus-4-8` ACMG narrative; the deterministic validator stays authoritative. The test suite is pinned offline (`TSC_OFFLINE=1`) so it never incurs API cost.
- **Documented upgrades still open (as the PRD anticipates):** the BAMSurgeon→Parabricks blinded-calling genomic substrate (RunPod GPU), and frontier-model cohort notes/imaging *generation*. The Phase-1 integration boundaries (Epic/LIMS/PACS/IAM/SIEM) remain out of v0.1 by design.

### Update — expanded build (June 2026)

Subsequent work has carried the engine well beyond the v0.1 walking skeleton; the file/test counts above are superseded. The repository is now **~86 application Python files (~5,700 LOC) with a 92-test suite (offline, deterministic)**, and the following are implemented and verified on the Spark:

- **Ontology-grounded phenome.** Every HPO term is validated against the official Human Phenotype Ontology release (19,389 terms): labels corrected to the canonical term, obsolete/alt IDs remapped, and codes that are not real HPO terms dropped rather than trusted (this caught and corrected two miscoded terms in the synthetic cohort itself).
- **A measured evaluation harness** (`/eval`, `scripts/run_eval.py`): variant-classification accuracy 1.00 (construct validity), a **+12-point diagnostic detection uplift** driven by 100% recovery of the six sub-threshold mosaics, zero truncating-variants-called-benign, and 100% provenance completeness — each carrying an explicit "construct validity on synthetic data, not prospective accuracy" caveat.
- **Read-level genomics (FR-VC-4).** VCFs carry strand-resolved allelic depths + a strand-balance, and every sample includes a low-VAF strand-biased artifact the curator must reject — so "zero false-positive pathogenic" is a real discrimination. The curator surfaces depth, alt-reads, strand balance, and an artifact assessment.
- **Multi-quantity trajectories (FR-TM-1/2/4/5, config-driven).** SEGA/AML size plus **renal function (eGFR)** and **seizure frequency**, with a Bayesian growth slope (population-prior shrinkage for SEGA, conjugate weak-prior otherwise), a survival-style threshold-crossing probability per horizon, a graded crossing (likely/possible/unlikely), and a surveillance-cadence recommendation that tightens vs. the ITSC floor.
- **Span-grounded clinical NLP (FR-CG-3, FR-PM-2/3, FR-TS-2/3).** A deterministic note generator embeds HPO + TAND markers at recorded character offsets with polarity/temporality; the Phenome Mapper extracts verbatim-validated spans (negated mentions logged, never admitted) and a present/absent discordance log; the TAND agent scores from spans with configurable weights, recency, recurrence, and an explicit negation filter.
- **Governance & operations.** A per-tier token + USD cost ledger with a spend cap (`/system/usage`, NFR-COST); an AI-labeled DRAFT report + molecular-geneticist sign-off endpoint (FR-VC-5); incremental dispatch — a new note re-runs only Phenome + TAND (FR-OR-3); a per-clinician sliding-window alert counter + recalibration (FR-SF-3).
- **Production substrate now real with runs-here fallbacks:** statsmodels mixed-effects, the LangGraph runtime (verified byte-equivalent), and PostgreSQL/Redis/MinIO/Milvus behind env flags, with SQLite/in-memory/filesystem defaults so the engine runs on a bare Spark.
- **Omniverse / OpenUSD digital-twin layer (Surface d).** The Spark authors four scene kinds from real projections — lesion-trajectory twin (envelope radii = the prediction intervals), mosaic "powers-of-ten" cell field (exactly the recovered VAF fraction glows), whole-child organ atlas (organs lit by the phenome), and a 50-patient population array (seven recovered mosaics ringed in gold) — with MDL materials for film-quality RTX. Authoring is CPU-side and dependency-free; rendering is a RunPod RTX step.

**The honest edges still hold.** All data remains synthetic and watermarked. The faithful BAMSurgeon→Parabricks blinded-calling substrate needs a burst GPU and is unbuilt; Synthea/FHIR + Epic/Clarity, the biobank LIMS, imaging-AI/PACS, IAM, and SIEM remain Phase-1 work; frontier-model note/imaging *generation* and patient-specific (segmented) anatomy are documented upgrades; the RTX render runs off-box. Nothing is FDA-cleared. Prospective validation on real Cincinnati specimens is the next step.

---

## Table of Contents

1. Overview, Goals, Scope & Users
2. System Architecture, Repo Layout & Deployment
3. Functional Requirements — Cohort Generation, Variant Curator, Phenome Mapper
4. Functional Requirements — Trajectory Modeler, TAND, Therapeutics, Orchestrator, Surfaces
5. Non-Functional Requirements, Data Specs, Evaluation Harness & Integration
6. Build Plan, Test Plan, Risk Register, Dependencies & Roadmap

---


## 1. Overview, Goals, Scope & Users

The TSC Intelligence Engine is Engine 7 of the HCLS AI Factory — the open-source precision-medicine platform whose headline is now "Seven Engines. Eight Intelligence Agents. One Platform." Where the other six engines (Genomic Foundation, Precision Intelligence, Therapeutic Discovery, Clinical Imaging, Precision Oncology, Cardiology Intelligence) each address a broad capability or organ system, Engine 7 is purpose-built for a single rare, multisystem disease: Tuberous Sclerosis Complex (TSC). It earns the "engine" designation rather than "agent" because it is not one model behind one surface. It is five coordinated agents plus a deterministic orchestrator, spanning genomics, longitudinal phenotyping, classical forecasting, clinical-NLP surveillance, and therapeutic strategy — cross-modal and coordinating, with shared event-sourced state feeding three clinician surfaces.

This document is the Product Requirements Document for the **v0.1 demonstration build**. It defines what we are building first, what we are explicitly not building, who it is for, and how we will know it worked. The remaining sections (beginning with §2 Architecture) specify the system in detail; this section sets the boundary conditions that every later section must respect.

A guiding principle runs through the whole PRD: **the build is the argument.** The audience for this work — clinical informaticians and physician-scientists who have seen many demos that did not survive contact with reality — is persuaded by a working system that is honest about its edges, not by claims. So this section is deliberately conservative about scope. We would rather ship a smaller thing that is exactly true than a larger thing that requires footnotes to defend.

### 1.1 Product summary

**TSC Intelligence Engine v0.1** is a working, single-node-plus-burst demonstration of a multi-agent precision-medicine engine for TSC, running on synthetic patient data. It takes a 50-patient synthetic TSC cohort — genomic substrate (BAM-level), longitudinal structured records, clinical notes, and imaging reports, all generated and version-controlled — and runs it through five agents and one orchestrator to produce clinician-facing decision-support surfaces:

- **Agent 1 — TSC-Variant Curator**: deterministic BAM→VCF calling (mosaic-aware, low-VAF ≥5%), deterministic annotation, evidence aggregation, and ACMG-AMP classification synthesis. Recovers low-VAF somatic mosaic variants that standard blood testing misses — the ~10–15% "no-mutation-identified" (NMI) cohort. Output is an AI-labeled draft molecular-genetics report for board-certified molecular-geneticist sign-off, never autonomous.
- **Agent 2 — TSC-Phenome Mapper**: extracts time-anchored, HPO-coded phenotypes from synthetic Epic Clarity-shaped structured data and clinical notes; produces a longitudinal HPO profile with evidence spans and an International TSC (ITSC) surveillance-gap report. It is the foundation the other agents build on.
- **Agent 3 — TSC-Trajectory Modeler**: classical statistics (mixed-effects, Gaussian process regression, survival analysis, Bayesian updating) — not LLM-driven — forecasting SEGA/AML growth, seizure burden, and renal function at 6/12/18 months with 50%/90% prediction intervals and threshold-crossing alerts.
- **Agent 4 — TAND Surveillance Agent**: surfaces under-recognized TSC-Associated Neuropsychiatric Disorders (TAND) from longitudinal notes using the Marshall-Hagedorn diagnostic-uncertainty discourse-marker taxonomy. Surfaces as pre-visit briefing material, never interruptive alerts, never diagnoses.
- **Agent 5 — TSC-Therapeutics Strategist**: integrates the four prior agents plus medications/adherence, adverse events, PubMed/PMC RAG, a ClinicalTrials.gov snapshot, and FDA actions into a six-section structured options brief — every claim source-attributed, framed as decision-support.
- **TSC-Orchestrator**: a deterministic LangGraph event router (not an LLM) over an append-only PostgreSQL event log with materialized current-state projections, coordinating agent execution and assembling surfaces on demand.

The engine reuses the HCLS AI Factory v1.3.0 substrate (LangGraph, Milvus RAG, tiered Claude models with a local Llama 3.1 70B fallback, PostgreSQL/Redis/MinIO, audit/provenance) and is scaffolded from the existing `precision_oncology_agent` template (Streamlit UI + FastAPI + Milvus + Pydantic + Claude). It is Apache 2.0, open, with non-commercial intent.

### 1.2 Goals

1. **Demonstrate mosaic recovery as a credible, auditable capability.** Show, live, that the Variant Curator recovers a low-VAF somatic mosaic variant in tissue (Patient A: 8.3% VAF *TSC2* frameshift in tuber tissue) that blood-based testing would report as NMI, classifies it correctly under ACMG-AMP (Likely Pathogenic, PVS1+PM2+PP4), recommends ddPCR validation, and exposes a complete provenance trail.
2. **Demonstrate cross-agent coordination on a single patient.** Show that one patient's genomic, phenotypic, trajectory, TAND, and therapeutic outputs assemble into coherent clinician surfaces (Patient B's four-quadrant in-visit dashboard; Patient C's therapeutics brief; a one-screen pre-visit briefing) through the deterministic orchestrator.
3. **Demonstrate that the TAND agent is a direct extension of published CCHMC research**, not an external graft — the Marshall-Hagedorn diagnostic-uncertainty methodology applied to TSC longitudinal notes.
4. **Demonstrate the cost and infrastructure thesis**: a clinically interesting multi-agent engine running first on a single DGX Spark with burst-to-RunPod for GPU-heavy steps, reproducibly, on version-controlled synthetic data.
5. **Produce a reusable engine, not a one-off demo.** The orchestrator, shared event-sourced state, synthetic-cohort pipeline, and agent contracts are written so that "swap the box labels, keep the wiring" replication to other institutions and adjacent diseases (NF1/NF2, Rett, Williams, mTORopathies) is a configuration exercise, not a rewrite.

### 1.3 Non-goals

These are deliberately out of scope for v0.1. Several are real, valuable, and architecturally described elsewhere in the platform — but they are **institutional Phase-1 work**, not demonstration-build work.

- **Not** live Epic Clarity/Caboodle, biobank LIMS, or any production data-plumbing integration. The synthetic cohort is shaped *like* Clarity output; no connector is built.
- **Not** an imaging-AI pipeline. Imaging is represented as frontier-model-generated radiology *reports*; there is no DICOM ingestion, segmentation, or volumetric model.
- **Not** clinical validation. The eval targets in §1.5 measure performance against *synthetic ground truth*, which is a software-correctness test, not evidence of clinical accuracy.
- **Not** FDA-cleared, and the SaMD posture is undetermined — that determination is institutional work.
- **Not** an Epic-embedded product. Surfaces are standalone web apps, persistently watermarked synthetic.
- **Not** autonomous clinical decision-making. Every agent output is decision-support behind a human gate.

### 1.4 Success criteria

v0.1 is successful if, at the demo (target early Q3 2026), the following are simultaneously true:

| # | Criterion | How verified |
|---|-----------|--------------|
| S1 | Mosaic recovery runs live and correctly on Patient A | Live run; ACMG class matches ground truth; ddPCR rec present; audit trail opened |
| S2 | All five agents produce coherent, source-attributed output for Patients B and C | Walkthrough of four-quadrant dashboard + therapeutics brief |
| S3 | Every agent output carries full provenance | Provenance record present and queryable for each output |
| S4 | The engine reproduces from a clean checkout on the Spark+RunPod | Cohort regen (~12 hr) + dry runs in W8 |
| S5 | The do-not-overclaim boundary holds throughout the demo | Every deferred capability is marked as Phase-1 in the narrative |
| S6 | Clinician reviewer (CCHMC TSC lead or designee) finds the outputs clinically legible | W8 clinician review pass |

Quantitative per-agent eval targets (vs synthetic ground truth) are owned by the evaluation section and summarized here for orientation: Variant Curator recovers all 7 mosaic variants at VAF ≥5% with correct ACMG class, no false-positive Pathogenic, <5 min/case; Phenome Mapper recall ≥90% / precision ≥85%, full cohort <1 hr; Trajectory Modeler forecasts Patient B's SEGA crossing threshold within a 12–18-month window with no false alarms; TAND detects embedded signals with no spurious flags; Therapeutics Strategist produces correct trial matches with appropriate hedging and full attribution, <3 min.

### 1.5 Scope: in v0.1 vs deferred to institutional Phase-1

| Capability | v0.1 (demo, now) | Institutional Phase-1 (later) |
|------------|------------------|-------------------------------|
| Patient data | 50-patient synthetic cohort, version-controlled | Real CCHMC cohort under IRB |
| Genomics entry point | BAM (BAMSurgeon-inserted variants on NA12878-derived BAMs) | Raw FASTQ from sequencing core |
| Variant calling | Parabricks-equivalent, mosaic-aware | Production validated pipeline |
| Structured data | Synthetic, Clarity-shaped | Live Epic Clarity/Caboodle + biobank LIMS |
| Notes / imaging | Frontier-model-generated, watermarked | Real notes; real imaging-AI pipeline |
| Surfaces | Standalone watermarked web apps | Integrated clinician workflow |
| Regulatory | Synthetic; no IRB required | IRB; SaMD posture determined |

The institutional context — the Winslow Research Pavilion as the infrastructure envelope, and the five CCHMC source areas (Discover Together Biobank, Biomedical Informatics, the TSC clinical/research program, and the Epic Clarity/Caboodle + LIMS data plumbing) that feed the engine — is the destination this build argues toward. v0.1 builds the engine; Phase-1 wires it to those sources. As the project mantra puts it: a biobank without an intelligence layer is a freezer full of tubes. This PRD specifies the intelligence layer; it does not pretend the wiring is already done.

### 1.6 Users and personas

**Adam M. Jones — builder (primary, now).** Architect of the HCLS AI Factory; designs, implements, runs, and demos the engine on the DGX Spark with RunPod burst. Needs reproducible builds, version-controlled synthetic data, clear agent contracts, and an orchestrator that fails conservatively (a down agent yields a "pending" surface tile, never silent missing output).

**CCHMC leadership and Dr. Philip A. Hagedorn's informatics team — demo audience (primary, near-term).** Skeptical, sophisticated clinical and informatics readers. They need to see that the TAND agent extends Hagedorn's own published methodology (Marshall et al. 2023; Nickels et al. 2024), that mosaic recovery addresses the real NMI gap, and that every claim is bounded and attributable. They will be persuaded by the audit trail and the explicit do-not-overclaim discipline, not by polish. Hagedorn has offered to engage his team and the faculty TSC lead.

**Open-source community — adopters (ongoing).** Engineers and informaticians cloning the Apache 2.0 repo. They need a clean scaffold from the `precision_oncology_agent` template, documented agent interfaces, the synthetic-cohort pipeline (Synthea + BAMSurgeon + frontier-model notes/reports), and a "swap the box labels" replication path.

**Future institutional clinicians — end users (Phase-1, not v0.1).** TSC clinicians who would eventually consume the pre-visit briefing, in-visit dashboard, and async alert surface against real data. They are designed *for* but not deployed *to* in v0.1; their needs shape the surface design (one-screen mobile-readable briefings, progressive disclosure, source navigation, ≤3 alerts/clinician/week discipline) without being a delivery target.

With users, goals, scope, and success criteria fixed, the next section specifies how the five agents, the deterministic orchestrator, the shared event-sourced state, and the clinician surfaces fit together into a single engine.


## 2. System Architecture, Repo Layout & Deployment

Section 1 established what the TSC Intelligence Engine is and why it earns the designation of Engine 7 rather than a single agent: five coordinated agents, a deterministic orchestrator, cross-modal inputs, and three clinician surfaces. This section is the engineering counterpart. It describes the components and how they wire together, the repository we will actually scaffold, the services and ports, the data model that holds the event-sourced state together, and the deployment topology that lets the whole thing run first on a single NVIDIA DGX Spark with RunPod GPUs added only where they pay for themselves.

The governing constraint is honesty about scope. The synthetic-data demo is what runs now. Epic Clarity / Caboodle extraction, the biobank LIMS feed, and any imaging-AI pipeline are described here architecturally because the architecture has to anticipate them, but they are explicitly not built in the demo. Where a component is institutional Phase-1 work rather than W1-W8 build work, it is marked. The build is the argument, and the build has edges.

### 2.1 Component architecture

The engine is a directed pipeline with a fan-in. The TSC-Phenome Mapper runs first because every other agent reads against the longitudinal HPO profile it produces. The Variant Curator runs in parallel off the genomic substrate. Trajectory Modeler, TAND Surveillance, and the Therapeutics Strategist consume upstream outputs. The orchestrator is not in the data path of any single agent; it is an event router that decides which agent runs when, persists state, and assembles surfaces on demand.

```
                         ┌──────────────────────────────────────────────┐
   DATA SUBSTRATE        │              TSC ENGINE (Engine 7)             │
   (sources, not parts)  │                                                │
                         │   ┌────────────────┐                           │
  Synthea + TSC modules ─┼──▶│ Phenome Mapper │──┐ longitudinal HPO       │
  (FHIR R4 + Clarity-    │   │   (Agent 2)    │  │ profile (foundation)   │
   shaped relational)    │   └────────────────┘  │                        │
                         │           │           ▼                        │
  BAMSurgeon BAMs ──────▶│   ┌────────────────┐  ┌────────────────┐       │
  (TSC1/TSC2, mosaic     │   │ Variant Curator│  │ Trajectory     │       │
   VAF 4-12%) →          │   │   (Agent 1)    │  │ Modeler (A3)   │       │
   Parabricks-equiv      │   └────────────────┘  └────────────────┘       │
                         │           │                    │               │
  Frontier-model notes ─▶│           │           ┌────────────────┐       │
   (~600-1000)           │           │           │ TAND Surveil.  │       │
                         │           │           │   (Agent 4)    │       │
  Frontier-model imaging │           │           └────────────────┘       │
   reports ─────────────▶│           │                    │               │
                         │           └──────┬─────────────┘               │
  ── PHASE-1 (not built) ─┤                  ▼                             │
  Epic Clarity/Caboodle  │          ┌──────────────────┐                  │
  Biobank LIMS           │          │ Therapeutics     │                  │
  Imaging-AI / DICOM     │          │ Strategist (A5)  │ + RAG + CT.gov   │
                         │          └──────────────────┘ + PubMed + FDA   │
                         │                  │                             │
                         │   ┌──────────────┴──────────────┐              │
                         │   │   TSC-Orchestrator           │              │
                         │   │   (deterministic LangGraph   │              │
                         │   │    event router; NOT an LLM) │              │
                         │   │   PostgreSQL event log +     │              │
                         │   │   projections, Redis TTL     │              │
                         │   └──────────────┬──────────────┘              │
                         └──────────────────┼─────────────────────────────┘
                                            ▼
                         ┌──────────────────────────────────────────────┐
   CLINICIAN SURFACES    │  (a) pre-visit briefing  (b) in-visit         │
   (standalone web apps, │      (1-screen mobile)       dashboard         │
    watermarked synthetic)│                          (c) async alerts     │
                         └──────────────────────────────────────────────┘
```

Three properties of this diagram matter for the engineering. First, the five CCHMC source areas sit outside the engine boundary. They feed it; they are not part of it. This is the "swap the box labels, keep the wiring" replication argument made concrete: another institution's biobank and EHR plug into the same left edge. Second, the orchestrator sits below the agents, not between them; agents do not call each other directly, they emit and consume events. Third, the surfaces pull from materialized state, never from a live agent invocation, so a slow or failed agent degrades to a "pending" panel rather than a broken page.

### 2.2 Reuse vs net-new

The engine is scaffolded from the `precision_oncology_agent` template (verified on disk: `api/` FastAPI app with a `routes/` package, `app/` Streamlit UI, `src/` with `utils/` and `workflows/`, `config/settings.py`, `tests/`, `docker-compose.yml`, `Dockerfile`, `requirements.txt`, `ruff.toml`). We inherit the HCLS AI Factory v1.3.0 substrate and add the multi-agent machinery on top.

| Layer | Status | Notes |
|---|---|---|
| LangGraph orchestration runtime | Reuse | New TSC graph definition; deterministic router config is net-new |
| Milvus RAG (BAAI/bge-large-en-v1.5 + BiomedBERT-derived clinical embeddings) | Reuse | New `tsc_literature` partition/corpus |
| Tiered Claude models (Haiku/Sonnet/Opus via API) + local Llama 3.1 70B Instruct via Ollama fallback | Reuse | New per-agent tier policy in config |
| PostgreSQL + Redis + MinIO | Reuse | New event-sourced schema, Redis key conventions, buckets |
| Audit/provenance envelope | Reuse | Extended with mosaic-flag and model-version fields |
| FastAPI / Streamlit / Pydantic scaffolding | Reuse | New routers, three surface apps |
| **Multi-agent orchestrator + shared event-sourced state** | **Net-new** | The coordination layer; the reason this is an Engine |
| **Synthetic cohort pipeline (Synthea + BAMSurgeon + frontier-model)** | **Net-new** | W1 deliverable |
| **The 5 TSC agents** | **Net-new** | Variant Curator, Phenome Mapper, Trajectory Modeler, TAND, Therapeutics |
| Epic Clarity/Caboodle + LIMS connectors | **Phase-1, not built** | Architecturally anticipated, marked at every boundary |
| Imaging-AI / DICOM pipeline | **Phase-1, not built** | Demo uses frontier-model imaging *reports*, not pixels |

### 2.3 Repository layout

The repo is `ai_agent_adds/tuberous_clerosis_engine/`. Today it is docs-only (`docs/`). The W1 scaffold lays down the following, mirroring the template and extending `src/` with the four net-new packages (`agents`, `orchestrator`, `cohort`, `rag`):

```
tuberous_clerosis_engine/
├── api/                              # FastAPI engine API (reuse template shape)
│   ├── main.py                       # app factory, lifespan, router mount
│   └── routes/
│       ├── events.py                 # POST /events (ingest), GET /events/{patient}
│       ├── surfaces.py               # GET /surfaces/{kind}/{patient}
│       ├── agents.py                 # per-agent run/status (debug + demo control)
│       └── provenance.py             # GET /provenance/{output_id}
├── app/                              # three Streamlit surfaces (standalone web apps)
│   ├── briefing_app.py               # (a) pre-visit, 1-screen mobile-readable
│   ├── dashboard_app.py              # (b) in-visit 4-quadrant
│   └── alerts_app.py                 # (c) async alert surface
├── src/
│   ├── agents/
│   │   ├── variant_curator/          # Agent 1: BAM→VCF→annotate→ACMG synthesis
│   │   ├── phenome_mapper/           # Agent 2: HPO extraction (runs first)
│   │   ├── trajectory_modeler/       # Agent 3: classical stats, NOT LLM
│   │   ├── tand_surveillance/        # Agent 4: Marshall-Hagedorn taxonomy
│   │   └── therapeutics_strategist/  # Agent 5: 6-section brief, Opus-class
│   ├── orchestrator/
│   │   ├── graph.py                  # LangGraph deterministic router
│   │   ├── events.py                 # 13 event-type definitions (enum + schemas)
│   │   ├── state.py                  # append-only writer + projection materializer
│   │   └── policies.py               # enrollment ordering, failure handling
│   ├── cohort/
│   │   ├── synthea_tsc/              # Synthea modules + post-processing → Clarity shape
│   │   ├── bamsurgeon/               # variant spike-in, mosaic VAF control
│   │   ├── notes_gen/                # frontier-model clinical notes
│   │   └── imaging_reports/          # frontier-model imaging reports (text)
│   ├── rag/                          # Milvus client, tsc corpus loaders, retrievers
│   └── utils/                        # provenance envelope, model router, hashing, audit
├── config/
│   ├── settings.py                   # pydantic-settings; ports, model tiers, URIs
│   ├── model_policy.yaml             # per-agent Haiku/Sonnet/Opus assignment
│   └── demo.yaml                     # demo cohort config, featured patients A/B/C
├── tests/                            # pytest; unit + eval harness (vs synthetic ground truth)
│   ├── eval/                         # recall/precision/latency targets per agent
│   └── unit/
├── scripts/
│   ├── regen_cohort.py               # full deterministic regen (~12 hr)
│   ├── load_rag.py                   # build tsc_literature Milvus partition
│   └── dry_run_demo.py               # 3-act demo rehearsal driver
├── docker-compose.yml                # Spark-local stack
├── docker-compose.runpod.yml         # burst overlay (Parabricks, parallel gen)
├── Dockerfile
├── requirements.txt
└── ruff.toml
```

The agent packages each follow the same internal contract: a `runner.py` exposing `run(patient_id, context) -> AgentOutput`, a `prompts/` directory with versioned templates, a `schema.py` of Pydantic models, and a `tier.py` reading from `model_policy.yaml`. This uniformity is what lets the orchestrator treat all five agents as interchangeable event consumers despite their wildly different internals (classical statistics for Trajectory, Opus-class synthesis for Therapeutics).

### 2.4 Services and ports

The stack runs as Docker Compose services. Engine ports sit in the 856x band to stay clear of the core platform and the existing intelligence agents documented in the platform port map.

| Service | Port | Tier | Notes |
|---|---|---|---|
| Engine API (FastAPI) | 8560 | net-new | event ingest, surface assembly, provenance |
| Pre-visit briefing (Streamlit) | 8561 | net-new | surface (a) |
| In-visit dashboard (Streamlit) | 8562 | net-new | surface (b) |
| Async alert surface (Streamlit) | 8563 | net-new | surface (c) |
| PostgreSQL | 5432 | reuse | event log + projections |
| Redis | 6379 | reuse | ephemeral state, TTL keys |
| MinIO | 9000 / 9001 | reuse | object store + console |
| Milvus | 19530 | reuse | `tsc_literature` collection |
| Ollama (Llama 3.1 70B fallback) | 11434 | reuse | local-LLM fallback path |

The three surfaces are deliberately separate Streamlit apps rather than tabs in one app. They have different latency budgets, different refresh cadences, and different audiences (a clinician glances at the briefing on a phone; the dashboard is a desk-during-visit tool). Splitting them also keeps the alert-discipline logic (recalibrate if more than roughly three alerts per clinician per week) isolated in one service.

### 2.5 Data model

The state is event-sourced. Every agent output, every ingest, and every orchestrator decision is an append-only event; the read side is a set of materialized projections rebuilt deterministically from the log. This is the property that makes the whole engine auditable: you can replay the log and reproduce any surface exactly, and provenance is not bolted on after the fact, it is the substrate.

#### 2.5.1 PostgreSQL — append-only event log

```sql
CREATE TABLE engine_events (
    event_id      BIGSERIAL PRIMARY KEY,
    patient_id    TEXT NOT NULL,
    event_type    TEXT NOT NULL,          -- one of 13 enumerated types
    payload       JSONB NOT NULL,         -- type-specific, schema-validated
    provenance    JSONB NOT NULL,         -- model id/version, prompt template version,
                                          -- RAG source URIs, input hash, latency_ms
    created_at    TIMESTAMPTZ NOT NULL DEFAULT now(),
    parent_event  BIGINT REFERENCES engine_events(event_id)  -- causal lineage
);
CREATE INDEX ix_events_patient ON engine_events(patient_id, event_id);
CREATE INDEX ix_events_type    ON engine_events(event_type, created_at);
```

The thirteen event types span the lifecycle: `cohort_loaded`, `patient_enrolled`, `genomic_substrate_ready`, `variant_curated`, `phenome_mapped`, `trajectory_forecast`, `tand_surveyed`, `therapeutics_briefed`, `surface_requested`, `surface_assembled`, `alert_emitted`, `agent_failed`, `provenance_logged`. The router (Section 2.6) keys all of its decisions off this enum.

#### 2.5.2 PostgreSQL — materialized projections

Projections are current-state views the surfaces read directly. They never carry information the log does not; they are a denormalized cache for read latency.

```sql
CREATE TABLE proj_patient_current (
    patient_id      TEXT PRIMARY KEY,
    hpo_profile     JSONB,        -- latest from phenome_mapped, with evidence spans
    variant_interp  JSONB,        -- latest variant_curated incl. mosaic flag
    trajectory      JSONB,        -- forecasts w/ 50/90% prediction intervals
    tand_briefing   JSONB,        -- pre-visit TAND material (never alerts)
    therapeutics    JSONB,        -- 6-section options brief
    staleness       JSONB,        -- per-section {source_event_id, ts, status}
    last_event_id   BIGINT NOT NULL
);
```

`staleness` is the conservative-failure mechanism in data form. If `tand_surveillance` failed for a patient, `tand_briefing` is left at its prior value and `staleness.tand = {status: "pending", ...}`, so the surface renders "pending" rather than silently dropping a panel.

#### 2.5.3 Redis — ephemeral state

Redis holds nothing durable. Keys are namespaced and TTL'd:

- `lock:agent:{patient_id}:{agent}` — run mutex, short TTL, prevents duplicate dispatch
- `inflight:{patient_id}` — set of agents currently running, for surface "computing" hints
- `surface:cache:{kind}:{patient_id}` — assembled surface payload, 60 s TTL
- `ratelimit:alerts:{clinician_id}` — sliding-window counter feeding alert discipline

#### 2.5.4 Milvus — collections

One collection, `tsc_literature`, partitioned by source (`pubmed_pmc`, `ctgov_snapshot`, `fda_actions`, `itsc_guidelines`). Embeddings: BAAI/bge-large-en-v1.5 for general literature and the BiomedBERT-derived clinical embedding model for clinically dense passages, stored as separate vector fields with metadata (`source_uri`, `pub_year`, `section`, `embedding_model`). The Therapeutics Strategist is the primary consumer; every retrieved chunk's `source_uri` flows into that agent's source attribution.

#### 2.5.5 MinIO — buckets

- `tsc-cohort` — version-controlled synthetic cohort artifacts (BAMs, VCFs, FHIR bundles, notes, imaging reports), content-addressed
- `tsc-reports` — generated draft molecular-genetics reports and option briefs (watermarked synthetic)
- `tsc-provenance` — large provenance blobs (full prompt renders, RAG context bundles) referenced by `event_id`

### 2.6 Orchestration in the deployment

The orchestrator is a deterministic LangGraph router, not an LLM, and not a service the surfaces talk to directly. Its three operative patterns map cleanly onto the data model. Dependency-ordered enrollment: on `patient_enrolled` it dispatches Phenome Mapper first and gates the others on `phenome_mapped`. Incremental-update minimization: an incoming event is diffed against the projection, and only agents whose inputs actually changed are re-dispatched, which keeps cohort-scale reruns cheap. Demand-driven surface assembly: a surface is built only on `surface_requested`, reading projections, with the 60 s Redis cache absorbing repeated reads during a visit. Conservative failure handling is the `agent_failed` event writing `staleness` and never a silent gap.

### 2.7 Deployment topology

The default and primary target is a single DGX Spark: GB10 Grace Blackwell, roughly 1,000 TOPS, 128 GB unified LPDDR5x, 4 TB NVMe, DGX OS. The full Compose stack, all three surfaces, PostgreSQL, Redis, MinIO, Milvus, and the Ollama fallback fit inside the unified-memory envelope for the 50-patient cohort. The agents' heavy reasoning is Claude via API (Haiku/Sonnet/Opus per `model_policy.yaml`), so the Spark's GPU is not the bottleneck for routine operation; it carries the local Llama 3.1 70B fallback, embedding generation, and the deterministic statistical models in the Trajectory Modeler.

RunPod is a burst overlay (`docker-compose.runpod.yml`), used only where a job is embarrassingly parallel or GPU-bound beyond the Spark's comfort:

| Workload | Where | Why |
|---|---|---|
| Routine agent operation, surfaces, demo | DGX Spark | Fits in 128 GB; API-driven LLM calls |
| Parabricks-equivalent variant calling (cohort build) | RunPod burst | GPU-accelerated alignment/calling at cohort scale |
| Parallel synthetic-cohort generation (~12 hr regen) | RunPod burst | Fan out Synthea/BAMSurgeon/notes across GPUs |
| Heavier local-LLM inference if API budget-capped | RunPod burst | Larger local models than the Spark holds comfortably |

The split is deliberate and is itself part of the cost argument made in the demo's third act: the steady state runs on a $4,699 box, and you rent GPUs by the hour only for the cohort-build and regeneration spikes. Constraints are features. RunPod is an unlock for throughput, not a dependency for the demo to function. If RunPod is unavailable, cohort regeneration runs slower on the Spark; the engine still serves.

Everything is Apache 2.0, open, and non-commercial in intent. The surfaces are standalone web apps, not Epic integrations, and they are persistently watermarked as synthetic. Not FDA-cleared; the SaMD posture is undetermined and is institutional work; IRB is not required for the synthetic demo and is required before any real-data Phase-1.

With the architecture, repo, data model, and deployment fixed, Section 3 turns to the functional requirements, beginning with the data substrate and the Variant Curator that recovers the low-VAF mosaic variants standard blood testing misses.


## 3. Functional Requirements — Cohort Generation, Variant Curator, Phenome Mapper

This section specifies the first three buildable components of the TSC Intelligence Engine to implementation grade: the synthetic cohort pipeline that produces the engine's only data, the TSC-Variant Curator (Agent 1), and the TSC-Phenome Mapper (Agent 2). Section §2 established the engine's architecture and the source-area mapping onto CCHMC. Here the contract tightens into numbered functional requirements (FRs), each with acceptance criteria, I/O schemas, model-tier assignments, deterministic-tool specs, endpoints, provenance hooks, and per-component evaluation targets measured against synthetic ground truth.

Two framing constraints govern everything below. First, the demo runs on **synthetic data only** — generated on Adam's DGX Spark with RunPod GPU burst, version-controlled, regenerable in ~12 hours. No Epic Clarity, Caboodle, or biobank LIMS feed is built; those are institutional Phase-1 work and are marked as such wherever the architecture touches them. Second, every numeric target in this section is a **demo eval target against synthetic ground truth, not a clinical-validation claim**. The cohort is the answer key; passing the eval means the engine reproduces what we deliberately injected, nothing more.

### 3.1 Synthetic Cohort Generation (FR-CG-*)

The cohort is generated once, committed, and treated as a fixed artifact for the duration of the build. Determinism is the load-bearing property: a fixed master seed (`COHORT_SEED=2026`) plus pinned tool versions must reproduce byte-identical genomic substrate and stably-hashed clinical text. The pipeline is four layers, each a discrete buildable task with its own outputs and gate.

The 50-patient composition is fixed by CANON and is the ground-truth manifest every downstream eval reads:

| Stratum | Count | Genotype | Notes |
|---|---|---|---|
| TSC2 germline | 30 (60%) | pathogenic/LP TSC2, VAF ~50% | bulk of cohort |
| TSC1 germline | 12 (24%) | pathogenic/LP TSC1, VAF ~50% | |
| TSC2 mosaic | 5 (10%) | TSC2, VAF 4–12% | mosaic-recovery targets |
| TSC1 mosaic | 2 (4%) | TSC1, VAF 4–12% | mosaic-recovery targets |
| NMI | 1 (2%) | tissue-only low-VAF (Patient A, 8.3% TSC2) | no-mutation-identified on blood |
| **Mosaic total** | **7** | | Variant Curator recall denominator |

#### FR-CG-1 — Layer 1: Synthea clinical skeleton

**Description.** Generate demographic and clinical skeletons for all 50 patients using Synthea (MIT-licensed) extended with custom TSC disease modules. Output is FHIR R4 bundles plus a Clarity-shaped relational projection (synthetic stand-in for the institutional schema, explicitly not the real Clarity feed). TSC modules encode the natural history from CANON disease facts: epilepsy onset (infantile spasms in infancy, ~85% lifetime, ~2/3 refractory), SEGA emergence near the foramen of Monro, renal AML (~80%, bleeding risk above ~4 cm), TAND features (~90% prevalence), and age/sex-appropriate surveillance encounters per ITSC 2021 cadence.

**Acceptance criteria.**
- 50 FHIR R4 bundles validate against the R4 profile; 0 validation errors.
- Each bundle carries a longitudinal encounter history of ≥ 18 months (the trajectory-forecast horizon) with TSC-appropriate problem lists, medications, and surveillance imaging encounters.
- The three featured patients (A: 4yo F; B: 12yo M; C: 18yo F) match their CANON phenotypes exactly.
- Re-running with `COHORT_SEED=2026` produces identical patient IDs, demographics, and encounter timelines.

**Implementation notes.** The TSC modules are authored as Synthea GMF (Generic Module Framework) JSON state machines — one module per organ axis (neuro, renal, pulmonary, cardiac, dermatologic, TAND) plus a master TSC-onset module that gates the others by genotype stratum. Module transitions are seeded from the per-stratum natural-history priors in CANON, so a TSC2-mosaic patient and a TSC2-germline patient traverse the same module graph but with stratum-conditioned probabilities (mosaic patients trend toward milder, later-onset, more variable phenotypes — the clinical reason they reach the NMI pathway). The Clarity-shaped projection is produced by a deterministic flattener that maps FHIR resources onto a small relational schema (`patient`, `encounter`, `diagnosis`, `medication`, `lab_result`, `imaging_order`) whose column names deliberately echo the institutional schema so that the Phenome Mapper's structured-normalization code (FR-PM-1) is written once and re-points at the real Caboodle views during Phase-1 without rework. That column-name parity is the only coupling to the real institution; no real data flows.

**Outputs.** `/cohort/layer1/fhir/{patient_id}.json`, `/cohort/layer1/clarity/*.parquet`, and a `cohort_manifest.json` (the ground-truth answer key).

```json
// cohort_manifest.json (per-patient entry, abbreviated)
{
  "patient_id": "TSC-0007",
  "stratum": "TSC2_mosaic",
  "sex": "F", "dob": "2022-03-14",
  "ground_truth_variant": {
    "gene": "TSC2", "hgvs_c": "c.4255C>T", "hgvs_p": "p.Arg1419Ter",
    "vaf": 0.083, "zygosity": "mosaic", "acmg_expected": "Likely Pathogenic",
    "acmg_criteria_expected": ["PVS1","PM2","PP4"]
  },
  "ground_truth_hpo": ["HP:0002359","HP:0009732","HP:0001250"],
  "ground_truth_tand_signals": [{"cluster":"behavioral","span_id":"note-014:p3"}]
}
```

#### FR-CG-2 — Layer 2: Genomic substrate (BAMSurgeon → variant calling → VCF)

**Description.** Build realistic per-patient BAMs by spiking TSC1/TSC2 variants into NA12878-derived alignments with BAMSurgeon, then call variants with a Parabricks-equivalent pipeline to produce VCFs the Variant Curator will consume. Germline variants are spiked at VAF ~50%; mosaic and NMI variants at VAF 4–12% (Patient A pinned at 8.3%). The pipeline **starts at BAM** — no raw FASTQ is produced or stored.

This is the GPU-heavy layer and the canonical RunPod-burst task: parallel BAMSurgeon spiking and Parabricks-equivalent calling across 50 samples on rented GPUs, results staged back to the Spark.

**Acceptance criteria.**
- All 50 VCFs produced; each spiked variant is recoverable in its source BAM at the intended VAF ± 1.5 percentage points.
- Mosaic variants show realistic read-level signal: low-VAF alt reads on both strands, mapping/base qualities consistent with true variants (not artifacts), so the Variant Curator's mosaic-vs-artifact discrimination is genuinely exercised.
- VCFs carry standard FORMAT fields (`AD`, `DP`, `AF`, `SB` for strand bias).
- Deterministic given `COHORT_SEED` and pinned BAMSurgeon/caller versions.

**Deterministic-tool spec.** BAMSurgeon (pinned commit) for spiking; calling mirrors the production Variant Curator stack (BWA-MEM/GATK HaplotypeCaller + Mutect2, mosaic-aware, low-VAF ≥ 5%) so demo and production paths agree. Outputs `/cohort/layer2/bam/{patient_id}.bam(.bai)`, `/cohort/layer2/vcf/{patient_id}.vcf.gz`.

**Why the substrate must be hard.** The single most important property of this layer is that the mosaic and NMI variants are genuinely difficult — not trivially separable from background error. If BAMSurgeon spiked clean, high-quality alt reads at 8% VAF the Variant Curator's recovery would be unimpressive and the demo would be hollow. We therefore tune BAMSurgeon to place mosaic alt reads with realistic base-quality distributions, balanced strand representation, and depth comparable to the surrounding germline regions (~300–500×), so that distinguishing the 8.3% TSC2 frameshift in Patient A from a sequencing artifact is a real discrimination task. The Layer-2 acceptance check explicitly seeds a small number of deliberate low-VAF artifacts (strand-biased, low-MAPQ) that are **not** in the manifest, so that "0 false-positive Pathogenic" in the Variant Curator eval is a meaningful result rather than a guarantee handed over by an easy substrate.

#### FR-CG-3 — Layer 3: Clinical notes

**Description.** Generate ~600–1,000 longitudinal clinical notes across the cohort using a frontier model (Opus-class), conditioned on each patient's Layer-1 timeline and seeded from published TSC note templates (neurology, nephrology, genetics, developmental-behavioral). Every TAND ground-truth signal in the manifest is **embedded at a specific note span** so the TAND agent (§4) has a recoverable target; uncertainty/discourse markers are placed per the Marshall-Hagedorn taxonomy. All notes carry a persistent synthetic watermark in text and metadata.

**Acceptance criteria.**
- 600 ≤ note count ≤ 1,000; every patient has ≥ 6 notes spanning the ≥ 18-month window.
- Each note carries metadata: `patient_id`, `encounter_date`, `author_role`, `note_type`, `synthetic: true`, `watermark_id`, and a generation provenance block.
- Every `ground_truth_hpo` and `ground_truth_tand_signals` entry traces to ≥ 1 concrete note span (`note_id:paragraph`), enabling recall scoring.
- A clinician-sampled review (≥ 10% of notes) passes plausibility with no edits that would change extractable phenotypes.

**Model tier.** Opus for generation (clinical realism is the gate); a deterministic post-processor injects watermarks and validates that each manifest signal is present in text. Outputs `/cohort/layer3/notes/{patient_id}/{note_id}.json`.

**Note on the TAND embedding discipline.** The TAND signals are the subtlest ground truth in the cohort and the hardest to validate, because the whole point of the TAND agent (§4) is to surface features that are *present but under-recognized* in documentation. The note generator therefore embeds TAND signals as the Marshall-Hagedorn taxonomy predicts they appear in real notes: hedged ("parents report some behavioral concerns at home, will monitor"), deferred ("recommend developmental-behavioral referral when scheduling allows"), third-party-attributed ("teacher notes attention difficulties"), conditional, and follow-up-without-formalization. Each such embedded marker is recorded in the manifest with its `span_id` and its taxonomy category, so that §4's eval can measure both detection and correct cluster assignment. Crucially, the generator also writes notes for patients with **no** TAND ground truth and notes containing benign behavioral language that should **not** trip the agent, giving the TAND eval a real specificity denominator.

#### FR-CG-4 — Layer 4: Imaging reports

**Description.** Generate frontier-model **textual** imaging reports (brain MRI, renal ultrasound, echocardiogram, ophthalmology) — explicitly **no DICOM, no pixel data**. The longitudinal SEGA series for Patient B is the marquee artifact: 0.8 → 1.1 → 1.3 cm at the foramen of Monro across the window (~2–4 mm/yr), feeding the Trajectory Modeler's threshold-crossing eval. Renal AML measurements track toward the ~4 cm bleeding-risk threshold for Patient C.

**Acceptance criteria.**
- Each imaging encounter from Layer 1 has a matching structured report with measured findings (lesion sizes in mm, locations, interval-change language).
- Patient B's SEGA series reproduces the CANON trajectory exactly (0.8/1.1/1.3 cm).
- Reports carry the same synthetic watermark and provenance metadata as notes.

**Outputs.** `/cohort/layer4/imaging/{patient_id}/{report_id}.json`.

#### FR-CG-5 — Cohort gate and regeneration

**Description.** A single orchestration target (`make cohort`) runs Layers 1–4 in order and emits a `cohort_build_report.json` summarizing counts, validation results, and a diff against the committed manifest. Full regeneration completes in ≤ ~12 hours wall-clock (Spark + RunPod burst for Layer 2).

**Acceptance criteria.** Clean run from empty workspace reproduces the committed cohort with 0 manifest diffs (genomic substrate byte-identical; clinical text stable-hashed). Build report flags any drift.

**Cohort eval targets.** Composition exactly 30/12/5/2/1; 7 mosaic targets present and recoverable; all featured-patient facts pinned; full regen ≤ 12 h; deterministic.

### 3.2 TSC-Variant Curator — Agent 1 (FR-VC-*)

The Variant Curator turns a per-patient BAM into a draft, provenance-complete molecular-genetics interpretation, with its defining capability being recovery of **low-VAF somatic mosaic variants (VAF ≥ 5%)** that standard blood testing misses — the ~10–15% NMI cohort (Tyburczy 2015; Giannikou 2016; Lim 2017). The pipeline is a deliberate alternation of deterministic tools and tiered LLM reasoning, and **no output is autonomous**: every result is a draft for board-certified molecular-geneticist sign-off.

**Stepwise model-tier assignment.**

| Step | Mechanism | Tier | Rationale |
|---|---|---|---|
| 1. Variant calling | BWA-MEM/GATK HaplotypeCaller + Mutect2, mosaic-aware, ≥ 5% VAF | Deterministic | reproducibility, no LLM in the numeric path |
| 2. Annotation | snpEff/VEP + ClinVar/gnomAD v4/LOVD-TSC/dbSNP joins | Deterministic | |
| 3. Evidence aggregation | structure per-variant evidence dossier | Sonnet | synthesis over heterogeneous annotations |
| 4. ACMG-AMP classification | criteria assignment + class synthesis | Opus, validated vs combinatorial rules | clinical stakes; cross-checked deterministically |
| 5. Mosaic adjudication | VAF / read depth / strand / artifact assessment | Sonnet (numbers deterministic; narrative Sonnet) | mosaic-vs-artifact call |
| 6. Draft report | AI-labeled molecular-genetics narrative | Opus | reviewer-facing prose |

#### FR-VC-1 — Deterministic variant calling and annotation

**Description.** Ingest a patient BAM; call variants with the mosaic-aware stack (germline HaplotypeCaller + somatic/low-VAF Mutect2), retaining calls with VAF ≥ 5%; annotate with snpEff/VEP and join ClinVar, gnomAD v4, LOVD-TSC, and dbSNP.

**Acceptance criteria.** Recovers all spiked TSC1/TSC2 variants at VAF ≥ 5% including all 7 mosaic targets; emits per-variant `vaf`, `depth`, `alt_reads`, `strand_balance`, `mapping_quality`; produces 0 false-positive Pathogenic-eligible calls outside ground truth. Deterministic given fixed inputs.

**Scope discipline.** Calling and annotation are scoped to the TSC1/TSC2 loci plus the small contiguous-deletion region (TSC2/PKD1) that CANON flags for the ~1–2% severe-PKD subset. The Curator is a TSC-focused interpreter, not a general clinical-exome pipeline; restricting the locus set keeps per-case latency under the 5-minute target and removes the incidental-findings surface area that would otherwise demand its own governance discussion. The annotation join is performed against locally-staged, version-pinned snapshots of ClinVar (`2026-04`), gnomAD v4, LOVD-TSC, and dbSNP so that the deterministic path has no network dependency and reproduces exactly months later.

#### FR-VC-2 — Evidence aggregation (Sonnet)

**Description.** For each candidate TSC1/TSC2 variant, assemble a structured evidence dossier: population frequency, ClinVar/LOVD assertions, in-silico predictors, gene-level constraint, and functional/literature context retrieved from the Milvus TSC corpus partition (BAAI/bge-large-en-v1.5 + BiomedBERT-derived clinical embeddings).

**Prompt-design notes.** Sonnet receives only deterministic, source-attributed fields and retrieved snippets with URIs; it is instructed to **aggregate and cite, not to assign ACMG criteria**. Each evidence item must carry its source. Temperature 0; prompt template versioned (`vc_evidence_v1`).

```json
// VC evidence dossier (per variant)
{
  "variant": {"gene":"TSC2","hgvs_c":"c.4255C>T","hgvs_p":"p.Arg1419Ter",
              "vaf":0.083,"depth":412,"alt_reads":34,"strand_balance":0.47},
  "evidence": [
    {"type":"population","value":"absent in gnomAD v4","source":"gnomad:v4","supports":"PM2"},
    {"type":"consequence","value":"nonsense, NMD-predicted","source":"vep","supports":"PVS1"},
    {"type":"clinvar","value":"no submission","source":"clinvar:2026-04"},
    {"type":"literature","value":"recurrent truncating TSC2 variants pathogenic",
     "source":"pmc:PMC1234567","supports":"PP4"}
  ]
}
```

#### FR-VC-3 — ACMG-AMP classification synthesis (Opus, validated)

**Description.** Opus assigns ACMG-AMP criteria and synthesizes a class (Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign). A **deterministic combinatorial-rule engine** independently computes the class from the assigned criteria; the two must agree or the case is flagged for mandatory reviewer attention.

**Acceptance criteria.** Class matches manifest `acmg_expected` for all ground-truth variants (e.g., Patient A → Likely Pathogenic via PVS1+PM2+PP4); 0 false-positive Pathogenic calls; Opus-assigned criteria and the rule engine agree on 100% of ground-truth variants (disagreements are surfaced, never silently resolved). Prompt template `vc_acmg_v1`, temperature 0, criteria definitions pinned to the ACMG-AMP reference.

**The two-path design is the credibility argument.** This is the most clinically consequential step in the engine and the one a skeptical molecular geneticist will probe hardest, so the design deliberately refuses to let an LLM be the sole arbiter of a classification. Opus reasons about *which criteria apply* (it is good at marshalling PVS1 loss-of-function logic, PM2 rarity, and PP4 phenotype-fit from the dossier), but the *class that follows from a criteria set* is a fixed combinatorial rule, and a small deterministic engine computes it independently. When the two paths agree — the expected case — the reviewer sees a classification backed by both an explainable narrative and a mechanical rule check. When they disagree, the case is flagged `rule_engine_agree:false` and routed to the front of the reviewer queue with both rationales shown side by side. The engine never picks a winner. This is what lets us claim the Curator is decision-support that *shows its work*, not a black box emitting a verdict.

#### FR-VC-4 — Mosaic adjudication and ddPCR recommendation

**Description.** For any variant with VAF below the germline expectation (~50%), emit a **mosaic flag** with the supporting numbers (VAF, alt-read count, strand balance, artifact assessment) and an orthogonal-validation recommendation (ddPCR). Numbers are deterministic; Sonnet writes the adjudication narrative only.

**Acceptance criteria.** All 7 mosaic targets flagged with correct VAF (± 1.5 pp), correctly distinguished from artifact, each carrying a ddPCR recommendation; 0 germline (~50% VAF) variants mis-flagged as mosaic.

#### FR-VC-5 — Draft molecular-genetics report and human gate

**Description.** Opus assembles an **AI-labeled, ClinVar-spec** draft interpretation: variant, classification with criteria, mosaic status, validation recommendation, and full provenance. The report is explicitly a **draft for review**; the UI requires molecular-geneticist sign-off before any downstream surface treats it as final.

```json
// VC output contract
{
  "patient_id":"TSC-0007",
  "status":"DRAFT_FOR_REVIEW",
  "variant":{"gene":"TSC2","hgvs_c":"c.4255C>T","hgvs_p":"p.Arg1419Ter"},
  "classification":{"class":"Likely Pathogenic","criteria":["PVS1","PM2","PP4"],
                    "rule_engine_agree":true},
  "mosaic":{"is_mosaic":true,"vaf":0.083,"alt_reads":34,"strand_balance":0.47,
            "artifact_assessment":"unlikely_artifact","validation_rec":"ddPCR"},
  "ai_labeled":true,
  "provenance":{"model_ids":{"evidence":"claude-sonnet-4-5","acmg":"claude-opus-4-8"},
    "prompt_versions":{"evidence":"vc_evidence_v1","acmg":"vc_acmg_v1"},
    "rag_sources":["gnomad:v4","clinvar:2026-04","pmc:PMC1234567"],
    "input_bam_sha256":"…","calling_pipeline":"parabricks-equiv@<ver>","latency_ms":211342}
}
```

**Acceptance criteria.** Output validates against the contract; `ai_labeled:true` and `status:"DRAFT_FOR_REVIEW"` always set; provenance block complete (model ids/versions, prompt versions, RAG URIs, input hash, latency); no surface consumes the report before sign-off.

**Endpoints.**

```
POST /api/v1/variant-curator/run        body: {patient_id} -> {job_id}
GET  /api/v1/variant-curator/result/{patient_id}   -> VC output contract
POST /api/v1/variant-curator/sign-off   body: {patient_id, reviewer_id, decision} 
GET  /api/v1/variant-curator/provenance/{patient_id}
```

**Variant Curator eval targets.** Recover all 7 mosaic variants at VAF ≥ 5%; correct ACMG class on all ground-truth variants; **0 false-positive Pathogenic**; rule-engine/Opus agreement 100% on ground truth; **< 5 min per case**; provenance complete on every output. Patient A (NMI, 8.3% VAF) is the live Act-One target.

### 3.3 TSC-Phenome Mapper — Agent 2 (FR-PM-*)

The Phenome Mapper is the foundation the other four agents build on: it converts synthetic Epic Clarity-shaped structured data and clinical notes into a **longitudinal, time-anchored HPO-coded phenotype profile with evidence spans**. It also produces a discordance log and an ITSC-surveillance-gap report. (The real Clarity/Caboodle feed is institutional Phase-1 work; the demo reads only the synthetic Clarity-shaped projection and Layer-3 notes.)

**Stepwise model-tier assignment.**

| Step | Mechanism | Tier | Rationale |
|---|---|---|---|
| 1. Structured normalization | ICD-10 / lab → HPO mapping | Haiku | high-volume, low-ambiguity |
| 2. Note phenotype extraction | per-note HPO extraction with spans | Sonnet | discourse-level reasoning |
| 3. Conflict resolution | rare cross-source contradictions | Opus | high-ambiguity adjudication |
| 4. Surveillance-gap report | ITSC 2021 cadence check | Deterministic | rule-based |

#### FR-PM-1 — Structured-data normalization (Haiku)

**Description.** Map ICD-10 diagnoses, problem-list entries, and lab results from the synthetic Clarity-shaped tables to HPO terms (via SNOMED-CT crosswalks), each as a time-anchored phenotype assertion with its source row as the evidence span.

**Acceptance criteria.** All ICD-10/lab rows with a defined HPO mapping are normalized; each assertion carries `hpo_id`, `onset_date`/`assertion_date`, and `source: {table, row_id}`; mappings are deterministic given a pinned crosswalk version. Haiku, temperature 0, template `pm_normalize_v1`.

#### FR-PM-2 — Note-based phenotype extraction (Sonnet)

**Description.** For each Layer-3 note, extract HPO-coded phenotypes with **character-offset evidence spans**, polarity (present/absent/uncertain), and temporality (current/historical). Sonnet runs per note (bounded context, parallelizable). Retrieval from the Milvus TSC partition supplies HPO definitions and TSC phenotype context.

**Prompt-design notes.** Sonnet must (a) return a span for every extracted term, (b) mark negation and uncertainty explicitly, and (c) never invent HPO IDs — terms are validated against the pinned HPO release post-hoc; unvalidated terms are dropped and logged. Template `pm_extract_v1`, temperature 0. The prompt does not ask the model to free-recall HPO codes (a known hallucination trap); instead it asks for the surface phenotype phrase and a span, and a deterministic post-step resolves phrases to HPO IDs via a pinned synonym index, leaving the LLM to do the linguistically-hard work (finding, negating, and temporally situating mentions) while the brittle code-assignment step stays mechanical and auditable. Span fidelity is enforced: every returned span must be a verbatim substring of the source note at the claimed offsets, or the extraction is rejected and re-requested once before being logged as a miss.

```json
// PM per-note extraction
{
  "note_id":"TSC-0012:note-014",
  "phenotypes":[
    {"hpo_id":"HP:0001250","label":"Seizure","polarity":"present",
     "temporality":"current","span":{"para":3,"start":118,"end":141},
     "evidence_text":"recurrent focal seizures"},
    {"hpo_id":"HP:0009732","label":"Hypomelanotic macule","polarity":"present",
     "temporality":"historical","span":{"para":1,"start":40,"end":72}}
  ]
}
```

#### FR-PM-3 — Longitudinal aggregation and discordance log

**Description.** Merge structured (FR-PM-1) and note-derived (FR-PM-2) phenotypes into one **longitudinal HPO profile** keyed by HPO ID with first-onset, evidence-span list, and source provenance. Where sources disagree (e.g., note asserts present, structured asserts absent), record a **discordance entry**; only genuinely rare/ambiguous conflicts are escalated to Opus (template `pm_conflict_v1`) for an adjudication note — never an autonomous resolution.

```json
// PM patient profile (abbreviated)
{
  "patient_id":"TSC-0012",
  "hpo_profile":[
    {"hpo_id":"HP:0001250","first_onset":"2014-06-01",
     "evidence":[{"note_id":"TSC-0012:note-014","span":{"para":3}},
                 {"source":{"table":"diagnosis","row_id":88}}]}
  ],
  "discordance_log":[
    {"hpo_id":"HP:0002353","note_says":"present","structured_says":"absent",
     "resolution":"escalated","adjudication_model":"claude-opus-4-8"}
  ],
  "provenance":{"models":{"normalize":"claude-haiku-4-5","extract":"claude-sonnet-4-5"},
    "hpo_release":"2026-03","crosswalk_version":"…","latency_ms":2841190}
}
```

**Acceptance criteria.** Profile validates against the contract; every HPO assertion carries ≥ 1 evidence span; discordances logged (not silently dropped); escalations recorded with the adjudicating model id.

#### FR-PM-4 — ITSC surveillance-gap report (deterministic)

**Description.** Compare each patient's encounter/imaging history against the ITSC 2021 consensus surveillance cadence (e.g., brain MRI intervals by age, renal imaging, neuropsych/TAND screening) and emit a deterministic gap report: which surveillance items are due or overdue given the patient's age and last-performed dates.

**Acceptance criteria.** Gap report deterministically reproduces overdue items for the featured patients given their Layer-1 timelines; no LLM in the cadence logic; each gap cites the ITSC 2021 rule it derives from.

**Endpoints.**

```
POST /api/v1/phenome-mapper/run         body: {patient_id} -> {job_id}
POST /api/v1/phenome-mapper/run-cohort  -> {job_id}   # full 50-patient pass
GET  /api/v1/phenome-mapper/profile/{patient_id}
GET  /api/v1/phenome-mapper/surveillance-gaps/{patient_id}
GET  /api/v1/phenome-mapper/provenance/{patient_id}
```

**Phenome Mapper eval targets.** Phenotype **recall ≥ 90%, precision ≥ 85%** against the manifest's `ground_truth_hpo` and embedded note spans; every assertion span-grounded; full **cohort pass < 1 hour**; provenance complete. Recall/precision are computed per patient against the answer key and aggregated cohort-wide; spans must overlap the embedded ground-truth location to count as a true positive.

The scoring is intentionally strict on grounding: a correct HPO ID with no overlapping span does **not** count as a true positive, because an ungrounded phenotype is exactly the kind of unverifiable assertion a clinical reviewer cannot trust. Negated and historical mentions are scored against their manifest polarity, so the agent is rewarded for correctly extracting an *absent* finding (a surveillance-relevant negative) and penalized for asserting it as present. The < 1-hour cohort budget is met by running the 50 per-patient Sonnet passes concurrently with bounded fan-out; the Haiku normalization and the deterministic surveillance-gap pass are negligible by comparison, so the wall-clock is dominated by — and tuned against — note-extraction concurrency on the Spark.

### 3.4 Cross-component provenance and orchestration hooks

All three components write to the engine's shared, append-only event-sourced state. Each component emits a typed completion event (`cohort.built`, `variant_curator.completed`, `phenome_mapper.completed`) consumed by the deterministic TSC-Orchestrator (§5); the Phenome Mapper's completion is the dependency the orchestrator waits on before enrolling downstream agents, per the dependency-ordered enrollment pattern. Every output object carries the standard provenance block — model ids/versions, prompt-template versions, RAG source URIs, input hash, latency — queryable through each component's `/provenance/{patient_id}` endpoint.

With the cohort fixed and the variant and phenotype substrates established and span-grounded, §4 specifies the second tier of agents — the Trajectory Modeler, TAND Surveillance Agent, and Therapeutics Strategist — which consume these outputs rather than the raw cohort.


## 4. Functional Requirements — Trajectory Modeler, TAND, Therapeutics, Orchestrator, Surfaces

Section 3 specified the two upstream agents the rest of the engine builds on: the TSC-Variant Curator (the molecular substrate, including low-VAF mosaic recovery for the NMI cohort) and the TSC-Phenome Mapper (the longitudinal HPO foundation). This section specifies everything that consumes those outputs: the TSC-Trajectory Modeler (Agent 3), the TAND Surveillance Agent (Agent 4), the TSC-Therapeutics Strategist (Agent 5), the deterministic TSC-Orchestrator that wires all five agents together over an event-sourced state, and the three clinician surfaces that present the result. Each requirement carries an ID, a description, and acceptance criteria stated against the synthetic cohort's ground truth — not against any claim of clinical validity, which is out of scope until the institutional Phase-1 work with real Epic Clarity and Discover Together Biobank data.

A standing reminder that applies to every requirement below: the synthetic-data demo is what runs now on the DGX Spark with RunPod overflow. The Epic Clarity/Caboodle and biobank-LIMS feeds that would supply these agents in production are described architecturally only and are explicitly not built. Where an FR's input would in production come from those feeds, it is satisfied in the demo by the version-controlled 50-patient synthetic cohort.

### 4.1 TSC-Trajectory Modeler (Agent 3) — FR-TM-*

The Trajectory Modeler is the one agent in the engine that is deliberately **not LLM-driven** for its core function. Forecasting SEGA growth, AML growth, seizure burden, and renal function is a classical longitudinal-statistics problem, and an LLM is the wrong tool for it: it cannot produce calibrated prediction intervals, it cannot be audited as a fitted model, and it would be the easiest thing in the demo for a skeptical reader to dismiss. Claude appears here only at the edges — Haiku for prose rendering of a computed forecast, Sonnet for interpreting an unusual trajectory that the statistical layer flags. The numbers come from fitted models with reproducible parameters.

**FR-TM-1 — Longitudinal measurement assembly.** The agent SHALL assemble, per patient and per tracked quantity, a time-ordered series of measurements from the synthetic imaging reports (SEGA max diameter, AML max diameter), the structured data (eGFR, seizure-frequency counts from the Phenome Mapper's HPO timeline), each point carrying value, units, measurement date, source document ID, and a measurement-uncertainty estimate.
- *Acceptance:* For Patient B the SEGA series resolves to the canonical `[0.8, 1.1, 1.3] cm` at the foramen of Monro across the cohort's defined imaging dates, with each point traceable to its source report.

**FR-TM-2 — Model fitting per quantity.** For each quantity with ≥3 longitudinal points, the agent SHALL fit the configured model class: a **linear or piecewise-linear mixed-effects model** for SEGA/AML diameter, **Gaussian process regression** (Matérn-5/2 kernel default) where curvature/non-monotonic behavior is plausible, **survival analysis** (time-to-threshold-crossing) for clinically defined thresholds, and **Bayesian updating** of the cohort-level prior as each patient's points arrive. Model class, hyperparameters, and random seed SHALL be read from `trajectory_config.yaml` so a fit is exactly reproducible.
- *Acceptance:* Re-running a fit with the same config and cohort version reproduces point forecasts and interval bounds to within floating-point tolerance; the chosen model class is recorded in the output provenance.

**FR-TM-3 — Forecast horizons and prediction intervals.** The agent SHALL emit forecasts at **6, 12, and 18 months** for each quantity, each with a point estimate and **50% and 90% prediction intervals**. Intervals SHALL widen with horizon and reflect the per-patient data density (sparser series → wider intervals).
- *Acceptance:* Coverage of the 90% interval over the synthetic cohort's held-out future points is within the demo tolerance band (target 0.85–0.95); intervals for a 2-point-sparse patient are visibly wider than for a dense series.

**FR-TM-4 — Threshold-crossing alerts.** The agent SHALL evaluate forecasts against the configured clinical thresholds (SEGA growth implying hydrocephalus risk at the foramen of Monro; AML ≥ ~4 cm bleeding-risk threshold; eGFR decline bands) and emit a threshold-crossing event when the **point estimate or the 50% interval** crosses within a horizon, classified as `likely` (point crosses) or `possible` (interval crosses).
- *Acceptance:* For Patient B the agent forecasts the SEGA crossing its watch threshold within the **12–18 month window** and emits a `possible`-then-`likely` graded alert; it raises **no false threshold alarm** on the cohort's stable-AML and stable-renal patients (zero false alarms is the eval target).

**FR-TM-5 — Surveillance-cadence recommendation.** From the fitted growth rate and interval width, the agent SHALL output a recommended next-imaging interval consistent with the ITSC 2021 surveillance cadence as a floor, tightening (never loosening below) it when growth velocity or interval width warrants. The recommendation is decision-support text, not an order.
- *Acceptance:* A patient with an accelerating SEGA slope receives a tighter-than-baseline recommended interval with the ITSC baseline cited; a stable patient receives the ITSC baseline.

**FR-TM-6 — LLM only at the edges.** Haiku SHALL render the computed forecast into a one-paragraph plain-language summary that introduces no number absent from the fitted output; Sonnet SHALL be invoked only when FR-TM-2's diagnostics flag an unusual trajectory (poor fit, regime change, outlier). No LLM SHALL produce a forecast value.
- *Acceptance:* A diff of every numeric token in the Haiku summary against the model output is exact; the Sonnet path fires only on flagged cases and is logged as such.

**I/O schema (output):**

```json
{
  "patient_id": "B",
  "generated_at": "2026-07-15T14:22:03Z",
  "forecasts": [{
    "quantity": "sega_max_diameter_cm",
    "model_class": "piecewise_linear_mixed_effects",
    "series": [{"date": "2025-01-10", "value": 0.8, "source_doc": "img_0142"}],
    "horizons": [{"months": 18, "point": 1.7, "pi50": [1.5, 1.9], "pi90": [1.2, 2.3]}],
    "threshold_events": [{"threshold": "foramen_monro_watch", "grade": "likely", "horizon_months": 18}],
    "surveillance_recommendation": {"next_interval_months": 6, "baseline_itsc_months": 12}
  }],
  "provenance": {"config_hash": "…", "seed": 42, "fit_diagnostics": {"r2": 0.97}}
}
```

**Endpoints:** `POST /agents/trajectory/fit` (cohort or single patient), `GET /agents/trajectory/{patient_id}/forecast`. Implementation is Python (statsmodels / scikit-learn GP / lifelines) behind the standard FastAPI agent service scaffolded from the precision_oncology_agent template. Build tasks land in **W4** (Patient B SEGA the gating case).

### 4.2 TAND Surveillance Agent (Agent 4) — FR-TS-*

The TAND Surveillance Agent surfaces under-recognized TSC-Associated Neuropsychiatric Disorders, the features the TOSCA registry showed are missed or unaddressed in 30–50% of patients despite affecting ~90%. Its method is a direct extension of Hagedorn's published clinical-NLP work on diagnostic-uncertainty language (Marshall/Nickels/Brady/Hagedorn 2023; Nickels et al. 2024) — it is an extension of his team's own research, not an external graft. The agent reads longitudinal notes for the linguistic *markers of uncertainty around* the six TAND clusters, not for diagnoses. It never diagnoses, never interrupts; it produces pre-visit briefing material.

**FR-TS-1 — Per-note discourse analysis (Sonnet).** For each clinical note in a patient's longitudinal record, the agent SHALL run a Sonnet pass that, for each of the **six TAND clusters** (behavioral, psychiatric, intellectual, academic, neuropsychological, psychosocial), extracts spans exhibiting the **Marshall-Hagedorn diagnostic-uncertainty discourse markers**: hedging, deferral, third-party attribution, conditional framing, and follow-up-without-formalization. Each extracted span carries the cluster label, the marker type(s), the verbatim text span with character offsets, and the note ID/date.
- *Acceptance:* On the synthetic notes with deliberately embedded TAND signals (Patient B's scattered under-recognized signals being the canonical set), every embedded signal is recovered with its correct cluster and at least one correct marker type.

**FR-TS-2 — Deterministic scoring and aggregation.** A **non-LLM** layer SHALL aggregate per-note spans into a per-cluster longitudinal signal score using configured weights per marker type and a recency/recurrence function (a marker recurring across visits weighs more than a one-off). Scoring is deterministic and reproducible from `tand_scoring.yaml`.
- *Acceptance:* Re-running aggregation over fixed per-note extractions yields identical scores; a marker appearing in three successive notes scores higher than the same marker once.

**FR-TS-3 — Spurious-flag suppression.** The scoring layer SHALL apply a minimum-evidence threshold and a negation/normal-finding filter so that explicitly negated or routinely-normal statements do not generate signal.
- *Acceptance:* On cohort patients with no embedded TAND signal the agent produces **no spurious flag** (the eval target); negated statements ("denies behavioral concerns") generate no signal.

**FR-TS-4 — Briefing summary (Opus).** For clusters crossing the briefing threshold, Opus SHALL compose a concise briefing paragraph per cluster that names the pattern, cites the contributing note spans by date, and frames it as something to consider asking about — never as a diagnosis or a directive. The TAND-L lifetime checklist is referenced as the structured backbone.
- *Acceptance:* Every claim in the briefing maps to a cited span; the language contains no diagnostic assertion; output reads as a question to raise, not a finding to act on.

**FR-TS-5 — Surface posture.** TAND output SHALL be delivered only as pre-visit briefing material and the in-visit dashboard's TAND quadrant. It SHALL NOT generate interruptive alerts.
- *Acceptance:* No TAND output path writes to the async alert surface.

**I/O schema (output):**

```json
{
  "patient_id": "B",
  "clusters": [{
    "cluster": "behavioral",
    "signal_score": 0.71,
    "spans": [{"note_id": "n_0231", "date": "2025-03-04", "marker_types": ["hedging", "deferral"],
               "text": "mother mentions he 'sometimes' has trouble settling; will revisit next visit",
               "offsets": [412, 498]}],
    "briefing": "Across two visits, settling/behavioral concerns are raised tentatively by the parent and deferred without formal assessment — worth asking about directly."
  }],
  "provenance": {"sonnet_prompt_version": "tand-v3", "scoring_config_hash": "…"}
}
```

**Endpoints:** `POST /agents/tand/analyze` (patient), `GET /agents/tand/{patient_id}/briefing`. Build in **W5**, validated against Patient B's embedded signals.

### 4.3 TSC-Therapeutics Strategist (Agent 5) — FR-TR-*

The Therapeutics Strategist is the integrative agent and is **Opus-class, non-negotiable** — it reasons across all four prior agents plus medications/adherence, adverse events, a PubMed/PMC RAG retrieval over the TSC corpus partition, a ClinicalTrials.gov snapshot, and FDA actions. Its output is a structured options brief framed as decision-support, never a recommendation, with every claim source-attributed.

**FR-TR-1 — Integrative context assembly.** The agent SHALL assemble a single context object from: Variant Curator interpretation, Phenome Mapper HPO profile, Trajectory forecasts/threshold events, TAND signals, and the patient's synthetic medication list, adherence notes, and adverse-event log.
- *Acceptance:* For Patient C the context includes the TSC1 interpretation, partial-everolimus-response history, the mucositis-driven dose reduction, the ~4 cm AML, and the refractory focal seizures.

**FR-TR-2 — RAG retrieval with attribution.** The agent SHALL retrieve from the Milvus TSC corpus partition (BAAI/bge-large-en-v1.5 + clinical embeddings) and from a static ClinicalTrials.gov + FDA-actions snapshot, returning passages with source URIs. Every factual claim in the brief SHALL carry an attribution to a retrieved source or to an upstream agent output.
- *Acceptance:* Each of the brief's claims links to a source URI or agent-output ID; an unattributable claim blocks brief completion (fails closed).

**FR-TR-3 — Six-section structured brief.** The agent SHALL emit exactly six sections: **Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions.** Trial Matching SHALL list ClinicalTrials.gov entries with match rationale and eligibility caveats; Emerging Evidence SHALL flag next-gen selective mTORC1 inhibitors and EXIST-3-class evidence with appropriate recency hedging.
- *Acceptance:* For Patient C the brief surfaces the **correct trial matches** against the snapshot, uses **appropriate hedging** on emerging evidence, and carries **full attribution**; generation completes in **< 3 min** (the eval target).

**FR-TR-4 — Decision-support framing and human gate.** The brief SHALL be framed throughout as options for clinician consideration. It SHALL NOT contain an imperative recommendation, dosing order, or autonomous action. A watermark and the synthetic-data banner SHALL be present.
- *Acceptance:* Linguistic check finds no imperative recommendation; the output is labeled decision-support.

**Output schema (abridged):**

```json
{
  "patient_id": "C",
  "sections": {
    "current_therapy": {"text": "...", "claims": [{"text": "...", "source": "agent:variant_curator"}]},
    "trial_matching": {"trials": [{"nct_id": "NCT…", "rationale": "...", "eligibility_caveats": "...",
                                   "source": "ctgov_snapshot_2026_06"}]},
    "open_questions": {"items": ["..."]}
  },
  "provenance": {"model": "claude-opus-4-8", "rag_sources": ["pmc://…"], "latency_ms": 168000}
}
```

**Endpoints:** `POST /agents/therapeutics/brief`, `GET /agents/therapeutics/{patient_id}/brief`. Build in **W5–W6**; Patient C is the gating case. Per the de-scope ladder, the **last** thing cut if W6 slips is reducing the brief from six sections to four (Combination and Emerging Evidence drop first).

### 4.4 TSC-Orchestrator — FR-OR-*

The orchestrator is a **deterministic LangGraph event router, not an LLM.** It owns the engine's coordination logic and its event-sourced state. Its job is to decide which agent runs when, to minimize recomputation, and to fail conservatively so a surface never silently drops an output.

**FR-OR-1 — Thirteen event types.** The orchestrator SHALL recognize and route these 13 event types:

| # | Event | Trigger | Primary effect |
|---|-------|---------|----------------|
| 1 | `patient_enrolled` | new patient added | dependency-ordered agent bootstrap |
| 2 | `genomic_data_ingested` | BAM/VCF available | enqueue Variant Curator |
| 3 | `variant_interpretation_ready` | Agent 1 done | unblock Therapeutics context |
| 4 | `clinical_note_ingested` | new note | enqueue Phenome Mapper + TAND (incremental) |
| 5 | `phenome_profile_updated` | Agent 2 done | unblock Trajectory, TAND, Therapeutics |
| 6 | `imaging_report_ingested` | new imaging report | enqueue Trajectory measurement assembly |
| 7 | `trajectory_forecast_ready` | Agent 3 done | evaluate threshold events |
| 8 | `threshold_crossing_detected` | FR-TM-4 | emit async alert candidate |
| 9 | `tand_signal_updated` | Agent 4 done | update briefing material |
| 10 | `therapeutics_brief_requested` | clinician/demand-driven | enqueue Agent 5 |
| 11 | `therapeutics_brief_ready` | Agent 5 done | mark in-visit/briefing surface ready |
| 12 | `surface_requested` | clinician opens a surface | demand-driven assembly |
| 13 | `agent_failed` | any agent error/timeout | conservative degradation (mark `pending`/`stale`) |

**FR-OR-2 — Dependency-ordered enrollment.** On `patient_enrolled`, the orchestrator SHALL run the **Phenome Mapper first** (it is the foundation), then fan out to the agents that depend on its profile.
- *Acceptance:* The event log for an enrollment shows Phenome Mapper started before Trajectory/TAND/Therapeutics for that patient.

**FR-OR-3 — Incremental-update minimization.** On a single new note or report, the orchestrator SHALL recompute only the affected agents and only the affected patient, never the whole cohort.
- *Acceptance:* `clinical_note_ingested` for one patient triggers Phenome Mapper + TAND for that patient only; the Variant Curator is not re-run.

**FR-OR-4 — Demand-driven surface assembly.** On `surface_requested` the orchestrator SHALL assemble the surface from the current-state projection, computing nothing that is already materialized and fresh.
- *Acceptance:* Opening an unchanged patient's dashboard issues zero agent recomputations.

**FR-OR-5 — Conservative failure handling.** On `agent_failed`, the orchestrator SHALL mark the affected surface region `pending` or `stale` with the failure reason and timestamp, and SHALL NOT silently omit the output or substitute a fabricated one.
- *Acceptance:* A forced Trajectory failure renders the dashboard's forecast quadrant as `pending` with reason text; the rest of the dashboard renders normally.

**FR-OR-6 — Event-sourced state.** State SHALL be an **append-only PostgreSQL event log** plus **materialized current-state projections** rebuildable by replay; Redis holds ephemeral working state with TTL; the demo's wiring is read from a YAML demo config. Every projection row records the event offset it was built from.
- *Acceptance:* Dropping and replaying the event log reconstructs identical current-state projections; no projection is writable except by replay.

**State model (event log row):**

```json
{
  "event_id": "evt_000142", "offset": 142, "type": "trajectory_forecast_ready",
  "patient_id": "B", "ts": "2026-07-15T14:22:05Z",
  "payload": {"forecast_id": "fc_77"}, "actor": "agent:trajectory", "input_hash": "…"
}
```

**Endpoints:** `POST /orchestrator/events` (append), `GET /orchestrator/state/{patient_id}` (projection), `POST /orchestrator/replay` (rebuild). LangGraph defines the routing graph; nodes are agent-invocation or projection-update functions, edges are event-type guards. Orchestrator integration is the **W7** focus.

### 4.5 Clinician Surfaces — FR-SF-*

The three surfaces are **standalone web apps** (Streamlit/FastAPI from the template), **not Epic embeds**, and are **persistently watermarked synthetic.** They share the orchestrator's current-state projection as their single source of truth and never call agents directly.

**FR-SF-1 — Pre-visit briefing surface.** A **single-screen, mobile-readable** layout: header (patient, visit context), a "what's new since last visit" block, **0–3 action items** maximum, a watchlist, and source links. It composes the TAND briefing, trajectory threshold events, and a therapeutics-brief link.
- *Acceptance:* For Patient B the briefing renders on a phone-width viewport with ≤3 action items, a "what's new" reflecting the latest SEGA reading, and the under-recognized TAND signals as watchlist items with source links.

**FR-SF-2 — In-visit dashboard surface.** A **four-quadrant** layout — (1) variant interpretation, (2) HPO timeline, (3) trajectory forecasts with 50%/90% intervals, (4) TAND + therapeutics — with **progressive disclosure** (summary first, drill to spans/sources on demand) and **source navigation** from any claim to its provenance.
- *Acceptance:* Patient B's dashboard shows all four quadrants populated; clicking a forecast opens its fit diagnostics and source measurements; clicking a TAND item opens the cited note span; a `pending` quadrant (FR-OR-5) renders as such.

**FR-SF-3 — Async alert surface.** Four alert categories with **strict discipline**: it SHALL track alert volume and recalibrate thresholds if a clinician would receive **more than ~3 alerts/week.** TAND never routes here (FR-TS-5). Threshold-crossing (FR-TM-4) is the primary source.
- *Acceptance:* Over a simulated cohort-week the surface stays at or below the ~3/clinician/week budget; exceeding it triggers a logged recalibration prompt.

**FR-SF-4 — Watermark and provenance affordance.** Every surface SHALL render a persistent synthetic-data watermark and expose, on every displayed output, the model id/version, prompt-template version, retrieved sources, input hash, and latency from the producing agent's provenance.
- *Acceptance:* The Act One audit-trail demo opens the full provenance for Patient A's mosaic interpretation directly from the surface.

**Build tasks.** Surfaces scaffold in **W6** (against Patient C therapeutics and the briefing) and complete in **W7** alongside orchestrator integration; **W8** is regeneration, dry runs, clinician review, and delivery. The de-scope ladder simplifies the alert surface before the TAND cluster set, and never touches Act One.

These functional requirements define behavior and interfaces; Section 5 turns to the non-functional requirements, the synthetic-cohort data specification, and the evaluation protocol that scores each FR's acceptance criteria against ground truth.


## 5. Non-Functional Requirements, Data Specs, Evaluation Harness & Integration

Section 4 specified what each of the five agents and the deterministic orchestrator produce. This section specifies the constraints those components run under, the data they consume, how we measure whether they work, and — critically for a skeptical CCHMC audience — which integration surfaces are real in v0.1 versus described architecturally as institutional Phase-1 work. The governing principle throughout: the build is the argument. Every requirement here is testable against the 50-patient synthetic cohort running on the DGX Spark, with RunPod GPUs attached for the genomics and cohort-generation bursts. Nothing in this section depends on access to real patient data, an Epic instance, or a biobank LIMS.

### 5.1 Non-Functional Requirements

NFRs carry stable IDs (`NFR-<domain>-<n>`) so the build plan in Section 6 and the evaluation harness in 5.3 can reference them directly. Each is phrased as a verifiable condition against the synthetic cohort, not an aspiration.

#### 5.1.1 Performance / latency

Latencies are wall-clock on the Spark (GB10 Grace Blackwell, 128 GB unified LPDDR5x) for LLM-light paths, with genomics offloaded to RunPod-attached GPUs for Parabricks-equivalent calling. They are budgets for the demo, not SLAs for production.

| ID | Requirement | Target | Notes |
|----|-------------|--------|-------|
| NFR-PERF-1 | Variant Curator end-to-end per case (BAM-in to draft report) | < 5 min/case | Excludes one-time cohort BAM staging; mosaic-aware calling dominates |
| NFR-PERF-2 | Phenome Mapper full 50-patient cohort | < 1 hr | Per-note Sonnet + Haiku normalization; parallelizable |
| NFR-PERF-3 | Phenome Mapper single incremental note (orchestrator-triggered) | < 30 s | Incremental-update path, not full re-run |
| NFR-PERF-4 | Trajectory Modeler per-patient forecast (all four endpoints) | < 20 s | Classical stats; CPU-bound, no GPU |
| NFR-PERF-5 | TAND Surveillance per-patient briefing pass | < 90 s | Per-note discourse analysis over 6 clusters |
| NFR-PERF-6 | Therapeutics Strategist six-section brief | < 3 min | Opus-class + RAG retrieval; the heaviest single call path |
| NFR-PERF-7 | Pre-visit briefing surface assembly (cached agent outputs) | < 2 s | Demand-driven read of materialized state; no agent invocation |
| NFR-PERF-8 | In-visit dashboard initial render | < 3 s | Progressive disclosure; quadrants hydrate independently |

NFR-PERF-1 and NFR-PERF-6 are the demo-critical ones: Act One depends on sub-5-minute mosaic recovery on Patient A, and the Patient C therapeutics walk in Act Two depends on a sub-3-minute brief.

#### 5.1.2 Reliability, reproducibility, determinism

The credibility of a draft molecular-genetics report rests on reproducibility. We distinguish three determinism tiers and hold each component to the tightest tier it can honor.

- **Tier D (fully deterministic):** variant calling, snpEff/VEP annotation, ACMG combinatorial rule evaluation, the Trajectory Modeler's statistical core, the orchestrator's event routing, and all scoring/aggregation layers. Same input bytes -> identical output bytes, verified by hashing.
- **Tier R (reproducible-under-pin):** all LLM calls. Temperature pinned (0.0 for classification/extraction synthesis; <=0.3 only for prose summaries), model id/version pinned, prompt template version pinned, retrieved RAG source set pinned per call. We do not claim bit-identical LLM output; we claim that the inputs to every LLM call are fully recorded and that re-running with the same pins lands in the same classification/extraction decision.
- **Tier C (controlled-stochastic):** none in the shipped path. No component is allowed free-temperature generation in an output that reaches a clinician surface.

| ID | Requirement |
|----|-------------|
| NFR-REL-1 | The 50-patient synthetic cohort regenerates deterministically from a pinned seed in <= 12 hr; regenerated artifacts are byte-identical for Tier-D layers and content-equivalent for note/imaging text |
| NFR-REL-2 | Every agent output carries a provenance record (5.1.5); no output without provenance reaches a surface |
| NFR-REL-3 | ACMG synthesis (Opus) is validated against the deterministic combinatorial-rule evaluation on every case; disagreement is surfaced, never silently overridden |
| NFR-REL-4 | Orchestrator failure handling is conservative: a failed agent yields a surface state of "pending"/"stale," never a silently missing or fabricated output |
| NFR-REL-5 | All state transitions are event-sourced (append-only PostgreSQL event log + materialized projections); current state is reconstructable by replaying events |
| NFR-REL-6 | RAG retrieval is reproducible: corpus partition is version-tagged; retrieved chunk URIs and scores are recorded per call (Tier R) |

#### 5.1.3 Cost

Cost discipline is part of the thesis — the platform runs on a $4,699 Spark plus burst GPU rental, not a cluster.

| ID | Requirement | Cap |
|----|-------------|-----|
| NFR-COST-1 | Per-patient full engine pass (all 5 agents, cold) Claude API spend | < $2.50 |
| NFR-COST-2 | Full 50-patient cohort cold pass Claude API spend | < $100 |
| NFR-COST-3 | One-time cohort generation RunPod GPU spend (genomics + parallel note/imaging generation) | < $150 |
| NFR-COST-4 | Steady-state demo operation (Spark only, agents idle/on-demand) marginal infra cost | ~$0 |

A per-tier token/cost ledger is recorded for every call (Haiku/Sonnet/Opus + local Llama 3.1 70B fallback via Ollama), mirroring the cost-ledger pattern already in the HCLS AI Factory substrate. The Llama fallback exists so the demo degrades gracefully — and remains *runnable offline* — if the API is unreachable; fallback invocations are flagged in provenance and excluded from clinical-quality claims.

#### 5.1.4 Licensing

| ID | Requirement |
|----|-------------|
| NFR-LIC-1 | All net-new TSC Engine code (5 agents, orchestrator, cohort pipeline, surfaces) ships under Apache 2.0 |
| NFR-LIC-2 | All external tools in the build path are license-compatible and attributed: Synthea (Apache 2.0), BAMSurgeon (MIT), snpEff (LGPL), VEP (Apache 2.0), Milvus (Apache 2.0), LangGraph (MIT) |
| NFR-LIC-3 | Reference databases are used within their terms: gnomAD (open), ClinVar/dbSNP (public domain), HPO (open), SNOMED-CT (license required for production; demo uses only HPO/ICD-10 mappings to avoid SNOMED redistribution), LOVD-TSC (open, cite) |
| NFR-LIC-4 | Claude (API) and any proprietary model weights are *runtime dependencies*, not redistributed; the repo contains no model weights |

SNOMED-CT carries a licensing constraint we honor conservatively: the synthetic demo normalizes to HPO and ICD-10 only, and any SNOMED dependency is deferred to institutional Phase-1 where CCHMC's existing license applies.

#### 5.1.5 Privacy, no-real-data, watermarking

| ID | Requirement |
|----|-------------|
| NFR-PRIV-1 | No real patient data, real DICOM, or real FASTQ enters the v0.1 system. The genomic substrate starts from BAMSurgeon edits over the public NA12878-derived BAM; clinical text and imaging reports are frontier-model-generated from published templates |
| NFR-PRIV-2 | Every generated clinical note, imaging report, and surface view is persistently watermarked "SYNTHETIC — NOT A REAL PATIENT," including in exported/printed views |
| NFR-PRIV-3 | No IRB is required for the synthetic demo; this is asserted in-product and in the demo deck. Real-data Phase-1 requires IRB and is explicitly out of v0.1 scope |
| NFR-PRIV-4 | The system is not FDA-cleared; SaMD posture is undetermined and is institutional work. Every clinician surface states this and frames output as decision-support with a human gate |
| NFR-PRIV-5 | Provenance records and audit logs contain no PHI by construction (synthetic only); the audit schema is nonetheless built PHI-safe so it ports to Phase-1 unchanged |

The watermark is not cosmetic. With a skeptical clinical audience, the fastest way to lose the room is ambiguity about whether a screen shows a real child. NFR-PRIV-2 makes that impossible to mistake.

### 5.2 Data Specifications

Three data classes feed the engine: reference databases (read-only, version-pinned), the RAG corpus (the Therapeutics Strategist's evidence substrate), and the synthetic cohort artifacts (generated, version-controlled).

#### 5.2.1 Reference databases

All references are pinned to a named release and stored under `data/reference/<db>/<version>/`. A `MANIFEST.json` records each version, download date, and checksum so NFR-REL-1 reproducibility holds.

| Database | Pinned version | Used by | Role |
|----------|----------------|---------|------|
| gnomAD | v4.0 | Variant Curator | Population allele frequency for PM2/BA1/BS1 |
| ClinVar | 2026-04 release | Variant Curator | Prior interpretations, ACMG evidence |
| dbSNP | b156 | Variant Curator | rsID annotation |
| LOVD-TSC | 2026-Q1 snapshot | Variant Curator | TSC1/TSC2 locus-specific variant evidence |
| HPO | 2026-04 release | Phenome Mapper, TAND | Phenotype ontology + ICD-10/lab->HPO maps |
| ClinicalTrials.gov | dated snapshot | Therapeutics Strategist | Trial matching (static snapshot, not live) |
| FDA actions | dated snapshot | Therapeutics Strategist | Label/safety actions for mTOR inhibitors |

ClinicalTrials.gov and FDA are *snapshots*, not live feeds — pinned for reproducibility (NFR-REL-6) and to make the demo offline-capable. Liveness is a Phase-1 concern.

#### 5.2.2 RAG corpus

The Therapeutics Strategist and the rare-conflict paths of other agents retrieve from a dedicated TSC corpus partition in Milvus, embedded with BAAI/bge-large-en-v1.5 plus BiomedBERT-derived clinical embeddings (reusing the HCLS AI Factory v1.3.0 substrate). Corpus contents:

- ITSC 2021 consensus surveillance and diagnostic guidelines (the spine of surveillance-gap reasoning).
- Somatic-mosaicism literature: Tyburczy 2015, Giannikou 2016, Lim 2017 — the evidentiary basis for the NMI/low-VAF recovery claim.
- EXIST-3 and the targeted-epilepsy-therapy literature; mTOR-inhibitor pharmacology (everolimus/sirolimus) including adverse-event profiles (mucositis, the Patient C dose-reduction scenario).
- TOSCA registry findings and the TAND consensus framework (TAND-L checklist).
- A curated PubMed/PMC TSC subset, ingested as full-text where open-access, abstract-only otherwise.

Ingestion is a deterministic pipeline: source -> normalize -> chunk (semantic, ~512-token target with overlap) -> embed -> upsert into the version-tagged `tsc_corpus_v<n>` partition. Each chunk retains a source URI, citation key (author/year), and license flag. NFR-REL-6 requires that retrieval records the partition version and returned chunk URIs/scores in provenance.

#### 5.2.3 Synthetic-data formats

The 50-patient cohort (30 TSC2 germline / 12 TSC1 germline / 5 TSC2 mosaic / 2 TSC1 mosaic / 1 NMI) is the single source of truth for the demo and for evaluation. Artifacts per patient, under `data/cohort/v<n>/patient_<id>/`:

```
patient_<id>/
  clarity/            # Epic Clarity-shaped relational extract (CSV/Parquet), Synthea-derived
  fhir/               # FHIR R4 bundle (Synthea native output)
  genomics/
    aligned.bam       # BAMSurgeon-edited, NA12878-derived; mosaic VAF 4-12%, germline ~50%
    variants.vcf.gz   # Parabricks-equivalent calling output
    truth.json        # ground-truth variant + intended VAF + ACMG class (eval only)
  notes/              # frontier-model clinical notes, watermarked, with HPO ground-truth spans
  imaging/            # frontier-model imaging *reports* (text); NO DICOM
  cohort_meta.json    # genotype class, featured-patient flag, seed lineage
```

The data the cohort deliberately does **not** contain is as load-bearing as what it does: no real imaging/DICOM, no raw FASTQ (the substrate starts at BAM), no neuropsych test scores, no pedigree beyond what appears in notes, no pharmacy/claims data, no patient-reported-outcome scores. Stating these gaps explicitly pre-empts the "you're hiding the hard part" objection — the hard parts are named, not papered over.

Ground-truth files (`truth.json`, HPO span annotations) are stored alongside but are read **only** by the evaluation harness, never by the agents. The agents see the same artifacts a Phase-1 clinician pipeline would.

### 5.3 Evaluation Harness

The harness is a standalone, scriptable runner (`eval/run_eval.py`) that executes each agent over the cohort, compares output to the embedded ground truth, and emits a pass/fail report against the NFRs and the per-agent eval targets from Section 4. It is explicitly **not clinical validation** — it measures behavior against *synthetic* ground truth we ourselves constructed. That limitation is stated in every report header.

#### 5.3.1 Ground-truth construction

Ground truth is generated *before* and *independently of* agent output, then frozen:

- **Variant truth:** the BAMSurgeon edit spec is the ground truth — we know exactly which TSC1/TSC2 variant was inserted at what intended VAF. ACMG class is pre-assigned by combinatorial-rule evaluation and a one-time clinician-sampled review.
- **Phenotype truth:** HPO terms are seeded into the Synthea modules and into the note-generation templates, so the intended HPO profile and the evidence spans are known by construction.
- **Trajectory truth:** SEGA/AML growth curves and seizure-burden trajectories are parameterized in the generator (e.g., Patient B's SEGA at ~2-4 mm/yr crossing threshold in the 12-18 month window), so threshold-crossing timing is known.
- **TAND truth:** TAND signals are deliberately embedded into specific notes at known severities; the harness checks recall of embedded signals and counts spurious flags.

#### 5.3.2 Metrics and gating

| Agent | Metric | Gate (demo) |
|-------|--------|-------------|
| Variant Curator | Mosaic-variant recovery at VAF >= 5% | 7/7 recovered |
| Variant Curator | ACMG class agreement vs truth | 100%, **zero** false-positive Pathogenic |
| Variant Curator | Per-case latency | < 5 min (NFR-PERF-1) |
| Phenome Mapper | HPO recall | >= 90% |
| Phenome Mapper | HPO precision | >= 85% |
| Phenome Mapper | Cohort latency | < 1 hr (NFR-PERF-2) |
| Trajectory Modeler | Patient B SEGA threshold-crossing forecast | within 12-18 mo window |
| Trajectory Modeler | False-alarm rate | zero on stable patients |
| TAND Surveillance | Embedded-signal detection | all embedded signals surfaced |
| TAND Surveillance | Spurious flags | zero |
| Therapeutics Strategist | Trial-match correctness | all valid matches, no invalid |
| Therapeutics Strategist | Attribution + hedging | every claim source-attributed; appropriate uncertainty language |
| Therapeutics Strategist | Latency | < 3 min (NFR-PERF-6) |

Gating policy: the **zero-false-positive-Pathogenic** gate on the Variant Curator is hard — a single FP fails the build, because in this domain an over-called Pathogenic is the most dangerous error. The TAND zero-spurious-flag gate is similarly hard, consistent with the agent's discipline of surfacing patterns as briefing material rather than firing alerts. Recall/precision and latency gates are graded thresholds. The harness writes a `eval_report_<run_id>.json` with per-patient detail and a one-page Markdown summary suitable for the demo appendix.

Calibration of the async alert surface is itself evaluated: if the simulated alert volume exceeds ~3 alerts per clinician per week over the cohort, thresholds are recalibrated — alert burden is a first-class failure mode (cf. Orenstein et al. 2021), not an afterthought.

### 5.4 Integration Requirements

This is the section that most directly governs credibility, so it states plainly what is real now versus what is institutional Phase-1. The rule: v0.1 runs entirely on synthetic data on the Spark + RunPod. Real institutional integrations are described architecturally and are **explicitly not built**.

| Integration | v0.1 status | Phase-1 (institutional) |
|-------------|-------------|-------------------------|
| Epic Clarity / Caboodle | **Stubbed** — Synthea-derived Clarity-shaped relational extracts stand in for Clarity table structures | Real read-only extract against CCHMC's Clarity/Caboodle, governed by data-use agreement |
| FHIR R4 | **Simulated** — Synthea-native FHIR bundles; the FHIR ingest path is exercised against synthetic bundles | Real FHIR endpoint integration; SMART-on-FHIR auth |
| Biobank LIMS | **Not built** — banked-tissue provenance is represented as fields in `cohort_meta.json` | Discover Together Biobank LIMS linkage: banked tuber/AML/SEGA specimen -> Variant Curator |
| Imaging / PACS | **Reports only** — frontier-model imaging *text reports*; no DICOM, no imaging-AI inference | PACS + imaging-AI pipeline (e.g., longitudinal SEGA volumetrics); a separate institutional effort |
| IAM / SSO | **Local** — demo auth only; standalone web apps, not embedded in Epic | Enterprise IAM/SSO, role-based access |
| Audit / SIEM | **Local append-only log** — PHI-safe by construction (NFR-PRIV-5) | Export to institutional SIEM; retention policy |

Two architectural commitments make the Phase-1 path credible without overclaiming it now. First, the **source-area model**: the five CCHMC areas (Discover Together Biobank, Biomedical Informatics/Hagedorn, the TSC clinical & research program, and the Clarity/Caboodle + LIMS data plumbing) are *sources that feed the engine*, not parts of it. The Winslow Research Pavilion is the physical envelope. In v0.1, each source is replaced by a synthetic stand-in with the same shape, so the wiring is real even when the box behind it is simulated. Replication to another site is "swap the box labels, keep the wiring" — a biobank without an intelligence layer is a freezer full of tubes.

Second, the **audit and provenance schema is built Phase-1-ready**: PHI-safe today (synthetic only), but structured so that turning on real Clarity/FHIR/LIMS feeds requires configuration and a data-use agreement, not a rewrite. The integration boundaries are interfaces, not assumptions baked through the agents.

With the constraints, data, evaluation gates, and integration boundaries fixed, Section 6 turns to the eight-week build plan that delivers against them — week by week, with the de-scope order that protects the demo-critical paths (Variant Curator mosaic recovery and the Therapeutics brief) when time slips.


## 6. Build Plan, Test Plan, Risk Register, Dependencies & Roadmap

This section is the operational contract for the eight-week build. It turns the architecture described earlier into a dated task breakdown, a test plan that culminates in the three-act demo runbook, a risk register I will actively manage, an honest dependency map, and a roadmap that connects the synthetic-data demo running now on the Spark to the institutional, real-data work that is explicitly out of scope for this phase. The recurring discipline throughout: the build is the argument. Every claim in the demo must be backed by something a skeptic on Hagedorn's team can open, inspect, and reproduce.

A note on ownership before the tables: the owner column reads "Adam" on every line, because this is a single-builder eight-week sprint on my DGX Spark with RunPod burst capacity. That is a deliberate constraint, not an oversight. A solo build with a fixed cohort and a fixed demo target is the cheapest honest way to test whether the engine is worth an institutional Phase 1. Where collaboration is assumed (clinician note review, the TSC faculty lead, Hagedorn's informatics team), it appears in the dependency map, not as build labor.

### 6.1 Build Plan: Task Breakdown by Week

The eight weeks map to the components in the order their dependencies allow. The Phenome Mapper foundation and the Variant Curator core land first because every downstream agent reads from them. The orchestrator integrates late but its event schema is fixed early (W1) so agents can be written against a stable contract. The roughly 61 tasks below are grouped by week; each is tagged with its primary dependency and a de-scope marker (•) where it sits on the cut list defined in 6.1.7.

#### W1 — Cohort Foundation & Substrate (T-01 .. T-08)

| ID | Task | Depends on |
|----|------|-----------|
| T-01 | Stand up TSC repo by scaffolding from `precision_oncology_agent` (Streamlit UI + FastAPI + Milvus client + Pydantic models + Claude tiered client) | HCLS v1.3.0 substrate |
| T-02 | Define the 13-event-type orchestrator event schema (Pydantic + JSON Schema) and freeze v1 | T-01 |
| T-03 | PostgreSQL append-only event log + materialized current-state projection tables; Redis TTL ephemeral store; MinIO buckets | T-01 |
| T-04 | Synthea + TSC disease modules: generate demographic/clinical skeletons for a 10-patient pilot slice (FHIR R4 + Clarity-shaped relational) | Synthea |
| T-05 | BAMSurgeon harness: insert TSC1/TSC2 variants into NA12878-derived BAMs; parameterize germline (~50% VAF) vs mosaic (4-12% VAF) | T-04 |
| T-06 | Parabricks variant-calling path on RunPod GPU (BWA-MEM + GATK HaplotypeCaller + Mutect2 mosaic-aware) producing VCFs from the pilot BAMs | T-05, RunPod |
| T-07 | Provenance envelope: shared record (model id/version, prompt template version, RAG source URIs, input hash, latency) wired into the FastAPI response wrapper | T-01 |
| T-08 | Synthetic watermark utility (persistent on every surface render and every generated artifact) | T-01 |

W1 exit gate: a 10-patient pilot slice exists end to end as BAMs, VCFs, and Clarity-shaped tables, and the event log records its first synthetic events.

#### W2 — Full Cohort & Variant Curator Core (T-09 .. T-17)

| ID | Task | Depends on |
|----|------|-----------|
| T-09 | Scale the cohort generator to the full 50: 30 TSC2 germline, 12 TSC1 germline, 5 TSC2 mosaic, 2 TSC1 mosaic, 1 NMI (7 mosaic total) | W1 gate |
| T-10 | Pin the cohort generation as deterministic + version-controlled (fixed seeds, ~12 hr full regen documented) | T-09 |
| T-11 | Variant Curator: deterministic snpEff/VEP annotation stage over cohort VCFs | T-06 |
| T-12 | Variant Curator: mosaic flag logic (VAF, read depth, strand bias, artifact heuristics) with ddPCR validation recommendation | T-11 |
| T-13 | Variant Curator: evidence-aggregation stage (Sonnet) consolidating ClinVar / gnomAD v4 / LOVD-TSC / dbSNP lookups | T-11, Milvus TSC partition |
| T-14 | Build the Milvus TSC corpus partition (BAAI/bge-large-en-v1.5 + BiomedBERT-derived clinical embeddings) over PubMed/PMC TSC literature | HCLS substrate |
| T-15 | Featured Patient A authored into the cohort: 4yo F, NMI, 8.3% VAF TSC2 frameshift in tuber tissue | T-09 |
| T-16 | Variant Curator FastAPI endpoint + per-case provenance output | T-07, T-12, T-13 |
| T-17 | Variant Curator unit tests against synthetic ground truth (see 6.2.1) | T-16 |

#### W3 — ACMG Synthesis & Phenome Mapper (T-18 .. T-26)

| ID | Task | Depends on |
|----|------|-----------|
| T-18 | ACMG-AMP classification synthesis (Opus) producing ClinVar-spec interpretation + AI-labeled draft molecular-genetics report | T-13 |
| T-19 | Validate Opus ACMG output against a deterministic combinatorial-rules reference implementation; log every disagreement | T-18 |
| T-20 | "Draft for molecular-geneticist sign-off" gating: no autonomous classification path; review state in the event log | T-18 |
| T-21 | Frontier-model clinical notes generator (~600-1000 notes) from published templates, watermarked | T-09 |
| T-22 | Phenome Mapper: per-note HPO extraction (Sonnet) with evidence spans | T-21, HPO/SNOMED-CT |
| T-23 | Phenome Mapper: ICD-10/lab → HPO normalization (Haiku) | T-22 |
| T-24 | Phenome Mapper: Opus rare-conflict resolution path | T-22 |
| T-25 | Phenome Mapper outputs: longitudinal HPO profile, discordance log, ITSC surveillance-gap report | T-22, T-23 |
| T-26 | Phenome Mapper recall/precision tests vs ground truth (≥90% / ≥85%) | T-25 |

#### W4 — Phenome Completion & Trajectory Modeler (T-27 .. T-34)

| ID | Task | Depends on |
|----|------|-----------|
| T-27 | Frontier-model imaging reports generator (brain MRI / renal US / echo / ophtho), longitudinal SEGA series ~2-4 mm/yr, watermarked | T-09 |
| T-28 | Featured Patient B authored: 12yo M, TSC2 c.3037C>T p.Arg1013Ter, SEGA 0.8→1.1→1.3 cm at foramen of Monro, bilateral 2.8 cm AMLs, well-controlled focal seizures, scattered under-recognized TAND signals | T-27 |
| T-29 | Trajectory Modeler: mixed-effects + Gaussian process regression for SEGA/AML growth (classical statistics, not LLM) | T-25, T-27 |
| T-30 | Trajectory Modeler: survival analysis + Bayesian updating for seizure burden and renal function; 50%/90% prediction intervals | T-29 |
| T-31 | Threshold-crossing alert logic + surveillance-cadence recommendations | T-30 |
| T-32 | Haiku prose-summary layer; Sonnet unusual-trajectory interpretation | T-30 |
| T-33 | Trajectory test: forecast Patient B SEGA crossing threshold in the 12-18 month window with no false alarms | T-31, T-28 |
| T-34 | Phenome Mapper cohort run completes < 1 hr end to end | T-25 |

#### W5 — TAND Surveillance & Therapeutics Start (T-35 .. T-42)

| ID | Task | Depends on |
|----|------|-----------|
| T-35 | Encode the Marshall-Hagedorn discourse-marker taxonomy (hedging, deferral, third-party attribution, conditional, follow-up-without-formalization) as a structured rubric | Marshall/Hagedorn 2023, Nickels 2024 |
| T-36 | TAND Agent: per-note Sonnet discourse analysis across the 6 TAND clusters (behavioral, psychiatric, intellectual, academic, neuropsychological, psychosocial) | T-35, T-21 |
| T-37 | TAND Agent: deterministic scoring/aggregation layer (non-LLM) | T-36 |
| T-38 | TAND Agent: Opus briefing-summary generation, framed as pre-visit briefing material, never interruptive, never diagnostic | T-37 |
| T-39 | TAND test: detect Patient B embedded signals, zero spurious flags | T-38, T-28 |
| T-40 | Therapeutics Strategist (Opus-class): integrate the four prior agents + meds/adherence + adverse events | T-25, T-30, T-38 |
| T-41 | Therapeutics RAG wiring: PubMed/PMC + ClinicalTrials.gov snapshot + FDA actions, all source-attributed | T-14, ClinicalTrials.gov, FDA |
| T-42 | Featured Patient C authored: 18yo F, TSC1, partial everolimus response with mucositis dose reduction, ~4 cm AML, refractory focal seizures | T-09 |

#### W6 — Therapeutics Completion & Clinician Surfaces (T-43 .. T-49)

| ID | Task | Depends on |
|----|------|-----------|
| T-43 | Therapeutics six-section options brief: Current Therapy, Optimization, Combination, Trial Matching, Emerging Evidence, Open Questions — every claim attributed | T-40, T-41 |
| T-44 | Therapeutics test: correct trial matches, appropriate hedging, full attribution, < 3 min/case | T-43, T-42 |
| T-45 | Pre-visit briefing surface (1-screen mobile-readable: header / what's-new / 0-3 action items / watchlist / links) | T-38, T-43 |
| T-46 | In-visit 4-quadrant dashboard (variant interp / HPO timeline / trajectory forecasts / TAND+therapeutics) with progressive disclosure + source navigation | T-18, T-25, T-30, T-43 |
| T-47 | Async alert surface (4 categories, strict discipline, recalibrate if > ~3 alerts/clinician/week) • | T-31, T-38 |
| T-48 | Surface watermarking + standalone web-app packaging (explicitly not Epic) | T-08, T-45, T-46 |
| T-49 | Source-navigation drill-through from any surface claim to its provenance record | T-07, T-46 |

#### W7 — Orchestrator Integration & End-to-End (T-50 .. T-56)

| ID | Task | Depends on |
|----|------|-----------|
| T-50 | TSC-Orchestrator: deterministic LangGraph event router over the 13 event types (not an LLM) | T-02, all agents |
| T-51 | Dependency-ordered enrollment (Phenome Mapper first) | T-50 |
| T-52 | Incremental-update minimization (only recompute affected agents on new events) | T-50 |
| T-53 | Demand-driven surface assembly (compose a surface only when requested) | T-50, T-45, T-46 |
| T-54 | Conservative failure handling: a failed agent yields "pending"/staleness on the surface, never a silent missing output | T-50 |
| T-55 | YAML demo-config wiring (which patients, which surfaces, which scripted events) | T-50 |
| T-56 | Full 50-patient end-to-end integration run through the orchestrator | T-50 .. T-55 |

#### W8 — Regen, Dry Runs, Clinician Review & Delivery (T-57 .. T-61)

| ID | Task | Depends on |
|----|------|-----------|
| T-57 | Clean cohort regeneration from seeds (proves determinism end to end, ~12 hr) | T-10, T-56 |
| T-58 | Clinician-sampled review of generated notes/imaging realism (TSC faculty lead) | dependency: clinician |
| T-59 | Three-act demo runbook authored + rehearsed; timing held to 30 min + 15 discussion | T-56 |
| T-60 | Three full dry runs against the regenerated cohort; defect triage between runs | T-57, T-59 |
| T-61 | Delivery package: repo (Apache 2.0), runbook, eval report vs ground truth, provenance samples, known-limitations doc | T-60 |

#### 6.1.7 De-scope order (if the schedule slips)

Cut in this order, top first. Each cut preserves the integrity of Act One (mosaic recovery), which is the load-bearing demonstration:

1. **Cohort 50 → 30** — drop the long tail of germline TSC2 patients; keep all 7 mosaic patients, Patients A/B/C, and the NMI case.
2. **Simplify the alert surface** (T-47) — collapse to one or two categories or defer it; the pre-visit briefing and in-visit dashboard carry the demo.
3. **Simplify the TAND cluster set** — run 3-4 of the 6 clusters rather than all 6.
4. **Last resort: Therapeutics six sections → four** — cut Combination and Emerging Evidence, keep Current Therapy, Optimization, Trial Matching, Open Questions.

Act One is never on the cut list. If the mosaic recovery does not work, there is no reason to do the demo.

### 6.2 Test Plan

Three layers, each feeding the next, ending in the demo runbook acceptance. Every test asserts against the synthetic ground truth baked into the cohort generator. This is demo validation, not clinical validation; that distinction is stated explicitly in the delivery package and is restated wherever the engine touches a real patient claim.

#### 6.2.1 Unit tests (per agent, vs synthetic ground truth)

- **Variant Curator** — recover all 7 mosaic variants at VAF ≥ 5%; assign the correct ACMG-AMP class for each cohort variant; zero false-positive Pathogenic calls; < 5 min/case. Patient A specifically must resolve to Likely Pathogenic (PVS1 + PM2 + PP4) with a ddPCR validation recommendation.
- **Phenome Mapper** — HPO extraction recall ≥ 90%, precision ≥ 85% against the per-note ground-truth HPO labels; evidence spans present on every extracted term; full cohort < 1 hr.
- **Trajectory Modeler** — Patient B SEGA forecast crosses the surveillance threshold within the 12-18 month window; prediction intervals calibrated (50%/90% coverage checked against held-out synthetic trajectories); no false alarms on stable patients.
- **TAND Agent** — detect the embedded signals authored into Patient B; zero spurious flags on patients with no authored TAND signal; output is briefing material, never a diagnosis string.
- **Therapeutics Strategist** — correct ClinicalTrials.gov matches for the cohort; every claim carries a source attribution; hedging language present where evidence is thin; < 3 min/case.
- **Orchestrator** — event router is deterministic (same event stream → same projection); incremental update touches only affected agents; a forced agent failure yields a "pending" surface state, never a silent gap.

#### 6.2.2 Integration tests

- ACMG synthesis (Opus) cross-checked against the deterministic combinatorial-rules reference for every cohort variant; any disagreement is logged and reviewed, not silently overridden.
- Phenome Mapper → Trajectory Modeler handoff: HPO timeline feeds the growth models correctly for the full cohort.
- All five agents → orchestrator → surface: a new synthetic imaging-report event for Patient B propagates through Trajectory, updates the in-visit dashboard quadrant, and appears in the pre-visit briefing "what's-new" within the orchestrator's incremental-update path.
- Provenance integrity: every surface claim drills through to a provenance record with model id/version, prompt template version, RAG source URIs, input hash, and latency.

#### 6.2.3 End-to-end dry runs & demo runbook acceptance

Three full dry runs against the freshly regenerated cohort (T-57, T-60). The runbook is accepted when all three of these hold on a clean regeneration:

1. **Act One** — Patient A mosaic recovery runs live, resolves to Likely Pathogenic, and the audit trail opens to inspection. (8 min.)
2. **Act Two** — Patient B 4-quadrant dashboard renders with current forecasts and TAND briefing signals; Patient C therapeutics brief renders with six sections and full attribution; the pre-visit briefing assembles on demand. (12 min.)
3. **Act Three** — infrastructure/cost/scaling narrative and the Apache 2.0 posture present cleanly. (5 min.)

Total runtime ≤ 30 min with 15 min discussion headroom, achieved on three consecutive runs without manual intervention beyond the scripted demo events.

### 6.3 Risk Register

| ID | Risk | Likelihood / Impact | Mitigation |
|----|------|--------------------|------------|
| R-01 | Synthetic notes/imaging are unconvincing to clinicians; the demo reads as toy data | Med / High | Generate from published TSC templates; clinician-sampled review by the TSC faculty lead in W8 (T-58); persistent synthetic watermark removes any pretense; de-scope cohort to 30 to raise per-patient quality if needed |
| R-02 | Mosaic recovery (Act One) fails at low VAF — the load-bearing demonstration breaks | Low / Critical | Mutect2 mosaic-aware calling validated against the 7 known authored variants in W2; BAMSurgeon VAF parameters tuned to 4-12%; Act One is never de-scoped; ddPCR-recommendation framing means the engine flags for confirmation rather than asserting certainty |
| R-03 | Opus ACMG synthesis disagrees with combinatorial rules or produces a false-positive Pathogenic | Med / High | Deterministic combinatorial-rules reference as the validator (T-19); every disagreement logged; output is explicitly a draft for molecular-geneticist sign-off, never autonomous (T-20) |
| R-04 | RunPod GPU availability or Parabricks throughput blocks variant calling | Med / Med | Cohort starts at BAM (no FASTQ alignment burden); calling is batchable overnight; Spark can run reduced-throughput calling locally as fallback; cohort de-scope 50→30 halves the workload |
| R-05 | TAND agent over-flags, eroding the "briefing not alert" discipline and clinician trust | Med / High | Deterministic scoring/aggregation layer downstream of the LLM discourse pass; zero-spurious-flag unit test (T-39); alert-surface recalibration trigger at > ~3 alerts/clinician/week; cluster-set de-scope available |
| R-06 | Eight weeks slips for a solo builder | Med / Med | Fixed de-scope order (6.1.7); Act One protected; weekly exit gates; featured patients authored early so the demo narrative is testable before the cohort is complete |
| R-07 | Overclaiming creep — institutional integrations (Epic/Clarity/LIMS, imaging AI) get described as built | Low / High | Every surface and the delivery doc explicitly marks these as architectural, not built; known-limitations doc in T-61; demo Act Three frames them as Phase-1 institutional work |
| R-08 | Clinical-collaboration dependency stalls (TSC faculty lead or Hagedorn's team unavailable in the window) | Med / Med | The synthetic demo is self-contained and requires no real data or IRB; clinician review (T-58) is a quality enhancer, not a blocker for the demo to run; engagement framed as offered, not assumed |

### 6.4 Dependencies

**Internal (reused HCLS AI Factory v1.3.0 substrate).** LangGraph orchestration, Milvus RAG (with the net-new TSC corpus partition), the tiered Claude client (Haiku/Sonnet/Opus via API) with Llama 3.1 70B Instruct via Ollama as local fallback, PostgreSQL + Redis + MinIO, the audit/provenance layer, and the `precision_oncology_agent` template that the TSC dir scaffolds from. Net-new and therefore on the critical path: the multi-agent orchestrator with shared event-sourced state, the synthetic cohort pipeline, and the five TSC agents.

**External tooling & data.** Synthea (MIT) + TSC modules; BAMSurgeon over NA12878-derived BAMs; NVIDIA Parabricks (RunPod GPU); snpEff/VEP. Reference resources: ClinVar, gnomAD v4, LOVD-TSC, dbSNP, HPO, SNOMED-CT, PubMed/PMC, ClinicalTrials.gov, FDA. Compute: my DGX Spark (GB10 Grace Blackwell, ~1,000 TOPS, 128 GB unified LPDDR5x, 4 TB NVMe) as the primary, with RunPod GPUs for burst variant calling, parallel cohort generation, and heavier local inference.

**Clinical collaboration (enhancers, not blockers).** Dr. Philip A. Hagedorn (CHIO, CCHMC) and the Division of Biomedical Informatics — the Marshall-Hagedorn discourse methodology that the TAND agent extends, plus output-surfacing review and sponsorship. The CCHMC TSC clinical & research faculty lead — note/imaging realism review and demo audience. None of these gate the synthetic demo running on the Spark; they shape its credibility and its path to a real Phase 1.

### 6.5 Roadmap

- **v0.1 — Synthetic demo (this build, target early Q3 2026).** The 50-patient (or de-scoped 30) deterministic cohort, all five agents, the deterministic orchestrator, the three clinician surfaces, and the three-act demo. Runs on Spark + RunPod. No real data, no IRB, no Epic/LIMS integration. This is what exists now.
- **v0.2 — Hardening & breadth.** Expand the synthetic cohort and the literature partition; widen the eval harness; add the deferred alert categories and the full TAND cluster set if de-scoped; tighten provenance drill-through and the regeneration tooling. Still fully synthetic.
- **v0.3 — Real-data validation prep (institutional, gated).** Define the IRB protocol; specify the Epic Clarity/Caboodle + biobank LIMS data-plumbing contracts (the Phenome Mapper / TAND / Trajectory ingestion that is explicitly not built in the demo); design a shadow-mode validation against curated real Variant Curator cases sourced from the Discover Together Biobank's banked tuber/AML/SEGA tissue. No autonomous output; molecular-geneticist sign-off retained. SaMD posture assessed here as institutional work, not asserted.
- **v1.0 — Institutional deployment.** Production ingestion within the CCHMC envelope (the Winslow Pavilion as infrastructure, the five CCHMC source areas feeding the engine), validated real-data agents under IRB, integrated clinician surfaces, and a replication pattern ("swap the box labels, keep the wiring") toward partner sites such as TGen and City of Hope. This is the destination, not the deliverable of this build.

### 6.6 Acceptance Criteria (AC-1 .. AC-9)

The build is accepted when all nine hold against a clean cohort regeneration:

| ID | Acceptance criterion |
|----|---------------------|
| AC-1 | Variant Curator recovers all 7 mosaic variants at VAF ≥ 5% with correct ACMG-AMP class and zero false-positive Pathogenic calls, < 5 min/case |
| AC-2 | Patient A resolves live to Likely Pathogenic (PVS1 + PM2 + PP4) with a ddPCR recommendation and an openable audit trail |
| AC-3 | Phenome Mapper achieves recall ≥ 90% and precision ≥ 85% with evidence spans, full cohort < 1 hr |
| AC-4 | Trajectory Modeler forecasts Patient B SEGA threshold crossing in the 12-18 month window with calibrated intervals and no false alarms |
| AC-5 | TAND agent surfaces Patient B's embedded signals as briefing material with zero spurious flags and no diagnostic language |
| AC-6 | Therapeutics Strategist produces a six-section brief with correct trial matches, appropriate hedging, and full source attribution, < 3 min/case |
| AC-7 | Orchestrator is deterministic, handles a forced agent failure with a "pending" surface state, and assembles all three surfaces on demand |
| AC-8 | Every surface claim drills through to a complete provenance record; synthetic watermark present on all surfaces and artifacts |
| AC-9 | Three consecutive end-to-end dry runs complete the three-act demo in ≤ 30 min with no manual intervention beyond scripted events, and the known-limitations doc correctly marks all unbuilt institutional integrations |

These nine criteria are the contract. Meeting them is what earns the conversation with Hagedorn's team that the final section frames; missing any of them is a reason to keep building rather than to demo. The closing section ties this build back to the institutional thesis it is meant to test.
