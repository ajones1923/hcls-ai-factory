---
name: demo-foundation-alignment
description: >-
  Anchors ALL future development of the HCLS AI Factory to the demonstration foundation — the
  D1–D7 demo portfolio and its capability-coverage matrix. Consult whenever planning, adding, or
  reviewing any engine, agent, model, or capability: every new capability must earn a demo home,
  update the coverage matrix, stay honestly scoped, and keep the shock-and-hope portfolio coherent.
  This is the bridge between the north star (why we build) and build-housekeeping-standards (how we
  build): it governs WHAT we build and how it is proven.
---

# Demo-Foundation Alignment — build to the demonstration portfolio

The demonstrations are not marketing afterthoughts — they are the **proving ground and the
coverage contract**. The master demo document is the single source of truth for the demonstration
program (`skills/tmp/HCLS_AI_Factory_Demos_Main_Document.docx`, local); this skill carries its
governing structure so development stays anchored to it even where the doc isn't present.

**The rule in one line:** nothing ships that cannot be shown, honestly, inside the demo portfolio.
Every new capability must (1) earn a demo home, (2) update the coverage matrix, (3) be honestly
labeled, and (4) keep the portfolio clinically real.

## The demo foundation (what everything maps to)

Coverage is achieved by a **portfolio**, never by one demo doing everything (that would break
clinical realism and violate the output-honesty gate). The portfolio:

| ID | Demonstration | Role | Primary reach |
|---|---|---|---|
| **D3** | **Tuberous Sclerosis Complex** | **Flagship** | The whole factory converging on one infant — every engine, most agents, the platform layer, one governed afternoon |
| D1 | Secondary Genomics → Novel Drug Molecule Generation | Core | Engines 1→2→3, MLOps, Workflow Composer, Chai-1 — "raw reads to a validated lead molecule" |
| D2 | Clinical Imaging → Cardiac CAC → Statin PGx → Cardiometabolic/Longevity | Core | Engines 4 + 6, Pharmacogenomics + Biomarker agents, GWAS — the cross-modality flagship |
| D4 | Pediatric Molecular Tumor Board | Completion | Engine 5, E1 germline, E8, Rare-Disease / Biomarker / Trial agents (fusion-first) |
| D5 | Oncology / CAR-T & Biologic Design | Completion | CAR-T agent, Engine 7 design suite, Chai-2, single-cell |
| D6 | Autoimmune Program | Completion | Precision Autoimmune agent, HLA, Pharmacogenomics |
| D7 | Neurology Suite (Parkinson's featured) | Completion | Neurology agent: Parkinson's, Alzheimer's, ALS, MS, stroke |

Coverage status: **Engines 1–8, all eight agents, every platform/infra item, and the frontier
models (Chai-1 in D1/D3/D5/D7, Chai-2 in D5, RFdiffusion + Evo 2 on-demand) each have ≥1 demo
home.** Completion demos exist to guarantee coverage and serve specialist audiences; they need not
all be shown live.

## The principles that govern development

1. **Portfolio, not one demo.** Spread coverage across demos. If a new capability would force one
   demo to "do everything," it belongs in a different (or new) demo — protecting realism and the
   honesty gate.
2. **Shock-and-hope.** The unit of impact is the **disease program** (one patient, whole factory),
   built in three beats — **weight → compression → hope**. Build toward capabilities that deepen
   that arc, not isolated party tricks.
3. **Honesty by construction.** Everything carries an honest status. A `live` capability is never
   `mock`-served; preclinical, roadmap, and research-use items are labeled as such; all clinical
   output is **decision support for a qualified clinician — never autonomous diagnosis or
   prescribing.** (See the Honesty ledger below.)
4. **Doorway principle.** Every focused demo is a door into the whole factory. A new capability
   should either open a doorway wider or let a spotlight recruit more of the platform.
5. **Clinical realism, verifiably.** Real tools, real guidelines, claims traceable to primary or
   authoritative sources. If a clinical claim can't be cited, it doesn't go on screen.

## The rule: every capability earns a demo home

When you plan, add, or modify any capability, answer all of these before it is "done":

- **Which demo(s) does it serve?** Extend an existing demo, or justify a new one. If it maps to no
  demo, either it's premature or a demo is missing — resolve that first.
- **Update the coverage matrix.** Add the capability → demo mapping so the "every capability has a
  demo home" invariant holds. Treat a coverage gap like a failing test.
- **Assign an honest status.** live / planned / preclinical / research-use / roadmap / gated — and
  add it to the Honesty ledger if it's anything other than plainly live.
- **Register + govern it.** Add it to the capability registry (the CI drift-guard
  `scripts/validate_registry.py` enforces this) and front it with the governance gates
  (`create_governed_app` / `install_governance`). See `build-housekeeping-standards`.
- **Check portfolio coherence.** Does it keep demos clinically real and appropriately scoped, or
  does it bloat one demo past believability?

## Alignment workflow (checklist)

1. **Start from a demo need**, not a capability in isolation — "which patient story does this make
   more real or more complete?"
2. Pick its **home demo(s)** from D1–D7 (or propose a new demo with a weight→compression→hope arc).
3. Build it to `build-housekeeping-standards` (governed app, tests, neutral, registered).
4. **Register the capability** and add its row to the coverage matrix.
5. **Label its honesty status**; if not plainly live, add it to the ledger and make sure any demo
   script says so out loud.
6. **Verify the clinical/scientific claims** against primary sources; keep the citation with the
   capability.
7. Update the master demo document (and, if kept in-repo, its markdown mirror) so the foundation
   stays the single source of truth.

## The honesty ledger (carry forward, state plainly wherever relevant)

- Gene therapy for TSC1/TSC2 (and gene-correction generally): **preclinical** — the factory is the
  open design/analysis bench, not a cure today.
- MAISI synthetic imaging: **research / augmentation / QA use, never a diagnostic source.**
- Single-cell atlas similarity and foundation-model cell embeddings: **roadmap.**
- Chai-2 de novo binder/antibody design: **gated / partnership access; integration contingent.**
- α-synuclein SAA, plasma p-tau217, NSD-ISS / SynNeurGe staging: **research- / trial-use biomarker
  inputs and frameworks the agents reason over — not routine clinical diagnostics.**
- Pediatric RNA fusion detection (Arriba / STAR-Fusion): **recommended near-term registry
  addition** to make the pediatric MTB (D4) fully airtight for a St. Jude / Cincinnati audience.
- "One box, elastic burst": heavy / ARM-incompatible models (Chai-1/2, RFdiffusion, Evo 2) burst
  to RunPod over a private Tailscale mesh — **say "elastic burst," never imply everything runs
  locally.**
- All clinical outputs: **decision support, not autonomous diagnosis or prescribing.**

## Do / Don't

**Do:** start development from a demo/patient-story need; map every capability to a demo home;
keep the coverage matrix green; label honesty status out loud; cite clinical claims; distinguish
live vs. representative vs. burst in any demo.

**Don't:** build capabilities with no demo home; overload one demo to look comprehensive; show a
roadmap/preclinical item as if it ships; imply "all on one box" when a model bursts to RunPod;
present any output as autonomous diagnosis; let a clinical claim on screen go uncited.

## Suggested arc for a room (reference)

Open on **D3** (whole factory, one infant — weight→compression→hope) → prove the engine with
**D1** (genome → novel molecule) → prove cross-modality with **D2** → go deep for the audience
(**D4** oncologists, **D7** neurologists, **D5** protein scientists) → close on the mission (hand a
discovered molecule to the trial agent; return to the child; their data makes the next child's
afternoon faster via the biobank).

## Related
- `hcls-core-vision-mission` — **why** we build (the north star the demos dramatize).
- `build-housekeeping-standards` — **how** we build (governance, registry, tests, neutrality).
- `09-ai-orchestration` — the registry + composer that make a capability demo-able.
- `08-inference-serving` — the live-never-mock honesty rule the demos depend on.
- `11-security-and-secrets` — the governance/honesty gates every demo run passes.
