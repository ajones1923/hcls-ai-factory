# Slide 19 — The logical architecture (walk-through diagram)

*Placed after Slide 18. A dense, technical "systems diagram" of the whole TSC AI Engine — modeled on
the HCLS AI Factory Pipeline MVP reference layout (left inputs → labeled processing lanes → right
consumers, with a spec/metrics badge strip along the bottom), rendered in the deck's dark house style.*

## 🎙️ Narration script (~85 sec — the full walk-through)

Let me put the entire engine on one diagram and walk you through it, left to right. [beat]

On the **far left** are the inputs — everything the engine reads. A child's genomic data, the tissue
specimen, their longitudinal records, the clinical notes, the imaging. And beneath those, the
reference knowledge it grounds against: the human phenotype ontology, the ClinVar and AlphaMissense
variant databases, the ACMG classification rules, and the trial and FDA sources.

That feeds the engine itself — four stages. **Stage one, Ground:** the Phenome Mapper runs first and
builds the patient's full, ontology-checked clinical picture; everything downstream depends on it.
**Stage two, Analyze** — three agents work in parallel: the Variant Curator recovers the mosaic
diagnosis from tissue and classifies it by the rulebook; the Trajectory Modeler — classical
statistics, not a chatbot — forecasts each growth with honest intervals; and the TAND agent surfaces
the hidden neuro-behavioral signals. **Stage three, Synthesize:** those three converge into the
Therapeutics Strategist, which assembles a single sourced options brief.

Running *beneath* all of it is the **orchestrator** — the deterministic conductor. It's event-sourced,
so every step is logged and replayable, and it routes each model call to the right tier — a small,
fast model for simple work, the most capable one for the hard reasoning, with a local fallback if the
cloud is unavailable.

**Stage four, Deliver:** the results flow out to the clinician surfaces on the right — the pre-visit
briefing, the in-visit dashboard, the disciplined alerts, and the three-D digital twin — and on to the
people who use them: the clinician, the geneticist, the care team, and the family.

And along the bottom, the foundation under everything: provenance on every output, a cost ledger, and
the validation scorecard — all on one NVIDIA DGX Spark, with the twin rendered on RunPod. [beat]
That's the whole engine, end to end.

---

## 🖼️ Nano Banana Pro prompt

**HOUSE STYLE — match the deck (this is the reference-architecture frame, index "19").** 16:9 landscape, but **dense and technical like an engineering systems diagram** — this is the most detailed frame in the deck. Background deep navy #0F1830; processing blocks #16213F with thin #26345C borders; primary accent sky-blue #7FB0FF; key numbers, ribbons & the foundation agent in gold #FFD166; body text off-white #E8EEF9. Small mono-style labels are welcome. Clean modern geometric sans-serif. Small pink "SYNTHETIC · not FDA-cleared" pill top-right; a thin footer band; a small gold **"19"** index in the bottom-right corner. Render ALL text crisply, legibly, and spelled exactly as written — no duplicated blocks or labels. This frame has many small labels: keep each one sharp; put exact numbers in text.

Create **FRAME 19 — THE TSC AI ENGINE LOGICAL ARCHITECTURE**: a comprehensive left-to-right **swimlane systems diagram** (in the spirit of a reference "pipeline MVP" architecture poster) — INPUT DATA SOURCES down the far left, FOUR labeled processing stages flowing across the middle, the ORCHESTRATOR as a horizontal spine beneath the stages, CLINICIAN SURFACES + CONSUMERS down the far right, and a HARDWARE / METRICS / OPEN-SOURCE badge strip along the very bottom. Thin glowing connector arrows show the data flow; a small color-coded **legend** sits top-right of the header.

**TOP HEADER BAND (full width):** Title "TSC AI Engine — Logical Architecture", subtitle "Engine 7 · HCLS AI Factory · Five Agents · One Deterministic Orchestrator · on a single NVIDIA DGX Spark". Small **LEGEND** chips: "▮ Ground (sky-blue)  ▮ Analyze (sky-blue)  ▮ Synthesize (gold)  ▮ Deliver (teal)  — — Data flow  ⟶ Orchestrated call". A small top-left use-case tag: "Synthetic 50-patient TSC cohort · featured cases A / B / C".

**FAR-LEFT COLUMN — "INPUT DATA SOURCES"** (a vertical stack of small icon+label tiles, each once):
- "Genomic VCF (mosaic-aware)" (DNA icon)
- "Tissue specimen (biobank)" (specimen-tube icon)
- "Longitudinal records — HPO" (timeline icon)
- "Clinical notes" (document icon)
- "Imaging reports" (scan icon)
- a thin divider, then **"REFERENCE KNOWLEDGE"**: "HPO ontology — 19,389 terms" · "ClinVar + AlphaMissense" · "ACMG-AMP rules" · "Clinical trials + FDA".

**CENTER — FOUR PROCESSING STAGES, left to right, each a titled lane:**
- **STAGE 1 · GROUND** (gold "RUNS FIRST · FOUNDATION" ribbon): block **"TSC-Phenome Mapper"** — "ontology-grounded HPO profile + surveillance-gap detection." Tag "Sonnet / Haiku".
- **STAGE 2 · ANALYZE** (three parallel blocks stacked): **"TSC-Variant Curator"** — "mosaic-aware ACMG-AMP; recovers low-VAF variants blood tests miss" tag "Opus + deterministic classifier"; **"TSC-Trajectory Modeler"** — "classical forecasting — SEGA / AML / eGFR / seizures, 50/90% intervals" tag "classical · not an LLM"; **"TAND Surveillance Agent"** — "neuropsychiatric signals via diagnostic-uncertainty discourse-marker NLP; briefing only" tag "Sonnet / Opus".
- **STAGE 3 · SYNTHESIZE** (gold block): **"TSC-Therapeutics Strategist"** — "six-section options brief: trials · evidence · FDA actions · all source-attributed" tag "Opus". (The three Analyze blocks converge with arrows into this one.)
- **STAGE 4 · DELIVER** (teal lane) — leads into the surfaces column.

**HORIZONTAL SPINE beneath stages 1–3 — "TSC ORCHESTRATOR":** "deterministic, NOT an LLM · event-sourced · 13 event types · dependency-ordered dispatch · materialized projections · conservative failure." A small sub-bar on the spine: **"MODEL TIER ROUTER — Haiku (fast) · Sonnet (analysis) · Opus (high-stakes) · Local model (offline fallback)."** Thin "⟶ orchestrated call" arrows rise from the spine to each agent.

**FAR-RIGHT COLUMN — "CLINICIAN SURFACES" then "CONSUMERS":** surfaces (teal-tinted tiles): "Pre-Visit Briefing" · "In-Visit Dashboard (variant · phenotype · trajectory · TAND)" · "Async Alerts (≤ 3 / clinician / week)" · gold-tagged **"Anatomical Digital Twin (3-D · OpenUSD · RTX)"**. Then arrows to consumer icons: "Clinician" · "Geneticist" · "Care team" · "Family".

**BOTTOM BADGE STRIP (full width, three segments):**
- **TRUST LAYER:** "Provenance on every output · Cost ledger · Validation harness".
- **VALIDATION SCORECARD (gold):** "100% classification accuracy · +12 pts detection · 6/6 mosaics recovered · 100% traceable".
- **PLATFORM:** "1× NVIDIA DGX Spark (~$4,699) · twin rendered on RunPod (RTX) · OpenUSD + NVIDIA Omniverse · Apache 2.0 open source".

**FOOTER (small, muted):** "SYNTHETIC demonstration data · decision support, clinician review required · construct validity, not yet clinical validity · not FDA-cleared."

A dense but beautifully organized systems diagram: inputs at far left, four labeled stages flowing right, the orchestrator spine beneath, surfaces and consumers at far right, the metrics/hardware badge strip across the bottom — every block, label, and number crisp and exactly as written.
