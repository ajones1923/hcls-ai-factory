# The TSC AI Engine — Master Study Guide

### A visual, module-by-module learning path to understand the whole system end to end

**Author:** Adam M. Jones — HCLS AI Factory, TSC Intelligence Engine (Engine 7)
**Purpose:** a self-study companion that ties together every infographic and document into one ordered
learning journey — so you can build a complete, correct mental model of the engine across biology,
software, statistics, visualization, validation, and strategy.
**License:** Apache 2.0 · 2026

> **How to use this guide.** It is organized as **eleven modules** that build on one another. Each module
> follows the same loop:
>
> 1. **🎯 Objective** — what you will understand by the end.
> 2. **📊 Look** — the infographic to study, and *what to look for in it.*
> 3. **📖 Read** — the in-depth explanation (this is the core of each module).
> 4. **🔑 Key terms** — the vocabulary to lock in.
> 5. **✅ Check yourself** — questions you should be able to answer before moving on.
> 6. **📚 Go deeper** — the full document to read if you want the fullest treatment.
>
> Work in order the first time. The infographics live in `docs/infographics/` (generate the images from
> the Nano Banana Pro prompts, or view ones you have already made); the deep documents live in `docs/`.
> Everything describes **synthetic** demonstration data — the engine is clinician decision support and is
> **not FDA-cleared.**

---

## Module map (the whole journey at a glance)

*(Infographics are the numbered prompt files in `docs/infographics/`, `01`–`18`.)*

| # | Module | Infographic(s) to study | Deep-dive doc |
| --- | --- | --- | --- |
| 0 | Orientation | `01` executive one-pager (start here) | — |
| 1 | The disease: TSC biology & clinical reality | `02` biology_mtor + `03` problem_solution | LEARN_GUIDE_ADVANCED, Parts 1–2 |
| 2 | The engine at a glance | `04` engine7_hero, `05` engine7_portrait | LEARN_GUIDE, Part 3 |
| 3 | The five agents in depth | `05` engine7_portrait | LEARN_GUIDE_ADVANCED, Parts 5–9 |
| 4 | Flagship: mosaic recovery | `09` mosaic_recovery | LEARN_GUIDE_ADVANCED, Part 6 |
| 5 | Forecasting with honest uncertainty | `12` omniverse_hero | LEARN_GUIDE_ADVANCED, Part 7 |
| 6 | The orchestrator & the architecture | `06` engine7_landscape + `07` ai_tiering | LEARN_GUIDE_ADVANCED, Part 4 |
| 7 | The 3-D digital twin (Omniverse) | `10` four_scenes, `11` pipeline | OMNIVERSE_VISUALS_GUIDE |
| 8 | Trust: how we know it works | `13` validation_scorecard + `14` trust_layer | LEARN_GUIDE_ADVANCED, Part 11 |
| 9 | The live demo, end to end | `15` engine7_demo_flow | TSC_DEMO_RUNBOOK |
| 10 | Strategy & the road ahead | `16` cincinnati_roadmap + `17` generalization + `18` small_cheap_open | CCHMC_HANDOVER; LEARN_GUIDE_ADVANCED, Parts 13–14 |
| 11 | Synthesis & mastery check | `08` patient_journey, `01` executive one-pager (+ all) | — |

---

## Module 1 — The disease: TSC biology and clinical reality

**🎯 Objective.** Understand what Tuberous Sclerosis Complex is at the molecular level, why it affects many
organs, what mosaicism is, and the six reasons caring for it is hard today.

**📊 Look — `biology_mtor_mechanism`, then `story_problem_solution`.** Start with the biology poster: trace
genes → the hamartin–tuberin *brake* → the normal-green vs disease-red mTORC1 comparison → the body map →
the mosaicism panel. Then read the *left column only* of `story_problem_solution` (the six struggles). The
first image teaches you *why the disease happens*; the second teaches you *why it is hard to manage.* Hold
both in mind as you go.

**📖 Read.** TSC is caused by loss-of-function changes in **TSC1** or **TSC2** — genes whose proteins
(**hamartin** and **tuberin**) form a complex that acts as a *brake* on a master growth pathway called
**mTORC1**. Break the gene, release the brake, and cells grow when they should not, producing benign
growths (**hamartomas**) across the body: brain tubers and the **SEGA** tumor, kidney **AMLs**, skin
lesions, cardiac and pulmonary growths. Because growths sit in specific places and change over time, TSC
care is dominated by *surveillance* (scheduled monitoring) and by *timing decisions* (when to intervene).

The single most important biological idea for this engine is **mosaicism**. A mutation that occurs *after*
fertilization ends up in only some of the body's cells. The fraction carrying it is the **variant allele
fraction (VAF)**. In blood, that fraction can be so low that standard testing sees nothing — which is why
roughly **10–15%** of clinically diagnosed TSC patients are reported "no mutation identified." Crucially,
the variant is usually *enriched in affected tissue* (e.g., a removed brain tuber), where it can be
recovered if someone looks. Keep that fact; it is the hinge of Module 4.

The six clinical struggles to internalize: (1) the **diagnostic odyssey** from missed mosaics; (2)
**fragmented** multidisciplinary care; (3) hard **timing** decisions; (4) under-recognized neuropsychiatric
disease (**TAND**); (5) **alert fatigue**; (6) **scattered, untrusted data**.

**🔑 Key terms.** TSC1/TSC2, hamartin/tuberin, mTORC1, hamartoma, SEGA, AML, mosaicism, VAF, NMI, TAND,
surveillance.

**✅ Check yourself.**
- In one sentence, why does a broken TSC gene cause growths?
- Why does a blood test miss a mosaic variant that tissue testing can find?
- Name four of the six clinical struggles.

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE_ADVANCED`, Parts 1–2.

---

## Module 2 — The engine at a glance

**🎯 Objective.** Form the big-picture mental model: five agents + one orchestrator + clinician surfaces,
running on one small computer.

**📊 Look — `engine7_hero`, then `engine7_portrait`.** The hero gives you the five agent names in a flow;
the portrait shows each agent with its one-line job, the orchestrator spine, the trust layer, and the
three surfaces. *Look for the left-to-right flow:* inputs → Phenome Mapper (first) → the other agents →
orchestrator → surfaces.

**📖 Read.** Picture a small team of five expert assistants (**agents**), a conductor that keeps them
organized (the **orchestrator**), and a few clean screens the clinician actually uses (**surfaces**). The
five agents are the **Phenome Mapper** (builds the phenotype map, runs first), the **Variant Curator**
(genetics + mosaic recovery), the **Trajectory Modeler** (forecasting), the **TAND Surveillance Agent**
(neuropsychiatric signals), and the **Therapeutics Strategist** (the integrated options brief). The whole
system runs on a single **NVIDIA DGX Spark** (~$4,699), is **open source (Apache 2.0)**, and is built to
be **honest by design** — traceable, uncertainty-aware, and clearly labeled as a demonstration.

The most important structural idea: **determinism where correctness matters, AI only at the edges.** The
genetic verdict, the coordination logic, and the numerical forecasting are fixed, auditable code. Large
language models are used only to read free text and write explanations — and never to overrule a
deterministic result.

**🔑 Key terms.** Agent, orchestrator, surface, deterministic vs. LLM, DGX Spark, Apache 2.0.

**✅ Check yourself.**
- Name the five agents and which one runs first (and why).
- Which parts of the engine are *not* AI, and why is that deliberate?

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE`, Parts 3–5.

---

## Module 3 — The five agents in depth

**🎯 Objective.** For each agent, be able to state its *job, the problem it solves, and how it works.*

**📊 Look — `engine7_portrait`.** Study the five agent cards one at a time. Note the small tier tags
(Haiku/Sonnet/Opus) and that the Trajectory Modeler is marked "not an LLM."

**📖 Read.**

- **Phenome Mapper.** Builds the ontology-grounded phenotype map. It assigns each finding an official
  **HPO** code and *validates every code against the real 19,389-term HPO* — correcting labels, remapping
  obsolete IDs, and dropping codes that are not real terms. It records each finding's exact text span,
  **polarity** (present vs. negated), and **temporality**, and logs **discordances** (e.g., a feature
  present at onset but negated now).
- **Variant Curator.** Recovers low-VAF mosaic variants from tissue and classifies them with the
  **ACMG-AMP** rulebook. A deterministic classifier makes the verdict; an LLM only explains it. It must
  *reject* strand-biased artifacts, which makes its "zero false positives" claim real.
- **Trajectory Modeler.** Classical statistics (not an LLM) forecasting SEGA, AML, kidney function, and
  seizures with **prediction intervals** (a "cone of uncertainty"), crossing probabilities, and
  surveillance-cadence advice.
- **TAND Surveillance Agent.** Finds hedged neuropsychiatric signals using a discourse-marker method, scores
  them with weights/recency/recurrence, ignores *negated* concerns, and surfaces findings as gentle
  briefing material — never an alarm.
- **Therapeutics Strategist.** Integrates everything plus trials, literature (via **RAG**), and FDA actions
  into a six-section, fully **source-attributed** options brief.

**🔑 Key terms.** HPO, polarity/temporality, discordance, ACMG-AMP, strand bias, prediction interval,
discourse marker, RAG, source attribution.

**✅ Check yourself.**
- Which agent recovers mosaic diagnoses, and what rulebook does it use?
- Why is the Trajectory Modeler deliberately *not* an LLM?
- How does the TAND agent avoid false positives?

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE_ADVANCED`, Parts 5–9.

---

## Module 4 — Flagship capability: mosaic recovery and the diagnostic odyssey

**🎯 Objective.** Deeply understand the engine's headline win — recovering diagnoses blood tests miss.

**📊 Look — `story_mosaic_recovery`.** Follow the gray→gold arc: blood test (negative) → mosaicism →
tissue → recovery → diagnosis. Note the central image: a cell field where ~1 in 12 cells glows.

**📖 Read.** A child clearly has TSC, but the blood test reports "no mutation identified." The reason is
**mosaicism** at a low **VAF** — too few variant-carrying cells in blood to detect. The variant is
concentrated in affected tissue (a resected tuber), where the VAF is much higher. The Variant Curator runs
**mosaic-aware** analysis on that tissue, finds the low-VAF variant (e.g., a *TSC2* change at 8.3% VAF —
about one cell in twelve), and classifies it under ACMG. Two rigor points: the verdict comes from a fixed
rulebook (with **PVS1** down-weighted for low-VAF mosaics, yielding "Likely Pathogenic"); and the agent
must reject a planted strand-biased **artifact**, so "zero false positives" is *earned*. Measured result on
the synthetic cohort: variant **detection 86% → 98% (+12 points)**, with **6 of 6** sub-threshold mosaics
recovered. This is the capability a tissue **biobank** unlocks — and the basis of the proposed first real
study (Module 10).

**🔑 Key terms.** Mosaic recovery, VAF, PVS1 strength modulation, artifact rejection, detection uplift.

**✅ Check yourself.**
- Walk the path from "negative blood test" to "confirmed diagnosis."
- What makes the "zero false positives" claim meaningful rather than automatic?
- What does "+12-point detection uplift" actually count?

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE_ADVANCED`, Part 6.

---

## Module 5 — Forecasting with honest uncertainty

**🎯 Objective.** Understand how the engine predicts the future *and how it expresses confidence.*

**📊 Look — `omniverse_hero`.** This single image is the best teacher here: a growing tumor wrapped in a
glass cloud whose *size is the uncertainty range.* Burn in the idea: bigger cloud = less certain.

**📖 Read.** Forecasting numbers with calibrated confidence is a job for **classical statistics**, not a
chatbot — which is why this agent is deterministic. With only a few past measurements, an individual's
trend is unreliable, so the engine borrows strength from a **mixed-effects population model** (how SEGAs
typically grow) and **shrinks** the individual estimate toward that prior in proportion to how trustworthy
each is — a **Bayesian** update. A **Gaussian process** then sets the *width* of the forecast, producing
**50% and 90% prediction intervals** — the "cone of uncertainty," wide when data is sparse, tight when it
is clear. From that distribution the engine computes a **probability of crossing** the clinical threshold
by each horizon, grades the risk **likely/possible/unlikely**, and recommends checking *more often* when
risk is elevated. It does this for four quantities at once (SEGA, AML, eGFR, seizures).

**🔑 Key terms.** Mixed-effects model, Bayesian shrinkage, Gaussian process, prediction interval, crossing
probability, surveillance cadence.

**✅ Check yourself.**
- Why does the engine borrow from a population model instead of using only the patient's own data?
- What does the *size* of the uncertainty cloud tell a clinician?
- Why is this agent not an LLM?

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE_ADVANCED`, Part 7.

---

## Module 6 — The orchestrator and the architecture

**🎯 Objective.** Understand how the pieces are coordinated and why the design is trustworthy.

**📊 Look — `engine7_landscape`, then `architecture_ai_tiering`.** The landscape view traces the overall
flow (inputs → agents in dependency order → orchestrator spine → surfaces). The AI-tiering view shows *how
the AI is used responsibly*: a deterministic core (green) that an LLM never overrules, and the LLM tiers
(Haiku → Sonnet → Opus → local Llama) routed by stakes, with a cost ledger.

**📖 Read.** The **orchestrator** is a deterministic conductor (not an AI). It runs agents in **dependency
order** — Phenome Mapper first because the others build on it — and records every step in an **append-only
event log** (event sourcing) with **13 event types**. The state a clinician sees is *computed by replaying
that log*, which makes the system reproducible and auditable: nothing is silently changed, and the whole
history can be reconstructed. If an agent fails, the orchestrator marks that section "pending" rather than
showing stale data — **conservative, fail-safe**. It also supports **incremental updates** (a new note
re-runs only the note-dependent agents). When LLMs are used, a **tiered router** sends each step to the
right-cost model (Haiku/Sonnet/Opus, with a local fallback) and a **cost ledger** tracks spend. Production
backends (database, cache, storage, vector search) each have a lightweight **fallback**, so the whole
system runs on a bare single computer offline.

**🔑 Key terms.** Event sourcing, append-only log, projection, dependency order, conservative failure,
incremental update, model tiering, cost ledger, runs-here fallback.

**✅ Check yourself.**
- What is event sourcing, and why does it make the engine auditable?
- What happens when an agent fails — and why is that the safe choice?

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE_ADVANCED`, Part 4.

---

## Module 7 — The 3-D digital twin (Omniverse)

**🎯 Objective.** Understand what each 3-D scene shows, how the visuals are built, and the honesty designed
into them.

**📊 Look — `omniverse_four_scenes`, then `omniverse_pipeline`.** From the four-scenes poster, learn what
each scene shows; from the pipeline, learn the **author-on-Spark / render-on-RunPod** split.

**📖 Read.** A **digital twin** here is a *deterministic 3-D rendering of the engine's validated results* —
not a simulation and not a patient's real scan (the anatomy is a *stylized atlas* scaled by real
measurements). Four scenes: the **lesion-trajectory twin** (a lesion grows along its forecast inside an
uncertainty envelope, crossing a threshold membrane); the **mosaic "powers-of-ten"** scene (a cell field
where exactly the recovered VAF fraction glows); the **whole-child atlas** (organs lit by the phenome); and
the **population view** (50 patients, the 7 recoveries ringed in gold). They are written in **OpenUSD**,
rendered with **RTX** in **NVIDIA Omniverse**, and dressed in **MDL** materials (glass envelopes, glowing
cells). The defining honesty choice: **the envelope's radii equal the engine's 50%/90% prediction
intervals** — so the picture cannot imply more confidence than the math supports. And the workflow is
split: the Spark *authors* the scene (CPU, free, deterministic); a rented cloud GPU (**RunPod**) *renders*
it.

**🔑 Key terms.** OpenUSD, prim, time-sample, point instancer, MDL, RTX, uncertainty envelope, atlas vs.
patient-specific, authoring/render split.

**✅ Check yourself.**
- What does each of the four scenes show?
- Why can the trajectory twin's render "not lie" about confidence?
- Which machine authors the scene, and which renders it?

**📚 Go deeper.** `TSC_ENGINE_OMNIVERSE_VISUALS_GUIDE` (the full visual companion).

---

## Module 8 — Trust: how we know it works

**🎯 Objective.** Understand the validation results *and* the crucial limit on what they prove.

**📊 Look — `story_validation_scorecard`, then `story_trust_layer`.** The scorecard gives the measured
numbers (read the caveat band twice — it is the most important part). The trust-layer poster shows the five
honesty disciplines that surround those numbers: provenance, visible uncertainty, the human sign-off gate,
conservative failure, and honest construct-validity labeling.

**📖 Read.** Because the demo cohort is synthetic, the correct answers are known in advance, so the engine
can be graded against them. Results: **100%** classification accuracy, **+12-point** detection uplift,
**6/6** mosaics recovered, **0** null-variants-called-benign, **100%** traceable. The engine also separates
two honest ideas: **detection** (did it find the variant at all?) versus **molecular diagnosis** (detection
*plus* a returnable Pathogenic/Likely-Pathogenic call). The decisive caveat: these prove **construct
validity** — the logic works on data with known answers — **not clinical validity**, which requires real
patients. The flight-simulator analogy: rigorous proof the systems work, but not the first real flight.
This honesty is not a weakness to hide; it is what makes the work credible and is ethically required.

**🔑 Key terms.** Construct validity, clinical validity, detection vs. molecular-diagnosis yield, provenance
completeness, ground truth.

**✅ Check yourself.**
- What is the difference between construct validity and clinical validity here?
- Why separate "detection" from "molecular diagnosis"?
- Why is stating the caveat a *strength*?

**📚 Go deeper.** `TSC_ENGINE_LEARN_GUIDE_ADVANCED`, Part 11.

---

## Module 9 — The live demonstration, end to end

**🎯 Objective.** Understand how the parts come together as a three-act story.

**📊 Look — `engine7_demo_flow`.** Trace the numbered nodes: pre-flight → Act One → Act Two → Act Three →
the close, with the author/render rail underneath.

**📖 Read.** **Act One (mosaic recovery):** show the recovered variant and the glowing cell field (the
diagnostic-odyssey win). **Act Two (the longitudinal twin):** scrub the lesion twin from −24 to +18 months,
watching it grow inside its uncertainty envelope toward the threshold. **Act Three (scale):** the
population view, the validation scorecard, the $0 cost — one Spark, a whole disease. The narrative mirrors
the *evidence*: a single dramatic recovery, a single calibrated forecast, then the whole cohort. The
guiding line: *lead with the tissue, show the glass, tell the truth.*

**🔑 Key terms.** Three-act arc, pre-flight, fly-in, the demo's honest close.

**✅ Check yourself.**
- What is shown in each act, and why that order?
- What is the one-line ethos of the demo?

**📚 Go deeper.** `TSC_DEMO_RUNBOOK` (the presenter's script, with live numbers).

---

## Module 10 — Strategy and the road ahead

**🎯 Objective.** Understand the two-phase plan, the Cincinnati touchpoints, and why the biobank is the
wedge.

**📊 Look — `story_cincinnati_roadmap`, then `story_generalization` and `story_small_cheap_open`.** The
roadmap shows the Cincinnati path (Phase 1 → bridge study → Phase 2's five touchpoints, biobank starred).
The generalization poster shows the *reusable-engine* thesis (TSC today → NF1/NF2/Rett/Williams tomorrow,
"swap the box labels"); the small/cheap/open poster shows the *positioning* thesis (one ~$4,699 Spark +
RunPod, Apache 2.0, vs. closed hyperscale).

**📖 Read.** **Phase 1** is the complete synthetic demo on the Spark (+ RunPod render) — the version that
earns the meeting. **Phase 2** connects to a real institution. For Cincinnati Children's, five touchpoints
map onto the engine: the **Discover Together Biobank** (banked tissue — the substrate, and the
highest-leverage piece), the **Winslow Pavilion** (infrastructure), the **Division of Biomedical
Informatics / Dr. Hagedorn** (methodology + sponsor), the **TSC clinic** (patients), and **Epic + the
LIMS** (data plumbing). The smart first real study is a **retrospective "no-mutation-identified"
re-analysis**: run the Variant Curator on *already-banked, already-consented* tissue from patients told
years ago that nothing was found. If it recovers missed diagnoses, that is a concrete, publishable result —
no new surgery, no new wait. It is the fastest credible route from construct validity to a real clinical
result. The same architecture also **generalizes** to other rare diseases (NF1, NF2, Rett, Williams,
mTORopathies) by swapping disease-specific pieces — "swap the box labels, keep the wiring."

**🔑 Key terms.** Phase 1 vs. Phase 2, the five touchpoints, retrospective NMI study, biobank wedge,
generalization.

**✅ Check yourself.**
- Why is the biobank the highest-leverage touchpoint?
- Why is a *retrospective* study the smart first step?
- What turns this from "a demo" into "a result"?

**📚 Go deeper.** `TSC_INTELLIGENCE_ENGINE_CCHMC_HANDOVER`; `LEARN_GUIDE_ADVANCED`, Parts 13–14.

---

## Module 11 — Synthesis and mastery check

**🎯 Objective.** Confirm you can explain the whole system, in your own words, to anyone.

**📊 Look — `story_patient_journey`, then `summary_executive_onepager`.** These are your two synthesis
images. The patient-journey poster walks one child's data through every agent to the clinical decision —
trace it and narrate each stage from memory. The executive one-pager compresses the entire system to a
single page (problem · what it is · the proof · where it goes); if you can reconstruct that one page
unaided, you have the whole thing.

**📖 Read & reflect.** You now hold the full chain: a **disease** driven by a released growth brake, whose
care is hard for six concrete reasons → a **five-agent engine** that grounds phenotype, recovers mosaic
diagnoses, forecasts with honest uncertainty, surfaces hidden neuropsychiatric signals, and integrates
options → coordinated by a **deterministic, auditable** orchestrator → presented in clean surfaces and a
**3-D twin** whose beauty cannot outrun its confidence → **measured** against known answers with the
honest caveat stated aloud → demonstrated as a **three-act story** → aimed at a **two-phase strategy**
whose wedge is a tissue biobank.

**✅ Final mastery checklist — can you, without notes:**
- [ ] Explain why a broken TSC gene causes growths, and what mosaicism is.
- [ ] Name the five agents and each one's job.
- [ ] Walk the mosaic-recovery path from negative blood test to confirmed diagnosis.
- [ ] Explain what the uncertainty envelope is and why it keeps the visuals honest.
- [ ] Describe event sourcing and why it makes the engine auditable.
- [ ] State the four headline numbers *and* the construct-validity caveat.
- [ ] Recount the three demo acts in order.
- [ ] Explain Phase 1 vs. Phase 2 and why the biobank is the wedge.
- [ ] Say, in one sentence each, what this engine *is* and what it is *not yet.*

**Teach it to lock it in.** The fastest way to confirm 100% understanding is to teach it. Use the
`STORYBOARD.md` recipes: pick an audience, choose the infographics in order, and narrate each from memory.
If you can run the "Cincinnati / Hagedorn" recipe out loud, you understand the system.

---

## Appendix — the full library

**Infographics** (`docs/infographics/`): 18 numbered Nano Banana Pro prompts, `01`–`18`, in walkthrough
order (each generated image saves alongside as `NN_<name>.png`). See that folder's `README.md` for the
ordered index and `STORYBOARD.md` for the audience deck-recipes.

**Documents** (`docs/`): this guide; `TSC_ENGINE_LEARN_GUIDE` (8–9th grade) and `_ADVANCED` (12th);
`TSC_ENGINE_OMNIVERSE_VISUALS_GUIDE`; `TSC_DEMO_RUNBOOK`; `TSC_INTELLIGENCE_ENGINE_RESEARCH_PAPER` / `_PRD`;
`TSC_DIGITAL_TWIN_RESEARCH_PAPER` / `_PRD`; `TSC_INTELLIGENCE_ENGINE_PEDIATRIC_IMPACT`;
`TSC_INTELLIGENCE_ENGINE_CCHMC_HANDOVER`.

*SYNTHETIC demonstration data throughout · decision support, clinician review required · not FDA-cleared ·
Apache 2.0 · Engine 7 of the HCLS AI Factory.*
