# The TSC AI Engine — A Foundational Learning Guide (Advanced Edition)

### The disease, the clinical problem, and the system — explained in depth

**Author:** Adam M. Jones — HCLS AI Factory, TSC Intelligence Engine (Engine 7)
**Reading level:** written for a strong high-school senior or first-year college reader (roughly 12th grade / AP). It assumes basic biology, a little statistics, and general computer literacy, and it defines specialized terms as they appear. A full glossary closes the guide.
**License:** Apache 2.0 (open source) · 2026

> **How to use this guide.** It builds in layers. Parts 1–2 establish the disease and the clinical
> problem. Part 3 states the design philosophy. Parts 4–10 dissect the system — architecture, the five
> agents, the surfaces, and the 3-D visualization. Parts 11–15 cover how the system is validated, who
> it helps, what it costs, where it goes next, and — importantly — its limits. **Read this first:**
> every result described here was produced on *synthetic* (computer-generated) patient data. No real
> patients are involved yet. The engine is clinician **decision support** — it does not act on its own,
> and it is **not FDA-cleared.** That honesty is not a disclaimer bolted on at the end; as you will see,
> it is built into the engineering.

---

## Part 1 — The biology of TSC: genes, the mTOR pathway, and why growths form

**Tuberous Sclerosis Complex (TSC)** is a rare, inherited, multisystem genetic disorder. To understand
why an AI engine for it looks the way it does, you first need to understand the disease at the level of
molecules.

### The two genes and the brake they build

TSC is caused by loss-of-function changes in one of two genes: **TSC1** (on chromosome 9) or **TSC2**
(on chromosome 16). These genes are the blueprints for two proteins — **hamartin** (from TSC1) and
**tuberin** (from TSC2). Inside the cell, hamartin and tuberin physically bind together to form a
complex. That complex acts as a **brake** on a master growth-control system called **mTOR** —
specifically a hub called **mTORC1** (mechanistic Target Of Rapamycin Complex 1).

Think of mTORC1 as the cell's "grow and divide" accelerator. Under normal conditions, the
hamartin–tuberin complex keeps a foot on the brake, switching mTORC1 off when growth is not needed. When
a TSC1 or TSC2 gene is broken, that brake fails. mTORC1 runs unchecked, and cells grow and proliferate
when they should not. The result is **hamartomas** — benign (non-cancerous) tumors made of disorganized
but normal-type cells — scattered through many organs.

### Why "multisystem" — the organs involved

Because mTOR signaling matters everywhere, TSC produces growths across the body:

- **Brain:** cortical **tubers** (the disease's namesake), **subependymal nodules** lining the fluid-
  filled ventricles, and the **subependymal giant-cell astrocytoma (SEGA)** — a tumor that tends to arise
  near a narrow channel called the *foramen of Monro*. If a SEGA grows there, it can block the normal
  flow of cerebrospinal fluid and cause dangerous pressure (hydrocephalus). Many patients also have
  **epilepsy** (recurrent seizures), often beginning in infancy.
- **Kidneys:** **angiomyolipomas (AMLs)** — fat-, vessel-, and muscle-containing growths. Once an AML
  approaches roughly **4 centimeters**, or develops abnormal vessels (aneurysms), the risk of
  life-threatening bleeding rises sharply.
- **Skin:** hypomelanotic ("ash-leaf") macules, facial angiofibromas, shagreen patches.
- **Heart:** rhabdomyomas (often present at birth and frequently regressing over time).
- **Lungs:** lymphangioleiomyomatosis (LAM), chiefly in women.
- **Eyes:** retinal hamartomas.

### Inheritance, the "two-hit" idea, and mosaicism

TSC is **autosomal dominant**, meaning a single broken copy of TSC1 or TSC2 is enough to cause the
disease. But roughly two-thirds of cases are **de novo** — a brand-new mutation not inherited from
either parent. At the tissue level, individual growths typically follow a **"two-hit" pattern** (the
Knudson hypothesis): the inherited or first mutation disables one gene copy in every cell, and a second,
local mutation knocks out the remaining good copy in the specific cells that then form a tumor.

The most important wrinkle for this engine is **mosaicism**. A **post-zygotic** mutation — one that
happens *after* fertilization, during early development — ends up present in only a *subset* of the
body's cells. A mosaic patient is, genetically, a patchwork. The fraction of cells carrying the variant
is the **variant allele fraction (VAF)**. In blood, that fraction may be very low (single-digit
percentages), and conventional sequencing simply does not see it. The consequence is stark: by standard
blood-based testing, roughly **10–15% of clinically diagnosed TSC patients are reported as "no mutation
identified" (NMI)** — not because they lack a causal variant, but because the test cannot detect a
low-level mosaic change. Those variants are often *enriched in affected tissue* (for example, a brain
tuber removed during epilepsy surgery), where the VAF is much higher. **This single biological fact —
mosaic variants hidden from blood but recoverable from tissue — is the clinical opening the engine's
flagship capability is built around.**

---

## Part 2 — The clinical reality: why managing TSC is hard

Biology sets the stage; the day-to-day difficulty of TSC care is where an assistive system earns its
keep. Six problems recur.

### 2.1 The diagnostic odyssey

As above, NMI patients may go years without a confirmed molecular diagnosis. A confirmed variant matters:
it can clarify prognosis, enable family testing and counseling, qualify patients for targeted therapy or
trials, and end the uncertainty that itself burdens families. The path to that answer often requires
recognizing that **tissue, not blood, is where the signal lives** — and applying mosaic-aware analysis to
it. That recognition is inconsistent in practice.

### 2.2 Fragmented, multidisciplinary care

A single patient may be followed by neurology, nephrology, genetics, dermatology, neurosurgery,
psychiatry, and more. Each specialist sees a slice. The longitudinal, whole-patient picture — the thing
needed to make good timing decisions — is assembled, if at all, by effortful human coordination.

### 2.3 A heavy surveillance burden

The **2021 International TSC Consensus** recommendations specify a schedule of brain MRIs, renal imaging,
EEGs, neuropsychiatric screening, skin checks, and more, at organ-specific intervals. Keeping every
patient on schedule — and intensifying surveillance for those at higher risk — is real cognitive and
logistical work, and gaps happen.

### 2.4 High-stakes timing decisions

The hardest TSC decisions are fundamentally about *when*, and *with what confidence*:

- **SEGA:** as it grows toward the size where it threatens CSF outflow, the team must decide between
  watchful waiting, starting an **mTOR inhibitor** (for example, **everolimus**, which directly targets
  the overactive pathway), or surgical resection.
- **AML:** as it approaches the ~4 cm bleeding-risk threshold, the choices are surveillance, an mTOR
  inhibitor, or pre-emptive embolization.
- **Epilepsy:** balancing seizure control against medication burden and developmental impact.

Each is a function of a *trajectory* (how fast is it growing?) and an *uncertainty* (how sure are we?).
Clinicians estimate both from a handful of prior measurements, mentally — skilled but invisible,
unshareable, and uncalibrated work.

### 2.5 Under-recognized neuropsychiatric disease (TAND)

**TSC-Associated Neuropsychiatric Disorders (TAND)** — spanning behavior, psychiatric, intellectual,
academic, neuropsychological, and psychosocial domains — affect roughly **90%** of people with TSC, yet
fewer than **20%** are formally assessed. A validated instrument exists (the **TAND Checklist**), but the
signals frequently surface first as soft, hedged language buried in clinic notes ("mom reports some
focus issues at school; will revisit") that never coalesce into a documented pattern.

### 2.6 Information overload, alert fatigue, and the trust deficit

Clinical data is fragmented across electronic health record (EHR) systems, imaging archives, and free
text. Software that responds with a flood of alerts trains clinicians to ignore alerts entirely (**alert
fatigue**). And clinicians rightly distrust any tool that cannot show *why* it said what it said. Any
credible system must be **disciplined** about when it speaks and **transparent** about where its claims
come from.

---

## Part 3 — Design philosophy: what kind of system this is (and is not)

Before the architecture, the principles — because they explain every later choice.

**It is decision support with a human in the loop, never an autonomous actor.** The engine drafts,
forecasts, and surfaces; a clinician decides. The molecular-genetics output, for instance, is an explicit
*draft* that requires a board-certified geneticist's sign-off.

**Determinism where correctness is non-negotiable; AI only at the edges.** The genetic verdict, the
orchestration logic, and the numerical forecasting are handled by *deterministic* code — fixed rules and
classical statistics that produce the same answer every time and can be audited. Large language models
(LLMs) are used only for what they are genuinely good at: extracting structure from free text and writing
explanatory narrative. Crucially, an LLM never *overrules* a deterministic verdict; it explains it.

**Honesty is an engineering property, not a disclaimer.** Three concrete commitments:
- **Provenance everywhere** — every output records which model produced it, at what tier, how long it
  took, and which prompt version was used.
- **Uncertainty made literal** — forecasts always carry calibrated intervals, and (as you will see) the
  3-D visualization renders that uncertainty as geometry that *cannot* imply more confidence than the
  math supports.
- **Construct-validity labeling** — results on synthetic data are labeled as exactly that, never dressed
  up as clinical proof.

**Small, open, and reusable.** The system runs on a single, affordable computer; its code is openly
licensed; and its structure is deliberately disease-agnostic so it can be retargeted to other rare
diseases by swapping disease-specific components.

---

## Part 4 — System architecture

At the highest level: **five agents**, a **deterministic orchestrator**, and a set of **clinician
surfaces**, all reading from a shared, event-sourced record.

### 4.1 Agents, orchestrator, surfaces

An **agent** is a focused software component that performs one analytical job and emits a structured
result. The five are the Phenome Mapper, Variant Curator, Trajectory Modeler, TAND Surveillance Agent,
and Therapeutics Strategist (detailed in Parts 5–9). The **orchestrator** coordinates them; the
**surfaces** present the assembled results.

### 4.2 Event sourcing: an append-only, replayable record

Rather than overwriting a database as it goes, the engine uses **event sourcing**. Every meaningful
occurrence — a patient enrolled, an agent finished, a surface requested — is appended as an immutable
**event** to a log that is never edited. There are **13 event types**. The current state a clinician sees
(the "projection") is *computed* by replaying that log. This design has powerful properties: nothing is
silently changed, the entire history can be reconstructed and audited, and the system is reproducible —
the same inputs always rebuild the same state.

### 4.3 Deterministic, dependency-ordered orchestration

The orchestrator is not an AI; it is a rule-following coordinator (think air-traffic control). It runs
agents in **dependency order** — the Phenome Mapper first, because the others build on the phenotype it
produces; then the Variant Curator, Trajectory Modeler, and TAND agent; then the Therapeutics Strategist,
which integrates everything. If an agent fails, the orchestrator records the failure and marks the
affected section **"pending"** rather than displaying stale or fabricated data — a conservative,
fail-safe stance. It also supports **incremental updates**: when a single new clinical note arrives, only
the note-dependent agents (Phenome Mapper and TAND) re-run, not the whole pipeline.

### 4.4 Model tiering and cost control

When LLMs are used, they are routed by **stakes**. A configuration policy maps each step to a tier:
**Haiku** (fast, inexpensive) for light extraction, **Sonnet** (balanced) for analysis, and **Opus**
(most capable) for high-stakes narrative such as the genetic interpretation — with a **local Llama**
model as an offline fallback. A **cost ledger** tracks token usage and dollar cost per tier and can
enforce a spending cap, so the system's economics are observable and bounded.

### 4.5 "Runs-here" fallbacks

A production deployment might use PostgreSQL (database), Redis (fast cache), MinIO (file storage), and
Milvus (vector search). The engine implements all of those, but each has a lightweight **fallback** —
SQLite, in-memory structures, the local filesystem, and an in-memory vector store — selected
automatically when the production service is absent. The effect: the entire system **runs on a bare
single computer, offline, with no external services**, which is exactly what makes a cheap-hardware demo
possible.

---

## Part 5 — Agent 1: the TSC-Phenome Mapper

**Goal:** construct an ontology-grounded, longitudinal map of the patient's **phenotype** (observable
traits and findings) over time — the foundation every other agent reads.

**Standardizing with the HPO.** Clinicians worldwide encode phenotypes using the **Human Phenotype
Ontology (HPO)** — a controlled vocabulary in which each concept (e.g., "Renal angiomyolipoma") has a
single stable identifier. The Phenome Mapper extracts findings from structured records and free-text
notes and assigns each one its HPO code.

**Why grounding matters — and how it is enforced.** An LLM asked for an HPO code can hallucinate a
plausible-looking but wrong identifier. The Phenome Mapper therefore **validates every term against the
real, official HPO release (19,389 terms)**: it corrects labels to the canonical term, remaps obsolete or
alternate IDs to their current form, *recovers* the correct ID when only the text name is right, and
**drops codes that are not real HPO terms entirely** rather than trusting them. (In testing, this check
even caught two miscoded entries in the synthetic data itself — a SEGA mistakenly coded as "optic nerve
glioma," and cortical tubers coded as a different white-matter finding.)

**Span-grounded extraction with polarity and temporality.** When reading notes, the agent does not just
say "this phenotype is present." It records the **exact character span** of the supporting text (and
verifies that the quoted span really matches the note), the **polarity** (is the finding asserted as
*present* or explicitly *negated*?), and the **temporality** (is it an *onset* finding, *current*, or
*historical*?). Negated findings are logged but never admitted as active phenotypes. When a finding is
present in one part of the record and negated in another — say, infantile spasms active in infancy but
resolved now — the agent records a **discordance**, surfacing it for clinician reconciliation rather than
silently picking one.

**Surveillance-gap analysis.** Finally, a deterministic rule engine compares the patient's findings and
timeline against a simplified ITSC surveillance schedule and flags items that are overdue.

---

## Part 6 — Agent 2: the TSC-Variant Curator

**Goal:** identify the causal genetic variant, classify its significance, and — the flagship capability —
**recover low-VAF mosaic variants that standard blood testing misses.**

### 6.1 From sequencing reads to a confident call

Genetic data arrives as a list of differences from a reference genome, each annotated with read-level
evidence: total **read depth**, the count of reads supporting the variant (**alt reads**), the **variant
allele fraction (VAF)**, and the **strand balance** (whether the supporting reads come evenly from both
DNA strands). These details are not decoration — they are how real from artifact is told apart.

To make that distinction concrete, the engine's synthetic data deliberately seeds, in every sample, a
**low-VAF, strand-biased artifact** — a look-alike "junk" call whose supporting reads pile up almost
entirely on one strand, a classic technical artifact. A naive tool would report it; the Variant Curator
must **reject** it. Because the curator demonstrably filters such artifacts, its claim of **zero
false-positive pathogenic calls** is a real discrimination, not a freebie.

### 6.2 Annotation and ACMG-AMP classification

For each genuine variant, the agent annotates the consequence (e.g., nonsense, frameshift, missense),
whether the gene acts through loss of function, population frequency (absence from large databases like
gnomAD is itself evidence), and any prior clinical classification. It then applies the international
standard for variant interpretation, **ACMG-AMP** (Richards et al., 2015). ACMG works by assembling
weighted **criteria** and combining them by fixed rules into a verdict:

- **PVS1** ("pathogenic, very strong") — a null variant in a gene where loss of function causes disease.
- **PM2** ("pathogenic, moderate") — absent or vanishingly rare in population databases.
- **PP4** ("pathogenic, supporting") — the patient's phenotype is highly specific for the disease.
- **PS1** ("pathogenic, strong") — the exact change is already established as pathogenic.

These combine into **Pathogenic**, **Likely Pathogenic**, or **Variant of Uncertain Significance (VUS)**,
among others. A crucial subtlety the engine models: for a **low-VAF mosaic** null variant, the strength of
PVS1 is conservatively **down-weighted** (from "very strong" toward "strong"), reflecting genuine
guideline practice for less-certain contexts. The featured demo patient — a *TSC2* frameshift recovered at
about 8.3% VAF in tissue — is therefore classified **Likely Pathogenic** rather than Pathogenic, exactly
as a careful geneticist would.

### 6.3 The deterministic verdict, the AI narrative, and the human gate

The classification itself is computed by a **deterministic, by-the-book combinatorial classifier** — the
authority. An LLM (the high-stakes Opus tier) is asked only to write the per-criterion *explanation* in
clinician-readable prose; it cannot change the verdict. The agent flags the variant as **mosaic** and
**recovered** when appropriate, recommends **orthogonal validation** (an independent assay such as ddPCR),
and produces an **AI-labeled DRAFT molecular-genetics report**. That draft is explicitly *pending a
board-certified molecular geneticist's sign-off* — a non-negotiable **human gate**, exposed as a formal
approve/reject/amend step.

---

## Part 7 — Agent 3: the TSC-Trajectory Modeler

**Goal:** forecast how the patient's monitored quantities will evolve over 6, 12, and 18 months, *with
honestly calibrated uncertainty.* This agent is deliberately **classical statistics, not an LLM** —
because forecasting numbers with calibrated intervals is precisely what statistical models do well and
chatbots do poorly.

### 7.1 Borrowing strength: the population prior and Bayesian shrinkage

A single patient usually has only a few past measurements — too few for a confident trend on their own.
The engine addresses this with a **mixed-effects population model**: it fits growth patterns across a
population of similar trajectories to derive a **prior** expectation for how, say, SEGAs typically grow.
The individual patient's simple trend (an ordinary least-squares slope) is then **shrunk toward that
prior**, weighted by how trustworthy each is — a **Bayesian** update. Sparse, noisy individual data leans
more on the population; clean individual data barely moves. For quantities without a population model, the
engine performs a conjugate-normal Bayesian update against a weak prior, so the method is principled in
every case.

### 7.2 Calibrated intervals and the "cone of uncertainty"

To express how *confident* the forecast is, a **Gaussian-process** model estimates the width of the
prediction at each horizon, producing **50% and 90% prediction intervals**. Conceptually this is the
**"cone of uncertainty"** from a hurricane forecast: the most likely path plus a widening band that
honestly reflects diminishing certainty further out. Tight data → narrow cone; sparse data → wide cone.

### 7.3 Survival-style crossing probability and graded risk

Beyond "how big," clinicians need "how likely, and when, will it cross a threshold." From the forecast's
distribution at each horizon, the agent computes a **probability of crossing** the clinical threshold (for
a SEGA, the size at which intervention is discussed; for an AML, ~4 cm) — a survival-analysis-style
treatment of time-to-event. It then **grades** the risk as **likely / possible / unlikely** and estimates
a **time-to-threshold.**

### 7.4 Multiple quantities and a cadence recommendation

The agent forecasts **four quantities at once**: SEGA size, AML size, **renal function (eGFR)**, and
**seizure frequency** — note that eGFR *declines* toward a lower threshold, the opposite direction from a
growing lesion, which the model handles explicitly. Finally, it recommends a **surveillance cadence**:
when crossing risk is elevated, it suggests checking the patient *more often* than the standard ITSC
floor; when stable, it holds the floor. An LLM is used only to write a one-sentence plain-language summary
— and is constrained to introduce no number the model did not produce.

---

## Part 8 — Agent 4: the TAND Surveillance Agent

**Goal:** surface the under-recognized neuropsychiatric (TAND) signals hidden in longitudinal notes —
gently, and never as an alarm.

**Discourse-marker detection.** Concerns about behavior, mood, attention, and learning often appear first
as *hedged, tentative* language. The agent uses a **diagnostic-uncertainty discourse-marker taxonomy** —
a method developed by researchers at Cincinnati Children's (the **Marshall-Hagedorn** approach) — to
detect exactly those markers (hedging, deferral, third-party attribution, conditional language) across the
six recognized **TAND clusters**: behavioral, psychiatric, intellectual, academic, neuropsychological, and
psychosocial.

**Weighted, longitudinal, negation-aware scoring.** A transparent, configurable scoring layer (not an LLM)
aggregates the signal: each marker type carries a **weight**, more **recent** mentions count for more (a
half-life decay), and **recurrence** across multiple notes adds confidence. Critically, the scoring is
**negation-aware**: a note stating "mom *denies* behavior concerns" contributes *nothing* — preventing the
classic false-positive where a denied concern is mistakenly counted as a finding. A cluster is flagged
only when its weighted longitudinal score crosses a configured threshold.

**Briefing, never an alert.** TAND findings are surfaced strictly as **pre-visit briefing material** for
the clinician to consider — never as an interruptive alert and never as a diagnosis. This is a deliberate
guard against alert fatigue and against over-medicalizing soft signals.

---

## Part 9 — Agent 5: the TSC-Therapeutics Strategist

**Goal:** integrate everything into a structured, fully sourced summary of therapeutic options.

The Therapeutics Strategist reads the variant interpretation, the phenotype, the forecasts, and the TAND
findings, and adds external context: matching **clinical trials**, relevant published evidence retrieved
via **RAG** (retrieval-augmented generation — a technique that fetches relevant source documents so the
model's statements are grounded in real text rather than invented), and applicable FDA actions. It
produces a **six-section options brief**. Its governing rule is **source attribution**: every claim names
where it came from. And its framing is strictly **decision support** — it lays out options and trade-offs;
it does not issue directives.

---

## Part 10 — The surfaces and the Omniverse digital twin

### 10.1 The clinician surfaces

The assembled results are presented through demand-driven **surfaces**, each computed from the projection:

- **Pre-Visit Briefing** — a single-screen, mobile-readable summary with at most three action items.
- **In-Visit Dashboard** — the working screen during a visit: four quadrants for variant, phenotype,
  trajectory, and TAND, with every claim traceable to its source and a browsable audit trail.
- **Async Alerts** — a *disciplined* alert surface governed by a per-clinician **sliding-window** counter
  that triggers recalibration if alerts would exceed roughly three per clinician per week; every alert
  carries a reason and a dismissal path.

### 10.2 The anatomical digital twin (OpenUSD + Omniverse)

The most visually arresting layer renders the engine's *actual outputs* as navigable 3-D scenes using
**OpenUSD** (a 3-D scene-description standard) and **NVIDIA Omniverse** (a high-fidelity rendering
platform). Four scene types exist:

- **The lesion-trajectory twin.** The lesion is a 3-D form whose radius is **animated along the forecast**
  as a timeline is scrubbed from the past into the future, wrapped in a translucent **uncertainty
  envelope.** Here is the elegant, honest design choice: **the radii of that envelope are set exactly
  equal to the engine's 50% and 90% prediction intervals.** The cloud is, geometrically, the forecast's
  uncertainty — so a confident forecast renders as a tight shell and an uncertain one as a large, soft
  cloud. *The picture cannot overstate the engine's confidence.* A red "membrane" marks the clinical
  threshold, and the lesion changes color as it approaches.
- **The mosaic "powers-of-ten" scene.** A field of cells in which **exactly the recovered VAF fraction**
  carries the variant — about **1 cell in 12** for an 8.3% VAF — rendered as glowing cells via efficient
  point-instancing. It makes "8.3% VAF" something you can literally count, dramatizing the signal a blood
  test reports as negative.
- **The whole-child organ atlas.** A stylized body with each affected organ system illuminated from the
  patient's phenome — the multisystem disease at a glance.
- **The population view.** All 50 synthetic patients as a grid of figures, colored by classification, with
  the **seven recovered mosaic patients ringed in gold halos.**

For film-quality rendering, the scenes use **MDL materials** (NVIDIA's material system) so envelopes
render as true glass and variant cells and halos glow — while a plain-color fallback keeps the scenes
viewable in any USD viewer. Notably, the lesion's *color* stays driven by the engine's animated values
rather than a fixed material, so its threshold-crossing color change survives.

**The authoring/render split — why it stays cheap.** Building these scenes is a deterministic, CPU-only
transformation that runs on the small computer with **no GPU and no special dependency**. The expensive,
RTX-accelerated *rendering* runs separately on a rented cloud GPU (via **RunPod**). In short: **the small
computer thinks and authors; the cloud GPU renders.** One honest caveat: the anatomy is *stylized,
atlas-based geometry scaled by real measurements* — not a reconstruction of a specific patient's actual
scans. Using segmented patient imaging is a future, institution-stage upgrade.

---

## Part 11 — Validation: measuring rather than asserting

A demonstration that merely *looks* convincing proves little. The harder discipline — and the one that
earns clinical credibility — is to **measure** whether the engine is *correct*. Because the demo cohort is
synthetic, the team knows the ground-truth answers in advance and can grade the engine against them, like
a practice exam with an answer key.

### 11.1 The scorecard

| Metric (on the synthetic cohort, with known ground truth) | Result |
| --- | --- |
| Variant-classification accuracy | 100% (50/50) |
| Diagnostic **detection** uplift (standard pipeline → engine) | 86% → 98% (**+12 points**) |
| Sub-threshold mosaic cases recovered | **6 / 6** (100% sensitivity) |
| Truncating (null) variants wrongly called benign | **0** |
| Forecast 90% interval coverage (held-out probe) | 100% |
| Outputs carrying full provenance | 100% |

### 11.2 Detection yield vs. molecular-diagnosis yield — an honest distinction

The engine carefully separates two different "yield" ideas. **Detection** asks: did the pipeline find the
variant *at all*? A standard, blood-focused pipeline reliably calls only higher-VAF variants, so it misses
the low-level mosaics — and the engine's mosaic-aware path recovers them (the +12-point detection uplift,
all true positives, with zero phantom variants invented in the genuinely variant-free patient).
**Molecular diagnosis** is the stricter bar of detection *plus* a returnable Pathogenic/Likely-Pathogenic
call. Reporting both — rather than conflating them — is part of the honesty discipline.

### 11.3 The most important caveat

> These results establish **construct validity** — proof that the engine's logic faithfully recovers the
> signal that was deliberately planted in synthetic data with known answers. **This is not the same as
> clinical validity** — proof that the engine performs on *real patients, prospectively, against orthogonal
> assays and expert adjudication.* That larger study has not been done.

Why emphasize a limitation? Because the honesty is precisely what makes the tool trustworthy — and it is
ethically required. An apt analogy is a **flight simulator**: it can rigorously demonstrate that a pilot's
skills and a plane's systems function, but it is not the first real flight. This engine is a very good
simulator; the real flight — running the Variant Curator on actual specimens — is the defining next step.

---

## Part 12 — Benefits mapped to problems, and who is helped

| Clinical problem (Part 2) | How the engine addresses it |
| --- | --- |
| Diagnostic odyssey / mosaic blind spot | Variant Curator recovers low-VAF mosaics from tissue (demo: 6/6; +12-point detection uplift), with artifact rejection making the result rigorous |
| Fragmented multidisciplinary care | Orchestrator integrates all agents into one shared, longitudinal projection and a unified dashboard |
| Hard timing decisions | Trajectory Modeler forecasts four quantities with calibrated intervals, crossing probabilities, graded risk, and cadence recommendations |
| Under-recognized TAND | TAND agent surfaces hedged, recurrent signals as gentle briefing material, negation-aware |
| Alert fatigue | Disciplined alert surface with a per-clinician sliding-window budget |
| Scattered, untrusted data | Provenance on every output; honest, geometric uncertainty; explicit construct-validity labeling |

**Stakeholders:** *patients and families* gain a shorter odyssey, earlier warning before emergencies, and
sooner recognition of learning and behavioral needs; *clinicians* gain the whole picture in one place with
the math done and sources cited, freeing attention for judgment; *care teams* gain a consistent shared
view across specialties; and the *research community* gains an open, inspectable platform to build on.

---

## Part 13 — Hardware, cost, openness, and generalization

**The hardware.** The full engine runs on a single **NVIDIA DGX Spark** — a compact desktop AI computer
(a GB10-class chip with a large pool of unified memory shared between CPU and GPU) priced around
**$4,699.** The deliberate, GPU-heavy steps (high-fidelity rendering; and, in future, accelerated genomic
analysis) burst to rented cloud GPUs rather than requiring a data center. The thesis is that a clinically
interesting, multi-agent precision-medicine system can run first on hardware a single lab or hospital
group can actually afford.

**Openness.** The code is **Apache 2.0**, meaning it can be inspected, trusted, reused, and extended
openly — the opposite of a closed commercial black box, and a posture that matters for clinical trust.

**Generalization.** The architecture is intentionally a **reusable pattern**, not a one-off. A different
disease supplies a different anatomical atlas and a different vocabulary of variants, lesions, and
phenotypes; the orchestrator, the event-sourced core, the uncertainty machinery, the provenance, and the
visualization are reused unchanged. The same structure could be retargeted to **NF1, NF2, Rett syndrome,
Williams syndrome, and other mTOR-pathway disorders**. That replicability — "swap the box labels, keep the
wiring" — is why it is called an *engine* rather than a project.

---

## Part 14 — The path to real impact: institutions, the biobank, and Phase 2

The engine today is a complete demonstration on synthetic data. Converting that into medical impact
follows a deliberate two-phase plan.

**Phase 1 — the demonstration (built):** the full five-agent engine on the DGX Spark with the synthetic
cohort, plus the Omniverse 3-D twins rendered via RunPod. This is the version that earns institutional
engagement.

**Phase 2 — the institution.** A real deployment connects to a hospital's components. For a beachhead at
**Cincinnati Children's**, five touchpoints map cleanly onto the engine: the **Discover Together Biobank**
(banked surgical tissue — the substrate that feeds the Variant Curator and the highest-leverage piece);
the research infrastructure that houses it; the **Division of Biomedical Informatics** (the methodological
and sponsorship home, and the source of the TAND discourse methodology); the **TSC clinical program** (the
patients); and the EHR plus laboratory information system (the data plumbing).

**The first real study — low-risk by design.** The natural opening is a **retrospective "no-mutation-
identified" re-analysis**: take patients historically reported NMI on blood, and run the mosaic-aware
Variant Curator on their *already-banked, already-consented* surgical tissue. If it recovers diagnoses
that were missed, that is a concrete, publishable result — obtained without new surgery, new collection,
or a long prospective wait. It is the fastest credible route from *construct validity* to a *real clinical
result*, and it converts every other touchpoint (informatics governance, the LIMS that links specimens to
patients, the clinic that consented the patients) from optional into necessary.

---

## Part 15 — Limitations, risks, and responsible use

A serious guide names its own edges plainly.

- **Synthetic data only, so far.** Every metric is construct validity. The engine has not yet diagnosed a
  real patient faster, managed a real lesion earlier, or changed a real decision. The biobank study is the
  inflection point.
- **Not FDA-cleared; not autonomous.** It is investigational decision support behind a clinician-review
  gate, with a non-negotiable molecular-geneticist sign-off on genetic findings.
- **Documented technical gaps.** The faithful, accelerated genomic-calling substrate (spiking variants
  into real read data and calling them blind) requires a burst GPU and is not built; live EHR, biobank
  LIMS, and imaging-AI/PACS integrations are institution-stage work; the clinical notes and imaging are
  synthetic; and the 3-D anatomy is stylized rather than segmented from a patient's real scans.
- **Equity and bias.** Any model trained or tuned on non-representative data risks unequal performance.
  Rare-disease populations are small and can be skewed; prospective validation must include attention to
  whether performance holds across ancestries, ages, and presentations — a reason the synthetic-data
  honesty matters even more, not less.
- **Human judgment is the point, not a backstop.** The system is engineered to *augment* clinicians'
  attention and consistency, not to replace their reasoning. Its disciplines — provenance, honest
  uncertainty, conservative failure, the human sign-off — exist to keep a person meaningfully in control.

Used within those boundaries, the engine is a credible, inspectable, affordable instrument for testing a
genuine clinical hypothesis: that mosaic-aware analysis and calibrated forecasting can measurably improve
TSC care. Used outside them — as if it were a finished, validated clinical product — it would overclaim.
The whole design is built to make the honest use the natural one.

---

## Glossary

- **ACMG-AMP** — the international standard (Richards et al., 2015) for classifying genetic-variant
  pathogenicity by combining weighted criteria.
- **Agent** — a focused software component performing one analytical job; the engine has five.
- **Allele / variant** — a specific version/change of a gene at a position in the genome.
- **AML (angiomyolipoma)** — a kidney growth in TSC; risk of hemorrhage rises near ~4 cm.
- **Bayesian update / shrinkage** — combining prior expectation with observed data, weighted by the
  reliability of each, to produce a posterior estimate.
- **Construct validity** — evidence that a system's logic works on data with known answers (achieved here).
- **Clinical validity** — evidence that it works prospectively on real patients (not yet done).
- **Deterministic** — producing the same output from the same input by fixed rules; not an AI guessing.
- **eGFR** — estimated glomerular filtration rate; a measure of kidney function.
- **Event sourcing** — recording state as an append-only log of immutable events that can be replayed.
- **Everolimus** — an mTOR-inhibitor drug used in TSC (e.g., for SEGA and AML).
- **Gaussian process** — a statistical model used here to estimate the width of forecast uncertainty.
- **Hamartin / tuberin** — the proteins made by TSC1 / TSC2 that together restrain mTORC1.
- **Hamartoma** — a benign, disorganized overgrowth of normal-type tissue.
- **HPO (Human Phenotype Ontology)** — the standardized vocabulary of phenotype codes.
- **LLM (large language model)** — an AI model for understanding and generating text.
- **Mixed-effects model** — a statistical model capturing both population-level and individual-level
  effects; used to build the growth prior.
- **Mosaicism** — presence of a genetic variant in only a subset of a person's cells (post-zygotic origin).
- **mTOR / mTORC1** — a master growth-signaling pathway/hub that the TSC proteins normally restrain.
- **NMI ("no mutation identified")** — a genetic test reporting no causal variant, often due to mosaicism.
- **OpenUSD** — a 3-D scene-description standard used for the digital twins.
- **Orchestrator** — the deterministic coordinator that runs agents in order and logs everything.
- **Phenotype** — the set of a person's observable traits and clinical findings.
- **Prediction interval** — a range expressing forecast uncertainty (here, 50% and 90%).
- **Provenance** — the recorded origin of an output (model, tier, latency, prompt version).
- **PVS1 / PM2 / PP4 / PS1** — specific ACMG evidence criteria of varying strength.
- **RAG (retrieval-augmented generation)** — fetching relevant documents so an AI can cite real sources.
- **SEGA** — a TSC brain tumor near the foramen of Monro that can obstruct CSF flow.
- **Strand balance** — whether variant-supporting reads come evenly from both DNA strands; imbalance
  suggests artifact.
- **Surveillance** — the scheduled monitoring (per the 2021 ITSC consensus) that catches TSC problems early.
- **Synthetic data** — realistic but computer-generated data with no real patients.
- **TAND** — TSC-Associated Neuropsychiatric Disorders (behavior, mood, cognition, learning, social).
- **TSC (Tuberous Sclerosis Complex)** — the multisystem genetic disorder this engine targets.
- **VAF (variant allele fraction)** — the fraction of sequencing reads (≈ cells) carrying a variant.

---

*SYNTHETIC demonstration data throughout. The TSC Intelligence Engine is investigational clinician
decision-support, requires clinician review, and is not FDA-cleared. Apache 2.0. Engine 7 of the HCLS AI
Factory.*
