# The TSC Anatomical Digital Twin

### Volumetric, Time-Resolved, Uncertainty-Aware Visualization of Multi-Agent Precision-Medicine Output — an Omniverse / OpenUSD Layer for Engine 7

**Author:** Adam M. Jones — HCLS AI Factory, TSC Intelligence Engine (Engine 7)
**Companion to:** *The TSC Intelligence Engine* (master research volume) and the *TSC Intelligence Engine PRD (v0.1)*
**Paired build document:** *TSC Anatomical Digital Twin — Visualization Layer PRD (v0.1)*
**License:** Apache 2.0 · **Status:** Design research preceding the build · 2026

> **What this document is.** A research-grade design rationale and technical foundation for a fourth clinician surface on the TSC Intelligence Engine: a volumetric, time-resolved *anatomical digital twin* that renders the engine's existing structured output — recovered mosaic variants, lesion-growth forecasts with prediction intervals, threshold-crossing windows, and HPO-coded multisystem phenotype — as a navigable 3D scene in NVIDIA Omniverse via OpenUSD. It is the visualization counterpart to the engine's analytic core. It introduces no new clinical claims; it is a new *way of seeing* claims the engine already makes and already validates.

---

## Implementation Status — Specification, Honestly Labeled (June 2026)

This is a **forward-looking design specification**, not a record of a completed build. It is written to the same standard as the engine's master volume so that the build, when it happens, is the argument — but the reader should hold the right expectation:

- **The engine it builds on is real.** The TSC Intelligence Engine v0.1 is built and verified: five agents, a deterministic orchestrator, three clinician surfaces, an event-sourced state store, an HPO-ontology-grounded phenome layer, and a validation harness that measures — against synthetic ground truth — variant-classification accuracy (1.00 construct validity), a +12-point diagnostic *detection* uplift driven by 100% recovery of sub-threshold mosaics, and 100% output provenance completeness. Those numbers are what this layer renders.
- **This layer is specified here and not yet implemented.** The companion PRD defines the build. Its first phase (Phase V0 — the CPU-side OpenUSD exporter and the lesion-trajectory scene) is buildable on the existing DGX Spark **today, with no new GPU and no new clinical data**, because authoring USD is a deterministic data transformation, not a render.
- **The honest hardware line is drawn early and kept.** The Spark authors USD and runs the engine; high-fidelity RTX rendering is a companion step on an RTX workstation or a burst RunPod instance. Native Omniverse Kit support on the Grace-Blackwell (GB10) platform should be verified rather than assumed; nothing in the core design depends on it.
- **The anatomy is stylized, and the document says so on every page that matters.** The synthetic cohort has no patient imaging. The twin is an *atlas-anchored* stylization scaled by the engine's real measurements — not a patient-specific reconstruction. The path to patient-specific anatomy (DICOM → MONAI segmentation → USD) is institutional Phase-1 and is specified as such.

---

## Table of Contents

1. Abstract and thesis
2. Part I — Motivation: the cognitive gap between tabular output and spatial-temporal disease
3. Part II — Scientific foundations: what a digital twin is here, and the data that already supports it
4. Part III — Technology and architecture: OpenUSD, Omniverse, and the authoring/render split
5. Part IV — The three signature scenes
6. Part V — Clinical integration, safety, and validation
7. Part VI — Generalization, roadmap, and references

---

## 1. Abstract and thesis

The TSC Intelligence Engine converts a child's genomic and longitudinal clinical record into decision-support output: a recovered mosaic variant with an ACMG-AMP classification and a provenance trail; a SEGA or angiomyolipoma (AML) growth forecast at 6, 12, and 18 months with 50% and 90% prediction intervals and a time-to-threshold estimate; an HPO-coded longitudinal phenotype with surveillance gaps; an under-recognized TAND signal; a source-attributed therapeutic options brief. Today that output is presented as cards, tables, and a two-dimensional growth chart. This is correct, auditable, and — for the specific decisions a TSC clinician makes — *under-dimensioned*.

The decisions are spatial and temporal. *When* will this SEGA reach the size at which the multidisciplinary team should discuss resection, and with *what* confidence? *Where*, relative to the foramen of Monro, is it, and is it the lesion that threatens CSF outflow? Is this 3.6 cm AML the one approaching the ~4 cm bleeding-risk threshold? A clinician answers these by building a mental model — a spatial, growing, uncertain object in time — from numbers on a page. The thesis of this document is simple: **the engine already computes that object; we should render it.**

We propose the **TSC Anatomical Digital Twin**: a fourth clinician surface that takes the engine's existing projections and authors them, deterministically, as an OpenUSD scene — anatomically anchored lesion geometry scaled to the engine's real measurements, animated across the observed history and into the forecast, wrapped in a translucent *uncertainty envelope* whose radii are exactly the engine's 50%/90% prediction intervals, and colored by proximity to the clinical threshold. The same provenance that travels with every analytic output travels into the scene as queryable metadata. Rendered in NVIDIA Omniverse with RTX path tracing, it is film-quality; authored in pure Python on the Spark, it is reproducible and free. It makes the engine's single most differentiating capability — forecasting the future of a lesion *with calibrated uncertainty* — something a clinician can see, scrub, and trust.

Crucially, the twin makes **no new claim**. Its uncertainty envelope is the antidote to the usual danger of medical 3D — that a beautiful render implies a false precision. Here the render's translucency *is* the prediction interval; the wider the engine's uncertainty, the larger and softer the cloud. Honesty is built into the geometry.

---

# Part I — Motivation

## 2. The cognitive gap

### 2.1 Disease is spatial and temporal; our output is neither

Tuberous Sclerosis Complex is a multisystem disease whose clinical management is dominated by *lesions that change over time in specific anatomical locations*: subependymal nodules and the subependymal giant-cell astrocytomas (SEGA) that can arise from them near the foramen of Monro; cortical tubers; renal angiomyolipomas that grow toward a bleeding-risk threshold; cardiac rhabdomyomas that often regress; pulmonary lymphangioleiomyomatosis; dermatologic features. The clinician's core surveillance task — codified in the 2021 International TSC Consensus recommendations — is fundamentally about *tracking objects in space over time and acting before a threshold is crossed*.

The engine's analytic output is faithful to this task in its *content* but not in its *form*. A growth forecast rendered as a table of numbers, or even as an excellent two-dimensional line chart with shaded prediction bands, asks the clinician to perform a mental transformation: to imagine the lesion, place it anatomically, inflate it along the trajectory, and hold the uncertainty in mind. Skilled clinicians do this constantly and well. But the transformation is cognitive work, it is error-prone under time pressure, and — critically — it is *not shareable*. A parent in the room, a trainee, a colleague in a different specialty cannot see the clinician's mental model. A rendered twin externalizes it.

### 2.2 Three decisions that volume, time, and uncertainty change

The case for the twin is not aesthetic; it is decision-specific. Three recurring TSC decisions are each a function of exactly the three dimensions the twin adds — spatial extent, temporal trajectory, and calibrated uncertainty:

| Decision | Spatial | Temporal | Uncertainty | Why 2D under-serves it |
| --- | --- | --- | --- | --- |
| SEGA resection / mTOR-inhibitor timing | Proximity to foramen of Monro and CSF pathways governs urgency | Growth rate and time-to-threshold drive *when* | The cost of acting early vs. late depends on the *width* of the forecast | A line chart shows size-vs-time but not anatomical risk, and shows bands but not their volumetric meaning |
| AML embolization / mTOR initiation | Lesion size and number; ~4 cm and vascular features raise bleeding risk | Slow growth across years; the question is the next window | Whether the lesion is *at* or merely *near* threshold | A number ("3.6 cm") flattens the difference between "stable below" and "accelerating toward" |
| Surveillance intensity and modality | Which organs are involved (from the HPO profile) sets the scan plan | Interval since last imaging vs. recommended cadence | Confidence that a quiet organ is truly quiet | A gap list is correct but does not convey *whole-child* burden at a glance |

In each row, the rightmost column is the gap the twin closes. The engine has already done the hard quantitative work; the twin is the display that matches the work to the decision.

### 2.3 Externalizing the model: trainees, families, multidisciplinary teams

TSC care is multidisciplinary by definition — neurology, nephrology, genetics, dermatology, neurosurgery, psychiatry, and the family — and the coordinating clinician spends real effort communicating a shared picture. A volumetric twin is a *communication artifact* as much as a clinical one. Showing a family the child's SEGA as a small object that the model expects to approach, but not necessarily cross, a discussion threshold over the next year — with the uncertainty visible as a soft cloud rather than a falsely precise point — is a more honest and more humane conversation than reciting millimeters. For trainees, the twin is a teaching object that ties the abstract surveillance schedule to concrete anatomy.

### 2.4 Prior art and why TSC is the right first twin

"Digital twin" is an overloaded term. In engineering it denotes a live, physics-simulating mirror of a physical asset. In medicine the term spans patient-specific computational physiology, organ-on-chip mirrors, and — most relevant here — *data-driven anatomical visualizations* that track a real patient's measured state over time. Radiology has long used 3D reconstruction; surgical planning uses patient-specific 3D models; cardiology uses anatomical twins for device planning. The TSC twin sits deliberately at the *data-driven visualization* end of this spectrum, not the physiology-simulation end (Part II makes this scoping explicit and honest).

TSC is an unusually good first target for three reasons. First, its management is lesion- and surveillance-centric, so the mapping from engine output to renderable geometry is direct. Second, the engine already produces *calibrated forecasts with uncertainty*, which is precisely the signal that 3D usually lacks and most needs. Third, the disease is multisystem, so a whole-child twin has immediate communicative value that a single-organ render would not. The same properties that made TSC the right first *engine* make it the right first *twin*.

---

# Part II — Scientific foundations

## 3. What "digital twin" means here — and what it does not

### 3.1 A data-driven scene, not a physiological simulation

We define the TSC Anatomical Digital Twin precisely, to forestall overclaim: it is a **deterministic, atlas-anchored, time-resolved 3D rendering of the engine's structured projections**. It is *not* a finite-element growth simulation, a biomechanical model, or a generative reconstruction of patient anatomy. It does not predict anything the engine has not already predicted; it *displays* the engine's predictions in three spatial dimensions plus time plus uncertainty. Every visible quantity traces to an analytic output that is itself validated and provenance-stamped.

This scoping is a feature, not a limitation. The danger of medical 3D is exactly the seductive implication that a render carries more information than its inputs. By defining the twin as a pure, deterministic *function of the projection* — same projection in, byte-identical scene out — we guarantee that the twin can make no claim the engine cannot defend. The render is a lens, not a source.

### 3.2 The data substrate already exists

The single most important technical observation in this document is that **the engine's projections are already geometry-ready**. The trajectory modeler emits, per lesion:

- a label and anatomical location (e.g., SEGA at the foramen of Monro; AML in the kidney);
- the observed size series in centimeters with their months;
- a forecast at 6/12/18 months, each with a mean and 50%/90% prediction intervals;
- a clinical threshold in centimeters, a time-to-threshold, and flags for "crosses in the 12–18-month window" and "at or above threshold."

The variant curator emits the recovered variant, its VAF, its mosaic status, and its ACMG classification with per-criterion evidence. The phenome mapper emits HPO-coded, ontology-validated phenotypes, each carrying an organ implication. Every one of these is exactly the input a scene needs. No new computation is required to author the twin — only a transformation from numbers to prims. The intermediate representation that performs this transformation (the *SceneSpec*, specified in the PRD) is engine-agnostic and small.

### 3.3 Uncertainty visualization as the core scientific contribution

The scientifically interesting part of the twin is not the lesion geometry; it is the *uncertainty geometry*. The engine's trajectory forecasts are produced by a classical pipeline — an OLS trend with a mixed-effects population prior shrinking the slope for sparse SEGA series, and a Gaussian-process posterior for interval width — and the resulting 50%/90% prediction intervals are the engine's honest statement of what it does and does not know.

Rendering that uncertainty volumetrically is a deliberate design choice grounded in the uncertainty-visualization literature, which consistently finds that explicit, perceptually-scaled uncertainty representations reduce overconfidence relative to point estimates. We map the prediction interval to a pair of nested, translucent shells around the forecast lesion: the inner shell at the 50% interval, the outer at the 90%. The shells are *not* decorative — their radii are computed from the same interval bounds the harness validates. When the engine is confident (tight interval, e.g., a well-sampled slow-growing AML), the cloud is thin and crisp. When the engine is uncertain (a sparse SEGA series leaning on the population prior), the cloud is large and soft. A clinician reading the scene reads the engine's confidence *directly off the geometry*, with no number to misinterpret.

This is the inversion that makes the twin safe: in most medical 3D, fidelity of rendering implies precision of knowledge, and the two are decoupled, which is dangerous. Here they are *coupled by construction* — the render is exactly as crisp as the engine is certain.

### 3.4 Anatomical grounding: the honest line

The twin's anatomy is **atlas-based and stylized**. We anchor each lesion to a canonical anatomical landmark in a stylized pediatric reference geometry (a low-polygon brain with ventricular landmarks including the foramen of Monro; stylized kidneys; a body shell for the whole-child view) and scale the lesion to the engine's measured size. We do *not* claim that the rendered brain is this child's brain. The synthetic cohort has no imaging, and even with imaging the v0.x twin renders measured-size lesions on reference anatomy, not segmented patient anatomy.

The line to patient-specific anatomy is clear and is Phase-1 institutional work: real DICOM imaging → automated segmentation (e.g., a MONAI-based pipeline) → per-patient organ and lesion meshes → USD. When that path is built, the twin's *scene schema does not change* — only the source of the geometry. The atlas meshes are swapped for segmented meshes; the lesion scaling, the uncertainty envelope, the threshold surface, the animation, and the provenance all remain. This is the same "swap the box labels, keep the wiring" discipline that governs the engine.

### 3.5 Provenance travels into the scene

The engine's defining discipline is that every output is traceable — model, tier, latency, prompt-template version, input hash. The twin extends this discipline into three dimensions rather than abandoning it. Each rendered element carries, as USD custom metadata, the provenance of the projection field that produced it: the lesion geometry references the trajectory event and its model record; the variant prim references the curation event and the ACMG rule; a dedicated provenance scope in the scene mirrors the audit trail. A reviewer can select a glowing lesion and read *why it is that size, from which model, at what latency, under which prompt version*. The twin is not a marketing render that hides its sources; it is an auditable artifact that surfaces them spatially.

---

# Part III — Technology and architecture

## 4. OpenUSD, Omniverse, and the authoring/render split

### 4.1 Why OpenUSD is the right substrate

OpenUSD (Universal Scene Description) is the interchange and scene-composition system originated by Pixar and now the backbone of NVIDIA Omniverse. It is the correct substrate for the twin for reasons that are technical, not fashionable:

- **It is a composition system, not just a file format.** USD's layering, references, sublayers, and variant sets let us compose a scene from independent pieces — a reusable anatomical atlas layer, a per-patient lesion layer, a provenance layer, a watermark layer — and override or swap any one without touching the others. The Phase-1 swap from atlas to segmented anatomy is a sublayer replacement.
- **It encodes animation natively as time-samples.** A USD attribute can hold values at multiple time codes; the renderer interpolates. The lesion's growth across observed months and into the forecast is encoded as time-sampled radius and color attributes — the trajectory *is* the animation, with no bespoke animation system.
- **It carries arbitrary metadata.** USD `customData` lets every prim hold the provenance dictionary, so the audit trail lives inside the scene.
- **It is the lingua franca of the target ecosystem.** Omniverse, USD Composer, Kit applications, and USDZ-for-AR all consume USD. Authoring USD targets all of them at once.

### 4.2 The authoring/render split — the key architectural decision

The central architectural decision is to **separate scene authoring from scene rendering**, and to put authoring where the engine already lives:

| Stage | Where | Needs a GPU? | Output |
| --- | --- | --- | --- |
| Engine analytics | DGX Spark (built) | Only for the genomic substrate burst | Projections (the engine's existing output) |
| Scene authoring | DGX Spark (this layer) | **No** | A deterministic `.usd`/`.usda` scene |
| RTX rendering | RTX workstation or RunPod | Yes (RTX) | Real-time interactive view; path-traced stills/films; USDZ for AR |

Authoring a USD scene is a pure data transformation — it constructs prims, sets attributes, writes time-samples — and runs comfortably on the Spark's CPU. The `usd-core` Python package (the `pxr` modules, CPU-only, cross-platform with ARM wheels) provides the full authoring API; where even that dependency is unwanted, USD's `.usda` form is a well-specified ASCII text format the layer can emit directly. Either way, **the Spark produces the scene with no new GPU and no Omniverse install.** Rendering — the GPU-bound, RTX-accelerated step — happens wherever an RTX GPU is available: an RTX workstation running USD Composer, or a burst RunPod instance. This mirrors the engine's own "single Spark, burst for GPU-heavy steps" thesis exactly.

This split is what makes the proposal credible rather than aspirational. The valuable, differentiating artifact — a correct, provenance-stamped, uncertainty-aware scene — is produced deterministically and freely on hardware Adam already owns. The render is the finish, not the foundation.

### 4.3 Omniverse, RTX, and interactivity

On the render side, NVIDIA Omniverse provides RTX rendering in two modes — a real-time rasterization-plus-ray-tracing path for interactive navigation and an interactive path-traced mode for reference-quality images and film. Omniverse's MDL material system gives physically-based translucency that makes the uncertainty envelope read as a true volumetric cloud and the threshold surface as a glowing membrane. Interactivity — the timeline scrubber that drives the lesion from −24 months into the +18-month forecast — is available at three levels of investment: the native USD timeline (free, immediate), an OmniGraph/Action-Graph behavior for richer in-scene controls, or a small Omniverse Kit extension that adds a purpose-built TSC panel (patient picker, scrub, organ toggles). The PRD phases these so that value arrives early (native timeline) and polish arrives later (Kit extension).

### 4.4 Hardware reality, stated plainly

Two honest constraints govern deployment. First, native Omniverse Kit support on the Grace-Blackwell (GB10) DGX Spark should be **verified, not assumed**; the architecture deliberately does not require it, because authoring is CPU-side and rendering can be off-box. Second, the twin's anatomy is reference geometry until the Phase-1 imaging pipeline lands. Neither constraint blocks the high-value first build, and both are surfaced rather than buried — the same discipline that lets the engine claim only what it can defend.

---

# Part IV — The three signature scenes

## 5. From projection to spectacle, without losing honesty

The layer ships three scene templates, each a deterministic function of the engine's output and each carrying the SYNTHETIC watermark and full provenance. They are ordered by build priority, which is also their order of clinical value.

### 5.1 Scene 1 — the volumetric lesion-trajectory twin (the centerpiece)

This is the scene that justifies the layer. For a selected patient and lesion (Patient B's SEGA at the foramen of Monro; Patient C's AML approaching 4 cm), the scene renders:

- the **lesion** as an atlas-anchored gprim (a sphere or ellipsoid) whose radius is scaled to the engine's measured size, positioned at the anatomical landmark in stylized reference anatomy;
- a **growth animation** driven by USD time-samples across the observed months (−24, −12, −6, now) and forward to the 6/12/18-month forecast means — the lesion literally grows along the engine's trajectory as the timeline scrubs;
- the **uncertainty envelope**: two nested translucent shells whose radii at each future time code are the engine's 50% and 90% prediction-interval bounds, so the forecast lesion sits inside a cloud that widens exactly as the engine's confidence narrows;
- the **threshold surface**: a glowing translucent shell at the clinical threshold radius (SEGA discussion threshold; AML ~4 cm bleeding-risk threshold), with the lesion's color transitioning through a colorblind-safe palette as it moves from *below* to *approaching* to *at/above* threshold;
- a **time-to-threshold marker** and the "12–18-month window" highlighted on the timeline, so the *preemptive intervention window* is a visible region, not a footnote.

Scrubbing this scene from the past into the forecast is the single most powerful demonstration the engine can give: a clinician watches a child's lesion approach a decision boundary, sees exactly how confident the model is by the size of the cloud, and reads the provenance of every element on selection. It is the two-dimensional growth chart made anatomical, volumetric, temporal, and — through the envelope — *honest*.

### 5.2 Scene 2 — the mosaic-recovery "powers of ten"

This scene dramatizes the engine's headline capability: recovery of low-VAF somatic mosaic variants that standard testing misses (the engine measures 100% recovery of the six sub-threshold mosaics in the synthetic cohort, a +12-point detection uplift). It is a continuous scale traversal:

- begin at the **organ** (the resected tuber tissue for Patient A);
- zoom to a **cellular field** in which only the recovered VAF fraction of cells carries the variant — for Patient A's 8.3% mosaic, roughly one cell in twelve glows, making "8.3% VAF" a thing you can *count* rather than a number to interpret;
- continue into the **genome** and onto the **variant** — the truncating change in *TSC2* — with the ACMG criteria (PVS1 applied Strong for the low-VAF mosaic, PM2, PP4) and the resulting Likely-Pathogenic call surfaced as annotation.

Where Scene 1 makes the forecast tangible, Scene 2 makes the *diagnostic odyssey* tangible: it shows, in one shot, the signal a blood test returns as negative and the engine recovers from tissue. It is the visual argument for the detection-uplift number.

### 5.3 Scene 3 — the whole-child organ atlas and the population view

This scene answers "is this a whole disease or a single lesion?" and "is this one patient or a population?":

- a **whole-child body shell** with each TSC-affected organ illuminated from the patient's HPO profile (brain tubers/SEGA, kidney AML, cardiac rhabdomyoma, skin features, pulmonary involvement), so the multisystem burden is legible at a glance and the surveillance schedule can be overlaid as a "what to scan, when" annotation keyed to the engine's ITSC gap analysis;
- a **population mode** in which the cohort's fifty patients appear as an array of ghosted bodies forming a lesion-burden and recovery heatmap — the engine's population/command-center view, rendered in three dimensions, with the seven recovered mosaics highlighted.

Scene 3 ties the twin back to the engine's cohort view and to the platform's scale thesis: one Spark, fifty patients, a whole disease.

---

# Part V — Clinical integration, safety, and validation

## 6. Entering the visit safely

### 6.1 Surface (d): the fly-in from the dashboard

The twin is **Surface (d)**, joining the pre-visit briefing, the in-visit dashboard, and the async alerts. Its natural entry is a "fly-in" from the in-visit dashboard: a clinician viewing Patient B's two-dimensional trajectory chart taps through to the volumetric twin of the same lesion, scrubs the forecast, and returns. The twin never originates a clinical claim; it visualizes a projection the dashboard already shows. This keeps the twin a *view*, consistent with the engine's principle that surfaces assemble from projections and never compute new conclusions.

### 6.2 Safety: the disciplines that do not relax in 3D

Every safety discipline the engine enforces applies in the twin, and one new discipline is added:

- **SYNTHETIC watermark, always.** A persistent in-scene watermark prim and a render-overlay HUD mark every frame as synthetic demonstration data; the watermark is part of the scene, not a removable layer.
- **Decision support, never autonomous; never a diagnosis from a render.** The twin shows forecasts and classifications the engine produced under clinician-review framing; it does not diagnose, and it carries the same "clinician review required, not FDA-cleared" framing as every other surface.
- **No implied patient-specific anatomy.** The stylized atlas is visually distinct from a patient reconstruction, and the scene is labeled as atlas-anchored, so no viewer mistakes reference geometry for this child's scan.
- **The new discipline — uncertainty is non-optional.** A forecast lesion may never be rendered *without* its uncertainty envelope. Showing the mean alone would imply false precision; the envelope is a required element of any forecast scene, enforced in the scene schema. This is the twin's most important safety property and the inversion described in §3.3.

### 6.3 What the twin must not do

To keep the layer trustworthy we enumerate anti-requirements: it must not render a forecast point estimate without its interval; it must not imply segmented patient anatomy it does not have; it must not introduce any quantity not present in the projection; it must not allow the watermark or provenance to be stripped in the authored scene; and it must not let visual drama outrun the engine's certainty — the palette and envelope are calibrated to the numbers, not chosen for effect.

### 6.4 Validation: the render must equal the math

Because the twin is a deterministic function of the projection, it is *testable as software*, and the validation harness extends naturally:

- **Geometry-scale fidelity.** The rendered lesion radius equals the engine's measured size under the scene's cm-to-scene-unit scale, to tolerance.
- **Envelope-equals-interval.** The 50%/90% shell radii at each forecast time code equal the engine's prediction-interval bounds — the visual uncertainty *is* the numeric uncertainty, asserted by test.
- **Threshold-state correctness.** The lesion's color state (below/approaching/at-or-above) matches the engine's threshold flags at every time code.
- **Provenance completeness in USD.** Every renderable element carries its provenance metadata, mirroring the engine's 100% provenance-completeness metric.
- **Determinism.** The same projection authors a byte-identical scene, so scenes are diffable and reviewable like code.

Beyond software tests, the layer's *clinical* validation is a forward question — does the twin change decision confidence, decision time, or communication quality relative to the 2D surfaces? — to be studied with clinicians as a usability and decision-impact evaluation. As with the engine, construct-level correctness is demonstrable now; prospective clinical-impact validation is institutional work, and the document does not blur the two.

---

# Part VI — Generalization, roadmap, and references

## 7. Generalization: one exporter, many diseases

The twin generalizes for the same reason the engine does: it is parameterized over a small, swappable surface. A different disease supplies a different anatomical atlas and a different lesion-and-organ vocabulary; the SceneSpec, the time-sampled animation, the uncertainty envelope, the threshold surface, the provenance, and the watermark are unchanged. NF1/NF2 (nerve-sheath and CNS tumors with their own growth-and-threshold decisions), Rett, Williams, and the broader mTORopathies are configuration targets, not rewrites — the same replication thesis the engine's master volume sets out for the analytic layer, now extended to the visualization layer. A platform that can *see* one rare disease in volumetric, uncertainty-aware time can see the next one by swapping the atlas.

## 8. Roadmap (summary; the PRD carries the detail)

- **Phase V0 — CPU-side authoring, buildable now.** The SceneSpec intermediate representation, the OpenUSD exporter (pxr where available, ASCII-USD fallback otherwise), and Scene 1 (the lesion-trajectory twin) with its uncertainty envelope, threshold surface, animation, provenance, and watermark — all authored and tested on the Spark with no new GPU, validated against the engine's projections.
- **Phase V1 — RTX rendering and the remaining scenes.** Omniverse/USD Composer setup on an RTX target; MDL materials; Scenes 2 and 3; the timeline-driven scrub as a polished experience.
- **Phase V2 — AR and portal preview.** USDZ packaging for tablet-based AR in clinic; a lightweight in-portal web preview (a pre-rendered turntable or a constrained web-3D view) that degrades honestly when Omniverse is not present.
- **Phase V3 — patient-specific anatomy (institutional).** The DICOM → MONAI segmentation → USD pipeline that replaces atlas geometry with segmented patient anatomy, changing the source of meshes while leaving the scene schema intact.

## 9. Closing

The TSC Intelligence Engine earned its credibility by being honest about its edges and by turning its claims into measured numbers. The Anatomical Digital Twin is the natural next surface precisely because it can inherit that discipline rather than dilute it: a render whose crispness equals the engine's certainty, whose every element traces to a validated projection, whose watermark and provenance cannot be stripped, and whose first and most valuable build runs for free on the hardware already in hand. It does not make the engine claim more. It makes what the engine already knows — the future of a child's lesion, with calibrated uncertainty — something a clinician, a trainee, and a family can see together. That is how a rigorous decision-support engine becomes, in the room, genuinely unforgettable.

## 10. References and notes

This document is a design rationale; its empirical anchors are the engine's own measured outputs (see the master research volume and the validation harness) and the established clinical framework for TSC surveillance (the 2021 International TSC Consensus recommendations). Its technical claims rest on the public OpenUSD specification and the NVIDIA Omniverse / RTX and OpenUSD documentation, and on the uncertainty-visualization literature's consistent finding that explicit, perceptually-scaled uncertainty representations mitigate overconfidence relative to point estimates. Specific software interfaces (the `pxr` USD authoring API, USD time-samples and `customData`, MDL materials, USDZ packaging, OmniGraph) are cited in the companion PRD at the level of the requirements that depend on them. All patient data referenced is synthetic; no real patient data, imaging, or institutional commitment is implied.
