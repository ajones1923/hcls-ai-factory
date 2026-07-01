# The TSC AI Engine in NVIDIA Omniverse — A Visual Foundations Guide (Advanced Edition)

### What the 3-D visuals look like, what they show, and how they are built — explained in depth

**Author:** Adam M. Jones — HCLS AI Factory, TSC Intelligence Engine (Engine 7)
**Reading level:** written for a strong high-school senior or first-year college reader (roughly 12th grade / AP). It assumes general computer literacy and a little geometry, and it defines specialized 3-D and rendering terms as they appear. A full glossary closes the guide.
**Companion to:** *The TSC AI Engine — A Foundational Learning Guide (Advanced Edition)*, the *TSC Anatomical Digital Twin* research paper, and the Digital Twin PRD.
**License:** Apache 2.0 (open source) · 2026

> **How to use this guide.** Parts 1–3 explain *why* the engine visualizes its results and *what kind*
> of 3-D object this is (and is not). Parts 4–8 walk scene by scene through what you actually see on
> screen and what each view shows clinically. Parts 9–11 cover how the images are built and rendered —
> the technology, the materials, and the computer pipeline. Parts 12–15 cover the demo experience, how the
> visuals are tested, the roadmap, and the honest limits. **Read this first:** every scene is generated
> from *synthetic* (computer-made) patient data, the anatomy is *stylized reference geometry* and not a
> picture of any real person's scans, and the system is clinician **decision support** that is **not
> FDA-cleared.** As you will see, that honesty is engineered directly into the pictures.

---

## Part 1 — Why visualize at all? The cognitive case

The TSC AI Engine already produces precise output — a recovered genetic variant, a growth forecast with
confidence ranges, a phenotype map, a list of overdue scans. Today that output is shown as cards, tables,
and a two-dimensional growth chart. That is correct and auditable. It is also, for the specific decisions
clinicians make, *under-dimensioned*.

The reason is that Tuberous Sclerosis Complex (TSC) is a fundamentally **spatial and temporal** disease.
Its growths sit in specific anatomical locations and change over time, and the core decisions are about
*when* to act and *with what confidence*. Consider three recurring decisions:

- **A brain tumor (SEGA):** *when* will it reach the size where the team discusses surgery, *where* is it
  relative to the narrow channel (the foramen of Monro) it can block, and *how sure* are we?
- **A kidney growth (AML):** is this lesion merely near, or genuinely crossing, the ~4-centimeter bleeding-
  risk threshold — and is it accelerating?
- **Surveillance intensity:** which organ systems are involved, and does the whole-body burden justify
  scanning more often than the default schedule?

Each of these is a function of three things at once — spatial extent, a time trajectory, and calibrated
uncertainty — and a clinician answers them by constructing a **mental model**: a spatial, growing,
uncertain object in time, assembled in the mind from numbers on a page. That mental work is skilled but
*invisible*, error-prone under time pressure, and — crucially — *unshareable*. A parent in the room, a
trainee, a colleague in another specialty cannot see the clinician's mental model.

A rendered 3-D twin **externalizes** that model. It turns the engine's already-computed object into
something a clinician, a care team, and a family can see, scrub through time, and discuss together. The
purpose of the visuals is not decoration; it is to match the *form* of the output to the *form* of the
decision.

---

## Part 2 — What a "digital twin" means here (and what it does not)

The phrase "digital twin" is overloaded, so the engine defines it precisely to avoid overclaiming.

**Here, a digital twin is a deterministic, atlas-anchored, time-resolved 3-D rendering of the engine's
structured results.** Three words in that sentence carry weight:

- **Deterministic** — the scene is a pure function of the engine's output. The same results always produce
  the byte-for-byte identical scene. There is no randomness and no separate "creative" step.
- **Atlas-anchored** — lesions are placed on a *stylized reference body* and scaled by the engine's real
  measurements. The brain you see is a low-detail standard brain, not a reconstruction of *this patient's*
  brain.
- **A rendering of results** — the twin displays only what the engine already computed and validated. It
  introduces no new quantity.

Equally important is what the twin is **not**:

- It is **not** a physics or biology *simulation*. It does not model how a tumor actually grows; it
  *displays* the engine's statistical forecast of growth.
- It is **not** a patient-specific anatomical *reconstruction*. The synthetic patients have no scans, and
  even with scans the current version renders measured-size lesions on reference geometry, not segmented
  patient organs.

This scoping is a feature, not a limitation. The central danger of medical 3-D is that a beautiful render
silently implies more knowledge than its inputs contain. By defining the twin as a deterministic function
of validated results — **a lens, not a source** — the engine guarantees the picture can make no claim the
analytics cannot defend. This is enforced as a rule the engineers call the **no-new-claim invariant**:
every number that appears in a scene must trace back to a field the engine produced.

---

## Part 3 — The foundational technologies, explained

To understand how the engine turns numbers into navigable 3-D, you need four building blocks. None require
prior 3-D experience.

### 3.1 OpenUSD — the scene as structured data

**OpenUSD (Universal Scene Description)**, originated at Pixar and now central to NVIDIA Omniverse, is the
format the scenes are written in. It is less a "picture file" and more a structured description of a 3-D
world. A few concepts matter:

- A **prim** (primitive) is an object in the scene — a sphere, a capsule, a group, a light. Prims are
  organized in a tree, like folders.
- An **attribute** is a property of a prim — its radius, its color, its position.
- **Time-samples** are how OpenUSD encodes *animation*: an attribute can hold different values at different
  **time codes** (think frame numbers), and the renderer smoothly interpolates between them. This is how a
  lesion can *grow* as you scrub a timeline — its radius attribute is time-sampled.
- **customData** is a slot for arbitrary metadata attached to any prim — used here to carry the watermark
  and the full provenance (which model produced a value, when, and how).
- A scene can be written as **.usda** (human-readable text) or a binary form; the engine emits text, which
  makes the scenes diff-able and reviewable like source code.

### 3.2 NVIDIA Omniverse and RTX rendering

**NVIDIA Omniverse** is the platform that opens and renders OpenUSD scenes with **RTX** — NVIDIA's
ray-tracing technology. Ray tracing simulates how light actually travels, bounces, and passes through
materials, producing photorealistic results: true glass, soft glows, realistic shadows. Omniverse offers a
fast interactive mode for navigation and a reference-quality **path-traced** mode for film-grade stills and
video. This is what makes the uncertainty cloud look like real translucent glass and the variant cells
look like they genuinely glow.

### 3.3 MDL — physically based materials

A **material** tells the renderer how a surface interacts with light — how shiny, transparent, or glowing
it is. NVIDIA's **MDL (Material Definition Language)** provides standard materials such as **OmniGlass**
(true glass) and **OmniPBR** (a general physically-based material that supports *emission*, i.e., glow).
The engine assigns MDL materials to specific elements so they render with cinematic fidelity.

### 3.4 The authoring/render split — why it stays cheap

The single most important practical idea is that **building a scene and rendering it are separate steps on
separate machines:**

- **Authoring** — converting the engine's results into an OpenUSD scene — is pure data manipulation. It
  runs on the small, affordable **NVIDIA DGX Spark** with **no GPU required** and no special software (the
  engine can even write the OpenUSD text directly, with zero third-party dependencies).
- **Rendering** — the GPU-intensive, RTX-accelerated step — runs separately on a rented cloud GPU via a
  service called **RunPod**.

In one sentence: **the Spark thinks and authors; the cloud GPU renders.** This keeps the differentiating,
trustworthy artifact — a correct, provenance-stamped scene — free to produce on hardware already on hand,
while the expensive beauty pass happens only when needed.

---

## Part 4 — The bridge: turning engine results into geometry

Before any scene is drawn, the engine's results are transformed into a small, neutral description called a
**SceneSpec** — an intermediate representation that holds everything a scene needs as plain data, with no
3-D or rendering concept in it. Two writers then consume the SceneSpec: one using the full OpenUSD software
library, and a dependency-free one that writes the OpenUSD text directly so the system works on a bare
machine. The mapping from results to geometry is direct and deterministic:

| Engine result | Becomes, in the scene |
| --- | --- |
| A lesion's measured size (cm) | A sphere's radius |
| Observed months and forecast horizons | Time codes on the timeline (e.g., −24 … +18) |
| The 50% / 90% prediction intervals | The radii of the two uncertainty shells |
| The clinical threshold | A translucent threshold "membrane" |
| The crossing risk (likely/possible/unlikely) | The lesion's color state |
| Phenotype (HPO) profile | Which organs light up in the body view |
| ACMG classification | A figure's body color in the population view |
| "Recovered mosaic" flag | A gold halo around the figure |
| Provenance (model, tier, latency) | customData metadata on the prim |

Because the mapping is fixed and the math comes straight from the engine, the scenes are reproducible and
testable like software — a point Part 13 returns to.

---

## Part 5 — Scene 1: the lesion-trajectory twin (the centerpiece)

This is the scene that justifies the whole layer. It takes one patient and one lesion and renders its past,
present, and forecasted future as a single, scrubbable 3-D object. The canonical example is **Patient B's
SEGA** at the foramen of Monro, measured at **0.8, 1.1, and 1.3 cm** over the past two years, forecast to
about **1.9 cm**, against a discussion threshold of **1.8 cm**.

**What you see on screen.** At the center sits the lesion as a sphere, anchored at its anatomical landmark
inside a faint, stylized brain for context. A timeline runs from frame −24 (two years ago) through 0 (now)
to +18 (eighteen months into the forecast). Three more elements surround the lesion:

- **The growth animation.** As you scrub the timeline forward, the lesion's radius is driven by time-samples
  — it is small at −24, grows through the observed measurements, and continues *along the engine's forecast*
  into the future. You literally watch it grow at the rate the model predicts.
- **The uncertainty envelope.** Two nested, translucent shells surround the forecast lesion — an inner shell
  and an outer shell. Here is the defining, honest design choice: **the radii of those shells are set
  exactly equal to the engine's 50% and 90% prediction intervals.** The cloud is, geometrically, the
  forecast's uncertainty. When the engine is confident, the shells hug the lesion tightly. When it is
  unsure, they balloon into a large, soft cloud. The render is *exactly as crisp as the engine is certain.*
- **The threshold membrane.** A glowing, semi-transparent shell marks the clinical threshold (1.8 cm for a
  SEGA; ~4 cm for a kidney AML). As the lesion grows toward it, the lesion's **color shifts** through a
  colorblind-safe palette — teal while comfortably below, amber as it approaches, vermillion at or beyond —
  so threshold proximity reads instantly, even without reading a number.

**What it shows clinically.** Scrubbing this scene from the past into the future is the most powerful single
demonstration the engine can give. A clinician watches a child's lesion approach a decision boundary, reads
the engine's confidence directly off the *size of the cloud*, and — because the forecast window (months 12
to 18) is highlighted on the timeline — sees the **preemptive-intervention window** as a visible region
rather than a footnote. The two-dimensional growth chart, made anatomical, volumetric, temporal, and honest.

The same scene type serves other patients and quantities. **Patient C** has a renal **AML** sitting right at
its 4 cm threshold (rendered at-or-beyond, vermillion), alongside **declining kidney function (eGFR)** — a
*decreasing* quantity whose threshold is a floor it falls toward, which the geometry handles by drawing the
"toward-threshold" interval bound — and **rising seizure frequency.** Each monitored quantity can be shown
as its own trajectory twin.

---

## Part 6 — Scene 2: the mosaic "powers-of-ten" scene

If Scene 1 makes the *forecast* tangible, Scene 2 dramatizes the engine's flagship diagnostic capability:
recovering the low-level mosaic variants that standard blood testing misses. The canonical example is
**Patient A** — a *TSC2* frameshift variant recovered at about **8.3% VAF** in tissue, reported "no mutation
identified" on blood.

**What you see on screen.** A field of cells — a grid of 144 in the current build — of which **exactly the
recovered VAF fraction carries the variant.** For 8.3%, that is **12 cells in 144 — about one in twelve —
glowing gold**, while the rest sit quiet and dim. The variant-carrying cells are rendered with an *emissive*
material so they genuinely glow under RTX. Riding along as metadata is the full variant call: the gene, the
specific change, the ACMG classification (Likely Pathogenic), the criteria applied, and the fact that it was
*recovered from affected tissue.* (Technically, the cell field uses **point-instancing** — an efficient
OpenUSD mechanism for placing many copies of a few prototype shapes — so even large fields render smoothly.)

**What it shows clinically.** It makes an abstract statistic *countable*. "8.3% VAF" stops being a number to
interpret and becomes a thing you can literally count on screen: roughly one cell in twelve. It shows, in a
single image, the exact signal a blood test returns as negative and the engine recovers from tissue — the
visual argument for the measured +12-point detection uplift and the end of a diagnostic odyssey. In the
demonstration's narrative, this is the Act One showpiece.

---

## Part 7 — Scene 3: the whole-child organ atlas

This scene answers a different question: *is this a single lesion, or a whole multisystem disease?*

**What you see on screen.** A stylized pediatric body (a translucent capsule torso with placed organ
markers — brain, heart, lungs, kidneys, skin). Each organ system is **illuminated based on the patient's
phenotype profile**: organs the engine's HPO map shows as involved glow brightly in their distinct color,
while uninvolved organs sit faint and translucent. For Patient B, for instance, the brain, kidney, and skin
light up. Any **overdue surveillance** items ride along as metadata, and the patient's classification and
organ burden are recorded on the figure.

**What it shows clinically.** The multisystem nature of TSC, legible in a single glance, rather than as a
list of phenotype codes. It is a natural communication object for a family conversation and a teaching
object for trainees — and it sets up the surveillance discussion (which organs, how often).

---

## Part 8 — Scene 4: the population view

The final scene answers the scale question: *is this one impressive patient, or a whole cohort?*

**What you see on screen.** All **50 synthetic patients** arranged as a grid of whole-child figures. Each
figure's **body color encodes its ACMG classification** (a distinct color for Pathogenic, Likely Pathogenic,
Variant of Uncertain Significance, and No-variant), its organs are lit by that patient's phenome, and —
the centerpiece — the **seven recovered mosaic patients are each ringed in a glowing gold halo**. A legend
records the convention, and the cohort's classification distribution rides along as metadata.

**What it shows clinically and strategically.** It is the engine's scale story rendered in three dimensions:
one small computer, fifty children, a whole disease — and the seven that a standard pipeline would have
missed, glowing gold across the room. It ties the cinematic single-patient scenes back to the cohort-level
evidence (the validation scorecard) and makes the "this scales" claim tangible. In the demonstration, this
is the Act Three closer.

---

## Part 9 — Materials and the look (MDL, RTX) — and one careful exception

The four scenes are designed to render at film quality under RTX, using MDL materials assigned to specific
elements:

- **The uncertainty envelope** → **OmniGlass** (true glass). The forecast cloud reads as a real translucent
  volume, not a flat tint.
- **The threshold membrane** → a red glass material, so the decision boundary glows as a surface.
- **The variant-carrying cells** (Scene 2) → **OmniPBR emissive** gold. They genuinely emit light, selling
  the "glowing one-in-twelve" effect.
- **The recovery halos** (Scene 4) → emissive gold, translucent, so each recovered patient is ringed in a
  soft glow.

**The careful exception — keeping the visuals honest.** One element is *deliberately left without an MDL
material*: the lesion in the trajectory twin. The lesion's *color* is time-sampled — it changes from teal to
amber to vermillion as the lesion crosses the threshold over the timeline — and binding a fixed material
would freeze that animation. So the lesion keeps its engine-driven, animated color. This is a small but
telling example of the project's discipline: the dramatic, threshold-crossing color change is information
from the engine, and it is protected from being overridden by a cosmetic material.

For portability, every element also keeps a plain display color as a fallback. The MDL modules resolve and
shine inside Omniverse; in a generic 3-D viewer without them, the scenes still show correct colors. Best of
both: cinematic in Omniverse, universally openable elsewhere.

---

## Part 10 — The honesty disciplines, built into the geometry

The visuals inherit — and physically embody — the engine's honesty principles. Five are worth naming because
each is a concrete property of the scenes, not a disclaimer.

**1. The uncertainty envelope is mandatory.** A forecast lesion may never be authored *without* its
envelope. Showing the mean alone would imply false precision; the envelope is a required element of any
forecast scene. This is the layer's single most important safety property — and it produces the inversion
that makes medical 3-D safe here: in most renders, visual polish implies knowledge, and the two are
dangerously decoupled. Here they are *coupled by construction* — the render is exactly as crisp as the
engine is certain.

**2. The watermark is part of every scene.** Each scene carries a persistent "SYNTHETIC demonstration data —
decision support — not FDA-cleared" mark in its metadata. It is authored into the scene, not a removable
overlay.

**3. Provenance travels into 3-D.** Each rendered element carries, as customData, the provenance of the
result that produced it — which model, which tier, how long it took. The audit trail does not stop at the
2-D dashboard; it follows the data into the scene, so a reviewer can select a glowing lesion and read *why*
it is that size, from which model.

**4. The anatomy is labeled as what it is.** Every scene states that the anatomy is stylized, atlas-anchored
reference geometry — *not* patient imaging — so no viewer mistakes the reference brain for this child's scan.

**5. Determinism makes the scenes testable.** Because a scene is a pure function of the engine's results, the
same results produce a byte-identical scene. Scenes can therefore be reviewed and *tested* like source code
— which leads directly to Part 13.

---

## Part 11 — The pipeline in practice: Spark → RunPod → Omniverse

Putting it together, the end-to-end workflow a demonstrator follows is:

1. **Author on the Spark.** A single command (or an API call) reads the engine's results and writes the
   OpenUSD scenes to disk — CPU-only, offline, deterministic, free. This is the entire authoring step.
2. **Move the scenes to a RunPod RTX pod.** The lightweight scene files (text OpenUSD) are copied to a
   cloud machine with an RTX GPU running NVIDIA Omniverse / USD Composer.
3. **Open and render.** The scene is opened, the renderer is set to RTX path tracing, and — for the
   trajectory twin — the timeline is scrubbed from −24 to +18 months. The glass envelopes, glowing cells,
   and gold halos come alive.

The division of labor is the whole point: the valuable, trustworthy artifact (a correct, provenance-stamped
scene) is produced for free on hardware already owned; the expensive RTX beauty pass is rented only when a
film-quality render is wanted. Nothing about the engine depends on Omniverse running on the Spark itself.

---

## Part 12 — The visual grammar: how to read any scene

The four scenes share a consistent **visual language**, so that once you learn to read one, you can read
all of them. Every visual property carries a fixed meaning.

**Color encodes state.** Color is never decorative; each color means something specific:

- **Lesion color (teal → amber → vermillion)** encodes threshold proximity: comfortably below, approaching,
  and at-or-beyond the clinical threshold.
- **Body color** in the population view encodes the ACMG classification (a distinct color for Pathogenic,
  Likely Pathogenic, Variant of Uncertain Significance, and No-variant).
- **Gold** always means "recovered" — a recovered mosaic variant: the glowing cells in Scene 2 and the
  halos in Scene 4.
- **Blue glass** is always the uncertainty envelope; **red glass** is always the threshold membrane.

**Size encodes magnitude and confidence.** A lesion's radius is its measured or forecast size. The
envelope's size is the *uncertainty*: a small, tight cloud means the engine is confident; a large, soft
cloud means it is not. You read the engine's confidence directly off how big the cloud is.

**Opacity encodes role.** Solid surfaces are the thing itself (the lesion, an involved organ). Translucent
surfaces are context or uncertainty (the reference anatomy, the envelope shells, the threshold membrane).
Faint, nearly invisible surfaces are *uninvolved* (organs the phenome does not flag).

**Position encodes anatomy and identity.** Lesions sit at their real anatomical anchors (a SEGA at the
foramen of Monro, an AML at the kidney). In the population grid, each figure holds a fixed position, so the
same patient is always in the same place.

**Accessibility is built in.** The palette is **colorblind-safe** — it uses a teal-to-amber-to-vermillion
progression rather than a pure red-green scale, which the most common forms of color blindness cannot
distinguish. And color is **never the only channel**: threshold state is also carried by the lesion's
position relative to the membrane, the envelope's width carries uncertainty independent of hue, and text
labels and metadata back up every visual. This matters because the design rests on a finding from the
uncertainty-visualization research literature: explicit, *perceptually-scaled* representations of
uncertainty (a cloud that literally gets bigger as confidence drops) reduce viewers' overconfidence far
more than a single point estimate does. The envelope is not just pretty; it is the mechanism that keeps the
viewer honest.

**A quick reading key:**

| You see | It means |
| --- | --- |
| A tight, small glass cloud | The engine is confident in this forecast |
| A large, soft glass cloud | The engine is uncertain — treat the forecast with caution |
| The lesion turning amber, then vermillion | It is approaching, then crossing, the clinical threshold |
| A gold glow (cells or a halo) | A mosaic variant was recovered |
| A faint, translucent organ | That organ system is *not* flagged for this patient |

---

## Part 13 — Anatomy of a single frame: scrubbing the SEGA twin

To make the trajectory twin concrete, here is what you actually see as you scrub **Patient B's SEGA** along
the timeline, frame by frame. Recall the data: observed sizes 0.8, 1.1, and 1.3 cm; forecast to ~1.9 cm; a
threshold membrane at 1.8 cm; a crossing graded "likely," in the 12–18-month window.

- **Frame −24 (two years ago):** a small **teal** sphere, 0.8 cm, sitting well *inside* (below) the red
  threshold membrane. No uncertainty cloud — this is an *observed* measurement, so it is shown as a precise
  fact, not a forecast.
- **Frames −12 and −6:** the sphere grows to 1.1 and then 1.3 cm, still teal, still without a cloud — three
  solid past measurements establishing the trend.
- **Frame 0 (now):** the lesion sits at 1.3 cm. This is the present moment; everything to the right is
  forecast.
- **Frame +6:** the lesion grows to roughly 1.6 cm, and — for the first time — the **uncertainty envelope
  appears** around it, modest in size because six months out the engine is still fairly confident. The
  lesion remains teal but is clearly climbing toward the membrane.
- **Frame +12:** the lesion reaches about 1.76 cm and turns **amber** — "approaching." The envelope has
  **widened**, honestly reflecting growing uncertainty further into the future. The highlighted
  12-to-18-month *intervention-discussion window* begins here.
- **Frame +18:** the lesion reaches about 1.9 cm, turns **vermillion**, and **touches or crosses the red
  threshold membrane.** The 90% uncertainty shell now extends out to roughly 2.37 cm — its widest — visibly
  expressing that an 18-month forecast is the least certain.

What a clinician *reads* from this single scrub: not just "it will cross," but *when* (in the highlighted
window), and *how sure* (the cloud is wide at +18, so the timing is a range, not a date). That combination
— the moment of crossing, framed by honest uncertainty — is exactly the information the SEGA-timing
decision requires, and it is delivered in one continuous motion instead of a table of numbers.

---

## Part 14 — Two views, one truth: the 2-D dashboard vs. the 3-D twin

The 3-D twin does not replace the engine's 2-D dashboard; the two are complementary tools for different
moments.

**The 2-D dashboard** is dense, fast, and precise. It is the working clinician's daily instrument: exact
numbers, the four-quadrant layout, the trajectory chart with its shaded bands, the browsable audit table.
When you need to *read values and check sources quickly*, the dashboard wins.

**The 3-D twin** is spatial, temporal, and intuitive. It is the *communication and decision-framing* tool:
the view you turn to when you want to build shared understanding — with a family, a trainee, a colleague in
another specialty — or when the decision itself is fundamentally about a growing object in space and time.
Watching a lesion swell toward a threshold inside a widening cloud conveys the *shape* of a situation in a
way a chart cannot, and the population view conveys *scale* at a glance.

The natural workflow ties them together: a clinician working in the dashboard, looking at Patient B's 2-D
trajectory chart, can **"fly in"** to the 3-D twin of the same lesion, scrub the forecast, and fly back.
The twin is a *lens* applied at specific moments — not a constant replacement for the dense, efficient 2-D
surface. Each is honest about what it is best at: the dashboard for precise, rapid reading; the twin for
spatial intuition, communication, and the demonstration.

---

## Part 15 — What the demo audience experiences (the three acts, visually)

The visuals are sequenced to tell a story that mirrors the clinical argument:

- **Act One — Mosaic recovery.** The mosaic "powers-of-ten" scene: a field of cells where one in twelve
  glows gold, dramatizing the diagnosis a blood test missed and the engine recovered from tissue. Paired
  with the measured +12-point detection uplift.
- **Act Two — The longitudinal twin.** The lesion-trajectory twin, scrubbed from past to future: the SEGA
  growing inside its glass uncertainty envelope toward the red threshold membrane. Paired with the line that
  defines the layer — *the envelope cannot lie; it widens exactly as the model is uncertain.*
- **Act Three — Scale.** The population view: fifty children, the seven recoveries ringed in gold. Paired
  with the validation scorecard and the "one affordable computer, a whole disease" thesis.

The arc moves from a single dramatic recovery, to a single calibrated forecast, to the whole cohort — the
same arc the engine's evidence follows.

---

## Part 16 — Validating the visuals (they are tested like software)

Because each scene is a deterministic function of the engine's results, the visuals are not merely *viewed*
— they are *tested*. The most important test asserts the headline invariant:

> **Envelope equals interval.** The radii of the 50% and 90% uncertainty shells, at every forecast time
> code, are asserted to equal the engine's prediction-interval bounds. The visual uncertainty *is* the
> numeric uncertainty.

Other automated checks assert that the lesion's rendered radius equals the engine's measured size; that the
lesion's color state matches the engine's threshold flags at every time code; that every renderable element
carries provenance metadata; that the watermark is present in every scene; that the mosaic cell fraction
equals the variant's VAF; that the population scene contains exactly the right number of figures and gold
halos; that the MDL materials are bound — *except* the lesion, whose animated color is verified to survive;
and that re-authoring a scene yields a byte-identical file. In short, the same rigor the engine applies to
its analytics applies to its pictures: the visuals make no claim the tests do not back.

What this proves, and what it does not, mirrors the engine's overall honesty: the tests establish that the
*rendering faithfully represents the engine's results.* They do not establish that those results are
clinically validated on real patients — that is the separate, larger study (see the companion guide).

---

## Part 17 — Roadmap of the visual layer

The visualization layer advances in phases, each adding visual capability:

- **Phase V0 — CPU-side authoring (built).** The SceneSpec, the dependency-free OpenUSD writer, and the
  lesion-trajectory twin with its uncertainty envelope, threshold membrane, animation, provenance, and
  watermark — all authored and tested on the Spark with no GPU.
- **Phase V1 — RTX rendering and the full scene set (built).** Omniverse rendering, MDL materials for glass
  and glow, and the mosaic, whole-child atlas, and population scenes.
- **Phase V2 — AR and an in-portal preview (next).** Packaging scenes as **USDZ** for tablet-based
  augmented reality in clinic (so a clinician could place the twin on a table and walk around it), plus a
  lightweight web preview that degrades honestly when Omniverse is not present.
- **Phase V3 — patient-specific anatomy (institutional).** Replacing the stylized atlas with anatomy
  *segmented from a patient's real medical images* (a DICOM-to-segmentation-to-OpenUSD pipeline). Critically,
  the scene structure does not change — only the *source* of the geometry. The lesion scaling, the
  uncertainty envelope, the threshold, the animation, and the provenance all remain; the reference brain is
  swapped for the patient's brain. This is the upgrade that turns "atlas-anchored" into "patient-specific."

---

## Part 18 — Limitations and honest edges of the visuals

A serious guide states its own boundaries plainly.

- **Stylized, not patient-specific.** The anatomy is reference geometry scaled by real measurements. The
  twin shows *this patient's measured lesion on a standard body*, not a reconstruction of their actual
  scans. Patient-specific anatomy is the Phase V3 upgrade and requires real imaging.
- **Synthetic data throughout.** The patients, variants, lesions, and notes are computer-generated. The
  visuals are honest renderings of *demonstration* results, not clinical findings.
- **Rendering is a separate, GPU-dependent step.** The Spark authors the scenes; film-quality RTX rendering
  requires Omniverse on a GPU (RunPod). The MDL "glass and glow" look specifically resolves inside
  Omniverse; a generic viewer shows correct but flatter colors.
- **The standing danger — and its antidote.** Any beautiful medical render risks implying certainty it does
  not have. The engine's countermeasure is structural, not cosmetic: the mandatory uncertainty envelope
  whose size *is* the prediction interval, the in-scene watermark, the provenance metadata, and the explicit
  "not patient imaging / not FDA-cleared" labeling. The visuals are built so that the honest reading is the
  natural one.

Used within those boundaries, the Omniverse layer does something genuinely valuable: it takes the engine's
already-validated, already-honest results and makes them *seeable* — a forecast you can scrub, an 8.3% VAF
you can count, a multisystem burden you can take in at a glance, and a cohort of recoveries glowing gold —
without ever letting the picture outrun the proof.

---

## Glossary

- **Atlas (anatomical)** — a stylized, standard reference body used to place lesions; not a specific
  patient's anatomy.
- **customData** — an OpenUSD metadata slot used here to carry the watermark and provenance on prims.
- **Determinism** — producing a byte-identical scene from the same results every time.
- **displayColor / displayOpacity** — basic OpenUSD color and transparency properties (the portable
  fallback when MDL is unavailable).
- **DICOM** — the standard format for medical images (CT, MRI); the input for future patient-specific
  anatomy.
- **Emissive** — a material property that makes a surface emit (glow with) light.
- **MDL (Material Definition Language)** — NVIDIA's material system; OmniGlass and OmniPBR are modules used
  here.
- **No-new-claim invariant** — the rule that every number shown in a scene must come from an engine result.
- **NVIDIA Omniverse** — the platform that renders OpenUSD scenes with RTX.
- **OpenUSD (Universal Scene Description)** — the structured 3-D scene format the twins are written in.
- **Path tracing** — a high-fidelity rendering method that simulates light realistically (RTX).
- **Point instancer** — an efficient OpenUSD mechanism for placing many copies of a few prototype shapes
  (used for the mosaic cell field).
- **Prediction interval** — a range expressing forecast uncertainty (50% and 90% here); rendered as the
  envelope shells.
- **Prim** — an object in an OpenUSD scene (sphere, capsule, group, etc.).
- **Provenance** — the recorded origin of a result (model, tier, latency), carried into the scene.
- **RunPod** — a cloud GPU rental service used for the RTX rendering step.
- **RTX** — NVIDIA's ray-tracing technology for realistic lighting, glass, and glow.
- **SceneSpec** — the engine-neutral intermediate description from which scenes are authored.
- **Threshold membrane** — the translucent surface marking a clinical decision boundary (e.g., 1.8 cm SEGA,
  4 cm AML).
- **Time-sample / time code** — OpenUSD's mechanism for animation: attribute values at points along a
  timeline.
- **Uncertainty envelope** — the two translucent shells around a forecast lesion whose radii equal the
  engine's prediction intervals.
- **USDZ** — a packaged OpenUSD file for augmented reality (e.g., on a tablet).
- **VAF (variant allele fraction)** — the fraction of cells carrying a variant; rendered as the glowing
  fraction of the mosaic cell field.

---

*SYNTHETIC demonstration data throughout. Anatomy is stylized, atlas-anchored reference geometry — not
patient imaging. The TSC Intelligence Engine is investigational clinician decision-support, requires
clinician review, and is not FDA-cleared. Apache 2.0. Engine 7 of the HCLS AI Factory.*
