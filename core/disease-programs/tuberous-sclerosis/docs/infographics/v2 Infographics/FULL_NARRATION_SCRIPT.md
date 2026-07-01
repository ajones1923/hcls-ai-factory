# TSC Intelligence Engine — Full Narration Script (v2)
### Continuous voiceover for all 18 slides · for ElevenLabs voice clone

**How to use this file.** This is the complete, in-order narration for the 18-slide deck — first
person, in your voice. Each block is one slide; advance the slide at the `▸ SLIDE NN` marker. Total
runtime ≈ **18–19 minutes** at an unhurried pace.

**Reading notes for ElevenLabs:**
- `[beat]` = a short pause (~0.7s). Leave it in if your voice model honors bracketed cues; otherwise
  delete all `[beat]` markers and rely on the sentence punctuation for pacing.
- *Italicized* words are light emphasis — say them with a touch more weight.
- Numbers are written as words where they should be spoken naturally (e.g. "a hundred percent"),
  so the model pronounces them the way you'd say them aloud.
- Record slide-by-slide (one take per block) for clean editing, **or** read straight through for a
  continuous master track and cut at the `▸ SLIDE` markers.
- No institution or person is named anywhere in this script — it is safe for any audience.

---

## ▸ SLIDE 01 — The one-line story  *(~70s)*

This is the TSC Intelligence Engine — a complete precision-medicine system for a rare genetic disease
called Tuberous Sclerosis Complex, and it runs on a single small computer: an NVIDIA DGX Spark. [beat]

Here's the whole story on one screen. The problem: TSC care is fragmented across many specialists, the
most important decisions are about *timing*, and about fifteen percent of patients are told "no
mutation identified" on a blood test — a diagnosis the test simply missed, which can leave a family
searching for years.

What we built is five coordinated AI agents and a deterministic orchestrator, feeding clean clinician
screens and a three-D digital twin. It recovers those missed diagnoses from tissue, forecasts how
growths will change with honest uncertainty, surfaces hidden neuro-behavioral signals, and pulls every
option into one place — and every output is traceable.

And we measured it. On data where we know the right answers: a hundred percent classification accuracy,
a twelve-point jump in diagnoses found, six out of six of the hardest cases recovered, and everything
traceable to its source. [beat] Those are real numbers — and I'll be just as clear about what they
don't yet prove.

In the next few minutes, I'll walk you through all of it — the disease, the engine, the visuals, the
proof, and where it goes next.

---

## ▸ SLIDE 02 — The biology  *(~70s)*

Before the engine, sixty seconds on the biology — because the whole design follows from it. [beat]

Every cell has a control system that decides when to grow. Picture it like a car: there's a gas pedal
that says "grow," and a brake that says "stop." In TSC, two genes — TSC1 and TSC2 — are the brake.
When one of them is broken, the brake fails. The growth signal, called mTOR, gets stuck *on*. And the
result is benign growths in all the wrong places — in the brain, the kidneys, the heart, the skin.

Two things follow from this. First, it explains why TSC touches so many organs at once — it's one
broken brake, showing up everywhere. That's why care is so fragmented. [beat]

Second, here's the subtle part — mosaicism. The break doesn't have to be in *every* cell. Think of the
body's DNA as a print run of a book. Sometimes the misprint is in only *some* of the copies. The change
happened after conception, so only a fraction of cells carry it. And that is exactly why a standard
blood test can come back clean even when the child clearly has the disease. [beat] Hold onto that idea —
it's the key to the engine's flagship capability.

---

## ▸ SLIDE 03 — The problem we set out to solve  *(~65s)*

So let's name the problem plainly, from the family's side of it. [beat]

A child with TSC doesn't have one doctor — they have a crowd of them. A neurologist, a nephrologist, a
cardiologist, a geneticist, a developmental specialist — each holding one piece of the picture, often
not talking to each other. The parents become the messengers, carrying information from one office to
the next.

And the hardest decisions are about *timing*. Is this brain growth stable, or is it the one that's
about to need surgery? Is the kidney lesion safe to watch, or not? These calls turn on subtle trends
over time — exactly the thing that's easy to miss when the data is scattered.

Meanwhile the behavioral and psychiatric side of TSC — which families often say is the hardest part —
goes under-recognized in busy visits. And underneath it all, that fifteen percent who were told "no
mutation found," still searching for an answer. [beat]

Fragmented care, decisions about timing, hidden signals, missed diagnoses. That's the problem. The
engine is built to take on all four.

---

## ▸ SLIDE 04 — Meet the engine  *(~55s)*

So here's the answer we built. [beat]

The TSC Intelligence Engine is five specialist AI agents working as a coordinated team — each one
taking on a piece of that problem. One maps the child's full clinical picture. One handles the
genetics, including the missed-diagnosis recovery. One forecasts how growths will change over time.
One watches for the hidden neuro-behavioral signals. And one assembles every treatment and trial
option into a single sourced brief.

Tying them together is an orchestrator — not an AI, but a precise, predictable conductor that runs the
agents in the right order and keeps a permanent record of everything it does.

And it all feeds three things a clinician actually uses: a briefing before the visit, a dashboard
during it, and a three-D digital twin. [beat] Five agents, one conductor, on one small machine. Let me
open it up.

---

## ▸ SLIDE 05 — The full architecture  *(~70s)*

Now the whole machine, in one view. [beat] On the left, the inputs — a synthetic fifty-patient cohort:
the genomic data, the longitudinal records, the clinical notes, the imaging reports.

That flows into the five agents. The Phenome Mapper runs *first*, because everyone builds on it — it
organizes the symptoms and checks every code against the official medical dictionary, nearly twenty
thousand terms. The Variant Curator recovers the mosaic variant and classifies it by the worldwide
rulebook. The Trajectory Modeler — which is classical statistics, not a chatbot — forecasts the
growths. The TAND agent reads the behavioral signals. And the Therapeutics Strategist pulls it all into
a sourced options brief.

Holding them together is the orchestrator. Think of it as an air-traffic controller — it runs the
agents in the right order and keeps a permanent, tamper-proof log of everything, so the whole process
can be replayed and checked.

And around all of it is the trust layer — provenance on every output — and a validation scorecard: a
hundred percent accuracy, a twelve-point detection lift, six of six mosaics, fully traceable.

---

## ▸ SLIDE 06 — How the data flows  *(~55s)*

Let's follow the data through it, left to right. [beat]

It starts with the synthetic cohort on the left. The first stop is always the Phenome Mapper — it's the
foundation, because everything downstream reads the phenotype it builds. From there the work fans out,
in parallel, to three agents: the Variant Curator on the genetics, the Trajectory Modeler on the
forecasts, and the TAND agent on the behavioral signals. Those three then converge into the
Therapeutics Strategist, which assembles the final options brief.

Underneath, the orchestrator coordinates all of it and routes the results to the clinician's surfaces —
the briefing, the dashboard, the disciplined alerts, and the digital twin. And running along the
bottom: provenance on every output, and the measured scorecard. [beat] It's a clean pipeline — data
in, a clear picture out.

---

## ▸ SLIDE 07 — The right model for each job  *(~60s)*

One design choice I'm especially proud of: this engine doesn't just throw the biggest AI at everything.
It picks the right tool for each job — which makes it both smarter and far cheaper to run. [beat]

For the heavy reasoning — weighing treatment options, judging a tricky variant — it uses the most
capable model. For mid-level work, a faster, lighter model. For simple, high-volume tasks, the
smallest and quickest. It's like a hospital: you don't put the chief of surgery on every routine blood
draw. You match the expertise to the task.

Two things make this safe. First, the forecasting and the variant classification don't use a language
model *at all* — they're classical statistics and a fixed rulebook, so they're predictable and
repeatable. And second, if a model is ever unavailable, the engine falls back to a local one running
right on the box, and tells you it did. [beat] Right-sized intelligence, with a safety net — that's how
it stays both affordable and trustworthy.

---

## ▸ SLIDE 08 — A patient's journey through the engine  *(~60s)*

Let me make this concrete by walking one child through the whole engine. [beat]

She arrives with a clinical diagnosis of TSC but — that familiar story — a blood test that found
nothing. The Phenome Mapper goes first, building her complete profile and flagging which routine
surveillance she's overdue for. Then the Variant Curator looks at *tissue* instead of blood, and
recovers the variant the blood test missed — a confirmed diagnosis, at last.

The Trajectory Modeler takes her history of scans and projects each growth forward, with honest
uncertainty bands. The TAND agent reads through her notes and surfaces a behavioral pattern that hadn't
been formally raised. And the Therapeutics Strategist gathers the relevant trials and treatment options
into one sourced brief.

By the end, the scattered crowd of specialists' worth of information is one coherent picture — genetics,
forecasts, behavior, and options — ready before the family even sits down. [beat] That's the journey,
start to finish.

---

## ▸ SLIDE 09 — The flagship: ending the diagnostic odyssey  *(~75s)*

This is the capability I'm most proud of — and it's the heart of the whole engine. [beat]

A child clearly has TSC. But the blood test comes back: "no mutation identified." Nothing found. For
the family, that can mean *years* of uncertainty — what doctors call a diagnostic odyssey.

So why does the test miss it? Mosaicism. Remember the misprint in only some copies of a book? The
genetic change is in only some of the child's cells, and in the blood those cells are too rare to see.
But — and this is the key — the variant is usually *concentrated in affected tissue*. In a piece of a
brain tumor already removed during surgery. It's been sitting in the freezer the whole time.

So we point the Variant Curator at the tissue instead. It runs mosaic-aware analysis, finds the
low-level variant — here, a TSC2 change at about eight percent, roughly one cell in twelve — and
classifies it Likely Pathogenic by the rulebook. And it has to *reject* look-alike junk to get there,
so the zero-false-alarms number is earned, not free.

The result: detection goes from eighty-six to ninety-eight percent. Six of six of the hardest cases,
recovered. [beat] A confirmed diagnosis — at last. And that's the capability a tissue biobank unlocks.

---

## ▸ SLIDE 10 — The visuals: four ways to see a child  *(~65s)*

Now the part that makes people lean in — the engine doesn't just produce reports, it produces a living
three-D digital twin. [beat]

There are four kinds of scene. The first is a single lesion over time — a brain or kidney growth you
can watch grow along its forecast. The second is the whole-child view — every organ TSC touches, lit up
in one body, so you grasp the full burden at a glance. The third is the uncertainty itself, rendered as
geometry — I'll come back to that, because it's the soul of the project. And the fourth follows a single
case through the timeline.

All of this is built in OpenUSD — the same 3-D format the visual-effects industry uses — and rendered
with RTX in NVIDIA Omniverse, with real glass and glowing materials. [beat] It turns a page of numbers
into something a family, and a care team, can actually *see*. And every bit of it stays tied to the
real data underneath.

---

## ▸ SLIDE 11 — How the visuals are built  *(~55s)*

You might assume rendering film-quality 3-D needs a data center. It doesn't — and the trick is a clean
split. [beat]

*Building* the scene is just turning the engine's results into a structured 3-D format called OpenUSD.
That's pure data work — it runs on the little DGX Spark, on the CPU, for free, no graphics card needed.
The lesion's size and color become an animation over time; the uncertainty envelope's radii are set
*equal* to the prediction intervals; the threshold and the provenance ride along as metadata.

Then the heavy, beautiful part — the RTX path tracing, the glass, the glow — happens separately, on a
rented cloud GPU through RunPod, in NVIDIA Omniverse. In short: the Spark does the thinking and
building; RunPod does the beautiful picture.

And because the envelope's size is literally the prediction interval, the render is exactly as confident
as the engine. [beat] It can't lie.

---

## ▸ SLIDE 12 — The honest-geometry hero  *(~55s)*

I want to sit on this one image for a moment, because it captures the whole philosophy of the project.
[beat]

This is the digital twin of a child's brain tumor. The glowing core is the lesion. The translucent
glass shells around it are the uncertainty — and their size *is* the engine's fifty- and ninety-percent
prediction interval. The red ring is the threshold where the care team would discuss acting.

Think of the weather: a hurricane forecast doesn't claim one exact landfall — it draws a cone that
widens the further out you look, because the future is less certain. This is that cone, in three
dimensions, around a child's tumor. When the engine is confident, the glass hugs the core. When it's
unsure, the cloud grows. [beat]

The whole point is right there in the geometry: the prettier this render gets, the more honest it stays.
The glass cloud *is* the forecast's uncertainty. A render that cannot lie about how sure it is.

---

## ▸ SLIDE 13 — How we know it works  *(~75s)*

Now the part that matters most to me — because anyone can build something that *looks* impressive. The
hard, rarer thing is to *measure* whether it's actually *right*. [beat]

Here's how we did that. Because the demonstration runs on synthetic patients, we know the correct
answers in advance — so we can grade the engine against them, like a practice exam with an answer key.
And the report card is strong: a hundred percent classification accuracy, a twelve-point lift in
diagnoses found, six out of six of the hardest mosaics recovered, zero disease-causing variants wrongly
called harmless, and a hundred percent of outputs traceable to their source.

But here is the single most important sentence in this whole story, and I'll always say it out loud.
[beat] These numbers prove the engine's *logic* works on practice data with known answers — that's
called construct validity. They do *not* prove it works on real patients in a real clinic — that's
clinical validity, and it still has to be done.

Think of it as a flight simulator: it can rigorously prove a pilot's skills and a plane's systems — but
it is not the first real flight. This engine is a very good simulator. [beat] The first real flight is
the next step.

---

## ▸ SLIDE 14 — Why you can trust it  *(~60s)*

In medicine, a confident-sounding AI that you can't check is worse than no AI at all. So trust was a
*design* requirement here, not a feature we bolted on. Five disciplines run through the whole engine.
[beat]

One — provenance. Every single output carries a paper trail back to its source: which record, which
guideline, which calculation. Nothing is "because the AI said so."

Two — the right tool for each job. The forecasting is classical statistics, not a chatbot. The
classification follows a deterministic rulebook. Language models are used only where judgment over text
actually helps — and even then they cite their sources.

Three — honest uncertainty. Every forecast carries its prediction interval, like that hurricane cone.
Four — disciplined alerts: no more than three a week per clinician, so the signal never drowns in noise.
And five — it degrades safely: if a model is unavailable, it falls back to a local one and tells you,
rather than guessing.

Provenance, the right tool, honest uncertainty, restraint, safe failure. [beat] That's what makes it
something a clinician can actually trust.

---

## ▸ SLIDE 15 — The live demo, in three acts  *(~65s)*

Let me show you how this actually plays out, the way I run it live. It's three acts. [beat]

**Act one — the briefing.** Before the visit, the engine has already done the prep. It opens on the
hardest case: the child whose blood test said "no mutation identified." You watch it point at the
tissue specimen instead, recover the variant, and turn a non-diagnosis into a confirmed one. That's the
gut-punch moment.

**Act two — the dashboard.** We move to a child being watched for a slow-growing brain tumor. One
screen: the genetics, the phenotype, and the forecast — the growth trending toward the line where the
team would discuss acting, with the honest uncertainty band around it. The whole picture, at a glance.

**Act three — the twin.** Then it lifts off the page into three dimensions. The same tumor, in NVIDIA
Omniverse, growing along its forecast inside the glass uncertainty cloud, crossing the threshold — the
render that's exactly as confident as the engine.

Briefing, dashboard, twin. [beat] Three acts, about ten minutes, one small computer.

---

## ▸ SLIDE 16 — Where it goes next  *(~65s)*

So where does this go from here? The path is deliberately concrete. [beat]

**Phase one is built** — what you've seen today: the full five-agent engine and the 3-D twin, running on
one small computer, validated against synthetic data with known answers.

**Phase two is the first real flight** — and it turns on one thing: tissue. A partner children's
hospital has something extraordinary sitting in its freezers — a biobank of specimens from past
surgeries, with years of follow-up already attached. Remember the flagship capability: recovering a
missed diagnosis from tissue. So the bridge study is almost poetic — take patients who were once told
"no mutation identified," point the engine at the specimen already on the shelf, and see how many
diagnoses we can recover *retrospectively*. The answers are, in a sense, already known — they're just
waiting to be read.

That study is the move from construct validity to clinical validity — from the simulator to the first
real flight. And it connects naturally to the people and infrastructure of a research institution: the
biobank, the informatics group, the clinical TSC program. [beat] One disease, proven end to end — then
it travels.

---

## ▸ SLIDE 17 — It's a blueprint, not a one-off  *(~55s)*

Here's the part that makes this bigger than one rare disease. [beat]

TSC was the proving ground — but almost nothing in the engine is *specific* to TSC. The architecture is
the real product: five coordinated agents, a deterministic orchestrator, honest uncertainty, provenance
on everything, an anatomical twin. That pattern is disease-agnostic.

The way I think about it: we built the wiring, and TSC is just the labels on the boxes. To move to
another condition, you swap the labels — the genes, the guidelines, the growth patterns — and keep the
wiring. The same engine could serve neurofibromatosis, Rett syndrome, other mTOR-related disorders, or
any rare disease where care is fragmented, timing is everything, and mosaicism hides diagnoses.

So this isn't a one-off tool for one disease. [beat] It's a blueprint for precision medicine in rare
disease — proven once, then reused.

---

## ▸ SLIDE 18 — Small, affordable, open  *(~60s)*

Let me close with what I think is the quietly radical part of this. [beat]

Everything you've seen — five AI agents, a deterministic orchestrator, classical forecasting, a
film-quality 3-D twin — all of it runs on a single desktop computer that costs about forty-seven
hundred dollars. Not a data center. Not a seven-figure cluster. A box that sits under a desk, with the
beautiful rendering rented by the hour from the cloud only when you need it.

And it's open source — Apache 2.0. Anyone can read it, run it, check it, and build on it. That matters
in medicine, because trust comes from being able to look inside.

So the whole thing comes down to three words. *Small* — it fits on one machine. *Affordable* — within
reach of a single lab or clinic, not just the largest institutions. And *open* — transparent and
reusable by anyone. [beat]

Precision medicine for a rare disease, recovering diagnoses that were missed, on a computer that fits on
a desk. That's the TSC Intelligence Engine. Thank you.

---

### Runtime summary

| Slides | Section | Approx. time |
|---|---|---|
| 01–03 | The story & the disease | ~3.5 min |
| 04–07 | The engine & architecture | ~4 min |
| 08–09 | Patient journey & the flagship | ~2.5 min |
| 10–12 | The Omniverse visuals | ~3 min |
| 13–14 | Proof & trust | ~2.5 min |
| 15–18 | Demo, roadmap, generalization, close | ~4 min |
| **Total** | **18 slides** | **≈ 18–19 min** |

*All narration is first-person, in your voice. No institution or individual is named — safe for any
audience. Pair each block with its matching Nano Banana Pro frame (`NN_*.md`).*
