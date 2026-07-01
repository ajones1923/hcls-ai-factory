# The TSC AI Engine — A Foundational Learning Guide

### What it is, the problems it tackles, and how it works — explained from the ground up

**Author:** Adam M. Jones — HCLS AI Factory, TSC Intelligence Engine (Engine 7)
**Reading level:** written for a general reader (about 8th–9th grade). Every technical word is explained the first time it appears, and there is a glossary at the end.
**License:** Apache 2.0 (open source) · 2026

> **How to read this guide.** Start at the beginning and go in order. Part 1 explains the
> disease. Part 2 explains why caring for it is so hard today. Part 3 onward explains the AI
> engine that was built to help — piece by piece — and what it does and does not do. You do not
> need a medical or computer background. **One important note up front:** everything in this
> demonstration uses *made-up* patient data (we call it "synthetic"). No real patients are
> involved yet. The engine is a decision-support tool for doctors — it never replaces a doctor,
> and it is not approved by the FDA. We will repeat that, because honesty is the whole point.

---

## Part 1 — First, what is TSC?

**TSC** stands for **Tuberous Sclerosis Complex**. It is a rare genetic disease. "Genetic" means
it is caused by a change in a person's genes. Think of your genes as the instruction manual that
tells your body how to build and run itself. TSC happens when there is a "spelling error" in one
of two genes, named **TSC1** or **TSC2**. That spelling error lets certain cells grow when they
should not.

Because those two genes matter all over the body, TSC is a **multisystem** disease — it affects
many organs at once. People with TSC can have:

- **Brain:** growths called **tubers**, and a special tumor called a **SEGA** (we will keep saying
  "SEGA" — just know it is a brain tumor that can block the normal flow of fluid in the brain if it
  gets big). Many people with TSC also have **seizures**.
- **Kidneys:** growths called **angiomyolipomas**, usually shortened to **AML**. If an AML gets
  large (around 4 centimeters), it can bleed, which is dangerous.
- **Skin:** light-colored patches and small bumps.
- **Heart and lungs:** other kinds of growths.
- **Brain and behavior:** a group of learning, mood, and behavior challenges that doctors call
  **TAND** (TSC-Associated Neuropsychiatric Disorders). This is extremely common but often missed.

Two more facts that matter a lot for this guide:

1. **TSC is lifelong and needs constant watching.** Doctors follow a checklist of regular scans and
   check-ups (this is called **surveillance**) to catch problems early.
2. **TSC is rare.** Far fewer people have it than have, say, diabetes. That is part of why it gets
   less attention from companies building medical tools — there is less money in rare diseases. That
   gap is exactly where this project aims.

---

## Part 2 — Why caring for TSC is so hard today (the struggles)

If TSC care were easy, we would not need a new tool. Here are the real problems families and
doctors run into right now. The AI engine is built to chip away at each one.

### Struggle 1 — The "diagnostic odyssey": some patients never get a clear answer

To confirm TSC, doctors usually run a **genetic test** on a blood sample, looking for that spelling
error in TSC1 or TSC2. But here is the catch: in about **15 out of every 100 patients**, the blood
test comes back saying **"no mutation identified"** — meaning *it found nothing*, even though the
person clearly has the disease.

Why does that happen? Because of something called **mosaicism**. Imagine a print run of a book where
only *some* of the copies have a misprint. In a mosaic patient, only *some* of their cells carry the
spelling error. If those cells are rare in the blood, the blood test misses them. The fraction of
cells carrying the change is called the **VAF** (variant allele fraction). A low VAF — say 8% — is
easy to miss in blood.

Families in this situation can spend years without a clear genetic answer. Doctors call this long,
frustrating search the **diagnostic odyssey**. The key insight: the answer is often hiding in
*tissue that was already removed during surgery* (like a piece of a brain tuber taken out to control
seizures). The cells with the change are concentrated there — if someone just looks.

### Struggle 2 — Care is split across many specialists

A person with TSC might see a brain doctor (neurologist), a kidney doctor (nephrologist), a genetics
doctor, a skin doctor (dermatologist), a brain surgeon (neurosurgeon), and a mental-health
specialist. Each one sees a slice of the patient. Nobody automatically holds the *whole picture*.
Important details fall through the cracks between specialists.

### Struggle 3 — Predicting the future is genuinely hard

The biggest decisions in TSC are about **timing**:

- *When* will this child's SEGA grow large enough that the team should talk about surgery?
- *When* will this kidney AML get close to the size where it might bleed?
- How *confident* are we in that prediction?

Doctors answer these by looking at a few past scans and building a mental guess. That is skilled
work, but it is invisible (nobody else can see the doctor's mental model) and it does not come with a
clear measure of *how sure* anyone is.

### Struggle 4 — The behavior and learning side (TAND) is under-recognized

About **9 out of 10** people with TSC have some TAND challenge — trouble with attention, anxiety,
learning, or behavior. But fewer than **2 out of 10** actually get formally checked for it. The
signals are often buried in clinic notes as soft, easy-to-skip phrases like *"mom mentions some
trouble focusing at school; will revisit next time."* Those quiet notes rarely get pulled together
into a clear pattern.

### Struggle 5 — Too much information, and alert fatigue

Medical software loves to fire alerts. When there are too many, doctors start ignoring all of them —
even the important ones. This is called **alert fatigue**. Any helpful tool has to be *disciplined*:
it should speak up only when it really matters.

### Struggle 6 — The data is scattered and hard to trust

Patient information lives in many disconnected systems — scans here, lab results there, free-text
notes somewhere else. And when a computer suggests something, doctors reasonably ask: *"Why should I
believe this? Where did it come from?"* If a tool cannot show its work, it will not be trusted.

---

## Part 3 — Meet the TSC AI Engine (the big picture)

The **TSC Intelligence Engine** is software built to help doctors with all six struggles above. It
is **"Engine 7"** in a larger open-source platform called the **HCLS AI Factory** (HCLS means
Healthcare and Life Sciences).

Here is the simplest way to picture it: **a small team of five expert assistants, plus a conductor
who keeps them organized, who together turn a patient's messy data into a few clear, trustworthy
screens for the doctor.**

The five expert assistants are called **agents**. An "agent" here is a focused piece of AI software
that does one job well. The conductor is called the **orchestrator**. And the screens the doctor
actually looks at are called **surfaces**.

A few things make this engine unusual:

- **It runs on one small, affordable computer.** The whole thing runs on a single **NVIDIA DGX
  Spark** — a desktop-sized AI computer that costs about **$4,699**. You do not need a giant data
  center.
- **It is open source (Apache 2.0).** The code is shared openly, not locked away.
- **It is honest by design.** Every answer it gives can be traced back to its source, it shows *how
  sure* it is, and it clearly labels everything as a demonstration on made-up data.

The rest of this guide walks through each part.

---

## Part 4 — The five agents, one by one

Each agent solves one of the struggles. They run in a set order, because some need the work of
others first. Below, each agent gets: *the job, the problem it solves,* and *how it works in plain
terms.*

### Agent 1 — The TSC-Phenome Mapper (runs first)

**The job.** Build a clean, organized list of everything going on with the patient's body over time —
their **phenotype** (the set of observable traits and symptoms). It runs *first* because every other
agent builds on its work.

**The problem it solves.** Patient information is scattered and written in everyday language. This
agent gathers it and tags each finding with a standard medical code so computers (and other doctors)
all mean the same thing.

**How it works, plainly.** Doctors around the world use a shared dictionary of symptom codes called
the **HPO** (Human Phenotype Ontology) — think of it as a giant, official list where "renal
angiomyolipoma" always has the exact same ID number. This agent reads the patient's records and
notes, finds the symptoms, and matches each one to its official HPO code. Importantly, it
**double-checks every code against the real, official HPO list (19,389 terms)**. If the AI guesses a
code that is not real, the agent throws it out instead of trusting it. (In testing, this safety check
even caught two miscoded entries in the practice data itself.) The agent also flags **surveillance
gaps** — scans or check-ups that are overdue.

### Agent 2 — The TSC-Variant Curator (the headline act)

**The job.** Look at the patient's genetic data, figure out which genetic change they have, and
decide how serious it is.

**The problem it solves.** This is the agent that tackles the **diagnostic odyssey** from Struggle 1.
It is built to *recover the low-VAF mosaic changes that standard blood tests miss* — by looking in
the right tissue.

**How it works, plainly.** Genetic data comes as a long list of differences between the patient's DNA
and a reference. This agent reads that list, finds the real change, and then classifies it using the
worldwide rulebook for genetic findings, called **ACMG-AMP**. ACMG gives a verdict like *"Pathogenic"*
(disease-causing), *"Likely Pathogenic,"* or *"Uncertain."*

Two details make this agent trustworthy:

1. **The verdict is decided by a fixed rulebook, not by the AI's opinion.** A reliable, by-the-book
   calculator makes the final call. The AI's job is to *explain* the reasoning in plain language, not
   to overrule the rules. (When connected to a real AI model, it produced a genuine, expert-level
   written explanation — while the rulebook stayed in charge.)
2. **It has to reject fakes.** The genetic data is built so that every sample also contains a
   look-alike "junk" signal (a technical artifact). A lazy tool might report the junk as real. This
   agent has to *throw it out*. That is why we can say the engine has **zero false alarms** on
   disease-causing calls — it *earned* that, it did not get it for free.

### Agent 3 — The TSC-Trajectory Modeler (forecasting with honesty)

**The job.** Predict how the patient's growths and other measures will change over the next 6, 12,
and 18 months — and say *how confident* it is.

**The problem it solves.** This is Struggle 3: the timing question. *When* will a SEGA or AML reach the
size where the team needs to act?

**How it works, plainly.** This agent is **not** an AI chatbot — it is **classical statistics** (math
and probability), which is the right tool for forecasting numbers. It takes a few past measurements
and projects them forward, and — this is the important part — it draws a **range of uncertainty**
around the prediction, not just a single line.

Here is the best way to picture it: think of a **hurricane forecast cone** on the news. The weather
service does not say "the storm will hit this exact spot." It draws a *cone* that gets wider the
farther out you look, because the future is less certain. This agent does the same for a lesion: it
shows the most likely size *and* a 50% and 90% "cone" around it. When the data is clear, the cone is
thin. When the data is sparse, the cone is wide and honest about it.

It does this for **four things at once**: brain SEGA size, kidney AML size, **kidney function
(eGFR)**, and **seizure frequency**. For each, it estimates *when* it might cross an important
threshold, grades that risk as **likely / possible / unlikely**, and even suggests whether the
patient should be checked *more often* than the standard schedule.

### Agent 4 — The TAND Surveillance Agent

**The job.** Find the easy-to-miss learning, mood, and behavior signals (TAND) hidden in the doctor's
notes, and gently surface them.

**The problem it solves.** Struggle 4 — TAND is everywhere but rarely caught.

**How it works, plainly.** People write about worries in soft, hedging language: *"seems a bit
distracted,"* *"mom reports some trouble,"* *"will reassess next visit."* This agent is trained to
notice exactly those **discourse markers** — the little uncertainty phrases — using a method developed
by researchers at Cincinnati Children's (the **Marshall-Hagedorn** approach). It then scores how
strong and how repeated the pattern is across many notes.

Two safety rules make it trustworthy and calm:

- **It respects "no."** If a note says *"mom denies any behavior concerns,"* the agent does **not**
  count that as a problem. (This sounds obvious, but many tools get it wrong.)
- **It never sets off an alarm.** TAND findings show up as gentle *briefing material* for the doctor
  to consider — never as an interruptive alert. This avoids the alert fatigue from Struggle 5.

### Agent 5 — The TSC-Therapeutics Strategist

**The job.** Pull everything together into a clear, organized summary of treatment options — and back
up every statement with a source.

**The problem it solves.** It fights Struggle 2 (scattered care) by combining the other four agents'
work into one place, and Struggle 6 (trust) by citing where each claim comes from.

**How it works, plainly.** It reads the variant, the phenotype, the forecasts, and the TAND findings,
then adds outside information — **clinical trials** the patient might qualify for, published research
(found using a search method called **RAG**, which fetches relevant documents), and official FDA
actions. It writes a **six-section options brief**. The most important rule: **every claim must name
its source.** It is decision-support — it lays out options for the doctor, it does not give orders.

---

## Part 5 — The conductor and the screens

### The orchestrator (the conductor)

The five agents need to run in the right order and share their results. The **orchestrator** is the
conductor that makes that happen. Like the agents, it is **deterministic** — meaning it follows fixed
rules and is *not* an AI guessing. It is more like an **air traffic controller**: it runs the Phenome
Mapper first, then the agents that depend on it, then the Therapeutics Strategist, and it keeps a
permanent, tamper-proof log of everything that happened (so the whole process can be replayed and
checked later).

If an agent ever fails, the orchestrator does the safe thing: it marks that section as "needs
attention" rather than quietly showing stale or wrong information. **Honest about failure, never
silent.**

### The surfaces (what the doctor actually sees)

All that work funnels into a few clean **surfaces**:

1. **Pre-Visit Briefing** — a one-screen summary for the night before a clinic visit, with at most
   three action items. Quick to read on a phone.
2. **In-Visit Dashboard** — the main screen during the appointment, split into four quadrants:
   the genetic variant, the phenotype timeline, the forecasts, and the TAND notes. Each finding can
   be traced back to its source.
3. **Async Alerts** — a *disciplined* alert screen. It is built to never overwhelm a doctor (aiming
   for no more than about three alerts per clinician per week), and every alert carries a clear reason
   and a way to dismiss it.
4. **The Anatomical Digital Twin** — a 3D view, explained next.

---

## Part 6 — The 3D digital twin (the part that wows people)

This is where the engine gets visually striking. A **digital twin** here means a **3D model**, almost
like a scene from a video game, built *automatically from the engine's real results*. It is rendered
using **NVIDIA Omniverse**, a tool for high-quality 3D graphics.

There are four kinds of scenes:

- **The lesion-trajectory twin.** It shows the patient's brain tumor as a 3D shape that **grows along
  the forecast** as you scrub a timeline from two years ago into the future. Around it is a glowing,
  see-through **"cloud" of uncertainty** — and here is the clever, honest part: *the size of that
  cloud is exactly the forecast's uncertainty range.* When the engine is confident, the cloud is
  small and tight. When it is unsure, the cloud is big and soft. **The picture literally cannot lie
  about how sure the engine is.** As the tumor nears the danger line, it changes color.
- **The mosaic "powers-of-ten" scene.** This dramatizes the diagnostic-odyssey win. It shows a field
  of cells where **exactly the right fraction glows** — for a patient with 8% VAF, about **1 cell in
  12 lights up gold**. It turns the abstract number "8% VAF" into something you can literally *count*.
  It shows the exact signal a blood test calls "negative."
- **The whole-child atlas.** A stylized body where each affected organ (brain, kidneys, skin, and so
  on) **lights up** based on that patient's symptoms — the multisystem disease in one glance.
- **The population view.** All 50 demo patients shown together, with the **recovered mosaic patients
  ringed in gold halos** — the "we found what others missed" story, at a glance.

**One honest note about the 3D:** these are *stylized, general* body shapes scaled by the real
measurements — **not** pictures of a specific patient's actual scans. Using real scan images is a
future step.

**Why the computer setup is clever.** Building these 3D scenes is done on the small, affordable DGX
Spark, with **no expensive graphics card needed** for that step. The fancy, movie-quality
*rendering* (the lighting, glass, and glow) is done separately on a rented cloud computer from a
service called **RunPod**. In short: **the Spark does the thinking and building; RunPod does the
beautiful picture.** This keeps costs low.

---

## Part 7 — How we know it works (and how we stay honest)

Anyone can build something that *looks* impressive. The harder, rarer thing is to **measure** whether
it is actually *right*. Because the demo uses made-up patients, the team knows the "correct answers"
in advance — so the engine can be graded against them, like a practice test with an answer key.

Here is the report card (these are real, measured results on the demo data):

| What was measured | Result |
| --- | --- |
| Got the genetic verdict right | **100%** (on the demo data) |
| Extra diagnoses found vs. a standard pipeline | **+12 out of 100 patients** |
| Hard-to-find mosaic cases recovered | **6 out of 6** |
| Disease-causing variants wrongly called harmless | **0** |
| Answers that can be traced to their source | **100%** |

Those are strong numbers. But here is the **most important sentence in this whole guide**, and the
team says it out loud every time:

> These results show the engine's *logic works correctly on practice data with known answers*. That
> is called **construct validity**. It is **not** the same as proving it works on *real patients in a
> real clinic*. That bigger proof — called **clinical validation** — still has to be done.

Why brag about a limit? Because that honesty is what makes the tool trustworthy. The engine never
pretends to be more than it is. Think of it like a **flight simulator**: it can prove a pilot's
skills and a plane's systems work, but it is not the same as the first real flight. This engine is a
very good simulator. The "first real flight" is the next step.

---

## Part 8 — What it solves and who it helps (the benefits)

Let us connect every struggle from Part 2 to a benefit.

| The struggle (today) | How the TSC AI Engine helps |
| --- | --- |
| Diagnostic odyssey — blood tests miss mosaic cases | The Variant Curator recovers low-VAF mosaic variants from tissue (in the demo: **6 of 6**, a **+12-point** jump in diagnoses) |
| Care split across many specialists | The orchestrator pulls everything into one shared, organized picture and a single dashboard |
| Hard to predict timing of growths | The Trajectory Modeler forecasts SEGA, AML, kidney function, and seizures *with an honest uncertainty range* |
| TAND is under-recognized | The TAND agent surfaces the quiet, easy-to-miss signals as gentle briefing material |
| Alert fatigue | A disciplined alert design that speaks up only when it matters |
| Scattered, hard-to-trust data | Every output is traceable to its source, and the engine is honest about how sure it is |

**Who benefits, in plain terms:**

- **Patients and families** — fewer years lost to the diagnostic odyssey, earlier warning before a
  problem becomes an emergency, and learning/behavior needs caught sooner.
- **Doctors** — the whole picture in one place, with the math done and the sources cited, so they can
  spend their attention on judgment instead of hunting for information.
- **Care teams** — a shared, consistent view across specialties.

---

## Part 9 — Why "small, cheap, and open" matters

It would be easy to build something like this only for giant, rich hospitals with huge data centers.
This project deliberately did the opposite:

- **Small and affordable:** it runs on a single ~$4,699 desktop AI computer. A children's hospital, a
  research lab, or even a motivated team could afford one.
- **Open source:** the code is shared (Apache 2.0), so others can inspect it, trust it, and build on
  it. Nothing is hidden.
- **A reusable pattern, not a one-off:** the same design can be pointed at *other* rare diseases —
  swap the disease-specific pieces, keep the structure. The engineers describe it as *"swap the box
  labels, keep the wiring."* What works for TSC could help with diseases like NF1, NF2, Rett
  syndrome, and others. That is why it is called an "engine" and not just a project.

---

## Part 10 — The road ahead

### From practice data to real patients

Right now, everything runs on synthetic (made-up) data. The single most important next step is to
test the engine on **real tissue** — and there is a smart, low-risk way to start.

Many hospitals keep a **biobank**: a carefully stored, consent-approved collection of tissue samples
from past surgeries. Cincinnati Children's has one called the **Discover Together Biobank**. The
proposed first real study is a **retrospective "no-mutation-identified" re-analysis**: take patients
who were told years ago that their blood test found nothing, and run the Variant Curator on their
*already-stored* tissue. If it recovers diagnoses that were missed, that is a real, publishable
result — *using tissue that is already in the freezer.* No new surgery, no new wait.

### The two phases

- **Phase 1 — the demonstration (done):** the full engine on the Spark with synthetic data, plus the
  Omniverse 3D twin rendered on RunPod. This is the version that earns the meeting.
- **Phase 2 — the institution:** connecting to a real hospital's pieces — the **biobank** (tissue),
  the research infrastructure, the clinical-informatics team and methods, the TSC clinic (patients),
  and the hospital's medical record system (the data plumbing).

### Why it could matter

If Phase 2 works, this stops being "a clever demo" and becomes "a tool that found a real child's
diagnosis that was missed for years." That is the line between *impressive* and *important*. The
engine is the apparatus built to make that test fast, cheap, and credible.

---

## Glossary (plain-language definitions)

- **Agent** — a focused piece of AI software that does one job. This engine has five.
- **ACMG-AMP** — the worldwide rulebook for deciding how serious a genetic change is.
- **Alert fatigue** — when too many computer alerts cause people to ignore them all.
- **AML (angiomyolipoma)** — a kidney growth common in TSC; dangerous if it gets large enough to bleed.
- **Biobank** — a hospital's stored, consent-approved collection of tissue samples.
- **Clinical validation** — proving a tool works on *real patients in a real clinic*. (Not yet done.)
- **Construct validity** — proving a tool's logic works on *practice data with known answers*. (Done.)
- **Deterministic** — follows fixed rules; not an AI making guesses.
- **Diagnostic odyssey** — a long, frustrating search for a diagnosis that standard tests keep missing.
- **Digital twin** — a 3D model built automatically from the engine's real results.
- **eGFR** — a number that measures how well the kidneys are working.
- **HPO (Human Phenotype Ontology)** — the official shared dictionary of symptom codes.
- **Mosaicism** — when only *some* of a person's cells carry a genetic change (like a misprint in only
  some copies of a book).
- **Orchestrator** — the "conductor" that runs the agents in order and keeps a complete log.
- **Phenotype** — the set of a person's observable traits and symptoms.
- **Provenance** — the paper trail showing where an answer came from (which model, when, how long it took).
- **RAG** — a search method that fetches relevant documents so the AI can cite real sources.
- **SEGA** — a brain tumor in TSC that can be dangerous if it grows.
- **Surveillance** — the schedule of regular scans and check-ups that catch TSC problems early.
- **Synthetic data** — made-up, realistic practice data (no real patients).
- **TAND** — TSC-related learning, mood, and behavior challenges; very common, often missed.
- **TSC (Tuberous Sclerosis Complex)** — the rare genetic disease this engine is built for.
- **VAF (variant allele fraction)** — the fraction of a person's cells that carry a genetic change.

---

*SYNTHETIC demonstration data throughout. The TSC Intelligence Engine is decision-support for
clinicians, requires clinician review, and is not FDA-cleared. Apache 2.0. Engine 7 of the HCLS AI
Factory.*
