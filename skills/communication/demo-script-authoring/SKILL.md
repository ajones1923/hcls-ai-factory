---
name: demo-script-authoring
description: >-
  How to write a demonstration walkthrough or script for the HCLS AI Factory. Use whenever
  scripting a live demo, a recorded walkthrough, a talk narrative, or a stage arc for any
  audience — clinicians, builders, protein scientists, the public, or families. Encodes the
  three-beat arc (weight → compression → hope), the portfolio-not-one-demo rule, the doorway
  principle, inline honesty flags, and the two-cut split, so a script lands emotionally AND
  survives expert scrutiny.
---

# Demo-Script Authoring — write the walkthrough that earns trust

A demo script is a **staged patient story**, not a feature tour. It opens on a human at stake,
compresses a years-long odyssey into a single governed afternoon inside the factory, and lands on
a multiplier that outlives the room. This skill is how you write one that thrills the public and
leaves the specialist nodding.

## The three-beat arc (build every script to this)

Every disease-program demo is one patient, the whole factory, three beats. Write to the arc, not
to a checklist of tools.

1. **WEIGHT — open on human stakes.** Start with a frightened family, a years-long diagnostic
   odyssey, an infant in the NICU — a concrete person, not a capability. One quiet, true sentence
   about a child hits harder than any montage and immunizes the piece against the "using sick kids
   for marketing" read. Establish what is at risk and how long the current path takes.
2. **COMPRESSION — the factory collapses the odyssey into a governed afternoon.** Walk the engines,
   agents, and tools **in sequence**, naming each as the story reaches it: genome interpreted →
   variants reasoned over → images and molecules cross-read → targets, structures, trial matches.
   Show the odyssey of years becoming an afternoon — governed, reproducible, and honest at every
   step. This is the body of the script.
3. **HOPE — land the multiplier.** Close on what makes this bigger than one patient: it is
   open-source (Apache-2.0, fork it and run it), it runs on one affordable box, and *the next
   patient gets the same afternoon* — this child's data makes the next child's afternoon faster via
   the biobank. Hope must be true hope the factory can actually deliver, never a cure implied.

> The unit of impact is the **disease program** — one patient, the whole factory — not an isolated
> capability. If a beat doesn't serve weight, compression, or hope, cut it.

## Portfolio, not one demo (the realism guardrail)

Do not force a single demo to touch every engine, agent, and model. That breaks clinical realism
and violates the honesty gate — a demo that "does everything" is a demo that overclaims. Coverage
is a property of the **portfolio (D1–D7)**, achieved across demos:

- Pick the **home demo** whose patient story the capability actually belongs to.
- Let flagship demos (D3 whole-factory, D2 cross-modality) carry breadth; let completion demos
  (D4/D5/D6/D7) go deep for a specialist audience.
- If a capability has no honest home in the story you're telling, it belongs in a different demo —
  not shoehorned into this one to look comprehensive.

## The doorway principle

Each spotlight is a **door into the whole factory**, not a dead end. When the script goes deep on
one thing — a fusion call, a docked molecule, a CAC score — narrate how that spotlight *recruits*
the rest of the platform: the registry that makes it composable, the agent that reasons over it,
the trial match it feeds. A focused demo should leave the room understanding the whole, reached
through one door.

## Honesty flags go INLINE (say them out loud, in the beat)

Honesty is load-bearing and disarming — experts convert from skeptic to advocate on it alone. Put
the limit **in the narration where it applies**, not in a disclaimer slide at the end:

- **"Decision support for a qualified clinician — not diagnosis or prescribing."** Say it the first
  time any clinical output appears.
- **Preclinical** (e.g. gene therapy for TSC1/TSC2, de novo binder design): "this is the open
  design bench, not a cure today."
- **Research-use / QA** (e.g. synthetic imaging augmentation): "for research and QA, never a
  diagnostic source."
- **Roadmap** (e.g. single-cell atlas similarity, cell-embedding models): name it as not-yet-live.
- **Elastic burst, not "all on one box."** Heavy or ARM-incompatible models (Chai-1/2, RFdiffusion,
  Evo 2) burst to remote GPUs over a private mesh. Say **"elastic burst"** — never imply everything
  runs locally when it doesn't.

## LIVE vs REPRESENTATIVE vs BURST (label what the room is seeing)

The audience must always know what kind of thing is on screen. A `live` capability is never
mock-served; if it isn't live, say what it is:

- **LIVE** — actually running now, in front of them, on real input. The strongest possible move:
  "you can run this yourself."
- **REPRESENTATIVE** — a pre-computed or curated result standing in for a long-running step (a
  genome alignment that takes hours, a cohort you can't process live). Call it representative; never
  pass it off as live.
- **BURST** — runs live, but on remote GPUs via elastic burst rather than the local box. Name the
  burst so "one box" stays honest.

## The two-cut principle

The beat that thrills the public is often the one that makes a specialist wince, and vice versa.
When a script must serve both, write **two cuts**, not one compromised script:

- **Technical cut** — for experts and builders: tools running live, reproducibility, honest limits,
  provenance. Density over drama.
- **Mission cut** — for the public and families: human stakes, understated, one legible "wow" plus
  the "and it's yours" hook.

Don't hand one script two contradictory jobs. When drama and accuracy conflict, **accuracy wins.**

## Suggested room arc (a full-session script)

A reference sequence that builds weight → proof → breadth → depth → mission:

1. **Open on D3** (flagship TSC) — the whole factory converging on one infant; run the full
   weight → compression → hope arc so the room feels the stakes first.
2. **Prove with D1** (genome → novel molecule) — "raw reads to a validated lead," so the engine is
   shown to be real, not narrated.
3. **Cross-modality with D2** — imaging + genomics + pharmacogenomics together; the platform
   reasoning across modalities, not one pipeline.
4. **Go deep for the audience** — D4 for oncologists (molecular tumor board), D7 for neurologists
   (Parkinson's/ALS/MS suite), D5 for protein scientists (CAR-T / biologic design). Pick the
   completion demo that matches the room.
5. **Close on the mission** — hand a discovered molecule to the trial agent, return to the child,
   and land the multiplier: the next patient gets the same afternoon.

## Do / Don't

**Do:** open on a real person; walk engines/agents/tools in story sequence; label LIVE vs
representative vs burst; put honesty flags inline in the beat they apply to; keep each demo in its
honest home; open a doorway to the whole factory; split into a technical cut and a mission cut when
audiences conflict; end on the open-source, one-box, next-patient multiplier.

**Don't:** script a feature tour with no patient; force one demo to touch everything; show a
representative result as live or imply "all on one box" during a burst; present clinical output as
diagnosis; bury limits in a closing disclaimer; let drama outrun the honesty ledger; overpromise a
cure in the hope beat.

## Authoring checklist

1. Name the **patient and the odyssey** (weight) before any technology.
2. Choose the **home demo** (D1–D7); confirm you are not overloading one demo.
3. Sequence the **compression** beat: which engines/agents/tools, in what order, telling the story.
4. Mark every step **LIVE / REPRESENTATIVE / BURST**; insert **inline honesty flags** where they
   apply (decision-support, preclinical, research-use, roadmap, elastic-burst).
5. Write the **hope** beat as a true multiplier (open-source, one box, next patient via biobank).
6. Decide if the room needs a **technical cut, a mission cut, or both**.
7. Read it back as a **skeptical expert**: is any claim ahead of the honesty ledger? Cut or qualify.

## Related
- `demo-foundation-alignment` — the D1–D7 portfolio, coverage matrix, and honesty ledger the script draws from.
- `broad-general-persona` — the five audiences, the two-cut principle, and "accuracy wins over drama."
- `hcls-core-vision-mission` — the mission the story dramatizes; why honesty is load-bearing.
