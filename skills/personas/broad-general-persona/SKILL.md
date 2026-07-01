---
name: broad-general-persona
description: >-
  Broad audience-awareness and messaging discipline for ANY outreach or communication artifact for
  the HCLS AI Factory — promo, README, talks, social, decks, docs, landing copy, release notes.
  Consult whenever writing or reviewing anything people outside the codebase will read or watch, so
  the artifact reaches the right audience AND leaves experts trusting what they find. Encodes the
  five audience personas, the universal guardrails, and the honesty discipline. Pairs with the
  north star (the mission it protects) and demo-foundation-alignment (what it shows).
---

# Broad / General Persona — audience-aware, honest outreach

Every outreach artifact amplifies **two reactions at once — awe and scrutiny**. The louder and
slicker it goes, the more of each it generates: the same footage or copy that thrills a general
audience is often exactly what makes a domain expert wince. Messaging is therefore never
one-size-fits-all — it is the deliberate management of that trade-off, persona by persona.

**The goal of outreach is not to fulfill the mission** (no artifact can). It is to make the right
people curious enough to reach the repository, and honest enough that experts trust what they find.
Optimizing for spectacle at the expense of honesty produces more attention and less of what
actually matters: forks, collaborators, and downstream impact.

## Two axes always in tension
- **Reach vs. Credibility** — breadth of attention against trust from the people whose respect the
  project needs.
- **Spectacle vs. Honesty** — emotional impact against accuracy about what is proven, preclinical,
  or roadmap.

The durable advantage is that the project is **real and open**. Openness is the argument that ends
most fights — a critic cannot credibly call it vaporware if they can fork it and run it themselves.

## The five personas

| Persona | Lead with | Avoid |
|---|---|---|
| **AI-curious public** | one legible "wow" + the *"and it's yours"* hook | letting them infer "cures disease now" |
| **Domain experts** (clinicians, geneticists, oncologists, pharma R&D) | precision, real tools shown running, limits stated plainly | blurring decision-support into "diagnosis"; spectacle outrunning evidence |
| **Builders / open-source** (home turf) | a forkable, documented repo; Linux-style framing; Apache-2.0 commercial-OK; clean provenance | polish over a messy or non-reproducible repo |
| **Patients / families** (highest-stakes) | dignified, understated, concrete true hope | hope the project can't yet deliver |
| **Skeptics / competitors** (guaranteed) | openness, reproducibility, "fork it and run it yourself" | any gap between claim and reality — it becomes their ammunition |

Notes that matter: **domain experts convert from skeptic to advocate on honesty alone**, and their
public pushback carries weight — losing them is expensive and hard to reverse. **Patients/families
are the most sensitive ground** — restraint there is not just ethics, it protects vulnerable
people and pre-empts the "false hope / using sick kids for marketing" critique. **Skeptics get
louder the louder the artifact is** — openness neutralizes them.

## Universal guardrails (apply to everything)

- **Show it running, not just its results.** *"You can run this yourself"* beats any claim and is
  the one thing a competitor cannot dismiss.
- **Put the limits on screen / in the text.** "Decision support, not diagnosis." "Preclinical."
  "Roadmap." Labeling limits inside outreach is rare and disarming — experts trust the person who
  does it.
- **Make openness the hook, not just a feature.** The shock is not "AI did this." It is *"AI did
  this — and it's yours."* That line is the project's alone.
- **Open-source is not non-commercial.** State plainly that Apache-2.0 welcomes commercial use; a
  business accelerating precision medicine is the goal succeeding, not leaking.
- **Handle the mission with dignity and understatement.** A quiet, true sentence about a child
  hits harder than a montage — and it immunizes against the manipulation read.
- **Never overpromise.** Overstatement betrays the mission more than any missing feature. If a
  claim can't be defended in front of a skeptical expert, cut or qualify it.

## The format failure mode (watch for this)
Polished outreach — especially video — pulls toward drama: swelling strings, "imagine a world
where…", a genome becoming a cure in thirty seconds. That compression is precisely where honesty
erodes. The worst outcome is not a bad artifact; it is **one more impressive than the honesty
ledger allows**, which detonates with experts and arms critics. The thing that makes the work
credible is the same thing most at risk in an edit bay or a headline.

> **Hard rule: when drama and accuracy conflict, accuracy wins.**

## The two-cut principle
The beat that thrills the general public is often the one that makes a specialist wince, and vice
versa. When one artifact must serve both, make **two versions** rather than one compromised one:
- a **technical cut** — for builders and experts: tools running, reproducibility, honest limits;
- a **mission cut** — for the broad audience and families: human stakes, understated.

Don't force one artifact to carry contradictory jobs.

## How to use this (workflow)
1. **Identify the persona(s)** the artifact primarily targets.
2. **Lead with what wins** that persona; **avoid what loses** them (table above).
3. **Apply the universal guardrails** to everything, regardless of persona.
4. **Flag any claim that risks the honesty line** *before it ships* — check it against the demo
   honesty ledger (`demo-foundation-alignment`) and cut or qualify anything unprovable.
5. If it must serve opposed audiences, **split into a technical cut and a mission cut.**

## The bottom line
Outreach converts a wave of attention into forks, collaborators, and eventually the downstream work
that helps someone — but only if it makes the right people curious *and* leaves experts trusting
what they find. **Hold the honesty discipline hardest exactly where it is hardest to hold — in the
places built to impress.**

## Related
- `hcls-core-vision-mission` — the mission this messaging protects (honesty is load-bearing).
- `demo-foundation-alignment` — what gets shown, and the honesty ledger to check claims against.
- `build-housekeeping-standards` — the "real, never mocked" discipline outreach depends on.
- *(Future: per-persona drill-down skills — `persona-domain-experts`, `persona-builders-oss`,
  `persona-patients-families`, `persona-skeptics`, `persona-ai-curious-public` — for deeper,
  audience-specific guidance.)*
