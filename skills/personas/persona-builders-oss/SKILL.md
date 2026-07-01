---
name: persona-builders-oss
description: >-
  Deep drill-down for reaching builders and the open-source community — developers, ML/bio
  engineers, HPC people, tinkerers. Consult when an outreach artifact primarily targets home turf:
  the README, repo landing, GitHub release notes, dev-facing posts, "clone and run" demos. Pairs
  with broad-general-persona and the north star's open-foundation principle.
---

# Builders / Open-Source — home turf

Developers, ML and bio engineers, HPC people, and tinkerers. This is the audience most like the
project itself, and their reaction is binary: they either **clone it and become the community**, or
they **find the cracks** and say so publicly. They are won or lost by the repo, not the pitch — a
slick promo over a messy or non-reproducible repo loses them faster than no promo at all.

## What wins them / What loses them

**Wins them:** a genuinely forkable, documented repo — clean provenance, a real license, one
affordable box, one-command bring-up, and honest "here's what runs, here's what bursts." The
**Linux-style framing** lands hard: a working foundation released for others to carry further.
Reproducibility they can verify is the whole game. And **commercial use explicitly welcomed** —
Apache-2.0, build a business on it, close your fork, nothing owed back.

**Loses them / the risk:** promising forkability the repo doesn't yet deliver. If the demo is
polished but `git clone` fails, the docs are stale, or "runs on one box" quietly means "plus a
cluster you don't mention," the gap becomes the story. Overselling reproducibility to this
audience is worse than underselling it.

## How to speak to them

- Lead with **"runs on one box, Apache-2.0, here's the GitHub"** — concrete, checkable, no
  marketing throat-clearing. Link the repo above the fold.
- Make the **"and it's yours"** hook literal: fork it, rebuild it, disagree with it, surpass it;
  the vision succeeds when the work no longer needs its author.
- Be explicit that **open-source here is not non-commercial** — Apache-2.0 welcomes products,
  clinics, and startups. A business accelerating precision medicine is the goal succeeding.
- **Show it running, honestly.** Name what is live vs. what bursts to remote GPUs ("elastic
  burst," never imply all-local). Clean provenance and neutral naming signal a repo worth trusting.
- Only promise the forkability the repo actually delivers *today*. If a path is rough, say so and
  file it as a known gap — builders forgive honest rough edges and punish hidden ones.

## Do / Don't

**Do:** put the repo and license first; make bring-up reproducible and documented; welcome
commercial use out loud; label burst vs. local; treat rough edges as honest TODOs. **Don't:**
ship promo polish over a repo that won't clone/run; overstate reproducibility; obscure the
architecture; imply one-box when it isn't; bury the license question.

## Related
- `broad-general-persona` — the five-persona map; builders are home turf.
- `hcls-core-vision-mission` — "take it and make it yours" and the open-foundation principle.
- `demo-foundation-alignment` — the live-not-mock honesty rule the repo's credibility rests on.
