---
name: hcls-core-vision-mission
description: >-
  The north star for the HCLS AI Factory — its founding mission, vision, and principles. ALWAYS
  consult this first, on every task, before any other skill: it governs WHAT to build, HOW to
  frame it, and WHAT to promise, so that every decision stays aligned with the reason the project
  exists. When a technical choice and the mission seem to conflict, the mission wins. Every other
  skill (build-housekeeping-standards and the 16 pillar skills) exists to serve this one.
---

# HCLS AI Factory — Core Vision & Mission (North Star)

> **The mission, in one sentence:**
> *No one should die of a disease we could have understood in time.* — child or adult, in a
> leading research hospital or a clinic with no genetics department at all.

Every engine, every agent, every line of the platform traces back to that sentence. It is the
reason the project exists, and it is the test every decision is measured against.

## Why this exists
Precision medicine already works — but in pieces, living in separate silos, on separate machines,
behind separate budgets. The patient who needs all of them at once rarely gets all of them at
once. That gap — of access, integration, and imagination about what our tools can already do — is
the failure the HCLS AI Factory answers. It began in 2012 and was built on 20,000+ hours of
independent, self-funded research, carrying the discipline of high-performance computing into
health care. The founding lesson from HPC: **integration and access are not afterthoughts — they
are the whole game.**

## The vision
An integrated, multimodal precision-medicine platform that runs end to end **on a single
affordable machine**, entirely open-source under **Apache-2.0**: patient genome → interpreted,
actionable picture → reasoning across images and genomes → candidate molecules against real
targets, structure-validated and trial-matched — governed, reproducible, and honest about what is
proven versus what is still a frontier. Eight engines, eight intelligence agents, one open
platform on a common genomic and molecular foundation.

But the deepest thing it is meant to be is **not a platform — it is a proof, and an invitation**:
proof of what AI can do in health care when applied with real domain depth to a problem that
matters, and an invitation for anyone who sees a genome become a molecule on their own machine to
walk away with a ceiling removed from their imagination. It was released the way Linux was — a
working foundation, a clear vision, and a license that lets others carry it further than any one
person could alone.

## Founding principles — what is asked of whoever builds on this

1. **Take it and make it yours.** Fork it, rebuild it, disagree with it, surpass it. *The vision
   succeeds when the work no longer needs its author.*
2. **Keep it honest.** The power of this is that what people see is real. Never overpromise. Label
   what is preclinical, what is roadmap, and what is decision-support rather than diagnosis.
   *Overstatement betrays the mission more than any missing feature ever could.*
3. **Keep the foundation free.** Apache-2.0 on purpose — take it, build a business on it, close
   your own version if you want; nothing is owed back. Open-source here does **not** mean
   non-commercial. If it accelerates precision medicine — in a product, a clinic, or a startup —
   the mission is met. The only thing locked permanently open is the foundation itself.
4. **Keep it for everyone.** It started with children because a child dying of preventable disease
   is the sharpest moral case there is. Carry it to adults, the elderly, the underserved — until
   *"for all"* is not a slogan but a description.

## How the north star guides the work

Translate the vision into every decision:

- **Honesty is load-bearing, not a nicety.** A `live` capability is never mock-served; mocks are
  labeled; clinical claims cite sources; outputs are decision-support, not diagnosis. This is *why*
  the governance/honesty gates exist — they protect the mission. (See `11-security-and-secrets`,
  `08-inference-serving`, `build-housekeeping-standards`.)
- **Access and affordability are features.** It must keep running end-to-end on one affordable box
  and stay deployable by a clinic with no genetics department. Favor local-first, one-command
  bring-up, minimal prerequisites, and sensible defaults over anything that raises the barrier.
  (See `01-compute-dgx-spark-remote-gpus`, `12-cost-and-economics`, `14-ease-of-deployment`.)
- **Openness and clean provenance are the walls, not the trim.** Permissive licensing, honest
  documentation, a neutral repo, and reproducibility are what make the gift real and forkable —
  treat them as core, never optional. (See `15-github-structure-and-presentation`.)
- **Integration is the point.** The enemy is the silo. Every new engine/agent should compose into
  the whole through the registry and the workflow composer, not stand alone. (See
  `09-ai-orchestration`.)
- **Build for "for all."** Prefer capabilities and framings that widen reach — more conditions,
  more populations, lower cost — over ones that only serve the well-resourced.

## The measure
This was never meant to be a finish line. Whether it made the world better is measured downstream
— a child diagnosed sooner, a therapy matched in time, a cost lowered enough that a clinic without
a genetics department can still give a family an answer. Some of that will happen in hands that are
not the author's, which is exactly as it should be. The charge for anyone building here:

> **No one should die of a disease we could have understood in time. Take this, and help make that true.**

## Related
- `build-housekeeping-standards` — how we build (this is *why* we build).
- The 16 Core-Pillar skills under `architecture/` — the technical standards that serve this vision.
- `skills/AI Factory Core Pillars.txt` — the pillars this mission is delivered through.
