---
title: HCLS AI Factory
template: home.html
hide:
  - navigation
  - toc
---

## The factory in under a minute

<video class="cap-video" controls preload="metadata" playsinline poster="/assets/factory-architecture.png" src="/assets/videos/factory-overview.mp4">
  Your browser can't play embedded video — <a href="/assets/videos/factory-overview.mp4">download the overview</a>.
</video>
/// caption
A narrated, captioned tour of the whole factory. Decision support for a qualified clinician, never diagnosis.
///

## Inside the factory

Eight compute engines, eight clinical intelligence agents, and one flagship disease program — on a
single [NVIDIA DGX Spark](https://www.nvidia.com/en-us/products/workstations/dgx-spark/), open under
Apache-2.0, bursting heavier work to remote GPUs over a private mesh only when it needs to.

<div class="grid cards" markdown>

- :material-engine-outline:{ .lg } __[8 Engines](factory/engines/index.md)__

    The horizontal compute muscle — genomics, clinical interpretation, therapeutic discovery,
    imaging, oncology, cardiology, protein design, and single-cell.

- :material-brain:{ .lg } __[8 Intelligence Agents](factory/agents/index.md)__

    The clinical reasoning layers that interpret what the engines produce — each one decision
    support for a qualified clinician.

- :material-heart-pulse:{ .lg } __[1 Disease Program](programs/tsc.md)__

    Tuberous Sclerosis Complex — the flagship vertical that composes the whole factory for one
    child, in one governed afternoon.

</div>

## Explore

- **[The Factory](factory/overview.md)** — the Platform → Engine → Agent → Program model, with a
  registry-generated page for every capability.
- **[Capability Maturity Matrix](honesty/maturity-matrix.md)** — every capability's real, audited
  status, generated from the registry. The site cannot claim ahead of the code.
- **[Run It Yourself](run/index.md)** — clone it, run the gate, verify a real capability.
- **[Demonstrations](demos/index.md)** — the D1–D7 portfolio, each with its honesty flags.
- **[Capability Brief](brief/README.md)** — the full picture in a **technical cut** and a
  **mission cut**.
- **[Source on GitHub](https://github.com/ajones1923/hcls-ai-factory)** — clone it, read it, run it.

## What it is — honestly

The honesty discipline is the load-bearing wall of this project, so the limits are stated plainly
everywhere they matter:

- **Decision support, not diagnosis.** Every clinical output supports a qualified clinician — never
  autonomous diagnosis or prescribing.
- **A bench, not a cure.** Where the factory points toward gene therapies, those are **preclinical**
  research — the promise of a direction, not a treatment available today.
- **"Elastic burst," not "all on one box."** Heavy or ARM-incompatible models burst to remote GPU
  over a private mesh; nothing pretends to fit on the desk that doesn't.
- **Real, never mocked.** A capability marked `live` is served by a real model against real input —
  and a [code-level audit](honesty/index.md#capability-audit) keeps that claim earned.

!!! note "Source of truth"
    This site reflects the canonical roster — **8 engines · 8 intelligence agents · 1 disease
    program (TSC)** — generated from the capability registry (`lib/hcls_common/capabilities.json`).
