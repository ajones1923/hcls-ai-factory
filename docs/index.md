# HCLS AI Factory

> **No one should die of a disease we could have understood in time.**

An open-source precision-medicine platform that takes a patient's raw DNA to a shortlist of
therapeutic candidates **in hours** — on a single [NVIDIA DGX Spark](https://www.nvidia.com/en-us/products/workstations/dgx-spark/),
open under Apache-2.0. One box you can fork and run yourself, bursting heavier work to remote GPUs
over a private mesh only when it needs to.

<div class="grid cards" markdown>

- :material-engine-outline:{ .lg } __8 Engines__

    The horizontal compute muscle — genomics, discovery, imaging, oncology, cardiology, protein
    design, and single-cell.

- :material-brain:{ .lg } __8 Intelligence Agents__

    The clinical reasoning layers that interpret what the engines produce — each one decision
    support for a qualified clinician.

- :material-heart-pulse:{ .lg } __1 Disease Program__

    Tuberous Sclerosis Complex — the flagship vertical that composes the whole factory for one
    child, in one governed afternoon.

</div>

## Start here

- **[Capability Brief](brief/README.md)** — the full picture: every engine and agent, the flagship
  TSC program, and the D1–D7 demo portfolio. Available as a **technical cut** (experts and builders)
  and a **mission cut** (the human "why").
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
- **Real, never mocked.** A capability marked `live` is served by a real model against real input.

!!! note "Source of truth"
    This site reflects the canonical roster — **8 engines · 8 intelligence agents · 1 disease
    program (TSC)** — reconciled with the capability registry (`lib/hcls_common/capabilities.json`).
