---
title: Frontier Models
description: The external models the factory builds on, in two roles — the reasoning layer (Anthropic's Claude) and the heavy domain-science models (NVIDIA BioNeMo, Chai Discovery, and the roadmap frontier) — each a registered capability with an honest status.
---

# Frontier Models

The factory builds on external models that no single box could hold alone — in **two distinct roles**:
the **reasoning layer** that interprets and orchestrates, and the **domain-science models** that
generate molecules and fold structures. Each connects **on demand** and registers as a native
capability with an **honest status**, so nothing here pretends to run locally, and nothing claims more
than it delivers.

## The reasoning layer

- **[Anthropic (Claude)](anthropic.md)** — the reasoning LLM behind the Precision Intelligence
  Engine and all eight intelligence agents. The "intelligence" in *intelligence agent*. Reached as a
  hosted API — the specific model is **configurable**, and the layer degrades honestly (never
  mock-served) when a key is absent.

## Domain-science models

Heavy or partner-gated models — too large for one box, or available only through a partner — reached
**on demand via [elastic burst](../../run/hardware.md)** to a remote GPU:

- **[NVIDIA BioNeMo](bionemo.md)** — GPU molecular microservices: **MolMIM** & **DiffDock**
  (<span class="cap-badge cap-live">live</span>) and **GenMol** (<span class="cap-badge cap-planned">planned</span>) for the Therapeutic Discovery Engine.
- **[Chai Discovery](chai-discovery.md)** — **Chai-1** (complex co-folding · <span class="cap-badge cap-planned">planned</span>)
  and **Chai-2** (de-novo binder & antibody design · <span class="cap-badge cap-gated">gated</span>):
  predict a molecular complex, or design a brand-new binder — and chain the two (design → validate).

### Also on the roadmap

Registered and connected the same way — heavy or containerized, honest status, elastic burst:

- **RFdiffusion** &nbsp; <span class="cap-badge cap-planned">planned</span> — de-novo protein backbone design (→ ProteinMPNN → Chai-1).
- **Evo 2** &nbsp; <span class="cap-badge cap-planned">planned</span> — genome-scale foundation model.

!!! note "Honest by construction"
    The statuses here are the registry's — see the live
    [Capability Maturity Matrix](../../honesty/maturity-matrix.md). The **domain-science models burst to
    a remote GPU** ([elastic burst](../../run/hardware.md)), never pretending to run on the box; the
    **reasoning layer** is a hosted API — configurable, and never mock-served. Any clinical use is
    **decision support for a qualified clinician**, never autonomous diagnosis.
