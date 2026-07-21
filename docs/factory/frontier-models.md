---
title: Frontier Models
description: Heavy or partner-gated models the factory connects to on demand via elastic burst — Chai Discovery's Chai-1 and Chai-2, and the roadmap frontier — each a registered capability with an honest status.
---

# Frontier Models

Some of the most capable models in biology are too heavy for one box, or available only through a
partner. The factory connects to them **on demand, via [elastic burst](../run/hardware.md)** — each
registers as a native capability in the registry with an **honest status**, so nothing here pretends to
run locally, and nothing claims more than it delivers.

## Connected services

- **[Anthropic (Claude)](frontier-anthropic.md)** — the reasoning LLM behind the Precision Intelligence
  Engine and all eight intelligence agents. The "intelligence" in *intelligence agent*.
- **[NVIDIA BioNeMo](frontier-bionemo.md)** — GPU molecular microservices: **MolMIM** & **DiffDock**
  (<span class="cap-badge cap-live">live</span>) and **GenMol** (<span class="cap-badge cap-planned">planned</span>) for the Therapeutic Discovery Engine.
- **[Chai Discovery](frontier-chai.md)** — **Chai-1** (complex co-folding · <span class="cap-badge cap-planned">planned</span>)
  and **Chai-2** (de-novo binder & antibody design · <span class="cap-badge cap-gated">gated</span>):
  predict a molecular complex, or design a brand-new binder — and chain the two (design → validate).

## Also on the roadmap

Registered and connected the same way — heavy or containerized, honest status, elastic burst:

- **RFdiffusion** &nbsp; <span class="cap-badge cap-planned">planned</span> — de-novo protein backbone design (→ ProteinMPNN → Chai-1).
- **Evo 2** &nbsp; <span class="cap-badge cap-planned">planned</span> — genome-scale foundation model.

!!! note "Honest by construction"
    The statuses here are the registry's — see the live
    [Capability Maturity Matrix](../honesty/maturity-matrix.md). Frontier models **burst to a remote
    GPU** ([elastic burst](../run/hardware.md)), never pretending to run on the box; and any clinical
    use is **decision support for a qualified clinician**, never autonomous diagnosis.
