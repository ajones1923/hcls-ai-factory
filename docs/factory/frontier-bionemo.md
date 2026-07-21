---
title: NVIDIA BioNeMo
description: NVIDIA BioNeMo provides GPU-accelerated molecular models as inference microservices (NIMs) — MolMIM, DiffDock, and GenMol — powering the Therapeutic Discovery Engine.
---

# NVIDIA BioNeMo

[NVIDIA BioNeMo](https://www.nvidia.com/en-us/clara/bionemo/) packages biology models as
GPU-accelerated inference microservices (**NIMs**) that the factory calls on demand. Its models power
the [Therapeutic Discovery Engine](engines/therapeutic-discovery-engine.md)'s *generate → dock → score*
loop, connected as one of the [Frontier Models](frontier-models.md) via elastic burst.

![NVIDIA BioNeMo in the factory: MolMIM, DiffDock, and GenMol feeding the Therapeutic Discovery Engine](../assets/infographics/frontier-bionemo.png)
/// caption
Illustrative. MolMIM and DiffDock are live clients (they need the NIM container); GenMol is planned.
///

## Where it plugs in

| Model | Role | Status |
|---|---|---|
| **MolMIM** | Generative small-molecule design from a seed compound. | <span class="cap-badge cap-live">live</span> |
| **DiffDock** | Protein–ligand docking and 3-D pose prediction. | <span class="cap-badge cap-live">live</span> |
| **GenMol** | SAFE-fragment molecule generation (a diversity complement to MolMIM). | <span class="cap-badge cap-planned">planned</span> |

All three serve the **Therapeutic Discovery Engine** — MolMIM and DiffDock as the generation and
docking stages, GenMol as a planned diversity generator.

## Honest limits

- **Real clients, real dependency.** MolMIM and DiffDock are `live` — real client code — but they need
  the **NVIDIA NIM container** (NGC + GPU) deployed to actually run; a clearly-labeled mock fallback is
  used only in development.
- **The default generator is open.** The engine's *default, verified* molecule generator is RDKit
  BRICS; MolMIM is the GPU-accelerated alternative when its NIM is deployed.
- **Elastic burst.** These are GPU microservices that run off the box, on demand — never pretending to
  run locally.
- **Preclinical design, not drugs.** Generated molecules and docking poses are preclinical design
  candidates — decision support for a scientist, not therapies.
