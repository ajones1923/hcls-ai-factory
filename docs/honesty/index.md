---
title: Honesty & Governance
description: The honesty discipline is the load-bearing wall of the HCLS AI Factory — the ledger, the decision-support posture, and how the site stays honest by construction.
---

# Honesty & Governance

The honesty discipline is the load-bearing wall of this project. The site is built so it **cannot
claim more than the code delivers**: statuses are rendered from the capability registry, and every
limit is stated plainly where it matters.

- **[Capability Maturity Matrix](maturity-matrix.md)** — every capability's real status, generated
  at build time from `lib/hcls_common/capabilities.json`.

## The four standing commitments

- **Decision support, not diagnosis.** Every clinical output supports a qualified clinician — never
  autonomous diagnosis or prescribing.
- **A bench, not a cure.** Where the factory points toward gene therapies, those are **preclinical**
  research — the promise of a direction, not a treatment available today.
- **"Elastic burst," not "all on one box."** Heavy or ARM-incompatible models burst to remote GPU
  over a private mesh; nothing pretends to fit on the desk that doesn't.
- **Real, never mocked.** A capability marked `live` is served by a real model against real input.

## The honesty ledger

Carried everywhere it is relevant, stated out loud:

| Item | Honest status |
|---|---|
| Gene therapy for TSC1/TSC2 (and gene correction generally) | **preclinical** — open design/analysis bench, not a cure today |
| MAISI synthetic imaging | **research / augmentation / QA use** — never a diagnostic source |
| Single-cell atlas similarity & cell-foundation embeddings | **roadmap** |
| Chai-2 de novo binder/antibody design | **gated** — partnership access; integration contingent |
| α-synuclein SAA · plasma p-tau217 · NSD-ISS / SynNeurGe staging | **research- / trial-use** inputs the agents reason over — not routine diagnostics |
| Frontier co-folding (Chai-1, RFdiffusion, Evo 2) | run via **elastic burst** to remote GPU, not locally |
| All clinical output | **decision support**, not autonomous diagnosis or prescribing |

## Capability audit

A **code-level audit of every `live` capability** (2026-07) checked the actual implementation behind
each registry claim — not the claim itself. The findings are reflected directly in the
[maturity matrix](maturity-matrix.md) and each capability's description:

- **The platform is genuinely implemented** — the library services (registry, MCP surface, workflow
  composer, MLOps, honesty gate), the genomics pipeline stages, all eight RAG intelligence agents,
  and the engines are real serving code with tests, not scaffolding.
- **`verified`** is reserved for capabilities proven against real recorded input: the variant store
  (real GIAB **HG002**, Ts/Tv ≈ 2.0) and **ProteinMPNN** design (vendored weights + a logged real
  design run + a real-run test).
- **Two capabilities were corrected `live` → `planned`** because their own code cannot currently run
  on the reference box (a required library isn't installed): MHCflurry immunogenicity and the ESM-2
  LoRA fine-tune.
- **Operational-dependency caveats are stated plainly** where the code runs but needs a runtime
  service: the RAG agents need a populated vector DB + an API key; the MolMIM/DiffDock NIMs need their
  NVIDIA container deployed (with a clearly-labeled mock fallback for development, never presented as
  real).

!!! note "Honest by construction"
    The maturity matrix is **generated at build time from the registry**, not hand-maintained prose —
    so it cannot drift from the code. That generation is the source of every status shown here; if the
    registry can't be read or is malformed, the strict build fails rather than publishing a stale page.
    A `live` capability is never mock-served, and `verified` cannot sit on anything that isn't actually
    running — both enforced in code.
