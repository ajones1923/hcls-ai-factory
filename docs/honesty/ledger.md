---
title: Honesty Ledger
description: The badge vocabulary and every non-live scope, stated plainly — the source of truth for what each status on the site means.
---

# Honesty Ledger

This page is the definition source for the badges you see across the site. Every badge on the
[Capability Maturity Matrix](maturity-matrix.md) and on each engine/agent page is a **deterministic
function of the capability registry plus this ledger** — never hand-typed prose.

## Badge vocabulary

| Badge | Meaning | Derivation |
|---|---|---|
| <span class="cap-badge cap-live">live</span> | Real and running — served by a real model against real input, never mock-served. | registry `status: live` |
| <span class="cap-badge cap-verified">verified</span> | Additionally **proven against real, recorded input** (e.g. GIAB HG002, pbmc3k, real ClinVar). | registry `status: live` **and** a recorded real-data run |
| <span class="cap-badge cap-planned">planned</span> | Not yet running (or a required dependency isn't deployed). | registry `status: planned` |
| <span class="cap-badge cap-verified">verified</span> | *(maturity tier)* proven against real recorded input. | registry `maturity: verified` |
| <span class="cap-badge cap-preclinical">preclinical</span> | A research bench, not a treatment or clinical use today. | registry maturity tier / honesty ledger |
| <span class="cap-badge cap-research-use">research-use</span> | Research / QA / trial-use input, not a routine clinical diagnostic. | registry maturity tier / honesty ledger |
| <span class="cap-badge cap-roadmap">roadmap</span> | A designed direction the factory is heading — not yet built. | registry maturity tier / honesty ledger |
| <span class="cap-badge cap-gated">gated</span> | Real but partnership- or license-gated access. | registry `maturity: gated` |

## The standing ledger

Carried everywhere it is relevant, stated out loud:

| Item | Honest status |
|---|---|
| Gene therapy for TSC1/TSC2 (and gene correction generally) | **preclinical** — an open design/analysis bench, not a cure today |
| MAISI synthetic imaging | **research / augmentation / QA use** — never a diagnostic source |
| Single-cell **atlas-similarity search & cell-foundation embeddings** (distinct from the `verified` single-cell compute engine) | **roadmap** |
| Chai-2 de novo binder / antibody design | **gated** — partnership access; integration contingent |
| α-synuclein SAA · plasma p-tau217 · NSD-ISS / SynNeurGe staging | **research- / trial-use** inputs the agents reason over — not routine diagnostics |
| Frontier co-folding (Chai-1, RFdiffusion, Evo 2) | when stood up, run via **elastic burst** to remote GPU, never locally — Chai-1 is `planned` |
| MolMIM / DiffDock NIMs | real clients; require the NVIDIA NIM container deployed (labeled mock fallback for dev) |
| All clinical output | **decision support**, not autonomous diagnosis or prescribing |

## Why it's trustworthy

Because the matrix is **generated at build time from the registry**, the site cannot show a status
the code does not declare, and a `verified` badge cannot sit on anything that isn't actually running —
both enforced in code. If the registry can't be read or is malformed, the strict build fails rather
than publishing a stale page. See the [capability audit](index.md#capability-audit) for how every
`live` claim was checked against the actual implementation.
