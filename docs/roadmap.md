---
title: Roadmap
description: What is live, planned, and preclinical in the HCLS AI Factory — and how the public site grows honestly alongside the platform.
---

# Roadmap

The factory and this site grow together, in public. The **[Capability Maturity
Matrix](honesty/maturity-matrix.md)** is the live scorecard; this page is the honest direction.

## Platform

- **Live today** — the platform library (capability registry, MCP tool-surface, workflow composer,
  governance gates) and the engines/agents the registry declares `live`.
- **Planned** — the frontier co-folding and generation models registered `planned` (e.g. Chai-1,
  Chai-2, GenMol) flip to `live` only when they are stood up and served for real.
- **Preclinical** — gene-therapy directions (including TSC1/TSC2) are an open design/analysis bench,
  not a treatment available today.

## This site

Built in four phases, each flipping matrix badges green as real capabilities come online:

| Phase | Focus |
|---|---|
| **0 — Foundation** *(in progress)* | design-system + honest status badges, OG/SEO + `llms.txt`, the registry-bound maturity matrix, the section skeleton, a strict docs build in CI |
| **1 — Depth** | per-engine / per-agent / TSC pages, published diagrams, badges flipping as capabilities verify |
| **2 — Proof** | "Run It Yourself" quickstart, verified-on-real-data proofs, citations, governance |
| **3 — Launch** | performance + accessibility pass, analytics for the right metrics, the launch moment |

## A note on the status vocabulary

Beyond the serving states **live** and **planned**, the registry now carries an optional, orthogonal
**maturity tier** — **verified · preclinical · research-use · roadmap · gated** — machine-checked
against the code (e.g. `verified` is rejected by the registry on anything not actually running). It
is populated **conservatively**, only where the registry documents a basis: today the queryable
variant store is `verified` (proven on real HG002) and Chai-2 is `gated` (partnership access).
Broadening coverage — and adding finer-grained capability entries for the feature-level caveats now
carried on the [Honesty & Governance ledger](honesty/index.md) (e.g. TSC gene therapy is preclinical)
— continues as capabilities are audited.
