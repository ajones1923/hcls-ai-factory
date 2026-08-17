---
title: Roadmap
description: What is live, planned, and preclinical in the HCLS AI Factory — and how the public site grows honestly alongside the platform.
---

# Roadmap

The factory and this site grow together, in public. The **[Capability Maturity
Matrix](honesty/maturity-matrix.md)** is the live scorecard; this page is the honest direction.

![The site goes green as the factory matures — roadmap to planned to live to verified](assets/infographics/roadmap.png)
/// caption
Every capability advertises its real status, generated from the registry — the site cannot claim ahead of the code.
///

## Platform

- **Live today** — the platform library (capability registry, MCP tool-surface, workflow composer,
  governance gates) and the engines/agents the registry declares `live`.
- **Planned** — the frontier co-folding and generation models registered `planned` (e.g. Chai-1,
  GenMol) flip to `live` only when they are stood up and served for real. (Chai-2 de-novo binder
  design is `gated` — partnership access — not `planned`.)
- **Preclinical** — gene-therapy directions (including TSC1/TSC2) are an open design/analysis bench,
  not a treatment available today.

## This site

Built in four phases, each flipping matrix badges green as real capabilities come online:

| Phase | Focus |
|---|---|
| **0 — Foundation** *(shipped)* | design-system + honest status badges, OG/SEO + `llms.txt`, the registry-bound maturity matrix, the section skeleton, a strict docs build in CI |
| **1 — Depth** *(shipped)* | per-engine / per-agent / TSC pages, published diagrams + narrated videos, badges flipping as capabilities verify |
| **2 — Proof** *(shipped)* | "Run It Yourself" quickstart, verified-on-real-data proofs, citations, governance |
| **3 — Launch** *(in progress)* | performance + accessibility pass, analytics for the right metrics, the launch moment |

## Changelog

The site build, in public — most-recent first.

- **2026-08** — **Security & governance:** authentication available on **12 of 12** clinical
  endpoints via one shared module, fail-closed once a key is set; the `X-HCLS-Governed` header now
  reports only gates that actually ran, and [Governance & Lineage](honesty/governance.md) states
  per-service gate coverage as a number rather than a claim.
- **2026-08** — **Proof of health:** one harness across all 17 subjects — 8 engines, 8 agents, the
  TSC program — **8,402 tests passing**; compose coverage 11/17 → 17/17; a single port convention
  (registry advertises the UI, API is UI + 1) now cross-checked against the process supervisor in CI.
- **2026-08** — **Site pass:** every page's narrated video re-recorded in one voice; **84
  infographics audited against the registry** and nine corrected where they claimed ahead of the
  code; per-capability subject guides published; the preprint added.

- **2026-07** — **Verification campaign:** ran capabilities on real public data and earned `verified`
  where they passed — **nine in total**: the variant store (HG002), single-cell compute (pbmc3k),
  ACMG secondary-findings (real ClinVar), the molecule generator + ADMET (real drugs), ESM-2 search +
  ProteinMPNN + ESMFold (real proteins), and the reproducibility capture.
- **2026-07** — **Fuller IA:** two-tier tab navigation, a page for every engine and agent, About /
  Mission, Honesty Ledger + Decision-Support Posture pages, Run It Yourself split
  (Quickstart / Hardware / Reproducibility), per-page social OG cards, `CITATION.cff`.
- **2026-07** — **Honesty audit:** code-level review of every `live` capability; corrected the real
  overclaims (`live → planned` where a dependency can't load); the maturity vocabulary shipped.
- **2026-07** — **Foundation (Phases 0–3):** the registry-bound maturity matrix, the honesty spine,
  citations & governance pages, `llms.txt` + SEO, the design-system status badges, a strict docs build
  gating every PR.

## A note on the status vocabulary

Beyond the serving states **live** and **planned**, the registry now carries an optional, orthogonal
**maturity tier** — **verified · preclinical · research-use · roadmap · gated** — machine-checked
against the code (e.g. `verified` is rejected by the registry on anything not actually running). It
is populated **conservatively**, only where the registry documents a basis: today **nine capabilities
are `verified`** (proven on real public data — the variant store on HG002, single-cell compute on
pbmc3k, ESMFold, ESM-2 search, ProteinMPNN, ACMG secondary findings, the molecule generator, ADMET,
and reproducibility), and Chai-2 is `gated` (partnership access).
Broadening coverage — and adding finer-grained capability entries for the feature-level caveats now
carried on the [Honesty & Governance ledger](honesty/index.md) (e.g. TSC gene therapy is preclinical)
— continues as capabilities are audited.
