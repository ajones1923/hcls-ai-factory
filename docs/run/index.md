---
title: Run It Yourself — Quickstart
description: Clone the HCLS AI Factory, run the build gate, and verify a real capability — the one claim a competitor cannot dismiss.
---

# Run It Yourself — Quickstart

!!! note "This page is for builders"
    It's for people who want to clone and run the factory themselves. If you'd rather just *see* what
    it produces, watch the [three-minute tour](../index.md) or the [Demonstrations](../demos/index.md).

Open-source means *you can run this yourself* — the one claim no one can wave away. Everything below
is real; heavier or ARM-incompatible models burst to remote GPU over a private mesh
([elastic burst](hardware.md)), and the site says so wherever that's the case. Apache-2.0,
**commercial use welcome.**

## Clone it

```bash
git clone https://github.com/ajones1923/hcls-ai-factory
cd hcls-ai-factory
```

## Run the platform gate

The same checks that gate every merge — real tests, no `|| true`:

```bash
ruff check --select E9,F82,F811,F706,F707 core lib scripts   # real-bug lint
( cd lib/hcls_common && pytest -q )                          # platform library tests
python scripts/validate_registry.py                          # capability-registry coverage
```

## What "live" means here

A capability marked `live` on the **[Capability Maturity Matrix](../honesty/maturity-matrix.md)** is
served by a real model against real input — never a mock. Where a capability is `planned`,
`preclinical`, or `gated`, it is labeled as such and not shown as if it ships. See the
[Honesty Ledger](../honesty/ledger.md) for what each badge means.

## Proof it's real

- **A real, passing test suite gates every merge** — the platform library's hundreds of tests, plus
  real-bug lint and registry validation, run for real (no `|| true`, no skipped gate).
- **The registry validates** — 41 registered capabilities, every engine and agent directory accounted
  for; a `live` capability can never be mock-served.
- **Verified on real data** — capabilities marked `verified` were run against real public reference
  data and are cited on the [Citations & Evidence](../honesty/citations.md) page: e.g. the variant
  store on **GIAB HG002** (Ts/Tv ≈ 2.0), single-cell on **pbmc3k**, ACMG secondary-findings on real
  **ClinVar**, ESM-2 search on real proteins, and ADMET on real drug molecules.
- **The site is honest by construction** — this page, the
  [Capability Maturity Matrix](../honesty/maturity-matrix.md), and every engine/agent page are
  generated from the registry, and the strict CI build fails before publishing anything stale.

*Next: the honest [hardware & elastic-burst story](hardware.md) · [reproducibility & datasets](reproducibility.md).*
