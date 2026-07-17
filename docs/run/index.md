---
title: Run It Yourself
description: Clone the HCLS AI Factory, run the build gate, and verify a real capability — the one claim a competitor cannot dismiss.
---

# Run It Yourself

Open-source means *you can run this yourself* — the one claim no one can wave away. Everything below
is real; heavier or ARM-incompatible models burst to remote GPU over a private mesh (**elastic
burst**), and the site says so wherever that's the case.

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
`preclinical`, or `gated`, it is labeled as such and not shown as if it ships.

!!! note "Hardware & honesty"
    The reference target is a single **NVIDIA DGX Spark** (ARM Grace + Blackwell). Heavy frontier
    models run on remote GPUs on demand — **"elastic burst," not "all on one box."** All patient data
    stays on the local box; only derived, non-identifying work is sent to a burst GPU.

## Proof it's real

- **A real, passing test suite gates every merge** — the platform library's hundreds of tests, plus
  real-bug lint and registry validation, run for real (no `|| true`, no skipped gate).
- **The registry validates** — 41 registered capabilities, every engine and agent directory accounted
  for; a `live` capability can never be mock-served.
- **Verified on real data** — the [queryable variant store](../factory/engines/genomics-engine.md) is
  marked `verified` because it was run on real **GIAB HG002** reads: PASS-SNV **Ts/Tv ≈ 2.0** (the
  expected genome-wide human range ≈ 2.0–2.1), behind a QC trust-gate that *withholds* interpretation
  on a bad call set. A real result on real data, not a mock — see [Citations & Evidence](../honesty/citations.md).
- **The site is honest by construction** — this page, the
  [Capability Maturity Matrix](../honesty/maturity-matrix.md), and every engine/agent page are
  generated from the registry, and the strict CI build fails before publishing anything stale.
