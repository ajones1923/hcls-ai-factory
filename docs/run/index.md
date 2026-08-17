---
title: Run It Yourself — Quickstart
description: Clone the HCLS AI Factory, set up the environment, run the real build gate, and execute a live capability on real public data — in about ten minutes, with no GPU and no gated software.
---

# Run It Yourself — Quickstart

!!! note "This page is for builders"
    It's for people who want to clone and run the factory themselves. If you'd rather just *see* what
    it produces, watch the [three-minute tour](../index.md) or the [Demonstrations](../demos/index.md).

Open-source means *you can run this yourself* — the one claim no one can wave away. Everything below
was run on a clean clone before it was published. Heavier or ARM-incompatible models burst to remote
GPU over a private mesh ([elastic burst](hardware.md)), and the site says so wherever that's the
case. Apache-2.0, **commercial use welcome.**

**No GPU and no gated software are needed for anything on this page.**

## 1. Clone it

```bash
git clone https://github.com/ajones1923/hcls-ai-factory
cd hcls-ai-factory
```

## 2. Set up the environment

This step is not optional. The shared library `hcls_common` must be installed before any gate will
run, and modern Linux marks the system Python "externally managed" (PEP 668), so it needs a venv:

```bash
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -e lib/hcls_common
.venv/bin/pip install duckdb statsmodels biopython peft     # two engines import these directly
```

Confirm it took:

```bash
.venv/bin/python -c "import hcls_common; print('ok')"
```

??? warning "If that prints `No module named 'hcls_common'`"
    A stale editable install can shadow it. Check
    `~/.local/lib/python3.12/site-packages/__editable___hcls_common_*_finder.py` — if its `MAPPING`
    dict is empty (`{}`), the install is corrupt. Delete both `__editable__*hcls_common*` files and
    repeat step 2. This exact failure blocked the project's own merge gate until it was found and
    documented.

## 3. Run the platform gate

The same checks that gate every merge — real tests, no `|| true`:

```bash
.venv/bin/python -m ruff check --select E9,F82,F811,F706,F707 core lib scripts
( cd lib/hcls_common && ../../.venv/bin/python -m pytest -q )
.venv/bin/python scripts/validate_registry.py
```

Expected:

```
All checks passed!
382 passed
registry: 42 capabilities, 9 typed 'engine'
OK — manifest valid and every engine/agent directory is registered.
```

`validate_registry.py` also enforces the [port convention](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md) — the registry
advertises the UI port and the API is UI + 1 — and cross-checks it against the process supervisor,
so a port change in one place and not the other fails the build.

## 4. Run all 17 subject suites

One command for the eight engines, eight intelligence agents and the TSC disease program:

```bash
.venv/bin/python scripts/run_all_tests.py
```

On a fresh clone:

```
17 subjects · passed 8194 · failed 0 · errors 0 · skipped 213
```

**The 213 skips are expected and correct.** The project deliberately never publishes `data/` —
`.gitignore` excludes it and the rule is *"data, weights and secrets stay local; only code and docs
publish."* Tests that assert on that local-only reference data skip rather than fail, so a clean
clone reports honestly instead of looking broken. On a machine that has the data, the same command
reports **8,402 passed**.

Use this script rather than a bare `pytest`: several subjects ship a module that shadows a Python
standard-library name, and two scientific packages register pytest plugins that abort collection.
The harness handles both, and explains why inline.

## 5. Run something real

The single-cell engine runs a genuine [scanpy](https://scanpy.readthedocs.io) workflow on the public
**PBMC 3k** dataset — on CPU, no gated model, no GPU:

```bash
.venv/bin/pip install scanpy anndata leidenalg igraph
.venv/bin/python scripts/run_demo.py E8
```

```
dataset       PBMC 3k via scanpy (local copy absent — data/ is not published)
loaded        2,700 cells x 32,738 genes
running       QC -> normalize -> HVG -> PCA -> neighbors -> Leiden -> marker DE
clusters      9
  cell type   B cells / CD14+ Monocytes / CD4 T cells / Dendritic cells
              FCGR3A+ Monocytes / Megakaryocytes / NK cells
RESULT        PASS — ran on real input
```

That is real computation on real public data — scanpy fetches PBMC 3k itself, so this reproduces on
any clone — recorded to `demo/transcripts/E8.txt` so you can diff it. See what else is ready:

```bash
.venv/bin/python scripts/run_demo.py --check-all
```

The runner **refuses to run a demo labelled LIVE whose service is unreachable** rather than quietly
returning a canned result — the same discipline that governs the maturity matrix.

## 6. Bring the platform up (optional)

```bash
cp .env.example .env      # then fill in the four required values
docker compose -f docker-compose.dgx-spark.yml up -d
```

!!! danger "Turn the security gate on"
    Every service ships with authentication **available but off**, preserving an open
    trusted-network posture. For any deployment reachable by others, set a key — it then fails
    closed on every route except health and docs:

    ```bash
    echo "HCLS_API_KEY=$(openssl rand -hex 24)" >> .env
    ```

Full cold-clone-to-running instructions, including seeding the agent corpora, are in the
[Build Guide](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/BUILD_GUIDE.md).

## What "live" means here

A capability marked `live` on the **[Capability Maturity Matrix](../honesty/maturity-matrix.md)** is
served by a real model against real input — never a mock. Where a capability is `planned`,
`preclinical`, or `gated`, it is labeled as such and not shown as if it ships. See the
[Honesty Ledger](../honesty/ledger.md) for what each badge means.

That rule is enforced, not aspirational: an internal audit found two capabilities registered `live`
with **nothing bound to their ports**. One was given a container; the other was demoted to `planned`
until an aggregator exists. Status follows evidence.

## Proof it's real

- **A real, passing gate on every merge** — 382 platform-library tests, plus real-bug lint and
  registry validation, with no `|| true` and no skipped gate. The 8,402 subject tests run in CI too.
- **The registry validates** — 42 registered capabilities, every engine and agent directory
  accounted for; a `live` capability can never be mock-served.
- **Verified on real public data** — the variant store on **GIAB HG002** (Ts/Tv ≈ 2.0), single-cell
  on **pbmc3k** (2,700 cells → 9 clusters, reproducible above), ACMG secondary-findings on real
  **ClinVar**, ESM-2 search on real proteins, ADMET on real drug molecules. Each is cited on
  [Citations & Evidence](../honesty/citations.md).
- **The site is honest by construction** — this page, the
  [Capability Maturity Matrix](../honesty/maturity-matrix.md) and every engine/agent page are
  generated from the registry, and the strict CI build fails before publishing anything stale.

*Next: the honest [hardware & elastic-burst story](hardware.md) · [reproducibility & datasets](reproducibility.md) · [what to install](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/ACQUISITION_MANIFEST.md).*
