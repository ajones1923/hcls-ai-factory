# HCLS AI Factory — Current State vs Required State

**Scope:** 8 Engines · 8 Intelligence Agents · 1 Disease Program (TSC) = 17 subjects.
**Method:** measured on the DGX Spark, 2026-08-15. Every number below was produced by running the
code, not by reading it. Where a first measurement was wrong, the correction is stated.

---

## 1. Headline

The codebase is **substantially more complete than it is operable**. 385,000 lines across 1,037
Python files and **8,402 passing tests** — but at the start of this session the platform's own merge
gate could not execute, and **not one engine or agent was running on the box**.

The gap is therefore not "write the features." It is **provisioning, deployment, and honesty of the
registry**.

| Dimension | State at start | State now | Remaining |
|---|---|---|---|
| Repo merge gate | **could not run** | passes | — |
| Subject test suites | unknown | **8,402 pass / 0 fail / 0 error** | — |
| Engines/agents running | **0 of 17** | 0 of 17 | deploy |
| Deployable via compose | 11 of 17 | 11 of 17 | **6 missing** |
| GPU usable from Python | **no** | no | CUDA torch |
| Registry accuracy | 2 claimed | 1 retracted | **1 open** |

---

## 2. What was broken and is now fixed

### 2.1 The merge gate could not run *(fixed)*

`CLAUDE.md` defines the gate as ruff + `pytest lib/hcls_common` + `validate_registry.py`. All three
failed with `No module named 'hcls_common'`.

Cause: a **stale editable install** in `~/.local/.../__editable___hcls_common_0_1_0_finder.py` whose
`MAPPING` dict was **empty** — the `.pth` loaded and resolved nothing. It could not be reinstalled
because PEP 668 marks system Python externally managed, and the repo had **no virtualenv**.

> I first reported `pyproject.toml` as misconfigured. That was wrong — `hcls_common = "."` is
> correctly set at line 83; my `sed` range had truncated the section. The stale artifact was the
> only cause.

**Fix:** created `.venv` at the repo root (`--system-site-packages`, so the heavy aarch64 wheels are
reused) and installed `lib/hcls_common` editable into it.

```
ruff --select E9,F82,F811,F706,F707 core lib scripts   All checks passed
validate_registry.py                                    42 capabilities — OK
pytest lib/hcls_common                                  372 passed
```

### 2.2 Neurology UI port *(changed)* — and a larger port defect it exposed *(OPEN)*

`config/settings.py` declared `STREAMLIT_PORT = 8534`, contradicting the agent's own README,
`neuro_ui.py`, `DESIGN_DOCUMENT.md`, `INDEX.md`, `DEMO_GUIDE.md`, its unit test and the capability
registry, which all say the UI is **8529**. 8534 is also a real service port —
`clinical-imaging/agent/docker-compose.yml:320` maps `8534:8000` for the **NV-Segment-CT NIM**.

Changed to 8529; two stale references in `docs/HCLS_AI_FACTORY_MINDMAP.md` corrected.

**But investigating it uncovered a systemic problem — see §3.7. Do not treat this as closed.**

### 2.3 Four missing Python dependencies *(fixed)*

`duckdb`, `statsmodels` (genomic-foundation variant store + GWAS), `biopython`, `peft`
(structural-biology). Ordinary PyPI packages, not gated. Installing them took two engines from
failing to green: genomic-foundation 142 → **156**, structural-biology 26 → **34**.

### 2.4 No unified test harness *(fixed)*

Every subject needed a different pytest invocation, so no single command reported platform health.
Added `scripts/run_all_tests.py`, which encodes the two traps that make a naive run wrong:

- **`src/collections.py` shadows the Python standard library** in 12 files across 11 subjects.
  Adding their `src/` to `PYTHONPATH` kills the interpreter *before collection*:
  `cannot import name 'namedtuple' from partially initialized module 'collections'`. The harness
  withholds `src/` for exactly those subjects.
- **`structural-biology/vendor_rfdiffusion/`** is vendored third-party code whose own tests need
  gated GPU packages. Excluded — it is not our code.

Both traps produced false readings during this analysis before being identified. See §5.

---

## 3. Current measured state

### 3.1 Test health — all 17 green

| Subject | Kind | Tests | Subject | Kind | Tests |
|---|---|---:|---|---|---:|
| cardiology | engine | 1,966 | pharmacogenomics | agent | 1,001 |
| clinical-imaging | engine | 1,365 | clinical-trial | agent | 769 |
| precision-oncology | engine | 556 | precision-biomarker | agent | 709 |
| precision-intelligence | engine | 157 | precision-autoimmune | agent | 455 |
| genomic-foundation | engine | 156 | cart | agent | 415 |
| therapeutic-discovery | engine | 124 | neurology | agent | 208 |
| structural-biology | engine | 34 | rare-disease-diagnostic | agent | 206 |
| single-cell | engine | 4 | single-cell | agent | 185 |
| tuberous-sclerosis | program | 92 | | | |

**Total 8,402 passed · 0 failed · 0 errors · 5 skipped.**

Test volume is wildly uneven — cardiology has 1,966 and the single-cell *engine* has 4. Coverage
depth, not pass rate, is the real signal.

### 3.2 Deployability — 6 of 17 cannot be started

`docker-compose.dgx-spark.yml` defines 17 services (plus 5 volumes). It covers **8 agents + 3
engines**. Missing entirely:

| Subject | Dockerfile | requirements | in compose |
|---|---|---|---|
| genomic-foundation | ✅ | ✅ | ❌ |
| precision-intelligence | ✅ | ✅ | ❌ |
| therapeutic-discovery | ✅ | ✅ | ❌ |
| **single-cell (engine)** | ❌ *(added)* | ❌ *(added)* | ❌ |
| **structural-biology** | ❌ **corrected** | ✅ | ❌ |
| tuberous-sclerosis | ✅ | ✅ | ❌ |

Five of these are one compose stanza away from running. **The flagship genomics engine and the
flagship TSC disease program are both un-startable by the platform's own launcher.**

The single-cell engine is the outlier: it is **real code** — a genuine scanpy pipeline (QC →
normalise → HVG → PCA → neighbours → Leiden → marker DE) with a canonical PBMC marker panel,
unit-tested pure logic, and a real 5.9 MB `pbmc3k_raw.h5ad` — but it has no Dockerfile, no
requirements, and no compose entry, so port 8573 is never bound.

> An earlier draft of this document called it a stub. It is not; it was checked before asserting.

### 3.3 Nothing is deployed

No engine or agent container exists on the box — not running, not even stopped. What *is* running is
the Marketing AI Factory (`intel-capture-*`), the cardiac-imaging Milvus stack, and monitoring
(Prometheus / Grafana / DCGM / node-exporter).

**Per-service virtualenvs exist for only 4 of ~20 supervised services.** `health-monitor.sh`
launches each service with `./venv/bin/python`; those venvs are present for
`precision-intelligence`, `therapeutic-discovery`, `tuberous-sclerosis` and
`genomic-foundation/web-portal`, and **absent for all 8 agents plus cardiology, imaging and
oncology**. `health-monitor.sh` does carry a `rebuild_venv()` path that can create them, so this is
self-healing in principle — it simply has not run. No container images have been built either.

> An earlier draft of this section said no per-service venvs existed at all. That was wrong: the
> search used `-maxdepth 3` and missed them. Four exist.

### 3.4 The GPU is unusable from Python

```
torch 2.10.0+cpu   torch.cuda.is_available() = False
```

On a DGX Spark whose entire value proposition is the GB10 and 128 GB of unified memory. Measured GPU
utilisation across this session: **0%, 13.7–15.1 W** (idle). Every PyTorch path — ESM-2 embeddings,
ESMFold, LoRA fine-tuning, any NIM run locally — is CPU-bound or non-functional.

This is the single highest-impact gap.

### 3.5 Dependency inventory

| Present | Missing |
|---|---|
| fastapi, streamlit, pymilvus, transformers, peft, biopython, duckdb, statsmodels | **CUDA torch**, scanpy, anndata, pysam, cyvcf2, rdkit, dgl, rfdiffusion, **nextflow** |

`nextflow` — the orchestrator named throughout the architecture — is **not installed**.

Code references by frequency: MolMIM 131 · DiffDock 121 · Parabricks 94 · BioNeMo 65 ·
DeepVariant 33 · ESMFold 28 · AlphaFold 7 · Chai-1 1.

### 3.6 Registry accuracy — 1 open defect (1 claim retracted)

`lib/hcls_common/capabilities.json`: 42 capabilities — 9 engine, 8 agent, 9 service, 10 model,
3 NIM, 3 stage. 35 `live`, 7 `planned`.

1. ~~`tuberous-sclerosis-engine` is typed `engine`, giving 9 engines.~~ **Not a defect — verified.**
   TSC is deliberately `type: engine` with `tags: ["engine", "tsc", "disease-program"]`, and
   `scripts/site_gen_matrix.py::_is_program()` filters it out of the engine group:
   *"the disease-program vertical (TSC) is tagged, not counted among the 8 engines."* The public
   site therefore shows 8 engines correctly. The `9 typed 'engine'` line printed by
   `validate_registry.py` is a raw type count, not a roster claim. **Do not "fix" this** — retyping
   TSC would drop it out of every `types=("engine","agent")` query, including the port drift-guard.

2. **`singlecell-compute` is marked `live` with `endpoint: localhost:8573`** but nothing can bind
   that port — no Dockerfile, no requirements, no compose entry. `CLAUDE.md` requires that a `live`
   capability is never mock-served; an unreachable endpoint is weaker still. **This is the one real
   registry defect.** The right fix is to make it genuinely live (the code is real), not to
   downgrade the label.

### 3.7 Port map: registry vs live supervisor — 8 of 13 disagreed *(RESOLVED 2026-08-15)*

The capability registry and `health-monitor.sh` (the cron-driven supervisor that is the **actual**
deployment mechanism) disagree on 8 of 13 subjects:

| Capability | registry | health-monitor | reading |
|---|---:|---:|---|
| precision-oncology-agent | 8526 | 8527 | UI / API |
| cart-intelligence-agent | 8521 | 8522 | UI / API |
| precision-biomarker-agent | 8528 | 8529 | UI / API |
| precision-autoimmune-agent | 8531 | 8532 | UI / API |
| neurology-intelligence-agent | 8529 | 8528 | **UI / API — inverted** |
| imaging-intelligence-agent | 8525 | 8524 | **UI / API — inverted** |
| rare-disease-diagnostic-agent | 8544 | 8134 | different scheme |
| pharmacogenomics-intelligence-agent | 8507 | 8107 | different scheme |

Most of this is **not** drift: the registry advertises the **UI** port and health-monitor supervises
the **API**, which is `UI + 1` for oncology, cart, biomarker and autoimmune. That convention is fine
but undocumented.

Two subjects invert it, and that produces a **real double collision**:

```
precision-biomarker :  UI 8528 (registry)   API 8529 (supervisor)
neurology           : API 8528 (supervisor)  UI 8529 (registry)
```

Both agents claim **both** ports. Neither can run alongside the other. The registry's drift-guard
does not catch it because it only compares registry entries with each other — the supervisor's ports
are invisible to it.

Separately, `health-monitor.sh` launches the neurology UI on **8534**, which belongs to the
NV-Segment-CT NIM.

**Resolved — Adam adopted *registry = UI, API = UI + 1* on 2026-08-15.** Three services were
re-seated (imaging 8525→8523, neurology 8529→8535, structural-biology 8579→8581) and the convention
is now gated in `validate_registry.py`. Canonical allocation: [`PORT_MAP.md`](PORT_MAP.md).

Original recommendation, for the record:
1. Adopt one documented convention: *registry endpoint = UI; API = UI + 1*.
2. Re-seat neurology and imaging onto it (neurology API 8530, imaging API 8526 — check availability).
3. Move the neurology UI off 8534 to its registry port 8529.
4. Extend `validate_registry.py` to parse `health-monitor.sh` so this class of conflict is gated.

~~Until 1–3 are decided, neurology and precision-biomarker cannot both be brought up.~~ Done.

---

## 4. Prioritised gap list

| # | Gap | Impact | Effort |
|---|---|---|---|
| G1 | CUDA-enabled torch for GB10/aarch64 | Unlocks every ML path | gated — NVIDIA channel |
| G2 | 5 engines + TSC absent from compose | Half the platform unreachable | low |
| G3 | single-cell engine has no container at all | Registered `live`, unreachable | low |
| G4 | Nothing deployed; no images built | Platform is cold | medium |
| G5 | nextflow not installed | Orchestrator missing | low |
| G7 | `singlecell-compute` mislabelled `live` | Honesty gate | trivial |
| G8 | scanpy/rdkit/pysam/cyvcf2 absent | Engine runtimes incomplete | low |
| G9 | `src/collections.py` shadows stdlib (12 files) | Latent catastrophic import failure | low |
| ~~G12~~ | ~~Neurology/biomarker double port collision~~ | **CLOSED** — convention adopted + gated | done |
| G13 | structural-biology has no Dockerfile of its own | Cannot be containerised | low |
| G10 | Test depth uneven (4 → 1,966) | Confidence uneven | high |
| G11 | BioNeMo NIMs / Parabricks / Chai-1 not provisioned | Demos cannot run end to end | gated |

G1 and G11 are **gated** — they need NGC credentials and licensed containers. They are covered in
the separate gated-software guide.

---

## 5. Measurement errors made during this analysis

Recorded because they would otherwise recur, and because two of them produced confidently wrong
numbers that looked plausible.

1. **First subject survey reported 6.7M LOC** for precision-intelligence and 9.4M for
   therapeutic-discovery — it was counting vendored virtualenvs. Real totals are 11,346 and 9,227.
2. **First test run reported all 17 suites as "collection failed."** The harness passed
   `--timeout=120` without `pytest-timeout` installed, so pytest exited on a usage error before
   running anything. cart alone actually has 415 passing tests.
3. **Second test run reported the 8 agents as "no tests collected."** The harness put `src/` on
   `PYTHONPATH`, hitting the stdlib-shadowing trap in §2.4.
4. **`pyproject.toml` was reported broken** — it was not; a truncated `sed` range hid the correct
   line.
5. **TSC's registry type was reported as a roster defect** — it is a deliberate, working design
   (typed `engine`, tagged `disease-program`, filtered by the site generator). Checked before
   changing it; nothing was changed.

The lesson for the PRD: **a uniform result across every subject is a signal that the harness is
wrong, not that the platform is uniformly broken.**
