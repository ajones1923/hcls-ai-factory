<div align="center">

# HCLS AI Factory

**Patient DNA → therapeutic candidates, in hours, on a single box.**

Eight Engines · Eight Intelligence Agents · One Platform — open-source (Apache-2.0),
running end-to-end on one NVIDIA DGX Spark ($4,699). No cloud lock-in.

[![CI](https://github.com/ajones1923/hcls-ai-factory/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/ajones1923/hcls-ai-factory/actions/workflows/ci.yml)
[![tests](https://img.shields.io/badge/tests-8%2C100%2B%20passing-brightgreen)](https://github.com/ajones1923/hcls-ai-factory/actions/workflows/ci.yml)
[![license](https://img.shields.io/badge/license-Apache--2.0-blue)](LICENSE)
[![python](https://img.shields.io/badge/python-3.11%20%7C%203.12-blue)](https://github.com/ajones1923/hcls-ai-factory/actions/workflows/ci.yml)
[![docs](https://img.shields.io/badge/docs-hcls--ai--factory.org-0b7285)](https://hcls-ai-factory.org)

<br>

<img src="docs/brief/architecture.svg" alt="HCLS AI Factory architecture — 8 engines, 8 intelligence agents, and the flagship Tuberous Sclerosis Complex disease program on one platform" width="900">

📄 **[Read the Capability Brief →](docs/brief/)** — the whole factory in one place, in a **technical cut** for experts & builders and a **mission cut** for the broad audience.

</div>

---

## What it is

An end-to-end precision-medicine platform: real GPU genomics, a clinical RAG/intelligence
layer, real protein and small-molecule modeling, real single-cell analysis, and an
AI workflow composer that lets you describe an experiment in plain language and have it
run — every step on real models, all on one machine.

## The pipeline (3 core engines)

| Stage | Engine | What it does |
|---|---|---|
| 1 | **Genomic Foundation** | GPU variant calling (Parabricks `fq2bam` → DeepVariant **or** HaplotypeCaller), FASTQ → VCF, + a queryable variant store (VCF → DuckDB, Ts/Tv QC) |
| 2 | **Precision Intelligence** | Variant annotation (VEP / ClinVar / AlphaMissense) + RAG over a clinical knowledge base → druggable targets |
| 3 | **Therapeutic Discovery** | Target → ranked candidates: generation (MolMIM + BRICS), docking (DiffDock), and **real ADMET/tox** (104 endpoints) |

…plus five more engines — **Clinical Imaging** (4), **Precision Oncology** (5),
**Cardiology** (6), **Large-Molecule / Structural Biology** (7), and **Single-Cell Analysis** (8) —
and the **Tuberous Sclerosis** disease-program built on top of them. **Eight engines in all.**

## Engines 7 & 8 — real compute, not just retrieval

- **Proteins** — structure prediction (ESMFold), ESM-2 embeddings + similarity search over your
  vector DB, and ProteinMPNN sequence design — all `verified` in the registry. Developability
  scoring and the developability-guided design optimizer are written and tested (34 tests) but
  their service is not yet stood up, so the registry marks them **`planned`**; the aggregate
  Structural Biology engine endpoint (:8581) is `planned` for the same reason.
- **Single-cell** — real scanpy analysis (QC → clustering → DE → cell-type annotation), with the
  clinical agent reasoning *on top of* computed results.

## Eight intelligence agents

CAR-T · Precision Biomarker · Pharmacogenomics · Precision Autoimmune · Neurology ·
Rare-Disease · Single-Cell · Clinical-Trial — RAG clinical-decision-support over a
shared vector database.

## The platform layer

- **Capability Registry** — one manifest of every engine/agent/model/service with typed I/O.
- **Assistant Tool-Surface (MCP)** — drive the whole factory from Claude / Cursor / any MCP client.
- **AI Workflow Composer** — natural-language goal → a *validated, executable, governed* pipeline
  (shape-based wiring, self-repair, pre-run validation, AI root-cause).
- **Single-box MLOps** — experiment tracking + model registry + run lineage (SQLite, no warehouse).
- **Governance** — a clinical-overclaim honesty gate, a 21 CFR Part 11 reproducibility manifest,
  OpenTelemetry tracing. Mock outputs are always labeled; a `live` capability can never be mock-served.

## Service / port map

| Port | Service | | Port | Service |
|---|---|---|---|---|
| 5000 | Genomics | | 8570 | ESMFold (protein folding) |
| 5001 | Precision Intelligence (RAG) | | 8571 | Protein sequence search |
| 8505 | Therapeutic Discovery | | 8572 | ADMET / toxicity |
| 8080 | Landing dashboard | | 8573 | Single-cell compute |
| 19530 | Vector DB (Milvus) | | 8574 | Molecule generation |
| 3000 | Grafana · 9099 Prometheus | | 8575 | Variant store |
| | (agents 8126–8545) | | 8576 | Protein developability + design *(planned)* |

> **UIs are loopback-only.** The Streamlit interfaces bind to `127.0.0.1`; LAN access goes
> through Caddy with basic auth on 8721–8725. The agent APIs are unchanged and require
> `X-API-Key`. See `.env.example`.

## Quickstart

### Run it on a laptop — no GPU, no entitlement, no data download

**This is the honest entry point, and it is exactly what CI does on a stock x86 Ubuntu runner on
every pull request.** No DGX, no NVIDIA account, no 500 GB download:

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git && cd hcls-ai-factory
python -m venv .venv && .venv/bin/pip install --upgrade pip
.venv/bin/pip install -e lib/hcls_common
.venv/bin/pip install pytest pytest-asyncio fastapi "uvicorn[standard]" streamlit httpx \
    duckdb statsmodels biopython peft pymilvus loguru anthropic prometheus-client apscheduler \
    python-multipart tqdm python-dotenv flask flask-cors nibabel pydicom highdicom tenacity \
    reportlab python-docx lxml plotly scanpy anndata sentence-transformers

.venv/bin/python scripts/run_all_tests.py       # 17 suites, 8,191 tests
.venv/bin/python scripts/validate_registry.py   # the capability manifest
```

That exercises the platform layer, all eight engines, all eight agents and the disease program —
the logic, the schemas, the governance gates and the workflow composer. What it does **not** do is
call a GPU, a gated model, or a clinical corpus. Those need the real thing.

### What a full end-to-end run actually needs

You cannot `docker compose up` your way to a working factory, and this section exists so you find
that out here rather than three hours in:

| Requirement | Detail |
|---|---|
| **GPU box** | An NVIDIA DGX Spark, or an equivalent CUDA machine. The repo is developed on aarch64 (GB10, 128 GB unified memory); several NVIDIA artefacts are x86-only and burst to a remote host instead. |
| **NGC entitlement** | Parabricks (Stage 1 variant calling) and the BioNeMo NIMs are gated behind an NVIDIA account. |
| **Anthropic API key** | Every agent's clinical synthesis. Without it the agents still retrieve evidence but return no prose — deliberately, rather than a stub. |
| **~500 GB of data** | Stage 1 (FASTQ + GRCh38) has a working downloader. **Stages 2 and 3 are manual** — no automated downloader ships for ClinVar or AlphaMissense. |
| **Licences** | Several datasets are **non-commercial**, most notably AlphaMissense (CC BY-NC-SA 4.0). Read [`DATA_LICENSES.md`](DATA_LICENSES.md) before a commercial deployment. |

```bash
# Stage 1 data (automated, idempotent — skips what it already has)
cd core/engines/genomic-foundation
./run.sh check && ./run.sh login && ./run.sh download && ./run.sh reference

# Stages 2 & 3 are manual — section 6 of the DGX Spark Deployment Guide has the exact commands
# docs/HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md

# bring the platform up
cd - && docker compose -f docker-compose.dgx-spark.yml up -d
./start-factory.sh                              # non-Docker services, health-checked
```

### Drive it

```bash
# from an assistant (MCP): point your client at
python -m hcls_common.mcp_server        # tools: list/describe/health/invoke/plan/compose_workflow

# or compose a pipeline in code
python - <<'PY'
from hcls_common import WorkflowComposer, get_registry, FactoryTools
c = WorkflowComposer(get_registry(), tools=FactoryTools())
pipe, meta = c.compose("predict ADMET for a candidate molecule")
print(meta["checklist"]); print(c.run(pipe))
PY
```

## Repository map

```
core/
├── engines/                     # 8 engines — horizontal capabilities
│   ├── genomic-foundation/      #  1  GPU variant calling + variant store
│   ├── precision-intelligence/  #  2  annotation + clinical RAG
│   ├── therapeutic-discovery/   #  3  generation, docking, real ADMET
│   ├── clinical-imaging/        #  4  DICOM analysis (VISTA-3D / MAISI / VILA-M3)
│   ├── precision-oncology/      #  5  MTB packets, therapy ranking, trial matching
│   ├── cardiology/              #  6  clinical workflows + risk calculators
│   ├── structural-biology/      #  7  ESMFold, ESM-2 search, ProteinMPNN, developability
│   └── single-cell/             #  8  scanpy compute → cell-type annotation
├── agents/                      # 8 intelligence agents — clinical decision support
│   ├── cart · precision-biomarker · pharmacogenomics · precision-autoimmune
│   └── neurology · clinical-trial · rare-disease-diagnostic · single-cell
└── disease-programs/            # verticals composing the engines + agents
    └── tuberous-sclerosis/      #  first clinical beachhead
lib/hcls_common/                 # Shared platform: Capability Registry, MCP, Composer, MLOps, governance
hcls-orchestrator/               # Nextflow + the cross-stage trigger fabric
monitoring · docs · scripts · demo · data   ·   docker-compose.dgx-spark.yml · Caddyfile
```

> **Structure:** engines and agents live under `core/engines/` and `core/agents/`, disease
> verticals under `core/disease-programs/`, and the shared platform layer in `lib/hcls_common/`.
> See `docs/STRUCTURE.md` for the full layout.

## How to cite

If you build on this, please cite it. Machine-readable metadata is in
[`CITATION.cff`](CITATION.cff); GitHub renders a "Cite this repository" button from it.

```
Jones, A. (2026). HCLS AI Factory: an open precision-medicine platform —
eight compute engines, eight clinical intelligence agents, and the Tuberous
Sclerosis disease program. https://github.com/ajones1923/hcls-ai-factory
```

## License

Apache-2.0. Built by Adam Jones. Copyright 2026 Adam Jones.

**Building on this is welcome — that is the point of the licence.** Apache-2.0 §4 asks
three things of anyone redistributing this work or a derivative: keep the copyright notice,
include a readable copy of [`NOTICE`](NOTICE), and state what you changed. That is the
difference between building on the work and claiming it.

**The code is Apache-2.0; the data it reads is not.** No third-party dataset is redistributed here
— you download each from its source under that source's terms. Several are **non-commercial**,
most notably **AlphaMissense (CC BY-NC-SA 4.0)**, the most-referenced external artefact in this
codebase. Every dataset the platform reads, and where to check its terms:
[`DATA_LICENSES.md`](DATA_LICENSES.md).
