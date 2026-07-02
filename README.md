<div align="center">

# HCLS AI Factory

**Patient DNA → therapeutic candidates, in hours, on a single box.**

Eight Engines · Eight Intelligence Agents · One Platform — open-source (Apache-2.0),
running end-to-end on one NVIDIA DGX Spark ($4,699). No cloud lock-in.

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
  vector DB, and developability scoring + a developability-guided design optimizer.
- **Single-cell** — real scanpy analysis (QC → clustering → DE → cell-type annotation), with the
  clinical agent reasoning *on top of* computed results.

## Eight intelligence agents

CAR-T · Precision Biomarker · Pharmacogenomics · Precision Autoimmune · Neurology ·
Rare-Disease Diagnostic · Single-Cell · Clinical-Trial — RAG clinical-decision-support over a
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
| | (agents 8521–8544) | | 8576 | Protein developability + design |

## Quickstart

```bash
# bring up the platform
docker compose -f docker-compose.dgx-spark.yml up -d

# drive it from an assistant (MCP): point your client at
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

## License

Apache-2.0. Built by Adam Jones.
