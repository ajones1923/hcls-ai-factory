# HCLS AI Factory — Repository Structure

The factory is **Eight Engines · Eight Intelligence Agents · One Platform**, plus disease-program
verticals. Engines and agents are horizontal capabilities under `core/`; the shared platform layer
is `lib/hcls_common/`.

```
core/
├── engines/                     # 8 engines — horizontal capabilities
│   ├── genomic-foundation/      #  1  GPU variant calling (Parabricks/DeepVariant) + variant store
│   ├── precision-intelligence/  #  2  annotation (ClinVar/AlphaMissense) + clinical RAG
│   ├── therapeutic-discovery/   #  3  MolMIM/BRICS generation, DiffDock docking, real ADMET
│   ├── clinical-imaging/        #  4  DICOM analysis (VISTA-3D / MAISI / VILA-M3)
│   ├── precision-oncology/      #  5  MTB packets, therapy ranking, trial matching
│   ├── cardiology/              #  6  clinical workflows + risk calculators
│   ├── structural-biology/      #  7  ESMFold, ESM-2 search, ProteinMPNN, developability
│   └── single-cell/             #  8  scanpy compute → cell-type annotation
├── agents/                      # 8 intelligence agents — clinical decision support
│   ├── cart/                    #     CAR-T cell therapy
│   ├── precision-biomarker/     #     biomarkers, PhenoAge/GrimAge
│   ├── pharmacogenomics/        #     star alleles, CPIC dosing
│   ├── precision-autoimmune/    #     autoantibody / HLA / flare
│   ├── neurology/               #     stroke triage, dementia, EDSS
│   ├── clinical-trial/          #     trial optimization + matching
│   ├── rare-disease-diagnostic/ #     HPO matching, ACMG classification
│   └── single-cell/             #     cell-type annotation, TME profiling (reasoning layer)
└── disease-programs/            # verticals composing the engines + agents for one condition
    └── tuberous-sclerosis/      #     first clinical beachhead
lib/hcls_common/                 # One Platform: capability registry, MCP, workflow composer,
                                 #                MLOps, governance gates, shared clients
hcls-orchestrator/               # Nextflow DSL2 + the cross-stage trigger fabric
docs/  monitoring/  scripts/  demo/  data/     docker-compose.dgx-spark.yml  Caddyfile
```

## Conventions

- **Directory names** are kebab-case; Python modules are snake_case.
- **`single-cell` appears in both `engines/` and `agents/`** by design — the *engine* is the scanpy
  compute layer (Engine 8), the *agent* is the reasoning layer on top. They are distinct
  capabilities, disambiguated by their registry IDs.
- **Two engines nest a deployable app under `agent/`** (`clinical-imaging/agent/`,
  `precision-oncology/agent/`) — historical, from when they were standalone apps. The engine
  directory still carries the engine-level README, docs, and (for oncology/imaging) the report
  tooling.
- Each engine/agent is self-contained: `src/`, `tests/`, `docs/`, `README.md`, and (for services)
  `api/`, `app/`, `config/`, `requirements.txt`, `Dockerfile`.

## Disease programs

Disease programs are the **vertical** solutions the horizontal engines + agents power. Each is a
self-contained folder under `core/disease-programs/` — its own engine, disease-specific agents,
config (which core capabilities it composes), and reference data — so it can be deployed or
replicated on its own. Tuberous Sclerosis Complex is the first; NF1/NF2/Rett/Williams and the
broader mTORopathies follow the same pattern.
