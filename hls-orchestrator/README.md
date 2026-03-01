# HLS Pipeline Orchestrator

Nextflow DSL2 orchestrator for the HCLS AI Factory. Coordinates the three-stage precision medicine pipeline: Genomics, RAG/Chat Evidence, and Drug Discovery.

## Architecture

```
                         ┌─────────────────────┐
                         │    ORCHESTRATOR      │
                         │   (main.nf / DSL2)   │
                         └──────────┬──────────┘
                                    │
              ┌─────────────────────┼─────────────────────┐
              │                     │                     │
              ▼                     ▼                     ▼
     ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
     │    GENOMICS      │  │    RAG/CHAT     │  │  DRUG DISCOVERY │
     │    MODULE        │  │    MODULE       │  │    MODULE       │
     │                  │  │                  │  │                  │
     │  FASTQ → VCF    │  │  VCF → Target   │  │  Target → Drugs │
     │  (Parabricks)   │  │  (Milvus+Claude)│  │  (BioNeMo NIMs) │
     └─────────────────┘  └─────────────────┘  └─────────────────┘
```

## Pipeline Modes

| Mode | Command | Stages | Use Case |
|---|---|---|---|
| **Full** | `--mode full` | 1 → 2 → 3 | Complete DNA-to-Drug pipeline |
| **Target** | `--mode target` | 2 → 3 | Start from existing VCF (skip genomics) |
| **Drug** | `--mode drug` | 3 | Start from validated target (skip genomics + RAG) |
| **Demo** | `--mode demo` | 3 | VCP/FTD demonstration with pre-loaded target |
| **Genomics Only** | `--mode genomics_only` | 1 | Variant calling only |

## Quick Start

### Run the Demo Pipeline

```bash
cd hls-orchestrator

# Using Nextflow
./nextflow run main.nf -profile dgx_spark --mode demo

# Using Python alternative (bypasses cgroup issues on some systems)
python3 run_pipeline.py --mode demo
```

### Run the Full Pipeline

```bash
./nextflow run main.nf \
  -profile dgx_spark \
  --mode full \
  --input /path/to/fastq_pairs.csv \
  --outdir results/
```

### Run from Existing VCF

```bash
./nextflow run main.nf \
  -profile dgx_spark \
  --mode target \
  --vcf /path/to/variants.vcf \
  --outdir results/
```

## Execution Profiles

| Profile | Executor | Resources | Use Case |
|---|---|---|---|
| `standard` | Local | Default limits | Development and testing |
| `docker` | Local + Docker | NVIDIA runtime | Standard deployment |
| `singularity` | Local + Singularity | GPU support | HPC without Docker |
| `dgx_spark` | Local | GPU + high memory | NVIDIA DGX Spark |
| `slurm` | SLURM | Cluster scheduling | HPC cluster |
| `test` | Local | Reduced resources | CI/CD and quick verification |

## Configuration

Key parameters in `nextflow.config`:

```groovy
params {
    // Pipeline mode
    mode = 'demo'

    // Stage 1: Genomics
    genome = 'GRCh38'

    // Stage 2: RAG/Chat
    llm_model = 'claude-sonnet-4-20250514'
    max_targets = 5
    confidence_threshold = 0.7

    // Stage 3: Drug Discovery
    num_molecules = 20
    diversity = 0.3
    max_mw = 550
    docking_poses = 10

    // NIM services
    molmim_url = 'http://localhost:8001'
    diffdock_url = 'http://localhost:8002'

    // Resources
    max_memory = '128.GB'
    max_cpus = 32
    max_gpus = 1
}
```

## Modules

| Module | File | Input | Output |
|---|---|---|---|
| Genomics | `modules/genomics.nf` | FASTQ pairs | BAM + VCF |
| RAG/Chat | `modules/rag_chat.nf` | VCF | Target hypotheses (JSON) |
| Drug Discovery | `modules/drug_discovery.nf` | Target JSON | Ranked molecules |
| Reporting | `modules/reporting.nf` | All outputs | HTML + JSON reports |

## Demo Target

The demo pipeline uses a pre-configured target for Frontotemporal Dementia:

| Field | Value |
|---|---|
| Gene | VCP (Valosin-containing protein, p97) |
| UniProt | P55072 |
| Disease | Frontotemporal Dementia, ALS, IBM |
| Mechanism | AAA+ ATPase inhibition |
| Seed Compound | CB-5083 (ATP-competitive inhibitor) |
| PDB Structures | 5FTK, 8OOI, 9DIL, 7K56 |

## Portal Dashboard

A Streamlit-based dashboard for monitoring pipeline execution:

```bash
# Start the portal on port 8510
./run_portal.sh

# Access at http://localhost:8510
```

## Outputs

Pipeline outputs are written to timestamped directories under `results/`:

```
results/demo_20260228_143000/
├── report.json     # Structured results (molecules, scores, metadata)
└── report.html     # Visual report with ranking tables
```

## Resume and Fault Tolerance

Nextflow provides built-in checkpoint/resume:

```bash
# Resume from last successful stage after a failure
./nextflow run main.nf -resume
```

Process error handling is configured with automatic retry on transient failures (out-of-memory, timeout).

## Container Images

| Process | Image |
|---|---|
| Genomics | `nfcore/sarek:3.4.0` |
| RAG/Chat | `hls-pipeline/rag-chat:latest` |
| Drug Discovery | `hls-pipeline/drug-discovery:latest` |
| MolMIM | `nvcr.io/nim/nvidia/molmim:1.0.0` |
| DiffDock | `nvcr.io/nim/mit/diffdock:2.2.0` |

## Python Alternative

For environments where Nextflow encounters cgroup issues (common on some ARM64 configurations), use the Python orchestrator:

```bash
python3 run_pipeline.py --mode demo
```

This provides the same pipeline stages with direct Python process management instead of Nextflow.
