---
title: Software Acquisition Manifest
description: Every package, container, model weight, dataset and credential the HCLS AI Factory needs — what is already present, what is a plain install, and what is gated behind an account.
---

# Software Acquisition Manifest

**Date:** 2026-08-16 · **Target:** NVIDIA DGX Spark — GB10, 128 GB unified LPDDR5x, 20 Grace cores,
**aarch64** · **Method:** scanned 22 `requirements*.txt`, every import, every container reference and
every binary invocation in the tree, then checked each against what is actually installed.

**Architecture matters throughout.** Several NVIDIA artefacts are x86-only. Where that is true it is
stated, because it changes the plan from "install locally" to "elastic burst".

---

## 0. How to read this

| Tier | Meaning |
|---|---|
| **A — free, no account** | `pip install` or `apt`. Do these first; they unblock real work today. |
| **B — free account** | Registration only (HuggingFace licence acceptance, NCBI key). |
| **C — gated** | NGC entitlement, licence, or partner access. The long pole. |

Current state: **51 of 82** pinned Python packages are installed; **31 are missing**. Nine system
binaries are referenced and **seven are absent**. No gated container is present.

---

## 1. Tier A — free installs, no account needed

### 1.1 Python packages (31 missing)

```bash
cd /home/adam/projects/hcls-ai-factory
.venv/bin/pip install \
  pysam cyvcf2 \
  rdkit py3dmol stmol \
  monai monai-deploy-app-sdk SimpleITK pyradiomics scikit-image \
  e3nn einops wandb \
  openai langgraph nemoguardrails \
  minio redis psycopg2-binary nvidia-ml-py \
  opentelemetry-api opentelemetry-sdk \
  plotly pyvis PyPDF2 streamlit-chat imageio-ffmpeg decorator \
  fhir.resources
```

| Group | Packages | What breaks without them |
|---|---|---|
| Genomics | `pysam`, `cyvcf2` | BAM/VCF access; fast variant parsing |
| Cheminformatics | `rdkit`, `py3dmol`, `stmol` | Molecular properties, 3-D viewers in the UIs |
| Medical imaging | `monai`, `monai-deploy-app-sdk`, `SimpleITK`, `pyradiomics`, `scikit-image` | The 9 MONAI Application Packages |
| Structural ML | `e3nn`, `einops`, `wandb` | SE3 equivariance; experiment tracking |
| Agents / LLM | `openai`, `langgraph`, `nemoguardrails` | OpenAI-compatible + vLLM clients, guardrails |
| Infrastructure | `minio`, `redis`, `psycopg2-binary`, `nvidia-ml-py` | Object store, cache, Postgres, GPU telemetry |
| Observability | `opentelemetry-api`, `opentelemetry-sdk` | Tracing |
| Docs / UI | `plotly`, `pyvis`, `PyPDF2`, `streamlit-chat`, `imageio-ffmpeg`, `decorator` | Charts, graphs, PDF reports |
| Interop | `fhir.resources` | FHIR resource models |

> **`mkdocs-material` belongs in `.venv-docs`, never `.venv`.** Installing the docs toolchain into
> the platform venv pulls in `zarr` and `fast-array-utils`, whose pytest plugins abort test
> collection for three subjects. That cost 373 tests before it was found.

> **`nvidia-ml-py` provides the `pynvml` module** — the package and import names differ.

### 1.2 Build-from-source

```bash
.venv/bin/pip install git+https://github.com/NVIDIA/dllogger.git
```

### 1.3 System binaries (7 of 9 missing)

```bash
sudo apt-get update && sudo apt-get install -y \
  samtools bcftools tabix   # tabix provides bgzip
```

| Binary | Refs | Present | Purpose |
|---|---:|---|---|
| `samtools` | 8 | **no** | BAM/SAM manipulation |
| `bcftools` | 5 | **no** | VCF manipulation |
| `tabix` / `bgzip` | 5 | **no** | Indexed compressed genomic files |
| `nextflow` | 5 | **no** | **The orchestrator named throughout the architecture** |
| `gatk` | 2 | **no** | Broad toolkit (CPU fallback for variant calling) |
| `ffmpeg` | 1 | **no** | Media encoding |
| `pbrun` | 6 | **no** | Parabricks CLI — **Tier C** |
| `docker`, `node`, `curl`, `caddy`, `nvidia-smi`, `bwa` | — | yes | — |

**Nextflow** (needs a JRE):

```bash
sudo apt-get install -y openjdk-17-jre-headless
curl -s https://get.nextflow.io | bash && sudo mv nextflow /usr/local/bin/
nextflow -version
```

**GATK** — download the release zip from <https://github.com/broadinstitute/gatk/releases>, unzip,
and put `gatk` on `PATH`.

---

## 2. Tier B — free, account or licence acceptance

| # | Item | Where | Needed for |
|---|---|---|---|
| B1 | **BAAI/bge-small-en-v1.5** (+ `bge-large`) | HuggingFace — no gate | The embedding model behind every RAG agent. **242 references.** |
| B2 | **facebook/esmfold_v1** | HuggingFace — accept the licence | ESMFold structure prediction `:8570` |
| B3 | **NCBI E-utilities API key** | <https://ncbiinsights.ncbi.nlm.nih.gov> | PubMed ingest — raises the rate limit from 3/s to 10/s |
| B4 | **Anthropic API key** | <https://console.anthropic.com> | Claude reasoning in every agent. **Required for 11 of 17 demos.** |

```bash
.venv/bin/pip install sentence-transformers
.venv/bin/python -c "from sentence_transformers import SentenceTransformer as S; S('BAAI/bge-small-en-v1.5')"
export HF_TOKEN=<token>      # keep in .env, never commit
```

### Reference data (free, large)

| Dataset | Size | Source | Referenced |
|---|---|---|---|
| **GRCh38** reference FASTA + index | ~3 GB | <https://ftp.ncbi.nlm.nih.gov/genomes/> | 134× |
| **ClinVar** VCF (`clinvar.vcf.gz` + `.tbi`) | ~200 MB | <https://ftp.ncbi.nlm.nih.gov/pub/clinvar/> | 66× |
| **AlphaMissense** (`AlphaMissense_hg38.tsv.gz`) | ~600 MB | Zenodo / DeepMind — **CC BY-NC-SA 4.0, non-commercial** | **455×** |
| **HG002** FASTQ (Genome in a Bottle) | ~100 GB | <https://ftp-trace.ncbi.nlm.nih.gov/giab/> | 316× |

> **AlphaMissense is the most-referenced external artefact in the codebase and is licensed
> non-commercially.** Worth confirming that fits your distribution intent, since the platform is
> Apache-2.0 and public-facing.

> **HG002 is publicly consented** — the correct sample for a public demo. Never a patient sample.

`hcls-ai-factory-core-data/` (1.4 GB) already holds VCF, ClinVar, AlphaMissense extracts and PDB
structures. Check there before re-downloading.

---

## 3. Tier C — gated

This is the long pole. Nothing here is obtainable without an account and, in places, an entitlement.

### 3.1 NGC account — the prerequisite for everything below

```bash
# ngc.nvidia.com -> Setup -> API Key -> Generate
echo "NGC_API_KEY=<paste>" >> .env
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
```

**Checkpoint:** must report *Login Succeeded*. Nothing below works otherwise.

### 3.2 The gated inventory

| # | Component | Image / source | Blocks | Arch risk |
|---|---|---|---|---|
| **G1** | **CUDA PyTorch** for GB10 (sm_121) | `nvcr.io/nvidia/pytorch:<tag>-py3` | **Every ML path.** Currently `torch 2.10.0+cpu`, `cuda_is_available() = False` | prefer the container over pip wheels |
| **G2** | **Parabricks** | `nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1` | Secondary analysis (BWA-MEM2 + DeepVariant). **94 refs.** The "under five hours" claim | licensed |
| **G3** | **MolMIM NIM** | `nvcr.io/nim/nvidia/molmim:1.0.0` | Molecule generation. **131 refs** | **may be x86-only** |
| **G4** | **DiffDock NIM** | `nvcr.io/nim/mit/diffdock:2.2.0` | Docking poses. **121 refs** | **may be x86-only** |
| G5 | BioNeMo framework | `nvcr.io/nvidia/bionemo/bionemo-framework` | Model training/finetune | |
| G6 | **NV-Segment-CT NIM** | `nvcr.io/nim/nvidia/nv-segment-ct:latest` | Imaging segmentation — **:8534** | |
| G7 | NV-Reason-CXR-3B NIM | `nvcr.io/nim/nvidia/nv-reason-cxr-3b:latest` | Chest X-ray reasoning | |
| G8 | Llama-3-8B NIM | `nvcr.io/.../meta-llama3-8b` | Local LLM alternative to Claude | |
| G9 | **Chai-1** | Chai Discovery — gated weights | Co-folding. Registry status `planned` — **keep it there** | partner access |
| G10 | RFdiffusion + DGL | GitHub + gated weights; **needs G1** | De novo backbone design (vendored, not installed) | GPU-only |
| G11 | MHCflurry models | PyPI + model download | Immunogenicity `:8577` (`planned`) | |

### 3.3 The architecture question you must answer early

If **MolMIM and DiffDock are x86-only**, they cannot run on this box at all. Then:

- run them on remote GPUs over the private mesh, and
- **say "elastic burst" in every demo and script** — never "all on one box".

That single fact changes the narration of the flagship drug-discovery demo. Determine it **before**
building demos around it.

---

## 4. Credentials summary

| Credential | Tier | Store in | Used by |
|---|---|---|---|
| `ANTHROPIC_API_KEY` | B | `.env` | Every agent · **11 of 17 demos** |
| `HCLS_API_KEY` | — | `.env` | **The security gate. Ships OFF — set it.** |
| `HCLS_MINIO_USER` / `_PASSWORD` | — | `.env` | Milvus object store |
| `GRAFANA_PASSWORD` | — | `.env` | Monitoring |
| `NGC_API_KEY` | C | `.env` | All gated containers |
| `HF_TOKEN` | B | `.env` | ESMFold weights |
| `NCBI_API_KEY` | B | `.env` | PubMed ingest rate limit |

`.env` is gitignored (`.gitignore:76-77,102`). `.env.example` is the tracked template. **Never
commit a real key** — the audit found zero committed secrets and it should stay that way.

---

## 5. Recommended order

```
Tier A packages + binaries          ← today, no account, unblocks real work
        │
   .env + HCLS_API_KEY              ← the security gate ships off
        │
   compose up + seed corpora        ← 11 of 17 demos become genuinely LIVE
        │
   Tier B (BGE, Claude, datasets)   ← free accounts
        │
   NGC entitlement ──► G1 CUDA ──► G2 Parabricks
                          │
                          └──► G3/G4 NIMs ──► the flagship demo
```

**Do Tier A and the bring-up before touching Tier C.** When the gated software lands you want to
slot it into a platform that already runs, or you will be debugging integration and bring-up at the
same time and cannot tell which is failing.

---

## 6. Verification after each tier

```bash
.venv/bin/python scripts/run_all_tests.py        # 17 subjects · 8,402 tests
.venv/bin/python scripts/validate_registry.py    # port convention + supervisor cross-check
.venv/bin/python scripts/run_demo.py --check-all # which demos are actually ready
```

And for the GPU, do not accept `cuda.is_available() == True` as sufficient — confirm real work:

```bash
nvidia-smi --query-gpu=utilization.gpu,power.draw --format=csv -l 1
# run a workload alongside; utilisation must leave 0% and power must rise above the ~13-16 W idle floor
```

**Status follows evidence.** A capability is promoted to `live` only after it answers a health
probe; a demo is labelled LIVE only after it has run on real input. Two capabilities were found
registered `live` with nothing bound to their ports — that is the failure mode this rule prevents.
