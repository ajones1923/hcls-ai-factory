# HCLS AI Factory — Project Bible

> **Purpose:** Complete implementation reference for building the HCLS AI Factory on NVIDIA DGX Spark. This platform transforms patient DNA into novel drug candidates in under 5 hours across three GPU-accelerated stages: Secondary Genomics, RAG-Grounded Target Identification, and AI-Driven Drug Discovery. Import this document into a Claude Code session as context for implementation.
>
> **License:** Apache 2.0 | **Author:** Adam Jones | **Date:** February 2026

---

## Table of Contents

1. [Project Overview & Goals](#1-project-overview--goals)
2. [DGX Spark Hardware Reference](#2-dgx-spark-hardware-reference)
3. [Repository Layout](#3-repository-layout)
4. [Docker Compose Services](#4-docker-compose-services)
5. [Stage 1: Genomics Pipeline](#5-stage-1-genomics-pipeline)
6. [Stage 2: RAG/Chat Pipeline](#6-stage-2-ragchat-pipeline)
7. [Milvus Vector Database Schema](#7-milvus-vector-database-schema)
8. [Variant Annotation Pipeline](#8-variant-annotation-pipeline)
9. [Knowledge Base — 201 Genes, 13 Therapeutic Areas](#9-knowledge-base--201-genes-13-therapeutic-areas)
10. [Anthropic Claude LLM Integration](#10-anthropic-claude-llm-integration)
11. [Stage 3: Drug Discovery Pipeline](#11-stage-3-drug-discovery-pipeline)
12. [BioNeMo NIM Services](#12-bionemo-nim-services)
13. [Drug-Likeness Scoring](#13-drug-likeness-scoring)
14. [Cryo-EM Structure Evidence](#14-cryo-em-structure-evidence)
15. [VCP/FTD Demo Walkthrough](#15-vcpftd-demo-walkthrough)
16. [Pydantic Data Models](#16-pydantic-data-models)
17. [Nextflow DSL2 Orchestration](#17-nextflow-dsl2-orchestration)
18. [Landing Page & Service Health](#18-landing-page--service-health)
19. [Monitoring Stack](#19-monitoring-stack)
20. [Cross-Modal Integration](#20-cross-modal-integration)
21. [Configuration Reference](#21-configuration-reference)
22. [Deployment Roadmap](#22-deployment-roadmap)
23. [Testing Strategy](#23-testing-strategy)
24. [Implementation Sequence](#24-implementation-sequence)

---

## 1. Project Overview & Goals

### What This Platform Does

The HCLS AI Factory is an end-to-end precision medicine platform that takes a patient's raw DNA sequencing data (FASTQ) and produces ranked novel drug candidates — all on a single NVIDIA DGX Spark desktop workstation. Three GPU-accelerated stages execute sequentially: variant calling, RAG-grounded target identification, and generative drug discovery.

### Three-Stage Pipeline

| Stage | Function | Duration | Key Output |
|---|---|---|---|
| 1 — Genomics | BWA-MEM2 alignment + DeepVariant calling | 120-240 min | VCF (~11.7M variants) |
| 2 — RAG/Chat | Annotation → Embedding → LLM reasoning | Interactive | Target gene + evidence |
| 3 — Drug Discovery | MolMIM generation → DiffDock docking → RDKit scoring | 8-16 min | 100 ranked drug candidates |

### End-to-End Flow

```
Patient DNA → Illumina Sequencer → FASTQ (~200 GB)
  → Parabricks fq2bam → BAM
  → DeepVariant → VCF (11.7M variants)
  → ClinVar + AlphaMissense + VEP annotation
  → Milvus vector indexing (3.5M embeddings)
  → Claude RAG reasoning → Target hypothesis (gene + evidence)
  → RCSB PDB structure retrieval
  → MolMIM molecule generation
  → DiffDock molecular docking
  → RDKit drug-likeness scoring
  → 100 ranked novel drug candidates + PDF report
```

### Design Principles

1. **GPU-first:** Every compute-intensive step runs on the GB10 GPU
2. **Clinically grounded:** ClinVar, AlphaMissense, and VEP provide evidence-based annotation
3. **Reproducible:** Nextflow DSL2 orchestration with containerized processes
4. **Open:** Apache 2.0 license, open-source tools, public reference databases
5. **Desktop-scale:** Runs entirely on a $3,999 DGX Spark

---

## 2. DGX Spark Hardware Reference

### Specifications

| Parameter | Value |
|---|---|
| CPU | NVIDIA Grace (ARM64 / aarch64), 144 cores |
| GPU | NVIDIA GB10, 1 GPU |
| Memory | 128 GB unified LPDDR5x (CPU + GPU shared pool) |
| System RAM | 512 GB |
| Storage | NVMe, high-throughput I/O |
| Storage Access | GPUDirect Storage (zero-copy GPU access) |
| Price | $3,999 |
| OS | Ubuntu-based (NVIDIA DGX OS) |

### Critical: ARM64 Architecture

**ALL containers must be ARM64-compatible.** The Grace CPU is aarch64, not x86_64. This affects:
- Base Docker images (must use ARM64 variants)
- Python wheel availability (most scientific packages have ARM64 wheels)
- NVIDIA container images (use NGC ARM64 variants)
- Any compiled C/C++ extensions (RDKit, BioPython)

### Unified Memory Model

The 128 GB LPDDR5x is **shared** between CPU and GPU — there is no separate GPU VRAM. This means:
- No explicit CPU→GPU data transfers needed for many operations
- Memory pressure from CPU workloads reduces GPU-available memory
- Monitor total system memory, not just "GPU memory"
- Parabricks fq2bam peaks at ~40 GB, DeepVariant peaks at ~60 GB

### Storage Requirements

| Dataset | Size | Notes |
|---|---|---|
| GRCh38 reference | 3.1 GB | Pre-indexed for BWA-MEM2 |
| FASTQ input (30× WGS) | ~200 GB | HG002 paired-end |
| BAM intermediate | ~100 GB | Temporary, deleted after VCF |
| ClinVar database | ~1.2 GB | 4.1M clinical variants |
| AlphaMissense database | ~4 GB | 71M predictions |
| Milvus index | ~2 GB | 3.5M × 384-dim vectors |
| BioNeMo model cache | ~10 GB | MolMIM + DiffDock weights |
| **Total minimum** | **~320 GB** | Plus OS and Docker layers |

---

## 3. Repository Layout

```
hcls-ai-factory-public/
├── README.md                           # Project overview
├── LICENSE                             # Apache 2.0
├── docker-compose.yml                  # All services
├── start-services.sh                   # Service startup orchestration
├── .env.example                        # Environment variable template
│
├── hls-orchestrator/                   # Nextflow pipeline orchestration
│   ├── main.nf                         # DSL2 entry point
│   ├── nextflow.config                 # Profiles and parameters
│   ├── run_pipeline.py                 # Python CLI launcher
│   ├── modules/
│   │   ├── genomics.nf                 # Stage 1 processes
│   │   ├── rag_chat.nf                 # Stage 2 processes
│   │   ├── drug_discovery.nf           # Stage 3 processes
│   │   └── reporting.nf                # Report generation
│   └── tests/
│
├── genomics-pipeline/                  # Stage 1: Parabricks
│   ├── README.md                       # Genomics documentation (48 KB)
│   ├── Dockerfile
│   ├── src/
│   │   ├── run_parabricks.py           # fq2bam + DeepVariant launcher
│   │   ├── vcf_stats.py                # VCF quality statistics
│   │   └── web_portal.py               # Flask portal (:5000)
│   ├── config/
│   │   └── parabricks.yaml             # GPU resource allocation
│   └── tests/
│
├── rag-chat-pipeline/                  # Stage 2: RAG + Claude
│   ├── README.md                       # RAG documentation (51 KB)
│   ├── Dockerfile
│   ├── src/
│   │   ├── rag_engine.py               # Core RAG orchestration (23 KB)
│   │   ├── milvus_client.py            # Milvus vector DB client (13 KB)
│   │   ├── annotator.py                # ClinVar + AlphaMissense + VEP (23 KB)
│   │   ├── knowledge.py                # 201 genes, 13 areas (88 KB)
│   │   ├── streamlit_chat.py           # Chat UI (:8501)
│   │   └── api.py                      # REST API (:5001)
│   ├── config/
│   │   └── milvus.yaml                 # Vector DB configuration
│   └── tests/
│
├── drug-discovery-pipeline/            # Stage 3: BioNeMo + RDKit
│   ├── README.md                       # Drug discovery documentation (56 KB)
│   ├── Dockerfile
│   ├── src/
│   │   ├── pipeline.py                 # 10-stage orchestration (18 KB)
│   │   ├── nim_clients.py              # MolMIM + DiffDock clients (15 KB)
│   │   ├── molecule_generator.py       # SMILES generation (11 KB)
│   │   ├── cryoem_evidence.py          # Cryo-EM structure scoring (6 KB)
│   │   ├── models.py                   # Pydantic data models (8 KB)
│   │   ├── streamlit_discovery.py      # Discovery UI (:8505)
│   │   └── portal.py                   # Discovery portal (:8510)
│   ├── config/
│   │   └── discovery.yaml              # Pipeline parameters
│   └── tests/
│
├── landing-page/                       # HCLS AI Factory entry point
│   ├── Dockerfile
│   └── src/
│       └── landing.py                  # Flask landing page (:8080)
│
├── monitoring/                         # Observability stack
│   ├── prometheus.yml                  # Scrape configuration
│   └── grafana/
│       └── dashboards/
│           └── hcls-factory.json       # GPU + pipeline dashboard
│
└── docs/                               # Documentation
    ├── PRODUCT_DOCUMENTATION.txt       # Full product docs (122 KB)
    ├── ARCHITECTURE_MINDMAP.md         # Architecture reference
    └── PIPELINE_REPORT.md              # Pipeline analysis (29 KB)
```

---

## 4. Docker Compose Services

### Port Allocation

| Service | Port | Protocol | Stage |
|---|---|---|---|
| Landing Page | 8080 | HTTP (Flask) | Orchestration |
| Genomics Portal | 5000 | HTTP (Flask) | Stage 1 |
| RAG REST API | 5001 | HTTP REST | Stage 2 |
| Milvus Vector DB | 19530 | gRPC | Stage 2 |
| Attu (Milvus UI) | 8000 | HTTP | Stage 2 |
| Streamlit Chat | 8501 | HTTP | Stage 2 |
| MolMIM NIM | 8001 | HTTP REST | Stage 3 |
| DiffDock NIM | 8002 | HTTP REST | Stage 3 |
| Discovery UI | 8505 | HTTP (Streamlit) | Stage 3 |
| Discovery Portal | 8510 | HTTP | Stage 3 |
| Grafana | 3000 | HTTP | Monitoring |
| Prometheus | 9099 | HTTP | Monitoring |
| Node Exporter | 9100 | HTTP | Monitoring |
| DCGM Exporter | 9400 | HTTP | Monitoring |

### Key Container Images

| Service | Image | Notes |
|---|---|---|
| Parabricks | `nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1` | GPU-accelerated genomics |
| Milvus | `milvusdb/milvus:v2.4-latest` | Vector database |
| MolMIM | `nvcr.io/nvidia/clara/bionemo-molmim:1.0` | Molecule generation NIM |
| DiffDock | `nvcr.io/nvidia/clara/diffdock:1.0` | Molecular docking NIM |
| Grafana | `grafana/grafana:10.2.2` | Dashboards |
| Prometheus | `prom/prometheus:v2.48.0` | Metrics TSDB |

### Service Startup Order

The `start-services.sh` script orchestrates startup in dependency order:

```bash
# 1. Infrastructure (Milvus, monitoring)
# 2. Stage 1 services (Parabricks, genomics portal)
# 3. Stage 2 services (RAG engine, Streamlit chat)
# 4. Stage 3 services (BioNeMo NIMs, discovery UI)
# 5. Landing page (health monitor for all 10 services)
```

### Health Monitoring

The landing page at port 8080 monitors 10 services with periodic health checks:

| Service | Health Endpoint | Check Interval |
|---|---|---|
| Parabricks | Port 5000 `/health` | 30s |
| Milvus | Port 19530 gRPC ping | 30s |
| RAG API | Port 5001 `/health` | 30s |
| Chat UI | Port 8501 `/healthz` | 30s |
| MolMIM NIM | Port 8001 `/v1/health/ready` | 30s |
| DiffDock NIM | Port 8002 `/v1/health/ready` | 30s |
| Discovery UI | Port 8505 `/healthz` | 30s |
| Grafana | Port 3000 `/api/health` | 30s |
| Prometheus | Port 9099 `/-/healthy` | 30s |
| DCGM Exporter | Port 9400 `/metrics` | 30s |

---

## 5. Stage 1: Genomics Pipeline

### Overview

Stage 1 takes raw FASTQ files from an Illumina sequencer and produces a Variant Call Format (VCF) file using NVIDIA Parabricks — a GPU-accelerated implementation of industry-standard bioinformatics tools.

### Input Specifications

| Parameter | Value |
|---|---|
| Sample | HG002 (GIAB reference standard) |
| Coverage | 30× whole-genome sequencing (WGS) |
| Read Length | 2×250 bp paired-end |
| File Size | ~200 GB (FASTQ pair) |
| Reference Genome | GRCh38 (3.1 GB, pre-indexed) |
| Format | FASTQ (gzip-compressed) |

### Pipeline Steps

#### Step 1: BWA-MEM2 Alignment (`fq2bam`)

```bash
pbrun fq2bam \
  --ref /reference/GRCh38.fa \
  --in-fq /data/HG002_R1.fastq.gz /data/HG002_R2.fastq.gz \
  --out-bam /output/HG002.bam \
  --num-gpus 1
```

| Metric | Value |
|---|---|
| Duration | 20-45 minutes |
| GPU Utilization | 70-90% |
| Peak Memory | ~40 GB |
| Output | Sorted BAM + BAI index |
| Algorithm | BWA-MEM2 (GPU-accelerated) |

#### Step 2: DeepVariant Variant Calling

```bash
pbrun deepvariant \
  --ref /reference/GRCh38.fa \
  --in-bam /output/HG002.bam \
  --out-variants /output/HG002.vcf.gz \
  --num-gpus 1
```

| Metric | Value |
|---|---|
| Duration | 10-35 minutes |
| GPU Utilization | 80-95% |
| Peak Memory | ~60 GB |
| Output | VCF (gzip-compressed + tabix index) |
| Algorithm | Google DeepVariant (CNN-based, >99% accuracy) |

### Output: VCF Statistics

| Metric | Count |
|---|---|
| Total Variants | ~11.7M |
| High-Quality (QUAL>30) | ~3.5M |
| SNPs | ~4.2M |
| Indels | ~1.0M |
| Coding Region Variants | ~35,000 |
| Multi-allelic Sites | ~150,000 |

### Parabricks Container

```
Image: nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
GPU: Required (CUDA)
Volumes: /reference, /data, /output
Port: 5000 (Flask web portal for run status)
```

### Genomics Portal (Port 5000)

The Flask-based portal provides:
- Real-time pipeline progress monitoring
- BAM quality statistics (mapping rate, duplication, coverage)
- VCF summary statistics and variant type distributions
- Run history and configuration logs

---

## 6. Stage 2: RAG/Chat Pipeline

### Overview

Stage 2 annotates the VCF variants with clinical and functional databases, indexes them in a Milvus vector database, and uses Anthropic Claude with RAG to identify druggable gene targets supported by evidence.

### Architecture

```
VCF (11.7M variants)
  → Quality filter (QUAL>30) → 3.5M variants
  → ClinVar annotation → clinical significance
  → AlphaMissense annotation → pathogenicity prediction
  → VEP annotation → functional consequences
  → BGE-small-en-v1.5 embedding → 384-dim vectors
  → Milvus IVF_FLAT indexing → 3.5M searchable embeddings
  → Claude RAG query → target hypothesis with evidence chain
```

### Annotation Funnel

| Stage | Variant Count | Filter |
|---|---|---|
| Raw VCF | ~11.7M | — |
| Quality filter | ~3.5M | QUAL > 30 |
| ClinVar match | ~35,616 | Clinical significance annotated |
| AlphaMissense match | ~6,831 | AI pathogenicity predicted |
| Coding + pathogenic | ~2,400 | Actionable subset |

### Embedding Model

| Parameter | Value |
|---|---|
| Model | BGE-small-en-v1.5 |
| Dimensions | 384 |
| Input | Text summary of annotated variant |
| Index Type | IVF_FLAT |
| Index Params | nlist=1024 |
| Search Params | nprobe=16 |
| Distance Metric | COSINE |
| Total Embeddings | ~3.5M |

### Query Flow

1. User asks a natural language question in the Streamlit chat
2. Query is expanded using 10 therapeutic area keyword maps
3. BGE-small-en-v1.5 embeds the expanded query
4. Milvus performs approximate nearest-neighbor search (top_k=20)
5. Retrieved variant contexts are assembled into a RAG prompt
6. Claude processes the prompt with knowledge base grounding
7. Response includes gene target, evidence chain, and confidence assessment

---

## 7. Milvus Vector Database Schema

### Collection: `genomic_evidence`

| Field | Type | Description |
|---|---|---|
| `id` | INT64 (PK, auto) | Primary key |
| `embedding` | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| `chrom` | VARCHAR(10) | Chromosome (chr1-chr22, chrX, chrY) |
| `pos` | INT64 | Genomic position |
| `ref` | VARCHAR(1000) | Reference allele |
| `alt` | VARCHAR(1000) | Alternate allele |
| `qual` | FLOAT | Variant quality score |
| `gene` | VARCHAR(100) | Gene symbol (e.g., VCP, EGFR) |
| `consequence` | VARCHAR(200) | Functional consequence (e.g., missense_variant) |
| `impact` | VARCHAR(20) | Impact level (HIGH, MODERATE, LOW, MODIFIER) |
| `genotype` | VARCHAR(10) | Sample genotype (0/1, 1/1, etc.) |
| `text_summary` | VARCHAR(2000) | Human-readable variant description |
| `clinical_significance` | VARCHAR(200) | ClinVar classification |
| `rsid` | VARCHAR(20) | dbSNP identifier (e.g., rs188935092) |
| `disease_associations` | VARCHAR(2000) | Associated diseases/conditions |
| `am_pathogenicity` | FLOAT | AlphaMissense pathogenicity score (0-1) |
| `am_class` | VARCHAR(20) | AlphaMissense class (pathogenic/ambiguous/benign) |

**Total: 17 fields**

### Index Configuration

```python
index_params = {
    "index_type": "IVF_FLAT",
    "metric_type": "COSINE",
    "params": {"nlist": 1024}
}

search_params = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16}
}
```

### Milvus Infrastructure

| Component | Port | Purpose |
|---|---|---|
| Milvus standalone | 19530 | gRPC vector operations |
| Attu UI | 8000 | Web-based Milvus management |
| etcd | 2379 | Metadata storage |
| MinIO | 9000 | Object storage for indexes |

---

## 8. Variant Annotation Pipeline

### ClinVar Integration

| Parameter | Value |
|---|---|
| Database | ClinVar (NCBI) |
| Total Variants | 4.1M clinical variants |
| Match Rate | ~35,616 / 3.5M variants (1.0%) |
| Classifications | Pathogenic, Likely pathogenic, VUS, Likely benign, Benign |
| Update Frequency | Monthly releases |
| Key Fields | Clinical significance, disease associations, review status |

### AlphaMissense Integration

| Parameter | Value |
|---|---|
| Database | AlphaMissense (DeepMind) |
| Total Predictions | 71,697,560 missense variant predictions |
| Match Rate | ~6,831 / 35,616 ClinVar variants (19.2%) |
| Model | AlphaFold-derived protein structure features |
| Output | Pathogenicity score (0.0-1.0) |

**AlphaMissense Thresholds:**

| Class | Score Range | Interpretation |
|---|---|---|
| Pathogenic | > 0.564 | Likely disease-causing |
| Ambiguous | 0.34 - 0.564 | Uncertain significance |
| Benign | < 0.34 | Likely neutral |

### Ensembl VEP Integration

| Parameter | Value |
|---|---|
| Tool | Ensembl Variant Effect Predictor (VEP) |
| Purpose | Functional consequence annotation |
| Output | Gene, transcript, consequence type, impact level |
| Impact Levels | HIGH, MODERATE, LOW, MODIFIER |
| Key Consequences | missense_variant, stop_gained, frameshift_variant, splice_donor_variant |

### Annotation Pipeline Code Pattern

```python
# From annotator.py — three-database annotation pipeline
def annotate_variants(vcf_path: str) -> List[AnnotatedVariant]:
    """
    Pipeline: VCF → ClinVar → AlphaMissense → VEP → Annotated variants
    """
    variants = parse_vcf(vcf_path, min_qual=30)        # ~3.5M pass filter
    variants = annotate_clinvar(variants)                # Clinical significance
    variants = annotate_alphamissense(variants)          # AI pathogenicity
    variants = annotate_vep(variants)                    # Functional consequences
    return variants
```

---

## 9. Knowledge Base — 201 Genes, 13 Therapeutic Areas

### Gene Distribution

| Therapeutic Area | Gene Count | Example Genes |
|---|---|---|
| Neurology | 36 | VCP, APP, PSEN1, MAPT, SOD1, FUS, C9orf72 |
| Oncology | 27 | EGFR, BRAF, KRAS, TP53, BRCA1, BRCA2, PIK3CA |
| Metabolic | 22 | GCK, PPARG, SLC2A2, ABCA1, PCSK9 |
| Infectious Disease | 21 | ACE2, CCR5, IFITM3, TLR4, TMPRSS2 |
| Respiratory | 13 | CFTR, SERPINA1, MUC5B, TERT |
| Rare Disease | 12 | VCP, HTT, SMN1, DMD, CFTR |
| Hematology | 12 | HBB, HBA1, F5, JAK2, CALR |
| GI/Hepatology | 12 | HFE, ATP7B, NOD2, SERPINA1 |
| Pharmacogenomics | 11 | CYP2D6, CYP2C19, CYP3A4, DPYD, TPMT |
| Ophthalmology | 11 | RHO, RPE65, RS1, ABCA4 |
| Cardiovascular | 10 | LDLR, PCSK9, SCN5A, MYBPC3, KCNQ1 |
| Immunology | 9 | HLA-B, TNF, IL6, JAK1, CTLA4 |
| Dermatology | 9 | FLG, MC1R, TYR, KRT14 |
| **Total** | **201** | **171 druggable (85% druggability)** |

### Knowledge Base Structure

Each gene entry in `knowledge.py` (88 KB) contains:

```python
{
    "gene": "VCP",
    "uniprot": "P55072",
    "therapeutic_area": "Neurology",
    "diseases": ["Frontotemporal Dementia", "ALS", "IBMPFD"],
    "druggability": "High",
    "drug_targets": ["D2 ATPase domain", "N-D1 interface"],
    "known_inhibitors": ["CB-5083", "NMS-873"],
    "variant_hotspots": ["R155H", "R191Q", "A232E"],
    "pathway": "Ubiquitin-proteasome system",
    "mechanism": "AAA+ ATPase, protein homeostasis",
}
```

### Query Expansion Maps

10 therapeutic area query expansion maps enrich user queries with domain-specific terminology:

```python
QUERY_EXPANSION = {
    "oncology": ["tumor", "cancer", "neoplasm", "carcinoma", "mutation driver",
                  "somatic", "germline", "tumor suppressor", "oncogene"],
    "neurology": ["neurodegeneration", "dementia", "ALS", "Parkinson",
                   "Alzheimer", "frontotemporal", "motor neuron"],
    # ... 8 more therapeutic areas
}
```

---

## 10. Anthropic Claude LLM Integration

### Configuration

| Parameter | Value |
|---|---|
| Model | `claude-sonnet-4-20250514` |
| Temperature | 0.3 |
| Max Tokens | 4096 |
| API | Anthropic Messages API |
| Role | RAG-grounded clinical reasoning |

### RAG Prompt Structure

```python
system_prompt = """You are a clinical genomics specialist analyzing patient
variant data. Ground all responses in the retrieved variant evidence and
knowledge base. Cite specific variants, genes, and clinical classifications.
When recommending drug targets, explain the evidence chain from variant
to disease mechanism to druggability assessment."""

user_prompt = f"""
## Retrieved Variant Evidence (top {top_k} matches)
{formatted_variants}

## Knowledge Base Context
{knowledge_context}

## User Question
{user_question}
"""
```

### Response Format

Claude generates structured target hypotheses:

```json
{
    "target_gene": "VCP",
    "confidence": "high",
    "evidence_chain": [
        "rs188935092 (chr9:35065263 G>A) — ClinVar: Pathogenic",
        "AlphaMissense: 0.87 (pathogenic, >0.564 threshold)",
        "VEP: missense_variant, HIGH impact",
        "Known drug target: CB-5083 (Phase I VCP inhibitor)",
        "Druggability: 0.92 (D2 ATPase domain, ~450Å³ pocket)"
    ],
    "therapeutic_area": "Neurology",
    "diseases": ["Frontotemporal Dementia", "ALS", "IBMPFD"],
    "recommended_action": "Proceed to drug discovery with VCP as primary target"
}
```

---

## 11. Stage 3: Drug Discovery Pipeline

### Overview

Stage 3 takes a target gene hypothesis from Stage 2 and produces 100 ranked novel drug candidates using BioNeMo generative chemistry, molecular docking, and drug-likeness scoring.

### 10-Stage Pipeline

| Stage | Process | Description |
|---|---|---|
| 1 | Initialize | Load target hypothesis, validate inputs |
| 2 | Normalize Target | Map gene → UniProt ID → PDB structures |
| 3 | Structure Discovery | Query RCSB PDB for Cryo-EM/X-ray structures |
| 4 | Structure Preparation | Score and rank structures, select best binding site |
| 5 | Molecule Generation | MolMIM generates novel SMILES from seed compound |
| 6 | Chemistry QC | RDKit validates chemical feasibility |
| 7 | Conformer Generation | RDKit 3D conformer embedding (ETKDG) |
| 8 | Molecular Docking | DiffDock predicts binding poses and affinities |
| 9 | Composite Ranking | Weighted scoring: 30% gen + 40% dock + 30% QED |
| 10 | Reporting | PDF report generation (ReportLab) |

### Pipeline Configuration

```python
# From pipeline.py
PIPELINE_CONFIG = {
    "num_candidates": 100,
    "molmim_endpoint": "http://localhost:8001/v1/generate",
    "diffdock_endpoint": "http://localhost:8002/v1/dock",
    "min_qed": 0.3,
    "min_dock_score": -6.0,         # kcal/mol
    "scoring_weights": {
        "generation": 0.30,
        "docking": 0.40,
        "qed": 0.30
    }
}
```

### UniProt Mappings

| Gene | UniProt ID | Function |
|---|---|---|
| VCP | P55072 | AAA+ ATPase, protein homeostasis |
| EGFR | P00533 | Receptor tyrosine kinase |
| BRAF | P15056 | Serine/threonine kinase |
| KRAS | P01116 | GTPase signaling |

---

## 12. BioNeMo NIM Services

### MolMIM (Port 8001) — Molecule Generation

| Parameter | Value |
|---|---|
| Endpoint | `POST http://localhost:8001/v1/generate` |
| Model | MolMIM (Molecular Masked Inverse Model) |
| Input | Seed SMILES string |
| Output | Novel SMILES candidates |
| Method | Masked language model on molecular tokens |
| Container | `nvcr.io/nvidia/clara/bionemo-molmim:1.0` |

**Request Format:**

```json
{
    "smiles": "CC(=O)Nc1ccc(O)cc1",
    "num_molecules": 100,
    "temperature": 0.7,
    "top_k": 50
}
```

**Response Format:**

```json
{
    "molecules": [
        {"smiles": "CC(=O)Nc1ccc(O)c(F)c1", "score": 0.85},
        {"smiles": "CC(=O)Nc1ccc(O)c(Cl)c1", "score": 0.82}
    ]
}
```

### DiffDock (Port 8002) — Molecular Docking

| Parameter | Value |
|---|---|
| Endpoint | `POST http://localhost:8002/v1/dock` |
| Model | DiffDock (diffusion-based docking) |
| Input | Ligand SMILES + protein PDB structure |
| Output | Binding pose + affinity score (kcal/mol) |
| Method | Score-based generative diffusion model |
| Container | `nvcr.io/nvidia/clara/diffdock:1.0` |

**Request Format:**

```json
{
    "ligand_smiles": "CC(=O)Nc1ccc(O)c(F)c1",
    "protein_pdb": "<PDB file content or path>",
    "num_poses": 5
}
```

**Response Format:**

```json
{
    "poses": [
        {"score": -8.7, "confidence": 0.92, "pose_pdb": "..."},
        {"score": -7.3, "confidence": 0.84, "pose_pdb": "..."}
    ]
}
```

### Docking Score Interpretation

| Score (kcal/mol) | Interpretation |
|---|---|
| -12 to -8 | Excellent binding affinity |
| -8 to -6 | Good binding affinity |
| -6 to -4 | Moderate binding affinity |
| > -4 | Weak binding affinity |

---

## 13. Drug-Likeness Scoring

### Lipinski's Rule of Five

| Rule | Threshold | Description |
|---|---|---|
| Molecular Weight | ≤ 500 Da | Oral absorption limit |
| LogP | ≤ 5 | Lipophilicity |
| H-Bond Donors | ≤ 5 | NH + OH groups |
| H-Bond Acceptors | ≤ 10 | N + O atoms |

### QED (Quantitative Estimate of Drug-likeness)

| Range | Interpretation |
|---|---|
| > 0.67 | Drug-like (favorable properties) |
| 0.49 - 0.67 | Moderate drug-likeness |
| < 0.49 | Less drug-like |

### TPSA (Topological Polar Surface Area)

| Range (Å²) | Interpretation |
|---|---|
| < 140 | Good oral bioavailability |
| 60-90 | Optimal range |
| > 140 | Poor oral absorption |

### Composite Scoring Formula

```python
def compute_composite_score(gen_score, dock_score, qed_score):
    """
    Weighted composite: 30% generation + 40% docking + 30% QED

    Docking normalization: scale raw kcal/mol to 0-1 range
    dock_normalized = max(0, min(1, (10 + dock_score) / 20))

    Example: dock_score = -8.5 → (10 + (-8.5)) / 20 = 0.075 → normalized
    """
    dock_normalized = max(0.0, min(1.0, (10.0 + dock_score) / 20.0))

    composite = (
        0.30 * gen_score +
        0.40 * dock_normalized +
        0.30 * qed_score
    )
    return composite
```

### RDKit Property Calculation

```python
from rdkit import Chem
from rdkit.Chem import Descriptors, QED

def calculate_properties(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "tpsa": Descriptors.TPSA(mol),
        "qed": QED.qed(mol),
        "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
        "lipinski_pass": all([
            Descriptors.MolWt(mol) <= 500,
            Descriptors.MolLogP(mol) <= 5,
            Descriptors.NumHDonors(mol) <= 5,
            Descriptors.NumHAcceptors(mol) <= 10,
        ])
    }
```

---

## 14. Cryo-EM Structure Evidence

### Structure Scoring Algorithm

The pipeline automatically retrieves and scores protein structures from RCSB PDB:

```python
def score_structure(structure: StructureInfo) -> float:
    """
    Score a PDB structure for suitability in drug discovery.

    Factors:
    - Resolution: lower is better (max 5 Å cutoff)
    - Inhibitor-bound: +3 bonus (binding site already defined)
    - Druggable pockets: +0.5 per pocket
    - Cryo-EM method: +0.5 (modern technique bonus)
    """
    score += max(0, 5.0 - resolution)               # Resolution: 0-5 scale
    if has_inhibitor_bound:
        score += 3.0
    score += num_druggable_pockets * 0.5
    if 'Cryo-EM' in method:
        score += 0.5
    return score
```

### VCP Structures (Demo)

| PDB ID | Resolution | Method | Description | Score |
|---|---|---|---|---|
| 8OOI | 2.9 Å | Cryo-EM | WT VCP hexamer | High |
| 9DIL | 3.2 Å | Cryo-EM | Mutant VCP | High |
| 7K56 | 2.5 Å | Cryo-EM | VCP complex | Highest |
| 5FTK | 2.3 Å | X-ray | VCP + CB-5083 inhibitor | Highest (inhibitor-bound) |

### VCP Binding Site

| Parameter | Value |
|---|---|
| Domain | D2 ATPase domain |
| Mechanism | ATP-competitive inhibition |
| Pocket Volume | ~450 Å³ |
| Druggability Score | 0.92 |
| Key Residues | ALA464, GLY479, ASP320, GLY215 |

---

## 15. VCP/FTD Demo Walkthrough

### Demo Target: Valosin-Containing Protein (VCP/p97)

| Parameter | Value |
|---|---|
| Gene | VCP |
| Protein | p97 / Valosin-Containing Protein |
| UniProt | P55072 |
| Function | AAA+ ATPase, ubiquitin-proteasome pathway |
| Diseases | Frontotemporal Dementia (FTD), ALS, IBMPFD |
| Variant | rs188935092 (chr9:35065263 G>A) |
| ClinVar | Pathogenic |
| AlphaMissense | 0.87 (pathogenic, >0.564 threshold) |
| Seed Compound | CB-5083 (Phase I clinical VCP inhibitor) |

### Demo Flow

**Stage 1 — Genomics (Demo Mode: ~20 min):**
1. Load pre-processed HG002 FASTQ subset
2. Run Parabricks fq2bam alignment
3. Run DeepVariant variant calling
4. Output VCF with ~11.7M variants including rs188935092

**Stage 2 — RAG/Chat (Interactive):**
1. VCF annotated: ClinVar flags rs188935092 as pathogenic in VCP
2. AlphaMissense scores the missense variant at 0.87 (pathogenic)
3. 3.5M variants embedded and indexed in Milvus
4. User queries: "What are the most promising drug targets in this patient's genome?"
5. Claude identifies VCP with full evidence chain
6. Target hypothesis: VCP → FTD → druggable D2 ATPase domain

**Stage 3 — Drug Discovery (~10 min):**
1. VCP → UniProt P55072 → PDB structure retrieval
2. Cryo-EM structures scored: 8OOI, 9DIL, 7K56, 5FTK
3. 5FTK selected (inhibitor-bound, highest score)
4. CB-5083 seed SMILES → MolMIM generates 100 novel analogs
5. RDKit validates Lipinski, QED, TPSA
6. DiffDock docks each candidate against VCP D2 domain
7. Composite ranking: 30% generation + 40% docking + 30% QED
8. Top candidates: novel VCP inhibitors with improved drug-likeness
9. PDF report generated via ReportLab

### Expected Demo Output

```
Pipeline: HCLS AI Factory — VCP/FTD Demo
Target: VCP (P55072) — Frontotemporal Dementia
Seed: CB-5083 (ATP-competitive VCP inhibitor)
Structure: 5FTK (2.3 Å, X-ray, inhibitor-bound)

Results:
- 100 novel VCP inhibitor candidates generated
- 87 pass Lipinski's Rule of Five
- 72 have QED > 0.67 (drug-like)
- Top 10 show docking scores -8.2 to -11.4 kcal/mol
- Composite scores range 0.68-0.89
```

---

## 16. Pydantic Data Models

### Core Models (from `models.py`)

```python
from pydantic import BaseModel, Field
from typing import List, Optional
from enum import Enum

class TargetHypothesis(BaseModel):
    """Output from Stage 2 — RAG-identified drug target"""
    gene: str                          # e.g., "VCP"
    uniprot_id: str                    # e.g., "P55072"
    confidence: str                    # high, medium, low
    evidence_chain: List[str]          # Supporting evidence items
    therapeutic_area: str              # e.g., "Neurology"
    diseases: List[str]               # Associated conditions
    druggability_score: float          # 0-1 scale

class StructureInfo(BaseModel):
    """PDB structure metadata"""
    pdb_id: str                        # e.g., "8OOI"
    resolution: float                  # Angstroms
    method: str                        # ELECTRON MICROSCOPY, X-RAY DIFFRACTION
    title: str
    has_inhibitor: bool
    num_pockets: int
    score: float                       # Computed suitability score

class StructureManifest(BaseModel):
    """Collection of scored PDB structures for a target"""
    target_gene: str
    uniprot_id: str
    structures: List[StructureInfo]
    best_structure: str                # PDB ID of highest-scored

class MoleculeProperties(BaseModel):
    """RDKit-computed molecular properties"""
    smiles: str
    molecular_weight: float
    logp: float
    hbd: int                           # H-bond donors
    hba: int                           # H-bond acceptors
    tpsa: float                        # Topological polar surface area
    qed: float                         # Quantitative drug-likeness
    rotatable_bonds: int
    lipinski_pass: bool

class GeneratedMolecule(BaseModel):
    """MolMIM output — a novel molecule candidate"""
    smiles: str
    generation_score: float            # MolMIM confidence
    properties: Optional[MoleculeProperties]

class DockingResult(BaseModel):
    """DiffDock output — binding prediction"""
    ligand_smiles: str
    dock_score: float                  # kcal/mol (negative = better)
    confidence: float                  # 0-1 model confidence
    pose_pdb: Optional[str]            # PDB-format binding pose

class RankedCandidate(BaseModel):
    """Final ranked drug candidate with composite score"""
    rank: int
    smiles: str
    generation_score: float
    dock_score: float
    qed: float
    composite_score: float             # 30% gen + 40% dock + 30% QED
    lipinski_pass: bool
    molecular_weight: float
    logp: float

class PipelineConfig(BaseModel):
    """Pipeline execution configuration"""
    mode: str                          # full, target, drug, demo, genomics_only
    num_candidates: int = 100
    min_qed: float = 0.3
    min_dock_score: float = -6.0
    molmim_url: str = "http://localhost:8001/v1/generate"
    diffdock_url: str = "http://localhost:8002/v1/dock"

class PipelineRun(BaseModel):
    """Complete pipeline execution record"""
    run_id: str
    mode: str
    target: Optional[TargetHypothesis]
    structures: Optional[StructureManifest]
    candidates: List[RankedCandidate]
    total_generated: int
    total_passed_qc: int
    total_docked: int
    duration_seconds: float
    status: str                        # running, completed, failed
```

---

## 17. Nextflow DSL2 Orchestration

### Pipeline Modes

| Mode | Stages | Description |
|---|---|---|
| `full` | 1 → 2 → 3 | Complete end-to-end pipeline |
| `target` | 2 → 3 | Skip genomics, use existing VCF |
| `drug` | 3 only | Skip to drug discovery with known target |
| `demo` | 1 → 2 → 3 | Pre-configured VCP/FTD demonstration |
| `genomics_only` | 1 only | Run only variant calling |

### Main Pipeline Entry (`main.nf`)

```groovy
#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { GENOMICS_PIPELINE } from './modules/genomics'
include { RAG_CHAT_PIPELINE } from './modules/rag_chat'
include { DRUG_DISCOVERY_PIPELINE } from './modules/drug_discovery'
include { REPORTING } from './modules/reporting'

workflow {
    if (params.mode in ['full', 'demo', 'genomics_only']) {
        GENOMICS_PIPELINE(
            params.fastq_r1,
            params.fastq_r2,
            params.reference
        )
    }

    if (params.mode in ['full', 'demo', 'target']) {
        RAG_CHAT_PIPELINE(
            params.mode == 'target' ? params.vcf : GENOMICS_PIPELINE.out.vcf
        )
    }

    if (params.mode in ['full', 'demo', 'target', 'drug']) {
        DRUG_DISCOVERY_PIPELINE(
            params.mode == 'drug' ? params.target_gene : RAG_CHAT_PIPELINE.out.target
        )
    }

    REPORTING(
        DRUG_DISCOVERY_PIPELINE.out.candidates
    )
}
```

### Nextflow Profiles

| Profile | Description |
|---|---|
| `standard` | Default local execution |
| `docker` | Docker container execution |
| `singularity` | Singularity container execution |
| `dgx_spark` | DGX Spark optimized (GPU resources) |
| `slurm` | HPC cluster submission |
| `test` | Minimal test data |

### Pipeline Launcher (`run_pipeline.py`)

```bash
# Full pipeline
python run_pipeline.py --mode full \
  --fastq-r1 /data/HG002_R1.fastq.gz \
  --fastq-r2 /data/HG002_R2.fastq.gz \
  --reference /reference/GRCh38.fa

# Demo mode (pre-configured VCP/FTD)
python run_pipeline.py --mode demo

# Drug discovery only (known target)
python run_pipeline.py --mode drug --target-gene VCP --seed-smiles "CC(=O)..."
```

---

## 18. Landing Page & Service Health

### Landing Page (Port 8080)

The Flask-based landing page serves as the entry point for the HCLS AI Factory:

- **URL:** `http://localhost:8080`
- **Framework:** Flask
- **Features:**
  - 10-service health status dashboard
  - Pipeline mode selector (full, target, drug, demo)
  - Quick-start links to all service UIs
  - Real-time service status with green/red indicators
  - Pipeline execution history

### Service Health Check Implementation

```python
SERVICES = [
    {"name": "Parabricks Portal", "url": "http://localhost:5000/health", "port": 5000},
    {"name": "Milvus Vector DB", "url": "http://localhost:19530", "port": 19530},
    {"name": "RAG API", "url": "http://localhost:5001/health", "port": 5001},
    {"name": "Streamlit Chat", "url": "http://localhost:8501/healthz", "port": 8501},
    {"name": "MolMIM NIM", "url": "http://localhost:8001/v1/health/ready", "port": 8001},
    {"name": "DiffDock NIM", "url": "http://localhost:8002/v1/health/ready", "port": 8002},
    {"name": "Discovery UI", "url": "http://localhost:8505/healthz", "port": 8505},
    {"name": "Grafana", "url": "http://localhost:3000/api/health", "port": 3000},
    {"name": "Prometheus", "url": "http://localhost:9099/-/healthy", "port": 9099},
    {"name": "DCGM Exporter", "url": "http://localhost:9400/metrics", "port": 9400},
]
```

---

## 19. Monitoring Stack

### Grafana (Port 3000)

| Parameter | Value |
|---|---|
| Image | `grafana/grafana:10.2.2` |
| Default User | admin / changeme |
| Dashboards | HCLS AI Factory (GPU, pipeline, services) |
| Data Source | Prometheus |

### Prometheus (Port 9099)

| Parameter | Value |
|---|---|
| Image | `prom/prometheus:v2.48.0` |
| Internal Port | 9090 → External 9099 |
| Retention | 30 days |
| Scrape Targets | Node Exporter, DCGM Exporter, service metrics |

### Node Exporter (Port 9100)

| Metric Category | Examples |
|---|---|
| CPU | Usage %, load average, core temperatures |
| Memory | Used/free/cached, swap usage |
| Disk | I/O throughput, NVMe utilization, space |
| Network | Bandwidth, packet rates, error rates |

### DCGM Exporter (Port 9400)

| Metric | Description |
|---|---|
| `DCGM_FI_DEV_GPU_UTIL` | GPU utilization percentage |
| `DCGM_FI_DEV_FB_USED` | GPU memory used (bytes) |
| `DCGM_FI_DEV_FB_FREE` | GPU memory free (bytes) |
| `DCGM_FI_DEV_GPU_TEMP` | GPU temperature (°C) |
| `DCGM_FI_DEV_POWER_USAGE` | GPU power draw (watts) |
| `DCGM_FI_DEV_SM_CLOCK` | SM clock frequency (MHz) |

### Key Dashboard Panels

1. **GPU Utilization Timeline** — fq2bam (70-90%) → DeepVariant (80-95%) → idle → MolMIM/DiffDock bursts
2. **Pipeline Stage Progress** — Stage 1/2/3 completion with timing
3. **Memory Pressure** — Unified memory usage across CPU + GPU workloads
4. **Service Health Grid** — Green/red status for all 10 services
5. **Variant Processing Rate** — Variants annotated per second
6. **Drug Discovery Throughput** — Molecules generated/docked per minute

---

## 20. Cross-Modal Integration

### HCLS AI Factory Ecosystem

The genomics-to-drug-discovery pipeline integrates with the broader HCLS AI Factory:

```
Imaging Intelligence Agent (CT/MRI/X-Ray)
    │
    ├── Lung-RADS 4B+ finding
    │       ↓
    │   FHIR ServiceRequest
    │       ↓
    ├── Trigger genomics analysis (Parabricks)
    │       ↓
    │   Tumor gene profiling
    │       ↓
    └── Drug candidates → Combined imaging + genomics report
```

### Cross-Modal Triggers

| Trigger | Source | Target | Action |
|---|---|---|---|
| Lung-RADS 4B+ | Imaging Agent | Genomics Pipeline | Initiate tumor profiling |
| Pathogenic Variant | Genomics Pipeline | Drug Discovery | Generate targeted therapies |
| Drug Candidates | Drug Discovery | Imaging Agent | Combined clinical report |

### NVIDIA FLARE — Federated Learning

For multi-site deployments (Phase 3), NVIDIA FLARE enables federated model training:
- Models train locally at each site
- Only model updates (not patient data) are shared
- Aggregation server combines updates
- Privacy-preserving: raw genomic data never leaves the institution

---

## 21. Configuration Reference

### Environment Variables

| Variable | Default | Description |
|---|---|---|
| `ANTHROPIC_API_KEY` | (required) | Anthropic API key for Claude |
| `NGC_API_KEY` | (required) | NVIDIA NGC key for BioNeMo NIMs |
| `REFERENCE_GENOME` | `/reference/GRCh38.fa` | Path to GRCh38 reference |
| `MILVUS_HOST` | `localhost` | Milvus server hostname |
| `MILVUS_PORT` | `19530` | Milvus gRPC port |
| `MOLMIM_URL` | `http://localhost:8001` | MolMIM NIM endpoint |
| `DIFFDOCK_URL` | `http://localhost:8002` | DiffDock NIM endpoint |
| `CLAUDE_MODEL` | `claude-sonnet-4-20250514` | Claude model identifier |
| `CLAUDE_TEMPERATURE` | `0.3` | LLM temperature |
| `PIPELINE_MODE` | `full` | Pipeline execution mode |
| `NUM_CANDIDATES` | `100` | Number of drug candidates to generate |
| `MIN_QED` | `0.3` | Minimum QED threshold |
| `MIN_DOCK_SCORE` | `-6.0` | Minimum docking score (kcal/mol) |
| `GRAFANA_USER` | `admin` | Grafana admin username |
| `GRAFANA_PASSWORD` | `changeme` | Grafana admin password |

### AlphaMissense Thresholds

```python
AM_PATHOGENIC_THRESHOLD = 0.564
AM_AMBIGUOUS_LOWER = 0.34
AM_AMBIGUOUS_UPPER = 0.564
AM_BENIGN_THRESHOLD = 0.34
```

### Scoring Weights

```python
SCORING_WEIGHTS = {
    "generation": 0.30,   # MolMIM generation confidence
    "docking": 0.40,      # DiffDock binding affinity
    "qed": 0.30           # RDKit drug-likeness
}
```

### Drug-Likeness Thresholds

```python
LIPINSKI = {
    "max_mw": 500,        # Molecular weight (Da)
    "max_logp": 5,        # Partition coefficient
    "max_hbd": 5,         # H-bond donors
    "max_hba": 10         # H-bond acceptors
}

QED_THRESHOLDS = {
    "drug_like": 0.67,    # QED > 0.67
    "moderate": 0.49,     # 0.49 < QED < 0.67
    "less_drug_like": 0   # QED < 0.49
}

DOCKING_THRESHOLDS = {
    "excellent": -8.0,    # kcal/mol
    "good": -6.0,
    "moderate": -4.0,
    "minimum": -6.0       # Pipeline cutoff
}
```

---

## 22. Deployment Roadmap

### Phase 1: Proof Build

| Parameter | Value |
|---|---|
| Hardware | NVIDIA DGX Spark ($3,999) |
| Orchestration | Docker Compose |
| Scale | Single patient, sequential processing |
| Timeline | Proof of concept |
| GPU | 1× GB10 |
| Memory | 128 GB unified |

### Phase 2: Departmental

| Parameter | Value |
|---|---|
| Hardware | 1-2× DGX B200 |
| Orchestration | Kubernetes |
| Scale | Multiple concurrent patients |
| GPU | 8× B200 per node |
| Memory | 1-2 TB HBM3e |
| Networking | InfiniBand |

### Phase 3: Enterprise / Multi-Site

| Parameter | Value |
|---|---|
| Hardware | DGX SuperPOD |
| Orchestration | Kubernetes + NVIDIA FLARE |
| Scale | Thousands of concurrent patients |
| GPU | Hundreds of B200 GPUs |
| Networking | InfiniBand fabric |
| Privacy | Federated learning (data stays local) |

### Scaling Considerations

| Bottleneck | Phase 1 Solution | Phase 2+ Solution |
|---|---|---|
| Genomics throughput | Sequential (1 sample) | Parallel Parabricks instances |
| Milvus query latency | Single-node Milvus | Milvus cluster with sharding |
| BioNeMo inference | Single NIM per model | Multiple NIM replicas |
| Storage I/O | NVMe direct | GPUDirect Storage + RAID |

---

## 23. Testing Strategy

### Unit Tests

| Component | Test Focus |
|---|---|
| VCF Parser | Variant extraction, quality filtering |
| Annotator | ClinVar/AlphaMissense/VEP lookup accuracy |
| Milvus Client | Index creation, search recall |
| MolMIM Client | SMILES generation, request format |
| DiffDock Client | Docking request/response parsing |
| RDKit Scoring | Lipinski, QED, TPSA calculations |
| Composite Scorer | Weight application, normalization |

### Integration Tests

| Test | Validates |
|---|---|
| VCF → Annotation → Milvus | End-to-end Stage 2 pipeline |
| Target → PDB → MolMIM → DiffDock | End-to-end Stage 3 pipeline |
| Health check endpoints | All 10 services responding |
| Nextflow modes | full, target, drug, demo execution |

### Demo Mode Validation

The `demo` pipeline mode uses pre-configured inputs to validate the complete pipeline:
- Input: HG002 FASTQ subset (smaller dataset for faster execution)
- Expected: VCP identified as target with rs188935092 evidence
- Expected: 100 novel VCP inhibitor candidates ranked
- Validation: Top candidates show improved QED vs CB-5083 seed

---

## 24. Implementation Sequence

### Recommended Build Order

1. **Infrastructure:** Docker Compose, Milvus, monitoring stack
2. **Stage 1 — Genomics:** Parabricks container, fq2bam, DeepVariant, VCF output
3. **Stage 2 — Annotation:** ClinVar + AlphaMissense + VEP pipeline
4. **Stage 2 — Vector DB:** Milvus schema, BGE embedding, IVF_FLAT index
5. **Stage 2 — RAG:** Claude integration, knowledge base, query expansion
6. **Stage 2 — Chat UI:** Streamlit interface, REST API
7. **Stage 3 — Structure:** RCSB PDB retrieval, Cryo-EM scoring
8. **Stage 3 — Generation:** MolMIM NIM, molecule generation
9. **Stage 3 — Docking:** DiffDock NIM, binding prediction
10. **Stage 3 — Scoring:** RDKit properties, composite ranking
11. **Stage 3 — Reporting:** PDF generation, Discovery UI
12. **Orchestration:** Nextflow DSL2, pipeline modes, landing page
13. **Testing:** Unit tests, integration tests, demo mode validation
14. **Monitoring:** Grafana dashboards, alerting rules

### Key Dependencies

```
GRCh38 reference → BWA-MEM2 index → fq2bam alignment
ClinVar + AlphaMissense databases → Annotation pipeline
Milvus running → Embedding indexing → RAG queries
BioNeMo NIMs running → Molecule generation + docking
All services healthy → Landing page green status
```

---

*This Project Bible is the authoritative technical reference for the HCLS AI Factory. All other documentation assets (White Paper, Demo Guide, Intelligence Report, Learning Guides) derive their technical details from this source.*
