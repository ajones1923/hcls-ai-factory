---
tags:
  - Genomics
  - Parabricks
  - DeepVariant
  - Docker
  - VCF
---

# Genomics Pipeline

[![NVIDIA Parabricks](https://img.shields.io/badge/NVIDIA-Parabricks%204.6-76B900?style=flat&logo=nvidia)](https://docs.nvidia.com/clara/parabricks/)
[![DGX Spark](https://img.shields.io/badge/NVIDIA-DGX%20Spark-76B900?style=flat&logo=nvidia)](https://www.nvidia.com/en-us/data-center/dgx-spark/)
[![Docker](https://img.shields.io/badge/Docker-Required-2496ED?style=flat&logo=docker)](https://www.docker.com/)
[![Python](https://img.shields.io/badge/Python-3.11+-3776AB?style=flat&logo=python)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](../../LICENSE)

**Stage 1 of the Precision Medicine to Drug Discovery AI Factory**

> GPU-accelerated germline variant calling pipeline using NVIDIA Parabricks. This pipeline transforms raw DNA sequencing data (FASTQ) into analysis-ready variant calls (VCF) in under 2 hours using GPU acceleration.

```
┌──────────────────────────────────────────────────────────────────────────────────────┐
│                    PRECISION MEDICINE TO DRUG DISCOVERY AI FACTORY                   │
├──────────────────────────────────────────────────────────────────────────────────────┤
│                                                                                      │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────────────┐   │
│  │  GENOMICS   │    │  RAG/CHAT   │    │   CRYO-EM   │    │ MOLECULE GENERATION │   │
│  │  PIPELINE   │───▶│  PIPELINE   │───▶│  EVIDENCE   │───▶│     (BioNeMo)       │   │
│  │ (This Repo) │    │             │    │             │    │                     │   │
│  └─────────────┘    └─────────────┘    └─────────────┘    └─────────────────────┘   │
│    FASTQ→VCF         VCF→Target        Target→Structure    Structure→Molecules      │
│    Parabricks        Milvus+Claude     PDB/EMDB            MolMIM+DiffDock          │
│                                                                                      │
└──────────────────────────────────────────────────────────────────────────────────────┘
```

---

## Table of Contents

- [The Biological Foundation](#the-biological-foundation)
- [From Biology to Digital Data](#from-biology-to-digital-data)
- [What This Pipeline Does](#what-this-pipeline-does)
- [Key Features](#key-features)
- [Architecture](#architecture)
- [System Requirements](#system-requirements)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps in Detail](#pipeline-steps-in-detail)
- [Understanding the Output VCF](#understanding-the-output-vcf)
- [Configuration](#configuration)
- [Directory Structure](#directory-structure)
- [Performance Benchmarks](#performance-benchmarks)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Related Pipelines](#related-pipelines)
- [References](#references)

---

## The Biological Foundation

### From DNA to Disease: Why This Matters

Every cell in the human body contains DNA—the blueprint of human biology. All of that DNA together is called the **genome**. DNA is organized into **chromosomes**, and chromosomes contain **genes**. Genes encode **proteins**, and proteins are the molecular machines that make the body work.

When **variants** occur in genes, they can change how proteins function, and those changes can disrupt normal biology and lead to disease. Understanding these changes—and doing it at scale—is the foundation of modern precision medicine.

```
DNA (Blueprint)
    │
    ▼
Chromosomes (23 pairs in humans)
    │
    ▼
Genes (~20,000 protein-coding genes)
    │
    ▼
Proteins (Molecular machines)
    │
    ▼
Function (Normal biology or disease)
```

### The Scale of Human Variation

A typical human genome contains approximately:
- **3 billion** base pairs (A, C, G, T)
- **4-5 million** variants compared to the reference genome
- **~20,000** protein-coding genes
- **~11.7 million** total variants (including non-coding regions) in a 30x whole-genome sequence

Most variants are harmless—normal human diversity. But some affect genes in ways that change proteins, and those changes are often where disease begins. This pipeline is the first step in finding those critical variants.

---

## From Biology to Digital Data

### The Sequencing Process

We start with data from **Genome in a Bottle (GIAB)**, a globally trusted reference initiative led by NIST that provides high-confidence human genome datasets used worldwide to validate sequencing and analysis. The genome we're using here, **HG002**, was generated with Illumina short-read sequencing—the clinical standard today.

**What happens during sequencing:**

1. **Sample Collection**: DNA is extracted from cells (blood, saliva, tissue)
2. **Library Preparation**: DNA is broken into millions of small fragments (~250bp each)
3. **Sequencing**: Each fragment is read as a string of DNA letters (A, C, G, T)
4. **Quality Scoring**: Each base call gets a confidence score (Phred quality)
5. **Output**: FASTQ files containing hundreds of millions of reads

**The challenge:** A single 30x whole-genome sequence produces:
- **~800 million** read pairs
- **~200 GB** of raw data
- Takes **24-48 hours** on traditional CPU pipelines

This is where GPU acceleration transforms what's possible.

### What is Secondary Analysis?

**Primary analysis** happens on the sequencer—converting light signals to base calls.

**Secondary analysis** (this pipeline) turns raw reads into meaningful variants:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         SECONDARY ANALYSIS PIPELINE                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  FASTQ Files              BAM File                 VCF File                 │
│  (Raw Reads)              (Aligned Reads)          (Variants)               │
│                                                                             │
│  ┌──────────┐    Align    ┌──────────┐    Call    ┌──────────┐            │
│  │ ACGTACGT │   ──────▶   │ chr7:pos │   ──────▶  │ G→A at   │            │
│  │ TGCATGCA │   BWA-MEM2  │ chr7:pos │  DeepVar   │ chr7:pos │            │
│  │ GCTAGCTA │             │ chr7:pos │            │ Quality  │            │
│  │   ...    │             │   ...    │            │ Score    │            │
│  └──────────┘             └──────────┘            └──────────┘            │
│                                                                             │
│  "What letters            "Where do they          "How does this           │
│   were read?"              belong?"                differ from             │
│                                                    reference?"             │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## What This Pipeline Does

This pipeline processes whole-genome sequencing (WGS) data through a complete bioinformatics workflow using **NVIDIA Parabricks**, which accelerates genome analysis on GPUs—turning what used to take days on CPUs into hours.

### The Four Steps

| Step | Tool | Input | Output | Time (GPU) |
|------|------|-------|--------|------------|
| **1. Alignment** | BWA-MEM2 | FASTQ reads | Mapped positions | 20-45 min |
| **2. Sorting & Dedup** | fq2bam | Mapped reads | Sorted BAM | (included) |
| **3. Indexing** | samtools | BAM file | BAM index + QC | 2-5 min |
| **4. Variant Calling** | DeepVariant | BAM file | VCF file | 10-35 min |

**Total: 120-240 minutes** (vs. 24-48 hours on CPU)

### Pipeline Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           NVIDIA PARABRICKS PIPELINE                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   Input: FASTQ Files (HG002_R1.fastq.gz, HG002_R2.fastq.gz)                │
│          ~200GB paired-end reads from Illumina sequencing                   │
│                                                                             │
│                              ┌─────────────────┐                            │
│                              │   Reference     │                            │
│                              │   GRCh38.fa     │                            │
│                              │   (3.1 GB)      │                            │
│                              └────────┬────────┘                            │
│                                       │                                     │
│   ┌────────────┐                      ▼                                     │
│   │   FASTQ    │    ┌─────────────────────────────────────────┐            │
│   │   Files    │───▶│              fq2bam                      │            │
│   │            │    │  • BWA-MEM2 alignment (GPU accelerated)  │            │
│   │  R1 + R2   │    │  • Coordinate sorting                    │            │
│   │  ~200 GB   │    │  • PCR duplicate marking                 │            │
│   └────────────┘    └──────────────────┬──────────────────────┘            │
│                                        │                                    │
│                                        ▼                                    │
│                     ┌─────────────────────────────────────────┐            │
│                     │            BAM File                      │            │
│                     │  • Aligned reads with coordinates        │            │
│                     │  • Quality scores preserved              │            │
│                     │  • ~100 GB output                        │            │
│                     └──────────────────┬──────────────────────┘            │
│                                        │                                    │
│                                        ▼                                    │
│                     ┌─────────────────────────────────────────┐            │
│                     │          samtools index                  │            │
│                     │  • Create BAM index (.bai)               │            │
│                     │  • Generate alignment statistics         │            │
│                     └──────────────────┬──────────────────────┘            │
│                                        │                                    │
│                                        ▼                                    │
│                     ┌─────────────────────────────────────────┐            │
│                     │           DeepVariant                    │            │
│                     │  • Deep learning variant caller          │            │
│                     │  • GPU-accelerated inference             │            │
│                     │  • Trained on millions of examples       │            │
│                     └──────────────────┬──────────────────────┘            │
│                                        │                                    │
│                                        ▼                                    │
│   Output: VCF File (HG002.genome.vcf.gz)                                   │
│           ~11.7 million variants identified                                 │
│           Ready for annotation in RAG/Chat Pipeline                         │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Key Features

### GPU Acceleration
- **10-50x faster** than CPU-only pipelines using NVIDIA Parabricks
- Full genome analysis in **120-240 minutes** vs. 24-48 hours
- Optimized for NVIDIA DGX Spark, A100, V100, and consumer GPUs

### Production-Ready Tools
- **BWA-MEM2**: Industry-standard aligner used in clinical labs worldwide
- **DeepVariant**: Google's deep learning variant caller (>99% accuracy)
- **GRCh38**: Latest human reference genome build

### Validated Output
- Uses **GIAB HG002** benchmark data with known truth sets
- Enables accuracy benchmarking against gold-standard calls
- Suitable for clinical genomics, research, and pharma applications

### Dual Interface
- **Command Line**: Scriptable pipeline for automation
- **Web Portal**: Visual interface with real-time monitoring

### Containerized & Reproducible
- Fully Dockerized with NVIDIA Container Runtime
- Version-controlled container images (Parabricks 4.6.0)
- Consistent results across different systems

---

## Architecture

### System Architecture

```
┌────────────────────────────────────────────────────────────────────────────┐
│                              HOST SYSTEM                                    │
├────────────────────────────────────────────────────────────────────────────┤
│                                                                            │
│   ┌────────────────────────────────────────────────────────────────────┐  │
│   │                      WEB PORTAL (Optional)                          │  │
│   │   Flask Server (Port 5000)                                          │  │
│   │   • Real-time pipeline monitoring     • Log streaming               │  │
│   │   • GPU utilization display           • Configuration management    │  │
│   │   • Progress tracking                 • One-click step execution    │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                    │                                       │
│   ┌────────────────────────────────────────────────────────────────────┐  │
│   │                      PIPELINE SCRIPTS                               │  │
│   │                                                                     │  │
│   │   00-setup-check.sh ──▶ 01-ngc-login.sh ──▶ 02-download-data.sh   │  │
│   │          │                                                          │  │
│   │          ▼                                                          │  │
│   │   03-setup-reference.sh ──▶ 04-run-chr20-test.sh                   │  │
│   │          │                                                          │  │
│   │          ▼                                                          │  │
│   │   05-run-full-genome.sh ──▶ Output: VCF File                       │  │
│   │                                                                     │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                    │                                       │
│   ┌────────────────────────────────────────────────────────────────────┐  │
│   │              DOCKER + NVIDIA CONTAINER RUNTIME                      │  │
│   │   ┌────────────────────────────────────────────────────────────┐   │  │
│   │   │         clara-parabricks:4.6.0-1 Container                 │   │  │
│   │   │                                                            │   │  │
│   │   │   ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐ │   │  │
│   │   │   │ BWA-MEM2 │  │ samtools │  │DeepVariant│ │ bcftools │ │   │  │
│   │   │   │  (GPU)   │  │          │  │  (GPU)   │  │          │ │   │  │
│   │   │   └──────────┘  └──────────┘  └──────────┘  └──────────┘ │   │  │
│   │   │                                                            │   │  │
│   │   └────────────────────────────────────────────────────────────┘   │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                    │                                       │
│   ┌────────────────────────────────────────────────────────────────────┐  │
│   │                         NVIDIA GPU                                  │  │
│   │                                                                     │  │
│   │   CUDA Cores ─────── Tensor Cores ─────── GPU Memory               │  │
│   │   (Alignment)        (DeepVariant)        (Data + Models)          │  │
│   │                                                                     │  │
│   │   Supported: DGX Spark (GB10), A100, V100, RTX 4090, RTX 3090     │  │
│   │                                                                     │  │
│   └────────────────────────────────────────────────────────────────────┘  │
│                                                                            │
└────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              DATA FLOW                                       │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│   GIAB FTP Server                    NVIDIA Sample Bundle                   │
│        │                                     │                              │
│        ▼                                     ▼                              │
│   ┌──────────┐                        ┌──────────────┐                     │
│   │  FASTQ   │                        │  Reference   │                     │
│   │  Files   │                        │  GRCh38.fa   │                     │
│   │  ~200GB  │                        │  + Indexes   │                     │
│   └────┬─────┘                        └──────┬───────┘                     │
│        │                                     │                              │
│        └──────────────┬──────────────────────┘                              │
│                       │                                                     │
│                       ▼                                                     │
│              ┌─────────────────┐                                           │
│              │    fq2bam       │                                           │
│              │  (Alignment)    │                                           │
│              └────────┬────────┘                                           │
│                       │                                                     │
│                       ▼                                                     │
│              ┌─────────────────┐                                           │
│              │   BAM File      │                                           │
│              │   (~100 GB)     │                                           │
│              └────────┬────────┘                                           │
│                       │                                                     │
│                       ▼                                                     │
│              ┌─────────────────┐                                           │
│              │  DeepVariant    │                                           │
│              │ (Variant Call)  │                                           │
│              └────────┬────────┘                                           │
│                       │                                                     │
│                       ▼                                                     │
│              ┌─────────────────┐                                           │
│              │   VCF File      │──────────▶  RAG/Chat Pipeline            │
│              │ ~11.7M variants │             (Stage 2)                     │
│              └─────────────────┘                                           │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## System Requirements

### Hardware Requirements

| Component | Minimum | Recommended (DGX Spark) | Notes |
|-----------|---------|-------------------------|-------|
| **GPU** | 8GB VRAM | 128GB (GB10) | More VRAM = faster processing |
| **GPU Architecture** | Volta (V100)+ | Blackwell (GB10) | Newer = better DeepVariant perf |
| **System RAM** | 32 GB | 128 GB (unified) | BAM sorting is memory-intensive |
| **Storage** | 500 GB SSD | 1 TB NVMe | Fast I/O critical for BAM files |
| **CPU** | 8 cores | ARM64 cores | Parallel I/O and preprocessing |

### Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| **Operating System** | Ubuntu 20.04+ / RHEL 8+ | Host OS |
| **Docker** | 20.10+ | Container runtime |
| **NVIDIA Driver** | 525+ | GPU driver |
| **nvidia-container-toolkit** | Latest | GPU container support |
| **Python** | 3.11+ | Web portal (optional) |
| **aria2c** | Latest | Parallel downloads (optional) |

### NGC Account (Required)

An NVIDIA NGC account is required to pull the Parabricks container:

1. Sign up at [ngc.nvidia.com](https://ngc.nvidia.com)
2. Generate an API key at [ngc.nvidia.com/setup/api-key](https://ngc.nvidia.com/setup/api-key)
3. The API key is used as your password (username is `$oauthtoken`)

---

## Quick Start

### Option 1: Command Line (Recommended)

```bash
# Clone the repository
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd genomics-pipeline

# Run the complete workflow
./run.sh check      # Verify prerequisites (Docker, GPU, disk space)
./run.sh login      # Authenticate with NGC
./run.sh download   # Download GIAB HG002 data (~200GB, 2-6 hours)
./run.sh reference  # Setup GRCh38 reference genome
./run.sh test       # Quick validation on chr20 (5-20 min)
./run.sh full       # Full genome analysis (120-240 min)
```

### Option 2: Web Portal

```bash
# Start the web portal
cd web-portal
./start-portal.sh

# Open browser to http://localhost:5000
# Click through each step in the visual interface
```

### Option 3: Using Your Own Data

```bash
# Copy your FASTQ files to the input directory
cp your_sample_R1.fastq.gz data/input/HG002_R1.fastq.gz
cp your_sample_R2.fastq.gz data/input/HG002_R2.fastq.gz

# Skip download, run the rest
./run.sh check
./run.sh login
./run.sh reference
./run.sh full
```

---

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd genomics-pipeline
```

### Step 2: Verify Prerequisites

```bash
./run.sh check
```

This verifies:
- Docker installation and daemon status
- NVIDIA Container Runtime configuration
- GPU availability and driver version
- Available disk space (minimum 500GB)

**Expected output:**
```
[CHECK] Docker installation... OK
[CHECK] Docker daemon running... OK
[CHECK] NVIDIA Container Runtime... OK
[CHECK] GPU detected: NVIDIA GB10 (128GB)
[CHECK] Driver version: 560.35.03
[CHECK] Available disk space: 1.2TB
[CHECK] All prerequisites met!
```

### Step 3: Authenticate with NGC

```bash
./run.sh login
```

When prompted:
- **Username:** `$oauthtoken` (enter this literally)
- **Password:** Your NGC API key (from ngc.nvidia.com)

### Step 4: Download Data

```bash
./run.sh download
```

This downloads the GIAB HG002 Illumina 2x250bp dataset:
- **Source:** NCBI GIAB FTP server
- **Size:** ~200GB (multiple lane files)
- **Features:** Parallel downloads, automatic resume on interruption

**Note:** Skip this step if using your own FASTQ files.

### Step 5: Setup Reference Genome

```bash
./run.sh reference
```

Downloads and indexes the GRCh38 human reference genome:
- Reference FASTA (GRCh38.fa) - 3.1 GB
- BWA index files (.bwt, .pac, .sa, .amb, .ann)
- FASTA index (.fai)
- Sequence dictionary (.dict)

---

## Usage

### Command Line Interface

The `run.sh` script provides a unified interface:

```bash
./run.sh <command>

Commands:
  check       Check prerequisites (Docker, GPU, disk space)
  login       Authenticate with NGC
  download    Download GIAB HG002 data (~200GB)
  reference   Setup GRCh38 reference genome
  test        Run chr20 fast test (~5-20 min)
  full        Run full genome pipeline (~120-240 min)
  clean       Clean output files (keeps input data)
  clean-all   Clean everything including downloaded data
  help        Show help message
```

### Running the Chr20 Test (Validation)

Before running the full genome, validate your setup with chromosome 20:

```bash
./run.sh test
```

**What this does:**
- Processes only chromosome 20 (~2% of genome)
- Validates GPU acceleration is working
- Completes in 5-20 minutes
- Produces: `HG002.chr20.bam`, `HG002.chr20.vcf.gz`

### Running Full Genome Analysis

After validation:

```bash
./run.sh full
```

**What this does:**
- Processes all chromosomes (chr1-22, chrX, chrY, chrM)
- Takes 120-240 minutes depending on GPU
- Produces: `HG002.genome.bam`, `HG002.genome.vcf.gz`

### Web Portal

The web portal provides visual pipeline management:

```bash
cd web-portal
./start-portal.sh
# Open http://localhost:5000
```

**Portal Features:**
- Click-to-run buttons for each pipeline step
- Real-time console output streaming
- Live GPU utilization monitoring
- Configuration management UI
- Historical log browser

---

## Pipeline Steps in Detail

### Step 0: Prerequisites Check

**Script:** `scripts/00-setup-check.sh`

Validates all system requirements before running the pipeline:

| Check | Requirement | Why It Matters |
|-------|-------------|----------------|
| Docker | Installed + running | Container runtime |
| NVIDIA Runtime | Configured | GPU access in containers |
| GPU | Detected | Acceleration |
| Driver | 525+ | CUDA compatibility |
| Disk Space | 500GB+ | BAM files are large |

### Step 1: NGC Authentication

**Script:** `scripts/01-ngc-login.sh`

Authenticates with NVIDIA NGC to pull the Parabricks container:

```bash
docker login nvcr.io
# Username: $oauthtoken
# Password: <your NGC API key>
```

### Step 2: Data Download

**Script:** `scripts/02-download-data.sh`

Downloads GIAB HG002 benchmark genome:

| File | Size | Description |
|------|------|-------------|
| HG002_R1.fastq.gz | ~100GB | Forward reads (Read 1) |
| HG002_R2.fastq.gz | ~100GB | Reverse reads (Read 2) |

**Features:**
- Parallel downloads using aria2c (8 connections)
- Automatic resume on network interruption
- Checksum validation
- Lane file merging

### Step 3: Reference Genome Setup

**Script:** `scripts/03-setup-reference.sh`

Downloads and prepares the GRCh38 human reference:

| File | Size | Purpose |
|------|------|---------|
| GRCh38.fa | 3.1 GB | Reference sequence |
| GRCh38.fa.fai | 3 KB | FASTA index |
| GRCh38.fa.bwt | 3 GB | BWA index |
| GRCh38.dict | 50 KB | Sequence dictionary |

### Step 4: Chr20 Test

**Script:** `scripts/04-run-chr20-test.sh`

Quick validation run using only chromosome 20:

```bash
# Alignment (fq2bam)
pbrun fq2bam \
  --ref data/ref/GRCh38.fa \
  --in-fq data/input/HG002_R1.fastq.gz data/input/HG002_R2.fastq.gz \
  --out-bam data/output/HG002.chr20.bam \
  --interval chr20

# Variant Calling (DeepVariant)
pbrun deepvariant \
  --ref data/ref/GRCh38.fa \
  --in-bam data/output/HG002.chr20.bam \
  --out-vcf data/output/HG002.chr20.vcf.gz
```

### Step 5: Full Genome Analysis

**Script:** `scripts/05-run-full-genome.sh`

Complete whole-genome analysis:

**Sub-step 5.1: fq2bam (Alignment)**
```bash
pbrun fq2bam \
  --ref data/ref/GRCh38.fa \
  --in-fq data/input/HG002_R1.fastq.gz data/input/HG002_R2.fastq.gz \
  --out-bam data/output/HG002.genome.bam \
  --num-gpus 1
```

What fq2bam does:
1. **Alignment**: Maps reads to reference using BWA-MEM2 (GPU)
2. **Sorting**: Orders reads by genomic coordinate
3. **Duplicate Marking**: Flags PCR duplicates

**Sub-step 5.2: BAM Indexing**
```bash
samtools index data/output/HG002.genome.bam
samtools flagstat data/output/HG002.genome.bam
```

**Sub-step 5.3: DeepVariant (Variant Calling)**
```bash
pbrun deepvariant \
  --ref data/ref/GRCh38.fa \
  --in-bam data/output/HG002.genome.bam \
  --out-vcf data/output/HG002.genome.vcf.gz \
  --num-gpus 1
```

What DeepVariant does:
1. **Pileup Images**: Creates images of read alignments at each position
2. **CNN Inference**: Deep learning model classifies variants (GPU)
3. **Variant Calling**: Outputs high-confidence variant calls

---

## Understanding the Output VCF

### What is a VCF?

A VCF (Variant Call Format) file is a structured summary of how one genome differs from the human reference. Instead of storing every DNA letter (3 billion bases), the VCF records only the meaningful differences.

### VCF Statistics (HG002 Full Genome)

| Metric | Value | Notes |
|--------|-------|-------|
| **Total Variants** | ~11.7 million | All chromosomes |
| **SNPs** | ~4.2 million | Single nucleotide changes |
| **Indels** | ~1.0 million | Insertions/deletions |
| **High Quality (QUAL>30)** | ~3.5 million | Confident calls |
| **In Coding Regions** | ~35,000 | Potentially functional |

### VCF Format Example

```
#CHROM  POS       ID           REF  ALT  QUAL   FILTER  INFO                    FORMAT  HG002
chr9    35065263  rs188935092  G    A    45.2   PASS    DP=32;AF=0.5           GT:DP   0/1:32
chr17   7674220   rs1042522    G    C    99.0   PASS    DP=45;AF=0.5           GT:DP   0/1:45
```

| Column | Meaning |
|--------|---------|
| CHROM | Chromosome (chr9) |
| POS | Position on chromosome (35065263) |
| ID | dbSNP identifier (rs188935092) |
| REF | Reference allele (G) |
| ALT | Alternate allele (A) |
| QUAL | Quality score (higher = more confident) |
| GT | Genotype (0/1 = heterozygous) |

### What Happens Next

The VCF file feeds into the **RAG/Chat Pipeline** (Stage 2), where:
1. Variants are annotated with ClinVar, AlphaMissense, VEP
2. Embedded into vector database (Milvus)
3. Connected to knowledge graph (Clinker)
4. Queried via natural language with Claude

---

## Configuration

### Configuration File

Edit `config/pipeline.env` to customize behavior:

```bash
# GPU Configuration
NUM_GPUS=1              # Number of GPUs to use (1-8)
LOW_MEMORY=0            # Set to 1 for GPUs with <16GB VRAM

# Sample Configuration
PATIENT_ID=HG002        # Sample identifier in output filenames

# Reference Genome
REF_BUILD=GRCh38        # Reference genome build

# Container Image
PB_IMG=nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1

# Performance Tuning
FQ2BAM_ARGS=""          # Additional fq2bam arguments
DEEPVARIANT_ARGS=""     # Additional DeepVariant arguments
```

### Configuration Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NUM_GPUS` | 1 | GPUs for parallel processing |
| `LOW_MEMORY` | 0 | Enable for GPUs with <16GB VRAM |
| `PATIENT_ID` | HG002 | Sample ID in output filenames |
| `PB_IMG` | 4.6.0-1 | Parabricks container version |

### Multi-GPU Configuration

For systems with multiple GPUs:

```bash
# config/pipeline.env
NUM_GPUS=4  # Use 4 GPUs in parallel
```

This distributes the workload across GPUs, reducing processing time proportionally.

---

## Directory Structure

```
genomics-pipeline/
├── run.sh                      # Main CLI interface
├── README.md                   # This documentation
├── QUICKSTART.md               # Quick start guide
├── WEB-PORTAL-GUIDE.md         # Web portal documentation
├── .gitignore                  # Git ignore patterns
│
├── config/
│   └── pipeline.env            # Pipeline configuration
│
├── scripts/
│   ├── 00-setup-check.sh       # Prerequisites validation
│   ├── 01-ngc-login.sh         # NGC authentication
│   ├── 02-download-data.sh     # GIAB data download
│   ├── 03-setup-reference.sh   # Reference genome setup
│   ├── 04-run-chr20-test.sh    # Chr20 validation test
│   └── 05-run-full-genome.sh   # Full genome pipeline
│
├── data/
│   ├── input/                  # Input FASTQ files
│   │   ├── HG002_R1.fastq.gz   # Forward reads (~100GB)
│   │   └── HG002_R2.fastq.gz   # Reverse reads (~100GB)
│   ├── ref/                    # Reference genome
│   │   ├── GRCh38.fa           # Reference FASTA (3.1GB)
│   │   ├── GRCh38.fa.fai       # FASTA index
│   │   ├── GRCh38.fa.bwt       # BWA index
│   │   └── GRCh38.dict         # Sequence dictionary
│   └── output/                 # Pipeline outputs
│       ├── logs/               # Execution logs
│       ├── HG002.chr20.bam     # Chr20 test BAM
│       ├── HG002.chr20.vcf.gz  # Chr20 test VCF
│       ├── HG002.genome.bam    # Full genome BAM (~100GB)
│       ├── HG002.genome.bam.bai # BAM index
│       ├── HG002.genome.vcf.gz # Full genome VCF (→ RAG Pipeline)
│       └── HG002.genome.vcf.gz.tbi # VCF index
│
├── web-portal/
│   ├── start-portal.sh         # Portal startup script
│   ├── requirements.txt        # Python dependencies
│   ├── app/
│   │   └── server.py           # Flask backend
│   ├── templates/
│   │   └── index.html          # Main UI
│   └── static/
│       ├── css/style.css       # Styles
│       └── js/app.js           # Frontend logic
│
└── docs/                       # Additional documentation
```

---

## Performance Benchmarks

### Expected Timings by GPU

| GPU Model | VRAM | Chr20 Test | Full Genome | Notes |
|-----------|------|------------|-------------|-------|
| **DGX Spark (GB10)** | 128GB | 5-10 min | 30-60 min | Recommended |
| **NVIDIA A100** | 80GB | 3-8 min | 25-45 min | Data center |
| **NVIDIA A100** | 40GB | 5-10 min | 35-55 min | Data center |
| **NVIDIA V100** | 32GB | 8-15 min | 50-75 min | Older data center |
| **NVIDIA RTX 4090** | 24GB | 6-12 min | 40-60 min | Consumer |
| **NVIDIA RTX 3090** | 24GB | 8-15 min | 50-75 min | Consumer |

*CPU baseline: 32-core Intel Xeon takes 24-48 hours*

### Step-by-Step Timing Breakdown

| Step | Chr20 | Full Genome | % of Total |
|------|-------|-------------|------------|
| fq2bam (alignment) | 2-8 min | 20-45 min | ~60% |
| BAM indexing | <1 min | 2-5 min | ~5% |
| DeepVariant | 2-6 min | 15-35 min | ~35% |
| **Total** | **5-15 min** | **37-85 min** | 100% |

### Resource Utilization

| Resource | During fq2bam | During DeepVariant |
|----------|---------------|-------------------|
| GPU Compute | 70-90% | 80-95% |
| GPU Memory | 8-12 GB | 12-20 GB |
| System RAM | 24-48 GB | 16-32 GB |
| Disk I/O | High (write) | Moderate (read) |

---

## Output Files

### Primary Outputs

| File | Description | Size | Next Step |
|------|-------------|------|-----------|
| `HG002.genome.bam` | Aligned, sorted, deduplicated reads | ~100 GB | Archive |
| `HG002.genome.bam.bai` | BAM index for random access | ~10 MB | With BAM |
| **`HG002.genome.vcf.gz`** | **Variant calls (main output)** | **~1-2 GB** | **→ RAG Pipeline** |
| `HG002.genome.vcf.gz.tbi` | VCF tabix index | ~2 MB | With VCF |

### Log Files

Located in `data/output/logs/`:

| Log | Contents |
|-----|----------|
| `genome_fq2bam.log` | Alignment metrics, timing |
| `genome_flagstat.log` | Alignment quality statistics |
| `genome_deepvariant.log` | Variant calling statistics |

### VCF Downstream Usage

```bash
# Count total variants
bcftools view -H data/output/HG002.genome.vcf.gz | wc -l
# Output: ~11,700,000

# Filter high-quality variants (QUAL > 30)
bcftools filter -i 'QUAL>30' data/output/HG002.genome.vcf.gz | wc -l
# Output: ~3,500,000

# Extract specific gene region (TP53)
bcftools view -r chr17:7668402-7687550 data/output/HG002.genome.vcf.gz

# Copy to RAG Pipeline
cp data/output/HG002.genome.vcf.gz* ../rag-chat-pipeline/data/input/
```

---

## Troubleshooting

### Common Issues

#### GPU Out of Memory

**Symptom:** `CUDA out of memory` error during fq2bam or DeepVariant

**Solution:**
```bash
# Enable low-memory mode
# Edit config/pipeline.env:
LOW_MEMORY=1

# Re-run pipeline
./run.sh full
```

#### Docker Permission Denied

**Symptom:** `Got permission denied while trying to connect to the Docker daemon`

**Solution:**
```bash
# Add user to docker group
sudo usermod -aG docker $USER

# Log out and back in, or:
newgrp docker
```

#### NVIDIA Container Runtime Not Found

**Symptom:** `could not select device driver "nvidia"`

**Solution:**
```bash
# Install nvidia-container-toolkit
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
```

#### NGC Authentication Failed

**Symptom:** `unauthorized: authentication required`

**Solution:**
```bash
# Verify credentials
docker logout nvcr.io
docker login nvcr.io
# Username: $oauthtoken (literally)
# Password: <your NGC API key>
```

#### Download Interrupted

**Symptom:** Partial FASTQ files after network interruption

**Solution:**
```bash
# aria2c automatically resumes - just re-run
./run.sh download
```

### Checking Logs

```bash
# Real-time monitoring
tail -f data/output/logs/genome_fq2bam.log

# GPU utilization during processing
watch -n 1 nvidia-smi

# Check all logs
ls -lh data/output/logs/
```

---

## Related Pipelines

This pipeline is **Stage 1** of the Precision Medicine to Drug Discovery AI Factory:

| Stage | Pipeline | Repository | Description |
|-------|----------|------------|-------------|
| **1** | **Genomics Pipeline** | [genomics-pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/genomics-pipeline) | FASTQ → VCF (This component) |
| **2** | RAG/Chat Pipeline | [rag-chat-pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/rag-chat-pipeline) | VCF → Target Hypothesis |
| **3** | Drug Discovery Pipeline | [drug-discovery-pipeline](https://github.com/ajones1923/hcls-ai-factory/tree/main/drug-discovery-pipeline) | Target → Molecule Candidates |

### Complete Demo Flow

```
Stage 1 (This Pipeline)
    │
    │  FASTQ files (200GB raw sequencing data)
    │        │
    │        ▼
    │  Parabricks fq2bam + DeepVariant
    │        │
    │        ▼
    │  VCF file (11.7M variants)
    │
    └───────────────────────────────────────────▶ Stage 2 (RAG/Chat)
                                                      │
                                                      │  Annotate with ClinVar,
                                                      │  AlphaMissense, VEP
                                                      │        │
                                                      │        ▼
                                                      │  Vector DB (Milvus)
                                                      │        │
                                                      │        ▼
                                                      │  Claude RAG + Clinker
                                                      │        │
                                                      │        ▼
                                                      │  "VCP is a druggable
                                                      │   target for FTD"
                                                      │
                                                      └──────▶ Stage 3 (Drug Discovery)
                                                                    │
                                                                    ▼
                                                              Generate drug
                                                              candidates for VCP
```

---

## References

### Documentation

- [NVIDIA Parabricks Documentation](https://docs.nvidia.com/clara/parabricks/)
- [BWA-MEM2 Algorithm](https://github.com/bwa-mem2/bwa-mem2)
- [DeepVariant Documentation](https://github.com/google/deepvariant)
- [VCF Specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)

### Data Sources

- [Genome in a Bottle (GIAB)](https://www.nist.gov/programs-projects/genome-bottle)
- [GRCh38 Reference Genome](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)
- [GIAB HG002 Dataset](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/)

### Related Projects

- [NVIDIA DGX Spark](https://www.nvidia.com/en-us/data-center/dgx-spark/)
- [NVIDIA NGC Catalog](https://catalog.ngc.nvidia.com/)
- [NVIDIA Clara](https://www.nvidia.com/en-us/clara/)

---

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](../../LICENSE) file for details.

---

## Acknowledgments

- **NVIDIA** for Parabricks and GPU-accelerated bioinformatics
- **NIST/GIAB** for the HG002 benchmark genome dataset
- **Google Health** for the DeepVariant variant caller
- **Broad Institute** for samtools and bcftools

---

**Note:** This pipeline uses public, de-identified genomics data (GIAB HG002) for demonstration and validation purposes. For clinical use, ensure compliance with relevant regulations and institutional guidelines.
