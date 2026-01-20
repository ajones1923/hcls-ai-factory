# HCLS AI Factory - Environment Guide

## Overview

The HCLS AI Factory is an end-to-end precision medicine platform that transforms patient DNA into drug candidates. The project is maintained in two environments:

1. **Development (Private)** - Where active development happens
2. **Public (Apache 2.0)** - The open-source release on GitHub

---

## Repository Structure

### Development Repositories (Private)

These are the active development repos where code changes are made:

| Repository | Local Path | GitHub URL |
|------------|------------|------------|
| **Genomics Pipeline** | `/home/adam/transfer/genomics-pipeline` | `https://github.com/ajones1923/genomics-pipeline-dev` |
| **RAG/Chat Pipeline** | `/home/adam/transfer/rag-chat-pipeline` | `https://github.com/ajones1923/rag-chat-pipeline-dev` |
| **Drug Discovery Pipeline** | `/home/adam/transfer/drug-discovery-pipeline` | `https://github.com/ajones1923/drug-discovery-pipeline-dev` |
| **Landing Page** | `/home/adam/transfer/landing-page` | `https://github.com/ajones1923/hcls-landing-page-dev` |
| **HLS Orchestrator** | `/home/adam/transfer/hls-orchestrator` | (local only) |
| **HCLS AI Factory (umbrella)** | `/home/adam/transfer/hcls-ai-factory` | `https://github.com/ajones1923/hcls-ai-factory-dev` |

### Public Repository (Apache 2.0)

The open-source release that combines all pipelines into a single repo:

| Repository | Local Path | GitHub URL |
|------------|------------|------------|
| **HCLS AI Factory** | `/home/adam/transfer/hcls-ai-factory-public` | `https://github.com/ajones1923/hcls-ai-factory` |

---

## Directory Structure

### Development Environment

```
/home/adam/transfer/
├── genomics-pipeline/           # Stage 1: FASTQ → VCF (Parabricks, DeepVariant) ~120-240 min
│   ├── scripts/                 # Pipeline scripts (00-setup through 05-run)
│   ├── config/                  # Pipeline configuration
│   ├── web-portal/              # Flask web interface
│   ├── data/                    # Input/output data (gitignored)
│   └── docs/                    # Pipeline documentation
│
├── rag-chat-pipeline/           # Stage 2: Evidence & AI Reasoning
│   ├── app/                     # Streamlit chat UI
│   ├── src/                     # Core modules (RAG engine, LLM client, etc.)
│   ├── data/                    # Knowledge base, targets
│   ├── portal/                  # Web portal
│   └── docker-compose.yml       # Milvus setup
│
├── drug-discovery-pipeline/     # Stage 3: Molecule Generation
│   ├── app/                     # Streamlit discovery UI
│   ├── src/                     # BioNeMo clients, molecule generator
│   ├── data/structures/         # PDB cache, structure images
│   ├── outputs/                 # Generated reports, molecules
│   └── monitoring/              # Grafana/Prometheus setup
│
├── landing-page/                # Main web entry point
│   ├── server.py                # Flask server
│   ├── templates/               # HTML templates
│   └── static/                  # CSS, JS, assets
│
├── hls-orchestrator/            # Nextflow orchestration
│   ├── main.nf                  # Main Nextflow workflow
│   ├── modules/                 # Pipeline modules
│   └── portal/                  # Orchestration portal
│
├── hcls-ai-factory/             # Umbrella repo (docs, scripts)
│   ├── start-services.sh        # Service startup script
│   ├── demo.sh                  # One-click demo launcher
│   └── *.md, *.txt              # Documentation files
│
├── hcls-ai-factory-public/      # PUBLIC RELEASE (synced from dev)
│   ├── genomics-pipeline/
│   ├── rag-chat-pipeline/
│   ├── drug-discovery-pipeline/
│   ├── hls-orchestrator/
│   ├── landing-page/
│   ├── docs/                    # Comprehensive documentation
│   ├── README.md                # Public readme with origin story
│   ├── LICENSE                  # Apache 2.0
│   └── ...
│
└── review/
    └── github-release/          # Release templates and sync tools
        ├── README.md
        ├── LICENSE
        ├── CONTRIBUTING.md
        ├── .gitignore
        ├── .env.example
        ├── docs/                # Documentation templates
        ├── scripts/
        │   ├── check-secrets.sh     # Pre-push secret scanner
        │   └── sync-to-public.sh    # Sync script
        └── ...
```

---

## The Sync Process

### How It Works

The sync script (`/home/adam/transfer/review/github-release/scripts/sync-to-public.sh`) copies code from the development repos to the public repo, excluding sensitive files and large data.

```
DEVELOPMENT                              PUBLIC
───────────────────                      ─────────────────────

genomics-pipeline/      ─┐
rag-chat-pipeline/       │               hcls-ai-factory-public/
drug-discovery-pipeline/ ├── rsync ──►   ├── genomics-pipeline/
hls-orchestrator/        │               ├── rag-chat-pipeline/
landing-page/           ─┘               ├── drug-discovery-pipeline/
                                         ├── hls-orchestrator/
github-release/         ─── copy ───►    ├── landing-page/
(docs, LICENSE, etc)                     ├── docs/
                                         ├── README.md
                                         └── LICENSE
```

### What Gets Excluded

The sync automatically excludes:

| Excluded | Reason |
|----------|--------|
| `.git/` | Each repo has its own git history |
| `.env` | Contains real API keys |
| `venv/` | Virtual environments are local |
| `data/` | Large data files (FASTQ, BAM, VCF) |
| `__pycache__/` | Python bytecode |
| `milvus_data/` | Vector database storage |
| `*.log` | Log files |
| `*.bam`, `*.vcf`, etc. | Large genomic files |

### Running the Sync

```bash
# Dry run (see what would happen)
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh --dry-run

# Sync and commit locally
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh

# Sync, commit, and push to GitHub
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh --push
```

### Sync Script Location

```
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh
```

### Configuration (in sync script)

```bash
# Public repo clone location
PUBLIC_REPO_PATH="/home/adam/transfer/hcls-ai-factory-public"

# Development repo locations
DEV_GENOMICS="/home/adam/transfer/genomics-pipeline"
DEV_RAG_CHAT="/home/adam/transfer/rag-chat-pipeline"
DEV_DRUG_DISCOVERY="/home/adam/transfer/drug-discovery-pipeline"
DEV_HLS_ORCHESTRATOR="/home/adam/transfer/hls-orchestrator"
DEV_LANDING_PAGE="/home/adam/transfer/landing-page"

# Documentation/templates location
RELEASE_DOCS="/home/adam/transfer/review/github-release"
```

---

## Credential Management

### API Keys Required

| Service | Purpose | Where to Get |
|---------|---------|--------------|
| **NVIDIA NGC** | Parabricks, BioNeMo NIMs | https://ngc.nvidia.com/setup/api-key |
| **Anthropic** | Claude AI for RAG pipeline | https://console.anthropic.com/ |
| **HuggingFace** | Local LLM models (optional) | https://huggingface.co/settings/tokens |

### Storage

- Development: Keys stored in `.env` files (gitignored)
- Public: Only `.env.example` with placeholders

### Secret Scanner

Before pushing to public, run:

```bash
/home/adam/transfer/review/github-release/scripts/check-secrets.sh
```

This scans for accidentally committed credentials.

---

## Service Ports

| Service | Port | Description |
|---------|------|-------------|
| Landing Page | 8080 | Main entry point |
| Genomics Portal | 5000 | Genomics pipeline web UI |
| RAG Portal | 5001 | RAG/Chat API portal |
| RAG Chat UI | 8501 | Streamlit chat interface |
| Drug Discovery UI | 8505 | Molecule generation interface |
| Drug Discovery Portal | 8510 | Pipeline orchestration |
| Milvus | 19530 | Vector database |
| Grafana | 3000 | Monitoring dashboards |

---

## Starting Services

### Quick Start (All Services)

```bash
cd /home/adam/transfer/hcls-ai-factory
./start-services.sh
```

### Demo Mode

```bash
cd /home/adam/transfer/hcls-ai-factory
./demo.sh
```

### Individual Pipelines

```bash
# Start only RAG/Chat
./start-services.sh --rag

# Start only Drug Discovery
./start-services.sh --drug

# Check status
./start-services.sh --status

# Stop all
./start-services.sh --stop
```

---

## Git Workflow

### Development Work

```bash
# Work in individual pipeline repos
cd /home/adam/transfer/genomics-pipeline
# Make changes...
git add .
git commit -m "Description of changes"
git push origin main
```

### Publishing to Public

```bash
# 1. Ensure all dev repos are committed
# 2. Run sync script
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh --push

# This will:
# - Copy code from dev repos to public repo
# - Exclude sensitive files and data
# - Run secret scanner
# - Commit and push to GitHub
```

---

## Key Files

### Public Repo Root Files

| File | Purpose |
|------|---------|
| `README.md` | Origin story + project overview |
| `LICENSE` | Apache 2.0 (Copyright Adam Jones) |
| `CONTRIBUTING.md` | Contribution guidelines |
| `.gitignore` | Excludes credentials, data, caches |
| `.env.example` | Configuration template |
| `start-services.sh` | Service startup script |
| `demo.sh` | One-click demo launcher |

### Documentation

| File | Purpose |
|------|---------|
| `docs/PRODUCT_DOCUMENTATION.txt` | Complete reference (3,300+ lines) |
| `docs/README.md` | Documentation index |
| `docs/genomics-pipeline/` | Stage 1 documentation |
| `docs/rag-chat-pipeline/` | Stage 2 documentation |
| `docs/drug-discovery-pipeline/` | Stage 3 documentation |
| `docs/diagrams/` | Mermaid architecture diagrams |

---

## Author & Licensing

**Author:** Adam Jones
**LinkedIn:** https://www.linkedin.com/in/socal-engineer
**License:** Apache 2.0

The public repository is released under Apache 2.0, which allows:
- Commercial use
- Modification
- Distribution
- Patent use

With requirements:
- License and copyright notice must be included
- Changes must be documented

---

## Quick Reference

### URLs

| Resource | URL |
|----------|-----|
| Public Repo | https://github.com/ajones1923/hcls-ai-factory |
| Dev Repo (umbrella) | https://github.com/ajones1923/hcls-ai-factory-dev |
| Genomics Dev | https://github.com/ajones1923/genomics-pipeline-dev |
| RAG/Chat Dev | https://github.com/ajones1923/rag-chat-pipeline-dev |
| Drug Discovery Dev | https://github.com/ajones1923/drug-discovery-pipeline-dev |
| Landing Page Dev | https://github.com/ajones1923/hcls-landing-page-dev |

### Local Paths

| Resource | Path |
|----------|------|
| Public Repo | `/home/adam/transfer/hcls-ai-factory-public` |
| Sync Script | `/home/adam/transfer/review/github-release/scripts/sync-to-public.sh` |
| Secret Scanner | `/home/adam/transfer/review/github-release/scripts/check-secrets.sh` |
| Release Templates | `/home/adam/transfer/review/github-release/` |

### Common Commands

```bash
# Sync to public (dry run)
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh --dry-run

# Sync to public (with push)
/home/adam/transfer/review/github-release/scripts/sync-to-public.sh --push

# Start all services
cd /home/adam/transfer/hcls-ai-factory && ./start-services.sh

# Check service status
cd /home/adam/transfer/hcls-ai-factory && ./start-services.sh --status

# Run secret scanner
/home/adam/transfer/review/github-release/scripts/check-secrets.sh
```

---

## Hardware

The platform is developed on:

- **NVIDIA DGX Spark** (personally owned)
- 128GB GPU memory
- 512GB RAM
- Ubuntu 22.04 LTS

---

## Project Background

In 2012, Adam Jones set out to use his HPC/supercomputing background to make a difference in healthcare. Starting with the conviction that no parent should ever lose a child to disease, he focused on Pediatric Neuroblastoma and spent 14+ years learning biology, genomics, and drug discovery.

The HCLS AI Factory is the result—a vendor-neutral, open-source platform released freely under Apache 2.0 for anyone to use, extend, and improve.

---

*Last updated: January 2026*
