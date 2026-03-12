# Precision Oncology Intelligence Agent -- Deployment Guide

**HCLS AI Factory / ai_agent_adds / precision_oncology_agent**

Version 1.0.0 | February 2026 | Author: Adam Jones

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Quick Start](#3-quick-start)
4. [Deployment Modes](#4-deployment-modes)
   - 4a. [Docker Lite (Milvus + API Only)](#4a-docker-lite-milvus--api-only)
   - 4b. [Docker Full Stack (All 6 Services)](#4b-docker-full-stack-all-6-services)
   - 4c. [DGX Spark Production](#4c-dgx-spark-production)
   - 4d. [Development Mode (Local Python)](#4d-development-mode-local-python)
5. [Configuration Reference](#5-configuration-reference)
6. [Collection Setup and Seeding](#6-collection-setup-and-seeding)
7. [Data Ingestion (CIViC, PubMed, ClinicalTrials.gov)](#7-data-ingestion-civic-pubmed-clinicaltrialsgov)
8. [Networking and Ports](#8-networking-and-ports)
9. [Storage and Persistence](#9-storage-and-persistence)
10. [Monitoring and Metrics](#10-monitoring-and-metrics)
11. [Security Hardening](#11-security-hardening)
12. [Health Checks and Troubleshooting](#12-health-checks-and-troubleshooting)
13. [Backup and Recovery](#13-backup-and-recovery)
14. [Scaling Considerations](#14-scaling-considerations)
15. [Integration with HCLS AI Factory](#15-integration-with-hcls-ai-factory)
16. [Updating and Maintenance](#16-updating-and-maintenance)
17. [Appendix A: Complete docker-compose.yml](#appendix-a-complete-docker-composeyml)
18. [Appendix B: Environment Variable Quick Reference](#appendix-b-environment-variable-quick-reference)

---

## 1. Overview

The Precision Oncology Intelligence Agent is a RAG-powered clinical decision
support system designed for molecular tumor boards (MTBs). It combines a
multi-collection Milvus vector database, BGE-small-en-v1.5 sentence embeddings,
and Claude LLM reasoning to deliver evidence-based therapy recommendations,
clinical trial matching, and resistance mechanism analysis.

### Core Capabilities

- **Multi-collection RAG search** across 11 knowledge domains (variants,
  literature, therapies, trials, biomarkers, pathways, guidelines, resistance
  mechanisms, outcomes, patient cases, and genomic evidence).
- **Clinical trial matching** with composite scoring (biomarker 40%, semantic
  25%, phase 20%, status 15%).
- **Therapy ranking** with evidence-level weighting and resistance awareness.
- **Cross-modal analysis** linking genomic, imaging, and clinical data.
- **FHIR-compatible** case management and export.
- **PDF report generation** for tumor board presentations.
- **Prometheus metrics** for operational monitoring.

### Architecture at a Glance

The agent runs as 6 Docker services on a bridge network (`onco-network`):

```
                    onco-network (bridge)
    +--------------------------------------------------+
    |                                                  |
    |  milvus-etcd -----> milvus-standalone <-----+    |
    |  milvus-minio ---/    :19530 :9091          |    |
    |                                             |    |
    |  onco-api (:8527) --------------------------+    |
    |  onco-streamlit (:8526) --------------------+    |
    |  onco-setup (one-shot) ---------------------+    |
    |                                                  |
    +--------------------------------------------------+
```

| Service | Container Name | Purpose |
|---------|---------------|---------|
| `milvus-etcd` | `onco-milvus-etcd` | Metadata key-value store for Milvus |
| `milvus-minio` | `onco-milvus-minio` | Object storage for Milvus log/index data |
| `milvus-standalone` | `onco-milvus-standalone` | Vector database (Milvus 2.4) |
| `onco-streamlit` | `onco-streamlit` | Streamlit clinical UI |
| `onco-api` | `onco-api` | FastAPI REST server |
| `onco-setup` | `onco-setup` | One-shot collection creation + data seeding |

### Software Components

The agent codebase is organized as follows:

```
agent/
  api/                  # FastAPI application (main.py + routes/)
    routes/             # Endpoint modules: meta_agent, cases, trials, reports, events
  app/                  # Streamlit UI (oncology_ui.py)
  config/               # Pydantic settings (settings.py)
  data/reference/       # 10 seed JSON files (~768 KB total)
  scripts/              # Setup, seed, ingest, and validation scripts
  src/                  # Core modules
    agent.py            # OncoIntelligenceAgent orchestrator
    case_manager.py     # FHIR-compatible case management
    collections.py      # 11 Milvus collection schemas + OncoCollectionManager
    cross_modal.py      # Cross-modal analysis trigger
    export.py           # PDF report generation
    knowledge.py        # In-memory knowledge graph (targets, therapies, resistance)
    metrics.py          # Prometheus metric definitions
    models.py           # Pydantic data models
    query_expansion.py  # Query rewriting and expansion
    rag_engine.py       # Multi-collection RAG search engine
    scheduler.py        # APScheduler-based periodic ingestion
    therapy_ranker.py   # Evidence-weighted therapy ranking
    trial_matcher.py    # Clinical trial matching with composite scoring
    ingest/             # Data source parsers (CIViC, OncoKB, PubMed, guidelines, etc.)
    utils/              # Utilities (pubmed_client.py)
    workflows/          # Multi-step agent workflows
  tests/                # Unit tests (6 test modules)
  Dockerfile            # Multi-stage build (Python 3.10-slim)
  docker-compose.yml    # All 6 services
  requirements.txt      # 27 dependency groups (~57 resolved packages)
```

---

## 2. Prerequisites

### 2.1 Hardware Requirements

| Tier | RAM | CPU | Disk | GPU | Use Case |
|------|-----|-----|------|-----|----------|
| **Minimum** | 16 GB | 4 cores | 50 GB | None | Development, demos |
| **Recommended** | 32 GB | 8 cores | 100 GB | NVIDIA GPU (8+ GB VRAM) | Staging, small production |
| **Production** | 128 GB | 20-core Grace CPU | 500 GB NVMe | RTX PRO 6000 | NVIDIA DGX Spark, full pipeline |

**Notes:**
- Milvus standalone requires approximately 4-8 GB RAM depending on index size.
- The BGE-small-en-v1.5 embedding model loads approximately 130 MB into memory.
- GPU is optional but significantly accelerates embedding generation during
  bulk seeding and live ingest operations.
- Disk requirements increase with ingested data volume. The base seed data
  is approximately 768 KB; a production deployment with full PubMed and
  ClinicalTrials.gov ingest may reach 10-50 GB in Milvus storage.

### 2.2 Software Requirements

| Component | Minimum Version | Purpose |
|-----------|----------------|---------|
| Docker Engine | 24.0+ | Container runtime |
| Docker Compose | 2.20+ | Multi-service orchestration |
| Python | 3.10+ | Local development (if not using Docker) |
| Git | 2.30+ | Source code management |
| curl | 7.0+ | Health checks and API testing |

**Verify Docker installation:**

```bash
docker --version          # Docker version 24.0+
docker compose version    # Docker Compose version v2.20+
```

### 2.3 API Keys

| Key | Required | Source | Purpose |
|-----|----------|--------|---------|
| `ANTHROPIC_API_KEY` | **Yes** (for LLM features) | [console.anthropic.com](https://console.anthropic.com) | Claude LLM for RAG answer generation |
| `NCBI_API_KEY` | No (recommended) | [ncbi.nlm.nih.gov/account](https://www.ncbi.nlm.nih.gov/account/) | PubMed E-utilities (higher rate limits) |

**Without `ANTHROPIC_API_KEY`:** Vector search, trial matching, and therapy
ranking still function. Only the LLM-generated narrative answers (the `/query`
endpoint) will fail.

**Without `NCBI_API_KEY`:** PubMed ingest works but is rate-limited to 3
requests/second instead of 10 requests/second.

### 2.4 Network Requirements

- Outbound HTTPS (port 443) access to:
  - `api.anthropic.com` (Claude API)
  - `eutils.ncbi.nlm.nih.gov` (PubMed)
  - `clinicaltrials.gov` (ClinicalTrials.gov API v2)
  - `civicdb.org` (CIViC variant database)
  - `huggingface.co` (initial model download for BGE-small-en-v1.5)
- No inbound ports are required beyond the service ports listed in
  [Section 8](#8-networking-and-ports).

---

## 3. Quick Start

Five commands to get the agent running with Docker Compose:

```bash
# 1. Clone or navigate to the agent directory
cd ai_agent_adds/precision_oncology_agent/agent

# 2. Create your environment file
cp .env.example .env
# Edit .env and set ANTHROPIC_API_KEY=sk-ant-...

# 3. Build the application images
docker compose build

# 4. Start all 6 services
docker compose up -d

# 5. Watch the setup/seed progress
docker compose logs -f onco-setup
```

**If `.env.example` does not exist**, create `.env` manually:

```bash
cat > .env << 'EOF'
ANTHROPIC_API_KEY=sk-ant-your-key-here
EOF
```

### Verify the Deployment

Once `onco-setup` exits with code 0 (typically 2-5 minutes), verify:

```bash
# Check all services are running
docker compose ps

# Milvus health
curl -s http://localhost:9091/healthz
# Expected: {"status":"OK"}

# API health
curl -s http://localhost:8527/health | python3 -m json.tool
# Expected: {"status": "healthy", "collections": {...}, ...}

# Open the Streamlit UI
open http://localhost:8526    # macOS
xdg-open http://localhost:8526  # Linux
```

### Quick Test Query

```bash
curl -s -X POST http://localhost:8527/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What targeted therapies are available for EGFR L858R in NSCLC?", "top_k": 5}' \
  | python3 -m json.tool
```

---

## 4. Deployment Modes

### 4a. Docker Lite (Milvus + API Only)

This mode runs only the vector database infrastructure and the FastAPI server,
omitting the Streamlit UI. Suitable for headless API-only deployments or when
the UI is served separately.

**Start Lite mode:**

```bash
docker compose up -d milvus-etcd milvus-minio milvus-standalone onco-api onco-setup
```

**Resource footprint:**
- RAM: ~8-12 GB
- CPU: 2-4 cores
- Disk: 20 GB minimum

**Access points:**
- API: `http://localhost:8527`
- API docs (Swagger): `http://localhost:8527/docs`
- Milvus gRPC: `localhost:19530`
- Milvus metrics: `http://localhost:9091/healthz`

**When to use:**
- Backend-only deployments where a custom frontend connects via the REST API.
- CI/CD pipelines that run integration tests against the API.
- Resource-constrained environments that cannot spare memory for Streamlit.

### 4b. Docker Full Stack (All 6 Services)

This is the default mode as described in the [Quick Start](#3-quick-start).
All six services run together on the `onco-network` bridge.

**Start Full Stack:**

```bash
docker compose up -d
```

**Service startup order** (enforced by `depends_on` with health checks):

1. `milvus-etcd` starts first (healthcheck: `etcdctl endpoint health`)
2. `milvus-minio` starts in parallel with etcd (healthcheck: MinIO `/minio/health/live`)
3. `milvus-standalone` waits for both etcd and MinIO to be healthy
4. `onco-streamlit`, `onco-api`, and `onco-setup` all wait for Milvus to be healthy
5. `onco-setup` runs once (creates collections, seeds data, then exits with code 0)

**Resource footprint:**
- RAM: ~12-16 GB
- CPU: 4 cores minimum
- Disk: 30-50 GB

**Access points:**
- Streamlit UI: `http://localhost:8526`
- API: `http://localhost:8527`
- API docs (Swagger): `http://localhost:8527/docs`
- MinIO Console: `http://localhost:9001` (credentials: `minioadmin` / `minioadmin`)
- Milvus gRPC: `localhost:19530`
- Milvus metrics/health: `http://localhost:9091/healthz`
- etcd: `localhost:2379`

### 4c. DGX Spark Production

For deployment on an NVIDIA DGX Spark (128 GB RAM, 20-core Grace CPU,
RTX PRO 6000), the agent integrates with the broader HCLS AI Factory pipeline.

#### Production Environment File

Create `/etc/onco-agent/.env`:

```bash
# === Required ===
ANTHROPIC_API_KEY=sk-ant-your-production-key

# === Milvus (use dedicated host if available) ===
ONCO_MILVUS_HOST=milvus-standalone
ONCO_MILVUS_PORT=19530

# === Embedding ===
ONCO_EMBEDDING_MODEL=BAAI/bge-small-en-v1.5
ONCO_EMBEDDING_DIM=384

# === LLM ===
ONCO_LLM_MODEL=claude-sonnet-4-20250514

# === API ===
ONCO_API_PORT=8527
ONCO_API_BASE_URL=http://onco-api:8527
ONCO_CORS_ORIGINS=http://localhost:8080,http://localhost:8526,http://localhost:8527

# === Security ===
ONCO_MAX_REQUEST_SIZE_MB=10

# === Monitoring ===
ONCO_METRICS_ENABLED=true

# === PubMed (recommended for production ingest) ===
NCBI_API_KEY=your-ncbi-key

# === Scheduler (weekly refresh) ===
ONCO_SCHEDULER_INTERVAL=168h
```

#### Production docker-compose Override

Create `docker-compose.prod.yml` to set resource limits and increase workers:

- `milvus-standalone`: memory limit 16G, reservation 8G, log level `warn`
- `onco-api`: memory limit 8G, reservation 4G, `--workers=4`
- `onco-streamlit`: memory limit 4G, reservation 2G

```bash
docker compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```

#### DGX Spark Systemd Service

For automatic startup on boot, create a systemd unit at
`/etc/systemd/system/onco-agent.service` with `Type=oneshot`,
`RemainAfterExit=yes`, `After=docker.service`, and `Requires=docker.service`.
Set `ExecStart` to `docker compose -f docker-compose.yml -f docker-compose.prod.yml up -d`
and `ExecStop` to `docker compose down`. Then:

```bash
sudo systemctl daemon-reload && sudo systemctl enable --now onco-agent
```

### 4d. Development Mode (Local Python)

Run the agent directly with Python for rapid iteration and debugging.

#### Step 1: Set Up Python Environment

```bash
cd ai_agent_adds/precision_oncology_agent/agent

# Create virtual environment
python3.10 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install --upgrade pip
pip install -r requirements.txt
```

#### Step 2: Start Milvus (Docker Required)

Even in development mode, Milvus still runs in Docker:

```bash
docker compose up -d milvus-etcd milvus-minio milvus-standalone
```

Wait for Milvus to be healthy:

```bash
# Poll until healthy
until curl -sf http://localhost:9091/healthz > /dev/null; do
  echo "Waiting for Milvus..."
  sleep 5
done
echo "Milvus is ready."
```

#### Step 3: Create Collections and Seed Data

```bash
# Set environment variables
export ONCO_MILVUS_HOST=localhost
export ONCO_MILVUS_PORT=19530
export ANTHROPIC_API_KEY=sk-ant-your-key

# Create all 11 collections and seed with reference data
python scripts/setup_collections.py --drop-existing --seed
```

Alternatively, run the setup and individual seed scripts separately:

```bash
# Create collections only
python scripts/setup_collections.py --drop-existing

# Seed individually
python scripts/seed_variants.py
python scripts/seed_literature.py
python scripts/seed_trials.py
python scripts/seed_therapies.py
python scripts/seed_biomarkers.py
python scripts/seed_pathways.py
python scripts/seed_guidelines.py
python scripts/seed_resistance.py
python scripts/seed_outcomes.py
python scripts/seed_cases.py
python scripts/seed_knowledge.py
```

#### Step 4: Run the API Server

```bash
# Option A: Direct uvicorn
uvicorn api.main:app --host 0.0.0.0 --port 8527 --reload

# Option B: Python module
python -m api.main
```

#### Step 5: Run the Streamlit UI

In a separate terminal:

```bash
source .venv/bin/activate
export ONCO_MILVUS_HOST=localhost
export ONCO_API_BASE_URL=http://localhost:8527
streamlit run app/oncology_ui.py --server.port 8526
```

#### Step 6: Run Tests

```bash
# All tests
pytest tests/ -v

# Individual test modules
pytest tests/test_collections.py -v
pytest tests/test_agent.py -v
pytest tests/test_case_manager.py -v
pytest tests/test_trial_matcher.py -v
pytest tests/test_therapy_ranker.py -v
pytest tests/test_knowledge.py -v
```

---

## 5. Configuration Reference

All configuration is managed through `config/settings.py` using Pydantic
`BaseSettings`. Every setting can be overridden via environment variables
with the `ONCO_` prefix (except `ANTHROPIC_API_KEY` and `NCBI_API_KEY`,
which use their standard names).

### 5.1 Connection Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_MILVUS_HOST` | `localhost` | Milvus server hostname. Set to `milvus-standalone` in Docker. |
| `ONCO_MILVUS_PORT` | `19530` | Milvus gRPC port. |
| `ONCO_API_HOST` | `0.0.0.0` | FastAPI bind address. |
| `ONCO_API_PORT` | `8527` | FastAPI listen port. |
| `ONCO_STREAMLIT_PORT` | `8526` | Streamlit server port. |
| `ONCO_API_BASE_URL` | `http://localhost:8527` | Base URL for API calls from the Streamlit UI. |

### 5.2 Embedding Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | HuggingFace model for sentence embeddings. |
| `ONCO_EMBEDDING_DIM` | `384` | Embedding vector dimension. Must match the model output. |
| `ONCO_EMBEDDING_BATCH_SIZE` | `32` | Batch size for embedding generation during ingest. |

### 5.3 LLM Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ANTHROPIC_API_KEY` | _(none)_ | **Required** for LLM features. Anthropic API key. |
| `ONCO_LLM_PROVIDER` | `anthropic` | LLM provider identifier. |
| `ONCO_LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model to use for RAG answer generation. |

### 5.4 RAG Search Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_TOP_K` | `5` | Number of top results to retrieve per collection. |
| `ONCO_SCORE_THRESHOLD` | `0.4` | Minimum cosine similarity score for result inclusion. |
| `ONCO_MIN_SUFFICIENT_HITS` | `3` | Minimum hits before the agent considers evidence sufficient. |
| `ONCO_MIN_COLLECTIONS_FOR_SUFFICIENT` | `2` | Minimum collections with hits for sufficient evidence. |
| `ONCO_MIN_SIMILARITY_SCORE` | `0.30` | Absolute minimum similarity for any hit. |

### 5.5 Collection Weight Settings

These weights control the relative importance of each collection in the
multi-collection RAG search. They should sum to approximately 1.0.

| Variable | Default | Collection |
|----------|---------|------------|
| `ONCO_WEIGHT_VARIANTS` | `0.18` | `onco_variants` |
| `ONCO_WEIGHT_LITERATURE` | `0.16` | `onco_literature` |
| `ONCO_WEIGHT_THERAPIES` | `0.14` | `onco_therapies` |
| `ONCO_WEIGHT_GUIDELINES` | `0.12` | `onco_guidelines` |
| `ONCO_WEIGHT_TRIALS` | `0.10` | `onco_trials` |
| `ONCO_WEIGHT_BIOMARKERS` | `0.08` | `onco_biomarkers` |
| `ONCO_WEIGHT_RESISTANCE` | `0.07` | `onco_resistance` |
| `ONCO_WEIGHT_PATHWAYS` | `0.06` | `onco_pathways` |
| `ONCO_WEIGHT_OUTCOMES` | `0.04` | `onco_outcomes` |
| `ONCO_WEIGHT_CASES` | `0.02` | `onco_cases` |
| `ONCO_WEIGHT_GENOMIC` | `0.03` | `genomic_evidence` |

### 5.6 Trial Matching Weights

These weights control the composite scoring for clinical trial matching.

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_TRIAL_WEIGHT_BIOMARKER` | `0.40` | Weight for biomarker match score. |
| `ONCO_TRIAL_WEIGHT_SEMANTIC` | `0.25` | Weight for semantic similarity score. |
| `ONCO_TRIAL_WEIGHT_PHASE` | `0.20` | Weight for trial phase preference (Phase 3 > Phase 1). |
| `ONCO_TRIAL_WEIGHT_STATUS` | `0.15` | Weight for trial recruitment status. |

### 5.7 Citation Thresholds

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_CITATION_STRONG_THRESHOLD` | `0.75` | Similarity score above which a citation is labeled "strong evidence." |
| `ONCO_CITATION_MODERATE_THRESHOLD` | `0.60` | Similarity score above which a citation is labeled "moderate evidence." |

### 5.8 Cross-Modal Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_CROSS_MODAL_ENABLED` | `true` | Enable cross-modal analysis (genomic + imaging linking). |
| `ONCO_CROSS_MODAL_THRESHOLD` | `0.40` | Minimum score to trigger cross-modal expansion. |
| `ONCO_GENOMIC_TOP_K` | `5` | Number of genomic evidence results to retrieve. |
| `ONCO_IMAGING_TOP_K` | `5` | Number of imaging results to retrieve. |

### 5.9 External API Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_API_KEY` | _(none)_ | NCBI E-utilities API key for PubMed ingest. Optional but recommended. |
| `ONCO_PUBMED_MAX_RESULTS` | `5000` | Maximum PubMed articles to fetch per ingest run. |
| `ONCO_CT_GOV_BASE_URL` | `https://clinicaltrials.gov/api/v2` | ClinicalTrials.gov API v2 base URL. |
| `ONCO_CIVIC_BASE_URL` | `https://civicdb.org/api` | CIViC database API base URL. |

### 5.10 Operational Settings

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_METRICS_ENABLED` | `true` | Enable Prometheus metrics collection. |
| `ONCO_SCHEDULER_INTERVAL` | `168h` | Interval for scheduled data refresh (default: 1 week). |
| `ONCO_CORS_ORIGINS` | `http://localhost:8080,http://localhost:8526,http://localhost:8527` | Comma-separated allowed CORS origins. |
| `ONCO_MAX_REQUEST_SIZE_MB` | `10` | Maximum HTTP request body size in megabytes. |
| `ONCO_CONVERSATION_MEMORY_DEPTH` | `3` | Number of conversation turns to retain for context. |

### 5.11 PDF Report Branding

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_PDF_BRAND_COLOR_R` | `118` | Brand color red component (0-255). |
| `ONCO_PDF_BRAND_COLOR_G` | `185` | Brand color green component (0-255). |
| `ONCO_PDF_BRAND_COLOR_B` | `0` | Brand color blue component (0-255). |

### 5.12 Collection Names

These are rarely changed but can be overridden if you need custom naming:

| Variable | Default |
|----------|---------|
| `ONCO_COLLECTION_LITERATURE` | `onco_literature` |
| `ONCO_COLLECTION_TRIALS` | `onco_trials` |
| `ONCO_COLLECTION_VARIANTS` | `onco_variants` |
| `ONCO_COLLECTION_BIOMARKERS` | `onco_biomarkers` |
| `ONCO_COLLECTION_THERAPIES` | `onco_therapies` |
| `ONCO_COLLECTION_PATHWAYS` | `onco_pathways` |
| `ONCO_COLLECTION_GUIDELINES` | `onco_guidelines` |
| `ONCO_COLLECTION_RESISTANCE` | `onco_resistance` |
| `ONCO_COLLECTION_OUTCOMES` | `onco_outcomes` |
| `ONCO_COLLECTION_CASES` | `onco_cases` |
| `ONCO_COLLECTION_GENOMIC` | `genomic_evidence` |

---

## 6. Collection Setup and Seeding

### 6.1 Collection Schemas

The agent manages 11 Milvus collections. Each uses IVF_FLAT indexing with
COSINE similarity on 384-dimensional BGE-small-en-v1.5 embeddings.

| # | Collection | Content | Key Fields |
|---|-----------|---------|------------|
| 1 | `onco_variants` | Actionable somatic/germline variants (CIViC, OncoKB) | gene, variant, cancer_type, evidence_level |
| 2 | `onco_literature` | PubMed/PMC/preprint chunks by cancer type | title, text_chunk, source_type, year, journal |
| 3 | `onco_therapies` | Approved and investigational therapies with MOA | drug_name, mechanism, cancer_type, approval_status |
| 4 | `onco_guidelines` | NCCN/ASCO/ESMO guideline recommendations | guideline_body, recommendation, cancer_type |
| 5 | `onco_trials` | ClinicalTrials.gov summaries with biomarker criteria | nct_id, title, phase, status, biomarkers |
| 6 | `onco_biomarkers` | Predictive/prognostic biomarkers and assays | biomarker_name, test_type, cancer_type |
| 7 | `onco_resistance` | Resistance mechanisms and bypass strategies | mechanism, gene, drug, bypass_strategy |
| 8 | `onco_pathways` | Signaling pathways, cross-talk, druggable nodes | pathway_name, genes, druggable_nodes |
| 9 | `onco_outcomes` | Real-world treatment outcome records | treatment, response, duration, cancer_type |
| 10 | `onco_cases` | De-identified patient case snapshots | diagnosis, mutations, treatment_history |
| 11 | `genomic_evidence` | Read-only VCF-derived evidence from Stage 1 pipeline | gene, variant, consequence, impact |

### 6.2 Setup Script

The `scripts/setup_collections.py` script creates all collection schemas and
optionally seeds them with reference data.

```bash
# Create collections only (preserves existing data)
python scripts/setup_collections.py

# Drop and recreate all collections
python scripts/setup_collections.py --drop-existing

# Drop, recreate, and seed with reference data
python scripts/setup_collections.py --drop-existing --seed

# Connect to a specific Milvus host
python scripts/setup_collections.py --host milvus-standalone --port 19530 --drop-existing --seed
```

### 6.3 Seed Scripts

There are 11 seed scripts, each populating a specific collection from
JSON files in `data/reference/`:

| Script | Source File | Collection |
|--------|-----------|------------|
| `seed_variants.py` | `variant_seed_data.json` | `onco_variants` |
| `seed_literature.py` | `literature_seed_data.json` | `onco_literature` |
| `seed_therapies.py` | `therapy_seed_data.json` | `onco_therapies` |
| `seed_guidelines.py` | `guideline_seed_data.json` | `onco_guidelines` |
| `seed_trials.py` | `trial_seed_data.json` | `onco_trials` |
| `seed_biomarkers.py` | `biomarker_seed_data.json` | `onco_biomarkers` |
| `seed_resistance.py` | `resistance_seed_data.json` | `onco_resistance` |
| `seed_pathways.py` | `pathway_seed_data.json` | `onco_pathways` |
| `seed_outcomes.py` | `outcome_seed_data.json` | `onco_outcomes` |
| `seed_cases.py` | `cases_seed_data.json` | `onco_cases` |
| `seed_knowledge.py` | _(in-memory maps)_ | Loads knowledge graph into module memory |

**Seed data files** are located in `data/reference/` and total approximately
768 KB (10 JSON files).

### 6.4 Docker-Based Seeding

In Docker Compose, the `onco-setup` service runs seeding automatically:

```bash
# Watch seed progress
docker compose logs -f onco-setup

# Re-run seeding (drop collections first)
docker compose rm -f onco-setup
docker compose up onco-setup
```

The `onco-setup` container executes this sequence:
1. `setup_collections.py --drop-existing` (create all 11 schemas)
2. `seed_variants.py` through `seed_outcomes.py` (9 seed scripts in order)

After completion, the container exits with code 0 (restart policy: `"no"`).

### 6.5 Verifying Seed Status

```bash
# Via API
curl -s http://localhost:8527/collections | python3 -m json.tool

# Via API health endpoint
curl -s http://localhost:8527/health | python3 -m json.tool

# Expected output includes collection counts:
# {
#   "status": "healthy",
#   "collections": {
#     "onco_variants": <count>,
#     "onco_literature": <count>,
#     "onco_therapies": <count>,
#     ...
#   },
#   "total_vectors": <total>,
#   ...
# }
```

---

## 7. Data Ingestion (CIViC, PubMed, ClinicalTrials.gov)

Beyond the static seed data, the agent includes 3 live ingest scripts that
fetch and parse data from external APIs.

### 7.1 CIViC Variant Ingest

Fetches actionable variants and evidence from the CIViC (Clinical
Interpretation of Variants in Cancer) database.

```bash
python scripts/ingest_civic.py
```

**Configuration:**
- `ONCO_CIVIC_BASE_URL` (default: `https://civicdb.org/api`)
- Uses `src/ingest/civic_parser.py` to parse CIViC API responses.
- Populates `onco_variants` collection with clinical evidence items.

**What it fetches:**
- Gene-variant-disease associations
- Evidence levels (A through E)
- Drug sensitivity/resistance annotations
- Supporting publications

### 7.2 PubMed Literature Ingest

Fetches recent oncology literature from PubMed via the NCBI E-utilities API.

```bash
# Without API key (3 req/sec rate limit)
python scripts/ingest_pubmed.py

# With API key (10 req/sec rate limit)
NCBI_API_KEY=your-key python scripts/ingest_pubmed.py
```

**Configuration:**
- `NCBI_API_KEY` (optional, recommended for higher rate limits)
- `ONCO_PUBMED_MAX_RESULTS` (default: `5000`)
- Uses `src/ingest/literature_parser.py` and `src/utils/pubmed_client.py`.
- Populates `onco_literature` collection.

**What it fetches:**
- Abstracts and metadata from recent oncology publications
- Chunks text for embedding and indexes by cancer type, gene, and variant

### 7.3 ClinicalTrials.gov Ingest

Fetches active oncology clinical trials from the ClinicalTrials.gov API v2.

```bash
python scripts/ingest_clinical_trials.py
```

**Configuration:**
- `ONCO_CT_GOV_BASE_URL` (default: `https://clinicaltrials.gov/api/v2`)
- Uses `src/ingest/clinical_trials_parser.py`.
- Populates `onco_trials` collection.

**What it fetches:**
- Trial summaries, eligibility criteria, and biomarker requirements
- Phase, status, and intervention details
- Geographic location data

### 7.4 Scheduled Ingest

The agent includes an APScheduler-based scheduler (`src/scheduler.py`) that
can periodically run ingest tasks.

**Configuration:**
- `ONCO_SCHEDULER_INTERVAL` (default: `168h` = 1 week)

To enable scheduled ingest in the API server, the scheduler is initialized
during the FastAPI lifespan and runs ingest jobs at the configured interval.

### 7.5 Custom Ingest Scripts

The `src/ingest/` directory contains parser modules that can be used to build
custom ingest pipelines:

| Parser Module | Purpose |
|--------------|---------|
| `base.py` | Base ingest class with common embedding and insertion logic |
| `civic_parser.py` | CIViC variant/evidence parsing |
| `clinical_trials_parser.py` | ClinicalTrials.gov API v2 parsing |
| `literature_parser.py` | PubMed abstract chunking and parsing |
| `guideline_parser.py` | NCCN/ASCO/ESMO guideline parsing |
| `oncokb_parser.py` | OncoKB variant annotation parsing |
| `outcome_parser.py` | Treatment outcome record parsing |
| `pathway_parser.py` | Signaling pathway parsing |
| `resistance_parser.py` | Resistance mechanism parsing |

---

## 8. Networking and Ports

### 8.1 Port Map

| Port | Protocol | Service | Container | Description |
|------|----------|---------|-----------|-------------|
| **8526** | HTTP | `onco-streamlit` | `onco-streamlit` | Streamlit clinical UI |
| **8527** | HTTP | `onco-api` | `onco-api` | FastAPI REST API + Swagger docs |
| **19530** | gRPC | `milvus-standalone` | `onco-milvus-standalone` | Milvus vector database client port |
| **9091** | HTTP | `milvus-standalone` | `onco-milvus-standalone` | Milvus metrics and health endpoint |
| **9000** | HTTP | `milvus-minio` | `onco-milvus-minio` | MinIO S3-compatible API |
| **9001** | HTTP | `milvus-minio` | `onco-milvus-minio` | MinIO web console |
| **2379** | HTTP | `milvus-etcd` | `onco-milvus-etcd` | etcd client endpoint (internal) |

### 8.2 Docker Network

All services communicate over the `onco-network` Docker bridge network.
Container-to-container communication uses service names as hostnames:

```
onco-api --> milvus-standalone:19530  (Milvus gRPC)
onco-streamlit --> milvus-standalone:19530  (Milvus gRPC)
milvus-standalone --> milvus-etcd:2379  (metadata)
milvus-standalone --> milvus-minio:9000  (object storage)
```

### 8.3 Exposing Ports

By default, only the following ports are mapped to the host:

- `8526` (Streamlit UI)
- `8527` (FastAPI API)
- `19530` (Milvus gRPC)
- `9091` (Milvus metrics)

MinIO ports (9000, 9001) and etcd (2379) are **not** exposed to the host by
default. To expose them, add port mappings in a `docker-compose.override.yml`:

```yaml
version: "3.8"
services:
  milvus-minio:
    ports:
      - "9000:9000"
      - "9001:9001"
  milvus-etcd:
    ports:
      - "2379:2379"
```

### 8.4 Changing Ports

Set `ONCO_API_PORT` and `ONCO_STREAMLIT_PORT` environment variables. You must
also update port mappings in `docker-compose.yml` or use an override file.

### 8.5 Firewall Rules (Production)

Allow ports 8526 and 8527; block direct Milvus access (19530, 9091) from
external networks. Example: `sudo ufw allow 8526/tcp && sudo ufw allow 8527/tcp`.

---

## 9. Storage and Persistence

### 9.1 Docker Volumes

The deployment uses three named Docker volumes for data persistence:

| Volume | Mount Point | Service | Content |
|--------|------------|---------|---------|
| `etcd_data` | `/etcd` | `milvus-etcd` | Milvus metadata (collection schemas, partitions) |
| `minio_data` | `/minio_data` | `milvus-minio` | Milvus index files and log segments |
| `milvus_data` | `/var/lib/milvus` | `milvus-standalone` | Milvus WAL, insert logs, and query cache |

### 9.2 Volume Lifecycle

Volumes persist across container restarts and `docker compose down`:

```bash
# Stop services (volumes preserved)
docker compose down

# Stop services AND delete volumes (DATA LOSS)
docker compose down -v

# List volumes
docker volume ls | grep -E "etcd_data|minio_data|milvus_data"

# Inspect a volume
docker volume inspect precision_oncology_agent_milvus_data
```

### 9.3 Disk Space Estimates

| Component | Base Size | With Full Ingest |
|-----------|----------|-----------------|
| Seed data (JSON) | ~768 KB | N/A |
| etcd metadata | ~50 MB | ~200 MB |
| MinIO (Milvus indexes) | ~500 MB | ~5-20 GB |
| Milvus data | ~500 MB | ~5-20 GB |
| Docker images | ~2 GB | ~2 GB |
| **Total** | **~3 GB** | **~12-42 GB** |

### 9.4 Cache Directory

The application creates a cache directory at `data/cache/` inside the
container (path: `/app/data/cache/`). This directory stores:

- Downloaded embedding model files (first run only, ~130 MB)
- Temporary processing artifacts

In the Dockerfile, this directory is pre-created and owned by the `oncouser`
non-root user.

### 9.5 Bind Mounts for Development

For development, you may want to bind-mount source code for live reloading:

```yaml
# docker-compose.override.yml
version: "3.8"
services:
  onco-api:
    volumes:
      - ./api:/app/api
      - ./src:/app/src
      - ./config:/app/config
    command:
      - uvicorn
      - api.main:app
      - --host=0.0.0.0
      - --port=8527
      - --reload
```

---

## 10. Monitoring and Metrics

### 10.1 Prometheus Metrics

The agent exposes Prometheus-compatible metrics at `GET /metrics` on the
FastAPI server (port 8527).

**Available metrics:**

| Metric | Type | Description |
|--------|------|-------------|
| `onco_agent_up` | Gauge | Service availability (1 = up, 0 = down) |
| `onco_collection_vectors` | Gauge | Vector count per collection (labeled by `collection`) |
| `onco_query_duration_seconds` | Histogram | RAG query latency |
| `onco_search_duration_seconds` | Histogram | Vector search latency |
| `onco_embedding_duration_seconds` | Histogram | Embedding generation latency |
| `onco_llm_tokens_total` | Counter | Total LLM tokens consumed |
| `onco_milvus_operations_total` | Counter | Milvus operation count (by operation type) |

**Scrape configuration for `prometheus.yml`:**

```yaml
scrape_configs:
  - job_name: "onco-agent"
    scrape_interval: 30s
    metrics_path: /metrics
    static_configs:
      - targets: ["onco-api:8527"]
        labels:
          service: "oncology-intelligence-agent"
```

### 10.2 Milvus Metrics

Milvus exposes its own metrics at `http://localhost:9091/metrics` in
Prometheus exposition format.

**Key Milvus metrics:**

| Metric | Description |
|--------|-------------|
| `milvus_datanode_flush_segments_total` | Segment flush count |
| `milvus_proxy_search_vectors_count` | Search request count |
| `milvus_querynode_search_latency` | Query latency histogram |

**Add to `prometheus.yml`:**

```yaml
  - job_name: "milvus"
    scrape_interval: 30s
    metrics_path: /metrics
    static_configs:
      - targets: ["milvus-standalone:9091"]
```

### 10.3 Grafana Dashboard

Add Prometheus as a data source (`http://prometheus:9090`), then create
panels using these PromQL queries:

| Panel | PromQL |
|-------|--------|
| Avg query latency (5m) | `rate(onco_query_duration_seconds_sum[5m]) / rate(onco_query_duration_seconds_count[5m])` |
| Total vectors | `sum(onco_collection_vectors)` |
| Request rate/min | `rate(onco_milvus_operations_total[1m]) * 60` |
| Service up | `onco_agent_up` |

### 10.4 Disabling Metrics

Set `ONCO_METRICS_ENABLED=false` to disable Prometheus metric collection.
The `/metrics` endpoint will still respond but will return minimal data.
When `prometheus_client` is not installed, the metrics module automatically
falls back to no-op stubs with zero overhead.

### 10.5 Health Dashboard

The `/health` endpoint returns a JSON summary with `status` (`"healthy"` or
`"degraded"`), per-collection vector counts, `total_vectors`, `version`, and
a `services` map showing boolean availability for `milvus`, `embedder`,
`rag_engine`, `intelligence_agent`, `case_manager`, `trial_matcher`, and
`therapy_ranker`. Use this for uptime monitoring and alerting.

---

## 11. Security Hardening

### 11.1 CORS Configuration

The API server enforces CORS (Cross-Origin Resource Sharing) restrictions.
By default, only requests from the Landing Page (`:8080`), Streamlit UI
(`:8526`), and the API itself (`:8527`) are allowed.

**Configure allowed origins:**

```bash
ONCO_CORS_ORIGINS="https://mtb.hospital.org,https://admin.hospital.org"
```

**In production, always restrict CORS to known origins.** The default
`localhost` values are suitable only for development.

### 11.2 API Key Protection

**Critical:** Never commit `ANTHROPIC_API_KEY` or `NCBI_API_KEY` to version
control.

Best practices:
- Use `.env` files (already in `.gitignore`)
- Use Docker secrets or a secrets manager (HashiCorp Vault, AWS Secrets Manager)
- Set environment variables directly on the host or in CI/CD pipelines
- Rotate API keys regularly

**Docker secrets example:**

```yaml
# docker-compose.yml addition
services:
  onco-api:
    secrets:
      - anthropic_key
    environment:
      ANTHROPIC_API_KEY_FILE: /run/secrets/anthropic_key

secrets:
  anthropic_key:
    file: ./secrets/anthropic_api_key.txt
```

### 11.3 Request Size Limits

The API enforces a maximum request body size via `ONCO_MAX_REQUEST_SIZE_MB`
(default: 10 MB). Requests exceeding this limit receive a `413 Payload Too
Large` response.

### 11.4 Non-Root Container User

The Dockerfile creates a dedicated `oncouser` with no shell access:

```dockerfile
RUN useradd -r -s /bin/false oncouser
USER oncouser
```

All application processes run as this non-root user inside the container.

### 11.5 TLS/HTTPS Configuration

The agent does not natively terminate TLS. For production, place a reverse
proxy (nginx, Caddy, or Traefik) in front of the services:

- Proxy `/` to `http://127.0.0.1:8526` (Streamlit -- requires WebSocket
  upgrade headers: `Upgrade`, `Connection`)
- Proxy `/api/` to `http://127.0.0.1:8527/` (FastAPI)
- Use TLS 1.2+ with valid certificates
- Set `X-Real-IP`, `X-Forwarded-For`, and `X-Forwarded-Proto` headers

### 11.6 MinIO Credentials

The default MinIO credentials (`minioadmin` / `minioadmin`) should be changed
in production:

```yaml
milvus-minio:
  environment:
    MINIO_ACCESS_KEY: ${MINIO_ACCESS_KEY}
    MINIO_SECRET_KEY: ${MINIO_SECRET_KEY}
```

### 11.7 Network Isolation

The default configuration exposes only ports 8526, 8527, 19530, and 9091 to
the host. For additional isolation, split into public and internal networks
(set `internal: true` on the Milvus infrastructure network) so that etcd
and MinIO are not reachable from the host.

---

## 12. Health Checks and Troubleshooting

### 12.1 Health Check Endpoints

| Service | Endpoint | Expected Response |
|---------|----------|-------------------|
| Milvus | `GET http://localhost:9091/healthz` | `{"status":"OK"}` |
| FastAPI | `GET http://localhost:8527/health` | JSON with `"status": "healthy"` |
| Streamlit | `GET http://localhost:8526/_stcore/health` | `ok` |
| MinIO | `GET http://localhost:9000/minio/health/live` | HTTP 200 |
| etcd | `etcdctl endpoint health` (in container) | `is healthy: true` |

### 12.2 Docker Health Check Configuration

Health checks are defined in both the Dockerfile and docker-compose.yml:

| Service | Check Command | Interval | Timeout | Retries | Start Period |
|---------|--------------|----------|---------|---------|--------------|
| Milvus | `curl -f http://localhost:9091/healthz` | 30s | 10s | 10 | 60s |
| FastAPI | `curl -f http://localhost:8527/health` | 30s | 10s | 3 | 30s |
| Streamlit | `curl -f http://localhost:8526/_stcore/health` | 30s | 10s | 3 | 40s |
| etcd | `etcdctl endpoint health` | 30s | 20s | 5 | -- |
| MinIO | `curl -f http://localhost:9000/minio/health/live` | 30s | 20s | 5 | -- |

### 12.3 Common Issues and Solutions

#### Issue: Milvus fails to start

**Symptoms:** `milvus-standalone` container restarts repeatedly.

**Diagnosis:**
```bash
docker compose logs milvus-standalone | tail -50
docker compose logs milvus-etcd | tail -20
docker compose logs milvus-minio | tail -20
```

**Common causes:**
1. **Insufficient memory:** Milvus needs at least 4 GB. Check with
   `docker stats`.
2. **etcd not ready:** etcd may take 30-60 seconds. The `depends_on` with
   health check should handle this, but check etcd logs.
3. **Port conflict:** Another service is using port 19530 or 9091.
   ```bash
   sudo lsof -i :19530
   sudo lsof -i :9091
   ```
4. **Corrupted data:** If volumes contain corrupted data from a previous run:
   ```bash
   docker compose down -v  # WARNING: deletes all data
   docker compose up -d
   ```

#### Issue: onco-setup exits with non-zero code

**Symptoms:** Seed scripts fail during collection creation or data loading.

**Diagnosis:**
```bash
docker compose logs onco-setup
```

**Common causes:**
1. **Milvus not ready:** The setup container started before Milvus was fully
   initialized. Re-run:
   ```bash
   docker compose rm -f onco-setup
   docker compose up onco-setup
   ```
2. **Collections already exist:** Use `--drop-existing` flag (already set in
   the default Docker Compose command).
3. **Missing seed data files:** Verify `data/reference/` contains all 10 JSON
   files.

#### Issue: API returns 503 "Service not initialised"

**Symptoms:** All API endpoints return 503 errors.

**Diagnosis:**
```bash
docker compose logs onco-api | tail -30
curl -s http://localhost:8527/health
```

**Common causes:**
1. **Milvus not reachable:** Verify `ONCO_MILVUS_HOST` and `ONCO_MILVUS_PORT`
   are correctly set.
2. **Embedding model download failed:** The first startup downloads
   `BAAI/bge-small-en-v1.5` (~130 MB). Ensure outbound HTTPS access to
   `huggingface.co`.
3. **Startup timeout:** The embedding model download can take several minutes
   on slow connections. Check the API logs for download progress.

#### Issue: LLM queries fail but search works

**Symptoms:** `/search` returns results but `/query` returns errors.

**Cause:** `ANTHROPIC_API_KEY` is not set or is invalid.

**Fix:**
```bash
# Check if the key is set in the container
docker compose exec onco-api env | grep ANTHROPIC

# Update the key
echo "ANTHROPIC_API_KEY=sk-ant-your-new-key" >> .env
docker compose up -d onco-api
```

#### Issue: Streamlit UI shows connection error

**Symptoms:** The Streamlit UI displays "Unable to connect to Milvus."

**Diagnosis:**
```bash
docker compose logs onco-streamlit | tail -20
```

**Fix:** Ensure `ONCO_MILVUS_HOST=milvus-standalone` (not `localhost`) in
the Docker environment, since the Streamlit container connects to Milvus
over the Docker network.

#### Issue: Slow embedding generation

**Symptoms:** Queries take 5+ seconds; seed scripts run slowly.

**Diagnosis:** Check if GPU acceleration is available:
```bash
python -c "import torch; print(torch.cuda.is_available())"
```

**Fix:** If GPU is available, ensure the Docker image has CUDA support.
For CPU-only deployments, increase `ONCO_EMBEDDING_BATCH_SIZE` to improve
throughput at the cost of memory.

### 12.4 Diagnostic Commands

| Command | Purpose |
|---------|---------|
| `docker compose ps` | Service status |
| `docker stats --no-stream` | Container resource usage |
| `curl -s localhost:8527/health \| python3 -m json.tool` | Full health check |
| `curl -s localhost:8527/collections \| python3 -m json.tool` | Collection stats |
| `curl -s localhost:8527/knowledge/stats \| python3 -m json.tool` | Knowledge stats |
| `curl -s localhost:9091/healthz` | Milvus health |
| `docker compose logs --tail=100 <service>` | View logs |
| `docker compose exec onco-api bash` | Shell into container |

---

## 13. Backup and Recovery

### 13.1 Volume Backup

Back up the three named Docker volumes to preserve all Milvus data:

```bash
#!/bin/bash
# backup-onco-volumes.sh
BACKUP_DIR="/backups/onco-agent/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$BACKUP_DIR"

# Stop services to ensure consistency
docker compose stop

# Backup each volume
for vol in etcd_data minio_data milvus_data; do
    echo "Backing up $vol..."
    docker run --rm \
        -v "precision_oncology_agent_${vol}:/data:ro" \
        -v "$BACKUP_DIR:/backup" \
        alpine tar czf "/backup/${vol}.tar.gz" -C /data .
done

# Restart services
docker compose start

echo "Backup complete: $BACKUP_DIR"
ls -lh "$BACKUP_DIR"
```

### 13.2 Volume Restore

```bash
#!/bin/bash
# restore-onco-volumes.sh <backup-dir>
BACKUP_DIR="$1"
[ -z "$BACKUP_DIR" ] && echo "Usage: $0 /path/to/backup" && exit 1

docker compose down
for vol in etcd_data minio_data milvus_data; do
    docker volume rm "precision_oncology_agent_${vol}" 2>/dev/null
    docker volume create "precision_oncology_agent_${vol}"
    docker run --rm \
        -v "precision_oncology_agent_${vol}:/data" \
        -v "$BACKUP_DIR:/backup:ro" \
        alpine tar xzf "/backup/${vol}.tar.gz" -C /data
done
docker compose up -d
```

### 13.3 Re-Seed from Scratch

If backups are unavailable, you can reconstruct the entire knowledge base
from the seed data and live ingest scripts:

```bash
# 1. Remove all data
docker compose down -v

# 2. Start infrastructure
docker compose up -d

# 3. Wait for setup to complete
docker compose logs -f onco-setup

# 4. Run live ingest (optional, for full data)
docker compose exec onco-api python scripts/ingest_civic.py
docker compose exec onco-api python scripts/ingest_pubmed.py
docker compose exec onco-api python scripts/ingest_clinical_trials.py
```

### 13.4 Backup Schedule

Schedule daily backups via cron (e.g., `0 2 * * *`) and weekly full backups
on Sundays. Log to `/var/log/onco-backup.log`.

### 13.5 Backup Verification

After restoring, verify with `curl -s localhost:8527/collections`, run a test
search query, and execute `python scripts/validate_e2e.py`.

---

## 14. Scaling Considerations

### 14.1 Vertical Scaling

The simplest scaling approach is to increase resources allocated to existing
containers.

**API server workers:**

```yaml
# docker-compose.override.yml
services:
  onco-api:
    command:
      - uvicorn
      - api.main:app
      - --host=0.0.0.0
      - --port=8527
      - --workers=8    # Increase from default 2
```

Each uvicorn worker loads its own copy of the embedding model (~130 MB).
Plan memory accordingly: `workers * 200 MB + base overhead`.

**Milvus resource allocation:**

```yaml
services:
  milvus-standalone:
    deploy:
      resources:
        limits:
          memory: 32G
        reservations:
          memory: 16G
```

### 14.2 Horizontal Scaling: API Layer

The FastAPI server is stateless (all state lives in Milvus). Run multiple
replicas behind a load balancer:

```yaml
services:
  onco-api:
    deploy:
      replicas: 3
```

Use nginx, HAProxy, or Traefik to distribute requests across the replicas.
Remove host port mapping from individual containers when using a load balancer.

### 14.3 Milvus Scaling

For larger deployments, consider migrating from Milvus standalone to
Milvus cluster mode:

| Deployment | Collections | Vectors | Concurrent Queries |
|-----------|------------|---------|-------------------|
| Standalone | 11 | Up to ~10M | 10-50 |
| Cluster | 11+ | 100M+ | 100+ |

Milvus cluster mode requires:
- Separate etcd cluster (3+ nodes)
- Separate MinIO or S3 storage
- Query nodes, data nodes, and index nodes

Refer to the [Milvus documentation](https://milvus.io/docs/install_cluster-milvusoperator.md)
for cluster deployment.

### 14.4 Embedding Model Scaling

- **GPU acceleration:** ~1000 embeddings/sec (RTX 4090) vs. ~100/sec (CPU).
- **Dedicated embedding service:** Run `sentence-transformers` as a standalone
  microservice with batching.
- **Model server:** Triton Inference Server for production embedding serving.

### 14.5 Performance Benchmarks

Typical latencies on DGX Spark hardware:

| Operation | Latency (p50) | Latency (p95) |
|-----------|--------------|--------------|
| Vector search (single collection) | ~15 ms | ~50 ms |
| Multi-collection RAG search | ~100 ms | ~300 ms |
| Full RAG query (with LLM) | ~2 s | ~5 s |
| Embedding generation (single text) | ~10 ms (GPU) / ~50 ms (CPU) | ~30 ms (GPU) / ~150 ms (CPU) |
| Trial matching (top 10) | ~200 ms | ~500 ms |

---

## 15. Integration with HCLS AI Factory

### 15.1 Pipeline Position

The Oncology Intelligence Agent operates as part of the HCLS AI Factory
three-stage pipeline:

```
Stage 1: Genomics Pipeline (Parabricks/DeepVariant)
    FASTQ --> VCF
        |
        v
Stage 2: RAG/Chat Pipeline (Milvus + Claude)
    VCF --> Variant Interpretation
        |
        v
    Oncology Intelligence Agent
    (MTB Decision Support)
        |
        v
Stage 3: Drug Discovery Pipeline (BioNeMo/DiffDock)
    Candidate Molecules --> Docking Scores
```

### 15.2 Shared Collections

The `genomic_evidence` collection is read-only and shared with the Stage 1
genomics pipeline. When the genomics pipeline processes a VCF file, it
writes annotated variants to this collection. The oncology agent reads from
it to provide cross-modal genomic context.

**Configuration:** The `ONCO_COLLECTION_GENOMIC` setting (default:
`genomic_evidence`) must match the collection name used by the genomics
pipeline.

### 15.3 Landing Page Integration

The HCLS AI Factory Landing Page (port 8080) provides a unified dashboard
with health monitoring for all agents. The oncology agent registers at:

- Health endpoint: `http://onco-api:8527/health`
- UI link: `http://localhost:8526`

Ensure `ONCO_CORS_ORIGINS` includes the landing page origin
(`http://localhost:8080`).

### 15.4 Docker Compose Integration

When running as part of the full HCLS AI Factory stack
(`docker-compose.dgx-spark.yml`), the oncology agent services are defined
within the main compose file. Ensure network connectivity between the
oncology services and the shared Milvus instance if using a centralized
vector database.

### 15.5 Cross-Agent Communication

The oncology agent can receive cross-modal triggers from other HCLS AI
Factory agents:

| Source Agent | Data Flow | Collection |
|-------------|-----------|------------|
| Biomarker Agent | Biomarker panel results | `onco_biomarkers` |
| CAR-T Agent | Immunotherapy candidates | `onco_therapies` |
| Imaging Agent | Radiomics features | Cross-modal trigger |
| Autoimmune Agent | Immune checkpoint data | `onco_resistance` |

### 15.6 Event Bus

The `api/routes/events.py` module provides an event endpoint for
inter-agent communication. Other agents can publish events that the
oncology agent processes for knowledge base updates.

---

## 16. Updating and Maintenance

### 16.1 Updating the Agent Code

```bash
# 1. Pull latest code
git pull origin main

# 2. Rebuild Docker images
docker compose build --no-cache

# 3. Restart services (preserves data volumes)
docker compose up -d

# 4. Re-run setup if schema changes occurred
docker compose rm -f onco-setup
docker compose up onco-setup
```

### 16.2 Updating Dependencies

```bash
# 1. Update requirements.txt with new versions
# 2. Rebuild the image
docker compose build --no-cache onco-api onco-streamlit onco-setup

# 3. Restart
docker compose up -d
```

### 16.3 Updating Seed Data

To refresh the seed data without losing live-ingested data:

```bash
# Run individual seed scripts (appends to existing collections)
docker compose exec onco-api python scripts/seed_variants.py

# Or drop and re-seed everything
docker compose exec onco-api python scripts/setup_collections.py --drop-existing --seed
```

### 16.4 Updating Milvus

To upgrade the Milvus version:

1. **Back up all volumes** (see [Section 13](#13-backup-and-recovery)).
2. Update the image tag in `docker-compose.yml`:
   ```yaml
   milvus-standalone:
     image: milvusdb/milvus:v2.5-latest  # Updated from v2.4
   ```
3. Review the [Milvus release notes](https://milvus.io/docs/release_notes.md)
   for breaking changes.
4. Restart:
   ```bash
   docker compose down
   docker compose up -d
   ```

### 16.5 Log Rotation

Configure Docker log rotation in `/etc/docker/daemon.json` with
`"max-size": "50m"` and `"max-file": "5"`. Restart Docker after changing.

### 16.6 Periodic Maintenance Tasks

| Task | Frequency | Command |
|------|-----------|---------|
| Backup volumes | Daily | `backup-onco-volumes.sh` |
| Refresh PubMed data | Weekly | `python scripts/ingest_pubmed.py` |
| Refresh ClinicalTrials.gov | Weekly | `python scripts/ingest_clinical_trials.py` |
| Refresh CIViC variants | Weekly | `python scripts/ingest_civic.py` |
| Compact Milvus | Monthly | Automatic (configured via etcd compaction) |
| Rotate API keys | Quarterly | Update `.env` and restart |
| Update Docker images | As needed | `docker compose build --no-cache` |
| Review metrics/alerts | Weekly | Grafana dashboard review |
| Test backup restore | Monthly | Restore to staging environment |

### 16.7 Version Pinning

For production stability, pin all image versions:

```yaml
services:
  milvus-etcd:
    image: quay.io/coreos/etcd:v3.5.5         # Pinned
  milvus-minio:
    image: minio/minio:RELEASE.2023-03-20T20-16-18Z  # Pinned
  milvus-standalone:
    image: milvusdb/milvus:v2.4.0              # Pinned to specific release
```

---

## Appendix A: Complete docker-compose.yml

The canonical `docker-compose.yml` is maintained at the repository root:

```
ai_agent_adds/precision_oncology_agent/agent/docker-compose.yml
```

It defines all 6 services (`milvus-etcd`, `milvus-minio`, `milvus-standalone`,
`onco-streamlit`, `onco-api`, `onco-setup`), 3 volumes (`etcd_data`,
`minio_data`, `milvus_data`), and the `onco-network` bridge network.

Key configuration facts for quick reference:

| Service | Image | Exposed Ports | Restart Policy |
|---------|-------|---------------|----------------|
| `milvus-etcd` | `quay.io/coreos/etcd:v3.5.5` | _(none)_ | `unless-stopped` |
| `milvus-minio` | `minio/minio:RELEASE.2023-03-20T20-16-18Z` | _(none)_ | `unless-stopped` |
| `milvus-standalone` | `milvusdb/milvus:v2.4-latest` | `19530`, `9091` | `unless-stopped` |
| `onco-streamlit` | Local `Dockerfile` | `8526` | `unless-stopped` |
| `onco-api` | Local `Dockerfile` | `8527` | `unless-stopped` |
| `onco-setup` | Local `Dockerfile` | _(none)_ | `no` |

**Startup order:** etcd + MinIO (parallel) --> Milvus (waits for both) -->
Streamlit + API + Setup (wait for Milvus healthy).

**etcd configuration:**
- Backend quota: 4 GB (`ETCD_QUOTA_BACKEND_BYTES=4294967296`)
- Auto-compaction: revision mode, retention 1000
- Snapshot count: 50000

**MinIO configuration:**
- Default credentials: `minioadmin` / `minioadmin`
- Console address: `:9001`

**Milvus configuration:**
- Security opt: `seccomp:unconfined`
- Health check: `start_period: 60s`, `retries: 10`

**onco-api:**
- Command: `uvicorn api.main:app --host=0.0.0.0 --port=8527 --workers=2`
- Health check: `start-period: 30s`, `retries: 3`

**onco-setup:**
- Runs `setup_collections.py --drop-existing` then 9 seed scripts sequentially
- Exits after completion (restart: `"no"`)

---

## Appendix B: Environment Variable Quick Reference

Complete list of all environment variables accepted by the Oncology
Intelligence Agent, sorted alphabetically within categories.

### Required

| Variable | Example | Description |
|----------|---------|-------------|
| `ANTHROPIC_API_KEY` | `sk-ant-api03-...` | Anthropic API key for Claude LLM |

### Connection

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_API_BASE_URL` | `http://localhost:8527` | API base URL for Streamlit UI |
| `ONCO_API_HOST` | `0.0.0.0` | FastAPI bind address |
| `ONCO_API_PORT` | `8527` | FastAPI listen port |
| `ONCO_MILVUS_HOST` | `localhost` | Milvus hostname |
| `ONCO_MILVUS_PORT` | `19530` | Milvus gRPC port |
| `ONCO_STREAMLIT_PORT` | `8526` | Streamlit server port |

### Embedding & LLM

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_EMBEDDING_BATCH_SIZE` | `32` | Embedding batch size |
| `ONCO_EMBEDDING_DIM` | `384` | Embedding vector dimension |
| `ONCO_EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | HuggingFace embedding model |
| `ONCO_LLM_MODEL` | `claude-sonnet-4-20250514` | Claude model name |
| `ONCO_LLM_PROVIDER` | `anthropic` | LLM provider |

### RAG Search

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_MIN_COLLECTIONS_FOR_SUFFICIENT` | `2` | Min collections with hits |
| `ONCO_MIN_SIMILARITY_SCORE` | `0.30` | Absolute minimum similarity |
| `ONCO_MIN_SUFFICIENT_HITS` | `3` | Min hits for sufficient evidence |
| `ONCO_SCORE_THRESHOLD` | `0.4` | Similarity score threshold |
| `ONCO_TOP_K` | `5` | Results per collection |

### Collection Weights

| Variable | Default | Variable | Default |
|----------|---------|----------|---------|
| `ONCO_WEIGHT_VARIANTS` | `0.18` | `ONCO_WEIGHT_BIOMARKERS` | `0.08` |
| `ONCO_WEIGHT_LITERATURE` | `0.16` | `ONCO_WEIGHT_RESISTANCE` | `0.07` |
| `ONCO_WEIGHT_THERAPIES` | `0.14` | `ONCO_WEIGHT_PATHWAYS` | `0.06` |
| `ONCO_WEIGHT_GUIDELINES` | `0.12` | `ONCO_WEIGHT_OUTCOMES` | `0.04` |
| `ONCO_WEIGHT_TRIALS` | `0.10` | `ONCO_WEIGHT_GENOMIC` | `0.03` |
| | | `ONCO_WEIGHT_CASES` | `0.02` |

### Trial Matching

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_TRIAL_WEIGHT_BIOMARKER` | `0.40` | Biomarker match weight |
| `ONCO_TRIAL_WEIGHT_PHASE` | `0.20` | Trial phase weight |
| `ONCO_TRIAL_WEIGHT_SEMANTIC` | `0.25` | Semantic similarity weight |
| `ONCO_TRIAL_WEIGHT_STATUS` | `0.15` | Recruitment status weight |

### Citation & Cross-Modal

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_CITATION_MODERATE_THRESHOLD` | `0.60` | Moderate evidence threshold |
| `ONCO_CITATION_STRONG_THRESHOLD` | `0.75` | Strong evidence threshold |
| `ONCO_CROSS_MODAL_ENABLED` | `true` | Enable cross-modal analysis |
| `ONCO_CROSS_MODAL_THRESHOLD` | `0.40` | Cross-modal trigger threshold |
| `ONCO_GENOMIC_TOP_K` | `5` | Genomic evidence result count |
| `ONCO_IMAGING_TOP_K` | `5` | Imaging result count |

### External APIs

| Variable | Default | Description |
|----------|---------|-------------|
| `NCBI_API_KEY` | _(none)_ | PubMed API key (optional) |
| `ONCO_CIVIC_BASE_URL` | `https://civicdb.org/api` | CIViC API base URL |
| `ONCO_CT_GOV_BASE_URL` | `https://clinicaltrials.gov/api/v2` | ClinicalTrials.gov API URL |
| `ONCO_PUBMED_MAX_RESULTS` | `5000` | Max PubMed fetch count |

### Operational

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_CONVERSATION_MEMORY_DEPTH` | `3` | Conversation turns to retain |
| `ONCO_CORS_ORIGINS` | `http://localhost:8080,...` | Allowed CORS origins |
| `ONCO_MAX_REQUEST_SIZE_MB` | `10` | Max request body size |
| `ONCO_METRICS_ENABLED` | `true` | Enable Prometheus metrics |
| `ONCO_SCHEDULER_INTERVAL` | `168h` | Data refresh interval |

### PDF Branding

| Variable | Default | Description |
|----------|---------|-------------|
| `ONCO_PDF_BRAND_COLOR_B` | `0` | Brand color blue (0-255) |
| `ONCO_PDF_BRAND_COLOR_G` | `185` | Brand color green (0-255) |
| `ONCO_PDF_BRAND_COLOR_R` | `118` | Brand color red (0-255) |

### Collection Names

| Variable | Default |
|----------|---------|
| `ONCO_COLLECTION_BIOMARKERS` | `onco_biomarkers` |
| `ONCO_COLLECTION_CASES` | `onco_cases` |
| `ONCO_COLLECTION_GENOMIC` | `genomic_evidence` |
| `ONCO_COLLECTION_GUIDELINES` | `onco_guidelines` |
| `ONCO_COLLECTION_LITERATURE` | `onco_literature` |
| `ONCO_COLLECTION_OUTCOMES` | `onco_outcomes` |
| `ONCO_COLLECTION_PATHWAYS` | `onco_pathways` |
| `ONCO_COLLECTION_RESISTANCE` | `onco_resistance` |
| `ONCO_COLLECTION_THERAPIES` | `onco_therapies` |
| `ONCO_COLLECTION_TRIALS` | `onco_trials` |
| `ONCO_COLLECTION_VARIANTS` | `onco_variants` |

---

*Document generated from verified codebase analysis of the Precision Oncology
Intelligence Agent. All defaults, file paths, and configuration values were
extracted from `config/settings.py`, `docker-compose.yml`, `Dockerfile`,
`requirements.txt`, and the `scripts/` directory.*
