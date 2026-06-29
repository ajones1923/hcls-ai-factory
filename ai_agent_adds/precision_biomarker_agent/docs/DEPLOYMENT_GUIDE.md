# Precision Biomarker Intelligence Agent -- Deployment Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**License:** Apache 2.0

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Quick Start (Docker)](#3-quick-start-docker)
4. [Local Development Setup](#4-local-development-setup)
5. [DGX Spark Production Deployment](#5-dgx-spark-production-deployment)
6. [Environment Variables Reference](#6-environment-variables-reference)
7. [Milvus Collection Setup](#7-milvus-collection-setup)
8. [Security Hardening](#8-security-hardening)
9. [Monitoring & Health Checks](#9-monitoring--health-checks)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Overview

The Precision Biomarker Intelligence Agent is a multi-collection RAG-powered
service that provides:

- **Biological age calculation** from clinical biomarker panels
- **Disease trajectory analysis** with genetic risk modifiers
- **Pharmacogenomic (PGx) drug-gene interaction mapping** (CPIC-aligned)
- **Genotype-adjusted reference ranges** for lab interpretation
- **Critical value alerting** with clinical decision support
- **Discordance detection** across related biomarker panels
- **Ashkenazi Jewish carrier screening** integration
- **Nutrition and lifestyle recommendations** linked to genotype

The agent searches across **13 specialized Milvus collections** plus 1 shared
read-only genomic collection, weighting results by clinical relevance before
synthesizing answers with Claude (Anthropic).

### Architecture

```
                     +-------------------+
                     |   Streamlit UI    |  :8528
                     |  (biomarker_ui)   |
                     +--------+----------+
                              |
                              v
                     +-------------------+
                     |   FastAPI REST    |  :8529
                     |   (api/main.py)   |
                     +--------+----------+
                              |
              +---------------+---------------+
              |               |               |
   +----------v--+   +-------v------+   +----v-----------+
   | RAG Engine  |   | Analysis     |   | Report         |
   | (rag_engine |   | Modules      |   | Generator      |
   |  .py)       |   | (bio_age,    |   | (export.py,    |
   |             |   |  trajectory, |   |  report_gen.py)|
   +------+------+   |  pgx, ...)   |   +----------------+
          |           +--------------+
          v
   +------+------+
   |  Milvus 2.4 |  :19530
   | (14 colls)  |
   +------+------+
          |
   +------+------+------+
   |      |             |
   v      v             v
  etcd   MinIO    Embedding Model
 (meta)  (index)  (BGE-small-en-v1.5)
```

### Service Summary

| Service              | Port  | Description                              |
|----------------------|-------|------------------------------------------|
| `biomarker-streamlit`| 8528  | Streamlit interactive UI                 |
| `biomarker-api`      | 8529  | FastAPI REST server (2 Uvicorn workers)  |
| `milvus-standalone`  | 19530 | Milvus 2.4 vector database (gRPC)        |
| `milvus-standalone`  | 9091  | Milvus metrics / health endpoint         |
| `milvus-etcd`        | 2379  | etcd metadata store (internal)           |
| `milvus-minio`       | 9000  | MinIO object storage (internal)          |
| `biomarker-setup`    | --    | One-shot collection creation + seeding   |

---

## 2. Prerequisites

### 2.1 Hardware Requirements

| Tier        | CPU     | RAM    | Disk   | GPU                | Use Case             |
|-------------|---------|--------|--------|--------------------|----------------------|
| Minimum     | 4 cores | 16 GB  | 50 GB  | None               | Dev / testing        |
| Recommended | 8 cores | 32 GB  | 100 GB | None               | Staging              |
| Production  | 16+ cores| 64 GB | 200 GB | NVIDIA DGX Spark   | Full deployment      |

**Disk breakdown (approximate):**

- Milvus data (14 collections, indexes): ~5 GB
- Embedding model cache (`BAAI/bge-small-en-v1.5`): ~130 MB
- Docker images (Python, Milvus, etcd, MinIO): ~8 GB
- Application code + reference data: ~500 MB
- Logs and temporary files: ~1 GB

### 2.2 Software Requirements

| Software         | Version     | Notes                                  |
|------------------|-------------|----------------------------------------|
| Docker           | 24.0+       | With Docker Compose V2 (`docker compose`) |
| Docker Compose   | 2.20+       | Included with Docker Desktop           |
| Python           | 3.10+       | For local development only             |
| Git              | 2.30+       | For cloning the repository             |
| NVIDIA Driver    | 535+        | Only if using GPU acceleration         |
| NVIDIA Container Toolkit | 1.14+ | Only for DGX Spark deployment       |

Verify Docker installation:

```bash
docker --version
docker compose version
```

### 2.3 API Keys

| Key                    | Required | How to Obtain                          |
|------------------------|----------|----------------------------------------|
| Anthropic API Key      | Yes      | https://console.anthropic.com          |

The Anthropic API key is required for LLM-powered analysis. Without it, the
RAG engine can still retrieve and rank evidence from Milvus, but synthesis
and natural-language responses will be unavailable.

### 2.4 Network Requirements

| Direction  | Port  | Protocol | Purpose                        |
|------------|-------|----------|--------------------------------|
| Inbound    | 8528  | HTTP     | Streamlit UI                   |
| Inbound    | 8529  | HTTP     | FastAPI REST API               |
| Outbound   | 443   | HTTPS    | Anthropic API (api.anthropic.com) |
| Internal   | 19530 | gRPC     | Milvus client connections      |
| Internal   | 9091  | HTTP     | Milvus health / metrics        |
| Internal   | 2379  | HTTP     | etcd metadata                  |
| Internal   | 9000  | HTTP     | MinIO object storage           |

---

## 3. Quick Start (Docker)

The fastest path from zero to running. All commands are run from the
`precision_biomarker_agent/` directory.

### 3.1 Clone and Navigate

```bash
git clone https://github.com/your-org/hcls-ai-factory.git
cd hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
```

### 3.2 Configure Environment

```bash
cp .env.example .env
```

Edit `.env` and set your Anthropic API key:

```bash
# .env
BIOMARKER_ANTHROPIC_API_KEY=sk-ant-api03-XXXXXXXXXXXXXXXXXXXXXXXX
```

### 3.3 Build and Launch

```bash
# Build images and start all services
docker compose up -d --build
```

This starts six services:

1. **milvus-etcd** -- metadata store (waits for healthy)
2. **milvus-minio** -- object storage (waits for healthy)
3. **milvus-standalone** -- vector database (waits for etcd + MinIO)
4. **biomarker-setup** -- creates 13 collections, seeds reference data, then exits
5. **biomarker-streamlit** -- Streamlit UI on port 8528
6. **biomarker-api** -- FastAPI REST server on port 8529

### 3.4 Watch Setup Progress

```bash
# Follow the one-shot setup container logs
docker compose logs -f biomarker-setup
```

You should see output like:

```
===== Biomarker Setup: Creating collections =====
[INFO] Creating collection: biomarker_reference (384-dim)
[INFO] Creating collection: biomarker_genetic_variants (384-dim)
...
[INFO] All 13 collections created successfully
===== Seeding all reference data =====
[INFO] Seeding biomarker_reference: 247 records
[INFO] Seeding biomarker_genetic_variants: 183 records
...
===== Biomarker Setup complete! =====
```

### 3.5 Verify Deployment

```bash
# Health check -- FastAPI
curl -s http://localhost:8529/health | python3 -m json.tool

# Expected output:
# {
#     "status": "healthy",
#     "collections": 14,
#     "total_vectors": 2847,
#     "agent_ready": true
# }

# Health check -- Milvus
curl -s http://localhost:9091/healthz

# Open Streamlit UI
echo "Open http://localhost:8528 in your browser"
```

### 3.6 Stop Services

```bash
# Stop all services (preserves data volumes)
docker compose down

# Stop and remove all data (full reset)
docker compose down -v
```

---

## 4. Local Development Setup

For development without Docker. Useful for debugging, writing tests, and
iterating on analysis modules.

### 4.1 Create Virtual Environment

```bash
cd ai_agent_adds/precision_biomarker_agent

python3 -m venv .venv
source .venv/bin/activate
```

### 4.2 Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

Key packages installed:

| Package               | Purpose                              |
|-----------------------|--------------------------------------|
| `pymilvus>=2.4`       | Milvus vector database client        |
| `sentence-transformers`| BGE-small-en-v1.5 embedding model   |
| `anthropic`           | Claude LLM client                    |
| `fastapi` + `uvicorn` | REST API server                      |
| `streamlit`           | Interactive UI                       |
| `pydantic-settings`   | Typed configuration with env vars    |
| `reportlab`           | PDF report generation                |
| `prometheus-client`   | Metrics collection                   |
| `pytest`              | Test framework                       |

### 4.3 Start Milvus (Docker)

Even for local development, Milvus runs in Docker:

```bash
# Start only the Milvus stack (etcd + MinIO + Milvus)
docker compose up -d milvus-etcd milvus-minio milvus-standalone
```

Wait for Milvus to become healthy:

```bash
# Poll until healthy (usually 30-60 seconds)
until curl -sf http://localhost:9091/healthz > /dev/null 2>&1; do
    echo "Waiting for Milvus..."
    sleep 5
done
echo "Milvus is ready"
```

### 4.4 Configure Environment

```bash
cp .env.example .env
```

Edit `.env`:

```bash
BIOMARKER_ANTHROPIC_API_KEY=sk-ant-api03-XXXXXXXXXXXXXXXXXXXXXXXX
BIOMARKER_MILVUS_HOST=localhost
BIOMARKER_MILVUS_PORT=19530
```

Alternatively, export environment variables directly:

```bash
export BIOMARKER_ANTHROPIC_API_KEY="sk-ant-api03-XXXXXXXXXXXXXXXXXXXXXXXX"
export BIOMARKER_MILVUS_HOST="localhost"
export BIOMARKER_MILVUS_PORT="19530"
```

### 4.5 Create Collections and Seed Data

```bash
# Create all 13 collections with proper schemas
python3 scripts/setup_collections.py --drop-existing

# Seed reference data into all collections
python3 scripts/seed_all.py
```

### 4.6 Start the FastAPI Server

```bash
# Development mode with auto-reload
uvicorn api.main:app --host 0.0.0.0 --port 8529 --reload

# Production-like mode with 2 workers
uvicorn api.main:app --host 0.0.0.0 --port 8529 --workers 2
```

### 4.7 Start the Streamlit UI

In a separate terminal:

```bash
source .venv/bin/activate
streamlit run app/biomarker_ui.py \
    --server.port 8528 \
    --server.address 0.0.0.0 \
    --server.headless true \
    --browser.gatherUsageStats false
```

### 4.8 Run the Test Suite

```bash
# Run all 709 tests with verbose output
python3 -m pytest tests/ -v

# Run tests with coverage
python3 -m pytest tests/ --cov=src --cov=api --cov-report=term-missing

# Run specific test categories
python3 -m pytest tests/ -v -m "not slow"           # Skip slow tests
python3 -m pytest tests/ -v -m "not integration"     # Skip integration tests
python3 -m pytest tests/ -v -m "security"            # Security tests only

# Run a single test file
python3 -m pytest tests/test_rag_engine.py -v
python3 -m pytest tests/test_pharmacogenomics.py -v
python3 -m pytest tests/test_critical_values.py -v
```

Test file inventory:

| Test File                       | Module Under Test               |
|---------------------------------|---------------------------------|
| `test_agent.py`                 | `src/agent.py`                  |
| `test_api.py`                   | `api/main.py` + routes          |
| `test_biological_age.py`        | `src/biological_age.py`         |
| `test_collections.py`           | `src/collections.py`            |
| `test_critical_values.py`       | `src/critical_values.py`        |
| `test_discordance_detector.py`  | `src/discordance_detector.py`   |
| `test_disease_trajectory.py`    | `src/disease_trajectory.py`     |
| `test_edge_cases.py`            | Cross-module edge cases         |
| `test_export.py`                | `src/export.py`                 |
| `test_genotype_adjustment.py`   | `src/genotype_adjustment.py`    |
| `test_integration.py`           | End-to-end integration          |
| `test_lab_range_interpreter.py` | `src/lab_range_interpreter.py`  |
| `test_longitudinal.py`          | Longitudinal tracking           |
| `test_models.py`                | `src/models.py`                 |
| `test_pharmacogenomics.py`      | `src/pharmacogenomics.py`       |
| `test_rag_engine.py`            | `src/rag_engine.py`             |
| `test_report_generator.py`      | `src/report_generator.py`       |
| `test_ui.py`                    | `app/biomarker_ui.py`           |

Coverage threshold is configured at **85%** in `pyproject.toml`.

### 4.9 Code Quality

```bash
# Lint with ruff (configured in pyproject.toml)
pip install ruff
ruff check src/ api/ app/

# Format
ruff format src/ api/ app/
```

---

## 5. DGX Spark Production Deployment

Production deployment on NVIDIA DGX Spark with external port mapping and
systemd service management.

### 5.1 Port Mapping

On DGX Spark, external-facing ports are remapped to avoid conflicts with
other HCLS AI Factory services:

| Service          | Container Port | DGX External Port | URL                           |
|------------------|----------------|--------------------|-------------------------------|
| Streamlit UI     | 8528           | **8502**           | `http://<dgx-ip>:8502`       |
| FastAPI REST     | 8529           | **8102**           | `http://<dgx-ip>:8102`       |
| Milvus gRPC      | 19530          | 19530 (internal)   | Not exposed externally        |
| Milvus health    | 9091           | 9091 (internal)    | Not exposed externally        |

### 5.2 Production Docker Compose Override

Create `docker-compose.dgx.yml` for DGX Spark-specific overrides:

```yaml
# docker-compose.dgx.yml -- DGX Spark production overrides
version: "3.8"

services:
  milvus-standalone:
    ports:
      - "127.0.0.1:19530:19530"   # Bind to localhost only
      - "127.0.0.1:9091:9091"

  biomarker-streamlit:
    ports:
      - "8502:8528"               # External: 8502 -> Container: 8528
    environment:
      BIOMARKER_CORS_ORIGINS: "http://localhost:8080,http://${DGX_HOST:-localhost}:8502"
    restart: always
    deploy:
      resources:
        limits:
          memory: 4G
        reservations:
          memory: 2G

  biomarker-api:
    ports:
      - "8102:8529"               # External: 8102 -> Container: 8529
    environment:
      BIOMARKER_API_KEY: "${BIOMARKER_API_KEY}"
      BIOMARKER_CORS_ORIGINS: "http://localhost:8080,http://${DGX_HOST:-localhost}:8502"
    restart: always
    deploy:
      resources:
        limits:
          memory: 8G
        reservations:
          memory: 4G
```

Launch with both compose files:

```bash
docker compose \
    -f docker-compose.yml \
    -f docker-compose.dgx.yml \
    up -d --build
```

### 5.3 Systemd Service (Auto-Start on Boot)

Create `/etc/systemd/system/biomarker-agent.service`:

```ini
[Unit]
Description=Precision Biomarker Intelligence Agent
After=docker.service
Requires=docker.service

[Service]
Type=oneshot
RemainAfterExit=yes
WorkingDirectory=/home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
ExecStart=/usr/bin/docker compose -f docker-compose.yml -f docker-compose.dgx.yml up -d
ExecStop=/usr/bin/docker compose down
TimeoutStartSec=300
User=adam
Group=adam

[Install]
WantedBy=multi-user.target
```

Enable and start:

```bash
sudo systemctl daemon-reload
sudo systemctl enable biomarker-agent.service
sudo systemctl start biomarker-agent.service

# Check status
sudo systemctl status biomarker-agent.service
```

### 5.4 Log Rotation

Create `/etc/logrotate.d/biomarker-agent`:

```
/home/adam/projects/hcls-ai-factory/logs/biomarker-*.log {
    daily
    rotate 14
    compress
    delaycompress
    missingok
    notifempty
    copytruncate
    maxsize 100M
}
```

### 5.5 Production Environment Variables

For DGX Spark production, create `.env.production`:

```bash
# .env.production -- DGX Spark production settings
BIOMARKER_ANTHROPIC_API_KEY=sk-ant-api03-XXXXXXXXXXXXXXXXXXXXXXXX
BIOMARKER_MILVUS_HOST=milvus-standalone
BIOMARKER_MILVUS_PORT=19530
BIOMARKER_API_HOST=0.0.0.0
BIOMARKER_API_PORT=8529
BIOMARKER_LLM_MODEL=claude-sonnet-4-6
BIOMARKER_METRICS_ENABLED=True
BIOMARKER_API_KEY=your-strong-api-key-here
BIOMARKER_CORS_ORIGINS=http://localhost:8080
BIOMARKER_MAX_REQUEST_SIZE_MB=10
BIOMARKER_REQUEST_TIMEOUT_SECONDS=60
BIOMARKER_MILVUS_TIMEOUT_SECONDS=10
BIOMARKER_LLM_MAX_RETRIES=3
BIOMARKER_CITATION_HIGH_THRESHOLD=0.75
BIOMARKER_CITATION_MEDIUM_THRESHOLD=0.60
```

Use it:

```bash
cp .env.production .env
docker compose -f docker-compose.yml -f docker-compose.dgx.yml up -d --build
```

### 5.6 Validate Production Deployment

```bash
DGX_IP="<your-dgx-spark-ip>"

# FastAPI health
curl -s http://${DGX_IP}:8102/health | python3 -m json.tool

# FastAPI health (with API key if configured)
curl -s -H "X-API-Key: your-strong-api-key-here" \
    http://${DGX_IP}:8102/v1/health | python3 -m json.tool

# Streamlit UI -- open in browser
echo "http://${DGX_IP}:8502"

# Collection status
curl -s http://${DGX_IP}:8102/collections | python3 -m json.tool

# Prometheus metrics
curl -s http://${DGX_IP}:8102/metrics
```

---

## 6. Environment Variables Reference

All environment variables use the `BIOMARKER_` prefix. The underlying
configuration class is `PrecisionBiomarkerSettings` in `config/settings.py`,
which uses Pydantic `BaseSettings` with `env_prefix="BIOMARKER_"`.

### 6.1 Connection Settings

| Variable                    | Default       | Description                           |
|-----------------------------|---------------|---------------------------------------|
| `BIOMARKER_MILVUS_HOST`    | `localhost`   | Milvus server hostname                |
| `BIOMARKER_MILVUS_PORT`    | `19530`       | Milvus gRPC port                      |

### 6.2 Embedding Settings

| Variable                          | Default                 | Description                    |
|-----------------------------------|-------------------------|--------------------------------|
| `BIOMARKER_EMBEDDING_MODEL`      | `BAAI/bge-small-en-v1.5`| HuggingFace model name         |
| `BIOMARKER_EMBEDDING_DIMENSION`  | `384`                   | Vector dimension               |
| `BIOMARKER_EMBEDDING_BATCH_SIZE` | `32`                    | Batch size for bulk embedding  |

### 6.3 LLM Settings

| Variable                       | Default       | Description                          |
|--------------------------------|---------------|--------------------------------------|
| `BIOMARKER_LLM_PROVIDER`     | `anthropic`   | LLM provider name                    |
| `BIOMARKER_LLM_MODEL`        | `claude-sonnet-4-6` | Model identifier              |
| `BIOMARKER_ANTHROPIC_API_KEY`| *(required)*  | Anthropic API key                    |
| `BIOMARKER_LLM_MAX_RETRIES`  | `3`           | Max retry attempts for LLM calls     |

### 6.4 RAG Search Settings

| Variable                           | Default | Description                            |
|------------------------------------|---------|----------------------------------------|
| `BIOMARKER_TOP_K_PER_COLLECTION`  | `5`     | Results retrieved per collection       |
| `BIOMARKER_SCORE_THRESHOLD`       | `0.4`   | Minimum similarity score to include    |

### 6.5 Collection Search Weights

Weights control relative importance of each collection in multi-collection
search. They should sum to approximately 1.0. The settings validator emits
a warning if the sum deviates by more than 0.05 from 1.0.

| Variable                                  | Default | Collection                     |
|-------------------------------------------|---------|--------------------------------|
| `BIOMARKER_WEIGHT_BIOMARKER_REF`         | `0.12`  | `biomarker_reference`          |
| `BIOMARKER_WEIGHT_GENETIC_VARIANTS`      | `0.11`  | `biomarker_genetic_variants`   |
| `BIOMARKER_WEIGHT_PGX_RULES`             | `0.10`  | `biomarker_pgx_rules`          |
| `BIOMARKER_WEIGHT_DISEASE_TRAJECTORIES`  | `0.10`  | `biomarker_disease_trajectories`|
| `BIOMARKER_WEIGHT_CLINICAL_EVIDENCE`     | `0.09`  | `biomarker_clinical_evidence`  |
| `BIOMARKER_WEIGHT_GENOMIC_EVIDENCE`      | `0.08`  | `genomic_evidence` (shared)    |
| `BIOMARKER_WEIGHT_DRUG_INTERACTIONS`     | `0.07`  | `biomarker_drug_interactions`  |
| `BIOMARKER_WEIGHT_AGING_MARKERS`         | `0.07`  | `biomarker_aging_markers`      |
| `BIOMARKER_WEIGHT_NUTRITION`             | `0.05`  | `biomarker_nutrition`          |
| `BIOMARKER_WEIGHT_GENOTYPE_ADJUSTMENTS`  | `0.05`  | `biomarker_genotype_adjustments`|
| `BIOMARKER_WEIGHT_MONITORING`            | `0.05`  | `biomarker_monitoring`         |
| `BIOMARKER_WEIGHT_CRITICAL_VALUES`       | `0.04`  | `biomarker_critical_values`    |
| `BIOMARKER_WEIGHT_DISCORDANCE_RULES`     | `0.04`  | `biomarker_discordance_rules`  |
| `BIOMARKER_WEIGHT_AJ_CARRIER_SCREENING` | `0.03`  | `biomarker_aj_carrier_screening`|

**Total:** 1.00

### 6.6 Server Settings

| Variable                    | Default    | Description                             |
|-----------------------------|------------|-----------------------------------------|
| `BIOMARKER_API_HOST`       | `0.0.0.0`  | FastAPI bind address                    |
| `BIOMARKER_API_PORT`       | `8529`     | FastAPI listen port                      |
| `BIOMARKER_STREAMLIT_PORT` | `8528`     | Streamlit listen port                    |
| `BIOMARKER_METRICS_ENABLED`| `True`     | Enable Prometheus metrics endpoint       |

### 6.7 Conversation & Citation Settings

| Variable                                  | Default | Description                          |
|-------------------------------------------|---------|--------------------------------------|
| `BIOMARKER_MAX_CONVERSATION_CONTEXT`     | `3`     | Prior turns kept in conversation     |
| `BIOMARKER_CITATION_HIGH_THRESHOLD`      | `0.75`  | Score above this = high confidence   |
| `BIOMARKER_CITATION_MEDIUM_THRESHOLD`    | `0.60`  | Score above this = medium confidence |

### 6.8 Security Settings

| Variable                         | Default                      | Description                      |
|----------------------------------|------------------------------|----------------------------------|
| `BIOMARKER_CORS_ORIGINS`        | `http://localhost:8080,...`   | Comma-separated allowed origins  |
| `BIOMARKER_MAX_REQUEST_SIZE_MB` | `10`                         | Max request body size (MB)       |
| `BIOMARKER_API_KEY`             | *(empty = no auth)*          | API key for request auth         |

**Important:** When `BIOMARKER_API_KEY` is empty (the default), authentication
is disabled entirely. In any deployment handling patient health information
(PHI/PII), always set a strong API key. Authenticated requests must include
the `X-API-Key` header:

```bash
curl -H "X-API-Key: your-key" http://localhost:8529/v1/analyze
```

### 6.9 Timeout Settings

| Variable                               | Default | Description                         |
|----------------------------------------|---------|-------------------------------------|
| `BIOMARKER_REQUEST_TIMEOUT_SECONDS`   | `60`    | Overall request timeout              |
| `BIOMARKER_MILVUS_TIMEOUT_SECONDS`    | `10`    | Milvus query timeout                 |
| `BIOMARKER_LLM_MAX_RETRIES`           | `3`     | LLM call retry attempts             |

### 6.10 Complete `.env` Example

```bash
# ── Required ──
BIOMARKER_ANTHROPIC_API_KEY=sk-ant-api03-XXXXXXXXXXXXXXXXXXXXXXXX

# ── Milvus ──
BIOMARKER_MILVUS_HOST=localhost
BIOMARKER_MILVUS_PORT=19530

# ── Embedding ──
BIOMARKER_EMBEDDING_MODEL=BAAI/bge-small-en-v1.5
BIOMARKER_EMBEDDING_DIMENSION=384
BIOMARKER_EMBEDDING_BATCH_SIZE=32

# ── LLM ──
BIOMARKER_LLM_PROVIDER=anthropic
BIOMARKER_LLM_MODEL=claude-sonnet-4-6
BIOMARKER_LLM_MAX_RETRIES=3

# ── RAG ──
BIOMARKER_TOP_K_PER_COLLECTION=5
BIOMARKER_SCORE_THRESHOLD=0.4

# ── Weights (sum = 1.00) ──
BIOMARKER_WEIGHT_BIOMARKER_REF=0.12
BIOMARKER_WEIGHT_GENETIC_VARIANTS=0.11
BIOMARKER_WEIGHT_PGX_RULES=0.10
BIOMARKER_WEIGHT_DISEASE_TRAJECTORIES=0.10
BIOMARKER_WEIGHT_CLINICAL_EVIDENCE=0.09
BIOMARKER_WEIGHT_GENOMIC_EVIDENCE=0.08
BIOMARKER_WEIGHT_DRUG_INTERACTIONS=0.07
BIOMARKER_WEIGHT_AGING_MARKERS=0.07
BIOMARKER_WEIGHT_NUTRITION=0.05
BIOMARKER_WEIGHT_GENOTYPE_ADJUSTMENTS=0.05
BIOMARKER_WEIGHT_MONITORING=0.05
BIOMARKER_WEIGHT_CRITICAL_VALUES=0.04
BIOMARKER_WEIGHT_DISCORDANCE_RULES=0.04
BIOMARKER_WEIGHT_AJ_CARRIER_SCREENING=0.03

# ── Server ──
BIOMARKER_API_HOST=0.0.0.0
BIOMARKER_API_PORT=8529
BIOMARKER_STREAMLIT_PORT=8528
BIOMARKER_METRICS_ENABLED=True

# ── Conversation ──
BIOMARKER_MAX_CONVERSATION_CONTEXT=3

# ── Citations ──
BIOMARKER_CITATION_HIGH_THRESHOLD=0.75
BIOMARKER_CITATION_MEDIUM_THRESHOLD=0.60

# ── Security ──
BIOMARKER_CORS_ORIGINS=http://localhost:8080,http://localhost:8528,http://localhost:8529
BIOMARKER_MAX_REQUEST_SIZE_MB=10
BIOMARKER_API_KEY=
BIOMARKER_REQUEST_TIMEOUT_SECONDS=60
BIOMARKER_MILVUS_TIMEOUT_SECONDS=10
```

---

## 7. Milvus Collection Setup

### 7.1 Collection Inventory

The agent uses **13 biomarker-specific collections** plus **1 shared
read-only collection** from the RAG pipeline. All collections use
384-dimensional float vectors (BGE-small-en-v1.5).

| #  | Collection Name                   | Purpose                                 | Read-Only |
|----|-----------------------------------|-----------------------------------------|-----------|
| 1  | `biomarker_reference`             | Reference biomarker definitions & ranges| No        |
| 2  | `biomarker_genetic_variants`      | Genetic variants affecting biomarkers   | No        |
| 3  | `biomarker_pgx_rules`             | CPIC pharmacogenomic dosing rules       | No        |
| 4  | `biomarker_disease_trajectories`  | Disease progression trajectories        | No        |
| 5  | `biomarker_clinical_evidence`     | Published clinical evidence             | No        |
| 6  | `biomarker_nutrition`             | Genotype-aware nutrition guidelines     | No        |
| 7  | `biomarker_drug_interactions`     | Gene-drug interactions                  | No        |
| 8  | `biomarker_aging_markers`         | Epigenetic aging clock markers          | No        |
| 9  | `biomarker_genotype_adjustments`  | Genotype-based range adjustments        | No        |
| 10 | `biomarker_monitoring`            | Condition-specific monitoring protocols | No        |
| 11 | `biomarker_critical_values`       | Critical value alert thresholds         | No        |
| 12 | `biomarker_discordance_rules`     | Discordance detection rules             | No        |
| 13 | `biomarker_aj_carrier_screening`  | Ashkenazi Jewish carrier screening      | No        |
| 14 | `genomic_evidence`                | Shared genomic variant evidence         | **Yes**   |

### 7.2 Schema Details

Every collection shares a common structure:

- **Primary key:** `id` (VARCHAR, max 100 chars)
- **Embedding field:** `embedding` (FLOAT_VECTOR, dim=384)
- **Metadata fields:** Collection-specific (VARCHAR, INT64, FLOAT, JSON)
- **Index:** IVF_FLAT on the embedding field, L2 distance metric

### 7.3 Automated Collection Setup

The recommended approach uses the provided scripts:

```bash
# Create all 13 collections with proper schemas
python3 scripts/setup_collections.py --drop-existing

# Seed reference data into all collections
python3 scripts/seed_all.py
```

The `setup_collections.py` script:
1. Connects to Milvus at `BIOMARKER_MILVUS_HOST:BIOMARKER_MILVUS_PORT`
2. Drops existing collections if `--drop-existing` is passed
3. Creates each collection with the schema defined in `src/collections.py`
4. Builds IVF_FLAT indexes on each embedding field
5. Loads all collections into memory

The `seed_all.py` script:
1. Reads reference data from `data/reference/` (JSON files)
2. Generates embeddings using BGE-small-en-v1.5
3. Inserts vectorized records into each collection
4. Reports per-collection record counts

### 7.4 Manual Collection Verification

```bash
# Python: verify all collections exist and have data
python3 -c "
from pymilvus import connections, utility
connections.connect(host='localhost', port=19530)
collections = utility.list_collections()
print(f'Total collections: {len(collections)}')
for name in sorted(collections):
    if name.startswith('biomarker_') or name == 'genomic_evidence':
        from pymilvus import Collection
        c = Collection(name)
        c.load()
        print(f'  {name}: {c.num_entities} entities')
connections.disconnect('default')
"
```

### 7.5 Backup and Restore

```bash
# Export collection data (requires pymilvus bulk_insert utilities)
python3 -c "
from pymilvus import connections, Collection
connections.connect(host='localhost', port=19530)

import json
collection_names = [
    'biomarker_reference', 'biomarker_genetic_variants',
    'biomarker_pgx_rules', 'biomarker_disease_trajectories',
    'biomarker_clinical_evidence', 'biomarker_nutrition',
    'biomarker_drug_interactions', 'biomarker_aging_markers',
    'biomarker_genotype_adjustments', 'biomarker_monitoring',
    'biomarker_critical_values', 'biomarker_discordance_rules',
    'biomarker_aj_carrier_screening',
]

for name in collection_names:
    c = Collection(name)
    c.load()
    print(f'{name}: {c.num_entities} entities')
"

# Full Milvus backup via Docker volume snapshot
docker run --rm \
    -v precision_biomarker_agent_milvus_data:/data \
    -v $(pwd)/backups:/backup \
    alpine tar czf /backup/milvus-data-$(date +%Y%m%d).tar.gz /data
```

### 7.6 Shared `genomic_evidence` Collection

The `genomic_evidence` collection is created and managed by the main
RAG Chat Pipeline (`rag-chat-pipeline/`). The Biomarker Agent connects
to it in **read-only** mode:

- The agent never inserts, updates, or deletes records in `genomic_evidence`
- It must exist before the agent starts (otherwise the agent operates
  without genomic context but remains functional)
- Search weight: 0.08 (8% of total relevance)

If the RAG pipeline is not deployed, the agent gracefully degrades --
queries still work across the 13 biomarker-specific collections.

---

## 8. Security Hardening

### 8.1 API Key Authentication

Enable request authentication by setting `BIOMARKER_API_KEY`:

```bash
# Generate a strong API key
BIOMARKER_API_KEY=$(python3 -c "import secrets; print(secrets.token_urlsafe(32))")
echo "BIOMARKER_API_KEY=${BIOMARKER_API_KEY}" >> .env
echo "Save this key: ${BIOMARKER_API_KEY}"
```

With authentication enabled:
- All endpoints except `/health` and `/metrics` require the `X-API-Key` header
- Requests without a valid key receive HTTP 401
- The `/health` endpoint remains unauthenticated for load balancer probes

Client usage:

```bash
# Authenticated request
curl -H "X-API-Key: ${BIOMARKER_API_KEY}" \
     -H "Content-Type: application/json" \
     -d '{"query": "What is HbA1c?"}' \
     http://localhost:8529/v1/query
```

### 8.2 CORS Configuration

Restrict allowed origins in production:

```bash
# Production: only allow the landing page and specific frontends
BIOMARKER_CORS_ORIGINS=http://dgx-spark:8080,https://your-domain.com
```

The default (`http://localhost:8080,http://localhost:8528,http://localhost:8529`)
is suitable for development only.

### 8.3 Request Size Limits

The default limit is 10 MB per request. Adjust for your use case:

```bash
# Restrict to 2 MB for API-only deployments
BIOMARKER_MAX_REQUEST_SIZE_MB=2
```

Requests exceeding this limit receive HTTP 413 (Payload Too Large).

### 8.4 Network Isolation

In production, internal services should not be exposed externally:

```yaml
# docker-compose.dgx.yml
services:
  milvus-standalone:
    ports:
      - "127.0.0.1:19530:19530"    # Localhost only
      - "127.0.0.1:9091:9091"

  milvus-etcd:
    # No port mapping -- internal only

  milvus-minio:
    # No port mapping -- internal only
```

### 8.5 Non-Root Container Execution

The Dockerfile creates a dedicated `biomarkeruser` and runs the application
as that user:

```dockerfile
RUN useradd -r -s /bin/false biomarkeruser \
    && mkdir -p /app/data/cache /app/data/reference \
    && chown -R biomarkeruser:biomarkeruser /app
USER biomarkeruser
```

Verify the container runs as non-root:

```bash
docker exec biomarker-api whoami
# Expected: biomarkeruser
```

### 8.6 Secret Management

**Never** commit API keys to version control. Follow these practices:

```bash
# 1. Ensure .env is in .gitignore
echo ".env" >> .gitignore

# 2. Use environment variables directly (no .env file)
export BIOMARKER_ANTHROPIC_API_KEY="sk-ant-api03-..."

# 3. Or use Docker secrets (Swarm mode)
echo "sk-ant-api03-..." | docker secret create biomarker_api_key -

# 4. Or mount from a secrets manager
docker run -e BIOMARKER_ANTHROPIC_API_KEY_FILE=/run/secrets/api_key \
    -v /path/to/secret:/run/secrets/api_key:ro \
    biomarker-api
```

### 8.7 TLS Termination

For HTTPS, use a reverse proxy (nginx, Caddy, or Traefik) in front of the
FastAPI server. Do not terminate TLS in the application itself:

```nginx
# /etc/nginx/sites-available/biomarker-api
server {
    listen 443 ssl http2;
    server_name biomarker.your-domain.com;

    ssl_certificate     /etc/ssl/certs/biomarker.crt;
    ssl_certificate_key /etc/ssl/private/biomarker.key;

    location / {
        proxy_pass http://127.0.0.1:8529;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

### 8.8 HIPAA Compliance Considerations

When deploying with patient health information (PHI):

1. **Enable API authentication** (`BIOMARKER_API_KEY`)
2. **Enable TLS** for all external endpoints
3. **Restrict CORS origins** to known frontends
4. **Bind Milvus to localhost** (no external access)
5. **Enable audit logging** (the agent logs all queries via `src/audit.py`)
6. **Encrypt Docker volumes** at the host OS level
7. **Rotate API keys** periodically (minimum quarterly)
8. **Monitor access logs** for anomalous patterns

---

## 9. Monitoring & Health Checks

### 9.1 Health Endpoints

The FastAPI server exposes two health check endpoints:

| Endpoint    | Method | Auth Required | Purpose                         |
|-------------|--------|---------------|---------------------------------|
| `/health`   | GET    | No            | Full health with collection stats|
| `/healthz`  | GET    | No            | Kubernetes-style liveness probe |

#### `/health` Response

```bash
curl -s http://localhost:8529/health | python3 -m json.tool
```

```json
{
    "status": "healthy",
    "collections": 14,
    "total_vectors": 2847,
    "agent_ready": true
}
```

Possible `status` values:

| Status     | Meaning                                            |
|------------|----------------------------------------------------|
| `healthy`  | All systems operational                            |
| `degraded` | Milvus unreachable; agent can't answer queries     |

#### HTTP Status Codes

| Code | Meaning                         | Action                     |
|------|---------------------------------|----------------------------|
| 200  | Healthy                         | None                       |
| 503  | Milvus unavailable              | Check Milvus container     |

### 9.2 Milvus Health

```bash
# Milvus built-in health endpoint
curl -s http://localhost:9091/healthz
# Expected: OK
```

### 9.3 Docker Health Checks

All services have Docker-level health checks configured:

```bash
# View health status of all containers
docker compose ps

# Expected output:
# NAME                        STATUS                    PORTS
# biomarker-milvus-etcd       Up (healthy)
# biomarker-milvus-minio      Up (healthy)
# biomarker-milvus-standalone Up (healthy)              0.0.0.0:19530->19530, 0.0.0.0:9091->9091
# biomarker-streamlit         Up (healthy)              0.0.0.0:8528->8528
# biomarker-api               Up (healthy)              0.0.0.0:8529->8529
# biomarker-setup             Exited (0)
```

### 9.4 Prometheus Metrics

The `/metrics` endpoint returns Prometheus-compatible metrics:

```bash
curl -s http://localhost:8529/metrics
```

```
# HELP biomarker_api_requests_total Total API requests
# TYPE biomarker_api_requests_total counter
biomarker_api_requests_total 142

# HELP biomarker_api_query_requests_total Total /query requests
# TYPE biomarker_api_query_requests_total counter
biomarker_api_query_requests_total 38

# HELP biomarker_api_search_requests_total Total /search requests
# TYPE biomarker_api_search_requests_total counter
biomarker_api_search_requests_total 12

# HELP biomarker_api_analyze_requests_total Total /analyze requests
# TYPE biomarker_api_analyze_requests_total counter
biomarker_api_analyze_requests_total 27

# HELP biomarker_api_bio_age_requests_total Total /biological-age requests
# TYPE biomarker_api_bio_age_requests_total counter
biomarker_api_bio_age_requests_total 15

# HELP biomarker_api_errors_total Total error responses
# TYPE biomarker_api_errors_total counter
biomarker_api_errors_total 3

# HELP biomarker_collection_vectors Number of vectors per collection
# TYPE biomarker_collection_vectors gauge
biomarker_collection_vectors{collection="biomarker_reference"} 247
biomarker_collection_vectors{collection="biomarker_genetic_variants"} 183
...
```

#### Grafana Dashboard Integration

Add the metrics endpoint to your Prometheus `scrape_configs`:

```yaml
# prometheus.yml
scrape_configs:
  - job_name: "biomarker-agent"
    scrape_interval: 30s
    static_configs:
      - targets: ["localhost:8529"]
    metrics_path: "/metrics"
```

### 9.5 Automated Health Monitoring Script

Create a monitoring cron job:

```bash
#!/bin/bash
# /usr/local/bin/biomarker-health-check.sh
set -euo pipefail

API_URL="${BIOMARKER_API_URL:-http://localhost:8529}"
ALERT_EMAIL="${ALERT_EMAIL:-admin@your-domain.com}"

# Check FastAPI health
HTTP_CODE=$(curl -s -o /dev/null -w "%{http_code}" "${API_URL}/health" || echo "000")

if [ "$HTTP_CODE" != "200" ]; then
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] ALERT: Biomarker API health check failed (HTTP ${HTTP_CODE})" \
        | tee -a /var/log/biomarker-health.log

    # Attempt container restart
    cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
    docker compose restart biomarker-api

    # Optional: send email alert
    # echo "Biomarker API returned HTTP ${HTTP_CODE}" | mail -s "Biomarker Alert" "${ALERT_EMAIL}"
fi

# Check Milvus health
MILVUS_STATUS=$(curl -sf http://localhost:9091/healthz || echo "FAIL")
if [ "$MILVUS_STATUS" != "OK" ]; then
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] ALERT: Milvus health check failed" \
        | tee -a /var/log/biomarker-health.log
    cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
    docker compose restart milvus-standalone
fi
```

Add to crontab:

```bash
# Run every 5 minutes
*/5 * * * * /usr/local/bin/biomarker-health-check.sh >> /var/log/biomarker-health.log 2>&1
```

### 9.6 Container Log Management

```bash
# Follow all service logs
docker compose logs -f

# Follow a specific service
docker compose logs -f biomarker-api

# Last 100 lines from the API server
docker compose logs --tail 100 biomarker-api

# Export logs for analysis
docker compose logs --no-color biomarker-api > /tmp/biomarker-api.log 2>&1
```

### 9.7 Knowledge Graph Statistics

```bash
curl -s http://localhost:8529/knowledge/stats | python3 -m json.tool
```

```json
{
    "disease_domains": 13,
    "total_biomarkers": 201,
    "total_genetic_modifiers": 171,
    "pharmacogenes": 42,
    "pgx_drug_interactions": 87,
    "phenoage_markers": 12,
    "cross_modal_links": 634
}
```

### 9.8 End-to-End Validation

```bash
# Run the built-in end-to-end validation script
python3 scripts/validate_e2e.py
```

This verifies:
1. Milvus connectivity and collection availability
2. Embedding model loading
3. LLM connectivity (Anthropic API)
4. Multi-collection search and ranking
5. Response generation pipeline

---

## 10. Troubleshooting

### 10.1 Milvus Will Not Start

**Symptom:** `milvus-standalone` container exits or stays unhealthy.

```bash
# Check Milvus logs
docker compose logs milvus-standalone

# Check etcd logs (Milvus depends on etcd)
docker compose logs milvus-etcd

# Check MinIO logs
docker compose logs milvus-minio
```

**Common causes:**

| Cause                          | Solution                                    |
|--------------------------------|---------------------------------------------|
| etcd not healthy               | `docker compose restart milvus-etcd`        |
| MinIO not healthy              | `docker compose restart milvus-minio`       |
| Port 19530 already in use      | `lsof -i :19530` and stop conflicting service|
| Insufficient memory            | Ensure at least 8 GB RAM available          |
| Corrupted data volume          | `docker compose down -v && docker compose up -d` |

```bash
# Full reset (destroys all data)
docker compose down -v
docker compose up -d
# Re-run setup after Milvus is healthy
docker compose logs -f biomarker-setup
```

### 10.2 Collections Not Created

**Symptom:** `/health` returns `collections: 0` or `biomarker-setup` exited
with non-zero code.

```bash
# Check setup logs
docker compose logs biomarker-setup

# Re-run setup manually
docker compose run --rm biomarker-setup

# Or run locally
python3 scripts/setup_collections.py --drop-existing
python3 scripts/seed_all.py
```

### 10.3 Embedding Model Download Fails

**Symptom:** `biomarker-api` or `biomarker-streamlit` hangs on startup, logs
show HuggingFace download errors.

```bash
# Pre-download the model before starting containers
python3 -c "
from sentence_transformers import SentenceTransformer
model = SentenceTransformer('BAAI/bge-small-en-v1.5')
print(f'Model loaded: {model.get_sentence_embedding_dimension()} dimensions')
"
```

**Common causes:**

| Cause                         | Solution                                      |
|-------------------------------|-----------------------------------------------|
| No internet access            | Pre-download model, mount as volume           |
| HuggingFace rate limit        | Set `HF_TOKEN` environment variable           |
| Disk full                     | Free space in Docker data root                |
| Proxy required                | Set `HTTP_PROXY` / `HTTPS_PROXY` env vars     |

**Air-gapped deployment:** Pre-download the model and mount it:

```bash
# On a machine with internet
python3 -c "
from sentence_transformers import SentenceTransformer
model = SentenceTransformer('BAAI/bge-small-en-v1.5')
model.save('/tmp/bge-small-en-v1.5')
"

# Copy /tmp/bge-small-en-v1.5 to the target machine
# Mount in docker-compose.yml:
#   volumes:
#     - /path/to/bge-small-en-v1.5:/app/data/cache/bge-small-en-v1.5
# Set:
#   BIOMARKER_EMBEDDING_MODEL=/app/data/cache/bge-small-en-v1.5
```

### 10.4 Anthropic API Errors

**Symptom:** Queries return errors about LLM being unavailable.

```bash
# Verify the API key is set
docker exec biomarker-api env | grep -i anthropic
# Should show ANTHROPIC_API_KEY=sk-ant-...

# Test connectivity directly
curl -s https://api.anthropic.com/v1/messages \
    -H "x-api-key: ${BIOMARKER_ANTHROPIC_API_KEY}" \
    -H "anthropic-version: 2023-06-01" \
    -H "content-type: application/json" \
    -d '{"model":"claude-sonnet-4-6","max_tokens":10,"messages":[{"role":"user","content":"hi"}]}'
```

**Common causes:**

| Cause                          | Solution                                    |
|--------------------------------|---------------------------------------------|
| API key not set                | Add `BIOMARKER_ANTHROPIC_API_KEY` to `.env` |
| API key invalid                | Regenerate at console.anthropic.com         |
| Rate limit exceeded            | Reduce request frequency, check plan limits |
| Network blocked                | Ensure outbound HTTPS to api.anthropic.com  |
| Wrong model name               | Verify `BIOMARKER_LLM_MODEL` value         |

**Note:** The agent degrades gracefully without LLM access. RAG search
results are still returned, but without synthesized natural-language answers.

### 10.5 Port Conflicts

**Symptom:** `docker compose up` fails with "address already in use."

```bash
# Find what is using the port
sudo lsof -i :8528
sudo lsof -i :8529
sudo lsof -i :19530

# Kill the conflicting process
sudo kill -9 <PID>

# Or change the port mapping in docker-compose.yml
# ports:
#   - "8530:8528"  # Use 8530 instead
```

### 10.6 Container Runs Out of Memory

**Symptom:** Container exits with OOM kill (exit code 137).

```bash
# Check Docker memory limits
docker stats --no-stream

# Increase memory limits in docker-compose.yml
# deploy:
#   resources:
#     limits:
#       memory: 8G
```

**Memory guidelines by service:**

| Service             | Minimum | Recommended |
|---------------------|---------|-------------|
| milvus-standalone   | 4 GB    | 8 GB        |
| biomarker-api       | 2 GB    | 4 GB        |
| biomarker-streamlit | 1 GB    | 2 GB        |
| milvus-etcd         | 512 MB  | 1 GB        |
| milvus-minio        | 512 MB  | 1 GB        |

### 10.7 Slow Query Performance

**Symptom:** Queries take more than 10 seconds.

```bash
# Check collection load status
python3 -c "
from pymilvus import connections, Collection
connections.connect(host='localhost', port=19530)
names = [
    'biomarker_reference', 'biomarker_genetic_variants',
    'biomarker_pgx_rules', 'biomarker_disease_trajectories',
    'biomarker_clinical_evidence', 'biomarker_nutrition',
    'biomarker_drug_interactions', 'biomarker_aging_markers',
    'biomarker_genotype_adjustments', 'biomarker_monitoring',
    'biomarker_critical_values', 'biomarker_discordance_rules',
    'biomarker_aj_carrier_screening', 'genomic_evidence',
]
for name in names:
    try:
        c = Collection(name)
        c.load()
        print(f'{name}: loaded, {c.num_entities} entities')
    except Exception as e:
        print(f'{name}: ERROR - {e}')
"
```

**Common causes:**

| Cause                            | Solution                                   |
|----------------------------------|--------------------------------------------|
| Collections not loaded to memory | Run `collection.load()` for each           |
| Missing indexes                  | Re-run `setup_collections.py`              |
| Too many concurrent requests     | Increase `--workers` or add replicas       |
| `TOP_K_PER_COLLECTION` too high | Reduce to 3 for faster results             |

### 10.8 Tests Failing

```bash
# Run tests with full output
python3 -m pytest tests/ -v --tb=long

# Run only unit tests (no external dependencies)
python3 -m pytest tests/ -v -m "not integration and not slow"

# Check test environment
python3 -c "
import pymilvus; print(f'pymilvus: {pymilvus.__version__}')
import anthropic; print(f'anthropic: {anthropic.__version__}')
import sentence_transformers; print(f'sentence-transformers: {sentence_transformers.__version__}')
import fastapi; print(f'fastapi: {fastapi.__version__}')
import streamlit; print(f'streamlit: {streamlit.__version__}')
"
```

### 10.9 Docker Build Failures

**Symptom:** `docker compose build` fails during pip install.

```bash
# Build with no cache (forces fresh dependency install)
docker compose build --no-cache

# Build a single service for faster debugging
docker compose build biomarker-api

# Check disk space (Docker build needs working space)
df -h /var/lib/docker
```

**Common causes:**

| Cause                          | Solution                                    |
|--------------------------------|---------------------------------------------|
| Pip version too old            | Dockerfile already runs `pip install --upgrade pip` |
| Missing system library         | Check Dockerfile `apt-get install` stage    |
| Network timeout during build   | Retry, or configure pip proxy               |
| Docker disk full               | `docker system prune -a` (removes unused images) |

### 10.10 Upgrading the Agent

```bash
# 1. Pull latest code
git pull origin main

# 2. Rebuild containers
docker compose build --no-cache

# 3. Stop old containers
docker compose down

# 4. Re-run setup (preserves Milvus data unless --drop-existing)
docker compose up -d

# 5. If collection schemas changed, re-create
docker compose run --rm biomarker-setup

# 6. Verify
curl -s http://localhost:8529/health | python3 -m json.tool
```

### 10.11 Complete Reset

When all else fails, perform a clean reset:

```bash
# Stop everything and remove all volumes
docker compose down -v

# Remove built images
docker compose down --rmi local

# Remove any orphaned containers
docker container prune -f

# Rebuild from scratch
docker compose up -d --build

# Watch setup
docker compose logs -f biomarker-setup

# Verify
curl -s http://localhost:8529/health | python3 -m json.tool
```

### 10.12 Diagnostic Commands Reference

Quick reference for common diagnostic tasks:

```bash
# ── Container Status ──
docker compose ps
docker stats --no-stream

# ── Logs ──
docker compose logs -f biomarker-api         # API server logs
docker compose logs -f biomarker-streamlit   # UI logs
docker compose logs -f milvus-standalone     # Milvus logs
docker compose logs biomarker-setup          # Setup/seed output

# ── Shell Access ──
docker exec -it biomarker-api bash           # API container shell
docker exec -it biomarker-milvus-standalone bash  # Milvus shell

# ── Network ──
docker network inspect biomarker-network     # Network details
docker exec biomarker-api curl -s http://milvus-standalone:19530  # Internal connectivity

# ── Health ──
curl -s http://localhost:8529/health         # API health
curl -s http://localhost:9091/healthz         # Milvus health
curl -s http://localhost:8529/collections     # Collection inventory
curl -s http://localhost:8529/metrics         # Prometheus metrics
curl -s http://localhost:8529/knowledge/stats # Knowledge stats

# ── Disk ──
docker system df                             # Docker disk usage
du -sh /var/lib/docker/volumes/              # Volume sizes

# ── Tests ──
python3 -m pytest tests/ -v                  # Full test suite (709 tests)
python3 -m pytest tests/ -v -x              # Stop on first failure
python3 -m pytest tests/ --cov=src --cov=api --cov-report=term-missing  # Coverage
```

---

## Appendix A: API Endpoint Reference

| Method | Path                 | Auth  | Description                          |
|--------|----------------------|-------|--------------------------------------|
| GET    | `/`                  | No    | Service info and links               |
| GET    | `/health`            | No    | Health with collection stats         |
| GET    | `/healthz`           | No    | Liveness probe                       |
| GET    | `/collections`       | Yes*  | Collection names and counts          |
| GET    | `/knowledge/stats`   | Yes*  | Knowledge graph statistics           |
| GET    | `/metrics`           | No    | Prometheus metrics                   |
| GET    | `/docs`              | Yes*  | OpenAPI interactive docs             |
| GET    | `/openapi.json`      | Yes*  | OpenAPI schema                       |
| POST   | `/v1/query`          | Yes*  | RAG Q&A query                        |
| POST   | `/v1/analyze`        | Yes*  | Full patient analysis                |
| POST   | `/v1/biological-age` | Yes*  | Biological age calculation           |
| POST   | `/v1/disease-risk`   | Yes*  | Disease trajectory analysis          |
| POST   | `/v1/pgx`            | Yes*  | Pharmacogenomic mapping              |
| GET    | `/v1/health`         | Yes*  | V1 versioned health check            |

*Auth required only when `BIOMARKER_API_KEY` is set (via `X-API-Key` header).

---

## Appendix B: File Structure

```
precision_biomarker_agent/
├── .env.example                 # Environment variable template
├── .streamlit/                  # Streamlit configuration
├── api/
│   ├── __init__.py
│   ├── main.py                  # FastAPI application entrypoint
│   └── routes/
│       ├── __init__.py
│       ├── analysis.py          # /v1/analyze, /v1/biological-age, etc.
│       ├── events.py            # Event streaming endpoints
│       └── reports.py           # Report generation endpoints
├── app/
│   ├── __init__.py
│   ├── biomarker_ui.py          # Main Streamlit UI
│   ├── patient_360.py           # Patient 360 view component
│   └── protein_viewer.py        # Protein structure viewer
├── config/
│   └── settings.py              # PrecisionBiomarkerSettings (Pydantic)
├── data/
│   ├── events/                  # Event data
│   └── reference/               # Reference data for seeding
├── docker-compose.yml           # Full stack compose file
├── Dockerfile                   # Multi-stage build (builder + runtime)
├── docs/
│   └── DEPLOYMENT_GUIDE.md      # This file
├── pyproject.toml               # Project metadata, pytest, ruff config
├── README.md                    # Project README
├── requirements.txt             # Python dependencies
├── scripts/
│   ├── demo_validation.py       # Demo validation script
│   ├── expand_biomarker_reference.py
│   ├── expand_variants_and_interactions.py
│   ├── gen_critical_values.py
│   ├── gen_lab_ranges_and_aj.py
│   ├── gen_patient_data.py
│   ├── seed_all.py              # Seed all collections with reference data
│   ├── setup_collections.py     # Create Milvus collections with schemas
│   └── validate_e2e.py          # End-to-end validation
├── src/
│   ├── __init__.py
│   ├── agent.py                 # PrecisionBiomarkerAgent orchestrator
│   ├── audit.py                 # Audit logging
│   ├── biological_age.py        # Biological age calculator
│   ├── collections.py           # BiomarkerCollectionManager (Milvus)
│   ├── critical_values.py       # Critical value alerting
│   ├── discordance_detector.py  # Biomarker discordance detection
│   ├── disease_trajectory.py    # Disease trajectory analyzer
│   ├── export.py                # PDF/data export
│   ├── genotype_adjustment.py   # Genotype-adjusted reference ranges
│   ├── knowledge.py             # Knowledge graph statistics
│   ├── lab_range_interpreter.py # Lab result interpretation
│   ├── models.py                # Pydantic data models
│   ├── pharmacogenomics.py      # PGx drug-gene mapping
│   ├── rag_engine.py            # Multi-collection RAG engine
│   ├── report_generator.py      # Report generation
│   └── translation.py          # Clinical language translation
└── tests/
    ├── conftest.py              # Shared fixtures
    ├── __init__.py
    ├── test_agent.py
    ├── test_api.py
    ├── test_biological_age.py
    ├── test_collections.py
    ├── test_critical_values.py
    ├── test_discordance_detector.py
    ├── test_disease_trajectory.py
    ├── test_edge_cases.py
    ├── test_export.py
    ├── test_genotype_adjustment.py
    ├── test_integration.py
    ├── test_lab_range_interpreter.py
    ├── test_longitudinal.py
    ├── test_models.py
    ├── test_pharmacogenomics.py
    ├── test_rag_engine.py
    ├── test_report_generator.py
    └── test_ui.py
```

---

## Appendix C: Version History

| Version | Date       | Changes                           |
|---------|------------|-----------------------------------|
| 1.0.0   | March 2026 | Initial deployment guide release  |

---

*Precision Biomarker Intelligence Agent -- HCLS AI Factory*
*Author: Adam Jones | License: Apache 2.0*
