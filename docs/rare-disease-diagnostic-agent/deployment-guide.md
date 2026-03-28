# Rare Disease Diagnostic Agent -- Deployment Guide

**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Quick Start (Docker Compose)](#2-quick-start-docker-compose)
3. [Integrated Deployment (HCLS AI Factory)](#3-integrated-deployment-hcls-ai-factory)
4. [Standalone Deployment](#4-standalone-deployment)
5. [Manual Setup](#5-manual-setup)
6. [Milvus Configuration](#6-milvus-configuration)
7. [Collection Setup](#7-collection-setup)
8. [Data Ingestion](#8-data-ingestion)
9. [Environment Variables](#9-environment-variables)
10. [Security Configuration](#10-security-configuration)
11. [Health Checks](#11-health-checks)
12. [Monitoring](#12-monitoring)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Prerequisites

### Hardware

| Component | Minimum | Recommended |
|---|---|---|
| CPU | 8 cores | 16+ cores |
| RAM | 16 GB | 32+ GB |
| GPU | NVIDIA GPU (CUDA 12.x) | NVIDIA DGX Spark |
| Disk | 20 GB | 50+ GB |

### Software

| Dependency | Version | Purpose |
|---|---|---|
| Python | 3.10+ | Runtime |
| Docker | 24.0+ | Containerization |
| Docker Compose | 2.20+ | Multi-container orchestration |
| CUDA | 12.x | GPU acceleration |
| pip | Latest | Package management |

### API Keys

| Key | Required | Purpose |
|---|---|---|
| ANTHROPIC_API_KEY | Optional | Claude LLM for RAG synthesis |
| ORPHANET_API_KEY | Optional | Orphanet data ingestion |
| NCBI_API_KEY | Optional | NCBI/OMIM data ingestion |

---

## 2. Quick Start (Docker Compose)

The fastest path to a running system:

```bash
# Clone and navigate
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/rare_disease_diagnostic_agent

# Set API key
export ANTHROPIC_API_KEY="your-key-here"

# Start all services
docker compose up -d

# Verify
curl http://localhost:8134/health
```

This starts:
- FastAPI API on port 8134
- Streamlit UI on port 8544
- Milvus on port 19530 (with etcd and MinIO)

### Docker Compose File

```yaml
version: "3.8"

services:
  # --- Milvus Dependencies ---
  etcd:
    image: quay.io/coreos/etcd:v3.5.5
    environment:
      ETCD_AUTO_COMPACTION_MODE: revision
      ETCD_AUTO_COMPACTION_RETENTION: "1000"
      ETCD_QUOTA_BACKEND_BYTES: "4294967296"
    command: >
      etcd
      --advertise-client-urls=http://127.0.0.1:2379
      --listen-client-urls=http://0.0.0.0:2379

  minio:
    image: minio/minio:RELEASE.2023-03-20T20-16-18Z
    environment:
      MINIO_ACCESS_KEY: minioadmin
      MINIO_SECRET_KEY: minioadmin
    command: minio server /minio_data

  # --- Milvus Vector Database ---
  milvus:
    image: milvusdb/milvus:v2.3.3
    environment:
      ETCD_ENDPOINTS: etcd:2379
      MINIO_ADDRESS: minio:9000
    ports:
      - "19530:19530"
    depends_on:
      - etcd
      - minio

  # --- Rare Disease API ---
  rare-disease-api:
    build: .
    ports:
      - "8134:8134"
    environment:
      RD_MILVUS_HOST: milvus
      RD_MILVUS_PORT: "19530"
      RD_ANTHROPIC_API_KEY: ${ANTHROPIC_API_KEY}
      RD_API_PORT: "8134"
    command: >
      uvicorn api.main:app --host 0.0.0.0 --port 8134
    depends_on:
      - milvus

  # --- Rare Disease UI ---
  rare-disease-ui:
    build: .
    ports:
      - "8544:8544"
    environment:
      RD_API_BASE: http://rare-disease-api:8134
    command: >
      streamlit run app/diagnostic_ui.py
      --server.port 8544
      --server.address 0.0.0.0
    depends_on:
      - rare-disease-api
```

---

## 3. Integrated Deployment (HCLS AI Factory)

When deployed as part of the full HCLS AI Factory platform:

```bash
# From the platform root
cd /home/adam/projects/hcls-ai-factory

# Start the full platform (includes all agents)
./start-factory.sh

# Or start just the rare disease agent
docker compose -f docker-compose.dgx-spark.yml up -d \
  rare-disease-agent-api rare-disease-agent-ui
```

### Port Map (Integrated)

| Service | Port | Notes |
|---|---|---|
| Rare Disease API | 8134 | Unique to this agent |
| Rare Disease UI | 8544 | Unique to this agent |
| Milvus | 19530 | Shared with other agents |
| Landing Page | 8080 | Platform hub |
| Genomics API | 8527 | Cross-agent integration |
| PGx Agent | 8107 | Cross-agent integration |
| Cardiology Agent | 8126 | Cross-agent integration |
| Biomarker Agent | 8529 | Cross-agent integration |
| Clinical Trial Agent | 8538 | Cross-agent integration |

---

## 4. Standalone Deployment

To deploy the agent independently of the full platform:

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/rare_disease_diagnostic_agent

# Install dependencies
pip install -r requirements.txt

# Start Milvus (standalone)
docker compose up -d milvus etcd minio

# Wait for Milvus to be ready
sleep 10

# Setup collections
python scripts/setup_collections.py

# Seed knowledge base
python scripts/seed_knowledge.py

# Start API
uvicorn api.main:app --host 0.0.0.0 --port 8134 &

# Start UI
streamlit run app/diagnostic_ui.py --server.port 8544 &
```

---

## 5. Manual Setup

### 5.1 Python Environment

```bash
# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "from src.knowledge import KNOWLEDGE_VERSION; print(KNOWLEDGE_VERSION)"
```

### 5.2 Configuration

Create a `.env` file in the project root:

```bash
# Required for LLM features
RD_ANTHROPIC_API_KEY=sk-ant-...

# Milvus connection
RD_MILVUS_HOST=localhost
RD_MILVUS_PORT=19530

# API configuration
RD_API_PORT=8134
RD_STREAMLIT_PORT=8544
RD_API_KEY=your-api-key-here

# Embedding model (auto-downloaded)
RD_EMBEDDING_MODEL=BAAI/bge-small-en-v1.5

# LLM model
RD_LLM_MODEL=claude-sonnet-4-6
```

### 5.3 Run Tests

```bash
# Run full test suite
python -m pytest tests/ -v

# Expected output: 193 passed, 0 failed (0.16s)
```

---

## 6. Milvus Configuration

### 6.1 Resource Allocation

| Parameter | Demo | Production |
|---|---|---|
| Memory | 4 GB | 16+ GB |
| Disk | 10 GB | 100+ GB |
| etcd Quota | 4 GB | 8 GB |

### 6.2 Index Tuning

All 14 collections use IVF_FLAT with COSINE metric. For larger deployments:

| Parameter | Demo | Production |
|---|---|---|
| nlist | 128 | 1024 |
| nprobe (search) | 10 | 32 |
| Index Type | IVF_FLAT | IVF_SQ8 or HNSW |

### 6.3 Persistence

Milvus persists data through MinIO (object storage) and etcd (metadata). Ensure these volumes are mounted to persistent storage:

```yaml
volumes:
  milvus_data:
    driver: local
  etcd_data:
    driver: local
  minio_data:
    driver: local
```

---

## 7. Collection Setup

### 7.1 Create Collections

```bash
python scripts/setup_collections.py
```

This creates all 14 collections with schemas:
- rd_phenotypes, rd_diseases, rd_genes, rd_variants
- rd_literature, rd_trials, rd_therapies, rd_case_reports
- rd_guidelines, rd_pathways, rd_registries
- rd_natural_history, rd_newborn_screening, genomic_evidence

### 7.2 Verify Collections

```bash
curl http://localhost:8134/collections
```

---

## 8. Data Ingestion

### 8.1 Seed Knowledge

```bash
# Seed curated knowledge base
python scripts/seed_knowledge.py
```

### 8.2 Full Ingest

```bash
# Run all ingest parsers (HPO, OMIM, Orphanet, gene therapies)
python scripts/run_ingest.py
```

### 8.3 Verify Data

```bash
curl http://localhost:8134/health
# Check vector_counts for each collection
```

---

## 9. Environment Variables

### 9.1 Complete Variable Reference

| Variable | Default | Description |
|---|---|---|
| RD_MILVUS_HOST | localhost | Milvus hostname |
| RD_MILVUS_PORT | 19530 | Milvus port |
| RD_API_HOST | 0.0.0.0 | API bind address |
| RD_API_PORT | 8134 | API port |
| RD_STREAMLIT_PORT | 8544 | Streamlit port |
| RD_ANTHROPIC_API_KEY | (none) | Claude API key |
| RD_EMBEDDING_MODEL | BAAI/bge-small-en-v1.5 | Embedding model |
| RD_EMBEDDING_DIMENSION | 384 | Embedding dimension |
| RD_EMBEDDING_BATCH_SIZE | 32 | Embedding batch size |
| RD_LLM_PROVIDER | anthropic | LLM provider |
| RD_LLM_MODEL | claude-sonnet-4-6 | LLM model |
| RD_SCORE_THRESHOLD | 0.4 | Minimum similarity score |
| RD_API_KEY | (empty) | API authentication key |
| RD_CORS_ORIGINS | localhost:8080,8134,8544 | CORS allowed origins |
| RD_CROSS_AGENT_TIMEOUT | 30 | Cross-agent HTTP timeout (seconds) |
| RD_INGEST_SCHEDULE_HOURS | 24 | Ingest schedule interval |
| RD_INGEST_ENABLED | false | Enable scheduled ingestion |
| RD_MAX_CONVERSATION_CONTEXT | 3 | Conversation memory depth |
| RD_MAX_REQUEST_SIZE_MB | 10 | Maximum request body size |
| RD_METRICS_ENABLED | true | Enable Prometheus metrics |

---

## 10. Security Configuration

### 10.1 API Authentication

Enable API key authentication:

```bash
export RD_API_KEY="your-secure-api-key"
```

Clients must include the key in requests:

```bash
curl -H "X-API-Key: your-secure-api-key" http://localhost:8134/v1/diagnostic/query
```

### 10.2 CORS

Configure allowed origins:

```bash
export RD_CORS_ORIGINS="http://localhost:8080,http://localhost:8544,https://your-domain.com"
```

### 10.3 Network Security

- Restrict Milvus port (19530) to internal network only
- Use reverse proxy (nginx/traefik) for production HTTPS termination
- Enable API key authentication for all external access
- Cross-agent communication should use internal Docker network

---

## 11. Health Checks

### 11.1 API Health

```bash
curl http://localhost:8134/health

# Response:
{
  "status": "healthy",
  "collections": ["rd_phenotypes", "rd_diseases", ...],
  "vector_counts": {"rd_phenotypes": 18000, ...},
  "uptime": "2h 15m"
}
```

### 11.2 Docker Health Check

```yaml
healthcheck:
  test: ["CMD", "curl", "-f", "http://localhost:8134/health"]
  interval: 30s
  timeout: 10s
  retries: 3
  start_period: 60s
```

---

## 12. Monitoring

### 12.1 Prometheus Metrics

Access metrics at:

```bash
curl http://localhost:8134/metrics
```

### 12.2 Key Metrics

| Metric | Type | Description |
|---|---|---|
| rd_queries_total | Counter | Total queries received |
| rd_query_duration_seconds | Histogram | Query processing time |
| rd_workflow_executions_total | Counter | Workflow executions |
| rd_errors_total | Counter | Error count |
| rd_active_connections | Gauge | Active connections |

---

## 13. Troubleshooting

### 13.1 Common Issues

| Issue | Cause | Solution |
|---|---|---|
| API returns 503 | Milvus not ready | Wait for Milvus startup, check `docker logs milvus` |
| Empty search results | Collections not populated | Run `python scripts/seed_knowledge.py` |
| LLM timeout | Anthropic API issue | Check API key, network connectivity |
| Port conflict | Port already in use | Change RD_API_PORT or RD_STREAMLIT_PORT |
| Import errors | Missing dependencies | Run `pip install -r requirements.txt` |

### 13.2 Log Locations

| Log | Location | Command |
|---|---|---|
| API logs | stdout | `docker logs rare-disease-api` |
| Milvus logs | stdout | `docker logs milvus` |
| Streamlit logs | stdout | `docker logs rare-disease-ui` |

### 13.3 Reset Collections

```bash
# Warning: destroys all data
python scripts/setup_collections.py --force
python scripts/seed_knowledge.py
```

---

*Apache 2.0 License -- HCLS AI Factory*
