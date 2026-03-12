# Precision Autoimmune Intelligence Agent -- Deployment Guide

**Author:** Adam Jones
**Date:** March 2026

---

## Table of Contents

1. [Deployment Overview](#1-deployment-overview)
2. [Prerequisites](#2-prerequisites)
3. [Local Development Setup](#3-local-development-setup)
4. [Docker Deployment](#4-docker-deployment)
5. [DGX Spark Deployment](#5-dgx-spark-deployment)
6. [Environment Variables](#6-environment-variables)
7. [Milvus Configuration](#7-milvus-configuration)
8. [Health Checks and Monitoring](#8-health-checks-and-monitoring)
9. [Security Configuration](#9-security-configuration)
10. [Production Hardening](#10-production-hardening)
11. [Scaling and High Availability](#11-scaling-and-high-availability)
12. [Troubleshooting](#12-troubleshooting)
13. [Quick Reference](#13-quick-reference)

---

## 1. Deployment Overview

The Precision Autoimmune Intelligence Agent is a multi-collection RAG (Retrieval-Augmented Generation) system purpose-built for autoimmune disease analysis. It combines 14 domain-specific Milvus vector collections with Claude LLM synthesis to provide evidence-based clinical decision support covering differential diagnosis, flare prediction, pharmacogenomics, and treatment optimization.

### Architecture Summary

```
                         +---------------------+
                         |   Streamlit UI      |
                         |   (port 8531)       |
                         +--------+------------+
                                  |
                         +--------v------------+
                         |   FastAPI Server    |
                         |   (port 8532)       |
                         +--------+------------+
                                  |
                    +-------------+-------------+
                    |             |              |
            +-------v--+  +------v-----+  +----v--------+
            | RAG      |  | Diagnostic |  | Document    |
            | Engine   |  | Engine     |  | Processor   |
            +-------+--+  +------+-----+  +----+--------+
                    |             |              |
                    +------+------+--------------+
                           |
                  +--------v---------+
                  |  Milvus Vector   |
                  |  DB (19530)      |
                  |  14 collections  |
                  |  384-dim BGE     |
                  +------------------+
```

### Service Topology

The agent deploys as 3 Docker services plus a shared Milvus dependency:

| Service | Container Name | Port | Description |
|---|---|---|---|
| **Streamlit UI** | `autoimmune-streamlit` | 8531 | Interactive clinical interface |
| **FastAPI API** | `autoimmune-api` | 8532 | REST API with RAG query, search, ingest, export |
| **Setup** | `autoimmune-setup` | -- | One-shot: creates collections + seeds knowledge |
| **Milvus** | `milvus-standalone` | 19530 | Shared vector database (external dependency) |

### Port Map

| Port | Service | Protocol |
|---|---|---|
| 8531 | Streamlit UI | HTTP |
| 8532 | FastAPI API | HTTP |
| 19530 | Milvus gRPC | gRPC |
| 9091 | Milvus Proxy | HTTP |

---

## 2. Prerequisites

### Hardware Requirements

**Minimum (Development):**
- 8 CPU cores
- 32 GB RAM (Milvus alone needs ~8 GB for 14 collections)
- 50 GB disk (model cache + collections + data)
- No GPU required (embedding model runs on CPU)

**Recommended (Production / DGX Spark):**
- NVIDIA DGX Spark: GB10 GPU, 128 GB unified LPDDR5x, 20 ARM cores (Grace CPU), NVLink-C2C
- 100 GB SSD for Milvus persistence and model cache
- GPU accelerates embedding generation for large ingestion workloads

### Software Requirements

| Component | Version |
|---|---|
| Python | 3.10+ |
| Docker | 24.0+ |
| Docker Compose | 2.20+ |
| Milvus | 2.4+ (standalone) |

**Python Dependencies** (from `requirements.txt`):

```
pydantic>=2.0
pydantic-settings>=2.7
pymilvus>=2.4.0
sentence-transformers>=2.2.0
anthropic>=0.18.0
streamlit>=1.30.0
fastapi>=0.109.0
uvicorn[standard]>=0.27.0
python-dotenv>=1.0.0
loguru>=0.7.0
numpy>=1.24.0
pandas>=2.0.0
plotly>=5.18.0
reportlab>=4.0.0
PyPDF2>=3.0.0
python-multipart>=0.0.6
prometheus-client>=0.20.0
requests>=2.31.0
httpx>=0.25.0
```

### Network Requirements

- Outbound HTTPS to `api.anthropic.com` (Claude API)
- Outbound HTTPS to `huggingface.co` (first-run model download for `BAAI/bge-small-en-v1.5`)
- Internal TCP to Milvus on port 19530
- Docker network `hcls-network` (external, shared across AI Factory services)

### API Keys

| Key | Required | Source |
|---|---|---|
| `ANTHROPIC_API_KEY` | Yes (for LLM features) | [console.anthropic.com](https://console.anthropic.com) |
| `AUTO_API_KEY` | No (optional auth) | Self-generated |

Without `ANTHROPIC_API_KEY`, the agent starts in degraded mode: vector search works but LLM synthesis is unavailable.

---

## 3. Local Development Setup

### 3.1 Clone and Configure

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_autoimmune_agent

# Create .env from example
cp .env.example .env

# Edit .env and set your Anthropic API key
nano .env
```

Minimum `.env` contents:

```bash
ANTHROPIC_API_KEY=sk-ant-api03-your-key-here
AUTO_MILVUS_HOST=localhost
AUTO_MILVUS_PORT=19530
```

### 3.2 Install Dependencies

The `run.sh` script auto-creates a virtual environment, but you can do it manually:

```bash
python3 -m venv venv
./venv/bin/pip install --upgrade pip
./venv/bin/pip install -r requirements.txt
```

The first run will download the `BAAI/bge-small-en-v1.5` embedding model (~134 MB) from Hugging Face. This is cached at `~/.cache/huggingface/` for subsequent runs.

### 3.3 Start Milvus

If Milvus is not already running from the main AI Factory stack:

```bash
# Standalone Milvus via Docker
docker run -d \
  --name milvus-standalone \
  -p 19530:19530 \
  -p 9091:9091 \
  -v milvus_data:/var/lib/milvus \
  milvusdb/milvus:v2.4.0 \
  milvus run standalone

# Verify Milvus is ready
curl -s http://localhost:9091/healthz
# Expected: {"status":"ok"}
```

### 3.4 Initialize Collections

Create all 14 vector collections and seed the knowledge base:

```bash
./run.sh --setup
```

This runs `scripts/setup_collections.py --seed`, which:
1. Connects to Milvus at `AUTO_MILVUS_HOST:AUTO_MILVUS_PORT`
2. Creates 13 autoimmune-specific collections (skips `genomic_evidence` if it does not exist, as it is a shared read-only collection)
3. Seeds HLA associations, autoantibody panels, biologic therapies, disease activity scores, flare patterns, classification criteria, and cross-disease patterns
4. Embeds all seed data using `BAAI/bge-small-en-v1.5` (384-dimensional vectors)

You can also run it directly with explicit host/port:

```bash
./venv/bin/python scripts/setup_collections.py \
    --host localhost \
    --port 19530 \
    --seed
```

To recreate collections from scratch (drops existing data):

```bash
./venv/bin/python scripts/setup_collections.py \
    --drop-existing \
    --seed
```

### 3.5 Load Demo Data

Demo patient data is stored in `demo_data/` and can be loaded via the API once the server is running:

```bash
# Start the API first (see 3.6), then:
curl -X POST http://localhost:8532/ingest/demo-data
```

Or load specific patient scenarios using the scripts:

```bash
./venv/bin/python scripts/patient_sarah.py
./venv/bin/python scripts/patient_maya.py
./venv/bin/python scripts/patient_emma.py
./venv/bin/python scripts/generate_demo_patients.py
```

### 3.6 Start Services

**Streamlit UI only** (default):

```bash
./run.sh
# Starts on port 8531
# Access: http://localhost:8531
```

**FastAPI API only:**

```bash
./run.sh --api
# Starts on port 8532 with 2 uvicorn workers
# Access: http://localhost:8532/docs (Swagger UI)
```

**Both UI and API** (recommended for development):

```bash
./run.sh --both
# UI on port 8531, API on port 8532
# Ctrl+C triggers graceful shutdown of both processes
```

The `--both` mode uses signal trapping (`SIGTERM`, `SIGINT`) and a 5-second graceful shutdown window before force-killing child processes.

---

## 4. Docker Deployment

### 4.1 Dockerfile Architecture

The Dockerfile uses a multi-stage build for minimal image size:

**Stage 1 -- Builder:**
- Base: `python:3.10-slim`
- Creates a virtualenv at `/opt/venv`
- Installs all dependencies from `requirements.txt` with `--no-cache-dir`

**Stage 2 -- Runtime:**
- Base: `python:3.10-slim`
- Installs `tini` as PID 1 init process (proper signal handling)
- Copies virtualenv from builder stage
- Creates non-root user `autouser`
- Copies application code: `config/`, `src/`, `api/`, `app/`, `scripts/`, `data/`, `.streamlit/`
- Sets ownership to `autouser`
- Exposes ports 8531 and 8532
- Default CMD: Streamlit UI

```dockerfile
# Key Dockerfile lines
FROM python:3.10-slim AS builder
# ... builds /opt/venv ...

FROM python:3.10-slim AS runtime
RUN apt-get update && apt-get install -y --no-install-recommends tini
COPY --from=builder /opt/venv /opt/venv
RUN useradd -m -s /bin/bash autouser
USER autouser
ENTRYPOINT ["tini", "--"]
CMD ["streamlit", "run", "app/autoimmune_ui.py", \
     "--server.port=8531", "--server.address=0.0.0.0", "--server.headless=true"]
```

### 4.2 Building the Image

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_autoimmune_agent

# Build with default tag
docker build -t autoimmune-agent:latest .

# Build with version tag
docker build -t autoimmune-agent:1.0.0 .

# Verify image size
docker images autoimmune-agent
```

### 4.3 Docker Compose Configuration

The `docker-compose.yml` defines three services:

```yaml
version: "3.8"

services:
  autoimmune-streamlit:
    build: .
    container_name: autoimmune-streamlit
    ports:
      - "8531:8531"
    environment:
      - AUTO_MILVUS_HOST=milvus-standalone
      - AUTO_MILVUS_PORT=19530
      - AUTO_STREAMLIT_PORT=8531
      - AUTO_API_PORT=8532
      - ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
    volumes:
      - ./demo_data:/app/demo_data:ro
      - ./data:/app/data
    networks:
      - hcls-network
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "python", "-c",
             "import urllib.request; urllib.request.urlopen('http://localhost:8531/_stcore/health')"]
      interval: 30s
      timeout: 5s
      start_period: 60s
      retries: 3

  autoimmune-api:
    build: .
    container_name: autoimmune-api
    command: ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8532", "--workers", "2"]
    ports:
      - "8532:8532"
    environment:
      - AUTO_MILVUS_HOST=milvus-standalone
      - AUTO_MILVUS_PORT=19530
      - AUTO_API_PORT=8532
      - ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
    volumes:
      - ./demo_data:/app/demo_data:ro
      - ./data:/app/data
    networks:
      - hcls-network
    restart: unless-stopped
    healthcheck:
      test: ["CMD", "python", "-c",
             "import urllib.request; urllib.request.urlopen('http://localhost:8532/healthz')"]
      interval: 30s
      timeout: 5s
      start_period: 30s
      retries: 3

  autoimmune-setup:
    build: .
    container_name: autoimmune-setup
    command: ["python", "scripts/setup_collections.py", "--seed"]
    environment:
      - AUTO_MILVUS_HOST=milvus-standalone
      - AUTO_MILVUS_PORT=19530
    networks:
      - hcls-network
    restart: "no"

networks:
  hcls-network:
    external: true
```

Key design decisions:
- **`demo_data` mounted read-only** (`:ro`) -- prevents accidental modification of demo scenarios
- **`data` mounted read-write** -- for cache, reference data, and ingested documents
- **`hcls-network` is external** -- must be created before starting, shared with Milvus and other AI Factory agents
- **`autoimmune-setup` runs once** (`restart: "no"`) -- creates collections and exits
- **Milvus host is `milvus-standalone`** in Docker (not `localhost`)

### 4.4 Starting Services

```bash
# Ensure the shared network exists
docker network create hcls-network 2>/dev/null || true

# Ensure Milvus is running on hcls-network
docker ps | grep milvus-standalone

# Start all services (setup runs first, then exits)
docker compose up -d

# Watch setup complete
docker logs -f autoimmune-setup

# Once setup exits, verify the long-running services
docker compose ps
```

### 4.5 Verifying Deployment

```bash
# Check service identity
curl -s http://localhost:8532/ | python3 -m json.tool
# Expected: {"service": "Precision Autoimmune Intelligence Agent", "version": "1.0.0", ...}

# Detailed health check
curl -s http://localhost:8532/health | python3 -m json.tool
# Expected: milvus_connected=true, collections=14, embedder_loaded=true, llm_available=true

# Lightweight probe
curl -s http://localhost:8532/healthz
# Expected: {"status": "ok"}

# Streamlit health
curl -s http://localhost:8531/_stcore/health
# Expected: "ok"

# List collections and vector counts
curl -s http://localhost:8532/collections | python3 -m json.tool

# Test a query
curl -s -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What are the early biomarkers for lupus flare?"}' | python3 -m json.tool

# Check Prometheus metrics
curl -s http://localhost:8532/metrics
```

---

## 5. DGX Spark Deployment

### 5.1 Hardware Specifications

The NVIDIA DGX Spark provides:

| Resource | Specification |
|---|---|
| GPU | GB10 (NVIDIA Blackwell architecture) |
| Memory | 128 GB unified LPDDR5x (shared CPU/GPU) |
| CPU | 20 ARM cores (NVIDIA Grace) |
| Interconnect | NVLink-C2C |
| Storage | NVMe SSD |
| Price | $3,999 |

The unified memory architecture eliminates PCIe bottlenecks for embedding model inference and allows Milvus to leverage the full 128 GB for large collection indexes.

### 5.2 DGX Spark-Specific Configuration

The autoimmune agent integrates with the main AI Factory stack via `docker-compose.dgx-spark.yml`. In this deployment model, the agent runs as a single container (port 8000 internally, mapped to 8105/8106 externally) rather than the standalone 3-service topology.

From the main compose file:

```yaml
precision-autoimmune-agent:
  build:
    context: ./ai_agent_adds/precision_autoimmune_agent
    dockerfile: Dockerfile
  restart: unless-stopped
  ports:
    - "8105:8531"
    - "8106:8000"
  environment:
    <<: *common-env
    AGENT_NAME: precision-autoimmune-agent
  depends_on:
    milvus:
      condition: service_healthy
```

The `*common-env` anchor provides:

```yaml
MILVUS_HOST: milvus
MILVUS_PORT: "19530"
ANTHROPIC_API_KEY: ${ANTHROPIC_API_KEY}
EMBEDDING_MODEL: BAAI/bge-small-en-v1.5
LOG_LEVEL: INFO
```

### 5.3 GPU Memory Considerations

On the DGX Spark, all services share 128 GB of unified memory:

| Component | Estimated Memory |
|---|---|
| Milvus (14 collections, IVF_FLAT indexes) | 4--8 GB |
| BGE-small-en-v1.5 embedding model | ~500 MB |
| Streamlit UI process | ~200 MB |
| FastAPI + 2 uvicorn workers | ~600 MB |
| Sentence-transformers runtime | ~1 GB |
| **Total autoimmune agent** | **~6--10 GB** |

This leaves ample headroom for the other 4 intelligence agents, the drug discovery pipeline, and Milvus itself.

### 5.4 Integration with Main docker-compose.dgx-spark.yml

To deploy as part of the full AI Factory stack:

```bash
cd /home/adam/projects/hcls-ai-factory

# Set API key
export ANTHROPIC_API_KEY=sk-ant-api03-your-key-here

# Start the full stack (includes Milvus, all 5 agents, monitoring, landing page)
docker compose -f docker-compose.dgx-spark.yml up -d

# Verify the autoimmune agent is healthy
curl -s http://localhost:8106/health | python3 -m json.tool
```

The landing page at `http://localhost:8080` provides a unified health dashboard showing all agents, including the autoimmune agent's status.

---

## 6. Environment Variables

All agent-specific variables use the `AUTO_` prefix and are managed by Pydantic Settings (`config/settings.py`). They can be set via environment variables, `.env` file, or Docker Compose `environment` blocks.

### Connection

| Variable | Default | Description |
|---|---|---|
| `ANTHROPIC_API_KEY` | *(empty)* | Anthropic API key for Claude. **Required for LLM features.** No `AUTO_` prefix. |
| `AUTO_MILVUS_HOST` | `localhost` | Milvus server hostname. Use `milvus-standalone` in Docker. |
| `AUTO_MILVUS_PORT` | `19530` | Milvus gRPC port. |
| `EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | Sentence-transformers model name. No `AUTO_` prefix. |

### Ports

| Variable | Default | Description |
|---|---|---|
| `AUTO_STREAMLIT_PORT` | `8531` | Streamlit UI listen port. |
| `AUTO_API_PORT` | `8532` | FastAPI listen port. |

### LLM

| Variable | Default | Description |
|---|---|---|
| `AUTO_LLM_MODEL` | `claude-sonnet-4-6` | Claude model ID for synthesis. |
| `AUTO_LLM_MAX_TOKENS` | `4096` | Maximum tokens in LLM response. |
| `AUTO_LLM_TEMPERATURE` | `0.2` | LLM temperature (lower = more deterministic). |
| `AUTO_LLM_MAX_RETRIES` | `3` | Retry count for failed LLM calls. |

### RAG Parameters

| Variable | Default | Description |
|---|---|---|
| `AUTO_TOP_K_PER_COLLECTION` | `5` | Number of top results retrieved per collection. |
| `AUTO_SCORE_THRESHOLD` | `0.40` | Minimum cosine similarity score to include a result. |
| `AUTO_MAX_EVIDENCE_ITEMS` | `30` | Maximum total evidence items across all collections. |
| `AUTO_CONVERSATION_MEMORY_SIZE` | `3` | Number of prior conversation turns retained for context. |

### Security

| Variable | Default | Description |
|---|---|---|
| `AUTO_API_KEY` | *(empty)* | API key for endpoint authentication. Empty = no auth required. |
| `AUTO_CORS_ORIGINS` | `http://localhost:8080,http://localhost:8531` | Comma-separated allowed CORS origins. |
| `AUTO_MAX_REQUEST_SIZE_MB` | `50` | Maximum upload size for PDF ingestion. |

### Thresholds

| Variable | Default | Description |
|---|---|---|
| `AUTO_CITATION_HIGH` | `0.80` | Cosine score threshold for high-confidence citations. |
| `AUTO_CITATION_MEDIUM` | `0.60` | Cosine score threshold for medium-confidence citations. |
| `AUTO_FLARE_RISK_IMMINENT` | `0.8` | Flare risk score threshold: imminent. |
| `AUTO_FLARE_RISK_HIGH` | `0.6` | Flare risk score threshold: high. |
| `AUTO_FLARE_RISK_MODERATE` | `0.4` | Flare risk score threshold: moderate. |

### Timeouts

| Variable | Default | Description |
|---|---|---|
| `AUTO_REQUEST_TIMEOUT_SECONDS` | `60` | Overall request timeout. |
| `AUTO_MILVUS_TIMEOUT_SECONDS` | `10` | Milvus connection/query timeout. |

### Streaming and Metrics

| Variable | Default | Description |
|---|---|---|
| `AUTO_STREAMING_ENABLED` | `True` | Enable SSE streaming for `/query/stream` endpoint. |
| `AUTO_METRICS_ENABLED` | `True` | Enable Prometheus metrics at `/metrics`. |

### Logging

| Variable | Default | Description |
|---|---|---|
| `AUTO_LOG_LEVEL` | `INFO` | Log level: DEBUG, INFO, WARNING, ERROR. |
| `AUTO_LOG_DIR` | `<project_root>/logs` | Directory for log files. |

### Document Processing

| Variable | Default | Description |
|---|---|---|
| `AUTO_MAX_CHUNK_SIZE` | `2500` | Maximum characters per text chunk when processing PDFs. |
| `AUTO_CHUNK_OVERLAP` | `200` | Character overlap between consecutive chunks. |
| `AUTO_PDF_DPI` | `200` | DPI for PDF rendering during OCR. |

---

## 7. Milvus Configuration

### 7.1 Collection Initialization

The agent manages 14 vector collections, each with a domain-specific schema:

| # | Collection Name | Description |
|---|---|---|
| 1 | `autoimmune_clinical_documents` | Ingested patient records (PDFs) |
| 2 | `autoimmune_patient_labs` | Lab results with flag analysis |
| 3 | `autoimmune_autoantibody_panels` | Autoantibody test result panels |
| 4 | `autoimmune_hla_associations` | HLA allele to disease risk mapping |
| 5 | `autoimmune_disease_criteria` | ACR/EULAR classification criteria |
| 6 | `autoimmune_disease_activity` | Activity scoring (DAS28, SLEDAI, BASDAI, etc.) |
| 7 | `autoimmune_flare_patterns` | Flare prediction biomarker patterns |
| 8 | `autoimmune_biologic_therapies` | Biologic drug database with PGx |
| 9 | `autoimmune_pgx_rules` | Pharmacogenomic dosing rules |
| 10 | `autoimmune_clinical_trials` | Autoimmune clinical trials |
| 11 | `autoimmune_literature` | Published literature and research |
| 12 | `autoimmune_patient_timelines` | Patient diagnostic timeline events |
| 13 | `autoimmune_cross_disease` | Cross-disease overlap syndromes |
| 14 | `genomic_evidence` | **Shared read-only** (from genomics pipeline) |

The `genomic_evidence` collection is shared across all 5 intelligence agents and is never created or dropped by the autoimmune agent. It is accessed read-only if it already exists.

### 7.2 Index Parameters

All collections use the same index and search configuration:

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

| Parameter | Value | Rationale |
|---|---|---|
| `metric_type` | `COSINE` | Standard for normalized text embeddings from BGE models |
| `index_type` | `IVF_FLAT` | Good recall/speed balance for collections under 1M vectors |
| `nlist` | `1024` | Number of Voronoi cells for IVF partitioning |
| `nprobe` | `16` | Cells searched at query time (higher = better recall, slower) |
| `embedding_dim` | `384` | Output dimension of `BAAI/bge-small-en-v1.5` |

### 7.3 Connection Settings

The collection manager uses a named Milvus connection alias (`autoimmune_agent`) with automatic reconnection:

```python
connections.connect(alias="autoimmune_agent", host=settings.MILVUS_HOST, port=settings.MILVUS_PORT)
```

The API startup includes a 2-attempt retry loop with a 2-second delay between attempts. If Milvus is unavailable, the agent starts in degraded mode with vector search disabled.

### 7.4 Shared Collections

The `genomic_evidence` collection is populated by the genomics pipeline (stage 1 of the AI Factory). The autoimmune agent reads from it with a low weight (0.02) in the multi-collection RAG search. It is weighted low because genomic variants are supplementary to the autoimmune-specific clinical data.

Collection search weights (sum to 1.0):

| Collection | Weight |
|---|---|
| `clinical_documents` | 0.18 |
| `patient_labs` | 0.14 |
| `autoantibody_panels` | 0.12 |
| `hla_associations` | 0.08 |
| `disease_criteria` | 0.08 |
| `disease_activity` | 0.07 |
| `flare_patterns` | 0.06 |
| `biologic_therapies` | 0.06 |
| `clinical_trials` | 0.05 |
| `literature` | 0.05 |
| `pgx_rules` | 0.04 |
| `patient_timelines` | 0.03 |
| `cross_disease` | 0.02 |
| `genomic_evidence` | 0.02 |

---

## 8. Health Checks and Monitoring

### 8.1 Health Endpoints

**`GET /`** -- Service identity:

```json
{
  "service": "Precision Autoimmune Intelligence Agent",
  "version": "1.0.0",
  "status": "running",
  "ports": {"api": 8532, "ui": 8531}
}
```

**`GET /health`** -- Detailed health (checks Milvus, collections, embedder, LLM):

```json
{
  "status": "healthy",
  "service": "autoimmune-agent",
  "milvus_connected": true,
  "collections": 14,
  "total_vectors": 1247,
  "embedder_loaded": true,
  "llm_available": true,
  "uptime_seconds": 3600
}
```

**`GET /healthz`** -- Lightweight probe (for orchestrators and landing page):

```json
{"status": "ok"}
```

**Streamlit health** (`GET /_stcore/health`): Returns `"ok"` as plain text.

**`GET /metrics`** -- Prometheus-compatible metrics (text format):

```
# HELP autoimmune_agent_up Whether the agent is running
# TYPE autoimmune_agent_up gauge
autoimmune_agent_up 1
# HELP autoimmune_collection_vectors Number of vectors per collection
# TYPE autoimmune_collection_vectors gauge
autoimmune_collection_vectors{collection="autoimmune_clinical_documents"} 342
autoimmune_collection_vectors{collection="autoimmune_hla_associations"} 87
...
# HELP autoimmune_agent_uptime_seconds Agent uptime
# TYPE autoimmune_agent_uptime_seconds gauge
autoimmune_agent_uptime_seconds 3600
```

### 8.2 Docker Health Checks

Both long-running containers include health checks:

**Streamlit container:**

```yaml
healthcheck:
  test: ["CMD", "python", "-c",
         "import urllib.request; urllib.request.urlopen('http://localhost:8531/_stcore/health')"]
  interval: 30s
  timeout: 5s
  start_period: 60s    # Streamlit needs more startup time
  retries: 3
```

**API container:**

```yaml
healthcheck:
  test: ["CMD", "python", "-c",
         "import urllib.request; urllib.request.urlopen('http://localhost:8532/healthz')"]
  interval: 30s
  timeout: 5s
  start_period: 30s
  retries: 3
```

### 8.3 Prometheus Metrics

The `/metrics` endpoint exposes:

| Metric | Type | Description |
|---|---|---|
| `autoimmune_agent_up` | Gauge | 1 if agent is running |
| `autoimmune_collection_vectors` | Gauge | Vector count per collection (labeled) |
| `autoimmune_agent_uptime_seconds` | Gauge | Agent uptime in seconds |

Prometheus scrape config (add to `prometheus.yml`):

```yaml
scrape_configs:
  - job_name: 'autoimmune-agent'
    scrape_interval: 30s
    static_configs:
      - targets: ['autoimmune-api:8532']
    metrics_path: /metrics
```

### 8.4 Grafana Dashboard Integration

The AI Factory monitoring stack includes Grafana at port 3000. To add the autoimmune agent:

1. Ensure Prometheus is scraping the `/metrics` endpoint (see above).
2. Import or create a dashboard with panels for:
   - Agent uptime
   - Collection vector counts (bar chart by collection)
   - Health status (using the `/health` endpoint via Grafana HTTP datasource)
3. Set alert rules for:
   - `autoimmune_agent_up == 0` (agent down)
   - `autoimmune_collection_vectors{collection="autoimmune_clinical_documents"} == 0` (empty clinical collection)

---

## 9. Security Configuration

### 9.1 API Key Authentication

When `AUTO_API_KEY` is set, all endpoints except `/`, `/health`, `/healthz`, and `/metrics` require authentication via:

- **Header:** `X-API-Key: <your-key>`
- **Query parameter:** `?api_key=<your-key>`

```bash
# Set the API key
export AUTO_API_KEY="your-secure-api-key-here"

# Authenticated request
curl -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-secure-api-key-here" \
  -d '{"question": "What biologics are indicated for lupus?"}'

# Without key (returns 401)
curl -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What biologics are indicated for lupus?"}'
# {"detail": "Invalid or missing API key"}
```

### 9.2 CORS Configuration

CORS origins are controlled by `AUTO_CORS_ORIGINS` (comma-separated):

```bash
# Default: allows landing page and Streamlit
AUTO_CORS_ORIGINS=http://localhost:8080,http://localhost:8531

# Production: restrict to specific domains
AUTO_CORS_ORIGINS=https://aifactory.example.com,https://portal.example.com
```

The CORS middleware allows all methods and headers (`allow_methods=["*"]`, `allow_headers=["*"]`) with credentials enabled.

### 9.3 Non-Root Docker User

The Dockerfile creates and switches to a non-root user:

```dockerfile
RUN useradd -m -s /bin/bash autouser
# ... copy application files ...
RUN chown -R autouser:autouser /app
USER autouser
```

This ensures the container process cannot modify system files or escalate privileges. The application runs entirely within `/app` owned by `autouser`.

### 9.4 Secret Management

**Current approach:** API keys are passed via environment variables or `.env` file.

**Production recommendations:**

- Never commit `.env` files to version control. The `.env.example` template uses placeholder values.
- Use Docker secrets or a secrets manager (HashiCorp Vault, AWS Secrets Manager) for `ANTHROPIC_API_KEY`.
- Rotate `AUTO_API_KEY` periodically.
- The `ANTHROPIC_API_KEY` is validated at startup; if missing, the agent logs a warning and enters degraded mode rather than failing.

**Request size limiting:** The `AUTO_MAX_REQUEST_SIZE_MB` setting (default: 50 MB) rejects uploads exceeding the limit with HTTP 413.

---

## 10. Production Hardening

### 10.1 Resource Limits

Add resource constraints to the Docker Compose services:

```yaml
autoimmune-api:
  # ... existing config ...
  deploy:
    resources:
      limits:
        memory: 4G
        cpus: "4.0"
      reservations:
        memory: 2G
        cpus: "1.0"

autoimmune-streamlit:
  deploy:
    resources:
      limits:
        memory: 2G
        cpus: "2.0"
      reservations:
        memory: 1G
        cpus: "0.5"
```

### 10.2 Restart Policies

The Docker Compose file already uses `unless-stopped` for long-running services and `"no"` for the one-shot setup container. For production, consider:

```yaml
autoimmune-api:
  restart: unless-stopped
  # Alternatively, for Kubernetes-like behavior:
  # restart: on-failure
  # deploy:
  #   restart_policy:
  #     condition: on-failure
  #     delay: 5s
  #     max_attempts: 10
  #     window: 120s
```

### 10.3 Log Management

**Loguru configuration** (from `config/logging.py`):

| Setting | Value |
|---|---|
| Console output | `stderr` with color |
| File output | `logs/autoimmune-agent.log` (API uses `autoimmune-api.log`) |
| Rotation | 10 MB per file |
| Retention | 5 rotated files |
| Thread safety | `enqueue=True` |
| Encoding | UTF-8 |

**Docker log driver** -- add to compose for centralized logging:

```yaml
autoimmune-api:
  logging:
    driver: "json-file"
    options:
      max-size: "10m"
      max-file: "5"
```

**Log level tuning** -- reduce verbosity in production:

```bash
AUTO_LOG_LEVEL=WARNING
```

### 10.4 Backup Strategy

**Milvus data:**

```bash
# Milvus data is stored in the milvus_data Docker volume
docker run --rm -v milvus_data:/data -v $(pwd)/backups:/backup \
  alpine tar czf /backup/milvus-$(date +%Y%m%d).tar.gz -C /data .
```

**Application data:**

```bash
# Back up ingested documents and cache
tar czf backups/autoimmune-data-$(date +%Y%m%d).tar.gz \
  data/ demo_data/ logs/
```

**Collection recreation:** Collections can be fully rebuilt from seed data and demo patients:

```bash
./run.sh --setup
curl -X POST http://localhost:8532/ingest/demo-data
```

### 10.5 Performance Tuning

**Embedding batch size** -- increase for bulk ingestion:

```bash
AUTO_EMBEDDING_BATCH_SIZE=64  # default: 32
```

**Milvus search parallelism** -- the `search_all` method uses `ThreadPoolExecutor` with `max_workers=6` by default. For higher throughput:

```python
# In rag_engine.py or via settings
results = cm.search_all(query_embedding, max_workers=10)
```

**Uvicorn workers** -- the API defaults to 2 workers. Increase for higher concurrency:

```yaml
autoimmune-api:
  command: ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8532", "--workers", "4"]
```

**Top-K tuning** -- reduce `AUTO_TOP_K_PER_COLLECTION` to 3 for faster responses at the cost of recall, or increase to 10 for maximum evidence coverage.

**IVF_FLAT nprobe** -- the current `nprobe=16` balances speed and recall. For latency-critical deployments, reduce to 8. For maximum recall, increase to 32 or 64.

---

## 11. Scaling and High Availability

### 11.1 Uvicorn Workers

The API runs 2 uvicorn workers by default. Each worker loads its own copy of the embedding model and Milvus connection. On the DGX Spark (20 cores), you can safely run 4--6 workers:

```bash
# Docker override
docker compose exec autoimmune-api \
  uvicorn api.main:app --host 0.0.0.0 --port 8532 --workers 4

# Or modify docker-compose.yml command
command: ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8532", "--workers", "4"]
```

Memory impact: each additional worker adds approximately 500 MB (embedding model + Milvus client + application state).

### 11.2 Milvus Scaling

For larger deployments (>1M vectors per collection), consider:

1. **Milvus Cluster mode** -- replaces standalone with distributed components (query node, data node, index node).
2. **HNSW index** -- replace IVF_FLAT with HNSW for better recall at high vector counts:

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "HNSW",
    "params": {"M": 16, "efConstruction": 256},
}
SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"ef": 64},
}
```

3. **Collection partitioning** -- partition `autoimmune_clinical_documents` by `patient_id` for faster patient-scoped queries.

### 11.3 Load Balancing Considerations

For multi-instance deployments:

- **Stateless API** -- the FastAPI server is stateless (no session affinity required). Place behind a reverse proxy (NGINX, Traefik) or cloud load balancer.
- **Milvus connection pooling** -- each API instance maintains its own Milvus connection via the `autoimmune_agent` alias. No connection pool sharing is needed.
- **Streamlit** -- Streamlit maintains WebSocket connections for UI interactivity. If scaling horizontally, use sticky sessions or deploy one Streamlit instance per user group.

```nginx
# NGINX example
upstream autoimmune_api {
    server autoimmune-api-1:8532;
    server autoimmune-api-2:8532;
}

server {
    listen 8532;
    location / {
        proxy_pass http://autoimmune_api;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

---

## 12. Troubleshooting

### 12.1 Milvus Connection Refused

**Symptom:** `MilvusException: connection refused` at startup.

**Solution:**
```bash
# Verify Milvus is running
docker ps | grep milvus
curl -s http://localhost:9091/healthz

# Check AUTO_MILVUS_HOST is correct
# Local: localhost | Docker: milvus-standalone
echo $AUTO_MILVUS_HOST

# Verify network connectivity (from inside container)
docker exec autoimmune-api python -c "
from pymilvus import connections
connections.connect(host='milvus-standalone', port=19530)
print('Connected')
"
```

### 12.2 Embedding Model Download Fails

**Symptom:** `OSError: Can't load tokenizer for 'BAAI/bge-small-en-v1.5'`.

**Solution:**
```bash
# Pre-download the model
python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('BAAI/bge-small-en-v1.5')"

# Or set a custom cache directory
export SENTENCE_TRANSFORMERS_HOME=/app/data/models
export HF_HOME=/app/data/huggingface
```

### 12.3 ANTHROPIC_API_KEY Not Detected

**Symptom:** `LLM features will be unavailable (demo mode)` in logs.

**Solution:**
```bash
# Verify the key is set
echo $ANTHROPIC_API_KEY | head -c 20

# In Docker, check it's passed through
docker exec autoimmune-api env | grep ANTHROPIC

# Common issue: .env file not loaded in Docker
# Fix: pass explicitly in docker-compose.yml environment block
```

### 12.4 Collections Not Created

**Symptom:** `/health` shows `"collections": 0`.

**Solution:**
```bash
# Run setup manually
docker exec autoimmune-api python scripts/setup_collections.py --seed

# Or check if setup container ran
docker logs autoimmune-setup
```

### 12.5 Streamlit Not Loading

**Symptom:** Browser shows connection refused on port 8531.

**Solution:**
```bash
# Check if Streamlit is actually running
docker logs autoimmune-streamlit

# Common issue: port conflict
lsof -i :8531

# Verify health
curl -s http://localhost:8531/_stcore/health
```

### 12.6 PDF Upload Returns 422

**Symptom:** `"No text extracted from PDF"` error.

**Solution:**
- Ensure the PDF contains selectable text (not scanned images without OCR).
- Check `AUTO_MAX_REQUEST_SIZE_MB` is large enough for the file.
- Verify the file has a `.pdf` extension (other formats return 400).

### 12.7 Slow Query Response

**Symptom:** Queries take >10 seconds.

**Solution:**
```bash
# Reduce collections searched
curl -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -d '{"question": "...", "collections_filter": ["autoimmune_clinical_documents", "autoimmune_patient_labs"]}'

# Reduce top_k
AUTO_TOP_K_PER_COLLECTION=3

# Check Milvus index status
docker exec autoimmune-api python -c "
from pymilvus import Collection, connections
connections.connect(host='milvus-standalone', port=19530)
c = Collection('autoimmune_clinical_documents')
print(c.indexes)
"
```

### 12.8 Out of Memory

**Symptom:** Container killed by OOM.

**Solution:**
```yaml
# Increase memory limit
deploy:
  resources:
    limits:
      memory: 6G

# Or reduce workers
command: ["uvicorn", "api.main:app", "--host", "0.0.0.0", "--port", "8532", "--workers", "1"]
```

### 12.9 Docker Network Not Found

**Symptom:** `network hcls-network declared as external, but could not be found`.

**Solution:**
```bash
docker network create hcls-network
```

### 12.10 Permission Denied in Container

**Symptom:** `PermissionError: [Errno 13]` writing to `/app/data` or `/app/logs`.

**Solution:**
```bash
# Fix host-side permissions for mounted volumes
chmod -R 777 data/ logs/

# Or match the autouser UID
docker exec autoimmune-api id autouser
# Then: chown -R <uid>:<gid> data/ logs/
```

### 12.11 Collection Weight Mismatch Warning

**Symptom:** `Collection weights sum to X.XXX, expected ~1.0` in logs.

**Solution:** This is a warning, not an error. If you override individual `AUTO_WEIGHT_*` variables, ensure all 14 weights sum to approximately 1.0. The tolerance is 0.05.

### 12.12 CORS Errors in Browser

**Symptom:** `Access-Control-Allow-Origin` errors in browser console.

**Solution:**
```bash
# Add your frontend origin
AUTO_CORS_ORIGINS=http://localhost:8080,http://localhost:8531,http://your-frontend:3000
```

---

## 13. Quick Reference

### Command Cheat Sheet

```bash
# ── Local Development ──
./run.sh              # Start Streamlit UI (port 8531)
./run.sh --api        # Start FastAPI (port 8532)
./run.sh --both       # Start both UI and API
./run.sh --setup      # Create collections + seed knowledge

# ── Docker ──
docker build -t autoimmune-agent:latest .
docker network create hcls-network
docker compose up -d
docker compose down
docker compose logs -f autoimmune-api
docker compose ps

# ── Collection Management ──
python scripts/setup_collections.py --seed
python scripts/setup_collections.py --drop-existing --seed
python scripts/setup_collections.py --host milvus-standalone --port 19530

# ── Health Checks ──
curl -s http://localhost:8532/
curl -s http://localhost:8532/health
curl -s http://localhost:8532/healthz
curl -s http://localhost:8532/metrics
curl -s http://localhost:8531/_stcore/health

# ── API Queries ──
curl -s -X POST http://localhost:8532/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What are early warning signs of a lupus flare?"}'

curl -s -X POST http://localhost:8532/search \
  -H "Content-Type: application/json" \
  -d '{"question": "rituximab pharmacogenomics"}'

curl -s http://localhost:8532/collections

# ── Full Stack (DGX Spark) ──
docker compose -f docker-compose.dgx-spark.yml up -d
docker compose -f docker-compose.dgx-spark.yml down
```

### Port Map

| Port | Service | Health Check |
|---|---|---|
| 8531 | Streamlit UI | `/_stcore/health` |
| 8532 | FastAPI API | `/healthz` |
| 19530 | Milvus gRPC | -- |
| 9091 | Milvus Proxy | `/healthz` |

### Health Check URLs

| URL | Purpose | Expected Response |
|---|---|---|
| `http://localhost:8532/` | Service identity | `{"service": "Precision Autoimmune Intelligence Agent", ...}` |
| `http://localhost:8532/health` | Full health | `{"status": "healthy", "milvus_connected": true, ...}` |
| `http://localhost:8532/healthz` | Lightweight probe | `{"status": "ok"}` |
| `http://localhost:8532/metrics` | Prometheus metrics | Text-format metrics |
| `http://localhost:8531/_stcore/health` | Streamlit probe | `ok` |

### API Endpoints

| Method | Path | Description |
|---|---|---|
| `GET` | `/` | Service identity |
| `GET` | `/health` | Detailed health check |
| `GET` | `/healthz` | Lightweight probe |
| `GET` | `/metrics` | Prometheus metrics |
| `POST` | `/query` | RAG query (retrieve + synthesize) |
| `POST` | `/query/stream` | Streaming RAG query (SSE) |
| `POST` | `/search` | Evidence-only search (no LLM) |
| `POST` | `/analyze` | Full patient analysis pipeline |
| `POST` | `/differential` | Differential diagnosis |
| `POST` | `/ingest/upload` | Upload and ingest a PDF |
| `POST` | `/ingest/demo-data` | Ingest all demo patients |
| `GET` | `/collections` | List collections with stats |
| `POST` | `/collections/create` | Create/recreate collections |
| `POST` | `/export` | Export report (markdown/FHIR/PDF) |

### Key File Paths

```
precision_autoimmune_agent/
  config/
    settings.py          # All AUTO_* configuration (Pydantic Settings)
    logging.py           # Loguru dual-sink setup
  src/
    collections.py       # 14 collection schemas + Milvus manager
    rag_engine.py        # Multi-collection RAG with weighted search
    agent.py             # Core autoimmune analysis agent
    diagnostic_engine.py # Differential diagnosis + criteria scoring
    knowledge.py         # Static knowledge base (HLA, antibodies, drugs)
    models.py            # Pydantic data models
    document_processor.py # PDF ingestion pipeline
    timeline_builder.py  # Patient diagnostic odyssey timeline
    export.py            # Report export (Markdown, FHIR R4, PDF)
  api/
    main.py              # FastAPI application with all endpoints
  app/
    autoimmune_ui.py     # Streamlit interactive UI
  scripts/
    setup_collections.py # Collection creation + knowledge seeding
  run.sh                 # Multi-mode startup script
  Dockerfile             # Multi-stage build
  docker-compose.yml     # 3-service compose
  requirements.txt       # Python dependencies
  .env.example           # Environment variable template
  logs/                  # Loguru output (autoimmune-agent.log)
  data/                  # Cache, reference data, ingested docs
  demo_data/             # Demo patient scenarios (read-only in Docker)
```
