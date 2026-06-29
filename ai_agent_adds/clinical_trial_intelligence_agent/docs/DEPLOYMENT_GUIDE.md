# Clinical Trial Intelligence Agent -- Deployment Guide

**Version:** 2.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Integrated Deployment (HCLS AI Factory)](#2-integrated-deployment)
3. [Standalone Docker Deployment](#3-standalone-docker-deployment)
4. [Manual Development Setup](#4-manual-development-setup)
5. [Milvus Configuration](#5-milvus-configuration)
6. [Collection Setup and Seeding](#6-collection-setup-and-seeding)
7. [Data Ingest](#7-data-ingest)
8. [Security Configuration](#8-security-configuration)
9. [Environment Variables](#9-environment-variables)
10. [Health Verification](#10-health-verification)
11. [Monitoring Setup](#11-monitoring-setup)
12. [Performance Tuning](#12-performance-tuning)
13. [Troubleshooting](#13-troubleshooting)
14. [Backup and Recovery](#14-backup-and-recovery)

---

## 1. Prerequisites

### Hardware Requirements

| Component | Minimum | Recommended |
|---|---|---|
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Storage | 20 GB | 100+ GB (with Milvus data) |
| GPU | Not required | NVIDIA GPU (for embedding acceleration) |

### Software Requirements

| Software | Version | Purpose |
|---|---|---|
| Python | 3.10+ | Runtime |
| Docker | 24.0+ | Container deployment |
| Docker Compose | 2.20+ | Multi-container orchestration |
| pip | 23.0+ | Package management |

### Network Requirements

| Port | Service | Direction |
|---|---|---|
| 8538 | FastAPI API | Inbound |
| 8128 | Streamlit UI | Inbound |
| 19530 | Milvus | Internal |
| 2379 | etcd | Internal |
| 9000 | MinIO | Internal |

---

## 2. Integrated Deployment

The Clinical Trial Intelligence Agent is deployed as part of the HCLS AI Factory stack using the master `docker-compose.dgx-spark.yml`:

### 2.1 Start the Full Stack

```bash
cd /home/adam/projects/hcls-ai-factory
./start-factory.sh
```

This starts all services including:
- Milvus standalone + etcd + MinIO
- Clinical Trial Intelligence Agent (API + UI)
- All peer intelligence agents
- Prometheus + Grafana monitoring
- Landing page

### 2.2 Start Only the Clinical Trial Agent

```bash
cd /home/adam/projects/hcls-ai-factory
docker-compose -f docker-compose.dgx-spark.yml up -d \
    milvus-standalone etcd minio \
    clinical-trial-agent-api clinical-trial-agent-ui
```

### 2.3 Verify Deployment

```bash
# API health
curl http://localhost:8538/health

# UI access
open http://localhost:8128

# Collection status
curl http://localhost:8538/collections
```

---

## 3. Standalone Docker Deployment

### 3.1 Build the Docker Image

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/clinical_trial_intelligence_agent

# Build
docker build -t clinical-trial-agent:latest .
```

### 3.2 Docker Compose (Standalone)

Create or use the existing `docker-compose.yml`:

```yaml
version: '3.8'

services:
  etcd:
    image: quay.io/coreos/etcd:v3.5.5
    environment:
      - ETCD_AUTO_COMPACTION_MODE=revision
      - ETCD_AUTO_COMPACTION_RETENTION=1000
    ports:
      - "2379:2379"

  minio:
    image: minio/minio:RELEASE.2023-03-20T20-16-18Z
    environment:
      MINIO_ACCESS_KEY: minioadmin
      MINIO_SECRET_KEY: minioadmin
    ports:
      - "9000:9000"
    command: minio server /data

  milvus:
    image: milvusdb/milvus:v2.3.3
    ports:
      - "19530:19530"
    depends_on:
      - etcd
      - minio
    environment:
      ETCD_ENDPOINTS: etcd:2379
      MINIO_ADDRESS: minio:9000

  trial-api:
    build: .
    ports:
      - "8538:8538"
    environment:
      - TRIAL_MILVUS_HOST=milvus
      - TRIAL_MILVUS_PORT=19530
      - TRIAL_ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}
      - TRIAL_API_KEY=${TRIAL_API_KEY:-}
    depends_on:
      - milvus
    command: uvicorn api.main:app --host 0.0.0.0 --port 8538

  trial-ui:
    build: .
    ports:
      - "8128:8128"
    environment:
      - TRIAL_API_BASE=http://trial-api:8538
    depends_on:
      - trial-api
    command: streamlit run app/trial_ui.py --server.port 8128 --server.address 0.0.0.0
```

### 3.3 Launch

```bash
# Set API key
export ANTHROPIC_API_KEY="sk-ant-..."

# Start
docker-compose up -d

# Check logs
docker-compose logs -f trial-api

# Setup collections
docker-compose exec trial-api python scripts/setup_collections.py

# Seed knowledge
docker-compose exec trial-api python scripts/seed_knowledge.py
```

---

## 4. Manual Development Setup

### 4.1 Virtual Environment

```bash
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/clinical_trial_intelligence_agent

# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 4.2 Environment Configuration

```bash
# Copy environment template
cat > .env << 'EOF'
TRIAL_MILVUS_HOST=localhost
TRIAL_MILVUS_PORT=19530
TRIAL_ANTHROPIC_API_KEY=sk-ant-...
TRIAL_API_PORT=8538
TRIAL_STREAMLIT_PORT=8128
TRIAL_API_KEY=
TRIAL_INGEST_ENABLED=false
EOF
```

### 4.3 Start Milvus (if not running)

```bash
# Using Docker
docker run -d --name milvus-standalone \
    -p 19530:19530 \
    -p 9091:9091 \
    milvusdb/milvus:v2.3.3 \
    milvus run standalone
```

### 4.4 Setup Collections and Seed Data

```bash
python scripts/setup_collections.py
python scripts/seed_knowledge.py
```

### 4.5 Start Services

```bash
# Terminal 1: API server
uvicorn api.main:app --host 0.0.0.0 --port 8538 --reload

# Terminal 2: Streamlit UI
streamlit run app/trial_ui.py --server.port 8128
```

### 4.6 Run Tests

```bash
python -m pytest tests/ -v
# Expected: 769 passed in 0.47s
```

---

## 5. Milvus Configuration

### 5.1 Collection Specifications

All 14 collections use identical index parameters:

| Parameter | Value |
|---|---|
| Embedding Dimension | 384 (BGE-small-en-v1.5) |
| Index Type | IVF_FLAT |
| Metric Type | COSINE |
| nlist | 128 |

### 5.2 Milvus Tuning for Production

```yaml
# milvus.yaml overrides for production
queryCoord:
  autoBalance: true
dataCoord:
  compaction:
    enabled: true
    autoCompaction: true
proxy:
  maxTaskNum: 1024
  maxNameLength: 255
```

### 5.3 Memory Allocation

For 14 collections with ~250K total records:

| Metric | Estimate |
|---|---|
| Raw vectors | ~250K * 384 * 4 bytes = ~384 MB |
| Index overhead | ~200 MB (IVF_FLAT) |
| Metadata | ~500 MB |
| **Total Milvus RAM** | **~2 GB** |

### 5.4 Scaling to Production

For larger deployments (millions of records):

1. Switch from Milvus Standalone to Milvus Cluster
2. Consider IVF_SQ8 or HNSW index for better search performance
3. Increase nlist to 256 or 512
4. Configure MinIO with multi-disk storage
5. Enable replica loading for read-heavy workloads

---

## 6. Collection Setup and Seeding

### 6.1 Create Collections

```bash
python scripts/setup_collections.py
```

This creates all 14 collections with proper schemas and indexes:

```
trial_protocols, trial_eligibility, trial_endpoints, trial_sites,
trial_investigators, trial_results, trial_regulatory, trial_literature,
trial_biomarkers, trial_safety, trial_rwe, trial_adaptive,
trial_guidelines, genomic_evidence
```

### 6.2 Seed Knowledge Base

```bash
python scripts/seed_knowledge.py
```

Seeds the knowledge base with:
- 40 landmark trials
- 13 therapeutic areas with success rates
- 9 regulatory agencies with pathways
- 9 endpoint types with statistical methods
- 9 adaptive designs with regulatory guidance
- 9 biomarker strategies
- 9 DCT components

### 6.3 Verify Collections

```bash
curl http://localhost:8538/collections | python -m json.tool
```

---

## 7. Data Ingest

### 7.1 Manual Ingest

```bash
# Ingest from ClinicalTrials.gov
python scripts/run_ingest.py --source clinicaltrials --query "NSCLC Phase 3"

# Ingest from PubMed
python scripts/run_ingest.py --source pubmed --query "adaptive trial design"

# Ingest regulatory documents
python scripts/run_ingest.py --source regulatory --agency FDA
```

### 7.2 Scheduled Ingest

Enable in `.env`:

```
TRIAL_INGEST_ENABLED=true
TRIAL_INGEST_SCHEDULE_HOURS=24
```

The scheduler runs ingest pipelines every 24 hours (default) via a daemon thread. It handles incremental updates and deduplication.

### 7.3 External API Keys

For increased API rate limits:

```
TRIAL_CLINICALTRIALS_API_KEY=your-key
TRIAL_NCBI_API_KEY=your-ncbi-key
```

---

## 8. Security Configuration

### 8.1 API Authentication

To enable API key authentication:

```bash
export TRIAL_API_KEY="your-secure-api-key-here"
```

Clients must include the key in requests:

```bash
curl -H "X-API-Key: your-secure-api-key-here" \
     http://localhost:8538/v1/trial/query \
     -d '{"question": "What are the key endpoints for NSCLC trials?"}'
```

### 8.2 CORS Configuration

Default CORS origins:

```
TRIAL_CORS_ORIGINS=http://localhost:8080,http://localhost:8538,http://localhost:8128
```

For production, restrict to specific domains:

```
TRIAL_CORS_ORIGINS=https://your-domain.com,https://api.your-domain.com
```

### 8.3 HTTPS (Production)

Configure a reverse proxy (nginx/traefik) for TLS termination:

```nginx
server {
    listen 443 ssl;
    server_name trial-api.your-domain.com;

    ssl_certificate /etc/ssl/certs/trial-api.crt;
    ssl_certificate_key /etc/ssl/private/trial-api.key;

    location / {
        proxy_pass http://localhost:8538;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```

### 8.4 Security Checklist

| Item | Status | Notes |
|---|---|---|
| API key authentication | Optional | Set TRIAL_API_KEY to enable |
| CORS restrictions | Configured | Default: localhost only |
| Input validation | Enforced | Pydantic models on all endpoints |
| Rate limiting | Active | 100 req/min per IP |
| Secret management | Env vars | Use vault/KMS in production |
| HTTPS | Proxy layer | Configure nginx/traefik |
| No SQL injection | N/A | Vector-only database |
| Logging | Enabled | Loguru with structured output |

---

## 9. Environment Variables

All variables use the `TRIAL_` prefix:

| Variable | Default | Required | Description |
|---|---|---|---|
| `TRIAL_MILVUS_HOST` | localhost | Yes | Milvus hostname |
| `TRIAL_MILVUS_PORT` | 19530 | Yes | Milvus port |
| `TRIAL_ANTHROPIC_API_KEY` | -- | For LLM | Anthropic API key |
| `TRIAL_API_PORT` | 8538 | No | API server port |
| `TRIAL_STREAMLIT_PORT` | 8128 | No | UI port |
| `TRIAL_API_KEY` | (empty) | No | API auth key |
| `TRIAL_CORS_ORIGINS` | localhost:* | No | CORS whitelist |
| `TRIAL_LLM_MODEL` | claude-sonnet-4-6 | No | LLM model name |
| `TRIAL_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | No | Embedding model |
| `TRIAL_EMBEDDING_DIMENSION` | 384 | No | Vector dimension |
| `TRIAL_TOP_K_PER_COLLECTION` | 5 | No | Search results per collection |
| `TRIAL_SCORE_THRESHOLD` | 0.4 | No | Minimum similarity score |
| `TRIAL_INGEST_ENABLED` | false | No | Enable scheduled ingest |
| `TRIAL_INGEST_SCHEDULE_HOURS` | 24 | No | Ingest interval |
| `TRIAL_CROSS_AGENT_TIMEOUT` | 30 | No | Cross-agent query timeout (s) |
| `TRIAL_CLINICALTRIALS_API_KEY` | -- | No | ClinicalTrials.gov API key |
| `TRIAL_NCBI_API_KEY` | -- | No | NCBI E-utilities API key |

---

## 10. Health Verification

### 10.1 API Health Check

```bash
curl -s http://localhost:8538/health | python -m json.tool
```

Expected response:

```json
{
    "status": "ok",
    "agent": "clinical-trial-intelligence",
    "version": "2.0.0",
    "milvus": "connected",
    "collections": 14,
    "knowledge_version": "2.0.0"
}
```

### 10.2 Collection Verification

```bash
curl -s http://localhost:8538/collections | python -m json.tool
```

### 10.3 Workflow Verification

```bash
curl -s http://localhost:8538/workflows | python -m json.tool
```

### 10.4 End-to-End Test

```bash
curl -X POST http://localhost:8538/v1/trial/query \
    -H "Content-Type: application/json" \
    -d '{"question": "What adaptive design should I use for a Phase 2/3 NSCLC trial?"}'
```

---

## 11. Monitoring Setup

### 11.1 Prometheus Scrape Configuration

Add to `prometheus.yml`:

```yaml
scrape_configs:
  - job_name: 'clinical-trial-agent'
    static_configs:
      - targets: ['localhost:8538']
    metrics_path: '/metrics'
    scrape_interval: 15s
```

### 11.2 Key Metrics to Monitor

| Metric | Type | Alert Threshold |
|---|---|---|
| `trial_queries_total` | Counter | Rate tracking |
| `trial_query_duration_seconds` | Histogram | p95 > 10s |
| `trial_query_errors_total` | Counter | Rate > 5/min |
| `trial_search_duration_seconds` | Histogram | p95 > 2s |
| `trial_collection_records_gauge` | Gauge | < expected count |

---

## 12. Performance Tuning

### 12.1 API Server

```bash
# Production: multiple workers
uvicorn api.main:app --host 0.0.0.0 --port 8538 --workers 4
```

### 12.2 Milvus

- Increase `nlist` for larger collections (256 for >100K records)
- Enable segment compaction for write-heavy workloads
- Configure memory limits: `queryNode.resource.limits.memory: 8Gi`

### 12.3 Embedding Generation

- Increase `TRIAL_EMBEDDING_BATCH_SIZE` for throughput (default: 32)
- Use GPU embedding if available on DGX Spark

### 12.4 Search Optimization

- Adjust `TRIAL_TOP_K_PER_COLLECTION` (lower for faster response)
- Increase `TRIAL_SCORE_THRESHOLD` to reduce noise (0.5-0.6 for precision)

---

## 13. Troubleshooting

### Common Issues

| Problem | Cause | Solution |
|---|---|---|
| API returns 503 | Milvus not running | Start Milvus: `docker start milvus-standalone` |
| Empty search results | Collections not seeded | Run `python scripts/seed_knowledge.py` |
| LLM synthesis fails | Missing API key | Set `TRIAL_ANTHROPIC_API_KEY` |
| Port conflict | Another service on 8538 | Change `TRIAL_API_PORT` |
| Import errors | Missing dependencies | Run `pip install -r requirements.txt` |
| Slow searches | Large collections, low nlist | Increase nlist, add replicas |
| Cross-agent timeout | Peer agent not running | Expected behavior; degrades gracefully |

### Log Locations

- API logs: stdout (Docker) or terminal (manual)
- Milvus logs: `docker logs milvus-standalone`
- Ingest logs: stdout from `run_ingest.py`

---

## 14. Backup and Recovery

### 14.1 Milvus Data Backup

```bash
# Backup MinIO data
docker cp minio:/data ./backup/milvus-data-$(date +%Y%m%d)

# Backup etcd metadata
docker exec etcd etcdctl snapshot save /tmp/etcd-backup.db
docker cp etcd:/tmp/etcd-backup.db ./backup/
```

### 14.2 Knowledge Base

The knowledge base is defined in code (`src/knowledge.py`) and can be re-seeded at any time:

```bash
python scripts/seed_knowledge.py
```

### 14.3 Configuration Backup

Back up the `.env` file and any custom Milvus configuration:

```bash
cp .env ./backup/.env.$(date +%Y%m%d)
```

---

*Clinical Trial Intelligence Agent v2.0.0 -- Deployment Guide -- March 2026*
