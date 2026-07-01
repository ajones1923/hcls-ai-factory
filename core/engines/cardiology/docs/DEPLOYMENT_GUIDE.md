# Cardiology Intelligence Agent -- Deployment Guide

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Docker Compose Deployment](#2-docker-compose-deployment)
3. [Manual Deployment](#3-manual-deployment)
4. [Milvus Tuning](#4-milvus-tuning)
5. [Port Map](#5-port-map)
6. [Environment Configuration](#6-environment-configuration)
7. [Security Checklist](#7-security-checklist)
8. [Monitoring with Prometheus](#8-monitoring-with-prometheus)
9. [Backup and Recovery](#9-backup-and-recovery)
10. [Scaling Considerations](#10-scaling-considerations)
11. [Updating and Maintenance](#11-updating-and-maintenance)

---

## 1. Prerequisites

### Hardware Requirements

| Component | Minimum | Recommended (DGX Spark) |
|-----------|---------|------------------------|
| CPU | 4 cores | 72 ARM cores (Grace) |
| RAM | 16 GB | 128 GB LPDDR5X |
| Storage | 50 GB SSD | 1 TB NVMe SSD |
| GPU | None (CPU inference) | NVIDIA Blackwell (128 GB) |

### Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| Docker | 24.0+ | Container runtime |
| Docker Compose | 2.20+ | Multi-service orchestration |
| Python | 3.12+ | Manual setup only |
| pip | 24.0+ | Manual setup only |

### Network Requirements

| Port | Direction | Purpose |
|------|-----------|---------|
| 8126 | Inbound | FastAPI REST API |
| 8536 | Inbound | Streamlit UI |
| 19530 | Internal | Milvus gRPC |
| 9091 | Internal | Milvus health |
| 9000/9001 | Internal | MinIO API/Console |
| 2379 | Internal | etcd |
| 443 | Outbound | Anthropic API (Claude) |
| 443 | Outbound | PubMed, ClinicalTrials.gov (ingest) |

### API Keys

| Key | Required | Purpose |
|-----|----------|---------|
| `ANTHROPIC_API_KEY` | Yes (for LLM synthesis) | Claude Sonnet 4.6 access |
| `NCBI_API_KEY` | Optional | PubMed ingest (higher rate limits) |

---

## 2. Docker Compose Deployment

### 2.1 Quick Start

```bash
# Clone the repository
cd ai_agent_adds/cardiology_intelligence_agent

# Configure environment
cp .env.example .env
# Edit .env:
#   CARDIO_ANTHROPIC_API_KEY=sk-ant-...
#   CARDIO_NCBI_API_KEY=... (optional)

# Start all services
docker compose up -d

# Watch setup logs
docker compose logs -f cardio-setup

# Verify all services are running
docker compose ps
```

### 2.2 Service Architecture

The `docker-compose.yml` defines 6 services:

```yaml
services:
  etcd:        # Milvus metadata store
  minio:       # Milvus object storage
  milvus:      # Vector database
  cardio-setup: # One-shot: create collections + seed data
  cardio-api:  # FastAPI server
  cardio-ui:   # Streamlit UI
```

### 2.3 Service Startup Order

1. **etcd** starts first (no dependencies)
2. **minio** starts first (no dependencies)
3. **milvus** starts after etcd + minio are healthy
4. **cardio-setup** starts after milvus is healthy; runs once and exits
5. **cardio-api** starts after cardio-setup completes
6. **cardio-ui** starts after cardio-api is healthy

### 2.4 Verifying Deployment

```bash
# Check all services
docker compose ps

# Expected output:
# NAME               STATUS     PORTS
# etcd               running    2379
# minio              running    9000, 9001
# milvus             running    19530, 9091
# cardio-setup       exited(0)  -
# cardio-api         running    8126
# cardio-ui          running    8536

# Health checks
curl http://localhost:8126/health
curl http://localhost:9091/healthz

# Open UI
open http://localhost:8536
```

### 2.5 Stopping and Restarting

```bash
# Stop all services (preserves data)
docker compose stop

# Start again
docker compose start

# Full teardown (removes containers, preserves volumes)
docker compose down

# Full teardown including data volumes
docker compose down -v
```

---

## 3. Manual Deployment

### 3.1 Milvus Setup

The Cardiology Intelligence Agent requires a running Milvus 2.4 instance. If you already have Milvus running (e.g., from the HCLS AI Factory), skip this step.

```bash
# Option 1: Use existing HCLS AI Factory Milvus
# Set CARDIO_MILVUS_HOST and CARDIO_MILVUS_PORT in .env

# Option 2: Start standalone Milvus
docker run -d --name milvus \
  -p 19530:19530 \
  -p 9091:9091 \
  -v milvus_data:/var/lib/milvus \
  milvusdb/milvus:v2.4-latest \
  milvus run standalone
```

### 3.2 Python Environment

```bash
cd ai_agent_adds/cardiology_intelligence_agent

# Create virtual environment
python3.12 -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Configure environment
cp .env.example .env
# Edit .env with your settings
```

### 3.3 Collection Setup

```bash
# Create all 12 cardiology collections
python scripts/setup_collections.py --drop-existing --seed

# Seed knowledge graph data
python scripts/seed_knowledge.py

# Optional: Run live data ingest
python scripts/run_ingest.py
```

### 3.4 Start Services

```bash
# Terminal 1: FastAPI server
uvicorn api.main:app --host 0.0.0.0 --port 8126 --workers 2

# Terminal 2: Streamlit UI
streamlit run app/cardio_ui.py --server.port 8536
```

---

## 4. Milvus Tuning

### 4.1 Index Configuration

All collections use IVF_FLAT with COSINE similarity. This provides high recall suitable for clinical queries at acceptable latency on DGX Spark.

```python
# Default index parameters (in src/collections.py)
index_params = {
    "index_type": "IVF_FLAT",
    "metric_type": "COSINE",
    "params": {"nlist": 128}
}
```

### 4.2 Search Parameters

```python
# Default search parameters
search_params = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16}
}
```

### 4.3 Performance Tuning

| Parameter | Default | Tuning Guidance |
|-----------|---------|----------------|
| `nlist` | 128 | Increase for larger collections (256 for >100K vectors) |
| `nprobe` | 16 | Increase for better recall (32 for clinical queries) at cost of latency |
| `top_k` | 5 | Per-collection; increase if coverage seems incomplete |
| `score_threshold` | 0.4 | Lower for broader results; raise for precision |

### 4.4 Memory Configuration

For DGX Spark (128GB RAM), configure Milvus for in-memory operation:

```yaml
# In milvus.yaml or docker-compose environment
MILVUS_QUERYNODE_CACHE_SIZE: 32GB
MILVUS_INDEXNODE_MEMORY_LIMIT: 16GB
```

### 4.5 Collection Compaction

Run periodic compaction to optimize storage:

```python
from pymilvus import utility
utility.compact("cardio_literature")
utility.wait_for_compaction("cardio_literature")
```

---

## 5. Port Map

| Service | Port | Protocol | Exposure |
|---------|------|----------|----------|
| FastAPI API | 8126 | HTTP | External (client-facing) |
| Streamlit UI | 8536 | HTTP | External (client-facing) |
| Milvus gRPC | 19530 | gRPC | Internal only |
| Milvus Health | 9091 | HTTP | Internal only |
| MinIO API | 9000 | HTTP | Internal only |
| MinIO Console | 9001 | HTTP | Internal (admin) |
| etcd | 2379 | HTTP | Internal only |
| Prometheus metrics | 8126/metrics | HTTP | Internal (monitoring) |

### Port Conflict Resolution

If ports conflict with existing services:

```bash
# Override via environment variables
export CARDIO_API_PORT=8127
export CARDIO_STREAMLIT_PORT=8528

# Or in docker-compose override
# docker-compose.override.yml
services:
  cardio-api:
    ports:
      - "8127:8126"
  cardio-ui:
    ports:
      - "8528:8536"
```

### Standalone Docker Compose (Offset Ports)

When running the Cardiology Intelligence Agent standalone (outside the main HCLS AI Factory stack), the standalone `docker-compose.yml` uses offset ports to avoid conflicts:

| Service | Standalone Port | Purpose |
|---------|----------------|---------|
| Milvus gRPC | 29530 | Offset from main Milvus (19530) |
| Milvus Health | 29091 | Offset from main Milvus health (9091) |

The agent is also integrated into the top-level `docker-compose.dgx-spark.yml` and the landing page health monitor for production deployments.

---

## 6. Environment Configuration

All settings use the `CARDIO_` prefix and are managed via Pydantic BaseSettings:

### Required Settings

```bash
# Anthropic API key (required for LLM synthesis)
CARDIO_ANTHROPIC_API_KEY=sk-ant-api03-...
```

### Optional Settings

```bash
# Milvus connection (defaults shown)
CARDIO_MILVUS_HOST=localhost
CARDIO_MILVUS_PORT=19530

# API server
CARDIO_API_HOST=0.0.0.0
CARDIO_API_PORT=8126

# Streamlit
CARDIO_STREAMLIT_PORT=8536

# Embedding model
CARDIO_EMBEDDING_MODEL=BAAI/bge-small-en-v1.5
CARDIO_EMBEDDING_DIMENSION=384

# LLM
CARDIO_LLM_MODEL=claude-sonnet-4-6

# RAG search
CARDIO_TOP_K_PER_COLLECTION=5
CARDIO_SCORE_THRESHOLD=0.4

# PubMed ingest
CARDIO_NCBI_API_KEY=your-ncbi-key
CARDIO_PUBMED_MAX_RESULTS=5000

# Monitoring
CARDIO_METRICS_ENABLED=true

# Scheduler
CARDIO_INGEST_ENABLED=false
CARDIO_INGEST_SCHEDULE_HOURS=168

# CORS (comma-separated origins)
CARDIO_CORS_ORIGINS=http://localhost:8080,http://localhost:8126,http://localhost:8536

# Citation scoring
CARDIO_CITATION_HIGH_THRESHOLD=0.75
CARDIO_CITATION_MEDIUM_THRESHOLD=0.60

# Conversation
CARDIO_MAX_CONVERSATION_CONTEXT=3
```

### Startup Validation

The `CardioSettings.validate()` method checks configuration at startup and logs warnings for:
- Missing Milvus host/port
- Missing Anthropic API key (search-only mode)
- Missing embedding model
- Port conflicts between API and Streamlit
- Collection weights not summing to ~1.0
- Non-existent RAG pipeline root directory

---

## 7. Security Checklist

### 7.1 API Key Management

- [ ] Store `ANTHROPIC_API_KEY` in `.env` file (not in code or docker-compose.yml)
- [ ] Add `.env` to `.gitignore`
- [ ] Use environment variables or secrets manager in production
- [ ] Rotate API keys periodically
- [ ] Set `NCBI_API_KEY` separately (different rotation schedule)

### 7.2 Network Security

- [ ] Restrict external access to ports 8126 (API) and 8536 (UI) only
- [x] API authentication via X-API-Key header (enabled by default)
- [ ] Keep Milvus (19530), etcd (2379), MinIO (9000) on internal network only
- [ ] Use reverse proxy (nginx/traefik) with TLS for production
- [ ] Enable CORS restrictions for production domains
- [x] Rate limit API endpoints (100 req/min per IP -- enabled by default)
- [x] API authentication via X-API-Key header

### 7.3 Data Security

- [ ] Enable Milvus authentication if available
- [ ] Encrypt MinIO storage at rest
- [ ] Audit log API access
- [ ] No PHI/PII should be stored in vector collections
- [ ] All patient data in queries is transient (not persisted)

### 7.4 Application Security

- [ ] Set `MAX_REQUEST_SIZE_MB` to prevent oversized payloads (default: 10)
- [ ] Validate all input via Pydantic models (automatic)
- [ ] Review CORS_ORIGINS for production (remove localhost)
- [ ] Disable debug mode in production

### 7.5 Container Security

- [ ] Use specific image tags (not `latest`) in docker-compose.yml
- [ ] Run containers as non-root user
- [ ] Limit container resource usage (CPU, memory limits)
- [ ] Scan images for vulnerabilities
- [ ] Keep base images updated

---

## 8. Monitoring with Prometheus

### 8.1 Metrics Endpoint

The FastAPI server exposes Prometheus metrics at:

```
GET http://localhost:8126/metrics
```

### 8.2 Available Metrics

| Metric | Type | Description |
|--------|------|-------------|
| `cardio_queries_total` | Counter | Total queries processed |
| `cardio_query_duration_seconds` | Histogram | Query processing latency |
| `cardio_risk_calculations_total` | Counter | Risk calculator invocations (by type) |
| `cardio_gdmt_optimizations_total` | Counter | GDMT optimization requests |
| `cardio_collection_search_duration` | Histogram | Per-collection search latency |
| `cardio_cross_modal_triggers_total` | Counter | Cross-modal triggers fired |
| `cardio_workflow_executions_total` | Counter | Workflow executions (by type) |
| `cardio_export_requests_total` | Counter | Export requests (by format) |
| `cardio_errors_total` | Counter | Error count (by type) |

### 8.3 Prometheus Configuration

Add to your Prometheus `prometheus.yml`:

```yaml
scrape_configs:
  - job_name: 'cardio-agent'
    scrape_interval: 15s
    static_configs:
      - targets: ['localhost:8126']
    metrics_path: '/metrics'
```

### 8.4 Grafana Dashboard

Create a Grafana dashboard with panels for:
- Query rate (queries per minute)
- Query latency (p50, p95, p99)
- Risk calculator usage distribution
- GDMT optimization rate
- Error rate
- Collection search latency by collection
- Cross-modal trigger rate

### 8.5 Alerting Rules

Suggested Prometheus alerting rules:

```yaml
groups:
  - name: cardio-agent
    rules:
      - alert: HighErrorRate
        expr: rate(cardio_errors_total[5m]) > 0.1
        for: 5m
        annotations:
          summary: "Cardio agent error rate above 10%"

      - alert: HighQueryLatency
        expr: histogram_quantile(0.95, rate(cardio_query_duration_seconds_bucket[5m])) > 10
        for: 5m
        annotations:
          summary: "Cardio agent p95 latency above 10 seconds"

      - alert: MilvusDown
        expr: up{job="milvus"} == 0
        for: 1m
        annotations:
          summary: "Milvus vector database is down"
```

---

## 9. Backup and Recovery

### 9.1 Milvus Data Backup

```bash
# Export collection data
python -c "
from pymilvus import utility
collections = utility.list_collections()
for col in collections:
    if col.startswith('cardio_'):
        print(f'Backing up {col}...')
        # Use Milvus backup utility or pymilvus export
"

# Volume-level backup (Docker)
docker run --rm -v milvus_data:/data -v $(pwd)/backup:/backup \
  alpine tar czf /backup/milvus-data-$(date +%Y%m%d).tar.gz /data
```

### 9.2 Configuration Backup

```bash
# Backup configuration and environment
cp .env .env.backup.$(date +%Y%m%d)
cp docker-compose.yml docker-compose.yml.backup.$(date +%Y%m%d)
```

### 9.3 Recovery Procedure

```bash
# 1. Stop services
docker compose down

# 2. Restore Milvus data (if needed)
docker run --rm -v milvus_data:/data -v $(pwd)/backup:/backup \
  alpine tar xzf /backup/milvus-data-YYYYMMDD.tar.gz -C /

# 3. Restore configuration
cp .env.backup.YYYYMMDD .env

# 4. Start services
docker compose up -d

# 5. Verify
curl http://localhost:8126/health
```

---

## 10. Scaling Considerations

### 10.1 Single-Node (DGX Spark)

The default deployment is optimized for single-node operation:
- Milvus standalone mode
- 2 Uvicorn workers for the API
- Single Streamlit instance
- All services on one machine

### 10.2 Multi-Worker API

To increase API throughput:

```bash
# Increase Uvicorn workers
uvicorn api.main:app --host 0.0.0.0 --port 8126 --workers 4
```

### 10.3 Milvus Cluster Mode

For larger deployments with millions of vectors per collection:

1. Deploy Milvus in cluster mode (separate query nodes, data nodes, index nodes)
2. Increase `nlist` and `nprobe` parameters
3. Consider HNSW index for sub-millisecond search at the cost of memory

### 10.4 Read Replicas

For high query volume:
- Deploy multiple API server instances behind a load balancer
- All instances share the same Milvus cluster
- Streamlit can be deployed separately for different user groups

---

## 11. Updating and Maintenance

### 11.1 Code Updates

```bash
# Pull latest code
git pull

# Rebuild and restart
docker compose build
docker compose up -d
```

### 11.2 Knowledge Base Updates

```bash
# Re-run knowledge seed (idempotent)
python scripts/seed_knowledge.py

# Run live ingest for latest data
python scripts/run_ingest.py
```

### 11.3 Dependency Updates

```bash
# Update Python dependencies
pip install -r requirements.txt --upgrade

# Update Docker images
docker compose pull
docker compose up -d
```

### 11.4 Milvus Upgrades

Follow the Milvus upgrade guide. Key steps:
1. Backup all collection data
2. Stop Milvus
3. Update the Milvus image version in docker-compose.yml
4. Start Milvus (automatic schema migration)
5. Verify collections and vector counts

### 11.5 Scheduled Ingest

To enable automatic weekly data updates:

```bash
# Enable in .env
CARDIO_INGEST_ENABLED=true
CARDIO_INGEST_SCHEDULE_HOURS=168  # Weekly

# Restart API server
docker compose restart cardio-api
```

The scheduler runs all 7 ingest parsers (PubMed, ClinicalTrials.gov, imaging, ECG, guidelines, devices, hemodynamics) on the configured interval.
