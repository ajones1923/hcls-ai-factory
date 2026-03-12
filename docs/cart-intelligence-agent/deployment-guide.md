# CAR-T Intelligence Agent — Deployment Guide

**HCLS AI Factory / ai_agent_adds / cart_intelligence_agent**

Version 1.0.0 | March 2026 | Author: Adam Jones

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
7. [Data Ingestion (PubMed, ClinicalTrials.gov)](#7-data-ingestion-pubmed-clinicaltrialsgov)
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

The CAR-T Intelligence Agent is a multi-collection RAG-powered clinical decision
support system for CAR-T cell therapy. It combines 11 Milvus vector collections
(10 agent-owned + 1 shared), a 71-node knowledge graph, BGE-small-en-v1.5
sentence embeddings (384-dim), and Claude LLM reasoning to deliver
evidence-based intelligence spanning literature, clinical trials, constructs,
assays, manufacturing, safety, biomarkers, regulatory, sequences, and
real-world evidence.

### Core Capabilities

- **Multi-collection RAG search** across 11 knowledge domains with weighted
  relevance scoring and antigen-specific filtering.
- **Knowledge graph augmentation** with 71 curated nodes (34 targets, 17
  toxicities, 20 manufacturing processes).
- **Query expansion** with 229 keyword categories mapping to 1,961 expansion
  terms across 12 domain-specific maps for improved recall.
- **Deep Research mode** with autonomous multi-pass agent pipeline.
- **Cross-collection entity linking** for holistic product intelligence.
- **PDF/Markdown/JSON report export** for clinical and research use.
- **Prometheus-compatible metrics** for operational monitoring.

### Architecture at a Glance

The agent runs as 6 Docker services on a bridge network (`cart-network`):

```
                    cart-network (bridge)
    +--------------------------------------------------+
    |                                                  |
    |  milvus-etcd -----> milvus-standalone <-----+    |
    |  milvus-minio ---/    :19530 :9091          |    |
    |                                             |    |
    |  cart-api (:8522) --------------------------+    |
    |  cart-streamlit (:8521) ---------------------+    |
    |  cart-setup (one-shot) ----------------------+    |
    |                                                  |
    +--------------------------------------------------+
```

| Service | Container Name | Purpose |
|---------|---------------|---------|
| `milvus-etcd` | `cart-milvus-etcd` | Metadata key-value store for Milvus |
| `milvus-minio` | `cart-milvus-minio` | Object storage for Milvus log/index data |
| `milvus-standalone` | `cart-milvus-standalone` | Vector database (Milvus 2.4) |
| `cart-streamlit` | `cart-streamlit` | Streamlit clinical chat UI (port 8521) |
| `cart-api` | `cart-api` | FastAPI REST API (port 8522) |
| `cart-setup` | `cart-setup` | One-shot: create collections + seed all data |

---

## 2. Prerequisites

### Hardware

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 4 cores | 8+ cores |
| RAM | 16 GB | 32 GB |
| Disk | 20 GB | 50 GB (with Milvus data) |
| GPU | Not required | NVIDIA GPU (for faster embedding) |

The agent runs entirely on CPU. GPU acceleration is beneficial for initial
embedding generation but not required for runtime.

### Software

| Dependency | Version | Purpose |
|------------|---------|---------|
| Docker Engine | 24.0+ | Container runtime |
| Docker Compose | 2.20+ | Multi-container orchestration |
| Python | 3.10+ | Development mode only |
| NVIDIA DGX Spark | Grace Blackwell | Production deployment target |

### API Keys

| Key | Required | Source |
|-----|----------|--------|
| `ANTHROPIC_API_KEY` | Yes (for LLM) | [console.anthropic.com](https://console.anthropic.com) |
| `NCBI_API_KEY` | Optional | [ncbi.nlm.nih.gov/account](https://www.ncbi.nlm.nih.gov/account/) |

The agent functions without an API key for embedding-only search (`/search`
endpoint), but requires `ANTHROPIC_API_KEY` for full RAG queries with LLM
synthesis.

---

## 3. Quick Start

```bash
# 1. Clone the repository
cd /path/to/hcls-ai-factory/ai_agent_adds/cart_intelligence_agent

# 2. Create environment file
cp .env.example .env
# Edit .env and add: ANTHROPIC_API_KEY=sk-ant-...

# 3. Start all services
docker compose up -d

# 4. Watch setup progress (one-shot container seeds all collections)
docker compose logs -f cart-setup

# 5. Access the interfaces
#    Streamlit UI: http://localhost:8521
#    FastAPI docs: http://localhost:8522/docs
#    Health check: http://localhost:8522/health
```

Startup sequence:
1. `milvus-etcd` + `milvus-minio` start first (health-checked)
2. `milvus-standalone` starts after etcd + minio are healthy (~60s)
3. `cart-setup` runs once: creates 11 collections and seeds 649 records from 13 JSON files
4. `cart-streamlit` + `cart-api` start after Milvus is healthy

Total cold-start time: ~3–5 minutes (including Milvus init and embedding model download).

---

## 4. Deployment Modes

### 4a. Docker Lite (Milvus + API Only)

For headless/API-only deployments without the Streamlit UI:

```bash
docker compose up -d milvus-etcd milvus-minio milvus-standalone cart-api cart-setup
```

This starts only the vector database infrastructure and the REST API. Useful for
backend integration where another frontend consumes the API.

### 4b. Docker Full Stack (All 6 Services)

The default `docker compose up -d` launches all 6 services. This is the
recommended deployment for demos and development.

```bash
docker compose up -d
docker compose ps   # Verify all services are running
```

### 4c. DGX Spark Production

On the NVIDIA DGX Spark ($3,999), the CAR-T agent runs as part of the full
HCLS AI Factory stack via the master `docker-compose.dgx-spark.yml`:

```bash
# From the HCLS AI Factory root
cd $HCLS_HOME
docker compose -f docker-compose.dgx-spark.yml up -d

# Or start just the CAR-T services
docker compose -f docker-compose.dgx-spark.yml up -d \
  cart-milvus-etcd cart-milvus-minio cart-milvus-standalone \
  cart-streamlit cart-api cart-setup
```

In production, the agent integrates with:
- **Landing Page** (`:8080`) — service health dashboard
- **Health Monitor** (`health-monitor.sh`) — auto-restart on failure
- **Prometheus + Grafana** — metrics collection and visualization
- **Cross-agent event bus** (`hcls_common.event_bus`) — inter-agent communication

### 4d. Development Mode (Local Python)

For development without Docker:

```bash
# 1. Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# 2. Install dependencies
pip install -r requirements.txt

# 3. Start Milvus (still needs Docker)
docker compose up -d milvus-etcd milvus-minio milvus-standalone

# 4. Set environment variables
export ANTHROPIC_API_KEY=sk-ant-...
export CART_MILVUS_HOST=localhost
export CART_MILVUS_PORT=19530

# 5. Seed collections
python3 scripts/setup_collections.py --drop-existing --seed-constructs
python3 scripts/seed_knowledge.py
python3 scripts/seed_assays.py
python3 scripts/seed_manufacturing.py
python3 scripts/seed_safety.py
python3 scripts/seed_biomarkers.py
python3 scripts/seed_regulatory.py
python3 scripts/seed_sequences.py
python3 scripts/seed_realworld.py

# 6. Launch Streamlit UI
streamlit run app/cart_ui.py --server.port 8521

# 7. (Optional) Launch FastAPI server
uvicorn api.main:app --host 0.0.0.0 --port 8522 --reload

# 8. Run tests
python3 -m pytest tests/ -v
```

---

## 5. Configuration Reference

All configuration is managed through `config/settings.py` using Pydantic
BaseSettings with `env_prefix="CART_"`. Every setting can be overridden via
environment variables prefixed with `CART_`.

### 5.1 Path Settings

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `PROJECT_ROOT` | Path | (auto-detected) | -- |
| `DATA_DIR` | Path | `{PROJECT_ROOT}/data` | -- |
| `CACHE_DIR` | Path | `{DATA_DIR}/cache` | -- |
| `REFERENCE_DIR` | Path | `{DATA_DIR}/reference` | -- |
| `RAG_PIPELINE_ROOT` | Path | `/app/rag-chat-pipeline` | `CART_RAG_PIPELINE_ROOT` |

### 5.2 Milvus Settings

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `MILVUS_HOST` | str | `localhost` | `CART_MILVUS_HOST` |
| `MILVUS_PORT` | int | `19530` | `CART_MILVUS_PORT` |

### 5.3 Collection Names

| Variable | Default Collection Name |
|----------|------------------------|
| `COLLECTION_LITERATURE` | `cart_literature` |
| `COLLECTION_TRIALS` | `cart_trials` |
| `COLLECTION_CONSTRUCTS` | `cart_constructs` |
| `COLLECTION_ASSAYS` | `cart_assays` |
| `COLLECTION_MANUFACTURING` | `cart_manufacturing` |
| `COLLECTION_GENOMIC` | `genomic_evidence` |
| `COLLECTION_SAFETY` | `cart_safety` |
| `COLLECTION_BIOMARKERS` | `cart_biomarkers` |
| `COLLECTION_REGULATORY` | `cart_regulatory` |
| `COLLECTION_SEQUENCES` | `cart_sequences` |
| `COLLECTION_REALWORLD` | `cart_realworld` |

### 5.4 Embedding Settings

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `EMBEDDING_MODEL` | str | `BAAI/bge-small-en-v1.5` | `CART_EMBEDDING_MODEL` |
| `EMBEDDING_DIMENSION` | int | `384` | `CART_EMBEDDING_DIMENSION` |
| `EMBEDDING_BATCH_SIZE` | int | `32` | `CART_EMBEDDING_BATCH_SIZE` |

### 5.5 LLM Settings

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `LLM_PROVIDER` | str | `anthropic` | `CART_LLM_PROVIDER` |
| `LLM_MODEL` | str | `claude-sonnet-4-6` | `CART_LLM_MODEL` |
| `ANTHROPIC_API_KEY` | str | `None` | `CART_ANTHROPIC_API_KEY` |

### 5.6 RAG Search Settings

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `TOP_K_PER_COLLECTION` | int | `5` | `CART_TOP_K_PER_COLLECTION` |
| `SCORE_THRESHOLD` | float | `0.4` | `CART_SCORE_THRESHOLD` |

**Collection Search Weights** (must sum to ~1.0):

| Collection | Weight |
|------------|--------|
| Literature | 0.20 |
| Trials | 0.16 |
| Constructs | 0.10 |
| Assays | 0.09 |
| Safety | 0.08 |
| Biomarkers | 0.08 |
| Manufacturing | 0.07 |
| Real-World | 0.07 |
| Regulatory | 0.06 |
| Sequences | 0.06 |
| Genomic | 0.04 |

### 5.7 API Settings

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `API_HOST` | str | `0.0.0.0` | `CART_API_HOST` |
| `API_PORT` | int | `8522` | `CART_API_PORT` |
| `STREAMLIT_PORT` | int | `8521` | `CART_STREAMLIT_PORT` |
| `CORS_ORIGINS` | str | `http://localhost:8080,http://localhost:8521,http://localhost:8522` | `CART_CORS_ORIGINS` |
| `MAX_REQUEST_SIZE_MB` | int | `10` | `CART_MAX_REQUEST_SIZE_MB` |

### 5.8 Monitoring & Scheduling

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `METRICS_ENABLED` | bool | `True` | `CART_METRICS_ENABLED` |
| `INGEST_SCHEDULE_HOURS` | int | `168` (weekly) | `CART_INGEST_SCHEDULE_HOURS` |
| `INGEST_ENABLED` | bool | `False` | `CART_INGEST_ENABLED` |

### 5.9 Conversation & Citation

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `MAX_CONVERSATION_CONTEXT` | int | `3` | `CART_MAX_CONVERSATION_CONTEXT` |
| `CITATION_HIGH_THRESHOLD` | float | `0.75` | `CART_CITATION_HIGH_THRESHOLD` |
| `CITATION_MEDIUM_THRESHOLD` | float | `0.60` | `CART_CITATION_MEDIUM_THRESHOLD` |

### 5.10 PubMed & ClinicalTrials.gov

| Variable | Type | Default | Env Override |
|----------|------|---------|-------------|
| `NCBI_API_KEY` | str | `None` | `CART_NCBI_API_KEY` |
| `PUBMED_MAX_RESULTS` | int | `5000` | `CART_PUBMED_MAX_RESULTS` |
| `CT_GOV_BASE_URL` | str | `https://clinicaltrials.gov/api/v2` | `CART_CT_GOV_BASE_URL` |

---

## 6. Collection Setup and Seeding

### 6.1 Collections

The agent operates across 11 Milvus collections:

| Collection | Description | Schema Fields |
|------------|-------------|---------------|
| `cart_literature` | Published research papers | title, abstract, target_antigen, year, journal |
| `cart_trials` | Clinical trial records | nct_id, title, phase, status, target_antigen |
| `cart_constructs` | CAR construct designs | name, target, costimulatory_domain, generation |
| `cart_assays` | Laboratory assay protocols | name, assay_type, target, sensitivity |
| `cart_manufacturing` | Manufacturing processes | process_name, scale, duration, yield |
| `cart_safety` | Safety/toxicity profiles | event_type, grade, management, product |
| `cart_biomarkers` | Predictive biomarkers | name, category, clinical_utility |
| `cart_regulatory` | Regulatory approvals | product, indication, agency, approval_date |
| `cart_sequences` | CAR/TCR sequences | sequence_type, target, framework |
| `cart_realworld` | Real-world evidence | study_type, endpoints, population |
| `genomic_evidence` | Shared genomic variants | (shared with genomics pipeline) |

All agent-owned collections use:
- **Embedding dimension:** 384 (BGE-small-en-v1.5)
- **Index type:** IVF_FLAT
- **Metric:** COSINE
- **nlist:** 128

### 6.2 Seed Data

The `data/seed/` directory contains 13 JSON files with 649 curated records:

| File | Records | Description |
|------|---------|-------------|
| `assay_seed_data.json` | 75 | Laboratory assay protocols |
| `biomarker_seed_data.json` | 60 | Predictive and prognostic biomarkers |
| `constructs_seed_data.json` | 41 | CAR construct designs |
| `immunogenicity_biomarker_seed.json` | 20 | Immunogenicity biomarkers |
| `immunogenicity_sequence_seed.json` | 18 | Immunogenicity sequence data |
| `literature_seed_data.json` | 60 | Published research papers |
| `manufacturing_seed_data.json` | 56 | Manufacturing process records |
| `patent_seed_data.json` | 45 | Patent records |
| `realworld_seed_data.json` | 54 | Real-world evidence studies |
| `regulatory_seed_data.json` | 40 | Regulatory approval records |
| `safety_seed_data.json` | 71 | Safety and toxicity profiles |
| `sequence_seed_data.json` | 40 | CAR/TCR sequence records |
| `trials_seed_data.json` | 69 | Clinical trial records |

### 6.3 Manual Seeding

If the `cart-setup` container didn't run or you need to re-seed:

```bash
# Re-run the setup container
docker compose run --rm cart-setup

# Or seed individually
docker compose exec cart-api python scripts/seed_knowledge.py
docker compose exec cart-api python scripts/seed_assays.py
docker compose exec cart-api python scripts/seed_manufacturing.py
# ... etc.
```

---

## 7. Data Ingestion (PubMed, ClinicalTrials.gov)

The agent includes automated ingestion pipelines for external data sources.

### 7.1 PubMed Ingestion

```bash
# One-time ingest
docker compose exec cart-api python scripts/ingest_pubmed.py \
  --query "CAR-T cell therapy" \
  --max-results 5000

# With NCBI API key (higher rate limit)
CART_NCBI_API_KEY=your_key docker compose exec cart-api \
  python scripts/ingest_pubmed.py --query "CAR-T cell therapy"
```

### 7.2 ClinicalTrials.gov Ingestion

```bash
docker compose exec cart-api python scripts/ingest_trials.py \
  --condition "CAR-T" \
  --status "RECRUITING,ACTIVE_NOT_RECRUITING"
```

### 7.3 Scheduled Ingestion

Enable automatic weekly ingestion:

```bash
# In .env
CART_INGEST_ENABLED=true
CART_INGEST_SCHEDULE_HOURS=168  # Weekly
```

The scheduler runs inside the API container and triggers ingestion pipelines
at the configured interval.

---

## 8. Networking and Ports

### 8.1 Port Map

| Port | Service | Protocol | Purpose |
|------|---------|----------|---------|
| 8521 | cart-streamlit | HTTP | Streamlit chat UI |
| 8522 | cart-api | HTTP | FastAPI REST API |
| 19530 | milvus-standalone | gRPC | Milvus client connections |
| 9091 | milvus-standalone | HTTP | Milvus health/metrics |
| 2379 | milvus-etcd | HTTP | etcd client (internal only) |
| 9000 | milvus-minio | HTTP | MinIO API (internal only) |
| 9001 | milvus-minio | HTTP | MinIO console (internal only) |

### 8.2 Network Configuration

All services communicate over the `cart-network` Docker bridge network.
Internal services (etcd, minio) are not exposed to the host by default.

To expose MinIO console for debugging:

```yaml
# In docker-compose.yml, add to milvus-minio:
ports:
  - "9001:9001"
```

### 8.3 CORS Configuration

CORS is restricted to 3 origins by default:

```
http://localhost:8080   # HCLS AI Factory Landing Page
http://localhost:8521   # CAR-T Streamlit UI
http://localhost:8522   # CAR-T API (self)
```

Override with `CART_CORS_ORIGINS` (comma-separated):

```bash
CART_CORS_ORIGINS="http://localhost:8080,http://your-domain.com"
```

---

## 9. Storage and Persistence

### 9.1 Docker Volumes

| Volume | Mount Point | Purpose | Typical Size |
|--------|-------------|---------|-------------|
| `etcd_data` | `/etcd` | Milvus metadata | ~100 MB |
| `minio_data` | `/minio_data` | Milvus index/log data | ~2 GB |
| `milvus_data` | `/var/lib/milvus` | Milvus WAL + data | ~5 GB |

### 9.2 Application Data

| Path | Purpose | Persistence |
|------|---------|-------------|
| `data/seed/` | Seed JSON files | Baked into image |
| `data/cache/` | Embedding cache | Ephemeral (recreated on restart) |
| `data/reference/` | Reference data | Baked into image |

### 9.3 Volume Backup

```bash
# Backup all Milvus volumes
docker compose stop
for vol in etcd_data minio_data milvus_data; do
  docker run --rm -v cart_intelligence_agent_${vol}:/data \
    -v $(pwd)/backups:/backup alpine \
    tar czf /backup/${vol}_$(date +%Y%m%d).tar.gz -C /data .
done
docker compose start
```

---

## 10. Monitoring and Metrics

### 10.1 Prometheus Endpoint

The API server exposes Prometheus-compatible metrics at `GET /metrics`:

```
# HELP cart_api_requests_total Total API requests
# TYPE cart_api_requests_total counter
cart_api_requests_total 42

# HELP cart_api_query_requests_total Total /query requests
# TYPE cart_api_query_requests_total counter
cart_api_query_requests_total 15

# HELP cart_collection_vectors Number of vectors per collection
# TYPE cart_collection_vectors gauge
cart_collection_vectors{collection="cart_literature"} 1200
cart_collection_vectors{collection="cart_trials"} 850
```

### 10.2 Health Endpoint

```bash
curl http://localhost:8522/health
# {"status":"healthy","collections":10,"total_vectors":6452}
```

### 10.3 Grafana Integration

On the DGX Spark, metrics flow to the central Prometheus instance and are
visualized on the HCLS AI Factory Grafana dashboard.

Prometheus scrape config:

```yaml
- job_name: 'cart-agent'
  static_configs:
    - targets: ['localhost:8522']
  metrics_path: '/metrics'
  scrape_interval: 30s
```

---

## 11. Security Hardening

### 11.1 Container Security

The Dockerfile implements security best practices:

- **Multi-stage build** — build dependencies excluded from runtime image
- **Non-root user** — `cartuser` with minimal permissions
- **Read-only filesystem** — application code is not writable at runtime
- **No shell access** — `cartuser` has `/bin/false` as shell

### 11.2 API Security

- **CORS restrictions** — limited to 3 explicit origins (not `["*"]`)
- **Request size limits** — `MAX_REQUEST_SIZE_MB=10` prevents payload abuse
- **Input validation** — Pydantic schemas validate all request bodies
- **No API key in responses** — API key is never reflected in output

### 11.3 Network Security

- Internal services (etcd, minio) are not exposed to the host
- Milvus gRPC port (19530) should be firewalled in production
- Use a reverse proxy (nginx, Traefik) for TLS termination

### 11.4 API Key Management

```bash
# Development: .env file
ANTHROPIC_API_KEY=sk-ant-...

# Production: Docker secrets or environment injection
docker compose run -e ANTHROPIC_API_KEY=$(vault read ...) cart-api
```

---

## 12. Health Checks and Troubleshooting

### 12.1 Service Health Checks

| Service | Health Check | Interval | Start Period |
|---------|-------------|----------|-------------|
| milvus-etcd | `etcdctl endpoint health` | 30s | -- |
| milvus-minio | `curl http://localhost:9000/minio/health/live` | 30s | -- |
| milvus-standalone | `curl http://localhost:9091/healthz` | 30s | 60s |
| cart-api | `curl http://localhost:8522/health` | 30s | 30s |
| cart-streamlit | `curl http://localhost:8521/_stcore/health` | 30s | 40s |

### 12.2 Common Issues

**Milvus fails to start:**
```bash
# Check etcd and minio are healthy first
docker compose ps
docker compose logs milvus-etcd
docker compose logs milvus-minio

# Common fix: remove stale lock files
docker compose down -v  # WARNING: deletes all data
docker compose up -d
```

**"Engine not initialized" on /health:**
```bash
# The embedding model is still downloading (~90 MB on first run)
docker compose logs cart-api | grep -i "model"
# Wait 1–2 minutes for sentence-transformers to download
```

**"LLM client not available":**
```bash
# Check API key is set
docker compose exec cart-api env | grep ANTHROPIC
# Verify key is valid
curl -H "x-api-key: $ANTHROPIC_API_KEY" https://api.anthropic.com/v1/messages \
  -d '{"model":"claude-sonnet-4-6","max_tokens":10,"messages":[{"role":"user","content":"hi"}]}'
```

**Empty search results:**
```bash
# Verify collections are seeded
curl http://localhost:8522/collections
# Re-run setup if needed
docker compose run --rm cart-setup
```

### 12.3 Log Inspection

```bash
# All services
docker compose logs -f --tail=100

# Specific service
docker compose logs -f cart-api
docker compose logs -f cart-streamlit

# Setup progress
docker compose logs -f cart-setup
```

---

## 13. Backup and Recovery

### 13.1 Milvus Data Backup

```bash
# Stop writes (optional, for consistency)
docker compose pause cart-api cart-streamlit

# Backup volumes
docker run --rm \
  -v cart_intelligence_agent_milvus_data:/source:ro \
  -v $(pwd)/backups:/backup \
  alpine tar czf /backup/milvus_$(date +%Y%m%d).tar.gz -C /source .

docker compose unpause cart-api cart-streamlit
```

### 13.2 Full Backup

```bash
# Backup everything: volumes + seed data + config
tar czf cart_backup_$(date +%Y%m%d).tar.gz \
  data/seed/ config/ .env docker-compose.yml
```

### 13.3 Restore

```bash
# Restore Milvus volume
docker compose down
docker volume rm cart_intelligence_agent_milvus_data
docker volume create cart_intelligence_agent_milvus_data
docker run --rm \
  -v cart_intelligence_agent_milvus_data:/target \
  -v $(pwd)/backups:/backup \
  alpine tar xzf /backup/milvus_YYYYMMDD.tar.gz -C /target
docker compose up -d
```

---

## 14. Scaling Considerations

### 14.1 Vertical Scaling

On the DGX Spark (128 GB unified memory), the agent can handle:
- Up to ~500K vectors per collection without performance degradation
- Concurrent query load of ~50 requests/second on the API
- Embedding batch sizes up to 256 for bulk ingestion

### 14.2 Horizontal Scaling

For higher throughput:

```yaml
# Scale API workers
cart-api:
  command:
    - uvicorn
    - api.main:app
    - --host=0.0.0.0
    - --port=8522
    - --workers=4  # Increase from 2
```

### 14.3 Milvus Scaling

For production workloads exceeding single-node Milvus capacity:
- Migrate from `milvus-standalone` to Milvus distributed mode
- Use external etcd cluster and S3-compatible object storage
- Refer to [Milvus documentation](https://milvus.io/docs) for cluster deployment

---

## 15. Integration with HCLS AI Factory

### 15.1 Landing Page Registration

The CAR-T agent registers with the HCLS AI Factory landing page (`:8080`)
for service health monitoring:

| Endpoint | Monitored By |
|----------|-------------|
| `http://localhost:8522/health` | `health-monitor.sh` |
| `http://localhost:8521/_stcore/health` | `health-monitor.sh` |

### 15.2 Cross-Agent Event Bus

The agent publishes events via `hcls_common.event_bus`:

```python
from hcls_common.event_bus import publish_event, EventType, PipelineStage

publish_event(
    EventType.CART_MANUFACTURING_READY,
    source_stage=PipelineStage.CART_ANALYSIS,
    payload={
        "target_antigen": "CD19",
        "evidence_count": 42,
        "mode": "deep_research",
    },
)
```

### 15.3 Shared Collections

The `genomic_evidence` collection is shared with the genomics pipeline
(Stage 1). The CAR-T agent reads from this collection but does not write to it.

### 15.4 Master Docker Compose

In production, the CAR-T services are defined in the master
`docker-compose.dgx-spark.yml` alongside all other HCLS AI Factory services.

---

## 16. Updating and Maintenance

### 16.1 Code Updates

```bash
# Pull latest code
git pull origin main

# Rebuild containers
docker compose build --no-cache
docker compose up -d
```

### 16.2 Dependency Updates

```bash
# Update requirements
pip install --upgrade -r requirements.txt

# Rebuild Docker image
docker compose build cart-api cart-streamlit
docker compose up -d cart-api cart-streamlit
```

### 16.3 Knowledge Updates

To refresh the knowledge graph and seed data:

```bash
# Re-seed with latest data
docker compose run --rm cart-setup

# Or update specific collections
docker compose exec cart-api python scripts/seed_knowledge.py
```

### 16.4 Milvus Updates

```bash
# Backup first
./scripts/backup_milvus.sh

# Update Milvus image
docker compose pull milvus-standalone
docker compose up -d milvus-standalone
```

---

## Appendix A: Complete docker-compose.yml

The complete `docker-compose.yml` is included at the project root. Key
parameters:

```yaml
version: "3.8"

services:
  milvus-etcd:
    image: quay.io/coreos/etcd:v3.5.5
    # Health-checked, 4 GB backend quota

  milvus-minio:
    image: minio/minio:RELEASE.2023-03-20T20-16-18Z
    # Health-checked, default credentials

  milvus-standalone:
    image: milvusdb/milvus:v2.4-latest
    ports: ["19530:19530", "9091:9091"]
    # Depends on etcd + minio healthy, 60s start period

  cart-streamlit:
    build: .
    ports: ["8521:8521"]
    # Default CMD: streamlit run app/cart_ui.py

  cart-api:
    build: .
    ports: ["8522:8522"]
    # CMD override: uvicorn api.main:app

  cart-setup:
    build: .
    restart: "no"
    # One-shot: create collections + seed all data

volumes:
  etcd_data:
  minio_data:
  milvus_data:

networks:
  cart-network:
    driver: bridge
```

---

## Appendix B: Environment Variable Quick Reference

All variables use the `CART_` prefix for Pydantic settings injection.

| Variable | Required | Default | Description |
|----------|----------|---------|-------------|
| `ANTHROPIC_API_KEY` | Yes* | -- | Claude API key (also accepts `CART_ANTHROPIC_API_KEY`) |
| `CART_MILVUS_HOST` | No | `localhost` | Milvus server hostname |
| `CART_MILVUS_PORT` | No | `19530` | Milvus server port |
| `CART_API_PORT` | No | `8522` | FastAPI server port |
| `CART_STREAMLIT_PORT` | No | `8521` | Streamlit UI port |
| `CART_LLM_MODEL` | No | `claude-sonnet-4-6` | Claude model ID |
| `CART_EMBEDDING_MODEL` | No | `BAAI/bge-small-en-v1.5` | Sentence embedding model |
| `CART_TOP_K_PER_COLLECTION` | No | `5` | Results per collection |
| `CART_SCORE_THRESHOLD` | No | `0.4` | Minimum similarity score |
| `CART_CORS_ORIGINS` | No | `localhost:8080,8521,8522` | Allowed CORS origins |
| `CART_METRICS_ENABLED` | No | `True` | Enable Prometheus metrics |
| `CART_INGEST_ENABLED` | No | `False` | Enable scheduled ingestion |
| `CART_INGEST_SCHEDULE_HOURS` | No | `168` | Ingestion interval (hours) |
| `CART_MAX_REQUEST_SIZE_MB` | No | `10` | Max request body size |
| `CART_RAG_PIPELINE_ROOT` | No | `/app/rag-chat-pipeline` | Path to RAG pipeline |
| `CART_NCBI_API_KEY` | No | -- | NCBI API key (PubMed) |
| `HCLS_LIB_PATH` | No | `/app/lib` | Path to hcls_common library |

\* Required for LLM-powered queries. Embedding-only search works without it.

---

**Author:** Adam Jones
**License:** Apache 2.0
**Repository:** [github.com/ajones1923/hcls-ai-factory](https://github.com/ajones1923/hcls-ai-factory)
