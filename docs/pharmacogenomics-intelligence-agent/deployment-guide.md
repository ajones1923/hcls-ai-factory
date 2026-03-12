# Pharmacogenomics Intelligence Agent -- Deployment Guide

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [Deployment Options](#1-deployment-options)
2. [Docker Compose Deployment](#2-docker-compose-deployment)
3. [Manual Deployment](#3-manual-deployment)
4. [Milvus Tuning](#4-milvus-tuning)
5. [GPU Considerations](#5-gpu-considerations)
6. [Production Security Checklist](#6-production-security-checklist)
7. [Monitoring Setup](#7-monitoring-setup)
8. [Backup and Disaster Recovery](#8-backup-and-disaster-recovery)
9. [Scaling Strategies](#9-scaling-strategies)
10. [Maintenance Runbook](#10-maintenance-runbook)

---

## 1. Deployment Options

| Option | Best For | Services Managed |
|--------|---------|-----------------|
| Docker Compose | Production / Demo | All 6 services (etcd, MinIO, Milvus, UI, API, setup) |
| Manual | Development / Debug | Each service started individually |
| HCLS AI Factory Stack | Full platform | Integrated with genomics, RAG, drug discovery pipelines |

### Port Assignments

| Port | Service | Protocol |
|------|---------|----------|
| 8507 | Streamlit UI | HTTP |
| 8107 | FastAPI REST API | HTTP |
| 19530 | Milvus gRPC | gRPC |
| 9091 | Milvus health/metrics | HTTP |

---

## 2. Docker Compose Deployment

### 2.1 Prerequisites

- Docker Engine 24.0+
- Docker Compose v2.20+
- 16GB+ RAM (32GB recommended for Milvus + embeddings)
- 20GB+ disk space
- Network access to Anthropic API

### 2.2 Quick Start

```bash
cd ai_agent_adds/pharmacogenomics_intelligence_agent

# Create environment file
cat > .env << 'EOF'
ANTHROPIC_API_KEY=sk-ant-your-key-here
PGX_MILVUS_HOST=milvus-standalone
PGX_MILVUS_PORT=19530
EOF

# Start all services
docker compose up -d

# Monitor setup progress
docker compose logs -f pgx-setup

# Verify all services are healthy
docker compose ps
```

### 2.3 Service Architecture

```
docker compose up -d
  |
  +-- milvus-etcd         (etcd v3.5.5)         Metadata store
  +-- milvus-minio        (MinIO)                Object storage
  +-- milvus-standalone   (Milvus v2.4)          Vector database
  |     depends_on: etcd, minio
  +-- pgx-streamlit       (Custom Python 3.12)   Streamlit UI (:8507)
  |     depends_on: milvus-standalone
  +-- pgx-api             (Custom Python 3.12)   FastAPI REST (:8107)
  |     depends_on: milvus-standalone
  +-- pgx-setup           (One-shot)             Create collections + seed
        depends_on: milvus-standalone
```

### 2.4 Health Checks

All services include Docker health checks:

| Service | Health Check | Interval | Retries |
|---------|-------------|----------|---------|
| milvus-etcd | `etcdctl endpoint health` | 30s | 5 |
| milvus-minio | `curl http://localhost:9000/minio/health/live` | 30s | 5 |
| milvus-standalone | `curl http://localhost:9091/healthz` | 30s | 10 |
| pgx-api | `curl http://localhost:8107/health` | 30s | 3 |

### 2.5 Volumes

Three named Docker volumes persist data:

| Volume | Mount | Purpose |
|--------|-------|---------|
| `etcd_data` | `/etcd` | Milvus metadata |
| `minio_data` | `/minio_data` | Milvus index/log objects |
| `milvus_data` | `/var/lib/milvus` | Milvus collection data |

### 2.6 Network

All services connect to the `pgx-network` bridge network. Inter-service communication uses container names as hostnames.

### 2.7 Stopping Services

```bash
# Graceful stop (preserves data)
docker compose down

# Stop and remove volumes (DESTRUCTIVE)
docker compose down -v
```

### 2.8 Environment Variable Reference

All environment variables use the `PGX_` prefix. The full list of configurable settings:

| Variable | Default | Description |
|----------|---------|-------------|
| `ANTHROPIC_API_KEY` | (required) | Anthropic API key for Claude Sonnet 4.6 |
| `PGX_MILVUS_HOST` | localhost | Milvus server hostname |
| `PGX_MILVUS_PORT` | 19530 | Milvus server port |
| `PGX_LLM_MODEL` | claude-sonnet-4-6 | LLM model identifier |
| `PGX_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model name |
| `PGX_EMBEDDING_DIMENSION` | 384 | Embedding vector dimension |
| `PGX_TOP_K_PER_COLLECTION` | 5 | Maximum results per collection |
| `PGX_SCORE_THRESHOLD` | 0.4 | Minimum similarity score |
| `PGX_API_PORT` | 8107 | FastAPI server port |
| `PGX_STREAMLIT_PORT` | 8507 | Streamlit UI port |
| `PGX_CORS_ORIGINS` | localhost:8080,8107,8507 | Allowed CORS origins |
| `PGX_MAX_REQUEST_SIZE_MB` | 10 | Maximum request body size |
| `PGX_INGEST_SCHEDULE_HOURS` | 168 | Ingest refresh interval (hours) |
| `PGX_INGEST_ENABLED` | false | Enable automated ingest scheduler |
| `PGX_MAX_CONVERSATION_CONTEXT` | 3 | Prior exchanges retained in context |
| `NCBI_API_KEY` | (optional) | NCBI API key for PubMed ingest rate limiting |

Collection search weights (all configurable, must sum to 1.0):

| Variable | Default | Collection |
|----------|---------|------------|
| `PGX_WEIGHT_DRUG_GUIDELINES` | 0.14 | pgx_drug_guidelines |
| `PGX_WEIGHT_DRUG_INTERACTIONS` | 0.12 | pgx_drug_interactions |
| `PGX_WEIGHT_GENE_REFERENCE` | 0.10 | pgx_gene_reference |
| `PGX_WEIGHT_HLA_HYPERSENSITIVITY` | 0.10 | pgx_hla_hypersensitivity |
| `PGX_WEIGHT_CLINICAL_EVIDENCE` | 0.08 | pgx_clinical_evidence |
| `PGX_WEIGHT_PHENOCONVERSION` | 0.08 | pgx_phenoconversion |
| `PGX_WEIGHT_DOSING_ALGORITHMS` | 0.07 | pgx_dosing_algorithms |
| `PGX_WEIGHT_FDA_LABELS` | 0.06 | pgx_fda_labels |
| `PGX_WEIGHT_POPULATION_DATA` | 0.06 | pgx_population_data |
| `PGX_WEIGHT_DRUG_ALTERNATIVES` | 0.05 | pgx_drug_alternatives |
| `PGX_WEIGHT_CLINICAL_TRIALS` | 0.04 | pgx_clinical_trials |
| `PGX_WEIGHT_GENOMIC_EVIDENCE` | 0.03 | genomic_evidence |
| `PGX_WEIGHT_PATIENT_PROFILES` | 0.03 | pgx_patient_profiles |
| `PGX_WEIGHT_IMPLEMENTATION` | 0.02 | pgx_implementation |
| `PGX_WEIGHT_EDUCATION` | 0.02 | pgx_education |

---

## 3. Manual Deployment

### 3.1 Python Environment

```bash
# Create virtual environment
python3.12 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 3.2 Milvus Setup

Ensure Milvus 2.4 is running on localhost:19530. For standalone:

```bash
# Using Docker for Milvus only
docker run -d --name milvus-standalone \
  -p 19530:19530 \
  -p 9091:9091 \
  -v milvus_data:/var/lib/milvus \
  milvusdb/milvus:v2.4-latest \
  milvus run standalone
```

### 3.3 Collection Setup and Seeding

```bash
# Create all 15 collections and seed 240 records
python scripts/setup_collections.py --drop-existing --seed

# Seed knowledge base references
python scripts/seed_knowledge.py
```

### 3.4 Start API Server

```bash
export ANTHROPIC_API_KEY=sk-ant-your-key-here
uvicorn api.main:app --host 0.0.0.0 --port 8107 --workers 2
```

### 3.5 Start Streamlit UI

```bash
export ANTHROPIC_API_KEY=sk-ant-your-key-here
streamlit run app/pgx_ui.py --server.port 8507 --server.address 0.0.0.0
```

---

## 4. Milvus Tuning

### 4.1 Index Configuration

Current settings optimized for the PGx collection sizes:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Index type | IVF_FLAT | Best quality for < 1M vectors per collection |
| Distance metric | COSINE | Standard for text embeddings |
| nlist | 1024 | Number of cluster centroids |
| nprobe | 16 | Clusters searched at query time (recall vs speed) |

### 4.2 Detailed Index Parameter Tuning

**nlist (Number of Inverted File Lists)**

The `nlist` parameter controls how many Voronoi cells the vector space is partitioned into. Each cell contains a subset of vectors, and at query time only `nprobe` cells are searched.

| Collection Size | Recommended nlist | Rationale |
|----------------|-------------------|-----------|
| < 1,000 vectors | 32-64 | Small collections need fewer partitions |
| 1,000-10,000 | 128-256 | Moderate partitioning for balanced search |
| 10,000-100,000 | 512-1024 | Current default; good for post-ingest sizes |
| 100,000-1,000,000 | 2048-4096 | Larger partitioning for scale |
| > 1,000,000 | Consider HNSW | IVF_FLAT recall degrades at this scale |

For the PGx agent with 240 seed records, `nlist=1024` is conservative. If you only use seed data, reducing to `nlist=64` will improve recall slightly. The default of 1024 is chosen for post-ingest scenarios where collections may grow to tens of thousands of records after PubMed and PharmGKB ingestion.

**nprobe (Number of Cells to Search)**

The `nprobe` parameter controls the recall-speed tradeoff at query time. Higher values improve recall at the cost of latency.

| nprobe | Approximate Recall | Latency Impact | Recommended For |
|--------|-------------------|----------------|----------------|
| 1 | ~10% | Fastest | Never (too low for clinical use) |
| 8 | ~60-70% | Very fast | Development/testing only |
| 16 | ~80-85% | Fast (current default) | Demo, standard queries |
| 32 | ~90-93% | Moderate | Production with quality emphasis |
| 64 | ~95-97% | Slower | Critical clinical queries |
| 128 | ~98-99% | Slowest | Exhaustive search scenarios |

To change nprobe at query time without rebuilding the index:

```python
search_params = {"metric_type": "COSINE", "params": {"nprobe": 32}}
collection.search(vectors, "embedding", search_params, limit=5)
```

**segment_row_limit**

This parameter controls the maximum number of rows per data segment in Milvus. Smaller segments enable faster index building but consume more metadata. Larger segments reduce overhead but slow down insertions and index rebuilds.

| segment_row_limit | Use Case |
|-------------------|----------|
| 1024 (default) | Small collections (< 10K records) |
| 4096 | Medium collections (10K-100K records) |
| 16384 | Large collections (100K-1M records) |
| 65536 | Very large collections (> 1M records) |

For the PGx agent's 15 collections with 240 seed records, the default of 1024 is appropriate. After full ingest, consider increasing to 4096 for collections that grow beyond 10,000 records (typically `pgx_clinical_evidence` and `pgx_clinical_trials` after PubMed and ClinicalTrials.gov ingestion).

### 4.3 Memory Configuration

For DGX Spark (128GB unified memory):

```yaml
# milvus.yaml overrides
queryNode:
  cacheSize: 8GB          # In-memory cache for loaded segments
dataNode:
  insertBufSize: 16MB     # Write buffer per collection
```

**Detailed memory allocation guidance:**

| Component | Memory Usage | Configuration |
|-----------|-------------|---------------|
| etcd metadata | 50-100MB | Minimal; scales with collection count |
| MinIO object storage | 100-500MB | Index files and logs; cached from disk |
| Milvus query node cache | 1-8GB | Set `cacheSize`; holds loaded collection segments |
| Milvus data node | 100-500MB | Set `insertBufSize`; buffers during ingest |
| BGE embedding model | ~500MB | Loaded into GPU or CPU memory |
| Python API/UI processes | 200-500MB each | FastAPI + Streamlit runtime |
| Operating system + overhead | 2-4GB | Kernel, Docker daemon, system services |

**Total recommended allocation:**

| Deployment Size | Total RAM | Milvus cacheSize | Notes |
|----------------|-----------|-------------------|-------|
| Development (8GB) | 8GB | 2GB | Tight; disable embedding model GPU offload |
| Demo (32GB) | 32GB | 8GB | Comfortable; all services with headroom |
| Production (64GB+) | 64GB+ | 16-24GB | Room for large post-ingest collections |
| DGX Spark (128GB) | 128GB | 32GB | Full capacity; GPU acceleration enabled |

### 4.4 Performance Tuning

| Scenario | Recommendation |
|----------|---------------|
| Slow first query | Pre-load collections: `collection.load()` during startup |
| High concurrent queries | Increase nprobe from 16 to 32 for better recall |
| Large collections (>100K vectors) | Consider HNSW index for better recall at scale |
| Memory pressure | Reduce nlist to 512, release unused collections |
| Ingest throughput | Increase `insertBufSize` to 64MB; batch insert with 1000-record batches |
| Index rebuild after large ingest | Call `collection.release()` then `collection.load()` to pick up new index |

### 4.5 HNSW Index Migration

For collections that grow beyond 100,000 vectors after full ingest, migrating from IVF_FLAT to HNSW provides better recall at comparable latency:

```python
# HNSW index parameters
index_params = {
    "index_type": "HNSW",
    "metric_type": "COSINE",
    "params": {
        "M": 16,          # Number of edges per node (higher = better recall, more memory)
        "efConstruction": 256  # Build-time search width (higher = better quality, slower build)
    }
}

# HNSW search parameters
search_params = {
    "metric_type": "COSINE",
    "params": {
        "ef": 64  # Query-time search width (higher = better recall, slower query)
    }
}
```

| Parameter | Recommended Range | Impact |
|-----------|------------------|--------|
| M | 8-64 (default: 16) | Memory per vector; 16 is a good balance |
| efConstruction | 128-512 (default: 256) | Index build quality; higher = slower build, better graph |
| ef (query time) | 32-256 (default: 64) | Recall vs speed tradeoff at query time |

**When to migrate to HNSW:**

- Collection exceeds 100K vectors
- IVF_FLAT recall at nprobe=32 drops below 90% for your query distribution
- You need consistent sub-100ms search latency per collection

**When to stay on IVF_FLAT:**

- Collections remain under 100K vectors
- You need the smallest possible memory footprint
- Index rebuild time must be minimal (IVF_FLAT rebuilds faster)

---

## 5. GPU Considerations

### 5.1 DGX Spark (GB10)

The NVIDIA DGX Spark's GB10 GPU accelerates:
- **Embedding generation**: BGE-small-en-v1.5 inference for ingest and query embedding
- **Milvus GPU index**: Optional GPU-accelerated IVF search (requires Milvus GPU build)

### 5.2 CPU-Only Deployment

The system runs fully on CPU. GPU is optional and only affects:
- Embedding generation speed (2-5x faster with GPU)
- Milvus index search (marginal improvement for current collection sizes)

### 5.3 Embedding Model on GPU

If using GPU for embeddings, ensure CUDA 12.x and torch are installed:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu121
```

### 5.4 GPU Memory Requirements

| Component | GPU Memory | Notes |
|-----------|-----------|-------|
| BGE-small-en-v1.5 model | ~300MB | Loaded once; shared across all queries |
| Embedding batch inference | 50-200MB | Depends on batch size; scales with `EMBEDDING_BATCH_SIZE` |
| Milvus GPU index (optional) | 500MB-4GB | Depends on collection sizes; only with Milvus GPU build |
| PyTorch runtime overhead | ~500MB | CUDA context, cuDNN, memory allocator |

**Total GPU memory for embedding inference:** ~1GB minimum, 2GB recommended.

On DGX Spark with 128GB unified memory (shared CPU/GPU), this is a negligible fraction of available resources. On discrete GPU systems, a GPU with at least 4GB VRAM is recommended.

### 5.5 GPU Acceleration for Ingest

During full ingest runs (PubMed, PharmGKB, ClinicalTrials.gov), GPU acceleration reduces embedding time significantly:

| Records | CPU Time (est.) | GPU Time (est.) | Speedup |
|---------|-----------------|-----------------|---------|
| 240 (seed) | 15 seconds | 5 seconds | 3x |
| 10,000 | 10 minutes | 2 minutes | 5x |
| 100,000 | 90 minutes | 20 minutes | 4.5x |

To maximize GPU throughput during ingest:

```bash
# Increase embedding batch size for GPU
export PGX_EMBEDDING_BATCH_SIZE=128  # Default is 32

# Run ingest
python scripts/setup_collections.py --drop-existing --seed
```

### 5.6 Multi-GPU Considerations

For systems with multiple GPUs, the embedding model can be pinned to a specific device:

```bash
export CUDA_VISIBLE_DEVICES=0  # Use only GPU 0 for embeddings
```

If running Milvus GPU index on the same system, assign different GPUs to avoid memory contention:

```bash
# Milvus GPU index on GPU 1
# Embedding model on GPU 0
```

---

## 6. Production Security Checklist

### 6.1 API Key Management

| Key | Purpose | Storage |
|-----|---------|---------|
| `ANTHROPIC_API_KEY` | Claude Sonnet 4.6 access | `.env` file (not committed to git) |
| `NCBI_API_KEY` | PubMed ingest rate limiting | `.env` file (optional) |

**Production recommendations:**

- **Never commit API keys to version control.** The `.env` file should be in `.gitignore`.
- **Use a secrets manager** for production deployments. Recommended options:
  - HashiCorp Vault (self-hosted)
  - AWS Secrets Manager (cloud)
  - Docker secrets (Docker Swarm)
  - Kubernetes Secrets (K8s deployments)

```bash
# Example: Using Docker secrets instead of .env
echo "sk-ant-your-key-here" | docker secret create anthropic_api_key -

# Reference in docker-compose.yml:
# services:
#   pgx-api:
#     secrets:
#       - anthropic_api_key
#     environment:
#       ANTHROPIC_API_KEY_FILE: /run/secrets/anthropic_api_key
```

- **Rotate API keys** on a regular schedule (monthly recommended).
- **Audit API key usage** via Anthropic's usage dashboard.
- **Use separate keys** for development, staging, and production environments.

### 6.2 Network Security

- **CORS**: Configured via `PGX_CORS_ORIGINS` (default: localhost:8080, 8107, 8507)
- **Milvus**: No authentication by default; restrict network access via firewall or Docker network
- **API**: No authentication by default; add API key middleware for production

### 6.3 Input Validation

- Milvus filter expressions are sanitized via `_SAFE_FILTER_RE` regex (alphanumeric + safe characters only)
- All API inputs validated via Pydantic models with field-level constraints
- Request size limited to 10MB via `PGX_MAX_REQUEST_SIZE_MB`

### 6.4 HTTPS / TLS Configuration

For production, deploy a reverse proxy with TLS termination in front of all HTTP services:

**Nginx example:**

```nginx
server {
    listen 443 ssl;
    server_name pgx.yourdomain.com;

    ssl_certificate /etc/ssl/certs/pgx.crt;
    ssl_certificate_key /etc/ssl/private/pgx.key;
    ssl_protocols TLSv1.2 TLSv1.3;
    ssl_ciphers HIGH:!aNULL:!MD5;

    # FastAPI backend
    location /api/ {
        proxy_pass http://localhost:8107/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Streamlit UI
    location / {
        proxy_pass http://localhost:8507/;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
    }
}

# Redirect HTTP to HTTPS
server {
    listen 80;
    server_name pgx.yourdomain.com;
    return 301 https://$host$request_uri;
}
```

**Traefik example (docker-compose integration):**

```yaml
services:
  traefik:
    image: traefik:v3.0
    command:
      - "--entrypoints.web.address=:80"
      - "--entrypoints.websecure.address=:443"
      - "--certificatesresolvers.letsencrypt.acme.tlschallenge=true"
      - "--certificatesresolvers.letsencrypt.acme.email=admin@yourdomain.com"
    ports:
      - "80:80"
      - "443:443"
    labels:
      - "traefik.http.routers.pgx-api.rule=Host(`pgx.yourdomain.com`) && PathPrefix(`/api`)"
      - "traefik.http.routers.pgx-api.tls.certresolver=letsencrypt"
```

### 6.5 API Authentication

For production deployments, add API key authentication middleware to the FastAPI application:

```python
# Example API key middleware (add to api/main.py)
from fastapi import Security, HTTPException
from fastapi.security import APIKeyHeader

api_key_header = APIKeyHeader(name="X-API-Key")

async def verify_api_key(api_key: str = Security(api_key_header)):
    if api_key != os.environ.get("PGX_API_KEY"):
        raise HTTPException(status_code=403, detail="Invalid API key")
    return api_key

# Apply to endpoints:
@app.get("/query", dependencies=[Security(verify_api_key)])
```

### 6.6 Rate Limiting

Implement rate limiting to prevent API abuse and control LLM API costs:

```python
# Example using slowapi
from slowapi import Limiter
from slowapi.util import get_remote_address

limiter = Limiter(key_func=get_remote_address)

@app.post("/query")
@limiter.limit("10/minute")  # 10 RAG queries per minute per IP
async def query(request: Request, body: QueryRequest):
    ...

@app.post("/search")
@limiter.limit("30/minute")  # 30 evidence-only searches per minute per IP
async def search(request: Request, body: SearchRequest):
    ...
```

Recommended rate limits:

| Endpoint | Limit | Rationale |
|----------|-------|-----------|
| `/query` (RAG) | 10/minute | LLM API cost control |
| `/search` (evidence only) | 30/minute | Lower cost, higher throughput |
| `/v1/pgx/*` (clinical) | 20/minute | Moderate computational cost |
| `/health` | 60/minute | Health checks should be responsive |
| `/metrics` | 12/minute | Prometheus scrape interval |

### 6.7 Milvus Security

By default, Milvus does not require authentication. For production:

```yaml
# milvus.yaml - enable authentication
common:
  security:
    authorizationEnabled: true

# Create admin user
# python -c "from pymilvus import utility; utility.create_credential('admin', 'secure_password')"

# Connect with credentials
# connections.connect(host='localhost', port='19530', user='admin', password='secure_password')
```

Additional Milvus security measures:

- **Restrict port 19530** to internal network only (Docker network or firewall rule)
- **Disable MinIO public access** (MinIO should only be accessible from the Milvus container)
- **Enable TLS for Milvus connections** in production (Milvus 2.4 supports TLS)

### 6.8 Full Production Security Checklist

- [ ] Store API keys in a secrets manager (not .env files)
- [ ] Enable HTTPS via reverse proxy (nginx, Traefik, or cloud load balancer)
- [ ] Restrict Milvus port (19530) to internal network only
- [ ] Restrict MinIO ports to internal network only
- [ ] Add API authentication (JWT or API key header)
- [ ] Enable Milvus authentication (username/password)
- [ ] Configure CORS to allow only trusted origins
- [ ] Set up log rotation for container logs (see Section 10.5)
- [ ] Implement rate limiting on public-facing endpoints
- [ ] Enable request size limits (default 10MB is appropriate)
- [ ] Validate all input via Pydantic models (already implemented)
- [ ] Sanitize Milvus filter expressions (already implemented via `_SAFE_FILTER_RE`)
- [ ] Disable debug logging in production (`LOG_LEVEL=INFO`)
- [ ] Set up network segmentation (API/UI on public network, Milvus/etcd/MinIO on private)
- [ ] Configure Docker container resource limits (memory, CPU)
- [ ] Implement audit logging for clinical decision support queries
- [ ] Review and comply with HIPAA/GDPR requirements if processing PHI/PII
- [ ] Enable Docker content trust for image signing
- [ ] Run containers as non-root user (already configured in Dockerfile)

---

## 7. Monitoring Setup

### 7.1 Prometheus Metrics

The API exposes 22 Prometheus metrics at `GET /metrics`:

```bash
# Scrape endpoint
curl http://localhost:8107/metrics
```

### 7.2 Prometheus Configuration

Add to `prometheus.yml`:

```yaml
global:
  scrape_interval: 15s
  evaluation_interval: 15s

scrape_configs:
  - job_name: 'pgx-agent'
    static_configs:
      - targets: ['localhost:8107']
    metrics_path: '/metrics'
    scrape_interval: 15s

  - job_name: 'milvus'
    static_configs:
      - targets: ['localhost:9091']
    metrics_path: '/metrics'
    scrape_interval: 30s
```

### 7.3 Docker Compose Integration with Prometheus and Grafana

Add monitoring services to your Docker Compose stack:

```yaml
services:
  prometheus:
    image: prom/prometheus:v2.51.0
    ports:
      - "9090:9090"
    volumes:
      - ./monitoring/prometheus.yml:/etc/prometheus/prometheus.yml
      - prometheus_data:/prometheus
    networks:
      - pgx-network
    restart: unless-stopped

  grafana:
    image: grafana/grafana:10.4.0
    ports:
      - "3000:3000"
    environment:
      - GF_SECURITY_ADMIN_PASSWORD=pgx_admin_2026
      - GF_USERS_ALLOW_SIGN_UP=false
    volumes:
      - grafana_data:/var/lib/grafana
      - ./monitoring/grafana/dashboards:/var/lib/grafana/dashboards
      - ./monitoring/grafana/provisioning:/etc/grafana/provisioning
    networks:
      - pgx-network
    restart: unless-stopped

volumes:
  prometheus_data:
  grafana_data:
```

### 7.4 Key Metrics to Monitor

| Metric | Alert Threshold | Action |
|--------|----------------|--------|
| `pgx_query_latency_seconds` | p95 > 10s | Check Milvus load, LLM API latency |
| `pgx_errors_total` | > 10/min | Check API logs |
| `pgx_collections_connected` | < 15 | Milvus collection health issue |
| `pgx_llm_api_latency_seconds` | p95 > 30s | Anthropic API throttling |
| `pgx_alerts_generated_total{severity="CONTRAINDICATED"}` | Any | Review clinical alert for accuracy |
| `pgx_total_vectors` | Drops to 0 | Milvus data loss; trigger recovery |
| `pgx_embedding_latency_seconds` | p95 > 2s | GPU/CPU overload for embedding model |
| `pgx_phenoconversion_detections_total` | Rising trend | May indicate medication list quality issue |

### 7.5 Prometheus Alerting Rules

```yaml
# monitoring/prometheus_rules.yml
groups:
  - name: pgx_alerts
    rules:
      - alert: PGxHighQueryLatency
        expr: histogram_quantile(0.95, rate(pgx_query_latency_seconds_bucket[5m])) > 10
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "PGx query latency p95 > 10s"
          description: "Check Milvus connectivity and LLM API responsiveness."

      - alert: PGxHighErrorRate
        expr: rate(pgx_errors_total[5m]) > 0.1
        for: 2m
        labels:
          severity: critical
        annotations:
          summary: "PGx error rate > 6/minute"
          description: "Check API logs for error details."

      - alert: PGxCollectionsDisconnected
        expr: pgx_collections_connected < 15
        for: 1m
        labels:
          severity: critical
        annotations:
          summary: "PGx collections disconnected ({{ $value }}/15)"
          description: "Milvus may be down or collections may need to be reloaded."

      - alert: PGxMilvusDown
        expr: pgx_total_vectors == 0
        for: 30s
        labels:
          severity: critical
        annotations:
          summary: "PGx Milvus reports 0 vectors"
          description: "Milvus data may be lost. Check volumes and trigger recovery."

      - alert: PGxLLMLatencyHigh
        expr: histogram_quantile(0.95, rate(pgx_llm_api_latency_seconds_bucket[5m])) > 30
        for: 5m
        labels:
          severity: warning
        annotations:
          summary: "Anthropic API latency p95 > 30s"
          description: "API may be rate-limited. Check usage dashboard."

      - alert: PGxContraindicatedAlert
        expr: increase(pgx_alerts_generated_total{severity="CONTRAINDICATED"}[1h]) > 0
        labels:
          severity: info
        annotations:
          summary: "Contraindicated drug-gene interaction detected"
          description: "Review the clinical alert for accuracy and follow up."
```

### 7.6 Grafana Dashboard

Import a Grafana dashboard with the following panels:

**Row 1: Overview**
- Query volume (queries/min, counter)
- Active collections (gauge, target 15)
- Total vectors (gauge)
- Error rate (queries vs errors)

**Row 2: Latency**
- Query latency distribution (histogram, p50/p95/p99)
- LLM API latency (histogram)
- Embedding latency (histogram)
- Clinical pipeline stage latency (histogram by stage)

**Row 3: Clinical**
- Clinical alerts by severity (stacked bar: CONTRAINDICATED, MAJOR, MODERATE, MINOR, INFORMATIONAL)
- Drug checks count (counter)
- HLA screens by result (pie: CONTRAINDICATED, HIGH_RISK, SAFE, UNKNOWN)
- Dosing calculations by algorithm (bar: warfarin, tacrolimus, fluoropyrimidine, thiopurine)

**Row 4: Collections**
- Evidence count per query (histogram)
- Collection hit distribution (heatmap: which collections contribute to answers)
- Phenoconversion detections (counter)
- Export operations by format (bar: Markdown, JSON, PDF, FHIR)

**Row 5: Infrastructure**
- Milvus memory usage (from Milvus metrics endpoint)
- Container CPU and memory (from Docker/cAdvisor)
- LLM token usage (estimated from response length)
- Ingest freshness (time since last ingest run)

### 7.7 Log Aggregation

For centralized log management, forward Docker container logs to a log aggregation system:

```yaml
# docker-compose.yml logging configuration
services:
  pgx-api:
    logging:
      driver: "json-file"
      options:
        max-size: "50m"
        max-file: "5"
        tag: "pgx-api"
```

For ELK stack integration:

```yaml
services:
  pgx-api:
    logging:
      driver: "fluentd"
      options:
        fluentd-address: "localhost:24224"
        tag: "pgx.api"
```

---

## 8. Backup and Disaster Recovery

### 8.1 Backup Strategy Overview

| Data Type | Backup Method | Frequency | Recovery Time |
|-----------|--------------|-----------|---------------|
| Milvus vectors + indexes | Volume snapshot | Daily | 5-30 minutes |
| Seed data (JSON files) | Version control | Every change | < 1 minute |
| Configuration (.env, settings) | Version control + secrets backup | Every change | < 1 minute |
| Prometheus metrics | Prometheus TSDB snapshot | Weekly | 10 minutes |
| Grafana dashboards | JSON export | Every change | 5 minutes |

### 8.2 Milvus Data Backup

```bash
# Option 1: Docker volume backup (cold)
docker compose stop milvus-standalone
docker run --rm -v pgx_milvus_data:/data -v $(pwd)/backups:/backup \
  alpine tar czf /backup/milvus_data_$(date +%Y%m%d).tar.gz -C /data .
docker run --rm -v pgx_etcd_data:/data -v $(pwd)/backups:/backup \
  alpine tar czf /backup/etcd_data_$(date +%Y%m%d).tar.gz -C /data .
docker run --rm -v pgx_minio_data:/data -v $(pwd)/backups:/backup \
  alpine tar czf /backup/minio_data_$(date +%Y%m%d).tar.gz -C /data .
docker compose start milvus-standalone

# Option 2: Milvus backup utility (hot - no downtime)
# Requires milvus-backup tool: https://github.com/zilliztech/milvus-backup
milvus-backup create \
  --milvus.address localhost:19530 \
  --minio.address localhost:9000 \
  --name pgx_backup_$(date +%Y%m%d)

# Option 3: Logical backup via seed data (always available)
# The seed data in data/reference/ (14 JSON files) can regenerate all collections.
# This is the simplest recovery method and should always be your baseline.
```

### 8.3 Automated Backup Script

```bash
#!/bin/bash
# scripts/backup.sh - Run daily via cron
BACKUP_DIR="/backups/pgx-agent"
DATE=$(date +%Y%m%d_%H%M%S)
RETENTION_DAYS=30

mkdir -p "$BACKUP_DIR"

echo "[$(date)] Starting PGx Agent backup..."

# 1. Backup Milvus data (requires brief stop)
docker compose stop milvus-standalone
for vol in pgx_milvus_data pgx_etcd_data pgx_minio_data; do
    docker run --rm -v ${vol}:/data -v ${BACKUP_DIR}:/backup \
        alpine tar czf /backup/${vol}_${DATE}.tar.gz -C /data .
    echo "[$(date)] Backed up volume: ${vol}"
done
docker compose start milvus-standalone

# 2. Backup configuration
tar czf ${BACKUP_DIR}/config_${DATE}.tar.gz \
    .env config/ data/reference/ docker-compose.yml requirements.txt

# 3. Verify backup integrity
for f in ${BACKUP_DIR}/*_${DATE}.tar.gz; do
    if tar tzf "$f" > /dev/null 2>&1; then
        echo "[$(date)] Verified: $f ($(du -h $f | cut -f1))"
    else
        echo "[$(date)] ERROR: Backup file corrupted: $f"
    fi
done

# 4. Cleanup old backups
find ${BACKUP_DIR} -name "*.tar.gz" -mtime +${RETENTION_DAYS} -delete
echo "[$(date)] Cleaned up backups older than ${RETENTION_DAYS} days"

echo "[$(date)] Backup complete."
```

Add to crontab:

```bash
# Run backup daily at 2:00 AM
0 2 * * * /path/to/pharmacogenomics_intelligence_agent/scripts/backup.sh >> /var/log/pgx-backup.log 2>&1
```

### 8.4 Seed Data Backup

The seed data in `data/reference/` (14 JSON files) is the source of truth. Keep this directory in version control. From seed data alone, you can fully reconstruct all 15 Milvus collections and 240 records.

### 8.5 Recovery Procedures

**Scenario 1: Complete rebuild from seed data (simplest)**

```bash
# 1. Stop and remove all data
docker compose down -v

# 2. Start fresh
docker compose up -d

# 3. Wait for setup to complete
docker compose logs -f pgx-setup  # Wait for "PGx Setup complete!"

# 4. Verify
curl http://localhost:8107/health
```

Recovery time: ~2 minutes. This is the recommended recovery method for most scenarios.

**Scenario 2: Restore from volume backup**

```bash
# 1. Stop all services
docker compose down

# 2. Restore volumes
for vol in pgx_milvus_data pgx_etcd_data pgx_minio_data; do
    docker volume create ${vol}
    docker run --rm -v ${vol}:/data -v /backups/pgx-agent:/backup \
        alpine sh -c "cd /data && tar xzf /backup/${vol}_YYYYMMDD_HHMMSS.tar.gz"
done

# 3. Start services
docker compose up -d

# 4. Verify
curl http://localhost:8107/health
```

Recovery time: 5-30 minutes depending on data volume size.

**Scenario 3: Partial collection recovery**

If a specific collection is corrupted or missing:

```bash
# Drop and re-seed a single collection
python -c "
from pymilvus import connections, utility
connections.connect(host='localhost', port='19530')
if utility.has_collection('pgx_drug_guidelines'):
    utility.drop_collection('pgx_drug_guidelines')
"

# Re-run setup (will only create missing collections)
python scripts/setup_collections.py --seed
```

### 8.6 Disaster Recovery Plan

| Scenario | RTO | RPO | Recovery Method |
|----------|-----|-----|----------------|
| Single collection corruption | 5 min | 0 (seed data) | Re-seed from JSON |
| Milvus crash (data intact) | 2 min | 0 | Restart container |
| Milvus data loss (volumes lost) | 5 min | 0 (seed data) | Rebuild from seed |
| Full system failure (VM/hardware) | 30 min | Last backup | Restore from volume backup on new host |
| Configuration loss | 1 min | Last commit | Restore from git |

RTO = Recovery Time Objective. RPO = Recovery Point Objective.

---

## 9. Scaling Strategies

### 9.1 Vertical Scaling

| Component | Scale Up | Impact |
|-----------|----------|--------|
| uvicorn workers | `--workers 4` (from 2) | 2x API throughput |
| Milvus cache | Increase `cacheSize` | Faster repeated queries |
| Embedding batch | Increase `EMBEDDING_BATCH_SIZE` from 32 to 128 | Faster ingest |
| nprobe | Increase from 16 to 32 | Better recall, slightly slower queries |
| Thread pool | Increase ThreadPoolExecutor max_workers | More parallel collection searches |

### 9.2 Horizontal Scaling

| Component | Strategy | Configuration |
|-----------|---------|---------------|
| API | Multiple uvicorn containers behind load balancer | Nginx upstream or Traefik service discovery |
| Milvus | Cluster mode with dedicated query/data/index nodes | See Milvus cluster deployment guide |
| Embedding | Dedicated GPU microservice with gRPC interface | Triton Inference Server or custom FastAPI |
| LLM | Request queue with rate-limited Anthropic API proxy | Redis queue + worker pool |

**Milvus Cluster Mode:**

For deployments serving more than 100 concurrent users or managing more than 10M vectors:

```yaml
# Milvus cluster components
services:
  milvus-rootcoord:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "rootcoord"]

  milvus-querycoord:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "querycoord"]

  milvus-querynode-1:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "querynode"]

  milvus-querynode-2:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "querynode"]

  milvus-datacoord:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "datacoord"]

  milvus-datanode:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "datanode"]

  milvus-indexcoord:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "indexcoord"]

  milvus-indexnode:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "indexnode"]

  milvus-proxy:
    image: milvusdb/milvus:v2.4-latest
    command: ["milvus", "run", "proxy"]
    ports:
      - "19530:19530"
```

**API Horizontal Scaling with Nginx:**

```nginx
upstream pgx_api {
    least_conn;
    server pgx-api-1:8107;
    server pgx-api-2:8107;
    server pgx-api-3:8107;
}

server {
    listen 8107;
    location / {
        proxy_pass http://pgx_api;
    }
}
```

### 9.3 Resource Requirements

| Deployment | RAM | CPU | Disk | GPU |
|-----------|-----|-----|------|-----|
| Development | 8GB | 4 cores | 10GB | None |
| Demo (DGX Spark) | 32GB | 10 cores | 20GB | GB10 (optional) |
| Production (single-node) | 64GB+ | 16+ cores | 100GB+ | Recommended |
| Production (clustered) | 128GB+ (distributed) | 32+ cores (distributed) | 500GB+ (distributed) | Required |

### 9.4 Load Testing

Before production deployment, conduct load testing to validate performance targets:

```bash
# Install load testing tool
pip install locust

# Create locustfile.py
cat > locustfile.py << 'PYEOF'
from locust import HttpUser, task, between

class PGxUser(HttpUser):
    wait_time = between(1, 5)

    @task(3)
    def health_check(self):
        self.client.get("/health")

    @task(1)
    def search_query(self):
        self.client.post("/search", json={
            "query": "CYP2D6 codeine poor metabolizer",
            "top_k": 10
        })

    @task(1)
    def drug_check(self):
        self.client.post("/v1/pgx/drug-check", json={
            "drug": "codeine",
            "gene": "CYP2D6",
            "phenotype": "poor_metabolizer"
        })
PYEOF

# Run load test
locust -f locustfile.py --host http://localhost:8107 --users 20 --spawn-rate 2
```

---

## 10. Maintenance Runbook

### 10.1 Daily Tasks

- Verify `/health` endpoint returns `status: healthy` and `collections: 15`
- Check Docker container status: `docker compose ps` (all services should show "Up" or "Exited (0)" for pgx-setup)
- Review error logs: `docker compose logs --tail 50 pgx-api | grep ERROR`

### 10.2 Weekly Tasks

- Check `/health` endpoint for collection status and vector counts
- Review Prometheus metrics for anomalies (query latency trends, error rate)
- Verify automated ingest scheduler is running (if enabled): check `pgx_last_ingest_timestamp` metric
- Review disk space usage: `docker system df`

### 10.3 Monthly Tasks

- Backup Milvus data volumes (see Section 8)
- Update seed data with new CPIC guideline releases
- Review and rotate API keys
- Check for updated Docker images: `docker compose pull`
- Review and apply security patches to base images
- Test disaster recovery procedure (restore from backup on a test environment)

### 10.4 Updating Seed Data

```bash
# 1. Update JSON files in data/reference/
# 2. Re-seed collections
python scripts/setup_collections.py --drop-existing --seed
python scripts/seed_knowledge.py

# Or via Docker:
docker compose run --rm pgx-setup
```

### 10.5 Log Rotation

Configure Docker log rotation to prevent disk exhaustion:

```json
// /etc/docker/daemon.json
{
  "log-driver": "json-file",
  "log-opts": {
    "max-size": "50m",
    "max-file": "5"
  }
}
```

After updating, restart Docker daemon:

```bash
sudo systemctl restart docker
```

For individual service override in docker-compose.yml:

```yaml
services:
  pgx-api:
    logging:
      driver: "json-file"
      options:
        max-size: "100m"
        max-file: "10"
```

### 10.6 Collection Compaction

Over time, Milvus collections accumulate deleted records and fragmented segments. Compact collections periodically to reclaim space and improve performance:

```python
from pymilvus import connections, Collection

connections.connect(host="localhost", port="19530")

collections = [
    "pgx_gene_reference", "pgx_drug_guidelines", "pgx_drug_interactions",
    "pgx_hla_hypersensitivity", "pgx_phenoconversion", "pgx_dosing_algorithms",
    "pgx_clinical_evidence", "pgx_population_data", "pgx_clinical_trials",
    "pgx_fda_labels", "pgx_drug_alternatives", "pgx_patient_profiles",
    "pgx_implementation", "pgx_education", "genomic_evidence"
]

for name in collections:
    col = Collection(name)
    col.compact()
    print(f"Compacted: {name}")
```

Run compaction monthly, or after large ingest runs that involve deletions and re-insertions.

### 10.7 Index Rebuild

After significant data changes (large ingest, seed data update), rebuild indexes to ensure optimal search performance:

```python
from pymilvus import connections, Collection

connections.connect(host="localhost", port="19530")

for name in collections:
    col = Collection(name)
    col.release()  # Release from memory
    col.drop_index()  # Drop existing index
    col.create_index(
        field_name="embedding",
        index_params={
            "index_type": "IVF_FLAT",
            "metric_type": "COSINE",
            "params": {"nlist": 1024}
        }
    )
    col.load()  # Reload into memory
    print(f"Rebuilt index: {name}")
```

### 10.8 Updating the Application

```bash
# 1. Pull latest code
git pull

# 2. Rebuild containers
docker compose build

# 3. Restart with new images
docker compose up -d

# 4. Verify health
curl http://localhost:8107/health
```

### 10.9 Troubleshooting Quick Reference

| Symptom | Probable Cause | Resolution |
|---------|---------------|------------|
| `/health` returns 0 collections | Milvus down or collections not created | `docker compose restart milvus-standalone` then re-run setup |
| Query latency > 10s | LLM API throttling or Milvus overloaded | Check `pgx_llm_api_latency_seconds`; increase nprobe or reduce TOP_K |
| "Connection refused" on 19530 | Milvus not running | `docker compose up -d milvus-standalone` |
| Embedding generation slow | CPU-only mode or GPU not detected | Verify `torch.cuda.is_available()`; install CUDA toolkit |
| Out of memory | Milvus cache too large | Reduce `cacheSize` in milvus.yaml; release unused collections |
| Stale search results after ingest | Index not rebuilt | Run `collection.release()` then `collection.load()` |
| Container restart loop | Missing environment variable or port conflict | Check `docker compose logs <service>` for error message |
| CORS errors in browser | Origin not in PGX_CORS_ORIGINS | Add the origin URL to the CORS configuration |
