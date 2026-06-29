# Neurology Intelligence Agent -- Deployment Guide

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## 1. Prerequisites

### Hardware Requirements

| Component | Minimum | Recommended |
|---|---|---|
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16+ GB |
| Disk | 20 GB free | 50+ GB SSD |
| GPU | Not required | NVIDIA GPU (for accelerated embedding) |

### Software Requirements

| Software | Version | Purpose |
|---|---|---|
| Docker | 24+ | Container runtime |
| Docker Compose | 2.x | Service orchestration |
| Python | 3.10+ | Local development |
| Git | 2.x | Source control |

### API Keys

| Key | Source | Required |
|---|---|---|
| `ANTHROPIC_API_KEY` | [console.anthropic.com](https://console.anthropic.com) | Yes (for LLM synthesis) |
| `NCBI_API_KEY` | [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/account/) | Optional (PubMed ingest) |
| `CLINICALTRIALS_API_KEY` | [clinicaltrials.gov](https://clinicaltrials.gov) | Optional (trial ingest) |

---

## 2. Standalone Docker Deployment

### 2.1 Clone and Configure

```bash
cd /path/to/hcls-ai-factory/ai_agent_adds/neurology_intelligence_agent

# Create environment file
cp .env.example .env

# Edit .env and set your Anthropic API key
# ANTHROPIC_API_KEY=sk-ant-...
```

### 2.2 Build and Start

```bash
# Build all images
docker compose build

# Start all services (detached)
docker compose up -d

# Watch setup progress (collection creation + seeding)
docker compose logs -f neuro-setup
```

### 2.3 Verify Deployment

```bash
# Check service health
curl http://localhost:8528/health

# Expected response:
# {
#   "status": "healthy",
#   "agent": "neurology-intelligence-agent",
#   "version": "1.0.0",
#   "components": {
#     "milvus": "connected",
#     "rag_engine": "ready",
#     "workflow_engine": "ready"
#   },
#   "collections": 14,
#   "total_vectors": ...,
#   "workflows": 9,
#   "scales": 10
# }

# List collections
curl http://localhost:8528/collections

# List available workflows
curl http://localhost:8528/workflows
```

### 2.4 Access the Application

| Interface | URL |
|---|---|
| REST API | http://localhost:8528 |
| API Documentation (Swagger) | http://localhost:8528/docs |
| Streamlit Chat UI | http://localhost:8529 |
| Milvus Health | http://localhost:59091/healthz |
| Prometheus Metrics | http://localhost:8528/metrics |

---

## 3. Integrated Deployment (HCLS AI Factory)

When deploying as part of the full HCLS AI Factory stack on DGX Spark, the agent uses the shared Milvus instance.

### 3.1 Configuration for Shared Milvus

Set the following environment variables to point to the shared Milvus:

```bash
NEURO_MILVUS_HOST=milvus-standalone   # Shared Milvus service name
NEURO_MILVUS_PORT=19530               # Default Milvus port
```

### 3.2 Add to Top-Level Docker Compose

Add the neurology agent services to `docker-compose.dgx-spark.yml`:

```yaml
neuro-api:
  build:
    context: ./ai_agent_adds/neurology_intelligence_agent
    dockerfile: Dockerfile
  container_name: neuro-api
  ports:
    - "8528:8528"
  environment:
    NEURO_MILVUS_HOST: milvus-standalone
    NEURO_MILVUS_PORT: "19530"
    ANTHROPIC_API_KEY: ${ANTHROPIC_API_KEY}
  command:
    - uvicorn
    - api.main:app
    - --host=0.0.0.0
    - --port=8528
    - --workers=2
  depends_on:
    milvus-standalone:
      condition: service_healthy
  networks:
    - hcls-network
  restart: unless-stopped

neuro-streamlit:
  build:
    context: ./ai_agent_adds/neurology_intelligence_agent
    dockerfile: Dockerfile
  container_name: neuro-streamlit
  ports:
    - "8529:8529"
  environment:
    NEURO_MILVUS_HOST: milvus-standalone
    NEURO_MILVUS_PORT: "19530"
    NEURO_API_BASE: "http://neuro-api:8528"
    ANTHROPIC_API_KEY: ${ANTHROPIC_API_KEY}
  depends_on:
    milvus-standalone:
      condition: service_healthy
  networks:
    - hcls-network
  restart: unless-stopped
```

### 3.3 Initialize Collections

After the shared Milvus is running:

```bash
# Run setup inside the API container
docker exec neuro-api python scripts/setup_collections.py --drop-existing --seed
docker exec neuro-api python scripts/seed_knowledge.py
```

---

## 4. Local Development

### 4.1 Python Environment

```bash
cd neurology_intelligence_agent

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Set environment variables
export ANTHROPIC_API_KEY=sk-ant-...
export NEURO_MILVUS_HOST=localhost
export NEURO_MILVUS_PORT=19530
```

### 4.2 Start Milvus Locally

```bash
# Start only the Milvus stack
docker compose up -d milvus-etcd milvus-minio milvus-standalone
```

### 4.3 Initialize Data

```bash
# Create collections
python scripts/setup_collections.py --drop-existing --seed

# Seed knowledge base
python scripts/seed_knowledge.py
```

### 4.4 Start the API Server

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8528 --reload
```

### 4.5 Start the Streamlit UI

```bash
# In a separate terminal
streamlit run app/neuro_ui.py --server.port 8529
```

### 4.6 Run Tests

```bash
pytest tests/ -v

# Run specific test module
pytest tests/test_clinical_scales.py -v

# Run with coverage
pytest tests/ --cov=src --cov-report=term-missing
```

---

## 5. Configuration Reference

All settings use the `NEURO_` environment variable prefix:

| Variable | Default | Description |
|---|---|---|
| `NEURO_MILVUS_HOST` | localhost | Milvus server hostname |
| `NEURO_MILVUS_PORT` | 19530 | Milvus server port |
| `NEURO_API_HOST` | 0.0.0.0 | API bind address |
| `NEURO_API_PORT` | 8528 | API listen port |
| `NEURO_STREAMLIT_PORT` | 8529 | Streamlit listen port |
| `NEURO_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model |
| `NEURO_LLM_MODEL` | claude-sonnet-4-6 | LLM model |
| `NEURO_SCORE_THRESHOLD` | 0.4 | Minimum cosine similarity |
| `NEURO_API_KEY` | (empty) | API key for auth (empty = no auth) |
| `NEURO_CORS_ORIGINS` | localhost:8080,8528,8529 | CORS allowlist |
| `NEURO_INGEST_ENABLED` | false | Enable scheduled ingestion |
| `NEURO_INGEST_SCHEDULE_HOURS` | 24 | Ingestion interval |
| `NEURO_MAX_REQUEST_SIZE_MB` | 10 | Max request body size |
| `NEURO_CROSS_AGENT_TIMEOUT` | 30 | Cross-agent request timeout (seconds) |

---

## 6. Service Management

### Start/Stop

```bash
# Start all services
docker compose up -d

# Stop all services
docker compose down

# Stop and remove volumes (WARNING: deletes all data)
docker compose down -v
```

### Logs

```bash
# All services
docker compose logs -f

# Specific service
docker compose logs -f neuro-api

# Last 100 lines
docker compose logs --tail=100 neuro-api
```

### Health Checks

```bash
# API health
curl -s http://localhost:8528/health | python -m json.tool

# Milvus health
curl -s http://localhost:59091/healthz

# Container status
docker compose ps
```

### Restart a Service

```bash
docker compose restart neuro-api
docker compose restart neuro-streamlit
```

---

## 7. Troubleshooting

| Issue | Cause | Solution |
|---|---|---|
| `"milvus": "unavailable"` | Milvus not started or not healthy | Check `docker compose ps`, wait for health check |
| `"rag_engine": "unavailable"` | Embedding model failed to load | Check disk space, verify `sentence-transformers` installed |
| 401 Unauthorized | API key mismatch | Verify `NEURO_API_KEY` matches `X-API-Key` header |
| 429 Rate Limited | Too many requests | Wait 60 seconds, or increase `_RATE_LIMIT_MAX` |
| Empty search results | Collections not seeded | Run `scripts/setup_collections.py --seed` |
| LLM timeout | Anthropic API unreachable | Check network, verify `ANTHROPIC_API_KEY` |
| Port conflict (8528) | Another service on port | Change `NEURO_API_PORT` or stop conflicting service |

---

## 8. Production Hardening

For production deployment beyond demo/development:

1. **TLS Termination:** Place nginx or Traefik reverse proxy in front of port 8528
2. **Milvus Cluster:** Replace standalone Milvus with distributed cluster for HA
3. **Redis Rate Limiting:** Replace in-memory rate limiter with Redis-backed solution
4. **Persistent Logging:** Configure log aggregation (ELK, Loki)
5. **Backup Strategy:** Schedule Milvus collection snapshots
6. **Monitoring:** Connect Prometheus scraping to Grafana dashboards
7. **Secrets Management:** Use Vault or AWS Secrets Manager for API keys

---

*Neurology Intelligence Agent -- Deployment Guide v1.3.0*
*HCLS AI Factory / GTC Europe 2026*
