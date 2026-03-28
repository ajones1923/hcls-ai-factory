# Single-Cell Intelligence Agent -- Deployment Guide

**Date:** 2026-03-22
**Author:** Adam Jones

---

## 1. Prerequisites

### 1.1 System Requirements

| Requirement | Minimum | Recommended |
|------------|---------|-------------|
| CPU | 4 cores | 8+ cores |
| RAM | 8 GB | 16 GB |
| Disk | 10 GB | 50 GB |
| GPU | None (CPU mode) | NVIDIA GPU with 8+ GB VRAM |
| Docker | 24.0+ | 25.0+ |
| Docker Compose | v2.20+ | v2.24+ |
| Python | 3.10+ | 3.10.x |
| Network | 100 Mbps | 1 Gbps |

### 1.2 Port Requirements

| Port | Service | Required |
|------|---------|----------|
| 8540 | FastAPI REST API | Yes |
| 8130 | Streamlit UI | Yes |
| 19530 | Milvus (shared mode) | Yes (integrated) |
| 69530 | Milvus (standalone mode) | Yes (standalone) |
| 69091 | Milvus health (standalone) | Optional |

### 1.3 API Keys

| Key | Required | Source |
|-----|----------|--------|
| ANTHROPIC_API_KEY | Yes (for LLM features) | https://console.anthropic.com |
| SC_API_KEY | Optional | Self-generated for API auth |
| SC_NCBI_API_KEY | Optional | https://www.ncbi.nlm.nih.gov/account/ |

---

## 2. Standalone Deployment (Docker Compose)

### 2.1 Clone and Configure

```bash
cd /path/to/hcls-ai-factory/ai_agent_adds/single_cell_intelligence_agent

# Create environment file
cp .env.example .env

# Edit .env and set your API key
echo "ANTHROPIC_API_KEY=sk-ant-your-key-here" >> .env
```

### 2.2 Start All Services

```bash
# Start Milvus infrastructure + agent services
docker compose up -d

# Watch the setup script progress
docker compose logs -f sc-setup

# Verify all services are healthy
docker compose ps
```

Expected output:

```
NAME                 STATUS          PORTS
sc-milvus-etcd       Up (healthy)
sc-milvus-minio      Up (healthy)
sc-milvus-standalone Up (healthy)    0.0.0.0:69530->19530
sc-api               Up (healthy)    0.0.0.0:8540->8540
sc-streamlit         Up              0.0.0.0:8130->8130
sc-setup             Exited (0)
```

### 2.3 Verify Deployment

```bash
# Health check
curl http://localhost:8540/health | python -m json.tool

# Expected response:
# {
#   "status": "healthy",
#   "agent": "single-cell-intelligence-agent",
#   "version": "1.0.0",
#   "components": {
#     "milvus": "connected",
#     "rag_engine": "ready",
#     "workflow_engine": "ready"
#   },
#   "collections": 12,
#   "total_vectors": 144,
#   "workflows": 10
# }

# List collections
curl http://localhost:8540/collections | python -m json.tool

# List workflows
curl http://localhost:8540/workflows | python -m json.tool
```

### 2.4 Access the UI

Open your browser to `http://localhost:8130` to access the Streamlit interface.

---

## 3. Integrated Deployment (DGX Spark)

### 3.1 Shared Milvus Configuration

When deployed alongside other HCLS AI Factory agents, the single-cell agent connects to the shared Milvus instance:

```bash
# Set environment variables for shared Milvus
export SC_MILVUS_HOST=milvus-standalone
export SC_MILVUS_PORT=19530
```

### 3.2 Docker Compose Override

Use the top-level `docker-compose.dgx-spark.yml` which includes the single-cell agent alongside all other services:

```bash
cd /path/to/hcls-ai-factory
docker compose -f docker-compose.dgx-spark.yml up -d sc-api sc-streamlit
```

### 3.3 Collection Setup

Run the setup script to create collections in the shared Milvus instance:

```bash
docker compose -f docker-compose.dgx-spark.yml run --rm sc-setup
```

---

## 4. Manual Installation (Development)

### 4.1 Python Environment

```bash
cd single_cell_intelligence_agent

# Create virtual environment
python3.10 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 4.2 Configure Environment

```bash
# Copy and edit environment file
cp .env.example .env
# Set ANTHROPIC_API_KEY in .env

# Or export directly
export ANTHROPIC_API_KEY=sk-ant-your-key-here
export SC_MILVUS_HOST=localhost
export SC_MILVUS_PORT=19530
```

### 4.3 Start Milvus (if not running)

```bash
# Start just the Milvus stack
docker compose up -d milvus-etcd milvus-minio milvus-standalone
```

### 4.4 Create Collections and Seed Data

```bash
python scripts/setup_collections.py --drop-existing --seed
python scripts/seed_knowledge.py
```

### 4.5 Start API Server

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8540 --reload
```

### 4.6 Start Streamlit UI

```bash
streamlit run app/sc_ui.py --server.port 8130 --server.address 0.0.0.0
```

---

## 5. Configuration Reference

### 5.1 Environment Variables

All configuration is via environment variables with `SC_` prefix:

| Variable | Default | Description |
|----------|---------|-------------|
| `SC_MILVUS_HOST` | localhost | Milvus server hostname |
| `SC_MILVUS_PORT` | 19530 | Milvus server port |
| `SC_API_PORT` | 8540 | FastAPI server port |
| `SC_STREAMLIT_PORT` | 8130 | Streamlit UI port |
| `SC_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Sentence-transformer model |
| `SC_LLM_MODEL` | claude-sonnet-4-6 | Anthropic model name |
| `SC_SCORE_THRESHOLD` | 0.4 | Minimum search similarity |
| `SC_API_KEY` | (empty) | API authentication key |
| `SC_MAX_REQUEST_SIZE_MB` | 10 | Max request body size |
| `SC_GPU_MEMORY_LIMIT_GB` | 120 | GPU memory allocation |
| `SC_INGEST_ENABLED` | false | Enable scheduled ingest |
| `SC_INGEST_SCHEDULE_HOURS` | 24 | Ingest interval (hours) |
| `SC_CROSS_AGENT_TIMEOUT` | 30 | Cross-agent call timeout (s) |
| `SC_CORS_ORIGINS` | localhost origins | CORS allowed origins |
| `ANTHROPIC_API_KEY` | (none) | Anthropic API key |

### 5.2 Collection Weights

Collection search weights can be overridden via environment variables:

| Variable | Default |
|----------|---------|
| `SC_WEIGHT_CELL_TYPES` | 0.14 |
| `SC_WEIGHT_MARKERS` | 0.12 |
| `SC_WEIGHT_SPATIAL` | 0.10 |
| `SC_WEIGHT_TME` | 0.10 |
| `SC_WEIGHT_DRUG_RESPONSE` | 0.09 |
| `SC_WEIGHT_LITERATURE` | 0.08 |
| `SC_WEIGHT_METHODS` | 0.07 |
| `SC_WEIGHT_DATASETS` | 0.06 |
| `SC_WEIGHT_TRAJECTORIES` | 0.07 |
| `SC_WEIGHT_PATHWAYS` | 0.07 |
| `SC_WEIGHT_CLINICAL` | 0.07 |
| `SC_WEIGHT_GENOMIC` | 0.03 |

Weights must sum to approximately 1.0 (tolerance: 0.05).

---

## 6. Monitoring

### 6.1 Health Check

```bash
# Continuous health monitoring
watch -n 10 'curl -s http://localhost:8540/health | python -m json.tool'
```

### 6.2 Prometheus Metrics

```bash
# Scrape metrics endpoint
curl http://localhost:8540/metrics

# Example output:
# sc_agent_requests_total 42
# sc_agent_query_requests_total 15
# sc_agent_errors_total 0
```

### 6.3 Docker Logs

```bash
# API server logs
docker compose logs -f sc-api

# Streamlit UI logs
docker compose logs -f sc-streamlit

# Milvus logs
docker compose logs -f milvus-standalone
```

---

## 7. Scaling

### 7.1 Horizontal Scaling (API Workers)

Increase uvicorn workers for higher throughput:

```yaml
# docker-compose.yml override
sc-api:
  command:
    - uvicorn
    - api.main:app
    - --host=0.0.0.0
    - --port=8540
    - --workers=4  # Increase from 2 to 4
```

### 7.2 Milvus Scaling

For large-scale deployments, replace Milvus standalone with Milvus cluster:

```bash
# See Milvus documentation for cluster deployment
# Requires: etcd cluster, Kafka/Pulsar, multiple query/data/index nodes
```

---

## 8. Backup and Recovery

### 8.1 Milvus Data Backup

```bash
# Backup Milvus volumes
docker compose stop milvus-standalone
docker run --rm -v sc_milvus_data:/data -v $(pwd)/backup:/backup \
    alpine tar czf /backup/milvus-$(date +%Y%m%d).tar.gz /data
docker compose start milvus-standalone
```

### 8.2 Conversation Data Backup

```bash
# Backup conversation cache
tar czf conversations-$(date +%Y%m%d).tar.gz data/cache/conversations/
```

---

## 9. Troubleshooting

### 9.1 Common Issues

| Issue | Cause | Solution |
|-------|-------|---------|
| "degraded" health status | Milvus not connected | Check `docker compose ps milvus-standalone` |
| "LLM unavailable" | Missing API key | Set `ANTHROPIC_API_KEY` env var |
| Port conflict on 8540 | Another service using port | Change `SC_API_PORT` |
| "No collections" | Setup not run | Run `docker compose run sc-setup` |
| Slow queries | Cold Milvus cache | First query loads collections; subsequent queries are fast |
| Rate limit errors (429) | Exceeded 100 req/min | Wait 60 seconds or reduce query rate |

### 9.2 Reset and Rebuild

```bash
# Full reset (destroys all data)
docker compose down -v
docker compose up -d
docker compose logs -f sc-setup
```

---

## 10. Security Checklist

- [ ] Set `ANTHROPIC_API_KEY` via environment variable (never in code)
- [ ] Set `SC_API_KEY` to a strong random value for production
- [ ] Deploy behind TLS-terminating reverse proxy (nginx/traefik)
- [ ] Restrict `SC_CORS_ORIGINS` to known frontend domains
- [ ] Review Docker network isolation (`sc-network`)
- [ ] Ensure non-root container user (`scuser`) is active
- [ ] Configure firewall rules for ports 8540, 8130

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Deployment Guide*

---

!!! warning "Clinical Decision Support Disclaimer"
    The Single-Cell Intelligence Agent is a clinical decision support research tool for single-cell genomic analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
