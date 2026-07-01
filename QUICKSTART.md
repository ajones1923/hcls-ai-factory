# HCLS AI Factory -- Quick Start

Patient DNA to drug candidates in under 5 hours on a single NVIDIA DGX Spark.

## Prerequisites

- **Docker** (v24+) and Docker Compose v2
- **NVIDIA GPU** with CUDA 12.x drivers and NVIDIA Container Toolkit
- **Disk space**: ~1.1 TB for models, reference genomes, and vector data
- **API keys**: `ANTHROPIC_API_KEY` (required), optional `GRAFANA_PASSWORD`

## Quick Start

```bash
# 1. Clone the repo
git clone https://github.com/your-org/hcls-ai-factory.git
cd hcls-ai-factory

# 2. Set your API key
export ANTHROPIC_API_KEY="sk-ant-..."

# 3. Launch the full stack
docker compose -f docker-compose.dgx-spark.yml up -d

# 4. Watch startup progress
docker compose -f docker-compose.dgx-spark.yml logs -f --tail=50
```

All services will be healthy within 2-5 minutes depending on build cache.

## Where to Find Things

| Service                        | URL                          |
|--------------------------------|------------------------------|
| Landing Page / Portal          | http://localhost:8080         |
| Precision Biomarker Agent      | http://localhost:8502 (UI) / http://localhost:8102 (API) |
| Precision Oncology Agent       | http://localhost:8503 (UI) / http://localhost:8103 (API) |
| CAR-T Intelligence Agent       | http://localhost:8504 (UI) / http://localhost:8104 (API) |
| Imaging Intelligence Agent     | http://localhost:8505 (UI) / http://localhost:8105 (API) |
| Autoimmune Agent               | http://localhost:8506 (UI) / http://localhost:8106 (API) |
| Cardiology Agent               | http://localhost:8536 (UI) / http://localhost:8126 (API) |
| Pharmacogenomics Agent         | http://localhost:8507 (UI) / http://localhost:8107 (API) |
| Clinical Trial Agent           | http://localhost:8128 (UI) / http://localhost:8538 (API) |
| Rare Disease Diagnostic Agent  | http://localhost:8544 (UI) / http://localhost:8134 (API) |
| Neurology Agent                | http://localhost:8529 (UI) / http://localhost:8528 (API) |
| Single Cell Agent              | http://localhost:8130 (UI) / http://localhost:8540 (API) |
| Milvus (vector DB)             | localhost:19530               |
| Prometheus                     | http://localhost:9099         |
| Grafana                        | http://localhost:3000         |

## First Demo

Check that the cardiology risk calculator is running:

```bash
curl -s http://localhost:8126/health | python3 -m json.tool
```

Run a sample risk assessment:

```bash
curl -s -X POST http://localhost:8126/assess-risk \
  -H "Content-Type: application/json" \
  -d '{"patient_id": "demo-001", "variants": ["rs1801133", "rs4343"]}' \
  | python3 -m json.tool
```

## Stopping

```bash
docker compose -f docker-compose.dgx-spark.yml down
```

## Full Documentation

See [https://hcls-ai-factory.org](https://hcls-ai-factory.org) for architecture details, pipeline configuration, and deployment guides.
