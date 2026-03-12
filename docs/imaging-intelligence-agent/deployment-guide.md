# Imaging Intelligence Agent — Deployment Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**License:** Apache 2.0

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Quick Start — Docker Lite (No GPU)](#3-quick-start--docker-lite-no-gpu)
4. [Full Stack Deployment (Docker with GPU + NIM Services)](#4-full-stack-deployment-docker-with-gpu--nim-services)
5. [DGX Spark Production Deployment](#5-dgx-spark-production-deployment)
6. [Environment Variables Reference](#6-environment-variables-reference)
7. [Milvus Collection Setup](#7-milvus-collection-setup)
8. [NIM Service Configuration](#8-nim-service-configuration)
9. [PACS Integration (Orthanc + OHIF)](#9-pacs-integration-orthanc--ohif)
10. [Security Hardening](#10-security-hardening)
11. [Monitoring and Health Checks](#11-monitoring-and-health-checks)
12. [Troubleshooting](#12-troubleshooting)

---

## 1. Overview

The Imaging Intelligence Agent is a clinical decision support system for radiology,
providing multi-collection RAG-based question answering, NVIDIA NIM inference
microservices, reference clinical workflows, and PACS integration. It is one of
five intelligence agents in the HCLS AI Factory platform.

### What This Guide Covers

This guide walks through three deployment profiles:

| Profile | Services | GPU Required | Use Case |
|---------|----------|-------------|----------|
| **Docker Lite** | 6 services | No | Demo, testing, CI/CD, development |
| **Docker Full Stack** | 11 services | Yes (NVIDIA) | Full inference with NIM containers |
| **DGX Spark Production** | 11 services | Yes (DGX Spark) | Production deployment on NVIDIA DGX Spark |

### Service Inventory

**Docker Full Stack — 11 services:**

| # | Service | Image | Port(s) |
|---|---------|-------|---------|
| 1 | Orthanc DICOM Server | `orthancteam/orthanc:24.12.1` | 8042 (REST), 4242 (DICOM) |
| 2 | OHIF Viewer | `ohif/app:v3.9.2` | 8526 |
| 3 | etcd | `quay.io/coreos/etcd:v3.5.5` | 2379 (internal) |
| 4 | MinIO | `minio/minio:RELEASE.2023-03-20T20-16-18Z` | 9000/9001 (internal) |
| 5 | Milvus | `milvusdb/milvus:v2.4.17` | 19530, 9091 |
| 6 | Imaging Streamlit UI | Built from Dockerfile | 8525 |
| 7 | Imaging FastAPI Server | Built from Dockerfile | 8524 |
| 8 | Imaging Setup (one-shot) | Built from Dockerfile | — |
| 9 | NIM LLM (Llama-3 8B) | `nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest` | 8520 |
| 10 | NIM VISTA-3D | `nvcr.io/nvidia/nim/vista3d:latest` | 8530 |
| 11 | NIM MAISI | `nvcr.io/nvidia/nim/maisi:latest` | 8531 |
| 12 | NIM VILA-M3 | `nvcr.io/nvidia/nim/vilam3:latest` | 8532 |

**Docker Lite Stack — 6 services:**

| # | Service | Notes |
|---|---------|-------|
| 1 | etcd | Milvus metadata store |
| 2 | MinIO | Milvus object storage |
| 3 | Milvus | Vector database |
| 4 | Imaging Streamlit UI | `NIM_MODE=mock` |
| 5 | Imaging FastAPI Server | `NIM_MODE=mock` |
| 6 | Imaging Setup (one-shot) | Collection creation + data seeding |

---

## 2. Prerequisites

### 2.1 Hardware Requirements

#### Lite Mode (No GPU)

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 4 cores | 8 cores |
| RAM | 8 GB | 16 GB |
| Disk | 10 GB | 20 GB |
| GPU | Not required | Not required |

#### Full Mode (With NIMs)

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 8 cores | 16 cores |
| RAM | 32 GB | 64 GB |
| Disk | 100 GB | 200 GB |
| GPU | 1x NVIDIA GPU (24 GB VRAM) | 1x NVIDIA A100/H100 (80 GB VRAM) |

#### DGX Spark Production

| Resource | Specification |
|----------|--------------|
| System | NVIDIA DGX Spark |
| GPU | NVIDIA Grace Blackwell |
| RAM | 128 GB unified memory |
| Storage | NVMe SSD |
| Network | 10 GbE minimum |

### 2.2 Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| Docker Engine | 24.0+ | Container runtime |
| Docker Compose | v2.20+ | Multi-container orchestration |
| NVIDIA Container Toolkit | 1.14+ | GPU passthrough (Full Mode only) |
| NVIDIA Driver | 535+ | GPU support (Full Mode only) |
| Python | 3.10+ | Local development only |
| Git | 2.30+ | Source code management |
| curl | Any | Health checks and testing |

#### Install Docker Engine

```bash
# Ubuntu / Debian
sudo apt-get update
sudo apt-get install -y docker.io docker-compose-plugin
sudo usermod -aG docker $USER
newgrp docker
```

#### Install NVIDIA Container Toolkit (Full Mode only)

```bash
# Add NVIDIA package repository
distribution=$(. /etc/os-release; echo $ID$VERSION_ID) \
  && curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey \
    | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg \
  && curl -s -L "https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list" \
    | sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' \
    | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt-get update
sudo apt-get install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

Verify GPU access from Docker:

```bash
docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
```

### 2.3 API Keys

| Key | Required For | How to Obtain |
|-----|-------------|---------------|
| `ANTHROPIC_API_KEY` | Claude LLM synthesis (RAG answers) | [console.anthropic.com](https://console.anthropic.com/) |
| `NGC_API_KEY` | Pulling NIM containers from NGC | [ngc.nvidia.com](https://ngc.nvidia.com/) > Setup > API Key |
| `NVIDIA_API_KEY` | Cloud NIM endpoints (optional) | [build.nvidia.com](https://build.nvidia.com/) |
| `NCBI_API_KEY` | PubMed ingest rate limit increase (optional) | [ncbi.nlm.nih.gov/account](https://www.ncbi.nlm.nih.gov/account/) |

#### NGC Authentication (Full Mode only)

```bash
# Log in to NGC container registry
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
```

---

## 3. Quick Start — Docker Lite (No GPU)

The Lite stack runs 6 services with all NIM-dependent features in mock mode, returning
synthetic responses. This is the fastest way to explore the RAG knowledge system,
Streamlit UI, and FastAPI endpoints without any GPU hardware.

### 3.1 Clone and Configure

```bash
# Navigate to the agent directory
cd /path/to/hcls-ai-factory/ai_agent_adds/imaging_intelligence_agent/agent

# Copy the environment template
cp .env.example .env
```

Edit `.env` and set your Anthropic API key:

```bash
# .env
ANTHROPIC_API_KEY=sk-ant-api03-your-key-here
IMAGING_NIM_MODE=mock
```

### 3.2 Launch the Lite Stack

```bash
# Start all 6 services
docker compose -f docker-compose.lite.yml up -d

# Watch the setup container seed data into Milvus
docker compose -f docker-compose.lite.yml logs -f imaging-setup
```

The setup container will:
1. Create all 10 Milvus collections with IVF_FLAT indexes
2. Seed findings, protocols, devices, anatomy, benchmarks, guidelines, report templates, and datasets
3. Exit with code 0 on success

Expected output (final lines):

```
imaging-setup  | ===== Seeding datasets =====
imaging-setup  | ===== Imaging Agent Setup complete! =====
imaging-setup exited with code 0
```

### 3.3 Verify Services

```bash
# Check all containers are running (setup will be Exited/0)
docker compose -f docker-compose.lite.yml ps

# Test FastAPI health endpoint
curl -s http://localhost:8524/health | python3 -m json.tool

# Test Milvus health
curl -s http://localhost:9091/healthz
```

### 3.4 Access the UI

| Service | URL |
|---------|-----|
| Streamlit Chat UI | [http://localhost:8525](http://localhost:8525) |
| FastAPI Docs (Swagger) | [http://localhost:8524/docs](http://localhost:8524/docs) |
| FastAPI Health | [http://localhost:8524/health](http://localhost:8524/health) |

### 3.5 Test a RAG Query

```bash
curl -s -X POST http://localhost:8524/query \
  -H "Content-Type: application/json" \
  -d '{
    "query": "What is the sensitivity of CT for detecting pulmonary embolism?",
    "top_k": 5
  }' | python3 -m json.tool
```

### 3.6 Stop the Lite Stack

```bash
docker compose -f docker-compose.lite.yml down

# To also remove volumes (deletes all Milvus data):
docker compose -f docker-compose.lite.yml down -v
```

---

## 4. Full Stack Deployment (Docker with GPU + NIM Services)

The Full Stack adds Orthanc DICOM server, OHIF viewer, and 4 NVIDIA NIM containers
for on-device medical imaging inference. This requires an NVIDIA GPU with the
Container Toolkit installed.

### 4.1 Configure Environment

```bash
cd /path/to/hcls-ai-factory/ai_agent_adds/imaging_intelligence_agent/agent

cp .env.example .env
```

Edit `.env`:

```bash
# .env — Full Stack Configuration
ANTHROPIC_API_KEY=sk-ant-api03-your-key-here
NGC_API_KEY=your-ngc-api-key-here

# NIM mode: "local" uses the Docker NIM containers
IMAGING_NIM_MODE=local

# NIM Service URLs (Docker service names resolve inside the network)
IMAGING_NIM_LLM_URL=http://nim-llm:8000
IMAGING_NIM_VISTA3D_URL=http://nim-vista3d:8000
IMAGING_NIM_MAISI_URL=http://nim-maisi:8000
IMAGING_NIM_VILAM3_URL=http://nim-vilam3:8000

# Milvus
IMAGING_MILVUS_HOST=milvus-standalone
IMAGING_MILVUS_PORT=19530

# Orthanc
IMAGING_ORTHANC_URL=http://orthanc:8042
IMAGING_ORTHANC_USERNAME=admin
IMAGING_ORTHANC_PASSWORD=orthanc

# Ports
IMAGING_API_PORT=8524
IMAGING_STREAMLIT_PORT=8525
```

### 4.2 Pull NIM Container Images

NIM images are large (10-30 GB each). Pull them before launching to avoid timeouts:

```bash
# Authenticate with NGC
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin

# Pull all NIM images (may take 30-60 minutes on first run)
docker pull nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest
docker pull nvcr.io/nvidia/nim/vista3d:latest
docker pull nvcr.io/nvidia/nim/maisi:latest
docker pull nvcr.io/nvidia/nim/vilam3:latest
```

### 4.3 Launch the Full Stack

```bash
# Start all 11 services
docker compose up -d

# Monitor startup progress
docker compose logs -f
```

NIM containers take 2-5 minutes to initialize and load models onto the GPU. Watch
for health check passes:

```bash
# Check NIM LLM readiness
docker compose logs nim-llm 2>&1 | tail -20

# Check all service status
docker compose ps
```

### 4.4 Verify Full Stack

```bash
# FastAPI health (includes NIM status)
curl -s http://localhost:8524/health | python3 -m json.tool

# NIM LLM health
curl -s http://localhost:8520/v1/health/ready

# NIM VISTA-3D health
curl -s http://localhost:8530/v1/health/ready

# NIM MAISI health
curl -s http://localhost:8531/v1/health/ready

# NIM VILA-M3 health
curl -s http://localhost:8532/v1/health/ready

# Orthanc health
curl -s http://localhost:8042/system | python3 -m json.tool

# Milvus health
curl -s http://localhost:9091/healthz
```

### 4.5 Access Services

| Service | URL |
|---------|-----|
| Streamlit Chat UI | [http://localhost:8525](http://localhost:8525) |
| FastAPI Docs (Swagger) | [http://localhost:8524/docs](http://localhost:8524/docs) |
| OHIF DICOM Viewer | [http://localhost:8526](http://localhost:8526) |
| Orthanc Explorer | [http://localhost:8042](http://localhost:8042) |
| NIM LLM | [http://localhost:8520](http://localhost:8520) |
| Milvus Metrics | [http://localhost:9091/metrics](http://localhost:9091/metrics) |

### 4.6 Test NIM Inference

```bash
# Test LLM generation via NIM
curl -s -X POST http://localhost:8520/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "meta-llama3-8b-instruct",
    "messages": [
      {"role": "user", "content": "What are the key findings on a chest CT for pulmonary embolism?"}
    ],
    "max_tokens": 256
  }' | python3 -m json.tool

# Test NIM status via FastAPI proxy
curl -s http://localhost:8524/nim/status | python3 -m json.tool
```

### 4.7 Upload a DICOM Study to Orthanc

```bash
# Upload a DICOM file via Orthanc REST API
curl -X POST http://localhost:8042/instances \
  -u admin:orthanc \
  -H "Content-Type: application/dicom" \
  --data-binary @/path/to/study.dcm

# List studies
curl -s http://localhost:8042/studies -u admin:orthanc | python3 -m json.tool

# Upload via DICOM C-STORE (from any DICOM toolkit)
# storescu localhost 4242 /path/to/study.dcm
```

### 4.8 Stop the Full Stack

```bash
# Graceful shutdown
docker compose down

# Shutdown and remove volumes
docker compose down -v

# Shutdown, remove volumes, and remove images
docker compose down -v --rmi all
```

---

## 5. DGX Spark Production Deployment

The DGX Spark deployment remaps external ports to avoid conflicts with other
HCLS AI Factory agents running on the same host.

### 5.1 Port Mapping

| Service | Internal Port | External Port (DGX Spark) |
|---------|--------------|---------------------------|
| Streamlit UI | 8525 | **8505** |
| FastAPI Server | 8524 | **8105** |
| Orthanc REST | 8042 | 8042 |
| Orthanc DICOM | 4242 | 4242 |
| OHIF Viewer | 80 (container) | 8526 |
| NIM LLM | 8000 (container) | 8520 |
| NIM VISTA-3D | 8000 (container) | 8530 |
| NIM MAISI | 8000 (container) | 8531 |
| NIM VILA-M3 | 8000 (container) | 8532 |
| Milvus gRPC | 19530 | 19530 |
| Milvus Metrics | 9091 | 9091 |

### 5.2 Create Production Override

Create a `docker-compose.dgx-spark.yml` override file:

```yaml
# docker-compose.dgx-spark.yml
# DGX Spark production overrides — remaps Streamlit and FastAPI ports

services:
  imaging-streamlit:
    ports:
      - "8505:8525"
    environment:
      IMAGING_NIM_MODE: local
      IMAGING_NIM_ALLOW_MOCK_FALLBACK: "false"
    restart: always
    deploy:
      resources:
        limits:
          memory: 4G

  imaging-api:
    ports:
      - "8105:8524"
    environment:
      IMAGING_NIM_MODE: local
      IMAGING_NIM_ALLOW_MOCK_FALLBACK: "false"
      IMAGING_CORS_ORIGINS: "http://localhost:8080,http://localhost:8505,http://localhost:8105"
    restart: always
    deploy:
      resources:
        limits:
          memory: 4G

  orthanc:
    restart: always

  milvus-standalone:
    restart: always

  nim-llm:
    restart: always
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
        limits:
          memory: 32G

  nim-vista3d:
    restart: always
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
        limits:
          memory: 16G

  nim-maisi:
    restart: always
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
        limits:
          memory: 16G

  nim-vilam3:
    restart: always
    deploy:
      resources:
        reservations:
          devices:
            - driver: nvidia
              count: 1
              capabilities: [gpu]
        limits:
          memory: 16G
```

### 5.3 Launch on DGX Spark

```bash
cd /path/to/hcls-ai-factory/ai_agent_adds/imaging_intelligence_agent/agent

# Configure environment
cp .env.example .env
# Edit .env with production API keys

# Launch with DGX Spark overrides
docker compose -f docker-compose.yml -f docker-compose.dgx-spark.yml up -d

# Monitor startup
docker compose -f docker-compose.yml -f docker-compose.dgx-spark.yml logs -f
```

### 5.4 Verify DGX Spark Deployment

```bash
# Verify remapped ports
curl -s http://localhost:8105/health | python3 -m json.tool
curl -s http://localhost:8505/_stcore/health

# Verify NIM services
curl -s http://localhost:8520/v1/health/ready
curl -s http://localhost:8530/v1/health/ready
curl -s http://localhost:8531/v1/health/ready
curl -s http://localhost:8532/v1/health/ready

# Verify GPU allocation
nvidia-smi
docker inspect imaging-nim-llm | grep -A 5 DeviceRequests
```

### 5.5 Systemd Service (Auto-Start on Boot)

Create `/etc/systemd/system/imaging-agent.service`:

```ini
[Unit]
Description=Imaging Intelligence Agent — HCLS AI Factory
After=docker.service nvidia-persistenced.service
Requires=docker.service

[Service]
Type=oneshot
RemainAfterExit=yes
WorkingDirectory=/path/to/hcls-ai-factory/ai_agent_adds/imaging_intelligence_agent/agent
ExecStart=/usr/bin/docker compose -f docker-compose.yml -f docker-compose.dgx-spark.yml up -d
ExecStop=/usr/bin/docker compose -f docker-compose.yml -f docker-compose.dgx-spark.yml down
TimeoutStartSec=300
TimeoutStopSec=120

[Install]
WantedBy=multi-user.target
```

Enable and start:

```bash
sudo systemctl daemon-reload
sudo systemctl enable imaging-agent.service
sudo systemctl start imaging-agent.service
sudo systemctl status imaging-agent.service
```

### 5.6 Log Rotation

Create `/etc/logrotate.d/imaging-agent`:

```
/var/lib/docker/containers/*/*.log {
    daily
    rotate 7
    compress
    delaycompress
    missingok
    notifempty
    maxsize 100M
    copytruncate
}
```

---

## 6. Environment Variables Reference

All settings use the `IMAGING_` prefix and are managed by Pydantic BaseSettings
in `config/settings.py`. Variables are loaded from environment variables first,
then from `.env` file.

### 6.1 Milvus / Vector Database

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_MILVUS_HOST` | `localhost` | Milvus server hostname |
| `IMAGING_MILVUS_PORT` | `19530` | Milvus gRPC port |
| `IMAGING_EMBEDDING_MODEL` | `BAAI/bge-small-en-v1.5` | HuggingFace embedding model name |
| `IMAGING_EMBEDDING_DIMENSION` | `384` | Embedding vector dimension |
| `IMAGING_EMBEDDING_BATCH_SIZE` | `32` | Batch size for embedding generation |
| `IMAGING_TOP_K_PER_COLLECTION` | `5` | Number of results per collection in RAG search |
| `IMAGING_SCORE_THRESHOLD` | `0.4` | Minimum cosine similarity score to include |

### 6.2 LLM Configuration

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_LLM_PROVIDER` | `anthropic` | LLM provider (`anthropic`) |
| `IMAGING_LLM_MODEL` | `claude-sonnet-4-6` | Claude model for RAG synthesis |
| `IMAGING_ANTHROPIC_API_KEY` | — | Anthropic API key (required) |

### 6.3 NIM Configuration

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_NIM_MODE` | `local` | NIM mode: `local`, `cloud`, or `mock` |
| `IMAGING_NIM_ALLOW_MOCK_FALLBACK` | `True` | Fall back to mock if NIM unavailable |
| `IMAGING_NIM_LLM_URL` | `http://localhost:8520/v1` | NIM LLM endpoint URL |
| `IMAGING_NIM_VISTA3D_URL` | `http://localhost:8530` | NIM VISTA-3D endpoint URL |
| `IMAGING_NIM_MAISI_URL` | `http://localhost:8531` | NIM MAISI endpoint URL |
| `IMAGING_NIM_VILAM3_URL` | `http://localhost:8532` | NIM VILA-M3 endpoint URL |
| `IMAGING_NGC_API_KEY` | — | NGC API key for NIM containers |
| `IMAGING_NIM_LOCAL_LLM_MODEL` | `meta/llama3-70b-instruct` | Local NIM LLM model name |

### 6.4 NVIDIA Cloud NIM

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_NVIDIA_API_KEY` | — | NVIDIA API key for cloud NIM endpoints |
| `IMAGING_NIM_CLOUD_URL` | `https://integrate.api.nvidia.com/v1` | Cloud NIM base URL |
| `IMAGING_NIM_CLOUD_LLM_MODEL` | `meta/llama-3.1-8b-instruct` | Cloud NIM LLM model |
| `IMAGING_NIM_CLOUD_VLM_MODEL` | `meta/llama-3.2-11b-vision-instruct` | Cloud NIM VLM model |

### 6.5 API Server

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_API_HOST` | `0.0.0.0` | FastAPI bind address |
| `IMAGING_API_PORT` | `8524` | FastAPI listen port |
| `IMAGING_STREAMLIT_PORT` | `8525` | Streamlit listen port |
| `IMAGING_API_BASE_URL` | `http://localhost:8524` | Internal URL Streamlit uses to reach FastAPI |
| `IMAGING_CORS_ORIGINS` | `http://localhost:8080,http://localhost:8524,http://localhost:8525` | Allowed CORS origins (comma-separated) |
| `IMAGING_MAX_REQUEST_SIZE_MB` | `10` | Maximum request body size in MB |

### 6.6 PACS / DICOM

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_ORTHANC_URL` | `http://localhost:8042` | Orthanc REST API URL |
| `IMAGING_ORTHANC_USERNAME` | `admin` | Orthanc HTTP basic auth username |
| `IMAGING_ORTHANC_PASSWORD` | *(empty)* | Orthanc HTTP basic auth password |
| `IMAGING_OHIF_URL` | `http://localhost:8526` | OHIF viewer URL |
| `IMAGING_DICOM_AUTO_INGEST` | `False` | Enable automatic DICOM study ingestion |
| `IMAGING_DICOM_WATCH_INTERVAL` | `5` | Seconds between Orthanc /changes polls |
| `IMAGING_DICOM_SERVER_URL` | `http://localhost:8042` | Legacy alias for Orthanc URL |

### 6.7 Preview Generation

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_PREVIEW_CACHE_DIR` | `data/cache/previews` | Directory for cached preview files |
| `IMAGING_PREVIEW_DEFAULT_FPS` | `8` | Default frames per second for video previews |
| `IMAGING_PREVIEW_DEFAULT_FORMAT` | `mp4` | Default preview format (`mp4`, `gif`) |
| `IMAGING_PREVIEW_MAX_FRAMES` | `200` | Maximum number of frames per preview |

### 6.8 Cross-Modal and Scheduling

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_CROSS_MODAL_ENABLED` | `False` | Enable cross-modal genomics enrichment |
| `IMAGING_INGEST_ENABLED` | `False` | Enable scheduled data ingestion |
| `IMAGING_INGEST_SCHEDULE_HOURS` | `168` | Ingest interval in hours (168 = weekly) |

### 6.9 RAG Search Weights

Collection weights control how results from each collection are scored. They should
sum to approximately 1.0.

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_WEIGHT_LITERATURE` | `0.18` | Weight for imaging_literature |
| `IMAGING_WEIGHT_TRIALS` | `0.12` | Weight for imaging_trials |
| `IMAGING_WEIGHT_FINDINGS` | `0.15` | Weight for imaging_findings |
| `IMAGING_WEIGHT_PROTOCOLS` | `0.08` | Weight for imaging_protocols |
| `IMAGING_WEIGHT_DEVICES` | `0.08` | Weight for imaging_devices |
| `IMAGING_WEIGHT_ANATOMY` | `0.06` | Weight for imaging_anatomy |
| `IMAGING_WEIGHT_BENCHMARKS` | `0.08` | Weight for imaging_benchmarks |
| `IMAGING_WEIGHT_GUIDELINES` | `0.10` | Weight for imaging_guidelines |
| `IMAGING_WEIGHT_REPORT_TEMPLATES` | `0.05` | Weight for imaging_report_templates |
| `IMAGING_WEIGHT_DATASETS` | `0.06` | Weight for imaging_datasets |
| `IMAGING_WEIGHT_GENOMIC` | `0.04` | Weight for genomic_evidence |

### 6.10 Monitoring and Observability

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_METRICS_ENABLED` | `True` | Enable Prometheus metrics endpoint |

### 6.11 Conversation and Citation

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_MAX_CONVERSATION_CONTEXT` | `3` | Number of prior exchanges to include |
| `IMAGING_CITATION_HIGH_THRESHOLD` | `0.75` | Cosine similarity for "high confidence" |
| `IMAGING_CITATION_MEDIUM_THRESHOLD` | `0.60` | Cosine similarity for "medium confidence" |

### 6.12 PubMed

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_NCBI_API_KEY` | — | NCBI API key (optional, increases rate limits) |
| `IMAGING_PUBMED_MAX_RESULTS` | `5000` | Max results for PubMed ingest |

### 6.13 ClinicalTrials.gov

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_CT_GOV_BASE_URL` | `https://clinicaltrials.gov/api/v2` | ClinicalTrials.gov API base URL |

---

## 7. Milvus Collection Setup

The Imaging Intelligence Agent uses 10 dedicated Milvus collections plus 1 shared
read-only collection from the Stage 2 RAG pipeline. All collections use
BGE-small-en-v1.5 embeddings (384 dimensions), IVF_FLAT index type, and COSINE
distance metric.

### 7.1 Collection Inventory

| # | Collection Name | Content | Source |
|---|----------------|---------|--------|
| 1 | `imaging_literature` | PubMed research papers and reviews | PubMed NCBI E-utilities |
| 2 | `imaging_trials` | ClinicalTrials.gov AI-in-imaging records | ClinicalTrials.gov API v2 |
| 3 | `imaging_findings` | Imaging finding templates and patterns | Seed data |
| 4 | `imaging_protocols` | Acquisition protocols and parameters | Seed data |
| 5 | `imaging_devices` | FDA-cleared AI/ML medical devices | Seed data |
| 6 | `imaging_anatomy` | Anatomical structure references | Seed data |
| 7 | `imaging_benchmarks` | Model performance benchmarks | Seed data |
| 8 | `imaging_guidelines` | Clinical practice guidelines (ACR, RSNA, NCCN) | Seed data |
| 9 | `imaging_report_templates` | Structured radiology report templates | Seed data |
| 10 | `imaging_datasets` | Public imaging datasets (TCIA, PhysioNet) | Seed data |
| 11 | `genomic_evidence` | Shared from Stage 2 RAG pipeline (read-only) | Pre-existing |

### 7.2 Automated Setup

The `imaging-setup` container runs automatically on `docker compose up` and handles
collection creation and data seeding:

```bash
# Watch setup progress
docker compose logs -f imaging-setup
```

### 7.3 Manual Setup (Local Development)

If running Milvus standalone (not via Docker Compose), create collections and seed
data manually:

```bash
# Ensure Milvus is running on localhost:19530

# Create all 10 collections (drops existing if --drop-existing flag used)
python scripts/setup_collections.py --drop-existing

# Seed each collection
python scripts/seed_findings.py
python scripts/seed_protocols.py
python scripts/seed_devices.py
python scripts/seed_anatomy.py
python scripts/seed_benchmarks.py
python scripts/seed_guidelines.py
python scripts/seed_report_templates.py
python scripts/seed_datasets.py
```

### 7.4 Ingest Live Data

After initial seeding, ingest real-world data from PubMed and ClinicalTrials.gov:

```bash
# Ingest PubMed literature (fetches up to 5000 papers)
python scripts/ingest_pubmed.py

# Ingest clinical trials
python scripts/ingest_clinical_trials.py
```

### 7.5 Verify Collections

```bash
# List all collections with record counts
curl -s http://localhost:8524/collections | python3 -m json.tool
```

Or via Python:

```python
from pymilvus import connections, utility

connections.connect(host="localhost", port=19530)

for name in utility.list_collections():
    if name.startswith("imaging_") or name == "genomic_evidence":
        from pymilvus import Collection
        c = Collection(name)
        c.load()
        print(f"{name}: {c.num_entities} records")
```

### 7.6 Backup and Restore

```bash
# Backup: flush all data to disk
curl -X POST http://localhost:9091/api/v1/persist

# The Milvus data is stored in the milvus_data Docker volume
docker volume inspect imaging_intelligence_agent_milvus_data

# Create a tarball backup
docker run --rm \
  -v imaging_intelligence_agent_milvus_data:/data \
  -v $(pwd)/backups:/backup \
  alpine tar czf /backup/milvus-backup-$(date +%Y%m%d).tar.gz -C /data .

# Restore from backup
docker compose down
docker volume rm imaging_intelligence_agent_milvus_data
docker volume create imaging_intelligence_agent_milvus_data
docker run --rm \
  -v imaging_intelligence_agent_milvus_data:/data \
  -v $(pwd)/backups:/backup \
  alpine tar xzf /backup/milvus-backup-20260311.tar.gz -C /data
docker compose up -d
```

---

## 8. NIM Service Configuration

The agent integrates four NVIDIA NIM microservices for on-device medical imaging
inference. Each NIM container exposes an OpenAI-compatible REST API on port 8000
internally, mapped to unique external ports.

### 8.1 NIM Service Overview

| NIM Service | Image | Internal Port | External Port | GPU Memory |
|-------------|-------|---------------|---------------|------------|
| **LLM** (Llama-3 8B Instruct) | `nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest` | 8000 | 8520 | ~16 GB |
| **VISTA-3D** (3D Segmentation) | `nvcr.io/nvidia/nim/vista3d:latest` | 8000 | 8530 | ~8 GB |
| **MAISI** (Synthetic CT) | `nvcr.io/nvidia/nim/maisi:latest` | 8000 | 8531 | ~8 GB |
| **VILA-M3** (Vision-Language) | `nvcr.io/nvidia/nim/vilam3:latest` | 8000 | 8532 | ~8 GB |

### 8.2 NIM Modes

The `IMAGING_NIM_MODE` variable controls how NIM services are accessed:

| Mode | Description | When to Use |
|------|-------------|-------------|
| `local` | Connect to local Docker NIM containers | Full Stack deployment with GPU |
| `cloud` | Connect to NVIDIA cloud NIM endpoints | No local GPU; pay-per-call |
| `mock` | Return synthetic responses | Testing, CI/CD, Lite mode |

```bash
# Local mode (default for Full Stack)
IMAGING_NIM_MODE=local

# Cloud mode (requires IMAGING_NVIDIA_API_KEY)
IMAGING_NIM_MODE=cloud
IMAGING_NVIDIA_API_KEY=nvapi-your-key-here
IMAGING_NIM_CLOUD_URL=https://integrate.api.nvidia.com/v1

# Mock mode (no GPU, no API key needed)
IMAGING_NIM_MODE=mock
```

### 8.3 Mock Fallback

When `IMAGING_NIM_ALLOW_MOCK_FALLBACK=True` (default), the agent gracefully degrades
to mock responses if a NIM service is unreachable. This prevents hard failures during
development. Disable in production:

```bash
IMAGING_NIM_ALLOW_MOCK_FALLBACK=false
```

### 8.4 NIM Health Checks

All NIM containers expose `/v1/health/ready`:

```bash
# Check each service
curl -s http://localhost:8520/v1/health/ready  # LLM
curl -s http://localhost:8530/v1/health/ready  # VISTA-3D
curl -s http://localhost:8531/v1/health/ready  # MAISI
curl -s http://localhost:8532/v1/health/ready  # VILA-M3

# Check all via FastAPI proxy
curl -s http://localhost:8524/nim/status | python3 -m json.tool
```

### 8.5 NIM Client Architecture

All NIM clients extend `BaseNIMClient` (defined in `src/nim/base.py`), which provides:

- **Cached health checks:** Results cached for 30 seconds to reduce polling overhead
- **Exponential backoff retry:** Via tenacity, retries on transient failures
- **Automatic mock fallback:** Returns synthetic results when service is unavailable
- **Request/response logging:** Via loguru for debugging

```
BaseNIMClient (ABC)
  +-- LLMClient         (src/nim/llm_client.py)
  +-- VISTA3DClient     (src/nim/vista3d_client.py)
  +-- MAISIClient       (src/nim/maisi_client.py)
  +-- VILAM3Client      (src/nim/vilam3_client.py)
```

### 8.6 Multi-GPU Allocation

On systems with multiple GPUs, you can pin each NIM to a specific GPU using
`NVIDIA_VISIBLE_DEVICES` or the `count`/`device_ids` fields in the compose file:

```yaml
# Example: assign specific GPUs to each NIM
nim-llm:
  environment:
    NVIDIA_VISIBLE_DEVICES: "0"
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            device_ids: ["0"]
            capabilities: [gpu]

nim-vista3d:
  environment:
    NVIDIA_VISIBLE_DEVICES: "1"
  deploy:
    resources:
      reservations:
        devices:
          - driver: nvidia
            device_ids: ["1"]
            capabilities: [gpu]
```

### 8.7 Testing NIM Services

```bash
# Run the NIM test script
python scripts/test_nim_services.py
```

---

## 9. PACS Integration (Orthanc + OHIF)

The Full Stack deployment includes Orthanc as a DICOM server and OHIF as a
zero-footprint web viewer for medical images.

### 9.1 Orthanc Configuration

Orthanc runs as the `orthanc` service with DICOMweb enabled:

| Setting | Value |
|---------|-------|
| Image | `orthancteam/orthanc:24.12.1` |
| REST API Port | 8042 |
| DICOM Port | 4242 |
| AE Title | `IMAGING_AI` |
| DICOMweb | Enabled at `/dicom-web/` |
| Auth | HTTP basic: `admin` / `orthanc` |
| Storage | Docker volume `orthanc_data` |

### 9.2 OHIF Viewer Configuration

OHIF connects to Orthanc via DICOMweb. The configuration file is mounted at
`config/ohif-config.js`:

| Setting | Value |
|---------|-------|
| Image | `ohif/app:v3.9.2` |
| Port | 8526 (maps to container port 80) |
| Data Source | Orthanc DICOMweb |
| WADO-RS | `http://<hostname>:8042/dicom-web` |
| QIDO-RS | `http://<hostname>:8042/dicom-web` |

The OHIF config uses `window.location.hostname` to dynamically resolve Orthanc URLs,
so it works on localhost, LAN IPs, and remote hosts.

### 9.3 Sending Studies to Orthanc

#### Via REST API (HTTP)

```bash
# Upload a single DICOM file
curl -X POST http://localhost:8042/instances \
  -u admin:orthanc \
  -H "Content-Type: application/dicom" \
  --data-binary @/path/to/image.dcm

# Upload a directory of DICOM files
for f in /path/to/dicom/directory/*.dcm; do
  curl -X POST http://localhost:8042/instances \
    -u admin:orthanc \
    -H "Content-Type: application/dicom" \
    --data-binary @"$f"
done

# Upload a ZIP archive
curl -X POST http://localhost:8042/instances \
  -u admin:orthanc \
  -H "Content-Type: application/zip" \
  --data-binary @/path/to/study.zip
```

#### Via DICOM C-STORE (Network)

```bash
# Using dcmtk's storescu
storescu localhost 4242 /path/to/image.dcm

# Using pynetdicom
python -c "
from pynetdicom import AE
ae = AE(ae_title='SENDING_SCU')
ae.add_requested_context('1.2.840.10008.5.1.4.1.1.2')  # CT
assoc = ae.associate('localhost', 4242, ae_title='IMAGING_AI')
if assoc.is_established:
    print('Connected to Orthanc')
    assoc.release()
"
```

### 9.4 Viewing Studies

1. Open OHIF at [http://localhost:8526](http://localhost:8526)
2. The study list shows all studies stored in Orthanc
3. Click a study to open the viewer with MPR, measurements, and annotations

### 9.5 DICOM Auto-Ingestion

When enabled, the agent polls Orthanc for new studies and automatically triggers
clinical workflows:

```bash
# Enable in .env
IMAGING_DICOM_AUTO_INGEST=True
IMAGING_DICOM_WATCH_INTERVAL=5
```

The DICOM watcher (`src/ingest/dicom_watcher.py`) polls Orthanc's `/changes` endpoint
every 5 seconds and dispatches new `StableStudy` events to the workflow engine via
the `/events/dicom-webhook` endpoint.

### 9.6 DICOM Event Webhook

External PACS systems can POST DICOM events directly:

```bash
curl -X POST http://localhost:8524/events/dicom-webhook \
  -H "Content-Type: application/json" \
  -d '{
    "event_type": "StableStudy",
    "study_id": "1.2.3.4.5.6.7.8.9",
    "modality": "CT",
    "body_part": "HEAD",
    "description": "CT Head without contrast"
  }'
```

---

## 10. Security Hardening

### 10.1 API Key Management

**Never commit API keys to version control.**

```bash
# Use .env file (already in .gitignore)
cp .env.example .env
chmod 600 .env

# Verify .env is gitignored
grep '.env' .gitignore
```

For production, use Docker secrets or environment variable injection from a vault:

```bash
# Docker secrets approach
echo "sk-ant-api03-..." | docker secret create anthropic_api_key -
echo "nvapi-..." | docker secret create ngc_api_key -
```

Or export from a secrets manager before launching:

```bash
# Example with environment variable export
export ANTHROPIC_API_KEY=$(vault kv get -field=api_key secret/imaging/anthropic)
export NGC_API_KEY=$(vault kv get -field=api_key secret/imaging/ngc)
docker compose up -d
```

### 10.2 CORS Configuration

By default, CORS allows requests from the landing page, FastAPI, and Streamlit:

```bash
IMAGING_CORS_ORIGINS=http://localhost:8080,http://localhost:8524,http://localhost:8525
```

For production, restrict to your actual domain:

```bash
IMAGING_CORS_ORIGINS=https://imaging.yourdomain.com,https://api.yourdomain.com
```

### 10.3 Non-Root Docker Containers

The Dockerfile creates and runs as a non-root user (`imaginguser`):

```dockerfile
RUN useradd -r -m -s /bin/false imaginguser \
    && mkdir -p /app/data/cache /app/data/reference /app/data/dicom \
    && mkdir -p /home/imaginguser/.cache/huggingface \
    && chown -R imaginguser:imaginguser /app /home/imaginguser
USER imaginguser
```

Verify the container is not running as root:

```bash
docker exec imaging-api whoami
# Expected: imaginguser

docker exec imaging-streamlit id
# Expected: uid=999(imaginguser) gid=999(imaginguser) groups=999(imaginguser)
```

### 10.4 Orthanc Authentication

Default credentials (`admin` / `orthanc`) must be changed for production:

```yaml
# In docker-compose.yml, update the orthanc service environment:
orthanc:
  environment:
    ORTHANC__REGISTERED_USERS: |
      {"your_admin_user": "your_strong_password_here"}
    ORTHANC__HTTP_AUTHENTICATION_ENABLED: "true"
```

Update the FastAPI environment to match:

```bash
IMAGING_ORTHANC_USERNAME=your_admin_user
IMAGING_ORTHANC_PASSWORD=your_strong_password_here
```

### 10.5 Network Isolation

All services communicate over the `imaging-network` bridge network. Only essential
ports are exposed to the host:

```bash
# Verify network isolation
docker network inspect imaging-network

# Only these ports should be bound to host interfaces:
# 8524 (FastAPI), 8525 (Streamlit), 8042 (Orthanc REST),
# 4242 (Orthanc DICOM), 8526 (OHIF), 19530 (Milvus gRPC),
# 9091 (Milvus metrics), 8520-8532 (NIMs)
```

For production, bind only to localhost or specific interfaces:

```yaml
# Restrict to localhost only
imaging-api:
  ports:
    - "127.0.0.1:8524:8524"
```

### 10.6 Request Size Limits

The FastAPI server enforces a maximum request body size:

```bash
IMAGING_MAX_REQUEST_SIZE_MB=10
```

Adjust upward if uploading large DICOM files directly through the API. For large
studies, use Orthanc's dedicated upload endpoints instead.

### 10.7 TLS Termination

For production, place a reverse proxy (nginx, Traefik, Caddy) in front of the
services with TLS:

```nginx
# /etc/nginx/sites-available/imaging-agent
server {
    listen 443 ssl http2;
    server_name imaging.yourdomain.com;

    ssl_certificate /etc/letsencrypt/live/imaging.yourdomain.com/fullchain.pem;
    ssl_certificate_key /etc/letsencrypt/live/imaging.yourdomain.com/privkey.pem;

    # FastAPI
    location /api/ {
        proxy_pass http://127.0.0.1:8524/;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    # Streamlit
    location / {
        proxy_pass http://127.0.0.1:8525/;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
        proxy_set_header Host $host;
    }

    # OHIF
    location /viewer/ {
        proxy_pass http://127.0.0.1:8526/;
    }
}
```

---

## 11. Monitoring and Health Checks

### 11.1 Health Endpoints

| Endpoint | Port | Description |
|----------|------|-------------|
| `GET /health` | 8524 | Full health report: Milvus, collections, NIM status |
| `GET /healthz` | 8524 | Simple liveness probe (returns 200 OK) |
| `GET /metrics` | 8524 | Prometheus-compatible metrics |
| `GET /_stcore/health` | 8525 | Streamlit internal health check |
| `GET /healthz` | 9091 | Milvus liveness check |
| `GET /metrics` | 9091 | Milvus Prometheus metrics |
| `GET /v1/health/ready` | 8520 | NIM LLM readiness |
| `GET /v1/health/ready` | 8530 | NIM VISTA-3D readiness |
| `GET /v1/health/ready` | 8531 | NIM MAISI readiness |
| `GET /v1/health/ready` | 8532 | NIM VILA-M3 readiness |
| `GET /system` | 8042 | Orthanc system info |

### 11.2 Health Check Script

```bash
#!/usr/bin/env bash
# imaging-health-check.sh — Quick health check for all services

set -e

echo "=== Imaging Intelligence Agent Health Check ==="
echo ""

# FastAPI
echo -n "FastAPI (8524):     "
curl -sf http://localhost:8524/healthz > /dev/null && echo "OK" || echo "FAIL"

# Streamlit
echo -n "Streamlit (8525):   "
curl -sf http://localhost:8525/_stcore/health > /dev/null && echo "OK" || echo "FAIL"

# Milvus
echo -n "Milvus (9091):      "
curl -sf http://localhost:9091/healthz > /dev/null && echo "OK" || echo "FAIL"

# Orthanc
echo -n "Orthanc (8042):     "
curl -sf http://localhost:8042/system > /dev/null && echo "OK" || echo "FAIL"

# OHIF
echo -n "OHIF (8526):        "
curl -sf http://localhost:8526 > /dev/null && echo "OK" || echo "FAIL"

# NIM services
for svc in "LLM:8520" "VISTA-3D:8530" "MAISI:8531" "VILA-M3:8532"; do
  name="${svc%%:*}"
  port="${svc##*:}"
  echo -n "NIM $name ($port):  "
  curl -sf "http://localhost:$port/v1/health/ready" > /dev/null 2>&1 && echo "OK" || echo "FAIL"
done

echo ""
echo "=== Done ==="
```

### 11.3 Docker Health Check Configuration

Each service in docker-compose.yml has a built-in health check. View status:

```bash
# Show health status for all containers
docker compose ps

# Detailed health info for a specific container
docker inspect --format='{{json .State.Health}}' imaging-api | python3 -m json.tool

# Watch health transitions
docker events --filter type=container --filter event=health_status
```

### 11.4 Prometheus Metrics

The FastAPI server exposes Prometheus metrics at `/metrics` (port 8524):

```bash
curl -s http://localhost:8524/metrics
```

Metrics include:
- `imaging_rag_queries_total` — Total RAG queries processed
- `imaging_rag_query_duration_seconds` — RAG query latency histogram
- `imaging_nim_requests_total` — NIM requests by service and status
- `imaging_nim_request_duration_seconds` — NIM request latency
- `imaging_collection_search_total` — Searches per collection
- `imaging_workflow_runs_total` — Workflow executions by name and status

### 11.5 Prometheus Scrape Configuration

Add to your `prometheus.yml`:

```yaml
scrape_configs:
  - job_name: 'imaging-agent-api'
    scrape_interval: 15s
    static_configs:
      - targets: ['localhost:8524']
        labels:
          service: 'imaging-api'

  - job_name: 'imaging-milvus'
    scrape_interval: 30s
    static_configs:
      - targets: ['localhost:9091']
        labels:
          service: 'imaging-milvus'
```

### 11.6 Grafana Dashboard

Import the Prometheus data source and create panels for:

- RAG query latency (p50, p95, p99)
- NIM service availability (up/down per service)
- Collection search distribution
- Workflow execution success rate
- Milvus memory and segment count

### 11.7 Log Monitoring

```bash
# Follow all service logs
docker compose logs -f

# Follow specific service
docker compose logs -f imaging-api
docker compose logs -f imaging-streamlit

# Filter by log level (loguru format)
docker compose logs imaging-api 2>&1 | grep "ERROR"
docker compose logs imaging-api 2>&1 | grep "WARNING"

# Export logs to file
docker compose logs --no-color > imaging-agent-logs-$(date +%Y%m%d).txt
```

---

## 12. Troubleshooting

### 12.1 Setup Container Exits with Non-Zero Code

**Symptom:** `imaging-setup exited with code 1`

**Diagnosis:**

```bash
docker compose logs imaging-setup
```

**Common causes:**

1. **Milvus not ready:** The setup container depends on Milvus health check, but
   Milvus may still be initializing indexes.

   ```bash
   # Wait for Milvus and re-run setup
   docker compose restart imaging-setup
   ```

2. **Embedding model download failed:** First run downloads BGE-small-en-v1.5 (~130 MB)
   from HuggingFace. If network is unrestricted, retry:

   ```bash
   docker compose restart imaging-setup
   ```

3. **Out of memory:** Milvus requires at least 4 GB RAM. Check:

   ```bash
   docker stats --no-stream
   ```

### 12.2 Milvus Fails to Start

**Symptom:** `imaging-milvus-standalone` stays unhealthy.

```bash
# Check Milvus logs
docker compose logs milvus-standalone

# Check etcd
docker compose logs milvus-etcd

# Check MinIO
docker compose logs milvus-minio
```

**Common causes:**

1. **etcd quota exceeded:** Reset etcd:
   ```bash
   docker compose down
   docker volume rm imaging_intelligence_agent_etcd_data
   docker compose up -d
   ```

2. **Port conflict:** Another service using 19530 or 9091:
   ```bash
   sudo lsof -i :19530
   sudo lsof -i :9091
   ```

3. **Disk space:** Milvus needs free disk for indexes:
   ```bash
   df -h
   docker system df
   ```

### 12.3 NIM Container Fails to Start

**Symptom:** NIM container restarts repeatedly or stays unhealthy.

```bash
docker compose logs nim-llm
```

**Common causes:**

1. **GPU not available:**
   ```bash
   nvidia-smi
   docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
   ```

2. **Insufficient GPU memory:** Each NIM needs 8-16 GB VRAM. Check usage:
   ```bash
   nvidia-smi --query-gpu=memory.used,memory.total --format=csv
   ```

3. **NGC authentication failed:**
   ```bash
   echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
   ```

4. **Image not pulled:** NIM images are large. Verify:
   ```bash
   docker images | grep nim
   ```

### 12.4 Streamlit UI Not Loading

**Symptom:** Browser shows connection refused at `http://localhost:8525`.

```bash
docker compose logs imaging-streamlit
docker compose ps imaging-streamlit
```

**Common causes:**

1. **Container not running:**
   ```bash
   docker compose up -d imaging-streamlit
   ```

2. **Waiting for Milvus:** Streamlit depends on Milvus health check:
   ```bash
   curl -s http://localhost:9091/healthz
   ```

3. **Port conflict:**
   ```bash
   sudo lsof -i :8525
   ```

### 12.5 FastAPI Returns 500 Errors

**Symptom:** API calls return HTTP 500.

```bash
# Check FastAPI logs
docker compose logs imaging-api 2>&1 | tail -50

# Test health endpoint
curl -s http://localhost:8524/health | python3 -m json.tool
```

**Common causes:**

1. **Anthropic API key missing or invalid:**
   ```bash
   # Verify key is set
   docker exec imaging-api env | grep ANTHROPIC
   ```

2. **Milvus connection failed:**
   ```bash
   docker exec imaging-api python -c "
   from pymilvus import connections
   connections.connect(host='milvus-standalone', port=19530)
   print('Connected')
   "
   ```

3. **Collections not created:** Re-run setup:
   ```bash
   docker compose restart imaging-setup
   docker compose logs -f imaging-setup
   ```

### 12.6 OHIF Viewer Shows No Studies

**Symptom:** OHIF loads but study list is empty.

1. Verify Orthanc has studies:
   ```bash
   curl -s http://localhost:8042/studies -u admin:orthanc | python3 -m json.tool
   ```

2. Verify DICOMweb is enabled:
   ```bash
   curl -s http://localhost:8042/dicom-web/studies -u admin:orthanc | python3 -m json.tool
   ```

3. Check OHIF config is mounted:
   ```bash
   docker exec imaging-ohif-viewer cat /usr/share/nginx/html/app-config.js
   ```

4. Check browser console for CORS errors. If accessing from a different host, ensure
   Orthanc allows the origin.

### 12.7 Cross-Modal Genomics Not Working

**Symptom:** Workflow results do not include genomic evidence enrichment.

1. Verify cross-modal is enabled:
   ```bash
   IMAGING_CROSS_MODAL_ENABLED=True
   ```

2. Verify `genomic_evidence` collection exists and has data:
   ```bash
   curl -s http://localhost:8524/collections | python3 -m json.tool
   # Look for genomic_evidence with non-zero count
   ```

3. The `genomic_evidence` collection must be pre-populated by the Stage 2 RAG
   pipeline. It is read-only from the imaging agent's perspective.

### 12.8 Embedding Model Download Fails

**Symptom:** Containers fail with HuggingFace download errors.

The BGE-small-en-v1.5 model (~130 MB) is downloaded on first run. If behind a
proxy or firewall:

```bash
# Pre-download the model
docker exec imaging-api python -c "
from sentence_transformers import SentenceTransformer
model = SentenceTransformer('BAAI/bge-small-en-v1.5')
print('Model loaded:', model.get_sentence_embedding_dimension())
"
```

Or set the HuggingFace cache directory and pre-populate it:

```bash
# On host, download model
pip install sentence-transformers
python -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('BAAI/bge-small-en-v1.5')"

# Mount host cache into container
# Add to docker-compose.yml:
# volumes:
#   - ~/.cache/huggingface:/home/imaginguser/.cache/huggingface:ro
```

### 12.9 Port Conflicts

If other HCLS AI Factory agents or services are using the same ports:

```bash
# Find what is using a port
sudo lsof -i :8524
sudo ss -tlnp | grep 8524

# Use the DGX Spark port mapping (Section 5) or override in .env:
IMAGING_API_PORT=8624
IMAGING_STREAMLIT_PORT=8625
```

### 12.10 Docker Compose Version Issues

If you see `version is obsolete` warnings:

```bash
# The lite compose file omits version (correct for Compose v2)
# The full compose file includes version: "3.8" for backward compatibility
# Both work with Docker Compose v2.20+

# Check your Docker Compose version
docker compose version
```

### 12.11 Running Tests

The project includes 620 tests:

```bash
# Run all tests
python3 -m pytest tests/ -v

# Run tests with coverage
python3 -m pytest tests/ -v --cov=src --cov=api --cov-report=term-missing

# Run specific test modules
python3 -m pytest tests/test_rag_engine.py -v
python3 -m pytest tests/test_nim_clients.py -v
python3 -m pytest tests/test_workflows.py -v

# Run end-to-end validation
python scripts/validate_e2e.py --quick
```

### 12.12 Complete Reset

To start completely fresh:

```bash
# Stop all containers and remove volumes, networks, and images
docker compose -f docker-compose.yml down -v --rmi all
docker compose -f docker-compose.lite.yml down -v --rmi all

# Remove orphan volumes
docker volume prune -f

# Remove the network
docker network rm imaging-network 2>/dev/null || true

# Rebuild and restart
docker compose -f docker-compose.lite.yml up -d --build
```

### 12.13 Getting Help

1. Check service logs: `docker compose logs <service-name>`
2. Check Docker events: `docker events --filter type=container`
3. Check resource usage: `docker stats --no-stream`
4. Review `docs/NIM_INTEGRATION_GUIDE.md` for NIM-specific issues
5. Review `docs/ARCHITECTURE_GUIDE.md` for system design details
6. Run the end-to-end validator: `python scripts/validate_e2e.py --quick`

---

## Appendix A: Docker Volumes

| Volume | Service | Contents |
|--------|---------|----------|
| `orthanc_data` | Orthanc | DICOM study storage (SQLite + files) |
| `etcd_data` | etcd | Milvus metadata key-value store |
| `minio_data` | MinIO | Milvus segment and index files |
| `milvus_data` | Milvus | Vector data and WAL |
| `nim_models` | NIM containers | Downloaded model weights |

## Appendix B: Network Configuration

| Network | Driver | Purpose |
|---------|--------|---------|
| `imaging-network` | bridge | All inter-service communication |

All services join `imaging-network`. Docker DNS resolves service names
(e.g., `milvus-standalone`, `nim-llm`, `orthanc`) to container IPs.

## Appendix C: Port Quick Reference

| Port | Service | Protocol |
|------|---------|----------|
| 4242 | Orthanc DICOM C-STORE | DICOM |
| 8042 | Orthanc REST API | HTTP |
| 8520 | NIM LLM | HTTP |
| 8524 | FastAPI REST Server | HTTP |
| 8525 | Streamlit Chat UI | HTTP |
| 8526 | OHIF Viewer | HTTP |
| 8530 | NIM VISTA-3D | HTTP |
| 8531 | NIM MAISI | HTTP |
| 8532 | NIM VILA-M3 | HTTP |
| 9091 | Milvus Metrics | HTTP |
| 19530 | Milvus gRPC | gRPC |

**DGX Spark external port overrides:**

| Internal | External | Service |
|----------|----------|---------|
| 8525 | 8505 | Streamlit UI |
| 8524 | 8105 | FastAPI Server |

---

*Last updated: March 2026*
*Imaging Intelligence Agent v1.0.0 — HCLS AI Factory*
