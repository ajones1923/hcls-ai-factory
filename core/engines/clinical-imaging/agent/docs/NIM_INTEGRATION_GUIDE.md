# Imaging Intelligence Agent -- NIM Integration Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0

---

## 1. What Are NVIDIA NIMs

NVIDIA NIM (NVIDIA Inference Microservice) packages are pre-built, optimized containers that expose AI models as REST API endpoints. Each NIM runs as a standalone Docker container with GPU acceleration, providing production-ready inference with health checks, auto-scaling support, and OpenAI-compatible APIs where applicable.

The Imaging Intelligence Agent integrates four NIMs:

| NIM | Purpose | Image |
|---|---|---|
| **Llama-3 8B Instruct** | Text generation, clinical reasoning | `nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest` |
| **VISTA-3D** | 3D medical image segmentation | `nvcr.io/nvidia/nim/vista3d:latest` |
| **MAISI** | Synthetic CT volume generation | `nvcr.io/nvidia/nim/maisi:latest` |
| **VILA-M3** | Vision-language medical image understanding | `nvcr.io/nvidia/nim/vilam3:latest` |

---

## 2. Prerequisites

### 2.1 NGC Account and API Key

1. Create a free account at [https://ngc.nvidia.com](https://ngc.nvidia.com)
2. Navigate to Setup > API Key
3. Generate a new API key
4. Save the key -- you will need it for container pulls and runtime authentication

### 2.2 Docker with NVIDIA Runtime

```bash
# Verify NVIDIA driver
nvidia-smi

# Verify NVIDIA Container Toolkit
docker run --rm --gpus all nvidia/cuda:12.0-base nvidia-smi

# If not installed, follow:
# https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html
```

### 2.3 NGC Docker Registry Login

```bash
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
```

---

## 3. GPU Memory Planning for DGX Spark

The NVIDIA DGX Spark provides 128 GB unified LPDDR5x memory shared between the Grace CPU and GB10 GPU via NVLink-C2C. This enables running all four NIMs simultaneously.

### Memory Budget

| Service | Estimated GPU Memory | Notes |
|---|---|---|
| **Llama-3 8B Instruct** | ~16 GB | FP16 weights + KV cache |
| **VISTA-3D** | ~8 GB | 3D segmentation model |
| **MAISI** | ~12 GB | Latent diffusion model |
| **VILA-M3** | ~16 GB | Vision-language model |
| **Milvus** | ~2 GB | Vector index (CPU, but shared memory) |
| **Application stack** | ~4 GB | FastAPI, Streamlit, embeddings |
| **Total estimated** | **~58 GB** | Well within 128 GB budget |
| **Headroom** | ~70 GB | Available for inference batch size, data loading |

### Running All NIMs Simultaneously

With 128 GB unified memory, DGX Spark can run all four NIMs plus the full application stack with significant headroom. This is not possible on typical workstation GPUs (24-48 GB).

### Reduced Configurations

If running on limited hardware:

| Configuration | NIMs | GPU Memory Required |
|---|---|---|
| Full | All 4 NIMs | ~52 GB |
| LLM + VISTA-3D | Llama-3 + VISTA-3D | ~24 GB |
| LLM only | Llama-3 | ~16 GB |
| Mock mode | None | 0 GB (CPU only) |

---

## 4. Pulling NIM Containers

### 4.1 Authenticate with NGC

```bash
export NGC_API_KEY=your_ngc_api_key_here
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin
```

### 4.2 Pull Individual NIMs

```bash
# Llama-3 8B Instruct
docker pull nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest

# VISTA-3D (3D medical image segmentation)
docker pull nvcr.io/nvidia/nim/vista3d:latest

# MAISI (medical AI synthetic imaging)
docker pull nvcr.io/nvidia/nim/maisi:latest

# VILA-M3 (vision-language model)
docker pull nvcr.io/nvidia/nim/vilam3:latest
```

### 4.3 Pull via Docker Compose

```bash
# Pulls all images defined in docker-compose.yml
docker compose pull
```

---

## 5. Configuration via .env File

Copy the example environment file and configure:

```bash
cp .env.example .env
```

### Key NIM Settings

```bash
# NVIDIA NGC API Key (required for NIM containers)
NGC_API_KEY=your_ngc_api_key_here

# NIM Mode: "local" uses real NIM containers, "mock" uses synthetic responses
IMAGING_NIM_MODE=local

# NIM Service URLs (defaults match docker-compose.yml port mappings)
IMAGING_NIM_LLM_URL=http://localhost:8520/v1
IMAGING_NIM_VISTA3D_URL=http://localhost:8530
IMAGING_NIM_MAISI_URL=http://localhost:8531
IMAGING_NIM_VILAM3_URL=http://localhost:8532

# Allow automatic fallback to mock when a NIM service is unavailable
IMAGING_NIM_ALLOW_MOCK_FALLBACK=true
```

### Mode Behavior

| `NIM_MODE` | `ALLOW_MOCK_FALLBACK` | Behavior |
|---|---|---|
| `local` | `true` | Try real NIM, fall back to mock on failure |
| `local` | `false` | Try real NIM, raise error on failure |
| `cloud` | `true` | Uses NVIDIA Cloud NIM endpoints (integrate.api.nvidia.com). Requires NVIDIA_API_KEY. LLM uses meta/llama-3.1-8b-instruct, VLM uses meta/llama-3.2-11b-vision-instruct. VISTA-3D and MAISI fall back to mock. |
| `cloud` | `false` | Same as above, but raises error if cloud endpoint is unreachable |
| `mock` | (ignored) | Always use mock responses |

### Feature Comparison by Mode

| Feature | local | cloud | mock |
|---|---|---|---|
| LLM inference | Local Llama-3 NIM (port 8520) | NVIDIA Cloud (llama-3.1-8b) | Template responses |
| VLM inference | Local VILA-M3 NIM (port 8532) | NVIDIA Cloud (llama-3.2-11b-vision) | Template responses |
| VISTA-3D segmentation | Local NIM (port 8530) | Mock fallback | Synthetic masks |
| MAISI generation | Local NIM (port 8531) | Mock fallback | Metadata only |
| GPU required | Yes (all 4 NIMs) | No (API only) | No |
| API key required | NGC_API_KEY | NVIDIA_API_KEY | None |

---

## NVIDIA Cloud NIM Endpoints

### Cloud NIM Configuration

When local NIM containers are unavailable or for development without GPU, the agent can use NVIDIA Cloud NIM endpoints. This requires an NVIDIA API key.

**Environment variables:**
```bash
IMAGING_NIM_MODE=cloud
IMAGING_NVIDIA_API_KEY=nvapi-xxxxx  # or NVIDIA_API_KEY
```

**Settings in `config/settings.py`:**
```python
NIM_CLOUD_URL: str = "https://integrate.api.nvidia.com/v1"
NIM_CLOUD_LLM_MODEL: str = "meta/llama-3.1-8b-instruct"
NIM_CLOUD_VLM_MODEL: str = "meta/llama-3.2-11b-vision-instruct"
```

**Behavior:**
- LLM requests route to `meta/llama-3.1-8b-instruct` via NVIDIA Cloud
- VLM requests route to `meta/llama-3.2-11b-vision-instruct` via NVIDIA Cloud
- VISTA-3D and MAISI remain in mock mode (no cloud equivalent -- these require local GPU)
- Cloud mode provides real LLM/VLM inference without local GPU

**Health check with cloud NIMs:**
```bash
curl -s http://localhost:8524/nim/status | python3 -m json.tool
# Expected:
# vista3d: mock, maisi: mock, vila_m3: cloud, llm: cloud
```

---

## 6. Running NIM Containers

### 6.1 Via Docker Compose (Recommended)

```bash
# Full stack with all NIMs
docker compose up -d

# Check NIM startup progress (NIMs take 60-120 seconds to load models)
docker compose logs -f nim-llm
docker compose logs -f nim-vista3d
docker compose logs -f nim-maisi
docker compose logs -f nim-vilam3
```

### 6.2 Running NIMs Individually

```bash
# Llama-3 8B Instruct
docker run -d --name nim-llm \
  --gpus all \
  -e NGC_API_KEY=$NGC_API_KEY \
  -p 8520:8000 \
  nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest

# VISTA-3D
docker run -d --name nim-vista3d \
  --gpus all \
  -e NGC_API_KEY=$NGC_API_KEY \
  -p 8530:8000 \
  nvcr.io/nvidia/nim/vista3d:latest

# MAISI
docker run -d --name nim-maisi \
  --gpus all \
  -e NGC_API_KEY=$NGC_API_KEY \
  -p 8531:8000 \
  nvcr.io/nvidia/nim/maisi:latest

# VILA-M3
docker run -d --name nim-vilam3 \
  --gpus all \
  -e NGC_API_KEY=$NGC_API_KEY \
  -p 8532:8000 \
  nvcr.io/nvidia/nim/vilam3:latest
```

### 6.3 Health Check

Each NIM exposes a `/v1/health/ready` endpoint:

```bash
# Check if Llama-3 NIM is ready
curl http://localhost:8520/v1/health/ready

# Check VISTA-3D
curl http://localhost:8530/v1/health/ready

# Check MAISI
curl http://localhost:8531/v1/health/ready

# Check VILA-M3
curl http://localhost:8532/v1/health/ready
```

Response: HTTP 200 when ready, 503 when still loading.

---

## 7. Mock Mode for Development

Mock mode enables full development and testing without any GPU or NIM containers. All NIM-dependent features return clinically realistic synthetic responses.

### Enabling Mock Mode

**Option 1: Environment variable**
```bash
export IMAGING_NIM_MODE=mock
```

**Option 2: In .env file**
```bash
IMAGING_NIM_MODE=mock
```

**Option 3: Docker Compose Lite (pre-configured)**
```bash
docker compose -f docker-compose.lite.yml up -d
```

### Mock Response Examples

**VISTA-3D mock** returns synthetic segmentation with realistic anatomical classes and volumes:
```json
{
  "classes_detected": ["liver", "spleen", "left_kidney", "right_kidney"],
  "volumes": {"liver": 1450.2, "spleen": 180.5, "left_kidney": 160.3, "right_kidney": 155.8},
  "num_classes": 4,
  "inference_time_ms": 45.2,
  "model": "vista3d",
  "is_mock": true
}
```

**MAISI mock** returns synthetic CT generation metadata:
```json
{
  "resolution": "512x512x512",
  "body_region": "chest",
  "num_classes_annotated": 104,
  "generation_time_ms": 120.5,
  "model": "maisi",
  "is_mock": true
}
```

---

## 8. Testing NIM Connectivity

### 8.1 Via API Health Endpoint

```bash
curl http://localhost:8524/health | python -m json.tool
```

The `/health` response includes NIM service status:

```json
{
  "status": "healthy",
  "collections": {"imaging_literature": 5000, "...": "..."},
  "total_vectors": 12345,
  "nim_services": {
    "vista3d": "available",
    "maisi": "mock",
    "vila_m3": "available",
    "llm": "available"
  }
}
```

### 8.2 Via NIM Status Endpoint

```bash
curl http://localhost:8524/nim/status | python -m json.tool
```

Returns detailed status for each NIM:

```json
{
  "services": [
    {"name": "vista3d", "status": "available", "url": "http://localhost:8530"},
    {"name": "maisi", "status": "mock", "url": "http://localhost:8531"},
    {"name": "vila_m3", "status": "available", "url": "http://localhost:8532"},
    {"name": "llm", "status": "available", "url": "http://localhost:8520/v1"}
  ],
  "available_count": 3,
  "mock_count": 1,
  "unavailable_count": 0
}
```

### 8.3 Direct NIM Test

```bash
# Test Llama-3 LLM
curl http://localhost:8520/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "meta-llama3-8b-instruct",
    "messages": [{"role": "user", "content": "What is intracranial hemorrhage?"}],
    "max_tokens": 200
  }'
```

---

## 9. API Endpoints for NIM Services

### 9.1 VISTA-3D Segmentation

```bash
curl -X POST http://localhost:8524/nim/vista3d/segment \
  -H "Content-Type: application/json" \
  -d '{
    "input_path": "",
    "target_classes": ["liver", "spleen", "kidney"]
  }'
```

### 9.2 MAISI Synthetic CT Generation

```bash
curl -X POST http://localhost:8524/nim/maisi/generate \
  -H "Content-Type: application/json" \
  -d '{
    "body_region": "chest",
    "resolution": "512x512x512",
    "num_classes": 104
  }'
```

### 9.3 VILA-M3 Image Analysis

```bash
curl -X POST http://localhost:8524/nim/vilam3/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Describe the findings in this chest X-ray",
    "input_path": ""
  }'
```

---

## 10. Troubleshooting

### NIM container fails to start

```bash
# Check container logs
docker logs imaging-nim-llm

# Common issues:
# - NGC_API_KEY not set or invalid
# - Insufficient GPU memory
# - NVIDIA driver version mismatch
# - Container image not pulled (authentication failure)
```

### NIM shows "unavailable" in health check

```bash
# NIM containers take 60-120 seconds to load models after start
# Check if the container is still initializing:
docker logs --tail 50 imaging-nim-vista3d

# Verify the container is running:
docker ps | grep nim
```

### Mock fallback is not working

Ensure `IMAGING_NIM_ALLOW_MOCK_FALLBACK=true` in your `.env` file. If `NIM_MODE=local` and `ALLOW_MOCK_FALLBACK=false`, failures will raise errors instead of falling back.

### GPU out of memory

If running on limited GPU hardware:
1. Reduce the number of simultaneous NIMs
2. Set unused NIM URLs to empty strings
3. Use `IMAGING_NIM_MODE=mock` for NIMs you cannot run
4. On DGX Spark (128 GB), all four NIMs fit comfortably

### Container pull fails with 401

```bash
# Re-authenticate with NGC
echo "$NGC_API_KEY" | docker login nvcr.io --username '$oauthtoken' --password-stdin

# Verify API key is valid at https://ngc.nvidia.com/setup/api-key
```

### Port conflicts

Default NIM ports (8520, 8530, 8531, 8532) may conflict with other services. To change, edit the port mappings in `docker-compose.yml` and update the corresponding `IMAGING_NIM_*_URL` variables in `.env`.

---

*For complete architecture details, see `ARCHITECTURE_GUIDE.md`. For the full implementation specification, see `PROJECT_BIBLE.md`.*
