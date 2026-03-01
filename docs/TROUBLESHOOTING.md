---
tags:
  - Troubleshooting
  - Operations
  - Support
---

# Troubleshooting Guide

Common issues and solutions for the HCLS AI Factory platform, organized by component.

---

## Quick Diagnostics

Before diving into specific issues, run these checks:

```bash
# 1. Check all service health from landing page
curl -s http://localhost:8080/api/check-services | python3 -m json.tool

# 2. Check Docker containers
docker compose ps

# 3. Check GPU availability
nvidia-smi

# 4. Check disk space
df -h /

# 5. Check memory
free -h
```

---

## Services Will Not Start

### Docker Compose fails with port conflicts

**Symptom:** `Bind for 0.0.0.0:8501 failed: port is already allocated`

**Solution:**
```bash
# Find what is using the port
lsof -i :8501

# Stop the conflicting process, or change the port in docker-compose.yml
```

### NVIDIA Container Runtime not found

**Symptom:** `docker: Error response from daemon: Unknown runtime specified nvidia`

**Solution:**
```bash
# Install NVIDIA Container Toolkit
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
distribution=$(. /etc/os-release; echo $ID$VERSION_ID)
curl -s -L "https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list" | \
  sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
  sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker
```

### Services start but show "unhealthy"

**Symptom:** Landing page shows services as red/offline despite containers running.

**Solution:**
```bash
# Check container logs for the specific service
docker compose logs <service-name> --tail=50

# Common causes:
# 1. Milvus not ready yet (takes 30-60s to initialize)
# 2. Missing API keys (check .env file)
# 3. Missing data files (run setup-data.sh)
```

---

## Stage 1: Genomics Pipeline

### Parabricks license error

**Symptom:** `NVIDIA Parabricks license validation failed`

**Solution:**
```bash
# Ensure NGC API key is set
echo $NGC_CLI_API_KEY

# Re-authenticate with NGC
docker login nvcr.io -u '$oauthtoken' -p $NGC_CLI_API_KEY
```

### Out of GPU memory during alignment

**Symptom:** `CUDA out of memory` during BWA-MEM2 or DeepVariant

**Solution:**

- Reduce batch size in Parabricks configuration
- For DGX Spark (128GB unified memory), this should not occur with default settings
- For GPUs with less VRAM, use the `--low-memory` flag if available, or process chromosomes individually

### FASTQ files not found

**Symptom:** `FileNotFoundError` when starting genomics pipeline

**Solution:**
```bash
# Verify data download completed
ls -la data/genomics/

# If missing, re-run data setup for genomics only
./setup-data.sh --genomics
```

---

## Stage 2: RAG/Chat Pipeline

### Milvus connection refused

**Symptom:** `ConnectionRefusedError: [Errno 111] Connection refused` on port 19530

**Solution:**
```bash
# Check if Milvus container is running
docker compose ps milvus

# Check Milvus logs
docker compose logs milvus --tail=50

# Restart Milvus
docker compose restart milvus

# Wait for initialization (30-60 seconds)
sleep 30
curl -s http://localhost:19530/v1/vector/collections
```

### Empty search results

**Symptom:** Queries return no results or empty evidence

**Solution:**
```bash
# Check collection counts via Attu or API
curl -s http://localhost:19530/v1/vector/collections

# If collections are empty, re-run data ingestion
cd rag-chat-pipeline
python3 ingest.py --all
```

### Claude API errors

**Symptom:** `AuthenticationError` or `RateLimitError` from Anthropic API

**Solution:**
```bash
# Verify API key is set
echo $ANTHROPIC_API_KEY

# Check API key validity
curl -s https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "content-type: application/json" \
  -H "anthropic-version: 2023-06-01" \
  -d '{"model":"claude-sonnet-4-20250514","max_tokens":10,"messages":[{"role":"user","content":"test"}]}'

# For rate limits: the pipeline automatically retries with exponential backoff
# If persistent, check your Anthropic plan limits
```

### Embedding model download fails

**Symptom:** `OSError` or timeout when loading `BAAI/bge-small-en-v1.5`

**Solution:**
```bash
# Pre-download the model
python3 -c "from sentence_transformers import SentenceTransformer; SentenceTransformer('BAAI/bge-small-en-v1.5')"

# If behind a firewall, set HuggingFace cache directory
export HF_HOME=/path/to/cache
export TRANSFORMERS_CACHE=/path/to/cache
```

---

## Stage 3: Drug Discovery Pipeline

### BioNeMo NIM services unavailable

**Symptom:** Connection refused on ports 8001 (MolMIM) or 8002 (DiffDock)

**Solution:**
```bash
# Check NIM mode
echo $NIM_MODE

# For cloud NIMs (recommended for DGX Spark ARM64):
export NIM_MODE=cloud
export NGC_CLI_API_KEY=your_key_here

# Verify cloud NIM access
curl -s https://health.api.nvidia.com/v1/health/ready

# For local NIMs (x86 only):
docker compose up -d molmim diffdock
```

### DiffDock docking fails

**Symptom:** `Error during docking` or empty results from DiffDock

**Solution:**

- Verify the PDB file is valid and contains the target protein
- Check that the ligand SMILES string is valid: `python3 -c "from rdkit import Chem; print(Chem.MolFromSmiles('your_smiles') is not None)"`
- For cloud NIM: verify NVCF asset staging completed (check logs for asset ID)
- Try reducing the number of docking poses

### RDKit import errors

**Symptom:** `ModuleNotFoundError: No module named 'rdkit'`

**Solution:**
```bash
# Install RDKit via conda (recommended)
conda install -c conda-forge rdkit

# Or via pip
pip install rdkit-pypi
```

---

## Intelligence Agents

### Agent UI not loading

**Symptom:** Streamlit UI returns connection error on agent port (8521, 8525, or 8526)

**Solution:**
```bash
# Check if the agent container is running
docker compose ps | grep agent

# Start the specific agent
docker compose up -d cart-agent    # Port 8521
docker compose up -d imaging-agent # Port 8525
docker compose up -d onco-agent   # Port 8526

# Check agent logs
docker compose logs cart-agent --tail=50
```

### Agent cannot connect to Milvus

**Symptom:** Agent health check shows `milvus: false`

**Solution:**
```bash
# Agents share the same Milvus instance as the core platform
# Verify Milvus is running and accessible
curl -s http://localhost:19530/v1/vector/collections

# Check the agent's Milvus host/port configuration
# Default: MILVUS_HOST=localhost, MILVUS_PORT=19530
# In Docker: MILVUS_HOST=milvus (container name)
```

### Cross-modal triggers not firing

**Symptom:** Agent queries do not pull evidence from shared genomic collections

**Solution:**

- Verify the shared `genomic_evidence` collection exists in Milvus
- Check that the cross-modal threshold is not set too high (default: 0.7)
- Ensure the agent has read access to shared collections

---

## Landing Page

### Landing page shows all services as offline

**Symptom:** All service tiles are red despite services running

**Solution:**
```bash
# Check if the landing page can reach services
# Services must be accessible from the landing page container/process
# If running in Docker, ensure services are on the same network

docker network ls
docker network inspect hcls-ai-factory_default
```

### Auto-start not working

**Symptom:** Genomics and RAG services do not auto-start from the landing page

**Solution:**

- Auto-start only works when the landing page runs on the same host as the services
- Check that the service directories exist at the expected paths
- Verify Python virtual environments are set up for each service

---

## Nextflow Orchestrator

### Nextflow cgroup errors on DGX Spark

**Symptom:** `Cannot get cgroup` or process resource errors

**Solution:**
```bash
# Use the Python orchestrator as an alternative
cd hls-orchestrator
python3 run_pipeline.py --mode demo

# Or run Nextflow with NXF_OPTS to bypass cgroup
export NXF_OPTS="-Xms512m -Xmx4g"
./nextflow run main.nf -profile dgx_spark --mode demo
```

### Pipeline hangs at a stage

**Symptom:** Nextflow shows a process running but no progress

**Solution:**
```bash
# Check Nextflow work directory for logs
ls -la work/

# Find the specific task directory
find work/ -name ".command.log" -newer work/ -exec tail -20 {} \;

# Resume from the last successful stage
./nextflow run main.nf -resume
```

---

## Data Setup

### setup-data.sh download failures

**Symptom:** Downloads fail or stall during `./setup-data.sh --all`

**Solution:**
```bash
# The script supports automatic retry — re-run safely
./setup-data.sh --all

# Download specific stages only
./setup-data.sh --genomics    # Reference genome, FASTQ files (~400GB)
./setup-data.sh --rag         # ClinVar, AlphaMissense, knowledge base (~2GB)
./setup-data.sh --drug        # PDB structures, seed compounds (~100MB)

# Verify checksums after download
./setup-data.sh --verify
```

### Insufficient disk space

**Symptom:** `No space left on device` during data download or pipeline execution

**Solution:**

| Component | Approximate Size |
|---|---|
| Reference genome (GRCh38) | 3.1 GB |
| FASTQ sequencing data | ~200 GB |
| ClinVar + AlphaMissense | ~2 GB |
| Milvus vector database | ~15 GB |
| Pipeline outputs (BAM, VCF) | ~120 GB |
| Docker images | ~30 GB |
| **Total recommended** | **500 GB minimum** |

```bash
# Check current usage
du -sh data/ genomics-pipeline/data/ rag-chat-pipeline/data/

# Clean up old pipeline outputs
rm -rf results/old_run_*/
docker system prune -f
```

---

## Monitoring

### Grafana dashboards empty

**Symptom:** Grafana loads but shows "No data" in panels

**Solution:**
```bash
# Check Prometheus is scraping targets
curl -s http://localhost:9099/api/v1/targets | python3 -m json.tool

# Verify Node Exporter is running
curl -s http://localhost:9100/metrics | head -5

# Verify DCGM Exporter is running (GPU metrics)
curl -s http://localhost:9400/metrics | head -5

# Default Grafana credentials: admin / admin
```

### Prometheus alerts firing

**Symptom:** Alert manager notifications for service down or GPU errors

**Solution:**
```bash
# Check active alerts
curl -s http://localhost:9099/api/v1/alerts | python3 -m json.tool

# Common alerts and resolutions:
# - ServiceDown: restart the affected service
# - GPUHighTemp: check GPU cooling, reduce workload
# - HighMemoryUsage: check for memory leaks, restart services
# - MilvusUnhealthy: restart Milvus container
```

---

## Network and Security

### CORS errors in browser

**Symptom:** Browser console shows `Access-Control-Allow-Origin` errors

**Solution:**

- Add your origin to the `CORS_ORIGINS` environment variable
- Default allows `http://localhost:*` patterns
- For production, set specific origins instead of wildcards

### API requests rejected (413)

**Symptom:** `Request body exceeds X MB limit`

**Solution:**

- Default limit is 50 MB per request
- Increase `MAX_REQUEST_SIZE_MB` in the service configuration
- For large file uploads (VCF, FASTQ), use the dedicated file upload endpoints

---

## Getting Help

If your issue is not covered here:

1. Check the [Deployment Guide](HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md) for detailed configuration
2. Review service logs: `docker compose logs <service> --tail=100`
3. Open an issue on [GitHub](https://github.com/ajones1923/hcls-ai-factory/issues) with:
   - Steps to reproduce
   - Relevant log output
   - Hardware and OS details
   - Docker and NVIDIA driver versions
