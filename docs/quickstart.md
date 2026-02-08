# Quick-Start Checklist

Deploy the HCLS AI Factory in 30 minutes. This checklist extracts the critical steps from the full [Deployment Guide](HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md).

---

## Prerequisites

- [ ] **Hardware**: NVIDIA DGX Spark (or equivalent with GB10 GPU, 128GB unified memory)
- [ ] **Storage**: 500GB+ free space for genomics data and models
- [ ] **Network**: Internet access for pulling containers and reference data
- [ ] **Software**: Docker 24+, Docker Compose 2.20+, Git
- [ ] **Download tools**: `aria2c`, `pigz` (`sudo apt-get install -y aria2 pigz`)

---

## Step 1: Clone the Repository

```bash
git clone https://github.com/NVIDIA/hcls-ai-factory.git
cd hcls-ai-factory
```

---

## Step 2: Environment Setup

```bash
# Copy environment template
cp .env.example .env

# Edit with your API keys
nano .env
```

**Required variables:**

| Variable | Description |
|----------|-------------|
| `ANTHROPIC_API_KEY` | Claude API key for RAG chat |
| `NGC_API_KEY` | NVIDIA NGC key for BioNeMo models |

---

## Step 3: Download Required Data

```bash
# Download all data (~500 GB, one-time)
./setup-data.sh --all

# Or download by stage
./setup-data.sh --stage2    # ClinVar + AlphaMissense (~2 GB, fast)
./setup-data.sh --stage1    # HG002 FASTQ + reference (~300 GB, 2-6 hours)

# Check status
./setup-data.sh --status
```

> **Note**: This is the most time-consuming step. See [DATA_SETUP.md](DATA_SETUP.md) for troubleshooting FASTQ checksum failures, disk space issues, and resuming interrupted downloads.

---

## Step 4: Start Core Services

```bash
# Start all services
docker compose up -d

# Verify services are running
docker compose ps
```

**Expected services:**

- `genomics-portal` — Parabricks + DeepVariant (port 5000)
- `rag-api` — RAG engine + Claude integration (port 5001)
- `streamlit-chat` — Chat UI (port 8501)
- `molmim` / `diffdock` — BioNeMo NIMs (ports 8001, 8002)
- `discovery-ui` — Drug discovery interface (port 8505)
- `milvus` / `etcd` / `minio` — Vector database stack (port 19530)
- `grafana` — Monitoring dashboard (port 3000)
- `landing-page` — Service health monitor (port 8080)

---

## Step 5: Verify GPU Access

```bash
# Check GPU visibility
docker compose exec genomics-portal nvidia-smi
```

You should see your GPU(s) listed with available memory.

---

## Step 6: Run a Test Pipeline

```bash
# Run the demo pipeline
python run_pipeline.py --mode demo

# Expected output: variant calls in output/demo/
```

---

## Step 7: Access the UI

Open your browser:

| Service | URL | Purpose |
|---------|-----|---------|
| Streamlit Chat | `http://localhost:8501` | Query variants with Claude |
| Grafana | `http://localhost:3000` | Monitor pipeline metrics |
| Landing Page | `http://localhost:8080` | Service health dashboard |
| RAG API | `http://localhost:5001` | REST API for variant queries |

---

## Troubleshooting

### Services won't start

```bash
# Check logs
docker compose logs -f

# Restart specific service
docker compose restart rag-api
```

### GPU not detected

```bash
# Verify NVIDIA runtime
docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
```

### Out of memory

Reduce batch sizes in `.env`:

```bash
PARABRICKS_BATCH_SIZE=1
MOLMIM_BATCH_SIZE=10
```

---

## Next Steps

- **Full deployment**: [Deployment Guide](HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md)
- **Run the demo**: [Demo Guide](HCLS_AI_FACTORY_DEMO_GUIDE.md)
- **Understand the architecture**: [White Paper](HCLS_AI_FACTORY_WHITE_PAPER_DGX_SPARK.md)

---

## Success Criteria

You're ready when:

- [ ] All Docker services show `Up` status
- [ ] GPU is visible in containers
- [ ] Streamlit chat responds to queries
- [ ] Grafana shows pipeline metrics

**Total time**: ~30 minutes (excluding data download — see Step 3 and [DATA_SETUP.md](DATA_SETUP.md))

---

*Need help? Open an issue on [GitHub](https://github.com/NVIDIA/hcls-ai-factory/issues).*
