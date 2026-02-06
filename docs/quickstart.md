# Quick-Start Checklist

Deploy the HCLS AI Factory in 30 minutes. This checklist extracts the critical steps from the full [Deployment Guide](HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md).

---

## Prerequisites

- [ ] **Hardware**: NVIDIA DGX Spark (or equivalent with GB10 GPU, 128GB unified memory)
- [ ] **Storage**: 500GB+ free space for genomics data and models
- [ ] **Network**: Internet access for pulling containers and reference data
- [ ] **Software**: Docker 24+, Docker Compose 2.20+, Git

---

## Step 1: Clone the Repository

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
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

## Step 3: Start Core Services

```bash
# Start all services
docker compose up -d

# Verify services are running
docker compose ps
```

**Expected services:**

- `genomics-pipeline` — Parabricks + DeepVariant
- `rag-service` — Milvus + Claude integration
- `drug-discovery` — BioNeMo MolMIM + DiffDock
- `grafana` — Monitoring dashboard (port 3000)
- `streamlit` — Chat UI (port 8501)

---

## Step 4: Verify GPU Access

```bash
# Check GPU visibility
docker compose exec genomics-pipeline nvidia-smi
```

You should see your GPU(s) listed with available memory.

---

## Step 5: Load Reference Data

```bash
# Download reference genome (one-time, ~15GB)
./scripts/download_reference.sh

# Load ClinVar and AlphaMissense annotations
./scripts/load_annotations.sh
```

---

## Step 6: Run a Test Pipeline

```bash
# Run with synthetic test data
./scripts/run_demo.sh --mock

# Expected output: variant calls in output/demo/
```

---

## Step 7: Access the UI

Open your browser:

| Service | URL | Purpose |
|---------|-----|---------|
| Streamlit Chat | `http://localhost:8501` | Query variants with Claude |
| Grafana | `http://localhost:3000` | Monitor pipeline metrics |
| API Docs | `http://localhost:8080/docs` | REST API reference |

---

## Troubleshooting

### Services won't start

```bash
# Check logs
docker compose logs -f

# Restart specific service
docker compose restart rag-service
```

### GPU not detected

```bash
# Verify NVIDIA runtime
docker run --rm --gpus all nvidia/cuda:12.0-base nvidia-smi
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

**Total time**: ~30 minutes (excluding reference data download)

---

*Need help? Open an issue on [GitHub](https://github.com/ajones1923/hcls-ai-factory/issues).*
