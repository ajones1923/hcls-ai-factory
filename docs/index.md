# HCLS AI Factory

**Transform patient DNA into therapeutic candidates in hours, not months.**

---

The Healthcare & Life Sciences (HCLS) AI Factory is an open-source precision medicine platform that unifies GPU-accelerated genomics, AI-powered evidence reasoning, and generative drug discovery into a single, continuous pipeline — running on a single NVIDIA DGX Spark workstation.

<div class="grid cards" markdown>

- :material-dna: **Stage 1 — Genomics**

    GPU-accelerated alignment and variant calling via NVIDIA Parabricks and DeepVariant. 200 GB of raw FASTQ data processed in under 1 hour.

    [:octicons-arrow-right-24: Genomics Pipeline](genomics-pipeline/README.md)

- :material-chat-processing: **Stage 2 — Evidence RAG**

    3.5 million variants annotated, embedded, and indexed in Milvus. Natural language queries grounded in ClinVar, AlphaMissense, and structural data via Claude.

    [:octicons-arrow-right-24: RAG & Chat Pipeline](rag-chat-pipeline/README.md)

- :material-molecule: **Stage 3 — Drug Discovery**

    BioNeMo MolMIM generates novel molecules. DiffDock predicts binding affinity. RDKit scores drug-likeness. 100 ranked candidates from a single target gene.

    [:octicons-arrow-right-24: Drug Discovery Pipeline](drug-discovery-pipeline/README.md)

</div>

---

## Key Numbers

| Metric | Value |
|---|---|
| Total variants called | ~11.7 million |
| High-quality variants (QUAL>30) | ~3.5 million |
| Genes covered | 201 across 13 therapeutic areas |
| Druggable targets | 171 (85%) |
| Drug candidates generated | 100 per target |
| End-to-end runtime | < 5 hours |
| Hardware | NVIDIA DGX Spark ($3,999) |

---

## Documentation

### Getting Started

- [**Demo Guide**](HCLS_AI_FACTORY_DEMO_GUIDE.md) — Step-by-step walkthrough for demonstrating the platform
- [**Deployment Guide**](HCLS_AI_FACTORY_DGX_SPARK_DEPLOYMENT_GUIDE.md) — Complete deployment and configuration guide for DGX Spark

### Architecture & Design

- [**Project Bible**](HCLS_AI_FACTORY_PROJECT_BIBLE.md) — Complete technical reference (scoring formulas, thresholds, configurations)
- [**White Paper**](HCLS_AI_FACTORY_WHITE_PAPER_DGX_SPARK.md) — Architecture overview and design rationale
- [**Architecture Mindmap**](HCLS_AI_Factory_Mindmap_Open.md) — Visual system map

### Learning

- [**Foundations Guide**](HCLS_AI_FACTORY_LEARNING_GUIDE_FOUNDATIONS.md) — Genomics, variant calling, and annotation fundamentals
- [**Advanced Guide**](HCLS_AI_FACTORY_LEARNING_GUIDE_ADVANCED.md) — RAG architecture, drug discovery, and BioNeMo deep dives

### Reports

- [**Intelligence Report**](HCLS_AI_FACTORY_INTELLIGENCE_REPORT.md) — Pipeline metrics and evidence analysis
- [**Executive Summary**](HCLS_AI_FACTORY_EXECUTIVE_BULLETS.md) — Key talking points

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd hcls-ai-factory

# Configure environment
cp .env.example .env
# Edit .env with your ANTHROPIC_API_KEY and NGC_API_KEY

# Start all services
./start-services.sh
```

Open the landing page at [http://localhost:8080](http://localhost:8080) to access all services.

---

## Technology Stack

| Layer | Technologies |
|---|---|
| Hardware | NVIDIA DGX Spark (GB10 GPU, 128 GB unified memory, 144 ARM64 cores) |
| Genomics | Parabricks 4.6, BWA-MEM2, DeepVariant |
| Vector DB | Milvus 2.4, BGE-small-en-v1.5 (384-dim) |
| LLM | Claude claude-sonnet-4-20250514 (Anthropic) |
| Drug Discovery | BioNeMo MolMIM, DiffDock, RDKit |
| Orchestration | Docker Compose, Nextflow |
| Monitoring | Grafana, Prometheus, DCGM Exporter |
| Frontend | Streamlit, Flask |

---

## License

Apache 2.0 — Free to use, modify, and distribute.

*Built by Adam Jones — February 2026*
