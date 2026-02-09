---
search:
  exclude: true
---

# HCLS AI Factory — Executive Bullets

> **One-page reference for executives, stakeholders, and demo audiences.**
>
> License: Apache 2.0 | Date: February 2026

---

## What It Is

The HCLS AI Factory transforms patient DNA into ranked novel drug candidates in under 5 hours on a single NVIDIA DGX Spark ($3,999). Three GPU-accelerated stages — genomics, RAG-grounded target identification, and AI-driven drug discovery — run end-to-end with no manual intervention.

---

## The Problem

- CPU-based genomics pipelines take **12-36 hours** for a single 30× WGS sample
- Variant annotation is **fragmented** across disconnected databases and manual curation
- The gap from identified variant to drug lead compound is **months of manual work**
- No integrated platform connects genomics, clinical reasoning, and drug discovery

---

## The Solution — Three Stages

### Stage 1: GPU-Accelerated Genomics (120-240 min)
- **NVIDIA Parabricks 4.6** — 10-20× faster than CPU
- BWA-MEM2 alignment: **20-45 min** (vs. 12-24 hours on CPU)
- Google DeepVariant: **10-35 min**, >99% accuracy
- Input: ~200 GB FASTQ (30× WGS, HG002)
- Output: VCF with **~11.7 million variants**

### Stage 2: RAG-Grounded Target Identification (Interactive)
- **3 annotation databases:** ClinVar (4.1M), AlphaMissense (71M), Ensembl VEP
- **3.5 million** high-quality variant embeddings in Milvus vector database
- **Anthropic Claude** (RAG-grounded reasoning) identifies druggable gene targets
- **201 genes** across **13 therapeutic areas**, **171 druggable targets** (85%)
- Output: Target gene with full evidence chain

### Stage 3: AI-Driven Drug Discovery (8-16 min)
- **BioNeMo MolMIM** — generative chemistry (novel molecule design)
- **BioNeMo DiffDock** — molecular docking (binding affinity prediction)
- **RDKit** — drug-likeness scoring (Lipinski, QED, TPSA)
- Composite scoring: 30% generation + 40% docking + 30% QED
- Output: **100 ranked novel drug candidates** + PDF report

---

## Key Numbers

| Metric | Value |
|---|---|
| Total Pipeline Time | < 5 hours |
| Input Data | ~200 GB FASTQ (30× WGS) |
| Variants Called | ~11.7 million |
| High-Quality Variants | ~3.5 million |
| Genes in Knowledge Base | 201 (13 therapeutic areas) |
| Druggable Targets | 171 (85%) |
| Drug Candidates Generated | 100 (ranked by composite score) |
| Hardware Cost | $3,999 (DGX Spark) |

---

## VCP/FTD Demo Highlights

- **Target:** VCP gene — Frontotemporal Dementia, ALS, IBMPFD
- **Variant:** rs188935092 — ClinVar Pathogenic, AlphaMissense 0.87
- **Seed:** CB-5083 (Phase I clinical VCP inhibitor)
- **Result:** Top candidate shows **+39% composite improvement** over seed
  - Docking: -11.4 kcal/mol (vs. -8.1 for CB-5083)
  - QED: 0.81 (vs. 0.62 for CB-5083)
  - All top 10 pass Lipinski's Rule of Five

---

## Technology Stack

| Layer | Technology |
|---|---|
| Hardware | NVIDIA DGX Spark (GB10 GPU, 128 GB unified, $3,999) |
| Genomics | NVIDIA Parabricks 4.6, DeepVariant (>99% accuracy) |
| Annotation | ClinVar, AlphaMissense, Ensembl VEP |
| Vector DB | Milvus 2.4, BGE-small-en-v1.5, IVF_FLAT |
| LLM | Anthropic Claude (RAG-grounded reasoning) |
| Drug Discovery | BioNeMo MolMIM, BioNeMo DiffDock, RDKit |
| Orchestration | Nextflow DSL2 (5 modes: full, target, drug, demo, genomics_only) |
| Monitoring | Grafana, Prometheus, DCGM Exporter |
| License | Apache 2.0 (fully open) |

---

## Deployment Roadmap

| Phase | Hardware | Scale | Cost |
|---|---|---|---|
| 1 — Proof Build | DGX Spark | 1 patient, Docker Compose | $3,999 |
| 2 — Departmental | DGX B200 | Multiple concurrent, Kubernetes | $500K-$1M |
| 3 — Enterprise | DGX SuperPOD | Thousands concurrent, FLARE federated | $7M-$60M+ |

---

## Cross-Modal Integration

- **Imaging → Genomics:** Lung-RADS 4B+ triggers tumor gene profiling
- **Genomics → Drug Discovery:** Pathogenic variants trigger molecule generation
- **Drug Discovery → Imaging:** Combined genomic + imaging clinical reports
- **NVIDIA FLARE:** Federated learning across institutions (data stays local)

---

## Competitive Differentiation

- **Only platform** running genomics-to-drug-candidates on a single desktop GPU
- **End-to-end:** No manual handoffs between stages
- **< 5 hours** total pipeline time (vs. weeks/months traditional)
- **$3,999** proof build cost (vs. $100K+ for equivalent CPU infrastructure)
- **Open project:** Apache 2.0, reproducible, auditable
- **Scalable:** Same Nextflow pipelines scale from DGX Spark to SuperPOD

---

*HCLS AI Factory — Apache 2.0 | February 2026*
