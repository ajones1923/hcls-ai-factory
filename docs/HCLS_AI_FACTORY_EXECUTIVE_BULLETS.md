---
search:
  exclude: true
---

# HCLS AI Factory — Executive Bullets

> **One-page reference for executives, stakeholders, and demo audiences.**
>
> License: Apache 2.0 | Date: March 2026

---

## What It Is

The HCLS AI Factory transforms patient DNA into ranked novel drug candidates in under 5 hours on a single NVIDIA DGX Spark ($4,699). Three GPU-accelerated engines -- Genomic Foundation, Precision Intelligence (11 agents), and Therapeutic Discovery -- run end-to-end with no manual intervention. Eleven domain-specialized intelligence agents provide comprehensive clinical decision support across oncology, cardiology, neurology, rare disease, pharmacogenomics, autoimmune disease, medical imaging, CAR-T therapy, biomarker analysis, single-cell genomics, and clinical trial operations.

---

## The Problem

- CPU-based genomics pipelines take **12-36 hours** for a single 30x WGS sample
- Variant annotation is **fragmented** across disconnected databases and manual curation
- The gap from identified variant to drug lead compound is **months of manual work**
- Clinical decision support is siloed by specialty -- no integrated platform connects genomics, clinical reasoning, and drug discovery
- Access requires $100K+ infrastructure and multiple specialist teams

---

## The Solution — Three Engines, 11 Agents

### Engine 1: Genomic Foundation Engine (120-240 min)
- **NVIDIA Parabricks 4.6** -- 10-20x faster than CPU
- BWA-MEM2 alignment: **20-45 min** (vs. 12-24 hours on CPU)
- Google DeepVariant: **10-35 min**, >99% accuracy
- Input: ~200 GB FASTQ (30x WGS, HG002)
- Output: **~11.7 million variants**, **3.56 million annotated variant** embeddings in Milvus

### Engine 2: Precision Intelligence Network (Interactive)
- **11 intelligence agents** sharing read-only access to 3.56M annotated variant vectors
- **139 Milvus collections** containing **~47,691 agent-owned vectors** across all domains
- **Anthropic Claude** (RAG-grounded reasoning) powers each agent
- **201 genes** across **13 therapeutic areas**, **171 druggable targets** (85%)
- Output: Validated target gene with full evidence chain, clinical reports (PDF, FHIR R4)

**The 11 Intelligence Agents:**

| Agent | Key Capabilities |
|---|---|
| **Precision Oncology** | Molecular tumor board, CIViC/OncoKB annotation, AMP/ASCO/CAP evidence tiers, therapy ranking |
| **Cardiology Intelligence** | 6 risk calculators (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II), GDMT optimizer, 8 workflows |
| **Neurology Intelligence** | 10 clinical scales (NIHSS, GCS, MoCA, MDS-UPDRS, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn-Yahr), 8 workflows |
| **Rare Disease Diagnostic** | 88 rare diseases across 13 categories, 23 ACMG criteria, HPO phenotype matching, GA4GH Phenopacket export |
| **Pharmacogenomics** | 25 pharmacogenes, CPIC/DPWG dosing, phenoconversion detection, HLA hypersensitivity screening |
| **Precision Autoimmune** | 13 autoimmune conditions, autoantibody panels, HLA typing, disease activity scoring, flare prediction |
| **Precision Biomarker** | Biological age estimation (PhenoAge/GrimAge), disease trajectory, pharmacogenomic profiling |
| **CAR-T Intelligence** | Construct comparison (4-1BB vs CD28), manufacturing intelligence, clinical trial matching |
| **Imaging Intelligence** | NVIDIA NIM (VISTA-3D, MAISI, VILA-M3), DICOM ingestion, Lung-RADS, cross-modal genomics triggers |
| **Single-Cell Intelligence** | 57 cell types, TME profiling, spatial niche mapping, drug response prediction, CAR-T target validation |
| **Clinical Trial Intelligence** | Protocol optimization, patient-trial matching, site selection, adaptive design, regulatory documents |

### Engine 3: Therapeutic Discovery Engine (8-16 min)
- **BioNeMo MolMIM** -- generative chemistry (novel molecule design)
- **BioNeMo DiffDock** -- molecular docking (binding affinity prediction)
- **RDKit** -- drug-likeness scoring (Lipinski, QED, TPSA)
- Composite scoring: 30% generation + 40% docking + 30% QED
- Output: **100 ranked novel drug candidates** + PDF report

---

## Key Numbers

| Metric | Value |
|---|---|
| Total Pipeline Time | < 5 hours |
| Input Data | ~200 GB FASTQ (30x WGS) |
| Variants Called | ~11.7 million |
| Annotated Variants | ~3.56 million |
| Intelligence Agents | 11 (spanning 11 medical specialties) |
| Milvus Collections | 139 (agent-owned) + shared genomic evidence |
| Agent Vectors | ~47,691 (domain-specific) |
| Services | 21 (engines + agents + infrastructure) |
| Genes in Knowledge Base | 201 (13 therapeutic areas) |
| Druggable Targets | 171 (85%) |
| Drug Candidates Generated | 100 (ranked by composite score) |
| Test Files | 158 (core + all 11 agents) |
| Hardware Cost | $4,699 (DGX Spark) |

---

## VCP/FTD Demo Highlights

- **Target:** VCP gene -- Frontotemporal Dementia, ALS, IBMPFD
- **Variant:** rs188935092 -- ClinVar Pathogenic, AlphaMissense 0.87
- **Seed:** CB-5083 (Phase I clinical VCP inhibitor)
- **Result:** Top candidate shows **+39% composite improvement** over seed
  - Docking: -11.4 kcal/mol (vs. -8.1 for CB-5083)
  - QED: 0.81 (vs. 0.62 for CB-5083)
  - All top 10 pass Lipinski's Rule of Five

---

## Technology Stack

| Layer | Technology |
|---|---|
| Hardware | NVIDIA DGX Spark (GB10 GPU, 128 GB unified, $4,699) |
| Genomics | NVIDIA Parabricks 4.6, DeepVariant (>99% accuracy) |
| Annotation | ClinVar (4.1M records), AlphaMissense (71M predictions), Ensembl VEP |
| Vector DB | Milvus 2.4, BGE-small-en-v1.5, IVF_FLAT, 139 collections |
| LLM | Anthropic Claude (RAG-grounded reasoning across all 11 agents) |
| Drug Discovery | BioNeMo MolMIM, BioNeMo DiffDock, RDKit |
| Orchestration | Nextflow DSL2 (5 modes: full, target, drug, demo, genomics_only) |
| Monitoring | Grafana, Prometheus, DCGM Exporter |
| License | Apache 2.0 (fully open) |

---

## Deployment Roadmap

| Phase | Hardware | Scale | Cost |
|---|---|---|---|
| 1 -- Proof Build | DGX Spark | 1 patient, Docker Compose, 21 services | $4,699 |
| 2 -- Departmental | DGX B200 | Multiple concurrent, Kubernetes | $500K-$1M |
| 3 -- Enterprise | DGX SuperPOD | Thousands concurrent, FLARE federated | $7M-$60M+ |

---

## Cross-Modal Integration

- **Imaging --> Genomics:** Lung-RADS 4B+ triggers tumor gene profiling (EGFR, ALK, ROS1, KRAS)
- **Genomics --> Drug Discovery:** Pathogenic variants trigger molecule generation
- **Single-Cell --> Oncology:** TME profiling informs immunotherapy selection
- **Pharmacogenomics --> All Agents:** Genotype-guided dosing across clinical domains
- **NVIDIA FLARE:** Federated learning across institutions (data stays local)

---

## Competitive Differentiation

- **Only platform** running genomics-to-drug-candidates with 11 clinical intelligence agents on a single desktop GPU
- **End-to-end:** No manual handoffs between engines
- **< 5 hours** total pipeline time (vs. weeks/months traditional)
- **$4,699** proof build cost (vs. $100K+ for equivalent CPU infrastructure)
- **11 agents** covering oncology, cardiology, neurology, rare disease, pharmacogenomics, autoimmune, biomarker, CAR-T, imaging, single-cell, and clinical trials
- **Open project:** Apache 2.0, reproducible, auditable, 158 test files
- **Scalable:** Same Nextflow pipelines scale from DGX Spark to SuperPOD

---

*HCLS AI Factory -- Apache 2.0 | March 2026*

---

!!! warning "Clinical Decision Support Disclaimer"
    The HCLS AI Factory platform and all intelligence agents described in this document are clinical decision support research tools. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
