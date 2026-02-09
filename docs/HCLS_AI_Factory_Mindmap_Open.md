---
search:
  exclude: true
---

# HCLS AI Factory — Architecture Mindmap (Open)

> **Purpose:** Comprehensive reference map of the HCLS AI Factory platform for Secondary Genomics to Novel Drug Discovery. Covers every component, data flow, technology, and integration point.
>
> License: Apache 2.0 | Date: February 2026

---

## 1. Patient Data Source

- **Input:** Illumina whole-genome sequencing (2×250 bp paired-end reads)
- **Sample:** GIAB HG002 (Genome in a Bottle, NIST reference standard)
- **Coverage:** 30× whole-genome sequencing
- **File Format:** FASTQ (gzipped, ~200 GB per sample)
  - `HG002_R1.fastq.gz` — Forward reads
  - `HG002_R2.fastq.gz` — Reverse reads
- **Reference Genome:** GRCh38 (Human Genome Build 38, 3.1 GB)
- **Sequencing Platform:** Illumina NovaSeq / HiSeq (compatible with any short-read platform)
- **Read Pairs:** ~800 million per 30× genome
- **Base Pairs:** 3 billion in the human genome

---

## 2. HCLS AI Factory Landing Page

- **Service:** Flask web server (port 8080)
- **Function:** Unified entry point for all three pipeline stages
- **Health Monitoring:** Real-time status of 10+ services
  - Port checking via `lsof`
  - 30-second refresh interval
  - Individual service probes
- **Service Dashboard:**
  - Landing Page (8080)
  - Genomics Portal (5000)
  - RAG/Chat API (5001)
  - Chat Interface (8501)
  - Drug Discovery UI (8505)
  - Discovery Portal (8510)
  - Milvus Vector DB (19530)
  - Grafana (3000)
  - Prometheus (9099)
  - Node Exporter (9100)
- **Startup Script:** `start-services.sh` (selective or full startup)
- **Commands:** `--all`, `--landing`, `--rag`, `--drug`, `--status`, `--stop`

---

## 3. Stage 1 — Genomics Pipeline (FASTQ → VCF)

### Compute Platform
- **Technology:** NVIDIA Parabricks 4.6.0-1 (GPU-accelerated bioinformatics)
- **Container:** `nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1`
- **Runtime:** NVIDIA Container Runtime (Docker)

### Pipeline Steps
1. **fq2bam (Alignment):** BWA-MEM2 GPU-accelerated read alignment
   - Maps ~800M read pairs to GRCh38 reference
   - Includes coordinate sorting and duplicate marking
   - Output: BAM file (~100 GB)
   - Time: 20–45 minutes (60% of total)
   - GPU utilization: 70–90% compute, 8–12 GB memory
2. **BAM Indexing:** samtools index
   - Output: BAM index file (~10 MB)
   - Time: 2–5 minutes
3. **DeepVariant (Variant Calling):** Google DeepVariant CNN
   - Accuracy: >99% (validated against GIAB truth set)
   - Output: VCF file (~1.5 GB)
   - Total variants: ~11.7 million
   - Time: 10–35 minutes (35% of total)
   - GPU utilization: 80–95% compute, 12–20 GB memory

### Output Statistics
- Total variants: ~11.7 million
- SNPs: ~4.2 million
- Indels: ~1.0 million
- High-quality (QUAL > 30): ~3.5 million
- In coding regions: ~35,000

### Performance
- **GPU time:** 120–240 minutes end-to-end
- **CPU baseline:** 24–48 hours (10–50× speedup)
- **Web Portal:** Flask (port 5000) — FASTQ upload, alignment monitoring

### Supported GPUs
- DGX Spark (GB10, 128 GB) — Recommended
- A100 (40/80 GB), V100 (32 GB), RTX 4090 (24 GB), RTX 3090 (24 GB)
- Minimum: 8 GB VRAM

---

## 4. Stage 2 — RAG/Chat Pipeline (VCF → Target Hypothesis)

### Variant Annotation Stack
- **ClinVar:** 4.1 million GRCh38 clinical variants (NCBI)
  - Clinical significance classifications (Pathogenic, Likely Pathogenic, VUS, Benign)
  - Disease associations (phenotype lists)
  - dbSNP rsID mapping
- **AlphaMissense:** 71,697,560 AI-predicted missense pathogenicity scores (Google DeepMind)
  - Pathogenic: score > 0.564
  - Ambiguous: 0.340–0.564
  - Benign: score < 0.340
  - High-confidence pathogenic: ≥ 0.80
  - High-confidence benign: ≤ 0.20
- **VEP (Ensembl Variant Effect Predictor):**
  - REST API: `https://rest.ensembl.org/vep/homo_sapiens/hgvs/`
  - Docker: `ensemblorg/ensembl-vep:release_110.1`
  - Functional consequence prediction (SIFT, PolyPhen)

### Annotation Funnel
- 11.7M total variants (from Stage 1)
- 3.5M high-quality variants (QUAL > 30)
- 35,616 ClinVar-annotated variants (0.3%)
- 6,831 AlphaMissense-annotated variants (0.2%)

### Vector Database — Milvus
- **Engine:** Milvus 2.4 (open-source vector database)
- **Port:** 19530 (gRPC/TCP)
- **Web UI (Attu):** Port 8000
- **Embedding Model:** BGE-small-en-v1.5 (BAAI, 384 dimensions)
- **Index Type:** IVF_FLAT (nlist = 1024)
- **Distance Metric:** Cosine similarity
- **Collection:** `genomic_evidence`
- **Indexed Records:** 3.5 million variant embeddings
- **Search Latency:** < 100 ms
- **Schema (17 fields):**
  - id (PK), embedding (384-dim), chrom, pos, ref, alt, qual
  - gene, consequence, impact, genotype, text_summary
  - clinical_significance, rsid, disease_associations
  - am_pathogenicity, am_class

### Knowledge Base (Clinker)
- **201 target genes** across 13 therapeutic areas
- **171 druggable targets** (85% druggability rate)
- **150+ disease conditions**
- **100+ FDA-approved drugs** referenced
- Per-gene data: protein name, function, pathway, diseases, drugs, drug status, PDB IDs, druggability

### 13 Therapeutic Areas
| Area | Gene Count | Key Examples |
|---|---|---|
| Neurology | 36 | VCP, GRN, C9orf72, MAPT, LRRK2, SNCA, HTT, APOE, TREM2 |
| Oncology | 27 | BRCA1, BRCA2, EGFR, KRAS, BRAF, ALK, ERBB2, TP53, PDCD1 |
| Metabolic/Endocrine | 22 | GLP1R, GIPR, GCGR, PPARG, DPP4, SGLT2, PCSK9, HMGCR |
| Infectious Disease | 21 | HIV1_RT, HCV_NS3, SARS2_MPRO, ACE2, RPOB, DHFR |
| Respiratory | 13 | ADRB2, IL5, IL4R, BMPR2, SERPINA1, PDE5A |
| Rare Disease | 12 | CFTR, SMN1, DMD, HBB, F8, GBA, PAH |
| Hematology | 12 | SYK, THPO, JAK2, BTK, F8, ADAMTS13 |
| GI/Hepatology | 12 | NOD2, S1PR1, PNPLA3, FXR, THR_BETA |
| Pharmacogenomics | 11 | CYP2D6, CYP2C19, CYP3A4, VKORC1, DPYD, TPMT |
| Ophthalmology | 11 | VEGFA, CFH, RPE65, RHO, ABCA4 |
| Cardiovascular | 10 | LDLR, PCSK9, TTR, MYBPC3, SCN5A, KCNH2 |
| Immunology | 9 | IL6, TNF, JAK1, JAK2, IL17A, IL23A |
| Dermatology | 9 | IL31RA, TYK2, IL13, COL7A1, FLG |

### LLM Integration
- **Provider:** Anthropic Claude (claude-sonnet-4-20250514)
- **Temperature:** 0.3 (factual consistency)
- **RAG grounding:** All responses cite retrieved variant evidence
- **System prompt:** Expert genomics assistant with variant interpretation, gene function, clinical significance, pharmacogenomics, drug target identification
- **Query expansion:** 10 therapeutic area gene maps (neurodegeneration, pharmacogenomics, metabolic, infectious, respiratory, ophthalmology, pain, GI, oncology, hematology)

### Interfaces
- **Streamlit Chat UI:** Port 8501 (interactive natural language queries)
- **REST API:** Port 5001 (programmatic access)
- **File Manager:** VCF upload, directory browsing, metadata viewing

---

## 5. Stage 3 — Drug Discovery Pipeline (Target → Drug Candidates)

### BioNeMo NIM Services
- **MolMIM (Molecule Generation):** Port 8001
  - Masked modeling for molecule generation
  - Seed compound → novel analogues
  - Temperature-controlled diversity
  - API: POST `/v1/generate`
- **DiffDock (Molecular Docking):** Port 8002
  - Diffusion-based molecular docking
  - Binding pose prediction
  - Affinity scoring
  - API: POST `/v1/dock`

### 10-Stage Pipeline
1. **Initialize:** Load target hypothesis
2. **Normalize Target:** Validate gene, map UniProt ID
3. **Structure Discovery:** Fetch PDB/Cryo-EM structures from RCSB
4. **Structure Prep:** Prepare structures for docking
5. **Molecule Generation:** MolMIM generates candidates from seed SMILES
6. **Chemistry QC:** Lipinski Rule of Five, molecular weight, LogP filters
7. **Conformer Generation:** RDKit 3D conformer + energy minimization
8. **Docking:** DiffDock binding prediction and scoring
9. **Ranking:** Composite score calculation and candidate selection
10. **Reporting:** PDF report generation with ranked candidates

### Scoring System
- **Composite Score:** 30% generation + 40% docking + 30% QED
- **Docking normalization:** `max(0, min(1, (10 + dock_score) / 20))`
- **Lipinski Rule of Five:**
  - MW ≤ 500 Da
  - LogP ≤ 5
  - H-bond donors ≤ 5
  - H-bond acceptors ≤ 10
- **QED (Quantitative Estimate of Drug-likeness):**
  - > 0.67: Drug-like
  - 0.49–0.67: Moderately drug-like
  - < 0.49: Less drug-like
- **Docking affinity:**
  - Good binders: -8 to -12 kcal/mol
  - Minimum binding: -6.0 kcal/mol

### Cryo-EM Structure Evidence
- **Source:** RCSB PDB / EMDB
- **Structure scoring:** Resolution (max 5 Å), inhibitor-bound (+3), druggable pockets (+0.5 each), Cryo-EM method (+0.5)
- **VCP structures:**
  - 8OOI — Wild-type hexamer (Cryo-EM, 2.9 Å)
  - 9DIL — Disease mutant (Cryo-EM, 3.2 Å)
  - 7K56 — Cofactor complex (Cryo-EM, 2.5 Å)
  - 5FTK — VCP + CB-5083 inhibitor (X-ray, 2.3 Å)

### Chemistry Tools
- **RDKit:** Lipinski scoring, QED calculation, TPSA, rotatable bonds, conformer generation
- **Properties computed:** MW, LogP, HBD, HBA, TPSA, rotatable bonds, Lipinski violations, QED, SA score
- **Visualization:** py3Dmol (3D molecular), Plotly (charts)
- **PDF reports:** ReportLab (NVIDIA branding #76B900)

### Interfaces
- **Discovery UI:** Port 8505 (Streamlit — molecule visualization, 3Dmol.js)
- **Discovery Portal:** Port 8510 (pipeline orchestration dashboard)

---

## 6. Nextflow Orchestration

- **Language:** Nextflow DSL2
- **Pipeline Name:** HLS-Pipeline v1.0.0
- **Main file:** `main.nf`
- **Modules:** `genomics.nf`, `rag_chat.nf`, `drug_discovery.nf`, `reporting.nf`
- **Execution modes:**
  - `full` — FASTQ → VCF → Target → Molecules (complete pipeline)
  - `target` — VCF → Target Hypothesis (skip genomics)
  - `drug` — Target → Drug Candidates (discovery only)
  - `demo` — VCP/FTD built-in demonstration
  - `genomics_only` — Stage 1 only
- **Profiles:** standard, docker, singularity, dgx_spark, slurm, test
- **Error handling:** Retry on exit codes 143, 137, 104, 134, 139 (max 2 retries)
- **Reports:** Pipeline HTML, timeline, DAG visualization, trace TSV
- **Container images:**
  - Genomics: `nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1`
  - RAG/Chat: `hls-pipeline/rag-chat:latest`
  - Drug Discovery: `hls-pipeline/drug-discovery:latest`
  - MolMIM: `nvcr.io/nvidia/clara/bionemo-molmim:1.0`
  - DiffDock: `nvcr.io/nvidia/clara/diffdock:1.0`

---

## 7. Monitoring and Infrastructure

### Observability Stack
- **Grafana:** Port 3000 — Real-time dashboards (GPU utilization, pipeline progress, service health)
- **Prometheus:** Port 9099 — Metrics collection and time-series storage
- **Node Exporter:** Port 9100 — System metrics (CPU, RAM, disk I/O)
- **DCGM Exporter:** Port 9400 — GPU metrics (temperature, power, memory, utilization)

### NVIDIA DGX Spark
- **GPU:** GB10 (Blackwell generation)
- **Memory:** 128 GB unified LPDDR5x (shared CPU/GPU)
- **CPU:** ARM64 (Grace) cores
- **Interconnect:** NVLink-C2C (~900 GB/s)
- **Storage:** NVMe SSD
- **CUDA:** 12.x
- **OS:** Ubuntu 22.04 LTS
- **Price:** $3,999

### Minimum Requirements
- GPU: 24 GB VRAM
- RAM: 64 GB
- Storage: 500 GB SSD
- CPU: 8 cores

---

## 8. VCP / Frontotemporal Dementia Demo Workflow

- **Gene:** VCP (Valosin-containing protein, also p97)
- **UniProt:** P55072
- **Function:** AAA+ ATPase for protein quality control (ERAD, autophagy, mitophagy)
- **Diseases:** Frontotemporal Dementia (FTD), ALS, IBMPFD (Inclusion Body Myopathy with Paget's disease and FTD)
- **Variant:** rs188935092 (chr9:35065263 G>A)
- **Total VCP variants found:** 13

### Demo Flow
1. Process HG002 whole-genome (Parabricks → 11.7M variants)
2. Annotate variants (ClinVar + AlphaMissense → identify 13 VCP variants)
3. RAG chat: "What variants are associated with frontotemporal dementia?"
4. Knowledge graph connects VCP variants to FTD via Clinker
5. Retrieve Cryo-EM structures (8OOI, 9DIL, 7K56, 5FTK)
6. Generate 100 novel VCP inhibitor candidates from CB-5083 seed
7. Dock candidates to VCP D2 ATPase domain (ATP-competitive binding pocket)
8. Rank by composite score, generate PDF report

### CB-5083 (Seed Compound)
- Phase I VCP inhibitor
- Binding site: D2 ATPase domain
- Mode: ATP-competitive
- Key residues: ALA464, GLY479, ASP320, GLY215
- Pocket volume: ~450 Å³
- Druggability score: 0.92

---

## 9. Cross-Modal Integration

- **Imaging → Genomics trigger:** Lung-RADS 4B+ finding in Imaging Intelligence Agent → FHIR ServiceRequest → Nextflow pipeline → Parabricks genomics analysis
- **Shared infrastructure:** PostgreSQL (common schema), FHIR ServiceRequest (trigger messages), Event bus (Kafka/Redis Streams at scale)
- **Unified monitoring:** Grafana dashboards spanning imaging and genomics services
- **Combined reports:** Imaging findings + genomic context + drug candidates in single clinical narrative

---

## 10. Technology Stack

| Component | Technology | License |
|---|---|---|
| Genomics alignment | BWA-MEM2 | MIT |
| Variant calling | Google DeepVariant | BSD-3-Clause |
| GPU acceleration | NVIDIA Parabricks | NVIDIA EULA (NVAIE) |
| Vector database | Milvus | Apache 2.0 |
| Embeddings | BGE-small-en-v1.5 (BAAI) | MIT |
| LLM | Anthropic Claude | Commercial API |
| Molecule generation | BioNeMo MolMIM | NVIDIA EULA (NIM) |
| Molecular docking | BioNeMo DiffDock | NVIDIA EULA (NIM) |
| Chemistry toolkit | RDKit | BSD-3-Clause |
| Web interface | Streamlit | Apache 2.0 |
| Pipeline orchestration | Nextflow | Apache 2.0 |
| Monitoring | Grafana + Prometheus | AGPLv3 / Apache 2.0 |
| Containers | Docker | Apache 2.0 |
| PDF generation | ReportLab | BSD |
| VCF parsing | cyvcf2 | MIT |
| 3D visualization | py3Dmol | MIT |

---

## 11. Deployment Roadmap

### Phase 1 — Proof Build (Current)
- Single DGX Spark ($3,999)
- Docker Compose orchestration
- 10+ services on one node
- VCP/FTD end-to-end demonstration
- All HCLS AI Factory code on GitHub (Apache 2.0) — see [Licensing & Cost Guide](licensing.md)

### Phase 2 — Departmental
- DGX B200 (8 GPUs) or small cluster
- Kubernetes orchestration
- Multi-user access with RBAC
- EHR integration (FHIR R4)
- Multiple simultaneous genome analyses
- Horizontal scaling of stateless services

### Phase 3 — Multi-Site / Enterprise
- DGX SuperPOD (256+ GPUs)
- Federated learning across sites (NVIDIA FLARE)
- Multi-site knowledge sharing
- Cross-modal triggers (imaging ↔ genomics ↔ drug discovery)
- Enterprise security (SMART on FHIR, OAuth2)
- Compliance frameworks (HIPAA, GDPR, 21 CFR Part 11)

---

## 12. Data Flow Summary

```
Patient DNA Sample
    ↓
Illumina Sequencer → FASTQ (~200 GB, 30× WGS)
    ↓
[Stage 1: Genomics Pipeline — 120-240 min]
    Parabricks fq2bam (BWA-MEM2) → BAM (~100 GB)
    DeepVariant → VCF (~11.7M variants)
    ↓
[Stage 2: RAG/Chat Pipeline — Interactive]
    ClinVar (4.1M) + AlphaMissense (71M) annotation
    → 3.5M high-quality variants embedded in Milvus
    → 35,616 ClinVar + 6,831 AlphaMissense annotated
    Claude RAG chat → Target hypothesis
    ↓
[Stage 3: Drug Discovery Pipeline — 8-16 min]
    Cryo-EM structure retrieval (RCSB PDB)
    MolMIM molecule generation → 100 candidates
    DiffDock molecular docking → binding predictions
    RDKit scoring (Lipinski, QED)
    → Ranked drug candidates + PDF report
```

**End-to-End:** < 5 hours (vs. 2–3 years traditional)

---

> *This document was created for the HCLS AI Factory — Secondary Genomics to Novel Drug Discovery.*
> *Apache 2.0 License | February 2026*
