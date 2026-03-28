---
tags:
  - Performance
  - Benchmarks
  - DGX Spark
  - Metrics
---

# Performance Benchmarks

Measured performance of the HCLS AI Factory on NVIDIA DGX Spark ($4,699). All timings represent end-to-end wall clock time under default configurations.

---

## Summary

| Metric | Value |
|---|---|
| **End-to-end time** (DNA to Drug Candidates) | < 5 hours |
| **Traditional approach** | 6-18 months |
| **Reduction** | ~99% |
| **Hardware cost** | $4,699 (single workstation) |
| **Traditional infrastructure cost** | $50K-500K+ (cluster + licenses) |
| **Intelligence agents** | 11 |
| **Milvus collections** | 139 |
| **Total vectors** | ~47,691 |
| **Services** | 21 |

---

## Hardware: NVIDIA DGX Spark

| Component | Specification |
|---|---|
| GPU | NVIDIA GB10 Grace Blackwell |
| Memory | 128 GB unified LPDDR5x |
| CPU | 20 ARM cores (Grace) |
| Interconnect | NVLink-C2C (GPU-CPU) |
| Storage | NVMe SSD |
| Power | Desktop form factor |
| Price | $4,699 |

---

## Engine 1: GPU Genomics (FASTQ to VCF)

**Pipeline:** NVIDIA Parabricks 4.6 with BWA-MEM2 + DeepVariant

| Step | Time | GPU Utilization |
|---|---|---|
| Alignment (BWA-MEM2) | 20-45 min | 85-95% |
| Sorting + deduplication | included | -- |
| Indexing (samtools) | 2-5 min | CPU |
| Variant calling (DeepVariant) | 10-35 min | 85-95% |
| **Total** | **120-240 min** | **85-95%** |

| Metric | Value |
|---|---|
| Input | ~200 GB paired-end FASTQ (HG002 WGS) |
| Output | ~11.7 million variant calls (VCF) |
| Accuracy | >99% concordance (DeepVariant) |
| Speedup vs. CPU | 10-50x |
| CPU baseline | 24-48 hours |

---

## Engine 2: Evidence RAG (VCF to Target Hypothesis)

**Pipeline:** Milvus + BGE-small-en-v1.5 + Claude

### Vector Database

| Collection | Records | Embedding Time |
|---|---|---|
| ClinVar variants | ~2.7M records | ~45 min (one-time) |
| AlphaMissense predictions | 71M records (sampled) | ~30 min (one-time) |
| Clinker knowledge base | 201 genes, 150+ diseases | ~5 min (one-time) |
| **Total annotated variants** | **3.56M** | -- |

### Query Performance

| Operation | Latency |
|---|---|
| Vector embedding (BGE-small-en-v1.5) | < 50 ms |
| Milvus similarity search (top-10) | < 100 ms |
| Claude evidence synthesis | 2-5 sec |
| **End-to-end query** | **< 5 sec** |

| Metric | Value |
|---|---|
| Embedding model | BGE-small-en-v1.5 (384 dimensions) |
| LLM | Claude (Anthropic) |
| Therapeutic areas covered | 13 |
| Target genes | 201 |
| Druggability rate | 85% |

---

## Engine 3: Drug Discovery (Target to Molecules)

**Pipeline:** BioNeMo MolMIM + DiffDock + RDKit

| Step | Time | Mode |
|---|---|---|
| Structure retrieval (RCSB PDB) | < 5 sec | API |
| Structure preparation | < 30 sec | CPU |
| Molecule generation (MolMIM) | 10-60 sec | Cloud NIM |
| 3D conformer generation (RDKit) | < 30 sec | CPU |
| Molecular docking (DiffDock) | 2-8 min | Cloud NIM |
| Scoring and ranking (QED + Lipinski) | < 10 sec | CPU |
| Report generation | < 30 sec | CPU |
| **Total** | **8-16 min** | -- |

| Metric | Value |
|---|---|
| Candidate molecules generated | 10-100 per run |
| Docking poses per molecule | 10 |
| Drug-likeness filter | Lipinski Rule of 5 + QED |
| Seed compound (demo) | CB-5083 (VCP inhibitor) |
| PDB structures used (demo) | 5FTK, 8OOI, 9DIL, 7K56 |

### NIM Execution Modes

| Mode | MolMIM | DiffDock | Best For |
|---|---|---|---|
| Cloud | health.api.nvidia.com | health.api.nvidia.com | DGX Spark (ARM64), no local GPU containers needed |
| Local | localhost:8001 | localhost:8002 | x86 workstations with dedicated GPU |
| Mock | Simulated output | Simulated output | Testing, CI/CD, demos without API keys |

---

## Intelligence Agent Benchmarks

### 1. Precision Biomarker Agent (Port 8502)

| Metric | Value |
|---|---|
| Collections | 11 (10 owned + 1 shared) |
| Biomarker interpretation | < 3 sec |
| Biological age estimation (PhenoAge/GrimAge) | < 2 sec |
| Pharmacogenomic profiling | < 3 sec |
| Disease trajectory detection | < 5 sec |
| PDF + FHIR R4 export | < 5 sec |

### 2. Precision Oncology Agent (Port 8503)

| Metric | Value |
|---|---|
| Collections | 11 (10 owned + 1 shared) |
| Case creation | < 2 sec |
| MTB packet generation | 10-30 sec |
| Trial matching | < 5 sec |
| Therapy ranking (CIViC/OncoKB) | < 5 sec |
| FHIR R4 bundle export | < 2 sec |
| Test suite | 516 tests, < 1 sec |

### 3. CAR-T Intelligence Agent (Port 8504)

| Metric | Value |
|---|---|
| Collections | 11 (10 owned + 1 shared) |
| Total vectors | 6,266+ |
| Query latency (evidence retrieval) | < 3 sec |
| Comparative analysis (e.g., 4-1BB vs CD28) | < 8 sec |
| Deep research mode | 10-30 sec |
| PDF export | < 5 sec |
| Test suite | 241 tests, < 1 sec |

### 4. Imaging Intelligence Agent (Port 8505)

| Metric | Value |
|---|---|
| Collections | 10 |
| NIM services | 4 (VISTA-3D, MAISI, VILA-M3, Llama-3) |
| Workflow demo execution | 5-15 sec per modality |
| Evidence query | < 5 sec |
| Comparative analysis | < 8 sec |
| FHIR R4 DiagnosticReport export | < 2 sec |
| Test suite | 539 tests, ~3 sec |

### 5. Precision Autoimmune Agent (Port 8506)

| Metric | Value |
|---|---|
| Collections | 10 |
| Autoimmune conditions covered | 13 |
| Diagnostic engine evaluation | < 3 sec |
| Disease timeline construction | < 5 sec |
| Cross-modal genomic enrichment | < 3 sec |
| PDF clinical report export | < 5 sec |
| Evidence query | < 5 sec |

### 6. Pharmacogenomics (PGx) Agent (Port 8507)

| Metric | Value |
|---|---|
| Collections | 15 |
| Pharmacogenes | 25 |
| Drugs covered | 100+ |
| CPIC dosing algorithms | 9 |
| HLA associations | 15 |
| Star allele to phenotype conversion | < 1 sec |
| Phenoconversion (inhibitor adjustment) | < 1 sec |
| HLA adverse reaction screening | < 1 sec |
| Full PGx pipeline (genotype to dosing) | < 3 sec |
| FHIR R4 PGx report export | < 2 sec |
| Test suite | 1,001+ tests |

### 7. Cardiology Intelligence Agent (Port 8527)

| Metric | Value |
|---|---|
| Collections | 13 (12 owned + 1 shared) |
| Risk calculators | 6 (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II) |
| ASCVD 10-year risk calculation | 14.6% example in < 1 sec |
| HEART score calculation | < 1 sec |
| CHA2DS2-VASc stroke risk | < 1 sec |
| GDMT optimization (4-pillar HFrEF) | < 3 sec |
| Clinical workflows | 8 (CAD, HF, valvular, arrhythmia, cardiac MRI, stress test, prevention, cardio-oncology) |
| Cross-modal genomic triggers | < 3 sec (cardiomyopathies, channelopathies, FH) |
| Evidence query | < 5 sec |
| FHIR R4 export | < 2 sec |
| Test suite | 1,927 tests, all passing |

### 8. Neurology Intelligence Agent (Port 8528)

| Metric | Value |
|---|---|
| Collections | 14 (13 owned + 1 shared) |
| Clinical scales | 10 (NIHSS, GCS, MoCA, MDS-UPDRS, EDSS, mRS, HIT-6, ALSFRS-R, ASPECTS, Hoehn & Yahr) |
| NIHSS stroke severity calculation | < 1 sec |
| GCS assessment | < 1 sec |
| MoCA cognitive screening | < 1 sec |
| MDS-UPDRS motor examination | < 1 sec |
| EDSS disability scoring (MS) | < 1 sec |
| ASPECTS early CT scoring | < 1 sec |
| Sub-domains covered | 8 (cerebrovascular, degenerative, epilepsy, neuro-oncology, MS, movement, headache, neuromuscular) |
| Clinical workflows | 8 (acute stroke, dementia eval, epilepsy focus, brain tumor, MS monitoring, Parkinson's, headache, neuromuscular) |
| Evidence query | < 5 sec |
| FHIR R4 export | < 2 sec |

### 9. Rare Disease Diagnostic Agent (Port 8526)

| Metric | Value |
|---|---|
| Collections | 14 (13 owned + 1 shared) |
| Rare diseases covered | 88 |
| ACMG variant classification criteria | 23 |
| HPO phenotype matching | < 2 sec |
| Phenotype-to-gene ranking | < 3 sec |
| WES/WGS variant interpretation | < 5 sec |
| Gene therapy eligibility assessment | < 3 sec |
| Newborn screening evaluation | < 2 sec |
| Clinical workflows | 10 (phenotype-driven, WES/WGS, metabolic screening, dysmorphology, neurogenetic, cardiac genetics, connective tissue, inborn errors, gene therapy, undiagnosed disease) |
| Evidence query | < 5 sec |
| FHIR R4 export | < 2 sec |

### 10. Clinical Trial Intelligence Agent (Port 8521)

| Metric | Value |
|---|---|
| Collections | 14 (13 owned + 1 shared) |
| Patient-trial matching | < 5 sec |
| Eligibility criteria analysis | < 3 sec |
| Protocol design assistance | 10-30 sec |
| Site selection ranking | < 5 sec |
| Safety signal detection | < 5 sec |
| Adaptive design evaluation | < 5 sec |
| Regulatory document intelligence | < 5 sec |
| Competitive intelligence query | < 5 sec |
| Clinical workflows | 11 (protocol design, patient matching, site selection, eligibility optimization, adaptive design, safety signal, regulatory docs, competitive intel, diversity assessment, decentralized planning, general) |
| FHIR R4 / DOCX export | < 5 sec |

### 11. Single-Cell Intelligence Agent (Port 8525)

| Metric | Value |
|---|---|
| Collections | 12 (11 owned + 1 shared) |
| Cell types annotated | 57 |
| Cell type identification query | < 3 sec |
| TME profiling (immune phenotype) | < 5 sec |
| Spatial niche analysis | < 5 sec |
| Drug response prediction | < 5 sec |
| Trajectory analysis | < 5 sec |
| Ligand-receptor interaction query | < 5 sec |
| CAR-T target validation | < 5 sec |
| Biomarker discovery query | < 5 sec |
| Clinical workflows | 10 (cell type annotation, TME profiling, drug response, subclonal architecture, spatial niche, trajectory analysis, ligand-receptor, biomarker discovery, CAR-T target, treatment monitoring) |
| Evidence query | < 5 sec |

### Combined Test Suite

| Agent | Tests | Time |
|---|---|---|
| Precision Biomarker | -- | -- |
| Precision Oncology | 516 | 0.40 sec |
| CAR-T Intelligence | 241 | 0.18 sec |
| Imaging Intelligence | 539 | 3.20 sec |
| Precision Autoimmune | -- | -- |
| Pharmacogenomics | 1,001+ | -- |
| Cardiology Intelligence | 1,927 | -- |
| Neurology Intelligence | -- | -- |
| Rare Disease Diagnostic | -- | -- |
| Clinical Trial Intelligence | -- | -- |
| Single-Cell Intelligence | -- | -- |
| **Total (measured)** | **4,224+** | -- |

---

## Infrastructure

### Service Startup

| Component | Cold Start | Warm Restart |
|---|---|---|
| Milvus | 30-60 sec | 10-15 sec |
| Landing page | < 5 sec | < 2 sec |
| RAG Chat UI | 5-10 sec | < 5 sec |
| Drug Discovery UI | 5-10 sec | < 5 sec |
| Agent UIs (11 agents) | 5-10 sec each | < 5 sec |
| Full platform (21 services) | 2-3 min | < 1 min |

### Resource Usage (Idle)

| Resource | Usage |
|---|---|
| CPU | < 5% (20 ARM cores) |
| Memory | ~8-12 GB (128 GB available) |
| GPU Memory | < 2 GB (128 GB unified) |
| Disk (platform + data) | ~400-500 GB |

### Resource Usage (Peak -- Genomics Pipeline Running)

| Resource | Usage |
|---|---|
| CPU | 40-60% |
| Memory | 30-50 GB |
| GPU | 85-95% utilization |
| Disk I/O | High (FASTQ read + BAM write) |

---

## Scalability

| Dimension | Current | Potential |
|---|---|---|
| Samples per day | 3-6 (sequential) | 10-20 (with pipeline parallelism) |
| Vector database | ~47,691 vectors across 139 collections | 100M+ (Milvus scales horizontally) |
| Knowledge base | 201 genes, 13 areas | Expandable with additional collections |
| Concurrent users | 5-10 (single workstation) | 50+ (with load balancing) |
| Intelligence agents | 11 | Additional agents via plugin architecture |

---

## Methodology

- All benchmarks measured on NVIDIA DGX Spark ($4,699) with Ubuntu 22.04 LTS
- Timings are wall clock measurements averaged over 3 runs
- GPU utilization measured via `nvidia-smi` and DCGM Exporter
- Query latencies measured end-to-end including network overhead
- Risk calculator and clinical scale timings measured as pure computation (no network)
- Test suite timings measured via `pytest` with default configuration
- "Traditional approach" estimates based on published literature for manual genomics + drug discovery workflows at academic medical centers
- Agent benchmarks represent typical single-query performance under normal load

---

!!! warning "Clinical Decision Support Disclaimer"
    The HCLS AI Factory platform and its components are clinical decision support research tools. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
