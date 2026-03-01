---
tags:
  - Performance
  - Benchmarks
  - DGX Spark
  - Metrics
---

# Performance Benchmarks

Measured performance of the HCLS AI Factory on NVIDIA DGX Spark ($3,999). All timings represent end-to-end wall clock time under default configurations.

---

## Summary

| Metric | Value |
|---|---|
| **End-to-end time** (DNA to Drug Candidates) | < 5 hours |
| **Traditional approach** | 6-18 months |
| **Reduction** | ~99% |
| **Hardware cost** | $3,999 (single workstation) |
| **Traditional infrastructure cost** | $50K-500K+ (cluster + licenses) |

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

---

## Stage 1: GPU Genomics (FASTQ to VCF)

**Pipeline:** NVIDIA Parabricks 4.6 with BWA-MEM2 + DeepVariant

| Step | Time | GPU Utilization |
|---|---|---|
| Alignment (BWA-MEM2) | 20-45 min | 85-95% |
| Sorting + deduplication | included | — |
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

## Stage 2: Evidence RAG (VCF to Target Hypothesis)

**Pipeline:** Milvus + BGE-small-en-v1.5 + Claude

### Vector Database

| Collection | Records | Embedding Time |
|---|---|---|
| ClinVar variants | ~2.7M records | ~45 min (one-time) |
| AlphaMissense predictions | 71M records (sampled) | ~30 min (one-time) |
| Clinker knowledge base | 201 genes, 150+ diseases | ~5 min (one-time) |
| **Total searchable vectors** | **3.56M** | — |

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

## Stage 3: Drug Discovery (Target to Molecules)

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
| **Total** | **8-16 min** | — |

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

### CAR-T Intelligence Agent (Port 8521)

| Metric | Value |
|---|---|
| Collections | 10 owned + 1 shared (read-only) |
| Total vectors | 6,266+ |
| Query latency (evidence retrieval) | < 3 sec |
| Comparative analysis | < 8 sec |
| Deep research mode | 10-30 sec |
| PDF export | < 5 sec |
| Test suite | 241 tests, < 1 sec |

### Imaging Intelligence Agent (Port 8525)

| Metric | Value |
|---|---|
| Collections | 10 |
| NIM services | 4 (VISTA-3D, MAISI, VILA-M3, Llama-3) |
| Workflow demo execution | 5-15 sec per modality |
| Evidence query | < 5 sec |
| Comparative analysis | < 8 sec |
| FHIR R4 DiagnosticReport export | < 2 sec |
| Test suite | 539 tests, ~3 sec |

### Precision Oncology Agent (Port 8526)

| Metric | Value |
|---|---|
| Collections | 11 (10 owned + 1 shared) |
| Case creation | < 2 sec |
| MTB packet generation | 10-30 sec |
| Trial matching | < 5 sec |
| Therapy ranking | < 5 sec |
| FHIR R4 bundle export | < 2 sec |
| Test suite | 516 tests, < 1 sec |

### Combined Test Suite

| Agent | Tests | Time |
|---|---|---|
| CAR-T Intelligence | 241 | 0.18 sec |
| Imaging Intelligence | 539 | 3.20 sec |
| Precision Oncology | 516 | 0.40 sec |
| **Total** | **1,296** | **3.78 sec** |

---

## Infrastructure

### Service Startup

| Component | Cold Start | Warm Restart |
|---|---|---|
| Milvus | 30-60 sec | 10-15 sec |
| Landing page | < 5 sec | < 2 sec |
| RAG Chat UI | 5-10 sec | < 5 sec |
| Drug Discovery UI | 5-10 sec | < 5 sec |
| Agent UIs | 5-10 sec each | < 5 sec |
| Full platform (all services) | 2-3 min | < 1 min |

### Resource Usage (Idle)

| Resource | Usage |
|---|---|
| CPU | < 5% (20 ARM cores) |
| Memory | ~8-12 GB (128 GB available) |
| GPU Memory | < 2 GB (128 GB unified) |
| Disk (platform + data) | ~400-500 GB |

### Resource Usage (Peak — Genomics Pipeline Running)

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
| Vector database | 3.56M vectors | 100M+ (Milvus scales horizontally) |
| Knowledge base | 201 genes, 13 areas | Expandable with additional collections |
| Concurrent users | 5-10 (single workstation) | 50+ (with load balancing) |
| Agent instances | 3 | Additional agents via plugin architecture |

---

## Methodology

- All benchmarks measured on NVIDIA DGX Spark with Ubuntu 22.04 LTS
- Timings are wall clock measurements averaged over 3 runs
- GPU utilization measured via `nvidia-smi` and DCGM Exporter
- Query latencies measured end-to-end including network overhead
- Test suite timings measured via `pytest` with default configuration
- "Traditional approach" estimates based on published literature for manual genomics + drug discovery workflows at academic medical centers
