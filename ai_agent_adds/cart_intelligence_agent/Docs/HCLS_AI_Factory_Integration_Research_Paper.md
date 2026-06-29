# HCLS AI Factory: An Integrated Precision Medicine Platform on VAST AI OS for Accelerated Genomics, Variant Interpretation, CAR-T Therapy Intelligence, and AI-Driven Drug Discovery

**Adam Jones**
Healthcare & Life Sciences AI Factory Project

**Date:** February 27, 2026

**Version:** 2.0

---

## Abstract

The Healthcare and Life Sciences (HCLS) AI Factory is an end-to-end precision medicine platform that transforms raw patient DNA sequencing data into actionable therapeutic candidates in under five hours on a single NVIDIA DGX Spark workstation ($4,699). The platform integrates four interconnected computational stages: (1) GPU-accelerated genomic variant calling via NVIDIA Parabricks 4.6, producing 11.7 million variants from whole-genome sequencing data in 2-4 hours; (2) Retrieval-Augmented Generation (RAG) for variant interpretation, leveraging 3.56 million vectors across clinical annotations (ClinVar, AlphaMissense) with a 201-gene knowledge graph spanning 13 therapeutic domains; (3) a CAR-T Intelligence Agent providing autonomous multi-collection reasoning across 11 specialized vector collections encompassing 6,266+ curated records on chimeric antigen receptor T-cell therapies; and (4) AI-driven drug discovery using BioNeMo NIM microservices (MolMIM and DiffDock) for de novo molecule generation, protein-ligand docking, and composite ranking in 8-16 minutes. This paper presents the complete architecture, details the technical implementation of each stage, identifies integration opportunities between the genomics pipeline, RAG engine, CAR-T agent, and drug discovery system, and proposes a phased research plan for transforming these currently semi-independent components into a unified, bidirectionally communicating precision medicine platform. The system is Apache 2.0 licensed and designed for deployment on commodity NVIDIA hardware, democratizing access to computational precision medicine workflows that traditionally require institutional-scale infrastructure.

An AI Factory of this scope must ingest extensive amounts of sequencing data, trigger GPU pipelines, embed millions of records into vector space, run RAG across clinical evidence, orchestrate autonomous drug design agents, and deliver it all at clinical-grade latency. That is not a storage problem. That is an operating system problem. This paper positions VAST AI OS---with its DataStore, DataBase, DataEngine, InsightEngine, AgentEngine, and ICMS subsystems---as the unified infrastructure foundation that collapses the middleware stack (Milvus, Docker Compose, Nextflow, custom RAG code, in-memory Python dictionaries) into a single, coherent data operating system. By migrating the HCLS AI Factory to VAST AI OS, we eliminate 17 GB of per-patient inter-service data movement, replace six independent infrastructure components with native platform services, and unlock event-driven pipeline orchestration, petabyte-scale genomic storage with zero-copy GPU access, unified SQL+vector search, and autonomous agent coordination---all on the same data plane.

**Keywords:** precision medicine, CAR-T cell therapy, retrieval-augmented generation, drug discovery, NVIDIA Parabricks, vector database, large language models, genomic variant interpretation, BioNeMo, DGX Spark, VAST AI OS, data operating system, DataEngine, InsightEngine, AgentEngine

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [System Architecture Overview](#2-system-architecture-overview)
3. [Hardware Platform: NVIDIA DGX Spark](#3-hardware-platform-nvidia-dgx-spark)
4. [Stage 1: GPU-Accelerated Genomics Pipeline](#4-stage-1-gpu-accelerated-genomics-pipeline)
5. [Stage 2: RAG-Powered Variant Interpretation](#5-stage-2-rag-powered-variant-interpretation)
6. [Stage 3: CAR-T Intelligence Agent](#6-stage-3-cart-intelligence-agent)
7. [Stage 4: AI-Driven Drug Discovery](#7-stage-4-ai-driven-drug-discovery)
8. [Orchestration Layer: Nextflow DSL2](#8-orchestration-layer-nextflow-dsl2)
9. [Shared Infrastructure and Services](#9-shared-infrastructure-and-services)
10. [Current Integration Analysis](#10-current-integration-analysis)
11. [Proposed Integration Research Plan](#11-proposed-integration-research-plan)
12. [Implementation Roadmap](#12-implementation-roadmap)
13. [Performance Benchmarks and Metrics](#13-performance-benchmarks-and-metrics)
14. [Security, Observability, and Production Hardening](#14-security-observability-and-production-hardening)
15. [Discussion](#15-discussion)
16. [Future Directions](#16-future-directions)
17. [Conclusion](#17-conclusion)
18. [References](#18-references)
19. [Appendices](#19-appendices)

*Note: Sections 1.5, 2.4, and 11.1 (Phase 0) introduce the VAST AI OS infrastructure foundation. Sections 9, 10, 15, 16, and 17 have been updated in Version 2.0 to reflect the VAST AI OS migration thesis.*

---

## 1. Introduction

### 1.1 The Precision Medicine Challenge

Precision medicine promises to tailor therapeutic interventions to individual patients based on their genomic profiles. However, the computational journey from raw sequencing data to actionable treatment recommendations remains fragmented across specialized tools, each requiring distinct expertise and infrastructure. A typical precision medicine workflow involves: (a) secondary analysis of whole-genome sequencing (WGS) data to identify variants, (b) clinical annotation and interpretation of those variants against databases such as ClinVar and AlphaMissense, (c) identification of therapeutic targets based on variant-driven evidence, and (d) computational exploration of candidate molecules or cell therapies directed at those targets. Each of these stages traditionally operates in isolation, with manual handoffs between bioinformaticians, clinical geneticists, immunologists, and medicinal chemists.

### 1.2 The Cell Therapy Imperative

Chimeric Antigen Receptor T-cell (CAR-T) therapies represent one of the most significant advances in cancer treatment. Six FDA-approved CAR-T products---Kymriah (tisagenlecleucel), Yescarta (axicabtagene ciloleucel), Tecartus (brexucabtagene autoleucel), Breyanzi (lisocabtagene maraleucel), Abecma (idecabtagene vicleucel), and Carvykti (ciltacabtagene autoleucel)---have demonstrated remarkable efficacy in hematological malignancies. Yet the field faces challenges in target antigen selection, resistance mechanism prediction, toxicity management (cytokine release syndrome, immune effector cell-associated neurotoxicity syndrome), and manufacturing optimization. These challenges require synthesis of evidence across clinical trials, molecular biology, regulatory submissions, real-world outcomes, and patient genomics---a cognitive task well-suited to AI augmentation.

### 1.3 Contribution of This Work

This paper presents the HCLS AI Factory, an integrated platform that unifies genomic analysis, evidence-based variant interpretation, CAR-T therapy intelligence, and computational drug discovery into a single, end-to-end system. Our specific contributions include:

1. **Architecture and implementation** of a four-stage pipeline executing on a single DGX Spark workstation, reducing both cost and complexity compared to distributed cloud solutions.
2. **A CAR-T Intelligence Agent** with 11 specialized Milvus vector collections, 25 target antigen profiles, 12 query expansion maps (8,000+ terms), and autonomous reasoning capabilities for cross-functional therapy evaluation.
3. **Integration analysis and research plan** identifying concrete opportunities for bidirectional communication between genomics, variant interpretation, CAR-T therapy selection, and drug discovery stages.
4. **Open-source release** under Apache 2.0, enabling academic and clinical research organizations to deploy the complete system on commodity hardware.

### 1.4 Demonstration Target

Throughout this paper, we use the VCP gene (Valosin-Containing Protein, p97) as our demonstration target. VCP is an AAA+ ATPase involved in protein quality control whose mutations cause Frontotemporal Dementia with Inclusion Body Myopathy and Paget Disease (IBMPFD). The seed compound CB-5083 (a Phase I VCP inhibitor) and PDB structures 5FTK, 8OOI, 9DIL, and 7K56 serve as reference inputs for the drug discovery stage.

### 1.5 The Infrastructure Operating System

The HCLS AI Factory, as described above, requires the simultaneous coordination of petabyte-scale genomic data ingestion, GPU-accelerated pipeline execution, embedding of millions of records into vector space, retrieval-augmented generation across clinical evidence, autonomous orchestration of drug design agents, and delivery of results at clinical-grade latency. Each of these demands has historically been addressed by a separate infrastructure component: Lustre or GPFS for parallel file access, Docker Compose and Nextflow for container orchestration, Milvus for vector search, PostgreSQL for relational queries, custom Python code for RAG pipelines, and bespoke agent frameworks for multi-step reasoning. The result is a middleware stack where data must be copied, serialized, and deserialized across process boundaries---incurring what we term the **egress/ingress tax**.

For a single patient genome processed through the full pipeline, we measure approximately 17 GB of inter-service data movement: 200 GB FASTQ ingested to local NVMe, 120 GB BAM written and re-read, 3 GB VCF parsed and re-embedded, 35,616 variant embeddings serialized to Milvus over TCP, RAG context assembled in Python memory and sent to an external LLM, and drug discovery results passed between MolMIM, DiffDock, and RDKit through HTTP APIs. This is not a storage problem. This is an operating system problem.

VAST AI OS addresses this by providing a unified data operating system with purpose-built subsystems that map directly to the HCLS AI Factory's requirements:

- **DataStore**: Petabyte-scale genomic data storage with zero-copy access, replacing Lustre/GPFS/S3 with a single namespace that GPUs and CPUs can access without serialization overhead.
- **DataSpace**: Secure multi-tenant data isolation across research teams and clinical workflows, replacing manual namespace management and ad hoc access controls.
- **DataEngine**: Event-driven pipeline triggers and serverless container execution on CNodes, replacing the Nextflow/Docker Compose orchestration layer with native data-aware scheduling.
- **DataBase**: Unified SQL and vector search in a single query engine---eliminating the need for separate Milvus, PostgreSQL, and Python in-memory dictionaries.
- **InsightEngine**: Knowledge extraction and RAG without manual ETL, replacing 3,000+ lines of custom embedding, annotation, and retrieval code.
- **AgentEngine**: Autonomous orchestration for drug discovery and CAR-T intelligence agents, replacing custom agent frameworks with a native agent execution environment.
- **ICMS (BlueField-4)**: 10-20x faster LLM time-to-first-token via petabyte-scale KV cache on DPU memory, replacing GPU HBM-only caching strategies.
- **SyncEngine**: Multi-site replication for multi-center clinical deployments, enabling federated genomics without custom synchronization code.
- **DASE Architecture**: Disaggregated Shared-Everything---stateless CNodes for compute elasticity, stateful DBoxes for data durability---enabling the AI Factory to scale from a single DGX Spark to a multi-rack deployment without architectural changes.

The remainder of this paper details the current HCLS AI Factory architecture and then maps each component to its VAST AI OS equivalent, presenting a migration path that eliminates the middleware stack while preserving all pipeline capabilities.

---

## 2. System Architecture Overview

### 2.1 Pipeline Topology

The HCLS AI Factory implements a linear four-stage pipeline with lateral connections through a shared vector database (Milvus 2.4+) and embedding model (BGE-small-en-v1.5, 384 dimensions):

```
FASTQ Reads (200GB)
    │
    ▼
┌──────────────────────────────────────────┐
│  STAGE 1: Genomics Pipeline (Port 5000)  │
│  Parabricks 4.6 · BWA-MEM2 · DeepVariant│
│  Input: FASTQ → Output: VCF (11.7M var)  │
│  Duration: 120-240 minutes               │
└──────────────────────────────────────────┘
    │ VCF
    ▼
┌──────────────────────────────────────────┐
│  STAGE 2: RAG Chat Pipeline (Port 8501)  │
│  Milvus · ClinVar · AlphaMissense · KG  │
│  35,616 annotated variants · 201 genes   │
│  Duration: Interactive                   │
└──────────────────────────────────────────┘
    │ Target Hypothesis JSON
    ▼
┌──────────────────────────────────────────┐  ┌───────────────────────────────┐
│  STAGE 3: Drug Discovery (Port 8505)     │  │  CAR-T Intelligence Agent     │
│  MolMIM · DiffDock · RDKit · Ranking     │  │  (Port 8521)                  │
│  10-stage pipeline · 20 molecules        │  │  11 collections · 6,266+ vec  │
│  Duration: 8-16 minutes                  │  │  25 antigens · Agent reasoning│
└──────────────────────────────────────────┘  └───────────────────────────────┘
    │                                              │
    ▼                                              ▼
┌──────────────────────────────────────────────────────┐
│              Unified Clinical Report                  │
│  Genomic Profile + Targets + Therapies + Molecules   │
└──────────────────────────────────────────────────────┘
```

### 2.2 Service Architecture

The platform deploys 15 services across distinct ports, monitored by a centralized landing page dashboard:

| Service | Port | Technology | Purpose |
|---------|------|------------|---------|
| Landing Page | 8080 | Flask | Service health dashboard |
| Genomics Portal | 5000 | Flask | Variant calling control |
| RAG API | 5001 | Flask/FastAPI | Programmatic variant queries |
| RAG Chat UI | 8501 | Streamlit | Interactive variant analysis |
| Drug Discovery UI | 8505 | Streamlit | Molecule generation |
| Discovery Portal | 8510 | Streamlit | Unified drug discovery |
| CAR-T Agent UI | 8521 | Streamlit | CAR-T therapy intelligence |
| Milvus | 19530 | Milvus 2.4 | Vector database |
| Attu | 8000 | Web UI | Milvus administration |
| MolMIM NIM | 8001 | BioNeMo | Molecule generation |
| DiffDock NIM | 8002 | BioNeMo | Molecular docking |
| Grafana | 3000 | Grafana | Metrics visualization |
| Prometheus | 9099 | Prometheus | Metrics collection |
| Node Exporter | 9100 | Prometheus | System metrics |
| DCGM Exporter | 9400 | NVIDIA | GPU metrics |

### 2.3 Shared Infrastructure

All stages share three critical components:

1. **Milvus 2.4+ Vector Database** (port 19530): Stores all embedding vectors using IVF_FLAT indexing with COSINE similarity metric, 1,024 partitions (nlist), and 16 probes per search (nprobe).

2. **BGE-small-en-v1.5 Embedding Model**: 384-dimensional embeddings from the BAAI/bge-small-en-v1.5 sentence-transformer model, used identically across the RAG pipeline and CAR-T Agent.

3. **Claude API** (Anthropic): Both the RAG pipeline and CAR-T Agent use Claude Sonnet/Opus for response synthesis, with domain-specific system prompts.

### 2.4 VAST AI OS Foundation

Under the VAST AI OS migration, the shared infrastructure described in Section 2.3 consolidates into native platform services. The following table maps each current component to its VAST AI OS equivalent:

| Current Component | VAST AI OS Service | Migration Impact |
|---|---|---|
| Milvus 2.4+ (port 19530) | **DataBase** native VECTOR columns | Eliminates standalone vector DB, TCP serialization, and IVF_FLAT index tuning. DataBase provides unified SQL+vector queries in a single engine. |
| BGE-small-en-v1.5 (custom Python) | **InsightEngine** built-in embedding | InsightEngine manages embedding generation, chunking, and indexing natively. Custom `embedder.py` code is replaced. |
| Claude API (custom RAG code) | **InsightEngine** RAG pipeline | InsightEngine orchestrates retrieval, context assembly, and LLM synthesis without 3,000+ lines of custom RAG code. |
| Docker Compose (15 services) | **DataEngine** serverless CNodes | DataEngine executes pipeline stages as event-driven functions on stateless compute nodes (CNodes), eliminating port management and Docker networking. |
| Nextflow DSL2 orchestrator | **DataEngine** event triggers | DataEngine's native event system (S3 object creation, data mutation) replaces Nextflow DAG scheduling with data-aware triggers. |
| Local NVMe (500GB) | **DataStore** petabyte namespace | DataStore provides a single, globally addressable namespace with zero-copy GPU access via GPUDirect, eliminating local storage bottlenecks. |
| Python in-memory dicts (ClinVar 2.7M, AlphaMissense 71M) | **DataBase** SQL tables | ClinVar and AlphaMissense data loads into DataBase tables, eliminating 12 GB of Python process RAM and enabling SQL joins across clinical annotations. |
| Manual namespace management | **DataSpace** multi-tenant isolation | DataSpace provides cryptographic data isolation per patient, research team, or clinical workflow---without filesystem-level ACLs. |
| Custom agent frameworks | **AgentEngine** autonomous orchestration | AgentEngine provides a native execution environment for the CAR-T Intelligence Agent and drug discovery pipeline, with built-in tool use, memory, and coordination. |
| GPU HBM-only KV cache | **ICMS (BlueField-4)** petabyte KV cache | ICMS extends LLM KV cache to DPU memory, achieving 10-20x faster time-to-first-token for RAG synthesis queries. |

This mapping preserves all existing pipeline capabilities while collapsing the middleware stack into a unified data operating system. The DASE (Disaggregated Shared-Everything) architecture ensures that compute (stateless CNodes) and data (stateful DBoxes) scale independently, enabling the HCLS AI Factory to grow from a single DGX Spark to a multi-rack clinical deployment without re-architecting.

---

## 3. Hardware Platform: NVIDIA DGX Spark

### 3.1 Specifications

The entire platform is designed to execute on a single NVIDIA DGX Spark workstation:

| Component | Specification |
|-----------|---------------|
| GPU | NVIDIA GB10 |
| CPU | 20 ARM cores (NVIDIA Grace) |
| Memory | 128GB unified LPDDR5x |
| Interconnect | NVLink-C2C (CPU-GPU) |
| Storage | NVMe SSD (500GB minimum) |
| Price | $4,699 |
| Power | Desktop form factor |

### 3.2 Design Rationale

The DGX Spark provides sufficient compute for each pipeline stage: Parabricks 4.6 requires a minimum of 8GB GPU memory for variant calling, the BGE-small-en-v1.5 model fits within 1GB, Milvus operates efficiently with 8-16GB RAM for million-scale vector collections, and the BioNeMo NIM services (MolMIM, DiffDock) run within the GPU memory budget. The 128GB unified memory architecture eliminates PCIe bottlenecks through NVLink-C2C, enabling seamless data movement between CPU and GPU memory spaces.

### 3.3 Resource Allocation Strategy

The Nextflow orchestrator configures resource allocation per stage:

| Process Label | CPUs | Memory | GPU | Use Case |
|---------------|------|--------|-----|----------|
| `process_low` | 2 | 4 GB | 0 | Lightweight tasks |
| `process_medium` | 8 | 32 GB | 0 | Standard processing |
| `process_high` | 16 | 64 GB | 0 | Memory-intensive |
| `gpu` | 8 | 32 GB | 1 | GPU-accelerated workloads |

Error recovery uses exit code-based retry with codes [143, 137, 104, 134, 139] triggering up to 2 retries.

---

## 4. Stage 1: GPU-Accelerated Genomics Pipeline

### 4.1 Overview

The genomics pipeline transforms raw whole-genome sequencing data (FASTQ format) into variant calls (VCF format) using NVIDIA Parabricks 4.6.0-1 with GPU acceleration. The pipeline processes the GIAB HG002 reference sample (Ashkenazi individual, Illumina 2x250bp paired-end, 30x depth) through four computational steps.

### 4.2 Processing Steps

**Step 1: Read Alignment (fq2bam)**

The `pbrun fq2bam` command executes GPU-accelerated BWA-MEM2 alignment:

```
pbrun fq2bam \
    --ref GRCh38.fa \
    --in-fq HG002_R1.fastq.gz HG002_R2.fastq.gz \
    --out-bam HG002.genome.bam \
    --num-gpus 1
```

BWA-MEM2 aligns ~800 million read pairs against the GRCh38 reference genome (3.1GB), producing a coordinate-sorted, deduplicated BAM file of approximately 100-120GB. GPU utilization peaks at 80-95% during seeding and extension, consuming 12-16GB of GPU memory. Duration: 60-180 minutes for full genome.

**Step 2: BAM Indexing and Quality Control**

Samtools generates a BAM index (`.bai`) for random access and produces alignment quality metrics via `flagstat`. Expected metrics for HG002: 98.7% mapping rate, 96.8% proper pairing, ~10M duplicate reads. Duration: 2-5 minutes.

**Step 3: Variant Calling (DeepVariant)**

The `pbrun deepvariant` command uses a convolutional neural network (CNN) trained on millions of validated variants to call SNPs and indels:

```
pbrun deepvariant \
    --ref GRCh38.fa \
    --in-bam HG002.genome.bam \
    --out-variants HG002.genome.vcf.gz \
    --num-gpus 1
```

DeepVariant achieves >99% concordance with GIAB truth sets. Output: ~11.7 million variants (4-5M SNPs, 600-800K indels, 100K+ multi-allelic sites) in a gzip-compressed, tabix-indexed VCF file of 2-3GB. Duration: 60-90 minutes.

**Step 4: VCF Indexing**

Tabix creates a block-based index enabling region-specific queries (e.g., `tabix file.vcf.gz chr1:100000-200000`). Duration: 1-2 minutes.

### 4.3 Web Portal

A Flask-based web portal (port 5000, 1,157 lines) provides:

- **Real-time monitoring** via Server-Sent Events (SSE) at 5 Hz update frequency
- **GPU metrics** via pynvml: utilization, memory, power draw, temperature
- **Pipeline control** endpoints: run, stop, reset with API key authentication
- **Rate limiting**: 60 requests per 60-second window per IP
- **Security**: shell metacharacter filtering, path traversal prevention, CSRF tokens, security headers (X-Frame-Options, CSP, X-XSS-Protection)
- **Error resilience**: 6-retry verification for BAM output, 3-retry DeepVariant with 30-second delays, and a custom nvidia-smi wrapper to handle DGX Spark GPU memory reporting quirks

### 4.4 Output Specification

The VCF output conforms to VCF 4.2 format with genotype fields (GT:GQ:DP:AD), enabling downstream annotation with ClinVar clinical significance, AlphaMissense pathogenicity scores, and Ensembl VEP consequence predictions.

---

## 5. Stage 2: RAG-Powered Variant Interpretation

### 5.1 Architecture

The RAG Chat Pipeline implements a Retrieval-Augmented Generation system for genomic variant interpretation. The pipeline ingests the VCF output from Stage 1, annotates variants with clinical and functional data, generates 384-dimensional embeddings, stores them in Milvus, and provides an interactive Streamlit chat interface for evidence-based variant analysis.

### 5.2 Annotation Pipeline

Variants undergo three annotation layers:

**ClinVar Annotation**: The ClinVarAnnotator loads the variant_summary.txt.gz database (~2.7 million GRCh38 variants) into memory, creating lookup keys of the form `chrom_pos_ref_alt`. Each matched variant receives: clinical significance (Pathogenic, Benign, VUS, etc.), phenotype associations, dbSNP rsID, HGVS nomenclature, and review status.

**AlphaMissense Annotation**: The AlphaMissenseAnnotator loads Google DeepMind's 71 million missense variant predictions (requiring 8-10GB RAM), assigning pathogenicity scores (0-1) and classifications: likely benign (<0.34), ambiguous (0.34-0.564), or likely pathogenic (>0.564). Variants scoring >=0.8 receive a `high_pathogenic` tag.

**VEP Annotation**: The LocalVEPAnnotator runs Ensembl VEP (release 110.1) via Docker in offline mode, adding gene symbol, biotype, consequence, impact level (HIGH/MODERATE/LOW/MODIFIER), SIFT, and PolyPhen predictions.

### 5.3 Vector Database Schema

The `genomic_evidence` collection stores 35,616 annotated variants with the following schema:

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(200) | chr_pos_ref_alt (primary key) |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| chrom | VARCHAR(10) | Chromosome |
| pos | INT64 | Genomic position |
| ref / alt | VARCHAR(500) | Reference/alternate alleles |
| qual | FLOAT | Quality score |
| gene | VARCHAR(50) | Gene symbol (HGNC) |
| consequence | VARCHAR(100) | VEP consequence |
| impact | VARCHAR(20) | Impact level |
| genotype | VARCHAR(10) | Sample genotype |
| text_summary | VARCHAR(2000) | Human-readable summary |
| clinical_significance | VARCHAR(200) | ClinVar annotation |
| rsid | VARCHAR(20) | dbSNP identifier |
| am_pathogenicity | FLOAT | AlphaMissense score |
| am_class | VARCHAR(30) | AlphaMissense classification |

Each variant's `text_summary` is a composite narrative: *"Variant at chr9:35065000 A>G in gene VCP (missense_variant) with HIGH impact. HGVS: p.Arg159His. dbSNP: rs188935092. Clinical: Pathogenic. AlphaMissense: 0.912 (likely_pathogenic)."*

### 5.4 Knowledge Graph

The system includes a Clinker-style knowledge graph mapping 201 genes across 13 therapeutic domains:

| Domain | Genes | Examples |
|--------|-------|---------|
| Neurodegeneration | 36 | VCP, C9orf72, GRN, MAPT, TBK1, FUS, TARDBP, SOD1 |
| Oncology | 27 | BRCA1, EGFR, ALK, BRAF, KRAS, ERBB2, PIK3CA |
| Metabolic | 22 | GLP1R, GIPR, INSR, PPARG, PCSK9, HMGCR |
| Infectious Disease | 21 | HIV1_RT, HCV_NS3, SARS2_MPRO, ACE2 |
| Pharmacogenomics | 11 | CYP2D6, CYP2C19, CYP3A4, VKORC1, DPYD |
| Cardiovascular | 10 | PCSK9, TTR, MYBPC3, SCN5A |
| Respiratory | 13 | ADRB2, IL4R, IL5, BMPR2 |
| Ophthalmology | 11 | VEGFA, CFH, RPE65 |
| Hematology | 12 | SYK, THPO, F10, ADAMTS13 |
| GI/Hepatology | 12 | NOD2, S1PR1, PNPLA3 |
| Immunology | 9 | IL6, TNF, JAK1, IL17A |
| Dermatology | 9 | IL31RA, TYK2, KIT |
| Rare Disease | 12 | Various monogenic targets |

Each gene entry maps: protein name, molecular function, biological pathway, disease associations, known drugs (with development status), PDB structure identifiers, and druggability classification. Of 201 genes, 171 (85%) are classified as druggable.

### 5.5 RAG Engine

The RAGEngine class implements the complete retrieval-generation pipeline:

1. **Query Embedding**: Text is prepended with the BGE instruction prefix (*"Represent this sentence for searching relevant passages:"*) before encoding to a 384-dimensional vector.

2. **Semantic Search**: The query vector is searched against the genomic_evidence collection with a default score threshold of 0.5, returning the top 10 results ranked by cosine similarity.

3. **Therapeutic Domain Expansion**: Query keywords are matched against 13 therapeutic area mappings, triggering gene-specific retrieval with a synthetic confidence score of 0.85. For example, the query *"What variants affect FTD treatment?"* expands to genes [VCP, C9orf72, GRN, MAPT, TBK1, FUS, TARDBP, SOD1].

4. **Knowledge Graph Augmentation**: Matched genes are enriched with pathway context, drug candidates, PDB structures, and disease associations from the KNOWLEDGE_CONNECTIONS dictionary.

5. **LLM Synthesis**: Evidence and knowledge context are assembled into a prompt for Claude, which generates a comprehensive response with variant-level citations.

### 5.6 Target Hypothesis Management

The system tracks therapeutic target candidates through a structured lifecycle:

```python
class TargetHypothesis:
    gene: str                    # VCP
    protein: str                 # p97/VCP ATPase
    rationale: str               # Scientific justification
    confidence: str              # low / medium / high
    priority: int                # 1-5
    status: str                  # hypothesis → validated → selected → rejected
    pdb_ids: List[str]           # ["5FTK", "8OOI"]
    reference_smiles: str        # CB-5083 SMILES
    druggability: str            # high / medium / low
```

Validated targets are exported as JSON for handoff to Stage 3 (Drug Discovery), preserving the complete chain of genomic evidence.

### 5.7 Multi-Provider LLM Support

The platform supports four LLM providers through a polymorphic client architecture:

| Provider | Models | Deployment |
|----------|--------|------------|
| Anthropic | Claude Sonnet 4, Claude Opus 4 | Cloud API |
| OpenAI | GPT-4 Turbo, GPT-4 | Cloud API |
| Ollama | Llama 3.1 8B, Llama 3.1 70B | Local (ARM-optimized) |
| vLLM | Llama-3.1-70B-Instruct | Local (x86) |

This flexibility enables deployment across cloud-connected and air-gapped environments.

---

## 6. Stage 3: CAR-T Intelligence Agent

### 6.1 Overview

The CAR-T Intelligence Agent is an autonomous reasoning system for CAR-T cell therapy evaluation. Unlike the RAG pipeline's single-collection approach, the CAR-T Agent searches across 11 specialized Milvus collections simultaneously, applies domain-specific query expansion, and uses an agent reasoning loop for multi-step evidence gathering.

### 6.2 Collection Architecture

The agent manages 11 purpose-built collections, each with a domain-specific schema:

| Collection | Records | Weight | Key Fields |
|------------|---------|--------|------------|
| cart_literature | ~1,200 | 1.2 | PMID, journal, cart_stage, keywords |
| cart_trials | ~800 | 1.3 | NCT number, phase, sponsor, enrollment |
| cart_constructs | ~400 | 1.1 | scFv origin, costimulatory domain, vector type |
| cart_assays | ~600 | 1.0 | assay type, cell line, E:T ratio, key metric |
| cart_manufacturing | ~500 | 0.9 | process step, parameter, target spec |
| cart_safety | ~700 | 1.0 | event type, severity, onset timing, incidence |
| cart_biomarkers | ~400 | 1.1 | biomarker name, assay method, cutoff |
| cart_regulatory | ~300 | 0.8 | product, regulatory event, agency, decision |
| cart_sequences | ~350 | 1.0 | scFv clone, binding affinity, species origin |
| cart_realworld | ~500 | 0.9 | study type, data source, population size |
| genomic_evidence | ~35,616 | 0.7 | (Read-only, shared with RAG pipeline) |

All collections use 384-dimensional BGE-small-en-v1.5 embeddings with IVF_FLAT indexing (COSINE metric, nlist=1024, nprobe=16). Collection weights influence ranked relevance scores: `weighted_score = raw_score × (1 + weight)`.

### 6.3 25 Target Antigen Profiles

The knowledge graph encodes detailed profiles for 25 CAR-T target antigens:

**CD19**: B-Lymphocyte Antigen (UniProt P15391), expressed on B-cell lineage (pro-B to mature B, excluding plasma cells). Approved products: Kymriah, Yescarta, Tecartus, Breyanzi. Key trials: ELIANA, ZUMA-1, ZUMA-2, TRANSFORM, TRANSCEND. Known resistance: CD19 loss/mutation, lineage switch, trogocytosis, alternative splicing (exon 2 deletion). Toxicity: CRS 30-90%, ICANS 20-65%, B-cell aplasia (expected). Normal tissue: B-cells (acceptable on-target/off-tumor).

Additional antigens include BCMA, CD22, CD20, CD30, CD33, CD38, CD123, GD2, HER2, GPC3, EGFR, EGFRvIII, Mesothelin, Claudin18.2, ROR1, PSMA, IL13RA2, CD5, CD7, FLT3, DLL3, B7-H3, and MUC1---each with equivalent depth of characterization.

### 6.4 Query Expansion System

Twelve expansion dictionaries with 1,200+ keywords and 8,000+ total expansion terms enable semantic coverage beyond literal query matching:

| Category | Keywords | Example |
|----------|----------|---------|
| Target Antigen | 25 | "CD19" → [B-ALL, DLBCL, Kymriah, FMC63, ...] |
| Disease | 15 | "DLBCL" → [r/r DLBCL, GCB subtype, ABC subtype, ...] |
| Toxicity | 15 | "CRS" → [tocilizumab, IL-6, ferritin, CRP, ...] |
| Manufacturing | 20 | "Transduction" → [MOI, VCN, RetroNectin, ...] |
| Mechanism | 15 | "Exhaustion" → [PD-1, LAG-3, TOX, ...] |
| Construct | 20 | "Bispecific" → [CD19/CD22, tandem, OR-gate, ...] |
| Safety | 8 | "Cardiac" → [cardiomyopathy, troponin, ...] |
| Biomarker | 10 | "Ferritin" → [serum ferritin, CRS prediction, ...] |
| Regulatory | 8 | "BLA" → [biologics license, NDA, ...] |
| Sequence | 8 | "CDR" → [hypervariable region, CDR3, paratope, ...] |
| RealWorld | 10 | "Registry" → [CIBMTR, EBMT, DESCAR-T, ...] |
| Immunogenicity | 11 | "HLA" → [MHC, epitope, ...] |

The expansion strategy operates in two modes: antigen-detected terms become Milvus field filters on target_antigen-capable collections, while non-antigen terms are re-embedded for semantic search with a 0.7-0.8 score discount factor. A maximum of 5 expansion terms are processed per query.

### 6.5 Agent Reasoning Loop

The CARTIntelligenceAgent implements a multi-step reasoning pipeline:

1. **Search Planning** (`search_plan`): Analyzes the query to identify target antigens (from 25 known), map to CAR-T development stages (TARGET_ID, CAR_DESIGN, VECTOR_ENG, TESTING, CLINICAL), and determine search strategy (broad, targeted, or comparative).

2. **Evidence Retrieval** (`retrieve`): Executes parallel multi-collection search via ThreadPoolExecutor, applies query expansion, merges and deduplicates results (capped at 30), and augments with full knowledge graph context.

3. **Evidence Evaluation**: Classifies retrieved evidence quality:
   - **Sufficient**: >=3 collections represented + >=10 total hits
   - **Partial**: >=2 collections + >=5 hits
   - **Insufficient**: below thresholds

4. **Adaptive Expansion**: When evidence is insufficient, the agent decomposes the query into sub-questions, retrieves evidence for each, and extends the hit pool.

5. **Comparative Analysis**: Queries containing "vs" or "versus" trigger dual-entity retrieval, producing separate evidence panels with a structured comparison prompt.

### 6.6 Citation and Relevance System

Evidence citations are classified into three tiers based on raw cosine similarity scores:

| Tier | Score Threshold | Display |
|------|----------------|---------|
| HIGH | >= 0.7 | Green badge, prioritized |
| MEDIUM | >= 0.5 | Orange badge |
| LOW | < 0.5 | Red badge |

Citations include clickable links to PubMed (PMID-based) and ClinicalTrials.gov (NCT number-based) for primary source verification.

### 6.7 Streamlit User Interface

The CAR-T Agent UI (port 8521) provides three main tabs:

**Chat Interface**: Supports both Quick RAG and Deep Research modes. Sidebar controls include target antigen filter (15 options + "All Targets"), development stage filter, date range sliders (2010-2030), collection selection with live record counts, and 13 pre-written demo queries.

**Knowledge Graph Explorer**: Interactive visualization of target, toxicity, manufacturing, biomarker, and regulatory knowledge using pyvis (with text-based fallback). Includes cross-collection entity search.

**Image Analysis**: Upload CAR-T-related figures (PNG, JPG, PDF) for Claude Vision extraction, claims parsing, and evidence verification against the 11 collections.

### 6.8 Export Engine

The agent produces three export formats:

- **Markdown**: Structured report with collection-specific evidence tables
- **JSON**: Pydantic-serialized data with full metadata
- **PDF**: Branded report using reportlab Platypus with NVIDIA-themed dark headers (#1B1B2F), green accent strips (#76B900), collection-specific accent colors (11 distinct colors), proper markdown-to-flowable conversion, clickable citation links, and page numbers with branded footers

---

## 7. Stage 4: AI-Driven Drug Discovery

### 7.1 Overview

The drug discovery pipeline implements a 10-stage workflow that transforms a validated therapeutic target into ranked drug candidates in 8-16 minutes. The pipeline uses NVIDIA BioNeMo NIM microservices for molecule generation (MolMIM) and protein-ligand docking (DiffDock), with RDKit for chemistry quality control and composite scoring.

### 7.2 Ten-Stage Pipeline

| Stage | Operation | Duration |
|-------|-----------|----------|
| 0 | Initialize: Load config, create target, check NIM health | <1 min |
| 1 | Normalize Target: Validate gene, enrich hypothesis | <1 min |
| 2 | Structure Discovery: Load PDB IDs, fetch from RCSB, select best | <1 min |
| 3 | Structure Prep: Remove water/ions, add hydrogens, identify binding site | <1 min |
| 4 | Molecule Generation: MolMIM de novo generation from seed SMILES | 2-5 min |
| 5 | Chemistry QC: Lipinski RoF5, QED scoring, filter drug-like | <1 min |
| 6 | Conformer Generation: SMILES → 3D via RDKit/MMFF force field | <1 min |
| 7 | Docking: DiffDock protein-ligand docking, 10 poses per molecule | 5-10 min |
| 8 | Ranking: Composite scoring (40% docking + 30% generation + 30% QED) | <1 min |
| 9 | Reporting: Serialize JSON, generate PDF report | <1 min |

### 7.3 MolMIM Integration

MolMIM uses masked language modeling on a learned molecular distribution to generate novel molecules from a seed compound:

**API Endpoint**: `POST http://localhost:8001/v1/generate`

**Parameters**:
- `smiles`: Seed molecule (e.g., CB-5083 for VCP)
- `num_molecules`: 10-50 candidates to generate
- `temperature`: 1.0 (sampling temperature)
- `num_samples_per_token`: 10 (MolMIM sampling strategy)
- `masked_ratio`: 0.1 (fraction of SMILES tokens to mask)

The cloud variant uses NVIDIA's hosted API at `health.api.nvidia.com` with CMA-ES (Covariance Matrix Adaptation Evolution Strategy) for property optimization (QED maximization).

### 7.4 Chemistry Quality Control

Generated molecules are filtered through Lipinski's Rule of Five:

| Property | Threshold | Calculation |
|----------|-----------|-------------|
| Molecular Weight | <= 550 Da | RDKit Descriptors.MolWt |
| LogP | <= 5.0 | Wildman-Crippen (MolLogP) |
| H-bond Donors | <= 5 | Lipinski.NumHDonors |
| H-bond Acceptors | <= 10 | Lipinski.NumHAcceptors |
| Max Violations | <= 1 | Sum of exceeded thresholds |

QED (Quantitative Estimate of Drug-likeness) provides a continuous 0-1 score incorporating MW, logP, HBD, HBA, rotatable bonds, and aromatic ring count.

### 7.5 DiffDock Integration

DiffDock uses a diffusion generative model for molecular docking:

**API Endpoint**: `POST http://localhost:8002/v1/dock`

**Parameters**:
- `protein`: PDB file content or path
- `ligand`: SMILES string
- `num_poses`: 10 conformations per molecule

**Output**: Pose-level data including docking score (kcal/mol, lower is better: -12 = excellent, 0 = no binding), confidence (0-1), hydrogen bond count, and contact residue list.

The cloud variant stages PDB and SDF files as NVCF assets before submission, with parameters `time_divisions=20` and `steps=18`.

### 7.6 Composite Ranking

Candidates are ranked by a weighted composite score:

```
dock_normalized = max(0.0, min(1.0, -dock_score / 12.0))

composite = (0.4 × dock_normalized) + (0.3 × generation_score) + (0.3 × qed_score)
```

The top N candidates (default: 10) are selected for the final report.

### 7.7 Service Mode Selection

The pipeline auto-selects NIM service mode:

1. **Cloud**: If `NIM_MODE=cloud` or API key has `nvapi-` prefix → NVIDIA hosted endpoints
2. **Local**: If Docker NIM containers detected → local endpoints (ports 8001/8002)
3. **Mock**: If `NIM_ALLOW_MOCK_FALLBACK=true` → deterministic mock services with VCP inhibitor analogues
4. **Error**: If no services available → RuntimeError

### 7.8 Checkpoint and Resume

A 10-stage checkpoint system serializes complete pipeline state (config, target, structures, molecules, docking results, rankings) to JSON after each stage, enabling resume-after-failure:

```
outputs/checkpoints/checkpoint_{run_id}_stage{stage}.json
```

### 7.9 VCP Demo Configuration

The demonstration target uses:

- **Gene**: VCP (p97, AAA+ ATPase)
- **Seed Compound**: CB-5083 (Phase I inhibitor)
- **SMILES**: `CC(C)C1=C(C=C(C=C1)NC2=NC3=C(C=N2)N(C=C3)C)C(=O)NC4=CC=C(C=C4)CN5CCOCC5`
- **PDB Structures**: 5FTK (2.3A, CB-5083 complex), 8OOI (2.9A, WT hexamer), 9DIL (3.2A, mutant), 7K56 (2.5A, complex)
- **Binding Site**: D2 ATPase domain (ATP-competitive)
- **Pre-designed Analogues**: 9 VCP inhibitor variants with systematic modifications (morpholine→piperazine, isopropyl→methyl, methyl→ethyl, etc.)

---

## 8. Orchestration Layer: Nextflow DSL2

### 8.1 Pipeline Configuration

The Nextflow DSL2 orchestrator (v1.0.3, `HLS-Pipeline`) coordinates the three primary stages through a modular workflow system:

```groovy
include { GENOMICS_PIPELINE } from './modules/genomics'
include { RAG_CHAT_PIPELINE } from './modules/rag_chat'
include { DRUG_DISCOVERY    } from './modules/drug_discovery'
include { GENERATE_REPORT   } from './modules/reporting'
```

### 8.2 Execution Modes

Five workflow modes support different entry points:

| Mode | Input | Pipeline | Output |
|------|-------|----------|--------|
| `full` | FASTQ file pairs | Genomics → RAG → Drug Discovery → Report | Complete analysis |
| `target` | Pre-computed VCF | RAG → Drug Discovery → Report | Skip genomics |
| `drug` | Target JSON | Drug Discovery → Report | Skip genomics + RAG |
| `demo` | Built-in VCP data | Drug Discovery → Report | Quick demonstration |
| `genomics_only` | FASTQ file pairs | Genomics only | VCF output |

### 8.3 Configuration Parameters

Key parameters from `nextflow.config`:

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `genome` | GRCh38 | Reference genome build |
| `rag_model` | claude-sonnet-4-20250514 | LLM for variant interpretation |
| `max_targets` | 5 | Maximum target hypotheses |
| `confidence_threshold` | 0.7 | Minimum variant call confidence |
| `num_molecules` | 20 | Candidate molecules to generate |
| `diversity` | 0.3 | Molecular diversity constraint |
| `max_mw` | 550 | Maximum molecular weight (Da) |
| `docking_poses` | 10 | Conformations per molecule |
| `max_memory` | 128.GB | DGX Spark memory ceiling |
| `max_cpus` | 32 | Grace CPU core count |
| `max_gpus` | 1 | GB10 GPU count |

### 8.4 Execution Profiles

Six profiles optimize for different deployment environments:

- **standard**: Local executor
- **docker**: Docker with user ID mapping
- **singularity**: Singularity with automounts (HPC)
- **dgx_spark**: Docker with `--gpus all`, GPU/memory optimized
- **slurm**: SLURM queue-based execution (HPC clusters)
- **test**: Demo mode with reduced molecule count

### 8.5 Container Assignments

| Process Pattern | Container Image |
|----------------|-----------------|
| `GENOMICS_.*` | nfcore/sarek:3.4.0 |
| `RAG_CHAT_.*` | hls-pipeline/rag-chat:latest |
| `DRUG_DISCOVERY_.*` | hls-pipeline/drug-discovery:latest |
| `MOLMIM_.*` | nvcr.io/nim/nvidia/molmim:1.0.0 |
| `DIFFDOCK_.*` | nvcr.io/nim/mit/diffdock:2.2.0 |

### 8.6 Reporting and Tracing

The orchestrator generates four reports per run:
- Pipeline execution report (HTML)
- Timeline visualization (HTML)
- Trace data (TSV: task_id, hash, status, duration, CPU, RSS, VMEM)
- DAG visualization (HTML)

---

## 9. Shared Infrastructure and Services

### 9.1 Milvus Vector Database (Current Architecture)

Milvus 2.4+ serves as the unified vector store for all RAG operations:

| Configuration | Value |
|---------------|-------|
| Port | 19530 |
| Index Type | IVF_FLAT |
| Metric | COSINE |
| nlist (partitions) | 1,024 |
| nprobe (search probes) | 16 |
| Embedding Dimension | 384 |

Total vector inventory:
- genomic_evidence: 35,616 variants
- CAR-T collections (10): 6,266+ records
- **Total: ~41,882 vectors**

### 9.2 Embedding Model

BGE-small-en-v1.5 (BAAI) generates all embeddings:
- Dimension: 384
- Normalization: L2-normalized
- Query prefix: "Represent this sentence for searching relevant passages: "
- Batch size: 32
- Speed: ~50-100 embeddings/second on GPU

### 9.3 Landing Page Dashboard

The Flask landing page (port 8080) provides centralized health monitoring for 11 services:

- **Health checking**: HTTP GET or TCP socket probing per service
- **Parallel execution**: ThreadPoolExecutor (10 workers) for simultaneous checks
- **Auto-start capability**: Genomics Portal (port 5000) and RAG API (port 5001) are auto-launched if their ports are closed
- **Readiness gate**: `/api/ready` returns HTTP 503 until critical services (genomics, rag-portal, drug-main) respond
- **Report freshness**: `/api/report-status` tracks PDF report generation timestamps

### 9.4 Monitoring Stack

- **Grafana** (port 3000): Dashboard visualization
- **Prometheus** (port 9099): Time-series metrics collection
- **Node Exporter** (port 9100): System-level metrics (CPU, memory, disk, network)
- **DCGM Exporter** (port 9400): GPU-specific metrics (utilization, temperature, power, memory bandwidth)

### 9.5 VAST AI OS Infrastructure Replacement

The infrastructure described in Sections 9.1--9.4 represents six independent services (Milvus, BGE embedder, Flask dashboard, Grafana, Prometheus, Docker Compose) that must be deployed, configured, monitored, and maintained separately. On VAST AI OS, this collapses into native platform capabilities:

**9.5.1 DataBase Replaces Milvus + Python Dictionaries**

VAST DataBase provides unified SQL and vector search in a single query engine. The 12 Milvus collections (1 genomic_evidence + 10 CAR-T + 1 clinvar) migrate to DataBase tables with native `VECTOR(384)` column types. Key advantages:

- **No standalone vector database**: Eliminates Milvus deployment (port 19530), Attu admin UI (port 8000), and the `pymilvus` client dependency.
- **SQL + vector in one query**: A single DataBase query can join ClinVar clinical_significance (SQL WHERE clause) with embedding similarity (vector distance function), replacing the two-step pattern of Milvus search followed by Python post-filtering.
- **Cross-collection queries**: The current architecture cannot query across CAR-T and genomic evidence collections in a single operation. DataBase enables `SELECT * FROM genomic_evidence ge JOIN cart_trials ct ON ge.gene = ct.target_antigen WHERE vector_distance(ge.embedding, query_vec) < 0.5`---a query that currently requires custom Python orchestration across two separate Milvus searches.
- **Eliminate 12 GB RAM overhead**: The ClinVar annotator loads 2.7 million variants into a Python dictionary (~4 GB). The AlphaMissense annotator loads 71 million predictions (~8 GB). On DataBase, these become indexed SQL tables queried on demand, freeing 12 GB of process memory.

**9.5.2 DataStore Replaces Local NVMe**

VAST DataStore provides a single, petabyte-scale namespace with zero-copy access for both CPUs and GPUs:

- **Zero-copy genomic data path**: FASTQ reads (200 GB), BAM alignments (120 GB), and VCF outputs (3 GB) reside in DataStore. Parabricks accesses them via GPUDirect without staging to local NVMe, eliminating the largest single data copy in the pipeline.
- **Single namespace**: All pipeline outputs---VCFs, annotated variants, target hypotheses, molecule candidates, docking results---live in one addressable namespace. No more scattered files across `/data/genomics/`, `/data/rag/`, `/data/drug_discovery/` with manual path management.
- **Immutable audit trail**: DataStore versioning provides a complete audit trail of every pipeline run, every VCF annotation, and every drug candidate generated---critical for clinical-grade reproducibility.

**9.5.3 DataEngine Replaces Docker Compose + Nextflow**

VAST DataEngine provides event-driven pipeline orchestration with serverless container execution on CNodes:

- **Event-driven triggers**: When a FASTQ file lands in DataStore, DataEngine automatically triggers the Parabricks pipeline. When the VCF is written, annotation begins. When annotations complete, RAG ingestion starts. This replaces Nextflow's DAG scheduling with native data-aware event triggers.
- **Serverless CNode execution**: Pipeline stages execute as serverless functions on stateless CNodes, eliminating Docker Compose service management (15 services, 15 ports, 15 health checks) and Nextflow process configuration (CPU, memory, GPU allocation per stage).
- **No port management**: The current architecture assigns 15 distinct ports and relies on TCP health checks. DataEngine routes requests through the data plane, eliminating port conflicts and service discovery overhead.

**9.5.4 InsightEngine Replaces Custom RAG Code**

VAST InsightEngine provides knowledge extraction and RAG without manual ETL:

- **Automated embedding**: InsightEngine handles document chunking, embedding generation, and index maintenance natively. This replaces the custom `embedder.py` (BGE-small-en-v1.5 wrapper), `annotation_pipeline.py` (three-layer annotation), and `rag_engine.py` (retrieval-generation pipeline)---approximately 3,000 lines of custom Python code.
- **Native knowledge extraction**: ClinVar, AlphaMissense, and PubMed documents are ingested directly by InsightEngine, which extracts entities, relationships, and embeddings without the current multi-stage ETL pipeline (download → parse → annotate → embed → insert).
- **Built-in RAG orchestration**: InsightEngine manages the complete retrieve-augment-generate cycle, including context window management, citation tracking, and source attribution---capabilities currently implemented manually in the RAGEngine class.

**9.5.5 AgentEngine Replaces Custom Agent Frameworks**

VAST AgentEngine provides a native execution environment for autonomous agents:

- **CAR-T Intelligence Agent**: The CARTIntelligenceAgent's multi-step reasoning loop (search planning → evidence retrieval → evaluation → adaptive expansion → synthesis) maps to AgentEngine's native agent orchestration, with built-in tool use, memory management, and coordination primitives.
- **Drug Discovery Orchestration**: The 10-stage drug discovery pipeline becomes an AgentEngine-managed workflow with automatic checkpoint/resume, eliminating the custom JSON checkpoint system.
- **Meta-Agent Coordination**: The proposed meta-agent (Section 11.3.2) that routes clinical questions across specialist pipelines is a native AgentEngine pattern, not a custom implementation.

---

## 10. Current Integration Analysis

### 10.1 Existing Integration Points

The platform currently has four integration points:

1. **VCF Handoff** (Genomics → RAG): The VCF file produced by Stage 1 is manually ingested by the RAG pipeline's VCF parser. This is a one-time batch process, not a continuous data flow.

2. **Target Export** (RAG → Drug Discovery): The RAG pipeline's TargetHypothesisManager exports validated targets as JSON, which the Drug Discovery pipeline's TargetImporter consumes.

3. **Shared Milvus** (RAG ↔ CAR-T): Both systems connect to the same Milvus instance. The CAR-T Agent has read-only access to the `genomic_evidence` collection populated by the RAG pipeline.

4. **Shared Embedder** (RAG ↔ CAR-T): Both use identical BGE-small-en-v1.5 models, ensuring vector compatibility.

### 10.2 Integration Gaps

Despite sharing infrastructure, the platform has significant integration gaps:

**Gap 1: No Bidirectional Communication.** The pipeline flows strictly downward (FASTQ → VCF → Targets → Molecules). The Drug Discovery pipeline cannot request additional variant context from the RAG system, and the CAR-T Agent cannot trigger drug design for identified targets.

**Gap 2: Hardcoded Gene Database.** The Nextflow orchestrator's `RAG_IDENTIFY_TARGETS` process contains a hardcoded 11-gene database (VCP, BRCA1, TP53, EGFR, etc.) rather than dynamically querying Milvus for variant-driven target identification.

**Gap 3: Duplicate Code.** The RAG pipeline (`src/milvus_client.py`) and CAR-T Agent (`src/milvus_manager.py`) each implement their own Milvus connection management, embedding wrappers, and Claude API clients. Bug fixes in one don't propagate to the other.

**Gap 4: CAR-T Agent Isolation.** The CAR-T Agent operates as a standalone application with no presence in the Nextflow orchestrator, landing page dashboard, or Unified Portal. It cannot receive patient-specific genomic context from Stage 1.

**Gap 5: No Cross-Collection Querying.** The RAG pipeline cannot access CAR-T-specific knowledge (clinical trials, toxicity profiles, construct designs) when a user asks about CAR-T therapies for an identified target. Conversely, the CAR-T Agent cannot leverage ClinVar annotations or the Clinker knowledge graph.

**Gap 6: No Unified Report.** Each stage produces its own output (VCF, chat response, molecule rankings, therapy recommendations) with no mechanism to synthesize a single patient-facing clinical report.

**Gap 7: No CAR-T Engineering Path.** The Drug Discovery pipeline designs small molecules. CAR-T therapies are biologics requiring different computational approaches (scFv design, protein-protein docking, immunogenicity prediction). There is no computational path from CAR-T target selection to engineered construct optimization.

### 10.3 VAST AI OS Gap Resolution

VAST AI OS natively resolves or substantially mitigates six of the seven integration gaps identified above:

| Gap | Current Problem | VAST AI OS Resolution |
|-----|----------------|----------------------|
| **Gap 1**: No Bidirectional Communication | Pipeline flows one direction; no stage can request data from another | **DataEngine** event triggers enable automatic stage transitions. When Drug Discovery needs variant context, a DataEngine function queries DataBase directly. When the CAR-T Agent identifies a target, a DataEngine event triggers the drug design pipeline. Bidirectional communication becomes native event routing. |
| **Gap 2**: Hardcoded Gene Database | Static 11-gene lookup in Nextflow | **DataBase** SQL queries replace hardcoded lookups. Target identification becomes `SELECT gene, COUNT(*) FROM genomic_evidence WHERE clinical_significance = 'Pathogenic' GROUP BY gene ORDER BY COUNT(*) DESC LIMIT 5`---a live query against the actual patient data. |
| **Gap 3**: Duplicate Code | Two independent Milvus clients, embedding wrappers, LLM clients | **DataBase** native vector operations eliminate both custom Milvus clients (`milvus_client.py` and `milvus_manager.py`). **InsightEngine** replaces both embedding wrappers. The duplicate code problem disappears because the infrastructure layer provides the functionality directly. |
| **Gap 4**: CAR-T Isolation | CAR-T Agent has no pipeline integration | **DataEngine + AgentEngine** integrate the CAR-T Intelligence Agent as a native pipeline stage. AgentEngine manages the CAR-T reasoning loop as a first-class workflow step, receiving patient context from DataStore and publishing therapy recommendations back to DataStore. |
| **Gap 5**: No Cross-Collection Querying | RAG and CAR-T collections are searched independently | **DataBase** unified query plane enables cross-collection SQL+vector queries natively. A single query can search genomic evidence and CAR-T clinical trials simultaneously: `SELECT * FROM genomic_evidence ge JOIN cart_trials ct ON ge.gene = ct.target_antigen WHERE vector_distance(ge.embedding, ?) < 0.5 AND ct.phase IN ('PHASE_2', 'PHASE_3')`. |
| **Gap 6**: No Unified Report | Each stage produces separate outputs | **DataStore** stores all outputs in a single namespace with consistent metadata. A DataEngine function assembles the unified report by reading VCF annotations, target hypotheses, CAR-T evaluations, and drug candidates from DataStore---all in one data access pattern, not across four separate file systems and APIs. |

**Gap 7 (CAR-T Engineering Path)** remains a research challenge that requires new computational biology methods (scFv design, protein-protein docking) regardless of infrastructure. However, AgentEngine provides the orchestration framework for integrating these methods as they become available as NIM microservices.

---

## 11. Proposed Integration Research Plan

### 11.0 Phase 0: VAST AI OS Migration

**Objective**: Migrate the HCLS AI Factory's infrastructure layer from the current middleware stack (Milvus, Docker Compose, Nextflow, custom Python) to VAST AI OS native services, establishing the unified data operating system as the foundation for all subsequent integration phases.

**11.0.1 Migrate Vector Collections to DataBase**

Migrate all 12 Milvus collections (1 genomic_evidence + 10 CAR-T + 1 clinvar_variants) to VAST DataBase tables with native `VECTOR(384)` columns. Each collection's schema (documented in Appendix A) maps directly to a DataBase table with the same fields, replacing `FLOAT_VECTOR(384)` with the DataBase `VECTOR(384)` type. The IVF_FLAT index configuration (nlist=1024, nprobe=16) is replaced by DataBase's native ANN indexing. Validation: vector search recall and latency must match or exceed current Milvus benchmarks (Section 13.2).

**11.0.2 Replace Python In-Memory Dictionaries with DataBase Tables**

The ClinVarAnnotator loads variant_summary.txt.gz (2.7 million GRCh38 variants, ~4 GB RAM) into a Python dictionary keyed by `chrom_pos_ref_alt`. The AlphaMissenseAnnotator loads 71 million missense predictions (~8 GB RAM). Together, these consume 12 GB of Python process memory. Replace both with DataBase tables:

- `clinvar_annotations(chrom VARCHAR, pos INT64, ref VARCHAR, alt VARCHAR, clinical_significance VARCHAR, phenotype VARCHAR, rsid VARCHAR, hgvs VARCHAR, review_status VARCHAR)` with a composite index on `(chrom, pos, ref, alt)`.
- `alphamissense_predictions(chrom VARCHAR, pos INT64, ref VARCHAR, alt VARCHAR, am_score FLOAT, am_class VARCHAR)` with the same composite index.

Annotation becomes a DataBase JOIN instead of a Python dictionary lookup, eliminating 12 GB of RAM overhead and enabling the pipeline to annotate arbitrarily large VCFs without memory constraints.

**11.0.3 Replace Docker Compose Orchestration with DataEngine**

The current architecture deploys 15 services via Docker Compose, each with dedicated port assignments (Section 2.2), health check endpoints, and manual startup sequencing. Replace with DataEngine event-driven functions:

- **FASTQ arrival trigger**: When FASTQ files are written to DataStore, DataEngine automatically launches the Parabricks container on a CNode with GPU allocation.
- **VCF completion trigger**: When the VCF file is written, DataEngine triggers the annotation pipeline (ClinVar + AlphaMissense + VEP), followed by embedding generation and DataBase ingestion.
- **Target export trigger**: When target hypotheses are written to DataStore, DataEngine triggers the Drug Discovery pipeline.
- **CAR-T evaluation trigger**: When pathogenic variants are identified in hematological targets, DataEngine triggers the CAR-T Intelligence Agent for therapy evaluation.

This eliminates the Nextflow configuration file (Section 8.3), Docker Compose service definitions, and the Flask landing page's health check infrastructure. DataEngine provides native service health monitoring and automatic restart.

**11.0.4 Replace Custom RAG Pipeline with InsightEngine**

The RAGEngine class (Section 5.5) implements a five-step pipeline: query embedding, semantic search, therapeutic domain expansion, knowledge graph augmentation, and LLM synthesis. Replace with InsightEngine:

- InsightEngine ingests ClinVar, AlphaMissense, and the Clinker knowledge graph as knowledge sources.
- InsightEngine handles embedding generation, semantic retrieval, and context assembly natively.
- The therapeutic domain expansion maps (13 areas, 201 genes) are configured as InsightEngine query expansion rules.
- LLM synthesis uses InsightEngine's built-in RAG orchestration with Claude.

This eliminates approximately 3,000 lines of custom Python code in the RAG pipeline while preserving all current query capabilities.

**11.0.5 Replace Manual Pipeline Triggers with DataEngine S3 Events**

The current pipeline requires manual initiation: a user uploads FASTQ files, clicks "Run" in the Genomics Portal, waits for completion, then manually triggers RAG ingestion. Replace with DataEngine S3-compatible event auto-execution:

- `s3://patient-data/{patient_id}/fastq/*.fastq.gz` → triggers genomics pipeline
- `s3://patient-data/{patient_id}/vcf/*.vcf.gz` → triggers annotation + RAG ingestion
- `s3://patient-data/{patient_id}/targets/*.json` → triggers drug discovery
- `s3://patient-data/{patient_id}/variants/pathogenic_heme.json` → triggers CAR-T evaluation

The complete patient journey from FASTQ upload to therapy recommendation becomes fully automated.

**11.0.6 Enable AgentEngine for Autonomous Drug Discovery**

The Drug Discovery pipeline's 10-stage workflow (Section 7.2) and the CAR-T Agent's reasoning loop (Section 6.5) both implement custom orchestration patterns that map to AgentEngine:

- Each drug discovery stage becomes an AgentEngine tool (initialize, normalize, structure_discovery, molecule_generation, chemistry_qc, conformer_generation, docking, ranking, reporting).
- The CAR-T Agent's search planning, evidence retrieval, evaluation, and adaptive expansion steps become AgentEngine reasoning steps with built-in memory and tool use.
- The meta-agent coordination pattern (Section 11.3.2) is a native AgentEngine capability, enabling autonomous cross-pipeline reasoning without custom framework code.

### 11.1 Phase 1: Foundation --- Shared Component Consolidation

**Objective**: Eliminate duplicate code and establish the CAR-T Agent as a first-class citizen of the platform.

**1.1.1 Unified Python Package (`hcls_common/`)**

Extract a shared package with four modules:
- `milvus.py`: Unified connection manager supporting both RAG's 6 collections and CAR-T's 11 collections, with connection pooling
- `embedder.py`: BGE-small-en-v1.5 wrapper with embedding caching (LRU)
- `llm.py`: Claude API client with retry logic and token tracking
- `export.py`: PDF/Markdown/JSON export engine (ported from CAR-T's 1,264-line implementation to serve both pipelines)

Install as editable package in both virtual environments. This consolidation ensures that improvements to Milvus connection handling, embedding performance, or export quality propagate automatically.

**1.1.2 Landing Page Integration**

Add the CAR-T Agent to the landing page's service registry:
```python
SERVICES['cart-agent'] = {
    'port': 8521, 'name': 'CAR-T Intelligence Agent',
    'path': '/healthz', 'type': 'Pipeline'
}
```
Add a CAR-T card to the Portal UI with direct navigation.

**1.1.3 Export Capabilities for RAG Pipeline**

The RAG Chat Pipeline currently has no export functionality. Port the shared export engine and wire download buttons (Markdown, JSON, PDF) into the Streamlit chat UI, adapting evidence table schemas for genomic collections.

### 11.2 Phase 2: Pipeline Integration --- Connecting the Stages

**Objective**: Enable data to flow between stages without manual intervention.

**11.2.1 Dynamic Target Identification**

Replace the hardcoded gene database in the Nextflow orchestrator with a Python script that:
1. Queries the `genomic_evidence` and `clinvar_variants` collections for pathogenic/likely pathogenic variants
2. Cross-references against the 201-gene knowledge graph for druggability
3. Returns prioritized target hypotheses with variant evidence

This transforms target identification from static lookup to variant-driven discovery.

**11.2.2 Genomics → CAR-T Patient Stratification**

Build a VCF ingestion module for the CAR-T Agent that:
1. Reads the Parabricks VCF output
2. Identifies patient-specific tumor antigens via ClinVar cross-reference
3. Infers HLA types relevant to CAR-T immunogenicity
4. Feeds patient context into target antigen selection

This enables personalized CAR-T therapy recommendations grounded in individual genomic profiles.

**11.2.3 Bidirectional RAG ↔ CAR-T Knowledge Sharing**

Implement a shared query router that determines which collections to search based on query intent:

- RAG pipeline queries CAR-T's `cart_clinical_trials`, `cart_toxicity`, and `cart_target_antigens` when users ask about cell therapies for identified targets
- CAR-T Agent queries RAG's `clinvar_variants` and `clinker_kg` when evaluating target viability and variant pathogenicity

The router uses keyword and intent classification to select the optimal collection subset.

**11.2.4 CAR-T → Drug Discovery: CAR-T Engineering Branch**

Propose a "Stage 3b" for biologics engineering:
- **scFv Design**: Use construct data (scFv sequences, hinge/TM domains, costimulatory domains) from the CAR-T Agent's collections
- **Antigen-Antibody Docking**: Adapt DiffDock or integrate a protein-protein docking NIM for scFv-target interaction modeling
- **CAR-T Fitness Scoring**: Replace RDKit ADMET metrics with CAR-T-specific scores: binding affinity prediction, immunogenicity risk, cytokine storm probability

### 11.3 Phase 3: Intelligence Layer --- Agent Coordination

**Objective**: Enable autonomous cross-pipeline reasoning.

**11.3.1 Nextflow CAR-T Workflow Mode**

Add two workflow modes to the orchestrator:
- `cart`: VCF → Target Identification → CAR-T Agent evaluation → Therapy Report
- `full_cart`: FASTQ → VCF → Targets → CAR-T evaluation + Drug Discovery → Unified Report

This requires adding a REST API to the CAR-T Agent (currently Streamlit-only) and a corresponding `CART_EVALUATE` Nextflow process.

**11.3.2 Meta-Agent Orchestration**

Build a lightweight meta-agent that:
1. Receives high-level clinical questions (e.g., *"What are the best treatment options for this patient's B-ALL?"*)
2. Routes sub-questions to the appropriate specialist: genomic context → RAG, therapy options → CAR-T, molecule design → Drug Discovery
3. Synthesizes a unified clinical narrative combining genomic findings, CAR-T recommendations, and drug candidates

This uses Claude's tool-use capability to orchestrate sub-agents, with each existing pipeline exposed as a callable tool.

**11.3.3 Unified Clinical Report**

Extend the PDF export engine to generate a comprehensive patient report:

- **Section 1: Genomic Profile** --- Key variants from VCF, ClinVar annotations, AlphaMissense scores
- **Section 2: Target Analysis** --- Identified therapeutic targets, evidence strength, druggability
- **Section 3: CAR-T Therapy Evaluation** --- Recommended constructs, toxicity profiles, clinical trial matches
- **Section 4: Small Molecule Candidates** --- MolMIM-generated molecules, DiffDock binding scores, QED
- **Section 5: Treatment Recommendation** --- Meta-agent synthesized recommendation with cross-pipeline evidence

---

## 12. Implementation Roadmap

### 12.1 Priority Matrix

| Phase | Item | Impact | Effort | Priority |
|-------|------|--------|--------|----------|
| 1 | Shared `hcls_common` package | High | Medium | **P0** |
| 1 | Landing page + Portal integration | Medium | Low | **P0** |
| 1 | Export to RAG pipeline | Medium | Low | **P1** |
| 2 | Dynamic target identification | High | Medium | **P1** |
| 2 | Genomics → CAR-T patient stratification | Very High | High | **P1** |
| 2 | Cross-collection query router | High | Medium | **P1** |
| 2 | CAR-T Engineering branch (Stage 3b) | Very High | Very High | **P2** |
| 3 | Nextflow CAR-T workflow mode + REST API | High | Medium | **P2** |
| 3 | Unified clinical report | High | High | **P2** |
| 3 | Meta-agent orchestration | Very High | Very High | **P3** |
| 4 | Real NIM service wiring | Medium | Low | **P1** |
| 4 | Connection pooling + embedding cache | Medium | Medium | **P2** |
| 4 | OpenTelemetry tracing | Medium | High | **P3** |

### 12.2 Sprint Plan

**Sprint 1: Foundation** (Items 1.1, 1.2, 4.1)
Get shared code right. Consolidate Milvus/embedding/LLM clients, add CAR-T to landing page, wire real NIM endpoints.

**Sprint 2: Cross-Pollination** (Items 1.3, 2.1, 2.3)
Export everywhere, dynamic target identification, cross-collection queries.

**Sprint 3: Patient-Specific Flow** (Items 2.2, 3.1)
VCF → CAR-T ingestion, Nextflow CAR-T workflow mode.

**Sprint 4: Engineering and Reporting** (Items 2.4, 3.3, 4.2)
CAR-T Engineering branch, unified report, performance optimization.

**Sprint 5: Intelligence** (Items 3.2, 4.3)
Meta-agent coordination, full observability stack.

---

## 13. Performance Benchmarks and Metrics

### 13.1 End-to-End Timing

| Stage | Duration | Data Volume |
|-------|----------|-------------|
| Genomics (fq2bam) | 60-180 min | 200GB FASTQ → 120GB BAM |
| Genomics (DeepVariant) | 60-90 min | 120GB BAM → 3GB VCF |
| VCF Annotation | 5-15 min | 11.7M variants annotated |
| Vector Ingestion | 10-20 min | 35,616 variants embedded |
| RAG Query | 1-5 sec | 10 results per query |
| CAR-T Agent Query | 2-8 sec | 30 results across 11 collections |
| Drug Discovery | 8-16 min | 20 molecules generated, docked, ranked |
| **Total** | **<5 hours** | Patient DNA → Drug Candidates |

### 13.2 Vector Search Performance

| Metric | Value |
|--------|-------|
| Single-collection search latency | <100ms |
| 11-collection parallel search | <500ms |
| Embedding generation (per query) | ~10ms |
| Top-k per collection | 5 |
| Max merged results | 30 |
| Score threshold | 0.3-0.5 (configurable) |

### 13.3 GPU Utilization Profile

| Stage | GPU Utilization | GPU Memory |
|-------|----------------|------------|
| Parabricks fq2bam | 80-95% | 12-16GB |
| Parabricks DeepVariant | 70-85% | 10-14GB |
| BGE Embedding | 10-20% | <1GB |
| MolMIM Generation | 40-60% | 4-8GB |
| DiffDock Docking | 50-70% | 4-8GB |

### 13.4 Traditional vs. HCLS AI Factory Comparison

| Task | Traditional | HCLS AI Factory | Speedup |
|------|-------------|-----------------|---------|
| Variant Calling (WGS) | 24-72 hours (CPU) | 2-4 hours (GPU) | 10-30x |
| Variant Interpretation | Days (manual curation) | Seconds (RAG) | 1000x+ |
| Lead Discovery | 6-12 months | 8-16 minutes | 30,000x+ |
| Lead Optimization | 1-2 years | Minutes (iterative) | 100,000x+ |
| CAR-T Target Evaluation | Weeks (literature review) | 2-8 seconds | 100,000x+ |

---

## 14. Security, Observability, and Production Hardening

### 14.1 Security Architecture

**Authentication**: Optional API key authentication (X-API-Key header) on the Genomics Portal. The CAR-T Agent and RAG Chat UI are currently session-based without authentication, suitable for single-user DGX Spark deployment.

**Input Validation**: Shell metacharacter filtering (`[$\`|;&]`), path traversal prevention via `secure_filename`, Milvus expression injection prevention via regex-validated gene names and chromosome identifiers.

**Security Headers**: X-Frame-Options: DENY, X-Content-Type-Options: nosniff, X-XSS-Protection, Content-Security-Policy, Permissions-Policy.

**Rate Limiting**: 60 requests per 60-second window per IP address on the Genomics Portal.

### 14.2 Observability

Current monitoring covers hardware metrics (GPU utilization, temperature, power, memory via DCGM and Node Exporter) visualized in Grafana dashboards. The proposed integration plan adds:

- **Application-level metrics**: Query latency histograms, collection hit rates, LLM token consumption, error rates per service
- **Distributed tracing**: OpenTelemetry spans across Nextflow → Genomics → RAG → CAR-T → Drug Discovery for end-to-end request tracking
- **Custom dashboards**: Pipeline throughput (queries/hour), per-agent performance, vector database growth curves

### 14.3 Error Handling

Each stage implements independent error recovery:

- **Genomics**: 6-retry BAM verification (60s total), 3-retry DeepVariant (30s delays), custom nvidia-smi wrapper for DGX Spark quirks
- **RAG Pipeline**: Exponential backoff on LLM API calls, graceful degradation to local Ollama models when Claude API is unavailable
- **CAR-T Agent**: ThreadPoolExecutor exception isolation per collection, adaptive evidence expansion on insufficient results
- **Drug Discovery**: 10-stage checkpoint/resume system, 5-retry NIM calls with exponential backoff (2s, 4s, 8s, 16s, 32s), 3 consecutive docking failures trigger abort
- **Orchestrator**: Exit code-based retry [143, 137, 104, 134, 139] with max 2 retries, optional email notification on completion/failure

---

## 15. Discussion

### 15.1 Architectural Decisions

**Single Machine vs. Distributed**: We chose single-machine deployment (DGX Spark) over distributed cloud architecture for three reasons: (a) data sovereignty---genomic data remains on-premises, (b) cost---$4,699 one-time vs. ongoing cloud compute charges, and (c) simplicity---no cluster management, network partitioning, or distributed consensus overhead. The 128GB unified memory and NVLink-C2C interconnect provide sufficient bandwidth for sequential stage execution.

**Milvus vs. Alternatives**: We selected Milvus over Pinecone, Weaviate, and Chroma because: (a) fully self-hosted (no cloud dependency), (b) IVF_FLAT indexing provides exact-in-partition search suitable for our collection sizes (6K-36K vectors), (c) rich scalar filtering (essential for antigen, phase, and stage filters), and (d) mature Python client with ThreadPoolExecutor-compatible API.

**BGE-small-en-v1.5 vs. Larger Models**: We chose the 384-dimensional BGE-small model over BGE-large (1024d) or domain-specific biomedical embeddings because: (a) 3x smaller vectors reduce Milvus memory footprint, (b) sufficient discrimination for our collection sizes, (c) fast inference on ARM CPU (DGX Spark), and (d) the query expansion system compensates for reduced embedding capacity through explicit term injection.

**VAST AI OS vs. Current Middleware Stack**: The V1.0 architecture validated the HCLS AI Factory's pipeline logic and clinical utility on commodity hardware. Version 2.0 recognizes that the middleware stack itself---Milvus, Docker Compose, Nextflow, custom RAG code, Python in-memory dictionaries---has become the primary scaling constraint. A unified data operating system is superior to the current middleware stack for three fundamental reasons:

1. **The Egress/Ingress Tax**: Processing a single patient genome through the full pipeline involves approximately 17 GB of inter-service data movement. FASTQ reads (200 GB) are copied from upload location to local NVMe. The BAM file (120 GB) is written by Parabricks and re-read by samtools. The VCF (3 GB) is parsed by Python, annotated, serialized to embeddings, and inserted into Milvus via TCP. RAG context is assembled in Python memory, serialized to JSON, and sent to the Claude API. Drug discovery results move between MolMIM, DiffDock, and RDKit through HTTP APIs. On VAST AI OS, this data movement goes to zero: DataStore provides zero-copy access, DataBase co-locates vector and SQL queries with the data, and DataEngine routes pipeline stages to the data rather than moving data to the pipeline.

2. **The Middleware Maintenance Burden**: The current architecture requires expertise in six distinct infrastructure technologies (Milvus administration, Docker networking, Nextflow DSL2 syntax, Flask deployment, Prometheus/Grafana configuration, and Python embedding pipelines). Each has its own versioning, configuration, failure modes, and upgrade paths. VAST AI OS consolidates these into a single platform with unified management, monitoring, and upgrade cycles. For a clinical deployment, this reduction in operational complexity is not a convenience---it is a regulatory requirement. FDA 21 CFR Part 11 compliance is substantially easier when the infrastructure layer is a single, validated platform rather than a manually assembled stack.

3. **The Scale Ceiling**: The DGX Spark deployment serves one patient at a time. Clinical-grade precision medicine requires processing hundreds of patients per day with multi-tenant isolation, audit trails, and guaranteed latency. The current architecture cannot scale to this level without fundamental re-architecture: Milvus requires sharding, Docker Compose requires migration to Kubernetes, Nextflow requires HPC queue integration, and file storage requires migration to a distributed file system. VAST AI OS's DASE architecture (stateless CNodes, stateful DBoxes) provides this scale path natively, with DataSpace providing multi-tenant isolation and SyncEngine enabling multi-site replication.

### 15.2 Limitations

1. **No Real-Time VCF Streaming**: The genomics pipeline produces a complete VCF before the RAG pipeline can begin annotation. Streaming variant ingestion would reduce end-to-end latency.

2. **Static Collection Data**: CAR-T collections are populated by batch ingestion scripts. There is no automated pipeline for continuous literature updates from PubMed, ClinicalTrials.gov, or FDA safety databases.

3. **Small Molecule Bias**: Stage 3 targets small molecule drug discovery. CAR-T therapies are biological products requiring fundamentally different computational approaches (protein engineering, not SMILES manipulation).

4. **Single-Patient Focus**: The system processes one patient at a time. Clinical deployment would require queue management and multi-tenant isolation.

5. **Validation Gap**: While DeepVariant achieves >99% accuracy on GIAB benchmarks, the end-to-end clinical validity of the target identification → therapy recommendation chain has not been validated in a clinical trial setting.

### 15.3 Comparison with Related Work

The HCLS AI Factory differs from existing platforms in its end-to-end scope. Systems like Terra/AnVIL focus on genomics orchestration, Watson for Genomics on variant interpretation, and OpenTargets on target identification---but none integrate GPU-accelerated variant calling, RAG-powered interpretation, domain-specific cell therapy intelligence, and de novo drug design into a single, deployable system. The closest comparable is Illumina's DRAGEN + Emedgene pipeline, which covers genomics through interpretation but lacks drug discovery and cell therapy evaluation.

---

## 16. Future Directions

### 16.1 Short-Term (3-6 Months): VAST AI OS Foundation

- **Migrate to DataBase**: Move all 12 Milvus collections to VAST DataBase native VECTOR columns; replace Python in-memory ClinVar/AlphaMissense dictionaries with DataBase SQL tables (eliminating 12 GB RAM overhead)
- **DataEngine event triggers**: Replace Docker Compose service management and Nextflow scheduling with DataEngine event-driven pipeline execution (FASTQ arrival → automatic genomics → automatic RAG → automatic drug discovery)
- **InsightEngine RAG**: Replace custom RAG pipeline code (~3,000 lines) with InsightEngine native knowledge extraction and retrieval-augmented generation
- **Shared component consolidation** (Phase 1 of research plan): Unified `hcls_common` Python package, landing page integration, export capabilities
- **Automated PubMed/ClinicalTrials.gov ingestion** via InsightEngine continuous knowledge ingestion for CAR-T collection updates

### 16.2 Medium-Term (6-12 Months): Intelligent Orchestration

- **AgentEngine orchestration**: Migrate the CAR-T Intelligence Agent reasoning loop and Drug Discovery 10-stage pipeline to AgentEngine native agent execution with built-in tool use, memory, and coordination
- **DataSpace multi-tenant isolation**: Enable per-patient, per-research-team, and per-clinical-workflow data isolation with cryptographic guarantees---supporting clinical deployment with multiple concurrent patients
- **ICMS acceleration**: Deploy BlueField-4 DPU-based petabyte-scale KV cache for 10-20x faster LLM time-to-first-token on RAG synthesis queries
- **Patient-specific CAR-T recommendations** grounded in individual VCF data via DataBase cross-table joins
- **Cross-collection query routing** via DataBase unified SQL+vector queries (replacing manual multi-collection Python orchestration)
- **CAR-T Engineering branch** (Stage 3b) for scFv design and protein-protein docking
- **Unified clinical report** assembled by DataEngine from DataStore outputs spanning all four stages

### 16.3 Long-Term (12-24 Months): Platform Scale

- **SyncEngine multi-site replication**: Enable multi-center clinical deployments with real-time data replication across sites---supporting federated genomics without custom synchronization code
- **DataSpace for federated genomics**: Secure data sharing across institutions with DataSpace tenant isolation, enabling collaborative research on pooled genomic datasets without data movement
- **BlueField-4 GPUDirect data path**: Zero-copy genomic data access from DataStore to GPU memory via GPUDirect RDMA, eliminating the PCIe/NVLink staging overhead for Parabricks and BioNeMo workloads
- **DASE elastic scaling**: Scale from single DGX Spark to multi-rack deployment by adding CNodes (compute) and DBoxes (data) independently, with DataEngine automatically distributing pipeline execution across available resources
- **Meta-agent orchestration** with Claude tool-use for autonomous cross-pipeline reasoning via AgentEngine
- **Federated learning** across multiple VAST AI OS deployments for collaborative model improvement
- **Clinical trial integration** for prospective validation of treatment recommendations
- **Imaging intelligence agent** for pathology slide analysis (H&E, IHC)
- **Precision biomarker agent** for longitudinal patient monitoring
- **Regulatory pathway automation** for IND/BLA submission support, with DataStore audit trails providing 21 CFR Part 11 compliance evidence

---

## 17. Conclusion

The HCLS AI Factory demonstrates that a comprehensive precision medicine platform---from raw sequencing data to therapeutic candidates---can operate on a single $4,699 workstation in under five hours. By integrating GPU-accelerated genomics (Parabricks), RAG-powered variant interpretation (Milvus + Claude), autonomous CAR-T therapy intelligence (11 specialized collections with agent reasoning), and AI-driven drug discovery (BioNeMo NIMs), the system bridges the traditional silos between bioinformatics, clinical genetics, immunology, and medicinal chemistry.

The Version 1.0 architecture validated the pipeline logic and clinical utility. Version 2.0 recognizes the deeper truth: an AI Factory that must ingest extensive amounts of sequencing data, trigger GPU pipelines, embed millions of records into vector space, run RAG across clinical evidence, orchestrate autonomous drug design agents, and deliver it all at clinical-grade latency is not solving a storage problem. It is solving an operating system problem.

VAST AI OS provides that operating system. DataStore replaces fragmented file systems with a single petabyte-scale namespace offering zero-copy GPU access. DataBase unifies vector search and SQL queries, eliminating standalone Milvus and 12 GB of Python in-memory dictionaries. DataEngine replaces Nextflow and Docker Compose with event-driven pipeline orchestration on stateless CNodes. InsightEngine replaces 3,000+ lines of custom RAG code with native knowledge extraction and retrieval-augmented generation. AgentEngine provides a native execution environment for the CAR-T Intelligence Agent and drug discovery orchestration. ICMS on BlueField-4 accelerates LLM inference with petabyte-scale KV caching. SyncEngine enables multi-site clinical deployment. The DASE architecture ensures that this foundation scales from a single DGX Spark to a multi-rack clinical deployment without re-architecture.

The proposed integration research plan identifies concrete, prioritized opportunities to transform these semi-independent components into a fully interconnected platform. Phase 0 migrates the infrastructure layer to VAST AI OS. Phase 1 consolidates shared application code. Phase 2 enables bidirectional data flow between stages. Phase 3 introduces agent-level coordination for autonomous cross-pipeline reasoning.

The HCLS AI Factory is being open-sourced at VAST for AI OS. The Apache 2.0 license and commodity hardware target ensure that this capability is accessible to academic medical centers, clinical research organizations, and biopharmaceutical companies. With VAST AI OS as the infrastructure foundation, these organizations gain not just the pipeline logic, but a production-grade data operating system that eliminates the middleware stack, the egress/ingress tax, and the operational complexity that have historically prevented computational precision medicine from reaching clinical deployment at scale.

---

## 18. References

1. NVIDIA Clara Parabricks Documentation. *Clara Parabricks 4.6.0 User Guide*. NVIDIA Corporation, 2025.

2. Zook, J.M. et al. "An open resource for accurately benchmarking small variant and reference calls." *Nature Biotechnology*, 37(5):561-566, 2019. (GIAB HG002)

3. Cheng, J. et al. "Accurate proteome-wide missense variant effect prediction with AlphaMissense." *Science*, 381(6664):eadg7492, 2023.

4. Landrum, M.J. et al. "ClinVar: improvements to accessing data." *Nucleic Acids Research*, 48(D1):D835-D844, 2020.

5. June, C.H. et al. "CAR T cell immunotherapy for human cancer." *Science*, 359(6382):1361-1365, 2018.

6. Xiao, S. et al. "C-Pack: Packaged Resources To Advance General Chinese Embedding." *arXiv:2309.07597*, 2023. (BGE-small-en-v1.5)

7. Wang, J. et al. "Milvus: A Purpose-Built Vector Data Management System." *SIGMOD*, 2021.

8. Corso, G. et al. "DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking." *ICLR*, 2023.

9. Di Tommaso, P. et al. "Nextflow enables reproducible computational workflows." *Nature Biotechnology*, 35(4):316-319, 2017.

10. Bickels, R.D. et al. "Quantitative Estimation of Drug-Likeness." *Nature Chemistry*, 4(2):90-98, 2012. (QED)

11. Lipinski, C.A. et al. "Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings." *Advanced Drug Delivery Reviews*, 23(1-3):3-25, 1997.

12. Poplin, R. et al. "A universal SNP and small-indel variant caller using deep neural networks." *Nature Biotechnology*, 36(10):983-987, 2018. (DeepVariant)

13. McLennan, S. et al. "Embedded AI for personalised medicine: bridging the translational gap." *Nature Medicine*, 29:1557-1558, 2023.

14. Lewis, P. et al. "Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks." *NeurIPS*, 2020.

15. Sterner, R.C. & Sterner, R.M. "CAR-T cell therapy: current limitations and potential strategies." *Blood Cancer Journal*, 11(4):69, 2021.

16. VAST Data, Inc. *VAST AI OS Documentation: DataStore, DataBase, DataEngine, InsightEngine, AgentEngine, ICMS, SyncEngine, DataSpace*. VAST Data Platform Documentation, 2025-2026. https://vastdata.com/platform

17. VAST Data, Inc. "The VAST Data Platform: A Disaggregated Shared-Everything (DASE) Architecture for AI-Native Data Infrastructure." *VAST Data Technical White Paper*, 2025.

---

## 19. Appendices

### Appendix A: Complete Collection Schemas

#### A.1 CAR-T Literature Collection (`cart_literature`)

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR | 100 | PMID or patent number |
| embedding | FLOAT_VECTOR | 384 | BGE-small embedding |
| title | VARCHAR | 500 | Paper/patent title |
| text_chunk | VARCHAR | 3000 | Text for embedding |
| source_type | VARCHAR | 20 | pubmed/pmc/patent/preprint/manual |
| year | INT64 | - | Publication year |
| cart_stage | VARCHAR | 30 | target_id/car_design/vector_eng/testing/clinical |
| target_antigen | VARCHAR | 100 | CD19, BCMA, etc. |
| disease | VARCHAR | 200 | Disease/indication |
| keywords | VARCHAR | 1000 | MeSH terms |
| journal | VARCHAR | 200 | Journal name |

#### A.2 Clinical Trials Collection (`cart_trials`)

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR | 20 | NCT number (^NCT\d{8}$) |
| embedding | FLOAT_VECTOR | 384 | BGE-small embedding |
| title | VARCHAR | 500 | Official trial title |
| text_summary | VARCHAR | 3000 | Summary for embedding |
| phase | VARCHAR | 30 | TrialPhase enum |
| status | VARCHAR | 30 | TrialStatus enum |
| sponsor | VARCHAR | 200 | Lead sponsor |
| target_antigen | VARCHAR | 100 | Target antigen |
| car_generation | VARCHAR | 20 | 1st/2nd/3rd/4th/armored/universal |
| costimulatory | VARCHAR | 50 | CD28, 4-1BB, dual |
| disease | VARCHAR | 200 | Indication |
| enrollment | INT64 | - | Target enrollment count |
| start_year | INT64 | - | Study start year |
| outcome_summary | VARCHAR | 2000 | Outcome summary |

#### A.3 CAR Constructs Collection (`cart_constructs`)

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR | 100 | Construct identifier |
| embedding | FLOAT_VECTOR | 384 | BGE-small embedding |
| name | VARCHAR | 200 | Product/construct name |
| text_summary | VARCHAR | 2000 | Description |
| target_antigen | VARCHAR | 100 | Target |
| scfv_origin | VARCHAR | 200 | Antibody clone/origin |
| costimulatory_domain | VARCHAR | 100 | CD28, 4-1BB, dual |
| signaling_domain | VARCHAR | 100 | CD3-zeta |
| generation | VARCHAR | 20 | CARGeneration enum |
| hinge_tm | VARCHAR | 200 | Hinge + transmembrane |
| vector_type | VARCHAR | 50 | lentiviral/retroviral/non-viral |
| fda_status | VARCHAR | 20 | FDAStatus enum |
| known_toxicities | VARCHAR | 500 | CRS, ICANS, etc. |

#### A.4 Safety Collection (`cart_safety`)

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR | 100 | Safety record ID |
| embedding | FLOAT_VECTOR | 384 | BGE-small embedding |
| text_summary | VARCHAR | 3000 | Safety summary |
| product | VARCHAR | 200 | Product name |
| event_type | VARCHAR | 30 | CRS/ICANS/CYTOPENIA/etc. |
| severity_grade | VARCHAR | 100 | Grade 1-5 |
| onset_timing | VARCHAR | 100 | Median day post-infusion |
| incidence_rate | VARCHAR | 200 | e.g., 42% any grade |
| management_protocol | VARCHAR | 500 | Management approach |
| outcome | VARCHAR | 100 | Recovery status |
| reporting_source | VARCHAR | 50 | FAERS/trial/registry/label |
| year | INT64 | - | Reporting year |

#### A.5 Remaining Collections

Assays (`cart_assays`), Manufacturing (`cart_manufacturing`), Biomarkers (`cart_biomarkers`), Regulatory (`cart_regulatory`), Sequences (`cart_sequences`), Real-World Evidence (`cart_realworld`), and Genomic Evidence (`genomic_evidence`) follow similar patterns with domain-specific fields as detailed in Section 6.2.

### Appendix B: Enumeration Types

```
CARTStage: TARGET_ID, CAR_DESIGN, VECTOR_ENG, TESTING, CLINICAL
SourceType: PUBMED, PMC, PATENT, PREPRINT, MANUAL
TrialPhase: EARLY_1, PHASE_1, PHASE_1_2, PHASE_2, PHASE_2_3, PHASE_3, PHASE_4, NA
TrialStatus: RECRUITING, ACTIVE, COMPLETED, TERMINATED, WITHDRAWN, SUSPENDED, NOT_YET, UNKNOWN
CARGeneration: FIRST, SECOND, THIRD, FOURTH, ARMORED, UNIVERSAL
AssayType: CYTOTOXICITY, CYTOKINE, FLOW_CYTOMETRY, PROLIFERATION, IN_VIVO, PERSISTENCE, EXHAUSTION
ProcessStep: TRANSDUCTION, EXPANSION, HARVEST, FORMULATION, RELEASE_TESTING, CRYOPRESERVATION
FDAStatus: APPROVED, BLA_FILED, PHASE_3, PHASE_2, PHASE_1, PRECLINICAL, DISCONTINUED
SafetyEventType: CRS, ICANS, CYTOPENIA, INFECTION, SECONDARY_MALIGNANCY, ORGAN_TOXICITY, NEUROLOGIC, CARDIAC
BiomarkerType: PREDICTIVE, PROGNOSTIC, PHARMACODYNAMIC, MONITORING, RESISTANCE
EvidenceLevel: VALIDATED, EMERGING, EXPLORATORY
RegulatoryEvent: BLA, BREAKTHROUGH_THERAPY, RMAT, ACCELERATED_APPROVAL, FULL_APPROVAL, LABEL_UPDATE, REMS, POST_MARKETING_REQ, COMPLETE_RESPONSE
RWEStudyType: RETROSPECTIVE, REGISTRY, CLAIMS, EHR_ANALYSIS, META_ANALYSIS
```

### Appendix C: Service Port Map

```
 8080 = Landing Page Dashboard (Flask)
 5000 = Genomics Portal (Flask)
 5001 = RAG/Chat API (Flask/FastAPI)
 8501 = RAG Chat UI (Streamlit)
 8505 = Drug Discovery UI (Streamlit)
 8510 = Discovery Portal (Streamlit)
 8521 = CAR-T Intelligence Agent (Streamlit)
19530 = Milvus Vector Database
 8000 = Attu (Milvus Admin UI)
 8001 = MolMIM NIM (BioNeMo)
 8002 = DiffDock NIM (BioNeMo)
 3000 = Grafana Dashboard
 9099 = Prometheus
 9100 = Node Exporter
 9400 = DCGM GPU Exporter
```

### Appendix D: Software Dependencies

**Shared Across Pipelines:**
- Python 3.10+, pydantic 2.12+, loguru 0.7+, python-dotenv 1.2+
- pymilvus 2.6+, sentence-transformers 5.2+, torch 2.9+
- anthropic 0.75+, streamlit 1.52+
- Flask 3.1+, Flask-CORS 6.0+

**Genomics-Specific:**
- Docker 20.10+, NVIDIA Container Runtime
- nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1
- psutil 7.2+, pynvml 13.0+

**RAG-Specific:**
- cyvcf2 0.31+, pysam 0.23+ (VCF processing)
- openai 2.15+ (multi-provider support)
- fastapi 0.128+, uvicorn 0.40+ (API server)

**CAR-T-Specific:**
- reportlab 4.4+ (PDF export)

**Drug Discovery-Specific:**
- rdkit 2025.9+ (chemistry toolkit)
- typer 0.21+, rich 14.2+ (CLI)
- stmol 0.0.9+, py3Dmol 2.5+ (3D molecular viewer)
- pillow 12.1+ (image processing)

---

*This paper describes the HCLS AI Factory as of February 2026 (Version 2.0). The system is open-source under Apache 2.0 license and is being released at VAST for AI OS.*

*Correspondence: Adam Jones, HCLS AI Factory Project*
