# Precision Oncology Agent — Project Bible

> **Purpose:** Complete implementation reference for building the HCLS Precision Oncology Agent on NVIDIA DGX Spark. Import this document into a Claude Code session as context for implementation.
>
> **License:** Apache 2.0 | **Author:** Adam Jones | **Date:** February 2026

---

## Table of Contents

1. [Project Overview & Goals](#1-project-overview--goals)
2. [DGX Spark Hardware Reference](#2-dgx-spark-hardware-reference)
3. [Repository Layout](#3-repository-layout)
4. [Docker Compose Services](#4-docker-compose-services)
5. [Milvus Collection Schemas](#5-milvus-collection-schemas)
6. [Pydantic Data Models](#6-pydantic-data-models)
7. [Configuration Reference](#7-configuration-reference)
8. [Embedding Strategy](#8-embedding-strategy)
9. [Multi-Collection RAG Engine](#9-multi-collection-rag-engine)
10. [Knowledge Graph](#10-knowledge-graph)
11. [Query Expansion](#11-query-expansion)
12. [Comparative Analysis Mode](#12-comparative-analysis-mode)
13. [Case Manager Workflow](#13-case-manager-workflow)
14. [Trial Matcher Implementation](#14-trial-matcher-implementation)
15. [Therapy Ranker Rules](#15-therapy-ranker-rules)
16. [Cross-Modal Integration](#16-cross-modal-integration)
17. [Data Ingest Pipelines](#17-data-ingest-pipelines)
18. [FastAPI REST Server](#18-fastapi-rest-server)
19. [Streamlit MTB Workbench](#19-streamlit-mtb-workbench)
20. [Export Formats](#20-export-formats)
21. [Monitoring & Metrics](#21-monitoring--metrics)
22. [Testing Strategy](#22-testing-strategy)
23. [HCLS AI Factory Integration](#23-hcls-ai-factory-integration)
24. [Implementation Sequence](#24-implementation-sequence)

---

## 1. Project Overview & Goals

### What This Agent Does

The Precision Oncology Agent is a closed-loop clinical decision support system that transforms raw genomic data (VCF files) into actionable Molecular Tumor Board (MTB) packets. It combines variant annotation, evidence retrieval, therapy ranking, clinical trial matching, and outcomes learning across **11 Milvus collections** to deliver grounded, citation-backed recommendations for precision cancer treatment.

### Clinical Workflow

| Step | Description |
|---|---|
| **VCF Upload** | Raw VCF text parsed, extracting PASS variants with gene, consequence, and position |
| **Variant Annotation** | Each variant classified against ~40 actionable targets using AMP/ASCO/CAP evidence tiers (A-D) |
| **Evidence Lookup** | RAG retrieval for each actionable variant across literature, therapies, and guidelines |
| **Therapy Ranking** | Evidence-level-sorted therapy recommendations with resistance flags and contraindication checks |
| **Trial Matching** | Hybrid deterministic + semantic search against oncology clinical trials |
| **Open Questions** | VUS variants, missing biomarkers, and evidence gaps flagged for MTB discussion |

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | **~1,490** across 10 owned collections + **3.5M** genomic (read-only) |
| Multi-collection search latency | **< 200 ms** (11 collections, top-5 each) |
| Comparative dual retrieval | **~400 ms** (2 x 11 collections, entity-filtered) |
| Full RAG query (search + Claude) | **~24 sec** end-to-end |
| MTB packet generation | **< 30 sec** (full workflow) |
| Cosine similarity scores | **0.72 - 0.92** on demo queries |
| Knowledge graph entities | **135+** (~40 targets, ~30 therapies, ~20 resistance, ~10 pathways, ~15 biomarkers, ~50 aliases) |
| Query expansion | **12 maps**, ~120 keywords -> ~700 terms |

### Pipeline Pattern

Every query follows the same pattern:

1. User question submitted (Streamlit UI or REST API)
2. Agent plans search strategy (identify genes, cancer types, decompose complex queries)
3. BGE-small embeds query with asymmetric instruction prefix
4. Parallel search across 11 Milvus collections (ThreadPoolExecutor)
5. Query expansion adds semantically related terms
6. Knowledge graph augments results with domain context
7. Evidence merged, deduplicated, weighted-ranked (cap at 30)
8. Claude Sonnet 4.6 generates grounded answer with citations
9. Response streamed to user with evidence panel

### HCLS AI Factory Integration

This agent is one node in the broader HCLS AI Factory. Cross-agent triggers include:

- **Oncology -> Genomics (Parabricks):** Shared Milvus instance, genomic_evidence collection (3.5M vectors, read-only)
- **Oncology -> Drug Discovery (BioNeMo):** Actionable targets feed molecule generation (MolMIM + DiffDock)
- **Oncology -> Imaging Agent:** Cross-modal triggers from imaging findings to genomic correlates
- **Oncology -> CAR-T Agent:** Shared target biology and biomarker intelligence

---

## 2. DGX Spark Hardware Reference

### Specifications

| Parameter | Value |
|---|---|
| GPU | NVIDIA GB10 Grace Blackwell Superchip |
| Architecture | Blackwell (GPU) + Grace (CPU) via NVLink-C2C |
| GPU Memory | 128 GB unified LPDDR5x (shared CPU+GPU) |
| CPU | 20 ARM Neoverse V2 cores (Grace) |
| Storage | NVMe SSD (local) |
| Networking | 10GbE, WiFi 7, Bluetooth 5.3 |
| Price | $3,999 |
| Power | Desktop form factor, standard AC |

### Why DGX Spark

- **128 GB unified memory** eliminates CPU<->GPU transfer bottleneck
- **NVLink-C2C** provides 900 GB/s bandwidth between Grace CPU and Blackwell GPU
- **ARM64 Grace CPU** runs all Python services natively
- **Single-box deployment** -- Milvus, embeddings, LLM inference, and UI on one machine
- **$3,999 price point** -- accessible for research labs and departmental pilots

---

## 3. Repository Layout

```
precision_oncology_agent/agent/
├── src/
│   ├── models.py                  # Pydantic data models (497 lines)
│   ├── collections.py             # 11 Milvus collection schemas + manager (606 lines)
│   ├── knowledge.py               # Knowledge graph: targets, therapies, resistance (1,194 lines)
│   ├── query_expansion.py         # 12 expansion maps (676 lines)
│   ├── rag_engine.py              # Multi-collection RAG + comparative analysis (780 lines)
│   ├── agent.py                   # Plan-search-synthesize pipeline (489 lines)
│   ├── case_manager.py            # VCF parsing + MTB packet generation (509 lines)
│   ├── trial_matcher.py           # Hybrid deterministic + semantic matching (393 lines)
│   ├── therapy_ranker.py          # Evidence-based therapy ranking (552 lines)
│   ├── cross_modal.py             # Cross-modal triggers to imaging + drug discovery (395 lines)
│   ├── export.py                  # Markdown, JSON, PDF, FHIR R4 export (876 lines)
│   ├── metrics.py                 # Prometheus metrics (362 lines)
│   ├── scheduler.py               # Data ingestion scheduler (263 lines)
│   ├── ingest/
│   │   ├── base.py                # Base ingest pipeline (249 lines)
│   │   ├── civic_parser.py        # CIViC actionable variant ingest (340 lines)
│   │   ├── oncokb_parser.py       # OncoKB data parser (104 lines)
│   │   ├── literature_parser.py   # PubMed E-utilities ingest (248 lines)
│   │   ├── clinical_trials_parser.py  # ClinicalTrials.gov API v2 (279 lines)
│   │   ├── guideline_parser.py    # NCCN/ASCO/ESMO guideline parser (168 lines)
│   │   ├── pathway_parser.py      # Signaling pathway parser (121 lines)
│   │   ├── resistance_parser.py   # Resistance mechanism parser (125 lines)
│   │   └── outcome_parser.py      # Outcome record parser (158 lines)
│   └── utils/
│       ├── vcf_parser.py          # VCF file parsing utilities (361 lines)
│       └── pubmed_client.py       # NCBI E-utilities HTTP client (296 lines)
├── app/
│   └── oncology_ui.py             # Streamlit MTB Workbench (703 lines)
├── api/
│   ├── main.py                    # FastAPI REST server (347 lines)
│   └── routes/
│       ├── meta_agent.py          # /api/ask, /api/deep-research (169 lines)
│       ├── cases.py               # /api/cases, /api/cases/{id}/mtb (234 lines)
│       ├── trials.py              # /api/trials/match (153 lines)
│       ├── reports.py             # /api/reports/{format} (236 lines)
│       └── events.py              # /api/events, cross-modal triggers (89 lines)
├── config/
│   └── settings.py                # Pydantic BaseSettings (109 lines)
├── scripts/
│   ├── setup_collections.py       # Create collections + seed data
│   ├── ingest_pubmed.py           # CLI: PubMed ingest
│   ├── ingest_clinical_trials.py  # CLI: ClinicalTrials.gov ingest
│   ├── ingest_civic.py            # CLI: CIViC variant ingest
│   └── validate_e2e.py            # End-to-end validation (5 tests)
├── tests/
│   └── conftest.py                # Test fixtures (214 lines)
├── requirements.txt
└── LICENSE                        # Apache 2.0
```

**39 Python files | ~12,301 lines of code | Apache 2.0**

---

## 4. Docker Compose Services

### Service Architecture

| Service | Image | Port | Role |
|---|---|---|---|
| `milvus-etcd` | quay.io/coreos/etcd:v3.5.5 | 2379 | Milvus metadata store |
| `milvus-minio` | minio/minio:v2023.03 | 9000, 9001 | Milvus object storage |
| `milvus-standalone` | milvusdb/milvus:v2.4 | 19530, 9091 | Vector database |
| `onco-streamlit` | Built from Dockerfile | 8526 | Streamlit MTB Workbench |
| `onco-api` | Built from Dockerfile | 8527 | FastAPI REST server |
| `onco-setup` | Built from Dockerfile | -- | One-shot collection setup + seed |

### Startup Sequence

```bash
# 1. Configure
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY

# 2. Launch all services
docker compose up -d

# 3. Watch setup completion
docker compose logs -f onco-setup

# 4. Verify
curl http://localhost:8527/health
# -> {"status": "healthy", "collections": {...}, "total_vectors": 1490}
```

### Dockerfile (Multi-Stage)

**Stage 1 -- Builder:**
- Base: `python:3.10-slim`
- System deps: `build-essential`, `gcc`, `g++`, `libxml2-dev`, `libxslt1-dev`
- Virtual environment with all pip dependencies

**Stage 2 -- Runtime:**
- Base: `python:3.10-slim`
- Copy venv from builder
- Copy source: `config/`, `src/`, `app/`, `api/`, `scripts/`
- Non-root user: `oncouser`
- Healthcheck: `curl -f http://localhost:8526/_stcore/health`

---

## 5. Milvus Collection Schemas

### 11 Collections Overview

| Collection | Records | Primary Fields | Source |
|---|---|---|---|
| `onco_variants` | ~300 | gene, variant_name, variant_type, evidence_level, drugs | CIViC, OncoKB |
| `onco_literature` | ~500 | PMID, title, text_chunk, year, cancer_type, gene | PubMed E-utilities |
| `onco_therapies` | ~120 | drug_name, category, targets, approved_indications, MoA | FDA-approved therapies |
| `onco_guidelines` | ~100 | org, cancer_type, version, year, key_recommendations | NCCN, ASCO, ESMO |
| `onco_trials` | ~200 | NCT ID, title, phase, status, biomarker_criteria | ClinicalTrials.gov API v2 |
| `onco_biomarkers` | ~80 | name, biomarker_type, testing_method, clinical_cutoff | TMB, MSI-H, PD-L1, HRD |
| `onco_resistance` | ~80 | primary_therapy, gene, mechanism, bypass_pathway | Resistance mechanisms |
| `onco_pathways` | ~50 | name, key_genes, therapeutic_targets, cross_talk | Signaling pathways |
| `onco_outcomes` | ~50 | therapy, cancer_type, response, duration_months | Treatment outcomes |
| `onco_cases` | ~10 | patient_id, cancer_type, stage, variants, biomarkers | Case snapshots |
| `genomic_evidence` | 3,561,170 | chrom, pos, gene, consequence, clinical_significance | Shared (read-only) |

### Collection Field Schemas

#### onco_literature

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key (PMID) |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `title` | VARCHAR | 500 | Article title |
| `text_chunk` | VARCHAR | 3000 | Chunked abstract text |
| `source_type` | VARCHAR | 20 | pubmed, pmc, preprint, manual |
| `year` | INT64 | -- | Publication year |
| `cancer_type` | VARCHAR | 50 | Classified cancer type |
| `gene` | VARCHAR | 50 | Primary gene discussed |
| `variant` | VARCHAR | 100 | Variant if mentioned |
| `keywords` | VARCHAR | 1000 | Extracted keywords |
| `journal` | VARCHAR | 200 | Journal name |

#### onco_trials

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(20) | 20 | Primary key (NCT ID) |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `title` | VARCHAR | 500 | Trial title |
| `text_summary` | VARCHAR | 3000 | Brief description |
| `phase` | VARCHAR | 30 | Phase 1, Phase 2, Phase 3 |
| `status` | VARCHAR | 30 | Recruiting, Completed, etc. |
| `sponsor` | VARCHAR | 200 | Lead sponsor |
| `cancer_types` | VARCHAR | 200 | Comma-separated cancer types |
| `biomarker_criteria` | VARCHAR | 500 | Biomarker eligibility criteria |
| `enrollment` | INT64 | -- | Target enrollment |
| `start_year` | INT64 | -- | Study start year |
| `outcome_summary` | VARCHAR | 2000 | Outcome results if available |

#### onco_variants

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `gene` | VARCHAR | 50 | Gene symbol |
| `variant_name` | VARCHAR | 100 | e.g. V600E, L858R |
| `variant_type` | VARCHAR | 30 | snv, indel, fusion, cnv_amplification |
| `cancer_type` | VARCHAR | 50 | Associated cancer type |
| `evidence_level` | VARCHAR | 20 | AMP/ASCO/CAP tier (A-E) |
| `drugs` | VARCHAR | 500 | Associated drugs (comma-separated) |
| `civic_id` | VARCHAR | 20 | CIViC variant ID |
| `vrs_id` | VARCHAR | 100 | GA4GH VRS identifier |
| `text_summary` | VARCHAR | 3000 | Clinical summary text |
| `clinical_significance` | VARCHAR | 200 | Pathogenic, likely pathogenic, etc. |
| `allele_frequency` | FLOAT | -- | Population allele frequency |

#### onco_biomarkers

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `name` | VARCHAR | 100 | Biomarker name |
| `biomarker_type` | VARCHAR | 30 | predictive, prognostic, diagnostic |
| `cancer_types` | VARCHAR | 200 | Associated cancer types |
| `predictive_value` | VARCHAR | 200 | Predictive interpretation |
| `testing_method` | VARCHAR | 100 | IHC, PCR, NGS, FISH |
| `clinical_cutoff` | VARCHAR | 100 | e.g. >= 10 mut/Mb |
| `text_summary` | VARCHAR | 3000 | Clinical summary |
| `evidence_level` | VARCHAR | 20 | Evidence tier |

#### onco_therapies

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `drug_name` | VARCHAR | 200 | Generic drug name |
| `category` | VARCHAR | 30 | targeted, immunotherapy, chemo |
| `targets` | VARCHAR | 200 | Molecular targets |
| `approved_indications` | VARCHAR | 500 | FDA-approved indications |
| `resistance_mechanisms` | VARCHAR | 500 | Known resistance patterns |
| `evidence_level` | VARCHAR | 20 | Evidence tier |
| `text_summary` | VARCHAR | 3000 | Drug summary |
| `mechanism_of_action` | VARCHAR | 500 | MoA description |

#### onco_pathways

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `name` | VARCHAR | 100 | Pathway name |
| `key_genes` | VARCHAR | 500 | Genes in pathway |
| `therapeutic_targets` | VARCHAR | 300 | Druggable nodes |
| `cross_talk` | VARCHAR | 500 | Cross-talk pathways |
| `text_summary` | VARCHAR | 3000 | Pathway summary |

#### onco_guidelines

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `org` | VARCHAR | 20 | NCCN, ASCO, ESMO, CAP/AMP |
| `cancer_type` | VARCHAR | 50 | Cancer type |
| `version` | VARCHAR | 20 | Guideline version |
| `year` | INT64 | -- | Publication year |
| `key_recommendations` | VARCHAR | 3000 | Recommendations |
| `text_summary` | VARCHAR | 3000 | Summary text |
| `evidence_level` | VARCHAR | 20 | Evidence tier |

#### onco_resistance

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `primary_therapy` | VARCHAR | 200 | Drug that resistance applies to |
| `gene` | VARCHAR | 50 | Gene involved |
| `mechanism` | VARCHAR | 500 | Resistance mechanism |
| `bypass_pathway` | VARCHAR | 200 | Bypass signaling pathway |
| `alternative_therapies` | VARCHAR | 500 | Alternative drugs |
| `text_summary` | VARCHAR | 3000 | Mechanism summary |

#### onco_outcomes

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `case_id` | VARCHAR | 100 | Linked case |
| `therapy` | VARCHAR | 200 | Drug used |
| `cancer_type` | VARCHAR | 50 | Cancer type |
| `response` | VARCHAR | 20 | CR, PR, SD, PD |
| `duration_months` | FLOAT | -- | Response duration |
| `toxicities` | VARCHAR | 500 | Observed toxicities |
| `biomarkers_at_baseline` | VARCHAR | 500 | Baseline biomarkers |
| `text_summary` | VARCHAR | 3000 | Outcome summary |

#### onco_cases

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(100) | 100 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `patient_id` | VARCHAR | 100 | De-identified patient ID |
| `cancer_type` | VARCHAR | 50 | Cancer type |
| `stage` | VARCHAR | 20 | AJCC stage |
| `variants` | VARCHAR | 1000 | Variant list (JSON) |
| `biomarkers` | VARCHAR | 1000 | Biomarker results (JSON) |
| `prior_therapies` | VARCHAR | 500 | Prior drug list |
| `text_summary` | VARCHAR | 3000 | Case summary |

#### genomic_evidence (read-only)

| Field | DataType | max_length | Notes |
|---|---|---|---|
| `id` | VARCHAR(200) | 200 | Primary key |
| `embedding` | FLOAT_VECTOR(384) | -- | BGE-small-en-v1.5 |
| `chrom` | VARCHAR | 10 | Chromosome |
| `pos` | INT64 | -- | Position |
| `ref` | VARCHAR | 500 | Reference allele |
| `alt` | VARCHAR | 500 | Alternate allele |
| `qual` | FLOAT | -- | Quality score |
| `gene` | VARCHAR | 50 | Gene symbol |
| `consequence` | VARCHAR | 100 | VEP consequence |
| `impact` | VARCHAR | 20 | HIGH, MODERATE, LOW |
| `genotype` | VARCHAR | 10 | e.g. 0/1, 1/1 |
| `text_summary` | VARCHAR | 2000 | Variant summary |
| `clinical_significance` | VARCHAR | 200 | ClinVar classification |
| `rsid` | VARCHAR | 20 | dbSNP ID |
| `disease_associations` | VARCHAR | 500 | Disease links |
| `am_pathogenicity` | FLOAT | -- | AlphaMissense score |
| `am_class` | VARCHAR | 30 | AlphaMissense class |

### Unified Index Configuration

```python
INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}

SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

### OncoCollectionManager

| Method | Purpose |
|---|---|
| `create_all_collections()` | Create all 11 collections with IVF_FLAT indexes |
| `get_collection_stats()` | Return entity count and field names per collection |
| `insert_batch()` | Bulk insert records (dict-of-lists transpose for pymilvus) |
| `search()` | Single-collection vector similarity search with optional filters |
| `search_all()` | **Parallel** multi-collection search via ThreadPoolExecutor (up to 6 threads) |

---

## 6. Pydantic Data Models

### Enumeration Types

```python
CancerType:       NSCLC, SCLC, BREAST, COLORECTAL, MELANOMA, PANCREATIC, OVARIAN,
                  PROSTATE, RENAL, BLADDER, HEAD_NECK, HEPATOCELLULAR, GASTRIC,
                  GLIOBLASTOMA, AML, CML, ALL, CLL, DLBCL, MULTIPLE_MYELOMA, OTHER

VariantType:      SNV, INDEL, CNV_AMP, CNV_DEL, FUSION, REARRANGEMENT, SV

EvidenceLevel:    A (FDA-approved), B (clinical evidence), C (case reports),
                  D (preclinical), E (computational)

TherapyCategory:  TARGETED, IMMUNOTHERAPY, CHEMOTHERAPY, HORMONAL, COMBINATION,
                  RADIOTHERAPY, CELL_THERAPY

TrialPhase:       EARLY_PHASE_1, PHASE_1, PHASE_1_2, PHASE_2, PHASE_2_3,
                  PHASE_3, PHASE_4, NA

TrialStatus:      NOT_YET_RECRUITING, RECRUITING, ENROLLING_BY_INVITATION,
                  ACTIVE_NOT_RECRUITING, SUSPENDED, TERMINATED, COMPLETED,
                  WITHDRAWN, UNKNOWN

ResponseCategory: CR, PR, SD, PD, NE

BiomarkerType:    PREDICTIVE, PROGNOSTIC, DIAGNOSTIC, MONITORING, RESISTANCE,
                  PHARMACODYNAMIC

PathwayName:      MAPK, PI3K_AKT_MTOR, DDR, CELL_CYCLE, APOPTOSIS,
                  WNT, NOTCH, HEDGEHOG, JAK_STAT, ANGIOGENESIS

GuidelineOrg:     NCCN, ESMO, ASCO, WHO, CAP_AMP

SourceType:       PUBMED, PMC, PREPRINT, MANUAL
```

### Collection Models

Each collection has a corresponding Pydantic model with a `to_embedding_text()` method that combines key fields for BGE encoding:

| Model | Collection | Key Fields |
|---|---|---|
| `OncologyLiterature` | onco_literature | PMID, title, text_chunk, source_type, year, cancer_type, gene, variant, keywords |
| `OncologyTrial` | onco_trials | NCT ID, title, text_summary, phase, status, cancer_types, biomarker_criteria, enrollment |
| `OncologyVariant` | onco_variants | gene, variant_name, variant_type, evidence_level, drugs, civic_id, clinical_significance |
| `OncologyBiomarker` | onco_biomarkers | name, biomarker_type, cancer_types, testing_method, clinical_cutoff, evidence_level |
| `OncologyTherapy` | onco_therapies | drug_name, category, targets, approved_indications, evidence_level, mechanism_of_action |
| `OncologyPathway` | onco_pathways | name, key_genes, therapeutic_targets, cross_talk |
| `OncologyGuideline` | onco_guidelines | org, cancer_type, version, year, key_recommendations, evidence_level |
| `ResistanceMechanism` | onco_resistance | primary_therapy, gene, mechanism, bypass_pathway, alternative_therapies |
| `OutcomeRecord` | onco_outcomes | case_id, therapy, cancer_type, response, duration_months, toxicities, biomarkers_at_baseline |
| `CaseSnapshot` | onco_cases | patient_id, cancer_type, stage, variants, biomarkers, prior_therapies |

### Search & Agent Models

```python
SearchHit:               # Single evidence item (collection, id, score, text, metadata)
CrossCollectionResult:   # Merged multi-collection results (query, hits, knowledge_context, search_time_ms)
ComparativeResult:       # Dual-entity comparison (entity_a, entity_b, evidence_a, evidence_b)
MTBPacket:               # Molecular Tumor Board packet (variant_table, evidence_table, therapy_ranking,
                         #   trial_matches, open_questions, citations, generated_at)
AgentQuery:              # User input (question, cancer_type, gene, include_genomic)
AgentResponse:           # Agent output (question, answer, evidence, knowledge_used, timestamp)
```

---

## 7. Configuration Reference

### OncoSettings (Pydantic BaseSettings)

```python
# Environment prefix: ONCO_
# Loaded from: .env file or environment variables

# Paths
PROJECT_ROOT: Path              # Auto-detected
DATA_DIR: Path                  # PROJECT_ROOT / "data"
CACHE_DIR: Path                 # PROJECT_ROOT / "cache"
REFERENCE_DIR: Path             # PROJECT_ROOT / "reference"
RAG_PIPELINE_ROOT: Path         # For shared resources

# Milvus
MILVUS_HOST: str = "localhost"
MILVUS_PORT: int = 19530

# Collection Names (11 total)
COLLECTION_LITERATURE: str = "onco_literature"
COLLECTION_TRIALS: str = "onco_trials"
COLLECTION_VARIANTS: str = "onco_variants"
COLLECTION_BIOMARKERS: str = "onco_biomarkers"
COLLECTION_THERAPIES: str = "onco_therapies"
COLLECTION_PATHWAYS: str = "onco_pathways"
COLLECTION_GUIDELINES: str = "onco_guidelines"
COLLECTION_RESISTANCE: str = "onco_resistance"
COLLECTION_OUTCOMES: str = "onco_outcomes"
COLLECTION_CASES: str = "onco_cases"
COLLECTION_GENOMIC: str = "genomic_evidence"

# Embeddings
EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
EMBEDDING_DIM: int = 384
EMBEDDING_BATCH_SIZE: int = 32

# LLM
LLM_PROVIDER: str = "anthropic"
LLM_MODEL: str = "claude-sonnet-4-20250514"
ANTHROPIC_API_KEY: Optional[str] = None

# RAG Search
TOP_K: int = 5
SCORE_THRESHOLD: float = 0.4

# Collection Weights (sum ~1.0)
WEIGHT_VARIANTS: float = 0.18
WEIGHT_LITERATURE: float = 0.16
WEIGHT_THERAPIES: float = 0.14
WEIGHT_GUIDELINES: float = 0.12
WEIGHT_TRIALS: float = 0.10
WEIGHT_BIOMARKERS: float = 0.08
WEIGHT_RESISTANCE: float = 0.07
WEIGHT_PATHWAYS: float = 0.06
WEIGHT_OUTCOMES: float = 0.04
WEIGHT_CASES: float = 0.02
WEIGHT_GENOMIC: float = 0.03

# External APIs
NCBI_API_KEY: Optional[str] = None
PUBMED_MAX_RESULTS: int = 5000
CT_GOV_BASE_URL: str = "https://clinicaltrials.gov/api/v2"
CIVIC_BASE_URL: str = "https://civicdb.org/api"

# Services
API_HOST: str = "0.0.0.0"
API_PORT: int = 8527
STREAMLIT_PORT: int = 8526

# Metrics & Scheduling
METRICS_ENABLED: bool = True
SCHEDULER_INTERVAL: str = "168h"        # Weekly refresh

# Conversation Memory
CONVERSATION_MEMORY_DEPTH: int = 3      # Prior exchanges

# Citation Thresholds
CITATION_STRONG_THRESHOLD: float = 0.75
CITATION_MODERATE_THRESHOLD: float = 0.60

# Cross-Modal
CROSS_MODAL_ENABLED: bool = True
```

---

## 8. Embedding Strategy

### BGE-small-en-v1.5

| Parameter | Value |
|---|---|
| Model | BAAI/bge-small-en-v1.5 |
| Dimensions | 384 |
| Parameters | ~33M |
| Encoding | Asymmetric (query vs document) |

### Asymmetric Encoding

```python
# Query embedding (with instruction prefix)
_BGE_INSTRUCTION = "Represent this sentence for searching relevant passages: "
query_embedding = embedder.encode(_BGE_INSTRUCTION + question)

# Document embedding (raw text, no prefix)
doc_embedding = embedder.encode(record.to_embedding_text())
```

This asymmetric approach improves retrieval relevance by 5-15% compared to symmetric encoding.

### to_embedding_text() Pattern

Each Pydantic model generates its embedding input by combining key fields:

```python
# OncologyVariant
def to_embedding_text(self) -> str:
    parts = [
        f"{self.gene} {self.variant_name}",
        self.text_summary,
        f"Type: {self.variant_type.value}",
        f"Evidence: {self.evidence_level.value}",
    ]
    if self.cancer_type:
        parts.append(f"Cancer: {self.cancer_type.value}")
    if self.drugs:
        parts.append(f"Drugs: {', '.join(self.drugs)}")
    return " | ".join(parts)
```

---

## 9. Multi-Collection RAG Engine

### System Prompt

```python
ONCO_SYSTEM_PROMPT = """You are a Precision Oncology Intelligence Agent — an expert AI assistant
purpose-built for clinical and translational oncology decision support.

Your core competencies include:
* Molecular profiling — somatic/germline variant interpretation, TMB, MSI, CNV, fusions
* Variant interpretation — CIViC/OncoKB evidence levels, AMP/ASCO/CAP classification
* Therapy selection — NCCN/ESMO guideline-concordant treatment recommendations
* Clinical trial matching — eligibility assessment against active registrations
* Resistance mechanisms — on-target mutations, bypass signaling, lineage plasticity
* Biomarker assessment — TMB, MSI, PD-L1, HRD scoring; companion diagnostics
* Outcomes monitoring — RECIST response criteria, survival endpoints, ctDNA dynamics
* Cross-modal integration — linking genomic findings to imaging, pathology, drug discovery

Behavioral instructions:
1. Cite evidence with clickable PubMed or ClinicalTrials.gov links
2. Think cross-functionally — connect variants to therapies, trials, resistance
3. Highlight resistance & contraindications proactively
4. Reference guidelines — cite NCCN, ESMO, or ASCO when recommending treatments
5. Acknowledge uncertainty — state evidence gaps and MTB review needs"""
```

### Collection Configuration

```python
COLLECTION_CONFIG = {
    "onco_variants":    {"weight": 0.18, "label": "Variant",    "filter_field": "gene"},
    "onco_literature":  {"weight": 0.16, "label": "Literature",  "filter_field": "gene",   "year_field": "year"},
    "onco_therapies":   {"weight": 0.14, "label": "Therapy",     "filter_field": None},
    "onco_guidelines":  {"weight": 0.12, "label": "Guideline",   "filter_field": None,     "year_field": "year"},
    "onco_trials":      {"weight": 0.10, "label": "Trial",       "filter_field": None,     "year_field": "start_year"},
    "onco_biomarkers":  {"weight": 0.08, "label": "Biomarker",   "filter_field": None},
    "onco_resistance":  {"weight": 0.07, "label": "Resistance",  "filter_field": "gene"},
    "onco_pathways":    {"weight": 0.06, "label": "Pathway",     "filter_field": None},
    "onco_outcomes":    {"weight": 0.04, "label": "Outcome",     "filter_field": None},
    "onco_cases":       {"weight": 0.02, "label": "Case",        "filter_field": None},
    "genomic_evidence": {"weight": 0.03, "label": "Genomic",     "filter_field": None},
}
```

### Core Methods

| Method | Purpose | Returns |
|---|---|---|
| `retrieve()` | Parallel search + expansion + knowledge augmentation | `CrossCollectionResult` |
| `query()` | Full RAG: retrieve + Claude generation | `str` (answer) |
| `query_stream()` | Streaming RAG: yields evidence then tokens | `Generator[str]` |
| `find_related()` | Cross-collection entity linking | `Dict[str, List[SearchHit]]` |
| `retrieve_comparative()` | Dual-entity comparative retrieval | `Dict` (entity_a/b, hits, shared) |

### Retrieve Pipeline

```
User Query
    |
1. Embed query (BGE asymmetric prefix)                        [< 5 ms]
    |
2. Parallel search across 11 collections (top-5 each)        [< 200 ms]
   (ThreadPoolExecutor, up to 8 concurrent threads)
    |
3. Query expansion: detect keywords -> expand -> re-search    [8-12 ms]
    |
4. Merge + deduplicate + weighted rank (cap at 30)            [< 1 ms]
    |
5. Score citations: high (>=0.75) / medium (>=0.60) / low     [< 1 ms]
    |
6. Knowledge graph augmentation (gene, therapy, resistance)   [< 1 ms]
    |
7. Build prompt: evidence + knowledge + question              [< 1 ms]
    |
8. Stream Claude Sonnet 4.6 response                          [~22-24 sec]
```

### Citation Formatting

```python
# Literature -> PubMed link
"[PubMed PMID](https://pubmed.ncbi.nlm.nih.gov/PMID/)"

# Trial -> ClinicalTrials.gov link
"[NCT...](https://clinicaltrials.gov/study/NCT...)"

# Other collections -> bracketed reference
"[Variant: EGFR_V600E]"
```

---

## 10. Knowledge Graph

### Components

| Component | Count | Description |
|---|---|---|
| **Actionable Targets** | ~40 | Full profiles with evidence levels, drugs, resistance, and pathways |
| **Therapy Mappings** | ~30 | Drug names, brand names, categories, MoA, indications |
| **Resistance Mechanisms** | ~20 | On-target mutations, bypass pathways, alternative therapies |
| **Signaling Pathways** | ~10 | Key genes, druggable nodes, cross-talk relationships |
| **Biomarker Panels** | ~15 | TMB, MSI-H, PD-L1, HRD, NTRK fusion, ctDNA, FGFR |
| **Entity Aliases** | ~50+ | Cancer type aliases, drug brand names, gene synonyms |

### Actionable Targets (~40 entries)

```
Lung:         EGFR, ALK, ROS1, KRAS (G12C), RET, MET, BRAF, HER2, NTRK
Breast:       HER2, ESR1, PIK3CA, BRCA1, BRCA2
Colorectal:   KRAS, NRAS, BRAF, MSI-H, HER2
Melanoma:     BRAF (V600E/K), NRAS, KIT, NTRK
Ovarian:      BRCA1, BRCA2, HRD, NTRK
Pan-tumor:    TMB-H, MSI-H/dMMR, NTRK fusion, RET fusion
Other:        IDH1, IDH2, FGFR, PTEN, TP53, CDK4/6, PDGFRA
```

Each target includes: evidence level, known variants (with per-variant evidence tiers), associated drugs, gene-level actionability flag, and default evidence level.

### Therapy Mappings

| Drug | Brand | Category | Targets | Evidence Level |
|---|---|---|---|---|
| vemurafenib | Zelboraf | Targeted | BRAF V600E | A |
| osimertinib | Tagrisso | Targeted | EGFR T790M, exon19/21 | A |
| pembrolizumab | Keytruda | Immunotherapy | PD-1, MSI-H, TMB-H | A |
| sotorasib | Lumakras | Targeted | KRAS G12C | A |
| lorlatinib | Lorbrena | Targeted | ALK | A |
| olaparib | Lynparza | PARP inhibitor | BRCA1/2, HRD | A |
| trastuzumab deruxtecan | Enhertu | ADC | HER2 | A |
| larotrectinib | Vitrakvi | Targeted | NTRK | A |
| entrectinib | Rozlytrek | Targeted | NTRK, ROS1 | A |

### Knowledge Graph API

```python
get_target_context("EGFR")            # Full EGFR knowledge block
get_therapy_context("osimertinib")     # Drug profile with MoA + indications
get_resistance_context("EGFR T790M")   # Resistance mechanism + alternatives
get_pathway_context("MAPK")            # Pathway genes, targets, cross-talk
get_biomarker_context("TMB-H")         # Biomarker threshold + drugs
resolve_comparison_entity("Keytruda")  # -> {"type": "drug", "canonical": "pembrolizumab"}
get_knowledge_stats()                  # {targets: ~40, therapies: ~30, ...}
```

---

## 11. Query Expansion

### 12 Expansion Map Categories

| Category | Keywords | Expanded Terms | Examples |
|---|---|---|---|
| Cancer Type | 20+ | ~120 | nsclc -> [non-small cell lung cancer, adenocarcinoma, squamous, ...] |
| Gene | 20+ | ~100 | egfr -> [epidermal growth factor receptor, ERBB1, HER1, ...] |
| Therapy | 15+ | ~80 | immunotherapy -> [checkpoint inhibitor, PD-1, PD-L1, CTLA-4, ...] |
| Biomarker | 12+ | ~60 | tmb -> [tumor mutational burden, mutations per megabase, ...] |
| Pathway | 10+ | ~50 | mapk -> [RAS-RAF-MEK-ERK, BRAF, KRAS, MEK inhibitor, ...] |
| Resistance | 10+ | ~40 | resistance -> [acquired resistance, bypass pathway, reversion, ...] |
| Clinical | 10+ | ~50 | stage -> [staging, AJCC, TNM, metastatic, locally advanced, ...] |
| Trial | 8+ | ~40 | clinical trial -> [phase, recruiting, eligibility, randomized, ...] |
| Immunotherapy | 10+ | ~50 | immune -> [T-cell, checkpoint, CAR-T, TIL, neoantigen, ...] |
| Surgery/Radiation | 8+ | ~30 | surgery -> [resection, lobectomy, SBRT, IMRT, proton, ...] |
| Toxicity | 8+ | ~30 | toxicity -> [adverse event, grade 3, dose limiting, ...] |
| Genomics | 8+ | ~50 | variant -> [mutation, SNV, indel, fusion, amplification, ...] |
| **Total** | **~120** | **~700** | |

### How It Works

```python
def expand_query(question: str) -> List[str]:
    """Detect keywords in user question, return expansion terms.

    Used by RAG engine to broaden recall via semantic re-search.
    Returns up to 10 related terms per matched keyword.
    """
```

---

## 12. Comparative Analysis Mode

### Auto-Detection

Comparative queries are detected by regex matching: `compare`, `vs`, `versus`, `difference between`, `head-to-head`.

### Pipeline

```
User: "Compare osimertinib vs erlotinib for EGFR-mutant NSCLC"
    |
1. _is_comparative() detects "vs"                              [< 1 ms]
    |
2. _parse_comparison_entities() extracts "osimertinib", "erlotinib"  [< 1 ms]
    |
3. resolve_comparison_entity() for each entity                  [< 1 ms]
   "osimertinib" -> {"type": "drug", "canonical": "osimertinib"}
   "erlotinib"   -> {"type": "drug", "canonical": "erlotinib"}
    |
4. Dual retrieve() -- one per entity                            [~400 ms]
   Entity A: retrieve(question, gene=entity_a.gene)
   Entity B: retrieve(question, gene=entity_b.gene)
    |
5. Identify shared / head-to-head evidence                      [< 1 ms]
    |
6. _build_comparative_prompt() -- structured comparison         [< 1 ms]
   8-point comparison template: MoA, efficacy, safety, biomarkers,
   resistance, guidelines, trial evidence, summary recommendation
    |
7. Stream Claude Sonnet 4.6 (max_tokens=3000)                  [~28-30 sec]
   Output: comparison table, head-to-head data, clinical context
```

### Supported Entity Types

| Type | Examples | Resolution |
|---|---|---|
| Genes / Targets | EGFR, BRAF, ALK, KRAS | Direct match in ACTIONABLE_TARGETS |
| Drugs | Osimertinib, Pembrolizumab | THERAPY_MAP -> canonical + targets |
| Drug Classes | PARP inhibitors, TKIs | ENTITY_ALIASES -> type: drug_class |
| Biomarkers | MSI-H, TMB-H, PD-L1 | BIOMARKER_PANELS -> type: biomarker |
| Cancer Types | NSCLC, Melanoma | CANCER_TYPE_ALIASES -> canonical |

### Fallback

If entity parsing fails, the query gracefully falls back to the standard single-query `retrieve()` path.

---

## 13. Case Manager Workflow

### OncologyCaseManager

The case manager orchestrates the complete MTB packet generation workflow from raw patient data to structured clinical decision support.

### Case Creation

```python
case = case_manager.create_case(
    patient_id="ONCO-2026-001",
    cancer_type="NSCLC",
    stage="IV",
    vcf_content_or_variants=vcf_text,     # Raw VCF text or pre-parsed variant list
    biomarkers={"MSI": "MSS", "TMB": 14.2, "PD-L1_TPS": 80},
    prior_therapies=["carboplatin", "pemetrexed"],
)
```

### VCF Parsing

```
VCF Text Input
    |
1. Skip header lines (starting with #)
    |
2. Split each data line into fixed fields (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)
    |
3. Filter: keep only PASS variants (FILTER = "PASS", ".", or "pass")
    |
4. Extract GENE from INFO field (regex: ANN|CSQ|GENE pattern)
    |
5. Extract CONSEQUENCE from INFO field (ANN|CSQ annotation)
    |
6. Build human-readable variant string: "EGFR chr7:55259515 T>G"
    |
7. Classify actionability against ACTIONABLE_TARGETS for each variant
```

### MTB Packet Generation

```python
packet = case_manager.generate_mtb_packet(case_id)
```

The MTB packet contains 6 sections:

| Section | Content | Source |
|---|---|---|
| **Variant Table** | All variants with actionability, evidence level, associated drugs | ACTIONABLE_TARGETS lookup |
| **Evidence Table** | RAG-retrieved evidence per actionable variant | onco_literature + onco_therapies retrieval |
| **Therapy Ranking** | Ranked therapies by evidence level, resistance check | ACTIONABLE_TARGETS + prior therapy history |
| **Trial Matches** | Matching clinical trials for cancer type + biomarkers | onco_trials semantic search |
| **Open Questions** | VUS variants, missing biomarkers, evidence gaps | Automated gap analysis |
| **Citations** | PubMed and ClinicalTrials.gov links | RAG retrieval metadata |

### Actionability Classification (AMP/ASCO/CAP)

| Level | Description | Action |
|---|---|---|
| **A** | Validated, FDA-approved companion diagnostic | Strong recommendation |
| **B** | Clinical evidence from well-powered studies | Consider in treatment plan |
| **C** | Case reports and small series | Discuss at MTB |
| **D** | Preclinical / in-vitro data | Research context only |
| **E** | Inferential / computational prediction | Lowest evidence tier |
| **VUS** | Variant of uncertain significance | Flag for further testing |

---

## 14. Trial Matcher Implementation

### TrialMatcher

The trial matcher uses a 4-step hybrid approach combining deterministic and semantic matching.

### Matching Algorithm

```
Patient Profile (cancer_type, biomarkers, stage, age)
    |
Step 1: Deterministic Filter                                [~50 ms]
    Filter on cancer_type + status="Recruiting"
    Returns up to 3x top_k candidates
    |
Step 2: Semantic Search                                     [~100 ms]
    Embed: "{cancer_type} clinical trial stage {stage} {biomarkers}"
    Vector search against onco_trials collection
    Returns up to 3x top_k candidates
    |
Step 3: Merge + Composite Scoring                           [< 5 ms]
    Union by trial_id (keep best score from both sources)
    Composite score = 0.40 * biomarker_match
                    + 0.25 * semantic_score
                    + 0.20 * phase_weight
                    + 0.15 * status_weight
    |
Step 4: Explanation Generation                              [< 5 ms]
    For each trial:
      - List matched criteria (cancer type, biomarkers, stage)
      - List unmatched/unconfirmed criteria
      - Phase and status context
```

### Phase Weights

| Phase | Weight |
|---|---|
| Phase 3 | 1.0 |
| Phase 2/3 | 0.9 |
| Phase 2 | 0.8 |
| Phase 1/2 | 0.7 |
| Phase 1 | 0.6 |
| Phase 4 | 0.5 |

### Status Weights

| Status | Weight |
|---|---|
| Recruiting | 1.0 |
| Enrolling by invitation | 0.8 |
| Active, not recruiting | 0.6 |
| Not yet recruiting | 0.4 |

### Biomarker Matching

Fuzzy text matching: for each patient biomarker key/value pair, check if it appears (case-insensitive) in the trial criteria text. Returns fraction matched (0.0 to 1.0).

---

## 15. Therapy Ranker Rules

### TherapyRanker

The therapy ranker implements a 6-step evidence-based ranking pipeline.

### Ranking Algorithm

```
Patient Profile (cancer_type, variants, biomarkers, prior_therapies)
    |
Step 1: Variant-Driven Therapy Identification               [< 5 ms]
    For each variant: check ACTIONABLE_TARGETS for matching drugs
    Extract evidence level from variant-specific or gene-level data
    |
Step 2: Biomarker-Driven Therapy Identification              [< 5 ms]
    Standard biomarker-therapy mappings:
      MSI-H / dMMR    -> pembrolizumab (Level A)
      TMB-H >= 10     -> pembrolizumab (Level A)
      HRD / BRCA      -> olaparib, rucaparib, niraparib (Level A/B)
      PD-L1 TPS >= 50 -> pembrolizumab first-line (Level A)
      NTRK fusion     -> larotrectinib, entrectinib (Level A)
    Also checks BIOMARKER_PANELS for additional mappings
    |
Step 3: Deduplicate + Evidence-Level Sort                    [< 1 ms]
    Combine variant + biomarker therapies
    Deduplicate by drug_name (keep stronger evidence)
    Sort: A > B > C > D > E
    |
Step 4: Resistance Check                                     [< 1 ms]
    For each candidate: check RESISTANCE_MAP against prior_therapies
    Flag if prior therapy triggers known resistance mechanism
    Include: resistance mechanism detail, alternative drugs
    |
Step 5: Contraindication Check                               [< 1 ms]
    Same drug previously used -> flag
    Same drug class (via THERAPY_MAP) as prior failed therapy -> flag
    |
Step 6: Supporting Evidence Retrieval                        [~200 ms]
    For each therapy: search onco_therapies + onco_literature
    Attach top-3 evidence items per drug with citations
    |
Final Ranking:
    Clean therapies first (sorted by evidence), then flagged (resistance/contraindication)
    Assign rank 1..N
```

### Evidence Level Ordering

```python
EVIDENCE_LEVEL_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "VUS": 5}
```

### Output Schema

Each ranked therapy contains:

```python
{
    "rank": 1,
    "drug_name": "osimertinib",
    "brand_name": "Tagrisso",
    "category": "targeted therapy",
    "targets": ["EGFR"],
    "evidence_level": "A",
    "guideline_recommendation": "NCCN first-line for EGFR-mutant NSCLC",
    "source": "variant",                    # variant | biomarker | biomarker_panel
    "source_gene": "EGFR",
    "resistance_flag": False,
    "resistance_detail": None,
    "contraindication_flag": False,
    "supporting_evidence": [
        {"collection": "onco_therapies", "source": "...", "text": "...", "score": 0.89}
    ],
}
```

---

## 16. Cross-Modal Integration

### OncoCrossModalTrigger

Connects oncology variants to imaging intelligence and drug discovery pipelines. When Level A or B actionable variants are detected, it queries:

1. **genomic_evidence** collection for variant context and literature
2. **imaging_*** collections for relevant imaging correlates (graceful failure if unavailable)

### Trigger Logic

```
Case Variants
    |
1. Filter: evidence level A or B only
    |
2. If no actionable A/B variants -> do not fire
    |
3. Build genomic queries:
   "{gene} {variant} targeted therapy evidence"
   "{gene} mutation clinical significance"
    |
4. Build imaging queries (if cancer_type available):
   "{gene} mutation {cancer_type} imaging findings"
    |
5. Query genomic_evidence (top-5 per query, threshold >= 0.40)
    |
6. Discover imaging_* collections and query (graceful failure)
    |
7. Build enrichment summary
```

### CrossModalResult

```python
@dataclass
class CrossModalResult:
    trigger_reason: str              # Why the trigger fired
    actionable_variants: List[Dict]  # Level A/B variants
    genomic_context: List[Dict]      # Genomic evidence hits
    imaging_context: List[Dict]      # Imaging findings (if available)
    genomic_hit_count: int
    imaging_hit_count: int
    enrichment_summary: str          # Human-readable multi-line summary
```

---

## 17. Data Ingest Pipelines

### BaseIngestPipeline Pattern

```python
class BaseIngestPipeline(ABC):
    def fetch(self, **kwargs) -> Any:          # Retrieve raw data
    def parse(self, raw_data) -> List[Model]:  # Validate into Pydantic models
    def embed_and_store(self, records) -> int:  # Embed + insert to Milvus
    def run(self, **kwargs) -> int:            # Orchestrate full pipeline
```

### Ingest Pipelines

| Pipeline | Source | Records | Time | Collection |
|---|---|---|---|---|
| PubMed | NCBI E-utilities (esearch + efetch) | ~500 | ~15 min | onco_literature |
| ClinicalTrials.gov | REST API v2 | ~200 | ~3 min | onco_trials |
| CIViC Variants | CIViC database API | ~300 | ~2 min | onco_variants |
| Seed Data | setup_collections.py --seed | ~290 | ~30 sec | therapies, guidelines, biomarkers, resistance, pathways, outcomes, cases |

### CIViC Variant Ingest

```bash
python3 scripts/ingest_civic.py
```

Fetches clinically actionable variants from CIViC, maps evidence levels to AMP/ASCO/CAP tiers (A-E), and stores in `onco_variants`.

### PubMed Literature Ingest

```bash
python3 scripts/ingest_pubmed.py --max-results 5000
```

Fetches oncology abstracts via NCBI E-utilities, classifies by cancer type and gene, embeds with BGE-small, and stores in `onco_literature`.

### ClinicalTrials.gov Ingest

```bash
python3 scripts/ingest_clinical_trials.py --max-results 1500
```

Fetches oncology trials via ClinicalTrials.gov API v2, extracts biomarker criteria, embeds, and stores in `onco_trials`.

### Seed Data

```bash
python3 scripts/setup_collections.py --seed
```

Seeds: `onco_therapies` (~120), `onco_guidelines` (~100), `onco_biomarkers` (~80), `onco_resistance` (~80), `onco_pathways` (~50), `onco_outcomes` (~50), `onco_cases` (~10).

---

## 18. FastAPI REST Server

### Service Configuration

- **Port:** 8527
- **Command:** `uvicorn api.main:app --host 0.0.0.0 --port 8527`

### Core Endpoints (api/main.py)

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/health` | Service health + collection counts + total vectors + service status |
| `GET` | `/collections` | List all collections with entity counts |
| `POST` | `/query` | Full RAG query (retrieve + Claude generation) |
| `POST` | `/search` | Evidence-only vector search (no LLM) |
| `POST` | `/find-related` | Cross-collection entity linking |
| `GET` | `/knowledge/stats` | Knowledge graph statistics |
| `GET` | `/metrics` | Prometheus-compatible metrics |

### Route Module Endpoints

| Method | Path | Module | Purpose |
|---|---|---|---|
| `POST` | `/api/ask` | meta_agent | Unified clinical Q&A (agent + RAG fallback) |
| `POST` | `/api/cases` | cases | Create oncology case from variants/VCF |
| `GET` | `/api/cases/{case_id}` | cases | Retrieve case details |
| `POST` | `/api/cases/{case_id}/mtb` | cases | Generate MTB packet |
| `GET` | `/api/cases/{case_id}/variants` | cases | List case variants |
| `POST` | `/api/trials/match` | trials | Match trials to patient profile |
| `POST` | `/api/trials/match-case/{case_id}` | trials | Match trials for existing case |
| `POST` | `/api/therapies/rank` | trials | Rank therapies by molecular profile |
| `POST` | `/api/reports/generate` | reports | Generate report (markdown/json/pdf) |
| `GET` | `/api/reports/{case_id}/{fmt}` | reports | Export case report (markdown/json/pdf/fhir) |

### Query Request Schema

```json
{
    "question": "What therapies target BRAF V600E in melanoma?",
    "cancer_type": "melanoma",
    "gene": "BRAF",
    "top_k": 10
}
```

### Ask Response Schema

```json
{
    "answer": "BRAF V600E melanoma has several FDA-approved targeted therapy options...",
    "sources": [
        {
            "collection": "onco_therapies",
            "text": "Vemurafenib is a selective BRAF V600E inhibitor...",
            "score": 0.91,
            "metadata": {"relevance": "high"}
        }
    ],
    "confidence": 0.89,
    "follow_up_questions": [
        "What resistance mechanisms emerge after BRAF inhibitor therapy?",
        "Should BRAF+MEK combination be used over BRAF monotherapy?",
        "What is the role of immunotherapy in BRAF-mutant melanoma?"
    ],
    "processing_time_ms": 24100
}
```

### Create Case Request Schema

```json
{
    "patient_id": "ONCO-2026-001",
    "cancer_type": "NSCLC",
    "stage": "IV",
    "variants": [
        {"gene": "EGFR", "variant": "L858R", "variant_type": "SNV"},
        {"gene": "TP53", "variant": "R273H", "variant_type": "SNV"}
    ],
    "biomarkers": {"MSI": "MSS", "TMB": 14.2, "PD-L1_TPS": 80},
    "prior_therapies": ["carboplatin", "pemetrexed"]
}
```

---

## 19. Streamlit MTB Workbench

### Service Configuration

- **Port:** 8526
- **Command:** `streamlit run app/oncology_ui.py --server.port 8526`
- **Theme:** NVIDIA dark (black/green)
- **API Backend:** `http://localhost:8527`

### 5-Tab Interface

| Tab | Description |
|---|---|
| **Case Management** | Create cases from VCF upload or variant entry; select cancer type (20+ types), AJCC stage, biomarkers (MSI, TMB, PD-L1, HRD), prior therapies |
| **Evidence Explorer** | RAG knowledge queries with streaming answers; Deep Research mode toggle; clickable citations; evidence panel grouped by collection |
| **Trial Matching** | Match patient profile to recruiting clinical trials; display ranked matches with match score, phase, and explanation |
| **Therapy Ranking** | Evidence-based therapy ranking with resistance flags; biomarker-driven recommendations; guideline concordance notes |
| **Outcomes Dashboard** | View treatment outcomes; comparative analysis; export reports (Markdown, JSON, PDF, FHIR R4) |

### Sidebar

- Collection statistics (10 owned + 1 read-only genomic)
- Deep Research mode toggle (activates agent reasoning pipeline)
- Cancer type filter dropdown (21 cancer types)
- Gene filter input
- Year range slider
- Export format selector

### Main Panel Features

- Chat message history (conversation memory: last 3 exchanges)
- Streaming response display
- Evidence panel (grouped by collection, color-coded relevance badges)
- Clickable PubMed and ClinicalTrials.gov links
- Comparative mode: auto-detected, produces structured tables

---

## 20. Export Formats

### Supported Formats

| Format | Method | Output |
|---|---|---|
| **Markdown** | `export_markdown()` | Structured report with headers, variant table, evidence, citations |
| **JSON** | `export_json()` | Machine-readable export with all metadata |
| **PDF** | `export_pdf()` | NVIDIA-themed report via ReportLab |
| **FHIR R4** | `export_fhir_r4()` | DiagnosticReport Bundle with SNOMED CT and LOINC coding |

### Report Structure

```
1. Query Metadata (timestamp, collections searched, evidence count)
2. Patient Summary (cancer type, stage, key biomarkers)
3. Variant Table (gene, variant, consequence, evidence level, drugs)
4. Evidence Sources (by collection, top 5 per collection)
5. Therapy Ranking (ranked therapies with resistance/contraindication flags)
6. Trial Matches (top matching trials with composite scores)
7. Open Questions (VUS, missing biomarkers, evidence gaps)
8. Citation Links (PubMed, ClinicalTrials.gov)
```

### PDF Export

- NVIDIA-themed styling (green accents on dark background)
- ReportLab-generated with professional layout
- Tables for variant and therapy data
- Clickable citation links embedded in PDF

### FHIR R4 Export

```json
{
    "resourceType": "Bundle",
    "type": "document",
    "entry": [
        {
            "resource": {
                "resourceType": "DiagnosticReport",
                "status": "final",
                "code": {
                    "coding": [
                        {
                            "system": "http://loinc.org",
                            "code": "51969-4",
                            "display": "Genetic analysis report"
                        }
                    ]
                },
                "conclusion": "..."
            }
        }
    ]
}
```

Coded with SNOMED CT (condition codes) and LOINC (observation codes) for interoperability.

---

## 21. Monitoring & Metrics

### Prometheus Metrics

```
onco_agent_up                           # Service availability (gauge)
onco_collection_vectors{collection}     # Vectors per collection (gauge)
```

### Endpoints

- **Metrics:** `GET /metrics` (Prometheus scrape format, text/plain)
- **Health:** `GET /health` (JSON with collection counts and service status)

---

## 22. Testing Strategy

### End-to-End Validation (`validate_e2e.py`)

5 tests:
1. Collection stats -- verify all 11 collections have expected record counts
2. Single-collection search -- verify vector similarity works
3. Multi-collection `search_all()` -- verify parallel search returns results from multiple collections
4. Filtered search -- verify gene filter works
5. Demo queries -- run all validated demo queries

### Demo Queries

```
1. "What therapies target BRAF V600E in melanoma?"
2. "Compare EGFR TKI generations for NSCLC"
3. "Resistance mechanisms to osimertinib"
4. "NCCN recommendations for HER2+ breast cancer"
5. "Match clinical trials for MSI-H colorectal cancer"
6. "What is the role of TMB as a predictive biomarker for immunotherapy?"
```

### conftest.py Fixtures

Test fixtures provide mock collection managers, embedders, knowledge stores, and sample data for unit testing without Milvus or LLM dependencies.

---

## 23. HCLS AI Factory Integration

### Shared Infrastructure

| Resource | Sharing |
|---|---|
| **Milvus 2.4** | Shared instance -- Oncology adds 10 owned collections alongside `genomic_evidence` |
| **BGE-small-en-v1.5** | Same embedding model as RAG/Chat pipeline |
| **ANTHROPIC_API_KEY** | Loaded from `rag-chat-pipeline/.env` if not set |
| **DGX Spark** | All services on single GB10 hardware ($3,999) |

### Architectural Insight

The platform is not disease-specific. By changing the knowledge graph, query expansion maps, and collection schemas, the same RAG architecture serves any therapeutic area. The Precision Oncology Agent demonstrates this with a clinical decision support domain (genomics-to-treatment) running on the same infrastructure as CAR-T cell therapy intelligence and imaging analysis.

---

## 24. Implementation Sequence

### Quick Start

```bash
# 1. Navigate
cd ai_agent_adds/precision_oncology_agent

# 2. Install
pip install -r requirements.txt

# 3. Configure
export ANTHROPIC_API_KEY=sk-ant-...

# 4. Create collections + seed data
python3 scripts/setup_collections.py --seed

# 5. Ingest PubMed (~15 min)
python3 scripts/ingest_pubmed.py --max-results 5000

# 6. Ingest ClinicalTrials.gov (~3 min)
python3 scripts/ingest_clinical_trials.py --max-results 1500

# 7. Ingest CIViC variants (~2 min)
python3 scripts/ingest_civic.py

# 8. Validate
python3 scripts/validate_e2e.py

# 9. Launch UI
streamlit run app/oncology_ui.py --server.port 8526

# 10. Launch API (separate terminal)
uvicorn api.main:app --host 0.0.0.0 --port 8527
```

### Docker Quick Start

```bash
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY
docker compose up -d
# UI: http://localhost:8526
# API: http://localhost:8527/docs
```

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| pydantic | >=2.0 | Data validation |
| pydantic-settings | >=2.0 | Environment config |
| pymilvus | >=2.4.0 | Vector DB client |
| sentence-transformers | >=2.2.0 | BGE embeddings |
| anthropic | >=0.18.0 | Claude API |
| streamlit | >=1.30.0 | Web UI |
| fastapi | >=0.109.0 | REST API |
| uvicorn[standard] | >=0.27.0 | ASGI server |
| lxml | >=5.0.0 | PubMed XML parsing |
| biopython | >=1.83 | NCBI E-utilities |
| reportlab | >=4.0.0 | PDF generation |
| prometheus-client | >=0.20.0 | Metrics export |
| apscheduler | >=3.10.0 | Job scheduling |

---

## Credits

- **Adam Jones** -- Architecture, implementation, knowledge curation
- **Apache 2.0 License**

---

*End of Project Bible*
