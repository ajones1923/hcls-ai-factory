# Pharmacogenomics Intelligence Agent -- Architecture Guide

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [Overview](#1-overview)
2. [Three-Tier Architecture](#2-three-tier-architecture)
3. [Data Flow](#3-data-flow)
4. [Ingest Tier](#4-ingest-tier)
5. [Vector Store Tier](#5-vector-store-tier)
6. [Inference Tier](#6-inference-tier)
7. [Collection Schema Design](#7-collection-schema-design)
8. [Embedding Strategy](#8-embedding-strategy)
9. [LLM Integration](#9-llm-integration)
10. [Clinical Pipeline Architecture](#10-clinical-pipeline-architecture)
11. [Agent Architecture](#11-agent-architecture)
12. [API Layer](#12-api-layer)
13. [UI Layer](#13-ui-layer)
14. [Observability](#14-observability)
15. [Scaling Considerations](#15-scaling-considerations)

---

## 1. Overview

The Pharmacogenomics Intelligence Agent is a multi-collection retrieval-augmented generation (RAG) system with integrated clinical decision support pipelines. It translates genetic data into actionable drug prescribing recommendations across 25 pharmacogenes and 100+ drugs.

**Key architectural decisions:**

- **Multi-collection over single-collection RAG**: 15 purpose-built collections with domain-specific schemas enable structured filtering alongside semantic search.
- **Knowledge graph augmentation**: Deterministic context from 9 structured dictionaries supplements probabilistic vector retrieval.
- **Clinical pipeline integration**: Star allele calling, phenoconversion, HLA screening, and dosing algorithms run as first-class components, not afterthoughts.
- **Parallel search**: ThreadPoolExecutor searches all 15 collections simultaneously, bounded at 30 merged results.

---

## 2. Three-Tier Architecture

```
+================================================================+
|                      TIER 1: INGEST                             |
|                                                                 |
|  +----------+  +----------+  +----------+  +----------+        |
|  |  CPIC    |  | PharmVar |  | PharmGKB |  |   FDA    |        |
|  |  Parser  |  |  Parser  |  |  Parser  |  |  Parser  |        |
|  +----+-----+  +----+-----+  +----+-----+  +----+-----+        |
|       |             |             |             |               |
|  +----+-----+  +----+-----+  +----+-----+  +----+-----+        |
|  |Population|  |  PubMed  |  | Clinical |  |   Base   |        |
|  |  Parser  |  |  Parser  |  |Trials Par|  |  Parser  |        |
|  +----+-----+  +----+-----+  +----+-----+  +----+-----+        |
|       |             |             |             |               |
|       +-------+-----+------+------+------+------+               |
|               |            |             |                      |
|               v            v             v                      |
|         [to_embedding_text()]  -->  [BGE-small-en-v1.5]         |
|                                         |                      |
|                                    384-dim vector               |
+=========================================|======================+
                                          |
                                          v
+================================================================+
|                    TIER 2: VECTOR STORE                         |
|                                                                 |
|                    Milvus 2.4 (port 19530)                      |
|                                                                 |
|  +------------------+  +------------------+  +--------------+  |
|  | pgx_gene_        |  | pgx_drug_        |  | pgx_drug_    |  |
|  | reference        |  | guidelines       |  | interactions |  |
|  | (W: 0.10)        |  | (W: 0.14)        |  | (W: 0.12)    |  |
|  +------------------+  +------------------+  +--------------+  |
|  +------------------+  +------------------+  +--------------+  |
|  | pgx_hla_         |  | pgx_pheno-       |  | pgx_dosing_  |  |
|  | hypersensitivity  |  | conversion       |  | algorithms   |  |
|  | (W: 0.10)        |  | (W: 0.08)        |  | (W: 0.07)    |  |
|  +------------------+  +------------------+  +--------------+  |
|  +------------------+  +------------------+  +--------------+  |
|  | pgx_clinical_    |  | pgx_population_  |  | pgx_clinical_|  |
|  | evidence         |  | data             |  | trials       |  |
|  | (W: 0.08)        |  | (W: 0.06)        |  | (W: 0.04)    |  |
|  +------------------+  +------------------+  +--------------+  |
|  +------------------+  +------------------+  +--------------+  |
|  | pgx_fda_labels   |  | pgx_drug_        |  | pgx_patient_ |  |
|  | (W: 0.06)        |  | alternatives     |  | profiles     |  |
|  |                   |  | (W: 0.05)        |  | (W: 0.03)    |  |
|  +------------------+  +------------------+  +--------------+  |
|  +------------------+  +------------------+  +--------------+  |
|  | pgx_implement-   |  | pgx_education    |  | genomic_     |  |
|  | ation            |  | (W: 0.02)        |  | evidence     |  |
|  | (W: 0.02)        |  |                   |  | (W: 0.03)    |  |
|  +------------------+  +------------------+  +--------------+  |
|                                                                 |
|  Index: IVF_FLAT  |  Metric: COSINE  |  nlist: 1024            |
+================================================================+
                            |
                            v
+================================================================+
|                    TIER 3: INFERENCE                            |
|                                                                 |
|  +--------------------------+  +---------------------------+   |
|  |      PGxRAGEngine        |  |   PGxIntelligenceAgent    |   |
|  |  - embed query            |  |   - search_plan()        |   |
|  |  - parallel search (15)   |  |   - evaluate_evidence()  |   |
|  |  - query expansion (14)   |  |   - generate_report()    |   |
|  |  - merge & rank (top 30)  |  |                          |   |
|  |  - knowledge augmentation |  |                          |   |
|  +-----------+---------------+  +---------------------------+   |
|              |                                                  |
|  +-----------v-------------------------------------------------+|
|  |                Clinical Pipelines                           ||
|  |  +-------------+ +----------+ +-------+ +----------------+ ||
|  |  |StarAllele   | |Pheno-    | |HLA    | |Dosing          | ||
|  |  |Caller       | |conversion| |Screener| |Calculator     | ||
|  |  |             | |Detector  | |       | | - IWPC warfarin| ||
|  |  |             | |          | |       | | - Tacrolimus   | ||
|  |  |             | |          | |       | | - Fluoropyrim. | ||
|  |  |             | |          | |       | | - Thiopurine   | ||
|  |  +-------------+ +----------+ +-------+ +----------------+ ||
|  +-------------------------------------------------------------+|
|              |                                                  |
|              v                                                  |
|  +-----------------------------------------------------------+ |
|  |  Claude Sonnet 4.6 (Anthropic API)                         | |
|  |  System prompt: 11 PGx expertise domains                   | |
|  |  Streaming response | Citation-linked output                | |
|  |  Max tokens: 2,048 (standard) / 3,000 (comparative)        | |
|  +-----------------------------------------------------------+ |
+================================================================+
```

---

## 3. Data Flow

### 3.1 Query Flow (Standard RAG)

```
User Question (e.g., "CYP2D6 poor metabolizer taking codeine")
     |
     v
[1] Query Embedding
     BGE-small-en-v1.5: "Represent this sentence..."
     Output: 384-dim float vector
     |
     v
[2] Parallel Collection Search
     ThreadPoolExecutor across 15 collections
     Top-K per collection: 5 (configurable)
     Score threshold: 0.4
     |
     v
[3] Query Expansion
     14 PGx expansion maps check for drug/gene/phenotype synonyms
     Expansion terms re-embedded and searched
     |
     v
[4] Merge & Rank
     Deduplicate by record ID
     Weighted scoring: similarity_score * collection_weight
     Cap at 30 results
     |
     v
[5] Knowledge Augmentation
     Extract matching PHARMACOGENES, DRUG_CATEGORIES, HLA_DRUG_ASSOCIATIONS
     Inject structured context (gene function, substrates, inhibitors)
     |
     v
[6] Clinical Pipeline (if applicable)
     Star allele interpretation, phenoconversion check,
     HLA screening, dosing calculation
     |
     v
[7] Prompt Assembly
     Evidence sections (top 5 per collection)
     + Knowledge context
     + Clinical pipeline results
     + Citation formatting instructions
     |
     v
[8] LLM Synthesis
     Claude Sonnet 4.6, streaming, 2,048 tokens
     |
     v
[9] Response
     Answer with clickable citations
     + Evidence panel
     + Clinical alerts (if any)
     + Export buttons (MD/JSON/PDF/FHIR)
```

### 3.2 Ingest Flow

```
External Source (CPIC API / PharmVar / PharmGKB / PubMed / ...)
     |
     v
[1] Parser (one of 8 ingest parsers)
     Fetch raw data via API or load from JSON seed file
     |
     v
[2] Record Construction
     Parse into domain-specific Pydantic model
     Validate all required fields
     |
     v
[3] Embedding Text Generation
     record.to_embedding_text() -> concatenated domain fields
     |
     v
[4] Vector Embedding
     BGE-small-en-v1.5 -> 384-dim vector
     Batch size: 32
     |
     v
[5] Milvus Insert
     Batch insert into target collection
     IVF_FLAT index auto-updated
```

---

## 4. Ingest Tier

### 4.1 Parser Architecture

All 8 parsers inherit from `BaseIngestParser` in `src/ingest/base.py`:

```
BaseIngestParser (abstract)
  |-- fetch()        # Retrieve raw data from source
  |-- parse()        # Transform into Pydantic records
  |-- embed()        # Generate 384-dim vectors
  |-- store()        # Insert into Milvus collection
  |-- run()          # Execute full pipeline
```

### 4.2 Parser Inventory

| Parser | Source | Target Collection | Cadence |
|--------|--------|------------------|---------|
| `cpic_parser.py` | CPIC API | pgx_drug_guidelines | Manual / weekly |
| `pharmvar_parser.py` | PharmVar API | pgx_gene_reference | Manual |
| `pharmgkb_parser.py` | PharmGKB | pgx_drug_interactions | Manual |
| `fda_label_parser.py` | FDA DailyMed | pgx_fda_labels | Manual |
| `population_parser.py` | Population DBs | pgx_population_data | Manual |
| `pubmed_parser.py` | PubMed/NCBI | pgx_clinical_evidence | Weekly (automated) |
| `clinical_trials_parser.py` | ClinicalTrials.gov | pgx_clinical_trials | Weekly (automated) |
| `base.py` | Local JSON files | All 14 PGx collections | One-shot seed |

---

## 5. Vector Store Tier

### 5.1 Milvus Configuration

- **Version**: Milvus 2.4 (standalone)
- **Storage backend**: etcd (metadata) + MinIO (object storage)
- **Index type**: IVF_FLAT (inverted file with flat quantization)
- **Distance metric**: COSINE similarity
- **nlist**: 1024 (number of cluster centroids)
- **nprobe**: 16 (clusters to search at query time)

### 5.2 Collection Design Principles

Each of the 15 collections follows these design principles:

1. **Primary key**: VARCHAR `id` field with domain-specific format (e.g., `CYP2D6_star4` for gene reference)
2. **Embedding field**: FLOAT_VECTOR `embedding` with dim=384
3. **Domain fields**: Collection-specific VARCHAR/FLOAT/INT fields for structured filtering
4. **Text field**: VARCHAR `text` storing the original text used for embedding

### 5.3 Search Weight Rationale

Weights are assigned based on clinical relevance:

- **Drug guidelines (0.14)**: Highest weight because CPIC/DPWG recommendations are the primary actionable output.
- **Drug interactions (0.12)**: PharmGKB annotations provide the evidence basis for gene-drug relationships.
- **Gene reference + HLA (0.10 each)**: Foundational pharmacogene data and safety-critical HLA screening.
- **Clinical evidence + phenoconversion (0.08 each)**: Published outcomes and drug-drug-gene interactions.
- **Dosing algorithms (0.07)**: Quantitative dosing formulas.
- **FDA labels + population data (0.06 each)**: Regulatory labeling and population context.
- **Drug alternatives (0.05)**: Therapeutic substitution options.
- **Clinical trials + genomic evidence (0.04 + 0.03)**: Research context and shared variant data.
- **Patient profiles + implementation + education (0.03 + 0.02 + 0.02)**: Supporting materials.

---

## 6. Inference Tier

### 6.1 PGxRAGEngine (796 lines)

The core RAG engine orchestrates the full retrieval-augmentation-generation pipeline:

1. **Parallel search**: `_search_all_collections()` uses `ThreadPoolExecutor` with up to 15 concurrent workers.
2. **Query expansion**: `_expanded_search()` enriches queries with drug/gene/phenotype synonyms from 14 expansion maps.
3. **Merge and rank**: `_merge_and_rank()` deduplicates and sorts by weighted score.
4. **Knowledge augmentation**: `_get_knowledge_context()` injects structured pharmacogene, drug, HLA, and phenoconversion data.
5. **Prompt assembly**: `_build_prompt()` constructs the evidence-rich prompt for LLM synthesis.
6. **LLM synthesis**: Calls Claude Sonnet 4.6 with streaming support.

### 6.2 Score Calculation

```
final_score = cosine_similarity(query_vector, record_vector) * collection_weight
```

Score thresholds:
- **Minimum**: 0.4 (configurable via `PGX_SCORE_THRESHOLD`)
- **High citation**: >= 0.75
- **Medium citation**: >= 0.60
- **Low citation**: < 0.60

---

## 7. Collection Schema Design

### 7.1 Example: pgx_gene_reference

```
Field Name               Type           Description
─────────────────────────────────────────────────────────────
id                       VARCHAR(100)   Gene-allele ID (e.g., CYP2D6_star1)
embedding                FLOAT_VECTOR   384-dim BGE embedding
gene                     VARCHAR(50)    Pharmacogene symbol
star_allele              VARCHAR(50)    Star allele designation
defining_variants        VARCHAR(1000)  Comma-separated rsIDs
activity_score           FLOAT          Enzyme activity score
function_status          VARCHAR(50)    no/decreased/normal/increased function
allele_frequency_global  FLOAT          Global allele frequency
allele_frequency_european FLOAT         European population frequency
allele_frequency_african FLOAT          African population frequency
allele_frequency_east_asian FLOAT       East Asian population frequency
text                     VARCHAR(4000)  Embedding source text
```

### 7.2 Example: pgx_drug_guidelines

```
Field Name        Type           Description
─────────────────────────────────────────────────────────────
id                VARCHAR(100)   Guideline ID
embedding         FLOAT_VECTOR   384-dim BGE embedding
gene              VARCHAR(50)    Pharmacogene symbol
drug              VARCHAR(200)   Drug name
guideline_body    VARCHAR(20)    CPIC, DPWG, FDA, CPNDS
cpic_level        VARCHAR(10)    A, A/B, B, C, D
phenotype         VARCHAR(100)   Metabolizer phenotype
recommendation    VARCHAR(2000)  Prescribing recommendation
text              VARCHAR(4000)  Embedding source text
```

---

## 8. Embedding Strategy

### 8.1 Model Selection

**Model:** BGE-small-en-v1.5 (BAAI)
- **Dimensions:** 384
- **Max sequence length:** 512 tokens
- **Size:** ~130MB
- **Selection rationale:** Optimal balance of embedding quality and inference speed on DGX Spark. Outperforms larger models on domain-specific retrieval when combined with PGx-optimized embedding text.

### 8.2 Embedding Text Construction

Each collection model implements `to_embedding_text()` to produce domain-optimized text:

```python
# pgx_gene_reference example:
def to_embedding_text(self) -> str:
    return (
        f"Pharmacogene {self.gene} star allele {self.star_allele} "
        f"with {self.function_status} function. "
        f"Activity score: {self.activity_score}. "
        f"Defining variants: {self.defining_variants}."
    )
```

### 8.3 Query Prefix

BGE retrieval requires a task-specific prefix:
```
"Represent this sentence for searching relevant passages: {query}"
```

---

## 9. LLM Integration

### 9.1 Model

**Claude Sonnet 4.6** (claude-sonnet-4-6) via Anthropic API.

### 9.2 System Prompt

11 domains of expertise defined in `PGX_SYSTEM_PROMPT`:
1. Pharmacogene Interpretation
2. Drug-Gene Interaction Analysis
3. Star Allele Nomenclature
4. Diplotype-to-Phenotype Translation
5. CPIC/DPWG/FDA Guideline Application
6. HLA-Mediated Hypersensitivity Screening
7. Phenoconversion Detection
8. Multi-Gene Interaction Modeling
9. Population Pharmacogenetics
10. Dosing Algorithm Application
11. Clinical Implementation

### 9.3 Parameters

| Parameter | Standard | Comparative |
|-----------|---------|-------------|
| Max tokens | 2,048 | 3,000 |
| Temperature | 0.7 | 0.7 |
| Streaming | Yes | Yes |

---

## 10. Clinical Pipeline Architecture

### 10.1 Pipeline Composition

```
VCF Input (optional)
     |
     v
+-------------------+
| StarAlleleCaller  |  pgx_pipeline.py (1,222 lines)
| VCF -> diplotypes |
+--------+----------+
         |
         v
+-------------------+     +----------------------+
|PhenotypeTranslator|     |PhenoconversionDetect | phenoconversion.py (494 lines)
|diplotype -> pheno |     |med list -> adjusted  |
+--------+----------+     +----------+-----------+
         |                           |
         +-------------+-------------+
                       |
                       v
              +--------+--------+
              |DrugGeneMatcher  |  pgx_pipeline.py
              |profile x meds  |
              |-> PGx alerts    |
              +--------+--------+
                       |
         +-------------+-------------+
         |                           |
         v                           v
+--------+--------+         +-------+-------+
|  HLAScreener    |         |DosingCalculator| dosing.py (692 lines)
|  12 HLA-drug    |         | 4 algorithms   |
|  associations   |         |                |
+-----------------+         +----------------+
```

### 10.2 Alert Severity Hierarchy

| Severity | Trigger | UI Display |
|----------|---------|-----------|
| CONTRAINDICATED | HLA positive for fatal ADR, PM + contraindicated drug | Red banner, blocks recommendation |
| MAJOR | PM/UM requiring dose change or drug switch | Orange alert with alternatives |
| MODERATE | IM requiring monitoring or minor adjustment | Yellow advisory |
| MINOR | Known interaction with limited clinical impact | Blue informational |
| INFORMATIONAL | Educational context, population data | Grey note |

---

## 11. Agent Architecture

**File:** `src/agent.py` (761 lines)

The `PGxIntelligenceAgent` implements autonomous reasoning:

```
Question -> search_plan() -> rag.retrieve() -> evaluate_evidence()
                                                     |
                                         sufficient? +--- yes --> rag.query() --> generate_report()
                                                     |
                                                     +--- no  --> sub_questions --> rag.retrieve() (retry)
```

**Workflow detection**: The agent identifies 8 PGx workflow types from keyword analysis and activates the corresponding clinical pipeline (star allele calling, phenoconversion, HLA screening, dosing).

---

## 12. API Layer

### 12.1 FastAPI Application

**File:** `api/main.py` (649 lines)
**Port:** 8107

Architecture:
- Lifespan-managed PGxCollectionManager and PGxRAGEngine initialization
- CORS middleware for cross-origin access from landing page (:8080) and UI (:8507)
- Pydantic request/response validation
- Route modules for PGx clinical, reports, and events

### 12.2 Route Modules

| Module | Path Prefix | Endpoints |
|--------|------------|-----------|
| `pgx_clinical.py` (803 lines) | `/v1/pgx` | 7 PGx-specific clinical endpoints |
| `reports.py` | `/v1/reports` | Report generation |
| `events.py` | `/v1/events` | Event audit trail |

---

## 13. UI Layer

**File:** `app/pgx_ui.py` (2,138 lines)
**Port:** 8507

Streamlit application with 10 specialized tabs, NVIDIA dark theme, cached engine initialization, conversation memory, and export buttons.

---

## 14. Observability

### 14.1 Prometheus Metrics

22 metrics with `pgx_` prefix: 10 histograms, 8 counters, 4 gauges. Exposed via `/metrics` endpoint in Prometheus text format.

### 14.2 Logging

Structured logging via `loguru` (src modules) and `logging` (clinical pipeline modules). Log levels: DEBUG for development, INFO for production.

### 14.3 Health Checks

- FastAPI: `GET /health` returns collection status, vector counts, Milvus connectivity
- Docker: Health check commands with intervals, timeouts, and retry counts for all services

---

## 15. Scaling Considerations

### 15.1 Current Design (DGX Spark)

The system is designed for single-node deployment on NVIDIA DGX Spark:
- **CPU**: 20 ARM cores handle parallel collection search and clinical pipeline logic
- **Memory**: 128GB unified LPDDR5x supports Milvus index + embeddings in-memory
- **GPU**: GB10 accelerates BGE embedding generation

### 15.2 Horizontal Scaling Path

For institutional deployment:
- **Milvus cluster mode**: Replace standalone with distributed Milvus (query nodes, data nodes, index nodes)
- **API replicas**: uvicorn workers (currently 2) can scale behind a load balancer
- **Embedding service**: Dedicated GPU-accelerated embedding microservice
- **LLM gateway**: Rate-limited Anthropic API proxy with request queuing

### 15.3 Data Growth

- **Current seed**: 240 records across 14 JSON files
- **Post-ingest target**: Thousands of records per collection (PubMed, PharmGKB, ClinicalTrials.gov)
- **IVF_FLAT scaling**: Efficient up to ~1M vectors per collection; beyond that, consider HNSW index
