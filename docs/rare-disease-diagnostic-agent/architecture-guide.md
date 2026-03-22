# Rare Disease Diagnostic Agent -- Architecture Guide

**Version:** 1.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Three-Tier Architecture](#2-three-tier-architecture)
3. [Component Map](#3-component-map)
4. [Data Flow](#4-data-flow)
5. [HPO Matching Pipeline](#5-hpo-matching-pipeline)
6. [ACMG Classification Pipeline](#6-acmg-classification-pipeline)
7. [RAG Engine Architecture](#7-rag-engine-architecture)
8. [Workflow Engine](#8-workflow-engine)
9. [Decision Support Layer](#9-decision-support-layer)
10. [Vector Database Design](#10-vector-database-design)
11. [Query Processing Pipeline](#11-query-processing-pipeline)
12. [Cross-Agent Communication](#12-cross-agent-communication)
13. [Deployment Architecture](#13-deployment-architecture)

---

## 1. Architecture Overview

The Rare Disease Diagnostic Agent follows a three-tier architecture designed for modular, gracefully-degrading operation. The system separates presentation (Streamlit UI), application logic (FastAPI + workflows + decision engines), and data (Milvus vector store + knowledge base) into independently scalable layers.

```
+==============================================================+
|                    RARE DISEASE DIAGNOSTIC AGENT               |
+==============================================================+
|                                                                |
|  +------------------+    +-----------------------------+       |
|  |  Streamlit UI    |    |     FastAPI REST API        |       |
|  |  Port 8544       +--->+     Port 8134               |       |
|  |  5 Tabs          |    |     20 Endpoints            |       |
|  |  NVIDIA Theme    |    |     Auth + CORS + Metrics   |       |
|  +------------------+    +----+--------+--------+------+       |
|                               |        |        |              |
|            +------------------+        |        +--------+     |
|            |                           |                 |     |
|  +---------v---------+  +-------------v----+  +---------v--+  |
|  | Workflow Engine   |  | Decision Support |  | RAG Engine |  |
|  | 10 Workflows      |  | 6 Engines        |  | Embed+Search|  |
|  | Template Method   |  | HPO, ACMG, Drug  |  | Fusion     |  |
|  +--------+----------+  | Family, Natural  |  +------+-----+  |
|           |              | History, Algo    |         |        |
|           |              +--------+---------+         |        |
|           |                       |                   |        |
|  +--------v-----------+-----------v---+    +----------v-----+ |
|  |    Knowledge Base                  |    |    Milvus      | |
|  |    13 disease categories           |    |    14 collections| |
|  |    97+ diseases, 45+ genes         |    |    384-dim BGE | |
|  |    12 gene therapies               |    |    IVF_FLAT    | |
|  |    28 ACMG criteria                |    |    COSINE      | |
|  +------------------------------------+    +----------------+ |
+================================================================+
```

---

## 2. Three-Tier Architecture

### Tier 1: Presentation Layer

| Component | Technology | Port | Responsibilities |
|---|---|---|---|
| Streamlit UI | Streamlit + NVIDIA CSS | 8544 | Patient intake, differential display, variant review, therapy matching, report export |

The UI communicates exclusively via HTTP REST calls to the API tier. It maintains no direct database connections, ensuring clean separation and the ability to replace or augment the frontend independently.

### Tier 2: Application Layer

| Component | Technology | Port | Responsibilities |
|---|---|---|---|
| FastAPI API | FastAPI + Uvicorn | 8134 | Request routing, authentication, rate limiting, CORS, metrics |
| Workflow Engine | Python classes | -- | 10 diagnostic workflow orchestration |
| Decision Support | Python classes | -- | 6 clinical decision support engines |
| RAG Engine | Python + sentence-transformers | -- | Embedding, search, context assembly, LLM synthesis |
| Query Expansion | Python module | -- | Entity alias resolution, HPO synonym mapping, term boosting |

### Tier 3: Data Layer

| Component | Technology | Port | Responsibilities |
|---|---|---|---|
| Milvus | Milvus 2.x + etcd + MinIO | 19530 | 14 vector collections, similarity search |
| Knowledge Base | Python dicts (in-memory) | -- | Disease catalogs, gene therapies, ACMG criteria, diagnostic algorithms |
| Ingest Pipeline | Python parsers | -- | HPO, OMIM, Orphanet, gene therapy data ingestion |

---

## 3. Component Map

```
src/
  agent.py              -- Autonomous agent pipeline (plan, search, evaluate, synthesize)
  clinical_workflows.py -- 10 diagnostic workflows (BaseRareDiseaseWorkflow subclasses)
  collections.py        -- 14 Milvus collection schemas and weight maps
  cross_modal.py        -- Cross-agent communication handlers
  decision_support.py   -- 6 decision support engines
  export.py             -- Report generation (Markdown, JSON, PDF)
  knowledge.py          -- Domain knowledge base (13 categories, 97+ diseases)
  metrics.py            -- Prometheus metric definitions
  models.py             -- Pydantic data models (11 enums, 8 models, 1 dataclass)
  query_expansion.py    -- Entity aliases, HPO synonyms, workflow-aware boosting
  rag_engine.py         -- Multi-collection RAG pipeline
  scheduler.py          -- Periodic ingest scheduling
  ingest/
    base.py             -- Abstract base ingest class
    hpo_parser.py       -- HPO ontology parser
    omim_parser.py      -- OMIM disease/gene parser
    orphanet_parser.py  -- Orphanet disease parser
    gene_therapy_parser.py -- Gene therapy catalog parser

api/
  main.py               -- FastAPI application setup, middleware, core endpoints
  routes/
    diagnostic_clinical.py -- Diagnostic and clinical endpoints
    reports.py           -- Report generation endpoints
    events.py            -- SSE event stream

app/
  diagnostic_ui.py      -- 5-tab Streamlit interface

config/
  settings.py           -- Pydantic BaseSettings configuration
```

---

## 4. Data Flow

### 4.1 Diagnostic Query Flow

```
1. User enters HPO terms + clinical notes
          |
2. Streamlit UI --> POST /v1/diagnostic/diagnose
          |
3. FastAPI validates PatientQuery (Pydantic)
          |
4. Workflow detection (auto or specified)
          |
5. Query expansion (aliases, HPO synonyms)
          |
6. Search plan generation (collections, weights, top_k)
          |
7. RAG Engine:
   a. Embed query with BGE-small-en-v1.5
   b. Parallel search across selected collections
   c. Weighted result fusion
   d. Top-K selection and context assembly
          |
8. Decision Support Engines (parallel):
   a. HPO-to-Gene Matcher: BMA similarity scoring
   b. ACMG Classifier: variant classification (if VCF provided)
   c. Orphan Drug Matcher: therapy matching
   d. Diagnostic Algorithm Recommender: test ordering
          |
9. LLM Synthesis (Claude):
   Build evidence context --> Generate diagnostic assessment
          |
10. Response: DiagnosticResult with candidates, variants, therapies
```

### 4.2 Variant Interpretation Flow

```
1. User submits variant data (gene, type, frequency, etc.)
          |
2. POST /v1/diagnostic/variants/interpret
          |
3. ACMG Variant Classifier:
   a. Check BA1 (standalone benign: freq > 5%)
   b. Evaluate pathogenic criteria (PVS1-PP5)
   c. Evaluate benign criteria (BS1-BP7)
   d. Calculate path_score and benign_score
   e. Apply classification thresholds
          |
4. Response: classification, criteria_met, evidence_summary
```

### 4.3 Gene Therapy Matching Flow

```
1. User submits disease + genotype
          |
2. POST /v1/diagnostic/therapy/search
          |
3. Orphan Drug Matcher:
   a. Exact disease match (approved indication)
   b. Pathway match (same gene, different disease)
   c. Repurposing candidates (mechanism-based)
   d. Genotype-specific filtering
          |
4. RAG search for clinical trials (rd_trials, rd_therapies)
          |
5. Response: ranked TherapyMatch list with eligibility
```

---

## 5. HPO Matching Pipeline

The HPO (Human Phenotype Ontology) matching pipeline is the core of phenotype-driven diagnosis. It uses Information Content (IC) scoring with Best-Match-Average (BMA) similarity to rank candidate genes/diseases.

```
Patient HPO Terms                Gene-HPO Associations
[HP:0001250, HP:0001252, ...]    {CFTR: [...], FBN1: [...], ...}
          |                                |
          v                                v
    +-----+----------+            +--------+---------+
    | IC Computation |            | IC Computation   |
    | IC = -log2(p)  |            | for gene terms   |
    +-----+----------+            +--------+---------+
          |                                |
          +---+---+----------+----+--------+
              |              |
     +--------v--------+   +v-----------+
     | Forward Match   |   | Reverse    |
     | avg maxIC P->G  |   | avg maxIC  |
     |                 |   | G->P       |
     +--------+--------+   +-----+------+
              |                   |
              +-------+-----------+
                      |
              +-------v--------+
              | BMA Similarity |
              | 0.5*(fwd+rev)  |
              +-------+--------+
                      |
              +-------v--------+
              | Freq Weighting |
              | 0.7*BMA+0.3*fw |
              +-------+--------+
                      |
              +-------v--------+
              | Ranked Genes   |
              | [SCN1A: 0.82,  |
              |  MECP2: 0.71]  |
              +----------------+
```

### IC Score Examples

| HPO Term | Frequency | IC Score | Discriminating Power |
|---|---|---|---|
| HP:0002816 (Genu recurvatum) | 0.005 | 7.64 | Very high -- specific to connective tissue |
| HP:0001083 (Ectopia lentis) | 0.008 | 6.97 | High -- Marfan/homocystinuria |
| HP:0001250 (Seizures) | 0.15 | 2.74 | Low -- common across many diseases |
| HP:0001252 (Hypotonia) | 0.12 | 3.06 | Low -- non-specific |

---

## 6. ACMG Classification Pipeline

The ACMG/AMP variant classification pipeline implements a simplified but complete scoring system covering 28 criteria.

```
Variant Input
{gene, type, frequency, de_novo, functional, computational, ...}
          |
          v
   +------+-------+
   | BA1 Check    |           freq > 5%? --> BENIGN (standalone)
   +------+-------+
          |
          v
   +------+-------+                    +--------+--------+
   | Pathogenic   |                    | Benign Criteria  |
   | Criteria     |                    |                  |
   | PVS1 (+8)   |                    | BS1: freq > 1%   |
   | PS1  (+4)   |                    | BS2: healthy      |
   | PS2  (+4)   |                    | BP1: missense     |
   | PS3  (+3)   |                    | BP3: repeat       |
   | PM1  (+2)   |                    | BP4: comp benign  |
   | PM2  (+2)   |                    | BP6: source       |
   | PM3  (+2)   |                    | BP7: synonymous   |
   | PM4  (+2)   |                    +--------+----------+
   | PM6  (+1)   |                             |
   | PP1  (+1)   |                             |
   | PP3  (+1)   |                      benign_score
   | PP4  (+1)   |                             |
   | PP5  (+1)   |                             |
   +------+------+                             |
          |                                    |
     path_score                                |
          |                                    |
          +-----+-----+----------+-----+------+
                |                       |
         +------v-----------------------v------+
         |        Classification Decision       |
         |                                      |
         | path >= 10  --> PATHOGENIC            |
         | path >= 6   --> LIKELY PATHOGENIC     |
         | path 1-5    --> VUS                   |
         | benign >= 4 --> LIKELY BENIGN         |
         | benign >= 6 --> BENIGN                |
         +--------------------------------------+
```

### LOF-Intolerant Genes (PVS1 Applicable)

SCN1A, MECP2, KCNQ1, KCNH2, SCN5A, MYH7, MYBPC3, FBN1, COL1A1, COL1A2, HTT, NSD1, NIPBL, KMT2D, CHD7, PTPN11, RAF1, BRAF, FGFR3, TCOF1

---

## 7. RAG Engine Architecture

```
                +-------------------+
                |   Query Input     |
                +--------+----------+
                         |
                +--------v----------+
                | Query Expansion   |
                | (aliases, HPO,    |
                |  workflow boost)  |
                +--------+----------+
                         |
                +--------v----------+
                | BGE-small-en-v1.5 |
                | 384-dim embedding |
                +--------+----------+
                         |
          +--------------+---------------+
          |              |               |
   +------v---+   +-----v----+   +------v---+
   | Coll #1  |   | Coll #2  |   | Coll #N  |   ... 14 collections
   | top_k=50 |   | top_k=30 |   | top_k=10 |   (parallel search)
   | w=0.22   |   | w=0.18   |   | w=0.06   |
   +------+---+   +-----+----+   +------+---+
          |              |               |
          +--------------+---------------+
                         |
                +--------v----------+
                | Weighted Fusion   |
                | score * weight    |
                | rank + deduplicate|
                +--------+----------+
                         |
                +--------v----------+
                | Context Assembly  |
                | top-K results     |
                | citation tracking |
                +--------+----------+
                         |
                +--------v----------+
                | LLM Synthesis     |
                | Claude API        |
                | Evidence-based    |
                | diagnostic report |
                +-------------------+
```

### Search Configuration by Workflow

The RAG engine dynamically adjusts search parameters based on the detected or specified workflow. Each workflow specifies:
- Which collections to search (all 14 by default)
- Per-collection top_k values (5-100)
- Per-collection relevance weights (sum to ~1.0)
- Score threshold for inclusion (default 0.4)

---

## 8. Workflow Engine

### Template Method Pattern

```
BaseRareDiseaseWorkflow (ABC)
    |
    +-- preprocess(inputs)    # Validate, extract HPO terms, detect urgency
    +-- execute(inputs)       # Core workflow logic (abstract)
    +-- postprocess(result)   # Add guideline refs, cross-agent triggers
    +-- run(inputs)           # Template method: preprocess -> execute -> postprocess
```

### Workflow Registration

All 10 workflows are registered in the `WorkflowEngine` for unified dispatch:

```
WorkflowEngine.dispatch(workflow_type, inputs) -> WorkflowResult
```

### Input Validation

Each workflow validates inputs during preprocessing:
- Required HPO terms or clinical notes
- VCF path accessibility (for variant workflows)
- Age/sex validation
- Urgency level assessment
- Validation warnings injected into findings

---

## 9. Decision Support Layer

```
+------------------------------------------------------------------+
|                     Decision Support Layer                         |
|                                                                    |
|  +-------------------+  +---------------------+  +--------------+ |
|  | HPOToGeneMatcher  |  | ACMGVariantClassifier|  | OrphanDrug   | |
|  | BMA + IC scoring  |  | 28-criteria scoring  |  | Matcher      | |
|  | 14 gene-HPO maps  |  | 20 LOF-intolerant    |  | 12+ drugs    | |
|  | 40+ HPO terms     |  | 4 hot spot genes     |  | 3 match types| |
|  +-------------------+  +---------------------+  +--------------+ |
|                                                                    |
|  +-------------------+  +---------------------+  +--------------+ |
|  | DiagnosticAlgo    |  | FamilySegregation   |  | NaturalHist  | |
|  | Recommender       |  | Analyzer            |  | Predictor    | |
|  | 6 cluster pathways|  | LOD score calc      |  | 6 diseases   | |
|  | Ordered test recs |  | AD/AR/XL support    |  | 24+ milestones| |
|  +-------------------+  +---------------------+  +--------------+ |
+------------------------------------------------------------------+
```

---

## 10. Vector Database Design

### Collection Organization

Collections are organized by information type, with each collection optimized for a specific domain query pattern:

```
Phenotype Layer:     rd_phenotypes (18K) -- HPO terms with IC scores
                     rd_case_reports (20K) -- Phenotype-genotype case data

Disease Layer:       rd_diseases (10K) -- OMIM/Orphanet disease entries
                     rd_natural_history (5K) -- Disease milestones

Genetic Layer:       rd_genes (5K) -- Gene-disease associations
                     rd_variants (500K) -- ACMG-classified variants
                     genomic_evidence (3.56M) -- Shared genomic data

Knowledge Layer:     rd_literature (50K) -- Published literature
                     rd_guidelines (3K) -- Clinical guidelines
                     rd_pathways (2K) -- Metabolic/signaling pathways

Clinical Layer:      rd_therapies (2K) -- Approved/investigational therapies
                     rd_trials (8K) -- Clinical trials
                     rd_registries (1.5K) -- Patient registries
                     rd_newborn_screening (80) -- NBS conditions
```

---

## 11. Query Processing Pipeline

```
1. Raw Query --> 2. Entity Detection --> 3. Alias Resolution
       |                                        |
4. HPO Term Extraction                   5. Workflow Detection
       |                                        |
6. Query Expansion (synonyms, related terms)
       |
7. Search Plan Generation
   - Select collections
   - Assign weights
   - Set per-collection top_k
       |
8. Parallel Collection Search
       |
9. Weighted Result Fusion
       |
10. Decision Support Enrichment
       |
11. LLM Context Assembly
       |
12. Claude Synthesis
       |
13. Structured Response (DiagnosticResult / WorkflowResult)
```

---

## 12. Cross-Agent Communication

```
Rare Disease Agent (:8134)
       |
       +---> Genomics Pipeline (:8527)    [VCF variant data]
       |
       +---> PGx Agent (:8107)            [Drug metabolism]
       |
       +---> Cardiology Agent (:8126)     [Cardiac genetics]
       |
       +---> Biomarker Agent (:8529)      [Biomarker stratification]
       |
       +---> Clinical Trial Agent (:8538) [Trial eligibility]
```

Cross-agent triggers are generated automatically by workflow logic:
- Detected cardiac phenotypes -> Cardiology Agent referral
- Identified gene therapy candidate -> Clinical Trial Agent query
- Pharmacogenomic implications -> PGx Agent consultation

Communication uses HTTP REST with configurable timeout (default 30s) and graceful failure handling.

---

## 13. Deployment Architecture

### Docker Compose Stack

```
+------------------------------------------------------------------+
|                         Docker Compose                            |
|                                                                    |
|  +-------------------+  +---------------------+                   |
|  | rare-disease-api  |  | rare-disease-ui     |                   |
|  | :8134             |  | :8544               |                   |
|  | FastAPI + Uvicorn  |  | Streamlit           |                   |
|  +--------+----------+  +----------+----------+                   |
|           |                        |                               |
|  +--------v------------------------v----------+                   |
|  |              Shared Network                 |                   |
|  +--------+-----------------------------------+                   |
|           |                                                        |
|  +--------v----------+  +----------+  +----------+               |
|  | Milvus Standalone |  | etcd     |  | MinIO    |               |
|  | :19530            |  | :2379    |  | :9000    |               |
|  +-------------------+  +----------+  +----------+               |
+------------------------------------------------------------------+
```

### Graceful Degradation Modes

| Mode | Available Features | Missing Dependency |
|---|---|---|
| Full | All features (RAG + LLM + search + workflows) | None |
| Search-Only | Vector search + workflows + decision support | Anthropic API key |
| Knowledge-Only | Knowledge base + decision support + workflows | Milvus |
| Offline | Decision support engines + knowledge lookups | Milvus + Anthropic |

---

*Apache 2.0 License -- HCLS AI Factory*
