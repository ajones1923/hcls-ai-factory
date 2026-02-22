# CAR-T Intelligence Agent — Project Bible

> **Purpose:** Complete implementation reference for building the HCLS CAR-T Intelligence Agent on NVIDIA DGX Spark. Import this document into a Claude Code session as context for implementation.
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
13. [Agent Reasoning Pipeline](#13-agent-reasoning-pipeline)
14. [Data Ingest Pipelines](#14-data-ingest-pipelines)
15. [FastAPI REST Server](#15-fastapi-rest-server)
16. [Streamlit Chat UI](#16-streamlit-chat-ui)
17. [Report Export](#17-report-export)
18. [Monitoring & Metrics](#18-monitoring--metrics)
19. [Testing Strategy](#19-testing-strategy)
20. [HCLS AI Factory Integration](#20-hcls-ai-factory-integration)
21. [Implementation Sequence](#21-implementation-sequence)

---

## 1. Project Overview & Goals

### What This Agent Does

The CAR-T Intelligence Agent is a cross-functional research intelligence system that unifies CAR-T cell therapy knowledge across **11 Milvus collections** spanning the entire development lifecycle — from target identification through post-market pharmacovigilance. It uses retrieval-augmented generation to scan curated biomedical collections and generate structured, traceable reports with clickable PubMed and ClinicalTrials.gov citations.

### Five Development Stages

| Stage | Scope |
|---|---|
| **Target Identification** | Antigen biology, expression profiling, disease association |
| **CAR Design** | scFv selection, costimulatory domains, signaling architecture |
| **Vector Engineering** | Transduction, viral vector production, manufacturing processes |
| **In Vitro / In Vivo Testing** | Cytotoxicity, cytokine assays, animal models, persistence |
| **Clinical Development** | Trial design, response rates, toxicity management |

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | **6,266** across 11 Milvus collections (10 owned + 1 read-only) |
| Multi-collection search latency | **12-16 ms** (11 collections, top-5 each, cached) |
| Comparative dual retrieval | **~365 ms** (2 × 11 collections, entity-filtered) |
| Full RAG query (search + Claude) | **~24 sec** end-to-end |
| Comparative RAG query | **~30 sec** end-to-end |
| Cosine similarity scores | **0.74 – 0.90** on demo queries |
| Knowledge graph entities | **70+** (25 targets, 8 toxicities, 10 manufacturing, 15 biomarkers, 6 regulatory) |
| Query expansion | **12 maps**, 169 keywords → 1,496 terms |
| Entity aliases | **39+** for comparative resolution |

### Pipeline Pattern

Every query follows the same pattern:

1. User question submitted (Streamlit UI or REST API)
2. Agent plans search strategy (identify antigens, stages, decompose complex queries)
3. BGE-small embeds query with asymmetric instruction prefix
4. Parallel search across 11 Milvus collections (ThreadPoolExecutor)
5. Query expansion adds semantically related terms
6. Knowledge graph augments results with domain context
7. Evidence merged, deduplicated, ranked (top 30)
8. Claude Sonnet 4 generates grounded answer with citations
9. Response streamed to user with evidence panel

### HCLS AI Factory Integration

This agent is one node in the broader HCLS AI Factory. Cross-agent triggers include:

- **CAR-T → Genomics (Parabricks):** Shared Milvus instance, genomic_evidence collection (3.56M vectors, read-only)
- **CAR-T → Drug Discovery (BioNeMo):** Target validation feeds molecule generation
- **CAR-T → Imaging Agent:** Cross-modal phenotype correlation (future)
- **CAR-T → Biomarker Agent:** CRS prediction, exhaustion monitoring (future)

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

- **128 GB unified memory** eliminates CPU↔GPU transfer bottleneck
- **NVLink-C2C** provides 900 GB/s bandwidth between Grace CPU and Blackwell GPU
- **ARM64 Grace CPU** runs all Python services natively
- **Single-box deployment** — Milvus, embeddings, LLM inference, and UI on one machine
- **$3,999 price point** — accessible for research labs and departmental pilots

---

## 3. Repository Layout

```
cart_intelligence_agent/
├── Docs/
│   ├── CART_Intelligence_Agent_Design.md         # Architecture design document
│   ├── CART_INTELLIGENCE_AGENT_PROJECT_BIBLE.md   # This document
│   ├── CART_Intelligence_Agent_DGX_Spark_Infographic_Prompt.md
│   └── CART_Intelligence_Agent_VAST_AI_OS_Infographic_Prompt.md
├── docs/
│   ├── ARCHITECTURE_GUIDE.md                     # Architecture walkthrough
│   ├── PROJECT_BIBLE.md                          # HCLS Factory project bible
│   ├── WHITE_PAPER.md                            # Technical white paper
│   ├── LEARNING_GUIDE_FOUNDATIONS.md             # Foundations guide
│   ├── LEARNING_GUIDE_ADVANCED.md                # Advanced guide
│   └── DEMO_GUIDE.md                            # Demo walkthrough
├── src/
│   ├── __init__.py
│   ├── models.py                    # Pydantic data models (474 lines)
│   ├── collections.py              # Milvus collection schemas + manager (1,004 lines)
│   ├── knowledge.py                # Knowledge graph (1,511 lines)
│   ├── query_expansion.py          # 12 expansion maps (1,257 lines)
│   ├── rag_engine.py               # Multi-collection RAG engine (686 lines)
│   ├── agent.py                    # Agent reasoning pipeline (262 lines)
│   ├── export.py                   # Report export: MD, JSON, PDF (903 lines)
│   ├── metrics.py                  # Prometheus metrics (189 lines)
│   ├── scheduler.py                # Data ingestion scheduler (226 lines)
│   ├── ingest/
│   │   ├── __init__.py
│   │   ├── base.py                 # Base ingest pipeline (184 lines)
│   │   ├── literature_parser.py    # PubMed E-utilities ingest (350 lines)
│   │   ├── clinical_trials_parser.py  # ClinicalTrials.gov API v2 (403 lines)
│   │   ├── construct_parser.py     # CAR construct parser + FDA seed (292 lines)
│   │   ├── assay_parser.py         # Assay data parser (163 lines)
│   │   └── manufacturing_parser.py # Manufacturing/CMC parser (111 lines)
│   └── utils/
│       ├── __init__.py
│       └── pubmed_client.py        # NCBI E-utilities HTTP client (390 lines)
├── app/
│   └── cart_ui.py                  # Streamlit chat + comparative UI (1,123 lines)
├── api/
│   └── main.py                     # FastAPI REST server (555 lines)
├── config/
│   └── settings.py                 # Pydantic BaseSettings (101 lines)
├── data/
│   ├── reference/
│   │   ├── assay_seed_data.json    # 45 curated assay records
│   │   └── manufacturing_seed_data.json  # 30 curated manufacturing records
│   └── cache/                      # Embedding cache
├── scripts/
│   ├── setup_collections.py        # Create collections + seed FDA constructs
│   ├── ingest_pubmed.py            # CLI: PubMed ingest
│   ├── ingest_clinical_trials.py   # CLI: ClinicalTrials.gov ingest
│   ├── seed_assays.py              # CLI: Seed assay data
│   ├── seed_manufacturing.py       # CLI: Seed manufacturing data
│   ├── validate_e2e.py             # End-to-end validation (5 tests)
│   ├── test_rag_pipeline.py        # Full RAG + LLM integration test
│   └── seed_knowledge.py           # Knowledge graph export
├── docker-compose.yml              # Full deployment (6 services)
├── Dockerfile                      # Multi-stage build (86 lines)
├── requirements.txt                # Python dependencies
├── .env.example                    # Environment variable template
├── LICENSE                         # Apache 2.0
└── README.md
```

**55 Python files | ~16,748 lines of code | Apache 2.0**

---

## 4. Docker Compose Services

### Service Architecture

| Service | Image | Port | Role |
|---|---|---|---|
| `milvus-etcd` | quay.io/coreos/etcd:v3.5.5 | 2379 | Milvus metadata store |
| `milvus-minio` | minio/minio:v2023.03 | 9000, 9001 | Milvus object storage |
| `milvus-standalone` | milvusdb/milvus:v2.4 | 19530, 9091 | Vector database |
| `cart-streamlit` | Built from Dockerfile | 8521 | Streamlit chat UI |
| `cart-api` | Built from Dockerfile | 8522 | FastAPI REST server |
| `cart-setup` | Built from Dockerfile | — | One-shot collection setup + seed |

### Startup Sequence

```bash
# 1. Configure
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY

# 2. Launch all services
docker compose up -d

# 3. Watch setup completion
docker compose logs -f cart-setup

# 4. Verify
curl http://localhost:8522/health
# → {"status": "healthy", "collections": 11, "total_vectors": 6266}
```

### Dockerfile (Multi-Stage)

**Stage 1 — Builder:**
- Base: `python:3.10-slim`
- System deps: `build-essential`, `gcc`, `g++`, `libxml2-dev`, `libxslt1-dev`
- Virtual environment with all pip dependencies

**Stage 2 — Runtime:**
- Base: `python:3.10-slim`
- Copy venv from builder
- Copy source: `config/`, `src/`, `app/`, `api/`, `scripts/`, `data/`
- Non-root user: `cartuser`
- Healthcheck: `curl -f http://localhost:8521/_stcore/health`

---

## 5. Milvus Collection Schemas

### 11 Collections Overview

| Collection | Records | Primary Fields | Source |
|---|---|---|---|
| `cart_literature` | 5,047 | PMID, title, text_chunk, year, cart_stage, target_antigen | PubMed E-utilities |
| `cart_trials` | 973 | NCT ID, title, phase, status, sponsor, target_antigen | ClinicalTrials.gov API v2 |
| `cart_constructs` | 6 | name, target_antigen, scfv_origin, costimulatory, generation | FDA product labels |
| `cart_assays` | 45 | assay_type, construct_id, target_antigen, key_metric, outcome | Landmark publications |
| `cart_manufacturing` | 30 | process_step, vector_type, parameter, target_spec, met_spec | Published CMC data |
| `cart_safety` | 40 | product, event_type, severity_grade, incidence_rate, management | Pharmacovigilance |
| `cart_biomarkers` | 43 | biomarker_name, biomarker_type, clinical_cutoff, evidence_level | CRS/response markers |
| `cart_regulatory` | 25 | product, regulatory_event, date, agency, indication, decision | FDA/EMA milestones |
| `cart_sequences` | 27 | construct_name, scfv_clone, binding_affinity_kd, species_origin | Molecular data |
| `cart_realworld` | 30 | study_type, data_source, product, population_size, outcome_value | Registry/RWE outcomes |
| `genomic_evidence` | 3,561,170 | chrom, pos, gene, consequence, clinical_significance | Shared (read-only) |

### Unified Schema Pattern

All collections share this structure:

```python
# Primary key
id: VARCHAR(100)              # PMID, NCT ID, construct ID, etc.

# Embedding vector (indexed)
embedding: FLOAT_VECTOR(384)  # BGE-small-en-v1.5
# Index: IVF_FLAT, metric: COSINE, nlist: 1024, nprobe: 16

# Text fields
title / text_summary / text_chunk: VARCHAR(500-3000)

# Collection-specific metadata
# VARCHAR for enums, INT64 for years/counts, FLOAT for metrics
```

### Index Configuration

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

### CARTCollectionManager

| Method | Purpose |
|---|---|
| `create_all_collections()` | Create all 11 collections with IVF_FLAT indexes |
| `get_collection_stats()` | Return record count per collection |
| `insert_batch()` | Bulk insert records with Enum→str conversion and UTF-8 truncation |
| `search()` | Single-collection vector similarity search with optional filters |
| `search_all()` | **Parallel** multi-collection search via ThreadPoolExecutor (11 threads) |

---

## 6. Pydantic Data Models

### Enumeration Types

```python
CARTStage:        TARGET_ID, CAR_DESIGN, VECTOR_ENG, TESTING, CLINICAL
TrialPhase:       EARLY_1, PHASE_1, PHASE_1_2, PHASE_2, PHASE_2_3, PHASE_3, PHASE_4, NA
TrialStatus:      RECRUITING, ACTIVE, COMPLETED, TERMINATED, WITHDRAWN, SUSPENDED, NOT_YET
CARGeneration:    1st, 2nd, 3rd, 4th, armored, universal
AssayType:        CYTOTOXICITY, CYTOKINE, FLOW_CYTOMETRY, PROLIFERATION, IN_VIVO, PERSISTENCE, EXHAUSTION
ProcessStep:      TRANSDUCTION, EXPANSION, HARVEST, FORMULATION, RELEASE_TESTING, CRYOPRESERVATION
SafetyEventType:  CRS, ICANS, CYTOPENIA, INFECTION, SECONDARY_MALIGNANCY, ORGAN_TOXICITY, NEUROLOGIC, CARDIAC
BiomarkerType:    PREDICTIVE, PROGNOSTIC, PHARMACODYNAMIC, MONITORING, RESISTANCE
EvidenceLevel:    VALIDATED, EMERGING, EXPLORATORY
RegulatoryEvent:  BLA, BREAKTHROUGH_THERAPY, RMAT, ACCELERATED_APPROVAL, FULL_APPROVAL, LABEL_UPDATE, REMS
RWEStudyType:     RETROSPECTIVE, REGISTRY, CLAIMS, EHR_ANALYSIS, META_ANALYSIS
```

### Collection Models

Each collection has a corresponding Pydantic model with a `to_embedding_text()` method that combines key fields for BGE encoding:

| Model | Collection | Key Fields |
|---|---|---|
| `CARTLiterature` | cart_literature | PMID, title, text_chunk, year, cart_stage, target_antigen, disease, journal |
| `ClinicalTrial` | cart_trials | NCT ID, title, phase, status, sponsor, target_antigen, car_generation, costimulatory |
| `CARConstruct` | cart_constructs | name, target_antigen, scfv_origin, costimulatory_domain, generation, fda_status |
| `AssayResult` | cart_assays | assay_type, construct_id, target_antigen, cell_line, key_metric, outcome |
| `ManufacturingRecord` | cart_manufacturing | process_step, vector_type, parameter, target_spec, met_spec |
| `SafetyRecord` | cart_safety | product, event_type, severity_grade, incidence_rate, management_protocol |
| `BiomarkerRecord` | cart_biomarkers | biomarker_name, biomarker_type, assay_method, clinical_cutoff, evidence_level |
| `RegulatoryRecord` | cart_regulatory | product, regulatory_event, date, agency, indication, decision |
| `SequenceRecord` | cart_sequences | construct_name, scfv_clone, binding_affinity_kd, species_origin |
| `RealWorldRecord` | cart_realworld | study_type, data_source, product, population_size, primary_endpoint |

### Search Result Models

```python
SearchHit:              # Single evidence item (collection, id, score, text, metadata)
CrossCollectionResult:  # Merged multi-collection results (hits, knowledge_context, search_time_ms)
ComparativeResult:      # Dual-entity comparison (entity_a, entity_b, evidence_a, evidence_b)
AgentQuery:             # User input (question, target_antigen, cart_stage, include_genomic)
AgentResponse:          # Agent output (question, answer, evidence, knowledge_used, timestamp)
```

---

## 7. Configuration Reference

### CARTSettings (Pydantic BaseSettings)

```python
# Environment prefix: CART_
# Loaded from: .env file or environment variables

# Paths
PROJECT_ROOT: Path          # Auto-detected
DATA_DIR: Path              # PROJECT_ROOT / "data"
RAG_PIPELINE_ROOT: Path     # For shared resources

# Milvus
MILVUS_HOST: str = "localhost"
MILVUS_PORT: int = 19530

# Embeddings
EMBEDDING_MODEL: str = "BAAI/bge-small-en-v1.5"
EMBEDDING_DIMENSION: int = 384
EMBEDDING_BATCH_SIZE: int = 32

# LLM
LLM_PROVIDER: str = "anthropic"
LLM_MODEL: str = "claude-sonnet-4-20250514"
ANTHROPIC_API_KEY: Optional[str] = None

# RAG Search
TOP_K_PER_COLLECTION: int = 5
SCORE_THRESHOLD: float = 0.4

# Collection Weights (sum ≈ 1.0)
WEIGHT_LITERATURE: float = 0.20
WEIGHT_TRIALS: float = 0.16
WEIGHT_CONSTRUCTS: float = 0.10
WEIGHT_ASSAYS: float = 0.09
WEIGHT_MANUFACTURING: float = 0.07
WEIGHT_SAFETY: float = 0.08
WEIGHT_BIOMARKERS: float = 0.08
WEIGHT_REGULATORY: float = 0.06
WEIGHT_SEQUENCES: float = 0.06
WEIGHT_REALWORLD: float = 0.07
WEIGHT_GENOMIC: float = 0.04

# External APIs
NCBI_API_KEY: Optional[str] = None
PUBMED_MAX_RESULTS: int = 5000
CT_GOV_BASE_URL: str = "https://clinicaltrials.gov/api/v2"

# Services
API_HOST: str = "0.0.0.0"
API_PORT: int = 8522
STREAMLIT_PORT: int = 8521

# Scheduling
INGEST_SCHEDULE_HOURS: int = 168     # Weekly refresh
INGEST_ENABLED: bool = False

# RAG Context
MAX_CONVERSATION_CONTEXT: int = 3    # Prior exchanges
CITATION_HIGH_THRESHOLD: float = 0.75
CITATION_MEDIUM_THRESHOLD: float = 0.60
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
prefix = "Represent this sentence for searching relevant passages: "
query_embedding = embedder.encode(prefix + question)

# Document embedding (raw text, no prefix)
doc_embedding = embedder.encode(record.to_embedding_text())
```

This asymmetric approach improves retrieval relevance by 5-15% compared to symmetric encoding.

### to_embedding_text() Pattern

Each Pydantic model generates its embedding input by combining key fields:

```python
# CARTLiterature
def to_embedding_text(self) -> str:
    parts = [self.title, self.text_chunk]
    if self.target_antigen:
        parts.append(f"Target: {self.target_antigen}")
    if self.disease:
        parts.append(f"Disease: {self.disease}")
    return " ".join(parts)
```

---

## 9. Multi-Collection RAG Engine

### System Prompt

```python
CART_SYSTEM_PROMPT = """You are a CAR-T cell therapy intelligence agent with deep expertise in:

1. Target Identification — antigen biology, expression profiling
2. CAR Design — scFv selection, costimulatory domains (CD28 vs 4-1BB)
3. Vector Engineering — lentiviral/retroviral production, VCN optimization
4. In Vitro & In Vivo Testing — cytotoxicity, cytokine profiling, persistence
5. Clinical Development — trial design, response rates, toxicity management
6. Manufacturing — leukapheresis, expansion, cryopreservation, release testing
7. Safety & Pharmacovigilance — REMS, FAERS, long-term follow-up
8. Biomarkers — CRS prediction (ferritin, CRP, IL-6), MRD monitoring
9. Regulatory Intelligence — FDA pathways, BLA timelines, breakthrough therapy
10. Molecular Design — scFv binding affinity, CDR sequences, humanization
11. Real-World Evidence — CIBMTR registry, community vs academic outcomes
12. Genomic Evidence — patient variants, ClinVar, AlphaMissense

When answering:
- Cite evidence using clickable markdown links
- Think cross-functionally — connect insights across stages
- Highlight failure modes and resistance mechanisms
- Be specific — cite trial names, products, quantitative data
- Acknowledge uncertainty — distinguish facts from emerging data"""
```

### Collection Configuration

```python
COLLECTION_CONFIG = {
    "cart_literature":    {"weight": 0.20, "label": "Literature",    "has_target_antigen": True},
    "cart_trials":        {"weight": 0.16, "label": "Trial",         "has_target_antigen": True},
    "cart_constructs":    {"weight": 0.10, "label": "Construct",     "has_target_antigen": True},
    "cart_assays":        {"weight": 0.09, "label": "Assay",         "has_target_antigen": True},
    "cart_manufacturing": {"weight": 0.07, "label": "Manufacturing", "has_target_antigen": False},
    "cart_safety":        {"weight": 0.08, "label": "Safety",        "has_target_antigen": False},
    "cart_biomarkers":    {"weight": 0.08, "label": "Biomarker",     "has_target_antigen": True},
    "cart_regulatory":    {"weight": 0.06, "label": "Regulatory",    "has_target_antigen": False},
    "cart_sequences":     {"weight": 0.06, "label": "Sequence",      "has_target_antigen": True},
    "cart_realworld":     {"weight": 0.07, "label": "RealWorld",     "has_target_antigen": False},
    "genomic_evidence":   {"weight": 0.04, "label": "Genomic",       "has_target_antigen": False},
}
```

### Core Methods

| Method | Purpose | Returns |
|---|---|---|
| `retrieve()` | Parallel search + expansion + knowledge augmentation | `CrossCollectionResult` |
| `query()` | Full RAG: retrieve + Claude generation | `str` (answer) |
| `query_stream()` | Streaming RAG: yields evidence then tokens | `Generator[Dict]` |
| `find_related()` | Cross-collection entity linking | `Dict[str, List[SearchHit]]` |
| `retrieve_comparative()` | Dual-entity comparative retrieval | `ComparativeResult` |

### Retrieve Pipeline

```
User Query
    ↓
1. Embed query (BGE asymmetric prefix)                        [< 5 ms]
    ↓
2. Parallel search across 11 collections (top-5 each)        [12-16 ms]
   (ThreadPoolExecutor, 11 concurrent threads)
    ↓
3. Query expansion: detect keywords → expand → re-search     [8-12 ms]
    ↓
4. Merge + deduplicate + weighted rank (cap at 30)            [< 1 ms]
    ↓
5. Score citations: high (≥0.75) / medium (≥0.60) / low      [< 1 ms]
    ↓
6. Knowledge graph augmentation (all domains)                 [< 1 ms]
    ↓
7. Build prompt: evidence + knowledge + question              [< 1 ms]
    ↓
8. Stream Claude Sonnet 4 response                            [~22-24 sec]
```

### Citation Formatting

```python
# Literature → PubMed link
"[Literature:PMID](https://pubmed.ncbi.nlm.nih.gov/PMID/)"

# Trial → ClinicalTrials.gov link
"[Trial:NCT...](https://clinicaltrials.gov/study/NCT...)"
```

---

## 10. Knowledge Graph

### Components

| Component | Count | Description |
|---|---|---|
| **Target Antigens** | 25 | Full profiles: protein, UniProt ID, expression, diseases, approved products, resistance, toxicity |
| **Toxicity Profiles** | 8 | CRS, ICANS, B-cell aplasia, HLH/MAS, cytopenias, TLS, GvHD, on-target/off-tumor |
| **Manufacturing Processes** | 10 | Transduction, expansion, harvest, formulation, release testing, cryopreservation |
| **Biomarker Database** | 15+ | CRS prediction (ferritin, CRP, IL-6), exhaustion (PD-1, LAG-3, TIM-3), response markers |
| **Regulatory Pathways** | 6 | FDA breakthrough, RMAT, accelerated approval, REMS, post-marketing |
| **Entity Aliases** | 39+ | Product names, generic names, costimulatory domains, biomarker terms, regulatory terms |

### Target Antigens (25 entries)

```
B-cell:    CD19, CD22, CD20, BCMA, GPRC5D, CD38, CD70
Myeloid:   CD33, CD123, FLT3, CLL-1
Solid:     HER2, GPC3, EGFR, EGFRvIII, Mesothelin, Claudin18.2, GD2,
           PSMA, ROR1, DLL3, B7-H3, MUC1, IL13RA2, EpCAM
```

Each target includes: protein name, UniProt ID, expression pattern, associated diseases, FDA-approved products, key trials, known resistance mechanisms, toxicity profile, and normal tissue expression.

### Toxicity Profiles

| Toxicity | Grading | Management |
|---|---|---|
| **CRS** | CTCAE Grade 1-4 | Tocilizumab, corticosteroids |
| **ICANS** | ICE score-based | Dexamethasone, supportive |
| **B-cell aplasia** | On-target depletion | IVIG replacement |
| **HLH/MAS** | Ferritin + fibrinogen | Anakinra, etoposide |
| **Cytopenias** | CBC monitoring | G-CSF, transfusions |
| **TLS** | Labs (uric acid, K+, PO4) | Rasburicase, hydration |
| **GvHD** | Allogeneic only | Steroids, ruxolitinib |
| **On-target/off-tumor** | Normal tissue expression | Affinity tuning, safety switches |

### API Functions

```python
get_target_context("CD19")              # Full CD19 knowledge block
get_toxicity_context("CRS")             # CRS grading + management
get_manufacturing_context("expansion")   # Expansion process specs
get_biomarker_context("ferritin")        # CRS prediction data
get_regulatory_context("Kymriah")        # FDA approval timeline
resolve_comparison_entity("Kymriah")     # → {"type": "product", "target": "CD19"}
get_comparison_context(entity_a, entity_b)  # Side-by-side knowledge
get_knowledge_stats()                    # {target_antigens: 25, ...}
```

### Entity Resolution Priority

```
CART_TARGETS (25) → ENTITY_ALIASES (39+) → CART_TOXICITIES (8) → CART_MANUFACTURING (10)
```

---

## 11. Query Expansion

### 12 Expansion Map Categories

| Category | Keywords | Expanded Terms | Examples |
|---|---|---|---|
| Target Antigen | 26 | 196 | CD19 → [B-ALL, DLBCL, tisagenlecleucel, axicabtagene, FMC63, ...] |
| Disease | 16 | 143 | multiple myeloma → [MM, RRMM, plasma cell neoplasm, ...] |
| Toxicity | 14 | 136 | CRS → [cytokine release syndrome, tocilizumab, IL-6, ...] |
| Manufacturing | 16 | 181 | transduction → [lentiviral, VCN, MOI, viral vector, ...] |
| Mechanism | 19 | 224 | resistance → [antigen loss, lineage switch, trogocytosis, ...] |
| Construct | 20 | 206 | bispecific → [dual-targeting, tandem, bicistronic, ...] |
| Safety | 15 | 135 | REMS → [risk evaluation, FAERS, adverse event, ...] |
| Biomarker | 14 | 125 | CRS prediction → [ferritin, CRP, IL-6, sIL-2R, ...] |
| Regulatory | 12 | 108 | BLA → [biologics license application, RMAT, ...] |
| Sequence | 10 | 95 | scFv → [single-chain variable fragment, VH, VL, CDR, ...] |
| Real-World | 12 | 110 | registry → [CIBMTR, real-world evidence, post-marketing, ...] |
| Immunogenicity | 10 | 92 | ADA → [anti-drug antibody, neutralizing antibody, ...] |
| **Total** | **169** | **1,496** | |

### How It Works

```python
def expand_query(question: str) -> List[str]:
    """Detect keywords in user question, return expansion terms.

    Used by RAG engine to broaden recall via semantic re-search.
    Returns up to 10-20 related terms per matched keyword.
    """
```

---

## 12. Comparative Analysis Mode

### Auto-Detection

Comparative queries are detected by regex matching: `compare`, `vs`, `versus`, `comparing`, `X vs Y`.

### Pipeline

```
User: "Compare 4-1BB vs CD28 costimulatory domains"
    ↓
1. _is_comparative() detects "vs"                              [< 1 ms]
    ↓
2. _parse_comparison_entities() extracts "4-1BB", "CD28"       [< 1 ms]
    ↓
3. resolve_comparison_entity() for each entity                  [< 1 ms]
   "4-1BB" → {"type": "costimulatory", "canonical": "4-1BB (CD137)"}
   "CD28"  → {"type": "costimulatory", "canonical": "CD28"}
    ↓
4. Dual retrieve() — one per entity                             [~365 ms]
   Entity A: retrieve(question, target_antigen=entity_a.target)
   Entity B: retrieve(question, target_antigen=entity_b.target)
    ↓
5. get_comparison_context() — side-by-side knowledge graph     [< 1 ms]
    ↓
6. _build_comparative_prompt() — structured comparison         [< 1 ms]
   Evidence grouped by entity (A section, B section)
   Instructions: table + advantages + limitations
    ↓
7. Stream Claude Sonnet 4 (max_tokens=3000)                    [~28-30 sec]
   Output: comparison table, pros/cons, clinical context
```

### Supported Entity Types

| Type | Examples | Resolution |
|---|---|---|
| Target Antigens | CD19, BCMA, CD22 | Direct match in CART_TARGETS (25) |
| FDA Products | Kymriah, Yescarta, Carvykti | ENTITY_ALIASES → canonical + target |
| Costimulatory Domains | 4-1BB, CD28 | ENTITY_ALIASES → type: costimulatory |
| Toxicity Profiles | CRS, ICANS | CART_TOXICITIES → type: toxicity |
| Manufacturing | Lentiviral, transduction | CART_MANUFACTURING → type: manufacturing |

### Fallback

If entity parsing fails, the query gracefully falls back to the standard single-query `retrieve()` path.

---

## 13. Agent Reasoning Pipeline

### CARTIntelligenceAgent

```python
class CARTIntelligenceAgent:
    """Autonomous reasoning agent: Plan → Search → Evaluate → Synthesize → Report"""
```

### SearchPlan

```python
@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str]
    target_antigens: List[str]
    relevant_stages: List[CARTStage]
    search_strategy: str           # "broad" | "targeted" | "comparative"
    sub_questions: List[str]       # Decomposed for complex queries
```

### Agent Pipeline

```
Phase 1: Plan
├── Extract target antigens (regex against 25 known antigens)
├── Identify development stages (keyword match)
├── Determine strategy (broad / targeted / comparative)
└── Decompose complex questions into sub-questions

Phase 2: Search
├── Build AgentQuery with plan parameters
└── Execute retrieve() via RAG engine

Phase 3: Evaluate
├── Sufficient:   ≥3 collections with hits + ≥10 total
├── Partial:      ≥2 collections + ≥5 total
└── Insufficient: try sub-questions for more evidence

Phase 4: Synthesize
└── LLM generates grounded answer with citations

Phase 5: Report
├── Query metadata (timestamp, collections, evidence count)
├── Analysis (the answer)
├── Evidence sources (by collection, top 5 each)
└── Knowledge graph items referenced
```

---

## 14. Data Ingest Pipelines

### BaseIngestPipeline Pattern

```python
class BaseIngestPipeline(ABC):
    def fetch(self, **kwargs) -> Any:          # Retrieve raw data
    def parse(self, raw_data) -> List[Model]:  # Validate into Pydantic models
    def embed_and_store(self, records, ...) -> int:  # Embed + insert to Milvus
    def run(self, **kwargs) -> int:            # Orchestrate full pipeline
```

### Five Ingest Pipelines

| Pipeline | Source | Records | Time | Collection |
|---|---|---|---|---|
| PubMed | NCBI E-utilities (esearch + efetch) | 5,047 | ~15 min | cart_literature |
| ClinicalTrials.gov | REST API v2 | 973 | ~3 min | cart_trials |
| FDA Constructs | Seed script (6 products) | 6 | ~5 sec | cart_constructs |
| Assays | Curated JSON (landmark papers) | 45 | ~30 sec | cart_assays |
| Manufacturing | Curated JSON (CMC data) | 30 | ~30 sec | cart_manufacturing |

Additional seed scripts populate: cart_safety (40), cart_biomarkers (43), cart_regulatory (25), cart_sequences (27), cart_realworld (30).

### Batch Processing Safety

```python
# Enum → string conversion for Milvus VARCHAR fields
for key, value in record_dict.items():
    if isinstance(value, Enum):
        record_dict[key] = value.value

# UTF-8 byte truncation for Milvus VARCHAR limits
for key, value in record_dict.items():
    if isinstance(value, str):
        encoded = value.encode("utf-8")
        if len(encoded) > 2990 and key in ("text_chunk", "text_summary"):
            record_dict[key] = encoded[:2990].decode("utf-8", errors="ignore")
```

---

## 15. FastAPI REST Server

### Service Configuration

- **Port:** 8522
- **Command:** `uvicorn api.main:app --host 0.0.0.0 --port 8522`

### Endpoints

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/health` | Service health + collection count + total vectors |
| `GET` | `/collections` | List all collections with record counts |
| `POST` | `/query` | Full RAG query (retrieve + Claude generation) |
| `POST` | `/search` | Evidence-only search (no LLM) |
| `POST` | `/find-related` | Cross-collection entity linking |
| `GET` | `/knowledge/stats` | Knowledge graph statistics |
| `GET` | `/metrics` | Prometheus-compatible metrics |

### Query Request Schema

```json
{
    "question": "Why do CD19 CAR-T therapies fail?",
    "target_antigen": "CD19",
    "collections": ["cart_literature", "cart_trials"],
    "year_min": 2015,
    "year_max": 2025
}
```

### Query Response Schema

```json
{
    "question": "...",
    "answer": "...",
    "evidence": [
        {
            "collection": "Literature",
            "id": "12345678",
            "score": 0.89,
            "text": "...",
            "metadata": {"relevance": "high"}
        }
    ],
    "knowledge_context": "...",
    "collections_searched": 11,
    "search_time_ms": 234.5
}
```

---

## 16. Streamlit Chat UI

### Service Configuration

- **Port:** 8521
- **Command:** `streamlit run app/cart_ui.py --server.port 8521`
- **Theme:** NVIDIA dark (black/green)

### Components

**Sidebar:**
- Collection statistics (10 owned + 1 read-only genomic)
- Deep Research mode toggle (activates agent reasoning pipeline)
- Collection filtering checkboxes
- Year range slider (1990–2030)
- Stage filter dropdown
- Target antigen multi-select

**Main Panel:**
- Chat message history (conversation memory: last 3 exchanges)
- Query input with demo query buttons
- Streaming response display
- Evidence panel (grouped by collection, color-coded badges)
- Citation relevance tags (high/medium/low)
- Clickable PubMed and ClinicalTrials.gov links

**Comparative Mode UI:**
- Entity A header (blue) — evidence cards for entity A
- "— VS —" divider (green)
- Entity B header (purple) — evidence cards for entity B

**Advanced Features:**
- Image upload for claim verification
- Knowledge graph visualization (PyVis network)
- Export buttons (PDF, Markdown, JSON)

---

## 17. Report Export

### Supported Formats

| Format | Method | Output |
|---|---|---|
| **Markdown** | `export_markdown()` | Structured report with headers, evidence table, citations |
| **JSON** | `export_json()` | Machine-readable export with all metadata |
| **PDF** | `export_pdf()` | Formatted report via ReportLab |

### Report Structure

```
1. Query Metadata (timestamp, collections searched, evidence count)
2. Analysis (LLM-generated answer)
3. Evidence Sources (by collection, top 5 per collection)
4. Knowledge Graph Items Referenced
5. Citation Links (PubMed, ClinicalTrials.gov)
```

---

## 18. Monitoring & Metrics

### Prometheus Metrics

```
cart_api_requests_total              # Total API requests (counter)
cart_collection_vectors{collection}  # Vectors per collection (gauge)
cart_query_duration_seconds          # Query latency histogram
cart_search_duration_seconds         # Search latency histogram
```

### Endpoints

- **Metrics:** `GET /metrics` (Prometheus scrape format)
- **Health:** `GET /health` (JSON)

---

## 19. Testing Strategy

### End-to-End Validation (`validate_e2e.py`)

5 tests:
1. Collection stats — verify all 11 collections have expected record counts
2. Single-collection search — verify vector similarity works
3. Multi-collection `search_all()` — verify parallel search returns results from multiple collections
4. Filtered search — verify `target_antigen == "CD19"` filter
5. Demo queries — run all 5 validated demo queries

### RAG Integration Test (`test_rag_pipeline.py`)

Tests full pipeline: embed → search_all → knowledge graph → Claude LLM response generation. Validates both synchronous and streaming modes.

### Demo Queries

```
1. "Why do CD19 CAR-T therapies fail in relapsed B-ALL?"
2. "Compare 4-1BB vs CD28 costimulatory domains for DLBCL"
3. "What manufacturing parameters predict clinical response?"
4. "BCMA CAR-T resistance mechanisms in multiple myeloma"
5. "How does T-cell exhaustion affect CAR-T persistence?"
```

---

## 20. HCLS AI Factory Integration

### Shared Infrastructure

| Resource | Sharing |
|---|---|
| **Milvus 2.4** | Shared instance — CAR-T adds 10 owned collections alongside `genomic_evidence` |
| **BGE-small-en-v1.5** | Same embedding model as RAG/Chat pipeline |
| **ANTHROPIC_API_KEY** | Loaded from `rag-chat-pipeline/.env` if not set |
| **DGX Spark** | All services on single GB10 hardware ($3,999) |

### Architectural Insight

The platform is not disease-specific. By changing the knowledge graph, query expansion maps, and collection schemas, the same RAG architecture serves any therapeutic area. The CAR-T agent proves this with a completely different domain (cell therapy manufacturing vs. small molecule drug discovery) running on the same infrastructure.

---

## 21. Implementation Sequence

### Quick Start

```bash
# 1. Clone
git clone https://github.com/ajones1923/cart-intelligence-agent.git
cd cart-intelligence-agent

# 2. Install
pip install -r requirements.txt

# 3. Configure
export ANTHROPIC_API_KEY=sk-ant-...

# 4. Create collections + seed FDA constructs
python3 scripts/setup_collections.py --seed-constructs

# 5. Ingest PubMed (~15 min)
python3 scripts/ingest_pubmed.py --max-results 5000

# 6. Ingest ClinicalTrials.gov (~3 min)
python3 scripts/ingest_clinical_trials.py --max-results 1500

# 7. Seed assays + manufacturing (~1 min)
python3 scripts/seed_assays.py
python3 scripts/seed_manufacturing.py

# 8. Validate
python3 scripts/validate_e2e.py

# 9. Launch UI
streamlit run app/cart_ui.py --server.port 8521

# 10. Launch API (separate terminal)
uvicorn api.main:app --host 0.0.0.0 --port 8522
```

### Docker Quick Start

```bash
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY
docker compose up -d
# UI: http://localhost:8521
# API: http://localhost:8522/docs
```

---

## Dependencies

| Package | Version | Purpose |
|---|---|---|
| pydantic | ≥2.0 | Data validation |
| pydantic-settings | ≥2.0 | Environment config |
| pymilvus | ≥2.4.0 | Vector DB client |
| sentence-transformers | ≥2.2.0 | BGE embeddings |
| anthropic | ≥0.18.0 | Claude API |
| streamlit | ≥1.30.0 | Web UI |
| fastapi | ≥0.109.0 | REST API |
| uvicorn[standard] | ≥0.27.0 | ASGI server |
| lxml | ≥5.0.0 | PubMed XML parsing |
| biopython | ≥1.83 | NCBI E-utilities |
| reportlab | ≥4.0.0 | PDF generation |
| pyvis | ≥0.3.0 | Network visualization |
| prometheus-client | ≥0.20.0 | Metrics export |
| apscheduler | ≥3.10.0 | Job scheduling |

---

*End of Project Bible*
