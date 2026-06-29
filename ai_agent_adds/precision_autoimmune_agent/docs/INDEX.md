# Precision Autoimmune Intelligence Agent

Autoimmune disease intelligence with autoantibody interpretation, HLA association analysis, disease activity scoring, flare prediction, and biologic therapy PGx recommendations. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

## Overview

The Precision Autoimmune Intelligence Agent transforms complex autoimmune clinical data into actionable diagnostic and therapeutic intelligence. It combines multi-collection RAG search across 14 Milvus collections with 6 deterministic clinical analysis engines to interpret autoantibody panels, evaluate HLA-disease associations, score disease activity across 20 validated instruments, predict flares from biomarker patterns, recommend biologic therapies with pharmacogenomic context, and analyze diagnostic odyssey timelines. The agent supports 13 autoimmune diseases and ingests patient clinical documents (PDFs) for longitudinal analysis.

| Collection | Records | Content |
|---|---|---|
| **Clinical Documents** | variable | Ingested patient records (PDFs) -- progress notes, labs, imaging, pathology |
| **Patient Labs** | variable | Lab results with flag analysis and reference range comparison |
| **Autoantibody Panels** | 24 | Autoantibody reference panels with disease associations, sensitivity, specificity |
| **HLA Associations** | 22 | HLA allele to disease risk mapping with odds ratios and PMIDs |
| **Disease Criteria** | 13 | ACR/EULAR classification criteria for 13 diseases |
| **Disease Activity** | 20 | Activity scoring reference (DAS28, SLEDAI, BASDAI, CDAI, etc.) |
| **Flare Patterns** | 13 | Flare prediction biomarker patterns and early warning signs |
| **Biologic Therapies** | 22 | Biologic drug database with PGx considerations |
| **PGx Rules** | variable | Pharmacogenomic dosing rules for autoimmune therapies |
| **Clinical Trials** | variable | Autoimmune disease clinical trials |
| **Literature** | variable | Published autoimmune literature and research |
| **Patient Timelines** | variable | Patient diagnostic timeline events for odyssey analysis |
| **Cross-Disease** | 9 | Cross-disease overlap syndromes and shared pathogenic mechanisms |
| **Genomic Evidence** | 30 | *(read-only)* Shared from Stage 2 RAG pipeline |

### 6 Clinical Analysis Engines

| Engine | Function |
|---|---|
| **AutoantibodyInterpreter** | Autoantibody panel interpretation with disease association scoring, titer analysis, and staining pattern recognition (ANA, anti-dsDNA, RF, anti-CCP, etc.) |
| **HLAAssociationAnalyzer** | HLA allele to disease risk mapping across 22 alleles with odds ratios, population context, and mechanistic annotations |
| **DiseaseActivityScorer** | 20 validated disease activity scoring systems across 13 diseases (DAS28-CRP, SLEDAI-2K, BASDAI, CDAI, EDSS, etc.) |
| **FlarePredictor** | Biomarker pattern-based flare prediction with risk stratification (low/moderate/high/imminent), contributing and protective factors |
| **BiologicTherapyAdvisor** | Biologic therapy recommendation with PGx filtering across 22 biologics (TNF inhibitors, IL-6 inhibitors, B-cell depleters, JAK inhibitors, etc.) |
| **DiagnosticOdysseyAnalyzer** | Diagnostic timeline analysis from ingested clinical documents -- tracks symptom onset to diagnosis journey, identifies delays and patterns |

### Example Queries

```
"Interpret ANA 1:640 homogeneous pattern with positive anti-dsDNA for this SLE patient"
"What is my flare risk given rising anti-dsDNA and falling C3/C4?"
"HLA-B27 positive with inflammatory back pain -- ankylosing spondylitis risk?"
"Compare rituximab vs belimumab for refractory lupus nephritis"
"Why did it take 4 years from first symptoms to my RA diagnosis?"
```

### Demo Guide

For a complete walkthrough of all 10 UI tabs with the 9 demo patients, see the **[Demo Guide](demo-guide.md)**.

## Architecture

```
Patient Clinical Documents (PDFs)
    |
    v
[Document Processor] ── Chunk + Embed ──> Milvus Collections
    |
    v
Patient Data Input (Autoantibodies + HLA + Labs + Genotypes)
    |
    v
[6 Clinical Analysis Engines -- Parallel Execution]
    |               |              |              |
    v               v              v              v
Autoantibody    HLA             Disease        Flare
Interpreter     Association     Activity       Predictor
(24 antibodies) Analyzer        Scorer         (13 patterns)
                (22 alleles)    (20 scores)
    |               |              |              |
    |               v              |              |
    |           Biologic        Diagnostic        |
    |           Therapy         Odyssey           |
    |           Advisor         Analyzer          |
    |           (22 drugs)      (timelines)       |
    +-------+-------+--------------+--------------+
            |
            v
    [Multi-Collection RAG Engine]
    Parallel search across 14 Milvus collections
    (ThreadPoolExecutor, configurable weights)
            |
            v
    [Claude Sonnet 4.6] -> Grounded response with citations
            |
            v
    [Export: FHIR R4 | PDF | Markdown]
```

Built on the HCLS AI Factory platform:

- **Vector DB:** Milvus 2.4 with IVF_FLAT/COSINE indexes (nlist=1024, nprobe=16)
- **Embeddings:** BGE-small-en-v1.5 (384-dim)
- **LLM:** Claude Sonnet 4.6 (Anthropic API)
- **UI:** Streamlit 10 tabs (port 8531) | **API:** FastAPI (port 8532)
- **Hardware target:** NVIDIA DGX Spark ($4,699)

### UI Tabs

| # | Tab | Content |
|---|---|---|
| 1 | Clinical Query | RAG-powered clinical questions with multi-collection search and streaming |
| 2 | Patient Analysis | Full autoimmune analysis pipeline (autoantibodies + HLA + activity + flare + therapy) |
| 3 | Document Ingest | PDF upload and ingestion for clinical documents (progress notes, labs, imaging) |
| 4 | Diagnostic Odyssey | Timeline visualization of symptom onset to diagnosis journey |
| 5 | Autoantibody Panel | Autoantibody panel interpretation with disease association scoring |
| 6 | HLA Analysis | HLA allele to disease risk mapping with odds ratios |
| 7 | Disease Activity | Disease activity scoring with 20 validated instruments |
| 8 | Flare Prediction | Flare risk assessment with biomarker pattern analysis |
| 9 | Therapy Advisor | Biologic therapy recommendations with PGx context |
| 10 | Knowledge Base | Knowledge base explorer with collection statistics and version info |

## Setup

### Prerequisites

- Python 3.10+
- Milvus 2.4 running on `localhost:19530`
- `AUTO_ANTHROPIC_API_KEY` environment variable (or `ANTHROPIC_API_KEY` in `.env`)

### Install

```bash
cd ai_agent_adds/precision_autoimmune_agent
pip install -r requirements.txt
```

### 1. Create Collections and Seed Knowledge

```bash
python3 scripts/setup_collections.py --seed
```

Creates 14 Milvus collections with IVF_FLAT indexes and seeds knowledge base v2.0.0.

### 2. Ingest Demo Patient Data

```bash
python3 scripts/generate_demo_patients.py
```

Generates and ingests 9 demo patient datasets with clinical documents, lab results, autoantibody panels, and diagnostic timelines.

### 3. Run Unit Tests

```bash
python3 -m pytest tests/ -v
```

455 tests covering all modules: collections (57), RAG engine (93), API (57), production readiness (59), diagnostic engine (51), timeline builder (49), export (34), and core agent (31).

### 4. Launch UI

```bash
streamlit run app/autoimmune_ui.py --server.port 8531
```

### 5. Launch API (separate terminal)

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8532
```

## Project Structure

```
precision_autoimmune_agent/
├── src/
│   ├── models.py                  # Pydantic data models (autoantibodies, HLA, activity scores)
│   ├── collections.py             # 14 Milvus collection schemas + AutoimmuneCollectionManager
│   ├── rag_engine.py              # Multi-collection RAG engine + Claude
│   ├── agent.py                   # AutoimmuneAgent orchestrator (6 engines)
│   ├── knowledge.py               # Knowledge base v2.0.0 (HLA, autoantibodies, biologics, flare patterns)
│   ├── diagnostic_engine.py       # ACR/EULAR criteria, diagnostic odyssey, differential diagnosis
│   ├── document_processor.py      # PDF ingestion, chunking, embedding
│   ├── timeline_builder.py        # Diagnostic timeline construction from clinical documents
│   └── export.py                  # FHIR R4, PDF, Markdown export
├── app/
│   └── autoimmune_ui.py           # Streamlit UI (10 tabs)
├── api/
│   └── main.py                    # FastAPI REST server (14 endpoints)
├── config/
│   ├── settings.py                # AutoimmuneSettings (AUTO_ prefix, 14 collection weights)
│   └── logging.py                 # Centralized logging configuration
├── demo_data/
│   ├── sarah_mitchell/            # SLE -- 27 clinical PDFs
│   ├── linda_chen/                # Sjogren's
│   ├── maya_rodriguez/            # POTS / hEDS / MCAS (dysautonomia)
│   ├── david_park/                # AS
│   ├── rachel_thompson/           # MCTD
│   ├── emma_williams/             # MS
│   ├── james_cooper/              # T1D + Celiac overlap
│   ├── karen_foster/              # SSc
│   └── michael_torres/            # Graves
├── scripts/
│   ├── setup_collections.py       # Create and seed Milvus collections
│   └── generate_demo_patients.py  # Generate 9 demo patient datasets
├── tests/
│   ├── test_autoimmune.py         # 31 tests -- core agent
│   ├── test_collections.py        # 57 tests -- collection schemas and manager
│   ├── test_rag_engine.py         # 93 tests -- RAG engine
│   ├── test_api.py                # 57 tests -- FastAPI endpoints
│   ├── test_production_readiness.py # 59 tests -- production readiness
│   ├── test_diagnostic_engine.py  # 51 tests -- diagnostic engine
│   ├── test_timeline_builder.py   # 49 tests -- timeline builder
│   └── test_export.py            # 34 tests -- export formats
├── docker-compose.yml
├── Dockerfile
├── requirements.txt
└── LICENSE                        # Apache 2.0
```

**44 Python files | ~10,500 lines | Apache 2.0**

## Performance

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified memory):

| Metric | Value |
|---|---|
| Vector search (14 collections, top-5 each) | 10-18 ms (cached) |
| Disease activity scoring (all instruments) | <100 ms |
| Flare risk prediction | <50 ms |
| Autoantibody panel interpretation | <50 ms |
| HLA association analysis | <20 ms |
| Full RAG query (search + Claude) | ~20-30 sec |
| PDF document ingestion (per page) | ~2 sec |
| FHIR R4 export + validation | <500 ms |
| Unit tests (455 tests) | ~60 sec |

## Key Stats

| Metric | Value |
|---|---|
| Unit tests | **455** |
| Diseases covered | **13** (RA, SLE, MS, T1D, IBD, Psoriasis, AS, Sjogren's, SSc, MG, Celiac, Graves, Hashimoto) |
| Biologic therapies | **22** |
| Autoantibodies | **24** |
| HLA alleles | **22** |
| Activity scoring systems | **20** |
| Demo patients | **9** |
| Knowledge version | **v2.0.0** |

## Status

- **Phase 1 (Scaffold)** -- Complete. Architecture, data models, 14 collection schemas, 6 clinical engines, RAG engine, agent orchestrator, and Streamlit UI.
- **Phase 2 (Knowledge)** -- Complete. Knowledge base v2.0.0 with 22 HLA alleles, 24 autoantibodies, 22 biologics, 20 activity scores, 13 flare patterns, and 9 overlap syndromes.
- **Phase 3 (Integration)** -- Complete. Full RAG pipeline with Claude generating grounded, disease-specific interpretations. FHIR R4 export with structural validation.
- **Phase 4 (Testing)** -- Complete. 455 unit tests passing across 8 test files. Production readiness validation.
- **Phase 5 (Demo Ready)** -- Complete. Production-quality Streamlit UI with 10 tabs. 9 demo patients with clinical document PDFs. Document ingestion pipeline.

## Credits

- **Adam Jones**
- **Apache 2.0 License**
