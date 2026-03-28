# Precision Autoimmune Intelligence Agent -- Architecture Guide

Architecture design document for the Precision Autoimmune Intelligence Agent,
part of the HCLS AI Factory on NVIDIA DGX Spark.

Author: Adam Jones
Date: March 2026

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [VAST AI OS Component Mapping](#vast-ai-os-component-mapping)
3. [System Architecture](#system-architecture)
4. [Milvus Collections](#milvus-collections)
5. [Clinical Engines](#clinical-engines)
6. [RAG Engine](#rag-engine)
7. [Export Pipeline](#export-pipeline)
8. [Infrastructure](#infrastructure)
9. [Demo Patients](#demo-patients)
10. [File Structure](#file-structure)

---

## Executive Summary

The Precision Autoimmune Intelligence Agent is a multi-collection
retrieval-augmented generation (RAG) system that provides clinical
decision support for autoimmune disease analysis. It combines domain-specific
knowledge bases with vector search and LLM synthesis to deliver
evidence-based recommendations.

### Key Results

| Metric | Value |
|--------|-------|
| Test suite | 455 tests across 8 files |
| Milvus collections | 14 domain-specific collections |
| Supported diseases | 13 autoimmune conditions |
| Biologic therapies | 22 drugs with PGx considerations |
| HLA alleles mapped | 22 alleles with disease associations |
| Autoantibodies mapped | 24 with sensitivity/specificity data |
| Disease activity scores | 20 scoring systems across all 13 diseases |
| Flare prediction patterns | 13 disease-specific biomarker patterns |
| Classification criteria | 10 ACR/EULAR and diagnostic criteria sets |
| Overlap syndromes | 9 cross-disease patterns |
| Lab test patterns | 45 extractable lab values |
| Demo patients | 9 with full clinical PDF datasets |
| API endpoints | 14 REST endpoints |
| Knowledge base version | 2.0.0 (last updated 2026-03-10) |

---

## VAST AI OS Component Mapping

The Autoimmune Agent maps to the VAST AI OS architecture as follows:

| VAST AI OS Layer | Agent Component |
|-----------------|-----------------|
| Data Layer | Milvus vector store (14 collections), demo_data/ PDFs |
| Model Layer | BAAI/bge-small-en-v1.5 (embedding), Claude (synthesis) |
| Inference Layer | FastAPI server (8532), RAG engine, clinical engines |
| Application Layer | Streamlit UI (8531), 10-tab interface |
| Orchestration | Docker Compose, health checks, CORS middleware |
| Integration | Cross-agent stubs (Biomarker, Imaging), event bus |

---

## System Architecture

```
+-------------------------------------------------------------------+
|                    Streamlit UI (:8531)                            |
|  +-------+ +--------+ +--------+ +--------+ +--------+ +------+  |
|  |Clinical| |Patient | |Document| |Diagnstc| |AutoAb  | | HLA  |  |
|  | Query  | |Analysis| |Ingest  | |Odyssey | |Panel   | |Anlys |  |
|  +-------+ +--------+ +--------+ +--------+ +--------+ +------+  |
|  +-------+ +--------+ +--------+ +--------+                      |
|  |Disease | |Flare   | |Therapy | |Knwldge |                      |
|  |Activty | |Predict | |Advisor | |Base    |                      |
|  +-------+ +--------+ +--------+ +--------+                      |
+-------------------------------------------------------------------+
        |                                         |
        v                                         v
+-------------------+                 +----------------------+
|  FastAPI (:8532)  |                 |  AutoimmuneAgent     |
|  14 endpoints     |                 |  analyze_patient()   |
|  Auth / CORS      |                 |  interpret_autoab()  |
|  Request limiter  |                 |  analyze_hla()       |
+-------------------+                 |  calc_activity()     |
        |                             |  predict_flares()    |
        v                             |  recommend_biologics |
+-------------------+                 +----------------------+
| AutoimmuneRAG     |                         |
| Engine            |                         v
|  retrieve()       |                 +----------------------+
|  query()          |                 | DiagnosticEngine     |
|  query_stream()   |                 |  eval_criteria()     |
|  find_related()   |                 |  gen_differential()  |
|  _embed()         |                 |  analyze_odyssey()   |
|  _build_knowledge |                 |  detect_overlap()    |
+-------------------+                 +----------------------+
        |                                     |
        v                                     v
+-------------------+                 +----------------------+
| CollectionManager |                 | TimelineBuilder      |
|  connect()        |                 |  extract_events()    |
|  search_all()     |<--------------->|  build_timeline()    |
|  insert_batch()   |                 |  classify_event()    |
+-------------------+                 +----------------------+
        |                                     |
        v                                     v
+-------------------+                 +----------------------+
| Milvus (:19530)   |                 | DocumentProcessor    |
| 14 collections    |                 |  process_pdf()       |
| IVF_FLAT / COSINE |                 |  chunk_text()        |
| 384-dim vectors   |                 |  classify_doc_type() |
+-------------------+                 |  extract_entities()  |
        |                             +----------------------+
        v                                     |
+-------------------+                         v
| BGE-small-en-v1.5 |                 +----------------------+
| Embedding Model   |                 | AutoimmuneExporter   |
| 384 dimensions    |                 |  to_markdown()       |
+-------------------+                 |  to_pdf_bytes()      |
                                      |  to_fhir_r4()        |
                                      +----------------------+
```

---

## Milvus Collections

The agent manages 14 specialized vector collections. All use COSINE similarity
with IVF_FLAT indexing (nlist=1024, nprobe=16) and 384-dimensional vectors
from BGE-small-en-v1.5.

| # | Collection Name | Description | Key Fields |
|---|----------------|-------------|------------|
| 1 | autoimmune_clinical_documents | Ingested patient records (PDFs) | patient_id, doc_type, specialty, provider, visit_date, source_file, page_number |
| 2 | autoimmune_patient_labs | Lab results with flag analysis | patient_id, test_name, value, unit, reference_range, flag, collection_date |
| 3 | autoimmune_autoantibody_panels | Autoantibody reference panels | antibody_name, associated_diseases, sensitivity, specificity, pattern, clinical_significance |
| 4 | autoimmune_hla_associations | HLA allele to disease risk mapping | allele, disease, odds_ratio, population, pmid, mechanism |
| 5 | autoimmune_disease_criteria | ACR/EULAR classification criteria | disease, criteria_set, criteria_type, year, required_score, criteria_items |
| 6 | autoimmune_disease_activity | Activity scoring reference | score_name, disease, components, thresholds, interpretation |
| 7 | autoimmune_flare_patterns | Flare prediction biomarker patterns | disease, biomarker_pattern, early_warning_signs, protective_factors |
| 8 | autoimmune_biologic_therapies | Biologic drug database with PGx | drug_name, drug_class, mechanism, indicated_diseases, pgx_considerations, contraindications |
| 9 | autoimmune_pgx_rules | Pharmacogenomic dosing rules | gene, variant, drug, phenotype, recommendation, evidence_level |
| 10 | autoimmune_clinical_trials | Autoimmune clinical trials | title, nct_id, phase, status, disease, intervention, enrollment |
| 11 | autoimmune_literature | Published literature | title, authors, journal, year, pmid, disease_focus, abstract_summary |
| 12 | autoimmune_patient_timelines | Patient diagnostic timelines | patient_id, event_type, event_date, description, provider, specialty, days_from_first_symptom |
| 13 | autoimmune_cross_disease | Cross-disease overlap syndromes | primary_disease, associated_conditions, shared_pathways, shared_biomarkers, co_occurrence_rate |
| 14 | genomic_evidence | Shared read-only collection | (managed by HCLS AI Factory core, read-only access) |

---

## Clinical Engines

The agent has 6 specialized engines that collaborate to deliver clinical
decision support.

### 1. AutoimmuneAgent (src/agent.py)

The main orchestrator. Runs the 5-step `analyze_patient` pipeline:

1. Interpret autoantibody panel against 24 antibody-disease associations
2. Analyze HLA profile against 22 allele-disease associations
3. Calculate disease activity scores using 20 scoring systems
4. Predict flare risk using 13 disease-specific biomarker patterns
5. Recommend biologic therapies from the 22-drug database with PGx filtering

Also provides cross-agent integration stubs for the Biomarker Agent
(inflammation context) and Imaging Agent (joint/organ assessment).

### 2. DiagnosticEngine (src/diagnostic_engine.py)

Clinical reasoning engine implementing:

- **Classification criteria evaluation:** 10 ACR/EULAR and diagnostic criteria
  sets with point-based scoring (e.g., 2010 ACR/EULAR RA requires >=6 points,
  2019 ACR/EULAR SLE requires >=10 points + ANA).
- **Differential diagnosis generation:** Ranks diseases by combining autoantibody
  specificity scores and HLA odds ratios (log2-scaled).
- **Diagnostic odyssey analysis:** Calculates time from first symptom to diagnosis,
  tracks specialists seen, identifies misdiagnoses and turning points.
- **Overlap syndrome detection:** 9 defined cross-disease patterns including MCTD,
  POTS/hEDS/MCAS triad, T1D-celiac overlap, and lupus-APS overlap.

### 3. AutoimmuneRAGEngine (src/rag_engine.py)

Multi-collection RAG engine providing:

- **14-collection parallel search** with ThreadPoolExecutor (max_workers=6)
- **Weighted scoring** using configurable per-collection weights (sum ~1.0)
- **Disease area detection** from query keywords (13 disease areas)
- **Knowledge augmentation** from the built-in knowledge base
- **Streaming and non-streaming synthesis** via Claude
- **Thread-safe conversation memory** (deque with configurable size)
- **Embedding cache** (256 entries, FIFO eviction)
- **Input sanitization** for Milvus filter expressions

### 4. DocumentProcessor (src/document_processor.py)

PDF ingestion pipeline:

- Text extraction via PyPDF2 (page-by-page)
- Sentence-boundary chunking (max 2500 chars, 200 char overlap)
- Document type classification (7 types: lab report, progress note, imaging,
  pathology, genetic, referral, medication list)
- Medical specialty detection (11 specialties)
- Entity extraction: 29 autoantibody names, 45 lab test patterns
- Provider and date extraction from text

### 5. TimelineBuilder (src/timeline_builder.py)

Diagnostic odyssey construction:

- Event classification from text using 12 event type patterns
  (symptom onset, diagnosis, misdiagnosis, lab result, imaging, biopsy,
  genetic test, treatment start, treatment change, flare, referral, ER visit)
- Date extraction from 4 date formats
- Chronological timeline building with days-from-first-symptom calculation
- Milvus-ready record generation for the patient timelines collection

### 6. AutoimmuneExporter (src/export.py)

Clinical report generation in 3 formats:

- **Markdown:** Structured report with tables for scores, HLA associations,
  and therapy recommendations
- **PDF:** Styled report using ReportLab with NVIDIA-themed colors
- **FHIR R4:** DiagnosticReport bundle with Patient, Observation, and
  DiagnosticReport resources

---

## RAG Engine

### Search Flow

```
Query --> Disease Area Detection --> Embed (BGE-small + instruction prefix)
  --> Parallel Search (14 collections, top_k=5 each, score_threshold=0.40)
  --> Weighted Scoring (per-collection weights)
  --> Deduplication (by ID and text content hash)
  --> Knowledge Augmentation (HLA, autoantibody, therapy, flare patterns)
  --> LLM Synthesis (Claude, system prompt, conversation memory)
  --> Streaming Response with Evidence Citations
```

### Collection Weights

| Collection | Weight | Label |
|-----------|--------|-------|
| autoimmune_clinical_documents | 0.18 | Clinical Document |
| autoimmune_patient_labs | 0.14 | Lab Result |
| autoimmune_autoantibody_panels | 0.12 | Autoantibody |
| autoimmune_hla_associations | 0.08 | HLA Association |
| autoimmune_disease_criteria | 0.08 | Classification Criteria |
| autoimmune_disease_activity | 0.07 | Disease Activity |
| autoimmune_flare_patterns | 0.06 | Flare Pattern |
| autoimmune_biologic_therapies | 0.06 | Biologic Therapy |
| autoimmune_clinical_trials | 0.05 | Clinical Trial |
| autoimmune_literature | 0.05 | Literature |
| autoimmune_pgx_rules | 0.04 | PGx Rule |
| autoimmune_patient_timelines | 0.03 | Timeline |
| autoimmune_cross_disease | 0.02 | Cross-Disease |
| genomic_evidence | 0.02 | Genomic Evidence |
| **Total** | **1.00** | |

### Citation Scoring

| Score Range | Relevance Level |
|------------|----------------|
| >= 0.80 | High |
| >= 0.60 | Medium |
| < 0.60 | Low |

### Embedding Configuration

- Model: BAAI/bge-small-en-v1.5
- Dimensions: 384
- Batch size: 32
- Instruction prefix: "Represent this sentence for searching relevant passages: "
- Asymmetric search: queries use the instruction prefix, documents do not

---

## Export Pipeline

### Markdown Export

Generates a structured clinical report with:
- Critical alerts section
- Disease activity score table
- Flare risk predictions with contributing/protective factors
- HLA-disease association table
- Biologic therapy recommendations with PGx and monitoring
- Evidence sources with relevance badges

### PDF Export

Uses ReportLab to generate styled PDF reports with:
- NVIDIA-themed colors (green #76B900 headers)
- Tabular disease activity scores
- Biologic therapy details with PGx considerations
- Clinical footer with knowledge base version

### FHIR R4 Export

Generates a FHIR R4 Bundle (type: collection) containing:
- Patient resource with identifier
- Observation resources for disease activity scores and flare risk predictions
- DiagnosticReport resource linking all observations
- HL7 FHIR StructureDefinition profile compliance

---

## Infrastructure

### Ports

| Service | Port | Protocol |
|---------|------|----------|
| Streamlit UI | 8531 | HTTP |
| FastAPI REST API | 8532 | HTTP |
| Milvus (shared) | 19530 | gRPC |

### Docker Compose Services

| Service | Container Name | Purpose |
|---------|---------------|---------|
| autoimmune-streamlit | autoimmune-streamlit | Streamlit 10-tab UI |
| autoimmune-api | autoimmune-api | FastAPI REST server (2 workers) |
| autoimmune-setup | autoimmune-setup | One-shot collection creation and knowledge seeding |

All services connect to the shared `hcls-network` and depend on the Milvus
stack from the main `docker-compose.dgx-spark.yml`.

Health checks:
- Streamlit: `/_stcore/health` every 30s (start period 60s)
- API: `/healthz` every 30s (start period 30s)

---

## Demo Patients

| # | Name | Conditions | Key Features |
|---|------|-----------|-------------|
| 1 | Sarah Mitchell, 34F | SLE (lupus nephritis) | ANA 1:640, anti-dsDNA, anti-Smith, low complement |
| 2 | Maya Rodriguez, 28F | POTS/hEDS/MCAS | Dysautonomia diagnostic odyssey |
| 3 | Linda Chen, 45F | Sjogren's | Anti-SSA/Ro+, anti-SSB/La+, Schirmer test, sicca syndrome |
| 4 | David Park, 45M | AS | HLA-B*27:05, uveitis, sacroiliitis, 3-year diagnostic odyssey |
| 5 | Rachel Thompson, 38F | MCTD | Mixed connective tissue disease, anti-U1 RNP |
| 6 | Emma Williams, 34F | MS (RRMS) | Optic neuritis, MRI lesions, oligoclonal bands, HLA-DRB1*15:01 |
| 7 | James Cooper, 19M | T1D + Celiac | Anti-GAD65, anti-IA-2, anti-ZnT8, tTG-IgA+, HLA-DQ2/DQ8 |
| 8 | Karen Foster, 48F | SSc (dcSSc) | Anti-Scl-70+, Raynaud's, ILD, modified Rodnan skin score |
| 9 | Michael Torres, 41M | Graves' Disease | Thyroid-stimulating immunoglobulin, Burch-Wartofsky scoring |

Each patient has a directory under `demo_data/` containing clinical PDFs
(progress notes, lab reports, imaging, pathology, genetic tests, referral
letters, and medication lists).

---

## File Structure

```
precision_autoimmune_agent/
|-- api/
|   |-- __init__.py
|   |-- main.py                  # FastAPI server (14 endpoints)
|
|-- app/
|   |-- autoimmune_ui.py         # Streamlit 10-tab UI
|
|-- config/
|   |-- settings.py              # Pydantic BaseSettings (AUTO_ prefix)
|   |-- logging.py               # Centralized logging config
|
|-- data/
|   |-- cache/                   # Embedding cache
|   |-- reference/               # Reference data
|
|-- demo_data/
|   |-- sarah_mitchell/          # SLE (lupus nephritis) patient PDFs
|   |-- maya_rodriguez/          # POTS/hEDS/MCAS patient PDFs
|   |-- linda_chen/              # Sjogren's patient PDFs
|   |-- david_park/              # AS patient PDFs
|   |-- rachel_thompson/         # MCTD patient PDFs
|   |-- emma_williams/           # MS (RRMS) patient PDFs
|   |-- james_cooper/            # T1D + Celiac patient PDFs
|   |-- karen_foster/            # SSc (dcSSc) patient PDFs
|   |-- michael_torres/          # Graves' Disease patient PDFs
|
|-- docs/
|   |-- PRECISION_AUTOIMMUNE_AGENT_RESEARCH_PAPER.md
|   |-- API_REFERENCE.md
|   |-- DEMO_GUIDE.md
|   |-- ARCHITECTURE_GUIDE.md    # This document
|   |-- PROJECT_BIBLE.md
|
|-- src/
|   |-- __init__.py
|   |-- agent.py                 # AutoimmuneAgent orchestrator
|   |-- collections.py           # Milvus collection manager (14 schemas)
|   |-- diagnostic_engine.py     # Classification criteria, differential, odyssey
|   |-- document_processor.py    # PDF ingestion pipeline
|   |-- export.py                # Markdown / PDF / FHIR R4 exporter
|   |-- knowledge.py             # Knowledge base (HLA, antibodies, therapies, flare)
|   |-- models.py                # Pydantic data models
|   |-- rag_engine.py            # Multi-collection RAG engine
|   |-- timeline_builder.py      # Diagnostic odyssey timeline builder
|
|-- tests/
|   |-- test_api.py              # 57 tests
|   |-- test_autoimmune.py       # 31 tests
|   |-- test_collections.py      # 57 tests
|   |-- test_diagnostic_engine.py # 51 tests
|   |-- test_export.py           # 34 tests
|   |-- test_production_readiness.py # 59 tests
|   |-- test_rag_engine.py       # 93 tests
|   |-- test_timeline_builder.py # 49 tests
|
|-- docker-compose.yml           # 3 services: streamlit, api, setup
|-- Dockerfile
|-- pyproject.toml
|-- requirements.txt             # 19 dependencies
|-- run.sh                       # Startup script
|-- README.md
```

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Autoimmune Agent is a clinical decision support research tool for autoimmune disease assessment. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
