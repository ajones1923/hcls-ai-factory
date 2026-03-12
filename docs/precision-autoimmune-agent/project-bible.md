# Precision Autoimmune Intelligence Agent -- Project Bible

Complete implementation reference for the Precision Autoimmune Intelligence
Agent, part of the HCLS AI Factory pipeline: Patient DNA -> Drug Candidates.

Author: Adam Jones
Date: March 2026

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Pipeline Pattern](#2-pipeline-pattern)
3. [DGX Spark Hardware](#3-dgx-spark-hardware)
4. [Repository Layout](#4-repository-layout)
5. [Docker Compose Services](#5-docker-compose-services)
6. [Milvus Collection Schemas](#6-milvus-collection-schemas)
7. [Pydantic Data Models](#7-pydantic-data-models)
8. [Configuration Reference](#8-configuration-reference)
9. [Embedding Strategy](#9-embedding-strategy)
10. [Autoantibody Engine](#10-autoantibody-engine)
11. [HLA Association Engine](#11-hla-association-engine)
12. [Disease Activity Engine](#12-disease-activity-engine)
13. [Flare Prediction Engine](#13-flare-prediction-engine)
14. [Biologic Therapy Engine](#14-biologic-therapy-engine)
15. [Diagnostic Engine](#15-diagnostic-engine)
16. [RAG Engine](#16-rag-engine)
17. [Agent Orchestrator](#17-agent-orchestrator)
18. [Export Pipeline](#18-export-pipeline)
19. [FastAPI REST Server](#19-fastapi-rest-server)
20. [Streamlit UI](#20-streamlit-ui)
21. [Demo Patients](#21-demo-patients)
22. [Cross-Agent Integration](#22-cross-agent-integration)
23. [Monitoring and Metrics](#23-monitoring-and-metrics)
24. [Testing](#24-testing)
25. [HCLS AI Factory Integration](#25-hcls-ai-factory-integration)
26. [Quick Start](#26-quick-start)
27. [Dependencies](#27-dependencies)

---

## 1. Project Overview

The Precision Autoimmune Intelligence Agent is a clinical decision-support
system specializing in autoimmune disease analysis. It provides:

- Autoantibody panel interpretation with sensitivity/specificity data
- HLA-disease association analysis with odds ratios
- Disease activity scoring across 20 validated scoring systems
- Biomarker-based flare prediction for 13 autoimmune conditions
- Biologic therapy recommendation with pharmacogenomic context
- Diagnostic odyssey analysis (timeline, misdiagnosis detection)
- RAG-powered clinical Q&A across 14 Milvus vector collections
- Clinical report export in Markdown, PDF, and FHIR R4 formats

The agent is deployed as part of the 5-agent HCLS AI Factory intelligence
layer alongside Biomarker, Imaging, Oncology, and CAR-T agents.

---

## 2. Pipeline Pattern

The autoimmune analysis follows a 5-step pipeline within the
`AutoimmuneAgent.analyze_patient()` method:

```
AutoimmunePatientProfile
  |
  |--> Step 1: Interpret autoantibody panel
  |      - Match positive antibodies against AUTOANTIBODY_DISEASE_MAP (24 entries)
  |      - Return disease associations with sensitivity/specificity
  |
  |--> Step 2: Analyze HLA associations
  |      - Match alleles against HLA_DISEASE_ASSOCIATIONS (22 alleles)
  |      - Support exact and broad allele matching (e.g., B*27:05 -> B*27)
  |      - Sort by odds ratio descending
  |
  |--> Step 3: Calculate disease activity scores
  |      - Select scoring system by diagnosed condition
  |      - Use CRP/ESR for simplified scoring
  |      - Assign level: remission/low/moderate/high/very_high
  |
  |--> Step 4: Predict flare risk
  |      - Base risk 0.3
  |      - Evaluate early warning biomarkers per disease pattern
  |      - Clamp score to [0.0, 1.0]
  |      - Configurable thresholds: imminent >=0.8, high >=0.6, moderate >=0.4
  |
  |--> Step 5: Recommend biologic therapies
  |      - Filter 22-drug database by indicated diseases
  |      - Include PGx considerations and contraindications
  |
  |--> Step 6: Generate critical alerts
         - High/very high disease activity
         - High/imminent flare risk
         - Strong HLA associations (OR > 5)

  ==> AutoimmuneAnalysisResult
```

---

## 3. DGX Spark Hardware

The HCLS AI Factory runs on NVIDIA DGX Spark:

| Component | Specification |
|-----------|--------------|
| GPU | NVIDIA Grace Hopper |
| CPU | ARM-based Grace CPU |
| Memory | Unified memory architecture |
| Storage | NVMe SSD |
| Network | High-speed interconnect |
| OS | Ubuntu / NVIDIA AI Enterprise |
| CUDA | 12.x |

The autoimmune agent uses GPU acceleration for:
- Embedding model inference (BGE-small-en-v1.5 via sentence-transformers)
- Milvus vector similarity search (IVF_FLAT with COSINE metric)

---

## 4. Repository Layout

```
precision_autoimmune_agent/
|
|-- api/
|   |-- __init__.py
|   |-- main.py                  # FastAPI server, 14 endpoints, lifespan mgmt
|
|-- app/
|   |-- autoimmune_ui.py         # Streamlit 10-tab interface
|
|-- config/
|   |-- settings.py              # AutoimmuneSettings (Pydantic BaseSettings)
|   |-- logging.py               # Loguru configuration
|
|-- data/
|   |-- cache/                   # Embedding and query cache
|   |-- reference/               # Reference datasets
|
|-- demo_data/
|   |-- sarah_mitchell/          # SLE (28F) -- 20+ clinical PDFs
|   |-- rachel_thompson/         # RA (52F)
|   |-- james_cooper/            # T1D + Celiac (19M)
|   |-- david_park/              # AS -- 29 clinical PDFs (3-year odyssey)
|   |-- emma_williams/           # MS
|   |-- linda_chen/              # Sjogren's
|   |-- maya_rodriguez/          # Psoriasis
|   |-- michael_torres/          # IBD
|   |-- karen_foster/            # SSc
|
|-- docs/
|   |-- PRECISION_AUTOIMMUNE_AGENT_RESEARCH_PAPER.md
|   |-- API_REFERENCE.md
|   |-- DEMO_GUIDE.md
|   |-- ARCHITECTURE_GUIDE.md
|   |-- PROJECT_BIBLE.md         # This document
|
|-- logs/                        # Runtime logs (loguru)
|
|-- scripts/
|   |-- setup_collections.py     # Collection creation and knowledge seeding
|
|-- src/
|   |-- __init__.py
|   |-- agent.py                 # AutoimmuneAgent (main orchestrator)
|   |-- collections.py           # AutoimmuneCollectionManager (14 schemas)
|   |-- diagnostic_engine.py     # DiagnosticEngine (criteria, differential, odyssey)
|   |-- document_processor.py    # DocumentProcessor (PDF ingestion)
|   |-- export.py                # AutoimmuneExporter (Markdown, PDF, FHIR R4)
|   |-- knowledge.py             # Knowledge base (v2.0.0)
|   |-- models.py                # Pydantic data models (13 models)
|   |-- rag_engine.py            # AutoimmuneRAGEngine (multi-collection RAG)
|   |-- timeline_builder.py      # TimelineBuilder (diagnostic odyssey)
|
|-- tests/
|   |-- __init__.py
|   |-- test_api.py              # 57 tests -- API endpoint testing
|   |-- test_autoimmune.py       # 31 tests -- Agent logic
|   |-- test_collections.py      # 57 tests -- Collection manager
|   |-- test_diagnostic_engine.py # 51 tests -- Diagnostic engine
|   |-- test_export.py           # 34 tests -- Export formats
|   |-- test_production_readiness.py # 59 tests -- Production checks
|   |-- test_rag_engine.py       # 93 tests -- RAG engine
|   |-- test_timeline_builder.py # 49 tests -- Timeline builder
|
|-- docker-compose.yml           # 3 services
|-- Dockerfile
|-- pyproject.toml
|-- requirements.txt             # 19 dependencies
|-- run.sh
|-- README.md
```

---

## 5. Docker Compose Services

Three services defined in `docker-compose.yml`:

### autoimmune-streamlit

```yaml
container: autoimmune-streamlit
port: 8531:8531
env:
  AUTO_MILVUS_HOST: milvus-standalone
  AUTO_MILVUS_PORT: 19530
  AUTO_STREAMLIT_PORT: 8531
  ANTHROPIC_API_KEY: ${ANTHROPIC_API_KEY}
volumes:
  - ./demo_data:/app/demo_data:ro
  - ./data:/app/data
healthcheck: /_stcore/health (30s interval, 60s start_period)
restart: unless-stopped
```

### autoimmune-api

```yaml
container: autoimmune-api
command: uvicorn api.main:app --host 0.0.0.0 --port 8532 --workers 2
port: 8532:8532
env:
  AUTO_MILVUS_HOST: milvus-standalone
  AUTO_MILVUS_PORT: 19530
  AUTO_API_PORT: 8532
  ANTHROPIC_API_KEY: ${ANTHROPIC_API_KEY}
healthcheck: /healthz (30s interval, 30s start_period)
restart: unless-stopped
```

### autoimmune-setup

```yaml
container: autoimmune-setup
command: python scripts/setup_collections.py --seed
restart: "no"
```

All services connect to the `hcls-network` (external Docker network) and
depend on the shared Milvus stack from `docker-compose.dgx-spark.yml`.

---

## 6. Milvus Collection Schemas

All collections share a common pattern: VARCHAR primary key `id`, FLOAT_VECTOR
`embedding` (384 dimensions), and VARCHAR `text_chunk` for the source text.
Indexing uses IVF_FLAT with nlist=1024, searched with nprobe=16 and COSINE
similarity.

### 1. autoimmune_clinical_documents

Ingested patient clinical records from PDF documents.

| Field | Type | Max Length | Description |
|-------|------|-----------|-------------|
| id | VARCHAR (PK) | 128 | Record identifier |
| embedding | FLOAT_VECTOR | 384 dim | BGE-small vector |
| text_chunk | VARCHAR | 3000 | Document text chunk |
| patient_id | VARCHAR | 64 | Patient identifier |
| doc_type | VARCHAR | 128 | Document type (progress_note, lab_report, imaging, pathology, genetic_report, referral_letter, medication_list) |
| specialty | VARCHAR | 128 | Medical specialty (rheumatology, neurology, etc.) |
| provider | VARCHAR | 256 | Provider name |
| visit_date | VARCHAR | 32 | Visit/collection date |
| source_file | VARCHAR | 512 | Source PDF filename |
| page_number | INT64 | -- | Page number in source PDF |
| chunk_index | INT64 | -- | Chunk index within page |

### 2. autoimmune_patient_labs

Laboratory results with reference range analysis.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| patient_id | VARCHAR | 64 |
| test_name | VARCHAR | 256 |
| value | FLOAT | -- |
| unit | VARCHAR | 64 |
| reference_range | VARCHAR | 128 |
| flag | VARCHAR | 32 |
| collection_date | VARCHAR | 32 |
| panel_name | VARCHAR | 256 |

### 3. autoimmune_autoantibody_panels

Autoantibody reference panels with disease associations.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| antibody_name | VARCHAR | 128 |
| associated_diseases | VARCHAR | 1024 |
| sensitivity | FLOAT | -- |
| specificity | FLOAT | -- |
| pattern | VARCHAR | 128 |
| clinical_significance | VARCHAR | 2000 |
| interpretation_guide | VARCHAR | 2000 |

### 4. autoimmune_hla_associations

HLA allele to disease risk mapping with odds ratios.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| allele | VARCHAR | 64 |
| disease | VARCHAR | 256 |
| odds_ratio | FLOAT | -- |
| population | VARCHAR | 128 |
| pmid | VARCHAR | 32 |
| mechanism | VARCHAR | 1024 |
| clinical_implication | VARCHAR | 2000 |

### 5. autoimmune_disease_criteria

ACR/EULAR classification and diagnostic criteria.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| disease | VARCHAR | 256 |
| criteria_set | VARCHAR | 256 |
| criteria_type | VARCHAR | 64 |
| year | INT64 | -- |
| required_score | VARCHAR | 128 |
| criteria_items | VARCHAR | 3000 |
| sensitivity_specificity | VARCHAR | 256 |

### 6. autoimmune_disease_activity

Disease activity scoring systems and interpretation.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| score_name | VARCHAR | 128 |
| disease | VARCHAR | 256 |
| components | VARCHAR | 2000 |
| thresholds | VARCHAR | 1024 |
| interpretation | VARCHAR | 2000 |
| monitoring_frequency | VARCHAR | 512 |

### 7. autoimmune_flare_patterns

Flare prediction biomarker patterns and early warning signs.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| disease | VARCHAR | 256 |
| biomarker_pattern | VARCHAR | 2000 |
| early_warning_signs | VARCHAR | 2000 |
| typical_timeline | VARCHAR | 512 |
| protective_factors | VARCHAR | 1024 |
| intervention_triggers | VARCHAR | 1024 |

### 8. autoimmune_biologic_therapies

Biologic therapy database with pharmacogenomic considerations.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| drug_name | VARCHAR | 128 |
| drug_class | VARCHAR | 128 |
| mechanism | VARCHAR | 512 |
| indicated_diseases | VARCHAR | 1024 |
| pgx_considerations | VARCHAR | 2000 |
| contraindications | VARCHAR | 1024 |
| monitoring | VARCHAR | 2000 |
| dosing | VARCHAR | 512 |
| evidence_level | VARCHAR | 64 |

### 9. autoimmune_pgx_rules

Pharmacogenomic dosing rules for autoimmune therapies.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| gene | VARCHAR | 64 |
| variant | VARCHAR | 128 |
| drug | VARCHAR | 128 |
| phenotype | VARCHAR | 128 |
| recommendation | VARCHAR | 2000 |
| evidence_level | VARCHAR | 64 |
| pmid | VARCHAR | 256 |

### 10. autoimmune_clinical_trials

Autoimmune disease clinical trials.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 256 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| title | VARCHAR | 1024 |
| nct_id | VARCHAR | 32 |
| phase | VARCHAR | 32 |
| status | VARCHAR | 64 |
| disease | VARCHAR | 256 |
| intervention | VARCHAR | 512 |
| biomarker_criteria | VARCHAR | 1024 |
| enrollment | INT64 | -- |
| start_year | INT64 | -- |
| sponsor | VARCHAR | 256 |

### 11. autoimmune_literature

Published autoimmune literature and research.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 256 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| title | VARCHAR | 1024 |
| authors | VARCHAR | 1024 |
| journal | VARCHAR | 256 |
| year | INT64 | -- |
| pmid | VARCHAR | 32 |
| disease_focus | VARCHAR | 256 |
| keywords | VARCHAR | 1024 |
| abstract_summary | VARCHAR | 3000 |

### 12. autoimmune_patient_timelines

Patient diagnostic timeline events for odyssey analysis.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| patient_id | VARCHAR | 64 |
| event_type | VARCHAR | 128 |
| event_date | VARCHAR | 32 |
| description | VARCHAR | 2000 |
| provider | VARCHAR | 256 |
| specialty | VARCHAR | 128 |
| days_from_first_symptom | INT64 | -- |

### 13. autoimmune_cross_disease

Cross-disease overlap syndromes and shared pathogenic mechanisms.

| Field | Type | Max Length |
|-------|------|-----------|
| id | VARCHAR (PK) | 128 |
| embedding | FLOAT_VECTOR | 384 dim |
| text_chunk | VARCHAR | 3000 |
| primary_disease | VARCHAR | 256 |
| associated_conditions | VARCHAR | 1024 |
| shared_pathways | VARCHAR | 1024 |
| shared_biomarkers | VARCHAR | 1024 |
| overlap_criteria | VARCHAR | 2000 |
| co_occurrence_rate | FLOAT | -- |

### 14. genomic_evidence

Shared read-only collection managed by the HCLS AI Factory core genomics
pipeline. The autoimmune agent reads from this collection but does not
create or modify it.

---

## 7. Pydantic Data Models

All data models are defined in `src/models.py`.

### Enumerations

**AutoimmuneDisease** -- 13 supported autoimmune conditions:

```
RHEUMATOID_ARTHRITIS       = "rheumatoid_arthritis"
SYSTEMIC_LUPUS             = "systemic_lupus_erythematosus"
MULTIPLE_SCLEROSIS         = "multiple_sclerosis"
TYPE_1_DIABETES            = "type_1_diabetes"
INFLAMMATORY_BOWEL         = "inflammatory_bowel_disease"
PSORIASIS                  = "psoriasis"
ANKYLOSING_SPONDYLITIS     = "ankylosing_spondylitis"
SJOGRENS_SYNDROME          = "sjogrens_syndrome"
SYSTEMIC_SCLEROSIS         = "systemic_sclerosis"
MYASTHENIA_GRAVIS          = "myasthenia_gravis"
CELIAC_DISEASE             = "celiac_disease"
GRAVES_DISEASE             = "graves_disease"
HASHIMOTO_THYROIDITIS      = "hashimoto_thyroiditis"
```

**DiseaseActivityLevel** -- 5 levels: REMISSION, LOW, MODERATE, HIGH, VERY_HIGH

**FlareRisk** -- 4 levels: LOW, MODERATE, HIGH, IMMINENT

### Core Models

**AutoantibodyResult** -- Single autoantibody test result:
- `antibody: str` -- Name (e.g., "ANA", "anti-dsDNA")
- `value: float` -- Measured value
- `unit: str` -- Unit of measurement
- `reference_range: str` -- Normal range
- `positive: bool` -- Positive/negative
- `titer: Optional[str]` -- Titer (e.g., "1:320")
- `pattern: Optional[str]` -- ANA staining pattern

**AutoantibodyPanel** -- Complete panel:
- `patient_id: str`
- `collection_date: str` -- ISO-8601
- `results: List[AutoantibodyResult]`
- Properties: `positive_antibodies`, `positive_count`

**HLAProfile** -- HLA typing results:
- `hla_a: List[str]` -- HLA-A alleles
- `hla_b: List[str]` -- HLA-B alleles
- `hla_c: List[str]` -- HLA-C alleles
- `hla_drb1: List[str]` -- HLA-DRB1 alleles
- `hla_dqb1: List[str]` -- HLA-DQB1 alleles
- Property: `all_alleles` -- Combined list

**DiseaseActivityScore** -- Activity score result:
- `disease: AutoimmuneDisease`
- `score_name: str` -- e.g., "DAS28-CRP", "SLEDAI-2K"
- `score_value: float`
- `level: DiseaseActivityLevel`
- `components: Dict[str, float]`
- `thresholds: Dict[str, float]`

**FlarePredictor** -- Flare prediction:
- `disease: AutoimmuneDisease`
- `current_activity: DiseaseActivityLevel`
- `predicted_risk: FlareRisk`
- `risk_score: float` -- 0.0 to 1.0
- `contributing_factors: List[str]`
- `protective_factors: List[str]`
- `recommended_monitoring: List[str]`
- `time_horizon_days: int` -- Default 90

**BiologicTherapy** -- Therapy recommendation:
- `drug_name: str`
- `drug_class: str`
- `mechanism: str`
- `indicated_diseases: List[AutoimmuneDisease]`
- `pgx_considerations: List[str]`
- `contraindications: List[str]`
- `monitoring_requirements: List[str]`
- `efficacy_evidence: str`

**AutoimmunePatientProfile** -- Input profile:
- `patient_id: str`
- `age: int` (0-150)
- `sex: str` (M/F)
- `autoantibody_panel: Optional[AutoantibodyPanel]`
- `hla_profile: Optional[HLAProfile]`
- `biomarkers: Dict[str, float]`
- `genotypes: Dict[str, str]`
- `diagnosed_conditions: List[AutoimmuneDisease]`
- `current_medications: List[str]`
- `symptom_duration_months: Optional[int]`
- `family_history: List[str]`
- Validator: at least one of autoantibody_panel, hla_profile, or biomarkers required

**AutoimmuneAnalysisResult** -- Output result:
- `patient_id: str`
- `disease_activity_scores: List[DiseaseActivityScore]`
- `flare_predictions: List[FlarePredictor]`
- `hla_associations: List[Dict[str, Any]]`
- `biologic_recommendations: List[BiologicTherapy]`
- `critical_alerts: List[str]`
- `cross_agent_findings: List[Dict[str, Any]]`

### Collection Record Models

**ClinicalDocumentRecord** -- For autoimmune_clinical_documents:
- Fields: id, text_chunk, patient_id, doc_type, specialty, provider,
  visit_date, source_file, page_number, chunk_index
- Method: `to_embedding_text()` -- combines text with doc_type and specialty

**LabResultRecord** -- For autoimmune_patient_labs:
- Fields: id, text_chunk, patient_id, test_name, value, unit,
  reference_range, flag, collection_date, panel_name
- Method: `to_embedding_text()` -- structured lab result format

**TimelineEventRecord** -- For autoimmune_patient_timelines:
- Fields: id, text_chunk, patient_id, event_type, event_date,
  description, provider, specialty, days_from_first_symptom
- Method: `to_embedding_text()` -- timestamped event format

---

## 8. Configuration Reference

All settings use the `AUTO_` environment variable prefix. Defined in
`config/settings.py` as `AutoimmuneSettings(BaseSettings)`.

### Paths

| Setting | Default | Description |
|---------|---------|-------------|
| PROJECT_ROOT | (auto-detected) | Project root directory |
| DATA_DIR | PROJECT_ROOT/data | Data directory (property) |
| CACHE_DIR | DATA_DIR/cache | Cache directory (property) |
| REFERENCE_DIR | DATA_DIR/reference | Reference data (property) |
| DEMO_DATA_DIR | PROJECT_ROOT/demo_data | Demo patient data (property) |

### Milvus

| Setting | Default | Env Var |
|---------|---------|---------|
| MILVUS_HOST | localhost | AUTO_MILVUS_HOST |
| MILVUS_PORT | 19530 | AUTO_MILVUS_PORT |

### Collection Names (14 collections)

| Setting | Default Value |
|---------|--------------|
| COLL_CLINICAL_DOCUMENTS | autoimmune_clinical_documents |
| COLL_PATIENT_LABS | autoimmune_patient_labs |
| COLL_AUTOANTIBODY_PANELS | autoimmune_autoantibody_panels |
| COLL_HLA_ASSOCIATIONS | autoimmune_hla_associations |
| COLL_DISEASE_CRITERIA | autoimmune_disease_criteria |
| COLL_DISEASE_ACTIVITY | autoimmune_disease_activity |
| COLL_FLARE_PATTERNS | autoimmune_flare_patterns |
| COLL_BIOLOGIC_THERAPIES | autoimmune_biologic_therapies |
| COLL_PGX_RULES | autoimmune_pgx_rules |
| COLL_CLINICAL_TRIALS | autoimmune_clinical_trials |
| COLL_LITERATURE | autoimmune_literature |
| COLL_PATIENT_TIMELINES | autoimmune_patient_timelines |
| COLL_CROSS_DISEASE | autoimmune_cross_disease |
| COLL_GENOMIC_EVIDENCE | genomic_evidence |

### Embedding

| Setting | Default | Env Var |
|---------|---------|---------|
| EMBEDDING_MODEL | BAAI/bge-small-en-v1.5 | AUTO_EMBEDDING_MODEL |
| EMBEDDING_DIM | 384 | AUTO_EMBEDDING_DIM |
| EMBEDDING_BATCH_SIZE | 32 | AUTO_EMBEDDING_BATCH_SIZE |
| BGE_INSTRUCTION | "Represent this sentence for searching relevant passages: " | AUTO_BGE_INSTRUCTION |

### LLM

| Setting | Default | Env Var |
|---------|---------|---------|
| ANTHROPIC_API_KEY | "" | AUTO_ANTHROPIC_API_KEY |
| LLM_MODEL | claude-sonnet-4-6 | AUTO_LLM_MODEL |
| LLM_MAX_TOKENS | 4096 | AUTO_LLM_MAX_TOKENS |
| LLM_TEMPERATURE | 0.2 | AUTO_LLM_TEMPERATURE |

### RAG Parameters

| Setting | Default | Env Var |
|---------|---------|---------|
| TOP_K_PER_COLLECTION | 5 | AUTO_TOP_K_PER_COLLECTION |
| SCORE_THRESHOLD | 0.40 | AUTO_SCORE_THRESHOLD |
| MAX_EVIDENCE_ITEMS | 30 | AUTO_MAX_EVIDENCE_ITEMS |
| CONVERSATION_MEMORY_SIZE | 3 | AUTO_CONVERSATION_MEMORY_SIZE |

### Collection Weights (sum = 1.00)

| Setting | Default | Env Var |
|---------|---------|---------|
| WEIGHT_CLINICAL_DOCUMENTS | 0.18 | AUTO_WEIGHT_CLINICAL_DOCUMENTS |
| WEIGHT_PATIENT_LABS | 0.14 | AUTO_WEIGHT_PATIENT_LABS |
| WEIGHT_AUTOANTIBODY_PANELS | 0.12 | AUTO_WEIGHT_AUTOANTIBODY_PANELS |
| WEIGHT_HLA_ASSOCIATIONS | 0.08 | AUTO_WEIGHT_HLA_ASSOCIATIONS |
| WEIGHT_DISEASE_CRITERIA | 0.08 | AUTO_WEIGHT_DISEASE_CRITERIA |
| WEIGHT_DISEASE_ACTIVITY | 0.07 | AUTO_WEIGHT_DISEASE_ACTIVITY |
| WEIGHT_FLARE_PATTERNS | 0.06 | AUTO_WEIGHT_FLARE_PATTERNS |
| WEIGHT_BIOLOGIC_THERAPIES | 0.06 | AUTO_WEIGHT_BIOLOGIC_THERAPIES |
| WEIGHT_CLINICAL_TRIALS | 0.05 | AUTO_WEIGHT_CLINICAL_TRIALS |
| WEIGHT_LITERATURE | 0.05 | AUTO_WEIGHT_LITERATURE |
| WEIGHT_PGX_RULES | 0.04 | AUTO_WEIGHT_PGX_RULES |
| WEIGHT_PATIENT_TIMELINES | 0.03 | AUTO_WEIGHT_PATIENT_TIMELINES |
| WEIGHT_CROSS_DISEASE | 0.02 | AUTO_WEIGHT_CROSS_DISEASE |
| WEIGHT_GENOMIC_EVIDENCE | 0.02 | AUTO_WEIGHT_GENOMIC_EVIDENCE |

A model_validator warns if weights do not sum to approximately 1.0.

### Ports

| Setting | Default | Env Var |
|---------|---------|---------|
| STREAMLIT_PORT | 8531 | AUTO_STREAMLIT_PORT |
| API_PORT | 8532 | AUTO_API_PORT |

### Authentication

| Setting | Default | Env Var |
|---------|---------|---------|
| API_KEY | "" | AUTO_API_KEY |
| CORS_ORIGINS | http://localhost:8080,http://localhost:8531 | AUTO_CORS_ORIGINS |
| MAX_REQUEST_SIZE_MB | 50 | AUTO_MAX_REQUEST_SIZE_MB |

### Document Processing

| Setting | Default | Env Var |
|---------|---------|---------|
| MAX_CHUNK_SIZE | 2500 | AUTO_MAX_CHUNK_SIZE |
| CHUNK_OVERLAP | 200 | AUTO_CHUNK_OVERLAP |
| PDF_DPI | 200 | AUTO_PDF_DPI |

### Relevance and Risk Thresholds

| Setting | Default | Env Var |
|---------|---------|---------|
| CITATION_HIGH | 0.80 | AUTO_CITATION_HIGH |
| CITATION_MEDIUM | 0.60 | AUTO_CITATION_MEDIUM |
| FLARE_RISK_IMMINENT | 0.8 | AUTO_FLARE_RISK_IMMINENT |
| FLARE_RISK_HIGH | 0.6 | AUTO_FLARE_RISK_HIGH |
| FLARE_RISK_MODERATE | 0.4 | AUTO_FLARE_RISK_MODERATE |

### Operational

| Setting | Default | Env Var |
|---------|---------|---------|
| MAX_EVIDENCE_TEXT_LENGTH | 1500 | AUTO_MAX_EVIDENCE_TEXT_LENGTH |
| MAX_KNOWLEDGE_CONTEXT_ITEMS | 25 | AUTO_MAX_KNOWLEDGE_CONTEXT_ITEMS |
| STREAMING_ENABLED | True | AUTO_STREAMING_ENABLED |
| REQUEST_TIMEOUT_SECONDS | 60 | AUTO_REQUEST_TIMEOUT_SECONDS |
| MILVUS_TIMEOUT_SECONDS | 10 | AUTO_MILVUS_TIMEOUT_SECONDS |
| LLM_MAX_RETRIES | 3 | AUTO_LLM_MAX_RETRIES |
| LOG_LEVEL | INFO | AUTO_LOG_LEVEL |
| LOG_DIR | "" | AUTO_LOG_DIR |
| METRICS_ENABLED | True | AUTO_METRICS_ENABLED |

---

## 9. Embedding Strategy

### Model

- **Model:** BAAI/bge-small-en-v1.5
- **Dimensions:** 384
- **Type:** Asymmetric bi-encoder
- **Library:** sentence-transformers

### Asymmetric Search

Queries and documents are embedded differently:

- **Queries:** Prefixed with instruction: `"Represent this sentence for searching relevant passages: "`
- **Documents:** Embedded without instruction prefix

This asymmetric approach improves retrieval quality by encoding the search
intent into the query vector.

### Caching

The RAG engine maintains an in-memory embedding cache:
- Maximum 256 entries
- Key: first 512 characters of input text
- Eviction: FIFO (oldest entry removed when cache is full)

### Batching

Document ingestion uses batch embedding (batch_size=32) via
`sentence-transformers` encode method for throughput.

---

## 10. Autoantibody Engine

Defined in `src/agent.py` method `interpret_autoantibodies()`.

Maps positive antibodies against `AUTOANTIBODY_DISEASE_MAP` (24 entries in
`src/knowledge.py`). For each positive result, returns:

- Antibody name, value, titer, and staining pattern
- Associated disease with sensitivity and specificity
- Clinical notes

### Autoantibody Database (24 entries)

| Antibody | Primary Disease | Sensitivity | Specificity |
|----------|----------------|-------------|-------------|
| ANA | SLE | 0.95 | 0.65 |
| anti-dsDNA | SLE | 0.70 | 0.95 |
| anti-Smith | SLE | 0.25 | 0.99 |
| RF | RA | 0.70 | 0.85 |
| anti-CCP | RA | 0.67 | 0.95 |
| anti-Scl-70 | SSc (diffuse) | 0.35 | 0.99 |
| anti-centromere | SSc (limited/CREST) | 0.40 | 0.98 |
| anti-SSA/Ro | Sjogren's | 0.70 | 0.90 |
| anti-SSB/La | Sjogren's | 0.40 | 0.95 |
| anti-Jo-1 | Antisynthetase syndrome | 0.30 | 0.99 |
| AChR antibody | Myasthenia gravis | 0.85 | 0.99 |
| anti-tTG IgA | Celiac disease | 0.93 | 0.97 |
| TSI | Graves disease | 0.90 | 0.95 |
| anti-TPO | Hashimoto thyroiditis | 0.90 | 0.85 |
| anti-RNP | MCTD | 0.95 | 0.85 |
| anti-histone | Drug-induced lupus | 0.95 | 0.50 |
| ANCA (c-ANCA/PR3) | GPA | 0.90 | 0.95 |
| ANCA (p-ANCA/MPO) | MPA | 0.70 | 0.90 |
| anti-Pm-Scl | SSc-myositis overlap | 0.10 | 0.98 |
| anti-RNA Pol III | SSc (diffuse, renal crisis) | 0.20 | 0.99 |
| anti-cardiolipin | APS | 0.80 | 0.80 |
| lupus anticoagulant | APS | 0.55 | 0.95 |
| anti-beta2-GP I | APS | 0.70 | 0.90 |
| anti-MuSK | MG (AChR-negative) | 0.40 | 0.99 |

---

## 11. HLA Association Engine

Defined in `src/agent.py` method `analyze_hla_associations()`.

Maps HLA alleles against `HLA_DISEASE_ASSOCIATIONS` (22 alleles in
`src/knowledge.py`). Supports exact match and broad allele group matching
(e.g., B*27:05 matches B*27). Results sorted by odds ratio descending.

### Key HLA Associations

| Allele | Disease | Odds Ratio |
|--------|---------|-----------|
| HLA-B*27:05 | Ankylosing Spondylitis | 87.4 |
| HLA-B*27:02 | Ankylosing Spondylitis | 50.0 |
| HLA-B*27:05 | Reactive Arthritis | 20.0 |
| HLA-C*06:02 | Psoriasis | 10.0 |
| HLA-DQB1*02:01 | Celiac Disease | 7.0 |
| HLA-DQA1*05:01 | Celiac Disease | 7.0 |
| HLA-DRB1*03:01 | Celiac Disease | 7.0 |
| HLA-DQB1*03:02 | Type 1 Diabetes | 6.5 |
| HLA-B*51:01 | Behcet's Disease | 5.9 |
| HLA-DRB1*04:01 | Rheumatoid Arthritis | 4.2 |

Protective alleles are also tracked:
- HLA-DRB1*13:01: OR=0.2 for T1D (strongest protection)
- HLA-DRB1*07:01: OR=0.3 for T1D

---

## 12. Disease Activity Engine

Defined in `src/agent.py` method `calculate_disease_activity()`.

Uses `DISEASE_ACTIVITY_THRESHOLDS` (20 scoring systems in `src/knowledge.py`)
to calculate standardized disease activity levels from biomarker data.

### Scoring Systems (20 total)

| Score | Disease | Range | Remission | Low | Moderate | High |
|-------|---------|-------|-----------|-----|----------|------|
| DAS28-CRP | RA | 0-10 | <2.6 | <3.2 | <5.1 | >=5.1 |
| DAS28-ESR | RA | 0-10 | <2.6 | <3.2 | <5.1 | >=5.1 |
| CDAI | RA | 0-76 | <2.8 | <10 | <22 | >=22 |
| SDAI | RA | 0-86 | <3.3 | <11 | <26 | >=26 |
| SLEDAI-2K | SLE | 0-105 | 0 | <4 | <8 | >=12 |
| BASDAI | AS | 0-10 | <2 | <3 | <4 | >=4 |
| ASDAS | AS | 0-6 | <1.3 | <2.1 | <3.5 | >=3.5 |
| PASI | Psoriasis | 0-72 | <1 | <5 | <10 | >=10 |
| DAPSA | Psoriatic Arthritis | 0-164 | <4 | <14 | <28 | >=28 |
| Mayo Score | IBD | 0-12 | <2 | <5 | <8 | >=8 |
| Harvey-Bradshaw | IBD | 0-30 | <4 | <7 | <16 | >=16 |
| ESSDAI | Sjogren's | 0-123 | <1 | <5 | <14 | >=14 |
| mRSS | SSc | 0-51 | <5 | <14 | <29 | >=29 |
| EDSS | MS | 0-10 | <1.5 | <3.5 | <6.0 | >=6.0 |
| QMGS | MG | 0-39 | <3 | <10 | <20 | >=20 |
| MG-ADL | MG | 0-24 | <1 | <5 | <10 | >=10 |
| Marsh Score | Celiac | 0-4 | 0 | <1 | <2 | >=3 |
| Burch-Wartofsky | Graves | 0-140 | <10 | <25 | <45 | >=45 |
| HbA1c-T1D | Type 1 Diabetes | 4-14 | <5.7 | <7.0 | <9.0 | >=9.0 |
| TSH-Hashimoto | Hashimoto's Thyroiditis | 0-100 | <4.5 | <10 | <20 | >=20 |

---

## 13. Flare Prediction Engine

Defined in `src/agent.py` method `predict_flares()`.

Uses `FLARE_BIOMARKER_PATTERNS` (13 disease patterns in `src/knowledge.py`)
to assess flare risk from current biomarker values.

### Algorithm

1. Start with base risk of 0.3
2. For each early warning biomarker in the disease pattern:
   - Elevated inflammatory markers (CRP, ESR, IL-6, calprotectin > 5): +0.15
   - Low complement (C3, C4 < 80): +0.15
   - Low albumin (< 3.5): +0.10
   - Normal values: added to protective factors
3. Clamp total score to [0.0, 1.0]
4. Apply thresholds: imminent >= 0.8, high >= 0.6, moderate >= 0.4
5. Generate monitoring recommendations from the first 3 early warning biomarkers

### Disease Patterns (13)

Each pattern contains:
- `early_warning_biomarkers` -- Ordered list of biomarkers to monitor
- `thresholds` -- Specific trigger values with time windows
- `protective_signals` -- Factors that reduce flare risk

Covered diseases: RA, SLE, IBD, AS, Psoriasis, Sjogren's, SSc, MS,
T1D, MG, Celiac, Graves, Hashimoto.

---

## 14. Biologic Therapy Engine

Defined in `src/agent.py` method `recommend_biologics()`.

Filters `BIOLOGIC_THERAPIES` (22 drugs in `src/knowledge.py`) by indicated
diseases and returns therapy recommendations with PGx context.

### Drug Classes (22 drugs)

| Class | Drugs | Indicated For |
|-------|-------|--------------|
| TNF inhibitor | Adalimumab, Etanercept, Infliximab, Golimumab, Certolizumab pegol | RA, Psoriasis, AS, IBD |
| Anti-CD20 | Rituximab | RA, SLE, MG |
| Anti-CD20 (humanized) | Ocrelizumab | MS |
| IL-6R inhibitor | Tocilizumab, Sarilumab | RA |
| IL-17A inhibitor | Secukinumab, Ixekizumab | Psoriasis, AS |
| BLyS inhibitor | Belimumab | SLE |
| IL-12/23 inhibitor | Ustekinumab | Psoriasis, IBD |
| IL-23 inhibitor (p19) | Guselkumab, Risankizumab | Psoriasis, IBD |
| T-cell co-stim modulator | Abatacept | RA |
| JAK inhibitor | Tofacitinib, Baricitinib, Upadacitinib | RA, AS, IBD, Psoriasis |
| Anti-VLA4 integrin | Natalizumab | MS |
| Integrin inhibitor | Vedolizumab | IBD |
| TYK2 inhibitor | Deucravacitinib | Psoriasis |

---

## 15. Diagnostic Engine

Defined in `src/diagnostic_engine.py` class `DiagnosticEngine`.

### Classification Criteria (10 diseases)

| Disease | Criteria Set | Threshold |
|---------|-------------|-----------|
| RA | 2010 ACR/EULAR RA Classification | >= 6 points |
| SLE | 2019 ACR/EULAR SLE Classification | >= 10 points + ANA |
| AS | ASAS Axial SpA Classification | Imaging arm OR clinical arm |
| SSc | 2013 ACR/EULAR SSc Classification | >= 9 points |
| Sjogren's | 2016 ACR/EULAR Sjogren's | >= 4 points |
| MS | 2017 McDonald Criteria | >= 2 (DIS + DIT) |
| MG | Clinical Diagnostic Criteria | >= 3 points |
| Celiac | ESPGHAN Diagnostic Criteria | >= 3 points |
| IBD | Montreal Classification | >= 3 points |
| Psoriasis | Clinical Diagnostic Criteria | >= 3 points |

### Differential Diagnosis

Generates ranked differential from:
- Autoantibody specificity (weighted by specificity * 2.0)
- HLA odds ratios (log2-scaled contribution * 0.5)

### Overlap Syndromes (9)

| Syndrome | Required Markers | Diseases |
|----------|-----------------|----------|
| Mixed Connective Tissue Disease | anti-RNP | SLE, SSc, RA |
| SLE-RA Overlap | anti-CCP, ANA, RF | SLE, RA |
| Sjogren's-SLE Overlap | anti-SSA/Ro, ANA | Sjogren's, SLE |
| POTS/hEDS/MCAS Triad | tilt table, Beighton >=5, tryptase | POTS, hEDS, MCAS |
| SSc-Myositis Overlap | anti-Pm-Scl | SSc + myositis |
| T1D-Celiac Overlap | anti-tTG, HLA-DQ2/DQ8 | T1D, Celiac |
| Thyroid-T1D Overlap (APS-2) | anti-TPO, anti-GAD | Hashimoto, T1D |
| RA-Sjogren's Overlap | RF, anti-SSA/Ro, ANA | RA, Sjogren's |
| Lupus-APS Overlap | anti-cardiolipin, LAC | SLE, APS |

---

## 16. RAG Engine

Defined in `src/rag_engine.py` class `AutoimmuneRAGEngine`.

### Search Flow

1. **Disease area detection:** Scan query for keywords matching 13 disease areas
   (plus POTS/EDS/MCAS triad)
2. **Embed query:** BGE-small-en-v1.5 with instruction prefix, cached (256 entries)
3. **Build filter expressions:** Patient ID filter with injection prevention
   (regex `^[A-Za-z0-9 _\-\.]+$`, max 64 chars); year range filter for literature
4. **Parallel search:** ThreadPoolExecutor (max_workers=6) across all/selected
   collections, top_k per collection (default 5), score threshold 0.40
5. **Weighted scoring:** `weighted_score = min(raw_score * (1 + weight), 1.0)`
6. **Deduplication:** By record ID and text content hash (MD5 of first 300 chars)
7. **Sort and cap:** Descending by raw score, capped at MAX_EVIDENCE_ITEMS (30)
8. **Knowledge augmentation:** Append relevant entries from HLA_DISEASE_ASSOCIATIONS,
   AUTOANTIBODY_DISEASE_MAP, BIOLOGIC_THERAPIES, FLARE_BIOMARKER_PATTERNS
   (capped at 25 items)
9. **LLM synthesis:** Build messages with system prompt, conversation history
   (last 3 exchanges, trimmed), evidence block, and instructions

### Collection Weight Table

See [Configuration Reference](#8-configuration-reference) for the 14 weights
summing to 1.00.

### Citation Scoring

| Score | Level |
|-------|-------|
| >= 0.80 | High relevance |
| >= 0.60 | Medium relevance |
| < 0.60 | Low relevance |

### System Prompt

The system prompt identifies the agent as the "Precision Autoimmune Intelligence
Agent" and instructs it to:
- Cite evidence using structured formats: `[AutoAb:name]`, `[HLA:allele]`,
  `[Activity:score_name]`, `[Therapy:drug]`, `[Literature:PMID]`, `[Trial:NCT_ID]`
- Distinguish confirmed diagnoses from differentials
- Flag critical alerts
- Consider PGx implications
- Note diagnostic odyssey timeline
- Provide actionable recommendations
- State when data is insufficient

### Conversation Memory

Thread-safe deque with configurable size (default 3). Each entry stores
the question (first 200 chars) and answer (first 800 chars) for context
in subsequent queries.

---

## 17. Agent Orchestrator

The `AutoimmuneAgent.analyze_patient()` method in `src/agent.py` orchestrates
the complete analysis pipeline:

```
Input: AutoimmunePatientProfile
  |
  |--> 1. interpret_autoantibodies(panel)
  |      Positive antibodies -> AUTOANTIBODY_DISEASE_MAP
  |      Returns: List[Dict] with antibody, disease, sensitivity, specificity, value, titer, pattern
  |
  |--> 2. analyze_hla_associations(hla_profile)
  |      All alleles -> HLA_DISEASE_ASSOCIATIONS (exact + broad match)
  |      Returns: List[Dict] sorted by odds_ratio descending
  |
  |--> 3. calculate_disease_activity(biomarkers, conditions)
  |      For each condition -> find scoring system -> CRP/ESR-based calculation
  |      Returns: List[DiseaseActivityScore] with level assignment
  |
  |--> 4. predict_flares(biomarkers, conditions)
  |      For each condition -> FLARE_BIOMARKER_PATTERNS -> risk scoring
  |      Returns: List[FlarePredictor] with contributing/protective factors
  |
  |--> 5. recommend_biologics(conditions, genotypes)
  |      Filter BIOLOGIC_THERAPIES by indicated diseases
  |      Returns: List[BiologicTherapy] with PGx considerations
  |
  |--> 6. _generate_alerts(result)
         High disease activity, imminent flare, strong HLA (OR>5)
         Returns: List[str] critical alerts

Output: AutoimmuneAnalysisResult
```

### Cross-Agent Enrichment

The `enrich_analysis_with_cross_agent()` method adds context from:
- **Biomarker Agent:** Longitudinal trends for contributing factors
  (stub: returns structured placeholder)
- **Imaging Agent:** Joint/organ assessment based on disease type
  (stub: returns structured placeholder)

Disease-to-region mapping for imaging requests:
- RA: hands, feet, knees
- AS: spine, sacroiliac joints
- SSc: lungs, heart
- SLE: kidneys, joints

### Event Publishing

The `publish_diagnosis_event()` method publishes diagnosis events to an
event bus (stub) for consumption by Biomarker, Imaging, and Oncology agents.
Event structure:

```json
{
  "event_type": "autoimmune_diagnosis",
  "source_agent": "precision_autoimmune",
  "patient_id": "...",
  "disease": "...",
  "confidence": 0.0,
  "supporting_evidence": [],
  "timestamp": null
}
```

---

## 18. Export Pipeline

Defined in `src/export.py` class `AutoimmuneExporter`.

### Markdown Export

Generates a structured report with sections:
- Header: patient ID, generation timestamp, knowledge base version
- Critical Alerts (if any)
- Disease Activity Scores table
- Flare Risk Predictions with contributing/protective factors
- HLA-Disease Associations table
- Biologic Therapy Recommendations with PGx and monitoring
- Evidence Sources with relevance badges
- Footer with clinical validation disclaimer

### PDF Export

Uses ReportLab to generate styled PDF reports:
- Page size: Letter
- Title style: NVIDIA green (#76B900), 18pt
- Heading style: NVIDIA green, 14pt
- Alert style: Red (#FF4444), 11pt
- Tables: Green header row, alternating row colors
- Content: Disease activity scores, biologic recommendations
- Footer: 8pt gray disclaimer

Falls back to UTF-8 encoded Markdown if ReportLab is not installed.

### FHIR R4 Bundle Structure

```json
{
  "resourceType": "Bundle",
  "type": "collection",
  "timestamp": "ISO-8601",
  "meta": {
    "profile": ["http://hl7.org/fhir/StructureDefinition/Bundle"],
    "source": "HCLS AI Factory - Precision Autoimmune Agent"
  },
  "entry": [
    { "resource": { "resourceType": "Patient", "id": "...", "identifier": [...] } },
    { "resource": { "resourceType": "Observation", "id": "activity-0", "code": { "text": "DAS28-CRP (...)" }, "valueQuantity": {...}, "interpretation": [...] } },
    { "resource": { "resourceType": "Observation", "id": "flare-risk-0", "code": { "text": "Flare Risk Prediction (...)" }, "valueQuantity": {...}, "interpretation": [...] } },
    { "resource": { "resourceType": "DiagnosticReport", "id": "autoimmune-report-...", "status": "final", "category": [{ "coding": [{ "code": "LAB" }] }], "code": { "text": "Autoimmune Disease Analysis" }, "conclusion": "...", "result": [{ "reference": "Observation/..." }] } }
  ]
}
```

Resources included:
- **Patient:** Identifier with system `urn:hcls-ai-factory:patient`
- **Observation (activity):** One per disease activity score, with valueQuantity
  and interpretation level
- **Observation (flare):** One per flare prediction, with probability value
  and risk level interpretation
- **DiagnosticReport:** Links all observations, includes concatenated conclusion
  with critical alerts

---

## 19. FastAPI REST Server

Defined in `api/main.py`. Runs on port 8532 with 2 uvicorn workers.

### Endpoints (14 total)

| Method | Path | Description |
|--------|------|-------------|
| GET | `/` | Service info (name, version, ports) |
| GET | `/health` | Full health check (Milvus, embedder, LLM, uptime, collections, vectors) |
| GET | `/healthz` | Lightweight probe for landing page |
| POST | `/query` | Full RAG query (retrieve + synthesize with Claude) |
| POST | `/query/stream` | Streaming RAG query via SSE (Server-Sent Events) |
| POST | `/search` | Evidence-only search (no LLM synthesis) |
| POST | `/analyze` | Full autoimmune analysis pipeline (5-step) |
| POST | `/differential` | Differential diagnosis from antibodies + HLA |
| POST | `/ingest/upload` | Upload and ingest a clinical PDF |
| POST | `/ingest/demo-data` | Ingest all demo patient data |
| GET | `/collections` | List collections with vector counts |
| POST | `/collections/create` | Create/recreate all 14 collections |
| POST | `/export` | Export report (Markdown, PDF, FHIR) |
| GET | `/metrics` | Prometheus-compatible metrics |

### Middleware Stack

1. **Request size limiter:** Rejects requests exceeding `MAX_REQUEST_SIZE_MB` (50 MB)
2. **API key authentication:** Optional; skipped if `AUTO_API_KEY` is empty.
   Checks `X-API-Key` header or `api_key` query parameter. Health and root
   endpoints are always accessible.
3. **Timing header:** Adds `X-Process-Time-Ms` to all responses

### Startup Lifecycle

The lifespan context manager initializes:
1. Centralized logging via loguru
2. Milvus connection (2 retry attempts, 2s between)
3. Embedding model (BAAI/bge-small-en-v1.5)
4. LLM client (Anthropic)
5. RAG engine
6. AutoimmuneAgent
7. DocumentProcessor
8. DiagnosticEngine
9. TimelineBuilder
10. Prints service status banner

Graceful shutdown disconnects Milvus.

### Request Models

**QueryRequest:**
- `question: str`
- `patient_id: Optional[str]`
- `patient_context: Optional[str]`
- `collections_filter: Optional[List[str]]`
- `top_k: Optional[int]`

**SearchRequest:**
- `question: str`
- `patient_id: Optional[str]`
- `collections_filter: Optional[List[str]]`
- `top_k: Optional[int]`

**DifferentialRequest:**
- `positive_antibodies: List[str]`
- `hla_alleles: Optional[List[str]]`
- `symptoms: Optional[List[str]]`

**ExportRequest:**
- `patient_id: str`
- `format: str` (markdown, fhir, pdf)
- `query_answer: Optional[str]`

---

## 20. Streamlit UI

Defined in `app/autoimmune_ui.py`. Runs on port 8531.

### 10 Tabs

| # | Tab Name | Functionality |
|---|----------|--------------|
| 1 | Clinical Query | RAG-powered Q&A with evidence citations and streaming |
| 2 | Patient Analysis | Full 5-step autoimmune analysis pipeline |
| 3 | Document Ingest | Upload clinical PDFs for patients |
| 4 | Diagnostic Odyssey | Timeline visualization and odyssey analysis |
| 5 | Autoantibody Panel | Interactive antibody interpretation |
| 6 | HLA Analysis | HLA-disease association lookup |
| 7 | Disease Activity | Activity scoring dashboards (20 scoring systems) |
| 8 | Flare Prediction | Biomarker-based flare risk assessment |
| 9 | Therapy Advisor | Biologic therapy recommendations with PGx |
| 10 | Knowledge Base | Collection stats, evidence explorer |

### Theme

- Background: #1a1a2e (dark navy)
- Tab inactive: #16213e
- Tab active: #0f3460 with #76B900 (NVIDIA green) bottom border
- Tab text: #76B900
- Metrics: #16213e background with 8px border radius

---

## 21. Demo Patients

All 9 demo patients with full clinical PDF datasets in `demo_data/`.

| # | Name | Age/Sex | Primary Dx | Key Antibodies | HLA | Key Biomarkers | PDF Count |
|---|------|---------|-----------|----------------|-----|---------------|-----------|
| 1 | Sarah Mitchell | 34F | SLE (lupus nephritis) | ANA 1:640 (homogeneous), anti-dsDNA, anti-Smith | DRB1*03:01 | Low C3/C4, elevated anti-dsDNA titer | ~20 |
| 2 | Maya Rodriguez | 28F | POTS/hEDS/MCAS | -- | -- | Dysautonomia diagnostic odyssey | ~8 |
| 3 | Linda Chen | 45F | Sjogren's | Anti-SSA/Ro+, anti-SSB/La+ | DRB1*03:01 | ESR, IgG, RF, C4 | ~10 |
| 4 | David Park | 45M | AS | -- | B*27:05 | CRP, ESR, IL-17 | 29 |
| 5 | Rachel Thompson | 38F | MCTD | Anti-U1 RNP | -- | Mixed connective tissue disease | ~15 |
| 6 | Emma Williams | 34F | MS (RRMS) | -- | DRB1*15:01 | NfL, IgG index, oligoclonal bands | ~8 |
| 7 | James Cooper | 19M | T1D + Celiac | Anti-GAD65, anti-IA-2, anti-ZnT8, tTG-IgA+ | DQ2 (DQB1*02:01), DQ8 (DQB1*03:02) | HbA1c, C-peptide, ferritin | ~12 |
| 8 | Karen Foster | 48F | SSc (dcSSc) | Anti-Scl-70+ | DRB1*01:01 | CRP, NT-proBNP, DLCO, creatinine | ~12 |
| 9 | Michael Torres | 41M | Graves' Disease | -- | -- | TSI, Burch-Wartofsky scoring | ~10 |

---

## 22. Cross-Agent Integration

The autoimmune agent integrates with other HCLS AI Factory agents through
stub interfaces and event publishing.

### Biomarker Agent Integration

Method: `request_biomarker_context(patient_id, biomarker_names)`

Requests longitudinal biomarker trends from the Biomarker Agent API.
Returns structured context with trend data for each requested marker.
Currently returns stub placeholder with status "Biomarker Agent integration pending".

### Imaging Agent Integration

Method: `request_imaging_context(patient_id, body_regions)`

Requests imaging findings correlated with HLA/antibody patterns.
Disease-to-region mapping:
- RA: hands, feet, knees
- AS: spine, sacroiliac joints
- SSc: lungs, heart
- SLE: kidneys, joints

Currently returns stub placeholder.

### Event Bus

Method: `publish_diagnosis_event(patient_id, disease, confidence, evidence)`

Publishes autoimmune diagnosis events for other agents. In production, this
would use Redis or Kafka. Currently returns a stub event structure.

### Enrichment Pipeline

Method: `enrich_analysis_with_cross_agent(result, patient_id)`

Combines biomarker and imaging context into the analysis result by:
1. Extracting biomarker names from flare prediction contributing factors
2. Requesting biomarker context from the Biomarker Agent
3. Requesting imaging context based on disease type
4. Appending findings to `cross_agent_findings` in the result

---

## 23. Monitoring and Metrics

### Prometheus Metrics Endpoint

`GET /metrics` returns Prometheus-compatible text format:

```
autoimmune_agent_up 1
autoimmune_collection_vectors{collection="autoimmune_clinical_documents"} 1234
autoimmune_collection_vectors{collection="autoimmune_patient_labs"} 567
...
autoimmune_agent_uptime_seconds 3600
```

### Health Endpoints

- `GET /health` -- Full health check returning: status, milvus_connected,
  collections count, total_vectors, embedder_loaded, llm_available, uptime_seconds
- `GET /healthz` -- Lightweight probe returning `{"status": "ok"}`

### Docker Health Checks

- Streamlit: checks `/_stcore/health` every 30s (start period 60s, 3 retries)
- API: checks `/healthz` every 30s (start period 30s, 3 retries)

### Logging

Uses loguru for structured logging. Configurable via:
- `AUTO_LOG_LEVEL` (default: INFO)
- `AUTO_LOG_DIR` (default: PROJECT_ROOT/logs)
- Service name: "autoimmune-api"

---

## 24. Testing

455 tests across 8 test files in `tests/`.

| File | Tests | Coverage Area |
|------|-------|--------------|
| test_rag_engine.py | 93 | RAG engine (search, query, streaming, caching, disease detection) |
| test_production_readiness.py | 59 | Production checks (security, config, error handling) |
| test_api.py | 57 | API endpoints (all 14 endpoints, auth, CORS, error cases) |
| test_collections.py | 57 | Collection manager (create, insert, search, stats) |
| test_diagnostic_engine.py | 51 | Diagnostic engine (criteria, differential, odyssey, overlap) |
| test_timeline_builder.py | 49 | Timeline builder (event extraction, date parsing, classification) |
| test_export.py | 34 | Export formats (Markdown, PDF, FHIR R4) |
| test_autoimmune.py | 31 | Core agent (analyze_patient, antibodies, HLA, flare, biologics) |

Run all tests:

```bash
cd precision_autoimmune_agent
python -m pytest tests/ -v
```

---

## 25. HCLS AI Factory Integration

The Precision Autoimmune Agent is one of 5 intelligence agents in the
HCLS AI Factory platform.

### Agent Family Port Map

| Agent | UI Port | API Port |
|-------|---------|----------|
| CAR-T Intelligence | 8521 | 8522 |
| Medical Imaging | 8523 | 8524 |
| Precision Oncology | 8525 | 8526 |
| Precision Biomarker | 8528 | 8529 |
| Precision Autoimmune | 8531 | 8532 |

### Shared Infrastructure

- **Milvus:** Shared vector database on port 19530 (with etcd and MinIO)
- **Docker network:** `hcls-network` (external)
- **Landing page:** Flask hub on port 8080 with real-time health monitoring
- **Monitoring:** Prometheus + Grafana dashboards
- **Genomic evidence:** Shared read-only collection `genomic_evidence`

### HCLS AI Factory Pipeline Position

The autoimmune agent operates in Stage 2 (RAG/Chat Pipeline) of the
3-stage HCLS AI Factory:

```
Stage 1: Genomics Pipeline (Parabricks/DeepVariant/BWA-MEM2)
  FASTQ -> VCF
    |
Stage 2: RAG/Chat Pipeline (Milvus + Claude AI)
  VCF -> Variant Interpretation -> Autoimmune Analysis
    |     [Autoimmune Agent operates here]
    |
Stage 3: Drug Discovery Pipeline (BioNeMo/MolMIM/DiffDock/RDKit)
  Targets -> Drug Candidates
```

The autoimmune agent receives genomic data (HLA typing, PGx variants)
from Stage 1 and feeds therapy recommendations to Stage 3 when applicable.

---

## 26. Quick Start

### Local Development

```bash
# 1. Clone and enter directory
cd precision_autoimmune_agent

# 2. Install dependencies
pip install -r requirements.txt

# 3. Set environment variables
export ANTHROPIC_API_KEY="your-key-here"
export AUTO_MILVUS_HOST="localhost"
export AUTO_MILVUS_PORT="19530"

# 4. Ensure Milvus is running
docker compose -f /path/to/docker-compose.dgx-spark.yml up -d milvus-standalone

# 5. Create collections and seed knowledge
python scripts/setup_collections.py --seed

# 6. Ingest demo data
curl -X POST http://localhost:8532/ingest/demo-data

# 7. Start API server
uvicorn api.main:app --host 0.0.0.0 --port 8532

# 8. Start Streamlit UI (separate terminal)
streamlit run app/autoimmune_ui.py --server.port 8531

# 9. Open browser
open http://localhost:8531
```

### Docker Deployment

```bash
# Build and start all services
docker compose up -d

# Verify
docker compose ps
curl http://localhost:8532/health
```

---

## 27. Dependencies

From `requirements.txt` (19 packages):

| Package | Version | Purpose |
|---------|---------|---------|
| pydantic | >=2.0 | Data models and validation |
| pydantic-settings | >=2.7 | Configuration management with env vars |
| pymilvus | >=2.4.0 | Milvus vector database client |
| sentence-transformers | >=2.2.0 | BGE-small-en-v1.5 embedding model |
| anthropic | >=0.18.0 | Claude LLM client |
| streamlit | >=1.30.0 | 10-tab UI framework |
| fastapi | >=0.109.0 | REST API framework |
| uvicorn[standard] | >=0.27.0 | ASGI server |
| python-dotenv | >=1.0.0 | .env file loading |
| loguru | >=0.7.0 | Structured logging |
| numpy | >=1.24.0 | Numerical operations |
| pandas | >=2.0.0 | Data manipulation |
| plotly | >=5.18.0 | Interactive charts |
| reportlab | >=4.0.0 | PDF report generation |
| PyPDF2 | >=3.0.0 | PDF text extraction |
| python-multipart | >=0.0.6 | File upload handling |
| prometheus-client | >=0.20.0 | Metrics collection |
| requests | >=2.31.0 | HTTP client |
| httpx | >=0.25.0 | Async HTTP client |
