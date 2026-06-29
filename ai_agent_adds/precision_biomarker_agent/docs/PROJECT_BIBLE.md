# Precision Biomarker Intelligence Agent -- Project Bible

Complete implementation reference for the Precision Biomarker Intelligence
Agent, part of the HCLS AI Factory pipeline: Patient DNA -> Drug Candidates.

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0
**Repository:** `hcls-ai-factory/ai_agent_adds/precision_biomarker_agent`

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
10. [Pharmacogenomic Engine](#10-pharmacogenomic-engine)
11. [Biological Age Engine](#11-biological-age-engine)
12. [Disease Trajectory Engine](#12-disease-trajectory-engine)
13. [Genotype Adjustment Engine](#13-genotype-adjustment-engine)
14. [Critical Value Engine](#14-critical-value-engine)
15. [Discordance Detector](#15-discordance-detector)
16. [Knowledge Graph](#16-knowledge-graph)
17. [RAG Engine](#17-rag-engine)
18. [Agent Orchestrator](#18-agent-orchestrator)
19. [Export Pipeline](#19-export-pipeline)
20. [Report Generator](#20-report-generator)
21. [FastAPI REST Server](#21-fastapi-rest-server)
22. [Streamlit UI](#22-streamlit-ui)
23. [Demo Patients](#23-demo-patients)
24. [Cross-Agent Integration](#24-cross-agent-integration)
25. [Monitoring and Metrics](#25-monitoring-and-metrics)
26. [Audit Logging](#26-audit-logging)
27. [Testing](#27-testing)
28. [Dependencies](#28-dependencies)
29. [Quick Start](#29-quick-start)

---

## 1. Project Overview

The Precision Biomarker Intelligence Agent is a genotype-aware biomarker
interpretation platform that transforms standard blood panel results and
genomic data into actionable precision health intelligence. It operates
within the HCLS AI Factory Precision Intelligence Network -- one of three
GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence
Network, Therapeutic Discovery Engine) -- delivering a complete Patient DNA
to Drug Candidates pipeline in under 5 hours on a single NVIDIA DGX Spark
desktop workstation ($4,699).

### Platform Positioning

The HCLS AI Factory comprises three engines:

1. **Genomic Foundation Engine** -- Parabricks/DeepVariant/BWA-MEM2
   (FASTQ to VCF)
2. **Precision Intelligence Network** -- 11 domain-specialized RAG agents
   including this biomarker agent
3. **Therapeutic Discovery Engine** -- BioNeMo MolMIM/DiffDock/RDKit
   (molecular generation and docking)

### The 11 Peer Agents

| # | Agent | Domain |
|---|-------|--------|
| 1 | Precision Oncology Agent | Molecular tumor board decision support |
| 2 | CAR-T Intelligence Agent | CAR-T cell therapy development lifecycle |
| 3 | **Precision Biomarker Agent** (this agent) | Genotype-aware biomarker interpretation |
| 4 | Clinical Trial Intelligence Agent | Trial landscape monitoring and matching |
| 5 | Cardiology Intelligence Agent | Cardiovascular risk and treatment optimization |
| 6 | Neurology Intelligence Agent | Neurological disease characterization |
| 7 | Pharmacogenomics (PGx) Agent | Drug-gene interaction and dosing guidance |
| 8 | Medical Imaging Intelligence Agent | Imaging AI with NVIDIA NIM microservices |
| 9 | Single-Cell Intelligence Agent | Single-cell omics analysis |
| 10 | Autoimmune Intelligence Agent | Autoimmune disease monitoring and treatment |
| 11 | Rare Disease Intelligence Agent | Rare disease diagnosis and therapeutic matching |

### Cross-Agent Integration

The biomarker agent calls 4 peer agents via `cross_modal/cross_agent.py` and
the `/integrated-assessment` endpoint: Oncology, CAR-T, PGx, and Clinical
Trial agents.

### Pediatric Biomarker Applications as Primary Use Case

Pediatric oncology biomarkers are a primary use case, with support for
CD19/CD22/CRLF2 markers in pediatric ALL, MRD detection for risk
stratification, and MYCN amplification for neuroblastoma prognostication.

**Core capabilities:**

- 14 pharmacogenes with CPIC-guided phenotyping (CYP2D6, CYP2C19, CYP2C9,
  CYP3A5, SLCO1B1, VKORC1, MTHFR, TPMT, DPYD)
- Biological age estimation using PhenoAge (Levine 2018) and GrimAge surrogate
- Disease trajectory prediction across 9 categories (diabetes, cardiovascular,
  liver, thyroid, iron, nutritional, kidney, bone health, cognitive)
- Ashkenazi Jewish carrier screening panel (10 genes)
- Genotype-based reference range adjustments for 7 modifier genes
- Critical value alerting with escalation targets
- Cross-biomarker discordance detection
- Longitudinal biomarker trend tracking
- 12-section clinical report generation with PDF and FHIR R4 export
- Multi-language report translation (English, Spanish, Chinese, Hindi, French,
  Arabic, Portuguese)

**Codebase statistics:**

| Metric                     | Value           |
| -------------------------- | --------------- |
| Total source lines         | ~20,200         |
| Source modules (src/)      | 16 Python files |
| API modules (api/)         | 4 Python files  |
| UI modules (app/)          | 3 Python files  |
| Scripts                    | 7 Python files  |
| Test files                 | 18 test files   |
| Tests                      | 709, all passing|
| Milvus collections         | 14 total        |
| Reference data files       | 18 JSON files   |

---

## 2. Pipeline Pattern

The agent follows the HCLS AI Factory three-stage pipeline:

```
Stage 1: Genomics Pipeline (Parabricks/DeepVariant)
  FASTQ -> VCF -> genomic_evidence collection (shared, read-only)
      |
Stage 2: RAG + Intelligence Agents
      |
      +-> Biomarker Agent (this project)
      |     Multi-collection RAG + 6 analysis engines
      |     14 Milvus collections (13 owned + 1 shared)
      |
      +-> CAR-T Agent, Oncology Agent, Imaging Agent, Autoimmune Agent
      |
Stage 3: Drug Discovery Pipeline (BioNeMo/DiffDock/RDKit)
  PGx findings -> drug candidate optimization
```

**Agent internal pipeline:**

```
Patient Profile
    |
    v
[1. Plan]  -> Identify topics, disease areas, relevant modules
    |
    v
[2. Analyze] -> Run analysis modules in parallel:
    |   - BiologicalAgeCalculator (PhenoAge/GrimAge)
    |   - DiseaseTrajectoryAnalyzer (9 categories)
    |   - PharmacogenomicMapper (14 pharmacogenes)
    |   - GenotypeAdjuster (7 modifier genes)
    |   - CriticalValueEngine (21 rules)
    |   - DiscordanceDetector (cross-biomarker patterns)
    |   - LabRangeInterpreter (standard vs optimal)
    |
    v
[3. Search] -> Multi-collection RAG across 14 collections
    |
    v
[4. Synthesize] -> LLM-grounded synthesis with citations
    |
    v
[5. Report] -> 12-section clinical report (MD/PDF/FHIR)
```

---

## 3. DGX Spark Hardware

| Spec              | Value                                |
| ----------------- | ------------------------------------ |
| GPU               | NVIDIA Grace Blackwell, 128 GB      |
| CPU               | ARM-based Grace CPU, 20 cores        |
| System RAM        | 128 GB unified memory                |
| Storage           | 4 TB NVMe SSD                        |
| CUDA              | 12.x                                 |
| Price             | $4,699                               |
| Power             | ~250 W (desktop form factor)         |

The biomarker agent runs alongside Milvus, embedding model, and the other four
intelligence agents on this single machine. The embedding model
(BGE-small-en-v1.5, 33M parameters) runs on CPU; the LLM (Claude) is called
via the Anthropic API.

---

## 4. Repository Layout

```
precision_biomarker_agent/
|-- api/                         # FastAPI REST server
|   |-- __init__.py
|   |-- main.py                  (465 lines) Entry point, lifespan, core endpoints
|   +-- routes/
|       |-- __init__.py
|       |-- analysis.py          (495 lines) /v1/analyze, /v1/biological-age, etc.
|       |-- events.py            (326 lines) /v1/events/cross-modal, /biomarker-alert
|       +-- reports.py           (296 lines) /v1/report/generate, /{id}/pdf, /fhir
|
|-- app/                         # Streamlit UI
|   |-- __init__.py
|   |-- biomarker_ui.py          (1,863 lines) 8-tab Streamlit application
|   |-- patient_360.py           (670 lines) Cross-agent Patient 360 dashboard
|   +-- protein_viewer.py        (168 lines) 3D protein structure viewer
|
|-- config/
|   +-- settings.py              (139 lines) PrecisionBiomarkerSettings (Pydantic)
|
|-- data/
|   |-- cache/                   # Runtime cache directory
|   +-- reference/               # 18 JSON seed files
|       |-- biomarker_reference.json
|       |-- biomarker_genetic_variants.json
|       |-- biomarker_pgx_rules.json
|       |-- biomarker_disease_trajectories.json
|       |-- biomarker_clinical_evidence.json
|       |-- biomarker_nutrition.json
|       |-- biomarker_drug_interactions.json
|       |-- biomarker_aging_markers.json
|       |-- biomarker_genotype_adjustments.json
|       |-- biomarker_monitoring.json
|       |-- biomarker_critical_values.json
|       |-- biomarker_discordance_rules.json
|       |-- biomarker_aj_carrier_screening.json
|       |-- biomarker_genomic_evidence.json
|       |-- biomarker_lab_ranges.json
|       |-- biomarker_longitudinal_tracking.json
|       |-- biomarker_sample_patients.json
|       +-- nutrition_guidelines_seed.json
|
|-- src/                         # Core engine modules
|   |-- __init__.py
|   |-- models.py                (786 lines) 14 collection + 8 analysis models + 3 enums
|   |-- collections.py           (1,391 lines) Milvus collection management
|   |-- rag_engine.py            (573 lines) Multi-collection RAG with parallel search
|   |-- agent.py                 (610 lines) Autonomous plan-analyze-search-synthesize
|   |-- knowledge.py             (1,326 lines) 6 disease domain knowledge graphs
|   |-- pharmacogenomics.py      (1,503 lines) CPIC-guided phenotyping, 14 pharmacogenes
|   |-- disease_trajectory.py    (1,421 lines) Pre-symptomatic detection, 9 categories
|   |-- biological_age.py        (408 lines) PhenoAge + GrimAge surrogate
|   |-- genotype_adjustment.py   (1,225 lines) Genotype-based reference range adjustments
|   |-- export.py                (1,392 lines) PDF + FHIR R4 + CSV + JSON + Markdown
|   |-- report_generator.py      (993 lines) 12-section clinical reports
|   |-- critical_values.py       (179 lines) Critical threshold detection (21 rules)
|   |-- discordance_detector.py  (299 lines) Cross-biomarker anomaly detection
|   |-- lab_range_interpreter.py (221 lines) Standard vs optimal range interpretation
|   |-- translation.py           (217 lines) Multi-language medical terminology
|   +-- audit.py                 (83 lines) HIPAA-compliant audit logging
|
|-- scripts/                     # Setup, seeding, and validation
|   |-- setup_collections.py     (57 lines) Create Milvus collections
|   |-- seed_all.py              (207 lines) Seed all 14 collections from JSON
|   |-- gen_patient_data.py      (345 lines) Generate sample patient data
|   |-- gen_critical_values.py   (510 lines) Generate critical value rules
|   |-- gen_lab_ranges_and_aj.py (527 lines) Generate lab ranges + AJ screening data
|   |-- expand_biomarker_reference.py  (580 lines) Expand reference data
|   |-- expand_variants_and_interactions.py (436 lines) Expand variant data
|   |-- demo_validation.py       (333 lines) End-to-end demo validation
|   +-- validate_e2e.py          (178 lines) E2E validation script
|
|-- tests/                       # 18 test files, 709 tests
|   |-- conftest.py              Shared fixtures (mock embedder, LLM, manager)
|   |-- test_models.py
|   |-- test_collections.py
|   |-- test_rag_engine.py
|   |-- test_agent.py
|   |-- test_pharmacogenomics.py
|   |-- test_biological_age.py
|   |-- test_disease_trajectory.py
|   |-- test_genotype_adjustment.py
|   |-- test_critical_values.py
|   |-- test_discordance_detector.py
|   |-- test_lab_range_interpreter.py
|   |-- test_report_generator.py
|   |-- test_export.py
|   |-- test_api.py
|   |-- test_integration.py
|   |-- test_edge_cases.py
|   |-- test_longitudinal.py
|   +-- test_ui.py
|
|-- docker-compose.yml           6-service stack
|-- Dockerfile                   Multi-stage Python 3.10-slim
|-- requirements.txt             Python dependencies
|-- pyproject.toml               Project metadata
+-- README.md                    Quick start guide
```

---

## 5. Docker Compose Services

File: `docker-compose.yml` -- 6 services on `biomarker-network` bridge.

| Service              | Image                                     | Port(s)       | Purpose                              |
| -------------------- | ----------------------------------------- | ------------- | ------------------------------------- |
| `milvus-etcd`        | `quay.io/coreos/etcd:v3.5.5`             | (internal)    | Milvus metadata store (4 GB quota)    |
| `milvus-minio`       | `minio/minio:RELEASE.2023-03-20T20-16-18Z` | (internal)  | Milvus object storage (log/index)     |
| `milvus-standalone`  | `milvusdb/milvus:v2.4-latest`            | 19530, 9091   | Vector database (gRPC + health)       |
| `biomarker-streamlit`| Built from `./Dockerfile`                 | 8528:8528     | Streamlit UI (8 tabs)                 |
| `biomarker-api`      | Built from `./Dockerfile`                 | 8529:8529     | FastAPI REST server (uvicorn, 2 workers) |
| `biomarker-setup`    | Built from `./Dockerfile`                 | --            | One-shot: create collections + seed   |

**DGX Spark external port mapping** (from `docker-compose.dgx-spark.yml`):

| Internal Port | External Port | Service          |
| ------------- | ------------- | ---------------- |
| 8528          | 8502          | Streamlit UI     |
| 8529          | 8102          | FastAPI API      |

**Dockerfile:** Multi-stage build on `python:3.10-slim`. Builder stage compiles
native extensions for sentence-transformers/numpy. Runtime stage copies the
virtualenv and runs as non-root `biomarkeruser`. Default CMD launches Streamlit;
the API container overrides with `uvicorn api.main:app`.

**Startup sequence:**

```bash
cp .env.example .env                     # Set ANTHROPIC_API_KEY
docker compose up -d                     # Start all services
docker compose logs -f biomarker-setup   # Watch seed progress
```

Setup runs: `setup_collections.py --drop-existing` then `seed_all.py`, which
reads 14 JSON files from `data/reference/`, embeds text chunks with
BGE-small-en-v1.5, and inserts into Milvus.

---

## 6. Milvus Collection Schemas

**14 collections total:** 13 biomarker-specific (owned) + 1 shared read-only.

**Embedding configuration:**

| Setting     | Value                   |
| ----------- | ----------------------- |
| Model       | BAAI/bge-small-en-v1.5  |
| Dimensions  | 384                     |
| Metric      | COSINE                  |
| Index       | IVF_FLAT                |

### 6.1 biomarker_reference

Reference biomarker definitions with standard and optimal ranges.

| Field                 | Type         | Notes                              |
| --------------------- | ------------ | ---------------------------------- |
| id                    | VARCHAR(100) | Primary key                        |
| embedding             | FLOAT_VECTOR | 384-dim BGE-small                  |
| name                  | VARCHAR(100) | Display name (e.g., "HbA1c")      |
| unit                  | VARCHAR(20)  | e.g., "%", "mg/dL", "ng/mL"       |
| category              | VARCHAR(30)  | CBC, CMP, Lipids, Thyroid, etc.    |
| ref_range_min         | FLOAT        | Standard lower bound               |
| ref_range_max         | FLOAT        | Standard upper bound               |
| text_chunk            | VARCHAR(3000)| Text for embedding                 |
| clinical_significance | VARCHAR(2000)| Clinical interpretation            |
| epigenetic_clock      | VARCHAR(50)  | PhenoAge/GrimAge coefficient       |
| genetic_modifiers     | VARCHAR(500) | Comma-separated modifier genes     |

### 6.2 biomarker_genetic_variants

Genetic variants affecting biomarker levels (MTHFR, APOE, PNPLA3, HFE, etc.).

| Field                | Type         | Notes                              |
| -------------------- | ------------ | ---------------------------------- |
| id                   | VARCHAR(100) | Primary key                        |
| embedding            | FLOAT_VECTOR | 384-dim                            |
| gene                 | VARCHAR(50)  | Gene symbol                        |
| rs_id                | VARCHAR(20)  | dbSNP rsID                         |
| risk_allele          | VARCHAR(5)   | Risk allele                        |
| protective_allele    | VARCHAR(5)   | Protective allele                  |
| effect_size          | VARCHAR(100) | Effect size description            |
| mechanism            | VARCHAR(2000)| Molecular mechanism                |
| disease_associations | VARCHAR(1000)| Disease associations               |
| text_chunk           | VARCHAR(3000)| Text for embedding                 |

### 6.3 biomarker_pgx_rules

CPIC pharmacogenomic dosing rules.

| Field          | Type         | Notes                              |
| -------------- | ------------ | ---------------------------------- |
| id             | VARCHAR(100) | Primary key                        |
| embedding      | FLOAT_VECTOR | 384-dim                            |
| gene           | VARCHAR(50)  | Pharmacogene (CYP2D6, etc.)       |
| star_alleles   | VARCHAR(100) | Star allele combination            |
| drug           | VARCHAR(100) | Drug name                          |
| phenotype      | VARCHAR(30)  | MetabolizerPhenotype enum          |
| cpic_level     | VARCHAR(5)   | CPIC evidence level (1A-3)         |
| recommendation | VARCHAR(2000)| Dosing recommendation              |
| evidence_url   | VARCHAR(500) | CPIC/PharmGKB URL                  |
| text_chunk     | VARCHAR(3000)| Text for embedding                 |

### 6.4 biomarker_disease_trajectories

Disease progression trajectories with intervention windows.

| Field               | Type         | Notes                              |
| ------------------- | ------------ | ---------------------------------- |
| id                  | VARCHAR(100) | Primary key                        |
| embedding           | FLOAT_VECTOR | 384-dim                            |
| disease             | VARCHAR(30)  | DiseaseCategory enum               |
| stage               | VARCHAR(30)  | Progression stage                  |
| biomarker_pattern   | VARCHAR(2000)| JSON thresholds for this stage     |
| years_to_diagnosis  | FLOAT        | Estimated years to diagnosis       |
| intervention_window | VARCHAR(500) | Intervention opportunity           |
| risk_reduction_pct  | FLOAT        | Potential risk reduction %         |
| text_chunk          | VARCHAR(3000)| Text for embedding                 |

### 6.5 biomarker_clinical_evidence

Published clinical evidence with PubMed references.

| Field         | Type         | Notes                              |
| ------------- | ------------ | ---------------------------------- |
| id            | VARCHAR(100) | Primary key                        |
| embedding     | FLOAT_VECTOR | 384-dim                            |
| pmid          | VARCHAR(20)  | PubMed ID                          |
| title         | VARCHAR(500) | Publication title                  |
| finding       | VARCHAR(3000)| Key finding                        |
| year          | INT16        | Publication year                   |
| disease_area  | VARCHAR(100) | Disease area                       |
| text_chunk    | VARCHAR(3000)| Text for embedding                 |

### 6.6 biomarker_nutrition

Genotype-aware nutrition guidelines.

| Field            | Type         | Notes                              |
| ---------------- | ------------ | ---------------------------------- |
| id               | VARCHAR(100) | Primary key                        |
| embedding        | FLOAT_VECTOR | 384-dim                            |
| nutrient         | VARCHAR(100) | Nutrient name                      |
| genetic_context  | VARCHAR(200) | e.g., "MTHFR C677T heterozygous"  |
| recommended_form | VARCHAR(200) | e.g., "methylfolate"              |
| dose_range       | VARCHAR(100) | e.g., "400-800 mcg/day"           |
| evidence_summary | VARCHAR(2000)| Evidence summary                   |
| text_chunk       | VARCHAR(3000)| Text for embedding                 |

### 6.7 biomarker_drug_interactions

Gene-drug interactions beyond PGx (substrate/inhibitor/inducer relationships).

| Field            | Type         | Notes                              |
| ---------------- | ------------ | ---------------------------------- |
| id               | VARCHAR(100) | Primary key                        |
| embedding        | FLOAT_VECTOR | 384-dim                            |
| drug             | VARCHAR(100) | Drug name                          |
| gene             | VARCHAR(50)  | Gene involved                      |
| interaction_type | VARCHAR(50)  | substrate, inhibitor, inducer      |
| text_chunk       | VARCHAR(3000)| Text for embedding                 |

### 6.8 biomarker_aging_markers

Epigenetic aging clock marker data (PhenoAge, GrimAge coefficients).

| Field      | Type         | Notes                              |
| ---------- | ------------ | ---------------------------------- |
| id         | VARCHAR(100) | Primary key                        |
| embedding  | FLOAT_VECTOR | 384-dim                            |
| text_chunk | VARCHAR(3000)| Text for embedding                 |

### 6.9 biomarker_genotype_adjustments

Genotype-based reference range adjustments for 7 modifier genes.

| Field      | Type         | Notes                              |
| ---------- | ------------ | ---------------------------------- |
| id         | VARCHAR(100) | Primary key                        |
| embedding  | FLOAT_VECTOR | 384-dim                            |
| text_chunk | VARCHAR(3000)| Text for embedding                 |

### 6.10 biomarker_monitoring

Condition-specific monitoring protocols (frequency, biomarkers to track).

| Field      | Type         | Notes                              |
| ---------- | ------------ | ---------------------------------- |
| id         | VARCHAR(100) | Primary key                        |
| embedding  | FLOAT_VECTOR | 384-dim                            |
| text_chunk | VARCHAR(3000)| Text for embedding                 |

### 6.11 biomarker_critical_values

Critical threshold definitions for life-threatening lab values.

| Field      | Type         | Notes                              |
| ---------- | ------------ | ---------------------------------- |
| id         | VARCHAR(100) | Primary key                        |
| embedding  | FLOAT_VECTOR | 384-dim                            |
| text_chunk | VARCHAR(3000)| Text for embedding                 |

### 6.12 biomarker_discordance_rules

Cross-biomarker discordance patterns (contradictory or unexpected relationships).

| Field      | Type         | Notes                              |
| ---------- | ------------ | ---------------------------------- |
| id         | VARCHAR(100) | Primary key                        |
| embedding  | FLOAT_VECTOR | 384-dim                            |
| text_chunk | VARCHAR(3000)| Text for embedding                 |

### 6.13 biomarker_aj_carrier_screening

Ashkenazi Jewish carrier screening panel (10 genes: BRCA1/2, GBA, HEXA, FANCC,
ASPA, BLM, SMPD1, IKBKAP/ELP1, MCOLN1).

| Field      | Type         | Notes                              |
| ---------- | ------------ | ---------------------------------- |
| id         | VARCHAR(100) | Primary key                        |
| embedding  | FLOAT_VECTOR | 384-dim                            |
| text_chunk | VARCHAR(3000)| Text for embedding                 |

### 6.14 genomic_evidence (shared, read-only)

Shared VCF-derived genomic variants from the Genomics Pipeline. Populated by
Stage 1 (Parabricks/DeepVariant). The biomarker agent reads but never writes
to this collection.

**Collection search weights** (must sum to ~1.0):

| Collection                 | Weight |
| -------------------------- | ------ |
| biomarker_reference        | 0.12   |
| biomarker_genetic_variants | 0.11   |
| biomarker_pgx_rules        | 0.10   |
| biomarker_disease_trajectories | 0.10 |
| biomarker_clinical_evidence | 0.09  |
| genomic_evidence           | 0.08   |
| biomarker_drug_interactions | 0.07  |
| biomarker_aging_markers    | 0.07   |
| biomarker_nutrition        | 0.05   |
| biomarker_genotype_adjustments | 0.05 |
| biomarker_monitoring       | 0.05   |
| biomarker_critical_values  | 0.04   |
| biomarker_discordance_rules | 0.04  |
| biomarker_aj_carrier_screening | 0.03 |
| **Total**                  | **1.00** |

---

## 7. Pydantic Data Models

File: `src/models.py` (786 lines)

### 7.1 Enums (7 total)

| Enum                 | Values                                                    |
| -------------------- | --------------------------------------------------------- |
| `RiskLevel`          | critical, high, moderate, low, normal                     |
| `ClockType`          | PhenoAge, GrimAge                                         |
| `DiseaseCategory`    | diabetes, cardiovascular, liver, thyroid, iron, nutritional, kidney, bone_health, cognitive |
| `MetabolizerPhenotype` | ultra_rapid, normal, intermediate, poor                 |
| `CPICLevel`          | 1A, 1B, 2A, 2B, 3                                        |
| `Zygosity`           | homozygous_ref, heterozygous, homozygous_alt              |

### 7.2 Collection Models (14)

Each collection model inherits `BaseModel` and provides:
- Schema fields matching the Milvus collection
- `to_embedding_text()` method for generating BGE-small embedding input
- `@model_validator` for cross-field validation

Models: `BiomarkerReference`, `GeneticVariant`, `PGxRule`, `DiseaseTrajectory`,
`ClinicalEvidence`, `NutritionGuideline`, `DrugInteraction`, `AgingMarker`,
`GenotypeAdjustment`, `MonitoringProtocol`, `CriticalValue`, `DiscordanceRule`,
`AJCarrierScreeningEntry`, (shared: `GenomicEvidence`).

### 7.3 Analysis Models (8+)

| Model                      | Purpose                                     |
| -------------------------- | ------------------------------------------- |
| `PatientProfile`           | Patient demographics, biomarkers, genotypes |
| `SearchHit`                | Single RAG search result with score         |
| `CrossCollectionResult`    | Merged results across all collections       |
| `AgentQuery`               | Structured query with filters               |
| `AgentResponse`            | Agent response with citations               |
| `AnalysisResult`           | Full analysis output (all modules)          |
| `BiologicalAgeResult`      | PhenoAge/GrimAge calculation result         |
| `DiseaseTrajectoryResult`  | Disease risk trajectory output              |
| `PGxResult`                | Pharmacogenomic mapping output              |
| `GenotypeAdjustmentResult` | Adjusted reference ranges                   |

---

## 8. Configuration Reference

File: `config/settings.py` -- `PrecisionBiomarkerSettings` extends Pydantic
`BaseSettings` with `env_prefix="BIOMARKER_"`.

### 8.1 Environment Variables

All settings can be overridden via environment variables prefixed with
`BIOMARKER_`. The `.env` file is also loaded automatically.

| Variable                        | Default            | Description                     |
| ------------------------------- | ------------------ | ------------------------------- |
| `BIOMARKER_MILVUS_HOST`        | `localhost`        | Milvus gRPC host               |
| `BIOMARKER_MILVUS_PORT`        | `19530`            | Milvus gRPC port               |
| `BIOMARKER_API_HOST`           | `0.0.0.0`          | API bind address                |
| `BIOMARKER_API_PORT`           | `8529`             | FastAPI port                    |
| `BIOMARKER_STREAMLIT_PORT`     | `8528`             | Streamlit port                  |
| `ANTHROPIC_API_KEY`            | (none)             | Claude API key (required)       |
| `BIOMARKER_LLM_PROVIDER`      | `anthropic`        | LLM provider                    |
| `BIOMARKER_LLM_MODEL`         | `claude-sonnet-4-6` | Claude model ID             |
| `BIOMARKER_EMBEDDING_MODEL`   | `BAAI/bge-small-en-v1.5` | Embedding model          |
| `BIOMARKER_EMBEDDING_DIMENSION`| `384`             | Embedding dimensions            |
| `BIOMARKER_EMBEDDING_BATCH_SIZE`| `32`             | Embedding batch size            |
| `BIOMARKER_TOP_K_PER_COLLECTION`| `5`              | Results per collection          |
| `BIOMARKER_SCORE_THRESHOLD`    | `0.4`              | Minimum cosine similarity       |
| `BIOMARKER_REQUEST_TIMEOUT_SECONDS` | `60`          | HTTP request timeout            |
| `BIOMARKER_MILVUS_TIMEOUT_SECONDS`  | `10`          | Milvus operation timeout        |
| `BIOMARKER_LLM_MAX_RETRIES`   | `3`                | LLM retry count                 |
| `BIOMARKER_MAX_CONVERSATION_CONTEXT` | `3`          | Conversation memory depth       |
| `BIOMARKER_CITATION_HIGH_THRESHOLD`  | `0.75`       | High-confidence citation cutoff |
| `BIOMARKER_CITATION_MEDIUM_THRESHOLD`| `0.60`       | Medium-confidence cutoff        |
| `BIOMARKER_CORS_ORIGINS`      | `http://localhost:8080,...` | CORS allowed origins     |
| `BIOMARKER_MAX_REQUEST_SIZE_MB`| `10`              | Maximum request body size       |
| `BIOMARKER_API_KEY`           | (empty)             | API auth key (empty = no auth)  |
| `BIOMARKER_METRICS_ENABLED`   | `true`              | Enable Prometheus metrics       |

### 8.2 Weight Overrides

Collection weights are individually overridable:

```bash
BIOMARKER_WEIGHT_BIOMARKER_REF=0.15
BIOMARKER_WEIGHT_GENETIC_VARIANTS=0.12
# ... etc.
```

A `@model_validator` warns if weights do not sum to approximately 1.0.

---

## 9. Embedding Strategy

| Parameter     | Value                   |
| ------------- | ----------------------- |
| Model         | BAAI/bge-small-en-v1.5  |
| Parameters    | 33M                     |
| Dimensions    | 384                     |
| Metric        | COSINE                  |
| Index type    | IVF_FLAT                |
| Batch size    | 32                      |
| Runtime       | CPU (no GPU required)   |

Each collection model provides a `to_embedding_text()` method that constructs
a domain-optimized text representation. Examples:

- **BiomarkerReference:** `"{name} ({unit}). {text_chunk}. Significance: {clinical_significance}. Category: {category}. Genetic modifiers: {genetic_modifiers}"`
- **PGxRule:** `"{gene} {star_alleles} -- {drug}. {text_chunk}. Recommendation: {recommendation}. Phenotype: {phenotype}. CPIC Level: {cpic_level}"`
- **GeneticVariant:** `"{gene} {rs_id}. {text_chunk}. Mechanism: {mechanism}. Diseases: {disease_associations}"`

The `seed_all.py` script reads JSON seed files, calls `to_embedding_text()` on
each record, embeds with SentenceTransformer, and inserts into Milvus.

---

## 10. Pharmacogenomic Engine

File: `src/pharmacogenomics.py` (1,503 lines)

Maps star alleles and genotypes to drug recommendations following CPIC Level 1A
guidelines. Pure computation -- no LLM or database calls required.

### 10.1 Covered Pharmacogenes (14)

| Gene        | Description                                | CPIC Level | Key Drugs                    |
| ----------- | ------------------------------------------ | ---------- | ---------------------------- |
| CYP2D6      | Metabolizes ~25% of drugs                  | 1A         | Codeine, tramadol, tamoxifen |
| CYP2C19     | Clopidogrel, PPIs, antidepressants         | 1A         | Clopidogrel, omeprazole      |
| CYP2C9      | Warfarin, NSAIDs, phenytoin                | 1A         | Warfarin, celecoxib          |
| CYP3A5      | Tacrolimus metabolism                       | 1A         | Tacrolimus                   |
| SLCO1B1     | Statin hepatic uptake transporter          | 1A         | Simvastatin, atorvastatin    |
| VKORC1      | Warfarin target sensitivity                | 1A         | Warfarin                     |
| MTHFR       | Folate metabolism, homocysteine             | Informational | Methotrexate             |
| TPMT        | Thiopurine metabolism                       | 1A         | Azathioprine, mercaptopurine |
| DPYD        | Fluoropyrimidine metabolism                 | 1A         | 5-FU, capecitabine           |

### 10.2 Phenotype Classification

Uses CPIC standard metabolizer terminology:
- CYP enzymes: Normal / Intermediate / Poor / Ultra-rapid / Rapid Metabolizer
- SLCO1B1: Normal / Intermediate / Poor Function (transporter activity)
- MTHFR: Normal / Intermediate / Reduced Activity (enzyme activity)
- VKORC1: Normal / Intermediate / High Sensitivity (drug sensitivity)

### 10.3 Drug Recommendation Actions

Each drug recommendation includes an `action` category:

| Action                | Meaning                                    |
| --------------------- | ------------------------------------------ |
| `STANDARD_DOSING`     | No change needed                           |
| `DOSE_REDUCTION`      | Reduce dose per recommendation             |
| `DOSE_ADJUSTMENT`     | Adjust dose based on context               |
| `CONSIDER_ALTERNATIVE`| Current drug may work, alternative preferred|
| `AVOID`               | Do not use this drug                       |
| `CONTRAINDICATED`     | Absolute contraindication (FDA/EMA)        |

Alert levels: `INFO` (routine), `WARNING` (clinical review), `CRITICAL`
(immediate action).

### 10.4 Usage

```python
from src.pharmacogenomics import PharmacogenomicMapper

mapper = PharmacogenomicMapper()
result = mapper.map_star_alleles({"CYP2D6": "*4/*4", "CYP2C19": "*1/*2"})
# Returns: PGxResult with phenotypes, drug recommendations, alerts
```

---

## 11. Biological Age Engine

File: `src/biological_age.py` (408 lines)

Implements PhenoAge (Levine et al. 2018, PMID:29676998) and GrimAge surrogate
estimation from standard blood biomarkers. Pure computation.

### 11.1 PhenoAge Algorithm

Uses 9 blood biomarkers with published coefficients from the dayoonkwon/BioAge
R package:

| Biomarker               | Coefficient | Unit (input) | Unit (internal SI) |
| ------------------------ | ----------- | ------------ | ------------------ |
| Albumin                  | -0.0336     | g/dL         | g/L (x10)          |
| Creatinine               | +0.0095     | mg/dL        | umol/L (x88.4)     |
| Glucose                  | +0.1953     | mg/dL        | mmol/L (/18.016)   |
| ln(CRP)                  | +0.0954     | mg/L         | ln(mg/L)           |
| Lymphocyte %             | -0.0120     | %            | %                  |
| MCV                      | +0.0268     | fL           | fL                 |
| RDW                      | +0.3306     | %            | %                  |
| Alkaline Phosphatase     | +0.0019     | U/L          | U/L                |
| WBC                      | +0.0554     | 10^3/uL      | 10^3/uL            |

Intercept: -19.9067. The module accepts standard US clinical units and converts
internally before applying coefficients.

### 11.2 Gompertz Mortality Model

Converts the linear predictor to biological age via Gompertz mortality
parameters:
- Mortality numerator: -1.51714
- Mortality denominator (gamma): 0.007692696
- BA intercept: 141.50225

Standard error: 4.9 years (from NHANES III validation).

### 11.3 GrimAge Surrogate

Uses available blood biomarkers to approximate GrimAge components when DNA
methylation data is unavailable. Returns estimated biological age with
confidence interval.

### 11.4 Usage

```python
from src.biological_age import BiologicalAgeCalculator

calc = BiologicalAgeCalculator()
result = calc.calculate(
    chronological_age=45,
    biomarkers={
        "albumin": 4.2,        # g/dL
        "creatinine": 0.9,     # mg/dL
        "glucose": 95,         # mg/dL
        "hs_crp": 1.2,        # mg/L
        "lymphocyte_pct": 30,  # %
        "mcv": 88,             # fL
        "rdw": 13.2,           # %
        "alkaline_phosphatase": 65,  # U/L
        "wbc": 6.5,            # 10^3/uL
    },
)
# Returns: BiologicalAgeResult with phenoage, grimage, delta, interpretation
```

---

## 12. Disease Trajectory Engine

File: `src/disease_trajectory.py` (1,421 lines)

Detects pre-symptomatic disease trajectories across 9 categories using
genotype-stratified biomarker thresholds. Pure computation.

### 12.1 Disease Categories (9)

| Category         | Display Name               | Key Biomarkers                        | Genetic Modifiers                |
| ---------------- | -------------------------- | ------------------------------------- | -------------------------------- |
| type2_diabetes   | Type 2 Diabetes            | HbA1c, fasting glucose, insulin, HOMA-IR | TCF7L2, PPARG, SLC30A8, KCNJ11, GCKR |
| cardiovascular   | Cardiovascular Disease     | Lp(a), LDL, ApoB, hs-CRP, HDL, TG   | APOE, PCSK9, LPA, IL6          |
| liver            | Liver Disease (NAFLD)      | ALT, AST, GGT, ferritin, platelets   | PNPLA3, TM6SF2, HSD17B13       |
| thyroid          | Thyroid Dysfunction        | TSH, free T4, free T3                | DIO2, DIO1                      |
| iron             | Iron Metabolism Disorder   | Ferritin, transferrin sat, serum iron | HFE C282Y, HFE H63D            |
| nutritional      | Nutritional Deficiency     | Omega-3, Vit D, B12, folate, Mg, Zn  | FADS1, FADS2, VDR, BCMO1, FUT2 |
| kidney           | Kidney Disease             | eGFR, cystatin C, BUN, urine ACR     | APOL1                           |
| bone_health      | Bone Health                | Calcium, PTH, Vit D, CTX, P1NP       | VDR, COL1A1                     |
| cognitive        | Cognitive Decline          | Homocysteine, B12, folate, hs-CRP    | APOE, MTHFR                    |

### 12.2 Genotype-Stratified Thresholds

Thresholds are adjusted based on genotype via `GENOTYPE_THRESHOLDS` in the
knowledge module:

```python
# TCF7L2 risk alleles lower the HbA1c threshold for concern:
#   0 risk alleles: 6.0%
#   1 risk allele:  5.8%
#   2 risk alleles: 5.5%

# PNPLA3 I148M genotype adjusts ALT upper limit:
#   CC (wild-type): 56 U/L
#   CG (heterozygous): 45 U/L
#   GG (homozygous): 35 U/L
```

### 12.3 Stage Progression

Each disease defines ordered stages. The engine identifies the current stage,
estimates years to clinical diagnosis, and calculates the intervention window
with potential risk reduction percentage.

---

## 13. Genotype Adjustment Engine

File: `src/genotype_adjustment.py` (1,225 lines)

Adjusts standard biomarker reference ranges based on individual genotype, age,
sex, and ancestry. Pure computation.

### 13.1 Modifier Genes (7)

| Gene    | Affected Biomarkers              | Adjustment Mechanism            |
| ------- | -------------------------------- | ------------------------------- |
| MTHFR   | Homocysteine, folate             | Reduced enzyme activity         |
| APOE    | LDL, total cholesterol           | Lipid metabolism variation      |
| PNPLA3  | ALT, AST                         | Hepatic lipid accumulation      |
| HFE     | Ferritin, transferrin saturation | Iron absorption dysregulation   |
| DIO2    | TSH, free T4                     | Impaired T4-to-T3 conversion   |
| VDR     | Vitamin D, calcium               | Vitamin D receptor sensitivity  |
| FADS1   | Omega-3 index                    | Fatty acid desaturation         |

### 13.2 Age-Stratified Brackets

Five age brackets for reference range stratification:
- 0-17, 18-39, 40-59, 60-79, 80+

Sources: NHANES III, Framingham Heart Study, KDIGO 2012, ATA 2017, ADA 2024,
ACC/AHA 2019, Endocrine Society guidelines.

### 13.3 Ancestry-Aware Adjustments

Population-specific biomarker adjustments from the knowledge module
(`ANCESTRY_ADJUSTMENTS`):
- **African:** Higher Lp(a) prevalence, lower TG, higher creatinine
- **South Asian:** Lower LDL/HbA1c thresholds due to higher CVD risk
- **East Asian:** Lower ALT limits, statin sensitivity
- **Hispanic:** Higher NAFLD prevalence, earlier diabetes screening

---

## 14. Critical Value Engine

File: `src/critical_values.py` (179 lines)

Evaluates patient biomarker values against 21 critical threshold rules loaded
from `biomarker_critical_values.json`. Real-time alerting for life-threatening
lab values.

### 14.1 Alert Severity Levels

| Level    | Meaning                           | Response Time     |
| -------- | --------------------------------- | ----------------- |
| CRITICAL | Immediately life-threatening      | Immediate         |
| URGENT   | Requires prompt clinical action   | Within hours      |
| WARNING  | Clinically significant deviation  | Next visit        |

### 14.2 Covered Biomarkers (15)

Platelet Count, Glucose, Potassium, INR, Sodium, Hemoglobin, Calcium,
Troponin I, WBC Count, Creatinine, Total Bilirubin, eGFR, Free T4, TSH,
Lactate.

Each rule includes: critical high/low thresholds, severity level, escalation
target (e.g., "Emergency Physician"), clinical action, cross-check biomarkers,
and LOINC code.

### 14.3 Biomarker Alias Resolution

The engine resolves multiple input names to canonical forms. For example,
`"platelets"`, `"platelet_count"`, and `"plt"` all resolve to `"Platelet Count"`.

---

## 15. Discordance Detector

File: `src/discordance_detector.py` (299 lines)

Detects cross-biomarker discordance patterns from
`biomarker_discordance_rules.json`. Identifies contradictory or unexpected
relationships between biomarker pairs.

### 15.1 Detection Output

Each `DiscordanceFinding` includes:
- `rule_name`: Descriptive rule name
- `biomarker_a`, `biomarker_b`: The discordant biomarker pair
- `condition`: The specific discordance condition detected
- `differential_diagnosis`: Possible explanations (list)
- `agent_handoff`: Recommended agents for follow-up (list)
- `priority`: high, moderate, low

---

## 16. Knowledge Graph

File: `src/knowledge.py` (1,326 lines)

Contains 4 major knowledge structures, versioned clinical thresholds, age/sex
reference ranges, ancestry adjustments, and plausible validation ranges.

### 16.1 BIOMARKER_DOMAINS (6 disease domains)

Each domain contains:
- Key biomarkers with units, normal/pre-disease/disease ranges, clinical notes
- Genetic modifiers with risk alleles and molecular mechanisms
- Intervention targets with evidence-based recommendations

Domains: Diabetes/Metabolic, Cardiovascular, Liver (NAFLD/Fibrosis), Thyroid,
Iron Metabolism, Nutritional Deficiency.

### 16.2 PGX_KNOWLEDGE (14 pharmacogenes)

Maps each pharmacogene to key drugs and CPIC guidance. Used by the RAG system
prompt and agent reasoning.

### 16.3 PHENOAGE_KNOWLEDGE

PhenoAge clock biomarker descriptions, coefficients, and clinical
interpretation. Used for biological age context in LLM responses.

### 16.4 CROSS_MODAL_LINKS (8 links)

Maps biomarker findings to triggers for other HCLS AI Factory agents:
- Elevated Lp(a) -> Imaging Agent (coronary calcium scoring)
- Iron overload -> Imaging Agent (liver MRI)
- PGx drug safety -> CAR-T/Oncology Agent
- VCF re-analysis -> Genomics Pipeline

### 16.5 Shared Clinical Thresholds

`GENOTYPE_THRESHOLDS` dictionary provides genotype-stratified thresholds used
by both the disease trajectory and genotype adjustment engines for consistency.

### 16.6 Age-Sex Reference Ranges

`AGE_SEX_REFERENCE_RANGES` provides clinically validated ranges stratified by
sex and age bracket (18-49, 50-69, 70+) for: creatinine, ALT, alkaline
phosphatase, ferritin, TSH, hemoglobin, BUN, cystatin C, homocysteine,
vitamin D.

### 16.7 Biomarker Plausible Ranges

`BIOMARKER_PLAUSIBLE_RANGES` provides validation bounds for 35+ biomarkers to
detect likely data entry errors before analysis.

### 16.8 Ancestry Adjustments

`ANCESTRY_ADJUSTMENTS` provides population-specific threshold multipliers for
African, South Asian, East Asian, and Hispanic populations with PMID citations.

### 16.9 Knowledge Versioning

```python
KNOWLEDGE_VERSION = {
    "version": "1.0.0",
    "cpic_version": "March 2025",
    "ada_standards": "2024",
    "esc_guidelines": "2021",
    "aasld_guidelines": "2023",
    "levine_phenoage": "2018",
    "lu_grimage": "2019",
}
```

---

## 17. RAG Engine

File: `src/rag_engine.py` (573 lines)

Multi-collection RAG engine that searches across all 14 collections
simultaneously using parallel `ThreadPoolExecutor`, applies weighted scoring,
and generates grounded LLM responses.

### 17.1 Search Pipeline

```
Query text
    |
    v
[Embed] -> BGE-small-en-v1.5 (384-dim vector)
    |
    v
[Parallel Search] -> ThreadPoolExecutor across 14 collections
    |                  Each returns top_k=5 results
    v
[Score & Weight] -> cosine_similarity * collection_weight
    |                 Filter by SCORE_THRESHOLD (0.4)
    v
[Merge & Rank] -> CrossCollectionResult (sorted by weighted score)
    |
    v
[Knowledge Augment] -> Inject domain knowledge, PGx context, PhenoAge info
    |
    v
[LLM Synthesis] -> Claude generates response with citations
```

### 17.2 Input Sanitization

Filter expressions are validated against `_SAFE_FILTER_RE = r"^[A-Za-z0-9 _\-]+$"`
to prevent Milvus filter injection.

### 17.3 System Prompt

The `BIOMARKER_SYSTEM_PROMPT` (75 lines) instructs Claude to:
- Cite evidence using collection labels (e.g., `[BiomarkerRef:marker-id]`)
- Always specify units when discussing values
- Provide genotype-specific interpretation
- Highlight critical findings (PGx drug safety, severe iron overload)
- Recommend actionable interventions grounded in CPIC evidence
- Explain pre-symptomatic disease trajectories with timelines
- Flag cross-modal triggers for other agents

---

## 18. Agent Orchestrator

File: `src/agent.py` (610 lines)

`PrecisionBiomarkerAgent` implements the plan-analyze-search-synthesize-report
pattern from the VAST AI OS AgentEngine model.

### 18.1 Initialization

```python
agent = PrecisionBiomarkerAgent(
    rag_engine=rag_engine,
    bio_age_calc=BiologicalAgeCalculator(),
    trajectory_analyzer=DiseaseTrajectoryAnalyzer(),
    pgx_mapper=PharmacogenomicMapper(),
    genotype_adjuster=GenotypeAdjuster(),
)
```

### 18.2 Agent Loop

1. **Plan** (`search_plan()`): Parse question, identify topics, disease areas,
   relevant analysis modules, and formulate sub-questions
2. **Analyze** (`analyze_patient()`): Run all applicable analysis modules:
   - BiologicalAgeCalculator
   - DiseaseTrajectoryAnalyzer
   - PharmacogenomicMapper
   - GenotypeAdjuster
   - CriticalValueEngine
   - DiscordanceDetector
   - LabRangeInterpreter
3. **Search** (`rag_engine.retrieve()`): Multi-collection RAG search
4. **Synthesize** (`evaluate_evidence()`): Merge analysis results with RAG
   evidence, generate LLM response with citations
5. **Report**: Pass to `ReportGenerator` for structured output

### 18.3 SearchPlan Dataclass

```python
@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str]
    disease_areas: List[str]
    relevant_modules: List[str]      # e.g., ["biological_age", "pgx", "trajectory"]
    search_strategy: str             # "broad", "targeted", "domain-specific"
    sub_questions: List[str]
```

### 18.4 Cross-Modal Triggers

The agent generates cross-modal triggers when findings warrant follow-up by
other HCLS AI Factory agents (imaging, oncology, genomics pipeline).

---

## 19. Export Pipeline

File: `src/export.py` (1,392 lines)

Exports analysis results in 5 formats:

| Format   | Function                          | Description                     |
| -------- | --------------------------------- | ------------------------------- |
| Markdown | `export_markdown()`              | Human-readable with evidence tables |
| JSON     | `export_json()`                  | Machine-readable structured data |
| PDF      | `export_pdf()`                   | Styled report via reportlab Platypus |
| CSV      | `export_csv()`                   | Tabular export for spreadsheets |
| FHIR R4  | `export_fhir_diagnostic_report()`| FHIR R4 DiagnosticReport JSON bundle |

### 19.1 PDF Generation

Uses `reportlab` Platypus engine with styled tables, headers, risk-level color
coding, and embedded charts. Generates timestamped filenames:
`biomarker_report_20260301T143025Z_a1b2.pdf`

### 19.2 FHIR R4 Export

Generates a compliant FHIR R4 `DiagnosticReport` resource bundle including:
- Patient reference
- Observation resources for each biomarker
- Condition resources for disease trajectories
- MedicationStatement resources for PGx recommendations

---

## 20. Report Generator

File: `src/report_generator.py` (993 lines)

`ReportGenerator` class produces structured 12-section clinical reports in
markdown format from an `AnalysisResult`.

### 20.1 Report Sections (12)

| #  | Section                         | Content                               |
| -- | ------------------------------- | ------------------------------------- |
| 1  | Biological Age Assessment       | PhenoAge/GrimAge, age acceleration    |
| 2  | Executive Findings              | Top 5 critical/high priority findings |
| 3  | Biomarker-Gene Correlation Map  | Genotype-biomarker interaction matrix |
| 4  | Disease Trajectory Analysis     | 9 disease categories with stages      |
| 5  | Pharmacogenomic Profile         | PGx phenotypes and drug recommendations |
| 6  | Nutritional Analysis            | Genotype-guided supplement protocols  |
| 7  | Interconnected Pathways         | Cross-domain pathway connections      |
| 8  | Prioritized Action Plan         | Ranked interventions by urgency       |
| 9  | Monitoring Schedule             | Follow-up timeline and biomarkers     |
| 10 | Supplement Protocol Summary     | Dosing schedules based on genotype    |
| 11 | Clinical Summary for MD         | Physician-oriented executive summary  |
| 12 | References                      | PMID citations and guideline versions |

Additional sections: Evidence Provenance, Clinical Validation.

### 20.2 Usage

```python
from src.report_generator import ReportGenerator

generator = ReportGenerator()
markdown = generator.generate(analysis_result, patient_profile)
```

---

## 21. FastAPI REST Server

File: `api/main.py` (465 lines) + route modules in `api/routes/`.

### 21.1 Application Setup

- Title: "Biomarker Intelligence Agent API"
- Version: 1.0.0
- Docs: `/docs` (Swagger UI), `/openapi.json`
- CORS: Configurable origins via `BIOMARKER_CORS_ORIGINS`
- Auth: Optional API key via `X-API-Key` header (skips `/health` and `/metrics`)
- Request size limit: 10 MB (configurable)

### 21.2 Lifespan Initialization

On startup, the lifespan handler initializes:
1. `BiomarkerCollectionManager` (Milvus connection)
2. `SentenceTransformer` embedder (BGE-small-en-v1.5)
3. `Anthropic` LLM client
4. All analysis modules (bio age, trajectory, PGx, genotype adjuster)
5. `BiomarkerRAGEngine` (wires everything together)
6. `PrecisionBiomarkerAgent` (autonomous reasoning)

On shutdown, Milvus connection is closed.

### 21.3 API Endpoints (19+)

**Core endpoints** (defined in `api/main.py`):

| Method | Path              | Description                          |
| ------ | ----------------- | ------------------------------------ |
| GET    | `/`               | Service info (name, docs, health)    |
| GET    | `/health`         | Collection count, vector count, agent status |
| GET    | `/collections`    | Collection names and record counts   |
| GET    | `/knowledge/stats`| Knowledge graph statistics           |
| GET    | `/metrics`        | Prometheus-compatible metrics        |

**Analysis router** (`/v1/`, defined in `api/routes/analysis.py`):

| Method | Path                  | Description                          |
| ------ | --------------------- | ------------------------------------ |
| POST   | `/v1/analyze`         | Full patient analysis (all modules)  |
| POST   | `/v1/biological-age`  | Biological age calculation           |
| POST   | `/v1/disease-risk`    | Disease trajectory analysis          |
| POST   | `/v1/pgx`             | Pharmacogenomic mapping              |
| POST   | `/v1/query`           | RAG Q&A query                        |
| POST   | `/v1/query/stream`    | Streaming RAG Q&A (SSE)             |
| GET    | `/v1/health`          | V1 router health check              |

**Events router** (`/v1/events/`, defined in `api/routes/events.py`):

| Method | Path                        | Description                      |
| ------ | --------------------------- | -------------------------------- |
| POST   | `/v1/events/cross-modal`    | Receive cross-modal event        |
| POST   | `/v1/events/biomarker-alert`| Send biomarker alert             |
| GET    | `/v1/events/cross-modal`    | List inbound events              |
| GET    | `/v1/events/biomarker-alert`| List outbound alerts             |

**Reports router** (`/v1/report/`, defined in `api/routes/reports.py`):

| Method | Path                          | Description                    |
| ------ | ----------------------------- | ------------------------------ |
| POST   | `/v1/report/generate`         | Generate 12-section report     |
| GET    | `/v1/report/{report_id}/pdf`  | Download report as PDF         |
| POST   | `/v1/report/fhir`             | Export as FHIR R4 bundle       |

### 21.4 Request Validation

All POST endpoints use Pydantic request models with:
- Field constraints (`ge=0`, `le=150`, `pattern="^(M|F)$"`)
- `@model_validator` requiring at least one data source (biomarkers, genotypes,
  or star alleles)

### 21.5 HIPAA Audit Logging

Analysis and report endpoints call `audit_log()` with action type and hashed
patient ID for every request.

---

## 22. Streamlit UI

File: `app/biomarker_ui.py` (1,863 lines)

### 22.1 Tabs (8)

| Tab               | Function                                           |
| ----------------- | -------------------------------------------------- |
| Biomarker Analysis| Full patient analysis with sample patient quick-load |
| Biological Age    | PhenoAge/GrimAge calculator with visualization      |
| Disease Risk      | Focused disease trajectory analysis                  |
| PGx Profile       | Pharmacogenomic drug interaction mapping             |
| Evidence Explorer | RAG Q&A with collection filtering                    |
| Reports           | PDF and FHIR R4 report generation and download       |
| Patient 360       | Unified cross-agent intelligence dashboard           |
| Longitudinal      | Biomarker trend tracking across multiple visits      |

### 22.2 Page Configuration

```python
st.set_page_config(
    page_title="Biomarker Intelligence Agent -- HCLS AI Factory",
    page_icon="🩸",
    layout="wide",
    initial_sidebar_state="expanded",
)
```

### 22.3 Engine Initialization

Uses `@st.cache_resource(ttl=300)` to cache the analysis engine across Streamlit
reruns (5-minute TTL). Initializes: BiomarkerCollectionManager,
BiologicalAgeCalculator, DiseaseTrajectoryAnalyzer, CriticalValueEngine,
DiscordanceDetector, LabRangeInterpreter, PharmacogenomicMapper.

### 22.4 Patient 360 Dashboard

File: `app/patient_360.py` (670 lines)

Cross-agent intelligence dashboard that aggregates findings from multiple HCLS
AI Factory agents into a unified patient view. Includes:
- Biomarker summary with risk indicators
- PGx drug interaction alerts
- Disease trajectory timelines
- Cross-agent trigger status
- Longitudinal trend charts

---

## 23. Demo Patients

Two pre-configured demo patients for validation and demonstration:

### 23.1 Patient HCLS-BIO-2026-00001

| Field         | Value                                      |
| ------------- | ------------------------------------------ |
| Patient ID    | HCLS-BIO-2026-00001                        |
| Age           | 45                                         |
| Sex           | Male                                       |
| Ethnicity     | Ashkenazi Jewish                           |
| Genome        | HG002 / NA24385 (Genome in a Bottle)       |
| Key genotypes | MTHFR C677T, APOE E3/E4, PNPLA3 I148M     |
| PGx alleles   | CYP2D6 *1/*4, CYP2C19 *1/*2               |

### 23.2 Patient HCLS-BIO-2026-00002

| Field         | Value                                      |
| ------------- | ------------------------------------------ |
| Patient ID    | HCLS-BIO-2026-00002                        |
| Age           | 38                                         |
| Sex           | Female                                     |
| Ethnicity     | Ashkenazi Jewish                           |

Demo patient data is stored in `data/reference/biomarker_sample_patients.json`
and can be loaded via the Streamlit UI "Biomarker Analysis" tab quick-load
buttons.

---

## 24. Cross-Agent Integration

The Biomarker Intelligence Agent communicates with other HCLS AI Factory agents
via the cross-modal event system and `cross_modal/cross_agent.py`. It calls 4
peer agents directly (Oncology, CAR-T, PGx, Clinical Trial) and the
`/integrated-assessment` endpoint for multi-agent clinical synthesis. The
platform comprises 11 agents total within the Precision Intelligence Network.

### 24.1 Outbound Triggers (Biomarker -> Other Agents)

| Finding                    | Target Agent       | Trigger                        |
| -------------------------- | ------------------ | ------------------------------ |
| Elevated Lp(a) > 125 nmol/L | Imaging Agent    | Coronary calcium scoring       |
| Iron overload (ferritin > 1000) | Imaging Agent | Liver MRI (T2*)               |
| DPYD Poor Metabolizer     | Oncology Agent     | 5-FU contraindication alert    |
| BRCA1/2 carrier (AJ panel)| Oncology Agent     | Cancer risk assessment         |
| Novel variant detected     | Genomics Pipeline  | VCF re-analysis trigger        |
| GBA + APOE E4 compound    | Autoimmune Agent   | Parkinson's risk assessment    |
| CD19/CD22 expression change| CAR-T Agent       | CAR-T resistance assessment    |
| MRD status change          | Oncology Agent     | Risk re-stratification         |
| PGx interaction detected   | PGx Agent          | Extended interaction analysis  |
| Biomarker-driven eligibility| Clinical Trial Agent | Trial matching trigger       |

### 24.2 Inbound Events (Other Agents -> Biomarker)

| Source Agent         | Event Type        | Action                           |
| -------------------- | ----------------- | -------------------------------- |
| Imaging Agent        | Imaging finding   | Correlate with biomarker patterns |
| Genomics Pipeline    | Genomic variant   | Update PGx profile, re-analyze   |
| Oncology Agent       | Drug alert        | Check PGx interactions           |

### 24.3 Event API

Events are exchanged via `/v1/events/cross-modal` (POST to send, GET to list).
In-memory store with configurable maximum (1,000 events). Production deployments
would use a message bus (e.g., Redis, Kafka).

---

## 25. Monitoring and Metrics

### 25.1 Prometheus Metrics

Endpoint: `GET /metrics` (Prometheus text format)

| Metric                              | Type    | Description                    |
| ----------------------------------- | ------- | ------------------------------ |
| `biomarker_api_requests_total`      | counter | Total API requests             |
| `biomarker_api_query_requests_total`| counter | Total /query requests          |
| `biomarker_api_search_requests_total`| counter| Total /search requests         |
| `biomarker_api_analyze_requests_total`| counter| Total /analyze requests       |
| `biomarker_api_bio_age_requests_total`| counter| Total /biological-age requests |
| `biomarker_api_errors_total`        | counter | Total error responses          |
| `biomarker_collection_vectors`      | gauge   | Vectors per collection (labeled) |

### 25.2 Health Check

Endpoint: `GET /health`

Returns: `{ status, collections, total_vectors, agent_ready }`
- `healthy`: All systems operational
- `degraded`: Milvus unavailable but API responsive

Docker health check: Every 30 seconds via Python stdlib urllib (no curl
dependency in container).

---

## 26. Audit Logging

File: `src/audit.py` (83 lines)

HIPAA-compliant audit logging for all patient data access. Uses loguru with a
dedicated `audit=True` binding for separate audit event routing.

### 26.1 Auditable Actions

| Action                 | When Logged                              |
| ---------------------- | ---------------------------------------- |
| `PATIENT_ANALYSIS`     | Full patient analysis request            |
| `BIOLOGICAL_AGE`       | Biological age calculation               |
| `DISEASE_RISK`         | Disease trajectory analysis              |
| `PGX_MAPPING`          | Pharmacogenomic mapping                  |
| `RAG_QUERY`            | RAG evidence query                       |
| `REPORT_GENERATED`     | Report generation                        |
| `REPORT_EXPORTED`      | Report downloaded (PDF/CSV)              |
| `FHIR_EXPORTED`        | FHIR R4 export                           |
| `PATIENT_DATA_ACCESSED`| Any patient data access                  |

### 26.2 Privacy

Patient IDs are SHA-256 hashed before logging. Full IDs are only stored in
encrypted production storage, never in log files.

---

## 27. Testing

**709 tests across 18 test files, all passing.**

### 27.1 Test Files

| File                         | Module Under Test          |
| ---------------------------- | -------------------------- |
| `test_models.py`            | Pydantic model validation  |
| `test_collections.py`       | Milvus collection manager  |
| `test_rag_engine.py`        | RAG engine search/synthesis |
| `test_agent.py`             | Agent orchestrator         |
| `test_pharmacogenomics.py`  | PGx phenotyping            |
| `test_biological_age.py`    | PhenoAge/GrimAge           |
| `test_disease_trajectory.py`| Disease trajectory engine  |
| `test_genotype_adjustment.py`| Genotype adjustments      |
| `test_critical_values.py`   | Critical value alerts      |
| `test_discordance_detector.py`| Discordance detection    |
| `test_lab_range_interpreter.py`| Lab range interpretation |
| `test_report_generator.py`  | Report generation          |
| `test_export.py`            | Export pipeline (MD/JSON/PDF/FHIR) |
| `test_api.py`               | FastAPI endpoint testing   |
| `test_integration.py`       | End-to-end integration     |
| `test_edge_cases.py`        | Edge cases and error handling |
| `test_longitudinal.py`      | Longitudinal tracking      |
| `test_ui.py`                | Streamlit UI components    |

### 27.2 Test Fixtures

Shared fixtures in `tests/conftest.py`:
- `mock_embedder`: Returns 384-dim zero vectors
- `mock_llm_client`: Returns "Mock response"
- `mock_collection_manager`: Simulates Milvus operations
- Sample patient profiles and biomarker data

All tests run without Milvus, Anthropic API, or external services.

### 27.3 Running Tests

```bash
# All tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=src --cov=api --cov-report=term-missing

# Single module
pytest tests/test_pharmacogenomics.py -v
```

---

## 28. Dependencies

File: `requirements.txt`

### 28.1 Core Dependencies

| Package              | Version Range     | Purpose                         |
| -------------------- | ----------------- | ------------------------------- |
| pydantic             | >=2.0,<3.0        | Data models and validation      |
| pydantic-settings    | >=2.7,<3.0        | Settings with env var support   |
| loguru               | >=0.7.0,<1.0      | Structured logging              |
| pymilvus             | >=2.4.0,<2.6      | Milvus vector database client   |
| sentence-transformers| >=2.2.0,<3.0      | BGE-small-en-v1.5 embedding     |
| anthropic            | >=0.18.0,<1.0     | Claude API client               |
| streamlit            | >=1.30.0,<2.0     | Web UI                          |
| fastapi              | >=0.109.0,<1.0    | REST API framework              |
| uvicorn[standard]    | >=0.27.0,<1.0     | ASGI server                     |
| python-multipart     | >=0.0.6,<1.0      | Form data parsing               |
| reportlab            | >=4.0.0,<5.0      | PDF generation                  |
| numpy                | >=1.24.0,<3.0     | Numerical computation           |
| pandas               | >=2.0.0,<3.0      | Data manipulation               |
| plotly                | >=5.18.0,<6.0     | Interactive charts              |
| tqdm                 | >=4.65.0,<5.0     | Progress bars                   |
| python-dotenv        | >=1.0.0,<2.0      | .env file loading               |
| prometheus-client    | >=0.20.0,<1.0     | Prometheus metrics              |

### 28.2 Testing Dependencies

| Package              | Version Range     | Purpose                         |
| -------------------- | ----------------- | ------------------------------- |
| pytest               | >=7.0,<8.0        | Test framework                  |
| pytest-asyncio       | >=0.21,<1.0       | Async test support              |
| pytest-cov           | >=4.0,<5.0        | Coverage reporting              |

---

## 29. Quick Start

### 29.1 Docker (Recommended)

```bash
cd ai_agent_adds/precision_biomarker_agent

# Configure API key
cp .env.example .env
# Edit .env: set ANTHROPIC_API_KEY=sk-ant-...

# Start all services
docker compose up -d

# Watch setup progress
docker compose logs -f biomarker-setup

# Access UI and API
open http://localhost:8528    # Streamlit UI
open http://localhost:8529/docs  # Swagger API docs
```

### 29.2 Local Development

```bash
cd ai_agent_adds/precision_biomarker_agent

# Create virtual environment
python -m venv .venv && source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Set API key
export ANTHROPIC_API_KEY=sk-ant-...

# Ensure Milvus is running (standalone or via Docker)
# docker run -d -p 19530:19530 milvusdb/milvus:v2.4-latest

# Create collections and seed data
python scripts/setup_collections.py --drop-existing
python scripts/seed_all.py

# Start Streamlit UI
streamlit run app/biomarker_ui.py --server.port 8528

# Start FastAPI server (separate terminal)
uvicorn api.main:app --host 0.0.0.0 --port 8529 --reload

# Run tests
pytest tests/ -v
```

### 29.3 Example API Call

```bash
curl -X POST http://localhost:8529/v1/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "age": 45,
    "sex": "M",
    "biomarkers": {
      "albumin": 4.2,
      "creatinine": 0.9,
      "glucose": 95,
      "hs_crp": 1.2,
      "hba1c": 5.8,
      "ldl": 142,
      "ferritin": 280,
      "alt": 38,
      "tsh": 2.8
    },
    "genotypes": {
      "rs1801133": "CT",
      "rs429358": "CT"
    },
    "star_alleles": {
      "CYP2D6": "*1/*4",
      "CYP2C19": "*1/*2"
    }
  }'
```

---

*Generated for the HCLS AI Factory. Import this document as context for Claude
Code sessions to get accurate, grounded implementation guidance.*
