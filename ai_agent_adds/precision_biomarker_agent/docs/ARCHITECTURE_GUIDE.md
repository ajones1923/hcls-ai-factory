# Precision Biomarker Intelligence Agent -- Architecture Guide

Architecture design document for the Precision Biomarker Intelligence Agent,
part of the HCLS AI Factory Precision Intelligence Network -- one of three
GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence
Network, Therapeutic Discovery Engine) that together deliver patient DNA to
drug candidates in under 5 hours on a single NVIDIA DGX Spark.

Author: Adam Jones
Date: March 2026

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [VAST AI OS Component Mapping](#vast-ai-os-component-mapping)
3. [System Architecture](#system-architecture)
4. [Component Architecture](#component-architecture)
5. [Data Flow](#data-flow)
6. [Milvus Collections](#milvus-collections)
7. [Embedding Strategy](#embedding-strategy)
8. [LLM Integration](#llm-integration)
9. [Knowledge Graph](#knowledge-graph)
10. [Clinical Pipelines](#clinical-pipelines)
11. [RAG Engine](#rag-engine)
12. [API Layer](#api-layer)
13. [UI Layer](#ui-layer)
14. [Export Pipeline](#export-pipeline)
15. [Scaling and Performance](#scaling-and-performance)
16. [Security](#security)
17. [File Structure](#file-structure)

---

## Executive Summary

The Precision Biomarker Intelligence Agent is a genotype-aware biomarker
interpretation platform that transforms standard blood panel results and
genomic data into actionable precision health intelligence. It combines
14 Milvus vector collections, 6 specialized clinical analysis engines,
and Claude LLM synthesis to deliver evidence-based recommendations
spanning pharmacogenomics, biological aging, disease trajectories,
genotype-adjusted reference ranges, critical value alerting, and
cross-biomarker discordance detection.

### Key Results

| Metric | Value |
|--------|-------|
| Source modules (src/) | 16 Python files |
| Total source lines | ~20,200 |
| Test suite | 709 tests across 18 files, all passing |
| Milvus collections | 14 (13 owned + 1 shared read-only) |
| Pharmacogenes mapped | 14 with CPIC-guided phenotyping |
| Disease trajectory categories | 9 |
| Genotype modifier genes | 7 |
| Critical value rules | 21 |
| Report sections | 12 |
| Export formats | 5 (Markdown, PDF, FHIR R4, CSV, JSON) |
| API endpoints | Multiple across 3 route modules |
| Demo patients | Included with sample data |
| Reference data files | 18 JSON seed files |
| Supported languages | 7 (English, Spanish, Chinese, Hindi, French, Arabic, Portuguese) |

---

## VAST AI OS Component Mapping

The Biomarker Agent maps to the VAST AI OS architecture as follows:

| VAST AI OS Layer | Agent Component |
|-----------------|-----------------|
| Data Layer | Milvus vector store (14 collections), 18 JSON reference files |
| Model Layer | BAAI/bge-small-en-v1.5 (embedding), Claude claude-sonnet-4-6 (synthesis) |
| Inference Layer | FastAPI server (8529), RAG engine, 6 clinical analysis engines |
| Application Layer | Streamlit UI (8528), 8-tab interface |
| Orchestration | Docker Compose (6 services), health checks, CORS middleware |
| Integration | Cross-agent stubs (Oncology, Autoimmune, CAR-T, Imaging), event bus |

---

## System Architecture

```
+-------------------------------------------------------------------+
|                    Streamlit UI (:8528)                            |
|  +-------+ +--------+ +--------+ +--------+ +--------+ +------+  |
|  |Clinical| |Patient | |PGx     | |Bio Age | |Disease | |Geno  |  |
|  | Query  | |Analysis| |Panel   | |Calc    | |Traject | |Adjust|  |
|  +-------+ +--------+ +--------+ +--------+ +--------+ +------+  |
|  +-------+ +--------+                                             |
|  |Critical| |Patient |                                            |
|  |Values  | | 360    |                                            |
|  +-------+ +--------+                                             |
+-------------------------------------------------------------------+
        |                                         |
        v                                         v
+-------------------+                 +----------------------------+
|  FastAPI (:8529)  |                 |  BiomarkerAgent            |
|  3 route modules  |                 |  plan()                    |
|  Auth / CORS      |                 |  analyze()                 |
|  Request limiter  |                 |  search()                  |
|  Prometheus /     |                 |  synthesize()              |
|    metrics        |                 |  report()                  |
+-------------------+                 +----------------------------+
        |                                     |
        v                                     v
+-------------------+            +----------------------------+
| BiomarkerRAG      |            | 6 Clinical Engines         |
| Engine            |            |                            |
|  retrieve()       |            | PharmacogenomicMapper      |
|  query()          |            |   14 pharmacogenes, CPIC   |
|  query_stream()   |            |                            |
|  _embed()         |            | BiologicalAgeCalculator    |
|  _build_context() |            |   PhenoAge + GrimAge       |
+-------------------+            |                            |
        |                        | DiseaseTrajectoryAnalyzer  |
        v                        |   9 disease categories     |
+-------------------+            |                            |
| CollectionManager |            | GenotypeAdjuster           |
|  connect()        |            |   7 modifier genes         |
|  search_all()     |            |                            |
|  insert_batch()   |            | CriticalValueEngine        |
+-------------------+            |   21 threshold rules       |
        |                        |                            |
        v                        | DiscordanceDetector        |
+-------------------+            |   cross-biomarker patterns |
| Milvus (:19530)   |            +----------------------------+
| 14 collections    |                    |
| IVF_FLAT / COSINE |                    v
| 384-dim vectors   |            +----------------------------+
+-------------------+            | ReportGenerator            |
        |                        |   12-section reports       |
        v                        | ExportPipeline             |
+-------------------+            |   MD/PDF/FHIR/CSV/JSON     |
| BGE-small-en-v1.5 |            | TranslationEngine          |
| Embedding Model   |            |   7 languages              |
| 384 dimensions    |            +----------------------------+
+-------------------+
```

---

## Component Architecture

The agent is organized into four layers with clear separation of concerns.

### Presentation Layer

- **Streamlit UI** (`app/biomarker_ui.py`, 1,863 lines): 8-tab application
  covering clinical query, patient analysis, pharmacogenomics, biological age,
  disease trajectory, genotype adjustments, critical values, and the cross-agent
  Patient 360 dashboard.
- **Patient 360 Dashboard** (`app/patient_360.py`, 670 lines): Cross-agent
  view aggregating findings from all five HCLS AI Factory intelligence agents.
- **Protein Viewer** (`app/protein_viewer.py`, 168 lines): 3D protein structure
  visualization for pharmacogenomically relevant variants.

### API Layer

- **FastAPI Server** (`api/main.py`, 465 lines): Entry point with lifespan
  management, core endpoints, health checks, and middleware configuration.
- **Analysis Routes** (`api/routes/analysis.py`, 495 lines): Endpoints for
  full analysis, biological age calculation, disease trajectory, PGx mapping,
  and genotype adjustments.
- **Event Routes** (`api/routes/events.py`, 326 lines): Cross-modal event
  endpoints and biomarker alert handling.
- **Report Routes** (`api/routes/reports.py`, 296 lines): Report generation,
  PDF export, and FHIR R4 bundle endpoints.

### Engine Layer

Six specialized clinical engines plus the RAG engine and agent orchestrator:

| Engine | Module | Lines | Purpose |
|--------|--------|-------|---------|
| PharmacogenomicMapper | `pharmacogenomics.py` | 1,503 | CPIC-guided phenotyping across 14 pharmacogenes |
| DiseaseTrajectoryAnalyzer | `disease_trajectory.py` | 1,421 | Pre-symptomatic detection across 9 disease categories |
| GenotypeAdjuster | `genotype_adjustment.py` | 1,225 | Genotype-based reference range adjustments for 7 genes |
| BiologicalAgeCalculator | `biological_age.py` | 408 | PhenoAge and GrimAge surrogate calculations |
| CriticalValueEngine | `critical_values.py` | 179 | Life-threatening threshold detection (21 rules) |
| DiscordanceDetector | `discordance_detector.py` | 299 | Cross-biomarker anomaly detection |
| BiomarkerRAGEngine | `rag_engine.py` | 573 | Multi-collection RAG with parallel search |
| BiomarkerAgent | `agent.py` | 610 | Autonomous plan-analyze-search-synthesize loop |

### Data Layer

- **Milvus** (14 vector collections, 384-dim BGE-small-en-v1.5 embeddings)
- **Knowledge Graph** (`knowledge.py`, 1,326 lines): 6 disease domain knowledge graphs
- **Reference Data** (18 JSON seed files in `data/reference/`)
- **Models** (`models.py`, 786 lines): 14 collection models, 8+ analysis models, 7 enums

---

## Data Flow

### End-to-End Pipeline

```
                    HCLS AI Factory Pipeline
                    =======================

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

### Agent Internal Pipeline

The agent orchestrator (`BiomarkerAgent`) executes a 5-step pipeline for
each patient analysis request:

```
Patient Profile (demographics, biomarkers, genotypes)
    |
    v
[1. PLAN]
    |  Identify relevant topics, disease areas, and modules
    |  Select analysis engines based on available patient data
    |  Determine which collections to prioritize
    |
    v
[2. ANALYZE]  (parallel execution)
    |  +-- BiologicalAgeCalculator: PhenoAge + GrimAge computation
    |  +-- DiseaseTrajectoryAnalyzer: 9 disease category risk assessment
    |  +-- PharmacogenomicMapper: 14 pharmacogene phenotyping
    |  +-- GenotypeAdjuster: 7 modifier gene range adjustments
    |  +-- CriticalValueEngine: 21-rule threshold check
    |  +-- DiscordanceDetector: cross-biomarker pattern analysis
    |  +-- LabRangeInterpreter: standard vs optimal range interpretation
    |
    v
[3. SEARCH]
    |  Multi-collection RAG across 14 collections
    |  ThreadPoolExecutor parallel search (max_workers=6)
    |  Weighted scoring with per-collection weights (sum = 1.0)
    |  Score threshold filtering (minimum 0.4 cosine similarity)
    |
    v
[4. SYNTHESIZE]
    |  Claude LLM synthesis with domain knowledge injection
    |  Grounded responses with evidence citations
    |  Conversation memory for follow-up questions
    |
    v
[5. REPORT]
    |  12-section clinical report generation
    |  Multi-format export: Markdown, PDF, FHIR R4, CSV, JSON
    |  Multi-language translation (7 languages)
    |
    v
AnalysisResult + Clinical Report
```

### Data Ingestion Flow

```
Reference Data (18 JSON files)
    |
    v
[seed_all.py]
    |  Read JSON -> Instantiate Pydantic models
    |  Call to_embedding_text() on each model
    |  Embed with BGE-small-en-v1.5 (batch_size=32)
    |  Insert into Milvus collections
    |
    v
14 Milvus Collections (indexed, searchable)
```

---

## Milvus Collections

The agent manages 14 specialized vector collections. All use COSINE similarity
with IVF_FLAT indexing (nlist=1024, nprobe=16) and 384-dimensional vectors
from BGE-small-en-v1.5.

| # | Collection Name | Description | Weight |
|---|----------------|-------------|--------|
| 1 | biomarker_reference | Reference biomarker definitions with standard/optimal ranges | 0.12 |
| 2 | biomarker_genetic_variants | Genetic variants affecting biomarker levels (MTHFR, APOE, etc.) | 0.11 |
| 3 | biomarker_pgx_rules | CPIC pharmacogenomic dosing rules | 0.10 |
| 4 | biomarker_disease_trajectories | Disease progression trajectories with intervention windows | 0.10 |
| 5 | biomarker_clinical_evidence | Published clinical evidence with PubMed references | 0.09 |
| 6 | genomic_evidence | Shared VCF-derived genomic variants (read-only) | 0.08 |
| 7 | biomarker_drug_interactions | Gene-drug interactions beyond PGx | 0.07 |
| 8 | biomarker_aging_markers | Epigenetic aging clock marker data (PhenoAge/GrimAge) | 0.07 |
| 9 | biomarker_nutrition | Genotype-aware nutrition guidelines | 0.05 |
| 10 | biomarker_genotype_adjustments | Genotype-based reference range adjustments | 0.05 |
| 11 | biomarker_monitoring | Condition-specific monitoring protocols | 0.05 |
| 12 | biomarker_critical_values | Critical threshold definitions | 0.04 |
| 13 | biomarker_discordance_rules | Cross-biomarker discordance patterns | 0.04 |
| 14 | biomarker_aj_carrier_screening | Ashkenazi Jewish carrier screening panel (10 genes) | 0.03 |
| | **Total** | | **1.00** |

### Collection Weight Distribution

```
biomarker_reference          ████████████          0.12
biomarker_genetic_variants   ███████████           0.11
biomarker_pgx_rules          ██████████            0.10
biomarker_disease_traject.   ██████████            0.10
biomarker_clinical_evidence  █████████             0.09
genomic_evidence             ████████              0.08
biomarker_drug_interactions  ███████               0.07
biomarker_aging_markers      ███████               0.07
biomarker_nutrition          █████                 0.05
biomarker_genotype_adjust.   █████                 0.05
biomarker_monitoring         █████                 0.05
biomarker_critical_values    ████                  0.04
biomarker_discordance_rules  ████                  0.04
biomarker_aj_carrier_screen  ███                   0.03
                                             Sum: 1.00
```

### Key Collection Schemas

**biomarker_reference** -- Core biomarker definitions:

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| name | VARCHAR(100) | Display name (e.g., "HbA1c") |
| unit | VARCHAR(20) | e.g., "%", "mg/dL" |
| category | VARCHAR(30) | CBC, CMP, Lipids, Thyroid |
| ref_range_min | FLOAT | Standard lower bound |
| ref_range_max | FLOAT | Standard upper bound |
| text_chunk | VARCHAR(3000) | Text for embedding |
| clinical_significance | VARCHAR(2000) | Clinical interpretation |
| epigenetic_clock | VARCHAR(50) | PhenoAge/GrimAge coefficient |
| genetic_modifiers | VARCHAR(500) | Comma-separated modifier genes |

**biomarker_pgx_rules** -- Pharmacogenomic dosing rules:

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| gene | VARCHAR(50) | Pharmacogene (CYP2D6, etc.) |
| star_alleles | VARCHAR(100) | Star allele combination |
| drug | VARCHAR(100) | Drug name |
| phenotype | VARCHAR(30) | MetabolizerPhenotype enum |
| cpic_level | VARCHAR(5) | CPIC evidence level (1A-3) |
| recommendation | VARCHAR(2000) | Dosing recommendation |
| evidence_url | VARCHAR(500) | CPIC/PharmGKB URL |
| text_chunk | VARCHAR(3000) | Text for embedding |

---

## Embedding Strategy

| Parameter | Value |
|-----------|-------|
| Model | BAAI/bge-small-en-v1.5 |
| Parameters | 33M |
| Dimensions | 384 |
| Metric | COSINE |
| Index type | IVF_FLAT (nlist=1024, nprobe=16) |
| Batch size | 32 |
| Runtime | CPU (no GPU required) |
| Instruction prefix | "Represent this sentence for searching relevant passages: " |
| Search mode | Asymmetric (queries use instruction prefix, documents do not) |
| Embedding cache | 256 entries, FIFO eviction |

### Domain-Optimized Embedding Text

Each collection model provides a `to_embedding_text()` method that constructs
a domain-optimized text representation for higher retrieval quality:

- **BiomarkerReference:** `"{name} ({unit}). {text_chunk}. Significance: {clinical_significance}. Category: {category}. Genetic modifiers: {genetic_modifiers}"`
- **PGxRule:** `"{gene} {star_alleles} -- {drug}. {text_chunk}. Recommendation: {recommendation}. Phenotype: {phenotype}. CPIC Level: {cpic_level}"`
- **GeneticVariant:** `"{gene} {rs_id}. {text_chunk}. Mechanism: {mechanism}. Diseases: {disease_associations}"`
- **DiseaseTrajectory:** `"{disease} {stage}. {text_chunk}. Intervention: {intervention_window}"`

This approach ensures that semantically critical fields (gene names, drug names,
phenotypes, clinical significance) are included in the embedding text even
when they appear only in structured metadata fields.

---

## LLM Integration

### Claude Configuration

| Setting | Value |
|---------|-------|
| Model | claude-sonnet-4-6 |
| Provider | Anthropic API |
| Environment variable | `ANTHROPIC_API_KEY` |
| Max retries | 3 |
| Conversation memory | 3 turns (configurable) |
| Streaming | Supported (via `query_stream()`) |

### System Prompt Design

The system prompt defines the agent as a **Precision Biomarker Intelligence Agent**
with expertise in:

1. Genotype-aware biomarker interpretation
2. Pharmacogenomic phenotyping (CPIC guidelines)
3. Biological age estimation (PhenoAge/GrimAge)
4. Disease trajectory prediction
5. Critical value alerting
6. Cross-biomarker discordance detection
7. Genotype-based reference range adjustments
8. Longitudinal trend analysis

Behavioral instructions enforce:
- Evidence-based responses with citations
- Acknowledgment of uncertainty and limitations
- Cross-referencing across biomarker domains
- Pharmacogenomic context in drug-related queries
- Flagging of critical values requiring immediate attention

### Knowledge Injection

Before each LLM call, the RAG engine injects structured domain knowledge from
the knowledge graph alongside retrieved evidence. This provides the LLM with
curated clinical context that may not be captured in the vector store alone,
including CPIC dosing tables, disease trajectory staging criteria, and genotype
adjustment coefficients.

---

## Knowledge Graph

**Module:** `src/knowledge.py` (1,326 lines)

The knowledge graph provides curated, structured domain knowledge organized
into 6 disease domains. This knowledge is injected into LLM prompts alongside
retrieved vector evidence to ensure comprehensive clinical context.

### Domain Structure

| Domain | Content | Purpose |
|--------|---------|---------|
| Pharmacogenomics | 14 pharmacogene profiles with star allele tables | PGx phenotyping context |
| Cardiovascular | Risk factors, biomarker thresholds, treatment targets | CV disease trajectory |
| Metabolic | Diabetes, thyroid, nutritional biomarker patterns | Metabolic disease detection |
| Hepatic/Renal | Liver function, kidney function, staging criteria | Organ function assessment |
| Hematology | CBC interpretation, iron studies, coagulation | Blood disorder detection |
| Aging | PhenoAge/GrimAge coefficients, aging biomarker norms | Biological age estimation |

### Helper Functions

- `lookup_biomarker(name)` -- Return full biomarker context including reference ranges, genetic modifiers, and clinical significance
- `lookup_pharmacogene(gene)` -- Return star allele tables, phenotype classifications, and affected drugs
- `lookup_disease_domain(category)` -- Return disease-specific biomarker patterns and trajectory staging
- `lookup_interaction(drug, gene)` -- Return gene-drug interaction context

---

## Clinical Pipelines

### 1. Pharmacogenomic Engine (src/pharmacogenomics.py, 1,503 lines)

Maps star alleles and genotypes to drug recommendations following CPIC Level 1A
guidelines. Pure computation with no LLM or database calls required.

**Covered pharmacogenes (14):**

| Gene | Description | CPIC Level | Key Drugs |
|------|-------------|------------|-----------|
| CYP2D6 | Metabolizes ~25% of drugs | 1A | Codeine, tramadol, tamoxifen |
| CYP2C19 | Clopidogrel, PPIs, antidepressants | 1A | Clopidogrel, omeprazole |
| CYP2C9 | Warfarin, NSAIDs, phenytoin | 1A | Warfarin, celecoxib |
| CYP3A5 | Tacrolimus metabolism | 1A | Tacrolimus |
| SLCO1B1 | Statin hepatic uptake transporter | 1A | Simvastatin, atorvastatin |
| VKORC1 | Warfarin target sensitivity | 1A | Warfarin |
| MTHFR | Folate metabolism, homocysteine | Informational | Methotrexate |
| TPMT | Thiopurine metabolism | 1A | Azathioprine, mercaptopurine |
| DPYD | Fluoropyrimidine metabolism | 1A | 5-FU, capecitabine |

**Phenotype classification:** Uses CPIC standard metabolizer terminology
(Normal/Intermediate/Poor/Ultra-rapid/Rapid Metabolizer) with gene-specific
variants for transporters (SLCO1B1), enzymes (MTHFR), and sensitivity (VKORC1).

**Drug recommendation actions:** STANDARD_DOSING, DOSE_REDUCTION, DOSE_ADJUSTMENT,
CONSIDER_ALTERNATIVE, AVOID, CONTRAINDICATED. Alert levels: INFO, WARNING, CRITICAL.

### 2. Biological Age Engine (src/biological_age.py, 408 lines)

Estimates biological age using two validated epigenetic aging clocks:

- **PhenoAge (Levine 2018):** Uses 9 blood biomarkers (albumin, creatinine,
  glucose, CRP, lymphocyte %, mean cell volume, red cell distribution width,
  alkaline phosphatase, white blood cell count) plus chronological age.
- **GrimAge surrogate:** Uses a subset of PhenoAge biomarkers with smoking
  pack-years and additional aging markers.

Output includes biological age estimate, acceleration (difference from
chronological age), percentile ranking, and contributing biomarkers sorted
by impact.

### 3. Disease Trajectory Engine (src/disease_trajectory.py, 1,421 lines)

Predicts disease progression across 9 categories:

| Category | Key Biomarkers | Stages |
|----------|---------------|--------|
| Diabetes | HbA1c, fasting glucose, HOMA-IR | Pre-diabetic -> T2DM -> Complications |
| Cardiovascular | LDL, CRP, Lp(a), BNP | Risk factors -> Subclinical -> Events |
| Liver | ALT, AST, GGT, albumin, bilirubin | Steatosis -> NASH -> Fibrosis -> Cirrhosis |
| Thyroid | TSH, free T4, free T3, TPO antibodies | Subclinical -> Overt -> Complications |
| Iron | Ferritin, transferrin saturation, TIBC | Depletion -> Deficiency -> Anemia |
| Nutritional | B12, folate, vitamin D, zinc | Suboptimal -> Deficient -> Clinical |
| Kidney | eGFR, creatinine, BUN, cystatin C | Stage 1 -> Stage 5 (CKD) |
| Bone Health | Calcium, phosphorus, PTH, vitamin D | Normal -> Osteopenia -> Osteoporosis |
| Cognitive | Homocysteine, B12, folate, CRP | Risk factors -> Mild -> Significant |

Each trajectory includes years-to-diagnosis estimates, intervention windows,
and potential risk reduction percentages.

### 4. Genotype Adjustment Engine (src/genotype_adjustment.py, 1,225 lines)

Adjusts standard reference ranges based on 7 modifier genes:

| Gene | Variant | Affected Biomarkers | Adjustment |
|------|---------|-------------------|------------|
| MTHFR | C677T | Homocysteine, folate | Narrowed ranges, higher folate targets |
| APOE | e4 allele | LDL, total cholesterol | Tighter lipid targets |
| PNPLA3 | I148M | ALT, AST, GGT | Adjusted liver enzyme thresholds |
| HFE | C282Y, H63D | Ferritin, transferrin sat | Iron overload-specific ranges |
| GBA | Various | GCase activity | Carrier-specific enzyme ranges |
| HEXA | Various | Hex A activity | Carrier-specific enzyme ranges |
| FTO | rs9939609 | BMI, glucose, HbA1c | Obesity-adjusted metabolic ranges |

### 5. Critical Value Engine (src/critical_values.py, 179 lines)

Monitors 21 critical threshold rules for life-threatening lab values that
require immediate clinical attention. Each rule defines:

- Biomarker name and threshold (high/low)
- Severity level (critical, urgent)
- Escalation target (emergency, physician, specialist)
- Required action text

### 6. Discordance Detector (src/discordance_detector.py, 299 lines)

Identifies cross-biomarker discordances -- situations where two or more
biomarker results are physiologically contradictory or clinically unexpected.
Examples include:

- Normal TSH with elevated free T4 (suggesting assay interference)
- Normal HbA1c with elevated fasting glucose (possible hemoglobin variant)
- Low ferritin with elevated transferrin saturation (mixed iron picture)

---

## RAG Engine

**Module:** `src/rag_engine.py` (573 lines)

### Search Flow

```
Query --> Embed (BGE-small + instruction prefix)
  --> Parallel Search (14 collections, top_k=5 each, score_threshold=0.40)
  --> Weighted Scoring (per-collection weights, sum=1.0)
  --> Deduplication (by ID and text content hash)
  --> Knowledge Augmentation (pharmacogenes, disease domains, interactions)
  --> LLM Synthesis (Claude, system prompt, conversation memory)
  --> Streaming Response with Evidence Citations
```

### Citation Scoring

| Score Range | Relevance Level |
|------------|----------------|
| >= 0.75 | High |
| >= 0.60 | Medium |
| < 0.60 | Low |

### Conversation Memory

The RAG engine maintains thread-safe conversation memory using a deque with
configurable maximum size (default: 3 turns). This allows follow-up questions
to reference previous context without re-executing the full analysis pipeline.

---

## API Layer

### FastAPI Server (port 8529)

The API is organized across 3 route modules with a shared lifespan manager.

**Core endpoints (api/main.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/healthz` | GET | Health check (returns service status) |
| `/readyz` | GET | Readiness check (Milvus connection verified) |
| `/metrics` | GET | Prometheus metrics endpoint |
| `/v1/query` | POST | RAG query with streaming support |

**Analysis endpoints (api/routes/analysis.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/analyze` | POST | Full patient analysis (all 6 engines) |
| `/v1/biological-age` | POST | PhenoAge/GrimAge calculation |
| `/v1/disease-trajectory` | POST | Disease trajectory prediction |
| `/v1/pgx-map` | POST | Pharmacogenomic mapping |
| `/v1/genotype-adjust` | POST | Genotype-based range adjustment |

**Event endpoints (api/routes/events.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/events/cross-modal` | POST | Cross-modal event handling |
| `/v1/biomarker-alert` | POST | Critical value alert processing |

**Report endpoints (api/routes/reports.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/report/generate` | POST | Generate clinical report |
| `/v1/report/{id}/pdf` | GET | Download report as PDF |
| `/v1/report/fhir` | POST | Generate FHIR R4 bundle |

### Authentication and Security

- API key authentication via `BIOMARKER_API_KEY` environment variable
  (empty string disables authentication for development)
- CORS middleware with configurable allowed origins
- Request size limiting (default 10 MB)
- Input sanitization for Milvus filter expressions

### Health Checks

- Docker health check: `/_stcore/health` (Streamlit), `/healthz` (API)
- Check interval: 30 seconds
- Start period: 30-60 seconds

---

## UI Layer

### Streamlit Application (port 8528)

**Module:** `app/biomarker_ui.py` (1,863 lines)

The UI provides an 8-tab interface for interactive biomarker analysis:

| Tab | Name | Functionality |
|-----|------|--------------|
| 1 | Clinical Query | RAG search with streaming responses and citation scoring |
| 2 | Patient Analysis | Full pipeline execution across all 6 engines |
| 3 | Pharmacogenomics | PGx panel with star allele input and drug recommendations |
| 4 | Biological Age | PhenoAge/GrimAge calculator with biomarker input |
| 5 | Disease Trajectory | 9-category disease risk assessment with intervention windows |
| 6 | Genotype Adjustments | Reference range adjustments for 7 modifier genes |
| 7 | Critical Values | Real-time critical value monitoring and alerting |
| 8 | Patient 360 | Cross-agent dashboard (Biomarker + Oncology + Autoimmune + CAR-T + Imaging) |

### Patient 360 Dashboard (app/patient_360.py, 670 lines)

Cross-agent view that aggregates findings from all five intelligence agents
into a unified patient profile. Displays:

- Biomarker summary with critical values highlighted
- PGx phenotype overview
- Oncology variant status (from Oncology Agent)
- Autoimmune disease activity (from Autoimmune Agent)
- Imaging findings (from Imaging Agent)
- CAR-T eligibility assessment (from CAR-T Agent)

---

## Export Pipeline

**Module:** `src/export.py` (1,392 lines)

Five export formats, each accepting the full AnalysisResult:

### Markdown Export

Structured clinical report with 12 sections including critical alerts,
PGx summary table, biological age assessment, disease trajectory findings,
genotype adjustments, discordance alerts, and evidence citations.

### PDF Export

Styled report using ReportLab with NVIDIA-themed colors (green #76B900 headers),
tabular PGx results, disease trajectory charts, and clinical footer with
knowledge base version.

### FHIR R4 Export

Generates a FHIR R4 Bundle (type: collection) containing:

| Resource | Content |
|----------|---------|
| Patient | Identifier with urn:hcls-ai-factory:patient |
| Observation (N) | One per biomarker result |
| Observation (PGx) | Pharmacogenomic phenotype observations |
| DiagnosticReport | Master biomarker intelligence report |

### CSV Export

Flat tabular export of biomarker results with reference ranges, genotype
adjustments, and clinical significance flags. Suitable for integration with
external analytics tools and EHR systems.

### JSON Export

Structured JSON export with full analysis results, metadata, and provenance
information. Includes the complete AnalysisResult model serialized via
Pydantic's `model_dump()`.

---

## Scaling and Performance

### Resource Footprint on DGX Spark

| Component | Resource Usage |
|-----------|---------------|
| Milvus standalone | ~2-4 GB RAM (14 collections) |
| BGE-small-en-v1.5 | ~200 MB RAM (CPU inference) |
| FastAPI server | ~500 MB RAM (2 uvicorn workers) |
| Streamlit UI | ~300 MB RAM |
| Total agent footprint | ~3-5 GB RAM |

### Performance Characteristics

| Operation | Typical Latency |
|-----------|----------------|
| Single collection search | 10-50 ms |
| 14-collection parallel search | 50-200 ms |
| BGE-small embedding (single) | 5-15 ms |
| Claude LLM synthesis | 2-8 seconds |
| Full patient analysis (6 engines) | 3-12 seconds |
| Report generation (Markdown) | < 100 ms |
| PDF generation | 200-500 ms |
| FHIR R4 export | < 100 ms |

### Parallelization

- Collection searches run in parallel via ThreadPoolExecutor (max_workers=6)
- Analysis engines run in parallel during the ANALYZE phase
- Embedding batching (batch_size=32) for bulk operations
- Connection pooling for Milvus gRPC connections

### Caching

- Embedding cache: 256 entries, FIFO eviction, thread-safe
- Knowledge graph: In-memory dictionaries, loaded at startup
- Reference data: Loaded from JSON at seed time, stored in Milvus

---

## Security

### Authentication

- API key authentication via `BIOMARKER_API_KEY` environment variable
- Empty API key disables authentication (development mode only)
- Bearer token format: `Authorization: Bearer <api_key>`

### Data Protection

- HIPAA-compliant audit logging (`src/audit.py`, 83 lines)
- All patient data encrypted at rest via Milvus storage encryption
- No PHI stored in application logs
- Anthropic API key stored in environment variables, never committed to source

### Input Validation

- Pydantic v2 model validation on all API inputs
- Milvus filter expression sanitization to prevent injection
- Request size limiting (configurable, default 10 MB)
- CORS middleware with explicit origin allowlist

### Network Security

- All services communicate via Docker bridge network (`biomarker-network`)
- Only Streamlit (8528) and FastAPI (8529) ports are exposed externally
- Milvus gRPC (19530) is internal to the Docker network
- Non-root container execution (`biomarkeruser`)

### Audit Trail

The audit module (`src/audit.py`) provides structured logging of:
- Patient data access events
- Analysis requests and results
- Report generation and export events
- API authentication events
- Configuration changes

---

## Cross-Agent Integration

### Platform Context: Precision Intelligence Network

The Precision Biomarker Intelligence Agent operates within the HCLS AI Factory
Precision Intelligence Network, one of three GPU-accelerated engines:

1. **Genomic Foundation Engine** -- Parabricks/DeepVariant/BWA-MEM2 for
   FASTQ-to-VCF processing (Stage 1)
2. **Precision Intelligence Network** -- 11 domain-specialized RAG agents
   providing clinical decision support (Stage 2)
3. **Therapeutic Discovery Engine** -- BioNeMo MolMIM/DiffDock/RDKit for
   molecular generation and docking (Stage 3)

### The 11 HCLS AI Factory Intelligence Agents

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

### Cross-Agent Calls via cross_modal/cross_agent.py

The Biomarker Agent calls 4 peer agents via `cross_modal/cross_agent.py` and
the `/integrated-assessment` endpoint for multi-agent clinical synthesis:

| Target Agent | Integration Purpose |
|-------------|-------------------|
| Precision Oncology Agent | DPYD poor metabolizer alerts, BRCA1/2 carrier cancer risk |
| CAR-T Intelligence Agent | CRS/ICANS biomarker correlation (ferritin, CRP, IL-6) |
| Pharmacogenomics (PGx) Agent | Extended PGx panel integration beyond the 14 core pharmacogenes |
| Clinical Trial Intelligence Agent | Biomarker-driven trial eligibility matching |

### Integrated Assessment Endpoint

The `/integrated-assessment` endpoint orchestrates multi-agent workflows by:

1. Collecting the biomarker agent's full patient analysis (PGx, biological age, trajectories)
2. Requesting variant-level therapy matching from the Oncology Agent
3. Querying the CAR-T Agent for CRS biomarker correlation when relevant
4. Checking extended PGx interactions via the PGx Agent
5. Aggregating findings into a unified Patient 360 assessment

---

## File Structure

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
|
|-- docs/
|   |-- ARCHITECTURE_GUIDE.md    # This document
|   |-- DEMO_GUIDE.md
|   |-- DEPLOYMENT_GUIDE.md
|   |-- DESIGN.md
|   |-- INDEX.md
|   |-- LEARNING_GUIDE_ADVANCED.md
|   |-- LEARNING_GUIDE_FOUNDATIONS.md
|   |-- PROJECT_BIBLE.md
|   +-- WHITE_PAPER.md
|
|-- src/                         # Core engine modules
|   |-- __init__.py
|   |-- agent.py                 (610 lines) Plan-analyze-search-synthesize loop
|   |-- audit.py                 (83 lines) HIPAA-compliant audit logging
|   |-- biological_age.py        (408 lines) PhenoAge + GrimAge surrogate
|   |-- collections.py           (1,391 lines) Milvus collection management
|   |-- critical_values.py       (179 lines) Critical threshold detection
|   |-- discordance_detector.py  (299 lines) Cross-biomarker anomaly detection
|   |-- disease_trajectory.py    (1,421 lines) 9-category disease trajectories
|   |-- export.py                (1,392 lines) PDF/FHIR/CSV/JSON/Markdown export
|   |-- genotype_adjustment.py   (1,225 lines) Genotype-based range adjustments
|   |-- knowledge.py             (1,326 lines) 6 disease domain knowledge graphs
|   |-- lab_range_interpreter.py (221 lines) Standard vs optimal ranges
|   |-- models.py                (786 lines) 14 collection + 8 analysis models
|   |-- pharmacogenomics.py      (1,503 lines) CPIC-guided 14 pharmacogenes
|   |-- rag_engine.py            (573 lines) Multi-collection RAG engine
|   |-- report_generator.py      (993 lines) 12-section clinical reports
|   +-- translation.py           (217 lines) Multi-language medical terminology
|
|-- scripts/                     # Setup, seeding, and validation
|   |-- setup_collections.py     Create Milvus collections
|   |-- seed_all.py              Seed all 14 collections from JSON
|   |-- gen_patient_data.py      Generate sample patient data
|   +-- demo_validation.py       End-to-end demo validation
|
|-- tests/                       # 18 test files, 709 tests
|   |-- conftest.py              Shared fixtures
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
|-- pyproject.toml               Project metadata
|-- requirements.txt             Python dependencies
+-- README.md                    Quick start guide
```
