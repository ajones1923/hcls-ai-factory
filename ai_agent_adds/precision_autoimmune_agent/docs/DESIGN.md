# Precision Autoimmune Intelligence Agent -- Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**Version:** 2.0.0
**License:** Apache 2.0

---

## 1. Executive Summary

The Precision Autoimmune Intelligence Agent extends the HCLS AI Factory platform to deliver comprehensive autoimmune disease intelligence. Unlike standard rheumatology workflows that require manual synthesis of autoantibody panels, HLA typing, disease activity scores, and therapy options across fragmented records, this agent integrates all data sources into a unified analysis pipeline. It ingests clinical documents (PDFs), interprets autoantibody patterns, evaluates HLA-disease associations, scores disease activity across 20 validated instruments, predicts flares from biomarker patterns, and recommends biologic therapies with pharmacogenomic context.

The agent combines **6 deterministic clinical analysis engines** with a **14-collection RAG pipeline** to answer questions like *"What is my flare risk given rising anti-dsDNA and falling C3/C4?"* -- simultaneously searching autoantibody databases, disease activity scoring references, flare prediction patterns, and biologic therapy databases, then synthesizing a grounded response through Claude.

### Key Results

| Metric | Value |
|---|---|
| Milvus collections | **14** (13 owned + 1 read-only) |
| Clinical analysis engines | **6** deterministic engines |
| Diseases covered | **13** (RA, SLE, MS, T1D, IBD, Psoriasis, AS, Sjogren's, SSc, MG, Celiac, Graves, Hashimoto) |
| Biologic therapies mapped | **22** across 7 drug classes |
| Autoantibodies characterized | **24** with disease associations |
| HLA alleles mapped | **22** with odds ratios |
| Activity scoring systems | **20** across 13 diseases |
| Demo patients | **9** with clinical document PDFs |
| Unit tests passing | **455** |
| Export formats | **3** (FHIR R4, PDF, Markdown) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Autoimmune Agent Role |
|---|---|
| **DataStore** | Clinical document PDFs, demo patient data, knowledge base (HLA, autoantibodies, biologics, flare patterns) |
| **DataEngine** | Document processor: PDF -> chunk -> BGE-small embedding -> Milvus vector insert; Knowledge seeder: JSON -> embeddings -> collections |
| **DataBase** | 14 Milvus collections (13 autoimmune-domain + 1 shared read-only) with 9 demo patients |
| **InsightEngine** | 6 clinical analysis engines + BGE-small embedding + multi-collection RAG with weighted parallel search |
| **AgentEngine** | AutoimmuneAgent orchestrator + Streamlit UI (10 tabs) + FastAPI REST (14 endpoints) |

### 2.2 System Diagram

```
                        +-----------------------------------------+
                        |    Streamlit UI (8531)                    |
                        |    10 tabs: Query | Analysis | Ingest |   |
                        |    Odyssey | Autoantibody | HLA |         |
                        |    Activity | Flare | Therapy | KB        |
                        +-------------------+---------------------+
                                            |
                        +-------------------v---------------------+
                        |  AutoimmuneAgent                         |
                        |  Orchestrates 6 analysis engines         |
                        |  + Document processor + RAG pipeline     |
                        +-------------------+---------------------+
                                            |
            +-------------------------------+-------------------------------+
            |                               |                               |
  +---------v-----------+       +-----------v-----------+       +-----------v-----------+
  | Deterministic        |       | RAG Pipeline           |       | Export                 |
  | Analysis Engines     |       |                        |       |                        |
  |                      |       | BGE-small-en-v1.5      |       | FHIR R4 Bundle         |
  | AutoantibodyInterp   |       | (384-dim embedding)    |       | PDF (reportlab)        |
  | HLAAssociationAnlyz  |       |         |              |       | Markdown               |
  | DiseaseActivityScor  |       |         v              |       |                        |
  | FlarePredictor       |       | Parallel Search        |       | + FHIR Validation      |
  | BiologicTherapyAdv   |       | 14 Milvus Collections  |       |                        |
  | DiagnosticOdyssey    |       | (ThreadPoolExecutor)   |       |                        |
  |                      |       |         |              |       |                        |
  +----------------------+       |         v              |       |                        |
            |                    | Claude Sonnet 4.6      |       |                        |
            |                    +------------------------+       +------------------------+
            |                               |
  +---------v-------------------------------v-------------------------------+
  |                  Milvus 2.4 -- 14 Collections                           |
  |                                                                         |
  |  autoimmune_clinical_documents       autoimmune_patient_labs            |
  |  autoimmune_autoantibody_panels      autoimmune_hla_associations        |
  |  autoimmune_disease_criteria         autoimmune_disease_activity        |
  |  autoimmune_flare_patterns           autoimmune_biologic_therapies      |
  |  autoimmune_pgx_rules               autoimmune_clinical_trials          |
  |  autoimmune_literature               autoimmune_patient_timelines       |
  |  autoimmune_cross_disease            genomic_evidence [RO]              |
  +-------------------------------------------------------------------------+
```

---

## 3. Data Collections -- Actual State

### 3.1 `autoimmune_clinical_documents` -- variable records

Ingested clinical documents (PDFs) chunked and embedded for RAG retrieval.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Chunked document text |
| patient_id | VARCHAR(64) | Patient identifier |
| doc_type | VARCHAR(128) | Document type (progress_note, lab_report, imaging, pathology, referral) |
| specialty | VARCHAR(128) | Medical specialty (rheumatology, neurology, nephrology, etc.) |
| provider | VARCHAR(256) | Provider name |
| visit_date | VARCHAR(32) | Date of visit |
| source_file | VARCHAR(512) | Source PDF filename |
| page_number | INT64 | Page number within document |
| chunk_index | INT64 | Chunk index within page |

### 3.2 `autoimmune_patient_labs` -- variable records

Laboratory results with reference range analysis and flag classification.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Lab result text |
| patient_id | VARCHAR(64) | Patient identifier |
| test_name | VARCHAR(256) | Test name |
| value | FLOAT | Measured value |
| unit | VARCHAR(64) | Unit of measurement |
| reference_range | VARCHAR(128) | Normal reference range |
| flag | VARCHAR(32) | Status: normal, high, low, critical |
| collection_date | VARCHAR(32) | Collection date |
| panel_name | VARCHAR(256) | Lab panel name |

### 3.3 `autoimmune_autoantibody_panels` -- 24 records

Autoantibody reference panels with disease associations, sensitivity, and specificity.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Autoantibody description |
| antibody_name | VARCHAR(128) | Antibody name (ANA, anti-dsDNA, RF, anti-CCP, etc.) |
| associated_diseases | VARCHAR(1024) | Diseases associated with this antibody |
| sensitivity | FLOAT | Diagnostic sensitivity |
| specificity | FLOAT | Diagnostic specificity |
| pattern | VARCHAR(128) | Staining pattern (homogeneous, speckled, nucleolar, centromere) |
| clinical_significance | VARCHAR(2000) | Clinical interpretation |
| interpretation_guide | VARCHAR(2000) | Interpretation guidelines |

### 3.4 `autoimmune_hla_associations` -- 22 records

HLA allele to autoimmune disease associations with odds ratios and population context.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Association description |
| allele | VARCHAR(64) | HLA allele (e.g., HLA-B*27:05, HLA-DRB1*04:01) |
| disease | VARCHAR(256) | Associated disease |
| odds_ratio | FLOAT | Odds ratio for disease risk |
| population | VARCHAR(128) | Population studied |
| pmid | VARCHAR(32) | PubMed identifier |
| mechanism | VARCHAR(1024) | Proposed mechanism |
| clinical_implication | VARCHAR(2000) | Clinical significance |

### 3.5 `autoimmune_disease_criteria` -- 13 records

ACR/EULAR classification and diagnostic criteria for 13 autoimmune diseases.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Criteria description |
| disease | VARCHAR(256) | Disease name |
| criteria_set | VARCHAR(256) | Criteria set (e.g., 2010 ACR/EULAR RA, 2019 ACR/EULAR SLE) |
| criteria_type | VARCHAR(64) | Classification or diagnostic |
| year | INT64 | Publication year |
| required_score | VARCHAR(128) | Minimum score for classification |
| criteria_items | VARCHAR(3000) | Individual criteria items with point values |
| sensitivity_specificity | VARCHAR(256) | Reported sensitivity and specificity |

### 3.6 `autoimmune_disease_activity` -- 20 records

Disease activity scoring systems and interpretation guides.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Scoring system description |
| score_name | VARCHAR(128) | Score name (DAS28-CRP, SLEDAI-2K, BASDAI, CDAI, EDSS, etc.) |
| disease | VARCHAR(256) | Associated disease |
| components | VARCHAR(2000) | Score components |
| thresholds | VARCHAR(1024) | Remission/low/moderate/high cutoffs (JSON) |
| interpretation | VARCHAR(2000) | Interpretation guidelines |
| monitoring_frequency | VARCHAR(512) | Recommended monitoring frequency |

### 3.7 `autoimmune_flare_patterns` -- 13 records

Flare prediction biomarker patterns and early warning signs.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Flare pattern description |
| disease | VARCHAR(256) | Disease name |
| biomarker_pattern | VARCHAR(2000) | Biomarker changes preceding flare |
| early_warning_signs | VARCHAR(2000) | Early warning signs |
| typical_timeline | VARCHAR(512) | Typical flare timeline |
| protective_factors | VARCHAR(1024) | Factors that reduce flare risk |
| intervention_triggers | VARCHAR(1024) | When to intervene |

### 3.8 `autoimmune_biologic_therapies` -- 22 records

Biologic therapy database with pharmacogenomic considerations.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Therapy description |
| drug_name | VARCHAR(128) | Drug name |
| drug_class | VARCHAR(128) | Drug class (TNF inhibitor, IL-6 inhibitor, B-cell depleter, etc.) |
| mechanism | VARCHAR(512) | Mechanism of action |
| indicated_diseases | VARCHAR(1024) | FDA-approved indications |
| pgx_considerations | VARCHAR(2000) | Pharmacogenomic factors |
| contraindications | VARCHAR(1024) | Contraindications |
| monitoring | VARCHAR(2000) | Required monitoring |
| dosing | VARCHAR(512) | Dosing regimen |
| evidence_level | VARCHAR(64) | Evidence level |

### 3.9 `autoimmune_pgx_rules` -- variable records

Pharmacogenomic dosing rules for autoimmune therapies.

### 3.10 `autoimmune_clinical_trials` -- variable records

Active and completed autoimmune disease clinical trials.

### 3.11 `autoimmune_literature` -- variable records

Published autoimmune research literature with abstracts and keywords.

### 3.12 `autoimmune_patient_timelines` -- variable records

Patient diagnostic timeline events for odyssey analysis.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Event description |
| patient_id | VARCHAR(64) | Patient identifier |
| event_type | VARCHAR(128) | Event type (symptom_onset, diagnosis, treatment_start, flare, etc.) |
| event_date | VARCHAR(32) | Date of event |
| description | VARCHAR(2000) | Event description |
| provider | VARCHAR(256) | Provider |
| specialty | VARCHAR(128) | Specialty |
| days_from_first_symptom | INT64 | Days from first symptom |

### 3.13 `autoimmune_cross_disease` -- 9 records

Cross-disease overlap syndromes and shared pathogenic mechanisms.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Overlap description |
| primary_disease | VARCHAR(256) | Primary disease |
| associated_conditions | VARCHAR(1024) | Associated conditions |
| shared_pathways | VARCHAR(1024) | Shared immunological pathways |
| shared_biomarkers | VARCHAR(1024) | Shared biomarkers |
| overlap_criteria | VARCHAR(2000) | Overlap classification criteria |
| co_occurrence_rate | FLOAT | Co-occurrence rate |

### 3.14 Index Configuration (all collections)

```
Algorithm:  IVF_FLAT
Metric:     COSINE
nlist:      1024
nprobe:     16
Dimension:  384 (BGE-small-en-v1.5)
```

---

## 4. Clinical Analysis Engines

### 4.1 AutoantibodyInterpreter

Interprets autoantibody panels by mapping antibody combinations to disease associations:

- 24 autoantibodies with sensitivity and specificity data
- Pattern recognition for ANA staining (homogeneous, speckled, nucleolar, centromere)
- Titer interpretation (1:40 through 1:1280+)
- Multi-antibody combination scoring (e.g., ANA + anti-dsDNA + anti-Smith = SLE)
- Positive antibody count and disease association ranking

### 4.2 HLAAssociationAnalyzer

Maps HLA alleles to disease susceptibility across 22 characterized alleles:

| Notable Associations | Allele | Disease | Odds Ratio |
|---|---|---|---|
| Strongest known | HLA-B*27:05 | Ankylosing Spondylitis | 87.4 |
| Shared epitope | HLA-DRB1*04:01 | Rheumatoid Arthritis | 4.2 |
| Celiac/T1D overlap | HLA-DQB1*02:01 | Celiac Disease | 7.0 |
| MS susceptibility | HLA-DRB1*15:01 | Multiple Sclerosis | 3.1 |
| Psoriasis PSORS1 | HLA-C*06:02 | Psoriasis | 10.0 |
| T1D highest risk | HLA-DQB1*03:02 | Type 1 Diabetes | 6.5 |

### 4.3 DiseaseActivityScorer

20 validated disease activity scoring systems across 13 diseases:

| Disease | Score | Thresholds |
|---|---|---|
| Rheumatoid Arthritis | DAS28-CRP | Remission <2.6, Low <3.2, Moderate <5.1, High >=5.1 |
| Rheumatoid Arthritis | CDAI | Remission <=2.8, Low <=10, Moderate <=22, High >22 |
| SLE | SLEDAI-2K | No activity 0, Mild 1-5, Moderate 6-10, High 11-19, Very high >=20 |
| Ankylosing Spondylitis | BASDAI | Inactive <4, Active >=4 |
| Multiple Sclerosis | EDSS | Minimal 0-2.5, Moderate 3-5.5, Severe >=6 |
| IBD (Crohn's) | CDAI (Harvey-Bradshaw) | Remission <150, Mild 150-220, Moderate 220-450, Severe >450 |

Each score returns: disease, score_name, score_value, level (remission/low/moderate/high/very_high), components, and thresholds.

### 4.4 FlarePredictor

Biomarker pattern-based flare prediction for 13 diseases:

- 4-tier risk stratification: LOW (<0.4), MODERATE (0.4-0.6), HIGH (0.6-0.8), IMMINENT (>=0.8)
- Contributing factors identification
- Protective factors identification
- Recommended monitoring schedules
- 90-day prediction window (configurable)

### 4.5 BiologicTherapyAdvisor

22 biologic therapies across 7 drug classes:

| Drug Class | Examples | Key Indications |
|---|---|---|
| TNF inhibitors | Adalimumab, Etanercept, Infliximab | RA, AS, Psoriasis, IBD |
| IL-6 inhibitors | Tocilizumab, Sarilumab | RA |
| B-cell depleters | Rituximab, Ocrelizumab | RA, SLE, MS |
| JAK inhibitors | Tofacitinib, Baricitinib, Upadacitinib | RA, Psoriasis, IBD |
| IL-17 inhibitors | Secukinumab, Ixekizumab | AS, Psoriasis |
| IL-23 inhibitors | Guselkumab, Risankizumab | Psoriasis, IBD |
| BLyS inhibitors | Belimumab | SLE |

Each recommendation includes PGx considerations, contraindications, monitoring requirements, and evidence level.

### 4.6 DiagnosticOdysseyAnalyzer

Analyzes diagnostic timelines from ingested clinical documents:

- Timeline construction from clinical document events
- Event classification: symptom_onset, diagnosis, treatment_start, flare, referral, lab_result
- Delay identification between symptom onset and diagnosis
- Pattern recognition for common diagnostic pitfalls
- Multi-specialty referral tracking

---

## 5. Multi-Collection RAG Engine

### 5.1 Search Flow

```
Query Text
    |
    v
BGE-small-en-v1.5 Embedding (384-dim)
    |
    v
ThreadPoolExecutor: Parallel search across 14 collections (max 6 workers)
    |
    v
Weighted merge (configurable per-collection weights)
    |
    v
Score threshold filtering (>= 0.40)
    |
    v
Knowledge context augmentation
    |
    v
Claude Sonnet 4.6 prompt with patient context
    |
    v
Grounded response with citations (streaming supported)
```

### 5.2 Collection Weights

| Collection | Weight | Rationale |
|---|---|---|
| clinical_documents | 0.18 | Primary patient clinical narrative |
| patient_labs | 0.14 | Lab results with flag analysis |
| autoantibody_panels | 0.12 | Autoantibody disease associations |
| hla_associations | 0.08 | HLA genetic susceptibility |
| disease_criteria | 0.08 | Classification criteria |
| disease_activity | 0.07 | Activity scoring reference |
| flare_patterns | 0.06 | Flare prediction patterns |
| biologic_therapies | 0.06 | Therapy database |
| clinical_trials | 0.05 | Clinical trial matching |
| literature | 0.05 | Published evidence |
| pgx_rules | 0.04 | Pharmacogenomic rules |
| patient_timelines | 0.03 | Diagnostic timelines |
| cross_disease | 0.02 | Overlap syndromes |
| genomic_evidence | 0.02 | Shared genomic context |
| **Total** | **1.00** | |

### 5.3 Embedding Strategy

**Model:** BGE-small-en-v1.5 (BAAI)

- Dimension: 384
- Metric: Cosine similarity
- Query prefix: `"Represent this sentence for searching relevant passages: "`
- Document embedding: Raw text (no prefix)

**`to_embedding_text()` pattern:** Each Pydantic model implements this method to produce an optimal embedding string combining key fields.

### 5.4 Citation Scoring

| Level | Threshold | Display |
|---|---|---|
| High confidence | >= 0.80 | Full citation with source link |
| Medium confidence | >= 0.60 | Citation with caveat |
| Below threshold | < 0.40 | Filtered out |

---

## 6. Export Pipeline

### 6.1 FHIR R4 DiagnosticReport

Produces a FHIR R4 Bundle containing:

- **Patient** resource with identifier, gender, birth year
- **DiagnosticReport** resource (main report)
- **Observation** resources for autoantibody results and disease activity scores
- **Condition** resources for diagnosed autoimmune diseases
- Reference integrity validation (all references resolve within the bundle)

### 6.2 PDF Export

Uses reportlab for clinical-grade PDF reports with autoimmune-specific formatting including autoantibody panel tables, HLA association summaries, disease activity score visualizations, and therapy recommendations.

### 6.3 Markdown

Structured plain-text reports for integration with clinical notes and downstream systems.

---

## 7. Sample Patient Data

Nine fully specified demo patients covering the disease spectrum:

| # | Patient | Disease | Key Features |
|---|---|---|---|
| 1 | Sarah Mitchell | SLE | ANA 1:640, anti-dsDNA+, lupus nephritis Class IV, 27 clinical PDFs, 3-year diagnostic odyssey |
| 2 | Linda Chen | Sjogren's | Anti-SSA/SSB+, HLA-DRB1*03:01, parotid involvement |
| 3 | Maya Rodriguez | POTS / hEDS / MCAS | Dysautonomia diagnostic odyssey, tilt table, hypermobility, mast cell mediators |
| 4 | David Park | AS | HLA-B*27:05 (OR 87.4), BASDAI scoring, sacroiliitis |
| 5 | Rachel Thompson | MCTD | Anti-U1 RNP+, Raynaud's, swollen fingers, overlap features (SLE/SSc/myositis) |
| 6 | Emma Williams | MS | HLA-DRB1*15:01, EDSS scoring, oligoclonal bands, DMT selection |
| 7 | James Cooper | T1D + Celiac | Anti-GAD65+, anti-tTG+, HLA-DQB1*02:01, overlap syndrome |
| 8 | Karen Foster | SSc | Anti-Scl-70+, modified Rodnan skin score, ILD monitoring |
| 9 | Michael Torres | Graves | TSH receptor Ab+, HLA-DRB1*03:01, antithyroid drug decision |

---

## 8. Performance Benchmarks

### 8.1 Clinical Engine Performance

| Engine | Latency |
|---|---|
| AutoantibodyInterpreter | <50 ms |
| HLAAssociationAnalyzer | <20 ms |
| DiseaseActivityScorer (all instruments) | <100 ms |
| FlarePredictor | <50 ms |
| BiologicTherapyAdvisor | <50 ms |
| DiagnosticOdysseyAnalyzer | <200 ms |
| **All 6 engines combined** | **<400 ms** |

### 8.2 RAG Pipeline Performance

| Operation | Latency |
|---|---|
| BGE-small query embedding | ~5 ms |
| 14-collection parallel search (top-5 each) | 10-18 ms |
| Claude Sonnet 4.6 generation | 15-25 sec |
| **Full RAG query end-to-end** | **~20-30 sec** |

### 8.3 Document Processing Performance

| Operation | Latency |
|---|---|
| PDF page extraction | ~200 ms/page |
| Text chunking (2,500 chars, 200 overlap) | <50 ms |
| BGE-small embedding per chunk | ~5 ms |
| Milvus insert per batch | <500 ms |
| **Full PDF ingestion per page** | **~2 sec** |

### 8.4 Export Performance

| Format | Latency |
|---|---|
| FHIR R4 Bundle + validation | <500 ms |
| PDF (reportlab) | <2 sec |
| Markdown | <100 ms |

---

## 9. Infrastructure

### 9.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4 |
| Embeddings | BGE-small-en-v1.5 (BAAI) -- 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit |
| REST API | FastAPI + Uvicorn |
| Configuration | Pydantic BaseSettings (AUTO_ prefix) |
| Testing | pytest |
| Export | FHIR R4, reportlab (PDF), Markdown |
| Document Processing | PyPDF2 |
| Containerization | Docker + Docker Compose |

### 9.2 Service Ports

| Service | Port |
|---|---|
| Streamlit UI | 8531 |
| FastAPI REST API | 8532 |
| Milvus (shared) | 19530 |

### 9.3 Dependencies on HCLS AI Factory

| Dependency | Type |
|---|---|
| Milvus 2.4 | Shared vector database (port 19530) |
| `genomic_evidence` collection | Read-only shared collection from Stage 2 RAG pipeline |
| BGE-small-en-v1.5 | Shared embedding model |
| Claude API key | Shared Anthropic API key |

---

## 10. Demo Scenarios

### 10.1 Validated Demo Queries

**Scenario 1 -- Autoantibody Interpretation:**
```
Patient: Sarah Mitchell, SLE
ANA: 1:640 homogeneous, anti-dsDNA: positive, anti-Smith: positive
Question: "What do these autoantibody results mean for my diagnosis?"
Expected: ANA 1:640 + anti-dsDNA + anti-Smith = pathognomonic for SLE. High specificity pattern.
```

**Scenario 2 -- HLA Risk Assessment:**
```
Patient: David Park, AS
HLA-B*27:05: positive
Question: "What is my risk for ankylosing spondylitis?"
Expected: HLA-B*27:05 OR 87.4 for AS. Strongest known HLA-disease association. Combined with inflammatory back pain = high clinical suspicion.
```

**Scenario 3 -- Flare Prediction:**
```
Patient: Sarah Mitchell, SLE
Anti-dsDNA: rising, C3/C4: falling, CRP: stable
Question: "Am I at risk for a lupus flare?"
Expected: Rising anti-dsDNA with falling complements = HIGH flare risk. Classic pre-flare pattern in SLE.
```

**Scenario 4 -- Overlap Therapy Selection:**
```
Patient: Rachel Thompson, MCTD
Status: Overlapping SLE/SSc/myositis features, anti-U1 RNP positive
Question: "What therapy approach is best for my overlap syndrome?"
Expected: Immunosuppressive therapy targeting predominant features. Consider organ-specific biologic options.
```

**Scenario 5 -- Overlap Syndrome:**
```
Patient: James Cooper, T1D + Celiac
HLA-DQB1*02:01: positive
Question: "Why do I have both Type 1 Diabetes and Celiac Disease?"
Expected: HLA-DQ2 shared risk. 10% co-occurrence rate. Polyautoimmune syndrome.
```

---

## 11. File Structure (Actual)

```
precision_autoimmune_agent/
├── src/                          # 9 core modules (~6,500 lines)
│   ├── models.py                 # Pydantic data models (autoantibodies, HLA, scores)
│   ├── collections.py            # 14 Milvus collection schemas + manager
│   ├── rag_engine.py             # Multi-collection RAG engine (parallel search + Claude)
│   ├── agent.py                  # AutoimmuneAgent orchestrator (6 engines)
│   ├── knowledge.py              # Knowledge base v2.0.0
│   ├── diagnostic_engine.py      # ACR/EULAR criteria, diagnostic odyssey, differentials
│   ├── document_processor.py     # PDF ingestion, chunking, embedding
│   ├── timeline_builder.py       # Diagnostic timeline construction
│   └── export.py                 # FHIR R4 + PDF + Markdown
├── app/
│   └── autoimmune_ui.py          # Streamlit (10 tabs)
├── api/
│   └── main.py                   # FastAPI REST server (14 endpoints)
├── config/
│   ├── settings.py               # AutoimmuneSettings (AUTO_ prefix)
│   └── logging.py                # Centralized logging
├── demo_data/                    # 9 patient directories with clinical PDFs
├── scripts/                      # setup_collections, generate_demo_patients
├── tests/                        # 455 tests (8 test files)
├── docker-compose.yml
├── Dockerfile
└── requirements.txt
```

**44 Python files | ~10,500 lines | Apache 2.0**

---

## 12. Implementation Status

- **Phase 1 (Architecture)** -- Complete. All data models, collection schemas, 6 clinical engines, RAG engine, and agent orchestrator implemented.
- **Phase 2 (Knowledge)** -- Complete. Knowledge base v2.0.0 with 22 HLA alleles, 24 autoantibodies, 22 biologics, 20 activity scores, 13 flare patterns.
- **Phase 3 (Integration)** -- Complete. Full RAG pipeline with Claude. FHIR R4 export with structural validation. Document ingestion pipeline.
- **Phase 4 (Testing)** -- Complete. 455 unit tests across 8 test files. Production readiness validation.
- **Phase 5 (Demo Ready)** -- Complete. Production-quality Streamlit UI with 10 tabs. 9 demo patients with clinical document PDFs.

---

## 13. Relationship to HCLS AI Factory

The Precision Autoimmune Intelligence Agent is the **fifth intelligence agent** in the HCLS AI Factory platform, joining:

1. **CAR-T Intelligence Agent** (8521/8522) -- Cross-functional CAR-T cell therapy intelligence
2. **Imaging Intelligence Agent** (8523/8524) -- Medical imaging detection, segmentation, and triage
3. **Precision Oncology Agent** (8525/8526) -- Tumor-specific treatment selection and clinical trial matching
4. **Precision Biomarker Agent** (8528/8529) -- Genomics-informed biomarker interpretation
5. **Precision Autoimmune Agent** (8531/8532) -- Autoimmune disease intelligence (this agent)

All agents share the same infrastructure (Milvus, BGE-small embeddings, Claude API) and can cross-reference the shared `genomic_evidence` collection.

---

## 14. Credits

- **Adam Jones**
- **Apache 2.0 License**
