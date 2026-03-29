# Precision Autoimmune Intelligence Agent -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [System Diagram](#1-system-diagram)
2. [Component Interactions](#2-component-interactions)
3. [Data Flow](#3-data-flow)
4. [Collection Design Rationale](#4-collection-design-rationale)
5. [Clinical Engines](#5-clinical-engines)
6. [Disease Activity Scoring](#6-disease-activity-scoring)
7. [Flare Prediction Engine](#7-flare-prediction-engine)
8. [Diagnostic Odyssey Engine](#8-diagnostic-odyssey-engine)
9. [Query Expansion](#9-query-expansion)
10. [RAG Pipeline](#10-rag-pipeline)
11. [Agent Orchestrator](#11-agent-orchestrator)
12. [Data Model Architecture](#12-data-model-architecture)

---

## 1. System Diagram

### 1.1 Full System Architecture

```
                          EXTERNAL USERS
                               |
                    +----------+----------+
                    |                     |
              +-----+------+      +------+-----+
              | Streamlit  |      | REST API   |
              | UI :8531   |      | :8532      |
              +-----+------+      +------+-----+
                    |                     |
                    +----------+----------+
                               |
                    +----------+----------+
                    |  AutoimmuneAgent   |
                    |  (src/agent.py)    |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Diagnostic  |    | Disease       |    | Biologic      |
   | Engine      |    | Activity      |    | Therapy       |
   |             |    | Scoring       |    | Advisor       |
   +------+------+    +-------+-------+    +-------+-------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Timeline    |    | Flare         |    | Document      |
   | Builder     |    | Prediction    |    | Processor     |
   +------+------+    +-------+-------+    +-------+-------+
          |                    |                    |
          +--------------------+--------------------+
                               |
                    +----------+----------+
                    |   RAG Engine       |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Knowledge   |    | Milvus        |    | LLM           |
   | Base        |    | Vector DB     |    | (Claude 4.6)  |
   | HLA, Ab,    |    | 14 Collections|    |               |
   | Biologics   |    |               |    |               |
   +-------------+    +-------+-------+    +---------------+
                               |
                    +----------+----------+
                    |  etcd    |  MinIO   |
                    +----------+----------+
```

### 1.2 PDF Ingestion Pipeline

```
  +------------------------------------------+
  |        Patient Clinical PDFs              |
  +------------------------------------------+
  | Progress Notes | Lab Reports | Imaging    |
  | Pathology | Genetic Tests | Referrals     |
  | Medication Lists                          |
  +---+--------+--------+--------+--------+--+
      |        |        |        |        |
  +---v-------------------------------------------+
  |         DocumentProcessor                      |
  |  - Text extraction (PyPDF2, page-by-page)      |
  |  - Sentence-boundary chunking                  |
  |    (max 2500 chars, 200 char overlap)          |
  |  - Document type classification (7 types)      |
  |  - Medical specialty detection (11 specialties)|
  |  - Entity extraction (29 autoantibodies,       |
  |    45 lab test patterns)                       |
  |  - Provider and date extraction                |
  +---+-------------------------------------------+
      |
  +---v------+
  | Embedding |
  | BGE-small |
  | 384-dim   |
  +---+------+
      |
  +---v--------------+
  | Milvus Insert     |
  | clinical_documents|
  | patient_labs      |
  | patient_timelines |
  +-----------------+
```

---

## 2. Component Interactions

### 2.1 Component Dependency Graph

```
AutoimmuneUI (Streamlit) ──> FastAPI Server ──> AutoimmuneAgent
                                                     |
                                 +-------------------+-------------------+
                                 |                   |                   |
                          DiagnosticEngine    DiseaseActivityScoring  BiologicAdvisor
                                 |                   |                   |
                          TimelineBuilder      FlarePrediction    DocumentProcessor
                                 |                   |                   |
                                 +-------------------+-------------------+
                                                     |
                                                RAGEngine
                                                     |
                                 +-------------------+-------------------+
                                 |                   |                   |
                            Milvus DB          Knowledge Base     genomic_evidence
                            (14 cols)          (HLA, Ab, Rx)      (shared col)
```

### 2.2 Module Responsibilities

| Module | File | Responsibilities |
|--------|------|-----------------|
| **AutoimmuneAgent** | `src/agent.py` | 5-step analyze_patient pipeline; autoantibody interpretation, HLA analysis, activity scoring, flare prediction, biologic recommendation |
| **DiagnosticEngine** | `src/diagnostic_engine.py` | 10 ACR/EULAR classification criteria evaluation, differential diagnosis, diagnostic odyssey analysis, 9 overlap syndrome detection |
| **AutoimmuneRAGEngine** | `src/rag_engine.py` | 14-collection parallel search, weighted scoring, disease area detection, knowledge augmentation, streaming LLM synthesis |
| **DocumentProcessor** | `src/document_processor.py` | PDF ingestion, chunking, document type classification (7 types), medical specialty detection (11 specialties), entity extraction |
| **TimelineBuilder** | `src/timeline_builder.py` | Event classification (12 types), date extraction, chronological timeline construction, days-from-first-symptom calculation |
| **AutoimmuneExporter** | `src/export.py` | Markdown, PDF (ReportLab), FHIR R4 DiagnosticReport Bundle |
| **Knowledge Base** | `src/knowledge.py` | HLA alleles, autoantibodies, biologic therapies, flare patterns, disease activity scoring |
| **Collections** | `src/collections.py` | 14 Milvus collection schemas, CRUD, parallel search |
| **Models** | `src/models.py` | Pydantic data models, enums, validation |

### 2.3 Interface Contracts

**AutoimmuneAgent inputs/outputs:**
```
Input:  PatientProfile(diseases, autoantibodies, hla_alleles, labs, medications, genotypes)
Output: AnalysisResult(activity_scores, flare_predictions, hla_associations, biologic_recommendations, alerts)
```

**DiagnosticEngine inputs/outputs:**
```
Input:  PatientProfile + clinical findings
Output: DiagnosticResult(criteria_evaluation, differential_diagnosis, odyssey_analysis, overlap_syndromes)
```

**Disease Activity inputs/outputs:**
```
Input:  Disease type + clinical measurements (e.g., joint counts, lab values)
Output: ActivityScore(score_name, score_value, interpretation, components, thresholds)
```

---

## 3. Data Flow

### 3.1 Patient Analysis Pipeline

```
Step 1: RECEIVE PATIENT PROFILE
  PatientProfile arrives via API or UI
  Demographics, diseases, autoantibodies, HLA, labs, medications
  |
Step 2: AUTOANTIBODY INTERPRETATION (agent.py)
  Match against 24 autoantibody-disease associations
  Calculate sensitivity/specificity for each antibody
  Identify disease-specific patterns (e.g., anti-dsDNA + low complement = active SLE)
  |
Step 3: HLA ANALYSIS (agent.py)
  Map 22 HLA alleles to disease associations
  Calculate odds ratios for identified alleles
  Flag high-risk allele-disease combinations
  (e.g., HLA-B*27:05 -> AS, OR 20-100x)
  |
Step 4: DISEASE ACTIVITY SCORING (agent.py)
  Execute 20 scoring systems across 13 diseases
  DAS28-ESR/CRP for RA (remission < 2.6, high > 5.1)
  SLEDAI-2K for SLE (inactive 0, mild 1-5, moderate 6-10, severe 11-19, very severe 20+)
  BASDAI for AS (active > 4.0)
  CDAI for IBD (remission < 150, severe > 450)
  EDSS for MS (0-10 disability scale)
  + 15 additional scoring systems
  |
Step 5: FLARE PREDICTION (agent.py)
  Evaluate 13 disease-specific biomarker patterns
  Identify early warning signs from lab trends
  Calculate flare risk (low/moderate/high)
  Identify protective and contributing factors
  |
Step 6: BIOLOGIC THERAPY RECOMMENDATION (agent.py)
  Query 22-drug biologic database
  Filter by disease indication
  Apply PGx considerations
  Check contraindications
  Rank by evidence level and guideline concordance
  |
Step 7: RAG SEARCH (rag_engine.py)
  Parallel search across 14 collections
  Weighted scoring with disease area boosting
  Knowledge augmentation from built-in knowledge base
  |
Step 8: LLM SYNTHESIS (rag_engine.py)
  Claude Sonnet 4.6 generates evidence-grounded response
  System prompt: autoimmune domain expert
  Inline citations with source references
  |
Step 9: RESPONSE ASSEMBLY (agent.py)
  Package AnalysisResult with:
  - Disease activity scores
  - Flare predictions
  - HLA associations
  - Biologic recommendations
  - Critical alerts
  - Evidence citations
```

### 3.2 Document Ingestion Flow

```
Patient PDF (progress note, lab report, imaging, pathology, etc.)
    |
    v
[Text Extraction] -- PyPDF2, page-by-page
    |
    v
[Chunking] -- Sentence boundaries, max 2500 chars, 200 char overlap
    |
    v
[Classification]
    |-- Document type: 7 types (lab report, progress note, imaging,
    |   pathology, genetic, referral, medication list)
    |-- Medical specialty: 11 specialties (rheumatology, neurology,
    |   endocrinology, dermatology, gastroenterology, etc.)
    |
    v
[Entity Extraction]
    |-- 29 autoantibody names detected
    |-- 45 lab test patterns extracted
    |-- Provider and date extraction
    |
    v
[Embedding] -- BGE-small-en-v1.5, 384-dim
    |
    v
[Milvus Insert]
    |-- autoimmune_clinical_documents (full chunks)
    |-- autoimmune_patient_labs (extracted lab values)
    |-- autoimmune_patient_timelines (timeline events)
```

---

## 4. Collection Design Rationale

### 4.1 Collection Inventory

The agent manages 14 specialized vector collections. All use COSINE similarity
with IVF_FLAT indexing (nlist=1024, nprobe=16) and 384-dimensional vectors
from BGE-small-en-v1.5.

| # | Collection Name | Description | Weight | Key Fields |
|---|----------------|-------------|--------|------------|
| 1 | autoimmune_clinical_documents | Ingested patient records (PDFs) | 0.18 | patient_id, doc_type, specialty, visit_date |
| 2 | autoimmune_patient_labs | Lab results with flag analysis | 0.14 | patient_id, test_name, value, unit, flag |
| 3 | autoimmune_autoantibody_panels | Autoantibody reference (24 antibodies) | 0.12 | antibody_name, associated_diseases, sensitivity, specificity |
| 4 | autoimmune_hla_associations | HLA allele-disease risk (22 alleles) | 0.08 | allele, disease, odds_ratio, population |
| 5 | autoimmune_disease_criteria | ACR/EULAR classification (10 sets) | 0.08 | disease, criteria_set, required_score |
| 6 | autoimmune_disease_activity | Activity scoring reference (20 systems) | 0.07 | score_name, disease, components, thresholds |
| 7 | autoimmune_flare_patterns | Flare prediction (13 patterns) | 0.06 | disease, biomarker_pattern, early_warning_signs |
| 8 | autoimmune_biologic_therapies | Biologic drug database (22 drugs) | 0.06 | drug_name, drug_class, mechanism, pgx_considerations |
| 9 | autoimmune_pgx_rules | Pharmacogenomic dosing rules | 0.04 | gene, variant, drug, phenotype, recommendation |
| 10 | autoimmune_clinical_trials | Autoimmune clinical trials | 0.05 | title, nct_id, phase, disease, intervention |
| 11 | autoimmune_literature | Published literature | 0.05 | title, journal, year, pmid, disease_focus |
| 12 | autoimmune_patient_timelines | Diagnostic odyssey timelines | 0.03 | patient_id, event_type, event_date, days_from_first_symptom |
| 13 | autoimmune_cross_disease | Overlap syndromes (9 patterns) | 0.02 | primary_disease, associated_conditions, shared_pathways |
| 14 | genomic_evidence | Shared VCF-derived genomic variants (read-only) | 0.02 | (managed by HCLS AI Factory core) |
| | **Total** | | **1.00** | |

### 4.2 Collection Weight Distribution

```
autoimmune_clinical_docs   ██████████████████████  0.18
autoimmune_patient_labs    █████████████████       0.14
autoimmune_autoantibody    ██████████████          0.12
autoimmune_hla_assoc       ██████████              0.08
autoimmune_disease_crit    ██████████              0.08
autoimmune_disease_act     █████████               0.07
autoimmune_flare_pat       ████████                0.06
autoimmune_biologic_rx     ████████                0.06
autoimmune_clinical_tri    ███████                 0.05
autoimmune_literature      ███████                 0.05
autoimmune_pgx_rules       █████                   0.04
autoimmune_timelines       ████                    0.03
autoimmune_cross_disease   ███                     0.02
genomic_evidence           ███                     0.02
                                              Sum: 1.00
```

### 4.3 Why 14 Collections

The multi-collection design serves autoimmune medicine's inherent complexity:

**1. Semantic Precision**
- "Positive": autoantibody titer (autoimmune_autoantibody_panels), disease activity (autoimmune_disease_activity), trial outcome (autoimmune_clinical_trials)
- "Score": DAS28 (activity), SLEDAI-2K (activity), ACR criteria (diagnosis)
- "Flare": biomarker-predicted (autoimmune_flare_patterns), documented clinical (autoimmune_clinical_documents)

**2. Disease-Specific Boosting**
When lupus is detected in a query, the engine boosts autoimmune_autoantibody_panels and autoimmune_disease_activity while reducing autoimmune_biologic_therapies weight.

**3. Independent Update Cycles**
- Biologic therapies: Updated when new FDA approvals or EULAR/ACR guideline revisions occur
- Clinical trials: Monthly automated refresh from ClinicalTrials.gov
- Patient documents: Real-time ingestion via PDF upload

---

## 5. Clinical Engines

### 5.1 AutoimmuneAgent (src/agent.py)

The main orchestrator runs the 5-step `analyze_patient` pipeline:

1. **Interpret autoantibody panel** against 24 antibody-disease associations with sensitivity/specificity data
2. **Analyze HLA profile** against 22 allele-disease associations with odds ratios
3. **Calculate disease activity scores** using 20 scoring systems across 13 diseases
4. **Predict flare risk** using 13 disease-specific biomarker patterns
5. **Recommend biologic therapies** from the 22-drug database with PGx filtering

Cross-agent integration stubs connect to:
- Biomarker Agent (inflammation context, CRP trends)
- Imaging Agent (joint/organ assessment, MRI findings)

### 5.2 DiagnosticEngine (src/diagnostic_engine.py)

Clinical reasoning engine implementing:

- **Classification criteria evaluation:** 10 ACR/EULAR and diagnostic criteria sets with point-based scoring
- **Differential diagnosis generation:** Ranks diseases by combining autoantibody specificity scores and HLA odds ratios (log2-scaled)
- **Diagnostic odyssey analysis:** Calculates time from first symptom to diagnosis, tracks specialists seen, identifies misdiagnoses and turning points
- **Overlap syndrome detection:** 9 defined cross-disease patterns including MCTD, POTS/hEDS/MCAS triad, T1D-celiac overlap, lupus-APS overlap

### 5.3 DocumentProcessor (src/document_processor.py)

PDF ingestion pipeline with:
- Text extraction via PyPDF2 (page-by-page)
- Sentence-boundary chunking (max 2500 chars, 200 char overlap)
- Document type classification (7 types)
- Medical specialty detection (11 specialties)
- Entity extraction: 29 autoantibody names, 45 lab test patterns
- Provider and date extraction from text

### 5.4 TimelineBuilder (src/timeline_builder.py)

Diagnostic odyssey construction:
- Event classification from text using 12 event type patterns (symptom onset, diagnosis, misdiagnosis, lab result, imaging, biopsy, genetic test, treatment start, treatment change, flare, referral, ER visit)
- Date extraction from 4 date formats
- Chronological timeline building with days-from-first-symptom calculation
- Milvus-ready record generation for the patient timelines collection

---

## 6. Disease Activity Scoring

### 6.1 Scoring Systems (20 Total)

The agent implements 20 validated disease activity scoring systems across all 13 supported autoimmune conditions:

| Disease | Score | Components | Thresholds |
|---------|-------|-----------|------------|
| Rheumatoid Arthritis | DAS28-ESR | Tender/swollen joints (28), ESR, patient VAS | Remission <2.6, Low 2.6-3.2, Moderate 3.2-5.1, High >5.1 |
| Rheumatoid Arthritis | DAS28-CRP | Tender/swollen joints (28), CRP, patient VAS | Same thresholds as DAS28-ESR |
| Rheumatoid Arthritis | CDAI | Tender/swollen joints (28), physician/patient VAS | Remission <=2.8, Low 2.8-10, Moderate 10-22, High >22 |
| SLE | SLEDAI-2K | 24 weighted items (seizure, psychosis, vasculitis, etc.) | Inactive 0, Mild 1-5, Moderate 6-10, Severe 11-19, Very Severe 20+ |
| SLE | BILAG-2004 | 9 organ systems, A-E per system | Grade A=severe, B=moderate, C=mild, D=previous, E=never |
| Ankylosing Spondylitis | BASDAI | 6 questions (fatigue, pain, stiffness, tenderness) | Active >4.0 |
| Ankylosing Spondylitis | ASDAS | Back pain, patient global, stiffness, CRP | Inactive <1.3, Low 1.3-2.1, High 2.1-3.5, Very High >3.5 |
| Psoriatic Arthritis | DAPSA | Tender/swollen joints, CRP, patient pain/global | Remission <=4, Low 4-14, Moderate 14-28, High >28 |
| IBD (Crohn's) | CDAI | 8 variables over 7 days | Remission <150, Mild 150-220, Moderate 220-450, Severe >450 |
| IBD (UC) | Mayo Score | Stool frequency, bleeding, endoscopy, physician | Remission 0-2, Mild 3-5, Moderate 6-10, Severe 11-12 |
| Multiple Sclerosis | EDSS | 8 functional systems | 0 (normal) to 10 (death), 0.5 increments |
| Scleroderma | mRSS | Skin thickness at 17 sites | 0-51, mild <14, moderate 14-29, severe >=30 |
| Sjogren's | ESSDAI | 12 organ domains | Low <5, Moderate 5-13, High >=14 |
| Myositis | MMT-8 | Manual muscle testing 8 groups | 0-80 per side |
| Vasculitis | BVAS | 9 organ systems | Active >0 |
| Graves' Disease | Burch-Wartofsky | Temperature, CNS, GI, HR, CHF, AF, bilirubin | <25 unlikely, 25-44 impending, >=45 thyroid storm |

### 6.2 Scoring Architecture

```
PatientProfile (disease, labs, clinical measurements)
    |
    v
[Disease Identification]
    |-- Match to AutoimmuneDisease enum (13 conditions)
    |
    v
[Score Selection]
    |-- Select applicable scoring systems for identified diseases
    |-- Multiple scores per disease (e.g., DAS28-ESR + DAS28-CRP + CDAI for RA)
    |
    v
[Component Extraction]
    |-- Extract required lab values and clinical measurements
    |-- Validate completeness (flag missing components)
    |
    v
[Score Calculation]
    |-- Apply published formulas and coefficients
    |-- Classify into severity categories
    |
    v
[Interpretation]
    ActivityScore(score_name, score_value, category, interpretation,
                  components_used, missing_components, recommendation)
```

---

## 7. Flare Prediction Engine

### 7.1 Disease-Specific Biomarker Patterns (13)

The flare prediction engine monitors biomarker trends against 13 disease-specific
patterns to provide early warning of disease flares:

| Disease | Early Warning Biomarkers | Protective Factors |
|---------|------------------------|-------------------|
| RA | Rising CRP, rising ESR, declining hemoglobin | Stable CRP, normal ESR, therapeutic drug levels |
| SLE | Rising anti-dsDNA, falling C3/C4, rising proteinuria | Stable complement, negative anti-dsDNA |
| AS | Rising CRP, rising ESR, elevated fecal calprotectin | Normal CRP, stable BASDAI |
| Sjogren's | Rising IgG, falling C4, rising RF | Stable IgG, normal C4 |
| SSc | Rising NT-proBNP, declining DLCO, rising CRP | Stable DLCO, normal NT-proBNP |
| MS | Rising neurofilament light chain, new MRI lesions | Stable NfL, no new lesions |
| T1D | Rising HbA1c, declining C-peptide, rising anti-GAD | Stable HbA1c, detectable C-peptide |
| IBD | Rising fecal calprotectin, rising CRP, declining albumin | Normal calprotectin, stable CRP |
| Psoriatic Arthritis | Rising CRP, rising ESR, new nail changes | Normal CRP, stable skin scores |
| Vasculitis | Rising ANCA titers, rising CRP, declining eGFR | Stable ANCA, normal CRP, stable eGFR |
| Myositis | Rising CK, rising aldolase, rising anti-Jo-1 | Stable CK, normal aldolase |
| Graves' Disease | Rising TSI, falling TSH, rising free T4 | Normalizing TSH, declining TSI |
| Celiac | Rising tTG-IgA, rising DGP-IgG | Declining tTG-IgA, resolving symptoms |

### 7.2 Flare Risk Classification

```
Biomarker Trend Analysis
    |
    v
[Pattern Matching]
    |-- Count matching early warning signs
    |-- Count matching protective factors
    |
    v
[Risk Calculation]
    |-- Low risk: 0-1 warning signs, 2+ protective factors
    |-- Moderate risk: 2 warning signs or declining trend
    |-- High risk: 3+ warning signs or critical threshold crossed
    |
    v
FlareRiskResult(risk_level, contributing_factors, protective_factors,
                recommended_actions, monitoring_frequency)
```

---

## 8. Diagnostic Odyssey Engine

### 8.1 Architecture

The diagnostic odyssey engine analyzes the full clinical timeline to identify
patterns common in autoimmune disease diagnosis:

```
Patient Timeline (from TimelineBuilder)
    |
    v
[Odyssey Metrics]
    |-- Time from first symptom to diagnosis (months/years)
    |-- Number of specialists seen
    |-- Number of misdiagnoses
    |-- Number of ER visits before diagnosis
    |
    v
[Turning Point Identification]
    |-- Key test that led to diagnosis
    |-- Provider who made the diagnosis
    |-- Triggering clinical event
    |
    v
[Overlap Syndrome Detection]
    |-- 9 defined patterns:
    |   1. MCTD (Mixed Connective Tissue Disease)
    |   2. POTS/hEDS/MCAS triad
    |   3. T1D + Celiac overlap
    |   4. SLE + APS (lupus-antiphospholipid) overlap
    |   5. RA + Sjogren's overlap (secondary Sjogren's)
    |   6. SSc + Myositis overlap (scleromyositis)
    |   7. SLE + SSc overlap
    |   8. IBD + SpA overlap
    |   9. Thyroid + T1D autoimmune polyglandular
    |
    v
OdysseyResult(duration_months, specialists_count, misdiagnoses,
              turning_point, overlap_syndromes, recommendations)
```

---

## 9. Query Expansion

### 9.1 Disease Area Detection

The RAG engine detects disease areas from query keywords across 13 autoimmune conditions:

| Disease | Trigger Keywords |
|---------|-----------------|
| RA | rheumatoid, DAS28, CDAI, joint erosion, RF, anti-CCP |
| SLE | lupus, SLEDAI, anti-dsDNA, complement, nephritis |
| AS | ankylosing, BASDAI, HLA-B27, sacroiliitis, spondylitis |
| Sjogren's | sjogren, sicca, anti-SSA, anti-SSB, Schirmer |
| SSc | scleroderma, sclerosis, mRSS, Raynaud, anti-Scl-70 |
| MS | multiple sclerosis, EDSS, oligoclonal, demyelination |
| T1D | type 1 diabetes, anti-GAD, anti-IA-2, C-peptide |
| IBD | crohn, colitis, calprotectin, CDAI, Mayo score |
| Psoriatic | psoriatic, DAPSA, dactylitis, enthesitis |
| Vasculitis | ANCA, vasculitis, BVAS, granulomatosis |
| Myositis | myositis, CK, aldolase, anti-Jo-1, MMT-8 |
| Graves' | graves, TSI, thyroid storm, Burch-Wartofsky |
| Celiac | celiac, tTG-IgA, DGP, villous atrophy |

### 9.2 Knowledge Augmentation

When disease areas are detected, the RAG engine injects structured knowledge:

- **HLA context:** Allele-disease associations with odds ratios and PMIDs
- **Autoantibody context:** Sensitivity/specificity data, disease associations
- **Therapy context:** Biologic drug profiles with PGx considerations
- **Flare pattern context:** Disease-specific early warning biomarkers
- **Activity scoring context:** Scoring system components and thresholds

---

## 10. RAG Pipeline

### 10.1 RAG Architecture

```
Query --> Disease Area Detection --> Embed (BGE-small + instruction prefix)
  --> Parallel Search (14 collections, top_k=5 each, score_threshold=0.40)
  --> Weighted Scoring (per-collection weights, sum=1.0)
  --> Deduplication (by ID and text content hash)
  --> Knowledge Augmentation (HLA, autoantibody, therapy, flare patterns)
  --> LLM Synthesis (Claude Sonnet 4.6, system prompt, conversation memory)
  --> Streaming Response with Evidence Citations
```

### 10.2 Embedding Configuration

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
| Search mode | Asymmetric |
| Embedding cache | 256 entries, FIFO eviction |

### 10.3 Citation Scoring

| Score Range | Relevance Level |
|------------|----------------|
| >= 0.80 | High |
| >= 0.60 | Medium |
| < 0.60 | Low |

### 10.4 Conversation Memory

Thread-safe conversation memory using a deque with configurable maximum size (default: 3 turns). Enables contextual follow-up questions like "What about adding rituximab?" after an initial SLE query.

---

## 11. Agent Orchestrator

### 11.1 Key Statistics

| Metric | Value |
|--------|-------|
| Test suite | 455 tests across 8 files |
| Milvus collections | 14 (13 owned + 1 read-only) |
| Supported diseases | 13 autoimmune conditions |
| Biologic therapies | 22 drugs with PGx considerations |
| HLA alleles mapped | 22 with disease associations |
| Autoantibodies mapped | 24 with sensitivity/specificity data |
| Disease activity scores | 20 scoring systems |
| Flare prediction patterns | 13 disease-specific biomarker patterns |
| Classification criteria | 10 ACR/EULAR criteria sets |
| Overlap syndromes | 9 cross-disease patterns |
| Lab test patterns | 45 extractable lab values |
| Demo patients | 9 with full clinical PDF datasets |
| API endpoints | 14 REST endpoints |

### 11.2 Error Handling Strategy

The orchestrator implements graceful degradation:

1. **Milvus unavailable**: Returns error with clear diagnostic message
2. **LLM unavailable**: Returns search results without synthesis
3. **Scoring engine error**: Logs warning, returns available scores, skips failed calculator
4. **Flare prediction error**: Logs warning, returns response without flare analysis
5. **Document processor error**: Logs error, returns partial ingestion results
6. **Timeout**: Returns partial results with timeout indicator

### 11.3 API Endpoints Summary

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/healthz` | GET | Health check |
| `/readyz` | GET | Readiness check (Milvus connection) |
| `/v1/query` | POST | RAG query with streaming support |
| `/v1/analyze` | POST | Full patient analysis (5-step pipeline) |
| `/v1/autoantibody` | POST | Autoantibody panel interpretation |
| `/v1/hla` | POST | HLA profile analysis |
| `/v1/activity-score` | POST | Disease activity score calculation |
| `/v1/flare-predict` | POST | Flare risk prediction |
| `/v1/biologic-recommend` | POST | Biologic therapy recommendation |
| `/v1/criteria-eval` | POST | ACR/EULAR classification criteria evaluation |
| `/v1/differential` | POST | Differential diagnosis generation |
| `/v1/odyssey` | POST | Diagnostic odyssey analysis |
| `/v1/ingest` | POST | Document ingestion (PDF upload) |
| `/v1/export` | POST | Report export (Markdown/PDF/FHIR) |

### 11.4 Port Configuration

| Service | Port | Protocol |
|---------|------|----------|
| Streamlit UI | 8531 | HTTP |
| FastAPI REST API | 8532 | HTTP |
| Milvus (shared) | 19530 | gRPC |

---

## 12. Data Model Architecture

### 12.1 Demo Patients

| # | Name | Conditions | Key Features |
|---|------|-----------|-------------|
| 1 | Sarah Mitchell, 34F | SLE (lupus nephritis) | ANA 1:640, anti-dsDNA, anti-Smith, low complement |
| 2 | Maya Rodriguez, 28F | POTS/hEDS/MCAS | Dysautonomia diagnostic odyssey |
| 3 | Linda Chen, 45F | Sjogren's | Anti-SSA/Ro+, anti-SSB/La+, Schirmer test |
| 4 | David Park, 45M | AS | HLA-B*27:05, uveitis, sacroiliitis, 3-year odyssey |
| 5 | Rachel Thompson, 38F | MCTD | Mixed connective tissue disease, anti-U1 RNP |
| 6 | Emma Williams, 34F | MS (RRMS) | Optic neuritis, MRI lesions, oligoclonal bands |
| 7 | James Cooper, 19M | T1D + Celiac | Anti-GAD65, anti-IA-2, tTG-IgA+, HLA-DQ2/DQ8 |
| 8 | Karen Foster, 48F | SSc (dcSSc) | Anti-Scl-70+, Raynaud's, ILD, mRSS |
| 9 | Michael Torres, 41M | Graves' Disease | TSI+, Burch-Wartofsky scoring |

### 12.2 Export Pipeline

Three export formats:

- **Markdown:** Structured report with critical alerts, disease activity tables, flare predictions, HLA associations, biologic recommendations, evidence citations
- **PDF:** ReportLab with NVIDIA-themed colors (green #76B900 headers), tabular scores, clinical footer
- **FHIR R4:** DiagnosticReport Bundle with Patient, Observation (disease activity scores, flare predictions), and DiagnosticReport resources

### 12.3 File Structure

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
|-- demo_data/                   # 9 patient directories with clinical PDFs
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
|-- tests/                       # 8 test files, 455 tests
|
|-- docker-compose.yml           # 3 services: streamlit, api, setup
|-- Dockerfile
|-- requirements.txt             # 19 dependencies
|-- README.md
```

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Autoimmune Agent is a clinical decision support research tool for autoimmune disease assessment. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
