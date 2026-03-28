# Precision Autoimmune Intelligence Agent -- Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Precision Autoimmune Intelligence Agent extends the HCLS AI Factory platform to deliver comprehensive autoimmune disease intelligence. Unlike point-of-care tools that address a single disease or scoring system, this agent integrates autoantibody interpretation, HLA-disease associations, disease activity scoring, flare prediction, biologic therapy selection with pharmacogenomic considerations, and diagnostic odyssey analysis across 13 autoimmune diseases -- all grounded in a 14-collection RAG pipeline.

The agent combines **6 deterministic clinical analysis engines** with a **14-collection RAG pipeline** to answer questions like *"Interpret ANA 1:640 homogeneous with positive anti-dsDNA and falling complement"* -- simultaneously searching autoantibody reference data, HLA associations, disease activity thresholds, flare patterns, biologic therapy databases, and clinical evidence, then synthesizing a grounded response through Claude.

### Key Results

| Metric | Value |
|---|---|
| Unit tests passing | **455** |
| Milvus collections | **14** (13 owned + 1 read-only) |
| Autoimmune diseases | **13** in AutoimmuneDisease enum |
| Biologic therapies | **22** with PGx considerations |
| Autoantibodies mapped | **24** to disease associations with sensitivity/specificity |
| HLA alleles | **22** with disease odds ratios and PMIDs |
| Disease activity scores | **20** scoring systems across 13 diseases |
| Demo patients | **9** with full clinical PDF records |
| Flare patterns | **13** disease-specific biomarker patterns |
| Classification criteria | **10** ACR/EULAR criteria sets |
| Overlap syndromes | **9** cross-disease patterns |
| Knowledge version | **v2.0.0** |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Autoimmune Agent Role |
|---|---|
| **DataStore** | Clinical PDFs, reference JSON files, knowledge base (HLA, antibodies, biologics, flare patterns) |
| **DataEngine** | Document processor: PDF -> chunks -> BGE-small embedding -> Milvus insert |
| **DataBase** | 14 Milvus collections (13 owned + 1 read-only) + 9 demo patients |
| **InsightEngine** | 6 clinical engines + BGE-small embedding + multi-collection RAG |
| **AgentEngine** | AutoimmuneAgent orchestrator + Streamlit UI (10 tabs) + FastAPI REST |

### 2.2 System Diagram

```
                        ┌─────────────────────────────────┐
                        │    Streamlit UI (8531)            │
                        │    10 tabs: Query | Analysis |    │
                        │    Ingest | Odyssey | Antibody |  │
                        │    HLA | Activity | Flare |       │
                        │    Therapy | Knowledge             │
                        └──────────────┬──────────────────┘
                                       │
                        ┌──────────────▼──────────────────┐
                        │  AutoimmuneAgent                  │
                        │  Orchestrates 6 analysis engines  │
                        │  + RAG pipeline + export           │
                        └──────────────┬──────────────────┘
                                       │
            ┌──────────────────────────┼───────────────────────────┐
            │                          │                           │
  ┌─────────▼──────────┐   ┌──────────▼──────────┐   ┌──────────▼──────────┐
  │ Deterministic       │   │ RAG Pipeline         │   │ Export               │
  │ Analysis Engines    │   │                      │   │                      │
  │                     │   │ BGE-small-en-v1.5    │   │ FHIR R4 Bundle       │
  │ AutoantibodyInterp  │   │ (384-dim embedding)  │   │ PDF (reportlab)      │
  │ HLAAssociation      │   │         │            │   │ Markdown             │
  │ DiseaseActivity     │   │         ▼            │   │                      │
  │ FlarePredictor      │   │ Parallel Search      │   │ + FHIR Validation    │
  │ BiologicTherapy     │   │ 14 Milvus Collections│   │                      │
  │ DiagnosticOdyssey   │   │ (ThreadPoolExecutor) │   │                      │
  │                     │   │         │            │   │                      │
  │ + DocumentProcessor │   │         ▼            │   │                      │
  │ + TimelineBuilder   │   │ Claude Sonnet 4      │   │                      │
  └─────────────────────┘   └──────────────────────┘   └──────────────────────┘
            │                          │
  ┌─────────▼──────────────────────────▼──────────────────────────────┐
  │                  Milvus 2.4 -- 14 Collections                      │
  │                                                                    │
  │  autoimmune_clinical_documents       autoimmune_patient_labs       │
  │  autoimmune_autoantibody_panels      autoimmune_hla_associations   │
  │  autoimmune_disease_criteria         autoimmune_disease_activity   │
  │  autoimmune_flare_patterns           autoimmune_biologic_therapies │
  │  autoimmune_pgx_rules               autoimmune_clinical_trials     │
  │  autoimmune_literature               autoimmune_patient_timelines  │
  │  autoimmune_cross_disease            genomic_evidence [read-only]  │
  └───────────────────────────────────────────────────────────────────┘
```

---

## 3. Data Collections -- Actual State

### 3.1 `autoimmune_clinical_documents`

Ingested clinical documents (PDFs) chunked for vector search.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Document text chunk |
| patient_id | VARCHAR(64) | Patient identifier |
| doc_type | VARCHAR(128) | progress_note, lab_report, imaging, pathology, etc. |
| specialty | VARCHAR(128) | rheumatology, neurology, nephrology, etc. |
| provider | VARCHAR(256) | Provider name |
| visit_date | VARCHAR(32) | ISO date |
| source_file | VARCHAR(512) | Original PDF filename |
| page_number | INT64 | Page number in source |
| chunk_index | INT64 | Chunk index within document |

### 3.2 `autoimmune_patient_labs`

Laboratory results with reference range analysis.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Lab result description |
| patient_id | VARCHAR(64) | Patient identifier |
| test_name | VARCHAR(256) | Lab test name |
| value | FLOAT | Measured value |
| unit | VARCHAR(64) | Unit of measurement |
| reference_range | VARCHAR(128) | Normal reference range |
| flag | VARCHAR(32) | normal, high, low, critical |
| collection_date | VARCHAR(32) | Date of collection |
| panel_name | VARCHAR(256) | Lab panel name |

### 3.3 `autoimmune_autoantibody_panels` -- 24 records

Autoantibody reference data with disease associations.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Autoantibody description |
| antibody_name | VARCHAR(128) | Antibody name (ANA, anti-dsDNA, RF, etc.) |
| associated_diseases | VARCHAR(1024) | Associated autoimmune diseases |
| sensitivity | FLOAT | Diagnostic sensitivity (0-1) |
| specificity | FLOAT | Diagnostic specificity (0-1) |
| pattern | VARCHAR(128) | Staining pattern (homogeneous, speckled, etc.) |
| clinical_significance | VARCHAR(2000) | Clinical interpretation |
| interpretation_guide | VARCHAR(2000) | Interpretation guidelines |

### 3.4 `autoimmune_hla_associations` -- 22 records

HLA allele to disease risk mapping.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Association description |
| allele | VARCHAR(64) | HLA allele (e.g., HLA-B*27:05) |
| disease | VARCHAR(256) | Associated disease |
| odds_ratio | FLOAT | Disease odds ratio |
| population | VARCHAR(128) | Population studied |
| pmid | VARCHAR(32) | PubMed ID |
| mechanism | VARCHAR(1024) | Pathogenic mechanism |
| clinical_implication | VARCHAR(2000) | Clinical implications |

### 3.5 `autoimmune_disease_criteria` -- 10 records

ACR/EULAR classification criteria for autoimmune diseases.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Criteria description |
| disease | VARCHAR(256) | Disease name |
| criteria_set | VARCHAR(256) | Criteria set name (e.g., 2010 ACR/EULAR RA) |
| criteria_type | VARCHAR(64) | classification or diagnostic |
| year | INT64 | Publication year |
| required_score | VARCHAR(128) | Score threshold |
| criteria_items | VARCHAR(3000) | Individual criteria items |
| sensitivity_specificity | VARCHAR(256) | Validation metrics |

### 3.6 `autoimmune_disease_activity` -- 20 records

Disease activity scoring systems with thresholds and components.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Scoring system description |
| score_name | VARCHAR(128) | Score name (DAS28-CRP, SLEDAI-2K, etc.) |
| disease | VARCHAR(256) | Applicable disease |
| components | VARCHAR(2000) | Component list |
| thresholds | VARCHAR(1024) | JSON thresholds for remission/low/moderate/high |
| interpretation | VARCHAR(2000) | Interpretation guide |
| monitoring_frequency | VARCHAR(512) | Recommended monitoring schedule |

### 3.7 `autoimmune_flare_patterns` -- 13 records

Flare prediction biomarker patterns for each disease.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Pattern description |
| disease | VARCHAR(256) | Disease name |
| biomarker_pattern | VARCHAR(2000) | Early warning biomarker list |
| early_warning_signs | VARCHAR(2000) | Early warning descriptions |
| typical_timeline | VARCHAR(512) | Typical flare timeline |
| protective_factors | VARCHAR(1024) | Factors reducing flare risk |
| intervention_triggers | VARCHAR(1024) | When to escalate treatment |

### 3.8 `autoimmune_biologic_therapies` -- 22 records

Biologic therapy database with PGx considerations.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Therapy description |
| drug_name | VARCHAR(128) | Drug name |
| drug_class | VARCHAR(128) | Drug class (TNF inhibitor, IL-6R inhibitor, etc.) |
| mechanism | VARCHAR(512) | Mechanism of action |
| indicated_diseases | VARCHAR(1024) | Approved indications |
| pgx_considerations | VARCHAR(2000) | Pharmacogenomic factors |
| contraindications | VARCHAR(1024) | Contraindications |
| monitoring | VARCHAR(2000) | Monitoring requirements |
| dosing | VARCHAR(512) | Dosing information |
| evidence_level | VARCHAR(64) | Evidence level |

### 3.9 `autoimmune_pgx_rules`

Pharmacogenomic dosing rules for autoimmune therapies.

### 3.10 `autoimmune_clinical_trials`

Autoimmune disease clinical trials with biomarker-based eligibility criteria.

### 3.11 `autoimmune_literature`

Published autoimmune research literature with abstracts and disease focus.

### 3.12 `autoimmune_patient_timelines`

Patient diagnostic timeline events for odyssey analysis.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Event description |
| patient_id | VARCHAR(64) | Patient identifier |
| event_type | VARCHAR(128) | symptom_onset, diagnosis, treatment_start, flare, etc. |
| event_date | VARCHAR(32) | ISO date |
| description | VARCHAR(2000) | Event description |
| provider | VARCHAR(256) | Provider name |
| specialty | VARCHAR(128) | Medical specialty |
| days_from_first_symptom | INT64 | Days since first symptom |

### 3.13 `autoimmune_cross_disease` -- 9 records

Cross-disease overlap syndromes and shared pathogenic mechanisms.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(128) | Primary key |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 |
| text_chunk | VARCHAR(3000) | Overlap description |
| primary_disease | VARCHAR(256) | Primary disease |
| associated_conditions | VARCHAR(1024) | Associated conditions |
| shared_pathways | VARCHAR(1024) | Shared pathogenic mechanisms |
| shared_biomarkers | VARCHAR(1024) | Shared biomarkers |
| overlap_criteria | VARCHAR(2000) | Overlap diagnostic criteria |
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

Maps 24 autoantibodies to disease associations with sensitivity and specificity data.

**24 Supported Autoantibodies:**

| Autoantibody | Primary Disease Association | Sensitivity | Specificity |
|---|---|---|---|
| ANA | SLE, Sjogren's, SSc | 0.95 (SLE) | 0.65 |
| anti-dsDNA | SLE | 0.70 | 0.95 |
| anti-Smith | SLE | 0.25 | 0.99 |
| RF | RA, Sjogren's | 0.70 (RA) | 0.85 |
| anti-CCP | RA | 0.67 | 0.95 |
| anti-Scl-70 | SSc (diffuse) | 0.35 | 0.99 |
| anti-centromere | SSc (limited/CREST) | 0.40 | 0.98 |
| anti-SSA/Ro | Sjogren's, SLE | 0.70 (SS) | 0.90 |
| anti-SSB/La | Sjogren's | 0.40 | 0.95 |
| anti-Jo-1 | Antisynthetase syndrome | 0.30 | 0.99 |
| AChR antibody | Myasthenia Gravis | 0.85 | 0.99 |
| anti-tTG IgA | Celiac Disease | 0.93 | 0.97 |
| TSI | Graves' Disease | 0.90 | 0.95 |
| anti-TPO | Hashimoto's, Graves' | 0.90 (Hashi) | 0.85 |
| anti-RNP | MCTD, SLE | 0.95 (MCTD) | 0.85 |
| anti-histone | Drug-induced lupus, SLE | 0.95 (DIL) | 0.50 |
| ANCA (c-ANCA/PR3) | GPA (Wegener's) | 0.90 | 0.95 |
| ANCA (p-ANCA/MPO) | MPA, EGPA | 0.70 | 0.90 |
| anti-Pm-Scl | SSc-myositis overlap | 0.10 | 0.98 |
| anti-RNA Pol III | SSc (diffuse, renal crisis) | 0.20 | 0.99 |
| anti-cardiolipin | APS | 0.80 | 0.80 |
| lupus anticoagulant | APS | 0.55 | 0.95 |
| anti-beta2-GP I | APS | 0.70 | 0.90 |
| anti-MuSK | MG (AChR-negative) | 0.40 | 0.99 |

**Interpretation logic:** For each positive antibody in the patient's panel, the engine returns all disease associations with their sensitivity, specificity, titer, and staining pattern.

### 4.2 HLAAssociationAnalyzer

Matches 22 HLA alleles against a curated disease association database with odds ratios and PubMed references.

**22 HLA Alleles:**

| Allele | Disease | Odds Ratio | PMID |
|---|---|---|---|
| HLA-B*27:05 | Ankylosing Spondylitis | 87.4 | 25603694 |
| HLA-B*27:02 | Ankylosing Spondylitis | 50.0 | 25603694 |
| HLA-C*06:02 | Psoriasis | 10.0 | 23143594 |
| HLA-DQB1*02:01 | Celiac Disease | 7.0 | 17554300 |
| HLA-DQB1*03:02 | Type 1 Diabetes | 6.5 | 17554300 |
| HLA-B*51:01 | Behcet's Disease | 5.9 | 22704706 |
| HLA-DRB1*04:01 | Rheumatoid Arthritis | 4.2 | 20301572 |
| HLA-DRB1*04:04 | Rheumatoid Arthritis | 3.8 | 20301572 |
| HLA-DRB1*03:01 | T1D, SLE, Sjogren's, Graves', Celiac | 2.2-7.0 | 17554300 |
| HLA-DRB1*15:01 | Multiple Sclerosis | 3.1 | 21833088 |
| HLA-B*08:01 | Myasthenia Gravis, SLE | 3.4 (MG) | 16710306 |
| HLA-DRB1*04:05 | Rheumatoid Arthritis | 3.5 | 20301572 |
| HLA-DRB1*01:01 | RA, SSc | 2.0-2.1 | 20301572 |
| HLA-DRB1*08:01 | SLE | 2.1 | 19864127 |
| HLA-DQA1*05:01 | Celiac Disease | 7.0 | 17554300 |
| HLA-DRB1*15:03 | MS (African-descent) | 2.8 | 21833088 |
| HLA-A*02:01 | Type 1 Diabetes | 1.5 | 17554300 |
| HLA-DRB1*07:01 | T1D (PROTECTIVE) | 0.3 | 17554300 |
| HLA-DPB1*05:01 | Systemic Sclerosis | 2.3 | 24098041 |
| HLA-B*44:03 | Type 1 Diabetes | 1.4 | 17554300 |
| HLA-DRB1*11:01 | Sjogren's Syndrome | 2.5 | 19864127 |
| HLA-DRB1*13:01 | T1D (PROTECTIVE) | 0.2 | 17554300 |

**Matching logic:** Supports exact allele match and broad allele group matching (e.g., B*27:05 matches B*27). Results sorted by odds ratio (highest risk first). Protective alleles (OR < 1.0) are identified.

### 4.3 DiseaseActivityScorer

Calculates disease activity levels using simplified CRP/ESR-based scoring against 20 validated scoring systems.

**20 Scoring Systems:**

| Score | Disease | Range | Remission | Low | Moderate | High | Reference |
|---|---|---|---|---|---|---|---|
| DAS28-CRP | RA | 0-10 | <2.6 | 3.2 | 5.1 | >5.1 | PMID:15593215 |
| DAS28-ESR | RA | 0-10 | <2.6 | 3.2 | 5.1 | >5.1 | PMID:15593215 |
| SLEDAI-2K | SLE | 0-105 | 0 | 4 | 8 | >=12 | PMID:12115176 |
| CDAI | RA | 0-76 | <2.8 | 10 | 22 | >22 | PMID:15641075 |
| BASDAI | AS | 0-10 | <2 | 3 | 4 | >4 | PMID:8003055 |
| SDAI | RA | 0-86 | <3.3 | 11 | 26 | >26 | PMID:14872836 |
| PASI | Psoriasis | 0-72 | <1 | 5 | 10 | >10 | PMID:15888150 |
| Mayo Score | IBD | 0-12 | <2 | 5 | 8 | >8 | PMID:3317057 |
| Harvey-Bradshaw | IBD | 0-30 | <4 | 7 | 16 | >16 | PMID:7014041 |
| ESSDAI | Sjogren's | 0-123 | <1 | 5 | 14 | >14 | PMID:20032223 |
| mRSS | SSc | 0-51 | <5 | 14 | 29 | >29 | PMID:8546527 |
| EDSS | MS | 0-10 | <1.5 | 3.5 | 6.0 | >6.0 | PMID:6685237 |
| QMGS | MG | 0-39 | <3 | 10 | 20 | >20 | PMID:10668691 |
| Marsh Score | Celiac | 0-4 | 0 | 1 | 2 | >=3 | PMID:1437871 |
| Burch-Wartofsky | Graves' | 0-140 | <10 | 25 | 45 | >45 | PMID:8432869 |
| ASDAS | AS | 0-6 | <1.3 | 2.1 | 3.5 | >3.5 | PMID:19139421 |
| MG-ADL | MG | 0-24 | <1 | 5 | 10 | >10 | PMID:10025780 |
| DAPSA | Psoriasis | 0-164 | <4 | 14 | 28 | >28 | PMID:22328740 |
| HbA1c-T1D | T1D | 4-14% | <6.5 | 7.0 | 8.5 | >8.5 | PMID:9742976 |
| TSH-Hashimoto | Hashimoto's | 0-100 | <2.5 | 5.0 | 10.0 | >10.0 | PMID:12487769 |

**Scoring logic:** Retrieves CRP and/or ESR from patient biomarkers, maps to the applicable scoring system for diagnosed conditions, and classifies activity level against thresholds. Returns level (REMISSION/LOW/MODERATE/HIGH/VERY_HIGH), components, and threshold context.

### 4.4 FlarePredictor

Predicts flare risk from biomarker patterns using 13 disease-specific configurations.

**Algorithm:**
1. Start with base risk score of 0.3
2. For each early warning biomarker present in patient data:
   - Elevated inflammatory markers (CRP, ESR, IL-6, calprotectin > 5) add +0.15
   - Low complement (C3, C4 < 80) adds +0.15
   - Low albumin (< 3.5) adds +0.10
   - Normal values become protective factors
3. Clamp score to [0.0, 1.0]
4. Classify: IMMINENT (>=0.8), HIGH (>=0.6), MODERATE (>=0.4), LOW (<0.4)
5. Generate recommended monitoring actions

**13 Disease Patterns:** RA, SLE, IBD, AS, Psoriasis, Sjogren's, SSc, MS, T1D, MG, Celiac, Graves', Hashimoto's. Each with disease-specific early warning biomarkers and protective signals.

### 4.5 BiologicTherapyAdvisor

Matches 22 biologic therapies to patient diagnoses with PGx filtering.

**22 Therapies by Class:**

| Class | Therapies | Key PGx Considerations |
|---|---|---|
| TNF inhibitors (5) | Adalimumab, Etanercept, Infliximab, Golimumab, Certolizumab | FCGR3A V158F, HLA-DRB1*03:01 (anti-drug antibodies), TNFA -308 |
| IL-6R inhibitors (2) | Tocilizumab, Sarilumab | IL6R Asp358Ala, CRP masked by IL-6R blockade |
| Anti-CD20 (2) | Rituximab, Ocrelizumab | FCGR3A V158F, TNFSF13B levels, hepatitis B risk |
| IL-17A inhibitors (2) | Secukinumab, Ixekizumab | HLA-C*06:02 response prediction, IBD worsening risk |
| IL-12/23 inhibitors (2) | Ustekinumab, Risankizumab | IL12B/IL23R variants |
| IL-23 p19 inhibitor (1) | Guselkumab | IL23R rs11209026, HLA-C*06:02 |
| BLyS inhibitor (1) | Belimumab | TNFSF13B genotype |
| T-cell modulator (1) | Abatacept | CTLA4 +49 A/G, shared epitope response |
| JAK inhibitors (3) | Tofacitinib, Baricitinib, Upadacitinib | CYP3A4/CYP2C19 metabolism, VTE risk |
| Integrin inhibitors (2) | Vedolizumab, Natalizumab | JC virus antibody index (PML risk), anti-drug antibodies |
| TYK2 inhibitor (1) | Deucravacitinib | TYK2 P1104A variant, minimal CYP interaction |

**Selection logic:** Filters therapies by `indicated_diseases` matching patient diagnoses, returns all matching therapies with PGx considerations, contraindications, and monitoring requirements.

### 4.6 DiagnosticOdysseyAnalyzer

Multi-function diagnostic reasoning engine implementing:

**Classification Criteria Evaluation:**
- 10 ACR/EULAR criteria sets (RA, SLE, AS, SSc, Sjogren's, MS, MG, Celiac, IBD, Psoriasis)
- Point-based scoring against disease-specific thresholds
- Returns met/unmet criteria, total points, and threshold comparison

**Overlap Syndrome Detection:**
- 9 overlap syndromes (MCTD, SLE-RA overlap, Sjogren's-SLE, POTS/hEDS/MCAS, SSc-myositis, T1D-celiac, thyroid-T1D, RA-Sjogren's, lupus-APS)
- Detects based on antibody profile and diagnosed conditions
- Reports matched markers, involved diseases, and confidence level

**Differential Diagnosis:**
- Scores diseases from positive antibodies (weighted by specificity) and HLA alleles (log2-scaled odds ratio)
- Returns ranked differential with evidence for each disease

**Diagnostic Odyssey Analysis:**
- Reconstructs timeline from event records
- Calculates diagnostic delay (days/months/years from first symptom to diagnosis)
- Counts specialists seen, misdiagnoses, and key diagnostic tests
- Identifies turning points and missed opportunities

---

## 5. Multi-Collection RAG Engine

### 5.1 Search Flow

```
Query Text
    │
    ▼
BGE-small-en-v1.5 Embedding (384-dim, asymmetric with query prefix)
    │
    ▼
ThreadPoolExecutor: Parallel search across 14 collections (max_workers=6)
    │
    ▼
Weighted merge + deduplication (content hash + ID dedup)
    │
    ▼
Knowledge base augmentation (HLA, autoantibody, therapy, flare pattern context)
    │
    ▼
Disease area detection (keyword-based routing)
    │
    ▼
Conversation history injection (3-turn memory)
    │
    ▼
Claude Sonnet 4 prompt with patient context + evidence block
    │
    ▼
Grounded response with citations ([AutoAb:name], [HLA:allele], [Therapy:drug], [Literature:PMID])
```

### 5.2 Collection Weights

| Collection | Weight | Rationale |
|---|---|---|
| autoimmune_clinical_documents | 0.18 | Primary patient records |
| autoimmune_patient_labs | 0.14 | Lab results with flags |
| autoimmune_autoantibody_panels | 0.12 | Autoantibody reference |
| autoimmune_hla_associations | 0.08 | HLA-disease mapping |
| autoimmune_disease_criteria | 0.08 | Classification criteria |
| autoimmune_disease_activity | 0.07 | Activity scoring |
| autoimmune_flare_patterns | 0.06 | Flare patterns |
| autoimmune_biologic_therapies | 0.06 | Therapy database |
| autoimmune_clinical_trials | 0.05 | Clinical trials |
| autoimmune_literature | 0.05 | Published research |
| autoimmune_pgx_rules | 0.04 | PGx rules |
| autoimmune_patient_timelines | 0.03 | Diagnostic timelines |
| autoimmune_cross_disease | 0.02 | Overlap syndromes |
| genomic_evidence | 0.02 | Shared genomic context |
| **Total** | **1.00** | |

### 5.3 Citation Scoring

| Level | Threshold | Display |
|---|---|---|
| High confidence | >= 0.80 | Full citation with source link |
| Medium confidence | >= 0.60 | Citation with caveat |
| Below threshold | < 0.40 | Filtered out |

---

## 6. Export Pipeline

### 6.1 FHIR R4 Bundle

Produces a FHIR R4 Bundle containing:

- **Patient** resource with identifier
- **DiagnosticReport** resource (main report)
- **Observation** resources for disease activity scores
- **Observation** resources for flare risk predictions
- Conclusion field with critical alerts

### 6.2 PDF Export

Uses reportlab for clinical-grade PDF reports with NVIDIA green branding. Includes critical alerts, disease activity score tables, biologic therapy recommendations with PGx, and clinical query responses.

### 6.3 Markdown

Structured text format with tables for disease activity, flare predictions, HLA associations, and biologic recommendations.

---

## 7. Infrastructure

### 7.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4 |
| Embeddings | BGE-small-en-v1.5 (BAAI) -- 384-dim |
| LLM | Claude Sonnet 4 (Anthropic API) |
| Web UI | Streamlit (10 tabs) |
| REST API | FastAPI + Uvicorn |
| Configuration | Pydantic BaseSettings (AUTO_ prefix) |
| Testing | pytest |
| PDF Processing | PyPDF2 (ingestion), reportlab (export) |
| Export | FHIR R4, PDF, Markdown |
| Containerization | Docker + Docker Compose |
| Monitoring | Prometheus metrics endpoint |

### 7.2 Service Ports

| Service | Port |
|---|---|
| Streamlit UI | 8531 |
| FastAPI REST API | 8532 |
| Milvus (shared) | 19530 |

### 7.3 Dependencies on HCLS AI Factory

| Dependency | Type |
|---|---|
| Milvus 2.4 | Shared vector database (port 19530) |
| `genomic_evidence` collection | Read-only shared collection from Stage 2 RAG pipeline |
| BGE-small-en-v1.5 | Shared embedding model |
| Claude API key | Shared Anthropic API key |

---

## 8. Demo Patients

### 9 Demo Patients

| # | Name | Age/Sex | Primary Disease | Key Features |
|---|---|---|---|---|
| 1 | Sarah Mitchell | 34F | SLE | ANA 1:640, anti-dsDNA+, anti-Smith+, lupus nephritis, 27+ PDFs |
| 2 | Maya Rodriguez | 28F | POTS/hEDS/MCAS | Dysautonomia diagnostic odyssey |
| 3 | Linda Chen | 45F | Sjogren's | anti-SSA/Ro+, ESSDAI scoring |
| 4 | David Park | 45M | AS | HLA-B*27+, BASDAI scoring |
| 5 | Rachel Thompson | 38F | MCTD | Anti-RNP+, overlap syndrome |
| 6 | Emma Williams | 34F | MS (RRMS) | EDSS scoring, relapse monitoring |
| 7 | James Cooper | 19M | T1D + Celiac | Overlap syndrome, HLA-DQ2/DQ8+, GAD65+ |
| 8 | Karen Foster | 48F | SSc (dcSSc) | anti-Scl-70+, mRSS scoring |
| 9 | Michael Torres | 41M | Graves' Disease | TSI+, Burch-Wartofsky scoring |

---

## 9. File Structure (Actual)

```
precision_autoimmune_agent/
├── src/                            # 9 core modules (~4,500 lines)
│   ├── models.py                   # Pydantic models (AutoimmunePatientProfile, etc.)
│   ├── collections.py              # 14 Milvus collection schemas + manager
│   ├── rag_engine.py               # AutoimmuneRAGEngine (parallel search + Claude)
│   ├── agent.py                    # AutoimmuneAgent orchestrator
│   ├── knowledge.py                # Knowledge v2.0.0 (HLA, antibodies, biologics, flare)
│   ├── diagnostic_engine.py        # Classification criteria, overlap, differential
│   ├── document_processor.py       # PDF ingestion and chunking
│   ├── timeline_builder.py         # Diagnostic odyssey timeline
│   └── export.py                   # FHIR R4 + PDF + Markdown
├── app/
│   └── autoimmune_ui.py            # Streamlit (10 tabs, ~1,100 lines)
├── api/
│   └── main.py                     # FastAPI REST server (14 endpoints, ~580 lines)
├── config/
│   └── settings.py                 # AutoimmuneSettings (AUTO_ prefix)
├── data/                           # Reference data + cache + events
├── demo_data/                      # 9 patient directories with clinical PDFs
├── scripts/                        # setup_collections, generate_demo_patients, pdf_engine
├── tests/                          # 8 test files (455 tests)
├── docker-compose.yml
├── Dockerfile
└── requirements.txt
```

**44 Python files | ~20,000 lines | Apache 2.0**

---

## 10. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Autoimmune Agent is a clinical decision support research tool for autoimmune disease assessment. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
