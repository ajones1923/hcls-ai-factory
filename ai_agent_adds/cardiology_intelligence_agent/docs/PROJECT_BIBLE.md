# Cardiology Intelligence Agent -- Project Bible

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0
**Repository:** hcls-ai-factory/ai_agent_adds/cardiology_intelligence_agent

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Vision & Mission](#2-vision--mission)
3. [System Overview](#3-system-overview)
4. [File Inventory](#4-file-inventory)
5. [Collections Catalog](#5-collections-catalog)
6. [Knowledge Graph](#6-knowledge-graph)
7. [Clinical Workflows](#7-clinical-workflows)
8. [Risk Calculators](#8-risk-calculators)
9. [GDMT Optimizer](#9-gdmt-optimizer)
10. [Cross-Modal Triggers](#10-cross-modal-triggers)
11. [RAG Engine](#11-rag-engine)
12. [Query Expansion](#12-query-expansion)
13. [Data Models](#13-data-models)
14. [Export System](#14-export-system)
15. [API Reference](#15-api-reference)
16. [UI Guide](#16-ui-guide)
17. [Ingest Parsers](#17-ingest-parsers)
18. [Metrics & Monitoring](#18-metrics--monitoring)
19. [Scheduler](#19-scheduler)
20. [Configuration](#20-configuration)
21. [Docker Deployment](#21-docker-deployment)
22. [Port Map](#22-port-map)
23. [Tech Stack](#23-tech-stack)
24. [Future Roadmap](#24-future-roadmap)

---

## 1. Executive Summary

The **Cardiology Intelligence Agent** is a domain-specialized retrieval-augmented generation (RAG) system that synthesizes cardiac imaging, electrophysiology, hemodynamics, heart failure management, valvular disease, preventive cardiology, interventional data, and cardio-oncology surveillance into guideline-aligned clinical recommendations using ACC/AHA/ESC evidence.

The system searches 12 cardiology-specific Milvus vector collections plus 1 shared genomic_evidence collection (3.5M variant vectors), implements 6 validated cardiovascular risk calculators (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II), optimizes guideline-directed medical therapy (GDMT) for heart failure (7 therapy classes including finerenone, omecamtiv, and sotagliflozin), and provides 11 clinical workflows covering the highest-impact cardiovascular use cases.

The Cardiology Intelligence Agent is built as part of the HCLS AI Factory Precision Intelligence Network -- one of three GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence Network, Therapeutic Discovery Engine) that compose the end-to-end precision medicine platform on NVIDIA DGX Spark. The Precision Intelligence Network comprises 11 specialized agents:

| # | Agent | Domain |
|---|-------|--------|
| 1 | Precision Biomarker | Biomarker interpretation and trends |
| 2 | Precision Oncology | Cancer genomics and treatment |
| 3 | CAR-T Intelligence | CAR-T cell therapy |
| 4 | Medical Imaging | Multi-modal imaging analysis |
| 5 | Precision Autoimmune | Autoimmune disease analysis |
| 6 | **Cardiology Intelligence** | **Cardiovascular decision support** |
| 7 | Clinical Trial Intelligence | Trial matching and eligibility |
| 8 | Neurology Intelligence | Neurological decision support |
| 9 | Rare Disease Intelligence | Rare disease diagnosis |
| 10 | Pharmacogenomics (PGx) | Drug-gene interactions |
| 11 | Pediatric Oncology | Pediatric cancer intelligence |

### Codebase at a Glance

| Metric | Value |
|--------|-------|
| Total Python LOC | 28,189 |
| Total files | 43 (37 Python, 1 TOML, 1 YAML, 2 Markdown, 1 requirements.txt, 1 .env template) |
| Milvus collections | 13 (12 cardiology-specific + 1 shared genomic_evidence) |
| Knowledge graph conditions | 45 cardiovascular conditions |
| Knowledge graph biomarkers | 29 cardiac biomarkers |
| Knowledge graph drug classes | 32 drug classes |
| Knowledge graph genes | 56 cardiovascular genes |
| Knowledge graph imaging modalities | 27 imaging protocols (12 echo, 7 CMR, 5 nuclear, 3 CT) |
| Guideline documents | 20 guideline documents |
| Guideline recommendations | 63 structured ACC/AHA/ESC recommendations |
| Entity aliases | 167 abbreviation mappings |
| Risk calculators | 6 validated scoring systems |
| Clinical workflows | 11 |
| Cross-modal imaging triggers | 18 genomic trigger patterns |
| Ingest parsers | 7 (PubMed, trials, imaging, ECG, guideline, device, hemodynamics) |
| Enums | 16 |
| Pydantic models | 13 |
| Dataclasses | 1 (SearchPlan) |
| LLM | Claude Sonnet 4.6 (claude-sonnet-4-6) |
| Embedding model | BGE-small-en-v1.5 (384-dim) |

---

## 2. Vision & Mission

### The Problem

Cardiovascular disease (CVD) remains the leading cause of death globally, claiming approximately 17.9 million lives annually -- 32% of all deaths worldwide. Despite exponential growth in cardiovascular AI research (over 4,500 publications in 2025 alone), critical barriers persist in clinical translation: fragmented evidence across modalities, siloed data systems, lack of integrated genomic-imaging correlation, and prohibitive infrastructure costs that limit advanced cardiac AI to elite academic centers.

Current tools fall short:

| Approach | Limitation |
|---|---|
| PubMed search | Keyword-based; no cross-modal integration |
| UpToDate / DynaMed | Expert-curated but static; no patient-specific reasoning |
| Commercial CVIS | Vendor-locked; $500K-$2M+ implementation |
| EHR-integrated CDS | Rule-based; cannot synthesize unstructured evidence |
| General AI assistants | No citation provenance; hallucination risk |
| Imaging-only AI | Single-modality; no genomic integration; cloud-dependent |

### Mission Statement

Democratize cardiovascular clinical decision support by providing unified, evidence-grounded intelligence across the entire cardiology landscape -- from imaging interpretation through risk stratification through GDMT optimization -- on hardware that any clinic or research lab can afford.

### Design Principles

1. **Guideline-Aligned**: Every recommendation traces to ACC/AHA/ESC guidelines with recommendation class and evidence level.
2. **Multi-Modal**: Integrates imaging, electrophysiology, hemodynamics, genomics, biomarkers, and clinical data.
3. **Quantitative**: Risk calculators produce specific scores, not vague guidance. GDMT optimizer recommends specific doses.
4. **Cross-Modal**: Imaging findings automatically trigger appropriate genomic workup panels.
5. **Transparent**: All responses cite source collections, similarity scores, and guideline references.

---

## 3. System Overview

### High-Level Architecture

```
+------------------------------------------------------------------+
|                    STREAMLIT UI (8536)                             |
|  [Dashboard] [CAD] [Heart Failure] [Valvular] [Arrhythmia]       |
|  [Cardiac MRI] [Stress Test] [Prevention] [Cardio-Onc] [Reports] |
|  Sidebar: Workflow filter, Collection stats, Risk calc, Export    |
+------------------------------------------------------------------+
         |                                          |
         v                                          v
+------------------+                    +---------------------+
|  FastAPI Server  |                    |  Clinical Engines   |
|  :8126           |                    |  - Risk Calculators |
|  POST /v1/cardio |                    |  - GDMT Optimizer   |
|  GET  /health    |                    |  - Cross-Modal      |
+------------------+                    |  - Workflows        |
         |                              +---------------------+
         v                                          |
+------------------+                                |
|  RAG Engine      |<-------------------------------+
|  - Query Expand  |
|  - Multi-Collect |
|  - Citation Score|
|  - LLM Synthesis |
+------------------+
         |
         v
+--------------------------------------------------+
|              Milvus Vector DB (:19530)            |
|  12 cardio collections + genomic_evidence         |
|  384-dim BGE-small-en-v1.5 / IVF_FLAT / COSINE  |
+--------------------------------------------------+
|  etcd (:2379)  |  MinIO (:9000)                  |
+--------------------------------------------------+
```

### Component Interaction Flow

1. **Query arrives** via Streamlit UI or FastAPI endpoint
2. **Query expansion** decomposes into sub-questions, identifies conditions/drugs/modalities
3. **Workflow routing** selects applicable clinical workflows (CAD, HF, valvular, etc.)
4. **Multi-collection search** queries relevant Milvus collections with weighted scoring
5. **Risk calculators** compute applicable scores (ASCVD, CHA2DS2-VASc, etc.)
6. **GDMT optimizer** evaluates heart failure medication status and recommends titration
7. **Cross-modal engine** identifies genomic workup triggers from imaging findings
8. **LLM synthesis** (Claude Sonnet 4.6) generates evidence-grounded clinical response
9. **Citation scoring** ranks and formats source citations
10. **Response assembly** packages answer, citations, risk scores, workflow results, triggers

---

## 4. File Inventory

### Source Files (`src/`)

| File | Lines | Description |
|------|-------|-------------|
| `agent.py` | 1,658 | Main agent orchestrator; coordinates all engines |
| `clinical_workflows.py` | 2,445 | 11 clinical workflow implementations |
| `collections.py` | 1,226 | Milvus collection CRUD and schema definitions |
| `cross_modal.py` | 1,734 | Imaging-genomics cross-modal trigger engine |
| `export.py` | 1,379 | PDF/CSV/FHIR export system |
| `gdmt_optimizer.py` | 2,457 | 4-pillar GDMT optimization engine |
| `knowledge.py` | 1,431 | Knowledge graph (conditions, biomarkers, drugs, genes, imaging, guidelines) |
| `metrics.py` | 537 | Prometheus metrics instrumentation |
| `models.py` | 717 | 16 enums, 13 Pydantic models, 1 dataclass |
| `query_expansion.py` | 2,025 | Query decomposition, synonym expansion, entity extraction |
| `rag_engine.py` | 1,589 | Multi-collection RAG retrieval and synthesis |
| `risk_calculators.py` | 2,397 | 6 validated risk calculator implementations |
| `scheduler.py` | 612 | APScheduler-based periodic ingest scheduler |

### Ingest Parsers (`src/ingest/`)

| File | Lines | Description |
|------|-------|-------------|
| `base.py` | 224 | Abstract base parser class |
| `pubmed_parser.py` | 477 | PubMed cardiovascular literature via NCBI E-utilities |
| `clinical_trials_parser.py` | 713 | ClinicalTrials.gov cardiovascular trials |
| `imaging_parser.py` | 519 | Cardiac imaging reports and protocols |
| `ecg_parser.py` | 424 | ECG interpretation data |
| `guideline_parser.py` | 555 | ACC/AHA/ESC guideline documents |
| `device_parser.py` | 428 | FDA-cleared cardiovascular devices |
| `hemodynamics_parser.py` | 402 | Catheterization and hemodynamic data |

### API Layer (`api/`)

| File | Lines | Description |
|------|-------|-------------|
| `main.py` | 433 | FastAPI application with CORS, error handling, lifespan |
| `routes/cardio_clinical.py` | 1,020 | Clinical query, risk calc, GDMT, workflow endpoints |
| `routes/reports.py` | 306 | Report generation and export endpoints |
| `routes/events.py` | 147 | SSE event streaming |

### Application (`app/`)

| File | Lines | Description |
|------|-------|-------------|
| `cardio_ui.py` | 1,182 | Streamlit 10-tab cardiovascular UI |

### Configuration (`config/`)

| File | Lines | Description |
|------|-------|-------------|
| `settings.py` | 181 | Pydantic BaseSettings with validation |

### Scripts (`scripts/`)

| File | Lines | Description |
|------|-------|-------------|
| `setup_collections.py` | 231 | Create/drop Milvus collections and seed |
| `seed_knowledge.py` | 412 | Seed knowledge graph into collections |
| `run_ingest.py` | 275 | Execute all ingest parsers |

### Tests (`tests/`)

| File | Lines | Description |
|------|-------|-------------|
| `conftest.py` | 15 | Pytest configuration and fixtures |

### Infrastructure

| File | Description |
|------|-------------|
| `docker-compose.yml` | Full stack: etcd, MinIO, Milvus, Streamlit, API, setup |
| `requirements.txt` | Python dependencies |
| `.streamlit/config.toml` | Streamlit theme and server configuration |
| `README.md` | Project overview and quickstart |

---

## 5. Collections Catalog

### 12 Cardiology-Specific Collections + 1 Shared

All collections use BGE-small-en-v1.5 embeddings (384-dim), IVF_FLAT index, COSINE similarity.

| # | Collection Name | Description | Search Weight |
|---|----------------|-------------|---------------|
| 1 | `cardio_literature` | Published cardiovascular research, reviews, meta-analyses | 0.10 |
| 2 | `cardio_trials` | Cardiovascular clinical trials and landmark trial results | 0.08 |
| 3 | `cardio_imaging` | Cardiac imaging protocols, findings, measurements (echo, CT, MRI, nuclear) | 0.10 |
| 4 | `cardio_electrophysiology` | ECG interpretation, arrhythmia classification, EP data | 0.08 |
| 5 | `cardio_heart_failure` | HF classification, GDMT protocols, management algorithms | 0.10 |
| 6 | `cardio_valvular` | Valvular heart disease assessment and intervention criteria | 0.08 |
| 7 | `cardio_prevention` | Preventive cardiology: risk stratification, lipid management | 0.10 |
| 8 | `cardio_interventional` | Interventional procedures, techniques, outcomes | 0.07 |
| 9 | `cardio_oncology` | Cardio-oncology surveillance, cardiotoxicity detection | 0.06 |
| 10 | `cardio_devices` | FDA-cleared cardiovascular AI devices, implantable devices | 0.04 |
| 11 | `cardio_guidelines` | ACC/AHA/ESC/HRS clinical practice guidelines | 0.10 |
| 12 | `cardio_hemodynamics` | Catheterization data, pressure tracings, derived calculations | 0.06 |
| 13 | `genomic_evidence` | Shared genomic evidence (read-only, 3.5M variants) | 0.03 |

**Total weights:** 1.00

### Collection Schema

Each collection follows a consistent schema:

| Field | Type | Description |
|-------|------|-------------|
| `id` | INT64 (PK, auto) | Unique identifier |
| `text` | VARCHAR(65535) | Source text content |
| `embedding` | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| `source` | VARCHAR(512) | Source identifier (DOI, URL, guideline ID) |
| `title` | VARCHAR(1024) | Document title |
| `metadata` | JSON | Additional structured metadata |
| `created_at` | VARCHAR(32) | Ingestion timestamp |

---

## 6. Knowledge Graph

The knowledge graph (`src/knowledge.py`, 1,431 lines) provides structured domain data across 7 categories:

### 6.1 Cardiovascular Conditions (45 entries)

Each condition includes: aliases, ICD-10 code, prevalence, diagnostic criteria, risk factors, imaging modalities, treatment options, associated genes, guideline references, and cross-modal trigger flag.

**Conditions by category:**

| Category | Conditions |
|----------|-----------|
| Cardiomyopathies (7) | HCM, DCM, ARVC, Takotsubo, Restrictive CM, Cardiac Amyloidosis, Cardiac Sarcoidosis |
| Coronary Disease (3) | CAD, ACS, Stable Angina |
| Heart Failure (2) | HFrEF, HFpEF |
| Arrhythmias (7) | AF, Atrial Flutter, VT, LQTS, Brugada, CPVT, WPW |
| Valvular Disease (5) | Aortic Stenosis, Aortic Regurgitation, Mitral Regurgitation, Mitral Stenosis, Tricuspid Regurgitation |
| Vascular (3) | Pulmonary Hypertension, Aortic Dissection, Thoracic Aortic Aneurysm |
| Inflammatory (3) | Myocarditis, Pericarditis, Infective Endocarditis |
| Metabolic/Other (2) | Familial Hypercholesterolemia, Cardiac Tamponade |

### 6.2 Cardiac Biomarkers (29 entries)

Each biomarker includes: full name, LOINC code, reference ranges (age/sex stratified where applicable), clinical use, kinetics, confounders, and guideline thresholds.

**Biomarkers:** hs-cTnI, hs-cTnT, NT-proBNP, BNP, CK-MB, D-dimer, hsCRP, LDL-C, HDL-C, Total Cholesterol, Triglycerides, Lp(a), ApoB, HbA1c, Creatinine/eGFR, Potassium, Magnesium, Ferritin/Iron/TSAT, Lactate, Procalcitonin, IL-6

### 6.3 Cardiac Drug Classes (32 entries)

Each drug class includes: specific drugs, mechanism, indications, contraindications, key trials, target doses, and monitoring parameters.

**Drug classes:** Beta-blockers, ACE inhibitors, ARBs, ARNI, SGLT2 inhibitors, MRA, Loop diuretics, Thiazide diuretics, DHP CCBs, Non-DHP CCBs, Statins, PCSK9 inhibitors, Ezetimibe, Antiplatelets, DOACs, Warfarin, Heparin/LMWH, Class I antiarrhythmics, Class III antiarrhythmics, Digoxin, Nitrates, Hydralazine/ISDN, Ivabradine, Mavacamten, GLP-1 receptor agonists, Inotropes

### 6.4 Cardiovascular Genes (56 entries)

Each gene includes: full name, chromosome locus, function, associated conditions, inheritance pattern, key variants, clinical significance, and genetic testing indication.

**Gene categories:**
- **Cardiomyopathy (11):** MYH7, MYBPC3, TNNT2, TNNI3, TPM1, ACTC1, MYL2, MYL3, TTN, LMNA, RBM20
- **Desmosomal/ARVC (6):** PKP2, DSP, DSG2, DSC2, JUP, TMEM43
- **Channelopathy (7):** SCN5A, KCNQ1, KCNH2, KCNJ2, RYR2, CASQ2, CALM1, ANK2
- **Aortopathy (7):** FBN1, TGFBR1, TGFBR2, SMAD3, ACTA2, MYH11, COL3A1
- **Lipid/Metabolic (5):** LDLR, PCSK9, APOB, APOE, LPA
- **Congenital/Other (4):** NKX2-5, GATA4, TBX5, GLA

### 6.5 Imaging Modalities (27 protocols: 12 echo, 7 CMR, 5 nuclear, 3 CT)

Each modality includes: full name, protocols, key measurements with reference values, indications, contraindications, and guideline society.

**Modalities:** TTE, TEE, Stress Echo, Cardiac CT Calcium Score, Coronary CTA, Cardiac MRI, Cardiac PET, SPECT MPI, MUGA, Right Heart Catheterization, Coronary Angiography, 12-Lead ECG, Holter Monitor, Event Monitor/ILR, Device Interrogation

### 6.6 Guideline Recommendations (63 entries across 20 guideline documents)

Structured guideline entries with: society, year, recommendation text, class (I/IIa/IIb/III), evidence level (A/B/C), condition, and source citation.

**Guideline coverage:**
- ACC/AHA 2022 Heart Failure (8 recommendations)
- ACC/AHA 2020 Valvular Heart Disease (6 recommendations)
- ACC/AHA 2018 Cholesterol (6 recommendations)
- ACC/AHA 2019 Primary Prevention (3 recommendations)
- ACC/AHA/HRS 2023 Atrial Fibrillation (5 recommendations)
- ACC/AHA 2024 HCM (6 recommendations)
- ESC 2023 ACS (4 recommendations)
- ESC 2022 Cardio-Oncology (4 recommendations)
- Additional key recommendations (10)

### 6.7 Entity Aliases (167 entries)

Bidirectional abbreviation-to-full-name mapping for cardiovascular terminology, enabling query expansion and entity recognition via 18 synonym maps. Covers abbreviations for conditions, procedures, medications, imaging modalities, biomarkers, professional societies, and clinical measurements.

### 6.8 Condition-Drug Map (11 entries)

Direct mapping from conditions to applicable drug classes for treatment recommendation generation.

### 6.9 Imaging Trigger Map (18 entries)

Maps imaging findings to recommended genetic testing panels. Expanded from 6 to 18 genomic trigger patterns, including:
- Unexplained LVH >=15mm --> HCM/Fabry/Danon gene panel
- Non-ischemic LGE with DCM --> DCM gene panel
- RV dilation with arrhythmia --> ARVC desmosomal panel
- Diffuse LGE + low voltage --> TTR (amyloidosis)
- Aortic root dilation --> Aortopathy panel
- Premature coronary calcification --> FH panel
- Plus 12 additional trigger patterns for channelopathies, congenital disease, and infiltrative cardiomyopathies

---

## 7. Clinical Workflows

Eleven clinical workflows (`src/clinical_workflows.py`) cover the highest-impact cardiovascular use cases:

| # | Workflow | Key Actions |
|---|----------|-------------|
| 1 | **CAD Assessment** | Calcium scoring interpretation, CAD-RADS classification, plaque characterization, FFR/iFR correlation, revascularization criteria |
| 2 | **Heart Failure Classification & GDMT** | NYHA/ACC staging, EF categorization (HFrEF/HFmrEF/HFpEF/HFimpEF), 4-pillar GDMT optimization, device eligibility (ICD/CRT) |
| 3 | **Valvular Heart Disease Quantification** | Severity grading per ASE criteria, intervention threshold assessment, TAVR vs SAVR decision support |
| 4 | **Arrhythmia Detection & Management** | ECG interpretation (15 criteria, 10 arrhythmia protocols), CHA2DS2-VASc scoring, anticoagulation recommendation, rhythm vs rate control strategy |
| 5 | **Cardiac MRI Tissue Characterization** | LGE pattern analysis (ischemic vs non-ischemic), T1/T2 mapping interpretation, ECV quantification, parametric analysis |
| 6 | **Stress Test Interpretation** | Duke Treadmill Score, SPECT/PET MPI defect analysis, stress echo wall motion, ischemic burden quantification |
| 7 | **Preventive Risk Stratification** | ASCVD risk calculation, statin eligibility algorithm, CAC reclassification, risk enhancer assessment |
| 8 | **Cardio-Oncology Surveillance** | Cardiotoxicity risk stratification, CTRCD detection (LVEF decline, GLS decline), biomarker monitoring, cardioprotection |
| 9 | **Acute Decompensated Heart Failure** | Hemodynamic profiling (14 parameters, 6 cath lab protocols), diuretic strategy, inotrope selection, mechanical support evaluation |
| 10 | **Post-Myocardial Infarction** | Post-MI risk stratification, secondary prevention, cardiac rehabilitation, follow-up imaging timing |
| 11 | **Myocarditis/Pericarditis** | Lake Louise criteria, inflammatory biomarker tracking, activity restriction guidance, colchicine/NSAID protocols |

Each workflow:
- Receives a `CardioQuery` with patient context
- Queries relevant Milvus collections with workflow-specific weights
- Executes applicable risk calculators
- Generates `WorkflowResult` with findings, recommendations, guideline references, severity, and cross-modal triggers

---

## 8. Risk Calculators

Six validated cardiovascular risk calculators (`src/risk_calculators.py`, 2,397 lines):

### 8.1 ASCVD Pooled Cohort Equations (Goff 2013)

- **Purpose:** 10-year atherosclerotic cardiovascular disease risk
- **Inputs:** Age, sex, race, total cholesterol, HDL, systolic BP, BP treatment status, diabetes, smoking
- **Coefficients:** 4 sex/race cohorts (white female, African American female, white male, African American male)
- **Output:** 10-year risk percentage + risk category (low <5%, borderline 5-7.5%, intermediate 7.5-20%, high >=20%)
- **Source:** Goff DC Jr, et al. J Am Coll Cardiol. 2014;63(25 Pt B):2935-2959

### 8.2 HEART Score

- **Purpose:** Chest pain risk stratification in the emergency department
- **Inputs:** History, ECG, Age, Risk factors, Troponin (0-2 points each)
- **Output:** Score 0-10; Low (0-3, 1.7% MACE), Moderate (4-6, 16.6% MACE), High (7-10, 50.1% MACE)
- **Source:** Six AJ, et al. Neth Heart J. 2008;16(6):191-196

### 8.3 CHA2DS2-VASc

- **Purpose:** Atrial fibrillation stroke risk assessment
- **Inputs:** CHF, Hypertension, Age (>=75 = 2 pts, 65-74 = 1 pt), Diabetes, Stroke/TIA history (2 pts), Vascular disease, Sex (female = 1 pt)
- **Output:** Score 0-9 with annual stroke rate lookup
- **Source:** Lip GYH, et al. Chest 2010;137(2):263-272

### 8.4 HAS-BLED

- **Purpose:** Bleeding risk assessment on anticoagulation
- **Inputs:** Hypertension, Abnormal renal/liver function, Stroke, Bleeding history, Labile INR, Elderly (>65), Drugs/alcohol
- **Output:** Score 0-9; Low (0-1), Moderate (2), High (>=3)
- **Source:** Pisters R, et al. Chest 2010;138(5):1093-1100

### 8.5 MAGGIC Heart Failure Risk Score

- **Purpose:** Heart failure 1-year and 3-year mortality prediction
- **Inputs:** Age, sex, LVEF, NYHA class, systolic BP, BMI, creatinine, diabetes, beta-blocker use, ACEi/ARB use, COPD, HF duration, smoking
- **Output:** Integer score with mortality percentage lookup (51 tiers for 1-year and 3-year)
- **Source:** Pocock SJ, et al. Eur Heart J. 2013;34(19):1404-1413

### 8.6 EuroSCORE II

- **Purpose:** Cardiac surgical mortality risk prediction
- **Inputs:** Patient factors (age, sex, renal function, mobility, prior surgery, lung disease, endocarditis, critical state, diabetes), cardiac factors (NYHA, angina, LVEF, recent MI, pulmonary hypertension), operation factors (urgency, procedure weight, thoracic aorta)
- **Coefficients:** Logistic regression with 28 beta coefficients
- **Output:** Predicted operative mortality percentage
- **Source:** Nashef SA, et al. Eur J Cardiothorac Surg. 2012;41(4):734-745

### Risk Calculator Engine

The `RiskCalculatorEngine` class provides:
- `calculate(score_type, input, extra)` -- compute a single score
- `calculate_all_applicable(input, extra)` -- compute all scores for which sufficient data exists
- Input validation with descriptive error messages via `RiskCalculatorError`
- Each calculator returns `RiskScoreResult` with score value, risk category, interpretation, recommendations, and guideline reference

---

## 9. GDMT Optimizer

The GDMT optimizer (`src/gdmt_optimizer.py`, 2,457 lines) implements evidence-based guideline-directed medical therapy optimization for heart failure with reduced ejection fraction (HFrEF).

### Four Pillars of GDMT

| Pillar | Drug Class | Target Dose | Key Trials |
|--------|-----------|-------------|------------|
| 1. Beta-Blocker | Carvedilol, Metoprolol Succinate, Bisoprolol | Carvedilol 25mg BID, Metoprolol 200mg QD, Bisoprolol 10mg QD | MERIT-HF, COPERNICUS, CIBIS-II |
| 2. ARNI/ACEi/ARB | Sacubitril/Valsartan (preferred), Enalapril, Valsartan | Sacubitril/Valsartan 97/103mg BID | PARADIGM-HF, PARAGON-HF |
| 3. MRA | Spironolactone, Eplerenone | Spironolactone 25-50mg QD, Eplerenone 50mg QD | RALES, EMPHASIS-HF |
| 4. SGLT2i | Dapagliflozin, Empagliflozin | 10mg QD (both) | DAPA-HF, EMPEROR-Reduced |

### Optimizer Features

- **EF-stratified recommendations:** Different GDMT algorithms for HFrEF, HFmrEF, HFpEF, and HFimpEF
- **Titration sequencing:** Recommends initiation order and uptitration steps
- **Contraindication checking:** Validates each pillar against patient-specific contraindications (K+, eGFR, BP, HR)
- **Lab monitoring guidance:** Generates lab check schedules for potassium, creatinine, BP after each titration
- **Device therapy eligibility:** Evaluates ICD and CRT criteria based on LVEF, QRS duration, LBBB, NYHA class
- **Additional therapies (7 total):** H-ISDN for African American patients, ivabradine for persistent tachycardia, IV iron for iron deficiency, finerenone (non-steroidal MRA), omecamtiv mecarbil (cardiac myosin activator), sotagliflozin (dual SGLT1/2 inhibitor), vericiguat (sGC stimulator)

---

## 10. Cross-Modal Triggers

The cross-modal engine (`src/cross_modal.py`, 1,734 lines) links cardiac imaging findings to genomic workup recommendations via the shared `genomic_evidence` collection.

### Trigger Patterns

| Imaging Finding | Triggered Gene Panel | Suspected Conditions |
|----------------|---------------------|---------------------|
| Unexplained LVH >=15mm | MYH7, MYBPC3, TNNT2, TNNI3, TPM1, GLA, LAMP2, PRKAG2 | HCM, Fabry Disease, Danon Disease |
| Non-ischemic LGE with DCM | TTN, LMNA, RBM20, MYH7, DSP, FLNC, BAG3 | DCM, Arrhythmogenic CM |
| RV dilation with arrhythmia | PKP2, DSP, DSG2, DSC2, JUP, TMEM43 | ARVC |
| Diffuse LGE + low voltage ECG | TTR | Cardiac Amyloidosis |
| Aortic root dilation (>=4.0cm) | FBN1, TGFBR1, TGFBR2, SMAD3, ACTA2, MYH11, COL3A1 | Marfan, Loeys-Dietz, Familial TAAD |
| Premature coronary calcification | LDLR, PCSK9, APOB | Familial Hypercholesterolemia |

### Cross-Modal Workflow

1. Imaging findings are extracted from the workflow result
2. Pattern matching identifies applicable trigger entries from `IMAGING_TRIGGER_MAP`
3. For each matched trigger, the genomic_evidence collection is queried for relevant variant data
4. Results are packaged as `CrossModalTrigger` objects with source, finding, gene panel, conditions, and clinical rationale
5. Triggers are included in the final `CardioResponse` for clinician review

---

## 11. RAG Engine

The RAG engine (`src/rag_engine.py`, 1,589 lines) implements multi-collection retrieval-augmented generation:

### Search Pipeline

1. **Embedding:** Input query is embedded using BGE-small-en-v1.5 (384-dim)
2. **Multi-collection search:** Queries all 13 collections in parallel with `top_k=5` per collection
3. **Score thresholding:** Filters results below `SCORE_THRESHOLD=0.4`
4. **Weighted ranking:** Applies collection-specific weights (see Collections Catalog)
5. **Deduplication:** Removes near-duplicate content across collections
6. **Citation scoring:** Classifies citations as high (>=0.75), medium (>=0.60), or standard confidence

### LLM Synthesis

- **Provider:** Anthropic (Claude Sonnet 4.6)
- **System prompt:** Cardiology domain expert with guideline citation requirements
- **Context window:** Retrieved passages assembled with source attribution
- **Output format:** Structured clinical answer with inline citations

### Conversation Memory

- Maintains up to 3 conversation turns for contextual follow-up questions
- Configurable via `MAX_CONVERSATION_CONTEXT` setting

---

## 12. Query Expansion

The query expansion module (`src/query_expansion.py`, 2,025 lines) enhances retrieval through:

### Expansion Strategies

1. **Synonym expansion:** Maps abbreviations to full names using 167 entity aliases across 18 synonym maps (e.g., "HCM" -> "hypertrophic cardiomyopathy")
2. **Condition identification:** Extracts cardiovascular conditions from the query using knowledge graph matching
3. **Drug identification:** Identifies mentioned drug classes and specific medications
4. **Imaging modality detection:** Recognizes imaging modality references for collection routing
5. **Sub-question decomposition:** Breaks complex queries into searchable sub-questions
6. **Search strategy selection:** Chooses between broad, targeted, comparative, and clinical strategies
7. **Topic identification:** Tags query with relevant topics for collection weight adjustment

### Output: SearchPlan

The query expansion module produces a `SearchPlan` dataclass containing:
- Original question
- Identified conditions, drugs, and imaging modalities
- Relevant workflow types
- Search strategy (broad/targeted/comparative/clinical)
- Decomposed sub-questions
- Identified topics

---

## 13. Data Models

### Enums (16)

| Enum | Values | Description |
|------|--------|-------------|
| `CardioWorkflowType` | 12 | Clinical workflow types (9 original + 3 new: acute_decompensated_hf, post_mi, myocarditis_pericarditis) |
| `RiskScoreType` | 6 | Risk calculator types |
| `SeverityLevel` | 5 | Clinical finding severity (critical/high/moderate/low/informational) |
| `ImagingModality` | 5 | Echo, CT, MRI, Nuclear, Angiography |
| `HeartFailureClass` | 4 | NYHA I-IV |
| `HeartFailureStage` | 4 | ACC/AHA Stage A-D |
| `EjectionFractionCategory` | 4 | HFrEF, HFmrEF, HFpEF, HFimpEF |
| `ValveSeverity` | 4 | Mild, Moderate, Severe, Critical |
| `CADRADSScore` | 7 | CAD-RADS 0-5 stenosis grading |
| `AnticoagulationRecommendation` | 4 | No anticoag to Contraindicated |
| `GDMTPillar` | 4 | Beta-blocker, ARNI/ACEi/ARB, MRA, SGLT2i |
| `GDMTStatus` | 6 | Not started through Intolerant |
| `CardiotoxicityRisk` | 4 | Low, Moderate, High, Very High |
| `LGEPattern` | 7 | Subendocardial through None |
| `GuidelineClass` | 5 | ACC/AHA recommendation classes I/IIa/IIb/III |
| `EvidenceLevel` | 5 | A, B-R, B-NR, C-LD, C-EO |

### Pydantic Models (13)

| Model | Fields | Purpose |
|-------|--------|---------|
| `CardioQuery` | 3 | Input query with optional workflow and patient context |
| `CardioSearchResult` | 4 | Single search result from any collection |
| `RiskScoreInput` | 30+ | Universal risk calculator input (demographics, vitals, labs, comorbidities, HF-specific, surgical) |
| `RiskScoreResult` | 6 | Calculator output (score, category, interpretation, recommendations) |
| `GDMTMedication` | 6 | Single GDMT pillar medication status |
| `GDMTRecommendation` | 5 | GDMT optimization output |
| `ValveAssessment` | 6 | Valve lesion assessment |
| `ECGInterpretation` | 6 | Structured ECG interpretation |
| `ImagingResult` | 5 | Cardiac imaging result with cross-modal triggers |
| `CardiotoxicityAssessment` | 8 | Cardio-oncology surveillance assessment |
| `CrossModalTrigger` | 5 | Imaging-to-genomics trigger |
| `WorkflowResult` | 7 | Single workflow execution output |
| `CardioResponse` | 6 | Top-level agent response |

### Dataclass (1)

| Class | Fields | Purpose |
|-------|--------|---------|
| `SearchPlan` | 8 | Pre-retrieval search planning output |

---

## 14. Export System

The export system (`src/export.py`, 1,379 lines) supports multiple output formats:

| Format | Description |
|--------|-------------|
| **PDF** | Structured clinical report with risk scores, GDMT status, recommendations, citations |
| **CSV** | Tabular export of risk scores, medications, findings |
| **JSON** | Full CardioResponse serialization |
| **FHIR** | HL7 FHIR R4-compatible DiagnosticReport and RiskAssessment resources |

---

## 15. API Reference

**Base URL:** `http://localhost:8126`

### Core Endpoints

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/health` | Health check with collection count and vector totals |
| `POST` | `/v1/cardio/query` | General cardiology query |
| `POST` | `/v1/cardio/risk/ascvd` | ASCVD risk calculation |
| `POST` | `/v1/cardio/risk/heart` | HEART score calculation |
| `POST` | `/v1/cardio/risk/cha2ds2-vasc` | CHA2DS2-VASc calculation |
| `POST` | `/v1/cardio/risk/has-bled` | HAS-BLED calculation |
| `POST` | `/v1/cardio/risk/maggic` | MAGGIC HF risk calculation |
| `POST` | `/v1/cardio/risk/euroscore-ii` | EuroSCORE II calculation |
| `POST` | `/v1/cardio/gdmt/optimize` | GDMT optimization |
| `POST` | `/v1/cardio/workflow/cad` | CAD assessment workflow |
| `POST` | `/v1/cardio/workflow/heart-failure` | Heart failure workflow |
| `POST` | `/v1/cardio/workflow/valvular` | Valvular disease workflow |
| `POST` | `/v1/cardio/workflow/arrhythmia` | Arrhythmia workflow |
| `POST` | `/v1/cardio/workflow/cardiac-mri` | Cardiac MRI workflow |
| `POST` | `/v1/cardio/workflow/stress-test` | Stress test workflow |
| `POST` | `/v1/cardio/workflow/prevention` | Prevention workflow |
| `POST` | `/v1/cardio/workflow/cardio-onc` | Cardio-oncology workflow |
| `POST` | `/v1/cardio/workflow/acute-dhf` | Acute decompensated HF workflow |
| `POST` | `/v1/cardio/workflow/post-mi` | Post-MI workflow |
| `POST` | `/v1/cardio/workflow/myocarditis-pericarditis` | Myocarditis/Pericarditis workflow |
| `GET` | `/v1/cardio/knowledge-version` | Knowledge base version info |
| `POST` | `/v1/cardio/report/pdf` | Generate PDF report |
| `POST` | `/v1/cardio/report/csv` | Generate CSV export |
| `POST` | `/v1/cardio/report/fhir` | Generate FHIR resources |
| `GET` | `/v1/cardio/collections` | Collection statistics |
| `GET` | `/metrics` | Prometheus metrics endpoint |

---

## 16. UI Guide

The Streamlit UI (`app/cardio_ui.py`, 1,182 lines) provides 10 tabs at port 8536:

| Tab | Description |
|-----|-------------|
| 1. **Dashboard** | System overview: collection stats, recent queries, health status |
| 2. **CAD Assessment** | Coronary artery disease evaluation with calcium scoring and CAD-RADS |
| 3. **Heart Failure** | HF classification, GDMT status, optimization recommendations |
| 4. **Valvular Disease** | Valve severity grading, intervention criteria assessment |
| 5. **Arrhythmia** | ECG interpretation, CHA2DS2-VASc, anticoagulation decision |
| 6. **Cardiac MRI** | Tissue characterization, LGE pattern analysis, parametric mapping |
| 7. **Stress Test** | Stress test interpretation, ischemic burden analysis |
| 8. **Prevention** | ASCVD risk, statin eligibility, CAC reclassification |
| 9. **Cardio-Oncology** | Cardiotoxicity surveillance, CTRCD monitoring |
| 10. **Reports** | Report generation and export (PDF, CSV, FHIR) |

**Sidebar:** Workflow filter, collection statistics, risk calculator quick access, export options.

---

## 17. Ingest Parsers

Seven domain-specific parsers (`src/ingest/`, 3,742 lines total) process cardiovascular data from external sources:

| Parser | Lines | Source | Target Collection |
|--------|-------|--------|-------------------|
| `pubmed_parser.py` | 477 | PubMed via NCBI E-utilities | cardio_literature |
| `clinical_trials_parser.py` | 713 | ClinicalTrials.gov API | cardio_trials |
| `imaging_parser.py` | 519 | Cardiac imaging report data | cardio_imaging |
| `ecg_parser.py` | 424 | ECG interpretation data | cardio_electrophysiology |
| `guideline_parser.py` | 555 | ACC/AHA/ESC guideline PDFs | cardio_guidelines |
| `device_parser.py` | 428 | FDA 510(k)/De Novo database | cardio_devices |
| `hemodynamics_parser.py` | 402 | Catheterization data | cardio_hemodynamics |

All parsers extend `BaseParser` (`base.py`, 224 lines) which provides:
- Chunking with configurable size and overlap
- Embedding generation via BGE-small-en-v1.5
- Batch Milvus insertion
- Deduplication by content hash
- Error handling and retry logic

---

## 18. Metrics & Monitoring

Prometheus metrics (`src/metrics.py`, 537 lines):

| Metric | Type | Description |
|--------|------|-------------|
| `cardio_queries_total` | Counter | Total queries processed |
| `cardio_query_duration_seconds` | Histogram | Query processing latency |
| `cardio_risk_calculations_total` | Counter | Risk calculator invocations by type |
| `cardio_gdmt_optimizations_total` | Counter | GDMT optimization requests |
| `cardio_collection_search_duration` | Histogram | Per-collection search latency |
| `cardio_cross_modal_triggers_total` | Counter | Cross-modal triggers fired |
| `cardio_workflow_executions_total` | Counter | Workflow executions by type |
| `cardio_export_requests_total` | Counter | Export requests by format |
| `cardio_errors_total` | Counter | Error count by type |

Endpoint: `GET /metrics` (Prometheus scrape target)

---

## 19. Scheduler

The APScheduler-based scheduler (`src/scheduler.py`, 612 lines) manages periodic data ingestion:

| Setting | Default | Description |
|---------|---------|-------------|
| `INGEST_SCHEDULE_HOURS` | 168 (weekly) | Hours between ingest runs |
| `INGEST_ENABLED` | False | Enable/disable automatic ingestion |

When enabled, the scheduler runs all 7 ingest parsers sequentially on the configured interval, updating collection data from PubMed, ClinicalTrials.gov, and other external sources.

---

## 20. Configuration

All configuration via `config/settings.py` (181 lines) using Pydantic BaseSettings with `CARDIO_` prefix environment variables:

### Key Configuration Groups

| Group | Settings |
|-------|----------|
| **Paths** | PROJECT_ROOT, DATA_DIR, CACHE_DIR, REFERENCE_DIR, RAG_PIPELINE_ROOT |
| **Milvus** | MILVUS_HOST (localhost), MILVUS_PORT (19530), 13 collection names |
| **Embeddings** | EMBEDDING_MODEL (BGE-small-en-v1.5), EMBEDDING_DIMENSION (384), EMBEDDING_BATCH_SIZE (32) |
| **LLM** | LLM_PROVIDER (anthropic), LLM_MODEL (claude-sonnet-4-6), ANTHROPIC_API_KEY |
| **RAG Search** | TOP_K_PER_COLLECTION (5), SCORE_THRESHOLD (0.4), 13 collection weights |
| **PubMed** | NCBI_API_KEY, PUBMED_MAX_RESULTS (5000) |
| **API Server** | API_HOST (0.0.0.0), API_PORT (8126) |
| **Streamlit** | STREAMLIT_PORT (8536) |
| **Metrics** | METRICS_ENABLED (True) |
| **Scheduler** | INGEST_SCHEDULE_HOURS (168), INGEST_ENABLED (False) |
| **Conversation** | MAX_CONVERSATION_CONTEXT (3) |
| **Citation** | CITATION_HIGH_THRESHOLD (0.75), CITATION_MEDIUM_THRESHOLD (0.60) |
| **CORS** | CORS_ORIGINS |
| **Limits** | MAX_REQUEST_SIZE_MB (10) |

### Startup Validation

The `validate()` method checks:
- Milvus host/port validity
- API key presence (warns if missing; search-only mode available)
- Embedding model presence
- Port conflict detection (API vs Streamlit)
- Collection weights sum to ~1.0 (tolerance 0.05)
- RAG pipeline root directory existence

---

## 21. Docker Deployment

`docker-compose.yml` defines the full stack:

| Service | Image | Port | Description |
|---------|-------|------|-------------|
| `etcd` | bitnami/etcd | 2379 | Milvus metadata store |
| `minio` | minio/minio | 9000, 9001 | Milvus object storage |
| `milvus` | milvusdb/milvus | 19530, 9091 | Vector database |
| `cardio-setup` | python:3.12-slim | - | Creates collections, seeds data (exits on completion) |
| `cardio-api` | python:3.12-slim | 8126 | FastAPI server |
| `cardio-ui` | python:3.12-slim | 8536 | Streamlit UI |

### Quickstart

```bash
# 1. Configure environment
cp .env.example .env
# Edit .env and set ANTHROPIC_API_KEY

# 2. Start all services
docker compose up -d

# 3. Watch setup progress
docker compose logs -f cardio-setup

# 4. Open the UI
open http://localhost:8536
```

---

## 22. Port Map

| Service | Port | Protocol |
|---------|------|----------|
| Streamlit UI | 8536 | HTTP |
| FastAPI API | 8126 | HTTP |
| Milvus gRPC | 19530 | gRPC |
| Milvus Health | 9091 | HTTP |
| MinIO API | 9000 | HTTP |
| MinIO Console | 9001 | HTTP |
| etcd | 2379 | HTTP |

---

## 23. Tech Stack

| Layer | Technology |
|-------|-----------|
| **Vector DB** | Milvus 2.4 (IVF_FLAT, COSINE) |
| **Embeddings** | BGE-small-en-v1.5 (384-dim, BAAI) |
| **LLM** | Claude Sonnet 4.6 (Anthropic) |
| **UI** | Streamlit |
| **API** | FastAPI + Uvicorn |
| **Configuration** | Pydantic BaseSettings |
| **Monitoring** | Prometheus metrics |
| **Scheduler** | APScheduler |
| **Containerization** | Docker Compose |
| **Compute** | NVIDIA DGX Spark |
| **Data Sources** | PubMed, ClinicalTrials.gov, ACC/AHA, ESC, FDA |

---

## 24. Cross-Agent Integration and Pediatric Oncology Focus

### Cross-Agent Endpoints

The Cardiology Intelligence Agent calls 4 sibling agents and exposes an `/v1/cardio/integrated-assessment` endpoint for coordinated multi-agent analysis:

- **Oncology Agent (8526)**: Cardiotoxicity correlation -- links CTRCD findings to active chemotherapy regimens, requests cumulative anthracycline dose data
- **Clinical Trial Agent (8537)**: Matches cardiology patients to active cardiovascular and cardio-oncology clinical trials
- **Biomarker Agent (8529)**: Longitudinal biomarker trending for troponin, BNP, GLS; threshold alert generation
- **Imaging Agent (8524)**: Cross-modal imaging correlation for echo, CMR, CT, and nuclear findings

### Pediatric Cardio-Oncology

The `pediatric_cardiotoxicity_assessment()` function provides specialized analysis for pediatric oncology patients:

- **Anthracycline cardiotoxicity**: 5-25% late-onset incidence in children; cumulative dose thresholds at 250 mg/m2 (surveillance) and 450 mg/m2 (high-risk)
- **Dexrazoxane cardioprotection**: Recommended for cumulative doxorubicin >=250 mg/m2 per COG guidelines
- **Pediatric echo monitoring**: Serial echocardiography with GLS tracking per COG LTFU guidelines (baseline, mid-treatment, end-of-treatment, 1/2/5/10-year follow-up)
- **Cross-agent coordination**: Oncology Agent provides treatment protocol context; Biomarker Agent tracks troponin/BNP trajectories; Imaging Agent correlates serial echo data; Trial Agent identifies pediatric cardioprotection studies

---

## 25. Future Roadmap

1. ~~**Test suite expansion**~~ -- 1,966 tests now passing across all modules
2. **NIM integration** -- On-device inference with NVIDIA NIM microservices
3. **FHIR interoperability** -- Full FHIR R4 integration for EHR connectivity
4. **Real-time ECG analysis** -- Streaming ECG interpretation pipeline
5. **Imaging AI integration** -- Direct DICOM analysis with on-device models
6. **Population analytics** -- Cohort-level cardiovascular risk dashboards
7. **Multi-language support** -- Clinical guidelines in multiple languages
8. **Wearable integration** -- Apple Watch, Fitbit cardiovascular data ingestion
