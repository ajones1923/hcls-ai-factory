# Cardiology Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Cardiology Intelligence Agent extends the HCLS AI Factory platform to deliver RAG-powered cardiovascular clinical decision support. It synthesizes cardiac imaging, electrophysiology, hemodynamics, heart failure management, valvular disease, preventive cardiology, interventional data, and cardio-oncology surveillance into guideline-aligned clinical recommendations using ACC/AHA/ESC evidence.

The system implements **6 validated cardiovascular risk calculators** (ASCVD, MAGGIC, EuroSCORE II, CHA2DS2-VASc, HAS-BLED, HEART), optimizes **guideline-directed medical therapy (GDMT)** for heart failure across 7 therapy classes, and provides **11 clinical workflows** covering the highest-impact cardiovascular use cases — all backed by 13 Milvus vector collections containing cardiac imaging protocols, ECG data, hemodynamic measurements, genomic-cardiac correlations, and clinical guidelines.

The platform enables cross-modal queries like *"This patient has an LVEF of 28%, LBBB on ECG, and a LMNA pathogenic variant — what is the optimal management strategy?"* that simultaneously search imaging protocols, guideline recommendations, genomic correlations, and clinical evidence.

### Key Results

| Metric | Value |
|---|---|
| Total Python LOC | **28,189** |
| Milvus collections | **13** (12 cardiology-specific + 1 shared genomic_evidence) |
| Cardiovascular conditions | **45** in knowledge graph |
| Cardiac biomarkers | **29** |
| Drug classes | **32** |
| Cardiovascular genes | **56** |
| Imaging protocols | **27** (12 echo, 7 CMR, 5 nuclear, 3 CT) |
| Guideline documents | **20** (ACC/AHA/ESC) |
| Guideline recommendations | **63** structured with class/level of evidence |
| Risk calculators | **6** validated scoring systems |
| Clinical workflows | **11** |
| Cross-modal imaging triggers | **18** genomic trigger patterns |
| GDMT therapy classes | **7** (including finerenone, omecamtiv, sotagliflozin) |
| Entity aliases | **167** abbreviation mappings |
| Test suite | **1,966** tests (100% pass, <1.2s runtime) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Cardiology Agent Role |
|---|---|
| **DataStore** | Raw files: PubMed XML, ClinicalTrials.gov JSON, imaging protocol specs, ECG templates, guideline PDFs |
| **DataEngine** | 7 ingest parsers (PubMed, trials, imaging, ECG, guideline, device, hemodynamics) |
| **DataBase** | 13 Milvus collections + knowledge graph (45 conditions, 56 genes, 29 biomarkers, 32 drug classes) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + 6 risk calculators + GDMT optimizer + cross-modal triggers |
| **AgentEngine** | CardiologyAgent orchestrator + Streamlit UI (10 tabs) + FastAPI REST |

### 2.2 System Diagram

```
                          EXTERNAL USERS
                               |
                    +----------+----------+
                    |                     |
              +-----+------+      +------+-----+
              | Streamlit  |      | REST API   |
              | UI :8536   |      | :8126      |
              +-----+------+      +------+-----+
                    |                     |
                    +----------+----------+
                               |
                    +----------+----------+
                    |  Agent Orchestrator  |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Query       |    | Workflow      |    | Clinical      |
   | Expansion   |    | Engine (11)   |    | Engines       |
   | 167 aliases |    |               |    |               |
   +------+------+    +-------+-------+    | 6 Risk Calcs  |
          |                    |           | GDMT Optimizer |
          |                    |           | Cross-Modal    |
          |                    |           +-------+-------+
          +--------------------+--------------------+
                               |
                    +----------+----------+
                    |    RAG Engine       |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Knowledge   |    | Milvus        |    | LLM           |
   | Graph       |    | Vector DB     |    | (Claude 4.6)  |
   | 45 conds    |    | 13 Collections|    |               |
   | 56 genes    |    |               |    |               |
   +-------------+    +---------------+    +---------------+
```

---

## 3. Data Collections — Actual State

### 3.1 Collection Catalog

| # | Collection | Est. Records | Weight | Primary Use |
|---|---|---|---|---|
| 1 | `cardio_literature` | 8,000 | 0.12 | Evidence synthesis, guideline grounding |
| 2 | `cardio_trials` | 3,000 | 0.10 | Trial evidence, outcome data |
| 3 | `cardio_imaging` | 2,500 | 0.10 | Echo, CMR, nuclear, CT protocols |
| 4 | `cardio_ecg` | 1,500 | 0.08 | ECG patterns, arrhythmia recognition |
| 5 | `cardio_hemodynamics` | 1,000 | 0.07 | Catheterization, pressure-volume data |
| 6 | `cardio_guidelines` | 1,200 | 0.10 | ACC/AHA/ESC recommendations |
| 7 | `cardio_drugs` | 800 | 0.08 | GDMT protocols, drug interactions |
| 8 | `cardio_devices` | 600 | 0.06 | ICD, CRT, LVAD, TAVR data |
| 9 | `cardio_genomics` | 1,500 | 0.08 | Cardiac gene panels, channelopathies |
| 10 | `cardio_biomarkers` | 500 | 0.06 | Troponin, BNP, hs-CRP interpretation |
| 11 | `cardio_risk_scores` | 200 | 0.05 | Validation cohort data for calculators |
| 12 | `cardio_rehab` | 400 | 0.04 | Cardiac rehabilitation protocols |
| 13 | `genomic_evidence` | 3,560,000 | 0.06 | Shared genomic variant context |

### 3.2 Index Configuration (all collections)

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 1024 (literature, trials), 256 (others) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Risk Calculators

### 4.1 ASCVD 10-Year Risk (Pooled Cohort Equations)

| Parameter | Input |
|---|---|
| Age | 40-79 years |
| Sex | Male / Female |
| Race | White / African American |
| Total cholesterol | mg/dL |
| HDL cholesterol | mg/dL |
| Systolic BP | mmHg |
| BP treatment | Yes / No |
| Diabetes | Yes / No |
| Smoking | Yes / No |
| **Output** | 10-year ASCVD risk (%), risk category (low/borderline/intermediate/high) |

### 4.2 MAGGIC Heart Failure Mortality

| Input Variables | Count |
|---|---|
| Age, sex, LVEF, NYHA class, SBP, BMI, creatinine, smoking, diabetes, COPD, HF duration, beta-blocker, ACE-I/ARB | 13 variables |
| **Output** | 1-year and 3-year mortality risk (%), risk stratification |

### 4.3 EuroSCORE II (Cardiac Surgery)

Logistic regression model with 18 variables predicting operative mortality for cardiac surgery. Covers patient factors (age, sex, renal function, extracardiac arteriopathy, mobility, cardiac surgery history, chronic lung disease, endocarditis, neurological status, diabetes on insulin) and cardiac factors (NYHA class, CCS angina, LVEF, recent MI, pulmonary hypertension, urgency, weight of procedure, surgery on thoracic aorta).

### 4.4 CHA2DS2-VASc (Stroke Risk in AF)

| Component | Points |
|---|---|
| Congestive heart failure | 1 |
| Hypertension | 1 |
| Age >= 75 | 2 |
| Diabetes mellitus | 1 |
| Stroke/TIA/thromboembolism | 2 |
| Vascular disease | 1 |
| Age 65-74 | 1 |
| Sex category (female) | 1 |
| **Output** | Score (0-9), annual stroke risk (%), anticoagulation recommendation |

### 4.5 HAS-BLED (Bleeding Risk)

| Component | Points |
|---|---|
| Hypertension (uncontrolled, SBP >160) | 1 |
| Abnormal renal/liver function | 1-2 |
| Stroke history | 1 |
| Bleeding predisposition | 1 |
| Labile INR | 1 |
| Elderly (>65) | 1 |
| Drugs/alcohol | 1-2 |
| **Output** | Score (0-9), annual major bleeding risk (%), risk category |

### 4.6 HEART Score (Chest Pain)

| Component | Points |
|---|---|
| History | 0-2 |
| ECG | 0-2 |
| Age | 0-2 |
| Risk factors | 0-2 |
| Troponin | 0-2 |
| **Output** | Score (0-10), 6-week MACE risk (%), disposition recommendation (discharge / observe / intervene) |

---

## 5. GDMT Optimizer

The GDMT Optimizer implements guideline-directed medical therapy titration for heart failure with reduced ejection fraction (HFrEF) across 7 therapy classes:

| Therapy Class | Target | Key Agents |
|---|---|---|
| **ACEi/ARB/ARNI** | Blood pressure, remodeling | Sacubitril/valsartan (ARNI preferred) |
| **Beta-blocker** | Heart rate, remodeling | Carvedilol, metoprolol succinate, bisoprolol |
| **MRA** | Potassium-sparing diuresis | Spironolactone, eplerenone |
| **SGLT2i** | Cardiovascular mortality | Dapagliflozin, empagliflozin |
| **Finerenone** | Non-steroidal MRA (CKD+HF) | Finerenone (FIDELIO/FIGARO evidence) |
| **Omecamtiv mecarbil** | Cardiac myosin activator | Omecamtiv (GALACTIC-HF evidence) |
| **Sotagliflozin** | Dual SGLT1/2 inhibitor | Sotagliflozin (SOLOIST/SCORED evidence) |

The optimizer checks contraindications, current doses vs. target doses, lab monitoring requirements (potassium, creatinine, eGFR), and generates a step-by-step titration plan.

---

## 6. Clinical Workflows

| # | Workflow | Clinical Question |
|---|---|---|
| 1 | Chest Pain Triage | "Risk-stratify this chest pain presentation using HEART score" |
| 2 | Heart Failure Management | "Optimize GDMT for this HFrEF patient" |
| 3 | Atrial Fibrillation | "CHA2DS2-VASc and HAS-BLED for anticoagulation decision" |
| 4 | Valvular Assessment | "Is this patient a TAVR or surgical AVR candidate?" |
| 5 | Cardiac Imaging | "Which imaging modality is optimal for this clinical scenario?" |
| 6 | Preventive Cardiology | "ASCVD risk and statin therapy recommendation" |
| 7 | Cardiomyopathy Genetics | "Evaluate this DCM/HCM genetic panel result" |
| 8 | Device Therapy | "Does this patient meet criteria for ICD or CRT?" |
| 9 | Perioperative Risk | "Estimate surgical risk with EuroSCORE II" |
| 10 | Cardio-Oncology | "Monitor for cardiotoxicity during anthracycline therapy" |
| 11 | Cardiac Rehabilitation | "Design a rehab protocol for post-CABG recovery" |

---

## 7. Cross-Modal Imaging Triggers

The agent detects 18 genomic patterns that trigger specific cardiac imaging recommendations:

| Genomic Finding | Triggered Imaging | Rationale |
|---|---|---|
| LMNA pathogenic variant | CMR with late gadolinium enhancement | DCM + conduction disease risk |
| MYH7/MYBPC3 variant | Echocardiogram + CMR | HCM screening |
| TTN truncating variant | Echocardiogram + CMR | DCM assessment |
| SCN5A variant | 12-lead ECG + signal-averaged ECG | Brugada syndrome |
| KCNQ1/KCNH2 variant | 12-lead ECG + exercise stress | Long QT syndrome |
| PKP2/DSP variant | CMR + signal-averaged ECG | ARVC evaluation |

---

## 8. Multi-Collection RAG Engine

### 8.1 Search Flow

```
User Query: "Optimize GDMT for 62M, LVEF 25%, NYHA III, CKD stage 3"
    │
    ├── 1. Embed query with BGE asymmetric prefix               [< 5 ms]
    │
    ├── 2. Parallel search across 13 collections (top-5 each)   [12-18 ms]
    │   ├── cardio_guidelines:  HFrEF GDMT recommendations      (score: 0.84-0.92)
    │   ├── cardio_drugs:       ARNI/BB/MRA/SGLT2i protocols     (score: 0.80-0.88)
    │   ├── cardio_literature:  DAPA-HF, EMPEROR-Reduced data    (score: 0.78-0.86)
    │   ├── cardio_biomarkers:  BNP/NT-proBNP monitoring         (score: 0.72-0.80)
    │   └── cardio_trials:      HFrEF landmark trials            (score: 0.70-0.82)
    │
    ├── 3. Query expansion + knowledge augmentation              [< 1 ms]
    │
    ├── 4. GDMT Optimizer: titration plan generation             [< 50 ms]
    │
    └── 5. Stream Claude Sonnet 4.6 response                    [~22-26 sec]
           GDMT plan with dose targets, monitoring schedule,
           CKD-specific adjustments, and guideline citations
```

---

## 9. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 9.1 Risk Calculator Performance

| Calculator | Latency | Validated Against |
|---|---|---|
| ASCVD (Pooled Cohort) | <10 ms | ACC/AHA risk calculator |
| MAGGIC | <15 ms | Original MAGGIC publication |
| EuroSCORE II | <20 ms | euroscore.org |
| CHA2DS2-VASc | <5 ms | ESC guidelines |
| HAS-BLED | <5 ms | ESC guidelines |
| HEART | <5 ms | HEART Pathway validation study |
| **All 6 calculators** | **<60 ms** | |

### 9.2 RAG Query Performance

| Operation | Latency |
|---|---|
| Full query (retrieve + Claude generate) | ~24 sec |
| Streaming query (time to first token) | ~3 sec |
| 13-collection parallel search | 12-18 ms |
| GDMT optimization | <50 ms |

---

## 10. Infrastructure

### 10.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4, localhost:19530 |
| Embeddings | BGE-small-en-v1.5 (BAAI) — 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit (port 8536, NVIDIA black/green theme) |
| REST API | FastAPI + Uvicorn (port 8126) |
| Configuration | Pydantic BaseSettings |
| Hardware target | NVIDIA DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 10.2 Service Ports

| Port | Service |
|---|---|
| 8126 | FastAPI REST API |
| 8536 | Streamlit Chat UI |
| 19530 | Milvus vector database (shared) |

### 10.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — adds 12 owned collections alongside existing `genomic_evidence` (3.56M vectors, read-only) |
| `ANTHROPIC_API_KEY` | Shared Anthropic API key |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 11. Knowledge Graph

### 11.1 Cardiovascular Conditions (45 entries)

Organized into 8 categories: heart failure (8), arrhythmia (7), coronary artery disease (5), valvular disease (6), cardiomyopathy (5), congenital heart disease (4), vascular disease (5), pericardial disease (5).

### 11.2 Cardiac Genes (56 entries)

| Category | Count | Key Genes |
|---|---|---|
| Cardiomyopathy | 18 | MYH7, MYBPC3, TTN, LMNA, DES, PLN, RBM20, FLNC |
| Arrhythmia (channelopathy) | 14 | SCN5A, KCNQ1, KCNH2, KCNJ2, RYR2, CASQ2, CALM1-3 |
| ARVC | 6 | PKP2, DSP, DSG2, DSC2, JUP, TMEM43 |
| Aortopathy | 5 | FBN1, TGFBR1, TGFBR2, SMAD3, ACTA2 |
| Lipid metabolism | 5 | LDLR, PCSK9, APOB, LDLRAP1, ABCG5 |
| Cardiac development | 4 | NKX2-5, GATA4, TBX5, TBX20 |
| Other | 4 | TNNT2, TNNI3, ACTC1, MYL2 |

### 11.3 Guidelines (51 structured recommendations)

All recommendations include ACC/AHA class (I, IIa, IIb, III) and level of evidence (A, B-R, B-NR, C-LD, C-EO).

---

## 12. Demo Scenarios

### 12.1 Validated Demo Queries

**1. "62-year-old male, LVEF 25%, NYHA Class III, on lisinopril and metoprolol — optimize GDMT"**
- GDMT Optimizer: Switch lisinopril to sacubitril/valsartan, add dapagliflozin, add spironolactone
- Titration plan with target doses, lab monitoring schedule

**2. "Calculate ASCVD risk: 55F, TC 240, HDL 45, SBP 138, on BP meds, non-diabetic, non-smoker, white"**
- ASCVD Calculator: 10-year risk with statin recommendation per ACC/AHA guidelines

**3. "New-onset AF, CHA2DS2-VASc 4, HAS-BLED 2 — anticoagulation strategy"**
- Dual calculator with anticoagulation recommendation, DOAC vs. warfarin guidance

**4. "LMNA p.R190W variant found in DCM patient — what cardiac workup is needed?"**
- Cross-modal trigger: CMR + ICD evaluation + family screening cascade
- Genomic knowledge: LMNA → high risk of sudden cardiac death, conduction disease

**5. "Chest pain, HEART score inputs: typical history, normal ECG, age 58, 3 risk factors, troponin normal"**
- HEART Score: Score calculation with 6-week MACE risk and disposition recommendation

---

## 13. File Structure (Actual)

```
cardiology_intelligence_agent/
├── src/
│   ├── agent.py                     # Agent orchestrator
│   ├── models.py                    # Pydantic models + 16 enums
│   ├── collections.py               # 13 Milvus collection schemas
│   ├── rag_engine.py                # Multi-collection RAG (1,589 LOC)
│   ├── clinical_workflows.py        # 11 workflows (2,445 LOC)
│   ├── risk_calculators.py          # 6 validated calculators (2,397 LOC)
│   ├── gdmt_optimizer.py            # GDMT titration engine (2,457 LOC)
│   ├── cross_modal.py               # 18 imaging trigger patterns (1,734 LOC)
│   ├── knowledge.py                 # Knowledge graph (1,431 LOC)
│   ├── query_expansion.py           # 167 aliases (2,025 LOC)
│   ├── metrics.py                   # Prometheus metrics
│   ├── export.py                    # Report generation
│   └── ingest/
│       ├── pubmed_parser.py
│       ├── trials_parser.py
│       ├── imaging_parser.py
│       ├── ecg_parser.py
│       ├── guideline_parser.py
│       ├── device_parser.py
│       └── hemodynamics_parser.py
├── app/
│   └── cardio_ui.py                # Streamlit (10 tabs, NVIDIA theme)
├── api/
│   └── main.py                     # FastAPI REST server
├── config/
│   └── settings.py                 # Pydantic BaseSettings
├── tests/                          # 1,966 tests
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

**43 files | ~28,189 lines of code | Apache 2.0**

---

## 14. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | 13 collections, knowledge graph, 6 risk calculators, GDMT optimizer, 11 workflows |
| **Phase 2: Data** | Complete | 45 conditions, 56 genes, 29 biomarkers, 32 drug classes, 27 imaging protocols, 63 recommendations |
| **Phase 3: RAG Integration** | Complete | Multi-collection parallel search, Claude Sonnet 4.6 streaming |
| **Phase 4: Testing** | Complete | 1,966 tests, 100% pass, <1.2s runtime |
| **Phase 5: UI + Demo** | Complete | 10-tab Streamlit UI, 5 demo scenarios validated |

---

## 15. Relationship to HCLS AI Factory

The Cardiology Intelligence Agent demonstrates the **multi-modal clinical extension** of the HCLS AI Factory architecture, integrating imaging, electrophysiology, hemodynamics, and genomics into a single decision support platform.

- **Same Milvus instance** — 12 new owned collections alongside existing `genomic_evidence` (3.56M vectors, read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($4,699)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The cross-modal trigger system uniquely bridges genomics (Stage 1) and cardiac imaging, enabling queries like *"This TTN truncating variant patient needs echocardiographic follow-up"* — connecting molecular findings to clinical imaging protocols automatically.

---

## 16. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Cardiology Intelligence Agent is a clinical decision support research tool for cardiovascular medicine. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Risk calculator outputs should be validated against institutional protocols. Apache 2.0 License.
