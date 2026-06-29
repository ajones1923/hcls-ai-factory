# Democratizing Cardiovascular AI: A Multi-Collection RAG Architecture and Product Requirements for the Cardiology Intelligence Agent

**Author:** Adam Jones
**Date:** March 2026
**Version:** 0.1.0 (Pre-Implementation)
**License:** Apache 2.0

Part of the HCLS AI Factory -- an end-to-end precision medicine platform.
https://github.com/ajones1923/hcls-ai-factory

---

## Abstract

Cardiovascular disease (CVD) remains the leading cause of death globally, claiming approximately 17.9 million lives annually -- representing 32% of all deaths worldwide. Despite exponential growth in cardiovascular AI research (over 4,500 publications in 2025 alone), critical barriers persist in clinical translation: fragmented evidence across modalities, siloed data systems, lack of integrated genomic-imaging correlation, and prohibitive infrastructure costs that limit advanced cardiac AI to elite academic centers.

This paper presents the architectural design, clinical rationale, and product requirements for the Cardiology Intelligence Agent -- a multi-collection retrieval-augmented generation (RAG) system purpose-built for cardiovascular medicine. The agent will unify 12 specialized Milvus vector collections spanning cardiac imaging (echocardiography, cardiac CT, cardiac MRI, nuclear cardiology), electrophysiology (12-lead ECG, Holter monitoring, device interrogation), hemodynamics, heart failure management, valvular heart disease, preventive cardiology, interventional cardiology, and cardio-oncology -- alongside a shared genomic_evidence collection containing 3.5 million variant vectors from the HCLS AI Factory genomics pipeline.

The system extends the proven multi-collection RAG architecture established by five existing intelligence agents in the HCLS AI Factory (Precision Biomarker, Precision Oncology, CAR-T, Imaging, and Autoimmune agents), adapting it with cardiology-specific clinical workflows, risk calculators, cross-modal imaging-genomics triggers, and structured reporting aligned with ACC/AHA guidelines. Eight reference clinical workflows will cover the highest-impact cardiovascular use cases: coronary artery disease assessment, heart failure classification, valvular disease quantification, arrhythmia detection, cardiac MRI tissue characterization, stress test interpretation, preventive risk stratification, and cardio-oncology surveillance.

The agent will deploy on a single NVIDIA DGX Spark ($4,699) using BGE-small-en-v1.5 embeddings (384-dimensional, IVF_FLAT, COSINE), Claude Sonnet 4.6 for evidence synthesis, and four NVIDIA NIM microservices for on-device inference. Licensed under Apache 2.0, the platform will democratize access to integrated cardiovascular intelligence that currently requires multi-million-dollar institutional investments in informatics infrastructure.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [The Cardiovascular Data Challenge](#2-the-cardiovascular-data-challenge)
3. [Clinical Landscape and Market Analysis](#3-clinical-landscape-and-market-analysis)
4. [Existing HCLS AI Factory Architecture](#4-existing-hcls-ai-factory-architecture)
5. [Cardiology Intelligence Agent Architecture](#5-cardiology-intelligence-agent-architecture)
6. [Milvus Collection Design](#6-milvus-collection-design)
7. [Clinical Workflows](#7-clinical-workflows)
8. [Cross-Modal Integration](#8-cross-modal-integration)
9. [NIM Integration Strategy](#9-nim-integration-strategy)
10. [Knowledge Graph Design](#10-knowledge-graph-design)
11. [Query Expansion and Retrieval Strategy](#11-query-expansion-and-retrieval-strategy)
12. [API and UI Design](#12-api-and-ui-design)
13. [Clinical Decision Support Engines](#13-clinical-decision-support-engines)
14. [Reporting and Interoperability](#14-reporting-and-interoperability)
15. [Product Requirements Document](#15-product-requirements-document)
16. [Data Acquisition Strategy](#16-data-acquisition-strategy)
17. [Validation and Testing Strategy](#17-validation-and-testing-strategy)
18. [Regulatory Considerations](#18-regulatory-considerations)
19. [DGX Compute Progression](#19-dgx-compute-progression)
20. [Implementation Roadmap](#20-implementation-roadmap)
21. [Risk Analysis](#21-risk-analysis)
22. [Competitive Landscape](#22-competitive-landscape)
23. [Discussion](#23-discussion)
24. [Conclusion](#24-conclusion)
25. [References](#25-references)

---

## 1. Introduction

### 1.1 The Cardiovascular Disease Burden

Cardiovascular disease represents the single largest cause of mortality and morbidity worldwide. According to the World Health Organization and the American Heart Association (AHA) 2025 Heart Disease and Stroke Statistics Update:

- **17.9 million deaths annually** -- 32% of all global deaths
- **523 million people** currently living with cardiovascular disease
- **$407 billion** in direct and indirect costs in the United States alone (2024)
- **Ischemic heart disease** and **stroke** are the top two causes of death globally, responsible for a combined 15.2 million deaths per year
- **Heart failure** affects 64.3 million people globally, with 5-year mortality rates exceeding 50%
- **Atrial fibrillation** prevalence has reached 59 million cases worldwide, with a projected 72 million by 2030

Despite these staggering numbers, cardiovascular AI adoption in clinical practice remains in early stages. A 2025 survey by the American College of Cardiology (ACC) found that while 89% of cardiologists believe AI will transform their practice, only 23% currently use AI-assisted tools in routine clinical decision-making. The gap between potential and adoption is driven by three factors: fragmented data, prohibitive cost, and lack of integrated systems.

### 1.2 The Opportunity for Integrated Cardiovascular Intelligence

The cardiovascular domain is uniquely suited for an integrated AI intelligence agent because:

1. **Multi-modal data convergence**: Cardiology routinely integrates imaging (echo, CT, MRI, nuclear), electrophysiology (ECG, Holter), hemodynamics (catheterization), biomarkers (troponin, BNP, lipids), genomics (familial cardiomyopathy, channelopathies), and clinical risk scores -- all of which must be synthesized for optimal patient care.

2. **Established clinical guidelines**: ACC/AHA guidelines provide structured frameworks for risk stratification, treatment decisions, and follow-up protocols that can be encoded into clinical decision support logic.

3. **Quantitative measurement-rich**: Unlike many medical specialties, cardiology generates highly quantitative data (ejection fraction percentages, valve gradients in mmHg, calcium scores in Agatston units, QTc intervals in milliseconds) that is well-suited for structured RAG retrieval and comparison.

4. **Strong cross-modal triggers**: A cardiac imaging finding (e.g., unexplained left ventricular hypertrophy) frequently triggers genomic workup (hypertrophic cardiomyopathy gene panel), which may lead to drug therapy optimization -- a natural fit for the HCLS AI Factory's three-stage pipeline.

5. **Large addressable market**: Every health system has a cardiology department. Cardiovascular AI is projected to reach $2.8 billion by 2028 (CAGR 42.3%), making this the highest-value vertical for an intelligence agent.

### 1.3 Our Contribution

This paper presents the complete architectural blueprint and product requirements for the Cardiology Intelligence Agent, the sixth domain-specific intelligence agent in the HCLS AI Factory platform. Our contributions include:

- A **12-collection Milvus vector schema** designed for the full spectrum of cardiovascular data: imaging, electrophysiology, hemodynamics, heart failure, valvular disease, preventive cardiology, interventional procedures, cardio-oncology, guidelines, trials, devices, and literature
- **Eight reference clinical workflows** covering coronary artery disease, heart failure, valvular disease, arrhythmia, cardiac MRI, stress testing, preventive risk stratification, and cardio-oncology surveillance
- A **cardiology knowledge graph** with structured data on 30+ cardiovascular conditions, 15+ imaging modalities/protocols, 20+ cardiac biomarkers, 25+ drug classes, and 50+ ACC/AHA guideline recommendations
- **Cross-modal triggers** linking cardiac imaging findings to genomic workup (familial hypercholesterolemia, cardiomyopathy gene panels, channelopathies) via the shared `genomic_evidence` collection
- **Clinical decision support engines** implementing validated risk calculators (ASCVD, HEART, CHA₂DS₂-VASc, HAS-BLED, MAGGIC, EuroSCORE II)
- A comprehensive **product requirements document** with user stories, acceptance criteria, and implementation prioritization
- **Deployment on a single NVIDIA DGX Spark** ($4,699), maintaining the platform's commitment to accessible AI

---

## 2. The Cardiovascular Data Challenge

### 2.1 Data Fragmentation in Cardiovascular Medicine

Cardiovascular clinical practice generates data across at least fifteen distinct categories, each with its own structure, vocabulary, source systems, and update cadences:

1. **Cardiac Imaging Literature** -- PubMed abstracts, JACC imaging supplements, Circulation reviews, conference proceedings (ACC, AHA, ESC, SCMR, SCCT, ASE). Over 4,500 publications in 2025 alone on cardiovascular AI.

2. **Echocardiography Data** -- Chamber dimensions (LVIDd, LVIDs, LA volume), systolic function (LVEF, GLS, TAPSE, S'), diastolic function (E/A ratio, E/e', deceleration time), valve assessments (gradients, regurgitant volumes, effective orifice area), strain imaging (longitudinal, circumferential, radial). Structured per ASE guidelines.

3. **Cardiac CT Data** -- Coronary artery calcium scores (Agatston, volume, mass scores), CTA stenosis grading (CAD-RADS 0-5), plaque characterization (calcified, non-calcified, mixed, vulnerable), fractional flow reserve from CT (FFR-CT), myocardial perfusion CT, cardiac structural assessment.

4. **Cardiac MRI Data** -- Volumes and function (LVEF, RVEF, indexed volumes), tissue characterization (T1 mapping, T2 mapping, ECV), late gadolinium enhancement (LGE) patterns (ischemic vs non-ischemic), perfusion (stress, rest, MPR), feature tracking strain, 4D flow, parametric mapping.

5. **Nuclear Cardiology** -- SPECT MPI (stress/rest perfusion, TID ratio, transient ischemic dilation), PET MPI (MBF quantification, CFR), MUGA (LVEF for cardio-oncology), cardiac amyloid scintigraphy (Tc-99m PYP, DPD), FDG-PET for sarcoidosis/endocarditis.

6. **Electrophysiology Data** -- 12-lead ECG interpretation (rhythm, intervals, morphology, axis, ischemic changes), Holter/event monitor data (arrhythmia burden, HR variability), device interrogation (ICD, pacemaker, CRT), EP study results (ablation lesion sets, activation mapping), QTc monitoring.

7. **Hemodynamic Data** -- Right heart catheterization (PA pressures, wedge pressure, cardiac output, PVR, SVR), left heart catheterization (LVEDP, aortic valve gradient), coronary angiography (TIMI flow, stenosis percentage, FFR/iFR), structural intervention hemodynamics.

8. **Heart Failure Data** -- NYHA classification, ACC/AHA staging (A-D), biomarkers (NT-proBNP, BNP, hs-troponin), GDMT titration (beta-blocker, ACEi/ARB/ARNI, MRA, SGLT2i), device therapy (ICD, CRT), transplant evaluation (MELD-XI, SHFM, HFSS), LVAD management.

9. **Valvular Heart Disease** -- Severity grading (mild/moderate/severe per ASE criteria), hemodynamic quantification (EOA, DVI, PHT, regurgitant fraction), intervention criteria (surgical AVR vs TAVR, mitral repair vs replacement), surveillance protocols.

10. **Preventive Cardiology** -- Lipid panels (LDL-C, HDL-C, TG, Lp(a), ApoB), risk calculators (PCE/ASCVD, HEART, MESA, Reynolds), statin eligibility, PCSK9i criteria, coronary calcium for risk reclassification, inflammatory markers (hsCRP, IL-6), metabolic parameters.

11. **Interventional Cardiology** -- PCI procedural data (stent type, vessel, lesion characteristics, TIMI flow), structural intervention data (TAVR, MitraClip, WATCHMAN, PFO closure), complications, antiplatelet/anticoagulant management, CTO techniques.

12. **Cardio-Oncology** -- Baseline cardiac assessment protocols (pre-anthracycline, pre-immunotherapy, pre-targeted therapy), surveillance schedules (GLS monitoring, troponin trends), cardiotoxicity detection (CTRCD definitions), risk scores (HFA-ICOS), cardioprotective strategies.

13. **Clinical Trials** -- ClinicalTrials.gov cardiovascular entries (16,000+ active/completed), landmark trial results (PARADIGM-HF, DAPA-HF, EMPEROR-Reduced, ISCHEMIA, REVIVED-BCIS2), outcome data, subgroup analyses.

14. **Cardiovascular Devices** -- FDA-cleared AI/ML cardiac devices (ECG interpretation, echo measurement, CT calcium scoring, arrhythmia detection), implantable devices (ICDs, pacemakers, CRT, LVAD), wearable cardiac monitors (smartwatch ECG, continuous glucose monitors for cardiometabolic risk).

15. **Cardiovascular Genomics** -- Familial hypercholesterolemia (LDLR, PCSK9, APOB), hypertrophic cardiomyopathy (MYH7, MYBPC3, TNNT2, TNNI3), dilated cardiomyopathy (TTN, LMNA, RBM20), arrhythmogenic cardiomyopathy (PKP2, DSP, DSG2), channelopathies (SCN5A, KCNQ1, KCNH2, RYR2), aortopathies (FBN1, TGFBR1/2, SMAD3, ACTA2).

### 2.2 Why Existing Tools Fall Short

Current approaches to cardiovascular intelligence fail to address this fragmentation:

| Approach | Limitation |
|---|---|
| **PubMed search** | Keyword-based; misses semantic connections; no cross-modal integration; no structured clinical data |
| **UpToDate / DynaMed** | Expert-curated but static; no patient-specific reasoning; no imaging or genomic data integration; subscription-based |
| **Commercial CVIS** (Solas, Lumedx) | Vendor-locked; limited to institutional data; no literature integration; $500K-$2M+ implementation |
| **EHR-integrated CDS** (Epic BestPractice) | Rule-based; cannot synthesize unstructured evidence; no imaging AI; limited to institutional data |
| **General AI assistants** | No citation provenance; hallucination risk; no structured cardiovascular data; not FDA-aligned |
| **Imaging-only AI** (Arterys, HeartFlow) | Single-modality; no genomic integration; no guideline reasoning; cloud-dependent; expensive per-scan pricing |

The Cardiology Intelligence Agent addresses all six limitations simultaneously by combining multi-collection vector search, cross-modal genomic triggers, validated clinical decision support logic, and guideline-grounded LLM synthesis -- all on a $4,699 desktop device.

### 2.3 The Case for Multi-Collection RAG in Cardiology

A cardiologist evaluating a patient with new-onset heart failure must simultaneously consider:

- **Imaging data**: Echocardiographic LVEF, cardiac MRI showing late gadolinium enhancement pattern (ischemic vs non-ischemic)
- **Electrophysiology**: ECG showing LBBB morphology (CRT candidacy), QTc for drug safety
- **Biomarkers**: NT-proBNP trend, hs-troponin for myocardial injury, iron studies for IV iron candidacy
- **Genomics**: If age < 50 or family history, cardiomyopathy gene panel (TTN, LMNA, MYH7)
- **Guidelines**: ACC/AHA heart failure guidelines for GDMT initiation sequence
- **Trials**: DAPA-HF, EMPEROR-Reduced for SGLT2i evidence; PARADIGM-HF for ARNI
- **Risk scores**: MAGGIC score for prognosis, Seattle Heart Failure Model for transplant timing

No existing tool synthesizes all seven data dimensions into a single clinical narrative. A multi-collection RAG architecture -- where each data dimension has its own optimized Milvus collection with domain-specific schema fields -- enables parallel retrieval across all dimensions with a single query, followed by LLM synthesis into a coherent clinical recommendation.

---

## 3. Clinical Landscape and Market Analysis

### 3.1 Cardiovascular AI Market

The global cardiovascular AI market demonstrates exceptional growth dynamics:

| Metric | Value | Source |
|---|---|---|
| Market size (2024) | $1.4 billion | Grand View Research |
| Projected size (2028) | $2.8 billion | Grand View Research |
| CAGR (2024-2028) | 42.3% | Grand View Research |
| FDA-cleared cardiac AI devices (cumulative) | 180+ | FDA AI/ML database |
| Active cardiovascular AI clinical trials | 450+ | ClinicalTrials.gov |
| Annual CV AI publications | 4,500+ | PubMed (2025) |
| US cardiology practices | 25,000+ | ACC Census |
| US cardiologists | 35,000+ | AHA Statistics |
| Global cardiologists | 200,000+ | WHO estimates |

### 3.2 Competitive Analysis

| Competitor | Strengths | Gaps |
|---|---|---|
| **HeartFlow** (FFR-CT) | FDA-cleared, clinical validation | Single modality (CT), no genomics, cloud-only, $1,100/scan |
| **Eko Health** (ECG + auscultation) | Point-of-care, stethoscope integration | Limited to murmur/arrhythmia screening, no imaging AI |
| **Viz.ai** (Stroke + cardiac) | Real-time triage, HIPAA-compliant | Narrow scope (LVO stroke, PE), SaaS pricing |
| **Ultromics** (EchoGo) | FDA-cleared echo AI, GLS analysis | Echo-only, no cross-modal integration |
| **Cleerly** (Coronary CTA) | Plaque quantification, FDA pathway | CT-only, cloud-only, per-patient pricing |
| **Caption Health** (Echo guidance) | AI-guided acquisition, GE partnership | Acquisition guidance only, no interpretation |
| **Tempus** (Cardiology platform) | Multi-modal data, genomics integration | Proprietary, expensive, cloud-dependent |

**Our differentiation**: The Cardiology Intelligence Agent is the only system that combines (1) multi-modal imaging AI, (2) genomic integration, (3) literature RAG, (4) guideline-aligned clinical decision support, (5) validated risk calculators, and (6) on-device deployment -- all in an open-source, $4,699 package. No competitor addresses more than two of these six dimensions.

### 3.3 Target Users

| User Segment | Use Case | Pain Point Addressed |
|---|---|---|
| **Community cardiologists** | Evidence-based decision support | Limited access to sub-specialty expertise |
| **Academic medical centers** | Research and education | Fragmented data across systems |
| **Heart failure programs** | GDMT optimization, transplant evaluation | Manual chart review for risk stratification |
| **Structural heart teams** | TAVR/MitraClip planning | Multi-modal data synthesis |
| **Cardio-oncology clinics** | Cardiotoxicity surveillance | No integrated monitoring platform |
| **Preventive cardiology** | Risk stratification, statin eligibility | Calculator fatigue, guideline complexity |
| **Clinical trial sites** | Patient screening, endpoint adjudication | Manual eligibility assessment |
| **Cardiovascular genomics** | Variant interpretation in cardiac context | Limited cardiac-specific annotation |

---

## 4. Existing HCLS AI Factory Architecture

### 4.1 Platform Overview

The HCLS AI Factory is a three-stage precision medicine platform running on NVIDIA DGX Spark:

```
Stage 1: Genomics Pipeline (Parabricks + DeepVariant)
    FASTQ → VCF → 3.56M annotated variants
         |
Stage 2: RAG/Chat Pipeline (Milvus + Claude)
    Variant interpretation, clinical significance
         |
Stage 3: Drug Discovery Pipeline (BioNeMo + DiffDock)
    Target → Lead compound → Docking → Drug-likeness
```

Five intelligence agents extend this core platform with domain-specific knowledge:

| Agent | Collections | Seed Vectors | Unique Capability |
|---|---|---|---|
| Precision Biomarker | 11 | 6,134 | Biological age calculators, biomarker panels |
| Precision Oncology | 10 | 609 | Molecular tumor board packets, trial matching |
| CAR-T Intelligence | 11 | 6,266 | CAR construct comparison, manufacturing optimization |
| Imaging Intelligence | 10 | 876 | NIM inference, DICOM workflows, 3D segmentation |
| Autoimmune Intelligence | 10 | ~500 | Autoantibody panels, flare prediction |

### 4.2 Shared Infrastructure

All agents share:

- **Milvus 2.4** vector database (IVF_FLAT, COSINE, 384-dim)
- **BGE-small-en-v1.5** embedding model (sentence-transformers)
- **Claude Sonnet 4.6** (Anthropic) primary LLM
- **`genomic_evidence`** collection (3,561,170 variants, read-only)
- **Docker Compose** orchestration
- **FastAPI** (REST) + **Streamlit** (UI) pattern
- **lib/hcls_common** shared library (23 modules)

### 4.3 Proven Patterns

The Cardiology Intelligence Agent will leverage battle-tested patterns from existing agents:

| Pattern | Proven In | Adaptation for Cardiology |
|---|---|---|
| Multi-collection parallel search | All 5 agents | 12 cardiology-specific collections |
| Knowledge graph augmentation | CAR-T, Biomarker | Cardiac conditions, drug classes, risk factors |
| Query expansion maps | CAR-T (12 maps), Biomarker | Cardiology terminology (e.g., "MI" → "myocardial infarction", "STEMI", "NSTEMI", "ACS") |
| Comparative analysis | CAR-T, Imaging | "TAVR vs SAVR", "Amiodarone vs Sotalol" |
| Cross-modal genomic triggers | Imaging (Lung-RADS → EGFR) | Imaging finding → cardiomyopathy gene panel |
| FHIR R4 export | Imaging | DiagnosticReport with cardiac SNOMED/LOINC codes |
| NIM inference workflows | Imaging (4 NIMs) | Cardiac-specific NIM models |
| Sidebar guided tour | Imaging | Cardiology demo flow |
| Pre-filled example queries | Imaging, Biomarker | Cardiology-specific starter questions |

---

## 5. Cardiology Intelligence Agent Architecture

### 5.1 System Diagram

```
+==========================================================================+
|  PRESENTATION:  Streamlit Cardiology Workbench (8536)                    |
|                 10 Tabs | Evidence | Workflows | Imaging | Risk Calcs    |
|                 FastAPI REST Server (8526)                                |
+==========================================================================+
                    |                            |
+==========================================================================+
|  INTELLIGENCE:   Cardiology RAG Engine                                   |
|                  12-collection parallel search                           |
|                  Knowledge graph (30 conditions, 20 biomarkers)          |
|                  Query expansion (15 maps, cardiology terminology)       |
|                  Comparative analysis ("TAVR vs SAVR")                   |
|                  Risk calculator engine (6 validated scores)             |
|                  GDMT optimization engine                                |
+==========================================================================+
                    |                            |
+==========================================================================+
|  INFERENCE:      NIM Services (VISTA-3D, MAISI, VILA-M3, Llama-3)       |
|                  8 Clinical Workflows:                                    |
|                  CAD | HF | Valve | Arrhythmia | CMR | Stress |         |
|                  Prevention | Cardio-Onc                                 |
+==========================================================================+
                    |                            |
+==========================================================================+
|  DATA:           Milvus 2.4 (12 cardiology collections + genomic)       |
|                  BGE-small-en-v1.5 (384-dim, IVF_FLAT, COSINE)          |
|                  PubMed, ClinicalTrials.gov, ACC/AHA guidelines          |
|                  Curated seed data (echo, CT, MRI, ECG, hemodynamics)    |
+==========================================================================+
```

### 5.2 Design Principles

1. **Guideline-first reasoning**: Every recommendation traces to ACC/AHA/ESC guideline evidence levels (Class I/IIa/IIb/III, LOE A/B/C)
2. **Quantitative precision**: Cardiac measurements retain units and reference ranges (LVEF 55-70%, E/e' < 14, QTc < 470ms)
3. **Cross-modal by default**: Every significant finding is checked against genomic context
4. **Risk-stratified output**: Severity badges and urgency routing aligned with clinical acuity
5. **Graceful degradation**: Full functionality in mock mode without GPU or live NIM services
6. **Familiar patterns**: Follows the same FastAPI + Streamlit + Milvus patterns as existing agents

### 5.3 Port Allocation

| Port | Service |
|---|---|
| 8526 | FastAPI REST Server |
| 8536 | Streamlit Cardiology Workbench |
| 19530 | Milvus (shared) |
| 8520 | NIM LLM (shared) |
| 8530 | NIM VISTA-3D (shared) |
| 8531 | NIM MAISI (shared) |
| 8532 | NIM VILA-M3 (shared) |

---

## 6. Milvus Collection Design

### 6.1 Index Configuration

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist / nprobe | 1024 / 16 |
| Dimension | 384 |
| Embedding model | BAAI/bge-small-en-v1.5 |

### 6.2 Collection Schemas

#### Collection 1: `cardio_literature` -- ~3,000 records

Published cardiovascular research papers, reviews, and meta-analyses.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | PubMed ID or unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| title | VARCHAR(500) | Paper title |
| text_chunk | VARCHAR(8000) | Abstract or text section |
| year | INT16 | Publication year |
| journal | VARCHAR(200) | Journal name (JACC, Circulation, EHJ, etc.) |
| cv_domain | VARCHAR(100) | Cardiovascular subdomain (imaging, EP, HF, valve, prevention) |
| modality | VARCHAR(50) | Imaging modality if applicable |
| study_type | VARCHAR(50) | RCT, meta-analysis, cohort, case-control, review |
| keywords | VARCHAR(500) | MeSH terms and author keywords |

**Source:** PubMed E-utilities with cardiovascular MeSH filters.

#### Collection 2: `cardio_trials` -- ~500 records

Cardiovascular clinical trials from ClinicalTrials.gov and landmark trial results.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | NCT number or trial identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| title | VARCHAR(500) | Official trial title |
| text_summary | VARCHAR(4000) | Trial summary including results |
| phase | VARCHAR(20) | Phase I-IV |
| status | VARCHAR(30) | Active, completed, recruiting |
| sponsor | VARCHAR(200) | Lead sponsor |
| cv_domain | VARCHAR(100) | CAD, HF, valve, arrhythmia, prevention |
| intervention | VARCHAR(300) | Drug, device, or procedure tested |
| primary_endpoint | VARCHAR(300) | Primary outcome measure |
| enrollment | INT32 | Number of participants |
| start_year | INT16 | Year trial began |
| outcome_summary | VARCHAR(2000) | Key results (if completed) |
| landmark | BOOL | Is this a landmark trial (PARADIGM-HF, DAPA-HF, etc.) |

**Source:** ClinicalTrials.gov V2 API with cardiovascular condition filters.

#### Collection 3: `cardio_imaging` -- ~200 records

Cardiac imaging protocols, findings, and measurements across all modalities.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Finding or protocol description |
| modality | VARCHAR(50) | echo, cardiac_ct, cardiac_mri, nuclear, cath |
| imaging_type | VARCHAR(100) | TTE, TEE, stress_echo, CTA, calcium_score, etc. |
| finding_category | VARCHAR(100) | Chamber size, systolic function, diastolic function, valve, plaque, LGE, perfusion |
| measurement_name | VARCHAR(100) | LVEF, GLS, E/e', calcium_score, stenosis_pct |
| measurement_value | VARCHAR(50) | Numeric value with units |
| reference_range | VARCHAR(100) | Normal range per guidelines |
| severity | VARCHAR(20) | Normal, mild, moderate, severe |
| guideline_source | VARCHAR(100) | ASE, SCCT, SCMR, ASNC |
| clinical_significance | VARCHAR(500) | Interpretation guidance |

**Source:** Curated from ASE/SCCT/SCMR/ASNC guideline documents and reference texts.

#### Collection 4: `cardio_electrophysiology` -- ~150 records

ECG interpretation, arrhythmia classification, and electrophysiology data.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | ECG pattern or arrhythmia description |
| ecg_pattern | VARCHAR(100) | LBBB, RBBB, LVH, ST_elevation, AF, VT, etc. |
| rhythm | VARCHAR(50) | Sinus, AF, flutter, SVT, VT, VF |
| rate_category | VARCHAR(30) | Bradycardia, normal, tachycardia |
| interval | VARCHAR(50) | PR, QRS, QT, QTc |
| interval_value_ms | FLOAT | Numeric interval value |
| clinical_significance | VARCHAR(500) | What this finding means |
| urgency | VARCHAR(20) | Routine, urgent, emergent |
| differential_diagnosis | VARCHAR(500) | DDx list |
| management | VARCHAR(500) | Recommended next steps |

**Source:** Curated from ACC/AHA/HRS guidelines and electrophysiology references.

#### Collection 5: `cardio_heart_failure` -- ~150 records

Heart failure classification, GDMT protocols, and management algorithms.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | HF management recommendation |
| hf_type | VARCHAR(30) | HFrEF, HFmrEF, HFpEF |
| acc_aha_stage | VARCHAR(5) | A, B, C, D |
| nyha_class | VARCHAR(5) | I, II, III, IV |
| drug_class | VARCHAR(100) | Beta-blocker, ARNI, MRA, SGLT2i, hydralazine/nitrate |
| drug_name | VARCHAR(100) | Specific medication |
| target_dose | VARCHAR(50) | Goal dose per guidelines |
| titration_protocol | VARCHAR(500) | How to uptitrate |
| evidence_level | VARCHAR(20) | Class I/IIa/IIb/III, LOE A/B/C |
| landmark_trial | VARCHAR(100) | Supporting trial name |
| contraindications | VARCHAR(500) | When not to use |

**Source:** ACC/AHA/HFSA 2022 Heart Failure Guidelines, ESC 2023 update.

#### Collection 6: `cardio_valvular` -- ~120 records

Valvular heart disease assessment, severity grading, and intervention criteria.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Valve assessment or recommendation |
| valve | VARCHAR(30) | Aortic, mitral, tricuspid, pulmonic |
| pathology | VARCHAR(50) | Stenosis, regurgitation, prolapse, endocarditis |
| severity | VARCHAR(20) | Mild, moderate, severe |
| quantitative_criteria | VARCHAR(500) | EOA, gradient, regurgitant volume, EROA, vena contracta |
| intervention_criteria | VARCHAR(500) | When to intervene per guidelines |
| intervention_type | VARCHAR(100) | SAVR, TAVR, mitral_repair, MitraClip, TMVr |
| sts_risk_threshold | VARCHAR(50) | Low (<4%), intermediate (4-8%), high (>8%) |
| evidence_level | VARCHAR(20) | Class/LOE |
| surveillance_protocol | VARCHAR(500) | Follow-up echo frequency |

**Source:** ACC/AHA 2020 Valvular Heart Disease Guidelines, ESC 2021 VHD Guidelines.

#### Collection 7: `cardio_prevention` -- ~150 records

Preventive cardiology: risk stratification, lipid management, and lifestyle intervention.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Prevention recommendation |
| risk_category | VARCHAR(50) | Low, borderline, intermediate, high, very_high |
| risk_calculator | VARCHAR(50) | PCE/ASCVD, MESA, Framingham, Reynolds |
| biomarker | VARCHAR(100) | LDL-C, Lp(a), ApoB, hsCRP, CAC |
| target_value | VARCHAR(50) | LDL < 70, LDL < 55 (very high risk) |
| therapy_class | VARCHAR(100) | Statin, ezetimibe, PCSK9i, bempedoic acid, inclisiran |
| therapy_name | VARCHAR(100) | Specific medication |
| evidence_level | VARCHAR(20) | Class/LOE |
| guideline_source | VARCHAR(100) | ACC/AHA 2018 Cholesterol, ESC 2019 Dyslipidemia |
| lifestyle_intervention | VARCHAR(500) | Diet, exercise, smoking, weight management |

**Source:** ACC/AHA 2019 Prevention Guidelines, 2018 Cholesterol Guidelines.

#### Collection 8: `cardio_interventional` -- ~100 records

Interventional cardiology procedures, techniques, and outcomes data.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Procedure description or outcome data |
| procedure_type | VARCHAR(100) | PCI, TAVR, MitraClip, WATCHMAN, PFO_closure, CTO_PCI |
| indication | VARCHAR(200) | Clinical indication |
| technique | VARCHAR(300) | Specific approach (radial, femoral, antegrade, retrograde) |
| device | VARCHAR(200) | Stent, valve, clip, occluder type |
| success_rate | VARCHAR(50) | Procedural success rate |
| complication_rate | VARCHAR(50) | Major adverse event rate |
| antithrombotic_protocol | VARCHAR(500) | Antiplatelet/anticoagulant regimen |
| evidence_level | VARCHAR(20) | Class/LOE |
| landmark_trial | VARCHAR(100) | Supporting trial |

**Source:** SCAI guidelines, ACC/AHA PCI and VHD guidelines.

#### Collection 9: `cardio_oncology` -- ~100 records

Cardio-oncology surveillance, cardiotoxicity detection, and cardioprotection.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Cardio-oncology recommendation |
| cancer_therapy | VARCHAR(100) | Anthracycline, trastuzumab, ICI, TKI, radiation |
| cardiotoxicity_type | VARCHAR(100) | CTRCD, myocarditis, QTc_prolongation, HTN, VTE |
| monitoring_protocol | VARCHAR(500) | Echo frequency, GLS thresholds, troponin schedule |
| gls_threshold | VARCHAR(50) | >15% relative decline = subclinical toxicity |
| risk_score | VARCHAR(50) | HFA-ICOS baseline risk |
| cardioprotective_strategy | VARCHAR(500) | Dexrazoxane, beta-blocker, ACEi, statin |
| when_to_hold | VARCHAR(500) | Criteria to hold/discontinue cancer therapy |
| evidence_level | VARCHAR(20) | Class/LOE |
| guideline_source | VARCHAR(100) | ESC 2022 Cardio-Oncology, ASCO |

**Source:** ESC 2022 Cardio-Oncology Guidelines, ASCO Cardiovascular Toxicity Guidelines.

#### Collection 10: `cardio_devices` -- ~80 records

FDA-cleared cardiovascular AI devices and implantable cardiac devices.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Device description and capabilities |
| device_name | VARCHAR(200) | Commercial name |
| manufacturer | VARCHAR(200) | Company |
| device_category | VARCHAR(50) | AI_software, implantable, wearable |
| modality | VARCHAR(50) | ECG, echo, CT, MRI, wearable |
| clinical_task | VARCHAR(100) | Detection, measurement, prediction |
| regulatory_status | VARCHAR(50) | FDA_cleared, CE_marked, investigational |
| clearance_date | VARCHAR(20) | Approval date |
| performance_summary | VARCHAR(500) | Sensitivity, specificity, AUC |
| intended_use | VARCHAR(500) | FDA-labeled intended use statement |

**Source:** FDA AI/ML-Enabled Medical Device database, company filings.

#### Collection 11: `cardio_guidelines` -- ~150 records

ACC/AHA/ESC/HRS clinical practice guidelines and focused updates.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Guideline recommendation |
| guideline_name | VARCHAR(300) | Full guideline title |
| organization | VARCHAR(50) | ACC/AHA, ESC, HRS, SCAI, ASE |
| year | INT16 | Publication year |
| cv_domain | VARCHAR(100) | CAD, HF, valve, arrhythmia, prevention |
| recommendation_class | VARCHAR(10) | I, IIa, IIb, III |
| level_of_evidence | VARCHAR(5) | A, B-R, B-NR, C-LD, C-EO |
| key_recommendation | VARCHAR(2000) | Specific recommendation text |
| clinical_scenario | VARCHAR(500) | When this recommendation applies |
| related_guidelines | VARCHAR(300) | Cross-references to other guidelines |

**Source:** ACC/AHA Practice Guidelines library, ESC Clinical Practice Guidelines.

#### Collection 12: `cardio_hemodynamics` -- ~80 records

Hemodynamic data: catheterization, pressure tracings, and derived calculations.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Hemodynamic scenario or reference |
| measurement_type | VARCHAR(100) | PCWP, PA_pressure, CO, SVR, PVR, LVEDP |
| normal_range | VARCHAR(50) | Expected normal values |
| abnormal_pattern | VARCHAR(200) | What elevation/depression indicates |
| associated_conditions | VARCHAR(500) | Clinical conditions causing this pattern |
| calculation_formula | VARCHAR(300) | How to derive (e.g., SVR = 80*(MAP-RAP)/CO) |
| clinical_significance | VARCHAR(500) | Interpretation guidance |
| catheterization_type | VARCHAR(50) | Right heart, left heart, coronary |

**Source:** Curated from hemodynamic references and catheterization guidelines.

#### Collection 13: `genomic_evidence` -- 3,561,170 records (read-only)

Shared genomic variant collection from HCLS AI Factory Stage 1+2.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Variant identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Variant annotation |
| gene | VARCHAR(50) | Gene symbol |
| clinical_significance | VARCHAR(50) | Pathogenic, likely pathogenic, VUS, benign |
| disease_associations | VARCHAR(1000) | Associated diseases |
| am_pathogenicity | FLOAT | AlphaMissense pathogenicity score |

**Purpose:** Cross-modal triggers query this collection for cardiovascular-relevant genes (LDLR, MYH7, TTN, SCN5A, etc.).

### 6.3 Collection Search Weights

| Collection | Weight | Rationale |
|---|---|---|
| Literature | 0.16 | Largest corpus, broadest evidence |
| Guidelines | 0.14 | Highest clinical authority |
| Trials | 0.12 | Primary evidence source |
| Imaging | 0.10 | Multi-modal imaging data |
| Heart Failure | 0.08 | Common, high-impact domain |
| Electrophysiology | 0.08 | Arrhythmia management |
| Valvular | 0.07 | Structural heart disease |
| Prevention | 0.07 | Primary/secondary prevention |
| Hemodynamics | 0.05 | Invasive assessment |
| Interventional | 0.05 | Procedural decision-making |
| Cardio-Oncology | 0.04 | Growing subspecialty |
| Devices | 0.04 | Technology landscape |

### 6.4 Estimated Vector Counts

| Collection | Seed Records | Post-Ingest Target |
|---|---|---|
| cardio_literature | 200 | 3,000+ (PubMed ingest) |
| cardio_trials | 50 | 500+ (ClinicalTrials.gov ingest) |
| cardio_imaging | 200 | 200 |
| cardio_electrophysiology | 150 | 150 |
| cardio_heart_failure | 150 | 150 |
| cardio_valvular | 120 | 120 |
| cardio_prevention | 150 | 150 |
| cardio_interventional | 100 | 100 |
| cardio_oncology | 100 | 100 |
| cardio_devices | 80 | 80 |
| cardio_guidelines | 150 | 150 |
| cardio_hemodynamics | 80 | 80 |
| **Total (owned)** | **~1,530** | **~4,780+** |
| genomic_evidence (read-only) | -- | 3,561,170 |

---

## 7. Clinical Workflows

### 7.1 Workflow Architecture

All workflows follow the established `BaseImagingWorkflow` pattern: `preprocess → infer → postprocess → WorkflowResult`. Each workflow supports full mock mode with clinically realistic synthetic results.

### 7.2 Eight Reference Workflows

#### Workflow 1: Coronary Artery Disease Assessment

| Attribute | Value |
|---|---|
| Workflow ID | `coronary_artery_disease` |
| Input | Coronary CTA, calcium score CT, or catheterization data |
| Target Latency | < 3 minutes |
| Models | VISTA-3D (coronary segmentation), vessel analysis pipeline |
| Key Outputs | CAD-RADS classification (0-5), Agatston calcium score, per-vessel stenosis grading, plaque characterization, FFR-CT estimate |
| Severity Routing | CAD-RADS 4-5 → Urgent cardiology consult |
| Cross-Modal Trigger | CAD-RADS 4+ with age < 55 → FH gene panel (LDLR, PCSK9, APOB) |
| Guideline Alignment | ACC/AHA 2021 Chest Pain Guidelines, SCCT CAD-RADS 2.0 |

**Clinical Decision Logic:**

```
Calcium Score:
  0        → Very low risk, no further imaging needed
  1-99     → Low risk, consider statin, lifestyle
  100-399  → Moderate risk, statin indicated, stress test if symptoms
  ≥400     → High risk, coronary CTA or functional testing

CAD-RADS:
  0        → No plaque, no stenosis
  1        → 1-24% stenosis, minimal
  2        → 25-49%, mild
  3        → 50-69%, moderate → consider functional testing
  4A       → 70-99% → catheterization or FFR-CT
  4B       → Left main ≥ 50% or 3-vessel ≥ 70% → surgical evaluation
  5        → Total occlusion → viability assessment

Plaque features → High-risk if:
  - Low-attenuation plaque (< 30 HU)
  - Positive remodeling (RI > 1.1)
  - Napkin-ring sign
  - Spotty calcification
```

#### Workflow 2: Heart Failure Classification and GDMT Optimization

| Attribute | Value |
|---|---|
| Workflow ID | `heart_failure_gdmt` |
| Input | Echo LVEF, BNP/NT-proBNP, clinical parameters |
| Target Latency | < 30 seconds |
| Key Outputs | HF classification (HFrEF/HFmrEF/HFpEF), ACC/AHA stage, NYHA class, GDMT recommendations with target doses, MAGGIC score |
| Severity Routing | Stage D or NYHA IV → Advanced HF evaluation |
| Cross-Modal Trigger | Age < 50 + non-ischemic DCM → cardiomyopathy gene panel (TTN, LMNA, MYH7, MYBPC3) |
| Guideline Alignment | ACC/AHA 2022 HF Guidelines, HFSA 2024 Update |

**GDMT Optimization Algorithm:**

```
HFrEF (LVEF ≤ 40%):
  Pillar 1: Beta-blocker (carvedilol, metoprolol succinate, bisoprolol)
  Pillar 2: ARNI (sacubitril/valsartan) or ACEi/ARB
  Pillar 3: MRA (spironolactone, eplerenone)
  Pillar 4: SGLT2i (dapagliflozin, empagliflozin)

  → If African American + NYHA III-IV: Add hydralazine/isosorbide dinitrate
  → If LVEF ≤ 35% + NYHA II-III + QRS ≥ 150ms LBBB: CRT
  → If LVEF ≤ 35% + ≥ 40 days post-MI or ≥ 90 days on optimal GDMT: ICD
  → If Stage D refractory: LVAD or transplant evaluation

HFpEF (LVEF ≥ 50%):
  SGLT2i (Class I, LOE A -- EMPEROR-Preserved, DELIVER)
  Diuretics for congestion
  Treat comorbidities (HTN, AF, obesity, sleep apnea, CAD)
  GLP-1 RA if obesity + HFpEF (STEP-HFpEF)
```

#### Workflow 3: Valvular Heart Disease Quantification

| Attribute | Value |
|---|---|
| Workflow ID | `valvular_disease` |
| Input | Echocardiographic measurements, clinical parameters |
| Target Latency | < 30 seconds |
| Key Outputs | Valve severity grading, intervention threshold assessment, SAVR vs TAVR recommendation for AS, STS-PROM risk score |
| Severity Routing | Severe symptomatic AS or MR → Heart team evaluation |
| Cross-Modal Trigger | Bicuspid AV + aortic dilation → aortopathy gene panel (FBN1, TGFBR1/2, SMAD3, ACTA2) |
| Guideline Alignment | ACC/AHA 2020 VHD Guidelines, ESC 2021 VHD Guidelines |

**Severity Grading Criteria (Aortic Stenosis):**

| Parameter | Mild | Moderate | Severe |
|---|---|---|---|
| Peak velocity (m/s) | 2.0-2.9 | 3.0-3.9 | ≥ 4.0 |
| Mean gradient (mmHg) | < 20 | 20-39 | ≥ 40 |
| AVA (cm²) | > 1.5 | 1.0-1.5 | < 1.0 |
| Indexed AVA (cm²/m²) | > 0.85 | 0.60-0.85 | < 0.60 |

#### Workflow 4: Arrhythmia Detection and Management

| Attribute | Value |
|---|---|
| Workflow ID | `arrhythmia_management` |
| Input | ECG data, clinical parameters, cardiac imaging |
| Target Latency | < 15 seconds |
| Key Outputs | Rhythm classification, CHA₂DS₂-VASc score (if AF), HAS-BLED score, rate vs rhythm control recommendation, ablation candidacy |
| Severity Routing | VT/VF or complete heart block → Emergent |
| Cross-Modal Trigger | Unexplained VT in young patient → channelopathy panel (SCN5A, KCNQ1, KCNH2, RYR2) |
| Guideline Alignment | ACC/AHA/HRS 2023 AF Guidelines, 2017 VA/SCD Guidelines |

**CHA₂DS₂-VASc Calculator:**

```
C  Congestive HF (or LVEF ≤ 40%)     +1
H  Hypertension                        +1
A₂ Age ≥ 75                            +2
D  Diabetes mellitus                   +1
S₂ Stroke/TIA/thromboembolism          +2
V  Vascular disease (PAD, prior MI)    +1
A  Age 65-74                           +1
Sc Sex category (female)               +1

Score 0 (males) or 1 (females): No anticoagulation
Score 1 (males): Consider anticoagulation
Score ≥ 2: Anticoagulation recommended (Class I)
```

#### Workflow 5: Cardiac MRI Tissue Characterization

| Attribute | Value |
|---|---|
| Workflow ID | `cardiac_mri_tissue` |
| Input | Cardiac MRI data (cine, T1/T2 mapping, LGE, perfusion) |
| Target Latency | < 5 minutes |
| Models | VISTA-3D (cardiac segmentation), parametric map analysis |
| Key Outputs | Biventricular volumes and function, LGE pattern (ischemic vs non-ischemic), T1/T2/ECV values, perfusion defect quantification, tissue diagnosis |
| Severity Routing | Active myocarditis or severe biventricular failure → Urgent |
| Cross-Modal Trigger | Non-ischemic LGE + LV dysfunction → cardiomyopathy gene panel |
| Guideline Alignment | SCMR 2020 Standardized Imaging Protocols, Lake Louise Criteria |

**LGE Pattern Interpretation:**

| Pattern | Location | Diagnosis |
|---|---|---|
| Subendocardial / transmural | Coronary territory | Ischemic cardiomyopathy (prior MI) |
| Mid-wall (septal) | Interventricular septum | Dilated cardiomyopathy, myocarditis |
| Epicardial / patchy | Lateral wall | Myocarditis (acute/chronic) |
| Diffuse subendocardial | Global | Cardiac amyloidosis |
| RV insertion point | Anterior/posterior | Pulmonary hypertension, HCM |
| Asymmetric septal | Basal septum | Hypertrophic cardiomyopathy |
| Basal inferolateral | Mid-wall | Cardiac sarcoidosis, Chagas |
| No LGE + elevated T1/ECV | Global | Early-stage amyloidosis or diffuse fibrosis |

#### Workflow 6: Stress Test Interpretation

| Attribute | Value |
|---|---|
| Workflow ID | `stress_test_interpretation` |
| Input | Exercise or pharmacologic stress test data (ECG, imaging) |
| Target Latency | < 30 seconds |
| Key Outputs | Duke Treadmill Score, ischemia assessment, functional capacity (METs), hemodynamic response, stress imaging interpretation |
| Severity Routing | High-risk DTS or large perfusion defect → Catheterization |
| Guideline Alignment | ACC/AHA 2021 Chest Pain Guidelines, ASNC Stress Testing Standards |

**Duke Treadmill Score:**

```
DTS = Exercise time (minutes) - (5 × max ST deviation mm) - (4 × angina index)
  Angina index: 0 = none, 1 = non-limiting, 2 = exercise-limiting

Low risk:     DTS ≥ +5    (annual mortality < 1%)
Moderate risk: DTS -10 to +4 (annual mortality 1-3%)
High risk:    DTS < -10   (annual mortality ≥ 5%)
```

#### Workflow 7: Preventive Risk Stratification

| Attribute | Value |
|---|---|
| Workflow ID | `preventive_risk_stratification` |
| Input | Demographics, lipids, BP, diabetes status, smoking, family history |
| Target Latency | < 15 seconds |
| Key Outputs | 10-year ASCVD risk, risk category, statin benefit group, statin intensity recommendation, risk-enhancing factors, CAC-based reclassification guidance |
| Cross-Modal Trigger | 10-year ASCVD ≥ 7.5% + family history of premature CAD → FH gene panel (LDLR, PCSK9, APOB) |
| Guideline Alignment | ACC/AHA 2019 Primary Prevention, 2018 Cholesterol Guidelines |

**Statin Decision Framework:**

```
1. Clinical ASCVD? → High-intensity statin (Class I)
   - If very high risk (recent ACS, multiple events): Add ezetimibe → if LDL still ≥ 55: Add PCSK9i

2. LDL ≥ 190 mg/dL? → High-intensity statin without risk calculation (Class I)
   - Consider FH genetic testing

3. Diabetes, age 40-75? → Moderate-intensity statin (Class I)
   - If multiple risk factors or 10y ASCVD ≥ 20%: High-intensity

4. 10-year ASCVD 7.5-19.9%? → Moderate-intensity statin (Class I)
   - If risk decision uncertain: Obtain CAC score
   - CAC = 0: Withhold statin (unless diabetes, FH, smoking)
   - CAC 1-99: Favors statin
   - CAC ≥ 100: Statin indicated

5. 10-year ASCVD 5-7.4%? → Consider if risk-enhancing factors present
   Risk-enhancing factors: Family hx premature ASCVD, Lp(a) ≥ 50 mg/dL,
   hsCRP ≥ 2.0, ABI < 0.9, metabolic syndrome, South Asian ancestry
```

#### Workflow 8: Cardio-Oncology Surveillance

| Attribute | Value |
|---|---|
| Workflow ID | `cardio_oncology_surveillance` |
| Input | Cancer therapy regimen, baseline cardiac function, risk factors |
| Target Latency | < 30 seconds |
| Key Outputs | HFA-ICOS baseline risk category, surveillance protocol (echo frequency, biomarker schedule), GLS monitoring thresholds, cardioprotective recommendations, criteria for therapy modification |
| Severity Routing | LVEF decline > 10% to < 50% or GLS decline > 15% → Urgent cardio-oncology consult |
| Guideline Alignment | ESC 2022 Cardio-Oncology Guidelines, ASCO CV Toxicity Guidelines |

**Cardiotoxicity Monitoring Protocol:**

```
Anthracycline-based chemotherapy:
  Baseline: Echo (LVEF + GLS) + troponin + NT-proBNP
  During: Echo every 2 cycles (or after cumulative dose milestones)
  After: Echo at 3, 6, 12 months post-completion, then annually
  GLS threshold: >15% relative decline from baseline = subclinical cardiotoxicity
  LVEF threshold: Decline >10% to below 50% = CTRCD

HER2-targeted therapy (trastuzumab):
  Echo every 3 months during therapy
  If LVEF 40-49%: Continue with closer monitoring + cardioprotection
  If LVEF < 40%: Hold therapy, refer cardio-oncology

Immune checkpoint inhibitors:
  Baseline ECG + troponin
  Troponin every cycle for first 4 cycles
  If troponin rise: Urgent CMR to rule out myocarditis
  ICI myocarditis mortality: 25-50% -- requires immediate therapy hold + high-dose steroids
```

---

## 8. Cross-Modal Integration

### 8.1 Cardiovascular Genomics Triggers

The Cardiology Intelligence Agent implements cross-modal triggers that automatically query the shared `genomic_evidence` collection (3.5M variants) when imaging or clinical findings suggest a genetic etiology:

| Trigger Condition | Gene Panel Queried | Clinical Rationale |
|---|---|---|
| Non-ischemic DCM, age < 50 | TTN, LMNA, RBM20, MYH7, MYBPC3, TNNT2, SCN5A, BAG3, PLN, FLNC, DSP | 25-40% of DCM has genetic basis; LMNA carriers need ICD regardless of LVEF |
| Unexplained LVH, septal hypertrophy | MYH7, MYBPC3, TNNT2, TNNI3, TPM1, ACTC1, MYL2, MYL3, GLA (Fabry) | HCM is most common inherited cardiac disease (1:500) |
| Unexplained VT or cardiac arrest | SCN5A, KCNQ1, KCNH2, RYR2, CASQ2, PKP2, DSP, DSG2, DSC2 | Channelopathy or ARVC screen; cascade screening for families |
| Premature CAD (< 55M / < 65F) or LDL > 190 | LDLR, PCSK9, APOB, LDLRAP1, LIPA | FH prevalence 1:250; early detection enables prevention |
| Aortic root dilation > 4.5 cm | FBN1, TGFBR1, TGFBR2, SMAD3, ACTA2, MYH11, PRKG1, LOX | Heritable thoracic aortic disease; surgical thresholds differ by genotype |
| Cardiac amyloid (LGE + elevated T1) | TTR (transthyretin), AL amyloid genes | TTR amyloidosis is treatable with tafamidis; gene-specific therapy |
| Bicuspid aortic valve + coarctation | NOTCH1, GATA4/5/6, SMAD6 | Congenital heart disease genetics |

### 8.2 Imaging ↔ Genomics ↔ Drug Discovery Pipeline

```
Cardiac Imaging Finding (Echo, CT, MRI)
    |
    v
[Cross-Modal Trigger] -- Clinical criteria met?
    |                         |
    YES                       NO
    |                         |
    v                         v
[Query genomic_evidence]   [Standard RAG response]
(3.5M variant vectors)
    |
    v
[Variant Annotation]
ClinVar pathogenicity, AlphaMissense score
    |
    v
[HCLS AI Factory Stage 2: RAG Target ID]
Identify druggable targets from genetic findings
    |
    v
[HCLS AI Factory Stage 3: Drug Discovery]
BioNeMo molecule generation, DiffDock binding
    |
    v
[Clinical Output]
FHIR R4 DiagnosticReport with genomic enrichment
```

### 8.3 Integration with Other Agents

| Integration | Direction | Data Flow |
|---|---|---|
| Imaging Agent | Bidirectional | Shares cardiac CT/MRI workflow results; receives DICOM routing for cardiac studies |
| Precision Biomarker Agent | Inbound | Receives cardiac biomarker reference ranges and age-stratified norms |
| Precision Oncology Agent | Inbound | Receives cancer therapy regimens for cardio-oncology surveillance |
| CAR-T Intelligence Agent | Inbound | Receives CRS-related cardiac toxicity data |
| Genomics Pipeline | Read-only | Queries `genomic_evidence` for cardiovascular gene variants |
| Drug Discovery Pipeline | Outbound | Sends confirmed genetic targets for compound screening |

---

## 9. NIM Integration Strategy

### 9.1 Shared NIM Services

The Cardiology Agent will reuse the four NIM services already deployed by the Imaging Intelligence Agent:

| NIM | Port | Cardiology Application |
|---|---|---|
| VISTA-3D | 8530 | Cardiac chamber segmentation, coronary artery segmentation, pericardial fat quantification |
| MAISI | 8531 | Synthetic cardiac CT generation for training, rare pathology simulation |
| VILA-M3 | 8532 | Cardiac image interpretation, echo measurement validation, ECG pattern recognition |
| Llama-3 8B | 8520 | Clinical reasoning fallback when Claude API unavailable |

### 9.2 VISTA-3D Cardiac Applications

VISTA-3D's 132 anatomical classes include cardiac structures suitable for:

- **Chamber volumetrics**: LV, RV, LA, RA segmentation for volume and EF calculation
- **Myocardial segmentation**: LV myocardium for mass calculation and wall motion assessment
- **Pericardial analysis**: Pericardial effusion detection and quantification
- **Great vessel assessment**: Aortic root and ascending aorta measurement
- **Coronary segmentation**: Combined with vessel centerline analysis for stenosis grading

### 9.3 Fallback Logic

```
Request → NIM available? → Yes → Real NIM inference
                        → No  → Cloud NIM? → Yes → NVIDIA Cloud API
                                            → No  → Mock enabled? → Yes → Synthetic result
                                                                   → No  → Error
```

---

## 10. Knowledge Graph Design

### 10.1 Graph Structure

The cardiology knowledge graph will contain structured entries across six dictionaries:

| Dictionary | Entries | Content |
|---|---|---|
| Cardiovascular Conditions | ~35 | ICD-10 codes, diagnostic criteria, severity classification, imaging characteristics, guideline-based management |
| Cardiac Biomarkers | ~25 | Reference ranges, clinical cutoffs, kinetics, diagnostic/prognostic value |
| Drug Classes | ~30 | Mechanism, target dose, titration, contraindications, landmark trials |
| Risk Calculators | ~8 | Input variables, formula, risk categories, guideline recommendations per category |
| Imaging Protocols | ~20 | Modality, protocol parameters, normal values, abnormal patterns |
| Cardiovascular Genes | ~50 | Gene symbol, chromosome, associated conditions, inheritance, prevalence, clinical actionability |

### 10.2 Example Knowledge Graph Entries

**Condition: Hypertrophic Cardiomyopathy (HCM)**

```python
{
    "name": "Hypertrophic Cardiomyopathy",
    "icd10": "I42.1",
    "aliases": ["HCM", "HOCM", "IHSS", "hypertrophic obstructive cardiomyopathy"],
    "prevalence": "1:500 (most common inherited cardiac disease)",
    "inheritance": "Autosomal dominant, variable penetrance",
    "diagnostic_criteria": {
        "wall_thickness": "≥ 15 mm in any segment (or ≥ 13 mm with family history)",
        "lvot_gradient": "≥ 30 mmHg at rest or provocation = obstructive",
        "mri_lge": "Patchy mid-wall, RV insertion points, extensive LGE = worse prognosis"
    },
    "risk_stratification": {
        "model": "ACC/AHA 2024 HCM Guidelines SCD Risk Calculator",
        "high_risk_features": [
            "Family history of SCD from HCM",
            "Unexplained syncope",
            "Massive LVH (≥ 30 mm)",
            "NSVT on Holter",
            "Abnormal BP response to exercise",
            "Extensive LGE (≥ 15% LV mass)",
            "LVEF < 50%",
            "Apical aneurysm"
        ]
    },
    "genes": ["MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1", "MYL2", "MYL3"],
    "treatment": {
        "medical": "Beta-blockers (first-line), verapamil, disopyramide, mavacamten",
        "interventional": "Septal myectomy (surgical) or alcohol septal ablation",
        "device": "ICD if high SCD risk"
    },
    "imaging_modalities": ["echo", "cardiac_mri", "cardiac_ct"],
    "cross_modal_trigger": True
}
```

**Biomarker: NT-proBNP**

```python
{
    "name": "NT-proBNP",
    "full_name": "N-terminal pro-B-type natriuretic peptide",
    "loinc_code": "33762-6",
    "reference_ranges": {
        "age_lt_50": {"normal": "< 300 pg/mL", "grey_zone": "300-450", "elevated": "> 450"},
        "age_50_75": {"normal": "< 300 pg/mL", "grey_zone": "300-900", "elevated": "> 900"},
        "age_gt_75": {"normal": "< 300 pg/mL", "grey_zone": "300-1800", "elevated": "> 1800"}
    },
    "clinical_use": "Heart failure diagnosis and prognosis, GDMT titration monitoring",
    "kinetics": "Half-life ~120 min, less affected by obesity than BNP, not cleared by neprilysin",
    "confounders": ["Renal dysfunction (↑)", "Obesity (↓)", "AF (↑)", "Pulmonary disease (↑)"],
    "guideline_thresholds": {
        "rule_out_hf": "< 300 pg/mL (NPV > 98%)",
        "hospitalized_hf": "Admission > 1000 associated with higher mortality",
        "gdmt_success": "> 30% reduction on optimal therapy"
    }
}
```

---

## 11. Query Expansion and Retrieval Strategy

### 11.1 Cardiology-Specific Query Expansion Maps

Fifteen domain-specific expansion maps will map cardiovascular terminology to related terms:

| Map | Keywords → Terms | Example |
|---|---|---|
| Coronary Disease | 25 → 180 | "heart attack" → myocardial infarction, MI, STEMI, NSTEMI, ACS, acute coronary syndrome, troponin elevation |
| Heart Failure | 20 → 150 | "weak heart" → heart failure, HFrEF, HFpEF, reduced ejection fraction, systolic dysfunction, LVEF |
| Arrhythmia | 20 → 140 | "irregular heartbeat" → atrial fibrillation, AF, AFib, arrhythmia, flutter, SVT, palpitations |
| Valvular | 15 → 100 | "heart valve" → aortic stenosis, AS, mitral regurgitation, MR, TAVR, valve replacement |
| Lipids | 15 → 100 | "cholesterol" → LDL-C, HDL-C, triglycerides, statin, PCSK9, hyperlipidemia, dyslipidemia |
| Imaging Echo | 15 → 90 | "heart ultrasound" → echocardiogram, echo, TTE, TEE, LVEF, GLS, diastolic function |
| Imaging CT | 10 → 70 | "heart scan" → cardiac CT, CTA, calcium score, CAD-RADS, coronary CTA |
| Imaging MRI | 10 → 70 | "cardiac MRI" → CMR, LGE, T1 mapping, ECV, tissue characterization, parametric mapping |
| Electrophysiology | 15 → 100 | "ECG" → electrocardiogram, 12-lead, rhythm strip, QTc, ST segment, LBBB |
| Heart Failure Drugs | 15 → 90 | "entresto" → sacubitril/valsartan, ARNI, neprilysin inhibitor, PARADIGM-HF |
| Structural | 10 → 60 | "TAVR" → transcatheter aortic valve replacement, TAVI, Edwards SAPIEN, Medtronic CoreValve |
| Prevention | 10 → 70 | "heart risk" → ASCVD, cardiovascular risk, PCE, risk calculator, CAC score |
| Cardio-Oncology | 10 → 60 | "chemo heart" → cardiotoxicity, CTRCD, anthracycline, trastuzumab, GLS monitoring |
| Hemodynamics | 10 → 60 | "heart pressures" → right heart cath, Swan-Ganz, PCWP, wedge pressure, PVR, cardiac output |
| Devices | 10 → 60 | "pacemaker" → ICD, CRT, CRT-D, cardiac resynchronization, LVAD, defibrillator |

**Total: 200 keywords → ~1,400 expanded terms**

### 11.2 Comparative Analysis Detection

The RAG engine auto-detects comparative queries using the same pattern established in the CAR-T and Imaging agents:

**Trigger patterns:** "vs", "versus", "compared to", "compare", "difference between", "better than"

**Cardiology-specific comparisons:**

| Comparison | Clinical Relevance |
|---|---|
| "TAVR vs SAVR" | Aortic stenosis intervention strategy by risk category |
| "Amiodarone vs Sotalol" | AF rate/rhythm control |
| "DOAC vs Warfarin" | Anticoagulation strategy for AF |
| "PCI vs CABG" | Revascularization strategy for multivessel CAD |
| "CT vs MRI for cardiomyopathy" | Imaging modality selection |
| "Metoprolol vs Carvedilol" | Beta-blocker selection in HF |
| "Entresto vs Enalapril" | RAAS inhibition strategy in HFrEF |
| "Ablation vs Drugs for AF" | Rhythm control strategy |

---

## 12. API and UI Design

### 12.1 FastAPI Endpoints (Port 8526)

| Method | Path | Purpose |
|---|---|---|
| GET | `/health` | Service health, collection stats, NIM status |
| GET | `/collections` | List all collections with vector counts |
| POST | `/query` | Full RAG query with evidence synthesis |
| POST | `/search` | Evidence-only search (no LLM) |
| POST | `/api/ask` | Chat-style question answering |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/workflows` | List available clinical workflows |
| POST | `/workflow/{name}/run` | Execute a clinical workflow |
| GET | `/demo-cases` | List pre-loaded demo cases |
| POST | `/demo-cases/{id}/run` | Run a demo case |
| POST | `/risk/calculate` | Calculate validated risk score |
| POST | `/gdmt/optimize` | GDMT optimization recommendation |
| POST | `/protocol/recommend` | Imaging protocol recommendation |
| POST | `/reports/generate` | Generate report (markdown, JSON, PDF, FHIR) |
| GET | `/knowledge/stats` | Knowledge graph statistics |
| GET | `/metrics` | Prometheus-compatible metrics |
| GET | `/nim/status` | NIM service availability |

### 12.2 Streamlit UI (Port 8536)

Ten-tab interface following the Imaging Agent's 9-tab pattern:

| Tab | Purpose |
|---|---|
| **Evidence Explorer** | Multi-collection RAG Q&A with evidence citations, pre-filled cardiology example queries, Plotly donut chart for collection contribution |
| **Workflow Runner** | 11 clinical workflows with pre-loaded demo cases, annotated cardiac images, 6-step pipeline animation, download buttons |
| **Cardiac Imaging** | Echo, CT, MRI, nuclear imaging gallery with AI annotations, before/after toggle, 3D volume viewer |
| **Risk Calculators** | Interactive ASCVD, HEART, CHA₂DS₂-VASc, HAS-BLED, MAGGIC, EuroSCORE II calculators with guideline recommendations |
| **GDMT Optimizer** | Heart failure medication optimizer: input LVEF, labs, vitals → get pillar-by-pillar titration plan |
| **Device & AI Ecosystem** | 80+ FDA-cleared cardiac AI devices, searchable by modality and clinical task |
| **Protocol Advisor** | Patient-specific cardiac imaging protocol recommendations |
| **Reports & Export** | Markdown, JSON, NVIDIA-branded PDF, FHIR R4 DiagnosticReport export |
| **Patient 360** | Cross-modal cardiac dashboard: imaging + genomics + biomarkers + risk scores with interactive Plotly network graph |
| **Guidelines & Trials** | ACC/AHA guideline browser, landmark trial summaries, evidence level filtering |

### 12.3 Sidebar

| Section | Controls |
|---|---|
| **Guided Tour** | 10-step demo flow with dismiss button |
| **NIM Services** | 2x2 status grid (VISTA-3D, MAISI, VILA-M3, LLM) |
| **Collection Stats** | Metric widgets for each of 12 collections |
| **Filters** | CV domain dropdown, modality dropdown, year range slider |
| **Collections to Search** | Individual collection toggle checkboxes |
| **Demo Mode** | Load demo patient button |

### 12.4 Demo Cases

| ID | Title | Workflow | Key Features |
|---|---|---|---|
| DEMO-001 | Acute Chest Pain: STEMI Evaluation | coronary_artery_disease | ECG ST-elevation, troponin kinetics, cath lab activation |
| DEMO-002 | New Heart Failure: GDMT Initiation | heart_failure_gdmt | LVEF 25%, NYHA III, 4-pillar GDMT plan |
| DEMO-003 | Severe Aortic Stenosis: TAVR Evaluation | valvular_disease | AVA 0.8 cm², mean gradient 48 mmHg, STS risk |
| DEMO-004 | New-Onset AF: Stroke Prevention | arrhythmia_management | CHA₂DS₂-VASc 4, DOAC recommendation, rate control |
| DEMO-005 | Anthracycline Cardiotoxicity Surveillance | cardio_oncology_surveillance | Pre-chemo baseline, GLS monitoring protocol |

---

## 13. Clinical Decision Support Engines

### 13.1 Validated Risk Calculators

The agent will implement six validated cardiovascular risk calculators:

| Calculator | Use Case | Inputs | Output |
|---|---|---|---|
| **ASCVD (PCE)** | 10-year atherosclerotic CVD risk | Age, sex, race, TC, HDL, SBP, DM, smoking, HTN Rx | Risk %, statin recommendation |
| **HEART Score** | Chest pain risk stratification | History, ECG, age, risk factors, troponin | Score 0-10, risk category |
| **CHA₂DS₂-VASc** | AF stroke risk | CHF, HTN, age, DM, stroke, vascular, sex | Score 0-9, anticoagulation recommendation |
| **HAS-BLED** | AF bleeding risk on anticoagulation | HTN, renal/liver, stroke, bleeding, labile INR, elderly, drugs/alcohol | Score 0-9, risk category |
| **MAGGIC** | HF prognosis | Age, sex, LVEF, NYHA, SBP, BMI, creatinine, DM, COPD, HF duration, smoking, meds | 1- and 3-year mortality |
| **EuroSCORE II** | Cardiac surgery risk | Age, sex, renal function, extracardiac arteriopathy, poor mobility, prior cardiac surgery, COPD, active endocarditis, critical preop state, DM on insulin, NYHA, CCS angina, LVEF, recent MI, PA pressure, urgency, weight of procedure, thoracic aorta surgery | Operative mortality % |

### 13.2 GDMT Optimization Engine

The GDMT optimizer implements the ACC/AHA 2022 Heart Failure Guidelines algorithm:

```python
class GDMTOptimizer:
    """Heart failure guideline-directed medical therapy optimizer.

    Implements the 4-pillar HFrEF GDMT framework with:
    - Current medication assessment
    - Target dose calculation
    - Titration schedule recommendation
    - Contraindication checking
    - Drug interaction screening
    - Lab monitoring requirements (K+, Cr, eGFR)
    """

    def optimize(
        self,
        lvef: float,
        nyha_class: int,
        current_meds: List[Medication],
        labs: LabValues,
        vitals: Vitals,
        comorbidities: List[str],
    ) -> GDMTRecommendation:
        ...
```

---

## 14. Reporting and Interoperability

### 14.1 Export Formats

| Format | Use Case | Standards |
|---|---|---|
| Markdown | Clinical narrative, consultation notes | -- |
| JSON | Programmatic consumption, dashboards | -- |
| PDF | NVIDIA-themed clinical documentation | ReportLab |
| FHIR R4 | EHR integration, interoperability | SNOMED CT, LOINC, DICOM |

### 14.2 FHIR R4 Cardiovascular Coding

| Element | Code System | Example Codes |
|---|---|---|
| Findings | SNOMED CT | 22298006 (MI), 84114007 (HF), 49436004 (AF), 60234000 (AS) |
| Observations | LOINC | 10230-1 (LVEF), 33762-6 (NT-proBNP), 2093-3 (Total cholesterol), 18262-6 (LDL) |
| Procedures | SNOMED CT | 232717009 (CABG), 26212005 (PCI), 448076004 (TAVR) |
| Medications | RxNorm | Sacubitril/valsartan, carvedilol, dapagliflozin, apixaban |
| Imaging Studies | DICOM | US (echo), CT, MR, NM (nuclear) |
| Risk Scores | Custom extension | ASCVD %, CHA₂DS₂-VASc, HEART, MAGGIC |

### 14.3 Structured Report Sections

Cardiovascular reports will follow ACC/AHA reporting standards:

1. **Clinical Context** -- Indication, relevant history, medications
2. **Findings** -- Organized by organ system (LV, RV, valves, great vessels, coronary)
3. **Measurements** -- Quantitative data with reference ranges and severity grading
4. **Impression** -- Synthesized clinical interpretation
5. **Recommendations** -- Guideline-based next steps with evidence levels
6. **Genomic Enrichment** -- Cross-modal genetic findings (if triggered)
7. **Risk Scores** -- Calculated risk stratification with category
8. **Citations** -- Evidence sources with relevance scores

---

## 15. Product Requirements Document

### 15.1 Product Vision

**Vision Statement:** Enable any cardiologist, anywhere, to access integrated cardiovascular intelligence combining imaging AI, genomic analysis, guideline-based decision support, and evidence synthesis -- on a single $4,699 device.

**Target Users:** Community cardiologists, academic cardiology fellows, heart failure programs, structural heart teams, cardio-oncology clinics, preventive cardiology practices, clinical trial sites.

### 15.2 User Stories

#### Epic 1: Evidence-Based Clinical Queries

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-001 | As a cardiologist, I want to ask clinical questions and receive evidence-grounded answers with citations, so that I can make informed decisions. | P0 | Query returns answer with ≥3 citations from ≥2 collections; response time < 30 sec |
| US-002 | As a cardiologist, I want comparative analysis ("TAVR vs SAVR") with structured tables, so that I can compare treatment options. | P0 | Comparative query auto-detected; side-by-side evidence display; structured comparison table |
| US-003 | As a fellow, I want pre-filled example queries for common scenarios, so that I can learn the system quickly. | P1 | ≥4 clickable example queries; each returns relevant results |
| US-004 | As a researcher, I want to filter by CV domain, modality, and year, so that I can narrow my evidence search. | P1 | Sidebar filters applied to all queries; results reflect applied filters |

#### Epic 2: Clinical Workflows

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-005 | As a cardiologist, I want to run CAD assessment with calcium score and CTA data, so that I get CAD-RADS classification and management recommendations. | P0 | Workflow returns CAD-RADS 0-5, per-vessel stenosis, plaque characterization, guideline recommendation |
| US-006 | As an HF specialist, I want GDMT optimization for a patient with HFrEF, so that I get a 4-pillar titration plan. | P0 | All 4 GDMT pillars assessed; current vs target dose displayed; titration schedule provided |
| US-007 | As a structural heart team member, I want valve severity grading from echo measurements, so that I can assess intervention criteria. | P0 | Severity grading matches ASE criteria; SAVR vs TAVR recommendation based on STS risk |
| US-008 | As an EP physician, I want arrhythmia classification with CHA₂DS₂-VASc, so that I can manage AF stroke risk. | P0 | CHA₂DS₂-VASc calculated correctly; anticoagulation recommendation per guidelines |
| US-009 | As a cardio-oncologist, I want surveillance protocol generation, so that I know when to order echo and biomarkers. | P1 | Protocol specific to cancer therapy; GLS and LVEF thresholds defined; cardioprotection recommendations |

#### Epic 3: Risk Calculators

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-010 | As a preventive cardiologist, I want to calculate 10-year ASCVD risk, so that I can determine statin eligibility. | P0 | Correct PCE calculation; risk category assignment; statin intensity recommendation |
| US-011 | As an ED physician, I want HEART score for chest pain, so that I can risk-stratify patients. | P0 | Score 0-10 calculated; risk category (low/mod/high); disposition recommendation |
| US-012 | As a cardiologist, I want interactive risk calculator forms, so that I can enter patient data directly in the UI. | P1 | Form inputs with validation; real-time score update; guideline recommendation display |

#### Epic 4: Cross-Modal Integration

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-013 | As a cardiomyopathy specialist, I want genetic triggers from cardiac imaging, so that DCM patients automatically get gene panel recommendations. | P0 | Non-ischemic DCM + age < 50 triggers gene panel; genomic hits displayed in Patient 360 |
| US-014 | As a cardiologist, I want a Patient 360 view combining imaging, genomics, biomarkers, and risk scores, so that I see the complete picture. | P1 | Interactive Plotly network graph; cross-modal connections displayed; drill-down to evidence |

#### Epic 5: Reporting and Export

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-015 | As a cardiologist, I want to export clinical reports as PDF, so that I can include them in patient records. | P0 | NVIDIA-branded PDF with clinical question, analysis, evidence, risk scores; download button |
| US-016 | As a health IT engineer, I want FHIR R4 DiagnosticReport export, so that findings integrate with our EHR. | P1 | Valid FHIR R4 Bundle; SNOMED CT + LOINC coded; passes FHIR validator |
| US-017 | As a cardiologist, I want 4 export formats (Markdown, JSON, PDF, FHIR), so that I can choose the right format for my use case. | P1 | All 4 formats functional; consistent content across formats |

#### Epic 6: Demo and Presentation

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-018 | As a demo presenter, I want 5 pre-loaded demo cases, so that I can run the demo without entering data. | P0 | 5 demo cases selectable from dropdown; each runs in < 30 seconds; realistic clinical output |
| US-019 | As a new user, I want a sidebar guided tour, so that I understand the 10-tab interface quickly. | P1 | Expandable tour with numbered steps; dismiss button; persists across session |
| US-020 | As a demo presenter, I want a cardiac imaging gallery, so that I can show impressive AI-annotated images. | P1 | Echo, CT, MRI images with AI overlays; before/after toggle; 3D slice viewer |

### 15.3 Non-Functional Requirements

| Requirement | Target | Rationale |
|---|---|---|
| RAG query latency | < 30 seconds end-to-end | Acceptable for clinical consultation workflow |
| Risk calculator latency | < 5 seconds | Near-instantaneous for interactive use |
| Workflow execution (mock) | < 10 seconds | Responsive demo experience |
| Availability | 99.5% uptime | Clinical support tool (not life-critical) |
| Memory footprint | < 16 GB (agent only) | Coexist with other agents on 128 GB DGX Spark |
| Seed data completeness | 1,500+ records across 12 collections | Sufficient for meaningful RAG retrieval |
| Unit test coverage | > 80% | Reliable development cycle |
| FHIR R4 compliance | Passes HL7 FHIR Validator | Interoperability requirement |

### 15.4 Prioritization Matrix

| Phase | Features | Timeline |
|---|---|---|
| **Phase 1 (MVP)** | RAG engine (12 collections), Evidence Explorer, 3 workflows (CAD, HF, arrhythmia), 3 risk calculators (ASCVD, CHA₂DS₂-VASc, HEART), 3 demo cases, PDF export | 4-6 weeks |
| **Phase 2 (Complete)** | All 11 workflows, all 6 calculators, GDMT optimizer, FHIR R4 export, cardiac imaging gallery, Patient 360, all 5 demo cases | 4-6 weeks |
| **Phase 3 (Polish)** | Guided tour, pipeline animation, cross-modal triggers, network graph, Guidelines & Trials tab, benchmark validation | 2-3 weeks |

---

## 16. Data Acquisition Strategy

### 16.1 Automated Ingest Pipelines

| Source | Collection(s) | Method | Update Cadence |
|---|---|---|---|
| PubMed (NCBI E-utilities) | cardio_literature | MeSH-filtered abstract retrieval | Weekly |
| ClinicalTrials.gov (V2 API) | cardio_trials | Cardiovascular condition filter | Weekly |
| ACC/AHA Guidelines PDFs | cardio_guidelines | Manual curation + embedding | Per guideline update |
| ASE/SCCT/SCMR references | cardio_imaging | Manual curation | Per guideline update |

### 16.2 Curated Seed Data

| Collection | Records | Curation Source |
|---|---|---|
| cardio_imaging | 200 | ASE, SCCT, SCMR, ASNC guidelines and reference texts |
| cardio_electrophysiology | 150 | ACC/AHA/HRS guidelines, EP textbooks |
| cardio_heart_failure | 150 | ACC/AHA 2022 HF Guidelines, landmark trial data |
| cardio_valvular | 120 | ACC/AHA 2020 VHD Guidelines |
| cardio_prevention | 150 | ACC/AHA 2019 Prevention, 2018 Cholesterol Guidelines |
| cardio_interventional | 100 | SCAI guidelines, procedural references |
| cardio_oncology | 100 | ESC 2022 Cardio-Oncology Guidelines |
| cardio_devices | 80 | FDA AI/ML database, manufacturer data |
| cardio_hemodynamics | 80 | Catheterization references and guidelines |

### 16.3 PubMed Search Strategy

```
MeSH Terms:
  "Cardiovascular Diseases"[MeSH] OR "Heart Diseases"[MeSH] OR
  "Coronary Artery Disease"[MeSH] OR "Heart Failure"[MeSH] OR
  "Atrial Fibrillation"[MeSH] OR "Cardiomyopathies"[MeSH]

AND ("Artificial Intelligence"[MeSH] OR "Machine Learning"[MeSH] OR
     "Deep Learning"[MeSH] OR "Neural Networks"[MeSH])

Filters: Published 2018-2026, English, Humans
Expected: 3,000-5,000 abstracts
```

---

## 17. Validation and Testing Strategy

### 17.1 Unit Tests

| Test Category | Target Count | Coverage |
|---|---|---|
| Collection schemas | 36 | 12 collections × 3 tests each |
| Risk calculators | 48 | 6 calculators × 8 test cases each |
| GDMT optimizer | 30 | Multiple HF phenotypes and medication combos |
| Workflow logic | 40 | 11 workflows × 5 test cases each |
| RAG engine | 20 | Query expansion, scoring, synthesis |
| Knowledge graph | 15 | Entity lookup, alias resolution |
| API endpoints | 30 | All REST endpoints |
| FHIR R4 export | 15 | Schema validation, coding accuracy |
| Cross-modal triggers | 12 | All trigger conditions |
| **Total** | **~250** | |

### 17.2 Clinical Validation

| Validation Type | Method | Target |
|---|---|---|
| Risk calculator accuracy | Compare against published validation cohorts | < 1% deviation from reference implementations |
| GDMT recommendations | Review by board-certified HF cardiologist | 95%+ guideline concordance |
| Severity grading | Compare against ASE/ACC criteria | 100% match for standard inputs |
| Guideline citation accuracy | Verify against source guidelines | 100% accurate citations |
| FHIR R4 compliance | HL7 FHIR Validator | Zero validation errors |

### 17.3 End-to-End Validation Checks

| Check | Criteria |
|---|---|
| Health endpoint | Returns status=healthy, all 12 collections with non-zero counts |
| RAG query | Returns answer with ≥3 citations in < 30 seconds |
| Risk calculator | All 6 calculators return correct results for known inputs |
| Workflow execution | All 11 workflows complete in mock mode |
| FHIR export | Valid R4 Bundle passes FHIR validator |
| PDF export | Generates downloadable PDF with NVIDIA branding |
| Demo cases | All 5 demo cases execute successfully |
| Cross-modal trigger | Genomic query fires for qualifying conditions |
| Comparative analysis | "TAVR vs SAVR" produces structured comparison |

---

## 18. Regulatory Considerations

### 18.1 Intended Use Classification

The Cardiology Intelligence Agent is a **clinical decision support (CDS) tool** intended for use by licensed healthcare professionals. It is designed to:

- Retrieve and synthesize published cardiovascular evidence
- Calculate validated clinical risk scores
- Generate guideline-based recommendations
- Assist with (not replace) clinical decision-making

### 18.2 FDA CDS Exemption Criteria (21st Century Cures Act)

Under the 21st Century Cures Act, software functions meeting all four criteria are exempt from FDA device regulation:

| Criterion | Assessment |
|---|---|
| Not intended to acquire, process, or analyze a medical image or signal | **Met** -- RAG and risk calculation; does not process raw medical images or signals |
| Intended for displaying, analyzing, or printing medical information | **Met** -- Displays evidence and calculations |
| Intended for use by healthcare professional | **Met** -- Designed for cardiologists |
| Healthcare professional does not primarily rely on the software | **Met** -- Provides recommendations for review, not autonomous decisions |

**Assessment:** The core RAG and risk calculator functions likely qualify for CDS exemption. Clinical workflows involving NIM-based image analysis would require separate regulatory consideration.

### 18.3 Disclaimers

All outputs will include standard disclaimers:

> *This tool is for clinical decision support only and does not replace professional medical judgment. All findings require verification by a qualified cardiologist. Risk calculations are estimates based on population-level data and may not apply to individual patients. Not FDA-cleared for autonomous clinical decision-making. For research and educational purposes only.*

---

## 19. DGX Compute Progression

| Phase | Hardware | Price | Scope |
|---|---|---|---|
| **Phase 1 -- Proof Build** | DGX Spark | $4,699 | All 11 workflows (mock/cloud NIM), 12 collections, 6 risk calculators |
| **Phase 2 -- Departmental** | 1-2x DGX B200 | $500K-$1M | Full NIM stack, live echo/CT/MRI processing, PACS integration |
| **Phase 3 -- Multi-Site** | 4-8x DGX B200 | $2M-$4M | NVIDIA FLARE federated learning across sites, population analytics |
| **Phase 4 -- AI Factory** | DGX SuperPOD | $7M-$60M+ | Thousands concurrent studies, real-time ICU monitoring, national registries |

---

## 20. Implementation Roadmap

### Phase 1: Foundation (Weeks 1-6)

| Week | Deliverable |
|---|---|
| 1-2 | Repository scaffolding, Pydantic models, settings, Docker Compose, 12 collection schemas |
| 3-4 | Seed data curation (1,530 records), ingest pipelines (PubMed, ClinicalTrials.gov), embedding generation |
| 5-6 | RAG engine (parallel search, weighted scoring, query expansion, Claude synthesis), Evidence Explorer tab |

### Phase 2: Clinical Intelligence (Weeks 7-12)

| Week | Deliverable |
|---|---|
| 7-8 | 3 priority workflows (CAD, HF, arrhythmia), risk calculators (ASCVD, CHA₂DS₂-VASc, HEART), demo cases |
| 9-10 | Remaining 5 workflows (valve, CMR, stress, prevention, cardio-onc), remaining calculators (HAS-BLED, MAGGIC, EuroSCORE II) |
| 11-12 | GDMT optimizer, cross-modal triggers, FHIR R4 export, PDF reports |

### Phase 3: UI and Polish (Weeks 13-16)

| Week | Deliverable |
|---|---|
| 13-14 | 10-tab Streamlit UI, cardiac imaging gallery, Patient 360 network graph, Guidelines & Trials browser |
| 15-16 | Sidebar guided tour, pipeline animation, pre-filled examples, 5 demo cases finalized, documentation |

### Phase 4: Integration and Validation (Weeks 17-18)

| Week | Deliverable |
|---|---|
| 17 | Integration with Imaging Agent (shared cardiac workflows), cross-agent triggers, Docker Compose integration |
| 18 | Clinical validation, end-to-end testing (250+ unit tests), documentation publication (demo guide, design doc, project bible) |

---

## 21. Risk Analysis

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Risk calculator implementation errors | Medium | High | Validate against published reference implementations; extensive unit testing |
| GDMT logic complexity | Medium | Medium | Start with HFrEF only (best-defined); expand to HFmrEF/HFpEF later |
| Insufficient seed data for niche collections | Medium | Medium | Focus on high-impact collections first; iteratively expand |
| NIM cardiac segmentation accuracy | Low | Medium | Mock mode fallback; VISTA-3D already supports cardiac structures |
| Guideline updates during development | Low | Low | Modular guideline collection; easy to update individual recommendations |
| Memory pressure with 12 collections + existing agents | Low | Medium | BGE-small (384-dim) is compact; Milvus handles multi-collection efficiently |
| Claude API rate limits during demo | Low | High | Llama-3 NIM fallback; response caching for demo scenarios |

---

## 22. Competitive Landscape

### 22.1 Positioning

The Cardiology Intelligence Agent occupies a unique position in the cardiovascular AI landscape:

```
                    Multi-Modal Integration
                           ↑
                           |
                    [Cardiology Agent]
                           |
         Cloud-only ←------+------→ On-device
                           |
                    [HeartFlow, Cleerly]
                           |
                           ↓
                    Single-Modality
```

No existing product combines:

1. Multi-modal cardiac imaging AI (echo + CT + MRI + nuclear + ECG)
2. Genomic integration with cross-modal triggers
3. RAG-grounded evidence synthesis with citations
4. Validated risk calculators (ASCVD, HEART, CHA₂DS₂-VASc)
5. GDMT optimization engine
6. On-device deployment ($4,699)
7. Open-source (Apache 2.0)

### 22.2 Defensibility

| Advantage | Defensibility |
|---|---|
| Multi-collection RAG architecture | High -- 5 agents prove the pattern; cardiology is the 6th instantiation |
| Cross-modal genomic triggers | High -- unique to HCLS AI Factory; requires integrated platform |
| On-device inference | High -- DGX Spark + NIM stack is NVIDIA-exclusive |
| Open source | Medium -- community contribution, institutional customization |
| Guideline-aligned CDS | Low -- guidelines are public; competitors can implement |
| Clinical validation | Medium -- requires domain expertise and clinical partnerships |

---

## 23. Discussion

### 23.1 Why Cardiology Is the Right Next Agent

Of all medical specialties, cardiology presents the strongest case for an integrated intelligence agent:

1. **Data richness**: No other specialty routinely integrates as many data modalities (imaging, EP, hemodynamics, biomarkers, genomics, risk scores) into a single clinical encounter.

2. **Quantitative precision**: Cardiovascular medicine is measurement-driven. LVEF, calcium scores, valve gradients, QTc intervals, and risk percentages are all inherently structured data -- ideal for RAG retrieval and clinical decision support.

3. **Guideline density**: The ACC/AHA and ESC produce more clinical practice guidelines than any other medical specialty. These guidelines are highly structured (Class I-III, LOE A-C) and directly mappable to clinical decision support rules.

4. **Clinical workflow fit**: Cardiologists already use structured reporting, risk calculators, and decision algorithms. The Cardiology Intelligence Agent enhances rather than disrupts existing workflow.

5. **Market size**: Cardiovascular AI is the largest segment of the medical AI market ($2.8B projected by 2028). Every health system has cardiologists. The addressable market is essentially universal.

6. **Existing foundation**: The HCLS AI Factory already has cardiac CT workflows (DEMO-003 in the Imaging Agent), coronary segmentation via VISTA-3D, and cardiovascular genomics in the shared variant collection. The Cardiology Agent builds on existing infrastructure rather than starting from scratch.

### 23.2 Limitations

1. **Mock mode inference**: Phase 1 uses clinically realistic mock results rather than live model inference. This is appropriate for proof-of-concept and demo but requires GPU deployment for clinical validation.

2. **Guideline currency**: Clinical guidelines are updated periodically. The agent requires manual curation to incorporate focused updates, although the modular collection architecture makes updates straightforward.

3. **Risk calculator applicability**: Validated risk calculators (PCE, HEART) were developed in specific populations and may not generalize equally across all demographic groups. The agent should display appropriate caveats.

4. **Not a replacement for clinical judgment**: The agent assists with evidence synthesis and calculation but cannot replace the nuanced clinical judgment of an experienced cardiologist, particularly in complex or atypical presentations.

### 23.3 Future Directions

1. **Real-time ICU integration**: Continuous hemodynamic monitoring with AI-driven trend analysis and early warning for cardiac decompensation.

2. **Wearable data integration**: Apple Watch ECG, continuous glucose monitoring, and remote blood pressure data feeding into longitudinal risk assessment.

3. **Population health analytics**: Cohort-level cardiovascular outcomes analysis across institutions using federated learning (NVIDIA FLARE).

4. **Clinical trial matching**: Automated screening of patients for cardiovascular clinical trial eligibility based on imaging, biomarker, and genomic criteria.

5. **Cardiac digital twins**: Patient-specific computational models combining imaging anatomy, electrophysiology, and hemodynamics for treatment simulation.

---

## 24. Conclusion

The Cardiology Intelligence Agent represents the natural next extension of the HCLS AI Factory platform into the highest-impact medical specialty. By leveraging proven multi-collection RAG patterns, shared infrastructure, and established NIM integration -- while adding cardiology-specific clinical workflows, validated risk calculators, GDMT optimization, and cross-modal genomic triggers -- the agent will deliver integrated cardiovascular intelligence that currently requires multi-million-dollar institutional investments.

The 12-collection Milvus architecture, 8 reference clinical workflows, 6 validated risk calculators, and GDMT optimization engine provide a comprehensive decision support platform that covers the full spectrum of cardiovascular practice: from preventive risk stratification through acute coronary syndromes to advanced heart failure management and cardio-oncology surveillance.

Deploying on a single NVIDIA DGX Spark ($4,699) ensures that this intelligence is accessible not just to elite academic centers but to community cardiologists, rural health systems, and resource-limited institutions worldwide. Combined with Apache 2.0 licensing and the HCLS AI Factory's proven open-source approach, the Cardiology Intelligence Agent will democratize cardiovascular AI in a way that no existing commercial product achieves.

---

## 25. References

1. Virani SS, et al. Heart Disease and Stroke Statistics -- 2025 Update. *Circulation*. 2025;151:e1-e375.
2. Heidenreich PA, et al. 2022 AHA/ACC/HFSA Guideline for the Management of Heart Failure. *Circulation*. 2022;145:e895-e1032.
3. Writing Committee, et al. 2021 ACC/AHA/SCAI Guideline for Coronary Artery Revascularization. *J Am Coll Cardiol*. 2022;79:e21-e129.
4. Otto CM, et al. 2020 ACC/AHA Guideline for the Management of Patients with Valvular Heart Disease. *Circulation*. 2021;143:e72-e227.
5. January CT, et al. 2019 AHA/ACC/HRS Focused Update of the 2014 Guideline for Management of Atrial Fibrillation. *Circulation*. 2019;140:e125-e151.
6. Arnett DK, et al. 2019 ACC/AHA Guideline on the Primary Prevention of Cardiovascular Disease. *Circulation*. 2019;140:e596-e646.
7. Grundy SM, et al. 2018 AHA/ACC Cholesterol Clinical Practice Guideline. *Circulation*. 2019;139:e1082-e1143.
8. Lyon AR, et al. 2022 ESC Guidelines on Cardio-Oncology. *Eur Heart J*. 2022;43:4229-4361.
9. Vahanian A, et al. 2021 ESC/EACTS Guidelines for the Management of Valvular Heart Disease. *Eur Heart J*. 2022;43:561-632.
10. Hindricks G, et al. 2020 ESC Guidelines for the Diagnosis and Management of Atrial Fibrillation. *Eur Heart J*. 2021;42:373-498.
11. Ommen SR, et al. 2024 AHA/ACC Guideline for the Diagnosis and Management of Hypertrophic Cardiomyopathy. *Circulation*. 2024;149:e1239-e1311.
12. Cury RC, et al. CAD-RADS 2.0 -- 2022 Coronary Artery Disease Reporting and Data System. *Radiology*. 2022;305:209-221.
13. Lang RM, et al. Recommendations for Cardiac Chamber Quantification by Echocardiography in Adults: An Update from the ASE. *J Am Soc Echocardiogr*. 2015;28:1-39.
14. Kramer CM, et al. Standardized Cardiovascular Magnetic Resonance Imaging (CMR) Protocols: 2020 Update. *J Cardiovasc Magn Reson*. 2020;22:17.
15. Goff DC, et al. 2013 ACC/AHA Guideline on the Assessment of Cardiovascular Risk. *Circulation*. 2014;129:S49-S73.
16. Lip GY, et al. Refining Clinical Risk Stratification for Predicting Stroke and Thromboembolism in Atrial Fibrillation Using a Novel Risk Factor-Based Approach: The Euro Heart Survey. *Chest*. 2010;137:263-272.
17. Six AJ, et al. The HEART Score for the Assessment of Patients with Chest Pain in the Emergency Department. *Neth Heart J*. 2008;16:191-196.
18. Pocock SJ, et al. Predicting Survival in Heart Failure: A Risk Score Based on 39,372 Patients from 30 Studies (MAGGIC). *Eur J Heart Fail*. 2013;15:1082-1094.
19. Nashef SA, et al. EuroSCORE II. *Eur J Cardiothorac Surg*. 2012;41:734-744.
20. World Health Organization. Cardiovascular Diseases (CVDs) Fact Sheet. WHO, 2024.
21. FDA. Artificial Intelligence and Machine Learning (AI/ML)-Enabled Medical Devices. FDA Database, 2025.
22. McMurray JJ, et al. Angiotensin-Neprilysin Inhibition versus Enalapril in Heart Failure (PARADIGM-HF). *N Engl J Med*. 2014;371:993-1004.
23. McMurray JJ, et al. Dapagliflozin in Patients with Heart Failure and Reduced Ejection Fraction (DAPA-HF). *N Engl J Med*. 2019;381:1995-2008.
24. Packer M, et al. Cardiovascular and Renal Outcomes with Empagliflozin in Heart Failure (EMPEROR-Reduced). *N Engl J Med*. 2020;383:1413-1424.
25. Solomon SD, et al. Dapagliflozin in Heart Failure with Mildly Reduced or Preserved Ejection Fraction (DELIVER). *N Engl J Med*. 2022;387:1089-1098.
26. Maron DJ, et al. Initial Invasive or Conservative Strategy for Stable Coronary Disease (ISCHEMIA). *N Engl J Med*. 2020;382:1395-1407.

---

*HCLS AI Factory -- Cardiology Intelligence Agent Research Paper and PRD*
*Apache 2.0 License | March 2026*
*Author: Adam Jones*
