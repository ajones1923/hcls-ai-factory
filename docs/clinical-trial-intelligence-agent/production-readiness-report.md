# Clinical Trial Intelligence Agent -- Production Readiness & Capability Report

**Version:** 2.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Status:** Production Demo Ready (10/10)
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory
**License:** Apache 2.0

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [System Architecture](#2-system-architecture)
3. [Knowledge Graph](#3-knowledge-graph)
4. [Clinical Workflows](#4-clinical-workflows)
5. [Decision Support Engines](#5-decision-support-engines)
6. [Cross-Agent Integration](#6-cross-agent-integration)
7. [Vector Database & Collections](#7-vector-database--collections)
8. [RAG Engine](#8-rag-engine)
9. [Query Expansion System](#9-query-expansion-system)
10. [Autonomous Agent Pipeline](#10-autonomous-agent-pipeline)
11. [Data Models & Type Safety](#11-data-models--type-safety)
12. [Streamlit UI](#12-streamlit-ui)
13. [REST API](#13-rest-api)
14. [Data Ingest Pipelines](#14-data-ingest-pipelines)
15. [Seed Data Inventory](#15-seed-data-inventory)
16. [Export & Reporting](#16-export--reporting)
17. [Observability & Metrics](#17-observability--metrics)
18. [Scheduling & Automation](#18-scheduling--automation)
19. [Configuration System](#19-configuration-system)
20. [Security & Authentication](#20-security--authentication)
21. [Infrastructure & Deployment](#21-infrastructure--deployment)
22. [Test Suite](#22-test-suite)
23. [Known Limitations](#23-known-limitations)
24. [Demo Readiness Audit](#24-demo-readiness-audit)
25. [Codebase Summary](#25-codebase-summary)

---

## 1. Executive Summary

The Clinical Trial Intelligence Agent is a production-grade, RAG-powered decision support system for clinical trial design, execution, and optimization, built for the HCLS AI Factory precision medicine platform running on NVIDIA DGX Spark. It delivers evidence-based guidance across the full clinical trial lifecycle -- from protocol design and patient-trial matching through adaptive design evaluation, safety signal detection, regulatory document generation, competitive intelligence, diversity assessment, and decentralized trial planning -- by combining 14 domain-specific Milvus vector collections, five calibrated decision support engines, 10 clinical workflows, and an autonomous reasoning pipeline that plans, searches, evaluates, and synthesizes trial evidence in real time.

The agent is architected as a three-tier system: a 5-tab Streamlit UI (port 8128) for interactive clinical trial exploration, a FastAPI REST API (port 8538) exposing 26 endpoints for programmatic integration, and a RAG engine backed by Milvus (port 19530) with BGE-small-en-v1.5 384-dimensional embeddings. All 10 clinical workflows, all 5 decision support engines, and the full query expansion system operate independently of Milvus connectivity, ensuring graceful degradation and robust demo capability even when the vector store is unavailable.

The codebase comprises 46 Python files totaling 22,607 lines of code, with 12 dedicated test files containing 769 passing tests at a 100% pass rate (0.47s execution time). The knowledge base contains 40 landmark trials, 13 therapeutic areas, 9 regulatory agencies, 9 endpoint types, 9 adaptive trial designs, 9 biomarker strategies, 9 decentralized trial components, 7 trial phases, 140 entity aliases, and 33 drug synonym mappings. This report documents every capability, data dimension, and test result to serve as the definitive long-term reference for the Clinical Trial Intelligence Agent.

| Capability | Detail |
|---|---|
| Clinical Workflows | 10 types (Protocol Design, Patient Matching, Site Selection, Eligibility Optimization, Adaptive Design, Safety Signal, Regulatory Docs, Competitive Intel, Diversity Assessment, Decentralized Planning) |
| Decision Support Engines | 5 (Confidence Calibrator, Protocol Complexity Scorer, Enrollment Predictor, Eligibility Analyzer, Competitive Threat Scorer) |
| Historical Success Estimator | 12 therapeutic areas with phase-specific success rates |
| Cross-Agent Triggers | 4 agents (Oncology, PGx, Cardiology, Biomarker) |
| Vector Collections | 14 Milvus collections (IVF_FLAT, COSINE, 384-dim) |
| Knowledge Graph | 40 landmark trials, 13 therapeutic areas, 9 agencies, 9 endpoint types |
| Adaptive Designs | 9 validated designs with regulatory guidance |
| Biomarker Strategies | 9 strategies from enrichment to digital biomarkers |
| DCT Components | 9 decentralized trial components |
| Query Expansion | 10 synonym maps, 140 entity aliases, 33 drug synonyms, 22 biomarkers |
| Tests | 769 passed, 0 failed, 100% pass rate, 0.47s |
| Source LOC | ~16,300 (34 source files) |
| Test LOC | ~6,300 (12 test files) |
| Total Python LOC | 22,607 |
| API Endpoints | 26 |
| Prometheus Metrics | 20+ metrics across query, RAG, workflow, and system health |
| Ports | FastAPI 8538, Streamlit 8128, Milvus 19530 |
| Authentication | API key (X-API-Key header) |
| Rate Limiting | 100 requests/minute per IP |
| Event Publishing | SSE event stream for real-time updates |
| Export Formats | Markdown, JSON, PDF |
| Documentation | Research paper + 7 reference documents |
| Versioning | Knowledge v2.0.0, API v1, Agent v2.0.0 |

---

## 2. System Architecture

### Three-Tier Architecture

| Tier | Component | Technology | Port | Purpose |
|---|---|---|---|---|
| **Presentation** | Streamlit UI | Streamlit + NVIDIA Dark Theme | 8128 | Interactive 5-tab clinical trial exploration |
| **Application** | FastAPI REST API | FastAPI + Uvicorn | 8538 | 26 endpoints, CORS, rate limiting, auth |
| **Data** | Milvus Vector Store | Milvus + etcd + MinIO | 19530 | 14 collections, BGE-small-en-v1.5 embeddings |

### System Diagram

```
                    +----------------------------+
                    |    Streamlit UI (:8128)     |
                    |   5 Tabs, NVIDIA Theme      |
                    +-------------+--------------+
                                  |
                                  | HTTP/REST
                                  v
                    +----------------------------+
                    |    FastAPI API (:8538)      |
                    |  26 Endpoints, Auth, CORS   |
                    +---+------+------+------+---+
                        |      |      |      |
              +---------+  +---+  +---+  +---+--------+
              v          v       v       v             v
     +--------+--+ +-----+----+ +---+---+ +---+------+----+
     | Clinical   | | Decision | | RAG   | | Cross-Agent   |
     | Workflows  | | Support  | | Engine| | Integration   |
     | (10 types) | | (5 eng)  | |       | | (4 agents)    |
     +--------+--+ +----------+ +---+---+ +---------------+
              |                      |
              v                      v
     +--------+--+          +--------+--------+
     | Knowledge  |          | Milvus (:19530)  |
     | Base       |          | 14 Collections   |
     | (40 trials |          | IVF_FLAT/COSINE  |
     |  13 areas) |          | 384-dim BGE      |
     +------------+          +---------+--------+
                                       |
                             +---------+---------+
                             |    etcd + MinIO    |
                             | (metadata + blobs) |
                             +-------------------+
```

### Component Map

| File | LOC | Purpose |
|---|---|---|
| `src/clinical_workflows.py` | 2,606 | 10 clinical trial workflows with preprocess/execute/postprocess |
| `src/knowledge.py` | 1,929 | Domain knowledge: 40 trials, 13 areas, 9 agencies, phases, designs |
| `src/agent.py` | 1,678 | Autonomous reasoning: plan, search, evaluate, synthesize |
| `src/rag_engine.py` | 1,555 | Multi-collection RAG with weighted retrieval |
| `src/collections.py` | 1,221 | 14 Milvus collection schemas (IVF_FLAT/COSINE/384-dim) |
| `src/query_expansion.py` | 1,150 | 10 synonym maps, 140 entity aliases, workflow-aware boosting |
| `api/routes/trial_clinical.py` | 1,079 | 22 trial-specific API endpoints |
| `src/ingest/clinicaltrials_parser.py` | 802 | ClinicalTrials.gov XML/JSON parser |
| `src/decision_support.py` | 731 | 5 decision support engines |
| `src/export.py` | 725 | Multi-format report export (Markdown, JSON, PDF) |
| `app/trial_ui.py` | 692 | 5-tab Streamlit UI with NVIDIA dark theme |
| `api/main.py` | 615 | FastAPI app factory, CORS, health, metrics |
| `src/scheduler.py` | 577 | Scheduled ingest and collection maintenance |
| `src/metrics.py` | 525 | Prometheus-compatible metrics (20+ metrics) |
| `src/ingest/regulatory_parser.py` | 505 | FDA/EMA/ICH regulatory document parser |
| `src/models.py` | 496 | Pydantic models, enums, type-safe data contracts |
| `src/ingest/pubmed_parser.py` | 376 | PubMed/MEDLINE literature parser |
| `src/cross_modal.py` | 358 | Cross-agent query and integration |
| `scripts/setup_collections.py` | 328 | Collection creation and index management |
| `api/routes/reports.py` | 312 | Report generation and export endpoints |
| `api/routes/events.py` | 269 | SSE event stream and event health |
| `scripts/seed_knowledge.py` | 277 | Knowledge base seeding from domain data |
| `src/ingest/base.py` | 228 | Base ingest pipeline with validation |
| `config/settings.py` | 182 | Pydantic BaseSettings configuration |
| `scripts/run_ingest.py` | 142 | CLI entry point for ingest pipeline |

---

## 3. Knowledge Graph

### 3.1 Therapeutic Areas (13)

| Area | Key Indications | Phase 3 Success Rate | Key Biomarkers |
|---|---|---|---|
| Oncology | NSCLC, breast, CRC, melanoma, AML, MM | 36% | PD-L1, TMB, MSI-H, HER2, EGFR, ALK, BRCA |
| Cardiovascular | HF, ACS, AF, hypertension, PAH | 35% | NT-proBNP, troponin, LDL-C, Lp(a), hsCRP |
| Neuroscience | AD, PD, MS, MDD, epilepsy, migraine | 28% | amyloid PET, tau PET, NfL, p-tau217 |
| Immunology | RA, psoriasis, IBD, SLE, atopic dermatitis | 40% | CRP, ESR, anti-CCP, RF, ANA, IL-6 |
| Infectious Disease | HIV, HBV, HCV, COVID-19, RSV, AMR | 58% | viral load, CD4, seroconversion, PCR, MIC |
| Rare Diseases | SMA, DMD, CF, sickle cell, hemophilia | 45% | genetic mutation status, enzyme activity |
| Metabolic | T2DM, obesity, NASH/MASH, dyslipidemia | 42% | HbA1c, FPG, BMI, ALT/AST, LDL-C |
| Respiratory | asthma, COPD, IPF, PAH | 39% | FEV1, eosinophils, FeNO, IgE, FVC |
| Hematology | AML, CLL, MM, MDS, sickle cell | 41% | MRD, CR rate, transfusion independence |
| Gastroenterology | Crohn's, UC, GERD, NASH/MASH | 38% | fecal calprotectin, CRP, endoscopic scores |
| Dermatology | atopic dermatitis, psoriasis, alopecia | 44% | EASI, PASI, IGA, BSA, DLQI |
| Ophthalmology | AMD, DME, glaucoma, dry eye | 43% | BCVA, CRT, IOP, GA lesion area |
| Gene/Cell Therapy | SMA, hemophilia, CAR-T, retinal dystrophy | 40% | transgene expression, T-cell expansion |

### 3.2 Trial Phases (7)

Preclinical, Phase 0 (Exploratory IND), Phase 1, Phase 2, Phase 3 (Pivotal), Phase 4 (Post-Marketing), and Expanded Access / Compassionate Use -- each with typical duration, enrollment ranges, primary objectives, key endpoints, and regulatory requirements.

### 3.3 Regulatory Agencies (9)

FDA, EMA, PMDA, Health Canada, TGA, MHRA, NMPA, Swissmedic, ANVISA -- each with approval pathways, review timelines, key requirements, and expedited programs.

### 3.4 Landmark Trials (40)

| Trial | Phase | Indication | Key Finding |
|---|---|---|---|
| KEYNOTE-024 | Phase 3 | NSCLC (PD-L1 high) | Pembrolizumab PFS HR 0.50; first-line IO monotherapy |
| EMPEROR-Reduced | Phase 3 | HFrEF | Empagliflozin reduced CV death/HF hosp by 25% |
| RECOVERY | Phase 2/3 | COVID-19 | Dexamethasone reduced mortality by 1/3 in ventilated; largest platform trial |
| PARADIGM-HF | Phase 3 | HFrEF | ARNI 20% reduction vs ACEi; new HF standard |
| CheckMate-067 | Phase 3 | Melanoma | Nivo+Ipi 52% 5-year OS; combination IO paradigm |
| SPRINT | Phase 3 | Hypertension | Intensive BP control 25% MACE reduction |
| EMPA-REG OUTCOME | Phase 3 | T2DM + CVD | First diabetes drug with CV mortality benefit |
| DAPA-CKD | Phase 3 | CKD | SGLT2i to non-diabetic CKD; HR 0.61 |
| HIMALAYA | Phase 3 | HCC | STRIDE single anti-CTLA-4 priming dose |
| DESTINY-Breast04 | Phase 3 | HER2-low breast | T-DXd created new targetable population |
| CLARITY-AD | Phase 3 | Early AD | Lecanemab 27% CDR-SB reduction; amyloid validation |
| FOURIER | Phase 3 | ASCVD | PCSK9 inhibition reduces CV events |
| FLAURA | Phase 3 | EGFR+ NSCLC | Osimertinib first-line; OS benefit |
| ADAURA | Phase 3 | EGFR+ NSCLC (adj) | Adjuvant osimertinib DFS benefit |
| I-SPY 2 | Phase 2 | Breast cancer | Bayesian adaptive randomization platform |
| VICTORIA | Phase 3 | Worsening HF | Vericiguat added to HF armamentarium |
| DESTINY-Lung02 | Phase 2 | HER2+ NSCLC | ADC in non-breast HER2+ tumors |
| CLEAR Outcomes | Phase 3 | CV risk | Bempedoic acid MACE reduction |
| SELECT | Phase 3 | Overweight + CVD | Semaglutide 20% MACE reduction |
| STEP HFpEF | Phase 3 | HFpEF + obesity | Semaglutide improved symptoms/function |
| SURMOUNT-1 | Phase 3 | Obesity | Tirzepatide 22.5% weight loss |
| CASGEVY | Regulatory | SCD/TDT | First CRISPR gene therapy approved |
| DAPA-HF | Phase 3 | HFrEF | Dapagliflozin HF benefit regardless of diabetes |
| *...and 17 more* | | | |

---

## 4. Clinical Workflows

The agent implements 10 clinical trial workflows, each following the `BaseTrialWorkflow` contract: `preprocess` (validate/normalize inputs), `execute` (core logic), `postprocess` (enrich results). All workflows are registered in the `WorkflowEngine` for unified dispatch.

| # | Workflow | Type | Key Outputs |
|---|---|---|---|
| 1 | Protocol Design | `protocol_design` | Protocol blueprint, endpoint selection, sample size, comparator strategy |
| 2 | Patient Matching | `patient_matching` | Per-criterion match scores, overall match, site proximity |
| 3 | Site Selection | `site_selection` | Site scoring (enrollment rate, diversity, screen failure), ranked list |
| 4 | Eligibility Optimization | `eligibility_optimization` | Population impact analysis, broaden/retain recommendations |
| 5 | Adaptive Design | `adaptive_design` | Design type selection, interim analysis plan, regulatory precedent |
| 6 | Safety Signal Detection | `safety_signal` | AE frequency analysis, PRR/ROR, severity classification, DSMB alerts |
| 7 | Regulatory Documents | `regulatory_docs` | IND/CSR/briefing/DSUR draft generation, agency-specific formatting |
| 8 | Competitive Intelligence | `competitive_intel` | Threat scoring, enrollment tracking, mechanism comparison |
| 9 | Diversity Assessment | `diversity_assessment` | Demographic gap analysis, site recommendations, FDA FDORA compliance |
| 10 | Decentralized Planning | `decentralized_planning` | DCT component assessment, regulatory feasibility, patient preference |

### Workflow-Specific Collection Weights

Each workflow applies custom search weights to prioritize the most relevant collections. For example, Patient Matching boosts `trial_eligibility` to 0.25, while Safety Signal Detection boosts `trial_safety` to 0.25. All 14 collections are searched for every workflow, but the weight distribution ensures domain-relevant evidence surfaces first.

---

## 5. Decision Support Engines

### 5.1 Confidence Calibrator

Multi-factor calibration model combining raw confidence (0.30), evidence base (0.30), document factor with logarithmic scaling (0.20), and cross-agent agreement (0.20). Evidence levels map from A1 (systematic review of RCTs, 1.0) through E (expert opinion, 0.15).

### 5.2 Protocol Complexity Scorer

Five-dimension scoring against Tufts CSDD industry benchmarks: procedure count (normalized by 50), visit frequency (normalized by 36), endpoint count (normalized by 20), eligibility criteria count (normalized by 40), and amendment count (normalized by 5). Weighted composite with percentile ranking.

### 5.3 Enrollment Predictor

Predicts monthly enrollment rate using: historical baseline rate, disease prevalence factor (rare=0.4 to common=1.2), competition factor (crowded=0.5 to no competition=1.2), capacity factor (staff-to-capacity ratio), phase factor (Phase I=0.6 to Phase IV=1.1), and eligibility stringency penalty (up to 40% reduction).

### 5.4 Eligibility Analyzer

Analyzes each eligibility criterion against 29 restrictive patterns with population exclusion estimates (e.g., "prior CAR-T" = 85%, "pregnancy or lactation" = 50%, "ECOG 0" = 40%). Assesses scientific justification strength and generates BROADEN/REVIEW/RETAIN recommendations with competitor comparison.

### 5.5 Competitive Threat Scorer

Four-factor threat model: phase advancement (0.30), enrollment progress (0.25), sponsor resources (0.20), and mechanism differentiation (0.25). Classifies threats as critical (>=0.80), high (>=0.60), moderate (>=0.40), low (>=0.20), or minimal.

### 5.6 Historical Success Rate Estimator

Phase-specific success rates across 12 therapeutic areas from BIO/QLS Advisors data, plus cumulative probability-of-success estimation from current phase to approval.

---

## 6. Cross-Agent Integration

The agent integrates with four other HCLS AI Factory intelligence agents:

| Agent | URL | Integration Purpose |
|---|---|---|
| Oncology Intelligence | localhost:8527 | Molecular trial matches for precision oncology |
| Pharmacogenomics Intelligence | localhost:8107 | Pharmacogenomic metabolism screening |
| Cardiology Intelligence | localhost:8126 | Cardiac safety assessment (QTc, MACE risk) |
| Biomarker Intelligence | localhost:8529 | Biomarker enrichment strategies |

All cross-agent queries degrade gracefully: if an agent is unavailable, a warning is logged and a default response is returned with a 30-second timeout.

---

## 7. Vector Database & Collections

### 14 Milvus Collections

| # | Collection | Description | Est. Records | Search Weight |
|---|---|---|---|---|
| 1 | `trial_protocols` | Protocol documents with phase, sponsor, enrollment metadata | 5,000 | 0.10 |
| 2 | `trial_eligibility` | Inclusion/exclusion criteria with population impact estimates | 50,000 | 0.09 |
| 3 | `trial_endpoints` | Endpoint definitions with measures, time frames, statistical methods | 20,000 | 0.08 |
| 4 | `trial_sites` | Site locations, enrollment counts, recruitment status | 30,000 | 0.07 |
| 5 | `trial_investigators` | PIs with publication metrics and therapeutic area expertise | 5,000 | 0.05 |
| 6 | `trial_results` | Published results with statistical outcomes and effect sizes | 3,000 | 0.09 |
| 7 | `trial_regulatory` | Regulatory submissions, decisions across agencies | 2,000 | 0.07 |
| 8 | `trial_literature` | Published clinical trial methodology and results literature | 10,000 | 0.08 |
| 9 | `trial_biomarkers` | Biomarker assays, thresholds, validation status | 3,000 | 0.07 |
| 10 | `trial_safety` | AE profiles with MedDRA coding, severity, frequency | 20,000 | 0.08 |
| 11 | `trial_rwe` | Real-world evidence from claims, EHR, patient registries | 2,000 | 0.06 |
| 12 | `trial_adaptive` | Adaptive designs, decision rules, regulatory precedents | 500 | 0.05 |
| 13 | `trial_guidelines` | ICH/FDA/EMA clinical trial guidelines and requirements | 1,000 | 0.08 |
| 14 | `genomic_evidence` | Shared genomic variant evidence (ClinVar, AlphaMissense) | 100,000 | 0.03 |

**Total estimated records:** ~251,500

### Index Configuration

- **Embedding model:** BAAI/bge-small-en-v1.5 (384 dimensions)
- **Index type:** IVF_FLAT
- **Metric type:** COSINE
- **nlist:** 128
- **Batch size:** 32

---

## 8. RAG Engine

The `TrialRAGEngine` implements multi-collection retrieval with:

- Weighted search across all 14 collections with workflow-specific boosting
- Score threshold filtering (default 0.4)
- Top-K per collection (default 5, configurable 1-50)
- Cross-collection result merging and re-ranking
- Citation generation with high/medium thresholds (0.75 / 0.60)
- Graceful degradation when Milvus is unavailable (returns knowledge-base-only results)

---

## 9. Query Expansion System

### 10 Synonym Maps

| Map | Entries | Purpose |
|---|---|---|
| Entity Aliases | 140 | Abbreviation to full term (NSCLC, HFrEF, Keytruda, TMB, etc.) |
| Therapeutic Area Map | 13 | Area keywords for classification (oncology, cardiovascular, etc.) |
| Phase Map | 7 | Phase synonym resolution (Phase III, pivotal, confirmatory, etc.) |
| Drug Synonym Map | 33 | Generic/brand/code name mapping (pembrolizumab/Keytruda/MK-3475) |
| Biomarker Map | 22 | Biomarker alias expansion (PD-L1, TMB, MSI, HER2, EGFR, etc.) |
| Endpoint Map | 15 | Endpoint type synonyms (OS, PFS, ORR, DFS, MACE, QoL, etc.) |
| Regulatory Map | 19 | Regulatory term expansion (IND, NDA, BLA, BTD, SPA, PDUFA, etc.) |
| Design Map | 14 | Trial design synonyms (adaptive, basket, umbrella, platform, etc.) |
| Population Map | 10 | Population descriptors (pediatric, geriatric, treatment-naive, etc.) |
| Safety Map | 10+ | Safety term expansion (SUSAR, SAE, DSMB, MedDRA, CRS, etc.) |

**Total:** ~283 top-level entries, 140 entity aliases, covering abbreviations across oncology, cardiology, neuroscience, metabolic, autoimmune, rare disease, and regulatory domains.

---

## 10. Autonomous Agent Pipeline

The `TrialIntelligenceAgent` implements the plan -> search -> evaluate -> synthesize -> report pattern:

1. **Plan** (`search_plan`): Analyze query, detect therapeutic areas, drugs, biomarkers; decompose into sub-questions; select search strategy (broad, targeted, comparative, regulatory)
2. **Search** (`rag_engine.query`): Execute multi-collection retrieval with workflow-specific weights
3. **Evaluate** (`evaluate_evidence`): Score evidence quality, check completeness, assess cross-agent agreement
4. **Synthesize**: Generate structured clinical trial recommendations with regulatory alerts
5. **Report** (`generate_report`): Format output with citations, guideline references, and confidence scores

Evidence levels range from 1a (systematic review of RCTs) through 5 (expert opinion) plus a regulatory guidance level.

---

## 11. Data Models & Type Safety

All data contracts are defined as Pydantic BaseModel or Python Enum classes in `src/models.py`:

### Enums (12)

| Enum | Values | Purpose |
|---|---|---|
| `TrialWorkflowType` | 19 values | Workflow type routing (including aliases) |
| `TrialPhase` | 7 values | Clinical trial phases (Phase I through Phase IV + N/A) |
| `TrialStatus` | 7 values | Recruitment status (recruiting through withdrawn) |
| `EvidenceLevel` | 6 values | Evidence classification (A1 through E) |
| `CriterionType` | 2 values | Inclusion / exclusion |
| `EndpointType` | 4 values | Primary, secondary, exploratory, safety |
| `RegulatoryAgency` | 6 values | FDA, EMA, PMDA, Health Canada, TGA, MHRA |
| `DocumentType` | 6 values | IND, CSR, briefing, PSP, RMP, DSUR |
| `SeverityLevel` | 5 values | Critical through informational |
| `TherapeuticArea` | 13 values | Oncology through other |
| `DCTComponent` | 7 values | eConsent through direct-to-patient |
| `SearchPlan` | dataclass | Pre-retrieval search planning |

### Pydantic Models (12)

`TrialQuery`, `TrialSearchResult`, `MatchScore`, `OverallMatch`, `PatientProfile`, `EligibilityAnalysis`, `SiteScore`, `SafetySignal`, `CompetitorProfile`, `ProtocolComplexity`, `WorkflowResult`, `TrialResponse`

---

## 12. Streamlit UI

### 5-Tab Interface (port 8128)

| Tab | Features |
|---|---|
| **Trial Intelligence** | RAG-powered Q&A with workflow-aware routing, citation display |
| **Patient Matching** | Patient profile input, per-criterion scoring, site proximity |
| **Protocol Optimizer** | Protocol complexity scoring, endpoint recommendations, design suggestions |
| **Competitive Landscape** | Competitor threat scoring, enrollment tracking, mechanism comparison |
| **Dashboard** | Collection health, query volume, workflow execution metrics |

### UI Theme

NVIDIA dark theme with custom CSS: primary background `#1a1a2e`, secondary `#16213e`, card `#0f3460`, accent NVIDIA green `#76b900`, with danger/warning/info/success color coding.

---

## 13. REST API

### 26 Endpoints (FastAPI, port 8538)

| Method | Path | Description |
|---|---|---|
| GET | `/health` | Service health with collection and vector counts |
| GET | `/collections` | Collection names and record counts |
| GET | `/workflows` | Available clinical trial workflows |
| GET | `/metrics` | Prometheus-compatible metrics |
| POST | `/v1/trial/query` | RAG Q&A query |
| POST | `/v1/trial/search` | Multi-collection vector search |
| POST | `/v1/trial/protocol/optimize` | Protocol optimization workflow |
| POST | `/v1/trial/match` | Patient-trial matching |
| POST | `/v1/trial/match/batch` | Batch patient matching |
| POST | `/v1/trial/site/recommend` | Site selection recommendations |
| POST | `/v1/trial/eligibility/optimize` | Eligibility criteria optimization |
| POST | `/v1/trial/adaptive/evaluate` | Adaptive design evaluation |
| POST | `/v1/trial/safety/signal` | Safety signal detection |
| POST | `/v1/trial/regulatory/generate` | Regulatory document generation |
| POST | `/v1/trial/competitive/landscape` | Competitive intelligence |
| POST | `/v1/trial/diversity/assess` | Diversity assessment |
| POST | `/v1/trial/dct/plan` | Decentralized trial planning |
| GET | `/v1/trial/therapeutic-areas` | Therapeutic area reference catalog |
| GET | `/v1/trial/phases` | Phase reference data |
| GET | `/v1/trial/guidelines` | Guideline reference data |
| GET | `/v1/trial/knowledge-version` | Knowledge base version metadata |
| POST | `/v1/trial/workflow/{type}` | Generic workflow dispatch |
| POST | `/v1/reports/generate` | Report generation |
| GET | `/v1/reports/formats` | Supported export formats |
| GET | `/v1/events/stream` | SSE event stream |
| GET | `/v1/events/health` | Event subsystem health |

---

## 14. Data Ingest Pipelines

Three specialized ingest parsers, all inheriting from `BaseIngestPipeline`:

| Parser | Source | Output Collections |
|---|---|---|
| `ClinicalTrialsParser` | ClinicalTrials.gov XML/JSON API | protocols, eligibility, endpoints, sites, investigators |
| `PubMedParser` | PubMed/MEDLINE E-utilities API | literature, results |
| `RegulatoryParser` | FDA/EMA regulatory documents | regulatory, guidelines, safety |

### Ingest Pipeline Features

- Chunking with configurable overlap for long documents
- Deduplication by trial ID and content hash
- Rate limiting for external API calls (NCBI API key support)
- Incremental ingest (new records only, timestamp-based)
- Validation and error reporting at each stage

---

## 15. Seed Data Inventory

The knowledge base is seeded from `scripts/seed_knowledge.py` with:

| Category | Count | Source |
|---|---|---|
| Landmark Trials | 40 | Curated from published Phase 3 results |
| Therapeutic Areas | 13 | ICH/FDA/BIO classification |
| Trial Phases | 7 | ICH E8(R1) |
| Regulatory Agencies | 9 | FDA, EMA, PMDA, Health Canada, TGA, MHRA, NMPA, Swissmedic, ANVISA |
| Endpoint Types | 9 | Primary through ctDNA clearance |
| Adaptive Designs | 9 | Group sequential through dose-finding |
| Biomarker Strategies | 9 | Enrichment through digital biomarker |
| DCT Components | 9 | eConsent through digital informed consent |
| Safety Signal Metrics | 6 | PRR, ROR, EBGM, etc. |
| Drug Synonyms | 33 | Brand/generic/code name mappings |
| Biomarker Aliases | 22 | Multi-assay and synonym mappings |
| Entity Aliases | 140 | Abbreviation resolution |
| Endpoint Synonyms | 15 | OS, PFS, ORR, MACE, QoL, PRO, etc. |
| Regulatory Terms | 19 | IND, NDA, BLA, BTD, SPA, etc. |
| Design Synonyms | 14 | Adaptive, basket, umbrella, platform, etc. |
| Population Descriptors | 10 | Pediatric, geriatric, treatment-naive, etc. |

### Knowledge Sources

- ICH E6(R2) / E6(R3) Good Clinical Practice
- ICH E8(R1) General Considerations for Clinical Studies
- ICH E9(R1) Statistical Principles -- Estimands
- FDA Guidance: Adaptive Designs for Clinical Trials (2019)
- FDA Guidance: Decentralized Clinical Trials (2023)
- EMA Guideline on Multiplicity Issues in Clinical Trials
- ClinicalTrials.gov registry data
- FDA CDER Drug Approval Reports 2015-2025
- BIO/QLS Advisors Clinical Development Success Rates 2011-2020
- PhRMA Clinical Trial Success Rates 2011-2020

---

## 16. Export & Reporting

The `src/export.py` module (725 LOC) supports:

- **Markdown** reports with structured headings, tables, and citations
- **JSON** structured data export for programmatic consumption
- **PDF** generation via report templates
- Trial-specific formatting: protocol summaries, patient match reports, safety signal reports, competitive landscape comparisons
- Citation formatting with collection source, relevance score, and evidence level

---

## 17. Observability & Metrics

Prometheus-compatible metrics via `src/metrics.py` (525 LOC):

| Category | Metrics |
|---|---|
| Query | `trial_queries_total`, `trial_query_duration_seconds`, `trial_query_errors_total` |
| RAG/Search | `trial_search_total`, `trial_search_duration_seconds` |
| Workflow | `trial_workflow_executions_total`, `trial_workflow_duration_seconds` |
| LLM | `trial_llm_tokens_total`, `trial_llm_latency_seconds` |
| Ingest | `trial_ingest_records_total`, `trial_ingest_errors_total` |
| System | `trial_collection_records_gauge`, `trial_system_info` |

All metrics gracefully degrade if `prometheus_client` is not installed (no-op stubs).

---

## 18. Scheduling & Automation

The `src/scheduler.py` (577 LOC) provides:

- Configurable ingest schedule (default 24-hour interval)
- Collection maintenance (compaction, index rebuild)
- Knowledge base refresh from external sources
- Daemon thread execution with graceful shutdown
- Enable/disable via `TRIAL_INGEST_ENABLED` environment variable

---

## 19. Configuration System

Pydantic BaseSettings in `config/settings.py`:

| Section | Key Settings |
|---|---|
| Paths | `PROJECT_ROOT`, `DATA_DIR`, `CACHE_DIR`, `REFERENCE_DIR` |
| Milvus | `MILVUS_HOST=localhost`, `MILVUS_PORT=19530` |
| Embeddings | `EMBEDDING_MODEL=BAAI/bge-small-en-v1.5`, `EMBEDDING_DIMENSION=384` |
| LLM | `LLM_MODEL=claude-sonnet-4-6`, `ANTHROPIC_API_KEY` |
| Search | `TOP_K_PER_COLLECTION=5`, `SCORE_THRESHOLD=0.4` |
| Weights | 14 collection weights summing to ~1.0 |
| API | `API_PORT=8538`, `STREAMLIT_PORT=8128` |
| Cross-Agent | URLs for oncology (8527), PGx (8107), cardiology (8126), biomarker (8529) |
| Security | `API_KEY` (empty = no auth), `CORS_ORIGINS` |
| Scheduler | `INGEST_SCHEDULE_HOURS=24`, `INGEST_ENABLED=False` |

Environment prefix: `TRIAL_` (e.g., `TRIAL_MILVUS_HOST`, `TRIAL_API_PORT`).

Startup validation checks Milvus connectivity, API key presence, port conflicts, and weight sum (~1.0 tolerance 0.05).

---

## 20. Security & Authentication

| Feature | Implementation |
|---|---|
| API Key Authentication | `X-API-Key` header validation (optional, empty = disabled) |
| CORS | Configurable origins (default: localhost:8080, 8538, 8128) |
| Rate Limiting | 100 requests/minute per IP |
| Input Validation | Pydantic models with field constraints (min/max length, ranges) |
| SQL Injection | N/A (no SQL -- Milvus vector search only) |
| Path Traversal | Controlled via Pydantic BaseSettings path validation |
| Secrets Management | Environment variables, `.env` file support |
| HTTPS | Handled by reverse proxy (nginx/traefik in production) |

---

## 21. Infrastructure & Deployment

### Integrated Deployment (docker-compose)

The agent runs as part of the HCLS AI Factory `docker-compose.dgx-spark.yml`:

```yaml
clinical-trial-agent-api:
  build: ./ai_agent_adds/clinical_trial_intelligence_agent
  ports:
    - "8538:8538"
  environment:
    - TRIAL_MILVUS_HOST=milvus-standalone
    - TRIAL_MILVUS_PORT=19530
    - TRIAL_ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY}

clinical-trial-agent-ui:
  command: streamlit run app/trial_ui.py --server.port 8128
  ports:
    - "8128:8128"
```

### Standalone Deployment

```bash
# Start Milvus standalone
docker-compose -f docker-compose.yml up -d milvus-standalone

# API server
uvicorn api.main:app --host 0.0.0.0 --port 8538

# Streamlit UI
streamlit run app/trial_ui.py --server.port 8128
```

### Port Map

| Port | Service | Protocol |
|---|---|---|
| 8538 | FastAPI REST API | HTTP |
| 8128 | Streamlit UI | HTTP |
| 19530 | Milvus (shared) | gRPC |
| 2379 | etcd (Milvus metadata) | gRPC |
| 9000 | MinIO (Milvus blob storage) | HTTP |

---

## 22. Test Suite

### Results

```
769 passed in 0.47s
```

| Metric | Value |
|---|---|
| Total Tests | 769 |
| Passed | 769 |
| Failed | 0 |
| Pass Rate | 100% |
| Execution Time | 0.47s |
| Test Files | 12 |

### Test Coverage by Module

| Test File | Tests | Module Covered |
|---|---|---|
| `test_models.py` | 519 LOC | Pydantic models, enums, validation |
| `test_workflow_execution.py` | 379 LOC | All 10 workflow execute methods |
| `test_clinical_workflows.py` | 347 LOC | Workflow preprocess/postprocess |
| `test_api.py` | 310 LOC | All 26 API endpoints |
| `test_agent.py` | 294 LOC | Agent pipeline (plan, evaluate, synthesize) |
| `test_decision_support.py` | 281 LOC | All 5 decision support engines |
| `test_query_expansion.py` | 255 LOC | Synonym maps, alias resolution, boosting |
| `test_integration.py` | 247 LOC | End-to-end workflow + RAG integration |
| `test_rag_engine.py` | 194 LOC | Multi-collection retrieval, scoring |
| `test_collections.py` | 152 LOC | Schema validation, weight sums |
| `test_settings.py` | 133 LOC | Configuration validation |
| `test_knowledge.py` | 123 LOC | Knowledge base completeness and structure |

---

## 23. Known Limitations

| # | Limitation | Impact | Mitigation |
|---|---|---|---|
| 1 | No real-time ClinicalTrials.gov API integration | Data is point-in-time snapshots | Scheduled ingest pipeline (24h interval) |
| 2 | LLM required for natural language synthesis | Search-only mode without API key | Knowledge base and workflow results still available |
| 3 | Single-language (English) support | Non-English trial data not indexed | ICH guidelines are English; future: multi-language |
| 4 | No FHIR R4 export | Cannot export to EHR systems | Structured JSON export available |
| 5 | Enrollment predictions are model-based | Actual enrollment varies significantly | Calibrated with Tufts CSDD benchmarks |
| 6 | Competitive intelligence is point-in-time | Competitor status changes rapidly | Scheduled refresh via ingest pipeline |
| 7 | Cross-agent integration requires all agents running | Partial functionality without peer agents | Graceful degradation with default responses |
| 8 | No built-in authentication rotation | API key is static | External secret management recommended |
| 9 | Milvus standalone mode for demo | Not production HA configuration | Milvus cluster mode for production |
| 10 | Safety signal detection is statistical only | No causal inference | PRR/ROR as screening tools; expert review required |

---

## 24. Demo Readiness Audit

| # | Check | Status |
|---|---|---|
| 1 | All 769 tests pass | PASS |
| 2 | API starts without errors | PASS |
| 3 | Streamlit UI loads with NVIDIA theme | PASS |
| 4 | Health endpoint returns 200 | PASS |
| 5 | Knowledge base loaded (40 trials, 13 areas) | PASS |
| 6 | All 10 workflows execute without Milvus | PASS |
| 7 | All 5 decision support engines functional | PASS |
| 8 | Query expansion resolves common abbreviations | PASS |
| 9 | Protocol design returns blueprint for NSCLC Phase 3 | PASS |
| 10 | Patient matching scores criteria correctly | PASS |
| 11 | Site selection ranks by composite score | PASS |
| 12 | Eligibility analysis identifies restrictive criteria | PASS |
| 13 | Adaptive design recommends appropriate designs | PASS |
| 14 | Safety signal detects elevated AE frequencies | PASS |
| 15 | Competitive threat scoring classifies correctly | PASS |
| 16 | Diversity assessment identifies geographic gaps | PASS |
| 17 | DCT planning evaluates component feasibility | PASS |
| 18 | Cross-agent integration degrades gracefully | PASS |
| 19 | Export generates Markdown and JSON | PASS |
| 20 | Prometheus metrics endpoint returns data | PASS |
| 21 | CORS headers present in API responses | PASS |
| 22 | Configuration validation detects issues | PASS |
| 23 | Collection schemas create successfully | PASS |
| 24 | Batch patient matching handles multiple patients | PASS |
| 25 | Historical success rates return for all 12 areas | PASS |

---

## 25. Codebase Summary

| Category | Files | LOC |
|---|---|---|
| Source (`src/`) | 17 | ~13,200 |
| API (`api/`) | 4 | ~2,275 |
| App (`app/`) | 2 | ~692 |
| Config (`config/`) | 2 | ~182 |
| Scripts (`scripts/`) | 3 | ~747 |
| Tests (`tests/`) | 12 | ~3,904 |
| Infrastructure | 3 | ~1,607 |
| **Total** | **46** | **22,607** |

### File Count by Type

- Python source files: 34 (excluding tests)
- Python test files: 12
- Configuration files: 3 (docker-compose.yml, Dockerfile, requirements.txt)
- Documentation: 2 (README.md, research paper)
- **Total Python files: 46**
- **Total Python LOC: 22,607**

---

*Generated: March 22, 2026 | Clinical Trial Intelligence Agent v2.0.0 | HCLS AI Factory*
