# Clinical Trial Intelligence Agent -- Project Bible

**Version:** 2.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [Overview](#1-overview)
2. [Architecture](#2-architecture)
3. [Collections Reference](#3-collections-reference)
4. [Workflow Reference](#4-workflow-reference)
5. [API Endpoint Reference](#5-api-endpoint-reference)
6. [Knowledge Base Reference](#6-knowledge-base-reference)
7. [Decision Support Engines](#7-decision-support-engines)
8. [Query Expansion Reference](#8-query-expansion-reference)
9. [Configuration Reference](#9-configuration-reference)
10. [Port Map](#10-port-map)
11. [Tech Stack](#11-tech-stack)
12. [Data Models](#12-data-models)
13. [Cross-Agent Integration](#13-cross-agent-integration)
14. [Ingest Pipeline Reference](#14-ingest-pipeline-reference)
15. [Test Reference](#15-test-reference)

---

## 1. Overview

### 1.1 Platform Positioning

The Clinical Trial Intelligence Agent is one of **11 intelligence agents** in the HCLS AI Factory, a three-engine precision medicine platform (Genomics, RAG/Chat, Drug Discovery) running on NVIDIA DGX Spark. It occupies the clinical trial intelligence niche, connecting genomic biomarker data from the Genomics Engine with trial design optimization and patient-trial matching.

**All 11 HCLS AI Factory Intelligence Agents:**

| # | Agent | Port | Focus |
|---|---|---|---|
| 1 | Biomarker Intelligence | :8529 | Biomarker discovery and stratification |
| 2 | Oncology Intelligence | :8527/:8528 | Cancer genomics and treatment |
| 3 | CAR-T Intelligence | -- | CAR-T cell therapy development |
| 4 | Imaging Intelligence | :8524 | Medical imaging AI |
| 5 | Autoimmune Intelligence | -- | Autoimmune disease genomics |
| 6 | Pharmacogenomics Intelligence | :8107 | Drug metabolism and dosing |
| 7 | **Clinical Trial Intelligence** | **:8538** | **Trial design and matching (this agent)** |
| 8 | Rare Disease Diagnostic | :8134 | Rare disease diagnosis |
| 9 | Single-Cell Intelligence | :8540 | Single-cell transcriptomics |
| 10 | Cardiology Intelligence | :8126 | Cardiac genetics |
| 11 | Neurology Intelligence | -- | Neurological genetics |

### 1.2 Agent Summary

The Clinical Trial Intelligence Agent is an AI-powered clinical trial decision support system that integrates RAG-based evidence retrieval across 14 Milvus vector collections, 10 clinical workflows, 5 decision support engines, and an autonomous reasoning pipeline. It serves pharmaceutical R&D teams, clinical operations, regulatory affairs, and medical affairs with evidence-based guidance across the entire clinical trial lifecycle.

### Key Numbers

| Metric | Value |
|---|---|
| Python files | 46 |
| Lines of code | 22,607 |
| Milvus collections | 14 |
| Clinical workflows | 10 |
| Decision support engines | 5 + 1 (Historical Success Estimator) |
| API endpoints | 26 |
| Landmark trials | 40 |
| Therapeutic areas | 13 |
| Regulatory agencies | 9 |
| Entity aliases | 140 |
| Drug synonym entries | 33 |
| Biomarker entries | 22 |
| Tests | 769 (100% pass, 0.47s) |
| Knowledge version | 2.0.0 |

---

## 2. Architecture

```
User --> Streamlit UI (:8128) --> FastAPI API (:8538)
                                       |
                      +----------------+----------------+
                      |                |                |
               Workflows(10)    Decision Engines(5)  RAG Engine
                      |                                |
                 Knowledge Base                  Milvus (:19530)
               (40 trials, 13 areas,             14 collections
                9 agencies, 9 designs)           384-dim BGE
                                                 IVF_FLAT/COSINE
```

### Tiers

- **Presentation:** Streamlit (5 tabs, NVIDIA dark theme, port 8128)
- **Application:** FastAPI (26 endpoints, CORS, auth, rate limiting, port 8538)
- **Data:** Milvus (14 collections, BGE-small-en-v1.5 embeddings, port 19530)

---

## 3. Collections Reference

### 3.1 Full Collection Catalog

| # | Name | Fields | Est. Records | Weight | Primary Use |
|---|---|---|---|---|---|
| 1 | `trial_protocols` | trial_id, title, phase, status, therapeutic_area, sponsor, start_date, enrollment_target, text_content | 5,000 | 0.10 | Protocol design, competitive intel |
| 2 | `trial_eligibility` | trial_id, criterion_type, criterion_text, logic_operator, population_impact | 50,000 | 0.09 | Patient matching, eligibility optimization |
| 3 | `trial_endpoints` | trial_id, endpoint_type, measure, time_frame, statistical_method | 20,000 | 0.08 | Protocol design, adaptive design |
| 4 | `trial_sites` | trial_id, site_id, facility_name, city, state, country, status, enrollment_count | 30,000 | 0.07 | Site selection, diversity assessment |
| 5 | `trial_investigators` | investigator_id, name, specialty, h_index, publication_count, therapeutic_areas | 5,000 | 0.05 | Site selection, competitive intel |
| 6 | `trial_results` | trial_id, outcome, p_value, effect_size, confidence_interval, publication_pmid | 3,000 | 0.09 | Protocol design, competitive intel |
| 7 | `trial_regulatory` | submission_id, agency, decision, document_type, drug_name, indication | 2,000 | 0.07 | Regulatory docs, competitive intel |
| 8 | `trial_literature` | pmid, title, journal, mesh_terms, publication_year, study_type | 10,000 | 0.08 | Evidence synthesis, protocol design |
| 9 | `trial_biomarkers` | biomarker, assay, threshold, validated, trial_context | 3,000 | 0.07 | Patient matching, biomarker strategy |
| 10 | `trial_safety` | trial_id, event_type, severity, frequency, soc_term | 20,000 | 0.08 | Safety signal, regulatory docs |
| 11 | `trial_rwe` | source, population, outcome, study_design, sample_size | 2,000 | 0.06 | Eligibility optimization, diversity |
| 12 | `trial_adaptive` | design_type, decision_rule, trigger_criteria, historical_precedent | 500 | 0.05 | Adaptive design evaluation |
| 13 | `trial_guidelines` | guideline_id, organization, version, recommendation_text, evidence_class | 1,000 | 0.08 | All workflows (regulatory reference) |
| 14 | `genomic_evidence` | gene, variant, clinical_significance, condition, evidence_summary, source | 100,000 | 0.03 | Cross-modal genomic queries |

### 3.2 Index Configuration

All collections use identical index parameters:

```python
EMBEDDING_DIM = 384       # BGE-small-en-v1.5
INDEX_TYPE = "IVF_FLAT"
METRIC_TYPE = "COSINE"
NLIST = 128
```

### 3.3 Workflow-Specific Weight Overrides

Each workflow boosts its primary collections. Example weights for the top-3 collections per workflow:

| Workflow | #1 Collection (Weight) | #2 Collection (Weight) | #3 Collection (Weight) |
|---|---|---|---|
| Protocol Design | trial_protocols (0.20) | trial_endpoints (0.15) | trial_eligibility (0.12) |
| Patient Matching | trial_eligibility (0.25) | trial_protocols (0.15) | trial_sites (0.12) |
| Site Selection | trial_sites (0.25) | trial_investigators (0.18) | trial_protocols (0.10) |
| Eligibility Optimization | trial_eligibility (0.25) | trial_protocols (0.12) | trial_rwe (0.10) |
| Adaptive Design | trial_adaptive (0.25) | trial_endpoints (0.15) | trial_guidelines (0.12) |
| Safety Signal | trial_safety (0.25) | trial_results (0.15) | trial_protocols (0.10) |
| Regulatory Docs | trial_regulatory (0.25) | trial_guidelines (0.18) | trial_results (0.12) |
| Competitive Intel | trial_protocols (0.20) | trial_results (0.15) | trial_endpoints (0.12) |
| Diversity Assessment | trial_sites (0.22) | trial_eligibility (0.18) | trial_protocols (0.12) |
| Decentralized Planning | trial_sites (0.18) | trial_protocols (0.15) | trial_guidelines (0.12) |
| General | trial_protocols (0.10) | trial_results (0.09) | trial_eligibility (0.09) |

---

## 4. Workflow Reference

### 4.1 Workflow Catalog

| Workflow | Enum Value | Key Inputs | Key Outputs |
|---|---|---|---|
| Protocol Design | `protocol_design` | indication, phase, comparator, mechanism | Protocol blueprint, endpoint recommendations, sample size |
| Patient Matching | `patient_matching` | PatientProfile (age, dx, biomarkers, variants) | MatchScore per criterion, OverallMatch, site distances |
| Site Selection | `site_selection` | therapeutic_area, phase, target_enrollment | SiteScore list, enrollment rate, diversity index |
| Eligibility Optimization | `eligibility_optimization` | criteria list with justifications | Population impact, BROADEN/REVIEW/RETAIN recs |
| Adaptive Design | `adaptive_design` | trial parameters, uncertainty profile | Design type, interim analysis schedule, regulatory guidance |
| Safety Signal | `safety_signal` | AE data, trial_id | SafetySignal list, PRR/ROR, severity, DSMB alerts |
| Regulatory Docs | `regulatory_docs` | trial data, agency, document_type | Structured document draft, agency-specific formatting |
| Competitive Intel | `competitive_intel` | indication, mechanism, target trial | CompetitorProfile list, threat scores, timeline comparison |
| Diversity Assessment | `diversity_assessment` | trial sites, target demographics | Demographic gaps, site recommendations, FDORA compliance |
| Decentralized Planning | `decentralized_planning` | trial type, patient population | DCT component assessment, regulatory feasibility |

### 4.2 Workflow Contract

All workflows inherit from `BaseTrialWorkflow`:

```python
class BaseTrialWorkflow(ABC):
    workflow_type: TrialWorkflowType

    def run(self, inputs: dict) -> WorkflowResult:
        processed = self.preprocess(inputs)
        result = self.execute(processed)
        return self.postprocess(result)

    def preprocess(self, inputs: dict) -> dict: ...
    @abstractmethod
    def execute(self, inputs: dict) -> WorkflowResult: ...
    def postprocess(self, result: WorkflowResult) -> WorkflowResult: ...
```

---

## 5. API Endpoint Reference

### Base URL: `http://localhost:8538`

### 5.1 System Endpoints

| Method | Path | Response | Auth Required |
|---|---|---|---|
| GET | `/health` | `{"status": "ok", "collections": {...}}` | No |
| GET | `/collections` | Collection names and record counts | No |
| GET | `/workflows` | Available workflow types | No |
| GET | `/metrics` | Prometheus text format | No |

### 5.2 Trial Endpoints (`/v1/trial/`)

| Method | Path | Request Body | Response Model |
|---|---|---|---|
| POST | `/query` | `{"question": "...", "workflow_type": "..."}` | QueryResponse |
| POST | `/search` | `{"query": "...", "collections": [...], "top_k": 5}` | SearchResponse |
| POST | `/protocol/optimize` | `{"indication": "...", "phase": "..."}` | ProtocolOptimizeResponse |
| POST | `/match` | `{"patient": {...}, "trial_ids": [...]}` | PatientMatchResponse |
| POST | `/match/batch` | `{"patients": [...]}` | BatchMatchResponse |
| POST | `/site/recommend` | `{"therapeutic_area": "...", "phase": "..."}` | SiteRecommendResponse |
| POST | `/eligibility/optimize` | `{"criteria": [...]}` | EligibilityOptimizeResponse |
| POST | `/adaptive/evaluate` | `{"trial_params": {...}}` | AdaptiveEvaluateResponse |
| POST | `/safety/signal` | `{"events": [...], "trial_id": "..."}` | SafetySignalResponse |
| POST | `/regulatory/generate` | `{"trial_data": {...}, "doc_type": "..."}` | RegulatoryGenerateResponse |
| POST | `/competitive/landscape` | `{"indication": "...", "mechanism": "..."}` | CompetitiveLandscapeResponse |
| POST | `/diversity/assess` | `{"sites": [...], "demographics": {...}}` | DiversityAssessResponse |
| POST | `/dct/plan` | `{"trial_type": "...", "population": "..."}` | DCTPlanResponse |
| GET | `/therapeutic-areas` | -- | Therapeutic area catalog |
| GET | `/phases` | -- | Phase reference |
| GET | `/guidelines` | -- | Guideline reference |
| GET | `/knowledge-version` | -- | Version metadata |
| POST | `/workflow/{type}` | `{"inputs": {...}}` | WorkflowResponse |

### 5.3 Report and Event Endpoints

| Method | Path | Description |
|---|---|---|
| POST | `/v1/reports/generate` | Generate structured report |
| GET | `/v1/reports/formats` | List supported export formats |
| GET | `/v1/events/stream` | SSE event stream |
| GET | `/v1/events/health` | Event subsystem health |

---

## 6. Knowledge Base Reference

### 6.1 Therapeutic Areas (13)

oncology, cardiovascular, neuroscience, immunology, infectious_disease, rare_diseases, metabolic, respiratory, hematology, gastroenterology, dermatology, ophthalmology, gene_cell_therapy

### 6.2 Trial Phases (7)

preclinical, phase_0, phase_1, phase_2, phase_3, phase_4, expanded_access

### 6.3 Regulatory Agencies (9)

FDA, EMA, PMDA, Health_Canada, TGA, MHRA, NMPA, Swissmedic, ANVISA

### 6.4 Endpoint Types (9)

primary, secondary, exploratory, safety, patient_reported, digital, composite, minimal_residual_disease, ctDNA_clearance

### 6.5 Adaptive Designs (9)

group_sequential, sample_size_reestimation, response_adaptive, biomarker_adaptive, platform_trial, seamless_phase, master_protocol, enrichment_adaptive, dose_finding

### 6.6 Biomarker Strategies (9)

enrichment, stratification, prognostic, predictive, pharmacodynamic, surrogate, companion_diagnostic, liquid_biopsy, digital_biomarker

### 6.7 DCT Components (9)

econsent, telemedicine, home_health, local_labs, wearables, epro_ecoa, direct_to_patient, remote_monitoring, digital_informed_consent

### 6.8 Landmark Trials (40)

KEYNOTE-024, EMPEROR-Reduced, RECOVERY, PARADIGM-HF, CheckMate-067, SPRINT, EMPA-REG_OUTCOME, DAPA-CKD, HIMALAYA, DESTINY-Breast04, ADVANCE, CLARITY-AD, FOURIER, FLAURA, ADAURA, MAGELLAN, I-SPY_2, VICTORIA, TOPAZ-1, KRYSTAL-1, ELARA, CREST, CheckMate-227, KEYNOTE-522, DESTINY-Lung02, CLEAR_Outcomes, SELECT, STEP_HFpEF, TRAILBLAZER-ALZ_2, EMERGE_ENGAGE, SPRINT_SMA, SUNFISH, CASGEVY, RINVOQ_SELECT, SURMOUNT-1, PURPOSE_1, PANORAMIC, EPIC-HR, IMpower110, DAPA-HF

---

## 7. Decision Support Engines

### 7.1 Engine Catalog

| Engine | Class | Purpose | Key Factors |
|---|---|---|---|
| Confidence Calibrator | `ConfidenceCalibrator` | Calibrate raw confidence scores | raw (0.3), evidence (0.3), docs (0.2), agreement (0.2) |
| Protocol Complexity | `ProtocolComplexityScorer` | Score protocol complexity | procedures, visits, endpoints, criteria, amendments |
| Enrollment Predictor | `EnrollmentPredictor` | Predict monthly enrollment | historical rate, prevalence, competition, capacity, phase |
| Eligibility Analyzer | `EligibilityAnalyzer` | Analyze criteria restrictiveness | 29 population impact patterns, justification scoring |
| Competitive Threat | `CompetitiveThreatScorer` | Score competitor threat level | phase (0.3), enrollment (0.25), sponsor (0.2), differentiation (0.25) |
| Success Estimator | `HistoricalSuccessEstimator` | Phase-specific success probability | 12 therapeutic areas, cumulative POS calculation |

### 7.2 Evidence Level Scoring

| Level | Score | Description |
|---|---|---|
| A1 | 1.00 | Systematic review of RCTs |
| A2 | 0.85 | High-quality RCT |
| B | 0.65 | Non-randomized controlled study |
| C | 0.45 | Observational study |
| D | 0.25 | Case series / case report |
| E | 0.15 | Expert opinion |

---

## 8. Query Expansion Reference

### 8.1 Synonym Maps

| Map | Key Count | Example Entry |
|---|---|---|
| ENTITY_ALIASES | 140 | `"NSCLC" -> "non-small cell lung cancer"` |
| THERAPEUTIC_AREA_MAP | 13 | `"oncology" -> ["cancer", "tumor", "neoplasm", ...]` |
| PHASE_MAP | 7 | `"phase 3" -> ["phase III", "pivotal", "confirmatory", ...]` |
| DRUG_SYNONYM_MAP | 33 | `"pembrolizumab" -> ["Keytruda", "MK-3475", ...]` |
| BIOMARKER_MAP | 22 | `"PD-L1" -> ["CD274", "TPS", "CPS", "SP263", ...]` |
| ENDPOINT_MAP | 15 | `"OS" -> ["overall survival", "mortality", ...]` |
| REGULATORY_MAP | 19 | `"IND" -> ["investigational new drug", ...]` |
| DESIGN_MAP | 14 | `"adaptive" -> ["interim analysis", "Bayesian adaptive", ...]` |
| POPULATION_MAP | 10 | `"pediatric" -> ["children", "adolescent", "PREA", ...]` |
| SAFETY_MAP | 10+ | `"SAE" -> ["serious adverse event", "hospitalization", ...]` |

### 8.2 QueryExpander Pipeline

1. Resolve entity aliases (NSCLC -> non-small cell lung cancer)
2. Detect therapeutic areas from query text
3. Expand drug names (brand -> generic, code names)
4. Expand biomarker references (PD-L1 -> all assay names)
5. Expand endpoint and regulatory terms
6. Apply workflow-aware term boosting

---

## 9. Configuration Reference

### Environment Variables

All variables use the `TRIAL_` prefix:

| Variable | Default | Description |
|---|---|---|
| `TRIAL_MILVUS_HOST` | localhost | Milvus server hostname |
| `TRIAL_MILVUS_PORT` | 19530 | Milvus server port |
| `TRIAL_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model name |
| `TRIAL_EMBEDDING_DIMENSION` | 384 | Embedding vector dimension |
| `TRIAL_LLM_MODEL` | claude-sonnet-4-6 | LLM model for synthesis |
| `TRIAL_ANTHROPIC_API_KEY` | (none) | Anthropic API key |
| `TRIAL_API_PORT` | 8538 | FastAPI server port |
| `TRIAL_STREAMLIT_PORT` | 8128 | Streamlit UI port |
| `TRIAL_API_KEY` | (empty) | API authentication key |
| `TRIAL_CORS_ORIGINS` | localhost:8080,8538,8128 | Allowed CORS origins |
| `TRIAL_TOP_K_PER_COLLECTION` | 5 | Results per collection |
| `TRIAL_SCORE_THRESHOLD` | 0.4 | Minimum similarity score |
| `TRIAL_INGEST_SCHEDULE_HOURS` | 24 | Ingest interval |
| `TRIAL_INGEST_ENABLED` | False | Enable scheduled ingest |
| `TRIAL_CROSS_AGENT_TIMEOUT` | 30 | Cross-agent query timeout (s) |

---

## 10. Port Map

| Port | Service | Protocol | Notes |
|---|---|---|---|
| 8538 | FastAPI REST API | HTTP | Clinical Trial Intelligence Agent |
| 8128 | Streamlit UI | HTTP | 5-tab clinical trial interface |
| 19530 | Milvus | gRPC | Shared vector store |
| 2379 | etcd | gRPC | Milvus metadata |
| 9000 | MinIO | HTTP | Milvus blob storage |
| 8527 | Oncology Agent | HTTP | Cross-agent integration |
| 8107 | PGx Agent | HTTP | Cross-agent integration |
| 8126 | Cardiology Agent | HTTP | Cross-agent integration |
| 8529 | Biomarker Agent | HTTP | Cross-agent integration |
| 8080 | Landing Page | HTTP | HCLS AI Factory hub |

---

## 11. Tech Stack

| Layer | Technology | Version/Details |
|---|---|---|
| Compute | NVIDIA DGX Spark | CUDA 12.x |
| LLM | Claude (Anthropic) | claude-sonnet-4-6 |
| Vector DB | Milvus | 2.x with etcd + MinIO |
| Embeddings | BGE-small-en-v1.5 | 384 dimensions, sentence-transformers |
| API Framework | FastAPI | Uvicorn ASGI server |
| UI Framework | Streamlit | NVIDIA dark theme |
| Data Models | Pydantic v2 | BaseModel + BaseSettings |
| Config | pydantic-settings | .env file support |
| Metrics | prometheus_client | Counter, Histogram, Gauge, Info |
| Container | Docker | Multi-stage build |
| Orchestration | docker-compose | DGX Spark stack |
| Testing | pytest | 769 tests, 0.47s |
| Python | 3.10+ | Type hints, dataclasses |

---

## 12. Data Models

### 12.1 Enums

| Enum | Count | Members |
|---|---|---|
| TrialWorkflowType | 19 | protocol_design, patient_matching, site_selection, eligibility_optimization, eligibility_analysis, endpoint_strategy, adaptive_design, safety_signal, safety_monitoring, regulatory_docs, regulatory_strategy, competitive_intel, competitive_intelligence, biomarker_strategy, rwe_analysis, recruitment_optimization, diversity_assessment, decentralized_planning, general |
| TrialPhase | 7 | phase_i through not_applicable |
| TrialStatus | 7 | recruiting through not_yet_recruiting |
| EvidenceLevel | 6 | a1 through e |
| TherapeuticArea | 13 | oncology through other |
| SeverityLevel | 5 | critical through informational |
| RegulatoryAgency | 6 | fda through mhra |
| DocumentType | 6 | ind through dsur |
| CriterionType | 2 | inclusion, exclusion |
| EndpointType | 4 | primary through safety |
| DCTComponent | 7 | econsent through direct_to_patient |

### 12.2 Pydantic Models

| Model | Fields | Purpose |
|---|---|---|
| TrialQuery | question, workflow_type, patient_context, top_k, include_guidelines | Input query |
| TrialSearchResult | collection, content, score, metadata | Single search result |
| PatientProfile | age, sex, diagnosis, biomarkers, medications, genomic_variants, comorbidities, geographic_location | Patient data |
| MatchScore | criterion_text, criterion_type, met, confidence, evidence | Per-criterion match |
| OverallMatch | trial_id, title, phase, status, inclusion_met/total, exclusion_clear/total, overall_score, confidence | Overall match |
| EligibilityAnalysis | criterion, population_impact, justification_score, competitor_comparison, recommendation | Criterion analysis |
| SiteScore | site_id, facility_name, city, country, enrollment_rate, screen_failure_rate, diversity_index, overall_score | Site evaluation |
| SafetySignal | event_type, severity, frequency, prr, ror, causality_assessment | Safety signal |
| CompetitorProfile | trial_id, sponsor, phase, indication, mechanism, enrollment_progress, estimated_completion, threat_level | Competitor data |
| ProtocolComplexity | procedure_count, visit_count, endpoint_count, eligibility_criteria_count, complexity_score, percentile_rank | Complexity assessment |
| WorkflowResult | workflow_type, findings, recommendations, guideline_references, severity, cross_agent_triggers, confidence | Workflow output |
| TrialResponse | answer, citations, workflow_results, matches, confidence | Top-level response |

---

## 13. Cross-Agent Integration

### 13.1 Direct Cross-Agent Endpoints (8 Agents)

| Agent | Endpoint | Purpose | Timeout |
|---|---|---|---|
| Oncology | `http://localhost:8527` | Molecular trial matches, biomarker-driven eligibility | 30s |
| PGx | `http://localhost:8107` | PGx-stratified enrollment, metabolizer phenotype data | 30s |
| Cardiology | `http://localhost:8126` | Cardiac safety endpoints, QTc risk assessment | 30s |
| Biomarker | `http://localhost:8529` | Enrichment strategies, companion diagnostic requirements | 30s |
| Rare Disease | `http://localhost:8134` | Orphan drug trial eligibility, gene therapy trials | 30s |
| Neurology | -- | CNS endpoint design, BBB penetration considerations | 30s |
| Single-Cell | `http://localhost:8540` | Cellular biomarker endpoints, MRD monitoring criteria | 30s |
| Imaging | `http://localhost:8524` | Imaging response criteria (RECIST, RANO), imaging eligibility | 30s |

All queries use graceful degradation: unavailable agents return default responses with warnings.

### 13.2 Pediatric Oncology Trial Focus

| Trial/Program | Identifier | Pediatric Context |
|---|---|---|
| COG AALL0932 | NCT00137111 | Standard-risk B-ALL; risk-adapted therapy de-escalation |
| COG AALL1231 | NCT02112916 | High-risk T-ALL; augmented BFM backbone with nelarabine |
| COG ANBL1232 | NCT01704716 | High-risk neuroblastoma; anti-GD2 immunotherapy (dinutuximab) |
| Pediatric MATCH | APEC1621 | Precision medicine basket trial for pediatric solid tumors and lymphomas |

**Pediatric-Specific Trial Design Features:**
- **Age Filters:** Neonatal (0-28d), infant (1-12mo), child (1-12yr), adolescent (12-18yr)
- **Rolling-6 Design:** COG Phase I dose escalation with continuous enrollment and safety monitoring
- **EFS Endpoints:** Event-free survival with pediatric-specific event definitions (relapse, progression, second malignancy, death)
- **Consent/Assent:** Parental consent for all minors; written assent >= 7 years; adolescent assent 12-17 years

---

## 14. Ingest Pipeline Reference

### 14.1 Parsers

| Parser | Source | API | Output Collections |
|---|---|---|---|
| ClinicalTrialsParser | ClinicalTrials.gov | XML/JSON REST API | protocols, eligibility, endpoints, sites, investigators |
| PubMedParser | PubMed/MEDLINE | E-utilities API | literature, results |
| RegulatoryParser | FDA/EMA docs | Document parsing | regulatory, guidelines, safety |

### 14.2 Pipeline Features

- Chunking with configurable overlap
- Content hash deduplication
- NCBI API key support for rate limit increase
- Incremental ingest (timestamp-based)
- Error reporting and validation

---

## 15. Test Reference

### 15.1 Test Files

| File | LOC | Coverage |
|---|---|---|
| test_models.py | 519 | All enums and Pydantic models |
| test_workflow_execution.py | 379 | All 10 workflow execute methods |
| test_clinical_workflows.py | 347 | Preprocess/postprocess logic |
| test_api.py | 310 | All 26 API endpoints |
| test_agent.py | 294 | Agent pipeline stages |
| test_decision_support.py | 281 | All decision support engines |
| test_query_expansion.py | 255 | All synonym maps and expansion |
| test_integration.py | 247 | End-to-end integration |
| test_rag_engine.py | 194 | RAG retrieval and scoring |
| test_collections.py | 152 | Schema validation |
| test_settings.py | 133 | Configuration validation |
| test_knowledge.py | 123 | Knowledge completeness |

### 15.2 Results

- **Total:** 769 tests
- **Passed:** 769
- **Failed:** 0
- **Pass rate:** 100%
- **Execution time:** 0.47s

---

*Clinical Trial Intelligence Agent v2.0.0 -- HCLS AI Factory -- March 2026*
