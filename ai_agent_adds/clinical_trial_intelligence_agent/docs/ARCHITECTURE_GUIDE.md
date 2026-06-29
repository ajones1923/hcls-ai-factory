# Clinical Trial Intelligence Agent -- Architecture Guide

**Version:** 2.0.0
**Date:** March 22, 2026
**Author:** Adam Jones
**Platform:** NVIDIA DGX Spark -- HCLS AI Factory

---

## Table of Contents

1. [System Overview](#1-system-overview)
2. [Three-Tier Architecture](#2-three-tier-architecture)
3. [Data Flow](#3-data-flow)
4. [Component Architecture](#4-component-architecture)
5. [RAG Engine Design](#5-rag-engine-design)
6. [Workflow Engine Design](#6-workflow-engine-design)
7. [Decision Support Architecture](#7-decision-support-architecture)
8. [Agent Pipeline](#8-agent-pipeline)
9. [Query Expansion System](#9-query-expansion-system)
10. [Collection Architecture](#10-collection-architecture)
11. [Cross-Agent Integration](#11-cross-agent-integration)
12. [Ingest Pipeline Architecture](#12-ingest-pipeline-architecture)
13. [Observability Architecture](#13-observability-architecture)
14. [Security Architecture](#14-security-architecture)
15. [Deployment Architecture](#15-deployment-architecture)

---

## 1. System Overview

### 1.1 Platform Context: HCLS AI Factory 3-Engine Architecture

The Clinical Trial Intelligence Agent operates within the HCLS AI Factory, a three-engine precision medicine platform on NVIDIA DGX Spark:

- **Stage 1 -- Genomics Engine:** Parabricks/DeepVariant/BWA-MEM2 for FASTQ-to-VCF variant calling (GPU-accelerated)
- **Stage 2 -- RAG/Chat Engine:** Milvus (3.56M vectors) + Claude AI for variant interpretation and evidence synthesis
- **Stage 3 -- Drug Discovery Engine:** BioNeMo MolMIM/DiffDock/RDKit for lead optimization across 171 druggable targets

The platform includes **11 intelligence agents**:

| # | Agent | Port | Domain |
|---|---|---|---|
| 1 | Biomarker Intelligence | :8529 | Biomarker discovery and stratification |
| 2 | Oncology Intelligence | :8527/:8528 | Cancer genomics and treatment selection |
| 3 | CAR-T Intelligence | -- | CAR-T cell therapy development |
| 4 | Imaging Intelligence | :8524 | Medical imaging AI and radiology |
| 5 | Autoimmune Intelligence | -- | Autoimmune disease genomics |
| 6 | Pharmacogenomics Intelligence | :8107 | Drug metabolism and dosing |
| 7 | **Clinical Trial Intelligence** | **:8538** | **Trial design and patient matching (this agent)** |
| 8 | Rare Disease Diagnostic | :8134 | Rare disease diagnosis |
| 9 | Single-Cell Intelligence | :8540 | Single-cell transcriptomics |
| 10 | Cardiology Intelligence | :8126 | Cardiac genetics and risk |
| 11 | Neurology Intelligence | -- | Neurological genetics |

### 1.2 Agent Design

The Clinical Trial Intelligence Agent is designed as a modular, layered system that separates concerns across three tiers (presentation, application, data) with well-defined interfaces between each layer. The architecture prioritizes graceful degradation: every component operates without Milvus, without the LLM, and without peer agents, ensuring the system delivers value at every connectivity level.

### Design Principles

1. **Graceful Degradation:** Each tier and component operates independently. Milvus down? Workflows and knowledge base still function. LLM unavailable? Search-only mode. Peer agents offline? Default responses logged.
2. **Workflow-First Design:** All clinical trial intelligence is channeled through 10 typed workflows following a common preprocess/execute/postprocess contract.
3. **Multi-Collection RAG:** Evidence is distributed across 14 specialized collections, with workflow-specific weights ensuring domain relevance.
4. **Type Safety:** Pydantic models enforce data contracts at every boundary (API input/output, workflow I/O, search results).
5. **Observable:** Prometheus metrics at every layer (query, search, workflow, LLM, ingest, system).

---

## 2. Three-Tier Architecture

```
+==============================================================+
|                    PRESENTATION TIER                           |
|  +------------------+                                         |
|  | Streamlit UI     |   5 tabs: Intelligence, Matching,       |
|  | Port 8128        |   Protocol, Competitive, Dashboard      |
|  | NVIDIA Theme     |                                         |
|  +--------+---------+                                         |
|           | HTTP/REST (JSON)                                  |
+===========|==================================================+
            v
+==============================================================+
|                    APPLICATION TIER                            |
|  +------------------+  +----------------+  +---------------+  |
|  | FastAPI API      |  | Workflow       |  | Decision      |  |
|  | Port 8538        |  | Engine (10)    |  | Support (5)   |  |
|  | 26 Endpoints     |  |                |  |               |  |
|  | CORS, Auth, Rate |  | Protocol       |  | Confidence    |  |
|  +--------+---------+  | Patient Match  |  | Complexity    |  |
|           |             | Site Select    |  | Enrollment    |  |
|  +--------+---------+  | Eligibility    |  | Eligibility   |  |
|  | RAG Engine       |  | Adaptive       |  | Competitive   |  |
|  | Multi-collection |  | Safety Signal  |  | Success Rate  |  |
|  | Weighted search  |  | Regulatory     |  +---------------+  |
|  +--------+---------+  | Competitive    |                     |
|           |             | Diversity      |  +---------------+  |
|  +--------+---------+  | DCT Planning   |  | Query         |  |
|  | Agent Pipeline   |  +----------------+  | Expansion     |  |
|  | Plan-Search-     |                      | 10 maps       |  |
|  | Evaluate-Synth   |  +----------------+  | 140 aliases   |  |
|  +--------+---------+  | Knowledge      |  +---------------+  |
|           |             | Base           |                     |
|           |             | 40 trials      |  +---------------+  |
|           |             | 13 areas       |  | Cross-Agent   |  |
|           |             | 9 agencies     |  | Integration   |  |
|           |             +----------------+  | 4 agents      |  |
|           |                                 +---------------+  |
+===========|==================================================+
            v
+==============================================================+
|                       DATA TIER                               |
|  +------------------+  +----------------+  +---------------+  |
|  | Milvus           |  | etcd           |  | MinIO         |  |
|  | Port 19530       |  | Port 2379      |  | Port 9000     |  |
|  | 14 Collections   |  | Metadata       |  | Blob Storage  |  |
|  | IVF_FLAT/COSINE  |  |                |  |               |  |
|  | 384-dim BGE      |  |                |  |               |  |
|  +------------------+  +----------------+  +---------------+  |
+==============================================================+
```

### Tier Responsibilities

| Tier | Responsibility | Failure Mode |
|---|---|---|
| Presentation | User interaction, visualization, form input | API errors shown as warnings |
| Application | Business logic, search, synthesis, decision support | Degrades per component |
| Data | Vector storage, indexing, metadata, blob storage | Workflows still run from knowledge base |

---

## 3. Data Flow

### 3.1 Query Flow

```
User Question
    |
    v
[Streamlit UI] --> HTTP POST /v1/trial/query
    |
    v
[FastAPI API]
    |
    v
[Query Expansion] --> Resolve aliases, expand synonyms
    |
    v
[Agent Pipeline]
    |
    +-> [Plan] --> Detect workflows, decompose sub-questions
    |
    +-> [Search] --> RAG Engine
    |       |
    |       +-> [Embed query] --> BGE-small-en-v1.5 (384-dim)
    |       |
    |       +-> [Search 14 collections] --> Milvus IVF_FLAT
    |       |
    |       +-> [Weight & merge] --> Workflow-specific weights
    |       |
    |       +-> [Score threshold] --> Filter < 0.4
    |
    +-> [Evaluate] --> Evidence quality, completeness check
    |
    +-> [Workflow Execute] --> One or more of 10 workflows
    |       |
    |       +-> [Decision Support] --> Confidence, complexity, etc.
    |       |
    |       +-> [Cross-Agent] --> Optional oncology/PGx/cardio/biomarker
    |
    +-> [Synthesize] --> LLM (Claude) generates natural language answer
    |
    +-> [Report] --> Format with citations, guidelines, confidence
    |
    v
[TrialResponse] --> JSON to Streamlit UI
```

### 3.2 Ingest Flow

```
[External Sources]
    |
    +-> ClinicalTrials.gov API --> ClinicalTrialsParser
    |                                    |
    |                                    +-> trial_protocols
    |                                    +-> trial_eligibility
    |                                    +-> trial_endpoints
    |                                    +-> trial_sites
    |                                    +-> trial_investigators
    |
    +-> PubMed E-utilities API --> PubMedParser
    |                                    |
    |                                    +-> trial_literature
    |                                    +-> trial_results
    |
    +-> FDA/EMA/ICH docs --> RegulatoryParser
                                    |
                                    +-> trial_regulatory
                                    +-> trial_guidelines
                                    +-> trial_safety
```

### 3.3 Cross-Agent Flow

```
[Clinical Trial Agent]
    |
    +-- query_oncology_agent() --> [:8527] Molecular matches
    |       |
    |       +-- Graceful degradation if unavailable
    |
    +-- query_pgx_agent() ------> [:8107] PGx screening
    |
    +-- query_cardiology_agent() -> [:8126] Cardiac safety
    |
    +-- query_biomarker_agent() --> [:8529] Biomarker enrichment
    |
    v
[integrate_cross_agent_results()] --> Unified assessment
```

---

## 4. Component Architecture

### 4.1 Component Dependency Graph

```
api/main.py (FastAPI)
    |
    +-> api/routes/trial_clinical.py (22 endpoints)
    |       |
    |       +-> src/agent.py (TrialIntelligenceAgent)
    |       |       |
    |       |       +-> src/rag_engine.py (TrialRAGEngine)
    |       |       |       |
    |       |       |       +-> src/collections.py (14 schemas)
    |       |       |       +-> src/query_expansion.py (10 maps)
    |       |       |
    |       |       +-> src/clinical_workflows.py (10 workflows)
    |       |       |       |
    |       |       |       +-> src/decision_support.py (5 engines)
    |       |       |       +-> src/knowledge.py (domain knowledge)
    |       |       |       +-> src/models.py (data contracts)
    |       |       |
    |       |       +-> src/cross_modal.py (4 agent integrations)
    |       |
    |       +-> src/export.py (report generation)
    |
    +-> api/routes/reports.py (2 endpoints)
    +-> api/routes/events.py (2 endpoints)
    +-> src/metrics.py (Prometheus)
    +-> config/settings.py (TrialSettings)
```

### 4.2 Module Roles

| Module | Role | Dependencies |
|---|---|---|
| `agent.py` | Orchestrator: plan, search, evaluate, synthesize | rag_engine, workflows, cross_modal |
| `rag_engine.py` | Multi-collection RAG with weighted retrieval | collections, query_expansion |
| `clinical_workflows.py` | 10 workflow implementations | models, knowledge, decision_support |
| `decision_support.py` | Quantitative scoring engines | models |
| `knowledge.py` | Static domain knowledge (40 trials, 13 areas, etc.) | -- (no dependencies) |
| `models.py` | Pydantic data contracts and enums | pydantic |
| `collections.py` | Milvus schema definitions | models, pymilvus |
| `query_expansion.py` | Synonym resolution and term expansion | models |
| `cross_modal.py` | Cross-agent HTTP integration | config/settings |
| `export.py` | Report formatting and export | models |
| `metrics.py` | Prometheus metric definitions | prometheus_client |
| `scheduler.py` | Timed ingest and maintenance | ingest, settings |

---

## 5. RAG Engine Design

### 5.1 Multi-Collection Search Strategy

The RAG engine searches all 14 collections in parallel, with each collection receiving a workflow-specific weight. Results are merged by weighted score and deduplicated.

```
Query
  |
  v
[Expand Query] --> query_expansion.py
  |
  v
[Embed] --> BGE-small-en-v1.5 (384-dim vector)
  |
  v
[Parallel Search] --> 14 collections x top_k results
  |
  +-> trial_protocols     (weight: varies by workflow)
  +-> trial_eligibility   (weight: varies by workflow)
  +-> ...
  +-> genomic_evidence    (weight: varies by workflow)
  |
  v
[Weight & Merge] --> score = raw_score * collection_weight
  |
  v
[Threshold Filter] --> remove if score < 0.4
  |
  v
[Re-rank] --> sort by weighted score descending
  |
  v
[Return top_k] --> TrialSearchResult list
```

### 5.2 Graceful Degradation Levels

| Level | Milvus | LLM | Peer Agents | Capability |
|---|---|---|---|---|
| Full | Available | Available | Available | Complete RAG + synthesis + cross-agent |
| Search-only | Available | Unavailable | Any | Vector search with structured results |
| Workflow-only | Unavailable | Available | Any | Knowledge-based workflows + LLM synthesis |
| Minimal | Unavailable | Unavailable | Unavailable | Decision engines + knowledge base queries |

---

## 6. Workflow Engine Design

### 6.1 Template Method Pattern

All workflows follow the `BaseTrialWorkflow` contract, which enforces a three-phase execution model:

```python
class BaseTrialWorkflow(ABC):
    workflow_type: TrialWorkflowType

    def run(self, inputs: dict) -> WorkflowResult:
        processed = self.preprocess(inputs)    # validate, normalize
        result = self.execute(processed)        # core logic
        result = self.postprocess(result)       # enrich, add warnings
        return result
```

### 6.2 Workflow Routing

The `WorkflowEngine` maps `TrialWorkflowType` enum values to workflow class instances, enabling both automatic detection from query text and explicit routing via the API:

```
Query Text --> [Detect Keywords] --> TrialWorkflowType
    |
    v
WorkflowEngine.dispatch(workflow_type, inputs)
    |
    v
[Specific Workflow Instance].run(inputs)
    |
    v
WorkflowResult
```

### 6.3 Collection Weight Boosting

Each workflow overrides the default collection weights to prioritize domain-relevant evidence. For example, the Safety Signal workflow boosts `trial_safety` to 0.25 (vs. default 0.08), ensuring adverse event data surfaces prominently.

---

## 7. Decision Support Architecture

### 7.1 Engine Independence

All five decision support engines are stateless, pure-function classes with no external dependencies (no Milvus, no LLM, no network calls). This ensures they operate at every degradation level.

```
Decision Engines (all stateless, self-contained)
    |
    +-> ConfidenceCalibrator     --> multi-factor calibration
    +-> ProtocolComplexityScorer --> Tufts CSDD benchmarks
    +-> EnrollmentPredictor      --> multi-factor prediction
    +-> EligibilityAnalyzer      --> 29 pattern matching
    +-> CompetitiveThreatScorer  --> 4-factor threat model
    +-> HistoricalSuccessEstimator -> 12 areas x 3 phases
```

### 7.2 Calibration Pipeline

The Confidence Calibrator is applied as a post-processing step after workflow execution and RAG retrieval, combining four signals into a single calibrated confidence:

```
raw_confidence (0.30)  --> from workflow
evidence_base (0.30)   --> from evidence level (A1=1.0, E=0.15)
doc_factor (0.20)      --> log(n_docs + 1) / log(12)
agreement (0.20)       --> cross-agent consensus
                           = calibrated confidence (0.0-1.0)
```

---

## 8. Agent Pipeline

### 8.1 Five-Stage Pipeline

The `TrialIntelligenceAgent` implements the VAST AI OS AgentEngine pattern:

| Stage | Method | Input | Output |
|---|---|---|---|
| Plan | `search_plan()` | Query text | SearchPlan (areas, drugs, biomarkers, sub-questions) |
| Search | `rag_engine.query()` | SearchPlan | TrialSearchResult list |
| Evaluate | `evaluate_evidence()` | Search results | Quality scores, completeness flags |
| Synthesize | `synthesize()` | Results + evaluation | Natural language answer (via Claude) |
| Report | `generate_report()` | All above | TrialResponse with citations |

### 8.2 Evidence Hierarchy

```
Level 1a: Systematic review of RCTs (highest)
Level 1b: Individual RCT
Level 2a: Systematic review of cohort studies
Level 2b: Individual cohort study
Level 3:  Case-control study
Level 4:  Case series
Level 5:  Expert opinion
Reg:      Regulatory guidance (FDA/EMA/ICH)
```

---

## 9. Query Expansion System

### 9.1 Expansion Pipeline

```
Raw Query: "What NSCLC trials use Keytruda with TMB-H?"
    |
    v
[Entity Alias Resolution]
    NSCLC -> "non-small cell lung cancer"
    Keytruda -> "pembrolizumab"
    TMB-H -> "tumor mutational burden high"
    |
    v
[Therapeutic Area Detection]
    "lung cancer" -> oncology
    |
    v
[Drug Expansion]
    "pembrolizumab" -> ["Keytruda", "MK-3475", "lambrolizumab", "anti-PD-1"]
    |
    v
[Biomarker Expansion]
    "TMB" -> ["tumor mutational burden", "mutational load", "TMB-high", ...]
    |
    v
[Expanded Query]
    "non-small cell lung cancer pembrolizumab anti-PD-1
     tumor mutational burden TMB-high oncology"
```

### 9.2 Map Statistics

- 10 synonym maps covering all clinical trial domains
- 140 entity aliases for instant abbreviation resolution
- 33 drug entries with brand/generic/code name mapping
- 22 biomarker entries with assay and synonym coverage

---

## 10. Collection Architecture

### 10.1 Schema Design Principles

1. **Shared embedding field:** All 14 collections use identical 384-dim FLOAT_VECTOR fields
2. **Shared index config:** IVF_FLAT with COSINE metric and nlist=128
3. **Domain-specific metadata:** Each collection has typed metadata fields (VARCHAR, INT32, FLOAT, BOOL) relevant to its domain
4. **Auto-generated primary keys:** INT64 auto_id for all collections
5. **Text content fields:** VARCHAR with appropriate max_length for full-text and chunk storage

### 10.2 Estimated Data Distribution

```
trial_eligibility  ████████████████████████████████████████████████████  50,000
trial_sites        ████████████████████████████████                      30,000
trial_endpoints    ████████████████████                                  20,000
trial_safety       ████████████████████                                  20,000
trial_literature   ██████████                                           10,000
trial_protocols    █████                                                 5,000
trial_investigators █████                                                5,000
trial_biomarkers   ███                                                   3,000
trial_results      ███                                                   3,000
trial_regulatory   ██                                                    2,000
trial_rwe          ██                                                    2,000
trial_guidelines   █                                                     1,000
trial_adaptive     ▌                                                       500
genomic_evidence   ████████████████████████████████████████████████ ~100,000
```

---

## 11. Cross-Agent Integration

### 11.1 Integration Architecture (8 Agents)

```
                         +----------------------------+
                         | Clinical Trial Intelligence |
                         |         Agent (:8538)       |
                         +--+--+--+--+--+--+--+--+---+
                            |  |  |  |  |  |  |  |
          +-----------------+  |  |  |  |  |  |  +----------------+
          v                    v  |  |  |  |  v                   v
  +-------+-------+  +--------++ |  |  |  | ++--------+  +-------+------+
  | Oncology      |  | PGx     | |  |  |  | | Rare    |  | Imaging      |
  | Intelligence  |  | Intelli-| |  |  |  | | Disease |  | Intelligence |
  | (:8527)       |  | gence   | |  |  |  | | (:8134) |  | (:8524)      |
  | Molecular     |  | (:8107) | |  |  |  | | Orphan  |  | Imaging      |
  | trial matches |  | PGx     | |  |  |  | | drug    |  | endpoints    |
  +---------------+  | screen  | |  |  |  | | trials  |  | for trial    |
                     +--------++ |  |  |  | +--------++  | eligibility  |
                                 v  |  |  v               +-------------+
                     +-----------++ |  | ++----------+
                     | Cardiology | |  | | Neurology |
                     | (:8126)    | |  | | Agent     |
                     | Cardiac    | |  | | Neuro     |
                     | safety     | |  | | trial     |
                     +------------+ |  | | matching  |
                                    v  | +-----------+
                        +-----------++ |
                        | Biomarker  | v
                        | (:8529)    | +-------------+
                        | Enrichment | | Single-Cell |
                        | strategies | | (:8540)     |
                        +------------+ | Cellular    |
                                       | biomarker   |
                                       | endpoints   |
                                       +-------------+
```

### 11.2 Cross-Agent Trigger Map

| Peer Agent | Port | Trigger | Data Exchanged |
|---|---|---|---|
| Oncology | :8527 | Oncology trial query | Molecular matches, biomarker-driven trial eligibility |
| PGx | :8107 | PGx-stratified trial | Metabolizer phenotype for enrollment criteria |
| Cardiology | :8126 | Cardiac safety endpoint | QTc risk, cardiotoxicity monitoring requirements |
| Biomarker | :8529 | Biomarker-enriched trial | Enrichment strategies, companion diagnostic requirements |
| Rare Disease | :8134 | Orphan drug trial | Gene therapy trial eligibility, rare disease patient matching |
| Neurology | -- | Neuro trial design | CNS endpoint design, blood-brain barrier considerations |
| Single-Cell | :8540 | Cellular biomarker trial | Single-cell biomarker endpoints, MRD monitoring criteria |
| Imaging | :8524 | Imaging endpoint trial | Imaging response criteria (RECIST, RANO), imaging-based eligibility |

### 11.3 Pediatric Oncology Trial Intelligence

The agent provides specialized support for pediatric oncology clinical trials:

- **COG Cooperative Group Trials:** Knowledge base includes landmark COG trials -- AALL0932 (standard-risk B-ALL), AALL1231 (high-risk T-ALL), ANBL1232 (high-risk neuroblastoma). Trial matching considers COG risk stratification and treatment arms.
- **Pediatric MATCH (APEC1621):** Precision medicine basket trial matching pediatric solid tumors to molecular-targeted therapies based on actionable genomic alterations.
- **Age-Specific Filters:** Pediatric trial matching applies age-appropriate filters (neonatal 0-28d, infant 1-12mo, child 1-12yr, adolescent 12-18yr) and accounts for age-specific eligibility windows.
- **Rolling-6 Design:** Phase I dose-escalation design used in COG trials where up to 6 patients can be enrolled simultaneously; the adaptive design evaluation workflow models rolling enrollment patterns.
- **EFS Endpoints:** Event-free survival as the standard primary endpoint for pediatric oncology trials, with appropriate event definitions (relapse, progression, second malignancy, death from any cause).
- **Consent/Assent:** Trial matching includes consent/assent requirements: parental consent required for all minors, written assent for children >= 7 years (age varies by IRB), and adolescent assent for 12-17 year olds.

### 11.4 Failure Isolation

Each cross-agent call is wrapped in a try/except with a 30-second timeout. Failure in one agent never blocks the clinical trial agent's response. The integration module (`cross_modal.py`) returns structured default responses when agents are unavailable, clearly flagged as degraded.

---

## 12. Ingest Pipeline Architecture

### 12.1 Pipeline Hierarchy

```
BaseIngestPipeline (abstract)
    |
    +-> ClinicalTrialsParser (ClinicalTrials.gov)
    |       |
    |       +-> XML/JSON parsing
    |       +-> Field extraction per collection
    |       +-> BGE embedding generation
    |       +-> Milvus upsert
    |
    +-> PubMedParser (PubMed/MEDLINE)
    |       |
    |       +-> E-utilities API calls
    |       +-> MeSH term extraction
    |       +-> Literature chunking
    |
    +-> RegulatoryParser (FDA/EMA/ICH)
            |
            +-> Document parsing
            +-> Guideline extraction
            +-> Safety signal parsing
```

### 12.2 Scheduler Integration

The `src/scheduler.py` runs ingest pipelines on a configurable interval (default 24 hours) via a daemon thread. It handles collection maintenance (compaction, index rebuild) and supports enable/disable via environment variable.

---

## 13. Observability Architecture

```
[Application Code]
    |
    +-> src/metrics.py (Prometheus client)
    |       |
    |       +-> Counters: queries_total, search_total, errors_total
    |       +-> Histograms: query_duration, search_duration
    |       +-> Gauges: collection_records
    |       +-> Info: system_info
    |
    +-> GET /metrics endpoint
            |
            v
    [Prometheus Server] --> [Grafana Dashboard]
```

All metrics use the `trial_` prefix for dashboard filtering. The metrics module gracefully degrades if `prometheus_client` is not installed (no-op stubs).

---

## 14. Security Architecture

### 14.1 Security Layers

```
[Client Request]
    |
    v
[CORS Middleware] --> Check origin against whitelist
    |
    v
[Rate Limiter] --> 100 req/min per IP
    |
    v
[Auth Middleware] --> Validate X-API-Key header (if configured)
    |
    v
[Pydantic Validation] --> Enforce field constraints, types, ranges
    |
    v
[Business Logic] --> No SQL injection risk (vector-only DB)
    |
    v
[Response] --> Sanitized JSON output
```

### 14.2 Secret Management

- API keys stored in environment variables or `.env` file
- `TRIAL_` prefix isolates agent-specific config
- No secrets in source code or version control
- HTTPS termination at reverse proxy layer

---

## 15. Deployment Architecture

### 15.1 Docker Composition

```
+-------------------------------------------------------+
|                  DGX Spark Host                        |
|                                                        |
|  +------------------+  +-------------------+           |
|  | clinical-trial-  |  | clinical-trial-   |           |
|  | agent-api        |  | agent-ui          |           |
|  | Port: 8538       |  | Port: 8128        |           |
|  | FastAPI + Uvicorn |  | Streamlit         |           |
|  +--------+---------+  +---------+---------+           |
|           |                      |                     |
|           +----------+-----------+                     |
|                      |                                 |
|              +-------+-------+                         |
|              | milvus-       |                         |
|              | standalone    |                         |
|              | Port: 19530   |                         |
|              +---+---+---+---+                         |
|                  |       |                             |
|          +-------+  +----+------+                      |
|          | etcd  |  | MinIO     |                      |
|          | :2379 |  | :9000     |                      |
|          +-------+  +-----------+                      |
|                                                        |
|  +------------------+  +-------------------+           |
|  | Other Agents     |  | Monitoring        |           |
|  | Oncology (:8527) |  | Prometheus        |           |
|  | PGx (:8107)      |  | Grafana           |           |
|  | Cardio (:8126)   |  |                   |           |
|  | Biomarker (:8529)|  |                   |           |
|  +------------------+  +-------------------+           |
+-------------------------------------------------------+
```

### 15.2 Resource Requirements

| Component | CPU | Memory | GPU | Storage |
|---|---|---|---|---|
| API Server | 2 cores | 2 GB | None | <100 MB |
| Streamlit UI | 1 core | 512 MB | None | <50 MB |
| Milvus Standalone | 4 cores | 8 GB | None | 10-50 GB |
| BGE Embedding | 2 cores | 1 GB | Optional | ~500 MB model |

---

*Clinical Trial Intelligence Agent v2.0.0 -- Architecture Guide -- March 2026*
