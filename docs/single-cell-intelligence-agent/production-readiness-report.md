# Single-Cell Intelligence Agent -- Production Readiness Report

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones
**Status:** Production Ready (Conditional)
**Classification:** Internal Engineering Document

---

## Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [System Overview](#2-system-overview)
3. [Architecture Assessment](#3-architecture-assessment)
4. [Milvus Collection Infrastructure](#4-milvus-collection-infrastructure)
5. [Knowledge Base Inventory](#5-knowledge-base-inventory)
6. [Workflow Engine Assessment](#6-workflow-engine-assessment)
7. [Decision Support Engine Assessment](#7-decision-support-engine-assessment)
8. [Data Model Completeness](#8-data-model-completeness)
9. [API Surface Assessment](#9-api-surface-assessment)
10. [Authentication and Security](#10-authentication-and-security)
11. [Embedding Pipeline](#11-embedding-pipeline)
12. [LLM Integration](#12-llm-integration)
13. [GPU Acceleration Readiness](#13-gpu-acceleration-readiness)
14. [Docker and Container Infrastructure](#14-docker-and-container-infrastructure)
15. [Seed Data and Ingest Pipeline](#15-seed-data-and-ingest-pipeline)
16. [Cross-Agent Integration](#16-cross-agent-integration)
17. [Test Coverage Analysis](#17-test-coverage-analysis)
18. [Performance Benchmarks](#18-performance-benchmarks)
19. [Monitoring and Observability](#19-monitoring-and-observability)
20. [Configuration Management](#20-configuration-management)
21. [Error Handling and Resilience](#21-error-handling-and-resilience)
22. [Streamlit UI Assessment](#22-streamlit-ui-assessment)
23. [Known Issues and Technical Debt](#23-known-issues-and-technical-debt)
24. [Risk Register](#24-risk-register)
25. [Go/No-Go Recommendation](#25-gono-go-recommendation)

---

## 1. Executive Summary

The Single-Cell Intelligence Agent is the twelfth intelligence agent in the HCLS AI Factory platform. It provides RAG-powered clinical decision support for single-cell transcriptomics, spatial biology, tumor microenvironment profiling, drug response prediction, and CAR-T target validation.

**Key metrics at a glance:**

| Metric | Value |
|--------|-------|
| Milvus collections | 12 (11 domain-specific + 1 shared genomic) |
| Analysis workflows | 10 (+ 1 general) |
| Decision support engines | 4 |
| Cell types in knowledge base | 44 |
| Drugs modeled | 30 |
| Marker genes | 75 |
| Immune signatures | 10 |
| Ligand-receptor pairs | 25 |
| Cancer TME atlas profiles | 12 |
| Clinical conditions | 36 |
| Agent cell type aliases | 232 |
| CellxGene seed records | 49 |
| Marker seed records | 75 |
| TME seed records | 20 |
| Biomarkers | 23 |
| Spatial interaction maps | 14 |
| Test files | 12 |
| Total test lines | 1,760 (185 test cases estimated) |
| Source code lines | 14,560 |
| API port | 8540 |
| Streamlit port | 8130 |

**Verdict:** PRODUCTION READY (Conditional) -- the agent is architecturally sound and feature-complete for the documented use cases. Conditions for full production clearance are listed in Section 25.

---

## 2. System Overview

### 2.1 Purpose

The Single-Cell Intelligence Agent resolves single-cell transcriptomics data into clinically actionable insights. It addresses the "resolution gap" in precision medicine by operating at individual-cell granularity rather than bulk-tissue averages.

### 2.2 Capabilities

- **Cell type annotation** with Cell Ontology mapping, marker-based scoring, and LLM-augmented consensus
- **Tumor microenvironment classification** into four immunophenotypes (hot-inflamed, cold-desert, excluded, immunosuppressive) with treatment recommendations
- **Drug response prediction** at cellular resolution using GDSC/DepMap signatures
- **Subclonal architecture detection** with escape risk scoring and timeline estimation
- **Spatial transcriptomics analysis** across Visium, MERFISH, Xenium, and CODEX platforms
- **Trajectory inference** for differentiation, activation, exhaustion, EMT, and stemness
- **Ligand-receptor interaction mapping** with CellPhoneDB/NicheNet-style analysis
- **Biomarker discovery** with cell-type specificity scoring and clinical correlation
- **CAR-T target validation** with on-tumor/off-tumor safety profiling
- **Treatment monitoring** through longitudinal clonal dynamics tracking

### 2.3 Integration Points

The agent integrates with four peer agents via REST API:
- Genomics Agent (port 8527) -- variant-level evidence
- Biomarker Agent (port 8529) -- cross-modal biomarker correlation
- Oncology Agent (port 8528) -- therapy line recommendations
- Clinical Trial Agent (port 8538) -- trial matching for novel targets

---

## 3. Architecture Assessment

### 3.1 Component Architecture

```
Streamlit UI (8130)
     |
FastAPI REST API (8540)
     |
+----+----+--------+--------+
|         |        |        |
RAG     Workflow  Decision  Knowledge
Engine   Engine   Support   Base
     |         |        |
     +----+----+
          |
     Milvus Vector DB (19530)
          |
     +----+----+
     |         |
   etcd      MinIO
```

### 3.2 Module Inventory

| Module | File | Lines | Purpose |
|--------|------|-------|---------|
| Agent | `src/agent.py` | 2,090 | Autonomous reasoning, system prompt, enums, workflow dispatch |
| Models | `src/models.py` | 820 | Pydantic data models (15 model classes, 12 enums) |
| Collections | `src/collections.py` | 1,210 | 12 Milvus collection schemas with field definitions |
| RAG Engine | `src/rag_engine.py` | 1,490 | Multi-collection search, conversation memory, LLM synthesis |
| Clinical Workflows | `src/clinical_workflows.py` | 1,792 | 10 analysis workflows with BaseSCWorkflow pattern |
| Decision Support | `src/decision_support.py` | 886 | 4 clinical engines (TME, subclonal, target, deconvolution) |
| Knowledge | `src/knowledge.py` | 1,816 | Domain knowledge base (cell types, drugs, markers, TME profiles) |
| Query Expansion | `src/query_expansion.py` | 893 | Synonym expansion for cell types, genes, diseases |
| Cross-Modal | `src/cross_modal.py` | 392 | Inter-agent communication layer |
| Metrics | `src/metrics.py` | 476 | Prometheus metrics collection |
| Export | `src/export.py` | 588 | Report generation (PDF/DOCX/JSON) |
| Scheduler | `src/scheduler.py` | 496 | APScheduler-based ingest scheduling |
| Settings | `config/settings.py` | 197 | Pydantic BaseSettings with validation |
| API Main | `api/main.py` | 615 | FastAPI application factory with middleware |
| Routes | `api/routes/` | ~800 | Versioned clinical, report, and event routes |
| UI | `app/sc_ui.py` | ~600 | 5-tab Streamlit interface |
| Ingest | `src/ingest/` | ~1,611 | CellxGene, marker, and TME parsers |

### 3.3 Architectural Strengths

- Clean separation between RAG engine, workflow engine, and decision support
- Consistent Pydantic model validation across all data boundaries
- Graceful degradation when dependencies (Milvus, LLM, embedder) are unavailable
- Thread-safe metrics collection with lock-protected counters
- Conversation persistence with 24-hour TTL and disk-backed storage

### 3.4 Architectural Concerns

- Dual `SCWorkflowType` enum definitions (one in `models.py`, one in `agent.py`) -- values slightly diverge
- RAG engine imports from `agent.py` rather than `models.py` for some types
- No circuit breaker pattern for cross-agent REST calls

---

## 4. Milvus Collection Infrastructure

### 4.1 Collection Inventory

| # | Collection | Fields | Est. Records | Search Weight | Purpose |
|---|-----------|--------|-------------|--------------|---------|
| 1 | `sc_cell_types` | 9 | 5,000 | 0.14 | Cell type annotations with CL ontology IDs |
| 2 | `sc_markers` | 9 | 50,000 | 0.12 | Gene markers with specificity scores |
| 3 | `sc_spatial` | 9 | 10,000 | 0.10 | Spatial transcriptomics niches |
| 4 | `sc_tme` | 9 | 8,000 | 0.10 | Tumor microenvironment profiles |
| 5 | `sc_drug_response` | 9 | 25,000 | 0.09 | Drug sensitivity predictions |
| 6 | `sc_literature` | 9 | 50,000 | 0.08 | Published scRNA-seq literature |
| 7 | `sc_methods` | 9 | 2,000 | 0.07 | Analytical tools and pipelines |
| 8 | `sc_datasets` | 10 | 15,000 | 0.06 | Reference atlases (CellxGene, HCA) |
| 9 | `sc_trajectories` | 9 | 8,000 | 0.07 | Pseudotime and trajectory data |
| 10 | `sc_pathways` | 7 | 20,000 | 0.07 | Signaling/metabolic pathways |
| 11 | `sc_clinical` | 9 | 12,000 | 0.07 | Clinical biomarker correlations |
| 12 | `genomic_evidence` | 6 | 3,560,000 | 0.03 | Shared genomic variants (read-only) |
| | **Total** | | **3,765,000** | **1.00** | |

### 4.2 Index Configuration

All collections use identical index parameters:

| Parameter | Value |
|-----------|-------|
| Embedding dimension | 384 (BGE-small-en-v1.5) |
| Index type | IVF_FLAT |
| Metric type | COSINE |
| nlist | 128 |
| Primary key | INT64, auto_id |

### 4.3 Workflow-Specific Weight Profiles

The system defines 11 weight profiles, one per workflow type. Each profile redistributes the 1.0 total weight across all 12 collections to prioritize domain-relevant collections. Example:

| Workflow | Primary Collection | Primary Weight | Secondary | Secondary Weight |
|----------|-------------------|----------------|-----------|-----------------|
| Cell Type Annotation | sc_cell_types | 0.25 | sc_markers | 0.22 |
| TME Profiling | sc_tme | 0.25 | sc_cell_types | 0.15 |
| Drug Response | sc_drug_response | 0.25 | sc_tme | 0.12 |
| Spatial Niche | sc_spatial | 0.28 | sc_cell_types | 0.12 |
| CAR-T Validation | sc_markers | 0.18 | sc_cell_types | 0.15 |
| Biomarker Discovery | sc_markers | 0.20 | sc_clinical | 0.18 |

### 4.4 Collection Schema Quality

- All collections follow consistent naming conventions (`sc_` prefix)
- All include auto-incrementing INT64 primary keys
- All embed the standard 384-dim FLOAT_VECTOR field
- VARCHAR max_length is appropriately sized per field (32-8192)
- All descriptions are populated and meaningful

**Assessment:** PASS -- collection infrastructure is well-designed and production-ready.

---

## 5. Knowledge Base Inventory

### 5.1 Cell Type Atlas

The knowledge base contains **44 cell types** organized across immune, stromal, epithelial, neural, and specialized compartments:

| Compartment | Cell Types | Count |
|-------------|-----------|-------|
| T lymphocytes | T_cell, CD8_T, CD4_T, Treg, gamma_delta_T | 5 |
| B lymphocytes | B_cell, Plasma | 2 |
| Innate lymphoid | NK | 1 |
| Myeloid | Monocyte, Macrophage, DC, pDC, Neutrophil, Mast_cell, Basophil, Eosinophil | 8 |
| Stromal | Fibroblast, Endothelial, Pericyte, Smooth_muscle, Adipocyte | 5 |
| Epithelial | Epithelial, Hepatocyte, Podocyte | 3 |
| Neural | Neuron, Astrocyte, Oligodendrocyte, Microglia | 4 |
| Cardiac | Cardiomyocyte | 1 |
| Stem/Progenitor | HSC, Erythroid_progenitor, Megakaryocyte | 3 |
| Cancer-specific | Cancer_stem_cell, Cycling_tumor, Senescent_tumor | 3 |
| Specialized | Melanocyte, Schwann_cell, Mesothelial, Satellite_cell | 4 |
| Other | MDSC, ILC, Platelet, RBC, Doublet | 5 |
| **Total** | | **44** |

Each cell type entry includes:
- Canonical marker genes (5 per type)
- Cell Ontology (CL) identifier
- Tissue distribution
- Subtypes list
- Description

### 5.2 Drug Database

**30 drugs** across 10 drug classes:

| Drug Class | Drugs | Count |
|-----------|-------|-------|
| Checkpoint inhibitors | Pembrolizumab, Nivolumab, Atezolizumab, Durvalumab, Ipilimumab | 5 |
| Tyrosine kinase inhibitors | Osimertinib, Sunitinib, Imatinib, Erlotinib | 4 |
| Targeted therapy | Olaparib, Trastuzumab, Venetoclax, Ibrutinib | 4 |
| Chemotherapy | Cisplatin, Doxorubicin, Paclitaxel, Gemcitabine, Temozolomide | 5 |
| Immunomodulators | Lenalidomide, Thalidomide | 2 |
| Cell therapy | Tisagenlecleucel, Axicabtagene ciloleucel, Brexucabtagene | 3 |
| Bispecifics | Blinatumomab, Teclistamab | 2 |
| ADCs | Trastuzumab deruxtecan, Enfortumab vedotin | 2 |
| Epigenetic | Azacitidine, Decitabine | 2 |
| Radiopharmaceutical | Lutetium-177 PSMA | 1 |

### 5.3 Marker Genes

**75 marker genes** spanning all 44 cell types. Each gene includes:
- Gene symbol and Ensembl ID (where available)
- Associated cell type
- Specificity score (0-1)
- Surface protein flag
- Evidence text

### 5.4 Immune Signatures

**10 curated immune signatures:**

1. Cytotoxic T cell activation
2. T cell exhaustion
3. Treg suppressive program
4. M1 macrophage polarization
5. M2 macrophage polarization
6. Interferon-gamma response
7. Type I interferon response
8. Antigen presentation
9. NK cell cytotoxicity
10. B cell activation/germinal center

### 5.5 Ligand-Receptor Pairs

**25 curated ligand-receptor interaction pairs** covering:
- Checkpoint interactions (PD-L1/PD-1, CTLA-4/CD80)
- Growth factor signaling (EGF/EGFR, HGF/MET)
- Chemokine axes (CXCL12/CXCR4, CCL2/CCR2)
- Notch signaling (DLL1/NOTCH1)
- Wnt signaling (WNT5A/FZD5)
- Hedgehog signaling (SHH/PTCH1)

### 5.6 Cancer TME Atlas

**12 cancer type-specific TME profiles:**

| Cancer Type | Typical TME Class | Key Features |
|------------|------------------|--------------|
| Melanoma | Hot inflamed | High CD8+ T, high PD-L1, checkpoint responsive |
| NSCLC | Variable | Smoking-associated neoantigen load, spatial heterogeneity |
| Breast (TNBC) | Hot/excluded | Stromal barrier, TIL-dependent prognosis |
| Colorectal (MSI-H) | Hot inflamed | High TMB, Lynch syndrome association |
| Colorectal (MSS) | Cold desert | Low immune infiltrate, Wnt-driven |
| Pancreatic | Immunosuppressive | Dense stroma, M2 macrophage dominant |
| Glioblastoma | Cold/suppressive | BBB exclusion, microglia-dominated |
| Renal cell | Excluded | High angiogenesis, anti-VEGF responsive |
| Ovarian | Excluded | Ascites microenvironment, CD8 prognostic |
| Hepatocellular | Variable | Viral vs. non-viral etiology affects TME |
| Head and neck | Hot inflamed | HPV+/HPV- dichotomy |
| Bladder | Variable | TMB-dependent, checkpoint responsive |

### 5.7 Agent Knowledge Statistics

| Category | Count |
|----------|-------|
| Cell types (detailed) | 44 |
| Cell type aliases | 232 |
| Biomarker genes | 23 |
| Spatial interaction maps | 14 |
| Conditions/diseases | 36 |
| Agent-defined cell types | 31 |
| CellxGene seed records | 49 |
| Marker seed records | 75 |
| TME seed records | 20 |

**Assessment:** PASS -- knowledge base is comprehensive and covers the major cell types, drugs, and cancer types encountered in clinical single-cell studies.

---

## 6. Workflow Engine Assessment

### 6.1 Workflow Inventory

| # | Workflow | Input | Output Model | Clinical Value |
|---|---------|-------|-------------|---------------|
| 1 | Cell Type Annotation | Cluster markers, reference | `CellTypeAnnotation` | Cell identity at CL ontology level |
| 2 | TME Profiling | Cell proportions, checkpoints | `TMEProfile` | Immunotherapy response prediction |
| 3 | Drug Response | Cell signatures, drug name | `DrugResponsePrediction` | Sensitivity/resistance at cell level |
| 4 | Subclonal Architecture | Clone frequencies, CNV | `SubclonalResult` | Escape risk and timeline |
| 5 | Spatial Niche | Spatial coordinates, genes | `SpatialNiche` | Tissue architecture mapping |
| 6 | Trajectory Analysis | Pseudotime, driver genes | `TrajectoryResult` | Differentiation path identification |
| 7 | Ligand-Receptor | L-R pairs, cell types | `LigandReceptorInteraction` | Cell communication networks |
| 8 | Biomarker Discovery | DE results, clinical data | `BiomarkerCandidate` | Diagnostic/prognostic markers |
| 9 | CAR-T Validation | Target expression, safety | `CARTTargetValidation` | On-tumor/off-tumor profiling |
| 10 | Treatment Monitoring | Longitudinal samples | `TreatmentResponse` | Resistance tracking |

### 6.2 Workflow Architecture

All workflows implement the `BaseSCWorkflow` abstract base class:

```python
class BaseSCWorkflow(ABC):
    def preprocess(self, query, context) -> dict
    def execute(self, preprocessed) -> WorkflowResult
    def postprocess(self, result) -> WorkflowResult
```

The `WorkflowEngine` class provides unified dispatch via `execute(workflow_type, data)`.

### 6.3 Workflow Quality Assessment

| Criterion | Status | Notes |
|-----------|--------|-------|
| All 10 workflows implemented | PASS | Each has preprocess/execute/postprocess |
| Structured output models | PASS | All return typed Pydantic models |
| Error handling | PASS | Graceful fallback to LLM when engine unavailable |
| Weight profile per workflow | PASS | 11 weight profiles in collections.py |
| Clinical recommendations | PASS | Treatment recs generated for TME, drug, CAR-T |

**Assessment:** PASS.

---

## 7. Decision Support Engine Assessment

### 7.1 Engine Inventory

| # | Engine | Input | Output | Clinical Application |
|---|--------|-------|--------|---------------------|
| 1 | TMEClassifier | Cell proportions, gene expression, PD-L1 TPS | TME class + treatment recs | Immunotherapy selection |
| 2 | SubclonalRiskScorer | Clone data, target antigen | Risk level + timeline | CAR-T escape monitoring |
| 3 | TargetExpressionValidator | Tumor/normal expression | Safety verdict + TI | CAR-T/ADC target selection |
| 4 | CellularDeconvolutionEngine | Bulk expression | Cell type proportions | Bulk-to-single-cell inference |

### 7.2 TMEClassifier Detail

**Classification logic:**
- Uses CD8+ T cell percentage, total immune score, suppressive fraction, stromal fraction
- Spatial context override: "absent" forces COLD_DESERT, "margin" forces EXCLUDED
- Score thresholds: CD8 >= 15% AND immune >= 25% for HOT_INFLAMED
- Checkpoint scoring: 6 checkpoint genes (CD274, PDCD1LG2, CTLA4, LAG3, HAVCR2, TIGIT)
- Suppressive scoring: 6 immunosuppressive genes (IDO1, TGFB1, IL10, VEGFA, ARG1, NOS2)
- Generates treatment recommendations per TME class (checkpoint inhibitors, priming strategies, stromal barrier targeting, Treg depletion)

**Evidence levels:** STRONG (spatial + PD-L1 TPS), MODERATE (one available), LIMITED (neither)

### 7.3 SubclonalRiskScorer Detail

**Risk scoring per clone:**
- Antigen-negative (expression < 0.1): +0.4
- Expanding clone: +0.2
- Proliferation index contribution: up to +0.2
- Resistance genes: +0.05 per gene (max +0.2)

**Overall risk:**
- HIGH: antigen-negative fraction > 10%
- MEDIUM: antigen-negative > 3% or any HIGH-risk clone
- LOW: all others

**Timeline estimation:** Exponential growth model -- `t = Td * log2(0.5 / neg_fraction)`

### 7.4 TargetExpressionValidator Detail

**Checks:**
- On-tumor coverage: percentage of tumor cells expressing target > threshold
- Off-tumor safety: expression in 8 vital organs (brain, heart, lung, liver, kidney, pancreas, bone_marrow, intestine)
- Therapeutic index: mean_on_tumor / (max_off_tumor + 0.01)
- Verdicts: FAVORABLE, CONDITIONAL, UNFAVORABLE

**Thresholds:**
- TI >= 10.0: favorable
- TI >= 3.0: acceptable
- TI < 3.0: unfavorable
- Tumor coverage >= 90%: excellent; >= 70%: adequate; >= 50%: marginal; < 50%: insufficient

### 7.5 CellularDeconvolutionEngine Detail

**Method:** Simplified NNLS (iterative proportional fitting)
- 10 reference cell types in default signature matrix
- 8 marker genes per cell type (80 total reference genes)
- Convergence threshold: 1e-4, max 100 iterations
- Quality metric: R-squared goodness of fit
- Minimum gene overlap: 5 genes required

**Assessment:** PASS -- all four engines are well-implemented with clear clinical logic. The deconvolution engine uses a simplified approach appropriate for real-time API response; production use for publication-grade deconvolution should integrate CIBERSORTx or MuSiC.

---

## 8. Data Model Completeness

### 8.1 Enum Coverage

| Enum | Values | Used In |
|------|--------|---------|
| SCWorkflowType | 11 | Query routing, weight selection |
| TMEClass | 4 | TME classification output |
| CellTypeConfidence | 3 | Annotation quality scoring |
| SpatialPlatform | 4 | Spatial data tagging |
| ResistanceRisk | 3 | Drug/CAR-T escape assessment |
| SeverityLevel | 5 | Clinical finding severity |
| EvidenceLevel | 5 | Evidence quality grading |
| TrajectoryType | 6 | Trajectory categorization |
| AssayType | 6 | Assay modality tagging |
| NormalizationMethod | 4 | Pipeline configuration |
| ClusteringMethod | 4 | Algorithm selection |
| AnalysisModality | 11 | Modality tracking (agent.py) |

### 8.2 Pydantic Model Coverage

| Model | Fields | Validation | Purpose |
|-------|--------|-----------|---------|
| SCQuery | 12 | max_length, ge/le, enum | Input query specification |
| SCSearchResult | 5 | ge/le on score | Search hit container |
| CellTypeAnnotation | 10 | ge/le, confidence enum | Annotation output |
| TMEProfile | 10 | ge/le, enum defaults | TME classification output |
| DrugResponsePrediction | 10 | ge/le, enum | Drug sensitivity output |
| SubclonalResult | 9 | ge constraints | Clone analysis output |
| SpatialNiche | 9 | ge/le | Spatial niche output |
| TrajectoryResult | 10 | ge constraints | Trajectory output |
| LigandReceptorInteraction | 10 | ge/le | Cell communication output |
| BiomarkerCandidate | 10 | ge/le, bool flags | Biomarker output |
| CARTTargetValidation | 10 | ge/le, enum | CAR-T safety output |
| TreatmentResponse | 9 | ge/le | Longitudinal monitoring output |
| WorkflowResult | 11 | enum, Optional | Aggregate workflow container |
| SCResponse | 9 | ge/le | Top-level API response |
| SearchPlan | 6 | dataclass | Search strategy container |

**Assessment:** PASS -- models are comprehensive with appropriate validation constraints.

---

## 9. API Surface Assessment

### 9.1 Endpoint Inventory

| Method | Path | Auth | Purpose |
|--------|------|------|---------|
| GET | `/health` | No | Service health with component status |
| GET | `/collections` | Yes | Collection names and counts |
| GET | `/workflows` | Yes | Available workflow definitions |
| GET | `/metrics` | No | Prometheus-compatible metrics |
| POST | `/v1/sc/query` | Yes | RAG Q&A query |
| POST | `/v1/sc/search` | Yes | Multi-collection vector search |
| POST | `/v1/sc/annotate` | Yes | Cell type annotation |
| POST | `/v1/sc/tme-profile` | Yes | TME profiling |
| POST | `/v1/sc/drug-response` | Yes | Drug response prediction |
| POST | `/v1/sc/subclonal` | Yes | Subclonal architecture |
| POST | `/v1/sc/spatial-niche` | Yes | Spatial niche mapping |
| POST | `/v1/sc/trajectory` | Yes | Trajectory inference |
| POST | `/v1/sc/ligand-receptor` | Yes | L-R interaction analysis |
| POST | `/v1/sc/biomarker` | Yes | Biomarker discovery |
| POST | `/v1/sc/cart-validate` | Yes | CAR-T target validation |
| POST | `/v1/sc/treatment-monitor` | Yes | Treatment monitoring |
| POST | `/v1/sc/workflow/{type}` | Yes | Generic workflow dispatch |
| GET | `/v1/sc/cell-types` | Yes | Cell type catalogue |
| GET | `/v1/sc/markers` | Yes | Marker gene reference |
| GET | `/v1/sc/tme-classes` | Yes | TME classification reference |
| GET | `/v1/sc/spatial-platforms` | Yes | Spatial platform reference |
| GET | `/v1/sc/knowledge-version` | Yes | Knowledge base version metadata |
| POST | `/v1/reports/generate` | Yes | Report generation |
| GET | `/v1/reports/formats` | Yes | Supported export formats |
| GET | `/v1/events/stream` | Yes | SSE event stream |

**Total endpoints:** 25

### 9.2 Middleware Stack

| Layer | Purpose | Config |
|-------|---------|--------|
| API key authentication | X-API-Key header validation | `settings.API_KEY` |
| Request size limiting | Reject bodies > configured MB | `settings.MAX_REQUEST_SIZE_MB` (10 MB) |
| Rate limiting | IP-based, in-memory | 100 req/60s window |
| Metrics counting | Thread-safe request/error counters | Always active |
| CORS | Configured origins | `settings.CORS_ORIGINS` |

### 9.3 API Quality Assessment

| Criterion | Status | Notes |
|-----------|--------|-------|
| OpenAPI/Swagger docs | PASS | FastAPI auto-generates at /docs |
| Input validation | PASS | Pydantic models on all POST endpoints |
| Error handling | PASS | Custom exception handlers return JSON |
| Health check | PASS | Component-level status reporting |
| Versioned routes | PASS | `/v1/` prefix on all domain routes |
| CORS configured | PASS | Non-wildcard origins |

**Assessment:** PASS.

---

## 10. Authentication and Security

### 10.1 Authentication

- **Mechanism:** API key via `X-API-Key` header or `api_key` query parameter
- **Configuration:** `settings.API_KEY` (empty string disables auth)
- **Skip paths:** `/health`, `/healthz`, `/metrics` bypass auth
- **Storage:** Environment variable `SC_API_KEY`

### 10.2 Security Measures

| Measure | Status | Notes |
|---------|--------|-------|
| API key authentication | IMPLEMENTED | Optional, header-based |
| Request size limiting | IMPLEMENTED | 10 MB default |
| Rate limiting | IMPLEMENTED | 100/min per IP |
| CORS restrictions | IMPLEMENTED | Configured origins only |
| SQL injection prevention | N/A | No SQL database |
| Non-root container user | IMPLEMENTED | `scuser` in Dockerfile |
| Input validation | IMPLEMENTED | Pydantic on all inputs |

### 10.3 Security Concerns

- API key transmitted in plaintext (requires TLS termination upstream)
- Rate limiting is in-memory (resets on restart, not shared across workers)
- No JWT/OAuth support for production multi-tenant scenarios
- `ANTHROPIC_API_KEY` passed via environment variable (standard practice but not encrypted at rest)

**Assessment:** CONDITIONAL PASS -- adequate for single-tenant deployment behind a reverse proxy with TLS. Multi-tenant production requires JWT/OAuth upgrade.

---

## 11. Embedding Pipeline

### 11.1 Configuration

| Parameter | Value |
|-----------|-------|
| Model | BAAI/bge-small-en-v1.5 |
| Dimension | 384 |
| Batch size | 32 |
| Framework | sentence-transformers 2.7.0 |

### 11.2 Embedding Quality

BGE-small-en-v1.5 is a well-validated embedding model for biomedical text:
- MTEB benchmark rank: competitive in the small model category
- Biomedical domain: adequate for gene names, cell types, and clinical terminology
- Dimension efficiency: 384-dim provides good compression vs. 768-dim models

### 11.3 Concerns

- No domain-adapted fine-tuning for single-cell terminology
- Consider upgrading to BGE-base or PubMedBERT embeddings for improved biomedical recall

**Assessment:** PASS -- adequate for current scale. Domain-specific fine-tuning is a future improvement.

---

## 12. LLM Integration

### 12.1 Configuration

| Parameter | Value |
|-----------|-------|
| Provider | Anthropic |
| Model | claude-sonnet-4-6 |
| Max tokens | 2,048 (default) |
| Temperature | 0.7 |
| System prompt | Single-cell genomics specialist |

### 12.2 LLM Usage Patterns

- **Primary:** RAG response synthesis from multi-collection search results
- **Secondary:** Workflow execution when dedicated engine is unavailable
- **Fallback:** Search-only mode when LLM is unavailable

### 12.3 Graceful Degradation

The system operates in three modes:
1. **Full mode:** Milvus + embedder + LLM all available
2. **Search-only mode:** Milvus + embedder available, no LLM
3. **Degraded mode:** Only static knowledge base available

**Assessment:** PASS -- graceful degradation is well-implemented.

---

## 13. GPU Acceleration Readiness

### 13.1 RAPIDS Integration Points

| Component | GPU Library | Status | Speedup |
|-----------|-----------|--------|---------|
| Dimensionality reduction | cuML UMAP | Planned | ~50x |
| Clustering | cuML Leiden/HDBSCAN | Planned | ~30x |
| Nearest neighbor | cuML kNN | Planned | ~100x |
| Graph analytics | cuGraph | Planned | ~40x |
| Sparse matrix ops | cuSPARSE | Planned | ~20x |

### 13.2 Configuration

- `GPU_MEMORY_LIMIT_GB`: 120 GB (DGX Spark default)
- RAPIDS integration is architecture-ready but not yet activated in the current codebase
- GPU-accelerated methods are flagged in the `sc_methods` collection schema (`gpu_accelerated` boolean field)

### 13.3 Foundation Model Integration

The knowledge base references three foundation models:
1. **scGPT** -- pre-trained on 33M cells, gene expression modeling
2. **Geneformer** -- attention-based, context-aware gene embeddings
3. **scFoundation** -- large-scale pre-trained model for cell representation

These are documented in the knowledge base for reference but not yet integrated as inference endpoints.

**Assessment:** PARTIAL -- GPU acceleration is architecturally prepared but not implemented. This is appropriate for v1.0 where the primary value is RAG-based intelligence.

---

## 14. Docker and Container Infrastructure

### 14.1 Docker Composition

| Service | Image | Port | Purpose |
|---------|-------|------|---------|
| `milvus-etcd` | quay.io/coreos/etcd:v3.5.5 | internal | Milvus metadata store |
| `milvus-minio` | minio/minio:RELEASE.2023-03-20 | internal | Milvus object storage |
| `milvus-standalone` | milvusdb/milvus:v2.4-latest | 69530, 69091 | Vector database |
| `sc-streamlit` | Custom (Dockerfile) | 8130 | Chat UI |
| `sc-api` | Custom (Dockerfile) | 8540 | REST API server |
| `sc-setup` | Custom (Dockerfile) | -- | One-shot seed script |

### 14.2 Dockerfile Assessment

| Criterion | Status | Notes |
|-----------|--------|-------|
| Multi-stage build | PASS | Builder + runtime stages |
| Non-root user | PASS | `scuser` created and used |
| Health check | PASS | curl-based on /health |
| PYTHONPATH set | PASS | `/app` |
| Minimal runtime packages | PASS | Only libgomp1 and XML libs |
| .env not copied | PASS | Via environment variables |
| Exposed ports | PASS | 8130, 8540 |

### 14.3 Docker Compose Assessment

| Criterion | Status | Notes |
|-----------|--------|-------|
| Service dependencies | PASS | `condition: service_healthy` |
| Health checks | PASS | All services have healthchecks |
| Network isolation | PASS | `sc-network` bridge |
| Named volumes | PASS | etcd_data, minio_data, milvus_data |
| Restart policies | PASS | `unless-stopped` (services), `no` (setup) |
| Environment variable pass-through | PASS | `${ANTHROPIC_API_KEY}` |

**Assessment:** PASS -- container infrastructure is production-grade.

---

## 15. Seed Data and Ingest Pipeline

### 15.1 Ingest Pipeline Components

| Parser | Source | Output Collection | Seed Records |
|--------|--------|------------------|-------------|
| `cellxgene_parser.py` | CellxGene/HCA references | sc_cell_types, sc_datasets | 49 |
| `marker_parser.py` | CellMarker/PanglaoDB | sc_markers | 75 |
| `tme_parser.py` | Cancer TME atlas | sc_tme | 20 |

### 15.2 Ingest Architecture

All parsers extend `BaseIngestParser`:
- `parse()` -- extract records from source
- `transform()` -- normalize to collection schema
- `load()` -- embed and insert into Milvus

### 15.3 Seed Script

The `seed_knowledge.py` script:
1. Creates all 12 collections via `setup_collections.py --drop-existing --seed`
2. Generates embeddings via sentence-transformers
3. Inserts seed records into Milvus
4. Gracefully degrades if pymilvus or sentence-transformers unavailable

### 15.4 Scheduled Ingest

- APScheduler integration via `src/scheduler.py`
- Configurable interval: `INGEST_SCHEDULE_HOURS` (default 24)
- Disabled by default: `INGEST_ENABLED = False`

**Assessment:** PASS -- seed pipeline is functional. Scheduled ingest is available but disabled pending production data source integration.

---

## 16. Cross-Agent Integration

### 16.1 Integration Map

| Peer Agent | URL | Timeout | Use Case |
|-----------|-----|---------|----------|
| Genomics Agent | `http://localhost:8527` | 30s | Variant-level evidence for mutations |
| Biomarker Agent | `http://localhost:8529` | 30s | Cross-modal biomarker correlation |
| Oncology Agent | `http://localhost:8528` | 30s | Therapy line recommendations |
| Clinical Trial Agent | `http://localhost:8538` | 30s | Trial matching for novel targets |

### 16.2 Integration Mechanism

- `src/cross_modal.py` handles inter-agent REST calls
- Timeout: 30 seconds per agent
- Failure mode: graceful degradation (result returned without cross-agent enrichment)
- Shared collection: `genomic_evidence` (read-only from genomics pipeline)

### 16.3 Concerns

- No circuit breaker -- repeated failures to a down agent will cause 30s timeouts per request
- No retry with backoff
- No service discovery -- hardcoded localhost URLs

**Assessment:** CONDITIONAL PASS -- functional for co-located deployment. Production distributed deployment requires service discovery and circuit breaker patterns.

---

## 17. Test Coverage Analysis

### 17.1 Test File Inventory

| Test File | Lines | Focus |
|-----------|-------|-------|
| `test_models.py` | 364 | Pydantic model validation (all 15 models) |
| `test_decision_support.py` | 230 | 4 decision engines with edge cases |
| `test_clinical_workflows.py` | 228 | 10 workflow execution paths |
| `test_knowledge.py` | 140 | Knowledge base completeness |
| `test_integration.py` | 168 | End-to-end RAG + workflow |
| `test_collections.py` | 115 | Schema validation, weight sums |
| `test_api.py` | 107 | API endpoint testing |
| `test_settings.py` | 86 | Configuration validation |
| `test_query_expansion.py` | 63 | Synonym expansion |
| `test_workflow_execution.py` | 152 | Workflow dispatch |
| `test_agent.py` | 49 | Agent plan-execute-report |
| `test_rag_engine.py` | 43 | RAG search and synthesis |
| `conftest.py` | 15 | Shared fixtures |

**Total:** 1,760 lines across 12 test files, approximately **185 test cases**.

### 17.2 Coverage Assessment

| Module | Test Coverage | Quality |
|--------|-------------|---------|
| Models | HIGH | All 15 models validated |
| Decision Support | HIGH | All 4 engines with edge cases |
| Clinical Workflows | HIGH | All 10 workflows tested |
| Collections | MEDIUM | Schema and weight validation |
| Knowledge Base | MEDIUM | Completeness checks |
| API | MEDIUM | Endpoint smoke tests |
| Settings | MEDIUM | Validation rule testing |
| RAG Engine | LOW | Basic search test only |
| Agent | LOW | Basic plan-execute test |
| Ingest | NOT TESTED | No dedicated tests |
| Cross-Modal | NOT TESTED | No dedicated tests |
| Scheduler | NOT TESTED | No dedicated tests |

### 17.3 Estimated Line Coverage

Based on test file analysis: approximately **65-70% line coverage** for core modules, **40-50% overall** including API/ingest/scheduler.

**Assessment:** CONDITIONAL PASS -- critical paths (models, decision support, workflows) are well-tested. RAG engine, ingest pipeline, and cross-agent integration need additional test coverage.

---

## 18. Performance Benchmarks

### 18.1 Expected Latency Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| Single-collection search | < 50ms | IVF_FLAT with nlist=128 |
| Multi-collection search (12) | < 500ms | ThreadPoolExecutor parallel |
| LLM synthesis | 2-5s | Claude Sonnet response time |
| Full RAG query | < 6s | Search + synthesis |
| Workflow execution | < 8s | Workflow + RAG + synthesis |
| Health check | < 100ms | Component status check |
| Report generation | < 15s | DOCX/PDF export |

### 18.2 Throughput Targets

| Metric | Target |
|--------|--------|
| Concurrent queries | 10 (2 uvicorn workers) |
| Queries per minute | 30-50 (rate limited to 100/min) |
| Seed records per minute | 500+ (batch embedding) |

### 18.3 Resource Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 4 cores | 8 cores |
| RAM | 8 GB | 16 GB |
| GPU VRAM | None (CPU inference) | 8 GB (RAPIDS acceleration) |
| Disk | 10 GB | 50 GB (with full collections) |
| Network | 100 Mbps | 1 Gbps |

**Assessment:** PASS -- targets are achievable on the DGX Spark platform.

---

## 19. Monitoring and Observability

### 19.1 Metrics Endpoints

- `/metrics` -- Prometheus-compatible text format
- Counters: requests_total, query_requests_total, search_requests_total, annotate_requests_total, workflow_requests_total, report_requests_total, errors_total
- `src/metrics.py` -- dedicated metrics module (476 lines) for extended metrics

### 19.2 Logging

- **Framework:** Loguru (structured logging)
- **Log levels:** DEBUG through CRITICAL
- **Output:** stdout (Docker logs compatible)
- **Correlation:** No request-ID tracing (improvement needed)

### 19.3 Health Monitoring

The `/health` endpoint reports:
- Overall status: healthy / degraded
- Component status: milvus, rag_engine, workflow_engine
- Collection count and total vector count
- Workflow count

### 19.4 Concerns

- No distributed tracing (OpenTelemetry)
- No request-ID propagation for cross-agent debugging
- No alerting integration (Prometheus + AlertManager recommended)

**Assessment:** CONDITIONAL PASS -- basic monitoring is in place. Production requires distributed tracing and alerting.

---

## 20. Configuration Management

### 20.1 Configuration Source

Pydantic BaseSettings with layered configuration:
1. Default values in `SingleCellSettings` class
2. Environment variables with `SC_` prefix
3. `.env` file support

### 20.2 Key Configuration Parameters

| Parameter | Default | Env Var |
|-----------|---------|---------|
| MILVUS_HOST | localhost | SC_MILVUS_HOST |
| MILVUS_PORT | 19530 | SC_MILVUS_PORT |
| API_PORT | 8540 | SC_API_PORT |
| STREAMLIT_PORT | 8130 | SC_STREAMLIT_PORT |
| EMBEDDING_MODEL | BAAI/bge-small-en-v1.5 | SC_EMBEDDING_MODEL |
| LLM_MODEL | claude-sonnet-4-6 | SC_LLM_MODEL |
| SCORE_THRESHOLD | 0.4 | SC_SCORE_THRESHOLD |
| GPU_MEMORY_LIMIT_GB | 120 | SC_GPU_MEMORY_LIMIT_GB |
| MAX_REQUEST_SIZE_MB | 10 | SC_MAX_REQUEST_SIZE_MB |
| INGEST_SCHEDULE_HOURS | 24 | SC_INGEST_SCHEDULE_HOURS |
| CROSS_AGENT_TIMEOUT | 30 | SC_CROSS_AGENT_TIMEOUT |

### 20.3 Startup Validation

The `SingleCellSettings.validate()` method checks:
- MILVUS_HOST is non-empty
- MILVUS_PORT is in valid range (1-65535)
- ANTHROPIC_API_KEY is set (warns if not)
- EMBEDDING_MODEL is non-empty
- API_PORT and STREAMLIT_PORT are valid and non-conflicting
- Collection search weights sum to ~1.0 (tolerance 0.05)

**Assessment:** PASS -- configuration management is well-structured with validation.

---

## 21. Error Handling and Resilience

### 21.1 Error Handling Patterns

| Pattern | Implementation | Location |
|---------|---------------|----------|
| HTTP exception handler | Custom JSON response with agent name | api/main.py |
| General exception handler | 500 with error logging | api/main.py |
| Milvus connection failure | Graceful degradation, manager = None | api/main.py lifespan |
| Embedding model missing | embedder = None, search disabled | api/main.py lifespan |
| LLM unavailable | search-only mode | api/main.py lifespan |
| Cross-agent timeout | 30s timeout, return without enrichment | src/cross_modal.py |
| Invalid collection name | KeyError with valid collection list | src/collections.py |
| Insufficient gene overlap | Warning return with quality metrics | decision_support.py |

### 21.2 Resilience Assessment

| Criterion | Status |
|-----------|--------|
| Graceful degradation | PASS |
| Component isolation | PASS |
| Error logging | PASS |
| No silent failures | PASS |
| Timeout handling | PARTIAL (no circuit breaker) |
| Retry with backoff | NOT IMPLEMENTED |

**Assessment:** CONDITIONAL PASS -- resilience is good for single-service deployment. Circuit breaker and retry patterns needed for distributed operation.

---

## 22. Streamlit UI Assessment

### 22.1 UI Features

The Streamlit interface provides 5 tabs:

1. **Chat** -- RAG-powered Q&A with conversation history
2. **TME Profiler** -- Interactive TME classification with cell proportion inputs
3. **Workflows** -- Workflow selection and execution
4. **Dashboard** -- Real-time health monitoring and metrics
5. **Reference** -- Cell type and marker gene catalogues

### 22.2 Theming

- NVIDIA dark theme with branded green accent (#76b900)
- Custom CSS for consistent styling
- Responsive layout using `st.columns()`

### 22.3 UI Quality

| Criterion | Status |
|-----------|--------|
| All 10 workflows accessible | PASS |
| API integration | PASS (via SC_API_BASE) |
| Error messaging | PASS |
| Loading states | PASS |
| Theme consistency | PASS |

**Assessment:** PASS.

---

## 23. Known Issues and Technical Debt

### 23.1 Critical Issues

| # | Issue | Impact | Mitigation |
|---|-------|--------|-----------|
| 1 | Dual SCWorkflowType enum definitions | Type confusion, import ambiguity | Consolidate to single source in models.py |
| 2 | No circuit breaker for cross-agent calls | 30s timeout cascading under agent failure | Implement tenacity or pybreaker |

### 23.2 High-Priority Issues

| # | Issue | Impact | Mitigation |
|---|-------|--------|-----------|
| 3 | Rate limiting is in-memory | Resets on restart, not shared across workers | Move to Redis-backed rate limiting |
| 4 | No request-ID tracing | Difficult cross-agent debugging | Add OpenTelemetry trace context |
| 5 | Ingest pipeline has no tests | Silent regressions in data loading | Add test_ingest.py |
| 6 | RAG engine test coverage is minimal | Core search logic undertested | Expand test_rag_engine.py |
| 7 | No schema migration strategy | Collection schema changes require drop/recreate | Implement migration versioning |

### 23.3 Medium-Priority Issues

| # | Issue | Impact | Mitigation |
|---|-------|--------|-----------|
| 8 | RAPIDS/GPU acceleration not implemented | Missing 30-100x speedup for large datasets | Implement cuML integration |
| 9 | Foundation model integration not implemented | scGPT/Geneformer not available for inference | Add NIM endpoints |
| 10 | Conversation cleanup lacks scheduling | Disk space growth over time | Add TTL-based cleanup cron |
| 11 | Deconvolution engine is simplified | Not suitable for publication-grade analysis | Integrate CIBERSORTx |
| 12 | No A/B testing framework | Cannot compare model/prompt changes | Add experiment tracking |

### 23.4 Low-Priority Issues

| # | Issue | Impact |
|---|-------|--------|
| 13 | No DOCX template customization | Reports use default styling |
| 14 | No batch query endpoint | Single-query only |
| 15 | No webhook/callback support | Polling required for long workflows |

---

## 24. Risk Register

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|-----------|
| LLM API outage | Medium | High -- no synthesis capability | Search-only fallback implemented |
| Milvus data corruption | Low | Critical -- all search disabled | Regular backup schedule, replicated deployment |
| Embedding model drift | Low | Medium -- search quality degradation | Version-pinned model, periodic evaluation |
| Cross-agent cascade failure | Medium | Medium -- degraded but functional | Circuit breaker (not yet implemented) |
| API key exposure | Low | High -- unauthorized access | Environment variable injection, no hardcoded keys |
| Knowledge base staleness | Medium | Medium -- outdated recommendations | Scheduled ingest (configurable interval) |
| Rate limit bypass | Low | Medium -- resource exhaustion | Upgrade to Redis-backed rate limiting |

---

## 25. Go/No-Go Recommendation

### 25.1 Summary Scores

| Category | Score | Status |
|----------|-------|--------|
| Architecture | 9/10 | PASS |
| Data model | 9/10 | PASS |
| API surface | 8/10 | PASS |
| Knowledge base | 9/10 | PASS |
| Decision support | 9/10 | PASS |
| Workflow engine | 8/10 | PASS |
| Container infrastructure | 9/10 | PASS |
| Security | 6/10 | CONDITIONAL |
| Test coverage | 6/10 | CONDITIONAL |
| Monitoring | 6/10 | CONDITIONAL |
| GPU acceleration | 4/10 | NOT READY |
| Cross-agent integration | 6/10 | CONDITIONAL |
| **Overall** | **7.4/10** | **CONDITIONAL PASS** |

### 25.2 Pre-Production Conditions

The following must be addressed before production deployment:

1. **MUST:** Consolidate dual SCWorkflowType enum definitions
2. **MUST:** Deploy behind TLS-terminating reverse proxy
3. **MUST:** Set non-empty API_KEY in production environment
4. **SHOULD:** Add circuit breaker for cross-agent calls
5. **SHOULD:** Increase RAG engine test coverage to > 80%
6. **SHOULD:** Add ingest pipeline tests
7. **SHOULD:** Implement request-ID tracing
8. **NICE:** Integrate Redis-backed rate limiting
9. **NICE:** Add OpenTelemetry distributed tracing
10. **NICE:** Implement RAPIDS GPU acceleration

### 25.3 Recommendation

**APPROVED FOR PRODUCTION** with the three MUST conditions above. The Single-Cell Intelligence Agent demonstrates strong architectural quality, comprehensive domain knowledge, and well-implemented clinical decision support engines. The conditional items should be addressed in the next sprint cycle.

---

*Report generated: 2026-03-22*
*HCLS AI Factory -- Single-Cell Intelligence Agent v1.0.0*
