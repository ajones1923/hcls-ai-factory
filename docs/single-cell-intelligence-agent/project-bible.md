# Single-Cell Intelligence Agent -- Project Bible

**Date:** 2026-03-22
**Author:** Adam Jones
**Classification:** Internal Reference Document

---

## 1. Mission Statement

The Single-Cell Intelligence Agent delivers RAG-powered clinical decision support at single-cell resolution, transforming raw transcriptomic data into clinically actionable insights for oncology, immunology, and cell therapy teams. It bridges the gap between single-cell research output and bedside treatment decisions by combining curated domain knowledge, vector-based evidence retrieval, and LLM-powered synthesis.

---

## 2. Problem Statement

### 2.1 The Resolution Gap

Precision medicine has advanced from tissue-level to gene-level analysis, but clinical decisions still rely on bulk-averaged measurements that mask critical cellular heterogeneity. A tumor biopsy reporting "PD-L1 positive" may contain 30% immune-hot niches and 70% immune-cold desert -- information invisible at bulk resolution but decisive for immunotherapy selection.

### 2.2 The Interpretation Bottleneck

Single-cell RNA-seq generates datasets with 10,000-500,000 cells, each expressing 20,000+ genes. A trained bioinformatician requires 2-4 weeks to fully annotate, profile, and interpret a single dataset. Clinical turnaround expectations are 24-72 hours.

### 2.3 The Knowledge Integration Challenge

Actionable single-cell interpretation requires simultaneous access to:
- 44+ cell type definitions with canonical marker genes
- 12 cancer-specific TME reference profiles
- 30+ drug sensitivity databases
- Spatial transcriptomics platform specifications
- Active clinical trial registries
- 75+ validated marker gene associations

No single analyst maintains current expertise across all these domains.

---

## 3. Solution Architecture

### 3.1 High-Level Design

The agent follows the **Plan-Search-Evaluate-Synthesize-Report** pattern:

1. **Plan:** Parse the query, classify the workflow type (1 of 11), and construct a search plan with collection-specific weights
2. **Search:** Execute parallel vector searches across 12 Milvus collections using BGE-small-en-v1.5 embeddings
3. **Evaluate:** Score evidence quality, check cross-collection corroboration, assess clinical relevance
4. **Synthesize:** Generate a grounded LLM response using Claude Sonnet with retrieved evidence
5. **Report:** Format structured output with citations, severity levels, and actionable recommendations

### 3.2 Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| 12 separate collections (not 1 monolithic) | Workflow-specific weight boosting, schema optimization per data type |
| BGE-small-en-v1.5 (384-dim) | Balance of quality and speed; 384-dim sufficient for biomedical text |
| IVF_FLAT index | Simple, accurate, adequate for < 100K records per collection |
| Pydantic models for all I/O | Type safety, validation, auto-documentation |
| 4 dedicated decision engines | Deterministic clinical logic separate from LLM stochasticity |
| Graceful degradation | Each component failure reduces capability, never crashes |
| COSINE similarity | Standard for normalized text embeddings |

### 3.3 What This Agent Does NOT Do

- Does not process raw FASTQ/BAM files (that is the genomics pipeline)
- Does not perform de novo clustering or dimensionality reduction (that is the computational pipeline)
- Does not store patient health records (no PHI storage)
- Does not replace pathologist review (decision support, not decision making)
- Does not run GPU-accelerated analysis (v1.0 -- RAPIDS integration planned for v2.0)

---

## 4. Stakeholder Map

| Stakeholder | Role | Interest |
|------------|------|----------|
| Clinical oncologist | End user | TME classification, drug response, treatment monitoring |
| Bioinformatician | Power user | Cell type annotation, trajectory inference, method selection |
| Cell therapy team | End user | CAR-T target validation, escape risk assessment |
| Spatial biology researcher | Power user | Spatial niche mapping, L-R interaction analysis |
| Clinical trial coordinator | Consumer | Biomarker discovery, trial-eligible target identification |
| Platform engineering | Operator | Deployment, monitoring, scaling |

---

## 5. Feature Inventory

### 5.1 Analysis Workflows (10)

| # | Workflow | Clinical Question Answered |
|---|---------|---------------------------|
| 1 | Cell Type Annotation | "What cell types are in this sample and at what proportions?" |
| 2 | TME Profiling | "Is this tumor hot, cold, excluded, or immunosuppressive?" |
| 3 | Drug Response | "Which drugs will this tumor respond to at the cellular level?" |
| 4 | Subclonal Architecture | "Are there resistant subclones that could cause relapse?" |
| 5 | Spatial Niche | "Where are the immune cells relative to tumor cells in tissue?" |
| 6 | Trajectory Analysis | "What differentiation or exhaustion trajectories are active?" |
| 7 | Ligand-Receptor | "Which cell-cell communication axes are driving tumor progression?" |
| 8 | Biomarker Discovery | "What cell-type-specific biomarkers predict treatment outcome?" |
| 9 | CAR-T Validation | "Is this target safe and effective for CAR-T therapy?" |
| 10 | Treatment Monitoring | "How has the tumor composition changed under treatment?" |

### 5.2 Decision Support Engines (4)

| Engine | Input | Output | Deterministic |
|--------|-------|--------|--------------|
| TMEClassifier | Cell proportions, gene expression | TME class + treatment recs | Yes |
| SubclonalRiskScorer | Clone data, target antigen | Risk level + timeline | Yes |
| TargetExpressionValidator | Tumor/normal expression | Safety verdict | Yes |
| CellularDeconvolutionEngine | Bulk expression | Cell type proportions | Yes |

### 5.3 Knowledge Resources

| Resource | Records | Update Frequency |
|----------|---------|-----------------|
| Cell Type Atlas | 44 types, 232 aliases | Static (v3.0.0) |
| Drug Database | 30 drugs, 10 classes | Semi-annual |
| Marker Genes | 75 genes | Static |
| Immune Signatures | 10 signatures | Static |
| L-R Pairs | 25 pairs | Static |
| Cancer TME Atlas | 12 cancer types | Semi-annual |
| CellxGene Seeds | 49 records | On ingest |
| Marker Seeds | 75 records | On ingest |
| TME Seeds | 20 records | On ingest |

---

## 6. Technical Stack

### 6.1 Core Technologies

| Layer | Technology | Version |
|-------|-----------|---------|
| Language | Python | 3.10 |
| API framework | FastAPI | 0.111.0 |
| UI framework | Streamlit | 1.33.0 |
| Vector database | Milvus | 2.4 |
| Embedding model | BGE-small-en-v1.5 | 384-dim |
| LLM | Claude Sonnet (Anthropic) | claude-sonnet-4-6 |
| Containerization | Docker | Multi-stage |
| Orchestration | Docker Compose | 3.8 |
| Validation | Pydantic | 2.7.4 |
| Metrics | Prometheus client | 0.20.0 |
| Scheduling | APScheduler | 3.10.4 |
| Export | python-docx | 1.1.0 |

### 6.2 Data Sources

| Source | Data Type | Integration |
|--------|----------|-------------|
| Human Cell Atlas | Cell type references | Seed data |
| CellMarker 2.0 | Marker-cell associations | Seed data |
| Cell Ontology (CL) | Ontology identifiers | Static mapping |
| PanglaoDB | Marker gene database | Seed data |
| CellxGene | Dataset metadata | API ingest |
| GDSC | Drug sensitivity | Knowledge base |
| DepMap | Cancer dependency | Knowledge base |
| ClinVar | Variant classification | Shared collection |
| TISCH2 | TME atlas | Knowledge base |
| CellPhoneDB | L-R interactions | Knowledge base |

---

## 7. Port Allocation

| Port | Service | Protocol |
|------|---------|----------|
| 8540 | FastAPI REST API | HTTP |
| 8130 | Streamlit UI | HTTP |
| 19530 | Milvus (shared) | gRPC |
| 69530 | Milvus (standalone) | gRPC |
| 69091 | Milvus health (standalone) | HTTP |

---

## 8. Data Flow

```
User Query
    |
    v
Query Classification (SCWorkflowType)
    |
    v
Search Plan Construction
    |-- Collection selection (12 collections)
    |-- Weight profile selection (11 profiles)
    |-- Filter expression generation
    |
    v
Parallel Vector Search (ThreadPoolExecutor)
    |-- sc_cell_types (weight: 0.14)
    |-- sc_markers (weight: 0.12)
    |-- ... (10 more collections)
    |
    v
Evidence Aggregation & Scoring
    |-- Cross-collection entity linking
    |-- Citation relevance scoring
    |-- Evidence level assessment
    |
    v
LLM Synthesis (Claude Sonnet)
    |-- System prompt: single-cell specialist
    |-- Context: top-K evidence from search
    |-- Conversation history (3-turn window)
    |
    v
Structured Response (SCResponse)
    |-- answer (natural language)
    |-- workflow_result (typed)
    |-- citations (formatted)
    |-- confidence (0-1)
```

---

## 9. Quality Gates

### 9.1 Code Quality

| Gate | Tool | Threshold |
|------|------|-----------|
| Type safety | Pydantic validation | All I/O models typed |
| Unit tests | pytest | 185+ test cases |
| Configuration validation | SingleCellSettings.validate() | 0 critical warnings |
| Weight sum validation | Collections.py | Sum within 0.05 of 1.0 |

### 9.2 Clinical Quality

| Gate | Mechanism |
|------|-----------|
| TME classification accuracy | Validated against TISCH2 reference profiles |
| Drug sensitivity correlation | Cross-referenced with GDSC IC50 data |
| CAR-T safety thresholds | Based on published vital organ expression data |
| Evidence grading | Four-tier evidence level system |
| Severity classification | Five-level clinical severity scale |

---

## 10. Deployment Models

### 10.1 Standalone (Docker Compose)

Includes dedicated Milvus instance (etcd + MinIO + standalone). Suitable for development, testing, and single-user deployment.

### 10.2 Integrated (DGX Spark)

Connects to shared Milvus instance via the top-level `docker-compose.dgx-spark.yml`. Reads from shared `genomic_evidence` collection. Suitable for production deployment alongside other HCLS AI Factory agents.

### 10.3 VAST AI OS

Deployed as a function within the VAST AI OS platform with automatic scaling and health monitoring. Uses the VAST AI OS AgentEngine model for lifecycle management.

---

## 11. Roadmap

### 11.1 v1.0 (Current)

- 12 Milvus collections with seed data
- 10 analysis workflows
- 4 decision support engines
- FastAPI + Streamlit deployment
- Cross-agent integration (4 peer agents)
- 185 test cases

### 11.2 v1.1 (Q2 2026)

- RAPIDS GPU acceleration for cuML UMAP/clustering
- scGPT foundation model integration via NIM
- CIBERSORTx-grade deconvolution
- OpenTelemetry distributed tracing
- Redis-backed rate limiting

### 11.3 v2.0 (Q3 2026)

- Multi-modal integration (scATAC-seq, CITE-seq, Multiome)
- Real-time spatial analysis pipeline
- Automated report generation with institutional templates
- Patient longitudinal tracking dashboard
- FDA 21 CFR Part 11 compliance features

---

## 12. Glossary

| Term | Definition |
|------|-----------|
| AnnData | Annotated data matrix format for single-cell data (.h5ad) |
| CAR-T | Chimeric Antigen Receptor T-cell therapy |
| CITE-seq | Cellular Indexing of Transcriptomes and Epitopes by Sequencing |
| CL | Cell Ontology (standardized cell type identifier system) |
| DE | Differential Expression analysis |
| GDSC | Genomics of Drug Sensitivity in Cancer |
| HCA | Human Cell Atlas consortium |
| L-R | Ligand-Receptor (cell-cell communication) |
| MERFISH | Multiplexed Error-Robust FISH (spatial transcriptomics) |
| MRD | Minimal Residual Disease |
| NNLS | Non-Negative Least Squares (deconvolution method) |
| PCA | Principal Component Analysis |
| RAG | Retrieval-Augmented Generation |
| scRNA-seq | Single-cell RNA sequencing |
| TME | Tumor Microenvironment |
| UMAP | Uniform Manifold Approximation and Projection |
| Visium | 10x Genomics spatial transcriptomics platform |
| Xenium | 10x Genomics in situ spatial platform |

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Project Bible*
