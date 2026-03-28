# Single-Cell Intelligence Agent -- Architecture Guide

**Date:** March 2026
**Author:** Adam Jones

---

## 1. Architecture Overview

The Single-Cell Intelligence Agent is built on a layered architecture that separates presentation, API, reasoning, search, clinical decision support, and data storage concerns.

```
+----------------------------------------------------------+
|                  PRESENTATION LAYER                       |
|  Streamlit UI (8130)           External Clients           |
+----------------------------------------------------------+
                           |
+----------------------------------------------------------+
|                     API LAYER                             |
|  FastAPI REST (8540)                                     |
|  - Auth middleware      - Rate limiting                   |
|  - Request validation   - CORS                            |
|  - Error handling       - Metrics collection              |
+----------------------------------------------------------+
                           |
+----------------------------------------------------------+
|                  REASONING LAYER                          |
|  SingleCellAgent                                         |
|  - Plan: Query classification, search strategy            |
|  - Execute: Multi-collection RAG search                   |
|  - Evaluate: Evidence quality scoring                     |
|  - Synthesize: LLM-powered response generation            |
|  - Report: Structured output with citations               |
+----------------------------------------------------------+
                           |
          +----------------+----------------+
          |                |                |
+---------+--+    +--------+---+    +------+------+
| RAG ENGINE |    | WORKFLOW   |    | DECISION    |
| Parallel   |    | ENGINE     |    | SUPPORT     |
| search     |    | 10 clinical|    | 4 engines   |
| synthesis  |    | workflows  |    |             |
+------+-----+    +------+-----+    +------+------+
       |                 |                |
+------+-----------------+----------------+------+
|                  DATA LAYER                     |
|  Milvus Vector DB (12 collections)              |
|  Knowledge Base (57 cell types, 30 drugs, ...)  |
|  Conversation Store (disk-backed, 24h TTL)      |
+-------------------------------------------------+
```

---

## 2. Module Architecture

### 2.1 Core Modules

```
single_cell_intelligence_agent/
|-- config/
|   |-- settings.py              # Pydantic BaseSettings (197 lines)
|-- src/
|   |-- agent.py                 # Autonomous reasoning engine (2,090 lines)
|   |-- models.py                # Pydantic data models (820 lines)
|   |-- collections.py           # 12 Milvus collection schemas (1,210 lines)
|   |-- rag_engine.py            # Multi-collection RAG search (1,490 lines)
|   |-- clinical_workflows.py    # 10 analysis workflows (1,792 lines)
|   |-- decision_support.py      # 4 clinical engines (886 lines)
|   |-- knowledge.py             # Domain knowledge base (1,816 lines)
|   |-- query_expansion.py       # Synonym expansion (893 lines)
|   |-- cross_modal.py           # Inter-agent communication (392 lines)
|   |-- metrics.py               # Prometheus metrics (476 lines)
|   |-- export.py                # Report generation (588 lines)
|   |-- scheduler.py             # APScheduler ingest (496 lines)
|   |-- ingest/
|       |-- base.py              # BaseIngestParser ABC (228 lines)
|       |-- cellxgene_parser.py  # CellxGene data parser (679 lines)
|       |-- marker_parser.py     # Marker gene parser (286 lines)
|       |-- tme_parser.py        # TME profile parser (418 lines)
|-- api/
|   |-- main.py                  # FastAPI application (615 lines)
|   |-- routes/
|       |-- sc_clinical.py       # Clinical endpoint routes
|       |-- reports.py           # Report generation routes
|       |-- events.py            # SSE event stream routes
|-- app/
|   |-- sc_ui.py                 # Streamlit 5-tab UI (~600 lines)
```

### 2.2 Dependency Graph

```
settings.py
    |
    v
models.py <---- collections.py
    |                |
    v                v
agent.py -------> rag_engine.py
    |                |
    +--------+-------+
    |        |       |
    v        v       v
clinical_  decision_ knowledge.py
workflows  support
    |        |
    v        v
    api/main.py
        |
        v
    app/sc_ui.py
```

---

## 3. GPU Acceleration Pipeline

### 3.1 RAPIDS Integration Architecture

The agent is architecturally prepared for GPU acceleration via NVIDIA RAPIDS. The integration targets five computational bottlenecks:

```
Single-Cell Data (AnnData .h5ad)
         |
    +----+----+
    |         |
  CPU Path   GPU Path (RAPIDS)
    |         |
    v         v
  scikit-   cuML
  learn     UMAP/PCA
    |         |
    v         v
  scanpy    cuGraph
  Leiden    Leiden/Louvain
    |         |
    v         v
  scipy     cuSPARSE
  sparse    operations
    |         |
    +----+----+
         |
         v
    Annotation + Analysis
```

### 3.2 RAPIDS Component Mapping

| CPU Library | RAPIDS GPU | Operation | Expected Speedup |
|------------|-----------|-----------|-----------------|
| sklearn.decomposition.PCA | cuml.PCA | Dimensionality reduction | 30-50x |
| umap-learn | cuml.UMAP | Manifold embedding | 50-100x |
| sklearn.neighbors.NearestNeighbors | cuml.NearestNeighbors | kNN graph | 80-120x |
| igraph/leidenalg | cugraph.leiden | Community detection | 20-40x |
| scipy.sparse | cuSPARSE | Sparse matrix ops | 15-25x |
| sklearn.cluster.KMeans | cuml.KMeans | K-means clustering | 40-60x |

### 3.3 GPU Memory Management

```python
# Configuration from settings.py
GPU_MEMORY_LIMIT_GB = 120  # DGX Spark default

# Planned memory allocation strategy:
# - 40% for RAPIDS cuML operations (48 GB)
# - 30% for Milvus GPU index (36 GB)
# - 20% for foundation model inference (24 GB)
# - 10% for system overhead (12 GB)
```

### 3.4 GPU-Accelerated Methods Registry

The `sc_methods` collection includes a `gpu_accelerated` boolean field to track which analytical methods support GPU execution:

| Method | GPU Support | Library |
|--------|-----------|---------|
| UMAP | Yes | cuml.UMAP |
| Leiden clustering | Yes | cugraph.leiden |
| PCA | Yes | cuml.PCA |
| t-SNE | Yes | cuml.TSNE |
| kNN graph | Yes | cuml.NearestNeighbors |
| Differential expression | Partial | RAPIDS cuDF |
| Trajectory inference | No | Monocle3 (R) / scVelo |
| CellChat | No | R-based |
| Scanpy preprocessing | Partial | rapids-singlecell |

---

## 4. TME Classification Architecture

### 4.1 Classification Pipeline

```
Input: Cell Type Proportions + Gene Expression
              |
              v
    +-------------------+
    | Immune Score       |  Sum of 8 immune cell type fractions
    | (CD8_T, CD4_T,    |  (CD8_T, CD4_T, NK, B_cell,
    |  NK, B_cell, ...)  |   Macrophage_M1, Dendritic,
    +-------------------+   Plasma, Neutrophil)
              |
              v
    +-------------------+
    | Suppressive Score  |  Weighted combination of:
    | (Treg, MDSC,      |  - Suppressive cell fraction (50%)
    |  M2 Macrophage)   |  - Suppressive gene score (50%)
    +-------------------+    (IDO1, TGFB1, IL10, VEGFA, ARG1, NOS2)
              |
              v
    +-------------------+
    | Checkpoint Score   |  6 checkpoint genes normalized:
    | (CD274, PDCD1LG2, |  CD274, PDCD1LG2, CTLA4,
    |  CTLA4, LAG3, ... |  LAG3, HAVCR2, TIGIT
    +-------------------+
              |
              v
    +-------------------+
    | Spatial Override   |  "absent" -> COLD_DESERT
    |                   |  "margin" -> EXCLUDED
    +-------------------+
              |
              v
    +-------------------+
    | Classification    |
    | Decision Tree     |
    +-------------------+
              |
    +---------+---------+-----------+
    |         |         |           |
    v         v         v           v
HOT       COLD      EXCLUDED  IMMUNO-
INFLAMED  DESERT               SUPPRESSIVE
```

### 4.2 Classification Decision Tree

```
IF spatial == "absent" AND immune < 0.05:
    -> COLD_DESERT

IF spatial == "margin" AND immune > 0.05:
    -> EXCLUDED

IF CD8 >= 0.15 AND immune >= 0.25:
    IF suppressive > 0.4:
        -> IMMUNOSUPPRESSIVE
    ELSE:
        -> HOT_INFLAMED

IF immune >= 0.10 AND stromal > 0.20:
    -> EXCLUDED

IF suppressive > 0.3 AND immune >= 0.10:
    -> IMMUNOSUPPRESSIVE

IF immune < 0.10:
    -> COLD_DESERT

IF PD-L1_high AND CD8 >= 0.05:
    -> HOT_INFLAMED

DEFAULT:
    -> COLD_DESERT
```

### 4.3 Treatment Recommendation Engine

Each TME class maps to a set of evidence-based treatment recommendations:

| TME Class | Primary Recommendation | Conditional Recommendations |
|-----------|----------------------|---------------------------|
| HOT_INFLAMED | Checkpoint inhibitor (anti-PD-1/PD-L1) | PD-L1 TPS >= 50%: pembrolizumab mono; LAG3+: relatlimab + nivolumab |
| COLD_DESERT | Priming strategies (oncolytic virus, STING agonist) | BiTE or adoptive cell therapy |
| EXCLUDED | Anti-TGFb or anti-VEGF to remove stromal barrier | Anti-CXCL12/CXCR4 for T-cell migration |
| IMMUNOSUPPRESSIVE | Dual checkpoint (anti-PD-1 + anti-CTLA-4) | Anti-CCR8 for Treg depletion; CSF1R inhibitor for M2 repolarization |

---

## 5. Spatial Deconvolution Architecture

### 5.1 Spatial Platform Support

| Platform | Resolution | Genes | Spatial Feature |
|----------|-----------|-------|----------------|
| Visium (10x) | 55 um spots | Whole transcriptome | H&E morphology overlay |
| MERFISH (Vizgen) | Subcellular | 100-500 panel | Single-molecule FISH |
| Xenium (10x) | Subcellular | 100-5000 panel | In situ sequencing |
| CODEX | Single-cell | 40-60 proteins | Protein co-detection |

### 5.2 Spatial Niche Detection Pipeline

```
Spatial Coordinates + Gene Expression
              |
              v
    +-------------------+
    | Cell Type          |  Assign cell types to spatial
    | Annotation         |  locations using marker genes
    +-------------------+
              |
              v
    +-------------------+
    | Spatial            |  Moran's I statistic for
    | Autocorrelation    |  spatially variable genes
    +-------------------+
              |
              v
    +-------------------+
    | Niche              |  k-NN graph on spatial
    | Construction       |  coordinates, community
    +-------------------+  detection on cell types
              |
              v
    +-------------------+
    | L-R Interaction    |  Spatially-aware ligand-
    | Analysis           |  receptor scoring
    +-------------------+
              |
              v
    +-------------------+
    | Clinical           |  Map niches to clinical
    | Interpretation     |  significance
    +-------------------+
```

### 5.3 Spatial Data Schema

The `sc_spatial` collection stores spatial niche data with the following key fields:
- `niche_label`: Descriptive niche name (e.g., "Tumor-immune interface")
- `platform`: Spatial technology (Visium, MERFISH, Xenium, CODEX)
- `cell_types`: Pipe-delimited cell types in the niche
- `signature_genes`: Spatially variable genes characterizing the niche
- `morans_i`: Spatial autocorrelation statistic (0-1)
- `clinical_relevance`: Clinical significance text

---

## 6. RAG Search Architecture

### 6.1 Multi-Collection Parallel Search

```
User Query
    |
    v
BGE-small-en-v1.5 Embedding (384-dim)
    |
    v
+-- ThreadPoolExecutor (max_workers=12) --+
|                                          |
|  sc_cell_types    (w=0.14, k=50)        |
|  sc_markers       (w=0.12, k=40)        |
|  sc_spatial       (w=0.10, k=30)        |
|  sc_tme           (w=0.10, k=30)        |
|  sc_drug_response (w=0.09, k=20)        |
|  sc_literature    (w=0.08, k=20)        |
|  sc_methods       (w=0.07, k=15)        |
|  sc_datasets      (w=0.06, k=15)        |
|  sc_trajectories  (w=0.07, k=20)        |
|  sc_pathways      (w=0.07, k=20)        |
|  sc_clinical      (w=0.07, k=15)        |
|  genomic_evidence (w=0.03, k=20)        |
|                                          |
+------------------------------------------+
    |
    v
Score Aggregation
    |-- Weighted score = cosine_similarity * collection_weight
    |-- Deduplication across collections
    |-- Score threshold filter (>= 0.4)
    |
    v
Evidence Ranking
    |-- Sort by weighted score
    |-- Citation relevance: HIGH (>0.75), MEDIUM (>0.60), LOW
    |
    v
Context Window Construction
    |-- Top-K evidence formatted for LLM
    |-- Conversation history (3-turn window)
    |
    v
Claude Sonnet Synthesis
    |
    v
SCResponse
```

### 6.2 Workflow-Specific Weight Boosting

When a query is classified as a specific workflow type, the default weights are replaced with a workflow-optimized profile. For example, a TME Profiling query boosts `sc_tme` from 0.10 to 0.25 and redistributes weight from less relevant collections.

### 6.3 Query Expansion

The `query_expansion.py` module expands queries with:
- Cell type synonyms (e.g., "T cell" -> "T lymphocyte", "CD3+ cell")
- Gene aliases (e.g., "PD-L1" -> "CD274", "B7-H1")
- Disease synonyms (e.g., "lung cancer" -> "NSCLC", "non-small cell lung carcinoma")
- 232 cell type aliases mapped to canonical names

---

## 7. Cross-Agent Communication

### 7.1 Integration Architecture

```
+------------------+      +------------------+
| Single-Cell      | <--> | Genomics Agent   |
| Intelligence     |      | (port 8527)      |
| Agent            |      +------------------+
| (port 8540)      |
|                  | <--> +------------------+
|                  |      | Biomarker Agent  |
|                  |      | (port 8529)      |
|                  |      +------------------+
|                  |
|                  | <--> +------------------+
|                  |      | Oncology Agent   |
|                  |      | (port 8528)      |
|                  |      +------------------+
|                  |
|                  | <--> +------------------+
|                  |      | Trial Agent      |
|                  |      | (port 8538)      |
+------------------+      +------------------+
         |
         v
+------------------+
| Shared Milvus    |
| genomic_evidence |
| (3.56M records)  |
+------------------+
```

### 7.2 Communication Protocol

- **Protocol:** HTTP REST
- **Timeout:** 30 seconds per agent
- **Failure mode:** Graceful degradation (response returned without cross-agent data)
- **Authentication:** Internal network trust (no auth between agents)

---

## 8. Data Storage Architecture

### 8.1 Milvus Vector Database

- **Index type:** IVF_FLAT
- **Metric:** COSINE similarity
- **nlist:** 128
- **Embedding dimension:** 384
- **Total estimated records:** 3,765,000 across 12 collections

### 8.2 Conversation Store

- **Backend:** Disk-backed JSON files
- **Location:** `data/cache/conversations/{session_id}.json`
- **TTL:** 24 hours
- **Format:** `{session_id, updated, messages[]}`

### 8.3 Knowledge Base

- **Backend:** In-memory Python dictionaries (loaded from `knowledge.py`)
- **Size:** 57 cell types, 30 drugs, 75 markers, 10 signatures, 25 L-R pairs, 12 TME profiles
- **Update:** Code deployment (static knowledge), scheduled ingest (dynamic data)

---

## 9. Security Architecture

```
External Client
       |
       v
TLS Termination (nginx/traefik)
       |
       v
API Key Authentication (X-API-Key header)
       |
       v
Rate Limiting (100 req/min per IP)
       |
       v
Request Size Limiting (10 MB)
       |
       v
CORS Validation (configured origins)
       |
       v
Pydantic Input Validation
       |
       v
Application Logic
       |
       v
Non-root Container (scuser)
```

---

## 10. Deployment Architecture

### 10.1 Standalone Deployment

```
docker-compose.yml
|
+-- milvus-etcd (etcd:v3.5.5)
+-- milvus-minio (minio)
+-- milvus-standalone (milvus:v2.4)
+-- sc-api (FastAPI, port 8540)
+-- sc-streamlit (Streamlit, port 8130)
+-- sc-setup (one-shot seed)
```

### 10.2 Integrated DGX Spark Deployment

```
docker-compose.dgx-spark.yml (top-level)
|
+-- Shared Milvus (port 19530)
+-- ... (other agents)
+-- sc-api (port 8540)
+-- sc-streamlit (port 8130)
```

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Architecture Guide*

---

!!! warning "Clinical Decision Support Disclaimer"
    The Single-Cell Intelligence Agent is a clinical decision support research tool for single-cell genomic analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
