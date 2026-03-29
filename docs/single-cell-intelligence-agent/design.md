# Single-Cell Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## 1. Executive Summary

The Single-Cell Intelligence Agent extends the HCLS AI Factory platform to deliver RAG-powered clinical decision support at single-cell resolution. It bridges the gap between high-dimensional transcriptomic data and bedside treatment decisions by combining curated domain knowledge, vector-based evidence retrieval, and LLM-powered synthesis. The agent transforms raw single-cell analysis outputs into clinically actionable insights for oncology tumor boards, immunology teams, and cell therapy programs.

The platform enables cross-functional queries like *"What is the TME composition of this melanoma sample and which immunotherapy is most likely to respond?"* that simultaneously search cell type references, TME profiles, drug sensitivity databases, spatial transcriptomics data, and clinical evidence — returning grounded answers with citations.

### Key Results

| Metric | Value |
|---|---|
| Total vectors indexed | **~850** across 12 Milvus collections (11 owned + 1 read-only) |
| Cell types in knowledge base | **57** with canonical marker genes |
| Analysis workflows | **10** clinical workflows |
| Decision support engines | **4** deterministic engines (TME classifier, drug response, CAR-T validator, spatial deconvolution) |
| Drugs modeled | **30** with cell-type-specific sensitivity profiles |
| Marker genes cataloged | **75** validated associations |
| Immune signatures | **10** TME-specific signatures |
| Cancer TME atlas profiles | **12** tumor-type reference profiles |
| Ligand-receptor pairs | **25** curated communication axes |
| API endpoints | **25** |
| Source code lines | **14,560** |
| Test cases | **~185** |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Single-Cell Agent Role |
|---|---|
| **DataStore** | Raw files: CellxGene datasets, marker gene references, TME profiles, spatial data JSON |
| **DataEngine** | Event-driven ingest pipelines (fetch → parse → embed → store) with CellxGene, marker, and TME parsers |
| **DataBase** | 12 Milvus collections (11 owned + 1 read-only) + knowledge base (57 cell types, 30 drugs, 75 markers) |
| **InsightEngine** | BGE-small embedding + multi-collection RAG + query expansion + 4 decision engines |
| **AgentEngine** | SingleCellAgent (Plan-Search-Evaluate-Synthesize-Report) + Streamlit UI |

### 2.2 System Diagram

```
                        ┌─────────────────────────────────┐
                        │    Streamlit Chat UI (8130)       │
                        │    Cell annotation | TME |        │
                        │    Drug response | Spatial        │
                        └──────────────┬──────────────────┘
                                       │
                        ┌──────────────▼──────────────────┐
                        │     FastAPI REST API (8540)       │
                        │     25 endpoints, CORS, Auth      │
                        │     Rate limiting, Metrics         │
                        └──────────────┬──────────────────┘
                                       │
                        ┌──────────────▼──────────────────┐
                        │     SingleCellAgent               │
                        │  Plan → Search → Evaluate →       │
                        │  Synthesize → Report               │
                        └──────────────┬──────────────────┘
                                       │
            ┌──────────────────────────┼───────────────────────────┐
            │                          │                           │
  ┌─────────▼──────────┐   ┌──────────▼──────────┐   ┌──────────▼──────────┐
  │ Decision Engines    │   │ RAG Pipeline         │   │ Workflow Engine      │
  │                     │   │                      │   │                      │
  │ TME Classifier      │   │ BGE-small-en-v1.5    │   │ 10 clinical          │
  │ Drug Response       │   │ (384-dim embedding)  │   │ workflows            │
  │ CAR-T Validator     │   │         │            │   │                      │
  │ Spatial Deconv.     │   │         ▼            │   │ Cell annotation      │
  │                     │   │ Parallel Search      │   │ TME profiling        │
  │                     │   │ 12 Milvus Collections│   │ Drug response        │
  │                     │   │ (ThreadPoolExecutor) │   │ Subclonal arch.      │
  │                     │   │         │            │   │ Spatial niche        │
  │                     │   │         ▼            │   │ Trajectory           │
  │                     │   │ Claude Sonnet 4.6    │   │ Ligand-receptor      │
  └─────────────────────┘   └──────────────────────┘   │ Biomarker discovery  │
            │                          │                │ CAR-T validation     │
            │                          │                │ Treatment monitor    │
  ┌─────────▼──────────────────────────▼───────────┐   └──────────────────────┘
  │                  Milvus 2.4 — 12 Collections    │
  │                                                  │
  │  sc_cell_types (57)        sc_marker_genes (75)  │
  │  sc_tme_profiles (12)      sc_drug_response (30) │
  │  sc_spatial_data (40)      sc_ligand_receptor(25)│
  │  sc_clinical_evidence(200) sc_immune_sigs (10)   │
  │  sc_trajectories (35)      sc_trials (80)        │
  │  sc_methods (50)           genomic_evidence [RO]  │
  └──────────────────────────────────────────────────┘
```

---

## 3. Data Collections — Actual State

All 12 collections (11 owned + 1 read-only) are populated and searchable.

### 3.1 `sc_cell_types` — 57 records

| Attribute | Value |
|---|---|
| **Source** | CellxGene Census, PanglaoDB, curated ontology |
| **Fields** | cell_type_id, name, lineage, canonical_markers, tissue_distribution, disease_associations, frequency_range |
| **Embedding** | FLOAT_VECTOR(384), BGE-small-en-v1.5 |
| **Index** | IVF_FLAT, COSINE, nlist=128, nprobe=16 |
| **Coverage** | T cells (12 subtypes), B cells (5), myeloid (10), stromal (8), epithelial (6), endothelial (4), neural (3), stem/progenitor (5), other (4) |

### 3.2 `sc_marker_genes` — 75 records

| Attribute | Value |
|---|---|
| **Source** | CellMarker 2.0, PanglaoDB, curated literature |
| **Fields** | gene_symbol, cell_type, specificity_score, expression_pattern, validation_status, tissue_context |
| **Coverage** | Canonical markers (CD3D, CD8A, CD19, CD14, EPCAM, etc.) plus emerging markers from scRNA-seq atlases |

### 3.3 `sc_tme_profiles` — 12 records

| Attribute | Value |
|---|---|
| **Source** | TCGA deconvolution studies, single-cell tumor atlases |
| **Tumor Types** | Melanoma, NSCLC, breast (TNBC/HR+), colorectal, pancreatic, glioblastoma, renal, ovarian, head/neck, hepatocellular, bladder, prostate |
| **Fields** | tumor_type, tme_class (hot/cold/excluded/immunosuppressive), immune_infiltrate_composition, response_predictors, resistance_mechanisms |

### 3.4 `sc_drug_response` — 30 records

| Attribute | Value |
|---|---|
| **Source** | GDSC, CCLE, curated clinical correlative studies |
| **Fields** | drug_name, mechanism, cell_type_sensitivity, resistance_markers, combination_synergies, clinical_evidence |
| **Drug Classes** | Checkpoint inhibitors (6), targeted therapies (8), chemotherapies (6), ADCs (4), CAR-T-related (3), novel agents (3) |

### 3.5 `sc_spatial_data` — 40 records

| Attribute | Value |
|---|---|
| **Source** | Curated from Visium, MERFISH, CosMx, Slide-seq publications |
| **Fields** | platform, tissue, spatial_resolution, niche_type, cell_interactions, clinical_relevance |
| **Platforms** | 10x Visium, MERFISH, CosMx SMI, Slide-seq V2, CODEX, IMC |

### 3.6 `sc_ligand_receptor` — 25 records

| Attribute | Value |
|---|---|
| **Source** | CellPhoneDB, CellChat, NicheNet |
| **Fields** | ligand, receptor, cell_type_source, cell_type_target, pathway, functional_role, therapeutic_target |

### 3.7 Index Configuration (all collections)

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist | 128 (all collections) |
| nprobe | 16 |
| Embedding dim | 384 (BGE-small-en-v1.5) |

---

## 4. Knowledge Base

### 4.1 Cell Type Taxonomy (57 entries)

Each entry includes: canonical markers, lineage hierarchy, tissue distribution, disease associations, functional states, and therapeutic relevance.

| Lineage | Count | Key Types |
|---|---|---|
| **T cells** | 12 | CD4+ naive, Th1, Th2, Th17, Treg, CD8+ naive, CD8+ effector, CD8+ exhausted, gamma-delta, NKT, MAIT, CD4+ Tfh |
| **B cells** | 5 | Naive B, memory B, germinal center B, plasmablast, plasma cell |
| **Myeloid** | 10 | Classical monocyte, non-classical monocyte, cDC1, cDC2, pDC, M1 macrophage, M2 macrophage, MDSC, mast cell, neutrophil |
| **Stromal** | 8 | Fibroblast, myofibroblast, CAF, pericyte, smooth muscle, mesenchymal stem, adipocyte, osteoblast |
| **Epithelial** | 6 | Basal, luminal, secretory, ciliated, alveolar type I, alveolar type II |
| **Endothelial** | 4 | Arterial, venous, lymphatic, tip cell |
| **Other** | 12 | Neural crest, melanocyte, hepatocyte, podocyte, Schwann cell, astrocyte, oligodendrocyte, microglia, erythrocyte precursor, megakaryocyte, stem cell, progenitor |

### 4.2 TME Classification Framework

| TME Class | Definition | Response Prediction |
|---|---|---|
| **Immune-hot** | >30% immune infiltrate, CD8+ T cell enriched, PD-L1+ | High checkpoint inhibitor response (40-60%) |
| **Immune-cold** | <5% immune infiltrate, minimal T cell presence | Poor immunotherapy response (<10%) |
| **Immune-excluded** | Immune cells at tumor margin but not penetrating | Moderate response, requires combination therapy |
| **Immunosuppressive** | Treg/MDSC-enriched, TGF-beta high | Requires TME reprogramming before immunotherapy |

### 4.3 CAR-T Target Validation Criteria

The decision engine evaluates proposed CAR-T targets against five safety and efficacy dimensions:

1. **Target expression** — Must be expressed on >80% of tumor cells (scRNA-seq quantification)
2. **Normal tissue expression** — Flags targets with expression on essential normal tissues
3. **Known escape mechanisms** — Cross-references antigen loss data from CAR-T clinical trials
4. **Immunosuppressive TME** — Warns if TME classification predicts poor T cell function
5. **Spatial accessibility** — Assesses target accessibility from spatial transcriptomics data

---

## 5. Multi-Collection RAG Engine

### 5.1 Search Flow

```
User Query: "What cell types dominate the TME in this melanoma sample?"
    │
    ├── 1. Embed query with BGE asymmetric prefix               [< 5 ms]
    │      "Represent this sentence for searching relevant passages: ..."
    │
    ├── 2. Parallel search across 12 collections (top-5 each)   [10-16 ms]
    │   ├── sc_cell_types:        Melanoma-associated cell types  (score: 0.80-0.88)
    │   ├── sc_tme_profiles:      Melanoma TME reference          (score: 0.82-0.90)
    │   ├── sc_marker_genes:      Melanoma markers                (score: 0.75-0.83)
    │   ├── sc_drug_response:     Melanoma drug sensitivities     (score: 0.72-0.80)
    │   └── sc_spatial_data:      Melanoma spatial niche data     (score: 0.70-0.78)
    │
    ├── 3. Query expansion: "melanoma TME" → [melanoma, BRAF,    [< 1 ms]
    │      immune checkpoint, PD-L1, T cell, macrophage, ...]
    │
    ├── 4. Weighted merge + deduplicate (cap at 30 results)      [< 1 ms]
    │      Weights: cell_types 0.15, tme 0.15, markers 0.12,
    │               drug 0.10, spatial 0.10, evidence 0.10, ...
    │
    ├── 5. Knowledge base augmentation                           [< 1 ms]
    │      Melanoma TME → immune-hot/cold classification
    │      + canonical markers + drug sensitivities
    │
    └── 6. Stream Claude Sonnet 4.6 response                    [~22-24 sec]
           Grounded answer with citations and cell type
           composition breakdown
```

**Total: ~24 sec** (dominated by LLM generation; retrieval is ~20 ms)

### 5.2 Collection Weights

| Collection | Weight | Rationale |
|---|---|---|
| sc_cell_types | 0.15 | Cell type definitions are the primary reference |
| sc_tme_profiles | 0.15 | TME classification drives treatment decisions |
| sc_marker_genes | 0.12 | Marker validation supports cell type annotation |
| sc_drug_response | 0.10 | Drug sensitivity is the clinical endpoint |
| sc_spatial_data | 0.10 | Spatial context adds tissue-level intelligence |
| sc_clinical_evidence | 0.10 | Published evidence grounds recommendations |
| sc_ligand_receptor | 0.07 | Cell communication explains TME dynamics |
| sc_immune_sigs | 0.06 | Immune signatures refine TME classification |
| sc_trajectories | 0.05 | Differentiation trajectories inform prognosis |
| sc_trials | 0.04 | Active trials provide translational context |
| sc_methods | 0.03 | Platform selection guides experimental design |
| genomic_evidence | 0.03 | Shared genomic context for cross-agent queries |
| **Total** | **1.00** | |

### 5.3 Embedding Strategy

BGE-small-en-v1.5 uses asymmetric encoding:

| Mode | Prefix | Usage |
|---|---|---|
| **Query** | `"Represent this sentence for searching relevant passages: "` | User questions via `_embed_query()` |
| **Document** | None (raw text) | Ingested records via `to_embedding_text()` |

---

## 6. Decision Support Engines

### 6.1 TME Classifier

Classifies tumor microenvironment into four categories (hot/cold/excluded/immunosuppressive) based on cell type proportions, spatial distribution, and immune signature scores.

Input: Cell type composition (dict of cell_type → fraction)
Output: TME classification with confidence score, treatment recommendations, and supporting evidence.

### 6.2 Drug Response Predictor

Maps cell-type-specific expression patterns to drug sensitivity profiles across 30 modeled drugs. Integrates resistance markers and combination synergy data.

### 6.3 CAR-T Target Validator

Evaluates proposed CAR-T targets against the five validation criteria (expression, normal tissue, escape, TME, spatial) and returns a go/no-go recommendation with risk factors.

### 6.4 Spatial Deconvolution Interpreter

Interprets spatial transcriptomics deconvolution results to identify tissue niches, immune exclusion zones, and cell-cell interaction hotspots.

---

## 7. Clinical Workflows

| # | Workflow | Clinical Question |
|---|---|---|
| 1 | Cell Type Annotation | "What cell types are in this sample and at what proportions?" |
| 2 | TME Profiling | "Is this tumor hot, cold, excluded, or immunosuppressive?" |
| 3 | Drug Response | "Which drugs will this tumor respond to at the cellular level?" |
| 4 | Subclonal Architecture | "Are there resistant subclones that could cause relapse?" |
| 5 | Spatial Niche | "Where are the immune cells relative to tumor cells in tissue?" |
| 6 | Trajectory Analysis | "What differentiation or exhaustion trajectories are active?" |
| 7 | Ligand-Receptor | "Which cell-cell communication axes are driving tumor progression?" |
| 8 | Biomarker Discovery | "What cell-type-specific biomarkers predict treatment outcome?" |
| 9 | CAR-T Validation | "Is this target safe and effective for CAR-T therapy?" |
| 10 | Treatment Monitoring | "How has the TME changed between treatment timepoints?" |

---

## 8. Performance Benchmarks

Measured on NVIDIA DGX Spark (GB10 GPU, 128GB unified LPDDR5x memory, 20 ARM cores).

### 8.1 Search Performance

| Operation | Latency | Notes |
|---|---|---|
| Single collection search (top-5) | 3-5 ms | Milvus IVF_FLAT with cached index |
| 12-collection parallel search (top-5 each) | 10-16 ms | ThreadPoolExecutor, 60 total results |
| Query expansion + filtered search | 6-10 ms | Up to 5 expanded terms |
| Full retrieve() pipeline | 18-28 ms | Embed + search + expand + merge + knowledge |

### 8.2 RAG Query Performance

| Operation | Latency | Notes |
|---|---|---|
| Full query (retrieve + Claude generate) | ~24 sec | Dominated by LLM generation |
| Streaming query (time to first token) | ~3 sec | Evidence returned immediately |
| Response length | 800-2000 chars | Grounded answer with citations |

### 8.3 Decision Engine Performance

| Engine | Latency | Notes |
|---|---|---|
| TME Classifier | <50 ms | Dictionary-based classification |
| Drug Response Predictor | <100 ms | 30-drug scan |
| CAR-T Target Validator | <100 ms | 5-criteria evaluation |
| Spatial Deconvolution | <50 ms | Niche identification |
| **All 4 engines combined** | **<300 ms** | |

---

## 9. Infrastructure

### 9.1 Technology Stack

| Component | Technology |
|---|---|
| Language | Python 3.10+ |
| Vector DB | Milvus 2.4, localhost:19530 |
| Embeddings | BGE-small-en-v1.5 (BAAI) — 384-dim |
| LLM | Claude Sonnet 4.6 (Anthropic API) |
| Web UI | Streamlit (port 8130, NVIDIA black/green theme) |
| REST API | FastAPI + Uvicorn (port 8540) |
| Configuration | Pydantic BaseSettings |
| Data models | Pydantic BaseModel + Field validation |
| Hardware target | NVIDIA DGX Spark (GB10 GPU, 128GB unified, $4,699) |

### 9.2 Service Ports

| Port | Service |
|---|---|
| 8540 | FastAPI REST API |
| 8130 | Streamlit Chat UI |
| 19530 | Milvus vector database (shared) |

### 9.3 Dependencies on HCLS AI Factory

| Dependency | Usage |
|---|---|
| Milvus 2.4 instance | Shared vector database — Single-Cell adds 11 owned collections alongside existing `genomic_evidence` (read-only) |
| `ANTHROPIC_API_KEY` | Loaded from `rag-chat-pipeline/.env` if not set in environment |
| BGE-small-en-v1.5 | Same embedding model as main RAG pipeline |

---

## 10. Demo Scenarios

### 10.1 Validated Demo Queries

**1. "Annotate the cell types in this NSCLC sample"**
- Searches: cell_types (lung-associated), marker_genes (EPCAM, CD45, CD3D), tme_profiles (NSCLC reference)
- Knowledge base: NSCLC TME → immune-excluded pattern, EGFR mutation context
- Expected: Epithelial (40%), T cells (15%), macrophages (20%), fibroblasts (15%), endothelial (10%)

**2. "Is this melanoma immune-hot or immune-cold?"**
- Searches: tme_profiles (melanoma), immune_sigs (checkpoint signatures), cell_types (immune subtypes)
- Expected: TME classification with CD8/Treg ratio, PD-L1 status, spatial distribution context

**3. "Validate CD19 as a CAR-T target for this B-ALL sample"**
- Searches: cell_types (B lineage), drug_response (CD19 CAR-T), spatial_data (bone marrow niche)
- CAR-T Validator: Expression check, normal tissue risk, escape mechanisms, TME assessment
- Expected: Go recommendation with CD19 loss risk warning (28% of relapses)

**4. "Which checkpoint inhibitor would work best for this triple-negative breast cancer?"**
- Searches: drug_response (TNBC sensitivities), tme_profiles (TNBC), immune_sigs (PD-L1, TMB)
- Expected: Pembrolizumab for CPS >= 10, atezolizumab + nab-paclitaxel for PD-L1+ IC

**5. "Map the spatial niches in this colorectal cancer Visium dataset"**
- Searches: spatial_data (CRC Visium), tme_profiles (CRC), ligand_receptor (tumor-stroma communication)
- Expected: Tumor core, immune margin, desmoplastic stroma, and tertiary lymphoid structure niches

---

## 11. File Structure (Actual)

```
single_cell_intelligence_agent/
├── src/
│   ├── __init__.py
│   ├── agent.py                     # Plan-Search-Evaluate-Synthesize-Report (2,090 lines)
│   ├── models.py                    # Pydantic data models (820 lines)
│   ├── collections.py               # 12 Milvus collection schemas (1,210 lines)
│   ├── rag_engine.py                # Multi-collection RAG search (1,490 lines)
│   ├── clinical_workflows.py        # 10 analysis workflows (1,792 lines)
│   ├── decision_support.py          # 4 clinical engines (886 lines)
│   ├── knowledge.py                 # Domain knowledge base (1,816 lines)
│   ├── query_expansion.py           # Synonym expansion (893 lines)
│   ├── cross_modal.py               # Inter-agent communication (392 lines)
│   ├── metrics.py                   # Prometheus metrics (476 lines)
│   ├── export.py                    # Report generation (588 lines)
│   ├── scheduler.py                 # APScheduler ingest (496 lines)
│   └── ingest/
│       ├── base.py                  # BaseIngestParser ABC (228 lines)
│       ├── cellxgene_parser.py      # CellxGene data parser (679 lines)
│       ├── marker_parser.py         # Marker gene parser (286 lines)
│       └── tme_parser.py            # TME profile parser (418 lines)
├── app/
│   └── single_cell_ui.py           # Streamlit chat UI
├── api/
│   └── main.py                     # FastAPI REST server (25 endpoints)
├── config/
│   └── settings.py                 # Pydantic BaseSettings (197 lines)
├── data/
│   ├── reference/                  # Seed data JSON files
│   └── cache/                      # Embedding cache
├── scripts/
│   ├── setup_collections.py        # Create Milvus schemas
│   ├── seed_knowledge.py           # Populate knowledge base
│   └── validate_e2e.py             # End-to-end validation
├── tests/                          # ~185 test cases
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
└── README.md
```

**~40 Python files | ~14,560 lines of code | Apache 2.0**

---

## 12. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | Data models, 12 collection schemas, knowledge base, 4 decision engines, RAG engine, agent orchestrator |
| **Phase 2: Data** | Complete | 57 cell types, 75 marker genes, 12 TME profiles, 30 drug records, 25 L-R pairs, spatial data |
| **Phase 3: RAG Integration** | Complete | Multi-collection parallel search, knowledge augmentation, Claude Sonnet 4.6 streaming |
| **Phase 4: Workflows** | Complete | 10 clinical workflows with workflow-specific collection weight boosting |
| **Phase 5: UI + Demo** | Complete | Streamlit UI on port 8130, NVIDIA theme, demo query buttons, streaming responses |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| RAPIDS GPU acceleration for deconvolution | Medium | 2-3 days |
| scGPT/Geneformer foundation model integration | Low | 1-2 weeks |
| Integration with HCLS AI Factory landing page | Low | 1 hour |

---

## 13. Relationship to HCLS AI Factory

The Single-Cell Intelligence Agent demonstrates the **resolution extension** of the HCLS AI Factory architecture. While the core platform operates at the variant level (VCF → drug candidates), this agent extends analysis to single-cell resolution — the frontier of precision medicine.

- **Same Milvus instance** — 11 new owned collections alongside existing `genomic_evidence` (3.56M vectors, read-only)
- **Same embedding model** — BGE-small-en-v1.5 (384-dim)
- **Same LLM** — Claude via Anthropic API
- **Same hardware** — NVIDIA DGX Spark ($4,699)
- **Same patterns** — Pydantic models, BaseIngestPipeline, knowledge graph, query expansion

The key insight: bulk genomics identifies *what* variants are present; single-cell analysis reveals *which cells* are affected and *how* the tumor microenvironment responds. This agent closes the resolution gap between genome-level and cell-level precision medicine.

---

## 14. Credits

- **Adam Jones**
- **Apache 2.0 License**

---

!!! warning "Clinical Decision Support Disclaimer"
    The Single-Cell Intelligence Agent is a clinical decision support research tool for single-cell transcriptomic analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
