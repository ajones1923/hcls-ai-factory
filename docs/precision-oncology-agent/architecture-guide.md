# Precision Oncology Intelligence Agent -- Architecture Guide

Architecture design document for the Precision Oncology Intelligence Agent,
part of the HCLS AI Factory on NVIDIA DGX Spark.

Author: Adam Jones
Date: March 2026

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [VAST AI OS Component Mapping](#vast-ai-os-component-mapping)
3. [System Architecture](#system-architecture)
4. [Component Architecture](#component-architecture)
5. [Data Flow](#data-flow)
6. [Milvus Collections](#milvus-collections)
7. [Embedding Strategy](#embedding-strategy)
8. [LLM Integration](#llm-integration)
9. [Knowledge Graph](#knowledge-graph)
10. [Clinical Pipelines](#clinical-pipelines)
11. [RAG Engine](#rag-engine)
12. [API Layer](#api-layer)
13. [UI Layer](#ui-layer)
14. [Export Pipeline](#export-pipeline)
15. [Scaling and Performance](#scaling-and-performance)
16. [Security](#security)
17. [File Structure](#file-structure)

---

## Executive Summary

The Precision Oncology Intelligence Agent is a RAG-powered clinical
decision-support system purpose-built for molecular tumor board (MTB)
workflows. It operates as Stage 2.5 of the HCLS AI Factory pipeline,
sitting between the genomics pipeline (Stage 1) and the drug discovery
pipeline (Stage 3) to provide real-time variant interpretation, therapy
ranking, clinical trial matching, and evidence synthesis.

The agent combines 11 Milvus vector collections, domain-aware query
expansion, a multi-step plan-search-evaluate-synthesize agent loop, and
Claude LLM synthesis to deliver guideline-concordant therapy recommendations
with full evidence provenance.

### Key Results

| Metric | Value |
|--------|-------|
| Python files | 66 |
| Total lines of code | ~20,490 |
| Milvus collections | 11 (10 owned + 1 read-only shared) |
| Live vectors | 609 (526 seed + 83 knowledge graph) |
| Actionable gene targets | 40+ |
| Therapy mappings | 30+ |
| Resistance mechanisms | 12+ |
| Oncogenic pathways | 10+ |
| Biomarker panels | 20+ |
| Test files / cases | 10 files, 556 test cases, all passing |
| Docker services | 6 |
| Seed data files | 10 JSON files, ~773 KB total |
| Export formats | 4 (Markdown, JSON, PDF, FHIR R4) |
| Cancer types supported | 26 |
| Enumerations | 13 |
| Domain models | 10 |

---

## VAST AI OS Component Mapping

The Oncology Agent maps to the VAST AI OS architecture as follows:

| VAST AI OS Layer | Agent Component |
|-----------------|-----------------|
| Data Layer | Milvus vector store (11 collections), 10 JSON seed files |
| Model Layer | BAAI/bge-small-en-v1.5 (embedding), Claude claude-sonnet-4-6 (synthesis) |
| Inference Layer | FastAPI server (8527), RAG engine, therapy ranker, trial matcher |
| Application Layer | Streamlit UI (8526), multi-tab MTB interface |
| Orchestration | Docker Compose (6 services), health checks, APScheduler |
| Integration | Cross-agent stubs, genomic_evidence shared collection, event bus |

---

## System Architecture

```
+-------------------------------------------------------------------+
|                    Streamlit UI (:8526)                            |
|  +-------+ +--------+ +--------+ +--------+ +--------+ +------+  |
|  |Clinical| |Case    | |Variant | |Therapy | |Trial   | |MTB   |  |
|  | Query  | |Manager | |Viewer  | |Ranking | |Match   | |Packet|  |
|  +-------+ +--------+ +--------+ +--------+ +--------+ +------+  |
|  +-------+ +--------+ +--------+                                  |
|  |Pathway | |Resist  | |Knwldge |                                 |
|  |Explorer| |Analys  | |Base    |                                 |
|  +-------+ +--------+ +--------+                                  |
+-------------------------------------------------------------------+
        |                                         |
        v                                         v
+-------------------+                 +----------------------------+
|  FastAPI (:8527)  |                 |  OncoIntelligenceAgent     |
|  REST endpoints   |                 |  plan()                    |
|  Auth / CORS      |                 |  search()                  |
|  Prometheus /     |                 |  evaluate()                |
|    metrics        |                 |  synthesize()              |
+-------------------+                 +----------------------------+
        |                                     |
        v                                     v
+-------------------+            +----------------------------+
| OncoRAGEngine     |            | Clinical Engines           |
|  cross_collection |            |                            |
|    _search()      |            | TherapyRanker              |
|  query()          |            |   7-step ranking algorithm  |
|  compare()        |            |                            |
|  _embed()         |            | TrialMatcher               |
|  _build_prompt()  |            |   hybrid deterministic +   |
+-------------------+            |   semantic matching        |
        |                        |                            |
        v                        | OncologyCaseManager        |
+-------------------+            |   VCF parsing, case CRUD,  |
| CollectionManager |            |   MTB packet generation    |
|  connect()        |            |                            |
|  search_all()     |            | QueryExpander              |
|  insert_batch()   |            |   12-category domain       |
+-------------------+            |   expansion                |
        |                        |                            |
        v                        | CrossModalIntegrator       |
+-------------------+            |   genomic-imaging-drug     |
| Milvus (:19530)   |            |   discovery events         |
| 11 collections    |            |                            |
| IVF_FLAT / COSINE |            | MetricsCollector           |
| 384-dim vectors   |            |   Prometheus + custom      |
+-------------------+            +----------------------------+
        |                                |
        v                                v
+-------------------+            +----------------------------+
| BGE-small-en-v1.5 |            | OncologyExporter           |
| Embedding Model   |            |   Markdown / JSON          |
| 384 dimensions    |            |   PDF / FHIR R4            |
+-------------------+            +----------------------------+
```

---

## Component Architecture

The agent is organized into four layers with clear separation of concerns.

### Presentation Layer

- **Streamlit UI** (`app/`): Multi-tab application for MTB workflows covering
  clinical query, case management, variant viewing, therapy ranking, trial
  matching, MTB packet generation, pathway exploration, resistance analysis,
  and knowledge base management.

### API Layer

- **FastAPI Server** (`api/`): RESTful API with endpoints for querying,
  case management, therapy ranking, trial matching, export, and health checks.
  Includes Prometheus metrics and APScheduler for background tasks.

### Engine Layer

Six specialized engines orchestrated by the OncoIntelligenceAgent:

| Engine | Module | Lines | Purpose |
|--------|--------|-------|---------|
| OncoIntelligenceAgent | `agent.py` | 553 | Plan-search-evaluate-synthesize loop |
| OncoRAGEngine | `rag_engine.py` | 908 | Multi-collection RAG with query expansion |
| TherapyRanker | `therapy_ranker.py` | 748 | 7-step evidence-based therapy ranking |
| TrialMatcher | `trial_matcher.py` | 513 | Hybrid deterministic + semantic trial matching |
| OncologyCaseManager | `case_manager.py` | 516 | Case CRUD, VCF parsing, MTB packets |
| QueryExpander | `query_expansion.py` | 812 | 12-category domain-aware query expansion |
| CrossModalIntegrator | `cross_modal.py` | -- | Cross-agent event handling |
| MetricsCollector | `metrics.py` | -- | Prometheus metrics and monitoring |

### Data Layer

- **Milvus** (11 vector collections, 384-dim BGE-small-en-v1.5 embeddings)
- **Knowledge Graph** (`knowledge.py`, 1,662 lines): 5 knowledge domains
- **Seed Data** (10 JSON files in `data/`)
- **Models** (`models.py`, 538 lines): 13 enums, 10 domain models, 4 search models, 2 agent I/O models

---

## Data Flow

### End-to-End Pipeline

```
                    HCLS AI Factory Pipeline
                    =======================

Stage 1: Genomics Pipeline (Parabricks/DeepVariant)
  FASTQ -> VCF -> genomic_evidence collection (shared, read-only)
      |
Stage 2.5: Oncology Intelligence Agent (this project)
      |
      +-> Variant Interpretation (CIViC/OncoKB evidence)
      +-> Therapy Ranking (NCCN/ESMO guideline-concordant)
      +-> Clinical Trial Matching (ClinicalTrials.gov)
      +-> MTB Packet Generation (Markdown/JSON/PDF/FHIR)
      |
Stage 3: Drug Discovery Pipeline (BioNeMo/DiffDock/RDKit)
  Therapy targets -> lead compound optimization
```

### Agent Internal Pipeline

The OncoIntelligenceAgent executes a 4-step plan-search-evaluate-synthesize
loop for each query:

```
User Question
    |
    v
[1. PLAN] --> SearchPlan
    |   - Identify topics from keyword matching (20+ triggers)
    |   - Extract target genes (30 KNOWN_GENES)
    |   - Detect cancer types (25 canonical + 70+ aliases)
    |   - Select strategy: broad | targeted | comparative
    |   - Decompose complex queries into sub-questions
    |
    v
[2. SEARCH] --> Cross-collection retrieval
    |   - Primary query + all sub-questions
    |   - Optional query expansion (12 categories)
    |   - Parallel search across 11 collections
    |   - Weighted scoring with per-collection weights
    |
    v
[3. EVALUATE] --> "sufficient" | "partial" | "insufficient"
    |   - sufficient: >= 3 hits from >= 2 collections
    |   - If insufficient and retries remain: broaden and retry
    |   - MAX_RETRIES = 2
    |   - Minimum similarity score: 0.30
    |
    v
[4. SYNTHESIZE] --> AgentResponse
    |   - Knowledge injection (genes, therapies, resistance, pathways, biomarkers)
    |   - Claude claude-sonnet-4-6 generates answer with citations
    |   - Markdown report attached
    |
    v
AgentResponse (answer, evidence, knowledge_used, report)
```

### Case Creation Flow

```
Patient Data (ID, cancer type, stage, VCF, biomarkers, prior therapies)
    |
    v
[VCF Parsing] -- If raw VCF text, parse via cyvcf2
    |   Extract: gene, variant, chrom, pos, ref, alt, consequence
    |
    v
[Variant Annotation] -- Cross-reference ACTIONABLE_TARGETS
    |   Classify actionability (A/B/C/D/E/VUS)
    |
    v
[CaseSnapshot Creation] -- Generate UUID, build text_summary
    |   Embed via BGE-small-en-v1.5
    |
    v
[Persist to onco_cases] -- Insert into Milvus collection
    |
    v
[MTB Packet Generation] -- Combine variant table, evidence,
    |   therapy ranking, trial matches, open questions, citations
    |
    v
MTBPacket + CaseSnapshot returned
```

---

## Milvus Collections

The agent manages 11 specialized vector collections. All use COSINE similarity
with IVF_FLAT indexing (nlist=1024, nprobe=16) and 384-dimensional vectors
from BGE-small-en-v1.5.

| # | Collection Name | Description | Weight | Seed Records |
|---|----------------|-------------|--------|-------------|
| 1 | onco_variants | Actionable somatic/germline variants (CIViC/OncoKB) | 0.18 | 90 |
| 2 | onco_literature | PubMed/PMC/preprint literature chunks | 0.16 | 60 |
| 3 | onco_therapies | Approved and investigational therapies | 0.14 | 64 |
| 4 | onco_guidelines | NCCN/ASCO/ESMO guideline recommendations | 0.12 | 45 |
| 5 | onco_trials | ClinicalTrials.gov summaries with biomarker criteria | 0.10 | 55 |
| 6 | onco_biomarkers | Predictive and prognostic biomarkers | 0.08 | 50 |
| 7 | onco_resistance | Resistance mechanisms and bypass strategies | 0.07 | 50 |
| 8 | onco_pathways | Signaling pathways, cross-talk, druggable nodes | 0.06 | 35 |
| 9 | onco_outcomes | Real-world treatment outcome records | 0.04 | 40 |
| 10 | onco_cases | De-identified patient case snapshots | 0.02 | 37 |
| 11 | genomic_evidence | Shared VCF-derived genomic variants (read-only) | 0.03 | -- |
| | **Total** | | **1.00** | **526** |

### Collection Weight Distribution

```
onco_variants      ██████████████████████  0.18
onco_literature    ████████████████████    0.16
onco_therapies     █████████████████       0.14
onco_guidelines    ██████████████          0.12
onco_trials        ████████████            0.10
onco_biomarkers    █████████               0.08
onco_resistance    ████████                0.07
onco_pathways      ███████                 0.06
onco_outcomes      █████                   0.04
genomic_evidence   ███                     0.03
onco_cases         ██                      0.02
                                      Sum: 1.00
```

### Key Collection Schemas

**onco_variants** -- Actionable somatic/germline variants:

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| gene | VARCHAR(50) | Gene symbol |
| variant_name | VARCHAR(100) | Variant designation |
| variant_type | VARCHAR(30) | SNV, INDEL, CNV_AMP, FUSION, etc. |
| cancer_type | VARCHAR(50) | Associated cancer type |
| evidence_level | VARCHAR(20) | A (FDA) through E (Computational) |
| drugs | VARCHAR(500) | Indicated therapies |
| civic_id | VARCHAR(20) | CIViC database ID |
| vrs_id | VARCHAR(100) | GA4GH VRS identifier |
| text_summary | VARCHAR(3000) | Clinical narrative for embedding |
| clinical_significance | VARCHAR(200) | Pathogenic, likely pathogenic, VUS |
| allele_frequency | FLOAT | Population allele frequency |

**onco_therapies** -- Approved and investigational therapies:

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| drug_name | VARCHAR(200) | Generic drug name |
| category | VARCHAR(30) | TARGETED, IMMUNOTHERAPY, CHEMO, etc. |
| targets | VARCHAR(200) | Molecular targets |
| approved_indications | VARCHAR(500) | FDA-approved indications |
| resistance_mechanisms | VARCHAR(500) | Known resistance mechanisms |
| evidence_level | VARCHAR(20) | Evidence tier |
| text_summary | VARCHAR(3000) | Clinical summary for embedding |
| mechanism_of_action | VARCHAR(500) | MOA description |

---

## Embedding Strategy

| Parameter | Value |
|-----------|-------|
| Model | BAAI/bge-small-en-v1.5 |
| Parameters | 33M |
| Dimensions | 384 |
| Metric | COSINE |
| Index type | IVF_FLAT (nlist=1024, nprobe=16) |
| Batch size | 32 |
| Runtime | CPU (no GPU required) |
| Instruction prefix | "Represent this sentence for searching relevant passages: " |
| Search mode | Asymmetric (queries use instruction prefix, documents do not) |

### Domain-Optimized Embedding Text

Each of the 10 domain models provides a `to_embedding_text()` method that
constructs a domain-optimized text representation:

- **OncologyVariant:** Combines gene, variant_name, cancer_type, evidence_level,
  drugs, clinical_significance, and text_summary into a single embedding string.
- **OncologyTherapy:** Combines drug_name, category, targets, mechanism_of_action,
  approved_indications, and text_summary.
- **OncologyTrial:** Combines title, phase, cancer_types, biomarker_criteria,
  outcome_summary, and text_summary.

This ensures that structured metadata fields (gene names, drug names, trial IDs)
contribute to embedding similarity even when they appear only in metadata.

---

## LLM Integration

### Claude Configuration

| Setting | Value |
|---------|-------|
| Model | claude-sonnet-4-6 |
| Provider | Anthropic API |
| Environment variable | `ANTHROPIC_API_KEY` |
| Configuration prefix | `ONCO_` |
| Streaming | Supported |
| Conversation memory | Configurable |

### System Prompt Design

The system prompt defines the agent as an **Oncology Intelligence Agent**
with core competencies in:

1. Molecular profiling (TMB, MSI, CNV, fusions)
2. Variant interpretation (CIViC/OncoKB evidence levels, AMP/ASCO/CAP)
3. Therapy selection (NCCN/ESMO guideline-concordant)
4. Clinical trial matching (ClinicalTrials.gov, basket/umbrella)
5. Resistance mechanisms (on-target, bypass, lineage plasticity)
6. Biomarker assessment (TMB, MSI, PD-L1, HRD, companion diagnostics)
7. Outcomes monitoring (RECIST, survival, MRD, ctDNA)
8. Cross-modal integration (genomic-imaging-drug discovery)

Behavioral instructions enforce citation, cross-functional thinking,
resistance flagging, guideline referencing, and uncertainty acknowledgment.

### Comparative Retrieval

The engine detects comparative questions via regex
(`compare|vs|versus|difference between|head.to.head`) and routes them to
a dual-entity retrieval pipeline:

1. Parse entity A and entity B from the question
2. Retrieve evidence independently for each entity
3. Identify shared/head-to-head evidence (intersection by ID)
4. Build a structured comparison prompt with 8 comparison axes
5. Generate comparative synthesis via LLM

### Citation Formatting

- PubMed IDs: `[PubMed 12345](https://pubmed.ncbi.nlm.nih.gov/12345/)`
- NCT IDs: `[NCT01234567](https://clinicaltrials.gov/study/NCT01234567)`
- All others: `[Label: record_id]`

---

## Knowledge Graph

**Module:** `src/knowledge.py` (1,662 lines)

The knowledge graph provides curated, structured domain knowledge that is
injected into LLM prompts alongside retrieved evidence. It consists of
five primary domains.

### 1. ACTIONABLE_TARGETS (~40 genes)

Each entry contains: gene, full_name, cancer_types, key_variants,
targeted_therapies, combination_therapies, resistance_mutations, pathway,
evidence_level, and a free-text clinical description.

**Representative genes:** BRAF, EGFR, ALK, ROS1, KRAS, HER2, NTRK, RET,
MET, FGFR, PIK3CA, IDH1, IDH2, BRCA1, BRCA2, TP53, PTEN, CDKN2A, STK11,
ESR1, ERBB2, NRAS, APC, VHL, KIT, PDGFRA, FLT3, NPM1, DNMT3A, and others.

### 2. THERAPY_MAP (~30 drugs)

Maps drug names (lowercase) to structured metadata: brand_name, category,
drug_class, and guideline reference (NCCN/ESMO recommendation).

### 3. RESISTANCE_MAP (~12 mechanisms)

Maps drug names to documented resistance mechanisms: mutation (e.g., EGFR T790M
for erlotinib), next_line (recommended subsequent therapies), and mechanism
(description of resistance biology).

### 4. PATHWAY_MAP (~10 pathways)

Maps oncogenic signaling pathways to key_genes, therapeutic_targets, cross_talk,
and pathway descriptions. Covered pathways: MAPK, PI3K/AKT/mTOR, DDR,
Cell Cycle, Apoptosis, WNT, NOTCH, Hedgehog, JAK/STAT, Angiogenesis,
Hippo, NF-kB, TGF-beta.

### 5. BIOMARKER_PANELS (~20 panels)

Maps biomarker identifiers to clinical decision rules: marker name, threshold,
positive_values, recommended drugs, evidence_level, and guideline reference text.

### Helper Functions

- `lookup_gene(query)` -- Return knowledge context for gene mentions
- `lookup_therapy(query)` -- Return knowledge context for therapy mentions
- `lookup_resistance(query)` -- Return resistance mechanism context
- `lookup_pathway(query)` -- Return pathway context
- `lookup_biomarker(query)` -- Return biomarker context
- `get_target_context(gene)` -- Return full ACTIONABLE_TARGETS entry
- `classify_variant_actionability(gene, variant)` -- Return evidence tier

---

## Clinical Pipelines

### 1. Therapy Ranker (src/therapy_ranker.py, 748 lines)

Seven-step evidence-based therapy ranking algorithm:

```
Patient Profile (cancer_type, variants, biomarkers, prior_therapies)
    |
    v
[Step 1: Variant-Driven Therapies]
    |   ACTIONABLE_TARGETS lookup for each gene/variant
    |   Evidence level from knowledge graph
    |
    v
[Step 2: Biomarker-Driven Therapies]
    |   MSI-H -> pembrolizumab, nivolumab, dostarlimab (Level A)
    |   TMB-H (>=10 mut/Mb) -> pembrolizumab (Level A)
    |   HRD/BRCA -> olaparib, rucaparib, niraparib, talazoparib (Level A/B)
    |   PD-L1 TPS >=50% -> pembrolizumab first-line (Level A)
    |   NTRK fusion -> larotrectinib, entrectinib (Level A)
    |   + BIOMARKER_PANELS registry check
    |
    v
[Step 3: Evidence Level Sort]
    |   A (FDA-approved) > B (Clinical) > C (Case reports) > D > E
    |
    v
[Step 4: Resistance Check]
    |   RESISTANCE_MAP: mutation-level resistance
    |   _DRUG_CLASS_GROUPS: same-mechanism class resistance
    |
    v
[Step 5: Contraindication Check]
    |   Same drug previously used -> flag
    |   Same drug_class as prior failed therapy -> flag
    |
    v
[Step 6: Supporting Evidence + Combination Therapy]
    |   Search onco_therapies + onco_literature for each drug
    |   Known FDA-approved combos (dabrafenib+trametinib, etc.)
    |
    v
[Step 7: Final Ranking]
    Clean therapies first (sorted by evidence level)
    Flagged therapies after (resistance/contraindication)
    Assign rank 1..N
```

### 2. Trial Matcher (src/trial_matcher.py, 513 lines)

Hybrid deterministic + semantic clinical trial matching:

```
Patient Profile (cancer_type, biomarkers, stage, age)
    |
    v
[Step 1: Deterministic Filter]
    |   Cancer type (fuzzy via 18+ alias groups)
    |   Open statuses: Recruiting, Active, Enrolling by invitation
    |   Milvus filter expressions
    |
    v
[Step 2: Semantic Search]
    |   Embed eligibility query -> vector similarity search
    |
    v
[Step 3: Merge and Deduplicate]
    |   Union by trial_id, keep best score
    |
    v
[Step 4: Composite Scoring]
    |   biomarker_match (0.40) + semantic_score (0.25)
    |   + phase_weight (0.20) + status_weight (0.15)
    |   * age_penalty (1.0 or 0.5)
    |
    v
[Step 5: Explanation Generation]
    Ranked trial list with match rationale
```

### 3. Case Manager (src/case_manager.py, 516 lines)

Handles case lifecycle: creation from VCF or structured input, variant
annotation against ACTIONABLE_TARGETS, persistence to Milvus, and MTB packet
generation with variant tables, evidence summaries, therapy rankings, trial
matches, open questions, and formatted citations.

### 4. Query Expansion (src/query_expansion.py, 812 lines)

Domain-aware query expansion across 12 categories:

| Category | Example Input | Example Expansions |
|----------|---------------|-------------------|
| Cancer types | NSCLC | lung adenocarcinoma, EGFR-mutant lung |
| Genes | EGFR | L858R, exon 19 deletion, T790M, C797S |
| Therapies | osimertinib | Tagrisso, 3rd-gen EGFR TKI |
| Biomarkers | TMB | tumor mutational burden, mut/Mb |
| Pathways | MAPK | RAS-RAF-MEK-ERK, RTK signaling |
| Resistance | T790M | gatekeeper mutation, osimertinib |
| Clinical terms | PFS | progression-free survival, HR |
| Trial terms | Phase 3 | randomized, pivotal, registration trial |
| Immunotherapy | checkpoint | PD-1, PD-L1, CTLA-4, pembrolizumab |
| Surgery/radiation | lobectomy | surgical resection, VATS |
| Toxicity | pneumonitis | ILD, interstitial lung disease |
| Genomics | ctDNA | circulating tumor DNA, liquid biopsy |

---

## RAG Engine

**Module:** `src/rag_engine.py` (908 lines)
**Class:** `OncoRAGEngine`

### Search Flow

```
User Question
      |
      v
[Query Embedding] -- BGE-small-en-v1.5 with instruction prefix
      |
      v
[Parallel Collection Search] -- ThreadPoolExecutor (max 8 workers)
      |   Each collection searched with per-collection weight,
      |   filter_field, year_field
      |
      v
[Query Expansion Search] -- (optional) expand_query() -> re-embed -> search
      |
      v
[Merge and Rank] -- De-duplicate by ID, weight-adjusted score, top 30
      |
      v
[Knowledge Injection] -- gene, therapy, resistance, pathway, biomarker
      |
      v
[Prompt Assembly] -- Domain Knowledge + Evidence + Question + Instructions
      |
      v
[LLM Synthesis] -- Claude claude-sonnet-4-6 -> answer with citations
      |
      v
AgentResponse
```

### Relevance Classification

| Score Range | Classification |
|-------------|----------------|
| >= 0.85 | High |
| >= 0.65 | Medium |
| < 0.65 | Low |

### Evidence Evaluation Thresholds

```python
MIN_SUFFICIENT_HITS = 3
MIN_COLLECTIONS_FOR_SUFFICIENT = 2
MIN_SIMILARITY_SCORE = 0.30
```

Evidence items with scores below 0.30 are discarded. The verdict is:
- `sufficient` -- >= 3 quality hits from >= 2 distinct collections
- `partial` -- some hits but below threshold
- `insufficient` -- zero quality hits; triggers fallback queries

### Fallback Queries

When evidence is insufficient, the agent generates broader fallback queries:
- `"{gene} oncology therapeutic implications"`
- `"{gene} mutation clinical significance"`
- `"{cancer_type} precision medicine current landscape"`

---

## API Layer

### FastAPI Server (port 8527)

The API provides RESTful endpoints for all agent capabilities.

**Core endpoints:**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/healthz` | GET | Health check |
| `/readyz` | GET | Readiness check (Milvus connection) |
| `/metrics` | GET | Prometheus metrics |
| `/v1/query` | POST | RAG query with optional streaming |

**Case management:**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/cases` | POST | Create case from VCF or structured input |
| `/v1/cases/{id}` | GET | Retrieve case by ID |
| `/v1/cases/{id}/mtb` | GET | Generate MTB packet for case |

**Clinical analysis:**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/therapy-rank` | POST | Rank therapies for a patient profile |
| `/v1/trial-match` | POST | Match patient to clinical trials |
| `/v1/compare` | POST | Comparative evidence retrieval |

**Export:**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/export/markdown` | POST | Generate Markdown report |
| `/v1/export/json` | POST | Generate JSON export |
| `/v1/export/pdf` | POST | Generate PDF report |
| `/v1/export/fhir` | POST | Generate FHIR R4 bundle |

### Authentication and Security

- API key authentication via `ONCO_API_KEY` environment variable
- CORS middleware with configurable allowed origins
- Request validation via Pydantic v2 models
- Input sanitization for Milvus filter expressions

### Health Checks

- Docker health check: `/healthz` every 30 seconds
- Start period: 30 seconds
- Milvus connectivity verified on startup and at readiness check

---

## UI Layer

### Streamlit Application (port 8526)

The UI provides a multi-tab interface for molecular tumor board workflows:

| Tab | Name | Functionality |
|-----|------|--------------|
| 1 | Clinical Query | RAG search with streaming, citations, comparative mode |
| 2 | Case Manager | Create/view cases, VCF upload, variant annotation |
| 3 | Variant Viewer | Interactive variant table with evidence levels |
| 4 | Therapy Ranking | Ranked therapies with resistance flags, combos |
| 5 | Trial Matching | Matched trials with composite scores, eligibility |
| 6 | MTB Packet | Full tumor board packet generation and export |
| 7 | Pathway Explorer | Oncogenic pathway visualization and druggable nodes |
| 8 | Resistance Analysis | Resistance mechanisms and alternative therapies |
| 9 | Knowledge Base | Collection stats, knowledge graph coverage |

---

## Export Pipeline

**Module:** `src/export.py` (1,055 lines)

Four export formats, each accepting dict, MTBPacket, or string input.

### Markdown Export

Sections: Header (patient/meta), Clinical Summary, Somatic Variant Profile
(table), Biomarker Summary, Evidence Summary, Therapy Ranking (table),
Clinical Trial Matches, Pathway Context, Known Resistance Mechanisms,
Open Questions, Disclaimer.

### JSON Export

Standardized schema with `meta` block (format, version, generated_at,
pipeline, author) plus all clinical sections. Suitable for downstream
programmatic consumption and integration.

### PDF Export

NVIDIA-themed PDF via ReportLab with:
- Green header bar (RGB 118, 185, 0) with white title
- Structured tables for variants, therapies, trials
- Alternating row colors (whitesmoke/white)
- Footer disclaimer in gray 7pt text

### FHIR R4 Export

Generates a FHIR R4 Bundle (type=collection) containing:

| Resource | Content |
|----------|---------|
| Patient | Identifier with urn:hcls-ai-factory:patient |
| Observation (N) | One per variant (LOINC 69548-6) |
| Observation (TMB) | Tumor mutation burden (LOINC 94076-7) |
| Observation (MSI) | Microsatellite instability (LOINC 81695-9) |
| Specimen | Tumor tissue (SNOMED 119376003) |
| Condition | Cancer diagnosis (SNOMED-coded, 22 cancer types) |
| MedicationRequest | Therapy recommendations (top 10) |
| DiagnosticReport | Master genomic report (LOINC 81247-9) |

---

## Scaling and Performance

### Resource Footprint on DGX Spark

| Component | Resource Usage |
|-----------|---------------|
| Milvus standalone | ~1-3 GB RAM (11 collections, 609 vectors) |
| BGE-small-en-v1.5 | ~200 MB RAM (CPU inference) |
| FastAPI server | ~400 MB RAM (uvicorn workers) |
| Streamlit UI | ~300 MB RAM |
| Total agent footprint | ~2-4 GB RAM |

### Performance Characteristics

| Operation | Typical Latency |
|-----------|----------------|
| Single collection search | 10-50 ms |
| 11-collection parallel search | 40-150 ms |
| BGE-small embedding (single) | 5-15 ms |
| Query expansion | < 5 ms |
| Claude LLM synthesis | 2-8 seconds |
| Therapy ranking (full pipeline) | 200-800 ms |
| Trial matching (hybrid) | 100-500 ms |
| MTB packet generation | 3-15 seconds |
| PDF generation | 200-500 ms |
| FHIR R4 export | < 100 ms |

### Parallelization

- Collection searches run in parallel via ThreadPoolExecutor (max_workers=8)
- Sub-question searches run in parallel during the SEARCH phase
- Embedding batching (batch_size=32) for bulk seed operations
- APScheduler for background metric collection and housekeeping

### Caching

- Embedding cache for frequently used queries
- Knowledge graph: In-memory dictionaries, loaded at startup
- Query expansion maps: In-memory, loaded at module import

---

## Security

### Authentication

- API key authentication via `ONCO_API_KEY` environment variable
- Empty API key disables authentication (development mode only)
- Bearer token format: `Authorization: Bearer <api_key>`

### Data Protection

- All patient data encrypted at rest via Milvus storage encryption
- De-identified patient case snapshots (no PHI in onco_cases)
- Anthropic API key stored in environment variables, never committed to source
- No PHI stored in application logs

### Input Validation

- Pydantic v2 model validation on all API inputs (13 enums enforce valid values)
- Milvus filter expression sanitization to prevent injection
- VCF file validation before parsing (cyvcf2 error handling)
- CORS middleware with explicit origin allowlist

### Network Security

- All services communicate via Docker bridge network
- Only Streamlit (8526) and FastAPI (8527) ports are exposed externally
- Milvus gRPC (19530) is internal to the Docker network
- Non-root container execution

---

## File Structure

```
precision_oncology_agent/
|-- api/                         # FastAPI REST server
|   |-- main.py                  # Entry point, endpoints
|
|-- app/                         # Streamlit UI
|   |-- oncology_ui.py           # Multi-tab MTB application
|
|-- config/                      # Configuration
|   |-- settings.py              # OncologySettings (Pydantic, ONCO_ prefix)
|
|-- data/                        # Seed data
|   |-- variant_seed_data.json         (90 records)
|   |-- literature_seed_data.json      (60 records)
|   |-- therapy_seed_data.json         (64 records)
|   |-- guideline_seed_data.json       (45 records)
|   |-- trial_seed_data.json           (55 records)
|   |-- biomarker_seed_data.json       (50 records)
|   |-- resistance_seed_data.json      (50 records)
|   |-- pathway_seed_data.json         (35 records)
|   |-- outcome_seed_data.json         (40 records)
|   +-- cases_seed_data.json           (37 records)
|
|-- docs/
|   |-- ARCHITECTURE_GUIDE.md    # This document
|   |-- DEMO_GUIDE.md
|   |-- DEPLOYMENT_GUIDE.md
|   |-- DESIGN.md
|   |-- INDEX.md
|   |-- LEARNING_GUIDE_ADVANCED.md
|   |-- LEARNING_GUIDE_FOUNDATIONS.md
|   |-- PROJECT_BIBLE.md
|   +-- WHITE_PAPER.md
|
|-- src/                         # Core engine modules
|   |-- __init__.py
|   |-- agent.py                 (553 lines) Plan-search-evaluate-synthesize
|   |-- case_manager.py          (516 lines) VCF parsing, case CRUD, MTB packets
|   |-- collections.py           Milvus collection management (11 schemas)
|   |-- cross_modal.py           Cross-agent event handling
|   |-- export.py                (1,055 lines) Markdown/JSON/PDF/FHIR R4
|   |-- knowledge.py             (1,662 lines) 5 knowledge domains
|   |-- metrics.py               Prometheus metrics
|   |-- models.py                (538 lines) 13 enums, 10 domain models
|   |-- query_expansion.py       (812 lines) 12-category domain expansion
|   |-- rag_engine.py            (908 lines) Multi-collection RAG engine
|   |-- scheduler.py             APScheduler background tasks
|   |-- therapy_ranker.py        (748 lines) 7-step therapy ranking
|   |-- trial_matcher.py         (513 lines) Hybrid trial matching
|   +-- utils/
|       +-- vcf_parser.py        VCF parsing via cyvcf2
|   +-- ingest/                  Data ingestion scripts
|   +-- workflows/               Workflow definitions
|
|-- tests/                       # 10 test files, 556 tests
|
|-- docker-compose.yml           6-service stack
|-- Dockerfile
|-- requirements.txt             23 dependencies
|-- pyproject.toml
|-- LICENSE
+-- README.md
```

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Oncology Agent is a clinical decision support research tool for oncology. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
