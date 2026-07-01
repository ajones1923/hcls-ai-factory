# Precision Oncology Intelligence Agent - Architecture Guide

## Technical Deep-Dive: System Design and Implementation

**Author:** Adam Jones
**Date:** March 2026
**Version:** 0.1.0

---

## Table of Contents

1. [Overview](#1-overview)
2. [High-Level Architecture](#2-high-level-architecture)
3. [Component Deep Dives](#3-component-deep-dives)
   - 3.1 [RAG Engine](#31-rag-engine)
   - 3.2 [Intelligence Agent](#32-intelligence-agent)
   - 3.3 [Knowledge Graph](#33-knowledge-graph)
   - 3.4 [Collections Manager](#34-collections-manager)
   - 3.5 [Case Manager](#35-case-manager)
   - 3.6 [Therapy Ranker](#36-therapy-ranker)
   - 3.7 [Trial Matcher](#37-trial-matcher)
   - 3.8 [Export Module](#38-export-module)
   - 3.9 [Query Expansion](#39-query-expansion)
   - 3.10 [Cross-Modal Triggers](#310-cross-modal-triggers)
   - 3.11 [Metrics and Monitoring](#311-metrics-and-monitoring)
   - 3.12 [Ingest Pipeline](#312-ingest-pipeline)
   - 3.13 [Data Models](#313-data-models)
4. [Data Flow Diagrams](#4-data-flow-diagrams)
5. [Collection Architecture](#5-collection-architecture)
6. [API Layer](#6-api-layer)
7. [UI Architecture](#7-ui-architecture)
8. [Docker and Deployment](#8-docker-and-deployment)
9. [Testing Architecture](#9-testing-architecture)
10. [Security Considerations](#10-security-considerations)
11. [Performance Characteristics](#11-performance-characteristics)
12. [Extension Points](#12-extension-points)

---

## 1. Overview

The Precision Oncology Intelligence Agent is a RAG-powered clinical decision
support system for Molecular Tumor Boards (MTBs). It operates as Stage 2b
of the HCLS AI Factory pipeline, sitting between the genomics pipeline
(Stage 1, Parabricks/DeepVariant) and the drug discovery pipeline (Stage 3,
BioNeMo MolMIM/DiffDock).

### Codebase Inventory

| Layer              | Files | Lines  |
|--------------------|-------|--------|
| Source modules     | 13    | 11,440 |
| Ingest parsers     | 9 + base class | 1,793 |
| Utility modules    | 2     | 668    |
| API routes + main  | 5 + main | 1,300 |
| Streamlit UI       | 1     | 758    |
| Test files         | 10 + conftest | 4,585 |
| Setup/seed scripts | 17    | 2,273  |

### Design Principles

1. **Multi-collection RAG** -- Eleven specialized vector collections, each
   with a tunable weight, contribute to every evidence retrieval. This
   prevents over-reliance on any single data source.

2. **Plan-search-evaluate-synthesize** -- The intelligence agent decomposes
   clinical queries into structured search plans, executes parallel searches,
   evaluates sufficiency, retries if needed, and synthesizes a cited answer.

3. **Deterministic + semantic hybrid** -- Hard filters (cancer type, gene,
   recruiting status) reduce the search space before semantic similarity
   ranking. This ensures clinically relevant results without sacrificing
   recall.

4. **Knowledge-augmented generation** -- A curated knowledge graph (40+
   actionable targets, 30 drugs, 10 pathways) injects domain context into
   LLM prompts alongside retrieved evidence.

5. **Export interoperability** -- Reports are exported in four formats:
   Markdown, JSON, PDF (NVIDIA-branded), and FHIR R4 Bundles for EHR
   integration.

---

## 2. High-Level Architecture

```
+--------------------------------------------------------------------+
|                     HCLS AI Factory Pipeline                       |
|                                                                    |
|   Stage 1               Stage 2b                  Stage 3          |
|   Genomics              Oncology Agent             Drug Discovery   |
|   (Parabricks)          (This System)              (BioNeMo)       |
|                                                                    |
|   FASTQ -> VCF  ------> Variant Interpretation --> MolMIM/DiffDock |
|                          Therapy Ranking                            |
|                          Trial Matching                             |
+--------------------------------------------------------------------+

+====================================================================+
|              Oncology Intelligence Agent Architecture               |
+====================================================================+
|                                                                    |
|  +--------------+    +----------------+    +--------------------+  |
|  |  Streamlit   |    |   FastAPI       |    |  Cross-Modal       |  |
|  |  UI (:8526)  |--->|   API (:8527)   |<-->|  Event Bus         |  |
|  |              |    |                |    |  (other agents)     |  |
|  |  5 tabs:     |    |  5 route files |    +--------------------+  |
|  |  - Cases     |    |  - meta_agent  |                            |
|  |  - Evidence  |    |  - cases       |                            |
|  |  - Trials    |    |  - trials      |    +--------------------+  |
|  |  - Therapy   |    |  - reports     |    |  Ingest Scheduler   |  |
|  |  - Outcomes  |    |  - events      |    |  (APScheduler)      |  |
|  +--------------+    +-------+--------+    +--------------------+  |
|                              |                                     |
|          +-------------------+--------------------+                |
|          |                   |                    |                 |
|  +-------v-------+  +-------v--------+  +-------v--------+       |
|  | Intelligence   |  | Case Manager   |  | Therapy Ranker |       |
|  | Agent          |  |                |  |                |       |
|  | (plan/search/  |  | VCF parsing    |  | Evidence-level |       |
|  |  evaluate/     |  | MTB packets    |  | ranking        |       |
|  |  synthesize)   |  | Actionability  |  | Resistance     |       |
|  +-------+--------+  +-------+--------+  +-------+--------+       |
|          |                   |                    |                 |
|  +-------v-------------------v--------------------v--------+       |
|  |                    RAG Engine                            |       |
|  |  - 11 collection configs with weights                   |       |
|  |  - ThreadPoolExecutor for parallel search               |       |
|  |  - Knowledge injection from knowledge.py                |       |
|  |  - Query expansion from query_expansion.py              |       |
|  |  - Citation formatting with relevance scoring           |       |
|  +---------------------------+------------------------------+       |
|                              |                                     |
|  +---------------------------v------------------------------+       |
|  |              Collections Manager (pymilvus)              |       |
|  |  11 collections, 384-dim FLOAT_VECTOR, IVF_FLAT/COSINE  |       |
|  +---------------------------+------------------------------+       |
|                              |                                     |
+------------------------------v-------------------------------------+
                               |
              +----------------v-----------------+
              |        Milvus 2.4 Standalone      |
              |  etcd (metadata) + MinIO (storage)|
              +----------------------------------+
```

### Service Initialization Order

At startup, the FastAPI lifespan handler creates nine services in sequence:

```
1. OncoSettings            (config/settings.py -- Pydantic BaseSettings)
2. OncoCollectionManager   (connect to Milvus)
3. EmbedderWrapper         (SentenceTransformer BGE-small-en-v1.5)
4. OncoRAGEngine           (collection_manager + embedder + settings)
5. OncoIntelligenceAgent   (rag_engine)
6. OncologyCaseManager     (collection_manager + embedder + knowledge + rag_engine)
7. TrialMatcher            (collection_manager + embedder)
8. TherapyRanker           (collection_manager + embedder + knowledge)
9. OncoCrossModalTrigger   (collection_manager + embedder + cross-modal settings)
```

All nine services are stored in a shared `_state` dict and accessed by route
handlers via `get_state()`.

---

## 3. Component Deep Dives

### 3.1 RAG Engine

**File:** `src/rag_engine.py` (908 lines)
**Class:** `OncoRAGEngine`

The RAG engine is the central retrieval-and-generation orchestrator. It
searches across all 11 Milvus collections in parallel, merges results using
weighted scoring, injects curated knowledge context, and generates LLM
answers with formatted citations.

#### Collection Configuration

Each collection has a weight that determines its contribution to the final
ranked result set. Weights sum to approximately 1.0 and are configurable
via environment variables with the `ONCO_` prefix:

| Collection         | Weight | Label       | Filter Field | Year Field |
|--------------------|--------|-------------|-------------|------------|
| onco_variants      | 0.18   | Variant     | gene        | --         |
| onco_literature    | 0.16   | Literature  | gene        | year       |
| onco_therapies     | 0.14   | Therapy     | gene        | --         |
| onco_guidelines    | 0.12   | Guideline   | gene        | --         |
| onco_trials        | 0.10   | Trial       | gene        | --         |
| onco_biomarkers    | 0.08   | Biomarker   | gene        | --         |
| onco_resistance    | 0.07   | Resistance  | gene        | --         |
| onco_pathways      | 0.06   | Pathway     | gene        | --         |
| onco_outcomes      | 0.04   | Outcome     | gene        | --         |
| onco_cases         | 0.02   | Case        | gene        | --         |
| genomic_evidence   | 0.03   | Genomic     | gene        | --         |

#### System Prompt

The engine uses a detailed system prompt defining 8 core competencies:

1. **Molecular profiling** -- somatic/germline variant interpretation, TMB,
   MSI, CNA, gene fusions
2. **Variant interpretation** -- CIViC/OncoKB evidence levels (Tier I-IV),
   AMP/ASCO/CAP classification
3. **Therapy selection** -- NCCN/ESMO guideline-concordant recommendations,
   FDA-approved indications
4. **Clinical trial matching** -- eligibility assessment against active
   ClinicalTrials.gov registrations, basket/umbrella trial awareness
5. **Resistance mechanisms** -- on-target mutations, bypass signaling,
   lineage plasticity, actionable resistance biomarkers
6. **Biomarker assessment** -- TMB, MSI, PD-L1, HRD scoring; companion
   diagnostic requirements
7. **Outcomes monitoring** -- RECIST criteria, survival endpoints, MRD
   tracking, ctDNA dynamics
8. **Cross-modal integration** -- linking genomic findings to imaging,
   pathology, and drug-discovery pipelines

The prompt also embeds 5 behavioral instructions: cite evidence, think
cross-functionally, highlight resistance and contraindications, reference
guidelines, and acknowledge uncertainty.

#### Key Methods

**`retrieve(question, top_k, filters)`**

Executes parallel vector searches across all 11 collections using a
`ThreadPoolExecutor`. Each collection search runs as a separate thread.
Results are merged by applying the collection-specific weight to each hit's
similarity score, then sorted by weighted score descending.

**`cross_collection_search(question, collections, top_k)`**

Targeted search across a specified subset of collections. Used by the
intelligence agent when a search plan identifies specific domains (e.g.,
only search variants and resistance for a resistance-mechanism query).

**`comparative_search(entity_a, entity_b, top_k)`**

Retrieves evidence for two entities (e.g., two drugs, two genes) and
structures the results for side-by-side comparison. Used for comparative
therapy analysis.

**`synthesize(question, evidence, context)`**

Constructs a prompt combining the system prompt, knowledge-graph context,
retrieved evidence with citations, and the user's question. Sends to the
Claude LLM (claude-sonnet-4-20250514 by default) and returns the generated
answer.

#### Citation Formatting

Citations are formatted with relevance tiers based on configurable
thresholds:

- **Strong** (score >= 0.75): High-confidence evidence
- **Moderate** (score >= 0.60): Supporting evidence
- **Weak** (score < 0.60): Background context

Each citation includes the source collection label, relevance score,
and a reference identifier linking back to the original document.

---

### 3.2 Intelligence Agent

**File:** `src/agent.py` (553 lines)
**Class:** `OncoIntelligenceAgent`

The intelligence agent implements the plan-search-evaluate-synthesize
pattern. It decomposes clinical oncology queries into structured search
plans, executes those plans via the RAG engine, evaluates evidence
sufficiency, and generates cited answers.

#### SearchPlan Dataclass

```python
@dataclass
class SearchPlan:
    question: str
    search_strategy: str          # "broad", "targeted", or "comparative"
    identified_topics: List[str]
    target_genes: List[str]
    relevant_cancer_types: List[str]
    sub_questions: List[str]
```

The search plan drives how the RAG engine is invoked. A "broad" strategy
searches all 11 collections. A "targeted" strategy focuses on the 2-3
most relevant collections. A "comparative" strategy retrieves evidence for
multiple entities for side-by-side analysis.

#### Oncology Vocabulary

The agent maintains curated vocabularies for entity recognition:

- **KNOWN_GENES (30 genes):** BRAF, EGFR, ALK, ROS1, KRAS, HER2, NTRK,
  RET, MET, FGFR, PIK3CA, IDH1, IDH2, BRCA, BRCA1, BRCA2, TP53, PTEN,
  CDKN2A, STK11, ESR1, ERBB2, NRAS, APC, VHL, KIT, PDGFRA, FLT3, NPM1, DNMT3A
- **KNOWN_CANCER_TYPES (25 types):** NSCLC through MESOTHELIOMA
- **_CANCER_ALIASES (30+ mappings):** Natural-language mentions mapped to
  canonical types (e.g., "lung" -> NSCLC, "tnbc" -> BREAST,
  "gbm" -> GLIOBLASTOMA, "crc" -> COLORECTAL)

#### Agent Workflow

```
User Query
    |
    v
1. PLAN: Extract genes, cancer types, topics from query text
         Match against KNOWN_GENES and _CANCER_ALIASES
         Choose strategy: broad / targeted / comparative
         Generate sub-questions if multi-faceted
    |
    v
2. SEARCH: Execute plan via RAG engine
           Parallel cross-collection vector search
           Apply gene/cancer_type filters where available
    |
    v
3. EVALUATE: Assess evidence sufficiency
             - Sufficient: >= 3 hits from >= 2 distinct collections
             - Partial: some evidence, but limited coverage
             - Insufficient: too few or too low-scoring results
    |
    v
4. RETRY (up to 2x): If insufficient evidence, broaden search:
                      - Relax filters
                      - Expand query terms
                      - Search additional collections
    |
    v
5. SYNTHESIZE: Inject knowledge context + evidence into LLM prompt
               Generate cited clinical answer
               Format with evidence tiers and references
```

#### Evidence Evaluation

The agent classifies retrieved evidence into three categories using
configurable thresholds (from `OncoSettings`):

| Category     | Criteria                                          |
|-------------|---------------------------------------------------|
| Sufficient  | >= 3 hits AND >= 2 distinct collections represented |
| Partial     | Some hits, but fewer than 3 or from a single collection |
| Insufficient| Too few hits or all below minimum similarity (0.30) |

When evidence is insufficient, the agent retries up to 2 times with
progressively broader search parameters.

---

### 3.3 Knowledge Graph

**File:** `src/knowledge.py` (1,662 lines)

The knowledge graph is a curated, Python-native domain knowledge store
that augments RAG retrieval with structured oncology intelligence. It is
not a database -- it is a set of Python dictionaries and helper functions
compiled directly into the application.

#### ACTIONABLE_TARGETS (40+ genes)

Each entry contains:
- `gene` -- canonical gene symbol
- `full_name` -- expanded gene name
- `cancer_types` -- list of relevant cancer types
- `key_variants` -- clinically significant variants
- `targeted_therapies` -- approved single-agent therapies
- `combination_therapies` -- approved combinations
- `resistance_mutations` -- known resistance mechanisms
- `pathway` -- primary signaling pathway
- `evidence_level` -- A/B/C/D/E tier
- `description` -- narrative summary

Example (BRAF):
```
Gene:         BRAF (B-Raf Proto-Oncogene)
Cancers:      melanoma, NSCLC, colorectal, thyroid, hairy cell leukemia
Variants:     V600E, V600K, V600D, V600R, class II, class III
Therapies:    vemurafenib, dabrafenib, encorafenib
Combinations: dabrafenib + trametinib, encorafenib + binimetinib
Resistance:   MEK1/2 mutations, NRAS activation, BRAF amplification
Pathway:      MAPK
Evidence:     A
```

#### THERAPY_MAP (30 drugs)

Maps drug names to structured records containing:
- Mechanism of action (MOA)
- FDA-approved indications
- Associated clinical trials
- Dosing and scheduling notes
- Key adverse events

#### RESISTANCE_MAP

Maps drug classes to resistance mechanisms (primary mutations, bypass
pathways, lineage plasticity) with recommended bypass strategies.

#### PATHWAY_MAP (10 pathways)

Curated definitions for MAPK, PI3K/AKT/mTOR, DDR, Cell Cycle,
WNT/beta-catenin, Hedgehog, JAK/STAT, Notch. Each entry includes
key nodes, druggable targets, and cross-talk relationships.

#### BIOMARKER_PANELS (20 biomarkers)

Each biomarker includes testing methods, scoring cutoffs, associated
therapies, and companion diagnostic requirements.

#### ENTITY_ALIASES (30+ aliases)

Maps alternate names and abbreviations to canonical identifiers for
entity resolution (e.g., "trastuzumab" and "Herceptin" resolve to
the same therapy).

#### Helper Functions

- `get_target_context(gene)` -- Returns formatted context string for a gene
- `get_therapy_context(drug)` -- Returns therapy details as context
- `get_resistance_context(drug_class)` -- Returns resistance information
- `get_pathway_context(pathway)` -- Returns pathway details
- `classify_variant_actionability(gene, variant)` -- Returns A/B/C/VUS
  classification based on ACTIONABLE_TARGETS lookup

---

### 3.4 Collections Manager

**File:** `src/collections.py` (665 lines)
**Class:** `OncoCollectionManager`

The collections manager wraps `pymilvus` operations for all 11 Milvus
collections. It defines schemas, manages connections, and provides
search/insert abstractions.

#### Collection Schemas

All 11 collections share these common properties:
- **Embedding dimension:** 384 (BGE-small-en-v1.5)
- **Vector field type:** FLOAT_VECTOR
- **Index type:** IVF_FLAT
- **Similarity metric:** COSINE
- **Index parameters:** nlist=1024, nprobe=16
- **Primary key:** VARCHAR id field (max 100 chars)

Each collection schema includes `id` (VARCHAR primary key), `embedding`
(FLOAT_VECTOR dim=384), `text_chunk` (VARCHAR 3000), and domain-specific
metadata fields. For example, `onco_literature` adds `title`, `source_type`,
`year`, `cancer_type`, `gene`, `variant`, `keywords`, `journal`.
`onco_trials` adds `nct_id`, `title`, `phase`, `status`, `cancer_type`,
`biomarkers`. The `genomic_evidence` collection is read-only, shared with
the Stage 1 genomics pipeline.

#### Key Methods

- `connect()` / `disconnect()` -- Milvus connection lifecycle
- `is_connected()` -- Connection health check
- `create_collection(name)` -- Creates collection with schema and index
- `list_collections()` -- Returns all oncology collection names
- `get_collection_count(name)` -- Entity count for a collection
- `search(collection, embedding, top_k, filters)` -- Vector similarity search
  with optional metadata filters
- `insert(collection, records)` -- Batch insert with embedding vectors
- `get_collection_stats()` -- Aggregate statistics across all collections

#### Index Configuration

IVF_FLAT with nlist=1024, nprobe=16, COSINE metric. This gives ~98.5%
recall vs brute-force while reducing query latency by ~60x on 10K+ vectors.

---

### 3.5 Case Manager

**File:** `src/case_manager.py` (516 lines)
**Class:** `OncologyCaseManager`

The case manager handles the patient case lifecycle from VCF ingestion to
MTB packet generation. It is the primary interface for clinical workflow
integration.

#### VCF Parsing

The case manager delegates VCF parsing to `src/utils/vcf_parser.py`
(371 lines), which supports multiple annotation formats:

- **SnpEff ANN=** format: Extracts gene from position 4 in pipe-delimited
  ANN field. Pattern: `ANN=Allele|Annotation|Impact|Gene|...`
- **VEP CSQ=** format: Extracts SYMBOL from the CSQ field.
  Pattern: `CSQ=Allele|Consequence|...|SYMBOL|...`
- **Simple tags**: Falls back to `GENE=` or `GENEINFO=` INFO field tags.

Supports VCF 4.x as produced by GATK, DeepVariant, Strelka, Mutect2,
and Parabricks pipelines.

#### Variant Actionability Classification

Each parsed variant is classified into one of four tiers by looking up
the gene and variant against the ACTIONABLE_TARGETS knowledge graph:

| Tier | Meaning                                              |
|------|------------------------------------------------------|
| A    | FDA-approved companion diagnostic, validated target  |
| B    | Clinical evidence from well-powered studies          |
| C    | Case reports, small series, or emerging data         |
| VUS  | Variant of uncertain significance, no matching target|

Classification logic:
1. Look up gene in ACTIONABLE_TARGETS
2. If found, check if the specific variant matches a `key_variants` entry
3. If variant matches, return the gene's `evidence_level`
4. If gene found but variant not in key_variants, return "C"
5. If gene not found, return "VUS"

#### MTB Packet Generation

The `generate_mtb_packet()` method produces a structured Molecular Tumor
Board packet containing:

```
MTBPacket:
  +-- Patient metadata (anonymized ID, cancer type, stage)
  +-- Variant summary (all parsed variants with annotations)
  +-- Actionable variants (Tier A/B/C, with evidence)
  +-- Therapy ranking (from TherapyRanker)
  +-- Trial matches (from TrialMatcher)
  +-- Resistance alerts (from RESISTANCE_MAP)
  +-- Biomarker status (MSI, TMB, PD-L1, HRD)
  +-- Open questions for MTB discussion
  +-- Evidence citations with source links
```

#### Input Validation

Filter values are validated against `_SAFE_FILTER_RE` pattern
(`^[A-Za-z0-9 _.\-/]+$`) to prevent Milvus expression injection attacks.

---

### 3.6 Therapy Ranker

**File:** `src/therapy_ranker.py` (748 lines)
**Class:** `TherapyRanker`

The therapy ranker produces evidence-ranked therapy recommendations with
resistance awareness and contraindication detection.

#### Ranking Strategy

```
Input: CaseSnapshot (variants, biomarkers, prior therapies)
    |
    v
Step 1: Variant-Driven Therapies
    - For each actionable variant, look up ACTIONABLE_TARGETS
    - Retrieve targeted_therapies and combination_therapies
    - Assign evidence level from the target entry
    |
    v
Step 2: Biomarker-Driven Therapies
    - MSI-H    -> pembrolizumab (immunotherapy)
    - TMB-H    -> pembrolizumab (immunotherapy)
    - HRD      -> PARP inhibitors (olaparib, rucaparib)
    - PD-L1+   -> pembrolizumab, nivolumab, atezolizumab
    - HER2+    -> trastuzumab, T-DXd
    |
    v
Step 3: Evidence-Level Ranking
    - Sort by evidence strength: A(0) > B(1) > C(2) > D(3) > E(4) > VUS(5)
    - Within the same level, sort by number of supporting references
    |
    v
Step 4: Resistance Checking
    - For each candidate therapy, look up RESISTANCE_MAP
    - If the patient has a known resistance mutation, flag the therapy
    - Include bypass pathway suggestions from the resistance entry
    |
    v
Step 5: Contraindication Detection
    - If patient previously received a therapy in the same drug class
      and progressed, flag as contraindicated
    - Example: prior erlotinib failure flags other 1st-gen EGFR TKIs
    |
    v
Step 6: Evidence Retrieval
    - For each ranked therapy, search onco_therapies and onco_literature
    - Attach supporting citations with relevance scores
    |
    v
Output: Ranked list of TherapyRecommendation objects
```

#### Evidence Level Ordering

```python
EVIDENCE_LEVEL_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "VUS": 5}
```

Level A therapies (FDA-approved with companion diagnostics) always rank
above Level B (well-powered clinical studies), and so on. This ordering
is deterministic -- the ranker never promotes a lower-evidence therapy
above a higher-evidence one unless the higher-evidence therapy is flagged
for resistance.

---

### 3.7 Trial Matcher

**File:** `src/trial_matcher.py` (513 lines)
**Class:** `TrialMatcher`

The trial matcher uses a hybrid deterministic + semantic approach to
match oncology patients to clinical trials.

#### Matching Strategy

```
Input: CaseSnapshot (cancer_type, biomarkers, stage, age)
    |
    v
1. DETERMINISTIC FILTER
    - Filter onco_trials by cancer_type (fuzzy match)
    - Retain only open statuses: Recruiting, Active not recruiting,
      Enrolling by invitation, Not yet recruiting
    |
    v
2. SEMANTIC SEARCH
    - Construct query: "{cancer_type} {biomarkers} {stage}"
    - Embed query with BGE-small-en-v1.5
    - Vector search against onco_trials collection
    |
    v
3. COMPOSITE SCORING
    Each trial receives a weighted composite score:

    composite = biomarker_score * 0.40
              + semantic_score  * 0.25
              + phase_score     * 0.20
              + status_score    * 0.15
    |
    v
4. AGE FILTERING
    - If trial has min/max age criteria, exclude ineligible patients
    |
    v
5. RANKING AND EXPLANATION
    - Sort by composite score descending
    - Generate structured explanation for each match
    - Include NCT ID, title, phase, status, and match rationale
```

#### Score Component Weights

| Component     | Weight | Description                              |
|---------------|--------|------------------------------------------|
| Biomarker     | 0.40   | Overlap between patient and trial biomarkers |
| Semantic      | 0.25   | Cosine similarity from vector search     |
| Phase         | 0.20   | Trial phase (Phase 3 > Phase 2 > Phase 1)|
| Status        | 0.15   | Recruitment status                       |

#### Phase Weights

| Phase       | Weight |
|------------|--------|
| Phase 3    | 1.0    |
| Phase 2/3  | 0.9    |
| Phase 2    | 0.8    |
| Phase 1/2  | 0.7    |
| Phase 1    | 0.6    |
| Phase 4    | 0.5    |

#### Status Weights

| Status                    | Weight |
|--------------------------|--------|
| Recruiting               | 1.0    |
| Enrolling by invitation  | 0.8    |
| Active, not recruiting   | 0.6    |
| Not yet recruiting       | 0.4    |

---

### 3.8 Export Module

**File:** `src/export.py` (1,055 lines)

The export module generates clinical reports in four formats for different
consumption contexts.

#### export_markdown()

Generates a structured Markdown document with clinical report sections:
- Header with patient ID, cancer type, date
- Actionable variant summary table
- Therapy recommendations with evidence levels
- Trial matches with NCT IDs
- Biomarker status panel
- Resistance alerts
- Evidence citations

#### export_json()

Produces a structured JSON document with full metadata:
- Machine-parseable variant records
- Nested therapy recommendation objects
- Trial match records with composite scores
- Timestamp and version metadata

#### export_pdf()

Generates a branded PDF using ReportLab:
- NVIDIA brand color: RGB (118, 185, 0)
- Configurable via `PDF_BRAND_COLOR_R/G/B` environment variables
- Header with logo area (height 50px)
- Page margins: 40px
- Evidence level labels with narrative descriptions:
  - Level 1: FDA-approved / Standard of Care
  - Level 2: Clinical Evidence / Consensus
  - Level 3: Case Reports / Early Trials
  - Level 4: Preclinical / Biological Rationale
  - Level R: Resistance Evidence

#### export_fhir_r4()

Generates a FHIR R4 Bundle for EHR integration:

**Resource types in the Bundle:**
- `Patient` -- Anonymized patient demographics
- `DiagnosticReport` -- Master genomic report (SNOMED coded)
- `Observation` -- Individual variant observations (LOINC genomic codes)

**LOINC Codes:**

| Observation          | LOINC Code | Description                    |
|---------------------|------------|--------------------------------|
| Genomic Report      | 81247-9    | Master HL7 genomic report      |
| Gene Studied        | 48018-6    | Gene studied                   |
| Variant             | 69548-6    | Genetic variant assessment     |
| Therapeutic Implication | 51969-4 | Genetic analysis summary      |
| Tumor Mutation Burden | 94076-7  | TMB                            |
| MSI                 | 81695-9    | Microsatellite instability     |

**SNOMED Cancer Codes (subset):**

| Cancer Type   | SNOMED Code | Display Name                        |
|--------------|-------------|-------------------------------------|
| NSCLC        | 254637007   | Non-small cell lung cancer          |
| Breast       | 254837009   | Malignant neoplasm of breast        |
| Colorectal   | 363406005   | Malignant tumor of colon            |
| Melanoma     | 372244006   | Malignant melanoma                  |

---

### 3.9 Query Expansion

**File:** `src/query_expansion.py` (812 lines)

The query expansion module improves RAG retrieval recall by expanding
oncology-specific keywords with semantically related terms.

#### Expansion Categories

The module maintains expansion dictionaries for 12 oncology domains:
cancer types, gene names, therapy names, biomarker terms, pathway terms,
resistance terms, clinical terms, trial terms, immunotherapy terms,
surgery/radiation terms, toxicity terms, and genomics terms.

Each domain maps canonical terms to 3-6 semantically related alternatives.
For example, NSCLC expands to "non-small cell lung cancer",
"lung adenocarcinoma", "lung squamous cell", "EGFR-mutant lung",
"ALK-positive lung".

The module also performs query rewriting -- ambiguous clinical shorthand
is expanded to explicit search terms. For example, "EGFR TKI resistance"
is rewritten to include "T790M" and "C797S" as specific variant searches.

---

### 3.10 Cross-Modal Triggers

**File:** `src/cross_modal.py` (383 lines)
**Class:** `OncoCrossModalTrigger`

The cross-modal trigger system connects the oncology agent to other agents
in the HCLS AI Factory (imaging intelligence, drug discovery) via event-
driven integration.

#### Trigger Logic

When actionable variants with evidence level A or B are identified, the
cross-modal trigger fires and queries across modalities:

```
Variant Detected (Level A/B)
    |
    v
1. Query genomic_evidence collection (Stage 1 pipeline data)
   - top_k: configurable (default 5 via GENOMIC_TOP_K)
   - Retrieves VCF-derived evidence from Parabricks output
    |
    v
2. Query imaging collections (if available)
   - top_k: configurable (default 5 via IMAGING_TOP_K)
   - Retrieves correlated imaging findings
    |
    v
3. Generate enrichment summary
   - Combines genomic and imaging context
   - Formats as human-readable summary for MTB review
```

#### CrossModalResult Dataclass

```python
@dataclass
class CrossModalResult:
    trigger_reason: str           # Why the trigger fired
    actionable_variants: List     # Variants that triggered it
    genomic_context: List         # Genomic evidence snippets
    imaging_context: List         # Imaging findings
    genomic_hit_count: int        # Number of genomic hits
    imaging_hit_count: int        # Number of imaging hits
    enrichment_summary: str       # Human-readable summary
```

#### Threshold Configuration

Cross-modal triggers fire only when variant similarity scores exceed
the configured threshold:

- `CROSS_MODAL_THRESHOLD`: 0.40 (minimum similarity for trigger)
- `CROSS_MODAL_ENABLED`: True (can be disabled for testing)

---

### 3.11 Metrics and Monitoring

**File:** `src/metrics.py` (362 lines)

The metrics module provides Prometheus-compatible instrumentation. It
gracefully degrades to no-op stubs when `prometheus_client` is not
installed.

#### Exposed Metrics

| Metric                      | Type      | Description                    |
|-----------------------------|-----------|--------------------------------|
| onco_agent_up               | Gauge     | Service availability (0/1)     |
| onco_collection_vectors     | Gauge     | Vector count per collection    |
| onco_query_latency          | Histogram | Query processing time          |
| onco_retrieval_latency      | Histogram | Evidence retrieval time        |
| onco_llm_latency            | Histogram | LLM response generation time  |
| onco_embedding_latency      | Histogram | Embedding computation time     |
| onco_milvus_operations      | Counter   | Total Milvus operations        |
| onco_circuit_breaker_state  | Gauge     | Circuit breaker open/closed    |

When `prometheus_client` is not installed, the module creates no-op stubs
that accept the same API but discard data. The `/metrics` endpoint in
`api/main.py` renders metrics in Prometheus text exposition format.

### Scheduler

**File:** `src/scheduler.py` (263 lines)
**Class:** `IngestScheduler`

The ingest scheduler wraps APScheduler's `BackgroundScheduler` to
periodically refresh data from PubMed, ClinicalTrials.gov, and CIViC
into the agent's Milvus collections. Falls back to a no-op stub when
`apscheduler` is not installed.

- Default interval: 168 hours (weekly)
- Configurable via `ONCO_SCHEDULER_INTERVAL`
- Records last ingest timestamp via `LAST_INGEST` metric

---

### 3.12 Ingest Pipeline

**Directory:** `src/ingest/` (1,793 lines across 10 files)

The ingest pipeline consists of 9 specialized parsers plus an abstract
base class. Each parser implements the `fetch -> parse -> embed_and_store`
pattern for a specific data source.

#### Base Class

**File:** `src/ingest/base.py` (249 lines)
**Class:** `BaseIngestPipeline` (ABC)

```
Orchestration pattern:
    run() -> fetch() -> parse() -> embed_and_store()

Parameters:
    collection_manager  Milvus collection manager
    embedder            Embedding model (BGE-small-en-v1.5)
    collection_name     Target Milvus collection
    batch_size          Records per batch (default 50)

Abstract methods (implemented by subclasses):
    fetch()   Retrieve raw data from external source
    parse()   Normalize into a list of record dicts

Concrete methods:
    embed_and_store()  Embed text fields, insert into Milvus in batches
    run()              Orchestrate the full fetch-parse-embed pipeline
```

#### Specialized Parsers

| Parser                      | Lines | Source             | Target Collection   |
|----------------------------|-------|--------------------|---------------------|
| civic_parser.py            | 340   | CIViC API          | onco_variants       |
| oncokb_parser.py           | 104   | OncoKB API         | onco_variants       |
| clinical_trials_parser.py  | 279   | ClinicalTrials.gov | onco_trials         |
| literature_parser.py       | 248   | PubMed/PMC         | onco_literature     |
| guideline_parser.py        | 168   | Curated guidelines | onco_guidelines     |
| pathway_parser.py          | 121   | Curated pathways   | onco_pathways       |
| resistance_parser.py       | 125   | Curated data       | onco_resistance     |
| outcome_parser.py          | 158   | Curated outcomes   | onco_outcomes       |

Each parser normalizes source-specific data formats into the collection
schema and delegates embedding and insertion to the base class.

#### Utility Modules

**`src/utils/pubmed_client.py` (296 lines):**
PubMed E-Utilities client for fetching abstracts and metadata. Supports
API key authentication via `NCBI_API_KEY`. Handles rate limiting and
pagination for bulk retrieval (up to `PUBMED_MAX_RESULTS`, default 5000).

**`src/utils/vcf_parser.py` (371 lines):**
VCF 4.x parser supporting SnpEff ANN, VEP CSQ, and simple GENE/GENEINFO
annotation formats. Extracts gene names, consequence terms, and variant
identifiers from INFO fields using compiled regex patterns.

---

### 3.13 Data Models

**File:** `src/models.py` (538 lines)

Pydantic models defining the agent's type system with 13 enums and 14
domain models.

#### Enums (13 total)

| Enum             | Values | Purpose                                |
|-----------------|--------|----------------------------------------|
| CancerType      | 26     | Supported cancer types (NSCLC..OTHER)  |
| VariantType     | 7      | SNV, INDEL, CNV_AMP, CNV_DEL, FUSION, REARRANGEMENT, SV |
| EvidenceLevel   | 5      | A (FDA), B (clinical), C (case), D (preclinical), E (inferential) |
| TherapyCategory | 9      | TARGETED, IMMUNOTHERAPY, CHEMO, HORMONAL, COMBINATION, RT, CELL_THERAPY, ADC, BISPECIFIC |
| TrialPhase      | 8      | Phase 1 through Phase 4 + sub-phases   |
| TrialStatus     | 9      | Recruiting through Withdrawn           |
| ResponseCategory| 5      | CR, PR, SD, PD, NE                     |
| BiomarkerType   | 8      | TMB, MSI, PD-L1, HRD, etc.            |
| PathwayName     | 13     | MAPK, PI3K, DDR, cell cycle, etc.      |
| GuidelineOrg    | 8      | NCCN, ASCO, ESMO, etc.                |
| SourceType      | 4      | Data provenance categories             |

#### Domain Models (14 total)

Core domain models (Pydantic `BaseModel` subclasses):
- `OncologyVariant` -- Variant with gene, position, evidence level
- `OncologyLiterature` -- Literature chunk with PubMed metadata
- `OncologyTrial` -- Clinical trial with NCT ID, phase, eligibility
- `OncologyTherapy` -- Drug with MOA, indications, adverse events
- `OncologyBiomarker` -- Biomarker with testing method and cutoffs
- `OncologyPathway` -- Signaling pathway with druggable nodes
- `OncologyGuideline` -- Guideline recommendation with organization
- `ResistanceMechanism` -- Resistance entry with bypass strategies
- `OutcomeRecord` -- Treatment outcome with response category
- `CaseSnapshot` -- Patient case with variants and biomarkers
- `MTBPacket` -- Complete MTB review packet

Search and agent I/O models:
- `SearchHit` -- Single search result with score and metadata
- `CrossCollectionResult` -- Aggregated results across collections
- `AgentQuery` / `AgentResponse` -- Agent input/output containers

---

## 4. Data Flow Diagrams

### 4.1 Query Flow (User Question to Answer)

```
User Question: "What are the treatment options for BRAF V600E melanoma?"
    |
    v
[Streamlit UI or API]
    |
    v
[OncoIntelligenceAgent.plan()]
    |   Extracted: gene=BRAF, variant=V600E, cancer_type=MELANOMA
    |   Strategy: targeted
    |   Sub-questions: therapy options, resistance, trials
    |
    v
[OncoRAGEngine.retrieve()]
    |   ThreadPoolExecutor: 11 parallel collection searches
    |   Embedding: BGE-small-en-v1.5 encodes query -> 384-dim vector
    |   Milvus: IVF_FLAT COSINE search with gene="BRAF" filter
    |
    v
[Evidence Merge + Weighting]
    |   onco_variants hits   * 0.18
    |   onco_therapies hits  * 0.14
    |   onco_resistance hits * 0.07
    |   ... (all 11 collections)
    |   Sort by weighted score, take top-k
    |
    v
[OncoIntelligenceAgent.evaluate()]
    |   Result: "sufficient" (12 hits from 5 collections)
    |
    v
[Knowledge Injection]
    |   ACTIONABLE_TARGETS["BRAF"] -> context string
    |   THERAPY_MAP entries for vemurafenib, dabrafenib -> context
    |   RESISTANCE_MAP for BRAF inhibitors -> context
    |
    v
[OncoRAGEngine.synthesize()]
    |   System prompt + knowledge context + evidence + question
    |   -> Claude (claude-sonnet-4-20250514) -> answer with citations
    |
    v
[Formatted Response]
    Cited answer with evidence tiers, therapy options,
    resistance warnings, and trial suggestions
```

### 4.2 Case Ingestion Flow

```
VCF File Upload
    |
    v
[OncologyCaseManager.create_case()]
    |
    v
[vcf_parser.py]
    |   Parse VCF lines
    |   Extract gene annotations (SnpEff ANN / VEP CSQ / GENE tag)
    |   Extract variant identifiers (CHROM, POS, REF, ALT)
    |
    v
[Variant Classification]
    |   For each variant:
    |     classify_variant_actionability(gene, variant)
    |     -> A / B / C / VUS
    |
    v
[CaseSnapshot Creation]
    |   Embed case summary text -> 384-dim vector
    |   Store in onco_cases collection
    |
    v
[MTB Packet Generation]
    |
    +---> TherapyRanker.rank(case) -> ranked therapy list
    +---> TrialMatcher.match(case) -> matched trials
    +---> Cross-modal check (Level A/B variants)
    |
    v
[MTBPacket]
    Actionable variants, ranked therapies, matched trials,
    resistance alerts, biomarker panel, open questions
```

### 4.3 Data Ingest Flow

```
External Sources
    |
    +---> CIViC API ---------> civic_parser.py --------> onco_variants
    +---> OncoKB API --------> oncokb_parser.py ------> onco_variants
    +---> PubMed E-Utils ----> literature_parser.py ---> onco_literature
    +---> ClinicalTrials.gov > clinical_trials_parser -> onco_trials
    +---> Curated JSON ------> guideline_parser.py ----> onco_guidelines
    +---> Curated JSON ------> pathway_parser.py ------> onco_pathways
    +---> Curated JSON ------> resistance_parser.py ---> onco_resistance
    +---> Curated JSON ------> outcome_parser.py ------> onco_outcomes
    |
    v
[BaseIngestPipeline.embed_and_store()]
    |   Batch size: 50 records
    |   Embed text fields with BGE-small-en-v1.5
    |   Insert into Milvus collection
    |
    v
[Milvus Collections]
    11 collections, 384-dim vectors, IVF_FLAT/COSINE
```

### 4.4 Cross-Modal Integration Flow

When Level A/B variants are detected, `OncoCrossModalTrigger.evaluate()`
queries `genomic_evidence` (Stage 1 pipeline) and imaging collections
(if available). Results are combined into a `CrossModalResult` with
enrichment summary for MTB review. Downstream, actionable findings can
trigger the Drug Discovery pipeline (Stage 3, MolMIM/DiffDock).

---

## 5. Collection Architecture

### 5.1 Collection Catalog

The system operates 11 Milvus collections: 10 owned by the oncology agent
and 1 shared read-only collection from the genomics pipeline.

```
+---------------------------+-------------------------------------------+
| Collection                | Purpose                                   |
+---------------------------+-------------------------------------------+
| onco_variants             | Actionable somatic/germline variants      |
|                           | from CIViC and OncoKB                     |
+---------------------------+-------------------------------------------+
| onco_literature           | PubMed/PMC/preprint text chunks tagged    |
|                           | by cancer type and gene                   |
+---------------------------+-------------------------------------------+
| onco_therapies            | Approved and investigational therapies    |
|                           | with mechanism of action                  |
+---------------------------+-------------------------------------------+
| onco_guidelines           | NCCN, ASCO, ESMO guideline recommendations|
+---------------------------+-------------------------------------------+
| onco_trials               | ClinicalTrials.gov summaries with         |
|                           | biomarker eligibility criteria            |
+---------------------------+-------------------------------------------+
| onco_biomarkers           | Predictive and prognostic biomarkers with |
|                           | testing methods and cutoffs               |
+---------------------------+-------------------------------------------+
| onco_resistance           | Resistance mechanisms with bypass         |
|                           | strategies                                |
+---------------------------+-------------------------------------------+
| onco_pathways             | Signaling pathways with druggable nodes   |
|                           | and cross-talk relationships              |
+---------------------------+-------------------------------------------+
| onco_outcomes             | Real-world treatment outcome records      |
+---------------------------+-------------------------------------------+
| onco_cases                | De-identified patient case snapshots      |
+---------------------------+-------------------------------------------+
| genomic_evidence          | Read-only VCF-derived evidence from       |
|                           | Stage 1 genomics pipeline                 |
+---------------------------+-------------------------------------------+
```

### 5.2 Embedding Model

All collections use the same embedding model:

- **Model:** BAAI/bge-small-en-v1.5
- **Dimension:** 384
- **Framework:** SentenceTransformer
- **Batch size:** 32 (configurable via `ONCO_EMBEDDING_BATCH_SIZE`)
- **Wrapper:** `EmbedderWrapper` in `api/main.py` provides both
  `.encode()` (batch) and `.embed()` (single text) interfaces

### 5.3 Index Configuration

All collections use IVF_FLAT with nlist=1024, nprobe=16, COSINE metric.
FLAT stores full vectors (no quantization loss), giving ~98.5% recall
vs brute-force with ~60x query speedup on 10K+ vectors.

### 5.4 Schema Design Principles

1. **Uniform embedding dimension** -- All 384-dim for cross-collection queries.
2. **Filterable metadata** -- `gene` and/or `cancer_type` fields for pre-filtering.
3. **Text chunk field** -- VARCHAR(3000) stores embedded text for display.
4. **Source provenance** -- External IDs (PubMed, NCT, CIViC) for citations.

---

## 6. API Layer

### 6.1 Application Setup

**File:** `api/main.py` (410 lines)

The FastAPI application is configured with:
- **Title:** Oncology Intelligence MTB API
- **Version:** 0.1.0
- **Docs:** Swagger UI at `/docs`, OpenAPI schema at `/openapi.json`
- **Lifespan:** Async context manager initializing 9 services at startup

### 6.2 Middleware

**CORS Middleware:**
```
Origins: http://localhost:8080 (landing page)
         http://localhost:8526 (Streamlit UI)
         http://localhost:8527 (self)
Methods: *
Headers: *
Credentials: allowed
```

**Request Size Limiter:**
Custom HTTP middleware rejects requests exceeding `MAX_REQUEST_SIZE_MB`
(default 10 MB). Checks `Content-Length` header before processing.

### 6.3 Core Endpoints (main.py)

| Method | Path              | Purpose                                 |
|--------|-------------------|-----------------------------------------|
| GET    | /                 | Service info (name, docs link, health)  |
| GET    | /health           | Health check with collection stats       |
| GET    | /collections      | List collections with entity counts      |
| POST   | /query            | Full RAG query (retrieve + LLM answer)  |
| POST   | /search           | Evidence-only vector search (no LLM)    |
| POST   | /find-related     | Cross-collection entity linking          |
| GET    | /knowledge/stats  | Aggregate knowledge-base statistics      |
| GET    | /metrics          | Prometheus text exposition format        |

### 6.4 Route Files

**meta_agent.py (169 lines):**
Intelligence agent endpoints for clinical Q&A. Accepts natural-language
oncology questions with optional cancer type and gene filters. Returns
cited answers with evidence tiers and search plan metadata.

**cases.py (238 lines):**
Case lifecycle endpoints:
- Create case from VCF text or pre-parsed variants
- Retrieve case by ID
- Generate MTB packet for a case
- List recent cases

**trials.py (156 lines):**
Clinical trial matching endpoints:
- Match patient profile to trials
- Search trials by cancer type and biomarkers
- Retrieve trial details by NCT ID

**reports.py (236 lines):**
Report generation endpoints:
- Generate Markdown report
- Generate JSON report
- Generate PDF report (binary download)
- Generate FHIR R4 Bundle

**events.py (89 lines):**
Cross-modal event endpoints:
- Receive events from other agents
- Publish variant discovery events
- Query event history

### 6.5 Health Check

The `/health` endpoint returns status ("healthy" or "degraded"), per-
collection vector counts, total vectors, version, and boolean status for
all 7 services (milvus, embedder, rag_engine, intelligence_agent,
case_manager, trial_matcher, therapy_ranker). Status is "healthy" only
when all 7 service checks pass.

---

## 7. UI Architecture

**File:** `app/oncology_ui.py` (758 lines)
**Framework:** Streamlit
**Port:** 8526

The Streamlit UI provides a Molecular Tumor Board workbench with five
tabs, communicating with the FastAPI backend over HTTP.

### 7.1 Five-Tab Layout

```
+-------+----------+-------+----------+----------+
| Case  | Evidence | Trial | Therapy  | Outcomes |
| Work  | Explorer | Finder| Ranker   | Dashboard|
| bench |          |       |          |          |
+-------+----------+-------+----------+----------+
```

**Tab 1: Case Workbench**
- VCF file upload or text paste
- Cancer type selection (20 types + Other)
- Stage selection (I through IVB)
- Biomarker checkboxes (15 options: EGFR+, ALK+, BRAF V600E, etc.)
- Prior therapy selection
- Case creation triggers MTB packet generation

**Tab 2: Evidence Explorer**
- Free-text oncology question input
- Gene and cancer type filters
- RAG-powered evidence retrieval
- Results with collection labels and relevance scores

**Tab 3: Trial Finder**
- Patient profile input (cancer type, biomarkers, stage, age)
- Hybrid deterministic + semantic trial matching
- Results ranked by composite score
- NCT ID links to ClinicalTrials.gov

**Tab 4: Therapy Ranker**
- Variant input (gene, variant, cancer type)
- Prior therapy history for resistance checking
- Evidence-ranked therapy recommendations
- Resistance warnings and contraindications

**Tab 5: Outcomes Dashboard**
- Treatment outcome visualization
- Response category distribution
- Outcome records by cancer type and therapy

### 7.2 Configuration

```python
API_BASE = os.environ.get("ONCO_API_BASE_URL", "http://localhost:8527")
PAGE_TITLE = "Oncology Intelligence MTB Workbench"
```

The UI provides 20 cancer type options, 13 stage values (I through IVB),
and 15 biomarker checkboxes (EGFR+, ALK+, ROS1+, BRAF V600E, KRAS G12C,
MSI-H, TMB-H, PD-L1>=50%, HER2+, BRCA+, NTRK fusion, RET fusion,
MET amplification, PIK3CA mutation, FGFR alteration).

---

## 8. Docker and Deployment

### 8.1 Docker Compose Architecture

**File:** `docker-compose.yml` (209 lines)

Six services on the `onco-network` bridge network:

```
+----------------------------------------------------------+
|  onco-network (bridge)                                    |
|                                                           |
|  +-------------+  +-------------+  +------------------+  |
|  | milvus-etcd |  | milvus-minio|  | milvus-standalone|  |
|  | (etcd v3.5) |  | (MinIO)     |  | (Milvus 2.4)    |  |
|  | metadata    |  | object store|  | :19530 (gRPC)    |  |
|  | store       |  |             |  | :9091 (metrics)  |  |
|  +------+------+  +------+------+  +--------+---------+  |
|         |                |                   |            |
|         +-------+--------+-------------------+            |
|                 |                                         |
|  +--------------v-+  +----------------+  +-------------+ |
|  | onco-streamlit |  | onco-api       |  | onco-setup  | |
|  | :8526          |  | :8527          |  | (one-shot)  | |
|  | Streamlit UI   |  | FastAPI + uvi  |  | collections | |
|  |                |  | corn (2 wkrs)  |  | + seed data | |
|  +----------------+  +----------------+  +-------------+ |
+----------------------------------------------------------+
```

### 8.2 Dockerfile

**File:** `Dockerfile` (42 lines)
**Base image:** python:3.10-slim
**Build strategy:** Multi-stage (builder + runtime)

Stage 1 (Builder) installs build deps and Python packages into a venv.
Stage 2 (Runtime) copies the venv, installs runtime libraries, copies
application code (config/, src/, api/, app/, scripts/, data/), creates
non-root user `oncouser`, exposes ports 8526+8527, and defaults to
Streamlit as CMD.

### 8.3 Service Dependencies

```
milvus-etcd ----+
                +--> milvus-standalone --+--> onco-streamlit
milvus-minio ---+                       +--> onco-api
                                        +--> onco-setup
```

All three application services depend on `milvus-standalone` reaching
`service_healthy`. Milvus depends on both etcd and MinIO being healthy.

### 8.4 Health Checks

All services use 30-second interval health checks. Infrastructure services
(etcd, MinIO, Milvus) use `etcdctl`, MinIO health API, and Milvus healthz
endpoints respectively. Application services use their own HTTP health
endpoints. Milvus has a 60-second `start_period` for index loading.

### 8.5 Setup Service (One-Shot)

The `onco-setup` container runs once (`restart: "no"`) and executes:

```bash
1. python scripts/setup_collections.py --drop-existing
2. python scripts/seed_variants.py
3. python scripts/seed_literature.py
4. python scripts/seed_trials.py
5. python scripts/seed_therapies.py
6. python scripts/seed_biomarkers.py
7. python scripts/seed_pathways.py
8. python scripts/seed_guidelines.py
9. python scripts/seed_resistance.py
10. python scripts/seed_outcomes.py
```

### 8.6 Configuration

All settings are configurable via environment variables with the `ONCO_`
prefix (Pydantic BaseSettings with `env_prefix="ONCO_"`):

```
ONCO_MILVUS_HOST          default: localhost
ONCO_MILVUS_PORT          default: 19530
ONCO_EMBEDDING_MODEL      default: BAAI/bge-small-en-v1.5
ONCO_LLM_MODEL            default: claude-sonnet-4-20250514
ONCO_API_PORT              default: 8527
ONCO_STREAMLIT_PORT        default: 8526
ONCO_TOP_K                 default: 5
ONCO_SCORE_THRESHOLD       default: 0.4
ANTHROPIC_API_KEY          required for LLM calls
```

The `ANTHROPIC_API_KEY` is passed through from the host `.env` file via
Docker Compose environment variable interpolation (`${ANTHROPIC_API_KEY}`).

---

## 9. Testing Architecture

### 9.1 Overview

The test suite contains 10 test files plus a shared `conftest.py`,
totaling 4,584 lines and 556 test cases (including parametrized expansions).

**Framework:** pytest
**Mocking:** unittest.mock (MagicMock)
**Strategy:** Unit tests with mocked infrastructure (no live Milvus
required for standard test runs)

### 9.2 Test File Inventory

| File                  | Lines | Focus Area                        |
|-----------------------|-------|-----------------------------------|
| conftest.py           | 214   | Shared fixtures (embedder, LLM, Milvus mocks) |
| test_agent.py         | 264   | Intelligence agent plan/search/evaluate |
| test_case_manager.py  | 332   | VCF parsing, case creation, MTB packets |
| test_collections.py   | 276   | Collection CRUD, schema validation |
| test_export.py        | 439   | Markdown, JSON, PDF, FHIR R4 export |
| test_integration.py   | 785   | End-to-end flows with mocked services |
| test_knowledge.py     | 603   | Knowledge graph completeness and helpers |
| test_models.py        | 644   | Pydantic model validation, enum coverage |
| test_rag_engine.py    | 301   | RAG retrieval, weighting, citation formatting |
| test_therapy_ranker.py| 363   | Therapy ranking, resistance, contraindications |
| test_trial_matcher.py | 363   | Trial matching, scoring, phase weights |

### 9.3 Shared Fixtures (conftest.py)

The conftest provides mock objects that simulate infrastructure without
requiring live services:

- **`mock_embedder`:** Returns 384-dim zero vectors. Supports `.embed_text()`,
  `.encode()`, and `.embed()` interfaces.
- **`mock_llm_client`:** Returns "Mock response" for all LLM calls.
- **`mock_collection_manager`:** Simulates all 11 collections with empty
  search results and 42-entity stats per collection.
- **`ALL_COLLECTION_NAMES`:** Constant list of all 11 collection names.

### 9.4 Test Categories

- **Unit Tests** (test_agent, test_rag_engine, test_therapy_ranker,
  test_trial_matcher): Individual component logic with mocked dependencies.
- **Model Tests** (test_models, test_knowledge): Data model completeness,
  enum validity, Pydantic serialization, knowledge graph field coverage.
- **Export Tests** (test_export): All four export formats produce valid
  output with correct headers, codes (LOINC, SNOMED), and structure.
- **Integration Tests** (test_integration): End-to-end flows with all
  components wired together (Milvus and LLM mocked).

### 9.5 Scripts for Pipeline Testing

Beyond pytest, two scripts provide end-to-end validation:

**`scripts/test_rag_pipeline.py` (200 lines):**
Sends sample queries through the full RAG pipeline against a live Milvus
instance. Measures retrieval latency, LLM response time, and citation
quality.

**`scripts/validate_e2e.py` (329 lines):**
Comprehensive end-to-end validation script. Creates a test case from
sample VCF data, generates an MTB packet, matches trials, ranks therapies,
exports reports in all four formats, and validates each output.

---

## 10. Security Considerations

### 10.1 Input Validation

**Filter value sanitization:**
Both `case_manager.py` and `trial_matcher.py` validate filter values
against a regex whitelist before constructing Milvus expressions:

```python
_SAFE_FILTER_RE = re.compile(r"^[A-Za-z0-9 _.\-/]+$")
```

This prevents Milvus expression injection by rejecting filter values
containing special characters like quotes, semicolons, or boolean operators.

**Request size limiting:**
The FastAPI middleware rejects HTTP requests exceeding `MAX_REQUEST_SIZE_MB`
(default 10 MB) based on the `Content-Length` header.

### 10.2 Authentication and Authorization

No authentication middleware is included. The API is designed for internal
network use behind a reverse proxy. CORS restricts browser requests to
known origins. For production, add API key or OAuth2 authentication.

### 10.3 Secrets Management

`ANTHROPIC_API_KEY` and `NCBI_API_KEY` are loaded from environment variables
(`.env` file or Docker Compose interpolation). Both are Optional -- the
agent starts without them, but LLM and PubMed calls fail. Ensure `.env`
is in `.gitignore` to prevent key exposure.

### 10.4 Container Security

The Dockerfile runs as non-root user `oncouser`. Milvus uses
`seccomp:unconfined` (required for mmap, relaxes isolation). MinIO default
credentials should be changed for production.

### 10.5 Data Privacy

Patient data in `onco_cases` should be de-identified before ingestion.
FHIR R4 exports use anonymized patient references. No PHI is stored in
vector embeddings.

---

## 11. Performance Characteristics

### 11.1 Query Latency Budget

```
End-to-end query processing (typical):

Embedding:        ~20 ms   (BGE-small-en-v1.5, single query)
Milvus search:    ~50 ms   (11 collections in parallel, IVF_FLAT)
Knowledge inject: ~2 ms    (Python dict lookup)
LLM generation:   ~2-5 s   (Claude, depends on response length)
Citation format:  ~5 ms    (string formatting)
                  -------
Total:            ~2-5 s   (dominated by LLM latency)
```

### 11.2 Throughput Considerations

- FastAPI runs with 2 Uvicorn workers (`onco-api` Docker service)
- `ThreadPoolExecutor` parallelizes 11 collection searches per query
- Embedding model is CPU-bound (BGE-small-en-v1.5 is lightweight)
- Milvus standalone supports concurrent queries via gRPC

### 11.3 Memory Footprint

- BGE-small-en-v1.5 model: ~130 MB resident
- IVF_FLAT index for 10K vectors at 384-dim: ~15 MB per collection
- Milvus standalone: ~2-4 GB resident (depends on loaded collections)
- Python application: ~200-400 MB (includes model + dependencies)

### 11.4 Embedding Batch Performance

```
BGE-small-en-v1.5 throughput (CPU):
  Single text:     ~20 ms
  Batch of 32:     ~200 ms  (~6 ms per text)
  Batch of 100:    ~550 ms  (~5.5 ms per text)
```

Batch embedding during data ingest uses `EMBEDDING_BATCH_SIZE=32` for
a good balance between memory usage and throughput.

### 11.5 Collection Weights and Relevance

The default weight distribution is tuned for clinical MTB decision support:
- Variants (0.18) and literature (0.16) receive the highest weights
  because they are the primary evidence sources for variant interpretation.
- Therapies (0.14) and guidelines (0.12) are weighted next because
  treatment decisions depend heavily on approved therapies and guideline
  concordance.
- Trials (0.10) are moderately weighted -- important for investigational
  options but secondary to established evidence.
- Biomarkers (0.08), resistance (0.07), and pathways (0.06) provide
  supporting context.
- Outcomes (0.04), genomic (0.03), and cases (0.02) are low-weighted
  supplementary sources.

All weights are configurable via environment variables
(`ONCO_WEIGHT_VARIANTS`, `ONCO_WEIGHT_LITERATURE`, etc.).

---

## 12. Extension Points

### 12.1 Adding a New Collection

To add a 12th collection (e.g., `onco_immunoprofiles`):

1. **Define schema** in `src/collections.py`:
   Add a `ONCO_IMMUNOPROFILES_FIELDS` list following the existing pattern.
   Include `id`, `embedding`, `text_chunk`, and domain-specific metadata.

2. **Add config** in `config/settings.py`:
   Add `COLLECTION_IMMUNOPROFILES: str = "onco_immunoprofiles"` and a
   corresponding weight `WEIGHT_IMMUNOPROFILES: float = 0.05` (adjust
   other weights to maintain sum ~1.0).

3. **Register in RAG engine** (`src/rag_engine.py`):
   Add an entry to `COLLECTION_CONFIG` with the weight, label, and
   filter fields.

4. **Create ingest parser** in `src/ingest/`:
   Subclass `BaseIngestPipeline` and implement `fetch()` and `parse()`.

5. **Add seed script** in `scripts/`:
   Create `seed_immunoprofiles.py` following the existing seed patterns.

6. **Update setup** in `scripts/setup_collections.py` and the Docker
   Compose `onco-setup` command.

7. **Update conftest.py**:
   Add the new collection name to `ALL_COLLECTION_NAMES`.

### 12.2 Adding a New Data Source

To add a new external data source (e.g., COSMIC):

1. **Create parser** in `src/ingest/cosmic_parser.py`:
   Subclass `BaseIngestPipeline`. Implement `fetch()` to pull data
   from the COSMIC API and `parse()` to normalize records.

2. **Create ingest script** in `scripts/ingest_cosmic.py`.

3. **Register with scheduler** if automated refresh is desired.

### 12.3 Adding a New Export Format

To add a new export format (e.g., CSV):

1. Add `export_csv()` function in `src/export.py`.
2. Add route in `api/routes/reports.py`.
3. Add test in `tests/test_export.py`.

### 12.4 Adding a New API Route

To add new endpoints:

1. Create a new route file in `api/routes/` (e.g., `biomarkers.py`).
2. Define an `APIRouter` with a tag and prefix.
3. Register in `api/main.py` with `app.include_router()`.

### 12.5 Cross-Modal Integration with New Agents

To integrate with a new HCLS AI Factory agent:

1. Add event handlers in `api/routes/events.py`.
2. Extend `OncoCrossModalTrigger` to query the new agent's collections.
3. Define trigger conditions (which evidence levels, which variant types).

### 12.6 Custom Knowledge Graph Entries

Add entries to `ACTIONABLE_TARGETS` in `src/knowledge.py` following the
existing structure (gene, full_name, cancer_types, key_variants,
targeted_therapies, combination_therapies, resistance_mutations, pathway,
evidence_level, description). Update `THERAPY_MAP`, `RESISTANCE_MAP`,
`BIOMARKER_PANELS`, and `ENTITY_ALIASES` as needed.

### 12.7 LLM Provider Swap

The LLM provider is configurable via `ONCO_LLM_PROVIDER` (default:
"anthropic") and `ONCO_LLM_MODEL` (default: claude-sonnet-4-20250514). To add a
new provider, modify the LLM client instantiation in the RAG engine's
`synthesize()` method.

### 12.8 Embedding Model Swap

Update `ONCO_EMBEDDING_MODEL` and `ONCO_EMBEDDING_DIM` in settings. This
requires re-creating all collections (new dimension = schema change) and
re-running all seed scripts. There is no migration path for existing vectors.

---

## Appendix: Port Map

| Port  | Service            | Protocol |
|-------|--------------------|----------|
| 2379  | etcd (internal)    | gRPC     |
| 9000  | MinIO API          | HTTP     |
| 9001  | MinIO Console      | HTTP     |
| 9091  | Milvus metrics     | HTTP     |
| 19530 | Milvus gRPC        | gRPC     |
| 8526  | Streamlit UI       | HTTP     |
| 8527  | FastAPI API        | HTTP     |

All settings are configurable via environment variables with the `ONCO_`
prefix. See `config/settings.py` (134 lines) for the complete list of
46 configuration parameters with defaults and descriptions.

---

*This document reflects the codebase as of March 2026, version 1.3.0.*
