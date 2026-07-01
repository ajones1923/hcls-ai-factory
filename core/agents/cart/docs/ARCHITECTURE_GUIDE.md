# CAR-T Intelligence Agent -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0
**Codebase:** `hcls-ai-factory/ai_agent_adds/cart_intelligence_agent/`

---

## Table of Contents

1. [Overview](#1-overview)
2. [High-Level Architecture](#2-high-level-architecture)
3. [Component Deep Dives](#3-component-deep-dives)
   - 3a. [RAG Engine](#3a-rag-engine)
   - 3b. [Collection Manager](#3b-collection-manager)
   - 3c. [Knowledge Graph](#3c-knowledge-graph)
   - 3d. [Query Expansion](#3d-query-expansion)
   - 3e. [Agent](#3e-agent)
   - 3f. [Export System](#3f-export-system)
   - 3g. [Models](#3g-models)
   - 3h. [Ingest Pipeline](#3h-ingest-pipeline)
   - 3i. [Metrics](#3i-metrics)
   - 3j. [Scheduler](#3j-scheduler)
4. [Data Flow Diagrams](#4-data-flow-diagrams)
5. [API Layer](#5-api-layer)
6. [UI Architecture](#6-ui-architecture)
7. [Docker and Deployment](#7-docker-and-deployment)
8. [Testing Architecture](#8-testing-architecture)
9. [Security Considerations](#9-security-considerations)
10. [Performance Characteristics](#10-performance-characteristics)
11. [Extension Points](#11-extension-points)

---

## 1. Overview

The CAR-T Intelligence Agent is part of the HCLS AI Factory Precision Intelligence Network -- one of three GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence Network, Therapeutic Discovery Engine) that together deliver patient DNA to drug candidates in under 5 hours on a single NVIDIA DGX Spark.

The agent is a multi-collection Retrieval-Augmented Generation (RAG) system purpose-built for CAR-T cell therapy research and development. It breaks down data silos across the entire CAR-T development lifecycle -- from target identification and molecular design through clinical development, manufacturing, regulatory approval, and post-market pharmacovigilance -- and provides unified, cross-functional intelligence grounded in evidence.

### What It Does

A researcher asks a natural-language question such as *"Why do CD19 CAR-T patients relapse with antigen-negative disease?"* The system:

1. Expands the query across 12 domain-specific terminology maps (229 keyword entries)
2. Embeds the question using BGE-small-en-v1.5 (384 dimensions)
3. Searches 11 Milvus vector collections simultaneously using parallel ThreadPoolExecutor
4. Applies collection-specific score weights to rank cross-domain evidence
5. Augments with structured knowledge from a 3-dictionary knowledge graph
6. Synthesizes a grounded response via Claude Sonnet 4.6 with clickable citations
7. Returns the answer with a scored evidence panel and export options (Markdown, JSON, PDF)

### Position in the HCLS AI Factory

The CAR-T Intelligence Agent extends the three-stage HCLS AI Factory pipeline (Genomics, RAG/Chat, Drug Discovery) with a domain-specialized fourth capability. It reads from the existing `genomic_evidence` collection populated by Stage 2 (rag-chat-pipeline) and adds 10 new CAR-T-specific collections.

```
HCLS AI Factory Pipeline
========================

Stage 1: Genomics       FASTQ -> VCF (Parabricks)
Stage 2: RAG/Chat       VCF -> Targets (Milvus + Claude)
Stage 3: Drug Discovery  Target -> Molecules (BioNeMo)

         +-- CAR-T Intelligence Agent --+
         |  11 collections (10 new + 1  |
         |  shared genomic_evidence)    |
         |  Knowledge graph + expansion |
         |  Autonomous agent reasoning  |
         +------------------------------+
```

### Key Numbers

| Metric | Value |
|--------|-------|
| Source files | ~30 Python modules |
| Total lines of code | ~21,259 |
| Milvus collections | 11 (10 owned + 1 read-only) |
| Knowledge graph entries | ~34 targets + 17 toxicities + 20 manufacturing + 23 biomarkers + 6 regulatory + 6 immunogenicity |
| Query expansion maps | 12 maps, 229 keyword entries |
| Seed records | 649 across 13 JSON files |
| Test cases | 415 across 7 test files |
| Embedding model | BGE-small-en-v1.5 (384-dim, COSINE) |
| LLM | Claude Sonnet 4.6 (Anthropic) |
| Docker services | 6 (Milvus stack + UI + API + setup) |

---

## 2. High-Level Architecture

```
                                    CAR-T Intelligence Agent
                                    =======================

    +------------------+     +------------------+
    |   Streamlit UI   |     |   FastAPI REST   |
    |   (port 8521)    |     |   (port 8522)    |
    +--------+---------+     +--------+---------+
             |                        |
             +----------+-------------+
                        |
                        v
              +---------+----------+
              |   CARTRAGEngine    |    <-- Core orchestrator
              |   (rag_engine.py)  |
              +---------+----------+
                        |
         +--------------+--------------+--------------+
         |              |              |              |
         v              v              v              v
  +------+------+ +-----+-----+ +-----+------+ +----+--------+
  | Collection  | | Knowledge | | Query      | | LLM Client  |
  | Manager     | | Graph     | | Expansion  | | (Claude)    |
  | (Milvus)    | | (3 dicts) | | (12 maps)  | |             |
  +------+------+ +-----------+ +------------+ +-------------+
         |
         v
  +------+-------------------------------------------------------+
  | Milvus 2.4 Standalone                                        |
  |                                                              |
  | cart_literature    cart_trials       cart_constructs          |
  | cart_assays        cart_manufacturing cart_safety             |
  | cart_biomarkers    cart_regulatory   cart_sequences           |
  | cart_realworld     genomic_evidence (read-only, shared)      |
  +--------------------------------------------------------------+
         |                        |
         v                        v
  +------+------+          +------+------+
  | etcd v3.5.5 |          | MinIO       |
  | (metadata)  |          | (object     |
  |             |          |  storage)   |
  +-------------+          +-------------+
```

### Request Lifecycle

```
User Question
     |
     v
[1] AgentQuery model (question, target_antigen, cart_stage, include_genomic)
     |
     v
[2] Query Expansion: scan 12 maps for keyword matches -> expansion terms
     |
     v
[3] Embed: BGE-small-en-v1.5 with "Represent this sentence..." prefix
     |
     v
[4] Parallel Search: ThreadPoolExecutor -> search_all() across 11 collections
     |
     v
[5] Score Weighting: raw_score * (1 + collection_weight)
     |
     v
[6] Merge & Rank: deduplicate by ID, sort by weighted score, cap at 30
     |
     v
[7] Knowledge Augmentation: match query against 6 knowledge dictionaries
     |
     v
[8] Build Prompt: evidence snippets + knowledge context + system prompt
     |
     v
[9] LLM Synthesis: Claude Sonnet 4.6 generates response with citations
     |
     v
[10] Citation Scoring: High (>=0.75), Medium (>=0.60), Low (<0.60)
     |
     v
[11] Response with evidence panel -> UI or API client
```

---

## 3. Component Deep Dives

### 3a. RAG Engine

**File:** `src/rag_engine.py` (754 lines)
**Class:** `CARTRAGEngine`

The RAG engine is the central orchestrator. It wires together the collection manager, embedder, LLM client, knowledge graph, and query expander into a unified retrieval-and-generation pipeline.

#### Constructor

```python
class CARTRAGEngine:
    def __init__(self, collection_manager, embedder, llm_client,
                 knowledge=None, query_expander=None):
        self.collections = collection_manager
        self.embedder = embedder
        self.llm = llm_client
        self.knowledge = knowledge
        self.expander = query_expander
```

All dependencies are injected, making the engine fully testable with mocks.

#### COLLECTION_CONFIG

Defines the weight, label, target-antigen capability, and year-field availability for each of the 11 collections. Weights are read from `CARTSettings` at import time:

```python
COLLECTION_CONFIG = {
    "cart_literature":    {"weight": 0.20, "label": "Literature",    "has_target_antigen": True,  "year_field": "year"},
    "cart_trials":        {"weight": 0.16, "label": "Trial",         "has_target_antigen": True,  "year_field": "start_year"},
    "cart_constructs":    {"weight": 0.10, "label": "Construct",     "has_target_antigen": True,  "year_field": None},
    "cart_assays":        {"weight": 0.09, "label": "Assay",         "has_target_antigen": True,  "year_field": None},
    "cart_safety":        {"weight": 0.08, "label": "Safety",        "has_target_antigen": False, "year_field": "year"},
    "cart_biomarkers":    {"weight": 0.08, "label": "Biomarker",     "has_target_antigen": True,  "year_field": None},
    "cart_manufacturing": {"weight": 0.07, "label": "Manufacturing", "has_target_antigen": False, "year_field": None},
    "cart_realworld":     {"weight": 0.07, "label": "RealWorld",     "has_target_antigen": False, "year_field": None},
    "cart_regulatory":    {"weight": 0.06, "label": "Regulatory",    "has_target_antigen": False, "year_field": None},
    "cart_sequences":     {"weight": 0.06, "label": "Sequence",      "has_target_antigen": True,  "year_field": None},
    "genomic_evidence":   {"weight": 0.04, "label": "Genomic",       "has_target_antigen": False, "year_field": None},
}
```

Weights sum to approximately 1.01. Literature and trials receive the heaviest weighting because they contain the densest evidence text.

#### Core Methods

| Method | Purpose | Returns |
|--------|---------|---------|
| `retrieve()` | Full evidence retrieval pipeline (embed, search, expand, merge, augment) | `CrossCollectionResult` |
| `query()` | retrieve() + LLM synthesis (blocking) | `str` (answer text) |
| `query_stream()` | retrieve() + streaming LLM (yields evidence event, then token chunks, then done event) | `Generator[Dict]` |
| `find_related()` | Cross-collection entity linking ("show me everything about Yescarta") | `Dict[str, List[SearchHit]]` |
| `retrieve_comparative()` | Two-entity comparison with separate evidence sets | `ComparativeResult` |

#### CART_SYSTEM_PROMPT

A 64-line system prompt that establishes the agent's 12 domains of expertise, citation formatting rules (clickable PubMed/ClinicalTrials.gov links), cross-functional reasoning guidelines, and the directive to acknowledge uncertainty.

#### Score Weighting Formula

```
weighted_score = raw_cosine_score * (1 + collection_weight)
```

For a literature hit with raw score 0.85: `0.85 * (1 + 0.20) = 1.02`. For a genomic hit with raw score 0.85: `0.85 * (1 + 0.04) = 0.884`. This ensures higher-weighted collections float to the top when scores are comparable.

#### Query Expansion Strategy

The engine handles expanded terms differently based on whether they are known target antigens:

- **Antigen terms** (e.g., "CD19"): Used as Milvus field filters (`target_antigen == "CD19"`) on collections that have the `target_antigen` field. Scores are discounted by 0.8x.
- **Non-antigen terms** (e.g., "cytokine release syndrome"): Re-embedded as new query vectors and searched semantically across all collections. Scores are discounted by 0.7x.

The discount factors prevent expansion hits from dominating the primary search results.

#### Comparative Analysis

When the engine detects comparison patterns ("compare", "vs", "versus") in a question, it:

1. Parses two entities from the question text using regex
2. Resolves each entity through the knowledge graph's `ENTITY_ALIASES` dictionary
3. Runs two separate `retrieve()` calls, one per entity
4. Builds a structured comparison context from the knowledge graph
5. Generates a response requiring: comparison table, advantages/limitations lists, and clinical context

---

### 3b. Collection Manager

**File:** `src/collections.py`
**Class:** `CARTCollectionManager`

Manages 11 Milvus collections (10 owned by this agent + 1 read-only `genomic_evidence` from rag-chat-pipeline).

#### Schema Summary

All collections share:
- **384-dim FLOAT_VECTOR** `embedding` field (BGE-small-en-v1.5)
- **IVF_FLAT** index with `nlist=1024`
- **COSINE** similarity metric
- **Search params:** `nprobe=16`

| Collection | Primary Key | Key Fields | VARCHAR Limits |
|------------|-------------|------------|----------------|
| `cart_literature` | id (PMID/patent) | title, text_chunk, source_type, year, cart_stage, target_antigen, disease, keywords, journal | text_chunk: 3000, title: 500, keywords: 1000 |
| `cart_trials` | id (NCT number) | title, text_summary, phase, status, sponsor, target_antigen, car_generation, costimulatory, disease, enrollment, start_year, outcome_summary | text_summary: 3000, outcome_summary: 2000 |
| `cart_constructs` | id | name, text_summary, target_antigen, scfv_origin, costimulatory_domain, signaling_domain, generation, hinge_tm, vector_type, fda_status, known_toxicities | text_summary: 2000, known_toxicities: 500 |
| `cart_assays` | id | text_summary, assay_type, construct_id, target_antigen, cell_line, effector_ratio, key_metric, metric_value (FLOAT), outcome, notes | text_summary: 2000, notes: 1000 |
| `cart_manufacturing` | id | text_summary, process_step, vector_type, parameter, parameter_value, target_spec, met_spec, batch_id, notes | text_summary: 2000, notes: 1000 |
| `cart_safety` | id | text_summary, product, event_type, severity_grade, onset_timing, incidence_rate, management_protocol, outcome, reporting_source, year | text_summary: 3000, management_protocol: 500 |
| `cart_biomarkers` | id | text_summary, biomarker_name, biomarker_type, assay_method, clinical_cutoff, predictive_value, associated_outcome, target_antigen, disease, evidence_level | text_summary: 3000, predictive_value: 200 |
| `cart_regulatory` | id | text_summary, product, regulatory_event, date, agency, indication, decision, conditions, pivotal_trial | text_summary: 3000, conditions: 500 |
| `cart_sequences` | id | text_summary, construct_name, target_antigen, scfv_clone, binding_affinity_kd, heavy_chain_vregion, light_chain_vregion, framework, species_origin, immunogenicity_risk, structural_notes | heavy/light_chain: 500, structural_notes: 1000 |
| `cart_realworld` | id | text_summary, study_type, data_source, product, indication, population_size (INT64), median_followup_months (FLOAT), primary_endpoint, outcome_value, setting, special_population | text_summary: 3000, special_population: 200 |
| `genomic_evidence` | id | chrom, pos, ref, alt, qual, gene, consequence, impact, genotype, text_summary, clinical_significance, rsid, disease_associations, am_pathogenicity, am_class | **Read-only** -- created by rag-chat-pipeline |

#### Parallel Search Architecture

```python
def search_all(self, query_embedding, top_k_per_collection=5,
               filter_exprs=None, score_threshold=0.0):
    collections = list(COLLECTION_SCHEMAS.keys())

    def _search_one(name):
        expr = (filter_exprs or {}).get(name)
        return name, self.search(name, query_embedding, top_k, expr, score_threshold)

    with ThreadPoolExecutor(max_workers=len(collections)) as executor:
        futures = {executor.submit(_search_one, name): name for name in collections}
        for future in as_completed(futures):
            name, hits = future.result()
            all_results[name] = hits

    return all_results
```

With 11 collections, this spawns 11 concurrent threads. Each thread independently loads its collection into memory, executes the ANN search, and returns results. The `as_completed` pattern collects results as they arrive, without waiting for the slowest collection.

#### Registry Pattern

Two dictionaries provide the single source of truth for collection metadata:

```python
COLLECTION_SCHEMAS: Dict[str, CollectionSchema]  # Schema definitions for creation
COLLECTION_MODELS: Dict[str, type]               # Pydantic models for validation
```

`genomic_evidence` has `None` in `COLLECTION_MODELS` because it is read-only -- this agent never inserts into it.

---

### 3c. Knowledge Graph

**File:** `src/knowledge.py`
**Size:** 3 dictionaries, ~2,249 lines

The knowledge graph provides structured, curated domain knowledge that supplements the vector search results. Unlike the Milvus collections which contain embedded free-text, the knowledge graph contains organized factual data.

#### Three Dictionaries

| Dictionary | Entries | Content |
|-----------|---------|---------|
| `CART_TARGETS` | 34 | Target antigens with protein name, UniProt ID, expression profile, diseases, approved products, key trials, resistance mechanisms, toxicity profile, normal tissue expression |
| `CART_TOXICITIES` | 17 | Toxicity profiles with mechanism, ASTCT grading, incidence, timing, management protocols, biomarkers, risk factors |
| `CART_MANUFACTURING` | 20 | Manufacturing processes with parameters, specifications, failure modes, release criteria |
| `CART_BIOMARKERS` | 23 | Biomarkers with type (predictive/prognostic/PD/monitoring/resistance), assay method, clinical cutoff, predictive value, evidence level, PMID references |
| `CART_REGULATORY` | 6 | FDA-approved products with approval dates, indications, pivotal trials, designations (BTD/RMAT), REMS, EMA approval, subsequent expansions |
| `CART_IMMUNOGENICITY` | 6 | HLA-restricted epitopes, humanization strategies, ADA clinical impact, allogeneic HLA considerations, immunogenicity testing paradigms |

#### Context Retrieval

The public API surfaces five retrieval functions that scan query text for entity mentions:

```python
get_target_context(antigen)        # "CD19" -> formatted target knowledge
get_toxicity_context(toxicity)     # "CRS"  -> grading, management, biomarkers
get_manufacturing_context(process) # "lentiviral_transduction" -> parameters
get_biomarker_context(biomarker)   # "ferritin" -> cutoffs, predictive value
get_regulatory_context(product)    # "Kymriah" -> approval history
get_immunogenicity_context(topic)  # "humanization" -> strategies, tools
```

The master function `get_all_context_for_query()` scans the query against all three dictionaries using keyword matching, returning combined context.

#### Entity Resolution for Comparisons

The `ENTITY_ALIASES` dictionary (54 entries) maps product names, generic names, costimulatory domains, vector types, biomarker names, and immunogenicity terms to canonical entities. This powers the comparative analysis pipeline:

```python
ENTITY_ALIASES = {
    "KYMRIAH":          {"type": "product", "canonical": "Kymriah (tisagenlecleucel)", "target": "CD19"},
    "TISAGENLECLEUCEL": {"type": "product", "canonical": "Kymriah (tisagenlecleucel)", "target": "CD19"},
    "4-1BB":            {"type": "costimulatory", "canonical": "4-1BB (CD137)"},
    "FERRITIN":         {"type": "biomarker", "canonical": "ferritin"},
    ...
}
```

`resolve_comparison_entity()` tries five resolution strategies in priority order: target antigens (exact), product aliases, toxicities, manufacturing processes, and biomarkers.

---

### 3d. Query Expansion

**File:** `src/query_expansion.py`
**Size:** 12 maps, 229 keyword entries, 1,592 lines

Query expansion improves recall by broadening search terms. When a user asks about "CRS," the system also searches for "cytokine release syndrome," "tocilizumab," "IL-6," "ferritin," "Lee grading," and 25+ related terms.

#### Twelve Expansion Maps

| Map | Keywords | Purpose |
|-----|----------|---------|
| `TARGET_ANTIGEN_EXPANSION` | 32 | Antigen -> diseases, products, clones |
| `DISEASE_EXPANSION` | 16 | Disease -> antigens, products, subtypes |
| `TOXICITY_EXPANSION` | 18 | Toxicity -> mechanisms, treatments, biomarkers |
| `MANUFACTURING_EXPANSION` | 21 | CMC -> platforms, parameters, reagents |
| `MECHANISM_EXPANSION` | 19 | Biology -> pathways, markers, interventions |
| `CONSTRUCT_EXPANSION` | 20 | Engineering -> domains, switches, generations |
| `SAFETY_EXPANSION` | 8 | Pharmacovigilance -> FAERS, REMS, signals |
| `BIOMARKER_EXPANSION` | 18 | Markers -> assays, cutoffs, outcomes |
| `REGULATORY_EXPANSION` | 8 | Regulatory -> designations, pathways |
| `SEQUENCE_EXPANSION` | 8 | Molecular -> scFv, CDR, affinity, nanobody |
| `REALWORLD_EXPANSION` | 10 | RWE -> registries, populations, disparities |
| `IMMUNOGENICITY_EXPANSION` | 12 | HLA/ADA -> epitopes, deimmunization, testing |

#### Expansion Algorithm

```python
def expand_query(query: str) -> List[str]:
    query_lower = query.lower()
    matched_terms: Set[str] = set()

    for category, mapping in ALL_EXPANSION_MAPS:
        for keyword, terms in mapping.items():
            if keyword in query_lower:
                matched_terms.update(terms)

    return sorted(matched_terms)
```

The algorithm is intentionally simple: substring matching against lowercased keywords. This trades precision for speed -- expansion runs in microseconds and the downstream vector search handles semantic relevance.

A category-grouped variant `expand_query_by_category()` is also available for cases where different collections should weight different expansion categories.

---

### 3e. Agent

**File:** `src/agent.py` (309 lines)
**Class:** `CARTIntelligenceAgent`

The agent wraps the RAG engine with planning and reasoning capabilities, implementing the **plan -> search -> synthesize -> report** pattern.

#### Agent Pipeline

```
                  CARTIntelligenceAgent.run(question)
                               |
                  Phase 1: search_plan()
                      |-- Identify target antigens
                      |-- Identify relevant CAR-T stages
                      |-- Determine strategy: broad/targeted/comparative
                      |-- Decompose into sub-questions
                               |
                  Phase 2: rag.retrieve(query)
                               |
                  Phase 3: evaluate_evidence()
                      |-- "sufficient"   (>=3 collections, >=10 hits)
                      |-- "partial"      (>=2 collections, >=5 hits)
                      |-- "insufficient" (anything less)
                               |
                  Phase 4: Sub-question expansion (if insufficient)
                      |-- Execute up to 2 sub-queries
                      |-- Merge additional hits into evidence
                               |
                  Phase 5: rag.query() -> LLM synthesis
                               |
                  Phase 6: Build AgentResponse
```

#### SearchPlan Dataclass

```python
@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str]
    target_antigens: List[str]         # e.g., ["CD19", "BCMA"]
    relevant_stages: List[CARTStage]   # e.g., [CLINICAL, TESTING]
    search_strategy: str               # "broad", "targeted", "comparative"
    sub_questions: List[str]           # Decomposed sub-queries
```

#### Strategy Selection Logic

- **Comparative:** Question contains "compare", "vs", or "versus"
- **Targeted:** Specific antigens mentioned AND at most one development stage
- **Broad:** Everything else (multi-stage, exploratory questions)

#### Sub-Question Decomposition

For "why did X fail?" questions, the agent generates:
1. "What are the resistance mechanisms for [antigen] therapy?"
2. "What manufacturing issues lead to CAR-T therapy failure?"
3. "What patient factors predict poor CAR-T response?"

This ensures that even when the primary search returns thin evidence, the agent can find relevant results by reframing the question.

---

### 3f. Export System

**File:** `src/export.py` (1,487 lines)
**Version:** 1.2.0

Three export formats, each with full support for both standard and comparative query results.

#### Public Functions

```python
export_markdown(query, response_text, evidence, comp_result, filters_applied) -> str
export_json(query, response_text, evidence, comp_result, filters_applied) -> str
export_pdf(query, response_text, evidence, comp_result, filters_applied) -> bytes
```

#### Markdown Export

Generates a structured report with:
- Header: query, timestamp, active filters
- Response section: the LLM-generated answer
- Evidence sources: collection-specific tables (see below)
- Knowledge graph context: if available
- Search metrics: total results, collections searched, latency
- Footer: version stamp

#### Collection-Specific Evidence Tables

Each collection renders a different set of columns optimized for its data type:

| Collection | Table Columns |
|-----------|---------------|
| Literature | #, ID, Score, Source, Title, Year, Target, Journal |
| Trial | #, NCT ID, Score, Source, Phase, Status, Sponsor, Enrollment |
| Construct | #, ID, Score, Name, Generation, Costimulatory, FDA Status |
| Assay | #, ID, Score, Type, Cell Line, Metric, Value, Outcome |
| Manufacturing | #, ID, Score, Process Step, Parameter, Batch |
| Safety | #, ID, Score, Product, Event, Severity, Onset, Source |
| Biomarker | #, ID, Score, Biomarker, Type, Method, Cutoff, Outcome |
| Regulatory | #, ID, Score, Product, Event, Date, Agency, Decision |
| Sequence | #, ID, Score, Construct, Target, Clone, Kd, Origin |
| RealWorld | #, ID, Score, Study Type, Product, Population, Endpoint, Outcome |
| Genomic | #, ID, Score, Gene, Consequence, Impact, Clinical Significance, AlphaMissense |

#### PDF Export

Uses reportlab Platypus with NVIDIA-themed styling:

- **Color palette:** NVIDIA Green (`#76B900`), Dark BG (`#1a1a1a`), Light Gray (`#666666`)
- **Layout:** Letter size, 0.6-inch margins
- **Features:**
  - Green header row on all tables with white text
  - Alternating row colors (`#f0f0f0`) for readability
  - Clickable PubMed and ClinicalTrials.gov links in evidence tables
  - Markdown-to-flowable conversion: handles headings, bold, bullet lists, block quotes, and inline tables from LLM responses
  - NVIDIA green horizontal rules for section separation

#### JSON Export

Uses Pydantic `model_dump()` for proper serialization:

```json
{
  "report_type": "cart_intelligence_query",
  "version": "1.2.0",
  "generated_at": "2026-02-20T14:30:25+00:00",
  "query": "...",
  "response": "...",
  "is_comparative": false,
  "filters_applied": {},
  "evidence": { ... },
  "search_metrics": {
    "total_results": 23,
    "collections_searched": 11,
    "search_time_ms": 142.3
  }
}
```

---

### 3g. Models

**File:** `src/models.py`
**Contents:** 10 collection models, 13 enums, 4 search/agent models

#### Enums (14)

| Enum | Values | Used By |
|------|--------|---------|
| `CARTStage` | target_id, car_design, vector_eng, testing, clinical | Agent planning, literature classification |
| `SourceType` | pubmed, pmc, patent, preprint, manual | Literature records |
| `TrialPhase` | Early Phase 1 through Phase 4, N/A | Clinical trial records |
| `TrialStatus` | Recruiting, Active, Completed, Terminated, Withdrawn, Suspended, Not yet, Unknown | Clinical trial records |
| `CARGeneration` | 1st, 2nd, 3rd, 4th, armored, universal | Construct and trial records |
| `AssayType` | cytotoxicity, cytokine, flow, proliferation, in_vivo, persistence, exhaustion, migration, trafficking, serial_killing | Assay records |
| `ProcessStep` | transduction, expansion, harvest, formulation, release, cryo, non_viral, mrna_electroporation, crispr_knock_in, ipsc_derived, automated | Manufacturing records |
| `FDAStatus` | approved, bla_filed, phase3, phase2, phase1, preclinical, discontinued | Construct records |
| `SafetyEventType` | CRS, ICANS, cytopenia, infection, secondary_malignancy, organ_toxicity, neurologic, cardiac, coagulopathy, renal | Safety records |
| `BiomarkerType` | predictive, prognostic, pharmacodynamic, monitoring, resistance | Biomarker records |
| `EvidenceLevel` | validated, emerging, exploratory | Biomarker records |
| `RegulatoryEvent` | BLA, breakthrough_therapy, RMAT, accelerated_approval, full_approval, label_update, REMS, post_marketing_requirement, complete_response | Regulatory records |
| `RWEStudyType` | retrospective, registry, claims, ehr_analysis, meta_analysis | Real-world evidence records |

#### Collection Models (10)

Each model maps to a Milvus collection schema and provides a `to_embedding_text()` method that generates the string input for the BGE-small-en-v1.5 embedding model:

```python
class CARTLiterature(BaseModel):
    id: str                     # PMID or patent number
    title: str                  # max 500 chars
    text_chunk: str             # max 3000 chars -- embedded text
    source_type: SourceType
    year: int                   # 1990-2030
    cart_stage: CARTStage
    target_antigen: str
    disease: str
    keywords: str
    journal: str

    def to_embedding_text(self) -> str:
        parts = [self.title]
        if self.text_chunk: parts.append(self.text_chunk)
        if self.target_antigen: parts.append(f"Target: {self.target_antigen}")
        if self.disease: parts.append(f"Disease: {self.disease}")
        return " ".join(parts)
```

All 10 models follow this pattern: Pydantic fields with validation, plus `to_embedding_text()` for ingest.

#### Search/Agent Models (4)

| Model | Purpose |
|-------|---------|
| `SearchHit` | Single search result: collection, id, score, text, metadata dict |
| `CrossCollectionResult` | Merged multi-collection results with `hits_by_collection()` grouping and `hit_count` property |
| `ComparativeResult` | Two-entity comparison: entity_a/b names, evidence_a/b, comparison_context |
| `AgentQuery` | Input: question, optional target_antigen, optional cart_stage, include_genomic flag |
| `AgentResponse` | Output: question, answer, evidence, knowledge_used list, timestamp |

---

### 3h. Ingest Pipeline

**Directory:** `src/ingest/`
**Files:** 15 parser modules + `__init__.py`

#### Base Class

```python
class BaseIngestPipeline(ABC):
    def __init__(self, collection_manager, embedder): ...

    @abstractmethod
    def fetch(self, **kwargs) -> Any: ...

    @abstractmethod
    def parse(self, raw_data) -> List[BaseModel]: ...

    def embed_and_store(self, records, collection_name, batch_size=32) -> int: ...
    def run(self, collection_name=None, batch_size=32, **fetch_kwargs) -> int: ...
```

The `embed_and_store()` method handles:
1. Calling `record.to_embedding_text()` on each Pydantic model
2. Batch encoding via the embedder (BGE-small-en-v1.5)
3. Enum-to-string conversion for Milvus compatibility
4. UTF-8 byte-length truncation to prevent VARCHAR overflow
5. Batch insertion via `collection_manager.insert_batch()`
6. Error handling with continue-on-failure per batch

#### Ingest Pipeline Flow

```
    fetch(**kwargs)                 parse(raw_data)              embed_and_store()
         |                              |                              |
         v                              v                              v
  +------+--------+            +--------+---------+           +--------+--------+
  | API call /    |            | XML/JSON -> list |           | to_embedding_   |
  | file read /   |  raw data  | of validated     |  records  | text() on each  |
  | web scrape    |----------->| Pydantic models  |---------->| Embed batch     |
  |               |            |                  |           | Insert to Milvus|
  +--------------+             +------------------+           +-----------------+
```

#### Parser Inventory

| Parser | Source | Collection | Notes |
|--------|--------|-----------|-------|
| `literature_parser.py` | PubMed E-utilities API | cart_literature | CAR-T-specific MeSH query, XML parsing via Biopython |
| `clinical_trials_parser.py` | ClinicalTrials.gov v2 API | cart_trials | JSON API with CAR-T condition filter |
| `construct_parser.py` | Local seed JSON | cart_constructs | FDA-approved products + investigational constructs |
| `assay_parser.py` | Local seed JSON | cart_assays | Cytotoxicity, cytokine, in vivo assay results |
| `manufacturing_parser.py` | Local seed JSON | cart_manufacturing | Process parameters and release criteria |
| `safety_parser.py` | Local seed JSON | cart_safety | Pharmacovigilance and toxicity data |
| `biomarker_parser.py` | Local seed JSON | cart_biomarkers | Predictive and PD biomarkers |
| `regulatory_parser.py` | Local seed JSON | cart_regulatory | FDA approval milestones |
| `sequence_parser.py` | Local seed JSON | cart_sequences | scFv sequences and binding affinity data |
| `realworld_parser.py` | Local seed JSON | cart_realworld | Registry and observational outcomes |
| `faers_parser.py` | FDA FAERS API | cart_safety | Live adverse event reports for CAR-T products |
| `dailymed_parser.py` | DailyMed SPL API | cart_safety, cart_regulatory | FDA drug label sections (boxed warnings, dosing) |
| `uniprot_parser.py` | UniProt REST API | cart_sequences | Protein sequences and annotations for target antigens |
| `cibmtr_parser.py` | CIBMTR public data | cart_realworld | Transplant/cell therapy registry outcomes |

The first 10 parsers work with static seed data (JSON files in `data/reference/`). The last 4 are live data fetchers that pull from external APIs at runtime or on a scheduled basis.

---

### 3i. Metrics

**File:** `src/metrics.py` (404 lines)

Prometheus metrics with graceful degradation when `prometheus_client` is not installed.

#### Metric Inventory

| Type | Name | Labels | Purpose |
|------|------|--------|---------|
| Histogram | `cart_query_latency_seconds` | query_type | Processing time distribution (buckets: 0.1s to 30s) |
| Histogram | `cart_evidence_count` | -- | Evidence items per query (buckets: 0 to 30) |
| Counter | `cart_queries_total` | query_type, status | Total queries (rag/agent/comparative/entity_link x success/error) |
| Counter | `cart_collection_hits_total` | collection | Hits per collection across all queries |
| Counter | `cart_llm_tokens_total` | direction | LLM token usage (prompt/completion) |
| Gauge | `cart_active_connections` | -- | Current active connections |
| Gauge | `cart_collection_size` | collection | Records per collection |
| Gauge | `cart_last_ingest_timestamp` | source | Unix timestamp of last ingest (pubmed/clinical_trials) |

#### Graceful Degradation

```python
try:
    from prometheus_client import Counter, Gauge, Histogram, generate_latest
    _PROMETHEUS_AVAILABLE = True
except ImportError:
    _PROMETHEUS_AVAILABLE = False

    class _NoOpLabeled:
        def labels(self, *args, **kwargs): return self
        def observe(self, *args, **kwargs): pass
        def inc(self, *args, **kwargs): pass
        def set(self, *args, **kwargs): pass

    QUERY_LATENCY = _NoOpLabeled()
    # ... all metrics become no-ops
```

The rest of the application imports metrics helpers unconditionally. If prometheus_client is missing, all calls silently do nothing. This avoids a hard dependency for development and testing environments.

#### Helper Functions

```python
record_query(query_type, latency, hit_count, status)  # One-call update for all query metrics
record_collection_hits(hits_by_collection)              # Per-collection hit counters
update_collection_sizes(stats)                          # Set current record counts
get_metrics_text() -> str                               # Prometheus exposition format
```

---

### 3j. Scheduler

**File:** `src/scheduler.py` (226 lines)
**Class:** `IngestScheduler`

Automated refresh of PubMed and ClinicalTrials.gov collections using APScheduler's `BackgroundScheduler`.

#### Architecture

```
  IngestScheduler
       |
       +-- BackgroundScheduler (daemon=True)
            |
            +-- Job: refresh_pubmed
            |     |-- PubMedIngestPipeline.run()
            |     |-- Update LAST_INGEST gauge
            |
            +-- Job: refresh_clinical_trials
                  |-- ClinicalTrialsIngestPipeline.run()
                  |-- Update LAST_INGEST gauge
```

#### Configuration

| Setting | Default | Description |
|---------|---------|-------------|
| `INGEST_SCHEDULE_HOURS` | 168 | Refresh interval (168h = weekly) |
| `INGEST_ENABLED` | False | Master switch for scheduler |

#### Graceful Degradation

Like metrics, the scheduler exports a no-op stub class when `apscheduler` is not installed:

```python
try:
    from apscheduler.schedulers.background import BackgroundScheduler
    _APSCHEDULER_AVAILABLE = True
except ImportError:
    _APSCHEDULER_AVAILABLE = False
```

#### Public API

```python
scheduler = IngestScheduler(collection_manager, embedder, interval_hours=168)
scheduler.start()        # Non-blocking, runs in daemon thread
scheduler.get_status()   # Returns next_run_time, last_run_time, job_count
scheduler.stop()         # Graceful shutdown
```

---

## 4. Data Flow Diagrams

### Query Flow

```
User Question: "What are the resistance mechanisms to BCMA-targeted CAR-T?"
     |
     v
+----+---------------------------+
| AgentQuery                     |
|   question: "What are the..."  |
|   target_antigen: None         |
|   include_genomic: True        |
+----+---------------------------+
     |
     v
+----+---------------------------+
| Query Expansion                |
|   "bcma" matched:              |
|     -> BCMA, B-cell maturation |
|        antigen, TNFRSF17,     |
|        multiple myeloma,       |
|        Abecma, Carvykti, ...  |
|   "resistance" matched:        |
|     -> antigen loss, escape,   |
|        lineage switch, ...    |
+----+---------------------------+
     |
     v
+----+---------------------------+
| Embed Query                    |
|   "Represent this sentence     |
|    for searching relevant      |
|    passages: What are the      |
|    resistance mechanisms..."   |
|   -> [384-dim float vector]    |
+----+---------------------------+
     |
     +----> Parallel Search (11 threads)
     |        |
     |        +-> cart_literature:    5 hits (weight 0.20)
     |        +-> cart_trials:        3 hits (weight 0.16)
     |        +-> cart_constructs:    4 hits (weight 0.10)
     |        +-> cart_assays:        2 hits (weight 0.09)
     |        +-> cart_safety:        5 hits (weight 0.08)
     |        +-> cart_biomarkers:    3 hits (weight 0.08)
     |        +-> cart_manufacturing: 1 hit  (weight 0.07)
     |        +-> cart_realworld:     2 hits (weight 0.07)
     |        +-> cart_regulatory:    2 hits (weight 0.06)
     |        +-> cart_sequences:     3 hits (weight 0.06)
     |        +-> genomic_evidence:   0 hits (weight 0.04)
     |
     +----> Expanded Search (antigen: field filter, non-antigen: re-embed)
     |        |
     |        +-> "BCMA" -> field filter on target_antigen-capable collections
     |        +-> "antigen loss" -> re-embed and search all collections
     |
     v
+----+---------------------------+
| Merge & Rank                   |
|   Deduplicate by ID            |
|   Sort by weighted score desc  |
|   Cap at 30 results            |
+----+---------------------------+
     |
     v
+----+---------------------------+
| Knowledge Augmentation         |
|   "BCMA" -> CART_TARGETS[BCMA] |
|     protein, expression,       |
|     approved products,         |
|     resistance mechanisms,     |
|     toxicity profile           |
+----+---------------------------+
     |
     v
+----+---------------------------+
| Build LLM Prompt               |
|   ## Retrieved Evidence         |
|   ### Evidence from Literature  |
|   1. [Literature:PMID ...]     |
|   ...                          |
|   ### Knowledge Graph Context   |
|   ## Target Antigen: BCMA      |
|   ---                          |
|   ## Question                   |
|   What are the resistance...   |
+----+---------------------------+
     |
     v
+----+---------------------------+
| Claude Sonnet 4.6              |
|   System: CART_SYSTEM_PROMPT   |
|   max_tokens: 2048             |
|   temperature: 0.7             |
+----+---------------------------+
     |
     v
+----+---------------------------+
| Response + Evidence Panel      |
|   Answer with citations        |
|   Evidence grouped by          |
|   collection with scores       |
+--------------------------------+
```

### Ingest Flow

```
  External Data Source                    Internal Pipeline
  ====================                    =================

  PubMed E-utilities          fetch()     Parse XML records
  ClinicalTrials.gov API  ------------>   Validate with Pydantic models
  FDA FAERS API                           Call to_embedding_text()
  DailyMed SPL API                        |
  UniProt REST API                        v
  CIBMTR Public Data           embed_and_store()
  Local JSON seed files    <----------    Batch encode (BGE-small, 32/batch)
                                          Insert into Milvus collection
                                          Log progress
                                          |
                                          v
                                     Milvus Collection
                                     (IVF_FLAT, COSINE)
```

### Export Flow

```
  CrossCollectionResult
  (or ComparativeResult)
         |
         +---> export_markdown()
         |        |
         |        +-> Header (query, timestamp, filters)
         |        +-> Response section
         |        +-> Collection-specific evidence tables
         |        +-> Knowledge graph context
         |        +-> Search metrics table
         |        +-> Footer with version
         |        |
         |        v
         |     Markdown string (.md)
         |
         +---> export_json()
         |        |
         |        +-> Pydantic model_dump()
         |        +-> json.dumps(indent=2)
         |        |
         |        v
         |     JSON string (.json)
         |
         +---> export_pdf()
                  |
                  +-> SimpleDocTemplate (letter, 0.6" margins)
                  +-> NVIDIA-themed styles (green headers, alt rows)
                  +-> _md_to_flowables() for response text
                  +-> _build_pdf_evidence_table() per collection
                  +-> Clickable PubMed/CT.gov links
                  |
                  v
               PDF bytes (.pdf)
```

---

## 5. API Layer

**Files:** `api/main.py` (588 lines), `api/routes/events.py` (122 lines), `api/routes/meta_agent.py` (142 lines), `api/routes/reports.py` (180 lines)
**Framework:** FastAPI with Pydantic request/response schemas

### Lifespan Management

The API uses FastAPI's `@asynccontextmanager` lifespan pattern for startup/shutdown:

```python
@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    _manager = CARTCollectionManager(host, port)
    _manager.connect()
    embedder = _Embedder()              # SentenceTransformer wrapper
    llm_client = _LLMClient()           # Anthropic client wrapper
    _engine = CARTRAGEngine(manager, embedder, llm_client, knowledge, query_expander)
    yield
    # Shutdown
    _manager.disconnect()
```

### Endpoints

| Method | Path | Tag | Description |
|--------|------|-----|-------------|
| GET | `/health` | status | Service health with collection count and total vector count. Returns 503 if engine or Milvus unavailable |
| GET | `/collections` | status | All collection names and record counts |
| POST | `/query` | rag | Full RAG: retrieve evidence + LLM synthesis. Requires both embedder and LLM |
| POST | `/search` | rag | Evidence-only retrieval (no LLM). Fast mode for when only evidence snippets are needed |
| POST | `/find-related` | rag | Cross-collection entity linking. "Show me everything about Yescarta" |
| GET | `/knowledge/stats` | knowledge | Knowledge graph statistics (target counts, toxicity profiles, etc.) |
| GET | `/metrics` | monitoring | Prometheus-compatible metrics exposition |
| GET | `/events` | events | List pipeline events with optional filters |
| GET | `/events/{event_id}` | events | Get a specific pipeline event by ID |
| POST | `/ask` | meta_agent | Meta-agent question answering endpoint |
| GET | `/reports/{patient_id}` | reports | Get report for a specific patient |
| GET | `/reports/{patient_id}/{fmt}` | reports | Get report in a specific format |

### Request/Response Schemas

**QueryRequest:**
```python
class QueryRequest(BaseModel):
    question: str                          # Required, min_length=1
    target_antigen: Optional[str] = None   # Filter by target (e.g., "CD19")
    collections: Optional[List[str]] = None # Restrict to specific collections
    year_min: Optional[int] = None         # Minimum publication year (1990-2030)
    year_max: Optional[int] = None         # Maximum publication year (1990-2030)
```

**QueryResponse:**
```python
class QueryResponse(BaseModel):
    question: str
    answer: str
    evidence: List[EvidenceItem]
    knowledge_context: str = ""
    collections_searched: int = 0
    search_time_ms: float = 0.0
```

### CORS Configuration

```python
_cors_origins = [o.strip() for o in settings.CORS_ORIGINS.split(",") if o.strip()]
app.add_middleware(
    CORSMiddleware,
    allow_origins=_cors_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

CORS is restricted to 3 origins configured in `config/settings.py`: `http://localhost:8080` (landing page), `http://localhost:8521` (Streamlit UI), and `http://localhost:8522` (FastAPI self).

### API Key Loading

The API loads the Anthropic API key from the rag-chat-pipeline `.env` file if `ANTHROPIC_API_KEY` is not already set in the environment. This allows the API to share credentials with the parent platform without duplication.

---

## 6. UI Architecture

**File:** `app/cart_ui.py` (1,162 lines)
**Framework:** Streamlit v1.30+
**Port:** 8521

### Engine Initialization

Uses `@st.cache_resource` to initialize the RAG engine once and share it across all Streamlit reruns:

```python
@st.cache_resource
def init_engine():
    manager = CARTCollectionManager()
    manager.connect()
    embedder = SimpleEmbedder()       # BGE-small-en-v1.5
    llm_client = SimpleLLMClient()    # Anthropic wrapper with streaming
    engine = CARTRAGEngine(manager, embedder, llm_client, knowledge, query_expander)
    return engine, manager
```

### Feature Set

| Feature | Description |
|---------|-------------|
| **Collection Stats** | All 11 collections displayed in sidebar with record counts |
| **Deep Research Mode** | Autonomous agent pipeline (plan-search-synthesize-report) |
| **Conversation Memory** | Maintains history for follow-up queries (configurable depth via `MAX_CONVERSATION_CONTEXT`) |
| **Collection Filtering** | Sidebar checkboxes to restrict search to specific collections |
| **Temporal Filtering** | Year range sliders for date-bounded queries |
| **Target Antigen Filter** | Dropdown for specific antigen focus |
| **CAR-T Stage Filter** | Filter by development stage (target_id, car_design, vector_eng, testing, clinical) |
| **Citation Relevance** | Evidence items tagged as high/medium/low relevance |
| **Image Upload** | Upload slides/images for claim verification against the knowledge base |
| **Knowledge Graph Viz** | Interactive PyVis network visualization of knowledge graph entities |
| **Streaming Response** | Token-by-token streaming via `query_stream()` |
| **Export Controls** | Download buttons for Markdown, JSON, and PDF reports |
| **Comparative Analysis** | Automatic detection and structured comparison rendering |

### Session State Management

Key session state variables:

```python
st.session_state.messages           # Conversation history
st.session_state.current_evidence   # Latest CrossCollectionResult
st.session_state.current_comp       # Latest ComparativeResult (if any)
st.session_state.engine             # Cached CARTRAGEngine reference
st.session_state.collection_filter  # Selected collections list
```

### CSS Theming

Configured via `.streamlit/config.toml` for dark theme with custom styling. NVIDIA green (`#76B900`) is used for accents and interactive elements.

---

## 7. Docker and Deployment

### Multi-Stage Dockerfile

```
Stage 1: builder (python:3.10-slim)
    |-- apt-get: build-essential, gcc, g++, libxml2-dev, libxslt1-dev
    |-- Create venv at /opt/venv
    |-- pip install requirements.txt
    |
Stage 2: runtime (python:3.10-slim)
    |-- apt-get: curl, libgomp1, libxml2, libxslt1.1  (minimal runtime libs)
    |-- COPY --from=builder /opt/venv /opt/venv
    |-- COPY application source
    |-- Create non-root user: cartuser
    |-- EXPOSE 8521, 8522
    |-- HEALTHCHECK against Streamlit /_stcore/health
    |-- Default CMD: streamlit run app/cart_ui.py
```

The builder stage installs compilation dependencies that are not needed at runtime, reducing the final image size.

### Docker Compose Topology

```
                    cart-network (bridge)
                          |
      +-------------------+-------------------+
      |                   |                   |
+-----+------+    +------+-------+    +------+------+
| milvus-etcd|    | milvus-minio |    | milvus-     |
| (etcd v3.5)|    | (MinIO)      |    | standalone  |
| metadata   |    | object store |    | (v2.4)      |
| store      |    |              |    | 19530, 9091 |
+------------+    +--------------+    +------+------+
                                             |
                          +------------------+------------------+
                          |                  |                  |
                   +------+------+    +------+------+    +-----+------+
                   | cart-       |    | cart-api    |    | cart-setup |
                   | streamlit   |    | (FastAPI)   |    | (one-shot) |
                   | (port 8521) |    | (port 8522) |    | collections|
                   |             |    | 2 workers   |    | + seed     |
                   +-------------+    +-------------+    +------------+
```

### Service Details

| Service | Image | Ports | Healthcheck | Restart |
|---------|-------|-------|-------------|---------|
| milvus-etcd | quay.io/coreos/etcd:v3.5.5 | -- | `etcdctl endpoint health` | unless-stopped |
| milvus-minio | minio/minio:RELEASE.2023-03-20 | -- | `curl localhost:9000/minio/health/live` | unless-stopped |
| milvus-standalone | milvusdb/milvus:v2.4-latest | 19530, 9091 | `curl localhost:9091/healthz` | unless-stopped |
| cart-streamlit | (built from Dockerfile) | 8521 | `curl localhost:8521/_stcore/health` | unless-stopped |
| cart-api | (built from Dockerfile) | 8522 | `curl localhost:8522/health` | unless-stopped |
| cart-setup | (built from Dockerfile) | -- | -- | no (one-shot) |

### Environment Variables

All settings use the `CART_` prefix for environment variable override (via Pydantic BaseSettings):

| Variable | Default | Description |
|----------|---------|-------------|
| `CART_MILVUS_HOST` | localhost | Milvus server hostname |
| `CART_MILVUS_PORT` | 19530 | Milvus gRPC port |
| `CART_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model identifier |
| `CART_EMBEDDING_DIMENSION` | 384 | Vector dimension |
| `CART_LLM_MODEL` | claude-sonnet-4-6 | Claude model identifier |
| `CART_TOP_K_PER_COLLECTION` | 5 | Max results per collection per query |
| `CART_SCORE_THRESHOLD` | 0.4 | Minimum cosine similarity threshold |
| `CART_WEIGHT_LITERATURE` | 0.20 | Literature collection weight |
| `CART_WEIGHT_TRIALS` | 0.16 | Trials collection weight |
| `CART_WEIGHT_CONSTRUCTS` | 0.10 | Constructs collection weight |
| `CART_WEIGHT_ASSAYS` | 0.09 | Assays collection weight |
| `CART_WEIGHT_SAFETY` | 0.08 | Safety collection weight |
| `CART_WEIGHT_BIOMARKERS` | 0.08 | Biomarkers collection weight |
| `CART_WEIGHT_MANUFACTURING` | 0.07 | Manufacturing collection weight |
| `CART_WEIGHT_REALWORLD` | 0.07 | Real-world evidence collection weight |
| `CART_WEIGHT_REGULATORY` | 0.06 | Regulatory collection weight |
| `CART_WEIGHT_SEQUENCES` | 0.06 | Sequences collection weight |
| `CART_WEIGHT_GENOMIC` | 0.04 | Genomic evidence collection weight |
| `CART_API_PORT` | 8522 | FastAPI server port |
| `CART_STREAMLIT_PORT` | 8521 | Streamlit server port |
| `CART_INGEST_SCHEDULE_HOURS` | 168 | Scheduler interval (hours) |
| `CART_INGEST_ENABLED` | False | Enable background ingest scheduler |
| `CART_CITATION_HIGH_THRESHOLD` | 0.75 | High relevance citation cutoff |
| `CART_CITATION_MEDIUM_THRESHOLD` | 0.60 | Medium relevance citation cutoff |
| `CART_MAX_CONVERSATION_CONTEXT` | 3 | Prior exchanges to inject for follow-ups |
| `ANTHROPIC_API_KEY` | -- | Claude API key (no CART_ prefix) |

### Volumes

| Volume | Mounted To | Purpose |
|--------|-----------|---------|
| etcd_data | /etcd | Milvus metadata persistence |
| minio_data | /minio_data | Milvus log and index object storage |
| milvus_data | /var/lib/milvus | Milvus vector data and WAL |

### Setup Sequence

The `cart-setup` service runs once after Milvus is healthy and executes 9 seed scripts in sequence:

1. `setup_collections.py --drop-existing --seed-constructs` -- Create all 10 collection schemas and seed construct data
2. `seed_knowledge.py` -- Load knowledge graph data (triggers any one-time knowledge base operations)
3. `seed_assays.py` -- Seed in-vitro/in-vivo assay records from `assay_seed_data.json`
4. `seed_manufacturing.py` -- Seed manufacturing/CMC records from `manufacturing_seed_data.json`
5. `seed_safety.py` -- Seed pharmacovigilance data from `safety_seed_data.json`
6. `seed_biomarkers.py` -- Seed biomarker records from `biomarker_seed_data.json`
7. `seed_regulatory.py` -- Seed FDA regulatory milestones from `regulatory_seed_data.json`
8. `seed_sequences.py` -- Seed molecular/structural data from `sequence_seed_data.json`
9. `seed_realworld.py` -- Seed real-world evidence from `realworld_seed_data.json`

---

## 8. Testing Architecture

**Directory:** `tests/`
**Files:** 7 test modules + conftest.py
**Total:** 415 tests, 4,321 lines

### Test Files

| File | Focus | Test Count |
|------|-------|-----------|
| `test_models.py` | All 10 collection models, 13 enums, 4 search/agent models, `to_embedding_text()` methods, validation constraints | ~80 |
| `test_rag_engine.py` | CARTRAGEngine: retrieve, query, query_stream, find_related, comparative, prompt building, citation formatting, score weighting | ~40 |
| `test_knowledge.py` | Knowledge graph: all 3 dictionaries, context retrieval, entity resolution, comparison context, `get_all_context_for_query()` | ~50 |
| `test_query_expansion.py` | All 12 expansion maps, `expand_query()`, `expand_query_by_category()`, stats, edge cases | ~30 |
| `test_agent.py` | CARTIntelligenceAgent: run, search_plan, evaluate_evidence, generate_report, sub-question decomposition | ~20 |
| `test_export.py` | Markdown, JSON, PDF export for standard and comparative results, evidence tables, citation links | ~21 |
| `test_integration.py` | Integration tests: end-to-end pipeline validation, cross-collection queries | ~37 |

### conftest.py Fixtures

The shared fixtures enable all tests to run without Milvus, embeddings, or LLM access:

```python
@pytest.fixture
def mock_embedder():
    """384-dim zero vectors for any embed_text() call."""
    embedder = MagicMock()
    embedder.embed_text.return_value = [0.0] * 384
    return embedder

@pytest.fixture
def mock_llm_client():
    """Always returns 'Mock response' for generate() and generate_stream()."""
    client = MagicMock()
    client.generate.return_value = "Mock response"
    client.generate_stream.return_value = iter(["Mock ", "response"])
    return client

@pytest.fixture
def mock_collection_manager():
    """search() -> [], search_all() -> {name: [] for all 10}, stats -> 42 each."""
    manager = MagicMock()
    manager.search.return_value = []
    manager.search_all.return_value = {name: [] for name in collection_names}
    manager.get_collection_stats.return_value = {name: 42 for name in collection_names}
    return manager

@pytest.fixture
def sample_search_hits():
    """5 SearchHit objects spanning Literature, Trial, Construct, Safety, Manufacturing."""
    return [
        SearchHit(collection="Literature", id="12345678", score=0.92, ...),
        SearchHit(collection="Trial", id="NCT03958656", score=0.87, ...),
        SearchHit(collection="Construct", id="construct-kymriah", score=0.83, ...),
        SearchHit(collection="Safety", id="safety-crs-001", score=0.78, ...),
        SearchHit(collection="Manufacturing", id="mfg-lenti-001", score=0.71, ...),
    ]

@pytest.fixture
def sample_evidence(sample_search_hits):
    """CrossCollectionResult with 5 hits, knowledge context, and search metrics."""
    return CrossCollectionResult(
        query="What is the efficacy of CD19 CAR-T therapy?",
        hits=sample_search_hits,
        knowledge_context="## Target Antigen: CD19\n...",
        total_collections_searched=10,
        search_time_ms=42.5,
    )
```

### Running Tests

```bash
cd cart_intelligence_agent
pytest tests/ -v                     # All 415 tests
pytest tests/test_rag_engine.py -v   # RAG engine tests only
pytest tests/ -k "comparative"       # Only comparison-related tests
```

---

## 9. Security Considerations

### Non-Root Docker User

The Dockerfile creates a dedicated `cartuser` and switches to it before exposing ports:

```dockerfile
RUN useradd -r -s /bin/false cartuser \
    && mkdir -p /app/data/cache /app/data/reference \
    && chown -R cartuser:cartuser /app
USER cartuser
```

### API Key Management

- The `ANTHROPIC_API_KEY` is never stored in code or configuration files
- It is loaded from environment variables or the parent pipeline's `.env` file at runtime
- The Pydantic `CARTSettings` class marks it as `Optional[str]` so the application can start without it (in degraded mode, with LLM features disabled)
- Docker Compose passes the key via `${ANTHROPIC_API_KEY}` environment variable interpolation

### CORS Configuration

CORS is configured via `config/settings.py` with `CORS_ORIGINS` defaulting to 3 origins: `http://localhost:8080,http://localhost:8521,http://localhost:8522`. For production, update the `CORS_ORIGINS` setting to include only the required origins.

### No Secrets in Code

All sensitive values are externalized:
- API keys via environment variables
- Milvus credentials via Docker Compose environment
- MinIO credentials via Docker Compose environment (default: minioadmin/minioadmin)

### Input Validation

- All API endpoints use Pydantic request schemas with `min_length`, `ge`, `le` constraints
- Collection names are validated against `COLLECTION_SCHEMAS.keys()`
- Year ranges are bounded to 1990-2030
- The `score_threshold` parameter prevents low-quality results from reaching the user

---

## 10. Performance Characteristics

### Parallel Search

The most critical performance optimization is parallel collection search. With 11 collections, sequential search would take 11x longer. The `ThreadPoolExecutor` with `max_workers=len(collections)` (11) runs all searches concurrently:

```
Sequential:  11 x ~50ms = ~550ms
Parallel:    max(~50ms each) = ~50-80ms  (limited by slowest collection)
```

### Score-Based Filtering

The `SCORE_THRESHOLD` (default 0.4) filters out low-relevance results at the Milvus level, reducing the amount of data transferred and processed:

```python
if score < score_threshold:
    continue
```

### Vector Index Optimization

**IVF_FLAT** with `nlist=1024` and `nprobe=16` provides:
- Good recall (searching 16 of 1024 partitions)
- Fast search time (~milliseconds per collection for <100K records)
- Exact distance computation within probed partitions (no quantization loss)

For collections growing beyond 1M records, consider switching to **IVF_PQ** or **HNSW** for better scalability.

### Embedding Batch Processing

Ingest pipelines batch embed 32 records at a time, leveraging the sentence-transformers library's batched encoding for GPU/CPU efficiency:

```python
embeddings = self.embedder.encode(texts)  # Batch of 32 texts at once
```

### Result Capping

The merge-and-rank step caps results at 30 to prevent excessive prompt lengths:

```python
unique.sort(key=lambda h: h.score, reverse=True)
return unique[:30]
```

This keeps the LLM prompt within a reasonable token budget while retaining the most relevant evidence.

### Streamlit Caching

The `@st.cache_resource` decorator ensures the embedding model, LLM client, and Milvus connection are initialized only once, even across multiple Streamlit reruns. This eliminates the ~5-10 second model loading time on each interaction.

---

## 11. Extension Points

### How to Add a New Collection

1. **Define the Milvus schema** in `src/collections.py`:

```python
NEW_DOMAIN_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, is_primary=True, max_length=100),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
    # ... domain-specific fields
]

NEW_DOMAIN_SCHEMA = CollectionSchema(
    fields=NEW_DOMAIN_FIELDS,
    description="Description of the new collection",
)
```

2. **Register in the two dictionaries**:

```python
COLLECTION_SCHEMAS["cart_new_domain"] = NEW_DOMAIN_SCHEMA
COLLECTION_MODELS["cart_new_domain"] = NewDomainRecord  # or None if read-only
```

3. **Create the Pydantic model** in `src/models.py`:

```python
class NewDomainRecord(BaseModel):
    id: str = Field(...)
    text_summary: str = Field(...)
    # ... fields matching the schema

    def to_embedding_text(self) -> str:
        return f"{self.text_summary}"
```

4. **Add to COLLECTION_CONFIG** in `src/rag_engine.py`:

```python
"cart_new_domain": {
    "weight": 0.05,
    "label": "NewDomain",
    "has_target_antigen": False,
    "year_field": None,
},
```

5. **Add a collection weight setting** in `config/settings.py`:

```python
WEIGHT_NEW_DOMAIN: float = 0.05
```

6. **Add collection-specific evidence table columns** in `src/export.py` (both `_format_evidence_table()` for Markdown and `_build_pdf_evidence_table()` for PDF).

7. **Add seed data** in `data/reference/new_domain_seed_data.json` and a corresponding `scripts/seed_new_domain.py`.

8. **Add to docker-compose.yml** `cart-setup` command sequence.

9. **Update tests** in `tests/conftest.py` to include the new collection in `mock_collection_manager`.

### How to Add a New Ingest Source

1. **Create a parser** in `src/ingest/new_source_parser.py`:

```python
from src.ingest.base import BaseIngestPipeline
from src.models import NewDomainRecord

class NewSourceIngestPipeline(BaseIngestPipeline):
    COLLECTION_NAME = "cart_new_domain"

    def fetch(self, **kwargs) -> Any:
        # API call, file read, or web scrape
        response = requests.get("https://api.example.com/data")
        return response.json()

    def parse(self, raw_data: Any) -> List[NewDomainRecord]:
        records = []
        for item in raw_data:
            records.append(NewDomainRecord(
                id=item["id"],
                text_summary=item["description"],
                # ... map fields
            ))
        return records

    def run(self, **kwargs) -> int:
        return super().run(
            collection_name=self.COLLECTION_NAME,
            **kwargs,
        )
```

2. **Optionally add to the scheduler** in `src/scheduler.py` if the source should be refreshed periodically.

### How to Add a New Knowledge Domain

1. **Add a new dictionary** in `src/knowledge.py`:

```python
CART_NEW_DOMAIN: Dict[str, Dict[str, Any]] = {
    "entry_key": {
        "field_1": "value",
        "field_2": ["list", "of", "values"],
    },
}
```

2. **Create a context retrieval function**:

```python
def get_new_domain_context(key: str) -> str:
    data = CART_NEW_DOMAIN.get(key)
    if not data:
        return ""
    lines = [f"## New Domain: {key}"]
    lines.append(f"- **Field 1:** {data['field_1']}")
    return "\n".join(lines)
```

3. **Wire into `get_all_context_for_query()`** to add keyword-based scanning.

4. **Wire into `CARTRAGEngine._get_knowledge_context()`** for query-time augmentation.

5. **Add entries to `ENTITY_ALIASES`** if the new domain should participate in comparative analysis.

6. **Update `get_knowledge_stats()`** to include the new dictionary count.

### How to Add a New Query Expansion Map

1. **Define the expansion dictionary** in `src/query_expansion.py`:

```python
NEW_DOMAIN_EXPANSION: Dict[str, List[str]] = {
    "keyword_1": ["term_1", "term_2", "term_3"],
    "keyword_2": ["term_4", "term_5"],
}
```

2. **Register in `ALL_EXPANSION_MAPS`**:

```python
ALL_EXPANSION_MAPS: List[tuple] = [
    # ... existing maps
    ("NewDomain", NEW_DOMAIN_EXPANSION),
]
```

The `expand_query()` function automatically picks up new maps from this list.

---

## Cross-Agent Integration

### The 11 HCLS AI Factory Intelligence Agents

| # | Agent | Domain |
|---|-------|--------|
| 1 | Precision Oncology Agent | Molecular tumor board decision support |
| 2 | **CAR-T Intelligence Agent** (this agent) | CAR-T cell therapy development lifecycle |
| 3 | Precision Biomarker Agent | Genotype-aware biomarker interpretation |
| 4 | Clinical Trial Intelligence Agent | Trial landscape monitoring and matching |
| 5 | Cardiology Intelligence Agent | Cardiovascular risk and treatment optimization |
| 6 | Neurology Intelligence Agent | Neurological disease characterization |
| 7 | Pharmacogenomics (PGx) Agent | Drug-gene interaction and dosing guidance |
| 8 | Medical Imaging Intelligence Agent | Imaging AI with NVIDIA NIM microservices |
| 9 | Single-Cell Intelligence Agent | Single-cell omics analysis |
| 10 | Autoimmune Intelligence Agent | Autoimmune disease monitoring and treatment |
| 11 | Rare Disease Intelligence Agent | Rare disease diagnosis and therapeutic matching |

### Cross-Agent Calls via cross_modal/cross_agent.py

The CAR-T Agent calls 5 peer agents via `cross_modal/cross_agent.py` and the `/integrated-assessment` endpoint:

| Target Agent | Integration Purpose |
|-------------|-------------------|
| Precision Biomarker Agent | CRS/ICANS biomarker correlation (ferritin, CRP, IL-6) |
| Precision Oncology Agent | Variant-level therapy matching for CAR-T-eligible histologies |
| Single-Cell Intelligence Agent | T-cell phenotype and exhaustion marker analysis |
| Cardiology Intelligence Agent | Cardiac monitoring during CRS management |
| Clinical Trial Intelligence Agent | CAR-T trial landscape and eligibility matching |

### Integrated Assessment Endpoint

The `/integrated-assessment` endpoint orchestrates multi-agent workflows for CAR-T patients by:

1. Collecting CAR-T eligibility assessment and product recommendation
2. Requesting biomarker correlation from the Biomarker Agent (CRS predictors)
3. Querying the Oncology Agent for prior therapy history and variant status
4. Checking cardiac risk via the Cardiology Agent (pre-lymphodepletion assessment)
5. Aggregating findings into a unified CAR-T treatment plan

---

## Appendix A: Directory Structure

```
cart_intelligence_agent/
├── api/
│   ├── __init__.py
│   ├── main.py                          # FastAPI REST API (588 lines)
│   └── routes/
│       ├── __init__.py
│       ├── events.py                    # Pipeline events (122 lines)
│       ├── meta_agent.py               # Meta-agent endpoint (142 lines)
│       └── reports.py                   # Report generation (180 lines)
├── app/
│   └── cart_ui.py                       # Streamlit UI (1,162 lines)
├── config/
│   └── settings.py                      # Pydantic BaseSettings (113 lines)
├── data/
│   └── reference/
│       ├── assay_seed_data.json         # 75 records
│       ├── biomarker_seed_data.json     # 60 records
│       ├── constructs_seed_data.json    # 41 records
│       ├── immunogenicity_biomarker_seed.json  # 20 records
│       ├── immunogenicity_sequence_seed.json   # 18 records
│       ├── literature_seed_data.json    # 60 records
│       ├── manufacturing_seed_data.json # 56 records
│       ├── patent_seed_data.json        # 45 records
│       ├── realworld_seed_data.json     # 54 records
│       ├── regulatory_seed_data.json    # 40 records
│       ├── safety_seed_data.json        # 71 records
│       ├── sequence_seed_data.json      # 40 records
│       └── trials_seed_data.json        # 69 records
├── docs/
│   └── ARCHITECTURE_GUIDE.md            # This document
├── scripts/                             # Setup and seed scripts (1,686 lines)
├── src/
│   ├── __init__.py
│   ├── agent.py                         # Autonomous agent (309 lines)
│   ├── collections.py                   # Milvus collection manager
│   ├── export.py                        # MD/JSON/PDF export (1,487 lines)
│   ├── knowledge.py                     # Knowledge graph (3 dicts, 71 entries)
│   ├── metrics.py                       # Prometheus metrics (404 lines)
│   ├── models.py                        # 10 models, 13 enums
│   ├── query_expansion.py              # 12 maps, 229 keywords
│   ├── rag_engine.py                    # Core RAG engine (754 lines)
│   ├── scheduler.py                     # APScheduler ingest (226 lines)
│   └── ingest/
│       ├── __init__.py
│       ├── base.py                      # BaseIngestPipeline ABC
│       ├── literature_parser.py         # PubMed
│       ├── clinical_trials_parser.py    # ClinicalTrials.gov
│       ├── construct_parser.py          # CAR constructs
│       ├── assay_parser.py              # Assay results
│       ├── manufacturing_parser.py      # Manufacturing/CMC
│       ├── safety_parser.py             # Pharmacovigilance
│       ├── biomarker_parser.py          # Biomarkers
│       ├── regulatory_parser.py         # FDA regulatory
│       ├── sequence_parser.py           # Molecular data
│       ├── realworld_parser.py          # Real-world evidence
│       ├── faers_parser.py              # FDA FAERS (live)
│       ├── dailymed_parser.py           # DailyMed labels (live)
│       ├── uniprot_parser.py            # UniProt sequences (live)
│       └── cibmtr_parser.py             # CIBMTR registry (live)
├── tests/
│   ├── __init__.py
│   ├── conftest.py                      # Shared fixtures
│   ├── test_agent.py
│   ├── test_export.py
│   ├── test_knowledge.py
│   ├── test_models.py
│   ├── test_query_expansion.py
│   └── test_rag_engine.py
├── .streamlit/
│   └── config.toml                      # Dark theme configuration
├── Dockerfile                           # Multi-stage build
├── docker-compose.yml                   # 6 services
└── requirements.txt                     # 22 packages
```

## Appendix B: Package Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| pydantic | >=2.0 | Data validation, model serialization |
| pydantic-settings | >=2.0 | Environment-driven configuration |
| loguru | >=0.7.0 | Structured logging |
| pymilvus | >=2.4.0 | Milvus vector database client |
| sentence-transformers | >=2.2.0 | BGE-small-en-v1.5 embedding model |
| anthropic | >=0.18.0 | Claude LLM client |
| streamlit | >=1.30.0 | Chat UI framework |
| fastapi | >=0.109.0 | REST API framework |
| uvicorn[standard] | >=0.27.0 | ASGI server |
| python-multipart | >=0.0.6 | File upload support |
| requests | >=2.31.0 | HTTP client for ingest pipelines |
| lxml | >=5.0.0 | PubMed XML parsing |
| biopython | >=1.83 | NCBI E-utilities interface |
| apscheduler | >=3.10.0 | Background job scheduling |
| prometheus-client | >=0.20.0 | Prometheus metrics exposition |
| reportlab | >=4.0.0 | PDF generation (Platypus) |
| pyvis | >=0.3.0 | Knowledge graph visualization |
| numpy | >=1.24.0 | Numerical operations |
| tqdm | >=4.65.0 | Progress bars for ingest |
| python-dotenv | >=1.0.0 | .env file loading |

---

*This document is the definitive architecture reference for the CAR-T Intelligence Agent. For questions or contributions, contact Adam Jones.*
