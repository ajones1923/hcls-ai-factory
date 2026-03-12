# Precision Oncology Intelligence Agent -- Advanced Learning Guide

**Audience:** Experienced developers contributing to or extending the agent.
**Prerequisite reading:** `LEARNING_GUIDE.md` (introductory guide).
**Codebase snapshot:** February 2026 -- ~14,000 lines across `src/`, `api/`, `app/`, `tests/`, `scripts/`.

---

## Prerequisites

Before diving in you should be comfortable with:

1. **Python 3.10+** -- dataclasses, async generators, `concurrent.futures`, Pydantic v2.
2. **Vector databases** -- IVF indexing, cosine similarity, embedding pipelines.
   Specifically Milvus (pymilvus SDK, collection schemas, filter expressions).
3. **Clinical genomics vocabulary** -- VCF format, SnpEff / VEP annotations,
   somatic vs. germline, AMP/ASCO/CAP tiering, CIViC evidence levels.
4. **RAG architecture** -- query embedding, multi-collection retrieval, weighted
   re-ranking, prompt construction with grounded citations.
5. **FastAPI + Streamlit** -- lifespan management, dependency injection,
   Pydantic request/response models, Streamlit session state.

---

## Codebase Map

```
agent/
├── config/
│   └── settings.py               OncoSettings (Pydantic BaseSettings, ONCO_ env prefix)
├── src/                           11,440 lines across 13+ modules
│   ├── agent.py          (553)    OncoIntelligenceAgent -- plan/search/evaluate/synthesize
│   ├── case_manager.py   (516)    VCF parsing, case lifecycle, MTB packet generation
│   ├── collections.py    (665)    11 collection schemas, OncoCollectionManager
│   ├── cross_modal.py    (383)    Cross-agent triggers (genomic + imaging enrichment)
│   ├── export.py       (1,055)    Markdown / JSON / PDF / FHIR R4 export
│   ├── knowledge.py    (1,662)    Knowledge graphs (targets, therapies, resistance, pathways)
│   ├── metrics.py        (362)    Prometheus instrumentation
│   ├── models.py         (538)    14 Pydantic models, 13 enums
│   ├── query_expansion.py(812)    Domain-aware query rewriting
│   ├── rag_engine.py     (908)    Multi-collection RAG with comparative retrieval
│   ├── scheduler.py               Background task scheduling
│   ├── therapy_ranker.py (748)    Evidence-based therapy ranking
│   ├── trial_matcher.py  (513)    Hybrid deterministic + semantic trial matching
│   ├── ingest/         (1,793)    9 parsers + abstract base class
│   │   ├── base.py                BaseIngestPipeline (fetch -> parse -> embed_and_store)
│   │   ├── civic_parser.py        CIViC variant evidence
│   │   ├── oncokb_parser.py       OncoKB annotations
│   │   ├── literature_parser.py   PubMed / PMC abstracts
│   │   ├── clinical_trials_parser.py  ClinicalTrials.gov v2 API
│   │   ├── guideline_parser.py    NCCN / ESMO / ASCO guidelines
│   │   ├── resistance_parser.py   Resistance mechanism seeds
│   │   ├── pathway_parser.py      Signaling pathway seeds
│   │   └── outcome_parser.py      Real-world outcome records
│   └── utils/            (668)    VCF parser, PubMed client
├── api/                 (1,300)
│   ├── main.py          (410)     FastAPI lifespan, CORS, health endpoint
│   └── routes/                    5 routers (meta_agent, cases, trials, reports, events)
├── app/                   (758)
│   └── oncology_ui.py             Streamlit 5-tab workbench
├── tests/               (4,370)   556 tests across 10 files + conftest
│   ├── test_models.py   (644)     Parametrized enum validation
│   ├── test_integration.py(785)   4 patient profiles, full pipeline
│   ├── test_knowledge.py(603)     Knowledge graph integrity
│   ├── test_export.py   (439)     All 4 export formats
│   ├── test_therapy_ranker.py(363) Ranking logic
│   ├── test_trial_matcher.py(363) Matching logic
│   ├── test_case_manager.py(332)  VCF and case lifecycle
│   ├── test_rag_engine.py(301)    Retrieval pipeline
│   ├── test_collections.py(276)   Schema validation
│   └── test_agent.py    (264)     Planning and evaluation
└── scripts/             (2,273)   17 scripts (ingest, seed, validate, benchmark)
```

---

## Data Models Overview (src/models.py)

The agent defines 13 enums and 14 Pydantic models that enforce type safety
across all modules:

### Enums (13 total)

| Enum              | Values                                               |
|-------------------|------------------------------------------------------|
| CancerType        | 25 types: NSCLC, SCLC, BREAST, COLORECTAL, ...      |
| VariantType       | SNV, INDEL, CNV_AMP, CNV_DEL, FUSION, REARRANGEMENT, SV |
| EvidenceLevel     | A (FDA-approved), B (clinical), C (case), D (preclinical), E (inferential) |
| TherapyCategory   | TARGETED, IMMUNOTHERAPY, CHEMOTHERAPY, HORMONAL, COMBINATION, RADIOTHERAPY, CELL_THERAPY, ADC, BISPECIFIC |
| TrialPhase        | Early Phase 1, Phase 1, Phase 1/2, Phase 2, Phase 2/3, Phase 3, Phase 4, N/A |
| TrialStatus       | 9 values: Recruiting, Active not recruiting, Completed, Terminated, ... |
| ResponseCategory  | CR (complete), PR (partial), SD (stable), PD (progressive), NE (not evaluable) |
| BiomarkerType     | PREDICTIVE, PROGNOSTIC, DIAGNOSTIC, MONITORING, RESISTANCE, PHARMACODYNAMIC, SCREENING, THERAPEUTIC_SELECTION |
| PathwayName       | 13 pathways: MAPK, PI3K_AKT_MTOR, DDR, CELL_CYCLE, APOPTOSIS, WNT, NOTCH, HEDGEHOG, JAK_STAT, ANGIOGENESIS, HIPPO, NF_KB, TGF_BETA |
| GuidelineOrg      | NCCN, ESMO, ASCO, WHO, CAP_AMP, FDA, EMA, AACR     |
| SourceType        | PUBMED, PMC, PREPRINT, MANUAL                       |

### Domain Models (14 total)

Each model implements `to_embedding_text()` to generate the text string
that gets embedded for vector storage. This method concatenates the most
semantically relevant fields with pipe separators:

| Model                | Key Fields                                           |
|----------------------|------------------------------------------------------|
| OncologyLiterature   | id, title, text_chunk, source_type, year, gene       |
| OncologyTrial        | id (NCT), title, phase, status, biomarker_criteria   |
| OncologyVariant      | id, gene, variant_name, evidence_level, drugs        |
| OncologyBiomarker    | id, name, biomarker_type, testing_method, cutoff     |
| OncologyTherapy      | id, drug_name, category, targets, mechanism_of_action|
| OncologyPathway      | id, name, key_genes, therapeutic_targets, cross_talk |
| OncologyGuideline    | id, org, cancer_type, version, year, recommendations |
| ResistanceMechanism  | id, primary_therapy, gene, mechanism, alternatives   |
| OutcomeRecord        | id, case_id, therapy, response, duration_months      |
| CaseSnapshot         | case_id, patient_id, cancer_type, variants, biomarkers|
| MTBPacket            | case_id, variant_table, therapy_ranking, trial_matches|
| AgentQuery           | question, gene, cancer_type, filters                 |
| SearchHit            | collection, id, score, text, metadata, label         |
| CrossCollectionResult| query, hits, total_collections_searched              |
| AgentResponse        | question, answer, evidence, plan, report             |

### Example: OncologyVariant.to_embedding_text()

```python
def to_embedding_text(self) -> str:
    parts = [
        f"{self.gene} {self.variant_name}",
        self.text_summary,
        f"Type: {self.variant_type.value}",
        f"Evidence: {self.evidence_level.value}",
    ]
    if self.cancer_type:
        parts.append(f"Cancer: {self.cancer_type.value}")
    if self.drugs:
        parts.append(f"Drugs: {', '.join(self.drugs)}")
    if self.clinical_significance:
        parts.append(f"Significance: {self.clinical_significance}")
    return " | ".join(parts)
```

This produces strings like:
`"EGFR L858R | Sensitizing mutation... | Type: snv | Evidence: A | Cancer: nsclc | Drugs: osimertinib, erlotinib"`

---

# Chapter 1: Deep Dive into the RAG Engine

File: `src/rag_engine.py` (908 lines)

The `OncoRAGEngine` is the central retrieval component. Every query passes
through it -- whether from the agent, the API, or the Streamlit UI.

## 1.1 Collection Configuration

The engine searches 11 Milvus collections in parallel. Each collection has a
configurable weight that scales raw cosine-similarity scores before merging:

```python
COLLECTION_CONFIG: Dict[str, Dict[str, Any]] = {
    "onco_variants":   {"weight": settings.WEIGHT_VARIANTS,   "label": "Variant",    "filter_field": "gene",  "year_field": None},
    "onco_literature": {"weight": settings.WEIGHT_LITERATURE,  "label": "Literature", "filter_field": "gene",  "year_field": "year"},
    "onco_therapies":  {"weight": settings.WEIGHT_THERAPIES,   "label": "Therapy",    "filter_field": None,    "year_field": None},
    "onco_guidelines": {"weight": settings.WEIGHT_GUIDELINES,  "label": "Guideline",  "filter_field": None,    "year_field": "year"},
    "onco_trials":     {"weight": settings.WEIGHT_TRIALS,      "label": "Trial",      "filter_field": None,    "year_field": "start_year"},
    "onco_biomarkers": {"weight": settings.WEIGHT_BIOMARKERS,  "label": "Biomarker",  "filter_field": None,    "year_field": None},
    "onco_resistance": {"weight": settings.WEIGHT_RESISTANCE,  "label": "Resistance", "filter_field": "gene",  "year_field": None},
    "onco_pathways":   {"weight": settings.WEIGHT_PATHWAYS,    "label": "Pathway",    "filter_field": None,    "year_field": None},
    "onco_outcomes":   {"weight": settings.WEIGHT_OUTCOMES,    "label": "Outcome",    "filter_field": None,    "year_field": None},
    "onco_cases":      {"weight": settings.WEIGHT_CASES,       "label": "Case",       "filter_field": None,    "year_field": None},
    "genomic_evidence":{"weight": settings.WEIGHT_GENOMIC,     "label": "Genomic",    "filter_field": None,    "year_field": None},
}
```

Default weights (from `config/settings.py`) sum to 1.0:

| Collection        | Default Weight | Purpose                            |
|-------------------|---------------:|------------------------------------|
| onco_variants     |           0.18 | CIViC / OncoKB variant evidence    |
| onco_literature   |           0.16 | PubMed / PMC literature chunks     |
| onco_therapies    |           0.14 | Approved & investigational drugs   |
| onco_guidelines   |           0.12 | NCCN / ESMO / ASCO recommendations |
| onco_trials       |           0.10 | ClinicalTrials.gov summaries       |
| onco_biomarkers   |           0.08 | Predictive / prognostic biomarkers |
| onco_resistance   |           0.07 | Resistance mechanisms              |
| onco_pathways     |           0.06 | Signaling pathway context          |
| onco_outcomes     |           0.04 | Real-world treatment outcomes      |
| genomic_evidence  |           0.03 | VCF-derived evidence (Stage 1)     |
| onco_cases        |           0.02 | De-identified patient snapshots    |

**Key insight:** `filter_field` and `year_field` allow the engine to apply
Milvus metadata filters (gene-level narrowing, publication-year ranges) on a
per-collection basis. Collections without a `filter_field` rely purely on
vector similarity.

## 1.2 The retrieve() Pipeline

The `retrieve()` method is the heart of the engine. Here is the step-by-step
execution flow:

```
retrieve(query, top_k, collections_filter, year_min, year_max, conversation_context)
  │
  ├─ 1. Embed query text
  │     └─ _embed_query(): prepends BGE instruction prefix
  │        "Represent this sentence for searching relevant passages: <query>"
  │
  ├─ 2. _search_all_collections() -- parallel ThreadPoolExecutor
  │     ├─ For each target collection (up to 8 workers):
  │     │   ├─ Build Milvus filters (gene, year range)
  │     │   ├─ collection_manager.search(collection, vector, top_k, filters)
  │     │   ├─ Scale raw _distance by collection weight
  │     │   ├─ Wrap in SearchHit with label, citation, relevance
  │     │   └─ Return hits
  │     └─ Merge all hits into a flat list
  │
  ├─ 3. Expanded search (if query_expander provided)
  │     ├─ query_expander(query.text) -> expansion terms
  │     ├─ Concatenate: "{query} {expansion_terms}"
  │     ├─ Re-embed expanded text
  │     └─ _search_all_collections() with top_k // 2
  │
  ├─ 4. _merge_and_rank()
  │     ├─ Deduplicate by record_id (first-seen wins)
  │     ├─ Sort descending by weighted score
  │     └─ Cap at _MAX_EVIDENCE = 30
  │
  └─ 5. Return CrossCollectionResult(query, hits, total_collections_searched)
```

## 1.3 Relevance Scoring

After weighting, each hit receives a human-readable relevance label:

```python
@staticmethod
def _score_relevance(score: float) -> str:
    if score >= 0.85:
        return "high"
    if score >= 0.65:
        return "medium"
    return "low"
```

These labels appear in the prompt sent to the LLM so it can calibrate
confidence in individual evidence items.

## 1.4 Citation Formatting

The engine automatically generates clickable citation links:

- **PubMed IDs** (`PMID:12345`) -> `[PubMed 12345](https://pubmed.ncbi.nlm.nih.gov/12345/)`
- **NCT IDs** (`NCT01234567`) -> `[NCT01234567](https://clinicaltrials.gov/study/NCT01234567)`
- **Everything else** -> `[Label: record_id]` using the collection label

## 1.5 Comparative Retrieval

When the query contains comparative language (`vs`, `versus`, `compare`,
`difference between`, `head to head`), the engine switches to a dual-entity
retrieval path:

```
_is_comparative(question)
  └─ regex: r"\b(compare|vs\.?|versus|difference between|head.to.head)\b"

retrieve_comparative(question, ...)
  ├─ _parse_comparison_entities() -> (entity_a, entity_b)
  │   Handles: "A vs B", "compare A and B", "difference between A and B"
  ├─ retrieve(entity_a) -> hits_a
  ├─ retrieve(entity_b) -> hits_b
  ├─ Compute shared_hits = intersection by record_id
  └─ Return {entity_a, entity_b, hits_a, hits_b, shared_hits}
```

The comparative prompt template instructs the LLM to structure its answer
across 8 comparison dimensions: mechanism of action, efficacy data, safety
profile, biomarker considerations, resistance mechanisms, guideline
recommendations, clinical trial evidence, and summary recommendation.

## 1.6 Knowledge Context Injection

Before building the LLM prompt, `_get_knowledge_context()` queries the
knowledge module for five types of domain knowledge:

1. Gene mentions -> `knowledge.lookup_gene()`
2. Therapy mentions -> `knowledge.lookup_therapy()`
3. Resistance mentions -> `knowledge.lookup_resistance()`
4. Pathway mentions -> `knowledge.lookup_pathway()`
5. Biomarker mentions -> `knowledge.lookup_biomarker()`

Each successful lookup is tagged with a section header (`[Gene Knowledge]`,
`[Therapy Knowledge]`, etc.) and prepended to the evidence in the prompt.

## 1.7 Prompt Construction

The final prompt has this structure:

```
=== Domain Knowledge ===
[Gene Knowledge] ...
[Therapy Knowledge] ...

=== Retrieved Evidence ===
1. Variant [high] (score 0.892) -- [PubMed 12345](...)
   EGFR L858R confers sensitivity to osimertinib...
2. Literature [medium] (score 0.734) -- [Variant: civic-123]
   ...

=== Question ===
<original question>

Using the evidence above, provide a thorough, well-cited answer...
```

## 1.8 System Prompt

The system prompt establishes 8 core competency areas:

1. Molecular profiling
2. Variant interpretation (CIViC/OncoKB evidence levels, AMP/ASCO/CAP)
3. Therapy selection (NCCN/ESMO guideline-concordant)
4. Clinical trial matching
5. Resistance mechanisms
6. Biomarker assessment
7. Outcomes monitoring (RECIST, MRD, ctDNA)
8. Cross-modal integration (imaging, drug discovery pipelines)

Five behavioral instructions enforce citation, cross-functional reasoning,
resistance/contraindication awareness, guideline references, and uncertainty
acknowledgment.

---

# Chapter 2: The OncoIntelligenceAgent

File: `src/agent.py` (553 lines)

The agent implements the **plan-search-evaluate-synthesize** loop -- the
highest-level orchestration pattern in the system.

## 2.1 SearchPlan Dataclass

Every query starts with a structured plan:

```python
@dataclass
class SearchPlan:
    question: str
    identified_topics: List[str] = field(default_factory=list)
    target_genes: List[str] = field(default_factory=list)
    relevant_cancer_types: List[str] = field(default_factory=list)
    search_strategy: str = "broad"   # "broad" | "targeted" | "comparative"
    sub_questions: List[str] = field(default_factory=list)
```

## 2.2 Gene and Cancer-Type Recognition

The planner uses two static vocabularies for entity extraction:

- **KNOWN_GENES** (30 entries): `BRAF`, `EGFR`, `ALK`, `ROS1`, `KRAS`,
  `HER2`, `NTRK`, `RET`, `MET`, `FGFR`, `PIK3CA`, `IDH1`, `IDH2`, `BRCA`,
  `BRCA1`, `BRCA2`, `TP53`, `PTEN`, `CDKN2A`, `STK11`, `ESR1`, `ERBB2`,
  `NRAS`, `APC`, `VHL`, `KIT`, `PDGFRA`, `FLT3`, `NPM1`, `DNMT3A`

- **KNOWN_CANCER_TYPES** (25 entries): `NSCLC`, `BREAST`, `MELANOMA`,
  `COLORECTAL`, `PANCREATIC`, `OVARIAN`, `PROSTATE`, `GLIOMA`,
  `GLIOBLASTOMA`, `AML`, `CML`, `CLL`, `DLBCL`, `BLADDER`, `RENAL`,
  `HEPATOCELLULAR`, `GASTRIC`, `ESOPHAGEAL`, `THYROID`, `ENDOMETRIAL`,
  `CERVICAL`, `HEAD_AND_NECK`, `SARCOMA`, `CHOLANGIOCARCINOMA`,
  `MESOTHELIOMA`

- **_CANCER_ALIASES** (50+ entries): maps natural language like
  `"lung cancer"` -> `"NSCLC"`, `"triple negative breast"` -> `"BREAST"`,
  `"gbm"` -> `"GLIOBLASTOMA"`, `"crpc"` -> `"PROSTATE"`, etc.

Gene detection is case-insensitive uppercase matching against `q_upper`.
Cancer type detection checks both canonical names and aliases.

## 2.3 Topic Detection

The planner scans for 20+ topic keywords and maps them to clinical concepts:

| Keyword            | Topic                          |
|--------------------|--------------------------------|
| `resistance`       | therapeutic resistance         |
| `biomarker`        | biomarker identification       |
| `immunotherapy`    | immunotherapy response         |
| `combination`      | combination therapy            |
| `tmb`              | TMB                            |
| `msi`              | MSI / microsatellite instab.   |
| `ctdna`            | liquid biopsy / ctDNA          |
| `pdl1`, `pd-l1`    | PD-L1 / immune checkpoint      |
| `fusion`           | gene fusion                    |
| `methylation`      | epigenetic regulation          |

## 2.4 Strategy Selection

```python
if any(sig in q_lower for sig in comparative_signals):
    search_strategy = "comparative"
elif target_genes and relevant_cancer_types:
    search_strategy = "targeted"
else:
    search_strategy = "broad"
```

- **comparative**: triggered by `compare`, `vs`, `versus`,
  `difference between`, `head to head`
- **targeted**: both gene(s) and cancer type(s) identified
- **broad**: fallback when context is ambiguous

## 2.5 Question Decomposition

Complex queries are broken into focused sub-questions:

- Multiple genes -> one sub-question per gene
  (e.g., "What is the role of EGFR in NSCLC?")
- Multiple cancer types -> one sub-question per type
  (e.g., "BRAF therapeutic landscape in melanoma")
- Topic-driven sub-questions for resistance, trials, biomarkers, combinations

## 2.6 The run() Pipeline

```python
def run(self, question: str, **kwargs) -> AgentResponse:
    # 1. Plan
    plan = self.search_plan(question)

    # 2. Search with adaptive retry
    for attempt in range(1, MAX_RETRIES + 2):     # up to 3 attempts
        for q in [plan.question] + plan.sub_questions:
            results = self.rag_engine.cross_collection_search(AgentQuery(question=q))
            all_evidence.extend(results)

        # 3. Evaluate
        verdict = self.evaluate_evidence(all_evidence)
        if verdict == "sufficient" or attempt > MAX_RETRIES:
            break
        # Broaden: targeted -> broad, generate fallback queries
        queries_to_run = self._generate_fallback_queries(plan)

    # 4. Synthesize
    response = self.rag_engine.synthesize(question, all_evidence, plan)
    response.report = self.generate_report(response)
    return response
```

## 2.7 Evidence Evaluation

The `evaluate_evidence()` method classifies evidence adequacy as one of three
verdicts:

| Verdict        | Criteria                                                     |
|----------------|--------------------------------------------------------------|
| `sufficient`   | >= 3 hits AND >= 2 collections represented                   |
| `partial`      | > 0 hits but insufficient diversity or count                 |
| `insufficient` | 0 usable hits (all below MIN_SIMILARITY_SCORE = 0.30)        |

Evidence items with scores below 0.30 are filtered out before evaluation.
An average score >= 0.50 is preferred but not required for `sufficient`.

## 2.8 Adaptive Retry

When evidence is insufficient:

1. Switch strategy from `targeted` to `broad`
2. Generate fallback queries:
   - Per gene: `"{gene} oncology therapeutic implications"`,
     `"{gene} mutation clinical significance"`
   - Per cancer type: `"{ct} precision medicine current landscape"`
   - Default: `"{question} precision oncology"`

Maximum retries: 2 (total 3 attempts including the initial search).

## 2.9 Report Generation

The agent produces a structured Markdown report with sections:
- Query (original question)
- Analysis (strategy, genes, cancer types, topics, sub-questions)
- Evidence Sources (grouped by collection, top 10 per collection)
- Knowledge Graph (if attached)
- Synthesis (LLM-generated answer)

---

# Chapter 3: Knowledge Graph Architecture

File: `src/knowledge.py` (1,662 lines)

The knowledge module is a curated, code-embedded domain knowledge graph
providing instant lookup without vector search latency.

## 3.1 ACTIONABLE_TARGETS (~40 genes)

Each entry has a consistent structure:

```python
ACTIONABLE_TARGETS["EGFR"] = {
    "gene": "EGFR",
    "full_name": "Epidermal Growth Factor Receptor",
    "cancer_types": ["NSCLC", "head and neck", "colorectal", "glioblastoma"],
    "key_variants": ["L858R", "exon 19 deletion", "T790M", "C797S",
                     "exon 20 insertion", "S768I", "L861Q", "G719X"],
    "targeted_therapies": ["osimertinib", "erlotinib", "gefitinib",
                           "afatinib", "dacomitinib", "amivantamab"],
    "combination_therapies": ["osimertinib + chemotherapy",
                              "amivantamab + lazertinib"],
    "resistance_mutations": ["T790M", "C797S", "MET amplification",
                             "HER2 amplification", "small cell transformation",
                             "BRAF V600E", "PIK3CA mutations"],
    "pathway": "MAPK",
    "evidence_level": "A",
    "description": "EGFR mutations occur in ~15-20% of NSCLC ..."
}
```

**Covered genes include:** BRAF, EGFR, ALK, ROS1, KRAS, HER2, NTRK, RET,
MET, FGFR, PIK3CA, IDH1, IDH2, and many more -- each with full_name,
cancer_types, key_variants, targeted_therapies, combination_therapies,
resistance_mutations, pathway, evidence_level, and description.

## 3.2 THERAPY_MAP

Maps drug names (lowercase keys) to structured therapy records:

```python
THERAPY_MAP["osimertinib"] = {
    "drug_name": "osimertinib",
    "brand_name": "Tagrisso",
    "category": "targeted therapy",
    "targets": ["EGFR"],
    "approved_indications": ["EGFR-mutant NSCLC (first-line)", ...],
    "mechanism": "3rd-gen EGFR TKI, active against T790M",
    "key_trials": ["FLAURA", "ADAURA", "FLAURA2"],
}
```

## 3.3 RESISTANCE_MAP

Maps primary therapies to their known resistance mechanisms:

```python
RESISTANCE_MAP["osimertinib"] = {
    "primary_therapy": "osimertinib",
    "resistance_triggers": [...],
    "mechanism": "C797S mutation, MET amplification, ...",
    "alternatives": ["amivantamab + lazertinib", ...],
}
```

The therapy ranker and case manager both consult this map to flag resistance
concerns and suggest next-line options.

## 3.4 PATHWAY_MAP

Maps signaling pathway names to their constituent genes, druggable targets,
and cross-talk connections:

- **MAPK**: BRAF, KRAS, NRAS, MEK1/2, ERK1/2
- **PI3K/AKT/mTOR**: PIK3CA, PTEN, AKT1, mTOR
- **DNA Damage Repair**: BRCA1, BRCA2, ATM, ATR, PALB2
- **Cell Cycle**: CDK4/6, CDKN2A, RB1, CCND1

Each pathway entry includes `cross_talk` describing how pathways interact
(e.g., MAPK <-> PI3K bypass signaling).

## 3.5 BIOMARKER_PANELS

Defines clinically validated biomarker panels with testing methods, cutoffs,
and associated therapies:

| Biomarker  | Type       | Testing Method    | Clinical Cutoff     | Evidence |
|------------|------------|-------------------|---------------------|----------|
| TMB        | Predictive | WGS / WES / Panel | >= 10 mut/Mb        | A        |
| MSI        | Predictive | IHC / PCR / NGS   | MSI-H               | A        |
| PD-L1 TPS  | Predictive | IHC (22C3/SP263)  | >= 50% (first-line) | A        |
| HRD        | Predictive | Myriad MyChoice   | HRD score >= 42     | A        |

## 3.6 ENTITY_ALIASES

Provides 30+ alias mappings for entity resolution in natural language queries:

```python
ENTITY_ALIASES = {
    "keytruda": "pembrolizumab",
    "opdivo": "nivolumab",
    "tagrisso": "osimertinib",
    "herceptin": "trastuzumab",
    ...
}
```

## 3.7 Lookup Functions

The module exposes helper functions consumed by the RAG engine:

- `lookup_gene(query)` -- searches ACTIONABLE_TARGETS for gene mentions
- `lookup_therapy(query)` -- searches THERAPY_MAP
- `lookup_resistance(query)` -- searches RESISTANCE_MAP
- `lookup_pathway(query)` -- searches PATHWAY_MAP
- `lookup_biomarker(query)` -- searches BIOMARKER_PANELS
- `get_target_context(gene)` -- formatted context string for a gene
- `classify_variant_actionability(gene, variant)` -- returns evidence level
  (A/B/C/D/VUS)

---

# Chapter 4: Query Expansion and Rewriting

File: `src/query_expansion.py` (812 lines)

## 4.1 Expansion Categories

The module contains 12 domain-specific expansion dictionaries:

1. **CANCER_TYPE_EXPANSIONS** -- maps abbreviations to full names and subtypes
   (e.g., `"NSCLC"` -> `["non-small cell lung cancer", "lung adenocarcinoma",
   "lung squamous cell", "EGFR-mutant lung", "ALK-positive lung"]`)

2. **GENE_EXPANSIONS** -- maps gene symbols to full names and common variants
   (e.g., `"EGFR"` -> `["epidermal growth factor receptor", "EGFR L858R",
   "EGFR exon 19 deletion", "EGFR T790M", "EGFR C797S"]`)

3. **THERAPY_EXPANSIONS** -- maps drug names to brand names, mechanisms,
   and related compounds

4. **BIOMARKER_EXPANSIONS** -- maps biomarker abbreviations to full terms

5. **PATHWAY_EXPANSIONS** -- signaling pathway synonyms

6. **RESISTANCE_EXPANSIONS** -- resistance mechanism terminology

7. **CLINICAL_TERM_EXPANSIONS** -- clinical outcomes and staging terms

8. **TRIAL_EXPANSIONS** -- trial-related terminology

9. **IMMUNOTHERAPY_EXPANSIONS** -- checkpoint inhibitor vocabulary

10. **SURGERY_RADIATION_EXPANSIONS** -- procedural terms

11. **TOXICITY_EXPANSIONS** -- adverse event terminology

12. **GENOMICS_EXPANSIONS** -- sequencing and variant calling terms

## 4.2 How Expansion Works

The RAG engine calls the expander as a callable:

```python
expansion_terms = self.query_expander(query.text)
expanded_text = f"{query.text} {' '.join(expansion_terms)}"
expanded_vector = self._embed_query(expanded_text)
```

The expander scans the input query for keywords matching any expansion
dictionary. For each match, it appends the associated expansion terms.
This widens the embedding to capture semantically related documents that
might use different terminology (e.g., a query mentioning "osimertinib"
also captures "Tagrisso" and "3rd-generation EGFR TKI").

## 4.3 Expansion Coverage

| Category         | Key Count | Example Key  | Example Expansions                              |
|------------------|-----------|--------------|-------------------------------------------------|
| Cancer types     |        16 | `NSCLC`      | non-small cell lung cancer, lung adenocarcinoma  |
| Genes            |        12 | `KRAS`       | KRAS G12C, KRAS G12D, Kirsten rat sarcoma        |
| Therapies        |       ~15 | `osimertinib`| Tagrisso, 3rd-gen EGFR TKI, T790M active         |
| Biomarkers       |        ~8 | `TMB`        | tumor mutational burden, TMB-high, TMB >= 10     |
| Pathways         |        ~6 | `MAPK`       | RAS/RAF/MEK/ERK, MAP kinase cascade              |
| Resistance       |        ~8 | `T790M`      | gatekeeper mutation, acquired resistance EGFR    |
| Clinical terms   |       ~10 | `PFS`        | progression-free survival, time to progression   |
| Trials           |        ~6 | `Phase 3`    | pivotal trial, registration study                |
| Immunotherapy    |        ~8 | `checkpoint`  | PD-1, PD-L1, CTLA-4, immune checkpoint           |

---

# Chapter 5: Collection Schemas and Indexing

File: `src/collections.py` (665 lines)

## 5.1 Shared Index Configuration

All 11 collections use identical vector index parameters:

```python
EMBEDDING_DIM = 384  # BGE-small-en-v1.5

INDEX_PARAMS = {
    "metric_type": "COSINE",
    "index_type": "IVF_FLAT",
    "params": {"nlist": 1024},
}

SEARCH_PARAMS = {
    "metric_type": "COSINE",
    "params": {"nprobe": 16},
}
```

**IVF_FLAT** partitions the vector space into 1024 Voronoi cells. At query
time, `nprobe=16` cells are searched (1.6% of partitions), balancing speed
and recall. For larger datasets (>1M vectors per collection), consider
increasing `nprobe` to 32-64 or switching to `IVF_PQ` for memory savings.

## 5.2 Collection Schema Reference

Each collection follows a pattern: VARCHAR primary key `id`, a 384-dim
`FLOAT_VECTOR` embedding field, and typed metadata fields.

### onco_variants

| Field                 | Type         | Max Length | Notes                |
|-----------------------|-------------|-----------|----------------------|
| id                    | VARCHAR (PK) | 100       |                      |
| embedding             | FLOAT_VECTOR | dim=384   |                      |
| gene                  | VARCHAR      | 50        | Filterable           |
| variant_name          | VARCHAR      | 100       |                      |
| variant_type          | VARCHAR      | 30        |                      |
| cancer_type           | VARCHAR      | 50        |                      |
| evidence_level        | VARCHAR      | 20        | A/B/C/D/E            |
| drugs                 | VARCHAR      | 500       |                      |
| civic_id              | VARCHAR      | 20        |                      |
| vrs_id                | VARCHAR      | 100       | VRS identifier       |
| text_summary          | VARCHAR      | 3000      |                      |
| clinical_significance | VARCHAR      | 200       |                      |
| allele_frequency      | FLOAT        |           |                      |

### onco_literature

| Field       | Type         | Max Length | Notes        |
|-------------|-------------|-----------|--------------|
| id          | VARCHAR (PK) | 100       |              |
| embedding   | FLOAT_VECTOR | dim=384   |              |
| title       | VARCHAR      | 500       |              |
| text_chunk  | VARCHAR      | 3000      |              |
| source_type | VARCHAR      | 20        | pubmed/pmc   |
| year        | INT64        |           | Filterable   |
| cancer_type | VARCHAR      | 50        |              |
| gene        | VARCHAR      | 50        | Filterable   |
| variant     | VARCHAR      | 100       |              |
| keywords    | VARCHAR      | 1000      |              |
| journal     | VARCHAR      | 200       |              |

### onco_trials

| Field              | Type         | Max Length | Notes        |
|--------------------|-------------|-----------|--------------|
| id                 | VARCHAR (PK) | 20        | NCT ID       |
| embedding          | FLOAT_VECTOR | dim=384   |              |
| title              | VARCHAR      | 500       |              |
| text_summary       | VARCHAR      | 3000      |              |
| phase              | VARCHAR      | 30        |              |
| status             | VARCHAR      | 30        |              |
| sponsor            | VARCHAR      | 200       |              |
| cancer_types       | VARCHAR      | 200       |              |
| biomarker_criteria | VARCHAR      | 500       |              |
| enrollment         | INT64        |           |              |
| start_year         | INT64        |           | Filterable   |
| outcome_summary    | VARCHAR      | 2000      |              |

### genomic_evidence (read-only)

| Field                 | Type         | Max Length | Notes                   |
|-----------------------|-------------|-----------|-------------------------|
| id                    | VARCHAR (PK) | 200       |                         |
| embedding             | FLOAT_VECTOR | dim=384   |                         |
| chrom                 | VARCHAR      | 10        |                         |
| pos                   | INT64        |           |                         |
| ref                   | VARCHAR      | 500       |                         |
| alt                   | VARCHAR      | 500       |                         |
| qual                  | FLOAT        |           |                         |
| gene                  | VARCHAR      | 50        |                         |
| consequence           | VARCHAR      | 100       |                         |
| impact                | VARCHAR      | 20        | HIGH/MODERATE/LOW/MOD   |
| genotype              | VARCHAR      | 10        | 0/1, 1/1, etc.         |
| text_summary          | VARCHAR      | 2000      |                         |
| clinical_significance | VARCHAR      | 200       |                         |
| rsid                  | VARCHAR      | 20        |                         |
| disease_associations  | VARCHAR      | 500       |                         |
| am_pathogenicity      | FLOAT        |           | AlphaMissense score     |
| am_class              | VARCHAR      | 30        | likely_pathogenic, etc. |

## 5.3 Schema and Model Registries

Two registries map collection names to schemas and Pydantic models:

```python
COLLECTION_SCHEMAS: Dict[str, CollectionSchema] = {
    "onco_literature": ONCO_LITERATURE_SCHEMA,
    "onco_trials": ONCO_TRIALS_SCHEMA,
    # ... all 11 collections
}

COLLECTION_MODELS: Dict[str, Optional[Type]] = {
    "onco_literature": OncologyLiterature,
    "onco_trials": OncologyTrial,
    # ...
    "genomic_evidence": None,  # read-only, populated by Stage 1
}
```

## 5.4 OncoCollectionManager

The manager wraps pymilvus operations:

- `connect()` / `disconnect()` -- Milvus connection lifecycle
- `create_collection(name)` -- creates from COLLECTION_SCHEMAS registry
- `get_collection(name)` -- returns cached `Collection` handle
- `get_collection_count(name)` -- entity count after flush
- `insert(collection_name, data)` -- single-record or batch insert
- `search(collection_name, query_vector, top_k, filters)` -- ANN search
- `query(collection_name, filter_expr, output_fields, limit)` -- filter query

Parallel search across all collections uses `ThreadPoolExecutor` with
`max_workers = min(len(collections), 8)`.

---

# Chapter 6: Therapy Ranking Engine

File: `src/therapy_ranker.py` (748 lines)

## 6.1 Evidence Level Ordering

```python
EVIDENCE_LEVEL_ORDER = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "VUS": 5}
```

Lower ordinal = stronger evidence. Level A (FDA-approved / companion
diagnostic) always ranks above Level B (clinical evidence), and so on.

## 6.2 The rank_therapies() Pipeline

```
rank_therapies(cancer_type, variants, biomarkers, prior_therapies)
  │
  ├─ Step 1: Identify variant-driven therapies
  │   └─ For each variant: _identify_variant_therapies(gene, variant, cancer_type)
  │       ├─ Look up gene in ACTIONABLE_TARGETS
  │       ├─ Check if variant matches key_variants
  │       ├─ Get drug list from targeted_therapies
  │       └─ Return therapy dict with evidence_level
  │
  ├─ Step 2: Identify biomarker-driven therapies
  │   └─ _identify_biomarker_therapies(biomarkers, cancer_type)
  │       ├─ MSI-H/dMMR -> pembrolizumab, nivolumab, dostarlimab (Level A)
  │       ├─ TMB >= 10 -> pembrolizumab (Level A), atezolizumab (Level B)
  │       ├─ HRD/BRCA -> olaparib, rucaparib, niraparib, talazoparib
  │       ├─ PTEN loss -> alpelisib (Level C)
  │       ├─ PD-L1 TPS >= 50% -> pembrolizumab first-line (Level A)
  │       ├─ NTRK fusion -> larotrectinib, entrectinib (Level A)
  │       └─ BIOMARKER_PANELS for additional mappings
  │
  ├─ Step 3: Deduplicate (keep strongest evidence per drug)
  │   └─ Sort by EVIDENCE_LEVEL_ORDER
  │
  ├─ Step 4: Check resistance via _check_resistance()
  │   ├─ Mutation-level: RESISTANCE_MAP lookup
  │   └─ Class-level: _DRUG_CLASS_GROUPS cross-resistance
  │
  ├─ Step 5: Check contraindications via _check_contraindication()
  │   ├─ Direct match: same drug used before
  │   └─ Same drug class via THERAPY_MAP category matching
  │
  ├─ Step 6: Retrieve supporting evidence from Milvus
  │   └─ Search onco_therapies + onco_literature per drug
  │
  ├─ Step 6.5: Identify combination regimens
  │   └─ _COMBO_REGIMENS: 6 FDA-approved combos
  │       dabrafenib+trametinib, encorafenib+binimetinib,
  │       encorafenib+cetuximab, ipilimumab+nivolumab,
  │       lenvatinib+pembrolizumab, trastuzumab+pertuzumab
  │
  └─ Step 7: _assign_final_ranks()
      ├─ Partition: clean vs. flagged (resistance/contraindication)
      ├─ Sort each group by evidence level
      ├─ Combine: clean first, then flagged
      └─ Assign rank 1..N
```

## 6.3 Drug Class Groups

The ranker defines 10 drug class groups for cross-resistance detection:

```python
_DRUG_CLASS_GROUPS = {
    "egfr_tki_1g":     ["erlotinib", "gefitinib"],
    "egfr_tki_2g":     ["afatinib", "dacomitinib"],
    "egfr_tki_3g":     ["osimertinib"],
    "alk_tki":         ["crizotinib", "ceritinib", "alectinib",
                        "brigatinib", "lorlatinib"],
    "braf_inhibitor":  ["vemurafenib", "dabrafenib", "encorafenib"],
    "mek_inhibitor":   ["trametinib", "cobimetinib", "binimetinib"],
    "anti_pd1":        ["pembrolizumab", "nivolumab", "dostarlimab",
                        "cemiplimab"],
    "anti_pdl1":       ["atezolizumab", "durvalumab", "avelumab"],
    "parp_inhibitor":  ["olaparib", "rucaparib", "niraparib",
                        "talazoparib"],
    "kras_g12c":       ["sotorasib", "adagrasib"],
}
```

If a patient previously received erlotinib (1st-gen EGFR TKI), gefitinib
is flagged for cross-resistance because they share the same drug class.

## 6.4 Resistance Check Logic

Two layers of resistance detection:

1. **Mutation-level** (RESISTANCE_MAP): Checks if the candidate drug has
   documented resistance mutations that overlap with prior therapy history.
   Returns mechanism details and suggested alternatives.

2. **Class-level** (_DRUG_CLASS_GROUPS): If the candidate belongs to the
   same drug class as a prior therapy, flags likely cross-resistance even
   without mutation-specific data.

## 6.5 Convenience API

```python
def rank_for_case(self, case: CaseSnapshot) -> List[Dict]:
    """Rank therapies directly from a CaseSnapshot."""
    return self.rank_therapies(
        cancer_type=case.cancer_type,
        variants=case.variants,
        biomarkers=case.biomarkers or {},
        prior_therapies=case.prior_therapies or [],
    )
```

---

# Chapter 7: Clinical Trial Matching

File: `src/trial_matcher.py` (513 lines)

## 7.1 Two-Stage Matching Strategy

### Stage 1: Deterministic Filter

```python
_deterministic_search(cancer_type, top_k=30)
```

- Resolve cancer type to aliases via `_CANCER_ALIASES` (18 cancer type groups)
- For each alias x each open status: build Milvus filter expression
- Filter on: `cancer_type == "{alias}" and status == "{status}"`
- Open statuses: `Recruiting`, `Active, not recruiting`,
  `Enrolling by invitation`, `Not yet recruiting`
- Input validation via `_SAFE_FILTER_RE = r"^[A-Za-z0-9 _.\-/,]+$"` to
  prevent Milvus injection

### Stage 2: Semantic Search

```python
_semantic_search(query_text, top_k=30)
```

- Build natural-language query: `"{cancer_type} clinical trial stage {stage}
  {marker1} {value1} {marker2} {value2}"`
- Embed and search `onco_trials` collection

### Merge

Results are merged by `trial_id` (union). If a trial appears in both sets,
the semantic score is preserved as `_semantic_score`.

## 7.2 Composite Scoring

```python
composite = (
    0.40 * biomarker_score      # fraction of patient biomarkers in criteria
  + 0.25 * semantic_score       # vector similarity
  + 0.20 * phase_weight         # trial phase
  + 0.15 * status_weight        # recruitment status
) * age_penalty                 # 1.0 or 0.5 if age outside range
```

### Phase Weights

| Phase       | Weight |
|-------------|--------|
| Phase 3     | 1.0    |
| Phase 2/3   | 0.9    |
| Phase 2     | 0.8    |
| Phase 1/2   | 0.7    |
| Phase 1     | 0.6    |
| Phase 4     | 0.5    |

### Status Weights

| Status                   | Weight |
|--------------------------|--------|
| Recruiting               | 1.0    |
| Enrolling by invitation  | 0.8    |
| Active, not recruiting   | 0.6    |
| Not yet recruiting       | 0.4    |

## 7.3 Biomarker Matching

`_score_biomarker_match()` performs case-insensitive fuzzy matching of each
patient biomarker key and value against the combined trial criteria text.
Returns fraction matched (0.0 to 1.0).

## 7.4 Age Penalty

`_compute_age_penalty()` parses age eligibility from criteria text using
regex patterns:

- `"Age >= 18"` / `"minimum age: 18"`
- `"Age <= 75"` / `"maximum age: 75"`
- `"18-75 years"` / `"18 to 75 years"`

Returns 1.0 (no penalty) if age is within range or unspecified; 0.5 if
out of range.

## 7.5 Match Explanation

Each matched trial gets a structured explanation:

```python
{
    "trial_id": "NCT04185831",
    "title": "Phase 3 Study of Osimertinib + Chemo in EGFR-mutant NSCLC",
    "phase": "Phase 3",
    "status": "Recruiting",
    "sponsor": "AstraZeneca",
    "match_score": 0.8234,
    "matched_criteria": ["Cancer type: NSCLC", "EGFR=L858R", "Age 62"],
    "unmatched_criteria": ["TMB=8.5 (not explicitly listed)"],
    "explanation": "Matched: Cancer type: NSCLC, EGFR=L858R, Age 62. ..."
}
```

---

# Chapter 8: Case Management and VCF Parsing

File: `src/case_manager.py` (516 lines)

## 8.1 Case Lifecycle

```
create_case(patient_id, cancer_type, stage, vcf_or_variants, biomarkers, prior_therapies)
  ├─ Parse VCF text (if string) via _parse_vcf_text()
  │   └─ Delegates to src/utils/vcf_parser.py:
  │       parse_vcf_text() -> filter_pass_variants() -> extract gene/consequence
  ├─ Classify actionability per variant
  │   └─ classify_variant_actionability(gene, variant) -> A/B/C/D/VUS
  ├─ Generate case_id (UUID)
  ├─ Build text_summary for embedding
  ├─ Create CaseSnapshot
  └─ _store_case() -> embed summary -> insert into onco_cases collection
```

## 8.2 VCF Parsing Details

The VCF parser (`src/utils/vcf_parser.py`) handles three annotation formats:

1. **SnpEff ANN field**: `ANN=A|missense_variant|MODERATE|EGFR|...`
2. **VEP CSQ field**: `CSQ=A|missense_variant|MODERATE|EGFR|...`
3. **GENE / GENEINFO fields**: `GENE=EGFR` or `GENEINFO=EGFR:1956`

Parsing pipeline:
- `parse_vcf_text()` -- split lines, skip headers, extract CHROM/POS/REF/ALT/FILTER/INFO
- `filter_pass_variants()` -- keep only FILTER == "PASS" or "."
- `extract_gene_from_info()` -- try ANN, then CSQ, then GENE/GENEINFO
- `extract_consequence_from_info()` -- extract SnpEff/VEP consequence term

## 8.3 Variant Actionability Classification

```python
def _classify_variant_actionability(self, gene: str, variant: str) -> str:
    return classify_variant_actionability(gene, variant)
    # from src.knowledge:
    # 1. Is gene in ACTIONABLE_TARGETS? No -> "VUS"
    # 2. Does variant match any key_variants? Yes -> evidence_level (A/B/C)
    # 3. Gene-level actionable? -> default_evidence_level
    # 4. Fallback -> "VUS"
```

## 8.4 MTB Packet Generation

`generate_mtb_packet(case_id_or_snapshot)` assembles 5 sections:

1. **variant_table**: all variants with actionability classification and
   associated drugs from ACTIONABLE_TARGETS

2. **evidence_table**: RAG-retrieved evidence for each actionable variant.
   Queries: `"{gene} {variant} {cancer_type} targeted therapy clinical
   evidence"` across `onco_literature` and `onco_therapies`.

3. **therapy_ranking**: delegates to the TherapyRanker (Chapter 6)

4. **trial_matches**: delegates to the TrialMatcher (Chapter 7)

5. **open_questions**: identifies gaps:
   - VUS variants that may need reclassification
   - Missing biomarker results (TMB, MSI, PD-L1 if not provided)
   - Uncertain evidence items needing tumor board discussion

## 8.5 Case Storage

Cases are stored in `onco_cases` with embedded text summaries:

```python
self.collection_manager.insert(
    collection_name="onco_cases",
    data={
        "id": str(snapshot.case_id)[:100],
        "patient_id": str(snapshot.patient_id)[:100],
        "cancer_type": str(snapshot.cancer_type)[:50],
        "stage": str(snapshot.stage or "")[:20],
        "variants": variants_str[:1000],      # serialized to CSV string
        "biomarkers": biomarkers_str[:1000],   # serialized to CSV string
        "prior_therapies": therapies_str[:500],
        "embedding": embedding,
        "text_summary": summary_text[:3000],
    },
)
```

Note the explicit length truncation -- Milvus VARCHAR fields enforce
`max_length` limits and will reject oversized values.

---

# Chapter 9: Export System

File: `src/export.py` (1,055 lines)

## 9.1 Input Normalization

All four export functions accept `MTBPacket`, `dict`, or `str` via
`_normalise_input()`:

```python
def _normalise_input(mtb_packet_or_response):
    if isinstance(..., dict): return it
    if isinstance(..., str):  try json.loads, else wrap in {"raw_text": ...}
    # Pydantic: try .model_dump(), .dict(), .__dict__
```

## 9.2 Markdown Export

`export_markdown(mtb_packet_or_response, title=None) -> str`

Sections generated:
- Header with timestamp, pipeline name, patient ID, cancer type
- Clinical Summary
- Somatic Variant Profile (Markdown table: Gene | Variant | Type | VAF | Consequence | Tier)
- Biomarker Summary (TMB, MSI, PD-L1 + any additional)
- Evidence Summary (per-gene, per-level, with source citations)
- Therapy Ranking (table: Rank | Therapy | Targets | Evidence | Line | Notes)
- Clinical Trial Matches (NCT ID, title, phase, status, match rationale)
- Pathway Context
- Known Resistance Mechanisms
- Open Questions / Follow-Up
- Disclaimer

## 9.3 JSON Export

`export_json(mtb_packet_or_response) -> dict`

Returns a structured dictionary with:

```python
{
    "meta": {
        "format": "hcls-ai-factory-oncology-report",
        "version": "1.0.0",
        "generated_at": "2026-02-15T...",
        "pipeline": "Oncology Intelligence Agent",
        "author": "HCLS AI Factory",
    },
    "patient_id": ...,
    "cancer_type": ...,
    "variants": [...],
    "biomarkers": {...},
    "evidence": [...],
    "therapy_ranking": [...],
    "clinical_trials": [...],
    "pathways": [...],
    "resistance_mechanisms": [...],
    "open_questions": [...],
}
```

## 9.4 PDF Export

`export_pdf(mtb_packet_or_response, output_path) -> str`

Requires ReportLab. Features:

- **NVIDIA branding**: header bar in RGB (118, 185, 0) with white title
- **Custom styles**: NVTitle (20pt white), NVHeading (14pt dark),
  NVBody (10pt), NVDisclaimer (7pt gray)
- **Page layout**: letter size, 40pt margins
- **Structured tables** for variants, therapies, trials using
  `reportlab.platypus.Table` with `TableStyle`
- **Disclaimer footer** on every page

Brand color is overridable via `ONCO_PDF_BRAND_COLOR_R/G/B` environment
variables.

## 9.5 FHIR R4 Export

`export_fhir_r4(mtb_packet_or_response) -> dict`

Generates a FHIR R4 Bundle resource containing:

| FHIR Resource     | Content                                | Coding System |
|-------------------|----------------------------------------|---------------|
| Patient           | Patient demographics                   |               |
| DiagnosticReport  | Genomic report                         | SNOMED        |
| Observation       | Genetic variant assessment             | LOINC 69548-6 |
| Observation       | Tumor mutation burden                  | LOINC 94076-7 |
| Observation       | Microsatellite instability             | LOINC 81695-9 |

### LOINC Codes Used

```python
FHIR_LOINC_CODES = {
    "genomic_report":          "81247-9",
    "gene_studied":            "48018-6",
    "variant":                 "69548-6",
    "therapeutic_implication":  "51969-4",
    "tumor_mutation_burden":   "94076-7",
    "microsatellite_instability": "81695-9",
}
```

### SNOMED Cancer Codes

22 cancer types are mapped to SNOMED CT codes:

```python
FHIR_SNOMED_CANCER_CODES = {
    "nsclc": ("254637007", "Non-small cell lung cancer"),
    "breast": ("254837009", "Malignant neoplasm of breast"),
    "colorectal": ("363406005", "Malignant tumor of colon"),
    "melanoma": ("372244006", "Malignant melanoma"),
    # ... 18 more
}
```

## 9.6 FHIR R4 Bundle Structure

The generated FHIR Bundle has `type: "collection"` and contains:

```json
{
  "resourceType": "Bundle",
  "id": "<uuid>",
  "type": "collection",
  "timestamp": "2026-02-15T10:30:00Z",
  "entry": [
    {
      "fullUrl": "urn:uuid:<patient-uuid>",
      "resource": {
        "resourceType": "Patient",
        "id": "<patient-uuid>",
        "identifier": [{"system": "urn:hcls-ai-factory:patient", "value": "<patient_id>"}],
        "active": true
      }
    },
    {
      "fullUrl": "urn:uuid:<obs-uuid>",
      "resource": {
        "resourceType": "Observation",
        "status": "final",
        "code": {"coding": [{"system": "http://loinc.org", "code": "69548-6"}]},
        "subject": {"reference": "urn:uuid:<patient-uuid>"},
        "valueCodeableConcept": {
          "coding": [{"system": "http://varnomen.hgvs.org", "code": "EGFR L858R"}]
        },
        "component": [
          {"code": {"text": "Gene"}, "valueString": "EGFR"},
          {"code": {"text": "Consequence"}, "valueString": "missense_variant"},
          {"code": {"text": "VAF"}, "valueQuantity": {"value": 0.35}}
        ]
      }
    },
    {
      "resource": {
        "resourceType": "DiagnosticReport",
        "code": {"coding": [{"system": "http://snomed.info/sct", "code": "254637007"}]},
        "result": [{"reference": "urn:uuid:<obs-uuid>"}]
      }
    }
  ]
}
```

### TMB and MSI Observations

When biomarkers include TMB or MSI values, dedicated Observation resources
are created:

- TMB: LOINC `94076-7`, `valueQuantity` with unit `1/1000000{Base}`
- MSI: LOINC `81695-9`, `valueCodeableConcept` with MSI-H/MSS/MSI-L

## 9.7 Export Format Comparison

| Feature             | Markdown | JSON  | PDF   | FHIR R4 |
|---------------------|----------|-------|-------|---------|
| Human-readable      | Yes      | No    | Yes   | No      |
| Machine-parseable   | Partial  | Yes   | No    | Yes     |
| Interoperable       | No       | No    | No    | Yes     |
| Print-ready         | No       | No    | Yes   | No      |
| Branding            | No       | No    | NVIDIA| No      |
| ReportLab required  | No       | No    | Yes   | No      |
| SNOMED/LOINC coded  | No       | No    | No    | Yes     |

---

# Chapter 10: Ingest Pipeline Architecture

Directory: `src/ingest/` (1,793 lines, 9 parsers + base)

## 10.1 BaseIngestPipeline

All parsers inherit from `BaseIngestPipeline` which provides the standard
three-step orchestration:

```python
class BaseIngestPipeline(ABC):
    def __init__(self, collection_manager, embedder, collection_name, batch_size=50):
        ...

    def run(self, query=None, max_results=None) -> int:
        raw_data = self.fetch(**kwargs)            # Step 1: Fetch
        parsed_records = self.parse(raw_data)      # Step 2: Parse
        count = self.embed_and_store(parsed_records) # Step 3: Embed & Store
        return count

    @abstractmethod
    def fetch(self, **kwargs) -> List[Dict]: ...

    @abstractmethod
    def parse(self, raw_data: List[Dict]) -> List[Dict]: ...

    def embed_and_store(self, records: List[Dict]) -> int:
        # Batch embed text fields using self.embedder
        # Insert into Milvus via self.collection_manager
        # Returns total records inserted
```

## 10.2 Parser Inventory

| Parser                    | Collection        | Data Source                    |
|---------------------------|-------------------|--------------------------------|
| `civic_parser.py`         | onco_variants     | CIViC GraphQL API              |
| `oncokb_parser.py`        | onco_variants     | OncoKB annotation files        |
| `literature_parser.py`    | onco_literature   | PubMed E-utilities API         |
| `clinical_trials_parser.py`| onco_trials      | ClinicalTrials.gov v2 API      |
| `guideline_parser.py`     | onco_guidelines   | Curated guideline seed files   |
| `resistance_parser.py`    | onco_resistance   | Curated resistance seed files  |
| `pathway_parser.py`       | onco_pathways     | Curated pathway seed files     |
| `outcome_parser.py`       | onco_outcomes     | Curated outcome records        |

## 10.3 Running Ingest

Each parser can be invoked independently:

```python
from src.ingest.literature_parser import LiteratureParser

parser = LiteratureParser(
    collection_manager=collection_manager,
    embedder=embedder,
    collection_name="onco_literature",
    batch_size=50,
)
count = parser.run(query="EGFR NSCLC targeted therapy", max_results=500)
```

The `scripts/` directory contains 17 scripts for automated ingestion,
seeding, validation, and benchmarking of all collections.

## 10.4 Embedding Strategy

All parsers use the same BGE-small-en-v1.5 model (384 dimensions). Each
domain model has a `to_embedding_text()` method that concatenates the most
semantically relevant fields:

```python
class OncologyLiterature(BaseModel):
    def to_embedding_text(self) -> str:
        parts = [self.title, self.text_chunk]
        if self.gene:    parts.append(f"Gene: {self.gene}")
        if self.variant: parts.append(f"Variant: {self.variant}")
        if self.cancer_type:
            parts.append(f"Cancer: {self.cancer_type.value}")
        if self.keywords:
            parts.append(f"Keywords: {', '.join(self.keywords)}")
        return " | ".join(parts)
```

---

# Chapter 11: Cross-Modal Integration

File: `src/cross_modal.py` (383 lines)

## 11.1 When Cross-Modal Triggers Fire

The `OncoCrossModalTrigger` fires when a case contains variants with
evidence level A or B in ACTIONABLE_TARGETS. It enriches the clinical
context by querying across modalities.

## 11.2 CrossModalResult Dataclass

```python
@dataclass
class CrossModalResult:
    trigger_reason: str               # why the trigger fired
    actionable_variants: List[Dict]   # A/B-level variants
    genomic_context: List[Dict]       # genomic evidence hits
    imaging_context: List[Dict]       # imaging findings (if available)
    genomic_hit_count: int
    imaging_hit_count: int
    enrichment_summary: str           # human-readable summary
```

## 11.3 Evaluation Flow

```
evaluate(case_or_variants)
  ├─ Extract variants from input
  ├─ Filter to actionability A or B
  │   (re-classify via ACTIONABLE_TARGETS if not pre-computed)
  ├─ If no actionable variants: return None (trigger not fired)
  │
  ├─ Build queries per actionable variant:
  │   Genomic: "{gene} {variant} targeted therapy evidence"
  │            "{gene} mutation clinical significance"
  │   Imaging: "{gene} mutation {cancer_type} imaging findings"
  │
  ├─ _query_genomics(queries) -> genomic hits
  │   Collection: "genomic_evidence"
  │   Threshold: DEFAULT_THRESHOLD = 0.40
  │   Top-K: genomic_top_k = 5
  │
  ├─ _query_imaging(queries) -> imaging hits
  │   Collection prefix: "imaging_"
  │   Graceful failure if imaging collections don't exist
  │
  └─ Build enrichment_summary and return CrossModalResult
```

## 11.4 Graceful Degradation

The imaging query is wrapped in try/except -- if the Imaging Intelligence
Agent is not deployed (no `imaging_*` collections in Milvus), the trigger
still returns genomic context without imaging data. This allows the
oncology agent to function standalone or in a multi-agent deployment.

---

# Chapter 12: Testing Strategy

Directory: `tests/` (4,370 lines, 556 tests)

## 12.1 Test File Map

| File                     | Lines | Tests | Focus                          |
|--------------------------|------:|------:|--------------------------------|
| test_models.py           |   644 |  ~120 | Enum membership, model fields  |
| test_integration.py      |   785 |   ~80 | 4 patient profiles, end-to-end |
| test_knowledge.py        |   603 |   ~70 | Knowledge graph integrity      |
| test_export.py           |   439 |   ~60 | All 4 export formats           |
| test_therapy_ranker.py   |   363 |   ~50 | Ranking, resistance, combos    |
| test_trial_matcher.py    |   363 |   ~50 | Matching, scoring, explanations|
| test_case_manager.py     |   332 |   ~40 | VCF parsing, case lifecycle    |
| test_rag_engine.py       |   301 |   ~35 | Retrieval pipeline             |
| test_collections.py      |   276 |   ~30 | Schema validation              |
| test_agent.py            |   264 |   ~21 | Planning, evaluation           |

## 12.2 Test Patterns

### Parametrized Enum Tests (test_models.py)

```python
@pytest.mark.parametrize("cancer_type", [
    CancerType.NSCLC, CancerType.BREAST, CancerType.MELANOMA, ...
])
def test_cancer_type_values(cancer_type):
    assert cancer_type.value == cancer_type.name.lower() or ...
```

### Integration Test Patient Profiles (test_integration.py)

Four synthetic patient profiles exercise the full pipeline:

1. **NSCLC with EGFR L858R** -- targeted therapy path
2. **Melanoma with BRAF V600E** -- combination therapy path
3. **CRC with MSI-H** -- biomarker-driven immunotherapy path
4. **Breast with BRCA2 mutation** -- PARP inhibitor path

Each profile tests: case creation -> VCF parsing -> actionability classification
-> therapy ranking -> trial matching -> export generation.

### Knowledge Graph Integrity (test_knowledge.py)

Validates structural consistency across all knowledge dictionaries:

- Every ACTIONABLE_TARGETS entry has required keys
- Every therapy in targeted_therapies appears in THERAPY_MAP
- Every resistance_mutations entry has corresponding RESISTANCE_MAP entries
- Pathway references are valid

### Mock Fixtures (conftest.py)

The test suite uses pytest fixtures for:
- `mock_collection_manager` -- in-memory Milvus simulation
- `mock_embedder` -- returns fixed-dimension zero vectors
- `mock_llm_client` -- returns canned responses
- `sample_vcf_text` -- synthetic VCF content with SnpEff annotations

### Export Format Tests (test_export.py)

Tests validate all four export formats against the same input data:

```python
def test_export_markdown_contains_sections(sample_mtb_packet):
    result = export_markdown(sample_mtb_packet)
    assert "# " in result                     # has heading
    assert "Somatic Variant Profile" in result # has variant section
    assert "Therapy Ranking" in result         # has therapy section
    assert "research use only" in result       # has disclaimer

def test_export_fhir_r4_valid_bundle(sample_mtb_packet):
    bundle = export_fhir_r4(sample_mtb_packet, patient_id="P001")
    assert bundle["resourceType"] == "Bundle"
    assert bundle["type"] == "collection"
    resources = [e["resource"]["resourceType"] for e in bundle["entry"]]
    assert "Patient" in resources
    assert "Observation" in resources
    assert "DiagnosticReport" in resources
```

### Therapy Ranker Tests (test_therapy_ranker.py)

Tests cover the full ranking pipeline including edge cases:

- Variant-driven therapy identification for each major gene
- Biomarker-driven therapy identification (MSI-H, TMB-H, HRD, PD-L1)
- Resistance flagging from prior therapy history
- Contraindication detection for same drug class
- Combination regimen identification
- Final rank ordering (clean before flagged)

### Trial Matcher Tests (test_trial_matcher.py)

Tests validate:

- Cancer type alias resolution
- Deterministic filter construction
- Biomarker matching scoring
- Phase and status weighting
- Age penalty computation
- Composite score calculation
- Match explanation generation

## 12.3 Running Tests

```bash
# All tests
pytest tests/ -v

# Specific module
pytest tests/test_therapy_ranker.py -v

# Integration tests only
pytest tests/test_integration.py -v -k "integration"

# With coverage
pytest tests/ --cov=src --cov-report=html
```

---

# Chapter 13: Performance Tuning

## 13.1 Collection Weight Tuning

Weights control how much each collection influences final ranking. To tune:

1. Run a set of benchmark queries with known-good answers
2. Adjust weights in `config/settings.py` (or via `ONCO_WEIGHT_*` env vars)
3. Measure retrieval quality (NDCG, MRR, or manual relevance assessment)

**Default rationale:**
- Variants (0.18) and literature (0.16) weighted highest because they carry
  the most direct clinical evidence
- Guidelines (0.12) and therapies (0.14) provide actionability context
- Cases (0.02) and genomic evidence (0.03) are supplementary

## 13.2 Milvus Index Tuning

Current defaults:

| Parameter   | Value    | Effect                                |
|-------------|----------|---------------------------------------|
| `nlist`     | 1024     | Number of IVF partitions              |
| `nprobe`    | 16       | Partitions searched per query         |
| `metric`    | COSINE   | Similarity measure                    |

Tuning guidelines:

| Scenario                    | Recommendation                      |
|-----------------------------|-------------------------------------|
| < 100K vectors/collection   | nlist=256, nprobe=16                |
| 100K-1M vectors/collection  | nlist=1024, nprobe=16-32            |
| > 1M vectors/collection     | nlist=2048, nprobe=32-64            |
| Memory constrained          | Switch to IVF_PQ (lossy)            |
| Maximum recall needed       | Use FLAT index (brute-force)        |

## 13.3 Parallel Search Workers

The RAG engine caps ThreadPoolExecutor at `min(len(collections), 8)`.
On systems with more cores, increase this ceiling. On memory-constrained
systems, reduce to 4 to limit concurrent Milvus connections.

## 13.4 Evidence Cap

`_MAX_EVIDENCE = 30` limits the number of evidence items sent to the LLM.
This balances context window utilization against prompt cost:

- Increase to 50 for models with large context windows (Claude Opus)
- Decrease to 15-20 for faster response times

## 13.5 Embedding Batch Size

`EMBEDDING_BATCH_SIZE = 32` (from settings). During ingest, texts are
batched before calling `embedder.encode()`. Increase to 64-128 on GPU
systems for throughput. During query time, only single texts are embedded.

## 13.6 Query Expansion Impact

Expanded searches use `top_k // 2` to avoid overwhelming the result set.
If expansion adds too much noise, consider:

- Reducing expansion terms per category
- Increasing the score threshold for expanded hits
- Disabling expansion for targeted strategies

## 13.7 Trial Matcher Optimization

The deterministic search issues `len(aliases) * len(statuses)` Milvus
queries per cancer type (typically 5 aliases x 4 statuses = 20 queries).
For latency-sensitive deployments:

- Pre-filter the alias list to the top 2-3 most relevant
- Cache deterministic results with a short TTL
- Increase `top_k` multiplier for semantic search, reduce for deterministic

---

# Chapter 14: Extending the Agent

## 14.1 Adding a New Collection

1. **Define schema** in `src/collections.py`:

```python
MY_NEW_FIELDS = [
    FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=100, is_primary=True),
    FieldSchema(name="embedding", dtype=DataType.FLOAT_VECTOR, dim=EMBEDDING_DIM),
    # ... your metadata fields
    FieldSchema(name="text_summary", dtype=DataType.VARCHAR, max_length=3000),
]

MY_NEW_SCHEMA = CollectionSchema(
    fields=MY_NEW_FIELDS,
    description="Your collection description",
)
```

2. **Register** in `COLLECTION_SCHEMAS` and `COLLECTION_MODELS`

3. **Add weight** in `config/settings.py`:

```python
WEIGHT_MY_NEW: float = 0.05  # adjust other weights to maintain sum ~1.0
```

4. **Add config entry** in `src/rag_engine.py` `COLLECTION_CONFIG`

5. **Create ingest parser** by subclassing `BaseIngestPipeline`

6. **Add Pydantic model** in `src/models.py` if needed

7. **Write tests** in a new `test_my_new.py` file

## 14.2 Adding a New Biomarker Rule

In `src/therapy_ranker.py`, add to `_identify_biomarker_therapies()`:

```python
# Example: TROP2 overexpression -> sacituzumab govitecan
trop2 = biomarkers.get("TROP2", "").upper()
if trop2 in ("POSITIVE", "OVEREXPRESSION", "HIGH"):
    therapies.append({
        "drug_name": "sacituzumab govitecan",
        "brand_name": "Trodelvy",
        "category": "ADC",
        "targets": ["TROP2"],
        "evidence_level": "A",
        "guideline_recommendation": "FDA-approved for TROP2+ TNBC and urothelial.",
        "source": "biomarker",
        "source_biomarker": "TROP2",
    })
```

Then add corresponding entries to:
- `BIOMARKER_PANELS` in `src/knowledge.py`
- Query expansion terms in `src/query_expansion.py`
- Test cases in `tests/test_therapy_ranker.py`

## 14.3 Adding a New Drug Class Group

In `src/therapy_ranker.py`, extend `_DRUG_CLASS_GROUPS`:

```python
_DRUG_CLASS_GROUPS = {
    ...
    "her2_adc": ["trastuzumab deruxtecan", "trastuzumab emtansine"],
    "trop2_adc": ["sacituzumab govitecan", "datopotamab deruxtecan"],
}
```

## 14.4 Adding a New Ingest Parser

1. Create `src/ingest/my_source_parser.py`:

```python
from src.ingest.base import BaseIngestPipeline

class MySourceParser(BaseIngestPipeline):
    def fetch(self, query=None, max_results=None):
        # Call external API or read seed files
        return [{"raw_field": "value", ...}, ...]

    def parse(self, raw_data):
        # Normalize to collection schema fields
        records = []
        for item in raw_data:
            records.append({
                "id": ...,
                "text_summary": ...,  # this gets embedded
                # ... other schema fields
            })
        return records
```

2. Register in `src/ingest/__init__.py`
3. Add a script in `scripts/` to run the parser

## 14.5 Adding a New Export Format

Follow the pattern in `src/export.py`:

```python
def export_my_format(mtb_packet_or_response: Any, **kwargs) -> Any:
    data = _normalise_input(mtb_packet_or_response)
    # Transform data into your target format
    return result
```

The `_normalise_input()` helper handles MTBPacket, dict, and str inputs.

## 14.6 Adding a New API Route

1. Create `api/routes/my_router.py`:

```python
from fastapi import APIRouter, Depends
router = APIRouter(prefix="/my-endpoint", tags=["my-feature"])

@router.post("/action")
async def my_action(request: MyRequest):
    state = get_state()
    # Use state["rag_engine"], state["collection_manager"], etc.
    return {"result": ...}
```

2. Import and include in `api/main.py`:

```python
from api.routes import my_router
app.include_router(my_router.router)
```

## 14.7 Adding a Cross-Modal Trigger

Extend `OncoCrossModalTrigger.evaluate()` in `src/cross_modal.py` to query
additional collection prefixes. The pattern is:

1. Define the new collection prefix (e.g., `DRUG_DISCOVERY_PREFIX = "drug_"`)
2. Build domain-specific queries for actionable variants
3. Query the collections with graceful failure handling
4. Add results to the `CrossModalResult`

## 14.8 Modifying the System Prompt

The system prompt in `src/rag_engine.py` (`ONCO_SYSTEM_PROMPT`) defines
the agent's persona and behavioral constraints. To modify:

1. Edit the `ONCO_SYSTEM_PROMPT` string in `src/rag_engine.py`
2. Keep the 8 competency areas unless deliberately removing capability
3. Preserve the 5 behavioral instructions (cite evidence, cross-functional
   thinking, resistance/contraindications, guideline references, uncertainty)
4. Test with `test_rag_engine.py` to verify prompt construction still works

## 14.9 Changing the Embedding Model

To switch from BGE-small-en-v1.5 to a different model:

1. Update `EMBEDDING_MODEL` and `EMBEDDING_DIM` in `config/settings.py`
2. Update `EMBEDDING_DIM` constant in `src/collections.py`
3. All existing collections must be dropped and recreated (dimension mismatch
   will cause Milvus errors)
4. Re-run all ingest pipelines to re-embed existing data
5. Update the BGE instruction prefix `_BGE_INSTRUCTION` in `src/rag_engine.py`
   (different models may use different instruction prefixes or none at all)

## 14.10 Working with the Metrics Module

The `src/metrics.py` module (362 lines) provides Prometheus instrumentation.
Key metrics tracked:

- `onco_query_total` -- total queries processed (counter)
- `onco_query_duration_seconds` -- query latency (histogram)
- `onco_evidence_hits` -- evidence items retrieved per query (histogram)
- `onco_collection_search_duration` -- per-collection search time (histogram)
- `onco_therapy_candidates` -- therapies identified per ranking (histogram)
- `onco_trial_matches` -- trials matched per patient (histogram)

Metrics are enabled by default (`METRICS_ENABLED=True`). Disable in
development with `ONCO_METRICS_ENABLED=false`.

---

# Appendix A: Complete API Reference

## FastAPI Application (`api/main.py`)

| Endpoint            | Method | Router       | Description                         |
|---------------------|--------|--------------|-------------------------------------|
| `/health`           | GET    | main         | Liveness check with component status|
| `/api/v1/query`     | POST   | meta_agent   | Full agent query (plan/search/synth)|
| `/api/v1/search`    | POST   | meta_agent   | Evidence-only search (no LLM)       |
| `/api/v1/compare`   | POST   | meta_agent   | Comparative retrieval               |
| `/api/v1/cases`     | POST   | cases        | Create case from VCF/variants       |
| `/api/v1/cases/{id}`| GET    | cases        | Retrieve case by ID                 |
| `/api/v1/cases/{id}/mtb` | GET | cases      | Generate MTB packet                 |
| `/api/v1/trials`    | POST   | trials       | Match trials for patient profile    |
| `/api/v1/reports`   | POST   | reports      | Export report (markdown/json/pdf/fhir)|
| `/api/v1/events`    | POST   | events       | Cross-modal trigger evaluation      |

## Startup Configuration

The `lifespan()` function initializes all components in order:

```
1. Load OncoSettings (ONCO_ env prefix, .env file)
2. Connect OncoCollectionManager to Milvus
3. Load SentenceTransformer (BGE-small-en-v1.5) via EmbedderWrapper
4. Initialize OncoRAGEngine
5. Initialize OncoIntelligenceAgent
6. Initialize OncologyCaseManager
7. Initialize TrialMatcher
8. Initialize TherapyRanker
9. Initialize OncoCrossModalTrigger
10. Store all in _state dict for dependency injection
```

## EmbedderWrapper

The `EmbedderWrapper` class in `api/main.py` adapts SentenceTransformer
to expose both `.encode()` and `.embed()` APIs:

```python
class EmbedderWrapper:
    def __init__(self, model: SentenceTransformer):
        self._model = model

    def encode(self, texts):
        """SentenceTransformer native API -- accepts str or list."""
        return self._model.encode(texts)

    def embed(self, text) -> list:
        """Single-text convenience. Returns list of floats (384-dim)."""
        if isinstance(text, str):
            return self._model.encode([text])[0].tolist()
        return self._model.encode(text).tolist()
```

The RAG engine uses `.encode()` (inherited from SentenceTransformer).
The therapy ranker, trial matcher, and case manager use `.embed()` for
single-text embedding. Both methods produce identical 384-dimensional
vectors.

## Request/Response Models

The API uses Pydantic models for request validation:

```python
# Query endpoint
class QueryRequest(BaseModel):
    question: str
    top_k: int = 5
    collections_filter: Optional[List[str]] = None
    year_min: Optional[int] = None
    year_max: Optional[int] = None

# Case creation endpoint
class CreateCaseRequest(BaseModel):
    patient_id: str
    cancer_type: str
    stage: str
    vcf_content: Optional[str] = None
    variants: Optional[List[Dict]] = None
    biomarkers: Optional[Dict[str, Any]] = None
    prior_therapies: Optional[List[str]] = None

# Trial matching endpoint
class TrialMatchRequest(BaseModel):
    cancer_type: str
    biomarkers: Dict[str, Any]
    stage: str
    age: Optional[int] = None
    top_k: int = 10
```

## CORS Configuration

The FastAPI app enables CORS for development:

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

For production deployments, restrict `allow_origins` to the Streamlit UI
host and any frontend applications.

---

# Appendix B: Collection Schema Reference

## Summary Table

| Collection        | PK Length | Text Field    | Text Max | Typed Metadata Fields           |
|-------------------|----------|---------------|----------|---------------------------------|
| onco_variants     | 100      | text_summary  | 3000     | gene, variant_name, evidence_level, drugs, civic_id |
| onco_literature   | 100      | text_chunk    | 3000     | title, source_type, year, gene, journal |
| onco_trials       | 20       | text_summary  | 3000     | phase, status, sponsor, cancer_types, enrollment |
| onco_biomarkers   | 100      | text_summary  | 3000     | name, biomarker_type, predictive_value, testing_method |
| onco_therapies    | 100      | text_summary  | 3000     | drug_name, category, targets, mechanism_of_action |
| onco_pathways     | 100      | text_summary  | 3000     | name, key_genes, therapeutic_targets, cross_talk |
| onco_guidelines   | 100      | text_summary  | 3000     | org, cancer_type, version, year, evidence_level |
| onco_resistance   | 100      | text_summary  | 3000     | primary_therapy, gene, mechanism, bypass_pathway |
| onco_outcomes     | 100      | text_summary  | 3000     | case_id, therapy, response, duration_months |
| onco_cases        | 100      | text_summary  | 3000     | patient_id, cancer_type, stage, variants, biomarkers |
| genomic_evidence  | 200      | text_summary  | 2000     | chrom, pos, gene, consequence, impact, am_pathogenicity |

## Common Index Configuration

All collections share:
- Embedding: `FLOAT_VECTOR`, dim=384 (BGE-small-en-v1.5)
- Index: `IVF_FLAT`, nlist=1024
- Metric: `COSINE`
- Search: nprobe=16

---

# Appendix C: Configuration Parameter Reference

All parameters use the `ONCO_` environment variable prefix (via Pydantic
BaseSettings). Override any value by setting `ONCO_<PARAM_NAME>` in your
environment or `.env` file.

## Paths

| Parameter        | Default                         | Description                    |
|------------------|---------------------------------|--------------------------------|
| PROJECT_ROOT     | `<auto-detected>`               | Repository root                |
| DATA_DIR         | `{PROJECT_ROOT}/data`           | Data file directory            |
| CACHE_DIR        | `{PROJECT_ROOT}/cache`          | Cache directory                |
| REFERENCE_DIR    | `{PROJECT_ROOT}/reference`      | Reference file directory       |
| RAG_PIPELINE_ROOT| `<parent>/rag-chat-pipeline`    | Shared RAG pipeline root       |

## Milvus

| Parameter   | Default     | Description               |
|-------------|-------------|---------------------------|
| MILVUS_HOST | `localhost` | Milvus server hostname    |
| MILVUS_PORT | `19530`     | Milvus server port        |

## Embeddings

| Parameter          | Default               | Description              |
|--------------------|-----------------------|--------------------------|
| EMBEDDING_MODEL    | `BAAI/bge-small-en-v1.5` | HuggingFace model ID  |
| EMBEDDING_DIM      | `384`                 | Vector dimensionality    |
| EMBEDDING_BATCH_SIZE| `32`                 | Batch size for ingest    |

## LLM

| Parameter       | Default                 | Description             |
|-----------------|-------------------------|-------------------------|
| LLM_PROVIDER    | `anthropic`             | LLM provider            |
| LLM_MODEL       | `claude-sonnet-4-20250514`  | Model identifier        |
| ANTHROPIC_API_KEY| `None`                 | Anthropic API key       |

## RAG Search

| Parameter       | Default | Description                              |
|-----------------|---------|------------------------------------------|
| TOP_K           | `5`     | Per-collection hit limit                 |
| SCORE_THRESHOLD | `0.4`   | Minimum similarity score                 |

## Collection Weights

| Parameter          | Default | Collection          |
|--------------------|--------:|---------------------|
| WEIGHT_VARIANTS    |    0.18 | onco_variants       |
| WEIGHT_LITERATURE  |    0.16 | onco_literature     |
| WEIGHT_THERAPIES   |    0.14 | onco_therapies      |
| WEIGHT_GUIDELINES  |    0.12 | onco_guidelines     |
| WEIGHT_TRIALS      |    0.10 | onco_trials         |
| WEIGHT_BIOMARKERS  |    0.08 | onco_biomarkers     |
| WEIGHT_RESISTANCE  |    0.07 | onco_resistance     |
| WEIGHT_PATHWAYS    |    0.06 | onco_pathways       |
| WEIGHT_OUTCOMES    |    0.04 | onco_outcomes       |
| WEIGHT_CASES       |    0.02 | onco_cases          |
| WEIGHT_GENOMIC     |    0.03 | genomic_evidence    |

## External APIs

| Parameter          | Default                                | Description              |
|--------------------|----------------------------------------|--------------------------|
| NCBI_API_KEY       | `None`                                 | NCBI E-utilities API key |
| PUBMED_MAX_RESULTS | `5000`                                 | Max PubMed fetch results |
| CT_GOV_BASE_URL    | `https://clinicaltrials.gov/api/v2`    | ClinicalTrials.gov API   |
| CIVIC_BASE_URL     | `https://civicdb.org/api`              | CIViC API endpoint       |

## Server

| Parameter      | Default     | Description                |
|----------------|-------------|----------------------------|
| API_HOST       | `0.0.0.0`  | FastAPI bind address       |
| API_PORT       | `8527`      | FastAPI port               |
| STREAMLIT_PORT | `8526`      | Streamlit UI port          |

## Operational

| Parameter                | Default | Description                      |
|--------------------------|---------|----------------------------------|
| METRICS_ENABLED          | `True`  | Enable Prometheus metrics        |
| SCHEDULER_INTERVAL       | `168h`  | Background task interval (7 days)|
| CONVERSATION_MEMORY_DEPTH| `3`     | Conversation turns to retain     |

## Scheduler Configuration

| Parameter           | Default  | Description                         |
|---------------------|----------|-------------------------------------|
| SCHEDULER_INTERVAL  | `168h`   | Ingest re-run interval (7 days)     |

The scheduler (`src/scheduler.py`) runs background ingest tasks at the
configured interval. Each run refreshes collection data from external
sources (PubMed, ClinicalTrials.gov, CIViC) without requiring manual
intervention.

## Evidence Level Labels (for Export)

The export module maps internal evidence levels to human-readable labels:

```python
EVIDENCE_LEVEL_LABELS = {
    "level_1": "Level 1 -- FDA-approved / Standard of Care",
    "level_2": "Level 2 -- Clinical Evidence / Consensus",
    "level_3": "Level 3 -- Case Reports / Early Trials",
    "level_4": "Level 4 -- Preclinical / Biological Rationale",
    "level_R": "Level R -- Resistance Evidence",
}
```

## Environment Variable Quick Reference

Set any parameter by prefixing with `ONCO_`:

```bash
# Example .env file
ONCO_MILVUS_HOST=milvus-server.internal
ONCO_MILVUS_PORT=19530
ONCO_EMBEDDING_MODEL=BAAI/bge-small-en-v1.5
ONCO_LLM_MODEL=claude-sonnet-4-20250514
ONCO_ANTHROPIC_API_KEY=sk-ant-...
ONCO_WEIGHT_VARIANTS=0.20
ONCO_WEIGHT_LITERATURE=0.18
ONCO_TOP_K=10
ONCO_SCORE_THRESHOLD=0.35
ONCO_API_PORT=8527
ONCO_METRICS_ENABLED=true
ONCO_CONVERSATION_MEMORY_DEPTH=5
```

---

*This guide reflects the codebase as of February 2026. For introductory
material, see `LEARNING_GUIDE.md`. For deployment instructions, see
`DEPLOYMENT.md`.*
