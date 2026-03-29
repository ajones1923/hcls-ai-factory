# Precision Oncology Intelligence Agent -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [System Diagram](#1-system-diagram)
2. [Component Interactions](#2-component-interactions)
3. [Data Flow](#3-data-flow)
4. [Collection Design Rationale](#4-collection-design-rationale)
5. [Therapy Ranker Engine](#5-therapy-ranker-engine)
6. [Trial Matcher Engine](#6-trial-matcher-engine)
7. [Resistance Analysis Engine](#7-resistance-analysis-engine)
8. [Cross-Modal Engine](#8-cross-modal-engine)
9. [Query Expansion](#9-query-expansion)
10. [RAG Pipeline](#10-rag-pipeline)
11. [Agent Orchestrator](#11-agent-orchestrator)
12. [Data Model Architecture](#12-data-model-architecture)

---

## 1. System Diagram

### 1.1 Full System Architecture

```
                          EXTERNAL USERS
                               |
                    +----------+----------+
                    |                     |
              +-----+------+      +------+-----+
              | Streamlit  |      | REST API   |
              | UI :8526   |      | :8527      |
              +-----+------+      +------+-----+
                    |                     |
                    +----------+----------+
                               |
                    +----------+----------+
                    |  Agent Orchestrator  |
                    |  (src/agent.py)      |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Query       |    | Workflow      |    | Clinical      |
   | Expansion   |    | Engine        |    | Engines       |
   | (812 LOC)   |    |               |    |               |
   +------+------+    +-------+-------+    | Therapy       |
          |                    |           | Ranker        |
          |                    |           | (748 LOC)     |
          |                    |           |               |
          |                    |           | Trial         |
          |                    |           | Matcher       |
          |                    |           | (513 LOC)     |
          |                    |           |               |
          |                    |           | Case          |
          |                    |           | Manager       |
          |                    |           | (516 LOC)     |
          |                    |           +-------+-------+
          |                    |                    |
          +--------------------+--------------------+
                               |
                    +----------+----------+
                    |    RAG Engine       |
                    |    (908 LOC)        |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Knowledge   |    | Milvus        |    | LLM           |
   | Graph       |    | Vector DB     |    | (Claude 4.6)  |
   | (1,662 LOC) |    | 11 Collections|    |               |
   +-------------+    +-------+-------+    +---------------+
                               |
                    +----------+----------+
                    |  etcd    |  MinIO   |
                    +----------+----------+
```

### 1.2 Ingest Pipeline Architecture

```
  +------------------------------------------+
  |           External Data Sources           |
  +------------------------------------------+
  | PubMed | ClinicalTrials | NCCN/ESMO | FDA|
  +---+--------+--------+--------+--------+--+
      |        |        |        |        |
  +---v---+ +--v---+ +-v----+ +-v-----+ +v--------+
  |PubMed | |Trial | |Guide | |Variant| |Therapy  |
  |Parser | |Parser| |Parser| |Parser | |Parser   |
  +---+---+ +--+---+ +-+----+ +-+-----+ ++--------+
      |        |        |        |        |
      +--------+--------+--------+--------+
                         |
                  +------v------+
                  | Base Parser |
                  | - Chunking  |
                  | - Embedding |
                  | - Insertion |
                  +------+------+
                         |
                  +------v------+
                  |   Milvus    |
                  | Collections |
                  +-------------+
```

---

## 2. Component Interactions

### 2.1 Component Dependency Graph

```
OncoUI (Streamlit) ──> FastAPI Server ──> OncoIntelligenceAgent
                                               |
                                 +-------------+-------------+
                                 |             |             |
                          QueryExpansion  CaseManager     TherapyRanker
                                 |             |             |
                                 +------+------+      TrialMatcher
                                        |                    |
                                   RAGEngine          CrossModalEngine
                                        |                    |
                                 +------+------+             |
                                 |             |             |
                            Milvus DB    Knowledge      genomic_evidence
                            (11 cols)     Graph          (shared col)
```

### 2.2 Module Responsibilities

| Module | File | LOC | Responsibilities |
|--------|------|-----|-----------------|
| **Agent Orchestrator** | `src/agent.py` | 553 | Plan-search-evaluate-synthesize loop; top-level coordination |
| **Query Expansion** | `src/query_expansion.py` | 812 | 12-category domain-aware expansion, entity extraction, synonym resolution |
| **RAG Engine** | `src/rag_engine.py` | 908 | Multi-collection search, citation scoring, LLM prompt assembly and synthesis |
| **Therapy Ranker** | `src/therapy_ranker.py` | 748 | 7-step evidence-based therapy ranking with 80+ therapies mapped |
| **Trial Matcher** | `src/trial_matcher.py` | 513 | Hybrid deterministic + semantic trial matching |
| **Case Manager** | `src/case_manager.py` | 516 | VCF parsing, case CRUD, MTB packet generation |
| **Knowledge Graph** | `src/knowledge.py` | 1,662 | 40+ actionable targets, 80+ therapies, 12+ resistance mechanisms |
| **Collections** | `src/collections.py` | -- | Milvus collection CRUD, schema definitions, index management |
| **Export** | `src/export.py` | 1,055 | PDF, CSV, JSON, FHIR R4 report generation |
| **Cross-Modal** | `src/cross_modal.py` | -- | Genomic-imaging-drug discovery event integration |
| **Metrics** | `src/metrics.py` | -- | Prometheus instrumentation |
| **Models** | `src/models.py` | 538 | 13 enums, 10 Pydantic domain models |

### 2.3 Interface Contracts

**Agent Orchestrator inputs/outputs:**
```
Input:  OncoQuery(question, patient_context?, case_id?)
Output: AgentResponse(answer, evidence, knowledge_used, report, confidence)
```

**Therapy Ranker inputs/outputs:**
```
Input:  PatientProfile(cancer_type, variants, biomarkers, prior_therapies)
Output: TherapyRanking(ranked_therapies, flagged_therapies, combination_options)
```

**Trial Matcher inputs/outputs:**
```
Input:  PatientProfile(cancer_type, biomarkers, stage, age)
Output: TrialMatchResult(matches, composite_scores, match_explanations)
```

---

## 3. Data Flow

### 3.1 Query Processing Pipeline

```
Step 1: RECEIVE QUERY
  OncoQuery arrives via API or UI
  |
Step 2: PLAN (agent.py)
  Identify topics from keyword matching (20+ triggers)
  Extract target genes (40+ ACTIONABLE_TARGETS)
  Detect cancer types (26 canonical + 70+ aliases)
  Select strategy: broad | targeted | comparative
  Decompose complex queries into sub-questions
  Produce SearchPlan
  |
Step 3: SEARCH (rag_engine.py)
  Embed expanded query via BGE-small-en-v1.5
  Search 11 collections in parallel (top_k=5 each)
  Optional query expansion (12 categories)
  Weighted scoring with per-collection weights
  Filter by score_threshold (0.30)
  Deduplicate results
  Merge and rank: top 30
  |
Step 4: EVALUATE (agent.py)
  sufficient: >= 3 hits from >= 2 collections
  partial: some hits but below threshold
  insufficient: zero quality hits -> broaden and retry
  MAX_RETRIES = 2
  |
Step 5: KNOWLEDGE INJECTION (knowledge.py)
  Inject gene context (40+ ACTIONABLE_TARGETS)
  Inject therapy context (80+ therapies, NCCN/ESMO)
  Inject resistance mechanisms (12+ documented mechanisms)
  Inject pathway context (10+ oncogenic pathways)
  Inject biomarker panels (20+ clinical panels)
  |
Step 6: LLM SYNTHESIS (rag_engine.py)
  Assemble context: domain knowledge + evidence + question
  Send to Claude Sonnet 4.6 with oncology system prompt
  Generate evidence-grounded answer with inline citations
  |
Step 7: RESPONSE ASSEMBLY (agent.py)
  Package AgentResponse with:
  - Synthesized answer
  - Ranked evidence with citations
  - Knowledge domains used
  - Markdown report
  - Confidence score
  |
Step 8: DELIVERY
  Return via API (JSON) or display in UI (formatted)
  Optional: Export as PDF/JSON/FHIR R4/Markdown
```

### 3.2 Case Creation Flow

```
Patient Data (ID, cancer type, stage, VCF, biomarkers, prior therapies)
    |
    v
[VCF Parsing] -- If raw VCF text, parse via cyvcf2
    |   Extract: gene, variant, chrom, pos, ref, alt, consequence
    |
    v
[Variant Annotation] -- Cross-reference ACTIONABLE_TARGETS (40+ genes)
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

## 4. Collection Design Rationale

### 4.1 Why 11 Collections (Not One Giant Collection)

The multi-collection architecture was chosen over a single monolithic collection for five reasons:

**1. Semantic Precision**

The same terms carry different meanings across oncology subspecialties:
- "Resistance": acquired mutations (T790M), pathway bypass (MET amplification), lineage plasticity
- "Response": RECIST criteria, molecular response (ctDNA), pathological complete response
- "Marker": predictive biomarker (PD-L1), prognostic biomarker (TMB), diagnostic biomarker (FISH)

Collection-specific embedding spaces preserve these semantic boundaries.

**2. Relevance Weighting**

Different queries require different emphasis:
- "EGFR T790M resistance" -> weight onco_variants (0.18), onco_resistance (0.07)
- "First-line NSCLC therapy" -> weight onco_therapies (0.14), onco_guidelines (0.12)
- "MSI-H trial eligibility" -> weight onco_trials (0.10), onco_biomarkers (0.08)

A single collection would return mixed results without the ability to boost domain-specific content.

**3. Independent Lifecycle Management**

| Collection | Update Frequency | Source |
|-----------|-----------------|--------|
| onco_literature | Weekly (PubMed) | Automated |
| onco_guidelines | Quarterly (NCCN) | Manual review |
| onco_trials | Monthly | Automated |
| onco_variants | Monthly (CIViC/OncoKB) | Automated |
| genomic_evidence | Shared, maintained by genomics pipeline | External |

**4. Source Attribution**

Clinicians need to know the provenance of every citation:
- "This comes from the NCCN NSCLC v4.2026 Guideline" (onco_guidelines)
- "This is based on the KEYNOTE-024 trial" (onco_trials)
- "This variant is documented in CIViC as Level A evidence" (onco_variants)

Collection names provide automatic source categorization.

**5. Scalability**

Individual collections can be independently reindexed, compacted, partitioned, or replicated.

### 4.2 Collection Inventory

| # | Collection Name | Description | Weight | Seed Records |
|---|----------------|-------------|--------|-------------|
| 1 | onco_variants | Actionable somatic/germline variants (CIViC/OncoKB) | 0.18 | 130 |
| 2 | onco_literature | PubMed/PMC/preprint literature chunks | 0.16 | 60 |
| 3 | onco_therapies | Approved and investigational therapies | 0.14 | 94 |
| 4 | onco_guidelines | NCCN/ASCO/ESMO guideline recommendations | 0.12 | 45 |
| 5 | onco_trials | ClinicalTrials.gov summaries with biomarker criteria | 0.10 | 55 |
| 6 | onco_biomarkers | Predictive and prognostic biomarkers | 0.08 | 50 |
| 7 | onco_resistance | Resistance mechanisms and bypass strategies | 0.07 | 50 |
| 8 | onco_pathways | Signaling pathways, cross-talk, druggable nodes | 0.06 | 35 |
| 9 | onco_outcomes | Real-world treatment outcome records | 0.04 | 40 |
| 10 | onco_cases | De-identified patient case snapshots | 0.02 | 37 |
| 11 | genomic_evidence | Shared VCF-derived genomic variants (read-only) | 0.03 | -- |
| | **Total** | | **1.00** | **596** |

### 4.3 Collection Weight Distribution

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

### 4.4 Key Collection Schemas

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

## 5. Therapy Ranker Engine

### 5.1 Architecture

```
Patient Profile (cancer_type, variants, biomarkers, prior_therapies)
    |
    v
[Step 1: Variant-Driven Therapies]
    |   ACTIONABLE_TARGETS lookup for each gene/variant (40+ genes)
    |   Evidence level from knowledge graph
    |
    v
[Step 2: Biomarker-Driven Therapies]
    |   MSI-H -> pembrolizumab, nivolumab, dostarlimab (Level A)
    |   TMB-H (>=10 mut/Mb) -> pembrolizumab (Level A)
    |   HRD/BRCA -> olaparib, rucaparib, niraparib, talazoparib (Level A/B)
    |   PD-L1 TPS >=50% -> pembrolizumab first-line (Level A)
    |   NTRK fusion -> larotrectinib, entrectinib (Level A)
    |   + BIOMARKER_PANELS registry check (20+ panels)
    |
    v
[Step 3: Evidence Level Sort]
    |   A (FDA-approved) > B (Clinical) > C (Case reports) > D > E
    |
    v
[Step 4: Resistance Check]
    |   RESISTANCE_MAP: 12+ mutation-level resistance mechanisms
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
    Total therapy mappings: 80+
```

### 5.2 Evidence Level Classification

| Level | Description | Source |
|-------|------------|--------|
| A | FDA-approved companion diagnostic or indication | FDA label, NCCN Category 1 |
| B | Well-powered clinical evidence | Phase 2/3 trials, ESMO MCBS |
| C | Case reports, small series | Published case reports |
| D | Preclinical or early clinical | Phase 1, in vitro data |
| E | Computational prediction | In silico, pathway inference |

---

## 6. Trial Matcher Engine

### 6.1 Architecture

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

### 6.2 Cancer Type Alias Resolution

The trial matcher resolves 70+ cancer type aliases to 26 canonical types to improve deterministic matching. Examples:

- "NSCLC" -> "non-small cell lung cancer"
- "CRC" -> "colorectal cancer"
- "TNBC" -> "triple-negative breast cancer"
- "HCC" -> "hepatocellular carcinoma"
- "ccRCC" -> "clear cell renal cell carcinoma"

---

## 7. Resistance Analysis Engine

### 7.1 Resistance Map (12+ Mechanisms)

The knowledge graph documents resistance mechanisms for major targeted therapies:

| Drug Class | Primary Resistance | Mechanism | Next-Line |
|-----------|-------------------|-----------|-----------|
| EGFR TKI (1st/2nd gen) | EGFR T790M | Gatekeeper mutation | Osimertinib |
| EGFR TKI (3rd gen) | EGFR C797S | Binding site mutation | Combination strategies |
| BRAF inhibitor | MAPK reactivation | MEK bypass, NRAS mutation | BRAF+MEK combo |
| ALK TKI (crizotinib) | ALK G1202R | Solvent front mutation | Lorlatinib |
| HER2 therapy | HER3 upregulation | Bypass signaling | HER3 antibodies |
| PARP inhibitor | BRCA reversion | Homologous recombination restored | Platinum rechallenge |
| Anti-PD-1/PD-L1 | B2M loss | Antigen presentation defect | Combination IO |
| KRAS G12C inhibitor | KRAS amplification | On-target amplification | Combination strategies |
| RET inhibitor | RET V804M | Gatekeeper mutation | Next-gen RET inhibitors |
| MET inhibitor | MET D1228N | Kinase domain mutation | Combination or next-gen |
| PIK3CA inhibitor | PTEN loss | PI3K pathway reactivation | Combination strategies |
| CDK4/6 inhibitor | RB1 loss | Cell cycle bypass | Chemotherapy |

### 7.2 Drug Class Grouping

Same-class resistance flagging groups drugs by mechanism of action, so that a patient who progressed on erlotinib is flagged for all first/second-generation EGFR TKIs, not just erlotinib.

---

## 8. Cross-Modal Engine

### 8.1 Genomic-Imaging-Drug Discovery Integration

```
Genomics Pipeline (Stage 1)
    |-- VCF -> genomic_evidence collection
    |
    v
Oncology Agent (Stage 2.5)
    |-- Variant interpretation
    |-- Therapy ranking
    |-- Trial matching
    |
    v
Drug Discovery Pipeline (Stage 3)
    |-- Target validation
    |-- Lead compound optimization
    |-- DiffDock binding prediction
```

### 8.2 Cross-Agent Event Publishing

The agent publishes events for consumption by other HCLS AI Factory agents:

- `onco.variant.actionable` -- New actionable variant identified
- `onco.therapy.recommended` -- Therapy recommendation generated
- `onco.trial.matched` -- Patient-trial match found
- `onco.resistance.detected` -- Resistance mechanism identified
- `onco.case.created` -- New case snapshot persisted

---

## 9. Query Expansion

### 9.1 Expansion Pipeline

```
Raw Query: "What targeted therapy options exist for EGFR exon 19 deletion in NSCLC?"
    |
    v
Step 1: Entity Extraction
  - Genes: [EGFR]
  - Variants: [exon 19 deletion]
  - Cancer types: [NSCLC]
  - Topics: [targeted therapy]
    |
    v
Step 2: Synonym Expansion (12 categories)
  - "NSCLC" -> "non-small cell lung cancer, lung adenocarcinoma"
  - "EGFR exon 19 deletion" -> "del19, E746_A750del"
  - "targeted therapy" -> "TKI, tyrosine kinase inhibitor"
    |
    v
Step 3: Sub-Question Decomposition
  - Q1: "What are FDA-approved EGFR TKIs for exon 19 deletion?"
  - Q2: "What resistance mechanisms occur with EGFR-targeted therapy?"
  - Q3: "What clinical trials are open for EGFR-mutant NSCLC?"
    |
    v
Step 4: Strategy Selection
  - Strategy: "targeted" (specific gene/variant/cancer combination)
    |
    v
Output: SearchPlan
```

### 9.2 Expansion Categories (12)

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

## 10. RAG Pipeline

### 10.1 RAG Architecture

```
SearchPlan (from Query Expansion)
        |
  Embedding (BGE-small-en-v1.5, 384-dim)
        |
  Multi-Collection Search
  +------+------+------+------+------+------+
  |var   |lit   |ther  |gdl   |trial |biom  |
  |top5  |top5  |top5  |top5  |top5  |top5  |
  +------+------+------+------+------+------+
  |resist|path  |outc  |cases |gen   |
  |top5  |top5  |top5  |top5  |top5  |
  +------+------+------+------+------+
        |
  Score Filtering (threshold 0.30)
        |
  Weight Application (per-collection weights)
        |
  Deduplication (content hash)
        |
  Citation Scoring
  - High confidence: score >= 0.85
  - Medium confidence: score >= 0.65
  - Standard: score < 0.65
        |
  Knowledge Injection
  - Gene context (ACTIONABLE_TARGETS)
  - Therapy context (THERAPY_MAP)
  - Resistance context (RESISTANCE_MAP)
  - Pathway context (PATHWAY_MAP)
  - Biomarker context (BIOMARKER_PANELS)
        |
  LLM Synthesis (Claude Sonnet 4.6)
  - System: Oncology domain expert (MTB-ready)
  - User: Clinical question + retrieved context + knowledge
  - Instructions: Cite sources, follow NCCN/ESMO, flag resistance
        |
  Response Parsing
  - Extract inline citations
  - Map PubMed IDs and NCT IDs to links
  - Calculate confidence score
```

### 10.2 Embedding Configuration

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

### 10.3 Comparative Retrieval

The engine detects comparative questions via regex
(`compare|vs|versus|difference between|head.to.head`) and routes them to
a dual-entity retrieval pipeline:

1. Parse entity A and entity B from the question
2. Retrieve evidence independently for each entity
3. Identify shared/head-to-head evidence (intersection by ID)
4. Build a structured comparison prompt with 8 comparison axes
5. Generate comparative synthesis via LLM

---

## 11. Agent Orchestrator

### 11.1 Orchestration Flow

The OncoIntelligenceAgent (`src/agent.py`, 553 lines) executes a 4-step
plan-search-evaluate-synthesize loop:

```python
class OncoIntelligenceAgent:
    def __init__(self):
        self.rag_engine = OncoRAGEngine()
        self.therapy_ranker = TherapyRanker()
        self.trial_matcher = TrialMatcher()
        self.case_manager = OncologyCaseManager()
        self.query_expander = QueryExpander()
        self.cross_modal = CrossModalIntegrator()
        self.export_system = OncologyExporter()

    async def process_query(self, query: OncoQuery) -> AgentResponse:
        # 1. Plan
        search_plan = self.query_expander.expand(query)

        # 2. Search
        evidence = await self.rag_engine.cross_collection_search(search_plan)

        # 3. Evaluate
        verdict = self._evaluate_evidence(evidence)
        if verdict == "insufficient" and retries < MAX_RETRIES:
            # Broaden and retry
            ...

        # 4. Synthesize
        answer = await self.rag_engine.synthesize(query, evidence, knowledge)

        return AgentResponse(
            answer=answer.text,
            evidence=evidence,
            knowledge_used=knowledge_domains,
            report=markdown_report,
            confidence=answer.confidence
        )
```

### 11.2 Error Handling Strategy

The orchestrator implements graceful degradation:

1. **Milvus unavailable**: Returns error with search-only mode recommendation
2. **LLM unavailable**: Returns search results without synthesis
3. **Therapy ranker error**: Logs warning, returns response without rankings
4. **Trial matcher error**: Logs warning, returns response without matches
5. **Cross-modal error**: Logs warning, returns response without triggers
6. **Timeout**: Returns partial results with timeout indicator

### 11.3 Key Statistics

| Metric | Value |
|--------|-------|
| Python files | 66 |
| Total lines of code | ~20,490 |
| Milvus collections | 11 (10 owned + 1 read-only shared) |
| Actionable gene targets | 40+ |
| Therapy mappings | 80+ |
| Resistance mechanisms | 12+ |
| Oncogenic pathways | 10+ |
| Biomarker panels | 20+ |
| Test files / cases | 10 files, 556 test cases, all passing |
| Docker services | 6 |
| Export formats | 4 (Markdown, JSON, PDF, FHIR R4) |
| Cancer types supported | 26 |

---

## 12. Data Model Architecture

### 12.1 Model Hierarchy

```
                    AgentResponse (top-level output)
                    /        |         \           \
              answer    evidence    knowledge_used   report
                |           |            |              |
              str    List[Dict]    List[str]         str

                    TherapyRanking (from TherapyRanker)
                    /        |         \           \
         ranked_therapies  flagged  combinations  evidence_level
                |           |            |              |
         List[Therapy]  List[Flag]  List[Combo]     str

                    MTBPacket (from CaseManager)
                    /        |         \           \
           case_snapshot  variant_table  rankings  trial_matches
                |              |            |          |
         CaseSnapshot    List[Variant]  TherapyRanking TrialMatchResult
```

### 12.2 Enum Design

Enums enforce type safety and valid value sets:

- **CancerType** (26 values): Maps to canonical cancer type names with 70+ aliases
- **EvidenceLevel** (5 values): A through E mapping to AMP/ASCO/CAP tiers
- **VariantType** (7 values): SNV, INDEL, CNV_AMP, CNV_DEL, FUSION, SPLICE, FRAMESHIFT
- **TherapyCategory** (6 values): TARGETED, IMMUNOTHERAPY, CHEMO, HORMONAL, COMBINATION, OTHER
- **TrialPhase** (5 values): Phase 1 through Phase 4 plus Early Phase 1
- **TrialStatus** (5 values): Recruiting, Active, Enrolling, Completed, Terminated
- **ActionabilityTier** (6 values): A, B, C, D, E, VUS
- **SearchStrategy** (3 values): broad, targeted, comparative
- **ExportFormat** (4 values): markdown, json, pdf, fhir

### 12.3 API Endpoints Summary

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/healthz` | GET | Health check |
| `/readyz` | GET | Readiness check (Milvus connection) |
| `/metrics` | GET | Prometheus metrics |
| `/v1/query` | POST | RAG query with optional streaming |
| `/v1/cases` | POST | Create case from VCF or structured input |
| `/v1/cases/{id}` | GET | Retrieve case by ID |
| `/v1/cases/{id}/mtb` | GET | Generate MTB packet |
| `/v1/therapy-rank` | POST | Rank therapies for a patient profile |
| `/v1/trial-match` | POST | Match patient to clinical trials |
| `/v1/compare` | POST | Comparative evidence retrieval |
| `/v1/export/markdown` | POST | Generate Markdown report |
| `/v1/export/json` | POST | Generate JSON export |
| `/v1/export/pdf` | POST | Generate PDF report |
| `/v1/export/fhir` | POST | Generate FHIR R4 bundle |

### 12.4 Port Configuration

| Service | Port | Protocol |
|---------|------|----------|
| Streamlit UI | 8526 | HTTP |
| FastAPI REST API | 8527 | HTTP |
| Milvus (shared) | 19530 | gRPC |

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Oncology Agent is a clinical decision support research tool for oncology. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
