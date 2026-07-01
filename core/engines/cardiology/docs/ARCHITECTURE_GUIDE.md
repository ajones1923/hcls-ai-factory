# Cardiology Intelligence Agent -- Architecture Guide

Part of the HCLS AI Factory Precision Intelligence Network -- one of three GPU-accelerated engines (Genomic Foundation Engine, Precision Intelligence Network, Therapeutic Discovery Engine) that compose the end-to-end precision medicine platform on NVIDIA DGX Spark.

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [System Diagram](#1-system-diagram)
2. [Component Interactions](#2-component-interactions)
3. [Data Flow](#3-data-flow)
4. [Collection Design Rationale](#4-collection-design-rationale)
5. [Workflow Engine](#5-workflow-engine)
6. [Risk Calculator Engine](#6-risk-calculator-engine)
7. [GDMT Optimizer](#7-gdmt-optimizer)
8. [Cross-Modal Engine](#8-cross-modal-engine)
9. [Query Expansion](#9-query-expansion)
10. [RAG Pipeline](#10-rag-pipeline)
11. [Agent Orchestrator](#11-agent-orchestrator)
12. [Cross-Agent Integration](#12-cross-agent-integration)
13. [Data Model Architecture](#13-data-model-architecture)

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
              | UI :8536   |      | :8126      |
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
   | (2,025 LOC) |    | (2,445 LOC)   |    |               |
   +------+------+    +-------+-------+    | Risk Calc     |
          |                    |           | (2,397 LOC)   |
          |                    |           |               |
          |                    |           | GDMT Opt      |
          |                    |           | (2,457 LOC)   |
          |                    |           |               |
          |                    |           | Cross-Modal   |
          |                    |           | (1,734 LOC)   |
          |                    |           +-------+-------+
          |                    |                    |
          +--------------------+--------------------+
                               |
                    +----------+----------+
                    |    RAG Engine       |
                    |    (1,589 LOC)      |
                    +----------+----------+
                               |
          +--------------------+--------------------+
          |                    |                    |
   +------+------+    +-------+-------+    +-------+-------+
   | Knowledge   |    | Milvus        |    | LLM           |
   | Graph       |    | Vector DB     |    | (Claude 4.6)  |
   | (1,431 LOC) |    | 13 Collections|    |               |
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
  | PubMed | Trials | ACC/AHA | FDA | SCMR  |
  +---+--------+--------+--------+--------+--+
      |        |        |        |        |
  +---v---+ +--v---+ +-v----+ +-v-----+ +v--------+
  |PubMed | |Trial | |Guide | |Device | |Imaging  |
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
CardioUI (Streamlit) ──> FastAPI Server ──> Agent Orchestrator
                                               |
                                 +-------------+-------------+
                                 |             |             |
                          QueryExpansion  WorkflowEngine  RiskCalculators
                                 |             |             |
                                 +------+------+      GDMTOptimizer
                                        |                    |
                                   RAGEngine          CrossModalEngine
                                        |                    |
                                 +------+------+             |
                                 |             |             |
                            Milvus DB    Knowledge      genomic_evidence
                            (12 cols)     Graph          (shared col)
```

### 2.2 Module Responsibilities

| Module | File | LOC | Responsibilities |
|--------|------|-----|-----------------|
| **Agent Orchestrator** | `src/agent.py` | 1,658 | Top-level coordination; receives queries, dispatches to sub-engines, assembles responses |
| **Query Expansion** | `src/query_expansion.py` | 2,025 | Entity extraction, synonym expansion, sub-question decomposition, strategy selection |
| **Workflow Engine** | `src/clinical_workflows.py` | 2,445 | 11 clinical workflow implementations with collection-specific search weights |
| **RAG Engine** | `src/rag_engine.py` | 1,589 | Multi-collection search, citation scoring, LLM prompt assembly and synthesis |
| **Risk Calculators** | `src/risk_calculators.py` | 2,397 | 6 validated scoring systems with published coefficients |
| **GDMT Optimizer** | `src/gdmt_optimizer.py` | 2,457 | 4-pillar HFrEF therapy optimization with contraindication checking |
| **Cross-Modal Engine** | `src/cross_modal.py` | 1,734 | Imaging-to-genomics trigger pattern matching and evidence retrieval |
| **Knowledge Graph** | `src/knowledge.py` | 1,431 | Static domain ontology (conditions, biomarkers, drugs, genes, imaging, guidelines) |
| **Collections** | `src/collections.py` | 1,226 | Milvus collection CRUD, schema definitions, index management |
| **Export** | `src/export.py` | 1,379 | PDF, CSV, JSON, FHIR report generation |
| **Metrics** | `src/metrics.py` | 537 | Prometheus instrumentation |
| **Scheduler** | `src/scheduler.py` | 612 | Periodic ingest scheduling |
| **Models** | `src/models.py` | 717 | 16 enums, 13 Pydantic models, 1 dataclass |

### 2.3 Interface Contracts

**Agent Orchestrator inputs/outputs:**
```
Input:  CardioQuery(question, workflow_type?, patient_context?)
Output: CardioResponse(answer, citations, risk_scores, workflow_results, cross_modal_triggers, confidence)
```

**Workflow Engine inputs/outputs:**
```
Input:  CardioQuery + SearchPlan
Output: List[WorkflowResult(workflow_type, findings, risk_scores, recommendations, guideline_references, severity, cross_modal_triggers)]
```

**Risk Calculator inputs/outputs:**
```
Input:  RiskScoreInput(score_type, demographics, vitals, labs, comorbidities) + Optional[extra_dict]
Output: RiskScoreResult(score_type, score_value, risk_category, interpretation, recommendations, guideline_reference)
```

**GDMT Optimizer inputs/outputs:**
```
Input:  Patient context (LVEF, NYHA, labs, current medications)
Output: GDMTRecommendation(ef_category, current_meds, recommendations, next_steps, guideline_references)
```

---

## 3. Data Flow

### 3.1 Query Processing Pipeline

```
Step 1: RECEIVE QUERY
  CardioQuery arrives via API or UI
  |
Step 2: QUERY EXPANSION (query_expansion.py)
  Extract entities (conditions, drugs, modalities)
  Expand synonyms via ENTITY_ALIASES
  Decompose into sub-questions
  Select search strategy (broad/targeted/comparative/clinical)
  Produce SearchPlan
  |
Step 3: WORKFLOW ROUTING (agent.py)
  Match SearchPlan.relevant_workflows to CardioWorkflowType enums
  If workflow_type specified in query, use that
  Otherwise, auto-detect from identified entities
  |
Step 4: PARALLEL EXECUTION (agent.py)
  |
  +---> 4a: MULTI-COLLECTION SEARCH (rag_engine.py)
  |     Embed expanded query
  |     Search 13 collections in parallel (top_k=5 each)
  |     Filter by score_threshold (0.4)
  |     Apply collection weights
  |     Deduplicate results
  |     Score citations (high/medium/standard)
  |
  +---> 4b: RISK CALCULATION (risk_calculators.py)
  |     Identify applicable calculators from patient_context
  |     Compute all applicable scores
  |     Generate risk categories and recommendations
  |
  +---> 4c: GDMT OPTIMIZATION (gdmt_optimizer.py)
  |     If HF workflow detected and patient_context includes LVEF
  |     Classify EF category
  |     Evaluate each GDMT pillar
  |     Check contraindications
  |     Generate titration recommendations
  |
  +---> 4d: WORKFLOW EXECUTION (clinical_workflows.py)
        Execute each matched workflow
        Workflow-specific collection weight adjustments
        Workflow-specific clinical logic
        Generate WorkflowResult per workflow
  |
Step 5: CROSS-MODAL TRIGGERS (cross_modal.py)
  Scan workflow results for imaging findings
  Match against IMAGING_TRIGGER_MAP
  Query genomic_evidence collection for matching variants
  Generate CrossModalTrigger objects
  |
Step 6: LLM SYNTHESIS (rag_engine.py)
  Assemble context from search results + risk scores + GDMT + workflows
  Send to Claude Sonnet 4.6 with cardiology system prompt
  Generate evidence-grounded clinical answer with inline citations
  |
Step 7: RESPONSE ASSEMBLY (agent.py)
  Package CardioResponse with:
  - Synthesized answer
  - Ranked citations
  - Risk score results
  - Workflow results
  - Cross-modal triggers
  - Confidence score
  |
Step 8: DELIVERY
  Return via API (JSON) or display in UI (formatted)
  Optional: Export as PDF/CSV/FHIR
```

### 3.2 Ingest Pipeline Flow

```
Step 1: SCHEDULER TRIGGER (scheduler.py)
  APScheduler fires at configured interval (default: weekly)
  Or manual trigger via scripts/run_ingest.py
  |
Step 2: PARSER SELECTION (src/ingest/__init__.py)
  Instantiate all 7 parsers
  |
Step 3: DATA ACQUISITION (per parser)
  PubMed: NCBI E-utilities API (cardiovascular query, max 5000 results)
  ClinicalTrials.gov: API (cardiovascular interventional trials)
  Guidelines: ACC/AHA/ESC document parsing
  Devices: FDA 510(k)/De Novo database
  Imaging/ECG/Hemodynamics: Structured data sources
  |
Step 4: PROCESSING (base.py)
  Chunk documents (configurable size + overlap)
  Generate BGE-small-en-v1.5 embeddings (batch_size=32)
  Deduplicate by content hash
  |
Step 5: INSERTION (base.py -> collections.py)
  Batch insert into target Milvus collection
  Include metadata (source, title, date, DOI)
  |
Step 6: VERIFICATION
  Log vector counts per collection
  Report errors/skipped documents
```

---

## 4. Collection Design Rationale

### 4.1 Why 12+1 Collections (Not One Giant Collection)

The multi-collection architecture was chosen over a single monolithic collection for five reasons:

**1. Semantic Precision**

The same words carry different meanings across cardiovascular subspecialties:
- "Gradient": pressure gradient (valvular), risk gradient (prevention), field gradient (MRI)
- "Flow": coronary flow reserve (cath), TIMI flow grade (angio), 4D flow (MRI)
- "Wall motion": stress echo finding, MRI cine assessment, VT substrate mapping

Collection-specific embedding spaces preserve these semantic boundaries.

**2. Relevance Weighting**

Different queries require different emphasis:
- "Manage HFrEF" -> weight cardio_heart_failure (0.25), cardio_guidelines (0.20)
- "Interpret LGE pattern" -> weight cardio_imaging (0.30), cardio_heart_failure (0.15)
- "Anticoagulation for AF" -> weight cardio_electrophysiology (0.25), cardio_guidelines (0.25)

A single collection would return mixed results without the ability to boost domain-specific content.

**3. Independent Lifecycle Management**

| Collection | Update Frequency | Source |
|-----------|-----------------|--------|
| cardio_literature | Weekly (PubMed) | Automated |
| cardio_guidelines | Every 2-5 years | Manual review |
| cardio_trials | Monthly | Automated |
| cardio_devices | Monthly | Automated |
| genomic_evidence | Shared, maintained by genomics pipeline | External |

Separate collections allow independent refresh, rollback, and quality control.

**4. Source Attribution**

Clinicians need to know the provenance of every citation:
- "This recommendation comes from the ACC/AHA 2022 HF Guideline" (cardio_guidelines)
- "This finding is based on the PARADIGM-HF trial" (cardio_trials)
- "This measurement reference comes from ASE guidelines" (cardio_imaging)

Collection names provide automatic source categorization.

**5. Scalability**

Individual collections can be independently:
- Reindexed (change from IVF_FLAT to HNSW for high-volume collections)
- Compacted (remove old/superseded entries)
- Partitioned (by year, by guideline society)
- Replicated (for read-heavy collections like cardio_guidelines)

### 4.2 Collection Weight Strategy

Weights sum to 1.0 and reflect clinical importance:

```
                    Higher Weight (0.10)          Lower Weight (0.04-0.06)
                    ┌────────────────────┐        ┌─────────────────┐
                    │ Literature   0.10  │        │ Devices    0.04 │
                    │ Imaging      0.10  │        │ Oncology   0.06 │
                    │ Heart Failure 0.10 │        │ Hemodynamics0.06│
                    │ Prevention   0.10  │        │ Genomic    0.03 │
                    │ Guidelines   0.10  │        └─────────────────┘
                    └────────────────────┘
                    Mid Weight (0.07-0.08)
                    ┌────────────────────┐
                    │ Trials       0.08  │
                    │ EP           0.08  │
                    │ Valvular     0.08  │
                    │ Interventional0.07 │
                    └────────────────────┘
```

**Dynamic weight adjustment**: Workflow detection can boost relevant collections by up to 2.5x while proportionally reducing others, maintaining the sum at 1.0.

---

## 5. Workflow Engine

### 5.1 Architecture

The workflow engine (`src/clinical_workflows.py`) implements a dispatcher pattern with 11 workflows:

```
CardioQuery + SearchPlan
        |
  WorkflowDispatcher.route(query, search_plan)
        |
  +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
  |     |     |     |     |     |     |     |     |     |     |     |
  CAD   HF   VHD   Arr   MRI   Stress Prev  Onc  ADHF  PMI  Myoc
  |     |     |     |     |     |     |     |     |     |     |     |
  +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
        |
  List[WorkflowResult]
```

### 5.2 Workflow Routing Logic

1. **Explicit routing**: If `CardioQuery.workflow_type` is set, route directly
2. **Auto-detection**: Match identified entities from SearchPlan against workflow triggers:
   - Mentions of "LVEF", "heart failure", "GDMT", "BNP" -> heart_failure workflow
   - Mentions of "ECG", "atrial fibrillation", "QTc", "antiarrhythmic" -> arrhythmia workflow
   - Mentions of "aortic stenosis", "mitral regurgitation", "valve" -> valvular_disease workflow
   - Mentions of "calcium score", "CAD-RADS", "stenosis", "angina" -> cad_assessment workflow
   - Mentions of "LGE", "T1 mapping", "cardiac MRI" -> cardiac_mri workflow
   - Mentions of "stress test", "Duke treadmill", "MPI" -> stress_test workflow
   - Mentions of "ASCVD", "statin", "cholesterol", "Lp(a)" -> preventive_risk workflow
   - Mentions of "cardiotoxicity", "doxorubicin", "trastuzumab", "GLS" -> cardio_oncology workflow
3. **Multi-workflow**: A single query can trigger multiple workflows (e.g., "HF patient with AF" triggers both heart_failure and arrhythmia)

### 5.3 Per-Workflow Collection Weights

Each workflow dynamically adjusts collection weights for optimal retrieval:

| Workflow | Boosted Collections | Reduced Collections |
|----------|-------------------|-------------------|
| CAD Assessment | imaging (0.20), prevention (0.18) | devices (0.02), oncology (0.02) |
| Heart Failure | heart_failure (0.25), guidelines (0.20) | imaging (0.05), interventional (0.03) |
| Valvular Disease | valvular (0.25), imaging (0.18) | EP (0.03), oncology (0.02) |
| Arrhythmia | EP (0.25), guidelines (0.20) | imaging (0.05), hemodynamics (0.03) |
| Cardiac MRI | imaging (0.30), heart_failure (0.15) | prevention (0.03), devices (0.02) |
| Stress Test | imaging (0.25), prevention (0.15) | hemodynamics (0.03), oncology (0.02) |
| Prevention | prevention (0.25), guidelines (0.20) | hemodynamics (0.02), devices (0.02) |
| Cardio-Oncology | oncology (0.30), guidelines (0.15) | interventional (0.02), hemodynamics (0.02) |
| Acute Decompensated HF | heart_failure (0.25), hemodynamics (0.25) | prevention (0.02), devices (0.02) |
| Post-MI | trials (0.20), guidelines (0.20), interventional (0.15) | oncology (0.02), hemodynamics (0.03) |
| Myocarditis/Pericarditis | imaging (0.25), literature (0.20), guidelines (0.15) | devices (0.02), prevention (0.02) |

---

## 6. Risk Calculator Engine

### 6.1 Calculator Architecture

```
RiskScoreInput + Optional[extra_dict]
        |
  RiskCalculatorEngine.calculate(score_type, input, extra)
        |
  +-----+-----+-----+-----+-----+-----+
  |     |     |     |     |     |     |
 ASCVD HEART CHA2  HAS  MAGGIC Euro
  |     |     |     |     |     |     |
  +-----+-----+-----+-----+-----+-----+
        |
  RiskScoreResult
```

### 6.2 Input Validation Strategy

Each calculator validates only the fields it needs:

| Calculator | Required Fields | Optional Fields |
|-----------|----------------|----------------|
| ASCVD | age, sex, total_cholesterol, hdl, systolic_bp, hypertension_treatment, diabetes, smoker | race |
| HEART | (via extra dict: history, ecg, age, risk_factors, troponin) | - |
| CHA2DS2-VASc | age, sex | chf, hypertension, diabetes, stroke, vascular_disease |
| HAS-BLED | - | hypertension, renal_disease, liver_disease, stroke, bleeding, labile_inr, age, alcohol, antiplatelet_nsaid |
| MAGGIC | age, sex, lvef, systolic_bp | nyha_class, bmi, creatinine, diabetes, bb_use, acei_arb_use, COPD, HF_duration, smoker |
| EuroSCORE II | age, sex | (28 factors via extra dict) |

### 6.3 Auto-Detection: calculate_all_applicable

The `calculate_all_applicable` method iterates through all 6 calculators and computes scores for any calculator where sufficient data is present. This enables a single API call to return multiple risk assessments simultaneously.

```python
engine = RiskCalculatorEngine()
results = engine.calculate_all_applicable(risk_input, extra=extra_dict)
# Returns: List[RiskScoreResult] with all calculable scores
```

---

## 7. GDMT Optimizer

### 7.1 Optimizer Architecture

```
Patient Context (LVEF, NYHA, labs, medications)
        |
  GDMTOptimizer.optimize(patient_context)
        |
  EF Classification (HFrEF/HFmrEF/HFpEF/HFimpEF)
        |
  +-----+-----+-----+-----+
  |     |     |     |     |
  BB   ARNI   MRA  SGLT2i
  |     |     |     |     |
  +-----+-----+-----+-----+
  Per-pillar evaluation:
  1. Current status (not_started/initiated/uptitrating/at_target/contraindicated/intolerant)
  2. Contraindication check (K+, eGFR, BP, HR, etc.)
  3. Titration recommendation (specific dose + timeline)
  4. Monitoring schedule (labs at 1 week, 4 weeks, quarterly)
        |
  GDMTRecommendation
```

### 7.2 Contraindication Logic

Each pillar has specific contraindication checks:

```
Beta-Blocker:
  IF HR < 50 -> contraindicated (severe bradycardia)
  IF SBP < 90 -> contraindicated (hypotension)
  IF acute_decompensation -> hold (do not initiate during acute HF)
  IF severe_reactive_airway -> contraindicated

ARNI:
  IF K+ > 5.5 -> contraindicated
  IF eGFR < 20 -> contraindicated
  IF SBP < 100 -> contraindicated
  IF angioedema_history -> contraindicated
  IF on_acei_within_36h -> contraindicated (washout required)

MRA:
  IF K+ > 5.0 -> hold
  IF K+ > 5.5 -> contraindicated
  IF eGFR < 30 -> contraindicated

SGLT2i:
  IF eGFR < 20 -> do not initiate (can continue if already started)
  IF type_1_diabetes -> contraindicated (DKA risk)
```

### 7.3 Titration Sequencing

The optimizer recommends a specific sequence for medication changes to minimize hemodynamic instability:

1. Start or uptitrate only one pillar at a time (2-week intervals)
2. SGLT2i can be started simultaneously with other pillars (no titration needed)
3. Switch ACEi to ARNI before uptitrating beta-blocker (if BP allows)
4. Add MRA after ARNI/ACEi is at target dose
5. Check labs (K+, creatinine) after each MRA initiation or dose change

---

## 8. Cross-Modal Engine

### 8.1 Trigger Architecture

```
WorkflowResult.findings
        |
  CrossModalEngine.scan(workflow_results)
        |
  Pattern Matching against IMAGING_TRIGGER_MAP (18 patterns)
        |
  For each matched trigger:
        |
  +-----+-----+
  |           |
  Generate    Query genomic_evidence
  gene_panel  collection for variants
  |           |
  +-----+-----+
        |
  CrossModalTrigger(trigger_source, finding, gene_panel, conditions, rationale)
```

### 8.2 Pattern Matching Strategy

The engine uses a two-stage matching approach:

**Stage 1: Rule-based keyword matching**
- Scans findings for trigger keywords (e.g., "LVH", "wall thickness", "LGE", "aortic root")
- Fast, deterministic, high precision

**Stage 2: Semantic similarity matching**
- For ambiguous findings, embeds the finding text and compares against trigger descriptions
- Uses same BGE-small-en-v1.5 embedding model
- Score threshold: 0.6 for trigger activation

### 8.3 Genomic Evidence Integration

When a trigger fires, the engine queries the shared `genomic_evidence` collection:

```python
# Pseudocode for genomic evidence lookup
for gene in trigger.gene_panel:
    results = milvus.search(
        collection="genomic_evidence",
        query=f"pathogenic variant {gene} {trigger.conditions}",
        top_k=3,
        score_threshold=0.5
    )
    trigger.genomic_results.extend(results)
```

This enables the system to report not just "genetic testing is recommended" but also "here are known pathogenic variants in MYH7 associated with HCM from the genomic evidence base."

---

## 9. Query Expansion

### 9.1 Expansion Pipeline

```
Raw Query: "What GDMT changes for an HFrEF patient with K+ of 5.3?"
        |
  Step 1: Entity Extraction
  - Conditions: [HFrEF]
  - Drugs: [GDMT]
  - Labs: [K+ 5.3]
        |
  Step 2: Synonym Expansion
  - "HFrEF" -> "heart failure with reduced ejection fraction, systolic heart failure"
  - "GDMT" -> "guideline-directed medical therapy"
  - "K+" -> "potassium, serum potassium"
        |
  Step 3: Sub-Question Decomposition
  - Q1: "What are the 4 pillars of GDMT for HFrEF?"
  - Q2: "How does hyperkalemia (K+ 5.3) affect GDMT?"
  - Q3: "Which GDMT medications should be held or adjusted for K+ >5.0?"
        |
  Step 4: Strategy Selection
  - Strategy: "clinical" (patient-specific medication management)
        |
  Step 5: Workflow Routing
  - Workflows: [heart_failure]
  - Boost: cardio_heart_failure, cardio_guidelines
        |
  Output: SearchPlan
```

### 9.2 Entity Alias Resolution

The expansion module leverages the 167-entry `ENTITY_ALIASES` dictionary (18 synonym maps):

```
Input aliases resolved:
  "AFib" -> "atrial fibrillation"
  "BB" -> "beta-blocker"
  "ARNI" -> "angiotensin receptor-neprilysin inhibitor"
  "EF" -> "ejection fraction"
  "CAC" -> "coronary artery calcium"
  "TAVR" -> "transcatheter aortic valve replacement"
  "GLS" -> "global longitudinal strain"
```

This ensures that queries using standard clinical abbreviations match content using full terminology, and vice versa.

---

## 10. RAG Pipeline

### 10.1 RAG Architecture

```
SearchPlan (from Query Expansion)
        |
  Embedding (BGE-small-en-v1.5, 384-dim)
        |
  Multi-Collection Search
  +-----+-----+-----+-----+-----+-----+-----+
  |lit  |trial|img  |ep   |hf   |valv |prev |
  |top5 |top5 |top5 |top5 |top5 |top5 |top5 |
  +-----+-----+-----+-----+-----+-----+-----+
  |intv |onc  |dev  |gdl  |hemo |gen  |
  |top5 |top5 |top5 |top5 |top5 |top5 |
  +-----+-----+-----+-----+-----+-----+
        |
  Score Filtering (threshold 0.4)
        |
  Weight Application (per-collection weights)
        |
  Deduplication (content hash)
        |
  Citation Scoring
  - High confidence: score >= 0.75
  - Medium confidence: score >= 0.60
  - Standard: score < 0.60
        |
  Context Assembly
  - Rank by weighted score
  - Include source metadata
  - Format for LLM prompt
        |
  LLM Synthesis (Claude Sonnet 4.6)
  - System: Cardiology domain expert
  - User: Clinical question + retrieved context
  - Instructions: Cite sources, follow guidelines, flag uncertainties
        |
  Response Parsing
  - Extract inline citations
  - Map to source documents
  - Calculate confidence score
```

### 10.2 LLM Prompt Structure

```
System Prompt:
  "You are a cardiovascular clinical decision support system. Your responses
   must be evidence-based, citing specific ACC/AHA/ESC guidelines with
   recommendation class and evidence level. Always include relevant risk
   scores, GDMT status, and cross-modal considerations. Flag any areas
   of clinical uncertainty."

User Message:
  "[Clinical Question]

   --- Retrieved Evidence ---
   [Source 1: cardio_guidelines, score 0.89]
   In patients with HFrEF, ARNI is recommended to reduce morbidity and mortality.
   Class I, Level A. (2022 AHA/ACC/HFSA HF Guideline)

   [Source 2: cardio_heart_failure, score 0.85]
   ...

   [Source N: cardio_trials, score 0.72]
   ...

   --- Risk Scores ---
   MAGGIC: 24 (22.3% 1-year mortality)

   --- GDMT Status ---
   BB: uptitrating (carvedilol 12.5mg BID -> target 25mg BID)
   ARNI: not started
   MRA: not started (K+ 5.3 -> hold until K+ <5.0)
   SGLT2i: not started

   --- Patient Context ---
   Age 62, Female, LVEF 28%, NYHA III"
```

---

## 11. Agent Orchestrator

### 11.1 Orchestration Flow

The agent orchestrator (`src/agent.py`, 1,658 lines) is the central coordinator:

```python
class CardiologyAgent:
    def __init__(self):
        self.query_expander = QueryExpander()
        self.workflow_engine = WorkflowEngine()
        self.rag_engine = RAGEngine()
        self.risk_calculator = RiskCalculatorEngine()
        self.gdmt_optimizer = GDMTOptimizer()
        self.cross_modal_engine = CrossModalEngine()
        self.export_system = ExportSystem()

    async def process_query(self, query: CardioQuery) -> CardioResponse:
        # 1. Expand query
        search_plan = self.query_expander.expand(query)

        # 2. Execute workflows
        workflow_results = await self.workflow_engine.execute(query, search_plan)

        # 3. Calculate risk scores
        risk_scores = self.risk_calculator.calculate_all_applicable(
            query.patient_context
        )

        # 4. Optimize GDMT (if applicable)
        gdmt_rec = self.gdmt_optimizer.optimize(query.patient_context)

        # 5. Check cross-modal triggers
        triggers = self.cross_modal_engine.scan(workflow_results)

        # 6. RAG synthesis
        answer = await self.rag_engine.synthesize(
            query, search_plan, workflow_results, risk_scores, gdmt_rec, triggers
        )

        # 7. Assemble response
        return CardioResponse(
            answer=answer.text,
            citations=answer.citations,
            risk_scores=risk_scores,
            workflow_results=workflow_results,
            cross_modal_triggers=triggers,
            confidence=answer.confidence
        )
```

### 11.2 Error Handling Strategy

The orchestrator implements graceful degradation:

1. **Milvus unavailable**: Returns error with search-only mode recommendation
2. **LLM unavailable**: Returns search results without synthesis (search-only mode)
3. **Risk calculator error**: Logs warning, returns available scores, skips failed calculator
4. **GDMT error**: Logs warning, returns response without GDMT section
5. **Cross-modal error**: Logs warning, returns response without triggers
6. **Timeout**: Returns partial results with timeout indicator

### 11.3 Conversation Memory

The agent maintains a sliding window of conversation context:

```
MAX_CONVERSATION_CONTEXT = 3  # Configurable

conversation_history = deque(maxlen=3)

# Each entry stores:
{
    "query": CardioQuery,
    "response": CardioResponse,
    "timestamp": datetime
}
```

This enables contextual follow-up questions like "What about adding spironolactone?" after an initial HF query.

### 11.4 New Capabilities (Post-Launch)

Several capabilities have been added since the initial architecture:

- **Rate Limiting**: 100 requests/minute per IP address to prevent abuse
- **API Authentication**: X-API-Key header required for all API endpoints
- **Conversation Memory Persistence**: File-based conversation memory with 24-hour TTL (persists across API restarts)
- **Cross-Agent Event Publishing**: Server-Sent Events (SSE) via `/v1/cardio/events` endpoint with 5 event types for integration with other HCLS AI Factory agents
- **Knowledge Base Versioning**: `/v1/cardio/knowledge-version` endpoint reports current knowledge base version and last update timestamp
- **Input Validation**: All 11 workflows now validate input parameters before execution
- **Lazy Initialization**: Risk calculators work without Milvus connectivity (pure computation)
- **1,966 tests passing** across all modules

---

## 12. Cross-Agent Integration

### 12.1 Platform Context: 11 Intelligence Agents

The Cardiology Intelligence Agent is one of 11 specialized intelligence agents within the HCLS AI Factory Precision Intelligence Network:

| # | Agent | UI Port | API Port | Domain |
|---|-------|---------|----------|--------|
| 1 | Precision Biomarker | 8528 | 8529 | Biomarker interpretation and trends |
| 2 | Precision Oncology | 8525 | 8526 | Cancer genomics and treatment |
| 3 | CAR-T Intelligence | 8521 | 8522 | CAR-T cell therapy |
| 4 | Medical Imaging | 8523 | 8524 | Multi-modal imaging analysis |
| 5 | Precision Autoimmune | 8531 | 8532 | Autoimmune disease analysis |
| 6 | Cardiology Intelligence | 8536 | 8126 | Cardiovascular decision support |
| 7 | Clinical Trial Intelligence | 8538 | 8537 | Trial matching and eligibility |
| 8 | Neurology Intelligence | 8529 | 8528 | Neurological decision support |
| 9 | Rare Disease Intelligence | 8135 | 8134 | Rare disease diagnosis |
| 10 | Pharmacogenomics (PGx) | 8541 | 8540 | Drug-gene interactions |
| 11 | Pediatric Oncology | 8543 | 8542 | Pediatric cancer intelligence |

### 12.2 Cross-Agent Calls

The Cardiology Intelligence Agent calls 4 sibling agents for integrated assessments:

```
Cardiology Agent
  |
  +--> Oncology Agent (8526)
  |      Cardiotoxicity correlation: links CTRCD findings to chemotherapy regimens
  |      Requests: active treatment protocols, cumulative anthracycline dose
  |
  +--> Clinical Trial Agent (8537)
  |      Matches cardiology patients to active cardiovascular trials
  |      Requests: trial eligibility for HF devices, novel GDMT, cardio-onc interventions
  |
  +--> Biomarker Agent (8529)
  |      Longitudinal biomarker trending for troponin, BNP, GLS
  |      Requests: biomarker trajectories, threshold alerts
  |
  +--> Imaging Agent (8524)
         Cross-modal imaging correlation for echo, CMR, CT, nuclear
         Requests: imaging findings, measurement trends, tissue characterization
```

### 12.3 Integrated Assessment Endpoint

```
POST /v1/cardio/integrated-assessment

Request:
  patient_id: str
  include_agents: List[str]  # ["oncology", "trial", "biomarker", "imaging"]
  clinical_context: dict

Response:
  cardiology_assessment: CardioResponse
  cross_agent_findings: List[CrossAgentResult]
  pediatric_cardiotoxicity: Optional[PediatricCardiotoxicityAssessment]
  integrated_recommendations: List[str]
```

### 12.4 Pediatric Cardiotoxicity Assessment

The `pediatric_cardiotoxicity_assessment()` function provides specialized analysis for pediatric oncology patients receiving cardiotoxic chemotherapy:

- **Anthracycline cardiotoxicity in children**: 5-25% late-onset incidence, with risk increasing years after treatment completion
- **Cumulative dose thresholds**: 250 mg/m2 (increased surveillance), 450 mg/m2 (high-risk threshold)
- **Dexrazoxane cardioprotection**: Recommended for cumulative doxorubicin >=250 mg/m2 in pediatric patients per COG guidelines
- **Pediatric echo monitoring**: Serial echocardiography with GLS tracking at baseline, mid-treatment, end-of-treatment, and long-term follow-up (1, 2, 5, 10 years)
- **COG LTFU guidelines**: Children's Oncology Group Long-Term Follow-Up guidelines for cardiac surveillance in childhood cancer survivors

---

## 13. Data Model Architecture

### 12.1 Model Hierarchy

```
                    CardioResponse (top-level output)
                    /        |         \           \
              answer    citations    risk_scores   cross_modal_triggers
                |           |            |              |
              str    List[Dict]   List[RiskScoreResult]  List[CrossModalTrigger]
                                        |
                            (from RiskScoreInput)

                    WorkflowResult (per-workflow output)
                    /        |         \           \
              findings  risk_scores  recommendations  severity
                |           |            |              |
           List[str]  List[RSResult]  List[str]    SeverityLevel

                    GDMTRecommendation
                    /        |         \
           ef_category  current_meds  recommendations
                |           |            |
         EFCategory  List[GDMTMedication]  List[str]
```

### 12.2 Enum Design Rationale

Enums enforce type safety and valid value sets:

- **CardioWorkflowType** (12 values): Prevents invalid workflow routing (9 original + 3 new: acute_decompensated_hf, post_mi, myocarditis_pericarditis)
- **RiskScoreType** (6 values): Maps to calculator implementations
- **SeverityLevel** (5 values): Standardized urgency classification
- **HeartFailureClass/Stage**: NYHA I-IV and ACC/AHA A-D
- **EjectionFractionCategory** (4 values): HFrEF/HFmrEF/HFpEF/HFimpEF
- **GDMTPillar** (4 values): Maps to GDMT optimization pillars
- **GDMTStatus** (6 values): Tracks medication titration state
- **GuidelineClass/EvidenceLevel**: ACC/AHA recommendation classification

### 12.3 Pydantic Validation

All input models use Pydantic field validators:
- `RiskScoreInput.age`: ge=18, le=120
- `RiskScoreInput.systolic_bp`: ge=60, le=300
- `RiskScoreInput.lvef`: ge=5.0, le=90.0
- `RiskScoreInput.sex`: validated to "male" or "female" via custom validator
- `CardioQuery.question`: min_length=1

These validators prevent invalid clinical data from reaching calculator logic, ensuring reproducible results.
