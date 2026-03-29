# Precision Biomarker Intelligence Agent -- Architecture Guide

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [System Diagram](#1-system-diagram)
2. [Component Interactions](#2-component-interactions)
3. [Data Flow](#3-data-flow)
4. [Collection Design Rationale](#4-collection-design-rationale)
5. [Pharmacogenomic Engine](#5-pharmacogenomic-engine)
6. [Biological Age Engine](#6-biological-age-engine)
7. [Disease Trajectory Engine](#7-disease-trajectory-engine)
8. [Genotype Adjustment Engine](#8-genotype-adjustment-engine)
9. [Critical Value Engine](#9-critical-value-engine)
10. [Discordance Detector](#10-discordance-detector)
11. [RAG Pipeline](#11-rag-pipeline)
12. [Agent Orchestrator](#12-agent-orchestrator)
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
              | UI :8533   |      | :8529      |
              +-----+------+      +------+-----+
                    |                     |
                    +----------+----------+
                               |
                    +----------+----------+
                    |  BiomarkerAgent    |
                    |  (src/agent.py)    |
                    +----------+----------+
                               |
       +-----------------------+------------------------+
       |            |          |          |              |
+------+------+ +--+------+ +-+------+ +-+--------+ +--+--------+
| Pharmacoge- | |Bio Age  | |Disease | |Genotype  | |Critical   |
| nomic       | |Calc     | |Traject | |Adjuster  | |Value Eng  |
| Mapper      | |(408 LOC)| |(1,421) | |(1,225)   | |(179 LOC)  |
| (1,503 LOC) | +---------+ +--------+ +----------+ +-----------+
+------+------+                                       |
       |            +----------------------------------+
       |            |
+------+------+ +---+--------+
| Discordance | | Lab Range  |
| Detector    | | Interpreter|
| (299 LOC)   | | (221 LOC)  |
+------+------+ +------------+
       |
       +--------------------+
                            |
                 +----------+----------+
                 |    RAG Engine       |
                 |    (573 LOC)        |
                 +----------+----------+
                            |
       +--------------------+--------------------+
       |                    |                    |
+------+------+    +-------+-------+    +-------+-------+
| Knowledge   |    | Milvus        |    | LLM           |
| Graph       |    | Vector DB     |    | (Claude 4.6)  |
| (1,326 LOC) |    | 14 Collections|    |               |
| 6 domains   |    |               |    |               |
+-------------+    +-------+-------+    +---------------+
                            |
                 +----------+----------+
                 |  etcd    |  MinIO   |
                 +----------+----------+
```

### 1.2 Multi-Format Export Pipeline

```
  +------------------------------------------+
  |         AnalysisResult                    |
  +------------------------------------------+
      |        |        |        |        |
  +---v---+ +--v---+ +-v----+ +-v-----+ +v--------+
  |Markdwn| |  PDF | | FHIR | |  CSV  | |  JSON   |
  |Report | |Report| |  R4  | |Export | |Export   |
  |12 sect| |ReportLab| Bundle| |Flat  | |Struct  |
  +---+---+ +--+---+ +-+----+ +-+-----+ ++--------+
      |        |        |        |        |
      +--------+--------+--------+--------+
                         |
                  +------v------+
                  | Translation |
                  | Engine      |
                  | 7 languages |
                  +-------------+
```

---

## 2. Component Interactions

### 2.1 Component Dependency Graph

```
BiomarkerUI (Streamlit) ──> FastAPI Server ──> BiomarkerAgent
                                                     |
                          +--------------------------+----------------------------+
                          |              |           |           |                |
                   PGxMapper    BioAgeCalc   DiseaseTrajectory  GenotypeAdjuster  CriticalValues
                          |              |           |           |                |
                          +------+-------+-----------+-----------+-------+--------+
                                 |                                       |
                            DiscordanceDetector                    LabRangeInterpreter
                                 |                                       |
                                 +-------------------+-------------------+
                                                     |
                                                RAGEngine
                                                     |
                                 +-------------------+-------------------+
                                 |                   |                   |
                            Milvus DB          Knowledge Graph     genomic_evidence
                            (14 cols)          (6 domains)         (shared col)
```

### 2.2 Module Responsibilities

| Module | File | LOC | Responsibilities |
|--------|------|-----|-----------------|
| **BiomarkerAgent** | `src/agent.py` | 610 | Plan-analyze-search-synthesize-report loop |
| **PharmacogenomicMapper** | `src/pharmacogenomics.py` | 1,503 | CPIC-guided phenotyping across 14 pharmacogenes |
| **DiseaseTrajectoryAnalyzer** | `src/disease_trajectory.py` | 1,421 | Pre-symptomatic detection across 9 disease categories |
| **GenotypeAdjuster** | `src/genotype_adjustment.py` | 1,225 | Genotype-based reference range adjustments for 7 genes |
| **BiologicalAgeCalculator** | `src/biological_age.py` | 408 | PhenoAge and GrimAge surrogate calculations |
| **CriticalValueEngine** | `src/critical_values.py` | 179 | Life-threatening threshold detection (21 rules) |
| **DiscordanceDetector** | `src/discordance_detector.py` | 299 | Cross-biomarker anomaly detection |
| **LabRangeInterpreter** | `src/lab_range_interpreter.py` | 221 | Standard vs optimal range interpretation |
| **BiomarkerRAGEngine** | `src/rag_engine.py` | 573 | Multi-collection RAG with parallel search |
| **ReportGenerator** | `src/report_generator.py` | 993 | 12-section clinical reports |
| **ExportPipeline** | `src/export.py` | 1,392 | PDF/FHIR/CSV/JSON/Markdown export |
| **TranslationEngine** | `src/translation.py` | 217 | Multi-language medical terminology (7 languages) |
| **Knowledge Graph** | `src/knowledge.py` | 1,326 | 6 disease domain knowledge graphs |
| **Collections** | `src/collections.py` | 1,391 | Milvus collection management (14 schemas) |
| **Models** | `src/models.py` | 786 | 14 collection models, 8+ analysis models, 7 enums |
| **Audit** | `src/audit.py` | 83 | HIPAA-compliant audit logging |

### 2.3 Interface Contracts

**BiomarkerAgent inputs/outputs:**
```
Input:  PatientProfile(demographics, biomarkers, genotypes, medications)
Output: AnalysisResult(pgx_results, biological_age, disease_trajectories, genotype_adjustments, critical_values, discordances, report)
```

**PharmacogenomicMapper inputs/outputs:**
```
Input:  GenotypePanelInput(gene_genotypes: Dict[gene, star_alleles])
Output: PGxResult(phenotype, drug_recommendations: List[DrugRecommendation], alert_level)
```

**BiologicalAgeCalculator inputs/outputs:**
```
Input:  BiomarkerPanel(chronological_age, albumin, creatinine, glucose, CRP, ...)
Output: BiologicalAgeResult(pheno_age, grim_age, acceleration, percentile, contributing_biomarkers)
```

---

## 3. Data Flow

### 3.1 Patient Analysis Pipeline

```
Step 1: RECEIVE PATIENT PROFILE
  PatientProfile arrives via API or UI
  Demographics, biomarker values, genotype data, medications
  |
Step 2: PLAN (agent.py)
  Identify relevant topics, disease areas, and modules
  Select analysis engines based on available patient data
  Determine which collections to prioritize
  |
Step 3: ANALYZE (parallel execution)
  +-- BiologicalAgeCalculator: PhenoAge + GrimAge computation
  +-- DiseaseTrajectoryAnalyzer: 9 disease category risk assessment
  +-- PharmacogenomicMapper: 14 pharmacogene phenotyping
  +-- GenotypeAdjuster: 7 modifier gene range adjustments
  +-- CriticalValueEngine: 21-rule threshold check
  +-- DiscordanceDetector: cross-biomarker pattern analysis
  +-- LabRangeInterpreter: standard vs optimal range interpretation
  |
Step 4: SEARCH (rag_engine.py)
  Multi-collection RAG across 14 collections
  ThreadPoolExecutor parallel search (max_workers=6)
  Weighted scoring with per-collection weights (sum = 1.0)
  Score threshold filtering (minimum 0.4 cosine similarity)
  |
Step 5: SYNTHESIZE (rag_engine.py)
  Claude LLM synthesis with domain knowledge injection
  Grounded responses with evidence citations
  Conversation memory for follow-up questions
  |
Step 6: REPORT (report_generator.py + export.py)
  12-section clinical report generation
  Multi-format export: Markdown, PDF, FHIR R4, CSV, JSON
  Multi-language translation (7 languages)
  |
  v
AnalysisResult + Clinical Report
```

### 3.2 Data Ingestion Flow

```
Reference Data (18 JSON files)
    |
    v
[seed_all.py]
    |  Read JSON -> Instantiate Pydantic models
    |  Call to_embedding_text() on each model
    |  Embed with BGE-small-en-v1.5 (batch_size=32)
    |  Insert into Milvus collections
    |
    v
14 Milvus Collections (indexed, searchable)
```

---

## 4. Collection Design Rationale

### 4.1 Collection Inventory

The agent manages 14 specialized vector collections. All use COSINE similarity
with IVF_FLAT indexing (nlist=1024, nprobe=16) and 384-dimensional vectors
from BGE-small-en-v1.5.

| # | Collection Name | Description | Weight | Seed Records |
|---|----------------|-------------|--------|-------------|
| 1 | biomarker_reference | Reference biomarker definitions with standard/optimal ranges | 0.12 | 624 |
| 2 | biomarker_genetic_variants | Genetic variants affecting biomarker levels | 0.11 | 42 |
| 3 | biomarker_pgx_rules | CPIC pharmacogenomic dosing rules | 0.10 | 80 |
| 4 | biomarker_disease_trajectories | Disease progression trajectories with intervention windows | 0.10 | 60 |
| 5 | biomarker_clinical_evidence | Published clinical evidence with PubMed references | 0.09 | 240 |
| 6 | genomic_evidence | Shared VCF-derived genomic variants (read-only) | 0.08 | -- |
| 7 | biomarker_drug_interactions | Gene-drug interactions beyond PGx | 0.07 | 45 |
| 8 | biomarker_aging_markers | Epigenetic aging clock marker data (PhenoAge/GrimAge) | 0.07 | 30 |
| 9 | biomarker_nutrition | Genotype-aware nutrition guidelines | 0.05 | 50 |
| 10 | biomarker_genotype_adjustments | Genotype-based reference range adjustments | 0.05 | 35 |
| 11 | biomarker_monitoring | Condition-specific monitoring protocols | 0.05 | 40 |
| 12 | biomarker_critical_values | Critical threshold definitions | 0.04 | 21 |
| 13 | biomarker_discordance_rules | Cross-biomarker discordance patterns | 0.04 | 20 |
| 14 | biomarker_aj_carrier_screening | Ashkenazi Jewish carrier screening panel (10 genes) | 0.03 | 10 |
| | **Total** | | **1.00** | **1,297** |

### 4.2 Collection Weight Distribution

```
biomarker_reference          ████████████████        0.12
biomarker_genetic_variants   ███████████████         0.11
biomarker_pgx_rules          ██████████████          0.10
biomarker_disease_traject.   ██████████████          0.10
biomarker_clinical_evidence  █████████████           0.09
genomic_evidence             ████████████            0.08
biomarker_drug_interactions  ███████████             0.07
biomarker_aging_markers      ███████████             0.07
biomarker_nutrition          ████████                0.05
biomarker_genotype_adjust.   ████████                0.05
biomarker_monitoring         ████████                0.05
biomarker_critical_values    ██████                  0.04
biomarker_discordance_rules  ██████                  0.04
biomarker_aj_carrier_screen  █████                   0.03
                                               Sum: 1.00
```

### 4.3 Key Collection Schemas

**biomarker_reference** -- Core biomarker definitions (624 records):

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| name | VARCHAR(100) | Display name (e.g., "HbA1c") |
| unit | VARCHAR(20) | e.g., "%", "mg/dL" |
| category | VARCHAR(30) | CBC, CMP, Lipids, Thyroid |
| ref_range_min | FLOAT | Standard lower bound |
| ref_range_max | FLOAT | Standard upper bound |
| text_chunk | VARCHAR(3000) | Text for embedding |
| clinical_significance | VARCHAR(2000) | Clinical interpretation |
| epigenetic_clock | VARCHAR(50) | PhenoAge/GrimAge coefficient |
| genetic_modifiers | VARCHAR(500) | Comma-separated modifier genes |

**biomarker_clinical_evidence** -- Published evidence (240 records):

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| title | VARCHAR(500) | Publication title |
| pmid | VARCHAR(20) | PubMed ID |
| journal | VARCHAR(200) | Journal name |
| year | INT64 | Publication year |
| biomarker_focus | VARCHAR(200) | Primary biomarker(s) discussed |
| evidence_level | VARCHAR(20) | Evidence grade |
| text_chunk | VARCHAR(3000) | Text for embedding |

**biomarker_pgx_rules** -- Pharmacogenomic dosing rules:

| Field | Type | Notes |
|-------|------|-------|
| id (PK) | VARCHAR(100) | Primary key |
| embedding | FLOAT_VECTOR | 384-dim |
| gene | VARCHAR(50) | Pharmacogene (CYP2D6, etc.) |
| star_alleles | VARCHAR(100) | Star allele combination |
| drug | VARCHAR(100) | Drug name |
| phenotype | VARCHAR(30) | MetabolizerPhenotype enum |
| cpic_level | VARCHAR(5) | CPIC evidence level (1A-3) |
| recommendation | VARCHAR(2000) | Dosing recommendation |
| text_chunk | VARCHAR(3000) | Text for embedding |

---

## 5. Pharmacogenomic Engine

### 5.1 Architecture

```
GenotypePanelInput (14 pharmacogenes x star alleles)
    |
    v
[Phenotype Classification]
    |   Map star allele combinations to metabolizer phenotypes
    |   CPIC standard terminology: Normal, Intermediate, Poor, Ultra-rapid, Rapid
    |
    v
[Drug Recommendation Generation]
    |   For each gene-phenotype combination:
    |   - Identify affected drugs from CPIC guidelines
    |   - Determine dosing action (standard, reduce, avoid, etc.)
    |   - Set alert level (info, warning, critical)
    |
    v
PGxResult (per-gene phenotype + per-drug recommendations)
```

### 5.2 Covered Pharmacogenes (14)

| Gene | Description | CPIC Level | Key Drugs |
|------|-------------|------------|-----------|
| CYP2D6 | Metabolizes ~25% of drugs | 1A | Codeine, tramadol, tamoxifen |
| CYP2C19 | Clopidogrel, PPIs, antidepressants | 1A | Clopidogrel, omeprazole |
| CYP2C9 | Warfarin, NSAIDs, phenytoin | 1A | Warfarin, celecoxib |
| CYP3A5 | Tacrolimus metabolism | 1A | Tacrolimus |
| SLCO1B1 | Statin hepatic uptake transporter | 1A | Simvastatin, atorvastatin |
| VKORC1 | Warfarin target sensitivity | 1A | Warfarin |
| MTHFR | Folate metabolism, homocysteine | Informational | Methotrexate |
| TPMT | Thiopurine metabolism | 1A | Azathioprine, mercaptopurine |
| DPYD | Fluoropyrimidine metabolism | 1A | 5-FU, capecitabine |
| HLA-B*57:01 | Abacavir hypersensitivity | 1A | Abacavir |
| G6PD | Oxidative drug hemolysis | 1A | Rasburicase, dapsone |
| HLA-B*58:01 | Allopurinol hypersensitivity | 1A | Allopurinol |
| UGT1A1 | Irinotecan, atazanavir | 1A | Irinotecan |
| CYP4F2 | Vitamin K recycling | Informational | Warfarin |

### 5.3 Drug Recommendation Actions

| Action | Description | Alert Level |
|--------|-------------|-------------|
| STANDARD_DOSING | No adjustment needed | INFO |
| DOSE_REDUCTION | Reduce dose per CPIC | WARNING |
| DOSE_ADJUSTMENT | Adjust based on phenotype | WARNING |
| CONSIDER_ALTERNATIVE | Prefer alternative drug | WARNING |
| AVOID | Do not prescribe | CRITICAL |
| CONTRAINDICATED | Absolute contraindication | CRITICAL |

---

## 6. Biological Age Engine

### 6.1 PhenoAge (Levine 2018)

Uses 9 blood biomarkers plus chronological age:

| Biomarker | Unit | PhenoAge Weight |
|-----------|------|----------------|
| Albumin | g/dL | Negative (protective) |
| Creatinine | mg/dL | Positive |
| Glucose | mg/dL | Positive |
| C-reactive protein (log) | mg/L | Positive |
| Lymphocyte % | % | Negative (protective) |
| Mean cell volume | fL | Positive |
| Red cell distribution width | % | Positive |
| Alkaline phosphatase | U/L | Positive |
| White blood cell count | 10^3/uL | Positive |
| Chronological age | years | Positive |

### 6.2 GrimAge Surrogate

Uses a subset of PhenoAge biomarkers with smoking pack-years and additional aging markers to estimate the GrimAge epigenetic clock.

### 6.3 Output

```
BiologicalAgeResult:
  - pheno_age: 58.3 years
  - grim_age: 56.1 years
  - chronological_age: 52 years
  - acceleration: +6.3 years (PhenoAge)
  - percentile: 78th (aging faster than 78% of age-matched cohort)
  - contributing_biomarkers: [CRP (high impact), RDW (moderate), glucose (moderate)]
  - recommendations: ["Reduce CRP through anti-inflammatory interventions", ...]
```

---

## 7. Disease Trajectory Engine

### 7.1 Covered Disease Categories (9)

| Category | Key Biomarkers | Stages |
|----------|---------------|--------|
| Diabetes | HbA1c, fasting glucose, HOMA-IR | Pre-diabetic -> T2DM -> Complications |
| Cardiovascular | LDL, CRP, Lp(a), BNP | Risk factors -> Subclinical -> Events |
| Liver | ALT, AST, GGT, albumin, bilirubin | Steatosis -> NASH -> Fibrosis -> Cirrhosis |
| Thyroid | TSH, free T4, free T3, TPO antibodies | Subclinical -> Overt -> Complications |
| Iron | Ferritin, transferrin saturation, TIBC | Depletion -> Deficiency -> Anemia |
| Nutritional | B12, folate, vitamin D, zinc | Suboptimal -> Deficient -> Clinical |
| Kidney | eGFR, creatinine, BUN, cystatin C | Stage 1 -> Stage 5 (CKD) |
| Bone Health | Calcium, phosphorus, PTH, vitamin D | Normal -> Osteopenia -> Osteoporosis |
| Cognitive | Homocysteine, B12, folate, CRP | Risk factors -> Mild -> Significant |

### 7.2 Trajectory Output

Each trajectory includes:
- Current stage classification
- Years-to-progression estimate (if untreated)
- Intervention window description
- Potential risk reduction percentage
- Recommended monitoring frequency
- Actionable biomarker targets

---

## 8. Genotype Adjustment Engine

### 8.1 Modifier Genes (7)

| Gene | Variant | Affected Biomarkers | Adjustment |
|------|---------|-------------------|------------|
| MTHFR | C677T | Homocysteine, folate | Narrowed ranges, higher folate targets |
| APOE | e4 allele | LDL, total cholesterol | Tighter lipid targets |
| PNPLA3 | I148M | ALT, AST, GGT | Adjusted liver enzyme thresholds |
| HFE | C282Y, H63D | Ferritin, transferrin sat | Iron overload-specific ranges |
| GBA | Various | GCase activity | Carrier-specific enzyme ranges |
| HEXA | Various | Hex A activity | Carrier-specific enzyme ranges |
| FTO | rs9939609 | BMI, glucose, HbA1c | Obesity-adjusted metabolic ranges |

### 8.2 Adjustment Architecture

```
Standard Reference Range + Genotype Data
    |
    v
[Modifier Gene Lookup]
    |-- Check patient genotype against 7 modifier genes
    |-- Identify applicable adjustments
    |
    v
[Range Modification]
    |-- Apply narrowing/widening factors
    |-- Flag genotype-specific optimal targets
    |
    v
AdjustedRange(biomarker, standard_min, standard_max,
              adjusted_min, adjusted_max, genotype, rationale)
```

---

## 9. Critical Value Engine

### 9.1 Threshold Rules (21)

The engine monitors 21 critical threshold rules for life-threatening lab values:

| Biomarker | Critical Low | Critical High | Severity | Escalation |
|-----------|-------------|--------------|----------|------------|
| Potassium | < 2.5 mEq/L | > 6.5 mEq/L | Critical | Emergency |
| Sodium | < 120 mEq/L | > 160 mEq/L | Critical | Emergency |
| Glucose | < 40 mg/dL | > 500 mg/dL | Critical | Emergency |
| Calcium | < 6.0 mg/dL | > 13.0 mg/dL | Critical | Emergency |
| Hemoglobin | < 7.0 g/dL | > 20.0 g/dL | Critical | Emergency |
| Platelets | < 20 K/uL | > 1000 K/uL | Critical | Emergency |
| INR | -- | > 5.0 | Urgent | Physician |
| Troponin | -- | > 0.04 ng/mL | Critical | Emergency |
| Lactate | -- | > 4.0 mmol/L | Critical | Emergency |
| Creatinine | -- | > 10.0 mg/dL | Critical | Emergency |

(Plus 11 additional rules for WBC, pH, pCO2, pO2, bicarbonate, magnesium, phosphorus, bilirubin, ammonia, TSH, and lithium.)

---

## 10. Discordance Detector

### 10.1 Cross-Biomarker Patterns

Identifies physiologically contradictory or clinically unexpected biomarker combinations:

| Pattern | Biomarker A | Biomarker B | Implication |
|---------|-----------|-----------|------------|
| Thyroid discordance | Normal TSH | Elevated free T4 | Assay interference or central thyroid disease |
| HbA1c-glucose | Normal HbA1c | Elevated fasting glucose | Possible hemoglobin variant or early diabetes |
| Iron paradox | Low ferritin | Elevated transferrin sat | Mixed iron picture, consider inflammation |
| Liver-albumin | Elevated ALT/AST | Normal albumin | Acute vs chronic liver disease differentiation |
| Renal-calcium | Elevated creatinine | Elevated calcium | Evaluate for primary hyperparathyroidism |
| B12-MMA | Normal B12 | Elevated methylmalonic acid | Functional B12 deficiency |

---

## 11. RAG Pipeline

### 11.1 RAG Architecture

```
Query --> Embed (BGE-small + instruction prefix)
  --> Parallel Search (14 collections, top_k=5 each, score_threshold=0.40)
  --> Weighted Scoring (per-collection weights, sum=1.0)
  --> Deduplication (by ID and text content hash)
  --> Knowledge Augmentation (pharmacogenes, disease domains, interactions)
  --> LLM Synthesis (Claude Sonnet 4.6, system prompt, conversation memory)
  --> Streaming Response with Evidence Citations
```

### 11.2 Embedding Configuration

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
| Search mode | Asymmetric |
| Embedding cache | 256 entries, FIFO eviction |

### 11.3 Citation Scoring

| Score Range | Relevance Level |
|------------|----------------|
| >= 0.75 | High |
| >= 0.60 | Medium |
| < 0.60 | Low |

### 11.4 Knowledge Graph (6 Domains)

| Domain | Content | Purpose |
|--------|---------|---------|
| Pharmacogenomics | 14 pharmacogene profiles with star allele tables | PGx phenotyping context |
| Cardiovascular | Risk factors, biomarker thresholds, treatment targets | CV disease trajectory |
| Metabolic | Diabetes, thyroid, nutritional biomarker patterns | Metabolic disease detection |
| Hepatic/Renal | Liver function, kidney function, staging criteria | Organ function assessment |
| Hematology | CBC interpretation, iron studies, coagulation | Blood disorder detection |
| Aging | PhenoAge/GrimAge coefficients, aging biomarker norms | Biological age estimation |

---

## 12. Agent Orchestrator

### 12.1 Key Statistics

| Metric | Value |
|--------|-------|
| Source modules (src/) | 16 Python files |
| Total source lines | ~20,200 |
| Test suite | 709 tests across 18 files, all passing |
| Milvus collections | 14 (13 owned + 1 shared read-only) |
| Total seed vectors | 1,297 |
| Pharmacogenes mapped | 14 with CPIC-guided phenotyping |
| Disease trajectory categories | 9 |
| Genotype modifier genes | 7 |
| Critical value rules | 21 |
| Report sections | 12 |
| Export formats | 5 (Markdown, PDF, FHIR R4, CSV, JSON) |
| Supported languages | 7 (English, Spanish, Chinese, Hindi, French, Arabic, Portuguese) |

### 12.2 Error Handling Strategy

The orchestrator implements graceful degradation:

1. **Milvus unavailable**: Analysis engines still function (pure computation); RAG returns error
2. **LLM unavailable**: Returns analysis results without LLM synthesis
3. **PGx engine error**: Logs warning, returns available engine results
4. **Biological age error**: Logs warning, returns response without age calculation
5. **Critical value detection**: Always runs first; alerts bypass other engine failures
6. **Timeout**: Returns partial results with timeout indicator

### 12.3 API Endpoints Summary

**Core endpoints (api/main.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/healthz` | GET | Health check |
| `/readyz` | GET | Readiness check (Milvus connection) |
| `/metrics` | GET | Prometheus metrics endpoint |
| `/v1/query` | POST | RAG query with streaming support |

**Analysis endpoints (api/routes/analysis.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/analyze` | POST | Full patient analysis (all engines) |
| `/v1/biological-age` | POST | PhenoAge/GrimAge calculation |
| `/v1/disease-trajectory` | POST | Disease trajectory prediction |
| `/v1/pgx-map` | POST | Pharmacogenomic mapping |
| `/v1/genotype-adjust` | POST | Genotype-based range adjustment |

**Event endpoints (api/routes/events.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/events/cross-modal` | POST | Cross-modal event handling |
| `/v1/biomarker-alert` | POST | Critical value alert processing |

**Report endpoints (api/routes/reports.py):**

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/v1/report/generate` | POST | Generate clinical report |
| `/v1/report/{id}/pdf` | GET | Download report as PDF |
| `/v1/report/fhir` | POST | Generate FHIR R4 bundle |

### 12.4 Port Configuration

| Service | Port | Protocol |
|---------|------|----------|
| FastAPI REST API | 8529 | HTTP |
| Streamlit UI | 8533 | HTTP |
| Milvus (shared) | 19530 | gRPC |

---

## 13. Data Model Architecture

### 13.1 Export Pipeline (5 Formats)

| Format | Module | Description |
|--------|--------|------------|
| Markdown | `export.py` | 12-section clinical report with tables |
| PDF | `export.py` | ReportLab with NVIDIA-themed colors (#76B900) |
| FHIR R4 | `export.py` | Bundle with Patient, Observation, DiagnosticReport |
| CSV | `export.py` | Flat tabular biomarker results with adjustments |
| JSON | `export.py` | Full AnalysisResult via Pydantic model_dump() |

### 13.2 Translation Engine (7 Languages)

| Language | Code | Medical Terminology Coverage |
|----------|------|------------------------------|
| English | en | Full (primary) |
| Spanish | es | Clinical terms, biomarker names |
| Chinese | zh | Clinical terms, biomarker names |
| Hindi | hi | Clinical terms, biomarker names |
| French | fr | Clinical terms, biomarker names |
| Arabic | ar | Clinical terms, biomarker names |
| Portuguese | pt | Clinical terms, biomarker names |

### 13.3 File Structure

```
precision_biomarker_agent/
|-- api/
|   |-- __init__.py
|   |-- main.py                  (465 lines) Entry point, core endpoints
|   +-- routes/
|       |-- analysis.py          (495 lines) /v1/analyze, /v1/biological-age, etc.
|       |-- events.py            (326 lines) /v1/events/cross-modal
|       +-- reports.py           (296 lines) /v1/report/generate
|
|-- app/
|   |-- biomarker_ui.py          (1,863 lines) 8-tab Streamlit application
|   |-- patient_360.py           (670 lines) Cross-agent Patient 360 dashboard
|   +-- protein_viewer.py        (168 lines) 3D protein structure viewer
|
|-- config/
|   +-- settings.py              (139 lines) PrecisionBiomarkerSettings
|
|-- data/
|   +-- reference/               # 18 JSON seed files
|
|-- src/
|   |-- agent.py                 (610 lines) Plan-analyze-search-synthesize
|   |-- audit.py                 (83 lines) HIPAA-compliant audit logging
|   |-- biological_age.py        (408 lines) PhenoAge + GrimAge
|   |-- collections.py           (1,391 lines) Milvus collection management
|   |-- critical_values.py       (179 lines) Critical threshold detection
|   |-- discordance_detector.py  (299 lines) Cross-biomarker anomalies
|   |-- disease_trajectory.py    (1,421 lines) 9-category trajectories
|   |-- export.py                (1,392 lines) 5-format export
|   |-- genotype_adjustment.py   (1,225 lines) 7 modifier genes
|   |-- knowledge.py             (1,326 lines) 6 disease domains
|   |-- lab_range_interpreter.py (221 lines) Standard vs optimal
|   |-- models.py                (786 lines) 14 collection + 8 analysis models
|   |-- pharmacogenomics.py      (1,503 lines) 14 pharmacogenes
|   |-- rag_engine.py            (573 lines) Multi-collection RAG
|   |-- report_generator.py      (993 lines) 12-section reports
|   +-- translation.py           (217 lines) 7 languages
|
|-- scripts/
|   |-- setup_collections.py
|   |-- seed_all.py
|   +-- demo_validation.py
|
|-- tests/                       # 18 test files, 709 tests
|
|-- docker-compose.yml           # 6-service stack
|-- Dockerfile
|-- requirements.txt
+-- README.md
```

### 13.4 Security

- API key authentication via `BIOMARKER_API_KEY` environment variable
- HIPAA-compliant audit logging (`src/audit.py`)
- Milvus filter expression sanitization
- CORS middleware with explicit origin allowlist
- Non-root container execution (`biomarkeruser`)
- Request size limiting (default 10 MB)

---

!!! warning "Clinical Decision Support Disclaimer"
    The Precision Biomarker Agent is a clinical decision support research tool for biomarker analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
