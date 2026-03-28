# Pharmacogenomics Intelligence Agent -- Demo Guide

**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Starting the System](#2-starting-the-system)
3. [Verifying Health](#3-verifying-health)
4. [Demo Scenarios: 10-Tab Walkthrough](#4-demo-scenarios-10-tab-walkthrough)
5. [Sample API Queries](#5-sample-api-queries)
6. [Demo Talking Points](#6-demo-talking-points)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Prerequisites

### Hardware

- NVIDIA DGX Spark (recommended) or any system with 16GB+ RAM
- GPU optional (accelerates embedding generation)

### Software

- Docker and Docker Compose v2+
- Python 3.12+ (for manual setup without Docker)
- Anthropic API key with Claude Sonnet 4.6 access

### Environment Setup

```bash
cd ai_agent_adds/pharmacogenomics_intelligence_agent

# Copy environment template and add your Anthropic API key
cp .env.example .env
# Edit .env and set ANTHROPIC_API_KEY=sk-ant-...
```

---

## 2. Starting the System

### Option A: Docker Compose (Recommended)

```bash
# Start all 6 services (etcd, MinIO, Milvus, Streamlit, API, setup)
docker compose up -d

# Watch the setup process (creates collections, seeds 240 records)
docker compose logs -f pgx-setup

# Wait for "PGx Setup complete!" message, then verify:
docker compose ps
```

Expected output: 5 running services (pgx-setup exits after completion).

### Option B: Manual Setup

```bash
# 1. Ensure Milvus is running on localhost:19530

# 2. Install dependencies
pip install -r requirements.txt

# 3. Create collections and seed data
python scripts/setup_collections.py --drop-existing --seed
python scripts/seed_knowledge.py

# 4. Start the API server (port 8107)
uvicorn api.main:app --host 0.0.0.0 --port 8107 &

# 5. Start the Streamlit UI (port 8507)
streamlit run app/pgx_ui.py --server.port 8507
```

---

## 3. Verifying Health

### API Health Check

```bash
curl http://localhost:8107/health | python -m json.tool
```

Expected response:
```json
{
  "status": "healthy",
  "collections": 15,
  "total_vectors": 240
}
```

### Collection Counts

```bash
curl http://localhost:8107/collections | python -m json.tool
```

### Knowledge Graph Stats

```bash
curl http://localhost:8107/knowledge/stats | python -m json.tool
```

Expected: 25 pharmacogenes, 12 drug categories, 15 HLA-drug associations.

### UI Access

Open http://localhost:8507 in your browser. The sidebar should show all 15 collections with their vector counts.

---

## 4. Demo Scenarios: 10-Tab Walkthrough

Each scenario maps to one of the 10 Streamlit tabs. Follow these in order for a complete demonstration.

### Tab 1: Dashboard Overview

**Purpose:** Show the system's scope and capabilities at a glance.

**Walkthrough:**

1. Open the **Dashboard** tab (selected by default).
2. Point out the 15 collection stats in the sidebar:
   ```
   +------------------------------------+
   |  SIDEBAR                           |
   |  Connected Collections: 15/15      |
   |  Total Vectors: 240                |
   |                                    |
   |  pgx_gene_reference: 20 vectors   |
   |  pgx_drug_guidelines: 25 vectors  |
   |  pgx_drug_interactions: 18 vectors|
   |  pgx_hla_hypersensitivity: 12 vec |
   |  pgx_phenoconversion: 15 vectors  |
   |  ... (all 15 collections listed)  |
   +------------------------------------+
   ```
3. Show the knowledge graph summary in the main panel:
   ```
   +---------------------------------------------+
   |  PGx INTELLIGENCE AGENT DASHBOARD           |
   |                                              |
   |  Pharmacogenes:  25                          |
   |  Drugs Covered:  100+                        |
   |  Drug Categories: 12                         |
   |  HLA Associations: 15                        |
   |  Dosing Algorithms: 9                        |
   |  Clinical Workflows: 8                       |
   |  Prometheus Metrics: 22                      |
   |                                              |
   |  [Knowledge Graph Explorer]                  |
   |  CYP2D6 | CYP2C19 | CYP2C9 | DPYD | ...   |
   +---------------------------------------------+
   ```
4. Click on a pharmacogene (e.g., CYP2D6) to see its details: chromosome location, key variants, substrate count, CPIC guidelines available.

**Talking points:**
- "This system covers the 25 most clinically actionable pharmacogenes as defined by CPIC."
- "All 15 collections are searched in parallel for every query -- response time is under 5 seconds."
- "The knowledge graph contains 2,657 lines of structured pharmacogenomic data that augments every response."

### Tab 2: Single Drug Check

**Purpose:** Demonstrate PGx-guided prescribing for a single medication.

**Walkthrough:**

1. Switch to the **Drug Check** tab.
2. Enter the following values:
   - Drug: `codeine`
   - Gene: `CYP2D6`
   - Phenotype: `poor_metabolizer`
3. Click **Check Drug**.
4. Review the response:
   ```
   +---------------------------------------------+
   |  CLINICAL ALERT: CRITICAL                    |
   |  Severity: CONTRAINDICATED                   |
   |                                              |
   |  Codeine is a CYP2D6 prodrug. Poor metab-   |
   |  olizers cannot convert codeine to its       |
   |  active metabolite morphine.                 |
   |                                              |
   |  Recommendation: AVOID codeine.              |
   |  Alternatives:                               |
   |  - Morphine (non-CYP2D6 pathway)            |
   |  - Oxycodone (CYP3A4 primary pathway)       |
   |  - Non-opioid analgesics                    |
   |                                              |
   |  Evidence: CPIC Level A                      |
   |  Reference: Crews et al., CPT 2021          |
   +---------------------------------------------+
   ```

**Follow-up query:** Change phenotype to `ultra_rapid_metabolizer`:
- System warns about excessive morphine production and respiratory depression risk.
- Deaths reported in children and breastfeeding infants of UM mothers.
- Recommendation: Also AVOID codeine for ultra-rapid metabolizers.

**Second follow-up:** Try drug = `clopidogrel`, gene = `CYP2C19`, phenotype = `poor_metabolizer`:
- System flags FDA Boxed Warning about reduced antiplatelet effect.
- Recommends prasugrel or ticagrelor as alternatives.
- Notes elevated risk of stent thrombosis and cardiovascular events.

### Tab 3: Polypharmacy Review

**Purpose:** Show multi-drug interaction analysis with phenoconversion detection.

**Walkthrough:**

1. Switch to the **Medication Review** tab.
2. Enter medications (comma-separated): `codeine, fluoxetine, omeprazole, simvastatin`
3. Click **Review Medications**.
4. Review the interaction analysis:
   ```
   +---------------------------------------------+
   |  MEDICATION REVIEW RESULTS                   |
   |                                              |
   |  PHENOCONVERSION ALERT (CYP2D6):            |
   |  Fluoxetine is a STRONG CYP2D6 inhibitor.   |
   |  Any CYP2D6 Normal/Intermediate metabolizer |
   |  is effectively converted to Poor Metab.    |
   |                                              |
   |  DRUG-DRUG-GENE INTERACTION:                 |
   |  Codeine + Fluoxetine = codeine cannot be   |
   |  activated (CYP2D6 inhibited by fluoxetine) |
   |  Action: Avoid codeine while on fluoxetine  |
   |                                              |
   |  STATIN NOTE:                                |
   |  Simvastatin: Check SLCO1B1 rs4149056       |
   |  status for myopathy risk.                  |
   |                                              |
   |  OMEPRAZOLE NOTE:                            |
   |  CYP2C19 metabolizer status affects PPI     |
   |  efficacy and H. pylori eradication rates.  |
   +---------------------------------------------+
   ```

**Talking points:**
- "The system detected that fluoxetine phenoconverts any CYP2D6 metabolizer to a poor metabolizer."
- "This is a drug-drug-gene interaction -- something most drug interaction checkers miss entirely."
- "Even if this patient has a normal CYP2D6 genotype, the fluoxetine makes them functionally a poor metabolizer for codeine."

### Tab 4: Warfarin Dosing

**Purpose:** Demonstrate the IWPC genotype-guided warfarin dosing algorithm.

**Walkthrough:**

1. Switch to the **Warfarin Dosing** tab.
2. Enter patient parameters:
   - Age: `65`
   - Height: `170` cm
   - Weight: `80` kg
   - Race: `Caucasian`
   - CYP2C9: `*1/*3`
   - VKORC1: `A/G`
   - Amiodarone: `No`
   - Enzyme inducer: `No`
3. Click **Calculate Dose**.
4. Review the output:
   ```
   +---------------------------------------------+
   |  IWPC WARFARIN DOSE CALCULATION              |
   |                                              |
   |  Calculated Weekly Dose: 28.4 mg/week       |
   |  Daily Dose Equivalent: 4.1 mg/day          |
   |  Population Average: 35 mg/week             |
   |                                              |
   |  GENOTYPE IMPACT:                            |
   |  CYP2C9 *1/*3: Decreased function allele    |
   |  --> Reduced warfarin clearance              |
   |  --> Lower dose requirement (-19%)           |
   |                                              |
   |  VKORC1 A/G: Intermediate sensitivity        |
   |  --> Moderately increased vitamin K cycle    |
   |     sensitivity --> Lower dose requirement   |
   |                                              |
   |  Reference: Klein et al., NEJM 2009         |
   +---------------------------------------------+
   ```

**Variant demonstrations:**

Try CYP2C9 = `*3/*3`, VKORC1 = `A/A`:
- Dose drops to approximately 10-15 mg/week (extreme low-dose genotype)
- System highlights very high over-anticoagulation risk

Try CYP2C9 = `*1/*1`, VKORC1 = `G/G`, enzyme inducer = `Yes`:
- Dose increases to approximately 45-55 mg/week
- System explains enzyme induction increases warfarin clearance

### Tab 5: Chemotherapy Toxicity Screening

**Purpose:** Show pre-treatment DPYD screening for fluoropyrimidine chemotherapy.

**Walkthrough:**

1. Switch to the **Chemo Safety** tab.
2. Enter: Gene = `DPYD`, Diplotype = `*1/*2A`
3. Click **Screen**.
4. Review the output:
   ```
   +---------------------------------------------+
   |  CHEMOTHERAPY SAFETY SCREEN                  |
   |                                              |
   |  ALERT: MAJOR                                |
   |  DPYD*2A (IVS14+1G>A) is a no-function     |
   |  allele. Heterozygous carriers (*1/*2A)     |
   |  have ~50% reduced DPD enzyme activity.     |
   |                                              |
   |  Activity Score: 1.0                         |
   |  (1.0 normal + 0.0 no-function)             |
   |                                              |
   |  Recommendation: REDUCE 5-fluorouracil or   |
   |  capecitabine dose by 50%.                  |
   |                                              |
   |  Without dose reduction: 73% risk of        |
   |  grade 3+ toxicity (severe mucositis,       |
   |  myelosuppression, hand-foot syndrome).     |
   |                                              |
   |  With 50% dose reduction: 28% risk          |
   |  (comparable to wild-type patients).        |
   |                                              |
   |  Reference: Amstutz et al., CPT 2018        |
   |  Reference: Henricks et al., Lancet Oncol   |
   +---------------------------------------------+
   ```

**Follow-up:** Enter Diplotype = `*2A/*2A`:
- System flags as CONTRAINDICATED (activity score 0.0)
- Complete DPD deficiency; 10-20% mortality at standard doses
- Recommendation: Avoid fluoropyrimidines entirely; use alternative regimens

### Tab 6: HLA Screening

**Purpose:** Demonstrate pre-prescription HLA allele screening.

**Walkthrough:**

1. Switch to the **HLA Screening** tab.
2. Enter: Drug = `carbamazepine`, HLA typing = `HLA-B*15:02 positive`
3. Click **Screen**.
4. Review the output:
   ```
   +---------------------------------------------+
   |  HLA SCREENING RESULT                        |
   |                                              |
   |  Status: CONTRAINDICATED                     |
   |  Severity: FATAL/SEVERE                      |
   |                                              |
   |  HLA-B*15:02 positive patients have a        |
   |  10-30% risk of Stevens-Johnson syndrome    |
   |  (SJS) / toxic epidermal necrolysis (TEN)   |
   |  with carbamazepine exposure.               |
   |                                              |
   |  FDA Boxed Warning applies.                  |
   |                                              |
   |  ALTERNATIVES:                               |
   |  - Lamotrigine (with slow titration)        |
   |  - Levetiracetam                            |
   |  - Valproic acid                            |
   |                                              |
   |  POPULATION NOTE:                            |
   |  HLA-B*15:02 prevalence:                    |
   |  Southeast Asian: 2-15%                     |
   |  East Asian: 2-6%                           |
   |  European: < 1%                             |
   |  African: < 1%                              |
   +---------------------------------------------+
   ```

**Follow-up screens:**

Drug = `abacavir`, HLA = `HLA-B*57:01 positive`:
- Status: CONTRAINDICATED
- FDA-mandated testing; hypersensitivity syndrome in 5-8% of carriers
- No abacavir prescribing for HLA-B*57:01 positive patients

Drug = `allopurinol`, HLA = `HLA-B*58:01 positive`:
- Status: CONTRAINDICATED
- SJS/TEN risk; more common in Southeast Asian and African populations
- Alternative: febuxostat

Drug = `carbamazepine`, HLA = `HLA-B*15:02 negative`:
- Status: SAFE
- Standard prescribing with routine monitoring

### Tab 7: Report Generation

**Purpose:** Show clinical report export capabilities.

**Walkthrough:**

1. Switch to the **Report Generator** tab.
2. Enter query: `Complete PGx profile for CYP2D6 *1/*4, CYP2C19 *1/*2 patient on codeine, clopidogrel, and omeprazole`
3. Click **Generate Report**.
4. Wait for the RAG pipeline to complete (evidence retrieval + LLM synthesis).
5. Review the generated report, which includes:
   - Patient pharmacogenomic profile summary
   - Drug-by-drug PGx analysis with alerts
   - Dosing recommendations
   - Evidence citations
6. Export in each format:
   ```
   +---------------------------------------------+
   |  EXPORT OPTIONS                               |
   |                                              |
   |  [Download Markdown]  Human-readable report  |
   |  [Download JSON]      Machine-readable data  |
   |  [Download PDF]       Styled PGx Passport    |
   |  [Download FHIR R4]   DiagnosticReport       |
   |                       Bundle (LOINC 69548-6) |
   +---------------------------------------------+
   ```

**Export format details:**

- **Markdown**: Alert table with severity color coding, drug interaction matrix, evidence citations with hyperlinks
- **JSON**: Pydantic-serialized response with structured alert objects, evidence scores, and knowledge graph context
- **PDF**: Styled report with header, patient identifier placeholder, PGx Passport format, alert summary table
- **FHIR R4**: DiagnosticReport Bundle containing PGx Observations coded with LOINC 69548-6 (Pharmacogenomic Analysis Report), ready for EHR ingestion

### Tab 8: Evidence Explorer

**Purpose:** Browse raw evidence across all 15 collections.

**Walkthrough:**

1. Switch to the **Evidence Explorer** tab.
2. Enter search: `CYP2D6 codeine metabolism`
3. Click **Search**.
4. Review evidence grouped by collection:
   ```
   +---------------------------------------------+
   |  EVIDENCE RESULTS (24 hits across 8         |
   |  collections)                                |
   |                                              |
   |  pgx_drug_guidelines (5 hits)               |
   |  [0.89] CPIC Level A: CYP2D6 UM - avoid   |
   |         codeine; risk of morphine toxicity  |
   |  [0.87] CPIC Level A: CYP2D6 PM - avoid   |
   |         codeine; no analgesic effect        |
   |  [0.82] CPIC Level A: CYP2D6 NM - use     |
   |         per standard dosing                 |
   |  ...                                        |
   |                                              |
   |  pgx_drug_interactions (4 hits)             |
   |  [0.85] CYP2D6 activates codeine to        |
   |         morphine; PK: metabolism             |
   |  ...                                        |
   |                                              |
   |  pgx_gene_reference (3 hits)                |
   |  [0.78] CYP2D6*4 - no function -           |
   |         rs3892097 (1846G>A, splice defect)  |
   |  ...                                        |
   +---------------------------------------------+
   ```
5. Use collection-specific filtering in the sidebar to narrow results.
6. Click citation links to navigate to PubMed/ClinicalTrials.gov sources.

**Talking points:**
- "Every piece of evidence has a relevance score based on cosine similarity weighted by collection importance."
- "High relevance (>= 0.75) citations are shown in green, medium (>= 0.60) in yellow, and low (< 0.60) in grey."
- "You can filter by collection to see only drug guidelines, only clinical trials, or any combination."

### Tab 9: Phenoconversion Modeler

**Purpose:** Interactive visualization of how concomitant medications alter metabolizer phenotype.

**Walkthrough:**

1. Switch to the **Phenoconversion Modeler** tab.
2. Set genetic phenotype: CYP2D6 = `Normal Metabolizer`
3. Add medications one at a time and observe the phenotype shift:

**Step 1: No medications**
```
Genetic Phenotype:   CYP2D6 Normal Metabolizer (AS = 2.0)
Effective Phenotype: CYP2D6 Normal Metabolizer (AS = 2.0)
Shift: None
```

**Step 2: Add fluoxetine**
```
Genetic Phenotype:   CYP2D6 Normal Metabolizer (AS = 2.0)
Medications:         fluoxetine (STRONG CYP2D6 inhibitor)
Effective Phenotype: CYP2D6 Poor Metabolizer
Shift: NM --> PM (strong inhibition)
```

**Step 3: Add bupropion**
```
Genetic Phenotype:   CYP2D6 Normal Metabolizer (AS = 2.0)
Medications:         fluoxetine (STRONG), bupropion (STRONG)
Effective Phenotype: CYP2D6 Poor Metabolizer
Shift: Still PM (additional strong inhibitor does not worsen PM status)
Note: Multiple strong inhibitors confirm the phenoconversion
```

**Step 4: Remove fluoxetine, keep bupropion**
```
Genetic Phenotype:   CYP2D6 Normal Metabolizer (AS = 2.0)
Medications:         bupropion (STRONG CYP2D6 inhibitor)
Effective Phenotype: CYP2D6 Poor Metabolizer
Shift: NM --> PM (bupropion alone is sufficient)
```

**Step 5: Change to moderate inhibitor**
```
Genetic Phenotype:   CYP2D6 Normal Metabolizer (AS = 2.0)
Medications:         duloxetine (MODERATE CYP2D6 inhibitor)
Effective Phenotype: CYP2D6 Intermediate Metabolizer
Shift: NM --> IM (moderate inhibition)
```

**Multi-enzyme demonstration:** Show that a single medication can affect multiple CYP enzymes. Enter fluvoxamine:
- CYP2C19: Strong inhibitor (NM --> PM)
- CYP1A2: Strong inhibitor (NM --> PM)
- CYP2C9: Moderate inhibitor (NM --> IM)
- Three enzymes affected by one drug

### Tab 10: Population Analytics

**Purpose:** Show population-specific allele frequency comparisons.

**Walkthrough:**

1. Switch to the **Population Analytics** tab.
2. Select gene: `CYP2D6`
3. View the population comparison:
   ```
   +---------------------------------------------+
   |  CYP2D6 ALLELE FREQUENCIES BY POPULATION    |
   |                                              |
   |  Allele    EUR    EAS    AFR    SAS    AMR  |
   |  *1 (NF)  35%    25%    35%    40%    35%  |
   |  *2 (NF)  25%    10%    15%    20%    20%  |
   |  *4 (NoF) 22%    1%     3%     8%     10%  |
   |  *5 (NoF)  3%    6%     5%     3%     4%  |
   |  *10 (DF)  2%    40%    4%     5%     5%  |
   |  *17 (DF) <1%   <1%    25%    <1%    <1%  |
   |  *41 (DF)  8%    3%    10%    15%    10%  |
   |                                              |
   |  NF=Normal Function, NoF=No Function,       |
   |  DF=Decreased Function                      |
   |                                              |
   |  METABOLIZER PHENOTYPE DISTRIBUTION:         |
   |  UM:  EUR 1-2% | EAS <1%  | AFR 2-5%       |
   |  NM:  EUR 70%  | EAS 55%  | AFR 65%        |
   |  IM:  EUR 12%  | EAS 35%  | AFR 25%        |
   |  PM:  EUR 7%   | EAS 1%   | AFR 3%         |
   +---------------------------------------------+
   ```

4. Highlight key observations:
   - CYP2D6*10 is the most common decreased-function allele in East Asian populations (~40%) but rare in Europeans (~2%).
   - CYP2D6*4 is the most common no-function allele in Europeans (~22%) but rare in East Asians (~1%).
   - CYP2D6*17 is common in African populations (~25%) but rare elsewhere.

5. Switch to gene `CYP2C19` and show:
   - *2 (no function): 25-35% in East Asian vs 12-15% in European
   - This explains why CYP2C19 poor metabolizer prevalence is 15-20% in East Asian populations vs 2-5% in European populations.

**Talking points:**
- "PGx testing panels developed primarily for European populations may miss clinically important alleles in other populations."
- "CYP2D6*17 was not characterized until long after *4 and *5, despite being the most common decreased-function allele in African populations."
- "Health equity requires population-aware PGx testing that includes alleles relevant to all ancestries."

---

## 5. Sample API Queries

### 5.1 Health and Status Endpoints

```bash
# Health check
curl http://localhost:8107/health | python -m json.tool

# Collection inventory
curl http://localhost:8107/collections | python -m json.tool

# Knowledge graph statistics
curl http://localhost:8107/knowledge/stats | python -m json.tool

# Prometheus metrics
curl http://localhost:8107/metrics
```

### 5.2 Full RAG Query

```bash
curl -X POST http://localhost:8107/query \
  -H "Content-Type: application/json" \
  -d '{"question": "What are the CYP2D6 implications for codeine prescribing?"}'
```

### 5.3 Evidence-Only Search (No LLM)

```bash
curl -X POST http://localhost:8107/search \
  -H "Content-Type: application/json" \
  -d '{"question": "warfarin CYP2C9 VKORC1 dosing"}'
```

### 5.4 Cross-Collection Entity Search

```bash
curl -X POST http://localhost:8107/find-related \
  -H "Content-Type: application/json" \
  -d '{"entity": "CYP2D6"}'
```

### 5.5 Single Drug PGx Check

```bash
# Codeine + CYP2D6 Poor Metabolizer
curl -X POST http://localhost:8107/v1/pgx/drug-check \
  -H "Content-Type: application/json" \
  -d '{"drug": "codeine", "gene": "CYP2D6", "phenotype": "poor_metabolizer"}'

# Clopidogrel + CYP2C19 Poor Metabolizer
curl -X POST http://localhost:8107/v1/pgx/drug-check \
  -H "Content-Type: application/json" \
  -d '{"drug": "clopidogrel", "gene": "CYP2C19", "phenotype": "poor_metabolizer"}'

# Simvastatin + SLCO1B1 Decreased Function
curl -X POST http://localhost:8107/v1/pgx/drug-check \
  -H "Content-Type: application/json" \
  -d '{"drug": "simvastatin", "gene": "SLCO1B1", "phenotype": "decreased_function"}'

# Tamoxifen + CYP2D6 Intermediate Metabolizer
curl -X POST http://localhost:8107/v1/pgx/drug-check \
  -H "Content-Type: application/json" \
  -d '{"drug": "tamoxifen", "gene": "CYP2D6", "phenotype": "intermediate_metabolizer"}'
```

### 5.6 Polypharmacy Medication Review

```bash
curl -X POST http://localhost:8107/v1/pgx/medication-review \
  -H "Content-Type: application/json" \
  -d '{"medications": ["codeine", "fluoxetine", "omeprazole", "simvastatin"]}'
```

### 5.7 IWPC Warfarin Dosing

```bash
# Standard patient
curl -X POST http://localhost:8107/v1/pgx/dosing/warfarin \
  -H "Content-Type: application/json" \
  -d '{
    "age": 65,
    "height_cm": 170,
    "weight_kg": 80,
    "race": "caucasian",
    "cyp2c9_genotype": "*1/*3",
    "vkorc1_genotype": "A/G",
    "amiodarone": false,
    "enzyme_inducer": false
  }'

# Extreme low-dose genotype
curl -X POST http://localhost:8107/v1/pgx/dosing/warfarin \
  -H "Content-Type: application/json" \
  -d '{
    "age": 75,
    "height_cm": 160,
    "weight_kg": 65,
    "race": "asian",
    "cyp2c9_genotype": "*3/*3",
    "vkorc1_genotype": "A/A",
    "amiodarone": true,
    "enzyme_inducer": false
  }'

# High-dose genotype with enzyme inducer
curl -X POST http://localhost:8107/v1/pgx/dosing/warfarin \
  -H "Content-Type: application/json" \
  -d '{
    "age": 45,
    "height_cm": 185,
    "weight_kg": 95,
    "race": "black",
    "cyp2c9_genotype": "*1/*1",
    "vkorc1_genotype": "G/G",
    "amiodarone": false,
    "enzyme_inducer": true
  }'
```

### 5.8 HLA Screening

```bash
# Abacavir + HLA-B*57:01
curl -X POST http://localhost:8107/v1/pgx/hla-screen \
  -H "Content-Type: application/json" \
  -d '{"drug": "abacavir", "hla_alleles": ["HLA-B*57:01"]}'

# Carbamazepine + HLA-B*15:02
curl -X POST http://localhost:8107/v1/pgx/hla-screen \
  -H "Content-Type: application/json" \
  -d '{"drug": "carbamazepine", "hla_alleles": ["HLA-B*15:02"]}'

# Allopurinol + HLA-B*58:01
curl -X POST http://localhost:8107/v1/pgx/hla-screen \
  -H "Content-Type: application/json" \
  -d '{"drug": "allopurinol", "hla_alleles": ["HLA-B*58:01"]}'

# Carbamazepine with NEGATIVE HLA (safe result)
curl -X POST http://localhost:8107/v1/pgx/hla-screen \
  -H "Content-Type: application/json" \
  -d '{"drug": "carbamazepine", "hla_alleles": []}'
```

### 5.9 Phenoconversion Analysis

```bash
curl -X POST http://localhost:8107/v1/pgx/phenoconversion \
  -H "Content-Type: application/json" \
  -d '{
    "baseline_phenotype": "normal_metabolizer",
    "enzyme": "CYP2D6",
    "concomitant_drugs": ["fluoxetine", "omeprazole"]
  }'
```

### 5.10 Gene Reference Lookup

```bash
# List all pharmacogenes
curl http://localhost:8107/v1/pgx/genes | python -m json.tool
```

### 5.11 Drug Guideline Lookup

```bash
# List all drugs with PGx guidelines
curl http://localhost:8107/v1/pgx/drugs | python -m json.tool
```

---

## 6. Demo Talking Points

### Opening Statement

"The Pharmacogenomics Intelligence Agent translates genetic data into actionable drug prescribing recommendations. It integrates 15 specialized knowledge collections, 9 validated dosing algorithms, and a 25-pharmacogene knowledge graph into a single conversational interface -- all running on a single NVIDIA DGX Spark."

### Key Messages by Audience

**For clinicians:**
- "This system checks for drug-drug-gene interactions that standard drug interaction checkers miss."
- "Every recommendation traces to a CPIC guideline, FDA label, or peer-reviewed publication."
- "The warfarin dosing calculator implements the exact IWPC algorithm from Klein et al., NEJM 2009."
- "HLA screening prevents life-threatening hypersensitivity reactions like Stevens-Johnson syndrome."

**For pharmacists:**
- "Phenoconversion modeling shows how concomitant medications change effective metabolizer status."
- "The system covers 30+ CYP inhibitors and inducers classified by FDA potency."
- "Multi-drug reviews identify interactions that emerge only when the full medication list is analyzed together."

**For IT/technical audiences:**
- "The multi-collection RAG architecture uses 15 Milvus vector collections with weighted search."
- "Parallel search across all 15 collections completes in under 1 second."
- "The system includes 1,001 tests covering all clinical logic, with 100% pass rate."
- "FHIR R4 export supports EHR integration via DiagnosticReport Bundles."

**For executives/administrators:**
- "Adverse drug reactions cost $528 billion per year globally. PGx-guided prescribing reduces ADRs by 30% (PREPARE study, Lancet 2023)."
- "The system runs on a single workstation -- no cloud infrastructure or per-patient licensing fees."
- "Open-source with 23,049 lines of Python, zero proprietary dependencies."

### Differentiators vs. Competitors

- "Unlike OneOme or GeneSight, our algorithms are open-source and transparent."
- "Unlike PharmCAT, we provide a conversational interface with LLM-powered synthesis."
- "We are the only system that integrates phenoconversion modeling, HLA screening, dosing algorithms, and RAG-based evidence retrieval in a single platform."

---

## 7. Troubleshooting

### Milvus Connection Failed

```
Error: Cannot connect to Milvus at localhost:19530
```

**Fix:** Ensure Milvus is running and healthy:
```bash
docker compose ps milvus-standalone
docker compose logs milvus-standalone | tail -20

# If Milvus is not running:
docker compose up -d milvus-standalone

# Wait for health check to pass:
curl http://localhost:9091/healthz
```

### No Collections Found

```
Error: 0 collections connected
```

**Fix:** Run the setup script:
```bash
docker compose run --rm pgx-setup
# or manually:
python scripts/setup_collections.py --drop-existing --seed
```

### LLM API Errors

```
Error: Anthropic API key not found
```

**Fix:** Ensure `ANTHROPIC_API_KEY` is set:
```bash
echo $ANTHROPIC_API_KEY  # Should show sk-ant-...
# or add to .env file
```

```
Error: Anthropic API rate limit exceeded
```

**Fix:** Wait 60 seconds and retry. For sustained use, check your Anthropic API usage tier and consider upgrading. Reduce concurrent queries if running multiple demo sessions.

### Streamlit Not Loading

```
Error: ModuleNotFoundError
```

**Fix:** Ensure project root is on PYTHONPATH:
```bash
export PYTHONPATH=/path/to/pharmacogenomics_intelligence_agent:$PYTHONPATH
streamlit run app/pgx_ui.py --server.port 8507
```

### Slow Query Response

If queries take longer than 5 seconds:
1. Check Milvus memory usage: `docker stats pgx-milvus-standalone`
2. Verify index is loaded: Check `/health` endpoint for collection status
3. Reduce `TOP_K_PER_COLLECTION` from 5 to 3 for faster results
4. Check LLM API latency in Prometheus metrics: `curl http://localhost:8107/metrics | grep llm_api_latency`

### Empty or Low-Quality Responses

If the LLM response seems thin or generic:
1. Verify collections have data: `curl http://localhost:8107/collections | python -m json.tool`
2. Check that seed data was loaded: total vectors should be 240
3. Try a more specific query (include gene name, drug name, and phenotype)
4. Check the Evidence Explorer tab to see what evidence is being retrieved

### Port Conflicts

| Port | Service | Resolution |
|------|---------|-----------|
| 8507 | Streamlit UI | Set `PGX_STREAMLIT_PORT` to alternate port |
| 8107 | FastAPI API | Set `PGX_API_PORT` to alternate port |
| 19530 | Milvus | Set `PGX_MILVUS_PORT` to alternate port |

### Docker Resource Issues

```
Error: Container killed by OOM
```

**Fix:** Increase Docker memory limit:
```bash
# Check current Docker resource limits
docker info | grep -i memory

# For Docker Desktop, increase memory in Settings > Resources
# For Linux, check /etc/docker/daemon.json or systemd configuration
```

### Demo Reset

To reset the system to a clean demo state:

```bash
# Stop everything and remove data
docker compose down -v

# Restart fresh
docker compose up -d

# Wait for setup
docker compose logs -f pgx-setup
```

This will recreate all 15 collections and re-seed all 240 records. Total reset time: approximately 2 minutes.
