# Precision Biomarker Intelligence Agent -- Demo Guide

**Author:** Adam Jones
**Date:** March 2026
**Version:** 1.0.0

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Demo Flow Overview](#3-demo-flow-overview)
4. [Tab-by-Tab Walkthrough](#4-tab-by-tab-walkthrough)
5. [API Demo](#5-api-demo)
6. [Talking Points for Each Feature](#6-talking-points-for-each-feature)
7. [Troubleshooting Common Demo Issues](#7-troubleshooting-common-demo-issues)

---

## 1. Overview

The Precision Biomarker Intelligence Agent transforms raw clinical biomarker data and
genomic variants into actionable clinical intelligence. It combines laboratory result
interpretation, pharmacogenomic profiling, biological age estimation, disease risk
trajectory modeling, and RAG-powered evidence retrieval into a single unified platform.

**What the agent does:**

- Ingests patient demographics, lab panels, medications, and genomic data (VCF)
- Applies genotype-aware reference range adjustments and critical value detection
- Calculates biological age using PhenoAge and GrimAge surrogate algorithms
- Predicts disease trajectories across 9 clinical categories
- Maps pharmacogenomic star alleles to metabolizer phenotypes with CPIC-guided dosing
- Searches 14 Milvus vector collections via RAG for evidence-backed interpretation
- Generates 12-section clinical reports with PDF export and FHIR R4 bundle output
- Provides longitudinal trend tracking and rate-of-change analysis

**Demo duration:** ~20 minutes for the full walkthrough. Individual tabs can be
demonstrated independently in 2--3 minutes each.

**Architecture at a glance:**

```
Streamlit UI (:8528)  -->  FastAPI Backend (:8529)  -->  Milvus Vector DB
                                  |
                      +-----------+-----------+
                      |           |           |
                  Biomarker   PGx Engine   RAG Engine
                  Analysis    (14 genes)   (14 collections)
                      |           |           |
                  Claude LLM  CPIC DB     BGE-small-en-v1.5
                                          (384-dim embeddings)
```

---

## 2. Prerequisites

### 2.1 Services That Must Be Running

| Service            | Port  | Health Check                          |
|--------------------|-------|---------------------------------------|
| Streamlit UI       | 8528  | `curl http://localhost:8528`          |
| FastAPI Backend    | 8529  | `curl http://localhost:8529/health`   |
| Milvus             | 19530 | `curl http://localhost:9091/healthz`  |
| etcd (Milvus dep)  | 2379  | Checked via Milvus health             |
| MinIO (Milvus dep) | 9000  | `curl http://localhost:9000/minio/health/live` |

### 2.2 Environment Variables

Confirm the following are set before launching:

```bash
# Required
export ANTHROPIC_API_KEY="sk-ant-..."

# Verify services
curl -s http://localhost:8529/health | python3 -m json.tool
```

### 2.3 Quick Start

```bash
# From the agent directory
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent

# Start backend
python -m uvicorn api.main:app --host 0.0.0.0 --port 8529 &

# Start UI
streamlit run app/biomarker_ui.py --server.port 8528 &
```

Alternatively, use Docker Compose:

```bash
docker compose up -d
```

### 2.4 Demo Data

Two demo patients are pre-loaded. No additional data ingestion is required.

---

## 3. Demo Flow Overview

The recommended demo path walks through all 8 tabs in sequence. Each tab builds on
context from the previous one, creating a narrative arc from raw biomarkers to
clinical decision support.

```
Tab 1: Biomarker Analysis     -- "What do the labs show?"
  |
Tab 2: Biological Age         -- "How old is this patient biologically?"
  |
Tab 3: Disease Risk            -- "Where is this patient heading?"
  |
Tab 4: PGx Profile             -- "How should we dose medications?"
  |
Tab 5: Evidence Explorer       -- "What does the literature say?"
  |
Tab 6: Reports                 -- "Package it for the clinician."
  |
Tab 7: Patient 360             -- "See everything in one place."
  |
Tab 8: Longitudinal            -- "How are things changing over time?"
```

**Narrative arc:** Start with the raw data (Tab 1), derive deeper insights (Tabs 2--4),
validate with evidence (Tab 5), package for clinical use (Tab 6), then show the
unified and temporal views (Tabs 7--8).

---

## 4. Tab-by-Tab Walkthrough

### Demo Patients Reference

| Field           | Patient 1                         | Patient 2                          |
|-----------------|-----------------------------------|------------------------------------|
| **ID**          | HCLS-BIO-2026-00001               | HCLS-BIO-2026-00002               |
| **Age/Sex**     | 45M                               | 38F                                |
| **Ethnicity**   | Ashkenazi Jewish                  | Ashkenazi Jewish                   |
| **Genome**      | HG002 (NA24385)                   | --                                 |
| **BMI**         | 23.7                              | 22.6                               |
| **Medications** | Atorvastatin 10mg, Metformin 500mg | OCP (ethinyl estradiol/norgestimate) |

---

### Tab 1: Biomarker Analysis

**Purpose:** Load a patient and review the full biomarker panel with intelligent
interpretation.

**Steps:**

1. **Select demo patient.** From the patient selector dropdown, choose
   `HCLS-BIO-2026-00001` (45M, Ashkenazi Jewish, HG002).

2. **Review the biomarker dashboard.** The system displays 67+ biomarker results
   organized by clinical category.

3. **Highlight genotype adjustments.** Point out biomarkers where the system has
   adjusted reference ranges based on the patient's Ashkenazi Jewish ancestry and
   known genomic variants. These appear with a gene icon indicator.

   > *Talking point:* "Traditional lab reports use population-wide reference ranges.
   > Our agent adjusts ranges based on the patient's actual genotype -- for example,
   > creatinine reference ranges shift when we know the patient's ancestry and body
   > composition."

4. **Show critical values.** Any biomarker flagged as critically high or low appears
   with a red alert banner. Click on a critical value to see the recommended clinical
   action.

   > *Talking point:* "Critical values trigger immediate alerts. The system doesn't
   > just flag them -- it provides context on why the value is critical for THIS
   > specific patient."

5. **Demonstrate discordance alerts.** Show cases where biomarker values conflict
   with each other or with genomic data. The discordance detector cross-references
   related biomarkers to catch inconsistencies.

   > *Talking point:* "Discordance detection is something a seasoned clinician does
   > intuitively -- checking whether related labs tell a consistent story. We automate
   > that pattern recognition across all 67+ biomarkers simultaneously."

6. **Optionally switch to Patient 2** (HCLS-BIO-2026-00002, 38F) to show how the
   same panel adjusts for sex-specific reference ranges and OCP medication effects.

---

### Tab 2: Biological Age

**Purpose:** Demonstrate the gap between chronological and biological age using
validated algorithms.

**Steps:**

1. **Review PhenoAge calculation.** The system uses 9 clinical biomarkers to compute
   PhenoAge:
   - Albumin
   - Creatinine
   - Glucose
   - C-Reactive Protein (CRP)
   - Lymphocyte Percent
   - Mean Cell Volume (MCV)
   - Red Cell Distribution Width (RDW)
   - Alkaline Phosphatase
   - White Blood Cell Count

2. **Show the GrimAge surrogate.** A secondary biological age estimate using surrogate
   markers when DNA methylation data is unavailable.

3. **Explain age acceleration/deceleration.** The delta between chronological age (45)
   and biological age is displayed as either acceleration (biological > chronological,
   higher risk) or deceleration (biological < chronological, protective).

   > *Talking point:* "PhenoAge was developed from NHANES III data and validated
   > against mortality outcomes. A patient who is chronologically 45 but biologically
   > 52 has the disease risk profile of a 52-year-old. This single number synthesizes
   > 9 biomarkers into an actionable aging metric."

4. **Point out which biomarkers are driving the result.** The breakdown shows each
   biomarker's contribution to the biological age estimate, highlighting which ones
   are aging the patient faster or slower.

---

### Tab 3: Disease Risk

**Purpose:** Show predicted disease trajectories based on current biomarker patterns.

**Steps:**

1. **Review the 9 disease trajectory categories:**

   | #  | Category                 | Key Biomarkers                          |
   |----|--------------------------|-----------------------------------------|
   | 1  | Cardiovascular           | Lipid panel, CRP, homocysteine          |
   | 2  | Metabolic / Diabetes     | Glucose, HbA1c, insulin, HOMA-IR        |
   | 3  | Liver                    | ALT, AST, bilirubin, albumin, GGT       |
   | 4  | Kidney                   | Creatinine, BUN, eGFR, cystatin C       |
   | 5  | Thyroid                  | TSH, free T4, free T3                   |
   | 6  | Iron Metabolism          | Ferritin, transferrin sat, serum iron    |
   | 7  | Nutritional              | Omega-3, Vit D, B12, folate, Mg, Zn     |
   | 8  | Cognitive                | Homocysteine, B12, folate, hs-CRP       |
   | 9  | Bone Health              | Calcium, PTH, Vit D, CTX, P1NP         |

2. **Walk through the cardiovascular risk panel for Patient 1.** Note that
   Atorvastatin 10mg is already prescribed -- the system factors current medications
   into risk projections.

3. **Show the metabolic/diabetes trajectory.** Patient 1 is on Metformin 500mg.
   The system models expected trajectory with and without continued treatment.

4. **Demonstrate sex-specific differences.** Switch to Patient 2 to show how the
   risk categories change for a 38F on OCP, including estrogen-related
   cardiovascular risk modifiers.

   > *Talking point:* "We don't just flag abnormal values -- we project trajectories.
   > A fasting glucose of 105 means something very different in a 45-year-old already
   > on Metformin versus a 38-year-old with no metabolic history."

---

### Tab 4: PGx Profile

**Purpose:** Demonstrate pharmacogenomic profiling from genomic data to drug dosing
recommendations.

**Steps:**

1. **Review the 14 pharmacogenes:** Show each gene with its detected star alleles,
   diplotype, and resulting metabolizer status.

   | Gene     | Example Phenotype      | Clinical Impact                        |
   |----------|------------------------|----------------------------------------|
   | CYP2D6   | Intermediate Metabolizer | Codeine, tamoxifen, SSRIs             |
   | CYP2C19  | Poor Metabolizer       | Clopidogrel, PPIs, voriconazole        |
   | CYP2C9   | Normal Metabolizer     | Warfarin, NSAIDs, phenytoin            |
   | CYP3A5   | Non-Expressor          | Tacrolimus, many statins               |
   | SLCO1B1  | Decreased Function     | Statin myopathy risk                   |
   | VKORC1   | High Sensitivity       | Warfarin dose reduction                |
   | MTHFR    | Reduced Activity       | Folate metabolism, methotrexate        |
   | TPMT     | Normal Metabolizer     | Thiopurines (azathioprine)             |
   | DPYD     | Normal Metabolizer     | Fluoropyrimidines (5-FU, capecitabine) |

2. **Highlight the star allele to phenotype mapping.** Click on CYP2D6 to show
   the full star allele detail (e.g., *1/*4 -> Intermediate Metabolizer).

3. **Show drug interaction with current medications.** For Patient 1 on Atorvastatin:
   - Check SLCO1B1 status (statin transport gene)
   - If decreased function is detected, the system recommends dose adjustment or
     alternative statin selection per CPIC guidelines

4. **Demonstrate CPIC guidance integration.** Each gene-drug pair links to the
   relevant CPIC guideline with strength of recommendation and evidence level.

   > *Talking point:* "Ashkenazi Jewish populations have elevated carrier frequencies
   > for several pharmacogenomic variants. Patient 1's HG002 genome gives us real
   > variant calls -- not imputed data -- so the star allele assignments are
   > high-confidence."

5. **Show the dosing recommendation panel.** For each active medication, the system
   provides a traffic-light recommendation: green (standard dose), yellow (consider
   adjustment), red (contraindicated or major dose change needed).

---

### Tab 5: Evidence Explorer

**Purpose:** Demonstrate RAG-powered search across the biomarker knowledge base.

**Steps:**

1. **Show the 14 collections available for search:**
   - 13 biomarker-specific collections (one per clinical category, plus general
     clinical chemistry)
   - 1 genomic collection (ClinVar + AlphaMissense annotations)

2. **Run a sample query.** Type a clinical question such as:
   > "What is the significance of elevated homocysteine in an Ashkenazi Jewish male
   > with MTHFR reduced activity?"

3. **Review the results.** The RAG engine returns:
   - Ranked passages with citation scores
   - Source collection and document metadata
   - Relevance confidence indicator

4. **Show citation scoring.** Each retrieved passage has a relevance score. Explain
   that the system uses BGE-small-en-v1.5 embeddings (384-dimensional vectors) to
   compute semantic similarity.

5. **Demonstrate comparative analysis.** Select two biomarkers and run a comparative
   query to show how the system synthesizes evidence across multiple collections.

   > *Talking point:* "The RAG engine searches across 14 vector collections
   > simultaneously. Each collection is embedded with BGE-small-en-v1.5, a compact
   > but high-quality biomedical embedding model. The 384-dimensional vectors give
   > us fast retrieval without sacrificing semantic precision."

6. **Show how evidence links back to patient context.** Click on a retrieved passage
   to see how the system connects the evidence to the specific patient's biomarker
   values and genotype.

---

### Tab 6: Reports

**Purpose:** Generate and export a comprehensive clinical report.

**Steps:**

1. **Trigger report generation.** Click "Generate Report" for Patient 1. The system
   compiles a 12-section clinical report:

   | #  | Section                          |
   |----|----------------------------------|
   | 1  | Patient Demographics             |
   | 2  | Executive Summary                |
   | 3  | Biomarker Analysis               |
   | 4  | Critical Values & Alerts         |
   | 5  | Genotype Adjustments             |
   | 6  | Biological Age Assessment        |
   | 7  | Disease Risk Trajectories        |
   | 8  | Pharmacogenomic Profile          |
   | 9  | Drug Interaction Analysis        |
   | 10 | Evidence Summary                 |
   | 11 | Recommendations                  |
   | 12 | Appendix (Raw Data & Methods)    |

2. **Walk through the executive summary.** This section distills the entire analysis
   into 3--5 key findings with clinical priority rankings.

3. **Show PDF export.** Click "Export PDF" to generate a formatted clinical document
   suitable for EHR attachment or patient handoff.

4. **Demonstrate FHIR R4 bundle export.** Click "Export FHIR" to generate a standards-
   compliant FHIR R4 Bundle resource containing:
   - Patient resource
   - Observation resources (one per biomarker)
   - DiagnosticReport resource
   - MedicationStatement resources

   > *Talking point:* "The FHIR R4 bundle means this report can be ingested directly
   > by any FHIR-compliant EHR system. We're not generating a PDF that sits in a
   > drawer -- we're producing structured, interoperable clinical data."

---

### Tab 7: Patient 360

**Purpose:** Show a unified view that brings together all patient data on a single
screen.

**Steps:**

1. **Load the Patient 360 view.** The dashboard consolidates:
   - Biomarker summary with sparkline trends
   - Genomic variant highlights (pathogenic/likely pathogenic)
   - Active medications with PGx annotations
   - Risk factor radar chart (cardiovascular, metabolic, liver, kidney, etc.)
   - Biological age gauge

2. **Demonstrate the integration story.** Point out how data from every other tab
   feeds into this single view.

   > *Talking point:* "A clinician has 15 minutes per patient visit. Patient 360
   > gives them the full picture without clicking through 8 separate systems. Every
   > insight from biomarkers, genomics, PGx, and disease risk is synthesized into
   > one actionable view."

3. **Show drill-down capability.** Click on any element in the 360 view to navigate
   to the detailed tab for that data domain.

4. **Highlight the risk factor radar chart.** The radar chart plots relative risk
   across all 9 disease categories, making it immediately visible where the patient's
   greatest vulnerabilities lie.

---

### Tab 8: Longitudinal

**Purpose:** Demonstrate biomarker trend tracking and rate-of-change analysis over
time.

**Steps:**

1. **Select a biomarker for trend analysis.** Choose HbA1c or LDL cholesterol for
   the most visually compelling demo.

2. **Review the trend chart.** The system plots historical values with:
   - Reference range bands (green/yellow/red zones)
   - Genotype-adjusted reference ranges (shifted bands)
   - Trend line with confidence interval
   - Medication start/stop annotations on the timeline

3. **Show rate-of-change analysis.** The system calculates:
   - Absolute change per unit time
   - Percentage change
   - Acceleration/deceleration of the trend
   - Projected future values (with confidence bands)

4. **Demonstrate clinical significance detection.** The system flags when a
   rate of change is clinically significant even if the absolute value is still
   within the reference range.

   > *Talking point:* "A single lab value is a snapshot. Longitudinal tracking turns
   > it into a movie. An LDL of 125 is 'normal' -- but if it was 95 six months ago,
   > that trajectory matters. Our rate-of-change analysis catches these trends before
   > they become clinical events."

---

## 5. API Demo

For technical audiences, demonstrate the FastAPI backend directly. All endpoints
are available at `http://localhost:8529`.

### 5.1 Biomarker Analysis

```bash
curl -s -X POST http://localhost:8529/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "include_genotype_adjustments": true,
    "include_critical_values": true
  }' | python3 -m json.tool
```

### 5.2 Biological Age

```bash
curl -s -X POST http://localhost:8529/biological-age \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "algorithms": ["phenoage", "grimage_surrogate"]
  }' | python3 -m json.tool
```

### 5.3 Disease Risk

```bash
curl -s -X POST http://localhost:8529/disease-risk \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "categories": [
      "cardiovascular",
      "metabolic_diabetes",
      "liver",
      "kidney",
      "thyroid",
      "iron",
      "nutritional",
      "cognitive",
      "inflammatory"
    ],
    "include_medication_context": true
  }' | python3 -m json.tool
```

### 5.4 Pharmacogenomic Profile

```bash
curl -s -X POST http://localhost:8529/pgx \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "HCLS-BIO-2026-00001",
    "genes": [
      "CYP2D6", "CYP2C19", "CYP2C9", "CYP3A5",
      "SLCO1B1", "VKORC1", "MTHFR", "TPMT", "DPYD"
    ],
    "include_cpic_guidance": true
  }' | python3 -m json.tool
```

### 5.5 Evidence Query (RAG)

```bash
curl -s -X POST http://localhost:8529/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the clinical significance of SLCO1B1 decreased function for statin therapy?",
    "collections": ["pharmacogenomics", "cardiovascular"],
    "top_k": 5,
    "patient_id": "HCLS-BIO-2026-00001"
  }' | python3 -m json.tool
```

### 5.6 Interactive API Docs

Direct the audience to `http://localhost:8529/docs` for the full Swagger/OpenAPI
interactive documentation.

---

## 6. Talking Points for Each Feature

### Biomarker Analysis
- "We analyze 67+ biomarkers with genotype-aware reference range adjustments --
  something no standard lab report provides."
- "Critical value detection with patient-specific context reduces alert fatigue
  by filtering out false positives from population-level thresholds."
- "Discordance detection catches inconsistencies across related biomarkers that
  a single-test view would miss."

### Biological Age
- "PhenoAge distills 9 clinical biomarkers into a single biological age estimate
  validated against mortality outcomes in large cohort studies."
- "Age acceleration is one of the strongest predictors of all-cause mortality --
  stronger than any individual biomarker."
- "This gives patients and clinicians a concrete, motivating number: 'You're
  biologically 3 years younger than your age -- here's what's driving that.'"

### Disease Risk
- "We model 9 disease trajectory categories, each driven by the biomarkers most
  relevant to that organ system."
- "Medication context matters: a fasting glucose of 110 in a patient on Metformin
  tells a different story than the same value in an untreated patient."
- "Sex-specific modeling ensures that hormonal and inflammatory risk factors are
  appropriately weighted."

### PGx Profile
- "14 pharmacogenes cover the majority of clinically actionable drug-gene
  interactions per CPIC guidelines."
- "Star allele calling from real genomic data (HG002 / NA24385) means these are
  confirmed variant calls, not imputed from population frequencies."
- "The traffic-light dosing system (green/yellow/red) makes PGx actionable for
  non-specialist physicians."

### Evidence Explorer
- "RAG search across 14 vector collections provides evidence-backed interpretation,
  not just pattern matching."
- "BGE-small-en-v1.5 embeddings at 384 dimensions give us biomedically-tuned
  semantic search in a compact, fast model."
- "Citation scoring lets the clinician assess evidence quality, not just relevance."

### Reports
- "The 12-section clinical report follows a structure that mirrors how clinicians
  think: overview first, details on demand, raw data in the appendix."
- "PDF export for clinical handoff. FHIR R4 for system integration. Two output
  formats, one generation step."
- "FHIR R4 compliance means this output can flow directly into Epic, Cerner, or
  any standards-compliant EHR."

### Patient 360
- "A 15-minute patient visit demands a single-screen summary. Patient 360 delivers
  biomarkers, genomics, medications, and risk in one view."
- "The risk radar chart immediately shows where clinical attention should focus."
- "Every element is clickable -- drill down into any domain without losing context."

### Longitudinal
- "Rate-of-change analysis catches clinically significant trends while values are
  still in the 'normal' range."
- "Medication annotations on the timeline show treatment effects -- did the statin
  actually bend the LDL curve?"
- "Projected values with confidence bands support proactive intervention before
  thresholds are crossed."

---

## 7. Troubleshooting Common Demo Issues

### 7.1 Streamlit UI Not Loading (Port 8528)

**Symptom:** Browser shows "connection refused" on `localhost:8528`.

**Fix:**
```bash
# Check if Streamlit is running
lsof -i :8528

# Restart Streamlit
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
streamlit run app/biomarker_ui.py --server.port 8528
```

### 7.2 FastAPI Backend Not Responding (Port 8529)

**Symptom:** API calls return connection errors.

**Fix:**
```bash
# Check if FastAPI is running
lsof -i :8529

# Restart FastAPI
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
python -m uvicorn api.main:app --host 0.0.0.0 --port 8529
```

### 7.3 Milvus Connection Failed

**Symptom:** Evidence Explorer returns empty results; RAG queries fail.

**Fix:**
```bash
# Check Milvus health
curl -s http://localhost:9091/healthz

# If Milvus is down, restart via Docker Compose (from project root)
cd /home/adam/projects/hcls-ai-factory
docker compose -f docker-compose.dgx-spark.yml up -d milvus-standalone etcd minio
```

### 7.4 ANTHROPIC_API_KEY Not Set

**Symptom:** Report generation or RAG queries return authentication errors.

**Fix:**
```bash
# Verify the key is set
echo $ANTHROPIC_API_KEY | head -c 10

# If not set, export it
export ANTHROPIC_API_KEY="sk-ant-..."
```

### 7.5 Demo Patient Not Found

**Symptom:** Patient selector dropdown is empty or patient data fails to load.

**Fix:**
```bash
# Verify demo data files exist
ls /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent/data/

# If data is missing, reload demo patients
cd /home/adam/projects/hcls-ai-factory/ai_agent_adds/precision_biomarker_agent
python scripts/load_demo_data.py
```

### 7.6 Slow Response Times

**Symptom:** API calls or UI interactions take more than 10 seconds.

**Fix:**
- Check GPU utilization: `nvidia-smi`
- Check Milvus collection load status (collections may need to be loaded into memory)
- Verify no competing workloads on the DGX Spark

### 7.7 PDF Export Fails

**Symptom:** "Export PDF" button produces an error or empty file.

**Fix:**
```bash
# Ensure WeasyPrint or report dependencies are installed
pip install weasyprint

# Check write permissions on the output directory
ls -la /tmp/biomarker_reports/
```

### 7.8 FHIR Export Validation Errors

**Symptom:** FHIR R4 bundle fails validation in external tools.

**Fix:**
- Ensure the FHIR export uses the `fhir.resources` library for validation
- Check that all required FHIR resource fields are populated
- Validate the bundle at `https://inferno.healthit.gov/validator`

---

## Appendix: Quick Reference

### Key Stats for Slides

| Metric                        | Value                    |
|-------------------------------|--------------------------|
| Biomarkers analyzed           | 67+                      |
| Vector collections (RAG)     | 14 (13 biomarker + 1 genomic) |
| Pharmacogenes profiled        | 14                       |
| Disease trajectory categories | 9                        |
| Report sections               | 12                       |
| Embedding model               | BGE-small-en-v1.5        |
| Embedding dimensions          | 384                      |
| Test suite                    | 709 tests, all passing   |
| Streamlit UI port             | 8528                     |
| FastAPI backend port          | 8529                     |

### Recommended Demo Script Timing

| Segment                     | Duration |
|-----------------------------|----------|
| Introduction & context      | 2 min    |
| Tab 1: Biomarker Analysis   | 3 min    |
| Tab 2: Biological Age       | 2 min    |
| Tab 3: Disease Risk         | 2 min    |
| Tab 4: PGx Profile          | 3 min    |
| Tab 5: Evidence Explorer    | 2 min    |
| Tab 6: Reports              | 2 min    |
| Tab 7: Patient 360          | 2 min    |
| Tab 8: Longitudinal         | 2 min    |
| API Demo (optional)         | 3 min    |
| Q&A                         | 5 min    |
| **Total**                   | **~25 min (20 min without API/Q&A)** |

---

*This guide accompanies the Precision Biomarker Intelligence Agent, part of the
HCLS AI Factory platform. For architecture documentation, see the main project
README. For deployment instructions, see the Dockerfile and docker-compose.yml
in the agent root directory.*
