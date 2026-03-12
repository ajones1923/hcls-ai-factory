# Precision Oncology Intelligence Agent - Live Demo Guide

**Total demo time:** ~30 minutes (7 scenarios + opening/closing)
**Audience:** Clinical informaticists, oncologists, pharma R&D, hospital IT leadership
**Presenter:** Ensure you have rehearsed each scenario at least once before the live session

---

## Table of Contents

1. [Overview](#1-overview)
2. [Pre-Demo Checklist](#2-pre-demo-checklist)
3. [Opening Hook (~2 min)](#3-opening-hook-2-min)
4. [Scenario 1: Case Creation & MTB Packet (~5 min)](#4-scenario-1-case-creation--mtb-packet-5-min)
5. [Scenario 2: Evidence Explorer (~5 min)](#5-scenario-2-evidence-explorer-5-min)
6. [Scenario 3: Comparative Analysis (~5 min)](#6-scenario-3-comparative-analysis-5-min)
7. [Scenario 4: Resistance Intelligence (~3 min)](#7-scenario-4-resistance-intelligence-3-min)
8. [Scenario 5: Trial Matching (~3 min)](#8-scenario-5-trial-matching-3-min)
9. [Scenario 6: FHIR Export (~3 min)](#9-scenario-6-fhir-export-3-min)
10. [Scenario 7: Cross-Agent Bridge (~2 min)](#10-scenario-7-cross-agent-bridge-2-min)
11. [Advanced Features (~5 min)](#11-advanced-features-5-min)
12. [Closing: The Big Picture](#12-closing-the-big-picture)
13. [Troubleshooting](#13-troubleshooting)
14. [Quick Reference Card](#14-quick-reference-card)

---

## 1. Overview

The Precision Oncology Intelligence Agent is a RAG-powered clinical decision support
system for Molecular Tumor Boards (MTBs). It ingests evidence from 11 Milvus vector
collections (609 vectors across literature, trials, variants, therapies, resistance
mechanisms, biomarkers, pathways, guidelines, outcomes, cases, and genomic evidence)
and synthesizes actionable recommendations using Claude AI.

### What the audience will see

- **Streamlit MTB Workbench** at `http://localhost:8526` with 5 interactive tabs
- **FastAPI backend** at `http://localhost:8527` with Swagger docs at `/docs`
- Real-time RAG retrieval across 11 collections with collection-colored source badges
- Therapy ranking with AMP/ASCO/CAP evidence tiering (Levels A through E)
- Clinical trial matching with NCT IDs and eligibility scoring
- Multi-format export: Markdown, JSON, PDF, and FHIR R4 bundles
- Cross-agent event propagation (ONCOLOGY_CASE_CREATED, THERAPY_RANKED)

### Demo patient profiles

| Patient | Cancer Type | Key Variant | Stage | Primary Therapy |
|---------|------------|-------------|-------|-----------------|
| PT-NSCLC-001 | Non-Small Cell Lung Cancer (NSCLC) | EGFR L858R | IV | Osimertinib |
| PT-BRCA-001 | Breast Cancer | BRCA1 germline | IIIA | Olaparib |
| PT-CRC-001 | Colorectal Cancer | KRAS G12C | IV | Sotorasib |
| PT-MEL-001 | Melanoma | BRAF V600E | IIIC | Dabrafenib + Trametinib |

---

## 2. Pre-Demo Checklist

Run these checks 15 minutes before the demo. Every command should succeed.

### 2.1 Verify services are running

```bash
# Check FastAPI health endpoint
curl -s http://localhost:8527/health | python3 -m json.tool
```

**Expected output:** `"status": "healthy"` with all services showing `true`:
- milvus
- embedder
- rag_engine
- intelligence_agent
- case_manager
- trial_matcher
- therapy_ranker

### 2.2 Verify Milvus collections

```bash
# List all collections and their vector counts
curl -s http://localhost:8527/collections | python3 -m json.tool
```

**Expected output:** 11 collections totaling 609 vectors:
- onco_literature
- onco_trials
- onco_variants
- onco_biomarkers
- onco_therapies
- onco_pathways
- onco_guidelines
- onco_resistance
- onco_outcomes
- onco_cases
- genomic_evidence

### 2.3 Verify Streamlit UI

```bash
# Confirm Streamlit is serving on port 8526
curl -s -o /dev/null -w "%{http_code}" http://localhost:8526
```

**Expected output:** `200`

### 2.4 Open browser tabs

**Open these tabs in order before the demo begins:**

1. **Tab 1:** `http://localhost:8526` -- Streamlit MTB Workbench (main demo surface)
2. **Tab 2:** `http://localhost:8527/docs` -- FastAPI Swagger UI (for API demos)
3. **Tab 3:** `http://localhost:8527/health` -- Health endpoint (quick status reference)

### 2.5 Test a quick RAG query

```bash
curl -s -X POST http://localhost:8527/api/ask \
  -H "Content-Type: application/json" \
  -d '{"question": "What is osimertinib?", "top_k": 3}' \
  | python3 -m json.tool | head -20
```

**Expected output:** A structured answer with sources from multiple collections and a
confidence score above 0.5.

### 2.6 Clear any stale session state

Open the Streamlit UI and refresh the page (Ctrl+R) to reset session state. Confirm
the sidebar shows **API: healthy** with green indicators for all services.

---

## 3. Opening Hook (~2 min)

> **Say this:**
> "Imagine a patient walks into your clinic with advanced non-small cell lung cancer.
> Their tumor has been sequenced, and the report shows an EGFR L858R mutation. Today,
> it takes your molecular tumor board days -- sometimes weeks -- to compile the
> evidence, match clinical trials, rank therapies, and produce a recommendation.
>
> What if that entire workflow could happen in under 5 seconds?
>
> That is what this agent does. It searches across 11 specialized knowledge collections
> -- literature, clinical trials, resistance mechanisms, treatment guidelines, and more --
> using vector similarity search, then synthesizes everything through Claude AI into
> an actionable MTB packet. Let me show you."

**Action:** Point to the Streamlit sidebar showing the green health indicators and the
total vector count (609 vectors across 11 collections).

---

## 4. Scenario 1: Case Creation & MTB Packet (~5 min)

**Goal:** Create an NSCLC patient case, generate an MTB packet, and show therapy
ranking with evidence-level badges.

### Step 1: Navigate to Case Workbench

**Action:** Click the **Case Workbench** tab in the Streamlit UI.

> **Say this:**
> "The Case Workbench is where clinicians enter patient data. It supports manual
> variant entry or direct VCF paste from annotation tools like SnpEff, VEP, or
> any VCF with GENEINFO fields."

### Step 2: Enter patient demographics

**Action:** Fill in the following fields:

| Field | Value |
|-------|-------|
| Patient ID | `PT-NSCLC-001` |
| Cancer Type | `Non-Small Cell Lung Cancer (NSCLC)` |
| Stage | `IV` |

### Step 3: Enter biomarkers

**Action:** Fill in the biomarker panel on the right side:

| Field | Value |
|-------|-------|
| TMB (mut/Mb) | `8.5` |
| MSI Status | `MSS` |
| PD-L1 TPS (%) | `60` |
| HRD Score | `0.0` |

> **Say this:**
> "The workbench captures four key biomarkers: tumor mutation burden, microsatellite
> instability, PD-L1 expression, and HRD score. These directly influence therapy
> ranking and trial eligibility."

### Step 4: Add variants

**Action:** In the Variants section, ensure **Manual Entry** is selected. Enter:

| Gene | Variant | Type |
|------|---------|------|
| EGFR | L858R | SNV |

**Action:** Click **+ Add Variant** and add a second variant:

| Gene | Variant | Type |
|------|---------|------|
| TP53 | R248W | SNV |

### Step 5: Select prior therapies

**Action:** From the **Prior Therapies** dropdown, select:
- `Platinum-based chemotherapy`

### Step 6: Create the case

**Action:** Click the green **Create Case** button.

**Expected output:** A success banner showing the generated case ID, plus a JSON
response with case details including `variant_count: 2` and the assigned `case_id`.

> **Say this:**
> "The case is now stored in Milvus as a vector embedding, which means future cases
> with similar molecular profiles can be retrieved by similarity search. This builds
> institutional memory over time."

### Step 7: Generate MTB Packet

**Action:** Click **Generate MTB Packet**.

**Expected output:** A structured MTB packet with four sections:

1. **Actionable Variants** -- Table showing EGFR L858R and TP53 R248W
2. **Evidence Summary** -- Citations from multiple collections
3. **Therapy Ranking** -- Osimertinib at the top with Evidence Level `A`, followed by
   erlotinib, gefitinib, and immunotherapy options
4. **Matched Clinical Trials** -- NCT IDs with match scores

> **Say this:**
> "Notice osimertinib ranks first with Level A evidence -- that means FDA-approved,
> standard of care for this exact molecular profile. The agent searched across
> variants, therapies, guidelines, and literature to reach this ranking.
> The resistance warnings below each therapy flag known escape mechanisms --
> for osimertinib, that includes C797S and MET amplification."

### Step 8: Point out cross-agent event

> **Say this:**
> "Behind the scenes, creating this case published an ONCOLOGY_CASE_CREATED event
> to the cross-agent event bus. Other agents in the HCLS AI Factory -- the biomarker
> agent, the CAR-T agent -- can subscribe to these events and trigger their own
> analyses automatically."

---

## 5. Scenario 2: Evidence Explorer (~5 min)

**Goal:** Demonstrate free-text clinical Q&A with multi-collection evidence retrieval,
source badges, confidence scoring, and follow-up suggestions.

### Step 1: Navigate to Evidence Explorer

**Action:** Click the **Evidence Explorer** tab.

### Step 2: Ask a clinical question

**Action:** Type the following question into the search box:

```
What is the evidence for pembrolizumab in MSI-H colorectal cancer?
```

**Action:** Click **Ask**.

### Step 3: Walk through the results

**Expected output:**

1. **Answer section** -- A synthesized paragraph covering KEYNOTE-177 trial data,
   FDA approval, response rates, and biomarker requirements
2. **Confidence bar** -- Visual progress bar showing confidence (expect > 70%)
3. **Evidence Sources** -- Multiple source cards, each with:
   - A **collection-colored badge** (blue for onco_targets, green for onco_therapies,
     red for onco_resistance, violet for onco_pathways, orange for onco_biomarkers,
     gray for onco_trials)
   - A **similarity score** (0.000 to 1.000)
   - A **text excerpt** from the matched vector
4. **Processing time** -- Displayed at the bottom in milliseconds
5. **Follow-up Questions** -- Clickable suggestions for deeper exploration

> **Say this:**
> "Each source is tagged with its collection origin -- notice the green badges are from
> the therapy collection, the gray badges from clinical trials. The agent searched all
> 11 collections in parallel, merged the results by cosine similarity, and then
> synthesized an answer using Claude with full citations.
>
> The follow-up questions are generated contextually. Clicking one immediately runs a
> new search. This creates a conversational exploration flow that mirrors how
> oncologists actually think through a case."

### Step 4: Click a follow-up question

**Action:** Click one of the suggested follow-up questions (typically something like
"What are the biomarker requirements for pembrolizumab eligibility?").

**Expected output:** A new answer with fresh sources, demonstrating the conversational
RAG flow.

### Step 5: Use filters

**Action:** Expand the **Filters** section. Set:
- Cancer Type Filter: `Colorectal Cancer`
- Gene Filter: `KRAS`

**Action:** Type a new question:

```
What are the approved therapies for KRAS G12C mutant colorectal cancer?
```

**Action:** Click **Ask**.

> **Say this:**
> "Filters narrow the search to specific cancer types or genes. This is useful when
> the knowledge base grows -- you can focus the RAG retrieval on exactly the domain
> you need."

---

## 6. Scenario 3: Comparative Analysis (~5 min)

**Goal:** Show the agent's ability to perform structured drug-vs-drug comparisons
using multi-collection evidence synthesis.

### Step 1: Stay in Evidence Explorer

**Action:** Clear any filters (set both to blank).

### Step 2: Ask a comparative question

**Action:** Type:

```
Compare osimertinib vs erlotinib for EGFR-mutant NSCLC
```

**Action:** Click **Ask**.

> **Say this:**
> "The agent detects comparative keywords -- 'compare,' 'versus,' 'vs' -- and routes
> the query through its intelligence agent rather than the simple RAG path. This
> triggers multi-step reasoning: it retrieves evidence for each drug separately,
> then synthesizes a structured comparison."

### Step 3: Walk through the comparison

**Expected output:** A structured comparison covering:

- **Efficacy data** -- PFS, OS, ORR for each drug
- **Mechanism differences** -- 3rd-gen vs 1st-gen TKI, T790M coverage
- **Resistance profiles** -- C797S for osimertinib vs T790M for erlotinib
- **Guideline preference** -- NCCN preferred first-line recommendation
- **Side effect profiles** -- Differences in rash, diarrhea, ILD rates

> **Say this:**
> "This is not a simple keyword search. The agent pulled from therapies, resistance
> mechanisms, guidelines, and literature collections to build a side-by-side
> comparison. Notice it cites specific trial data -- FLAURA for osimertinib,
> EURTAC and ENSURE for erlotinib."

### Step 4: Try a second comparison

**Action:** Type:

```
Compare dabrafenib plus trametinib vs vemurafenib for BRAF V600E melanoma
```

**Action:** Click **Ask**.

> **Say this:**
> "Combination therapies are compared against monotherapy alternatives. The agent
> understands that dabrafenib plus trametinib is a BRAF+MEK inhibitor combination
> and evaluates the added benefit of MEK inhibition for resistance prevention."

---

## 7. Scenario 4: Resistance Intelligence (~3 min)

**Goal:** Show how the agent surfaces known resistance mechanisms and suggests
alternative treatment strategies.

### Step 1: Stay in Evidence Explorer

### Step 2: Ask about resistance

**Action:** Type:

```
What resistance mechanisms emerge after osimertinib therapy?
```

**Action:** Click **Ask**.

### Step 3: Walk through resistance data

**Expected output:** A detailed answer covering:

- **T790M bypass** -- The primary mechanism for 1st/2nd-gen TKI resistance, already
  addressed by osimertinib
- **C797S mutation** -- The dominant on-target resistance mechanism for osimertinib,
  occurring in the cysteine residue that osimertinib covalently binds
- **MET amplification** -- Off-target bypass pathway, addressable with capmatinib or
  tepotinib
- **HER2 amplification** -- Alternative bypass signaling
- **Small cell transformation** -- Histologic transformation requiring platinum + etoposide
- **PIK3CA mutations** -- Parallel pathway activation

> **Say this:**
> "Resistance intelligence is critical for treatment planning. The agent pulls from the
> dedicated onco_resistance collection, which maps each therapy to its known escape
> mechanisms, bypass pathways, and alternative treatment strategies.
>
> For a molecular tumor board, this means you can proactively discuss what to do
> when -- not if -- resistance emerges. C797S after osimertinib? Consider amivantamab
> plus lazertinib. MET amplification? Add capmatinib. The agent surfaces these
> strategies with evidence levels."

---

## 8. Scenario 5: Trial Matching (~3 min)

**Goal:** Match a BRCA1-positive breast cancer patient to relevant clinical trials
with eligibility scoring.

### Step 1: Navigate to Trial Finder

**Action:** Click the **Trial Finder** tab.

### Step 2: Enter patient profile

**Action:** Fill in:

| Field | Value |
|-------|-------|
| Cancer Type | `Breast Cancer` |
| Stage | `IIIA` |
| Patient Age | `52` |

### Step 3: Select biomarkers

**Action:** Check the following biomarker checkboxes on the right side:
- `BRCA+`
- `HER2+`

> **Say this:**
> "The Trial Finder uses a weighted scoring algorithm: 40% biomarker match, 25%
> semantic similarity to trial descriptions, 20% phase weighting (later phases score
> higher), and 15% trial status (recruiting trials preferred). These weights are
> configurable."

### Step 4: Find trials

**Action:** Click **Find Trials**.

### Step 5: Walk through results

**Expected output:** A list of matched clinical trials, each showing:

- **NCT ID** (e.g., NCT04821999)
- **Phase** (e.g., Phase III)
- **Status** (e.g., Recruiting)
- **Title** of the trial
- **Match Score** (0.00 to 1.00)
- **Eligibility badge** -- green for Eligible, orange for Potentially Eligible,
  red for Not Eligible
- **Explanation** of why the trial matched

> **Say this:**
> "Each trial result shows an eligibility badge based on how well the patient profile
> matches the trial's inclusion criteria. Green means the biomarker criteria align
> directly. Orange means some criteria match but there are unknowns. The NCT IDs
> link directly to ClinicalTrials.gov for verification."

---

## 9. Scenario 6: FHIR Export (~3 min)

**Goal:** Generate an MTB packet and export it as a FHIR R4 bundle, showing
interoperability with hospital EHR systems.

### Step 1: Use the API for this demo

**Action:** Switch to the **FastAPI Swagger UI** tab (`http://localhost:8527/docs`).

> **Say this:**
> "For interoperability, the agent supports four export formats: Markdown for human
> review, JSON for system integration, PDF for formal documentation, and FHIR R4
> for EHR interoperability. Let me show the FHIR export."

### Step 2: Export as FHIR R4

**Action:** In Swagger UI, find the endpoint:
```
GET /api/reports/{case_id}/fhir
```

**Action:** Enter the case ID from Scenario 1 (the one created for PT-NSCLC-001)
and click **Execute**.

### Step 3: Walk through the FHIR bundle

**Expected output:** A FHIR R4 Bundle (type: collection) containing:

- **Patient** resource with the patient identifier
- **Condition** resource with SNOMED CT coding for the cancer type
  (e.g., code `254637007` for NSCLC)
- **Observation** resources -- one per variant, plus TMB and MSI observations,
  each coded with LOINC:
  - `69548-6` -- Genetic variant assessment
  - `48018-6` -- Gene studied
  - `94076-7` -- Tumor mutation burden
  - `81695-9` -- Microsatellite instability
- **Specimen** resource with tissue type
- **MedicationRequest** resources (status: draft, intent: proposal) for each
  ranked therapy with evidence level annotations
- **DiagnosticReport** tying all observations together under LOINC `81247-9`
  (Master HL7 genetic variant reporting panel)

> **Say this:**
> "This is a fully compliant FHIR R4 genomics reporting bundle. The DiagnosticReport
> references all variant observations and uses SNOMED CT for the cancer diagnosis,
> LOINC for all lab observations, and HGVS notation for variants. This can be
> imported directly into any FHIR-compliant EHR -- Epic, Cerner, or any SMART on
> FHIR application.
>
> The MedicationRequest resources represent the therapy recommendations as draft
> proposals, ready for clinician review and approval within the EHR workflow."

### Step 4: Show other export formats

**Action:** Quickly demonstrate the other export endpoints:

```bash
# Markdown export
curl -s http://localhost:8527/api/reports/{case_id}/markdown | head -30

# JSON export
curl -s http://localhost:8527/api/reports/{case_id}/json | python3 -m json.tool | head -30

# PDF export (downloads a file)
curl -s http://localhost:8527/api/reports/{case_id}/pdf -o /tmp/mtb_report.pdf
```

> **Say this:**
> "Four formats, one API. Markdown for the tumor board review, JSON for downstream
> analytics, PDF for the patient record, FHIR R4 for EHR integration."

---

## 10. Scenario 7: Cross-Agent Bridge (~2 min)

**Goal:** Show how the oncology agent communicates with other agents in the
HCLS AI Factory through the event bus.

### Step 1: Navigate to Outcomes Dashboard

**Action:** Click the **Outcomes Dashboard** tab in the Streamlit UI.

### Step 2: Show the knowledge base statistics

**Expected output:** Five metric cards at the top:

- **Targets** -- Number of actionable gene targets (40+)
- **Therapies** -- Number of therapy entries
- **Resistance** -- Number of resistance mechanism records
- **Pathways** -- Number of signaling pathway entries
- **Biomarkers** -- Number of biomarker records

Below the metrics, a horizontal bar chart showing vector counts per collection.

### Step 3: Show the event log

**Action:** Scroll down to the **Recent Events** section.

**Expected output:** Event entries from the earlier demo scenarios:

- `case_created` -- From Scenario 1, with case ID, patient ID, cancer type
- `mtb_generated` -- From Scenario 1, with therapy count, trial count, processing time

> **Say this:**
> "Every significant action publishes an event to the HCLS AI Factory event bus.
> When we created the NSCLC case, an ONCOLOGY_CASE_CREATED event was published.
> When we generated the MTB packet, a THERAPY_RANKED event followed.
>
> Other agents in the factory -- the biomarker discovery agent, the CAR-T therapy
> agent, the imaging analysis agent -- can subscribe to these events and trigger
> their own workflows. For example, the biomarker agent might see a new EGFR L858R
> case and automatically queue a resistance monitoring panel. This is how the
> five agents in the factory collaborate without tight coupling."

### Step 4: Show the cross-agent API

**Action:** Switch briefly to the Swagger UI and point to:

```
GET /api/events?limit=20
```

> **Say this:**
> "The event log is accessible via API, so any downstream system -- a LIMS, a
> clinical trial management system, a pharmacy system -- can poll for events and
> act on them."

---

## 11. Advanced Features (~5 min)

Use this section to fill time or address specific audience interests.

### 11.1 VCF Upload

**Action:** In the Case Workbench, switch the Input Mode to **Paste VCF**.

**Action:** Paste the following minimal VCF content:

```
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr7	55259515	.	T	G	100	PASS	GENEINFO=EGFR;ANN=L858R|missense_variant
chr17	7578406	.	C	T	95	PASS	GENEINFO=TP53;ANN=R248W|missense_variant
```

> **Say this:**
> "The agent parses VCF files with annotations from SnpEff, VEP, or any tool that
> populates the GENEINFO field. This means you can go directly from your genomics
> pipeline output to an MTB recommendation -- no manual data entry required."

### 11.2 API-Driven Workflow

**Action:** In a terminal, demonstrate the full programmatic workflow:

```bash
# Step 1: Create a case via API
CASE_RESPONSE=$(curl -s -X POST http://localhost:8527/api/cases \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "PT-API-001",
    "cancer_type": "Melanoma",
    "stage": "IIIC",
    "variants": [
      {"gene": "BRAF", "variant": "V600E", "variant_type": "SNV"}
    ],
    "biomarkers": {"tmb": 12.0, "msi_status": "MSS", "pdl1_tps": 30},
    "prior_therapies": []
  }')
echo "$CASE_RESPONSE" | python3 -m json.tool

# Step 2: Extract case ID
CASE_ID=$(echo "$CASE_RESPONSE" | python3 -c "import sys,json; print(json.load(sys.stdin)['case_id'])")

# Step 3: Generate MTB packet
curl -s -X POST "http://localhost:8527/api/cases/${CASE_ID}/mtb" \
  -H "Content-Type: application/json" \
  -d '{}' | python3 -m json.tool

# Step 4: Rank therapies
curl -s -X POST http://localhost:8527/api/therapies/rank \
  -H "Content-Type: application/json" \
  -d '{
    "cancer_type": "Melanoma",
    "variants": [{"gene": "BRAF", "variant": "V600E", "variant_type": "SNV"}],
    "biomarkers": {"tmb": 12.0, "msi_status": "MSS", "pdl1_tps": 30},
    "prior_therapies": []
  }' | python3 -m json.tool

# Step 5: Export as FHIR R4
curl -s "http://localhost:8527/api/reports/${CASE_ID}/fhir" | python3 -m json.tool
```

> **Say this:**
> "Everything in the UI is backed by a REST API. This means the oncology agent can
> be integrated into any clinical workflow -- an EHR plugin, a LIMS webhook, or a
> Nextflow pipeline step. The entire flow from case creation to FHIR export is five
> API calls."

### 11.3 Knowledge Base Statistics

**Action:** Run:

```bash
curl -s http://localhost:8527/knowledge/stats | python3 -m json.tool
```

> **Say this:**
> "The knowledge base currently spans 11 collections with 609 vectors covering 40+
> actionable gene targets, 26 cancer types, and 9 therapy categories. The collections
> are weighted for search relevance: variants at 18%, literature at 16%, therapies
> at 14%, guidelines at 12%, and so on. These weights are tunable per institution."

### 11.4 Prometheus Metrics

**Action:** Open `http://localhost:8527/metrics` in a browser tab.

> **Say this:**
> "The agent exposes Prometheus-compatible metrics for monitoring. Each collection's
> vector count is tracked as a gauge, so you can set alerts if ingestion pipelines
> fail or data drifts."

### 11.5 Demo Mode Shortcut

**Action:** In the Streamlit sidebar, click the **Load Demo Patient** button.

> **Say this:**
> "For quick demonstrations, the sidebar has a Load Demo Patient button that
> pre-populates the Case Workbench with our reference NSCLC patient. This uses
> shared demo data from the HCLS common library, ensuring consistency across all
> five agents in the factory."

---

## 12. Closing: The Big Picture

> **Say this:**
> "Let me step back and put this in context. What you have seen today is one agent
> in a three-stage precision medicine pipeline that runs entirely on a single
> NVIDIA DGX Spark.
>
> Stage 1 takes raw patient DNA -- FASTQ files -- and runs GPU-accelerated variant
> calling with Parabricks and DeepVariant. In 2 to 4 hours, you have an annotated
> VCF file with 11.7 million variants.
>
> Stage 2 -- the RAG pipeline -- embeds those variants into Milvus alongside 3.56
> million searchable vectors from ClinVar, AlphaMissense, and the oncology knowledge
> base. That is what powers the evidence retrieval you just saw.
>
> Stage 3 -- drug discovery -- takes the top therapeutic targets and runs molecular
> docking with BioNeMo and DiffDock to generate novel drug candidates.
>
> From patient DNA to drug candidates in under 5 hours. On a single machine. That
> is the HCLS AI Factory.
>
> The oncology agent you saw today is the clinical intelligence layer -- it turns
> genomic data into actionable treatment recommendations for molecular tumor boards.
> But it does not work alone. It shares events with four other agents: biomarker
> discovery, CAR-T therapy design, imaging analysis, and autoimmune intelligence.
> Together, they cover the full spectrum of precision medicine.
>
> Questions?"

---

## 13. Troubleshooting

### Problem: Streamlit shows "Cannot connect to API"

**Cause:** The FastAPI backend is not running on port 8527.

**Fix:**
```bash
# Check if the API process is running
lsof -i :8527

# Start the API if needed
cd $HCLS_HOME/ai_agent_adds/precision_oncology_agent/agent
python -m api.main
```

### Problem: Health endpoint shows "degraded" status

**Cause:** One or more services failed to initialize (usually Milvus).

**Fix:**
```bash
# Verify Milvus is running
curl -s http://localhost:19530/v1/vector/collections

# Check Milvus container status
docker ps | grep milvus

# Restart Milvus if needed
docker restart milvus-standalone
```

### Problem: Collections show 0 vectors

**Cause:** Seed data has not been ingested.

**Fix:**
```bash
cd $HCLS_HOME/ai_agent_adds/precision_oncology_agent/agent

# Run all seed scripts
python scripts/setup_collections.py
python scripts/seed_knowledge.py
python scripts/seed_variants.py
python scripts/seed_therapies.py
python scripts/seed_biomarkers.py
python scripts/seed_trials.py
python scripts/seed_resistance.py
python scripts/seed_pathways.py
python scripts/seed_guidelines.py
python scripts/seed_outcomes.py
python scripts/seed_literature.py
python scripts/seed_cases.py
```

### Problem: RAG query returns "No answer generated"

**Cause:** The ANTHROPIC_API_KEY environment variable is not set, or the Claude API
is unreachable.

**Fix:**
```bash
# Check if the API key is set
echo $ONCO_ANTHROPIC_API_KEY

# Set it if missing (replace with your actual key)
export ONCO_ANTHROPIC_API_KEY="sk-ant-..."
```

### Problem: Evidence sources show only one collection

**Cause:** Most collections are empty. Run `curl http://localhost:8527/collections`
to check vector counts per collection.

**Fix:** Re-run the seed scripts for the empty collections (see "Collections show
0 vectors" above).

### Problem: PDF export fails with ImportError

**Cause:** ReportLab is not installed.

**Fix:**
```bash
pip install reportlab
```

### Problem: FHIR export returns a minimal stub instead of full bundle

**Cause:** The export module failed to load. Check the API logs for import errors.

**Fix:**
```bash
# Verify the export module loads
cd $HCLS_HOME/ai_agent_adds/precision_oncology_agent/agent
python -c "from src.export import export_fhir_r4; print('OK')"
```

### Problem: Cross-agent events are not appearing

**Cause:** The HCLS common event bus library is not on the Python path.

**Fix:**
```bash
export PYTHONPATH="$HCLS_HOME/lib:$PYTHONPATH"
```

### Problem: Streamlit session state is corrupted

**Cause:** Rapid clicking or browser back-navigation can desync Streamlit state.

**Fix:** Refresh the page (Ctrl+R) or append `?reset=1` to the URL.

---

## 14. Quick Reference Card

Print this card and keep it next to your laptop during the demo.

### Service URLs

| Service | URL | Purpose |
|---------|-----|---------|
| Streamlit UI | `http://localhost:8526` | MTB Workbench (5 tabs) |
| FastAPI | `http://localhost:8527` | REST API backend |
| Swagger Docs | `http://localhost:8527/docs` | Interactive API docs |
| Health Check | `http://localhost:8527/health` | Service status |
| Metrics | `http://localhost:8527/metrics` | Prometheus metrics |
| Milvus Attu | `http://localhost:8000` | Vector DB admin UI |
| Grafana | `http://localhost:3000` | Monitoring dashboards |
| HCLS Portal | `http://localhost:8510` | Landing page hub |

### API Endpoints

| Method | Endpoint | Purpose |
|--------|----------|---------|
| POST | `/api/ask` | RAG-powered clinical Q&A |
| POST | `/api/cases` | Create patient case |
| GET | `/api/cases/{id}` | Retrieve case details |
| POST | `/api/cases/{id}/mtb` | Generate MTB packet |
| GET | `/api/cases/{id}/variants` | List case variants |
| POST | `/api/trials/match` | Match clinical trials |
| POST | `/api/trials/match-case/{id}` | Match trials for existing case |
| POST | `/api/therapies/rank` | Rank therapy options |
| POST | `/api/reports/generate` | Generate report from question |
| GET | `/api/reports/{id}/{fmt}` | Export case report (md/json/pdf/fhir) |
| GET | `/api/events` | List recent events |
| GET | `/health` | Service health |
| GET | `/collections` | List Milvus collections |
| GET | `/knowledge/stats` | Knowledge base statistics |
| GET | `/metrics` | Prometheus metrics |
| POST | `/query` | Direct RAG query |
| POST | `/search` | Vector search (no LLM) |
| POST | `/find-related` | Cross-collection entity linking |

### Demo Patient Quick-Entry

**NSCLC (Scenario 1):**
- Patient ID: `PT-NSCLC-001`
- Cancer: Non-Small Cell Lung Cancer (NSCLC)
- Stage: IV
- TMB: 8.5 | MSI: MSS | PD-L1: 60% | HRD: 0
- Variants: EGFR L858R (SNV), TP53 R248W (SNV)
- Prior: Platinum-based chemotherapy

**Breast Cancer (Scenario 5):**
- Cancer: Breast Cancer
- Stage: IIIA
- Age: 52
- Biomarkers: BRCA+, HER2+

**Melanoma (Advanced):**
- Patient ID: `PT-MEL-001`
- Cancer: Melanoma
- Stage: IIIC
- TMB: 12.0 | MSI: MSS | PD-L1: 30%
- Variants: BRAF V600E (SNV)

### Key Talking Points

- **11 Milvus collections**, 609 vectors, BGE-small-en-v1.5 embeddings (dim=384)
- **40+ actionable gene targets** with AMP/ASCO/CAP evidence tiering (Levels A-E)
- **26 cancer types** supported across 20 options in the UI dropdown
- **9 therapy categories**: TKIs, immunotherapy, chemotherapy, PARP inhibitors,
  antibody-drug conjugates, MEK inhibitors, anti-angiogenics, hormone therapy,
  combination regimens
- **4 export formats**: Markdown, JSON, PDF (NVIDIA-branded), FHIR R4
- **Cross-agent events**: ONCOLOGY_CASE_CREATED, THERAPY_RANKED
- **Collection search weights**: variants 18%, literature 16%, therapies 14%,
  guidelines 12%, trials 10%, biomarkers 8%, resistance 7%, pathways 6%,
  outcomes 4%, genomic 3%, cases 2%
- **Trial matching weights**: biomarker 40%, semantic 25%, phase 20%, status 15%

### Scenario Timing

| # | Scenario | Time | Tab |
|---|----------|------|-----|
| -- | Opening Hook | 2 min | Sidebar |
| 1 | Case Creation & MTB | 5 min | Case Workbench |
| 2 | Evidence Explorer | 5 min | Evidence Explorer |
| 3 | Comparative Analysis | 5 min | Evidence Explorer |
| 4 | Resistance Intelligence | 3 min | Evidence Explorer |
| 5 | Trial Matching | 3 min | Trial Finder |
| 6 | FHIR Export | 3 min | Swagger UI |
| 7 | Cross-Agent Bridge | 2 min | Outcomes Dashboard |
| -- | Advanced Features | 5 min | Various |
| -- | Closing | 2 min | -- |
| | **Total** | **~35 min** | |

### Evidence Level Color Badges (Therapy Ranker)

| Level | Color | Meaning |
|-------|-------|---------|
| 1 / 1A / 1B | Green | FDA-approved, standard of care |
| 2 / 2A / 2B | Blue | Clinical evidence, consensus |
| 3 / 3A / 3B | Orange | Case reports, early trials |
| 4 | Red | Preclinical, biological rationale |
| R1 / R2 | Red | Resistance evidence |

### Collection Badge Colors (Evidence Explorer)

| Collection | Badge Color |
|------------|-------------|
| onco_targets / onco_variants | Blue |
| onco_therapies | Green |
| onco_resistance | Red |
| onco_pathways | Violet |
| onco_biomarkers | Orange |
| onco_trials | Gray |

---

*Last updated: March 2026*
*Author: Adam Jones*
*Pipeline: HCLS AI Factory - Precision Oncology Intelligence Agent v0.1.0*
