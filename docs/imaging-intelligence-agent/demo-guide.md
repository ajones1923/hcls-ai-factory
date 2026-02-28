---
search:
  boost: 2
tags:
  - Demo
  - Walkthrough
  - Imaging
  - Radiology
  - DICOM
  - NIM
---

# Imaging Intelligence Agent — Demo Guide

> **Step-by-step walkthrough for demonstrating the Imaging Intelligence Agent on DGX Spark.**
>
> License: Apache 2.0 | Date: February 2026

---

## Demo Overview

| Parameter | Value |
|---|---|
| **Route A Duration** | 20 minutes (standalone agent) |
| **Route B Duration** | 30-35 minutes (cross-platform integration) |
| **Hardware** | NVIDIA DGX Spark (GB10, 128 GB unified) |
| **Clinical Workflows** | 4 (CT hemorrhage, CT lung nodule, CXR findings, MRI MS lesion) |
| **NIM Services** | 4 (VISTA-3D, MAISI, VILA-M3, Llama-3 8B) |
| **Knowledge Base** | 3.56M vectors across 11 Milvus collections |
| **Export Formats** | Markdown, JSON, PDF, FHIR R4 DiagnosticReport |

### What the Audience Will See

1. A health dashboard showing 11 vector collections with 3.56 million indexed vectors
2. Four NVIDIA NIM microservices running on-device with automatic fallback logic
3. A CT head scan triaged in seconds — hemorrhage volume, midline shift, and urgency routing
4. A chest CT with lung nodules classified by ACR Lung-RADS with automated follow-up scheduling
5. A chest X-ray analyzed by DenseNet-121 with multi-label findings and GradCAM heatmaps
6. An MRI brain scan with MS lesion tracking and disease activity classification
7. RAG-grounded clinical answers from 2,678 PubMed papers with clickable citations
8. VISTA-3D zero-shot segmentation across 132 anatomical classes
9. Reports exported as Markdown, JSON, PDF, and FHIR R4 DiagnosticReport Bundles with SNOMED CT, LOINC, and DICOM coding
10. A DICOM study arriving via Orthanc and automatically routed to the correct clinical workflow
11. *(Route B)* A lung nodule triggering cross-modal genomic queries against 3.5 million variant vectors
12. *(Route B)* The complete pipeline: DICOM study → imaging AI → genomic trigger → target identification → drug candidates

---

## Pre-Demo Setup

### Step 1: Verify Hardware

```bash
# Verify DGX Spark GPU
nvidia-smi
# Expected: GB10 GPU, 128 GB unified memory

# Verify ARM64 architecture
uname -m
# Expected: aarch64
```

### Step 2: Set Environment Variables

```bash
cp .env.example .env

# Required variables:
# ANTHROPIC_API_KEY=sk-ant-...     (for Claude RAG synthesis)
# IMAGING_NVIDIA_API_KEY=nvapi-... (for NVIDIA Cloud NIM fallback)
# IMAGING_NGC_API_KEY=...          (for local NIM containers)
```

### Step 3: Start Services

```bash
cd ai_agent_adds/imaging_intelligence_agent/agent
docker compose up -d

# This starts 13 services:
# Milvus (etcd + MinIO + standalone)
# 4 NIM containers (VISTA-3D, MAISI, VILA-M3, Llama-3)
# FastAPI server (port 8524)
# Streamlit UI (port 8525)
# Orthanc DICOM server (ports 8042/4242)
# Prometheus + Grafana
```

### Step 4: Verify All Services Healthy

```bash
curl -s http://localhost:8524/health | python3 -m json.tool
```

Expected response:

```json
{
  "status": "healthy",
  "collections": {
    "imaging_literature": 2678,
    "imaging_trials": 12,
    "imaging_findings": 25,
    "imaging_protocols": 15,
    "imaging_devices": 15,
    "imaging_anatomy": 20,
    "imaging_benchmarks": 15,
    "imaging_guidelines": 12,
    "imaging_report_templates": 10,
    "imaging_datasets": 12,
    "genomic_evidence": 3561170
  },
  "total_vectors": 3563984,
  "nim_services": {
    "vista3d": "available",
    "maisi": "mock",
    "vila_m3": "cloud",
    "llm": "anthropic"
  }
}
```

All 11 collections should have non-zero counts. Total vectors should be ~3.56M.

### Step 5: Verify Orthanc DICOM Server

```bash
curl -s http://localhost:8042/system | python3 -m json.tool
# Expected: Orthanc version, DICOM AET, storage info
```

### Step 6: Open the Imaging Agent UI

Open **http://localhost:8525** in your browser. Confirm the Streamlit UI loads with five tabs (Ask, Comparative, Workflow Demo, Reports, Settings) and the sidebar shows NIM service status indicators and collection stats.

---

## Route A: Standalone Agent Demo (20 minutes)

### Route A Timeline

| Step | Time | Tab / Location | Action |
|------|------|---------------|--------|
| Opening | 1 min | Landing page :8080 | Show health grid, Imaging + NIM services green |
| NIM Status | 1 min | Sidebar :8525 | Show NIM service status indicators |
| Workflow 1 | 3 min | Workflow Demo tab | Select "CT Head -- Hemorrhage Detection", click "Run Demo" |
| Workflow 2 | 3 min | Workflow Demo tab | Select "CT Chest -- Lung Nodule Analysis", click "Run Demo" |
| Workflow 3 | 2 min | Workflow Demo tab | Select "CXR -- Rapid Findings Triage", click "Run Demo" |
| Workflow 4 | 2 min | Workflow Demo tab | Select "MRI Brain -- MS Lesion Quantification", click "Run Demo" |
| Evidence Q&A | 3 min | Ask tab | Type query about AI-assisted CT sensitivity |
| Comparative | 2 min | Comparative tab | Compare CT vs MRI for brain hemorrhage detection |
| Reports | 1 min | Reports tab | Export PDF report |
| Closing | 2 min | — | Talking points |

---

### Opening (1 minute)

**Show:** Landing page at http://localhost:8080 — highlight the service health grid with Imaging Agent and NIM services showing green.

**Talking points:**

- "This is the Imaging Intelligence Agent — it processes CT, MRI, and chest X-ray studies using NVIDIA NIM microservices and a 3.56 million vector knowledge base."
- "Four clinical workflows run real pretrained model weights: SegResNet, RetinaNet, DenseNet-121, and UNEST."
- "Everything runs on a single DGX Spark — a $3,999 desktop workstation with 128 GB unified memory."

---

### NIM Service Status (1 minute)

**Show:** Imaging Agent UI at http://localhost:8525

**Look at:** The sidebar on the left. The **NIM Services Status** section displays four services with color-coded indicators:

- **VISTA-3D** — green (available locally)
- **MAISI** — yellow (mock mode)
- **VILA-M3** — green (cloud fallback)
- **Llama-3** — green (Anthropic Claude)

**Expected result:** Four NIM services listed with green/yellow indicators. No red (unavailable) indicators.

**Talking points:**

- "Four NVIDIA NIM microservices: VISTA-3D for 3D segmentation, MAISI for synthetic CT generation, VILA-M3 for vision-language understanding, and Llama-3 for text generation."
- "The system supports three modes — local GPU, NVIDIA Cloud, or clinically realistic mock — with automatic fallback. No single point of failure."

**Look at:** Below the NIM status, the **Collection Stats** section shows all 10 imaging collections plus the shared `genomic_evidence` collection. Point out the total vector count of ~3.56M.

---

### Workflow 1: CT Head Hemorrhage Triage (3 minutes)

**Click:** The **Workflow Demo** tab in the main content area.

**Select:** From the workflow selectbox, choose **"CT Head -- Hemorrhage Detection"**.

**Expected result:** The metrics display shows:

- **Modality:** CT
- **Body Region:** Head
- **Target Latency:** < 90 sec

**Click:** The **"Run Demo"** button.

**Expected result:** The result panel displays with an **urgent** severity indicator (red) and the following findings and measurements:

```json
{
  "workflow_name": "ct_head_hemorrhage",
  "status": "completed",
  "findings": [
    {
      "category": "hemorrhage",
      "description": "Intraparenchymal hemorrhage in right basal ganglia, volume 12.5 mL, midline shift 3.2 mm",
      "severity": "urgent",
      "recommendation": "Surgical consultation recommended"
    }
  ],
  "measurements": {
    "volume_ml": 12.5,
    "midline_shift_mm": 3.2,
    "max_thickness_mm": 8.1,
    "hounsfield_mean": 62.0,
    "hounsfield_max": 78.0,
    "surrounding_edema_ml": 4.3
  },
  "classification": "urgent_hemorrhage",
  "severity": "urgent",
  "inference_time_ms": 42.5,
  "nim_services_used": ["vista3d"]
}
```

**Talking points:**

- "A CT head scan is triaged in under 90 seconds. The SegResNet model segments the hemorrhage and measures volume, midline shift, and maximum thickness."
- "12.5 mL with 3.2 mm midline shift — this is classified as 'urgent' using Brain Trauma Foundation criteria."
- "Critical cases (>30 mL or >5 mm shift) get immediate surgical routing. This one gets urgent neurosurgical consultation."

**Key metrics:**

| Measurement | Value | Clinical Significance |
|---|---|---|
| Hemorrhage volume | 12.5 mL | Above 5 mL threshold --> urgent |
| Midline shift | 3.2 mm | Below 5 mm critical threshold |
| Max thickness | 8.1 mm | Localized hemorrhage |
| Hounsfield (mean) | 62 HU | Consistent with acute blood |

---

### Workflow 2: CT Chest Lung Nodule Tracking (3 minutes)

**Select:** From the workflow selectbox, choose **"CT Chest -- Lung Nodule Analysis"**.

**Click:** The **"Run Demo"** button.

**Expected result:** The result panel shows nodule detection with Lung-RADS classification, volume measurements, and follow-up scheduling.

**Talking points:**

- "RetinaNet detects lung nodules and SegResNet segments each one for precise volume measurement."
- "ACR Lung-RADS v2022 classification — the same system radiologists use — automatically applied."
- "Lung-RADS 4A or higher triggers cross-modal genomic queries. We'll see that in Route B."
- "Volume doubling time is calculated for longitudinal tracking using diffeomorphic registration."

---

### Workflow 3: CXR Rapid Findings (2 minutes)

**Select:** From the workflow selectbox, choose **"CXR -- Rapid Findings Triage"**.

**Click:** The **"Run Demo"** button.

**Expected result:** The result panel shows multi-label findings with confidence scores and GradCAM attention regions.

**Talking points:**

- "DenseNet-121 pretrained on CheXpert — multi-label classification in under 30 seconds."
- "Findings include consolidation, effusion, pneumothorax, cardiomegaly, atelectasis, edema, and nodules."
- "Each finding includes a confidence score. GradCAM heatmaps show exactly where the model is looking."
- "Critical findings like tension pneumothorax get immediate urgency routing."

---

### Workflow 4: MRI Brain MS Lesion Tracking (2 minutes)

**Select:** From the workflow selectbox, choose **"MRI Brain -- MS Lesion Quantification"**.

**Click:** The **"Run Demo"** button.

**Expected result:** The result panel shows lesion counts, total lesion volume, disease activity classification, and longitudinal comparison.

**Talking points:**

- "UNEST segments white matter lesions from FLAIR MRI sequences."
- "Longitudinal matching tracks individual lesions across timepoints using ANTsPy diffeomorphic registration."
- "Disease activity is classified as Stable, Active, or Highly Active based on new and enlarging lesion counts."
- "Highly Active disease triggers neurological genomics queries — HLA-DRB1 and demyelination markers."

---

### RAG Knowledge Query (3 minutes)

**Click:** The **Ask** tab.

**Type this query:**
> "What is the sensitivity of AI-assisted CT for pulmonary nodule detection?"

**Expected result:** The main panel displays Claude's synthesized answer with grounded citations from PubMed and clinical guidelines. The **Evidence sidebar** on the right shows a collection breakdown — Literature, Guidelines, Trials, Benchmarks — each with relevance scores.

**Talking points:**

- "The RAG engine searches across 11 Milvus collections simultaneously — 2,678 PubMed papers, clinical guidelines, trial results, device clearances, and benchmarks."
- "BGE-small-en-v1.5 embeddings with asymmetric query encoding — retrieval completes in under 20 milliseconds."
- "Claude synthesizes a grounded answer with clickable PubMed and ClinicalTrials.gov citations."
- "Follow-up questions are automatically generated based on the topic area."

**Show:** Point out the evidence sidebar — highlight the collection labels (Literature, Guidelines, Trials, Benchmarks) and the relevance scores next to each retrieved chunk.

---

### Comparative Analysis (2 minutes)

**Click:** The **Comparative** tab.

**Type in Entity A:**
> "CT"

**Type in Entity B:**
> "MRI"

**Type in the comparison question field:**
> "brain hemorrhage detection"

**Click:** The **"Compare"** button.

**Expected result:** A side-by-side evidence display comparing CT and MRI for brain hemorrhage detection, with evidence pulled from relevant collections and relevance scores for each side.

**Talking points:**

- "Side-by-side comparative analysis grounded in the same 3.56 million vector knowledge base."
- "CT remains the gold standard for acute hemorrhage detection — faster acquisition, higher sensitivity for acute blood."
- "MRI offers superior contrast resolution for subacute and chronic hemorrhage, and detects microbleeds CT misses."
- "Each side shows its own evidence panel with collection labels and relevance scores."

---

### Multi-Format Report Export (1 minute)

**Click:** The **Reports** tab.

**Click:** The **"Export PDF"** button.

**Expected result:** A PDF report is generated containing the workflow results from this session. A **download button** appears.

**Click:** The download button to save the PDF locally.

**Talking points:**

- "Four export formats: Markdown for human reading, JSON for dashboards, PDF for clinical documentation, and FHIR R4 for EHR integration."
- "The FHIR R4 DiagnosticReport Bundle includes Patient, ImagingStudy, Observation, and DiagnosticReport resources."
- "Every finding is coded with SNOMED CT, LOINC, and DICOM standard terminologies — ready for Epic, Cerner, or any FHIR-compliant EHR."
- "Measurements include UCUM units — milliliters for volume, millimeters for shift."

**Key export formats:**

| Format | Use Case | Standards |
|---|---|---|
| Markdown | Human-readable clinical report | — |
| JSON | Programmatic consumption, dashboards | — |
| PDF | Clinical documentation, patient records | ReportLab styled |
| FHIR R4 | EHR integration, data exchange | SNOMED CT, LOINC, DICOM, HL7 v3, UCUM |

---

### Closing Route A (2 minutes)

**Talking points:**

- "Four clinical workflows, four NVIDIA NIMs, 11 knowledge collections, 3.56 million vectors — all on a single $3,999 DGX Spark."
- "From DICOM ingestion to clinical report in under 90 seconds for the most urgent cases."
- "Mock mode works identically to live mode — swap in real DICOM data and the same pipeline runs real model weights."

---

## Route B: Cross-Platform Integration Demo (30-35 minutes)

> **Prerequisite:** Complete Route A first, or start from a fresh session with all services running.
>
> This route demonstrates how the Imaging Intelligence Agent connects to the full HCLS AI Factory — genomics, RAG-grounded target identification, and AI-driven drug discovery — creating a closed-loop precision medicine workflow.

### Route B Timeline

| Step | Time | Tab / Location | Action |
|------|------|---------------|--------|
| Platform Overview | 2 min | Landing page :8080 | Show health grid, all platform services |
| Shared Data Layer | 2 min | Sidebar :8525 | Show collection stats, highlight genomic_evidence |
| Lung Nodule + Cross-Modal | 7 min | Workflow Demo tab :8525 | Select "CT Chest -- Lung Nodule Analysis", show cross-modal trigger |
| FHIR Export | 2 min | Reports tab :8525 | Export report with genomic context |
| Bridge to Stage 2 | 5 min | RAG Chat :8501 | Type genomic bridge query |
| Bridge to Stage 3 | 4 min | Drug Discovery :8505 | Show EGFR target pipeline |
| DICOM Auto-Routing | 3 min | — | Talking points about Orthanc webhook architecture |
| Complete Pipeline | 3 min | — | Walk through full pipeline diagram |
| Closing | 2 min | — | Talking points |

---

### Platform Overview (2 minutes)

**Show:** HCLS AI Factory landing page at http://localhost:8080

**Talking points:**

- "The HCLS AI Factory is a 3-stage precision medicine platform: GPU-accelerated genomics, RAG-grounded target identification, and AI-driven drug discovery."
- "The Imaging Intelligence Agent is one of several AI agents that extend this platform."
- "All agents share the same infrastructure — Milvus vector database, BGE embeddings, Claude LLM — creating a unified intelligence layer."

**Service overview:**

| Service | Port | Role |
|---|---|---|
| Landing Page | 8080 | Platform health dashboard |
| Imaging Agent API | 8524 | This demo's primary endpoint |
| Imaging Agent UI | 8525 | Streamlit chat interface |
| RAG Chat (Stage 2) | 8501 | Genomic variant analysis |
| Drug Discovery (Stage 3) | 8505 | Molecule generation + docking |
| Milvus | 19530 | Shared vector database |
| Orthanc DICOM | 8042 | Medical image server |

---

### Shared Data Layer (2 minutes)

**Show:** Imaging Agent UI at http://localhost:8525

**Look at:** The sidebar **Collection Stats** section. Point out each of the 10 imaging collections and the shared `genomic_evidence` collection.

**Talking points:**

- "Look at the `genomic_evidence` collection — 3,561,170 vectors. These are real patient variants from the genomics pipeline."
- "The Imaging Agent reads this collection in read-only mode. The same vectors are used by Stage 2 RAG, the CAR-T Agent, and now imaging."
- "This is the connective tissue of the platform — every agent sees the same genomic truth."

---

### CT Chest Lung Nodule with Cross-Modal Trigger (7 minutes)

#### Step 1: Run the lung nodule workflow

**Click:** The **Workflow Demo** tab.

**Select:** From the workflow selectbox, choose **"CT Chest -- Lung Nodule Analysis"**.

**Click:** The **"Run Demo"** button.

**Expected result:** The result panel shows lung nodule detection with Lung-RADS classification, volume, and follow-up scheduling.

**Talking points:**

- "This CT chest study shows a suspicious lung nodule. RetinaNet detects it, SegResNet segments it."
- "The Lung-RADS classification comes back as 4A — that's a high-risk nodule requiring tissue sampling."

#### Step 2: Cross-modal trigger fires automatically

**Expected result:** The workflow result includes a `cross_modal` section in the output. Point this out in the result display.

```json
{
  "cross_modal": {
    "trigger_reason": "Lung-RADS 4A — high-risk lung nodule",
    "genomic_context": [
      "EGFR exon 19 deletion (del19) — sensitive to erlotinib, gefitinib",
      "KRAS G12C — sotorasib responsive",
      "ALK fusion (EML4-ALK) — alectinib, lorlatinib first-line"
    ],
    "genomic_hit_count": 12,
    "query_count": 3,
    "enrichment_summary": "Patient variants suggest actionable lung cancer targets with approved targeted therapies"
  }
}
```

**Talking points:**

- "Here's where it gets interesting. Lung-RADS 4A or higher automatically triggers genomic queries."
- "The system queries the `genomic_evidence` collection — 3.5 million real variant vectors — for EGFR, ALK, ROS1, and KRAS mutations."
- "12 genomic hits across 3 queries. The system found EGFR del19, KRAS G12C, and ALK fusion variants."
- "Each variant is linked to specific targeted therapies — this isn't generic information, it's precision medicine."
- "The imaging finding triggered genomic analysis without any human intervention."

---

### FHIR R4 Export with Genomic Context (2 minutes)

**Click:** The **Reports** tab.

**Click:** The **"Export JSON"** button.

**Expected result:** A FHIR R4 DiagnosticReport Bundle is generated containing both the imaging findings and the genomic enrichment context. A download button appears.

**Talking points:**

- "The combined imaging + genomic findings are exported as a FHIR R4 DiagnosticReport Bundle."
- "The report includes SNOMED CT codes for the lung nodule, LOINC codes for the study type, DICOM modality codes, and the genomic enrichment context."
- "This is ready to send to Epic, Cerner, or any FHIR-compliant EHR."

---

### Bridge to Stage 2: RAG Target Identification (5 minutes)

**Show:** Open Streamlit chat at http://localhost:8501

**Type this query:**
> "The imaging study identified a Lung-RADS 4A nodule. Genomic analysis found EGFR exon 19 deletion and KRAS G12C. What are the most effective targeted therapies and their resistance mechanisms?"

**Talking points:**

- "We've transitioned from imaging to the RAG pipeline. Claude is now reasoning over the same genomic vectors the cross-modal trigger queried."
- "The answer is grounded in 3.56 million variant annotations — ClinVar pathogenicity, AlphaMissense predictions, and clinical evidence."
- "Claude identifies EGFR as the primary drug target — with specific inhibitor recommendations and known resistance mechanisms like T790M."

---

### Bridge to Stage 3: Drug Discovery (4 minutes)

**Show:** Drug Discovery UI at http://localhost:8505

**Talking points:**

- "Now we take the confirmed EGFR target and feed it into Stage 3."
- "BioNeMo MolMIM generates 100 novel inhibitor analogs. DiffDock simulates binding to the EGFR kinase domain."
- "RDKit scores each candidate for drug-likeness. The top candidates outperform the seed compound."

---

### DICOM Auto-Routing: Architecture Overview (3 minutes)

**Talking points:**

- "In a live clinical environment, DICOM studies arrive in Orthanc via standard DICOM networking (C-STORE on port 4242)."
- "When a study completes, an Orthanc webhook fires automatically and the system routes the study to the correct workflow based on modality and body region."
- "CT + head gets routed to hemorrhage triage. CT + chest goes to lung nodule tracking. CR + chest gets rapid findings. MR + brain gets MS lesion tracking."
- "No manual intervention. Studies are processed the moment they arrive."

**Expected result:** Describe the auto-routing concept. In production, this is a fully automated pipeline triggered by Orthanc's Lua scripting engine and webhook callbacks.

**Routing table:**

| Modality + Region | Workflow | Target Latency |
|---|---|---|
| CT + head | ct_head_hemorrhage | < 90 sec |
| CT + chest | ct_chest_lung_nodule | < 5 min |
| CR/DX + chest | cxr_rapid_findings | < 30 sec |
| MR + brain | mri_brain_ms_lesion | < 5 min |

---

### The Complete Pipeline (3 minutes)

**Talking points:**

Walk through the full precision medicine loop on a whiteboard or slide:

```
DICOM Study (Orthanc, port 8042)
    |
Imaging Intelligence Agent (port 8524/8525)
    | Lung-RADS 4A detected
Cross-Modal Trigger --> Genomic Evidence (3.5M vectors, Milvus 19530)
    | EGFR del19 + KRAS G12C found
Stage 2: RAG Target Identification (port 8501)
    | EGFR confirmed as drug target
Stage 3: Drug Discovery (port 8505)
    | 100 novel EGFR inhibitors generated
Clinical Output
    |
FHIR R4 DiagnosticReport --> EHR Integration
PDF Report --> Clinical Documentation
```

- "From DICOM image to drug candidates. On a $3,999 desktop."
- "Every step is grounded in evidence — imaging models, genomic annotations, clinical literature."
- "This is what precision medicine looks like when you remove the barriers."

---

### Closing Route B (2 minutes)

**Talking points:**

- "You've just seen a medical image trigger genomic analysis, target identification, and drug candidate generation — all automatically."
- "This runs on a single DGX Spark. The same pipelines scale to DGX B200 for departments and DGX SuperPOD for enterprises."
- "NVIDIA FLARE enables federated learning across institutions — the models get better without sharing patient data."
- "All HCLS AI Factory code is Apache 2.0. NVIDIA NIM components are free for development on DGX Spark."

**Scaling story:**

| Phase | Hardware | Scale |
|---|---|---|
| Phase 1 | DGX Spark ($3,999) | Proof build — what you just saw |
| Phase 2 | DGX B200 | Department — multiple concurrent studies |
| Phase 3 | DGX SuperPOD | Enterprise — thousands concurrent, federated learning |

---

## Troubleshooting

### NIM Services Not Available

If NIM services show "unavailable" instead of "available" or "mock":

```bash
# Check NIM container status
docker compose ps | grep nim

# Check NIM health endpoints
curl -s http://localhost:8530/v1/health/ready  # VISTA-3D
curl -s http://localhost:8520/v1/health/ready  # Llama-3

# Restart NIM services
docker compose restart nim-vista3d nim-llm
```

NIM services require 8-16 GB GPU memory each. If GPU memory is insufficient, enable mock mode:

```bash
export IMAGING_NIM_MODE=mock
```

### Milvus Connection Issues

```bash
# Check Milvus status
curl -s http://localhost:19530/v1/vector/collections

# Check etcd and MinIO dependencies
docker compose logs milvus-etcd
docker compose logs milvus-minio

# If collections show 0 records, re-run ingestion
python3 scripts/setup_collections.py
python3 scripts/ingest_pubmed.py --max-results 5000
```

### Orthanc Not Responding

```bash
# Check Orthanc status
curl -s http://localhost:8042/system

# Check Orthanc logs
docker compose logs orthanc

# Restart Orthanc
docker compose restart orthanc
```

### Workflow Returns Empty Findings

Ensure `IMAGING_NIM_ALLOW_MOCK_FALLBACK=true` is set in your `.env` file. Without this, workflows will fail if NIM services are unavailable.

### FHIR Export Returns Errors

Verify the workflow completed successfully before requesting FHIR export. The FHIR exporter requires a valid `WorkflowResult` object with findings and measurements.

### Streamlit UI Not Loading

```bash
# Check Streamlit container status
docker compose ps | grep streamlit

# Check Streamlit logs
docker compose logs imaging-streamlit

# Restart the Streamlit service
docker compose restart imaging-streamlit
```

Verify port 8525 is not in use by another process:

```bash
lsof -i :8525
```

---

## Quick Reference

### UI Endpoints

| Service | URL |
|---|---|
| Imaging Agent UI | http://localhost:8525 |
| Imaging Agent API docs | http://localhost:8524/docs |
| Orthanc DICOM | http://localhost:8042 |
| RAG Chat (Stage 2) | http://localhost:8501 |
| Drug Discovery UI | http://localhost:8505 |
| Landing Page | http://localhost:8080 |
| Grafana Monitoring | http://localhost:3000 |

### Imaging Agent UI Tabs

| Tab | Purpose |
|---|---|
| **Ask** | Chat-based RAG Q&A with evidence sidebar |
| **Comparative** | Side-by-side entity comparison with evidence |
| **Workflow Demo** | Pre-built clinical workflow demos with Run Demo button |
| **Reports** | Export Markdown, JSON, PDF reports with download |
| **Settings** | Results per collection, NIM mode, relevance thresholds |

### API Endpoints

| Action | Endpoint |
|---|---|
| Health check | `GET http://localhost:8524/health` |
| NIM status | `GET http://localhost:8524/nim/status` |
| List collections | `GET http://localhost:8524/collections` |
| List workflows | `GET http://localhost:8524/workflows` |
| Run workflow | `POST http://localhost:8524/workflow/{name}/run` |
| RAG query | `POST http://localhost:8524/api/ask` |
| Evidence search | `POST http://localhost:8524/search` |
| VISTA-3D segment | `POST http://localhost:8524/nim/vista3d/segment` |
| Generate report | `POST http://localhost:8524/reports/generate` |
| DICOM webhook | `POST http://localhost:8524/events/dicom-webhook` |
| Event history | `GET http://localhost:8524/events/history` |
| Event status | `GET http://localhost:8524/events/status` |

---

## Appendix: API Reference

The following curl commands exercise the same functionality shown in the Streamlit UI during the live demo. Use these for scripted testing, CI/CD verification, or programmatic access.

### Health Check

```bash
curl -s http://localhost:8524/health | python3 -m json.tool
```

### NIM Service Status

```bash
curl -s http://localhost:8524/nim/status | python3 -m json.tool
```

Expected:

```json
{
  "services": [
    {"name": "vista3d", "status": "available", "url": "http://localhost:8530"},
    {"name": "maisi", "status": "mock", "url": ""},
    {"name": "vila_m3", "status": "cloud", "url": "integrate.api.nvidia.com"},
    {"name": "llm", "status": "anthropic", "url": "api.anthropic.com"}
  ],
  "available_count": 2,
  "mock_count": 1,
  "unavailable_count": 0
}
```

### List Collections

```bash
curl -s http://localhost:8524/collections | python3 -m json.tool
```

### Run Workflow: CT Head Hemorrhage

```bash
curl -s -X POST http://localhost:8524/workflow/ct_head_hemorrhage/run \
  -H "Content-Type: application/json" \
  -d '{"input_path": "", "mock_mode": true}' | python3 -m json.tool
```

### Run Workflow: CT Chest Lung Nodule

```bash
curl -s -X POST http://localhost:8524/workflow/ct_chest_lung_nodule/run \
  -H "Content-Type: application/json" \
  -d '{"input_path": "", "mock_mode": true}' | python3 -m json.tool
```

### Run Workflow: CXR Rapid Findings

```bash
curl -s -X POST http://localhost:8524/workflow/cxr_rapid_findings/run \
  -H "Content-Type: application/json" \
  -d '{"input_path": "", "mock_mode": true}' | python3 -m json.tool
```

### Run Workflow: MRI Brain MS Lesion

```bash
curl -s -X POST http://localhost:8524/workflow/mri_brain_ms_lesion/run \
  -H "Content-Type: application/json" \
  -d '{"input_path": "", "mock_mode": true}' | python3 -m json.tool
```

### RAG Knowledge Query

```bash
curl -s -X POST http://localhost:8524/api/ask \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What are the latest AI models for lung nodule detection and how do they compare?",
    "modality": "ct",
    "body_region": "chest",
    "top_k": 5
  }' | python3 -m json.tool
```

### VISTA-3D Segmentation

```bash
curl -s -X POST http://localhost:8524/nim/vista3d/segment \
  -H "Content-Type: application/json" \
  -d '{
    "input_path": "",
    "target_classes": ["liver", "spleen", "left_kidney", "right_kidney"]
  }' | python3 -m json.tool
```

Expected:

```json
{
  "classes_detected": ["liver", "spleen", "left_kidney", "right_kidney"],
  "volumes": {
    "liver": 1450.2,
    "spleen": 180.5,
    "left_kidney": 160.3,
    "right_kidney": 155.8
  },
  "num_classes": 4,
  "inference_time_ms": 45.2,
  "model": "vista3d"
}
```

### Generate Report (Markdown)

```bash
curl -s -X POST http://localhost:8524/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Summarize current guidelines for incidental lung nodule management",
    "modality": "ct",
    "body_region": "chest",
    "format": "markdown"
  }' | python3 -m json.tool
```

### Generate Report (PDF)

```bash
curl -s -X POST http://localhost:8524/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Summarize current guidelines for incidental lung nodule management",
    "format": "pdf"
  }' --output lung_nodule_report.pdf
```

### DICOM Webhook (Simulated Study Arrival)

```bash
curl -s -X POST http://localhost:8524/events/dicom-webhook \
  -H "Content-Type: application/json" \
  -d '{
    "event_type": "study.complete",
    "study_uid": "1.2.840.113619.2.55.3.604688119",
    "patient_id": "DEMO-PT-001",
    "modality": "CT",
    "body_region": "head",
    "series_count": 3,
    "instance_count": 245
  }' | python3 -m json.tool
```

Expected:

```json
{
  "study_uid": "1.2.840.113619.2.55.3.604688119",
  "patient_id": "DEMO-PT-001",
  "modality": "CT",
  "workflow_triggered": "ct_head_hemorrhage",
  "workflow_status": "completed",
  "workflow_result": {
    "findings_count": 1,
    "classification": "urgent_hemorrhage",
    "severity": "urgent",
    "inference_time_ms": 42.5,
    "nim_services_used": ["vista3d"]
  },
  "processed_at": "2026-02-28T14:30:45.123Z"
}
```

### Event History

```bash
curl -s http://localhost:8524/events/history?limit=5 | python3 -m json.tool
```

### Event Status / Routing Configuration

```bash
curl -s http://localhost:8524/events/status | python3 -m json.tool
```

---

*HCLS AI Factory — Apache 2.0 | February 2026*
