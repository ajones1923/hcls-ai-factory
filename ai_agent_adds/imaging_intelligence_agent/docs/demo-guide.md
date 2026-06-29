---
search:
  boost: 2
tags:
  - Demo
  - Walkthrough
  - Medical Imaging
  - Radiology AI
  - Clinical Decision Support
---

# Imaging Intelligence Agent -- Demo Guide

> **UI-driven walkthrough for demonstrating the Imaging Intelligence Agent on DGX Spark.**
>
> All live demo interaction uses the Streamlit Imaging Workbench -- no terminal commands during the presentation.
>
> License: Apache 2.0 | Date: March 2026

---

## Demo Overview

| Parameter | Value |
|---|---|
| **Total Duration** | 25 minutes |
| **Hardware** | NVIDIA DGX Spark (GB10, 128 GB unified) |
| **Total Vectors** | 876 seed vectors across 10 owned collections |
| **Streamlit URL** | `http://localhost:8525` |
| **API URL** | `http://localhost:8524` |
| **LLM** | Claude Sonnet 4.6 (Anthropic) |
| **Milvus** | `http://localhost:19530` |
| **Collections** | 10 owned + 1 shared `genomic_evidence` |
| **Workflows** | 6 (CT Head Hemorrhage, CT Chest Lung Nodule, CT Coronary Angiography, CXR Rapid Findings, MRI Brain MS Lesion, MRI Prostate PI-RADS) |

### What the Audience Will See

1. Four pre-loaded clinical demo cases run through imaging AI workflows with mock-realistic results
2. Cross-modal genomic enrichment linking imaging findings to genetic risk factors via the shared `genomic_evidence` collection
3. Image Gallery with AI-annotated CXR showcase, cross-modality gallery, and 3D volume slice viewer
4. Interactive Plotly network graph in Patient 360 showing gene-finding relationships
5. Multi-collection RAG queries across 10 imaging-specific knowledge bases with Claude-synthesized answers and citations, with pre-filled example queries in Evidence Explorer and Protocol Advisor
6. Automatic comparative analysis triggered by natural language -- "CT vs MRI" produces structured comparison tables
7. Patient-specific imaging protocol recommendations with AI-optimized dose reduction
8. A landscape of 50 FDA-cleared and research AI devices filterable by modality and clinical task
9. Dose intelligence showing 36% average radiation dose reduction across 20 protocols via DLIR
10. 6-step pipeline animation in Workflow Runner with annotated AI images
11. Sidebar guided tour with 9-step demo flow
12. Professional reports exported as Markdown, JSON, NVIDIA-themed PDF, and FHIR R4 DiagnosticReport Bundles

---

## Pre-Demo Setup

### Step 1: Verify Milvus Running

```bash
curl -s http://localhost:19530/v1/vector/collections
# Expected: JSON response listing imaging collections
```

### Step 2: Verify Docker Containers

```bash
docker ps | grep imaging
# Expected: imaging-intelligence-agent containers running (API + Streamlit)
```

### Step 3: Health Check

```bash
curl -s http://localhost:8524/health | python3 -m json.tool
```

Expected response:

```json
{
  "status": "healthy",
  "collections": {
    "imaging_literature": 100,
    "imaging_trials": 80,
    "imaging_findings": 100,
    "imaging_protocols": 80,
    "imaging_devices": 96,
    "imaging_anatomy": 60,
    "imaging_benchmarks": 80,
    "imaging_guidelines": 80,
    "imaging_report_templates": 100,
    "imaging_datasets": 100
  },
  "total_vectors": 876,
  "nim_services": {
    "vista3d": "mock",
    "maisi": "mock",
    "vila_m3": "mock",
    "llm": "anthropic"
  }
}
```

All 10 collections should have non-zero counts. NIM services will show `mock` unless NVIDIA API keys are configured; the LLM should show `anthropic`.

### Step 4: Verify Streamlit

```bash
curl -s http://localhost:8525/_stcore/health
# Expected: "ok"
```

### Step 5: Verify ANTHROPIC_API_KEY

```bash
echo $ANTHROPIC_API_KEY | head -c 10
# Expected: sk-ant-api (first 10 characters of your key)
```

### Step 6: Open Browser Tabs

Before starting the live demo, open the following tabs in your browser so you can switch between them without typing URLs:

| Tab | URL | Used In |
|---|---|---|
| Imaging Agent UI | http://localhost:8525 | All steps |
| Imaging Agent API docs | http://localhost:8524/docs | Reference only |
| Landing Page | http://localhost:8080 | Opening context |

---

## Route A: Guided Demo Cases (15 minutes)

> **These four cases demonstrate the full imaging AI workflow: detection, classification, measurement, and cross-modal genomic enrichment.**
>
> All cases use pre-loaded mock data that runs without live DICOM images or NIM GPU services.

### Step 1: Emergency Stroke (5 minutes)

**Click:** the **Workflow Runner** tab (the second tab) in the Streamlit UI at `http://localhost:8525`.

**Select:** "DEMO-001: Emergency Stroke: Acute Intracranial Hemorrhage" from the demo case dropdown.

**Expected result:** The clinical scenario populates:

> 62-year-old male presents to ED with sudden onset severe headache, left-sided weakness, and slurred speech. GCS 12. History of uncontrolled hypertension and type 2 diabetes. BP 195/110 on arrival. Non-contrast CT head ordered stat for acute stroke workup.

**Click:** the **Run** button.

**Expected result:** The workflow executes the `ct_head_hemorrhage` pipeline and displays:

- **Hemorrhage Detection:** Intraparenchymal hemorrhage detected
- **Location:** Right basal ganglia extending to internal capsule
- **Volume:** 28.5 mL (ABC/2 method)
- **Midline Shift:** 4.8 mm
- **Max Thickness:** 42.0 mm
- **Hounsfield Units:** Mean 65, Max 82
- **Surrounding Edema:** 12.1 mL
- **Intraventricular Extension:** Present
- **Fisher Grade:** 3
- **Severity Classification:** Critical

**Click:** the **Patient 360** tab (the eighth tab).

**Show:** Cross-modal genomic enrichment results linking the hemorrhage finding to genetic risk factors:

- **APOE** -- e4 allele associated with lobar hemorrhage risk and worse outcomes
- **COL3A1** -- Variants linked to vascular fragility
- **ACE** -- Insertion/deletion polymorphism affects cerebrovascular risk

**Talking points:**

- "AI detected the hemorrhage in under 90 seconds -- faster than radiologist page response time."
- "Automatic midline shift measurement eliminates subjective assessment. 4.8 mm shift is clinically significant."
- "Cross-modal genomic enrichment identifies genetic risk factors for hemorrhage recurrence -- APOE, COL3A1, and ACE polymorphisms."
- "Structured report with ICH Score automatically generated for neurosurgery consult."

---

### Step 2: Lung Cancer Screening (5 minutes)

**Click:** the **Workflow Runner** tab.

**Select:** "DEMO-002: Lung Cancer Screening: Suspicious Nodule Detection" from the demo case dropdown.

**Expected result:** The clinical scenario populates:

> 58-year-old female, 30 pack-year smoking history, undergoing annual low-dose CT lung cancer screening per USPSTF guidelines. Prior screening CT 12 months ago showed a 6mm ground-glass nodule in the right upper lobe (Lung-RADS 3). Follow-up scan ordered.

**Click:** the **Run** button.

**Expected result:** The workflow executes the `ct_chest_lung_nodule` pipeline and displays two nodules:

- **Nodule 1 (Primary):**
    - Location: Right upper lobe, posterior segment
    - Type: Part-solid, spiculated margin
    - Size: 18 mm (long axis) x 14 mm (short axis)
    - Solid component: 10 mm
    - Volume: 1,890 mm3 (prior: 680 mm3)
    - Volume doubling time: 245 days
    - Density: -450 HU
    - **Lung-RADS: 4B** -- Tissue sampling recommended

- **Nodule 2 (Incidental):**
    - Location: Left lower lobe, superior segment
    - Type: Ground-glass, smooth margin
    - Size: 5 mm
    - Volume: 52 mm3
    - Density: -650 HU
    - **Lung-RADS: 2** -- Benign appearance, no intervention

**Show:** Genomic enrichment panel linking to driver mutations:

- **EGFR** -- Mutations found in 15-50% of lung adenocarcinomas
- **ALK** -- Rearrangements in 3-7% of NSCLC
- **ROS1, KRAS, BRAF, MET** -- Additional targetable driver mutations

**Talking points:**

- "AI tracked nodule growth from 6 mm to 18 mm over 12 months -- volume doubling time of 245 days is concerning for malignancy."
- "Automatic Lung-RADS 4B classification triggers tissue sampling recommendation per ACR guidelines."
- "The second nodule (5 mm GGN, Lung-RADS 2) is appropriately classified as benign -- no intervention needed."
- "Cross-modal genomic query identifies targetable driver mutations: EGFR, ALK, ROS1, KRAS."
- "Complete screening-to-genomics pipeline: low-dose CT -> AI detection -> Lung-RADS classification -> molecular profiling -> treatment planning."

---

### Step 3: Cardiac Workup (5 minutes)

**Click:** the **Workflow Runner** tab.

**Select:** "DEMO-003: Cardiac Workup: Coronary Artery Disease Assessment" from the demo case dropdown.

**Expected result:** The clinical scenario populates:

> 55-year-old male presents with exertional chest pain and dyspnea on exertion for 3 weeks. Family history of premature CAD (father MI at age 50). BMI 28, total cholesterol 265 mg/dL, LDL 185 mg/dL. Stress test equivocal. Coronary CT angiography ordered for definitive evaluation.

**Click:** the **Run** button.

**Expected result:** The workflow executes the `ct_coronary_angiography` pipeline and displays:

- **Calcium Score:** 385 Agatston (92nd percentile for age/sex)
- **CAD-RADS Classification:** 4A
- **Vessel Assessment:**
    - LAD (proximal): 72% stenosis, mixed plaque -- **significant**
    - LAD (mid): 40% stenosis, calcified plaque
    - LCx (proximal): 30% stenosis, calcified plaque
    - RCA (proximal): 15% stenosis, no plaque
    - Left Main: 0% stenosis, no plaque
- **High-Risk Plaque Features:** Low-attenuation plaque, positive remodeling
- **Ejection Fraction Estimate:** 55%
- **Severity Classification:** Urgent

**Show:** Genomic enrichment panel linking to familial hypercholesterolemia genes:

- **LDLR** -- Variants cause 60-80% of familial hypercholesterolemia cases
- **PCSK9** -- Gain-of-function variants increase LDL cholesterol
- **APOB** -- Defective ligand binding causes FH
- **LPA** -- Elevated Lp(a) is an independent cardiovascular risk factor
- **9p21** -- Strongest common CAD risk variant locus

**Talking points:**

- "Calcium score of 385 (92nd percentile for age and sex) -- high atherosclerotic burden quantified in seconds."
- "AI detected 72% LAD stenosis with high-risk plaque features: low-attenuation plaque and positive remodeling."
- "CAD-RADS 4A classification automatically generated with guideline-concordant follow-up recommendations."
- "Cross-modal genomic enrichment identifies familial hypercholesterolemia genes -- LDLR, PCSK9, APOB -- family cascade screening is indicated."
- "Complete cardiac loop: CTA -> calcium score -> stenosis grading -> genomics -> risk management."

---

### Step 3b: Emergency CXR (3 minutes)

**Select:** "DEMO-004: Emergency CXR: Bilateral Pneumonia with Sepsis Risk" from the demo case dropdown.

**Expected result:** The clinical scenario populates:

> 45-year-old female presents with high fever (39.8°C), productive cough, bilateral crackles on auscultation, and hypoxia (SpO2 88% on room air). WBC 18,500, Procalcitonin 8.2 ng/mL. CXR ordered stat for pneumonia evaluation and sepsis workup.

**Click:** the **Run** button.

**Expected result:** The workflow executes the `cxr_rapid_findings` pipeline and displays:

- **Primary Finding:** Bilateral consolidation with air bronchograms
- **Distribution:** Bilateral, predominantly lower lobes
- **Pattern:** Multifocal consolidation
- **Associated Findings:** Small bilateral pleural effusions
- **Sepsis Alert:** Triggered based on clinical markers
- **Severity Classification:** Critical

**Show:** Genomic enrichment panel linking to infection susceptibility genes:

- **TLR4** -- Toll-like receptor variants affecting innate immune response
- **MBL2** -- Mannose-binding lectin deficiency increases infection risk
- **CFTR** -- Cystic fibrosis carrier status affects pulmonary defense
- **IL6, TNF** -- Cytokine storm risk variants
- **SERPINA1** -- Alpha-1 antitrypsin deficiency

**Talking points:**

- "AI detected bilateral consolidation and automatically triggered a sepsis alert based on imaging pattern plus clinical markers."
- "Cross-modal genomic enrichment identifies infection susceptibility genes -- TLR4, MBL2, and CFTR variants."
- "From CXR to sepsis risk stratification in under 30 seconds."

---

### Step 4: Image Gallery (3 minutes)

**Click:** the **Image Gallery** tab (the third tab).

**Expected result:** A hero banner displays "Medical Imaging AI Detection Showcase" with NVIDIA DGX Spark branding.

**Show:** The CXR AI Detection Showcase with 5 pathology cards:
- Normal (no findings)
- Consolidation (pneumonia)
- Pleural Effusion
- Cardiomegaly
- Pneumothorax

Each card shows the AI-annotated image with bounding boxes, severity badges, and measurement annotations. Toggle the "Show AI Annotations" switch to compare raw vs AI-analyzed views.

**Show:** The Cross-Modality AI Showcase section with CT and MRI examples alongside the CXR cases.

**Show:** The 3D Volume Slice Viewer. Select a volume (CT Head, CT Chest, or MRI Brain) and use the slice slider to scroll through the 3D volume. Adjust window center and width for different tissue contrasts.

**Talking points:**

- "These are AI-annotated medical images showing detection overlays, bounding boxes, and measurement annotations."
- "The before/after toggle demonstrates what AI adds to the raw image -- from nothing visible to precise detection with measurements."
- "The 3D slice viewer renders NIfTI volumes with adjustable HU windowing -- brain window, lung window, or custom settings."

---

## Route B: Interactive Exploration (10 minutes)

> **This route demonstrates the breadth of the knowledge base, protocol intelligence, dose optimization, and export capabilities.**

### Step 5: Evidence Explorer (3 minutes)

**Click:** the **Evidence Explorer** tab (the first tab).

**Show:** The 4 clickable pre-filled example queries: "What AI models detect hemorrhage?", "Compare CT vs MRI for stroke", "Lung nodule follow-up guidelines", "FDA-cleared AI devices for mammography". Click any to auto-populate the query box.

**Type this query** in the "Ask a question" text input:

> What AI models are used for intracranial hemorrhage detection?

**Click:** the **Ask** button.

**Expected result:**

- A markdown answer appears identifying specific models and approaches: CNN-based classifiers, U-Net segmentation, attention-based architectures, and FDA-cleared commercial solutions.
- Evidence sources display badges from multiple collections: `imaging_literature`, `imaging_devices`, `imaging_benchmarks`, `imaging_findings`.
- Each source shows its cosine similarity score and a text excerpt.
- Processing time caption appears at the bottom.

**Show:** The evidence count -- 20+ sources retrieved across multiple collections in a single query.

**Type this comparative query** in the text input:

> CT vs MRI for stroke detection

**Click:** the **Ask** button.

**Expected result:** Claude detects the comparative keywords ("vs") and produces a structured comparison table:

| Dimension | CT | MRI |
|---|---|---|
| Speed | Minutes | 30-60 minutes |
| Sensitivity (acute hemorrhage) | >95% | Lower for acute blood |
| Sensitivity (ischemic stroke) | Limited (<6h) | High (DWI within minutes) |
| Availability | Ubiquitous, 24/7 | Limited off-hours |
| Cost | Lower | Higher |
| Guideline Role | First-line for acute stroke | Follow-up and characterization |

**Talking points:**

- "One question searched all 10 collections in parallel -- literature, devices, benchmarks, and findings all contributed evidence."
- "The comparative query was auto-detected. Claude parsed 'CT' and 'MRI' as two entities and ran dual retrieval."
- "Every claim is backed by a citation with cosine similarity scores."
- "The Plotly donut chart shows collection contribution -- which knowledge bases provided evidence for this query."

---

### Step 6: Protocol Advisor (2 minutes)

**Click:** the **Protocol Advisor** tab (the fourth tab).

**Show:** The 4 clickable example indications for quick selection.

**Type this value** in the clinical indication input:

> acute chest pain, rule out pulmonary embolism

**Enter:** Patient age = `45`, Patient weight = `80` kg.

**Click:** the **Recommend Protocol** button.

**Expected result:** A protocol recommendation card displays:

- **Protocol:** CT Pulmonary Angiography (CTPA)
- **Modality:** CT
- **Contrast Agent:** Iodinated IV contrast (100 mL, bolus tracking)
- **Estimated Dose:** Optimized CTDIvol with AI dose reduction percentage
- **Duration:** Estimated scan time
- **AI Optimization Notes:** DLIR (Deep Learning Image Reconstruction) recommendations, kVp optimization for patient weight
- **Alternatives:** V/Q scan (if contrast allergy), D-dimer follow-up (if low pretest probability)

**Talking points:**

- "Patient-specific dose adjustment based on age and weight -- the protocol adapts to the individual."
- "AI optimization notes recommend DLIR for noise reduction at lower dose levels."
- "Alternatives are provided for patients with contrast allergies -- following the ALARA principle."

---

### Step 7: Dose Intelligence (2 minutes)

**Click:** the **Dose Intelligence** tab (the sixth tab).

**Expected result:** The dose comparison dashboard displays:

- **Summary Statistics:** 20 protocols analyzed, average dose reduction, max/min reduction by modality.
- **Standard vs AI-Optimized Dose Chart:** Side-by-side bar chart comparing CTDIvol (mGy) for standard protocols vs AI-optimized protocols.
- **Protocol Table:** All 20 protocols with standard dose, AI-optimized dose, reduction percentage, technique used, and image quality assessment.

**Show:** The average dose reduction metric.

**Expected result:** ~36% average dose reduction across all 20 protocols. Individual protocols range from 15-55% reduction depending on body region and technique.

**Talking points:**

- "DLIR (Deep Learning Image Reconstruction) reduces radiation dose by an average of 36% without sacrificing diagnostic image quality."
- "The highest reductions are in pediatric protocols and repeat imaging studies."
- "Image quality is maintained or improved -- AI reconstruction compensates for lower photon counts."
- "This directly supports the ALARA principle: As Low As Reasonably Achievable."

---

### Step 8: Device & AI Ecosystem (1 minute)

**Click:** the **Device & AI Ecosystem** tab (the fifth tab).

**Expected result:** A searchable, filterable catalog of AI devices displays with columns for device name, manufacturer, modality, clinical task, regulatory status, and performance metrics.

**Select:** Modality filter = `CT`, Task filter = `detection`.

**Expected result:** The table filters to show CT-specific AI detection devices -- triage tools, hemorrhage detectors, pulmonary embolism detectors, nodule detection, and fracture detection solutions.

**Show:** The total device count: 50 FDA-cleared and research AI devices across all modalities.

**Talking points:**

- "A landscape of 50+ AI devices spanning CT, MRI, X-ray, ultrasound, mammography, and PET."
- "Each device is catalogued with regulatory status, clinical task, and published performance benchmarks."
- "This collection enables evidence-based evaluation of AI tools for clinical deployment."

---

### Step 9: Reports & Export (2 minutes)

**Click:** the **Reports & Export** tab (the seventh tab).

**Click:** the **Generate Report** button to create a report from the last query.

**Expected result:** A formatted clinical report renders inline with sections:

- **Clinical Question** -- the original query
- **Analysis** -- Claude-synthesized answer with structured findings
- **Evidence Citations** -- numbered list with collection badges, document IDs, and relevance scores

**Note:** The Workflow Runner now includes a 6-step pipeline animation and displays annotated AI images alongside clinical metrics.

**Show:** The four export format buttons:

1. **Markdown** -- Copy-paste ready for clinical notes
2. **JSON** -- Structured data for programmatic consumption
3. **PDF** -- NVIDIA-themed report with green headers and professional formatting
4. **FHIR R4** -- DiagnosticReport Bundle with SNOMED CT and LOINC coding

**Click:** the **PDF** button.

**Expected result:** A PDF downloads with NVIDIA-branded formatting, clinical question, analysis, evidence citations with relevance scores, and a disclaimer footer.

**Click:** the **FHIR R4** button.

**Expected result:** A FHIR R4 DiagnosticReport Bundle renders as JSON, including `resourceType: "Bundle"`, `DiagnosticReport` resource with LOINC-coded sections, and `Observation` resources for each finding.

**Talking points:**

- "Four export formats for different integration needs: Markdown for sharing, JSON for APIs, PDF for clinical documentation, FHIR R4 for EHR interoperability."
- "The FHIR R4 bundle uses SNOMED CT and LOINC coding for standards-based interoperability."
- "PDF reports include all evidence citations with relevance scores -- fully auditable AI-assisted findings."

---

## Troubleshooting

### No Collections Found

If the health check shows zero vectors or collections are missing:

```bash
cd ai_agent_adds/imaging_intelligence_agent/agent
python3 scripts/setup_collections.py --drop-existing
python3 scripts/seed_literature.py
python3 scripts/seed_trials.py
python3 scripts/seed_findings.py
python3 scripts/seed_protocols.py
python3 scripts/seed_devices.py
python3 scripts/seed_anatomy.py
python3 scripts/seed_benchmarks.py
python3 scripts/seed_guidelines.py
python3 scripts/seed_report_templates.py
python3 scripts/seed_datasets.py
```

### RAG Returns Empty Answers

If queries return evidence but the answer is empty or generic:

```bash
# Verify ANTHROPIC_API_KEY is set and valid
echo $ANTHROPIC_API_KEY | head -c 10
# Expected: sk-ant-api

# Test Claude directly
curl -s https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "anthropic-version: 2023-06-01" \
  -H "content-type: application/json" \
  -d '{"model": "claude-sonnet-4-6", "max_tokens": 10, "messages": [{"role": "user", "content": "Hi"}]}'
```

### Workflow Failed

If a demo case workflow fails to execute:

```bash
# Check Docker containers are healthy
docker ps | grep imaging
docker compose logs imaging-api

# Verify workflow registry
curl -s http://localhost:8524/demo-cases | python3 -m json.tool
# Expected: 4 demo cases listed
```

### Port In Use

If services fail to start due to port conflicts:

```bash
# Check what is using the API port
ss -tlnp | grep 8524

# Check what is using the Streamlit port
ss -tlnp | grep 8525

# Check running Docker containers
docker ps --format "table {{.Names}}\t{{.Ports}}"
```

### Streamlit UI Not Loading

```bash
# Check Streamlit container status
docker compose ps

# View Streamlit logs
docker compose logs streamlit-imaging

# Restart Streamlit service
docker compose restart streamlit-imaging
```

### PDF Export Fails

PDF generation requires ReportLab. If PDF export returns a 501 error:

```bash
pip install reportlab
```

---

## Quick Reference

### URLs

| Resource | URL |
|---|---|
| Imaging Agent UI (Streamlit) | http://localhost:8525 |
| Imaging Agent API docs (Swagger) | http://localhost:8524/docs |
| Landing Page | http://localhost:8080 |
| Milvus | http://localhost:19530 |
| Attu (Milvus UI) | http://localhost:8000 |
| Grafana Monitoring | http://localhost:3000 |

### Collections Reference

| Collection | Purpose | Approx. Vectors |
|---|---|---|
| `imaging_literature` | Published research papers and reviews | 100 |
| `imaging_trials` | ClinicalTrials.gov AI-in-imaging records | 80 |
| `imaging_findings` | Imaging finding templates and patterns | 100 |
| `imaging_protocols` | Acquisition protocols and parameters | 80 |
| `imaging_devices` | FDA-cleared AI/ML medical devices | 96 |
| `imaging_anatomy` | Anatomical structure references | 60 |
| `imaging_benchmarks` | Model performance benchmarks | 80 |
| `imaging_guidelines` | Clinical practice guidelines (ACR, RSNA) | 80 |
| `imaging_report_templates` | Structured radiology report templates | 100 |
| `imaging_datasets` | Public imaging datasets (TCIA, PhysioNet) | 100 |
| `genomic_evidence` | Shared genomic variants (read-only, Stage 1) | 3,561,170 |

### Key curl Examples

```bash
# Health check
curl -s http://localhost:8524/health | python3 -m json.tool

# List collections
curl -s http://localhost:8524/collections | python3 -m json.tool

# List demo cases
curl -s http://localhost:8524/demo-cases | python3 -m json.tool

# Run a demo case
curl -s -X POST http://localhost:8524/demo-cases/DEMO-001/run | python3 -m json.tool

# Protocol recommendation
curl -s -X POST http://localhost:8524/protocol/recommend \
  -H "Content-Type: application/json" \
  -d '{
    "indication": "acute chest pain, rule out pulmonary embolism",
    "patient_age": 45,
    "patient_weight_kg": 80
  }' | python3 -m json.tool

# Dose summary
curl -s http://localhost:8524/dose/summary | python3 -m json.tool
```

---

## Appendix: API Reference

All endpoints below are served by the FastAPI server at `http://localhost:8524`. These are provided for programmatic integration, automated testing, and scripting -- **not for live demo use**. Use the Streamlit UI at port 8525 for demos.

Interactive API documentation is available at `http://localhost:8524/docs` (Swagger UI).

### Health and Status

```bash
# Health check -- all collections, vectors, and NIM service status
curl -s http://localhost:8524/health | python3 -m json.tool

# List all collections with vector counts and labels
curl -s http://localhost:8524/collections | python3 -m json.tool

# Knowledge graph statistics
curl -s http://localhost:8524/knowledge/stats | python3 -m json.tool

# Prometheus metrics
curl -s http://localhost:8524/metrics
```

### Demo Cases

```bash
# List all demo cases
curl -s http://localhost:8524/demo-cases | python3 -m json.tool

# Run a specific demo case (DEMO-001, DEMO-002, DEMO-003, or DEMO-004)
curl -s -X POST http://localhost:8524/demo-cases/DEMO-001/run | python3 -m json.tool

# Response includes: workflow_result, genomic_context, talking_points
```

### RAG Queries

```bash
# Full RAG query (multi-collection retrieval + LLM synthesis)
curl -s -X POST http://localhost:8524/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What AI models are used for intracranial hemorrhage detection?",
    "modality": "ct",
    "body_region": "head",
    "top_k": 5,
    "include_genomic": true
  }' | python3 -m json.tool

# Evidence-only search (retrieval only, no LLM synthesis)
curl -s -X POST http://localhost:8524/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "lung nodule detection deep learning",
    "modality": "ct",
    "top_k": 10
  }' | python3 -m json.tool

# Cross-collection entity linking
curl -s -X POST http://localhost:8524/find-related \
  -H "Content-Type: application/json" \
  -d '{
    "entity": "pulmonary embolism",
    "top_k": 5
  }' | python3 -m json.tool
```

### Protocol Optimization

```bash
# Recommend imaging protocol
curl -s -X POST http://localhost:8524/protocol/recommend \
  -H "Content-Type: application/json" \
  -d '{
    "indication": "acute chest pain, rule out pulmonary embolism",
    "patient_age": 45,
    "patient_weight_kg": 80,
    "patient_sex": "male",
    "contrast_allergy": false,
    "pregnancy": false
  }' | python3 -m json.tool
```

### Dose Intelligence

```bash
# Get all dose comparison data (20 protocols)
curl -s http://localhost:8524/dose/reference | python3 -m json.tool

# Get dose comparison for a specific protocol
curl -s http://localhost:8524/dose/comparison/chest | python3 -m json.tool

# Get summary statistics (avg/max/min reduction, by modality)
curl -s http://localhost:8524/dose/summary | python3 -m json.tool
```

### Report Export

```bash
# Markdown report
curl -s -X POST http://localhost:8524/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the optimal imaging protocol for suspected stroke?",
    "modality": "ct",
    "format": "markdown"
  }' | python3 -m json.tool

# JSON report
curl -s -X POST http://localhost:8524/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the optimal imaging protocol for suspected stroke?",
    "format": "json"
  }' | python3 -m json.tool

# PDF report (binary download)
curl -s -X POST http://localhost:8524/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the optimal imaging protocol for suspected stroke?",
    "format": "pdf"
  }' --output imaging_report.pdf
```

---

*HCLS AI Factory -- Apache 2.0 | March 2026*
