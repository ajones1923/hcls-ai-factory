# Napkin AI Pro — Imaging Intelligence Agent MVP on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the MVP proof build running on a single NVIDIA DGX Spark. No the platform, no enterprise-scale — this is what runs on the desk during the proof build.

---

## OVERALL LAYOUT AND STYLE

Create a dense, professional technical architecture infographic in landscape orientation (16:9 aspect ratio). The visual style should feel like a polished technical white paper diagram — clean, structured, authoritative. Match the aesthetic of a polished enterprise technical white paper (bold sans-serif headers, generous use of white space within structured sections, thin-line icons, card-based layouts with subtle borders) combined with the dense technical poster feel of the "HCLS AI Factory on NVIDIA DGX Spark" reference infographic.

**Canvas:** White background (#FFFFFF). Dense but organized — every section carries information. Clean visual hierarchy with clear section boundaries. The diagram should feel like a reference architecture poster a solutions architect pins to their wall during a proof build.

**Typography:**
- Title: Large, bold, sans-serif (Inter, Helvetica, or similar), deep navy (#1B2333)
- Subtitle: Smaller, medium gray (#666666), directly below title
- Section headers: Bold, dark navy (#1B2333) on white, with a thin NVIDIA green (#76B900) left-border accent or underline
- Sub-headers: Bold, teal (#1AAFCC)
- Body text: Small (8-10pt equivalent), clean sans-serif, dark gray (#333333)
- Metric callouts: Bold, slightly larger than body, inside small rounded green (#76B900) or teal (#1AAFCC) pill badges with white text

**Color Palette (exact):**
- NVIDIA Green: #76B900 — primary accent for all NVIDIA components, pipeline headers, infrastructure bar, metric badges
- Deep Navy: #1B2333 — title text, dark section bars, footer
- Teal: #1AAFCC — secondary accent for data flow lines, sub-headers, connection lines
- Light Gray: #F5F5F5 — card backgrounds, pipeline row backgrounds
- White: #FFFFFF — canvas, text on dark backgrounds, card interiors
- Amber/Orange: #F5A623 — clinical output badges, elevated/watch indicators
- Red: #DC2626 — critical/urgent finding indicators, bold alert arrows
- Purple: #7B2D8E — clinical destination badges (EHR, Clinician, PACS)
- Medium Gray: #666666 — metadata text, secondary labels
- Emerald Green: #059669 — canonical data / normal status indicators

**Visual Elements:**
- Rounded-corner rectangles (8px radius) for every component/process box
- Thin-line icons (16x16 to 24x24) next to data sources and outputs — simple, monochrome line icons in a clean enterprise white paper style (not emoji)
- Directional arrows: solid medium gray (#999999) for primary data flow, dashed teal (#1AAFCC) for cross-modal triggers, bold red (#DC2626) for critical alert paths
- Color-coded pipeline rows with distinct light background tints
- Metric badges: small rounded pills with white text on green or teal background
- NVIDIA logo mark in the title bar and infrastructure footer
- No the platform logo, no the platform branding, no the platform component names anywhere

---

## CANVAS STRUCTURE (Top to Bottom, 7 horizontal bands)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Small rounded badge in NVIDIA green (#76B900) with white text: "Imaging Intelligence Agent" — plus a second smaller badge below it: "MVP Proof Build" in navy (#1B2333) with white text

**Center (the dominant visual element of the title bar — large, centered, commanding):**
- **Title (large, bold, navy #1B2333, centered):** "Imaging Intelligence Agent"
- **Subtitle line 1 (medium, gray #666666, centered):** "CT / MRI / X-Ray — MVP on NVIDIA DGX Spark"
- **Subtitle line 2 (smaller, gray #666666, centered):** "GB10 Grace Blackwell Superchip | 128 GB Unified Memory"
- **Date / Author line (smallest, gray #999999, centered):** "February 2026 | HCLS AI Factory"

**Right side — Key/Legend box** (small, top-right corner, thin gray border, white background):
```
Key
——————————————————
● CT Head — Hemorrhage Triage Pipeline
● CT Chest — Lung Nodule Tracking Pipeline
● CXR — Rapid Findings Pipeline
● MRI Brain — MS Lesion Tracking Pipeline
→ Data Flow (solid)
- → Cross-Modal Trigger (dashed)
```
Use small colored dots: green for CT Head, blue for CT Chest, amber for CXR, purple for MRI Brain. These four colors are used throughout the diagram to identify each pipeline.

---

### ━━━ BAND 2: IMAGING DATA SOURCES (left column, spanning vertically alongside Bands 3-4) ━━━

**Position:** Left edge of canvas, vertically stacked column of input cards. Each card has a thin-line icon, bold label, and 1-2 lines of detail. White background (#FFFFFF) with thin teal (#1AAFCC) left-border accent.

**Cards (top to bottom):**

1. **PACS** [hospital icon]
   System of Record
   Viewing / Worklists

2. **VNA** [archive icon]
   Vendor-Neutral Archive
   DICOMweb STOW-RS

3. **RIS** [clipboard icon]
   Radiology Information System
   Orders / Demographics

4. **CT Head** [brain icon]
   Non-contrast CT
   Acute headache / trauma

5. **CT Chest** [lungs icon]
   Low-dose CT
   Lung nodule follow-up

6. **CXR** [x-ray icon]
   Chest X-Ray
   ED / Inpatient

7. **MRI Brain** [magnet icon]
   FLAIR / T1 / T2
   MS surveillance

**Ingestion arrows** flow rightward from each card into Band 3, labeled:
- "DICOMweb STOW-RS" (from PACS/VNA)
- "DIMSE C-STORE" (alternative path)
- "Event-triggered pull" (from PACS/VNA)

---

### ━━━ BAND 3: DGX SPARK COMPUTE PLATFORM (center of canvas, largest section) ━━━

This is the core of the diagram. It occupies roughly 55-60% of the canvas width (center) and is organized as a layered architecture within the DGX Spark. The entire section is enclosed in a subtle rounded-corner container with a thin NVIDIA green border and a header bar.

**Container header bar (NVIDIA green #76B900 background, white bold text):**
"NVIDIA DGX Spark — GB10 Grace Blackwell Superchip | 128 GB Unified Memory"

#### ── Layer 3A: DATA INGESTION & STORAGE ──

**Left label (vertical text, navy background #1B2333, white text):** "① Data Layer"

**Contents (left to right):**

**Box: DICOM Ingestion**
- Thin green-bordered rounded rectangle
- Header: "DICOM Ingestion" in bold
- Body text:
  - DICOMweb STOW-RS receiver
  - DIMSE C-STORE SCP
  - Study-level event generation
  - DICOM validation and routing

**Box: Local NVMe Storage**
- Thin green-bordered rounded rectangle
- Header: "Local NVMe Storage" in bold
- Body text:
  - Immutable DICOM archive
  - Derived artifacts (masks, measurements, heatmaps)
  - Provenance bundles (model ID, version, params, timestamps)
  - Canonical data model
- Small badge: "128 GB unified memory" in green pill

**Box: Canonical Imaging Data (larger, prominent, emerald-bordered #059669)**
- Visually distinct — slightly larger, emerald green (#059669) border, light green (#D1FAE5) background tint
- Header: "Canonical Imaging Data" in bold
- 5 sub-items in a row of small cards:

  | Original DICOM | Derived Artifacts | Semantic Embeddings | Structured Findings | Provenance Bundles |
  |---|---|---|---|---|
  | CT / MRI / X-Ray | Segmentation masks | 384-dim vectors | DICOM SR compatible | Model ID + version |
  | Immutable evidence | Measurements | Study / Series / Lesion | Confidence scores | Inference params |
  | STOW-RS / C-STORE | GradCAM heatmaps | Cohort matching | Lung-RADS mappings | DICOM UID lineage |
  | | Spatial registrations | | | Timestamps |

#### ── Layer 3B: INFERENCE EXECUTION ──

**Left label (vertical text, NVIDIA green #76B900 background, white text):** "② Execution Layer"

**Contents (left to right):**

**Box: Pipeline Orchestrator**
- Green-bordered rounded rectangle
- Header: "Pipeline Orchestrator" in bold
- Body text:
  - study.complete event trigger
  - Multi-stage pipeline DAGs
  - Container-based execution
  - Automatic output persistence

**Then, four color-coded pipeline flow strips running left-to-right:**

**Pipeline 1 (green-tinted row #E8F5E9): CT Head — Hemorrhage Triage**
- Flow: `3D U-Net Segmentation` → `Volumetric Quantification` → `Midline Shift` → `Urgency Classification`
- Metric badges: "< 90 sec target" | "> 95% sensitivity (> 5 mL)"
- Output arrow: "Critical / Urgent / Routine"

**Pipeline 2 (blue-tinted row #E3F2FD): CT Chest — Lung Nodule Tracking**
- Flow: `RetinaNet Detection` → `SegResNet Segmentation` → `Volumetrics` → `Lung-RADS Assignment` → `Prior Matching` → `VDT Calculation` → `Risk Scoring`
- Metric badges: "< 5 min target" | "> 90% detection (≥ 4 mm)"
- Dashed arrow from "Lung-RADS 4B+" to an external trigger box labeled "→ Genomics Pipeline (Parabricks)"

**Pipeline 3 (amber-tinted row #FFF8E1): CXR — Rapid Findings**
- Flow: `DenseNet-121 Classification` → `GradCAM Heatmap` → `Confidence Scoring`
- Below Classification: Pneumothorax | Consolidation | Pleural Effusion | Cardiomegaly | Fracture
- Metric badges: "< 30 sec target" | "> 95% pneumothorax sensitivity"

**Pipeline 4 (purple-tinted row #F3E8FF): MRI Brain — MS Lesion Tracking**
- Flow: `3D U-Net on FLAIR` → `Lesion Count + Volume` → `Spatial Registration` → `New/Enlarging ID` → `Disease Activity`
- Below Disease Activity: "Stable / Active / Highly Active"

**All four pipelines output arrows flow rightward** into Band 4.

**Box: NVIDIA AI Software Stack (below or alongside pipelines)**
- Green-bordered, NVIDIA green header accent
- Sub-items as small badges in a row:
  - `MONAI Deploy MAPs` — Containerized inference packages
  - `MONAI Model Zoo` — 3D U-Net | RetinaNet | SegResNet | DenseNet-121
  - `NVIDIA NIM` — Inference microservices | Versioned | Auto-scaling
  - `NVIDIA FLARE` — Federated learning across sites
- Arrow labeled "Orchestrated by Pipeline Orchestrator"

#### ── Layer 3C: REASONING & QUERY ──

**Left label (vertical text, navy #1B2333 background, white text):** "③ Reasoning Layer"

**Contents (left to right):**

**Box: Structured Query**
- Teal-bordered rounded rectangle
- Header: "Structured Query" in bold
- Body text:
  - SQL queries on findings tables
  - "All Lung-RADS 4A+ patients"
  - Confidence score filtering
  - Guideline-mapped classifications

**Box: Vector Search**
- Teal-bordered rounded rectangle
- Header: "Vector Search" in bold
- Body text:
  - Embedding similarity (384-dim)
  - "10 most similar CT studies"
  - Cohort matching
  - Population-scale retrieval

**Box: Evidence-Grounded Reasoning**
- Teal-bordered rounded rectangle
- Header: "Evidence-Grounded Reasoning" in bold
- Body text:
  - RAG: Evidence retrieval + grounding
  - Longitudinal delta analysis
  - Cross-modal enrichment
  - ACR guidelines retrieval
  - LLM inference via NIM
- Bidirectional arrow with Structured Query and Vector Search
- Upward arrow into pipelines labeled "Evidence + Priors + Cohort context"

**Box: Longitudinal Comparison**
- Smaller teal-bordered box
- Prior study retrieval | Delta analysis | Volume doubling time

**Box: Cohort Retrieval**
- Smaller teal-bordered box
- Embedding similarity | Patients-like-this | Outcomes matching

---

### ━━━ BAND 4: OUTPUT ENCODING (right-center column) ━━━

**Position:** Right of the DGX Spark container, receiving arrows from all four pipelines.

**Four output boxes stacked vertically, amber/orange-bordered (#F5A623), white background:**

1. **DICOM SR (TID 1500)**
   Structured findings + measurements
   Confidence scores + classifications

2. **GSPS + Secondary Capture**
   Graphic annotation contours (GSPS)
   GradCAM heatmap images (Secondary Capture)
   Pushed via DICOMweb STOW-RS

3. **FHIR DiagnosticReport (R4)**
   SNOMED CT + LOINC coded
   Narrative summaries
   Critical finding alerts

4. **DICOM SEG**
   Volumetric segmentation masks
   Per-lesion measurements

---

### ━━━ BAND 5: CLINICAL DESTINATIONS (right edge of canvas) ━━━

**Position:** Far right column, receiving arrows from Band 4 outputs.

**Three destination cards, purple-bordered (#7B2D8E), with thin-line icons:**

1. **Radiologist / Clinician** [doctor icon]
   - Clinician-in-the-Loop
   - Decision support, not autonomous Dx
   - FDA AI/ML SaMD aligned
   - Arrow from DICOM SR + GSPS + Alerts

2. **PACS Worklist** [hospital icon]
   - Triage: P1 Stat → P4 Routine
   - On-call notification
   - Critical finding escalation
   - Arrow from DICOM SR + GSPS

3. **EHR** [computer icon]
   - FHIR DiagnosticReport
   - Structured clinical summaries
   - Arrow from FHIR output

---

### ━━━ BAND 6: HCLS AI FACTORY CROSS-MODAL (bottom strip, above infrastructure) ━━━

**Position:** Horizontal strip below the DGX Spark container. Light indigo background (#E0E7FF). This shows the Imaging Agent as one node in the broader HCLS AI Factory — all running on the same DGX Spark during MVP.

**Header bar (navy #1B2333 background, white text):** "HCLS AI Factory — Cross-Modal Integration (MVP on DGX Spark)"

**Five connected boxes in a horizontal row, with dashed arrows from the imaging pipelines above:**

1. **Imaging → Genomics**
   Parabricks
   Lung-RADS 4B+ triggers tumor profiling
   Somatic + germline variant calling
   *(dashed arrow from CT Chest pipeline)*

2. **Imaging → RAG / Clinical Chat**
   NIM-served LLMs
   Evidence retrieval + guideline grounding
   Cross-modal reasoning
   *(arrow from Evidence-Grounded Reasoning)*

3. **Imaging → Drug Discovery**
   BioNeMo
   Tumor volume endpoints
   Treatment stratification
   *(dashed arrow from Genomics box)*

4. **Imaging → Longitudinal Care**
   Continuous monitoring
   Automated change detection
   Population-scale cohort analysis
   *(arrow from Longitudinal Comparison)*

5. **Imaging → Outcomes**
   Cohort retrieval (patients-like-this)
   Imaging trajectories → outcomes
   Care pathway optimization
   *(arrow from Cohort Retrieval)*

---

### ━━━ BAND 7: NVIDIA DGX SPARK INFRASTRUCTURE BAR (very bottom of canvas) ━━━

**Position:** Full-width horizontal bar at the very bottom. NVIDIA green (#76B900) background. White text throughout. NVIDIA logo mark on left side.

**Single row of hardware + software components:**

| DGX Spark Hardware | MONAI Deploy + Model Zoo | NVIDIA NIM | GPU Compute | Storage & Memory |
|---|---|---|---|---|
| GB10 Grace Blackwell Superchip | MAPs: containerized inference | Inference microservices | Blackwell GPU architecture | 128 GB unified LPDDR5x |
| Desktop AI supercomputer | 3D U-Net / RetinaNet / SegResNet / DenseNet-121 | Versioned + auto-scaling | CUDA / cuDNN / TensorRT | NVMe local storage |
| Proof → Department → Enterprise | NVIDIA FLARE: federated learning | NIM-served LLMs | GPU-accelerated inference | GPUDirect memory access |

---

## ADDITIONAL DETAIL AND ANNOTATIONS

**Performance metric badges** scattered throughout (small rounded green pills with white text):
- "< 90 sec" on CT Head pipeline
- "< 5 min" on CT Chest pipeline
- "< 30 sec" on CXR pipeline
- "> 95% sensitivity" on CT Head and CXR
- "> 90% detection ≥ 4 mm" on CT Chest
- "384-dim vectors" on Semantic Embeddings
- "128 GB unified memory" on DGX Spark
- All latency metrics labeled as "target"

**Provenance annotation** (small dashed box, bottom-right of the DGX Spark section):
- "Provenance by Default"
- Complete audit trail: Model ID + Version + Params + Timestamps + DICOM UIDs
- FDA AI/ML SaMD aligned | Predetermined change control
- Immutable in local storage

**MVP scope annotation** (small dashed box, top-left of DGX Spark section):
- "MVP Proof Build Scope"
- Single DGX Spark (GB10 Grace Blackwell Superchip)
- 1-4 reference workflows validated
- MONAI Deploy MAPs on local compute
- Canonical data model proven
- Path: Proof → Department (DGX cluster) → Enterprise (DGX SuperPOD)

**Data flow arrows style guide:**
- **Solid arrows:** Primary data flow (ingestion → processing → output). Medium gray (#999999) with arrowheads.
- **Dashed arrows:** Cross-modal triggers and enrichment paths. Teal (#1AAFCC) dashes.
- **Bold arrows:** Critical finding alert paths. Red (#DC2626).

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **This runs on a single NVIDIA DGX Spark** — a desktop AI supercomputer, not a data center
2. **Four imaging modalities** (CT Head, CT Chest, CXR, MRI Brain) enter from the left
3. **Three-layer architecture** (Data / Execution / Reasoning) runs entirely on DGX Spark
4. **Four distinct pipelines** execute with specific MONAI models and target performance metrics
5. **Clinical outputs** (DICOM SR, GSPS, FHIR, DICOM SEG) flow to the right
6. **Clinical consumers** (Radiologist, PACS Worklist, EHR) receive results on the far right
7. **HCLS AI Factory cross-modal links** show imaging connecting to genomics, drug discovery, and outcomes
8. **NVIDIA infrastructure** anchors everything at the bottom — this is built entirely on NVIDIA silicon and software
9. **No the platform references anywhere** — this is the DGX Spark proof build, pure NVIDIA stack
10. **Clear path from MVP to enterprise** — annotation shows Proof → Department → Enterprise scaling trajectory

The overall impression should be: a complete, technically precise, production-quality imaging AI agent running on a single DGX Spark desktop supercomputer — dense enough to be a reference poster, clean enough to present to NVIDIA executives.
