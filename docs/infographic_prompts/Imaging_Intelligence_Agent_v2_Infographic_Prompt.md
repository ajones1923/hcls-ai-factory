# Nano Banana Pro — Imaging Intelligence Agent on NVIDIA DGX Spark v2.0

## IMPORTANT: Read this entire prompt before generating. This is an UPDATE to the Imaging Intelligence Agent infographic (v1.0). The original showed the imaging agent as an MVP proof build with 4 pipelines (CT Head, CT Chest, CXR, MRI Brain) running on DGX Spark. v2.0 preserves that exact visual style and architecture but updates for: renamed engine context (now "Precision Intelligence Engine"), cross-agent integration with 10 peer agents, pediatric oncology imaging capabilities (brain tumor staging, trilateral retinoblastoma screening, treatment response monitoring), 9-tab Streamlit UI, and the agent's role in the 6 pediatric demo workflows. ONE canvas, landscape 16:9, same reference architecture poster density as v1.0.

---

## REFERENCE: What v1.0 looked like (preserve this style EXACTLY)

The v1.0 had these visual characteristics that MUST be maintained:
- White canvas, landscape 16:9, dense architecture poster
- "MVP Proof Build" badge top-left (now changes to "v2.0 Production")
- Title top-center: "Imaging Intelligence Agent" with subtitle
- Key/Legend box top-right (CT Head green, CT Chest blue, CXR amber, MRI Brain purple)
- LEFT COLUMN: 7 data source cards (PACS, VNA, RIS, CT Head, CT Chest, CXR, MRI Brain) with thin teal left borders
- CENTER: Large DGX Spark container with green border and 3 internal layers:
  - Layer 1 (Data Layer): DICOM Ingestion → Local NVMe → Canonical Imaging Data
  - Layer 2 (Execution Layer): Pipeline Orchestrator → 4 color-coded pipeline strips → NVIDIA AI Software Stack
  - Layer 3 (Reasoning Layer): Structured Query + Vector Search + Evidence-Grounded Reasoning + Longitudinal + Cohort
- RIGHT-CENTER: Output Encoding column (DICOM SR, GSPS, FHIR, DICOM SEG)
- FAR RIGHT: Clinical Destinations (Radiologist, PACS Worklist, EHR)
- BOTTOM STRIP: HCLS AI Factory Cross-Modal Integration (5 connected boxes)
- VERY BOTTOM: NVIDIA infrastructure bar (green background)
- Color-coded pipeline rows (green CT Head, blue CT Chest, amber CXR, purple MRI Brain)
- Layer labels as vertical text on left edge of each layer
- Metric badges as green pills throughout
- Provenance and MVP scope annotation boxes

## WHAT CHANGES FROM v1.0 TO v2.0:

1. **Badge** → "v2.0 Production" (was "MVP Proof Build")
2. **Subtitle** → "Part of the Precision Intelligence Engine — HCLS AI Factory" (was standalone)
3. **NEW: Pediatric Imaging Pipelines** added alongside existing 4 pipelines
4. **NEW: Cross-Agent Integration section** inside Reasoning Layer showing connections to 10 peer agents
5. **NEW: 9-Tab Streamlit UI** reference in output section
6. **Cross-Modal strip** → Updated with all 8 agents and pediatric oncology focus
7. **Scope annotation** → Updated from "MVP 1-4 workflows" to "Production: 4 core + 3 pediatric workflows"
8. **Clinical Destinations** → Expanded with Pediatric Oncologist and Tumor Board
9. **Metrics** → Updated with current numbers

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background. Dense reference architecture poster. IDENTICAL visual language to v1.0.

**Typography (matching v1.0 exactly):**
- Title: Large, bold, sans-serif (Inter/Helvetica), navy (#1B2333)
- Subtitle: Smaller, gray (#666666)
- Section headers: Bold, navy, thin green left-border accent
- Sub-headers: Bold, teal (#1AAFCC)
- Body text: Small (8-10pt), sans-serif, dark gray (#333333) — SHORT PHRASES, not paragraphs
- Metric badges: Bold white text on green (#76B900) or teal (#1AAFCC) pills
- Layer labels: Vertical text, white on colored backgrounds

**Color Palette (matching v1.0 exactly):**
- NVIDIA Green: #76B900 — pipeline headers, infrastructure bar, metric badges
- Navy: #1B2333 — title, dark section bars, footer
- Teal: #1AAFCC — data flow lines, sub-headers, connections
- Light Gray: #F5F5F5 — card backgrounds
- White: #FFFFFF — canvas, card interiors
- Amber: #F5A623 — CXR pipeline, clinical output badges
- Red: #DC2626 — critical/urgent indicators
- Purple: #7B2D8E — MRI Brain pipeline, clinical destination badges
- Gray: #666666 — metadata text
- Emerald: #059669 — canonical data, normal status, pediatric verified
- Light Green: #E8F5E9 — CT Head pipeline row tint
- Light Blue: #E3F2FD — CT Chest pipeline row tint
- Light Amber: #FFF8E1 — CXR pipeline row tint
- Light Purple: #F3E8FF — MRI Brain pipeline row tint
- Light Teal: #E0F7FA — Pediatric imaging row tint [NEW]

**Visual Elements (matching v1.0):**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome)
- Solid gray arrows for primary data flow
- Dashed teal arrows for cross-modal triggers
- Bold red arrows for critical alert paths
- Color-coded pipeline rows
- Metric badge pills
- NVIDIA logo in title bar and infrastructure footer

**CRITICAL TEXT RULES:**
- All text must be legible at 1920x1080
- Use 3-5 word labels, not sentences
- Metric badges instead of numbers in sentences
- Same density as v1.0 but no garbled text

---

## CANVAS STRUCTURE (matching v1.0 layout with additions)

### ━━━ BAND 1: TITLE BAR (top of canvas) ━━━

**Left side:** Rounded badge: "Imaging Intelligence Agent" (green) + "v2.0 Production" (navy) — [CHANGED from "MVP Proof Build"]

**Center:**
- Title (large, bold, navy): **Imaging Intelligence Agent**
- Subtitle line 1 (gray): "CT / MRI / X-Ray — Production on NVIDIA DGX Spark (GB10 Grace Blackwell)"
- Subtitle line 2 (teal): "Part of the Precision Intelligence Engine — HCLS AI Factory" [NEW]
- Subtitle line 3 (gray, small): "128 GB Unified Memory | 4 Core + 3 Pediatric Imaging Pipelines" [UPDATED]

**Right side — Key/Legend box (matching v1.0, expanded):**
```
Key
━━━━━━━━━━━━━━━━━━━━
● CT Head — Hemorrhage Triage
● CT Chest — Lung Nodule Tracking
● CXR — Rapid Findings
● MRI Brain — MS Lesion Tracking
● Pediatric Brain — Tumor Staging [NEW]
● Pediatric Body — ALL/Neuroblastoma [NEW]
→  Data Flow (solid)
-→ Cross-Agent Trigger (dashed)
!→ Critical Alert (red)
```

---

### ━━━ LEFT COLUMN: IMAGING DATA SOURCES (matching v1.0, with additions) ━━━

**Matching v1.0 exactly — stacked input cards, thin teal left borders:**

1. **PACS** [hospital icon] — System of Record, Viewing/Worklists
2. **VNA** [archive icon] — Vendor-Neutral Archive, DICOMweb STOW-RS
3. **RIS** [clipboard icon] — Orders/Demographics
4. **CT Head** [brain icon] — Non-contrast CT, Acute headache/trauma
5. **CT Chest** [lungs icon] — Low-dose CT, Lung nodule follow-up
6. **CXR** [x-ray icon] — Chest X-Ray, ED/Inpatient
7. **MRI Brain** [magnet icon] — FLAIR/T1/T2, MS surveillance

**NEW ADDITIONS (teal left border, light teal background):**
8. **Pediatric Brain MRI** [child+brain icon] — Posterior fossa, Medulloblastoma staging, Trilateral screening
9. **Pediatric Body CT/MRI** [child icon] — ALL staging, Neuroblastoma, Retinoblastoma orbit

Ingestion arrows flow rightward into Band 2, labeled: "DICOMweb STOW-RS" and "DIMSE C-STORE"

---

### ━━━ BAND 2: DGX SPARK COMPUTE PLATFORM (center, largest section — matching v1.0) ━━━

Enclosed in rounded container with green (#76B900) border. Header bar: "NVIDIA DGX Spark — GB10 Grace Blackwell Superchip | 128 GB Unified Memory"

#### Layer 1: DATA LAYER (matching v1.0)

**Vertical label (navy background, white text):** "① Data Layer"

```
DICOM Ingestion → Local NVMe Storage → Canonical Imaging Data
DICOMweb STOW-RS    Immutable archive     Original DICOM | Derived Artifacts
DIMSE C-STORE SCP   Derived artifacts     Semantic Embeddings | Structured Findings
Event generation    Provenance bundles    Provenance Bundles
                    [128 GB unified]
```

#### Layer 2: EXECUTION LAYER (matching v1.0, with pediatric additions)

**Vertical label (green background, white text):** "② Execution Layer"

**Pipeline Orchestrator box** (matching v1.0):
```
Pipeline Orchestrator
study.complete trigger | Multi-stage DAGs | Container execution
```

**4 Original Pipeline Strips (matching v1.0 exactly):**

**Pipeline 1 (green row #E8F5E9): CT Head — Hemorrhage Triage**
```
3D U-Net Segmentation → Volumetric Quantification → Midline Shift → Urgency Classification
[<90 sec] [>95% sensitivity (>5 mL)]
→ Critical / Urgent / Routine
```

**Pipeline 2 (blue row #E3F2FD): CT Chest — Lung Nodule Tracking**
```
RetinaNet Detection → SegResNet Segmentation → Volumetrics → Lung-RADS → Prior Match → VDT → Risk
[<5 min] [>90% detection ≥4 mm]
- → Lung-RADS 4B+ triggers Genomics Pipeline (Parabricks) [dashed arrow]
```

**Pipeline 3 (amber row #FFF8E1): CXR — Rapid Findings**
```
DenseNet-121 Classification → GradCAM Heatmap → Confidence Scoring
Pneumothorax | Consolidation | Pleural Effusion | Cardiomegaly | Fracture
[<30 sec] [>95% pneumothorax sensitivity]
```

**Pipeline 4 (purple row #F3E8FF): MRI Brain — MS Lesion Tracking**
```
3D U-Net on FLAIR → Lesion Count + Volume → Spatial Registration → New/Enlarging ID → Disease Activity
Stable / Active / Highly Active
```

**NEW Pipeline Strips (teal tinted rows #E0F7FA):** [NEW in v2]

**Pipeline 5 (teal row): Pediatric Brain — Tumor Staging & Surveillance**
```
3D Segmentation → Tumor Volume → Posterior Fossa Assessment → Spinal Drop Mets → Response (RANO)
Medulloblastoma | Ependymoma | DIPG/DMG | Craniopharyngioma
[Trilateral screening: MRI q6mo until age 5 for RB1 carriers]
- → Oncology Agent (SHH subtype) [dashed teal arrow]
- → Neurology Agent (posterior fossa syndrome risk) [dashed teal arrow]
```
Metric badges: [Pediatric Protocol] [RANO criteria]

**Pipeline 6 (teal row): Pediatric Body — Cancer Staging**
```
Lesion Detection → Volumetrics → RECIST/Lugano Response → Staging Classification
ALL mediastinal mass | Neuroblastoma adrenal | Wilms renal | Retinoblastoma orbit
- → Oncology Agent (staging context) [dashed teal arrow]
- → Clinical Trial Agent (imaging eligibility) [dashed teal arrow]
```
Metric badges: [COG staging] [Sedation-aware protocols]

**NVIDIA AI Software Stack box (matching v1.0):**
```
MONAI Deploy MAPs | MONAI Model Zoo (3D U-Net, RetinaNet, SegResNet, DenseNet-121)
NVIDIA NIM (inference microservices, versioned) | NVIDIA FLARE (federated learning)
+ Vista3D | MAISI | ViLaM3 [NEW NIMs for v2]
```

#### Layer 3: REASONING LAYER (matching v1.0, with cross-agent addition)

**Vertical label (navy background, white text):** "③ Reasoning Layer"

**Matching v1.0 boxes:**
```
Structured Query          Vector Search            Evidence-Grounded Reasoning
SQL on findings tables    Embedding similarity     RAG: evidence + grounding
Confidence filtering      384-dim cohort match     ACR guidelines retrieval
Guideline classifications Population retrieval     LLM inference via NIM
```

**Longitudinal Comparison** | **Cohort Retrieval** (smaller boxes, matching v1.0)

**NEW: Cross-Agent Integration box (teal #1AAFCC border, prominent):** [NEW in v2]
```
Cross-Agent Integration (/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Oncology Agent (:8527)     Staging requirements, response criteria
→ Clinical Trial Agent (:8538)  Trial-specific imaging requirements
→ Cardiology Agent (:8126)   Cardiac imaging for cardiotoxicity monitoring
→ Neurology Agent (:8528)    Brain/spine MRI for CNS involvement

Timeout: 30s | Graceful degradation if peer unavailable
```

---

### ━━━ RIGHT COLUMN: OUTPUT ENCODING (matching v1.0) ━━━

**Amber-bordered boxes stacked vertically:**

1. **DICOM SR (TID 1500)** — Structured findings + measurements, Confidence scores
2. **GSPS + Secondary Capture** — GradCAM heatmaps, Annotation contours
3. **FHIR DiagnosticReport (R4)** — SNOMED CT + LOINC coded, Critical alerts
4. **DICOM SEG** — Volumetric masks, Per-lesion measurements

---

### ━━━ FAR RIGHT: CLINICAL DESTINATIONS (matching v1.0, expanded) ━━━

**Purple-bordered cards:**

1. **Radiologist / Clinician** [doctor icon] — Decision support, FDA AI/ML SaMD aligned
2. **PACS Worklist** [hospital icon] — Triage P1-P4, Critical escalation
3. **EHR** [computer icon] — FHIR DiagnosticReport
4. **Pediatric Oncologist** [child icon] [NEW] — Tumor staging, Treatment response, Trilateral screening
5. **Tumor Board** [group icon] [NEW] — Cross-modal integration with Oncology, Neurology agents

**NEW: 9-Tab Streamlit UI reference** (small box below destinations):
```
Streamlit UI (:8525)
9 Tabs: Evidence | Workflows | Gallery
Protocol | Devices | Dose | Reports
Patient 360 | Benchmarks
[9-step guided demo tour]
```

---

### ━━━ BOTTOM STRIP: HCLS AI FACTORY CROSS-MODAL (updated from v1.0) ━━━

**Header (navy background, white text):** "HCLS AI Factory — Cross-Agent Integration (11 Agents on DGX Spark)"

**Updated from 5 boxes to show broader integration:**

1. **Imaging → Genomics** — Parabricks, Lung-RADS 4B+ triggers tumor profiling [same as v1.0]
2. **Imaging → Oncology** [UPDATED] — SHH subtype analysis, Staging context, Response criteria for 11 pediatric cancers
3. **Imaging → Neurology** [NEW] — Posterior fossa syndrome, Brain MRI for CNS involvement, Neurotoxicity monitoring
4. **Imaging → Cardiology** [NEW] — Cardiac imaging for anthracycline cardiotoxicity, Echo coordination
5. **Imaging → Clinical Trial** [NEW] — Trial-specific imaging requirements, COG staging protocols
6. **Imaging → Rare Disease** [NEW] — Trilateral retinoblastoma screening for RB1 carriers

Dashed teal arrows from imaging pipelines above down to each integration box.

---

### ━━━ NVIDIA INFRASTRUCTURE BAR (matching v1.0 exactly) ━━━

**Green (#76B900) background, white text, NVIDIA logo left:**

```
DGX Spark Hardware | MONAI Deploy + Model Zoo | NVIDIA NIM | GPU Compute | Storage & Memory
GB10 Grace Blackwell | MAPs + 3D U-Net/RetinaNet/SegResNet/DenseNet-121 | Vista3D/MAISI/ViLaM3 | Blackwell GPU | 128 GB LPDDR5x
Desktop AI supercomputer | + Vista3D, MAISI, ViLaM3 [NEW] | Versioned + auto-scaling | CUDA/cuDNN/TensorRT | NVMe + GPUDirect
Proof → Department → Enterprise | NVIDIA FLARE: federated learning | NIM-served LLMs | GPU-accelerated | 
```

---

### ━━━ ANNOTATIONS (matching v1.0 style) ━━━

**Provenance annotation (dashed box, bottom-right of DGX Spark section):**
```
Provenance by Default
Complete audit trail: Model ID + Version + Params + Timestamps + DICOM UIDs
FDA AI/ML SaMD aligned | Predetermined change control
Immutable in local storage
```

**Scope annotation (dashed box, top-left of DGX Spark section):** [UPDATED]
```
v2.0 Production Scope
━━━━━━━━━━━━━━━━━━━━━
Single DGX Spark (GB10)
4 core + 3 pediatric workflows validated
6 demo workflows exercised
Cross-agent integration with 4 peer agents
MONAI Deploy MAPs + 3 NVIDIA NIMs
9-tab Streamlit UI with guided tour
Path: Proof → Department → Enterprise
```

**Performance metric badges throughout:**
- "<90 sec" on CT Head
- "<5 min" on CT Chest
- "<30 sec" on CXR
- ">95% sensitivity" on CT Head and CXR
- ">90% detection ≥4 mm" on CT Chest
- "384-dim vectors" on embeddings
- "128 GB unified memory" on DGX Spark
- "3 NVIDIA NIMs" on software stack [NEW]
- "4 peer agents" on cross-agent integration [NEW]
- "9-tab UI" on Streamlit reference [NEW]

---

## WHAT THIS v2.0 MUST COMMUNICATE vs v1.0

| v1.0 Said | v2.0 Says |
|-----------|-----------|
| MVP Proof Build | v2.0 Production |
| Standalone agent | Part of Precision Intelligence Engine with 10 peer agents |
| 4 imaging pipelines | 4 core + 3 pediatric imaging pipelines |
| No cross-agent integration | /integrated-assessment with Oncology, Neurology, Cardiology, Clinical Trial |
| No pediatric imaging | Medulloblastoma staging, trilateral screening, ALL staging, treatment response |
| Generic clinical destinations | + Pediatric Oncologist, + Tumor Board |
| 5 cross-modal connections | 6 cross-modal connections (added Neurology, Cardiology, Clinical Trial, Rare Disease) |
| No UI reference | 9-tab Streamlit UI with guided demo tour |
| No NIM models | Vista3D + MAISI + ViLaM3 prominently shown |
| "1-4 reference workflows" scope | "4 core + 3 pediatric workflows, 6 demos, 4 peer agents" scope |

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **Same architecture as v1.0** — recognizably the same poster, naturally evolved
2. **Now part of a larger platform** — "Precision Intelligence Engine" subtitle, cross-agent connections
3. **Pediatric imaging is a first-class capability** — 3 new pipeline rows for pediatric brain and body imaging
4. **Cross-agent integration is real** — visible connections to Oncology, Neurology, Cardiology, Clinical Trial agents
5. **4 core imaging pipelines preserved** — CT Head, CT Chest, CXR, MRI Brain unchanged
6. **3 NVIDIA NIMs added** — Vista3D, MAISI, ViLaM3 in the software stack
7. **Pediatric-specific protocols** — trilateral screening, sedation-aware, COG staging, RANO criteria
8. **9-tab Streamlit UI** — production-ready interface, not a prototype
9. **Still runs on one DGX Spark** — same machine, expanded capability
10. **Dense enough to be a radiology department reference poster, clear enough to present to NVIDIA**
