# Imaging Intelligence Agent — Architecture Design Document

**Author:** Adam Jones
**Date:** February 2026
**Version:** 1.0.0
**License:** Apache 2.0

---

## 1. Executive Summary

The Imaging Intelligence Agent extends the HCLS AI Factory platform to support automated detection, segmentation, longitudinal tracking, and clinical triage of medical imaging studies. The agent processes CT, MRI, and chest X-ray studies using NVIDIA MONAI Deploy Application Packages (MAPs) on DGX Spark hardware — from DICOM ingestion through agentic inference to structured clinical output.

Four reference workflows cover the highest-impact radiology use cases:

1. **CT Head Hemorrhage Triage** — 3D U-Net segmentation, volume estimation, midline shift measurement, urgency routing
2. **CT Chest Lung Nodule Tracking** — RetinaNet detection + SegResNet segmentation, volume doubling time, Lung-RADS classification
3. **CXR Rapid Findings** — DenseNet-121 multi-label classification with GradCAM heatmap localization
4. **MRI Brain MS Lesion Tracking** — 3D U-Net on FLAIR, SyN diffeomorphic registration, longitudinal lesion matching

The platform enables cross-modal triggers — Lung-RADS 4B+ findings can automatically trigger the Parabricks genomics pipeline for tumor profiling, and quantitative imaging endpoints feed into the Drug Discovery pipeline for treatment-response tracking.

### Key Results

| Metric | Value |
|---|---|
| Reference workflows | **4** (CT head, CT chest, CXR, MRI brain) |
| PostgreSQL + pgvector tables | **16** (studies, series, findings, measurements, embeddings, provenance, worklist_entries + views) |
| Embedding model | **BiomedCLIP-PubMedBERT** (384-dim, HNSW index) |
| LLM | **Meta-Llama3-8B-Instruct** via NIM (local, on-device) |
| Docker Compose services | **12** (Orthanc, PostgreSQL, NIM LLM, Embedding, DICOM Listener, FHIR Publisher, Agent, Portal, Prometheus, Grafana, DCGM, MAPs) |
| CT head hemorrhage latency | **< 90 seconds** end-to-end |
| CXR rapid findings latency | **< 30 seconds** end-to-end |
| CT chest lung nodule latency | **< 5 minutes** (includes prior retrieval + registration) |
| MRI brain MS lesion latency | **< 5 minutes** (includes prior retrieval + registration) |
| Clinical output formats | **3** (DICOM SR TID 1500, FHIR DiagnosticReport R4, Priority Worklist) |

---

## 2. Architecture Overview

### 2.1 Mapping to VAST AI OS

| VAST AI OS Component | Imaging Agent Role |
|---|---|
| **DataStore** | Orthanc DICOM server: DICOM studies (CT, MRI, CXR), DICOM SR, DICOM SEG |
| **DataEngine** | Event-driven pipelines: study.complete webhook → MONAI MAP inference → post-processing → persistence |
| **DataBase** | PostgreSQL + pgvector: 16 tables (structured findings + 384-dim embedding vectors), HNSW indexes |
| **InsightEngine** | BiomedCLIP embedding + pgvector similarity search + NIM LLM RAG-grounded clinical reports |
| **AgentEngine** | LangGraph reasoning agent (triage → longitudinal → population → report) + Streamlit portal |

### 2.2 System Diagram

```
                      ┌──────────────────────────────────┐
                      │     Streamlit Portal (8525)       │
                      │  Worklist | Studies | Monitoring  │
                      └──────────────┬───────────────────┘
                                     │
                      ┌──────────────▼───────────────────┐
                      │    LangGraph Agent API (8524)     │
                      │  triage → longitudinal →          │
                      │  population → report              │
                      └──────────────┬───────────────────┘
                                     │
                ┌───────────── Severity? ──────────────┐
                │ Critical/Urgent                Routine│
                ▼                                       ▼
      ┌──────────────────┐                  ┌──────────────────┐
      │ Longitudinal     │                  │ Direct Report    │
      │ Prior retrieval  │                  │ NIM LLM RAG     │
      │ Registration     │                  │ generation       │
      │ Delta analysis   │                  │                  │
      └────────┬─────────┘                  └────────┬─────────┘
               │                                     │
               └──────────────┬──────────────────────┘
                              │
          ┌───────────────────┼──────────────────────┐
          │                   │                       │
  ┌───────▼────────┐  ┌──────▼───────────┐  ┌───────▼──────────┐
  │  NIM LLM (8520)│  │ Embedding (8521) │  │ FHIR Pub. (8523) │
  │  Llama3-8B     │  │ BiomedCLIP       │  │ DiagnosticReport │
  │  Clinical RAG  │  │ 384-dim          │  │ R4 Bundle        │
  └────────────────┘  └──────────────────┘  └──────────────────┘
                              │
          ┌───────────────────┼──────────────────────┐
          │                   │                       │
  ┌───────▼────────┐  ┌──────▼───────────┐  ┌───────▼──────────┐
  │ PostgreSQL +   │  │ Orthanc (8042)   │  │ Prometheus (9099)│
  │ pgvector (5432)│  │ DICOM Server     │  │ Grafana (3000)   │
  │ 16 tables      │  │ DICOMweb + DIMSE │  │ DCGM (9400)      │
  └───────┬────────┘  └──────┬───────────┘  └──────────────────┘
          │                  │
          │           ┌──────▼───────────┐
          │           │ DICOM Listener   │
          │           │ (8522) webhook   │
          │           │ study.complete   │
          │           └──────┬───────────┘
          │                  │
          │    ┌─────────────┼──────────────────┐
          │    │             │                   │
          │  ┌─▼──────┐ ┌───▼────────┐  ┌──────▼──────┐
          │  │CT Head │ │CT Chest    │  │CXR Rapid   │
          │  │3D U-Net│ │RetinaNet + │  │DenseNet-121│
          │  │< 90s   │ │SegResNet   │  │< 30s       │
          │  └────────┘ └────────────┘  └────────────┘
          │                  │
          │           ┌──────▼───────────┐
          │           │MRI Brain MS      │
          │           │3D U-Net + SyN    │
          │           │< 5 min           │
          │           └──────────────────┘
          │
   ┌──────▼──────────────┐
   │ Nextflow DSL2       │
   │ Pipeline Orchestrator│
   │ dgx_spark profile   │
   └─────────────────────┘
```

---

## 3. Reference Workflows — Architecture Details

### 3.1 CT Head Hemorrhage Triage

| Attribute | Value |
|---|---|
| **Target latency** | < 90 seconds end-to-end |
| **Sensitivity target** | > 95% for hemorrhage > 5 mL |
| **Validation dataset** | RSNA ICH Dataset (752K slices) |
| **Model** | 3D U-Net (MONAI): spatial_dims=3, in_channels=1, out_channels=2, channels=(16,32,64,128,256), strides=(2,2,2,2), num_res_units=2, batch norm |
| **Preprocessing** | LoadImaged → EnsureChannelFirst → Orientationd (RAS) → Spacingd (1mm iso) → ScaleIntensityRanged (a_min=0, a_max=80 — CT blood window) → CropForegroundd → EnsureTyped (float32) |
| **Post-processing** | Binary segmentation → volume estimation (voxel count × voxel volume) → midline shift (center of mass vs. falx cerebri) → max thickness |
| **Urgency classification** | Brain Trauma Foundation thresholds: volume > 30 mL OR shift > 5 mm OR thickness > 10 mm → **Critical (P1)**; volume > 5 mL → **Urgent (P2)**; else **Routine (P4)** |
| **Output** | Finding + Measurements (volume_ml, shift_mm, thickness_mm) + WorklistEntry + DICOM SR (TID 1500) |

#### Pipeline Stages

```
DICOM CT Head → Orthanc (STOW-RS / C-STORE)
    │
    ├── 1. Preprocess: Reorient RAS, resample 1mm iso, blood window          [~5 sec]
    │
    ├── 2. Inference: 3D U-Net binary segmentation (hemorrhage vs. normal)    [~30 sec]
    │
    ├── 3. Volume estimation: positive voxels × voxel volume → mL            [< 1 sec]
    │
    ├── 4. Midline shift: center of mass vs. axial center × voxel spacing    [< 1 sec]
    │
    ├── 5. Urgency classification: BTF thresholds → P1/P2/P4                 [< 1 sec]
    │
    ├── 6. Persist: Finding + Measurements + WorklistEntry → PostgreSQL       [< 1 sec]
    │
    ├── 7. DICOM SR: highdicom TID 1500 measurement report → Orthanc STOW    [< 1 sec]
    │
    └── 8. Triage: P1 → Neurosurgery alert, P2 → Urgent queue                [< 1 sec]
```

**Total: < 90 seconds**

### 3.2 CT Chest Lung Nodule Tracking

| Attribute | Value |
|---|---|
| **Target latency** | < 5 minutes |
| **Detection target** | > 90% for nodules >= 4 mm |
| **Detection model** | RetinaNet (MONAI): spatial_dims=3, num_classes=1, num_anchors=6, resnet50 FPN backbone, score_thresh=0.5, nms_thresh=0.3 |
| **Segmentation model** | SegResNet (MONAI): spatial_dims=3, in_channels=1, out_channels=2, init_filters=16, blocks_down=(1,2,2,4), blocks_up=(1,1,1), batch norm |
| **Preprocessing** | Reorient RAS → resample 1mm iso → lung window (W:1500 L:-600) → normalize |
| **Longitudinal** | Prior CT chest retrieval from PostgreSQL → Orthanc WADO-RS → SyN diffeomorphic registration (ANTsPy) → volume doubling time: `VDT = (Δt × ln2) / ln(V2/V1)` |
| **Classification** | ACR Lung-RADS v2022 — rule-based: solid, ground-glass, part-solid morphology + size + growth rate → categories 1–4X |
| **Cross-modal trigger** | Lung-RADS 4B+ → POST to Nextflow API → Parabricks genomics pipeline for tumor profiling |
| **Nodule types** | Solid, ground-glass, part-solid |
| **Output** | Per-nodule Finding (type, size, volume, Lung-RADS, malignancy risk) + Measurements (long_axis_mm, short_axis_mm, volume_mm3, doubling_time_days) |

#### Pipeline Stages

```
DICOM CT Chest → Orthanc
    │
    ├── 1. Preprocess: RAS, 1mm iso, lung window                             [~5 sec]
    │
    ├── 2. Detection: RetinaNet → candidate bounding boxes                    [~30 sec]
    │
    ├── 3. Segmentation: SegResNet per-nodule within each bbox                [~60 sec]
    │
    ├── 4. Volumetrics: voxel counting × voxel spacing product                [< 1 sec]
    │
    ├── 5. Prior retrieval: Query PostgreSQL → retrieve from Orthanc          [~10 sec]
    │
    ├── 6. Registration: SyN diffeomorphic (ANTsPy) current ↔ prior          [~60 sec]
    │
    ├── 7. Volume doubling time: VDT = (Δt × ln2) / ln(V2/V1)               [< 1 sec]
    │
    ├── 8. Lung-RADS: ACR v2022 rule-based assignment                        [< 1 sec]
    │
    ├── 9. Persist + DICOM SR + FHIR                                         [< 2 sec]
    │
    └── 10. Cross-modal: If Lung-RADS 4B+ → trigger Parabricks               [< 1 sec]
```

**Total: < 5 minutes**

### 3.3 CXR Rapid Findings

| Attribute | Value |
|---|---|
| **Target latency** | < 30 seconds |
| **Sensitivity target** | > 95% for pneumothorax |
| **Validation datasets** | CheXpert + MIMIC-CXR (601K images) |
| **Model** | DenseNet-121 (MONAI): spatial_dims=2, in_channels=1 (grayscale), out_channels=5 |
| **Classification labels** | Pneumothorax, consolidation, pleural effusion, cardiomegaly, fracture |
| **Confidence thresholds** | Pneumothorax: 0.50 (lower — high-risk), consolidation: 0.60, pleural effusion: 0.55, cardiomegaly: 0.60, fracture: 0.55 |
| **Preprocessing** | LoadImaged → EnsureChannelFirst → Resized (224×224, bilinear) → ScaleIntensityd [0,1] → EnsureTyped (float32) |
| **Localization** | GradCAM heatmap (MONAI GradCAM, target layer: `class_layers.relu`) per positive finding |
| **Output** | Multi-label Findings + GradCAM images (DICOM Secondary Capture) + GSPS overlays |

#### Pipeline Stages

```
DICOM CXR → Orthanc
    │
    ├── 1. Preprocess: Resize 224×224, normalize [0,1]                        [~2 sec]
    │
    ├── 2. Inference: DenseNet-121 multi-label (5 heads)                      [~5 sec]
    │
    ├── 3. Thresholding: Per-class confidence thresholds                      [< 1 sec]
    │
    ├── 4. GradCAM: Heatmap generation per positive finding                   [~5 sec]
    │
    ├── 5. Persist: Findings → PostgreSQL                                     [< 1 sec]
    │
    ├── 6. Output: DICOM Secondary Capture (GradCAM) + GSPS                   [~2 sec]
    │
    └── 7. Triage: Pneumothorax → Emergency Medicine / Thoracic Surgery       [< 1 sec]
```

**Total: < 30 seconds**

### 3.4 MRI Brain MS Lesion Tracking

| Attribute | Value |
|---|---|
| **Target latency** | < 5 minutes |
| **Validation dataset** | ISBI MS Challenge + institutional data (1,200 MRIs) |
| **Model** | 3D U-Net (MONAI): spatial_dims=3, in_channels=1 (FLAIR), out_channels=2, channels=(32,64,128,256), strides=(2,2,2), num_res_units=2, batch norm |
| **Preprocessing** | LoadImaged → EnsureChannelFirst → Orientationd (RAS) → Spacingd (1mm iso) → NormalizeIntensityd (z-score, nonzero, channel-wise) → CropForegroundd → EnsureTyped (float32) |
| **Connected components** | scipy.ndimage.label → per-lesion volume + centroid |
| **Registration** | ANTsPy SyNRA (Rigid + Affine + SyN diffeomorphic) → warp prior lesion masks to current space |
| **Lesion matching** | Dice overlap analysis: > 0.3 overlap → matched; > 50 mm³ growth → enlarging; no match → new |
| **Disease activity** | Stable: 0 new + 0 enlarging + 0 enhancing; Active: < 3 total active; Highly Active: >= 3 active OR >= 2 enhancing |
| **Output** | MSLesionResult (total count, total volume, new/enlarging/enhancing counts, disease activity) + Measurements |

#### Pipeline Stages

```
DICOM MRI FLAIR → Orthanc
    │
    ├── 1. Preprocess: RAS, 1mm iso, z-score normalize                       [~5 sec]
    │
    ├── 2. Segmentation: 3D U-Net on FLAIR → binary lesion mask              [~30 sec]
    │
    ├── 3. Connected components: Individual lesion labeling + volume           [~5 sec]
    │
    ├── 4. Prior retrieval: Query PostgreSQL → retrieve from Orthanc          [~10 sec]
    │
    ├── 5. Registration: Affine + SyN (ANTsPy) current ↔ prior               [~90 sec]
    │
    ├── 6. Lesion matching: Overlap analysis → new/stable/enlarging           [~5 sec]
    │
    ├── 7. Disease activity: stable / active / highly_active                  [< 1 sec]
    │
    ├── 8. Persist + DICOM SR + FHIR                                         [< 2 sec]
    │
    └── 9. Routing: Highly Active → Neurology / MS Clinic alert               [< 1 sec]
```

**Total: < 5 minutes**

---

## 4. Database Schema

### 4.1 PostgreSQL + pgvector — 16 Tables

The database uses PostgreSQL 16 with the pgvector extension for hybrid structured + vector queries.

| Table | Purpose | Key Columns |
|---|---|---|
| **studies** | DICOM study metadata | study_instance_uid, patient_id, modality (CT/MR/CR/DX), body_part, status, orthanc_id |
| **series** | DICOM series within studies | series_instance_uid, study_id (FK), modality, num_instances |
| **findings** | Clinical findings from workflows | study_id (FK), workflow, finding_type, finding_code (SNOMED CT), location, laterality, severity, confidence, details (JSONB) |
| **measurements** | Quantitative measurements | finding_id (FK), measurement_type (volume/diameter/shift/count/doubling_time), value, unit, prior_value, delta_percent |
| **embeddings** | 384-dim BiomedCLIP vectors | study_id (FK), finding_id (FK), level (study/series/lesion), embedding vector(384), HNSW index |
| **provenance** | Audit trail for inference runs | study_id (FK), workflow, model_id, model_version, model_arch, duration_ms, gpu_memory_mb, status |
| **worklist_entries** | Priority-routed triage | study_id (FK), finding_id (FK), urgency, priority (P1–P4), notification, routing, acknowledged |

### 4.2 Helper Views

| View | Purpose |
|---|---|
| **active_worklist** | Unacknowledged entries ordered by priority (P1→P4), joined with studies and findings |
| **study_summary** | Aggregated finding counts per study: total, critical, urgent, max confidence |

### 4.3 Indexes

| Index | Type | Purpose |
|---|---|---|
| `idx_embeddings_hnsw` | HNSW (vector_cosine_ops, m=16, ef_construction=64) | Fast approximate nearest neighbor embedding search |
| `idx_findings_details` | GIN (JSONB) | Query workflow-specific structured data (e.g., `details->>'lung_rads'`) |
| `idx_studies_patient` | B-tree | Patient lookup for longitudinal prior retrieval |
| `idx_studies_modality` | B-tree | Modality-based filtering |
| `idx_findings_severity` | B-tree | Severity-based worklist queries |
| `idx_provenance_workflow` | B-tree | Audit trail queries by workflow type |

### 4.4 Example Hybrid Queries

```sql
-- All Lung-RADS 4A+ findings with volumes
SELECT f.*, m.value AS volume_mm3, m.delta_percent
FROM findings f
JOIN measurements m ON m.finding_id = f.id AND m.measurement_type = 'volume'
WHERE f.workflow = 'ct_chest_nodule'
  AND f.details->>'lung_rads' IN ('4A', '4B', '4X')
ORDER BY m.value DESC;

-- 10 most similar CT chest studies (vector search)
SELECT s.study_instance_uid, s.patient_id, s.study_date,
       e.embedding <=> $1::vector AS distance
FROM embeddings e
JOIN studies s ON e.study_id = s.id
WHERE e.level = 'study' AND s.modality = 'CT' AND s.body_part = 'CHEST'
ORDER BY e.embedding <=> $1::vector
LIMIT 10;

-- Growing nodules AND similar phenotype (hybrid structured + vector)
WITH growing_nodules AS (
    SELECT f.study_id, f.id AS finding_id
    FROM findings f
    JOIN measurements m ON m.finding_id = f.id
    WHERE f.workflow = 'ct_chest_nodule'
      AND m.measurement_type = 'doubling_time'
      AND m.value < 400
)
SELECT s.patient_id, s.study_date,
       e.embedding <=> $1::vector AS phenotype_distance
FROM growing_nodules gn
JOIN studies s ON gn.study_id = s.id
JOIN embeddings e ON e.study_id = s.id AND e.level = 'study'
ORDER BY e.embedding <=> $1::vector
LIMIT 10;
```

---

## 5. LangGraph Agent Architecture

### 5.1 Agent State

The LangGraph agent uses a `TypedDict` state shared across all nodes:

| Field | Type | Purpose |
|---|---|---|
| messages | list (add_messages) | Conversation history |
| study_id | int | Current study ID |
| study_uid | str | DICOM Study Instance UID |
| patient_id | str | Patient identifier |
| modality | str | CT, MR, CR, DX |
| findings | list[dict] | Clinical findings from workflow |
| measurements | list[dict] | Quantitative measurements |
| worklist_entries | list[dict] | Triage entries |
| prior_studies | list[dict] | Prior study data for longitudinal comparison |
| similar_studies | list[dict] | Similar cases from embedding search |
| rag_context | str | RAG-retrieved context |
| report_text | str | Generated clinical report |
| triage_complete | bool | Triage node completed |
| longitudinal_complete | bool | Longitudinal node completed |

### 5.2 Graph Definition

```
START → triage → [severity routing] → longitudinal → population → report → END
                         │                                            ▲
                         │ (routine)                                   │
                         └────────────────────────────────────────────┘
```

| Node | Function | Description |
|---|---|---|
| **triage** | `triage_node()` | Classify findings by urgency, create worklist entries (P1–P4), route to specialists |
| **longitudinal** | `longitudinal_node()` | Retrieve prior measurements from PostgreSQL for delta analysis (only for critical/urgent) |
| **population** | `population_node()` | Find similar studies via pgvector embedding search (top-10 nearest neighbors) |
| **report** | `report_node()` | Generate evidence-grounded clinical summary via NIM LLM RAG pipeline |

### 5.3 Conditional Routing

After triage, findings are routed based on severity:

- **Critical or Urgent** → longitudinal → population → report (full analysis path)
- **Routine** → report (direct generation, skip expensive longitudinal retrieval)

### 5.4 MCP Tool Definitions

| Tool | Purpose | Query |
|---|---|---|
| `query_findings()` | Structured finding lookup by workflow + severity | PostgreSQL SELECT with optional filters |
| `search_similar_studies()` | Vector similarity search | pgvector HNSW nearest neighbor (cosine distance) |
| `get_prior_measurements()` | Longitudinal comparison data | JOIN studies → findings → measurements for same patient + modality |

### 5.5 Specialist Routing

| Finding Type | Routing |
|---|---|
| Hemorrhage | Neurosurgery |
| Nodule | Pulmonology / Interventional Radiology |
| Pneumothorax | Emergency Medicine / Thoracic Surgery |
| MS Lesion | Neurology / MS Clinic |

---

## 6. Embedding and Similarity Search

### 6.1 Embedding Model

| Attribute | Value |
|---|---|
| Model | BiomedCLIP-PubMedBERT_256-vit_base_patch16_224 |
| Source | Microsoft |
| Dimension | 384 |
| Modality | Vision (image features) |
| Normalization | L2 unit vector |
| Index | HNSW (m=16, ef_construction=64, vector_cosine_ops) |

### 6.2 Embedding Levels

| Level | Granularity | Use Case |
|---|---|---|
| **study** | Entire imaging study | Find patients with similar overall presentation |
| **series** | Individual series within a study | Match specific sequences (e.g., FLAIR, T1, lung window) |
| **lesion** | Individual finding/lesion | Find cases with morphologically similar lesions |

### 6.3 Embedding Service API

```
POST /embed
  Body: { study_id, level, image_path }
  Response: { study_id, level, embedding: [384 floats] }

GET /health
  Response: { status: "healthy" }
```

The embedding service runs as a standalone microservice on port 8521, processing DICOM or NIfTI images through BiomedCLIP and storing normalized 384-dim vectors in pgvector.

---

## 7. Clinical Output Standards

### 7.1 DICOM Structured Report (TID 1500)

Generated using highdicom for standards-compliant measurement reports:

| Component | Implementation |
|---|---|
| Template | TID 1500 (Measurement Report) |
| SR Type | Comprehensive 3D SR |
| Coding | SNOMED CT (anatomic sites), UCUM (units), LOINC (measurements), DCM (device observer) |
| Measurements | NumericMeasurement (volume in mL, shift in mm, diameter in mm) |
| Observer | Device: "ImagingIntelligenceAgent" |
| Manufacturer | "HCLS AI Factory" |
| Delivery | STOW-RS push to Orthanc → PACS |

### 7.2 FHIR DiagnosticReport (R4)

| Component | Implementation |
|---|---|
| Resource type | DiagnosticReport + Observation (per finding) |
| Bundle type | Transaction |
| Coding systems | SNOMED CT (findings), LOINC (measurements), HL7 v3 (interpretation) |
| Interpretation | Abnormal (A) for critical/urgent, Normal (N) for routine |
| Quantitative data | Observation.component with valueQuantity (UCUM units) |
| Delivery | POST to FHIR server via fhir_publisher service (port 8523) |

### 7.3 Priority Worklist

| Priority | Urgency | SLA | Routing |
|---|---|---|---|
| **P1** | Critical | Stat — immediate read | Specialist alert |
| **P2** | Urgent | < 1 hour | Urgent queue |
| **P3** | Moderate | < 4 hours | Standard queue |
| **P4** | Routine | Standard turnaround | Normal workflow |

---

## 8. NIM LLM Integration

### 8.1 Deployment

| Attribute | Value |
|---|---|
| Model | Meta-Llama3-8B-Instruct |
| Container | `nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark` (ARM64) |
| Port | 8520 (mapped from container 8000) |
| API | OpenAI-compatible `/v1/chat/completions` |
| Authentication | NGC_API_KEY from environment |
| Temperature | 0.1 (low for clinical consistency) |

### 8.2 RAG Pipeline

```
Study findings + measurements
    │
    ├── 1. Retrieve relevant guidelines from pgvector                         [< 10 ms]
    │
    ├── 2. Retrieve prior measurements for longitudinal context               [< 10 ms]
    │
    ├── 3. Retrieve similar cases for outcome context                         [< 10 ms]
    │
    ├── 4. Construct context: findings + guidelines + priors + similar         [< 1 ms]
    │
    └── 5. Generate grounded clinical summary via NIM LLM                     [~2-5 sec]
           System: "You are a radiology AI assistant. Provide structured
           clinical summaries grounded in imaging evidence and ACR guidelines."
```

### 8.3 System Prompt

The agent uses a specialized system prompt instructing the NIM LLM to:
1. **Provide structured clinical summaries** grounded in imaging evidence
2. **Cite ACR guidelines** (Lung-RADS, BI-RADS, Brain Trauma Foundation)
3. **Reference specific measurements** and confidence scores
4. **Include recommendations** based on finding severity and longitudinal changes
5. **Distinguish AI findings from clinical judgment** — all findings are "AI-assisted" pending radiologist review

---

## 9. Infrastructure

### 9.1 Technology Stack

| Component | Technology | Version/Detail |
|---|---|---|
| DICOM server | Orthanc | 24.1.2, DICOMweb + DIMSE, port 4242/8042 |
| Vector database | PostgreSQL + pgvector | pg16, HNSW indexes, port 5432 |
| Inference framework | MONAI Deploy | MAPs (Application Packages), MONAI >= 1.3.0 |
| Embedding model | BiomedCLIP-PubMedBERT | 384-dim, Microsoft, port 8521 |
| LLM | Meta-Llama3-8B-Instruct | NIM container (ARM64), OpenAI-compatible, port 8520 |
| Agent framework | LangGraph | StateGraph with conditional routing |
| LLM client | LangChain (ChatOpenAI) | Points to local NIM endpoint |
| UI framework | Streamlit | Port 8525, NVIDIA black/green theme |
| Orchestration | Nextflow DSL2 | `dgx_spark` profile, GPU-aware |
| DICOM parsing | pydicom + highdicom | SR, SEG, GSPS generation |
| FHIR | fhir.resources | Pydantic FHIR R4 models |
| Registration | ANTsPy | SyN diffeomorphic registration |
| Monitoring | Prometheus + Grafana + DCGM | GPU metrics, pipeline traces, alerts |
| Data models | Pydantic | BaseModel + Field validation + enums |
| Hardware target | NVIDIA DGX Spark | GB10 GPU, 128GB unified LPDDR5x, $3,999 |

### 9.2 Service Ports

| Port | Service |
|---|---|
| 4242 | Orthanc DIMSE (DICOM network) |
| 5432 | PostgreSQL + pgvector |
| 8042 | Orthanc DICOMweb REST |
| 8520 | NIM LLM (Llama3-8B-Instruct) |
| 8521 | Embedding Service (BiomedCLIP) |
| 8522 | DICOM Listener (webhook) |
| 8523 | FHIR Publisher |
| 8524 | Agent API (LangGraph) |
| 8525 | Streamlit Portal |
| 3000 | Grafana |
| 9099 | Prometheus |
| 9400 | DCGM Exporter |

### 9.3 Docker Compose Services (12)

| Service | Image | GPU | Depends On |
|---|---|---|---|
| orthanc | orthancteam/orthanc:24.1.2 | No | — |
| postgres | pgvector/pgvector:pg16 | No | — |
| nim-llm | nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark | Yes (1) | — |
| embedding-service | Custom (BiomedCLIP) | Yes (1) | postgres |
| dicom-listener | Custom (webhook) | No | orthanc, postgres |
| fhir-publisher | Custom (FHIR R4) | No | postgres |
| agent | Custom (LangGraph) | No | nim-llm, postgres |
| portal | Custom (Streamlit) | No | agent, orthanc, postgres |
| dcgm-exporter | nvcr.io/nvidia/k8s/dcgm-exporter:3.3.5 | Yes (all) | — |
| prometheus | prom/prometheus:v2.48.0 | No | dcgm-exporter |
| grafana | grafana/grafana:10.2.2 | No | prometheus |
| MAPs (4) | Custom per workflow | Yes (1) | orthanc, postgres |

---

## 10. MONAI Deploy Application Packages (MAPs)

### 10.1 MAP Pattern

Every workflow is packaged as a MONAI Deploy Application Package with three operators:

```
PreprocessOperator → InferenceOperator → PostprocessOperator
```

| Operator | Responsibility |
|---|---|
| **PreprocessOperator** | Load DICOM → apply MONAI transforms → output preprocessed tensor |
| **InferenceOperator** | Load TorchScript model from mounted weights → GPU inference → output predictions |
| **PostprocessOperator** | Threshold → measurements → classification → persist to PostgreSQL → generate DICOM SR |

### 10.2 I/O Conventions

| Path | Purpose |
|---|---|
| `/var/holoscan/input/` | DICOM input directory (mounted by orchestrator) |
| `/var/holoscan/output/` | Output directory (SR, SEG, measurements JSON) |
| `/models/` | Model weights directory (mounted volume) |

### 10.3 Dockerfile Pattern

All MAPs use the same ARM64-compatible base:

```dockerfile
FROM nvcr.io/nvidia/pytorch:24.01-py3    # ARM64-compatible
WORKDIR /app
RUN pip install monai-deploy-app-sdk>=0.6.0 monai>=1.3.0 pydicom>=2.4.0 highdicom>=0.22.0
COPY app.py operators.py ./
ENV MONAI_DEPLOY_MODEL_PATH=/models
ENTRYPOINT ["python", "app.py"]
```

### 10.4 Resource Requirements

```python
@resource(cpu=4, gpu=1, memory="16Gi")
class WorkflowApp(Application):
    ...
```

---

## 11. Nextflow Pipeline Orchestration

### 11.1 Pipeline Structure

```groovy
// main.nf — Nextflow DSL2 entry point
nextflow.enable.dsl = 2

include { CT_HEAD_HEMORRHAGE }   from './modules/ct_head_hemorrhage'
include { CT_CHEST_LUNG_NODULE } from './modules/ct_chest_lung_nodule'
include { CXR_RAPID_FINDINGS }   from './modules/cxr_rapid_findings'
include { MRI_BRAIN_MS_LESION }  from './modules/mri_brain_ms_lesion'

workflow {
    // Route to selected workflow via --workflow parameter
    // ct_head | ct_chest | cxr | mri_brain
}
```

### 11.2 Workflow Module Pattern (CT Head Example)

Each workflow module follows the same three-process pattern:

```
PREPROCESS(study_uid) → INFERENCE(preprocessed) → POSTPROCESS(results)
```

| Process | Label | Container | I/O |
|---|---|---|---|
| PREPROCESS | `gpu` | imaging-agent/ct-head-hemorrhage:latest | study_uid → preprocessed/ |
| INFERENCE | `gpu` | imaging-agent/ct-head-hemorrhage:latest | preprocessed → results/ |
| POSTPROCESS | `process_low` | imaging-agent/ct-head-hemorrhage:latest | results → output/ |

### 11.3 DGX Spark Profile

```groovy
profiles {
    dgx_spark {
        docker.enabled = true
        docker.runOptions = '--gpus all'
        process {
            withLabel: 'gpu' {
                containerOptions = '--gpus all'
                accelerator = 1
                memory = '64 GB'
            }
            withLabel: 'process_low' {
                cpus = 2
                memory = '4 GB'
            }
        }
        params {
            max_gpus = 1
            max_memory = '128.GB'
            max_cpus = 12
        }
    }
}
```

### 11.4 Execution

```bash
nextflow run main.nf -profile dgx_spark --workflow ct_head --study_uid <uid>
```

---

## 12. Monitoring Stack

### 12.1 Components

| Component | Purpose | Port |
|---|---|---|
| DCGM Exporter | GPU utilization, temperature, power, memory | 9400 |
| Prometheus | Metrics scraping + alert rules | 9099 |
| Grafana | Dashboard visualization | 3000 |

### 12.2 Key DCGM Metrics

| Metric | Description |
|---|---|
| `DCGM_FI_DEV_GPU_UTIL` | GPU utilization % |
| `DCGM_FI_DEV_GPU_TEMP` | GPU temperature (C) |
| `DCGM_FI_DEV_POWER_USAGE` | Power draw (W) |
| `DCGM_FI_DEV_FB_USED` | GPU memory used (MB) |
| `DCGM_FI_DEV_FB_FREE` | GPU memory free (MB) |

### 12.3 Scrape Targets

Prometheus scrapes 5 targets: itself, DCGM exporter (5s interval), Orthanc (30s), embedding service, and agent API.

### 12.4 Alert Rules

| Alert | Condition | Severity |
|---|---|---|
| HighGPUMemory | GPU memory > 90% for 5 min | Warning |
| InferenceFailureRate | Failure rate > 0.1/sec for 2 min | Critical |

---

## 13. Pydantic Data Models

### 13.1 Enums

| Enum | Values |
|---|---|
| Modality | CT, MR, CR, DX |
| BodyPart | HEAD, CHEST, BRAIN, ABDOMEN, MSK |
| Severity | critical, urgent, moderate, routine |
| Priority | P1, P2, P3, P4 |
| WorkflowType | ct_head_hemorrhage, ct_chest_nodule, cxr_findings, mri_ms_lesion |
| StudyStatus | received, processing, completed, failed |
| DiseaseActivity | stable, active, highly_active |
| LungRADS | 1, 2, 3, 4A, 4B, 4X |

### 13.2 Core Models

| Model | Purpose | Key Fields |
|---|---|---|
| Study | DICOM study metadata | study_instance_uid, patient_id, modality, body_part, status |
| Finding | Clinical finding from a workflow | study_id, workflow, finding_type, finding_code (SNOMED), severity, confidence, details (JSONB) |
| Measurement | Quantitative measurement | finding_id, measurement_type, value, unit, prior_value, delta_percent |
| EmbeddingRecord | Vector for similarity search | study_id, finding_id, level, embedding (384-dim), model_name |
| ProvenanceBundle | Audit trail | study_id, workflow, model_id, model_version, duration_ms, gpu_memory_mb |
| WorklistEntry | Triage entry | study_id, finding_id, urgency, priority, routing |

### 13.3 Workflow-Specific Models

| Model | Workflow | Key Fields |
|---|---|---|
| HemorrhageResult | CT Head | detected, hemorrhage_type (subdural/epidural/subarachnoid/intraparenchymal), volume_ml, midline_shift_mm, urgency |
| NoduleResult | CT Chest | nodule_id, nodule_type (solid/ground-glass/part-solid), long_axis_mm, volume_mm3, lung_rads, malignancy_risk, doubling_time_days |
| CXRFindingResult | CXR | finding_name, detected, confidence, gradcam_region |
| MSLesionResult | MRI Brain | total_lesion_count, total_lesion_volume_ml, new/enlarging/enhancing counts, disease_activity |

---

## 14. Orthanc DICOM Server Configuration

### 14.1 Configuration

| Parameter | Value |
|---|---|
| AE Title | IMAGING_AGENT |
| DIMSE Port | 4242 |
| HTTP Port | 8042 |
| DICOMweb | Enabled (WADO, STOW, QIDO) |
| Lua Scripts | `on-stable-study.lua` — fires study.complete webhook after 10s stability |
| PACS Modality | Configurable (AET, host, port) |

### 14.2 Event-Driven Pipeline

```
DICOM Study → Orthanc C-STORE/STOW-RS
    │
    ├── StableAge timer (10 seconds, no new instances)
    │
    ├── Lua: OnStableStudy() fires HTTP POST to dicom-listener:8000
    │
    └── DICOM Listener: Parse study metadata → insert study row → route to workflow
```

---

## 15. ARM64 Compatibility

### 15.1 DGX Spark Requirements

The NVIDIA DGX Spark uses the Grace CPU (ARM64 / aarch64). **All containers must be ARM64-compatible.**

### 15.2 Verified ARM64 Base Images

| Image | ARM64 Support |
|---|---|
| `nvcr.io/nvidia/pytorch:24.01-py3` | Yes (NGC ARM64 tags) |
| `pgvector/pgvector:pg16` | Yes (multi-arch) |
| `orthancteam/orthanc:24.1.2` | Yes (multi-arch) |
| `prom/prometheus:v2.48.0` | Yes (multi-arch) |
| `grafana/grafana:10.2.2` | Yes (multi-arch) |
| `python:3.11-slim` | Yes (multi-arch) |

### 15.3 NIM Container Variants

For DGX Spark, append `-dgx-spark` to standard NIM image tags:

```
Standard (x86_64): nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest
DGX Spark (ARM64): nvcr.io/nvidia/nim/meta-llama3-8b-instruct:latest-dgx-spark
```

---

## 16. DGX Compute Progression

| Phase | Hardware | Price | Scope |
|---|---|---|---|
| **1 — Proof Build** | DGX Spark | $3,999 | 1-2 workflows, single user |
| **2 — Departmental** | 1-2× DGX B200 | $500K–$1M | All 4 workflows, PACS integration |
| **3 — Multi-Site** | 4-8× DGX B200 + InfiniBand | $2M–$4M | FLARE federated learning |
| **4 — AI Factory** | DGX SuperPOD | $7M–$60M+ | Thousands of concurrent studies |

---

## 17. File Structure

```
hls-imaging-agent/
├── main.nf                          # Nextflow DSL2 entry point
├── nextflow.config                  # Profiles: docker, dgx_spark
├── docker-compose.yml               # 12+ services
├── docker-compose.dev.yml           # Dev overrides (mock NIM)
├── start-services.sh                # Service startup script
├── demo.sh                          # Demo launcher
├── .env.example                     # Environment variable template
├── requirements.txt                 # Python dependencies
│
├── modules/                         # Nextflow workflow modules
│   ├── ct_head_hemorrhage.nf
│   ├── ct_chest_lung_nodule.nf
│   ├── cxr_rapid_findings.nf
│   └── mri_brain_ms_lesion.nf
│
├── maps/                            # MONAI Deploy Application Packages
│   ├── ct_head_hemorrhage/          # 3D U-Net hemorrhage triage
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   ├── operators.py
│   │   └── requirements.txt
│   ├── ct_chest_lung_nodule/        # RetinaNet + SegResNet nodule tracking
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   ├── operators.py
│   │   └── requirements.txt
│   ├── cxr_rapid_findings/          # DenseNet-121 multi-label CXR
│   │   ├── Dockerfile
│   │   ├── app.py
│   │   ├── operators.py
│   │   └── requirements.txt
│   └── mri_brain_ms_lesion/         # 3D U-Net MS lesion tracking
│       ├── Dockerfile
│       ├── app.py
│       ├── operators.py
│       └── requirements.txt
│
├── agent/                           # LangGraph clinical reasoning agent
│   ├── __init__.py
│   ├── graph.py                     # StateGraph definition
│   ├── nodes.py                     # Agent node implementations
│   ├── tools.py                     # MCP tool definitions
│   ├── state.py                     # AgentState TypedDict
│   └── prompts.py                   # Agent persona system prompts
│
├── services/                        # Microservices
│   ├── dicom_listener/              # Orthanc webhook → workflow router
│   ├── fhir_publisher/              # FHIR DiagnosticReport output
│   ├── embedding_service/           # BiomedCLIP 384-dim embeddings
│   └── portal/                      # Streamlit dashboard
│
├── src/                             # Shared library code
│   ├── models.py                    # Pydantic data models
│   ├── db.py                        # PostgreSQL + pgvector client
│   ├── dicom_utils.py               # DICOM SR, GSPS, SEG helpers
│   ├── fhir_utils.py                # FHIR resource builders
│   └── config.py                    # Environment-based configuration
│
├── db/
│   ├── init.sql                     # PostgreSQL + pgvector schema (16 tables)
│   └── migrations/                  # Schema migration scripts
│
├── models/
│   ├── README.md                    # Model download instructions
│   └── download_models.sh           # MONAI Model Zoo download script
│
├── config/
│   ├── orthanc.json                 # Orthanc DICOM server config
│   ├── prometheus.yml               # Monitoring config
│   └── grafana/dashboards/          # Grafana dashboard JSON
│
├── tests/
│   ├── conftest.py                  # Shared fixtures (testcontainers)
│   ├── test_ct_head.py
│   ├── test_ct_chest.py
│   ├── test_cxr.py
│   ├── test_mri_brain.py
│   ├── test_db.py
│   ├── test_dicom_sr.py
│   └── test_fhir.py
│
├── scripts/
│   ├── download_models.sh           # MONAI Model Zoo downloads
│   └── seed_test_data.sh            # Synthetic DICOM test data
│
└── docs/
    └── diagrams/
        └── architecture.mmd         # Mermaid architecture diagram
```

---

## 18. Testing Strategy

### 18.1 Synthetic DICOM Generation

Test data is generated using pydicom with configurable parameters:
- **CT Head:** 512×512 × 64 slices, 2.5mm slice thickness, MONOCHROME2
- **CT Chest:** 512×512 × 256 slices, 1.25mm slice thickness
- **CXR:** 2048×2048 single frame, MONOCHROME2
- **MRI Brain:** 256×256 × 128 slices, FLAIR sequence

### 18.2 Test Fixtures

| Fixture | Scope | Purpose |
|---|---|---|
| `postgres_container` | session | Spin up PostgreSQL + pgvector via testcontainers, run init.sql |
| `db_conn` | function | Per-test database connection with rollback |
| `synthetic_ct_head` | function | Generate synthetic CT head DICOM in temp directory |

### 18.3 Test Coverage

| Test File | Coverage |
|---|---|
| test_ct_head.py | Urgency classification thresholds, midline shift measurement, volume estimation |
| test_ct_chest.py | Lung-RADS assignment, VDT calculation, nodule type classification |
| test_cxr.py | Multi-label thresholding, GradCAM output |
| test_mri_brain.py | Connected component analysis, lesion matching, disease activity classification |
| test_db.py | Schema creation, CRUD operations, embedding search |
| test_dicom_sr.py | TID 1500 SR generation, SNOMED/UCUM coding |
| test_fhir.py | DiagnosticReport bundle, Observation components |

---

## 19. Demo Scenarios

### 19.1 CT Head Hemorrhage Demo

**Scenario:** 68-year-old patient arrives with acute headache. CT head shows a left subdural hematoma.

1. Upload synthetic CT head DICOM to Orthanc via STOW-RS
2. Orthanc fires study.complete webhook → DICOM listener
3. CT Head Hemorrhage MAP processes: 3D U-Net detects 25 mL subdural hemorrhage
4. Midline shift: 3.2 mm rightward
5. Classification: **Urgent (P2)** — volume > 5 mL
6. DICOM SR pushed to Orthanc with volume and shift measurements
7. Worklist entry created: P2 → Neurosurgery
8. Portal shows urgent finding in worklist with NVIDIA green/orange theming

### 19.2 Lung Nodule Tracking Demo

**Scenario:** Follow-up CT chest for a 55-year-old with a known 6mm solid nodule.

1. Upload current and prior CT chest studies
2. RetinaNet detects 3 nodules; SegResNet segments each
3. Largest nodule: 8.2 mm solid, volume 287 mm³
4. Prior comparison: was 6.1 mm / 119 mm³ → VDT = 215 days
5. Lung-RADS: **4A** (8-15mm solid + VDT < 400 → upgraded from 3)
6. Worklist: P2 → Pulmonology / Interventional Radiology
7. FHIR DiagnosticReport created with SNOMED-coded findings

### 19.3 CXR Rapid Triage Demo

**Scenario:** Trauma patient with chest X-ray.

1. Upload CXR DICOM
2. DenseNet-121 classifies: Pneumothorax (0.92), Fracture (0.71)
3. GradCAM heatmaps highlight left apex (pneumothorax) and left 7th rib (fracture)
4. Worklist: P1 → Emergency Medicine / Thoracic Surgery
5. DICOM Secondary Capture with GradCAM overlay pushed to Orthanc

### 19.4 MS Lesion Tracking Demo

**Scenario:** Follow-up MRI brain FLAIR for a 32-year-old MS patient.

1. Upload current and prior MRI brain FLAIR
2. 3D U-Net segments white matter lesions: 14 total, 2,450 mm³
3. SyN registration aligns prior scan
4. Lesion matching: 11 stable, 2 enlarging, 1 new
5. Disease activity: **Highly Active** (3 active lesions)
6. Worklist: P2 → Neurology / MS Clinic
7. Clinical report generated via NIM LLM with longitudinal context

---

## 20. Implementation Status

| Phase | Status | Details |
|---|---|---|
| **Phase 1: Architecture** | Complete | All 4 workflow designs, database schema (16 tables), Pydantic models, Docker Compose (12 services), Nextflow orchestration |
| **Phase 2: CT Head + CXR** | Complete | 3D U-Net hemorrhage MAP, DenseNet-121 CXR MAP, PostgreSQL + pgvector, Orthanc integration |
| **Phase 3: CT Chest + MRI** | Complete | RetinaNet + SegResNet nodule MAP, 3D U-Net MS lesion MAP, ANTsPy registration, longitudinal tracking |
| **Phase 4: Agent + RAG** | Complete | LangGraph agent (4 nodes), NIM LLM integration, RAG pipeline, MCP tools |
| **Phase 5: Output + Portal** | Complete | DICOM SR (TID 1500), FHIR DiagnosticReport R4, Streamlit portal, monitoring stack |

### Remaining Work

| Item | Priority | Effort |
|---|---|---|
| NVIDIA FLARE federated learning integration | Medium | 2-3 days |
| Additional MONAI model zoo models (cardiac, abdominal) | Medium | 1-2 days per workflow |
| BI-RADS mammography workflow | Low | 2-3 days |
| Full PACS integration testing (dcm4chee) | Low | 1-2 days |
| Multi-GPU inference scaling (DGX B200) | Low | 2-3 days |
| Synthetic benchmark suite (1000+ studies) | Low | 1-2 days |

---

## 21. Cross-Modal Integration

This agent demonstrates the **cross-modal intelligence** capability of the HCLS AI Factory. The same DGX Spark hardware running the Imaging Agent also supports:

- **Imaging → Genomics (Parabricks):** Lung-RADS 4B+ finding triggers tumor genomic profiling via the Stage 1 pipeline
- **Imaging → Drug Discovery (BioNeMo):** Quantitative imaging endpoints (tumor volume, response metrics) feed into the Stage 3 pipeline for treatment-response tracking
- **Imaging → Clinical Reasoning (NIM LLM):** RAG-grounded clinical reports combine imaging findings with retrieved guideline context
- **Imaging → Biomarker Agent:** Genomic variants + imaging phenotypes fused for cross-modal risk stratification
- **Imaging → CAR-T Agent:** Lymphoma staging imaging supports CAR-T treatment planning

The key architectural insight: **different data modalities share the same infrastructure patterns**. PostgreSQL + pgvector replaces Milvus for this agent (relational data + vector search in one system), but the embedding → retrieve → augment → generate pattern is identical to the CAR-T Intelligence Agent. The platform scales across therapeutic areas and imaging modalities on the same $3,999 DGX Spark.

---

## 22. Credits

- **Adam Jones**
- **Apache 2.0 License**
