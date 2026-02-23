# Imaging Intelligence Agent

Automated detection, segmentation, longitudinal tracking, and clinical triage of CT, MRI, and chest X-ray studies on NVIDIA DGX Spark. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

![HCLS AI Factory Imaging AI Agent on NVIDIA DGX Spark](HCLS%20AI%20Factory%20Imaging%20AI%20Agent%20on%20NVIDIA%20DGX%20Spark%20Infographic.png)

## Overview

The Imaging Intelligence Agent processes medical imaging studies using NVIDIA MONAI models on DGX Spark hardware. It automates the full pipeline from DICOM ingestion through AI inference to structured clinical output—DICOM SR, FHIR DiagnosticReport, and priority-routed worklist entries that push directly back to PACS and EHR systems.

## Four Reference Workflows

| Workflow | Modality | Target Latency | Key Metric |
|---|---|---|---|
| **Hemorrhage Triage** | CT Head | < 90 seconds | > 95% sensitivity for bleeds > 5 mL |
| **Lung Nodule Tracking** | CT Chest | < 5 minutes | > 90% detection for nodules >= 4 mm |
| **Rapid Findings** | CXR | < 30 seconds | > 95% pneumothorax sensitivity |
| **MS Lesion Tracking** | MRI Brain | < 5 minutes | Longitudinal lesion matching + disease activity |

## Architecture

```
DICOM Study Arrives (Orthanc)
    |
    v
[Workflow Router] ──── CT Head? / CT Chest? / CXR? / MRI Brain?
    |
    v
[MONAI Deploy MAP]
(3D U-Net / RetinaNet / SegResNet / DenseNet-121)
    |
    v
[Post-Processing]
Volume estimation, midline shift, Lung-RADS, GradCAM heatmaps
    |
    v
[PostgreSQL + pgvector]
Structured findings + 384-dim BiomedCLIP embeddings
    |
    v
[LangGraph Reasoning Agent + NIM LLM]
Evidence-grounded interpretation, longitudinal comparison
    |
    v
[Clinical Output]
DICOM SR (TID 1500) | FHIR DiagnosticReport R4 | Priority Worklist
```

Built on the HCLS AI Factory platform:

- **Inference:** MONAI Deploy Application Packages (MAPs) on GB10 GPU
- **Embeddings:** BiomedCLIP-PubMedBERT (384-dim)
- **Database:** PostgreSQL + pgvector (16 tables, HNSW indexes)
- **DICOM Server:** Orthanc (DICOMweb + DIMSE)
- **LLM:** Meta-Llama3-8B-Instruct via NIM
- **Orchestration:** Nextflow DSL2
- **UI:** Streamlit portal (port 8525)
- **Hardware target:** NVIDIA DGX Spark ($3,999)

## Key Capabilities

| Capability | Detail |
|---|---|
| **CT Head Hemorrhage** | 3D U-Net segmentation, volume estimation, midline shift measurement, urgency routing (Critical/Urgent/Routine) |
| **CT Chest Lung Nodule** | RetinaNet detection + SegResNet segmentation, volume doubling time, Lung-RADS 1–4B classification |
| **CXR Rapid Findings** | DenseNet-121 multi-label classification, GradCAM heatmaps, pneumothorax/effusion/consolidation/cardiomegaly |
| **MRI Brain MS Lesion** | 3D U-Net on FLAIR, SyN diffeomorphic registration, lesion matching, disease activity (Stable/Active/Highly Active) |

## Cross-Modal Integration

The Imaging Agent connects into the broader HCLS AI Factory:

- **Lung-RADS 4B+** triggers Parabricks genomics pipeline for tumor profiling
- **Imaging phenotypes** feed into the Precision Biomarker Agent for cross-modal risk stratification
- **Quantitative imaging endpoints** support Drug Discovery pipeline for treatment-response tracking

## Clinical Output Standards

| Output | Format | Usage |
|---|---|---|
| Structured Report | DICOM SR (TID 1500) | PACS viewing, radiologist review |
| Segmentation Masks | DICOM SEG | Overlay on source images |
| Heatmaps | GSPS + Secondary Capture | GradCAM localization |
| Clinical Report | FHIR DiagnosticReport R4 | EHR integration (SNOMED CT + LOINC coded) |
| Worklist | Priority-routed entries | P1 Stat → P4 Routine triage |

## Services

| Port | Service |
|---|---|
| 8520 | NIM LLM (Llama3-8B-Instruct) |
| 8521 | Embedding Service (BiomedCLIP) |
| 8522 | DICOM Listener |
| 8523 | FHIR Publisher |
| 8524 | Agent API (LangGraph) |
| 8525 | Streamlit Portal |
| 4242 | Orthanc DIMSE |
| 8042 | Orthanc REST |

## Status

- **Phase 1 (Proof Build)** — In active development. Architecture defined, all 4 reference workflow implementations complete, documentation comprehensive.

## Credits

- **Adam Jones**
- **Apache 2.0 License**
