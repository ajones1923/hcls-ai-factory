# Clinical Imaging Engine

![Clinical Imaging Engine Architecture Infographic](infographic.jpg)

**Engine 4 of the HCLS AI Factory** | **Apache 2.0** | **NVIDIA DGX Spark ($4,699)**

> The most comprehensive open-source medical imaging AI platform available. 9 clinical workflows, 7 standardized scoring systems, 20 NVIDIA technologies, cross-modal genomic integration — all on a single device for $4,699.

## Overview

The Clinical Imaging Engine (Engine 4) processes medical imaging studies across CT, MRI, X-ray, mammography, and ultrasound using 20 NVIDIA technologies on DGX Spark hardware. It delivers:

- **9 clinical workflows** with 7 standardized scoring systems (Lung-RADS, BI-RADS, TI-RADS, LI-RADS, CAD-RADS, PI-RADS, ASPECTS)
- **38,028 indexed vectors** across 13 Milvus collections, including 1,938 real PubMed research papers
- **8 cross-modal genomic triggers** that automatically connect imaging findings to 35,678 genomic variant vectors
- **9 NIM clients** with 3-tier fallback (local → cloud → mock)
- **~1,500 radiomics features** per segmented region via PyRadiomics-CUDA
- **Agentic reasoning** (Plan/Execute/Reflect/Refine) via NVIDIA AIQ Toolkit
- **Clinical safety guardrails** via NeMo Guardrails (PII detection, evidence grounding, disclaimer injection)
- **5 export formats** (Markdown, JSON, PDF, FHIR R4 with 103 SNOMED codes, DICOM SR)
- **Real-time streaming** for ultrasound and endoscopy via NVIDIA Holoscan (30fps)
- **9 MONAI Deploy Application Packages** (MAPs) for portable clinical AI deployment
- **Interactive annotation** via MONAI Label with NVIDIA FLARE federated learning bridge
- **React portal** (clinical-grade UI) + Streamlit (developer workbench)
- **1,324 tests passing** with comprehensive mock mode for GPU-free development
- **3-tier deployment:** Community (free), Enterprise (+ AI Enterprise NIMs), Research (+ noncommercial models)

## Nine Clinical Workflows

| Workflow | Modality | Scoring System | Key Output | Cross-Modal Trigger |
|---|---|---|---|---|
| **CT Head Hemorrhage** | CT | ASPECTS / BTF | Volume, midline shift, urgency | → APOE, COL3A1, ACE |
| **CT Chest Lung Nodule** | CT | Lung-RADS v2022 | Nodule size, VDT, classification | → EGFR, ALK, ROS1, KRAS |
| **CT Coronary Angiography** | CT | CAD-RADS 2.0 | Stenosis %, calcium score, plaque | → LDLR, PCSK9, APOB |
| **CXR Rapid Findings** | CXR | Multi-label | Consolidation, effusion, pneumothorax | → TLR4, MBL2 |
| **MRI Brain MS Lesion** | MRI | MS Activity | Lesion count, volume, activity | → HLA-DRB1, IL7R |
| **MRI Prostate PI-RADS** | MRI | PI-RADS v2.1 | Lesion PI-RADS score, ADC, zone | → BRCA2, HOXB13, ATM |
| **Breast BI-RADS** | Mammography | BI-RADS 5th Ed | Mass/calcification, BI-RADS 0-6 | → BRCA1, BRCA2, PALB2 |
| **Thyroid TI-RADS** | Ultrasound | ACR TI-RADS | Nodule TI-RADS TR1-5, FNA rec | → BRAF, RAS, RET/PTC |
| **Liver LI-RADS** | CT/MRI | LI-RADS v2018 | APHE, washout, capsule, LR-1 to LR-5 | → TP53, CTNNB1, TERT |

## Architecture

```
DICOM Study Arrives (Orthanc 8042/4242)
    |
    v
[Webhook Router] ── 9 workflow routing rules by modality + body region
    |
    v
[Clinical Workflow] (9 workflows with standardized scoring)
(SegResNet / RetinaNet / DenseNet-121 / UNEST / BI-RADS / TI-RADS / LI-RADS)
    |
    v
[Post-Processing + Cross-Modal Trigger]
8 triggers: Lung-RADS 4A+ → EGFR | BI-RADS 4+ → BRCA | TI-RADS 4+ → BRAF | LI-RADS 4+ → TP53 | ...
    |
    v
[RAG Engine + Agentic Reasoning (AIQ)]
13 collections (38,028 vectors) + genomic_evidence (35,678 vectors)
Plan → Execute → Reflect → Refine cycle with 6 registered tools
    |
    v
[Safety Guardrails (NeMo Guardrails)]
PII detection | Evidence grounding | Clinical disclaimer | Contraindication check
    |
    v
[Clinical Output — 5 Formats]
Markdown | JSON | PDF | FHIR R4 (103 SNOMED codes) | DICOM SR (TID 1500)
    |
    v
[React Portal (8550) + Streamlit (8525)]
10 pages: Dashboard, Workflows, Evidence, Protocol, Dose, Analytics, Reports, Benchmarks, Compare
```

## 20 NVIDIA Technologies (All Free)

| # | Technology | Category | License |
|---|---|---|---|
| 1 | MONAI Core | Medical AI Framework | Apache 2.0 |
| 2 | MONAI Deploy SDK (MAPs) | Clinical Packaging | Apache 2.0 |
| 3 | MONAI Label | Interactive Annotation | Apache 2.0 |
| 4 | NVIDIA FLARE | Federated Learning | Apache 2.0 |
| 5 | NVIDIA AIQ Toolkit | Agentic AI | Open source |
| 6 | NeMo Guardrails | Clinical Safety | Apache 2.0 |
| 7 | RAPIDS cuDF/cuML | GPU Analytics | Apache 2.0 |
| 8 | cuVS (CAGRA) | Vector Search | Apache 2.0 |
| 9 | cuCIM | Image Processing | Apache 2.0 |
| 10 | NVIDIA DALI | Data Loading | Apache 2.0 |
| 11 | Triton Inference Server | Model Serving | BSD 3-Clause |
| 12 | TensorRT | Inference Optimization | Free SDK |
| 13 | NVIDIA Dynamo | LLM Serving | Apache 2.0 |
| 14 | Holoscan SDK | Real-time Streaming | Apache 2.0 |
| 15 | NV-Segment-CT | Segmentation (132 classes) | Open Model |
| 16 | Llama-3 8B | On-device LLM | Meta Community |
| 17 | PyRadiomics-CUDA | Radiomics | BSD |
| 18 | torchxrayvision | CXR Classification | Apache 2.0 |
| 19 | Nemotron Nano | Edge LLM | Open weights |
| 20 | NV-Generate CT/MR | Synthetic Data | Check release |

## Knowledge Base

| Collection | Vectors | Description |
|---|---|---|
| imaging_literature | 1,938 | Real PubMed research papers (AI + medical imaging) |
| imaging_trials | 47 | ClinicalTrials.gov imaging AI studies |
| imaging_findings | 80 | Imaging finding templates and patterns |
| imaging_protocols | 55 | Acquisition protocols and parameters |
| imaging_devices | 50 | FDA-cleared AI/ML medical devices |
| imaging_anatomy | 40 | Anatomical structure references |
| imaging_benchmarks | 45 | Model performance benchmarks |
| imaging_guidelines | 35 | ACR, RSNA, ESR clinical guidelines |
| imaging_report_templates | 30 | Structured radiology report templates |
| imaging_datasets | 30 | Public datasets (TCIA, PhysioNet) |
| imaging_radiomics | — | PyRadiomics feature vectors |
| imaging_reports | — | Parsed radiology reports |
| genomic_evidence | 35,678 | Shared genomic variant evidence (read-only) |
| **Total** | **38,028** | |

## Demo Guides

| Demo | Duration | Description | Download |
|---|---|---|---|
| **Demo 1: Imaging Engine** | 21 min | Engine 4 standalone — 9 workflows, 3D visualization, cross-modal triggers, evidence RAG, protocol optimization | [:material-download: Demo Guide 1](Clinical_Imaging_Engine_Demo_Guide_1_(1_eng).docx) |
| **Demo 2: Closed Loop** | 27 min | CT scan → genomic analysis → 100 drug candidates — all 4 engines on one device | [:material-download: Demo Guide 2](Clinical_Imaging_Engine_Demo_Guide_2_(mult_eng).docx) |

## Services

| Port | Service | Description |
|---|---|---|
| 8550 | React Portal | Clinical-grade UI (10 pages) |
| 8525 | Streamlit UI | Developer workbench (10 tabs) |
| 8524 | FastAPI REST API | 33+ endpoints |
| 8520 | Llama-3 8B NIM | On-device LLM |
| 8527 | MONAI Label | Interactive annotation |
| 8530 | VISTA-3D NIM | 3D segmentation |
| 8531 | MAISI NIM | Synthetic CT generation |
| 8532 | VILA-M3 NIM | Vision-language model |
| 8534 | NV-Segment-CT | 132-class CT segmentation |
| 19530 | Milvus gRPC | Vector database |
| 8042 | Orthanc REST | DICOM server |
| 4242 | Orthanc DICOM | C-STORE/C-FIND |
| 3000 | Grafana | Monitoring dashboards |
| 9099 | Prometheus | Metrics |

## Key Numbers

| Metric | Value |
|---|---|
| NVIDIA technologies | 20 (all free) |
| Clinical workflows | 9 |
| Scoring systems | 7 |
| Cross-modal triggers | 8 |
| Vector collections | 13 |
| Indexed vectors | 38,028 |
| PubMed papers | 1,938 |
| SNOMED CT codes | 103 |
| Export formats | 5 |
| Tests passing | 1,324 |
| Demo cases | 9 |
| Hardware | DGX Spark ($4,699) |
| Software cost | $0 |
| License | Apache 2.0 |

## Credits

- **Adam Jones** — Author
- **Apache 2.0 License**
- Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory)

---

!!! warning "Clinical Decision Support Disclaimer"
    The Clinical Imaging Engine is a clinical decision support research tool. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
