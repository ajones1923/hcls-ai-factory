# Imaging Intelligence Agent

Automated detection, segmentation, longitudinal tracking, and clinical triage of CT, MRI, and chest X-ray studies on NVIDIA DGX Spark. Part of the [HCLS AI Factory](https://github.com/ajones1923/hcls-ai-factory).

![HCLS AI Factory Imaging AI Agent on NVIDIA DGX Spark](HCLS%20AI%20Factory%20Imaging%20AI%20Agent%20on%20NVIDIA%20DGX%20Spark%20Infographic.png)

**Source:** [github.com/ajones1923/imaging-intelligence-agent](https://github.com/ajones1923/imaging-intelligence-agent)

## Overview

The Imaging Intelligence Agent processes medical imaging studies using MONAI models and NVIDIA NIM microservices on DGX Spark hardware. A multi-collection RAG engine backed by 11 Milvus vector collections provides evidence-grounded clinical reasoning, while six reference workflows cover the highest-impact radiology use cases. Cross-modal triggers connect imaging findings to 3.5M genomic variant vectors for precision medicine enrichment. Output is exported as Markdown, JSON, PDF, or FHIR R4 DiagnosticReport Bundles with SNOMED CT, LOINC, and DICOM coding.

## Six Reference Workflows

| Workflow | Modality | Model (Pretrained Weights) | Key Output |
|---|---|---|---|
| **Hemorrhage Triage** | CT Head | SegResNet (MONAI `wholeBody_ct_segmentation`) | Volume (mL), midline shift (mm), urgency routing |
| **Lung Nodule Tracking** | CT Chest | RetinaNet + SegResNet (MONAI `lung_nodule_ct_detection`) | Lung-RADS 1--4B classification |
| **Coronary Angiography** | CT Cardiac | Calcium scoring + stenosis grading | CAD-RADS classification |
| **Rapid Findings** | CXR | DenseNet-121 (torchxrayvision `densenet121-res224-all`, CheXpert) | Multi-label classification + GradCAM heatmaps |
| **MS Lesion Tracking** | MRI Brain | UNEST (MONAI `wholeBrainSeg_Large_UNEST_segmentation`) | Lesion count, disease activity (Stable/Active/Highly Active) |
| **Prostate PI-RADS** | MRI Prostate | PI-RADS scoring and lesion characterization | PI-RADS classification |

## Architecture

```
DICOM Study Arrives (Orthanc 8042/4242)
    |
    v
[Webhook Router] -- CT+head / CT+chest / CT+cardiac / CR+chest / MR+brain / MR+prostate
    |
    v
[Clinical Workflow]
(SegResNet / RetinaNet / DenseNet-121 / UNEST)
    |
    v
[Post-Processing + Cross-Modal Trigger]
Volume, midline shift, Lung-RADS, CAD-RADS, GradCAM, PI-RADS
Lung-RADS 4A+ -> genomic variant queries (3.5M vectors)
    |
    v
[RAG Engine + Claude LLM]
11 Milvus collections + Claude Sonnet 4.6/Llama-3 synthesis
    |
    v
[Clinical Output]
Markdown | JSON | PDF | FHIR R4 DiagnosticReport Bundle
```

Built on the HCLS AI Factory platform:

- **RAG Engine:** Multi-collection Milvus vector search + Claude/Llama-3 LLM synthesis
- **Embeddings:** BGE-small-en-v1.5 (384-dim, IVF_FLAT, COSINE)
- **Database:** Milvus 2.4 (11 collections -- 10 imaging-specific + `genomic_evidence` read-only)
- **LLM:** Claude Sonnet 4.6 (Anthropic) primary, Llama-3 8B NIM fallback
- **NIM Services:** VISTA-3D (segmentation), MAISI (synthetic CT), VILA-M3 (VLM), Llama-3 8B (LLM)
- **DICOM Server:** Orthanc (webhook auto-routing to workflows)
- **UI:** 9-tab Streamlit Imaging Workbench (port 8525)
- **API:** FastAPI (port 8524)
- **Export:** Markdown, JSON, PDF, FHIR R4 DiagnosticReport Bundle
- **Hardware target:** NVIDIA DGX Spark ($4,699)

## Knowledge Base

| Source | Records |
|---|---|
| Seed vectors across 10 collections | 876 |
| Genomic evidence vectors (read-only) | 3,561,170 |
| **Total vectors** | **3,562,046** |

## Cross-Modal Integration

The Imaging Agent connects into the broader HCLS AI Factory genomics pipeline:

- **Lung-RADS 4A+** triggers EGFR/ALK/ROS1/KRAS genomic variant queries against 3.5M genomic vectors
- **CXR urgent findings** (consolidation, critical severity) trigger infection genomics queries
- **Brain lesion high activity** triggers neurological genomics queries (HLA-DRB1, demyelination markers)
- All cross-modal results are included in FHIR R4 export with SNOMED CT, LOINC, and DICOM coding

## Clinical Output

| Output | Format | Usage |
|---|---|---|
| Clinical Report | Markdown | Human-readable structured report |
| Structured Data | JSON | Programmatic consumption, dashboards |
| Printable Report | PDF | Clinical documentation, patient records |
| Interoperability | FHIR R4 DiagnosticReport Bundle | EHR integration (SNOMED CT + LOINC + DICOM coded) |

## Services

| Port | Service |
|---|---|
| 8524 | FastAPI REST Server |
| 8525 | Streamlit Imaging Workbench |
| 8520 | NIM LLM (Llama-3 8B) |
| 8530 | NIM VISTA-3D |
| 8531 | NIM MAISI |
| 8532 | NIM VILA-M3 |
| 19530 | Milvus (gRPC) |
| 8042 | Orthanc REST API |
| 4242 | Orthanc DICOM C-STORE |

## Demo Highlights

- **4 pre-loaded demo cases** (DEMO-001 through DEMO-004) with mock-realistic clinical results
- **Image Gallery** with 5 CXR pathologies, cross-modality gallery, and 3D volume slice viewer
- **Annotated AI images** in Workflow Runner with 6-step pipeline animation
- **Pre-filled example queries** in Evidence Explorer and Protocol Advisor
- **Plotly donut chart** for collection statistics
- **Interactive network graph** in Patient 360 showing gene-finding relationships
- **Sidebar guided tour** with 9-step demo flow
- **876 seed vectors** across 10 imaging-specific collections

## Credits

- **Adam Jones**
- **Apache 2.0 License**
