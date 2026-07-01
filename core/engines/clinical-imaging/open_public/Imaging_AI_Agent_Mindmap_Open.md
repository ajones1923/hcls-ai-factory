# Imaging AI Agent (CT / MRI / X-Ray) — Open Architecture

## Canonical Imaging State

### Original DICOM (Immutable)
- CT / MRI / X-Ray studies preserved as evidence
- Ingested via DICOMweb STOW-RS
- Ingested via DIMSE C-STORE
- Event-triggered pull from PACS/VNA

### Derived Artifacts
- Segmentation masks (volumetric)
- Measurements (diameters, volumes, shifts)
- Heatmaps (GradCAM localization)
- Spatial registrations

### Semantic Embeddings
- Study-level vectors (384-dim)
- Series-level vectors
- Lesion-level vectors
- Cohort retrieval and case matching

### Structured Findings
- Queryable clinical state
- DICOM SR compatible
- Confidence scores and thresholds
- Guideline mappings (Lung-RADS, etc.)

### Provenance Bundles
- Model ID and version
- Inference parameters and durations
- Input data lineage (DICOM UIDs)
- Timestamps and operator approvals
- Predetermined change control plans

## DGX Spark Software Architecture

### Data Layer

#### Local NVMe Storage (GPUDirect)
- DGX Spark 128 GB unified memory + NVMe
- GPUDirect Storage for zero-copy GPU access
- Immutable DICOM archive on local flash
- Derived artifacts persisted alongside originals
- Linear I/O scaling with NVMe capacity

#### Local Filesystem Namespace
- Single-site namespace for MVP proof build
- Prior study retrieval within local archive
- Longitudinal memory within institution
- Path to multi-site via shared storage at scale

### Execution Layer

#### MONAI Deploy + Container Orchestration
- MONAI Deploy Application Packages (MAPs)
- Containerized, portable inference pipelines
- Nextflow / Airflow for pipeline DAG orchestration
- study.complete event triggers pipeline execution
- Outputs persisted as canonical artifacts

#### Agent Framework (LangChain / LangGraph)
- Multi-step clinical reasoning workflows
- MCP (Model Context Protocol) tool integration
- Agent personas with distinct roles
  - Triage agent
  - Longitudinal tracker
  - Population analyst
- Checkpointing and observability

### Reasoning Layer

#### PostgreSQL + pgvector
- Structured finding queries ("all Lung-RADS 4A+")
- Vector similarity search ("10 most similar CT studies")
- Hybrid queries ("growing nodules AND similar phenotype AND APOE4 carrier")
- HNSW and IVFFlat indexing for fast retrieval
- Open-source, no proprietary database required

#### RAG Pipeline (NVIDIA NIM LLM Serving)
- Retrieval-augmented generation
- Evidence grounding — ACR guidelines, prior measurements, similar outcomes
- Longitudinal comparison and delta analysis
- Cross-modal context enrichment (imaging + genomics + biomarkers)
- Evidence-grounded reporting

## NVIDIA Accelerated Computing

### DGX Compute Progression
- DGX Spark — proof build ($4,699, GB10 Grace Blackwell, 128 GB unified memory)
- DGX B200 — departmental ($500K-$515K per 8-GPU system)
- DGX SuperPOD — enterprise ($7M-$60M+, thousands of concurrent studies)
- DGX Cloud — managed GPU subscription (hyperscaler availability)

### MONAI Deploy
- Application Packages (MAPs)
- Containerized inference pipelines
- Portable and versioned deployment
- Orchestrated by Nextflow / container runtime

### MONAI Model Zoo
- 3D U-Net (segmentation)
- RetinaNet (detection)
- SegResNet (segmentation)
- DenseNet-121 (classification)
- SwinUNETR (transformer-based)
- Fine-tuning starting points

### NVIDIA NIM
- Standardized inference microservices
- Versioned model serving
- Health checks and auto-scaling
- Production-grade deployment
- Licensed via NVIDIA AI Enterprise ($4,500/GPU/year)

### NVIDIA FLARE
- Federated learning
- Multi-site model improvement
- No centralization of patient data
- Privacy-preserving training
- Free (Apache 2.0) — drives DGX pull-through

## Reference Agent Workflows

### CT Head — Hemorrhage Triage
- Target: < 90 seconds end-to-end
- Sensitivity: > 95% for hemorrhage > 5 mL
- Detection -> Segmentation -> Volume estimation
- Midline shift measurement
- Urgency classification (critical / urgent / routine)
- Worklist prioritization and notification

### CT Chest — Lung Nodule Tracking
- Target: < 5 minutes multi-stage
- Detection sensitivity: > 90% for nodules >= 4 mm
- Per-nodule segmentation and volumetrics
- Lung-RADS category assignment
- Malignancy risk scoring
- Longitudinal matching to prior CT
- Volume doubling time calculation
- Triggers genomics pipeline if Lung-RADS 4B+

### CXR — Rapid Findings
- Target: < 30 seconds end-to-end
- Pneumothorax sensitivity: > 95%
- Multi-label classification
  - Pneumothorax
  - Consolidation
  - Pleural effusion
  - Cardiomegaly
  - Fracture
- GradCAM heatmap localization
- Immediate high-risk flagging

### MRI Brain — MS Lesion Tracking
- Target: < 5 minutes multi-stage
- 3D U-Net lesion segmentation on FLAIR
- Lesion count and volume measurement
- Spatial registration to prior MRI
- New / enlarging lesion identification
- Disease activity assessment (stable / active / highly active)
- Longitudinal trajectory tracking

## Clinical Integration

### Inputs from Clinical Systems
- PACS (system of record for viewing/worklists)
- VNA (vendor-neutral archive)
- RIS (radiology information system)

### Output to PACS
- DICOM SR (structured findings)
- GSPS overlays (visual annotations)
- Pushed via DICOMweb STOW-RS

### Output to EHR
- FHIR DiagnosticReport
- Structured clinical summaries

### Triage and Routing
- Worklist prioritization by urgency
- On-call notification by policy thresholds
- Escalation for critical findings

### Clinician-in-the-Loop
- Decision support, not autonomous diagnosis
- FDA AI/ML SaMD framework aligned
- All outputs reviewable and configurable
- Clinicians remain accountable

## Trust, Governance, and Compliance

### Provenance by Default
- Every output traceable to exact model + data + config
- Immutable audit trail

### Reproducibility
- Deterministic re-runs on canonical data
- Previous version outputs preserved alongside new

### Regulatory Alignment
- FDA AI/ML SaMD framework
- Predetermined change control plans
- Controlled rollouts for model updates

### Security and Multi-Tenancy
- Least-privilege access enforced at the data layer
- Tenant isolation
- Patient data remains within institutional control
- Agent personas with distinct security credentials

### Observability
- Pipeline traces (durations, failures, throughput)
- Model performance monitoring
- Drift detection

## HCLS AI Factory Integration

### Imaging to Genomics (Parabricks)
- Lung nodule confirmed malignant -> triggers tumor profiling (somatic / germline)
- Parabricks: 30x WGS from ~30 hours (CPU) to ~10 minutes on DGX (8-GPU)
- Up to 50% lower compute cost (NVIDIA Parabricks documentation)
- Imaging phenotype to molecular characterization

### Imaging to Clinical Reasoning (NIM LLM Serving)
- Evidence retrieval for clinical context
- Guidelines and literature grounding
- Cross-modal reasoning via RAG pipeline

### Imaging to Drug Discovery (BioNeMo)
- 200+ adopters (Eli Lilly, Astellas, Insilico, Recursion)
- Drug candidate scoring by imaging phenotype
- Treatment stratification
- Quantitative imaging endpoints for clinical trials (RECIST)

### Imaging to Biomarker Intelligence Agent
- Cross-agent data flows on shared DGX Spark
- Genomic + imaging biomarker fusion
- Combined phenotype profiling

### Imaging to Longitudinal Care
- Continuous monitoring across time points
- Automated detection of meaningful change
- Population-scale cohort analysis

### Imaging to Outcomes
- Cohort retrieval (patients like this)
- Imaging trajectories linked to outcomes
- Care pathway optimization

## Open-Source Technology Stack

### Core Frameworks (Apache 2.0 / BSD / MIT)
- MONAI — Medical imaging AI framework (Apache 2.0)
- MONAI Deploy — Containerized inference packaging (Apache 2.0)
- NVIDIA FLARE — Federated learning (Apache 2.0)
- LangChain / LangGraph — Agent orchestration (MIT)
- Nextflow — Pipeline orchestration (Apache 2.0)
- PostgreSQL — Relational database (PostgreSQL License)
- pgvector — Vector similarity search (PostgreSQL License)
- Orthanc — DICOM server (GPLv3)
- dcm4chee — DICOM archive (MPL 1.1 / GPL 2.0 / LGPL 2.1 triple license)
- pydicom — DICOM parsing (MIT)
- HAPI FHIR — FHIR server (Apache 2.0)

### NVIDIA Software (NVAIE Licensed)
- NVIDIA NIM — Inference microservices ($4,500/GPU/year)
- NVIDIA Parabricks — Genomics acceleration ($4,500/GPU/year)
- NVIDIA BioNeMo — Drug discovery ($4,500/GPU/year)

### Hardware
- NVIDIA DGX Spark — GB10 Grace Blackwell, 128 GB unified memory ($4,699)
- NVIDIA DGX B200 — 8x B200 GPUs ($500K-$515K)
- NVIDIA DGX SuperPOD — Thousands of GPUs ($7M-$60M+)

## Deployment Roadmap

### Phase 1 — Proof Build
- Single DGX Spark ($4,699)
- 1-2 workflows (CT head, CXR)
- MONAI Deploy MAPs validated
- Canonical data model proven
- Zero NVAIE software cost (desktop-class)

### Phase 2 — Departmental
- 1-2x DGX B200 cluster ($500K-$1M)
- Shared filesystem namespace
- Multi-user, multi-modality
- PACS integration via DICOMweb
- Clinical validation with radiologists
- NVAIE: $36K-$72K/year (8-16 GPUs)

### Phase 3 — Multi-Site Enterprise
- 4-8x DGX B200 + InfiniBand ($2M-$4M)
- Cross-site namespace via shared storage
- Continuous reprocessing
- Population-scale cohort retrieval
- NVIDIA FLARE federated learning
- NVAIE: $144K-$288K/year (32-64 GPUs)

### Phase 4 — AI Factory at Scale
- DGX SuperPOD ($7M-$60M+)
- Unified multimodal agent fabric
- Imaging + Genomics + Outcomes + Therapy
- Complete HCLS AI Factory operating as one platform
- Thousands of concurrent studies
- NVAIE: $576K-$1.15M/year (128-256 GPUs)
