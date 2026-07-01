# Imaging AI Agent (CT / MRI / X-Ray)

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

## the AI platform (Unified Substrate)
### Data Plane
#### Object Storage
- Sustained multi-pass throughput
- Repeated inference over large datasets
- No copy-and-stage pipeline sprawl
#### namespace
- Global namespace (edge / core / cloud)
- Prior study retrieval
- Longitudinal memory across sites
- Multi-site / multi-department addressing
### Execution Plane
#### data engine
- Event-driven orchestration
- study.complete triggers pipeline
- Containerized execution close to data
- Outputs persisted as canonical artifacts
#### agent engine
- Low-code AI agent deployment runtime
- Runs natively within data engine
- Studio (low-code environment) for workflow design
- MCP (Model Context Protocol) Tool Server
- Agent personas with distinct security credentials
- Fault-tolerant queuing with checkpointing
- Observability and feedback loops
### Reasoning Plane
#### Vector Database
- Unified SQL + vector query plane
- Structured finding queries
- Vector similarity search
- Single plane replaces separate databases
#### retrieval engine
- Retrieval-augmented reasoning
- Evidence grounding
- Longitudinal comparison and delta analysis
- Cross-modal context enrichment

## NVIDIA Accelerated Computing
### DGX Spark to SuperPOD
- Proof build on DGX Spark
- Department scale on DGX cluster
- Enterprise scale on DGX SuperPOD
- Thousands of concurrent studies
### MONAI Deploy
- Application Packages (MAPs)
- Containerized inference pipelines
- Portable and versioned deployment
- Orchestrated by data engine
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
### NVIDIA FLARE
- Federated learning
- Multi-site model improvement
- No centralization of patient data
- Privacy-preserving training

## Reference Agent Workflows
### CT Head — Hemorrhage Triage
- Target: < 90 seconds end-to-end
- Sensitivity: > 95% for hemorrhage > 5 mL
- Detection → Segmentation → Volume estimation
- Midline shift measurement
- Urgency classification (critical / urgent / routine)
- Worklist prioritization and notification
### CT Chest — Lung Nodule Tracking
- Target: < 5 minutes multi-stage
- Detection sensitivity: > 90% for nodules ≥ 4 mm
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
- Least-privilege access
- Tenant isolation
- Patient data never leaves the AI platform
### Observability
- Pipeline traces (durations, failures, throughput)
- Model performance monitoring
- Drift detection

## HCLS AI Factory Integration
### Imaging → Genomics (Parabricks)
- Lung nodule confirmed malignant
- Triggers tumor profiling (somatic / germline)
- Imaging phenotype to molecular characterization
### Imaging → RAG / Clinical Chat (NIM-served LLMs)
- Evidence retrieval for clinical context
- Guidelines and literature grounding
- Cross-modal reasoning
### Imaging → Drug Discovery (BioNeMo)
- Drug candidate scoring
- Treatment stratification by imaging phenotype
- Quantitative endpoints for trials
### Imaging → Longitudinal Care
- Continuous monitoring across time points
- Automated detection of meaningful change
- Population-scale cohort analysis
### Imaging → Outcomes
- Cohort retrieval (patients like this)
- Imaging trajectories linked to outcomes
- Care pathway optimization

## Competitive Differentiation
### the AI platform (Native Substrate)
- Storage + Execution + Query + Reasoning in one platform
- No external glue layers required
- Canonical data model built in
- Single operational surface
### DDN / WEKA (Assembled Stack)
- High-performance file storage
- Requires external workflow orchestrator
- Requires external relational database
- Requires external vector database
- Requires external reasoning engine
- Higher integration complexity

## Deployment Roadmap
### Phase 1 — Proof Build
- Single DGX Spark
- 1–2 workflows (CT head, CXR)
- MONAI Deploy MAPs validated
- Canonical data model proven
### Phase 2 — Departmental
- Small DGX cluster
- Shared namespace namespace
- Multi-user, multi-modality
- PACS integration via DICOMweb
- Clinical validation with radiologists
### Phase 3 — Multi-Site Enterprise
- Global namespace across edge/core/cloud
- Continuous reprocessing
- Population-scale cohort retrieval
- NVIDIA FLARE federated learning
### Phase 4 — AI Factory at Scale
- DGX SuperPOD
- Unified multimodal agent fabric
- Imaging + Genomics + Outcomes + Therapy
- Complete HCLS AI Factory operating as one substrate
