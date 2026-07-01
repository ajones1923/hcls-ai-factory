# Imaging AI Agent (CT / MRI / X-Ray) — v2

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
- disaggregated storage architecture
- Stateless compute nodes + persistent storage nodes
- Universal hash space with zero east-west traffic
- Sustained multi-pass throughput for repeated inference
- Linear scaling — add compute nodes for proportional I/O bandwidth
- Scales to exabytes of flash

#### namespace
- Global namespace (edge / core / cloud)
- Prior study retrieval across sites
- Longitudinal memory across institutions
- Policy-driven placement
- Multi-site archives become a single namespace

### Execution Plane

#### data engine
- Event-driven orchestration
- study.complete triggers pipeline
- Up to 2 million events per second per compute node
- Serverless, event-driven — no external orchestrator
- Pipeline DAGs with containerized execution
- Outputs persisted as canonical artifacts

#### agent engine
- Low-code AI agent deployment runtime
- Runs natively within data engine
- Studio (low-code environment) for workflow design
- MCP (Model Context Protocol) Tool Server
- Agent personas with distinct security credentials
  - Triage agent
  - Longitudinal tracker
  - Population analyst
- Fault-tolerant queuing with checkpointing
- Full observability and feedback loops

### Reasoning Plane

#### Vector Database
- Unified SQL + vector + analytics in one query plane
- Structured finding queries ("all Lung-RADS 4A+")
- Vector similarity search ("10 most similar CT studies")
- Hybrid queries ("growing nodules AND similar phenotype AND APOE4 carrier")
- GPU-accelerated CAGRA indexing
- Sub-millisecond vector search
- Replaces separate relational, vector, and analytics databases

#### retrieval engine
- Retrieval-augmented generation (RAG)
- Evidence grounding — ACR guidelines, prior measurements, similar outcomes
- Longitudinal comparison and delta analysis
- Cross-modal context enrichment (imaging + genomics + biomarkers)
- Evidence-grounded reporting

## NVIDIA Accelerated Computing

### DGX Compute Progression
- DGX Spark — proof build ($4,699, GB10 Grace Blackwell, 128 GB unified memory)
- DGX B200 — departmental ($500K–$515K per 8-GPU system)
- DGX SuperPOD — enterprise ($7M–$60M+, thousands of concurrent studies)
- DGX Cloud — managed GPU subscription (starts ~$37K/mo for 8-GPU)

### Hardware Acceleration

#### BlueField-4 DPU
- 64-core Arm Neoverse V2
- ConnectX-9 at 800 Gbps
- the platform compute node runs directly on DPU
- Every DGX node becomes storage access point and GPU compute endpoint
- DICOM data flows from NVMe flash through DPU directly to GPU memory at line rate
- Eliminates external storage controllers

#### ICMS / KV Cache
- NVIDIA Inference Context Memory Storage
- Hierarchical KV Block Manager
  - Tier 1: GPU memory (fastest)
  - Tier 2: CPU memory
  - Tier 3: NVMe flash (managed by the platform disaggregated storage)
- KV cache persists across inference turns
- Patient with four prior studies maintains conversational context
- No re-ingestion of historical findings between turns

#### NVIDIA Dynamo
- 90% improvement in GPU efficiency (the platform + Dynamo benchmark, 2025)
- 20x improvement in time-to-first-token vs GPU-memory-only caching
- Optimized LLM serving layer for retrieval engine RAG reasoning
- Faster clinical triage and higher concurrent throughput

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
- Detection → Segmentation → Volume estimation
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
- Least-privilege access enforced at the data plane
- Tenant isolation
- Patient data never leaves the AI platform
- Agent personas with distinct security credentials

### Observability
- Pipeline traces (durations, failures, throughput)
- Model performance monitoring
- Drift detection

## HCLS AI Factory Integration

### Imaging to Genomics (Parabricks)
- Lung nodule confirmed malignant → triggers tumor profiling (somatic / germline)
- Parabricks: 30x WGS from ~30 hours (CPU) to ~10 minutes on single GPU
- Up to 50% lower compute cost (NVIDIA Parabricks documentation)
- Imaging phenotype to molecular characterization

### Imaging to Clinical Reasoning (Med42-70B on NIM)
- Evidence retrieval for clinical context
- Guidelines and literature grounding
- Cross-modal reasoning via retrieval engine RAG

### Imaging to Drug Discovery (BioNeMo)
- 200+ adopters (Eli Lilly, Astellas, Insilico, Recursion)
- Drug candidate scoring by imaging phenotype
- Treatment stratification
- Quantitative imaging endpoints for clinical trials (RECIST)

### Imaging to Biomarker Intelligence Agent
- Cross-agent data flows on shared the platform substrate
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

## Competitive Differentiation

### the AI platform (Native Substrate)
- Storage + Execution + Query + Reasoning in one platform
- No external glue layers required
- Canonical data model built in
- Single operational surface
- agent engine: native agent runtime with MCP, personas, observability

### DDN / WEKA (Assembled Stack)
- High-performance file storage
- MCP-ready storage (WEKA); full agent runtime requires external frameworks
- Requires external workflow orchestrator (Airflow/Argo for complex pipelines)
- Requires external relational database
- Requires external vector database (Milvus)
- Reference architectures for RAG (DDN CORE RAG, WEKA WARRP) but LLM serving external
- Active BF-3/BF-4 integrations but DPU runs as add-on, not native data plane
- Per-cluster namespace; cross-site via snapshot replication

## Deployment Roadmap

### Phase 1 — Proof Build
- Single DGX Spark ($4,699)
- 1–2 workflows (CT head, CXR)
- MONAI Deploy MAPs validated
- Canonical data model proven
- Zero NVAIE software cost (desktop-class)

### Phase 2 — Departmental
- 1–2x DGX B200 cluster ($500K–$1M)
- Shared namespace namespace
- Multi-user, multi-modality
- PACS integration via DICOMweb
- Clinical validation with radiologists
- NVAIE: $36K–$72K/year (8–16 GPUs)

### Phase 3 — Multi-Site Enterprise
- 4–8x DGX B200 + InfiniBand ($2M–$4M)
- Global namespace across edge/core/cloud
- Continuous reprocessing
- Population-scale cohort retrieval
- NVIDIA FLARE federated learning
- NVAIE: $144K–$288K/year (32–64 GPUs)

### Phase 4 — AI Factory at Scale
- DGX SuperPOD ($7M–$60M+)
- Unified multimodal agent fabric
- Imaging + Genomics + Outcomes + Therapy
- Complete HCLS AI Factory operating as one substrate
- Thousands of concurrent studies
- NVAIE: $576K–$1.15M/year (128–256 GPUs)

## Revenue Opportunity (5-Year)

### NVIDIA Revenue — Imaging Agent
- 200 US large hospitals (>250 beds)
- Hardware: $225M (DGX Spark → B200 → multi-site)
- Software (NVAIE): $35M
- US total: $260M
- Global base (2.5x): $650M

### the platform Revenue — Imaging Agent
- Same 200 hospitals
- Platform revenue: $97M (US)
- Global base (2.5x): $243M

### Combined HCLS AI Factory (Imaging + Genomics)
- Conservative (imaging 2x + 50 pharma direct): $1.13B
- Base (imaging 2.5x + 50 pharma direct): $1.31B
- Aggressive (imaging 3.5x + DGX Cloud + co-innovation): $1.8B–$2.5B
