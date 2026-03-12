# Democratizing Medical Imaging AI: A Multi-Collection RAG Architecture with NVIDIA NIM Inference on a Single DGX Spark

**Author:** Adam Jones
**Date:** March 2026
**Version:** 2.0.0
**License:** Apache 2.0

Part of the HCLS AI Factory -- an end-to-end precision medicine platform.
https://github.com/ajones1923/hcls-ai-factory

---

## Abstract

Medical imaging AI has reached clinical maturity, with over 800 FDA-cleared AI/ML-enabled devices as of 2025, yet deployment remains concentrated in well-resourced academic medical centers. Community hospitals, rural health systems, and research institutions in low- and middle-income countries lack the GPU infrastructure, engineering expertise, and capital required to deploy, integrate, and maintain these systems. This paper presents the Imaging Intelligence Agent, an open-source medical imaging AI platform that combines multi-collection retrieval-augmented generation (RAG) with NVIDIA NIM on-device inference to deliver clinical decision support for radiology on a single NVIDIA DGX Spark ($3,999).

The system maintains 10 imaging-specific Milvus vector collections and one read-only genomic evidence collection (11 total), indexed with 384-dimensional BAAI/bge-small-en-v1.5 embeddings under IVF_FLAT indexing (nlist=1024, nprobe=16) with cosine similarity. A 3-domain knowledge graph spanning 25+ pathologies, 9 imaging modalities, and 21 anatomical structures (with SNOMED CT and FMA codes) augments vector retrieval with structured clinical context. Twelve query expansion maps improve recall by mapping domain keywords to semantically related terms. Four NVIDIA NIM microservices -- Llama-3 8B Instruct (port 8520), VISTA-3D (port 8530), MAISI (port 8531), and VILA-M3 (port 8532) -- provide on-device 3D segmentation, synthetic CT generation, vision-language understanding, and clinical text synthesis. Six reference clinical workflows implement end-to-end analysis pipelines for CT head hemorrhage triage, CT chest lung nodule tracking, CT coronary angiography, chest X-ray rapid findings, MRI brain MS lesion quantification, and MRI prostate PI-RADS scoring. Seven standardized scoring systems (Lung-RADS v2022, BI-RADS, TI-RADS, LI-RADS, CAD-RADS, PI-RADS, ASPECTS) ensure reproducible severity classification. Cross-modal genomic integration automatically triggers precision medicine queries when imaging findings exceed severity thresholds -- for example, a Lung-RADS 4A+ lung nodule triggers EGFR/ALK/ROS1 variant queries against 3.5 million genomic evidence vectors. The system exports to Markdown, JSON, PDF (ReportLab), and FHIR R4, and includes 620 tests across 12 modules. A 9-tab Streamlit interface and 19+ FastAPI endpoints across 8 routers provide both interactive and programmatic access. The entire stack deploys on a single DGX Spark with 128 GB unified LPDDR5x memory, requiring no cloud dependency and preserving full data sovereignty.

---

## 1. Introduction

### 1.1 The Challenge

Medical imaging generates approximately 3.6 billion diagnostic imaging examinations annually worldwide. AI has demonstrated radiologist-level performance on tasks ranging from chest X-ray triage to brain tumor segmentation. The FDA has cleared over 800 AI/ML-enabled medical devices, and the market is projected to exceed $20 billion by 2030. Yet deployment remains starkly unequal: the vast majority of clinical AI implementations are concentrated in large academic medical centers in high-income countries, while the hospitals that serve most of the world's patients have no access to these tools.

Three barriers explain this deployment gap:

1. **Infrastructure cost.** Production-grade medical imaging AI requires GPU servers costing $30,000-$200,000, plus networking, storage, and cooling infrastructure. Cloud inference solves the capital problem but introduces recurring costs and data sovereignty concerns.
2. **Integration complexity.** Deploying a single AI model requires PACS integration, HL7/FHIR interfacing, DICOM routing, result visualization, and clinical workflow embedding -- each demanding specialized engineering expertise that most community hospitals lack.
3. **Knowledge fragmentation.** Radiology AI knowledge is scattered across PubMed literature, ClinicalTrials.gov records, FDA device databases, ACR/RSNA guidelines, performance benchmarks, and institutional protocols. No single system synthesizes this information into actionable clinical intelligence.

### 1.2 Democratizing Medical Imaging AI

The NVIDIA DGX Spark addresses the first barrier with a $3,999 system featuring a GB10 GPU (Blackwell architecture), 128 GB unified LPDDR5x memory, and 20 ARM cores (Grace architecture) -- sufficient to run multiple AI inference models simultaneously with no cloud dependency. The Imaging Intelligence Agent addresses the second and third barriers by providing a complete, open-source platform that combines evidence retrieval, inference orchestration, and clinical workflow automation in a single deployable package.

### 1.3 Design Goals

The system was designed around six principles:

1. **Single-device deployment.** The entire stack -- vector database, embedding model, four NIM inference services, PACS server, web viewer, and application layer -- runs on one DGX Spark.
2. **No cloud dependency.** All inference happens on-device. No patient data, imaging studies, or clinical queries leave the device. This architecture supports deployment in HIPAA- and GDPR-regulated environments where cloud-based inference faces compliance barriers.
3. **Graceful degradation.** The system operates in three tiers: full mode (all NIMs active), lite mode (RAG-only, no GPU required), and mock mode (synthetic results for demonstration and testing). Degradation is automatic and transparent.
4. **Clinical relevance.** Six reference workflows address high-impact radiology use cases with standardized scoring systems and evidence-based severity classification.
5. **Open source.** Apache 2.0 licensing enables institutional customization, academic research, and commercial deployment without licensing barriers.
6. **Cross-modal integration.** The architecture connects imaging findings to genomic evidence, enabling precision medicine workflows that bridge radiology and molecular diagnostics.

### 1.4 Our Contribution

This paper makes five contributions:

1. A **multi-collection RAG architecture** with 10 imaging-specific vector collections and parallel search with weighted scoring, delivering cross-domain evidence synthesis for radiology queries.
2. A **3-domain knowledge graph** with 25+ pathologies, 9 modalities, and 21 anatomical structures that augments vector retrieval with structured clinical context including SNOMED CT, FMA, and ICD-10 codes.
3. A **cross-modal genomic bridge** that automatically triggers precision medicine queries when imaging findings exceed clinical severity thresholds.
4. A **reference deployment architecture** demonstrating that enterprise-grade medical imaging AI can run on a single $3,999 device with 128 GB unified memory.
5. A **complete open-source implementation** comprising 620 tests across 12 modules, 19+ API endpoints, and 6 clinical workflows released under Apache 2.0.

---

## 2. Background

### 2.1 Medical Imaging AI Landscape

The trajectory of medical imaging AI has followed three distinct phases. The first phase (2012-2017) established deep learning as a viable approach for medical image analysis, with landmark results in diabetic retinopathy screening, skin lesion classification, and chest X-ray interpretation. The second phase (2017-2023) produced the first wave of FDA-cleared products, concentrated initially in mammography, chest X-ray triage, and stroke detection. The third phase (2023-present) is characterized by foundation models (VISTA-3D, SAM-Med, BiomedCLIP), vision-language models (VILA-M3, LLaVA-Med), and the convergence of imaging AI with broader clinical AI systems.

Despite this maturation, the field faces persistent challenges:

- **Validation gap.** Models trained on data from large academic centers often underperform in community hospital settings with different scanners, protocols, and patient populations.
- **Integration burden.** Each AI tool requires DICOM routing, PACS integration, result display, and workflow embedding -- multiplied across potentially dozens of AI applications.
- **Knowledge disconnect.** Radiologists must manually search literature, guidelines, and device databases to contextualize AI outputs, a time-consuming process that limits adoption.
- **Cost concentration.** GPU infrastructure, software licensing, and integration engineering costs price out the majority of healthcare facilities worldwide.

### 2.2 The Multi-Modal Data Problem

Medical imaging intelligence requires synthesis across fundamentally different data types:

1. **Published literature** -- PubMed abstracts, conference proceedings, systematic reviews documenting evidence for imaging AI applications.
2. **Clinical trials** -- ClinicalTrials.gov records documenting ongoing and completed studies of imaging AI interventions.
3. **Imaging findings** -- Structured descriptions of pathological findings with modality-specific characteristics, severity criteria, and measurement standards.
4. **Acquisition protocols** -- Scanner parameters, contrast timing, reconstruction kernels, and quality metrics for reproducible imaging.
5. **FDA-cleared devices** -- Regulatory status, clinical indications, performance claims, and integration requirements for AI/ML medical devices.
6. **Anatomical structures** -- Standardized anatomical references with ontology codes, imaging characteristics, and segmentation labels.
7. **Performance benchmarks** -- Model accuracy, sensitivity, specificity, and AUC metrics across datasets and clinical settings.
8. **Clinical guidelines** -- ACR, RSNA, ESR, and specialty society recommendations for imaging utilization, reporting, and quality.
9. **Report templates** -- Structured radiology report templates following ACR and institutional standards.
10. **Public datasets** -- TCIA, PhysioNet, and other open imaging datasets for model development and validation.

Each data type has a different optimal schema, different source APIs, and different relevance to different clinical questions. A single-collection vector database cannot adequately represent this heterogeneity. A multi-collection architecture, where each collection has a purpose-built schema with domain-specific metadata fields, is essential.

---

## 3. System Architecture

### 3.1 Overview

The Imaging Intelligence Agent is organized into four functional layers:

**Presentation Layer.** A 9-tab Streamlit interface (port 8525) provides interactive access through Evidence Explorer, Workflow Runner, Patient 360, Image Gallery, Protocol Advisor, Device Ecosystem, Dose Intelligence, Reports & Export, and Benchmarks tabs. A FastAPI REST server (port 8524) exposes 19+ endpoints across 8 routers for programmatic access. An OHIF viewer (port 8526) provides DICOM visualization, and an Orthanc PACS server (ports 8042/4242) handles DICOM storage and retrieval.

**Intelligence Layer.** A multi-collection RAG engine searches 10 imaging-specific Milvus collections plus one read-only genomic evidence collection in parallel, applies weighted scoring, augments results with a 3-domain knowledge graph, and synthesizes answers via Claude Sonnet (Anthropic) or Llama-3 (NIM fallback). The engine supports query expansion, comparative analysis (auto-detected "X vs Y" queries), and multi-turn conversation memory.

**Inference Layer.** Four NVIDIA NIM microservices provide GPU-accelerated inference for 3D segmentation (VISTA-3D), synthetic CT generation (MAISI), vision-language understanding (VILA-M3), and clinical text synthesis (Llama-3 8B Instruct). A unified NIM client layer provides health checking, retry logic, and automatic mock fallback. Six reference clinical workflows implement complete analysis pipelines with standardized scoring.

**Data Layer.** Milvus 2.4 stores 384-dimensional BGE-small-en-v1.5 embeddings across all collections with IVF_FLAT indexing (nlist=1024, nprobe=16) and cosine similarity. Ingest pipelines fetch from PubMed (NCBI E-utilities) and ClinicalTrials.gov (API v2), with seed data scripts for all domain-specific collections. Four demo cases provide immediate demonstration capability.

### 3.2 The 10+1 Collections

The agent maintains 10 imaging-specific Milvus collections with purpose-built schemas, plus read-only access to the genomic evidence collection populated by the HCLS AI Factory's Stage 2 RAG pipeline:

| # | Collection | Content Domain | Seed Source | Weight |
|---|-----------|----------------|-------------|--------|
| 1 | `imaging_literature` | Published research papers and reviews | PubMed (NCBI E-utilities) | 0.18 |
| 2 | `imaging_trials` | Clinical trials of imaging AI | ClinicalTrials.gov API v2 | 0.12 |
| 3 | `imaging_findings` | Imaging finding patterns and templates | Curated seed data | 0.15 |
| 4 | `imaging_protocols` | Acquisition protocols and parameters | Curated seed data | 0.08 |
| 5 | `imaging_devices` | FDA-cleared AI/ML medical devices | Curated seed data | 0.08 |
| 6 | `imaging_anatomy` | Anatomical structure references | Curated seed data | 0.06 |
| 7 | `imaging_benchmarks` | Model performance benchmarks | Curated seed data | 0.08 |
| 8 | `imaging_guidelines` | Clinical guidelines (ACR, RSNA, ESR) | Curated seed data | 0.10 |
| 9 | `imaging_report_templates` | Structured radiology report templates | Curated seed data | 0.05 |
| 10 | `imaging_datasets` | Public datasets (TCIA, PhysioNet) | Curated seed data | 0.06 |
| 11 | `genomic_evidence` | Genomic variant evidence (read-only) | HCLS AI Factory Stage 2 | 0.04 |

The weights sum to 1.00 and reflect the relative clinical importance of each domain for typical radiology queries. Literature and findings receive the highest weights (0.18 and 0.15) because they most directly answer clinical questions. The genomic evidence collection receives the lowest weight (0.04) because it serves primarily as a cross-modal enrichment source triggered by specific clinical findings rather than as a primary search target.

### 3.3 Parallel Multi-Collection Search

The RAG pipeline implements a five-stage search process:

1. **Query expansion.** Twelve domain-specific expansion maps broaden the search query with related medical terminology, abbreviations, and synonyms. For example, "ct head hemorrhage" expands to include "computed tomography," "intracranial," "bleeding," "hematoma," and "ich."
2. **Parallel retrieval.** All 11 collections are searched simultaneously using `ThreadPoolExecutor`, with `top_k_per_collection=5` results per collection and a `score_threshold=0.4` minimum cosine similarity filter. This architecture ensures that response latency scales with the slowest single-collection search rather than with the number of collections.
3. **Weighted merge.** Results from all collections are scored by multiplying the cosine similarity score by the collection weight, then sorted by weighted score. This ensures that a moderately relevant literature result (cosine 0.75, weight 0.18) can outrank a highly relevant dataset result (cosine 0.85, weight 0.06) when the query is clinical rather than methodological.
4. **Knowledge augmentation.** Matched pathology, modality, and anatomy entries from the 3-domain knowledge graph are injected as structured context, providing ICD-10 codes, severity criteria, key measurements, and guideline references that may not appear in the vector-retrieved documents.
5. **LLM synthesis.** The query, retrieved evidence with citation metadata, and knowledge context are assembled into a prompt for Claude Sonnet (Anthropic API) as the primary LLM, with Llama-3 8B (NIM) as the fallback. The synthesized response includes evidence citations with relevance scores and source provenance.

### 3.4 Embedding Strategy

All collections use the same embedding configuration:

| Parameter | Value |
|-----------|-------|
| Model | BAAI/bge-small-en-v1.5 |
| Dimensions | 384 |
| Index type | IVF_FLAT |
| nlist (index) | 1024 |
| nprobe (search) | 16 |
| Similarity metric | COSINE |
| Query prefix | "Represent this sentence for searching relevant passages:" |

The choice of BGE-small-en-v1.5 balances retrieval quality against memory footprint. At 384 dimensions, it is one-third the size of BGE-large (1024 dimensions) while retaining over 95% of retrieval accuracy on the MTEB benchmark. This is critical for a deployment model where the embedding model shares memory with four NIM inference services on a 128 GB device. The asymmetric query prefix, specific to the BGE family, improves retrieval precision by distinguishing query embeddings from document embeddings in the shared vector space.

---

## 4. Knowledge Augmentation

### 4.1 The 3-Domain Knowledge Graph

Vector retrieval alone cannot provide the structured clinical context that radiologists require for decision-making. A cosine similarity search for "lung nodule" may retrieve relevant literature but will not surface the specific Lung-RADS classification criteria, measurement thresholds, or follow-up intervals that drive clinical management. The 3-domain knowledge graph addresses this gap with curated, structured entries spanning three domains:

**Pathologies (~25 entries).** Each pathology entry contains:
- ICD-10 code and display name
- Applicable imaging modalities and body regions
- Subtypes (e.g., epidural, subdural, subarachnoid, intraparenchymal, intraventricular for intracranial hemorrhage)
- CT and MRI imaging characteristics (e.g., "Hyperdense (acute 50-70 HU), hypodense (chronic). Blood window 0-80 HU.")
- Severity criteria with quantitative thresholds (critical, urgent, routine)
- Key measurements (e.g., volume_ml, midline_shift_mm, max_thickness_mm)
- Clinical guidelines and applicable AI models
- Associated NIM workflow identifiers

**Modalities (~9 entries).** Covering CT, MRI, X-ray, CXR, PET, PET/CT, ultrasound, mammography, and fluoroscopy, each entry includes physics principles, typical protocols, clinical indications, and relative strengths and limitations.

**Anatomy (~21 entries).** Each anatomical structure includes SNOMED CT codes, Foundational Model of Anatomy (FMA) identifiers, imaging characteristics across modalities, common associated pathologies, and VISTA-3D segmentation label identifiers for direct integration with the 3D segmentation NIM.

### 4.2 Entity Alias Resolution

Clinical queries use inconsistent terminology. A query about "PE" could mean pulmonary embolism, pleural effusion, or pericardial effusion depending on context. The knowledge graph implements alias resolution that maps common abbreviations and synonyms to canonical pathology, modality, and anatomy entries. This ensures that queries using informal or abbreviated terminology still receive the correct structured knowledge augmentation.

### 4.3 Context Functions

Five context functions provide programmatic access to the knowledge graph:

- `get_pathology_context(name)` -- Returns structured pathology information including severity criteria and clinical guidelines
- `get_modality_context(name)` -- Returns modality characteristics and clinical indications
- `get_anatomy_context(name)` -- Returns anatomical references with ontology codes
- `get_nim_recommendation(query)` -- Analyzes a query and recommends applicable NIM services
- `get_comparison_context(entity_a, entity_b)` -- Returns structured comparison data for "X vs Y" queries

These functions are called automatically during the knowledge augmentation stage of the RAG pipeline, injecting relevant structured context alongside vector-retrieved documents.

---

## 5. Query Expansion

### 5.1 The Recall Problem

Medical imaging terminology is highly varied. A radiologist might describe the same finding as a "pulmonary nodule," "lung nodule," "SPN" (solitary pulmonary nodule), or simply "nodule." BGE-small-en-v1.5 captures some of this semantic similarity, but embedding models trained on general text corpora do not reliably map domain-specific abbreviations and jargon to their full forms. Without query expansion, a search for "ptx" (pneumothorax) would miss documents that use only the full term, and vice versa.

### 5.2 Expansion Architecture

The system implements 12 domain-specific expansion maps:

| Map | Coverage | Example |
|-----|----------|---------|
| Modality | 9 modalities | "ct" -> "computed tomography," "cat scan," "helical ct," "spiral ct," "mdct" |
| Body Region | 12 regions | "chest" -> "thorax," "thoracic," "lung," "pulmonary," "mediastinal," "pleural" |
| Pathology | 18+ pathologies | "pneumothorax" -> "ptx," "collapsed lung," "air leak," "tension pneumothorax" |
| AI Task | 8 tasks | "segmentation" -> "delineation," "contouring," "roi extraction" |
| Image Quality | 6 categories | "artifact" -> "motion artifact," "beam hardening," "metal artifact" |
| Contrast Phase | 5 phases | "arterial" -> "arterial phase," "early enhancement," "bolus timing" |
| Scoring System | 7 systems | "lung-rads" -> "Lung-RADS," "lung cancer screening," "ACR classification" |
| Severity | 4 levels | "critical" -> "emergent," "stat," "urgent," "priority 1" |
| Anatomy | 15+ structures | "liver" -> "hepatic," "hepatobiliary," "right lobe," "left lobe" |
| Device Category | 6 categories | "cad" -> "computer-aided detection," "AI triage," "automated detection" |
| Guideline Source | 5 sources | "acr" -> "American College of Radiology," "ACR Appropriateness" |
| Clinical Action | 6 actions | "follow-up" -> "surveillance," "interval imaging," "short-term follow-up" |

Expansion is applied before embedding: the expanded terms are concatenated with the original query to produce a richer embedding vector that has higher cosine similarity with relevant documents regardless of their specific terminology. Expansion adds less than 1 ms of latency and has no impact on search throughput.

---

## 6. NVIDIA NIM Integration

### 6.1 The 4 NIM Services

The Imaging Intelligence Agent integrates four NVIDIA NIM inference microservices, each running as an independent Docker container with an HTTP API:

| NIM Service | Port | Model | Purpose |
|------------|------|-------|---------|
| Llama-3 8B Instruct | 8520 | meta/llama3-8b-instruct | Clinical text synthesis, report generation, question answering |
| VISTA-3D | 8530 | nvidia/vista-3d | Zero-shot 3D segmentation across 132 anatomical classes |
| MAISI | 8531 | nvidia/maisi | High-resolution synthetic CT generation (up to 512x512x512) |
| VILA-M3 | 8532 | nvidia/vila-m3 | Vision-language medical image understanding |

**Llama-3 8B Instruct** serves as the fallback LLM when Claude (Anthropic API) is unavailable. It provides on-device text synthesis for RAG response generation, clinical report drafting, and multi-turn conversation. The model runs in INT8 quantization on the DGX Spark's GB10 GPU.

**VISTA-3D** (Versatile Imaging Segmentation and Triage Accelerator) provides zero-shot 3D segmentation across 132 anatomical classes without task-specific fine-tuning. In the agent, VISTA-3D supports automated organ segmentation for quantitative volume analysis, hemorrhage segmentation in the CT head hemorrhage workflow, and lesion segmentation in the MRI brain MS lesion workflow.

**MAISI** (Medical AI Synthetic Imaging) uses a latent diffusion model to generate high-resolution synthetic CT volumes with paired segmentation masks for up to 127 anatomical classes. Applications include training data augmentation for rare pathologies, algorithm validation without patient data, and educational demonstration of imaging findings.

**VILA-M3** provides natural language understanding of medical images, combining visual feature extraction with medical domain knowledge. The model accepts an image and a natural language question, returning findings with a confidence score. Applications include automated preliminary report drafting, interactive image-based question answering, and quality assurance for existing reports.

### 6.2 Fallback Chain

The NIM client layer implements a three-tier fallback chain for maximum availability:

```
Local NIM Container -> NVIDIA Cloud NIM -> Mock Mode
```

1. **Local NIM.** The primary path queries the on-device NIM container at its configured localhost port. Health is checked via the standard `/v1/health/ready` endpoint with 30-second cached polling intervals.
2. **NVIDIA Cloud NIM.** When local containers are unavailable and an NVIDIA API key is configured, the agent falls back to NVIDIA Cloud NIM endpoints at `integrate.api.nvidia.com` for LLM (meta/llama-3.1-8b-instruct) and VLM (meta/llama-3.2-11b-vision-instruct) inference.
3. **Mock Mode.** When both local and cloud paths are unavailable and `NIM_ALLOW_MOCK_FALLBACK=true`, the client returns clinically realistic synthetic responses. Mock mode enables demonstration, testing, and development without GPU hardware.

All NIM clients inherit from `BaseNIMClient`, which provides:
- Cached health status polling (30-second interval)
- Exponential-backoff retry via tenacity (3 attempts, 1-10 second wait)
- Unified error handling and structured logging
- Consistent status reporting across all four services

### 6.3 Mock Mode

Mock mode is a first-class design concern, not an afterthought. Each NIM client includes hand-crafted mock responses that are clinically realistic and internally consistent. For example, the VISTA-3D mock returns segmentation volumes and label counts consistent with typical anatomical proportions; the VILA-M3 mock returns findings and confidence scores that match the simulated pathology. This design enables three critical use cases:

1. **Demonstration without GPU.** The full system can be demonstrated on any laptop using Docker Compose in lite mode with mock NIM responses.
2. **Automated testing.** The 620-test suite runs entirely in mock mode, enabling CI/CD without GPU infrastructure.
3. **Development velocity.** Engineers can develop and test new workflows, UI features, and API endpoints without waiting for NIM container startup (60-120 seconds per service).

### 6.4 Memory Budget on DGX Spark

The 128 GB unified LPDDR5x memory of the DGX Spark enables simultaneous deployment of all four NIMs alongside the full application stack:

| Component | Estimated Memory |
|-----------|-----------------|
| Llama-3 8B Instruct | ~16 GB |
| VISTA-3D | ~8 GB |
| MAISI | ~12 GB |
| VILA-M3 | ~16 GB |
| Milvus + etcd + MinIO | ~4 GB |
| Application stack (FastAPI + Streamlit) | ~2 GB |
| Orthanc PACS + OHIF Viewer | ~1 GB |
| Embedding model (BGE-small) | ~0.5 GB |
| **Total** | **~59.5 GB** |
| **Available headroom** | **~68.5 GB** |

This deployment model is not feasible on typical workstation GPUs (24-48 GB VRAM) but fits comfortably within the DGX Spark's unified memory architecture, which shares the full 128 GB pool between CPU and GPU workloads via NVLink-C2C interconnect.

---

## 7. Clinical Workflows

### 7.1 The 6 Reference Implementations

The agent includes six reference clinical workflows, each implementing a complete analysis pipeline from image ingestion to structured report generation:

**1. CT Head Hemorrhage Triage (`ct_head_hemorrhage`).** An emergency radiology workflow for automated intracranial hemorrhage detection and triage. A 3D U-Net (MONAI) performs binary hemorrhage segmentation on CT head volumes preprocessed to RAS orientation, 1 mm isotropic spacing, and blood window (0-80 HU). Postprocessing computes hemorrhage volume (mL), midline shift (mm), and maximum thickness (mm). Urgency classification follows Brain Trauma Foundation thresholds: volume > 30 mL or midline shift > 5 mm or thickness > 10 mm triggers Critical (P1); volume > 5 mL triggers Urgent (P2). Target latency: < 90 seconds.

**2. CT Chest Lung Nodule Tracking (`ct_chest_lung_nodule`).** A lung cancer screening workflow implementing ACR Lung-RADS v2022 classification. RetinaNet (MONAI) detects candidate nodules, SegResNet performs per-nodule segmentation for volumetric analysis, and SyN diffeomorphic registration enables longitudinal tracking with volume doubling time (VDT) computation. Rule-based Lung-RADS assignment considers nodule morphology (solid, ground-glass, part-solid), size, and growth rate. Cross-modal trigger routes Lung-RADS 4A+ findings to genomic pipeline for EGFR/ALK/ROS1 queries. Target latency: < 5 minutes.

**3. CT Coronary Angiography (`ct_coronary_angiography`).** A cardiac imaging workflow implementing CAD-RADS classification for coronary artery disease assessment. The pipeline evaluates coronary artery stenosis, calcium scoring, and plaque characterization. CAD-RADS >= 3 triggers cross-modal queries for cardiovascular genomics (LDLR, PCSK9, APOB variants). Target latency: < 5 minutes.

**4. Chest X-Ray Rapid Findings (`cxr_rapid_findings`).** A triage workflow using DenseNet-121 for multi-label CXR classification. The model detects consolidation, effusion, pneumothorax, cardiomegaly, atelectasis, edema, and nodules with GradCAM heatmap localization for visual explanation. CXR critical findings trigger infection-related genomic queries (TLR4, MBL2). Target latency: < 30 seconds.

**5. MRI Brain MS Lesion Tracking (`mri_brain_ms_lesion`).** A multiple sclerosis monitoring workflow using 3D U-Net segmentation on FLAIR sequences. The pipeline includes SyN diffeomorphic registration for longitudinal comparison, lesion matching to classify new, enlarged, stable, and resolved lesions, and quantitative metrics (total lesion volume, lesion count, lesion change rate). Brain lesions with high activity trigger MS-related genomic queries. Target latency: < 5 minutes.

**6. MRI Prostate PI-RADS (`mri_prostate_pirads`).** A prostate cancer screening workflow implementing PI-RADS v2.1 scoring for multi-parametric MRI. The pipeline evaluates T2-weighted, DWI/ADC, and DCE sequences for transition zone and peripheral zone lesions. PI-RADS >= 4 triggers cross-modal queries for prostate cancer genomics (BRCA2, HOXB13). Target latency: < 5 minutes.

### 7.2 Workflow Implementation Pattern

All workflows inherit from `BaseImagingWorkflow` and implement a three-stage pipeline:

1. **Preprocess.** Load DICOM or NIfTI input, reorient to standard orientation, resample to isotropic spacing, apply modality-specific windowing and normalization.
2. **Infer.** Run model inference via MONAI or NIM (live mode) or return clinically realistic synthetic data (mock mode). Mock mode is the default and requires no GPU or medical image data.
3. **Postprocess.** Extract quantitative measurements, apply classification rules (Lung-RADS, CAD-RADS, PI-RADS, ASPECTS), determine severity, and generate structured findings.

Each workflow returns a `WorkflowResult` object containing the workflow identifier, findings, measurements, severity classification, classification score (e.g., "Lung-RADS 4A"), processing time, and an optional cross-modal trigger flag indicating whether genomic enrichment should be invoked.

### 7.3 Standardized Scoring Systems

The agent implements seven standardized scoring systems used across its clinical workflows:

| Scoring System | Domain | Categories |
|---------------|--------|------------|
| Lung-RADS v2022 | Lung cancer screening CT | 1, 2, 3, 4A, 4B, 4X, S |
| BI-RADS | Breast imaging | 0-6 |
| TI-RADS | Thyroid ultrasound | TR1-TR5 |
| LI-RADS | Liver imaging | LR-1 to LR-5, LR-M, LR-TIV |
| CAD-RADS | Coronary CT angiography | 0-5, N |
| PI-RADS | Prostate MRI | 1-5 |
| ASPECTS | Acute stroke CT | 0-10 |

These scoring systems are implemented as rule-based classifiers operating on quantitative measurements from the inference stage. Each classifier maps measurements to the appropriate category, determines the recommended clinical action, and identifies whether the finding meets the threshold for cross-modal genomic enrichment.

### 7.4 Cross-Modal Genomic Integration

The cross-modal trigger mechanism is the architectural bridge between the Imaging Intelligence Agent and the HCLS AI Factory's genomics pipeline. When an imaging workflow produces findings that exceed clinical severity thresholds, the system automatically queries the genomic evidence collection (3.5+ million vectors) for relevant molecular context:

| Imaging Trigger | Threshold | Genomic Queries |
|----------------|-----------|-----------------|
| Lung nodule | Lung-RADS 4A+ | EGFR, ALK, ROS1, KRAS driver mutations |
| CXR critical finding | Consolidation (critical/urgent) | TLR4, MBL2 infection susceptibility |
| Brain lesion | High activity score | HLA-DRB1, MS susceptibility genes |
| Coronary artery disease | CAD-RADS >= 3 | LDLR, PCSK9, APOB cardiovascular variants |
| Prostate lesion | PI-RADS >= 4 | BRCA2, HOXB13 cancer susceptibility |

This integration enables a precision medicine workflow entirely on one device: an imaging study reveals a suspicious lung nodule, the workflow classifies it as Lung-RADS 4A, the cross-modal trigger queries the genomic evidence collection for EGFR/ALK/ROS1 variant evidence, and the synthesized report includes both imaging findings and relevant genomic context -- without any data leaving the DGX Spark.

### 7.5 Demo Cases

Four pre-configured demo cases enable immediate demonstration:

| Case ID | Pathology | Patient | Workflow |
|---------|-----------|---------|----------|
| DEMO-001 | Intracranial Hemorrhage | 62-year-old male | ct_head_hemorrhage |
| DEMO-002 | Lung Nodule | 58-year-old female | ct_chest_lung_nodule |
| DEMO-003 | Coronary Artery Disease | 55-year-old male | ct_coronary_angiography |
| DEMO-004 | Pneumonia | 45-year-old female | cxr_rapid_findings |

Each demo case includes pre-configured clinical context, synthetic imaging data, and expected workflow outputs, enabling complete end-to-end demonstration without real patient data.

---

## 8. Evidence Quality and Citation Provenance

### 8.1 Citation Scoring

Every evidence item returned by the RAG engine includes a cosine similarity score, a collection source identifier, and citation metadata (PubMed ID, ClinicalTrials.gov NCT number, or seed data reference). The system classifies citation relevance into three tiers:

| Tier | Score Range | Interpretation |
|------|-----------|----------------|
| High | >= 0.75 | Directly relevant to the query |
| Medium | 0.60 - 0.75 | Contextually relevant, may address a related aspect |
| Low | 0.40 - 0.60 | Tangentially relevant, included for completeness |

Results below the 0.40 threshold are discarded. This three-tier system helps clinicians quickly assess the strength of the evidence base supporting each synthesized answer.

### 8.2 Source Provenance

For literature results, PubMed IDs are rendered as clickable links to the original abstract. For clinical trial results, NCT numbers link to ClinicalTrials.gov records. For guideline results, source organizations (ACR, RSNA, ESR) and publication years are included. This provenance chain enables clinicians to verify any claim made by the system against the original source, maintaining the traceability that evidence-based medicine demands.

### 8.3 Export Formats

The agent supports four export formats for clinical and research use:

- **Markdown.** Human-readable reports with structured sections, findings tables, and citation lists. Suitable for documentation and collaboration.
- **JSON.** Machine-readable structured output containing all findings, measurements, scores, and metadata. Suitable for integration with downstream systems.
- **PDF (ReportLab).** Publication-quality reports with institutional headers, formatted tables, and embedded visualizations. Suitable for clinical documentation and archival.
- **FHIR R4.** Interoperable diagnostic reports conforming to the HL7 FHIR R4 specification, containing Patient, ImagingStudy, Observation, and DiagnosticReport resources. Findings are coded with SNOMED CT (10+ imaging finding codes), categories use LOINC (LP29684-5 Radiology), and modalities use standard DICOM codes. FHIR bundles can be submitted directly to compliant EHR systems.

---

## 9. Hardware Democratization

### 9.1 The DGX Spark

The NVIDIA DGX Spark represents a paradigm shift in AI infrastructure accessibility:

| Component | Specification |
|-----------|--------------|
| GPU | NVIDIA GB10 (Blackwell architecture) |
| Memory | 128 GB unified LPDDR5x |
| CPU | 20 ARM cores (Grace architecture) |
| Interconnect | NVLink-C2C (CPU-GPU) |
| Form factor | Desktop (comparable to Mac Studio) |
| Price | $3,999 |

The unified memory architecture is the critical enabler: unlike discrete GPU systems where VRAM and system RAM are separate pools connected by PCIe, the DGX Spark shares its full 128 GB between CPU and GPU workloads via NVLink-C2C. This means the four NIM inference services (~52 GB), Milvus vector database (~4 GB), and application stack (~3.5 GB) coexist in a single memory pool without the data transfer bottlenecks that characterize PCIe-connected systems.

### 9.2 Local vs Cloud

The deployment model prioritizes data sovereignty:

| Attribute | Local (DGX Spark) | Cloud |
|-----------|-------------------|-------|
| Data sovereignty | Complete -- no data leaves device | Dependent on cloud provider |
| Latency | Deterministic, no network variability | Variable, network-dependent |
| Recurring cost | None (one-time $3,999) | $500-5,000+/month |
| Compliance | HIPAA/GDPR straightforward | Requires BAA, data residency |
| Availability | 100% uptime (local) | Dependent on cloud availability |
| Scalability | Fixed (1 device) | Elastic |

For organizations processing fewer than 1,000 imaging studies per day -- which includes the vast majority of community hospitals, outpatient imaging centers, and research institutions -- the DGX Spark provides sufficient throughput at a fraction of the total cost of ownership of cloud-based solutions.

### 9.3 Deployment Models

**Full Stack (11 services).** Docker Compose deploys the complete system including Milvus (with etcd and MinIO), FastAPI server, Streamlit UI, four NIM containers, Orthanc PACS, and OHIF viewer. Requires DGX Spark or equivalent GPU system with 128+ GB memory.

**Lite Stack (6 services).** Deploys the RAG knowledge system (Milvus + application stack) without NIM containers. NIM-dependent features operate in mock mode. Requires no GPU and runs on any machine with Docker, making it accessible for development, testing, and demonstration on standard hardware.

---

## 10. Results and Capabilities

### 10.1 Performance Metrics

**RAG Performance:**

| Metric | Value |
|--------|-------|
| Vector search (11 collections, top-5 each) | 12-16 ms (warm cache) |
| Query expansion (12 maps) | < 1 ms |
| Knowledge graph augmentation | < 1 ms |
| Full RAG query (search + Claude synthesis) | ~20-30 seconds |
| Comparative RAG query (dual retrieval) | ~25-35 seconds |
| Cosine similarity scores (typical) | 0.70-0.90 |

**Workflow Target Latencies:**

| Workflow | Target Latency | Mock Latency |
|----------|---------------|-------------|
| CT Head Hemorrhage Triage | < 90 seconds | < 100 ms |
| CT Chest Lung Nodule Tracking | < 5 minutes | < 100 ms |
| CT Coronary Angiography | < 5 minutes | < 100 ms |
| CXR Rapid Findings | < 30 seconds | < 100 ms |
| MRI Brain MS Lesion Tracking | < 5 minutes | < 100 ms |
| MRI Prostate PI-RADS | < 5 minutes | < 100 ms |

**NIM Startup Times:**

| NIM Service | Model Load Time |
|------------|----------------|
| Llama-3 8B Instruct | 60-90 seconds |
| VISTA-3D | 60-90 seconds |
| MAISI | 90-120 seconds |
| VILA-M3 | 60-90 seconds |

### 10.2 Test Coverage

The system includes 620 tests across 12 modules:

| Module | Scope |
|--------|-------|
| Models | Pydantic data models for all collections and workflow results |
| Collections | Milvus collection creation, schema validation, CRUD operations |
| Knowledge | Knowledge graph lookup, alias resolution, comparison context |
| Query Expansion | Expansion map coverage, term deduplication, edge cases |
| RAG Engine | Multi-collection search, weighted scoring, LLM synthesis |
| NIM Clients | Health checking, retry logic, mock fallback, all four services |
| Workflows | All six workflow pipelines in mock mode |
| Cross-Modal | Trigger rules, genomic query generation, threshold evaluation |
| Export | Markdown, JSON, PDF, and FHIR R4 generation |
| API | All FastAPI endpoints across 8 routers |
| DICOM | Orthanc integration, DICOM parsing, study routing |
| Preview | Image preview generation and caching |

All 620 tests pass in mock mode without GPU or external service dependencies.

### 10.3 Query Capabilities

The system handles five categories of clinical queries:

1. **Evidence queries.** "What is the evidence for AI-assisted lung nodule detection in low-dose CT screening?" -- searches literature, trials, benchmarks, and guidelines with knowledge graph augmentation.
2. **Protocol queries.** "What is the optimal CT protocol for pulmonary embolism?" -- prioritizes protocol and guideline collections with modality-specific expansion.
3. **Device queries.** "Which FDA-cleared AI devices are available for chest X-ray triage?" -- searches device collection with regulatory metadata filtering.
4. **Comparative queries.** "Compare CT vs MRI for brain hemorrhage detection" -- auto-detected "X vs Y" pattern triggers dual retrieval and structured comparison output.
5. **Workflow queries.** "Analyze this CT head for hemorrhage" -- routes to the appropriate clinical workflow, executes inference, applies scoring, and optionally triggers cross-modal genomic enrichment.

---

## 11. Integration with the HCLS AI Factory

### 11.1 The 3-Stage Pipeline

The Imaging Intelligence Agent operates as a specialized intelligence module within the HCLS AI Factory, an end-to-end precision medicine platform that processes patient DNA to drug candidates in under 5 hours:

**Stage 1: Genomics Pipeline.** Parabricks, DeepVariant, and BWA-MEM2 process raw FASTQ sequencing data into annotated VCF files. GPU-accelerated alignment and variant calling reduce processing time from 24-48 hours (CPU) to 120-240 minutes.

**Stage 2: RAG/Chat Pipeline.** Milvus stores 3.56 million searchable vectors from ClinVar (4.1M records), AlphaMissense (71M predictions), and curated annotations. Claude AI synthesizes variant interpretations with evidence citations. This is the pipeline that populates the shared `genomic_evidence` collection.

**Stage 3: Drug Discovery Pipeline.** BioNeMo MolMIM, DiffDock, and RDKit generate and evaluate drug candidates targeting identified molecular variants. The pipeline covers 201 genes across 13 therapeutic areas with 171 druggable targets.

The Imaging Intelligence Agent connects to this pipeline at two points: it reads from the Stage 2 genomic evidence collection for cross-modal enrichment, and it can trigger Stage 1 genomic analysis when imaging findings suggest the need for molecular profiling (e.g., a Lung-RADS 4B lung nodule prompting tumor genomic profiling).

### 11.2 The Genomic Evidence Bridge

The cross-modal bridge is implemented through the shared `genomic_evidence` Milvus collection. This collection, populated by the Stage 2 RAG pipeline with 3.5+ million vectors from ClinVar, AlphaMissense, and curated variant annotations, is mounted as a read-only collection in the Imaging Intelligence Agent. When a cross-modal trigger fires, the agent:

1. Generates domain-specific genomic query strings (e.g., "lung cancer driver mutations EGFR ALK ROS1 KRAS")
2. Embeds these queries using the shared BGE-small-en-v1.5 model
3. Searches the genomic_evidence collection with the same IVF_FLAT/COSINE configuration
4. Returns a `CrossModalResult` containing the genomic evidence alongside the imaging findings

This architecture avoids data duplication: the genomic evidence is maintained by the Stage 2 pipeline and consumed read-only by the imaging agent. Updates to the genomic evidence (new ClinVar releases, new AlphaMissense predictions) are automatically available to the imaging agent without any re-indexing.

### 11.3 The Five Intelligence Agents

The HCLS AI Factory includes five specialized intelligence agents, each addressing a distinct clinical domain:

1. **Biomarker Intelligence Agent** -- Biomarker discovery and validation
2. **Oncology Intelligence Agent** -- Cancer genomics and treatment selection
3. **CAR-T Intelligence Agent** -- CAR-T cell therapy development lifecycle
4. **Imaging Intelligence Agent** -- Medical imaging AI and clinical decision support
5. **Autoimmune Intelligence Agent** -- Autoimmune disease genomics and therapy

All five agents share the same architectural pattern (multi-collection RAG, knowledge graph augmentation, query expansion, Claude synthesis) and the same infrastructure (Milvus, BGE-small-en-v1.5, DGX Spark). The `genomic_evidence` collection serves as a shared knowledge substrate across agents, enabling cross-domain intelligence without data silos.

---

## 12. Discussion

### 12.1 Implications for Clinical Radiology

The Imaging Intelligence Agent demonstrates that a comprehensive medical imaging AI platform -- encompassing evidence retrieval, clinical decision support, inference orchestration, and structured reporting -- can deploy on a single $3,999 device. This has three significant implications:

First, it removes the infrastructure barrier that has confined medical imaging AI to well-resourced institutions. A community hospital can deploy the complete stack with a single hardware purchase and no recurring cloud costs, gaining access to the same AI capabilities that previously required six-figure infrastructure investments.

Second, the cross-modal genomic integration demonstrates that precision medicine workflows can execute entirely on a local device. When imaging findings trigger genomic queries, the entire chain -- from image analysis through molecular context enrichment to report generation -- happens without any data leaving the device. This is particularly significant for regulated environments where data sovereignty concerns have been a persistent barrier to AI adoption.

Third, the multi-collection RAG architecture with knowledge graph augmentation provides a model for building domain-specific AI systems that are grounded in evidence rather than parametric knowledge alone. Every answer is traceable to specific documents, guidelines, or structured knowledge entries, maintaining the evidentiary standards that clinical practice demands.

### 12.2 Democratization Beyond Radiology

The architectural pattern demonstrated here -- multi-collection vector retrieval, domain knowledge graph, query expansion, NIM inference, cross-modal integration -- is applicable to any clinical domain where evidence synthesis is needed. The same approach could support pathology AI (whole slide image analysis with molecular context), cardiology AI (ECG interpretation with genomic risk scoring), or ophthalmology AI (retinal imaging with systemic disease correlation). The HCLS AI Factory's five intelligence agents already demonstrate this generalizability across biomarkers, oncology, CAR-T therapy, imaging, and autoimmune disease.

### 12.3 Limitations

Several limitations should be acknowledged:

1. **Research prototype status.** The system is a research prototype demonstrating architectural feasibility, not a cleared medical device. Clinical deployment would require FDA 510(k) or De Novo clearance, clinical validation studies, and institutional review.
2. **Mock inference.** The current public release operates primarily in mock mode, with synthetic inference results. Real inference requires NVIDIA NIM container access (NGC account) and DGX Spark hardware. Mock results are clinically realistic but not clinically validated.
3. **Single-device throughput.** The DGX Spark is designed for low-to-moderate throughput. High-volume imaging centers processing thousands of studies per day would require multi-device deployment or cloud augmentation.
4. **Embedding model limitations.** BGE-small-en-v1.5, while efficient, is a general-purpose embedding model. A domain-specific model fine-tuned on radiology text (e.g., RadBERT embeddings) could improve retrieval precision for highly specialized queries.
5. **Knowledge currency.** The seed data collections contain curated content that requires periodic updates to reflect evolving guidelines, new device clearances, and emerging evidence. The PubMed and ClinicalTrials.gov ingest pipelines provide automated updates for literature and trials, but other collections depend on manual curation.
6. **Cross-modal validation.** The genomic trigger thresholds (e.g., Lung-RADS 4A+ for lung cancer genomics) are based on clinical guidelines but have not been validated in a prospective clinical study to assess their impact on patient outcomes.

---

## 13. Conclusion

### 13.1 Key Contributions

The Imaging Intelligence Agent makes five contributions to the field of medical imaging AI:

1. **Multi-collection RAG for radiology.** A 10+1 collection architecture with parallel search, weighted scoring, and knowledge graph augmentation that synthesizes evidence across literature, trials, findings, protocols, devices, anatomy, benchmarks, guidelines, report templates, datasets, and genomic evidence into unified clinical intelligence.

2. **Cross-modal precision medicine.** An automated bridge between imaging findings and genomic evidence that triggers molecular context enrichment when imaging severity thresholds are exceeded, enabling precision medicine workflows entirely on a local device.

3. **Reference clinical workflows.** Six end-to-end analysis pipelines implementing standardized scoring systems (Lung-RADS, CAD-RADS, PI-RADS, ASPECTS, BI-RADS, TI-RADS, LI-RADS) with graceful degradation from live inference to mock mode.

4. **Hardware democratization.** A complete deployment architecture -- four NIM services, vector database, PACS server, web viewer, and application stack -- on a single NVIDIA DGX Spark ($3,999) with 128 GB unified memory, requiring no cloud dependency.

5. **Open-source implementation.** 620 tests across 12 modules, 19+ API endpoints across 8 routers, 9 Streamlit UI tabs, and 6 clinical workflows released under Apache 2.0 with comprehensive documentation.

### 13.2 Future Directions

**DICOM Structured Reporting.** Generation of DICOM SR (Structured Report) objects via highdicom TID 1500 measurement reports, enabling AI findings to be stored alongside source images in the DICOM archive.

**LangGraph Multi-Step Agent.** An advanced reasoning agent that chains triage, longitudinal analysis, population comparison, and report generation in a multi-step workflow with human-in-the-loop checkpoints.

**Federated Learning.** Three NVIDIA FLARE job configurations are in place (CXR classification with DenseNet-121, CT organ segmentation with SegResNet, lung nodule detection with RetinaNet). Next steps include multi-institution pilot deployment and integration of differential privacy mechanisms.

**Population Analytics.** Cohort-level imaging trends, disease prevalence monitoring, and outcomes tracking across institutional imaging archives.

**Domain-Specific Embeddings.** Fine-tuning the embedding model on radiology-specific corpora (RadLex, ACR guidelines, radiology reports) to improve retrieval precision for specialized terminology.

**Multi-Institutional Validation.** Prospective evaluation of the cross-modal trigger mechanism to assess whether automated genomic enrichment of imaging findings improves clinical decision-making and patient outcomes.

---

## 14. References

1. NVIDIA DGX Spark Technical Specifications, NVIDIA Corporation, 2025.
2. ACR Lung-RADS v2022: Assessment Categories and Management Recommendations, American College of Radiology, 2022.
3. Brain Trauma Foundation Guidelines for the Management of Severe Traumatic Brain Injury, 4th Edition, 2016.
4. Wang J, Yi X, Guo R, et al. Milvus: A Purpose-Built Vector Data Management System. *Proceedings of the ACM SIGMOD International Conference on Management of Data*, 2021.
5. Xiao S, Liu Z, Zhang P, Muennighoff N. C-Pack: Packaged Resources to Advance General Chinese Embedding (BGE). *arXiv preprint arXiv:2309.07597*, 2023.
6. Project MONAI Consortium. MONAI: Medical Open Network for Artificial Intelligence. https://monai.io, 2024.
7. NVIDIA. VISTA-3D: Versatile Imaging Segmentation and Triage Accelerator. NVIDIA NGC Catalog, 2024.
8. NVIDIA. MAISI: Medical AI Synthetic Imaging. NVIDIA NGC Catalog, 2024.
9. NVIDIA. VILA-M3: Vision-Language Model for Medical Imaging. NVIDIA NGC Catalog, 2024.
10. Touvron H, et al. Llama 2: Open Foundation and Fine-Tuned Chat Models. *arXiv preprint arXiv:2307.09288*, 2023.
11. Irvin J, et al. CheXpert: A Large Chest Radiograph Dataset with Uncertainty Labels and Expert Comparison. *Proceedings of the AAAI Conference on Artificial Intelligence*, 2019.
12. Rajpurkar P, et al. CheXNet: Radiologist-Level Pneumonia Detection on Chest X-Rays with Deep Learning. *arXiv preprint arXiv:1711.05225*, 2017.
13. HL7 International. FHIR R4 Specification. https://hl7.org/fhir/R4, 2019.
14. ACR. BI-RADS Atlas, 5th Edition. American College of Radiology, 2013.
15. Tessler FN, et al. ACR TI-RADS: American College of Radiology Thyroid Imaging, Reporting, and Data System. *Journal of the American College of Radiology*, 14(5):587-595, 2017.
16. Cury RC, et al. CAD-RADS 2.0: Coronary Artery Disease Reporting and Data System. *Radiology: Cardiothoracic Imaging*, 4(5), 2022.
17. Turkbey B, et al. PI-RADS Prostate Imaging Reporting and Data System: 2019, Version 2.1. *European Urology*, 76(3):340-351, 2019.
18. Barber PA, et al. Validity and Reliability of a Quantitative Computed Tomography Score in Predicting Outcome of Hyperacute Stroke Before Thrombolytic Therapy (ASPECTS). *The Lancet*, 355(9216):1670-1674, 2000.
19. Mitchell TM, et al. LI-RADS: Liver Imaging Reporting and Data System. *American College of Radiology*, 2018.
20. Clark K, et al. The Cancer Imaging Archive (TCIA): Maintaining and Operating a Public Information Repository. *Journal of Digital Imaging*, 26(6):1045-1057, 2013.

---

*Apache 2.0 License. HCLS AI Factory -- Imaging Intelligence Agent.*
*Part of the HCLS AI Factory precision medicine platform: Patient DNA to Drug Candidates in < 5 hours.*
