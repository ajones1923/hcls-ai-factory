# Learning Guide -- Foundations

**Imaging Intelligence Agent | HCLS AI Factory**

Author: Adam Jones
Date: March 2026
License: Apache 2.0

---

## Welcome

You are reading the foundational learning guide for the Imaging Intelligence Agent, an AI-powered research assistant that unifies scattered knowledge across the medical imaging landscape -- from peer-reviewed literature and clinical trials to device clearances, imaging protocols, scoring guidelines, and AI benchmarks. This system is part of the HCLS AI Factory, an end-to-end precision medicine platform that runs on a single NVIDIA DGX Spark ($4,699).

### Who this guide is for

This guide is written for three audiences:

- **Radiologists and clinicians** who know medical imaging but are new to AI/ML, vector databases, and retrieval-augmented generation.
- **Data scientists and ML engineers** who know embeddings and LLMs but are new to radiology, DICOM, imaging physics, and clinical scoring systems.
- **Software developers** who want to understand how the system works end to end -- from the Streamlit UI with its 9 specialized tabs, through Milvus vector search, to NIM-powered inference with VISTA-3D and MAISI.

You do not need to be an expert in all three areas. The whole point of this guide is to bring you up to speed on whichever parts are new to you.

### What you will learn

By the end of this guide, you will understand:

1. What medical imaging AI is and why it matters
2. Why data fragmentation is the central challenge in imaging AI
3. How Retrieval-Augmented Generation (RAG) works, from first principles
4. How this system searches 11 collections simultaneously
5. How to use the UI to ask questions, run workflows, and interpret results
6. What each of the 11 data collections contains and why it exists
7. How the knowledge graph covers 3 domains: pathologies, modalities, and anatomy
8. How query expansion improves search recall with synonyms and aliases
9. How to set up and run the system locally
10. How to use the REST API with its 19+ endpoints

### Prerequisites

- Basic Python knowledge (you can read a Python function and understand what it does)
- A computer with a terminal
- Curiosity about either medical imaging or AI-assisted research (or both)

No prior knowledge of radiology, DICOM, vector databases, or large language models is required. We will build every concept from the ground up.

---

## Chapter 1: What Is Medical Imaging AI?

### The core idea

Medical imaging is the practice of creating visual representations of the interior of the human body for clinical analysis and medical intervention. Every day, millions of CT scans, MRIs, X-rays, ultrasounds, and other studies are performed worldwide. Radiologists -- physicians who specialize in interpreting these images -- read each study, identify abnormalities, and write a report describing their findings.

The problem is scale. A single radiologist may read 50 to 100 studies per day, each containing dozens to hundreds of images. The volume of imaging studies grows 5-10% annually, outpacing the supply of trained radiologists. Errors of perception (missing a finding) and errors of interpretation (misclassifying a finding) are inevitable under this workload.

**Medical imaging AI** uses deep learning and computer vision to assist radiologists by automatically detecting, measuring, and classifying findings in medical images. The goal is not to replace radiologists, but to augment their capabilities -- catching findings they might miss, reducing measurement variability, and prioritizing urgent cases.

### The nine imaging modalities

This system covers nine primary imaging modalities:

| Modality | Full Name | Physics | Primary Use |
|----------|-----------|---------|-------------|
| **CT** | Computed Tomography | X-ray attenuation (Hounsfield units) | Trauma, cancer staging, pulmonary embolism, stroke |
| **MRI** | Magnetic Resonance Imaging | Hydrogen proton spin (T1, T2, DWI, FLAIR) | Brain, spine, joints, soft tissue, cardiac |
| **X-ray** | Radiography | X-ray transmission | Fractures, pneumonia, chest screening |
| **CXR** | Chest X-ray | X-ray transmission (chest-specific) | Pneumonia, pneumothorax, pleural effusion, cardiomegaly |
| **PET** | Positron Emission Tomography | Radiotracer uptake (FDG metabolism) | Oncology staging, cardiac viability, neurodegeneration |
| **PET/CT** | Combined PET-CT | PET + CT co-registration | Cancer staging with anatomical localization |
| **Ultrasound** | Ultrasonography | Sound wave reflection | OB/GYN, vascular, thyroid, breast, cardiac |
| **Mammography** | Breast Imaging | Low-dose X-ray | Breast cancer screening, tomosynthesis (3D) |
| **Fluoroscopy** | Fluoroscopic Imaging | Real-time X-ray | GI studies, catheter placement, joint injection |

Each modality produces different types of images, with different spatial resolutions, contrast characteristics, and clinical applications. A CT scan of the chest might contain 300-500 axial slices at sub-millimeter resolution, while a chest X-ray is a single 2D projection. An MRI brain study might include six different contrast sequences (T1, T2, FLAIR, DWI, ADC, post-contrast T1), each highlighting different tissue properties.

### What AI can do with medical images

AI applications in medical imaging fall into several categories:

```
DETECTION
  Finding abnormalities that might be missed by the human eye.
  Example: Detecting a 4 mm lung nodule on a chest CT among 400 slices.

CLASSIFICATION
  Categorizing findings according to standardized scoring systems.
  Example: Assigning a Lung-RADS 4A score to a part-solid nodule.

SEGMENTATION
  Outlining the precise boundaries of anatomical structures or lesions.
  Example: Measuring tumor volume by segmenting a liver lesion in 3D.

MEASUREMENT
  Quantifying characteristics like size, volume, density, or change over time.
  Example: Calculating coronary artery calcium score from a non-contrast CT.

TRIAGE
  Prioritizing urgent cases for immediate radiologist review.
  Example: Flagging a CT head with suspected intracranial hemorrhage.

REPORT GENERATION
  Drafting structured radiology reports from image findings.
  Example: Generating a Lung-RADS-compliant CT lung screening report.
```

### The FDA landscape

As of 2026, the U.S. FDA has cleared over 900 AI/ML-enabled medical devices. Radiology leads all medical specialties in FDA-cleared AI products, with roughly 75% of all cleared devices. Key regulatory pathways include:

- **510(k) clearance:** The most common pathway. The device must demonstrate substantial equivalence to an existing (predicate) device. Most radiology AI products use this pathway.
- **De Novo classification:** For novel devices without a predicate. Creates a new regulatory classification that future devices can reference.
- **PMA (Premarket Approval):** The most rigorous pathway, required for Class III devices. Requires clinical evidence of safety and effectiveness.
- **Breakthrough Device Designation:** An expedited pathway for devices that provide more effective treatment or diagnosis of life-threatening conditions.

The rapid growth in FDA-cleared imaging AI devices creates its own data challenge: how does a clinician or researcher keep track of which devices are available, what modalities and body regions they cover, how they perform on benchmark datasets, and what the current evidence says about their clinical impact?

That is one of the problems this system solves.

### The seven scoring systems

Standardized scoring systems are the backbone of structured radiology reporting. This system covers seven major scoring systems:

| System | Full Name | Modality | Body Region | Scale |
|--------|-----------|----------|-------------|-------|
| **Lung-RADS** | Lung Imaging Reporting and Data System | CT | Chest | 0, 1, 2, 3, 4A, 4B, 4X |
| **BI-RADS** | Breast Imaging Reporting and Data System | Mammography, US, MRI | Breast | 0-6 |
| **TI-RADS** | Thyroid Imaging Reporting and Data System | Ultrasound | Neck | TR1-TR5 |
| **LI-RADS** | Liver Imaging Reporting and Data System | CT, MRI | Abdomen | LR-1 to LR-5, LR-M, LR-TIV |
| **CAD-RADS** | Coronary Artery Disease Reporting and Data System | CT angiography | Cardiac | 0-5 |
| **PI-RADS** | Prostate Imaging Reporting and Data System | MRI | Pelvis | 1-5 |
| **ASPECTS** | Alberta Stroke Program Early CT Score | CT | Head | 0-10 |

Each scoring system standardizes how findings are described and communicated. For example, a Lung-RADS 4B nodule on chest CT means the patient needs tissue sampling or PET/CT, while a Lung-RADS 2 nodule means routine annual screening. These scores directly influence clinical decisions, and the ability to query across scoring systems, guidelines, and evidence is a core capability of this system.

---

## Chapter 2: The Data Challenge

### Where imaging AI knowledge lives

A radiologist evaluating a new AI tool, or a researcher developing one, needs to consider:

- **Published literature:** Tens of thousands of papers on PubMed describe imaging AI methods, validation studies, and clinical outcomes.
- **Clinical trials:** Hundreds of imaging AI trials are registered on ClinicalTrials.gov, spanning multiple modalities, body regions, and AI tasks.
- **FDA-cleared devices:** Over 900 FDA-cleared AI/ML medical devices, each with a specific intended use, cleared modality, and performance data.
- **Imaging protocols:** Standardized acquisition parameters (slice thickness, contrast agents, radiation dose) that vary by clinical indication and institution.
- **Clinical guidelines:** ACR Appropriateness Criteria, Fleischner Society guidelines, NCCN imaging recommendations -- each updated on different schedules by different organizations.
- **Benchmark datasets:** RSNA competitions, MICCAI challenges, NIH Clinical Center releases -- each with specific modalities, annotation types, and licensing terms.
- **Anatomy atlases:** Standardized anatomical terminology (SNOMED CT, FMA), segmentation label maps, and cross-sectional imaging characteristics.
- **Radiology findings:** Structured findings with severity levels, differential diagnoses, measurement types, and recommended follow-up actions.
- **Report templates:** Structured reporting templates aligned with ACR guidelines and specific coding systems (CPT, ICD-10, RadLex).
- **Genomic data:** Variant-level data from whole-genome sequencing that can reveal imaging-genomic correlations (radiogenomics).

### The problem: data silos

Each of these data sources lives in a different database, uses a different format, and is searched with a different interface. A researcher who wants to answer a cross-functional question -- for example, "What is the best AI model for detecting intracranial hemorrhage on non-contrast CT, and what are the current clinical guidelines for managing it?" -- would need to:

1. Search PubMed for AI performance studies on ICH detection
2. Search ClinicalTrials.gov for active trials validating ICH detection algorithms
3. Check the FDA device database for cleared ICH detection products
4. Review ACR and AHA/ASA guidelines for ICH management
5. Look up benchmark results from the RSNA 2019 intracranial hemorrhage challenge
6. Find the optimal CT acquisition protocol for non-contrast head imaging
7. Verify the anatomical regions involved and their segmentation labels

This is slow, error-prone, and nearly impossible to do comprehensively. A researcher might spend days manually searching across these sources and still miss a recently cleared device or an updated guideline recommendation.

### The solution: one unified intelligence platform

The Imaging Intelligence Agent solves this by:

1. **Ingesting** data from all of these sources into a single system
2. **Embedding** every piece of text as a 384-dimensional vector (a list of 384 numbers that captures the meaning of the text)
3. **Storing** these vectors in 11 purpose-built Milvus collections
4. **Searching** all 11 collections simultaneously when you ask a question
5. **Augmenting** the search results with structured knowledge from a hand-curated knowledge graph covering 25 pathologies, 9 modalities, and 21 anatomy structures
6. **Synthesizing** a grounded answer using an LLM, with citations back to the original sources

The result: you ask one question, and the system searches all 11 collections in parallel, augments the results with domain knowledge and query expansion, and produces a cross-functional answer that cites published papers, clinical trials, FDA-cleared devices, guidelines, protocols, and benchmark data -- all in one response.

---

## Chapter 3: What Is RAG?

### The limitation of language models

Large Language Models (LLMs) like Llama 3 or Claude are trained on vast amounts of text. They can write coherently, reason about complex topics, and follow instructions. But they have a fundamental limitation: **they do not have access to your specific data.**

If you ask an LLM about the best AI model for detecting lung nodules on CT, it can give you a general answer based on its training data. But it cannot tell you what is in your Milvus database. It cannot cite a specific benchmark result. It cannot look up whether a particular device has received FDA clearance since its training cutoff.

RAG solves this.

### RAG = Retrieval-Augmented Generation

RAG is a three-step pattern:

```
1. RETRIEVAL
   You ask a question. The system finds the most relevant documents
   from your database.

2. AUGMENTATION
   The retrieved documents are added to the prompt as context.
   The LLM now has your specific data in front of it.

3. GENERATION
   The LLM reads the evidence and generates an answer that is
   grounded in your data, with citations.
```

Think of it like this: imagine you are taking an open-book exam. The LLM is the student. RAG is the process of finding the right pages in the textbook (retrieval), putting them on the desk in front of the student (augmentation), and then asking the student to write an answer using those pages (generation).

Without RAG, the student is answering from memory (which may be outdated or vague). With RAG, the student is answering from evidence (which is specific, current, and citable).

### How retrieval works: embeddings and vector similarity

The retrieval step is the most technically interesting part, so let us break it down.

#### What is an embedding?

An embedding is a way of representing text as a list of numbers (a "vector") such that texts with similar meanings have similar numbers.

This system uses a model called **BGE-small-en-v1.5**, which converts any piece of text into a vector of **384 numbers** (384 dimensions). For example:

```
"Intracranial hemorrhage detection on non-contrast CT"
    --> [0.023, -0.156, 0.891, 0.044, ..., -0.312]  (384 numbers)

"CT head for acute brain bleeding without contrast"
    --> [0.019, -0.148, 0.883, 0.051, ..., -0.298]  (384 numbers)
```

These two texts have very similar vectors because they express similar meanings, even though the words are different.

Conversely:

```
"MRI prostate protocol with multiparametric sequences"
    --> [0.512, 0.078, -0.234, 0.667, ..., 0.112]  (384 numbers)
```

This vector looks very different from the hemorrhage vectors because the topic is different.

#### An analogy: GPS coordinates for meaning

Think of embeddings like GPS coordinates for meaning. Just as GPS coordinates place a physical location in a two-dimensional space (latitude and longitude), embeddings place a piece of text in a 384-dimensional meaning space. Texts about similar topics end up at nearby coordinates. Texts about unrelated topics end up far apart.

You cannot visualize 384 dimensions (nobody can), but the math works the same way as it does in two dimensions. To find texts similar to your question, you measure the "distance" between vectors.

#### How similarity search works

When you type a question into the Imaging Intelligence Agent:

1. Your question is embedded into a 384-dimensional vector.
2. That vector is compared to every vector in the database using **cosine similarity** (a measure of how close two vectors are in direction, regardless of length).
3. The vectors with the highest similarity scores are returned as the most relevant results.

Cosine similarity ranges from 0 (completely unrelated) to 1 (identical meaning). In practice, scores above 0.75 indicate high relevance, and scores between 0.60 and 0.75 indicate moderate relevance.

The database that stores these vectors and performs fast similarity searches is **Milvus**, a purpose-built vector database. Milvus uses an indexing algorithm called **IVF_FLAT** (Inverted File with Flat quantization) to search millions of vectors in milliseconds rather than scanning them one by one.

### How the full RAG pipeline works in this system

Here is the complete pipeline, from question to answer:

```
User types: "What is the best AI model for detecting intracranial hemorrhage?"
  |
  v
[1] EMBED THE QUESTION
    BGE-small-en-v1.5 converts the question to a 384-dim vector.
    The model prepends a special instruction prefix:
    "Represent this sentence for searching relevant passages: ..."
    This asymmetric prefix improves retrieval quality.
  |
  v
[2] SEARCH ALL 11 COLLECTIONS (in parallel)
    The query vector is sent to Milvus, which searches all 11
    collections simultaneously using ThreadPoolExecutor:
      - imaging_literature
      - imaging_trials
      - imaging_findings
      - imaging_protocols
      - imaging_devices
      - imaging_anatomy
      - imaging_benchmarks
      - imaging_guidelines
      - imaging_report_templates
      - imaging_datasets
      - genomic_evidence (read-only, shared)

    Each collection returns its top results.
  |
  v
[3] QUERY EXPANSION
    The system detects "hemorrhage" in the question and expands
    the search to include related terms: "bleeding", "hematoma",
    "ICH", "intracranial hemorrhage", "haemorrhage", etc.
    Additional searches are run with these expanded terms.
  |
  v
[4] KNOWLEDGE GRAPH AUGMENTATION
    The system detects "intracranial_hemorrhage" and adds
    structured knowledge:
      - ICD-10: I62.9
      - Subtypes: epidural, subdural, subarachnoid,
        intraparenchymal, intraventricular
      - CT characteristics: Hyperdense (acute 50-70 HU)
      - Severity criteria: critical (>30 mL or midline shift >5 mm)
      - AI models: 3D U-Net (MONAI), VISTA-3D
      - NIM workflow: ct_head_hemorrhage
  |
  v
[5] MERGE, DEDUPLICATE, AND RANK
    All results are merged, duplicates removed, and ranked by
    weighted score. Citation relevance is assigned:
      - Score >= 0.75: high relevance
      - Score >= 0.60: medium relevance
      - Score < 0.60: low relevance
  |
  v
[6] BUILD THE PROMPT
    The top-ranked evidence is formatted into a structured prompt:
      - Section per collection (Literature, Trial, Device, etc.)
      - Each evidence item includes a citation reference
      - Knowledge graph context is appended
      - The user's original question is stated
  |
  v
[7] LLM SYNTHESIS
    The LLM (Llama 3 8B via NIM, or Claude via API) receives the
    prompt with the system instruction and generates a comprehensive,
    citation-rich answer. Response is streamed token by token.
  |
  v
[8] DISPLAY
    The Streamlit UI shows:
      - The streaming LLM response with citations
      - An expandable evidence panel with collection badges
      - Citation relevance indicators (high/medium/low)
      - Links to related workflows (e.g., ct_head_hemorrhage)
```

---

## Chapter 4: System Overview

### Architecture at a high level

The Imaging Intelligence Agent has five main layers:

```
┌───────────────────────────────────────────────────────────────┐
│                       USER INTERFACE                          │
│              Streamlit (port 8525) + FastAPI (8524)            │
│  9 tabs: Evidence | Workflows | Patient 360 | Gallery |       │
│          Protocol | Devices | Dose | Reports | Benchmarks     │
└──────────────────────────┬────────────────────────────────────┘
                           |
                           v
┌───────────────────────────────────────────────────────────────┐
│                        RAG ENGINE                             │
│  ImagingRAGEngine (src/rag_engine.py, 690 lines)              │
│  - Query expansion (12 domain maps)                           │
│  - Parallel multi-collection search                           │
│  - Comparative analysis mode (CT vs MRI, etc.)                │
│  - Citation relevance scoring                                 │
│  - Cross-modal imaging -> genomic queries                     │
└──────────────────────────┬────────────────────────────────────┘
                           |
               ┌───────────┴───────────┐
               v                       v
┌─────────────────────────┐  ┌──────────────────────────┐
│   MILVUS VECTOR DB      │  │   KNOWLEDGE GRAPH        │
│   11 Collections        │  │   3 Domains:             │
│   IVF_FLAT / COSINE     │  │   - 25 Pathologies       │
│   384 dimensions        │  │   - 9 Modalities         │
│   (BGE-small-en-v1.5)   │  │   - 21 Anatomy structs   │
└─────────────────────────┘  └──────────────────────────┘
               |
               v
┌───────────────────────────────────────────────────────────────┐
│                     NIM SERVICES                              │
│  4 services:                                                  │
│  - LLM (Llama 3 8B) -- text generation                       │
│  - VISTA-3D -- 3D medical image segmentation                  │
│  - MAISI -- medical AI synthetic imaging                      │
│  - VILA-M3 -- multimodal medical vision-language              │
└───────────────────────────────────────────────────────────────┘
```

### How the pieces connect

1. **The user** interacts with the Streamlit UI (port 8525) or sends requests to the FastAPI REST API (port 8524).
2. **The RAG engine** embeds the question using BGE-small-en-v1.5, searches all 11 Milvus collections in parallel, expands the query using 12 domain-specific expansion maps, and retrieves knowledge graph context.
3. **Milvus** (port 19530) performs fast cosine-similarity search using IVF_FLAT indexes, returning the most relevant evidence from each collection.
4. **The knowledge graph** adds structured facts (pathology characteristics, modality physics, anatomical relationships, scoring system criteria) that complement the vector search results.
5. **NIM services** provide inference capabilities: Llama 3 8B generates text responses, VISTA-3D segments 3D medical images, MAISI generates synthetic imaging data, and VILA-M3 handles multimodal vision-language tasks.

### The 11 collections as "specialized libraries"

Think of each collection as a specialized library shelf:

- **Literature** is the research library -- published papers on imaging AI methods, validation studies, and clinical outcomes
- **Trials** is the clinical registry -- active and completed clinical trials evaluating imaging AI devices and techniques
- **Findings** is the radiologist's notebook -- structured descriptions of imaging findings with severity, differential diagnosis, and follow-up recommendations
- **Protocols** is the protocol binder -- standardized acquisition parameters for each modality, body region, and clinical indication
- **Devices** is the technology catalog -- FDA-cleared AI/ML devices with performance data and regulatory status
- **Anatomy** is the anatomy atlas -- standardized anatomical structures with SNOMED codes, FMA identifiers, and segmentation labels
- **Benchmarks** is the leaderboard -- AI model performance on standardized datasets with metrics, hardware, and inference times
- **Guidelines** is the clinical standards shelf -- ACR, Fleischner, NCCN, and other society guidelines with evidence levels and key recommendations
- **Report Templates** is the structured reporting library -- standardized templates for radiology reports organized by modality and finding type
- **Datasets** is the data catalog -- public and semi-public imaging datasets with download links, licenses, and annotation details
- **Genomic Evidence** is the genome center -- shared variant-level data enabling cross-modal imaging-to-genomic queries (read-only)

When you ask a question, the system does not just search one library -- it searches all eleven simultaneously and then cross-references the findings. A question about lung nodule detection might pull evidence from literature (published AI studies), devices (FDA-cleared nodule detection products), benchmarks (model performance on LIDC-IDRI), guidelines (Fleischner Society recommendations), protocols (low-dose CT acquisition parameters), and datasets (available training data).

---

## Chapter 5: Your First Query

This chapter walks you through the experience of using the system for the first time.

### Opening the UI

Once the system is running (see Chapter 9 for setup), open your browser and navigate to:

```
http://localhost:8525
```

You will see the Imaging Intelligence Agent interface with a dark theme (NVIDIA black and green accents) and nine tabs across the top of the page.

### The nine tabs

| Tab | Purpose |
|-----|---------|
| **Evidence Explorer** | The primary search interface. Ask questions, view evidence cards, get LLM-synthesized answers with citations. |
| **Workflow Runner** | Run predefined clinical workflows (e.g., CT head hemorrhage, chest lung nodule) with step-by-step results. |
| **Patient 360** | Load a demo patient case and see all available imaging data, findings, and cross-modal correlations. |
| **Image Gallery** | Browse and preview medical images with annotation overlays and NIM-generated segmentations. |
| **Protocol Advisor** | Get protocol recommendations for specific clinical indications, including acquisition parameters. |
| **Device Ecosystem** | Explore FDA-cleared AI devices filtered by modality, body region, manufacturer, and regulatory status. |
| **Dose Intelligence** | Compare radiation dose protocols, evaluate dose reduction strategies, and review ALARA compliance. |
| **Reports & Export** | Generate structured radiology reports, export results in multiple formats (Markdown, JSON, PDF). |
| **Benchmarks** | Compare AI model performance across datasets, tasks, and hardware configurations. |

### Asking a question

Click the **Evidence Explorer** tab and type your question in the search input. For your first query, try:

```
What is the best AI approach for detecting intracranial hemorrhage on non-contrast CT?
```

Press Enter. Here is what happens:

1. **Search status** appears, showing "Searching across imaging knowledge sources..."
2. The system reports how many results it found and from which collections.
3. An **Evidence Sources** panel appears, showing the raw evidence cards organized by collection.
4. The **LLM response** streams in token by token, with markdown formatting and citations.

### Understanding the response

The response will contain several types of content:

**Cited evidence:** References to published papers, clinical trials, FDA-cleared devices, and benchmark results. Each citation links back to the evidence card in the evidence panel.

**Cross-functional insights:** The response connects evidence from different domains. For example, it might explain how a benchmark result (from imaging_benchmarks) relates to an FDA-cleared device (from imaging_devices) and a clinical guideline (from imaging_guidelines).

**Structured analysis:** For complex topics, the response often includes categorized sections (e.g., "Detection Methods," "Clinical Validation," "Regulatory Status," "Implementation Considerations").

### Understanding the evidence panel

The evidence panel shows the raw evidence that the LLM used to generate its response. Each evidence card displays:

- **Collection badge** (color-coded): Literature (blue), Trial (purple), Finding (yellow), Protocol (orange), Device (green), Anatomy (teal), Benchmark (red), Guideline (indigo), Template (pink), Dataset (brown), Genomic (cyan)
- **Record ID**: The document identifier
- **Similarity score**: How closely the evidence matches your question (0.0 to 1.0)
- **Relevance tag**: `[high]`, `[medium]`, or `[low]` based on the score
- **Text snippet**: The first portion of the evidence text

### Reading citation scores

| Score Range | Relevance | What it means |
|-------------|-----------|--------------|
| 0.75 - 1.00 | **High** | Strong semantic match to your question. This evidence is directly relevant. |
| 0.60 - 0.74 | **Medium** | Partial match. The evidence is related but may address a subtopic or adjacent concept. |
| 0.40 - 0.59 | **Low** | Weak match. The evidence has some thematic overlap but may not directly answer your question. |
| Below 0.40 | Filtered out | Not returned. Below the minimum score threshold. |

### Using sidebar controls

The sidebar gives you fine-grained control over searches:

- **Modality Filter**: Restrict results to a specific modality (e.g., CT, MRI, ultrasound). This adds a Milvus field filter on collections that have a `modality` field.
- **Body Region Filter**: Focus on a specific anatomical region (e.g., head, chest, abdomen).
- **Date Range**: Filter evidence by publication year or trial start year.
- **Collection toggles**: Enable or disable specific collections. Each toggle shows the live record count.
- **Demo Cases**: Load a prebuilt patient scenario for guided exploration.

### The six predefined workflows

From the **Workflow Runner** tab, you can launch predefined clinical workflows:

| Workflow | Description |
|----------|-------------|
| `ct_head_hemorrhage` | Detect and classify intracranial hemorrhage subtypes on non-contrast CT |
| `ct_chest_lung_nodule` | Detect, measure, and classify lung nodules with Lung-RADS scoring |
| `ct_coronary_angiography` | Evaluate coronary artery disease with CAD-RADS scoring |
| `cxr_rapid_findings` | Rapid triage of chest X-ray findings (pneumonia, pneumothorax, effusion) |
| `mri_brain_ms_lesion` | Detect and quantify multiple sclerosis lesions on brain MRI |
| `mri_prostate_pirads` | Evaluate prostate lesions with PI-RADS scoring on multiparametric MRI |

Each workflow orchestrates a multi-step pipeline: retrieve relevant evidence, apply the appropriate scoring system, consult guidelines, and generate a structured report.

### The four demo cases

The system includes four demo patient cases for guided exploration:

| Case ID | Description | Key Finding |
|---------|-------------|-------------|
| **DEMO-001** | 62-year-old male | Intracranial hemorrhage (ICH) |
| **DEMO-002** | 58-year-old female | Lung nodule |
| **DEMO-003** | 55-year-old male | Coronary artery disease (CAD) |
| **DEMO-004** | 45-year-old female | Pneumonia |

Load any demo case from the Patient 360 tab or the sidebar to see how the system integrates findings, protocols, guidelines, and device recommendations for a specific clinical scenario.

### Downloading results

After each query or workflow run, you can export results:

- **Markdown**: A formatted `.md` file with the query, response, evidence, and metadata
- **JSON**: A structured `.json` file suitable for programmatic processing
- **PDF**: A formatted PDF report suitable for clinical documentation

---

## Chapter 6: Understanding Collections

### The 11 collections at a glance

| # | Collection | Purpose | Source |
|---|-----------|---------|--------|
| 1 | `imaging_literature` | Published research on imaging AI | PubMed, journals |
| 2 | `imaging_trials` | Clinical trials evaluating imaging AI | ClinicalTrials.gov |
| 3 | `imaging_findings` | Structured radiology findings | Clinical knowledge bases |
| 4 | `imaging_protocols` | Imaging acquisition protocols | Institutional protocols, ACR |
| 5 | `imaging_devices` | FDA-cleared AI/ML devices | FDA database, manufacturer data |
| 6 | `imaging_anatomy` | Anatomical structures and labels | SNOMED CT, FMA, VISTA-3D |
| 7 | `imaging_benchmarks` | AI model performance data | RSNA, MICCAI, published studies |
| 8 | `imaging_guidelines` | Clinical imaging guidelines | ACR, Fleischner, NCCN, AHA |
| 9 | `imaging_report_templates` | Structured reporting templates | ACR, RadLex, institutional |
| 10 | `imaging_datasets` | Public imaging datasets | NIH, TCIA, Kaggle, grand challenges |
| 11 | `genomic_evidence` | Patient variant data (shared) | VCF data from rag-chat-pipeline |

**10 imaging-specific collections + 1 shared genomic collection (read-only) = 11 total**

### Collection details

#### 1. imaging_literature

**What it contains:** Published research papers on medical imaging AI -- methods papers, validation studies, clinical outcome analyses, systematic reviews, and meta-analyses. Each record captures the title, text content, source type, publication year, imaging modality, body region, AI task, disease focus, keywords, and journal.

**Why it matters:** Published literature is the foundation of evidence-based radiology. When you ask "What deep learning architectures perform best for lung nodule detection?", the literature collection contains the peer-reviewed evidence.

**Example questions it helps answer:**
- "What is the sensitivity of AI-assisted mammography screening?"
- "How does federated learning improve multi-site imaging AI models?"
- "What are the latest advances in diffusion MRI analysis?"

**Key fields:** `title`, `text_chunk`, `source_type`, `year`, `modality`, `body_region`, `ai_task`, `disease`, `keywords`, `journal`

#### 2. imaging_trials

**What it contains:** Clinical trial records for imaging AI studies, including the trial title, summary, phase, recruitment status, sponsor, imaging modality, body region, AI task, disease focus, enrollment size, start year, and outcome summary.

**Why it matters:** Clinical trials represent the bridge between algorithm development and clinical deployment. This collection lets you find active trials for specific imaging AI applications, compare study designs, and track emerging evidence.

**Example questions it helps answer:**
- "Are there any Phase 3 trials for AI-assisted breast cancer screening?"
- "What imaging AI trials is Siemens Healthineers sponsoring?"
- "How large are the enrollment cohorts for CT lung screening AI trials?"

**Key fields:** `title`, `text_summary`, `phase`, `status`, `sponsor`, `modality`, `body_region`, `ai_task`, `disease`, `enrollment`, `start_year`, `outcome_summary`

#### 3. imaging_findings

**What it contains:** Structured descriptions of radiology findings -- the observations that radiologists document when reading imaging studies. Each record includes the finding text, category, severity level, associated modality and body region, clinical significance, differential diagnosis, recommended follow-up, measurement type and value, classification system, and classification score.

**Why it matters:** Findings are the core output of radiology. This collection encodes expert knowledge about what findings look like, how severe they are, what they might represent (differential diagnosis), and what should be done next. It supports 37 finding categories across 16 body regions at 5 severity levels.

**Example questions it helps answer:**
- "What is the differential diagnosis for a ground-glass opacity on chest CT?"
- "What follow-up is recommended for a PI-RADS 4 prostate lesion?"
- "How is the severity of intracranial hemorrhage classified?"

**Key fields:** `text_summary`, `finding_category`, `severity`, `modality`, `body_region`, `clinical_significance`, `differential_diagnosis`, `recommended_followup`, `measurement_type`, `measurement_value`, `classification_system`, `classification_score`

#### 4. imaging_protocols

**What it contains:** Imaging acquisition protocols -- the technical parameters used to acquire imaging studies. Each record documents the protocol name, modality, body region, contrast agent (if any), slice thickness, radiation dose, scan duration, clinical indication, and preprocessing steps.

**Why it matters:** Protocol optimization is critical for image quality, radiation safety, and AI model performance. The wrong protocol can render an AI algorithm ineffective. This collection helps researchers understand which protocols are standard for each clinical indication and how acquisition parameters affect downstream AI analysis.

**Example questions it helps answer:**
- "What is the recommended slice thickness for low-dose CT lung screening?"
- "Which contrast agents are used for cardiac CT angiography?"
- "How can radiation dose be reduced for pediatric abdominal CT?"

**Key fields:** `text_summary`, `protocol_name`, `modality`, `body_region`, `contrast_agent`, `slice_thickness_mm`, `radiation_dose`, `scan_duration`, `clinical_indication`, `preprocessing_steps`

#### 5. imaging_devices

**What it contains:** FDA-cleared AI/ML medical devices for radiology. Each record documents the device name, manufacturer, regulatory status, clearance date, imaging modality, body region, AI task, intended use, performance summary, and model architecture.

**Why it matters:** Understanding which AI devices are commercially available, how they are regulated, and how they perform is essential for clinical adoption decisions. This collection bridges the gap between research benchmarks and real-world device performance.

**Example questions it helps answer:**
- "Which FDA-cleared devices detect pulmonary embolism on CT?"
- "What is the performance of Aidoc's intracranial hemorrhage detection?"
- "How many AI devices have received De Novo classification?"

**Key fields:** `text_summary`, `device_name`, `manufacturer`, `regulatory_status`, `clearance_date`, `modality`, `body_region`, `ai_task`, `intended_use`, `performance_summary`, `model_architecture`

#### 6. imaging_anatomy

**What it contains:** Anatomical structure definitions relevant to medical imaging. Each record documents the structure name, body region, organ system, SNOMED CT code, FMA identifier, imaging characteristics (how the structure appears on different modalities), common pathologies affecting that structure, and the VISTA-3D segmentation label ID.

**Why it matters:** Anatomy is the foundation of radiology. Understanding where structures are, how they appear on imaging, and what pathologies affect them is essential for interpreting AI segmentation results and correlating findings across modalities. The system covers 21 anatomy structures.

**Example questions it helps answer:**
- "What is the SNOMED code for the left ventricle?"
- "How does the hippocampus appear on T2-weighted MRI?"
- "What pathologies commonly affect the pancreas?"

**Key fields:** `text_summary`, `structure_name`, `body_region`, `system`, `snomed_code`, `fma_id`, `imaging_characteristics`, `common_pathologies`, `segmentation_label_id`

#### 7. imaging_benchmarks

**What it contains:** AI model performance data from benchmark datasets and competitions. Each record documents the model name, architecture, AI task, modality, body region, dataset, metric name and value, training data size, inference time in milliseconds, and hardware used.

**Why it matters:** Benchmarks provide standardized comparisons between AI approaches. This collection lets you find the state-of-the-art model for a specific task and understand the tradeoffs between accuracy, speed, and hardware requirements.

**Example questions it helps answer:**
- "What is the best-performing model on the RSNA pneumonia detection challenge?"
- "How does nnU-Net compare to VISTA-3D for liver segmentation?"
- "What inference time can I expect for a 3D segmentation model on an A100 GPU?"

**Key fields:** `text_summary`, `model_name`, `model_architecture`, `ai_task`, `modality`, `body_region`, `dataset_name`, `metric_name`, `metric_value`, `training_data_size`, `inference_time_ms`, `hardware`

#### 8. imaging_guidelines

**What it contains:** Clinical imaging guidelines from professional societies. Each record documents the guideline name, issuing organization, publication year, modality, body region, clinical indication, classification system, key recommendation, and evidence level.

**Why it matters:** Guidelines represent expert consensus on best practices. They define when to image, how to image, how to interpret findings, and what to do next. AI systems that align with clinical guidelines are more likely to be adopted in clinical practice.

**Example questions it helps answer:**
- "What does the Fleischner Society recommend for incidental pulmonary nodules?"
- "What is the ACR Appropriateness Criteria for low back pain imaging?"
- "What evidence level supports screening mammography for women aged 40-49?"

**Key fields:** `text_summary`, `guideline_name`, `organization`, `year`, `modality`, `body_region`, `clinical_indication`, `classification_system`, `key_recommendation`, `evidence_level`

#### 9. imaging_report_templates

**What it contains:** Structured reporting templates for radiology reports. Each record documents the template name, modality, body region, finding type, structured fields (the data points to include), an example report, and the coding system used (CPT, ICD-10, RadLex).

**Why it matters:** Structured reporting reduces variability, improves communication with referring clinicians, and enables data mining of radiology reports. AI-generated reports should follow these templates to integrate smoothly into clinical workflows.

**Example questions it helps answer:**
- "What structured fields should a CT lung screening report include?"
- "What does a standard BI-RADS mammography report template look like?"
- "What ICD-10 codes are associated with intracranial hemorrhage findings?"

**Key fields:** `text_summary`, `template_name`, `modality`, `body_region`, `finding_type`, `structured_fields`, `example_report`, `coding_system`

#### 10. imaging_datasets

**What it contains:** Public and semi-public imaging datasets available for AI research and model training. Each record documents the dataset name, source institution, modality, body region, number of studies, number of images, disease labels, annotation type, license, and download URL.

**Why it matters:** Access to high-quality, well-annotated imaging data is the primary bottleneck in medical imaging AI research. This collection helps researchers find the right dataset for their task, understand licensing restrictions, and plan training data strategies.

**Example questions it helps answer:**
- "What chest X-ray datasets have bounding box annotations?"
- "Where can I download a large CT dataset with liver segmentation labels?"
- "Which imaging datasets have permissive (MIT/Apache) licenses?"

**Key fields:** `text_summary`, `dataset_name`, `source`, `modality`, `body_region`, `num_studies`, `num_images`, `disease_labels`, `annotation_type`, `license_type`, `download_url`

#### 11. genomic_evidence (read-only)

**What it contains:** Patient-level variant data from whole-genome sequencing, shared from the HCLS AI Factory's rag-chat-pipeline. This collection is read-only -- the Imaging Intelligence Agent can query it but does not write to it.

**Why it matters:** Radiogenomics -- the correlation between imaging features and genomic markers -- is an emerging field. A tumor that appears spiculated on CT may have a different mutational profile than one with smooth margins. Cross-modal queries that bridge imaging and genomics are a unique capability of this system.

**Example questions it helps answer:**
- "Are there genomic markers associated with ground-glass opacity morphology in lung adenocarcinoma?"
- "What imaging features correlate with BRCA1 mutation status in breast cancer?"

**Key fields:** Inherited from the rag-chat-pipeline (variant ID, gene, consequence, clinical significance, etc.)

### Cross-reference guide

| If you are interested in... | Primary collection | Cross-reference with |
|----------------------------|-------------------|---------------------|
| AI model performance | Benchmarks | Literature, Devices, Datasets |
| Clinical best practices | Guidelines | Literature, Protocols, Findings |
| FDA device selection | Devices | Benchmarks, Literature, Trials |
| Protocol optimization | Protocols | Guidelines, Devices, Anatomy |
| Structured reporting | Report Templates | Findings, Guidelines |
| Dataset discovery | Datasets | Benchmarks, Literature |
| Anatomical context | Anatomy | Findings, Protocols |
| Clinical trial landscape | Trials | Literature, Devices |
| Imaging-genomic correlations | Genomic Evidence | Literature, Findings |

---

## Chapter 7: The Knowledge Graph

### What it is

The knowledge graph is a hand-curated database of structured facts about medical imaging. Unlike the vector database (which stores free-text evidence as embeddings), the knowledge graph stores explicit, structured relationships.

Think of the difference this way:

- **Vector database:** "There is a document that mentions intracranial hemorrhage, CT, and non-contrast protocol. It is probably relevant to your question."
- **Knowledge graph:** "Intracranial hemorrhage (ICD-10: I62.9) has 5 subtypes: epidural, subdural, subarachnoid, intraparenchymal, intraventricular. On CT, acute hemorrhage appears hyperdense (50-70 HU). Critical severity: volume > 30 mL or midline shift > 5 mm. The NIM workflow ct_head_hemorrhage can be triggered."

The knowledge graph provides precise, factual context that helps the LLM generate more accurate and specific answers.

### The three knowledge domains

#### 1. Pathologies (25 entries)

For each of 25 pathologies, the knowledge graph stores:

- ICD-10 code and display name
- Applicable imaging modalities
- Body region
- Subtypes (e.g., hemorrhage subtypes: epidural, subdural, subarachnoid, intraparenchymal, intraventricular)
- CT characteristics (Hounsfield unit ranges, appearance descriptors)
- MRI characteristics (T1/T2 signal, DWI behavior, enhancement patterns)
- Severity criteria (critical, urgent, routine -- with specific measurement thresholds)
- Key measurements (e.g., volume_ml, midline_shift_mm, max_thickness_mm)
- Clinical guidelines (e.g., Brain Trauma Foundation, AHA/ASA Stroke Guidelines)
- AI models recommended for the pathology
- NIM workflow to invoke

**Examples of pathology entries:**

| Pathology | ICD-10 | Modalities | Body Region | NIM Workflow |
|-----------|--------|------------|-------------|--------------|
| Intracranial Hemorrhage | I62.9 | CT, MRI | Head | ct_head_hemorrhage |
| Pulmonary Nodule | R91.1 | CT | Chest | ct_chest_lung_nodule |
| Coronary Artery Disease | I25.10 | CT | Cardiac | ct_coronary_angiography |
| Pneumonia | J18.9 | CXR, CT | Chest | cxr_rapid_findings |
| MS Lesions | G35 | MRI | Head | mri_brain_ms_lesion |
| Prostate Lesion | N42.9 | MRI | Pelvis | mri_prostate_pirads |

The pathology knowledge graph encodes clinical decision-making logic. For lung nodules, the severity criteria map directly to Lung-RADS categories:

```
lung_nodule:
  severity_criteria:
    critical: "Lung-RADS 4B or 4X, size > 15 mm"
    urgent:   "Lung-RADS 4A, size 8-15 mm or growing"
    routine:  "Lung-RADS 1-3, size < 6 mm"
  key_measurements:
    - long_axis_mm
    - short_axis_mm
    - volume_mm3
    - doubling_time_days
```

#### 2. Modalities (9 entries)

For each of 9 imaging modalities, the knowledge graph stores:

- Full name and physics description
- Spatial resolution range
- Typical scan duration
- Radiation dose (where applicable)
- Strengths and limitations
- Common AI applications
- Key acquisition parameters
- Contrast agents used
- Body regions covered

This domain helps the system understand the capabilities and constraints of each modality. When a user asks "Should I use CT or MRI for evaluating a liver lesion?", the knowledge graph provides the factual basis for a modality comparison.

#### 3. Anatomy (21 entries)

For each of 21 anatomical structures, the knowledge graph stores:

- Structure name and body region
- Organ system
- SNOMED CT code and FMA identifier
- Imaging characteristics on CT, MRI, and ultrasound
- Common pathologies affecting the structure
- VISTA-3D segmentation label ID
- Adjacent structures and anatomical relationships

**Examples of anatomy entries:**

| Structure | Body Region | System | SNOMED Code | Common Pathologies |
|-----------|-------------|--------|-------------|-------------------|
| Cerebral Cortex | Head | Nervous | 40733004 | Stroke, tumor, atrophy |
| Lung (Right) | Chest | Respiratory | 3341006 | Nodule, pneumonia, PE |
| Liver | Abdomen | Digestive | 10200004 | HCC, metastasis, cirrhosis |
| Left Ventricle | Cardiac | Cardiovascular | 87878005 | MI, cardiomyopathy |
| Prostate | Pelvis | Reproductive | 41216001 | Cancer, BPH |

### How the knowledge graph differs from the vector database

| Feature | Vector Database (Milvus) | Knowledge Graph (Python dicts) |
|---------|-------------------------|-------------------------------|
| Data type | Free-text embeddings | Structured key-value facts |
| Size | Varies by collection | ~55 entities (25 + 9 + 21) |
| Search method | Cosine similarity | Keyword matching |
| Strength | Finding relevant text you did not know existed | Providing precise facts about known entities |
| Weakness | May return tangentially related text | Only covers pre-curated entities |
| Role in RAG | Retrieval (finding evidence) | Augmentation (adding structured context) |

Both are essential. The vector database finds evidence. The knowledge graph ensures the LLM has the right factual context.

---

## Chapter 8: Query Expansion

### The problem

When you search for "hemorrhage," you want results that also mention:

- "bleeding" (a common synonym)
- "hematoma" (a related condition)
- "ICH" (an abbreviation for intracranial hemorrhage)
- "haemorrhage" (British spelling)
- "blood" (the substance involved)

A pure vector similarity search will catch some of these (because the embeddings capture semantic relationships), but it may miss others -- especially when the related term uses completely different words or abbreviations.

Query expansion solves this by explicitly mapping keywords to their related terms.

### The 12 expansion maps

The system contains 12 hand-curated expansion dictionaries, covering imaging-specific domains:

| # | Map | Example Keywords | Domain |
|---|-----|-----------------|--------|
| 1 | Modality | ct, mri, xray, cxr, pet, ultrasound, mammography, fluoroscopy | Imaging modality synonyms |
| 2 | Body Region | head, neck, chest, abdomen, pelvis, spine, extremity, cardiac, breast, shoulder, musculoskeletal | Anatomical region synonyms |
| 3 | Pathology | hemorrhage, nodule, mass, fracture, pneumonia, effusion, pneumothorax, stroke, embolism, dissection, stenosis, fibrosis, emphysema, lymphadenopathy, breast mass, thyroid nodule, pancreatic mass, bowel obstruction | Disease and finding synonyms |
| 4 | AI Task | detection, segmentation, classification, registration, reconstruction | ML task synonyms |
| 5 | Scoring System | lung-rads, bi-rads, ti-rads, li-rads, cad-rads, pi-rads, aspects | Standardized scoring terms |
| 6 | Anatomy | brain, liver, kidney, lung, heart, spine, prostate | Structure-specific terms |
| 7 | Device | fda, clearance, 510k, de novo, pma | Regulatory terms |
| 8 | Protocol | contrast, dose, slice, sequence, reconstruction | Acquisition parameter terms |
| 9 | Guideline | acr, fleischner, nccn, aha, eso | Society guideline terms |
| 10 | Metric | sensitivity, specificity, auc, dice, accuracy | Performance metric terms |
| 11 | Dataset | rsna, miccai, nih, tcia, lidc, chestxray | Dataset name terms |
| 12 | Architecture | unet, resnet, transformer, cnn, vision transformer | ML architecture terms |

### How expansion works

When you ask a question, the `expand_query()` function:

1. Converts the question to lowercase
2. Scans it for any keyword that appears in any of the 12 maps
3. Collects all related terms for every matched keyword
4. Deduplicates and returns the expansion terms

For example, given the query "What CT protocol is best for detecting lung nodules?":

- **"ct"** matches the Modality map, expanding to: "computed tomography", "ct scan", "cat scan", "helical ct", "spiral ct", "mdct", "multi-detector ct"
- **"lung"** matches the Body Region map (via "chest"), expanding to: "thorax", "thoracic", "pulmonary", "mediastinal", "pleural"
- **"nodule"** matches the Pathology map, expanding to: "nodular", "pulmonary nodule", "lung nodule", "solitary pulmonary nodule", "spn"

### How expanded terms are used

Expanded terms are used in two ways:

1. **Modality and body region terms** can be used as **field filters** on collections that have `modality` or `body_region` fields. This is a precise filter: `modality == "ct"`.

2. **Synonym and alias terms** are **re-embedded** as separate queries and searched across all collections. This semantic search finds additional evidence the original query might have missed.

The results from expanded searches are merged with the original results, deduplicated, and scored at a lower weight to prioritize the original direct-match results.

### Modality and region context injection

Beyond synonyms, query expansion also injects modality and body region context from the knowledge graph. If the system detects that the query is about "CT of the chest," it can add structured context about CT physics (Hounsfield unit ranges, typical slice thickness, radiation dose ranges) and chest anatomy (lung lobes, mediastinal structures, pleural spaces) to the prompt. This helps the LLM generate more technically precise responses.

---

## Chapter 9: Setting Up Locally

### Prerequisites

You will need:

- **Python 3.10+** (check with `python3 --version`)
- **Milvus 2.4** running on `localhost:19530` (see below for Docker setup)
- **Approximately 4 GB of disk space** for the vector database and models
- **Internet access** for the initial data ingestion from PubMed and ClinicalTrials.gov

Optional but recommended:

- **NVIDIA GPU** with CUDA 12.x for NIM services (VISTA-3D, MAISI, VILA-M3)
- **Docker** for containerized deployment
- **An LLM API key** (for Anthropic Claude or a local Llama 3 8B NIM instance)

### Step 1: Clone the repository

```bash
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd hcls-ai-factory/ai_agent_adds/imaging_intelligence_agent
```

### Step 2: Install Python dependencies

```bash
cd agent
pip install -r requirements.txt
```

This installs the core stack:
- `pymilvus` (vector database client)
- `sentence-transformers` (BGE-small-en-v1.5 embedding model)
- `anthropic` (Claude API client, optional)
- `streamlit` (UI framework)
- `fastapi` + `uvicorn` (REST API)
- `pydantic` + `pydantic-settings` (data models and configuration)
- `biopython` + `lxml` (PubMed data parsing)
- `reportlab` (PDF export)
- `loguru` (structured logging)
- `apscheduler` (automated data refresh)
- `pydicom` (DICOM file handling)

### Step 3: Start Milvus

If you have Docker installed, the simplest approach is:

```bash
# Start Milvus standalone (creates a local vector database on port 19530)
docker compose up -d milvus-standalone
```

Or if you are using the HCLS AI Factory's full docker-compose:

```bash
cd /home/adam/projects/hcls-ai-factory
docker compose up -d milvus-standalone milvus-etcd milvus-minio
```

Verify Milvus is running:

```bash
curl http://localhost:19530/v1/vector/collections
```

### Step 4: Create collections and seed initial data

```bash
python3 scripts/setup_collections.py
```

This creates all 10 imaging-specific Milvus collections with IVF_FLAT indexes (384 dimensions, COSINE metric). The shared `genomic_evidence` collection is not created here -- it is managed by the rag-chat-pipeline.

### Step 5: Seed domain-specific data

```bash
# Seed literature from PubMed
python3 scripts/seed_literature.py

# Seed clinical trials
python3 scripts/seed_trials.py

# Seed structured findings
python3 scripts/seed_findings.py

# Seed imaging protocols
python3 scripts/seed_protocols.py

# Seed FDA-cleared devices
python3 scripts/seed_devices.py

# Seed anatomical structures
python3 scripts/seed_anatomy.py

# Seed AI benchmarks
python3 scripts/seed_benchmarks.py

# Seed clinical guidelines
python3 scripts/seed_guidelines.py

# Seed report templates
python3 scripts/seed_report_templates.py

# Seed public datasets
python3 scripts/seed_datasets.py
```

### Step 6: Ingest PubMed literature (optional, for live data)

```bash
python3 scripts/ingest_pubmed.py --max-results 5000
```

This fetches imaging AI abstracts from PubMed via NCBI E-utilities, classifies each paper by modality and AI task, embeds the text with BGE-small-en-v1.5, and stores the results in `imaging_literature`.

### Step 7: Ingest clinical trials (optional, for live data)

```bash
python3 scripts/ingest_clinical_trials.py --max-results 1500
```

This fetches imaging AI trials from the ClinicalTrials.gov API v2, extracts phase, status, modality, and other metadata, embeds the trial summaries, and stores them in `imaging_trials`.

### Step 8: Validate the setup

```bash
python3 scripts/validate_e2e.py
```

This runs end-to-end validation: collection statistics, single-collection search, multi-collection search, filtered search, and all demo queries.

### Step 9: Launch the Streamlit UI

```bash
streamlit run app/imaging_ui.py --server.port 8525
```

Open `http://localhost:8525` in your browser. You should see the Imaging Intelligence Agent interface with live collection statistics in the sidebar and all nine tabs functional.

### Step 10: (Optional) Launch the REST API

```bash
uvicorn api.main:app --host 0.0.0.0 --port 8524 --reload
```

The API documentation is automatically available at `http://localhost:8524/docs` (Swagger UI).

### Step 11: (Optional) Start NIM services

If you have an NVIDIA GPU, you can start the four NIM services:

```bash
python3 scripts/test_nim_services.py
```

This verifies connectivity to the LLM (Llama 3 8B), VISTA-3D, MAISI, and VILA-M3 NIM endpoints.

### Docker Lite: minimal setup

For a minimal setup without NIM services or GPU:

```bash
cd /home/adam/projects/hcls-ai-factory
docker compose up -d milvus-standalone milvus-etcd milvus-minio
cd ai_agent_adds/imaging_intelligence_agent/agent
pip install -r requirements.txt
python3 scripts/setup_collections.py
python3 scripts/seed_findings.py
python3 scripts/seed_protocols.py
python3 scripts/seed_devices.py
python3 scripts/seed_anatomy.py
python3 scripts/seed_guidelines.py
streamlit run app/imaging_ui.py --server.port 8525
```

This gives you a working system with core collections seeded, without requiring an LLM API key or GPU. Evidence retrieval works fully; LLM synthesis will be unavailable until you configure an LLM endpoint.

### Environment variables

All configuration is managed through environment variables with the `IMAGING_` prefix. Key variables:

| Variable | Default | Description |
|----------|---------|-------------|
| `IMAGING_MILVUS_HOST` | localhost | Milvus server hostname |
| `IMAGING_MILVUS_PORT` | 19530 | Milvus server port |
| `IMAGING_EMBEDDING_MODEL` | BAAI/bge-small-en-v1.5 | Embedding model name |
| `IMAGING_LLM_MODEL` | claude-sonnet-4-6 | Default LLM (Claude Sonnet 4.6) |
| `IMAGING_STREAMLIT_PORT` | 8525 | Streamlit UI port |
| `IMAGING_API_PORT` | 8524 | FastAPI REST API port |

---

## Chapter 10: Exploring the API

### What the REST API is

The REST API wraps the same RAG engine that powers the Streamlit UI, but exposes it as a set of HTTP endpoints. This enables:

- Integration with other applications (PACS viewers, EHR systems, research notebooks)
- Programmatic access from any programming language
- Automated testing and monitoring
- Pipeline orchestration with Nextflow or Airflow

### The 19+ endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/health` | Service health with collection count |
| GET | `/collections` | All collection names and their record counts |
| POST | `/query` | Full RAG query (retrieve + LLM synthesis) |
| POST | `/search` | Evidence-only retrieval (no LLM, fast) |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/knowledge/stats` | Knowledge graph statistics (pathologies, modalities, anatomy) |
| GET | `/metrics` | Prometheus-compatible metrics |
| POST | `/api/ask` | Alternative query endpoint |
| GET | `/workflows` | List all available workflows |
| GET | `/workflows/{name}` | Get details for a specific workflow |
| POST | `/workflows/{name}/run` | Execute a clinical workflow |
| POST | `/reports/generate` | Generate a structured radiology report |
| GET | `/nim` | List NIM service status |
| GET | `/nim/{service}` | Get status of a specific NIM service |
| POST | `/nim/check-all` | Health check all NIM services |
| POST | `/preview/generate` | Generate an image preview with annotations |
| GET | `/protocol/recommend` | Get protocol recommendation for a clinical indication |
| POST | `/dose/compare` | Compare radiation dose between protocols |
| GET | `/demo-cases` | List available demo patient cases |

### Testing with curl

**Health check:**

```bash
curl http://localhost:8524/health
```

Expected response:

```json
{
  "status": "healthy",
  "collections": 11,
  "version": "1.0.0"
}
```

**List collections:**

```bash
curl http://localhost:8524/collections
```

Expected response:

```json
{
  "collections": [
    {"name": "imaging_literature", "record_count": 5047},
    {"name": "imaging_trials", "record_count": 973},
    {"name": "imaging_findings", "record_count": 240},
    {"name": "imaging_protocols", "record_count": 85},
    {"name": "imaging_devices", "record_count": 120},
    {"name": "imaging_anatomy", "record_count": 21},
    {"name": "imaging_benchmarks", "record_count": 95},
    {"name": "imaging_guidelines", "record_count": 60},
    {"name": "imaging_report_templates", "record_count": 30},
    {"name": "imaging_datasets", "record_count": 45}
  ],
  "total": 11
}
```

**Evidence-only search (fast, no LLM):**

```bash
curl -X POST http://localhost:8524/search \
  -H "Content-Type: application/json" \
  -d '{
    "question": "lung nodule detection on CT",
    "modality": "ct",
    "body_region": "chest"
  }'
```

**Full RAG query (with LLM synthesis):**

```bash
curl -X POST http://localhost:8524/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the current evidence for AI-assisted lung nodule detection?",
    "modality": "ct",
    "body_region": "chest"
  }'
```

**Run a clinical workflow:**

```bash
curl -X POST http://localhost:8524/workflows/ct_head_hemorrhage/run \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "DEMO-001"
  }'
```

**Get protocol recommendation:**

```bash
curl "http://localhost:8524/protocol/recommend?modality=ct&body_region=chest&indication=lung_screening"
```

**Compare radiation doses:**

```bash
curl -X POST http://localhost:8524/dose/compare \
  -H "Content-Type: application/json" \
  -d '{
    "protocol_a": "standard_chest_ct",
    "protocol_b": "low_dose_chest_ct"
  }'
```

**Generate a structured report:**

```bash
curl -X POST http://localhost:8524/reports/generate \
  -H "Content-Type: application/json" \
  -d '{
    "modality": "ct",
    "body_region": "chest",
    "findings": ["6mm solid nodule in RUL"],
    "template": "lung_screening"
  }'
```

**Check NIM service status:**

```bash
curl http://localhost:8524/nim
```

Expected response:

```json
{
  "services": [
    {"name": "llm", "model": "meta/llama-3-8b-instruct", "status": "available"},
    {"name": "vista3d", "model": "nvidia/vista-3d", "status": "available"},
    {"name": "maisi", "model": "nvidia/maisi", "status": "available"},
    {"name": "vilam3", "model": "nvidia/vila-m3", "status": "available"}
  ]
}
```

**List demo cases:**

```bash
curl http://localhost:8524/demo-cases
```

Expected response:

```json
{
  "cases": [
    {"id": "DEMO-001", "description": "62M, intracranial hemorrhage"},
    {"id": "DEMO-002", "description": "58F, lung nodule"},
    {"id": "DEMO-003", "description": "55M, coronary artery disease"},
    {"id": "DEMO-004", "description": "45F, pneumonia"}
  ]
}
```

**Cross-collection entity search:**

```bash
curl -X POST http://localhost:8524/find-related \
  -H "Content-Type: application/json" \
  -d '{
    "entity": "intracranial hemorrhage",
    "top_k": 5
  }'
```

This returns everything the system knows about intracranial hemorrhage across all 11 collections -- literature, trials, findings, devices, guidelines, benchmarks, anatomy, and more.

### API documentation

The FastAPI framework automatically generates interactive API documentation. Once the API is running, open:

- **Swagger UI:** `http://localhost:8524/docs`
- **ReDoc:** `http://localhost:8524/redoc`

These pages let you explore endpoints, view request/response schemas, and test API calls directly from the browser.

---

## Chapter 11: Understanding the Codebase

### Project structure

```
imaging_intelligence_agent/agent/
├── src/                              # Core source code (the engine)
│   ├── __init__.py
│   ├── models.py                     # Pydantic data models (728 lines)
│   ├── collections.py                # 11 Milvus collection schemas + manager (1,420 lines)
│   ├── knowledge.py                  # Knowledge graph: 3 domains (1,843 lines)
│   ├── query_expansion.py            # 12 expansion maps
│   ├── rag_engine.py                 # Multi-collection RAG engine (690 lines)
│   ├── agent.py                      # Autonomous agent: plan-search-synthesize (283 lines)
│   ├── cross_modal.py                # Imaging <-> genomic cross-modal queries (401 lines)
│   ├── export.py                     # Markdown, JSON, PDF export (691 lines)
│   ├── demo_cases.py                 # 4 demo patient cases
│   ├── protocol_optimizer.py         # Protocol recommendation engine
│   ├── scheduler.py                  # Automated data refresh
│   ├── imaging/
│   │   └── preview_generator.py      # Image preview with annotations
│   ├── nim/                          # NVIDIA NIM service clients
│   │   ├── base.py                   # Base NIM client
│   │   ├── llm_client.py             # Llama 3 8B client
│   │   ├── vista3d_client.py         # VISTA-3D segmentation client
│   │   ├── maisi_client.py           # MAISI synthetic imaging client
│   │   ├── vilam3_client.py          # VILA-M3 vision-language client
│   │   └── service_manager.py        # NIM service orchestration
│   ├── ingest/                       # Data ingestion parsers
│   │   ├── base.py                   # Base ingest pipeline
│   │   ├── literature_parser.py      # PubMed NCBI E-utilities
│   │   ├── clinical_trials_parser.py # ClinicalTrials.gov API v2
│   │   ├── finding_parser.py         # Radiology findings
│   │   ├── protocol_parser.py        # Imaging protocols
│   │   ├── device_parser.py          # FDA device data
│   │   ├── anatomy_parser.py         # Anatomical structures
│   │   ├── benchmark_parser.py       # AI benchmark results
│   │   ├── guideline_parser.py       # Clinical guidelines
│   │   ├── report_template_parser.py # Report templates
│   │   ├── dataset_parser.py         # Public datasets
│   │   └── dicom_watcher.py          # DICOM file monitoring
│   ├── workflows/                    # 6 clinical workflows
│   │   ├── base.py                   # Base workflow class
│   │   ├── ct_head_hemorrhage.py     # ICH detection and classification
│   │   ├── ct_chest_lung_nodule.py   # Lung nodule with Lung-RADS
│   │   ├── ct_coronary_angiography.py # CAD-RADS scoring
│   │   ├── cxr_rapid_findings.py     # CXR triage
│   │   ├── mri_brain_ms_lesion.py    # MS lesion detection
│   │   └── mri_prostate_pirads.py    # PI-RADS scoring
│   └── utils/
│       └── pubmed_client.py          # NCBI E-utilities HTTP client
├── app/
│   ├── imaging_ui.py                 # Streamlit entry point (9 tabs)
│   └── tabs/                         # Individual UI tab implementations
│       ├── evidence_explorer.py      # Evidence search tab
│       ├── workflow_runner.py        # Clinical workflow tab
│       ├── patient_360.py            # Patient overview tab
│       ├── image_gallery.py          # Image browser tab
│       ├── protocol_advisor.py       # Protocol recommendation tab
│       ├── device_ecosystem.py       # Device explorer tab
│       ├── dose_intelligence.py      # Dose comparison tab
│       ├── reports_export.py         # Report generation tab
│       └── benchmarks.py            # Benchmark comparison tab
├── api/
│   ├── __init__.py
│   ├── main.py                       # FastAPI REST API (547 lines, 19+ endpoints)
│   └── routes/                       # API route modules
│       ├── workflows.py              # Workflow endpoints
│       ├── nim.py                    # NIM service endpoints
│       ├── preview.py               # Image preview endpoints
│       ├── protocol.py              # Protocol recommendation endpoints
│       ├── dose.py                  # Dose comparison endpoints
│       ├── reports.py               # Report generation endpoints
│       ├── demo_cases.py            # Demo case endpoints
│       ├── events.py                # Server-sent event endpoints
│       └── meta_agent.py           # Meta-agent orchestration endpoints
├── config/
│   └── settings.py                   # Pydantic BaseSettings configuration (158 lines)
├── scripts/                          # CLI scripts for setup and data ingestion
│   ├── setup_collections.py
│   ├── seed_literature.py
│   ├── seed_trials.py
│   ├── seed_findings.py
│   ├── seed_protocols.py
│   ├── seed_devices.py
│   ├── seed_anatomy.py
│   ├── seed_benchmarks.py
│   ├── seed_guidelines.py
│   ├── seed_report_templates.py
│   ├── seed_datasets.py
│   ├── ingest_pubmed.py
│   ├── ingest_clinical_trials.py
│   ├── test_nim_services.py
│   ├── test_rag_pipeline.py
│   ├── validate_e2e.py
│   ├── generate_sample_dicoms.py
│   ├── generate_annotated_images.py
│   ├── download_real_data.py
│   └── validate_real_data.py
├── tests/                            # 620 tests across 12 modules
│   ├── conftest.py
│   ├── test_knowledge.py
│   ├── test_query_expansion.py
│   ├── test_cross_modal.py
│   ├── test_export.py
│   ├── test_workflows.py
│   ├── test_nim_clients.py
│   ├── test_dicom_ingestion.py
│   ├── test_preview_api.py
│   └── test_preview_generator.py
├── docs/                             # Documentation
├── requirements.txt
├── Dockerfile
├── docker-compose.yml
├── .streamlit/
└── LICENSE                           # Apache 2.0
```

### File-by-file walkthrough of the core files

#### config/settings.py (158 lines)

The single source of truth for all configuration. Uses Pydantic `BaseSettings` so every value can be overridden via environment variables (prefixed with `IMAGING_`) or a `.env` file.

Key settings:
- `MILVUS_HOST` / `MILVUS_PORT`: Vector database connection (default: localhost:19530)
- `EMBEDDING_MODEL`: BGE-small-en-v1.5 (384 dimensions)
- `LLM_MODEL`: Claude Sonnet 4.6 (via Anthropic API), with Llama 3 8B (NIM) as fallback
- `TOP_K_PER_COLLECTION`: Maximum results per collection per query (default: 5)
- `SCORE_THRESHOLD`: Minimum cosine similarity score (default: 0.4)
- `WEIGHT_*`: Per-collection search weights
- `CITATION_HIGH_THRESHOLD` / `CITATION_MEDIUM_THRESHOLD`: Score thresholds for relevance tagging (0.75 / 0.60)
- NIM service endpoints and ports

#### src/models.py (728 lines)

All Pydantic data models used throughout the system. Contains:

- **Enums** for modality types, body regions, AI tasks, severity levels, regulatory statuses, trial phases, and evidence levels
- **10 collection models** (one per owned collection): `ImagingLiterature`, `ImagingTrial`, `ImagingFinding`, `ImagingProtocol`, `ImagingDevice`, `ImagingAnatomy`, `ImagingBenchmark`, `ImagingGuideline`, `ImagingReportTemplate`, `ImagingDataset`
- **Search result models:** `SearchHit`, `CrossCollectionResult`, `ComparativeResult`
- **Agent models:** `AgentQuery`, `AgentResponse`
- **Workflow models:** workflow input and output schemas

Each collection model has a `to_embedding_text()` method that generates the text string used for BGE-small embedding. The quality of the embedding depends on how the text is constructed from the structured fields.

#### src/collections.py (1,420 lines)

Manages the 11 Milvus collections. Contains:

- **Schema definitions** for all 11 collections, defining every field's name, data type, and constraints
- **COLLECTION_SCHEMAS** registry mapping collection names to their schemas
- **ImagingCollectionManager** class with methods for:
  - `connect()` / `disconnect()`: Milvus connection lifecycle
  - `create_collection()` / `create_all_collections()`: Schema creation with IVF_FLAT index (384 dim, COSINE)
  - `get_collection_stats()`: Record counts per collection
  - `insert_batch()`: Bulk data insertion with embedding
  - `search()`: Single-collection vector similarity search
  - `search_all()`: Parallel search across ALL collections using `ThreadPoolExecutor`

The `search_all()` method is the workhorse of the system. It launches concurrent searches across all 11 collections and merges the results.

#### src/rag_engine.py (690 lines)

The multi-collection RAG engine. This is the central piece of the system. It orchestrates:

1. **Query embedding** (with BGE asymmetric query prefix)
2. **Parallel collection search** (via the collection manager)
3. **Query expansion** (via the expansion module)
4. **Knowledge graph augmentation** (from all 3 domains: pathologies, modalities, anatomy)
5. **Result merging and ranking** (deduplication, weighted scoring, relevance tagging)
6. **Prompt construction** (formatting evidence, knowledge context, and the question into a structured LLM prompt)
7. **LLM generation** (synchronous and streaming modes)
8. **Comparative analysis** (detecting "vs" queries like "CT vs MRI," parsing entities, running dual retrievals)
9. **Cross-collection entity linking** (`find_related()` method)

#### src/knowledge.py (1,843 lines)

The knowledge graph. Contains three large Python dictionaries:

- `IMAGING_PATHOLOGIES`: 25 pathology entries with ICD-10 codes, imaging characteristics, severity criteria, and NIM workflows
- `IMAGING_MODALITIES`: 9 modality entries with physics, acquisition parameters, and AI applications
- `IMAGING_ANATOMY`: 21 anatomy entries with SNOMED codes, FMA identifiers, and VISTA-3D segmentation labels

Public API functions:
- `get_pathology_context()`: Formatted text for a pathology
- `get_modality_context()`: Formatted text for an imaging modality
- `get_anatomy_context()`: Formatted text for an anatomical structure
- `get_nim_recommendation()`: Suggests which NIM service to invoke for a given pathology
- `resolve_comparison_entity()`: Resolves a raw text string to a known entity for comparative analysis
- `get_comparison_context()`: Builds side-by-side knowledge for two entities
- `get_knowledge_stats()`: Returns counts of pathologies, modalities, and anatomy entries

#### src/agent.py (283 lines)

The autonomous Imaging Intelligence Agent. Implements the plan-search-synthesize pattern:

1. **`search_plan()`**: Analyzes the question to identify modalities, body regions, pathologies, and search strategy (broad, targeted, or comparative).
2. **`evaluate_evidence()`**: Assesses evidence quality as "sufficient," "partial," or "insufficient" based on collection coverage and hit count.
3. **`run()`**: Orchestrates the full pipeline: plan, search, evaluate, optionally expand with sub-questions, generate answer.
4. **`generate_report()`**: Creates a formatted markdown report from the results.

#### src/cross_modal.py (401 lines)

Cross-modal query engine that bridges imaging and genomic evidence. Enables questions like "Are there genomic markers associated with the imaging phenotype of ground-glass opacity in lung adenocarcinoma?" by:

1. Extracting imaging features from the question
2. Querying the 10 imaging collections for relevant evidence
3. Translating imaging features into genomic search terms
4. Querying the shared `genomic_evidence` collection
5. Merging and correlating the cross-modal results

#### src/export.py (691 lines)

Export engine for generating downloadable reports in multiple formats:

- **Markdown export**: Formatted `.md` file with query, response, evidence, and metadata
- **JSON export**: Structured `.json` file for programmatic processing
- **PDF export**: Formatted PDF report using ReportLab with headers, tables, and citations

#### app/imaging_ui.py (entry point)

The Streamlit application entry point. Configures the page layout, loads the dark theme with NVIDIA styling, initializes the RAG engine and NIM service manager, and renders the nine tabs by delegating to the tab modules in `app/tabs/`.

### How the pieces connect: a trace through the code

When you type "What Lung-RADS score should be assigned to a 10mm part-solid nodule?" in the Evidence Explorer tab:

1. `app/tabs/evidence_explorer.py` receives the input
2. The tab calls `engine.retrieve()` which:
   - Embeds the query via BGE-small-en-v1.5
   - Searches all 11 collections in parallel via `collections.search_all()`
   - Runs `query_expansion.expand_query()` which detects "lung-rads" and expands with Lung-RADS scoring terms, and detects "nodule" expanding with "pulmonary nodule", "lung nodule", "SPN", etc.
3. The RAG engine detects "lung_nodule" in the knowledge graph and adds structured context:
   - Severity criteria: "Lung-RADS 4A, size 8-15 mm or growing" for a 10mm part-solid nodule
   - Key measurements: long_axis_mm, short_axis_mm, volume_mm3, doubling_time_days
   - CT characteristics: "Part-solid: mixed"
4. Evidence from `imaging_guidelines` (Fleischner Society, ACR), `imaging_findings` (structured finding descriptions), and `imaging_literature` (published Lung-RADS validation studies) is merged and ranked
5. The LLM receives the combined evidence and knowledge context and generates a grounded answer citing specific guideline recommendations and literature evidence
6. The response streams to the UI with citation links and relevance indicators

---

## Chapter 12: Next Steps

### Where to go from here

Now that you understand the foundations, here are pathways for further learning:

**If you are a radiologist / clinician:**
- Try asking the system questions from your daily practice (e.g., "What follow-up is recommended for a BI-RADS 4A lesion on screening mammography?")
- Explore the Workflow Runner tab to see how the system handles end-to-end clinical scenarios
- Use the Protocol Advisor tab to compare acquisition protocols for your clinical indications
- Use the Device Ecosystem tab to discover FDA-cleared AI tools for your subspecialty
- Load the demo cases in Patient 360 to see how cross-collection evidence is integrated for a single patient

**If you are a data scientist / ML engineer:**
- Read `src/rag_engine.py` to understand the retrieval-augmentation-generation pipeline in detail
- Experiment with collection weights in `config/settings.py` to see how they affect result ranking
- Look at the query expansion maps in `src/query_expansion.py` and consider what additional domain maps might improve recall
- Explore the NIM client implementations in `src/nim/` to understand how VISTA-3D, MAISI, and VILA-M3 are integrated
- Use the Benchmarks tab to compare model performance and identify gaps in current AI capabilities
- Study `src/cross_modal.py` to understand how imaging-to-genomic queries are constructed

**If you are a software developer:**
- Read `src/collections.py` to understand Milvus schema design and IVF_FLAT index configuration
- Look at `api/main.py` and the route modules in `api/routes/` to see how the 19+ endpoints are implemented
- Explore the ingest parsers in `src/ingest/` to understand how data flows from PubMed, ClinicalTrials.gov, and DICOM files into the vector database
- Study the workflow implementations in `src/workflows/` to understand how clinical pipelines are orchestrated
- Run the test suite: `pytest tests/ -v` (620 tests across 12 modules)
- Read `src/export.py` to understand PDF report generation with ReportLab

**If you are a healthcare IT administrator:**
- Review `config/settings.py` to understand all configurable parameters
- Study `docker-compose.yml` for containerized deployment options
- Look at the `IMAGING_` environment variable prefix for deployment configuration
- Review the port assignments: Streamlit 8525, FastAPI 8524, Milvus 19530

### How to contribute

The Imaging Intelligence Agent is open-source under the Apache 2.0 license. Contributions are welcome in these areas:

- **New data sources**: Additional parsers for FDA device databases, DICOM SR, RSNA competition results
- **New collections**: Specialized collections for contrast agents, radiation physics, or quality metrics
- **Knowledge graph expansion**: Adding new pathology entries, modality details, or anatomy structures
- **Query expansion refinement**: Improving term coverage for subspecialty domains (neuroradiology, MSK, cardiac)
- **New workflows**: Additional clinical workflow pipelines for other modalities and body regions
- **NIM integration**: Connecting additional NVIDIA NIM services as they become available
- **UI enhancements**: Better visualization, DICOM viewer integration, measurement tools
- **Performance optimization**: Faster search, better caching, reduced latency
- **Testing**: Additional tests to improve coverage (currently 620 tests across 12 modules)

### Resources for learning more

**Medical imaging AI:**
- Rajpurkar et al., "AI in Health and Medicine," Nature Medicine 2022 (comprehensive review of clinical AI)
- Liu et al., "How to Read Articles That Use Machine Learning," JAMA 2019 (critical appraisal guide for clinicians)
- ACR Data Science Institute: https://www.acrdsi.org/ (FDA-cleared AI device directory)
- The MONAI Project: https://monai.io/ (open-source framework for medical imaging AI)

**Radiology scoring systems:**
- Lung-RADS: ACR Lung-RADS Assessment Categories v2022
- BI-RADS: ACR BI-RADS Atlas, 5th Edition
- PI-RADS: PI-RADS v2.1 (2019)
- LI-RADS: ACR LI-RADS v2018
- CAD-RADS: CAD-RADS 2.0 (2022)
- TI-RADS: ACR TI-RADS (2017)
- ASPECTS: Barber et al., Lancet 2000

**Vector databases and embeddings:**
- Milvus documentation: https://milvus.io/docs
- BGE embedding models: https://huggingface.co/BAAI/bge-small-en-v1.5
- "Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks" (Lewis et al., 2020) -- the original RAG paper

**NVIDIA NIM:**
- NVIDIA NIM documentation: https://developer.nvidia.com/nim
- VISTA-3D: Interactive 3D medical image segmentation
- MAISI: Medical AI Synthetic Imaging
- VILA-M3: Vision-Language model for medical imaging

**This project:**
- HCLS AI Factory GitHub: https://github.com/ajones1923/hcls-ai-factory
- Architecture guide: `docs/ARCHITECTURE_GUIDE.md` (in this repository)
- NIM integration guide: `docs/NIM_INTEGRATION_GUIDE.md` (in this repository)
- Project bible: `docs/PROJECT_BIBLE.md` (in this repository)

---

## Glossary

| Term | Definition |
|------|-----------|
| **510(k)** | An FDA premarket notification pathway. A device manufacturer demonstrates that their product is substantially equivalent to an existing legally marketed device (the predicate device). Most radiology AI products are cleared via 510(k). |
| **ALARA** | As Low As Reasonably Achievable. A radiation safety principle requiring that radiation exposure be minimized while maintaining diagnostic image quality. |
| **ASPECTS** | Alberta Stroke Program Early CT Score. A 10-point scoring system for quantifying early ischemic changes on non-contrast CT of the brain. A score of 10 is normal; lower scores indicate more extensive ischemia. Used to guide thrombolysis and thrombectomy decisions. |
| **BGE-small-en-v1.5** | The embedding model used by this system. Produced by BAAI (Beijing Academy of Artificial Intelligence). Converts text to 384-dimensional vectors optimized for semantic similarity search. |
| **BI-RADS** | Breast Imaging Reporting and Data System. A standardized classification system for breast imaging findings (mammography, ultrasound, MRI). Categories range from 0 (incomplete) to 6 (known malignancy). Developed by the ACR. |
| **Body region** | An anatomical area of the body targeted by an imaging study. This system covers 16 body regions including head, neck, chest, abdomen, pelvis, spine, extremity, cardiac, breast, shoulder, and musculoskeletal. |
| **CAD-RADS** | Coronary Artery Disease Reporting and Data System. A standardized reporting system for coronary CT angiography. Categories range from 0 (no stenosis) to 5 (total occlusion). |
| **Classification** | In imaging AI, the task of assigning a category or score to an imaging finding (e.g., benign vs. malignant, Lung-RADS 1 vs. 4A). |
| **Collection** | In Milvus, a collection is a table that stores vectors and their associated metadata. This system has 11 collections (10 imaging-specific + 1 shared genomic). |
| **Contrast agent** | A substance administered to a patient (intravenous, oral, or rectal) to enhance the visibility of specific structures or pathologies on imaging. Common agents include iodinated contrast (for CT) and gadolinium-based contrast (for MRI). |
| **Cosine similarity** | A measure of similarity between two vectors based on the angle between them. Ranges from 0 (completely different) to 1 (identical direction). Used to find relevant evidence in the vector database. |
| **Cross-modal query** | A query that bridges two data domains, such as searching imaging collections and genomic collections simultaneously. Enables radiogenomics research. |
| **CT** | Computed Tomography. An imaging modality that uses X-ray beams rotated around the patient to create cross-sectional images. Tissue density is measured in Hounsfield Units (HU). |
| **CXR** | Chest X-Ray. A plain radiograph of the chest, the most commonly performed imaging study worldwide. Used for pneumonia, pneumothorax, pleural effusion, and cardiac silhouette assessment. |
| **De Novo** | An FDA regulatory pathway for novel medical devices that have no substantially equivalent predicate. Creates a new device classification. |
| **Detection** | In imaging AI, the task of identifying the presence and location of an abnormality in a medical image (e.g., detecting a lung nodule on CT). |
| **DGX Spark** | An NVIDIA desktop computer with a GB10 GPU, 128 GB unified LPDDR5x memory, and 20 ARM cores. The hardware target for this system ($4,699). |
| **DICOM** | Digital Imaging and Communications in Medicine. The international standard for storing, transmitting, and displaying medical images. Every clinical imaging study produces DICOM files. |
| **Embedding** | A numerical representation of text as a vector (list of numbers) in a high-dimensional space. Similar texts have similar embeddings. This system uses 384-dimensional embeddings. |
| **FastAPI** | A modern Python web framework for building REST APIs. Used for the Imaging Intelligence Agent's API server on port 8524. |
| **FDA** | U.S. Food and Drug Administration. The regulatory agency that clears or approves medical devices, including AI/ML-based radiology software. |
| **Finding** | An observation documented by a radiologist when interpreting a medical image. Examples: "6mm solid nodule in right upper lobe," "subdural hematoma with 3mm midline shift." |
| **Finding category** | A classification of a radiology finding. This system supports 37 finding categories such as mass, nodule, fracture, hemorrhage, effusion, and pneumothorax. |
| **Fleischner Society** | An international society for thoracic radiology that publishes widely-used guidelines for managing incidental pulmonary nodules. |
| **FMA** | Foundational Model of Anatomy. A reference ontology for anatomy that assigns unique identifiers to anatomical structures. Used in the imaging_anatomy collection. |
| **Ground-glass opacity (GGO)** | A hazy area of increased attenuation on CT that does not obscure underlying structures. Can indicate inflammation, infection, or early malignancy. |
| **Hounsfield Unit (HU)** | The unit of measurement for CT density. Water = 0 HU, air = -1000 HU, bone = +1000 HU. Acute blood on CT is typically 50-70 HU. |
| **IVF_FLAT** | Inverted File with Flat quantization. The Milvus indexing algorithm used by this system. It partitions vectors into clusters (nlist) and searches the nearest clusters (nprobe) for fast approximate nearest neighbor search. |
| **Knowledge graph** | A structured database of facts and relationships. In this system, the knowledge graph contains curated data about 25 pathologies, 9 modalities, and 21 anatomy structures. |
| **LI-RADS** | Liver Imaging Reporting and Data System. A standardized system for interpreting and reporting liver imaging in patients at risk for hepatocellular carcinoma. Categories range from LR-1 (definitely benign) to LR-5 (definitely HCC). |
| **LLM** | Large Language Model. An AI model trained on vast amounts of text that can generate human-like responses. This system uses Llama 3 8B (via NVIDIA NIM) as its primary LLM. |
| **Lung-RADS** | Lung Imaging Reporting and Data System. A standardized classification system for low-dose CT lung cancer screening. Categories range from 0 (incomplete) to 4X (suspicious with additional features). |
| **MAISI** | Medical AI Synthetic Imaging. An NVIDIA NIM service for generating synthetic medical images for training and augmentation. |
| **Milvus** | An open-source vector database designed for similarity search on large-scale embedding datasets. This system uses Milvus 2.4 on port 19530. |
| **Modality** | An imaging technique or technology used to create medical images. This system covers 9 modalities: CT, MRI, X-ray, CXR, PET, PET/CT, ultrasound, mammography, and fluoroscopy. |
| **NIM** | NVIDIA Inference Microservices. Containerized AI model inference services. This system uses 4 NIM services: LLM (Llama 3 8B), VISTA-3D, MAISI, and VILA-M3. |
| **PI-RADS** | Prostate Imaging Reporting and Data System. A standardized scoring system for multiparametric MRI of the prostate. Categories range from 1 (very low suspicion) to 5 (very high suspicion for clinically significant cancer). |
| **PMA** | Premarket Approval. The most rigorous FDA regulatory pathway, required for Class III medical devices. Requires clinical evidence of safety and effectiveness. |
| **Protocol** | A standardized set of imaging acquisition parameters (e.g., slice thickness, kVp, mAs, contrast timing) optimized for a specific clinical indication. |
| **Pydantic** | A Python library for data validation using type annotations. Used throughout this system for data models, API schemas, and configuration. |
| **Query expansion** | The process of augmenting a search query with related terms to improve recall. For example, expanding "hemorrhage" to also search for "bleeding," "hematoma," and "ICH." |
| **RAG** | Retrieval-Augmented Generation. A technique that combines information retrieval (finding relevant documents) with language model generation (writing an answer based on those documents). |
| **Radiogenomics** | The study of correlations between imaging features (radiomics) and genomic data. Enabled by the cross-modal query capability bridging imaging and genomic evidence collections. |
| **Scoring system** | A standardized classification scheme used in radiology to categorize findings and guide clinical management. This system covers 7 scoring systems: Lung-RADS, BI-RADS, TI-RADS, LI-RADS, CAD-RADS, PI-RADS, and ASPECTS. |
| **Segmentation** | In imaging AI, the task of delineating the boundaries of anatomical structures or lesions in a medical image, typically at the pixel or voxel level. |
| **Severity** | A classification of finding urgency. This system uses 5 severity levels: critical, urgent, moderate, routine, and incidental. |
| **SNOMED CT** | Systematized Nomenclature of Medicine -- Clinical Terms. A comprehensive clinical terminology system. Used in the imaging_anatomy collection to standardize anatomical structure identification. |
| **Streamlit** | A Python framework for building data applications and dashboards. Used for the Imaging Intelligence Agent's web UI on port 8525. |
| **TI-RADS** | Thyroid Imaging Reporting and Data System. A standardized classification system for thyroid nodules detected on ultrasound. Categories range from TR1 (benign) to TR5 (highly suspicious). |
| **Triage** | In imaging AI, the automated prioritization of imaging studies based on urgency. For example, flagging a CT head with suspected ICH for immediate radiologist review. |
| **Vector** | In the context of this system, a list of numbers (384 floating-point values) that represents the meaning of a piece of text. |
| **Vector database** | A database optimized for storing and searching high-dimensional vectors. Enables fast similarity search across millions of embeddings. |
| **VILA-M3** | A multimodal vision-language model from NVIDIA designed for medical imaging. One of the 4 NIM services used by this system. |
| **VISTA-3D** | An NVIDIA NIM service for interactive 3D medical image segmentation. Supports multiple organ and lesion segmentation tasks. |
| **Workflow** | A predefined multi-step clinical pipeline that orchestrates evidence retrieval, scoring system application, guideline consultation, and report generation for a specific clinical scenario. This system includes 6 workflows. |

---

*Imaging Intelligence Agent | HCLS AI Factory | Apache 2.0 | Adam Jones | March 2026*
