# Ending the Diagnostic Odyssey: A Multi-Collection RAG Architecture for Clinical Document Intelligence and Genomic-Autoimmune Correlation

**Author:** Adam Jones
**Date:** March 2026
**Version:** 0.1.0 (Pre-Implementation)
**License:** Apache 2.0

Part of the HCLS AI Factory -- an end-to-end precision medicine platform.
https://github.com/ajones1923/hcls-ai-factory

---

## Abstract

Autoimmune diseases affect approximately 50 million Americans and 800 million people worldwide, yet the average time from symptom onset to definitive diagnosis remains staggeringly long: 4.6 years for systemic lupus erythematosus (SLE), 5-7 years for postural orthostatic tachycardia syndrome (POTS), 6-10 years for celiac disease, 4-5 years for multiple sclerosis, and over 10 years for conditions like Ehlers-Danlos syndrome and systemic mastocytosis. This "diagnostic odyssey" results from a fundamental structural problem: autoimmune diseases are multi-system disorders that generate fragmented clinical data across dozens of specialist encounters, hundreds of laboratory tests, imaging studies, procedure reports, and genetic evaluations -- yet no existing tool synthesizes this longitudinal patient journey into a coherent diagnostic picture.

This paper presents the architectural design, clinical rationale, and product requirements for the Precision Autoimmune Agent -- a clinical document intelligence system built on multi-collection retrieval-augmented generation (RAG) that ingests thousands of patient clinical documents from the complete medical journey, identifies patterns across laboratory results, specialist assessments, imaging findings, and genomic data, and cross-references these patterns against known autoimmune disease signatures to dramatically accelerate diagnosis. The agent will unify 14 specialized Milvus vector collections spanning clinical documents (progress notes, lab reports, imaging reports, pathology, procedure notes), autoimmune reference knowledge (literature, clinical trials, autoantibody databases, HLA-disease associations, disease activity indices), genomic data (3.5 million variants from the HCLS AI Factory genomics pipeline), and longitudinal biomarker tracking -- enabling queries like "What patterns in this patient's 3-year clinical history suggest an underlying autoimmune etiology?" or "Are there genomic variants in this patient that increase susceptibility to lupus nephritis?"

The system extends the proven multi-collection RAG architecture established by five existing intelligence agents in the HCLS AI Factory (Precision Biomarker, Precision Oncology, CAR-T, Imaging, and Autoimmune prototype), adapting it with clinical document ingestion pipelines capable of processing thousands of patient records, NLP-based entity extraction for laboratory values and clinical findings, longitudinal biomarker trend analysis, HLA-disease association scoring with 50+ known associations, autoantibody panel interpretation across 14 antibody types, disease activity scoring (DAS28, SLEDAI-2K, CDAI, BASDAI), flare prediction algorithms, and pharmacogenomic-guided biologic therapy recommendations for 8 drug classes. Eight reference clinical workflows will cover the highest-impact diagnostic and monitoring scenarios: diagnostic odyssey acceleration, lupus nephritis surveillance, POTS/dysautonomia evaluation, inflammatory arthritis differentiation, overlap syndrome detection, biologic therapy optimization, flare prediction and prevention, and genomic-autoimmune risk profiling.

The agent will deploy on a single NVIDIA DGX Spark ($3,999) using BGE-small-en-v1.5 embeddings (384-dimensional, IVF_FLAT, COSINE), Claude Sonnet 4.6 for evidence synthesis, and shared NVIDIA NIM microservices for on-device inference. Licensed under Apache 2.0, the platform will democratize access to integrated autoimmune intelligence that currently requires multi-million-dollar institutional investments in informatics infrastructure -- bringing the diagnostic power of world-class rheumatology and immunology centers to any clinic worldwide.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [The Autoimmune Diagnostic Crisis](#2-the-autoimmune-diagnostic-crisis)
3. [Clinical Landscape and Market Analysis](#3-clinical-landscape-and-market-analysis)
4. [Existing HCLS AI Factory Architecture](#4-existing-hcls-ai-factory-architecture)
5. [Precision Autoimmune Agent Architecture](#5-precision-autoimmune-agent-architecture)
6. [Clinical Document Ingestion Pipeline](#6-clinical-document-ingestion-pipeline)
7. [Milvus Collection Design](#7-milvus-collection-design)
8. [Clinical Workflows](#8-clinical-workflows)
9. [Cross-Modal Integration and Genomic Correlation](#9-cross-modal-integration-and-genomic-correlation)
10. [NIM Integration Strategy](#10-nim-integration-strategy)
11. [Knowledge Graph Design](#11-knowledge-graph-design)
12. [Query Expansion and Retrieval Strategy](#12-query-expansion-and-retrieval-strategy)
13. [API and UI Design](#13-api-and-ui-design)
14. [Clinical Decision Support Engines](#14-clinical-decision-support-engines)
15. [Reporting and Interoperability](#15-reporting-and-interoperability)
16. [Product Requirements Document](#16-product-requirements-document)
17. [Data Acquisition Strategy](#17-data-acquisition-strategy)
18. [Validation and Testing Strategy](#18-validation-and-testing-strategy)
19. [Regulatory Considerations](#19-regulatory-considerations)
20. [DGX Compute Progression](#20-dgx-compute-progression)
21. [Implementation Roadmap](#21-implementation-roadmap)
22. [Risk Analysis](#22-risk-analysis)
23. [Competitive Landscape](#23-competitive-landscape)
24. [Discussion](#24-discussion)
25. [Conclusion](#25-conclusion)
26. [References](#26-references)

---

## 1. Introduction

### 1.1 The Autoimmune Disease Burden

Autoimmune diseases represent one of the most significant and underappreciated categories of chronic illness worldwide. The National Institutes of Health estimates that over 80 distinct autoimmune conditions collectively affect approximately 50 million Americans -- more than heart disease (30 million) and cancer (18 million) combined. Globally, the World Health Organization estimates 800 million people live with autoimmune conditions, with prevalence rising at 3-9% annually in developed nations.

The clinical and economic burden is staggering:

- **Annual U.S. healthcare costs:** $100+ billion in direct medical expenses
- **Disability:** Autoimmune diseases are among the top 10 leading causes of death in women under 65
- **Workforce impact:** 25% of autoimmune patients report being unable to work during flares
- **Mental health:** 30-40% prevalence of comorbid depression and anxiety
- **Mortality:** SLE carries a standardized mortality ratio of 2.6; systemic sclerosis 3.5; systemic vasculitis 2.7

Yet despite this enormous burden, autoimmune diseases remain among the most difficult conditions to diagnose. The fundamental challenge is that autoimmune diseases are systemic -- they affect multiple organ systems simultaneously, producing constellations of symptoms that span the expertise of dozens of medical specialties. A lupus patient may see a rheumatologist for joint pain, a dermatologist for rashes, a nephrologist for kidney involvement, a hematologist for cytopenias, a neurologist for cognitive dysfunction, and a cardiologist for pericarditis -- each specialist generating clinical documentation in isolated electronic health record (EHR) silos that no existing system synthesizes into a unified diagnostic picture.

### 1.2 The Diagnostic Odyssey

The term "diagnostic odyssey" describes the years-long journey that autoimmune patients endure before receiving a correct diagnosis. Published data quantify this crisis:

| Condition | Average Time to Diagnosis | Specialists Seen | Source |
|-----------|--------------------------|-----------------|--------|
| Systemic lupus erythematosus (SLE) | 4.6 years | 3-5 | Lupus Foundation of America |
| Postural orthostatic tachycardia syndrome (POTS) | 5.9 years | 7+ | Dysautonomia International |
| Ehlers-Danlos syndrome (hEDS) | 10-12 years | 8+ | Ehlers-Danlos Society |
| Celiac disease | 6-10 years | 3-4 | Celiac Disease Foundation |
| Multiple sclerosis | 4-5 years | 2-4 | National MS Society |
| Ankylosing spondylitis | 6-8 years | 3-5 | Spondylitis Association |
| Systemic mastocytosis | 8-10 years | 5+ | The Mastocytosis Society |
| Myasthenia gravis | 2-3 years | 3-4 | Myasthenia Gravis Foundation |
| Sjogren's syndrome | 3-4 years | 2-4 | Sjogren's Foundation |
| Autoimmune encephalitis | 1-3 years | 4-6 | Autoimmune Encephalitis Alliance |

The consequences of delayed diagnosis are not merely inconvenient -- they are medically devastating. Each year of undiagnosed lupus nephritis increases the risk of irreversible kidney damage. Each year of untreated inflammatory arthritis produces measurable joint erosion. Each year of undiagnosed celiac disease compounds the risk of lymphoma, osteoporosis, and infertility. The diagnostic odyssey is not a quality-of-life issue; it is a patient safety crisis.

### 1.3 Why Clinical Document Intelligence Is the Solution

The data needed to diagnose most autoimmune patients already exists -- it is scattered across their medical records. A patient who has seen 7 specialists over 5 years may have generated:

- 50-200 progress notes documenting symptom patterns, physical exam findings, and clinical impressions
- 100-500 laboratory results tracking inflammatory markers, autoantibodies, complement levels, and organ function
- 10-30 imaging studies (X-rays, MRIs, CT scans, ultrasounds) with radiologist interpretations
- 5-15 procedure reports (biopsies, endoscopies, nerve conduction studies, tilt table tests)
- Genetic testing results (HLA typing, whole exome sequencing, pharmacogenomic panels)
- Medication histories documenting treatment responses and failures

No human clinician can read thousands of pages of clinical documentation, identify subtle patterns across years of data, and correlate those patterns with the 80+ known autoimmune diseases and their genomic risk factors. This is precisely the task that a multi-collection RAG system with clinical document ingestion, NLP entity extraction, and genomic cross-referencing is designed to solve.

### 1.4 Our Contribution

The Precision Autoimmune Agent addresses the diagnostic odyssey through a clinical document intelligence architecture that:

- Ingests **thousands of patient clinical documents** (progress notes, lab reports, imaging reports, pathology, procedures) via OCR and NLP pipelines into patient-specific vector collections
- Extracts **structured entities** from unstructured clinical text -- laboratory values, autoantibody results, imaging findings, medication changes, symptom patterns -- using medical NLP models
- Cross-references extracted patterns against **14 specialized Milvus collections** containing autoimmune reference knowledge, HLA-disease associations (50+ conditions), autoantibody databases (14 antibody-disease maps), and 3.5 million genomic variants
- Performs **longitudinal biomarker trend analysis** to detect rising inflammatory markers, falling complement levels, and other flare-predictive patterns that human review would miss
- Calculates **validated disease activity scores** (DAS28-CRP, DAS28-ESR, SLEDAI-2K, CDAI, BASDAI) from extracted clinical data
- Provides **pharmacogenomic-guided biologic therapy recommendations** for 8 drug classes based on HLA typing and genotype data
- Generates **diagnostic hypothesis reports** that synthesize thousands of clinical data points into prioritized differential diagnoses with supporting evidence citations
- Runs on a **single NVIDIA DGX Spark** ($3,999), democratizing access to diagnostic intelligence that currently requires multi-million-dollar institutional informatics platforms

---

## 2. The Autoimmune Diagnostic Crisis

### 2.1 Data Fragmentation as the Root Cause

The diagnostic odyssey is not primarily a knowledge gap -- rheumatologists and immunologists understand autoimmune diseases well. The crisis is fundamentally a **data synthesis problem**. Patient data is fragmented across:

1. **Multiple EHR systems:** Patients who see specialists at different health systems generate records in incompatible EHRs (Epic, Cerner, Allscripts, Meditech). Even within a single system, clinical notes, laboratory results, and imaging reports are stored in separate modules.

2. **Specialty-specific documentation:** Each specialist documents findings relevant to their domain. A dermatologist notes a malar rash but does not order complement levels. A nephrologist documents proteinuria but does not review the dermatology notes describing the rash. The pattern (rash + proteinuria + low complement = lupus) exists in the aggregate record but is never seen by any single provider.

3. **Temporal dispersion:** Autoimmune symptoms evolve over months to years. A laboratory finding from 18 months ago (mildly positive ANA at 1:160) takes on completely different significance when combined with a physical exam finding from last week (Raynaud's phenomenon) and a lab result from last month (elevated anti-centromere antibody). No existing tool performs this temporal pattern recognition across the full clinical timeline.

4. **Laboratory data in unstructured formats:** Many laboratory results arrive as free-text reports (especially autoantibody panels, HLA typing results, pathology reports, and genetic testing). Extracting structured values from these reports requires medical NLP.

5. **Genomic data in separate systems:** When patients undergo HLA typing, whole exome sequencing, or pharmacogenomic testing, these results typically reside in genetic testing portals (Invitae, GeneDx, 23andMe, Color Health) disconnected from the clinical EHR.

### 2.2 The Overlapping Disease Problem

Autoimmune diseases frequently coexist, creating diagnostic complexity that exceeds the capacity of single-disease frameworks:

- **Overlap syndromes:** Mixed connective tissue disease (MCTD), overlap of SLE with systemic sclerosis or myositis, "rhupus" (RA + SLE overlap)
- **POTS as a comorbid condition:** POTS frequently co-occurs with Ehlers-Danlos syndrome (hEDS), mast cell activation syndrome (MCAS), small fiber neuropathy, and Sjogren's syndrome -- forming the POTS/hEDS/MCAS triad
- **Autoimmune polyendocrine syndromes:** Type 1 (Addison's + hypoparathyroidism + chronic mucocutaneous candidiasis), Type 2 (Addison's + autoimmune thyroid disease + Type 1 diabetes)
- **Familial clustering:** First-degree relatives of autoimmune patients have 2-10x higher risk of developing any autoimmune disease, not just the same disease

These overlapping patterns require a system that can simultaneously evaluate evidence for multiple autoimmune conditions and identify syndromic clusters that human pattern recognition would miss.

### 2.3 Why Existing Tools Fall Short

Current approaches to autoimmune diagnostics fail because they address only one dimension of the problem:

1. **Classification criteria tools** (ACR/EULAR calculators) require structured input that clinicians must manually extract from records -- they cannot ingest raw clinical documents.
2. **Laboratory reference systems** (UpToDate, DynaMed) provide disease-level reference information but cannot analyze individual patient data.
3. **EHR analytics** (Epic Cogito, Cerner HealtheIntent) offer population-level dashboards but lack the medical NLP and genomic integration needed for individual diagnostic reasoning.
4. **General AI assistants** lack the structured autoimmune knowledge graphs, validated scoring systems, and clinical document ingestion pipelines required for rigorous diagnostic support.
5. **Commercial autoimmune panels** (Exagen AVISE, Labcorp ARUP) test for specific autoantibodies but do not integrate results with clinical history, genomic data, or longitudinal biomarker trends.

### 2.4 The Case for Clinical Document Intelligence

The Precision Autoimmune Agent represents a fundamentally different approach: rather than asking clinicians to manually input structured data into diagnostic calculators, the system ingests the patient's raw clinical documentation and extracts the patterns itself. This approach has three critical advantages:

1. **Completeness:** The system analyzes all available data, not just what a single provider remembers or chooses to enter. A positive ANA from 3 years ago, a borderline complement level from 18 months ago, and a new arthritis presentation from last week are all considered simultaneously.

2. **Pattern recognition at scale:** A human clinician reviewing 200 clinical documents might miss the subtle temporal pattern of rising anti-dsDNA titers preceding each lupus flare. The system identifies these patterns algorithmically.

3. **Genomic correlation:** By cross-referencing extracted clinical findings with the 3.5 million genomic variants in the shared `genomic_evidence` collection, the system can identify genetic risk factors (HLA alleles, STAT4 variants, IRF5 polymorphisms) that strengthen or weaken specific diagnostic hypotheses.

---

## 3. Clinical Landscape and Market Analysis

### 3.1 Autoimmune AI Market

The autoimmune diagnostics and AI market is experiencing rapid growth driven by rising prevalence, diagnostic delay awareness, and precision medicine adoption:

| Metric | Value | Source |
|--------|-------|--------|
| Global autoimmune diagnostics market (2025) | $5.2 billion | Grand View Research |
| Projected market (2030) | $8.7 billion | Grand View Research |
| CAGR | 10.8% | Grand View Research |
| AI in clinical diagnostics (2025) | $2.1 billion | MarketsandMarkets |
| AI in clinical diagnostics (2030) | $6.8 billion | MarketsandMarkets |
| U.S. autoimmune disease prevalence | 50 million | NIH/AARDA |
| Global autoimmune prevalence | 800 million | WHO |
| Average diagnostic delay (all autoimmune) | 4.5 years | AARDA |

### 3.2 Competitive Analysis

| Solution | Approach | Limitations |
|----------|----------|-------------|
| Exagen AVISE | Lupus-specific autoantibody panel with cell-bound complement | Single-disease focus; no clinical document ingestion; no genomics |
| DxTerity AutoImmune Profile | At-home RNA expression testing | Screening only; no longitudinal integration; no imaging/clinical notes |
| IBM Watson for Genomics | Genomic variant interpretation | Discontinued; no autoimmune-specific workflows; no clinical document NLP |
| Google Health DeepMind | Medical image analysis | Imaging only; no laboratory, clinical note, or genomic integration |
| Epic Cognitive Computing | EHR-integrated NLP | Single-EHR vendor; limited autoimmune-specific knowledge; no genomic correlation |
| **Precision Autoimmune Agent** | **Multi-collection RAG with clinical document ingestion, genomic correlation, 14 collections** | **Open-source; $3,999 hardware; cross-system document ingestion** |

### 3.3 Target Users

| User Segment | Estimated Size | Primary Value |
|-------------|---------------|---------------|
| Academic rheumatology centers | 150+ U.S. programs | Diagnostic acceleration, research data extraction |
| Community rheumatologists | 5,500 U.S. practitioners | Complex case support, genomic interpretation |
| Primary care physicians | 250,000+ U.S. PCPs | Early autoimmune detection, appropriate referral |
| Immunology researchers | 15,000+ globally | Literature synthesis, clinical-genomic correlation |
| Patient advocacy organizations | 100+ autoimmune foundations | Patient empowerment, second opinion support |
| Rare disease diagnostic centers | 50+ NIH UDP, UDNI | Undiagnosed autoimmune case evaluation |
| Pharma/biotech (autoimmune pipeline) | 200+ companies | Clinical trial design, biomarker identification |

---

## 4. Existing HCLS AI Factory Architecture

### 4.1 Platform Overview

The HCLS AI Factory is a three-stage precision medicine pipeline that processes a patient sample from raw DNA sequencing data to drug candidate molecules in under 5 hours on a single NVIDIA DGX Spark:

| Stage | Duration | Process | Output |
|-------|----------|---------|--------|
| **1. Genomics** | 120-240 min | FASTQ to VCF via Parabricks 4.6, BWA-MEM2, DeepVariant | 11.7 million variants |
| **2. RAG/Chat** | Interactive | ClinVar (~2.7M) + AlphaMissense (71M) + Milvus (3.5M vectors) + Claude | Variant interpretation, gene-drug associations |
| **3. Drug Discovery** | 8-16 min | MolMIM generation + DiffDock docking + RDKit scoring | Candidate molecules with docking scores |

Five intelligence agents currently extend Stage 2:

| Agent | Collections | Vectors | Port (API/UI) | Focus |
|-------|------------|---------|---------------|-------|
| Precision Biomarker | 10 | ~890 | 8510/8511 | Biomarker interpretation |
| Precision Oncology | 10 | ~950 | 8514/8515 | Cancer variant analysis |
| CAR-T Intelligence | 11 | 3,567,436 | 8521/8522 | Cell therapy development |
| Imaging Intelligence | 10 | 876 | 8524/8525 | Medical imaging AI |
| Autoimmune (prototype) | 0 | 0 | — | Agent class only (331 lines) |

The Precision Autoimmune Agent will be the sixth full intelligence agent, with the unique addition of a clinical document ingestion pipeline that enables patient-specific analysis at scale.

### 4.2 Shared Infrastructure

All intelligence agents share core infrastructure:

- **Vector database:** Milvus 2.4 on localhost:19530 with etcd and MinIO backing
- **Embedding model:** BGE-small-en-v1.5 (384-dimensional, IVF_FLAT index, COSINE similarity)
- **LLM:** Claude Sonnet 4.6 (primary) with Llama-3 8B NIM fallback
- **Shared collection:** `genomic_evidence` (3,561,170 variants) -- read-only access from all agents
- **Monitoring:** Prometheus + Grafana dashboards
- **Configuration:** Pydantic BaseSettings pattern

### 4.3 Proven Patterns Adapted for Autoimmune

The Precision Autoimmune Agent extends established patterns:

- **Multi-collection parallel search:** ThreadPoolExecutor dispatching to 14 collections simultaneously (adapted from CAR-T's 11-collection architecture)
- **Knowledge graph augmentation:** 7-dictionary knowledge graph (HLA associations, autoantibodies, disease activity, biologics, flare biomarkers, overlap syndromes, dysautonomia) -- adapted from CAR-T's 6-dictionary pattern
- **Query expansion:** 18 domain-specific expansion maps (250+ keywords to 2,000+ terms) -- adapted from CAR-T's 12-map pattern
- **Comparative analysis:** Auto-detection and resolution of "X vs Y" queries (e.g., "Compare RA vs lupus arthritis")
- **Multi-format export:** Markdown, JSON, PDF with NVIDIA branding, FHIR R4 DiagnosticReport

**New capability unique to this agent:** Clinical document ingestion pipeline with medical NLP entity extraction, enabling patient-specific analysis from raw clinical records rather than manually structured input.

---

## 5. Precision Autoimmune Agent Architecture

### 5.1 System Diagram

```
Patient Clinical Documents (PDFs, HL7 FHIR, CCDA, free text)
    |
    v
[Clinical Document Ingestion Pipeline]
  - OCR (Tesseract / Azure Document Intelligence)
  - Section segmentation (clinical notes, labs, imaging, pathology)
  - Medical NLP entity extraction (laboratory values, diagnoses, medications)
  - Temporal normalization (date extraction and timeline construction)
    |
    v
[Patient Document Collection] (patient-specific vectors in Milvus)
    |
    v
[BGE-small-en-v1.5 Embedding]
(384-dim, asymmetric query prefix)
    |
    v
[Parallel Search: 14 Milvus Collections]
(ThreadPoolExecutor, IVF_FLAT / COSINE)
  - 3 patient-specific collections (clinical_documents, patient_labs, patient_timeline)
  - 8 reference collections (literature, trials, autoantibodies, hla_associations,
    disease_activity, biologics, overlap_syndromes, guidelines)
  - 2 shared collections (genomic_evidence, biomarker_reference)
  - 1 cross-agent collection (imaging_findings)
    |
    v
[Query Expansion] (18 maps, 250+ keywords -> 2,000+ terms)
    |
    v
[Knowledge Graph Augmentation]
(50+ HLA associations, 14 autoantibody maps, 5 disease activity indices,
 8 biologic therapies, 3 flare biomarker patterns, 12 overlap syndromes,
 10 dysautonomia conditions)
    |
    v
[Longitudinal Pattern Analysis]
  - Biomarker trend detection (rising CRP, falling complement)
  - Autoantibody seroconversion tracking
  - Symptom pattern recognition across clinical notes
  - Medication response correlation
    |
    v
[Score-Weighted Merge & Rank]
(citation relevance: high >= 0.75, medium >= 0.60)
    |
    v
[Claude Sonnet 4.6] --> Grounded diagnostic hypotheses with evidence citations
    |
    v
[Export] --> Markdown | JSON | PDF | FHIR R4 DiagnosticReport
```

### 5.2 Design Principles

1. **Patient-centric:** Collections are organized around the patient's complete clinical journey, not around individual tests or encounters
2. **Longitudinal:** Temporal relationships between clinical events are preserved and queryable
3. **Genomic-first:** Every diagnostic hypothesis is automatically cross-referenced against the patient's genomic data (if available)
4. **Evidence-grounded:** Every claim cites specific clinical documents, laboratory values, or reference literature
5. **Multi-disease aware:** The system evaluates evidence for all 13+ supported autoimmune conditions simultaneously, detecting overlaps and polyautoimmunity
6. **Privacy-preserving:** Patient documents remain local on the DGX Spark; no clinical data is sent to cloud APIs (only anonymized queries reach the LLM)

### 5.3 Port Allocation

| Service | Port | Protocol |
|---------|------|----------|
| FastAPI (REST endpoints) | 8530 | HTTP |
| Streamlit UI | 8531 | HTTP |
| Document ingestion webhook | 8532 | HTTP |
| Shared Milvus | 19530 | gRPC |
| Shared etcd | 2379 | gRPC |
| Shared MinIO | 9000 | HTTP |

---

## 6. Clinical Document Ingestion Pipeline

### 6.1 The Core Innovation

The clinical document ingestion pipeline is the defining capability of the Precision Autoimmune Agent. Unlike other agents that work with curated reference data, this agent ingests a patient's actual clinical records -- potentially thousands of documents spanning years of care -- and transforms them into searchable, structured vectors that can be queried alongside reference knowledge.

### 6.2 Ingestion Architecture

```
Input Sources:
  - PDF clinical documents (scanned or digital)
  - HL7 FHIR R4 Bundles (from patient portals, EHR exports)
  - C-CDA documents (Consolidated Clinical Document Architecture)
  - Free-text clinical notes (typed or dictated)
  - Laboratory result files (HL7 ORU messages, CSV exports)
  - Genetic testing reports (PDF from Invitae, GeneDx, Color Health)

Processing Pipeline:
  Step 1: Document Classification
    - Classify each document by type: progress_note, lab_report, imaging_report,
      pathology_report, procedure_note, genetic_report, medication_list, referral
    - Use document header patterns, section markers, and content heuristics

  Step 2: OCR and Text Extraction
    - Digital PDFs: direct text extraction via PyMuPDF
    - Scanned PDFs: OCR via Tesseract with medical vocabulary enhancement
    - HL7/FHIR: structured field extraction via fhir.resources library
    - C-CDA: XML parsing with section-level extraction

  Step 3: Section Segmentation
    - Identify clinical note sections: Chief Complaint, HPI, Review of Systems,
      Physical Exam, Assessment/Plan, Laboratory Results, Imaging
    - Segment laboratory reports into individual test results
    - Extract structured findings from imaging report impressions

  Step 4: Medical NLP Entity Extraction
    - Laboratory values: extract test name, numeric value, unit, reference range, flag
    - Medications: extract drug name, dose, frequency, route, start/stop dates
    - Diagnoses: extract ICD-10 codes, problem descriptions, laterality
    - Symptoms: extract symptom descriptions, severity, duration, anatomic location
    - Vital signs: extract HR, BP, temperature, weight, orthostatic measurements
    - Autoantibodies: extract antibody name, titer, pattern, interpretation

  Step 5: Temporal Normalization
    - Extract document dates and encounter dates
    - Construct patient timeline with all clinical events
    - Normalize date formats across different source systems
    - Calculate intervals between events (onset-to-diagnosis, treatment-to-response)

  Step 6: Embedding and Indexing
    - Generate embedding text from structured entities + context
    - Embed with BGE-small-en-v1.5 (384-dim)
    - Insert into patient-specific Milvus collections with full metadata
```

### 6.3 Entity Extraction Models

The NLP pipeline uses a combination of approaches optimized for medical text:

| Component | Model/Approach | Purpose |
|-----------|---------------|---------|
| Named entity recognition | SciSpaCy (en_core_sci_lg) | Clinical entity detection |
| Negation detection | NegEx algorithm | Filter negated findings ("no rash") |
| Laboratory value extraction | Regex + rule-based parser | Structured lab result extraction |
| Medication extraction | MedEx-UIMA adapted patterns | Drug/dose/frequency parsing |
| Temporal expression | SUTime-based parser | Date and duration normalization |
| Section segmentation | SecTag algorithm patterns | Clinical note section identification |
| Abbreviation expansion | UMLS abbreviation dictionary | Medical abbreviation resolution |

### 6.4 Privacy and Security

Clinical document ingestion introduces significant privacy considerations:

- **All processing is local:** Document ingestion, NLP extraction, embedding, and Milvus indexing occur entirely on the DGX Spark. No patient data leaves the device.
- **LLM query anonymization:** When the system sends a query to Claude for synthesis, the prompt contains only the question and retrieved evidence snippets -- not raw clinical documents. Patient identifiers are stripped from evidence text before LLM submission.
- **Collection isolation:** Each patient's clinical documents are stored in a patient-specific Milvus partition, ensuring query isolation between patients.
- **Encryption at rest:** Milvus data on the DGX Spark NVMe SSD uses LUKS full-disk encryption.
- **Access control:** API endpoints for document ingestion and patient queries require JWT authentication with role-based access control (clinician, researcher, administrator).
- **Audit logging:** All document ingestion events and patient queries are logged with timestamps, user identity, and document identifiers for HIPAA compliance.

---

## 7. Milvus Collection Design

### 7.1 Index Configuration

All collections share the same index configuration:

| Parameter | Value |
|-----------|-------|
| Embedding model | BGE-small-en-v1.5 |
| Dimensions | 384 |
| Index type | IVF_FLAT |
| nlist | 1024 |
| Metric | COSINE |
| nprobe (search) | 16 |

### 7.2 Collection Schemas

**Collection 1: `autoimmune_clinical_documents`** -- Patient clinical document chunks

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique chunk identifier |
| patient_id | VARCHAR(32) | Patient identifier |
| document_type | VARCHAR(32) | progress_note, lab_report, imaging_report, pathology, procedure, genetic_report |
| document_date | VARCHAR(10) | ISO-8601 date |
| provider_specialty | VARCHAR(32) | Specialty of documenting provider |
| section | VARCHAR(32) | Note section (hpi, ros, exam, assessment, plan, labs, imaging) |
| text_chunk | VARCHAR(4096) | Document text chunk |
| extracted_entities | VARCHAR(2048) | JSON-encoded extracted entities |
| source_system | VARCHAR(32) | EHR/source system identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 2: `autoimmune_patient_labs`** -- Structured laboratory results extracted from clinical documents

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique result identifier |
| patient_id | VARCHAR(32) | Patient identifier |
| test_name | VARCHAR(64) | Laboratory test name (standardized) |
| loinc_code | VARCHAR(16) | LOINC code for standardized identification |
| value | FLOAT | Numeric result value |
| unit | VARCHAR(16) | Unit of measurement |
| reference_low | FLOAT | Lower reference range |
| reference_high | FLOAT | Upper reference range |
| flag | VARCHAR(8) | H (high), L (low), A (abnormal), N (normal) |
| collection_date | VARCHAR(10) | ISO-8601 date |
| text_context | VARCHAR(512) | Surrounding clinical context |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 3: `autoimmune_patient_timeline`** -- Longitudinal patient events for temporal analysis

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique event identifier |
| patient_id | VARCHAR(32) | Patient identifier |
| event_type | VARCHAR(32) | symptom_onset, diagnosis, lab_result, medication_change, procedure, flare, hospitalization |
| event_date | VARCHAR(10) | ISO-8601 date |
| event_description | VARCHAR(512) | Event description |
| severity | VARCHAR(16) | mild, moderate, severe, critical |
| associated_diagnoses | VARCHAR(256) | Related diagnoses |
| associated_labs | VARCHAR(512) | Related laboratory values |
| text_summary | VARCHAR(1024) | Narrative summary |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 4: `autoimmune_literature`** -- Published autoimmune research

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | PMID or unique identifier |
| title | VARCHAR(256) | Article title |
| text_chunk | VARCHAR(4096) | Abstract or text chunk |
| disease_category | VARCHAR(64) | Primary autoimmune disease |
| study_type | VARCHAR(32) | meta_analysis, rct, cohort, case_report, review |
| year | INT16 | Publication year |
| journal | VARCHAR(64) | Journal name |
| keywords | VARCHAR(512) | Author keywords |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 5: `autoimmune_trials`** -- Clinical trial registrations

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | NCT number |
| title | VARCHAR(256) | Trial title |
| text_summary | VARCHAR(2048) | Brief summary |
| phase | VARCHAR(8) | Phase 1, 2, 3, 4 |
| status | VARCHAR(32) | Recruiting, completed, active |
| disease | VARCHAR(64) | Target disease |
| intervention | VARCHAR(128) | Drug or intervention |
| sponsor | VARCHAR(64) | Lead sponsor |
| enrollment | INT32 | Target enrollment |
| start_year | INT16 | Year trial started |
| primary_endpoint | VARCHAR(256) | Primary outcome measure |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 6: `autoimmune_autoantibodies`** -- Autoantibody reference database

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| antibody_name | VARCHAR(32) | Autoantibody name (ANA, anti-dsDNA, RF, anti-CCP, etc.) |
| text_summary | VARCHAR(2048) | Clinical description |
| associated_diseases | VARCHAR(256) | Associated autoimmune conditions |
| sensitivity | FLOAT | Diagnostic sensitivity |
| specificity | FLOAT | Diagnostic specificity |
| pattern | VARCHAR(32) | Staining pattern (for ANA) |
| clinical_significance | VARCHAR(512) | Clinical interpretation guidance |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 7: `autoimmune_hla_associations`** -- HLA-disease association database

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| hla_allele | VARCHAR(16) | HLA allele (e.g., B*27:05, DRB1*04:01) |
| text_summary | VARCHAR(2048) | Association description |
| disease | VARCHAR(64) | Associated disease |
| odds_ratio | FLOAT | Odds ratio for disease risk |
| population | VARCHAR(32) | Study population |
| pmid | VARCHAR(16) | PubMed reference ID |
| mechanism | VARCHAR(512) | Proposed molecular mechanism |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 8: `autoimmune_disease_activity`** -- Disease activity scoring reference

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| score_name | VARCHAR(32) | DAS28-CRP, SLEDAI-2K, CDAI, BASDAI, etc. |
| text_summary | VARCHAR(2048) | Score description and interpretation |
| disease | VARCHAR(64) | Applicable disease |
| components | VARCHAR(512) | Score components |
| thresholds | VARCHAR(256) | Activity level thresholds (remission, low, moderate, high) |
| reference | VARCHAR(16) | PMID reference |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 9: `autoimmune_biologics`** -- Biologic therapy reference

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| drug_name | VARCHAR(32) | Drug name |
| drug_class | VARCHAR(32) | TNF inhibitor, IL-6 inhibitor, JAK inhibitor, etc. |
| text_summary | VARCHAR(2048) | Drug description and mechanism |
| mechanism | VARCHAR(256) | Mechanism of action |
| indicated_diseases | VARCHAR(256) | FDA-approved indications |
| pgx_considerations | VARCHAR(512) | Pharmacogenomic considerations |
| contraindications | VARCHAR(256) | Contraindications |
| monitoring | VARCHAR(256) | Required monitoring |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 10: `autoimmune_overlap_syndromes`** -- Overlap and polyautoimmunity reference

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| syndrome_name | VARCHAR(64) | Overlap syndrome name |
| text_summary | VARCHAR(2048) | Clinical description |
| component_diseases | VARCHAR(256) | Component diseases |
| diagnostic_criteria | VARCHAR(512) | Diagnostic criteria |
| key_autoantibodies | VARCHAR(128) | Characteristic antibodies |
| prevalence | VARCHAR(32) | Estimated prevalence |
| management | VARCHAR(512) | Treatment approach |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 11: `autoimmune_guidelines`** -- Clinical practice guidelines

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| guideline_title | VARCHAR(256) | Guideline title |
| text_chunk | VARCHAR(4096) | Guideline text chunk |
| issuing_body | VARCHAR(64) | ACR, EULAR, AGA, AAN, etc. |
| disease | VARCHAR(64) | Target disease |
| year | INT16 | Publication year |
| recommendation_strength | VARCHAR(16) | Strong, conditional, expert consensus |
| evidence_level | VARCHAR(8) | High, moderate, low, very low |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 12: `autoimmune_dysautonomia`** -- POTS/dysautonomia reference knowledge

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| condition | VARCHAR(64) | POTS, neurocardiogenic syncope, MSA, AAG, etc. |
| text_summary | VARCHAR(2048) | Clinical description |
| diagnostic_criteria | VARCHAR(512) | Diagnostic criteria |
| autonomic_tests | VARCHAR(256) | Recommended autonomic testing |
| comorbidities | VARCHAR(256) | Common comorbidities (hEDS, MCAS, SFN) |
| treatment_options | VARCHAR(512) | Treatment approaches |
| autoimmune_associations | VARCHAR(256) | Known autoimmune associations |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 13: `genomic_evidence`** -- Shared genomic variant data (read-only)

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Variant identifier |
| text_summary | VARCHAR(2048) | Variant description |
| chrom | VARCHAR(2) | Chromosome |
| pos | INT64 | Position |
| ref | VARCHAR(512) | Reference allele |
| alt | VARCHAR(512) | Alternate allele |
| gene | VARCHAR(16) | Gene symbol |
| consequence | VARCHAR(32) | Variant consequence |
| clinical_significance | VARCHAR(32) | ClinVar classification |
| disease_associations | VARCHAR(512) | Associated diseases |
| am_pathogenicity | FLOAT | AlphaMissense score |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

**Collection 14: `autoimmune_biomarker_trends`** -- Longitudinal biomarker trend data

| Field | Type | Description |
|-------|------|-------------|
| id | VARCHAR(64) | Unique identifier |
| patient_id | VARCHAR(32) | Patient identifier |
| biomarker | VARCHAR(32) | Biomarker name |
| trend_type | VARCHAR(16) | rising, falling, stable, fluctuating |
| values_json | VARCHAR(1024) | JSON array of {date, value} pairs |
| trend_slope | FLOAT | Calculated trend slope |
| clinical_significance | VARCHAR(256) | Clinical interpretation |
| associated_disease | VARCHAR(64) | Disease context |
| alert_level | VARCHAR(16) | normal, watch, warning, critical |
| text_summary | VARCHAR(512) | Trend narrative |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 vector |

### 7.3 Collection Search Weights

| Collection | Weight | Rationale |
|-----------|--------|-----------|
| Clinical Documents | 0.18 | Patient-specific primary evidence |
| Patient Labs | 0.16 | Quantitative diagnostic data |
| Patient Timeline | 0.10 | Longitudinal pattern context |
| Literature | 0.10 | Published evidence base |
| Autoantibodies | 0.08 | Diagnostic specificity data |
| HLA Associations | 0.07 | Genetic susceptibility evidence |
| Disease Activity | 0.06 | Scoring system reference |
| Guidelines | 0.06 | Clinical practice standards |
| Overlap Syndromes | 0.05 | Multi-disease pattern recognition |
| Biologics | 0.04 | Treatment reference |
| Dysautonomia | 0.04 | POTS/autonomic dysfunction reference |
| Genomic Evidence | 0.03 | Variant-level genetic data |
| Biomarker Trends | 0.02 | Longitudinal biomarker patterns |
| Trials | 0.01 | Clinical trial context |

### 7.4 Estimated Vector Counts

| Collection | Estimated Vectors | Source |
|-----------|------------------|--------|
| Clinical Documents | 5,000-50,000 per patient | Patient clinical records |
| Patient Labs | 1,000-10,000 per patient | Extracted laboratory results |
| Patient Timeline | 100-1,000 per patient | Longitudinal events |
| Literature | 8,000 | PubMed autoimmune corpus |
| Trials | 2,500 | ClinicalTrials.gov |
| Autoantibodies | 200 | Curated reference data |
| HLA Associations | 300 | Curated from GWAS/literature |
| Disease Activity | 150 | Scoring system reference |
| Biologics | 100 | Drug reference data |
| Overlap Syndromes | 80 | Curated reference data |
| Guidelines | 500 | ACR/EULAR/AGA guidelines |
| Dysautonomia | 150 | Curated reference data |
| Genomic Evidence | 3,561,170 | HCLS AI Factory Stage 1+2 (shared, read-only) |
| Biomarker Trends | 100-5,000 per patient | Computed from Patient Labs |
| **Reference total** | **~3,573,150** | |
| **Per-patient addition** | **~6,100-66,000** | |

---

## 8. Clinical Workflows

### 8.1 Workflow Architecture

Each workflow follows the pattern: **Trigger -> Ingest -> Extract -> Correlate -> Score -> Report**. Workflows can be triggered manually by clinician query, automatically upon document ingestion, or on a scheduled basis for longitudinal monitoring.

### 8.2 Eight Reference Workflows

#### Workflow 1: Diagnostic Odyssey Accelerator

**Purpose:** Ingest a patient's complete clinical record set and generate a prioritized differential diagnosis for suspected autoimmune disease.

**Trigger:** Clinician uploads patient's clinical document package or connects to EHR export.

**Process:**
1. Ingest all clinical documents through the document ingestion pipeline (Section 6)
2. Extract all laboratory values, autoantibody results, imaging findings, symptoms, and diagnoses
3. Construct patient timeline with all clinical events
4. Query all 14 collections simultaneously to identify patterns consistent with autoimmune diseases
5. Score each candidate diagnosis using classification criteria:
   - SLE: 2019 ACR/EULAR criteria (score >= 10 with ANA entry criterion)
   - RA: 2010 ACR/EULAR criteria (score >= 6)
   - Sjogren's: 2016 ACR/EULAR criteria (score >= 4)
   - Systemic sclerosis: 2013 ACR/EULAR criteria (score >= 9)
   - AS: Modified New York criteria or ASAS criteria
6. Cross-reference with genomic evidence (HLA alleles, autoimmune risk variants)
7. Generate diagnostic hypothesis report with confidence scores, supporting evidence, and recommended next steps

**Output:** Ranked differential diagnosis with evidence citations, classification criteria scores, recommended confirmatory tests, and urgency assessment.

**Example query:** "This 34-year-old woman has been evaluated by 5 specialists over 3 years. What autoimmune diagnoses are supported by her clinical record?"

#### Workflow 2: Lupus Nephritis Surveillance

**Purpose:** Monitor SLE patients for early signs of renal involvement through longitudinal biomarker tracking.

**Trigger:** Automated quarterly review or upon ingestion of new laboratory results.

**Process:**
1. Track longitudinal biomarkers: anti-dsDNA titers, complement C3/C4, urinalysis (protein, casts), serum creatinine, urine protein-to-creatinine ratio
2. Detect flare-predictive patterns: rising anti-dsDNA + falling complement precedes clinical flare by 4-12 weeks
3. Calculate SLEDAI-2K renal domain score from available data
4. Cross-reference with ISN/RPS 2003 lupus nephritis classification if biopsy data available
5. Assess current immunosuppressive regimen against ACR 2024 lupus nephritis guidelines
6. Generate surveillance report with trend visualizations and action items

**Output:** Lupus nephritis risk assessment, biomarker trend analysis, SLEDAI-2K score, treatment adequacy evaluation, and recommended monitoring schedule.

**Example query:** "Has this lupus patient's complement been trending down? What is her current nephritis risk?"

#### Workflow 3: POTS/Dysautonomia Evaluation

**Purpose:** Identify POTS and associated conditions (hEDS, MCAS, small fiber neuropathy) from clinical records, reducing the typical 5-7 year diagnostic delay.

**Trigger:** Clinician query or detection of orthostatic vital sign patterns in ingested records.

**Process:**
1. Extract orthostatic vital signs from clinical documents (supine vs. standing HR and BP)
2. Identify POTS criteria: sustained HR increase >= 30 bpm (>= 40 bpm if 12-19 years old) within 10 minutes of standing, without orthostatic hypotension
3. Screen for comorbid conditions:
   - Ehlers-Danlos syndrome: joint hypermobility scores, skin findings, family history
   - Mast cell activation syndrome: tryptase levels, prostaglandin D2, histamine metabolites
   - Small fiber neuropathy: IENFD biopsy results, QSART data, sudomotor testing
4. Cross-reference with autoimmune etiologies: Sjogren's-associated autonomic neuropathy, autoimmune autonomic ganglionopathy (ganglionic AChR antibodies)
5. Check genomic evidence for relevant variants (TPSAB1 for hereditary alpha-tryptasemia, COL5A1/COL3A1 for EDS)
6. Generate comprehensive dysautonomia evaluation report

**Output:** POTS diagnostic assessment, comorbidity screening results, autoimmune etiology evaluation, genetic risk factors, and management recommendations.

**Example query:** "This patient has had tachycardia, fatigue, and syncope for 2 years. Do her records support a POTS diagnosis? Are there autoimmune associations?"

#### Workflow 4: Inflammatory Arthritis Differentiation

**Purpose:** Differentiate between rheumatoid arthritis, psoriatic arthritis, reactive arthritis, crystal arthropathies, and lupus arthritis using clinical document analysis.

**Trigger:** Clinician query for a patient presenting with inflammatory joint symptoms.

**Process:**
1. Extract joint examination findings from clinical documents (tender joints, swollen joints, distribution pattern)
2. Retrieve autoantibody results: RF, anti-CCP, ANA, HLA-B27
3. Extract imaging findings: erosions, joint space narrowing, enthesitis, dactylitis, sacroiliitis
4. Calculate DAS28-CRP and CDAI if components available
5. Apply classification criteria for RA (2010 ACR/EULAR), PsA (CASPAR), SpA (ASAS), and gout (2015 ACR/EULAR)
6. Cross-reference with HLA associations: B27 for SpA, DRB1 shared epitope for RA, C06:02 for psoriasis
7. Generate differential diagnosis with evidence

**Output:** Inflammatory arthritis differential diagnosis, classification criteria scores, HLA risk profile, disease activity assessment, and treatment pathway recommendations.

#### Workflow 5: Overlap Syndrome Detection

**Purpose:** Identify multi-disease autoimmune overlap syndromes that single-disease evaluation would miss.

**Trigger:** Detection of autoantibody patterns or symptom constellations spanning multiple autoimmune categories.

**Process:**
1. Analyze complete autoantibody profile for overlap patterns:
   - Anti-U1-RNP at high titer: mixed connective tissue disease (MCTD)
   - ANA + anti-dsDNA + RF + anti-CCP: rhupus (RA + SLE overlap)
   - Anti-SSA + anti-centromere: Sjogren's + limited systemic sclerosis overlap
   - Anti-Jo-1 + anti-SSA: antisynthetase syndrome with Sjogren's features
2. Screen for polyautoimmune syndromes (APS-1, APS-2, IPEX)
3. Cross-reference with genomic evidence for shared susceptibility loci (STAT4, IRF5, PTPN22, TNFAIP3)
4. Evaluate for immune dysregulation syndromes if multiple autoimmune conditions present
5. Generate overlap syndrome report with component disease assessments

**Output:** Overlap syndrome assessment, polyautoimmunity risk evaluation, shared genetic susceptibility analysis, and integrated management recommendations.

#### Workflow 6: Biologic Therapy Optimization

**Purpose:** Recommend biologic therapies based on disease profile, prior treatment responses, and pharmacogenomic considerations.

**Trigger:** Clinician query about treatment options or detection of inadequate disease control.

**Process:**
1. Assess current disease activity from most recent clinical data
2. Review treatment history: prior biologics, duration, response, reason for discontinuation
3. Apply 8-drug biologic database with indication matching:
   - TNF inhibitors: adalimumab, etanercept, infliximab, certolizumab, golimumab
   - IL-6 inhibitors: tocilizumab, sarilumab
   - IL-17 inhibitors: secukinumab, ixekizumab
   - IL-12/23 inhibitors: ustekinumab
   - IL-23 inhibitors: guselkumab, risankizumab
   - B-cell depleters: rituximab, obinutuzumab
   - BLyS inhibitors: belimumab
   - JAK inhibitors: tofacitinib, baricitinib, upadacitinib
   - T-cell co-stimulation modulators: abatacept
4. Apply pharmacogenomic filters:
   - HLA-DRB1*03:01: increased risk of anti-drug antibodies with adalimumab
   - FCGR3A V158F: affects rituximab ADCC efficacy
   - CYP3A4/CYP2C19: affects tofacitinib metabolism
   - IL6R Asp358Ala (rs2228145): affects tocilizumab response
   - HLA-C*06:02: predicts better PASI response to secukinumab
5. Check contraindications against patient profile (TB status, hepatitis B, heart failure, IBD)
6. Generate personalized treatment recommendation with evidence

**Output:** Ranked biologic therapy recommendations, pharmacogenomic considerations, contraindication alerts, monitoring requirements, and switching rationale from current therapy.

#### Workflow 7: Flare Prediction and Prevention

**Purpose:** Predict autoimmune disease flares by detecting early biomarker patterns in longitudinal laboratory data.

**Trigger:** Automated analysis upon ingestion of new laboratory results, or scheduled weekly/monthly review.

**Process:**
1. Extract longitudinal biomarker data from patient labs collection
2. Calculate trend analysis for disease-specific biomarkers:
   - RA: CRP trend, ESR trend, IL-6, MMP-3, 14-3-3eta
   - SLE: anti-dsDNA titer trend, complement C3/C4 trends, lymphocyte count, proteinuria
   - IBD: fecal calprotectin trend, CRP, lactoferrin, albumin
3. Apply flare prediction algorithm:
   - Base risk: 0.3 (30% baseline)
   - Each elevated inflammatory marker: +0.15
   - Each falling protective marker (complement, albumin): +0.15
   - Stable markers: protective factor
4. Classify risk: Low (<0.4), Moderate (0.4-0.6), High (0.6-0.8), Imminent (>0.8)
5. Generate monitoring recommendations and preemptive intervention suggestions
6. Alert clinician if risk crosses threshold

**Output:** Flare risk assessment with contributing factors, protective factors, recommended biomarker monitoring schedule, and intervention recommendations.

#### Workflow 8: Genomic-Autoimmune Risk Profiling

**Purpose:** Analyze a patient's genomic data for autoimmune disease risk variants, combining HLA typing with non-HLA susceptibility loci.

**Trigger:** Availability of genomic data (from HCLS AI Factory Stage 1 or external genetic testing reports).

**Process:**
1. Query `genomic_evidence` collection for autoimmune-associated variants:
   - HLA alleles: B*27 (AS), DRB1*04 (RA), DRB1*03 (SLE, Sjogren's, T1D), DRB1*15 (MS), DQB1*02 (celiac), C*06 (psoriasis), B*51 (Behcet's), B*08 (MG)
   - Non-HLA risk genes: PTPN22 R620W (multiple autoimmune), STAT4 rs7574865 (SLE, RA), IRF5 rs2004640 (SLE), TNFAIP3 rs2230926 (SLE, RA), IL23R rs11209026 (IBD, psoriasis, AS), CTLA4 +49 A/G (autoimmune thyroid, T1D), TPSAB1 copy number (hereditary alpha-tryptasemia / MCAS)
2. Calculate polygenic risk scores for each autoimmune condition
3. Cross-reference genetic risk with clinical findings from patient documents
4. Identify pharmacogenomically actionable variants
5. Generate comprehensive genomic-autoimmune risk report

**Output:** Autoimmune genetic risk profile, HLA-disease associations with odds ratios, non-HLA risk variant assessment, pharmacogenomic actionable findings, and recommended genetic counseling topics.

---

## 9. Cross-Modal Integration and Genomic Correlation

### 9.1 The Genomic-Autoimmune Bridge

The shared `genomic_evidence` collection (3,561,170 variants) enables a transformative capability: automatic correlation between a patient's genomic profile and their autoimmune clinical presentation. This bridge operates through several mechanisms:

**HLA Typing from Genomic Data:**
When a patient's genome is processed through Stage 1 (Parabricks), HLA alleles can be extracted using HLA typing tools (OptiType for Class I, HLA-HD for Class I+II). These alleles are then automatically cross-referenced against the `autoimmune_hla_associations` collection to generate disease susceptibility profiles.

**Non-HLA Autoimmune Risk Variants:**
The system queries the `genomic_evidence` collection for known autoimmune susceptibility loci:

| Gene | Variant | Associated Diseases | OR Range |
|------|---------|-------------------|----------|
| PTPN22 | R620W (rs2476601) | RA, SLE, T1D, Hashimoto's, Graves' | 1.5-2.0 |
| STAT4 | rs7574865 | SLE, RA, Sjogren's | 1.3-1.6 |
| IRF5 | rs2004640 | SLE | 1.4-1.8 |
| TNFAIP3 | rs2230926 | SLE, RA | 1.7-2.3 |
| IL23R | rs11209026 | Crohn's, UC, AS, psoriasis | 0.4 (protective) |
| CTLA4 | rs231775 | T1D, autoimmune thyroid, RA | 1.2-1.5 |
| IL2RA | rs2104286 | MS, T1D | 1.1-1.3 |
| BANK1 | rs10516487 | SLE | 1.3-1.4 |
| BLK | rs13277113 | SLE | 1.3-1.4 |
| TPSAB1 | CNV (duplication) | Hereditary alpha-tryptasemia / MCAS | 4-6% prevalence |

### 9.2 Autoimmune Trigger Conditions

The system defines 12 genomic-autoimmune trigger conditions that automatically activate cross-collection queries:

1. **HLA-B*27 detected** -> Query AS, reactive arthritis, anterior uveitis, IBD-associated spondyloarthritis collections
2. **HLA-DRB1*04 shared epitope detected** -> Query RA literature, anti-CCP interpretation, ACPA-positive RA management
3. **HLA-DRB1*03 + DQB1*02 detected** -> Query celiac, T1D, SLE, Sjogren's, autoimmune thyroid collections
4. **HLA-DRB1*15 detected** -> Query MS literature, NMO differential, imaging (brain MRI white matter lesions)
5. **PTPN22 R620W detected** -> Query polyautoimmunity risk, T-cell signaling pathway, multiple autoimmune screening
6. **HLA-C*06:02 detected** -> Query psoriasis, PsA, IL-17 inhibitor response prediction
7. **IL23R protective variant detected** -> Note protective factor for IBD, psoriasis, AS in diagnostic scoring
8. **STAT4 risk variant detected** -> Increase SLE, RA susceptibility weighting in diagnostic algorithm
9. **Anti-dsDNA + low complement + TNFAIP3 variant** -> High-priority lupus nephritis surveillance trigger
10. **HLA-B*51 detected** -> Query Behcet's disease, evaluate oral/genital ulcers, pathergy testing
11. **TPSAB1 duplication detected** -> Query MCAS, hereditary alpha-tryptasemia, POTS/hEDS/MCAS triad
12. **Multiple autoimmune risk variants detected** -> Activate polyautoimmunity screening workflow

### 9.3 Integration with Other Agents

The Precision Autoimmune Agent integrates with sibling agents in the HCLS AI Factory:

| Agent | Integration | Example |
|-------|------------|---------|
| Precision Biomarker | Inflammation biomarker interpretation | CRP, ESR, ferritin, calprotectin trends shared with autoimmune flare prediction |
| Precision Oncology | Autoimmune paraneoplastic screening | Anti-NMDA receptor encephalitis -> ovarian teratoma screening |
| Imaging Intelligence | Joint/organ imaging correlation | Joint erosions on hand MRI -> RA disease activity assessment |
| CAR-T Intelligence | Autoimmune complications of CAR-T | CRS, autoimmune cytopenias post-CAR-T therapy |
| Cardiology (future) | Cardiac autoimmune manifestations | Lupus pericarditis, myocarditis, autoimmune POTS |
| Neurology (future) | Neurological autoimmune manifestations | MS, NMO, autoimmune encephalitis, autoimmune neuropathy |

---

## 10. NIM Integration Strategy

### 10.1 Shared NIM Services

The Precision Autoimmune Agent leverages NVIDIA NIM microservices already deployed for the HCLS AI Factory:

| NIM Service | Port | Autoimmune Application |
|------------|------|----------------------|
| Llama-3 8B | 8520 | Local LLM fallback for evidence synthesis when Claude API unavailable |
| VISTA-3D | 8530 (shared with Imaging) | Joint imaging segmentation for arthritis assessment |
| VILA-M3 | 8532 (shared with Imaging) | Clinical document image understanding (scanned documents) |

### 10.2 Future NIM Extensions

| NIM Service | Application | Status |
|------------|-------------|--------|
| BioNeMo ESMFold | Autoantibody structure prediction for epitope analysis | Planned |
| NeMo Guardrails | Safety guardrails for clinical recommendation generation | Planned |
| NVIDIA FLARE | Federated learning for multi-institutional autoimmune pattern discovery | Research |

---

## 11. Knowledge Graph Design

### 11.1 Graph Structure

The Precision Autoimmune Agent employs a 7-dictionary knowledge graph containing structured clinical data that complements vector retrieval:

| Dictionary | Entries | Content |
|-----------|---------|---------|
| `HLA_DISEASE_ASSOCIATIONS` | 50+ | HLA alleles mapped to autoimmune diseases with odds ratios, PMIDs, and mechanism notes |
| `AUTOANTIBODY_DISEASE_MAP` | 14 antibodies x 1-4 diseases | Autoantibody-disease associations with sensitivity, specificity, and staining patterns |
| `DISEASE_ACTIVITY_THRESHOLDS` | 5 scoring systems | DAS28-CRP, DAS28-ESR, SLEDAI-2K, CDAI, BASDAI with component definitions and level thresholds |
| `BIOLOGIC_THERAPIES` | 8 drugs (expandable to 20+) | Drug class, mechanism, indications, PGx considerations, contraindications, monitoring |
| `FLARE_BIOMARKER_PATTERNS` | 3 diseases (expandable to 13) | Early warning biomarkers, threshold patterns, protective signals for RA, SLE, IBD |
| `OVERLAP_SYNDROMES` | 12 syndromes | Component diseases, diagnostic criteria, key autoantibodies, prevalence, management |
| `DYSAUTONOMIA_CONDITIONS` | 10 conditions | POTS, NCS, MSA, AAG, hEDS, MCAS, SFN, PAF, baroreflex failure, familial dysautonomia |

### 11.2 Example Knowledge Graph Entries

**HLA-Disease Association (HLA-B*27:05):**
```json
{
  "allele": "HLA-B*27:05",
  "associations": [
    {
      "disease": "ankylosing_spondylitis",
      "odds_ratio": 87.4,
      "pmid": "25603694",
      "note": "Strongest known HLA-disease association. Arthritogenic peptide hypothesis: B27 presents self-peptides to autoreactive CD8+ T-cells. Also: B27 misfolding triggers UPR and IL-23 production.",
      "population": "European",
      "mechanism": "Arthritogenic peptide presentation + misfolding/UPR + IL-23 axis"
    },
    {
      "disease": "reactive_arthritis",
      "odds_ratio": 20.0,
      "pmid": "25603694"
    },
    {
      "disease": "anterior_uveitis",
      "odds_ratio": 10.5,
      "pmid": "25603694"
    },
    {
      "disease": "ibd_spondyloarthritis",
      "odds_ratio": 8.0,
      "pmid": "25603694"
    }
  ]
}
```

**Autoantibody Map (anti-dsDNA):**
```json
{
  "antibody": "anti-dsDNA",
  "associations": [
    {
      "disease": "systemic_lupus_erythematosus",
      "sensitivity": 0.70,
      "specificity": 0.95,
      "note": "Titers correlate with disease activity, especially lupus nephritis. Rising titers precede clinical flare by 4-12 weeks. Part of 2019 ACR/EULAR SLE criteria (6 points).",
      "assay_methods": ["Farr assay (gold standard)", "ELISA", "CLIFT (Crithidia luciliae)"],
      "monitoring_frequency": "Every 3-6 months in active SLE; with complement levels"
    }
  ]
}
```

**POTS/Dysautonomia Entry:**
```json
{
  "condition": "postural_orthostatic_tachycardia_syndrome",
  "abbreviation": "POTS",
  "diagnostic_criteria": {
    "heart_rate_increase": ">= 30 bpm within 10 min of standing (adults); >= 40 bpm (12-19 years)",
    "absence_of": "Orthostatic hypotension (BP drop > 20/10 mmHg)",
    "duration": "Symptoms present >= 6 months",
    "exclusion": "No other cause of tachycardia (anemia, hyperthyroidism, deconditioning)"
  },
  "subtypes": ["Neuropathic POTS", "Hyperadrenergic POTS", "Hypovolemic POTS", "Autoimmune POTS"],
  "comorbidities": {
    "ehlers_danlos_heds": {"prevalence": "up to 80%", "mechanism": "Connective tissue laxity -> venous pooling"},
    "mcas": {"prevalence": "up to 65%", "mechanism": "Mast cell mediator release -> vasodilation"},
    "small_fiber_neuropathy": {"prevalence": "up to 50%", "mechanism": "Sudomotor/vasomotor nerve loss"},
    "sjogrens_syndrome": {"prevalence": "15-25%", "mechanism": "Autoimmune autonomic neuropathy"},
    "autoimmune_ganglionopathy": {"antibody": "ganglionic AChR", "prevalence": "10-15%"}
  },
  "autoimmune_associations": [
    "Ganglionic AChR antibodies (autoimmune autonomic ganglionopathy)",
    "Anti-adrenergic receptor antibodies",
    "Anti-muscarinic receptor antibodies",
    "Sjogren's-associated autonomic neuropathy",
    "Post-viral autoimmune autonomic dysfunction"
  ]
}
```

**Overlap Syndrome Entry (Mixed Connective Tissue Disease):**
```json
{
  "syndrome": "mixed_connective_tissue_disease",
  "abbreviation": "MCTD",
  "component_diseases": ["SLE", "Systemic sclerosis", "Polymyositis"],
  "key_autoantibody": "anti-U1-RNP (high titer, required for diagnosis)",
  "diagnostic_criteria": "Alarcon-Segovia criteria: anti-U1-RNP >= 1:1600 + 3 of 5 clinical criteria (edema of hands, synovitis, myositis, Raynaud's, acrosclerosis)",
  "prevalence": "1.9-3.8 per 100,000",
  "management": "Treat dominant clinical feature; corticosteroids for myositis/serositis; immunosuppressants for organ involvement",
  "prognosis": "May evolve into definite SLE, SSc, or PM over time (differentiated MCTD)"
}
```

---

## 12. Query Expansion and Retrieval Strategy

### 12.1 Autoimmune-Specific Query Expansion Maps

The system implements 18 domain-specific query expansion maps:

| # | Category | Keywords | Expanded Terms | Example |
|---|----------|----------|---------------|---------|
| 1 | Autoimmune diseases | 30 | 280 | "lupus" -> SLE, lupus nephritis, anti-dsDNA, complement C3/C4, hydroxychloroquine, belimumab, ACR/EULAR criteria |
| 2 | Autoantibodies | 20 | 180 | "ANA" -> antinuclear antibody, ANA pattern, homogeneous, speckled, nucleolar, centromere, IIF, HEp-2 |
| 3 | HLA alleles | 15 | 120 | "B27" -> HLA-B*27, ankylosing spondylitis, spondyloarthritis, sacroiliitis, uveitis, reactive arthritis |
| 4 | Inflammatory markers | 12 | 90 | "CRP" -> C-reactive protein, inflammation, acute phase, IL-6 driven, liver synthesis |
| 5 | Biologics | 20 | 160 | "rituximab" -> anti-CD20, B-cell depletion, ADCC, IV infusion, PML risk, hepatitis B screening |
| 6 | Disease activity | 8 | 60 | "DAS28" -> disease activity score, 28 joints, CRP, ESR, tender, swollen, patient global, remission |
| 7 | Dysautonomia | 15 | 140 | "POTS" -> postural orthostatic tachycardia, tilt table test, standing HR, orthostatic intolerance, dysautonomia |
| 8 | Overlap syndromes | 10 | 80 | "MCTD" -> mixed connective tissue disease, anti-U1-RNP, overlap syndrome, Raynaud's, edema of hands |
| 9 | Immunology mechanisms | 18 | 150 | "Th17" -> T-helper 17, IL-17, RORgammaT, autoimmune, mucosal immunity, psoriasis, SpA, IBD |
| 10 | Genomics | 15 | 120 | "PTPN22" -> protein tyrosine phosphatase, R620W, rs2476601, T-cell signaling, autoimmune risk variant |
| 11 | Laboratory panels | 12 | 90 | "ENA panel" -> extractable nuclear antigens, SSA, SSB, Sm, RNP, Scl-70, Jo-1 |
| 12 | Imaging | 10 | 70 | "joint erosion" -> marginal erosion, bone erosion, MRI, X-ray, Sharp score, modified Sharp, progression |
| 13 | Flare patterns | 8 | 60 | "flare" -> disease flare, exacerbation, relapse, breakthrough, loss of response, secondary failure |
| 14 | Treatments | 18 | 140 | "methotrexate" -> MTX, csDMARD, folic acid, hepatotoxicity, pneumonitis, weekly dosing, anchor drug |
| 15 | Pregnancy | 8 | 60 | "pregnancy lupus" -> neonatal lupus, anti-SSA, congenital heart block, hydroxychloroquine continuation |
| 16 | Pediatric | 10 | 80 | "JIA" -> juvenile idiopathic arthritis, oligoarticular, polyarticular, systemic, enthesitis-related, RF positive |
| 17 | Infections | 10 | 70 | "TB screening" -> tuberculosis, QuantiFERON, T-SPOT, latent TB, biologic contraindication, isoniazid |
| 18 | Comorbidities | 12 | 90 | "cardiovascular risk" -> accelerated atherosclerosis, lupus vasculitis, anti-phospholipid, thrombosis |
| | **Total** | **251** | **2,040** | |

### 12.2 Comparative Analysis Detection

The system auto-detects comparative queries and routes them to dual-retrieval:

**Supported comparison types:**
- Disease vs. disease: "Compare RA vs lupus arthritis," "POTS vs neurocardiogenic syncope"
- Drug vs. drug: "Adalimumab vs rituximab for RA," "Tofacitinib vs baricitinib"
- Antibody vs. antibody: "RF vs anti-CCP for RA diagnosis"
- Scoring systems: "DAS28 vs CDAI for RA monitoring"
- HLA alleles: "B*27:05 vs B*27:02 for AS risk"

---

## 13. API and UI Design

### 13.1 FastAPI Endpoints (Port 8530)

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/query` | Multi-collection RAG query with evidence synthesis |
| POST | `/compare` | Comparative analysis (X vs Y) |
| POST | `/documents/ingest` | Ingest patient clinical documents |
| POST | `/documents/batch` | Batch ingest multiple documents |
| GET | `/documents/{patient_id}/status` | Document ingestion status |
| GET | `/patient/{patient_id}/timeline` | Patient clinical timeline |
| GET | `/patient/{patient_id}/labs` | Patient laboratory trends |
| POST | `/patient/{patient_id}/diagnostic` | Run diagnostic odyssey workflow |
| POST | `/patient/{patient_id}/flare-risk` | Calculate flare risk assessment |
| POST | `/patient/{patient_id}/genomic-risk` | Genomic-autoimmune risk profile |
| POST | `/patient/{patient_id}/biologic-rec` | Biologic therapy recommendation |
| POST | `/reports/generate` | Generate clinical report (Markdown, JSON, PDF) |
| POST | `/reports/fhir` | Generate FHIR R4 DiagnosticReport |
| GET | `/collections/stats` | Collection vector counts |
| GET | `/health` | Health check |
| GET | `/metrics` | Prometheus metrics |

### 13.2 Streamlit UI (Port 8531) -- 10 Tabs

| Tab | Name | Purpose |
|-----|------|---------|
| 1 | Evidence Explorer | Multi-collection RAG search across autoimmune knowledge |
| 2 | Document Ingest | Upload and process patient clinical documents |
| 3 | Diagnostic Workup | Run diagnostic odyssey workflow with classification criteria scoring |
| 4 | Patient Timeline | Interactive timeline visualization of clinical events |
| 5 | Lab Trends | Longitudinal biomarker trend analysis with flare prediction |
| 6 | Genomic Risk | Autoimmune genetic risk profiling (HLA + non-HLA variants) |
| 7 | Biologic Advisor | Pharmacogenomic-guided therapy recommendations |
| 8 | Overlap Detector | Multi-disease overlap syndrome assessment |
| 9 | Reports & Export | Generate PDF/FHIR reports with NVIDIA branding |
| 10 | Benchmarks | Diagnostic accuracy validation against classification criteria |

### 13.3 Demo Cases

| # | Case | Key Findings | Expected Diagnosis |
|---|------|-------------|-------------------|
| 1 | 34F, 3-year diagnostic odyssey | Malar rash, arthritis, proteinuria, ANA 1:640 homogeneous, anti-dsDNA+, low C3/C4 | SLE with lupus nephritis |
| 2 | 28F, chronic fatigue and syncope | HR increase 42 bpm on standing, joint hypermobility score 7/9, elevated tryptase | POTS / hEDS / MCAS triad |
| 3 | 45M, inflammatory back pain | HLA-B*27+, bilateral sacroiliitis on MRI, CRP 24 mg/L, morning stiffness >30 min | Ankylosing spondylitis |
| 4 | 52F, dry eyes/mouth + fatigue | ANA 1:320 speckled, anti-SSA+, anti-SSB+, Schirmer test <5mm, lip biopsy focus score 3 | Sjogren's syndrome |
| 5 | 38F, multiple autoimmune features | Anti-U1-RNP 1:5120, Raynaud's, puffy fingers, myositis, arthritis | Mixed connective tissue disease |

---

## 14. Clinical Decision Support Engines

### 14.1 Validated Disease Activity and Diagnostic Scores

| Score | Disease | Components | Thresholds |
|-------|---------|-----------|-----------|
| DAS28-CRP | Rheumatoid arthritis | Tender joints (28), swollen joints (28), CRP, patient global VAS | Remission <2.6, Low <3.2, Moderate <5.1, High >=5.1 |
| DAS28-ESR | Rheumatoid arthritis | Tender joints (28), swollen joints (28), ESR, patient global VAS | Remission <2.6, Low <3.2, Moderate <5.1, High >=5.1 |
| SLEDAI-2K | SLE | 16 weighted items (seizure, psychosis, vasculitis, arthritis, etc.) | Inactive 0, Mild 1-4, Moderate 5-11, High 12+, Very high 20+ |
| CDAI | Rheumatoid arthritis | Tender joints (28), swollen joints (28), patient global, evaluator global | Remission <=2.8, Low <=10, Moderate <=22, High >22 |
| BASDAI | Ankylosing spondylitis | 6 questions (fatigue, spinal pain, joint pain, enthesitis, stiffness) | Inactive <2, Low 2-3, Moderate 3-4, Active >=4 |
| ACR/EULAR SLE 2019 | SLE (diagnosis) | Entry: ANA+; 7 domains, 22 criteria with weights | Score >=10: classify as SLE |
| ACR/EULAR RA 2010 | RA (diagnosis) | Joint involvement, serology (RF/anti-CCP), acute phase, duration | Score >=6: classify as RA |
| Beighton score | Joint hypermobility | 9-point scale assessing bilateral flexibility | >=5/9: generalized joint hypermobility |
| COMPASS-31 | Dysautonomia | 31 items across 6 autonomic domains | Higher scores = more severe autonomic dysfunction |
| Sheldon POTS criteria | POTS | Standing HR increase, absence of orthostatic hypotension, duration | Meet all criteria = POTS diagnosis |

### 14.2 Classification Criteria Engine

The system implements automated scoring for major classification criteria:

- **2019 ACR/EULAR SLE criteria:** Entry criterion (ANA >= 1:80) + additive weighted criteria across 7 clinical domains and 4 immunology domains. Score >= 10 classifies as SLE. Each criterion is scored only if not better explained by another diagnosis.
- **2010 ACR/EULAR RA criteria:** Four domains (joint involvement 0-5, serology 0-3, acute phase reactants 0-1, symptom duration 0-1). Score >= 6 classifies as RA.
- **2016 ACR/EULAR Sjogren's criteria:** Weighted items including labial gland biopsy focus score >= 1 (3 points), anti-SSA+ (3 points), ocular staining score >= 5 (1 point), Schirmer test <= 5mm (1 point), unstimulated salivary flow <= 0.1 mL/min (1 point). Score >= 4 classifies as Sjogren's.
- **CASPAR PsA criteria:** Inflammatory musculoskeletal disease + score >= 3 from: current psoriasis (2), personal/family history psoriasis (1), nail dystrophy (1), negative RF (1), dactylitis (1), juxta-articular bone formation (1).

---

## 15. Reporting and Interoperability

### 15.1 Export Formats

| Format | Use Case | Implementation |
|--------|----------|---------------|
| Markdown | In-app display, clinician review | Template-based rendering with evidence tables |
| JSON | API integration, downstream analysis | Pydantic `.model_dump()` serialization |
| PDF | Formal clinical reports, sharing | ReportLab with NVIDIA branding (#76B900 green) |
| FHIR R4 | EHR integration, interoperability | DiagnosticReport resource with coded observations |

### 15.2 FHIR R4 Autoimmune Coding

| FHIR Resource | Coding System | Example |
|--------------|---------------|---------|
| DiagnosticReport | LOINC | 51967-8 (Genetic analysis summary) |
| Condition | SNOMED CT | 55464009 (Systemic lupus erythematosus) |
| Condition | SNOMED CT | 69896004 (Rheumatoid arthritis) |
| Condition | ICD-10-CM | M32.14 (Lupus nephritis) |
| Observation | LOINC | 33935-8 (Anti-dsDNA Ab, quantitative) |
| Observation | LOINC | 14585-3 (ANA by IIF) |
| MedicationRequest | RxNorm | 327361 (Adalimumab 40mg/0.8mL) |
| AllergyIntolerance | SNOMED CT | 294468003 (Rituximab adverse reaction) |
| Observation | LOINC | 30522-7 (C-reactive protein, high sensitivity) |

---

## 16. Product Requirements Document

### 16.1 Product Vision

**For** rheumatologists, immunologists, and primary care physicians **who** need to diagnose complex autoimmune diseases faster and more accurately, **the** Precision Autoimmune Agent **is a** clinical document intelligence system **that** ingests a patient's complete clinical record, identifies patterns across years of clinical data, cross-references with genomic and autoimmune reference knowledge, and generates prioritized diagnostic hypotheses with evidence citations. **Unlike** existing autoimmune diagnostic tools that require manual data entry and evaluate one disease at a time, **our product** ingests raw clinical documents, analyzes all supported autoimmune conditions simultaneously, and incorporates genomic correlation -- reducing the average diagnostic odyssey from 4+ years to weeks.

### 16.2 User Stories

**Epic 1: Clinical Document Ingestion**

| ID | Story | Priority |
|----|-------|----------|
| US-1.1 | As a rheumatologist, I want to upload a patient's PDF medical records so that the system can extract and index all relevant clinical data | P0 |
| US-1.2 | As a clinician, I want to import FHIR R4 bundles from the patient's EHR portal so that structured data is automatically ingested | P0 |
| US-1.3 | As a clinician, I want to see the status of document processing (in progress, completed, errors) so that I know when analysis is ready | P0 |
| US-1.4 | As a clinician, I want scanned documents to be OCR-processed so that handwritten or faxed records are included in the analysis | P1 |

**Epic 2: Diagnostic Intelligence**

| ID | Story | Priority |
|----|-------|----------|
| US-2.1 | As a rheumatologist, I want the system to score my patient against ACR/EULAR classification criteria for SLE, RA, and Sjogren's so that I can see which diagnoses are supported by the data | P0 |
| US-2.2 | As a clinician, I want a ranked differential diagnosis with confidence scores so that I can prioritize my diagnostic workup | P0 |
| US-2.3 | As a PCP, I want the system to flag autoimmune red flags in my patient's records so that I can make appropriate specialist referrals | P0 |
| US-2.4 | As a rheumatologist, I want the system to detect overlap syndromes (MCTD, rhupus) so that I don't miss multi-disease patterns | P1 |

**Epic 3: Genomic-Autoimmune Correlation**

| ID | Story | Priority |
|----|-------|----------|
| US-3.1 | As a clinician, I want to see HLA-disease associations with odds ratios when HLA typing is available so that I can assess genetic susceptibility | P0 |
| US-3.2 | As a geneticist, I want to query the patient's genomic data for known autoimmune risk variants (PTPN22, STAT4, IRF5) so that I can generate a genetic risk profile | P1 |
| US-3.3 | As a rheumatologist, I want pharmacogenomic data integrated with biologic therapy recommendations so that I can select the most effective treatment | P1 |

**Epic 4: Longitudinal Monitoring**

| ID | Story | Priority |
|----|-------|----------|
| US-4.1 | As a rheumatologist, I want to track disease activity scores (DAS28, SLEDAI) over time so that I can assess treatment response | P0 |
| US-4.2 | As a clinician, I want automated flare risk predictions based on biomarker trends so that I can intervene before clinical flares occur | P1 |
| US-4.3 | As a clinician, I want a visual patient timeline showing all clinical events, labs, and diagnoses so that I can see the complete clinical picture at a glance | P1 |

**Epic 5: POTS/Dysautonomia**

| ID | Story | Priority |
|----|-------|----------|
| US-5.1 | As a clinician, I want the system to identify POTS criteria from orthostatic vital signs in the clinical record so that POTS is not missed | P1 |
| US-5.2 | As a clinician, I want the system to screen for POTS comorbidities (hEDS, MCAS, SFN) so that the complete syndrome is identified | P1 |
| US-5.3 | As a clinician, I want autoimmune etiologies of POTS (Sjogren's, AAG) evaluated automatically so that treatable causes are identified | P2 |

**Epic 6: Reporting and Export**

| ID | Story | Priority |
|----|-------|----------|
| US-6.1 | As a clinician, I want to generate a PDF diagnostic report with evidence citations so that I can share findings with colleagues | P0 |
| US-6.2 | As a health system, I want FHIR R4 DiagnosticReport export so that findings can be integrated back into the EHR | P1 |
| US-6.3 | As a researcher, I want JSON export of all analysis results so that I can perform downstream statistical analysis | P2 |

### 16.3 Non-Functional Requirements

| Requirement | Target | Rationale |
|------------|--------|-----------|
| Document ingestion throughput | 100 documents/minute | Support batch upload of complete patient records |
| Entity extraction accuracy | >90% F1 for laboratory values | Clinical reliability requirement |
| Query response time (reference) | <30 seconds | Acceptable for clinical decision support |
| Query response time (patient) | <60 seconds (with 10K patient vectors) | Larger patient-specific search space |
| Concurrent patients | 50 per DGX Spark | Multi-clinician usage scenario |
| Data retention | Configurable per institution | HIPAA compliance |
| Uptime | 99.5% | Clinical workflow reliability |
| HIPAA compliance | Full | Required for patient data handling |

### 16.4 Prioritization Matrix

| Phase | Scope | Duration |
|-------|-------|----------|
| **Phase 1 (MVP)** | Reference collections (literature, trials, autoantibodies, HLA, guidelines) + basic RAG query + autoantibody interpretation + HLA analysis + disease activity scoring | 6 weeks |
| **Phase 2 (Document Intelligence)** | Clinical document ingestion pipeline + patient-specific collections + timeline construction + laboratory extraction + NLP entity extraction | 6 weeks |
| **Phase 3 (Advanced Analytics)** | Genomic correlation + flare prediction + overlap syndrome detection + biologic therapy optimization + POTS/dysautonomia evaluation + FHIR export + comparative analysis | 6 weeks |

---

## 17. Data Acquisition Strategy

### 17.1 Automated Ingest Pipelines

| Source | Method | Target Collection | Refresh Cadence |
|--------|--------|------------------|----------------|
| PubMed | E-utilities API (esearch + efetch) | autoimmune_literature | Weekly |
| ClinicalTrials.gov | V2 API | autoimmune_trials | Weekly |
| ACR/EULAR guidelines | Manual curation + PDF ingestion | autoimmune_guidelines | Quarterly |
| HLA-disease GWAS | Literature curation + GWAS Catalog | autoimmune_hla_associations | Monthly |
| Autoantibody reference | Expert curation from ACR/EULAR criteria | autoimmune_autoantibodies | Quarterly |
| Biologic drug database | DailyMed + FDA labels | autoimmune_biologics | Monthly |
| Dysautonomia reference | Expert curation from consensus statements | autoimmune_dysautonomia | Quarterly |

### 17.2 Patient Document Sources

| Source | Format | Ingestion Method |
|--------|--------|-----------------|
| Patient portal exports | PDF, CCDA | Upload via Streamlit UI or API |
| EHR FHIR endpoints | FHIR R4 Bundle | SMART on FHIR integration |
| Genetic testing portals | PDF reports | OCR + genetic report parser |
| Laboratory result files | HL7 ORU, CSV | Structured parser |
| Outside records (faxed/scanned) | Image PDF | Tesseract OCR + NLP |
| Patient-provided records | Various | Document classifier + appropriate parser |

### 17.3 PubMed Search Strategy

```
("autoimmune disease" OR "autoimmunity" OR "systemic lupus erythematosus" OR
 "rheumatoid arthritis" OR "multiple sclerosis" OR "inflammatory bowel disease" OR
 "psoriasis" OR "ankylosing spondylitis" OR "sjogren" OR "scleroderma" OR
 "myasthenia gravis" OR "celiac disease" OR "type 1 diabetes" OR "POTS" OR
 "dysautonomia" OR "postural orthostatic tachycardia" OR "Ehlers-Danlos" OR
 "mast cell activation" OR "autoantibody" OR "HLA association" OR
 "biologic therapy" OR "immunosuppressive") AND
("diagnosis" OR "biomarker" OR "classification criteria" OR "disease activity" OR
 "genomic" OR "genetic risk" OR "pharmacogenomics" OR "flare prediction" OR
 "overlap syndrome" OR "treatment response")
```

Estimated yield: 8,000-12,000 abstracts per refresh cycle.

---

## 18. Validation and Testing Strategy

### 18.1 Test Architecture

| Test Type | Scope | Target Count |
|-----------|-------|-------------|
| Unit tests | Knowledge graph, NLP extractors, scoring engines, models | 150+ |
| Integration tests | Collection operations, ingestion pipelines, API endpoints | 80+ |
| NLP accuracy tests | Entity extraction precision/recall against annotated clinical notes | 50+ |
| Clinical validation | Classification criteria scoring against expert-adjudicated cases | 30+ |
| End-to-end tests | Full workflow execution with synthetic patient records | 15+ |
| Performance tests | Query latency, ingestion throughput, concurrent access | 10+ |

### 18.2 NLP Validation Approach

The clinical document NLP pipeline will be validated against:

1. **i2b2/VA NLP challenge datasets:** De-identified clinical notes with annotated entities (medications, diagnoses, laboratory values) -- standard benchmark for clinical NLP
2. **MIMIC-III discharge summaries:** Large corpus of ICU discharge summaries with structured data for validation
3. **Synthetic autoimmune cases:** Custom-generated clinical documents with known entities for regression testing
4. **Expert review:** Board-certified rheumatologist review of entity extraction accuracy on 100 sampled documents

### 18.3 Diagnostic Accuracy Validation

Classification criteria scoring engines will be validated by:

1. Creating 50+ synthetic patient profiles with known ACR/EULAR classification criteria scores
2. Comparing system-calculated scores against expert-calculated scores
3. Measuring sensitivity, specificity, and concordance for each classification system
4. Target: >95% concordance with expert scoring for structured data, >85% for NLP-extracted data

---

## 19. Regulatory Considerations

### 19.1 FDA Classification

The Precision Autoimmune Agent is designed as a **Clinical Decision Support (CDS)** tool that meets the criteria for exemption from FDA device regulation under the 21st Century Cures Act, Section 3060(a):

1. **Not intended to replace clinical judgment:** The system provides diagnostic hypotheses and evidence for clinician review, not autonomous diagnosis
2. **Displays underlying evidence:** All recommendations include source citations, classification criteria scores, and confidence levels
3. **Clinician can independently verify:** Every evidence citation links to its primary source (PubMed, ClinicalTrials.gov, clinical documents)
4. **Intended for qualified professionals:** The system is designed for use by licensed clinicians, not patients

### 19.2 HIPAA Compliance

| Requirement | Implementation |
|------------|---------------|
| PHI at rest | LUKS full-disk encryption on DGX Spark NVMe |
| PHI in transit | TLS 1.3 for all API communications |
| Access control | JWT authentication with RBAC |
| Audit logging | All patient data access logged with timestamp, user, action |
| Minimum necessary | LLM queries contain anonymized evidence snippets, not raw PHI |
| BAA | Required for cloud LLM provider (Anthropic) if PHI included in queries |
| De-identification | Patient identifiers stripped before LLM prompt construction |

### 19.3 Data Privacy Architecture

The privacy architecture ensures that patient data never leaves the DGX Spark:

```
Patient Documents -> [Local OCR/NLP] -> [Local Milvus] -> [Local Query]
                                                              |
                                                              v
                                              [Anonymized evidence snippets]
                                                              |
                                                              v
                                              [Cloud LLM (Claude Sonnet 4.6)]
                                                              |
                                                              v
                                              [Synthesized response returned]
```

Only the anonymized evidence snippets (with patient identifiers removed) are sent to the cloud LLM. All raw clinical documents, extracted entities, and patient-specific collections remain on the local device.

---

## 20. DGX Compute Progression

### 20.1 DGX Spark (Current Target -- $3,999)

| Component | Specification | Autoimmune Agent Usage |
|-----------|-------------|----------------------|
| GPU | GB10 (Blackwell) | Embedding generation, NLP inference |
| Memory | 128 GB unified LPDDR5x | Milvus index (~6 GB), embedding model (~130 MB), NLP models (~500 MB) |
| CPU | 20 ARM cores (Grace) | Document ingestion, OCR, entity extraction |
| Storage | NVMe SSD | Milvus data, patient documents, audit logs |

**Memory budget for 10 concurrent patients (worst case):**
- Milvus reference collections: ~6 GB
- 10 patients x 50,000 vectors x 384 dims x 4 bytes: ~740 MB
- BGE-small-en-v1.5: 130 MB
- SciSpaCy NLP models: 500 MB
- LLM inference (Llama-3 8B NIM): ~16 GB
- Operating system and overhead: ~8 GB
- **Total: ~31.4 GB of 128 GB** (75% headroom)

### 20.2 Future: DGX Spark Cluster / DGX Station

For institutions processing hundreds of patients concurrently, horizontal scaling via multiple DGX Sparks or vertical scaling to DGX Station provides additional capacity. The Milvus distributed architecture (with external etcd and MinIO) supports seamless cluster scaling.

---

## 21. Implementation Roadmap

### 21.1 Three-Phase, 18-Week Plan

**Phase 1: Reference Knowledge Foundation (Weeks 1-6)**

| Week | Deliverable |
|------|------------|
| 1-2 | Milvus collection schemas (autoimmune_literature, autoimmune_trials, autoimmune_autoantibodies, autoimmune_hla_associations, autoimmune_guidelines) + PubMed ingest pipeline |
| 3-4 | Knowledge graph (7 dictionaries), query expansion maps (18 maps), entity alias resolution, disease activity scoring engines |
| 5-6 | FastAPI endpoints (query, compare, collections), Streamlit Evidence Explorer tab, basic RAG engine with parallel search, 150+ unit tests |

**Phase 2: Clinical Document Intelligence (Weeks 7-12)**

| Week | Deliverable |
|------|------------|
| 7-8 | Document ingestion pipeline: PDF extraction (PyMuPDF), OCR (Tesseract), document classification, section segmentation |
| 9-10 | Medical NLP entity extraction: laboratory values, medications, diagnoses, symptoms, vital signs. Patient-specific Milvus collections (clinical_documents, patient_labs, patient_timeline) |
| 11-12 | Patient timeline construction, lab trend analysis, Streamlit Document Ingest + Patient Timeline + Lab Trends tabs, integration tests |

**Phase 3: Advanced Analytics and Genomic Integration (Weeks 13-18)**

| Week | Deliverable |
|------|------------|
| 13-14 | Genomic-autoimmune correlation engine, HLA typing integration, non-HLA risk variant queries, autoimmune trigger conditions |
| 15-16 | Classification criteria scoring (SLE, RA, Sjogren's, PsA), overlap syndrome detection, POTS/dysautonomia evaluation, flare prediction engine |
| 17-18 | Biologic therapy optimizer with PGx, PDF/FHIR report generation, demo cases, end-to-end testing, documentation, docker-compose deployment |

---

## 22. Risk Analysis

### 22.1 Technical Risks

| Risk | Severity | Mitigation |
|------|----------|-----------|
| NLP entity extraction accuracy insufficient for clinical use | High | Validate against i2b2 benchmarks; implement confidence scores; flag low-confidence extractions for human review |
| OCR quality for scanned/faxed documents | Medium | Pre-process with image enhancement; support manual correction; track OCR confidence scores |
| Patient document volume exceeds Milvus capacity on single node | Medium | Implement patient partition pruning; archive inactive patient collections; scale to Milvus cluster if needed |
| Classification criteria scoring discordance with expert assessment | High | Validate against expert-adjudicated cases; implement "uncertain" category for borderline scores; always display raw criteria components |
| LLM hallucination of diagnostic conclusions | High | Ground all conclusions in cited evidence; implement confidence thresholds; append disclaimers; require clinician verification |

### 22.2 Clinical Risks

| Risk | Severity | Mitigation |
|------|----------|-----------|
| False positive autoimmune diagnosis leads to unnecessary immunosuppression | Critical | System provides diagnostic hypotheses, not diagnoses; always recommend confirmatory testing; clinician makes final determination |
| False negative: system misses autoimmune diagnosis | High | Multi-disease simultaneous evaluation reduces single-disease bias; overlap syndrome detection captures atypical presentations |
| PHI exposure via LLM API | Critical | De-identify all patient data before LLM submission; implement PHI detection guardrails; audit all outbound API calls |
| Patient self-diagnosis from system output | Medium | System designed for clinician use only; add role-based access control; include prominent disclaimers |

### 22.3 Operational Risks

| Risk | Severity | Mitigation |
|------|----------|-----------|
| HIPAA violation from inadequate audit logging | Critical | Comprehensive audit trail for all patient data access; regular compliance review; encryption at rest and in transit |
| Document ingestion pipeline overwhelmed during batch upload | Medium | Queue-based architecture with progress tracking; configurable batch sizes; async processing |
| Knowledge graph becomes stale | Low | Automated weekly literature refresh; quarterly guideline review; version tracking for all knowledge base entries |

---

## 23. Competitive Landscape

### 23.1 Why the Precision Autoimmune Agent Is Unique

No existing product combines all four of these capabilities:

1. **Clinical document ingestion at patient scale:** The system ingests thousands of patient documents and extracts structured data -- not just a single lab panel or questionnaire
2. **Multi-disease simultaneous evaluation:** The system evaluates evidence for 13+ autoimmune conditions in parallel, detecting overlaps that single-disease tools miss
3. **Genomic correlation:** HLA typing, non-HLA risk variants, and pharmacogenomic data are integrated directly into diagnostic and therapeutic reasoning
4. **$3,999 hardware:** The system runs on a desktop workstation, not a cloud platform with per-query pricing or enterprise licensing

### 23.2 Competitive Matrix

| Capability | Precision Autoimmune Agent | Exagen AVISE | Epic Cogito | Google Health | Commercial AI Platforms |
|-----------|--------------------------|-------------|-------------|---------------|----------------------|
| Clinical document ingestion | Yes (thousands) | No | Limited (single EHR) | No | No |
| Multi-disease evaluation | 13+ diseases | Lupus only | Custom build required | No | 1-3 diseases |
| Genomic correlation | HLA + non-HLA + PGx | No | No | Limited | No |
| Autoantibody interpretation | 14 antibody types | 10 lupus-specific | No | No | Variable |
| Disease activity scoring | 5 validated indices | SLE only | Custom | No | Variable |
| Longitudinal trend analysis | Full timeline | No | Basic charts | No | Limited |
| POTS/dysautonomia | Full evaluation | No | No | No | No |
| Overlap syndrome detection | 12 syndromes | No | No | No | No |
| Open-source | Apache 2.0 | Proprietary | Proprietary | Proprietary | Proprietary |
| Hardware cost | $3,999 (one-time) | Per-test fee | EHR license | Cloud pricing | $50K-500K/year |

---

## 24. Discussion

### 24.1 Technical Feasibility

The clinical document ingestion approach described in this paper is technically feasible with current technology. Medical NLP has matured significantly: SciSpaCy achieves >90% F1 on clinical entity extraction benchmarks, Tesseract OCR with medical vocabulary enhancement handles most clinical document formats, and FHIR R4 provides a standardized interoperability layer for structured data exchange. The multi-collection RAG architecture has been proven across five existing HCLS AI Factory agents, with the CAR-T Intelligence Agent demonstrating that 3.5+ million vectors can be searched in parallel across 11 collections in under 30 seconds on a single DGX Spark.

The primary technical challenge is NLP accuracy on the wide variety of clinical document formats encountered in real-world autoimmune patient records. Handwritten notes, faxed outside records, and non-standard laboratory report formats will require robust fallback handling and confidence scoring. The validation strategy (Section 18) addresses this through multi-benchmark evaluation and expert review.

### 24.2 Clinical Impact

The potential clinical impact of reducing the autoimmune diagnostic odyssey is substantial. For lupus alone, compressing the average 4.6-year diagnostic delay could prevent irreversible organ damage (nephritis, cardiovascular disease, neuropsychiatric complications) in tens of thousands of patients annually. For POTS patients, reducing the 5-7 year diagnostic delay addresses a condition that affects an estimated 1-3 million Americans, many of whom are unable to work or attend school during their diagnostic journey.

The system's ability to detect overlap syndromes -- conditions that by definition span multiple specialty domains -- addresses a diagnostic blind spot that no single-specialty tool can solve. A patient with MCTD features seen by a rheumatologist, a dermatologist, and a pulmonologist separately may never receive a unifying diagnosis; the Precision Autoimmune Agent evaluates all evidence simultaneously.

### 24.3 The Genomic-Autoimmune Frontier

The integration of genomic data with clinical document analysis represents a frontier capability. Most autoimmune patients never undergo HLA typing or whole exome sequencing unless a specific indication arises. As the cost of genomic testing continues to fall (whole exome sequencing now under $250), routine genomic profiling for autoimmune risk assessment becomes increasingly practical. The Precision Autoimmune Agent is designed to leverage this data when available, providing a bridge between clinical phenotyping and genetic risk assessment that few existing tools offer.

The POTS/hEDS/MCAS triad is a particularly compelling use case for genomic integration. TPSAB1 duplication (hereditary alpha-tryptasemia) is present in 4-6% of the general population and substantially increases MCAS risk. COL5A1 and COL3A1 variants are associated with Ehlers-Danlos subtypes. Identifying these genetic factors alongside clinical autonomic testing data creates a diagnostic framework that is far more powerful than either approach alone.

### 24.4 Limitations

1. **NLP accuracy ceiling:** Clinical NLP on real-world documents will not achieve 100% accuracy. Handwritten notes, non-standard abbreviations, and context-dependent interpretations will produce extraction errors that must be flagged and reviewed.
2. **Diagnosis vs. diagnostic support:** The system generates diagnostic hypotheses, not diagnoses. Classification criteria scores are calculated from available data, which may be incomplete. Clinician judgment remains essential.
3. **Limited to supported diseases:** The initial 13 autoimmune diseases and 10 dysautonomia conditions cover the most common entities but do not encompass all autoimmune conditions. Rare diseases (autoimmune retinopathy, susac syndrome, IgG4-related disease) will require future collection expansion.
4. **Cloud LLM dependency:** The synthesis step requires Claude API access. On-device LLM deployment via NIM would eliminate this dependency but with reduced generation quality.
5. **Single-patient optimization:** The current architecture is optimized for individual patient analysis, not population-level epidemiology. Cohort analysis features are a future extension.

---

## 25. Conclusion

### 25.1 Key Contributions

This paper has presented the architectural design and product requirements for the Precision Autoimmune Agent, a clinical document intelligence system that addresses the autoimmune diagnostic odyssey through six key innovations:

1. **Clinical document ingestion at scale:** A pipeline that ingests thousands of patient clinical documents (PDFs, FHIR bundles, C-CDA, free text), extracts structured clinical entities using medical NLP, and indexes them in patient-specific Milvus vector collections for semantic search.

2. **Multi-disease simultaneous evaluation:** 14 specialized Milvus collections with parallel search across autoimmune reference knowledge, patient-specific clinical data, and genomic evidence -- enabling evaluation of 13+ autoimmune conditions in a single query.

3. **Genomic-autoimmune correlation:** Automatic cross-referencing of clinical findings with HLA-disease associations (50+ conditions), non-HLA autoimmune risk variants (PTPN22, STAT4, IRF5, TNFAIP3, IL23R, CTLA4, TPSAB1), and pharmacogenomic data for therapy optimization.

4. **Longitudinal pattern recognition:** Biomarker trend analysis across years of laboratory data to detect flare-predictive patterns (rising anti-dsDNA, falling complement, rising calprotectin) that human review of individual lab reports would miss.

5. **POTS/dysautonomia integration:** Purpose-built evaluation for POTS, hEDS, MCAS, and associated autonomic dysfunction -- conditions that are systematically underdiagnosed due to their multi-system nature and overlap with autoimmune disease.

6. **Hardware democratization:** The complete system runs on a single NVIDIA DGX Spark ($3,999), with all patient data processed locally for HIPAA compliance -- bringing world-class autoimmune diagnostic intelligence to any clinic worldwide.

### 25.2 Future Directions

- **Multi-institutional federated learning** via NVIDIA FLARE for autoimmune pattern discovery across health systems without sharing patient data
- **Patient-reported outcomes integration** via smartphone apps for symptom tracking between clinic visits
- **Wearable device integration** for continuous heart rate monitoring (POTS) and activity tracking (RA, AS)
- **Natural language patient history intake** where patients describe symptoms conversationally and the system structures the data
- **Autoimmune disease prediction** using pre-symptomatic biomarker patterns and genetic risk scores to identify individuals at risk before clinical onset
- **Drug repurposing** by analyzing treatment response patterns across the autoimmune disease spectrum
- **International guideline support** for NICE (UK), DGRh (Germany), JCR (Japan), and APLAR (Asia-Pacific) guidelines

### 25.3 Closing Remarks

The autoimmune diagnostic odyssey is not inevitable. The data needed to diagnose most autoimmune patients already exists in their medical records -- it is simply fragmented across too many systems, too many specialists, and too many years for any human to synthesize. The Precision Autoimmune Agent demonstrates that a carefully designed clinical document intelligence system, built on multi-collection RAG with genomic correlation and running on accessible hardware, can transform this fragmented data landscape into actionable diagnostic intelligence.

For the 50 million Americans and 800 million people worldwide living with autoimmune diseases -- many still undiagnosed -- this technology has the potential to compress years of suffering into weeks of systematic evaluation. By open-sourcing this system under the Apache 2.0 license, we aim to make this capability available to every clinic, every rheumatologist, and every patient who needs it.

---

## 26. References

### Autoimmune Disease Epidemiology and Burden

1. Jacobson DL, Gange SJ, Rose NR, Graham NM. Epidemiology and estimated population burden of selected autoimmune diseases in the United States. *Clin Immunol Immunopathol*. 1997;84(3):223-243. doi:10.1006/clin.1997.4412

2. American Autoimmune Related Diseases Association. Autoimmune Disease Statistics. https://autoimmune.org/resource-center/autoimmune-statistics/

3. GBD 2019 Diseases and Injuries Collaborators. Global burden of 369 diseases and injuries in 204 countries and territories: a systematic analysis for the Global Burden of Disease Study 2019. *Lancet*. 2020;396(10258):1204-1222.

### Diagnostic Delay in Autoimmune Disease

4. Lupus Foundation of America. Lupus Facts and Statistics. Median diagnostic delay 4.6 years. https://www.lupus.org/resources/lupus-facts-and-statistics

5. Dysautonomia International. Diagnostic Delay Survey Results (2019). Median 5.9 years to POTS diagnosis.

6. Rubio-Tapia A, Kyle RA, Kaplan EL, et al. Increased prevalence and mortality in undiagnosed celiac disease. *Gastroenterology*. 2009;137(1):88-93. doi:10.1053/j.gastro.2009.03.059

7. Feldtkeller E, Khan MA, van der Heijde D, et al. Age at disease onset and diagnosis delay in HLA-B27 negative vs. positive patients with ankylosing spondylitis. *Rheumatol Int*. 2003;23(2):61-66.

### Classification Criteria

8. Aringer M, Costenbader K, Daikh D, et al. 2019 European League Against Rheumatism/American College of Rheumatology Classification Criteria for Systemic Lupus Erythematosus. *Arthritis Rheumatol*. 2019;71(9):1400-1412. doi:10.1002/art.40930

9. Aletaha D, Neogi T, Silman AJ, et al. 2010 Rheumatoid arthritis classification criteria: an American College of Rheumatology/European League Against Rheumatism collaborative initiative. *Arthritis Rheum*. 2010;62(9):2569-2581. doi:10.1002/art.27584

10. Shiboski CH, Shiboski SC, Seror R, et al. 2016 American College of Rheumatology/European League Against Rheumatism Classification Criteria for Primary Sjogren's Syndrome. *Arthritis Rheumatol*. 2017;69(1):35-45. doi:10.1002/art.39859

11. van den Hoogen F, Khanna D, Fransen J, et al. 2013 classification criteria for systemic sclerosis: an American College of Rheumatology/European League Against Rheumatism collaborative initiative. *Arthritis Rheum*. 2013;65(11):2737-2747.

### HLA and Autoimmune Genetics

12. Trowsdale J, Knight JC. Major histocompatibility complex genomics and human disease. *Annu Rev Genomics Hum Genet*. 2013;14:301-323. doi:10.1146/annurev-genom-091212-153455

13. Brown MA, Kenna T, Wordsworth BP. Genetics of ankylosing spondylitis -- insights into pathogenesis. *Nat Rev Rheumatol*. 2016;12(2):81-91. doi:10.1038/nrrheum.2015.133

14. Raychaudhuri S, Sandor C, Stahl EA, et al. Five amino acids in three HLA proteins explain most of the association between MHC and seropositive rheumatoid arthritis. *Nat Genet*. 2012;44(3):291-296. doi:10.1038/ng.1076

### POTS and Dysautonomia

15. Sheldon RS, Grubb BP 2nd, Olshansky B, et al. 2015 Heart Rhythm Society Expert Consensus Statement on the Diagnosis and Treatment of Postural Tachycardia Syndrome. *Heart Rhythm*. 2015;12(6):e41-e63. doi:10.1016/j.hrthm.2015.03.029

16. Vernino S, Bourne KM, Stiles LE, et al. Postural Orthostatic Tachycardia Syndrome (POTS): State of the Science and Clinical Care from a 2019 National Institutes of Health Expert Consensus Meeting. *Auton Neurosci*. 2021;235:102828.

17. Dahan S, Tomljenovic L, Shoenfeld Y. Postural Orthostatic Tachycardia Syndrome (POTS) -- A Novel Member of the Autoimmune Family. *Lupus*. 2016;25(4):339-342.

### Autoimmune Biomarkers and Disease Activity

18. van Gestel AM, Prevoo ML, van 't Hof MA, et al. Development and validation of the European League Against Rheumatism response criteria for rheumatoid arthritis. *Arthritis Rheum*. 1996;39(1):34-40. (DAS28)

19. Gladman DD, Ibanez D, Urowitz MB. Systemic lupus erythematosus disease activity index 2000. *J Rheumatol*. 2002;29(2):288-291. (SLEDAI-2K)

20. Garrett S, Jenkinson T, Kennedy LG, et al. A new approach to defining disease status in ankylosing spondylitis: the Bath Ankylosing Spondylitis Disease Activity Index. *J Rheumatol*. 1994;21(12):2286-2291. (BASDAI)

### Biologic Therapy and Pharmacogenomics

21. Smolen JS, Landewe RBM, Bijlsma JWJ, et al. EULAR recommendations for the management of rheumatoid arthritis with synthetic and biological disease-modifying antirheumatic drugs: 2019 update. *Ann Rheum Dis*. 2020;79(6):685-699.

22. Fanouriakis A, Kostopoulou M, Alunno A, et al. 2019 update of the EULAR recommendations for the management of systemic lupus erythematosus. *Ann Rheum Dis*. 2019;78(6):736-745.

### Clinical NLP and RAG Architecture

23. Neumann M, King D, Beltagy I, Ammar W. ScispaCy: Fast and Robust Models for Biomedical Natural Language Processing. *Proceedings of BioNLP*. 2019. (SciSpaCy)

24. Lewis P, Perez E, Piktus A, et al. Retrieval-Augmented Generation for Knowledge-Intensive NLP Tasks. *Advances in Neural Information Processing Systems*. 2020;33:9459-9474.

25. Wang J, Yi X, Guo R, et al. Milvus: A Purpose-Built Vector Data Management System. *Proceedings of the 2021 International Conference on Management of Data*. 2021:2614-2627. doi:10.1145/3448016.3457550

26. Xiao S, Liu Z, Zhang P, Muennighoff N. C-Pack: Packaged Resources To Advance General Chinese Embedding. 2023. arXiv:2309.07597 (BGE embedding model)

---

*Precision Autoimmune Agent -- HCLS AI Factory v0.1.0 (Pre-Implementation)*
*Apache 2.0 License | https://github.com/ajones1923/hcls-ai-factory*
