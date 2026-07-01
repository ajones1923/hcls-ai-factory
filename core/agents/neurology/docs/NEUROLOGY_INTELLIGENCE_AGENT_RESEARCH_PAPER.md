# Advancing Neurological AI: A Multi-Collection RAG Architecture and Product Requirements for the Neurology Intelligence Agent

**Author:** Adam Jones
**Date:** March 2026
**Version:** 0.1.0 (Pre-Implementation)
**License:** Apache 2.0

Part of the HCLS AI Factory -- an end-to-end precision medicine platform.
https://github.com/ajones1923/hcls-ai-factory

---

## Abstract

Neurological disorders represent the leading cause of disability and the second leading cause of death worldwide, affecting over 3.4 billion people -- 43% of the global population. The World Health Organization estimates that neurological conditions account for 443 million disability-adjusted life years (DALYs) annually, with Alzheimer's disease and other dementias, stroke, migraine, epilepsy, Parkinson's disease, and multiple sclerosis among the most prevalent and devastating. Despite rapid advances in neuroimaging, electrophysiology, genomics, and computational neuroscience, clinical neurology remains hampered by profound data fragmentation: imaging findings exist in PACS, EEG data in separate neurophysiology systems, genomic results in laboratory databases, cognitive assessments in paper-based or siloed electronic records, and treatment evidence scattered across thousands of publications and clinical trials.

This paper presents the architectural design, clinical rationale, and product requirements for the Neurology Intelligence Agent -- a multi-collection retrieval-augmented generation (RAG) system purpose-built for clinical neuroscience. The agent will unify 13 specialized Milvus vector collections spanning neuroimaging (structural MRI, functional MRI, diffusion tensor imaging, PET, SPECT), electrophysiology (EEG, EMG/NCS, evoked potentials), neurodegenerative disease management (Alzheimer's, Parkinson's, ALS, Huntington's), cerebrovascular disease, epilepsy, neuro-oncology, multiple sclerosis, movement disorders, headache medicine, neuromuscular disease, neurogenetics, clinical trials, and neuroradiology literature -- alongside a shared genomic_evidence collection containing 3.5 million variant vectors from the HCLS AI Factory genomics pipeline.

The system extends the proven multi-collection RAG architecture established by six existing intelligence agents in the HCLS AI Factory (Precision Biomarker, Precision Oncology, CAR-T, Imaging, Autoimmune, and the planned Cardiology agent), adapting it with neurology-specific clinical workflows, validated assessment scales, cross-modal neuroimaging-genomics triggers, and structured reporting aligned with AAN practice guidelines. Eight reference clinical workflows will cover the highest-impact neurological use cases: acute stroke triage, dementia evaluation, epilepsy focus localization, brain tumor grading, multiple sclerosis monitoring, Parkinson's disease assessment, headache classification, and neuromuscular disease evaluation.

The agent will deploy on a single NVIDIA DGX Spark ($4,699) using BGE-small-en-v1.5 embeddings (384-dimensional, IVF_FLAT, COSINE), Claude Sonnet 4.6 for evidence synthesis, and four NVIDIA NIM microservices for on-device inference. Licensed under Apache 2.0, the platform will democratize access to integrated neurological intelligence that currently requires tertiary academic medical center infrastructure and subspecialty expertise concentrated in a small number of institutions globally.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [The Neurological Data Challenge](#2-the-neurological-data-challenge)
3. [Clinical Landscape and Market Analysis](#3-clinical-landscape-and-market-analysis)
4. [Existing HCLS AI Factory Architecture](#4-existing-hcls-ai-factory-architecture)
5. [Neurology Intelligence Agent Architecture](#5-neurology-intelligence-agent-architecture)
6. [Milvus Collection Design](#6-milvus-collection-design)
7. [Clinical Workflows](#7-clinical-workflows)
8. [Cross-Modal Integration](#8-cross-modal-integration)
9. [NIM Integration Strategy](#9-nim-integration-strategy)
10. [Knowledge Graph Design](#10-knowledge-graph-design)
11. [Query Expansion and Retrieval Strategy](#11-query-expansion-and-retrieval-strategy)
12. [API and UI Design](#12-api-and-ui-design)
13. [Clinical Decision Support Engines](#13-clinical-decision-support-engines)
14. [Reporting and Interoperability](#14-reporting-and-interoperability)
15. [Product Requirements Document](#15-product-requirements-document)
16. [Data Acquisition Strategy](#16-data-acquisition-strategy)
17. [Validation and Testing Strategy](#17-validation-and-testing-strategy)
18. [Regulatory Considerations](#18-regulatory-considerations)
19. [DGX Compute Progression](#19-dgx-compute-progression)
20. [Implementation Roadmap](#20-implementation-roadmap)
21. [Risk Analysis](#21-risk-analysis)
22. [Competitive Landscape](#22-competitive-landscape)
23. [Discussion](#23-discussion)
24. [Conclusion](#24-conclusion)
25. [References](#25-references)

---

## 1. Introduction

### 1.1 The Neurological Disease Burden

Neurological disorders constitute an unprecedented global health challenge. According to the Global Burden of Disease (GBD) 2021 study and the WHO's 2024 update on neurological conditions:

- **3.4 billion people** (43% of the global population) are affected by neurological conditions
- **443 million DALYs** lost annually to neurological disorders
- **11 million deaths per year** directly attributable to neurological disease (second leading cause of death after cardiovascular disease)
- **Alzheimer's disease and other dementias** affect 55 million people globally, projected to reach 139 million by 2050
- **Stroke** kills 6.6 million people annually and is the leading cause of acquired disability in adults
- **Epilepsy** affects 50 million people worldwide, with 80% in low- and middle-income countries
- **Parkinson's disease** prevalence has doubled in the past 25 years to 8.5 million, the fastest-growing neurological disorder
- **Multiple sclerosis** affects 2.8 million people globally, with incidence rising in every world region
- **Migraine** is the second leading cause of years lived with disability globally, affecting 1.1 billion people
- **Brain tumors** account for 308,000 new diagnoses annually, with glioblastoma carrying a median survival of 15 months

The economic burden is equally staggering. In the United States alone, neurological diseases cost an estimated $789 billion annually in direct medical costs and lost productivity. Alzheimer's disease alone accounts for $345 billion in annual costs, projected to exceed $1 trillion by 2050. Globally, the economic impact of neurological disorders exceeds $2.5 trillion per year.

### 1.2 The Crisis in Neurological Expertise

Unlike cardiology, where advanced imaging and intervention are widely available, neurology faces a severe workforce shortage that compounds the data challenge:

- The United States has approximately **19,500 practicing neurologists** for a population of 340 million -- a ratio of 1:17,400
- The AAN projects a **shortfall of 11,000 neurologists** by 2031
- Average wait time for a new neurology appointment is **4-8 weeks**, reaching **6+ months** in rural areas
- Only **~200 Level 4 epilepsy centers** exist in the US, leaving most patients without access to specialized evaluation
- Subspecialty expertise (neuro-oncology, movement disorders, neurogenetics) is concentrated in fewer than 100 academic medical centers

This shortage makes AI-assisted clinical decision support not merely convenient but essential. A community neurologist managing 2,500+ patients cannot maintain current knowledge across the breadth of neurology -- from stroke protocols to epilepsy genetics to immunotherapy for MS. An intelligence agent that synthesizes evidence across the full neurological spectrum can function as an always-available subspecialty consultant.

### 1.3 The Opportunity for Integrated Neurological Intelligence

The neurological domain presents unique characteristics that make it exceptionally well-suited for an integrated AI intelligence agent:

1. **Multi-modal data convergence**: Neurology integrates structural imaging (MRI, CT), functional imaging (fMRI, PET, SPECT), electrophysiology (EEG, EMG/NCS, evoked potentials), cognitive/behavioral assessment (neuropsychological testing, rating scales), genomics (neurogenetic panels), cerebrospinal fluid analysis, and biomarkers -- all of which must be synthesized for accurate diagnosis and treatment.

2. **Pattern recognition dependency**: Many neurological diagnoses depend on recognizing complex spatiotemporal patterns -- EEG seizure localization, white matter lesion distribution in MS, dopaminergic deficit patterns in Parkinson's, brain tumor morphology and enhancement characteristics -- that are well-suited for AI augmentation.

3. **Longitudinal disease tracking**: Neurodegenerative diseases unfold over years to decades. Tracking hippocampal volume loss, lesion burden evolution, motor score progression, and cognitive decline trajectories requires systematic longitudinal analysis that AI excels at.

4. **Genomic revolution**: Neurogenetics has transformed from an academic curiosity to a clinical imperative. Over 1,000 genes are now associated with neurological disease, with genetic diagnosis changing management in epilepsy (SCN1A/Dravet), movement disorders (GBA/Parkinson's), dementia (APOE, PSEN1/2, MAPT), neuromuscular disease (SMN1/SMA, DMD), and brain tumors (IDH1/2, 1p/19q, MGMT).

5. **Therapeutic pipeline explosion**: The neuroscience drug pipeline has expanded dramatically with anti-amyloid antibodies (lecanemab, donanemab), antisense oligonucleotides (nusinersen, tofersen), gene therapies (onasemnogene for SMA, atidarsagene for MLD), and emerging CRISPR-based approaches -- creating urgent need for precision treatment selection.

6. **Strong cross-modal triggers**: A neuroimaging finding (e.g., hippocampal atrophy pattern suggestive of genetic frontotemporal dementia) frequently triggers genomic workup (FTD gene panel: MAPT, GRN, C9orf72), which may guide therapy selection and family counseling.

### 1.4 Our Contribution

This paper presents the complete architectural blueprint and product requirements for the Neurology Intelligence Agent, the seventh domain-specific intelligence agent in the HCLS AI Factory platform. Our contributions include:

- A **13-collection Milvus vector schema** designed for the full spectrum of neurological data: neuroimaging, electrophysiology, neurodegenerative disease, cerebrovascular disease, epilepsy, neuro-oncology, multiple sclerosis, movement disorders, headache, neuromuscular disease, neurogenetics, clinical trials, and literature
- **Eight reference clinical workflows** covering acute stroke triage, dementia evaluation, epilepsy focus localization, brain tumor grading, MS monitoring, Parkinson's assessment, headache classification, and neuromuscular evaluation
- A **neurology knowledge graph** with structured data on 40+ neurological conditions, 20+ neuroimaging protocols, 30+ validated clinical scales, 25+ drug classes, and 50+ AAN/EAN guideline recommendations
- **Cross-modal triggers** linking neuroimaging findings to genomic workup (epilepsy gene panels, dementia gene panels, movement disorder genetics, neuromuscular genetics) via the shared `genomic_evidence` collection
- **Clinical decision support engines** implementing validated neurological scales (NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6) and diagnostic algorithms
- A comprehensive **product requirements document** with user stories, acceptance criteria, and implementation prioritization
- **Deployment on a single NVIDIA DGX Spark** ($4,699), maintaining the platform's commitment to accessible AI

---

## 2. The Neurological Data Challenge

### 2.1 Data Fragmentation in Clinical Neuroscience

Clinical neuroscience generates data across at least sixteen distinct categories, each with its own structure, vocabulary, source systems, and interpretive frameworks:

1. **Structural Neuroimaging** -- Brain MRI sequences (T1, T2, FLAIR, T2*, SWI, post-contrast T1), spinal MRI, CT head (non-contrast, CTA, CTP), measurements (hippocampal volume, cortical thickness, ventricular size, white matter hyperintensity volume, lesion counts and volumes, brain parenchymal fraction).

2. **Functional Neuroimaging** -- fMRI (task-based and resting-state connectivity), PET (FDG for metabolism, amyloid PET with florbetapir/flutemetamol/florbetaben, tau PET with flortaucipir, dopamine transporter DaT-SPECT), SPECT (ictal/interictal perfusion for epilepsy).

3. **Diffusion Imaging** -- DTI (fractional anisotropy, mean diffusivity, axial/radial diffusivity), tractography (corticospinal tract, arcuate fasciculus, optic radiations), DWI for acute stroke (ADC maps, mismatch assessment), connectomics.

4. **Electrophysiology -- EEG** -- Scalp EEG (routine, ambulatory, continuous ICU monitoring, long-term video-EEG monitoring), intracranial EEG (subdural grids, stereo-EEG depth electrodes), quantitative EEG, event-related potentials, spectral analysis, source localization (dipole modeling, LORETA), seizure detection and classification.

5. **Electrophysiology -- EMG/NCS** -- Nerve conduction studies (motor, sensory, F-waves, H-reflexes), needle electromyography (insertional activity, spontaneous activity, MUAP analysis, recruitment), repetitive nerve stimulation (NMJ disorders), single-fiber EMG.

6. **Neurodegenerative Disease Data** -- Cognitive assessments (MoCA, MMSE, CDR, ADAS-Cog), motor scales (UPDRS, Hoehn & Yahr, ALSFRS-R), biomarkers (CSF Abeta-42, p-tau-181, neurofilament light chain, alpha-synuclein seed amplification), amyloid/tau PET status, genetic risk factors (APOE, LRRK2, GBA, C9orf72, SOD1).

7. **Cerebrovascular Disease Data** -- Stroke scales (NIHSS, mRS, ASPECTS), vascular imaging (CTA, MRA, digital subtraction angiography), perfusion imaging (CTP, MR perfusion), hemorrhage grading (ICH Score, Fisher Scale, Hunt-Hess), thrombolysis eligibility criteria, thrombectomy selection (DAWN, DEFUSE-3 criteria), secondary prevention protocols.

8. **Epilepsy Data** -- Seizure classification (ILAE 2017), epilepsy syndrome diagnosis, EEG findings (interictal epileptiform discharges, ictal patterns, HFOs), MRI findings (hippocampal sclerosis, cortical dysplasia, tumors, vascular malformations), ASM efficacy and side effects, surgical evaluation data (Wada testing, neuropsychological lateralization, concordance matrix), stimulation device data (VNS, RNS, DBS parameters).

9. **Neuro-Oncology Data** -- WHO 2021 CNS tumor classification (molecular-integrated), grading (Grade 1-4), molecular markers (IDH1/2, 1p/19q codeletion, MGMT promoter methylation, H3K27M, TERT, ATRX, BRAF, CDKN2A), treatment protocols (Stupp protocol, CCNU/lomustine, bevacizumab, tumor treating fields), RANO response criteria, survival data.

10. **Multiple Sclerosis Data** -- McDonald criteria (2017), MS phenotype (CIS, RRMS, SPMS, PPMS), EDSS scoring, MRI lesion metrics (new T2, enhancing, PRL, cortical lesions, spinal cord lesions, brain atrophy), OCT (RNFL thickness, ganglion cell layer), CSF (oligoclonal bands, IgG index, kappa free light chains), DMT efficacy data (platform therapies through BTK inhibitors), NEDA-3 status.

11. **Movement Disorders Data** -- Parkinson's disease (UPDRS parts I-IV, Hoehn & Yahr, dopaminergic medication LEDDs, DBS programming, DaT-SPECT), essential tremor, dystonia (Burke-Fahn-Marsden), Huntington's disease (UHDRS, CAG repeat length, prodromal markers), ataxia (SARA, BARS), functional movement disorders.

12. **Headache Medicine** -- ICHD-3 classification, migraine subtypes (with/without aura, chronic, vestibular, hemiplegic), medication overuse headache, cluster headache, trigeminal autonomic cephalalgias, secondary headache red flags, preventive efficacy data (CGRP mAbs, gepants, neuromodulation), disability scales (HIT-6, MIDAS).

13. **Neuromuscular Disease Data** -- NCS/EMG patterns (demyelinating vs axonal, generalized vs focal), myasthenia gravis (AChR/MuSK antibodies, QMG score, MGFA classification), GBS (Hughes scale, treatment criteria), CIDP (EFNS/PNS criteria), muscular dystrophies (CK levels, genetic confirmation, functional scales), motor neuron disease (El Escorial criteria, ALSFRS-R, FVC).

14. **Neurogenetics Data** -- Mendelian neurological disorders (>1,000 genes), pharmacogenomics (HLA-B*15:02 for carbamazepine, CYP2C19 for clopidogrel), polygenic risk scores (Alzheimer's, Parkinson's, epilepsy), variant interpretation (ACMG criteria adapted for neurological variants), genetic counseling considerations, gene therapy eligibility.

15. **Clinical Trials** -- ClinicalTrials.gov neurology entries (12,000+ active/completed), landmark trial results (EMERGE/CLARITY for lecanemab, TRAILBLAZER for donanemab, DAWN/DEFUSE-3 for thrombectomy, ARISE for antisense oligonucleotides), outcome measures, biomarker endpoints.

16. **Neurorehabilitation Data** -- Functional assessments (FIM, Barthel Index, mRS at 90 days), rehabilitation intensity and protocols, plasticity biomarkers, assistive technology, cognitive rehabilitation outcomes.

### 2.2 Why Existing Tools Fall Short

| Approach | Limitation |
|---|---|
| **PubMed search** | Keyword-based; misses semantic connections across neurology subspecialties; no imaging or electrophysiology data integration |
| **UpToDate / DynaMed** | Expert-curated but static; no patient-specific reasoning; no imaging AI; no genomic integration; subscription-based |
| **PACS with AI** (Aidoc, Viz.ai) | Single-modality (CT for stroke); no EEG, EMG, genomics, or longitudinal tracking; cloud-dependent |
| **EEG platforms** (Persyst, BESA) | Electrophysiology-only; no imaging correlation; no literature synthesis; expensive licensing |
| **Genetic platforms** (Invitae, GeneDx) | Genetics-only; no imaging or electrophysiology correlation; report delivery, not decision support |
| **General AI assistants** | No citation provenance; hallucination risk in high-stakes neurological decisions; no structured data; not guideline-aligned |

### 2.3 The Case for Multi-Collection RAG in Neurology

A neurologist evaluating a patient with new-onset seizures must simultaneously consider:

- **Structural imaging**: MRI showing hippocampal sclerosis, cortical dysplasia, or tumor
- **Electrophysiology**: EEG showing temporal intermittent rhythmic delta activity or epileptiform discharges
- **Functional imaging**: PET showing temporal hypometabolism concordant with EEG focus
- **Genomics**: If family history or drug-resistant epilepsy, genetic panel (SCN1A, KCNQ2, CDKL5, TSC1/2)
- **Guidelines**: ILAE seizure classification, AAN/AES treatment guidelines for first seizure
- **Trials**: Evidence for specific ASMs by seizure type and epilepsy syndrome
- **Clinical scales**: Seizure frequency, side effect burden, quality of life measures

No existing tool synthesizes all seven dimensions. A multi-collection RAG architecture enables parallel retrieval across all dimensions with a single query, followed by LLM synthesis into coherent clinical guidance.

---

## 3. Clinical Landscape and Market Analysis

### 3.1 Neurology AI Market

| Metric | Value | Source |
|---|---|---|
| Global neuro AI market (2024) | $1.7 billion | Markets and Markets |
| Projected market (2029) | $4.9 billion | Markets and Markets |
| CAGR (2024-2029) | 23.7% | Markets and Markets |
| FDA-cleared neuro AI devices (cumulative) | 120+ | FDA AI/ML database |
| Active neuro AI clinical trials | 380+ | ClinicalTrials.gov |
| Annual neuro AI publications | 3,800+ | PubMed (2025) |
| US neurologists | ~19,500 | AAN Census |
| Global neurologists | ~85,000 | WHO estimates |
| US epilepsy centers (Level 4) | ~200 | NAEC |
| US memory disorder clinics | ~500 | Alzheimer's Association |
| Annual US neurological costs | $789 billion | AAN Economic Burden Study |
| Global Alzheimer's drug pipeline | 140+ candidates | Alzheimer's Drug Discovery Foundation |

### 3.2 Competitive Analysis

| Competitor | Strengths | Gaps |
|---|---|---|
| **Viz.ai** (LVO stroke) | FDA-cleared, real-time CT triage, hospital networks | Single use case (stroke), no EEG/EMG, no genomics, no outpatient neuro, SaaS pricing |
| **Aidoc** (CT head) | Multi-pathology CT triage (hemorrhage, PE, C-spine) | Imaging-only, no longitudinal tracking, no treatment guidance |
| **Brainomix** (e-ASPECTS, CTA) | Stroke imaging AI, European market presence | Stroke-only, no other neurological conditions |
| **Persyst** (EEG) | Automated seizure detection, spike detection | EEG-only, no imaging correlation, no genomics, expensive licensing ($30K+) |
| **Natus/Nihon Kohden** (EEG) | Comprehensive EEG platforms, ICU monitoring | Hardware-dependent, no AI-driven interpretation, no cross-modal |
| **Biogen Digital Health** | MS monitoring, Parkinson's wearables | Disease-specific, not comprehensive, pharma-biased |
| **QMENTA** (Neuroimaging platform) | Brain volumetrics, cloud neuroimaging analysis | Research-focused, no clinical workflow integration, cloud-only |
| **Tempus** (Neuro-oncology) | Molecular profiling, clinical data | Proprietary, expensive, limited to neuro-oncology |

**Our differentiation**: The Neurology Intelligence Agent will be the only system combining (1) multi-modal neuroimaging AI, (2) EEG/EMG integration, (3) genomic cross-modal triggers, (4) literature RAG with citations, (5) validated neurological scales, (6) guideline-aligned decision support, and (7) on-device deployment -- all open-source at $4,699.

### 3.3 Target Users

| User Segment | Use Case | Pain Point Addressed |
|---|---|---|
| **Community neurologists** | Subspecialty-level decision support | 4-8 week wait for subspecialty referral |
| **Academic neurology** | Research and education | Fragmented data across systems |
| **Epilepsy centers** | Pre-surgical evaluation, surgical candidacy | Complex multi-modal concordance analysis |
| **Memory clinics** | Dementia differential diagnosis | Distinguishing AD from FTD, DLB, vascular |
| **Stroke centers** | Acute triage and secondary prevention | Time-critical decision support |
| **Neuro-oncology programs** | Molecular-integrated tumor classification | Rapid molecular result integration |
| **MS centers** | Treatment escalation decisions | NEDA monitoring, DMT selection |
| **Movement disorder clinics** | Parkinson's vs atypical parkinsonism | Complex differential diagnosis |
| **Clinical trial sites** | Patient screening, biomarker endpoints | Manual eligibility assessment |
| **Neurogenetics clinics** | Variant interpretation in neurological context | Limited neuro-specific annotation |

---

## 4. Existing HCLS AI Factory Architecture

### 4.1 Platform Overview

The HCLS AI Factory is a three-stage precision medicine platform running on NVIDIA DGX Spark:

```
Stage 1: Genomics Pipeline (Parabricks + DeepVariant)
    FASTQ -> VCF -> 3.56M annotated variants
         |
Stage 2: RAG/Chat Pipeline (Milvus + Claude)
    Variant interpretation, clinical significance
         |
Stage 3: Drug Discovery Pipeline (BioNeMo + DiffDock)
    Target -> Lead compound -> Docking -> Drug-likeness
```

Six existing intelligence agents extend this platform:

| Agent | Collections | Seed Vectors | Unique Capability |
|---|---|---|---|
| Precision Biomarker | 11 | 6,134 | Biological age calculators, biomarker panels |
| Precision Oncology | 10 | 609 | Molecular tumor board packets, trial matching |
| CAR-T Intelligence | 11 | 6,266 | CAR construct comparison, manufacturing optimization |
| Imaging Intelligence | 10 | 876 | NIM inference, DICOM workflows, 3D segmentation |
| Autoimmune Intelligence | 10 | ~500 | Autoantibody panels, flare prediction |
| Cardiology (planned) | 12 | ~1,530 | Risk calculators, GDMT optimization |

### 4.2 Shared Infrastructure

All agents share:

- **Milvus 2.4** vector database (IVF_FLAT, COSINE, 384-dim)
- **BGE-small-en-v1.5** embedding model (sentence-transformers)
- **Claude Sonnet 4.6** (Anthropic) primary LLM
- **`genomic_evidence`** collection (3,561,170 variants, read-only)
- **Docker Compose** orchestration
- **FastAPI** (REST) + **Streamlit** (UI) pattern
- **lib/hcls_common** shared library (23 modules)

### 4.3 Proven Patterns Adapted for Neurology

| Pattern | Proven In | Adaptation for Neurology |
|---|---|---|
| Multi-collection parallel search | All agents | 13 neurology-specific collections |
| Knowledge graph augmentation | CAR-T, Biomarker, Cardiology | Neurological conditions, drug classes, genetic disorders |
| Query expansion maps | CAR-T (12 maps), Cardiology (15 maps) | Neurology terminology (e.g., "seizure" -> epilepsy, convulsion, ictal, epileptiform, paroxysmal) |
| Comparative analysis | CAR-T, Imaging | "Lecanemab vs Donanemab", "Carbamazepine vs Lamotrigine" |
| Cross-modal genomic triggers | Imaging (Lung-RADS -> EGFR), Cardiology (DCM -> TTN) | Neuroimaging -> epilepsy/dementia gene panels |
| FHIR R4 export | Imaging, Cardiology | DiagnosticReport with neuro SNOMED/LOINC codes |
| NIM inference workflows | Imaging (4 NIMs) | Brain segmentation, volumetrics, lesion tracking |
| Validated clinical scales | Cardiology (ASCVD, CHA2DS2-VASc) | NIHSS, GCS, MoCA, UPDRS, EDSS |
| Sidebar guided tour | Imaging, Cardiology | Neurology demo flow |

---

## 5. Neurology Intelligence Agent Architecture

### 5.1 System Diagram

```
+==========================================================================+
|  PRESENTATION:  Streamlit Neurology Workbench (8529)                     |
|                 10 Tabs | Evidence | Workflows | Imaging | Scales        |
|                 FastAPI REST Server (8528)                                |
+==========================================================================+
                    |                            |
+==========================================================================+
|  INTELLIGENCE:   Neurology RAG Engine                                    |
|                  13-collection parallel search                           |
|                  Knowledge graph (40 conditions, 30 scales, 50 genes)    |
|                  Query expansion (16 maps, neuroscience terminology)     |
|                  Comparative analysis ("Lecanemab vs Donanemab")         |
|                  Clinical scale calculators (NIHSS, GCS, MoCA, etc.)    |
|                  Diagnostic algorithm engine                             |
+==========================================================================+
                    |                            |
+==========================================================================+
|  INFERENCE:      NIM Services (VISTA-3D, MAISI, VILA-M3, Llama-3)       |
|                  8 Clinical Workflows:                                    |
|                  Stroke | Dementia | Epilepsy | Neuro-Onc | MS |        |
|                  Parkinson's | Headache | Neuromuscular                   |
+==========================================================================+
                    |                            |
+==========================================================================+
|  DATA:           Milvus 2.4 (13 neuro collections + genomic)            |
|                  BGE-small-en-v1.5 (384-dim, IVF_FLAT, COSINE)          |
|                  PubMed, ClinicalTrials.gov, AAN/EAN guidelines          |
|                  Curated seed data (imaging, EEG, scales, genetics)      |
+==========================================================================+
```

### 5.2 Design Principles

1. **Diagnostic rigor**: Neurological diagnoses carry life-altering implications (epilepsy, dementia, brain tumor). Every output includes differential diagnosis, red flags, and confidence qualifiers.
2. **Longitudinal awareness**: Neurodegenerative diseases require tracking over months to years. Collections and workflows support temporal comparisons (baseline vs follow-up volumetrics, lesion evolution, scale score trajectories).
3. **Multi-modal concordance**: Epilepsy and dementia evaluation require concordance across imaging, electrophysiology, and clinical data. The architecture supports concordance matrices.
4. **Genomic integration by default**: Every significant neurological finding is checked against the neurogenetics knowledge base for actionable genetic associations.
5. **Graceful degradation**: Full functionality in mock mode without GPU or live NIM services.
6. **Pattern consistency**: Follows the same FastAPI + Streamlit + Milvus patterns as all other HCLS AI Factory agents.

### 5.3 Port Allocation

| Port | Service |
|---|---|
| 8528 | FastAPI REST Server |
| 8529 | Streamlit Neurology Workbench |
| 19530 | Milvus (shared) |
| 8520 | NIM LLM (shared) |
| 8530 | NIM VISTA-3D (shared) |
| 8531 | NIM MAISI (shared) |
| 8532 | NIM VILA-M3 (shared) |

---

## 6. Milvus Collection Design

### 6.1 Index Configuration

| Parameter | Value |
|---|---|
| Index type | IVF_FLAT |
| Metric | COSINE |
| nlist / nprobe | 1024 / 16 |
| Dimension | 384 |
| Embedding model | BAAI/bge-small-en-v1.5 |

### 6.2 Collection Schemas

#### Collection 1: `neuro_literature` -- ~3,500 records

Published neuroscience and neurology research papers, reviews, and meta-analyses.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | PubMed ID or unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| title | VARCHAR(500) | Paper title |
| text_chunk | VARCHAR(8000) | Abstract or text section |
| year | INT16 | Publication year |
| journal | VARCHAR(200) | Journal name (Neurology, Ann Neurol, Brain, Lancet Neurol, JAMA Neurol, Epilepsia) |
| neuro_domain | VARCHAR(100) | Neurological subdomain |
| modality | VARCHAR(50) | Imaging or diagnostic modality if applicable |
| study_type | VARCHAR(50) | RCT, meta-analysis, cohort, case-control, review |
| keywords | VARCHAR(500) | MeSH terms and author keywords |

**Source:** PubMed E-utilities with neurology/neuroscience MeSH filters.

#### Collection 2: `neuro_trials` -- ~600 records

Neurology clinical trials from ClinicalTrials.gov and landmark trial results.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | NCT number or trial identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| title | VARCHAR(500) | Official trial title |
| text_summary | VARCHAR(4000) | Trial summary including results |
| phase | VARCHAR(20) | Phase I-IV |
| status | VARCHAR(30) | Active, completed, recruiting |
| sponsor | VARCHAR(200) | Lead sponsor |
| neuro_domain | VARCHAR(100) | Stroke, epilepsy, dementia, MS, PD, neuro-onc, headache |
| intervention | VARCHAR(300) | Drug, device, or procedure tested |
| primary_endpoint | VARCHAR(300) | Primary outcome measure |
| enrollment | INT32 | Number of participants |
| start_year | INT16 | Year trial began |
| outcome_summary | VARCHAR(2000) | Key results (if completed) |
| landmark | BOOL | Is this a landmark trial |

**Source:** ClinicalTrials.gov V2 API with neurological condition filters.

#### Collection 3: `neuro_imaging` -- ~250 records

Neuroimaging protocols, findings, and measurements across all modalities.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Finding or protocol description |
| modality | VARCHAR(50) | mri, ct, pet, spect, dti, fmri |
| sequence_type | VARCHAR(100) | T1, T2, FLAIR, SWI, DWI, ADC, T1_post, DTI_FA |
| finding_category | VARCHAR(100) | Atrophy, lesion, enhancement, hemorrhage, infarction, mass |
| brain_region | VARCHAR(100) | Temporal, frontal, parietal, occipital, cerebellum, brainstem, spinal_cord |
| measurement_name | VARCHAR(100) | Hippocampal_volume, cortical_thickness, WMH_volume, lesion_count |
| measurement_value | VARCHAR(50) | Numeric value with units |
| reference_range | VARCHAR(100) | Normal range per age/sex norms |
| clinical_significance | VARCHAR(500) | Interpretation guidance |
| differential_diagnosis | VARCHAR(500) | What this finding suggests |

**Source:** Curated from AAN/ASNR guidelines, neuroradiology references, and brain atlas data.

#### Collection 4: `neuro_electrophysiology` -- ~200 records

EEG patterns, EMG/NCS findings, and electrophysiological data.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | EEG or EMG pattern description |
| test_type | VARCHAR(50) | EEG, EMG, NCS, SSEP, VEP, BAEP |
| pattern_name | VARCHAR(200) | Specific pattern (e.g., "3-Hz generalized spike-and-wave") |
| lateralization | VARCHAR(50) | Left, right, bilateral, generalized, multifocal |
| brain_region | VARCHAR(100) | Temporal, frontal, parietal, occipital, generalized |
| clinical_correlation | VARCHAR(500) | Associated conditions |
| urgency | VARCHAR(20) | Routine, urgent, emergent |
| classification | VARCHAR(100) | ACNS standardized terminology |
| differential_diagnosis | VARCHAR(500) | DDx list |

**Source:** Curated from ACNS guidelines, ABEM standards, and EP references.

#### Collection 5: `neuro_degenerative` -- ~200 records

Neurodegenerative disease management: Alzheimer's, Parkinson's, ALS, FTD, Huntington's.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Disease management recommendation |
| disease | VARCHAR(100) | Alzheimer's, Parkinson's, ALS, FTD, Huntington's, DLB, MSA, PSP, CBD |
| disease_stage | VARCHAR(50) | Preclinical, prodromal, mild, moderate, severe |
| biomarker | VARCHAR(100) | Amyloid_PET, tau_PET, NfL, CSF_Abeta42, DaT_SPECT |
| biomarker_status | VARCHAR(50) | Positive, negative, borderline |
| drug_class | VARCHAR(100) | Cholinesterase inhibitor, anti-amyloid, dopaminergic, riluzole |
| drug_name | VARCHAR(100) | Specific medication |
| clinical_scale | VARCHAR(100) | MoCA, CDR, UPDRS, ALSFRS-R, UHDRS |
| evidence_level | VARCHAR(20) | Level of evidence |
| guideline_source | VARCHAR(100) | AAN, NIA-AA, MDS, EFNS |
| genetic_association | VARCHAR(200) | Associated genes |

**Source:** AAN practice guidelines, NIA-AA criteria, MDS criteria.

#### Collection 6: `neuro_cerebrovascular` -- ~180 records

Stroke, cerebrovascular disease, and neurovascular conditions.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Stroke management or vascular neurology recommendation |
| stroke_type | VARCHAR(50) | Ischemic, hemorrhagic, SAH, TIA, CVT |
| vascular_territory | VARCHAR(100) | MCA, ACA, PCA, basilar, PICA, AChA |
| imaging_finding | VARCHAR(200) | DWI restriction, CTA occlusion, perfusion mismatch, hemorrhage |
| scale_name | VARCHAR(50) | NIHSS, ASPECTS, ICH_Score, Fisher, Hunt_Hess, mRS |
| scale_value | VARCHAR(50) | Score value and interpretation |
| treatment | VARCHAR(300) | tPA, thrombectomy, EVD, surgical evacuation |
| time_window | VARCHAR(100) | Eligibility time window |
| eligibility_criteria | VARCHAR(500) | Treatment selection criteria (DAWN, DEFUSE-3) |
| secondary_prevention | VARCHAR(500) | Long-term prevention strategy |

**Source:** AHA/ASA stroke guidelines, AAN practice parameters.

#### Collection 7: `neuro_epilepsy` -- ~200 records

Epilepsy classification, treatment, and surgical evaluation data.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Epilepsy management recommendation |
| seizure_type | VARCHAR(100) | Focal_aware, focal_impaired, focal_to_bilateral, generalized_tonic_clonic, absence, myoclonic |
| epilepsy_syndrome | VARCHAR(200) | JME, CAE, Dravet, West, Lennox-Gastaut, TLE-HS, BECTS |
| eeg_finding | VARCHAR(300) | Interictal and ictal EEG patterns |
| mri_finding | VARCHAR(300) | Hippocampal sclerosis, FCD, tumor, vascular malformation |
| asm_name | VARCHAR(100) | Anti-seizure medication name |
| asm_class | VARCHAR(100) | Na+ channel, SV2A, GABA, broad-spectrum |
| asm_first_line | BOOL | Is this first-line for this seizure type |
| surgical_candidate | VARCHAR(50) | Yes, no, pending_evaluation |
| genetic_etiology | VARCHAR(200) | SCN1A, KCNQ2, CDKL5, TSC1/2, SLC2A1 |
| evidence_level | VARCHAR(20) | Level of evidence |

**Source:** ILAE classification, AAN/AES practice guidelines, epilepsy genetics databases.

#### Collection 8: `neuro_oncology` -- ~150 records

Brain tumor classification, treatment, and molecular data.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Neuro-oncology recommendation |
| tumor_type | VARCHAR(200) | WHO 2021 integrated diagnosis |
| who_grade | VARCHAR(10) | 1, 2, 3, 4 |
| molecular_markers | VARCHAR(500) | IDH1/2, 1p/19q, MGMT, H3K27M, TERT, ATRX, BRAF |
| location | VARCHAR(100) | Supratentorial, infratentorial, spinal, sellar |
| imaging_characteristics | VARCHAR(500) | Enhancement pattern, diffusion, perfusion, spectroscopy |
| treatment_protocol | VARCHAR(500) | Surgery, radiation, chemotherapy regimen |
| response_criteria | VARCHAR(200) | RANO criteria, iRANO for immunotherapy |
| prognosis | VARCHAR(300) | Expected survival, prognostic factors |

**Source:** WHO 2021 CNS tumor classification, NCCN guidelines, RANO working group.

#### Collection 9: `neuro_ms` -- ~180 records

Multiple sclerosis diagnosis, monitoring, and treatment.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | MS management recommendation |
| ms_phenotype | VARCHAR(30) | CIS, RRMS, SPMS, PPMS |
| diagnostic_criteria | VARCHAR(300) | McDonald 2017 criteria met |
| mri_metric | VARCHAR(100) | New_T2, enhancing, PRL, brain_atrophy, spinal_lesions |
| mri_value | VARCHAR(50) | Numeric measurement |
| edss_score | FLOAT | Expanded Disability Status Scale (0-10) |
| dmt_name | VARCHAR(100) | Disease-modifying therapy name |
| dmt_category | VARCHAR(50) | Platform, moderate_efficacy, high_efficacy |
| neda_status | VARCHAR(20) | NEDA-3 achieved, not achieved |
| escalation_criteria | VARCHAR(500) | When to escalate therapy |
| monitoring_protocol | VARCHAR(500) | MRI frequency, JCV testing, labs |

**Source:** AAN MS practice guidelines, ECTRIMS/EAN guidelines, McDonald 2017.

#### Collection 10: `neuro_movement` -- ~150 records

Movement disorders: Parkinson's disease, essential tremor, dystonia, ataxia, Huntington's.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Movement disorder recommendation |
| disorder | VARCHAR(100) | Parkinson's, essential_tremor, dystonia, Huntington's, ataxia, MSA, PSP, CBD, DLB |
| clinical_feature | VARCHAR(200) | Specific motor or non-motor feature |
| scale_name | VARCHAR(100) | UPDRS, H&Y, UHDRS, SARA, Burke_Fahn_Marsden |
| scale_value | VARCHAR(50) | Score and interpretation |
| medication | VARCHAR(100) | Drug name |
| medication_ledd | FLOAT | Levodopa equivalent daily dose (mg) |
| surgical_option | VARCHAR(200) | DBS target (STN, GPi, VIM), FUS, LITT |
| genetic_association | VARCHAR(200) | LRRK2, GBA, SNCA, PARK2, PINK1, DJ-1 |
| dat_spect_result | VARCHAR(100) | Normal, reduced (pattern description) |

**Source:** MDS clinical criteria, AAN treatment guidelines.

#### Collection 11: `neuro_headache` -- ~120 records

Headache classification, treatment, and prevention.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Headache management recommendation |
| headache_type | VARCHAR(100) | Migraine_without_aura, migraine_with_aura, chronic_migraine, cluster, TTH, MOH |
| ichd3_code | VARCHAR(20) | ICHD-3 classification code |
| red_flags | VARCHAR(500) | SNOOP mnemonic: Systemic, Neurologic, Onset, Older, Previous_hx |
| acute_treatment | VARCHAR(300) | Triptans, NSAIDs, gepants, ditans, ergots |
| preventive_treatment | VARCHAR(300) | CGRP mAbs, topiramate, amitriptyline, propranolol, OnabotulinumtoxinA |
| neuromodulation | VARCHAR(200) | sTMS, eTNS, nVNS, REN |
| disability_score | VARCHAR(50) | HIT-6, MIDAS score and grade |
| frequency | VARCHAR(50) | Monthly headache days |

**Source:** ICHD-3 classification, AAN/AHS practice guidelines.

#### Collection 12: `neuro_neuromuscular` -- ~130 records

Neuromuscular disease: neuropathy, myopathy, NMJ disorders, motor neuron disease.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Neuromuscular disease recommendation |
| disease_category | VARCHAR(100) | Neuropathy, myopathy, NMJ_disorder, motor_neuron_disease |
| disease_name | VARCHAR(200) | GBS, CIDP, CMT, MG, ALS, SMA, DMD, myotonic_dystrophy |
| emg_ncs_pattern | VARCHAR(300) | Demyelinating, axonal, mixed, myopathic, NMJ_decrement |
| antibody | VARCHAR(100) | AChR, MuSK, anti-MAG, anti-GM1, anti-GQ1b |
| ck_level | VARCHAR(50) | Normal, elevated (range) |
| genetic_test | VARCHAR(200) | SMN1, DMD, PMP22, MFN2, TTR |
| treatment | VARCHAR(300) | IVIg, PLEX, rituximab, nusinersen, gene_therapy |
| functional_scale | VARCHAR(100) | ALSFRS-R, QMG, INCAT, CMTNS, FVC |
| prognosis | VARCHAR(300) | Expected course and milestones |

**Source:** AAN practice parameters, EFNS/PNS guidelines.

#### Collection 13: `neuro_guidelines` -- ~200 records

AAN, EAN, and subspecialty clinical practice guidelines and practice parameters.

| Field | Type | Description |
|---|---|---|
| id | VARCHAR(64) | Unique identifier |
| embedding | FLOAT_VECTOR(384) | BGE-small-en-v1.5 embedding |
| text_summary | VARCHAR(4000) | Guideline recommendation |
| guideline_name | VARCHAR(300) | Full guideline title |
| organization | VARCHAR(50) | AAN, EAN, AES, ILAE, MDS, AHA/ASA, NCCN |
| year | INT16 | Publication year |
| neuro_domain | VARCHAR(100) | Stroke, epilepsy, MS, dementia, movement, headache, NMD |
| recommendation_level | VARCHAR(30) | Level A/B/C/U (AAN) or Class I/IIa/IIb/III |
| key_recommendation | VARCHAR(2000) | Specific recommendation text |
| clinical_scenario | VARCHAR(500) | When this recommendation applies |

**Source:** AAN Practice Guideline library, EAN guidelines, ILAE/AES guidelines.

#### Collection 14: `genomic_evidence` -- 3,561,170 records (read-only)

Shared genomic variant collection from HCLS AI Factory Stage 1+2.

**Purpose:** Cross-modal triggers query this collection for neurologically relevant genes (APOE, PSEN1/2, APP, MAPT, GRN, C9orf72, SCN1A, LRRK2, GBA, SMN1, DMD, HTT, etc.).

### 6.3 Collection Search Weights

| Collection | Weight | Rationale |
|---|---|---|
| Literature | 0.14 | Largest corpus, broadest evidence base |
| Guidelines | 0.12 | Highest clinical authority |
| Trials | 0.10 | Primary evidence source |
| Neuroimaging | 0.10 | Central to neurological diagnosis |
| Neurodegenerative | 0.08 | High-impact chronic diseases |
| Cerebrovascular | 0.08 | Acute high-acuity conditions |
| Epilepsy | 0.07 | Complex diagnostic workup |
| Neuro-Oncology | 0.06 | Molecular-integrated diagnosis |
| MS | 0.06 | Growing therapeutic complexity |
| Movement Disorders | 0.05 | Subspecialty differential diagnosis |
| Electrophysiology | 0.05 | Diagnostic correlate |
| Headache | 0.04 | High-volume clinical need |
| Neuromuscular | 0.04 | Subspecialty domain |
| Genomic Evidence | 0.01 | Supplementary variant context |

### 6.4 Estimated Vector Counts

| Collection | Seed Records | Post-Ingest Target |
|---|---|---|
| neuro_literature | 250 | 3,500+ (PubMed ingest) |
| neuro_trials | 60 | 600+ (ClinicalTrials.gov) |
| neuro_imaging | 250 | 250 |
| neuro_electrophysiology | 200 | 200 |
| neuro_degenerative | 200 | 200 |
| neuro_cerebrovascular | 180 | 180 |
| neuro_epilepsy | 200 | 200 |
| neuro_oncology | 150 | 150 |
| neuro_ms | 180 | 180 |
| neuro_movement | 150 | 150 |
| neuro_headache | 120 | 120 |
| neuro_neuromuscular | 130 | 130 |
| neuro_guidelines | 200 | 200 |
| **Total (owned)** | **~2,270** | **~6,060+** |
| genomic_evidence (read-only) | -- | 3,561,170 |

---

## 7. Clinical Workflows

### 7.1 Workflow Architecture

All workflows follow the established `BaseImagingWorkflow` pattern: `preprocess -> infer -> postprocess -> WorkflowResult`. Each workflow supports full mock mode with clinically realistic synthetic results.

### 7.2 Eight Reference Workflows

#### Workflow 1: Acute Stroke Triage

| Attribute | Value |
|---|---|
| Workflow ID | `acute_stroke_triage` |
| Input | CT head (non-contrast), CTA, CTP, clinical data |
| Target Latency | < 90 seconds |
| Models | VISTA-3D (hemorrhage segmentation, infarct core), perfusion analysis |
| Key Outputs | Stroke type (ischemic vs hemorrhagic), ASPECTS score, LVO detection, perfusion mismatch ratio, NIHSS estimate, tPA/thrombectomy eligibility |
| Severity Routing | LVO detected or NIHSS >= 6 -> Emergent stroke alert |
| Cross-Modal Trigger | Cryptogenic stroke + age < 50 -> thrombophilia/CADASIL panel (NOTCH3, COL4A1/2) |
| Guideline Alignment | AHA/ASA 2019 Acute Ischemic Stroke Guidelines, 2022 Update |

**Clinical Decision Logic:**

```
Acute Ischemic Stroke Triage:

1. CT Non-Contrast:
   - Hemorrhage? -> YES -> ICH pathway (ICH Score, surgical evaluation)
   - ASPECTS score -> < 6: poor candidate for intervention
                    -> >= 6: proceed to vascular imaging

2. CTA:
   - LVO present? (ICA, M1, M2, basilar)
   -> YES + NIHSS >= 6 + ASPECTS >= 6:
      -> Last known well < 6 hours: Direct to thrombectomy
      -> Last known well 6-24 hours: Apply DAWN/DEFUSE-3 criteria
         DAWN: Clinical-core mismatch (NIHSS >= 10, core < 31 mL for age < 80)
         DEFUSE-3: Perfusion mismatch ratio >= 1.8, core < 70 mL
   -> NO: Medical management (tPA if < 4.5 hours, antiplatelet)

3. tPA Eligibility (< 4.5 hours from last known well):
   - Age >= 18
   - Measurable neurological deficit (NIHSS >= 4 typically)
   - No contraindications (recent surgery, active bleeding, INR > 1.7)
   - BP controllable to < 185/110

4. Hemorrhagic Stroke:
   ICH Score (0-6):
     GCS 3-4: +2, GCS 5-12: +1, GCS 13-15: +0
     ICH volume >= 30 mL: +1
     IVH present: +1
     Infratentorial: +1
     Age >= 80: +1
   30-day mortality: Score 0=0%, 1=13%, 2=26%, 3=72%, 4=97%, 5-6=100%
```

#### Workflow 2: Dementia Evaluation

| Attribute | Value |
|---|---|
| Workflow ID | `dementia_evaluation` |
| Input | Brain MRI, cognitive testing scores, biomarkers |
| Target Latency | < 3 minutes |
| Models | VISTA-3D (hippocampal volumetry, cortical thickness, ventricular volume, WMH quantification) |
| Key Outputs | Atrophy pattern classification, hippocampal volume percentile, Fazekas WMH score, differential diagnosis (AD vs FTD vs DLB vs vascular), NIA-AA ATN classification, treatment eligibility |
| Severity Routing | Rapid progression or unusual pattern -> Urgent neuro referral |
| Cross-Modal Trigger | Pattern suggestive of genetic FTD -> FTD gene panel (MAPT, GRN, C9orf72); Early-onset AD (<65) -> AD gene panel (PSEN1, PSEN2, APP); DLB features -> GBA |
| Guideline Alignment | NIA-AA 2024 Revised Criteria, AAN Mild Cognitive Impairment Guidelines |

**Atrophy Pattern Recognition:**

| Pattern | Brain Region | Suggestive Diagnosis |
|---|---|---|
| Medial temporal atrophy (MTA) | Hippocampus, entorhinal cortex | Alzheimer's disease |
| Frontal > temporal atrophy | Frontal lobes, anterior temporal | Behavioral variant FTD |
| Asymmetric left temporal | Left temporal pole, fusiform | Semantic variant PPA |
| Asymmetric left perisylvian | Left inferior frontal, insula | Nonfluent variant PPA |
| Posterior cortical atrophy | Parieto-occipital | Posterior cortical atrophy (AD variant) |
| Midbrain atrophy ("hummingbird") | Midbrain tegmentum | PSP |
| Pontine/cerebellar atrophy ("hot cross bun") | Pons, cerebellum | MSA-C |
| Caudate atrophy | Caudate heads | Huntington's disease |
| Diffuse cortical + hippocampal | Global | DLB (less MTA than AD) |

**NIA-AA ATN Biomarker Framework:**

```
A (Amyloid): Amyloid PET or CSF Abeta-42/40 ratio
  A+: Abnormal amyloid -> Alzheimer's continuum
  A-: Normal amyloid -> Non-Alzheimer's pathway

T (Tau): Tau PET or CSF p-tau-181/217
  T+: Abnormal tau -> AD pathologic change
  T-: Normal tau

N (Neurodegeneration): MRI atrophy, FDG-PET, CSF NfL, CSF total-tau
  N+: Evidence of neurodegeneration
  N-: No neurodegeneration

Classification:
  A-T-N-: Normal biomarkers
  A+T-N-: Preclinical Alzheimer's (asymptomatic)
  A+T+N-: Alzheimer's pathologic change
  A+T+N+: Alzheimer's disease
  A-T-N+: Non-AD neurodegeneration (suspected SNAP)
  A-T+N+: Non-AD tauopathy (consider FTLD-tau)
```

#### Workflow 3: Epilepsy Focus Localization

| Attribute | Value |
|---|---|
| Workflow ID | `epilepsy_focus_localization` |
| Input | Brain MRI, EEG data, clinical seizure semiology, PET/SPECT if available |
| Target Latency | < 3 minutes |
| Models | VISTA-3D (hippocampal volumetry, lesion detection), volumetric comparison |
| Key Outputs | ILAE seizure classification, epilepsy syndrome identification, MRI-EEG concordance assessment, surgical candidacy evaluation, ASM recommendation by seizure type |
| Severity Routing | Status epilepticus or drug-resistant epilepsy -> Epilepsy center referral |
| Cross-Modal Trigger | Drug-resistant epilepsy + age < 25 -> epilepsy gene panel (SCN1A, KCNQ2, CDKL5, TSC1/2, DEPDC5, SLC2A1, PCDH19) |
| Guideline Alignment | ILAE 2017 Classification, AAN/AES 2018 First Seizure Guidelines |

**Concordance Matrix for Surgical Evaluation:**

```
Data Source          | Finding               | Lateralization | Region
---------------------|----------------------|----------------|--------
Seizure semiology    | Aura type, motor sx   | L / R / Bilat  | Temporal/Frontal/etc
Interictal EEG       | Spike location         | L / R / Bilat  | Temporal/Frontal/etc
Ictal EEG            | Seizure onset          | L / R / Bilat  | Temporal/Frontal/etc
MRI                  | Structural lesion      | L / R / Bilat  | Temporal/Frontal/etc
PET                  | Hypometabolism         | L / R / Bilat  | Temporal/Frontal/etc
Neuropsych           | Memory lateralization  | L / R           | Temporal/Frontal/etc

Concordance Score: 6/6 = Excellent surgical candidate
                   4-5/6 = Good candidate, consider invasive monitoring
                   < 4/6 = Invasive monitoring required (sEEG/grids)
```

#### Workflow 4: Brain Tumor Grading and Molecular Classification

| Attribute | Value |
|---|---|
| Workflow ID | `brain_tumor_grading` |
| Input | Brain MRI (T1, T1+C, T2, FLAIR, DWI, perfusion), molecular data if available |
| Target Latency | < 3 minutes |
| Models | VISTA-3D (tumor segmentation -- enhancing, non-enhancing, edema), perfusion analysis |
| Key Outputs | WHO 2021 integrated diagnosis prediction, tumor volume (enhancing + FLAIR), eloquent cortex proximity, molecular marker predictions (IDH, MGMT, 1p/19q), treatment protocol recommendation, RANO baseline measurements |
| Severity Routing | Large mass effect or midline shift -> Emergent neurosurgery |
| Cross-Modal Trigger | Suspected glioma -> molecular panel (IDH1/2, ATRX, TERT, CDKN2A, H3K27M) |
| Guideline Alignment | WHO 2021 CNS Tumor Classification, NCCN CNS Cancers Guidelines |

**WHO 2021 Integrated Classification:**

| Histology | IDH | 1p/19q | ATRX | Grade | Integrated Diagnosis | Prognosis |
|---|---|---|---|---|---|---|
| Astrocytic | Mutant | Intact | Lost | 2-4 | Astrocytoma, IDH-mutant | Better (mOS 5-10y) |
| Oligodendroglial | Mutant | Codeleted | Retained | 2-3 | Oligodendroglioma, IDH-mutant, 1p/19q-codeleted | Best (mOS 12-15y) |
| Astrocytic/GBM features | Wildtype | Intact | Retained | 4 | Glioblastoma, IDH-wildtype | Worst (mOS 15 mo) |
| Diffuse midline | Wildtype | -- | -- | 4 | Diffuse midline glioma, H3K27M-altered | Poor (mOS 9-11 mo) |

#### Workflow 5: Multiple Sclerosis Monitoring

| Attribute | Value |
|---|---|
| Workflow ID | `ms_monitoring` |
| Input | Brain and spinal cord MRI (current + prior), clinical data, EDSS, labs |
| Target Latency | < 3 minutes |
| Models | VISTA-3D (lesion segmentation, lesion comparison, brain volume change) |
| Key Outputs | New T2 lesion count, enhancing lesion count, total lesion volume change, brain volume change (annualized), NEDA-3 status, DMT escalation recommendation, PML risk (JCV index) |
| Severity Routing | Multiple new enhancing lesions on DMT -> Treatment failure, urgent MS specialist |
| Cross-Modal Trigger | Atypical MS features (progressive from onset, bilateral optic neuritis) -> NMOSD/MOGAD antibody panel + AQP4/MOG genetics |
| Guideline Alignment | AAN DMT Guidelines 2018, McDonald 2017 Criteria |

**NEDA-3 Assessment:**

```
No Evidence of Disease Activity (NEDA-3):
  1. No relapses in past 12 months
  2. No new or enlarging T2 lesions on MRI
  3. No confirmed disability progression (EDSS stable or improved)

All three criteria must be met = NEDA-3 achieved
Any criterion failed = NEDA-3 not achieved -> Consider DMT escalation

DMT Escalation Ladder:
  Platform therapies: Interferons, glatiramer acetate, teriflunomide, DMF
  Moderate efficacy: Fingolimod, cladribine
  High efficacy: Natalizumab (JCV risk), ocrelizumab, ofatumumab, alemtuzumab
  Emerging: BTK inhibitors (tolebrutinib, fenebrutinib -- trials ongoing)

Escalation triggers:
  - >= 1 relapse on current DMT
  - >= 2 new T2 or >= 1 enhancing lesion
  - Confirmed EDSS worsening (>= 1.0 if EDSS <= 5.5, >= 0.5 if EDSS > 5.5)
  - Suboptimal brain volume loss (> 0.4%/year)
```

#### Workflow 6: Parkinson's Disease Assessment

| Attribute | Value |
|---|---|
| Workflow ID | `parkinsons_assessment` |
| Input | Clinical motor exam, DaT-SPECT, brain MRI, cognitive testing |
| Target Latency | < 2 minutes |
| Key Outputs | MDS clinical diagnostic criteria assessment, UPDRS Part III motor score, Hoehn & Yahr stage, tremor-dominant vs PIGD classification, DaT-SPECT interpretation, medication optimization (LEDD calculation), DBS candidacy assessment, red flags for atypical parkinsonism |
| Severity Routing | Red flags for MSA/PSP/CBD -> Movement disorder specialist referral |
| Cross-Modal Trigger | Early-onset PD (<50) or Jewish ancestry -> PD gene panel (LRRK2, GBA, SNCA, PARK2/Parkin, PINK1, DJ-1) |
| Guideline Alignment | MDS Clinical Diagnostic Criteria, AAN PD Treatment Guidelines |

**Atypical Parkinsonism Red Flags:**

| Red Flag | Suggests |
|---|---|
| Poor levodopa response | MSA, PSP, CBD, vascular parkinsonism |
| Early falls (within 1 year) | PSP |
| Vertical supranuclear gaze palsy | PSP |
| Cerebellar ataxia | MSA-C |
| Severe autonomic failure | MSA |
| Asymmetric cortical signs (apraxia, alien limb) | CBD |
| Rapid cognitive decline | DLB, CBD, CJD |
| Wheelchair-bound within 5 years | PSP, MSA |
| Midbrain atrophy ("hummingbird sign") | PSP |
| Pontine "hot cross bun" sign | MSA-C |
| Putaminal rim sign (T2*) | MSA-P |

#### Workflow 7: Headache Classification and Management

| Attribute | Value |
|---|---|
| Workflow ID | `headache_classification` |
| Input | Headache characteristics, associated symptoms, exam findings, imaging if indicated |
| Target Latency | < 30 seconds |
| Key Outputs | ICHD-3 classification, red flag assessment (SNOOP), imaging recommendation, acute treatment plan, preventive treatment candidacy, disability score (HIT-6, MIDAS), CGRP therapy eligibility |
| Severity Routing | Red flags (thunderclap, new neurological deficit, papilledema) -> Emergent imaging |
| Guideline Alignment | ICHD-3, AAN Migraine Prevention Guidelines, AHS Position Statements |

**SNOOP Red Flag Mnemonic:**

```
S - Systemic symptoms/signs (fever, weight loss, cancer, HIV, pregnancy)
N - Neurologic symptoms/signs (confusion, focal deficit, papilledema, seizure)
O - Onset: Sudden/thunderclap (< 1 minute to peak -> SAH until proven otherwise)
O - Older: New headache onset after age 50 (temporal arteritis, mass, subdural)
P - Previous headache history: First or worst headache, change in character
    + Pattern change, Positional component, Precipitated by Valsalva

Any SNOOP flag present -> Imaging and further workup required before primary headache diagnosis
```

#### Workflow 8: Neuromuscular Disease Evaluation

| Attribute | Value |
|---|---|
| Workflow ID | `neuromuscular_evaluation` |
| Input | Clinical presentation, EMG/NCS data, serological data, genetic data |
| Target Latency | < 2 minutes |
| Key Outputs | EMG/NCS pattern classification (axonal vs demyelinating, motor vs sensory, proximal vs distal), differential diagnosis ranked by probability, antibody testing recommendations, genetic testing recommendations, treatment plan, functional assessment |
| Severity Routing | Rapidly progressive weakness (GBS) or respiratory compromise -> ICU admission |
| Cross-Modal Trigger | Hereditary neuropathy pattern -> CMT gene panel (PMP22, MFN2, GJB1, MPZ); Suspected SMA -> SMN1 testing; Suspected muscular dystrophy -> dystrophy panel (DMD, DMPK, CNBP) |
| Guideline Alignment | AAN Practice Parameters for GBS, CIDP, MG, ALS |

**EMG/NCS Pattern Recognition:**

| Pattern | NCS Motor | NCS Sensory | EMG | Diagnosis |
|---|---|---|---|---|
| Diffuse demyelinating motor + sensory | Slow CV, prolonged F-waves, temporal dispersion | Slow CV, low amplitude | Reduced recruitment | CIDP, GBS (AIDP) |
| Length-dependent axonal sensorimotor | Low amplitude distally | Low amplitude distally | Fibrillations, large MUAPs | Diabetic/metabolic neuropathy |
| Motor-predominant axonal | Low amplitude, normal CV | Normal | Fibrillations, fasciculations, large MUAPs, reduced recruitment | ALS, MMN |
| NMJ decrement | Decrement >10% at 3Hz RNS | Normal | Unstable MUAPs, variable firing | Myasthenia gravis, LEMS |
| Myopathic | Normal | Normal | Short, small, polyphasic MUAPs, early recruitment | Inflammatory myopathy, dystrophy |
| Uniform demyelinating motor + sensory | Uniformly slow CV | Uniformly slow | Reduced recruitment | CMT1A (hereditary) |

---

## 8. Cross-Modal Integration

### 8.1 Neurogenetics Triggers

The Neurology Intelligence Agent implements cross-modal triggers that automatically query the shared `genomic_evidence` collection (3.5M variants) when clinical or imaging findings suggest a genetic etiology:

| Trigger Condition | Gene Panel Queried | Clinical Rationale |
|---|---|---|
| Early-onset AD (<65) or family history | PSEN1, PSEN2, APP, APOE | 1-5% of AD is autosomal dominant; APOE e4 major risk factor |
| FTD features (behavioral/language) | MAPT, GRN, C9orf72 | 30-50% of FTD is genetic; C9orf72 is most common cause |
| Early-onset PD (<50) or family history | LRRK2, GBA, SNCA, PARK2, PINK1, DJ-1 | GBA variants in 5-15% of PD; affects prognosis and therapy |
| Drug-resistant epilepsy, age < 25 | SCN1A, KCNQ2, CDKL5, TSC1/2, DEPDC5, SLC2A1, PCDH19 | Genetic diagnosis changes ASM selection (e.g., avoid Na+ channel blockers in SCN1A) |
| Unexplained neuropathy (young onset) | PMP22, MFN2, GJB1, MPZ, TTR | CMT affects 1:2,500; TTR amyloid is treatable |
| Suspected SMA/muscular dystrophy | SMN1, DMD, DMPK, CNBP | Gene therapies available (nusinersen, onasemnogene) |
| Huntington's features (chorea + cognitive) | HTT (CAG repeat) | Definitive diagnosis; predictive testing for family members |
| Cryptogenic stroke, age < 50 | NOTCH3 (CADASIL), COL4A1/2, MTHFR, Factor V Leiden | Monogenic small vessel disease, thrombophilia |
| Unexplained leukoencephalopathy | CSF1R, EIF2B1-5, LMNB1, AARS2, DARS2 | Adult-onset leukodystrophies are underdiagnosed |
| Atypical MS / NMOSD features | AQP4, MOG (serological), HLA-DRB1 | NMOSD requires different treatment than MS |

### 8.2 Neuroimaging -> Genomics -> Therapeutics Pipeline

```
Neuroimaging Finding (MRI, CT, PET, DaT-SPECT)
    |
    v
[Cross-Modal Trigger] -- Clinical criteria met?
    |                         |
    YES                       NO
    |                         |
    v                         v
[Query genomic_evidence]   [Standard RAG response]
(3.5M variant vectors)
    |
    v
[Variant Annotation]
ClinVar pathogenicity, AlphaMissense score, ACMG classification
    |
    v
[Precision Therapeutics Selection]
Gene-specific therapy matching:
  SCN1A+ -> Avoid carbamazepine, phenytoin; use stiripentol + VPA + CLB
  GBA+ PD -> Consider GBA-targeted therapies (venglustat, ambroxol trials)
  SMN1 del -> Nusinersen, onasemnogene, risdiplam
  TTR+ amyloid neuropathy -> Tafamidis, patisiran, vutrisiran
    |
    v
[HCLS AI Factory Stage 3: Drug Discovery]
Novel compound generation for undruggable targets
    |
    v
[Clinical Output]
FHIR R4 DiagnosticReport with genomic enrichment
```

### 8.3 Integration with Other Agents

| Integration | Direction | Data Flow |
|---|---|---|
| Imaging Agent | Bidirectional | Shares brain MRI/CT workflows; receives DICOM routing for neuro studies |
| Cardiology Agent | Bidirectional | Stroke-cardiology interface (AF detection -> stroke prevention); cardio-embolic stroke workup |
| Precision Biomarker Agent | Inbound | Receives NfL, tau, amyloid biomarker reference ranges |
| Precision Oncology Agent | Bidirectional | Neuro-oncology molecular data; CNS metastasis management |
| Genomics Pipeline | Read-only | Queries `genomic_evidence` for neurological gene variants |
| Drug Discovery Pipeline | Outbound | Sends confirmed genetic targets for compound screening |

---

## 9. NIM Integration Strategy

### 9.1 Shared NIM Services

| NIM | Port | Neurology Application |
|---|---|---|
| VISTA-3D | 8530 | Brain segmentation (hippocampus, ventricles, cortex, white matter, tumors, lesions), volumetrics |
| MAISI | 8531 | Synthetic brain MRI generation for training, rare pathology simulation |
| VILA-M3 | 8532 | Brain image interpretation, EEG pattern recognition, radiology report assistance |
| Llama-3 8B | 8520 | Clinical reasoning fallback when Claude API unavailable |

### 9.2 VISTA-3D Neurology Applications

- **Hippocampal volumetry**: Bilateral hippocampal segmentation for dementia evaluation; age/sex-percentile comparison
- **White matter lesion quantification**: FLAIR lesion segmentation for MS monitoring, small vessel disease grading
- **Brain tumor segmentation**: Enhancing tumor, non-enhancing tumor, peritumoral edema volumes for RANO criteria
- **Brain parenchymal fraction**: Global and regional atrophy quantification for neurodegeneration tracking
- **Ventricular volume**: Hydrocephalus assessment, NPH evaluation
- **Cortical thickness mapping**: Regional cortical atrophy patterns for differential diagnosis

---

## 10. Knowledge Graph Design

### 10.1 Graph Structure

| Dictionary | Entries | Content |
|---|---|---|
| Neurological Conditions | ~45 | ICD-10, diagnostic criteria, imaging patterns, scales, genetic associations |
| Clinical Scales | ~30 | Input variables, scoring, interpretation, clinical use |
| Drug Classes | ~35 | Mechanism, dosing, interactions, specific neuro indications/contraindications |
| Neuroimaging Patterns | ~25 | Modality, sequence, pattern description, differential diagnosis |
| Neurological Genes | ~60 | Gene symbol, inheritance, phenotype, genetic testing indications, treatment implications |
| EEG/EMG Patterns | ~20 | Pattern name, morphology, clinical correlation, urgency |

### 10.2 Example Knowledge Graph Entries

**Condition: Dravet Syndrome**

```python
{
    "name": "Dravet Syndrome",
    "icd10": "G40.409",
    "aliases": ["SMEI", "severe myoclonic epilepsy of infancy"],
    "prevalence": "1:15,000 - 1:40,000",
    "inheritance": "De novo (90%), autosomal dominant",
    "age_of_onset": "4-8 months, first febrile seizure",
    "diagnostic_criteria": {
        "clinical": "Prolonged febrile seizures before age 1, followed by afebrile seizures",
        "eeg": "Normal initially, then generalized spike-wave, photosensitivity",
        "mri": "Normal initially, later hippocampal sclerosis possible",
        "genetic": "SCN1A pathogenic variant (>80%)"
    },
    "genes": ["SCN1A"],
    "treatment": {
        "effective": ["Valproic acid", "Clobazam", "Stiripentol", "Cannabidiol (Epidiolex)", "Fenfluramine (Fintepla)"],
        "contraindicated": ["Carbamazepine", "Oxcarbazepine", "Phenytoin", "Lamotrigine", "Vigabatrin"],
        "contraindication_reason": "Sodium channel blockers worsen seizures in SCN1A loss-of-function"
    },
    "prognosis": "Drug-resistant epilepsy, intellectual disability, increased SUDEP risk",
    "cross_modal_trigger": True,
    "critical_pharmacogenomic": True
}
```

**Clinical Scale: NIHSS (National Institutes of Health Stroke Scale)**

```python
{
    "name": "NIHSS",
    "full_name": "National Institutes of Health Stroke Scale",
    "clinical_use": "Quantify stroke severity, guide treatment decisions, predict outcome",
    "scoring": {
        "range": "0-42",
        "categories": {
            "0": "No stroke symptoms",
            "1-4": "Minor stroke",
            "5-15": "Moderate stroke",
            "16-20": "Moderate to severe stroke",
            "21-42": "Severe stroke"
        }
    },
    "domains_assessed": [
        "1a. Level of consciousness",
        "1b. LOC questions (month, age)",
        "1c. LOC commands (open/close eyes, grip/release)",
        "2. Best gaze (horizontal eye movement)",
        "3. Visual fields",
        "4. Facial palsy",
        "5. Motor arm (L and R)",
        "6. Motor leg (L and R)",
        "7. Limb ataxia",
        "8. Sensory",
        "9. Best language",
        "10. Dysarthria",
        "11. Extinction/inattention"
    ],
    "clinical_decision_thresholds": {
        "tpa_consideration": ">= 4 (measurable deficit)",
        "thrombectomy_consideration": ">= 6 with LVO",
        "minor_stroke_no_intervention": "0-3 (unless disabling deficit)"
    },
    "guideline_source": "AHA/ASA 2019 Acute Ischemic Stroke Guidelines"
}
```

---

## 11. Query Expansion and Retrieval Strategy

### 11.1 Neurology-Specific Query Expansion Maps

Sixteen domain-specific expansion maps:

| Map | Keywords -> Terms | Example |
|---|---|---|
| Stroke | 25 -> 200 | "brain attack" -> stroke, CVA, cerebrovascular accident, ischemic stroke, hemorrhagic stroke, TIA, thrombectomy, tPA |
| Dementia | 25 -> 180 | "memory loss" -> dementia, Alzheimer's, cognitive decline, MCI, neurodegeneration, amyloid, tau |
| Epilepsy | 20 -> 160 | "fits" -> seizure, epilepsy, convulsion, ictal, interictal, epileptiform, paroxysmal, ASM |
| MS | 20 -> 140 | "relapse" -> multiple sclerosis, demyelination, lesion, RRMS, SPMS, DMT, NEDA |
| Parkinson's | 20 -> 140 | "tremor" -> Parkinson's disease, parkinsonism, dopaminergic, bradykinesia, rigidity, levodopa |
| Brain Tumor | 15 -> 100 | "brain cancer" -> glioma, glioblastoma, meningioma, IDH, MGMT, neuro-oncology, WHO grade |
| Headache | 15 -> 100 | "bad headache" -> migraine, cluster, tension, CGRP, triptan, aura, photophobia |
| Neuromuscular | 15 -> 90 | "weak muscles" -> neuropathy, myopathy, myasthenia, ALS, GBS, CIDP, motor neuron |
| EEG | 15 -> 100 | "brain waves" -> EEG, electroencephalogram, epileptiform, spike, seizure monitoring, VEEG |
| Neuroimaging | 15 -> 100 | "brain scan" -> MRI, CT, PET, FLAIR, DWI, brain volumetrics, white matter |
| Neurogenetics | 15 -> 100 | "genetic brain disease" -> neurogenetic, channelopathy, trinucleotide repeat, gene panel |
| Movement | 10 -> 70 | "shaking" -> tremor, dystonia, chorea, ataxia, tics, myoclonus, DBS |
| Sleep Neurology | 10 -> 60 | "sleep problems" -> insomnia, narcolepsy, RBD, sleep apnea, restless legs |
| Neuroimmunology | 10 -> 80 | "brain inflammation" -> encephalitis, autoimmune, NMDA, LGI1, CASPR2, MOG, AQP4 |
| Neurorehab | 10 -> 60 | "brain recovery" -> neuroplasticity, rehabilitation, functional recovery, constraint therapy |
| CSF Analysis | 10 -> 70 | "spinal tap" -> lumbar puncture, CSF, oligoclonal bands, protein, glucose, cytology |

**Total: ~250 keywords -> ~1,750 expanded terms**

### 11.2 Comparative Analysis Detection

**Cardiology-specific comparisons:**

| Comparison | Clinical Relevance |
|---|---|
| "Lecanemab vs Donanemab" | Anti-amyloid therapy selection in early AD |
| "Carbamazepine vs Lamotrigine" | First-line ASM for focal epilepsy |
| "Natalizumab vs Ocrelizumab" | High-efficacy DMT selection in MS |
| "Levodopa vs Dopamine Agonist" | Initial PD therapy |
| "tPA vs Thrombectomy" | Acute stroke intervention strategy |
| "TMZ vs CCNU" | Chemotherapy for recurrent glioma |
| "DBS STN vs GPi" | DBS target selection in PD |
| "Erenumab vs Fremanezumab" | CGRP mAb selection for migraine prevention |

---

## 12. API and UI Design

### 12.1 FastAPI Endpoints (Port 8528)

| Method | Path | Purpose |
|---|---|---|
| GET | `/health` | Service health, collection stats, NIM status |
| GET | `/collections` | List all collections with vector counts |
| POST | `/query` | Full RAG query with evidence synthesis |
| POST | `/search` | Evidence-only search (no LLM) |
| POST | `/api/ask` | Chat-style question answering |
| POST | `/find-related` | Cross-collection entity linking |
| GET | `/workflows` | List available clinical workflows |
| POST | `/workflow/{name}/run` | Execute a clinical workflow |
| GET | `/demo-cases` | List pre-loaded demo cases |
| POST | `/demo-cases/{id}/run` | Run a demo case |
| POST | `/scale/calculate` | Calculate validated neurological scale |
| POST | `/diagnostic/algorithm` | Run diagnostic algorithm (e.g., dementia differential) |
| POST | `/protocol/recommend` | Imaging protocol recommendation |
| POST | `/reports/generate` | Generate report (markdown, JSON, PDF, FHIR) |
| GET | `/knowledge/stats` | Knowledge graph statistics |
| GET | `/metrics` | Prometheus-compatible metrics |

### 12.2 Streamlit UI (Port 8529) -- 10 Tabs

| Tab | Purpose |
|---|---|
| **Evidence Explorer** | Multi-collection RAG Q&A with citations, pre-filled neuro queries, Plotly donut chart |
| **Workflow Runner** | 8 clinical workflows with pre-loaded demo cases, AI-annotated brain images, pipeline animation |
| **Neuroimaging Gallery** | Brain MRI/CT/PET showcase with AI annotations, 3D volume viewer, before/after toggle |
| **Clinical Scales** | Interactive NIHSS, GCS, MoCA, UPDRS, EDSS, mRS, HIT-6, ALSFRS-R calculators |
| **Diagnostic Algorithms** | Step-through diagnostic pathways for dementia DDx, seizure classification, neuropathy workup |
| **Device & AI Ecosystem** | 120+ FDA-cleared neuro AI devices, searchable by modality and condition |
| **Protocol Advisor** | Patient-specific neuroimaging protocol recommendations |
| **Reports & Export** | Markdown, JSON, NVIDIA-branded PDF, FHIR R4 DiagnosticReport |
| **Patient 360** | Cross-modal neuro dashboard: imaging + EEG + genomics + scales with Plotly network graph |
| **Guidelines & Trials** | AAN/EAN/ILAE guideline browser, landmark trial summaries |

### 12.3 Demo Cases

| ID | Title | Workflow | Key Features |
|---|---|---|---|
| DEMO-001 | Acute Stroke: LVO with Thrombectomy Decision | acute_stroke_triage | NIHSS 14, ASPECTS 8, LVO M1, DAWN criteria |
| DEMO-002 | Memory Clinic: Early-Onset Alzheimer's | dementia_evaluation | MoCA 18, hippocampal atrophy, amyloid PET+, PSEN1 trigger |
| DEMO-003 | Drug-Resistant Epilepsy: Surgical Evaluation | epilepsy_focus_localization | Concordance matrix, hippocampal sclerosis, SCN1A trigger |
| DEMO-004 | New Brain Mass: Suspected Glioblastoma | brain_tumor_grading | Ring-enhancing mass, tumor volumes, IDH/MGMT prediction |
| DEMO-005 | Young Woman with Optic Neuritis: MS vs NMOSD | ms_monitoring | McDonald criteria, OCB+, NMOSD differential |

---

## 13. Clinical Decision Support Engines

### 13.1 Validated Neurological Scales

| Scale | Use Case | Range | Interpretation |
|---|---|---|---|
| **NIHSS** | Stroke severity | 0-42 | 0=no symptoms, 1-4=minor, 5-15=moderate, 16-20=mod-severe, 21+=severe |
| **GCS** | Consciousness level | 3-15 | 3-8=severe (coma), 9-12=moderate, 13-15=mild |
| **MoCA** | Cognitive screening | 0-30 | >=26 normal, <26 cognitive impairment, <18 probable dementia |
| **UPDRS Part III** | PD motor severity | 0-132 | Higher = more severe motor impairment |
| **EDSS** | MS disability | 0-10 | 0=normal, 4=walking limited, 6=unilateral assistance, 7=wheelchair, 10=death |
| **mRS** | Global disability (stroke) | 0-6 | 0=no symptoms, 2=slight disability, 4=mod-severe, 6=dead |
| **HIT-6** | Headache impact | 36-78 | <=49=little/no impact, 50-55=some, 56-59=substantial, >=60=severe |
| **ALSFRS-R** | ALS function | 0-48 | 48=normal, decline ~1 pt/month typical, faster=worse prognosis |
| **CDR** | Dementia staging | 0-3 | 0=normal, 0.5=MCI, 1=mild, 2=moderate, 3=severe dementia |
| **ASPECTS** | CT stroke scoring | 0-10 | 10=normal, >=6=good candidate for intervention, <6=large infarct |

---

## 14. Reporting and Interoperability

### 14.1 Export Formats

| Format | Use Case | Standards |
|---|---|---|
| Markdown | Clinical narrative, consultation notes | -- |
| JSON | Programmatic consumption, dashboards | -- |
| PDF | NVIDIA-themed clinical documentation | ReportLab |
| FHIR R4 | EHR integration, interoperability | SNOMED CT, LOINC, DICOM |

### 14.2 FHIR R4 Neurological Coding

| Element | Code System | Example Codes |
|---|---|---|
| Findings | SNOMED CT | 230690007 (Stroke), 26929004 (AD), 84757009 (Epilepsy), 24700007 (MS), 49049000 (PD) |
| Observations | LOINC | 72172-0 (NIHSS), 9269-2 (GCS), 72133-2 (MoCA), LP6464-2 (MRI Brain) |
| Procedures | SNOMED CT | 429064006 (Thrombectomy), 445185007 (Brain MRI), 54550000 (EEG) |
| Medications | RxNorm | Lecanemab, levetiracetam, ocrelizumab, levodopa/carbidopa, sumatriptan |

---

## 15. Product Requirements Document

### 15.1 Product Vision

**Vision Statement:** Enable any neurologist, anywhere, to access integrated neurological intelligence combining neuroimaging AI, electrophysiology correlation, genomic analysis, validated clinical scales, and evidence synthesis -- on a single $4,699 device.

### 15.2 User Stories

#### Epic 1: Evidence-Based Clinical Queries

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-001 | As a neurologist, I want to ask clinical questions and receive evidence-grounded answers with citations. | P0 | Query returns answer with >=3 citations from >=2 collections; <30 sec |
| US-002 | As a neurologist, I want comparative analysis ("Lecanemab vs Donanemab") with structured tables. | P0 | Auto-detected; side-by-side evidence; structured comparison |
| US-003 | As a fellow, I want pre-filled example queries for common neuro scenarios. | P1 | >=4 clickable examples; each returns relevant results |
| US-004 | As a researcher, I want to filter by neurological domain, imaging modality, and year. | P1 | Sidebar filters applied; results reflect filters |

#### Epic 2: Clinical Workflows

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-005 | As a stroke neurologist, I want acute stroke triage with NIHSS, ASPECTS, LVO detection, and thrombectomy eligibility. | P0 | Correct NIHSS scoring; ASPECTS calculation; DAWN/DEFUSE-3 criteria applied |
| US-006 | As a memory clinic physician, I want dementia differential diagnosis with ATN classification. | P0 | Atrophy pattern -> differential; ATN staging; genetic triggers for early-onset |
| US-007 | As an epileptologist, I want seizure classification and concordance assessment for surgical evaluation. | P0 | ILAE classification; concordance matrix; ASM recommendations by type |
| US-008 | As a neuro-oncologist, I want WHO 2021 integrated tumor classification prediction from imaging. | P0 | Tumor segmentation volumes; molecular prediction; treatment protocol |
| US-009 | As an MS specialist, I want NEDA-3 assessment with DMT escalation recommendations. | P1 | Lesion comparison; brain atrophy rate; NEDA-3 status; escalation ladder |

#### Epic 3: Clinical Scales

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-010 | As a neurologist, I want interactive NIHSS, GCS, and MoCA calculators with interpretation. | P0 | All items scored; total with severity category; guideline recommendation |
| US-011 | As an MS specialist, I want EDSS calculation with progression tracking over time. | P1 | EDSS 0-10; functional system scores; progression confirmation criteria |
| US-012 | As a movement disorder specialist, I want UPDRS scoring with LEDD calculation. | P1 | UPDRS Part III; H&Y staging; LEDD from medication list |

#### Epic 4: Cross-Modal Integration

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-013 | As a neurogenetics specialist, I want automatic genomic triggers from neuroimaging findings. | P0 | Qualifying imaging pattern triggers gene panel query; genomic hits displayed |
| US-014 | As a neurologist, I want Patient 360 combining imaging, EEG, genomics, and scales. | P1 | Interactive network graph; cross-modal connections; drill-down to evidence |

#### Epic 5: Reporting and Export

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-015 | As a neurologist, I want PDF clinical reports with neurological findings, scales, and evidence. | P0 | NVIDIA-branded PDF; all sections populated; download button |
| US-016 | As a health IT engineer, I want FHIR R4 DiagnosticReport with neuro SNOMED/LOINC codes. | P1 | Valid FHIR R4; neuro-specific coding; passes validator |

#### Epic 6: Demo and Presentation

| ID | User Story | Priority | Acceptance Criteria |
|---|---|---|---|
| US-017 | As a demo presenter, I want 5 pre-loaded demo cases covering the breadth of neurology. | P0 | 5 cases executable; <30 seconds each; realistic output |
| US-018 | As a new user, I want a sidebar guided tour for the 10-tab interface. | P1 | Expandable tour; numbered steps; dismiss button |

### 15.3 Non-Functional Requirements

| Requirement | Target |
|---|---|
| RAG query latency | < 30 seconds end-to-end |
| Scale calculator latency | < 5 seconds |
| Workflow execution (mock) | < 10 seconds |
| Memory footprint | < 16 GB (agent only) |
| Seed data completeness | 2,270+ records across 13 collections |
| Unit test coverage | > 80% |
| FHIR R4 compliance | Passes HL7 FHIR Validator |

### 15.4 Prioritization Matrix

| Phase | Features | Timeline |
|---|---|---|
| **Phase 1 (MVP)** | RAG engine (13 collections), Evidence Explorer, 3 workflows (stroke, dementia, epilepsy), 4 scales (NIHSS, GCS, MoCA, mRS), 3 demo cases, PDF export | 5-7 weeks |
| **Phase 2 (Complete)** | All 8 workflows, all 10 scales, diagnostic algorithms, FHIR R4, neuroimaging gallery, Patient 360, all 5 demo cases | 5-7 weeks |
| **Phase 3 (Polish)** | Guided tour, pipeline animation, cross-modal triggers, network graph, Guidelines & Trials tab, EEG gallery | 2-3 weeks |

---

## 16. Data Acquisition Strategy

### 16.1 Automated Ingest Pipelines

| Source | Collection(s) | Method | Update Cadence |
|---|---|---|---|
| PubMed (NCBI E-utilities) | neuro_literature | MeSH-filtered neuroscience retrieval | Weekly |
| ClinicalTrials.gov (V2 API) | neuro_trials | Neurological condition filter | Weekly |
| AAN Practice Guidelines | neuro_guidelines | Manual curation + embedding | Per guideline update |

### 16.2 PubMed Search Strategy

```
MeSH Terms:
  "Nervous System Diseases"[MeSH] OR "Neurodegenerative Diseases"[MeSH] OR
  "Epilepsy"[MeSH] OR "Stroke"[MeSH] OR "Multiple Sclerosis"[MeSH] OR
  "Brain Neoplasms"[MeSH] OR "Parkinson Disease"[MeSH] OR "Alzheimer Disease"[MeSH]

AND ("Artificial Intelligence"[MeSH] OR "Machine Learning"[MeSH] OR
     "Deep Learning"[MeSH] OR "Neural Networks"[MeSH] OR
     "Neuroimaging"[MeSH])

Filters: Published 2018-2026, English, Humans
Expected: 3,500-6,000 abstracts
```

---

## 17. Validation and Testing Strategy

### 17.1 Unit Tests

| Test Category | Target Count |
|---|---|
| Collection schemas | 39 (13 x 3) |
| Clinical scales | 50 (10 scales x 5 cases) |
| Workflow logic | 40 (8 x 5) |
| RAG engine | 20 |
| Knowledge graph | 15 |
| API endpoints | 30 |
| FHIR R4 export | 15 |
| Cross-modal triggers | 15 |
| Diagnostic algorithms | 20 |
| **Total** | **~245** |

### 17.2 Clinical Validation

| Validation Type | Method | Target |
|---|---|---|
| Scale accuracy | Compare against published validation data | < 1% deviation |
| Diagnostic algorithms | Board-certified neurologist review | 95%+ guideline concordance |
| Seizure classification | Compare against ILAE criteria | 100% match for standard inputs |
| FHIR R4 compliance | HL7 FHIR Validator | Zero validation errors |

---

## 18. Regulatory Considerations

### 18.1 FDA CDS Exemption (21st Century Cures Act)

| Criterion | Assessment |
|---|---|
| Not intended to acquire, process, or analyze a medical image or signal | **Met** -- RAG and scale calculation; does not process raw images/EEG |
| Intended for displaying, analyzing, or printing medical information | **Met** |
| Intended for use by healthcare professional | **Met** -- Designed for neurologists |
| Healthcare professional does not primarily rely on the software | **Met** -- Provides recommendations for review |

**Assessment:** Core RAG and scale calculator functions likely qualify for CDS exemption. NIM-based image analysis workflows would require separate regulatory consideration.

### 18.2 Disclaimers

> *This tool is for clinical decision support only and does not replace professional medical judgment. All findings require verification by a qualified neurologist. Scale calculations are estimates and must be confirmed by clinical examination. Not FDA-cleared for autonomous clinical decision-making. For research and educational purposes only.*

---

## 19. DGX Compute Progression

| Phase | Hardware | Price | Scope |
|---|---|---|---|
| **Phase 1 -- Proof Build** | DGX Spark | $4,699 | All 8 workflows (mock/cloud NIM), 13 collections, 10 clinical scales |
| **Phase 2 -- Departmental** | 1-2x DGX B200 | $500K-$1M | Full NIM stack, live MRI/CT/EEG processing, PACS integration |
| **Phase 3 -- Multi-Site** | 4-8x DGX B200 | $2M-$4M | NVIDIA FLARE federated learning across neurology centers |
| **Phase 4 -- AI Factory** | DGX SuperPOD | $7M-$60M+ | National brain health surveillance, real-time ICU neuro monitoring |

---

## 20. Implementation Roadmap

### Phase 1: Foundation (Weeks 1-7)

| Week | Deliverable |
|---|---|
| 1-2 | Repository scaffolding, Pydantic models, settings, Docker Compose, 13 collection schemas |
| 3-4 | Seed data curation (2,270 records), ingest pipelines (PubMed, ClinicalTrials.gov) |
| 5-7 | RAG engine, Evidence Explorer tab, 3 priority workflows (stroke, dementia, epilepsy), 4 clinical scales |

### Phase 2: Clinical Intelligence (Weeks 8-14)

| Week | Deliverable |
|---|---|
| 8-9 | Remaining 5 workflows (brain tumor, MS, Parkinson's, headache, neuromuscular) |
| 10-11 | Remaining 6 clinical scales, diagnostic algorithm engine |
| 12-14 | Cross-modal triggers, FHIR R4 export, PDF reports |

### Phase 3: UI and Polish (Weeks 15-18)

| Week | Deliverable |
|---|---|
| 15-16 | 10-tab Streamlit UI, neuroimaging gallery, Patient 360 network graph |
| 17-18 | Guided tour, demo cases finalized, documentation, integration testing |

---

## 21. Risk Analysis

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Clinical scale implementation errors | Medium | High | Validate against published reference implementations |
| Diagnostic algorithm complexity (dementia DDx) | Medium | Medium | Start with most common differentials; expand iteratively |
| Insufficient EEG/EMG seed data | Medium | Medium | Focus on patterns with highest clinical impact |
| NIM brain segmentation accuracy for volumetrics | Low | Medium | Age/sex-normalized percentiles with published norms |
| Guideline updates during development | Low | Low | Modular guideline collection |
| Memory pressure with 13 collections + other agents | Low | Medium | BGE-small (384-dim) is compact |
| Neurogenetics variant interpretation complexity | Medium | Medium | Leverage ClinVar pathogenicity + AlphaMissense scores from shared collection |

---

## 22. Competitive Landscape

### 22.1 Positioning

```
                    Multi-Condition Coverage
                           ^
                           |
                    [Neurology Agent]
                           |
         Cloud-only <------+------> On-device
                           |
                    [Viz.ai, Aidoc]
                           |
                           v
                    Single-Condition
```

No existing product combines:
1. Multi-condition neurological coverage (stroke + dementia + epilepsy + tumors + MS + PD + headache + NMD)
2. Cross-modal integration (imaging + EEG + genomics)
3. RAG-grounded evidence synthesis with citations
4. Validated neurological scales (NIHSS, GCS, MoCA, UPDRS, EDSS)
5. Diagnostic algorithms
6. On-device deployment ($4,699)
7. Open-source (Apache 2.0)

### 22.2 Defensibility

| Advantage | Defensibility |
|---|---|
| Multi-collection RAG architecture | High -- 6 agents prove the pattern |
| Cross-modal genomic triggers | High -- unique to HCLS AI Factory |
| On-device inference | High -- DGX Spark + NIM is NVIDIA-exclusive |
| Breadth of neurological coverage | High -- no competitor covers 8 subspecialties |
| Open source | Medium -- community adoption |
| Clinical validation | Medium -- requires domain expertise |

---

## 23. Discussion

### 23.1 Why Neurology Is a High-Impact Agent

1. **Workforce crisis**: With a projected shortfall of 11,000 neurologists by 2031, AI-assisted decision support is not optional -- it is a clinical necessity. A community neurologist covering stroke, epilepsy, MS, dementia, and headache cannot maintain subspecialty-level expertise across all domains.

2. **Diagnostic complexity**: Many neurological conditions require multi-modal data synthesis (imaging + EEG + genomics + clinical scales) that exceeds the cognitive capacity of a single clinician reviewing siloed data. AI excels at this pattern.

3. **Genomic transformation**: Over 1,000 neurological genes with clinical actionability. SCN1A changes epilepsy treatment. GBA changes Parkinson's prognosis. APOE/PSEN1 changes Alzheimer's management. A genomic-integrated intelligence agent enables precision neurology.

4. **Therapeutic pipeline explosion**: Anti-amyloid antibodies, gene therapies, antisense oligonucleotides, and BTK inhibitors are transforming previously untreatable conditions. Selecting the right therapy for the right patient at the right time requires evidence synthesis beyond PubMed searches.

5. **Existing foundation**: The HCLS AI Factory already has stroke triage (DEMO-001 in Imaging Agent), brain MRI MS lesion tracking, hippocampal volumetry via VISTA-3D, and neurological genes in the genomic_evidence collection.

6. **Differentiation**: No competitor addresses neurology's full breadth. Viz.ai does stroke. Persyst does EEG. QMENTA does imaging. None integrate all three with genomics, scales, and literature.

### 23.2 Limitations

1. **Mock mode inference**: Phase 1 uses synthetic results; GPU deployment required for clinical validation.
2. **EEG signal processing**: Raw EEG analysis requires specialized preprocessing not covered by standard NIM services. Initial version will work with pre-interpreted EEG findings.
3. **Scale automation**: Some scales (NIHSS, UPDRS) require bedside examination and cannot be fully automated from data alone.
4. **Rare disease coverage**: With >7,000 rare neurological conditions, complete coverage is infeasible. Focus is on the 50 most common/impactful conditions.

### 23.3 Future Directions

1. **Real-time ICU neuromonitoring**: Continuous EEG seizure detection with automated NIHSS trending
2. **Brain-computer interface integration**: Neurofeedback and BCI data for rehabilitation
3. **Longitudinal brain health platform**: Population-level brain aging trajectories
4. **Digital biomarkers**: Smartphone-based gait analysis, speech markers for cognitive decline
5. **Federated learning across epilepsy centers**: NVIDIA FLARE for multi-site seizure detection models

---

## 24. Conclusion

The Neurology Intelligence Agent addresses a critical unmet need at the intersection of the global neurological disease burden, the neurologist workforce shortage, and the fragmentation of neurological data. By leveraging the HCLS AI Factory's proven multi-collection RAG architecture -- adapted with 13 neurology-specific collections, 8 clinical workflows, 10 validated scales, and cross-modal neuroimaging-genomics triggers -- the agent will deliver subspecialty-level neurological intelligence to any neurologist, anywhere.

The breadth of coverage (stroke, dementia, epilepsy, brain tumors, MS, Parkinson's, headache, neuromuscular disease) is unmatched by any existing commercial product, all of which address only one or two subspecialties. Combined with genomic integration that enables precision neurology (SCN1A-guided epilepsy treatment, GBA-stratified Parkinson's management, PSEN1/2 Alzheimer's diagnosis), the agent represents a new paradigm in neurological clinical decision support.

Deploying on a single NVIDIA DGX Spark ($4,699) ensures accessibility for community neurologists, rural health systems, and resource-limited institutions globally -- precisely the settings where the neurologist shortage is most acute and subspecialty expertise least available.

---

## 25. References

1. GBD 2021 Nervous System Disorders Collaborators. Global, regional, and national burden of disorders affecting the nervous system, 1990-2021. *Lancet Neurol*. 2024;23:344-381.
2. World Health Organization. *Intersectoral Global Action Plan on Epilepsy and Other Neurological Disorders 2022-2031*. WHO, 2023.
3. Feigin VL, et al. The global burden of neurological disorders. *Lancet Neurol*. 2024;23:2-3.
4. Dorsey ER, et al. Global, regional, and national burden of Parkinson's disease, 1990-2016. *Lancet Neurol*. 2018;17:939-953.
5. 2023 Alzheimer's Disease Facts and Figures. *Alzheimers Dement*. 2023;19:1598-1695.
6. GBD 2016 Stroke Collaborators. Global, regional, and national burden of stroke, 1990-2016. *Lancet Neurol*. 2019;18:439-458.
7. Powers WJ, et al. Guidelines for the Early Management of Patients with Acute Ischemic Stroke: 2019 Update. *Stroke*. 2019;50:e344-e418.
8. Nogueira RG, et al. Thrombectomy 6 to 24 Hours after Stroke with a Mismatch between Deficit and Infarct (DAWN). *N Engl J Med*. 2018;378:11-21.
9. Albers GW, et al. Thrombectomy for Stroke at 6 to 16 Hours with Selection by Perfusion Imaging (DEFUSE-3). *N Engl J Med*. 2018;378:708-718.
10. Jack CR, et al. NIA-AA Research Framework: Toward a biological definition of Alzheimer's disease. *Alzheimers Dement*. 2018;14:535-562.
11. van Dyck CH, et al. Lecanemab in Early Alzheimer's Disease (CLARITY AD). *N Engl J Med*. 2023;388:9-21.
12. Sims JR, et al. Donanemab in Early Symptomatic Alzheimer Disease (TRAILBLAZER-ALZ 2). *JAMA*. 2023;330:512-527.
13. Fisher RS, et al. ILAE Official Report: A practical clinical definition of epilepsy. *Epilepsia*. 2014;55:475-482.
14. Scheffer IE, et al. ILAE Classification of the Epilepsies: Position paper. *Epilepsia*. 2017;58:512-521.
15. Thompson AJ, et al. Diagnosis of multiple sclerosis: 2017 revisions of the McDonald criteria. *Lancet Neurol*. 2018;17:162-173.
16. Rae-Grant A, et al. Practice guideline recommendations summary: Disease-modifying therapies for adults with multiple sclerosis (AAN). *Neurology*. 2018;90:777-788.
17. Postuma RB, et al. MDS clinical diagnostic criteria for Parkinson's disease. *Mov Disord*. 2015;30:1591-1601.
18. Headache Classification Committee of the IHS. The International Classification of Headache Disorders, 3rd edition (ICHD-3). *Cephalalgia*. 2018;38:1-211.
19. Louis DN, et al. The 2021 WHO Classification of Tumors of the Central Nervous System. *Neuro Oncol*. 2021;23:1231-1251.
20. Hemphill JC, et al. Guidelines for the Management of Spontaneous Intracerebral Hemorrhage. *Stroke*. 2015;46:2032-2060.
21. FDA. Artificial Intelligence and Machine Learning (AI/ML)-Enabled Medical Devices. FDA Database, 2025.
22. American Academy of Neurology. *2023 Neurology Workforce Study*. AAN, 2023.
23. Wirrell EC, et al. Dravet syndrome: Diagnosis, treatment, and long-term management. *Epilepsia*. 2022;63:S131-S147.
24. Hogan DB, et al. Nusinersen and the treatment of spinal muscular atrophy. *CMAJ*. 2018;190:E1-E3.
25. Grossman M, et al. Frontotemporal lobar degeneration. *Nat Rev Dis Primers*. 2023;9:40.
26. Filippi M, et al. Assessment of lesions on magnetic resonance imaging in multiple sclerosis. *Brain*. 2019;142:1858-1875.

---

*HCLS AI Factory -- Neurology Intelligence Agent Research Paper and PRD*
*Apache 2.0 License | March 2026*
*Author: Adam Jones*
