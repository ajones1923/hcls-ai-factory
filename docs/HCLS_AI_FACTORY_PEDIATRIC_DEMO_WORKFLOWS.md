# HCLS AI Factory: Pediatric Oncology Demo Workflows

## A Comprehensive Guide to Demonstrating Precision Medicine for Children with Cancer

---

**Author:** Adam Jones

**Date:** March 2026

**Version:** 1.0

**Classification:** Technical Reference / Demo Guide

---

**Document Information**

| Field | Value |
|---|---|
| Document Title | HCLS AI Factory: Pediatric Oncology Demo Workflows |
| Subtitle | A Comprehensive Guide to Demonstrating Precision Medicine for Children with Cancer |
| Author | Adam Jones |
| Organization | HCLS AI Factory Project |
| Publication Date | March 2026 |
| Version | 1.0 |
| Platform | NVIDIA DGX Spark |
| Target Audience | Clinical researchers, oncologists, demonstration engineers, healthcare IT professionals |

---

**Abstract**

This document presents six comprehensive demonstration workflows for the HCLS AI Factory
platform, each specifically designed around pediatric oncology use cases. The workflows span
the full precision medicine pipeline--from raw genomic sequencing data (FASTQ) through variant
calling, multi-agent clinical intelligence, and therapeutic candidate discovery. Each workflow
is designed with dual entry points: a FASTQ path for audiences interested in the complete
genomic processing pipeline, and a VCF path for audiences focused on clinical intelligence
and drug discovery capabilities. The six demos collectively exercise every major component
of the platform: the Genomic Foundation Engine (Parabricks, DeepVariant, BWA-MEM2), the
Precision Intelligence Engine (seven specialized AI agents with 80+ Milvus vector
collections), and the Therapeutic Discovery Engine (BioNeMo MolMIM, DiffDock, RDKit).
Together, they tell the story of how an integrated AI platform can compress the timeline
from a child's diagnosis to an actionable treatment plan from weeks or months to hours.

---

**Table of Contents**

1. [Title Page](#hcls-ai-factory-pediatric-oncology-demo-workflows)
2. [Executive Summary](#executive-summary)
3. [Introduction: Why Pediatric Oncology](#introduction-why-pediatric-oncology)
    - 3.1 [The Scale of Pediatric Cancer](#31-the-scale-of-pediatric-cancer)
    - 3.2 [Progress and Plateaus](#32-progress-and-plateaus)
    - 3.3 [The Hardest Cancers in Childhood](#33-the-hardest-cancers-in-childhood)
    - 3.4 [Current Clinical Workflow Challenges](#34-current-clinical-workflow-challenges)
    - 3.5 [How the HCLS AI Factory Addresses Each Challenge](#35-how-the-hcls-ai-factory-addresses-each-challenge)
4. [Platform Overview](#platform-overview)
    - 4.1 [Three-Engine Architecture](#41-three-engine-architecture)
    - 4.2 [The Agent Constellation](#42-the-agent-constellation)
    - 4.3 [Cross-Agent Integration](#43-cross-agent-integration)
    - 4.4 [Shared Genomic Evidence Layer](#44-shared-genomic-evidence-layer)
    - 4.5 [Dual Entry Points](#45-dual-entry-points)
    - 4.6 [Service Port Map](#46-service-port-map)
5. [Demo 1: End-to-End Precision Medicine Pipeline](#section-5-demo-1----end-to-end-precision-medicine-pipeline-the-foundation-demo)
6. [Demo 2: Pediatric ALL Multi-Agent Tumor Board](#section-6-demo-2----pediatric-all-multi-agent-tumor-board)
7. [Demo 3: Cardiotoxicity Prevention in Pediatric Chemotherapy](#section-7-demo-3----cardiotoxicity-prevention-in-pediatric-chemotherapy)
8. [Demo 4: Rare Disease with Cancer Predisposition](#section-8-demo-4----rare-disease-with-cancer-predisposition)
9. [Demo 5: CAR-T Therapy Decision and Monitoring](#section-9-demo-5----car-t-therapy-decision-and-monitoring)
10. [Demo 6: Medulloblastoma Precision Treatment with Novel Drug Discovery](#section-10-demo-6----medulloblastoma-precision-treatment-with-novel-drug-discovery)
11. [Demo Execution Guide](#section-11-demo-execution-guide)
12. [Clinical Validation References](#section-12-clinical-validation-references)
13. [Appendix — Sample API Payloads](#section-13-appendix--sample-api-payloads)

---

## Executive Summary

Pediatric oncology stands at an inflection point. While five-year survival rates for childhood
cancers have risen dramatically over the past six decades--from roughly 10% in the 1960s to
approximately 85% today--progress has plateaued for many of the most aggressive malignancies.
Diffuse midline gliomas carry a five-year survival rate below 10%. Relapsed acute
lymphoblastic leukemia (ALL), despite initial treatment success rates exceeding 90%, sees
survival drop to 30-40% upon relapse. High-risk neuroblastoma, metastatic Ewing sarcoma,
and relapsed medulloblastoma remain among the most challenging diagnoses in all of medicine.
For these children, the traditional clinical workflow--manual tumor board review, fragmented
data systems, and therapies derived primarily from adult trials--is insufficient.

The HCLS AI Factory is a precision medicine platform built on the NVIDIA DGX Spark that
compresses the journey from a patient's raw DNA sequencing data to actionable drug candidates
from weeks or months to under five hours. This document describes six demonstration workflows
that showcase the platform's capabilities through the lens of pediatric oncology, the domain
where speed, precision, and integration matter most.

### The Six Demo Workflows

Each demonstration workflow targets a specific pediatric cancer type and exercises a distinct
combination of platform capabilities:

1. **Pediatric Leukemia Molecular Tumor Board** -- Demonstrates the full pipeline from FASTQ
   input through variant calling to a structured molecular tumor board report. This workflow
   exercises the Genomic Foundation Engine (Parabricks alignment, DeepVariant variant calling),
   the Precision Oncology Agent (CIViC/OncoKB annotation, therapy ranking), and the Precision
   Biomarker Agent (pharmacogenomic profiling). The demo produces a tiered evidence report
   suitable for clinical review, with matched clinical trials and FDA-approved therapy
   recommendations specific to pediatric ALL.

2. **Neuroblastoma Risk Stratification** -- Focuses on multi-agent integration for complex
   risk assessment. Neuroblastoma is uniquely heterogeneous: outcomes range from spontaneous
   regression to rapid metastatic progression. This workflow demonstrates how the Precision
   Oncology Agent, Biomarker Agent, and Imaging Intelligence Agent coordinate through
   cross-agent triggers to produce an integrated risk stratification that incorporates
   MYCN amplification status, segmental chromosomal aberrations, tumor imaging features,
   and age-adjusted prognostic factors.

3. **Diffuse Midline Glioma Therapeutic Discovery** -- Showcases the Therapeutic Discovery
   Engine for the hardest pediatric brain tumors. Starting from a VCF with the characteristic
   H3K27M mutation, this workflow moves through variant interpretation to novel molecule
   generation using BioNeMo MolMIM, protein-ligand docking via DiffDock, and ADMET
   property scoring via RDKit. The demo produces a ranked list of candidate compounds
   targeting the altered histone modification pathway.

4. **CAR-T Therapy Design for Relapsed ALL** -- Highlights the CAR-T Intelligence Agent's
   specialized knowledge base of 6,266+ vectors covering CAR construct designs, manufacturing
   protocols, clinical trial outcomes, and toxicity management. This workflow demonstrates
   comparative analysis of costimulatory domain choices (4-1BB vs. CD28), target antigen
   selection (CD19, CD22, dual-targeting), and manufacturing optimization for pediatric
   patients.

5. **Pharmacogenomic-Guided Chemotherapy Dosing** -- Demonstrates the Pharmacogenomics
   Agent's ability to translate germline variants into actionable dosing recommendations.
   Pediatric patients are especially vulnerable to chemotherapy toxicity; this workflow
   shows how TPMT/NUDT15 variants guide mercaptopurine dosing in ALL maintenance therapy,
   DPYD variants inform fluoropyrimidine dosing, and UGT1A1 variants affect irinotecan
   metabolism--all critical for preventing life-threatening toxicity in growing children.

6. **Integrated Multi-Agent Tumor Board** -- The capstone demonstration bringing all agents
   together through the /integrated-assessment endpoint. A single patient case flows
   through genomic analysis, oncology interpretation, biomarker profiling, pharmacogenomic
   screening, imaging correlation, and therapeutic discovery in a unified workflow,
   producing a comprehensive report that mirrors--and accelerates--the output of a
   multidisciplinary tumor board.

### Dual Entry Points

Every workflow supports two entry points to accommodate different audiences and time
constraints:

- **FASTQ Entry Point** (Full Pipeline, approximately 5 hours): Begins with raw sequencing
  reads. Exercises the complete three-engine pipeline. Best for demonstrating end-to-end
  capability, genomics audiences, and full platform validation.

- **VCF Entry Point** (Intelligence + Discovery, approximately 30 minutes): Begins with
  pre-called variants. Skips the computationally intensive alignment and variant calling
  steps. Best for clinical audiences, time-constrained demonstrations, and iterative
  workflow development.

### Platform at a Glance

The HCLS AI Factory platform runs entirely on a single NVIDIA DGX Spark workstation with
an approximate 1.1 TB footprint. It comprises three major engines, seven operational
intelligence agents (with two additional agents in design), over 80 Milvus vector
collections containing 3.56 million annotated variant embeddings, and integration with
4.1 million ClinVar records and 71 million AlphaMissense predictions. The platform is
fully open-source, containerized via Docker Compose, and orchestrated through Nextflow
DSL2 with a Python fallback orchestrator.

This document provides step-by-step instructions for each demo, including prerequisites,
configuration, execution commands, expected outputs, troubleshooting guidance, and
talking points for presenters. It is designed to be used by demonstration engineers,
clinical researchers evaluating the platform, and healthcare IT professionals assessing
deployment feasibility.

---

## Introduction: Why Pediatric Oncology

### 3.1 The Scale of Pediatric Cancer

Cancer is the leading cause of death by disease among children past infancy in the United
States. According to the American Cancer Society's 2025 estimates, approximately 15,590
children and adolescents aged 0 to 19 will be diagnosed with cancer in the US annually,
and roughly 1,780 will die from the disease. Globally, the World Health Organization
estimates that approximately 400,000 children and adolescents aged 0 to 19 develop cancer
each year, with approximately 300,000 new diagnoses in the 0 to 14 age group alone. In
low- and middle-income countries, where access to specialized pediatric oncology care is
limited, survival rates can be as low as 20-30%, compared to the 85% five-year survival
seen in high-income nations.

The spectrum of childhood cancers differs fundamentally from adult malignancies. While
adult cancers are dominated by carcinomas arising from accumulated somatic mutations over
decades of environmental exposure, pediatric cancers are predominantly embryonal and
hematologic in origin. The most common pediatric cancers include:

- **Leukemias** (28% of pediatric cancers): Acute lymphoblastic leukemia (ALL) accounts
  for approximately 75% of pediatric leukemia cases, with acute myeloid leukemia (AML)
  comprising most of the remainder. ALL is the single most common childhood malignancy,
  with approximately 3,000 new cases per year in the US.

- **Central nervous system tumors** (26%): Including medulloblastoma, diffuse midline
  glioma (formerly DIPG), ependymoma, low-grade glioma, and atypical teratoid/rhabdoid
  tumor (AT/RT). Brain tumors are the leading cause of cancer-related death in children.

- **Neuroblastoma** (6%): An embryonal tumor of the sympathetic nervous system, almost
  exclusively a disease of early childhood. Approximately 700-800 new cases are diagnosed
  annually in the US. Its clinical behavior is remarkably heterogeneous--some tumors
  spontaneously regress, while others are rapidly fatal despite aggressive therapy.

- **Wilms tumor (nephroblastoma)** (5%): The most common renal malignancy of childhood,
  with approximately 500 new cases per year in the US. Overall survival exceeds 90%,
  but certain histologic subtypes (diffuse anaplasia) and relapsed disease carry
  significantly worse prognoses.

- **Lymphomas** (10%): Both Hodgkin lymphoma and non-Hodgkin lymphoma, including Burkitt
  lymphoma, diffuse large B-cell lymphoma, and anaplastic large cell lymphoma.

- **Soft tissue sarcomas** (7%): Rhabdomyosarcoma is the most common, with approximately
  350 new cases per year in the US.

- **Bone tumors** (5%): Osteosarcoma and Ewing sarcoma, primarily affecting adolescents.

These cancers are driven by fundamentally different molecular mechanisms than their adult
counterparts. Pediatric cancers typically have far fewer somatic mutations--often just a
handful of driver events, compared to the hundreds or thousands seen in adult carcinomas.
Instead, they are frequently driven by chromosomal translocations (e.g., ETV6-RUNX1 in
ALL, EWSR1-FLI1 in Ewing sarcoma), epigenetic dysregulation (e.g., H3K27M mutations in
diffuse midline glioma), and developmental signaling pathway disruption (e.g., MYCN
amplification in neuroblastoma, aberrant Hedgehog signaling in medulloblastoma). This
molecular distinctiveness means that adult-derived therapeutic approaches and clinical
decision-support tools are often poorly suited to pediatric cases.

### 3.2 Progress and Plateaus

The history of pediatric oncology is one of medicine's great success stories--and one of
its most frustrating plateaus. In the 1960s, a diagnosis of childhood cancer was
essentially a death sentence, with overall five-year survival rates of approximately 10%.
The introduction of combination chemotherapy regimens in the 1960s and 1970s, pioneered
by researchers including Sidney Farber, James Holland, and Donald Pinkel, transformed
outcomes for childhood leukemia. Pinkel's "Total Therapy" protocols at St. Jude Children's
Research Hospital demonstrated that ALL could be cured, not merely palliated.

By the 1980s, five-year survival for childhood ALL had risen above 70%. The development
of risk-stratified treatment protocols--adjusting therapy intensity based on age, white
blood cell count, cytogenetics, and minimal residual disease (MRD)--pushed survival rates
progressively higher. The Children's Oncology Group (COG) and its international
counterparts refined protocols through successive clinical trials enrolling tens of
thousands of children.

Today, the overall five-year survival rate for all childhood cancers combined stands at
approximately 85% in high-income countries, according to the National Cancer Institute's
Surveillance, Epidemiology, and End Results (SEER) program. For standard-risk ALL, five-year
event-free survival now exceeds 95%. Hodgkin lymphoma in children has a five-year survival
rate above 98%. Standard-risk Wilms tumor achieves five-year survival above 90%.

However, this aggregate success masks a critical reality: **survival rates have plateaued
for the cancers that kill the most children**. The improvements of the 1970s through 2000s
were largely achieved through intensification of cytotoxic chemotherapy, and we have reached
the limits of what further intensification can achieve without unacceptable toxicity. For
the past 15 to 20 years, overall survival curves for the most aggressive pediatric cancers
have been largely flat:

- **Diffuse midline glioma (DMG)**: Median survival remains approximately 9 to 11 months
  from diagnosis. Five-year overall survival is less than 10%. More than 300 clinical
  trials over four decades have failed to improve outcomes beyond supportive care and
  radiation therapy. The blood-brain barrier, diffuse infiltrative growth pattern, and
  the tumor's location in critical brainstem structures make surgical resection impossible
  and drug delivery extraordinarily challenging.

- **Metastatic Ewing sarcoma**: Five-year survival remains approximately 20-30%, a figure
  that has been largely unchanged for two decades despite intensification of multimodal
  therapy (surgery, radiation, and aggressive combination chemotherapy).

- **High-risk neuroblastoma**: Despite the introduction of anti-GD2 immunotherapy
  (dinutuximab), which improved event-free survival from approximately 46% to 66% in
  the landmark COG ANBL0032 trial, long-term survival for MYCN-amplified, stage 4
  neuroblastoma diagnosed after 18 months of age remains in the range of 40-50%.

- **Infant ALL with MLL (KMT2A) rearrangements**: Five-year event-free survival remains
  approximately 35-45%, markedly worse than other ALL subtypes. These leukemias have a
  distinct biology driven by epigenetic dysregulation rather than the typical ALL
  mechanisms.

- **Relapsed medulloblastoma**: Outcomes after relapse are dismal, with five-year
  post-relapse survival below 25% for Group 3 and Group 4 tumors. Second-line therapies
  are limited, and re-irradiation carries substantial neurocognitive risks.

The plateau is not for lack of effort. The pediatric oncology community has been
extraordinarily collaborative--more than 90% of children with cancer in the US are treated
at COG-affiliated institutions, and clinical trial enrollment rates far exceed those in
adult oncology (approximately 60% vs. less than 5%). The challenge is biological: these
are fundamentally different diseases than the ones we have learned to cure, and they
require fundamentally different approaches.

### 3.3 The Hardest Cancers in Childhood

To understand why the HCLS AI Factory focuses its demonstration workflows on pediatric
oncology, it is helpful to examine in more detail the specific cancers that remain most
challenging. These are the diseases where precision medicine and AI-accelerated workflows
have the greatest potential to change outcomes.

**Diffuse Midline Glioma (DMG)**, previously classified as diffuse intrinsic pontine
glioma (DIPG), is perhaps the most devastating diagnosis in all of pediatric oncology.
Approximately 200 to 300 children are diagnosed with DMG annually in the US. The median
age at diagnosis is 6 to 7 years. The disease is defined by the presence of histone
H3K27M mutations (found in approximately 80% of cases), which drive a global loss of
the repressive histone mark H3K27me3 and widespread epigenetic reprogramming. The
discovery of this mechanism by Sturm, Wu, and colleagues in 2012 opened new therapeutic
avenues--including ONC201 (dordaviprone), a dopamine receptor D2 antagonist that has
shown unprecedented responses in H3K27M-mutant tumors in recent trials--but the overall
prognosis remains grim. DMG exemplifies the need for rapid molecular characterization
and accelerated therapeutic discovery.

**Relapsed Acute Lymphoblastic Leukemia** is the fifth most common childhood cancer
when counted as a separate entity. Approximately 15-20% of children with ALL will relapse,
and for these patients, overall survival drops from >90% to 30-60% depending on the
timing and site of relapse. Early bone marrow relapse (within 18 months of diagnosis)
carries the worst prognosis, with five-year survival below 30%. The introduction of
CAR-T cell therapy (tisagenlecleucel, approved by the FDA in 2017 for relapsed/refractory
B-cell ALL in patients up to 25 years of age) has provided a new treatment option, but
durable remissions are achieved in only approximately 40-60% of treated patients.
Understanding resistance mechanisms (antigen loss, T-cell exhaustion, the tumor
microenvironment) requires integrated multi-agent analysis that spans genomics,
immunology, and therapeutic design.

**High-Risk Neuroblastoma** accounts for approximately 15% of pediatric cancer deaths
despite comprising only 6% of pediatric cancers. The disease is biologically extraordinary
in its heterogeneity: tumors with favorable biology (hyperdiploid, no MYCN amplification,
localized) may spontaneously regress without treatment, while tumors with unfavorable
biology (MYCN-amplified, segmental chromosomal aberrations, metastatic) are among the
most treatment-resistant cancers known. The current standard of care for high-risk
neuroblastoma involves induction chemotherapy, surgical resection, high-dose
chemotherapy with autologous stem cell rescue, radiation therapy, and maintenance
immunotherapy--a treatment course spanning 12 to 18 months with substantial acute and
long-term toxicity. Despite this intensity, relapse occurs in approximately 50% of
high-risk patients.

### 3.4 Current Clinical Workflow Challenges

The clinical workflows currently used to manage pediatric oncology cases face four
systemic challenges that limit their ability to deliver precision medicine at the speed
and scale these children need.

**Challenge 1: Manual Tumor Board Review Timelines.** The molecular tumor board (MTB)
is the standard mechanism by which genomic and molecular data are translated into
treatment recommendations. In current practice, a typical MTB workflow proceeds as
follows: tumor or blood sample collection (day 0), DNA extraction and library preparation
(days 1-3), sequencing (days 3-7), bioinformatic analysis and variant calling (days
7-14), manual variant annotation and interpretation (days 14-21), MTB meeting presentation
and discussion (days 21-28), and clinical action (days 28+). For a child with an
aggressive cancer, a four-week turnaround from biopsy to actionable recommendation is
an eternity. Studies published in Nature Medicine and the Journal of Clinical Oncology
have documented that delays in initiating molecularly targeted therapy are associated
with disease progression and loss of treatment windows in pediatric oncology.

**Challenge 2: Fragmented Data Systems.** The information needed to make a precision
treatment decision for a child with cancer is scattered across multiple disconnected
systems: the electronic health record (EHR) for clinical data, the laboratory information
system (LIS) for genomic results, clinical trial registries (ClinicalTrials.gov) for
trial matching, drug databases (DrugBank, ChEMBL) for therapeutic options, the published
literature (PubMed, with over 36 million citations) for evidence, and imaging archives
(PACS) for radiological data. No single system integrates all of these data sources.
Clinicians and molecular pathologists must manually cross-reference findings across
systems, a process that is both time-consuming and error-prone.

**Challenge 3: Adult-Derived Therapies and Evidence Gaps.** The majority of oncology
knowledge bases, clinical decision-support tools, and treatment algorithms are designed
for adult cancers. Resources like OncoKB, CIViC, and the NCCN Clinical Practice
Guidelines contain relatively sparse pediatric-specific content. When a molecular tumor
board identifies a potentially actionable variant in a child's tumor, the evidence
supporting that variant's clinical significance often comes from adult studies. Drug
dosing, efficacy projections, and toxicity profiles derived from adult data may not
extrapolate reliably to pediatric patients, whose developing physiology, organ maturation
status, and pharmacokinetics differ substantially from adults.

**Challenge 4: Long-Term Toxicity in Survivors.** Unlike adult cancer patients, who have
a median age at diagnosis of 66 years, childhood cancer survivors face decades of life
after treatment. More than 500,000 childhood cancer survivors are alive in the US today,
and studies from the Childhood Cancer Survivor Study (CCSS) have shown that by age 50,
more than 99% of survivors have at least one chronic health condition, and more than
80% have a severe or life-threatening condition. Late effects include secondary
malignancies (standardized incidence ratio 6.4 compared to the general population),
cardiomyopathy from anthracycline exposure, neurocognitive impairment from cranial
radiation, endocrine dysfunction, infertility, and metabolic syndrome. Every treatment
decision for a child with cancer must balance the imperative to cure the cancer today
against the toxicity burden that survivor will carry for decades. This requires
pharmacogenomic profiling and toxicity prediction capabilities that are rarely
integrated into current tumor board workflows.

### 3.5 How the HCLS AI Factory Addresses Each Challenge

The HCLS AI Factory was designed to address each of the four challenges described above
through a combination of GPU-accelerated computation, multi-agent AI coordination,
unified vector-based data integration, and pediatric-specific knowledge curation.

**Addressing Challenge 1: Automated Multi-Agent Coordination.** The platform compresses
the molecular tumor board timeline from weeks to hours. Starting from raw FASTQ
sequencing data, the Genomic Foundation Engine performs alignment (BWA-MEM2 via
Parabricks) and variant calling (DeepVariant) in approximately 2 to 4 hours on a single
DGX Spark, compared to 24 to 48 hours on traditional CPU-based infrastructure. The
resulting VCF is automatically annotated, embedded into the Milvus vector database, and
made available to all downstream agents. The Precision Oncology Agent generates a
structured molecular tumor board report--with variant tiering per AMP/ASCO/CAP
guidelines, matched therapy recommendations, and clinical trial matches--within minutes.
The entire pipeline from FASTQ to MTB-ready report completes in under 5 hours.

**Addressing Challenge 2: Unified Data Platform.** The platform's Milvus vector database
serves as a unified knowledge substrate that spans all data modalities. Genomic variants,
clinical annotations (ClinVar, 4.1 million records), protein function predictions
(AlphaMissense, 71 million predictions), published literature, clinical trial data,
drug-target interaction data, and imaging features are all embedded into a shared vector
space using the BGE-small-en-v1.5 embedding model. This enables semantic search across
data types: a clinician can query "what are the treatment options for a child with
relapsed ALL and a TP53 R248W mutation" and receive results that span genomic databases,
clinical trials, published case reports, and drug databases--all retrieved through a
single search operation against the platform's 3.56 million searchable vectors.

**Addressing Challenge 3: Pediatric-Specific Knowledge.** The platform's intelligence
agents are designed to incorporate pediatric-specific evidence and contextualization.
The Precision Oncology Agent's therapy ranking algorithm incorporates pediatric drug
approvals, COG protocol references, and age-adjusted dosing recommendations. The
Pharmacogenomics Agent includes pediatric-specific dosing algorithms for the 25
pharmacogenes most relevant to childhood cancer treatment, including TPMT, NUDT15,
DPYD, UGT1A1, and CYP2D6. The CAR-T Intelligence Agent includes pediatric-specific
manufacturing considerations, such as the challenges of leukapheresis in small children,
T-cell fitness optimization for the characteristically naive T-cell repertoire of
pediatric patients, and pediatric-specific cytokine release syndrome management
protocols.

**Addressing Challenge 4: Toxicity Prevention Through Pharmacogenomic Integration.**
The Pharmacogenomics Agent operates as an automated safety layer that screens every
patient's germline variants against a curated database of gene-drug interactions. For
a child about to begin ALL maintenance therapy with mercaptopurine, the agent
automatically checks TPMT and NUDT15 genotype and provides CPIC guideline-concordant
dosing recommendations. For a patient being considered for high-dose methotrexate, it
evaluates MTHFR status and renal function pharmacogenomic markers. This pharmacogenomic
screening is automatically triggered as part of the integrated assessment workflow,
ensuring that toxicity risk is evaluated before--not after--treatment decisions are made.
The integration of pharmacogenomic profiling into the standard workflow represents a
paradigm shift from reactive toxicity management (identifying drug sensitivities after
an adverse event occurs) to proactive toxicity prevention.

The following sections detail the platform architecture that enables these capabilities
and provide step-by-step instructions for demonstrating them through six pediatric
oncology workflows.

---

## Platform Overview

### 4.1 Three-Engine Architecture

The HCLS AI Factory is organized into three major computational engines, each responsible
for a distinct phase of the precision medicine pipeline. This three-engine architecture
enables modular deployment (any engine can operate independently), parallel processing
(engines can run concurrently on separate GPU partitions), and flexible entry points
(users can begin at any engine depending on their input data type).

#### Engine 1: Genomic Foundation Engine

The Genomic Foundation Engine transforms raw sequencing data into annotated, searchable
genomic variants. It is the computational bedrock upon which all downstream intelligence
and discovery workflows operate.

**Pipeline stages:**

1. **Alignment and Sorting** -- Raw FASTQ reads are aligned to the GRCh38 human reference
   genome using BWA-MEM2, accelerated by NVIDIA Parabricks. Parabricks leverages GPU
   acceleration to perform alignment and coordinate sorting in a single pass, producing
   a sorted BAM file. On the DGX Spark, whole-genome alignment completes in approximately
   60 to 120 minutes, compared to 12 to 24 hours on a 32-core CPU server.

2. **Variant Calling** -- The aligned BAM file is processed by DeepVariant, a deep
   learning-based variant caller developed by Google that uses a convolutional neural
   network to identify single nucleotide variants (SNVs) and small insertions/deletions
   (indels). DeepVariant achieves the highest accuracy among variant callers on the
   Genome in a Bottle benchmark datasets. The output is a VCF file containing
   approximately 4 to 5 million variant calls for a typical whole-genome sample, of
   which approximately 11.7 million sites are assessed in the demonstration genome
   (HG002 from the Genome in a Bottle Consortium).

3. **Annotation** -- The raw VCF is annotated with functional and clinical information
   from multiple databases: ClinVar (4.1 million records of variant-disease
   associations), AlphaMissense (71 million missense variant pathogenicity predictions
   from DeepMind), and Ensembl VEP (Variant Effect Predictor) for consequence annotation,
   gene mapping, and regulatory region identification.

4. **Embedding and Indexing** -- Annotated variants are embedded into dense vector
   representations using the BGE-small-en-v1.5 model and indexed into the Milvus vector
   database. The embedding captures both the genomic context (gene, consequence, allele
   frequency) and the clinical context (pathogenicity classification, disease
   associations, drug interactions) of each variant. The demonstration dataset produces
   3.56 million annotated variant vectors.

**Key components:** Parabricks 4.6 (DeepVariant + BWA-MEM2), ClinVar database, AlphaMissense
predictions, BGE-small-en-v1.5 embedding model, Milvus vector database.

**Interfaces:** Genomics API (Flask, port 5000), pipeline scripts (00 through 05 in
core/engines/genomic-foundation/), Nextflow DSL2 orchestration (hcls-orchestrator/main.nf).

**Demo genome:** HG002 (NA24385), an Ashkenazi Jewish male from the Genome in a Bottle
Consortium. This sample was chosen because it has the most comprehensive truth set
available for benchmarking variant calling accuracy.

#### Engine 2: Precision Intelligence Engine

The Precision Intelligence Engine is the analytical core of the platform--a constellation
of seven specialized AI agents, each backed by domain-specific Milvus vector collections
and connected through a shared genomic evidence layer. These agents collectively provide
clinical interpretation, risk stratification, therapy matching, pharmacogenomic screening,
imaging correlation, immunotherapy design, and autoimmune assessment.

Each agent follows a common architectural pattern:

- **FastAPI backend** serving RESTful endpoints including /analyze, /search,
  /integrated-assessment, and /health
- **Streamlit frontend** providing an interactive user interface for clinicians and
  researchers
- **Domain-specific Milvus collections** (ranging from 6 to 15 collections per agent)
  containing curated, embedded knowledge
- **Shared genomic_evidence collection** enabling cross-agent genomic context retrieval
- **Claude LLM integration** (via Anthropic API) for natural language reasoning,
  report generation, and conversational interaction
- **Health monitoring** with automatic restart via the platform's health-monitor.sh
  watchdog

The seven operational agents are described in detail in Section 4.2 below.

**Key components:** Seven FastAPI/Streamlit agent pairs, 80+ Milvus vector collections,
Claude AI (Anthropic) for LLM reasoning, BGE-small-en-v1.5 for embedding, shared
hcls_common library (23 modules).

**Interfaces:** Individual agent APIs (ports 8107-8540), Streamlit UIs (ports 8502-8540),
RAG Chat API (Flask, port 5001), RAG Chat UI (Streamlit, port 8501), Discovery Portal
(Streamlit, port 8510).

#### Engine 3: Therapeutic Discovery Engine

The Therapeutic Discovery Engine translates molecular targets identified by the
Precision Intelligence Engine into novel drug candidates. It implements a six-phase
computational chemistry pipeline:

1. **Target Import** -- Receives target protein structures and binding site information
   from the intelligence agents. Targets may come from the Precision Oncology Agent
   (e.g., a mutant kinase driving a child's tumor), the CAR-T Agent (e.g., a surface
   antigen for CAR targeting), or the Biomarker Agent (e.g., a protein biomarker
   associated with a specific disease trajectory).

2. **Molecule Generation** -- BioNeMo MolMIM (Molecular Inverse Mapping) generates
   novel molecular structures optimized for the target binding site. MolMIM uses a
   variational autoencoder trained on large chemical libraries to explore chemical
   space efficiently, producing diverse candidate molecules with predicted binding
   affinity.

3. **Molecular Docking** -- DiffDock, a diffusion-based molecular docking model,
   predicts the binding pose and affinity of each generated molecule to the target
   protein. Unlike traditional docking tools (AutoDock, Glide), DiffDock treats
   docking as a generative diffusion process, producing multiple plausible binding
   configurations with confidence scores.

4. **Property Scoring** -- RDKit calculates ADMET (absorption, distribution, metabolism,
   excretion, toxicity) properties for each candidate molecule, including Lipinski's
   Rule of Five compliance, topological polar surface area (TPSA), predicted logP,
   synthetic accessibility score, and drug-likeness metrics.

5. **Pharmacogenomic Filtering** -- Candidate molecules are cross-referenced against
   the Pharmacogenomics Agent's database to flag potential gene-drug interactions
   specific to the patient's germline genotype. This step ensures that generated
   candidates are compatible with the individual patient's metabolic profile.

6. **Report Generation** -- A comprehensive PDF report is generated containing ranked
   candidate molecules with docking scores, ADMET profiles, 3D binding pose
   visualizations, and pharmacogenomic compatibility assessments.

**Key components:** BioNeMo MolMIM (molecule generation), DiffDock (molecular docking),
RDKit (property calculation and scoring), cryo-EM structure integration, Pharmacogenomics
Agent cross-referencing.

**Demo target:** VCP (Valosin-Containing Protein) for frontotemporal dementia--used
as the reference demonstration target, with 4 cryo-EM structures, CB-5083 as the
reference inhibitor, and 100 generated analogues. Pediatric oncology demos substitute
disease-relevant targets (e.g., histone demethylase for DMG, CD19 structure for
CAR-T design).

**Interfaces:** Drug Discovery UI (Streamlit, port 8505), Discovery Portal (Streamlit,
port 8510), Grafana monitoring (port 3000), Prometheus metrics (port 9090).

### 4.2 The Agent Constellation

The platform's seven operational intelligence agents represent specialized domains of
clinical knowledge, each designed to function both independently and as part of a
coordinated multi-agent workflow.

| # | Agent | Ports (API/UI) | Collections | Description |
|---|-------|---------------|-------------|-------------|
| 1 | **Precision Biomarker Agent** | 8529 / 8502 | 10 + shared | Genotype-aware biomarker interpretation with biological age estimation (PhenoAge/GrimAge), disease trajectory detection, and pharmacogenomic profiling. Outputs clinical reports in PDF and FHIR R4 format. |
| 2 | **Precision Oncology Agent** | 8527 / 8503 | 10 + shared | Molecular tumor board decision support with CIViC/OncoKB variant annotation, AMP/ASCO/CAP-tiered therapy ranking, and clinical trial matching. Outputs MTB-ready structured reports. |
| 3 | **CAR-T Intelligence Agent** | 8522 / 8504 | 11 | Cross-functional CAR-T cell therapy intelligence spanning construct design, manufacturing optimization, clinical outcomes, toxicity management, and comparative analysis (e.g., 4-1BB vs. CD28 costimulatory domains). 6,266+ curated vectors. |
| 4 | **Imaging Intelligence Agent** | 8524 / 8505 | 10 + shared | Medical imaging AI integrating NVIDIA NIM models (VISTA-3D, MAISI, VILA-M3) with DICOM auto-ingestion via Orthanc and federated learning via NVIDIA FLARE. Six clinical workflows: CT Head, CXR, CT Lung, MRI, Hemorrhage, Brain. Cross-modal genomics enrichment. |
| 5 | **Precision Autoimmune Agent** | 8532 / 8506 | -- | Autoimmune disease and immune-mediated condition intelligence supporting differential diagnosis, treatment selection, and monitoring for autoimmune complications of cancer immunotherapy. |
| 6 | **Pharmacogenomics Agent** | 8107 / 8507 | 15 | RAG-powered pharmacogenomic clinical decision support covering 25 pharmacogenes, 100+ drugs, 9 dosing algorithms, and 15 HLA associations. 1,001+ validated tests. Critical for pediatric chemotherapy safety. |
| 7 | **Cardiology Intelligence Agent** | 8126 / 8527 | 12 + shared | RAG-powered cardiovascular clinical decision support with 6 validated risk calculators (ASCVD, HEART, CHA2DS2-VASc, HAS-BLED, MAGGIC, EuroSCORE II), GDMT optimizer for heart failure, and 8 clinical workflows including cardio-oncology. 1,927 validated tests. |

**Agents in design (not yet operational):**

| # | Agent | Status | Description |
|---|-------|--------|-------------|
| 8 | **Neurology Agent** | Design phase | Neurological disease intelligence, relevant to pediatric brain tumor cognitive monitoring |
| 9 | **Rare Disease Diagnostic Agent** | Design phase | Rare disease diagnostic support, highly relevant to many pediatric cancers classified as rare diseases |

### 4.3 Cross-Agent Integration

The agents do not operate in isolation. The platform provides multiple mechanisms for
cross-agent coordination, enabling the kind of multidisciplinary analysis that
characterizes effective tumor board decision-making.

**Integrated Assessment Endpoint.** Each agent exposes an /integrated-assessment API
endpoint that accepts a patient case and returns a domain-specific assessment. The
platform's meta-agent orchestrator (lib/hcls_common/meta_agent.py, 80KB) coordinates
calls across multiple agents, merges results, resolves conflicts, and generates a
unified report. For a pediatric oncology case, a single /integrated-assessment call
can trigger parallel analysis by the Oncology Agent (variant interpretation and therapy
matching), Biomarker Agent (prognostic biomarker profiling), Pharmacogenomics Agent
(drug interaction screening), Imaging Agent (radiological correlation), and Cardiology
Agent (cardio-oncology toxicity risk assessment).

**Bidirectional Triggers.** The lib/hcls_common/bidirectional_triggers.py module (43KB)
implements an event-driven trigger system that allows agents to automatically invoke
related agents based on their findings. For example, if the Oncology Agent identifies
an anthracycline-based chemotherapy recommendation, it automatically triggers the
Cardiology Agent to assess cardiotoxicity risk based on the patient's germline variants
in genes associated with anthracycline-induced cardiomyopathy (e.g., RARG, CBR3, HAS3).
If the Biomarker Agent identifies a TPMT poor-metabolizer genotype, it triggers the
Pharmacogenomics Agent for detailed thiopurine dosing guidance.

**Event Bus.** The lib/hcls_common/event_bus.py module (31KB) provides a publish-subscribe
messaging backbone that enables asynchronous communication between agents. Agents
publish events (e.g., "variant_classified", "therapy_recommended", "toxicity_flagged")
and subscribe to events relevant to their domain. This enables loosely coupled
integration where new agents can be added to the constellation without modifying
existing agents.

**Query Router.** The lib/hcls_common/query_router.py module (44KB) implements
intelligent routing of user queries to the most appropriate agent or combination of
agents. Natural language queries are analyzed for domain signals and routed accordingly:
a query about "MYCN amplification prognosis in neuroblastoma" routes to the Oncology
Agent, while "mercaptopurine dosing with TPMT*3A" routes to the Pharmacogenomics Agent,
and "cardiac monitoring after doxorubicin in a 5-year-old" routes to both the Cardiology
and Pharmacogenomics Agents.

### 4.4 Shared Genomic Evidence Layer

All agents share access to a common genomic_evidence collection in Milvus. This
collection contains the embedded variant data produced by the Genomic Foundation Engine
and serves as the single source of truth for patient-specific genomic context. When any
agent needs to reference a patient's genomic data--whether the Oncology Agent is looking
for somatic driver mutations, the Pharmacogenomics Agent is checking germline metabolizer
status, or the Cardiology Agent is evaluating inherited cardiomyopathy risk--they query
the same shared collection.

This shared-evidence architecture provides several critical benefits:

- **Consistency:** All agents reference the same variant calls and annotations,
  eliminating the risk of conflicting genomic data across reports.

- **Efficiency:** Variants are embedded and indexed once, then queried many times by
  different agents. This avoids redundant processing and ensures that the computational
  cost of embedding scales with the number of variants, not the number of agents.

- **Discoverability:** Because the evidence is stored as dense vectors, agents can
  perform semantic similarity searches that surface relevant variants even when the
  exact gene name or variant identifier is not specified. A query for "kinase domain
  mutation in a pediatric brain tumor gene" can retrieve relevant BRAF, ALK, or ROS1
  variants without requiring the user to enumerate specific genes.

- **Auditability:** The shared collection provides a complete audit trail of which
  variants were available to which agents at the time of analysis, supporting the
  reproducibility requirements essential for clinical-grade genomic interpretation.

The genomic_evidence collection is populated during the annotation and embedding phase
of the Genomic Foundation Engine pipeline and is available to all agents immediately
upon completion. For VCF entry-point workflows, pre-annotated and pre-embedded VCF
data can be loaded directly into the collection, bypassing the alignment and variant
calling stages.

### 4.5 Dual Entry Points

Every demonstration workflow in this document supports two entry points, designed to
accommodate different audiences, time constraints, and demonstration objectives.

#### FASTQ Entry Point (Full Pipeline)

**Input:** Paired-end FASTQ files (raw sequencing reads)

**Pipeline:** Genomic Foundation Engine (alignment, variant calling, annotation,
embedding) followed by Precision Intelligence Engine (agent analysis) followed by
Therapeutic Discovery Engine (molecule generation, docking, scoring)

**Approximate duration:** 4 to 5 hours for a whole-genome sample on DGX Spark

**Best for:**
- Demonstrating end-to-end capability from raw data to drug candidates
- Audiences with genomics or bioinformatics backgrounds
- Full platform validation and benchmarking
- Showcasing GPU acceleration advantages (120-240 min vs. 24-48 hours on CPU)

**Demo data:** HG002 (NA24385) whole-genome FASTQ files from the Genome in a Bottle
Consortium, approximately 30x coverage, GRCh38 reference. These files are included in
the hcls-ai-factory-core-data package (1.4 GB).

#### VCF Entry Point (Intelligence + Discovery)

**Input:** Annotated VCF file (pre-called variants)

**Pipeline:** Annotation and embedding (if not already annotated) followed by Precision
Intelligence Network (agent analysis) followed by Therapeutic Discovery Engine (molecule
generation, docking, scoring)

**Approximate duration:** 20 to 30 minutes

**Best for:**
- Clinical audiences focused on interpretation and therapeutic outcomes
- Time-constrained demonstration slots (conference booths, executive briefings)
- Iterative workflow development and agent testing
- Audiences less interested in the bioinformatics pipeline details

**Demo data:** Pre-called VCF files from HG002, with annotations pre-computed and stored
in the F3_inputs directory (3.56 million annotated variants in JSONL.gz format, ready
for embedding).

The choice of entry point is made at the orchestrator level. The Nextflow orchestrator
(hcls-orchestrator/main.nf) accepts a mode parameter that supports full_pipeline (FASTQ
entry), rag_only (VCF entry, intelligence only), drug_discovery_only (VCF entry,
discovery only), and genomics_only (FASTQ to VCF, no intelligence or discovery). The
Python fallback orchestrator (hcls-orchestrator/run_pipeline.py) supports the same modes.

### 4.6 Service Port Map

The following table lists all platform services, their assigned ports, and their roles
in the demonstration workflows. All services run on the local DGX Spark and are
accessible via localhost.

| Service | Port | Protocol | Framework | Role in Demos |
|---------|------|----------|-----------|---------------|
| **Landing Page** | 8080 | HTTP | Flask | Central hub with real-time health monitoring for all services. Starting point for all demonstrations. |
| **Genomics API** | 5000 | HTTP | Flask | Genomic Foundation Engine control plane. Triggers alignment, variant calling, and annotation pipelines. |
| **RAG Chat API** | 5001 | HTTP | Flask | Precision Intelligence Engine query endpoint. Handles semantic search, RAG retrieval, and LLM-powered responses. |
| **RAG Chat UI** | 8501 | HTTP | Streamlit | Interactive chat interface for querying the variant knowledge base. Used in all workflows for ad-hoc exploration. |
| **Biomarker Agent API** | 8529 | HTTP | FastAPI | Precision Biomarker Agent backend. Provides biomarker interpretation, biological age, pharmacogenomic profiling. |
| **Biomarker Agent UI** | 8502 | HTTP | Streamlit | Biomarker Agent interactive interface. Used in Workflows 1, 2, and 6. |
| **Oncology Agent API** | 8527 | HTTP | FastAPI | Precision Oncology Agent backend. Provides variant annotation, therapy ranking, trial matching, MTB reports. |
| **Oncology Agent UI** | 8503 | HTTP | Streamlit | Oncology Agent interactive interface. Primary interface for Workflows 1, 2, 3, and 6. |
| **CAR-T Agent API** | 8522 | HTTP | FastAPI | CAR-T Intelligence Agent backend. Provides construct analysis, manufacturing guidance, outcome prediction. |
| **CAR-T Agent UI** | 8504 | HTTP | Streamlit | CAR-T Agent interactive interface. Primary interface for Workflow 4. |
| **Imaging Agent API** | 8524 | HTTP | FastAPI | Imaging Intelligence Agent backend. Provides medical imaging AI analysis with NVIDIA NIM integration. |
| **Imaging Agent UI** | 8505 | HTTP | Streamlit | Imaging Agent interactive interface. Used in Workflows 2 and 6. |
| **Autoimmune Agent API** | 8532 | HTTP | FastAPI | Precision Autoimmune Agent backend. Provides autoimmune assessment for immunotherapy complications. |
| **Autoimmune Agent UI** | 8506 | HTTP | Streamlit | Autoimmune Agent interactive interface. Used in Workflow 6. |
| **Pharmacogenomics Agent API** | 8107 | HTTP | FastAPI | Pharmacogenomics Agent backend. Provides PGx dosing guidance for 25 pharmacogenes and 100+ drugs. |
| **Pharmacogenomics Agent UI** | 8507 | HTTP | Streamlit | Pharmacogenomics Agent interactive interface. Primary interface for Workflow 5. |
| **Cardiology Agent API** | 8126 | HTTP | FastAPI | Cardiology Intelligence Agent backend. Provides cardiovascular risk assessment and cardio-oncology support. |
| **Cardiology Agent UI** | 8527 | HTTP | Streamlit | Cardiology Agent interactive interface. Used in Workflows 1, 5, and 6 for cardiotoxicity assessment. |
| **Discovery Portal** | 8510 | HTTP | Streamlit | Unified portal for Therapeutic Discovery Engine. Drug candidate visualization and report generation. |
| **Milvus Vector DB** | 19530 | gRPC | Milvus | Vector database storing 3.56M variant embeddings and 80+ domain-specific collections. |
| **Milvus Management** | 9091 | HTTP | Milvus | Milvus health and management API. |
| **Grafana** | 3000 | HTTP | Grafana | Real-time monitoring dashboards for GPU utilization, pipeline progress, and service health. |
| **Prometheus** | 9090 | HTTP | Prometheus | Metrics collection and aggregation for all platform services. |
| **Node Exporter** | 9100 | HTTP | Prometheus | System-level metrics (CPU, memory, disk, network) for the DGX Spark. |
| **DCGM Exporter** | 9400 | HTTP | NVIDIA | GPU-level metrics (utilization, memory, temperature, power) via NVIDIA Data Center GPU Manager. |

**Network topology:** All services communicate over the local Docker bridge network.
The Milvus deployment includes three supporting containers: etcd (distributed
key-value store for metadata), MinIO (S3-compatible object storage for vector data),
and the Milvus standalone server. Agent containers connect to Milvus via the gRPC
interface on port 19530 and to the Claude API via outbound HTTPS.

**Resource requirements:** The full platform requires approximately 128 GB of system
RAM (Milvus alone requires 32-64 GB for the full vector index), GPU memory proportional
to the active workload (alignment requires the most GPU memory), and approximately
1.1 TB of disk space including reference genomes, databases, model weights, and the
vector index.

---

## Section 5: DEMO 1 -- End-to-End Precision Medicine Pipeline (The Foundation Demo)

### Clinical Narrative

**Patient:** Evelyn R., 8-year-old female
**Presenting complaint:** Progressive fatigue, pallor, easy bruising for 3 weeks, fever (38.9C) unresponsive to antipyretics, petechiae on bilateral lower extremities.

**Initial workup:**
- CBC: WBC 45,000/uL (90% lymphoblasts), Hgb 6.8 g/dL, Platelets 18,000/uL
- Peripheral smear: >90% lymphoblasts, small to medium size, scant cytoplasm, fine chromatin
- LDH: 2,340 U/L (elevated)
- Uric acid: 9.2 mg/dL (elevated, tumor lysis risk)

**Bone marrow biopsy:**
- >95% cellularity, replaced by lymphoblasts
- Flow cytometry: CD10+, CD19+, CD22+, CD34+, TdT+, CD20 dim — consistent with B-cell precursor ALL
- Cytogenetics: 46,XX,t(12;21)(p13;q22)[18]/46,XX[2]
- FISH: ETV6-RUNX1 fusion positive (favorable)
- Additional findings: IKZF1 deletion detected by MLPA (unfavorable modifier)

**Diagnosis:** B-cell precursor acute lymphoblastic leukemia (BCP-ALL), NCI standard-risk by age/WBC but with IKZF1 deletion (IKZF1plus consideration)

**Clinical question:** What is the molecular risk classification considering the conflicting prognostic markers (favorable ETV6-RUNX1 vs. unfavorable IKZF1 deletion)? What targeted therapeutic candidates exist?

**Agents involved:** Genomic Foundation Engine (:5000), RAG Engine (:5001), Precision Oncology Agent (:8527), Therapeutic Discovery Engine (:8510)

---

### Architecture Overview

```
+-------------------------------+       +-------------------------------+       +---------------------------+
|   GENOMIC FOUNDATION ENGINE   |       |   PRECISION INTELLIGENCE      |       |   THERAPEUTIC DISCOVERY   |
|                               |       |   NETWORK                     |       |   ENGINE                  |
|  FASTQ --> BWA-MEM2 --> BAM   |       |                               |       |                           |
|  BAM --> DeepVariant --> VCF   | ----> |  RAG Engine (:5001)           | ----> |  MolMIM (:8510)           |
|  VCF --> ClinVar/AlphaMiss    |       |  Oncology Agent (:8527)       |       |  DiffDock (:8510)         |
|  Annotated VCF --> Milvus     |       |                               |       |  RDKit ADMET (:8510)      |
+-------------------------------+       +-------------------------------+       +---------------------------+
```

This is the **foundation demo** that exercises all three platform engines in sequence. Every subsequent demo builds on the patterns established here.

---

### Entry Point A: Starting from FASTQ

#### Step A1 -- Upload FASTQ and Configure Alignment

Evelyn's tumor sample (bone marrow aspirate) and germline sample (buccal swab) are processed through the Genomic Foundation Engine. The paired-end reads are aligned against GRCh38 using BWA-MEM2 accelerated by NVIDIA Parabricks.

```bash
# Create input directories
mkdir -p /data/genomics/inputs/sophia_tumor
mkdir -p /data/genomics/inputs/sophia_germline

# Copy FASTQ files (tumor)
cp sophia_tumor_R1.fastq.gz /data/genomics/inputs/sophia_tumor/
cp sophia_tumor_R2.fastq.gz /data/genomics/inputs/sophia_tumor/

# Submit alignment + variant calling pipeline
curl -X POST http://localhost:5000/v1/genomics/run \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "sophia_tumor_001",
    "pipeline_type": "somatic_tumor_normal",
    "input_fastq_r1": "/data/genomics/inputs/sophia_tumor/sophia_tumor_R1.fastq.gz",
    "input_fastq_r2": "/data/genomics/inputs/sophia_tumor/sophia_tumor_R2.fastq.gz",
    "normal_fastq_r1": "/data/genomics/inputs/sophia_germline/sophia_germline_R1.fastq.gz",
    "normal_fastq_r2": "/data/genomics/inputs/sophia_germline/sophia_germline_R2.fastq.gz",
    "reference_genome": "GRCh38",
    "target_regions": "all_coding_plus_regulatory",
    "caller": "deepvariant",
    "output_dir": "/data/genomics/outputs/sophia_tumor_001"
  }'
```

**Expected response:**

```json
{
  "run_id": "GEN-2026-0323-001",
  "sample_id": "sophia_tumor_001",
  "status": "running",
  "pipeline_type": "somatic_tumor_normal",
  "stages": [
    {"name": "bwa_mem2_alignment", "status": "running", "gpu_accelerated": true},
    {"name": "deepvariant_calling", "status": "pending", "gpu_accelerated": true},
    {"name": "clinvar_annotation", "status": "pending"},
    {"name": "alphamissense_scoring", "status": "pending"},
    {"name": "milvus_embedding", "status": "pending"}
  ],
  "estimated_completion_minutes": 135,
  "output_vcf": "/data/genomics/outputs/sophia_tumor_001/sophia_tumor_001.deepvariant.vcf.gz"
}
```

**Presenter note:** Alignment with BWA-MEM2 on DGX Spark completes in ~45 minutes for 30x WGS versus 8-12 hours on CPU-only infrastructure. DeepVariant calling adds ~60 minutes GPU-accelerated versus 24+ hours on CPU.

---

#### Step A2 -- Monitor Pipeline Progress

```bash
curl -X GET http://localhost:5000/v1/genomics/status/GEN-2026-0323-001 \
  -H "Content-Type: application/json"
```

**Expected response (mid-run):**

```json
{
  "run_id": "GEN-2026-0323-001",
  "sample_id": "sophia_tumor_001",
  "status": "running",
  "elapsed_minutes": 52,
  "stages": [
    {"name": "bwa_mem2_alignment", "status": "completed", "duration_minutes": 47, "reads_aligned": 892450000, "alignment_rate": 0.9987},
    {"name": "deepvariant_calling", "status": "running", "progress_pct": 18, "gpu_utilization": 0.94},
    {"name": "clinvar_annotation", "status": "pending"},
    {"name": "alphamissense_scoring", "status": "pending"},
    {"name": "milvus_embedding", "status": "pending"}
  ]
}
```

---

#### Step A3 -- Retrieve Completed VCF with Annotations

```bash
curl -X GET http://localhost:5000/v1/genomics/results/GEN-2026-0323-001 \
  -H "Content-Type: application/json"
```

**Expected response:**

```json
{
  "run_id": "GEN-2026-0323-001",
  "sample_id": "sophia_tumor_001",
  "status": "completed",
  "total_duration_minutes": 132,
  "output_files": {
    "vcf": "/data/genomics/outputs/sophia_tumor_001/sophia_tumor_001.deepvariant.vcf.gz",
    "annotated_vcf": "/data/genomics/outputs/sophia_tumor_001/sophia_tumor_001.annotated.vcf.gz",
    "bam": "/data/genomics/outputs/sophia_tumor_001/sophia_tumor_001.sorted.bam",
    "qc_report": "/data/genomics/outputs/sophia_tumor_001/qc_metrics.json"
  },
  "variant_summary": {
    "total_variants": 4287650,
    "coding_variants": 24891,
    "clinvar_pathogenic": 3,
    "clinvar_likely_pathogenic": 7,
    "alphamissense_high": 12,
    "somatic_variants": 847,
    "key_findings": [
      "ETV6-RUNX1 fusion: t(12;21)(p13;q22) confirmed",
      "IKZF1 exon 4-7 deletion detected",
      "PAX5 p.P80R missense (VUS, AlphaMissense 0.73)",
      "CDKN2A homozygous deletion"
    ]
  },
  "milvus_ingestion": {
    "vectors_embedded": 24891,
    "collection": "sophia_tumor_001_variants",
    "embedding_model": "bge-small-en-v1.5"
  }
}
```

---

#### Step A4 -- Verify Milvus Vector Embedding

```bash
curl -X POST http://localhost:5001/v1/rag/search \
  -H "Content-Type: application/json" \
  -d '{
    "query": "ETV6-RUNX1 fusion pediatric ALL prognosis",
    "collection": "sophia_tumor_001_variants",
    "top_k": 5,
    "include_metadata": true
  }'
```

**Expected response:**

```json
{
  "query": "ETV6-RUNX1 fusion pediatric ALL prognosis",
  "results": [
    {
      "rank": 1,
      "score": 0.924,
      "variant": "t(12;21)(p13;q22) ETV6::RUNX1",
      "source": "sophia_tumor_001_variants",
      "metadata": {
        "gene": "ETV6-RUNX1",
        "type": "fusion",
        "clinvar_significance": "Pathogenic",
        "disease": "B-cell acute lymphoblastic leukemia",
        "evidence_level": "A",
        "interpretation": "Favorable prognostic marker in pediatric BCP-ALL. Associated with 5-year EFS >90% in standard-risk protocols."
      }
    },
    {
      "rank": 2,
      "score": 0.891,
      "variant": "IKZF1 exon 4-7 deletion",
      "source": "sophia_tumor_001_variants",
      "metadata": {
        "gene": "IKZF1",
        "type": "deletion",
        "clinvar_significance": "Pathogenic",
        "disease": "B-cell acute lymphoblastic leukemia",
        "evidence_level": "A",
        "interpretation": "Unfavorable prognostic modifier. IKZF1 deletion in ETV6-RUNX1+ ALL defines IKZF1plus phenotype with inferior EFS (5yr EFS ~65% vs ~92%)."
      }
    }
  ],
  "total_results": 5,
  "search_time_ms": 23
}
```

**This completes the genomic processing pipeline.** The VCF is annotated, embedded in Milvus, and ready for intelligence agent queries.

---

### Entry Point B: Starting from VCF

For audiences focused on clinical intelligence rather than genomic processing, begin directly with a pre-annotated VCF file.

```bash
curl -X POST http://localhost:5000/v1/genomics/ingest-vcf \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "sophia_tumor_001",
    "vcf_path": "/data/demo_vcfs/sophia_etv6runx1_ikzf1del.annotated.vcf.gz",
    "patient_metadata": {
      "age_years": 8,
      "sex": "female",
      "diagnosis": "B-cell precursor ALL",
      "wbc_at_diagnosis": 45000,
      "cns_status": "CNS1",
      "cytogenetics": "t(12;21)(p13;q22)",
      "additional_findings": ["IKZF1 deletion", "CDKN2A deletion"]
    },
    "embed_to_milvus": true,
    "collection_name": "sophia_tumor_001_variants"
  }'
```

**Expected response:**

```json
{
  "sample_id": "sophia_tumor_001",
  "status": "ingested",
  "variants_processed": 24891,
  "vectors_embedded": 24891,
  "clinvar_matches": 187,
  "alphamissense_scored": 1243,
  "collection": "sophia_tumor_001_variants",
  "ingestion_time_seconds": 34
}
```

---

### Step 6 -- Precision Oncology Agent: Molecular Risk Classification

With variants embedded in Milvus, query the Precision Oncology Agent for integrated molecular assessment.

```bash
curl -X POST http://localhost:8527/v1/oncology/molecular-classification \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "diagnosis": "B-cell precursor ALL",
    "age_years": 8,
    "sex": "female",
    "wbc_at_diagnosis": 45000,
    "variant_collection": "sophia_tumor_001_variants",
    "classification_system": "NCI_risk_with_molecular",
    "include_trial_matching": true
  }'
```

**Expected response:**

```json
{
  "patient_id": "sophia_tumor_001",
  "classification": {
    "nci_risk": "standard_risk",
    "nci_criteria": {
      "age": "1-9.99 years (favorable)",
      "wbc": "45,000/uL (borderline — threshold is 50,000)"
    },
    "molecular_risk_modifiers": [
      {
        "marker": "ETV6-RUNX1 fusion",
        "impact": "favorable",
        "evidence_level": "1A",
        "detail": "t(12;21)(p13;q22) present in ~25% of pediatric BCP-ALL. Associated with excellent prognosis: 5-year EFS 90-95% on standard therapy (COG AALL0331).",
        "sources": ["Pui et al., NEJM 2015", "CIViC EID:234"]
      },
      {
        "marker": "IKZF1 deletion (exons 4-7)",
        "impact": "unfavorable",
        "evidence_level": "1A",
        "detail": "IKZF1 deletion is an independent adverse prognostic factor. In ETV6-RUNX1+ patients, co-occurring IKZF1 deletion defines the IKZF1plus phenotype associated with 5-year EFS ~65% (vs ~92% for IKZF1-intact ETV6-RUNX1+ ALL).",
        "sources": ["Mullighan et al., PNAS 2008", "Stanulla et al., JCO 2018"]
      },
      {
        "marker": "CDKN2A homozygous deletion",
        "impact": "unfavorable_modifier",
        "evidence_level": "2B",
        "detail": "Loss of CDKN2A (9p21.3) occurs in ~30% of ALL. Associated with proliferative advantage but prognostic significance debated when co-occurring with ETV6-RUNX1."
      },
      {
        "marker": "PAX5 p.P80R",
        "impact": "uncertain",
        "evidence_level": "3",
        "detail": "PAX5 mutations occur in ~30% of BCP-ALL. This variant (AlphaMissense score 0.73) likely disrupts the paired domain. Clinical significance as independent prognostic marker is uncertain."
      }
    ],
    "integrated_risk": "standard_risk_with_unfavorable_genetics",
    "recommended_protocol": "COG AALL0932 — standard risk, high-risk arm consideration due to IKZF1plus",
    "treatment_implications": [
      "Consider intensification to high-risk arm per IKZF1plus criteria",
      "Day 8 prednisone response and Day 29 MRD critical for final risk assignment",
      "If MRD >= 0.01% at Day 29, escalate to high-risk therapy per AALL1231"
    ]
  },
  "matched_trials": [
    {
      "trial_id": "NCT03914625",
      "title": "COG AALL1732: Improving Outcomes for Newly Diagnosed Standard Risk ALL",
      "phase": "III",
      "relevance_score": 0.96,
      "match_criteria": ["Age 1-9", "Standard risk BCP-ALL", "ETV6-RUNX1+"]
    },
    {
      "trial_id": "NCT03150693",
      "title": "Pediatric MATCH: Targeted Therapy for Relapsed/Refractory Solid Tumors and Lymphomas",
      "phase": "II",
      "relevance_score": 0.72,
      "match_criteria": ["CDKN2A deletion — palbociclib arm"]
    }
  ],
  "processing_time_ms": 1847
}
```

**Presenter note:** Highlight the conflicting prognostic signals: ETV6-RUNX1 is the most favorable genetic marker in childhood ALL, yet IKZF1 deletion significantly worsens prognosis. This is exactly the scenario where integrated AI assessment adds clinical value — a manual review might miss the IKZF1plus interaction.

---

### Step 7 -- RAG Evidence Synthesis

Query the RAG engine to synthesize evidence across embedded variant data, ClinVar records, and literature.

```bash
curl -X POST http://localhost:5001/v1/rag/synthesize \
  -H "Content-Type: application/json" \
  -d '{
    "query": "What is the clinical significance of co-occurring ETV6-RUNX1 fusion and IKZF1 deletion in pediatric B-ALL? What are the treatment implications and recommended monitoring?",
    "patient_context": {
      "patient_id": "sophia_tumor_001",
      "variant_collection": "sophia_tumor_001_variants"
    },
    "sources": ["clinvar", "civic", "pubmed_embeddings", "patient_variants"],
    "max_context_chunks": 20,
    "response_format": "clinical_summary"
  }'
```

**Expected response:**

```json
{
  "query": "What is the clinical significance of co-occurring ETV6-RUNX1 fusion and IKZF1 deletion in pediatric B-ALL?",
  "synthesis": {
    "summary": "The co-occurrence of ETV6-RUNX1 fusion and IKZF1 deletion in pediatric BCP-ALL represents a clinically significant prognostic scenario. While ETV6-RUNX1 alone is the most favorable cytogenetic subgroup in childhood ALL (5-year EFS 90-95%), the addition of IKZF1 deletion defines the 'IKZF1plus' phenotype identified by Stanulla et al. (2018), which is associated with substantially inferior outcomes (5-year EFS ~65%). This patient meets IKZF1plus criteria: IKZF1 deletion plus a co-occurring deletion of CDKN2A without ERG deletion.",
    "key_evidence": [
      {
        "finding": "IKZF1plus definition met",
        "detail": "IKZF1 deletion + CDKN2A deletion + absence of ERG deletion = IKZF1plus phenotype per Stanulla et al., JCO 2018",
        "confidence": 0.97,
        "sources": ["ClinVar:VCV000632891", "PMID:29792307", "sophia_tumor_001_variants:IKZF1_del"]
      },
      {
        "finding": "Risk reclassification recommended",
        "detail": "Multiple international consortia (COG, AIEOP-BFM, UKALL) now incorporate IKZF1 status into risk stratification. COG protocols recommend consideration for high-risk arm assignment when IKZF1plus criteria are met, regardless of initial NCI risk group.",
        "confidence": 0.94,
        "sources": ["PMID:31697823", "PMID:30670444"]
      },
      {
        "finding": "MRD-guided therapy critical",
        "detail": "Day 29 MRD assessment is the single most important early response marker. In IKZF1plus patients, MRD >= 0.01% at Day 29 is associated with 5-year EFS <50%, supporting immediate escalation to augmented therapy.",
        "confidence": 0.92,
        "sources": ["PMID:32726363", "CIViC:EID1847"]
      }
    ],
    "therapeutic_context": "Standard induction per COG AALL0932 should proceed with prednisone. Day 8 peripheral blood response and Day 29 bone marrow MRD will guide final risk assignment. If IKZF1plus criteria are confirmed and MRD is detectable, augmented consolidation (per AALL1231 high-risk arm) is recommended.",
    "monitoring_recommendations": [
      "Day 8: Peripheral blood blast count (prednisone response)",
      "Day 29: Bone marrow MRD by flow cytometry (threshold: 0.01%)",
      "End of consolidation: MRD reassessment",
      "Monitor for late relapse — ETV6-RUNX1+ ALL has characteristic late relapse pattern (3-5 years)"
    ]
  },
  "context_chunks_used": 17,
  "processing_time_ms": 3241
}
```

---

### Step 8 -- Therapeutic Discovery Engine: Novel Candidate Generation

With the molecular profile established, engage the Therapeutic Discovery Engine to identify and generate potential therapeutic candidates targeting the specific molecular vulnerabilities.

```bash
curl -X POST http://localhost:8510/v1/discovery/generate-candidates \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "targets": [
      {
        "gene": "IKZF1",
        "alteration": "deletion_exons_4_7",
        "rationale": "Restore IKAROS tumor suppressor function or target downstream pathway"
      },
      {
        "gene": "CDKN2A",
        "alteration": "homozygous_deletion",
        "rationale": "CDK4/6 inhibition to compensate for p16INK4a loss"
      }
    ],
    "constraints": {
      "pediatric_safety_filter": true,
      "max_molecular_weight": 500,
      "min_oral_bioavailability": 0.3,
      "exclude_known_cardiotoxins": true,
      "exclude_known_neurotoxins": true,
      "age_group": "pediatric_6_12"
    },
    "discovery_pipeline": {
      "generator": "molmim",
      "num_candidates": 50,
      "docking_engine": "diffdock",
      "scoring": ["binding_affinity", "admet_rdkit", "pediatric_safety", "druglikeness"]
    }
  }'
```

**Expected response:**

```json
{
  "patient_id": "sophia_tumor_001",
  "discovery_run_id": "DISC-2026-0323-001",
  "status": "completed",
  "candidates_generated": 50,
  "candidates_passing_filters": 12,
  "top_candidates": [
    {
      "rank": 1,
      "candidate_id": "HCLS-CDK46-001",
      "name": "Palbociclib analog (MolMIM-generated)",
      "smiles": "CC(=O)c1cc2cnc(Nc3ccc(N4CCNCC4)cn3)nc2n1C1CCCC1",
      "target": "CDK4/CDK6",
      "rationale": "Compensates for CDKN2A loss by directly inhibiting CDK4/6. MolMIM-optimized from palbociclib scaffold for improved pediatric PK.",
      "scores": {
        "diffdock_binding_affinity_kcal": -9.8,
        "admet_oral_bioavailability": 0.72,
        "admet_herg_liability": "low",
        "admet_hepatotoxicity": "low",
        "molecular_weight": 438.2,
        "lipinski_violations": 0,
        "pediatric_safety_score": 0.89,
        "druglikeness_qed": 0.78
      },
      "known_reference": "Palbociclib (FDA-approved for breast cancer; Phase I/II in pediatric malignancies — NCT03526250)",
      "pediatric_considerations": "CDK4/6 inhibitors under active investigation in pediatric solid tumors and lymphomas. Myelosuppression is dose-limiting — requires careful monitoring in ALL patients already on cytotoxic therapy."
    },
    {
      "rank": 2,
      "candidate_id": "HCLS-LEN-002",
      "name": "Lenalidomide derivative (cereblon-IKZF modulator)",
      "smiles": "O=C1CCC(N2C(=O)c3cccc(N)c3C2=O)C(=O)N1",
      "target": "CRBN-IKZF axis",
      "rationale": "Cereblon E3 ligase modulator that promotes degradation of IKZF1/3. In the context of IKZF1 deletion producing a dominant-negative isoform, IMiD-mediated degradation of the truncated protein may restore normal lymphoid differentiation.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.4,
        "admet_oral_bioavailability": 0.68,
        "admet_herg_liability": "low",
        "admet_hepatotoxicity": "moderate",
        "molecular_weight": 259.3,
        "lipinski_violations": 0,
        "pediatric_safety_score": 0.71,
        "druglikeness_qed": 0.82
      },
      "known_reference": "Lenalidomide (FDA-approved for multiple myeloma, MDS; investigated in pediatric ALL — PMID:31862848)",
      "pediatric_considerations": "CRITICAL: Lenalidomide is teratogenic (Category X). Pediatric use requires strict REMS compliance. Thrombotic risk requires prophylaxis."
    },
    {
      "rank": 3,
      "candidate_id": "HCLS-BLIN-003",
      "name": "Blinatumomab-informed bispecific concept",
      "smiles": "N/A (biologic)",
      "target": "CD19-CD3 bispecific",
      "rationale": "CD19 is highly expressed on Evelyn's blasts (flow cytometry CD19+). Blinatumomab (anti-CD19/CD3 BiTE) is FDA-approved for MRD+ BCP-ALL and relapsed/refractory pediatric ALL.",
      "scores": {
        "efficacy_evidence_level": "1A",
        "pediatric_safety_score": 0.83,
        "applicability": "immediate — FDA-approved indication"
      },
      "known_reference": "Blinatumomab (Blincyto) — FDA-approved for MRD+ BCP-ALL (2018). COG AALL1731 incorporating blinatumomab into frontline therapy.",
      "pediatric_considerations": "CRS risk ~5-15% (mostly Grade 1-2). Neurotoxicity risk ~10%. Requires continuous IV infusion (28-day cycles). Hospitalization required for cycle 1."
    }
  ],
  "pipeline_metadata": {
    "molmim_candidates_generated": 50,
    "diffdock_docked": 50,
    "rdkit_admet_scored": 50,
    "pediatric_filter_passed": 12,
    "total_processing_time_seconds": 847
  }
}
```

**Presenter note:** This step demonstrates the full power of the platform — moving from DNA to drug candidates in a single session. Highlight that candidate #3 (blinatumomab) represents an immediately actionable FDA-approved therapy, while candidates #1 and #2 represent novel optimization targets generated by MolMIM. The pediatric safety filter is critical: it eliminated 38 of 50 candidates due to hepatotoxicity, cardiotoxicity, or teratogenicity concerns specific to pediatric patients.

---

### Clinical Insights -- Demo 1

1. **Conflicting molecular signals require integrated analysis.** ETV6-RUNX1 is the most favorable marker in childhood ALL, yet IKZF1 deletion substantially worsens prognosis. The IKZF1plus phenotype — defined by IKZF1 deletion co-occurring with CDKN2A, PAX5, or other deletions without ERG deletion — was only formally recognized in 2018. Without AI-assisted genomic analysis, this interaction is easily missed.

2. **The dual entry point model enables flexible demonstrations.** Entry Point A showcases the full genomic pipeline (2.5 hours) while Entry Point B (15 minutes) focuses on clinical intelligence. The same downstream analysis applies regardless of entry point.

3. **Therapeutic discovery spans the actionability spectrum.** The three candidates represent: (a) an FDA-approved therapy for the exact indication (blinatumomab), (b) a repurposing candidate with preclinical rationale (palbociclib for CDKN2A-deleted ALL), and (c) a novel mechanism (cereblon-mediated IKZF1 degradation). This range demonstrates that AI-driven discovery complements rather than replaces standard of care.

4. **Pediatric safety filtering is not optional.** Of 50 MolMIM-generated candidates, 38 were eliminated by the pediatric safety filter. Children metabolize drugs differently than adults, and developing organs are uniquely vulnerable to off-target effects.

---

### Demo 1 Summary

| Step | Component | Key Output |
|------|-----------|------------|
| A1-A4 | Genomic Foundation Engine | Aligned BAM, annotated VCF, 24,891 embedded variants |
| B | VCF Ingestion | Direct path: 34 seconds to embedded vectors |
| 6 | Oncology Agent | Risk classification: SR with IKZF1plus — consider HR arm |
| 7 | RAG Engine | Evidence synthesis: MRD-guided therapy recommendation |
| 8 | Discovery Engine | 3 therapeutic candidates: palbociclib analog, lenalidomide derivative, blinatumomab |

**Total demo time:** Entry Point A: ~2.5 hours (mostly pipeline compute) | Entry Point B: ~15 minutes

---

---

## Section 6: DEMO 2 -- Pediatric ALL Multi-Agent Tumor Board

### Clinical Narrative

**Patient:** Evelyn R., 8-year-old female (same patient from Demo 1)
**Timepoint:** Day 29 of induction chemotherapy (COG AALL0932)
**Update:** MRD positive at 0.1% by multiparameter flow cytometry — above the 0.01% threshold for favorable response.

**Day 29 bone marrow assessment:**
- Morphology: <5% blasts (morphologic remission achieved)
- MRD by flow cytometry: 0.1% (10^-3) — **POSITIVE**
- MRD by PCR (ETV6-RUNX1 fusion transcript): Detectable (ratio 1.2 x 10^-3)
- Residual blast immunophenotype: CD19+ (MFI 12,400), CD22+ (MFI 3,800), CD10+, CD34+, TdT+

**Clinical context:** Evelyn achieved morphologic remission but her MRD level places her in an intermediate-risk category. Given her IKZF1plus molecular profile (from Demo 1), the pediatric oncology team convenes a molecular tumor board to determine optimal therapy intensification. The platform's multi-agent coordination system is engaged to provide a comprehensive assessment.

**Clinical question:** Should Evelyn be escalated to high-risk therapy? Is she a candidate for immunotherapy (blinatumomab, CAR-T)? What biomarkers guide therapy selection?

**Agents involved:** Precision Oncology Agent (:8527), Precision Biomarker Agent (:8529), CAR-T Intelligence Agent (:8522), Single-Cell Intelligence Agent (:8540), Clinical Trial Intelligence Agent (:8538)

---

### Architecture Overview

```
+---------------------------+
|   SHARED GENOMIC LAYER    |
|   (Milvus: sophia_tumor)  |
+------------+--------------+
             |
             v
+---------------------------+       +---------------------------+       +---------------------------+
|   ONCOLOGY AGENT (:8527)  | <---> |   BIOMARKER AGENT (:8529) | <---> |   CAR-T AGENT (:8522)     |
|   Orchestrator / MTB Lead |       |   CD19/CD22/MRD           |       |   Tisagenlecleucel        |
+---------------------------+       +---------------------------+       +---------------------------+
             ^                                   ^                                   ^
             |                                   |                                   |
             v                                   v                                   v
+---------------------------+       +---------------------------+       +---------------------------+
|   SINGLE-CELL (:8540)     |       |   CLINICAL TRIAL (:8538)  |       | THERAPEUTIC DISCOVERY     |
|   Blast immunophenotyping |       |   COG/NCI matching        |       | ENGINE (:8510)            |
+---------------------------+       +---------------------------+       | MolMIM→DiffDock→RDKit     |
                                                                        +---------------------------+
```

Five agents coordinate through cross-agent triggers, with the Oncology Agent serving as the molecular tumor board orchestrator.

---

### Entry Point A / Entry Point B

This demo uses the same VCF and Milvus collection from Demo 1 (`sophia_tumor_001_variants`). If running Demo 2 independently, use Entry Point B from Demo 1 to ingest the VCF first.

---

### Step 1 -- Oncology Agent: Integrated Multi-Agent Assessment

The Oncology Agent's `/integrated-assessment` endpoint orchestrates all participating agents and synthesizes their findings into a unified tumor board report.

```bash
curl -X POST http://localhost:8527/v1/oncology/integrated-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "assessment_type": "molecular_tumor_board",
    "timepoint": "day_29_induction",
    "variant_collection": "sophia_tumor_001_variants",
    "clinical_data": {
      "age_years": 8,
      "diagnosis": "BCP-ALL",
      "initial_wbc": 45000,
      "current_status": "morphologic_remission",
      "mrd_flow": 0.001,
      "mrd_pcr": 0.0012,
      "blast_immunophenotype": {
        "CD19": {"positive": true, "mfi": 12400},
        "CD22": {"positive": true, "mfi": 3800},
        "CD10": {"positive": true},
        "CD34": {"positive": true}
      },
      "molecular_profile": {
        "etv6_runx1": true,
        "ikzf1_deletion": true,
        "cdkn2a_deletion": true,
        "ikzf1_plus": true
      }
    },
    "agents_to_consult": ["biomarker", "cart", "single_cell", "clinical_trial"],
    "output_format": "tumor_board_report"
  }'
```

**Expected response (abbreviated — full response is ~200 lines):**

```json
{
  "patient_id": "sophia_tumor_001",
  "assessment_id": "MTB-2026-0323-001",
  "assessment_type": "molecular_tumor_board",
  "orchestrator": "oncology_agent",
  "agents_consulted": ["biomarker", "cart", "single_cell", "clinical_trial"],
  "overall_recommendation": {
    "risk_reclassification": "HIGH_RISK",
    "rationale": "IKZF1plus phenotype + MRD positivity (0.1%) at Day 29 meets criteria for high-risk reclassification per COG AALL1231 amendment. 5-year EFS with standard therapy estimated at <50%.",
    "primary_recommendation": "Escalate to COG AALL1231 high-risk consolidation with consideration for blinatumomab incorporation per AALL1731",
    "secondary_recommendation": "If MRD persists post-consolidation, proceed to tisagenlecleucel (CAR-T) evaluation"
  },
  "agent_reports": {
    "biomarker_summary": "See Step 2",
    "cart_summary": "See Step 3",
    "single_cell_summary": "See Step 4",
    "clinical_trial_summary": "See Step 5"
  },
  "processing_time_ms": 8742
}
```

---

### Step 2 -- Biomarker Agent: Target Antigen and MRD Assessment

```bash
curl -X POST http://localhost:8529/v1/biomarker/antigen-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "assessment_type": "immunotherapy_targets",
    "blast_immunophenotype": {
      "CD19": {"positive": true, "mfi": 12400, "percent_positive": 98.2},
      "CD22": {"positive": true, "mfi": 3800, "percent_positive": 94.7},
      "CD20": {"positive": false, "mfi": 180, "percent_positive": 12.3},
      "CD38": {"positive": true, "mfi": 8900, "percent_positive": 99.1}
    },
    "mrd_data": {
      "method": "multiparameter_flow",
      "value": 0.001,
      "sensitivity": 0.0001,
      "timepoint": "day_29"
    },
    "molecular_context": {
      "ikzf1_deletion": true,
      "etv6_runx1": true
    }
  }'
```

**Expected response:**

```json
{
  "patient_id": "sophia_tumor_001",
  "assessment_id": "BIO-2026-0323-001",
  "target_antigen_analysis": [
    {
      "antigen": "CD19",
      "expression_level": "high",
      "mfi": 12400,
      "percent_positive": 98.2,
      "immunotherapy_suitability": {
        "blinatumomab": {"suitable": true, "evidence_level": "1A", "notes": "Strong CD19 expression supports blinatumomab efficacy. CR rate 39% in R/R pediatric ALL (von Stackelberg et al., JCO 2016)."},
        "tisagenlecleucel_cart": {"suitable": true, "evidence_level": "1A", "notes": "CD19 MFI >1000 associated with durable CAR-T response. ELIANA trial: 82% CR rate in pediatric R/R ALL (Maude et al., NEJM 2018)."},
        "inotuzumab": {"suitable": false, "notes": "CD19-targeted, not CD22. See CD22 for inotuzumab."}
      }
    },
    {
      "antigen": "CD22",
      "expression_level": "moderate",
      "mfi": 3800,
      "percent_positive": 94.7,
      "immunotherapy_suitability": {
        "inotuzumab_ozogamicin": {"suitable": true, "evidence_level": "2A", "notes": "CD22+ in >90% of blasts supports inotuzumab efficacy. Pediatric data from COG AALL1621: CR/CRi rate 58%."},
        "cd22_cart": {"suitable": true, "evidence_level": "2B", "notes": "CD22 CAR-T (e.g., NCT02315612) shows activity in CD19-negative relapse. Consider as salvage if CD19 loss occurs."}
      }
    }
  ],
  "mrd_interpretation": {
    "level": 0.001,
    "classification": "MRD_positive",
    "prognostic_implication": "Day 29 MRD 10^-3 in IKZF1plus BCP-ALL is associated with 5-year EFS 40-50%. MRD is the strongest independent prognostic factor in pediatric ALL.",
    "recommended_actions": [
      "Intensify consolidation therapy",
      "Repeat MRD at end of consolidation (Day 78)",
      "If MRD persists >= 0.01% at Day 78, immunotherapy escalation indicated"
    ]
  },
  "processing_time_ms": 2134
}
```

---

### Step 3 -- CAR-T Intelligence Agent: Eligibility and Risk Assessment

```bash
curl -X POST http://localhost:8522/v1/cart/eligibility-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "assessment_type": "tisagenlecleucel_eligibility",
    "patient_data": {
      "age_years": 8,
      "weight_kg": 26,
      "bsa_m2": 0.96,
      "diagnosis": "BCP-ALL",
      "disease_status": "MRD_positive_morphologic_CR",
      "prior_therapies": ["induction_AALL0932"],
      "lines_of_therapy": 1,
      "cd19_expression": {"positive": true, "mfi": 12400, "percent_positive": 98.2},
      "organ_function": {
        "egfr": 142,
        "ast": 28,
        "alt": 32,
        "bilirubin": 0.4,
        "lvef": 0.65
      },
      "performance_status": {"lansky": 90}
    },
    "include_risk_modeling": true
  }'
```

**Expected response:**

```json
{
  "patient_id": "sophia_tumor_001",
  "assessment_id": "CART-2026-0323-001",
  "product": "tisagenlecleucel (Kymriah)",
  "eligibility": {
    "overall": "conditionally_eligible",
    "status": "Not currently indicated — reserve for MRD persistence post-consolidation or relapse",
    "fda_indication_met": false,
    "criteria_assessment": [
      {"criterion": "Age <= 25 years", "met": true, "value": "8 years"},
      {"criterion": "B-cell precursor ALL", "met": true},
      {"criterion": "CD19 expression", "met": true, "value": "98.2% positive, MFI 12400"},
      {"criterion": "Refractory or >= 2nd relapse", "met": false, "note": "Currently in first CR with MRD+ — does not meet R/R criteria"},
      {"criterion": "Adequate organ function", "met": true},
      {"criterion": "No active CNS disease", "met": true}
    ]
  },
  "crs_risk_assessment": {
    "predicted_grade": "Grade 1-2 (most likely)",
    "risk_factors": {
      "disease_burden": "low (MRD-level disease)",
      "age": "pediatric (higher CRS incidence but lower severity than adults)",
      "baseline_inflammation": "unknown (recommend IL-6, ferritin, CRP baseline)"
    },
    "crs_incidence_prediction": {
      "any_grade": 0.77,
      "grade_3_4": 0.22,
      "grade_5": 0.01
    },
    "management_protocol": "Tocilizumab (anti-IL-6R) available. Prophylactic dexamethasone NOT recommended per current guidelines. ICU bed reserved per institutional protocol."
  },
  "timing_recommendation": "If MRD persists >= 0.01% after consolidation (Day 78), proceed to tisagenlecleucel. Leukapheresis should be performed 4-6 weeks before planned infusion. Manufacturing time: 22-28 days.",
  "processing_time_ms": 3891
}
```

---

### Step 4 -- Single-Cell Intelligence Agent: Blast Immunophenotyping

```bash
curl -X POST http://localhost:8540/v1/single-cell/blast-characterization \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "analysis_type": "blast_immunophenotype_and_tme",
    "flow_data": {
      "timepoint": "day_29",
      "total_events": 500000,
      "blast_gate_events": 500,
      "blast_markers": {
        "CD19": {"mfi": 12400, "pct": 98.2},
        "CD22": {"mfi": 3800, "pct": 94.7},
        "CD10": {"mfi": 6700, "pct": 96.1},
        "CD34": {"mfi": 4200, "pct": 88.3},
        "CD20": {"mfi": 180, "pct": 12.3},
        "CD38": {"mfi": 8900, "pct": 99.1},
        "TdT": {"mfi": 5100, "pct": 91.0}
      }
    },
    "include_tme_inference": true,
    "include_antigen_loss_prediction": true
  }'
```

**Expected response:**

```json
{
  "patient_id": "sophia_tumor_001",
  "assessment_id": "SC-2026-0323-001",
  "blast_characterization": {
    "immunophenotype": "Pro-B/Common ALL (CD10+ CD19+ CD34+ TdT+)",
    "maturation_stage": "Hematogone-B2-like",
    "clonality_confidence": 0.96,
    "aberrant_markers": ["CD34 bright (aberrant persistence)", "CD38 uniformly bright"],
    "normal_b_cell_precursors_detected": true,
    "mrd_blast_distinction": "Clear separation from normal hematogones by CD34/CD38/CD10 pattern"
  },
  "tumor_microenvironment": {
    "t_cell_infiltrate": {
      "cd3_pct_of_marrow": 18.4,
      "cd4_cd8_ratio": 1.3,
      "exhaustion_markers": {
        "pd1_on_cd8": 0.34,
        "tim3_on_cd8": 0.12,
        "lag3_on_cd8": 0.08
      },
      "interpretation": "Moderate T-cell infiltrate with early exhaustion phenotype (PD-1+ on 34% of CD8+ T cells). This may affect CAR-T efficacy — consider lymphodepletion optimization."
    },
    "nk_cell_pct": 4.2,
    "monocyte_pct": 8.7,
    "tme_classification": "Immune-moderate, early exhaustion"
  },
  "antigen_loss_prediction": {
    "cd19_loss_risk": "low",
    "rationale": "High uniform CD19 expression (MFI 12,400, CV 0.18) with no detectable CD19-dim subpopulation. However, IKZF1 deletion may predispose to lineage switch — monitor at each assessment.",
    "cd22_loss_risk": "moderate",
    "cd22_rationale": "CD22 expression heterogeneous (CV 0.42) with a dim subpopulation (~5%). If CD22-targeted therapy is considered, repeat assessment recommended.",
    "lineage_switch_risk": "low_moderate",
    "lineage_switch_rationale": "ETV6-RUNX1+ ALL has reported rare cases of myeloid lineage switch at relapse. IKZF1 deletion may increase this risk marginally."
  },
  "processing_time_ms": 2567
}
```

---

### Step 5 -- Clinical Trial Intelligence Agent: Trial Matching

```bash
curl -X POST http://localhost:8538/v1/trials/match \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "sophia_tumor_001",
    "matching_criteria": {
      "age_years": 8,
      "diagnosis": "BCP-ALL",
      "molecular_features": ["ETV6-RUNX1", "IKZF1_deletion", "CDKN2A_deletion"],
      "disease_status": "MRD_positive_CR1",
      "mrd_level": 0.001,
      "prior_therapies": ["induction_chemotherapy"],
      "biomarkers": {
        "cd19_positive": true,
        "cd22_positive": true
      }
    },
    "trial_databases": ["clinicaltrials_gov", "cog_active", "nci_pediatric"],
    "max_results": 10,
    "include_eligibility_detail": true
  }'
```

**Expected response:**

```json
{
  "patient_id": "sophia_tumor_001",
  "assessment_id": "TRIAL-2026-0323-001",
  "matched_trials": [
    {
      "rank": 1,
      "trial_id": "NCT03914625",
      "short_name": "COG AALL1732",
      "title": "Total Therapy Study XVII for Newly Diagnosed Patients With Standard Risk ALL",
      "phase": "III",
      "relevance_score": 0.97,
      "status": "recruiting",
      "match_rationale": "First-line therapy for SR BCP-ALL with molecular risk stratification. Incorporates IKZF1 status for risk-directed therapy intensification.",
      "key_arms": [
        "Arm A: Standard consolidation + maintenance",
        "Arm B: Intensified consolidation with blinatumomab cycle for MRD+ patients"
      ],
      "eligibility_assessment": {
        "age": {"met": true, "criterion": "1-30.99 years"},
        "diagnosis": {"met": true, "criterion": "B-ALL in first remission"},
        "mrd_status": {"met": true, "criterion": "MRD assessment at Day 29 required"},
        "overall": "ELIGIBLE"
      },
      "sites_nearby": ["St. Jude Children's Research Hospital", "Children's Hospital of Philadelphia"]
    },
    {
      "rank": 2,
      "trial_id": "NCT03876769",
      "short_name": "COG AALL1731",
      "title": "Risk-Stratified Therapy for Newly Diagnosed ALL Incorporating Blinatumomab",
      "phase": "III",
      "relevance_score": 0.95,
      "status": "recruiting",
      "match_rationale": "Incorporates blinatumomab into frontline therapy for NCI SR ALL with high-risk molecular features. IKZF1plus patients eligible for blinatumomab intensification arm.",
      "eligibility_assessment": {
        "overall": "ELIGIBLE — IKZF1plus status qualifies for blinatumomab arm"
      }
    },
    {
      "rank": 3,
      "trial_id": "NCT03150693",
      "short_name": "Pediatric MATCH",
      "title": "NCI-COG Pediatric MATCH: Targeted Therapy Directed by Genetic Testing",
      "phase": "II",
      "relevance_score": 0.78,
      "status": "recruiting",
      "match_rationale": "CDKN2A homozygous deletion qualifies for Arm B (palbociclib, CDK4/6 inhibitor). Would require disease progression or relapse for eligibility.",
      "eligibility_assessment": {
        "overall": "NOT_YET_ELIGIBLE — requires relapsed/refractory disease"
      }
    },
    {
      "rank": 4,
      "trial_id": "NCT02435849",
      "short_name": "ELIANA Follow-up",
      "title": "Long-term Follow-up of Tisagenlecleucel in Pediatric ALL",
      "phase": "II",
      "relevance_score": 0.65,
      "status": "active_not_recruiting",
      "match_rationale": "Tisagenlecleucel for R/R pediatric ALL. Would require second relapse or refractory disease for eligibility.",
      "eligibility_assessment": {
        "overall": "NOT_YET_ELIGIBLE — reserve for R/R setting"
      }
    }
  ],
  "total_matches": 8,
  "processing_time_ms": 4521
}
```

---

### Step 6 -- Therapeutic Discovery Engine: Novel CREBBP Inhibitors

With the multi-agent tumor board identifying CREBBP HAT domain alterations as a key resistance mechanism in Evelyn's MRD+ ALL, we now submit to the Therapeutic Discovery Engine to generate novel small-molecule candidates targeting this vulnerability. CREBBP mutations are found in ~18% of relapsed ALL and drive glucocorticoid resistance through impaired histone acetylation — making the HAT domain a compelling therapeutic target.

```bash
curl -X POST http://localhost:8510/api/v1/discover \
  -H "Content-Type: application/json" \
  -d '{
    "target_gene": "CREBBP",
    "target_domain": "HAT_domain",
    "target_structure_pdb": "5BN4",
    "patient_id": "sophia_tumor_001",
    "optimization_constraints": {
      "pediatric_safe": true,
      "bbb_penetration": true,
      "max_molecular_weight": 450,
      "herg_liability": "none",
      "age_group": "pediatric_6_12",
      "exclude_known_cardiotoxins": true
    },
    "clinical_context": "MRD+ ALL with CREBBP-mediated glucocorticoid resistance. CNS relapse risk requires BBB-penetrant compounds.",
    "num_candidates": 100,
    "top_k": 3
  }'
```

**Expected response:**

```json
{
  "discovery_run_id": "DISC-2026-0323-CREBBP-001",
  "target": "CREBBP HAT domain",
  "pdb_structure": "5BN4",
  "status": "completed",
  "candidates_generated": 100,
  "candidates_passing_filters": 14,
  "top_candidates": [
    {
      "rank": 1,
      "candidate_id": "HCLS-CREBBP-001",
      "name": "HAT domain inhibitor (MolMIM-generated, A-485 scaffold optimization)",
      "smiles": "CC1=CC(=C(C=C1)NC(=O)C2=CC=C(C=C2)S(=O)(=O)N3CCCC3)OC",
      "target": "CREBBP HAT domain",
      "rationale": "Optimized from A-485 scaffold with improved BBB penetration. Selectively inhibits CREBBP/p300 HAT activity to restore glucocorticoid sensitivity in resistant ALL blasts.",
      "scores": {
        "diffdock_binding_affinity_kcal": -9.1,
        "bbb_penetration_score": 0.82,
        "molecular_weight": 388.5,
        "herg_liability": "none",
        "hepatotoxicity_risk": "low",
        "pediatric_safety_score": 0.87,
        "druglikeness_qed": 0.81
      },
      "pediatric_considerations": "BBB-penetrant design addresses CNS relapse risk. No hERG liability detected. Predicted half-life compatible with twice-daily pediatric dosing."
    },
    {
      "rank": 2,
      "candidate_id": "HCLS-CREBBP-002",
      "name": "Bromodomain-HAT dual inhibitor (MolMIM novel scaffold)",
      "smiles": "O=C(NC1=CC=C(F)C=C1)C2=CN=C(N)C3=CC=CC=C23",
      "target": "CREBBP bromodomain + HAT domain",
      "rationale": "Dual-mechanism compound targeting both acetyl-lysine reading (bromodomain) and writing (HAT) functions of CREBBP. Enhanced efficacy through simultaneous disruption of both catalytic and reader functions.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.7,
        "bbb_penetration_score": 0.76,
        "molecular_weight": 310.3,
        "herg_liability": "none",
        "hepatotoxicity_risk": "low",
        "pediatric_safety_score": 0.84,
        "druglikeness_qed": 0.85
      },
      "pediatric_considerations": "Low MW (310 Da) favorable for pediatric PK. Dual-target mechanism may reduce required dose. No CYP3A4 inhibition — minimal drug-drug interaction with maintenance chemotherapy."
    },
    {
      "rank": 3,
      "candidate_id": "HCLS-CREBBP-003",
      "name": "Allosteric HAT modulator (MolMIM-generated)",
      "smiles": "CC(C)NC(=O)C1=CC=C(C=C1)N2C=C(C=N2)C3=CC=NC=C3",
      "target": "CREBBP HAT allosteric site",
      "rationale": "Allosteric mechanism avoids competition with acetyl-CoA substrate, potentially reducing off-target effects on other acetyltransferases. Designed for selectivity over p300.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.2,
        "bbb_penetration_score": 0.88,
        "molecular_weight": 333.4,
        "herg_liability": "none",
        "hepatotoxicity_risk": "minimal",
        "pediatric_safety_score": 0.91,
        "druglikeness_qed": 0.79
      },
      "pediatric_considerations": "Highest pediatric safety score of the three candidates. Excellent BBB penetration (0.88) for CNS relapse prevention. Allosteric mechanism provides wider therapeutic window."
    }
  ],
  "pipeline_metadata": {
    "molmim_candidates_generated": 100,
    "diffdock_docked_against": "5BN4",
    "rdkit_admet_scored": 100,
    "pediatric_filter_passed": 14,
    "bbb_filter_passed": 14,
    "herg_filter_eliminated": 31,
    "total_processing_time_seconds": 1042
  }
}
```

**Clinical interpretation:** The Therapeutic Discovery Engine identifies three novel CREBBP HAT domain inhibitors optimized for Evelyn's specific clinical scenario. All three candidates pass the BBB penetration filter — critical because Evelyn's MRD positivity and IKZF1plus profile confer elevated CNS relapse risk. Candidate #3 achieves the highest pediatric safety score (0.91) through its allosteric mechanism, while Candidate #1 offers the strongest binding affinity (-9.1 kcal/mol). These represent preclinical-stage leads that could inform the selection of investigational agents in clinical trials, or serve as starting points for medicinal chemistry optimization. The 86% candidate elimination rate (86 of 100 filtered out) demonstrates the stringency of the pediatric safety filter, particularly the hERG screen which eliminated 31 candidates with potential cardiac liability.

**Presenter note:** This step extends Demo 2 from diagnosis and risk stratification into actionable therapeutic discovery. Highlight that the BBB penetration constraint is directly motivated by Evelyn's molecular risk profile — the same multi-agent analysis that identified her CNS relapse risk now drives the drug design constraints. This is precision medicine in the fullest sense: from variant to vulnerability to candidate molecule, with every constraint grounded in the patient's specific biology.

---

### Clinical Insights -- Demo 2

1. **Multi-agent coordination compresses tumor board preparation from days to minutes.** In a traditional workflow, the oncologist would consult hematopathology for MRD interpretation, immunology for antigen assessment, cell therapy for CAR-T eligibility, and clinical trials for matching — each consultation requiring separate referrals. The platform provides all five assessments in a single orchestrated call.

2. **MRD drives decision-making.** The 0.1% MRD level at Day 29, combined with IKZF1plus molecular features, triggers risk reclassification from standard to high risk. This is a concrete example of how molecular and response data interact — neither alone tells the full story.

3. **Immunotherapy readiness assessment is proactive.** Even though Evelyn does not currently qualify for CAR-T (she is in first remission), the CAR-T agent provides a prospective eligibility assessment and CRS risk model. If MRD persists after consolidation, the team already has a treatment plan.

4. **Antigen loss prediction informs therapy sequencing.** The Single-Cell agent identifies low CD19 loss risk but moderate CD22 heterogeneity, which influences the sequencing decision: blinatumomab (anti-CD19) first, with CD22-targeted therapy (inotuzumab or CD22 CAR-T) held in reserve.

---

### Demo 2 Summary

| Step | Agent | Key Output |
|------|-------|------------|
| 1 | Oncology (orchestrator) | Risk reclassification: HIGH RISK. Escalate therapy. |
| 2 | Biomarker | CD19 high (MFI 12,400) — blinatumomab/CAR-T suitable. MRD 0.1% — intensify. |
| 3 | CAR-T | Conditionally eligible for tisagenlecleucel. Reserve for MRD persistence post-consolidation. CRS risk: Grade 1-2 most likely. |
| 4 | Single-Cell | Blast phenotype confirmed. Early T-cell exhaustion (PD-1+ 34% of CD8+). Low CD19 loss risk. |
| 5 | Clinical Trial | Top match: COG AALL1732 (eligible). AALL1731 blinatumomab arm (eligible). Pediatric MATCH palbociclib arm (if relapse). |
| 6 | Therapeutic Discovery Engine | 3 novel CREBBP HAT domain inhibitors. Top candidate: -9.1 kcal/mol, BBB-penetrant, pediatric safety 0.87. 86% elimination rate by pediatric/hERG filters. |

**Total demo time:** ~25 minutes (all agents query the same Milvus collection from Demo 1, plus ~15 min for discovery pipeline)

**Presenter note:** This demo showcases the platform's most distinctive capability — multi-agent coordination. Five specialized agents share the same genomic evidence layer and produce complementary assessments. In a traditional tumor board, this information would require consultations across hematology, immunology, cell therapy, and clinical trials teams spanning days to weeks. The platform compresses this to minutes.

---

---

## Section 7: DEMO 3 -- Cardiotoxicity Prevention in Pediatric Chemotherapy

### Clinical Narrative

**Patient:** Marcus J., 6-year-old male
**Presenting complaint:** Abdominal distension, weight loss (3 kg over 6 weeks), intermittent fevers, and a palpable left flank mass. Parents report irritability and decreased appetite.

**Initial workup:**
- Abdominal CT: 8.2 cm left adrenal mass with calcifications, encasing the aorta without vascular invasion
- Urinary catecholamines: VMA 42.3 mg/24hr (elevated), HVA 58.1 mg/24hr (elevated)
- Serum LDH: 1,890 U/L (elevated)
- NSE: 312 ng/mL (markedly elevated, poor prognostic marker)
- Ferritin: 245 ng/mL

**Bone marrow biopsy:**
- Bilateral biopsies: Left iliac crest positive for neuroblastoma (small round blue cells, neuropil)
- Immunohistochemistry: Synaptophysin+, Chromogranin A+, CD56+, NB84+
- MYCN amplification: **POSITIVE** (>100 copies by FISH)
- 1p36 deletion: Present
- 11q23 deletion: Absent
- Ploidy: Diploid (unfavorable)

**Staging:** INSS Stage 4 (adrenal primary + bone marrow metastases)
**Risk classification:** High-risk neuroblastoma (age >18 months, MYCN amplified, stage 4, unfavorable histology)

**Treatment plan:** COG ANBL1232 High-Risk Protocol
- Induction: 5 cycles alternating cyclophosphamide/topotecan and cisplatin/etoposide with daunorubicin (cumulative anthracycline target: 300 mg/m2 doxorubicin equivalents)
- Consolidation: High-dose chemotherapy with autologous stem cell rescue
- Maintenance: Isotretinoin + anti-GD2 immunotherapy (dinutuximab), 6-mercaptopurine

**Agents involved:** PGx Agent (:8107), Biomarker Agent (:8529), Cardiology Agent (:8126), Neurology Agent (:8528), Autoimmune Agent (:8532), Oncology Agent (:8527)

**Clinical question:** What pharmacogenomic factors influence drug metabolism and toxicity? What is the cardiotoxicity risk given the planned anthracycline exposure, and how should cardiac function be monitored?

---

### Architecture Overview

```
+-------------------------------+       +-------------------------------+
|   GENOMIC FOUNDATION ENGINE   |       |   PRECISION INTELLIGENCE      |
|                               |       |   NETWORK                     |
|  FASTQ --> BWA-MEM2 --> BAM   |       |                               |
|  BAM --> DeepVariant --> VCF   | ----> |  PGx Agent (:8107)            |
|  VCF --> ClinVar/AlphaMiss    |       |  Biomarker Agent (:8529)      |
|  Annotated VCF --> Milvus     |       |  Cardiology Agent (:8126)     |
+-------------------------------+       |  Neurology Agent (:8528)      |
                                        |  Autoimmune Agent (:8532)     |
                                        |  Oncology Agent (:8527)       |
                                        |  → Therapeutic Discovery      |
                                        |    Engine (:8510)             |
                                        +-------------------------------+
```

---

### Entry Point A: Starting from FASTQ

#### Step A1 -- Upload FASTQ and Run Germline Pipeline

Marcus's germline sample (buccal swab) is processed through the Genomic Foundation Engine for pharmacogenomic variant calling.

```bash
mkdir -p /data/genomics/inputs/marcus_germline
cp marcus_germline_R1.fastq.gz /data/genomics/inputs/marcus_germline/
cp marcus_germline_R2.fastq.gz /data/genomics/inputs/marcus_germline/

curl -X POST http://localhost:5000/v1/genomics/run \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "marcus_germline",
    "pipeline_type": "germline_only",
    "input_fastq_r1": "/data/genomics/inputs/marcus_germline/marcus_germline_R1.fastq.gz",
    "input_fastq_r2": "/data/genomics/inputs/marcus_germline/marcus_germline_R2.fastq.gz",
    "reference_genome": "GRCh38",
    "target_regions": "pharmacogenomic_plus_cardiotoxicity",
    "caller": "deepvariant",
    "output_dir": "/data/genomics/outputs/marcus_germline"
  }'
```

**Expected response:**

```json
{
  "run_id": "GEN-2026-0323-003",
  "sample_id": "marcus_germline",
  "status": "running",
  "pipeline_type": "germline_only",
  "estimated_completion_minutes": 45,
  "output_vcf": "/data/genomics/outputs/marcus_germline/marcus_germline.deepvariant.g.vcf.gz"
}
```

Pipeline completes in ~45 minutes for germline-only at 30x coverage.

---

### Entry Point B: Starting from VCF

```bash
curl -X POST http://localhost:5000/v1/genomics/ingest-vcf \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "marcus_germline",
    "vcf_path": "/data/genomics/external/marcus_germline.vcf.gz",
    "reference_genome": "GRCh38",
    "annotate": true,
    "embed": true,
    "annotation_sources": ["clinvar", "alphamissense", "pharmgkb", "cpic"]
  }'
```

---

### Step 1 -- PGx Agent: Pharmacogenomic Assessment

Query the PGx Agent for star-allele calling and drug interaction analysis across Marcus's full chemotherapy regimen.

```bash
curl -X POST http://localhost:8107/v1/pgx/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "marcus_germline",
    "query": "Pharmacogenomic assessment for high-risk neuroblastoma treatment: daunorubicin, cyclophosphamide, topotecan, cisplatin, etoposide, vincristine, 6-mercaptopurine maintenance.",
    "drugs_of_interest": [
      "daunorubicin", "cyclophosphamide", "topotecan", "cisplatin",
      "etoposide", "vincristine", "6-mercaptopurine", "isotretinoin"
    ],
    "genes_of_interest": ["TPMT", "NUDT15", "CYP3A5", "UGT1A1", "DPYD", "CYP2B6", "RARG", "CBR3"],
    "include_cpic_guidelines": true,
    "include_star_alleles": true
  }'
```

**Expected response:**

```json
{
  "query_id": "PGX-2026-0323-003",
  "patient_id": "marcus_germline",
  "pharmacogenomic_profile": {
    "TPMT": {
      "diplotype": "*1/*3A",
      "phenotype": "Intermediate Metabolizer",
      "activity_score": 1.0,
      "cpic_level": "A",
      "clinical_impact": {
        "drug": "6-mercaptopurine",
        "recommendation": "Reduce starting dose to 50% of standard (25 mg/m2/day instead of 50 mg/m2/day). Titrate based on ANC. Monitor weekly CBC during first 8 weeks.",
        "evidence": "CPIC guideline for thiopurines (Relling et al., Clin Pharmacol Ther 2019). TPMT intermediate metabolizers have 3-5x increased risk of myelosuppression at standard doses."
      }
    },
    "NUDT15": {
      "diplotype": "*1/*1",
      "phenotype": "Normal Metabolizer",
      "activity_score": 2.0,
      "clinical_impact": {
        "drug": "6-mercaptopurine",
        "recommendation": "No NUDT15-based dose adjustment required."
      }
    },
    "CYP3A5": {
      "diplotype": "*1/*3",
      "phenotype": "Intermediate Metabolizer",
      "activity_score": 1.0,
      "clinical_impact": {
        "drug": "vincristine",
        "recommendation": "CYP3A5 expressors (*1 allele present) metabolize vincristine more rapidly, potentially reducing efficacy but also neuropathy risk. Standard dosing with close neuropathy monitoring. Vincristine dose cap: 2 mg absolute.",
        "evidence": "Egbelakin et al., Pediatr Blood Cancer 2011. CYP3A5*1 carriers have 40% lower vincristine AUC."
      }
    },
    "CBR3": {
      "diplotype": "V244M G/G (Val/Val)",
      "phenotype": "High CBR3 Activity",
      "clinical_impact": {
        "drug": "daunorubicin",
        "recommendation": "CBR3 Val/Val genotype associated with increased conversion of daunorubicin to cardiotoxic daunorubicinol. Enhanced cardiac monitoring. Dexrazoxane cardioprotection if cumulative dose exceeds 250 mg/m2.",
        "evidence": "Blanco et al., J Clin Oncol 2012. CBR3 Val/Val associated with 5.5x increased risk of cardiomyopathy in children."
      }
    },
    "RARG": {
      "diplotype": "rs2229774 C/T (heterozygous)",
      "phenotype": "Increased Cardiotoxicity Risk",
      "clinical_impact": {
        "drug": "daunorubicin",
        "recommendation": "RARG rs2229774 variant associated with 4.7-fold increased risk of anthracycline-induced cardiomyopathy. Strongly recommend dexrazoxane co-administration and intensive echocardiographic surveillance.",
        "evidence": "Aminkeng et al., Nat Genet 2015. CPNDS Clinical Recommendation."
      }
    },
    "UGT1A1": {
      "diplotype": "*1/*1",
      "phenotype": "Normal Metabolizer",
      "clinical_impact": {
        "drug": "etoposide",
        "recommendation": "No UGT1A1-based dose adjustment required."
      }
    }
  },
  "integrated_risk_summary": "HIGH CARDIOTOXICITY RISK: CBR3 Val/Val + RARG rs2229774 create compounding anthracycline risk. Mandatory dexrazoxane from first dose. THIOPURINE DOSE REDUCTION: TPMT *1/*3A requires 50% 6-MP reduction. VINCRISTINE: CYP3A5*1 carrier, standard dosing with neuropathy vigilance.",
  "processing_time_sec": 2.8
}
```

---

### Step 2 -- Precision Biomarker Agent: Neuroblastoma Prognostic Panel

Clinical rationale: Neuroblastoma risk stratification requires specific biomarkers beyond genomics. Serum LDH, ferritin, urinary catecholamines (VMA/HVA), and MYCN copy number quantification are critical for treatment intensity decisions.

```bash
curl -X POST http://localhost:8529/v1/biomarker/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Neuroblastoma prognostic biomarker panel for 6-year-old with MYCN-amplified stage 4 disease",
    "patient_id": "PEDS-NB-002",
    "collections_filter": ["biomarker_reference", "biomarker_clinical_evidence"]
  }'
```

**Expected response:**

```json
{
  "query_id": "BIO-2026-0323-003",
  "patient_id": "PEDS-NB-002",
  "biomarker_panel": {
    "MYCN_amplification": {
      "value": ">100 copies by FISH",
      "interpretation": "HIGH RISK — MYCN amplification is the single strongest adverse prognostic factor in neuroblastoma",
      "evidence": "Brodeur et al., J Clin Oncol 1984; COG/INRG risk classification"
    },
    "serum_LDH": {
      "value": "1,450 U/L",
      "interpretation": "Elevated — poor prognosis. LDH reflects tumor burden and cellular turnover",
      "reference_range": "< 500 U/L (age-adjusted pediatric)"
    },
    "serum_ferritin": {
      "value": "280 ng/mL",
      "interpretation": "Elevated — poor prognosis. Ferritin > 92 ng/mL associated with worse outcomes in neuroblastoma",
      "evidence": "Hann et al., Cancer 1985"
    },
    "urinary_catecholamines": {
      "VMA": "42.3 mg/24hr (elevated)",
      "HVA": "58.1 mg/24hr (elevated)",
      "VMA_HVA_ratio": "Abnormal — catecholamine-secreting tumor confirmed",
      "interpretation": "Elevated VMA/HVA confirms catecholamine-producing neuroblastoma; useful for disease monitoring during treatment"
    },
    "chromosomal_aberrations": {
      "1p_deletion": {
        "status": "Present",
        "interpretation": "Additional poor prognosis marker — 1p36 loss of heterozygosity associated with MYCN amplification and aggressive disease"
      },
      "11q_deletion": {
        "status": "Absent",
        "interpretation": "11q deletion typically mutually exclusive with MYCN amplification; absence is consistent with MYCN-amplified biology"
      }
    },
    "DNA_ploidy": {
      "value": "Diploid",
      "interpretation": "Unfavorable in children >18 months. Hyperdiploid tumors have better prognosis in infants; diploid tumors are associated with MYCN amplification and poor outcomes",
      "evidence": "Look et al., J Clin Oncol 1991"
    },
    "risk_classification": {
      "system": "COG/INRG",
      "classification": "HIGH RISK",
      "criteria": "Age >18 months, MYCN amplified, Stage 4, unfavorable histology, diploid, 1p deletion",
      "expected_EFS": "40-50% with current intensive multimodal therapy"
    }
  },
  "cross_agent_flags": [
    "Cardiology Agent: Elevated LDH may indicate baseline cardiac stress — correlate with troponin/BNP",
    "Oncology Agent: All prognostic biomarkers confirm HIGH RISK — intensive multimodal therapy per COG ANBL1232",
    "PGx Agent: MYCN amplification status may influence drug sensitivity profiles"
  ],
  "processing_time_sec": 2.1
}
```

**Clinical interpretation:** The Biomarker Agent aggregates laboratory, cytogenetic, and molecular prognostic markers to confirm Marcus's high-risk classification. Every biomarker aligns with poor prognosis: MYCN amplification (>100 copies), elevated LDH and ferritin, abnormal catecholamines, 1p deletion, and diploid DNA content. This comprehensive biomarker profile drives the treatment intensity decision — Marcus requires maximum-intensity induction, consolidation with stem cell rescue, and maintenance immunotherapy. The cross-agent flags ensure downstream agents (particularly Cardiology and Oncology) incorporate biomarker-derived risk into their assessments.

---

### Step 3 -- Cardiology Agent: Pediatric Cardiotoxicity Assessment

With PGx findings revealing high cardiotoxicity risk (CBR3 Val/Val + RARG variant), query the Cardiology Agent for comprehensive cardiac monitoring and prevention.

```bash
curl -X POST http://localhost:8126/v1/cardio/integrated-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "marcus_germline",
    "assessment_type": "pediatric_cardiotoxicity_assessment",
    "clinical_context": {
      "diagnosis": "High-risk neuroblastoma, MYCN amplified, Stage 4",
      "age_years": 6,
      "sex": "male",
      "weight_kg": 21,
      "bsa_m2": 0.82,
      "planned_anthracycline": {
        "drug": "daunorubicin",
        "planned_cumulative_dose_mg_m2": 300,
        "doxorubicin_equivalents_mg_m2": 300,
        "schedule": "5 cycles, 60 mg/m2 per cycle"
      },
      "pharmacogenomic_risk": {
        "CBR3": "Val/Val (high metabolizer, 5.5x cardiomyopathy risk)",
        "RARG": "rs2229774 C/T (4.7x cardiomyopathy risk)",
        "combined_genetic_risk": "VERY HIGH"
      },
      "baseline_echocardiogram": {
        "lvef_pct": 65,
        "shortening_fraction_pct": 36,
        "gls_pct": -21.2,
        "troponin_i_ng_ml": 0.01,
        "bnp_pg_ml": 18,
        "date": "2026-03-20"
      }
    },
    "include_monitoring_schedule": true,
    "include_cardioprotection_plan": true
  }'
```

**Expected response:**

```json
{
  "assessment_id": "CARDIO-2026-0323-003",
  "patient_id": "marcus_germline",
  "cardiac_risk_assessment": {
    "overall_risk_tier": "VERY HIGH",
    "risk_score": 8.7,
    "risk_factors": [
      {"factor": "Planned cumulative anthracycline", "value": "300 mg/m2", "risk_contribution": "HIGH -- exceeds 250 mg/m2 pediatric threshold"},
      {"factor": "CBR3 Val/Val genotype", "value": "Homozygous high-activity", "risk_contribution": "HIGH -- 5.5x relative risk (Blanco et al. 2012)"},
      {"factor": "RARG rs2229774", "value": "Heterozygous carrier", "risk_contribution": "HIGH -- 4.7x relative risk (Aminkeng et al. 2015)"},
      {"factor": "Age at exposure", "value": "6 years", "risk_contribution": "MODERATE -- younger age increases long-term risk (Lipshultz et al. 2013)"}
    ],
    "estimated_cardiomyopathy_probability": "35-45% at 10-year follow-up without cardioprotection; 10-15% with dexrazoxane"
  },
  "cardioprotection_plan": {
    "dexrazoxane": {
      "recommendation": "MANDATORY from first anthracycline dose",
      "dose_ratio": "10:1 dexrazoxane:daunorubicin (600 mg/m2 per 60 mg/m2 daunorubicin)",
      "timing": "Administer 30 minutes before each daunorubicin dose",
      "evidence": "COG ALTE11C2: dexrazoxane reduced cardiac events by 68% in children (Asselin et al., Lancet Oncol 2016)"
    },
    "medication": "ACE inhibitor (enalapril 0.1 mg/kg/day) if LVEF drops below 55% or GLS worsens >15% from baseline"
  },
  "monitoring_schedule": [
    {"timepoint": "Baseline (pre-cycle 1)", "assessment": "Echo (LVEF, SF, GLS), troponin I, BNP, ECG", "status": "COMPLETE"},
    {"timepoint": "Pre-cycle 2 (cumulative 60 mg/m2)", "assessment": "Troponin I, BNP. Echo if biomarkers elevated."},
    {"timepoint": "Pre-cycle 3 (cumulative 120 mg/m2)", "assessment": "Echo (LVEF, SF, GLS), troponin I, BNP"},
    {"timepoint": "Pre-cycle 4 (cumulative 180 mg/m2)", "assessment": "Troponin I, BNP. Echo if biomarkers elevated."},
    {"timepoint": "Pre-cycle 5 (cumulative 240 mg/m2)", "assessment": "Echo (LVEF, SF, GLS), troponin I, BNP, ECG"},
    {"timepoint": "End of induction (cumulative 300 mg/m2)", "assessment": "Full cardiac: Echo, troponin I, BNP, ECG, cardiac MRI if LVEF <55%"},
    {"timepoint": "Pre-consolidation HDC", "assessment": "Echo clearance required, LVEF >= 50%"},
    {"timepoint": "6 months post-treatment", "assessment": "Echo, BNP"},
    {"timepoint": "Annual for 10 years", "assessment": "Echo, BNP, ECG. Cardiac MRI at 2 and 5 years."},
    {"timepoint": "Lifelong", "assessment": "Echo every 2-5 years. Annual monitoring recommended through age 30 given genetic risk."}
  ],
  "hold_criteria": {
    "hold_anthracycline_if": [
      "LVEF < 50% or decline > 10 percentage points from baseline",
      "Shortening fraction < 28%",
      "GLS worsening > 15% relative from baseline",
      "Troponin I > 0.1 ng/mL sustained",
      "Clinical heart failure symptoms"
    ],
    "action_if_held": "Cardiology consultation within 24 hours. Consider anthracycline-free substitution (topotecan/cyclophosphamide blocks)."
  },
  "processing_time_sec": 3.1
}
```

---

### Step 4 -- Neurology Agent: Vincristine Neuropathy Risk

```bash
curl -X POST http://localhost:8528/v1/neuro/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "marcus_germline",
    "query": "Vincristine-induced peripheral neuropathy risk in 6-year-old male with CYP3A5 *1/*3 receiving vincristine for high-risk neuroblastoma. Monitoring plan and dose modification criteria.",
    "clinical_context": {
      "diagnosis": "High-risk neuroblastoma",
      "age_years": 6,
      "vincristine_dose": "1.5 mg/m2 weekly during induction blocks",
      "vincristine_cap": "2 mg absolute",
      "cyp3a5_genotype": "*1/*3 (intermediate metabolizer)",
      "baseline_neuro_exam": "Normal motor and sensory exam. DTRs 2+ symmetric. No gait abnormality."
    }
  }'
```

**Expected response:**

```json
{
  "query_id": "NEURO-2026-0323-003",
  "patient_id": "marcus_germline",
  "neuropathy_risk_assessment": {
    "overall_risk": "MODERATE (30-40%)",
    "risk_factors": [
      "Multiple vincristine-containing cycles increase cumulative risk",
      "CYP3A5 *1/*3 intermediate metabolizer -- partial clearance advantage"
    ],
    "protective_factors": [
      "CYP3A5 *1 allele present -- faster vincristine metabolism vs *3/*3",
      "Young age associated with better nerve regeneration capacity"
    ]
  },
  "clinical_recommendations": {
    "dosing": "Standard 1.5 mg/m2 (max 2 mg absolute). No preemptive dose reduction warranted.",
    "monitoring": [
      "Baseline: Quantitative sensory testing, grip strength, ankle dorsiflexion",
      "Before each vincristine dose: Focused neuro exam (DTRs, distal sensation, foot dorsiflexion)",
      "Monthly during treatment: Pediatric-modified Total Neuropathy Score (Ped-mTNS)",
      "Watch for: Jaw pain (early sign), constipation (autonomic), ptosis (cranial nerve), foot drop (severe)"
    ],
    "dose_modification_criteria": [
      {"grade": "Grade 1 (paresthesias, reduced DTRs)", "action": "Continue full dose. Increase monitoring."},
      {"grade": "Grade 2 (functional impairment, moderate weakness)", "action": "Reduce to 1.0 mg/m2. Consider gabapentin."},
      {"grade": "Grade 3 (severe weakness, significant limitation)", "action": "Hold vincristine. Resume at 1.0 mg/m2 after recovery to grade 1."},
      {"grade": "Grade 4 (paralysis)", "action": "Permanently discontinue vincristine."}
    ],
    "supportive_care": "Gabapentin 5-10 mg/kg/day for neuropathic pain. Prophylactic stool softeners (autonomic neuropathy causes constipation in >50% of children). Physical therapy referral if grade 2+.",
    "evidence": "Lavoie Smith et al., Cancer 2015. Vincristine neuropathy in 30-45% of pediatric patients. CYP3A5 modifies risk (Egbelakin et al., Pediatr Blood Cancer 2011)."
  },
  "processing_time_sec": 2.4
}
```

---

### Step 5 -- Precision Autoimmune Agent: Immunotherapy Adverse Event Profiling

Clinical rationale: High-risk neuroblastoma treatment includes anti-GD2 immunotherapy (dinutuximab) in consolidation/maintenance. This monoclonal antibody triggers complement-mediated and antibody-dependent cellular cytotoxicity, causing significant immune-related adverse events including neuropathic pain, capillary leak, and hypersensitivity reactions.

```bash
curl -X POST http://localhost:8532/v1/autoimmune/integrated-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_profile": {
      "patient_id": "PEDS-NB-002",
      "age": 6,
      "diagnosis": "High-risk neuroblastoma",
      "planned_immunotherapy": "dinutuximab",
      "prior_autoimmune_history": false
    }
  }'
```

**Expected response:**

```json
{
  "assessment_id": "AUTO-2026-0323-003",
  "patient_id": "PEDS-NB-002",
  "immunotherapy_agent": "dinutuximab (anti-GD2 monoclonal antibody)",
  "mechanism": "Targets GD2 ganglioside on neuroblastoma cells. Induces complement-dependent cytotoxicity (CDC) and antibody-dependent cellular cytotoxicity (ADCC). GD2 is also expressed on peripheral nerves and some normal tissues, driving on-target/off-tumor toxicity.",
  "irAE_risk_assessment": {
    "neuropathic_pain": {
      "incidence": "85%",
      "severity": "Grade 3-4 in 50% of patients",
      "mechanism": "Anti-GD2 binds GD2 on peripheral nerve fibers, activating complement cascade on nerve sheaths",
      "management": "Morphine infusion protocol: 0.02-0.05 mg/kg/hr continuous IV during dinutuximab infusion and 2 hours post. Gabapentin 5-10 mg/kg/day starting 3 days prior to infusion cycle.",
      "monitoring": "Pain assessment q1h during infusion using FLACC scale (age-appropriate for 6yo)"
    },
    "capillary_leak_syndrome": {
      "incidence": "25%",
      "severity": "Grade 3 in 5-8%",
      "mechanism": "Complement activation causes endothelial damage and increased vascular permeability",
      "management": "Monitor daily weights, serum albumin, and hemodynamic status during infusion days. Hold dinutuximab if systolic BP < 70 mmHg (age-adjusted) or albumin < 2.0 g/dL. IV albumin replacement PRN.",
      "monitoring": "Daily weights, albumin, vital signs q4h during infusion cycles"
    },
    "hypersensitivity_anaphylaxis": {
      "incidence": "15%",
      "severity": "Anaphylaxis in 1-2%",
      "mechanism": "IgE-mediated and complement-mediated hypersensitivity",
      "management": "Premedicate with diphenhydramine 1 mg/kg IV, acetaminophen 15 mg/kg PO, and hydrocortisone 1 mg/kg IV 30 minutes before each infusion. Epinephrine and resuscitation equipment at bedside.",
      "monitoring": "Continuous pulse oximetry, vital signs q15min during first hour of each infusion"
    },
    "hemolytic_uremic_syndrome": {
      "incidence": "Rare (<1%)",
      "severity": "Potentially life-threatening",
      "mechanism": "Complement-mediated thrombotic microangiopathy",
      "management": "Monitor CBC with schistocyte review, LDH, haptoglobin, creatinine, and urinalysis before each cycle. Discontinue dinutuximab if HUS develops.",
      "monitoring": "Pre-cycle labs: CBC, reticulocyte count, LDH, haptoglobin, BUN/creatinine, urinalysis"
    }
  },
  "recommended_monitoring_protocol": {
    "pre_infusion": "CBC, CMP, LDH, complement levels (C3, C4, CH50), urinalysis",
    "during_infusion": "Continuous SpO2, vital signs q15min x1hr then q1h, pain assessment q1h (FLACC scale)",
    "daily_during_cycle": "Weight, intake/output, serum albumin, vital signs q4h",
    "post_cycle": "CBC, CMP, LDH, complement levels — repeat at 1 week post-infusion"
  },
  "cross_agent_flags": [
    "Cardiology Agent: Capillary leak syndrome requires careful fluid management — coordinate with cardiac monitoring protocol given anthracycline cardiotoxicity risk (CBR3 Val/Val + RARG)",
    "Neurology Agent: Dinutuximab neuropathic pain is distinct from vincristine neuropathy — additive risk during concurrent administration. Coordinate pain management protocols.",
    "Oncology Agent: If grade 4 irAE occurs, consider reduced dinutuximab dosing (14 mg/m2/day) or switch to naxitamab (humanized anti-GD2 with potentially lower complement activation)"
  ],
  "processing_time_sec": 2.6
}
```

**Clinical interpretation:** The Autoimmune Agent profiles dinutuximab's immune-mediated toxicity spectrum, which is particularly relevant for Marcus. The 85% incidence of neuropathic pain — caused by anti-GD2 binding to peripheral nerve GD2 gangliosides — compounds the vincristine neuropathy risk already flagged by the Neurology Agent. The capillary leak syndrome risk (25%) creates a critical intersection with the Cardiology Agent's monitoring protocol: Marcus's anthracycline-stressed heart must tolerate fluid shifts from capillary leak. The cross-agent flags highlight these intersections, ensuring no monitoring gap between the cardiac, neurologic, and immunologic surveillance plans.

---

### Step 6 -- Oncology Agent: Integrated Dose-Adjusted Protocol

Synthesize all agent findings into a unified treatment protocol.

```bash
curl -X POST http://localhost:8527/v1/onco/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "marcus_germline",
    "query": "Generate integrated pharmacogenomically-adjusted treatment protocol for 6yo male with high-risk MYCN-amplified neuroblastoma Stage 4. Incorporate: TPMT *1/*3A (50% 6-MP reduction), CBR3 Val/Val + RARG rs2229774 (mandatory dexrazoxane), CYP3A5 *1/*3 (vincristine neuropathy monitoring).",
    "clinical_context": {
      "diagnosis": "High-risk neuroblastoma, MYCN amplified, Stage 4, unfavorable histology",
      "age_years": 6,
      "weight_kg": 21,
      "bsa_m2": 0.82,
      "protocol": "COG ANBL1232",
      "pgx_findings": {
        "TPMT": "*1/*3A -- intermediate metabolizer, 50% 6-MP reduction",
        "CBR3": "Val/Val -- high cardiotoxicity risk",
        "RARG": "rs2229774 C/T -- high cardiotoxicity risk",
        "CYP3A5": "*1/*3 -- intermediate metabolizer, vincristine monitoring"
      }
    },
    "output_format": "protocol_summary"
  }'
```

**Expected response:**

```json
{
  "query_id": "ONCO-2026-0323-003",
  "patient_id": "marcus_germline",
  "integrated_protocol": {
    "protocol_name": "COG ANBL1232 -- Pharmacogenomically Adjusted",
    "induction": {
      "cycles_1_3_5": {
        "cyclophosphamide": "400 mg/m2/day x5 days",
        "topotecan": "1.2 mg/m2/day x5 days",
        "modifications": "No PGx-based adjustments for this block"
      },
      "cycles_2_4": {
        "cisplatin": "50 mg/m2/day x4 days",
        "etoposide": "200 mg/m2/day x3 days",
        "daunorubicin": "60 mg/m2 x1 day",
        "dexrazoxane": "600 mg/m2 IV 30 min before daunorubicin (MANDATORY -- CBR3 Val/Val + RARG)",
        "modifications": "Dexrazoxane added per COG ALTE11C2 data and genetic high-risk profile"
      },
      "all_cycles": {
        "vincristine": "1.5 mg/m2 (max 2 mg) per protocol schedule",
        "vincristine_monitoring": "Ped-mTNS before each dose. Reduce to 1.0 mg/m2 if grade 2 neuropathy.",
        "cardiac_monitoring": "Troponin I + BNP each cycle. Echo every 2 cycles."
      }
    },
    "consolidation": {
      "high_dose_chemotherapy": "Busulfan/melphalan with autologous stem cell rescue",
      "cardiac_clearance": "LVEF >= 50% required"
    },
    "maintenance": {
      "isotretinoin": "160 mg/m2/day x14 days per 28-day cycle x6 cycles",
      "dinutuximab": "17.5 mg/m2/day x4 days per cycle x5 cycles",
      "6_mercaptopurine": {
        "standard_dose": "50 mg/m2/day",
        "adjusted_dose": "25 mg/m2/day (50% reduction for TPMT *1/*3A)",
        "monitoring": "Weekly CBC x8 weeks then biweekly. Target ANC 750-1500/uL.",
        "rationale": "TPMT intermediate metabolizers accumulate thioguanine nucleotides at 3-5x normal rate (Relling et al., Clin Pharmacol Ther 2019)"
      }
    }
  },
  "summary": "Protocol adjusted at 3 points: (1) Dexrazoxane added to all anthracycline cycles for CBR3/RARG risk, (2) 6-MP reduced 50% for TPMT *1/*3A, (3) Enhanced vincristine neuropathy monitoring for CYP3A5 *1/*3.",
  "processing_time_sec": 3.8
}
```

---

### Step 7 -- Therapeutic Discovery Engine: Novel ALK Inhibitors with Cardiac Safety

Marcus's genomic profiling revealed ALK F1174L as the primary oncogenic driver in his MYCN-amplified neuroblastoma. This mutation confers resistance to crizotinib but remains sensitive to next-generation ALK inhibitors. Given Marcus's elevated cardiotoxicity risk (CBR3 Val/Val + RARG rs2229774), the Therapeutic Discovery Engine is tasked with generating ALK inhibitors specifically optimized for cardiac safety — the central theme of this demo.

```bash
curl -X POST http://localhost:8510/api/v1/discover \
  -H "Content-Type: application/json" \
  -d '{
    "target_gene": "ALK",
    "target_mutation": "F1174L",
    "target_domain": "kinase_domain",
    "target_structure_pdb": "2XP2",
    "patient_id": "marcus_neuro_001",
    "optimization_constraints": {
      "pediatric_safe": true,
      "herg_liability": "strict_exclusion",
      "cardiotoxicity_filter": "stringent",
      "hepatotoxicity_filter": "elevated_sensitivity",
      "max_molecular_weight": 500,
      "age_group": "pediatric_3_8",
      "exclude_known_cardiotoxins": true,
      "cardiac_safety_priority": "highest"
    },
    "clinical_context": "MYCN-amplified neuroblastoma with ALK F1174L. Patient has CBR3 Val/Val and RARG rs2229774 conferring 35-45% baseline anthracycline cardiomyopathy risk. Cardiac safety is paramount.",
    "num_candidates": 100,
    "top_k": 3
  }'
```

**Expected response:**

```json
{
  "discovery_run_id": "DISC-2026-0323-ALK-001",
  "target": "ALK kinase domain (F1174L mutant)",
  "pdb_structure": "2XP2",
  "status": "completed",
  "candidates_generated": 100,
  "candidates_passing_filters": 11,
  "top_candidates": [
    {
      "rank": 1,
      "candidate_id": "HCLS-ALK-001",
      "name": "Lorlatinib analog (MolMIM cardiac-optimized)",
      "smiles": "CC(C1=C(C=CC(=C1)NC2=NC=C3C(=N2)N(C=C3)C4CCC(CC4)N5CCOCC5)F)O",
      "target": "ALK F1174L kinase domain",
      "rationale": "MolMIM-optimized from lorlatinib scaffold with modifications to eliminate hERG channel binding. Lorlatinib is the most active ALK inhibitor against F1174L, but standard formulation carries QTc prolongation risk. This analog maintains ALK potency while reducing cardiac ion channel interactions.",
      "scores": {
        "diffdock_binding_affinity_kcal": -9.1,
        "herg_ic50_predicted_uM": 48.2,
        "herg_liability": "none",
        "cardiac_safety_index": 0.93,
        "hepatotoxicity_risk": "low",
        "molecular_weight": 463.5,
        "pediatric_safety_score": 0.86,
        "druglikeness_qed": 0.74
      },
      "pediatric_considerations": "Designed for a patient with pre-existing cardiac risk factors. hERG IC50 >40 uM (safe margin). No QTc prolongation predicted. Compatible with concurrent dexrazoxane cardioprotection. Liver-friendly metabolism via non-CYP3A4 pathway to protect immature hepatic function."
    },
    {
      "rank": 2,
      "candidate_id": "HCLS-ALK-002",
      "name": "Macrocyclic ALK inhibitor (MolMIM novel scaffold)",
      "smiles": "C1CC(=O)N(C1)C2=NC3=CC=C(C=C3N2)NC(=O)C4=CC=C(C=C4)C#N",
      "target": "ALK F1174L kinase domain",
      "rationale": "Novel macrocyclic scaffold designed for high selectivity against ALK F1174L over wild-type ALK and other kinases. Macrocyclic constraint reduces off-target kinase inhibition that drives cardiotoxicity in first-generation ALK inhibitors.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.6,
        "herg_ic50_predicted_uM": 62.7,
        "herg_liability": "none",
        "cardiac_safety_index": 0.96,
        "hepatotoxicity_risk": "low",
        "molecular_weight": 371.4,
        "pediatric_safety_score": 0.89,
        "druglikeness_qed": 0.82
      },
      "pediatric_considerations": "Highest cardiac safety index (0.96) of all candidates. Low MW (371 Da) favorable for pediatric dosing. Macrocyclic structure limits off-target effects. Predicted oral bioavailability >60% enabling once-daily dosing for treatment adherence in a 6-year-old."
    },
    {
      "rank": 3,
      "candidate_id": "HCLS-ALK-003",
      "name": "Covalent ALK inhibitor (MolMIM selective warhead)",
      "smiles": "C=CC(=O)NC1=CC=C(C=C1)C2=CN=C(N=C2N)NC3=CC(=CC=C3)OC",
      "target": "ALK F1174L Cys1156 covalent site",
      "rationale": "Covalent binding to Cys1156 near the F1174L mutation site provides sustained target engagement, allowing lower dosing to reduce systemic exposure and cardiac load. Acrylamide warhead has favorable safety profile compared to chloroacetamide-based covalent inhibitors.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.9,
        "herg_ic50_predicted_uM": 55.3,
        "herg_liability": "none",
        "cardiac_safety_index": 0.91,
        "hepatotoxicity_risk": "moderate",
        "molecular_weight": 387.4,
        "pediatric_safety_score": 0.82,
        "druglikeness_qed": 0.77
      },
      "pediatric_considerations": "Covalent mechanism enables lower doses, reducing total cardiac exposure. However, moderate hepatotoxicity risk requires monitoring — particularly important given Marcus's immature hepatic metabolism (age 6). Recommended only if Candidates #1 and #2 show insufficient efficacy."
    }
  ],
  "pipeline_metadata": {
    "molmim_candidates_generated": 100,
    "diffdock_docked_against": "2XP2",
    "rdkit_admet_scored": 100,
    "pediatric_filter_passed": 22,
    "cardiac_safety_filter_passed": 11,
    "herg_filter_eliminated": 44,
    "cardiotoxicity_filter_eliminated": 11,
    "total_processing_time_seconds": 1156
  }
}
```

**Clinical interpretation:** The cardiac safety constraint is the defining feature of this discovery run. Of 100 MolMIM-generated candidates, 44 were eliminated by the hERG screen alone, and an additional 11 by the broader cardiotoxicity filter — reflecting the reality that kinase inhibitors frequently carry cardiac liability. For Marcus, whose CBR3 Val/Val and RARG rs2229774 genotypes already confer 35-45% anthracycline cardiomyopathy risk, any additional cardiac insult could be catastrophic. Candidate #2 achieves the highest cardiac safety index (0.96) with its macrocyclic scaffold that minimizes off-target kinase inhibition. The pipeline demonstrates how pharmacogenomic risk factors (identified by the PGx and Cardiology agents in earlier steps) directly constrain the drug design space, creating a feedback loop from diagnosis to discovery.

**Presenter note:** This step ties together the entire cardiac safety narrative of Demo 3. The same CBR3/RARG variants that triggered dexrazoxane addition and enhanced cardiac monitoring now constrain the drug discovery pipeline. Point out the 55% cardiac-related elimination rate (55 of 100 candidates) — this is dramatically higher than the ~30% rate seen in non-cardiac-risk patients (Demo 1), illustrating how patient-specific PGx profiles reshape the therapeutic search space.

---

### Clinical Insights -- Demo 3

1. **Pharmacogenomics transforms cardiac risk prediction.** Standard monitoring uses cumulative dose thresholds alone. By incorporating CBR3 and RARG genotypes, the system identifies Marcus as having 35-45% baseline cardiomyopathy risk, necessitating dexrazoxane from the first dose rather than waiting until 300 mg/m2.

2. **Biomarker integration confirms risk stratification.** The Biomarker Agent aggregates LDH, ferritin, catecholamines, cytogenetics, and ploidy into a unified prognostic panel. Every marker aligns with high-risk classification, driving maximum treatment intensity and informing downstream agent assessments.

3. **TPMT genotyping prevents maintenance toxicity.** Without PGx testing, Marcus would receive standard-dose 6-mercaptopurine and likely develop severe myelosuppression. The 50% dose reduction maintains efficacy while avoiding toxicity.

4. **Immunotherapy toxicity profiling prevents compounding harm.** The Autoimmune Agent identifies dinutuximab's 85% neuropathic pain incidence and 25% capillary leak risk — both of which intersect critically with Marcus's existing vincristine neuropathy risk and anthracycline cardiac vulnerability. Cross-agent flags ensure coordinated monitoring across neurologic, cardiac, and immunologic surveillance.

5. **Multi-agent coordination creates a safety net.** PGx identifies risks, Biomarker confirms risk stratification, Cardiology creates the monitoring plan, Neurology sets dose-modification criteria, Autoimmune profiles immunotherapy toxicity, and Oncology synthesizes everything into one actionable protocol. Seven agents contribute to a single unified care plan.

6. **Genetic risk compounding is non-obvious.** CBR3 Val/Val (5.5x) plus RARG rs2229774 (4.7x) likely create multiplicative rather than additive risk. Multi-agent analysis surfaces these compounding interactions that single-gene analysis misses.

---
---

## Section 8: DEMO 4 -- Rare Disease with Cancer Predisposition

### Clinical Narrative

**Patient:** Aurora T., 4-year-old female
**Presenting complaint:** Strabismus (esotropia) of the right eye detected at well-child visit. Fundoscopic exam reveals bilateral white pupillary reflex (leukocoria). Parents report no family history of eye cancer.

**Ophthalmic workup:**
- Right eye: 8 mm x 7 mm x 6 mm endophytic retinal mass in the macula, vitreous seeding (group D per International Classification of Retinoblastoma)
- Left eye: Two smaller foci (3 mm and 2 mm) in the peripheral retina, no vitreous seeds (group B)
- **Bilateral retinoblastoma** -- highly suspicious for hereditary (germline RB1 mutation)

**Systemic workup:**
- Brain MRI: No pinealoblastoma (trilateral retinoblastoma screening negative at diagnosis)
- Metastatic workup: CSF cytology negative, bone marrow negative, bone scan negative
- CBC, LFTs, renal function: All within normal limits

**Molecular workup ordered:** Whole-genome sequencing of blood (germline) to identify the causative RB1 mutation. FASTQ files generated on Illumina NovaSeq 6000, 150bp paired-end, ~30x germline coverage.

**Agents involved:** Rare Disease Agent (:8134), Oncology Agent (:8527), Imaging Agent (:8524), Clinical Trial Agent (:8538)

**Clinical questions:**
1. What is the germline RB1 mutation and its ACMG classification?
2. What is the management plan for bilateral retinoblastoma?
3. What cancer surveillance is needed for hereditary RB1 carriers?
4. Are Aurora's siblings at risk and should they undergo cascade genetic testing?

---

### Architecture Overview

```
+-------------------------------+       +-------------------------------+
|   GENOMIC FOUNDATION ENGINE   |       |   PRECISION INTELLIGENCE      |
|                               |       |   NETWORK                     |
|  FASTQ --> BWA-MEM2 --> BAM   |       |                               |
|  BAM --> DeepVariant --> VCF   | ----> |  Rare Disease Agent (:8134)   |
|  VCF --> ClinVar/AlphaMiss    |       |  Oncology Agent (:8527)       |
|  Focus: RB1 (13q14.2)        |       |  Imaging Agent (:8524)        |
+-------------------------------+       |  Clinical Trial Agent (:8538) |
                                        |  → Therapeutic Discovery      |
                                        |    Engine (:8510)             |
                                        +-------------------------------+
```

---

### Entry Point A: Starting from FASTQ

#### Step A1 -- Upload FASTQ and Run Germline Pipeline

```bash
mkdir -p /data/genomics/inputs/lily_germline

cp lily_germline_R1.fastq.gz /data/genomics/inputs/lily_germline/
cp lily_germline_R2.fastq.gz /data/genomics/inputs/lily_germline/

curl -X POST http://localhost:5000/v1/genomics/run \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "lily_germline",
    "pipeline_type": "germline_only",
    "input_fastq_r1": "/data/genomics/inputs/lily_germline/lily_germline_R1.fastq.gz",
    "input_fastq_r2": "/data/genomics/inputs/lily_germline/lily_germline_R2.fastq.gz",
    "reference_genome": "GRCh38",
    "target_regions": "cancer_predisposition_genes",
    "caller": "deepvariant",
    "output_dir": "/data/genomics/outputs/lily_germline"
  }'
```

**Expected response:**

```json
{
  "run_id": "GEN-2026-0323-004",
  "sample_id": "lily_germline",
  "status": "running",
  "pipeline_type": "germline_only",
  "estimated_completion_minutes": 45,
  "output_vcf": "/data/genomics/outputs/lily_germline/lily_germline.deepvariant.g.vcf.gz"
}
```

---

### Entry Point B: Starting from VCF

```bash
curl -X POST http://localhost:5000/v1/genomics/ingest-vcf \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "lily_germline",
    "vcf_path": "/data/genomics/external/lily_germline.vcf.gz",
    "reference_genome": "GRCh38",
    "annotate": true,
    "embed": true,
    "annotation_sources": ["clinvar", "alphamissense", "gnomad", "omim"]
  }'
```

---

### Step 1 -- Rare Disease Agent: RB1 Variant Classification

Query the Rare Disease Agent with HPO terms for ACMG-guided classification of the RB1 variant.

```bash
curl -X POST http://localhost:8134/v1/rare/diagnose \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "lily_germline",
    "query": "Classify the identified RB1 germline variant using ACMG/AMP criteria. Patient is a 4-year-old female with bilateral retinoblastoma.",
    "hpo_terms": [
      "HP:0009919",
      "HP:0030084",
      "HP:0000486",
      "HP:0007700",
      "HP:0000555"
    ],
    "hpo_labels": [
      "Retinoblastoma",
      "Bilateral retinoblastoma",
      "Strabismus",
      "Leukocoria",
      "Leukocoria"
    ],
    "gene_focus": ["RB1"],
    "locus": "13q14.2",
    "variant_of_interest": {
      "gene": "RB1",
      "hgvs_c": "c.958C>T",
      "hgvs_p": "p.Arg320Ter",
      "consequence": "stop_gained",
      "exon": "10/27",
      "zygosity": "heterozygous",
      "gnomad_af": 0.0,
      "clinvar_id": "VCV000012345",
      "clinvar_significance": "Pathogenic"
    },
    "acmg_classification_requested": true,
    "include_evidence_codes": true
  }'
```

**Expected response:**

```json
{
  "query_id": "RARE-2026-0323-004",
  "patient_id": "lily_germline",
  "diagnosis": {
    "disease": "Hereditary Retinoblastoma (OMIM #180200)",
    "gene": "RB1",
    "variant": "c.958C>T (p.Arg320Ter)",
    "inheritance": "Autosomal dominant with high penetrance (~90%)"
  },
  "acmg_classification": {
    "final_classification": "PATHOGENIC",
    "criteria_met": [
      {
        "code": "PVS1",
        "strength": "Very Strong",
        "description": "Null variant (nonsense) in RB1, a gene where loss of function is an established mechanism of disease. Truncating variants in RB1 are definitively pathogenic for retinoblastoma."
      },
      {
        "code": "PM2",
        "strength": "Moderate",
        "description": "Absent from gnomAD (0/251,450 alleles). Not observed in any population database."
      },
      {
        "code": "PP4",
        "strength": "Supporting",
        "description": "Patient phenotype (bilateral retinoblastoma at age 4) is highly specific for germline RB1 mutations. Bilateral disease has >95% positive predictive value for hereditary retinoblastoma."
      },
      {
        "code": "PP5",
        "strength": "Supporting",
        "description": "ClinVar classified as Pathogenic by multiple submitters (3-star review status)."
      }
    ],
    "classification_logic": "PVS1 + PM2 + PP4 + PP5 = Pathogenic (exceeds threshold: 1 Very Strong + 1 Moderate = Pathogenic)",
    "confidence": "HIGH"
  },
  "penetrance_data": {
    "retinoblastoma_risk": "~90% by age 5 (bilateral in ~60% of germline carriers)",
    "second_primary_cancer_risk": {
      "osteosarcoma": "7% cumulative incidence by age 20 (Kleinerman et al., J Clin Oncol 2012)",
      "soft_tissue_sarcoma": "3.5% by age 50",
      "melanoma": "2-4% lifetime risk",
      "lung_cancer": "increased risk, especially with radiation exposure",
      "bladder_cancer": "increased risk after cyclophosphamide exposure"
    }
  },
  "processing_time_sec": 3.2
}
```

---

### Step 2 -- Rare Disease Agent: Integrated Assessment

Cross-agent coordination for comprehensive care planning.

```bash
curl -X POST http://localhost:8134/v1/rare/integrated-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "lily_germline",
    "assessment_type": "hereditary_cancer_predisposition",
    "clinical_context": {
      "diagnosis": "Bilateral retinoblastoma, hereditary RB1 (c.958C>T, p.Arg320Ter)",
      "age_years": 4,
      "sex": "female",
      "right_eye_classification": "Group D (large tumor, vitreous seeding)",
      "left_eye_classification": "Group B (small peripheral foci)",
      "trilateral_screen": "Negative at diagnosis",
      "metastatic_workup": "Negative"
    },
    "agents_requested": ["oncology", "imaging", "clinical_trial"],
    "include_cascade_testing": true,
    "include_surveillance_plan": true
  }'
```

**Expected response:**

```json
{
  "assessment_id": "RARE-INTEG-2026-0323-004",
  "patient_id": "lily_germline",
  "integrated_findings": {
    "diagnosis_confirmed": "Hereditary retinoblastoma, RB1 p.Arg320Ter (Pathogenic)",
    "management_summary": "Globe-sparing therapy for both eyes is the primary goal. Right eye (group D) requires intensive intra-arterial chemotherapy +/- intravitreal melphalan. Left eye (group B) amenable to focal therapy (laser/cryotherapy) after systemic chemotherapy.",
    "cross_agent_flags": [
      "ONCOLOGY: Second primary cancer surveillance is critical and lifelong",
      "IMAGING: Trilateral retinoblastoma screening with brain MRI every 6 months until age 5",
      "CLINICAL TRIAL: Active COG ARET0321 and Rb-focused trials available",
      "GENETIC COUNSELING: Cascade testing mandatory for siblings and future pregnancies"
    ]
  },
  "processing_time_sec": 4.8
}
```

---

### Step 3 -- Oncology Agent: Retinoblastoma Management

```bash
curl -X POST http://localhost:8527/v1/onco/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "lily_germline",
    "query": "Management plan for bilateral hereditary retinoblastoma in a 4-year-old: right eye group D (vitreous seeding), left eye group B. Include treatment sequencing, second primary cancer surveillance protocol, and radiation avoidance rationale.",
    "clinical_context": {
      "diagnosis": "Bilateral hereditary retinoblastoma, RB1 p.Arg320Ter",
      "age_years": 4,
      "sex": "female",
      "right_eye": "Group D, 8mm macular mass, vitreous seeding",
      "left_eye": "Group B, two peripheral foci (3mm, 2mm)",
      "metastatic_status": "Negative"
    }
  }'
```

**Expected response:**

```json
{
  "query_id": "ONCO-2026-0323-004",
  "patient_id": "lily_germline",
  "treatment_plan": {
    "primary_goal": "Globe preservation bilateral. Avoid enucleation if possible. AVOID EXTERNAL BEAM RADIATION -- germline RB1 carriers have markedly increased second cancer risk in irradiated fields.",
    "right_eye_group_D": {
      "first_line": "Intra-arterial chemotherapy (IAC) with melphalan 5 mg via ophthalmic artery, 3 cycles at 4-week intervals",
      "adjuvant": "Intravitreal melphalan (20-30 mcg) for vitreous seeds, typically 6-8 injections at 1-2 week intervals",
      "monitoring": "Exam under anesthesia (EUA) with RetCam imaging every 3-4 weeks during active treatment",
      "salvage": "If progressive disease despite IAC + intravitreal: enucleation with long optic nerve section (>10 mm)"
    },
    "left_eye_group_B": {
      "first_line": "Systemic chemotherapy (vincristine 0.05 mg/kg, etoposide 5 mg/kg, carboplatin 18.6 mg/kg -- VEC protocol) x 6 cycles at 3-week intervals",
      "consolidation": "Focal therapy: Transpupillary thermotherapy (TTT) or cryotherapy to residual tumors after chemoreduction",
      "monitoring": "EUA with RetCam monthly during treatment, then every 3 months"
    },
    "radiation_avoidance": "External beam radiation therapy (EBRT) is CONTRAINDICATED as first-line in hereditary retinoblastoma. Germline RB1 carriers exposed to radiation have 38% cumulative incidence of second cancers by age 50 vs 6% without radiation (Kleinerman et al., J Clin Oncol 2005). Proton beam may be considered as a last resort if globe-sparing chemotherapy fails.",
    "second_primary_cancer_surveillance": {
      "schedule": [
        {"age_range": "0-5 years", "assessment": "Brain MRI every 6 months (trilateral screening). Annual physical exam."},
        {"age_range": "5-10 years", "assessment": "Annual brain MRI. Annual whole-body MRI starting age 8."},
        {"age_range": "10-20 years", "assessment": "Annual whole-body MRI. Annual dermatologic exam. Bone pain workup promptly."},
        {"age_range": "20+ years", "assessment": "Annual whole-body MRI. Annual dermatologic exam. Consider breast MRI for female carriers starting age 30. Annual CBC."}
      ],
      "radiation_exposure_minimization": "Avoid CT scans whenever possible. Use MRI and ultrasound preferentially. Dental radiographs with lead shielding."
    }
  },
  "processing_time_sec": 3.6
}
```

---

### Step 4 -- Imaging Agent: Trilateral Screening and Staging

```bash
curl -X POST http://localhost:8524/v1/imaging/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "lily_germline",
    "query": "Imaging surveillance protocol for hereditary retinoblastoma with germline RB1 mutation. Include trilateral screening (brain MRI for pinealoblastoma), orbital imaging, and systemic staging protocol.",
    "clinical_context": {
      "diagnosis": "Bilateral hereditary retinoblastoma, RB1 p.Arg320Ter",
      "age_years": 4,
      "current_trilateral_status": "Negative at diagnosis",
      "radiation_contraindication": true
    }
  }'
```

**Expected response:**

```json
{
  "query_id": "IMG-2026-0323-004",
  "patient_id": "lily_germline",
  "imaging_protocol": {
    "trilateral_screening": {
      "rationale": "Germline RB1 carriers have 5-15% risk of primitive neuroectodermal tumor (PNET) of the pineal or suprasellar region (trilateral retinoblastoma). Peak incidence is 20-36 months but can occur up to age 5.",
      "protocol": "Brain MRI with gadolinium contrast every 6 months until age 5",
      "focus_regions": ["Pineal gland", "Suprasellar region", "Posterior fossa"],
      "sequences": "T1 pre/post gadolinium, T2 FLAIR, DWI. Thin-slice (1 mm) through pineal region.",
      "current_status": "Negative at diagnosis (age 4). Continue q6mo for 1 more year."
    },
    "orbital_imaging": {
      "at_diagnosis": "Orbital MRI with fat suppression to assess extraocular extension",
      "during_treatment": "Orbital MRI pre-IAC to assess tumor response (every 2 cycles)",
      "post_treatment": "Orbital MRI every 6 months for 2 years, then annually for 5 years"
    },
    "staging_protocol": {
      "at_diagnosis": [
        "Brain + orbits MRI with gadolinium",
        "Lumbar puncture with CSF cytology (if extraocular extension suspected)",
        "Bone marrow aspirate/biopsy (if extraocular extension suspected)",
        "Bone scan or whole-body MRI (if metastatic disease suspected)"
      ],
      "current_status": "M0 -- no extraocular, no CNS, no metastatic disease"
    },
    "radiation_avoidance_note": "ALL imaging should be MRI-based. No CT scans unless emergent. Germline RB1 carriers have heightened radiation-induced second cancer risk (radiosensitivity due to haploinsufficiency of the G1/S checkpoint)."
  },
  "processing_time_sec": 2.6
}
```

---

### Step 5 -- Clinical Trial Agent: Trial Matching

```bash
curl -X POST http://localhost:8538/v1/trials/match \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "lily_germline",
    "diagnosis": "Bilateral hereditary retinoblastoma, group D (RE) / group B (LE)",
    "age_years": 4,
    "biomarkers": {
      "RB1": "c.958C>T, p.Arg320Ter (germline, pathogenic)"
    },
    "prior_therapy": "None (treatment-naive)",
    "search_criteria": {
      "disease_terms": ["retinoblastoma", "hereditary retinoblastoma", "RB1"],
      "therapy_types": ["intra-arterial chemotherapy", "intravitreal chemotherapy", "globe-sparing", "immunotherapy"],
      "pediatric_only": true,
      "phase": ["2", "3"],
      "status": "recruiting"
    }
  }'
```

**Expected response:**

```json
{
  "query_id": "TRIAL-2026-0323-004",
  "patient_id": "lily_germline",
  "matched_trials": [
    {
      "nct_id": "NCT03042429",
      "title": "COG ARET0321: A Study of Vincristine, Carboplatin, and Etoposide with Local Ophthalmic Therapy for Intraocular Retinoblastoma",
      "phase": "Phase 3",
      "status": "Active, not recruiting (follow-up ongoing)",
      "relevance": "Standard-of-care reference protocol for systemic chemoreduction. Aurora's left eye treatment follows this backbone.",
      "eligibility_match": "PARTIAL -- enrollment closed but treatment approach is standard-of-care"
    },
    {
      "nct_id": "NCT05765045",
      "title": "Intra-Arterial Melphalan Combined with Intravitreal Melphalan for Advanced Intraocular Retinoblastoma",
      "phase": "Phase 2",
      "status": "Recruiting",
      "relevance": "Directly applicable to Aurora's right eye (group D with vitreous seeding). Evaluates combined IAC + IVT melphalan for globe salvage.",
      "eligibility_match": "HIGH -- age 4, group D, treatment-naive"
    },
    {
      "nct_id": "NCT04587544",
      "title": "Immune Checkpoint Inhibition for Recurrent/Refractory Retinoblastoma",
      "phase": "Phase 1/2",
      "status": "Recruiting",
      "relevance": "Salvage option if primary globe-sparing therapy fails. Novel PD-L1 inhibitor approach for retinoblastoma.",
      "eligibility_match": "NOT YET -- reserve for refractory disease"
    }
  ],
  "genetic_counseling_trials": [
    {
      "nct_id": "NCT04218201",
      "title": "Long-Term Outcomes and Second Cancer Risk in Hereditary Retinoblastoma Survivors",
      "phase": "Observational",
      "status": "Enrolling",
      "relevance": "Registry study for long-term follow-up. All hereditary Rb patients eligible."
    }
  ],
  "processing_time_sec": 2.9
}
```

---

### Step 6 -- Genetic Counseling: Cascade Testing

```bash
curl -X POST http://localhost:8134/v1/rare/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "lily_germline",
    "query": "Genetic counseling plan for cascade testing in the family of a child with hereditary retinoblastoma (RB1 c.958C>T, p.Arg320Ter). Patient has a 2-year-old sibling and parents are considering future pregnancies.",
    "clinical_context": {
      "proband": "Aurora, age 4, bilateral retinoblastoma",
      "variant": "RB1 c.958C>T, p.Arg320Ter, heterozygous",
      "family_history": "No known family history of retinoblastoma or cancer",
      "siblings": [
        {"relation": "brother", "age_years": 2, "status": "asymptomatic, not yet tested"}
      ],
      "parental_status": "Not yet tested"
    }
  }'
```

**Expected response:**

```json
{
  "query_id": "RARE-COUNSEL-2026-0323-004",
  "patient_id": "lily_germline",
  "cascade_testing_plan": {
    "urgency": "HIGH -- sibling at 50% risk if a parent carries the variant; even if de novo in Aurora, somatic mosaicism in a parent cannot be excluded",
    "parental_testing": {
      "recommendation": "Both parents should undergo targeted RB1 testing for c.958C>T immediately",
      "scenarios": [
        {
          "result": "One parent positive",
          "implication": "Sibling has 50% risk. That parent needs retinal exam and second-cancer surveillance. All future pregnancies have 50% risk.",
          "sibling_action": "Immediate targeted testing. If positive: EUA within 2 weeks, brain MRI for trilateral screening."
        },
        {
          "result": "Both parents negative (apparent de novo)",
          "implication": "Likely de novo in Aurora, but parental gonadal mosaicism cannot be excluded (estimated 3-5% risk). Sibling risk is lower (~1-3%) but NOT zero.",
          "sibling_action": "Still recommend targeted testing. If negative: clinical screening with fundoscopy every 6 months until age 5 (given residual mosaicism risk)."
        }
      ]
    },
    "sibling_testing": {
      "brother_age_2": {
        "urgency": "IMMEDIATE -- age 2 is within peak retinoblastoma risk window",
        "test": "Targeted RB1 c.958C>T sequencing on blood sample",
        "turnaround": "5-7 business days at most genetic testing laboratories",
        "if_positive": "EUA within 1-2 weeks. Fundoscopy under anesthesia. Brain MRI. Begin screening protocol: EUA every 3-4 weeks until age 3, then every 3 months to age 5.",
        "if_negative": "Reassurance, but maintain clinical vigilance with annual eye exams through age 5 (given mosaicism possibility)"
      }
    },
    "future_pregnancies": {
      "preconception": "Offer preimplantation genetic testing for monogenic disorders (PGT-M) for RB1 c.958C>T if either parent carries the variant",
      "prenatal": "Chorionic villus sampling (CVS) at 10-12 weeks or amniocentesis at 15-18 weeks for targeted RB1 testing",
      "postnatal_if_positive": "First EUA at 1-2 weeks of life. Brain MRI at 3 months. Intensive screening through age 5."
    },
    "psychosocial": "Refer family to retinoblastoma support network. Address guilt/anxiety in parents (especially if one is a carrier). Genetic counseling for Aurora at age-appropriate intervals regarding her own future reproductive risks."
  },
  "processing_time_sec": 3.0
}
```

---

### Step 7 -- Therapeutic Discovery Engine: Novel CDK4/6 Inhibitors for Ocular Delivery

Aurora's hereditary retinoblastoma is driven by biallelic RB1 loss, which removes the critical cell cycle checkpoint at the G1/S transition. With RB1 protein absent, CDK4/6 activity is unchecked, driving uncontrolled proliferation. CDK4/6 inhibition represents a rational therapeutic strategy to pharmacologically restore the checkpoint that RB1 loss eliminates. Critically, Aurora's bilateral disease and hereditary cancer predisposition make radiation contraindicated (38% second cancer risk), creating an urgent need for novel local therapies. The Therapeutic Discovery Engine generates CDK4/6 inhibitors optimized for intravitreal or periocular delivery to treat ocular disease while minimizing systemic exposure in a 4-year-old.

```bash
curl -X POST http://localhost:8510/api/v1/discover \
  -H "Content-Type: application/json" \
  -d '{
    "target_gene": "CDK4/CDK6",
    "target_domain": "ATP_binding_pocket",
    "target_structure_pdb": "1JOW",
    "patient_id": "lily_retino_001",
    "optimization_constraints": {
      "pediatric_safe": true,
      "delivery_route": "intravitreal_periocular",
      "max_molecular_weight": 400,
      "ocular_penetration": true,
      "minimize_systemic_exposure": true,
      "teratogenicity": "strict_exclusion",
      "age_group": "pediatric_2_6",
      "exclude_known_cardiotoxins": true
    },
    "clinical_context": "Bilateral hereditary retinoblastoma (RB1 p.Arg320Ter). Group D right eye, Group B left eye. Radiation contraindicated due to second cancer risk. Need locally delivered CDK4/6 inhibitors for globe preservation.",
    "num_candidates": 100,
    "top_k": 3
  }'
```

**Expected response:**

```json
{
  "discovery_run_id": "DISC-2026-0323-CDK46-001",
  "target": "CDK4/CDK6 ATP-binding pocket",
  "pdb_structure": "1JOW",
  "status": "completed",
  "candidates_generated": 100,
  "candidates_passing_filters": 9,
  "top_candidates": [
    {
      "rank": 1,
      "candidate_id": "HCLS-CDK46-OC-001",
      "name": "Ocular-optimized CDK4/6 inhibitor (MolMIM-generated, ribociclib scaffold)",
      "smiles": "CC(=O)NC1=NC=C2CN=C(NC3=CC(=CC=C3)N4CCNCC4)N=C2N1",
      "target": "CDK4/CDK6",
      "rationale": "MolMIM-optimized from ribociclib scaffold for enhanced vitreous humor solubility and retinal penetration. Low MW (348 Da) enables diffusion through vitreous to reach retinal tumors. Designed for intravitreal injection formulation.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.9,
        "ocular_penetration_score": 0.87,
        "vitreous_solubility_mg_ml": 2.4,
        "systemic_exposure_index": 0.08,
        "molecular_weight": 348.4,
        "teratogenicity_risk": "none_detected",
        "pediatric_safety_score": 0.90,
        "druglikeness_qed": 0.83
      },
      "pediatric_considerations": "Designed for local delivery to minimize systemic exposure in a 4-year-old (systemic exposure index 0.08). MW 348 Da enables vitreous diffusion. No teratogenicity signal — important for a child with lifetime cancer surveillance ahead. Compatible with concurrent systemic chemotherapy (vincristine-carboplatin-etoposide)."
    },
    {
      "rank": 2,
      "candidate_id": "HCLS-CDK46-OC-002",
      "name": "Periocular depot CDK4/6 inhibitor (MolMIM slow-release design)",
      "smiles": "O=C1N=C(NC2=CC=C(N3CCCC3)C=C2)C2=CN=CC=C2N1C1CC1",
      "target": "CDK4/CDK6",
      "rationale": "Designed for periocular (subconjunctival) depot injection with slow scleral penetration. Higher lipophilicity enables sustained release from periocular depot while maintaining CDK4/6 selectivity. Targets the smaller Group B lesions in Aurora's left eye where periocular delivery may suffice.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.5,
        "ocular_penetration_score": 0.72,
        "scleral_permeability": 0.68,
        "systemic_exposure_index": 0.05,
        "molecular_weight": 361.4,
        "teratogenicity_risk": "none_detected",
        "pediatric_safety_score": 0.88,
        "druglikeness_qed": 0.79
      },
      "pediatric_considerations": "Periocular route avoids repeated intravitreal injections — reducing procedural sedation burden in a 4-year-old. Lowest systemic exposure (0.05) of all candidates. Suitable for the less advanced left eye (Group B) where periocular penetration may be adequate."
    },
    {
      "rank": 3,
      "candidate_id": "HCLS-CDK46-OC-003",
      "name": "Dual CDK4/6-HDAC inhibitor for retinoblastoma (MolMIM novel scaffold)",
      "smiles": "OC(=O)CCCCCNC(=O)C1=CC=C(C=C1)C2=NC3=CC=CC=C3N2",
      "target": "CDK4/CDK6 + HDAC",
      "rationale": "Dual-target compound combining CDK4/6 inhibition with HDAC inhibition. Rationale: RB1-null retinoblastoma cells show aberrant histone acetylation patterns; HDAC inhibition synergizes with CDK4/6 blockade to induce differentiation and apoptosis. Hydroxamic acid moiety enables HDAC binding.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.1,
        "ocular_penetration_score": 0.79,
        "systemic_exposure_index": 0.12,
        "molecular_weight": 395.5,
        "teratogenicity_risk": "none_detected",
        "pediatric_safety_score": 0.82,
        "druglikeness_qed": 0.71
      },
      "pediatric_considerations": "Dual mechanism may enable lower individual target doses, reducing toxicity. However, HDAC inhibition carries thrombocytopenia risk — requires platelet monitoring. Reserve for refractory lesions not responding to single-target CDK4/6 inhibition. Higher systemic exposure (0.12) than Candidates #1-2 due to hydroxamic acid absorption."
    }
  ],
  "pipeline_metadata": {
    "molmim_candidates_generated": 100,
    "diffdock_docked_against": "1JOW",
    "rdkit_admet_scored": 100,
    "pediatric_filter_passed": 24,
    "ocular_delivery_filter_passed": 15,
    "teratogenicity_filter_eliminated": 6,
    "systemic_exposure_filter_eliminated": 9,
    "total_processing_time_seconds": 978
  }
}
```

**Clinical interpretation:** The Therapeutic Discovery Engine addresses the core therapeutic challenge in hereditary retinoblastoma: how to deliver effective therapy locally while avoiding the systemic and long-term risks that are amplified in a child with a germline cancer predisposition. All three candidates are designed for local ocular delivery with minimal systemic exposure (index 0.05-0.12, compared to 0.6-0.8 for oral CDK4/6 inhibitors). Candidate #1 is optimized for intravitreal delivery to the advanced right eye (Group D), while Candidate #2 targets the left eye's smaller lesions via the less invasive periocular route. The teratogenicity filter is particularly relevant — while Aurora is only 4, she carries a lifetime cancer risk, and any therapeutic agent she receives should carry no teratogenic legacy. The 91% elimination rate reflects the compound difficulty of simultaneously satisfying ocular penetration, minimal systemic exposure, and pediatric safety constraints.

**Presenter note:** This step demonstrates how the Therapeutic Discovery Engine adapts to fundamentally different therapeutic modalities. Unlike Demos 1-3 which seek oral systemic agents, Aurora's case demands local delivery optimization — a different pharmacological challenge entirely. Highlight that Candidate #2's periocular depot design reduces procedural sedation burden, a real clinical consideration when treating a 4-year-old who may need monthly treatments for years.

---

### Clinical Insights -- Demo 4

1. **ACMG classification provides definitive diagnosis.** The RB1 p.Arg320Ter variant meets Pathogenic criteria through PVS1 (null variant) + PM2 (absent from population databases) + PP4 (specific phenotype). This is not a variant of uncertain significance -- it is definitively actionable.

2. **Radiation avoidance is life-saving in hereditary Rb.** The system correctly flags that external beam radiation increases second cancer risk from 6% to 38% by age 50. This critical constraint shapes the entire treatment plan toward globe-sparing chemotherapy.

3. **Cascade testing has immediate clinical urgency.** Aurora's 2-year-old brother is within the peak retinoblastoma risk window. If he carries the variant, early detection through screening can catch tumors when they are small (group A/B) and curable with focal therapy alone, rather than waiting for symptoms when tumors may be advanced.

4. **Multi-agent integration addresses the full disease spectrum.** Hereditary retinoblastoma is not just an eye cancer -- it is a lifelong cancer predisposition syndrome. The Rare Disease Agent identifies the variant, Oncology plans treatment, Imaging creates the surveillance protocol, and Clinical Trial matching identifies active studies. No single specialist typically covers all four domains.

5. **The platform models real genetic counseling workflows.** By generating cascade testing plans with scenario-based recommendations, the system supports the actual decision tree that genetic counselors use, including the nuance of gonadal mosaicism risk when both parents test negative.

---
---

## Section 9: DEMO 5 -- CAR-T Therapy Decision and Monitoring

### 10.1 Clinical Narrative

**Patient:** Ethan M., 12-year-old male
**Diagnosis:** Relapsed/refractory B-cell acute lymphoblastic leukemia (B-ALL)
**Weight:** 42 kg | **BSA:** 1.32 m² | **Karnofsky/Lansky:** 80%

**Clinical History:**
Ethan was diagnosed with National Cancer Institute (NCI) high-risk B-ALL at age 10. Initial
induction chemotherapy per COG AALL1231 achieved morphologic remission, but day 29
minimal residual disease (MRD) by flow cytometry was 5% — a devastating prognostic
indicator associated with 5-year event-free survival (EFS) below 30%. Reinduction with
FLAG-IDA (fludarabine, cytarabine, G-CSF, idarubicin) failed to achieve MRD negativity.
Bone marrow at day 42 of reinduction showed 12% lymphoblasts with preserved CD19
expression (mean fluorescence intensity 8,500).

His parents, Sarah and Michael, are desperate for options. The pediatric oncology team
presents two paths: chimeric antigen receptor T-cell (CAR-T) therapy with tisagenlecleucel
(Kymriah), or proceeding directly to allogeneic hematopoietic stem cell transplantation
(HSCT) from his HLA-matched sibling. The family wants a data-driven understanding of
risks, timelines, and outcomes for both approaches.

**Key Clinical Parameters:**
- Cytogenetics: Hyperdiploidy (trisomy 4, 10, 17), no Philadelphia chromosome
- Immunophenotype: CD19+ (>95%), CD22+ (85%), CD10+, TdT+
- Prior anthracycline exposure: Daunorubicin cumulative dose 250 mg/m²
- Cardiac function: LVEF 58% by echocardiography
- Hepatic function: ALT 28 U/L, AST 32 U/L, total bilirubin 0.6 mg/dL
- Renal function: Creatinine 0.5 mg/dL, eGFR >120 mL/min/1.73m²
- CNS status: CNS-1 (no blasts in CSF)
- Infections: None active; CMV IgG positive, EBV IgG positive

### 10.2 Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────────────────┐
│                    CAR-T THERAPY DECISION PIPELINE                          │
│                                                                              │
│  Entry Point A (FASTQ)              Entry Point B (VCF)                     │
│  Tumor DNA + Germline DNA           Pre-analyzed molecular data             │
│         │                                    │                               │
│         ▼                                    │                               │
│  ┌─────────────────────┐                     │                               │
│  │ Genomic Foundation  │                     │                               │
│  │ Engine (Parabricks/ │                     │                               │
│  │ DeepVariant)        │                     │                               │
│  └────────┬────────────┘                     │                               │
│           │ VCF output                       │                               │
│           ▼                                  ▼                               │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 1: Single-Cell Agent /v1/sc/query                      │            │
│  │   CD19 expression validation, TME analysis, MRD tracking    │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 2: CAR-T Agent /v1/cart/integrated-assessment          │            │
│  │   Eligibility, CRS/ICANS risk, manufacturing timeline       │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 3: Cardiology Agent /v1/cardio/query                   │            │
│  │   Pre-lymphodepletion cardiac assessment                    │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 4: Autoimmune Agent /v1/autoimmune/integrated-assessment│           │
│  │   irAE risk profiling, cytokine storm monitoring plan       │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 5: PGx Agent /v1/pgx/query                             │            │
│  │   Tocilizumab metabolism, corticosteroid pharmacogenomics   │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 6: Clinical Trial Agent /v1/trial/match                │            │
│  │   CAR-T trials, novel constructs, combination studies       │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 7: Therapeutic Discovery Engine /api/v1/discover        │            │
│  │   CD22 bridging therapy candidates (MolMIM→DiffDock→RDKit)  │            │
│  └──────────────────────────────────────────────────────────────┘            │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 10.3 Entry Point A: FASTQ Processing

Tumor and germline FASTQ files are submitted to the Genomic Foundation Engine:

```bash
# Submit paired tumor + germline FASTQ for somatic analysis
curl -X POST http://localhost:5000/v1/genomics/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "ETHAN_CART_001",
    "tumor_fastq_r1": "/data/fastq/ethan_tumor_R1.fastq.gz",
    "tumor_fastq_r2": "/data/fastq/ethan_tumor_R2.fastq.gz",
    "germline_fastq_r1": "/data/fastq/ethan_germline_R1.fastq.gz",
    "germline_fastq_r2": "/data/fastq/ethan_germline_R2.fastq.gz",
    "reference": "GRCh38",
    "analysis_type": "somatic_paired",
    "pipeline": "parabricks_deepvariant"
  }'
```

**Expected output:** VCF with somatic variants, germline variants for PGx, and BAM files
for downstream analysis. Typical runtime on DGX Spark: 2-4 hours for paired tumor/normal
at 60x/30x coverage.

### 10.4 Entry Point B: VCF Direct Input

For sites with existing molecular profiling:

```bash
curl -X POST http://localhost:5000/v1/genomics/annotate \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "ETHAN_CART_001",
    "vcf_path": "/data/vcf/ethan_somatic.vcf.gz",
    "germline_vcf_path": "/data/vcf/ethan_germline.vcf.gz",
    "annotations": ["clinvar", "cosmic", "alphamissense", "gnomad"]
  }'
```

### 10.5 Step 1: Single-Cell Analysis — CD19 Expression Validation

The single-cell agent processes flow cytometry and scRNA-seq data to validate CAR-T
target expression and characterize the tumor microenvironment.

```bash
curl -X POST http://localhost:8540/v1/sc/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ETHAN_CART_001",
    "query": "Validate CD19 expression on B-ALL blasts, characterize tumor microenvironment, and assess MRD clone architecture",
    "data_sources": {
      "flow_cytometry": "/data/flow/ethan_mrd_panel.fcs",
      "scrna_seq": "/data/scrna/ethan_bm_10x.h5"
    },
    "analysis_type": "cart_target_validation"
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "analysis": {
    "cd19_expression": {
      "blast_percentage_cd19_positive": 97.2,
      "mean_fluorescence_intensity": 8500,
      "cd19_expression_heterogeneity": "low",
      "dim_cd19_population": "2.8% of blasts with MFI < 1000",
      "clinical_interpretation": "High uniform CD19 expression supports anti-CD19 CAR-T. Minimal dim-CD19 population reduces antigen escape risk."
    },
    "cd22_expression": {
      "blast_percentage_cd22_positive": 85.3,
      "mean_fluorescence_intensity": 4200,
      "clinical_interpretation": "Adequate CD22 expression provides viable backup target if CD19 loss occurs post-CAR-T."
    },
    "tumor_microenvironment": {
      "total_cells_analyzed": 12847,
      "blast_percentage": 78.4,
      "cd8_t_cells_percentage": 2.1,
      "cd4_t_cells_percentage": 3.8,
      "nk_cells_percentage": 1.2,
      "monocytes_percentage": 8.3,
      "tme_classification": "immune_desert",
      "exhaustion_markers": {
        "pd1_on_cd8": "68% positive",
        "tim3_on_cd8": "42% positive",
        "lag3_on_cd8": "31% positive"
      },
      "clinical_interpretation": "Immune desert phenotype with severely depleted and exhausted endogenous T-cells. This does NOT preclude CAR-T efficacy as manufactured CAR-T cells bypass endogenous immune suppression."
    },
    "mrd_clone_tracking": {
      "dominant_clone_vdj": "IGHV3-23/IGHD2-2/IGHJ4",
      "clone_frequency": 0.92,
      "subclonal_populations": 3,
      "subclone_with_cd19_loss": "none detected",
      "clinical_interpretation": "Monoclonal disease with stable CD19 expression across all subclones. Low risk of pre-existing CD19-negative escape."
    }
  }
}
```

**Clinical Significance:** The >95% CD19 positivity with high MFI confirms Ethan is an
excellent candidate for anti-CD19 CAR-T. The immune desert TME, while prognostically
unfavorable for checkpoint inhibitors, does not impair CAR-T efficacy. The absence of
pre-existing CD19-negative subclones is reassuring, though post-infusion antigen loss
remains a 10-20% risk per published ELIANA data.

### 10.6 Step 2: CAR-T Agent — Integrated Assessment

```bash
curl -X POST http://localhost:8522/v1/cart/integrated-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ETHAN_CART_001",
    "indication": "relapsed_refractory_b_all",
    "age_years": 12,
    "weight_kg": 42,
    "bsa_m2": 1.32,
    "prior_therapies": ["AALL1231_induction", "FLAG_IDA"],
    "blast_percentage": 12,
    "cd19_status": {"positive_pct": 97.2, "mfi": 8500},
    "organ_function": {
      "lvef_pct": 58,
      "alt_ul": 28,
      "ast_ul": 32,
      "bilirubin_mgdl": 0.6,
      "creatinine_mgdl": 0.5,
      "egfr_ml_min": 120
    },
    "cns_status": "CNS1",
    "hla_matched_sibling": true,
    "assess": ["eligibility", "crs_risk", "icans_risk", "manufacturing", "alternatives"]
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "eligibility_assessment": {
    "tisagenlecleucel": {
      "eligible": true,
      "indication_met": "Relapsed/refractory B-ALL in patient ≤25 years",
      "contraindications_present": [],
      "organ_function_adequate": true,
      "rationale": "FDA-approved for r/r B-ALL. Age 12, adequate organ function, no active infections, CNS-1 status. All ELIANA-equivalent eligibility criteria met."
    }
  },
  "risk_assessment": {
    "crs": {
      "overall_risk": "moderate",
      "grade_any_probability": 0.77,
      "grade_3_4_probability": 0.47,
      "risk_factors": [
        "High disease burden (12% blasts) — strongest CRS predictor",
        "Age 12 — moderate pediatric risk",
        "Prior FLAG-IDA — enhanced inflammatory priming"
      ],
      "mitigation": "Tocilizumab (12 mg/kg IV) available bedside. PICU bed reserved. IL-6 and ferritin monitoring q6h for first 14 days.",
      "evidence": "ELIANA trial: 77% any-grade CRS, 47% grade ≥3. Median onset day 3 (range 1-22). Disease burden >5% blasts is primary risk factor (Maude et al., NEJM 2018)."
    },
    "icans": {
      "overall_risk": "moderate",
      "any_grade_probability": 0.40,
      "grade_3_4_probability": 0.13,
      "risk_factors": [
        "High disease burden",
        "Anticipated CRS severity (correlated with ICANS)"
      ],
      "monitoring": "ICE score q8h days 1-28. Neurology consultation for grade ≥2.",
      "evidence": "ELIANA: 40% any-grade neurotoxicity, 13% grade ≥3. Typically follows CRS by 1-5 days."
    }
  },
  "manufacturing_plan": {
    "product": "Tisagenlecleucel (Kymriah)",
    "manufacturer": "Novartis/Cellular Therapy",
    "estimated_manufacturing_time_days": 22,
    "leukapheresis_window": "Within 7 days of enrollment",
    "bridging_therapy": {
      "recommended": true,
      "rationale": "22-day manufacturing window with 12% blasts requires disease control",
      "options": [
        {"regimen": "Dexamethasone 10 mg/m²/day x 14 days", "rationale": "Avoid lymphotoxic agents that impair T-cell collection"},
        {"regimen": "Vincristine 1.5 mg/m² IV weekly x 3", "rationale": "Non-lymphotoxic debulking"},
        {"regimen": "Low-dose cytarabine 100 mg/m²/day x 5", "rationale": "If rapid debulking needed"}
      ],
      "avoid": "Cyclophosphamide or clofarabine — may impair collected T-cell viability"
    },
    "lymphodepletion": {
      "regimen": "Fludarabine 30 mg/m²/day x 4 days + Cyclophosphamide 500 mg/m²/day x 2 days",
      "timing": "Complete 2-6 days before CAR-T infusion",
      "rationale": "Creates lymphopenic niche for CAR-T expansion"
    }
  },
  "cd22_backup_strategy": {
    "rationale": "10-20% of patients develop CD19-negative relapse post-CAR-T",
    "cd22_expression_confirmed": true,
    "options": [
      "Inotuzumab ozogamicin (anti-CD22 ADC, FDA-approved for r/r ALL)",
      "Anti-CD22 CAR-T (clinical trials: NCT02315612)",
      "Dual CD19/CD22 CAR-T (clinical trials: NCT03448393)"
    ]
  },
  "comparison_with_hsct": {
    "cart_advantages": [
      "No graft-vs-host disease risk",
      "No prolonged immunosuppression",
      "ELIANA: 82% CR rate, 66% 12-month EFS",
      "Outpatient monitoring possible after initial period"
    ],
    "hsct_advantages": [
      "HLA-matched sibling donor available — best HSCT scenario",
      "Proven long-term curative potential (55-65% EFS with matched sibling)",
      "No risk of B-cell aplasia (engraftment restores B-cells)",
      "Established long-term follow-up data (decades)"
    ],
    "recommendation": "CAR-T first, HSCT as consolidation or salvage. ELIANA data supports CAR-T as bridge-to-transplant or definitive therapy depending on response depth and MRD status at day 28 post-infusion."
  }
}
```

**Clinical Significance:** The integrated assessment confirms tisagenlecleucel eligibility
and provides a comprehensive risk-benefit analysis. The 77% CRS rate requires PICU
availability and tocilizumab on hand. The 22-day manufacturing window necessitates
bridging therapy to control Ethan's 12% blasts — dexamethasone is preferred as it
avoids lymphotoxicity that could compromise the leukapheresis product.

### 10.7 Step 3: Cardiology Agent — Pre-Lymphodepletion Cardiac Assessment

```bash
curl -X POST http://localhost:8126/v1/cardio/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ETHAN_CART_001",
    "query": "Pre-lymphodepletion cardiac risk assessment for CAR-T therapy",
    "clinical_context": {
      "age_years": 12,
      "prior_anthracycline_cumulative_mg_m2": 250,
      "current_lvef_pct": 58,
      "troponin_ng_ml": 0.01,
      "bnp_pg_ml": 45,
      "ecg": "normal_sinus_rhythm_80bpm",
      "planned_regimen": "fludarabine_cyclophosphamide_lymphodepletion"
    }
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "cardiac_assessment": {
    "baseline_function": {
      "lvef": "58% — normal (lower limit 55%)",
      "gls": "-19.2% — normal (>-18% abnormal)",
      "troponin": "0.01 ng/mL — normal",
      "bnp": "45 pg/mL — normal (<100)",
      "interpretation": "Preserved systolic function. No evidence of subclinical anthracycline cardiotoxicity despite 250 mg/m² cumulative daunorubicin."
    },
    "anthracycline_risk": {
      "cumulative_dose": "250 mg/m² daunorubicin (doxorubicin equivalent ~210 mg/m²)",
      "risk_category": "moderate",
      "lifetime_cardiomyopathy_risk": "5-10% at current cumulative dose",
      "recommendation": "Avoid further anthracyclines. If HSCT is pursued, non-anthracycline conditioning preferred."
    },
    "lymphodepletion_cardiac_risk": {
      "fludarabine": "Minimal direct cardiac toxicity",
      "cyclophosphamide_500_mg_m2": "Low risk at this dose (cardiotoxicity typically >150 mg/kg)",
      "overall_risk": "LOW",
      "recommendation": "Proceed with Flu/Cy lymphodepletion. No cardiac contraindication."
    },
    "crs_cardiac_monitoring": {
      "risk": "CRS-associated cardiac dysfunction in 5-15% of severe CRS",
      "monitoring_plan": [
        "Troponin and BNP at baseline, day 3, 7, 14, 28",
        "Echocardiography if grade ≥3 CRS",
        "Continuous telemetry during CRS (typically days 1-14)",
        "Urgent echo if troponin rise >3x baseline or new arrhythmia"
      ],
      "interventions": {
        "grade_1_2_crs": "Supportive care, continue monitoring",
        "grade_3_crs_with_cardiac": "Tocilizumab + consider low-dose vasopressors",
        "grade_4_crs_with_cardiac": "Tocilizumab + corticosteroids + PICU + inotropes PRN"
      }
    }
  }
}
```

### 10.8 Step 4: Autoimmune Agent — Immune-Related Adverse Event Profiling

```bash
curl -X POST http://localhost:8532/v1/autoimmune/integrated-assessment \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ETHAN_CART_001",
    "context": "post_cart_immune_toxicity_risk",
    "therapy": "tisagenlecleucel",
    "assess": ["cytokine_storm", "autoimmune_cytopenias", "immune_myocarditis", "hematologic_toxicity"]
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "irae_risk_profile": {
    "cytokine_release_syndrome": {
      "mechanism": "CAR-T activation → massive IL-6, IFN-gamma, IL-10, TNF-alpha release",
      "key_cytokines_to_monitor": {
        "il6": {"expected_peak": "day 3-7", "threshold_pg_ml": 1000, "action": "tocilizumab if >1000 + clinical grade ≥2"},
        "ifn_gamma": {"expected_peak": "day 2-5", "threshold_pg_ml": 500},
        "ferritin": {"expected_peak": "day 5-10", "threshold_ng_ml": 10000, "action": "suggests macrophage activation syndrome overlap"},
        "crp": {"expected_peak": "day 3-7", "threshold_mg_l": 200}
      },
      "macrophage_activation_syndrome": {
        "risk": "5-10% of severe CRS cases",
        "distinguishing_features": "Ferritin >10,000, triglycerides elevated, fibrinogen decreased, hepatosplenomegaly",
        "treatment": "Anakinra (IL-1 receptor antagonist) if refractory to tocilizumab"
      }
    },
    "autoimmune_cytopenias": {
      "incidence": "30-40% post-CAR-T",
      "types": {
        "neutropenia": {"incidence_pct": 35, "median_duration_days": 21, "management": "G-CSF if ANC <500 after day 28"},
        "thrombocytopenia": {"incidence_pct": 25, "median_duration_days": 28, "management": "Platelet transfusion if <10,000 or bleeding"},
        "anemia": {"incidence_pct": 20, "median_duration_days": 14, "management": "Transfusion support"}
      },
      "mechanism": "Immune-mediated destruction + bone marrow suppression from lymphodepletion",
      "monitoring": "CBC with differential daily x 28 days, then weekly x 3 months"
    },
    "immune_myocarditis": {
      "incidence": "rare (<1%)",
      "risk_factors": ["Prior anthracycline exposure", "Severe CRS with high IL-6"],
      "presentation": "New arrhythmia, troponin elevation, reduced LVEF during CRS",
      "management": "High-dose corticosteroids, cardiology consultation, consider IVIG",
      "monitoring": "Serial troponin (already in CRS monitoring plan)"
    },
    "b_cell_aplasia": {
      "expected": true,
      "mechanism": "On-target off-tumor: CAR-T eliminates normal CD19+ B-cells",
      "duration": "Months to years (indicates CAR-T persistence)",
      "management": "IVIG replacement (400 mg/kg IV q4 weeks) when IgG <400 mg/dL",
      "infection_prophylaxis": {
        "bacterial": "Prophylactic antibiotics per institutional protocol",
        "viral": "Acyclovir 400 mg BID for HSV/VZV prophylaxis",
        "fungal": "Fluconazole during neutropenic period",
        "pjp": "TMP-SMX or pentamidine for Pneumocystis prophylaxis"
      }
    }
  }
}
```

### 10.9 Step 5: PGx Agent — Supportive Medication Pharmacogenomics

```bash
curl -X POST http://localhost:8107/v1/pgx/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ETHAN_CART_001",
    "germline_vcf": "/data/vcf/ethan_germline.vcf.gz",
    "medications": ["tocilizumab", "dexamethasone", "fludarabine", "cyclophosphamide", "acyclovir", "fluconazole"],
    "query": "Pharmacogenomic assessment for CAR-T supportive medications"
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "pgx_assessment": {
    "tocilizumab": {
      "drug_class": "IL-6 receptor monoclonal antibody",
      "metabolism": "Proteolytic degradation (not CYP-mediated)",
      "pgx_relevance": "None — monoclonal antibodies not metabolized by CYP enzymes",
      "dosing": "12 mg/kg IV (weight <30 kg) or 8 mg/kg IV (weight ≥30 kg). Ethan 42 kg → 336 mg IV. Maximum 2 doses in 24 hours.",
      "interactions": "IL-6 suppression normalizes CYP3A4 activity — may affect concurrent CYP3A4 substrates"
    },
    "dexamethasone": {
      "metabolism": "CYP3A4 (primary)",
      "pgx_variants_detected": {
        "CYP3A4": "*1/*1 (normal metabolizer)",
        "NR3C1": "No clinically relevant variants detected"
      },
      "dosing_adjustment": "None required. Standard dosing: 10 mg/m² BID for bridging.",
      "note": "During CRS, elevated IL-6 may inhibit CYP3A4 → transiently increased dexamethasone exposure. Clinical significance minimal at standard CRS doses."
    },
    "fludarabine": {
      "metabolism": "Intracellular phosphorylation to active F-ara-ATP. Renal excretion (40-60%).",
      "pgx_relevance": "Limited PGx data. DCK (deoxycytidine kinase) polymorphisms may affect activation.",
      "dosing": "30 mg/m²/day x 4 days for lymphodepletion. No PGx-based adjustment needed.",
      "renal_function": "eGFR >120 — no dose reduction required"
    },
    "cyclophosphamide": {
      "metabolism": "CYP2B6 (primary activation), CYP3A4, CYP2C9",
      "pgx_variants_detected": {
        "CYP2B6": "*1/*1 (normal metabolizer)",
        "CYP2C9": "*1/*1 (normal metabolizer)"
      },
      "dosing_adjustment": "None required. 500 mg/m²/day x 2 days for lymphodepletion.",
      "hemorrhagic_cystitis_risk": "Low at this dose. Mesna not required per most protocols."
    },
    "acyclovir": {
      "metabolism": "Renal excretion (60-90% unchanged)",
      "pgx_relevance": "Minimal. No CYP-mediated metabolism.",
      "dosing": "400 mg PO BID for HSV/VZV prophylaxis. Adequate renal function."
    },
    "fluconazole": {
      "metabolism": "CYP2C9 inhibitor, CYP3A4 inhibitor",
      "pgx_variants_detected": {
        "CYP2C9": "*1/*1 — no increased sensitivity"
      },
      "interaction_alert": "Fluconazole inhibits CYP3A4 — monitor dexamethasone levels if co-administered. Unlikely to be clinically significant at prophylactic fluconazole doses (200 mg daily)."
    }
  }
}
```

### 10.10 Step 6: Clinical Trial Agent — CAR-T Trial Matching

```bash
curl -X POST http://localhost:8538/v1/trial/match \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "ETHAN_CART_001",
    "diagnosis": "relapsed_refractory_b_all",
    "age_years": 12,
    "prior_therapies": ["standard_induction", "FLAG_IDA"],
    "biomarkers": {"cd19_positive": true, "cd22_positive": true, "ph_negative": true},
    "therapy_interest": ["car_t", "bispecific_antibody", "novel_immunotherapy"]
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "matched_trials": [
    {
      "trial_id": "NCT02435849",
      "name": "ELIANA — Tisagenlecleucel for Pediatric/Young Adult r/r B-ALL",
      "phase": "II (registration trial, now post-marketing)",
      "status": "FDA approved — commercial access available",
      "relevance_score": 0.98,
      "key_data": {
        "overall_remission_rate": "82%",
        "complete_remission_rate": "63%",
        "twelve_month_efs": "66%",
        "twelve_month_os": "76%"
      },
      "note": "Commercial product. Not a trial enrollment but provides benchmark data."
    },
    {
      "trial_id": "NCT02625480",
      "name": "ZUMA-4 — Axicabtagene Ciloleucel for Pediatric r/r B-ALL",
      "phase": "I/II",
      "status": "Recruiting",
      "relevance_score": 0.92,
      "eligibility_match": ["Age 2-21", "r/r B-ALL", "CD19+", "Adequate organ function"],
      "sites": ["NCI", "MD Anderson", "CHOP", "Seattle Children's"],
      "key_difference": "CD28 costimulatory domain (vs 4-1BB in tisagenlecleucel). Faster expansion, potentially higher CRS."
    },
    {
      "trial_id": "NCT03448393",
      "name": "Dual CD19/CD22 CAR-T for B-ALL",
      "phase": "I",
      "status": "Recruiting",
      "relevance_score": 0.88,
      "eligibility_match": ["r/r B-ALL", "CD19+", "CD22+", "Age ≥1"],
      "rationale": "Dual targeting reduces antigen escape risk. Particularly relevant given Ethan's CD22+ status.",
      "sites": ["Stanford", "NIH", "CHOP"]
    },
    {
      "trial_id": "NCT03876769",
      "name": "Blinatumomab + CAR-T Sequencing Study",
      "phase": "II",
      "status": "Recruiting",
      "relevance_score": 0.82,
      "rationale": "Blinatumomab (CD19/CD3 BiTE) as bridging before CAR-T. May reduce disease burden and improve CAR-T outcomes."
    },
    {
      "trial_id": "NCT04049383",
      "name": "CAR-T Followed by HSCT Consolidation — Pediatric ALL",
      "phase": "II",
      "status": "Recruiting",
      "relevance_score": 0.79,
      "rationale": "Relevant to Ethan's scenario: CAR-T first for remission, followed by matched sibling HSCT for consolidation."
    }
  ]
}
```

### 10.11 Step 7: Therapeutic Discovery Engine — Bridging Therapy Candidates for CAR-T Manufacturing Window

Ethan's CAR-T assessment reveals a 22-day manufacturing window during which his relapsed ALL must be controlled without compromising the T-cell product. The Single-Cell agent identified a CD19-dim escape clone (8.3% of blasts), making CD22 the backup target. The Therapeutic Discovery Engine now generates small-molecule CD22 modulators designed specifically for bridging therapy — compounds that control disease during manufacturing without suppressing the T-cells that will become his CAR-T product.

```bash
curl -X POST http://localhost:8510/api/v1/discover \
  -H "Content-Type: application/json" \
  -d '{
    "target_gene": "CD22",
    "target_domain": "extracellular_Ig_domains",
    "target_structure_pdb": "5VKJ",
    "patient_id": "ethan_cart_001",
    "optimization_constraints": {
      "pediatric_safe": true,
      "t_cell_sparing": true,
      "short_half_life_preferred": true,
      "max_half_life_hours": 24,
      "max_molecular_weight": 500,
      "age_group": "pediatric_10_16",
      "exclude_known_cardiotoxins": true,
      "immunosuppression": "strict_exclusion"
    },
    "clinical_context": "Relapsed/refractory B-ALL with CD19-dim escape clone. CD22 85% positive. 22-day CAR-T manufacturing window requires bridging therapy. MUST NOT suppress T-cells (would compromise CAR-T product). Short half-life preferred for washout before CAR-T infusion.",
    "num_candidates": 100,
    "top_k": 3
  }'
```

**Expected response:**

```json
{
  "discovery_run_id": "DISC-2026-0323-CD22-001",
  "target": "CD22 extracellular domain",
  "pdb_structure": "5VKJ",
  "status": "completed",
  "candidates_generated": 100,
  "candidates_passing_filters": 8,
  "top_candidates": [
    {
      "rank": 1,
      "candidate_id": "HCLS-CD22-001",
      "name": "CD22 small-molecule modulator (MolMIM B-cell selective)",
      "smiles": "CC1=CC(=CC(=C1)NC(=O)C2=CC=C(C=C2)S(=O)(=O)NC3=CC=CC=N3)OC(F)(F)F",
      "target": "CD22 ligand-binding domain",
      "rationale": "Targets CD22 sialic acid-binding domain to disrupt B-cell survival signaling without affecting T-cell populations. Designed for selective B-cell cytotoxicity via CD22-mediated internalization and apoptosis induction.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.4,
        "t_cell_impact_score": 0.02,
        "predicted_half_life_hours": 8.5,
        "molecular_weight": 465.4,
        "herg_liability": "none",
        "hepatotoxicity_risk": "low",
        "pediatric_safety_score": 0.85,
        "druglikeness_qed": 0.76
      },
      "pharmacokinetics": {
        "half_life_hours": 8.5,
        "washout_time_5_half_lives_hours": 42.5,
        "days_before_cart_to_stop": 2,
        "steady_state_days": 1.5
      },
      "pediatric_considerations": "T-cell sparing confirmed (impact score 0.02 — minimal). 8.5-hour half-life enables 2-day washout before CAR-T infusion, fitting within the manufacturing timeline. Rapid onset (steady state in 1.5 days) provides immediate disease control. No immunosuppressive properties that would compromise leukapheresis product."
    },
    {
      "rank": 2,
      "candidate_id": "HCLS-CD22-002",
      "name": "Sialic acid-CD22 interaction disruptor (MolMIM glycomimetic)",
      "smiles": "OC1C(O)C(OC(C1O)CO)NC(=O)C2=CC=C(C=C2)C3=CC=NO3",
      "target": "CD22 Siglec-2 domain",
      "rationale": "Glycomimetic compound that competes with natural sialic acid ligands for CD22 binding, disrupting the CD22-mediated B-cell receptor inhibitory checkpoint. By blocking CD22 ligand engagement, renders B-ALL blasts susceptible to immune-mediated killing while preserving T-cell function.",
      "scores": {
        "diffdock_binding_affinity_kcal": -7.9,
        "t_cell_impact_score": 0.01,
        "predicted_half_life_hours": 6.2,
        "molecular_weight": 382.3,
        "herg_liability": "none",
        "hepatotoxicity_risk": "minimal",
        "pediatric_safety_score": 0.88,
        "druglikeness_qed": 0.72
      },
      "pharmacokinetics": {
        "half_life_hours": 6.2,
        "washout_time_5_half_lives_hours": 31.0,
        "days_before_cart_to_stop": 1.5,
        "steady_state_days": 1.0
      },
      "pediatric_considerations": "Shortest half-life (6.2 hours) allows the fastest washout — therapy can continue until 1.5 days before CAR-T infusion, maximizing bridging coverage. Glycomimetic scaffold has favorable safety profile with minimal off-target effects. Lowest T-cell impact (0.01) of all candidates."
    },
    {
      "rank": 3,
      "candidate_id": "HCLS-CD22-003",
      "name": "CD22-targeted degrader (MolMIM PROTAC-inspired)",
      "smiles": "CC(NC(=O)COC1=CC=C(C=C1)CNC(=O)C2=CC=C3C(=O)NC(=O)C3=C2)C4=CC=CC=C4",
      "target": "CD22 surface protein (degradation)",
      "rationale": "PROTAC-inspired molecular glue that recruits cereblon E3 ligase to CD22, inducing proteasomal degradation of CD22 on blast surfaces. Surface CD22 loss may sensitize blasts to standard bridging agents (dexamethasone) while maintaining the CD22 intracellular epitope for potential CD22 CAR-T backup.",
      "scores": {
        "diffdock_binding_affinity_kcal": -8.1,
        "t_cell_impact_score": 0.05,
        "predicted_half_life_hours": 14.3,
        "molecular_weight": 487.5,
        "herg_liability": "none",
        "hepatotoxicity_risk": "moderate",
        "pediatric_safety_score": 0.79,
        "druglikeness_qed": 0.68
      },
      "pharmacokinetics": {
        "half_life_hours": 14.3,
        "washout_time_5_half_lives_hours": 71.5,
        "days_before_cart_to_stop": 3,
        "steady_state_days": 2.5
      },
      "pediatric_considerations": "Longest half-life (14.3 hours) requires 3-day washout — still feasible within 22-day manufacturing window. Moderate hepatotoxicity risk requires LFT monitoring. T-cell impact slightly higher (0.05) than Candidates #1-2 due to potential cereblon engagement in T-cells. Reserve for cases where Candidates #1-2 show insufficient blast control."
    }
  ],
  "pipeline_metadata": {
    "molmim_candidates_generated": 100,
    "diffdock_docked_against": "5VKJ",
    "rdkit_admet_scored": 100,
    "pediatric_filter_passed": 26,
    "t_cell_sparing_filter_passed": 12,
    "immunosuppression_filter_eliminated": 38,
    "half_life_filter_eliminated": 4,
    "total_processing_time_seconds": 1089
  }
}
```

**Clinical interpretation:** The defining constraint in Ethan's therapeutic discovery is the T-cell sparing requirement. Of 100 candidates, 38 were eliminated solely because they showed immunosuppressive properties that could compromise the leukapheresis product or the final CAR-T infusion. This is a unique pharmacological challenge: the bridging therapy must kill B-cell blasts while protecting T-cells — the very cell lineage being engineered into his treatment. Candidate #2 achieves the optimal balance with the lowest T-cell impact (0.01), shortest half-life (6.2 hours for rapid washout), and highest pediatric safety score (0.88). The pharmacokinetic data is particularly actionable: each candidate includes explicit washout timing to guide the clinical team on when to discontinue bridging therapy relative to CAR-T infusion day.

**Presenter note:** This step highlights a constraint unique to immuno-oncology: the therapy must selectively kill cancer cells while preserving the patient's own immune cells that are being repurposed as treatment. Point out the pharmacokinetic washout calculations — these directly inform clinical logistics around the 22-day manufacturing window. Also note that Candidate #3's PROTAC-inspired degrader represents cutting-edge targeted protein degradation applied to a pediatric context, showing the platform's awareness of emerging modalities.

---

### 10.11 Clinical Decision Summary for Ethan

The seven-step pipeline (six agents plus Therapeutic Discovery Engine) produces a unified clinical decision framework:

| Decision Point | Recommendation | Evidence Level |
|---|---|---|
| CAR-T vs HSCT first | CAR-T first (tisagenlecleucel) | ELIANA Phase II |
| Target validation | CD19 confirmed (97.2%, MFI 8,500) | Flow cytometry |
| CRS risk | Moderate (77% any grade, 47% severe) | ELIANA published data |
| Cardiac clearance | Approved for lymphodepletion | LVEF 58%, troponin normal |
| Bridging therapy | Dexamethasone 10 mg/m²/day x 14d | COG guidelines |
| PGx concerns | None for tocilizumab (non-CYP) | PharmGKB Level 1 |
| Backup target | CD22 (85% positive) available | Institutional data |
| Bridging candidates | 3 CD22 modulators (T-cell sparing, 6-14hr half-life) | Therapeutic Discovery Engine |
| Trial option | Dual CD19/CD22 (NCT03448393) | Phase I recruiting |

**Family communication:** Ethan's parents are told: "We recommend CAR-T therapy with
tisagenlecleucel. There is an 82% chance of complete remission. The main risks are
cytokine release syndrome (77% chance, usually manageable with tocilizumab) and
neurological symptoms (40% chance, usually reversible). Manufacturing takes about 3
weeks, during which Ethan will receive dexamethasone to control his leukemia. His heart
is healthy enough for the conditioning chemotherapy. If the leukemia comes back after
CAR-T, his matched sibling donor remains available for transplant, and we have dual
CD19/CD22 CAR-T trials as another option."

---

## Section 10: DEMO 6 -- Medulloblastoma Precision Treatment with Novel Drug Discovery

### 11.1 Clinical Narrative

**Patient:** Aiden K., 10-year-old male
**Diagnosis:** Newly diagnosed medulloblastoma, Sonic Hedgehog (SHH) subtype
**Weight:** 32 kg | **BSA:** 1.12 m² | **Lansky Performance:** 70%

**Clinical History:**
Aiden presented with 6 weeks of progressive morning headaches, vomiting, and unsteady
gait. Brain MRI revealed a 4.2 cm enhancing posterior fossa mass centered in the
cerebellar vermis with compression of the fourth ventricle and obstructive hydrocephalus.
Emergent external ventricular drain (EVD) placement was followed by near-total resection
(>95%). Post-operative MRI at 48 hours shows a 3 mm residual enhancing nodule along
the floor of the fourth ventricle.

Histopathology confirmed medulloblastoma with desmoplastic/nodular histology.
Molecular profiling revealed SHH subtype with a PTCH1 nonsense mutation
(c.1729C>T, p.Arg577*), TP53 wild-type, no MYCN amplification, and no chromosome 14q
loss. Spinal MRI and CSF cytology are negative for leptomeningeal dissemination (M0).

**Risk Stratification:**
- WHO Grade IV medulloblastoma
- SHH subtype with PTCH1 mutation — favorable within SHH group
- TP53 wild-type — critical favorable factor (TP53-mutant SHH has dismal prognosis)
- Desmoplastic/nodular histology — favorable
- M0 staging — no metastatic disease
- Near-total resection — <1.5 cm² residual
- **Overall risk: Average risk per current COG stratification**

Aiden's parents, Jennifer and David, understand standard treatment involves craniospinal
irradiation (CSI) at 23.4 Gy with posterior fossa boost to 54 Gy, followed by adjuvant
chemotherapy. They are deeply concerned about long-term neurocognitive effects,
endocrine dysfunction, and secondary malignancies. They want to explore whether targeted
SHH pathway inhibitors could reduce or replace radiation.

### 11.2 Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────────────────┐
│            MEDULLOBLASTOMA PRECISION TREATMENT PIPELINE                      │
│                                                                              │
│  Entry Point A (FASTQ)              Entry Point B (VCF)                     │
│  Tumor DNA + Germline DNA           Molecular profiling results             │
│         │                                    │                               │
│         ▼                                    │                               │
│  ┌─────────────────────┐                     │                               │
│  │ Genomic Foundation  │                     │                               │
│  │ Engine              │                     │                               │
│  └────────┬────────────┘                     │                               │
│           │ VCF output                       │                               │
│           ▼                                  ▼                               │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 2: Oncology Agent /v1/onco/query                       │            │
│  │   SHH-subtype interpretation, vismodegib contraindication   │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 3: Neurology Agent /v1/neuro/query                     │            │
│  │   Posterior fossa syndrome risk, radiation neurotoxicity     │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 4: Imaging Agent /v1/imaging/query                     │            │
│  │   Staging, post-op residual, surveillance planning          │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 5: Therapeutic Discovery Engine (MolMIM/DiffDock)      │            │
│  │   Novel SMO antagonists with pediatric BBB optimization     │            │
│  └────────────────────────┬─────────────────────────────────────┘            │
│                           ▼                                                  │
│  ┌──────────────────────────────────────────────────────────────┐            │
│  │ Step 6: Clinical Trial Agent /v1/trial/match                │            │
│  │   Pediatric brain tumor trials, SHH-targeted studies        │            │
│  └──────────────────────────────────────────────────────────────┘            │
└──────────────────────────────────────────────────────────────────────────────┘
```

### 11.3 Entry Point A: FASTQ Processing

```bash
curl -X POST http://localhost:5000/v1/genomics/analyze \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "AIDEN_MDB_001",
    "tumor_fastq_r1": "/data/fastq/aiden_tumor_R1.fastq.gz",
    "tumor_fastq_r2": "/data/fastq/aiden_tumor_R2.fastq.gz",
    "germline_fastq_r1": "/data/fastq/aiden_germline_R1.fastq.gz",
    "germline_fastq_r2": "/data/fastq/aiden_germline_R2.fastq.gz",
    "reference": "GRCh38",
    "analysis_type": "somatic_paired",
    "pipeline": "parabricks_deepvariant",
    "tumor_type": "medulloblastoma"
  }'
```

### 11.4 Entry Point B: VCF Direct Input

```bash
curl -X POST http://localhost:5000/v1/genomics/annotate \
  -H "Content-Type: application/json" \
  -d '{
    "sample_id": "AIDEN_MDB_001",
    "vcf_path": "/data/vcf/aiden_somatic.vcf.gz",
    "germline_vcf_path": "/data/vcf/aiden_germline.vcf.gz",
    "annotations": ["clinvar", "cosmic", "alphamissense", "gnomad", "cancerhotspots"]
  }'
```

**Key Variants Identified:**

| Gene | Variant | Type | Classification | Clinical Significance |
|---|---|---|---|---|
| PTCH1 | c.1729C>T, p.Arg577* | Nonsense | Pathogenic | Loss of SHH pathway suppressor |
| TP53 | Wild-type | — | — | Favorable prognostic factor |
| MYCN | No amplification | CNV | — | No high-risk amplification |
| SUFU | Wild-type | — | — | Rules out SUFU-driven SHH |
| SMO | Wild-type | — | — | Upstream activation (PTCH1 loss) |
| TERT | Promoter wild-type | — | — | No telomerase activation |
| KMT2D | c.8390G>A, p.Arg2797Gln | Missense | VUS | Chromatin modifier — significance uncertain |

### 11.5 Step 2: Oncology Agent — SHH-Subtype Medulloblastoma

```bash
curl -X POST http://localhost:8527/v1/onco/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "AIDEN_MDB_001",
    "query": "Treatment recommendations for SHH-subtype medulloblastoma with PTCH1 loss in a 10-year-old",
    "diagnosis": {
      "tumor_type": "medulloblastoma",
      "molecular_subtype": "SHH",
      "histology": "desmoplastic_nodular",
      "key_mutations": ["PTCH1_p.Arg577*"],
      "tp53_status": "wild_type",
      "mycn_status": "not_amplified",
      "m_stage": "M0",
      "extent_of_resection": "near_total"
    },
    "patient": {"age_years": 10, "sex": "male"}
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "molecular_interpretation": {
    "shh_pathway_activation": {
      "mechanism": "PTCH1 nonsense mutation (p.Arg577*) → loss of Patched1 tumor suppressor → constitutive release of SMO inhibition → unchecked GLI1/GLI2 transcription",
      "pathway_level": "Upstream (receptor level)",
      "therapeutic_implication": "SMO inhibitors (vismodegib, sonidegib) act DOWNSTREAM of PTCH1 and should be effective against PTCH1-loss tumors"
    },
    "prognostic_factors": {
      "favorable": [
        "TP53 wild-type (5-year OS >80% vs <40% if TP53-mutant)",
        "Desmoplastic/nodular histology",
        "No MYCN amplification",
        "M0 staging",
        "Age >5 years (infant SHH has different biology)",
        "Near-total resection"
      ],
      "unfavorable": [],
      "overall_prognosis": "5-year EFS 75-85% with standard therapy"
    }
  },
  "treatment_recommendations": {
    "standard_of_care": {
      "surgery": "Completed — near-total resection with 3mm residual",
      "radiation": {
        "protocol": "Craniospinal irradiation 23.4 Gy + posterior fossa boost to 54 Gy",
        "duration": "6 weeks",
        "proton_therapy_benefit": "STRONGLY RECOMMENDED — reduces exit dose to cochlea (hearing preservation), hippocampus (neurocognition), hypothalamus (endocrine), and reduces secondary malignancy risk by 50-60%"
      },
      "chemotherapy": {
        "protocol": "Modified Packer regimen per COG ACNS0331/ACNS0332",
        "during_radiation": "Vincristine 1.5 mg/m² IV weekly x 6-8 weeks",
        "maintenance": "Cisplatin 75 mg/m² + CCNU 75 mg/m² + Vincristine 1.5 mg/m², 6-week cycles x 8"
      }
    },
    "shh_targeted_therapy": {
      "vismodegib": {
        "mechanism": "SMO antagonist — blocks Smoothened receptor downstream of PTCH1",
        "efficacy_in_shh_mdb": "40-50% response rate in recurrent SHH medulloblastoma (Robinson et al., JCO 2015)",
        "CRITICAL_PEDIATRIC_CONCERN": "CONTRAINDICATED in growing children. Causes IRREVERSIBLE premature growth plate closure (physeal fusion). Published case reports of growth arrest in children treated with vismodegib. FDA label: not approved for pediatric use.",
        "age_threshold": "Generally avoided in patients with open growth plates (males typically until age 16-18)",
        "aiden_recommendation": "DO NOT USE — Aiden is 10 years old with open growth plates"
      },
      "sonidegib": {
        "mechanism": "SMO antagonist (different binding site than vismodegib)",
        "same_growth_plate_concern": true,
        "aiden_recommendation": "DO NOT USE — same class effect on growth plates"
      },
      "alternative_shh_targets": {
        "gli_inhibitors": "GANT61 — preclinical, targets GLI transcription factors downstream of SMO. Avoids growth plate concerns theoretically, but no pediatric clinical data.",
        "itraconazole": "Antifungal with SMO inhibitory activity. Weaker than vismodegib. Some pediatric safety data from antifungal use. Being explored in combination approaches.",
        "arsenic_trioxide": "GLI1 inhibitor. Pediatric safety established in APL. Being studied in SHH medulloblastoma (NCT04299113)."
      }
    },
    "clinical_bottom_line": "Standard therapy (CSI + chemotherapy) remains the backbone for average-risk SHH medulloblastoma in a 10-year-old. Proton therapy should be pursued to minimize late effects. FDA-approved SHH inhibitors (vismodegib, sonidegib) are CONTRAINDICATED due to growth plate toxicity. Novel SHH-targeting approaches (GLI inhibitors, drug discovery) are the appropriate avenue for targeted therapy in pediatric patients."
  }
}
```

**Clinical Significance:** The oncology agent identifies the critical pediatric safety
concern — vismodegib causes irreversible growth plate closure and cannot be used in a
10-year-old with open physes. This drives the pipeline toward the Therapeutic Discovery
Engine to identify novel SMO antagonists without growth plate toxicity.

### 11.6 Step 3: Neurology Agent — Posterior Fossa Syndrome and Late Effects

```bash
curl -X POST http://localhost:8528/v1/neuro/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "AIDEN_MDB_001",
    "query": "Posterior fossa syndrome risk assessment and neurocognitive late effects of craniospinal irradiation",
    "clinical_context": {
      "tumor_location": "cerebellar_vermis_4th_ventricle",
      "surgery": "near_total_resection",
      "planned_csi_gy": 23.4,
      "posterior_fossa_boost_gy": 54,
      "age_at_radiation": 10,
      "proton_therapy_planned": true
    }
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "posterior_fossa_syndrome": {
    "incidence": "25% of posterior fossa tumor resections",
    "risk_factors_present": [
      "Vermian tumor location (highest risk)",
      "Near-total resection (brainstem manipulation)",
      "Fourth ventricle involvement"
    ],
    "predicted_risk": "HIGH (30-35% given vermian location)",
    "clinical_features": {
      "cerebellar_mutism": "Inability to speak, onset 1-5 days post-operatively",
      "ataxia": "Truncal and appendicular, may persist months",
      "emotional_lability": "Common accompaniment",
      "cranial_nerve_palsies": "6th and 7th nerve, typically transient"
    },
    "current_status": "Aiden should be assessed — if surgery was recent, posterior fossa syndrome may still be evolving",
    "recovery_timeline": {
      "speech": "Typically recovers over 1-6 months, may have residual dysarthria",
      "ataxia": "Improves over 3-12 months, some patients have permanent mild imbalance",
      "overall": "80% have functional recovery, 20% have lasting deficits"
    },
    "management": "Speech therapy, occupational therapy, physical therapy. No pharmacological intervention proven effective."
  },
  "neurocognitive_late_effects": {
    "craniospinal_irradiation_impact": {
      "iq_decline": {
        "estimated_annual_loss": "1.5-3 IQ points per year post-radiation",
        "cumulative_by_5_years": "7-15 IQ point decline",
        "primary_deficit": "Processing speed and working memory, with relatively preserved verbal comprehension",
        "age_factor": "Greater decline when treated at younger ages (<7 most vulnerable). Age 10 is intermediate risk."
      },
      "white_matter_changes": {
        "incidence": "70-90% within 2 years on MRI",
        "clinical_correlation": "Correlates with processing speed deficits",
        "proton_benefit": "Reduced dose to hippocampi and supratentorial white matter with proton therapy"
      },
      "academic_impact": "40-60% require educational accommodations, 20% require special education"
    },
    "proton_vs_photon_comparison": {
      "hippocampal_dose_reduction": "60-70% lower with proton therapy",
      "cochlear_dose_reduction": "50-70% lower → reduced hearing loss",
      "hypothalamic_dose_reduction": "40-60% lower → reduced endocrine dysfunction",
      "estimated_iq_preservation": "3-5 IQ points preserved vs photon at 5 years (Yock et al., Lancet Oncol 2016)",
      "secondary_malignancy_reduction": "50-60% lower integral body dose → reduced secondary cancer risk"
    },
    "endocrine_late_effects": {
      "growth_hormone_deficiency": {"incidence_pct": 80, "onset": "1-3 years post-radiation", "management": "GH replacement therapy"},
      "hypothyroidism": {"incidence_pct": 40, "management": "Levothyroxine replacement"},
      "precocious_puberty": {"incidence_pct": 15, "management": "GnRH agonist if indicated"},
      "adrenal_insufficiency": {"incidence_pct": 5, "management": "Hydrocortisone replacement"},
      "gonadal_dysfunction": {"incidence_pct": 30, "management": "Fertility preservation counseling before treatment"}
    },
    "monitoring_plan": [
      "Neuropsychological testing at baseline, 1, 2, and 5 years post-treatment",
      "Annual endocrine panel (GH stimulation, TSH/fT4, cortisol, LH/FSH, IGF-1)",
      "Audiology every 6 months during cisplatin and annually thereafter",
      "Annual ophthalmology for posterior subcapsular cataracts",
      "Educational liaison — 504 plan or IEP as needed"
    ]
  }
}
```

### 11.7 Step 4: Imaging Agent — Staging and Surveillance

```bash
curl -X POST http://localhost:8524/v1/imaging/query \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "AIDEN_MDB_001",
    "query": "Complete staging assessment and post-treatment surveillance plan for medulloblastoma",
    "imaging_data": {
      "pre_op_brain_mri": "/data/imaging/aiden_preop_brain.nii.gz",
      "post_op_brain_mri": "/data/imaging/aiden_postop_brain.nii.gz",
      "spinal_mri": "/data/imaging/aiden_spine.nii.gz"
    }
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "staging_assessment": {
    "pre_operative_brain_mri": {
      "tumor_location": "Cerebellar vermis with fourth ventricle extension",
      "maximum_diameter_cm": 4.2,
      "enhancement_pattern": "Heterogeneous enhancement with cystic components — consistent with desmoplastic/nodular histology",
      "hydrocephalus": "Obstructive hydrocephalus, lateral ventricle dilation",
      "brainstem_invasion": "Tumor abuts but does not invade brainstem",
      "leptomeningeal_enhancement": "None detected on pre-operative scan"
    },
    "post_operative_brain_mri": {
      "timing": "48 hours post-resection (within recommended 24-72 hour window)",
      "residual_tumor": {
        "present": true,
        "maximum_diameter_cm": 0.3,
        "location": "Floor of fourth ventricle",
        "volume_cm3": 0.014,
        "classification": "Near-total resection (residual <1.5 cm²)"
      },
      "post_surgical_changes": "Expected post-operative blood products in surgical cavity. No diffusion restriction to suggest ischemia.",
      "hydrocephalus_status": "Improved with EVD in situ"
    },
    "spinal_mri": {
      "leptomeningeal_metastases": "NONE DETECTED",
      "drop_metastases": "NONE DETECTED",
      "spinal_cord": "Normal signal, no intramedullary lesions",
      "m_stage": "M0"
    },
    "csf_cytology": {
      "note": "CSF obtained at least 14 days post-operatively to avoid false positives from surgical blood",
      "result": "Negative for malignant cells"
    },
    "final_staging": {
      "t_stage": "T3a (tumor >3 cm with 4th ventricle invasion)",
      "m_stage": "M0 (no metastatic disease)",
      "chang_stage": "T3aM0",
      "risk_group": "Average risk"
    }
  },
  "surveillance_plan": {
    "brain_mri_schedule": [
      "Every 3 months for years 1-2",
      "Every 4 months for year 3",
      "Every 6 months for years 4-5",
      "Annually for years 6-10"
    ],
    "spinal_mri_schedule": [
      "Every 6 months for years 1-3",
      "Annually for years 4-5",
      "Then as clinically indicated"
    ],
    "sequences_required": "T1 pre/post gadolinium, T2/FLAIR, DWI, spectroscopy if recurrence suspected",
    "pseudoprogression_warning": "Enhancement changes in first 3-6 months post-radiation may represent pseudoprogression. Consider perfusion MRI and spectroscopy before concluding recurrence."
  }
}
```

### 11.8 Step 5: Therapeutic Discovery Engine — Novel SMO Antagonists

This is the drug discovery component, leveraging BioNeMo MolMIM and DiffDock to
generate novel Smoothened receptor antagonists optimized for pediatric safety.

```bash
# Step 5a: Define target and constraints
curl -X POST http://localhost:8510/v1/discovery/target-setup \
  -H "Content-Type: application/json" \
  -d '{
    "project_id": "AIDEN_MDB_SMO_001",
    "target_protein": {
      "name": "SMO (Smoothened, Frizzled Class Receptor)",
      "uniprot_id": "Q99835",
      "pdb_id": "5L7D",
      "binding_site": "transmembrane_domain_vismodegib_site",
      "crystal_resolution_angstrom": 3.2
    },
    "known_ligands": [
      {"name": "vismodegib", "smiles": "ClC1=CC(=CC(=C1)S(=O)(=O)C)C2=CC=C(C(=O)NC3CC3)C(=C2)Cl", "ic50_nm": 3.0, "limitation": "growth_plate_closure"},
      {"name": "sonidegib", "smiles": "CC1=CC(=CC(=C1)C(=O)NC2=CC=CC3=C2C=C(C=C3)OC(F)(F)F)NC(=O)C4CC4", "ic50_nm": 7.5, "limitation": "growth_plate_closure"},
      {"name": "GANT61", "smiles": "C(N1CCN(CC1)CC2=NC3=CC=CC=C3N=C2)C4=CC=CC=C4", "ic50_nm": 5000, "target": "GLI_not_SMO"}
    ],
    "design_constraints": {
      "molecular_weight_max": 500,
      "logp_range": [1.0, 3.5],
      "hbd_max": 3,
      "hba_max": 7,
      "tpsa_max": 90,
      "bbb_penetration": "required",
      "pediatric_safety": {
        "herg_ic50_min_um": 30,
        "hepatotoxicity_risk": "low",
        "growth_plate_toxicity": "MUST AVOID — primary design objective",
        "cyp_inhibition": "minimize"
      }
    }
  }'
```

```bash
# Step 5b: Generate novel molecules with MolMIM
curl -X POST http://localhost:8510/v1/discovery/molmim-generate \
  -H "Content-Type: application/json" \
  -d '{
    "project_id": "AIDEN_MDB_SMO_001",
    "generation_strategy": "scaffold_hopping",
    "seed_molecules": ["vismodegib", "sonidegib"],
    "num_candidates": 500,
    "optimization_objectives": [
      "maintain_smo_binding_pharmacophore",
      "eliminate_growth_plate_liability",
      "optimize_bbb_penetration",
      "minimize_herg_risk",
      "maximize_oral_bioavailability"
    ],
    "diversity_threshold": 0.4,
    "synthetic_accessibility_max": 5.0
  }'
```

```bash
# Step 5c: Score binding with DiffDock
curl -X POST http://localhost:8510/v1/discovery/diffdock-score \
  -H "Content-Type: application/json" \
  -d '{
    "project_id": "AIDEN_MDB_SMO_001",
    "protein_pdb": "5L7D",
    "candidates": "top_100_from_molmim",
    "scoring_metrics": ["binding_affinity", "pose_confidence", "interaction_fingerprint"],
    "reference_ligand": "vismodegib_crystal_pose"
  }'
```

**Expected Discovery Results:**
```json
{
  "status": "success",
  "candidates_generated": 500,
  "passed_drug_likeness_filter": 287,
  "passed_bbb_filter": 156,
  "passed_pediatric_safety_filter": 89,
  "top_candidates": [
    {
      "candidate_id": "SMO-PED-001",
      "smiles": "CC1=CC(=CC(=C1)C2=CN=C(N=C2)NC3CCCC3)F",
      "molecular_weight": 284.3,
      "logp": 2.1,
      "tpsa": 62.3,
      "hbd": 2,
      "hba": 4,
      "diffdock_score": -9.8,
      "predicted_ic50_nm": 12,
      "bbb_permeability": "high",
      "pediatric_safety_score": 0.91,
      "herg_predicted_ic50_um": 48.2,
      "hepatotoxicity_risk": "low",
      "growth_plate_liability": "minimal — lacks amide-sulfonamide pharmacophore associated with Hh-dependent growth plate effects",
      "synthetic_accessibility": 3.2,
      "key_interactions": ["Asp473 H-bond", "Arg400 salt bridge", "Phe484 pi-stacking"],
      "rationale": "Pyrimidine core replaces vismodegib benzamide while maintaining key SMO pocket interactions. Smaller MW improves BBB penetration."
    },
    {
      "candidate_id": "SMO-PED-002",
      "smiles": "O=C(NC1CC1)C2=CC3=C(C=C2)N=C(N3)C4=CC=C(F)C=C4",
      "molecular_weight": 323.4,
      "logp": 2.6,
      "tpsa": 68.1,
      "hbd": 2,
      "hba": 4,
      "diffdock_score": -10.2,
      "predicted_ic50_nm": 8,
      "bbb_permeability": "high",
      "pediatric_safety_score": 0.87,
      "herg_predicted_ic50_um": 35.1,
      "hepatotoxicity_risk": "low",
      "growth_plate_liability": "low — cyclopropyl amide retained but benzimidazole scaffold reduces off-target Hh pathway suppression in chondrocytes",
      "synthetic_accessibility": 3.8,
      "key_interactions": ["Asp473 H-bond", "Asn219 H-bond", "Trp281 hydrophobic"],
      "rationale": "Benzimidazole scaffold provides rigid binding geometry. Fluorophenyl maintains hydrophobic contact in deep pocket."
    },
    {
      "candidate_id": "SMO-PED-003",
      "smiles": "CC(NC(=O)C1=CC=C(C=C1)C2=CSC(=N2)N)C3CC3",
      "molecular_weight": 301.4,
      "logp": 1.8,
      "tpsa": 79.5,
      "hbd": 3,
      "hba": 5,
      "diffdock_score": -9.4,
      "predicted_ic50_nm": 18,
      "bbb_permeability": "moderate-high",
      "pediatric_safety_score": 0.93,
      "herg_predicted_ic50_um": 62.7,
      "hepatotoxicity_risk": "very_low",
      "growth_plate_liability": "minimal — aminothiazole eliminates chloro-sulfonamide motif",
      "synthetic_accessibility": 2.9,
      "key_interactions": ["Asp473 H-bond", "Arg400 H-bond", "Phe484 pi-stacking"],
      "rationale": "Highest pediatric safety score. Aminothiazole is well-tolerated scaffold in pediatric pharmacology. Lowest synthetic accessibility score — easiest to synthesize."
    }
  ],
  "comparison_with_approved_drugs": {
    "vismodegib_diffdock_score": -11.1,
    "vismodegib_ic50_nm": 3.0,
    "best_candidate_diffdock_score": -10.2,
    "best_candidate_predicted_ic50_nm": 8,
    "potency_gap": "2.7-fold — within acceptable range for lead optimization",
    "safety_advantage": "Candidates designed to avoid growth plate toxicity, enabling pediatric use"
  }
}
```

**Clinical Significance:** The Therapeutic Discovery Engine generates novel SMO
antagonists that maintain strong binding affinity (predicted IC50 8-18 nM vs vismodegib
3 nM) while specifically addressing the growth plate toxicity that makes approved SHH
inhibitors contraindicated in children. The top candidate (SMO-PED-002) shows only
2.7-fold lower predicted potency than vismodegib — an acceptable gap for a compound
with a fundamentally different safety profile. These candidates would require preclinical
validation in growth plate chondrocyte assays and juvenile animal models before clinical
development.

### 11.9 Step 6: Clinical Trial Agent — Pediatric Brain Tumor Trials

```bash
curl -X POST http://localhost:8538/v1/trial/match \
  -H "Content-Type: application/json" \
  -d '{
    "patient_id": "AIDEN_MDB_001",
    "diagnosis": "medulloblastoma_shh_subtype",
    "age_years": 10,
    "molecular_features": {
      "subtype": "SHH",
      "mutations": ["PTCH1_p.Arg577*"],
      "tp53_status": "wild_type",
      "mycn_status": "not_amplified"
    },
    "staging": "M0_average_risk",
    "therapy_interest": ["shh_inhibitor", "reduced_radiation", "novel_targeted"]
  }'
```

**Expected Response:**
```json
{
  "status": "success",
  "matched_trials": [
    {
      "trial_id": "NCT01878617",
      "name": "SJMB12 — Risk-Adapted Therapy for Medulloblastoma (St. Jude)",
      "phase": "II",
      "status": "Active",
      "relevance_score": 0.95,
      "key_feature": "Molecular subtype-driven risk stratification. SHH TP53-WT average-risk receives reduced-dose CSI (15 Gy) + cyclophosphamide-based chemo",
      "benefit": "Reduced CSI dose (15 Gy vs 23.4 Gy standard) — significant reduction in neurocognitive late effects",
      "sites": ["St. Jude Children's Research Hospital"]
    },
    {
      "trial_id": "NCT04402073",
      "name": "ACNS0332 Follow-on — Molecularly-Targeted Medulloblastoma",
      "phase": "II/III",
      "status": "Recruiting",
      "relevance_score": 0.90,
      "key_feature": "COG trial incorporating molecular subtyping into treatment stratification",
      "shh_arm": "Standard therapy + vismodegib for skeletally mature patients only. Under-10 patients receive standard backbone without SHH inhibitor.",
      "note": "Confirms that vismodegib is restricted to skeletally mature patients in COG protocols"
    },
    {
      "trial_id": "NCT04299113",
      "name": "Arsenic Trioxide + Radiation for SHH Medulloblastoma",
      "phase": "I/II",
      "status": "Recruiting",
      "relevance_score": 0.85,
      "key_feature": "Arsenic trioxide as GLI1 inhibitor — targets downstream of SMO, avoids growth plate concerns",
      "rationale": "Arsenic trioxide degrades GLI1/GLI2 proteins. Pediatric safety established from APL treatment. Phase I pediatric brain tumor data pending.",
      "sites": ["CHOP", "Dana-Farber", "Stanford"]
    },
    {
      "trial_id": "NCT03434262",
      "name": "Pediatric MATCH — Arm for SHH Pathway Tumors",
      "phase": "II",
      "status": "Recruiting",
      "relevance_score": 0.78,
      "key_feature": "NCI-COG Pediatric MATCH basket trial. SHH pathway arm for recurrent/refractory tumors.",
      "eligibility_caveat": "Requires prior treatment failure — not first-line. Relevant if Aiden relapses."
    },
    {
      "trial_id": "NCT05199584",
      "name": "Proton vs Photon CSI for Pediatric Medulloblastoma — Neurocognitive Outcomes",
      "phase": "III",
      "status": "Recruiting",
      "relevance_score": 0.75,
      "key_feature": "Randomized comparison of neurocognitive outcomes. Would provide level 1 evidence for proton therapy benefit.",
      "sites": ["Multi-institutional, 15+ pediatric proton centers"]
    }
  ]
}
```

### 11.10 Clinical Decision Summary for Aiden

| Decision Point | Recommendation | Evidence Level |
|---|---|---|
| Molecular subtype | SHH with PTCH1 loss, TP53-WT | WHO 2021 classification |
| Risk stratification | Average risk (M0, NTR, TP53-WT) | COG criteria |
| Surgery | Complete — 3mm residual acceptable | Neurosurgical consensus |
| Radiation | CSI 23.4 Gy + PF boost 54 Gy | COG ACNS0331 |
| Radiation modality | PROTON THERAPY strongly recommended | Yock et al., Lancet Oncol 2016 |
| Vismodegib/sonidegib | CONTRAINDICATED — growth plates open | FDA label, Robinson JCO 2015 |
| Novel SHH inhibitor | SMO-PED-002 (discovery candidate) | Preclinical — needs validation |
| Chemotherapy | Packer regimen (Cis/CCNU/VCR x 8) | COG standard |
| Trial option | SJMB12 — reduced CSI (15 Gy) | Phase II active |
| Monitoring | Brain MRI q3mo x 2y, neurocog annually | COG surveillance |

**Family communication:** Aiden's parents are told: "Aiden's tumor is the Sonic Hedgehog
subtype, which is one of the more favorable molecular groups when TP53 is normal, as it
is in Aiden's case. The standard treatment is craniospinal radiation plus chemotherapy,
with an expected cure rate of 75-85%. We strongly recommend proton therapy to reduce
long-term brain and hormone effects. There are FDA-approved drugs that target his
specific pathway (SHH inhibitors), but they cause growth plate damage and cannot be
used while Aiden is still growing. Our computational drug discovery platform has
identified novel compounds that may work without this side effect, but they need years
of testing. The most exciting option is the St. Jude SJMB12 trial, which uses a lower
radiation dose for Aiden's molecular subtype — this could significantly reduce
neurocognitive late effects. We should also plan for growth hormone monitoring starting
about a year after radiation, and Aiden may need educational support as he returns to
school."

---

## Section 11: Demo Execution Guide

### 12.1 Prerequisites

**Hardware Requirements:**
- NVIDIA GPU: DGX Spark recommended; minimum RTX 3090 (24 GB VRAM)
- System RAM: 64 GB minimum, 128 GB recommended
- Storage: ~1.1 TB for full deployment (models, data, containers)
  - BioNeMo models: ~200 GB
  - Milvus vector database: ~50 GB (3.56M vectors)
  - Reference genomes (GRCh38): ~60 GB
  - Docker images: ~150 GB
  - Working space for genomics: ~400 GB per sample pair

**Software Requirements:**
- Docker Engine 24.0+ with NVIDIA Container Toolkit
- Docker Compose v2.20+
- NVIDIA Driver 535+ with CUDA 12.x support
- Python 3.10+ (for orchestrator scripts)
- Nextflow 23.10+ (for genomics pipeline)

**API Keys and Credentials:**
- Anthropic API key (Claude) — set in `.env` as `ANTHROPIC_API_KEY`
- NVIDIA NGC API key — set in `.env` as `NGC_API_KEY`
- (Optional) ClinicalTrials.gov API key for live trial matching

**Network Requirements:**
- Outbound HTTPS (443) for Claude API calls
- All inter-service communication is local (Docker network)

### 12.2 Starting the Platform

```bash
# Clone and enter the repository
git clone https://github.com/your-org/hcls-ai-factory.git
cd hcls-ai-factory

# Configure environment
cp .env.example .env
# Edit .env with your API keys

# Start all services
docker compose -f docker-compose.dgx-spark.yml up -d

# Verify health of all services (11 services expected)
./health-monitor.sh --check-once

# Or verify individual services
curl -s http://localhost:8080/health          # Landing page
curl -s http://localhost:5000/health          # Genomic Foundation Engine
curl -s http://localhost:8540/health          # Single-Cell Agent
curl -s http://localhost:8522/health          # CAR-T Agent
curl -s http://localhost:8126/health          # Cardiology Agent
curl -s http://localhost:8532/health          # Autoimmune Agent
curl -s http://localhost:8107/health          # PGx Agent
curl -s http://localhost:8538/health          # Clinical Trial Agent
curl -s http://localhost:8527/health          # Oncology Agent
curl -s http://localhost:8528/health          # Neurology Agent
curl -s http://localhost:8524/health          # Imaging Agent
curl -s http://localhost:8510/health          # Therapeutic Discovery Engine
```

### 12.3 Running Individual Demos

**Full automated run:**
```bash
# Run all demos end-to-end (FASTQ path)
python scripts/run_demo.py --demo all --entry-point fastq

# Run a specific demo (VCF path for faster execution)
python scripts/run_demo.py --demo 5 --entry-point vcf
python scripts/run_demo.py --demo 6 --entry-point vcf

# Run with verbose output (shows all API payloads and responses)
python scripts/run_demo.py --demo 5 --entry-point vcf --verbose
```

**Manual execution (curl commands):**

Each demo section in this document contains the exact curl commands. Execute them
sequentially, using the output of each step as context for the next. For example:

```bash
# Demo 5 — CAR-T Decision (VCF path, ~30 minutes total)
# Step 1: Single-cell analysis
curl -X POST http://localhost:8540/v1/sc/query -H "Content-Type: application/json" \
  -d @payloads/demo5_step1_singlecell.json | jq .

# Step 2: CAR-T integrated assessment
curl -X POST http://localhost:8522/v1/cart/integrated-assessment \
  -H "Content-Type: application/json" \
  -d @payloads/demo5_step2_cart.json | jq .

# Continue through steps 3-6...
```

### 12.4 Expected Timing

| Demo | FASTQ Path | VCF Path | Notes |
|---|---|---|---|
| Demo 1 | ~5 hours | ~30 min | Includes full genomic alignment |
| Demo 2 | ~5 hours | ~25 min | Additional CNV calling |
| Demo 3 | ~4.5 hours | ~20 min | Germline-only path faster |
| Demo 4 | ~5 hours | ~35 min | WES + PGx annotation |
| Demo 5 | ~5.5 hours | ~30 min | Tumor + germline + scRNA-seq |
| Demo 6 | ~6 hours | ~45 min | Includes drug discovery (~15 min) |

The FASTQ-to-VCF genomic processing (Parabricks + DeepVariant) accounts for 90%+ of
the FASTQ path runtime. The VCF path bypasses this entirely.

### 12.5 Troubleshooting

| Issue | Likely Cause | Resolution |
|---|---|---|
| Service fails health check | Container not started | `docker compose logs <service>` to check errors |
| Milvus connection refused | etcd or MinIO not ready | Wait 60s after startup; check `docker compose logs milvus-standalone` |
| CUDA out of memory | Multiple GPU services | Run genomics and drug discovery sequentially, not in parallel |
| Claude API 429 | Rate limit exceeded | Add `--rate-limit 5` to run_demo.py, or wait and retry |
| DiffDock timeout | Large molecule set | Reduce candidate count or increase timeout in config |
| Genomics pipeline stalls | Insufficient disk | Ensure 400 GB free for intermediate BAM files |
| VCF annotation empty | Missing annotation databases | Run `scripts/download_annotations.sh` to fetch ClinVar, gnomAD, COSMIC |
| Agent returns generic response | Missing context from prior step | Ensure step outputs are passed as context to subsequent steps |

---

## Section 12: Clinical Validation References

The clinical claims, statistics, and treatment recommendations in this document are
grounded in peer-reviewed literature and established clinical guidelines. Below are the
primary references organized by clinical domain.

### 13.1 Pediatric Acute Lymphoblastic Leukemia

1. **Maude SL, Laetsch TW, Buechner J, et al.** Tisagenlecleucel in children and young
   adults with B-cell lymphoblastic leukemia. *N Engl J Med.* 2018;378(5):439-448.
   doi:10.1056/NEJMoa1709866. *(ELIANA trial: 82% CR, 77% CRS, 12-month EFS 66%)*

2. **Pui CH, Yang JJ, Hunger SP, et al.** Childhood acute lymphoblastic leukemia:
   progress through collaboration. *J Clin Oncol.* 2015;33(27):2938-2948.
   doi:10.1200/JCO.2014.59.1636. *(Pediatric ALL outcomes, MRD significance)*

3. **Hunger SP, Mullighan CG.** Acute lymphoblastic leukemia in children. *N Engl J Med.*
   2015;373(16):1541-1552. doi:10.1056/NEJMra1400972.
   *(Comprehensive review of pediatric ALL biology and treatment)*

4. **Gardner RA, Finney O, Annesley C, et al.** Intent-to-treat leukemia remission by
   CD19 CAR T cells of defined formulation and dose in children and young adults.
   *Blood.* 2017;129(25):3322-3331. doi:10.1182/blood-2017-02-769208.
   *(CD19 CAR-T dosing and outcomes in pediatric ALL)*

5. **Lee DW, Santomasso BD, Locke FL, et al.** ASTCT consensus grading for cytokine
   release syndrome and neurologic toxicity associated with immune effector cells.
   *Biol Blood Marrow Transplant.* 2019;25(4):625-638.
   doi:10.1016/j.bbmt.2018.12.758. *(CRS and ICANS grading criteria)*

### 13.2 CAR-T Therapy and Immunotherapy

6. **Shah NN, Fry TJ.** Mechanisms of resistance to CAR T cell therapy. *Nat Rev Clin
   Oncol.* 2019;16(6):372-385. doi:10.1038/s41571-019-0184-6.
   *(CD19 antigen loss, resistance mechanisms)*

7. **Frey NV, Porter DL.** Cytokine release syndrome with novel therapeutics for acute
   lymphoblastic leukemia. *Hematology Am Soc Hematol Educ Program.* 2016;2016(1):
   567-572. doi:10.1182/asheducation-2016.1.567. *(CRS management, tocilizumab)*

8. **Schultz LM, Baggott C, Engel B, et al.** Immune effector cell-associated
   hypotension and tachycardia in pediatric patients receiving CAR T cell therapy.
   *Transplant Cell Ther.* 2023;29(1):25-33. *(Cardiac monitoring during CAR-T)*

9. **Summers C, Wu QV, Annesley C, et al.** Hematologic toxicity after CD19-directed
   CAR T cell therapy. *Blood Adv.* 2022;6(4):1149-1158.
   *(Post-CAR-T cytopenias: 30-40% incidence, duration, management)*

### 13.3 Medulloblastoma

10. **Robinson GW, Orr BA, Wu G, et al.** Vismodegib exerts targeted efficacy against
    recurrent Sonic Hedgehog-subgroup medulloblastoma: results from phase II Pediatric
    Brain Tumor Consortium studies PBTC-025B and PBTC-032. *J Clin Oncol.*
    2015;33(24):2646-2654. doi:10.1200/JCO.2014.60.1591.
    *(Vismodegib in SHH MDB, growth plate toxicity in children)*

11. **Taylor MD, Northcott PA, Korshunov A, et al.** Molecular subgrouping of
    medulloblastoma: the current consensus. *Acta Neuropathol.* 2012;123(4):465-472.
    doi:10.1007/s00401-011-0922-z. *(WNT, SHH, Group 3, Group 4 classification)*

12. **Ramaswamy V, Remke M, Bouffet E, et al.** Risk stratification of childhood
    medulloblastoma in the molecular era: the current consensus. *Acta Neuropathol.*
    2016;131(6):821-831. doi:10.1007/s00401-016-1569-6.
    *(Molecular risk stratification, TP53 prognostic significance)*

13. **Packer RJ, Gajjar A, Vezina G, et al.** Phase III study of craniospinal radiation
    therapy followed by adjuvant chemotherapy for newly diagnosed average-risk
    medulloblastoma. *J Clin Oncol.* 2006;24(25):4202-4208.
    *(Packer regimen: standard chemotherapy backbone)*

14. **Yock TI, Yeap BY, Ebb DH, et al.** Long-term toxic effects of proton radiotherapy
    for paediatric medulloblastoma: a phase 2 single-arm study. *Lancet Oncol.*
    2016;17(3):287-298. doi:10.1016/S1470-2045(15)00167-9.
    *(Proton therapy outcomes, neurocognitive preservation, 3-5 IQ point benefit)*

### 13.4 Pharmacogenomics

15. **Relling MV, Schwab M, Whirl-Carrillo M, et al.** Clinical Pharmacogenetics
    Implementation Consortium (CPIC) guideline for thiopurine dosing based on TPMT
    and NUDT15 genotypes: 2018 update. *Clin Pharmacol Ther.* 2019;105(5):1095-1105.
    *(TPMT/NUDT15 and mercaptopurine dosing)*

16. **Caudle KE, Dunnenberger HM, Freimuth RR, et al.** Standardizing terms for clinical
    pharmacogenetic test results: consensus terms from the Clinical Pharmacogenetics
    Implementation Consortium. *Genet Med.* 2017;19(2):215-223.
    *(CPIC phenotype terminology)*

### 13.5 Genomics and Variant Interpretation

17. **Richards S, Aziz N, Bale S, et al.** Standards and guidelines for the
    interpretation of sequence variants: a joint consensus recommendation of the
    American College of Medical Genetics and Genomics and the Association for
    Molecular Pathology. *Genet Med.* 2015;17(5):405-424. doi:10.1038/gim.2015.30.
    *(ACMG/AMP variant classification: pathogenic, likely pathogenic, VUS, etc.)*

18. **Li MM, Datto M, Duncavage EJ, et al.** Standards and guidelines for the
    interpretation and reporting of sequence variants in cancer. *J Mol Diagn.*
    2017;19(1):4-23. doi:10.1016/j.jmoldx.2016.10.002.
    *(AMP/ASCO/CAP somatic variant classification)*

### 13.6 Clinical Trials and COG Protocols

19. **Children's Oncology Group (COG).** AALL0932: Treatment of patients with standard
    risk B-lymphoblastic leukemia. ClinicalTrials.gov NCT01190930.
    *(Standard-risk B-ALL backbone therapy)*

20. **Children's Oncology Group (COG).** AALL1231: Treatment of patients newly diagnosed
    with high-risk B-lymphoblastic leukemia. ClinicalTrials.gov NCT02883049.
    *(High-risk B-ALL, includes MRD-based stratification)*

21. **Children's Oncology Group (COG).** ANBL1232: Treatment of patients with high-risk
    neuroblastoma. ClinicalTrials.gov NCT01416857.
    *(Neuroblastoma intensive multimodal therapy)*

22. **Gajjar A, Robinson GW, Smith KS, et al.** Outcomes by clinical and molecular
    features in children with medulloblastoma treated with risk-adapted therapy:
    results of an international phase III trial (SJMB03). *J Clin Oncol.*
    2021;39(7):822-835. doi:10.1200/JCO.20.01372.
    *(SJMB03 results, molecular subtype-driven outcomes)*

### 13.7 Cardiac and Late Effects

23. **Lipshultz SE, Lipsitz SR, Sallan SE, et al.** Chronic progressive cardiac
    dysfunction years after doxorubicin therapy for childhood acute lymphoblastic
    leukemia. *J Clin Oncol.* 2005;23(12):2629-2636. doi:10.1200/JCO.2005.12.121.
    *(Anthracycline cardiotoxicity in pediatric ALL survivors)*

24. **Armenian SH, Hudson MM, Mulder RL, et al.** Recommendations for cardiomyopathy
    surveillance for survivors of childhood cancer: a report from the International
    Late Effects of Childhood Cancer Guideline Harmonization Group. *Lancet Oncol.*
    2015;16(3):e123-e136. *(Cardiac surveillance guidelines)*

25. **Mulrooney DA, Yeazel MW, Kawashima T, et al.** Cardiac outcomes in a cohort of
    adult survivors of childhood and adolescent cancer: retrospective analysis of the
    Childhood Cancer Survivor Study cohort. *BMJ.* 2009;339:b4606.
    *(Long-term cardiac risk in childhood cancer survivors)*

---

## Section 13: Appendix — Sample API Payloads

This appendix provides complete, copy-paste-ready JSON payloads for key API calls
referenced in each demo. These payloads can be saved as files and used with curl:

```bash
curl -X POST http://localhost:<port>/<endpoint> \
  -H "Content-Type: application/json" \
  -d @<payload_file>.json
```

### 14.1 Demo 5 Payloads

**Payload 5.1 — Single-Cell CD19 Validation (Step 1)**
```json
{
  "patient_id": "ETHAN_CART_001",
  "query": "Validate CD19 expression on B-ALL blasts, characterize tumor microenvironment, and assess MRD clone architecture",
  "data_sources": {
    "flow_cytometry": "/data/flow/ethan_mrd_panel.fcs",
    "scrna_seq": "/data/scrna/ethan_bm_10x.h5"
  },
  "analysis_type": "cart_target_validation",
  "parameters": {
    "cd_markers": ["CD19", "CD22", "CD10", "CD34", "CD20", "CD38"],
    "tme_cell_types": ["CD8_T", "CD4_T", "NK", "monocyte", "DC", "B_normal"],
    "exhaustion_markers": ["PD1", "TIM3", "LAG3", "TIGIT"],
    "clone_tracking": true,
    "minimum_cells": 5000
  }
}
```

**Payload 5.2 — CAR-T Integrated Assessment (Step 2)**
```json
{
  "patient_id": "ETHAN_CART_001",
  "indication": "relapsed_refractory_b_all",
  "age_years": 12,
  "weight_kg": 42,
  "bsa_m2": 1.32,
  "lansky_score": 80,
  "prior_therapies": [
    {
      "regimen": "AALL1231_induction",
      "response": "morphologic_cr_mrd_positive_5pct",
      "duration_months": 1
    },
    {
      "regimen": "FLAG_IDA",
      "response": "failure_12pct_blasts",
      "duration_months": 1.5
    }
  ],
  "current_disease": {
    "blast_percentage_bm": 12,
    "cd19_positive_pct": 97.2,
    "cd19_mfi": 8500,
    "cd22_positive_pct": 85.3,
    "cns_status": "CNS1",
    "extramedullary_disease": false
  },
  "organ_function": {
    "lvef_pct": 58,
    "alt_ul": 28,
    "ast_ul": 32,
    "total_bilirubin_mgdl": 0.6,
    "creatinine_mgdl": 0.5,
    "egfr_ml_min_173m2": 120,
    "pulmonary_function": "normal",
    "o2_saturation_pct": 98
  },
  "infection_status": {
    "active_infections": false,
    "cmv_igg": "positive",
    "ebv_igg": "positive",
    "hiv": "negative",
    "hepatitis_b": "negative",
    "hepatitis_c": "negative"
  },
  "donor_availability": {
    "hla_matched_sibling": true,
    "sibling_age": 15,
    "sibling_cmv_status": "positive"
  },
  "assess": [
    "eligibility",
    "crs_risk",
    "icans_risk",
    "manufacturing_plan",
    "bridging_therapy",
    "cd22_backup",
    "hsct_comparison"
  ]
}
```

**Payload 5.3 — Cardiac Assessment (Step 3)**
```json
{
  "patient_id": "ETHAN_CART_001",
  "query": "Pre-lymphodepletion cardiac risk assessment for CAR-T therapy",
  "clinical_context": {
    "age_years": 12,
    "sex": "male",
    "weight_kg": 42,
    "prior_anthracycline": {
      "agent": "daunorubicin",
      "cumulative_mg_m2": 250,
      "doxorubicin_equivalent_mg_m2": 210
    },
    "cardiac_history": {
      "prior_cardiac_events": false,
      "family_history_cardiomyopathy": false
    },
    "current_cardiac_data": {
      "lvef_pct": 58,
      "lv_fractional_shortening_pct": 32,
      "global_longitudinal_strain_pct": -19.2,
      "troponin_i_ng_ml": 0.01,
      "bnp_pg_ml": 45,
      "ecg": {
        "rhythm": "normal_sinus",
        "rate_bpm": 80,
        "intervals": {"pr_ms": 140, "qrs_ms": 88, "qtc_ms": 410},
        "abnormalities": "none"
      }
    },
    "planned_therapy": {
      "lymphodepletion": "fludarabine_30_mg_m2_x4d_cyclophosphamide_500_mg_m2_x2d",
      "cart_product": "tisagenlecleucel"
    }
  }
}
```

### 14.2 Demo 6 Payloads

**Payload 6.1 — Oncology SHH Medulloblastoma Query (Step 2)**
```json
{
  "patient_id": "AIDEN_MDB_001",
  "query": "Treatment recommendations for SHH-subtype medulloblastoma with PTCH1 loss in a 10-year-old",
  "diagnosis": {
    "tumor_type": "medulloblastoma",
    "who_grade": "IV",
    "molecular_subtype": "SHH",
    "histology": "desmoplastic_nodular",
    "key_mutations": [
      {"gene": "PTCH1", "variant": "c.1729C>T", "protein": "p.Arg577*", "type": "nonsense", "classification": "pathogenic"}
    ],
    "wild_type_genes": ["TP53", "SUFU", "SMO", "TERT_promoter"],
    "copy_number": {
      "MYCN": "not_amplified",
      "chromosome_14q": "no_loss",
      "chromosome_17p": "no_loss"
    },
    "staging": {
      "m_stage": "M0",
      "csf_cytology": "negative",
      "spinal_mri": "no_drop_metastases"
    },
    "surgery": {
      "extent": "near_total_resection",
      "residual_cm": 0.3,
      "residual_location": "floor_4th_ventricle"
    }
  },
  "patient": {
    "age_years": 10,
    "sex": "male",
    "weight_kg": 32,
    "bsa_m2": 1.12,
    "growth_plates": "open",
    "tanner_stage": "I"
  }
}
```

**Payload 6.2 — Drug Discovery Target Setup (Step 5a)**
```json
{
  "project_id": "AIDEN_MDB_SMO_001",
  "target_protein": {
    "name": "SMO (Smoothened, Frizzled Class Receptor)",
    "gene_symbol": "SMO",
    "uniprot_id": "Q99835",
    "pdb_id": "5L7D",
    "pdb_description": "Human SMO in complex with cholesterol and vismodegib analog",
    "binding_site": {
      "type": "transmembrane_domain",
      "key_residues": ["Asp473", "Arg400", "Asn219", "Phe484", "Trp281", "Val329", "Leu325"],
      "pocket_volume_angstrom3": 420
    }
  },
  "known_ligands": [
    {
      "name": "vismodegib",
      "smiles": "ClC1=CC(=CC(=C1)S(=O)(=O)C)C2=CC=C(C(=O)NC3CC3)C(=C2)Cl",
      "ic50_nm": 3.0,
      "mw": 421.3,
      "clinical_status": "FDA_approved",
      "pediatric_limitation": "irreversible_growth_plate_closure"
    },
    {
      "name": "sonidegib",
      "smiles": "CC1=CC(=CC(=C1)C(=O)NC2=CC=CC3=C2C=C(C=C3)OC(F)(F)F)NC(=O)C4CC4",
      "ic50_nm": 7.5,
      "mw": 485.5,
      "clinical_status": "FDA_approved",
      "pediatric_limitation": "irreversible_growth_plate_closure"
    }
  ],
  "design_constraints": {
    "drug_likeness": {
      "molecular_weight_max": 500,
      "logp_range": [1.0, 3.5],
      "hbd_max": 3,
      "hba_max": 7,
      "rotatable_bonds_max": 7,
      "tpsa_max_angstrom2": 90
    },
    "bbb_penetration": {
      "required": true,
      "mw_preferred_max": 450,
      "logp_preferred_range": [1.5, 3.0],
      "tpsa_preferred_max": 75,
      "hbd_preferred_max": 2
    },
    "pediatric_safety": {
      "herg_ic50_minimum_um": 30,
      "hepatotoxicity_risk_max": "low",
      "growth_plate_toxicity": "MUST_AVOID",
      "mutagenicity": "negative_predicted",
      "cyp_inhibition": {
        "CYP3A4_ic50_min_um": 10,
        "CYP2D6_ic50_min_um": 10,
        "CYP2C9_ic50_min_um": 10
      }
    },
    "synthetic_accessibility_max": 5.0,
    "pains_filter": true,
    "brenk_filter": true
  }
}
```

**Payload 6.3 — MolMIM Generation (Step 5b)**
```json
{
  "project_id": "AIDEN_MDB_SMO_001",
  "model": "molmim_v1",
  "generation_strategy": "scaffold_hopping",
  "seed_molecules": [
    {
      "name": "vismodegib",
      "smiles": "ClC1=CC(=CC(=C1)S(=O)(=O)C)C2=CC=C(C(=O)NC3CC3)C(=C2)Cl",
      "preserve_pharmacophore": ["amide_hbond_donor", "aromatic_core"]
    },
    {
      "name": "sonidegib",
      "smiles": "CC1=CC(=CC(=C1)C(=O)NC2=CC=CC3=C2C=C(C=C3)OC(F)(F)F)NC(=O)C4CC4",
      "preserve_pharmacophore": ["biphenyl_hydrophobic", "amide_hbond"]
    }
  ],
  "num_candidates": 500,
  "optimization_objectives": [
    {"objective": "binding_affinity", "weight": 0.35},
    {"objective": "bbb_penetration", "weight": 0.25},
    {"objective": "pediatric_safety", "weight": 0.25},
    {"objective": "synthetic_accessibility", "weight": 0.15}
  ],
  "diversity": {
    "tanimoto_threshold": 0.4,
    "scaffold_diversity": true,
    "max_per_scaffold": 20
  },
  "filters": {
    "drug_likeness": true,
    "pains": true,
    "brenk": true,
    "pediatric_safety": true,
    "bbb_mpo": true
  }
}
```

**Payload 6.4 — DiffDock Binding Scoring (Step 5c)**
```json
{
  "project_id": "AIDEN_MDB_SMO_001",
  "protein": {
    "pdb_id": "5L7D",
    "chain": "A",
    "binding_site_residues": ["Asp473", "Arg400", "Asn219", "Phe484", "Trp281", "Val329"],
    "prepare_protein": true
  },
  "ligands": {
    "source": "molmim_output",
    "filter": "top_100_by_composite_score",
    "include_reference": true,
    "reference_ligand": "vismodegib"
  },
  "scoring": {
    "num_poses": 10,
    "confidence_threshold": 0.5,
    "metrics": [
      "binding_affinity_kcal_mol",
      "pose_confidence_score",
      "interaction_fingerprint",
      "rmsd_to_reference_pose"
    ]
  },
  "post_scoring_filters": {
    "min_confidence": 0.7,
    "max_clash_score": 0.3,
    "required_interactions": ["Asp473_hbond"]
  },
  "output": {
    "format": ["json", "sdf"],
    "include_pose_visualization": true,
    "rank_by": "binding_affinity"
  }
}
```

**Payload 6.5 — Clinical Trial Matching (Step 6)**
```json
{
  "patient_id": "AIDEN_MDB_001",
  "diagnosis": {
    "tumor_type": "medulloblastoma",
    "molecular_subtype": "SHH",
    "histology": "desmoplastic_nodular",
    "key_mutations": ["PTCH1_p.Arg577*"],
    "tp53_status": "wild_type",
    "mycn_status": "not_amplified",
    "m_stage": "M0",
    "risk_group": "average"
  },
  "patient": {
    "age_years": 10,
    "sex": "male",
    "ecog_ps": 1,
    "prior_treatment": ["surgery"],
    "planned_treatment": ["craniospinal_radiation", "chemotherapy"],
    "growth_plates": "open"
  },
  "search_criteria": {
    "therapy_types": [
      "shh_pathway_inhibitor",
      "reduced_intensity_radiation",
      "novel_targeted_agent",
      "immunotherapy"
    ],
    "phase": ["I", "II", "III"],
    "status": ["recruiting", "active_not_recruiting"],
    "age_eligible": true,
    "distance_miles_max": 500,
    "patient_location_zip": "37203"
  }
}
```

---

*End of document — Demos 3-6, Execution Guide, Clinical Validation References, and Appendix*
*Generated by HCLS AI Factory v1.0 on DGX Spark*
