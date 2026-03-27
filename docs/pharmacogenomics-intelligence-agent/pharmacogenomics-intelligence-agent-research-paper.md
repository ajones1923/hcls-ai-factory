# From Genome to Safe Prescription: A Multi-Collection RAG Architecture for Clinical Pharmacogenomic Decision Support

**Author:** Adam Jones
**Date:** March 2026
**Version:** 0.1.0 (Pre-Implementation)
**License:** Apache 2.0

Part of the HCLS AI Factory -- an end-to-end precision medicine platform.
https://github.com/ajones1923/hcls-ai-factory

---

## Abstract

Adverse drug reactions (ADRs) are the fourth leading cause of death in the United States, responsible for approximately 106,000 deaths and 2.2 million hospitalizations annually, at a cost exceeding $136 billion per year. Yet an estimated 95-99% of patients carry at least one actionable pharmacogenomic (PGx) variant that would alter prescribing decisions -- and fewer than 1% of prescriptions in the U.S. are informed by genetic testing. This catastrophic gap between genomic knowledge and clinical practice persists because pharmacogenomic data is complex, rapidly evolving, fragmented across multiple guideline bodies (CPIC, DPWG, FDA, CPNDS), and difficult for non-specialist clinicians to interpret and act upon at the point of prescribing.

This paper presents the architectural design, clinical rationale, and product requirements for the Pharmacogenomics Intelligence Agent -- a clinical decision support system built on multi-collection retrieval-augmented generation (RAG) that transforms raw genomic data from the HCLS AI Factory's genomics pipeline (VCF output) into actionable, patient-specific prescribing guidance for over 400 drug-gene interactions across 25+ pharmacogenes. The agent will unify 14 specialized Milvus vector collections spanning pharmacogene reference data (star allele definitions, diplotype-to-phenotype mappings, activity scores), clinical guideline knowledge (CPIC Level A/B guidelines, DPWG recommendations, FDA Table of Pharmacogenomic Biomarkers), drug interaction intelligence (PharmGKB annotations, drug-drug-gene interactions, phenoconversion modeling), HLA-mediated hypersensitivity screening (12 HLA-drug associations with absolute contraindications), population pharmacokinetics (ethnicity-adjusted allele frequencies, dosing nomograms), clinical evidence (published PGx implementation studies, outcomes data), and the shared genomic evidence collection (3.5 million variants) -- enabling queries like "What are the prescribing implications of this patient's CYP2D6 *4/*41 genotype for their current medication list?" or "Does this patient carry any HLA alleles that contraindicate specific drugs before surgery?"

The system extends the proven multi-collection RAG architecture established by six existing intelligence agents in the HCLS AI Factory (Precision Biomarker, Precision Oncology, CAR-T, Imaging, Autoimmune, and Cardiology), adapting it with a genomic variant-to-drug mapping pipeline capable of processing whole-genome VCF files, star allele calling for CYP450 enzymes and transporters, diplotype-to-phenotype translation using CPIC standardized terms, multi-gene interaction modeling (e.g., CYP2C9 + VKORC1 for warfarin), phenoconversion detection (drug-induced CYP inhibition altering metabolizer status), HLA typing from NGS data for hypersensitivity screening, and real-time medication list cross-referencing against the patient's complete PGx profile. Eight reference clinical workflows will cover the highest-impact prescribing scenarios: pre-emptive PGx panel interpretation, opioid prescribing safety (CYP2D6/codeine/tramadol), anticoagulant optimization (CYP2C9/VKORC1/warfarin), antidepressant selection (CYP2D6/CYP2C19/SSRIs/TCAs), statin myopathy risk (SLCO1B1), chemotherapy toxicity prevention (DPYD/5-FU, TPMT/NUDT15/thiopurines), HLA-mediated hypersensitivity screening (abacavir/carbamazepine/allopurinol), and polypharmacy drug-drug-gene interaction resolution.

The agent will deploy on a single NVIDIA DGX Spark ($4,699) using BGE-small-en-v1.5 embeddings (384-dimensional, IVF_FLAT, COSINE), Claude Sonnet 4.6 for evidence synthesis, and shared NVIDIA NIM microservices for on-device inference. Licensed under Apache 2.0, the platform will democratize access to pharmacogenomic intelligence that currently requires multi-million-dollar institutional PGx implementation programs -- bringing the prescribing safety of world-class pharmacogenomics centers to any clinic, pharmacy, or emergency department worldwide.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [The Pharmacogenomic Implementation Gap](#2-the-pharmacogenomic-implementation-gap)
3. [Clinical Landscape and Market Analysis](#3-clinical-landscape-and-market-analysis)
4. [Existing HCLS AI Factory Architecture](#4-existing-hcls-ai-factory-architecture)
5. [Pharmacogenomics Agent Architecture](#5-pharmacogenomics-agent-architecture)
6. [Genomic Variant-to-Drug Mapping Pipeline](#6-genomic-variant-to-drug-mapping-pipeline)
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

### 1.1 The Adverse Drug Reaction Crisis

Adverse drug reactions represent one of the most significant -- and most preventable -- causes of morbidity and mortality in modern medicine. The scope of the problem is staggering:

- **Deaths:** 106,000 Americans die annually from ADRs, making them the 4th leading cause of death (Lazarou et al., JAMA)
- **Hospitalizations:** 2.2 million ADR-related hospitalizations per year in the U.S. alone
- **Cost:** $136 billion annually in direct healthcare costs (more than cardiovascular disease or diabetes management)
- **Emergency visits:** ADRs account for 27% of all emergency department drug-related visits
- **ICU admissions:** 10-20% of ICU admissions are ADR-related
- **Global impact:** The WHO estimates ADRs cause 197,000 deaths annually in the European Union

What makes this crisis particularly tragic is that a substantial proportion of ADRs are genetically predictable. Pharmacogenomic variants -- inherited differences in genes encoding drug-metabolizing enzymes, transporters, receptors, and immune molecules -- directly influence how individuals respond to medications. A patient who is a CYP2D6 poor metabolizer cannot activate codeine into morphine and will receive no pain relief. A patient who is a CYP2D6 ultra-rapid metabolizer will convert codeine too quickly, potentially causing fatal respiratory depression. A patient carrying HLA-B*57:01 who receives abacavir will develop a life-threatening hypersensitivity reaction in approximately 50% of cases. These are not rare edge cases -- they are common genetic variants present in significant proportions of the population.

### 1.2 The Promise of Pharmacogenomics

Pharmacogenomics -- the study of how genetic variation affects drug response -- has been one of the most successful translational applications of the Human Genome Project. Over the past two decades, the field has produced:

- **CPIC guidelines:** 27 gene-drug guidelines covering 80+ drug-gene pairs with actionable prescribing recommendations, graded by evidence strength (Level A = strong evidence, Level B = moderate evidence)
- **FDA labeling:** 450+ drugs carry pharmacogenomic information in their FDA-approved labeling, including 80+ with boxed warnings or contraindications based on genetic status
- **PharmGKB:** Over 780 clinical annotations, 180 drug label annotations, 150 clinical guideline annotations, and 700+ variant annotations
- **DPWG (Dutch Pharmacogenetics Working Group):** 120+ gene-drug therapeutic recommendations used across European healthcare systems
- **Economic evidence:** Multiple studies demonstrate cost-effectiveness of PGx testing, with return on investment (ROI) of $4-$13 per $1 invested across various healthcare systems

The clinical genes with the strongest evidence include:

| Gene | Function | Key Drugs Affected | Population Impact |
|------|----------|-------------------|-------------------|
| CYP2D6 | Phase I metabolism | Codeine, tramadol, tamoxifen, SSRIs, TCAs, antipsychotics | ~7% poor metabolizers, ~5% ultra-rapid |
| CYP2C19 | Phase I metabolism | Clopidogrel, PPIs, SSRIs, voriconazole | ~2-15% poor metabolizers (varies by ethnicity) |
| CYP2C9 | Phase I metabolism | Warfarin, phenytoin, NSAIDs | ~1-3% poor metabolizers |
| VKORC1 | Warfarin target | Warfarin | ~37% carry dose-reduction variant |
| SLCO1B1 | Hepatic transporter | Statins (simvastatin, atorvastatin) | ~15% carry myopathy risk variant |
| DPYD | Pyrimidine catabolism | 5-FU, capecitabine | ~3-5% carry partial deficiency |
| TPMT/NUDT15 | Thiopurine metabolism | Azathioprine, 6-MP, thioguanine | ~10% intermediate, ~0.3% poor |
| HLA-B*57:01 | Immune presentation | Abacavir | ~6-8% in Caucasians |
| HLA-B*58:01 | Immune presentation | Allopurinol | ~1-6% (varies by ethnicity) |
| HLA-A*31:01 | Immune presentation | Carbamazepine | ~2-5% in Caucasians |
| HLA-B*15:02 | Immune presentation | Carbamazepine, phenytoin | ~8% in Southeast Asian populations |
| UGT1A1 | Phase II conjugation | Irinotecan, atazanavir | ~10% poor metabolizers |
| G6PD | Redox protection | Rasburicase, primaquine, dapsone | ~8% globally (up to 25% in some populations) |
| CYP3A5 | Phase I metabolism | Tacrolimus | ~70-90% non-expressers (varies by ethnicity) |
| IFNL3 (IL28B) | Immune response | PEG-IFN/ribavirin | ~70% favorable genotype in Europeans |

### 1.3 Why Pharmacogenomics Remains Unused

Despite overwhelming evidence, PGx adoption remains dismally low:

- **<1% of prescriptions** in the U.S. are guided by pharmacogenomic testing
- **Only 15%** of physicians report feeling comfortable interpreting PGx results
- **42%** of physicians surveyed stated they had never ordered a PGx test
- **88%** of medical schools provide fewer than 8 hours of PGx education across 4 years
- **<5%** of health systems have integrated PGx into their electronic health records

The barriers are well-characterized:

1. **Knowledge gap:** Clinicians lack training in PGx interpretation. A CYP2D6 *4/*41 diplotype means nothing to most prescribers.
2. **Complexity:** Star allele nomenclature, activity scores, diplotype-to-phenotype translations, gene-drug-drug interactions, and phenoconversion are intellectually demanding even for specialists.
3. **Fragmentation:** CPIC, DPWG, FDA, and institutional guidelines sometimes disagree. Clinicians don't know which to follow.
4. **EHR integration:** Most EHRs cannot store, display, or trigger alerts from structured PGx data. Results are often returned as PDFs that sit unread in the chart.
5. **Point-of-care timing:** PGx results must be available at the moment of prescribing, not days later from a reference lab.
6. **Multi-gene complexity:** Real patients have variants in multiple PGx genes simultaneously. A patient on warfarin needs CYP2C9 + VKORC1 + CYP4F2 interpreted together. A patient on psychiatric medications may need CYP2D6 + CYP2C19 + CYP1A2 + CYP3A4.
7. **Population diversity:** Allele frequencies vary dramatically across ethnic groups. The most common CYP2D6 poor metabolizer allele in Europeans (*4) is rare in East Asians, where *10 predominates.

### 1.4 The HCLS AI Factory Opportunity

The HCLS AI Factory's existing genomics pipeline already generates the raw data needed for pharmacogenomic analysis. Every patient genome processed through the pipeline produces a VCF file containing all pharmacogenomic variants -- but this data is not currently translated into prescribing guidance. The Pharmacogenomics Intelligence Agent closes this gap by:

1. **Extracting PGx variants** from VCF files produced by the genomics pipeline (Parabricks/DeepVariant)
2. **Calling star alleles** using standardized nomenclature (PharmVar database)
3. **Translating diplotypes to phenotypes** using CPIC activity score algorithms
4. **Cross-referencing** the patient's PGx profile against their current and potential medication list
5. **Generating actionable prescribing recommendations** grounded in CPIC/DPWG guidelines
6. **Storing results** in persistent Milvus collections for lifetime re-querying as new drugs are prescribed

The agent doesn't just report what variants a patient carries -- it answers the question every prescriber actually needs answered: **"Is this drug safe for this patient, and if not, what should I prescribe instead?"**

---

## 2. The Pharmacogenomic Implementation Gap

### 2.1 Preventable Deaths: The Scale of the Problem

To understand why the Pharmacogenomics Intelligence Agent matters, consider these real-world scenarios that occur daily in hospitals and clinics worldwide:

**Scenario 1: Codeine and CYP2D6**
A 3-year-old child undergoes tonsillectomy. Prescribed codeine for post-operative pain. The child is a CYP2D6 ultra-rapid metabolizer (alleles *1/*1xN), converting codeine to morphine at 4-8x the normal rate. The child develops fatal respiratory depression and dies. This exact scenario led to an FDA Black Box Warning in 2013 -- yet codeine is still prescribed to children without CYP2D6 testing.

**Scenario 2: Clopidogrel and CYP2C19**
A 58-year-old man receives a coronary stent and is prescribed clopidogrel (Plavix) as standard antiplatelet therapy. He is a CYP2C19 poor metabolizer (*2/*2), producing functionally no active metabolite. One month later, he suffers stent thrombosis and dies. An alternative antiplatelet (prasugrel or ticagrelor) -- not dependent on CYP2C19 -- would have prevented this death. The FDA label has carried a boxed warning about CYP2C19 since 2010.

**Scenario 3: Abacavir and HLA-B*57:01**
A 34-year-old woman newly diagnosed with HIV is prescribed an abacavir-containing regimen. She is HLA-B*57:01 positive. Two weeks later, she develops abacavir hypersensitivity syndrome -- fever, rash, malaise progressing to hypotension and organ failure. Before mandatory HLA-B*57:01 screening (implemented ~2008), this reaction occurred in ~5% of abacavir-treated patients and was fatal in some cases. This is one of the few PGx tests with near-universal adoption, proving that when guidelines are clear and testing is mandated, PGx saves lives.

**Scenario 4: 5-Fluorouracil and DPYD**
A 62-year-old woman with colon cancer begins adjuvant 5-FU chemotherapy. She carries DPYD*2A (splice site variant), rendering her unable to metabolize 5-FU. She develops grade 4 mucositis, pancytopenia, and sepsis. She dies from treatment toxicity, not from cancer. Pre-treatment DPYD testing with dose adjustment would have prevented this death. The European Medicines Agency now mandates DPYD testing before fluoropyrimidine therapy; the FDA does not.

### 2.2 The Implementation Gap by the Numbers

The disconnect between PGx knowledge and clinical practice can be quantified:

| Metric | Current State | Ideal State | Gap |
|--------|--------------|------------|-----|
| Prescriptions guided by PGx | <1% | 30-50% (for PGx-relevant drugs) | 30-50x |
| Health systems with PGx CDS | <5% | 100% | 20x |
| Time to PGx result availability | 3-14 days (reference lab) | Pre-emptive (already in chart) | Paradigm shift |
| Physicians comfortable with PGx | 15% | >80% | 5x |
| Pharmacogenes routinely tested | 1-2 (reactive) | 12-25 (pre-emptive panel) | 10x |
| PGx variants detected per genome | ~0 (not analyzed) | 15-50 actionable findings | Infinite gap |
| Drug-gene alerts at prescribing | Near zero | Every PGx-relevant prescription | Total gap |

### 2.3 The Pre-emptive vs. Reactive Testing Paradigm

The field is shifting from reactive PGx testing (test one gene when prescribing one drug) to pre-emptive PGx testing (test a panel of pharmacogenes BEFORE any drug is needed, store results for life). This paradigm shift is essential because:

1. **Genomic data is stable:** Unlike lab values, a patient's PGx genotype never changes. Test once, use forever.
2. **Reactive testing is too slow:** When a patient needs pain medication in the ER, there is no time to wait 7 days for CYP2D6 results.
3. **Pre-emptive panels are cost-effective:** A $250-$500 multi-gene panel tested once replaces dozens of single-gene tests ($100-$300 each) over a lifetime.
4. **EHR alerts require pre-existing data:** Clinical decision support can only fire at the moment of prescribing if PGx results are already in the system.

The Pharmacogenomics Intelligence Agent is architected for the pre-emptive paradigm: it processes a patient's entire genome once, extracts all pharmacogenomic variants, and stores the complete PGx profile in persistent Milvus collections. Every subsequent prescribing query -- whether today or in 20 years -- can be answered instantly.

---

## 3. Clinical Landscape and Market Analysis

### 3.1 Market Size and Growth

The pharmacogenomics market is experiencing rapid growth driven by declining sequencing costs, regulatory mandates, and health system adoption:

- **Global PGx market (2025):** $4.1 billion
- **Projected (2030):** $11.2 billion (CAGR 22.1%)
- **PGx testing volume growth:** 30-40% annually
- **DTC PGx testing:** $800 million (2025), growing 25% annually
- **PGx clinical decision support software:** $420 million (2025), growing 35% annually

### 3.2 Key Market Drivers

1. **Regulatory mandates:** EMA requires DPYD testing before fluoropyrimidines (2020). FDA adding PGx to more drug labels annually. CMS considering PGx coverage expansion.
2. **Health system initiatives:** 80+ U.S. health systems have PGx implementation programs (IGNITE Network, CPIC institutions). Mayo Clinic, St. Jude, Vanderbilt, University of Florida leading adoption.
3. **Payer coverage:** UnitedHealthcare, Cigna, and Aetna now cover multi-gene PGx panels for select indications. Medicare MAC coverage expanding.
4. **Pharmacist-led models:** Clinical pharmacists are driving PGx adoption through medication therapy management programs.
5. **Legal liability:** Failure to test PGx before prescribing drugs with known genetic interactions is becoming a malpractice concern (Marchetti v. United States, 2022).

### 3.3 Competitive Landscape Overview

| Company | Product | Approach | Limitations |
|---------|---------|----------|-------------|
| OneOme | RightMed | Panel test + CDS portal | No genomics integration, proprietary |
| Myriad Genetics | GeneSight | Psychiatric PGx panel | Narrow scope (psych only), criticized evidence base |
| Invitae/Tempus | PGx panels | Testing + basic CDS | No RAG, limited multi-gene interaction |
| Clinical Pharmacogenomics Implementation Consortium | CPIC guidelines | Gold-standard guidelines (free) | Text-based, no CDS integration |
| Translational Software | PGx CDS | Standalone CDS engine | No genomic pipeline integration |
| Color Health | PGx program | Employer-based testing | Limited clinical depth |

**Gap our agent fills:** No existing product combines (1) whole-genome pharmacogenomic extraction from a genomics pipeline, (2) multi-collection RAG over CPIC/DPWG/FDA guidelines, (3) multi-gene interaction modeling, (4) phenoconversion detection, (5) HLA-mediated hypersensitivity screening, and (6) natural language clinical queries -- all running on a $4,699 local device with no cloud data exposure.

---

## 4. Existing HCLS AI Factory Architecture

### 4.1 Three-Stage Pipeline

The HCLS AI Factory is an end-to-end precision medicine platform deployed on NVIDIA DGX Spark. Its three-stage pipeline provides the foundational infrastructure that the Pharmacogenomics Intelligence Agent builds upon:

**Stage 1: Genomics Pipeline** (`genomics-pipeline/`)
- Input: FASTQ (raw sequencing data) or pre-aligned BAM
- Processing: BWA-MEM2 alignment → Parabricks/DeepVariant variant calling
- Output: Annotated VCF files with 11.7 million variants per genome
- Relevance to PGx: **VCF output contains all pharmacogenomic variants** but they are not currently extracted, interpreted, or translated to prescribing guidance

**Stage 2: Precision Intelligence Network** (`rag-chat-pipeline/`)
- Milvus vector database (19530) with BGE-small-en-v1.5 embeddings
- Claude Sonnet 4.6 for evidence synthesis
- 3.5 million searchable genomic variant vectors in shared `genomic_evidence` collection
- Multi-collection architecture proven across 11 agents

**Stage 3: Therapeutic Discovery Engine** (`drug-discovery-pipeline/`)
- BioNeMo MolMIM for molecular generation
- DiffDock for binding pose prediction
- RDKit for ADMET property calculation
- Relevance to PGx: Drug metabolism predictions complement PGx by modeling how genetic variants affect drug pharmacokinetics at the molecular level

### 4.2 Existing Intelligence Agents

Six intelligence agents currently operate within the HCLS AI Factory, each demonstrating the multi-collection RAG architecture that the Pharmacogenomics Agent will extend:

| Agent | Collections | Port (UI/API) | Key Capability |
|-------|-------------|---------------|----------------|
| Precision Biomarker | 13 + shared | 8502/8102 | Biomarker interpretation, includes existing PGx module (71K lines) |
| Precision Oncology | 10 + shared | 8503/8103 | Tumor genomics, therapy selection |
| CAR-T Intelligence | 10 + shared | 8504/8104 | CAR-T construct design, manufacturing |
| Imaging Intelligence | 10 + shared | 8505/8105 | Medical imaging AI, radiology CDS |
| Precision Autoimmune | 14 + shared | 8506/8106 | Diagnostic odyssey acceleration, clinical document intelligence |
| Cardiology Intelligence | 12 + shared | 8526/8527 | Cardiac risk, ECG/imaging interpretation |

The Pharmacogenomics Intelligence Agent will be assigned ports **8507 (Streamlit UI)** and **8107 (FastAPI API)**.

### 4.3 Relationship to Existing Biomarker Agent PGx Module

The Precision Biomarker Agent already contains a substantial pharmacogenomics module (`src/pharmacogenomics.py`, 71,246 bytes) implementing CPIC guidelines for 13 genes. The Pharmacogenomics Intelligence Agent is NOT a replacement -- it is a **dedicated, deep-dive expansion** that provides:

| Capability | Biomarker Agent PGx Module | Pharmacogenomics Intelligence Agent |
|-----------|---------------------------|-------------------------------------|
| Genes covered | 13 (CPIC Level 1A) | 25+ (CPIC Level A, B, and C + DPWG + FDA) |
| Drug-gene pairs | ~80 | 400+ |
| Multi-gene interactions | None | Full modeling (e.g., CYP2C9+VKORC1+CYP4F2 for warfarin) |
| Phenoconversion | None | Drug-drug-gene interaction detection |
| HLA screening | 3 alleles | 12+ alleles with population-specific frequencies |
| Star allele calling | Pre-computed lookup | VCF-to-star-allele pipeline |
| Population adjustments | None | Ethnicity-adjusted allele frequencies and dosing |
| Milvus collections | 1 (pgx_rules) | 14 specialized collections |
| Clinical workflows | 1 (basic PGx query) | 8 comprehensive workflows |
| Dosing calculators | None | Warfarin dosing algorithm, tacrolimus dosing, 5-FU dose adjustment |

The two agents complement each other: the Biomarker Agent provides quick PGx screening as part of a broader biomarker analysis, while the Pharmacogenomics Agent provides deep clinical pharmacogenomic consultation when detailed PGx guidance is needed.

---

## 5. Pharmacogenomics Agent Architecture

### 5.1 System Design

```
┌─────────────────────────────────────────────────────────┐
│                    PHARMACOGENOMICS                       │
│                  INTELLIGENCE AGENT                       │
│                                                           │
│  ┌──────────┐  ┌──────────────┐  ┌────────────────────┐ │
│  │ Streamlit │  │ FastAPI      │  │ VCF-to-PGx         │ │
│  │ UI :8507  │  │ API :8107    │  │ Pipeline           │ │
│  └────┬─────┘  └──────┬───────┘  └────────┬───────────┘ │
│       │               │                    │             │
│  ┌────▼───────────────▼────────────────────▼──────────┐ │
│  │              PGx Intelligence Core                   │ │
│  │                                                       │ │
│  │  ┌─────────────┐ ┌──────────────┐ ┌──────────────┐  │ │
│  │  │ Star Allele  │ │ Phenotype    │ │ Drug-Gene    │  │ │
│  │  │ Caller       │ │ Translator   │ │ Matcher      │  │ │
│  │  └──────┬──────┘ └──────┬───────┘ └──────┬───────┘  │ │
│  │         │               │                │           │ │
│  │  ┌──────▼──────┐ ┌──────▼───────┐ ┌──────▼───────┐  │ │
│  │  │ Multi-Gene   │ │ Phenoconv.   │ │ HLA          │  │ │
│  │  │ Interaction  │ │ Detector     │ │ Screener     │  │ │
│  │  └──────┬──────┘ └──────┬───────┘ └──────┬───────┘  │ │
│  │         └───────────────┼────────────────┘           │ │
│  │                         │                             │ │
│  │  ┌──────────────────────▼────────────────────────┐   │ │
│  │  │         Multi-Collection RAG Engine            │   │ │
│  │  │    (14 PGx collections + shared genomic)       │   │ │
│  │  └──────────────────────┬────────────────────────┘   │ │
│  │                         │                             │ │
│  │  ┌──────────────────────▼────────────────────────┐   │ │
│  │  │      Claude Sonnet 4.6 Synthesis Engine        │   │ │
│  │  │  (Grounded PGx recommendations with evidence)  │   │ │
│  │  └────────────────────────────────────────────────┘   │ │
│  └───────────────────────────────────────────────────────┘ │
│                                                           │
│  ┌───────────────────────────────────────────────────────┐ │
│  │                  Milvus (19530)                        │ │
│  │  14 PGx collections + shared genomic_evidence         │ │
│  └───────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────┘
```

### 5.2 Core Processing Modules

The agent contains six core processing modules:

**1. Star Allele Caller**
Extracts pharmacogenomic variants from VCF files and resolves them to star allele nomenclature. Uses PharmVar-defined haplotype tables for each gene. Handles:
- SNV-based star alleles (e.g., CYP2D6 *4 = 1846G>A)
- Structural variants (CYP2D6 gene deletion *5, gene duplication *1xN, *2xN)
- Hybrid alleles (CYP2D6/CYP2D7 hybrids like *36, *13)
- Suballeles (CYP2D6 *4.001 vs. *4.002)

**2. Phenotype Translator**
Converts diplotypes to standardized phenotype terms using CPIC's activity score system:
- CYP enzymes: Ultra-rapid, Rapid, Normal, Intermediate, Poor Metabolizer
- Transporters: Increased, Normal, Decreased, Poor Function
- HLA: Positive (carrier) vs. Negative
- Enzyme deficiency: Normal, Intermediate, Deficient (G6PD, DPYD, TPMT)

**3. Drug-Gene Matcher**
Cross-references the patient's phenotype profile against their current medication list. For each drug-gene interaction, returns:
- CPIC recommendation level (A, B, C, D)
- Clinical action (standard dosing, dose adjustment, alternative drug, avoid, contraindicated)
- Alert level (INFO, WARNING, CRITICAL)
- Alternative drug suggestions with rationale

**4. Multi-Gene Interaction Engine**
Models complex scenarios where multiple pharmacogenes affect a single drug or therapeutic area:
- Warfarin: CYP2C9 + VKORC1 + CYP4F2 → personalized dose calculator
- Psychiatric medications: CYP2D6 + CYP2C19 + CYP1A2 + CYP3A4 → comprehensive metabolizer profile
- Immunosuppressants: CYP3A5 + CYP3A4 + ABCB1 → tacrolimus dosing

**5. Phenoconversion Detector**
Identifies when a patient's medication list contains CYP inhibitors or inducers that change their effective metabolizer status:
- Example: A CYP2D6 normal metabolizer taking fluoxetine (strong CYP2D6 inhibitor) is phenoconverted to a CYP2D6 poor metabolizer. If codeine is then prescribed, it will be ineffective despite a "normal" genotype.
- Example: A CYP3A4 normal metabolizer taking carbamazepine (strong CYP3A4 inducer) metabolizes cyclosporine too rapidly, potentially causing organ rejection.

**6. HLA Screener**
Screens for HLA alleles associated with drug hypersensitivity reactions:

| HLA Allele | Drug | Reaction | Risk if Positive | CPIC Level |
|-----------|------|----------|------------------|------------|
| HLA-B*57:01 | Abacavir | Hypersensitivity syndrome | ~50% | A (mandatory screening) |
| HLA-B*58:01 | Allopurinol | SJS/TEN | OR 80-580 | A |
| HLA-B*15:02 | Carbamazepine | SJS/TEN | OR 1,357 | A |
| HLA-A*31:01 | Carbamazepine | DRESS/maculopapular | OR 25.9 | A |
| HLA-B*15:02 | Phenytoin | SJS/TEN | OR 6.7 | A |
| HLA-B*15:02 | Oxcarbazepine | SJS/TEN | OR 27.9 | B |
| HLA-B*13:01 | Dapsone | Hypersensitivity | OR 20.5 | B |
| HLA-A*33:03 | Ticlopidine | Hepatotoxicity | OR 13 | B |
| HLA-DRB1*07:01 | Lapatinib | Hepatotoxicity | OR 2.5 | B |
| HLA-B*35:01 | Minocycline | DRESS | Under study | C |
| HLA-B*38:02 | Sulfasalazine | DRESS | Under study | C |
| HLA-DPB1*03:01 | Aspirin | AERD | Under study | C |

---

## 6. Genomic Variant-to-Drug Mapping Pipeline

### 6.1 Pipeline Overview

```
VCF File (11.7M variants)
    │
    ▼
┌──────────────────────────┐
│ 1. PGx Variant Extraction │
│    Filter to ~2,500 PGx   │
│    relevant positions      │
└─────────┬────────────────┘
          │
          ▼
┌──────────────────────────┐
│ 2. Star Allele Resolution │
│    PharmVar haplotype      │
│    tables for each gene    │
└─────────┬────────────────┘
          │
          ▼
┌──────────────────────────┐
│ 3. Diplotype Assembly     │
│    Phase-aware diplotype   │
│    determination           │
└─────────┬────────────────┘
          │
          ▼
┌──────────────────────────┐
│ 4. Phenotype Translation  │
│    CPIC activity scores    │
│    → standardized terms    │
└─────────┬────────────────┘
          │
          ▼
┌──────────────────────────┐
│ 5. Drug Matching          │
│    Cross-reference med     │
│    list → recommendations  │
└─────────┬────────────────┘
          │
          ▼
┌──────────────────────────┐
│ 6. Report Generation      │
│    Clinical action items   │
│    with evidence grading   │
└──────────────────────────┘
```

### 6.2 Variant Extraction Detail

From a typical whole-genome VCF with 11.7 million variants, the pipeline extracts approximately 2,500 pharmacogenomically relevant positions across 25+ genes. Extraction uses a curated BED file of PGx-relevant coordinates derived from:

- PharmVar database (all defined star allele positions)
- CPIC gene-specific tables (defining variants for each guideline)
- ClinVar PGx-classified variants
- PharmGKB high-confidence variant annotations

### 6.3 Star Allele Calling Algorithm

Star allele calling is the most complex step and varies by gene:

**Simple genes (SNV-only):** CYP2C19, CYP2C9, VKORC1, DPYD, TPMT, NUDT15
- Haplotype matching against PharmVar-defined allele tables
- Each star allele defined by 1-5 SNVs
- Diplotype by standard phasing (statistical or read-backed)

**Complex genes (SNV + structural):** CYP2D6
- CYP2D6 is the most pharmacogenomically important and most difficult gene to genotype
- Requires: SNV calling, copy number determination (gene deletions *5, duplications *1xN, *2xN), hybrid allele detection (CYP2D6/CYP2D7 fusions)
- The agent implements a tiered approach:
  - Tier 1: SNV-based star allele assignment from VCF
  - Tier 2: Read depth analysis for copy number (requires BAM)
  - Tier 3: Hybrid allele detection (requires specialized alignment)
- When structural variant data is unavailable, the agent flags this limitation and provides recommendations based on SNV data alone with appropriate caveats

**HLA genes:** HLA-A, HLA-B, HLA-DRB1, HLA-DPB1
- HLA typing from WGS data using established tools (OptiType, HLA-HD)
- Four-digit resolution sufficient for most PGx associations
- Population-specific frequency annotations

### 6.4 Activity Score System

CPIC uses activity scores to standardize phenotype assignment across genes. Example for CYP2D6:

| Allele | Activity Score | Classification |
|--------|---------------|----------------|
| *1 (normal function) | 1.0 | Functional |
| *2 (normal function) | 1.0 | Functional |
| *9 (decreased function) | 0.5 | Reduced |
| *10 (decreased function) | 0.25 | Reduced |
| *17 (decreased function) | 0.5 | Reduced |
| *41 (decreased function) | 0.5 | Reduced |
| *4 (no function) | 0 | Non-functional |
| *5 (gene deletion) | 0 | Non-functional |
| *6 (no function) | 0 | Non-functional |

**Diplotype to Phenotype:**
- Activity Score ≥2.25: Ultra-rapid Metabolizer (UM)
- Activity Score 1.25-2.25: Normal Metabolizer (NM)
- Activity Score 0.25-1.25: Intermediate Metabolizer (IM)
- Activity Score 0: Poor Metabolizer (PM)

Example: CYP2D6 *4/*41 → Activity Score 0 + 0.5 = 0.5 → **Intermediate Metabolizer**

---

## 7. Milvus Collection Design

### 7.1 Collection Overview

The Pharmacogenomics Intelligence Agent maintains 14 specialized Milvus collections plus access to the shared genomic evidence collection:

| # | Collection Name | Record Estimate | Purpose |
|---|----------------|----------------|---------|
| 1 | `pgx_gene_reference` | ~5,000 | Star allele definitions, activity scores, allele frequencies |
| 2 | `pgx_drug_guidelines` | ~8,000 | CPIC/DPWG/FDA guideline recommendations |
| 3 | `pgx_drug_interactions` | ~12,000 | Drug-gene interaction annotations (PharmGKB) |
| 4 | `pgx_hla_hypersensitivity` | ~2,000 | HLA-drug hypersensitivity associations |
| 5 | `pgx_phenoconversion` | ~3,000 | CYP inhibitor/inducer drug interactions |
| 6 | `pgx_dosing_algorithms` | ~1,500 | Population PK models, dosing nomograms |
| 7 | `pgx_clinical_evidence` | ~15,000 | Published PGx implementation studies, outcomes |
| 8 | `pgx_population_data` | ~4,000 | Ethnicity-specific allele frequencies, dosing adjustments |
| 9 | `pgx_clinical_trials` | ~6,000 | PGx-related clinical trial data |
| 10 | `pgx_fda_labels` | ~3,500 | FDA pharmacogenomic labeling information |
| 11 | `pgx_drug_alternatives` | ~5,000 | Alternative drug recommendations by metabolizer status |
| 12 | `pgx_patient_profiles` | Variable | Patient-specific PGx profiles (per-patient) |
| 13 | `pgx_implementation` | ~4,000 | Health system PGx implementation protocols |
| 14 | `pgx_education` | ~2,500 | Clinician-facing PGx education materials |
| S | `genomic_evidence` | 3,500,000 | Shared genomic variant evidence (read-only) |

### 7.2 Collection Schemas

**Collection 1: `pgx_gene_reference`**
Core pharmacogene reference data including star allele definitions and allele function.

```python
FieldSchema(name="id", dtype=VARCHAR, is_primary=True, max_length=100)
FieldSchema(name="embedding", dtype=FLOAT_VECTOR, dim=384)
FieldSchema(name="gene", dtype=VARCHAR, max_length=20)          # CYP2D6, SLCO1B1, etc.
FieldSchema(name="star_allele", dtype=VARCHAR, max_length=30)    # *1, *4, *10, *41
FieldSchema(name="defining_variants", dtype=VARCHAR, max_length=500)  # rs numbers
FieldSchema(name="activity_score", dtype=FLOAT)                  # 0, 0.25, 0.5, 1.0
FieldSchema(name="function_status", dtype=VARCHAR, max_length=30) # Normal, Decreased, No function
FieldSchema(name="allele_frequency_global", dtype=FLOAT)
FieldSchema(name="allele_frequency_european", dtype=FLOAT)
FieldSchema(name="allele_frequency_african", dtype=FLOAT)
FieldSchema(name="allele_frequency_east_asian", dtype=FLOAT)
FieldSchema(name="allele_frequency_south_asian", dtype=FLOAT)
FieldSchema(name="allele_frequency_latino", dtype=FLOAT)
FieldSchema(name="pharmvar_id", dtype=VARCHAR, max_length=30)
FieldSchema(name="text_chunk", dtype=VARCHAR, max_length=3000)
FieldSchema(name="source", dtype=VARCHAR, max_length=100)
```

**Collection 2: `pgx_drug_guidelines`**
CPIC, DPWG, and FDA guideline recommendations for drug-gene pairs.

```python
FieldSchema(name="id", dtype=VARCHAR, is_primary=True, max_length=100)
FieldSchema(name="embedding", dtype=FLOAT_VECTOR, dim=384)
FieldSchema(name="gene", dtype=VARCHAR, max_length=20)
FieldSchema(name="drug", dtype=VARCHAR, max_length=50)
FieldSchema(name="phenotype", dtype=VARCHAR, max_length=40)      # Poor Metabolizer, etc.
FieldSchema(name="guideline_body", dtype=VARCHAR, max_length=10) # CPIC, DPWG, FDA
FieldSchema(name="cpic_level", dtype=VARCHAR, max_length=5)      # A, A/B, B, C, D
FieldSchema(name="recommendation", dtype=VARCHAR, max_length=1000)
FieldSchema(name="clinical_action", dtype=VARCHAR, max_length=30) # STANDARD, DOSE_ADJUST, AVOID, etc.
FieldSchema(name="alert_level", dtype=VARCHAR, max_length=10)    # INFO, WARNING, CRITICAL
FieldSchema(name="alternative_drugs", dtype=VARCHAR, max_length=500)
FieldSchema(name="dose_adjustment", dtype=VARCHAR, max_length=200)
FieldSchema(name="evidence_pmids", dtype=VARCHAR, max_length=300)
FieldSchema(name="guideline_version", dtype=VARCHAR, max_length=20)
FieldSchema(name="last_updated", dtype=VARCHAR, max_length=10)
FieldSchema(name="text_chunk", dtype=VARCHAR, max_length=3000)
```

**Collection 3: `pgx_drug_interactions`**
Comprehensive drug-gene interaction annotations from PharmGKB.

```python
FieldSchema(name="id", dtype=VARCHAR, is_primary=True, max_length=100)
FieldSchema(name="embedding", dtype=FLOAT_VECTOR, dim=384)
FieldSchema(name="drug", dtype=VARCHAR, max_length=50)
FieldSchema(name="gene", dtype=VARCHAR, max_length=20)
FieldSchema(name="variant_rsid", dtype=VARCHAR, max_length=20)
FieldSchema(name="interaction_type", dtype=VARCHAR, max_length=30) # PK, PD, efficacy, toxicity
FieldSchema(name="effect_description", dtype=VARCHAR, max_length=500)
FieldSchema(name="evidence_level", dtype=VARCHAR, max_length=5)  # 1A, 1B, 2A, 2B, 3, 4
FieldSchema(name="clinical_significance", dtype=VARCHAR, max_length=20) # Actionable, Informative
FieldSchema(name="pharmgkb_id", dtype=VARCHAR, max_length=30)
FieldSchema(name="affected_phenotype", dtype=VARCHAR, max_length=100)
FieldSchema(name="text_chunk", dtype=VARCHAR, max_length=3000)
```

**Collection 4: `pgx_hla_hypersensitivity`**
HLA-mediated drug hypersensitivity associations with population-specific risks.

```python
FieldSchema(name="id", dtype=VARCHAR, is_primary=True, max_length=100)
FieldSchema(name="embedding", dtype=FLOAT_VECTOR, dim=384)
FieldSchema(name="hla_allele", dtype=VARCHAR, max_length=20)
FieldSchema(name="drug", dtype=VARCHAR, max_length=50)
FieldSchema(name="reaction_type", dtype=VARCHAR, max_length=50)  # SJS/TEN, DRESS, HSR, hepatotoxicity
FieldSchema(name="risk_if_positive", dtype=VARCHAR, max_length=100) # OR, absolute risk %
FieldSchema(name="severity", dtype=VARCHAR, max_length=20)       # Life-threatening, Severe, Moderate
FieldSchema(name="cpic_level", dtype=VARCHAR, max_length=5)
FieldSchema(name="recommendation", dtype=VARCHAR, max_length=500)
FieldSchema(name="screening_mandatory", dtype=BOOL)
FieldSchema(name="prevalence_european", dtype=FLOAT)
FieldSchema(name="prevalence_african", dtype=FLOAT)
FieldSchema(name="prevalence_east_asian", dtype=FLOAT)
FieldSchema(name="prevalence_south_asian", dtype=FLOAT)
FieldSchema(name="prevalence_latino", dtype=FLOAT)
FieldSchema(name="alternative_drugs", dtype=VARCHAR, max_length=300)
FieldSchema(name="text_chunk", dtype=VARCHAR, max_length=3000)
```

**Collection 5: `pgx_phenoconversion`**
Drug-drug-gene interactions where concomitant medications alter metabolizer phenotype.

```python
FieldSchema(name="id", dtype=VARCHAR, is_primary=True, max_length=100)
FieldSchema(name="embedding", dtype=FLOAT_VECTOR, dim=384)
FieldSchema(name="affected_enzyme", dtype=VARCHAR, max_length=20)  # CYP2D6, CYP3A4, etc.
FieldSchema(name="precipitant_drug", dtype=VARCHAR, max_length=50) # The inhibitor or inducer
FieldSchema(name="interaction_type", dtype=VARCHAR, max_length=15) # Strong/Moderate/Weak inhibitor or inducer
FieldSchema(name="effect_on_phenotype", dtype=VARCHAR, max_length=100) # "Converts NM to PM"
FieldSchema(name="clinical_significance", dtype=VARCHAR, max_length=500)
FieldSchema(name="affected_substrate_drugs", dtype=VARCHAR, max_length=500)
FieldSchema(name="time_to_onset", dtype=VARCHAR, max_length=50)
FieldSchema(name="reversibility", dtype=VARCHAR, max_length=50)
FieldSchema(name="evidence_level", dtype=VARCHAR, max_length=5)
FieldSchema(name="text_chunk", dtype=VARCHAR, max_length=3000)
```

**Collection 6: `pgx_dosing_algorithms`**
Pharmacokinetic models and dosing nomograms for dose-critical drugs.

```python
FieldSchema(name="id", dtype=VARCHAR, is_primary=True, max_length=100)
FieldSchema(name="embedding", dtype=FLOAT_VECTOR, dim=384)
FieldSchema(name="drug", dtype=VARCHAR, max_length=50)
FieldSchema(name="genes_involved", dtype=VARCHAR, max_length=100)
FieldSchema(name="algorithm_name", dtype=VARCHAR, max_length=100) # e.g., "IWPC Warfarin Dosing"
FieldSchema(name="input_variables", dtype=VARCHAR, max_length=500)
FieldSchema(name="formula_description", dtype=VARCHAR, max_length=1000)
FieldSchema(name="validation_cohort", dtype=VARCHAR, max_length=200)
FieldSchema(name="accuracy_metrics", dtype=VARCHAR, max_length=200) # R², MAE, etc.
FieldSchema(name="clinical_context", dtype=VARCHAR, max_length=500)
FieldSchema(name="text_chunk", dtype=VARCHAR, max_length=3000)
```

**Collections 7-14** follow similar patterns for clinical evidence, population data, clinical trials, FDA labels, drug alternatives, patient profiles, implementation protocols, and education materials. Each uses the standard 384-dimension BGE-small-en-v1.5 embedding with IVF_FLAT COSINE indexing.

### 7.3 Collection Search Weights

```python
WEIGHT_GENE_REFERENCE = 0.10
WEIGHT_DRUG_GUIDELINES = 0.14      # Highest -- clinical guidelines are primary
WEIGHT_DRUG_INTERACTIONS = 0.12
WEIGHT_HLA_HYPERSENSITIVITY = 0.10
WEIGHT_PHENOCONVERSION = 0.08
WEIGHT_DOSING_ALGORITHMS = 0.07
WEIGHT_CLINICAL_EVIDENCE = 0.08
WEIGHT_POPULATION_DATA = 0.06
WEIGHT_CLINICAL_TRIALS = 0.04
WEIGHT_FDA_LABELS = 0.06
WEIGHT_DRUG_ALTERNATIVES = 0.05
WEIGHT_PATIENT_PROFILES = 0.03
WEIGHT_IMPLEMENTATION = 0.02
WEIGHT_EDUCATION = 0.02
WEIGHT_GENOMIC_EVIDENCE = 0.03
```

---

## 8. Clinical Workflows

### 8.1 Workflow 1: Pre-emptive PGx Panel Interpretation

**Trigger:** New patient genome processed through HCLS AI Factory genomics pipeline.

**Process:**
1. VCF file ingested by PGx variant extraction module
2. Star alleles called for all 25+ pharmacogenes
3. Diplotypes assembled and phenotypes translated
4. Complete PGx profile stored in `pgx_patient_profiles` collection
5. Profile cross-referenced against common medication classes
6. "PGx Passport" report generated with:
   - All metabolizer phenotypes (table format)
   - High-risk drug-gene interactions (CRITICAL alerts)
   - Pre-emptive recommendations for common drug classes
   - HLA hypersensitivity risks
   - Population-adjusted allele frequency context

**Demo Query:** *"Generate a complete pharmacogenomic profile for this patient and identify any high-priority drug-gene interactions."*

**Output Example:**
```
PHARMACOGENOMIC PROFILE SUMMARY

Gene        Diplotype    Phenotype              Key Implications
CYP2D6      *4/*41       Intermediate Met.      Reduced codeine activation, consider alternatives
CYP2C19     *1/*2        Intermediate Met.      Reduced clopidogrel activation -- use prasugrel
CYP2C9      *1/*3        Intermediate Met.      Reduced warfarin metabolism -- lower starting dose
VKORC1      -1639 G/A    Intermediate Sens.     Combined with CYP2C9 IM → warfarin ~3 mg/day
SLCO1B1     *1/*5        Intermediate Function  Simvastatin myopathy risk ↑ -- use pravastatin
DPYD        *1/*1        Normal                 Standard fluoropyrimidine dosing
TPMT        *1/*3A       Intermediate Met.      Reduce thiopurine dose 30-50%
HLA-B*57:01 Negative                            Abacavir safe
HLA-B*58:01 Negative                            Allopurinol safe
HLA-B*15:02 Negative                            Carbamazepine SJS risk low

CRITICAL ALERTS (2):
⚠ CYP2C19 *1/*2: AVOID clopidogrel -- use prasugrel or ticagrelor
⚠ CYP2D6 *4/*41: AVOID codeine/tramadol -- use morphine or oxycodone
```

### 8.2 Workflow 2: Opioid Prescribing Safety

**Trigger:** Opioid prescription initiated or planned for patient with PGx data.

**Process:**
1. Retrieve patient's CYP2D6 diplotype and phenotype
2. Check for CYP2D6 phenoconversion (concomitant inhibitors: fluoxetine, paroxetine, bupropion, duloxetine, terbinafine)
3. Determine effective phenotype (genetic + phenoconversion)
4. Map to opioid-specific recommendations:
   - Codeine: PM → avoid (no activation), UM → avoid (rapid activation, toxicity risk)
   - Tramadol: PM → avoid (no activation), UM → avoid (seizure/serotonin risk)
   - Hydrocodone: PM → reduced efficacy, UM → increased effect
   - Oxycodone: Minimal CYP2D6 dependence → safer alternative for PM/UM
5. Generate recommendation with alternative analgesics

**Demo Query:** *"This patient needs post-surgical pain management. What opioids are safe given their CYP2D6 status?"*

### 8.3 Workflow 3: Anticoagulant Optimization

**Trigger:** Warfarin initiation or dose adjustment for patient with PGx data.

**Process:**
1. Retrieve CYP2C9 + VKORC1 + CYP4F2 genotypes
2. Apply IWPC (International Warfarin Pharmacogenetics Consortium) dosing algorithm:
   - Inputs: CYP2C9 genotype, VKORC1 genotype, age, weight, height, race, amiodarone use, smoker status
   - Output: Predicted therapeutic dose (mg/day)
3. Cross-reference with current INR values and clinical context
4. Flag drug-drug interactions that affect warfarin metabolism
5. Generate personalized dose recommendation with evidence level

**Dosing Algorithm (IWPC):**
```
Predicted weekly dose (mg) =
    5.6044
    - 0.2614 × age (decades)
    + 0.0087 × height (cm)
    + 0.0128 × weight (kg)
    - 0.8677 × VKORC1 A/G
    - 1.6974 × VKORC1 A/A
    - 0.5211 × CYP2C9 *1/*2
    - 0.9357 × CYP2C9 *1/*3
    - 1.0616 × CYP2C9 *2/*2
    - 1.9206 × CYP2C9 *2/*3
    - 2.3312 × CYP2C9 *3/*3
    - 0.2188 × CYP4F2 *1/*3
    - 0.2760 × CYP4F2 *3/*3
    + 1.1816 × race (African American)
    - 0.1070 × race (Asian)
    - 0.2029 × amiodarone (yes)
    + 0.2107 × smoker (yes)
```

### 8.4 Workflow 4: Antidepressant Selection

**Trigger:** SSRI/SNRI/TCA initiation for depression or anxiety.

**Process:**
1. Retrieve CYP2D6 and CYP2C19 phenotypes (primary metabolizers for most antidepressants)
2. Check CYP1A2 status if fluvoxamine or clozapine considered
3. Map each candidate antidepressant to patient's metabolizer profile:
   - SSRIs metabolized primarily by CYP2D6: fluoxetine, paroxetine, fluvoxamine
   - SSRIs metabolized primarily by CYP2C19: citalopram, escitalopram, sertraline
   - TCAs metabolized by both: amitriptyline, nortriptyline, imipramine
4. Flag drug-drug interactions with other psychiatric medications
5. Rank antidepressants by suitability for patient's genotype
6. Generate recommendation with dose adjustments if needed

**Demo Query:** *"Patient is CYP2D6 poor metabolizer, CYP2C19 normal metabolizer. Which SSRI should I start for generalized anxiety?"*

**Expected Output:** *"Avoid paroxetine and fluoxetine (primarily CYP2D6 metabolized -- toxicity risk in PM). Sertraline or escitalopram (primarily CYP2C19 metabolized) are preferred. Standard dosing appropriate given CYP2C19 NM status. Note: fluoxetine is also a strong CYP2D6 inhibitor, which would compound the poor metabolizer phenotype."*

### 8.5 Workflow 5: Statin Myopathy Risk Assessment

**Trigger:** Statin prescription initiated or patient reports muscle symptoms.

**Process:**
1. Retrieve SLCO1B1 genotype (rs4149056, T>C)
2. Assess statin-specific myopathy risk:
   - SLCO1B1 *5/*5 (CC): Simvastatin myopathy risk 18% (vs. 0.6% in *1/*1)
   - SLCO1B1 *1/*5 (TC): Simvastatin myopathy risk 3%
3. Cross-reference with concurrent medications affecting statin levels
4. Recommend statin selection and dosing:
   - High risk: Avoid simvastatin >20 mg; consider pravastatin or rosuvastatin (not SLCO1B1 dependent)
   - Intermediate risk: Simvastatin ≤20 mg or alternative statin
5. If patient on existing statin with muscle symptoms, assess whether PGx explains the ADR

### 8.6 Workflow 6: Chemotherapy Toxicity Prevention

**Trigger:** Fluoropyrimidine (5-FU, capecitabine) or thiopurine (azathioprine, 6-MP) initiation.

**Process:**
1. **Fluoropyrimidines (DPYD):**
   - Retrieve DPYD genotype for 4 key variants: *2A (splice), *13 (missense), c.2846A>T, HapB3
   - Activity score calculation: each variant assigned 0 (no function) or 0.5 (decreased function)
   - Dose recommendation:
     - Activity score 2.0 (NM): Full dose
     - Activity score 1.5 (IM): Reduce dose 50%
     - Activity score 1.0 (IM): Reduce dose 50%, consider therapeutic drug monitoring
     - Activity score 0.5: Strongly reduce dose or avoid
     - Activity score 0 (PM): AVOID fluoropyrimidines -- alternative regimen required

2. **Thiopurines (TPMT + NUDT15):**
   - Retrieve TPMT and NUDT15 diplotypes
   - Combined phenotype assessment (both genes contribute independently)
   - Dose recommendation per CPIC:
     - Both NM: Full dose
     - TPMT IM + NUDT15 NM: Reduce dose 30-50%
     - TPMT PM or NUDT15 PM: Reduce dose 90% or avoid

### 8.7 Workflow 7: HLA-Mediated Hypersensitivity Screening

**Trigger:** Prescription of a drug with known HLA-mediated hypersensitivity risk.

**Process:**
1. Drug triggers HLA screening alert
2. Retrieve patient's HLA typing from genomic data
3. Cross-reference against HLA-drug hypersensitivity database
4. Generate risk assessment:
   - **CONTRAINDICATED** if positive for mandatory-screening alleles (HLA-B*57:01/abacavir)
   - **HIGH RISK** with alternative recommendation if positive for strongly associated alleles
   - **LOW RISK** if negative for relevant alleles
5. Population-specific frequency context (e.g., HLA-B*15:02 is rare in Europeans but common in Southeast Asians -- screening is ethnicity-guided for carbamazepine)

### 8.8 Workflow 8: Polypharmacy Drug-Drug-Gene Interaction Resolution

**Trigger:** Patient on 5+ medications, new drug being added, or comprehensive medication review.

**Process:**
1. Retrieve patient's complete PGx profile (all 25+ genes)
2. Cross-reference entire medication list against PGx profile
3. Identify drug-drug-gene interactions:
   - Drug A is a CYP2D6 inhibitor + Drug B is a CYP2D6 substrate + Patient is CYP2D6 IM → effective PM phenotype → Drug B toxicity risk
   - Drug C is a CYP3A4 inducer + Drug D is a CYP3A4 substrate → Drug D subtherapeutic levels
4. Model phenoconversion cascade effects
5. Prioritize interactions by clinical severity
6. Generate comprehensive medication safety report with actionable recommendations

**Demo Query:** *"This patient is on 12 medications. Review the complete medication list against their PGx profile and identify all drug-drug-gene interactions."*

---

## 9. Cross-Modal Integration and Genomic Correlation

### 9.1 VCF-to-Clinical Integration

The Pharmacogenomics Agent uniquely bridges the gap between genomic data and clinical action. Unlike other agents that operate primarily on text-based clinical knowledge, this agent directly processes structured genomic data (VCF format) and translates it to structured clinical recommendations.

Integration points with other HCLS AI Factory components:

- **Genomics Pipeline → PGx Agent:** VCF output feeds directly into the PGx variant extraction module. Every genome processed through the pipeline automatically generates a PGx profile.
- **PGx Agent → Biomarker Agent:** PGx findings inform biomarker interpretation (e.g., a CYP2D6 PM patient on tamoxifen will have subtherapeutic endoxifen levels, affecting breast cancer biomarker interpretation).
- **PGx Agent → Oncology Agent:** Chemotherapy drug selection informed by DPYD, TPMT, NUDT15, UGT1A1 status.
- **PGx Agent → Autoimmune Agent:** Biologic therapy metabolism (CYP3A4/5 for some biologics), thiopurine dosing for autoimmune conditions (azathioprine for lupus, IBD).
- **PGx Agent → Cardiology Agent:** Warfarin dosing, clopidogrel selection, statin safety.

### 9.2 Temporal Considerations

Unlike most pharmacogenomic data (which is static -- genotype doesn't change), the PGx agent must handle temporal dynamics:

1. **Phenoconversion is temporal:** A patient's effective metabolizer status changes when CYP inhibitors/inducers are started or stopped.
2. **Guidelines evolve:** CPIC updates guidelines every 2-3 years. New drug-gene associations are published continuously.
3. **Medication lists change:** The same PGx profile produces different recommendations as drugs are added or removed.
4. **Age-related changes:** Some drug metabolism changes with age (reduced CYP activity in elderly), modifying PGx-based dose recommendations.

---

## 10. NIM Integration Strategy

### 10.1 On-Device Inference

The agent leverages NVIDIA NIM microservices for computationally intensive on-device tasks:

- **BGE-small-en-v1.5 NIM:** Text embedding for all collection searches (384-dim, optimized for NVIDIA GPU)
- **Re-ranking NIM:** Cross-encoder re-ranking of search results for improved retrieval accuracy
- **Genomics NIMs:** Parabricks-based variant calling provides input VCF data

### 10.2 Cloud LLM Integration

Claude Sonnet 4.6 handles evidence synthesis and natural language response generation:

- System prompt includes PGx-specific clinical reasoning framework
- Structured output format for drug recommendations (gene, phenotype, drug, action, evidence level)
- Citation grounding to CPIC/DPWG guideline versions
- Uncertainty quantification for novel drug-gene combinations

**Privacy architecture:** Patient genomic data and PGx profiles remain on the local DGX Spark. Only anonymized query text and retrieved evidence snippets are sent to the cloud LLM. No patient identifiers, genotypes, or medication lists leave the local device.

---

## 11. Knowledge Graph Design

### 11.1 Core Entity Dictionaries

The PGx knowledge graph contains 8 core entity types:

**1. PHARMACOGENES (25+ entries)**
```python
PHARMACOGENES = {
    "CYP2D6": {
        "full_name": "Cytochrome P450 2D6",
        "chromosome": "22q13.2",
        "function": "Phase I oxidative metabolism",
        "substrates_count": 80,
        "percent_drugs_metabolized": 25,
        "star_alleles_defined": 140,
        "key_variants": ["*1", "*2", "*3", "*4", "*5", "*6", "*9", "*10", "*17", "*41"],
        "structural_variation": True,
        "complexity_level": "Very High",
        "cpic_guidelines": ["codeine", "tramadol", "tamoxifen", "ondansetron", "SSRIs", "TCAs"],
    },
    "CYP2C19": { ... },
    "CYP2C9": { ... },
    # ... 22 more genes
}
```

**2. METABOLIZER_PHENOTYPES**
```python
METABOLIZER_PHENOTYPES = {
    "Ultra-rapid Metabolizer": {
        "abbreviation": "UM",
        "clinical_meaning": "Metabolizes drug faster than normal. May need higher dose or different drug.",
        "risk": "Subtherapeutic drug levels, treatment failure (prodrugs: toxicity from rapid activation)",
    },
    "Normal Metabolizer": {
        "abbreviation": "NM",
        "clinical_meaning": "Standard drug metabolism. Standard dosing appropriate.",
        "risk": "None -- standard of care",
    },
    "Intermediate Metabolizer": {
        "abbreviation": "IM",
        "clinical_meaning": "Reduced drug metabolism. May need dose reduction.",
        "risk": "Elevated drug levels, increased ADR risk",
    },
    "Poor Metabolizer": {
        "abbreviation": "PM",
        "clinical_meaning": "Severely reduced or absent metabolism. Drug accumulates or prodrug fails to activate.",
        "risk": "Toxicity (active drugs) or treatment failure (prodrugs)",
    },
}
```

**3. DRUG_CATEGORIES (12 therapeutic areas)**
```python
DRUG_CATEGORIES = {
    "opioids": ["codeine", "tramadol", "hydrocodone", "oxycodone", "morphine"],
    "anticoagulants": ["warfarin", "clopidogrel", "prasugrel", "ticagrelor"],
    "antidepressants": ["fluoxetine", "paroxetine", "sertraline", "citalopram", "escitalopram",
                        "amitriptyline", "nortriptyline", "imipramine", "venlafaxine", "duloxetine"],
    "antipsychotics": ["aripiprazole", "haloperidol", "risperidone", "clozapine"],
    "statins": ["simvastatin", "atorvastatin", "rosuvastatin", "pravastatin", "lovastatin"],
    "chemotherapy": ["5-fluorouracil", "capecitabine", "irinotecan", "tamoxifen",
                     "azathioprine", "6-mercaptopurine", "thioguanine", "cisplatin"],
    "anticonvulsants": ["carbamazepine", "oxcarbazepine", "phenytoin", "lamotrigine", "valproate"],
    "antivirals": ["abacavir", "efavirenz", "atazanavir"],
    "immunosuppressants": ["tacrolimus", "cyclosporine", "mycophenolate", "azathioprine"],
    "cardiovascular": ["metoprolol", "propranolol", "verapamil", "amiodarone"],
    "proton_pump_inhibitors": ["omeprazole", "pantoprazole", "lansoprazole"],
    "anti_gout": ["allopurinol", "febuxostat"],
}
```

**4. CYP_INHIBITORS_INDUCERS**
```python
CYP_INHIBITORS = {
    "CYP2D6": {
        "strong": ["fluoxetine", "paroxetine", "bupropion", "quinidine", "terbinafine"],
        "moderate": ["duloxetine", "sertraline", "diphenhydramine", "abiraterone"],
        "weak": ["citalopram", "escitalopram", "amiodarone"],
    },
    "CYP3A4": {
        "strong": ["ketoconazole", "itraconazole", "clarithromycin", "ritonavir", "cobicistat"],
        "moderate": ["fluconazole", "erythromycin", "diltiazem", "verapamil", "grapefruit"],
        "weak": ["cimetidine"],
    },
    "CYP2C19": {
        "strong": ["fluoxetine", "fluvoxamine", "ticlopidine"],
        "moderate": ["omeprazole", "esomeprazole", "voriconazole"],
        "weak": ["cimetidine"],
    },
    "CYP1A2": {
        "strong": ["fluvoxamine", "ciprofloxacin", "enoxacin"],
        "moderate": ["oral contraceptives", "mexiletine"],
    },
}

CYP_INDUCERS = {
    "CYP3A4": {
        "strong": ["rifampin", "carbamazepine", "phenytoin", "St. John's wort", "phenobarbital"],
        "moderate": ["efavirenz", "bosentan", "modafinil"],
    },
    "CYP1A2": {
        "strong": ["smoking (tobacco)", "charcoal-grilled meats"],
        "moderate": ["omeprazole (high-dose)"],
    },
    "CYP2C19": {
        "strong": ["rifampin"],
        "moderate": ["carbamazepine", "efavirenz"],
    },
}
```

**5-8.** Additional dictionaries for POPULATION_ALLELE_FREQUENCIES, DRUG_ALTERNATIVE_MAPS, DOSING_PARAMETERS, and EVIDENCE_GRADING_CRITERIA follow similar structured patterns.

### 11.2 Query Expansion Maps

```python
QUERY_EXPANSION = {
    "warfarin": ["coumadin", "anticoagulant", "INR", "blood thinner", "CYP2C9", "VKORC1",
                 "vitamin K antagonist", "bleeding risk", "dose adjustment"],
    "codeine": ["opioid", "pain", "CYP2D6", "morphine", "prodrug", "ultra-rapid",
                "poor metabolizer", "respiratory depression", "analgesic"],
    "clopidogrel": ["plavix", "antiplatelet", "stent", "CYP2C19", "prasugrel", "ticagrelor",
                    "stent thrombosis", "ACS", "PCI"],
    "statin": ["simvastatin", "atorvastatin", "SLCO1B1", "myopathy", "rhabdomyolysis",
               "cholesterol", "pravastatin", "rosuvastatin", "muscle pain"],
    "tamoxifen": ["breast cancer", "CYP2D6", "endoxifen", "ER-positive", "SERM",
                  "aromatase inhibitor", "poor metabolizer"],
    "5-FU": ["fluorouracil", "capecitabine", "DPYD", "DPD deficiency", "mucositis",
             "neutropenia", "chemotherapy toxicity", "dose reduction"],
    "abacavir": ["HIV", "HLA-B*57:01", "hypersensitivity", "antiretroviral",
                 "immune-mediated reaction"],
    "carbamazepine": ["tegretol", "epilepsy", "HLA-B*15:02", "HLA-A*31:01", "SJS",
                      "TEN", "DRESS", "Stevens-Johnson"],
    "tacrolimus": ["immunosuppressant", "transplant", "CYP3A5", "organ rejection",
                   "calcineurin inhibitor", "dose adjustment"],
    # ... 15+ more drug expansion maps
}
```

---

## 12. Query Expansion and Retrieval Strategy

### 12.1 Multi-Stage Retrieval

The RAG engine uses a four-stage retrieval strategy optimized for pharmacogenomic queries:

**Stage 1: Query Classification**
Classify incoming query into PGx workflow type:
- Gene-specific query ("What does CYP2D6 *4/*41 mean?")
- Drug-specific query ("Is codeine safe for this patient?")
- Patient profile query ("Generate PGx report for this genome")
- Interaction query ("Any drug-drug-gene interactions in this med list?")
- Dosing query ("What warfarin dose for this genotype?")

**Stage 2: Targeted Collection Search**
Based on classification, prioritize relevant collections:
- Gene query → `pgx_gene_reference` (weight boost), `pgx_drug_guidelines`
- Drug query → `pgx_drug_guidelines` (weight boost), `pgx_drug_interactions`, `pgx_drug_alternatives`
- Profile query → All collections, emphasizing `pgx_gene_reference` and `pgx_drug_guidelines`
- Interaction query → `pgx_phenoconversion` (weight boost), `pgx_drug_interactions`

**Stage 3: Evidence Merging and Re-ranking**
- Merge results across collections with weighted scoring
- Cross-encoder re-ranking for relevance
- Deduplication of overlapping evidence
- Evidence grading annotation (CPIC Level A vs. B vs. emerging)

**Stage 4: LLM Synthesis**
- Structured prompt with retrieved evidence
- Clinical reasoning framework (gene → phenotype → drug → recommendation → evidence)
- Output format enforced: actionable recommendations with evidence grading

---

## 13. API and UI Design

### 13.1 FastAPI Endpoints (Port 8107)

```
Health and Status:
GET  /health                      → Service health, collection counts
GET  /collections                 → Collection names and record counts
GET  /metrics                     → Prometheus-compatible metrics

Core PGx Queries:
POST /v1/pgx/profile             → Generate complete PGx profile from VCF
POST /v1/pgx/query               → Natural language PGx query with RAG
POST /v1/pgx/drug-check          → Check single drug against PGx profile
POST /v1/pgx/medication-review   → Full medication list review
POST /v1/pgx/dosing              → Genotype-guided dosing calculation
POST /v1/pgx/hla-screen          → HLA hypersensitivity screening

Interaction Analysis:
POST /v1/pgx/interactions        → Drug-drug-gene interaction analysis
POST /v1/pgx/phenoconversion     → Phenoconversion detection

Reporting:
POST /v1/pgx/report              → Generate clinical PGx report (PDF/JSON)
POST /v1/pgx/passport            → Generate PGx Passport card
GET  /v1/pgx/profile/{patient_id} → Retrieve stored patient PGx profile

Data Management:
POST /v1/pgx/ingest-vcf          → Ingest VCF and extract PGx variants
POST /v1/pgx/update-guidelines   → Update CPIC/DPWG guideline data
GET  /v1/pgx/gene/{gene_name}    → Gene-specific reference information
GET  /v1/pgx/drug/{drug_name}    → Drug-specific PGx information
```

### 13.2 Streamlit UI Design (Port 8507)

**Tab 1: PGx Dashboard**
- Patient PGx profile overview (metabolizer status table)
- Active medication list with drug-gene interaction alerts
- Risk severity heatmap (genes × drugs)

**Tab 2: Drug Check**
- Enter drug name → instant PGx safety assessment
- Visual traffic light system: Green (safe), Yellow (adjust), Red (avoid)
- Alternative drug suggestions with rationale

**Tab 3: Medication Review**
- Paste or upload complete medication list
- Comprehensive drug-drug-gene interaction matrix
- Phenoconversion detection and cascade effects
- Prioritized action items

**Tab 4: Warfarin Dosing**
- Interactive dosing calculator (IWPC algorithm)
- Input: genotypes, demographics, concurrent medications
- Output: predicted therapeutic dose with confidence interval

**Tab 5: Chemotherapy Safety**
- DPYD screening results and 5-FU dose recommendation
- TPMT/NUDT15 results and thiopurine dose recommendation
- UGT1A1 results and irinotecan dose recommendation

**Tab 6: HLA Screening**
- Complete HLA typing results
- Drug hypersensitivity risk table
- Population-specific frequency context

**Tab 7: PGx Report Generator**
- Generate clinical PGx report (PDF)
- Generate PGx Passport (wallet card format)
- Generate provider letter (for specialist communication)

**Tab 8: Evidence Explorer**
- Search CPIC/DPWG guidelines by gene or drug
- View published PGx implementation evidence
- Clinical trial results for PGx-guided therapy

**Tab 9: Phenoconversion Modeler**
- Interactive tool: add/remove drugs and see phenotype changes
- Visual diagram of CYP enzyme inhibition/induction cascades
- "What if" scenarios for medication changes

**Tab 10: Population Analytics**
- Allele frequency explorer by ethnicity
- Metabolizer phenotype distribution charts
- Health equity analysis (differential PGx risk by population)

---

## 14. Clinical Decision Support Engines

### 14.1 Dosing Calculators

**Warfarin Dose Calculator**
- Algorithm: IWPC (International Warfarin Pharmacogenetics Consortium)
- Inputs: CYP2C9, VKORC1, CYP4F2 genotypes + demographics
- Output: Predicted weekly dose (mg) with 95% CI
- Validation: R² = 0.43 (vs. 0.17 for clinical-only algorithm)

**Tacrolimus Dose Calculator**
- Algorithm: CYP3A5-guided dosing
- CYP3A5 expressers (*1/*1, *1/*3): Standard dose → higher trough levels expected
- CYP3A5 non-expressers (*3/*3): May need 1.5-2x dose for target trough

**5-FU Dose Adjustment Calculator**
- Algorithm: DPYD activity score-based
- Input: DPYD diplotype → activity score
- Output: Percent dose reduction recommendation per EMA/CPIC guidelines

### 14.2 Alert Classification Engine

Three-tier alert system for clinical decision support:

| Alert Level | Criteria | Clinical Action | Example |
|------------|---------|----------------|---------|
| **CRITICAL** | CPIC Level A, avoid/contraindicate | Immediate prescriber notification, hard stop recommended | CYP2D6 UM + codeine, HLA-B*57:01+ + abacavir |
| **WARNING** | CPIC Level A/B, dose adjustment needed | Prescriber review required before dispensing | CYP2C19 IM + clopidogrel, SLCO1B1 *5/*5 + simvastatin |
| **INFO** | CPIC Level B/C, informational | FYI for clinical record, no immediate action required | CYP2D6 IM + ondansetron (may need higher dose) |

### 14.3 Drug Alternative Recommendation Engine

When a drug is flagged as AVOID or CONTRAINDICATED, the engine recommends alternatives:

```python
DRUG_ALTERNATIVES = {
    "codeine": {
        "CYP2D6_PM": ["morphine (not CYP2D6 dependent)", "oxycodone (minimal CYP2D6)",
                       "acetaminophen (non-opioid)", "NSAIDs (if appropriate)"],
        "CYP2D6_UM": ["morphine (direct-acting, no activation needed)",
                       "fentanyl (CYP3A4 metabolized)", "hydromorphone"],
    },
    "clopidogrel": {
        "CYP2C19_PM": ["prasugrel (not CYP2C19 dependent)",
                        "ticagrelor (not CYP2C19 dependent)"],
        "CYP2C19_IM": ["prasugrel", "ticagrelor",
                        "clopidogrel at increased dose (75 mg → 150 mg, limited evidence)"],
    },
    "simvastatin": {
        "SLCO1B1_poor": ["pravastatin (not SLCO1B1 dependent)",
                          "rosuvastatin (minimal SLCO1B1 dependence)",
                          "fluvastatin (not SLCO1B1 dependent)"],
    },
    # ... 30+ drug alternative maps
}
```

---

## 15. Reporting and Interoperability

### 15.1 Output Formats

- **Clinical PGx Report (PDF):** Comprehensive multi-page report for medical record
- **PGx Passport (PDF/card):** Wallet-sized card with critical PGx information for patient to carry
- **HL7 FHIR PGx Report:** Structured data output compatible with EHR integration (DiagnosticReport + Observation resources)
- **CDS Hooks:** Integration with EHR clinical decision support via CDS Hooks specification (when available)
- **JSON/API:** Structured data for programmatic consumption by other agents

### 15.2 FHIR Integration

PGx results will be represented using HL7 FHIR Genomics specifications:
- `DiagnosticReport` for overall PGx panel results
- `Observation` (component type `69548-6` Genetic variant assessment) for each gene
- `MedicationRequest` extensions for PGx-informed prescribing
- `Task` resources for recommended clinical actions

---

## 16. Product Requirements Document

### 16.1 User Stories

**Epic 1: Genome-to-PGx Profile (Foundation)**

| ID | Story | Acceptance Criteria | Priority |
|----|-------|-------------------|----------|
| PGX-001 | As a clinician, I want to upload a patient's VCF file and receive a complete PGx profile so I can make informed prescribing decisions | VCF processed, star alleles called for 25+ genes, phenotypes assigned, profile stored in Milvus | P0 |
| PGX-002 | As a clinician, I want to see which of my patient's current medications have PGx implications so I can address high-risk interactions | Medication list cross-referenced against PGx profile, alerts generated with CPIC evidence levels | P0 |
| PGX-003 | As a clinician, I want HLA alleles extracted and screened for drug hypersensitivity risk so I can prevent life-threatening reactions | HLA typing from WGS, screening against 12+ drug-HLA associations, mandatory screening compliance | P0 |

**Epic 2: Clinical Decision Support**

| ID | Story | Acceptance Criteria | Priority |
|----|-------|-------------------|----------|
| PGX-004 | As a prescriber, I want a traffic-light alert when I query a specific drug so I can quickly assess safety | Green/Yellow/Red classification with clinical recommendation and alternative drugs | P0 |
| PGX-005 | As a pharmacist, I want phenoconversion detection when reviewing medication lists so I can identify hidden drug-drug-gene interactions | All CYP inhibitors/inducers in med list flagged, effective phenotype calculated, cascade effects modeled | P1 |
| PGX-006 | As a cardiologist, I want genotype-guided warfarin dosing using the IWPC algorithm so I can optimize anticoagulation faster | CYP2C9+VKORC1+CYP4F2 genotypes + demographics → predicted dose with CI | P1 |

**Epic 3: Natural Language PGx Consultation**

| ID | Story | Acceptance Criteria | Priority |
|----|-------|-------------------|----------|
| PGX-007 | As a clinician with no PGx training, I want to ask questions in plain English and get clear recommendations so I don't need to interpret star alleles myself | Natural language query → grounded RAG response with CPIC citations and evidence levels | P0 |
| PGX-008 | As a clinician, I want to ask "What should I use instead?" when a drug is flagged and get ranked alternatives with rationale | Alternative drug list ranked by PGx suitability with explanation of why each is preferred | P1 |
| PGX-009 | As an oncologist, I want to check chemotherapy safety (DPYD, TPMT, NUDT15, UGT1A1) before starting treatment so I can prevent fatal toxicity | All relevant chemotherapy PGx genes assessed, dose adjustment calculated, CRITICAL alerts for deficiency | P0 |

**Epic 4: Reporting and Communication**

| ID | Story | Acceptance Criteria | Priority |
|----|-------|-------------------|----------|
| PGX-010 | As a clinician, I want to generate a PDF PGx report for the medical record so results are documented and accessible to other providers | Multi-page clinical report with gene table, alerts, recommendations, evidence citations | P1 |
| PGX-011 | As a patient, I want a PGx Passport card I can carry so I can inform emergency providers about my drug sensitivities | Wallet-card format with critical AVOID drugs, HLA alerts, and metabolizer status for key genes | P1 |
| PGX-012 | As a healthcare system, I want FHIR-formatted PGx data so I can integrate results into our EHR | DiagnosticReport + Observation FHIR resources with proper LOINC/SNOMED coding | P2 |

**Epic 5: Population Health and Analytics**

| ID | Story | Acceptance Criteria | Priority |
|----|-------|-------------------|----------|
| PGX-013 | As a health equity researcher, I want population-specific PGx analytics so I can identify disparities in PGx-relevant prescribing | Allele frequency data by ethnicity, metabolizer phenotype distributions, differential risk analysis | P2 |
| PGX-014 | As a pharmacy director, I want aggregate PGx data across our patient population so I can prioritize formulary changes | Population-level metabolizer phenotype distributions for drugs on formulary | P2 |

**Epic 6: Multi-Gene and Complex Scenarios**

| ID | Story | Acceptance Criteria | Priority |
|----|-------|-------------------|----------|
| PGX-015 | As a psychiatrist prescribing multiple psychotropic medications, I want a comprehensive CYP profile (2D6+2C19+1A2+3A4) so I can avoid multi-drug interactions | All psychiatric-relevant CYP phenotypes assessed, drug-drug-gene interactions modeled, recommendations for the combination | P1 |
| PGX-016 | As an internist managing a polypharmacy patient, I want the system to model phenoconversion cascades when adding a new drug so I can anticipate downstream effects | Interactive "what if" modeling: add Drug X → see phenotype changes → see downstream drug level predictions | P1 |
| PGX-017 | As a transplant physician, I want CYP3A5-guided tacrolimus dosing so I can reach therapeutic levels faster and reduce rejection risk | CYP3A5 phenotype → starting dose recommendation with expected trough level range | P1 |
| PGX-018 | As a pain specialist, I want comprehensive opioid safety assessment (CYP2D6 genetic + phenoconversion) so I can choose the safest analgesic | CYP2D6 genetic phenotype + effective phenotype (after phenoconversion) → opioid recommendation matrix | P0 |

### 16.2 Non-Functional Requirements

| ID | Requirement | Target |
|----|------------|--------|
| NFR-001 | VCF processing time (PGx extraction + star allele calling) | <60 seconds |
| NFR-002 | Single drug query response time | <5 seconds |
| NFR-003 | Full medication list review (20 drugs) | <15 seconds |
| NFR-004 | PGx report generation (PDF) | <10 seconds |
| NFR-005 | Concurrent user support | 10 simultaneous |
| NFR-006 | PGx profile storage capacity | 10,000+ patients |
| NFR-007 | Guideline update propagation | <24 hours after CPIC publication |
| NFR-008 | Patient data remains on-device | 100% (HIPAA/GDPR compliance) |

---

## 17. Data Acquisition Strategy

### 17.1 Primary Data Sources

| Source | Data Type | Access | Cost |
|--------|----------|--------|------|
| CPIC Guidelines | Gene-drug guidelines, star allele tables | Open access (cpicpgx.org) | Free |
| PharmGKB | Drug-gene annotations, clinical annotations | Academic license | Free |
| PharmVar | Star allele definitions, haplotype tables | Open access (pharmvar.org) | Free |
| DPWG | Dutch PGx therapeutic recommendations | Open access | Free |
| FDA Table of PGx Biomarkers | Drug label PGx information | Open access (FDA.gov) | Free |
| ClinVar | PGx variant classifications | Open access (NCBI) | Free |
| gnomAD | Population allele frequencies | Open access | Free |
| 1000 Genomes | Population reference genotypes | Open access | Free |
| ClinicalTrials.gov | PGx clinical trials | Open access | Free |
| PubMed/PMC | PGx literature | Open access | Free |

All primary data sources for the Pharmacogenomics Agent are **free and open access**, making this one of the most data-accessible agents in the HCLS AI Factory.

### 17.2 Data Ingestion Pipeline

1. **CPIC Guidelines:** Parse structured guideline PDFs and supplementary data tables → embed in `pgx_drug_guidelines`
2. **PharmVar:** Download allele definition tables → populate `pgx_gene_reference`
3. **PharmGKB:** API access for clinical annotations → populate `pgx_drug_interactions`
4. **FDA Labels:** Parse pharmacogenomic labeling sections → populate `pgx_fda_labels`
5. **gnomAD:** Extract PGx-relevant allele frequencies by population → populate `pgx_population_data`
6. **PubMed:** Systematic search for PGx implementation studies → populate `pgx_clinical_evidence`

---

## 18. Validation and Testing Strategy

### 18.1 Validation Tiers

**Tier 1: Variant Concordance (Unit)**
- Gold standard: GeT-RM (Genetic Testing Reference Materials) samples with known PGx genotypes
- Test: VCF → star allele caller → compare to GeT-RM consensus genotype
- Target: >99% concordance for simple genes, >95% for CYP2D6

**Tier 2: Phenotype Accuracy (Integration)**
- Gold standard: CPIC allele functionality tables
- Test: All possible diplotype combinations → phenotype → compare to CPIC assignment
- Target: 100% concordance with CPIC standardized terms

**Tier 3: Recommendation Accuracy (Clinical)**
- Gold standard: CPIC guideline recommendations + clinical pharmacist review
- Test: Simulated patient cases with known genotypes + medication lists → recommendations
- Target: 100% concordance with CPIC Level A guidelines for CRITICAL alerts

**Tier 4: Clinical Utility (Outcomes)**
- Post-implementation: Track ADR rates, time to therapeutic dose, prescriber adoption
- Compare PGx-guided vs. standard prescribing outcomes
- Measure cost avoidance from prevented ADRs

### 18.2 Test Case Library

50+ synthetic test cases covering:
- All 25+ pharmacogenes with multiple diplotype combinations
- All CPIC Level A drug-gene pairs
- Multi-gene scenarios (warfarin with CYP2C9 + VKORC1 + CYP4F2)
- Phenoconversion scenarios
- HLA screening edge cases
- Population-specific allele frequency considerations
- Polypharmacy with ≥10 medications

---

## 19. Regulatory Considerations

### 19.1 FDA Classification

The Pharmacogenomics Intelligence Agent operates as a **clinical decision support (CDS) tool** that:
- Presents PGx information from validated sources (CPIC, DPWG, FDA labels)
- Does NOT make autonomous prescribing decisions
- Requires clinician review and judgment for all recommendations
- Falls under FDA 21st Century Cures Act CDS exemptions (criteria i-iv) when:
  - Not intended to replace clinical judgment
  - Displays underlying evidence for clinician review
  - Does not generate alarms that trigger automatic action

### 19.2 Compliance Framework

- **HIPAA:** Patient genomic data stays on local DGX Spark. No PHI transmitted to cloud.
- **CLIA:** Star allele calling from research-grade WGS is informational. Clinical-grade results require CLIA-certified laboratory confirmation.
- **GDPR:** On-device processing eliminates cross-border data transfer concerns.
- **State PGx laws:** Some states require genetic counseling before PGx testing. The agent provides educational content but does not replace genetic counseling.

---

## 20. DGX Compute Progression

### 20.1 DGX Spark (Current Target)

- **GPU:** NVIDIA GPU with 128 GB unified memory
- **CPU:** NVIDIA Grace (72 ARM Neoverse cores)
- **RAM:** 128 GB unified (CPU+GPU shared)
- **Storage:** 4 TB NVMe SSD
- **Price:** $4,699
- **PGx Agent footprint:** ~15 GB (collections + models + reference data)
- **Concurrent capacity:** 10 simultaneous users

### 20.2 Scaling Path

As PGx adoption grows (more patients, more queries):
- **DGX Spark cluster (2-4 nodes):** Support 50+ concurrent users, faster batch VCF processing
- **DGX Station / Workstation:** Support health system deployment with 100+ users
- **DGX SuperPOD:** Population-scale PGx analytics across millions of patients

---

## 21. Implementation Roadmap

### Phase 1: Core PGx Engine (Weeks 1-6)

| Week | Deliverable |
|------|------------|
| 1-2 | VCF-to-PGx variant extraction pipeline; star allele caller for 15 simple genes |
| 3-4 | CYP2D6 complex allele handling; diplotype-to-phenotype translator |
| 5-6 | Milvus collections setup; CPIC guideline ingestion; basic drug-gene matching |

### Phase 2: Clinical Workflows (Weeks 7-12)

| Week | Deliverable |
|------|------------|
| 7-8 | Pre-emptive panel workflow; opioid safety workflow; HLA screening |
| 9-10 | Warfarin dosing calculator; antidepressant selection workflow |
| 11-12 | Chemotherapy safety (DPYD, TPMT, NUDT15); statin myopathy workflow |

### Phase 3: Advanced Features (Weeks 13-18)

| Week | Deliverable |
|------|------------|
| 13-14 | Phenoconversion detection engine; polypharmacy drug-drug-gene interaction analysis |
| 15-16 | Streamlit UI (10 tabs); PDF/FHIR report generation; PGx Passport |
| 17-18 | Population analytics; validation suite (50+ test cases); performance optimization |

---

## 22. Risk Analysis

### 22.1 Critical Risks

| Risk | Severity | Mitigation |
|------|---------|------------|
| Incorrect star allele calling leads to wrong phenotype and wrong drug recommendation | **Critical** -- patient harm | Multi-tier validation against GeT-RM; CLIA disclaimer for research-grade WGS; mandatory clinician review |
| CYP2D6 structural variants missed from short-read WGS | **High** -- incomplete profile | Clear documentation of limitations; flag when structural variant data unavailable; recommend confirmatory testing |
| Guideline updates not propagated, stale recommendations served | **High** -- outdated guidance | Automated CPIC RSS monitoring; version stamping on all recommendations; manual review process for updates |
| Clinician over-reliance on PGx recommendations without clinical context | **High** -- inappropriate automation | CDS design emphasizes "decision SUPPORT" not "decision MAKING"; all outputs include "clinician review required" |
| Phenoconversion modeling errors in complex polypharmacy | **Medium** -- missed interactions | Conservative alerting (flag potential phenoconversion even with weak inhibitors); clear uncertainty language |

### 22.2 Ethical Considerations

- **Health equity:** PGx guidelines are predominantly derived from European-ancestry populations. Allele frequency databases for African, Asian, and Latino populations are less complete. The agent must clearly communicate when evidence is population-limited.
- **Access equity:** The $4,699 DGX Spark + open-source model makes PGx CDS accessible to under-resourced clinics, but genomic sequencing itself remains costly.
- **Genetic determinism:** PGx phenotype is one factor in drug response. Environment, adherence, comorbidities, age, and other medications all contribute. The agent must contextualize genetic findings appropriately.
- **Incidental findings:** WGS-based PGx extraction may reveal disease-risk variants (e.g., BRCA1/2) incidentally. The agent must have a clear policy for managing incidental findings.

---

## 23. Competitive Landscape

### 23.1 Detailed Competitive Analysis

| Feature | Our Agent | OneOme RightMed | Myriad GeneSight | Invitae PGx | CPIC Guidelines |
|---------|-----------|----------------|-----------------|------------|----------------|
| Genes covered | 25+ | 24 | 12 (psych only) | 14 | 27 (guidelines) |
| Drug-gene pairs | 400+ | ~300 | ~60 | ~150 | ~80 |
| WGS integration | Yes (VCF pipeline) | No (panel only) | No | No | N/A |
| Multi-gene modeling | Yes | Limited | No | No | Manual |
| Phenoconversion | Yes | No | No | No | Mentioned |
| HLA screening | 12+ alleles | 3 | None | 2 | 6 |
| Natural language query | Yes (RAG + Claude) | No | No | No | No |
| On-device / HIPAA | Yes (DGX Spark) | Cloud-based | Cloud-based | Cloud-based | N/A |
| Cost | $4,699 (one-time HW) | Per-test subscription | Per-test ($400+) | Per-test | Free (text) |
| Open source | Yes (Apache 2.0) | No | No | No | Yes (guidelines) |

### 23.2 Unique Differentiators

1. **Genome-first:** Only solution that starts from WGS/WES VCF data rather than a targeted panel
2. **Multi-collection RAG:** Searches across 14 specialized knowledge bases simultaneously
3. **Phenoconversion modeling:** Detects when concurrent medications alter the patient's effective metabolizer status -- no competitor does this
4. **Natural language interface:** Clinicians ask questions in English, not star allele codes
5. **Local deployment:** Patient genomic data never leaves the clinic
6. **Zero per-test cost:** After initial hardware investment, unlimited PGx consultations

---

## 24. Discussion

### 24.1 Transformative Impact

The Pharmacogenomics Intelligence Agent addresses one of the most actionable gaps in modern medicine: the translation of genomic data into safe prescribing decisions. Unlike many precision medicine applications that remain aspirational, pharmacogenomics has Level A evidence from CPIC for dozens of drug-gene pairs, FDA boxed warnings for genetic testing, and demonstrated cost-effectiveness across multiple healthcare systems. The barriers to adoption are not scientific -- they are practical: complexity, integration, and clinical workflow friction.

By embedding PGx intelligence into a natural language interface backed by multi-collection RAG, the agent eliminates the need for clinicians to understand star allele nomenclature, activity scores, or guideline interpretation. A physician can simply ask: "Is this drug safe for my patient?" and receive a grounded, evidence-cited answer in seconds.

### 24.2 Integration with Existing Agents

The Pharmacogenomics Agent completes a critical circle in the HCLS AI Factory:
- The **Genomics Pipeline** generates the variants
- The **Biomarker Agent** provides initial PGx screening
- The **Pharmacogenomics Agent** provides deep clinical PGx consultation
- The **Oncology Agent** uses PGx for chemotherapy safety
- The **Autoimmune Agent** uses PGx for immunosuppressant dosing
- The **Cardiology Agent** uses PGx for anticoagulant and statin optimization

### 24.3 The Vision: PGx as Standard of Care

The ultimate goal is a future where every patient's genome is processed once and their PGx profile is available for every prescribing decision for life. The Pharmacogenomics Intelligence Agent makes this vision technically feasible on a $4,699 device -- democratizing access to pharmacogenomic intelligence that currently requires multi-million-dollar institutional programs.

With 106,000 ADR deaths per year in the U.S. alone, and 95-99% of patients carrying actionable PGx variants, the potential impact is measured not in efficiency gains but in lives saved.

---

## 25. Conclusion

This paper has presented the comprehensive architecture, clinical rationale, and product requirements for the Pharmacogenomics Intelligence Agent -- a multi-collection RAG system that transforms raw genomic data into actionable, evidence-grounded prescribing guidance for over 400 drug-gene interactions across 25+ pharmacogenes. The agent addresses the most preventable cause of drug-related death by bridging the gap between genomic knowledge and clinical practice through natural language clinical decision support.

Key architectural innovations include:
- **VCF-to-prescribing pipeline:** Direct integration with the HCLS AI Factory genomics pipeline
- **14 specialized Milvus collections** covering CPIC/DPWG/FDA guidelines, PharmGKB annotations, HLA hypersensitivity screening, phenoconversion modeling, dosing algorithms, and population-specific data
- **Multi-gene interaction engine** for complex scenarios (warfarin, psychiatric polypharmacy, transplant immunosuppression)
- **Phenoconversion detector** that identifies when concurrent medications alter effective metabolizer status -- a capability no competing product offers
- **8 clinical workflows** covering the highest-impact prescribing scenarios where PGx testing has CPIC Level A evidence
- **Privacy-first architecture** with all patient genomic data remaining on the local DGX Spark device

The agent represents the seventh intelligence module in the HCLS AI Factory platform, bringing the total agent portfolio to: Precision Biomarker, Precision Oncology, CAR-T Intelligence, Imaging Intelligence, Precision Autoimmune, Cardiology Intelligence, and Pharmacogenomics Intelligence -- a comprehensive precision medicine platform that takes patient care from genome to safe prescription on a single $4,699 device.

---

## 26. References

1. Lazarou J, Pomeranz BH, Corey PN. Incidence of adverse drug reactions in hospitalized patients: a meta-analysis of prospective studies. JAMA. 1998;279(15):1200-1205. PMID: 9555760
2. Relling MV, Klein TE. CPIC: Clinical Pharmacogenetics Implementation Consortium of the Pharmacogenomics Research Network. Clin Pharmacol Ther. 2011;89(3):464-467. PMID: 21270786
3. Caudle KE, et al. Standardizing CYP2D6 Genotype to Phenotype Translation: Consensus Recommendations from the Clinical Pharmacogenetics Implementation Consortium and Dutch Pharmacogenetics Working Group. Clin Transl Sci. 2020;13(1):116-124. PMID: 31647186
4. Crews KR, et al. Clinical Pharmacogenetics Implementation Consortium Guidelines for Cytochrome P450 2D6 Genotype and Codeine Therapy: 2014 Update. Clin Pharmacol Ther. 2014;95(4):376-382. PMID: 24458010
5. Scott SA, et al. Clinical Pharmacogenetics Implementation Consortium Guidelines for CYP2C19 Genotype and Clopidogrel Therapy: 2013 Update. Clin Pharmacol Ther. 2013;94(3):317-323. PMID: 23698643
6. Johnson JA, et al. Clinical Pharmacogenetics Implementation Consortium Guidelines for CYP2C9 and VKORC1 Genotypes and Warfarin Dosing: 2017 Update. Clin Pharmacol Ther. 2017;102(3):397-404. PMID: 28198005
7. SEARCH Collaborative Group. SLCO1B1 variants and statin-induced myopathy -- a genomewide study. N Engl J Med. 2008;359(8):789-799. PMID: 18650507
8. Amstutz U, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for Dihydropyrimidine Dehydrogenase Genotype and Fluoropyrimidine Dosing: 2017 Update. Clin Pharmacol Ther. 2018;103(2):210-216. PMID: 29152729
9. Mallal S, et al. HLA-B*5701 screening for hypersensitivity to abacavir. N Engl J Med. 2008;358(6):568-579. PMID: 18256392
10. Hershfield MS, et al. Clinical Pharmacogenetics Implementation Consortium Guidelines for Human Leukocyte Antigen-B Genotype and Allopurinol Dosing. Clin Pharmacol Ther. 2013;93(2):153-158. PMID: 23232549
11. Leckband SG, et al. Clinical Pharmacogenetics Implementation Consortium Guidelines for HLA-B Genotype and Carbamazepine Dosing. Clin Pharmacol Ther. 2013;94(3):324-328. PMID: 23695185
12. International Warfarin Pharmacogenetics Consortium. Estimation of the warfarin dose with clinical and pharmacogenetic data. N Engl J Med. 2009;360(8):753-764. PMID: 19228618
13. Birdwell KA, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guidelines for CYP3A5 Genotype and Tacrolimus Dosing. Clin Pharmacol Ther. 2015;98(1):19-24. PMID: 25801146
14. Relling MV, et al. Clinical Pharmacogenetics Implementation Consortium Guideline for Thiopurine Dosing Based on TPMT and NUDT15 Genotypes: 2018 Update. Clin Pharmacol Ther. 2019;105(5):1095-1105. PMID: 30447069
15. Hicks JK, et al. Clinical Pharmacogenetics Implementation Consortium Guideline (CPIC) for CYP2D6 and CYP2C19 Genotypes and Dosing of Tricyclic Antidepressants: 2016 Update. Clin Pharmacol Ther. 2017;102(1):37-44. PMID: 27997040
16. Hicks JK, et al. Clinical Pharmacogenetics Implementation Consortium (CPIC) Guideline for CYP2D6 and CYP2C19 Genotypes and Dosing of Selective Serotonin Reuptake Inhibitors. Clin Pharmacol Ther. 2015;98(2):127-134. PMID: 25974703
17. PharmGKB. https://www.pharmgkb.org/
18. PharmVar. https://www.pharmvar.org/
19. CPIC. https://cpicpgx.org/
20. Dunnenberger HM, et al. Preemptive clinical pharmacogenetics implementation: current programs in five US medical centers. Annu Rev Pharmacol Toxicol. 2015;55:89-106. PMID: 25292429
21. Weitzel KW, et al. The IGNITE network: a model for genomic medicine implementation and research. BMC Med Genomics. 2016;9:1. PMID: 26729011
22. Swen JJ, et al. Pharmacogenetics: from bench to byte -- an update of guidelines. Clin Pharmacol Ther. 2011;89(5):662-673. PMID: 21412232
23. FDA Table of Pharmacogenomic Biomarkers in Drug Labeling. https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling
24. Bousman CA, et al. Pharmacogenetic tests and depressive symptom remission: a meta-analysis of randomized controlled trials. Pharmacogenomics. 2019;20(1):37-47. PMID: 30520364
25. Empey PE, et al. Multisite investigation of outcomes with implementation of CYP2C19 genotype-guided antiplatelet therapy after percutaneous coronary intervention. JACC Cardiovasc Interv. 2018;11(2):181-191. PMID: 29348010
26. Samwald M, et al. Incidence of exposure of patients in the United States to multiple drugs for which pharmacogenomic guidelines are available. PLoS One. 2016;11(10):e0164972. PMID: 27736993
