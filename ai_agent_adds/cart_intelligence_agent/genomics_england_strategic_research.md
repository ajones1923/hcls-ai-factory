# Genomics England -- Strategic Research Brief

**Prepared for:** VAST Data Sales/Partnerships Teams
**Date:** March 3, 2026
**Classification:** Internal -- Strategic Planning

---

## 1. Organization Overview

**Genomics England** is a company wholly owned by the **UK Department of Health and Social Care (DHSC)**. It was established on **5 July 2013** by then-Secretary of State for Health Jeremy Hunt, with the launch timed to coincide with the 65th birthday of the NHS.

- **Legal Status:** Government-owned limited company (Company No. 08493132)
- **Headquarters:** London, UK
- **Headcount:** ~250 employees (2024 Gender Pay Gap Report), growing ~2% YoY
- **Mission:** To enable genomic healthcare delivery and genomic research at scale, working with the NHS to help doctors diagnose, treat, and prevent illnesses -- particularly rare diseases and cancer
- **Dual Role:** Operates both as a **clinical service provider** for the NHS and as a **research data platform** enabling linked genomic and clinical data to be used for research

**Key fact:** The UK became the **first nation in the world** to apply whole genome sequencing at scale in direct healthcare, driven largely by Genomics England's infrastructure.

---

## 2. The 100,000 Genomes Project

The foundational programme that established Genomics England's global position.

| Metric | Detail |
|--------|--------|
| **Goal** | Sequence 100,000 whole genomes from NHS patients |
| **Completion** | 100,000th genome sequenced in **December 2018** (5 years) |
| **Participants** | ~85,000-97,000 patients and family members |
| **Genomes sequenced** | 106,292 uniquely sequenced samples from 88,505 participants (Data Release v19, Oct 2024) |
| **Disease split** | ~70,000 rare disease genomes, ~30,000 cancer genomes |
| **Data generated** | **21-50 petabytes** (figures vary by reporting period and what is included) |
| **Clinical data items** | Over **2 billion** linked data items (hospital events, deaths, cancer diagnostics/pathology/therapeutics) |
| **Cancer tumors profiled** | 13,880 solid tumors across 33 cancer types |
| **Recruiting hospitals** | 84 hospitals across England |

### Clinical Impact
- A landmark study in the **New England Journal of Medicine** demonstrated that WGS led to a **new diagnosis for 25% of rare disease participants**
- **14%** of diagnoses involved variants in regions of the genome that would be missed by other testing methods
- Over **2,000 researchers** have been enabled to analyze data from the project
- Created the **largest whole genome dataset linked to lifelong health records** anywhere in the world

---

## 3. Current Programs (2024-2026)

### 3.1 Generation Study (Newborn Genomes Programme)
**Funding:** GBP 105 million

The NHS's pioneering study to evaluate whole genome sequencing for newborn screening.

- **Target:** Sequence genomes of **100,000 newborn babies**
- **Conditions screened:** Over **200 rare diseases** caused by variants in over **500 genes**
- **Recruitment started:** March 2024
- **Recruitment end (expected):** December 2026; results returned through March 2027
- **Current enrollment:** 25,000+ babies enrolled (milestone announced); 36,056 participants from 72+ hospitals as of January 2026
- **Early finding:** WGS identifies rare, treatable conditions in approximately **1 in every 200 babies**
- **National expansion:** In June 2025, NHS announced a plan to offer WGS to **every newborn within the next decade**, with national rollout beginning in 2026 -- a **GBP 650 million initiative**
- **Scale implication:** If fully rolled out, this would cover **~600,000 babies per year** born in England and Wales

### 3.2 Cancer 2.0 Programme
**Funding:** GBP 26 million (initial)

Two interlocking sub-programmes targeting the **300,000+ people diagnosed with cancer per year** in England:

**Long-Read Sequencing Programme:**
- Sequences long strands of DNA without fragmenting, revealing structural features invisible to short-read sequencing
- Proof-of-concept for building analytical pipeline alongside existing short-read sequencing
- Long-read formats contain **~5x more data** than short-read formats

**Multi-Modal Programme:**
- Combines genomic, pathology imaging (whole slide images), radiology (MRI/CT), and clinical follow-up data
- Partnership with **National Pathology Imaging Co-operative (NPIC)** to digitize all pathology slides from 100KGP cancer participants
- Onboarding **250,000 curated digital pathology and radiology images**
- Machine learning used to identify cancer features driving prognosis or treatment response
- **15,000+ cancer participants** from 84 hospitals across 20 solid cancer subspecialties

### 3.3 Diverse Data Initiative
**Funding:** GBP 22 million

- **Target:** Sequence genomes of up to **25,000 non-European ancestry participants**
- **Problem addressed:** 80% of genome-wide association studies involve people of European ancestry; polygenic risk scores are ~4x more accurate for European ancestry populations
- **Focus areas:** Sickle cell genomic dataset, maternal health/preterm birth inequalities, equity in genomic medicine, emerging technologies
- **Approach:** Community engagement-driven, designed with underrepresented communities

### 3.4 Pharmacogenomics (PROGRESS Study)
- **20 GP practices** across England participating in the PROGRESS trial
- **9 practices** using PROGRESSRx (integrates pharmacogenomic results directly into prescribing systems with guidance pop-ups at point of prescribing)
- CYP2C19 testing pilot launched at 4 sites
- NHS 10 Year Plan commits to pharmacogenomics as routine care

### 3.5 Adult Population Genomics Programme
- **Pilot study** of up to **150,000 participants** recruited from the general population
- "Genomics, Healthcare and You" engagement programme running through 2025/26
- Assessing delivery approaches and health benefits from population-level genomic screening

### 3.6 NHS Genomic Medicine Service (GMS) Partnership
- Genomics England provides analytical and IT systems for the NHS GMS
- **7 supra-regional Genomic Laboratory Hubs (GLHs)** deliver testing nationally
- **National Genomic Research Library** -- partnership between NHS England and Genomics England
- **Liquid biopsy** testing available for all eligible lung and breast cancer patients in NHS hospitals across England (2025) -- first national "liquid biopsy first" approach
- **GeNotes** educational platform: 12 specialties, ~600 resources covering 150+ rare diseases
- **28 rare disease collaborative networks** supporting treatment pathway development
- Digitised National Genomic Test Directory being introduced for clinical use nationally during 2026

### 3.7 COVID-19 Response
- Sequenced entire genomes of thousands of severely affected patients
- Compared with mild-symptom patients to investigate how genetic variation influenced susceptibility/severity
- Complemented COG-UK consortium (1.8M SARS-CoV-2 viral genomes sequenced by Dec 2021)

---

## 4. Data Infrastructure & Technology Stack

### 4.1 Hybrid Architecture: On-Premises + AWS Cloud

Genomics England operates a **dual-tier architecture** spanning on-premises high-performance storage and AWS cloud compute.

#### On-Premises Storage (WEKA)
| Component | Specification |
|-----------|--------------|
| **Primary tier** | 1.3 PB NVMe-based flash storage |
| **Secondary tier** | 40 PB object storage |
| **Total namespace** | 41 PB unified single namespace with automated data tiering |
| **Throughput** | 135+ GB/s from NVMe tier |
| **Performance gain** | 10x+ improvement over legacy NFS-based NAS |
| **Disaster recovery** | Object tier geo-distributed across **3 sites, each 50 miles apart** |
| **Data tiering** | Automated flash-to-disk tiering for optimal economics |
| **Legacy limit** | Previous storage system capped at 25 PB |

#### AWS Cloud Infrastructure
| Component | Detail |
|-----------|--------|
| **Research Environment** | Linux Virtual Desktop Interface (VDI) hosted by AWS |
| **HPC Cluster** | High Performance Compute cluster for genomic analyses |
| **Compute** | Amazon EC2 (on-demand capacity), exploring AWS Graviton Processors |
| **Storage** | Amazon S3 for genomic data storage |
| **ML/AI** | Amazon SageMaker for cancer sub-typing and survival prediction pipelines |
| **Bioinformatics** | AWS HealthOmics (Nextflow, WDL, CWL pipeline execution) |
| **Modernization** | Migrating to microservices: AWS Lambda, AWS Fargate (serverless) |
| **Cost optimization** | Saved GBP 1M+/year through FinOps with The Server Labs (67% from non-production instance scheduling) |
| **Performance** | 99% reduction in researcher task completion time (25 hours to 23 seconds for common tasks) |

#### Database Layer
- **MongoDB Enterprise Advanced** -- powers data science platform; reduced complex query times from hours to milliseconds; selected for end-to-end encryption, fine-grained access control

#### Platform Partner
- **Lifebit** -- Selected as AWS Partner to develop the Trusted Research Environment platform; standardizing data to OMOP common data model for interoperability

### 4.2 Trusted Research Environment (TRE)

The secure research access platform enabling external researchers to work with genomic data:

- **Access model:** Virtual Desktop Interface -- researchers access only de-identified data
- **Oversight:** Access Review Committee approves all researcher access
- **Data export:** Only summary data exported through "Airlock" after manual review
- **Researchers enabled:** 2,000+ approved researchers
- **Data available:** WGS data enriched with clinical/family information + linked NHS health records
- **Tools:** Portfolio for both coding and non-coding researchers (Python, R, command-line bioinformatics)
- **Researcher cost model:** Small budget for compute costs + computer access (pay-as-you-go)

### 4.3 Data Scale

| Metric | Value |
|--------|-------|
| **Total data (100KGP)** | 21-50 PB |
| **Genomes sequenced** | 100,000+ (100KGP) expanding toward 5 million |
| **Clinical data items** | 2 billion+ |
| **Projected growth** | From ~100K to **5 million genomes** (UK-wide ambition) |
| **Newborn programme potential** | ~600,000 genomes/year if fully rolled out |
| **Long-read data impact** | ~5x more data per genome than short-read |
| **Imaging data** | 250,000 curated pathology/radiology images being onboarded |

---

## 5. Partnerships & Collaborations

### 5.1 Illumina (Primary Sequencing Partner)
- **Contract:** Up to **300,000 whole genome equivalents** over 5 years, option to increase to **500,000**
- **Historical partnership:** GBP 524 million (~$524M) partnership for original 100KGP
- **Sequencing lab:** Illumina Laboratory Services, Cambridge, UK (ISO 15189 accredited)
- **Technology:** NovaSeq 6000, BaseSpace Sequence Hub, Illumina analysis tools

### 5.2 GENE Consortium (Pharma)
The **Genomics Expert Network for Enterprises** consortium brings major pharmaceutical companies to mine genomic data:
- **Members:** GSK, AstraZeneca, Roche, Biogen, AbbVie, Takeda, Alexion Pharmaceuticals, UCB, Helomics, Dimension Therapeutics
- **Entry fee:** GBP 250,000 per company
- **Purpose:** New diagnostic avenues and therapeutic targets from NHS patient data
- **Model:** Data accessed within the secure Research Environment; no data leaves the TRE

### 5.3 Technology Partners
| Partner | Role |
|---------|------|
| **AWS** | Primary cloud infrastructure provider |
| **WEKA** | On-premises high-performance storage |
| **Lifebit** | TRE platform development and data standardization |
| **The Server Labs** | AWS FinOps optimization, infrastructure modernization |
| **MongoDB** | Database platform for data science |
| **National Pathology Imaging Co-operative (NPIC)** | Cancer pathology slide digitization |

### 5.4 Academic & Research Partners
- Wellcome Sanger Institute
- European Bioinformatics Institute
- University of Edinburgh
- King's College London (Chief Scientific Officer's home institution)
- UK Biobank (partner in 5 million genomes ambition)
- NIHR Biomedical Research Centres

### 5.5 Funding Partners
- UK Research and Innovation (UKRI)
- LifeArc
- National Institute for Health and Care Research (NIHR)
- Wellcome Trust

---

## 6. Leadership Team

### Executive Leadership
| Role | Name | Background |
|------|------|------------|
| **Chair** | Baroness Nicola Blackwood | |
| **CEO** | Dr Richard Scott | Joined Genomics England in 2015; appointed CEO March 2024 |
| **Chief Scientific Officer** | Professor Matt Brown | Professor of Medicine, King's College London; Director, NIHR Guy's & St Thomas' BRC |
| **Chief Technology & Product Officer** | Julian Thomas | 20+ years leading technology teams; former NHS Digital and UK Health Security Agency |
| **Chief Medical Officer** | Dr Ellen Thomas | |
| **Chief Bioinformatician** | Dr Augusto Rendon | |
| **Chief Ethics & Engagement Officer** | Dr Natalie Banner | |
| **Chief Financial & Performance Officer** | Catherine Byers | |
| **Chief People Officer** | Jackie Kinsey | |
| **General Counsel & DPO** | Nick Maltby | |
| **Chief of Staff & Director of Strategy** | Chris Schonewald | |

### Board Non-Executive Directors
- Professor Ewan Birney
- Dr Annalisa Jenkins
- Dr Gail Marzetti
- Nicola Perrin MBE
- Dr Keith Stewart
- Dr Sergei Yakneen
- Dr Vikram Bajaj

### Previous CEOs
- Chris Wigley (2018-2024)
- Sir Mark Caulfield (Interim CEO, also served as Chief Scientist)

---

## 7. Funding Model

### Primary Funding: Government (DHSC)
The vast majority of Genomics England's funding comes from the Department of Health and Social Care.

### Major Government Investments
| Programme | Amount (GBP) |
|-----------|-------------|
| Generation Study (Newborn) | 105 million |
| Cancer 2.0 | 26 million (initial) |
| Diverse Data Initiative | 22 million |
| Cutting-edge genomics research | 22.4 million |
| **Total GBP 175M package** (Dec 2022) | 175 million |
| National newborn WGS rollout (announced June 2025) | 650 million |
| NHS Secure Data Environments | 13.5 million (GE share of GBP 200M NHS data programme) |

### Commercial Revenue
- Data access fees from industry researchers via the Research Network
- GENE Consortium pharma membership fees (GBP 250K per company)
- Revenue reinvested in the national healthcare system

### Research Grants
- UK Research and Innovation (UKRI)
- LifeArc
- National Institute for Health and Care Research (NIHR)
- Wellcome Trust
- Medical Research Council (MRC)

---

## 8. AI/ML Initiatives

### Amazon SageMaker Cancer Pipelines
- Automatic **cancer sub-typing and survival detection** pipeline
- Multi-modal ML combining gene expression, mutation data, copy number variants, and histopathology whole slide images
- Proof of concept on breast cancer (BRCA) and gastrointestinal cancer types (Pan-GI)
- Can predict genomic features directly from pathology images

### Multi-Modal Machine Learning
- Combining genomic, pathology imaging, radiology imaging, and clinical data
- Machine learning identifies cancer features driving prognosis or treatment response
- "Pathogenomic" approach improves classification vs. unimodal methods

### Clinical AI
- GBP 180 million Healthcare AI Solutions framework expected to open for bids summer 2025, going live early 2026
- Liquid biopsy for lung and breast cancer using AI-enhanced genomic testing
- 1,600+ lung cancer patients and ~600 advanced breast cancer patients tested since April 2025

### Infrastructure for AI
- AWS SageMaker for model training and deployment
- Fully managed ML infrastructure within the TRE
- Graviton processors being explored for efficiency
- Serverless compute (Lambda, Fargate) for scalable pipeline execution

---

## 9. Recent Developments (2024-2026)

| Date | Development |
|------|-------------|
| **March 2024** | Dr Rich Scott appointed CEO |
| **March 2024** | Generation Study recruitment begins |
| **Oct 2024** | Data Release v19 (106,292 genomes, 90,173 participants) |
| **Nov 2024** | 1 million genomes sequenced milestone |
| **Jan 2026** | Generation Study reaches 36,056 participants from 72+ hospitals |
| **2025** | Julian Thomas appointed Chief Technology and Product Officer |
| **2025** | National liquid biopsy first approach for lung and breast cancer |
| **2025** | PROGRESS pharmacogenomics trial expands to 20 GP practices |
| **June 2025** | GBP 650M national newborn WGS rollout plan announced |
| **2025/26** | "Genomics, Healthcare and You" adult genomics engagement programme |
| **2026** | Digitised National Genomic Test Directory going live nationally |
| **2026** | Cancer 2.0 long-read sequencing pipeline development |
| **2026** | Rare disease network of excellence launching with proteomics and metabolomics pathways |
| **2026** | Adult population genomics pilot study (up to 150,000 participants) |

---

## 10. Strategic Assessment: Why VAST AI OS Matters to Genomics England

### Pain Points Visible in Their Current Infrastructure

1. **Hybrid complexity:** Operating WEKA on-premises (41 PB) alongside AWS cloud creates management overhead, data movement costs, and consistency challenges. Their infrastructure is split between on-prem NVMe flash + object storage and AWS S3/EC2.

2. **Explosive data growth:** From 50 PB to a projected need for hundreds of petabytes as they scale from 100K to 5 million genomes. Long-read sequencing adds 5x more data per genome. The legacy system already hit its 25 PB ceiling once.

3. **Multi-modal data integration:** Cancer 2.0 is onboarding 250,000 pathology/radiology images alongside genomic and clinical data -- fundamentally different data types requiring unified access at speed.

4. **AI/ML pipeline performance:** SageMaker cancer pipelines, multi-modal ML, and future generative genomics all demand high-throughput, low-latency access to massive datasets. Their 135 GB/s throughput from WEKA NVMe will face pressure at 5M-genome scale.

5. **Cost pressure:** Already engaged FinOps optimization (saved GBP 1M/year). Government-funded organization under constant pressure to demonstrate value. Pay-as-you-go model for researchers means infrastructure cost directly affects accessibility.

6. **Secure multi-tenancy:** 2,000+ researchers, 10+ pharma companies, NHS clinical services, and internal teams all accessing the same data through the TRE with different permission levels. Airlock-based data export requires granular control.

7. **Geo-distributed DR:** Currently maintaining 3 sites 50 miles apart for disaster recovery. VAST's native geo-distribution could simplify this.

8. **Newborn screening scale:** If the GBP 650M national rollout proceeds, ~600,000 genomes/year will need to be processed, stored, and made available for both clinical and research use -- a step-change in throughput requirements.

### Where VAST AI OS Could Add Value

- **Unified data platform** replacing the WEKA on-prem + AWS S3 hybrid with a single namespace across flash/object/cloud tiers
- **Multi-protocol access** supporting the diverse workloads (genomic pipelines, ML training, imaging, researcher desktops)
- **Scalability to exabyte-class** as Genomics England grows toward 5 million genomes and national newborn screening
- **Performance at scale** maintaining or exceeding 135 GB/s as dataset grows 10-50x
- **Native S3 compatibility** for seamless integration with existing AWS-based SageMaker/HealthOmics pipelines
- **Cost efficiency** through intelligent tiering and data reduction, directly benefiting the government-funded cost model
- **Security and compliance** for NHS data sovereignty requirements and TRE multi-tenancy
- **AI-native architecture** aligned with their growing ML/AI ambitions (cancer prediction, multi-modal learning, pharmacogenomics)

### Competitive Context
- **WEKA** is the current on-premises storage incumbent (selected for 5M Genomes Project)
- **AWS** is the primary cloud partner (S3, EC2, SageMaker, HealthOmics)
- **Lifebit** is the platform partner for TRE
- **MongoDB** handles the database layer
- **The Server Labs** is the AWS FinOps/modernization partner
- **Illumina** provides the sequencing infrastructure

### Entry Points for Engagement
1. **CTPO Julian Thomas** -- new in role, 20+ years tech leadership including NHS Digital; likely evaluating infrastructure strategy
2. **Chief Bioinformatician Dr Augusto Rendon** -- owns pipeline performance and data architecture decisions
3. **Cancer 2.0 and Multi-Modal teams** -- facing the most acute data diversity and performance challenges
4. **Generation Study expansion planning** -- 600,000 genomes/year requires fundamental infrastructure rethinking
5. **Adult Population Genomics Pilot** -- 150,000 participants is a new programme that needs infrastructure planning from scratch

---

## Sources

- [Genomics England -- About Us](https://www.genomicsengland.co.uk/about-us)
- [Genomics England -- Origins](https://www.genomicsengland.co.uk/about-us/origins)
- [Genomics England -- Governance](https://www.genomicsengland.co.uk/about-us/governance)
- [100,000 Genomes Project](https://www.genomicsengland.co.uk/initiatives/100000-genomes-project)
- [Genomics England -- Wikipedia](https://en.wikipedia.org/wiki/Genomics_England)
- [Cancer 2.0 Programme](https://www.genomicsengland.co.uk/initiatives/cancer)
- [Newborn Genomes Programme](https://www.genomicsengland.co.uk/initiatives/newborns)
- [Diverse Data Initiative](https://www.genomicsengland.co.uk/initiatives/diverse-data)
- [NHS Genomic Medicine Service](https://www.genomicsengland.co.uk/genomic-medicine/nhs-gms)
- [Genomics England Research Environment](https://www.genomicsengland.co.uk/research/research-environment)
- [Research Environment Technical Docs](https://re-docs.genomicsengland.co.uk/)
- [Data Release v19](https://re-docs.genomicsengland.co.uk/release19/)
- [AWS Case Study -- Genomics England](https://aws.amazon.com/solutions/case-studies/genomics-england/)
- [AWS + The Server Labs Partnership](https://aws.amazon.com/partners/success/genomics-england-the-server-labs/)
- [AWS SageMaker Cancer Prediction](https://aws.amazon.com/blogs/machine-learning/genomics-england-uses-amazon-sagemaker-to-predict-cancer-subtypes-and-patient-survival-from-multi-modal-data/)
- [AWS Multimodal Whole Slide Images](https://aws.amazon.com/blogs/publicsector/aws-helps-genomics-englands-multimodal-programme-accelerate-research-with-whole-slide-images/)
- [WEKA -- Genomics England Case Study](https://www.weka.io/customers/genomics-england/)
- [WEKA Selection for 5M Genomes Project](https://www.weka.io/company/weka-newsroom/press-releases/genomics-england-selects-wekaio-to-accelerate-genomics-research-for-5-million-genomes-project/)
- [MongoDB -- Genomics England](https://www.mongodb.com/company/newsroom/press-releases/genomics-england-uses-mongodb-to-power-the-data-science-behind-the-100000-genomes-project)
- [Lifebit -- Genomics England Data Standardisation](https://lifebit.ai/blog/lifebit-genomics-england-data-standardisation/)
- [Lifebit -- Genomics England 2025 Revolution](https://lifebit.ai/blog/genomics-england/)
- [Illumina -- GE Whole Genome Partnership](https://www.genomicsengland.co.uk/news/genomics-england-and-illumina-sequence-whole-genomes-for-nhs-gms)
- [Illumina + GE $524M Partnership](https://www.genengnews.com/topics/translational-medicine/illumina-genomics-england-launch-524m-partnership-for-100k-genomes-project/)
- [GENE Consortium](https://www.fiercebiotech.com/partnering/big-pharma-joins-u-k-s-genomics-project-eyes-on-drug-discovery)
- [Industry Research Partnerships](https://www.genomicsengland.co.uk/research/partnerships)
- [GBP 175M Funding Announcement](https://www.genomicsengland.co.uk/news/funding-genomics-research)
- [GBP 175M Government Announcement](https://www.gov.uk/government/news/over-175-million-for-cutting-edge-genomics-research)
- [GBP 650M Newborn WGS Rollout](https://www.theanalyticalscientist.com/issues/2025/articles/july/uk-plans-nationwide-dna-sequencing-of-newborns/)
- [Generation Study 25,000 Milestone](https://www.genomicsengland.co.uk/news/25-000-babies-join-groundbreaking-generation-study)
- [Generation Study Protocol v8.0](https://www.genomicsengland.co.uk/assets/documents/Protocol-V8.0-02-December-2025-For-External-Use-1.pdf)
- [Dr Rich Scott -- CEO Appointment](https://www.genomicsengland.co.uk/news/genomics-england-appoints-dr-rich-scott-as-chief-executive-officer)
- [Julian Thomas -- CTPO Appointment](https://www.genomicsengland.co.uk/news/genomics-england-welcomes-new-chief-technology-and-product-officer)
- [Professor Matt Brown -- CSO Appointment](https://www.genomicsengland.co.uk/news/genomics-england-appoints-chief-scientific-officer)
- [NHS GMS Achievements 2024](https://www.england.nhs.uk/long-read/the-nhs-genomic-medicine-service-achievements-in-2024/)
- [NHS GMS Structure](https://www.england.nhs.uk/genomics/the-structure-of-the-nhs-genomic-medicine-service/)
- [Genome UK 2022-2025 Implementation Plan](https://www.gov.uk/government/publications/genome-uk-2022-to-2025-implementation-plan-for-england/genome-uk-2022-to-2025-implementation-plan-for-england)
- [England Rare Diseases Action Plan 2026](https://www.gov.uk/government/publications/england-rare-diseases-action-plan-2026/england-rare-diseases-action-plan-2026-main-report)
- [5 Million Genomes Ambition](https://www.genomicsengland.co.uk/news/matt-hancock-announces-sequencing-ambition)
- [PROGRESS Pharmacogenomics Trial](https://centralsouthgenomics.nhs.uk/2025/02/05/national-rollout-of-phase-ii-progress-pharmacogenomics-project/)
- [Pharmacogenomics NHS Implementation](https://bpspubs.onlinelibrary.wiley.com/doi/10.1002/bcp.70109)
- [Genomics, Healthcare and You](https://www.genomicsengland.co.uk/news/genomics-healthcare-and-you-engagement-programme-launched-to-shape-future-adult-genomics-study)
- [Adult Population Genomics -- Festival of Genomics 2026](https://festivalofgenomics.com/london/2026-agenda/sponsored-slot-available-abzu)
- [COVID-19 Study](https://www.genomicsengland.co.uk/initiatives/covid-19)
- [Cancer Long-Read Sequencing](https://www.genomicsengland.co.uk/news/cancer-long-read-sequencing)
- [Nature Medicine -- 13,880 Tumors Study](https://www.nature.com/articles/s41591-023-02682-0)
- [AI in NHS Genomics Report (Oct 2024)](https://genomicainetwork.nhs.uk/wp-content/uploads/2025/01/AI-in-genomics-report-final-251124.pdf)
- [Three Trends Reshaping UK Precision Diagnostics 2026](https://med4nexus.com/insights/three-transformative-trends-reshaping-uk-precision-diagnostics-in-2026/)
- [Gender Pay Gap Report 2024](https://www.genomicsengland.co.uk/assets/documents/People/Genomics_England_Gender-pay-gap_2024_v4.pdf)
