---
search:
  exclude: true
tags:
  - Launch
  - Social Media
  - Announcement
---

# Launch Announcement Drafts

> Ready-to-post drafts for LinkedIn, Twitter/X, and Hacker News.
> Edit to your voice before posting.

---

## LinkedIn Post

---

**I've been working on this for 14 years. Today I'm giving it away.**

In 2012, I set out to use my high-performance computing background for something that mattered. I started with one conviction: no parent should ever have to lose a child to disease.

That conviction led me to Pediatric Neuroblastoma. I taught myself biology, genomics, molecular pathways, drug discovery -- whatever the work required. I made one commitment early: I would not profit from this.

Today I'm releasing the HCLS AI Factory -- an open-source platform that takes raw patient DNA and produces ranked drug candidates in under 5 hours. On a single desktop workstation.

What it does:

- Takes raw DNA sequencing data (~200 GB) and runs GPU-accelerated variant calling (NVIDIA Parabricks) -- 11.7 million variants in 2-4 hours instead of days
- Annotates variants against ClinVar, AlphaMissense, and a curated knowledge base of 201 genes across 13 therapeutic areas
- Uses RAG-grounded AI reasoning (Milvus + Claude) to identify druggable gene targets
- Generates 100 novel drug candidates via BioNeMo (MolMIM + DiffDock) in under 16 minutes
- Three domain-specific intelligence agents (CAR-T, Medical Imaging, Precision Oncology) with 1,296 tests

The hardware: an NVIDIA DGX Spark. $3,999. Sits on a desk.

The traditional path from patient DNA to drug candidates takes 6-18 months and $50K-500K+ in infrastructure. This does it in an afternoon.

Everything is Apache 2.0. No vendor lock-in. No strings attached.

My hope is simple: a family sitting in a hospital room gets access to the same computational tools that used to require a research institution. A graduate student at any university can run a complete precision medicine pipeline on hardware they can actually afford. A startup can prototype therapeutic hypotheses without a six-figure compute budget.

If you work in genomics, drug discovery, or clinical research -- take it. Use it. Make it better.

GitHub: https://github.com/ajones1923/hcls-ai-factory
Docs: https://hcls-ai-factory.org

#PrecisionMedicine #OpenSource #Genomics #DrugDiscovery #NVIDIA #DGXSpark #AI #Bioinformatics #Apache2

---

## Twitter/X Thread

---

**Tweet 1 (hook):**

I spent 14 years building this. Today I'm releasing it for free.

An open-source platform that goes from raw patient DNA to ranked drug candidates in under 5 hours.

On a $3,999 desktop workstation.

Apache 2.0. No strings attached.

Thread:

**Tweet 2 (what it does):**

The HCLS AI Factory is a 3-stage pipeline:

Stage 1: GPU-accelerated genomics (Parabricks) -- 11.7M variants in 2-4 hrs
Stage 2: RAG-grounded target ID (Milvus + Claude) -- <5 sec per query
Stage 3: AI drug discovery (MolMIM + DiffDock) -- 100 candidates in 16 min

All on one NVIDIA DGX Spark.

**Tweet 3 (the numbers):**

Key numbers:

- 3.56M searchable vectors in Milvus
- 201 genes across 13 therapeutic areas
- 171 druggable targets (85%)
- 1,296 agent tests in 3.78 sec
- Top drug candidate: +39% improvement over seed compound

Traditional approach: 6-18 months, $50K-500K+
This: <5 hours, $3,999

**Tweet 4 (intelligence agents):**

It includes 3 domain-specific intelligence agents:

- CAR-T: cell therapy evidence + comparative analysis
- Imaging: NVIDIA NIM workflows (VISTA-3D, MAISI) + FHIR R4
- Precision Oncology: MTB packets, trial matching, therapy ranking

Each connects to a shared genomic evidence base (3.5M vectors).

**Tweet 5 (why):**

I started this because of Pediatric Neuroblastoma.

My hope: a family in a hospital room gets access to the same tools that used to require a research institution. A grad student can run precision medicine on hardware they can afford.

Take it. Use it. Make it better.

**Tweet 6 (links):**

GitHub: https://github.com/ajones1923/hcls-ai-factory
Docs: https://hcls-ai-factory.org
Demo video: [link]

Apache 2.0 -- free for commercial and non-commercial use.

Built with: NVIDIA Parabricks, BioNeMo, Milvus, Anthropic Claude, RDKit

---

## Hacker News (Show HN)

---

**Title:** Show HN: Open-source platform -- Patient DNA to drug candidates in 5 hours on a $4K workstation

**Body:**

I've been working on this for 14 years, starting from a HPC background and teaching myself genomics along the way. The HCLS AI Factory is an end-to-end precision medicine platform that processes raw DNA sequencing data into ranked novel drug candidates -- entirely on a single NVIDIA DGX Spark ($3,999).

Three stages:

1. GPU genomics (NVIDIA Parabricks): 200GB FASTQ -> 11.7M variant calls in 2-4 hours (vs 1-2 days on CPU)

2. RAG-grounded target identification (Milvus + Claude): 3.56M vectors, 201 genes across 13 therapeutic areas, <5 sec per query

3. AI drug discovery (BioNeMo MolMIM + DiffDock + RDKit): 100 ranked drug candidates in 8-16 minutes

The platform also includes three domain-specific intelligence agents (CAR-T cell therapy, medical imaging, precision oncology) with a combined 1,296 tests.

Traditional path: 6-18 months and $50K-500K+ in infrastructure. This does it in an afternoon on a desktop.

Everything is Apache 2.0. The demo targets VCP/Frontotemporal Dementia and produces candidates with +39% composite improvement over the clinical seed compound.

GitHub: https://github.com/ajones1923/hcls-ai-factory

Tech stack: NVIDIA Parabricks 4.6, Milvus 2.4, BGE-small-en-v1.5, Anthropic Claude, BioNeMo MolMIM/DiffDock, RDKit, Nextflow, Docker, Streamlit

Happy to answer questions about the architecture, benchmarks, or the biology.

---

## Reddit (r/bioinformatics, r/genomics, r/MachineLearning)

---

**Title:** I open-sourced a platform that goes from raw patient DNA to drug candidates in <5 hours on a $4K workstation (Apache 2.0)

**Body:**

After 14 years of work, I'm releasing the HCLS AI Factory -- an end-to-end precision medicine platform that runs entirely on a single NVIDIA DGX Spark ($3,999).

**What it does:**

- Stage 1: GPU-accelerated variant calling via Parabricks (11.7M variants from 30x WGS in 2-4 hours)
- Stage 2: RAG pipeline with Milvus (3.56M vectors), ClinVar, AlphaMissense, and Claude for druggable target identification
- Stage 3: Generative drug discovery via BioNeMo MolMIM + DiffDock + RDKit (100 ranked candidates in 8-16 min)

**Also includes:**

- 3 domain-specific agents (CAR-T, Imaging, Oncology) with 1,296 tests
- 201 genes across 13 therapeutic areas, 85% druggable
- Synthetic demo VCF for quick evaluation
- One-command quickstart: `./quickstart.sh`

**Why it matters:**

The traditional path from patient DNA to a drug lead takes 6-18 months and requires significant infrastructure investment. This brings the entire workflow to a desktop workstation that anyone can afford.

Apache 2.0 -- use it for research, commercial work, whatever you need.

GitHub: https://github.com/ajones1923/hcls-ai-factory
Docs: https://hcls-ai-factory.org

Happy to discuss architecture decisions, benchmark methodology, or answer questions.
