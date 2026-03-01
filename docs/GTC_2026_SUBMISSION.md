---
tags: [GTC, NVIDIA, Conference, Submission]
---

# NVIDIA GTC 2026 — Talk Submission

## Session Information

| Field | Value |
|---|---|
| Session Type | Technical Talk (50 min) |
| Track | Healthcare & Life Sciences |
| Sub-track | Drug Discovery / Genomics |
| Difficulty | Intermediate |
| Target Audience | Computational biologists, bioinformaticians, pharma R&D, clinical researchers |

## Title

From Patient DNA to Drug Candidates in 5 Hours on DGX Spark

## Abstract

What if a family facing a devastating diagnosis could go from a DNA sample to ranked drug candidates in under five hours — on a single $3,999 machine sitting on a desk?

The HCLS AI Factory is an open-source, Apache 2.0-licensed precision medicine platform that runs end-to-end on the NVIDIA DGX Spark. It chains three GPU-accelerated stages into a seamless pipeline: (1) Genomics — Parabricks 4.6 with BWA-MEM2 and DeepVariant calls 11.7 million variants from 30x whole-genome sequencing in 120-240 minutes; (2) RAG/Chat — a retrieval-augmented generation engine backed by Milvus (3.56 million vectors across ClinVar, AlphaMissense, and curated literature) and Claude delivers clinician-ready variant interpretation in under five seconds per query; (3) Drug Discovery — BioNeMo NIMs (MolMIM and DiffDock) generate and dock 100 ranked small-molecule candidates in 8-16 minutes, with the top candidate showing a 39% improvement over the established seed compound CB-5083.

In this session, attendees will see the full pipeline demonstrated live on DGX Spark hardware, walk through the architecture decisions that make it possible on a single node with 128 GB of unified memory, and learn how to fork the repository and deploy the platform in their own labs. We will cover three intelligence agents — CAR-T Cell Therapy, Imaging, and Precision Oncology — that extend the platform into specialized clinical domains, validated across 1,296 tests. The talk closes with a patient advocacy perspective: how democratizing this technology shifts power from institutions to families.

## Description

Precision medicine has a bottleneck problem. The science exists to connect a patient's genomic profile to targeted therapies, but the infrastructure to do so has historically required months of calendar time, teams of specialists, and budgets measured in hundreds of thousands of dollars. A family confronting frontotemporal dementia, a rare cancer, or an inherited cardiomyopathy does not have months. They need answers now, and they need those answers to be actionable.

The HCLS AI Factory was built to close that gap. It is an end-to-end precision medicine platform — fully open source under the Apache 2.0 license — that takes raw patient DNA (FASTQ files from 30x whole-genome sequencing) and produces 100 ranked drug candidates in under five hours. The entire pipeline runs on a single NVIDIA DGX Spark: a GB10 GPU with 128 GB of unified LPDDR5x memory, 20 ARM cores, and a price tag of $3,999. No cluster. No cloud bill. No vendor lock-in.

This talk walks through the platform in three acts, mirroring its three-stage architecture.

**Act I — Genomics (Stage 1).** We begin with approximately 200 GB of raw FASTQ data from the Genome in a Bottle HG002 reference sample. Using NVIDIA Parabricks 4.6 — BWA-MEM2 for alignment, DeepVariant for variant calling — the DGX Spark produces a fully annotated VCF containing 11.7 million variants in 120-240 minutes. We will discuss memory management on unified LPDDR5x, how NVLink-C2C eliminates the traditional CPU-GPU transfer penalty, and why this matters for desktop-scale genomics.

**Act II — RAG/Chat (Stage 2).** The VCF feeds into a retrieval-augmented generation engine. Milvus indexes 3.56 million vectors derived from ClinVar (~2.7 million variant records), AlphaMissense (71 million missense predictions), and curated clinical literature. The BGE-small-en-v1.5 embedding model and Claude work together to answer natural-language questions about the patient's variants — pathogenicity, drug interactions, clinical trial eligibility — in under five seconds per query. We cover 201 genes across 13 therapeutic areas. The live demo will show a clinician-style interrogation of variants in the VCP gene, the target for our frontotemporal dementia case study.

**Act III — Drug Discovery (Stage 3).** Actionable variants become drug targets. BioNeMo NIMs — specifically MolMIM for molecular generation and DiffDock for protein-ligand docking — produce 100 candidate small molecules in 8-16 minutes. RDKit scores each candidate on drug-likeness, synthetic accessibility, and predicted binding affinity. In the VCP/frontotemporal dementia demo, the top-ranked candidate shows a 39% improvement in composite score over CB-5083, an established VCP inhibitor, docked against PDB structures 5FTK, 8OOI, 9DIL, and 7K56.

Beyond the core pipeline, we introduce three intelligence agents that extend the platform into specialized clinical domains: a CAR-T Cell Therapy agent for immunotherapy target selection, an Imaging Intelligence agent for radiogenomic correlation, and a Precision Oncology agent for tumor-specific treatment recommendations. Together, these agents have been validated across 1,296 tests.

The session concludes with the human story. This platform was not built for pharmaceutical companies with seven-figure IT budgets. It was built so that a parent, a caregiver, or a rural oncologist can plug in a DGX Spark, load a genome, and have a conversation with their patient's DNA — then walk into a treatment meeting with data. That is what democratization looks like.

Attendees will leave with the GitHub repository URL, a five-minute quickstart guide, and the knowledge to deploy the HCLS AI Factory in their own environment.

## Key Takeaways

- **Desktop-scale genomics is here.** A single DGX Spark ($3,999) can run a complete precision medicine pipeline — from raw FASTQ to ranked drug candidates — in under five hours, replacing workflows that traditionally require clusters, cloud infrastructure, and weeks of calendar time.

- **The NVIDIA life sciences stack is production-ready end-to-end.** Parabricks 4.6 for variant calling, BioNeMo NIMs for molecular generation and docking, and DCGM for GPU monitoring compose a coherent, GPU-accelerated pipeline with no gaps.

- **RAG transforms variant interpretation.** By indexing 3.56 million vectors from ClinVar, AlphaMissense, and curated literature in Milvus, clinicians and researchers can interrogate a patient's 11.7 million variants through natural language — no bioinformatics expertise required.

- **Open source lowers the barrier to zero.** The entire platform is Apache 2.0 licensed. Attendees can fork the repository, load their own genomes, and be running within minutes on any DGX Spark.

- **Precision medicine is a patient advocacy tool.** When families can generate their own genomic insights and drug candidates, the power dynamic in treatment decisions shifts from institutional gatekeepers to informed advocates.

## Speaker Bio

Adam Jones brings more than 14 years of hands-on experience in genomic research, bioinformatics, and precision medicine. He is the creator of the HCLS AI Factory, an open-source platform that compresses the precision medicine pipeline — from raw patient DNA to ranked drug candidates — into under five hours on a single NVIDIA DGX Spark.

Adam's career spans the full arc of modern genomics: from early next-generation sequencing workflows to large-scale variant interpretation, clinical annotation, and GPU-accelerated drug discovery. He has built production pipelines that bridge the gap between raw sequencing data and clinical decision-making, with deep expertise in variant calling (Parabricks, GATK, DeepVariant), vector-based retrieval systems (Milvus, ClinVar, AlphaMissense), and AI-driven molecular generation (BioNeMo NIMs).

What sets Adam apart is his focus on accessibility. The HCLS AI Factory was not designed for institutions with unlimited compute budgets — it was designed for the smallest credible unit of precision medicine: one patient, one machine, one afternoon. His work is motivated by a belief that families facing rare and devastating diagnoses deserve the same computational firepower that pharmaceutical giants take for granted. Adam is a practical builder who ships open-source tools that work, documented to the level where a motivated clinician or researcher can deploy them without a DevOps team.

## Demo Plan

The live demonstration follows the three-stage pipeline in real time, with pre-computed checkpoints to ensure the session stays within the 50-minute window.

**Screen 1 — Landing Page & Service Health (Port 8080).** Open the Flask-based service dashboard. All services green: Genomics (5000), RAG API (5001), Chat UI (8501), Drug Discovery UI (8505), Portal (8510), Milvus (19530), MolMIM (8001), DiffDock (8002). Show Grafana (3000) with GPU utilization via DCGM (9400).

**Screen 2 — Genomics Pipeline (Port 5000).** Launch a Parabricks 4.6 run against HG002 30x WGS data (~200 GB FASTQ). Show the BWA-MEM2 alignment phase, then cut to a pre-computed VCF with 11.7 million variants. Highlight GPU memory utilization on the DGX Spark's unified 128 GB pool.

**Screen 3 — RAG/Chat Interface (Port 8501).** Load the Streamlit chat UI. Ask natural-language questions about VCP gene variants: "What pathogenic variants exist in VCP?", "What drugs target p97?", "Is this patient eligible for any clinical trials?" Show sub-5-second response times. Demonstrate vector similarity search in Attu (8000) — visualize the 3.56 million indexed vectors.

**Screen 4 — Drug Discovery UI (Port 8505).** Trigger MolMIM molecular generation from the VCP target. Show DiffDock docking against PDB 5FTK. Display the ranked results table: 100 candidates scored by binding affinity, drug-likeness (QED), and synthetic accessibility. Highlight the top candidate's 39% improvement over CB-5083.

**Screen 5 — Intelligence Agents.** Briefly demonstrate the Precision Oncology agent (Port 8510) with a tumor-specific query. Show the agent selecting relevant genes from the 201-gene panel and cross-referencing 13 therapeutic areas.

**Fallback plan:** All demo segments have pre-recorded video backups. If the live DGX Spark encounters any issue, the pre-recorded segments maintain session flow without interruption.

## Technical Requirements

| Requirement | Details |
|---|---|
| **On-stage hardware** | 1x NVIDIA DGX Spark (GB10 GPU, 128 GB unified memory, 20 ARM cores) — preferred live; pre-recorded demo acceptable as fallback |
| **Network** | Internet access required for Claude API calls during the RAG/Chat demo segment; all other pipeline components run locally |
| **Display** | Minimum 1920x1080 resolution; 4K preferred for code and terminal readability |
| **Projector/screen** | Single screen sufficient; dual screen preferred (one for slides, one for live demo terminal) |
| **Audio/video** | Standard conference A/V; lapel microphone for mobility during demo |
| **Backup media** | USB drive with pre-recorded demo segments (1080p MP4) for all five demo screens |
| **Software on demo machine** | Docker, Nextflow, Python 3.10+, all HCLS AI Factory services pre-deployed and tested |
| **Setup time** | 30 minutes pre-session for hardware check, service verification, and network test |

## Why This Talk Matters for GTC

The DGX Spark was announced as NVIDIA's entry point for personal AI computing — this talk proves it is a serious life sciences instrument, not just a developer workstation. By running a complete precision medicine pipeline end-to-end on a single Spark, we demonstrate that NVIDIA's hardware and software stack (Parabricks, BioNeMo NIMs, DCGM, NVLink-C2C) can deliver clinical-grade results at a price point that puts genomic analysis within reach of individual researchers, small labs, and patient advocacy groups. The platform is fully open source under Apache 2.0, giving every NVIDIA customer in healthcare and life sciences a ready-to-fork reference architecture. Most importantly, the talk carries a human story — families using a $3,999 machine to advocate for their loved ones' treatments — that resonates far beyond the technical audience and reinforces NVIDIA's mission to bring accelerated computing to the problems that matter most.

## Session Outline (50 min)

| Time | Section | Content |
|---|---|---|
| 0-5 min | Introduction | The precision medicine gap: why 6-18 months and $50K-500K+ is unacceptable. A family's story. The thesis — what if it took 5 hours and $3,999? |
| 5-15 min | Architecture | Three-stage pipeline overview. DGX Spark hardware deep dive: GB10 GPU, 128 GB unified LPDDR5x, NVLink-C2C, 20 ARM cores. Why unified memory changes everything for genomics. |
| 15-20 min | Live Demo Part 1 | Genomics pipeline launch (Parabricks 4.6, BWA-MEM2, DeepVariant). RAG/Chat interface — natural-language variant interrogation of VCP gene. Sub-5-second query responses. |
| 20-30 min | Deep Dive | RAG engine architecture: Milvus vector store (3.56M vectors), BGE-small-en-v1.5 embeddings, ClinVar + AlphaMissense data integration, Claude for clinical reasoning. 201 genes, 13 therapeutic areas. |
| 30-35 min | Live Demo Part 2 | Drug Discovery pipeline: MolMIM molecular generation, DiffDock protein-ligand docking, RDKit scoring. 100 ranked candidates. Top hit: +39% over CB-5083 against VCP/p97. |
| 35-42 min | Intelligence Agents | CAR-T Cell Therapy agent for immunotherapy targets. Imaging Intelligence agent for radiogenomic correlation. Precision Oncology agent for tumor-specific recommendations. 1,296 validated tests across all three agents. |
| 42-47 min | Impact & Future | Patient advocacy and democratization. Roadmap: multi-patient batching, federated learning, clinical validation pathways. The $3,999 precision medicine lab. How to get started — GitHub repo, quickstart guide, community. |
| 47-50 min | Q&A | Open questions. GitHub URL on screen. |

## Supplementary Materials

| Material | Description | License |
|---|---|---|
| **GitHub Repository** | Complete HCLS AI Factory source code, Dockerfiles, Nextflow orchestrator, and deployment scripts | Apache 2.0 |
| **White Paper** | _HCLS AI Factory: Desktop-Scale Precision Medicine on DGX Spark_ — architecture, benchmarks, and clinical validation methodology | Apache 2.0 |
| **Demo Guide** | Step-by-step walkthrough of the five-minute quickstart flow, from `docker compose up` to ranked drug candidates | Apache 2.0 |
| **Deployment Guide** | DGX Spark-specific deployment instructions, memory tuning, and service configuration | Apache 2.0 |
| **Architecture Diagrams** | Pipeline flowcharts, data flow diagrams, and system architecture visuals (Draw.io source files included) | Apache 2.0 |
| **Intelligence Report** | Comprehensive validation results for all three intelligence agents (1,296 tests) | Apache 2.0 |
| **Video Script** | Production-ready script for the pre-recorded demo backup and promotional materials | Apache 2.0 |
