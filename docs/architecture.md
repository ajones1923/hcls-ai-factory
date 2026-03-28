---
title: Architecture
hide:
  - toc
search:
  boost: 2
tags:
  - Architecture
  - DGX Spark
  - GB10
  - NVLink
  - Pipeline
---

# Architecture Guide

The HCLS AI Factory runs three GPU-accelerated AI pipelines on a single NVIDIA DGX Spark workstation — from raw FASTQ sequencing data to ranked drug candidates in under 5 hours.

### Processing Stages

**Stage 0 — Data Acquisition (one-time):** Before the pipeline can run, all required data must be downloaded: HG002 FASTQ sequencing files (~200 GB), the GRCh38 reference genome with BWA-MEM2 index (~11 GB), ClinVar clinical annotations (4.1M variants), and AlphaMissense pathogenicity predictions (71M variants). The [`setup-data.sh`](DATA_SETUP.md) script automates this entire process with checksum verification, parallel downloads, and idempotent resumption. This is a one-time step (~500 GB total).

**Stage 1 — GPU Genomics (120–240 min):** Raw FASTQ sequencing files are aligned to the human reference genome using BWA-MEM2, then variant-called with Google DeepVariant — both accelerated through NVIDIA Parabricks. The output is a clinical-grade VCF containing 11.7 million variants at >99% accuracy.

**Stage 2 — Evidence RAG (interactive):** Variants are annotated against ClinVar (4.1M clinical records) and AlphaMissense (71M pathogenicity predictions), embedded with BGE-small-en-v1.5, and indexed into a Milvus vector database (3.56M vectors). A conversational RAG interface powered by Anthropic Claude lets researchers query variants in natural language, with every answer grounded in retrieved clinical evidence.

**Stage 3 — Drug Discovery (8–16 min):** For validated targets, NVIDIA BioNeMo MolMIM generates novel molecular candidates from a seed compound, DiffDock predicts protein-ligand binding poses, and RDKit scores drug-likeness (QED, Lipinski, synthetic accessibility). The result: 100 ranked candidates per target with full structural and pharmacological profiles.

All three stages run on the DGX Spark's GB10 Grace Blackwell Superchip with 128 GB unified memory — connected via NVLink-C2C so the GPU and CPU share the same memory pool without transfer bottlenecks.

---

## Architectural Infographic

<div class="architecture-showcase">
  <a href="../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Architectural%20Infographic%20v1.0.png" class="lightbox-link" target="_blank">
    <img src="../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Architectural%20Infographic%20v1.0.png" alt="HCLS AI Factory on NVIDIA DGX Spark — Architectural Infographic" class="architecture-img" loading="lazy">
  </a>
  <p class="architecture-caption">Click image to expand</p>
</div>

---

## From Patient DNA to New Medicine Infographic

<div class="architecture-showcase">
  <a href="../diagrams/dgx_spark/From%20Patient%20DNA%20to%20New%20Medicine%20Infographic.png" class="lightbox-link" target="_blank">
    <img src="../diagrams/dgx_spark/From%20Patient%20DNA%20to%20New%20Medicine%20Infographic.png" alt="From Patient DNA to New Medicine Infographic" class="architecture-img" loading="lazy">
  </a>
  <p class="architecture-caption">Click image to expand</p>
</div>

---

## Pipeline Mindmap

[:material-file-pdf-box: View Pipeline Mindmap (PDF)](../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Pipeline%20Mindmap%20v1.0.pdf)

---

## Hardware Platform

| Component | Specification |
|---|---|
| **System** | NVIDIA DGX Spark |
| **GPU** | GB10 Grace Blackwell Superchip |
| **Memory** | 128 GB unified LPDDR5x |
| **CPU** | ARM64 cores (Grace) |
| **Storage** | NVMe SSD |
| **Price** | $4,699 |

---

## Deep Dives

- [**Project Bible**](HCLS_AI_FACTORY_PROJECT_BIBLE.md) — Complete technical reference with scoring formulas, thresholds, and configurations
- [**White Paper**](HCLS_AI_FACTORY_WHITE_PAPER_DGX_SPARK.md) — Architecture overview and design rationale

---

!!! warning "Clinical Decision Support Disclaimer"
    The HCLS AI Factory platform and its components are clinical decision support research tools. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
