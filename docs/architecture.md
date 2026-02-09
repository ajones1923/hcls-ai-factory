---
title: Architecture
hide:
  - toc
---

# Architecture

The HCLS AI Factory runs three GPU-accelerated AI pipelines on a single NVIDIA DGX Spark workstation — from raw FASTQ sequencing data to ranked drug candidates in under 5 hours.

### Stage 0 + Three Processing Stages

**Stage 0 — Data Acquisition (one-time):** Before the pipeline can run, all required data must be downloaded: HG002 FASTQ sequencing files (~200 GB), the GRCh38 reference genome with BWA-MEM2 index (~11 GB), ClinVar clinical annotations (4.1M variants), and AlphaMissense pathogenicity predictions (71M variants). The [`setup-data.sh`](DATA_SETUP.md) script automates this entire process with checksum verification, parallel downloads, and idempotent resumption. This is a one-time step (~500 GB total).

**Stage 1 — GPU Genomics (120–240 min):** Raw FASTQ sequencing files are aligned to the human reference genome using BWA-MEM2, then variant-called with Google DeepVariant — both accelerated through NVIDIA Parabricks. The output is a clinical-grade VCF containing 11.7 million variants at >99% accuracy.

**Stage 2 — Evidence RAG (interactive):** Variants are annotated against ClinVar (4.1M clinical records) and AlphaMissense (71M pathogenicity predictions), embedded with BGE-small-en-v1.5, and indexed into a Milvus vector database (3.5M vectors). A conversational RAG interface powered by Anthropic Claude lets researchers query variants in natural language, with every answer grounded in retrieved clinical evidence.

**Stage 3 — Drug Discovery (8–16 min):** For validated targets, NVIDIA BioNeMo MolMIM generates novel molecular candidates from a seed compound, DiffDock predicts protein-ligand binding poses, and RDKit scores drug-likeness (QED, Lipinski, synthetic accessibility). The result: 100 ranked candidates per target with full structural and pharmacological profiles.

All three stages run on the DGX Spark's GB10 Grace Blackwell Superchip with 128 GB unified memory — connected via NVLink-C2C so the GPU and CPU share the same memory pool without transfer bottlenecks.

---

## System Architecture

<div class="architecture-showcase" markdown="1">

![HCLS AI Factory on NVIDIA DGX Spark — Architectural Infographic](../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Architectural%20Infographic%20v1.0.png){ loading=lazy }

</div>

---

## Architectural Infographic (Alt View)

<div class="architecture-showcase" markdown="1">

![HCLS AI Factory on NVIDIA DGX Spark — Architectural Infographic Alt View](../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Architectural%20Infographic%20v1.0%20%28Alt%20View%29.png){ loading=lazy }

</div>

---

## Pipeline Logical Diagram

<div class="architecture-showcase" markdown="1">

![HCLS AI Factory on NVIDIA DGX Spark — Pipeline Logical Diagram](../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Pipeline%20Logical%20Diagram%20v1.0.png){ loading=lazy }

</div>

---

## draw.io Diagrams

### High Level

<div class="architecture-showcase" markdown="1">

![HCLS AI Factory — draw.io High Level](../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20drawio%20Diagram-High%20Level%20v1.0.png){ loading=lazy }

</div>

### Medium Level

<div class="architecture-showcase" markdown="1">

![HCLS AI Factory — draw.io Medium Level](../diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20drawio%20Diagram-Medium%20Level%20v1.0.png){ loading=lazy }

</div>

---

## Interactive Diagrams

Explore the architecture in detail with interactive draw.io diagrams:

- [**Detailed Architecture** — DGX Spark](../diagrams/dgx_spark/HCLS_AI_Factory_DGX_Spark.drawio.html) — Full component-level view
- [**High-Level Architecture** — DGX Spark](../diagrams/dgx_spark/HCLS_AI_Factory_DGX_Spark_HighLevel.drawio.html) — Simplified pipeline overview

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
| **CPU** | 144 ARM64 cores (Grace) |
| **Storage** | NVMe SSD |
| **Price** | $3,999 |

---

## Deep Dives

- [**Project Bible**](HCLS_AI_FACTORY_PROJECT_BIBLE.md) — Complete technical reference with scoring formulas, thresholds, and configurations
- [**White Paper**](HCLS_AI_FACTORY_WHITE_PAPER_DGX_SPARK.md) — Architecture overview and design rationale
- [**Architecture Mindmap**](HCLS_AI_Factory_Mindmap_Open.md) — Visual system map
