---
title: Architecture
hide:
  - toc
---

# Architecture

The HCLS AI Factory runs three GPU-accelerated AI pipelines on a single NVIDIA DGX Spark workstation — from raw FASTQ sequencing data to ranked drug candidates in under 5 hours.

---

## System Architecture

<div class="architecture-showcase">

![HCLS AI Factory on NVIDIA DGX Spark — Architectural Infographic](diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Architectural%20Infographic.png){ loading=lazy }

</div>

---

## Pipeline Logical Diagram

<div class="architecture-showcase">

![HCLS AI Factory on NVIDIA DGX Spark — Pipeline Logical Diagram](diagrams/dgx_spark/HCLS%20AI%20Factory%20on%20NVIDIA%20DGX%20Spark%20Pipeline%20Logical%20Diagram.png){ loading=lazy }

</div>

---

## Interactive Diagrams

Explore the architecture in detail with interactive draw.io diagrams:

- [**Detailed Architecture** — DGX Spark](diagrams/dgx_spark/HCLS_AI_Factory_DGX_Spark.drawio.html) — Full component-level view
- [**High-Level Architecture** — DGX Spark](diagrams/dgx_spark/HCLS_AI_Factory_DGX_Spark_HighLevel.drawio.html) — Simplified pipeline overview

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
