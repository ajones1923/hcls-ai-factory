---
title: Architecture
hide:
  - toc
---

# Architecture

The HCLS AI Factory runs three GPU-accelerated AI pipelines on a single NVIDIA DGX Spark workstation — from raw FASTQ sequencing data to ranked drug candidates in under 5 hours.

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
