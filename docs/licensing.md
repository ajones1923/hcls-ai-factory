# Licensing and Cost Guide

> What's free, what requires API keys, and what changes at enterprise scale.

---

## HCLS AI Factory Code

All code in this repository — pipeline scripts, orchestrator, portal, documentation — is released under the **Apache 2.0 license**. You can use, modify, and distribute it freely.

---

## Third-Party Component Licenses

The HCLS AI Factory integrates open-source tools, NVIDIA software, and a commercial API. Each has its own license terms.

### Open-Source Components (No Cost)

| Component | Technology | License |
|---|---|---|
| Genomics alignment | BWA-MEM2 | MIT |
| Variant calling | Google DeepVariant | BSD-3-Clause |
| Vector database | Milvus 2.4 | Apache 2.0 |
| Embeddings | BGE-small-en-v1.5 (BAAI) | MIT |
| Chemistry toolkit | RDKit | BSD-3-Clause |
| Pipeline orchestration | Nextflow | Apache 2.0 |
| Web interface | Streamlit | Apache 2.0 |
| Monitoring | Grafana / Prometheus | AGPLv3 / Apache 2.0 |
| Containers | Docker | Apache 2.0 |
| PDF generation | ReportLab | BSD |
| VCF parsing | cyvcf2 | MIT |
| 3D visualization | py3Dmol | MIT |

These are fully open-source. No API keys, no license fees, no restrictions on commercial use.

### NVIDIA Software

| Component | License | Dev/Test Cost | Production Cost |
|---|---|---|---|
| **Parabricks 4.x** | NVIDIA EULA | **Free** — no license required | **Free** — enterprise support optional via NVAIE |
| **BioNeMo MolMIM** (NIM) | NVIDIA AI Product Agreement | **Free** on DGX Spark | NVAIE subscription required |
| **BioNeMo DiffDock** (NIM) | NVIDIA AI Product Agreement | **Free** on DGX Spark | NVAIE subscription required |

**Key details:**

- **Parabricks** is free for all users. No license is required for Parabricks 4.x and later. Enterprise support is available through NVIDIA AI Enterprise (NVAIE) but is not required to run the software.

- **BioNeMo NIMs** (MolMIM and DiffDock) are available at no cost for research, development, and testing through the [NVIDIA Developer Program](https://developer.nvidia.com/nim). On DGX Spark and other Grace Blackwell client systems, NIMs designated as Free SDKs may be used without a subscription when used on a personal workstation and not serving multiple users. Moving to **production** (hospital deployments, multi-user systems, or enterprise scale) requires an [NVIDIA AI Enterprise subscription](https://www.nvidia.com/en-us/data-center/products/ai-enterprise/) — currently $4,500/GPU/year.

- **NGC API Key** is required to pull NIM containers. Free to obtain through the [NVIDIA Developer Program](https://developer.nvidia.com).

### Commercial API

| Component | Provider | License | Cost |
|---|---|---|---|
| **Claude** (RAG reasoning) | Anthropic | Commercial API | Pay-per-token ([pricing](https://www.anthropic.com/pricing)) |

An Anthropic API key is required for Stage 2 RAG chat. Claude is a commercial service — not open-source, not free. API costs depend on usage volume.

---

## Cost Summary by Deployment Phase

| Phase | Hardware | Software Licenses | API Keys |
|---|---|---|---|
| **Phase 1 — Proof Build** | DGX Spark ($3,999) | None required | Anthropic + NGC (free) |
| **Phase 2 — Departmental** | DGX B200 ($500K-$1M) | NVAIE ($4,500/GPU/year) | Anthropic + NGC |
| **Phase 3 — Enterprise** | DGX SuperPOD ($7M-$60M+) | NVAIE ($4,500/GPU/year) | Anthropic + NGC |

For a Phase 1 proof build on DGX Spark, the only recurring costs are Anthropic API usage. All NVIDIA software components — Parabricks, MolMIM, and DiffDock — are available at no additional cost for research and development.

---

## What "Open Source" Means

The HCLS AI Factory is an open-source project: all code authored for this platform is Apache 2.0, published on GitHub, and free to use.

The platform **depends on** components with different license models. The open-source tools (Milvus, RDKit, Nextflow, DeepVariant, BWA-MEM2) are fully free. The NVIDIA NIM containers are free for development on DGX Spark but require enterprise licensing at scale. The Anthropic Claude API is a paid commercial service.

This is a common pattern in modern AI platforms — open-source orchestration and code, with commercial or proprietary components at the inference layer.

---

## References

- [NVIDIA Parabricks — Free Download](https://www.nvidia.com/en-us/industries/healthcare-life-sciences/genomics/free-trial-vs-buy/)
- [NVIDIA NIM for Developers](https://developer.nvidia.com/nim)
- [NVIDIA AI Enterprise Pricing](https://docs.nvidia.com/ai-enterprise/planning-resource/licensing-guide/latest/pricing.html)
- [NVIDIA AI Product Terms](https://www.nvidia.com/en-us/agreements/enterprise-software/product-specific-terms-for-ai-products/)
- [Anthropic API Pricing](https://www.anthropic.com/pricing)
- [Apache 2.0 License](https://www.apache.org/licenses/LICENSE-2.0)
