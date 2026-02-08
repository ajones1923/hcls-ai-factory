# HCLS AI Factory Documentation

## Quick Navigation

| Document | Description |
|----------|-------------|
| [DATA_SETUP.md](DATA_SETUP.md) | **Start here** — Download and verify all required data (~500 GB) |
| [PRODUCT_DOCUMENTATION.txt](PRODUCT_DOCUMENTATION.txt) | **Complete reference** (3,300+ lines) — Installation, configuration, API reference, troubleshooting |
| [COMPLETE_PIPELINE_OVERVIEW.txt](COMPLETE_PIPELINE_OVERVIEW.txt) | Visual pipeline diagrams and data flow |
| [ARCHITECTURE_MINDMAP.md](ARCHITECTURE_MINDMAP.md) | System architecture overview |
| [ARCHITECTURE.mmd](ARCHITECTURE.mmd) | Mermaid diagram source |
| [PIPELINE_REPORT.md](PIPELINE_REPORT.md) | Pipeline integration report |

---

## Pipeline-Specific Documentation

### Genomics Pipeline (Stage 1)

| Document | Description |
|----------|-------------|
| [genomics-pipeline/README.md](genomics-pipeline/README.md) | Overview and quick start |
| [genomics-pipeline/QUICKSTART.md](genomics-pipeline/QUICKSTART.md) | Step-by-step setup guide |
| [genomics-pipeline/WEB-PORTAL-GUIDE.md](genomics-pipeline/WEB-PORTAL-GUIDE.md) | Web interface documentation |
| [genomics-pipeline/PIPELINE_OVERVIEW.md](genomics-pipeline/PIPELINE_OVERVIEW.md) | Technical pipeline details |
| [genomics-pipeline/PHASES_OVERVIEW.md](genomics-pipeline/PHASES_OVERVIEW.md) | Pipeline phases explained |

### RAG/Chat Pipeline (Stage 2)

| Document | Description |
|----------|-------------|
| [rag-chat-pipeline/README.md](rag-chat-pipeline/README.md) | Overview and configuration |
| [rag-chat-pipeline/PHASE4_SUMMARY.md](rag-chat-pipeline/PHASE4_SUMMARY.md) | Implementation summary |
| [rag-chat-pipeline/PHASE4_BUILD_PLAN.md](rag-chat-pipeline/PHASE4_BUILD_PLAN.md) | Build and deployment plan |

### Drug Discovery Pipeline (Stage 3)

| Document | Description |
|----------|-------------|
| [drug-discovery-pipeline/README.md](drug-discovery-pipeline/README.md) | BioNeMo integration, molecule generation, docking |

---

## Diagrams

All Mermaid (.mmd) diagrams are in the [diagrams/](diagrams/) folder:

| Diagram | Description |
|---------|-------------|
| `genomics-pipeline-diagram.mmd` | Genomics pipeline flow |
| `phase4-pipeline-diagram.mmd` | RAG/Chat pipeline flow |
| `drug-discovery-pipeline-diagram.mmd` | Drug discovery flow |
| `system-architecture.mmd` | Overall system architecture |
| `data-flow.mmd` | Data flow between stages |
| `sequence-diagram.mmd` | Process sequence diagram |

To render Mermaid diagrams:
- GitHub renders them automatically in markdown
- Use [Mermaid Live Editor](https://mermaid.live/)
- VS Code with Mermaid extension

---

## Key Sections in Product Documentation

The main `PRODUCT_DOCUMENTATION.txt` contains:

1. **Executive Summary** — Key metrics and use case
2. **Product Overview** — Vision, users, use cases
3. **System Architecture** — Diagrams, topology, technology stack
4. **Genomics Pipeline** — BWA-MEM2, DeepVariant, benchmarks
5. **RAG/Chat Pipeline** — Milvus, ClinVar, AlphaMissense, Claude AI
6. **Drug Discovery Pipeline** — BioNeMo MolMIM, DiffDock, scoring
7. **HLS Orchestrator** — Service management, Nextflow
8. **Installation & Deployment** — Requirements, setup steps
9. **Configuration Reference** — Environment variables, Docker services
10. **User Guide** — Web interfaces, CLI, workflow examples
11. **Monitoring & Observability** — Grafana, Prometheus, GPU metrics
12. **Security Considerations** — API keys, network, data privacy
13. **Troubleshooting** — Common issues, diagnostics
14. **Appendices** — Glossary (50+ terms), 201 gene list, file formats

---

## Getting Started

1. Start with the main [README.md](../README.md) in the repository root
2. Configure environment: `cp .env.example .env` and add your API keys
3. **Download all required data** (~500 GB, one-time): `./setup-data.sh --all` — see [DATA_SETUP.md](DATA_SETUP.md) for details and troubleshooting
4. Run `./start-services.sh` to launch all services
5. Access the landing page at http://localhost:8080

---

## File Formats Reference

| Format | Extension | Description |
|--------|-----------|-------------|
| FASTQ | `.fastq.gz` | Raw sequencing reads |
| BAM | `.bam` | Aligned sequences |
| VCF | `.vcf.gz` | Variant calls |
| SMILES | `.smi` | Molecular structures |
| PDB | `.pdb` | Protein structures |

---

## Support

- Open an issue on GitHub for bugs or questions
- See [Troubleshooting](PRODUCT_DOCUMENTATION.txt) section 13 for common issues

---

*HCLS AI Factory — From Patient DNA to Drug Candidates*
