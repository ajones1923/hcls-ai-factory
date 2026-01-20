# HCLS AI Factory

**Transform patient DNA into therapeutic candidates in hours, not months.**

---

## Origin

In 2012, I set out to use my high-performance computing background for something that mattered. I started with one conviction: no parent should ever have to lose a child to disease.

That conviction led me to Pediatric Neuroblastoma. I taught myself biology, genomics, molecular pathways, drug discovery—whatever the work required. I made one commitment early: I would not profit from this. Whatever I built, I would give away freely, so others could build on it and move faster than any one person ever could alone.

Thousands of hours later, this is the result.

— **Adam Jones**

---

## What It Is

The Healthcare & Life Sciences (HCLS) AI Factory unifies three production-grade AI workflows into a single, continuous system—designed to take raw patient DNA and produce viable drug candidates without the fragmentation and delays that define traditional approaches.

Raw FASTQ files flow through NVIDIA Parabricks for GPU-accelerated genomics—alignment, variant calling, clinical-grade accuracy via DeepVariant—completing in hours instead of days. Outputs feed directly into an evidence layer where millions of variants can be queried in natural language, grounded in ClinVar, AlphaMissense, structural data, and curated biomedical knowledge. Validated targets then move into generative drug discovery via NVIDIA BioNeMo, where novel molecules are created, docked, scored, and ranked.

No batch jobs. No manual handoffs. Full lineage from patient DNA to candidate therapeutic.

One workstation. One workflow. Hours, not months.

Take it. Use it. Make it better.

---

## What It Does

| Stage | Pipeline | Input | Output | Key Technology |
|-------|----------|-------|--------|----------------|
| 1 | **Genomics** | FASTQ (raw sequences) | VCF (variant calls) | NVIDIA Parabricks, DeepVariant |
| 2 | **RAG/Chat** | VCF + natural language query | Target hypothesis | Milvus, Claude AI, ClinVar |
| 3 | **Drug Discovery** | Protein target | Ranked drug candidates | BioNeMo MolMIM, DiffDock |

### Performance

- **Genome analysis**: 30-90 minutes (vs. 24-48 hours on CPU)
- **Variant database**: 3.5 million searchable with semantic search
- **Query response**: <5 seconds for complex natural language
- **Molecule generation**: 10-100 candidates in seconds

---

## Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           HCLS AI FACTORY                                   │
├─────────────────┬─────────────────┬─────────────────────────────────────────┤
│                 │                 │                                         │
│   GENOMICS      │    RAG/CHAT     │         DRUG DISCOVERY                  │
│   PIPELINE      │    PIPELINE     │            PIPELINE                     │
│                 │                 │                                         │
│  FASTQ → VCF    │  VCF → Target   │      Target → Molecules                 │
│                 │                 │                                         │
│  • Parabricks   │  • Milvus       │      • BioNeMo MolMIM                   │
│  • BWA-MEM2     │  • ClinVar      │      • BioNeMo DiffDock                 │
│  • DeepVariant  │  • Claude AI    │      • RDKit scoring                    │
│                 │                 │                                         │
└────────┬────────┴────────┬────────┴──────────────────┬──────────────────────┘
         │                 │                           │
         ▼                 ▼                           ▼
   Web Portal:8080   Chat UI:8501              Discovery UI:8505
```

---

## Knowledge Coverage

- **201 genes** across 13 therapeutic areas
- **171 druggable targets** (85% druggability rate)
- **4.1M variants** from ClinVar (clinical evidence)
- **71M predictions** from AlphaMissense (AI pathogenicity)

### Therapeutic Areas

Oncology | Neurology | Rare Disease | Cardiovascular | Immunology | Pharmacogenomics | Metabolic | Ophthalmology | Dermatology | Pulmonary | Infectious Disease | Hematology | Musculoskeletal

---

## Demo: VCP and Frontotemporal Dementia

The platform includes a complete demonstration using VCP (Valosin-containing protein) as a target for Frontotemporal Dementia:

1. **Genomic Analysis** — Process HG002 whole-genome sequencing data
2. **Variant Discovery** — Identify 13 VCP variants including rs188935092
3. **Evidence Synthesis** — Connect variants to FTD through knowledge graph
4. **Structure Retrieval** — Access Cryo-EM structures (8OOI, 9DIL, 7K56, 5FTK)
5. **Molecule Generation** — Create novel VCP inhibitor candidates
6. **Binding Prediction** — Dock molecules to ATP-binding pocket
7. **Report Generation** — Export PDF with ranked candidates

---

## Requirements

### Hardware

| Component | Recommended | Minimum |
|-----------|-------------|---------|
| GPU | NVIDIA DGX Spark (128GB) | NVIDIA GPU with 24GB+ VRAM |
| RAM | 512GB | 64GB |
| Storage | High-performance NVMe | 2TB SSD |

### Software

- Ubuntu 22.04 LTS
- Docker with NVIDIA Container Runtime
- NVIDIA CUDA 12.x
- NGC CLI (for Parabricks and BioNeMo)

### API Keys

- **NVIDIA NGC** — Required for Parabricks and BioNeMo NIMs
- **Anthropic** — Required for Claude AI in RAG pipeline

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/ajones1923/hcls-ai-factory.git
cd hcls-ai-factory

# Configure environment
cp .env.example .env
# Edit .env with your NGC and Anthropic API keys

# Start all services
./start-services.sh

# Access the landing page
open http://localhost:8080
```

### Service Ports

| Service | Port | Description |
|---------|------|-------------|
| Landing Page | 8080 | Main entry point |
| Chat UI | 8501 | RAG-powered variant queries |
| Discovery UI | 8505 | Drug candidate generation |
| Portal | 8510 | Pipeline orchestration |
| Grafana | 3000 | Monitoring dashboards |

---

## Documentation

- [Installation Guide](docs/INSTALLATION.md)
- [Genomics Pipeline](docs/GENOMICS_PIPELINE.md)
- [RAG/Chat Pipeline](docs/RAG_CHAT_PIPELINE.md)
- [Drug Discovery Pipeline](docs/DRUG_DISCOVERY_PIPELINE.md)
- [Configuration Reference](docs/CONFIGURATION.md)
- [API Reference](docs/API_REFERENCE.md)

---

## Technology Stack

**Compute**: NVIDIA DGX Spark, CUDA 12.x

**Genomics**: NVIDIA Parabricks 4.6, BWA-MEM2, DeepVariant

**AI/ML**: Anthropic Claude, NVIDIA BioNeMo NIMs, HuggingFace Transformers

**Databases**: Milvus (vectors), ClinVar, AlphaMissense, RCSB PDB

**Chemistry**: RDKit, BioNeMo MolMIM, BioNeMo DiffDock

**Infrastructure**: Docker, Nextflow, Grafana, Prometheus

---

## Use Cases

**Pharmaceutical R&D**
- Accelerate target identification from months to days
- Generate novel chemical matter for hit-to-lead programs
- Rapid hypothesis testing across genomic datasets

**Research Institutions**
- Process patient cohorts at scale
- Natural language queries over genomic evidence
- Publication-ready analyses and visualizations

**Healthcare Organizations**
- Clinical-grade variant calling (>99% accuracy)
- Actionable therapeutic insights from patient genomes
- Compliance-ready documentation

---

## License

This project is released under the [Apache License 2.0](LICENSE).

You are free to use, modify, and distribute this software. Attribution is appreciated but not required.

---

## Contributing

Contributions are welcome. Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

---

## Citation

If you use HCLS AI Factory in academic work, please consider citing:

```
HCLS AI Factory: An End-to-End Precision Medicine Platform
https://github.com/ajones1923/hcls-ai-factory
```

---

## Acknowledgments

Built with technologies from:

- **NVIDIA** — Parabricks, BioNeMo, DGX Spark
- **Anthropic** — Claude AI
- **ClinVar** — Clinical variant database
- **AlphaMissense** — AI pathogenicity predictions

---

## About

HCLS AI Factory was created by **Adam Jones** as a vendor-neutral baseline platform, designed for organizations to adopt, extend, and optimize for their own infrastructure. It represents 14+ years of genomic research experience distilled into an accessible, end-to-end workflow.

The platform intentionally contains no vendor-specific dependencies beyond the core NVIDIA computing stack, allowing storage vendors, cloud providers, system integrators, and life sciences organizations to build tailored solutions on top of this foundation.

---

## Author

**Adam Jones**
Creator and Maintainer
[LinkedIn](https://www.linkedin.com/in/socal-engineer)

---

*Transform genomic data into therapeutic insights.*
