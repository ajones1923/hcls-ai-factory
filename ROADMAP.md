# Roadmap

This document outlines the development trajectory of the HCLS AI Factory, from its initial release through planned expansions in clinical intelligence, regulatory readiness, and multi-institutional deployment.

---

## Completed

### v1.0.0 — Foundation Release (February 2026)

The initial release establishes the complete end-to-end pipeline: Patient DNA to Drug Candidates on a single workstation.

- **Stage 1: GPU Genomics** — FASTQ to VCF via NVIDIA Parabricks 4.6 (120-240 min vs. 24-48 hrs on CPU)
- **Stage 2: Evidence RAG** — VCF to Target Hypothesis via Milvus + ClinVar + AlphaMissense + Claude (3.56M vectors, 201 genes, 13 therapeutic areas)
- **Stage 3: Drug Discovery** — Target to Ranked Drug Candidates via BioNeMo MolMIM + DiffDock + RDKit
- **Nextflow DSL2 Orchestrator** — 5 execution modes (full, target, drug, demo, genomics_only)
- **Landing Page** — Flask service health dashboard monitoring 12 services
- **Monitoring** — Grafana dashboards + Prometheus alerting + DCGM GPU metrics
- **Documentation** — MkDocs Material site, deployment guide (2,690 lines), project bible, learning guides
- **CI/CD** — GitHub Actions (lint, test, security scan, docs build), Codecov, Dependabot
- **Demo Pipeline** — VCP (p97) for Frontotemporal Dementia with CB-5083 seed compound

### v1.0.1–v1.0.3 — Hardening (February 2026)

- Cloud NIM API support (health.api.nvidia.com) with 3-tier fallback: cloud → local → mock
- ARM64/DGX Spark compatibility via cloud NIM bypass
- Landing page containerization and test suite
- Real NIM inference report with 27 docked molecules

---

## Current: Intelligence Agents (v1.1.0)

Three domain-specific intelligence agents extending the core platform with specialized RAG, cross-modal analysis, and clinical decision support.

### CAR-T Intelligence Agent

Cross-functional intelligence across the CAR-T cell therapy development lifecycle.

| Capability | Detail |
|---|---|
| Knowledge Base | 10 collections, 6,266+ vectors (literature, trials, constructs, manufacturing, safety, biomarkers) |
| Shared Data | Read-only access to Stage 2 genomic evidence (3.5M vectors) |
| Comparative Analysis | Auto-detected entity comparison with structured tables |
| Deep Research Mode | Multi-hop reasoning across all collections |
| Export | PDF report generation |
| Test Coverage | 241 tests |
| Port | 8521 (API), 8521 (Streamlit UI) |

### Imaging Intelligence Agent

AI-powered medical imaging analysis with NVIDIA NIM microservices.

| Capability | Detail |
|---|---|
| NIM Services | VISTA-3D (segmentation), MAISI (synthesis), VILA-M3 (VLM), Llama-3 (reasoning) |
| Workflow Demos | CT Head hemorrhage, CT Chest lung nodules, CXR triage, MRI brain MS lesions |
| Knowledge Base | 10 collections, cross-modal triggers |
| FHIR R4 | DiagnosticReport export |
| Test Coverage | 539 tests |
| Port | 8525 (API), 8525 (Streamlit UI) |

### Precision Oncology Agent

Clinical decision support for molecular tumor boards.

| Capability | Detail |
|---|---|
| Case Management | VCF ingest, variant annotation, case lifecycle |
| MTB Packet Generation | Evidence synthesis, therapy ranking, trial matching, open questions |
| Trial Matching | Biomarker-based clinical trial matching with eligibility scoring |
| Therapy Ranking | Evidence-tiered recommendations with resistance warnings |
| Knowledge Base | 11 collections (variants, therapies, resistance, pathways, biomarkers, trials, guidelines, outcomes, imaging, case history, shared genomic evidence) |
| FHIR R4 | DiagnosticReport bundle export with SNOMED CT / LOINC coding |
| Test Coverage | 516 tests |
| Port | 8526 (API), 8526 (Streamlit UI) |

**Combined Agent Test Coverage: 1,296 tests, all passing.**

---

## Planned

### Multi-Institutional Deployment

- Kubernetes Helm charts for production cluster deployment
- Singularity/Apptainer profiles for HPC environments (SLURM)
- Cloud deployment guides (AWS, GCP, Azure) with GPU instance recommendations
- Multi-tenant isolation with namespace-based resource partitioning

### Data Expansion

- COSMIC (Catalogue of Somatic Mutations in Cancer) integration
- PharmGKB pharmacogenomics knowledge base
- ChEMBL bioactivity data for drug-target interactions
- OncoKB precision oncology knowledge base
- CIViC (Clinical Interpretation of Variants in Cancer)
- gnomAD population frequency annotations

### Clinical Integration

- HL7 FHIR R4 interoperability for EHR integration
- DICOM ingestion pipeline for medical imaging
- Specimen tracking and biobank integration
- Clinical report templates for regulatory submission

### Platform Capabilities

- Multi-GPU distributed inference for large cohort processing
- Federated learning support for multi-site collaboration without data sharing
- Automated variant classification per ACMG/AMP guidelines
- Longitudinal patient tracking and outcome correlation
- Cohort analysis tools for population-level insights

### Observability and Governance

- OpenTelemetry distributed tracing across all pipeline stages
- Audit logging with tamper-evident storage
- Role-based access control (RBAC) framework
- Data lineage tracking from FASTQ through drug candidates
- Model versioning and reproducibility guarantees

---

## Contributing

We welcome contributions at any stage. See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines. If you are interested in contributing to a specific roadmap item, open an issue to coordinate.

---

## Versioning

This project follows [Semantic Versioning](https://semver.org/). See [CHANGELOG.md](CHANGELOG.md) for release history.
