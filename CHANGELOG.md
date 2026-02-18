# Changelog

All notable changes to the HCLS AI Factory will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.3] - 2026-02-18

### Added
- Cloud NIM API support for MolMIM and DiffDock (health.api.nvidia.com)
- 3-tier NIM strategy: cloud → local → mock fallback
- NIM_MODE environment variable (cloud/local) with auto-detection
- CloudMolMIMClient and CloudDiffDockClient in nim_clients.py
- NVCF asset staging for DiffDock cloud docking
- Batched MolMIM generation (CMA-ES popsize=20 limit handling)
- Real NIM inference report generator (run_cloud_nim_report.py)
- Regenerated VCP Drug Candidate Report with real NIM data (27 molecules docked)

### Fixed
- ARM64/DGX Spark compatibility — cloud NIM bypasses x86-only container limitation
- DiffDock response parsing for array-format position_confidence
- Path(protein_pdb).exists() OSError on PDB content strings
- Docker image names: nvcr.io/nim/nvidia/molmim:1.0.0 and nvcr.io/nim/mit/diffdock:2.2.0

## [1.0.2] - 2026-02-17

### Changed
- Project relocated from /transfer/ to /projects/hcls-ai-factory/
- Updated all internal path references

## [1.0.1] - 2026-02-14

### Added
- Landing page Dockerfile and containerization
- Web portal test suite
- Genomics pipeline quickstart guide improvements

## [1.0.0] - 2026-02-10

### Added
- End-to-end 10-stage drug discovery pipeline (Initialize -> Reporting)
- RAG/Chat pipeline with ClinVar (4.1M variants) and AlphaMissense (71M predictions)
- Genomics pipeline with Parabricks 4.6 GPU-accelerated analysis
- Pipeline checkpoint/resume capability for fault tolerance
- Grafana dashboards for pipeline health and GPU monitoring
- Prometheus alerting rules for service availability and GPU health
- Docker Compose orchestration for all services
- OpenAPI specifications for Genomics and RAG APIs
- Milvus vector database integration with BGE-small-en-v1.5 embeddings
- Streamlit UIs for Chat (port 8501) and Drug Discovery (port 8505)
- Flask landing page with service health dashboard (port 8080)
- Nextflow DSL2 orchestrator for pipeline coordination
- Environment preflight check script
- Comprehensive test suite (255+ tests across all services)
- CI/CD with lint, test, security scanning, and docs build
- Codecov coverage reporting with per-service flags
- Multi-stage Docker builds for all services
- MkDocs Material documentation site
- VCP/FTD demo pipeline with CB-5083 seed compound
- Performance benchmarking with per-stage timing
- Apache 2.0 license
