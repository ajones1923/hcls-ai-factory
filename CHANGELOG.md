# Changelog

All notable changes to the HCLS AI Factory will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2025-02-10

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
- Comprehensive test suite (251+ tests across all services)
- CI/CD with lint, test, security scanning, and docs build
- Codecov coverage reporting with per-service flags
- Multi-stage Docker builds for all services
- MkDocs Material documentation site
- VCP/FTD demo pipeline with CB-5083 seed compound
- Performance benchmarking with per-stage timing
- Apache 2.0 license
