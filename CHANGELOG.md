# Changelog

All notable changes to the HCLS AI Factory will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-03-01

### Added

- `quickstart.sh` — One-command setup (prerequisites, venvs, deps, Milvus, all services)
- `make quickstart` Makefile target
- `data/demo/demo_variants.vcf` — 90 synthetic variants across 13 therapeutic areas for Stage 2 demos
- `docs/HCLS_AI_FACTORY_ARXIV_PAPER.md` — 10-section academic paper with 23 references
- `docs/GTC_2026_SUBMISSION.md` — NVIDIA GTC 2026 talk submission (50-min session outline)
- `docs/QUICKSTART_VIDEO_SCRIPT.md` — 5-minute video storyboard with 29 shots
- `docs/LAUNCH_ANNOUNCEMENT.md` — LinkedIn, Twitter/X, Hacker News, and Reddit launch drafts
- `CITATION.cff` — GitHub "Cite this repository" metadata
- README badges: Python 3.10+, MkDocs Docs, NVIDIA DGX Spark
- MkDocs "Publications" navigation section

### Changed

- All 3 agent demo guides rewritten to UI-driven flow (zero curl in live demo sections, curl preserved in API Reference appendix)

### Fixed

- CAR-T: double retrieval in `agent.py` — `run()` called `rag.retrieve()` then `rag.query()` which internally called `retrieve()` again
- CAR-T: `identified_topics` never populated in `search_plan()`
- CAR-T: stale docstrings ("5 collections" hardcoded instead of dynamic count)
- Imaging: `SyntheticCTResult` missing `volume_path`, `segmentation_mask_path`, `voxel_spacing_mm`, `model_name` fields
- Imaging: `knowledge_used` logic checked same boolean for all 3 categories
- Imaging: `body_region` filter accepted but never passed through to query
- Imaging: typo "solitiary" → "solitary"
- Imaging: deleted 140-line dead code `metrics.py` (never imported)
- Oncology: `target_collections` NameError in `retrieve()` method
- Oncology: monkey-patched attributes on `SearchHit` — added `label`, `citation`, `relevance` fields and `record_id` property
- Oncology: `variants` vs `key_variants` structure mismatch in `therapy_ranker.py`

### Security

- CORS origins restricted from wildcard to explicit localhost origins in all 3 agents
- All agents migrated from deprecated `class Config:` to `model_config = SettingsConfigDict(...)`
- Cleared hardcoded Orthanc default password in Imaging agent settings
- Imaging Dockerfile: removed unnecessary `curl` from builder stage

## [1.1.0] - 2026-02-28

### Added
- **CAR-T Intelligence Agent** — Cross-functional intelligence across the CAR-T cell therapy development lifecycle (10 collections, 6,266+ vectors, comparative analysis, deep research mode, PDF export, 241 tests)
- **Imaging Intelligence Agent** — AI-powered medical imaging analysis with NVIDIA NIM microservices: VISTA-3D, MAISI, VILA-M3, Llama-3 (10 collections, 4 workflow demos, FHIR R4 export, 539 tests)
- **Precision Oncology Agent** — Clinical decision support for molecular tumor boards (11 collections, case management, MTB packet generation, trial matching, therapy ranking, FHIR R4 bundle export, 516 tests)
- Cross-modal triggers connecting agents to shared genomic evidence (3.5M vectors)
- Agent Streamlit UIs on ports 8521 (CAR-T), 8525 (Imaging), 8526 (Precision Oncology)
- UI-driven demo guides for all three agents with pre-demo setup and live demo scripts
- Agent design documents and project bibles
- ROADMAP.md — Development trajectory and planned features
- TROUBLESHOOTING.md — Centralized troubleshooting guide for all components
- PERFORMANCE.md — Benchmark data for all pipeline stages and agents on DGX Spark
- HLS Orchestrator README.md — Pipeline orchestrator documentation
- Combined test suite: 1,296 tests across all agents (all passing)
- Agents integrated into platform infrastructure (landing page, docker-compose, CI/CD, startup scripts, .env.example)
- `pyproject.toml` with ruff, mypy, pytest, and coverage configuration
- `Makefile` with standard targets (lint, test, format, security, docs, clean)
- `.pre-commit-config.yaml` with ruff, gitleaks, and pre-commit-hooks
- `.gitleaks.toml` for automated secrets scanning in CI
- `.editorconfig` for cross-editor formatting consistency
- GitHub issue templates (bug report, feature request) with agent options
- Gitleaks secrets scanning and mypy type checking jobs in CI
- CODEOWNERS for automated PR review assignment

### Changed
- Codecov project coverage threshold raised from 2% to 50%
- CORS default restricted from wildcard to explicit localhost origins
- Milvus and monitoring ports bound to localhost only
- Docker Compose services given memory and CPU resource limits
- Streamlit version aligned across services (1.54.0)
- Pytest markers consolidated across all service pytest.ini files
- CI agent paths updated with graceful skip when directories not present
- Orphaned documentation pages added to MkDocs nav
- PR template updated with agent and CI/CD checkboxes

### Fixed
- Milvus filter expression injection prevention across all agents
- Async/sync method call correctness in FastAPI route handlers
- UTC ISO-8601 timestamp standardization (replacing naive datetime.utcnow)
- Missing RAG engine methods (cross_collection_search, search, synthesize)
- AgentQuery/AgentResponse model completeness for agent-RAG integration
- FHIR R4 export parameter handling
- Bare exception blocks replaced with specific error handling and logging

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
