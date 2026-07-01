# HCLS AI Factory — Comprehensive Improvement Recommendations

**Date:** 2026-02-27
**Scope:** End-to-end audit of all components across the HCLS AI Factory
**Method:** Deep code review of every source file across all pipelines, orchestrator, shared library, landing page, infrastructure, tests, CI/CD, Docker configs, and documentation

---

## Executive Summary

The HCLS AI Factory is a mature, well-architected precision medicine platform with strong fundamentals. The audit identified **~130 findings** across all components. The breakdown by severity:

| Severity | Count | Key Themes |
|----------|-------|------------|
| **CRITICAL** | 4 | API key exposure (2 repos), undefined Nextflow channels, silent download failures |
| **HIGH** | 25 | Security gaps, missing tests for critical paths, docking score inconsistencies, no parallelism in docking |
| **MEDIUM** | 55 | Code duplication, RAG quality improvements, Docker inconsistencies, dependency drift |
| **LOW** | ~45 | Style issues, documentation mismatches, minor optimizations |

**Overall Rating: 8.5/10** — Production-quality architecture with specific areas needing hardening before enterprise deployment.

---

## Priority 1: CRITICAL — Fix Immediately

### 1.1 Rotate Exposed API Keys
- **Files:** `rag-chat-pipeline/.env` (Anthropic key), `drug-discovery-pipeline/.env` (NVIDIA key)
- **Risk:** Keys are in plaintext in working directories. If synced to any repo, they're compromised.
- **Action:**
  1. Rotate both API keys immediately (Anthropic console + NVIDIA NGC)
  2. Run `git rm --cached .env` in any repo that tracked them
  3. Verify `.gitignore` includes `.env` in every pipeline directory
  4. Consider using a secrets manager or environment-only injection

### 1.2 Fix Undefined Nextflow Channels
- **File:** `hls-orchestrator/main.nf`, line 75
- **Risk:** Adding new modes or edge cases causes runtime crash (`No such variable`)
- **Action:** Pre-declare all channels at the top of the workflow block:
  ```groovy
  ch_versions  = Channel.empty()
  ch_molecules = Channel.empty()
  ch_targets   = Channel.empty()
  ```

### 1.3 Fix Silent Download Failures in Genomics
- **File:** `genomics-pipeline/scripts/02-download-data.sh`, lines 110, 130
- **Risk:** `|| true` swallows aria2c failures; corrupt merged FASTQ produced silently
- **Action:** Remove `|| true`, add per-file checksum verification before merge, fail fast on download errors

### 1.4 Fix `genomics_only` Mode — Declared But Never Handled
- **File:** `hls-orchestrator/main.nf` and `nextflow.config`
- **Risk:** Using `--mode genomics_only` falls through all conditionals, producing empty outputs
- **Action:** Add an `else if (params.mode == 'genomics_only')` block in `main.nf`

---

## Priority 2: HIGH — Fix Before Next Release

### 2.1 Security Hardening

| # | Issue | File | Action |
|---|-------|------|--------|
| A | Chat UI filter injection | `rag-chat-pipeline/app/chat_ui.py:981-997` | Sanitize gene/chrom inputs using `MilvusClient._sanitize_gene()` before constructing Milvus filter expressions |
| B | Portal executes arbitrary shell commands | `rag-chat-pipeline/portal/app/server.py:322-366` | Require `PORTAL_API_KEY` (fail closed when unset) or use parameterized commands |
| C | File upload/delete has no auth | `genomics-pipeline/web-portal/app/server.py:760,804` | Add `@require_api_key` to upload and delete endpoints |
| D | Security headers defined but never applied | `genomics-pipeline/web-portal/app/security.py` | Add `app.after_request(add_security_headers)` in `server.py` |
| E | `stop_all` kills arbitrary host processes | `genomics-pipeline/web-portal/app/server.py:463-499` | Track child PIDs explicitly instead of pattern-matching all system processes |
| F | `allow_pickle=True` without integrity check | `rag-chat-pipeline/src/embedder.py:133` | Add HMAC checksum to cache files or use pickle-free serialization |
| G | Grafana default credentials in .env | `drug-discovery-pipeline/.env:24-25` | Change from `admin/admin` and don't store in version control |
| H | `unsafe_allow_html=True` throughout Chat UI | `rag-chat-pipeline/app/chat_ui.py` | HTML-escape all dynamic values with `html.escape()` before injection |

### 2.2 Fix Milvus Backup/Restore (Will Fail)
- **File:** `rag-chat-pipeline/scripts/milvus_backup.py`
- **Issues:**
  1. `expr="id >= 0"` but `id` is VARCHAR — invalid query (line 59)
  2. Restore has wrong field order and omits `id` primary key (lines 121-138)
- **Action:** Fix query to `expr='id != ""'`; restore data must match schema order and include `id`

### 2.3 Docking Score Normalization Inconsistency
- **Files:** `drug-discovery-pipeline/src/pipeline.py:514-518` vs `run_cloud_nim_report.py:209-214`
- **Issue:** Pipeline uses `-dock/12.0`, cloud report uses `-dock/6.0` — non-comparable results
- **Action:** Create a shared `normalize_docking_score()` utility; document the actual DiffDock score range

### 2.4 Add Parallel Docking (Major Performance Win)
- **File:** `drug-discovery-pipeline/src/pipeline.py:461-489`
- **Issue:** Sequential docking of 50 molecules takes 4-25 minutes
- **Action:** Use `concurrent.futures.ThreadPoolExecutor(max_workers=4)` or `dock_batch()` method

### 2.5 Cache Cloud DiffDock Protein Asset
- **File:** `drug-discovery-pipeline/src/nim_clients.py:494-503`
- **Issue:** Protein PDB re-uploaded for every molecule (100+ NVCF API calls for 50 molecules)
- **Action:** Upload protein once, cache asset ID, reuse for all docking calls in the run

### 2.6 Add Tests for Critical Untested Code

| Module | Location | Priority |
|--------|----------|----------|
| `hcls_common` library (14 modules, 0 tests) | `lib/hcls_common/` | **Highest** — shared across all pipelines |
| Cloud NIM clients | `drug-discovery-pipeline/src/nim_clients.py` | **High** — production code path for DGX Spark |
| Annotator module | `rag-chat-pipeline/src/annotator.py` | **High** — ClinVar/AlphaMissense parsing |
| Embedder module | `rag-chat-pipeline/src/embedder.py` | **High** — cache hit/miss logic |
| CryoEM evidence | `drug-discovery-pipeline/src/cryoem_evidence.py` | **Medium** |
| Shell scripts (18 scripts, 0 tests) | `genomics-pipeline/scripts/` | **Medium** — add BATS or shellcheck CI |

### 2.7 Version Mismatch in Nextflow Pipeline
- **File:** `hls-orchestrator/main.nf:17` says `1.0.0`, `nextflow.config:164` says `1.0.3`
- **Action:** Update `main.nf` to `1.0.3`; create a single `VERSION` file read everywhere

### 2.8 Add Schema Validation at Stage Boundaries
- **File:** `drug-discovery-pipeline/src/target_import.py:57-83`
- **Issue:** Silently accepts malformed target JSON with defaults for all fields
- **Action:** Use Pydantic model validation; make `gene` required (not defaulting to "Unknown")

---

## Priority 3: MEDIUM — Fix Soon (Next Sprint)

### 3.1 RAG Quality Improvements

| # | Issue | Action |
|---|-------|--------|
| A | Query expansion uses substring match ("false" matches "als") | Use word boundary regex: `re.search(r'\b' + re.escape(keyword) + r'\b', query_lower)` |
| B | No token budget management — prompts can exceed context window | Estimate tokens (~4 chars/token), truncate evidence if exceeding 80% of context |
| C | No conversation history — each query is independent | Add `messages` parameter to LLM `generate()`; pass recent chat history |
| D | Hardcoded retrieval score (0.85) for gene expansion results | Lower to 0.60 or derive from query relevance |
| E | No LLM retry logic — single failure crashes query | Add `tenacity` retry with exponential backoff for rate limits, timeouts |
| F | `to_dict()` missing `disease_associations` field | Add to `VariantEvidence.to_dict()` |

### 3.2 Code Quality Improvements

| # | Issue | Action |
|---|-------|--------|
| A | Massive shell script duplication (7 scripts with identical helpers) | Extract into `scripts/common.sh` and `source` from each |
| B | Duplicate `GeneratedMolecule` class (Pydantic vs dataclass) | Consolidate to single Pydantic model |
| C | `validation.py` never imported in genomics `server.py` | Import and use the validation module |
| D | `sys.path.insert` everywhere instead of proper packaging | Add `pyproject.toml` with `pip install -e .` for each pipeline |
| E | Pydantic v1 `.dict()` mixed with v2 `.model_dump()` | Replace all `.dict()` with `.model_dump()` |
| F | Hardcoded path in Drug Discovery UI | `discovery_ui.py:22` — use env var or relative path |
| G | UI uses fallback `MoleculeGenerator` instead of NIM clients | Use `NIMServiceManager` in UI for real molecule generation |
| H | Multiple competing download/merge scripts (7 scripts) | Consolidate into one with flags; archive unused |
| I | OllamaClient and VLLMClient are near-identical | Extract `OpenAICompatibleClient` base class |

### 3.3 Dependency Management

| # | Issue | Action |
|---|-------|--------|
| A | numpy version divergence (2.4.0 / 2.4.1 / 2.4.2) | Pin to single version across all pipelines |
| B | Working repo vs public repo version drift | Sync working dirs with `hcls-ai-factory-public/` |
| C | No dependency lock files | Add `pip-compile` workflow or constraints files |
| D | Orchestrator portal uses unpinned `>=` ranges | Pin to exact versions |
| E | Missing opentelemetry/reportlab in public repo | Sync with working dir |

### 3.4 Docker & Infrastructure

| # | Issue | Action |
|---|-------|--------|
| A | Inconsistent SHA pinning (2 of 4 Dockerfiles) | Pin all base images with SHA256 digest |
| B | Landing page `COPY . .` copies everything | Add `.dockerignore` or selective COPY |
| C | Both NIM containers claim 1 GPU (conflicts on single-GPU) | Document cloud mode for DGX Spark; add GPU sharing config |
| D | NVCF assets never cleaned up | Add cleanup method after docking completes |
| E | CORS allows `*` by default in genomics portal | Default to localhost; require explicit opt-in |
| F | Remove deprecated `version: "3.8"` from docker-compose files | Delete the line |

### 3.5 Performance Optimizations

| # | Issue | Impact | Action |
|---|-------|--------|--------|
| A | AlphaMissense loads 71M records into RAM (~10GB) | Memory pressure | Consider SQLite or memory-mapped file |
| B | `get_cpu_utilization()` blocks for 0.2s per request | API latency | Use `psutil.cpu_percent(interval=None)` |
| C | NVML init/shutdown on every metrics request | API latency | Init once at startup |
| D | Genomics download decompresses/recompresses 200GB unnecessarily | 30-60 min wasted | Use `cat` for gzip concatenation |
| E | BED file gene lookup is O(n) linear scan | Ingestion speed | Use `bisect` binary search or interval tree |
| F | Embedding batch size of 32 is conservative | Ingestion speed | Increase to 128-256 for BGE-small |
| G | `get_variant_stats` makes 9 separate Milvus queries | Page load speed | Single query with `IN` expression |

### 3.6 Composite Score Validation
- **File:** `drug-discovery-pipeline/src/models.py`
- **Action:** Add Pydantic `@model_validator` ensuring `docking_weight + generation_weight + qed_weight == 1.0`

### 3.7 CI/CD Improvements
- Add pre-commit hooks (ruff, mypy)
- Add root-level `pyproject.toml` with `[tool.ruff]` settings
- Add `hcls_common` to CI matrix
- Add CAR-T agent to CI matrix
- Update OpenAPI specs from 1.0.0 to 1.0.3

---

## Priority 4: LOW — Address When Convenient

### 4.1 Minor Code Issues
- Remove debug `echo "SCRIPT STARTED"` messages from shell scripts
- Add lock file (`flock`) to prevent concurrent `run.sh` executions
- Fix step numbering in `03-setup-reference.sh` (STEP X/4 vs X/5)
- Add cleanup trap in `03-setup-reference.sh` for heartbeat processes
- Fix `get_file_type()` to handle compound extensions (`.vcf.gz`)
- Cap portal log list to prevent unbounded memory growth
- Clean up NVCF checkpoint files (add rotation/cleanup)
- Remove unused `llm_models` volume from rag-chat docker-compose
- Use `hmac.compare_digest` for API key comparison (timing attack prevention)
- Add Milvus port (19530) to README port table
- Use `loguru` consistently (not `print()` in `target_import.py`)

### 4.2 Documentation Updates
- Re-run test suite and update audit report counts (356 → 548+)
- Clarify port 8510 documentation (which service actually runs there)
- Document Prometheus host port (9099 vs container 9090)
- Document AlphaMissense memory requirements (~10GB)
- Document single-sample VCF assumption in parser

---

## Improvement Roadmap for Future Development

### Phase 1: Foundation Hardening (Before New Agent Development)
1. Rotate API keys and fix `.gitignore`
2. Fix Nextflow orchestrator bugs (undefined channels, genomics_only mode)
3. Add tests for `hcls_common` library
4. Sync dependency versions across all pipelines
5. Add pre-commit hooks and linting config

### Phase 2: Performance & Quality (Next 2-4 Weeks)
1. Add parallel docking in drug discovery
2. Cache cloud DiffDock protein assets
3. Improve RAG query expansion (word boundaries, token budget)
4. Add conversation history to LLM calls
5. Add LLM retry logic with exponential backoff
6. Consolidate duplicate code (shell helpers, GeneratedMolecule, LLM clients)

### Phase 3: Enterprise Readiness (Ongoing)
1. Schema validation at all stage boundaries
2. Complete test coverage for cloud NIM clients
3. Dependency lock files
4. Consistent Docker SHA pinning
5. Security hardening (CORS, auth, HTML escaping)
6. Add integration tests with real Milvus

### Phase 4: Agent Ecosystem Preparation
When adding new agents (CAR-T already production-ready, Imaging and Biomarker in progress):
1. Use `hcls_common` library for all shared functionality (once tested)
2. Follow the established patterns: Streamlit UI + Flask API + Docker + tests
3. Add new agents to CI matrix
4. Maintain consistent dependency versions
5. Use the query router in `hcls_common` for cross-collection searches
6. Register new agents in the landing page service health monitor
7. Add Nextflow modules for agent-specific pipeline stages

---

## Summary: Top 10 Highest-Impact Actions

| # | Action | Impact | Effort |
|---|--------|--------|--------|
| 1 | Rotate exposed API keys | Eliminates security risk | 15 min |
| 2 | Fix Nextflow undefined channels | Prevents runtime crashes | 30 min |
| 3 | Add parallel docking | 3-5x speed improvement in Stage 3 | 2 hours |
| 4 | Cache cloud protein assets | 50% fewer API calls in docking | 1 hour |
| 5 | Fix query expansion substring matching | Eliminates false positive gene expansions | 30 min |
| 6 | Add LLM retry logic | Prevents query failures from transient errors | 1 hour |
| 7 | Add tests for `hcls_common` | De-risks shared library used by all pipelines | 1 day |
| 8 | Sync dependency versions | Prevents cross-pipeline ABI issues | 2 hours |
| 9 | Add conversation history to RAG | Enables multi-turn clinical reasoning | 3 hours |
| 10 | Consolidate shell script helpers | Eliminates 500+ lines of duplication | 2 hours |

---

*This document should be treated as a living reference. Check off items as they're completed and re-audit quarterly.*
