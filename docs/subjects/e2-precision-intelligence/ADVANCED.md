# Precision Intelligence — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/precision-intelligence` · 30 Python files · 11,346 LOC · 11 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `precision-intelligence-engine` | engine | **live** | localhost:5001 |

**:5001** — a single-service portal, so there is no separate API port

## Principal modules

### `src/knowledge.py`

`get_gene_reference_data`, `get_knowledge_for_genes`, `get_knowledge_for_evidence`, `format_knowledge_for_prompt`, `get_druggable_genes`, `get_gene_drugs`

- **`get_gene_reference_data`** — Get reference data for a gene including UniProt ID, seed compound SMILES,
- **`get_knowledge_for_genes`** — Get knowledge connections for a list of genes.
- **`get_knowledge_for_evidence`** — Extract genes from evidence and return their knowledge connections.
- **`format_knowledge_for_prompt`** — Format knowledge connections as context for Claude's prompt.

### `app/chat_ui.py`

`get_variant_stats`, `load_shared_model`, `get_provider_for_model`, `get_rag_engine`, `get_target_manager`, `get_vcf_preview`

- **`get_variant_stats`** — Get variant counts from Milvus database.
- **`load_shared_model`** — Load model and provider from shared file (set by portal).
- **`get_provider_for_model`** — Get the provider for a given model.
- **`get_rag_engine`** — Initialize RAG engine with specified model.

### `portal/app/server.py`

`RateLimiter`, `require_api_key`, `load_config`, `save_config`, `check_file_exists`, `get_file_size`

- **`RateLimiter`** — Simple in-memory rate limiter.
- **`require_api_key`** — Require X-API-Key header for dangerous endpoints.
- **`load_config`** — Load pipeline configuration from .env
- **`save_config`** — Save configuration to .env

### `src/rag_engine.py`

`RAGEngine`, `create_rag_engine`

- **`RAGEngine`** — RAG Engine for genomic evidence retrieval and question answering.
- **`create_rag_engine`** — Factory function to create a fully configured RAG engine.

### `src/annotator.py`

`AnnotationResult`, `VariantAnnotator`, `LocalVEPAnnotator`, `ClinVarAnnotator`, `AlphaMissenseAnnotator`

- **`AnnotationResult`** — Result of variant annotation.
- **`VariantAnnotator`** — Annotate variants with gene names, consequences, and clinical information.
- **`LocalVEPAnnotator`** — Annotate using local VEP installation (Docker or native).
- **`ClinVarAnnotator`** — Fast local ClinVar annotation using pre-downloaded variant_summary.txt.gz.


## Dependencies

`anthropic==0.75.0`, `cyvcf2==0.31.4`, `fastapi==0.128.0`, `flask-cors==6.0.2`, `flask==3.1.2`, `loguru==0.7.3`, `numpy==2.4.0`, `openai==2.15.0`, `opentelemetry-api>=1.29.0`, `opentelemetry-sdk>=1.29.0`, `pandas==2.3.3`, `psutil==7.2.1`, `pydantic-settings==2.12.0`, `pydantic==2.12.5`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py precision-intelligence
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Retrieval quality is bounded by what has been indexed. An unseeded collection returns nothing — that is a data problem, not a reasoning failure.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/precision-intelligence`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
