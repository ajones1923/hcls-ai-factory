# Single-Cell Intelligence — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/single-cell` · 47 Python files · 20,940 LOC · 14 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `single-cell-intelligence-agent` | agent | **live** | localhost:8540 |

**UI :8540 · API :8541** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/agent.py`

`SCWorkflowType`, `EvidenceLevel`, `SeverityLevel`, `AnalysisModality`, `CellOntologyDomain`, `SCResponse`

- **`SCWorkflowType`** — Types of single-cell analysis workflows.
- **`EvidenceLevel`** — Clinical evidence hierarchy for single-cell findings.
- **`SeverityLevel`** — Finding severity classification.
- **`AnalysisModality`** — Single-cell analysis modalities.

### `src/rag_engine.py`

`SCSearchResult`, `get_all_collection_names`, `SingleCellRAGEngine`

- **`SCSearchResult`** — A single search result from a Milvus collection.
- **`get_all_collection_names`** — Return all collection names.
- **`SingleCellRAGEngine`** — Multi-collection RAG engine for single-cell intelligence.

### `src/knowledge.py`

_no public symbols_


### `src/clinical_workflows.py`

`BaseSCWorkflow`, `CellTypeAnnotationWorkflow`, `TMEProfilingWorkflow`, `DrugResponseWorkflow`, `SubclonalArchitectureWorkflow`, `SpatialNicheWorkflow`

- **`BaseSCWorkflow`** — Abstract base for all single-cell clinical workflows.
- **`CellTypeAnnotationWorkflow`** — Multi-strategy consensus cell type annotation.
- **`TMEProfilingWorkflow`** — Tumor microenvironment profiling: classify hot/cold/excluded/
- **`DrugResponseWorkflow`** — Cell-type-specific drug sensitivity from DepMap/CCLE data.

### `api/routes/sc_clinical.py`

`integrated_assessment`, `QueryRequest`, `QueryResponse`, `SearchRequest`, `SearchResult`, `SearchResponse`

- **`integrated_assessment`** — Multi-agent integrated assessment combining insights from across the HCLS AI Factory.
- **`QueryRequest`** — Free-text RAG query with optional domain and patient context.
- **`SearchRequest`** — Multi-collection semantic search.


## Dependencies

`anthropic==0.25.0`, `apscheduler==3.10.4`, `fastapi==0.111.0`, `httpx==0.27.0`, `loguru==0.7.2`, `lxml==5.2.1`, `numpy==1.26.4`, `prometheus-client==0.20.0`, `pydantic-settings==2.2.1`, `pydantic==2.7.4`, `pymilvus==2.4.1`, `python-docx==1.1.0`, `python-dotenv==1.0.1`, `python-multipart==0.0.9`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py single-cell
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Interpretation rests on the quality of the upstream clustering and the marker panel used.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/single-cell`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
