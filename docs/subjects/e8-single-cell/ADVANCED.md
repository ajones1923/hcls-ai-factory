# Single-Cell Analysis — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/single-cell` · 3 Python files · 174 LOC · 1 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `singlecell-compute` | engine | **live** | localhost:8573 |

**:8573** — a single-service portal, so there is no separate API port

## Principal modules

### `src/single_cell_compute.py`

`annotate_cluster`, `summarize_clusters`, `SingleCellAnalysis`

- **`annotate_cluster`** — Annotate a cluster by overlap of its top DE genes with each cell type's markers.
- **`summarize_clusters`** — clusters: [{"cluster": id, "n_cells": int, "top_genes": [..]}] -> annotated summary.

### `src/single_cell_service.py`

`create_app`, `SingleCellClient`



## Dependencies

`anndata>=0.10,<1.0`, `fastapi>=0.110,<1.0`, `httpx>=0.26,<1.0`, `igraph>=0.11,<1.0`, `leidenalg>=0.10,<1.0`, `numpy>=1.24,<3.0`, `pandas>=2.0,<3.0`, `pydantic>=2.0,<3.0`, `pynndescent>=0.5,<1.0`, `scanpy>=1.10,<2.0`, `scipy>=1.11,<2.0`, `uvicorn[standard]>=0.27,<1.0`

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

Cluster annotation is by marker-gene overlap against a canonical panel. It is a strong heuristic, not a ground truth.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/single-cell`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
