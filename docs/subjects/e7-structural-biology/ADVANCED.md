# Structural Biology — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/structural-biology` · 21 Python files · 1,233 LOC · 7 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `esmfold-model` | model | **live** | localhost:8570 |
| `esm2-search` | model | **live** | localhost:8571 |
| `protein-developability` | model | **planned** | localhost:8576 |
| `mhcflurry-immunogenicity` | model | **planned** | localhost:8577 |
| `proteinmpnn-design` | model | **live** | localhost:8578 |
| `esm2-finetune` | model | **planned** | — |
| `chai1-structure` | model | **planned** | — |

_No service port — this subject exposes no single registered endpoint._

## Principal modules

### `src/esmfold_service.py`

`SequenceError`, `validate_sequence`, `parse_pdb_stats`, `FoldResult`, `ESMFoldModel`, `create_app`

- **`SequenceError`** — Raised for an invalid protein sequence.
- **`validate_sequence`** — Normalize + validate a one-letter amino-acid sequence.
- **`parse_pdb_stats`** — Extract structure stats from a PDB string. ESMFold writes per-residue pLDDT

### `src/protein_search.py`

`VectorBackend`, `InMemoryBackend`, `MilvusBackend`, `ProteinSearchIndex`

- **`InMemoryBackend`** — Exact cosine search over normalized vectors (dot product). Test/fallback backend.
- **`MilvusBackend`** — Real Milvus-backed protein index (reuses the platform vector DB).
- **`ProteinSearchIndex`** — Embed proteins and search them by sequence similarity.

### `src/esm2_finetune.py`

`ESM2FineTuner`


### `src/client.py`

`ESMFoldClient`, `MockESMFoldClient`, `create_esmfold_client`, `ProteinSearchClient`

- **`ESMFoldClient`** — Calls a running ESMFold service (real model).
- **`MockESMFoldClient`** — Offline placeholder — clearly labeled, NEVER a real structure.
- **`create_esmfold_client`** — Factory mirroring NIMServiceManager: ESMFOLD_MODE = local | mock (default local).
- **`ProteinSearchClient`** — Calls a running protein embedding + search service (B2).

### `src/proteinmpnn_design.py`

`design`

- **`design`** — Run ProteinMPNN on a backbone (PDB text or a path) → designed sequences.


## Dependencies

`accelerate>=1.0`, `einops>=0.8`, `fastapi>=0.110`, `httpx>=0.27`, `pydantic>=2.0`, `torch>=2.4`, `transformers>=5.0`, `uvicorn>=0.30`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py structural-biology
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

ESMFold prediction and ProteinMPNN design need CUDA; this box runs a CPU-only PyTorch build. Structures shown are deposited PDB entries, not predictions made live.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/structural-biology`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
