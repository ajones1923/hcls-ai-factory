# Tuberous Sclerosis Complex — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/disease-programs/tuberous-sclerosis` · 115 Python files · 7,598 LOC · 24 test files

## Registered capabilities

_No capability is registered for this subject in `lib/hcls_common/capabilities.json`._

_No service port — this subject exposes no single registered endpoint._

## Principal modules

### `src/viz/usd_ascii.py`

`write_usda`, `write_mosaic_usda`, `write_atlas_usda`, `write_population_usda`

- **`write_mosaic_usda`** — Author the mosaic 'powers-of-ten' cell field as a USD PointInstancer: exactly the

### `api/main.py`

`lifespan`, `create_app`


### `src/eval/harness.py`

`evaluate`

- **`evaluate`** — Run the full scorecard. `featured` defaults to the engine's featured map if present.

### `src/agents/trajectory_modeler/runner.py`

`TrajectoryModeler`


### `src/orchestrator/graph.py`

`Orchestrator`



## Dependencies

`anthropic>=0.40`, `fastapi>=0.109`, `langgraph>=0.2`, `minio>=7.2`, `numpy>=1.26`, `psycopg2-binary>=2.9`, `pydantic-settings>=2.0`, `pydantic>=2.0`, `pymilvus>=2.4`, `pytest>=8.0`, `pyyaml>=6.0`, `redis>=5.0`, `sentence-transformers>=2.2`, `statsmodels>=0.14`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py tuberous-sclerosis
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Gene-therapy correction of TSC1/TSC2 is preclinical — an open design bench, not a cure today. Pediatric caution at full force. Decision support, never prescribing.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/disease-programs/tuberous-sclerosis`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
