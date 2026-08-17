# Precision Biomarker — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/precision-biomarker` · 58 Python files · 29,398 LOC · 20 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `precision-biomarker-agent` | agent | **live** | localhost:8528 |

**UI :8528 · API :8529** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `scripts/expand_biomarker_reference.py`

`main`


### `app/biomarker_ui.py`

`init_engine`, `risk_badge`, `get_pgx_phenotype`

- **`init_engine`** — Initialize the Biomarker Intelligence analysis engine (cached across reruns).
- **`risk_badge`** — Return an HTML span with risk-level coloring.
- **`get_pgx_phenotype`** — Determine metabolizer phenotype from star alleles (simplified).

### `src/pharmacogenomics.py`

`PharmacogenomicMapper`

- **`PharmacogenomicMapper`** — Maps star alleles and genotypes to drug recommendations.

### `src/disease_trajectory.py`

`DiseaseTrajectoryAnalyzer`

- **`DiseaseTrajectoryAnalyzer`** — Detects pre-symptomatic disease trajectories using genotype-stratified

### `src/knowledge.py`

`get_domain_context`, `get_pgx_context`, `get_biomarker_context`, `get_cross_modal_context`, `get_knowledge_stats`

- **`get_domain_context`** — Return formatted knowledge context for a disease domain.
- **`get_pgx_context`** — Return formatted PGx knowledge context for a pharmacogene.
- **`get_biomarker_context`** — Return formatted knowledge context for a specific biomarker.
- **`get_cross_modal_context`** — Return formatted cross-modal link context.


## Dependencies

`anthropic>=0.18.0,<1.0`, `fastapi>=0.109.0,<1.0`, `loguru>=0.7.0,<1.0`, `numpy>=1.24.0,<3.0`, `pandas>=2.0.0,<3.0`, `plotly>=5.18.0,<6.0`, `prometheus-client>=0.20.0,<1.0`, `pydantic-settings>=2.7,<3.0`, `pydantic>=2.0,<3.0`, `pymilvus>=2.4.0,<2.6`, `pytest-asyncio>=0.21,<1.0`, `pytest-cov>=4.0,<5.0`, `pytest>=7.0,<8.0`, `python-dotenv>=1.0.0,<2.0`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py precision-biomarker
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Many biomarker inputs are research- or trial-use, not validated for routine clinical practice.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/precision-biomarker`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
