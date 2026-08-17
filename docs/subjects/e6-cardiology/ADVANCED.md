# Cardiology — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/cardiology` · 59 Python files · 45,345 LOC · 20 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `cardiology-intelligence-agent` | engine | **live** | localhost:8126 |

**UI :8126 · API :8127** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/clinical_workflows.py`

`BaseCardioWorkflow`, `CADAssessmentWorkflow`, `HeartFailureWorkflow`, `ValvularDiseaseWorkflow`, `ArrhythmiaWorkflow`, `CardiacMRIWorkflow`

- **`BaseCardioWorkflow`** — Abstract base for all cardiology clinical workflows.
- **`CADAssessmentWorkflow`** — Coronary artery disease assessment integrating calcium scoring,
- **`HeartFailureWorkflow`** — Heart failure classification, GDMT assessment, and device
- **`ValvularDiseaseWorkflow`** — Valvular heart disease severity grading and intervention criteria

### `src/knowledge.py`

_no public symbols_


### `src/gdmt_optimizer.py`

`GDMTOptimizer`

- **`GDMTOptimizer`** — GDMT optimization engine based on ACC/AHA 2022 HF Guidelines.

### `src/query_expansion.py`

`QueryExpander`

- **`QueryExpander`** — Expand cardiology queries with related clinical terms.

### `src/risk_calculators.py`

`RiskCalculatorError`, `calculate_ascvd`, `calculate_heart_score`, `calculate_cha2ds2_vasc`, `calculate_has_bled`, `calculate_maggic`

- **`RiskCalculatorError`** — Raised when input validation fails for a risk calculator.
- **`calculate_ascvd`** — Calculate 10-year ASCVD risk using the Pooled Cohort Equations.
- **`calculate_heart_score`** — Calculate the HEART Score for emergency department chest pain triage.
- **`calculate_cha2ds2_vasc`** — Calculate CHA2DS2-VASc score for stroke risk in atrial fibrillation.


## Dependencies

`anthropic>=0.18.0`, `apscheduler>=3.10.0`, `biopython>=1.83`, `fastapi>=0.109.0`, `loguru>=0.7.0`, `lxml>=5.0.0`, `numpy>=1.24.0`, `opentelemetry-api>=1.29.0`, `opentelemetry-sdk>=1.29.0`, `prometheus-client>=0.20.0`, `pydantic-settings>=2.7`, `pydantic>=2.0`, `pymilvus>=2.4.0`, `python-dotenv>=1.0.0`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py cardiology
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Risk scores are population statistics applied to an individual. They inform, they do not determine.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/cardiology`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
