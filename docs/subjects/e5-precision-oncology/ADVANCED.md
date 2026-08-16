# Precision Oncology — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/precision-oncology/agent` · 67 Python files · 21,723 LOC · 13 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `precision-oncology-agent` | engine | **live** | localhost:8526 |

**UI :8526 · API :8527** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/knowledge.py`

`get_pediatric_dosing_context`, `classify_variant_actionability`, `get_target_context`, `get_therapy_context`, `get_resistance_context`, `get_pathway_context`

- **`get_pediatric_dosing_context`** — Return formatted pediatric dosing context for a given drug or principle.
- **`classify_variant_actionability`** — Classify a variant's actionability against known ACTIONABLE_TARGETS.
- **`get_target_context`** — Return formatted knowledge context for an actionable target gene.
- **`get_therapy_context`** — Return formatted knowledge context for a therapy.

### `src/export.py`

`export_markdown`, `export_json`, `export_pdf`, `export_fhir_r4`, `case_to_markdown`, `markdown_to_pdf`

- **`export_markdown`** — Generate a formatted Markdown report suitable for Molecular Tumor Board
- **`export_json`** — Export MTB packet or agent response as structured JSON.
- **`export_pdf`** — Generate an NVIDIA-themed PDF report via ReportLab.
- **`export_fhir_r4`** — Export an MTB packet as a FHIR R4 Bundle.

### `src/rag_engine.py`

`OncoRAGEngine`

- **`OncoRAGEngine`** — Multi-collection RAG engine for precision oncology queries.

### `src/therapy_ranker.py`

`TherapyRanker`

- **`TherapyRanker`** — Evidence-based therapy ranking engine for precision oncology.

### `src/query_expansion.py`

`expand_query`

- **`expand_query`** — Return up to 10 expansion terms for *query* by checking all maps.


## Dependencies

`anthropic>=0.18.0`, `apscheduler>=3.10.0`, `biopython>=1.83`, `cyvcf2>=0.30.0`, `fastapi>=0.109.0`, `fhir.resources>=7.0.0`, `loguru>=0.7.0`, `lxml>=5.0.0`, `numpy>=1.24.0`, `opentelemetry-api>=1.29.0`, `opentelemetry-sdk>=1.29.0`, `prometheus-client>=0.20.0`, `pydantic-settings>=2.7`, `pydantic>=2.0`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py precision-oncology
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Trial matching is informational. Therapy selection remains the treating oncologist's decision.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/precision-oncology/agent`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
