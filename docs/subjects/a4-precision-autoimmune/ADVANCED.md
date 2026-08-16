# Precision Autoimmune — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/precision-autoimmune` · 45 Python files · 20,427 LOC · 9 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `precision-autoimmune-agent` | agent | **live** | localhost:8531 |

**UI :8531 · API :8532** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `app/autoimmune_ui.py`

`init_services`, `render_sidebar`

- **`init_services`** — Initialize all backend services (cached).

### `scripts/patient_maya.py`

`generate`

- **`generate`** — Generate all clinical documents for Maya Rodriguez.

### `scripts/patient_sarah.py`

`generate`

- **`generate`** — Generate all clinical documents for Sarah Mitchell.

### `src/knowledge.py`

_no public symbols_


### `scripts/patients_dismissals.py`

`generate_maya_dismissals`, `generate_sarah_dismissals`, `generate_linda_dismissals`, `generate_david_dismissals`

- **`generate_maya_dismissals`** — Generate dismissal documents for Maya Rodriguez.
- **`generate_sarah_dismissals`** — Generate dismissal documents for Sarah Mitchell.
- **`generate_linda_dismissals`** — Generate dismissal documents for Linda Chen.
- **`generate_david_dismissals`** — Generate dismissal documents for David Park.


## Dependencies

`PyPDF2>=3.0.0`, `anthropic>=0.18.0`, `fastapi>=0.109.0`, `httpx>=0.25.0`, `loguru>=0.7.0`, `numpy>=1.24.0`, `pandas>=2.0.0`, `plotly>=5.18.0`, `prometheus-client>=0.20.0`, `pydantic-settings>=2.7`, `pydantic>=2.0`, `pymilvus>=2.4.0`, `python-dotenv>=1.0.0`, `python-multipart>=0.0.6`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py precision-autoimmune
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Decision support only; autoimmune management is highly individual and specialist-led.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/precision-autoimmune`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
