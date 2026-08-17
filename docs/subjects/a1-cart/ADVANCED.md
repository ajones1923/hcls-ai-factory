# CAR-T Intelligence — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/cart` · 62 Python files · 22,155 LOC · 11 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `cart-intelligence-agent` | agent | **live** | localhost:8521 |

**UI :8521 · API :8522** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/knowledge.py`

`get_target_context`, `get_toxicity_context`, `get_manufacturing_context`, `get_all_context_for_query`, `get_knowledge_stats`, `get_pediatric_cart_context`

- **`get_target_context`** — Get formatted knowledge context for a target antigen.
- **`get_toxicity_context`** — Get formatted knowledge context for a toxicity profile.
- **`get_manufacturing_context`** — Get formatted knowledge context for a manufacturing process.
- **`get_all_context_for_query`** — Extract all relevant knowledge context from a query string.

### `src/export.py`

`generate_filename`, `export_markdown`, `export_json`, `export_pdf`, `export_fhir_r4`

- **`generate_filename`** — Generate a timestamped filename for export.
- **`export_markdown`** — Export a query result as a Markdown report.
- **`export_json`** — Export a query result as structured JSON.
- **`export_pdf`** — Export a query result as a professionally styled PDF report.

### `src/query_expansion.py`

`expand_query`, `expand_query_by_category`, `get_expansion_stats`

- **`expand_query`** — Extract expansion terms from a user query.
- **`expand_query_by_category`** — Like expand_query but returns terms grouped by expansion category.
- **`get_expansion_stats`** — Return the number of keywords and total terms per expansion map.

### `app/cart_ui.py`

`init_engine`, `init_agent`, `render_evidence_cards`, `build_conversation_context`

- **`init_engine`** — Initialize the CAR-T RAG engine (cached across reruns).
- **`init_agent`** — Initialize the autonomous CAR-T Intelligence Agent.
- **`render_evidence_cards`** — Render evidence cards with relevance indicators.
- **`build_conversation_context`** — Build conversation context from recent exchanges for follow-up queries.

### `src/collections.py`

`CARTCollectionManager`

- **`CARTCollectionManager`** — Manages 11 CAR-T Milvus collections (10 owned + 1 read-only genomic).


## Dependencies

`anthropic>=0.18.0`, `apscheduler>=3.10.0`, `biopython>=1.83`, `fastapi>=0.109.0`, `loguru>=0.7.0`, `lxml>=5.0.0`, `numpy>=1.24.0`, `opentelemetry-api>=1.29.0`, `opentelemetry-sdk>=1.29.0`, `prometheus-client>=0.20.0`, `pydantic-settings>=2.7`, `pydantic>=2.0`, `pymilvus>=2.4.0`, `python-dotenv>=1.0.0`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py cart
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

On-target off-tumour toxicity and cytokine release are real. Any construct suggestion carries its monitoring requirement.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/cart`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
