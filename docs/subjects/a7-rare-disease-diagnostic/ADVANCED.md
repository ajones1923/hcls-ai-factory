# Rare Disease Diagnostic — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/rare-disease-diagnostic` · 48 Python files · 22,454 LOC · 14 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `rare-disease-diagnostic-agent` | agent | **live** | localhost:8544 |

**UI :8544 · API :8545** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/knowledge.py`

_no public symbols_


### `src/clinical_workflows.py`

`BaseRareDiseaseWorkflow`, `PhenotypeDrivenWorkflow`, `WESWGSInterpretationWorkflow`, `MetabolicScreeningWorkflow`, `DysmorphologyWorkflow`, `NeurogeneticWorkflow`

- **`BaseRareDiseaseWorkflow`** — Abstract base for all rare disease diagnostic workflows.
- **`PhenotypeDrivenWorkflow`** — Phenotype-driven differential diagnosis using HPO terms.
- **`WESWGSInterpretationWorkflow`** — WES/WGS variant interpretation workflow.
- **`MetabolicScreeningWorkflow`** — Metabolic screening workflow for inborn errors of metabolism.

### `src/agent.py`

`EvidenceLevel`, `ACMGClassification`, `InheritancePattern`, `SeverityLevel`, `DiagnosticResult`, `SearchPlan`

- **`EvidenceLevel`** — Clinical evidence hierarchy for rare disease diagnostics.
- **`ACMGClassification`** — ACMG/AMP variant classification.
- **`InheritancePattern`** — Mendelian inheritance patterns.
- **`SeverityLevel`** — Finding severity classification.

### `src/rag_engine.py`

`RareDiseaseSearchResult`, `get_all_collection_names`, `RareDiseaseRAGEngine`

- **`RareDiseaseSearchResult`** — A single search result from a Milvus collection.
- **`get_all_collection_names`** — Return all collection names.
- **`RareDiseaseRAGEngine`** — Multi-collection RAG engine for rare disease diagnostics.

### `api/routes/diagnostic_clinical.py`

`integrated_assessment`, `QueryRequest`, `EvidenceItem`, `QueryResponse`, `SearchRequest`, `SearchResult`

- **`integrated_assessment`** — Multi-agent integrated assessment combining insights from across the HCLS AI Factory.
- **`QueryRequest`** — Free-text RAG query with optional workflow and patient context.


## Dependencies

`anthropic==0.25.0`, `apscheduler==3.10.4`, `fastapi==0.111.0`, `httpx==0.27.0`, `loguru==0.7.2`, `lxml==5.2.1`, `numpy==1.26.4`, `prometheus-client==0.20.0`, `pydantic-settings==2.2.1`, `pydantic==2.7.4`, `pymilvus==2.4.1`, `python-docx==1.1.0`, `python-dotenv==1.0.1`, `python-multipart==0.0.9`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py rare-disease-diagnostic
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

A ranked differential is a starting point for investigation, never a diagnosis.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/rare-disease-diagnostic`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
