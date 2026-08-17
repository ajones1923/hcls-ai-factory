# Neurology Intelligence — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/neurology` · 48 Python files · 22,440 LOC · 14 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `neurology-intelligence-agent` | agent | **live** | localhost:8535 |

**UI :8535 · API :8536** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/clinical_workflows.py`

`BaseNeuroWorkflow`, `AcuteStrokeTriageWorkflow`, `DementiaEvaluationWorkflow`, `EpilepsyFocusWorkflow`, `BrainTumorGradingWorkflow`, `MSMonitoringWorkflow`

- **`BaseNeuroWorkflow`** — Abstract base for all neurology clinical workflows.
- **`AcuteStrokeTriageWorkflow`** — Acute stroke triage: NIHSS, ASPECTS, tPA eligibility, thrombectomy
- **`DementiaEvaluationWorkflow`** — Dementia evaluation: atrophy pattern -> differential (AD vs FTD vs DLB
- **`EpilepsyFocusWorkflow`** — Epilepsy focus localization: ILAE seizure classification, syndrome

### `src/agent.py`

`NeuroWorkflowType`, `EvidenceLevel`, `SeverityLevel`, `ImagingModality`, `ElectrophysiologyType`, `NeuroResponse`

- **`NeuroWorkflowType`** — Types of neurology query workflows.
- **`EvidenceLevel`** — Clinical evidence hierarchy for neurological recommendations.
- **`SeverityLevel`** — Finding severity classification.
- **`ImagingModality`** — Neuroimaging modalities.

### `api/routes/neuro_clinical.py`

`integrated_assessment`, `QueryRequest`, `QueryResponse`, `SearchRequest`, `SearchResult`, `SearchResponse`

- **`integrated_assessment`** — Multi-agent integrated assessment combining insights from across the HCLS AI Factory.
- **`QueryRequest`** — Free-text RAG query with optional domain and patient context.
- **`SearchRequest`** — Multi-collection semantic search.

### `src/knowledge.py`

`get_domain_count`, `get_drug_count`, `get_gene_count`, `get_scale_count`, `get_drugs_by_domain`, `get_genes_by_domain`

- **`get_domain_count`** — Return the number of neurological disease domains.
- **`get_drug_count`** — Return the number of curated neurology drugs.
- **`get_gene_count`** — Return the number of curated neurology genes.
- **`get_scale_count`** — Return the number of clinical scales.

### `src/rag_engine.py`

`NeuroSearchResult`, `get_all_collection_names`, `NeuroRAGEngine`

- **`NeuroSearchResult`** — A single search result from a Milvus collection.
- **`get_all_collection_names`** — Return all collection names.
- **`NeuroRAGEngine`** — Multi-collection RAG engine for neurology intelligence.


## Dependencies

`anthropic==0.25.0`, `apscheduler==3.10.4`, `fastapi==0.111.0`, `httpx==0.27.0`, `loguru==0.7.2`, `lxml==5.2.1`, `numpy==1.26.4`, `prometheus-client==0.20.0`, `pydantic-settings==2.2.1`, `pydantic==2.7.4`, `pymilvus==2.4.1`, `python-docx==1.1.0`, `python-dotenv==1.0.1`, `python-multipart==0.0.9`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py neurology
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Scales such as NIHSS are clinician-administered instruments. The agent supports their interpretation, it does not perform the examination.

Before changing a port, read [`../../build/PORT_MAP.md`](https://github.com/ajones1923/hcls-ai-factory/blob/main/docs/build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/neurology`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
