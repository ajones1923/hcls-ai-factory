# Clinical Trial Intelligence — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/clinical-trial` · 47 Python files · 23,396 LOC · 14 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `clinical-trial-intelligence-agent` | agent | **live** | localhost:8538 |

**UI :8538 · API :8539** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/clinical_workflows.py`

`BaseTrialWorkflow`, `ProtocolDesignWorkflow`, `PatientMatchingWorkflow`, `SiteSelectionWorkflow`, `EligibilityOptimizationWorkflow`, `AdaptiveDesignWorkflow`

- **`BaseTrialWorkflow`** — Abstract base for all clinical trial workflows.
- **`ProtocolDesignWorkflow`** — Protocol design workflow that generates evidence-based protocol
- **`PatientMatchingWorkflow`** — Patient-trial matching workflow that evaluates a patient profile
- **`SiteSelectionWorkflow`** — Site selection and ranking workflow based on weighted scoring

### `src/knowledge.py`

_no public symbols_


### `src/agent.py`

`EvidenceLevel`, `TrialPhase`, `RegulatoryBody`, `SeverityLevel`, `TrialResponse`, `SearchPlan`

- **`EvidenceLevel`** — Clinical evidence hierarchy for trial recommendations.
- **`TrialPhase`** — Clinical trial phases.
- **`RegulatoryBody`** — Regulatory bodies referenced in trial intelligence.
- **`SeverityLevel`** — Finding severity classification.

### `src/rag_engine.py`

`TrialSearchResult`, `get_all_collection_names`, `TrialRAGEngine`

- **`TrialSearchResult`** — A single search result from a Milvus collection.
- **`get_all_collection_names`** — Return all collection names.
- **`TrialRAGEngine`** — Multi-collection RAG engine for clinical trial intelligence.

### `src/query_expansion.py`

`QueryExpander`

- **`QueryExpander`** — Expand clinical trial queries with related domain terms.


## Dependencies

`anthropic==0.25.0`, `apscheduler==3.10.4`, `fastapi==0.111.0`, `httpx==0.27.0`, `loguru==0.7.2`, `lxml==5.2.1`, `numpy>=1.24.0`, `prometheus-client==0.20.0`, `pydantic-settings==2.2.1`, `pydantic==2.7.4`, `pymilvus==2.4.1`, `python-docx==1.1.0`, `python-dotenv>=1.0.0`, `python-multipart==0.0.9`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py clinical-trial
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Matching is informational. Eligibility is determined by the trial investigators, not by this agent.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/clinical-trial`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
