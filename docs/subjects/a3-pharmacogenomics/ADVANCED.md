# Pharmacogenomics — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/agents/pharmacogenomics` · 56 Python files · 27,751 LOC · 18 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `pharmacogenomics-intelligence-agent` | agent | **live** | localhost:8507 |

**UI :8507 · API :8508** (platform convention: registry endpoint is the UI, API is UI + 1)

## Principal modules

### `src/knowledge.py`

`get_gene_context`, `get_phenotype_context`, `get_drug_category_context`, `get_hla_context`, `get_inhibitor_context`, `get_alternative_drugs`

- **`get_gene_context`** — Return formatted pharmacogene knowledge for a given gene.
- **`get_phenotype_context`** — Return formatted metabolizer phenotype knowledge.
- **`get_drug_category_context`** — Return formatted drug category knowledge.
- **`get_hla_context`** — Return formatted HLA-drug association knowledge.

### `app/pgx_ui.py`

`init_engine`, `init_agent`, `init_pipeline_components`, `load_knowledge`, `render_evidence_cards`, `run_query`

- **`init_engine`** — Initialize the PGx RAG engine (cached across reruns).
- **`init_agent`** — Initialize the autonomous PGx Intelligence Agent.
- **`init_pipeline_components`** — Initialize PGx pipeline components (cached).
- **`load_knowledge`** — Load PGx knowledge graph data for UI display.

### `src/pgx_pipeline.py`

`AlertSeverity`, `PGxAlert`, `PGxPosition`, `StarAlleleCaller`, `PhenotypeTranslator`, `DrugGeneMatcher`

- **`AlertSeverity`** — CPIC-aligned clinical action levels.
- **`PGxAlert`** — Alert generated when a patient's PGx profile impacts a medication.
- **`PGxPosition`** — A single pharmacogenomic variant position.
- **`StarAlleleCaller`** — Resolves VCF variants to pharmacogene star allele nomenclature.

### `src/dosing.py`

`DosingCalculator`

- **`DosingCalculator`** — Genotype-guided dosing calculators for dose-critical drugs.

### `src/export.py`

`generate_filename`, `export_markdown`, `export_json`, `export_pdf`, `export_fhir_r4`

- **`generate_filename`** — Generate a timestamped filename for export.
- **`export_markdown`** — Export a query result as a Markdown report.
- **`export_json`** — Export a query result as structured JSON.
- **`export_pdf`** — Export a query result as a professionally styled PDF report.


## Dependencies

`anthropic>=0.18.0`, `apscheduler>=3.10.0`, `biopython>=1.83`, `cyvcf2>=0.31.0`, `fastapi>=0.109.0`, `loguru>=0.7.0`, `lxml>=5.0.0`, `numpy>=1.24.0`, `opentelemetry-api>=1.29.0`, `opentelemetry-sdk>=1.29.0`, `prometheus-client>=0.20.0`, `pydantic-settings>=2.7`, `pydantic>=2.0`, `pymilvus>=2.4.0`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py pharmacogenomics
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

CPIC guidance covers specific gene-drug pairs. Substrate relationships outside those pairs are pharmacology, not guideline-backed dosing advice.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/agents/pharmacogenomics`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
