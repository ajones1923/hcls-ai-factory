# Therapeutic Discovery — Advanced Learning Guide

For engineers extending or operating this subject.

**Source:** `core/engines/therapeutic-discovery` · 33 Python files · 9,227 LOC · 13 test files

## Registered capabilities

| Capability | Type | Status | Endpoint |
|---|---|---|---|
| `therapeutic-discovery-engine` | engine | **live** | localhost:8505 |
| `molmim-nim` | nim | **live** | localhost:8001 |
| `diffdock-nim` | nim | **live** | localhost:8002 |
| `genmol-nim` | nim | **planned** | — |
| `chemprop-admet` | model | **live** | localhost:8572 |
| `molecule-generator` | model | **live** | localhost:8574 |

**:8505** — a single-service portal, so there is no separate API port

## Principal modules

### `generate_vcp_report_enhanced.py`

`Colors`, `GradientRect`, `VCPReportGeneratorEnhanced`, `main`

- **`Colors`** — Dark theme color palette.
- **`GradientRect`** — Custom flowable for gradient header bars.
- **`VCPReportGeneratorEnhanced`** — Generate stunning VCP Drug Candidate PDF Report.
- **`main`** — Generate enhanced VCP Drug Candidate Report.

### `app/discovery_ui.py`

`load_targets_from_export`, `get_active_target`, `render_header`, `render_sidebar`, `render_target_hypothesis`, `render_structural_evidence`

- **`load_targets_from_export`** — Load targets from Stage 2 RAG/Chat export file.
- **`get_active_target`** — Get the currently active target from session state or load from export.
- **`render_header`** — Render the main header.
- **`render_sidebar`** — Render the sidebar.

### `src/nim_clients.py`

`NIMServiceConfig`, `MolMIMClient`, `DiffDockClient`, `NIMServiceManager`, `CloudMolMIMClient`, `CloudDiffDockClient`

- **`NIMServiceConfig`** — Configuration for NIM service connection.
- **`MolMIMClient`** — Client for MolMIM molecule generation service.
- **`DiffDockClient`** — Client for DiffDock molecular docking service.
- **`NIMServiceManager`** — Manages NIM service connections and provides fallback behavior.

### `generate_vcp_report.py`

`VCPReportGenerator`, `main`

- **`VCPReportGenerator`** — Generate VCP → Drug Candidate PDF Report.
- **`main`** — Generate VCP Drug Candidate Report.

### `src/pipeline.py`

`pediatric_safety_assessment`, `DrugDiscoveryPipeline`, `run_vcp_demo_pipeline`

- **`pediatric_safety_assessment`** — Assess drug candidate safety for pediatric patients.
- **`DrugDiscoveryPipeline`** — Main pipeline orchestrator for drug discovery.
- **`run_vcp_demo_pipeline`** — Run the VCP FTD demo pipeline.


## Dependencies

`loguru==0.7.3`, `numpy==2.4.1`, `opentelemetry-api>=1.29.0`, `opentelemetry-sdk>=1.29.0`, `pillow==12.1.0`, `py3Dmol==2.5.3`, `pydantic==2.12.5`, `pytest-cov==7.0.0`, `pytest==9.0.2`, `rdkit==2025.9.3`, `reportlab==4.4.0`, `requests==2.32.5`, `rich==14.2.0`, `stmol==0.0.9`

## Running the tests

```bash
.venv/bin/python scripts/run_all_tests.py therapeutic-discovery
```

Two traps the shared harness handles, which a hand-rolled `pytest` invocation will not:

1. Several subjects ship `src/collections.py`, which **shadows the Python standard library**. Putting
   their `src/` on `PYTHONPATH` kills the interpreter before collection.
2. `structural-biology/vendor_rfdiffusion/` is vendored third-party code needing gated GPU packages
   and is excluded.

## Operational notes

Molecule generation (MolMIM) and docking (DiffDock) are gated NVIDIA NIMs and are not installed. Candidates shown are pre-computed. This is the flagship demo and the easiest to overclaim.

Before changing a port, read [`../../build/PORT_MAP.md`](../../build/PORT_MAP.md). The convention is
enforced by `scripts/validate_registry.py`, which also cross-checks `health-monitor.sh` — a port
change in one place and not the other fails the build.

## Extending it

1. Add or change code under `core/engines/therapeutic-discovery`.
2. Keep the capability entry in `lib/hcls_common/capabilities.json` truthful — **a `live` status must
   answer a health probe**. Two capabilities were found registered `live` with nothing bound to their
   ports; do not add a third.
3. Run the gate: `ruff`, `pytest lib/hcls_common`, `validate_registry.py`, `run_all_tests.py`.
