# Single-Cell Compute Engine

Single-cell transcriptomics compute for the HCLS AI Factory — the compute layer that
feeds the Single-Cell Intelligence Agent. Default port **8573**.

> **Engine vs. Agent — not a duplicate.** "Single-cell" is the one capability that exists as
> both an **engine** and an **agent**, on purpose. This **engine** (`core/engines/single-cell/`,
> :8573, registry `singlecell-compute`) is deterministic scanpy *compute* — it turns a raw
> expression matrix into cell-type clusters. The **agent**
> (`core/agents/single-cell/`, :8130, registry `single-cell-intelligence-agent`) is the RAG
> *reasoning* layer that clinically interprets those clusters (TME profiling, drug response,
> subclonal architecture). **The engine computes; the agent interprets.** The engine is also a
> shared service other components call (e.g. oncology TME, the multi-omics join) — which is
> exactly why it's a horizontal engine, not folded into the agent.

## What it does
Runs the standard [scanpy](https://scanpy.readthedocs.io/) workflow over an `.h5ad`
expression matrix: **QC → normalize → HVG → PCA/neighbors → Leiden clustering →
differential expression → marker-based cell-type annotation**.

- `src/single_cell_compute.py` — the pipeline. The marker-annotation logic
  (`annotate_cluster`, `summarize_clusters`) is pure and unit-tested; the heavy scanpy
  import is lazy so the module loads without the ML stack installed.
- `src/single_cell_service.py` — FastAPI service (`POST /analyze`).
- Ships a canonical PBMC marker panel (Seurat/scanpy tutorial cell types) for annotation.

## Why this engine is small
Unlike the structural-biology engine, single-cell has **no vendored third-party repos** —
it's a thin wrapper over `scanpy` (a pip dependency), so the whole engine is a few hundred
lines of source. Its `data/` directory (example matrices) is gitignored and stays local.

## Run
```bash
pip install scanpy anndata fastapi uvicorn leidenalg
uvicorn single_cell_service:create_app --factory --port 8573   # from src/
# POST /analyze  {"h5ad_path": "...", "resolution": 1.0}  -> clusters + annotations
```

## Roadmap (this engine)
Foundation-model cell embeddings (Geneformer / scGPT) + atlas similarity search are on the
roadmap (D3–D5) — a learned cell encoder on top of the scanpy compute described above.
