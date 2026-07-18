## In plain terms

The Single-Cell Analysis Engine takes a **single-cell expression matrix** — thousands of individual
cells, each with its own gene-activity profile — and works out **which cell types are present**. It is
a real, standard **scanpy** workflow: the deterministic compute layer that turns a raw matrix into a
labeled map of cells.

## Why it matters

Knowing exactly which cells are in a sample — and, for a tumor, what makes up its microenvironment —
underlies oncology, immunology, and cell therapy. This engine is the honest, reproducible compute step
that every higher-level single-cell reasoning depends on.

## How it works

![Inside the Single-Cell Analysis Engine — QC, normalize, cluster, annotate cell types](../../assets/infographics/pages/singlecell-compute-how.png)
/// caption
A real scanpy workflow: QC → normalize → cluster → annotate. Verified on PBMC 3k. Illustrative.
///

1. **QC & normalize** — filter cells, normalize counts, select highly variable genes.
2. **Reduce & cluster** — **PCA**, then **Leiden** clustering groups cells that behave alike.
3. **Find markers** — differential expression surfaces each cluster's marker genes.
4. **Annotate** — marker genes map clusters to cell types. **Verified on the public PBMC 3k dataset:**
   2,700 cells → 9 clusters → CD4 T, B, NK, CD14+ and FCGR3A+ Monocytes, Dendritic, and Megakaryocytes.

## What goes in, what comes out

- **In:** an **expression matrix** and a clustering **resolution** (`POST /analyze`).
- **Out:** **cell-type clusters** — a labeled cell map.

## Where it fits

![Where the Single-Cell Analysis Engine sits — compute base under the intelligence agent and a shared service](../../assets/infographics/pages/singlecell-compute-fits.png)
/// caption
The engine computes; the Single-Cell Intelligence Agent interprets. A shared service others call. Illustrative.
///

It is the compute base beneath the [Single-Cell Intelligence Agent](../agents/single-cell-intelligence-agent.md)
— the engine computes, the agent interprets — and a **shared service** other capabilities call (e.g.
the oncology tumor-microenvironment view and the multi-omics join).

## Honest limits

- **Verified on public data.** Marked `verified` because it was run end-to-end on the public **PBMC 3k**
  benchmark with the expected, canonical cell types — real-data evidence, not a claim.
- **Deterministic compute, not reasoning.** It computes clusters; interpretation is the agent's job.
- **CPU today, GPU is the upgrade.** It is CPU-served now; a GPU (`rapids-singlecell`) path is the
  planned acceleration.
