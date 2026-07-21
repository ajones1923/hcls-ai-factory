## In plain terms

Every tissue is a mix of thousands of individual cells that are *not* all alike — an immune cell, a
tumor cell, and a healthy cell can sit side by side. The Single-Cell Analysis Engine takes a dataset
that measures gene activity **one cell at a time** (a *single-cell expression matrix*) and works out
**which cell types are present**. It's a real, standard scientific workflow, built on the widely-used
**scanpy** toolkit — the reliable compute step that turns a raw table of numbers into a labeled map of
cells.

## Why it matters

Knowing exactly which cells are in a sample — and, for a tumor, which immune cells surround it — is the
foundation of cancer, immunology, and cell-therapy work. This engine is the honest, reproducible
compute step that every higher-level single-cell interpretation depends on.

*For a patient: knowing exactly which cells make up a tumor — including the immune cells around it —
helps a clinician choose therapies aimed at the right cells.*

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

- **Verified on public data.** Run end-to-end on the public **PBMC 3k** benchmark — about 2,700 human
  immune cells that researchers use as a standard, known-answer test — and it recovered exactly the
  cell types it should. Real-data evidence, not a claim.
- **Deterministic compute, not reasoning.** It computes clusters; interpretation is the agent's job.
- **CPU today, GPU is the upgrade.** It is CPU-served now; a GPU (`rapids-singlecell`) path is the
  planned acceleration.
