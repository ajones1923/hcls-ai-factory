"""
Single-cell analysis compute backend (D1).

Turns the single-cell agent from RAG-over-literature into real computation: ingest an
.h5ad expression matrix and run the standard scanpy workflow (QC -> normalize -> HVG ->
PCA -> neighbors -> Leiden clustering -> marker DE), then annotate clusters to cell types
by marker-gene overlap. The clinical RAG agent sits *on top* of these real results.

The marker annotation logic is pure and unit-testable; the heavy scanpy import is lazy.
"""
from __future__ import annotations

from typing import Any

# Canonical PBMC marker panel (Seurat/scanpy tutorial cell types) for annotation.
PBMC_MARKERS: dict[str, list[str]] = {
    "CD4 T cells": ["IL7R", "CD3D", "CD3E", "CCR7"],
    "CD8 T cells": ["CD8A", "CD8B", "CD3D", "GZMK"],
    "B cells": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "NK cells": ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "CD14+ Monocytes": ["CD14", "LYZ", "S100A8", "S100A9"],
    "FCGR3A+ Monocytes": ["FCGR3A", "MS4A7", "LYZ"],
    "Dendritic cells": ["FCER1A", "CST3", "CLEC10A"],
    "Megakaryocytes": ["PPBP", "PF4", "ITGA2B"],
}


# --------------------------------------------------------------------------- #
# Pure, testable logic
# --------------------------------------------------------------------------- #
def annotate_cluster(top_genes: list[str], markers: dict[str, list[str]] = PBMC_MARKERS) -> dict[str, Any]:
    """Annotate a cluster by overlap of its top DE genes with each cell type's markers."""
    top = set(top_genes)
    scores = {ct: len(top & set(genes)) for ct, genes in markers.items()}
    best = max(scores, key=scores.get) if scores else "Unknown"
    n = scores.get(best, 0)
    return {
        "cell_type": best if n > 0 else "Unknown",
        "marker_overlap": n,
        "evidence": sorted(top & set(markers.get(best, []))) if n > 0 else [],
    }


def summarize_clusters(clusters: list[dict[str, Any]],
                       markers: dict[str, list[str]] = PBMC_MARKERS) -> dict[str, Any]:
    """clusters: [{"cluster": id, "n_cells": int, "top_genes": [..]}] -> annotated summary."""
    out = []
    for c in clusters:
        ann = annotate_cluster(c.get("top_genes", []), markers)
        out.append({"cluster": c["cluster"], "n_cells": c.get("n_cells"),
                    "top_genes": c.get("top_genes", [])[:8], **ann})
    types = sorted({c["cell_type"] for c in out if c["cell_type"] != "Unknown"})
    return {"n_clusters": len(out), "cell_types": types, "clusters": out}


# --------------------------------------------------------------------------- #
# Real scanpy pipeline (heavy; lazy import)
# --------------------------------------------------------------------------- #
class SingleCellAnalysis:
    def __init__(self, markers: dict[str, list[str]] = PBMC_MARKERS) -> None:
        self.markers = markers

    def run(self, adata, n_top_genes: int = 2000, resolution: float = 1.0,
            min_genes: int = 200, min_cells: int = 3) -> dict[str, Any]:
        """Run the standard scRNA-seq workflow on an AnnData (or .h5ad path) and annotate."""
        import scanpy as sc
        if isinstance(adata, str):
            adata = sc.read_h5ad(adata)

        # QC + filtering
        sc.pp.filter_cells(adata, min_genes=min_genes)
        sc.pp.filter_genes(adata, min_cells=min_cells)
        # normalize + log
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        # HVG + scale + PCA
        sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars - 1))
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1))
        # neighbors + cluster + DE
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2, directed=False)
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

        # per-cluster top genes + sizes
        names = adata.uns["rank_genes_groups"]["names"]
        sizes = adata.obs["leiden"].value_counts().to_dict()
        clusters = []
        for cl in adata.obs["leiden"].cat.categories:
            top = [names[i][cl] for i in range(min(10, len(names)))]
            clusters.append({"cluster": str(cl), "n_cells": int(sizes.get(cl, 0)), "top_genes": top})

        summary = summarize_clusters(clusters, self.markers)
        summary["n_cells"] = int(adata.n_obs)
        summary["n_genes"] = int(adata.raw.n_vars if adata.raw is not None else adata.n_vars)
        return summary
