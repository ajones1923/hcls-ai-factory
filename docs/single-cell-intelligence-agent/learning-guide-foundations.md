# Single-Cell Intelligence Agent -- Learning Guide: Foundations

**Date:** March 2026
**Author:** Adam Jones

---

## 1. Why Single-Cell Analysis Matters

### 1.1 The Bulk Averaging Problem

Traditional gene expression analysis (microarrays, bulk RNA-seq) measures the average expression of 20,000+ genes across millions of cells in a tissue sample. This averaging obscures critical biological information:

- A tissue containing 50% Type A cells (gene X ON) and 50% Type B cells (gene X OFF) reports gene X at moderate expression -- indistinguishable from a tissue where all cells express gene X at moderate levels
- Rare cell populations (< 5% of total) are invisible in bulk measurements
- Cell state transitions (e.g., T cell activation to exhaustion) appear as continuous gradients rather than discrete states

### 1.2 What Single-Cell Sequencing Reveals

Single-cell RNA sequencing (scRNA-seq) profiles each cell individually, enabling:

- **Cell type identification:** Distinguish 20-50 distinct cell types within a single biopsy
- **State characterization:** Identify activated, resting, exhausted, or cycling states within the same cell type
- **Trajectory inference:** Map developmental and activation pathways
- **Rare cell detection:** Find populations as rare as 0.1% of total cells
- **Heterogeneity quantification:** Measure the diversity within a tumor or tissue

### 1.3 Clinical Impact

| Application | Bulk Resolution | Single-Cell Resolution |
|-------------|----------------|----------------------|
| Immunotherapy selection | PD-L1 TPS (binary: positive/negative) | TME classification (4 phenotypes with treatment recs) |
| Drug resistance | "Tumor resistant" | "5% subclone with KRAS G12C drives resistance" |
| Cell therapy targets | "CD19 expressed on tumor" | "CD19 on 85% of tumor, 1.8 TPM in bone marrow" |
| Biomarker discovery | Tissue-level markers | Cell-type-specific markers with higher specificity |
| Treatment monitoring | "Disease progression" | "Resistant clone expanding from 2% to 15% over 4 weeks" |

---

## 2. Single-Cell Technologies

### 2.1 Droplet-Based (10x Genomics Chromium)

**How it works:** Individual cells are encapsulated in oil droplets with barcoded gel beads. Each cell's mRNA is tagged with a unique cell barcode and unique molecular identifier (UMI), then pooled and sequenced.

**Characteristics:**
- Throughput: 500 - 10,000 cells per run
- Genes detected: 2,000 - 5,000 per cell
- Cost: $1-3 per cell
- Strengths: High throughput, standardized protocols
- Limitations: 3' bias, droplet doublets (2-5%)

**Use cases:** Immune profiling, tumor characterization, atlas building

### 2.2 Plate-Based (Smart-seq2/Smart-seq3)

**How it works:** Individual cells are sorted into 96- or 384-well plates using FACS. Full-length cDNA is generated from each cell and sequenced independently.

**Characteristics:**
- Throughput: 96 - 384 cells per run
- Genes detected: 5,000 - 10,000 per cell
- Cost: $10-50 per cell
- Strengths: Full-length transcripts, high sensitivity, isoform detection
- Limitations: Low throughput, higher cost

**Use cases:** Rare cell characterization, splicing analysis, transcript isoform studies

### 2.3 Spatial Transcriptomics

**How it works:** Gene expression is measured in situ, preserving spatial location within the tissue. Technologies vary in resolution and gene coverage:

| Platform | Company | Resolution | Genes | Method |
|----------|---------|-----------|-------|--------|
| Visium | 10x Genomics | 55 um spots | Whole transcriptome | Spatial barcoding on tissue section |
| MERFISH | Vizgen | Subcellular | 100-500 panel | Multiplexed FISH with error correction |
| Xenium | 10x Genomics | Subcellular | 100-5000 | Padlock probe in situ sequencing |
| CODEX | Akoya | Single cell | 40-60 proteins | Sequential antibody staining |

### 2.4 Multi-Modal Technologies

| Technology | Measurements | Platform |
|-----------|-------------|----------|
| CITE-seq | RNA + 200+ surface proteins | 10x Chromium |
| Multiome | RNA + chromatin accessibility (ATAC) | 10x Chromium |
| scTCR/BCR-seq | RNA + T/B cell receptor sequences | 10x Chromium |
| SHARE-seq | RNA + chromatin accessibility | Custom |

---

## 3. Data Formats

### 3.1 AnnData (.h5ad)

AnnData is the standard data format for single-cell analysis in the Python ecosystem (Scanpy):

```
AnnData object
|-- .X         # Expression matrix (cells x genes), sparse or dense
|-- .obs       # Cell metadata (DataFrame: cell_type, batch, patient, ...)
|-- .var       # Gene metadata (DataFrame: gene_name, highly_variable, ...)
|-- .obsm      # Cell embeddings (PCA, UMAP, tSNE coordinates)
|-- .obsp      # Cell-cell graphs (kNN connectivities, distances)
|-- .uns       # Unstructured data (clustering params, colors, ...)
|-- .layers    # Additional expression matrices (raw counts, normalized)
```

**Example:**
```python
import scanpy as sc

# Load an h5ad file
adata = sc.read_h5ad("sample.h5ad")

print(adata)
# AnnData object with n_obs x n_vars = 10000 x 20000
#     obs: 'cell_type', 'patient', 'tissue'
#     var: 'gene_name', 'highly_variable'
#     obsm: 'X_pca', 'X_umap'

print(adata.obs['cell_type'].value_counts())
# CD8+ T cell      2500
# Macrophage       2000
# B cell           1500
# ...
```

### 3.2 Seurat Object (R)

The R ecosystem uses Seurat objects with similar structure:

```r
# Seurat object slots:
# @assays$RNA       # Expression data (counts, data, scale.data)
# @meta.data        # Cell metadata
# @reductions       # PCA, UMAP, tSNE
# @graphs           # kNN graphs
# @active.ident     # Current cell identity
```

### 3.3 Count Matrix Formats

| Format | Extension | Size | Read Speed |
|--------|-----------|------|-----------|
| Dense matrix | .csv/.tsv | Large | Slow |
| Sparse matrix (MTX) | .mtx + .genes + .barcodes | Medium | Medium |
| HDF5 | .h5/.h5ad | Compact | Fast |
| Loom | .loom | Compact | Fast |
| Zarr | .zarr/ | Compact | Fast (cloud-native) |

---

## 4. The Single-Cell Analysis Pipeline

### 4.1 Overview

```
Raw Data (FASTQ)
      |
      v
Alignment & Quantification (Cell Ranger / STARsolo)
      |
      v
Count Matrix (cells x genes)
      |
      v
+-----+-----+-----+-----+-----+-----+
|     |     |     |     |     |     |
QC   Filter Norm  HVG   PCA  Batch
                              Correction
      |
      v
kNN Graph Construction
      |
      v
+-----+-----+
|           |
UMAP/tSNE   Clustering
             (Leiden/Louvain)
      |
      v
+-----+-----+-----+
|           |     |
Cell Type   DE    Trajectory
Annotation  Genes Inference
```

### 4.2 Quality Control (QC)

QC identifies and removes low-quality cells and artifacts:

| QC Metric | Low-Quality Indicator | Threshold (typical) |
|-----------|----------------------|-------------------|
| nUMI (total counts) | Very low or very high | < 500 or > 50,000 |
| nGenes (detected genes) | Very low | < 200 |
| Mitochondrial % | High (dying/stressed cell) | > 20% |
| Ribosomal % | Very high (lysing cell) | > 50% |
| Doublet score | High (two cells in one droplet) | > 0.25 (Scrublet) |

```python
# Scanpy QC example
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[adata.obs['pct_counts_mt'] < 20]
adata = adata[adata.obs['n_genes_by_counts'] > 200]
adata = adata[adata.obs['n_genes_by_counts'] < 8000]
```

### 4.3 Normalization

Normalization corrects for differences in sequencing depth between cells:

| Method | Library | Description |
|--------|---------|-------------|
| Log-normalize | Scanpy/Seurat | Divide by total counts, multiply by scale factor (10,000), log1p |
| scran | scran (R) | Pool-based size factor estimation, deconvolution |
| SCTransform | Seurat | Regularized negative binomial regression |
| Pearson residuals | Scanpy | Analytic Pearson residuals from expected counts |

```python
# Standard log-normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

### 4.4 Highly Variable Gene (HVG) Selection

Not all 20,000+ genes are informative. HVG selection identifies 2,000-5,000 genes with the highest cell-to-cell variation:

```python
# Select top 2000 highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']]
```

### 4.5 Dimensionality Reduction: PCA

Principal Component Analysis reduces the gene space from 2,000+ dimensions to 30-50 principal components that capture the major axes of variation:

```python
sc.pp.scale(adata)  # Zero-center and unit variance
sc.tl.pca(adata, n_comps=50)
```

### 4.6 Batch Correction

When combining data from multiple patients, technologies, or experiments, batch effects must be corrected:

| Method | Library | Approach |
|--------|---------|----------|
| Harmony | harmonypy | Iterative soft clustering in PCA space |
| scVI | scvi-tools | Variational autoencoder |
| BBKNN | bbknn | Batch-balanced kNN graph |
| Scanorama | scanorama | Mutual nearest neighbor panorama stitching |

### 4.7 kNN Graph and Clustering

Cells are connected in a k-nearest-neighbor graph based on PCA distances, then clustered:

```python
# Build kNN graph (k=30 neighbors)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

# Leiden clustering (resolution controls granularity)
sc.tl.leiden(adata, resolution=0.5)

# UMAP for visualization
sc.tl.umap(adata)
```

**Clustering algorithms:**

| Algorithm | Description | When to Use |
|-----------|-------------|-------------|
| Leiden | Community detection on kNN graph | Default choice, modularity-optimized |
| Louvain | Similar to Leiden, older algorithm | Legacy compatibility |
| k-means | Partition-based, requires k | When k is known a priori |
| HDBSCAN | Density-based, handles noise | Irregular cluster shapes |

### 4.8 Cell Type Annotation

Assigning biological identity to clusters:

| Strategy | Method | Automation |
|----------|--------|-----------|
| Manual annotation | Expert review of marker genes | Low |
| Reference-based | Transfer labels from annotated atlas | High |
| Marker-based | Score known marker gene sets per cluster | Medium |
| Automated | CellTypist, SingleR, scGPT | High |

```python
# Check canonical markers per cluster
marker_genes = {
    'T cells': ['CD3D', 'CD3E'],
    'CD8+ T': ['CD8A', 'CD8B', 'GZMB'],
    'B cells': ['CD19', 'MS4A1', 'CD79A'],
    'Monocytes': ['CD14', 'LYZ', 'S100A9'],
    'NK cells': ['NKG7', 'GNLY', 'NCAM1'],
}
sc.pl.dotplot(adata, marker_genes, groupby='leiden')
```

### 4.9 Differential Expression (DE)

Identify genes that distinguish one cell type/condition from others:

```python
# Wilcoxon rank-sum test (recommended for scRNA-seq)
sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=10)
```

| DE Method | Strengths | Limitations |
|-----------|----------|-------------|
| Wilcoxon rank-sum | Non-parametric, robust to outliers | Does not model count distribution |
| t-test | Fast, parametric | Assumes normality |
| MAST | Models dropout (zero-inflation) | Slower, R-based |
| pseudobulk DESeq2 | Controls for patient effects | Requires replicates |

---

## 5. Key Concepts for Single-Cell Intelligence

### 5.1 Cell Ontology (CL)

The Cell Ontology is a standardized vocabulary for cell types:

| CL ID | Cell Type |
|-------|----------|
| CL:0000084 | T cell |
| CL:0000625 | CD8-positive, alpha-beta T cell |
| CL:0000624 | CD4-positive, alpha-beta T cell |
| CL:0000815 | Regulatory T cell |
| CL:0000236 | B cell |
| CL:0000623 | Natural killer cell |
| CL:0000576 | Monocyte |
| CL:0000235 | Macrophage |

The Single-Cell Intelligence Agent maps all 57 cell types to CL identifiers for standardized communication.

### 5.2 Gene Markers

Canonical marker genes identify cell types:

| Cell Type | Key Markers | Surface Protein |
|----------|-------------|----------------|
| CD8+ T cell | CD8A, CD8B, GZMB, PRF1 | CD8 |
| CD4+ T cell | CD4, IL7R, CCR7 | CD4 |
| Treg | FOXP3, IL2RA, CTLA4 | CD25 |
| B cell | CD19, MS4A1, CD79A | CD19, CD20 |
| NK cell | NKG7, GNLY, NCAM1, KLRD1 | CD56 |
| Monocyte | CD14, LYZ, S100A9 | CD14 |
| Macrophage | CD68, CD163, MARCO | CD68 |
| Dendritic cell | CLEC9A, CD1C, FCER1A | CD11c |
| Fibroblast | COL1A1, DCN, PDGFRA | PDGFRA |
| Endothelial | PECAM1, VWF, CDH5 | CD31 |
| Epithelial | EPCAM, KRT18, CDH1 | EpCAM |

### 5.3 Tumor Microenvironment (TME)

The TME describes the cellular ecosystem surrounding and within a tumor:

| TME Class | Characteristics | Immunotherapy Response |
|-----------|----------------|----------------------|
| Hot-inflamed | High CD8+ T cell infiltration, active immune response | Checkpoint inhibitor responsive |
| Cold-desert | Minimal immune infiltrate, low neoantigen load | Requires immune priming |
| Excluded | Immune cells at tumor margin, blocked from infiltrating | Target stromal barriers |
| Immunosuppressive | Immune cells present but suppressed (Tregs, M2, MDSCs) | Dual checkpoint or depletion |

### 5.4 Embeddings and Vector Search

The Single-Cell Intelligence Agent stores knowledge as 384-dimensional vectors:

1. Text (cell type description, drug mechanism, clinical evidence) is encoded by BGE-small-en-v1.5
2. Vectors are stored in Milvus with metadata fields
3. User queries are embedded with the same model
4. Cosine similarity identifies the most relevant records
5. Top-K results are synthesized by Claude into a coherent response

---

## 6. Tools and Libraries

### 6.1 Python Ecosystem

| Tool | Purpose | Link |
|------|---------|------|
| Scanpy | Core single-cell analysis | scanpy.readthedocs.io |
| AnnData | Data structure | anndata.readthedocs.io |
| scvi-tools | Probabilistic models | scvi-tools.org |
| CellTypist | Automated cell typing | celltypist.org |
| Squidpy | Spatial analysis | squidpy.readthedocs.io |
| scVelo | RNA velocity | scvelo.org |
| CellPhoneDB | Cell-cell interaction | cellphonedb.org |

### 6.2 R Ecosystem

| Tool | Purpose |
|------|---------|
| Seurat | Core single-cell analysis |
| SingleR | Reference-based annotation |
| Monocle3 | Trajectory inference |
| CellChat | Cell-cell communication |
| Harmony | Batch correction |

### 6.3 Reference Databases

| Database | Content | URL |
|----------|---------|-----|
| Human Cell Atlas | Cell type references | humancellatlas.org |
| CellxGene | Public datasets | cellxgene.cziscience.com |
| CellMarker 2.0 | Marker-cell type associations | bio-bigdata.hrbmu.edu.cn/CellMarker |
| PanglaoDB | Marker gene database | panglaodb.se |
| Cell Ontology | Standardized cell vocabulary | obofoundry.org/ontology/cl |

---

## 7. Recommended Learning Path

### 7.1 Beginner (Week 1-2)

1. Read the Scanpy tutorial: "Preprocessing and clustering 3k PBMCs"
2. Load a public dataset from CellxGene and explore cell type annotations
3. Run the Single-Cell Intelligence Agent Demo 1 (Cell Type Annotation)
4. Study the 57 cell types in the knowledge base (`src/knowledge.py`)

### 7.2 Intermediate (Week 3-4)

1. Complete the Scanpy PBMC tutorial end-to-end (QC through DE)
2. Run Demo 2 (TME Profiling) and understand classification logic
3. Study the TMEClassifier decision tree in `src/decision_support.py`
4. Explore the 12 Milvus collection schemas in `src/collections.py`

### 7.3 Advanced (Week 5-8)

1. Work through the Advanced Learning Guide
2. Run all 5 demos with custom parameters
3. Study the 10 clinical workflows in `src/clinical_workflows.py`
4. Explore spatial transcriptomics analysis with Squidpy

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Learning Guide: Foundations*

---

!!! warning "Clinical Decision Support Disclaimer"
    The Single-Cell Intelligence Agent is a clinical decision support research tool for single-cell genomic analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
