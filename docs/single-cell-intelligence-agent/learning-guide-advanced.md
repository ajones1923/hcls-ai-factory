# Single-Cell Intelligence Agent -- Learning Guide: Advanced Topics

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## 1. TME Classification: Deep Dive

### 1.1 The Four Immunophenotypes

Tumor microenvironment classification is central to immunotherapy patient selection. The field has converged on four canonical phenotypes, each with distinct cellular composition, spatial organization, and therapeutic implications.

#### Hot-Inflamed TME

**Cellular hallmarks:**
- CD8+ T cells > 15% of total cells
- Active cytotoxic gene program (GZMB+, PRF1+, IFNG+)
- PD-L1 expression on tumor and immune cells
- Tertiary lymphoid structures may be present

**Molecular signatures:**
- Interferon-gamma signaling: STAT1, IRF1, CXCL9, CXCL10, CXCL11
- Cytotoxic effector program: GZMA, GZMB, PRF1, GNLY, NKG7
- Antigen presentation: HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2

**Clinical implication:** Strong candidate for checkpoint inhibitor monotherapy. Response rates: 30-50% for anti-PD-1/PD-L1.

#### Cold-Desert TME

**Cellular hallmarks:**
- Total immune infiltrate < 10% of cells
- Minimal T cell presence (< 2% CD8+)
- Low neoantigen load (low TMB)
- No tertiary lymphoid structures

**Molecular signatures:**
- Absent IFN-gamma signaling
- Low MHC class I expression
- Active Wnt/beta-catenin signaling (immune exclusion)
- PTEN loss (PI3K pathway activation)

**Clinical implication:** Checkpoint inhibitors alone are ineffective. Requires immune priming:
- Oncolytic virus (T-VEC) to induce immunogenic cell death
- STING agonist to activate innate immunity
- Radiation therapy for abscopal effect
- Bispecific T-cell engagers (BiTEs) to bypass recruitment defect

#### Excluded TME

**Cellular hallmarks:**
- Immune cells present at tumor margin but excluded from core
- Dense stromal barrier (fibroblasts, CAFs, myofibroblasts)
- Immune cells "stuck" at the invasive front
- Angiogenic vasculature without immune migration signals

**Molecular signatures:**
- TGF-beta signaling: TGFB1, TGFB2, SMAD2, SMAD3
- Stromal activation: COL1A1, COL1A2, FN1, POSTN
- CXCL12/CXCR4 axis (immune cell trapping)
- VEGF-driven angiogenesis without CXCL9/10 chemokines

**Clinical implication:** Target the stromal barrier:
- Anti-TGF-beta (bintrafusp alfa) to reduce fibrosis
- Anti-VEGF + anti-PD-L1 combination (atezolizumab + bevacizumab)
- FAK inhibitors to disrupt stromal architecture
- Anti-CXCR4 to release trapped immune cells

#### Immunosuppressive TME

**Cellular hallmarks:**
- Immune cells present and infiltrating but functionally suppressed
- High regulatory T cell (Treg) fraction (> 10%)
- M2-polarized macrophages dominant
- Myeloid-derived suppressor cells (MDSCs) present

**Molecular signatures:**
- Immunosuppressive cytokines: IL-10, TGF-beta, IL-35
- Metabolic suppression: IDO1, ARG1, NOS2 (tryptophan/arginine depletion)
- Exhaustion markers on T cells: LAG3, TIM3/HAVCR2, TIGIT, TOX
- Checkpoint overexpression: CTLA-4, PD-1 on T cells; PD-L1, PD-L2 on myeloid

**Clinical implication:** Multi-pronged approach:
- Dual checkpoint (anti-PD-1 + anti-CTLA-4)
- Treg depletion (anti-CCR8, low-dose cyclophosphamide)
- Macrophage reprogramming (CSF1R inhibitor)
- MDSC differentiation (ATRA, HDAC inhibitor)

### 1.2 Classification Algorithm

The Single-Cell Intelligence Agent's TMEClassifier implements a hierarchical decision tree:

```
Step 1: Spatial override
  - "absent" + immune < 0.05  -->  COLD_DESERT
  - "margin" + immune > 0.05  -->  EXCLUDED

Step 2: Hot-inflamed check
  - CD8 >= 15% AND immune >= 25%
    - Suppressive > 0.4  -->  IMMUNOSUPPRESSIVE
    - Otherwise           -->  HOT_INFLAMED

Step 3: Excluded check
  - Immune >= 10% AND stromal > 20%  -->  EXCLUDED

Step 4: Immunosuppressive check
  - Suppressive > 0.3 AND immune >= 10%  -->  IMMUNOSUPPRESSIVE

Step 5: Cold check
  - Immune < 10%  -->  COLD_DESERT

Step 6: PD-L1 rescue
  - PD-L1 high AND CD8 >= 5%  -->  HOT_INFLAMED

Default: COLD_DESERT
```

The suppressive score is a weighted combination:
- 50%: suppressive cell fraction (Treg + MDSC + M2 macrophage) / 0.2
- 50%: suppressive gene score (IDO1, TGFB1, IL10, VEGFA, ARG1, NOS2)

### 1.3 Evidence Levels for TME Classification

| Evidence Available | Level | Confidence |
|-------------------|-------|-----------|
| Spatial context + PD-L1 TPS + scRNA-seq | STRONG | High |
| PD-L1 TPS + scRNA-seq (no spatial) | MODERATE | Medium |
| scRNA-seq only (no PD-L1, no spatial) | LIMITED | Low |

---

## 2. Subclonal Architecture and Clonal Dynamics

### 2.1 Why Subclones Matter

Cancer is not a monolithic disease. A single tumor contains multiple subclonal populations, each with distinct:
- Somatic mutation profiles (driver and passenger)
- Copy number aberrations (gains, losses, LOH)
- Transcriptomic programs (proliferation, invasion, immune evasion)
- Drug sensitivity profiles

Under therapeutic selective pressure, resistant subclones expand:

```
Before Treatment:
  Clone A (80%): Drug-sensitive, antigen+
  Clone B (15%): Moderate sensitivity, antigen+
  Clone C (5%):  Resistant, antigen-negative

After 8 Weeks of CAR-T:
  Clone A (5%):  Depleted by CAR-T
  Clone B (20%): Partially depleted
  Clone C (75%): Expanded (antigen escape)
```

### 2.2 Single-Cell Subclonal Detection

Methods for inferring subclonal architecture from scRNA-seq:

| Method | Input | Output | Mechanism |
|--------|-------|--------|-----------|
| inferCNV | scRNA-seq expression | Clone-specific CNV profiles | Expression deviation from normal reference |
| CopyKAT | scRNA-seq expression | Aneuploid/diploid classification | Bayesian segmentation |
| Numbat | scRNA-seq + genotype | Haplotype-aware CNV + clone tree | Allele-specific expression |
| clonealign | scRNA-seq + scDNA-seq | Clone-to-transcriptome mapping | Statistical alignment |

### 2.3 Escape Risk Scoring

The SubclonalRiskScorer evaluates four risk factors per clone:

| Factor | Weight | Threshold |
|--------|--------|-----------|
| Antigen-negative (expression < 0.1) | +0.4 | Binary flag |
| Clone expanding | +0.2 | Boolean (serial samples) |
| High proliferation index | up to +0.2 | Proportional to MKI67/TOP2A |
| Resistance genes present | +0.05/gene (max +0.2) | Count of resistance-associated genes |

**Overall risk classification:**
- HIGH: antigen-negative fraction > 10%
- MEDIUM: antigen-negative > 3% or any individual clone at HIGH risk
- LOW: all clones below thresholds

**Timeline estimation:**
Using exponential growth: `t = T_doubling * log2(0.5 / current_fraction)`

Example: If antigen-negative fraction is 5% and tumor doubling time is 14 days:
- `t = 14 * log2(0.5 / 0.05) = 14 * 3.32 = 46.5 days` to reach 50% dominance

---

## 3. Spatial Transcriptomics

### 3.1 Technology Landscape

Spatial transcriptomics preserves the physical location of gene expression measurements within tissue:

#### Visium (10x Genomics)

- **Resolution:** 55-micron spots (5-10 cells per spot)
- **Coverage:** Whole transcriptome (~20,000 genes)
- **Tissue:** Fresh-frozen or FFPE
- **Workflow:** Tissue on barcoded slide -> permeabilization -> mRNA capture -> sequencing
- **Analysis:** Requires computational deconvolution (cell2location, RCTD) to resolve cell types within spots

#### MERFISH (Vizgen)

- **Resolution:** Subcellular (individual transcripts)
- **Coverage:** 100-500 gene panel (custom design)
- **Tissue:** Fresh-frozen
- **Workflow:** Tissue on slide -> sequential rounds of hybridization + imaging
- **Analysis:** Direct cell segmentation and gene assignment

#### Xenium (10x Genomics)

- **Resolution:** Subcellular
- **Coverage:** 100-5,000 gene panel (expanding)
- **Tissue:** Fresh-frozen or FFPE
- **Workflow:** In situ padlock probe hybridization + rolling circle amplification
- **Analysis:** Cell segmentation -> direct transcript counting per cell

#### CODEX (Akoya)

- **Resolution:** Single cell
- **Coverage:** 40-60 proteins (antibody panel)
- **Tissue:** FFPE or fresh-frozen
- **Workflow:** Sequential antibody staining + fluorescence imaging
- **Analysis:** Protein co-expression -> cell typing

### 3.2 Spatial Analysis Methods

| Analysis | Method | What It Reveals |
|----------|--------|----------------|
| Spatial autocorrelation | Moran's I | Genes with spatially structured expression |
| Niche identification | Cell neighborhood analysis | Co-occurring cell type combinations |
| Cell-cell proximity | Pairwise distance analysis | Which cell types are physically adjacent |
| Spatial deconvolution | cell2location, RCTD | Cell type composition of Visium spots |
| Tissue segmentation | Histological features + expression | Tumor vs. stroma vs. necrosis regions |
| Spatial communication | MISTy, SpaTalk | Location-aware ligand-receptor analysis |

### 3.3 Spatial Niches in Oncology

Clinically relevant spatial patterns:

| Spatial Niche | Cell Types | Clinical Significance |
|--------------|-----------|---------------------|
| Tumor-immune interface | CD8+ T, tumor, DC | Active immune surveillance, checkpoint response |
| Tertiary lymphoid structure | B cell, T cell, FDC | Positive prognosis, improved immunotherapy response |
| Fibrotic barrier | CAF, myofibroblast | Immune exclusion, anti-TGFb target |
| Hypoxic core | Tumor, few immune | Radioresistance, angiogenesis driver |
| Perivascular niche | Endothelial, pericyte, tumor | Metastatic dissemination route |
| Necrotic zone | Dead/dying cells | Antigen release, DAMP signaling |

---

## 4. Trajectory Inference

### 4.1 What Are Cellular Trajectories?

Single-cell snapshots capture cells at different stages of continuous processes (differentiation, activation, exhaustion). Trajectory inference algorithms order cells along these continuous paths in "pseudotime."

### 4.2 Trajectory Types

| Type | Start State | End State | Clinical Relevance |
|------|-----------|----------|-------------------|
| Differentiation | Progenitor/stem | Mature cell | HSC transplant engraftment |
| Activation | Naive T cell | Effector T cell | Immune response quality |
| Exhaustion | Effector T cell | Exhausted T cell (TOX+) | Checkpoint inhibitor response |
| EMT | Epithelial | Mesenchymal | Metastatic potential |
| Stemness | Differentiated tumor | Cancer stem cell | Treatment resistance |
| Cell cycle | G1 | G2/M | Proliferation rate, chemo sensitivity |

### 4.3 Trajectory Inference Methods

| Method | Approach | Strengths | Key Paper |
|--------|---------|----------|-----------|
| Monocle3 | Principal graph | Handles branching, scalable | Cao et al., Nature 2019 |
| PAGA | Partition-based | Robust, preserves topology | Wolf et al., Genome Biology 2019 |
| RNA velocity (scVelo) | Spliced/unspliced ratios | Directionality without time series | Bergen et al., Nature Biotech 2020 |
| Palantir | Diffusion maps | Probabilistic fate assignment | Setty et al., Nature Biotech 2019 |
| CytoTRACE | Gene counts as proxy | Simple, no assumptions | Gulati et al., Science 2020 |

### 4.4 RNA Velocity

RNA velocity infers the direction and speed of gene expression change by comparing unspliced (nascent) and spliced (mature) mRNA:

- **Positive velocity (unspliced > expected):** Gene is being upregulated
- **Negative velocity (unspliced < expected):** Gene is being downregulated
- **Zero velocity (equilibrium):** Gene is at steady state

```python
import scvelo as scv

# Load data with spliced/unspliced counts
adata = scv.read("sample.h5ad")

# Compute velocity
scv.pp.moments(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

# Visualize on UMAP
scv.pl.velocity_embedding_stream(adata)
```

---

## 5. Foundation Models for Single-Cell Biology

### 5.1 scGPT

**Architecture:** Transformer-based generative pre-trained model for single-cell data.

**Pre-training:** 33 million cells from CellxGene, trained on gene expression prediction using masked token modeling.

**Capabilities:**
- Zero-shot cell type annotation
- Gene expression imputation
- Perturbation response prediction
- Multi-batch integration
- Gene regulatory network inference

**Performance benchmarks (from Cui et al., Nature Methods 2024):**
- Cell type annotation: 93.5% accuracy (zero-shot on held-out datasets)
- Batch integration: superior to scVI on 6/8 benchmarks
- Perturbation prediction: R=0.85 correlation with observed perturbation effects

### 5.2 Geneformer

**Architecture:** BERT-style transformer trained on gene expression rank order.

**Pre-training:** 30 million cells from public data, using attention-based gene embeddings.

**Key innovation:** Represents cells as ordered sequences of genes (ranked by expression), enabling transfer learning across tissues and species.

**Capabilities:**
- Context-aware gene function prediction
- Disease state classification
- Therapeutic target nomination
- Dosage sensitivity prediction

**Performance (from Theodoris et al., Nature 2023):**
- Transfer learning accuracy: 85-95% across tissue types
- Network biology prediction: improved over expression-based methods
- Chromatin dynamics prediction: validated experimentally

### 5.3 scFoundation

**Architecture:** Large-scale pre-trained model (100M+ parameters) for cell representation learning.

**Pre-training:** 50 million+ cells from diverse tissues and species.

**Capabilities:**
- Universal cell embeddings for cross-dataset integration
- Drug response prediction
- Cell fate prediction

### 5.4 Integration with the Single-Cell Intelligence Agent

Foundation models can serve as:
1. **Embedding backbone:** Replace BGE-small-en-v1.5 with scGPT cell embeddings for cell-level vector search
2. **Annotation engine:** Zero-shot cell type prediction via scGPT
3. **Perturbation simulator:** Predict drug response at single-cell resolution
4. **Integration layer:** Cross-dataset harmonization via Geneformer embeddings

The agent's knowledge base documents these models and their capabilities. NIM endpoint integration is planned for v2.0.

---

## 6. GPU Benchmarks for Single-Cell Analysis

### 6.1 RAPIDS vs. CPU Benchmarks

| Operation | Dataset Size | CPU (seconds) | GPU (seconds) | Speedup |
|-----------|-------------|--------------|--------------|---------|
| PCA (50 comps) | 50K cells | 45 | 0.9 | 50x |
| PCA (50 comps) | 500K cells | 480 | 4.2 | 114x |
| UMAP | 50K cells | 120 | 2.4 | 50x |
| UMAP | 500K cells | 1,800 | 12 | 150x |
| kNN (k=30) | 50K cells | 90 | 0.8 | 112x |
| kNN (k=30) | 500K cells | 960 | 3.5 | 274x |
| Leiden (res=0.5) | 50K cells | 30 | 1.0 | 30x |
| Leiden (res=0.5) | 500K cells | 350 | 5.0 | 70x |
| Full pipeline | 50K cells | 345 | 7.5 | 46x |
| Full pipeline | 500K cells | 3,590 | 24.7 | 145x |

*Benchmarks on NVIDIA A100 80GB. CPU benchmarks on AMD EPYC 7742 64-core.*

### 6.2 Memory Requirements

| Dataset Size | CPU RAM | GPU VRAM |
|-------------|---------|---------|
| 10K cells | 2 GB | 1 GB |
| 50K cells | 8 GB | 4 GB |
| 100K cells | 16 GB | 8 GB |
| 500K cells | 64 GB | 32 GB |
| 1M cells | 128 GB | 64 GB |

### 6.3 rapids-singlecell

The `rapids-singlecell` package provides GPU-accelerated Scanpy-compatible functions:

```python
import rapids_singlecell as rsc

# GPU-accelerated preprocessing
rsc.pp.normalize_total(adata)
rsc.pp.log1p(adata)
rsc.pp.highly_variable_genes(adata)
rsc.pp.pca(adata)

# GPU-accelerated analysis
rsc.pp.neighbors(adata)
rsc.tl.leiden(adata)
rsc.tl.umap(adata)

# Results are identical to Scanpy, 50-150x faster
```

---

## 7. Cell-Cell Communication Analysis

### 7.1 Ligand-Receptor Databases

| Database | Interactions | Source |
|----------|------------|--------|
| CellPhoneDB | 2,500+ | Curated from literature |
| CellTalkDB | 3,000+ | Curated + predicted |
| NicheNet | 6,000+ | Ligand-target predicted |
| CellChatDB | 2,000+ | Curated with pathway context |
| LIANA | Meta-database | Consensus of multiple databases |

### 7.2 Analysis Methods

| Method | Approach | Output |
|--------|---------|--------|
| CellPhoneDB | Permutation test on L-R co-expression | P-values per L-R pair per cell type pair |
| CellChat | Quantitative mass-action model | Interaction strength, pathway activity |
| NicheNet | Ligand activity prediction from target genes | Ligand prioritization by downstream effect |
| LIANA | Consensus of multiple methods | Aggregated interaction scores |

### 7.3 The Single-Cell Intelligence Agent's L-R Knowledge

The agent curates 25 ligand-receptor pairs across clinically actionable pathways:

| Pathway | Ligand | Receptor | Clinical Relevance |
|---------|--------|----------|-------------------|
| Checkpoint | CD274 (PD-L1) | PDCD1 (PD-1) | Checkpoint inhibitor target |
| Checkpoint | CD80 | CTLA4 | Ipilimumab target |
| Chemokine | CXCL12 | CXCR4 | Immune cell migration/trapping |
| Chemokine | CCL2 | CCR2 | Monocyte/macrophage recruitment |
| Growth factor | EGF | EGFR | TKI target (erlotinib, osimertinib) |
| Growth factor | HGF | MET | MET inhibitor target |
| Notch | DLL1 | NOTCH1 | Cancer stem cell maintenance |
| Wnt | WNT5A | FZD5 | Immune exclusion, beta-catenin |
| Hedgehog | SHH | PTCH1 | Stromal activation |
| Angiogenesis | VEGFA | KDR (VEGFR2) | Anti-VEGF target (bevacizumab) |

---

## 8. Biomarker Discovery at Single-Cell Resolution

### 8.1 Advantages Over Bulk Discovery

| Feature | Bulk Discovery | Single-Cell Discovery |
|---------|---------------|---------------------|
| Specificity | Tissue-level | Cell-type-specific (AUROC > 0.9) |
| Confounders | Cell composition changes confound DE | Direct cell-type DE |
| Sensitivity | Rare cell markers diluted | Detectable at 0.1% frequency |
| Actionability | Unknown cellular source | Known cell type enables targeted therapy |

### 8.2 Discovery Workflow

```
scRNA-seq data (disease vs. control)
         |
         v
Cell type annotation (57 cell types)
         |
         v
Per-cell-type differential expression
         |
    +----+----+----+
    |    |    |    |
    v    v    v    v
  CD8  Treg  Mac  ...
  DE   DE    DE
    |
    v
Specificity scoring (AUROC per gene per cell type)
    |
    v
Surface protein filter (is_surface = True)
    |
    v
Clinical validation check (existing assay, clinical trial)
    |
    v
BiomarkerCandidate output
```

### 8.3 Biomarker Types

| Type | Definition | Example |
|------|-----------|---------|
| Diagnostic | Distinguishes disease from normal | CD19 for B-ALL detection |
| Prognostic | Predicts outcome regardless of treatment | TOX+ exhausted CD8 fraction predicts poor OS |
| Predictive | Predicts treatment response | PD-L1 on tumor cells predicts anti-PD-1 response |
| Pharmacodynamic | Measures treatment effect | CD8/Treg ratio change under immunotherapy |

---

## 9. CAR-T Target Validation

### 9.1 The Ideal CAR-T Target

| Property | Ideal | Acceptable | Unacceptable |
|----------|-------|-----------|--------------|
| On-tumor coverage | > 95% | > 70% | < 50% |
| Off-tumor vital organs | 0 hits | Low-level (< 0.5 TPM) | High in heart, brain, lung |
| Therapeutic index | > 10 | > 3 | < 3 |
| Heterogeneity | Low (uniform expression) | Moderate | High (bimodal) |
| Escape risk | Low (essential gene) | Medium | High (dispensable antigen) |

### 9.2 The Agent's Target Validation Pipeline

```
Target Gene (e.g., CD19, MSLN, HER2)
         |
    +----+----+
    |         |
    v         v
On-Tumor    Off-Tumor
Analysis    Safety Check
    |         |
    v         v
Coverage    8 vital organs:
percentage  brain, heart, lung,
Mean expr.  liver, kidney, pancreas,
            bone_marrow, intestine
    |         |
    +----+----+
         |
         v
Therapeutic Index = mean_on_tumor / (max_off_tumor + 0.01)
         |
         v
    +----+----+----+
    |         |    |
    v         v    v
FAVORABLE  COND.  UNFAVORABLE
(safe +    (risk   (safety or
 effective) mitig.) efficacy fail)
```

### 9.3 Safety Switch Integration

For CONDITIONAL targets, the agent recommends:
- **iCasp9 (inducible caspase 9):** Dimerizer-activated suicide switch
- **EGFRt:** Truncated EGFR enabling cetuximab-mediated depletion
- **Affinity-tuned CAR:** Reduced scFv affinity discriminates high-expression tumor from low-expression normal tissue

---

## 10. Advanced Study Resources

### 10.1 Key Papers

| Year | Paper | Impact |
|------|-------|--------|
| 2017 | Zheng et al., "Massively parallel digital transcriptional profiling" | 10x Chromium technology paper |
| 2018 | Wolf et al., "SCANPY: large-scale single-cell gene expression data analysis" | Standard Python toolkit |
| 2019 | Stuart et al., "Comprehensive Integration of Single-Cell Data" | Seurat v3, integration methods |
| 2020 | Bergen et al., "Generalizing RNA velocity" | RNA velocity dynamical model |
| 2021 | Stahl et al., "Visualization and analysis of gene expression in tissue sections by spatial transcriptomics" | Visium technology |
| 2023 | Theodoris et al., "Transfer learning enables predictions in network biology" | Geneformer foundation model |
| 2024 | Cui et al., "scGPT: toward building a foundation model for single-cell multi-omics" | scGPT foundation model |

### 10.2 Online Courses

- **Single Cell Genomics** (Wellcome Sanger Institute) -- Comprehensive bioinformatics training
- **Analysis of Single Cell RNA-seq Data** (Cambridge University) -- Scanpy/Seurat tutorials
- **NVIDIA RAPIDS for Single-Cell** -- GPU acceleration training

### 10.3 Practice Datasets

| Dataset | Cells | Tissue | Modality | Access |
|---------|-------|--------|----------|--------|
| PBMC 3K | 2,700 | Blood | scRNA-seq | 10x Genomics |
| PBMC 68K | 68,000 | Blood | scRNA-seq | 10x Genomics |
| Tabula Sapiens | 500,000 | Multi-tissue | scRNA-seq | CellxGene |
| Human Lung Cell Atlas | 580,000 | Lung | scRNA-seq + spatial | CellxGene |
| TCGA Pan-Cancer scRNA | 1M+ | Multi-cancer | scRNA-seq | TISCH2 |

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Learning Guide: Advanced Topics v1.0.0*
