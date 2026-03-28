# Single-Cell Intelligence Agent -- Demo Guide

**Date:** March 2026
**Author:** Adam Jones

---

## Prerequisites

Before running any demo, ensure the agent is deployed and healthy:

```bash
curl -s http://localhost:8540/health | python -m json.tool
# Verify: "status": "healthy"
```

Access the Streamlit UI at `http://localhost:8130` or use the API directly via `curl` or the Swagger UI at `http://localhost:8540/docs`.

---

## Demo 1: Cell Type Annotation

### Scenario

A researcher has performed scRNA-seq on a PBMC sample from a melanoma patient and needs to identify the immune cell populations present.

### Streamlit UI

1. Open the **Chat** tab
2. Enter the query:

> "I have a PBMC scRNA-seq dataset from a melanoma patient. The top clusters show these marker genes: Cluster 0: CD3D, CD3E, CD8A, GZMB, PRF1; Cluster 1: CD14, LYZ, S100A9, VCAN, FCN1; Cluster 2: CD19, MS4A1, CD79A, PAX5; Cluster 3: NCAM1, NKG7, GNLY, KLRD1; Cluster 4: FOXP3, IL2RA, CTLA4, CD4. What are these cell types and what is their clinical significance?"

3. Review the response:
   - Each cluster annotated with cell type and Cell Ontology ID
   - Marker gene evidence for each annotation
   - Clinical significance in melanoma context

### API Call

```bash
curl -X POST http://localhost:8540/v1/sc/query \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "query": "Identify cell types from these PBMC cluster markers: CD3D/CD3E/CD8A/GZMB/PRF1, CD14/LYZ/S100A9, CD19/MS4A1/CD79A, NCAM1/NKG7/GNLY, FOXP3/IL2RA/CTLA4",
    "tissue_type": "PBMC",
    "disease_context": "melanoma",
    "workflow_type": "cell_type_annotation"
  }'
```

### Expected Output

| Cluster | Cell Type | CL ID | Confidence | Key Markers |
|---------|----------|-------|-----------|-------------|
| 0 | CD8+ cytotoxic T cell | CL:0000625 | High | CD8A, GZMB, PRF1 |
| 1 | Classical monocyte | CL:0000576 | High | CD14, LYZ, S100A9 |
| 2 | B cell | CL:0000236 | High | CD19, MS4A1, CD79A |
| 3 | Natural killer cell | CL:0000623 | High | NCAM1, NKG7, GNLY |
| 4 | Regulatory T cell | CL:0000815 | High | FOXP3, IL2RA, CTLA4 |

### Talking Points

- The agent maps marker gene signatures to canonical cell types using the 44-entry Cell Type Atlas
- Cell Ontology (CL) identifiers enable interoperability with Human Cell Atlas references
- Clinical significance contextualized to melanoma (e.g., Treg fraction relates to immunotherapy resistance)

---

## Demo 2: TME Profiling

### Scenario

An oncologist wants to classify the tumor microenvironment of a lung adenocarcinoma biopsy to guide immunotherapy selection.

### Streamlit UI

1. Open the **TME Profiler** tab
2. Enter cell type proportions:
   - CD8_T: 0.18
   - CD4_T: 0.12
   - Treg: 0.04
   - NK: 0.03
   - B_cell: 0.05
   - Macrophage_M1: 0.06
   - Macrophage_M2: 0.08
   - Fibroblast: 0.15
   - Epithelial: 0.25
   - Endothelial: 0.04

3. Set PD-L1 TPS to 55 (optional)
4. Click **Classify TME**

### API Call

```bash
curl -X POST http://localhost:8540/v1/sc/tme-profile \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "cell_type_proportions": {
      "CD8_T": 0.18, "CD4_T": 0.12, "Treg": 0.04,
      "NK": 0.03, "B_cell": 0.05, "Macrophage_M1": 0.06,
      "Macrophage_M2": 0.08, "Fibroblast": 0.15,
      "Epithelial": 0.25, "Endothelial": 0.04
    },
    "gene_expression": {
      "CD274": 2.5, "CTLA4": 1.8, "LAG3": 1.2,
      "IDO1": 0.5, "TGFB1": 1.1
    },
    "pdl1_tps": 55,
    "cancer_type": "lung_adenocarcinoma"
  }'
```

### Expected Output

```json
{
  "tme_class": "hot_inflamed",
  "immune_score": 0.48,
  "stromal_score": 0.15,
  "pdl1_status": "high",
  "treatment_recommendations": [
    "Immune-hot TME: strong candidate for checkpoint inhibitor therapy (anti-PD-1/PD-L1)",
    "PD-L1 TPS >= 50%: consider first-line pembrolizumab monotherapy",
    "LAG-3 co-expressed: consider relatlimab + nivolumab (Opdualag)"
  ],
  "evidence_level": "strong",
  "severity": "informational"
}
```

### Talking Points

- Classification uses deterministic TME Classifier engine (not LLM)
- PD-L1 TPS combined with scRNA-seq data provides stronger evidence than either alone
- LAG-3 co-expression detected, enabling dual-checkpoint recommendation
- Immune score of 0.48 exceeds the 0.25 threshold for hot-inflamed classification

---

## Demo 3: Drug Response Prediction

### Scenario

A clinical team needs to understand which drugs a triple-negative breast cancer (TNBC) patient is likely to respond to, based on single-cell profiling.

### Streamlit UI

1. Open the **Chat** tab
2. Enter the query:

> "A TNBC patient has scRNA-seq showing 40% tumor cells with high BRCA1 pathway activity, 15% exhausted CD8+ T cells, and 20% M2 macrophages. What drugs would you recommend and what resistance mechanisms should we watch for?"

### API Call

```bash
curl -X POST http://localhost:8540/v1/sc/drug-response \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "query": "Drug response prediction for TNBC with BRCA1 activity, exhausted CD8+ T cells, and M2 macrophage infiltration",
    "disease_context": "triple-negative breast cancer",
    "cell_types_of_interest": ["tumor", "CD8_T_exhausted", "Macrophage_M2"],
    "genes_of_interest": ["BRCA1", "PDCD1", "LAG3", "HAVCR2", "CD163"]
  }'
```

### Expected Output Highlights

- **PARP inhibitor (Olaparib):** High sensitivity predicted based on BRCA1 pathway activity; watch for BRCA1 reversion mutations in expanding subclones
- **Checkpoint inhibitor (Atezolizumab):** Moderate sensitivity; exhausted T cells may benefit from anti-PD-L1 but M2 macrophage suppression limits efficacy
- **CSF1R inhibitor combination:** Recommended to reprogram M2 macrophages and enhance checkpoint blockade

### Talking Points

- Drug response prediction integrates cell type composition with drug mechanism databases
- Resistance mechanisms identified at the subpopulation level
- Combination strategies emerge from multi-compartment analysis

---

## Demo 4: CAR-T Target Validation

### Scenario

A cell therapy team is evaluating CD19 as a target for CAR-T therapy in a B-cell lymphoma patient and wants to assess safety and efficacy.

### Streamlit UI

1. Open the **Workflows** tab
2. Select **CAR-T Target Validation**
3. Enter:
   - Target gene: CD19
   - Cancer type: B-cell lymphoma

### API Call

```bash
curl -X POST http://localhost:8540/v1/sc/cart-validate \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "target_gene": "CD19",
    "tumor_expression_values": [3.2, 4.1, 2.8, 3.5, 4.0, 3.8, 3.1, 4.2, 0.05, 0.02, 3.9, 4.1, 3.7, 3.3, 2.9, 4.0, 3.6, 3.8, 0.08, 4.1],
    "normal_tissue_expression": {
      "brain": 0.01, "heart": 0.02, "lung": 0.15,
      "liver": 0.05, "kidney": 0.03, "pancreas": 0.01,
      "bone_marrow": 1.8, "intestine": 0.08
    },
    "cancer_type": "B-cell lymphoma"
  }'
```

### Expected Output

```json
{
  "target_gene": "CD19",
  "on_tumor_pct": 0.85,
  "mean_tumor_expression": 3.14,
  "off_tumor_hits": {
    "bone_marrow": 1.8
  },
  "therapeutic_index": 1.74,
  "safety_assessment": "moderate_risk",
  "safety_detail": "Low-level expression in: bone_marrow. Monitor for on-target off-tumour toxicity",
  "efficacy_assessment": "adequate",
  "overall_verdict": "CONDITIONAL",
  "recommendations": [
    "CD19: viable with risk mitigation. Consider dose escalation protocol and enhanced monitoring",
    "Incorporate safety switch (iCasp9 or EGFRt) in CAR design",
    "Consider tandem or dual-target CAR to improve coverage"
  ]
}
```

### Talking Points

- 85% on-tumor coverage is adequate but 15% antigen-negative cells pose escape risk
- Bone marrow expression (1.8) creates on-target off-tumor toxicity concern (expected: B cell aplasia)
- Therapeutic index of 1.74 is below the 3.0 "acceptable" threshold
- Safety switch (iCasp9) recommended as a mitigation strategy
- This matches the known clinical profile of CD19 CAR-T (B cell aplasia is an accepted on-target toxicity)

---

## Demo 5: Spatial Analysis

### Scenario

A researcher has Visium spatial transcriptomics data from a colorectal cancer biopsy and wants to identify tissue architecture patterns.

### Streamlit UI

1. Open the **Chat** tab
2. Enter the query:

> "In our Visium spatial transcriptomics data from a colorectal cancer sample, we see three distinct spatial niches: (1) an immune-rich region at the tumor invasive margin with high CD8A, GZMB, and CXCL9 expression; (2) a central tumor core with high MKI67, TOP2A, and low immune markers; (3) a fibrotic stroma with COL1A1, FAP, and TGFB1. What is the clinical significance of these spatial patterns?"

### API Call

```bash
curl -X POST http://localhost:8540/v1/sc/spatial-niche \
  -H "Content-Type: application/json" \
  -H "X-API-Key: your-api-key" \
  -d '{
    "query": "Spatial niche analysis of colorectal cancer with immune margin, proliferative core, and fibrotic stroma",
    "tissue_type": "colorectal",
    "disease_context": "colorectal cancer",
    "spatial_platform": "visium",
    "genes_of_interest": ["CD8A", "GZMB", "CXCL9", "MKI67", "TOP2A", "COL1A1", "FAP", "TGFB1"]
  }'
```

### Expected Output Highlights

- **Niche 1 (Immune margin):** Tertiary lymphoid structure-like organization; positive prognostic indicator; suggests immune-excluded TME with potential for anti-TGFb to enable infiltration
- **Niche 2 (Proliferative core):** High cycling tumor with no immune infiltrate; cold desert phenotype in this region
- **Niche 3 (Fibrotic stroma):** CAF-rich barrier with TGFb signaling; primary mechanism of immune exclusion

### Talking Points

- Spatial context reveals that this tumor is EXCLUDED (not COLD) -- immune cells are present but confined to the margin
- The fibrotic stroma (FAP+, TGFb1+) physically excludes T cells from the tumor core
- Anti-TGFb combination (e.g., bintrafusp alfa) could enable T cell infiltration
- This spatial architecture is invisible in dissociated scRNA-seq and bulk sequencing

---

## Demo Flow Recommendations

### 5-Minute Quick Demo

Run Demo 2 (TME Profiling) via the Streamlit TME Profiler tab. It is the most visual and clinically impactful demonstration.

### 15-Minute Executive Demo

1. Demo 1 (Cell Type Annotation) -- 3 minutes
2. Demo 2 (TME Profiling) -- 5 minutes
3. Demo 4 (CAR-T Validation) -- 5 minutes
4. Health dashboard overview -- 2 minutes

### 30-Minute Technical Deep Dive

Run all 5 demos sequentially, with discussion of:
- Multi-collection RAG architecture
- Decision support engine determinism
- Cross-agent integration points
- GPU acceleration roadmap

### Post-Demo API Exploration

Direct attendees to `http://localhost:8540/docs` for interactive Swagger documentation where they can construct custom queries against any of the 25 API endpoints.

---

## Troubleshooting Demo Issues

| Issue | Fix |
|-------|-----|
| "degraded" health | Wait 30s for Milvus startup; run `docker compose ps` |
| Empty search results | Run `docker compose run sc-setup` to seed data |
| LLM timeout | Check ANTHROPIC_API_KEY; Claude API may be rate-limited |
| Streamlit not loading | Check port 8130 is not blocked; try `docker compose restart sc-streamlit` |
| Slow first query | Normal -- first query loads embedding model and Milvus collections into memory |

---

*HCLS AI Factory -- Single-Cell Intelligence Agent Demo Guide*

---

!!! warning "Clinical Decision Support Disclaimer"
    The Single-Cell Intelligence Agent is a clinical decision support research tool for single-cell genomic analysis. It is not FDA-cleared and is not intended as a standalone diagnostic device. All recommendations should be reviewed by qualified healthcare professionals. Apache 2.0 License.
