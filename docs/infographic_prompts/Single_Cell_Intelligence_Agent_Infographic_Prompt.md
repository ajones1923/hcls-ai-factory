# Nano Banana Pro — Single-Cell Intelligence Agent

## RENDERING RULES
- Landscape 16:9, ONE canvas, white background
- NO text smaller than 12pt. Metric numbers: 18pt bold. Headers: 16pt bold.
- NO dense code blocks. Render everything as clean visual cards.
- Every character legible at 1920x1080

## STYLE
- Colors: NVIDIA Green #76B900, Navy #1B2333, Teal #1AAFCC, Red #DC2626 (hot TME), Blue #2196F3 (cold TME), Amber #F5A623 (excluded TME), Purple #7B2D8E (immunosuppressive TME), Emerald #059669, Orange #FF9800
- Rounded rectangles, thin-line icons, metric badge pills
- Reference architecture poster density

---

## TITLE BAR

**Left badges:** "Single-Cell Intelligence Agent" (green) | "57 Cell Types | 10 Workflows" (teal)

**Center:** Single-Cell Intelligence Agent — Cellular-Resolution Tumor Intelligence — Part of the HCLS AI Factory

**Right badge column:**
- 57 Cell Types
- 75 Marker Genes
- 30 Drugs (IC50)
- 4 TME Classifications
- 4 Decision Engines

---

## LEFT COLUMN: CELL TYPE ATLAS

**Header:** "Cell Type Atlas — 57 Types"

Three compartment cards with colored left borders:

**Immune** (blue): CD8+ T, CD4+ T, Tregs, NK, B Cells, Macrophages

**Stromal** (emerald): Fibroblasts/CAFs, Endothelial, Pericytes

**Tumor** (red): Malignant, Cycling, Hypoxic, Stem-like

Note: "+42 additional types in full atlas"

Arrow right: "12 Collections | BGE 384-dim"

---

## CENTER-TOP: MILVUS COLLECTIONS

**Header:** "Milvus 2.4 | 12 Collections"

Two rows of 6 colored pill badges:
- sc_cell_types, sc_markers, sc_spatial, sc_tme, sc_drug_response, sc_literature
- sc_methods, sc_datasets, sc_trajectories, sc_pathways, sc_clinical, genomic_evidence (shared)

---

## CENTER: TME CLASSIFICATION

**Header:** "4 TME Phenotypes — Immunotherapy Prediction"

Four large horizontal color cards:

| Class | Color | Key Feature | Therapy |
|-------|-------|-------------|---------|
| Hot-Inflamed | RED | High CD8+, IFN-gamma | Checkpoint inhibitors |
| Excluded | AMBER | Immune at margin only | Combination therapy |
| Cold-Desert | BLUE | No immune infiltration | CAR-T / vaccines |
| Immunosuppressive | PURPLE | Tregs, M2 dominant | Treg depletion + combo |

Small callout: "TMEClassifier — deterministic, no LLM"

---

## CENTER: WORKFLOWS + ENGINES

**Left: 10 Workflows** (numbered list, name only)

1. Cell Type Annotation
2. TME Profiling
3. Drug Response Prediction
4. Subclonal Architecture
5. Spatial Niche Mapping
6. Trajectory Inference
7. Ligand-Receptor Interaction
8. Biomarker Discovery
9. CAR-T Target Validation
10. Treatment Monitoring

**Right: 4 Decision Engines** (cards)

- TMEClassifier — immune/stromal scoring
- SubclonalRiskScorer — resistance prediction
- TargetExpressionValidator — CAR-T feasibility
- CellularDeconvolutionEngine — bulk to single-cell

---

## CENTER: TME ATLAS + DRUGS

**Left: Cancer TME Atlas (top 6)**

| Cancer | TME | Response |
|--------|-----|----------|
| Melanoma | Hot | High (40%) |
| NSCLC | Mixed | Moderate |
| Breast TNBC | Excluded | Emerging |
| Colorectal | Cold/Hot | MSI-H responds |
| Pancreatic | Cold | Very low |
| RCC | Hot | High (42%) |

**Right: Drug Sensitivity (top 6 of 30)**

| Drug | Target | IC50 |
|------|--------|------|
| Venetoclax | BCL2 | 10-100 nM |
| Bortezomib | Proteasome | 1-10 nM |
| Imatinib | BCR-ABL | 100-500 nM |
| Dabrafenib | BRAF | 0.5-5 nM |
| Osimertinib | EGFR | 5-15 nM |
| Vincristine | Tubulin | 1-10 nM |

---

## BOTTOM: PEDIATRIC + CROSS-AGENT

**Pediatric card** (emerald border):
- ALL Blast Typing: Pre-B (CD19+/CD10+), T-ALL (CD3+/CD7+)
- CAR-T Targets: CD19 >95%, CD22 80-90%, GD2 >90%
- MRD sensitivity: 10^-4

**Cross-Agent card** (teal border):
- Oncology (:8527), CAR-T (:8522), Biomarker (:8529)

---

## RIGHT COLUMN: OUTPUTS

Six output cards (icon + name):
1. TME Profile Report
2. Cell Type Annotation
3. Drug Response Prediction
4. CAR-T Validation
5. Integrated Assessment
6. Reports (MD/JSON/PDF/FHIR)

**UI:** Streamlit :8130 (5 tabs) | FastAPI :8540

---

## INFRASTRUCTURE BAR

Green full-width bar: NVIDIA DGX Spark | Milvus 2.4 | BGE 384-dim | Claude Sonnet 4 | Ports 8130, 8540

---

## FOOTER

Adam Jones | Apache 2.0 | hcls-ai-factory.org | March 2026
