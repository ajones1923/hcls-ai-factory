# Nano Banana Pro — Precision Biomarker Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Precision Biomarker Agent — the most analytically deep agent in the HCLS AI Factory, featuring 14 Milvus collections, 14 specialized analysis modules (including PhenoAge biological age calculation, genotype-adjusted lab interpretation, pharmacogenomics profiling, and critical value detection), and an 8-tab Streamlit UI. It connects to 4 peer agents and is the platform's primary biomarker intelligence engine for both adult wellness and pediatric oncology (CD19/CD22 expression, MRD detection, MYCN amplification). Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual language as the Imaging, CAR-T, and Oncology Agent infographics in this series — data sources left, processing center, outputs right, cross-agent strip at bottom.

**Canvas:** Dense but organized. Every section carries information. The diagram should feel like a clinical laboratory reference poster — comprehensive enough for a clinical pathologist, informative enough for a precision medicine clinician, clean enough for an executive.

**Typography:**
- Title: Large, bold, Inter/Helvetica, navy (#1B2333)
- Subtitle: Teal (#1AAFCC), smaller
- Section headers: Bold, navy, green left-border accent
- Sub-headers: Bold, teal (#1AAFCC)
- Body: Dark gray (#333333), 8-10pt, SHORT PHRASES only
- Metric badges: White text on colored pills
- All text MUST be legible at 1920x1080

**Color Palette:**
- NVIDIA Green: #76B900 — badges, infrastructure bar
- Navy: #1B2333 — title, headers
- Teal: #1AAFCC — sub-headers, data flow, agent connections
- Emerald: #059669 — normal/healthy ranges, favorable findings
- Red: #DC2626 — critical values, high-risk alerts
- Amber: #F5A623 — moderate-risk, watch indicators
- Purple: #7B2D8E — genomic/PGx collections
- Blue: #2196F3 — clinical evidence collections
- Orange: #FF9800 — nutrition/drug interaction collections
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome)
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for critical alerts
- Color-coded collection and module badges
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Precision Biomarker Agent" (green #76B900, white text)
- "14 Analysis Modules" (teal #1AAFCC, white text)

**Center:**
- Line 1 (large, bold, navy): **Precision Biomarker Agent**
- Line 2 (gray): "Multi-Omic Biomarker Intelligence on NVIDIA DGX Spark"
- Line 3 (teal): "Part of the Precision Intelligence Network — HCLS AI Factory"
- Line 4 (small, gray): "GB10 Superchip | 128 GB Unified Memory | 14 Collections | 14 Analysis Modules | 8-Tab UI"

**Right side — Capability summary badge:**
```
ANALYTICAL DEPTH
━━━━━━━━━━━━━━━━━━━
14 Milvus Collections
14 Analysis Modules
PhenoAge Biological Age
Genotype-Adjusted Ranges
Cross-Modal Triggers to 8 Agents
Demo 2 + Demo 3
```

---

### ━━━ LEFT COLUMN: DATA DOMAINS & KNOWLEDGE SOURCES ━━━

**Stacked input cards representing the biomarker knowledge domains:**

1. **Reference Ranges** [chart icon]
   Age/sex-stratified normals
   Genotype-adjusted thresholds
   CPIC March 2025 guidelines

2. **Genetic Variants** [DNA icon]
   Pharmacogenes (CYP2D6, CYP2C19, TPMT, DPYD)
   Disease-risk alleles (APOL1, MTHFR, PNPLA3)
   Carrier screening panels

3. **PGx Rules** [pill+DNA icon]
   Drug-gene interaction rules
   Metabolizer phenotype mapping
   Dosing adjustment algorithms

4. **Disease Trajectories** [trend icon]
   Diabetes progression (TCF7L2, HbA1c)
   Liver disease (PNPLA3, ALT)
   Thyroid (DIO2, TSH)
   Cardiac risk (Lp(a), LDL)

5. **Clinical Evidence** [journal icon]
   Published biomarker studies
   Validation cohort data
   Guideline references (ADA, ESC, AASLD)

6. **Aging Markers** [clock icon]
   PhenoAge coefficients (Levine 2018)
   GrimAge surrogates (Lu 2019)
   Epigenetic clock biomarkers

7. **Critical Values** [alert icon]
   Life-threatening thresholds
   Immediate notification triggers
   Panic value protocols

8. **Pediatric Oncology** [child icon] [NEW]
   CD19/CD22 expression (B-ALL)
   MRD detection (flow/PCR)
   MYCN amplification (neuroblastoma)
   Tumor biomarker panels

Arrows flow rightward labeled: "BGE-small-en-v1.5 | 384-dim | Milvus Ingest"

---

### ━━━ CENTER-TOP: MILVUS VECTOR DATABASE — 14 COLLECTIONS ━━━

```
Milvus 2.4 | Vector Database — 14 Collections | IVF_FLAT / COSINE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

**Collection badges in colored rectangles (3 rows):**

Row 1 (core clinical):
[biomarker_reference] [biomarker_genetic_variants] [biomarker_pgx_rules] [biomarker_disease_trajectories] [biomarker_clinical_evidence]
   (blue)                (purple)                      (purple)               (blue)                          (blue)

Row 2 (specialized):
[biomarker_nutrition] [biomarker_drug_interactions] [biomarker_aging_markers] [biomarker_genotype_adjustments] [biomarker_monitoring]
   (orange)              (orange)                      (teal)                    (purple)                         (green)

Row 3 (safety + shared):
[biomarker_critical_values] [biomarker_discordance_rules] [biomarker_aj_carrier_screening] [genomic_evidence]
   (red)                       (amber)                        (purple)                         (emerald — shared)

**Note:** 14 collections — the most of any agent in the platform.

---

### ━━━ CENTER: 14 ANALYSIS MODULES (the unique depth of this agent) ━━━

**Section header (navy, green underline):** "14 Specialized Analysis Modules — Pure Computation, No LLM"

**Visual: Grid of 14 module cards arranged in a structured layout (7×2 or similar). Each card is a small rounded rectangle with an icon, module name, and 1-2 line description.**

**Row 1 — Core Analysis:**

```
[chart] Lab Range          [DNA] Genotype         [age] Biological        [trend] Disease
        Interpreter              Adjustment              Age (PhenoAge)          Trajectory

Age/sex reference     Genotype-specific     Levine 2018 algorithm   Longitudinal risk
ranges per analyte    threshold adjustment  9 blood biomarkers      prediction models
                      TCF7L2, PNPLA3,       → biological vs         Diabetes, liver,
                      DIO2, APOL1, MTHFR    chronological age       cardiac, thyroid
```

**Row 2 — Pharmacogenomics & Safety:**

```
[pill] Pharmaco-          [alert] Critical        [puzzle] Discordance    [shield] Audit
       genomics                   Values                   Detector               Trail

CYP2D6, CYP2C19,    Life-threatening       Conflicting results    Complete provenance
TPMT, DPYD, UGT1A1  thresholds             between related        Model ID + version
Metabolizer phenotype Panic value alerts    biomarkers flagged     Immutable logging
Star allele → dosing                        (e.g., low iron +
                                            high ferritin)
```

**Row 3 — Advanced & Export:**

```
[network] Cross-Modal     [document] Report        [translate] Translation   [export] Export
          Triggers                   Generator                                        Engine

8 trigger rules to    Markdown, JSON, PDF    Clinical → patient-     FHIR R4
peer agents           FHIR R4                friendly language       PDF branded
(imaging, oncology,   Phenopacket            Bilingual support       Phenopacket v2
genomics, cardio)     Branded NVIDIA theme                           CSV data
```

**Row 4 — Specialized:**

```
[dna] Carrier              [nutrition] Nutrition
      Screening                        Analysis

AJ panel (BRCA1/2,    Nutrient-biomarker
Tay-Sachs, CF,        interactions
Gaucher, etc.)         Diet modification
                       recommendations
```

**Key callout (emerald border, prominent):**
```
DETERMINISTIC COMPUTATION
Every module uses validated algorithms
No LLM involvement in calculations
PhenoAge: Levine et al. 2018 coefficients
GrimAge: Lu et al. 2019 surrogates
CPIC: March 2025 guideline versions
```

---

### ━━━ CENTER-BOTTOM: RAG + CROSS-AGENT ━━━

**Claude LLM box:**
```
Claude Sonnet 4
Anthropic API | Streaming RAG
Evidence synthesis + natural language explanation
Interprets module outputs — NEVER computes biomarker values
```

**Cross-Agent Integration box (teal border):**
```
Cross-Agent Integration (/biomarker/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Oncology Agent (:8527)      Correlate biomarkers with tumor molecular profile
→ CAR-T Agent (:8522)         Validate CAR-T target suitability from expression data
→ PGx Agent (:8107)           Pharmacogenomic implications of biomarker findings
→ Clinical Trial Agent (:8538) Match biomarker-driven clinical trials

8 Cross-Modal Triggers:
elevated_lpa → Cardiology imaging | high_psa → Oncology workup | low_egfr → Nephrology
abnormal_cbc → Hematology | elevated_afp → Hepatology/Oncology | thyroid_panel → Endocrine
mrd_positive → Oncology intensification | cd19_loss → CAR-T CD22 backup
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

**Output cards:**

1. **Biomarker Analysis Report** [document icon]
   Complete panel interpretation
   Genotype-adjusted ranges highlighted
   Critical values flagged in red

2. **Biological Age Assessment** [clock icon]
   PhenoAge calculation
   Chronological vs biological
   Aging acceleration/deceleration rate
   Health trajectory projection

3. **PGx Profile** [pill+DNA icon]
   Metabolizer phenotype per gene
   Drug dosing recommendations
   CPIC guideline citations

4. **Disease Risk Projections** [trend icon]
   Multi-domain risk assessment
   Longitudinal trajectories
   Intervention recommendations

5. **Integrated Assessment** [network icon]
   /biomarker/integrated-assessment
   4-agent coordinated output

6. **Cross-Modal Alerts** [alert icon]
   Automatic triggers to peer agents
   8 clinical escalation pathways

7. **Reports** [export icon]
   Markdown | JSON | PDF | FHIR R4
   Phenopacket v2

**Streamlit UI (prominent box):**
```
Streamlit UI (:8528) — 8 Tabs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. Biomarker Analysis (full panel)
2. Biological Age (PhenoAge calculator)
3. Disease Risk (trajectory analysis)
4. PGx Profile (drug-gene mapping)
5. Evidence Explorer (RAG Q&A)
6. Reports (multi-format export)
7. Patient 360 (cross-agent dashboard)
8. Longitudinal (trend tracking)

Demo patient quick-load button
```

**FastAPI REST:**
```
Port :8529
/analyze | /biological-age | /disease-risk
/pgx | /query | /query/stream
/biomarker/integrated-assessment
```

---

### ━━━ BOTTOM STRIP: DEMO PARTICIPATION + CROSS-MODAL ━━━

**Header (navy background, white text):** "Precision Biomarker Agent in HCLS AI Factory"

**Left half — Demo participation:**
```
Demo 2 ✓ "The 30-Second Tumor Board"
  → CD19/CD22 expression quantification
  → MRD kinetics assessment
  → Immunotherapy eligibility panel

Demo 3 ✓ "Protecting the Survivor"
  → MYCN amplification confirmation
  → LDH, ferritin, VMA/HVA panel
  → COG/INRG risk classification
```

**Right half — Cross-Modal Trigger Map (visual showing 8 trigger paths):**
```
Biomarker Finding          →    Triggered Agent          Action
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Lp(a) > 125 nmol/L        →    Cardiology/Imaging       CCTA referral
PSA > 4.0 ng/mL           →    Oncology                 Prostate workup
eGFR < 60 mL/min          →    Nephrology               CKD staging
CBC abnormalities          →    Hematology               Smear review
AFP > 400 ng/mL           →    Oncology/Imaging         HCC surveillance
Thyroid panel abnormal     →    Endocrine                Thyroid workup
MRD ≥ 0.01%               →    Oncology                 Treatment intensification
CD19 loss post-CAR-T       →    CAR-T Agent              CD22 backup strategy
```

---

### ━━━ INFRASTRUCTURE BAR ━━━

**Green (#76B900) background, white text:**

```
DGX Spark (GB10)    Milvus 2.4           BGE-small-en-v1.5    Claude Sonnet 4     Docker Compose
$4,699              14 collections        384-dim vectors       Streaming RAG       Ports 8528, 8529
128 GB unified      IVF_FLAT/COSINE      Async embedding       + Biomarker prompt  14 analysis modules
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 Open Source | hcls-ai-factory.org | HCLS AI Factory v2.0 | March 2026
PhenoAge (Levine 2018) | GrimAge (Lu 2019) | CPIC March 2025 | ADA 2024 | ESC/EAS 2021 | AASLD 2023
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **14 collections — the most of any agent.** This is the deepest analytical agent in the platform. The collection count must be visually prominent.
2. **14 analysis modules — all deterministic.** PhenoAge, genotype adjustment, critical values, discordance detection, PGx — none use LLM for computation. The "DETERMINISTIC COMPUTATION" callout must be prominent.
3. **PhenoAge biological age** is a unique capability — calculating biological vs chronological age from 9 blood biomarkers. This is a headline feature.
4. **Genotype-adjusted lab interpretation** — TCF7L2 adjusts HbA1c thresholds, PNPLA3 adjusts ALT thresholds, APOL1 adjusts eGFR. This is precision medicine applied to routine lab work.
5. **8 cross-modal triggers** — biomarker findings automatically escalate to peer agents (Cardiology, Oncology, etc.). The trigger map must be visible.
6. **Pediatric oncology biomarkers** — CD19/CD22 expression, MRD detection, MYCN amplification. Not just adult wellness.
7. **8-tab Streamlit UI** — the most tabs of any agent. Full patient analysis from a single interface.
8. **4 peer agent connections** — Oncology, CAR-T, PGx, Clinical Trial via /integrated-assessment.
9. **Participates in Demo 2 + Demo 3** — tumor board biomarkers and neuroblastoma risk classification.
10. **Dense enough for a clinical laboratory wall, clear enough to present to a chief medical officer.**
