# Nano Banana Pro — Rare Disease Diagnostic Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Rare Disease Diagnostic Agent — the HCLS AI Factory's most diagnostically complex agent, covering 88 rare diseases across 13 categories, ACMG/AMP variant classification with 23 of the 28 criteria, 9 diagnostic algorithms, gene therapies, HPO phenotype matching, 14 Milvus collections, and 6 decision support engines. It connects to 4 peer agents and is the centerpiece of Demo 4 ("One Gene, One Family" — Aurora's retinoblastoma triggering sibling cascade testing). Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual series as all other HCLS AI Factory agent infographics.

**Canvas:** Dense but organized. The diagram should feel like a clinical genetics department reference poster — comprehensive enough for a medical geneticist, precise enough for a molecular diagnostics lab director, human enough for the family receiving the diagnosis.

**Typography — READABILITY IS CRITICAL:**
- Title: 24pt bold, Inter/Helvetica, navy (#1B2333)
- Subtitle: 16pt, teal (#1AAFCC)
- Section headers: 14pt bold, navy, green left-border accent
- Sub-headers: 12pt bold, teal (#1AAFCC)
- Body text: 10pt minimum everywhere — dark gray (#333333), SHORT PHRASES only
- Metric badges: 14pt minimum, white text on colored pills, LARGE AND PROMINENT
- ACMG classifications: 14pt bold color-coded blocks (red/amber/emerald)
- NO TEXT SMALLER THAN 10pt ANYWHERE ON THE CANVAS
- All text MUST be legible at 1920x1080
- RENDER ALL STRUCTURED DATA AS VISUAL CARDS AT 10pt MINIMUM — no raw code blocks

**Color Palette:**
- NVIDIA Green: #76B900 — badges, infrastructure bar
- Navy: #1B2333 — title, headers
- Teal: #1AAFCC — sub-headers, data flow, agent connections
- Red: #DC2626 — Pathogenic classification, critical findings
- Amber: #F5A623 — VUS, uncertain significance
- Emerald: #059669 — Benign, favorable, gene therapy available
- Purple: #7B2D8E — genetic/genomic collections, rare disease categories
- Blue: #2196F3 — phenotype/HPO collections
- Orange: #FF9800 — metabolic diseases, newborn screening
- Pink: #EC4899 — family cascade, genetic counseling
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome) — DNA, family tree, puzzle piece, magnifying glass motifs
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for pathogenic alerts
- Color-coded ACMG classification blocks (5-tier)
- Metric badge pills: LARGE TEXT, 14pt minimum, high contrast

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges (LARGE TEXT, 14pt minimum):
- "Rare Disease Diagnostic Agent" (green #76B900, white text)
- "88 Diseases | 23 ACMG Criteria" (purple #7B2D8E, white text)

**Center:**
- Line 1 (24pt bold, navy): **Rare Disease Diagnostic Agent**
- Line 2 (14pt, gray): "Phenotype-Driven Diagnosis & ACMG Variant Classification"
- Line 3 (12pt, teal): "Part of the Precision Intelligence Network — HCLS AI Factory"

**Right side — Capability badge (RENDER AS VISUAL CARD AT 10pt MINIMUM):**
Seven metric pills, each 14pt bold white text on colored background:
- "88 Rare Diseases" (purple pill)
- "23 ACMG Criteria" (red pill)
- "9 Algorithms" (navy pill)
- "6 Engines" (teal pill)
- "12 Gene Therapies" (emerald pill)
- "Phenopacket v2" (blue pill)
- "Demo 4" (pink pill)

---

### ━━━ LEFT COLUMN: DISEASE CATEGORIES + DATA SOURCES ━━━

**Section header (14pt bold):** "13 Disease Categories — 88 Diseases"

**RENDER AS VISUAL CARDS AT 10pt MINIMUM — show category name and disease count only, no individual disease names:**

| Color Border | Category | Count |
|---|---|---|
| Orange | Metabolic | 28 |
| Purple | Neurological | 23 |
| Red | Hematologic | 15 |
| Blue | Immunologic | 13 |
| Teal | Connective Tissue | 10 |
| Pink | Cancer Predisposition | 8 |
| Gray | Skeletal Dysplasias | — |
| Gray | Renal | — |
| Gray | Cardiac | — |
| Gray | Pulmonary | — |
| Gray | Endocrine | — |
| Gray | Dermatologic | — |
| Gray | Multi-System | — |

Categories 7-13 render as a single compact row of labels beneath the six primary cards.

**Data Sources (5 compact pills beneath categories, 10pt):**
HPO | OMIM/Orphanet | ClinVar | gnomAD | HGMD

Arrow rightward: "14 Collections | BGE-small-en-v1.5 | 384-dim"

---

### ━━━ CENTER-TOP: MILVUS COLLECTIONS — 14 ━━━

**Section header (14pt bold):** "Milvus 2.4 — 14 Collections — IVF_FLAT / COSINE"

**Collection badges in 3 rows (each badge is a colored rounded pill, 10pt label):**

Row 1 — Clinical:
[rd_phenotypes] blue | [rd_diseases] purple | [rd_genes] purple | [rd_variants] red | [rd_literature] blue

Row 2 — Therapeutic:
[rd_trials] blue | [rd_therapies] emerald | [rd_case_reports] amber | [rd_guidelines] purple | [rd_pathways] orange

Row 3 — Specialized:
[rd_registries] teal | [rd_natural_history] teal | [rd_newborn_screening] orange | [genomic_evidence] emerald (shared)

---

### ━━━ CENTER: 6 DECISION SUPPORT ENGINES ━━━

**Section header (14pt bold, navy, green underline):** "6 Decision Support Engines — Deterministic Clinical Logic"

**RENDER AS 6 VISUAL CARDS IN A 3x2 GRID AT 10pt MINIMUM. Each card: icon + name (12pt bold) + single-line description (10pt):**

Row 1:
- [puzzle icon] **HPO-to-Gene Matcher** — Semantic similarity scoring, HPO terms to ranked candidate genes
- [checklist icon] **ACMG Variant Classifier** — 23 criteria, point-based 5-tier pathogenicity classification
- [pill icon] **Orphan Drug Matcher** — FDA orphan drug + gene therapy eligibility matching

Row 2:
- [algorithm icon] **Algorithm Recommender** — 9 algorithms, phenotype cluster to recommended testing pathway
- [family icon] **Family Segregation Analyzer** — Inheritance pattern analysis (AD/AR/XL/XR/MT), co-segregation scoring
- [timeline icon] **Natural History Predictor** — Disease progression modeling, age-of-onset, treatment windows

---

### ━━━ ACMG 5-TIER CLASSIFICATION (prominent center callout) ━━━

**RENDER AS 5 LARGE HORIZONTAL COLOR BLOCKS, each with 14pt bold text. This must be one of the most visually prominent elements on the canvas:**

- RED block (#DC2626): **Pathogenic** — 10+ points — Clinical action required
- RED block (lighter): **Likely Pathogenic** — 6-9 points — Clinical action recommended
- AMBER block (#F5A623): **VUS** — -1 to 5 points — Insufficient evidence
- EMERALD block (lighter): **Likely Benign** — -6 to -2 points — Low significance
- EMERALD block (#059669): **Benign** — -7 or below — No clinical significance

Below the blocks, single line (10pt): "23 of 28 ACMG/AMP evidence codes implemented (PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7)"

---

### ━━━ CENTER-MIDDLE: 9 DIAGNOSTIC ALGORITHMS + GENE THERAPY ━━━

**Left half — 9 Diagnostic Algorithms (RENDER AS VISUAL CARDS AT 10pt MINIMUM, name only with 1-line description):**

**Section header (14pt bold):** "9 Diagnostic Algorithms"

1. **Phenotype-Driven Diagnosis** — HPO terms to ranked differentials
2. **WES/WGS Variant Interpretation** — VCF to ACMG classification
3. **Metabolic Screening** — Newborn screening to IEM diagnosis
4. **Dysmorphology Assessment** — Features to syndrome matching
5. **Neurogenetic Evaluation** — Seizures/regression gene panel pathway
6. **Cardiac Genetics** — Cardiomyopathy/arrhythmia cascade testing
7. **Connective Tissue** — Marfan/EDS criteria to gene panel
8. **IEM Emergency** — Acute metabolic crisis rapid diagnosis
9. **Gene Therapy Eligibility** — Confirmed diagnosis to therapy/trial match

**Right half — Top 6 Approved Gene Therapies (RENDER AS VISUAL CARDS AT 10pt MINIMUM):**

**Section header (14pt bold):** "Gene Therapies (6 of 12 Approved)"

| Therapy | Disease | Gene | Year |
|---|---|---|---|
| Zolgensma | SMA Type 1 | SMN1 | 2019 |
| Luxturna | Leber Congenital | RPE65 | 2017 |
| Casgevy | Sickle Cell/Thal | BCL11A | 2023 |
| Elevidys | Duchenne MD | DMD | 2023 |
| Hemgenix | Hemophilia B | F9 | 2022 |
| Skysona | CALD | ABCD1 | 2022 |

Below table, single line (10pt): "+6 more approved therapies tracked | Eligibility check integrated"

---

### ━━━ CENTER-BOTTOM: PEDIATRIC + CROSS-AGENT ━━━

**Pediatric / Cancer Predisposition box (pink #EC4899 border, RENDER AS VISUAL CARD AT 10pt MINIMUM):**

**Header (12pt bold):** "Hereditary Cancer Predisposition in Children"

Four lines only:
- **Li-Fraumeni (TP53)** — Lifetime cancer risk >90%, childhood adrenocortical carcinoma
- **Retinoblastoma (RB1)** — Demo 4: Aurora. Bilateral = germline. Sibling 50% risk. Cascade testing CRITICAL
- **Lynch (MLH1/MSH2/MSH6/PMS2)** — Lifetime CRC risk 40-80%
- **FAP (APC)** — Polyps by teens, CRC by 40 without colectomy

Single line beneath (10pt): "Cascade: Proband > Parents > Siblings > Extended Family | PGT-M available"

**Cross-Agent Integration box (teal border, RENDER AS VISUAL CARD AT 10pt MINIMUM):**

**Header (12pt bold):** "Cross-Agent Integration"

Four compact lines:
- Cardiology Agent (:8126) — Inherited cardiac conditions
- Biomarker Agent (:8529) — Metabolic profiles for IEM
- PGx Agent (:8107) — Post-diagnosis dosing optimization
- Imaging Agent (:8524) — Trilateral screening, skeletal surveys

Single line beneath (10pt): "4 peer agents | Graceful degradation | Pediatric-first"

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

**6 output cards (RENDER AS VISUAL CARDS AT 10pt MINIMUM, each card: icon + name + 1-line description):**

1. [ranked list icon] **Differential Diagnosis** — Phenotype-driven ranked candidates with confidence scores
2. [classification icon] **ACMG Classification** — 5-tier variant classification with 23 criteria checklist
3. [therapy icon] **Therapeutic Options** — Gene therapy eligibility, ERT, orphan drugs, trial matching
4. [family tree icon] **Family Cascade Plan** — Sibling/parent testing, inheritance analysis, PGT-M
5. [network icon] **Integrated Assessment** — 4-agent coordinated output via /v1/diagnostic/integrated-assessment
6. [export icon] **Reports** — Markdown | JSON | PDF | FHIR R4 | **GA4GH Phenopacket v2** (unique to this agent)

**Streamlit UI (:8544) — 5 tabs as compact pills (10pt):**
Patient Intake | Diagnostic Dashboard | Variant Review | Therapeutic Options | Report Generator

**FastAPI (:8134) — Key endpoints as single line (10pt):**
/diagnose | /variants/interpret | /phenotype/match | /therapy/search | /v1/diagnostic/integrated-assessment

---

### ━━━ BOTTOM STRIP: DEMO + DIAGNOSTIC JOURNEY ━━━

**Left — Demo 4 walkthrough (RENDER AS VISUAL CARD AT 10pt MINIMUM):**

**Header (14pt bold, pink):** "Demo 4: One Gene, One Family — Aurora, 4yo, Bilateral Retinoblastoma"

Six steps as compact single lines (10pt):
1. HPO terms entered: Retinoblastoma + Microphthalmia
2. ACMG classifies RB1 c.958C>T as PATHOGENIC (12 points)
3. Oncology: Globe-sparing therapy, avoid radiation
4. Imaging: Trilateral screening — brain MRI q6mo until age 5
5. Trial matched: COG ARET0321
6. Sibling (age 2): 50% risk — cascade testing CRITICAL — PGT-M for future pregnancies

**Right — Rare Disease Statistics (RENDER AS LARGE METRIC BADGES, 14pt minimum per number):**

**Header (14pt bold):** "The Rare Disease Challenge"

Six large metric badges (number in 18pt bold, label in 10pt):
- **30M** Americans affected
- **7,000** Known rare diseases
- **5-7 yrs** Avg diagnostic odyssey
- **7+** Specialists before diagnosis
- **30%** Never receive diagnosis
- **50%** Are children

Single line (12pt bold, teal): "This agent compresses the diagnostic odyssey from years to minutes."

---

### ━━━ INFRASTRUCTURE BAR ━━━

Green bar (#76B900), white text, 12pt:
NVIDIA DGX Spark (GB10) | Milvus 2.4, 14 collections | BGE-small 384-dim | Claude Sonnet 4 | Ports 8544, 8134 | 206 tests

---

### ━━━ FOOTER ━━━

Single line, 10pt, gray:
Adam Jones | Apache 2.0 | hcls-ai-factory.org | v2.0 | March 2026 | ACMG: Richards et al. 2015 | HPO: Kohler et al. 2021 | GA4GH Phenopacket v2.0

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **88 rare diseases across 13 categories** — the broadest rare disease coverage of any open-source platform.
2. **23 of the 28 ACMG/AMP criteria** — 5-tier classification must be visually prominent with large color blocks.
3. **6 decision support engines** — deterministic clinical logic, not LLM guesses.
4. **9 diagnostic algorithms** — distinct clinical pathways from phenotype to treatment.
5. **Gene therapies tracked** — diagnosis connected to treatment options.
6. **GA4GH Phenopacket v2** — the only agent with this export format.
7. **Family cascade testing** — one diagnosis changes an entire family.
8. **Demo 4: Aurora's story** — the most emotionally powerful rare disease demo.
9. **Every character readable at 1920x1080** — no text below 10pt, key metrics at 14pt+.
