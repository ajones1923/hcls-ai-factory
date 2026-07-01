# Nano Banana Pro — Pharmacogenomics Intelligence Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Pharmacogenomics Intelligence Agent — the HCLS AI Factory's right-drug-right-dose engine featuring 25 pharmacogenes (CYP2D6 with 150 star alleles, CYP2C19/38, CYP2C9/75, CYP3A4/34, and more), 12 HLA-drug hypersensitivity pairs, 12 drug categories, 5 metabolizer phenotype classifications, activity score algorithms (CYP2D6, DPYD), warfarin dosing calculator, phenoconversion modeling, CYP inhibitor/inducer databases, 15 Milvus collections, and a 10-tab Streamlit UI. It connects to 4 peer agents and is critical in Demo 3 ("Protecting the Survivor" — TPMT dosing for Marcus) and Demo 5 ("Last Line of Defense" — tocilizumab clearance for Ethan). Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual series as all other HCLS AI Factory agent infographics.

**Canvas:** Dense but organized. The diagram should feel like a clinical pharmacogenomics laboratory reference poster — comprehensive enough for a PharmD with PGx board certification, practical enough for a prescribing oncologist, clean enough for a hospital P&T committee.

**Typography:**
- Title: Large, bold, Inter/Helvetica, navy (#1B2333)
- Subtitle: Teal (#1AAFCC), smaller
- Section headers: Bold, navy, green left-border accent
- Sub-headers: Bold, teal (#1AAFCC)
- Body: Dark gray (#333333), 8-10pt, SHORT PHRASES only
- Metric badges: White text on colored pills
- Metabolizer phenotypes: Color-coded (red=poor, amber=intermediate, emerald=normal, blue=rapid, purple=ultra-rapid)
- All text MUST be legible at 1920x1080

**Color Palette:**
- NVIDIA Green: #76B900 — badges, infrastructure bar
- Navy: #1B2333 — title, headers
- Teal: #1AAFCC — sub-headers, data flow, agent connections
- Red: #DC2626 — Poor Metabolizer, contraindicated drug, HLA positive (do not prescribe)
- Amber: #F5A623 — Intermediate Metabolizer, dose adjustment needed
- Emerald: #059669 — Normal Metabolizer, standard dosing, HLA negative (safe)
- Blue: #2196F3 — Rapid Metabolizer, consider dose increase
- Purple: #7B2D8E — Ultra-Rapid Metabolizer, toxicity risk from prodrug conversion
- Orange: #FF9800 — HLA hypersensitivity, drug interaction warnings
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome) — pill, DNA, enzyme, warning triangle, dosing syringe motifs
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for contraindication alerts
- Color-coded 5-tier metabolizer phenotype spectrum
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Pharmacogenomics Intelligence Agent" (green #76B900, white text)
- "25 Pharmacogenes | 12 HLA Pairs" (purple #7B2D8E, white text)

**Center:**
- Line 1 (large, bold, navy): **Pharmacogenomics Intelligence Agent**
- Line 2 (gray): "Right Drug, Right Dose, Right Patient — Guided by Genotype"
- Line 3 (teal): "Part of the Precision Intelligence Engine — HCLS AI Factory"
- Line 4 (small, gray): "GB10 Superchip | 128 GB | 25 Genes | 308+ Star Alleles | 12 HLA Pairs | 12 Drug Categories | 15 Collections"

**Right side — Capability badge:**
```
PHARMACOGENOMIC DEPTH
━━━━━━━━━━━━━━━━━━━━━━━
25 Pharmacogenes Profiled
308+ Star Alleles Mapped
12 HLA-Drug Pairs
5 Metabolizer Phenotypes
Warfarin Dosing Algorithm
Phenoconversion Modeling
CPIC March 2025 Guidelines
Demo 3 + Demo 5
```

---

### ━━━ LEFT COLUMN: PHARMACOGENE DATABASE ━━━

**25 pharmacogenes as stacked cards, grouped by function:**

**Phase I Metabolism (CYP450 Enzymes):**

1. **CYP2D6** [enzyme icon] — "The most complex pharmacogene"
   150 star alleles defined
   Codeine, tamoxifen, atomoxetine, opioids
   Gene deletions, duplications, tandems, CYP2D7 hybrids
   [150 alleles] (purple badge) [VERY HIGH complexity]

2. **CYP2C19** [enzyme icon]
   38 star alleles
   Clopidogrel, PPIs, voriconazole, escitalopram
   [38 alleles] (blue badge)

3. **CYP2C9** [enzyme icon]
   75 star alleles
   Warfarin, phenytoin, NSAIDs
   [75 alleles] (blue badge)

4. **CYP3A4** [enzyme icon]
   34 star alleles
   ~30% of all drugs metabolized
   Tacrolimus, cyclosporine, statins
   [34 alleles] (teal badge)

5. **CYP3A5** [enzyme icon]
   11 star alleles
   Tacrolimus, vincristine (pediatric)
   *1/*1 expressers metabolize faster
   [11 alleles] (teal badge)

6. **CYP1A2 | CYP2B6 | CYP4F2 | CYP2C8** [enzyme icons]
   Additional CYP enzymes
   Clozapine, efavirenz, warfarin, paclitaxel

**Phase II Metabolism:**

7. **TPMT** [enzyme icon] — "Critical for pediatric ALL"
   Thiopurine methyltransferase
   6-MP, azathioprine | PM: 10% dose, IM: 50% dose
   [Pediatric critical] (red badge)

8. **DPYD** [enzyme icon]
   Dihydropyrimidine dehydrogenase
   5-FU, capecitabine | PM: contraindicated
   [CPIC Level A] (red badge)

9. **UGT1A1** [enzyme icon]
   Irinotecan, atazanavir
   *28/*28: 30% dose reduction
   [pediatric solid tumors] (amber badge)

10. **NUDT15** [enzyme icon]
    Thiopurine metabolism
    *3/*3: 10% dose (East Asian prevalence)
    [Pediatric critical] (red badge)

**Drug Transporters:**

11. **SLCO1B1** — Statin myopathy risk
12. **ABCG2** — Drug efflux transporter

**HLA Immunogenetics:**

13. **HLA-B / HLA-A** [immune icon]
    12 HLA-drug hypersensitivity pairs
    B*57:01 (abacavir), B*15:02 (carbamazepine SJS)
    B*58:01 (allopurinol SJS), A*31:01 (CBZ DRESS)
    [12 HLA pairs] (orange badge)

**Other:**
14. **VKORC1** — Warfarin sensitivity
15. **IFNL3** — Interferon response
16. **G6PD** — Hemolytic anemia risk
17-25. **Additional pharmacogenes**

Arrows rightward: "15 Collections | BGE-small-en-v1.5 | 384-dim"

---

### ━━━ CENTER-TOP: MILVUS COLLECTIONS — 15 ━━━

```
Milvus 2.4 | 15 Collections | IVF_FLAT / COSINE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
```

**Collection badges (3 rows):**

Row 1:
[pgx_gene_reference] [pgx_drug_guidelines] [pgx_drug_interactions] [pgx_hla_hypersensitivity] [pgx_phenoconversion]
  (purple)              (blue)                 (orange)                (red)                       (teal)

Row 2:
[pgx_dosing_algorithms] [pgx_clinical_evidence] [pgx_population_data] [pgx_clinical_trials] [pgx_fda_labels]
  (emerald)               (blue)                  (teal)                 (blue)                 (purple)

Row 3:
[pgx_drug_alternatives] [pgx_patient_profiles] [pgx_implementation] [pgx_education] [genomic_evidence]
  (emerald)               (amber)                (teal)               (blue)           (emerald — shared)

**15 collections — most of any agent in the platform (tied with Biomarker)**

---

### ━━━ CENTER: METABOLIZER PHENOTYPE ENGINE (the core) ━━━

**Section header (navy, green underline):** "Metabolizer Phenotype Classification — Genotype to Clinical Action"

**Visual: 5-tier metabolizer spectrum as a horizontal color bar (prominent, central):**

```
METABOLIZER PHENOTYPE SPECTRUM
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
[PURPLE]          [RED]              [AMBER]            [EMERALD]          [BLUE]
Ultra-Rapid       Poor               Intermediate       Normal             Rapid
Metabolizer       Metabolizer        Metabolizer        Metabolizer        Metabolizer

Excess active     No enzyme          Reduced enzyme     Standard enzyme    Increased
metabolite        activity           activity           activity           activity

Prodrug toxicity  Drug accumulation  Dose reduction     Standard dosing    Consider dose
(codeine→morphine (active drug       needed             recommended        increase
  overdose)       toxicity)

Example:          Example:           Example:           Example:           Example:
CYP2D6 UM +      TPMT *3A/*3A +     TPMT *1/*3A +      CYP2D6 *1/*1 +    CYP2C19 *1/*17 +
codeine → DEATH   6-MP → 10% dose    6-MP → 50% dose    tamoxifen → std    omeprazole → ↑dose
```

---

### ━━━ CENTER-MIDDLE: SPECIALIZED ENGINES ━━━

**Four processing engine cards:**

```
[dosing]                     [hla]                       [phenoconv]                [interactions]
Activity Score               HLA Screening               Phenoconversion            CYP Inhibitor/
Algorithms                   Engine                      Modeler                    Inducer Database
━━━━━━━━━━━━━━               ━━━━━━━━━━━━━               ━━━━━━━━━━━━━              ━━━━━━━━━━━━━━━

CYP2D6 activity scores       12 HLA-drug pairs           Drug-drug conversion       4 CYP enzymes
(CPIC 2019 system)           Screen BEFORE prescribing   of metabolizer status      Strong/moderate/weak
AS 0 → PM                    B*57:01 → NO abacavir       Paroxetine (CYP2D6         inhibitors mapped
AS 0.5-1.0 → IM             B*15:02 → NO carbamazepine    strong inhibitor) →
AS 1.5-2.0 → NM             B*58:01 → NO allopurinol      converts NM → PM         3 CYP enzymes
AS >2.0 → UM                A*31:01 → NO carbamazepine   phenotype for co-meds      strong/moderate
                                                                                    inducers mapped
DPYD activity scores         Population-specific          Real-time interaction
(CPIC 2018 system)           prevalence data              modeling
AS 0 → PM: NO 5-FU          [pre-test probability]       [medication list input]
AS 1.0 → IM: 50% dose

+ Warfarin Dosing Algorithm
  CYP2C9 + VKORC1 → predicted dose
  Clinical + genetic model (IWPC)
```

---

### ━━━ CENTER: 12 DRUG CATEGORIES ━━━

```
12 Drug Categories with PGx Relevance
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Category              Primary Genes              Clinical Impact
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Opioids               CYP2D6, CYP3A4            Codeine→morphine toxicity in UMs
Anticoagulants        CYP2C9, VKORC1, CYP2C19   Warfarin dosing, clopidogrel response
Antidepressants       CYP2D6, CYP2C19           SSRI/TCA metabolism and efficacy
Antipsychotics        CYP2D6, CYP1A2            Clozapine, aripiprazole metabolism
Statins               SLCO1B1                    Myopathy risk stratification
Chemotherapy          DPYD, UGT1A1, TPMT        5-FU toxicity, irinotecan, thiopurines
Antiepileptics        HLA-B, CYP2C9             Carbamazepine SJS/TEN, phenytoin
Antivirals            HLA-B, CYP2B6, UGT1A1     Abacavir hypersensitivity, efavirenz
Immunosuppressants    CYP3A5, TPMT              Tacrolimus dosing, azathioprine
Cardiovascular        CYP2D6, CYP2C9            Beta-blockers, antiarrhythmics
PPIs                  CYP2C19                    Omeprazole/pantoprazole efficacy
Gout                  HLA-B, G6PD               Allopurinol SJS/TEN risk
```

---

### ━━━ CENTER-BOTTOM: PEDIATRIC PGx + CROSS-AGENT ━━━

**Pediatric PGx box (emerald border):**
```
PEDIATRIC PHARMACOGENOMICS
━━━━━━━━━━━━━━━━━━━━━━━━━━

TPMT for 6-Mercaptopurine in ALL (Demo 3 — Marcus):
• *1/*3A (IM): Start at 50% dose → titrate by TGN levels
• *3A/*3A (PM): Start at 10% dose → weekly CBC monitoring
• NUDT15 *3/*3: 10% dose (higher prevalence East Asian)

CYP3A5 for Vincristine:
• *1/*1 expressers: rapid clearance → monitor efficacy
• *3/*3 non-expressers: standard dosing

UGT1A1 for Irinotecan (Pediatric Solid Tumors):
• *28/*28: 30% dose reduction → monitor severe diarrhea

MTHFR for Methotrexate:
• C677T: reduced folate metabolism → increased MTX toxicity risk

Asparaginase: No established PGx markers
• Anti-asparaginase antibodies affect efficacy → switch to Erwinia
```

**Cross-Agent box (teal border):**
```
Cross-Agent (/v1/pgx/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Oncology Agent (:8527)      Planned therapy for drug-gene interaction context
→ Cardiology Agent (:8126)    Cardiac drug interactions (QT prolongation, warfarin)
→ Neurology Agent (:8528)     Neurotoxic drug interactions (vincristine, methotrexate)
→ Clinical Trial Agent (:8538) PGx-guided trial matching

4 peer agents | Dosing integrated before first infusion
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

1. **Drug Check Report** [pill+check icon]
   Single drug PGx assessment
   Metabolizer phenotype + dosing recommendation
   Alternative drug suggestions

2. **Medication Review** [list icon]
   Full medication list analysis
   Drug-drug-gene interactions
   Phenoconversion alerts

3. **Warfarin Dosing** [dosing icon]
   CYP2C9 + VKORC1 algorithm
   Predicted weekly dose (mg)
   IWPC clinical + genetic model

4. **HLA Screening** [shield icon]
   Pre-prescription safety screen
   12 HLA-drug pairs checked
   GO / DO NOT PRESCRIBE result

5. **Integrated Assessment** [network icon]
   /v1/pgx/integrated-assessment
   4-agent coordinated output

6. **Reports** [export icon]
   Markdown | JSON | PDF | FHIR R4

**Streamlit UI (prominent — 10 tabs):**
```
Streamlit UI (:8507) — 10 Tabs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 1. PGx Dashboard (gene overview)
 2. Drug Check (single drug query)
 3. Medication Review (full list)
 4. Warfarin Dosing (CYP2C9+VKORC1)
 5. Chemotherapy Safety (DPYD/TPMT/UGT1A1)
 6. HLA Screening (12 allele-drug pairs)
 7. PGx Report Generator (multi-format)
 8. Evidence Explorer (RAG Q&A)
 9. Phenoconversion Modeler (DDI impact)
10. Population Analytics (allele frequencies)

CYP inhibitor/inducer reference tables
```

**FastAPI (:8107):**
```
/drug-check | /medication-review | /dosing/warfarin
/hla-screen | /phenoconversion | /genes | /drugs
/v1/pgx/integrated-assessment
```

---

### ━━━ BOTTOM STRIP: DEMO + CLINICAL DECISION TABLE ━━━

**Left — Demo participation:**
```
Demo 3 ✓ "Protecting the Survivor" — Marcus, 6yo, Neuroblastoma
  → TPMT *1/*3A: Intermediate Metabolizer → 50% 6-MP dose
  → CYP3A5: Vincristine metabolism screening
  → CBR3/RARG: Anthracycline cardiotoxicity genotype
  → Dosing adjustments BEFORE first infusion

Demo 5 ✓ "Last Line of Defense" — Ethan, 12yo, CAR-T
  → Tocilizumab: No CYP interactions (monoclonal antibody)
  → Corticosteroids: CYP3A4 metabolism cleared
  → All supportive medications PGx-screened
```

**Right — Clinical Decision Framework:**
```
The PGx Clinical Decision:
━━━━━━━━━━━━━━━━━━━━━━━━━

Genotype → Star Allele → Activity Score → Metabolizer Phenotype → Clinical Action

CYP2D6 *1/*4   → AS 1.0  → IM  → Reduce codeine 50%, consider alternative
CYP2C19 *2/*2  → AS 0    → PM  → Switch clopidogrel → ticagrelor
TPMT *3A/*3A   → AS 0    → PM  → 6-MP at 10% dose → weekly CBC
DPYD *2A/wt    → AS 1.0  → IM  → 5-FU at 50% dose → monitor DPD
HLA-B*57:01+   → N/A     → N/A → DO NOT prescribe abacavir → tenofovir

The cost of NOT testing: adverse drug reactions cause
~100,000 deaths/year in the US (Lazarou 1998, JAMA)
This agent prevents them.
```

---

### ━━━ INFRASTRUCTURE BAR ━━━

Green (#76B900), white text:
```
DGX Spark (GB10) $4,699 | Milvus 2.4, 15 collections | BGE-small 384-dim | Claude Sonnet 4 | Ports 8507, 8107 | 200+ tests
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 | hcls-ai-factory.org | HCLS AI Factory v2.0 | March 2026
CPIC: March 2025 Guidelines | PharmGKB: Variant annotations | FDA: Pharmacogenomic biomarker table
Star alleles: PharmVar (CYP2D6 150, CYP2C19 38, CYP2C9 75) | HLA: IPD-IMGT/HLA Database
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **25 pharmacogenes with 308+ star alleles** — from CYP2D6 (150 alleles, the most complex) to NUDT15. The gene database must be visually dominant.
2. **5-tier metabolizer phenotype spectrum** — Ultra-Rapid through Poor, each with clinical consequence and example. The color-coded horizontal bar must be the visual centerpiece.
3. **12 HLA-drug hypersensitivity pairs** — screen BEFORE prescribing. HLA-B*57:01/abacavir, B*15:02/carbamazepine SJS, B*58:01/allopurinol SJS. These prevent fatal reactions.
4. **12 drug categories** — opioids through gout, each with primary genes and clinical impact.
5. **15 Milvus collections** — tied for most in the platform. From gene reference to FDA labels to population data.
6. **Activity score algorithms** — CYP2D6 and DPYD CPIC systems with clear AS→phenotype→action mapping.
7. **Warfarin dosing calculator** — CYP2C9 + VKORC1 clinical + genetic model (IWPC).
8. **Phenoconversion modeling** — unique capability: drug-drug interactions that change metabolizer phenotype in real-time.
9. **Pediatric PGx** — TPMT/NUDT15 for 6-MP in ALL, CYP3A5 for vincristine, UGT1A1 for irinotecan, MTHFR for MTX. Dosing before the first infusion.
10. **Dense enough for a pharmacogenomics laboratory, practical enough for a prescribing oncologist, and life-saving enough to prevent the 100,000 annual US deaths from adverse drug reactions.**
