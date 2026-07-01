# Nano Banana Pro — Cardiology Intelligence Agent on NVIDIA DGX Spark

## IMPORTANT: Read this entire prompt before generating. This describes a single, dense, technical architecture infographic — NOT a slide deck. Every element described below appears on ONE canvas. This diagram represents the Cardiology Intelligence Agent — the most clinically validated agent in the HCLS AI Factory, featuring 6 published risk calculators (ASCVD, MAGGIC, EuroSCORE II, CHA₂DS₂-VASc, HAS-BLED, HEART Score), 11 clinical workflows, a GDMT optimization engine (ACC/AHA 2022), 13 Milvus collections, 45 cardiovascular conditions, 56 cardiac genes, 51 guideline references, and a 10-tab Streamlit UI. It connects to 5 peer agents and plays a critical role in Demo 3 ("Protecting the Survivor") for pediatric cardiotoxicity prevention. Landscape 16:9, reference architecture poster density.

---

## OVERALL LAYOUT AND STYLE

Landscape 16:9. White background (#FFFFFF). Dense reference architecture poster. Same visual series as the Imaging, CAR-T, Oncology, and Biomarker Agent infographics.

**Canvas:** Dense but organized. The diagram should feel like a cardiology department reference poster — comprehensive enough for a heart failure specialist, detailed enough for a clinical informaticist, clean enough for a hospital CTO.

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
- Red: #DC2626 — critical values, high-risk indicators, EuroSCORE high
- Amber: #F5A623 — moderate-risk, elevated markers
- Emerald: #059669 — low-risk, normal values, favorable outcomes
- Blue: #2196F3 — imaging, electrophysiology collections
- Purple: #7B2D8E — genomic, guideline collections
- Orange: #FF9800 — interventional, device collections
- Light Gray: #F5F5F5 — card backgrounds
- Gray: #666666 — secondary text
- Crimson: #B91C1C — cardiac emergency indicators

**Visual Elements:**
- Rounded rectangles (8px radius)
- Thin-line icons (16-24px, monochrome) — heart, ECG, stethoscope, pill, chart motifs
- Directional arrows: solid gray for data flow, dashed teal for cross-agent, bold red for critical alerts
- Color-coded risk stratification (green/amber/red)
- Metric badge pills throughout

---

## CANVAS STRUCTURE

### ━━━ BAND 1: TITLE BAR ━━━

**Left side:** Two stacked badges:
- "Cardiology Intelligence Agent" (green #76B900, white text)
- "6 Validated Risk Calculators" (red #DC2626, white text — red signifies clinical gravity)

**Center:**
- Line 1 (large, bold, navy): **Cardiology Intelligence Agent**
- Line 2 (gray): "Comprehensive Cardiovascular Intelligence on NVIDIA DGX Spark"
- Line 3 (teal): "Part of the Precision Intelligence Network — HCLS AI Factory"
- Line 4 (small, gray): "GB10 Superchip | 128 GB | 6 Calculators | 11 Workflows | 45 Conditions | 56 Genes | 63 Guidelines"

**Right side — Clinical capability badge:**
```
CLINICAL VALIDATION
━━━━━━━━━━━━━━━━━━━━━
6 Risk Calculators (published formulas)
GDMT Optimization (ACC/AHA 2022)
11 Clinical Workflows
Pediatric Cardio-Oncology
Deterministic — LLM never computes scores
Demo 3: "Protecting the Survivor"
Demo 5: "Last Line of Defense"
```

---

### ━━━ LEFT COLUMN: CLINICAL DATA DOMAINS ━━━

**Stacked input cards representing the 12 cardiology domains + shared genomic:**

1. **Literature** [journal icon]
   Cardiology research corpus
   Heart failure, CAD, arrhythmia
   [cardio_literature] [weight: 0.10]

2. **Clinical Trials** [trial icon]
   Cardiovascular trial data
   DAPA-HF, PARADIGM-HF, EMPEROR
   [cardio_trials] [weight: 0.08]

3. **Cardiac Imaging** [echo icon]
   Echocardiography, cardiac MRI
   CT angiography, nuclear
   [cardio_imaging] [weight: 0.10]

4. **Electrophysiology** [ECG icon]
   Arrhythmia, EP studies
   Device programming, ablation
   [cardio_electrophysiology] [weight: 0.08]

5. **Heart Failure** [heart-broken icon]
   HFrEF, HFpEF, HFmrEF
   GDMT protocols, staging
   [cardio_heart_failure] [weight: 0.10]

6. **Valvular Disease** [valve icon]
   Stenosis, regurgitation
   Surgical/TAVR criteria
   [cardio_valvular] [weight: 0.08]

7. **Prevention** [shield icon]
   Primary/secondary prevention
   Risk factor modification
   [cardio_prevention] [weight: 0.10]

8. **Interventional** [catheter icon]
   PCI, CABG, structural
   Revascularization criteria
   [cardio_interventional] [weight: 0.07]

9. **Cardio-Oncology** [heart+ribbon icon]
   Anthracycline cardiotoxicity
   Checkpoint inhibitor myocarditis
   [cardio_oncology] [weight: 0.06]

10. **Devices** [pacemaker icon]
    ICD, CRT, LVAD
    Programming, monitoring
    [cardio_devices] [weight: 0.04]

11. **Guidelines** [book icon]
    ACC/AHA, ESC, HRS
    51 guideline references
    [cardio_guidelines] [weight: 0.10]

12. **Hemodynamics** [pressure icon]
    Invasive hemodynamics
    Right heart cath, wedge
    [cardio_hemodynamics] [weight: 0.06]

13. **Genomic Evidence** [DNA icon]
    Shared collection
    56 cardiac genes
    [genomic_evidence] [weight: 0.03]

Arrows flow rightward: "BGE-small-en-v1.5 | 384-dim | 13 Milvus Collections"

---

### ━━━ CENTER-TOP: 6 RISK CALCULATORS (the headline feature) ━━━

**Section header (navy, green underline):** "6 Validated Risk Calculators — Deterministic, Published Formulas"

**Visual: 6 calculator cards arranged in a 3×2 grid. Each card has the calculator name, key inputs, output, and source citation. Color-coded risk output (green/amber/red scale).**

**Row 1:**

```
[heart+chart]                    [heart+broken]                [heart+surgery]
ASCVD 10-Year Risk               MAGGIC Mortality              EuroSCORE II
━━━━━━━━━━━━━━━━                 ━━━━━━━━━━━━━━                ━━━━━━━━━━━━
Inputs: Age, sex, race,          Inputs: Age, EF, NYHA,        Inputs: Age, sex, EF,
TC, HDL, SBP, DM, smoking,      SBP, BMI, creatinine,         creatinine, NYHA,
Rx status                        COPD, DM, HF duration         procedure type

Output: % 10-year risk           Output: 1yr/3yr mortality     Output: % operative
[Low <5%] [Moderate 5-20%]       [Low] [Moderate] [High]       mortality
[High >20%]                                                    [Low <2%] [High >5%]

Source: ACC/AHA 2013              Source: Pocock 2013           Source: Nashef 2012
Pooled Cohort Equations                                        STS/EACTS
```

**Row 2:**

```
[heart+lightning]                [heart+drop]                  [heart+flame]
CHA₂DS₂-VASc                    HAS-BLED                      HEART Score
━━━━━━━━━━━━━━━                  ━━━━━━━━━                     ━━━━━━━━━━━
Inputs: Age, sex, CHF,           Inputs: Hypertension,         Inputs: History, ECG,
hypertension, DM, stroke,        renal/liver disease,          age, risk factors,
vascular disease                 stroke, bleeding, age,        troponin
                                 drugs, alcohol, labile INR
Output: Annual stroke risk       Output: Annual bleed risk     Output: ACS risk
Score 0: 0.2% (males)           [Low 0-1] [Moderate 2]        [Low 0-3] [Moderate 4-6]
Score ≥2: anticoagulate          [High ≥3]                     [High 7-10]

Source: Lip 2010                 Source: Pisters 2010          Source: Six 2008
Chest 2010 (verified)
```

**Critical callout (emerald border):**
```
ALL CALCULATORS ARE DETERMINISTIC
Published formulas implemented in Python
CHA₂DS₂-VASc stroke rates: Lip GYH et al. Chest 2010
Score 0 includes sex-specific guidance (males: 0.2% annual risk)
LLM NEVER computes risk scores — only explains results
```

---

### ━━━ CENTER: GDMT OPTIMIZATION ENGINE + 11 WORKFLOWS ━━━

**Two side-by-side sections:**

#### Left half: GDMT Optimization Engine

```
GDMT Optimization Engine
ACC/AHA 2022 Heart Failure Guidelines
(Heidenreich et al., Circulation 2022)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Four-Pillar HFrEF GDMT:
[pill] Beta-Blocker (carvedilol, metoprolol, bisoprolol)
[pill] ARNI/ACEi/ARB (sacubitril-valsartan, enalapril)
[pill] MRA (spironolactone, eplerenone)
[pill] SGLT2i (dapagliflozin, empagliflozin)

Additional:
Hydralazine-ISDN | Ivabradine | Digoxin | Vericiguat

Device Therapy:
ICD eligibility | CRT criteria | CRT-D assessment

HFpEF:
SGLT2i | GLP-1 RA | Diuretics

[titration planning] [contraindication check] [interaction screen]
```

#### Right half: 11 Clinical Workflows

```
11 Clinical Workflows
━━━━━━━━━━━━━━━━━━━━

 1. CAD Assessment
 2. Heart Failure Management
 3. Valvular Disease Evaluation
 4. Arrhythmia Classification
 5. Cardiac MRI Interpretation
 6. Stress Test Analysis
 7. Preventive Risk Stratification
 8. Cardio-Oncology Surveillance
 9. Acute Decompensated HF
10. Post-MI Management
11. Myocarditis / Pericarditis

Each: preprocess → validate → execute → postprocess
Input validation on all 11 workflows
Template-method pattern (BaseCardioWorkflow)
```

---

### ━━━ CENTER-BOTTOM: PEDIATRIC CARDIO-ONCOLOGY + CROSS-AGENT ━━━

**Pediatric Cardio-Oncology box (emerald #059669 border, prominent — this is the Demo 3 feature):**

```
PEDIATRIC CARDIO-ONCOLOGY [NEW]
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Anthracycline Cardiotoxicity in Children:
• Acute: 5% incidence during therapy
• Late-onset: 5-25% (ages 10-30 post-treatment)
• Cumulative dose: >250 mg/m² moderate, >450 mg/m² high risk

Dexrazoxane Cardioprotection:
• 10:1 ratio to doxorubicin
• Recommended at >300 mg/m² cumulative dose (COG guidelines)

Pediatric Monitoring Protocol:
• Baseline echo: LVEF + fractional shortening before anthracyclines
• During therapy: echo every 3 months
• Post-treatment: annually for 5 years, then every 5 years lifetime

Pediatric Normal Values:
• LVEF >55% | FS >28% | Troponin I <0.04 ng/mL

pediatric_cardiotoxicity_assessment() function:
therapy_plan → cumulative dose check → dexrazoxane decision → monitoring schedule
```

**Cross-Agent Integration box (teal border):**
```
Cross-Agent Integration (/v1/cardio/integrated-assessment)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
→ Oncology Agent (:8527)       Planned chemotherapy, cumulative anthracycline dose
→ Clinical Trial Agent (:8538)  Trial-specific cardiac monitoring requirements
→ Biomarker Agent (:8529)      Troponin/BNP trends for cardiotoxicity detection
→ Imaging Agent (:8524)        Baseline echo coordination, cardiac MRI
→ Neurology Agent (:8528)      Shared toxicity monitoring (vincristine + cardiac)

+ pediatric_cardiotoxicity_assessment() — orchestrates all agents for pediatric chemo patients
```

---

### ━━━ RIGHT COLUMN: OUTPUTS & UI ━━━

**Output cards:**

1. **Risk Assessment Report** [chart icon]
   6-calculator output panel
   Color-coded risk stratification
   Guideline-linked recommendations

2. **GDMT Optimization Plan** [pill icon]
   Four-pillar therapy plan
   Titration schedule
   Contraindication alerts
   Device therapy eligibility

3. **Workflow Results** [workflow icon]
   11 workflow outputs
   Evidence-cited recommendations
   Clinical pathway navigation

4. **Integrated Assessment** [network icon]
   /v1/cardio/integrated-assessment
   5-agent coordinated output
   Pediatric cardiotoxicity report

5. **Reports** [export icon]
   Markdown | JSON | PDF | FHIR R4
   NVIDIA-branded PDF

**Streamlit UI (prominent box):**
```
Streamlit UI (:8536) — 10 Tabs
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
 1. Dashboard (health, metrics)
 2. Clinical Query (RAG Q&A)
 3. Risk Calculator (6 calculators)
 4. Heart Failure (GDMT optimizer)
 5. CAD Assessment (workflow)
 6. Arrhythmia (workflow)
 7. Imaging (echo/MRI review)
 8. Cardio-Oncology (surveillance)
 9. Evidence Explorer (collections)
10. Report Generator (multi-format)

Agent capabilities display:
[45 conditions] [32 drug classes]
[56 cardiac genes] [51 guidelines]
```

**FastAPI REST:**
```
Port :8126
/risk/ascvd | /risk/maggic | /risk/euroscore
/risk/cha2ds2-vasc | /risk/has-bled | /risk/heart-score
/gdmt/optimize | /workflow/* (11 endpoints)
/v1/cardio/integrated-assessment
```

---

### ━━━ BOTTOM STRIP: KNOWLEDGE STATS + DEMO PARTICIPATION ━━━

**Header (navy background, white text):** "Cardiology Intelligence Agent — By the Numbers"

**Left half — Knowledge statistics as a visual grid:**
```
45                  32                  56                  29                  63
Conditions          Drug Classes        Cardiac Genes       Biomarkers          Guideline Refs

Covering:
HCM • DCM • ARVC • Long QT • Brugada • AF • VT • STEMI • NSTEMI
Aortic Stenosis • Mitral Regurgitation • Heart Failure (HFrEF/HFpEF/HFmrEF)
Pericarditis • Myocarditis • Endocarditis • PAH • Cardiac Amyloidosis...
```
Each number as a large green metric badge.

**Right half — Demo participation:**
```
Demo 3 ✓ "Protecting the Survivor" — Marcus, 6yo, Neuroblastoma
  → Anthracycline cardiotoxicity prevention
  → Dexrazoxane at 300 mg/m² trigger
  → Baseline echo + monitoring schedule
  → pediatric_cardiotoxicity_assessment()

Demo 5 ✓ "Last Line of Defense" — Ethan, 12yo, CAR-T
  → Pre-lymphodepletion cardiac clearance
  → LVEF 58%, prior anthracycline 250 mg/m²
  → CRS cardiac monitoring plan
```

---

### ━━━ INFRASTRUCTURE BAR ━━━

**Green (#76B900) background, white text:**

```
DGX Spark (GB10)    Milvus 2.4          BGE-small-en-v1.5   Claude Sonnet 4     Docker Compose
$4,699              13 collections       384-dim vectors      Streaming RAG       Ports 8536, 8126
128 GB unified      IVF_FLAT/COSINE     Async embedding      + Cardio prompt     1,966 tests passing
```

---

### ━━━ FOOTER ━━━

```
Created by Adam Jones | Apache 2.0 Open Source | hcls-ai-factory.org | HCLS AI Factory v2.0 | March 2026
Risk Calculators: ACC/AHA 2013 (ASCVD) | Lip 2010 (CHA₂DS₂-VASc) | Pocock 2013 (MAGGIC) | Nashef 2012 (EuroSCORE II) | Pisters 2010 (HAS-BLED) | Six 2008 (HEART)
GDMT: Heidenreich et al., Circulation 2022;145:e895-e1032 | Pediatric: COG LTFU Guidelines
```

---

## WHAT THIS DIAGRAM MUST COMMUNICATE AT A GLANCE

1. **6 risk calculators — all from published clinical literature.** ASCVD, MAGGIC, EuroSCORE II, CHA₂DS₂-VASc, HAS-BLED, HEART Score. Each with named source citation. This is the most clinically validated agent.
2. **GDMT optimization engine** — ACC/AHA 2022 four-pillar HFrEF therapy with titration planning, contraindication checking, and device eligibility. Not a generic recommendation — a guideline-linked optimizer.
3. **11 clinical workflows** — from CAD assessment to myocarditis/pericarditis. Each with input validation and evidence-cited output.
4. **13 Milvus collections** — covering every cardiology subdomain (EP, HF, valvular, interventional, imaging, oncology, devices, hemodynamics, guidelines).
5. **45 conditions, 56 genes, 51 guidelines** — these numbers are specific and verified. Not rounded estimates.
6. **Pediatric cardio-oncology** — anthracycline cardiotoxicity prevention, dexrazoxane protocols, pediatric echo monitoring, COG LTFU guidelines. The pediatric_cardiotoxicity_assessment() function is unique to this agent.
7. **5 peer agent connections** — Oncology (chemo plans), Trial (monitoring requirements), Biomarker (troponin/BNP), Imaging (echo coordination), Neurology (shared toxicity).
8. **CHA₂DS₂-VASc verified against Lip 2010** — sex-specific guidance for score 0. This was audited and corrected during the certification process.
9. **10-tab Streamlit UI** — the most comprehensive UI of all specialty agents. Every cardiology subdomain has its own tab.
10. **Dense enough for a cardiology catheterization lab wall, precise enough for an ACC/AHA guideline committee, and human enough to protect a 6-year-old from heart failure 20 years after his cancer treatment.**
