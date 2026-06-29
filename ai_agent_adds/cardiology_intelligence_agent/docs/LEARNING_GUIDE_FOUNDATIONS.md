# Cardiology Intelligence Agent -- Learning Guide: Foundations

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

A primer on the cardiology fundamentals underlying the Cardiology Intelligence Agent. Designed for engineers, data scientists, and non-cardiologists working with the system.

---

## Table of Contents

1. [Cardiac Anatomy and Physiology Basics](#1-cardiac-anatomy-and-physiology-basics)
2. [Heart Failure Pathophysiology](#2-heart-failure-pathophysiology)
3. [Coronary Artery Disease](#3-coronary-artery-disease)
4. [ECG Interpretation Basics](#4-ecg-interpretation-basics)
5. [Echocardiography Basics](#5-echocardiography-basics)
6. [Risk Score Fundamentals](#6-risk-score-fundamentals)
7. [GDMT Four-Pillar Framework](#7-gdmt-four-pillar-framework)
8. [RAG Architecture Introduction](#8-rag-architecture-introduction)

---

## 1. Cardiac Anatomy and Physiology Basics

### 1.1 The Four Chambers

The heart consists of four chambers that work as two synchronized pumps:

```
                    Superior Vena Cava
                          |
              +-----------+-----------+
              |     RIGHT ATRIUM      |
              |   (receives venous    |
              |    blood from body)   |
              +-----------+-----------+
                          |
                    Tricuspid Valve
                          |
              +-----------+-----------+
              |   RIGHT VENTRICLE     |
              |   (pumps blood to     |
              |    lungs for O2)      |
              +-----------+-----------+
                          |
                    Pulmonic Valve
                          |
                    Pulmonary Artery --> LUNGS --> Pulmonary Veins
                                                         |
                                              +----------+----------+
                                              |     LEFT ATRIUM     |
                                              |   (receives         |
                                              |    oxygenated blood)|
                                              +----------+----------+
                                                         |
                                                   Mitral Valve
                                                         |
                                              +----------+----------+
                                              |    LEFT VENTRICLE   |
                                              |   (pumps blood to   |
                                              |    entire body)     |
                                              +----------+----------+
                                                         |
                                                   Aortic Valve
                                                         |
                                                       Aorta --> Body
```

**Key measurements in the system:**
- **LVIDd (LV internal dimension, diastole):** 42-58 mm (measures LV size)
- **LA volume index:** <34 mL/m2 (left atrial size, elevated in diastolic dysfunction)
- **RV diameter:** <42 mm at base (right ventricular size)

### 1.2 The Cardiac Cycle

The cardiac cycle has two main phases:

| Phase | Description | Key Events |
|-------|-------------|-----------|
| **Systole** | Contraction | Ventricles contract, eject blood through aortic and pulmonic valves |
| **Diastole** | Relaxation | Ventricles relax and fill from atria through mitral and tricuspid valves |

**Ejection Fraction (EF):** The percentage of blood pumped out of the left ventricle with each heartbeat.

```
EF = (End-diastolic volume - End-systolic volume) / End-diastolic volume x 100

Normal LVEF: 55-70%
```

EF is the single most important measurement in the system, as it determines:
- Heart failure classification (HFrEF vs HFmrEF vs HFpEF)
- GDMT eligibility
- Device therapy eligibility (ICD, CRT)
- Prognosis

### 1.3 The Coronary Arteries

Three major coronary arteries supply blood to the heart muscle:

| Artery | Abbreviation | Territory |
|--------|-------------|-----------|
| Left Anterior Descending | LAD | Anterior wall, septum, apex |
| Left Circumflex | LCx | Lateral wall, posterior wall (some) |
| Right Coronary Artery | RCA | Inferior wall, RV, SA/AV nodes (usually) |
| Left Main | LM | Bifurcates into LAD and LCx |

**Why this matters for the system:** The CAD Assessment workflow evaluates stenosis severity using CAD-RADS scoring, which grades stenosis percentage in each vessel.

### 1.4 The Conduction System

Electrical signals control the heart rhythm:

```
SA Node (60-100 bpm, "natural pacemaker")
    |
    v
Atrial Depolarization (P wave on ECG)
    |
    v
AV Node (delay ~120-200ms = PR interval)
    |
    v
Bundle of His
    |
    +---> Left Bundle Branch (LBBB if blocked)
    |         |
    |         +---> Left Anterior Fascicle
    |         +---> Left Posterior Fascicle
    |
    +---> Right Bundle Branch (RBBB if blocked)
    |
    v
Purkinje Fibers
    |
    v
Ventricular Depolarization (QRS complex on ECG, <120ms normally)
```

**Why this matters for the system:** The Arrhythmia workflow interprets ECG findings and the GDMT optimizer evaluates CRT eligibility (requires LBBB + QRS >=150ms).

### 1.5 Key Cardiac Biomarkers

| Biomarker | What It Measures | Why It Matters |
|-----------|-----------------|---------------|
| **Troponin (hs-cTnI/T)** | Myocardial cell death/injury | Heart attack diagnosis |
| **BNP / NT-proBNP** | Ventricular wall stress | Heart failure diagnosis and monitoring |
| **LDL cholesterol** | "Bad cholesterol" | ASCVD risk, statin eligibility |
| **Potassium (K+)** | Electrolyte level | GDMT safety (MRA contraindication if >5.0) |
| **Creatinine / eGFR** | Kidney function | GDMT dosing, DOAC adjustment |

---

## 2. Heart Failure Pathophysiology

### 2.1 What Is Heart Failure?

Heart failure is a clinical syndrome where the heart cannot pump blood efficiently enough to meet the body's needs. It is classified by:

**By Ejection Fraction:**

| Category | LVEF | Agent Enum | Key Feature |
|----------|------|-----------|-------------|
| **HFrEF** | <=40% | `EjectionFractionCategory.HFrEF` | Systolic dysfunction; full GDMT indicated |
| **HFmrEF** | 41-49% | `EjectionFractionCategory.HFmrEF` | "Mildly reduced"; SGLT2i + consider GDMT |
| **HFpEF** | >=50% | `EjectionFractionCategory.HFpEF` | Diastolic dysfunction; SGLT2i + comorbidity management |
| **HFimpEF** | Was <=40%, now >40% | `EjectionFractionCategory.HFimpEF` | Improved; continue all GDMT indefinitely |

**By NYHA Functional Class:**

| Class | Agent Enum | Symptoms |
|-------|-----------|----------|
| I | `HeartFailureClass.NYHA_I` | No limitation of physical activity |
| II | `HeartFailureClass.NYHA_II` | Slight limitation; comfortable at rest |
| III | `HeartFailureClass.NYHA_III` | Marked limitation; comfortable only at rest |
| IV | `HeartFailureClass.NYHA_IV` | Unable to carry on any physical activity without discomfort |

**By ACC/AHA Stage:**

| Stage | Agent Enum | Description |
|-------|-----------|-------------|
| A | `HeartFailureStage.STAGE_A` | At risk (HTN, DM, CAD) but no structural disease |
| B | `HeartFailureStage.STAGE_B` | Structural disease (reduced EF, LVH, valve disease) but no symptoms |
| C | `HeartFailureStage.STAGE_C` | Structural disease WITH symptoms |
| D | `HeartFailureStage.STAGE_D` | Advanced HF refractory to standard therapy |

### 2.2 The Neurohormonal Model

Heart failure activates harmful neurohormonal pathways that GDMT medications target:

```
Heart failure (low cardiac output)
    |
    +---> Sympathetic activation --> Elevated heart rate, vasoconstriction
    |     TARGET: Beta-blockers (reduce HR, remodeling)
    |
    +---> RAAS activation --> Angiotensin II, aldosterone
    |     TARGET: ARNI/ACEi/ARB (block angiotensin)
    |              MRA (block aldosterone)
    |
    +---> Sodium/water retention --> Volume overload, congestion
    |     TARGET: SGLT2i (natriuresis), Loop diuretics (volume removal)
    |
    +---> Cardiac remodeling --> Progressive dilation, fibrosis
          TARGET: All 4 pillars reduce remodeling
```

### 2.3 Key Heart Failure Biomarkers

| Biomarker | Rule-Out Threshold | Elevated Threshold | Clinical Use |
|-----------|-------------------|-------------------|-------------|
| BNP | <100 pg/mL | >400 pg/mL | HF diagnosis |
| NT-proBNP | <300 pg/mL | Age-adjusted: >450/<900/<1800 | HF diagnosis, GDMT monitoring |
| hs-Troponin | 99th percentile URL | Rising pattern | Myocardial injury detection |

---

## 3. Coronary Artery Disease

### 3.1 Atherosclerosis

Coronary artery disease is caused by atherosclerotic plaque buildup in the coronary arteries. The process:

1. **Fatty streak**: Lipid deposition in arterial wall (begins in youth)
2. **Fibroatheroma**: Plaque with lipid core and fibrous cap
3. **Progressive stenosis**: Gradual narrowing of the vessel lumen
4. **Plaque rupture**: Unstable plaque ruptures, triggering thrombosis -> heart attack

### 3.2 CAD-RADS Classification

The CAD Assessment workflow uses CAD-RADS scoring from coronary CT angiography:

| CAD-RADS | Agent Enum | Stenosis | Management |
|----------|-----------|----------|-----------|
| 0 | `CADRADSScore.CAD_RADS_0` | 0% (no plaque) | No further workup |
| 1 | `CADRADSScore.CAD_RADS_1` | 1-24% (minimal) | Preventive therapy |
| 2 | `CADRADSScore.CAD_RADS_2` | 25-49% (mild) | Preventive therapy |
| 3 | `CADRADSScore.CAD_RADS_3` | 50-69% (moderate) | Consider functional testing |
| 4A | `CADRADSScore.CAD_RADS_4A` | 70-99% (severe, 1-2 vessels) | Consider invasive angiography |
| 4B | `CADRADSScore.CAD_RADS_4B` | 70-99% (severe, 3 vessels or LM >=50%) | Invasive angiography recommended |
| 5 | `CADRADSScore.CAD_RADS_5` | 100% (total occlusion) | Invasive angiography |

### 3.3 Functional Significance

Not all anatomic stenosis causes ischemia. Functional testing determines significance:

| Test | Positive Result | Significance |
|------|----------------|-------------|
| **FFR** (catheterization) | <=0.80 | Hemodynamically significant |
| **iFR** (catheterization) | <=0.89 | Hemodynamically significant |
| **Stress echo** | New wall motion abnormality | Inducible ischemia |
| **SPECT MPI** | Reversible perfusion defect | Inducible ischemia |
| **Cardiac PET** | Reduced myocardial blood flow | Inducible ischemia |

### 3.4 Coronary Artery Calcium Scoring

Coronary artery calcium (CAC) quantifies calcified plaque using non-contrast CT:

| Agatston Score | Risk Category | Statin Decision |
|---------------|--------------|----------------|
| 0 | Very low | May defer statin if intermediate risk |
| 1-99 | Low | Favors statin if intermediate risk |
| 100-399 | Moderate | Statin recommended |
| >=400 | High | Statin recommended; consider further evaluation |

**Why this matters for the system:** The Prevention workflow uses CAC scoring for risk reclassification in patients with intermediate 10-year ASCVD risk (7.5-20%).

---

## 4. ECG Interpretation Basics

### 4.1 The 12-Lead ECG

The ECG records the heart's electrical activity from 12 perspectives (leads):

| Lead Group | Leads | View |
|-----------|-------|------|
| Limb leads | I, II, III | Frontal plane |
| Augmented leads | aVR, aVL, aVF | Frontal plane |
| Precordial leads | V1-V6 | Horizontal plane |

### 4.2 Key Intervals and Normal Values

| Interval | Normal | Clinical Significance |
|----------|--------|---------------------|
| **Heart Rate** | 60-100 bpm | Bradycardia <60, Tachycardia >100 |
| **PR Interval** | 120-200 ms | >200 = first-degree AV block; <120 = pre-excitation (WPW) |
| **QRS Duration** | <120 ms | 120-149 = incomplete bundle branch block; >=150 = complete BBB |
| **QTc Interval** | <450 ms (male), <460 ms (female) | >500 ms = high risk for torsades de pointes |
| **Axis** | -30 to +90 degrees | LAD (<-30), RAD (>+90) |

### 4.3 Common ECG Findings in the System

| Finding | ECG Features | Clinical Implication |
|---------|-------------|---------------------|
| **Atrial Fibrillation** | Irregular rhythm, no P waves, fibrillatory baseline | CHA2DS2-VASc calculation, anticoagulation |
| **ST Elevation** | ST elevation >=1mm in 2 contiguous leads | Acute MI (STEMI) -- emergent cath lab |
| **Left Bundle Branch Block** | QRS >=120ms, broad R in I/aVL/V5-V6 | CRT eligibility if HFrEF |
| **Long QT** | QTc >480ms | Risk of torsades; check for LQTS genes |
| **Brugada Pattern** | Coved ST elevation in V1-V3 | Risk of sudden death; consider ICD |
| **LVH** | Voltage criteria + repolarization abnormalities | Screen for HCM if unexplained |

### 4.4 How the System Uses ECG Data

The `ECGInterpretation` Pydantic model captures structured ECG data:
- **rhythm**: Free text rhythm interpretation
- **rate**: Ventricular rate in bpm
- **intervals**: Dictionary with PR, QRS, QTc values in milliseconds
- **axis**: Electrical axis description
- **findings**: List of identified abnormalities
- **urgency**: SeverityLevel classification

---

## 5. Echocardiography Basics

### 5.1 What Is Echocardiography?

Echocardiography (echo) uses ultrasound to image the heart in real-time. It is the most commonly ordered cardiac imaging test and provides:

- **Chamber sizes**: Is the heart dilated?
- **Wall thickness**: Is there hypertrophy?
- **Systolic function**: How well does the heart contract? (LVEF)
- **Diastolic function**: How well does the heart relax and fill?
- **Valve function**: Are valves stenotic (narrowed) or regurgitant (leaking)?
- **Wall motion**: Are specific segments abnormal? (ischemia, infarction)

### 5.2 Key Echo Measurements

| Measurement | Normal Range | Clinical Significance |
|-------------|-------------|---------------------|
| **LVEF** | 55-70% | <40% = HFrEF; 41-49% = HFmrEF; >=50% = normal/HFpEF |
| **LVIDd** | 42-58 mm | >58 mm = LV dilation (DCM, chronic volume overload) |
| **IVS/PW thickness** | 6-11 mm | >=15 mm = consider HCM |
| **LA volume index** | <34 mL/m2 | Elevated in diastolic dysfunction, mitral disease, AF |
| **E/e'** | <14 | >=14 = elevated filling pressures (diastolic dysfunction) |
| **TAPSE** | >17 mm | <17 mm = RV systolic dysfunction |
| **GLS** | <- 18% | Less negative = subclinical dysfunction; used in cardio-oncology |
| **TR velocity** | <2.8 m/s | >2.8 = elevated PASP |

### 5.3 Valve Assessment

Each valve is assessed for stenosis (narrowing) and regurgitation (leaking):

**Aortic Stenosis Severity:**

| Parameter | Mild | Moderate | Severe |
|-----------|------|----------|--------|
| Vmax | 2.0-2.9 m/s | 3.0-3.9 m/s | >=4.0 m/s |
| Mean gradient | <20 mmHg | 20-39 mmHg | >=40 mmHg |
| AVA | >1.5 cm2 | 1.0-1.5 cm2 | <1.0 cm2 |

These correspond to the `ValveSeverity` enum: MILD, MODERATE, SEVERE, CRITICAL.

### 5.4 Global Longitudinal Strain (GLS)

GLS measures myocardial deformation using speckle tracking:

- **Normal:** < -18% (more negative = better function)
- **Subclinical dysfunction:** -14% to -18%
- **Reduced:** > -14% (less negative)

**Why this matters for the system:** The Cardio-Oncology workflow uses GLS as an early marker of cardiotoxicity. A >15% relative decline from baseline triggers a CTRCD alert, even if LVEF remains normal.

---

## 6. Risk Score Fundamentals

### 6.1 What Are Clinical Risk Scores?

Clinical risk scores are validated mathematical models that predict patient outcomes. They convert multiple clinical variables into a single number that guides treatment decisions.

### 6.2 ASCVD Pooled Cohort Equations

**Purpose:** Estimate 10-year risk of a first atherosclerotic cardiovascular event (heart attack or stroke).

**Inputs:** Age, sex, race, total cholesterol, HDL, systolic BP, BP treatment, diabetes, smoking.

**How it works:** Uses log-transformed variables with sex/race-specific beta coefficients in a survival function:

```
Risk = 1 - S0^exp(individual_sum - mean_coefficient)
```

Where S0 is the 10-year baseline survival for the patient's sex/race cohort.

**Clinical decision thresholds:**

| Risk Level | 10-Year Risk | Action |
|-----------|-------------|--------|
| Low | <5% | Lifestyle modification |
| Borderline | 5-7.5% | Consider risk enhancers |
| Intermediate | 7.5-20% | Moderate-intensity statin; consider CAC |
| High | >=20% | High-intensity statin |

### 6.3 CHA2DS2-VASc

**Purpose:** Estimate annual stroke risk in atrial fibrillation patients.

**Scoring:**

| Factor | Points |
|--------|--------|
| **C** - Congestive heart failure | 1 |
| **H** - Hypertension | 1 |
| **A2** - Age >=75 | 2 |
| **D** - Diabetes | 1 |
| **S2** - Stroke/TIA history | 2 |
| **V** - Vascular disease | 1 |
| **A** - Age 65-74 | 1 |
| **Sc** - Sex category (female) | 1 |

**Decision rule:**
- Score 0 (males) or 1 (females): No anticoagulation needed
- Score 1 (males) or 2 (females): Consider anticoagulation
- Score >=2 (males) or >=3 (females): Anticoagulation recommended

### 6.4 HAS-BLED

**Purpose:** Assess bleeding risk on anticoagulation (used alongside CHA2DS2-VASc).

**Scoring:**

| Factor | Points |
|--------|--------|
| **H** - Hypertension (uncontrolled, SBP >160) | 1 |
| **A** - Abnormal renal or liver function | 1-2 |
| **S** - Stroke history | 1 |
| **B** - Bleeding history or predisposition | 1 |
| **L** - Labile INR (TTR <60%) | 1 |
| **E** - Elderly (age >65) | 1 |
| **D** - Drugs (antiplatelets/NSAIDs) or alcohol | 1-2 |

**Interpretation:** Score >=3 = high bleeding risk. Does NOT mean "don't anticoagulate" -- means closer monitoring required.

### 6.5 HEART Score

**Purpose:** Risk-stratify chest pain patients in the emergency department for major adverse cardiac events (MACE).

**Scoring (0-2 each):**

| Factor | 0 | 1 | 2 |
|--------|---|---|---|
| **H** - History | Slightly suspicious | Moderately suspicious | Highly suspicious |
| **E** - ECG | Normal | Non-specific changes | Significant ST deviation |
| **A** - Age | <45 | 45-64 | >=65 |
| **R** - Risk factors | None | 1-2 | >=3 or h/o ASCVD |
| **T** - Troponin | Normal | 1-3x URL | >3x URL |

**Risk categories:** Low (0-3): 1.7% MACE, Moderate (4-6): 16.6% MACE, High (7-10): 50.1% MACE

### 6.6 MAGGIC

**Purpose:** Predict 1-year and 3-year mortality in heart failure patients.

**Key variables:** Age, sex, LVEF, NYHA class, SBP, BMI, creatinine, diabetes, beta-blocker use, ACEi/ARB use.

**Output:** Integer score (0-50) mapped to mortality percentage via published lookup table.

### 6.7 EuroSCORE II

**Purpose:** Predict operative mortality for cardiac surgery.

**Key variables:** 28 factors across patient, cardiac, and operation categories.

**Output:** Predicted operative mortality percentage. Used in the system's Valvular Disease workflow to assess surgical risk for TAVR vs SAVR decisions.

---

## 7. GDMT Four-Pillar Framework

### 7.1 What Is GDMT?

Guideline-Directed Medical Therapy (GDMT) refers to the four medication classes that have been proven in randomized controlled trials to reduce mortality and hospitalization in heart failure with reduced ejection fraction (HFrEF, LVEF <=40%).

### 7.2 The Four Pillars

| Pillar | Drug Class | Agent Enum | Mechanism | Key Evidence |
|--------|-----------|-----------|-----------|-------------|
| 1 | Beta-blocker | `GDMTPillar.BETA_BLOCKER` | Blocks sympathetic overstimulation | MERIT-HF, COPERNICUS, CIBIS-II |
| 2 | ARNI (or ACEi/ARB) | `GDMTPillar.ARNI_ACEI_ARB` | Blocks RAAS + enhances natriuretic peptides | PARADIGM-HF |
| 3 | MRA | `GDMTPillar.MRA` | Blocks aldosterone; anti-fibrotic | RALES, EMPHASIS-HF |
| 4 | SGLT2 inhibitor | `GDMTPillar.SGLT2I` | Natriuresis, osmotic diuresis, cardiac remodeling | DAPA-HF, EMPEROR-Reduced |

### 7.3 GDMT Status Tracking

The `GDMTStatus` enum tracks each pillar's medication status:

| Status | Agent Enum | Description |
|--------|-----------|-------------|
| Not started | `GDMTStatus.NOT_STARTED` | Medication not yet initiated |
| Initiated | `GDMTStatus.INITIATED` | Started at low dose |
| Uptitrating | `GDMTStatus.UPTITRATING` | Dose being increased toward target |
| At target | `GDMTStatus.AT_TARGET` | At guideline-recommended target dose |
| Contraindicated | `GDMTStatus.CONTRAINDICATED` | Cannot use due to clinical contraindication |
| Intolerant | `GDMTStatus.INTOLERANT` | Cannot tolerate (side effects) |

### 7.4 Titration Example

```
Patient: 60-year-old male, LVEF 30%, NYHA II

Step 1 (Day 0):
  Start carvedilol 3.125mg BID + sacubitril/valsartan 24/26mg BID + dapagliflozin 10mg
  (SGLT2i needs no titration)

Step 2 (Week 2):
  If HR >55 and SBP >100: uptitrate carvedilol to 6.25mg BID
  Check K+, creatinine

Step 3 (Week 4):
  Uptitrate sacubitril/valsartan to 49/51mg BID
  If K+ <5.0 and eGFR >30: add spironolactone 12.5mg
  Check K+ at 1 week after MRA start

Step 4 (Week 6):
  Uptitrate carvedilol to 12.5mg BID

Step 5 (Week 8):
  Uptitrate sacubitril/valsartan to 97/103mg BID (target)
  Uptitrate spironolactone to 25mg (if K+ still <5.0)

Step 6 (Week 10-12):
  Uptitrate carvedilol to 25mg BID (target)

At target: All 4 pillars at guideline-recommended doses.
```

### 7.5 Why GDMT Matters

The mortality reduction from each pillar is roughly additive:

| Therapy | Relative Risk Reduction (Mortality) |
|---------|-----------------------------------|
| Beta-blocker alone | ~34% |
| ACEi alone | ~23% |
| MRA added to ACEi + BB | ~30% additional |
| ARNI vs ACEi | ~20% additional |
| SGLT2i added to background GDMT | ~18% additional |
| **Combined 4-pillar vs placebo** | **~60-70% estimated** |

Despite this evidence, real-world data shows that fewer than 25% of eligible patients receive target-dose GDMT for all four pillars. The GDMT optimizer addresses this gap.

---

## 8. RAG Architecture Introduction

### 8.1 What Is RAG?

Retrieval-Augmented Generation (RAG) is an AI architecture that combines:

1. **Retrieval**: Search a knowledge base for relevant information
2. **Augmentation**: Provide retrieved information as context to an LLM
3. **Generation**: LLM generates a response grounded in the retrieved evidence

```
          User Query
              |
         [RETRIEVAL]
              |
    Search vector database
    for similar content
              |
    Top-K most relevant
    documents retrieved
              |
        [AUGMENTATION]
              |
    "Here is the question + here is
     the relevant evidence from our
     knowledge base"
              |
        [GENERATION]
              |
    LLM synthesizes an answer
    grounded in the evidence
    with citations
              |
         Response
```

### 8.2 Why RAG for Cardiology?

RAG solves critical problems for clinical AI:

| Problem | How RAG Solves It |
|---------|------------------|
| **Hallucination** | LLM only uses retrieved evidence, not fabricated knowledge |
| **Stale knowledge** | Vector DB is updated regularly; LLM doesn't need retraining |
| **Citation provenance** | Every claim traces to a specific source document |
| **Domain specificity** | Vector DB contains only curated cardiovascular content |
| **Data privacy** | Patient data stays in the prompt; never used for training |

### 8.3 Vector Embeddings

Text is converted to numerical vectors (embeddings) that capture semantic meaning:

```
"heart failure with reduced ejection fraction"
    --> [0.023, -0.118, 0.456, ..., 0.089]  (384 numbers)

"systolic dysfunction with low EF"
    --> [0.025, -0.115, 0.449, ..., 0.091]  (384 numbers)

These two vectors are very similar (high cosine similarity)
because they describe the same clinical concept.
```

The system uses **BGE-small-en-v1.5** to generate 384-dimensional embeddings. All queries and documents are embedded into the same 384-dimensional space, enabling semantic search.

### 8.4 Multi-Collection RAG

The Cardiology Intelligence Agent extends basic RAG with multi-collection search:

```
User Query: "Should this HFrEF patient start an MRA?"
    |
    Embed query --> 384-dim vector
    |
    Search ALL 13 collections in parallel:
    |
    cardio_heart_failure: [result1 (score 0.91), result2 (score 0.87), ...]
    cardio_guidelines:    [result1 (score 0.88), result2 (score 0.82), ...]
    cardio_trials:        [result1 (score 0.75), result2 (score 0.71), ...]
    cardio_literature:    [result1 (score 0.68), ...]
    ... (9 more collections)
    |
    Apply collection weights:
    heart_failure results * 0.10 (or boosted to 0.25 for HF workflow)
    guidelines results * 0.10 (or boosted to 0.20)
    trials results * 0.08
    |
    Combine, rank, deduplicate
    |
    Top results assembled as LLM context
```

### 8.5 From RAG to Clinical Decision Support

The Cardiology Intelligence Agent goes beyond basic RAG by adding clinical engines:

```
Standard RAG:
  Query --> Search --> LLM --> Answer

Cardiology Intelligence Agent:
  Query --> Search ----+
                       |
         Risk Calc ----+--> LLM --> Answer + Risk Scores
                       |              + GDMT Recommendations
         GDMT Opt -----+              + Cross-Modal Triggers
                       |              + Workflow Results
         Cross-Modal --+              + Citations with Confidence
                       |
         Workflows ----+
```

This layered approach ensures that the system provides not just text answers but actionable clinical data: computed risk scores, specific medication titration plans, and genomic testing recommendations.
