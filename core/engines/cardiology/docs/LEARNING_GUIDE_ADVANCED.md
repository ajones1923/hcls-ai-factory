# Cardiology Intelligence Agent -- Learning Guide: Advanced Topics

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

Advanced cardiovascular concepts for users who have completed the Foundations guide and want deeper understanding of the clinical logic embedded in the Cardiology Intelligence Agent.

---

## Table of Contents

1. [Cardiac MRI Tissue Characterization](#1-cardiac-mri-tissue-characterization)
2. [Hemodynamic Assessment](#2-hemodynamic-assessment)
3. [Valvular Quantification](#3-valvular-quantification)
4. [Channelopathy Genetics](#4-channelopathy-genetics)
5. [Cross-Modal Imaging-Genomics Integration](#5-cross-modal-imaging-genomics-integration)
6. [Advanced Risk Stratification](#6-advanced-risk-stratification)
7. [Cardio-Oncology](#7-cardio-oncology)
8. [Multi-Collection RAG Optimization](#8-multi-collection-rag-optimization)

---

## 1. Cardiac MRI Tissue Characterization

### 1.1 Why Cardiac MRI Is Transformative

Cardiac MRI (CMR) is the only imaging modality that provides simultaneous assessment of:
- **Structure**: Chamber volumes, wall thickness, mass
- **Function**: Ejection fraction, regional wall motion, strain
- **Tissue composition**: Edema (T2), fibrosis (LGE, T1, ECV), iron overload (T2*)
- **Perfusion**: Myocardial blood flow (stress/rest)
- **Flow**: 4D flow for valvular assessment and shunt quantification

### 1.2 Late Gadolinium Enhancement (LGE)

LGE is the gold standard for detecting myocardial fibrosis and scar. Gadolinium contrast accumulates in areas of expanded extracellular space (fibrosis, necrosis, infiltration).

**LGE Patterns and Differential Diagnosis:**

| Pattern | Agent Enum | Associated Conditions | Cross-Modal Trigger? |
|---------|-----------|----------------------|---------------------|
| **Subendocardial** | `LGEPattern.SUBENDOCARDIAL` | Ischemic heart disease (coronary territory) | No (known etiology) |
| **Transmural** | `LGEPattern.TRANSMURAL` | Completed MI (full-thickness infarction) | No |
| **Mid-wall** | `LGEPattern.MID_WALL` | DCM (especially LMNA, TTN), myocarditis, sarcoidosis | Yes: DCM gene panel |
| **Epicardial** | `LGEPattern.EPICARDIAL` | Myocarditis, sarcoidosis, Chagas disease | Sometimes |
| **RV insertion point** | `LGEPattern.RV_INSERTION` | Pulmonary hypertension (mechanical stress), HCM | Depends on context |
| **Patchy** | `LGEPattern.PATCHY` | HCM, sarcoidosis, Anderson-Fabry | Yes: HCM or Fabry panel |
| **None** | `LGEPattern.NONE` | Normal, HFpEF, early disease | N/A |

**Key clinical principle:** The LGE pattern distinguishes ischemic from non-ischemic etiologies:
- **Subendocardial/transmural in a coronary territory** = ischemic (follows coronary distribution)
- **Mid-wall, epicardial, or patchy** = non-ischemic (does not follow coronary distribution)

This distinction is critical because non-ischemic patterns often warrant genetic testing (cross-modal trigger), while ischemic patterns warrant coronary evaluation.

### 1.3 T1 Mapping and Extracellular Volume (ECV)

T1 mapping quantifies the T1 relaxation time of myocardial tissue, which reflects tissue composition:

| Parameter | Normal (1.5T) | Elevated In | Reduced In |
|-----------|--------------|------------|-----------|
| **Native T1** | 950-1050 ms | Edema, fibrosis, amyloid, iron overload | Anderson-Fabry (lipid storage) |
| **Post-contrast T1** | 400-500 ms | Fibrosis (shorter post-contrast T1) | - |
| **ECV** | 25-30% | Diffuse fibrosis, amyloid, edema | Athletic heart (lower normal) |

**Why this matters for the system:**
- **Cardiac amyloidosis**: Elevated native T1 (>1100ms) + very high ECV (>40%) + diffuse LGE + clinical features triggers the amyloidosis cross-modal trigger (TTR gene testing)
- **Anderson-Fabry disease**: Paradoxically LOW native T1 (<900ms at 1.5T) due to lipid storage in myocytes, combined with LVH, triggers GLA gene testing
- **Myocarditis**: Elevated T2 (edema) + elevated T1/ECV + non-ischemic LGE = Lake Louise Criteria

### 1.4 T2 Mapping

T2 mapping detects myocardial edema (acute inflammation):

| Parameter | Normal (1.5T) | Interpretation |
|-----------|--------------|---------------|
| **T2** | ~52 ms | >55-60 ms suggests active edema/inflammation |

Used in combination with T1 and LGE for the **2018 Updated Lake Louise Criteria** for myocarditis diagnosis:
- **At least 1 T2-based criterion** (edema): Elevated T2, T2 mapping
- **Plus at least 1 T1-based criterion** (necrosis/fibrosis): Elevated T1, elevated ECV, non-ischemic LGE
- **Supportive:** Pericardial effusion, LV systolic dysfunction

### 1.5 Parametric Mapping Clinical Decision Tree

```
Unexplained Cardiomyopathy on CMR
    |
    +-- T1 elevated + ECV very high (>40%) + diffuse LGE
    |   --> AMYLOIDOSIS (trigger: TTR gene, bone scan)
    |
    +-- T1 LOW (<900ms at 1.5T) + LVH
    |   --> FABRY DISEASE (trigger: GLA gene)
    |
    +-- T2 elevated + T1 elevated + non-ischemic LGE
    |   --> MYOCARDITIS (no genetic trigger; follow Lake Louise)
    |
    +-- Mid-wall LGE + dilated LV + reduced EF
    |   --> DCM (trigger: TTN, LMNA, RBM20 gene panel)
    |
    +-- RV dilation + RV fatty replacement + arrhythmias
    |   --> ARVC (trigger: PKP2, DSP, DSG2 gene panel)
    |
    +-- Patchy LGE + FDG-PET uptake
    |   --> SARCOIDOSIS (no genetic trigger; biopsy)
    |
    +-- Asymmetric LVH + patchy LGE at RV insertion
        --> HCM (trigger: MYH7, MYBPC3 gene panel)
```

---

## 2. Hemodynamic Assessment

### 2.1 Right Heart Catheterization

Right heart catheterization (RHC) provides definitive hemodynamic data:

| Measurement | Normal Value | Clinical Significance |
|-------------|-------------|---------------------|
| **RA pressure** | 0-8 mmHg | Elevated in RV failure, tamponade, constrictive pericarditis |
| **RV pressure** | 15-30/0-8 mmHg | Elevated in PH, PS, RV failure |
| **PA pressure** | 15-30/4-12 mmHg | mPAP >20 = pulmonary hypertension (2022 definition) |
| **PCWP** | 4-12 mmHg | >15 = elevated left-sided filling pressures |
| **Cardiac output** | 4-8 L/min | <4 = low output (cardiogenic shock if + hypoperfusion) |
| **Cardiac index** | >2.2 L/min/m2 | <2.2 = reduced; <1.8 = severely reduced |
| **PVR** | <3 Wood units | >3 WU = elevated; >5 WU = severe PH |
| **SVR** | 800-1200 dyn*s/cm5 | Elevated in cardiogenic shock; reduced in septic shock |

### 2.2 Hemodynamic Profiles in Heart Failure

The hemodynamic assessment workflow classifies patients into profiles:

```
                           PCWP
                    Low (<15)       High (>15)
                 +------------+----------------+
    CI >2.2      | "Warm-Dry" |   "Warm-Wet"   |
    (adequate    | Profile A  |   Profile B     |
     perfusion)  | (compensated)| (congested)   |
                 +------------+----------------+
    CI <2.2      | "Cold-Dry" |   "Cold-Wet"   |
    (poor        | Profile L  |   Profile C     |
     perfusion)  | (hypovolemic)| (most severe) |
                 +------------+----------------+
```

**Profile-specific management:**
- **Warm-Dry (A):** Optimize oral GDMT; no acute intervention needed
- **Warm-Wet (B):** IV diuretics to decongest; continue/optimize GDMT
- **Cold-Dry (L):** Cautious volume challenge; consider inotropes
- **Cold-Wet (C):** Inotropes + diuretics; consider mechanical support (IABP, Impella)

### 2.3 Pulmonary Hypertension Classification

The 2022 ESC/ERS definition uses hemodynamic criteria:

| Type | Hemodynamic Definition | Common Causes |
|------|----------------------|---------------|
| **Pre-capillary PH** | mPAP >20, PCWP <=15, PVR >2 WU | PAH, CTEPH, lung disease |
| **Post-capillary PH** | mPAP >20, PCWP >15 | Left heart disease (HFrEF, HFpEF, VHD) |
| **Combined pre+post** | mPAP >20, PCWP >15, PVR >2 WU | Advanced HF with reactive component |

---

## 3. Valvular Quantification

### 3.1 Aortic Stenosis Grading

The Valvular Disease workflow uses ASE guidelines for severity grading:

| Parameter | Mild | Moderate | Severe |
|-----------|------|----------|--------|
| **Vmax** | 2.0-2.9 m/s | 3.0-3.9 m/s | >=4.0 m/s |
| **Mean gradient** | <20 mmHg | 20-39 mmHg | >=40 mmHg |
| **AVA** | >1.5 cm2 | 1.0-1.5 cm2 | <1.0 cm2 |
| **DVI** | >0.35 | 0.25-0.35 | <0.25 |

**Discordant grading:** When parameters disagree (e.g., low gradient but small AVA), the system considers:
- **Low-flow, low-gradient severe AS:** LVEF <50% + AVA <1.0 + mean gradient <40 -> dobutamine stress echo
- **Paradoxical low-flow:** LVEF >=50% but small stroke volume (stroke volume index <35 mL/m2)

### 3.2 Mitral Regurgitation Quantification

**Primary MR (degenerative):**

| Parameter | Mild | Moderate | Severe |
|-----------|------|----------|--------|
| **Vena contracta** | <0.3 cm | 0.3-0.69 cm | >=0.7 cm |
| **Regurgitant volume** | <30 mL | 30-59 mL | >=60 mL |
| **Regurgitant fraction** | <30% | 30-49% | >=50% |
| **ERO** | <0.20 cm2 | 0.20-0.39 cm2 | >=0.40 cm2 |

**Secondary MR (functional):** Different thresholds apply:
- Severe secondary MR: ERO >=0.20 cm2, RVol >=30 mL (lower thresholds than primary)

### 3.3 Intervention Decision Algorithms

The system implements ACC/AHA 2020 VHD guideline decision logic:

**Aortic Stenosis:**
```
Severe AS (Vmax >=4.0, MG >=40, AVA <1.0)?
    |
    Yes --> Symptomatic?
    |           |
    |           Yes --> AVR indicated (Class I)
    |           |       Age >=65 or high surgical risk --> TAVR
    |           |       Age <65 and low surgical risk --> SAVR
    |           |
    |           No --> LVEF <50%?
    |                   |
    |                   Yes --> AVR indicated (Class I)
    |                   |
    |                   No --> Very severe (Vmax >=5.0)?
    |                           |
    |                           Yes --> AVR reasonable (Class IIa)
    |                           No --> Surveillance
    |
    No --> Not severe; surveillance
```

**Mitral Regurgitation:**
```
Severe Primary MR?
    |
    Yes --> Symptomatic?
    |           |
    |           Yes --> MV surgery indicated (Class I)
    |           |       Repair preferred over replacement
    |           |
    |           No --> LVEF <60% or LVESD >40mm?
    |                   |
    |                   Yes --> MV surgery indicated (Class I)
    |                   No --> Surveillance (consider if repair likelihood >95%)
    |
Severe Secondary MR?
    |
    Yes --> Optimize GDMT + CRT first
            Persistent symptoms? --> Consider TEER (MitraClip)
```

---

## 4. Channelopathy Genetics

### 4.1 Overview

Channelopathies are inherited disorders of cardiac ion channels that predispose to life-threatening arrhythmias in structurally normal hearts. The knowledge graph contains detailed entries for 7 channelopathy genes.

### 4.2 Long QT Syndrome (LQTS)

Three major types, each with distinct triggers and gene-specific therapy:

| Type | Gene | Current | Trigger | Specific Therapy |
|------|------|---------|---------|-----------------|
| **LQT1** | KCNQ1 | IKs | Exercise, swimming | Beta-blocker (nadolol); avoid competitive sports |
| **LQT2** | KCNH2 | IKr | Auditory stimuli, emotional stress, postpartum | Beta-blocker; avoid QT-prolonging drugs; supplemental K+ |
| **LQT3** | SCN5A | INa (gain of function) | Rest, sleep | Mexiletine (sodium channel blocker); consider ICD |

**Why gene-specific therapy matters:** LQT3 patients respond to mexiletine (reduces late sodium current) but not beta-blockers. The cross-modal trigger system identifies LQTS gene testing candidates from ECG findings (QTc >480ms).

### 4.3 Brugada Syndrome

- **Gene:** SCN5A (loss of function) -- only ~20% of cases are genotype-positive
- **ECG hallmark:** Type 1 Brugada pattern (coved ST elevation >=2mm in V1-V3)
- **Triggers:** Fever, sodium channel blockers, vagal stimulation
- **Risk stratification:** Prior cardiac arrest, spontaneous type 1 pattern, syncope
- **Treatment:** ICD for high-risk; quinidine for recurrent VF storms; fever management

### 4.4 CPVT

- **Genes:** RYR2 (AD, ~60%) or CASQ2 (AR)
- **Hallmark:** Bidirectional or polymorphic VT with exercise/catecholamines in a structurally normal heart
- **Treatment:** Nadolol (preferred beta-blocker), flecainide as add-on, strict exercise restriction
- **Key point:** Normal resting ECG; exercise stress test needed for diagnosis

### 4.5 Genotype-Phenotype Correlations in the Knowledge Graph

The knowledge graph cross-references genes to conditions:

```
SCN5A --> Brugada (loss of function)
      --> LQT3 (gain of function)
      --> Progressive cardiac conduction disease
      --> Sick sinus syndrome
      --> DCM
      --> AF

Different variant types in the SAME gene cause different diseases.
The system's cross-modal engine identifies which gene panel to recommend
based on the specific clinical phenotype.
```

---

## 5. Cross-Modal Imaging-Genomics Integration

### 5.1 The Cross-Modal Paradigm

The Cardiology Intelligence Agent implements a novel approach to cardiovascular precision medicine: automated imaging-to-genomics triggers. This mirrors the clinical workflow where a cardiologist identifies an imaging finding that warrants genetic evaluation, but automates the pattern recognition.

### 5.2 Detailed Trigger Workflows

**Trigger 1: Unexplained LVH (Wall Thickness >=15mm)**

```
Echocardiogram or CMR shows wall thickness >=15mm
    |
    Exclude common causes: Hypertension? Aortic stenosis? Athletic heart?
    |
    If unexplained:
    |
    +-- Gene panel: MYH7, MYBPC3, TNNT2, TNNI3, TPM1 (HCM sarcomeric)
    |               GLA (Fabry disease -- X-linked, check alpha-gal A enzyme)
    |               LAMP2 (Danon disease)
    |               PRKAG2 (glycogen storage with WPW)
    |
    +-- Query genomic_evidence for known pathogenic variants
    |
    +-- Clinical action:
        - If sarcomeric gene positive: cascade family screening
        - If GLA positive: enzyme replacement or chaperone therapy
        - If LAMP2 positive: evaluate for transplant (aggressive course)
```

**Trigger 2: Non-Ischemic LGE with DCM**

```
Cardiac MRI shows mid-wall LGE + dilated LV + reduced EF
    |
    Pattern is NON-ischemic (not following coronary territory)
    |
    +-- Gene panel: TTN, LMNA, RBM20, MYH7, DSP, FLNC, BAG3
    |
    +-- Special attention to LMNA:
    |   - Mid-wall septal LGE is highly specific for LMNA mutations
    |   - LMNA cardiomyopathy has high SCD risk
    |   - Early ICD consideration (lower LVEF threshold than standard DCM)
    |   - Conduction disease often present (AV block, sinus node dysfunction)
    |
    +-- If TTN truncating variant (TTNtv) in A-band:
        - Most common DCM gene (~20-25%)
        - Generally better prognosis than LMNA
        - Standard GDMT + standard ICD criteria
```

**Trigger 3: Aortic Root Dilation in Young Patient**

```
Aortic root >=4.0cm (or Z-score >=2) in patient <50 years
    |
    +-- Gene panel: FBN1, TGFBR1, TGFBR2, SMAD3, ACTA2, MYH11, COL3A1
    |
    +-- Gene-specific surgical thresholds:
    |   FBN1 (Marfan): Surgery at 5.0 cm
    |   TGFBR1/2 (Loeys-Dietz): Surgery at 4.0-4.2 cm (more aggressive)
    |   COL3A1 (vascular EDS): Avoid elective surgery if possible (tissue fragility)
    |   ACTA2: Surgery at 4.5-5.0 cm depending on variant
    |
    +-- This demonstrates the clinical impact of genomic data on surgical timing:
        Knowing the gene changes the surgical threshold by up to 1.5 cm.
```

### 5.3 Integration with HCLS AI Factory Genomics Pipeline

The cross-modal engine queries the shared `genomic_evidence` collection, which is populated by the HCLS AI Factory's genomics pipeline:

```
HCLS AI Factory Genomics Pipeline
    |
    FASTQ --> BWA-MEM2 --> BAM --> DeepVariant --> VCF
    |
    VCF --> Annotation (ClinVar, AlphaMissense)
    |
    Annotated variants --> BGE-small-en-v1.5 embedding
    |
    3.5M variant vectors in genomic_evidence collection
    |
    <-- Queried by Cardiology Intelligence Agent cross-modal engine
```

This creates a seamless bridge between the genomics pipeline and the cardiology agent, enabling true precision cardiovascular medicine.

---

## 6. Advanced Risk Stratification

### 6.1 Beyond ASCVD: Risk Enhancers

The 2018 ACC/AHA Cholesterol Guidelines introduced "risk-enhancing factors" for patients in the borderline or intermediate risk category:

| Risk Enhancer | Threshold | Agent Integration |
|---------------|-----------|------------------|
| Family history of premature ASCVD | Male <55, Female <65 in first-degree relative | Knowledge graph cross-reference |
| Persistently elevated LDL-C | >=160 mg/dL | Biomarker knowledge (LDL-C entry) |
| Elevated Lp(a) | >=50 mg/dL or >=125 nmol/L | Biomarker knowledge + LPA gene |
| Elevated hsCRP | >=2.0 mg/L | Biomarker knowledge (hsCRP entry) |
| Elevated ApoB | >=130 mg/dL | Biomarker knowledge (ApoB entry) |
| Metabolic syndrome | 3+ criteria | Patient context evaluation |
| CKD | eGFR 15-59 | Biomarker knowledge (Creatinine/eGFR) |
| Chronic inflammatory conditions | RA, lupus, psoriasis, HIV | Patient history |
| South Asian ancestry | - | Demographics |
| Preeclampsia/premature menopause | - | Patient history |
| Ankle-brachial index | <0.9 | Patient context |
| CAC score | >=100 Agatston units | Imaging data |

### 6.2 MAGGIC Score Deep Dive

The MAGGIC score integrates 13 variables with interaction terms:

**Age scoring:**
```
Age points = (age - 55) / 5 * coefficient

Interaction: Age effect is LARGER when EF is higher (paradox)
- At EF 20%: each 5-year increment adds ~2 points
- At EF 40%: each 5-year increment adds ~3 points

This reflects that in severely reduced EF, age contributes
less to incremental risk than the EF itself.
```

**LVEF scoring:**
```
EF points = (30 - EF) / 5 * coefficient (if EF <30)
           OR more gradual scoring if EF 30-40

Lower EF = dramatically more points
The relationship is non-linear: going from EF 30 to 20 adds
more risk than going from 40 to 30.
```

### 6.3 HCM Sudden Cardiac Death Risk Stratification

The system evaluates HCM patients for SCD risk using ACC/AHA criteria:

| Risk Factor | ICD Consideration |
|-------------|------------------|
| Prior cardiac arrest or sustained VT | ICD recommended (Class I) |
| Family history of SCD in first-degree relative | Factor in risk discussion |
| Unexplained syncope | Major risk factor |
| Massive LVH (>=30mm) | Major risk factor |
| NSVT on Holter monitoring | Major risk factor |
| Abnormal BP response to exercise | Additional factor |
| Extensive LGE on CMR (>=15% LV mass) | Additional factor |
| LV apical aneurysm | Additional factor |
| LVEF <50% (end-stage HCM) | ICD recommended |

### 6.4 SCAI Cardiogenic Shock Classification

For critically ill patients, the system references the SCAI shock classification:

| Stage | Hemodynamics | Lactate | Management |
|-------|-------------|---------|-----------|
| A (At risk) | Normal | Normal | Monitor, optimize |
| B (Beginning) | SBP <90 or MAP <60 | Normal | Vasopressors, IV fluids |
| C (Classic) | CI <2.2, PCWP >15 | 2-5 mmol/L | Inotropes, MCS evaluation |
| D (Deteriorating) | Worsening on initial interventions | Rising >5 | Escalate MCS (Impella, ECMO) |
| E (Extremis) | PEA, refractory VT, multiorgan failure | >8 | ECPR, emergent MCS |

---

## 7. Cardio-Oncology

### 7.1 The Cardio-Oncology Challenge

Cancer therapies can cause cardiac damage through multiple mechanisms. The Cardio-Oncology workflow monitors for cancer therapy-related cardiac dysfunction (CTRCD).

### 7.2 Cardiotoxic Agent Classification

| Agent | Risk Level (Agent Enum) | Mechanism | Monitoring |
|-------|------------------------|-----------|-----------|
| **Doxorubicin** (anthracycline) | `CardiotoxicityRisk.HIGH` | Dose-dependent myocardial cell death, ROS generation | Echo + GLS q2 cycles; troponin each cycle |
| **Trastuzumab** (HER2-targeted) | `CardiotoxicityRisk.MODERATE` | Type II dysfunction (reversible); no cell death | Echo + GLS q3 months |
| **Pembrolizumab** (ICI) | `CardiotoxicityRisk.MODERATE` | Immune-mediated myocarditis (rare but fulminant) | Troponin; ECG; echo if symptomatic |
| **5-FU / Capecitabine** | `CardiotoxicityRisk.MODERATE` | Coronary vasospasm | ECG; symptoms (chest pain) |
| **Ibrutinib** (BTK inhibitor) | `CardiotoxicityRisk.MODERATE` | AF (~10%), VT, hypertension | ECG; BP monitoring |
| **VEGF inhibitors** | `CardiotoxicityRisk.MODERATE` | Hypertension, arterial thromboembolism | BP monitoring q2 weeks initially |
| **Radiation (chest)** | `CardiotoxicityRisk.HIGH` | Pericarditis, CAD, valvular disease (delayed) | Long-term surveillance (5-10+ years post-RT) |

### 7.3 CTRCD Definitions (ESC 2022)

The system uses the ESC 2022 cardio-oncology guideline definitions:

| Category | Definition | System Response |
|----------|-----------|----------------|
| **Symptomatic CTRCD** | New HF symptoms + LVEF decline | Flag as HIGH severity; cardiology consultation |
| **Asymptomatic CTRCD** | LVEF decline >10% to below 50% | Flag as MODERATE; consider cardioprotection |
| **Subclinical CTRCD** | GLS decline >15% from baseline (relative) | Flag as MODERATE; early intervention opportunity |
| **Biomarker CTRCD** | Persistent troponin elevation above URL | Flag as MODERATE; intensify imaging surveillance |

### 7.4 The GLS Advantage

Global Longitudinal Strain (GLS) detects subclinical myocardial dysfunction before LVEF decline:

```
Timeline of Cardiotoxicity Detection:

      Molecular damage  GLS decline    LVEF decline    HF symptoms
           |               |               |               |
     Day 0          Week 4-8        Month 3-6        Month 6-12+
           |               |               |               |
    Troponin rise   15% relative     >10% absolute    Clinical HF
                    GLS decline      LVEF decline
           |               |               |               |
    EARLIEST -----> EARLY ---------> MODERATE -------> LATE
    detection       detection        detection         detection

The Cardio-Oncology workflow monitors GLS as the primary
early detection biomarker, flagging a >15% relative decline
as subclinical CTRCD.
```

### 7.5 Cardioprotective Strategies

The system recommends cardioprotection based on risk level:

| Risk Level | Pre-Treatment | During Treatment | Post-Treatment |
|-----------|---------------|-----------------|---------------|
| **Low** | Baseline echo | Symptom-based monitoring | Follow-up echo at 12 months |
| **Moderate** | Baseline echo + GLS + biomarkers | Echo + GLS q3 cycles; troponin q1 cycle | Echo + GLS at 3, 12 months |
| **High** | Baseline echo + GLS + biomarkers; consider ACEi/BB prophylaxis | Echo + GLS q2 cycles; troponin each cycle; dexrazoxane if anthracycline dose >300 mg/m2 | Echo + GLS at 3, 6, 12 months; lifelong surveillance |

---

## 8. Multi-Collection RAG Optimization

### 8.1 Collection Weight Tuning

The default collection weights were calibrated through iterative testing against clinical vignettes. Advanced users may need to tune weights for specific use cases:

**Tuning principles:**
1. **Higher weight = more influence** on final ranked results
2. **Workflow-specific boosts** are applied multiplicatively (default weights * boost factor)
3. **Weights must sum to ~1.0** (tolerance 0.05, enforced by settings validation)
4. **Minimum weight = 0.02** to prevent complete exclusion of any collection

### 8.2 Score Threshold Optimization

The `SCORE_THRESHOLD` setting (default 0.4) controls the minimum similarity score for a result to be included:

| Threshold | Behavior | Use Case |
|-----------|----------|----------|
| 0.3 | Broad retrieval; may include less relevant results | Exploratory queries, rare conditions |
| 0.4 | Balanced (default) | General clinical queries |
| 0.5 | Precise retrieval; fewer but more relevant results | Specific guideline lookups |
| 0.6 | Very strict; may miss relevant content | When precision is critical |

### 8.3 Query Decomposition Strategies

The query expansion module selects from four search strategies:

| Strategy | When Selected | Behavior |
|----------|--------------|----------|
| **Broad** | General questions, multi-topic queries | Search all collections with default weights |
| **Targeted** | Specific condition or drug query | Boost relevant collections 2-3x |
| **Comparative** | "Drug A vs Drug B" or "TAVR vs SAVR" | Search for both entities, present side-by-side |
| **Clinical** | Patient-specific management question | Boost guidelines + heart_failure/relevant workflow |

### 8.4 Embedding Space Analysis

All collections share the same BGE-small-en-v1.5 embedding space (384 dimensions). This creates a unified semantic space where:

- Queries about "heart failure with reduced ejection fraction" are close to documents about "HFrEF GDMT"
- But also close to documents about "systolic dysfunction management" (different terminology, same concept)
- Entity aliases in the query expansion module bridge terminology gaps before embedding

**Embedding quality considerations:**
- BGE-small-en-v1.5 was trained on general English text; it handles medical terminology well but may underperform on highly specialized abbreviations
- The entity alias expansion step mitigates this by expanding abbreviations before embedding
- For production deployments with large collections, consider fine-tuning the embedding model on cardiovascular text

### 8.5 Citation Confidence Calibration

The citation scoring system uses two thresholds:

| Level | Score Range | Display | Interpretation |
|-------|-----------|---------|---------------|
| **High confidence** | >= 0.75 | Strong recommendation | Source directly addresses the query |
| **Medium confidence** | >= 0.60 | Supporting evidence | Source is relevant but may not be exact match |
| **Standard** | 0.40-0.59 | Additional context | Source provides background; should not be primary citation |

**Calibration guidance:**
- If too many results are marked "high confidence," raise `CITATION_HIGH_THRESHOLD`
- If too few results pass citation scoring, lower `CITATION_MEDIUM_THRESHOLD`
- Monitor the `cardio_queries_total` and citation distributions via Prometheus metrics

### 8.6 Conversation Context Management

The sliding window of 3 conversation turns enables contextual follow-up:

```
Turn 1: "What GDMT does this patient need?"
  -> Full GDMT analysis

Turn 2: "What about adding an MRA?"
  -> System remembers patient context from Turn 1
  -> Focuses on MRA specifically (K+, eGFR, contraindications)

Turn 3: "Are there any relevant clinical trials?"
  -> System searches cardio_trials collection
  -> Contextualizes results to the HFrEF + MRA context from Turns 1-2

Turn 4: (New question)
  -> Turn 1 falls out of the window
  -> Only Turns 2-3-4 retained
```

**Configuration:** `MAX_CONVERSATION_CONTEXT = 3` in settings. Increase for longer conversations at the cost of larger LLM prompts. Decrease for faster responses with less context.

### 8.7 Performance Optimization Checklist

For production deployments:

1. **Index type**: IVF_FLAT is the default; switch to HNSW for collections with >500K vectors (faster search, more memory)
2. **Batch embedding**: Use `EMBEDDING_BATCH_SIZE=32` (default); increase to 64 for faster bulk ingest
3. **Parallel search**: The RAG engine searches all 13 collections concurrently; ensure Milvus has sufficient query node capacity
4. **LLM caching**: Consider caching LLM responses for identical queries (risk calculator results are deterministic; LLM responses are not)
5. **Compaction**: Run Milvus compaction monthly on actively ingested collections to reclaim space from deleted/updated segments
