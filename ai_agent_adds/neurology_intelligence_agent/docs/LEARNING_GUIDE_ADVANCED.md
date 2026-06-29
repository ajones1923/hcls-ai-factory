# Neurology Intelligence Agent -- Learning Guide: Advanced Topics

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## Purpose

This advanced guide covers the specialized clinical frameworks, trial criteria, assessment methodologies, and pattern recognition skills implemented within the Neurology Intelligence Agent. It is intended for those who have completed the Foundations guide and want to understand the clinical reasoning encoded in the agent's workflows and scale calculators.

---

## 1. The ATN Biomarker Framework in Detail

### 1.1 Biological Definition of Alzheimer's Disease

The 2018 NIA-AA Research Framework fundamentally shifted Alzheimer's from a clinical diagnosis to a biological one. AD is now defined by the presence of biomarkers, not symptoms alone.

### 1.2 Biomarker Modalities

| Marker | Fluid-Based | Imaging-Based | Cut-offs (Approximate) |
|---|---|---|---|
| **A (Amyloid)** | CSF Abeta42 < 500 pg/mL; Abeta42/40 ratio < 0.058 | Amyloid PET (Florbetapir, Florbetaben, Flutemetamol): visual read positive | Plasma p-tau217 emerging as blood-based screen |
| **T (Tau)** | CSF p-tau181 > 27 pg/mL; p-tau217 > 0.38 pg/mL | Tau PET (Flortaucipir): Braak stage I-VI uptake pattern | Regional tau correlates with symptom domain |
| **N (Neurodegeneration)** | CSF total tau > 400 pg/mL; serum NfL > 19 pg/mL | FDG-PET (temporo-parietal hypometabolism); MRI volumetrics (hippocampal volume) | NfL is non-specific (elevated in ALS, MS, stroke) |

### 1.3 Clinical Staging Integration

| ATN Profile | Biological Interpretation | Clinical Implication |
|---|---|---|
| A-T-N- | Normal AD biomarkers | No AD pathology; evaluate other causes |
| A+T-N- | Alzheimer's pathologic change | Preclinical AD; biomarker-only stage |
| A+T+N- | AD (early biological) | Biological AD; may be cognitively normal or MCI |
| A+T+N+ | AD (full biological) | Full AD continuum; correlates with clinical dementia |
| A+T-N+ | Alzheimer's + non-tau neurodegeneration | AD amyloid with concomitant pathology (vascular, Lewy body) |
| A-T+N- | Non-AD pathologic change | Primary tauopathy (PSP, CBD, FTD-tau); SNAP |
| A-T-N+ | Non-AD neurodegeneration | Vascular, Lewy body, hippocampal sclerosis |
| A-T+N+ | Non-AD with neurodegeneration | Primary tauopathy with advanced disease |

### 1.4 Anti-Amyloid Therapy Eligibility

Current eligibility criteria (based on CLARITY AD and TRAILBLAZER-ALZ 2):
- Age 50-90 (some variation by trial)
- MoCA typically 18-26 (MCI to mild dementia)
- Amyloid-positive on PET or CSF
- No more than 4 microhemorrhages on MRI
- No macrohemorrhage history
- APOE e4/e4 carriers have higher ARIA risk -- requires informed consent

### 1.5 ARIA Monitoring Protocol

- Baseline MRI before treatment initiation
- MRI at weeks 7, 14, 52 (minimum for lecanemab)
- ARIA-E (edema): Usually asymptomatic; hold treatment until resolved
- ARIA-H (hemorrhage): Microhemorrhages may or may not require treatment hold
- Symptomatic ARIA: Hold treatment, close follow-up

---

## 2. Extended Stroke Treatment Windows: DAWN and DEFUSE-3

### 2.1 The Paradigm Shift

Before 2018, mechanical thrombectomy was limited to 6 hours from symptom onset. DAWN and DEFUSE-3 trials extended the window to 24 hours using imaging selection, saving thousands of additional patients annually.

### 2.2 DAWN Trial Criteria (6-24 hours)

**Key principle:** Clinical-imaging mismatch -- the clinical deficit is disproportionately worse than the infarct core size.

| Patient Group | NIHSS Requirement | Infarct Core Limit |
|---|---|---|
| Age >= 80 | NIHSS >= 10 | Core volume < 21 mL |
| Age < 80, moderate deficit | NIHSS >= 10 | Core volume < 31 mL |
| Age < 80, severe deficit | NIHSS >= 20 | Core volume 31-51 mL |

**Additional requirements:**
- ICA or M1 MCA occlusion
- Pre-stroke mRS 0-1
- Core volume measured by CTP or DWI-MRI

**Results:** 49% vs 13% functional independence (mRS 0-2) at 90 days. NNT = 2.8.

### 2.3 DEFUSE-3 Trial Criteria (6-16 hours)

**Key principle:** Perfusion-diffusion mismatch -- significant salvageable tissue (penumbra) exists beyond the infarct core.

| Criterion | Requirement |
|---|---|
| Vessel occlusion | ICA or M1 MCA |
| Age | 18-90 |
| NIHSS | >= 6 |
| Pre-stroke mRS | 0-2 |
| Infarct core volume | < 70 mL |
| Mismatch ratio | >= 1.8 |
| Mismatch volume | >= 15 mL |
| Tmax > 6 seconds | Defines penumbra on CTP/MR perfusion |

**Results:** 45% vs 17% functional independence. Mismatch volume predicted treatment benefit.

### 2.4 Practical Integration in the Agent

The stroke triage workflow evaluates:
1. Time from last known well
2. NIHSS score (calculated or provided)
3. ASPECTS (from CT) or core volume (from CTP/DWI)
4. Vessel occlusion status (from CTA)
5. Pre-stroke functional status (mRS)
6. Maps patient to appropriate protocol (standard window, DAWN, DEFUSE-3, or not eligible)

---

## 3. NEDA-3 and MS Treatment Targets

### 3.1 NEDA-3 Definition

"No Evidence of Disease Activity" is the composite treatment target for relapsing MS:

| Component | Measurement | Definition of "Activity" |
|---|---|---|
| **Clinical relapses** | Clinician-assessed | Any new or recurrent neurological symptom lasting > 24 hours |
| **MRI activity** | Annual brain MRI | Any new or enlarging T2 lesion OR Gd-enhancing lesion |
| **Disability progression** | EDSS | >= 1.0 point increase sustained at 3-6 months (baseline <= 5.5) OR >= 0.5 point increase (baseline >= 6.0) |

### 3.2 NEDA-4 Extension

Adds **brain volume loss** (annualized brain volume loss > 0.4% is abnormal). More sensitive but less routinely measured.

### 3.3 DMT Escalation Decision Framework

```
Patient on DMT
  |
  +-- Meeting NEDA-3? ---> Yes: Continue current DMT
  |
  +-- No: Evidence of disease activity
       |
       +-- Relapses? ------> 1 mild relapse: Monitor closely
       |                      2+ relapses or severe: Escalate
       |
       +-- MRI activity? --> 1-2 new T2 lesions: May monitor
       |                      3+ new lesions or Gd+: Escalate
       |
       +-- EDSS progression? -> Confirmed: Escalate
       |
       +-- Combined: Any 2+ of above: Strong escalation indication
       |
       +-- Escalation path:
            Platform --> Moderate efficacy --> High efficacy
            Consider: JCV status, pregnancy planning, comorbidities
```

### 3.4 JCV/PML Risk Stratification (Natalizumab)

| JCV Antibody Index | Prior Immunosuppression | PML Risk (per 1000) |
|---|---|---|
| Negative | No | ~0.07 |
| Low (< 0.9) | No | ~0.1 |
| Moderate (0.9-1.5) | No | ~1.0 |
| High (> 1.5) | No | ~6.0 |
| High (> 1.5) | Yes | ~11.0 |

Risk increases significantly after 24 months of treatment. Extended interval dosing (every 6 weeks instead of 4) may reduce risk.

---

## 4. ACLS Stroke Algorithm Integration

### 4.1 Prehospital Stroke Recognition

**BE-FAST mnemonic:**
- **B**alance: Sudden loss of balance
- **E**yes: Sudden vision change
- **F**ace: Facial droop
- **A**rm: Arm weakness
- **S**peech: Speech difficulty
- **T**ime: Time to call 911

### 4.2 Emergency Department Stroke Protocol

```
Stroke Alert Activation
  |
  +-- Time zero: Last known well (LKW) documented
  |
  +-- 10 min: CT head (rule out hemorrhage)
  |
  +-- 15 min: NIHSS scored
  |
  +-- 20 min: Labs drawn (glucose, coags, CBC)
  |
  +-- 25 min: CT angiography (if NIHSS >= 6 or LVO suspected)
  |
  +-- tPA Decision:
  |     LKW < 4.5h + NIHSS >= 1 + no contraindications
  |     --> IV alteplase 0.9 mg/kg (max 90 mg)
  |     --> 10% bolus, 90% over 60 minutes
  |     --> Target door-to-needle < 60 minutes
  |
  +-- Thrombectomy Decision:
  |     LVO confirmed (ICA/M1) + NIHSS >= 6 + ASPECTS >= 6
  |     --> Mechanical thrombectomy
  |     --> Target door-to-groin < 90 minutes
  |
  +-- Extended Window (6-24h):
        DAWN or DEFUSE-3 criteria met
        --> CTP or MRI perfusion for core/penumbra
        --> Thrombectomy if eligible
```

### 4.3 Post-tPA Monitoring

- Blood pressure: < 180/105 for 24 hours
- Neurological checks: Every 15 minutes for 2 hours, then q30min for 6 hours, then hourly for 16 hours
- No antiplatelets or anticoagulants for 24 hours
- CT head at 24 hours before starting antithrombotic therapy
- Watch for signs of hemorrhagic transformation: sudden neurological worsening, headache, nausea/vomiting

---

## 5. EMG/NCS Pattern Recognition

### 5.1 Fundamentals

**Nerve Conduction Studies (NCS):** Electrical stimulation of a nerve with recording of the response. Measures:
- **Amplitude:** Reflects number of functioning axons (reduced in axonal disease)
- **Conduction velocity:** Reflects myelin integrity (slowed in demyelinating disease)
- **Distal latency:** Time from stimulus to response
- **F-wave latency:** Proximal nerve conduction (tests nerve roots)

**Electromyography (EMG):** Needle electrode inserted into muscle. Evaluates:
- **Insertional activity:** Normal, increased (denervation), decreased (end-stage)
- **Spontaneous activity:** Fibrillation potentials and positive sharp waves (active denervation)
- **Motor unit potential (MUP) morphology:** Large/polyphasic (neurogenic), small/short (myopathic)
- **Recruitment pattern:** Reduced (neurogenic), early (myopathic)

### 5.2 Pattern Recognition Guide

| Pattern | NCS Findings | EMG Findings | Differential |
|---|---|---|---|
| **Axonal neuropathy** | Low amplitudes, normal velocities | Fibrillations in distal muscles, large MUAPs | Diabetic, toxic, hereditary (CMT2) |
| **Demyelinating neuropathy** | Slow velocities, prolonged distal latencies, conduction block | May have secondary axonal changes | CIDP, GBS (AIDP), CMT1 |
| **Motor neuron disease** | Normal sensory NCS; may have low CMAPs | Diffuse active denervation, fasciculations, large MUAPs | ALS, SMA, Kennedy disease |
| **NMJ (postsynaptic)** | Decremental response on RNS at 2-3 Hz | Normal or minimal changes | Myasthenia gravis (AChR+, MuSK+) |
| **NMJ (presynaptic)** | Low CMAP amplitudes; incremental response on rapid RNS | -- | Lambert-Eaton (VGCC antibodies) |
| **Myopathy** | Normal NCS or mildly reduced CMAPs | Small, short, polyphasic MUAPs; early recruitment | Inflammatory myopathy, dystrophy |
| **Radiculopathy** | Normal NCS (or low CMAP if severe) | Fibrillations in myotomal distribution, paraspinal muscles | Disc herniation, spondylosis |

### 5.3 Critical Diagnostic Distinctions

**GBS vs. CIDP:**
- GBS: Acute onset (< 4 weeks to nadir), monophasic, post-infectious
- CIDP: Chronic (> 8 weeks), relapsing or progressive, no trigger
- NCS: Both demyelinating, but CIDP shows more uniform slowing

**ALS Diagnosis (El Escorial Criteria):**
- Definite: UMN + LMN signs in 3+ regions (bulbar, cervical, thoracic, lumbosacral)
- Probable: UMN + LMN signs in 2+ regions with UMN rostral to LMN
- EMG: Active denervation in 3+ regions with normal sensory NCS is highly suggestive

---

## 6. DBS Candidacy Assessment

### 6.1 Indications for DBS in Parkinson's Disease

DBS is considered when medical therapy is inadequate despite optimization:

**CAPSIT-PD (Core Assessment Program for Surgical Interventional Therapies) criteria:**

| Criterion | Requirement |
|---|---|
| Diagnosis | Idiopathic PD (not atypical parkinsonism) |
| Disease duration | >= 5 years |
| Levodopa response | >= 30% improvement in UPDRS Part III on vs. off |
| Motor complications | Disabling dyskinesias OR motor fluctuations |
| Cognition | MoCA >= 26 (no dementia) |
| Psychiatric | No active untreated psychiatric disease |
| Age | Typically < 70-75 (relative, not absolute) |
| MRI | No structural lesions contraindicating surgery |

### 6.2 Target Selection

| Target | Best For | Advantages | Considerations |
|---|---|---|---|
| **STN (subthalamic nucleus)** | Broadly effective PD | Reduces medication needs by ~50% | May affect mood, speech |
| **GPi (globus pallidus interna)** | Dyskinesia-predominant | Directly suppresses dyskinesias | Less medication reduction |
| **VIM thalamus** | Essential tremor, PD tremor-dominant | Most effective for tremor | Does not help bradykinesia |

### 6.3 Emerging Alternatives

- **Focused ultrasound thalamotomy:** Non-invasive, unilateral tremor control. MRI-guided.
- **Adaptive DBS:** Closed-loop stimulation that adjusts in real-time based on brain signals.

---

## 7. Epilepsy Surgical Evaluation

### 7.1 Pre-Surgical Workup

```
Drug-Resistant Epilepsy (failed 2+ ASMs)
  |
  +-- Phase I: Non-invasive evaluation
  |     +-- Video-EEG monitoring (capture habitual seizures)
  |     +-- 3T MRI with epilepsy protocol
  |     +-- Neuropsychological testing
  |     +-- PET (FDG): Interictal hypometabolism
  |     +-- ictal SPECT (if available)
  |     +-- MEG (if available)
  |
  +-- Concordance analysis:
  |     EEG focus + MRI lesion + PET hypometabolism + semiology
  |     All point to same region? --> Proceed to surgery
  |
  +-- Discordant or non-lesional?
        +-- Phase II: Invasive monitoring
              +-- Stereo-EEG (SEEG): Depth electrodes
              +-- Subdural grid electrodes (less common now)
              +-- Functional mapping (language, motor)
              +-- Identify seizure onset zone
```

### 7.2 Surgical Options and Outcomes

| Procedure | Indication | Seizure Freedom Rate | Engel Class I |
|---|---|---|---|
| Anterior temporal lobectomy | TLE with MTS | 60-80% | 65-75% |
| Selective amygdalohippocampectomy | TLE, sparing lateral temporal | 50-70% | 55-65% |
| Lesionectomy | Focal lesion (FCD, cavernoma, tumor) | 60-80% | 60-75% |
| LITT (laser ablation) | MTS, focal lesions, hypothalamic hamartoma | 50-65% | 50-60% |
| Corpus callosotomy | Generalized/drop attacks (LGS) | Palliative; reduces drop attacks 60-80% | Rarely class I |
| Hemispherectomy | Hemispheric pathology (Rasmussen, Sturge-Weber) | 60-80% | 60-70% |

### 7.3 Engel Classification (Surgical Outcomes)

| Class | Description |
|---|---|
| I | Free of disabling seizures |
| II | Rare disabling seizures |
| III | Worthwhile improvement (> 75% reduction) |
| IV | No worthwhile improvement |

---

## 8. Neuro-Oncology Molecular Integration

### 8.1 WHO 2021 CNS Tumor Classification

The 2021 WHO classification requires integrated molecular diagnosis for many tumor types:

| Tumor Type | Required Molecular | Treatment Implications |
|---|---|---|
| **Glioblastoma, IDH-wildtype** | IDH1/2 wt, +TERT mutation OR +EGFR amp OR +7/−10 | Stupp protocol: maximal resection + RT 60Gy + TMZ |
| **Astrocytoma, IDH-mutant** | IDH1/2 mut, ATRX loss, no 1p19q codel | Grade-dependent; RT + TMZ for grade 3-4 |
| **Oligodendroglioma, IDH-mutant, 1p19q-codeleted** | IDH1/2 mut + 1p19q codel | PCV or TMZ + RT; better prognosis |
| **H3 K27-altered diffuse midline glioma** | H3K27M mutation | Radiation; poor prognosis |
| **Medulloblastoma** | WNT, SHH, Group 3, Group 4 | Subgroup-specific therapy |

### 8.2 MGMT Methylation

O6-methylguanine-DNA methyltransferase (MGMT) promoter methylation is the strongest predictive biomarker for temozolomide response in glioblastoma:
- **Methylated (~40%):** Better response to TMZ, median OS ~21 months
- **Unmethylated (~60%):** Poorer response, median OS ~14 months
- Testing: Methylation-specific PCR or pyrosequencing

### 8.3 RANO Criteria (Response Assessment in Neuro-Oncology)

| Response | Criteria |
|---|---|
| Complete response | No enhancing disease, stable/improved FLAIR, no new lesions |
| Partial response | >= 50% decrease in enhancing tumor area |
| Stable disease | < 50% decrease and < 25% increase |
| Progressive disease | >= 25% increase or new enhancing lesion |

---

## 9. Advanced EEG Pattern Recognition

### 9.1 ACNS Standardized Terminology (2021)

| Pattern | Description | Clinical Significance |
|---|---|---|
| **LPDs (lateralized periodic discharges)** | Periodic sharp waves, one hemisphere | Acute structural lesion, may be ictal |
| **GPDs (generalized periodic discharges)** | Bilateral periodic discharges | Metabolic encephalopathy, CJD, anoxia |
| **LRDA (lateralized rhythmic delta)** | Rhythmic delta, one hemisphere | Strongly associated with seizures |
| **GRDA (generalized rhythmic delta)** | Bilateral rhythmic delta | Encephalopathy, increased ICP |
| **Burst suppression** | Bursts of activity alternating with suppression | Severe encephalopathy, anesthesia, hypothermia |
| **Triphasic waves** | Anterior-predominant, positive-negative-positive | Hepatic encephalopathy (classic), also metabolic, CJD |
| **3 Hz spike-and-wave** | Generalized, rhythmic | Typical absence seizures |
| **Hypsarrhythmia** | Chaotic, high-amplitude, multifocal spikes | West syndrome (infantile spasms) |
| **PLEDs (now LPDs)** | Periodic lateralized epileptiform discharges | Herpes encephalitis, acute stroke, tumor |
| **Alpha coma** | Unreactive alpha, diffuse, non-posterior | Brainstem lesion, anoxia (poor prognosis) |

### 9.2 Ictal vs. Interictal Patterns

**Interictal:** Epileptiform discharges (spikes, sharp waves) between seizures. Localizing but not diagnostic of ongoing seizure.

**Ictal:** Evolving rhythmic activity representing an electrographic seizure. Requires:
- Definite evolution in frequency, morphology, and distribution
- Duration typically > 10 seconds
- Clinical correlation (may be subclinical/NCSE)

### 9.3 Non-Convulsive Status Epilepticus (NCSE)

NCSE is a medical emergency that is underdiagnosed because the patient may not convulse. Suspect in:
- Unexplained altered mental status
- Subtle motor signs (eye deviation, nystagmus, facial twitching)
- Failure to improve after convulsive status epilepticus

**Salzburg Criteria for NCSE:**
- >= 2.5 Hz epileptiform discharges, OR
- Epileptiform discharges < 2.5 Hz + clinical improvement with IV benzodiazepine, OR
- Subtle clinical signs + EEG correlate

---

## 10. Neuromuscular Junction: Advanced Diagnostics

### 10.1 Repetitive Nerve Stimulation (RNS)

| Disease | Stimulation Rate | Response |
|---|---|---|
| **Myasthenia gravis** | 2-3 Hz (low rate) | Decremental response > 10% |
| **Lambert-Eaton** | 2-3 Hz then 50 Hz | Decrement at low rate; > 100% increment at high rate |
| **Botulism** | 2-3 Hz then 50 Hz | Similar to LEMS pattern |

### 10.2 Single-Fiber EMG (SFEMG)

Most sensitive test for NMJ disorders:
- **Jitter:** Variability in time between two muscle fiber potentials from same motor unit
- **Blocking:** Failure of one fiber to fire
- Abnormal jitter + blocking = NMJ dysfunction
- 95% sensitivity for MG (even in ocular MG)

### 10.3 Antibody Panels

| Antibody | Disease | Frequency | Treatment Implication |
|---|---|---|---|
| **AChR** | Generalized MG | 85% | Cholinesterase inhibitors, immunosuppression, thymectomy |
| **MuSK** | MG (often bulbar-predominant) | 5-10% | Rituximab preferred; avoid cholinesterase inhibitors in some |
| **LRP4** | MG | 1-3% | Similar to AChR-MG management |
| **VGCC** | Lambert-Eaton | >90% | 3,4-DAP; screen for SCLC |
| **Ganglioside (GM1, GQ1b)** | GBS subtypes | Variable | IVIg or PLEX |
| **SOD1** | Familial ALS | ~20% of fALS | Tofersen (ASO therapy) |

---

## 11. Sleep-Movement Disorder Connection

### 11.1 RBD as a Prodromal Marker

REM sleep behavior disorder (dream enactment behavior with loss of REM atonia) is the strongest prodromal marker for synucleinopathies:
- **Conversion rate:** ~75-90% will develop PD, DLB, or MSA within 10-15 years
- **Confirmed by polysomnography:** REM without atonia (RWA)
- **Clinical significance:** Patients with isolated RBD should be counseled and monitored
- **Alpha-synuclein seed amplification assay (SAA):** Can detect alpha-synuclein pathology in CSF of RBD patients before motor onset

### 11.2 Sleep Disorders in Parkinson's Disease

| Disorder | Prevalence in PD | Mechanism |
|---|---|---|
| RBD | 25-50% | Brainstem nuclei degeneration |
| Excessive daytime sleepiness | 30-50% | Medication effect + hypothalamic involvement |
| Insomnia | 30-80% | Multifactorial (motor, nocturia, depression) |
| RLS/PLMD | 15-30% | Dopaminergic dysfunction |
| Sleep apnea | 20-60% | Upper airway dysfunction, autonomic |

---

## 12. Precision Medicine in Neurology

### 12.1 Gene-Specific Therapies (Approved and Pipeline)

| Gene/Target | Disease | Therapy | Mechanism |
|---|---|---|---|
| **SOD1** | Familial ALS | Tofersen (Qalsody) | Antisense oligonucleotide |
| **SMN1/2** | SMA | Nusinersen (Spinraza), Risdiplam (Evrysdi), Onasemnogene (Zolgensma) | ASO, SMN2 splicing modifier, gene replacement |
| **SCN1A** | Dravet syndrome | Fenfluramine, cannabidiol | Serotonergic, multi-mechanism |
| **HTT** | Huntington disease | Tominersen (suspended), other ASOs in development | Huntingtin lowering |
| **APOE** | Alzheimer's | Gene editing approaches (preclinical) | Risk modification |
| **LRRK2** | Parkinson's | LRRK2 kinase inhibitors (Phase II) | Kinase inhibition |
| **GBA1** | Parkinson's/GD | Venglustat, ambroxol (Phase II/III) | GCase enhancement |

### 12.2 Biomarker-Driven Treatment Selection

The agent encodes biomarker-to-treatment pathways:
- Amyloid PET positive --> lecanemab/donanemab eligibility
- SOD1 mutation confirmed --> tofersen initiation
- MGMT methylated GBM --> temozolomide benefit predicted
- JCV antibody positive --> natalizumab risk stratification
- SMN2 copy number --> SMA therapy selection and prognosis

---

*Neurology Intelligence Agent -- Learning Guide: Advanced Topics v1.3.0*
*HCLS AI Factory / GTC Europe 2026*
