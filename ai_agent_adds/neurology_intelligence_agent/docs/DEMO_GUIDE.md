# Neurology Intelligence Agent -- Demo Guide

**Version:** 1.0.0
**Date:** 2026-03-22
**Author:** Adam Jones

---

## Pre-Demo Setup

### Verify Services

```bash
# Confirm all services running
docker compose ps

# Verify health
curl -s http://localhost:8528/health | python -m json.tool

# Expected: "status": "healthy", "collections": 14, "workflows": 9, "scales": 10

# Open Streamlit UI
open http://localhost:8529

# Open API docs
open http://localhost:8528/docs
```

### Demo Environment Checklist

- [ ] All 14 collections loaded (check `/collections`)
- [ ] Streamlit UI responsive at port 8529
- [ ] API health returns "healthy"
- [ ] Anthropic API key configured and LLM responding
- [ ] Browser tabs ready for Streamlit UI and API docs
- [ ] Fallback slides prepared

---

## Demo Scenario 1: Acute Stroke Triage

**Clinical Setup:** Emergency department. 72-year-old woman presents with sudden-onset left-sided weakness and speech difficulty. Last known well 2 hours ago. CT head negative for hemorrhage. CTA shows right M1 MCA occlusion.

### Step 1: Calculate NIHSS

Via API:
```bash
curl -X POST http://localhost:8528/v1/neuro/scale/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "scale_type": "nihss",
    "scores": {
      "1a_loc": 0,
      "1b_loc_questions": 1,
      "1c_loc_commands": 0,
      "2_gaze": 1,
      "3_visual": 2,
      "4_facial": 2,
      "5a_left_arm": 4,
      "5b_right_arm": 0,
      "6a_left_leg": 3,
      "6b_right_leg": 0,
      "7_ataxia": 0,
      "8_sensory": 1,
      "9_language": 2,
      "10_dysarthria": 2,
      "11_extinction": 1
    }
  }'
```

**Expected output:** NIHSS = 19 (Moderate-to-severe stroke). Recommendations: IV alteplase, CTA for LVO evaluation, consider mechanical thrombectomy.

### Step 2: Calculate ASPECTS

```bash
curl -X POST http://localhost:8528/v1/neuro/scale/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "scale_type": "aspects",
    "affected_regions": ["I", "M2"]
  }'
```

**Expected output:** ASPECTS = 8/10 (Favorable for intervention). Two regions with early ischemic change.

### Step 3: Run Stroke Triage Workflow

Via Streamlit: Type the following query:

> "72-year-old woman with acute left hemiparesis and aphasia, NIHSS 19, last known well 2 hours ago, CT negative for hemorrhage, CTA showing right M1 MCA occlusion. ASPECTS 8. Is she eligible for tPA and thrombectomy?"

**Demo talking points:**
- NIHSS 19 indicates moderate-to-severe stroke -- tPA strongly recommended
- Within 4.5-hour window -- tPA eligible (no contraindications mentioned)
- M1 occlusion with ASPECTS >= 6 -- thrombectomy indicated
- Door-to-needle target: < 60 minutes for tPA
- Door-to-groin target: < 90 minutes for thrombectomy
- Agent cites AHA/ASA 2019 guidelines and DAWN/DEFUSE-3 trial criteria

---

## Demo Scenario 2: Memory Clinic -- Alzheimer's Evaluation

**Clinical Setup:** Outpatient memory clinic. 68-year-old man with 18 months of progressive memory loss. Retired professor. Wife reports word-finding difficulty and getting lost while driving.

### Step 1: Calculate MoCA

```bash
curl -X POST http://localhost:8528/v1/neuro/scale/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "scale_type": "moca",
    "domain_scores": {
      "visuospatial": 3,
      "naming": 3,
      "attention": 4,
      "language": 2,
      "abstraction": 1,
      "delayed_recall": 1,
      "orientation": 6
    },
    "education_years": 20
  }'
```

**Expected output:** MoCA = 20/30 (Mild cognitive impairment). Notably impaired: delayed_recall, abstraction. Recommendations: neuropsychological testing, structural MRI, consider amyloid PET or CSF biomarkers.

### Step 2: Query Dementia Evaluation

Via Streamlit:

> "68-year-old man with 18 months progressive memory loss, MoCA 20/30 with impaired delayed recall and abstraction. APOE e3/e4 genotype. MRI shows mild bilateral hippocampal atrophy. What is the ATN staging and is he a candidate for anti-amyloid therapy?"

**Demo talking points:**
- MoCA 20 suggests MCI-to-mild dementia range
- APOE e3/e4 increases AD risk (3x vs e3/e3)
- Hippocampal atrophy suggests neurodegeneration (N+)
- Recommend amyloid PET or CSF Abeta42/p-tau to determine A and T status
- If A+T+, patient may be candidate for lecanemab (CLARITY AD trial criteria: MoCA 18-26, amyloid-positive)
- Agent generates differential: AD most likely, but consider DLB, vascular, FTD

---

## Demo Scenario 3: Drug-Resistant Epilepsy

**Clinical Setup:** Epilepsy clinic. 28-year-old woman with temporal lobe epilepsy since age 14. Failed levetiracetam, lamotrigine, and carbamazepine. Seizure frequency: 4-6 focal impaired awareness seizures per month.

### Step 1: Query Epilepsy Classification

Via Streamlit:

> "28-year-old woman with drug-resistant temporal lobe epilepsy. EEG shows right anterior temporal sharp waves and TIRDA. MRI shows right mesial temporal sclerosis. Failed levetiracetam, lamotrigine, and carbamazepine. 4-6 focal impaired awareness seizures per month with deja vu aura and oral automatisms. Is she a surgical candidate?"

**Demo talking points:**
- Meets ILAE definition of drug-resistant epilepsy (failed 2+ appropriate ASMs)
- TLE with hippocampal sclerosis is the most surgically remediable epilepsy syndrome
- EEG-MRI concordance (right temporal spikes + right MTS) is favorable for surgery
- Anterior temporal lobectomy offers 60-80% seizure freedom rate
- Additional pre-surgical workup: video-EEG monitoring, neuropsychological testing, Wada test
- Alternative: LITT (laser interstitial thermal therapy) or RNS
- Agent identifies contraindicated ASMs for other syndromes (educational contrast)

### Step 2: Scale Calculators (Show Flexibility)

Demonstrate other scale calculators to show breadth:

```bash
# ALSFRS-R for ALS patient
curl -X POST http://localhost:8528/v1/neuro/scale/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "scale_type": "alsfrs_r",
    "scores": {
      "speech": 3, "salivation": 3, "swallowing": 3,
      "handwriting": 2, "cutting_food": 2, "dressing_hygiene": 2,
      "turning_in_bed": 3, "walking": 2, "climbing_stairs": 1,
      "dyspnea": 3, "orthopnea": 3, "respiratory_insufficiency": 3
    },
    "months_since_onset": 12
  }'
```

---

## Demo Scenario 4: New Brain Mass

**Clinical Setup:** Neurosurgery referral. 55-year-old man with 3-week history of progressive headache, new-onset seizure. MRI shows 4 cm ring-enhancing left temporal lobe mass with surrounding edema.

### Step 1: Query Brain Tumor Evaluation

Via Streamlit:

> "55-year-old man with new-onset seizure, progressive headache. MRI shows 4 cm ring-enhancing mass in left temporal lobe with significant vasogenic edema and midline shift. What is the differential diagnosis and what molecular workup is needed?"

**Demo talking points:**
- Top differential: glioblastoma (IDH-wildtype), brain metastasis, CNS lymphoma, abscess
- Critical molecular markers to determine: IDH mutation status, MGMT methylation, TERT promoter, EGFR amplification
- WHO 2021 classification requires integrated molecular diagnosis
- If GBM (IDH-wt): Stupp protocol (maximal safe resection + RT 60 Gy + temozolomide + TTFields)
- If IDH-mutant astrocytoma: different prognosis and treatment approach
- MGMT methylation predicts temozolomide response
- Agent cites NCCN CNS guidelines and WHO 2021 classification

---

## Demo Scenario 5: MS Disease Monitoring

**Clinical Setup:** MS clinic follow-up. 34-year-old woman with RRMS diagnosed 3 years ago. On dimethyl fumarate (Tecfidera). New MRI shows 3 new T2 lesions and 1 Gd-enhancing lesion. EDSS 3.0.

### Step 1: Calculate EDSS

```bash
curl -X POST http://localhost:8528/v1/neuro/scale/calculate \
  -H "Content-Type: application/json" \
  -d '{
    "scale_type": "edss",
    "fs_scores": {
      "visual": 0,
      "brainstem": 1,
      "pyramidal": 3,
      "cerebellar": 1,
      "sensory": 2,
      "bowel_bladder": 0,
      "cerebral": 0
    },
    "edss_step": 3.0
  }'
```

### Step 2: Query MS Monitoring

Via Streamlit:

> "34-year-old woman with RRMS on dimethyl fumarate for 2 years. New MRI shows 3 new T2 lesions and 1 gadolinium-enhancing lesion. EDSS 3.0. JCV antibody index 0.9. NfL elevated at 28 pg/mL. Should we escalate DMT? What are the options?"

**Demo talking points:**
- NOT meeting NEDA-3: new T2 lesions = disease activity on current DMT
- Elevated NfL supports active disease
- DMT escalation indicated: switch from moderate to high-efficacy
- Options: ocrelizumab (Ocrevus), natalizumab (Tysabri), ofatumumab (Kesimpta)
- JCV index 0.9 = moderate PML risk with natalizumab (requires monitoring q6 months)
- Ocrelizumab may be preferred given JCV status
- Agent explains NEDA-3 criteria and DMT tier classification
- Agent cites McDonald 2017 criteria and AAN DMT guidelines

---

## Demo Tips

1. **Start with stroke** -- most dramatic, time-critical, clear decision points
2. **Show scale calculators** -- they respond instantly and demonstrate clinical accuracy
3. **Use the Streamlit UI** for natural language queries (more engaging for audience)
4. **Use the API docs** (Swagger) to show the full endpoint catalog
5. **Highlight cross-collection search** -- mention the 14 collections being searched in parallel
6. **Point out citations** -- the agent cites specific guidelines and trials
7. **If LLM is slow**, use the scale calculators as a bridge (they are instant)
8. **Fallback**: if Anthropic API is down, switch to search-only mode and show raw retrieval results

---

## Quick Reference: Scale Calculator API

| Scale | POST Path | Required Fields |
|---|---|---|
| NIHSS | `/v1/neuro/scale/calculate` | `scale_type: "nihss"`, `scores: {}` |
| GCS | `/v1/neuro/scale/calculate` | `scale_type: "gcs"`, `eye`, `verbal`, `motor` |
| MoCA | `/v1/neuro/scale/calculate` | `scale_type: "moca"`, `domain_scores: {}` |
| UPDRS | `/v1/neuro/scale/calculate` | `scale_type: "updrs_part_iii"`, `scores: {}` |
| EDSS | `/v1/neuro/scale/calculate` | `scale_type: "edss"`, `fs_scores: {}`, `edss_step` |
| mRS | `/v1/neuro/scale/calculate` | `scale_type: "mrs"`, `score` |
| HIT-6 | `/v1/neuro/scale/calculate` | `scale_type: "hit6"`, `responses: []` |
| ALSFRS-R | `/v1/neuro/scale/calculate` | `scale_type: "alsfrs_r"`, `scores: {}` |
| ASPECTS | `/v1/neuro/scale/calculate` | `scale_type: "aspects"`, `affected_regions: []` |
| Hoehn-Yahr | `/v1/neuro/scale/calculate` | `scale_type: "hoehn_yahr"`, `stage` |

---

*Neurology Intelligence Agent -- Demo Guide v1.3.0*
*HCLS AI Factory / GTC Europe 2026*
