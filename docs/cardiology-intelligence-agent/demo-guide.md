# Cardiology Intelligence Agent -- Demo Guide

**Version:** 1.0.0
**Author:** Adam Jones
**Date:** March 2026
**License:** Apache 2.0

---

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Starting the System](#2-starting-the-system)
3. [Verifying Health](#3-verifying-health)
4. [Demo Scenarios: 10-Tab Walkthrough](#4-demo-scenarios-10-tab-walkthrough)
5. [Sample API Queries](#5-sample-api-queries)
6. [Full Demo Scenarios](#6-full-demo-scenarios)
7. [Demo Talking Points](#7-demo-talking-points)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Prerequisites

### Hardware

- NVIDIA DGX Spark (recommended) or any system with 16GB+ RAM
- GPU optional (accelerates embedding generation)

### Software

- Docker and Docker Compose v2+
- Python 3.12+ (for manual setup without Docker)
- Anthropic API key with Claude Sonnet 4.6 access

### Environment Setup

```bash
cd ai_agent_adds/cardiology_intelligence_agent

# Copy environment template and add your Anthropic API key
cp .env.example .env
# Edit .env and set ANTHROPIC_API_KEY=sk-ant-...
```

---

## 2. Starting the System

### Option A: Docker Compose (Recommended)

```bash
# Start all 6 services (etcd, MinIO, Milvus, Streamlit, API, setup)
docker compose up -d

# Watch the setup process (creates 12 collections, seeds knowledge graph)
docker compose logs -f cardio-setup

# Wait for "Cardio Setup complete!" message, then verify:
docker compose ps
```

Expected output: 5 running services (cardio-setup exits after completion).

### Option B: Manual Setup

```bash
# 1. Ensure Milvus is running on localhost:19530

# 2. Install dependencies
pip install -r requirements.txt

# 3. Create collections and seed data
python scripts/setup_collections.py --drop-existing --seed
python scripts/seed_knowledge.py

# 4. Start the API server (port 8126)
uvicorn api.main:app --host 0.0.0.0 --port 8126 &

# 5. Start the Streamlit UI (port 8536)
streamlit run app/cardio_ui.py --server.port 8536
```

---

## 3. Verifying Health

### API Health Check

```bash
curl http://localhost:8126/health | python -m json.tool
```

Expected response:
```json
{
  "status": "healthy",
  "collections": 13,
  "total_vectors": 0,
  "api_port": 8126,
  "streamlit_port": 8536,
  "rate_limit": "100/min per IP",
  "auth": "X-API-Key required"
}
```

### Collection Status

```bash
curl http://localhost:8126/v1/cardio/collections | python -m json.tool
```

### Streamlit UI

Open http://localhost:8536 in a browser. The Dashboard tab should show collection statistics and system health.

---

## 4. Demo Scenarios: 10-Tab Walkthrough

### Tab 1: Dashboard

**What it shows:** System overview with collection statistics, recent query history, service health status, and quick links to all workflows.

**Demo action:** Point out the 13 collections, vector counts per collection, and overall system health indicators.

**Key talking point:** "This dashboard gives a real-time view of the entire cardiovascular knowledge base -- 12 specialized collections plus 3.5 million genomic variants, all searchable in under 5 seconds."

### Tab 2: CAD Assessment

**Sample input:**
```
Patient is a 58-year-old male smoker with exertional chest pain.
Coronary CTA shows CAD-RADS 3 with mixed plaque in the LAD.
Calcium score: Agatston 342.
```

**Expected output:**
- CAD-RADS 3 interpretation (moderate stenosis, 50-69%)
- Calcium score percentile for age/sex
- Recommendation for functional testing (FFR-CT or stress MPI)
- Guideline citation: ACC/AHA 2021 Chest Pain Guidelines
- Risk enhancement factors identified
- Statin therapy recommendation

### Tab 3: Heart Failure

**Sample input:**
```
62-year-old female with LVEF 28%, NYHA class III symptoms.
Currently on metoprolol succinate 50mg daily and lisinopril 10mg daily.
Not on MRA or SGLT2i. Labs: K+ 4.2, eGFR 55, NT-proBNP 2,400.
```

**Expected output:**
- EF classification: HFrEF (LVEF <=40%)
- NYHA/ACC staging: Stage C, NYHA III
- GDMT optimization:
  - Switch lisinopril to sacubitril/valsartan (start 24/26mg BID, titrate to 97/103mg BID)
  - Add spironolactone 12.5-25mg (check K+ at 1 week)
  - Add dapagliflozin 10mg daily
  - Uptitrate metoprolol to 200mg daily target
- Device evaluation: Consider ICD (LVEF <=35% on GDMT >=90 days)
- MAGGIC risk score calculation

### Tab 4: Valvular Disease

**Sample input:**
```
78-year-old male with severe aortic stenosis.
Echo: Vmax 4.5 m/s, mean gradient 52 mmHg, AVA 0.7 cm2.
Symptomatic with exertional dyspnea and presyncope.
STS predicted mortality 3.2%.
```

**Expected output:**
- Severity classification: Severe AS (all criteria met)
- Intervention criteria: Met (symptomatic severe AS)
- Recommendation: AVR indicated (Class I, Level A)
- TAVR vs SAVR discussion based on age and surgical risk
- EuroSCORE II calculation
- Guideline citation: ACC/AHA 2020 VHD Guidelines

### Tab 5: Arrhythmia

**Sample input:**
```
72-year-old female with newly diagnosed atrial fibrillation.
History of hypertension and diabetes. No prior stroke.
No CHF, no vascular disease.
```

**Expected output:**
- CHA2DS2-VASc calculation:
  - CHF: 0, Hypertension: 1, Age >=75: 1, Diabetes: 1, Stroke: 0, Vascular: 0, Age 65-74: 0 (already counted >=75), Female: 1
  - Score: 4 (annual stroke rate ~4.0%)
- Anticoagulation recommended (score >=3 for females)
- DOAC preferred over warfarin (Class I, Level A)
- HAS-BLED assessment for bleeding risk
- Rate vs rhythm control discussion
- Guideline citation: ACC/AHA/HRS 2023 AF Guidelines

### Tab 6: Cardiac MRI

**Sample input:**
```
35-year-old male with unexplained heart failure.
Cardiac MRI shows LVEF 32%, dilated LV (LVIDd 68mm).
Mid-wall linear LGE in the interventricular septum.
Native T1 elevated at 1120ms (1.5T). No ischemic LGE pattern.
```

**Expected output:**
- LGE pattern: Mid-wall (non-ischemic pattern)
- Differential diagnosis: Dilated cardiomyopathy, LMNA-related cardiomyopathy
- Cross-modal trigger fired: DCM gene panel recommended
  - Genes: TTN, LMNA, RBM20, MYH7, DSP, FLNC, BAG3
- LMNA highlighted: mid-wall septal LGE strongly associated with LMNA mutations
- If LMNA positive: early ICD consideration regardless of LVEF threshold
- GDMT initiation recommendations
- Guideline citation: ESC 2023 Cardiomyopathy Guidelines

### Tab 7: Stress Test

**Sample input:**
```
55-year-old male. Exercise stress echo performed.
Achieved 10 METs, peak HR 158 (96% predicted).
New hypokinesis of the mid-anterior and apical segments at peak stress.
Resting LVEF 60%, no significant change with stress.
```

**Expected output:**
- Duke Treadmill Score calculation
- Wall motion abnormality interpretation: Inducible ischemia in LAD territory
- Risk stratification: Intermediate-to-high risk based on extent of ischemia
- Recommendation: Consider coronary angiography for definitive evaluation
- Medical therapy optimization discussion
- Guideline citation: ACC/AHA 2021 Chest Pain Guidelines

### Tab 8: Prevention

**Sample input:**
```
52-year-old White male. Non-smoker, no diabetes.
Total cholesterol 240, HDL 42, LDL 158.
SBP 138, on amlodipine. No prior ASCVD.
Father had MI at age 51.
```

**Expected output:**
- ASCVD 10-year risk calculation (using Pooled Cohort Equations)
- Risk enhancers identified: Family history of premature ASCVD, low HDL
- Statin eligibility assessment
- If intermediate risk: Consider CAC scoring for shared decision-making
- Lifestyle modification recommendations
- LDL target discussion based on risk category
- Guideline citation: ACC/AHA 2018 Cholesterol Guidelines

### Tab 9: Cardio-Oncology

**Sample input:**
```
48-year-old female with breast cancer, planned doxorubicin + trastuzumab.
Baseline LVEF 62%, GLS -21%.
No prior cardiac history. HTN controlled on lisinopril.
Cumulative doxorubicin dose planned: 240 mg/m2.
```

**Expected output:**
- Baseline cardiac assessment: Adequate (LVEF >50%, GLS normal)
- Cardiotoxicity risk stratification: High risk (anthracycline + trastuzumab)
- Monitoring schedule:
  - Echo with GLS: baseline, after 2 cycles, after 4 cycles, 3 months post-completion
  - Troponin: before each cycle
  - NT-proBNP: baseline and after 4 cycles
- Cardioprotective measures:
  - Continue ACEi (lisinopril) -- cardioprotective
  - Consider dexrazoxane if cumulative dose >300 mg/m2
- CTRCD definitions:
  - LVEF decline >10% to below 50%
  - GLS decline >15% from baseline
- Guideline citation: ESC 2022 Cardio-Oncology Guidelines

### Tab 10: Reports

**What it shows:** Report generation interface with export options.

**Demo action:**
1. Select a previous query result
2. Choose export format (PDF, CSV, FHIR)
3. Generate and download the report

**Key talking point:** "Every clinical analysis can be exported as a structured PDF report with full citations, as CSV data for further analysis, or as FHIR R4 resources for EHR integration."

---

## 5. Sample API Queries

### 5.1 General Cardiology Query

```bash
curl -X POST http://localhost:8126/v1/cardio/query \
  -H "Content-Type: application/json" \
  -d '{
    "question": "What is the recommended GDMT for a patient with HFrEF and LVEF of 25%?",
    "patient_context": {
      "age": 65,
      "sex": "male",
      "lvef": 25,
      "nyha_class": "nyha_iii",
      "medications": ["metoprolol 25mg", "lisinopril 5mg"],
      "labs": {"potassium": 4.1, "creatinine": 1.2, "egfr": 58}
    }
  }' | python -m json.tool
```

### 5.2 ASCVD Risk Calculation

```bash
curl -X POST http://localhost:8126/v1/cardio/risk/ascvd \
  -H "Content-Type: application/json" \
  -d '{
    "score_type": "ascvd",
    "age": 55,
    "sex": "male",
    "race": "white",
    "total_cholesterol": 213,
    "hdl": 50,
    "systolic_bp": 120,
    "hypertension_treatment": true,
    "diabetes": false,
    "smoker": false
  }' | python -m json.tool
```

Expected response:
```json
{
  "score_type": "ascvd",
  "score_value": 8.4,
  "risk_category": "intermediate",
  "interpretation": "10-year ASCVD risk of 8.4% (intermediate risk, 7.5-20%). Consider risk enhancers and shared decision-making for statin therapy.",
  "recommendations": [
    "Discuss statin therapy initiation",
    "Consider CAC scoring for risk reclassification",
    "Lifestyle modification: diet, exercise, weight management",
    "Assess risk-enhancing factors: family history, Lp(a), hsCRP, ABI"
  ],
  "guideline_reference": "2018 AHA/ACC Cholesterol Guideline; Goff 2013 PCE"
}
```

### 5.3 CHA2DS2-VASc Score

```bash
curl -X POST http://localhost:8126/v1/cardio/risk/cha2ds2-vasc \
  -H "Content-Type: application/json" \
  -d '{
    "score_type": "cha2ds2_vasc",
    "age": 72,
    "sex": "female",
    "congestive_heart_failure": false,
    "hypertension_treatment": true,
    "diabetes": true,
    "history_of_stroke": false,
    "vascular_disease": false
  }' | python -m json.tool
```

### 5.4 HAS-BLED Score

```bash
curl -X POST http://localhost:8126/v1/cardio/risk/has-bled \
  -H "Content-Type: application/json" \
  -d '{
    "score_type": "has_bled",
    "age": 72,
    "systolic_bp": 155,
    "renal_disease": false,
    "liver_disease": false,
    "history_of_stroke": false,
    "history_of_bleeding": false,
    "labile_inr": false,
    "alcohol_excess": false,
    "antiplatelet_nsaid": true
  }' | python -m json.tool
```

### 5.5 GDMT Optimization

```bash
curl -X POST http://localhost:8126/v1/cardio/gdmt/optimize \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Optimize GDMT for this patient",
    "patient_context": {
      "age": 58,
      "sex": "male",
      "lvef": 30,
      "nyha_class": "nyha_ii",
      "systolic_bp": 118,
      "heart_rate": 72,
      "labs": {
        "potassium": 4.3,
        "creatinine": 1.1,
        "egfr": 65,
        "nt_pro_bnp": 1800
      },
      "current_medications": [
        {"pillar": "beta_blocker", "drug_name": "carvedilol", "current_dose": "12.5mg BID", "status": "uptitrating"},
        {"pillar": "arni_acei_arb", "drug_name": "lisinopril", "current_dose": "10mg daily", "status": "initiated"},
        {"pillar": "mra", "drug_name": "", "status": "not_started"},
        {"pillar": "sglt2i", "drug_name": "", "status": "not_started"}
      ]
    }
  }' | python -m json.tool
```

### 5.6 CAD Assessment Workflow

```bash
curl -X POST http://localhost:8126/v1/cardio/workflow/cad \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Evaluate this patient with chest pain and positive CTA",
    "workflow_type": "cad_assessment",
    "patient_context": {
      "age": 62,
      "sex": "male",
      "symptoms": "exertional chest tightness for 3 months",
      "cta_findings": "CAD-RADS 4A, 80% stenosis proximal LAD, mixed plaque",
      "calcium_score": 520,
      "lvef": 55
    }
  }' | python -m json.tool
```

### 5.7 Heart Failure Workflow

```bash
curl -X POST http://localhost:8126/v1/cardio/workflow/heart-failure \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Classify and manage this heart failure patient",
    "workflow_type": "heart_failure",
    "patient_context": {
      "age": 68,
      "sex": "female",
      "lvef": 35,
      "nyha_class": "nyha_iii",
      "bnp": 850,
      "labs": {"potassium": 4.5, "creatinine": 1.4, "egfr": 42}
    }
  }' | python -m json.tool
```

### 5.8 Cardio-Oncology Workflow

```bash
curl -X POST http://localhost:8126/v1/cardio/workflow/cardio-onc \
  -H "Content-Type: application/json" \
  -d '{
    "question": "Assess cardiotoxicity risk and monitoring plan",
    "workflow_type": "cardio_oncology",
    "patient_context": {
      "age": 55,
      "sex": "female",
      "cancer_therapy": "doxorubicin 60mg/m2 x 4 cycles + trastuzumab",
      "baseline_lvef": 60,
      "current_lvef": 52,
      "baseline_gls": -20.5,
      "current_gls": -16.8,
      "troponin_trend": "rising from 0.01 to 0.08 ng/mL"
    }
  }' | python -m json.tool
```

### 5.9 PDF Report Generation

```bash
curl -X POST http://localhost:8126/v1/cardio/report/pdf \
  -H "Content-Type: application/json" \
  -d '{
    "query_id": "last",
    "include_risk_scores": true,
    "include_gdmt": true,
    "include_citations": true
  }' --output cardio_report.pdf
```

---

## 6. Full Demo Scenarios

### Scenario 1: CAD Patient Journey (5 minutes)

**Setup:** 58-year-old male smoker with exertional chest pain.

1. **Prevention Tab**: Calculate ASCVD risk (expected: elevated)
2. **CAD Tab**: Enter CTA findings (CAD-RADS 3, calcium score 342)
3. **API call**: POST /v1/cardio/risk/heart (HEART score for chest pain)
4. **Talking point**: "The system synthesizes risk scores, imaging findings, and guidelines into a unified recommendation -- something that typically requires consulting multiple reference sources."

### Scenario 2: HF Patient Needing GDMT (5 minutes)

**Setup:** 62-year-old female, LVEF 28%, NYHA III, on suboptimal GDMT.

1. **Heart Failure Tab**: Enter patient data (LVEF 28%, on low-dose metoprolol + lisinopril only)
2. **Observe GDMT optimization**: System recommends switching to ARNI, adding MRA + SGLT2i, uptitrating BB
3. **Risk Tab**: Calculate MAGGIC score for mortality risk
4. **Talking point**: "Studies show fewer than 25% of eligible patients are on target-dose GDMT. This optimizer identifies every gap and generates specific titration guidance with safety checks."

### Scenario 3: AF Patient Needing Anticoagulation (3 minutes)

**Setup:** 72-year-old female with newly diagnosed AF.

1. **Arrhythmia Tab**: Enter patient demographics and comorbidities
2. **CHA2DS2-VASc**: Automatic calculation (score 4, stroke rate 4.0%)
3. **HAS-BLED**: Automatic calculation for bleeding risk assessment
4. **Anticoagulation recommendation**: DOAC recommended with specific drug suggestions
5. **Talking point**: "Both scores are calculated simultaneously, and the recommendation weighs stroke risk against bleeding risk -- exactly as the 2023 AF guidelines instruct."

### Scenario 4: Cardio-Oncology Surveillance (5 minutes)

**Setup:** 48-year-old female with breast cancer starting doxorubicin + trastuzumab.

1. **Cardio-Oncology Tab**: Enter baseline cardiac assessment (LVEF 62%, GLS -21%)
2. **Risk stratification**: System identifies high cardiotoxicity risk
3. **Monitoring schedule**: Auto-generated per ESC 2022 guidelines
4. **Later visit**: Enter follow-up data showing GLS decline to -16.8% (>15% decline)
5. **CTRCD alert**: System flags subclinical cardiotoxicity and recommends cardiology consultation
6. **Talking point**: "The system tracks GLS decline -- a more sensitive marker than LVEF alone. A 15% relative decline triggers an alert before overt cardiac dysfunction occurs."

### Scenario 5: Cross-Modal Genomic Integration (3 minutes)

**Setup:** 35-year-old with unexplained LVH on echocardiogram.

1. **Cardiac MRI Tab**: Enter MRI findings (wall thickness 19mm, patchy LGE)
2. **Cross-modal trigger**: System recommends HCM gene panel + Fabry screening (18 genomic triggers available)
3. **Show genomic_evidence integration**: System queries 3.5M variant collection
4. **Talking point**: "This is the power of the HCLS AI Factory's three-stage pipeline. An imaging finding on the cardiology agent automatically triggers a query against the genomics pipeline's variant database -- connecting imaging to genomics in seconds."

### Scenario 6: Acute Decompensated Heart Failure (3 minutes)

**Setup:** 72-year-old male presenting with acute dyspnea, orthopnea, and lower extremity edema.

1. **Enter via API**: POST /v1/cardio/workflow/acute-dhf with hemodynamic data (PCWP 28, CI 1.8)
2. **Hemodynamic profiling**: System classifies as "Cold-Wet" (Profile C) using 14 hemodynamic parameters
3. **Management**: IV diuretic dosing, inotrope selection, mechanical support evaluation
4. **Talking point**: "The acute decompensated HF workflow integrates 14 hemodynamic parameters and 6 cath lab protocols to classify shock profiles and guide escalation decisions."

### Scenario 7: Post-MI Risk Management (3 minutes)

**Setup:** 55-year-old male, 48 hours post-STEMI with successful PCI to LAD.

1. **Enter via API**: POST /v1/cardio/workflow/post-mi with post-MI context
2. **Risk stratification**: LVEF assessment, arrhythmia risk, mechanical complications
3. **Secondary prevention**: Dual antiplatelet, high-intensity statin, beta-blocker, ACEi/ARB
4. **Talking point**: "The post-MI workflow generates a complete secondary prevention plan with 39 landmark trial references and specific follow-up imaging timing."

### Scenario 8: Myocarditis Workup (3 minutes)

**Setup:** 28-year-old male with chest pain and troponin elevation 2 weeks after viral illness.

1. **Enter via API**: POST /v1/cardio/workflow/myocarditis-pericarditis
2. **Lake Louise criteria**: System evaluates CMR findings against diagnostic criteria
3. **Activity restriction**: Evidence-based return-to-activity guidance
4. **Talking point**: "The myocarditis workflow applies the updated Lake Louise criteria and provides structured return-to-activity protocols based on current ACC/AHA recommendations."

---

## 7. Demo Talking Points

### Architecture Highlights

- "13 specialized Milvus vector collections, each tuned for a specific cardiovascular data type"
- "6 validated risk calculators using published coefficients -- not approximations"
- "4-pillar GDMT optimizer that checks every contraindication and generates specific titration steps"
- "Cross-modal triggers that bridge cardiac imaging to genomic testing"
- "All running on a $3,999 DGX Spark -- no cloud dependency, no data leaving the building"

### Clinical Value Propositions

- "End-to-end from patient question to guideline-grounded recommendation in under 5 seconds"
- "Every recommendation cites specific ACC/AHA/ESC guidelines with class and evidence level"
- "Connects to the HCLS AI Factory's 3.5 million genomic variants for integrated precision medicine"
- "Knowledge graph covers 45 conditions, 56 genes, 32 drug classes, 29 biomarkers"
- "Open source under Apache 2.0 -- any institution can deploy and customize"

### Differentiators from Existing Tools

- "Unlike PubMed: semantic search across structured clinical data, not just keywords"
- "Unlike UpToDate: patient-specific reasoning with computed risk scores"
- "Unlike commercial CVIS: $3,999 not $500K, with genomic integration"
- "Unlike general AI: validated calculators, citation provenance, no hallucination of clinical values"

---

## 8. Troubleshooting

### Milvus Connection Failed

```bash
# Check if Milvus is running
curl http://localhost:9091/healthz

# If not running:
docker compose restart milvus

# Wait 30 seconds, then check again
curl http://localhost:9091/healthz
```

### API Not Responding

```bash
# Check API logs
docker compose logs cardio-api

# Restart API service
docker compose restart cardio-api
```

### Streamlit UI Not Loading

```bash
# Check Streamlit logs
docker compose logs cardio-ui

# Verify port is not in use
lsof -i :8536

# Restart Streamlit
docker compose restart cardio-ui
```

### LLM Responses Not Working

```bash
# Verify API key is set
echo $CARDIO_ANTHROPIC_API_KEY

# Test Anthropic connectivity
curl https://api.anthropic.com/v1/messages \
  -H "x-api-key: $CARDIO_ANTHROPIC_API_KEY" \
  -H "anthropic-version: 2023-06-01" \
  -H "content-type: application/json" \
  -d '{"model": "claude-sonnet-4-6", "max_tokens": 10, "messages": [{"role": "user", "content": "Hello"}]}'
```

Note: Without an API key, the system operates in search-only mode -- vector search and risk calculators work, but LLM synthesis is unavailable.

### Collections Not Created

```bash
# Re-run setup manually
python scripts/setup_collections.py --drop-existing --seed
python scripts/seed_knowledge.py
```
